# Independent reference implementation for the complete-data static binary
# switcher/stayer estimators and influence functions.
#
# Used exclusively by tests to verify R/static_complete_core.R against:
#   1. Full stacked stayer regression (`stats::lm.wfit` and explicit normal
#      matrix `M = (1/n) sum_c sum_i w_i C_{ic} a_{ic} a_{ic}'`) rather than
#      compact QR/intercept-partialling helpers.
#   2. Direct weighted logit (`stats::glm.fit`) and direct matrix `solve()`
#      for `H_c^{-1}` and `M^{-1}` in the first-step, cell-level (DR, IPW, OR),
#      and switcher-frequency aggregate influence functions.
#   3. Centered empirical unit-weight directional derivatives
#      `n * (T(1 + eps * d_i) - T(1 - eps * d_i)) / (2 * eps)` by refitting the
#      weighted logit and weighted full-design least-squares models under
#      perturbed original-unit weights.

ref_static_complete_eval <- function(cells,
                                     or_blocks,
                                     method,
                                     unit_weights = NULL,
                                     aggregate_cells = NULL,
                                     aggregate_signs = NULL,
                                     cluster = NULL,
                                     glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
                                     compute_if = TRUE,
                                     polish = FALSE) {
  cell_ids <- unname(vapply(cells, function(cl) cl$cell_id, character(1L)))
  names(cells) <- cell_ids
  n_cells <- length(cells)
  n <- length(cells[[1L]]$dy)
  n_num <- as.double(n)

  w <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  stopifnot(length(w) == n, all(is.finite(w)), all(w > 0))

  # Parse partition into block map
  if (is.list(or_blocks)) {
    block_ids <- names(or_blocks)
    if (is.null(block_ids)) block_ids <- paste0("block_", seq_along(or_blocks))
    cell_block <- stats::setNames(character(n_cells), cell_ids)
    for (b_idx in seq_along(or_blocks)) {
      mems <- or_blocks[[b_idx]]
      if (is.numeric(mems)) mems <- cell_ids[as.integer(mems)]
      cell_block[mems] <- block_ids[b_idx]
    }
  } else {
    ob_chr <- as.character(or_blocks)
    if (!is.null(names(or_blocks))) {
      ob_chr <- ob_chr[match(cell_ids, names(or_blocks))]
    }
    cell_block <- stats::setNames(ob_chr, cell_ids)
    block_ids <- unique(ob_chr)
  }
  blocks <- lapply(block_ids, function(bid) cell_ids[cell_block == bid])
  names(blocks) <- block_ids

  # Determine full design vector layout a_{ic} in R^{n_cells + sum_P p_P}:
  # First `n_cells` coordinates are local intercepts alpha_c; remaining coordinates
  # are block slope vectors b_P for P in block_ids.
  get_X_tilde <- function(cl) {
    if (!is.null(cl$X_tilde)) cl$X_tilde else cbind("(Intercept)" = rep(1.0, n), cl$Z)
  }
  get_Z_or <- function(cl) {
    if (!is.null(cl$Z_or)) {
      cl$Z_or
    } else if (!is.null(cl$X_or_tilde)) {
      cl$X_or_tilde[, -1L, drop = FALSE]
    } else if (!is.null(cl$Z)) {
      cl$Z
    } else {
      cl$X_tilde[, -1L, drop = FALSE]
    }
  }

  block_p <- vapply(block_ids, function(bid) {
    ncol(get_Z_or(cells[[blocks[[bid]][1L]]]))
  }, integer(1L))
  names(block_p) <- block_ids

  slope_offsets <- integer(length(block_ids))
  names(slope_offsets) <- block_ids
  cur_off <- n_cells
  for (bid in block_ids) {
    slope_offsets[[bid]] <- cur_off
    cur_off <- cur_off + block_p[[bid]]
  }
  dim_beta <- cur_off

  # Build per-cell full design matrix A_c (n x dim_beta) where row i is a_{ic}'.
  A_by_cell <- vector("list", n_cells)
  names(A_by_cell) <- cell_ids
  for (j in seq_len(n_cells)) {
    cid <- cell_ids[j]
    cl <- cells[[cid]]
    Z_or_c <- get_Z_or(cl)
    bid <- cell_block[[cid]]
    p_b <- block_p[[bid]]
    off_b <- slope_offsets[[bid]]

    A_c <- matrix(0.0, nrow = n, ncol = dim_beta)
    A_c[, j] <- 1.0
    if (p_b > 0L) {
      A_c[, off_b + seq_len(p_b)] <- Z_or_c
    }
    A_by_cell[[cid]] <- A_c
  }

  # 1. Weighted cell-specific logit fits on R_{ic} = 1
  eps_clip <- 1e-6
  cap_prob <- 1 - eps_clip
  kbar <- log(cap_prob / eps_clip)
  cap_odds <- cap_prob / (1 - cap_prob)

  canon_method <- if (identical(method, "or")) "reg" else method
  use_logit <- canon_method %in% c("ipw", "dr")
  use_or <- canon_method %in% c("reg", "dr")

  gamma_hat <- vector("list", n_cells)
  names(gamma_hat) <- cell_ids
  H_hat <- vector("list", n_cells)
  names(H_hat) <- cell_ids
  ell_gamma <- vector("list", n_cells)
  names(ell_gamma) <- cell_ids
  prob_by_cell <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  odds_by_cell <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  active_by_cell <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))

  if (use_logit) {
    for (cid in cell_ids) {
      cl <- cells[[cid]]
      X_c <- get_X_tilde(cl)
      r_idx <- which(cl$R == 1L)
      fit_g <- suppressWarnings(stats::glm.fit(
        x = X_c[r_idx, , drop = FALSE],
        y = as.double(cl$G[r_idx]),
        weights = w[r_idx],
        start = rep(0, ncol(X_c)),
        family = stats::binomial(link = "logit"),
        control = glm_control,
        intercept = FALSE
      ))
      g_c <- as.double(fit_g$coefficients)
      if (polish) {
        # Polish to machine-precision root of the weighted logit score so finite-difference
        # perturbations do not suffer from IRLS stopping-step quantization.
        X_r <- X_c[r_idx, , drop = FALSE]
        y_r <- as.double(cl$G[r_idx])
        w_r <- w[r_idx]
        for (polish_iter in 1:3) {
          u_tmp <- drop(X_r %*% g_c)
          p_tmp <- 1 / (1 + exp(-u_tmp))
          grad_tmp <- colSums(X_r * (w_r * (y_r - p_tmp)))
          if (max(abs(grad_tmp)) <= 1e-14) break
          hess_tmp <- crossprod(X_r, X_r * (w_r * p_tmp * (1 - p_tmp)))
          g_c <- g_c + drop(solve(hess_tmp, grad_tmp))
        }
      }
      names(g_c) <- colnames(X_c)
      gamma_hat[[cid]] <- g_c

      u_r <- drop(X_c[r_idx, , drop = FALSE] %*% g_c)
      p_r <- 1 / (1 + exp(-u_r))
      act_r <- as.double(p_r < cap_prob)
      r_r <- ifelse(act_r == 1, exp(u_r), cap_odds)

      prob_by_cell[r_idx, cid] <- p_r
      odds_by_cell[r_idx, cid] <- r_r
      active_by_cell[r_idx, cid] <- act_r

      if (compute_if) {
        w_hess_r <- w[r_idx] * p_r * (1 - p_r)
        H_c <- crossprod(X_c[r_idx, , drop = FALSE], X_c[r_idx, , drop = FALSE] * w_hess_r) / n_num
        H_hat[[cid]] <- H_c

        score_c <- matrix(0.0, nrow = n, ncol = ncol(X_c), dimnames = list(NULL, colnames(X_c)))
        score_c[r_idx, ] <- X_c[r_idx, , drop = FALSE] * (w[r_idx] * (cl$G[r_idx] - p_r))
        ell_gamma[[cid]] <- t(solve(H_c, t(score_c)))
      }
    }
  }

  # 2. Full stacked stayer least-squares outcome regression
  beta_hat <- rep(0.0, dim_beta)
  M_hat <- NULL
  ell_beta <- matrix(0.0, nrow = n, ncol = dim_beta)
  e_by_cell <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))

  if (use_or) {
    A_stack_list <- vector("list", n_cells)
    y_stack_list <- vector("list", n_cells)
    w_stack_list <- vector("list", n_cells)
    M_acc <- matrix(0.0, nrow = dim_beta, ncol = dim_beta)

    for (cid in cell_ids) {
      cl <- cells[[cid]]
      A_c <- A_by_cell[[cid]]
      c_idx <- which(cl$C == 1L)
      A_st <- A_c[c_idx, , drop = FALSE]
      y_st <- cl$dy[c_idx]
      w_st <- w[c_idx]

      A_stack_list[[cid]] <- A_st
      y_stack_list[[cid]] <- y_st
      w_stack_list[[cid]] <- w_st
      M_acc <- M_acc + crossprod(A_st, A_st * w_st) / n_num
    }
    A_stack <- do.call(rbind, A_stack_list)
    y_stack <- unlist(y_stack_list, use.names = FALSE)
    w_stack <- unlist(w_stack_list, use.names = FALSE)

    wfit <- stats::lm.wfit(x = A_stack, y = y_stack, w = w_stack)
    beta_hat <- as.double(wfit$coefficients)
    M_hat <- M_acc

    # Sum stayer regression score across all cells c in C for each original unit i:
    #   score_beta_i = sum_{c in C} w_i C_{ic} a_{ic} e_hat_{ic}
    score_beta <- matrix(0.0, nrow = n, ncol = dim_beta)
    for (cid in cell_ids) {
      cl <- cells[[cid]]
      A_c <- A_by_cell[[cid]]
      e_c <- cl$dy - drop(A_c %*% beta_hat)
      e_by_cell[, cid] <- e_c
      if (compute_if) {
        score_beta <- score_beta + A_c * (w * cl$C * e_c)
      }
    }
    if (compute_if) {
      ell_beta <- t(solve(M_hat, t(score_beta)))
    }
  }

  # 3. Cell contrasts and complete original-unit cell influence functions
  theta_hat <- stats::setNames(numeric(n_cells), cell_ids)
  tau_hat <- stats::setNames(numeric(n_cells), cell_ids)
  sigma_vec <- stats::setNames(numeric(n_cells), cell_ids)
  pi_hat <- stats::setNames(numeric(n_cells), cell_ids)
  Q_hat <- stats::setNames(rep(NA_real_, n_cells), cell_ids)
  mu_1 <- stats::setNames(numeric(n_cells), cell_ids)
  mu_0 <- stats::setNames(rep(NA_real_, n_cells), cell_ids)
  phi_theta <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  phi_tau <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))

  for (cid in cell_ids) {
    cl <- cells[[cid]]
    G_c <- as.double(cl$G)
    C_c <- as.double(cl$C)
    dy_c <- as.double(cl$dy)
    sigma_c <- if (!is.null(cl$sigma)) as.double(cl$sigma) else as.double(1L - 2L * as.integer(cl$d))
    sigma_vec[[cid]] <- sigma_c

    pi_c <- sum(w * G_c) / n_num
    pi_hat[[cid]] <- pi_c
    A_c <- A_by_cell[[cid]]
    X_c <- get_X_tilde(cl)

    if (identical(canon_method, "ipw")) {
      r_c <- odds_by_cell[, cid]
      act_c <- active_by_cell[, cid]
      Q_c <- sum(w * C_c * r_c) / n_num
      Q_hat[[cid]] <- Q_c

      mu_1c <- sum(w * G_c * dy_c) / (n_num * pi_c)
      mu_0c <- sum(w * C_c * r_c * dy_c) / (n_num * Q_c)
      theta_c <- mu_1c - mu_0c

      if (compute_if) {
        dot_r_c <- X_c * (act_c * r_c)
        k_gamma_c <- colSums(dot_r_c * (w * C_c * (dy_c - mu_0c))) / (n_num * Q_c)

        phi_c <- (w * G_c / pi_c) * (dy_c - mu_1c) -
          (w * C_c * r_c / Q_c) * (dy_c - mu_0c) -
          drop(ell_gamma[[cid]] %*% k_gamma_c)
      } else {
        phi_c <- numeric(n)
      }
    } else if (identical(canon_method, "reg")) {
      e_c <- e_by_cell[, cid]
      mu_1c <- sum(w * G_c * e_c) / (n_num * pi_c)
      mu_0c <- 0.0
      theta_c <- mu_1c

      if (compute_if) {
        k_beta_or_c <- colSums(A_c * (w * G_c)) / (n_num * pi_c)
        phi_c <- (w * G_c / pi_c) * (e_c - mu_1c) -
          drop(ell_beta %*% k_beta_or_c)
      } else {
        phi_c <- numeric(n)
      }
    } else {
      r_c <- odds_by_cell[, cid]
      act_c <- active_by_cell[, cid]
      Q_c <- sum(w * C_c * r_c) / n_num
      Q_hat[[cid]] <- Q_c
      e_c <- e_by_cell[, cid]

      mu_1c <- sum(w * G_c * e_c) / (n_num * pi_c)
      mu_0c <- sum(w * C_c * r_c * e_c) / (n_num * Q_c)
      theta_c <- mu_1c - mu_0c

      if (compute_if) {
        dot_r_c <- X_c * (act_c * r_c)
        k_gamma_c <- colSums(dot_r_c * (w * C_c * (e_c - mu_0c))) / (n_num * Q_c)
        k_beta_c <- colSums(A_c * (w * G_c)) / (n_num * pi_c) -
          colSums(A_c * (w * C_c * r_c)) / (n_num * Q_c)

        phi_c <- (w * G_c / pi_c) * (e_c - mu_1c) -
          (w * C_c * r_c / Q_c) * (e_c - mu_0c) -
          drop(ell_gamma[[cid]] %*% k_gamma_c) -
          drop(ell_beta %*% k_beta_c)
      } else {
        phi_c <- numeric(n)
      }
    }

    theta_hat[[cid]] <- theta_c
    tau_hat[[cid]] <- sigma_c * theta_c
    mu_1[[cid]] <- mu_1c
    mu_0[[cid]] <- mu_0c
    phi_theta[, cid] <- phi_c
    phi_tau[, cid] <- sigma_c * phi_c
  }

  cov_from_phi <- function(Phi, clus) {
    if (is.null(clus)) {
      crossprod(Phi) / (n_num^2)
    } else {
      u_c <- unique(clus)
      idx_c <- match(clus, u_c)
      Phi_c <- t(vapply(seq_along(u_c), function(j) {
        colSums(Phi[idx_c == j, , drop = FALSE])
      }, numeric(ncol(Phi))))
      if (ncol(Phi) == 1L) Phi_c <- matrix(Phi_c, ncol = 1L)
      crossprod(Phi_c) / (n_num^2)
    }
  }

  cov_theta <- cov_from_phi(phi_theta, cluster)
  cov_tau <- cov_from_phi(phi_tau, cluster)

  agg_res <- NULL
  if (!is.null(aggregate_cells)) {
    sel <- if (is.numeric(aggregate_cells)) cell_ids[as.integer(aggregate_cells)] else as.character(aggregate_cells)
    s_vec <- as.double(aggregate_signs)
    if (!is.null(names(aggregate_signs))) {
      s_vec <- unname(s_vec[match(sel, names(aggregate_signs))])
    }
    names(s_vec) <- sel

    pi_A <- pi_hat[sel]
    Pi_A <- sum(pi_A)
    omega_A <- pi_A / Pi_A
    theta_A <- theta_hat[sel]
    delta_hat <- sum(omega_A * s_vec * theta_A)

    G_A <- vapply(sel, function(cid) as.double(cells[[cid]]$G), numeric(n))
    if (n == 1L) G_A <- matrix(G_A, nrow = 1L)
    cell_term <- drop(phi_theta[, sel, drop = FALSE] %*% (omega_A * s_vec))
    G_cent <- sweep(G_A * w, 2L, pi_A, FUN = "-")
    share_term <- drop(G_cent %*% ((s_vec * theta_A - delta_hat) / Pi_A))
    phi_agg <- cell_term + share_term
    var_agg <- as.double(cov_from_phi(matrix(phi_agg, ncol = 1L), cluster)[1L, 1L])

    agg_res <- list(
      delta_hat = delta_hat,
      omega_hat = omega_A,
      Pi_hat = Pi_A,
      phi = phi_agg,
      cell_if_term = cell_term,
      share_if_term = share_term,
      variance = var_agg,
      vhat = var_agg * n_num,
      se = sqrt(var_agg)
    )
  }

  # Extract alpha_hat and block slopes from beta_hat
  alpha_hat <- stats::setNames(beta_hat[seq_len(n_cells)], cell_ids)
  slope_hat <- lapply(block_ids, function(bid) {
    p_b <- block_p[[bid]]
    off_b <- slope_offsets[[bid]]
    if (p_b > 0L) beta_hat[off_b + seq_len(p_b)] else numeric(0L)
  })
  names(slope_hat) <- block_ids

  list(
    theta_hat = theta_hat,
    tau_hat = tau_hat,
    pi_hat = pi_hat,
    Q_hat = Q_hat,
    mu_1 = mu_1,
    mu_0 = mu_0,
    gamma_hat = gamma_hat,
    H_hat = H_hat,
    beta_hat = beta_hat,
    alpha_hat = alpha_hat,
    slope_hat = slope_hat,
    M_hat = M_hat,
    ell_gamma = ell_gamma,
    ell_beta = ell_beta,
    phi_theta = phi_theta,
    phi_tau = phi_tau,
    cov_theta = cov_theta,
    cov_tau = cov_tau,
    aggregate = agg_res
  )
}

# Centered finite-difference directional derivative of the independent reference
# estimator under original-unit empirical weight perturbations:
#   w_+(i) = 1 + eps * (e_i - (1/n) 1_n)
#   w_-(i) = 1 - eps * (e_i - (1/n) 1_n)
#   FD_i = n * (T(w_+(i)) - T(w_-(i))) / (2 * eps)
ref_static_complete_fd_if <- function(cells,
                                      or_blocks,
                                      method,
                                      aggregate_cells = NULL,
                                      aggregate_signs = NULL,
                                      eps = 1e-6,
                                      glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE)) {
  cell_ids <- unname(vapply(cells, function(cl) cl$cell_id, character(1L)))
  n_cells <- length(cells)
  n <- length(cells[[1L]]$dy)
  n_num <- as.double(n)

  fd_theta <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  fd_tau <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  fd_agg <- if (!is.null(aggregate_cells)) numeric(n) else NULL

  base_w <- rep(1.0, n)
  for (i in seq_len(n)) {
    dir_i <- rep(-1.0 / n_num, n)
    dir_i[i] <- dir_i[i] + 1.0
    w_plus <- base_w + eps * dir_i
    w_minus <- base_w - eps * dir_i

    fit_plus <- ref_static_complete_eval(
      cells = cells,
      or_blocks = or_blocks,
      method = method,
      unit_weights = w_plus,
      aggregate_cells = aggregate_cells,
      aggregate_signs = aggregate_signs,
      glm_control = glm_control,
      compute_if = FALSE,
      polish = TRUE
    )
    fit_minus <- ref_static_complete_eval(
      cells = cells,
      or_blocks = or_blocks,
      method = method,
      unit_weights = w_minus,
      aggregate_cells = aggregate_cells,
      aggregate_signs = aggregate_signs,
      glm_control = glm_control,
      compute_if = FALSE,
      polish = TRUE
    )

    fd_theta[i, ] <- n_num * (fit_plus$theta_hat - fit_minus$theta_hat) / (2.0 * eps)
    fd_tau[i, ] <- n_num * (fit_plus$tau_hat - fit_minus$tau_hat) / (2.0 * eps)
    if (!is.null(aggregate_cells)) {
      fd_agg[i] <- n_num * (fit_plus$aggregate$delta_hat - fit_minus$aggregate$delta_hat) / (2.0 * eps)
    }
  }

  list(
    phi_theta = fd_theta,
    phi_tau = fd_tau,
    phi_agg = fd_agg
  )
}
