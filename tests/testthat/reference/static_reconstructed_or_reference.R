# Independent reference implementation for deterministic stayer graph
# reconstruction of source-adjacent outcome regression, its complete
# original-unit influence function, and fixed switch-step link centering.
#
# Calls NO production functions (`static_reconstructed_or_*`, `or_wls_*`,
# `or_pooled_wls_*`, `static_obs_*`, `static_complete_*`).

ref_recon_or_incidence <- function(edges,
                                   t,
                                   baseline_date = t - 1,
                                   dates = NULL,
                                   reference_dates = NULL) {
  edges <- as.matrix(edges)
  storage.mode(edges) <- "double"
  s <- edges[, 1L]
  u <- edges[, 2L]
  verts <- if (!is.null(dates)) sort(unique(as.double(dates))) else sort(unique(c(s, u)))
  n_v <- length(verts)
  K <- nrow(edges)

  adj <- diag(TRUE, n_v)
  for (k in seq_len(K)) {
    si <- match(s[k], verts)
    ui <- match(u[k], verts)
    adj[si, ui] <- TRUE
    adj[ui, si] <- TRUE
  }
  repeat {
    adj_next <- (adj %*% adj) > 0
    if (identical(adj_next, adj)) break
    adj <- adj_next
  }

  base_idx <- match(as.double(baseline_date), verts)
  assigned <- logical(n_v)
  comp1_idx <- which(adj[base_idx, ])
  assigned[comp1_idx] <- TRUE
  ref_dates <- as.double(baseline_date)
  ref_arg <- if (!is.null(reference_dates)) as.double(reference_dates) else NULL

  for (v_i in seq_len(n_v)) {
    if (!assigned[v_i]) {
      cm_idx <- which(adj[v_i, ])
      assigned[cm_idx] <- TRUE
      cm_dates <- verts[cm_idx]
      rm_val <- if (!is.null(ref_arg) && any(cm_dates %in% ref_arg)) {
        cm_dates[cm_dates %in% ref_arg][1L]
      } else {
        min(cm_dates)
      }
      ref_dates <- c(ref_dates, rm_val)
    }
  }

  non_ground <- verts[!(verts %in% ref_dates)]
  V <- length(non_ground)
  D <- matrix(0.0, nrow = K, ncol = V,
              dimnames = list(paste0(s, "->", u), as.character(non_ground)))
  for (k in seq_len(K)) {
    if (!(s[k] %in% ref_dates)) D[k, match(s[k], non_ground)] <- -1.0
    if (!(u[k] %in% ref_dates)) D[k, match(u[k], non_ground)] <- 1.0
  }
  ell <- numeric(V)
  names(ell) <- as.character(non_ground)
  ell[match(t, non_ground)] <- 1.0
  u_pot <- drop(solve(crossprod(D), ell))
  # Enforce exact zero on columns outside the endpoint component
  comp1_non_ground <- verts[comp1_idx][verts[comp1_idx] != as.double(baseline_date)]
  u_pot[!(non_ground %in% comp1_non_ground)] <- 0.0
  varpi <- drop(D %*% u_pot)
  names(varpi) <- rownames(D)
  list(
    D = D,
    ell = ell,
    varpi = varpi,
    vertices = verts,
    non_ground = non_ground,
    reference_dates = ref_dates
  )
}

# Independent weighted edge mean and covariance evaluator under unit weights w.
# Supports routes: "baseline", "weighting", "imputation", "augmented".
ref_recon_or_edge_moment <- function(kind = c("mean", "cov"),
                                     A,
                                     route,
                                     O = NULL,
                                     dy = NULL,
                                     q = NULL,
                                     j = NULL,
                                     h = NULL,
                                     Z = NULL,
                                     unit_weights = NULL,
                                     q_jac = list(),
                                     j_jac = list(),
                                     h_jac = list(),
                                     coef_if = list()) {
  kind <- match.arg(kind)
  n <- length(A)
  n_num <- as.double(n)
  w_u <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  A_num <- as.double(A)
  eps_floor <- 1e-6

  if (kind == "mean") {
    h_vec <- as.double(h)
  } else {
    h_vec <- rep(1.0, n)
  }

  O_num <- if (!is.null(O)) as.double(O) else rep(0.0, n)
  phi_vec <- numeric(n)
  dphi_dq <- numeric(n)
  if (route %in% c("weighting", "augmented")) {
    idx_o <- which(A_num == 1 & O_num == 1)
    q_o <- as.double(q[idx_o])
    act_o <- as.double(q_o > eps_floor)
    qf_o <- pmax(q_o, eps_floor)
    phi_vec[idx_o] <- 1.0 / qf_o
    dphi_dq[idx_o] <- -act_o / (qf_o * qf_o)
  }

  w_route <- numeric(n)
  Y_pseudo <- numeric(n)
  if (route == "baseline") {
    idx <- which(A_num == 1 & O_num == 1)
    w_route[idx] <- w_u[idx] * h_vec[idx]
    Y_pseudo[idx] <- as.double(dy[idx])
  } else if (route == "weighting") {
    idx <- which(A_num == 1 & O_num == 1)
    w_route[idx] <- w_u[idx] * h_vec[idx] * phi_vec[idx]
    Y_pseudo[idx] <- as.double(dy[idx])
  } else if (route == "imputation") {
    idx <- which(A_num == 1)
    w_route[idx] <- w_u[idx] * h_vec[idx]
    Y_pseudo[idx] <- as.double(j[idx])
  } else if (route == "augmented") {
    idx <- which(A_num == 1)
    w_route[idx] <- w_u[idx] * h_vec[idx]
    j_idx <- as.double(j[idx])
    corr <- numeric(length(idx))
    obs_local <- which(O_num[idx] == 1)
    if (length(obs_local) > 0L) {
      obs_rows <- idx[obs_local]
      corr[obs_local] <- phi_vec[obs_rows] * (as.double(dy[obs_rows]) - j_idx[obs_local])
    }
    Y_pseudo[idx] <- j_idx + corr
  }

  S <- sum(w_route)
  M <- sum(w_route * Y_pseudo) / S

  if (kind == "mean") {
    est <- M
    Fw <- matrix(Y_pseudo - M, ncol = 1L)
    Gy <- matrix(1.0, nrow = n, ncol = 1L)
    direct <- n_num * w_route * (Y_pseudo - M) / S
  } else {
    Zm <- as.matrix(Z)
    storage.mode(Zm) <- "double"
    p <- ncol(Zm)
    Zbar <- colSums(Zm * w_route) / S
    Zc <- sweep(Zm, 2L, Zbar, FUN = "-")
    est <- drop(crossprod(Zc, w_route * Y_pseudo)) / S
    names(est) <- colnames(Zm)
    Fw <- Zc * (Y_pseudo - M) - matrix(est, nrow = n, ncol = p, byrow = TRUE)
    Gy <- Zc
    direct <- matrix(0.0, nrow = n, ncol = p, dimnames = list(rownames(Zm), colnames(Zm)))
    direct[ idx, ] <- (n_num * w_route[idx] / S) * Fw[idx, , drop = FALSE]
  }

  # Nuisance loadings
  loadings <- list()
  add_loading <- function(load_list, key, val) {
    if (is.null(load_list[[key]])) load_list[[key]] <- val else load_list[[key]] <- load_list[[key]] + val
    load_list
  }

  if (kind == "mean" && length(h_jac) > 0L) {
    dw_dh <- if (route == "weighting") A_num * O_num * phi_vec else if (route == "baseline") A_num * O_num else A_num
    coef_h <- (w_u * dw_dh) * Fw
    for (k in names(h_jac)) {
      Jk <- drop(crossprod(h_jac[[k]][idx, , drop = FALSE], coef_h[idx, , drop = FALSE])) / S
      loadings <- add_loading(loadings, k, Jk)
    }
  }
  if (route %in% c("weighting", "augmented") && length(q_jac) > 0L) {
    idx_o <- which(A_num == 1 & O_num == 1)
    if (route == "weighting") {
      coef_q <- (w_u * A_num * O_num * h_vec * dphi_dq) * Fw
    } else {
      resid_o <- numeric(n)
      resid_o[idx_o] <- as.double(dy[idx_o]) - as.double(j[idx_o])
      coef_q <- (w_u * A_num * O_num * h_vec * resid_o * dphi_dq) * Gy
    }
    for (k in names(q_jac)) {
      val <- crossprod(q_jac[[k]][idx_o, , drop = FALSE], coef_q[idx_o, , drop = FALSE]) / S
      Jk <- if (kind == "mean") drop(val[, 1L]) else t(val)
      loadings <- add_loading(loadings, k, Jk)
    }
  }
  if (route %in% c("imputation", "augmented") && length(j_jac) > 0L) {
    mult_j <- if (route == "imputation") rep(1.0, n) else (1.0 - O_num * phi_vec)
    coef_j <- (w_u * A_num * h_vec * mult_j) * Gy
    for (k in names(j_jac)) {
      val <- crossprod(j_jac[[k]][idx, , drop = FALSE], coef_j[idx, , drop = FALSE]) / S
      Jk <- if (kind == "mean") drop(val[, 1L]) else t(val)
      loadings <- add_loading(loadings, k, Jk)
    }
  }

  inf <- direct
  if (length(loadings) > 0L && length(coef_if) > 0L) {
    for (k in names(loadings)) {
      Lk <- coef_if[[k]]
      Jk <- loadings[[k]]
      if (kind == "mean") {
        inf <- inf + drop(Lk %*% Jk)
      } else {
        inf <- inf + Lk %*% t(Jk)
      }
    }
  }

  list(
    kind = kind,
    route = route,
    estimate = est,
    direct_if = direct,
    loadings = loadings,
    influence = inf
  )
}

# Independent evaluation of full stacked normal equations, per-block reconstructed
# OR, complete original-unit IF pieces, and fixed switch-step link centering.
ref_reconstructed_or_eval <- function(cells,
                                      or_blocks,
                                      unit_weights = NULL,
                                      compute_if = TRUE) {
  cell_ids <- names(cells)
  if (is.null(cell_ids)) cell_ids <- vapply(cells, `[[`, character(1L), "cell_id")
  names(cells) <- cell_ids
  n_cells <- length(cell_ids)
  n <- length(cells[[1L]]$C)
  n_num <- as.double(n)
  w_u <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  W_sum <- sum(w_u)

  # Parse partition
  if (is.list(or_blocks)) {
    block_ids <- names(or_blocks)
    if (is.null(block_ids)) block_ids <- paste0("block_", seq_along(or_blocks))
    cell_block <- stats::setNames(character(n_cells), cell_ids)
    for (bi in seq_along(or_blocks)) {
      mems <- or_blocks[[bi]]
      if (is.numeric(mems)) mems <- cell_ids[as.integer(mems)]
      cell_block[mems] <- block_ids[bi]
    }
  } else {
    ob_chr <- as.character(or_blocks)
    if (!is.null(names(or_blocks))) ob_chr <- ob_chr[match(cell_ids, names(or_blocks))]
    cell_block <- stats::setNames(ob_chr, cell_ids)
    block_ids <- unique(ob_chr)
  }
  blocks <- lapply(block_ids, function(bid) cell_ids[cell_block == bid])
  names(blocks) <- block_ids

  pi_hat <- stats::setNames(numeric(n_cells), cell_ids)
  zbar_hat <- vector("list", n_cells)
  names(zbar_hat) <- cell_ids
  S_hat <- vector("list", n_cells)
  names(S_hat) <- cell_ids
  varpi_hat <- vector("list", n_cells)
  names(varpi_hat) <- cell_ids
  mC_hat <- stats::setNames(numeric(n_cells), cell_ids)
  Cov_hat <- vector("list", n_cells)
  names(Cov_hat) <- cell_ids

  for (cid in cell_ids) {
    cl <- cells[[cid]]
    C_num <- as.double(cl$C)
    w_C <- w_u * C_num
    S_C <- sum(w_C)
    pi_c <- S_C / W_sum
    pi_hat[[cid]] <- pi_c

    Zm <- as.matrix(cl$Z)
    storage.mode(Zm) <- "double"
    p_c <- ncol(Zm)
    z_names <- colnames(Zm)
    if (is.null(z_names)) z_names <- if (p_c > 0L) paste0("z", seq_len(p_c)) else character(0L)
    st_idx <- which(C_num == 1)
    if (p_c > 0L) {
      Z_st <- Zm[st_idx, , drop = FALSE]
      w_st <- w_C[st_idx]
      zbar_c <- stats::setNames(colSums(Z_st * w_st) / S_C, z_names)
      Zc_st <- sweep(Z_st, 2L, zbar_c, FUN = "-")
      S_c <- crossprod(Zc_st, Zc_st * w_st) / S_C
      dimnames(S_c) <- list(z_names, z_names)
    } else {
      zbar_c <- stats::setNames(numeric(0L), character(0L))
      S_c <- matrix(0.0, 0L, 0L)
    }
    zbar_hat[[cid]] <- zbar_c
    S_hat[[cid]] <- S_c

    if (!is.null(cl$flow)) {
      vp <- as.double(cl$flow)
    } else if (!is.null(cl$D) && !is.null(cl$ell)) {
      D_c <- as.matrix(cl$D)
      ell_c <- as.double(cl$ell)
      vp <- drop(D_c %*% solve(crossprod(D_c), ell_c))
    } else {
      inc <- ref_recon_or_incidence(
        edges = cl$edges,
        t = cl$t,
        baseline_date = if (!is.null(cl$baseline_date)) cl$baseline_date else cl$t - 1,
        dates = cl$dates,
        reference_dates = cl$reference_dates
      )
      vp <- inc$varpi
    }
    varpi_hat[[cid]] <- vp

    em_vec <- vapply(cl$edge_means, `[[`, numeric(1L), "estimate")
    mC_hat[[cid]] <- sum(vp * em_vec)
    if (p_c > 0L) {
      ec_mat <- do.call(rbind, lapply(cl$edge_covs, `[[`, "estimate"))
      Cov_c <- stats::setNames(drop(crossprod(vp, ec_mat)), z_names)
    } else {
      Cov_c <- stats::setNames(numeric(0L), character(0L))
    }
    Cov_hat[[cid]] <- Cov_c
  }

  # Block-wise solve AND full stacked normal equations solve
  A_hat <- vector("list", length(block_ids))
  names(A_hat) <- block_ids
  B_hat <- vector("list", length(block_ids))
  names(B_hat) <- block_ids
  b_hat <- vector("list", length(block_ids))
  names(b_hat) <- block_ids
  alpha_hat <- stats::setNames(numeric(n_cells), cell_ids)

  block_p <- integer(length(block_ids))
  names(block_p) <- block_ids
  for (bid in block_ids) {
    cids <- blocks[[bid]]
    p_b <- ncol(as.matrix(cells[[cids[1L]]]$Z))
    block_p[[bid]] <- p_b
    zn_b <- names(zbar_hat[[cids[1L]]])
    if (p_b > 0L) {
      A_P <- matrix(0.0, p_b, p_b, dimnames = list(zn_b, zn_b))
      B_P <- stats::setNames(numeric(p_b), zn_b)
      for (cid in cids) {
        A_P <- A_P + pi_hat[[cid]] * S_hat[[cid]]
        B_P <- B_P + pi_hat[[cid]] * Cov_hat[[cid]]
      }
      b_P <- stats::setNames(drop(solve(A_P, B_P)), zn_b)
    } else {
      A_P <- matrix(0.0, 0L, 0L)
      B_P <- numeric(0L)
      b_P <- numeric(0L)
    }
    A_hat[[bid]] <- A_P
    B_hat[[bid]] <- B_P
    b_hat[[bid]] <- b_P
    for (cid in cids) {
      alpha_hat[[cid]] <- if (p_b > 0L) mC_hat[[cid]] - sum(zbar_hat[[cid]] * b_P) else mC_hat[[cid]]
    }
  }

  # Full stacked normal system M_full * beta_full = rhs_full:
  # E_n[sum_{c in P} C_c a_c (Delta_c Y - a_c' beta)] = 0 with E_n[C_c Delta_c Y]
  # replaced by pi_c mC_c and E_n[C_c Z_c Delta_c Y] replaced by pi_c (Cov_c + zbar_c mC_c).
  slope_offsets <- integer(length(block_ids))
  names(slope_offsets) <- block_ids
  cur_off <- n_cells
  for (bid in block_ids) {
    slope_offsets[[bid]] <- cur_off
    cur_off <- cur_off + block_p[[bid]]
  }
  dim_beta <- cur_off
  M_full <- matrix(0.0, dim_beta, dim_beta)
  rhs_full <- numeric(dim_beta)
  for (j in seq_len(n_cells)) {
    cid <- cell_ids[j]
    cl <- cells[[cid]]
    bid <- cell_block[[cid]]
    p_b <- block_p[[bid]]
    off_b <- slope_offsets[[bid]]
    C_num <- as.double(cl$C)
    st_idx <- which(C_num == 1)
    w_st <- w_u[st_idx]
    A_st <- matrix(0.0, nrow = length(st_idx), ncol = dim_beta)
    A_st[, j] <- 1.0
    if (p_b > 0L) {
      A_st[, off_b + seq_len(p_b)] <- as.matrix(cl$Z)[st_idx, , drop = FALSE]
    }
    M_full <- M_full + crossprod(A_st, A_st * w_st) / W_sum
    rhs_full[j] <- pi_hat[[cid]] * mC_hat[[cid]]
    if (p_b > 0L) {
      idx_b <- off_b + seq_len(p_b)
      rhs_full[idx_b] <- rhs_full[idx_b] +
        pi_hat[[cid]] * (Cov_hat[[cid]] + zbar_hat[[cid]] * mC_hat[[cid]])
    }
  }
  beta_full <- drop(solve(M_full, rhs_full))
  alpha_full <- stats::setNames(beta_full[seq_len(n_cells)], cell_ids)
  b_full <- lapply(block_ids, function(bid) {
    p_b <- block_p[[bid]]
    if (p_b == 0L) return(numeric(0L))
    off_b <- slope_offsets[[bid]]
    stats::setNames(beta_full[off_b + seq_len(p_b)], names(b_hat[[bid]]))
  })
  names(b_full) <- block_ids

  if (!compute_if) {
    return(list(
      pi_hat = pi_hat,
      zbar_hat = zbar_hat,
      S_hat = S_hat,
      varpi_hat = varpi_hat,
      mC_hat = mC_hat,
      Cov_hat = Cov_hat,
      A_hat = A_hat,
      B_hat = B_hat,
      b_hat = b_hat,
      alpha_hat = alpha_hat,
      M_full = M_full,
      rhs_full = rhs_full,
      alpha_full = alpha_full,
      b_full = b_full
    ))
  }

  # Influence functions by direct unit-loop over i = 1..n
  chi_pi <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  chi_mC <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  chi_alpha <- matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  chi_z <- vector("list", n_cells)
  names(chi_z) <- cell_ids
  chi_Cov <- vector("list", n_cells)
  names(chi_Cov) <- cell_ids
  chi_A_cell <- vector("list", n_cells)
  names(chi_A_cell) <- cell_ids
  chi_B_cell <- vector("list", n_cells)
  names(chi_B_cell) <- cell_ids
  chi_A_block <- vector("list", length(block_ids))
  names(chi_A_block) <- block_ids
  chi_A_b_block <- vector("list", length(block_ids))
  names(chi_A_b_block) <- block_ids
  chi_B_block <- vector("list", length(block_ids))
  names(chi_B_block) <- block_ids
  chi_b <- vector("list", length(block_ids))
  names(chi_b) <- block_ids

  for (bid in block_ids) {
    cids <- blocks[[bid]]
    p_b <- block_p[[bid]]
    zn_b <- names(b_hat[[bid]])
    b_P <- b_hat[[bid]]
    A_P <- A_hat[[bid]]

    chi_A_P <- array(0.0, dim = c(n, p_b, p_b), dimnames = list(NULL, zn_b, zn_b))
    chi_A_b_P <- matrix(0.0, nrow = n, ncol = p_b, dimnames = list(NULL, zn_b))
    chi_B_P <- matrix(0.0, nrow = n, ncol = p_b, dimnames = list(NULL, zn_b))

    for (cid in cids) {
      cl <- cells[[cid]]
      C_num <- as.double(cl$C)
      pi_c <- pi_hat[[cid]]
      zbar_c <- zbar_hat[[cid]]
      S_c <- S_hat[[cid]]
      Cov_c <- Cov_hat[[cid]]
      vp <- varpi_hat[[cid]]
      K_c <- length(vp)

      chi_pi[, cid] <- C_num - pi_c
      m_if <- numeric(n)
      for (k in seq_len(K_c)) {
        em_k <- cl$edge_means[[k]]
        if_k <- if (!is.null(em_k[["influence"]])) em_k[["influence"]] else em_k[["eval_if"]](seq_len(n))
        m_if <- m_if + vp[k] * as.double(if_k)
      }
      chi_mC[, cid] <- m_if

      if (p_b > 0L) {
        c_if <- matrix(0.0, nrow = n, ncol = p_b, dimnames = list(NULL, zn_b))
        for (k in seq_len(K_c)) {
          ec_k <- cl$edge_covs[[k]]
          if_ck <- if (!is.null(ec_k[["influence"]])) ec_k[["influence"]] else ec_k[["eval_if"]](seq_len(n))
          c_if <- c_if + vp[k] * as.matrix(if_ck)
        }
        chi_Cov[[cid]] <- c_if

        Zm <- as.matrix(cl$Z)
        chi_z_c <- matrix(0.0, nrow = n, ncol = p_b, dimnames = list(NULL, zn_b))
        chi_A_c <- array(0.0, dim = c(n, p_b, p_b), dimnames = list(NULL, zn_b, zn_b))
        chi_B_c <- matrix(0.0, nrow = n, ncol = p_b, dimnames = list(NULL, zn_b))
        for (i in seq_len(n)) {
          if (C_num[i] == 1) {
            dz_i <- as.double(Zm[i, ]) - zbar_c
            chi_z_c[i, ] <- dz_i / pi_c
            chi_A_c[i, , ] <- tcrossprod(dz_i) - pi_c * S_c
          } else {
            chi_A_c[i, , ] <- -pi_c * S_c
          }
          chi_B_c[i, ] <- (C_num[i] - pi_c) * Cov_c + pi_c * c_if[i, ]
          chi_A_b_P[i, ] <- chi_A_b_P[i, ] + drop(chi_A_c[i, , ] %*% b_P)
        }
        chi_z[[cid]] <- chi_z_c
        chi_A_cell[[cid]] <- chi_A_c
        chi_B_cell[[cid]] <- chi_B_c
        chi_A_P <- chi_A_P + chi_A_c
        chi_B_P <- chi_B_P + chi_B_c
      } else {
        chi_Cov[[cid]] <- matrix(0.0, nrow = n, ncol = 0L)
        chi_z[[cid]] <- matrix(0.0, nrow = n, ncol = 0L)
        chi_A_cell[[cid]] <- array(0.0, dim = c(n, 0L, 0L))
        chi_B_cell[[cid]] <- matrix(0.0, nrow = n, ncol = 0L)
      }
    }

    chi_b_P <- if (p_b > 0L) {
      sol <- t(solve(A_P, t(chi_B_P - chi_A_b_P)))
      if (p_b == 1L) sol <- matrix(sol, nrow = n, ncol = 1L)
      colnames(sol) <- zn_b
      sol
    } else {
      matrix(0.0, nrow = n, ncol = 0L)
    }

    chi_A_block[[bid]] <- chi_A_P
    chi_A_b_block[[bid]] <- chi_A_b_P
    chi_B_block[[bid]] <- chi_B_P
    chi_b[[bid]] <- chi_b_P

    for (cid in cids) {
      if (p_b > 0L) {
        chi_alpha[, cid] <- chi_mC[, cid] -
          drop(chi_z[[cid]] %*% b_P) -
          drop(chi_b_P %*% zbar_hat[[cid]])
      } else {
        chi_alpha[, cid] <- chi_mC[, cid]
      }
    }
  }

  list(
    pi_hat = pi_hat,
    zbar_hat = zbar_hat,
    S_hat = S_hat,
    varpi_hat = varpi_hat,
    mC_hat = mC_hat,
    Cov_hat = Cov_hat,
    A_hat = A_hat,
    B_hat = B_hat,
    b_hat = b_hat,
    alpha_hat = alpha_hat,
    M_full = M_full,
    rhs_full = rhs_full,
    alpha_full = alpha_full,
    b_full = b_full,
    chi_pi = chi_pi,
    chi_z = chi_z,
    chi_mC = chi_mC,
    chi_Cov = chi_Cov,
    chi_A_cell = chi_A_cell,
    chi_B_cell = chi_B_cell,
    chi_A = chi_A_block,
    chi_A_b = chi_A_b_block,
    chi_B = chi_B_block,
    chi_b = chi_b,
    chi_alpha = chi_alpha
  )
}

# Independent weighted logit fitter (polished Newton-Raphson) and IF evaluator
ref_weighted_logit_fit <- function(A, y, X, unit_weights = NULL) {
  n <- length(A)
  n_num <- as.double(n)
  w_u <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  W_sum <- sum(w_u)
  A_num <- as.double(A)
  idx <- which(A_num == 1)
  X_idx <- as.matrix(X)[idx, , drop = FALSE]
  y_idx <- as.double(y[idx])
  w_idx <- w_u[idx]

  fit0 <- suppressWarnings(stats::glm.fit(
    x = X_idx, y = y_idx, weights = w_idx,
    family = stats::binomial(link = "logit"),
    control = list(epsilon = 1e-13, maxit = 100L, trace = FALSE),
    intercept = FALSE
  ))
  beta <- as.double(fit0$coefficients)
  for (it in 1:6) {
    u <- drop(X_idx %*% beta)
    p <- 1.0 / (1.0 + exp(-u))
    g <- colSums(X_idx * (w_idx * (y_idx - p)))
    if (max(abs(g)) <= 1e-15) break
    H <- crossprod(X_idx, X_idx * (w_idx * p * (1.0 - p)))
    beta <- beta + drop(solve(H, g))
  }
  names(beta) <- colnames(X)

  u_all <- numeric(n)
  p_all <- numeric(n)
  u_all[idx] <- drop(X_idx %*% beta)
  p_all[idx] <- 1.0 / (1.0 + exp(-u_all[idx]))

  H_bar <- crossprod(X_idx, X_idx * (w_idx * p_all[idx] * (1.0 - p_all[idx]))) / W_sum
  score_mat <- matrix(0.0, nrow = n, ncol = ncol(X), dimnames = list(rownames(X), colnames(X)))
  score_mat[idx, ] <- X_idx * (w_idx * (y_idx - p_all[idx]))
  coef_if <- t(solve(H_bar, t(score_mat)))
  colnames(coef_if) <- colnames(X)

  jac <- matrix(0.0, nrow = n, ncol = ncol(X), dimnames = list(rownames(X), colnames(X)))
  jac[idx, ] <- X_idx * (p_all[idx] * (1.0 - p_all[idx]))

  eps_clip <- 1e-6
  cap_prob <- 1.0 - eps_clip
  cap_odds <- cap_prob / (1.0 - cap_prob)
  act_all <- numeric(n)
  odds_all <- numeric(n)
  act_all[idx] <- as.double(p_all[idx] < cap_prob)
  odds_all[idx] <- ifelse(act_all[idx] == 1, exp(u_all[idx]), cap_odds)
  odds_jac <- matrix(0.0, nrow = n, ncol = ncol(X), dimnames = list(rownames(X), colnames(X)))
  odds_jac[idx, ] <- X_idx * (act_all[idx] * odds_all[idx])

  list(
    beta = beta,
    prob = p_all,
    prob_jac = jac,
    odds = odds_all,
    active = act_all,
    odds_jac = odds_jac,
    coef_if = coef_if
  )
}

# Independent weighted OLS imputation fitter and IF evaluator
ref_weighted_imputation_fit <- function(A, O, dy, X_j, unit_weights = NULL) {
  n <- length(A)
  w_u <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  W_sum <- sum(w_u)
  A_num <- as.double(A)
  O_num <- as.double(O)
  idx_obs <- which(A_num == 1 & O_num == 1)
  idx_role <- which(A_num == 1)
  X_obs <- as.matrix(X_j)[idx_obs, , drop = FALSE]
  y_obs <- as.double(dy[idx_obs])
  w_obs <- w_u[idx_obs]

  G_bar <- crossprod(X_obs, X_obs * w_obs) / W_sum
  rhs <- colSums(X_obs * (w_obs * y_obs)) / W_sum
  theta <- drop(solve(G_bar, rhs))
  names(theta) <- colnames(X_j)

  j_pred <- numeric(n)
  j_pred[idx_role] <- drop(as.matrix(X_j)[idx_role, , drop = FALSE] %*% theta)
  j_jac <- matrix(0.0, nrow = n, ncol = ncol(X_j), dimnames = list(rownames(X_j), colnames(X_j)))
  j_jac[idx_role, ] <- as.matrix(X_j)[idx_role, , drop = FALSE]

  score_mat <- matrix(0.0, nrow = n, ncol = ncol(X_j), dimnames = list(rownames(X_j), colnames(X_j)))
  res_obs <- y_obs - drop(X_obs %*% theta)
  score_mat[idx_obs, ] <- X_obs * (w_obs * res_obs)
  coef_if <- t(solve(G_bar, t(score_mat)))
  colnames(coef_if) <- colnames(X_j)

  list(
    theta = theta,
    j = j_pred,
    j_jac = j_jac,
    coef_if = coef_if
  )
}

# Independent weighted role mean Z evaluator
ref_role_mean_z <- function(A, Z, r = NULL, X_tilde = NULL, active = NULL,
                            ell_gamma = NULL, unit_weights = NULL) {
  n <- length(A)
  w_u <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  W_sum <- sum(w_u)
  A_num <- as.double(A)
  idx <- which(A_num == 1)
  r_vec <- if (is.null(r)) rep(1.0, n) else as.double(r)
  w_eff <- numeric(n)
  w_eff[idx] <- w_u[idx] * r_vec[idx]
  S_w <- sum(w_eff)
  Q <- S_w / W_sum

  Zm <- as.matrix(Z)
  storage.mode(Zm) <- "double"
  p <- ncol(Zm)
  if (p == 0L) {
    return(list(estimate = numeric(0L), influence = matrix(0.0, nrow = n, ncol = 0L)))
  }
  zbar <- stats::setNames(colSums(Zm[idx, , drop = FALSE] * w_eff[idx]) / S_w, colnames(Zm))
  Zc <- sweep(Zm, 2L, zbar, FUN = "-")
  inf <- matrix(0.0, nrow = n, ncol = p, dimnames = list(rownames(Zm), colnames(Zm)))
  inf[idx, ] <- (Zc[idx, , drop = FALSE] * w_eff[idx]) / Q
  if (!is.null(ell_gamma)) {
    act_idx <- as.double(active[idx])
    X_idx <- as.matrix(X_tilde)[idx, , drop = FALSE]
    J_g <- crossprod(Zc[idx, , drop = FALSE], X_idx * (w_eff[idx] * act_idx)) / S_w
    inf <- inf + ell_gamma %*% t(J_g)
  }
  list(estimate = zbar, influence = inf)
}

# Independent fixed switch-step link centering reference
ref_center_cell_links <- function(ref_or,
                                  cell_id,
                                  block_id,
                                  method,
                                  t,
                                  G_edges,
                                  G_links,
                                  C_edges = NULL,
                                  C_links = NULL,
                                  zG = NULL,
                                  zCr = NULL,
                                  date_shift_G = NULL,
                                  date_shift_C = NULL) {
  pm <- if (identical(method, "reg")) "or" else method
  alpha_c <- ref_or$alpha_hat[[cell_id]]
  b_P <- ref_or$b_hat[[block_id]]
  p <- length(b_P)
  n <- nrow(ref_or$chi_alpha)

  G_edges <- as.matrix(G_edges)
  step_G <- as.double(G_edges[, 2L] >= t) - as.double(G_edges[, 1L] >= t)
  raw_G <- vapply(G_links, `[[`, numeric(1L), "estimate")
  chi_raw_G <- do.call(cbind, lapply(G_links, `[[`, "influence"))

  eval_shift <- function(edges, ds) {
    if (is.null(ds)) return(list(est = numeric(nrow(edges)), inf = matrix(0.0, n, nrow(edges))))
    u_pos <- match(edges[, 2L], ds$dates, nomatch = 0L)
    s_pos <- match(edges[, 1L], ds$dates, nomatch = 0L)
    u_val <- ifelse(u_pos > 0L, ds$values[u_pos], 0.0)
    s_val <- ifelse(s_pos > 0L, ds$values[s_pos], 0.0)
    inf <- matrix(0.0, n, nrow(edges))
    if (!is.null(ds$if_mat)) {
      for (k in seq_len(nrow(edges))) {
        if (u_pos[k] > 0L) inf[, k] <- inf[, k] + ds$if_mat[, u_pos[k]]
        if (s_pos[k] > 0L) inf[, k] <- inf[, k] - ds$if_mat[, s_pos[k]]
      }
    }
    list(est = u_val - s_val, inf = inf)
  }
  sh_G <- eval_shift(G_edges, date_shift_G)

  if (pm == "ipw") {
    C_edges <- as.matrix(C_edges)
    step_C <- as.double(C_edges[, 2L] >= t) - as.double(C_edges[, 1L] >= t)
    raw_C <- vapply(C_links, `[[`, numeric(1L), "estimate")
    chi_raw_C <- do.call(cbind, lapply(C_links, `[[`, "influence"))
    return(list(
      center_G = 0.0,
      center_C = 0.0,
      adjusted_G = raw_G,
      adjusted_C = raw_C,
      chi_center_G = numeric(n),
      chi_center_C = numeric(n),
      chi_adjusted_G = chi_raw_G,
      chi_adjusted_C = chi_raw_C
    ))
  }

  chi_alpha_c <- ref_or$chi_alpha[, cell_id]
  chi_b_P <- ref_or$chi_b[[block_id]]
  center_G <- if (p > 0L) alpha_c + sum(zG$estimate * b_P) else alpha_c
  chi_center_G <- if (p > 0L) {
    chi_alpha_c + drop(chi_b_P %*% zG$estimate) + drop(zG$influence %*% b_P)
  } else {
    chi_alpha_c
  }
  adj_G <- raw_G - step_G * center_G - sh_G$est
  chi_adj_G <- chi_raw_G - outer(chi_center_G, step_G) - sh_G$inf

  if (pm == "or") {
    return(list(
      center_G = center_G,
      center_C = NA_real_,
      adjusted_G = adj_G,
      adjusted_C = NULL,
      chi_center_G = chi_center_G,
      chi_center_C = NULL,
      chi_adjusted_G = chi_adj_G,
      chi_adjusted_C = NULL
    ))
  }

  C_edges <- as.matrix(C_edges)
  step_C <- as.double(C_edges[, 2L] >= t) - as.double(C_edges[, 1L] >= t)
  raw_C <- vapply(C_links, `[[`, numeric(1L), "estimate")
  chi_raw_C <- do.call(cbind, lapply(C_links, `[[`, "influence"))
  sh_C <- eval_shift(C_edges, date_shift_C)
  center_C <- if (p > 0L) alpha_c + sum(zCr$estimate * b_P) else alpha_c
  chi_center_C <- if (p > 0L) {
    chi_alpha_c + drop(chi_b_P %*% zCr$estimate) + drop(zCr$influence %*% b_P)
  } else {
    chi_alpha_c
  }
  adj_C <- raw_C - step_C * center_C - sh_C$est
  chi_adj_C <- chi_raw_C - outer(chi_center_C, step_C) - sh_C$inf

  list(
    center_G = center_G,
    center_C = center_C,
    adjusted_G = adj_G,
    adjusted_C = adj_C,
    chi_center_G = chi_center_G,
    chi_center_C = chi_center_C,
    chi_adjusted_G = chi_adj_G,
    chi_adjusted_C = chi_adj_C
  )
}

