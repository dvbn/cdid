# ==============================================================================
# Independent Tiny Dense Linear-Algebra Reference for Static Role-Graph Joint GMM
# ==============================================================================
#
# Used ONLY by tests (`tests/testthat/test-static-graph-gmm.R`) to verify
# `solve_static_graph_gmm()` against direct textbook matrix algebra:
#   hat_nu    = (D' W D)^{-1} D' W mhat
#   hat_tau   = C' hat_nu = Lambda' mhat
#   Lambda    = W D (D' W D)^{-1} C
#   psi       = Phi %*% Lambda
#   V_sand    = normalization_factor * sum (psi_c)(psi_c)'
#
# Supports explicit formula mutations (`mutate_*`) to verify that tests fail
# when the adapter contract is violated.
# ==============================================================================

dense_reference_static_graph_gmm <- function(
  blocks,
  method = c("ipw", "dr", "or"),
  mhat,
  phi,
  cluster_ids = NULL,
  weighting = c("optimal", "two_step", "identity"),
  centered = TRUE,
  normalization_factor = 1.0 / nrow(phi),
  condition_trigger = 1e-10,
  ridge_relative = 1e-4,
  cell_ids = NULL,
  mutate_target_sign = FALSE,
  mutate_block_diagonal_V = FALSE,
  mutate_ground_first_date = FALSE,
  mutate_inverse_info_vcov = FALSE,
  gauge_rule = c("lowest", "highest", "custom")
) {
  method <- match.arg(method)
  weighting <- match.arg(weighting)
  gauge_rule <- match.arg(gauge_rule)
  w_canonical <- if (weighting == "two_step") "optimal" else weighting

  n_blocks <- length(blocks)
  block_cells <- vapply(blocks, `[[`, character(1L), "cell_id")
  if (is.null(cell_ids)) {
    cell_ids <- unique(block_cells)
  }

  # 1. Build dense grounded incidence matrix D (m_R x rank) directly from block
  #    date/edge tables using independent boolean adjacency closure per block
  #    (one coordinate gauge per connected component, with t-1 in the endpoint component).
  d_blocks <- vector("list", n_blocks)
  edge_block_idx <- integer(0L)
  link_ids <- character(0L)
  param_ids <- character(0L)
  target_param_col <- integer(n_blocks)
  col_offset <- 0L

  for (b in seq_len(n_blocks)) {
    blk <- blocks[[b]]
    dates <- as.double(blk$dates)
    nv_b <- length(dates)
    nl_b <- nrow(blk$edges)
    s_vec <- as.double(blk$edges$s)
    u_vec <- as.double(blk$edges$u)

    adj_b <- diag(TRUE, nv_b)
    for (e in seq_len(nl_b)) {
      si <- match(s_vec[e], dates)
      ui <- match(u_vec[e], dates)
      adj_b[si, ui] <- TRUE
      adj_b[ui, si] <- TRUE
    }
    repeat {
      adj_next <- (adj_b %*% adj_b) > 0
      if (identical(adj_next, adj_b)) break
      adj_b <- adj_next
    }

    base_idx <- match(as.double(blk$baseline_time), dates)
    assigned <- logical(nv_b)
    ground_dates_b <- numeric(0L)

    # Endpoint component containing baseline_time = t - 1
    comp1_idx <- which(adj_b[base_idx, ])
    assigned[comp1_idx] <- TRUE
    g1 <- if (isTRUE(mutate_ground_first_date)) min(dates[comp1_idx]) else as.double(blk$baseline_time)
    ground_dates_b <- c(ground_dates_b, g1)

    # Additional connected components (including any isolated non-target dates)
    custom_refs <- if (!is.null(blk$custom_reference_dates)) as.double(blk$custom_reference_dates) else NULL
    for (v_i in seq_len(nv_b)) {
      if (!assigned[v_i]) {
        comp_m_idx <- which(adj_b[v_i, ])
        assigned[comp_m_idx] <- TRUE
        comp_m_dates <- dates[comp_m_idx]
        gm <- if (!is.null(custom_refs) && any(comp_m_dates %in% custom_refs)) {
          comp_m_dates[comp_m_dates %in% custom_refs][1L]
        } else if (gauge_rule == "highest") {
          max(comp_m_dates)
        } else {
          min(comp_m_dates)
        }
        ground_dates_b <- c(ground_dates_b, gm)
      }
    }

    free_dates <- dates[!(dates %in% ground_dates_b)]
    np_b <- length(free_dates)

    D_b <- matrix(0.0, nrow = nl_b, ncol = np_b)
    for (e in seq_len(nl_b)) {
      s_col <- match(s_vec[e], free_dates)
      u_col <- match(u_vec[e], free_dates)
      if (!is.na(s_col)) D_b[e, s_col] <- D_b[e, s_col] - 1.0
      if (!is.na(u_col)) D_b[e, u_col] <- D_b[e, u_col] + 1.0
    }
    d_blocks[[b]] <- D_b
    edge_block_idx <- c(edge_block_idx, rep(b, nl_b))
    l_ids_b <- if ("link_id" %in% names(blk$edges)) {
      as.character(blk$edges$link_id)
    } else {
      paste0(blk$block_id, ":", blk$edges$s, "->", blk$edges$u)
    }
    link_ids <- c(link_ids, l_ids_b)
    param_ids <- c(param_ids, paste0("nu(block=", blk$block_id, ",cell=", blk$cell_id, ",role=", blk$role, ",t=", free_dates, ")"))
    target_param_col[b] <- col_offset + match(as.double(blk$target_time), free_dates)
    col_offset <- col_offset + np_b
  }

  m_R <- length(link_ids)
  p <- length(param_ids)
  D <- matrix(0.0, nrow = m_R, ncol = p, dimnames = list(link_ids, param_ids))
  r_off <- 0L
  c_off <- 0L
  for (b in seq_len(n_blocks)) {
    nr <- nrow(d_blocks[[b]])
    nc <- ncol(d_blocks[[b]])
    D[r_off + seq_len(nr), c_off + seq_len(nc)] <- d_blocks[[b]]
    r_off <- r_off + nr
    c_off <- c_off + nc
  }

  # 2. Build dense target contrast matrix C (p x J)
  J <- length(cell_ids)
  C_mat <- matrix(0.0, nrow = p, ncol = J, dimnames = list(param_ids, cell_ids))
  block_roles <- vapply(blocks, function(x) if (x$role %in% c("G", "switcher")) "G" else "C", character(1L))
  for (j in seq_len(J)) {
    cid <- cell_ids[j]
    b_g <- which(block_cells == cid & block_roles == "G")
    C_mat[target_param_col[b_g], j] <- 1.0
    if (method != "or") {
      b_c <- which(block_cells == cid & block_roles == "C")
      C_mat[target_param_col[b_c], j] <- if (isTRUE(mutate_target_sign)) 1.0 else -1.0
    }
  }

  # 3. Moment covariance V from complete-unit IF matrix phi (N x m_R)
  N <- nrow(phi)
  phi_use <- if (isTRUE(centered)) sweep(phi, 2L, colMeans(phi), "-") else phi
  V <- as.double(normalization_factor) * crossprod(phi_use)
  V <- 0.5 * (V + t(V))
  if (isTRUE(mutate_block_diagonal_V)) {
    for (r in seq_len(m_R)) {
      for (c_idx in seq_len(m_R)) {
        if (edge_block_idx[r] != edge_block_idx[c_idx]) {
          V[r, c_idx] <- 0.0
        }
      }
    }
  }

  # 4. Weight matrix W (with Tikhonov ridge when triggered under optimal/two_step)
  ridge_lambda <- 0.0
  handler_used <- "none"
  if (w_canonical == "identity") {
    W <- diag(1.0, m_R)
    V_reg <- W
  } else {
    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    lam_max <- max(ev)
    lam_min <- min(ev)
    rho <- if (lam_max > 0.0) lam_min / lam_max else -Inf
    if (lam_max <= 0.0 || lam_min <= 0.0 || rho < condition_trigger) {
      ridge_lambda <- ridge_relative * (sum(diag(V)) / as.double(m_R))
      V_reg <- V + diag(ridge_lambda, m_R)
      W <- solve(V_reg)
      handler_used <- "tikhonov"
    } else {
      V_reg <- V
      W <- solve(V)
    }
  }

  # 5. Whiten by V_reg = L L', then factor A P = Q R with column pivots.
  # W is diagnostic only: neither W D nor design normal equations enter the solve.
  covariance_chol <- chol(V_reg)
  A <- forwardsolve(t(covariance_chol), D)
  design_qr <- qr(A, LAPACK = TRUE)
  design_R <- qr.R(design_qr)
  pivot <- design_qr$pivot
  whitened_mhat <- forwardsolve(t(covariance_chol), as.double(mhat))
  qty <- qr.qty(design_qr, whitened_mhat)
  nu_hat <- numeric(p)
  nu_hat[pivot] <- as.double(backsolve(design_R, qty[seq_len(p)]))
  names(nu_hat) <- param_ids

  # Lambda = L^{-T} Q R^{-T} P' C; undo both information-matrix pivots.
  target_coords <- forwardsolve(t(design_R), C_mat[pivot, , drop = FALSE])
  Lambda <- backsolve(covariance_chol, qr.qy(
    design_qr, rbind(target_coords, matrix(0.0, m_R - p, J))
  ))
  DtWD_inv <- matrix(0.0, p, p, dimnames = list(param_ids, param_ids))
  DtWD_inv[pivot, pivot] <- chol2inv(design_R)
  dimnames(Lambda) <- list(link_ids, cell_ids)
  estimates <- as.double(crossprod(C_mat, nu_hat))
  names(estimates) <- cell_ids

  # 6. Full-unit target IF and empirical sandwich covariance
  target_if <- phi %*% Lambda
  dimnames(target_if) <- list(rownames(phi), cell_ids)

  if (isTRUE(mutate_inverse_info_vcov)) {
    sandwich_cov <- crossprod(C_mat, DtWD_inv %*% C_mat)
  } else if (is.null(cluster_ids)) {
    psi_c <- if (isTRUE(centered)) sweep(target_if, 2L, colMeans(target_if), "-") else target_if
    sandwich_cov <- as.double(normalization_factor) * crossprod(psi_c)
  } else {
    u_cl <- unique(as.character(cluster_ids))
    G <- length(u_cl)
    cl_sums <- matrix(0.0, nrow = G, ncol = J, dimnames = list(u_cl, cell_ids))
    for (g_i in seq_len(G)) {
      rows_g <- which(as.character(cluster_ids) == u_cl[g_i])
      cl_sums[g_i, ] <- colSums(target_if[rows_g, , drop = FALSE])
    }
    cl_c <- if (isTRUE(centered)) sweep(cl_sums, 2L, colMeans(cl_sums), "-") else cl_sums
    sandwich_cov <- as.double(normalization_factor) * crossprod(cl_c)
  }
  sandwich_cov <- 0.5 * (sandwich_cov + t(sandwich_cov))
  dimnames(sandwich_cov) <- list(cell_ids, cell_ids)

  raw_S <- sandwich_cov / as.double(normalization_factor)
  asymptotic_cov <- raw_S / as.double(N)
  vcov_n2 <- raw_S / (as.double(N)^2)

  list(
    D = D,
    C = C_mat,
    V = V,
    W = W,
    ridge_lambda = ridge_lambda,
    handler_used = handler_used,
    nu_hat = nu_hat,
    Lambda = Lambda,
    estimates = estimates,
    se = sqrt(diag(sandwich_cov)),
    se_n2 = sqrt(diag(vcov_n2)),
    sandwich_cov = sandwich_cov,
    asymptotic_cov = asymptotic_cov,
    vcov_n2 = vcov_n2,
    target_if = target_if
  )
}
