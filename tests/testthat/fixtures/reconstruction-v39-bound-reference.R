# Exact pre-V39 resource reference from6ecfb78; no generated expected values.
# Only the two bound functions are copied; raw fixtures stay in existing files.
imputation_retained_core_bound <- function(cat, recipe, trt_recipe, dense_adjoints = FALSE) {
  E <- function(x) unname(cp04_payload(x)["entries"])
  # A names attribute costs 9 + E(names). A matrix costs its entries + 8,
  # the one-slot attributes list (9), and its two-element dim vector (10).
  A <- function(n, nm = NULL) n + 8 + if (is.null(nm)) 0 else 9 + E(nm)
  M <- function(nr, nc, rn = NULL, cn = NULL, dn = FALSE) {
    nr * nc + 27 + if (dn) 1 + E(list(rn, cn)) else 0
  }
  R <- function(fields, cls = NULL, nm = names(fields)) {
    length(fields) + 8 + sum(fields) +
      if (is.null(nm) && is.null(cls)) 0 else
        8 + as.integer(!is.null(nm)) + as.integer(!is.null(cls)) + E(nm) + E(cls)
  }
  V <- function(x) length(x) + 8 + sum(x)
  p <- as.double(recipe$p)
  z <- recipe$slope_columns
  q <- as.double(trt_recipe$p)
  nc <- length(cat$cell_specs)
  ng <- length(cat$group_specs)
  has_or <- length(cat$or_blocks) > 0L
  has_trt <- length(cat$treatment_specs) > 0L

  # A pooled summary and fit each own block records and accounting. The fit
  # additionally owns its slope_fit (including a second R and coefficients).
  wls <- function(bn) {
    block <- R(c(r11 = A(1), r1d = A(p, z), r1y = A(1), n = A(1),
      weight_sum = A(1), mean_d = A(p, z), mean_y = A(1),
      augmented_rows = A(1), residual_rows = A(1)))
    blocks <- R(rep(block, length(bn)), nm = bn)
    accounting <- R(c(n = A(1), p = A(1), n_blocks = A(1), residual_dim = A(2),
      max_augmented_rows = A(1), max_residual_rows = A(1),
      max_chunk_merge_rows = A(1), max_residual_merge_rows = A(1), n_compressions = A(1)))
    simple_summary <- R(c(factor = M(p + 1, p + 1), n = A(1), columns = E(z)),
      "cdid_or_wls_summary")
    simple_fit <- R(c(coefficients = A(p, z), rss = A(1), R = M(p, p),
      pivot = A(p), n = A(1), rank = A(1), rank_tol = A(1), cutoff = A(1)), "cdid_or_wls_fit")
    if (recipe$intercept) {
      summary <- R(c(residual = simple_summary, blocks = blocks, n = A(1), p = A(1),
        columns = E(z), accounting = accounting), "cdid_or_pooled_wls_summary")
      fit <- R(c(coefficients = A(p, z), intercepts = A(length(bn), bn), rss = A(1),
        R = M(p, p), pivot = A(p), n = A(1), p = A(1), n_blocks = A(1),
        rank = A(1), rank_tol = A(1), cutoff = A(1), blocks = blocks,
        slope_fit = if (p > 0) simple_fit else 0, accounting = accounting), "cdid_or_pooled_wls_fit")
    } else {
      summary <- simple_summary
      fit <- simple_fit
    }
    basis <- R(c(slope_columns = E(z), kept_columns = E(z), kept_index = A(p),
      proposed_columns = E(z), aliased_columns = E(z), nontrivial_candidates_retained = E(z),
      alias_map = M(p, p, z, z, TRUE), block_constants = M(length(bn), p, bn, z, TRUE),
      reconstruction_control_scaled = A(p), reconstruction_treated_scaled = A(p),
      reconstruction_roundoff_threshold = A(1), rows_checked = A(1),
      proposal = E("R stats::lm.fit / dqrdc2 pivoted QR at rank_tol on the pooled (partialled) factor"),
      rank_check = E("unchanged relative-diagonal rule of the accepted kernel on the reduced factor"),
      frozen_during_differentiation = A(1), numerical_certificate = A(1)))
    c(fit = fit, summary = summary, basis = basis, rss = A(1), objective = A(1))
  }

  cells <- V(vapply(cat$cell_specs, function(cs) E(cs) + R(c(
    group_id = A(1), local_id = A(1), n_pair = A(1), weight_sum = A(1), mass = A(1),
    alpha = A(1), mean_d = A(p, z), mean_y = A(1), a_k = A(1), M_one = A(1),
    Y_one = A(1), M_z = A(p, z), Y_z = A(p, z), M_d_kept = A(p, z))), numeric(1)))
  groups <- V(vapply(cat$group_specs, function(gs) {
    bn <- vapply(cat$cell_specs[gs$cell_ids], `[[`, character(1), "name")
    E(gs) + R(c(keep = A(p), p = A(1), kept_columns = E(z), kept_cols_F = A(p),
      beta = A(p, z), alpha = A(length(bn), bn), wls(bn)))
  }, numeric(1)))
  shared <- if (has_or) {
    bn <- vapply(cat$or_blocks, `[[`, character(1), "name")
    blocks <- V(vapply(cat$or_blocks, function(b) {
      k <- length(b$control_cell_ids)
      E(b) + R(c(alpha = A(1), mass = A(1), weight_sum = A(1), mean_d = A(p, z),
        a_force_alpha = A(k), a_force_beta = M(p, k), j_force_alpha_alpha = A(k),
        j_force_alpha_beta = k + 8 + k * A(p), j_force_beta_alpha = M(p, k),
        j_force_beta_beta = k + 8 + k * M(p, p, z, z, TRUE)))
    }, numeric(1)))
    nt <- length(bn) + p
    R(c(keep = A(p), p = A(1), kept_columns = E(z), kept_cols_F = A(p),
      intercept = A(1), n_blocks = A(1), alpha = A(length(bn), bn), beta = A(p, z),
      theta = A(nt, paste0("OR.theta", seq_len(nt) - 1L)), blocks = blocks, wls(bn)))
  } else 0
  treatments <- if (has_trt) V(vapply(cat$treatment_specs, function(ts) {
    E(ts) + R(c(coefficients = A(q, paste0("gamma", seq_len(q) - 1L)),
      column_names = E(trt_recipe$columns),
      H_chol = M(q, q, trt_recipe$columns, trt_recipe$columns, TRUE),
      # The current lfit$H partially matches H_chol and retains that alias too.
      H = M(q, q, trt_recipe$columns, trt_recipe$columns, TRUE),
      n_risk_F = A(1), a_forcing = M(q, length(ts$all_cohort_indices)),
      clipping = E(ipw_clip_summary(prob = 0.5, controls = FALSE, weights = 1))))
  }, numeric(1))) else 0
  dense_adjoint <- R(c(tilde_lambda = A(nc), tilde_nu = ng + 8 + ng * A(p)))
  adjoint <- function(k) {
    g <- min(k, ng)
    compact <- R(c(kind = E("compact_reconstruction_adjoint"), cell_ids = A(k),
      tilde_lambda = A(k), group_ids = A(g), tilde_nu = g + 8 + g * A(p)))
    if (dense_adjoints) dense_adjoint else min(dense_adjoint, compact)
  }
  links <- V(vapply(cat$links, function(lk) {
    k <- length(lk$control_cell_ids)
    E(lk) + R(c(mu = A(1), mu_ipw = A(1), mu_or = A(1), mu_reg = A(1), mu_dr = A(1),
      d_plus = A(1), num_T = A(1), T_gl = A(1), num_Bg = A(1), B_gl = A(1),
      d_minus = A(1), num_C = A(1), C_ol = A(1), num_Bo = A(1), B_ol = A(1),
      d_bar_plus_j = A(p, z), d_bar_plus_or = A(p, z), d_bar_minus_or = A(p, z),
      d_bar_diff_or = A(p, z), a_coeff_T_g = A(1), a_coeff_Bg_g = A(1),
      a_coeff_C_k = A(k), a_coeff_Bo_k = A(k), m_minus_k = A(k),
      d_bar_minus_j_k = k + 8 + k * A(p), denminus_gamma = A(q),
      g_C_gamma = A(q, trt_recipe$columns), g_Bo_gamma = A(q, trt_recipe$columns),
      g_dr_gamma = A(q, trt_recipe$columns), adjoint_T = adjoint(1L),
      if (has_trt) c(adjoint_C = adjoint(k))))
  }, numeric(1)))
  # cp04_payload charges logical aliases separately, even if R shares storage.
  c(cells = cells, groups = groups, shared_or = shared,
    shared = if (has_or) shared else E(list(p = 0L, keep = integer(), beta = numeric())),
    treatments = treatments, links = links)
}

augmented_retained_core_bound <- function(cat, recipe, trt_recipe, theta_ids, dense_adjoints = FALSE) {
  E <- function(x) unname(cp04_payload(x)["entries"])
  A <- function(n, nm = NULL) n + 8 + if (is.null(nm)) 0 else 9 + E(nm)
  M <- function(nr, nc, rn = NULL, cn = NULL, dn = FALSE) {
    nr * nc + 27 + if (dn) 1 + E(list(rn, cn)) else 0
  }
  R <- function(fields, cls = NULL, nm = names(fields)) {
    length(fields) + 8 + sum(fields) +
      if (is.null(nm) && is.null(cls)) 0 else
        8 + as.integer(!is.null(nm)) + as.integer(!is.null(cls)) + E(nm) + E(cls)
  }
  V <- function(x) length(x) + 8 + sum(x)
  p <- as.double(recipe$p)
  z <- recipe$slope_columns
  q <- as.double(trt_recipe$p)
  nc <- length(cat$cell_specs)
  ng <- length(cat$group_specs)
  dq <- vapply(theta_ids, length, integer(1L))
  pure <- imputation_retained_core_bound(cat, recipe, trt_recipe, dense_adjoints = dense_adjoints)

  cells <- V(vapply(seq_len(nc), function(cid) {
    cs <- cat$cell_specs[[cid]]
    E(cs) + R(c(group_id = A(1), local_id = A(1), n_pair = A(1), weight_sum = A(1), mass = A(1),
      alpha = A(1), mean_d = A(p, z), mean_y = A(1), a_k = A(1), M_one = A(1),
      Y_one = A(1), M_z = A(p, z), Y_z = A(p, z), M_d_kept = A(p, z))) +
      R(c(q_id = A(1), theta_ids = A(dq[cid]), Yaug_one = A(1), Yaug_z = A(p, z),
        MW_one = A(1), MW_d = A(p, z), MD_one = A(1), MD_d = A(p, z), Rq_one = A(dq[cid]),
        clipping = E(ipw_clip_summary(q = 0.5, weights = 1))))
  }, numeric(1L)))

  block_q <- vapply(cat$or_blocks, function(b) {
    d <- sum(dq[b$control_cell_ids])
    2 * A(d) + 64
  }, numeric(1L))
  d_global <- sum(dq[sort(unique(unlist(lapply(cat$or_blocks, `[[`, "control_cell_ids"))))])
  shared <- unname(pure["shared_or"]) + sum(block_q) + A(d_global) + M(p, d_global) + 64

  dense_adjoint <- R(c(tilde_lambda = A(nc), tilde_nu = ng + 8 + ng * A(p)))
  adjoint <- function(k) {
    g <- min(k, ng)
    compact <- R(c(kind = E("compact_reconstruction_adjoint"), cell_ids = A(k),
      tilde_lambda = A(k), group_ids = A(g), tilde_nu = g + 8 + g * A(p)))
    if (dense_adjoints) dense_adjoint else min(dense_adjoint, compact)
  }
  links <- V(vapply(cat$links, function(lk) {
    k <- length(lk$control_cell_ids)
    dg <- dq[lk$target_cell_id]
    dc <- sum(dq[lk$control_cell_ids])
    E(lk) + R(c(mu = A(1), mu_dr = A(1),
      d_plus = A(1), num_T = A(1), T_gl = A(1), num_Bg = A(1), B_gl = A(1),
      d_minus = A(1), num_C = A(1), C_ol = A(1), num_Bo = A(1), B_ol = A(1),
      d_bar_plus_or = A(p, z), d_bar_minus_or = A(p, z), d_bar_diff_or = A(p, z),
      a_coeff_T_g = A(1), a_coeff_Bg_g = A(1), a_coeff_C_k = A(k), a_coeff_Bo_k = A(k),
      m_minus_k = A(k), lam_T = A(1), nu_T = A(p, z), lam_C_k = A(k),
      nu_C_k = k + 8 + k * A(p, z), denminus_gamma = A(q),
      g_C_gamma = A(q, trt_recipe$columns), g_Bo_gamma = A(q, trt_recipe$columns),
      g_dr_gamma = A(q, trt_recipe$columns),
      sel_T_ids = A(dg), sel_T = A(dg), sel_C_ids = A(dc), sel_C = A(dc),
      adjoint_T = adjoint(1L), adjoint_C = adjoint(k)))
  }, numeric(1L)))

  # The pure-imputation convention charges the shared-OR alias field twice.
  c(cells = cells, groups = unname(pure["groups"]), shared_or = shared, shared = shared,
    treatments = unname(pure["treatments"]), links = links)
}
