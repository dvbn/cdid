# SPDX-License-Identifier: GPL-3.0-only
#
# Deterministic backend tests for CP03-B (`pt_assumption = "all"` under explicit
# `pt_start` and D34 no-extra-window handling) across:
#   - `R/chained_design.R`
#   - `R/chained_gmm_core.R`
#   - `R/chained_ipw_shared.R`
#   - `R/chained_ipw_pipeline.R`
#
# Directly exercises the accepted tiny references in `spec/pt_all_reference.py`,
# `spec/pt_all_reference_check.R`, and `spec/pt_all_reference_results.json`
# using literal rational / small-matrix expectations derived independently of
# the package adapter. No raw DGP fits, Monte Carlo simulations, or bootstrap
# campaigns are executed.

if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}

make_backend_test_gmm_controls <- function(
  method = "conventional_chol_qr",
  weighting = "optimal",
  require_matching_scale = FALSE,
  covariance_strategy = NULL
) {
  is_opt <- identical(weighting, "optimal")
  acc <- gmm_accuracy_control(
    method = method,
    target_error_unit = "abs",
    max_target_error_abs = if (identical(method, "conventional_chol_qr")) .Machine$double.xmax else 1e-6,
    max_target_error_se = NA_real_,
    max_orthogonality_defect = 1e-8,
    min_certified_lambda_rel = if (is_opt) 1e-10 else NA_real_,
    max_factor_residual_rel = 1e-10,
    min_rcond_R = if (is_opt) 1e-10 else NA_real_,
    min_rcond_U = 1e-10,
    min_sigma_rel_A = 1e-10,
    max_kkt_residual_rel = 1e-10,
    max_loading_residual_rel = 1e-10,
    require_matching_scale = require_matching_scale,
    max_loading_error_abs = if (identical(method, "conventional_chol_qr")) NULL else 1e-6,
    covariance_strategy = covariance_strategy
  )
  bnds <- gmm_bounds_control(
    max_matrix_entries = 100000L,
    max_batch_rows = 32L,
    max_targets = 32L,
    max_open_clusters = 64L,
    allow_materialize_if = TRUE,
    max_materialized_if_entries = 100000L
  )
  list(accuracy = acc, bounds = bnds)
}

# Build a deterministic centered 2M x M IF matrix whose sample covariance
# crossprod(Phi) / (2M) equals a supplied M x M SPD/PSD matrix `Omega`.
make_exact_if_replay_from_omega <- function(Omega, reduction = list(mode = "none")) {
  M <- nrow(Omega)
  ev <- eigen(Omega, symmetric = TRUE)
  vals <- pmax(0.0, ev$values)
  L_mat <- ev$vectors %*% diag(sqrt(vals), nrow = M, ncol = M)
  base_rows <- sqrt(as.double(M)) * t(L_mat)
  Phi <- rbind(base_rows, -base_rows)
  n_units <- nrow(Phi)
  colnames(Phi) <- colnames(Omega)
  uid <- paste0("u", seq_len(n_units))
  wi <- rep(1.0, n_units)
  cb <- function(consume) {
    consume(list(
      unit_id = uid,
      cluster_id = NULL,
      phi = Phi,
      weight = wi,
      nuisance_fit_id = "frozen_backend_test"
    ))
    invisible(NULL)
  }
  replay_spec <- gmm_if_replay_spec(
    cb, "unit", FALSE, n_units, n_units, n_units,
    TRUE, FALSE, "unweighted", 1 / n_units,
    "backend_test_scale", "frozen_backend_test", uid, NULL
  )
  prov <- gmm_covariance_provenance(
    row_level = "unit",
    n_units = n_units,
    n_clusters = NA,
    n_positive_contributing_rows = n_units,
    centered = TRUE,
    weighted = FALSE,
    weight_type = "unweighted",
    normalization_factor = 1 / n_units,
    scale_id = "backend_test_scale",
    nuisance_fit_id = "frozen_backend_test",
    reduction = reduction
  )
  prov$complete_unit_if_gram <- ipw_seal(
    list(covariance = Omega, n_units = prov$n_units, nuisance_fit_id = prov$nuisance_fit_id),
    "complete_unit_if_gram"
  )
  list(Phi = Phi, n_units = n_units, replay_spec = replay_spec, provenance = prov)
}

run_backend_gmm_solve <- function(
  design,
  mhat,
  V,
  weighting = "optimal",
  method = "conventional_chol_qr",
  targets = NULL,
  gauge_rule = "lowest",
  covariance_strategy = NULL,
  reduction = list(mode = "none"),
  request_param_cov = TRUE,
  materialize_target_if = TRUE
) {
  lids <- design$incidence$link_ids
  mhat_named <- stats::setNames(as.double(mhat), lids)
  V_named <- as.matrix(V)
  dimnames(V_named) <- list(lids, lids)
  if (is.null(targets)) {
    pt <- design$post_targets
    targets <- lapply(seq_len(nrow(pt)), function(i) {
      gmm_target_spec(pt$param_id[i], pt$param_index[i], 1.0)
    })
  }
  ctrls <- make_backend_test_gmm_controls(
    method = method,
    weighting = if (identical(weighting, "two_step")) "optimal" else weighting,
    covariance_strategy = covariance_strategy
  )
  if_fx <- make_exact_if_replay_from_omega(V_named, reduction = reduction)
  solve_chained_gmm_core(
    design = design,
    mhat = mhat_named,
    V = V_named,
    provenance = if_fx$provenance,
    targets = targets,
    weighting = weighting,
    if_replay = if_fx$replay_spec,
    accuracy = ctrls$accuracy,
    bounds = ctrls$bounds,
    gauge_rule = gauge_rule,
    request_param_cov = request_param_cov,
    materialize_target_if = materialize_target_if
  )
}

test_that("Fixture 1: SPD zero-loading row strictly improves efficiency (variance 5/3 vs 2 if discarded) and matches Lagrange/Schur references", {
  # Accepted Fixture 1 (`spec/pt_all_reference.py`, `spec/pt_all_reference_check.R`):
  # Cohort g = 4, calendar T = {1, 2, 3, 4}, delta = 0, b_g = 3, pt_start c = 2.
  # Links: l1 = (1, 2), l2 = (2, 3), l3 = (3, 4).
  links_f1 <- data.frame(
    link_id = c("l12", "l23", "l34"),
    g = c(4, 4, 4),
    s = c(1, 2, 3),
    t = c(2, 3, 4),
    functional_key = rep("NYT:{Inf}", 3L),
    control_pool = rep("{Inf}", 3L),
    s_qualified = c(FALSE, FALSE, TRUE),
    t_qualified = c(FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  cal <- c(1, 2, 3, 4)
  cohorts <- c(4, Inf)

  # 1. Bitwise post-mode invariance when pt_assumption = "post"
  cs_default <- validate_calendar_spec(cal, cohorts, anticipation = 0)
  cs_post <- validate_calendar_spec(cal, cohorts, anticipation = 0, pt_assumption = "post", pt_start = NULL)
  expect_identical(cs_default, cs_post)

  des_default <- build_calendar_chained_design(links_f1, calendar = cal, cohorts = cohorts, anticipation = 0)
  des_post <- build_calendar_chained_design(
    links_f1, calendar = cal, cohorts = cohorts, anticipation = 0,
    pt_assumption = "post", pt_start = NULL
  )
  expect_identical(des_default, des_post)

  # 2. PT-all restriction at c = 2 via both direct builder and adapter
  des_all <- build_calendar_chained_design(
    links_f1, calendar = cal, cohorts = cohorts, anticipation = 0,
    pt_assumption = "all", pt_start = 2
  )
  des_all_adapted <- restrict_chained_design_pt_all(des_post, pt_start = 2)
  expect_identical(des_all, des_all_adapted)

  # Vertex catalog is preserved (ordered by role: reference, pre, post); only free parameter indexing changes
  expect_identical(nrow(des_all$vertices), 4L)
  expect_identical(des_all$vertices$vertex_id, des_post$vertices$vertex_id)
  expect_identical(des_all$vertices$time, c(3, 1, 2, 4))
  expect_identical(des_all$vertices$role, c("reference", "pre", "pre", "post"))
  expect_identical(des_all$vertices$param_index, c(NA_integer_, 1L, NA_integer_, 2L))
  expect_identical(des_all$vertices$anchored_value, c(0, NA_real_, 0, NA_real_))

  # Link (2,3) has both endpoints anchored -> (0L, 0L) in incidence
  expect_identical(des_all$incidence$n_links, 3L)
  expect_identical(des_all$incidence$n_parameters, 2L)
  expect_identical(des_all$incidence$start, c(1L, 0L, 0L))
  expect_identical(des_all$incidence$end, c(0L, 0L, 2L))
  expect_identical(des_all$pt_all_spec$zero_loading_link_indices, 2L)
  expect_identical(des_all$pt_all_spec$free_post_param_indices, c(1L, 3L))
  expect_identical(des_all$pt_all_spec$restricted_post_param_indices, 2L)

  # Exact matrix identity W_all = W_post %*% H1
  W_post <- materialize_chained_design(des_post, max_entries = 100L)$W
  W_all <- materialize_chained_design(des_all, max_entries = 100L)$W
  H1 <- rbind(c(1L, 0L), c(0L, 0L), c(0L, 1L))
  expect_identical(unname(W_post), rbind(c(-1L, 1L, 0L), c(0L, -1L, 0L), c(0L, 0L, 1L)))
  expect_identical(unname(W_all), rbind(c(-1L, 0L), c(0L, 0L), c(0L, 1L)))
  expect_identical(unname(W_all * 1.0), unname(W_post %*% H1))

  # 3. Solve unrestricted post GLS vs PT-all restricted GLS on Fixture 1 SPD V1 and m1
  V1 <- matrix(c(
     4, -1,  0,
    -1,  3, -1,
     0, -1,  2
  ), nrow = 3L, byrow = TRUE)
  m1 <- c(3 / 2, 1 / 2, 5 / 2)
  targets_both <- list(
    gmm_target_spec("beta_1", 1L, 1.0),
    gmm_target_spec("tau_4", 2L, 1.0)
  )

  for (meth in c("conventional_chol_qr", "certified_eigen_chol_kkt")) {
    fit_post <- run_backend_gmm_solve(des_post, m1, V1, weighting = "optimal", method = meth)
    fit_all <- run_backend_gmm_solve(des_all, m1, V1, weighting = "optimal", method = meth, targets = targets_both)

    # Exact Fixture 1 rational reference:
    # eta_hat = (-5/3, 8/3), cov_eta = [[11/3, 1/3], [1/3, 5/3]],
    # Lambda_eta = [[-1, 0], [-1/3, 1/3], [0, 1]] (nonzero loading [-1/3, 1/3] on zero row l23!)
    expect_equal(unname(fit_all$internal_gauge$theta_gauge_coords), c(-5 / 3, 8 / 3), tolerance = 1e-13)
    expect_equal(unname(fit_all$targets$estimate), c(-5 / 3, 8 / 3), tolerance = 1e-13)
    expect_equal(unname(fit_all$targets$se_sandwich^2), c(11 / 3, 5 / 3), tolerance = 1e-13)
    expect_equal(
      unname(fit_all$internal_gauge$param_cov_gauge),
      matrix(c(11 / 3, 1 / 3, 1 / 3, 5 / 3), nrow = 2L, byrow = TRUE),
      tolerance = 1e-13
    )
    expect_equal(
      unname(fit_all$internal_gauge$Lambda),
      rbind(c(-1, 0), c(-1 / 3, 1 / 3), c(0, 1)),
      tolerance = 1e-13
    )
    psi_all <- fit_all$if_projection$materialized_target_if
    expect_equal(unname(crossprod(psi_all) / nrow(psi_all)), unname(fit_all$internal_gauge$param_cov_gauge), tolerance = 1e-13)

    # Explicit Lagrange multiplier & Schur-complement equality against unrestricted post GLS:
    Sigma_unr <- fit_post$internal_gauge$param_cov_gauge
    th_unr <- fit_post$internal_gauge$theta_gauge_coords
    th_lagrange <- th_unr - (Sigma_unr[, 2L] / Sigma_unr[2L, 2L]) * th_unr[2L]
    Sigma_schur <- Sigma_unr - tcrossprod(Sigma_unr[, 2L]) / Sigma_unr[2L, 2L]
    expect_equal(unname(th_lagrange), c(-5 / 3, 0, 8 / 3), tolerance = 1e-13)
    expect_equal(unname(th_lagrange[c(1L, 3L)]), unname(fit_all$internal_gauge$theta_gauge_coords), tolerance = 1e-13)
    expect_equal(unname(Sigma_schur[c(1L, 3L), c(1L, 3L)]), unname(fit_all$internal_gauge$param_cov_gauge), tolerance = 1e-13)
  }

  # 4. Analytic counterexample: discarding zero-loading row l23 changes eta_hat from (-5/3, 8/3)
  #    to (-3/2, 5/2) and inflates ATT(4,4) variance from 5/3 to 2 (strictly worse by 1/3!)
  keep_active <- c(1L, 3L)
  W_drop <- W_all[keep_active, , drop = FALSE]
  V_drop <- V1[keep_active, keep_active, drop = FALSE]
  m_drop <- m1[keep_active]
  cov_drop <- solve(t(W_drop) %*% solve(V_drop, W_drop))
  eta_drop <- as.double(cov_drop %*% t(W_drop) %*% solve(V_drop, m_drop))
  expect_equal(eta_drop, c(-3 / 2, 5 / 2), tolerance = 1e-13)
  expect_equal(unname(cov_drop), diag(c(4, 2)), tolerance = 1e-13)
  expect_equal(cov_drop[2L, 2L] - (5 / 3), 1 / 3, tolerance = 1e-13)

  # 5. Analytic counterexample: zeroing row/col 2 of V1 is singular and refused under optimal GMM
  V1_zeroed <- V1
  V1_zeroed[2L, ] <- 0
  V1_zeroed[, 2L] <- 0
  expect_error(
    run_backend_gmm_solve(des_all, m1, V1_zeroed, weighting = "optimal", method = "conventional_chol_qr"),
    class = "cdid_gmm_pd_error"
  )

  # 6. Target-specific non-improvement when c_A is orthogonal to W_A' tilde(Omega)_AA^{-1} Omega_AZ
  links_4 <- data.frame(
    link_id = c("l13", "l23", "l24", "l34"),
    g = c(4, 4, 4, 4),
    s = c(1, 2, 2, 3),
    t = c(3, 3, 4, 4),
    functional_key = rep("NYT:{Inf}", 4L),
    control_pool = rep("{Inf}", 4L),
    s_qualified = c(FALSE, FALSE, FALSE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_4 <- build_calendar_chained_design(
    links_4, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 2
  )
  Omega_orth <- matrix(c(
    2,  1, 0,  0,
    1,  2, 1, -1,
    0,  1, 2,  0,
    0, -1, 0,  2
  ), nrow = 4L, byrow = TRUE)
  m_orth <- c(1, 2, 3, 4)
  fit_orth <- run_backend_gmm_solve(des_4, m_orth, Omega_orth, weighting = "optimal")
  # Full Sigma = [[1, 0], [0, 1]]; without zero-loading row l23, Sigma_drop = [[2, 0], [0, 1]]:
  # ATT(4,4) variance is 1 in both, while free pre parameter variance improves from 2 to 1!
  expect_equal(unname(fit_orth$internal_gauge$param_cov_gauge), diag(c(1, 1)), tolerance = 1e-13)
  W_4 <- materialize_chained_design(des_4, 100L)$W
  Sigma_orth_drop <- solve(t(W_4[c(1L, 3L, 4L), ]) %*% solve(Omega_orth[c(1L, 3L, 4L), c(1L, 3L, 4L)], W_4[c(1L, 3L, 4L), ]))
  expect_equal(unname(Sigma_orth_drop), diag(c(2, 1)), tolerance = 1e-13)
})

test_that("Fixture 2: multi-baseline star vs chain basis invariance, identity non-invariance, and joint 2-post-target dominance", {
  # Part 2A (`spec/pt_all_reference.py`):
  # Star basis (1,4), (2,4), (3,4) vs Chain basis (1,2), (2,3), (3,4) under c = 1
  links_star <- data.frame(
    link_id = c("s14", "s24", "s34"),
    g = c(4, 4, 4),
    s = c(1, 2, 3),
    t = c(4, 4, 4),
    functional_key = rep("NYT:{Inf}", 3L),
    control_pool = rep("{Inf}", 3L),
    s_qualified = c(FALSE, FALSE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  links_chain <- data.frame(
    link_id = c("c12", "c23", "c34"),
    g = c(4, 4, 4),
    s = c(1, 2, 3),
    t = c(2, 3, 4),
    functional_key = rep("NYT:{Inf}", 3L),
    control_pool = rep("{Inf}", 3L),
    s_qualified = c(FALSE, FALSE, TRUE),
    t_qualified = c(FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_star <- build_calendar_chained_design(
    links_star, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 1
  )
  des_chain <- build_calendar_chained_design(
    links_chain, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 1
  )

  # Star has WH = (1, 1, 1)'; Chain has WH = (0, 0, 1)' = A_chain %*% WH_star (two zero rows!)
  expect_identical(des_star$incidence$start, c(0L, 0L, 0L))
  expect_identical(des_star$incidence$end, c(1L, 1L, 1L))
  expect_identical(des_chain$incidence$start, c(0L, 0L, 0L))
  expect_identical(des_chain$incidence$end, c(0L, 0L, 1L))
  expect_identical(des_chain$pt_all_spec$zero_loading_link_indices, c(1L, 2L))

  V_star <- matrix(c(
    3, 1, 1,
    1, 2, 1,
    1, 1, 2
  ), nrow = 3L, byrow = TRUE)
  d_hat <- c(1.4, 1.1, 1.0)
  A_chain <- rbind(
    c(1, -1,  0),
    c(0,  1, -1),
    c(0,  0,  1)
  )
  V_chain <- A_chain %*% V_star %*% t(A_chain)
  m_chain <- as.double(A_chain %*% d_hat)

  fit_star_opt <- run_backend_gmm_solve(des_star, d_hat, V_star, weighting = "optimal")
  fit_chain_opt <- run_backend_gmm_solve(des_chain, m_chain, V_chain, weighting = "optimal")

  # Exact rational reference: tau_star = tau_chain = 28/25, var = 7/5,
  # weights_star = (1/5, 2/5, 2/5), Lambda_chain = (1/5, 3/5, 1)
  expect_equal(unname(fit_star_opt$targets$estimate), 28 / 25, tolerance = 1e-13)
  expect_equal(unname(fit_chain_opt$targets$estimate), 28 / 25, tolerance = 1e-13)
  expect_equal(unname(fit_star_opt$targets$se_sandwich^2), 7 / 5, tolerance = 1e-13)
  expect_equal(unname(fit_chain_opt$targets$se_sandwich^2), 7 / 5, tolerance = 1e-13)
  expect_equal(unname(as.double(fit_star_opt$internal_gauge$Lambda)), c(1 / 5, 2 / 5, 2 / 5), tolerance = 1e-13)
  expect_equal(unname(as.double(fit_chain_opt$internal_gauge$Lambda)), c(1 / 5, 3 / 5, 1), tolerance = 1e-13)
  expect_equal(unname(as.double(crossprod(A_chain, fit_chain_opt$internal_gauge$Lambda))), c(1 / 5, 2 / 5, 2 / 5), tolerance = 1e-13)

  # Identity weighting counterexample: Star gives mean(d_hat) = 7/6, while Chain gives m_(3,4) = 1.0
  fit_star_id <- run_backend_gmm_solve(des_star, d_hat, V_star, weighting = "identity")
  fit_chain_id <- run_backend_gmm_solve(des_chain, m_chain, V_chain, weighting = "identity")
  expect_equal(unname(fit_star_id$targets$estimate), 7 / 6, tolerance = 1e-13)
  expect_equal(unname(fit_chain_id$targets$estimate), 1.0, tolerance = 1e-13)
  expect_gt(abs(fit_star_id$targets$estimate - fit_chain_id$targets$estimate), 0.16)

  # Part 2B (`spec/pt_all_reference.py`): Joint restricted GMM for two post periods (tau_4, tau_5)
  # vs separate per-cell averaging on links (2,4), (3,4), (3,5), (4,5) under c = 2
  links_2post <- data.frame(
    link_id = c("l24", "l34", "l35", "l45"),
    g = rep(4, 4L),
    s = c(2, 3, 3, 4),
    t = c(4, 4, 5, 5),
    functional_key = rep("NYT:{Inf}", 4L),
    control_pool = rep("{Inf}", 4L),
    s_qualified = c(FALSE, TRUE, TRUE, TRUE),
    t_qualified = rep(TRUE, 4L),
    stringsAsFactors = FALSE
  )
  des_2post <- build_calendar_chained_design(
    links_2post, calendar = c(2, 3, 4, 5), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 2
  )
  V_2post <- matrix(c(
     4,  1, 1, -1,
     1,  3, 1, -1,
     1,  1, 4,  1,
    -1, -1, 1,  3
  ), nrow = 4L, byrow = TRUE)
  m_2post <- c(6 / 5, 1, 2, 9 / 10)

  fit_2post <- run_backend_gmm_solve(des_2post, m_2post, V_2post, weighting = "optimal")
  # Exact rational reference: eta_joint = (173/160, 159/80), cov_joint = [[35/16, 9/8], [9/8, 11/4]],
  # strictly dominating per-cell estimates (27/25, 2) with variances (11/5, 4)
  expect_equal(unname(fit_2post$targets$estimate), c(173 / 160, 159 / 80), tolerance = 1e-13)
  expect_equal(
    unname(fit_2post$internal_gauge$param_cov_gauge),
    matrix(c(35 / 16, 9 / 8, 9 / 8, 11 / 4), nrow = 2L, byrow = TRUE),
    tolerance = 1e-13
  )
  expect_lt(fit_2post$targets$se_sandwich[1L]^2, 11 / 5)
  expect_lt(fit_2post$targets$se_sandwich[2L]^2, 4.0)
})

test_that("Fixture 3: early cancellation bridge (u < c < b_g) preserves operator separation and avoids pre-c zeroing or merging bias", {
  # Accepted Fixture 3 (`spec/pt_all_reference.py`):
  # Cohort g = 4, calendar T = {1, 2, 3, 4, 5}, delta = 0, b_4 = 3, pt_start c = 2.
  # Operators K1 = NYT:{5,Inf} and K2 = NYT:{Inf}.
  # Links: l1=(1,2) K1, l2=(1,2) K2, l3=(1,4) K2, l4=(2,3) K1, l5=(2,4) K1.
  links_f3 <- data.frame(
    link_id = c("l12_K1", "l12_K2", "l14_K2", "l23_K1", "l24_K1"),
    g = rep(4, 5L),
    s = c(1, 1, 1, 2, 2),
    t = c(2, 2, 4, 3, 4),
    functional_key = c("NYT:{5,Inf}", "NYT:{Inf}", "NYT:{Inf}", "NYT:{5,Inf}", "NYT:{5,Inf}"),
    control_pool = c("{5,Inf}", "{Inf}", "{Inf}", "{5,Inf}", "{5,Inf}"),
    s_qualified = rep(FALSE, 5L),
    t_qualified = c(FALSE, FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  cal <- c(1, 2, 3, 4, 5)
  cohorts <- c(4, 5, Inf)

  des_f3 <- build_calendar_chained_design(
    links_f3, calendar = cal, cohorts = cohorts, anticipation = 0,
    pt_assumption = "all", pt_start = 2
  )

  # Pre-start vertices at u = 1 keep distinct private parameters for K1 and K2;
  # in-window vertices at u = 2 (both K1 and K2) and reference u = 3 are anchored at 0!
  v_u1 <- des_f3$vertices[des_f3$vertices$time == 1, , drop = FALSE]
  expect_identical(nrow(v_u1), 2L)
  expect_identical(v_u1$functional_key, c("NYT:{5,Inf}", "NYT:{Inf}"))
  expect_identical(v_u1$param_index, c(1L, 2L))

  v_u2 <- des_f3$vertices[des_f3$vertices$time == 2, , drop = FALSE]
  expect_identical(nrow(v_u2), 2L)
  expect_identical(v_u2$functional_key, c("NYT:{5,Inf}", "NYT:{Inf}"))
  expect_true(all(is.na(v_u2$param_index)))
  expect_identical(v_u2$anchored_value, c(0, 0))

  WH_f3 <- materialize_chained_design(des_f3, 100L)$W
  expect_identical(
    unname(WH_f3),
    rbind(
      c(-1L,  0L, 0L),
      c( 0L, -1L, 0L),
      c( 0L, -1L, 1L),
      c( 0L,  0L, 0L),
      c( 0L,  0L, 1L)
    )
  )

  mu_pop <- c(3 / 2, 1 / 2, 7 / 2, 0, 3)
  V_f3 <- matrix(c(
     3, 1, 0, -1,  0,
     1, 3, 1,  0,  0,
     0, 1, 4,  0,  1,
    -1, 0, 0,  2, -1,
     0, 0, 1, -1,  3
  ), nrow = 5L, byrow = TRUE)

  fit_f3 <- run_backend_gmm_solve(des_f3, mu_pop, V_f3, weighting = "optimal")
  # Exact Fixture 3 reference: eta_correct = (-3/2, -1/2, 3), zero bias on tau_4 = 3!
  expect_equal(unname(fit_f3$internal_gauge$theta_gauge_coords), c(-3 / 2, -1 / 2, 3), tolerance = 1e-13)
  expect_equal(unname(fit_f3$targets$estimate), 3.0, tolerance = 1e-13)

  # Analytic counterexample A: falsely zeroing before c = 2 (i.e. setting pt_start = 1) gives bias +7/18
  des_false_zero <- build_calendar_chained_design(
    links_f3, calendar = cal, cohorts = cohorts, anticipation = 0,
    pt_assumption = "all", pt_start = 1
  )
  fit_false_zero <- run_backend_gmm_solve(des_false_zero, mu_pop, V_f3, weighting = "optimal")
  expect_equal(unname(fit_false_zero$targets$estimate - 3.0), 7 / 18, tolerance = 1e-13)

  # Analytic counterexample B: falsely merging K1 and K2 at u = 1 < c gives bias +1/34
  WH_false_merge <- rbind(
    c(-1, 0),
    c(-1, 0),
    c(-1, 1),
    c( 0, 0),
    c( 0, 1)
  )
  eta_false_merge <- as.double(solve(t(WH_false_merge) %*% solve(V_f3, WH_false_merge), t(WH_false_merge) %*% solve(V_f3, mu_pop)))
  expect_equal(eta_false_merge[2L] - 3.0, 1 / 34, tolerance = 1e-13)
})

test_that("Fixture 4 & D34: irregular calendar with anticipation preserves b_g and free anticipation while handling global/partial no-extra-window and D33 checks", {
  # Accepted Fixture 4 (`spec/pt_all_reference.py`):
  # Calendar {2010, 2012, 2015, 2017, 2020}, g = 2017, delta = 2 => b_g = 2012.
  cal <- c(2010, 2012, 2015, 2017, 2020)
  cohorts <- c(2017, Inf)
  links_irreg <- data.frame(
    link_id = c("l_10_12", "l_12_15", "l_12_17", "l_17_20"),
    g = rep(2017, 4L),
    s = c(2010, 2012, 2012, 2017),
    t = c(2012, 2015, 2017, 2020),
    functional_key = rep("NT:{Inf}", 4L),
    control_pool = rep("{Inf}", 4L),
    s_qualified = c(FALSE, TRUE, TRUE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE, TRUE),
    s_coordinate_role = c("pre_pt_all_zero", "reference", "reference", "post_att"),
    t_coordinate_role = c("reference", "anticipation_att", "post_att", "post_att"),
    stringsAsFactors = FALSE
  )

  des_2010 <- build_calendar_chained_design(
    links_irreg, calendar = cal, cohorts = cohorts, anticipation = 2,
    pt_assumption = "all", pt_start = 2010
  )
  expect_identical(des_2010$links$ref, rep(2012, 4L))
  V <- des_2010$vertices
  expect_identical(V$time, c(2012, 2010, 2015, 2017, 2020))
  expect_identical(V$role, c("reference", "pre", "anticipation", "post", "post"))
  expect_identical(V$param_index, c(NA_integer_, NA_integer_, 1L, 2L, 3L))
  expect_identical(V$anchored_value, c(0, 0, NA_real_, NA_real_, NA_real_))
  expect_identical(des_2010$incidence$start, c(0L, 0L, 0L, 2L))
  expect_identical(des_2010$incidence$end, c(0L, 1L, 2L, 3L))

  # In Fixture 4, pt_start = 2012 and pt_start = 2015 have c >= b_g = 2012 (D34 global no-extra-window!)
  # and must reduce algebraically to PT-post while keeping anticipation 2015 free!
  links_irreg_no_role <- links_irreg[, setdiff(names(links_irreg), c("s_coordinate_role", "t_coordinate_role"))]
  des_post_irreg <- build_calendar_chained_design(
    links_irreg_no_role, calendar = cal, cohorts = cohorts, anticipation = 2,
    pt_assumption = "post", pt_start = NULL
  )
  for (c_vac in c(2012, 2015)) {
    des_vac <- build_calendar_chained_design(
      links_irreg_no_role, calendar = cal, cohorts = cohorts, anticipation = 2,
      pt_assumption = "all", pt_start = c_vac
    )
    expect_true(des_vac$pt_all_spec$no_extra_calendar_window)
    expect_true(des_vac$pt_all_spec$no_extra_represented_restriction)
    expect_true(des_vac$pt_all_spec$no_extra_restrictions)
    expect_true(des_vac$pt_all_spec$restriction_is_identity)
    expect_identical(des_vac$parameters, des_post_irreg$parameters)
    expect_identical(des_vac$incidence, des_post_irreg$incidence)
  }

  # D33 refusals: missing/NULL/nonfinite/vector/non-numeric/off-grid pt_start under "all",
  # non-NULL pt_start under "post", and unknown pt_assumption
  for (bad_start in list(NULL, NA_real_, NaN, Inf, -Inf, "2010", c(2010, 2012), 2011, 2014)) {
    expect_error(
      validate_calendar_spec(cal, cohorts, anticipation = 2, pt_assumption = "all", pt_start = bad_start),
      class = "cdid_design_calendar_error"
    )
    expect_error(
      build_calendar_chained_design(
        links_irreg_no_role, calendar = cal, cohorts = cohorts, anticipation = 2,
        pt_assumption = "all", pt_start = bad_start
      ),
      class = "cdid_design_calendar_error"
    )
  }
  expect_error(
    validate_calendar_spec(cal, cohorts, anticipation = 2, pt_assumption = "post", pt_start = 2010),
    class = "cdid_design_calendar_error"
  )
  expect_error(
    validate_calendar_spec(cal, cohorts, anticipation = 2, pt_assumption = "pre", pt_start = 2010),
    class = "cdid_design_calendar_error"
  )

  # Calendar extra date exists (2 in [2, 3) for g = 4) but is unrepresented in links (1,3), (3,4)
  links_unrep <- data.frame(
    link_id = c("l13", "l34"),
    g = c(4, 4),
    s = c(1, 3),
    t = c(3, 4),
    functional_key = c("NT:{Inf}", "NT:{Inf}"),
    control_pool = c("{Inf}", "{Inf}"),
    s_qualified = c(FALSE, TRUE),
    t_qualified = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_unrep <- build_calendar_chained_design(
    links_unrep, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf), anticipation = 0,
    pt_assumption = "all", pt_start = 2
  )
  expect_identical(des_unrep$pt_all_spec$by_cohort$n_extra_calendar_dates, 1L)
  expect_false(des_unrep$pt_all_spec$no_extra_calendar_window)
  expect_identical(des_unrep$pt_all_spec$by_cohort$n_extra_represented_dates, 0L)
  expect_identical(des_unrep$pt_all_spec$n_restricted_pre_parameters, 0L)
  expect_true(des_unrep$pt_all_spec$no_extra_represented_restriction)
  expect_true(des_unrep$pt_all_spec$no_extra_restrictions)
  expect_true(des_unrep$pt_all_spec$restriction_is_identity)

  # Partial no-extra-window: g = 2 (no extra window at c = 2; u = 2 is post ATT) and g = 4 (u = 2 is anchored pre)
  links_part <- data.frame(
    link_id = c("g2_12", "g4_23", "g4_34"),
    g = c(2, 4, 4),
    s = c(1, 2, 3),
    t = c(2, 3, 4),
    functional_key = rep("NT:{Inf}", 3L),
    control_pool = rep("{Inf}", 3L),
    s_qualified = c(TRUE, FALSE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_part <- build_calendar_chained_design(
    links_part, calendar = c(1, 2, 3, 4), cohorts = c(2, 4, Inf), anticipation = 0,
    pt_assumption = "all", pt_start = 2
  )
  expect_identical(des_part$pt_all_spec$by_cohort$g, c(2, 4))
  expect_identical(des_part$pt_all_spec$by_cohort$no_extra_calendar_window, c(TRUE, FALSE))
  expect_identical(des_part$pt_all_spec$by_cohort$no_extra_restrictions, c(TRUE, FALSE))
  v_g2_u2 <- des_part$vertices[des_part$vertices$g == 2 & des_part$vertices$time == 2, , drop = FALSE]
  v_g4_u2 <- des_part$vertices[des_part$vertices$g == 4 & des_part$vertices$time == 2, , drop = FALSE]
  expect_identical(v_g2_u2$role, "post")
  expect_false(is.na(v_g2_u2$param_index))
  expect_true(is.na(v_g2_u2$anchored_value))
  expect_identical(v_g4_u2$role, "pre")
  expect_true(is.na(v_g4_u2$param_index))
  expect_identical(v_g4_u2$anchored_value, 0)
})

test_that("Fixture 5: newly identified post ATT from rank-deficient post design is estimable under PT-all without unrestricted post fit", {
  # Accepted Fixture 5 (`spec/pt_all_reference.py`):
  # Links l1 = (1, 3) and l2 = (2, 4) for g = 4, b_4 = 3, calendar = {1, 2, 3, 4}.
  links_f5 <- data.frame(
    link_id = c("l13", "l24"),
    g = c(4, 4),
    s = c(1, 2),
    t = c(3, 4),
    functional_key = c("NYT:{Inf}", "NYT:{Inf}"),
    control_pool = c("{Inf}", "{Inf}"),
    s_qualified = c(FALSE, FALSE),
    t_qualified = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  cal <- c(1, 2, 3, 4)
  cohorts <- c(4, Inf)

  des_post <- build_calendar_chained_design(
    links_f5, calendar = cal, cohorts = cohorts, anticipation = 0,
    pt_assumption = "post", pt_start = NULL
  )
  des_all <- build_calendar_chained_design(
    links_f5, calendar = cal, cohorts = cohorts, anticipation = 0,
    pt_assumption = "all", pt_start = 2
  )

  # Under post: rank = 2 < 3 = n_parameters (nullity = 1); component {2,4} is unanchored; ATT(4,4) is disconnected_component
  expect_identical(des_post$rank_diagnostics$rank, 2L)
  expect_identical(des_post$rank_diagnostics$n_parameters, 3L)
  expect_identical(des_post$rank_diagnostics$nullity, 1L)
  expect_false(des_post$post_targets$level_estimable)
  expect_identical(des_post$post_targets$status, "disconnected_component")

  # Under PT-all (c = 2): vertex (4,2) is anchored, anchoring component {2,4}!
  # Both components ({1,3} and {2,4}) now contain an anchor, so nullity = 0 and rank = 2 == n_parameters!
  expect_identical(des_all$rank_diagnostics$rank, 2L)
  expect_identical(des_all$rank_diagnostics$n_parameters, 2L)
  expect_identical(des_all$rank_diagnostics$nullity, 0L)
  expect_identical(des_all$rank_diagnostics$unanchored_components, integer(0L))
  expect_true(des_all$post_targets$level_estimable)
  expect_identical(des_all$post_targets$status, "estimable")

  # Witness for ATT(4,4) uses link 2 (+1) terminating at anchored pre vertex (4,2)
  wit_44 <- check_level_estimable(des_all, des_all$post_targets$param_index[1L])
  expect_true(wit_44$estimable)
  expect_identical(wit_44$witness_links$link_index, 2L)
  expect_identical(wit_44$witness_links$sign, 1L)
  W_all <- materialize_chained_design(des_all, 100L)$W
  expect_identical(unname(W_all), rbind(c(-1L, 0L), c(0L, 1L)))
  expect_identical(as.double(crossprod(W_all, c(0, 1))), c(0, 1))

  # Solve GMM on des_post vs des_all: under des_post, ATT(4,4) is non-estimable (NA),
  # whereas under des_all, ATT(4,4) is estimable without any unrestricted post fit!
  Omega <- matrix(c(3, 1, 1, 4), nrow = 2L)
  m_vec <- c(-2, 5)

  fit_post <- run_backend_gmm_solve(des_post, m_vec, Omega, weighting = "optimal", materialize_target_if = FALSE)
  expect_false(fit_post$targets$estimable)
  expect_true(is.na(fit_post$targets$estimate))

  fit_all <- run_backend_gmm_solve(des_all, m_vec, Omega, weighting = "optimal")
  expect_true(fit_all$targets$estimable)
  expect_equal(unname(fit_all$internal_gauge$theta_gauge_coords), c(2, 5), tolerance = 1e-13)
  expect_equal(unname(fit_all$targets$estimate), 5, tolerance = 1e-13)
  expect_equal(unname(fit_all$targets$se_sandwich^2), 4, tolerance = 1e-13)

  # Pipeline ledger support & core target plan recognize ATT(4,4) as available under des_all
  fake_meta <- list(
    catalog = list(
      rows = data.frame(g = c(4, 4), s = c(1, 2), t = c(3, 4), stringsAsFactors = FALSE),
      target_qualification = data.frame(
        g = 4, t = 4, ref = 3, role = "post", qualified = TRUE,
        has_eligible_controls = TRUE, stringsAsFactors = FALSE
      )
    )
  )
  ledger_init <- ipw_pipeline_request_ledger(
    target_cohorts = 4, calendar = cal, cohort_universe = cohorts,
    anticipation = 0, control_group = "nevertreated"
  )
  led_post <- ipw_pipeline_ledger_support(ledger_init, fake_meta, des_post)
  led_all <- ipw_pipeline_ledger_support(ledger_init, fake_meta, des_all)
  expect_false(led_post$available)
  expect_identical(led_post$support_status, "disconnected_component")
  expect_true(led_all$available)
  expect_identical(led_all$support_status, "estimable")
  plan_all <- ipw_pipeline_core_target_plan(led_all, des_all)
  expect_identical(plan_all$J_est_expected, 2L)
  expect_identical(plan_all$n_basis, 0L)
  expect_identical(plan_all$core_targets$ledger_index, c(NA_integer_, 1L))
})

test_that("multi-anchor and all-anchor components preserve graph connectivity, non-chaining witnesses, and unanchored post contrast basis", {
  # 1. Multi-anchor component: links (1,2), (2,3), (2,4), (3,5) with g=4, b_4=3, c=2
  links_ma <- data.frame(
    link_id = c("l12", "l23", "l24", "l35"),
    g = rep(4, 4L),
    s = c(1, 2, 2, 3),
    t = c(2, 3, 4, 5),
    functional_key = rep("NYT:{Inf}", 4L),
    control_pool = rep("{Inf}", 4L),
    s_qualified = c(FALSE, FALSE, FALSE, TRUE),
    t_qualified = c(FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_ma <- build_calendar_chained_design(
    links_ma, calendar = c(1, 2, 3, 4, 5), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 2
  )

  # Graph has 5 distinct vertices in 1 connected component with 2 anchors ((4,3) and (4,2))
  expect_identical(nrow(des_ma$vertices), 5L)
  expect_identical(des_ma$forest$n_components, 1L)
  expect_identical(des_ma$rank_diagnostics$nullity, 0L)
  expect_identical(des_ma$rank_diagnostics$unanchored_components, integer(0L))
  expect_identical(des_ma$rank_diagnostics$rank, 3L)
  expect_identical(des_ma$rank_diagnostics$n_parameters, 3L)

  # Both anchors have parent_vertex == 0L and 0-row root_path_certificate
  anchor_rows <- which(is.na(des_ma$vertices$param_index))
  expect_identical(sort(des_ma$vertices$time[anchor_rows]), c(2, 3))
  expect_identical(des_ma$forest$parent_vertex[anchor_rows], c(0L, 0L))
  for (v in anchor_rows) {
    expect_identical(nrow(root_path_certificate(des_ma, v)), 0L)
  }

  # Free vertex (4,4) connects to anchor (4,2) via link 3 (l24, 1 hop), NOT chaining through (4,3)
  v_44 <- which(des_ma$vertices$time == 4)
  cert_44 <- root_path_certificate(des_ma, v_44)
  expect_identical(cert_44$link_index, 3L)
  expect_identical(cert_44$sign, 1L)

  # Every free parameter j has a witness satisfying W_all' w = e_j
  W_ma <- materialize_chained_design(des_ma, 100L)$W
  free_rows <- which(!is.na(des_ma$vertices$param_index))
  for (v in free_rows) {
    p_idx <- des_ma$vertices$param_index[v]
    cert_v <- root_path_certificate(des_ma, v)
    w_v <- numeric(nrow(W_ma))
    w_v[cert_v$link_index] <- cert_v$sign
    e_j <- numeric(ncol(W_ma))
    e_j[p_idx] <- 1.0
    expect_identical(as.double(crossprod(W_ma, w_v)), e_j)
  }

  # 2. All-anchor component {1,2,3} + disconnected post component {4,5} under c = 1
  links_dc <- data.frame(
    link_id = c("l12", "l23", "l45"),
    g = c(4, 4, 4),
    s = c(1, 2, 4),
    t = c(2, 3, 5),
    functional_key = rep("NYT:{Inf}", 3L),
    control_pool = rep("{Inf}", 3L),
    s_qualified = c(FALSE, FALSE, TRUE),
    t_qualified = c(FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_dc <- build_calendar_chained_design(
    links_dc, calendar = c(1, 2, 3, 4, 5), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 1
  )
  expect_identical(des_dc$forest$n_components, 2L)
  expect_identical(des_dc$rank_diagnostics$n_parameters, 2L)
  expect_identical(des_dc$rank_diagnostics$rank, 1L)
  expect_identical(des_dc$rank_diagnostics$nullity, 1L)
  expect_identical(des_dc$post_targets$status, c("disconnected_component", "disconnected_component"))
  mat_dc <- materialize_chained_design(des_dc, 100L)
  expect_equal(unname(mat_dc$rank_diagnostics$nullspace_basis), matrix(c(1L, 1L), nrow = 2L))

  # Contrast ATT(4,5) - ATT(4,4) (params 1 and 2) is estimable while individual levels are not
  c_chk <- check_contrast_estimable(des_dc, 1L, 2L)
  expect_true(c_chk$estimable)
  l_chk <- check_level_estimable(des_dc, 1L)
  expect_false(l_chk$estimable)

  # Pipeline ledger support + post_basis + core_target_plan work with shifted free indices
  fake_meta_dc <- list(
    catalog = list(
      rows = data.frame(g = c(4, 4, 4), s = c(1, 2, 4), t = c(2, 3, 5), stringsAsFactors = FALSE),
      target_qualification = data.frame(
        g = c(4, 4), t = c(4, 5), ref = c(3, 3), role = c("post", "post"),
        qualified = c(TRUE, TRUE), has_eligible_controls = c(TRUE, TRUE),
        stringsAsFactors = FALSE
      )
    )
  )
  led_dc <- ipw_pipeline_request_ledger(
    target_cohorts = 4, calendar = c(1, 2, 3, 4, 5), cohort_universe = c(4, Inf),
    anticipation = 0, control_group = "nevertreated"
  )
  led_dc <- ipw_pipeline_ledger_support(led_dc, fake_meta_dc, des_dc)
  expect_identical(led_dc$support_status, c("disconnected_component", "disconnected_component"))
  plan_dc <- ipw_pipeline_core_target_plan(led_dc, des_dc)
  expect_identical(plan_dc$n_basis, 1L)
  expect_identical(plan_dc$basis$basis_id, "contrast:ATT(g=4,t=5)-ATT(g=4,t=4)")
  expect_identical(plan_dc$J_est_expected, 1L)

  # Solve GMM for the estimable post contrast while keeping zero-loading pre rows (1,2) and (2,3)
  Omega_dc <- matrix(c(
    2, 0, 1,
    0, 2, 0,
    1, 0, 2
  ), nrow = 3L, byrow = TRUE)
  m_dc <- c(2, -1, 4)
  fit_dc <- run_backend_gmm_solve(
    des_dc, m_dc, Omega_dc, weighting = "optimal",
    targets = plan_dc$target_specs
  )
  expect_identical(fit_dc$targets$estimable, c(FALSE, FALSE, TRUE))
  # Contrast l45 = 4 - 1/2 * 2 = 3, with variance 2 - 1/2 = 3/2
  expect_equal(unname(fit_dc$targets$estimate[3L]), 3.0, tolerance = 1e-13)
  expect_equal(unname(fit_dc$targets$se_sandwich[3L]^2), 3 / 2, tolerance = 1e-13)
})

test_that("D28 full vs reduced Tikhonov criterion, structural class certificates, and forged-map refusals hold under PT-all", {
  # 1. D28 full vs reduced Tikhonov equivalence with zero-loading rows:
  #    4-link shared-class design: (1,2), (1,3), (2,3), (3,4) for g = 4, b_4 = 3.
  links_full <- data.frame(
    link_id = c("l12", "l13", "l23", "l34"),
    g = rep(4, 4L),
    s = c(1, 1, 2, 3),
    t = c(2, 3, 3, 4),
    functional_key = rep("NYT:{Inf}", 4L),
    control_pool = rep("{Inf}", 4L),
    s_qualified = c(FALSE, FALSE, FALSE, TRUE),
    t_qualified = c(FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  # Canonical BFS spanning forest retains (1,2), (1,3), (3,4) (indices 1, 2, 4) and drops (2,3) (index 3)
  # where l23 = -l12 + l13.
  kept <- c(1L, 2L, 4L)
  links_red <- links_full[kept, , drop = FALSE]
  rownames(links_red) <- NULL

  E_mat <- rbind(
    c( 1, 0, 0),
    c( 0, 1, 0),
    c(-1, 1, 0),
    c( 0, 0, 1)
  )
  G_mat <- crossprod(E_mat)
  # Singular reduced covariance V_r_12_23 (det = 0) and reduced moments m_r_12_23 in the {(1,2), (2,3), (3,4)} basis:
  E_12_23 <- rbind(
    c(1, 0, 0),
    c(1, 1, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  V_r_12_23 <- matrix(c(
    2, 0, 1,
    0, 2, 1,
    1, 1, 1
  ), nrow = 3L, byrow = TRUE)
  m_r_12_23 <- c(1, 1, 2)
  Omega_full <- E_12_23 %*% V_r_12_23 %*% t(E_12_23)
  m_full <- as.double(E_12_23 %*% m_r_12_23)

  V_red <- Omega_full[kept, kept, drop = FALSE]
  m_red <- m_full[kept]
  expect_identical(unname(E_mat %*% V_red %*% t(E_mat)), unname(Omega_full))
  expect_identical(as.double(E_mat %*% m_red), unname(m_full))

  des_full_c2 <- build_calendar_chained_design(
    links_full, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 2
  )
  des_red_c2 <- build_calendar_chained_design(
    links_red, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 2
  )
  # Notice in des_full_c2:
  # link 1 (1,2): start=1, end=0 -> (-1, 0)
  # link 2 (1,3): start=1, end=0 -> (-1, 0)
  # link 3 (2,3): start=0, end=0 -> ( 0, 0)
  # link 4 (3,4): start=0, end=2 -> ( 0, 1)

  tik_strat <- gmm_covariance_strategy(
    method = "tikhonov", condition_trigger = 1e-10,
    ridge_relative = 1 / 4, overrides = "ridge_relative",
    weighting_requested = "optimal", weighting_canonical = "optimal"
  )
  red_full <- list(mode = "none")
  red_cert <- list(mode = "sparse_class_certificate", G = G_mat, n_links_all = 4L)

  targets_both <- list(
    gmm_target_spec("beta_4_1", 1L, 1.0),
    gmm_target_spec("ATT_4_4", 2L, 1.0)
  )

  fit_full_tik <- run_backend_gmm_solve(
    des_full_c2, m_full, Omega_full, weighting = "optimal",
    method = "conventional_chol_qr", targets = targets_both,
    covariance_strategy = tik_strat, reduction = red_full
  )
  fit_red_tik <- run_backend_gmm_solve(
    des_red_c2, m_red, V_red, weighting = "optimal",
    method = "conventional_chol_qr", targets = targets_both,
    covariance_strategy = tik_strat, reduction = red_cert
  )

  # Trace scale tr(Omega_full) = tr(V_red %*% G_mat) = 9, lambda = (1/4) * (9/4) = 9/16
  expect_equal(fit_full_tik$covariance_handling$lambda, 9 / 16, tolerance = 1e-13)
  expect_equal(fit_red_tik$covariance_handling$lambda, 9 / 16, tolerance = 1e-13)

  # Direct dense full-catalog reference: W_full' (Omega_full + (9/16) I_4)^{-1} W_full
  W_full_c2 <- materialize_chained_design(des_full_c2, 100L)$W
  Q_full_tik <- solve(Omega_full + diag(9 / 16, 4L))
  bread_tik <- solve(t(W_full_c2) %*% Q_full_tik %*% W_full_c2)
  A_tik <- bread_tik %*% t(W_full_c2) %*% Q_full_tik
  theta_tik_ref <- as.double(A_tik %*% m_full)
  sand_tik_ref <- A_tik %*% Omega_full %*% t(A_tik)

  expect_equal(unname(fit_full_tik$internal_gauge$theta_gauge_coords), theta_tik_ref, tolerance = 1e-13)
  expect_equal(unname(fit_red_tik$internal_gauge$theta_gauge_coords), theta_tik_ref, tolerance = 1e-13)
  expect_equal(unname(fit_full_tik$sandwich_cov), unname(sand_tik_ref), tolerance = 1e-13)
  expect_equal(unname(fit_red_tik$sandwich_cov), unname(sand_tik_ref), tolerance = 1e-13)
  expect_equal(
    unname(fit_red_tik$internal_gauge$Lambda),
    unname(crossprod(E_mat, fit_full_tik$internal_gauge$Lambda)),
    tolerance = 1e-13
  )

  # Analytic counterexample: naive reduced ridge (V_red + lambda_r I_3)^{-1} without G_mat disagrees with full-catalog ridge
  W_red_c2 <- materialize_chained_design(des_red_c2, 100L)$W
  lam_naive <- (1 / 4) * (sum(diag(V_red)) / 3)
  Q_naive <- solve(V_red + diag(lam_naive, 3L))
  th_naive <- as.double(solve(t(W_red_c2) %*% Q_naive %*% W_red_c2, t(W_red_c2) %*% Q_naive %*% m_red))
  expect_gt(max(abs(th_naive - fit_red_tik$internal_gauge$theta_gauge_coords)), 0.02)

  # Also verify Fixture 2A's 6-link vs 3-link D28 Tikhonov match (`tau_tikhonov_full6_float = 1.1200017499343775`)
  links_6 <- data.frame(
    link_id = c("l12", "l13", "l23", "l14", "l24", "l34"),
    g = rep(4, 6L),
    s = c(1, 1, 2, 1, 2, 3),
    t = c(2, 3, 3, 4, 4, 4),
    functional_key = rep("NYT:{Inf}", 6L),
    control_pool = rep("{Inf}", 6L),
    s_qualified = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    t_qualified = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_6_c1 <- build_calendar_chained_design(
    links_6, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 1
  )
  # Row order matches `links_6`: 1:(1,2), 2:(1,3), 3:(2,3), 4:(1,4), 5:(2,4), 6:(3,4)
  E_6 <- rbind(
    c(1, 0, 0), # (1,2)
    c(1, 1, 0), # (1,3)
    c(0, 1, 0), # (2,3)
    c(1, 1, 1), # (1,4)
    c(0, 1, 1), # (2,4)
    c(0, 0, 1)  # (3,4)
  )
  V_star_2a <- matrix(c(3, 1, 1, 1, 2, 1, 1, 1, 2), nrow = 3L, byrow = TRUE)
  A_chain_2a <- rbind(c(1, -1, 0), c(0, 1, -1), c(0, 0, 1))
  V_chain_2a <- A_chain_2a %*% V_star_2a %*% t(A_chain_2a)
  m_chain_2a <- as.double(A_chain_2a %*% c(1.4, 1.1, 1.0))
  V_full_6 <- E_6 %*% V_chain_2a %*% t(E_6)
  m_full_6 <- as.double(E_6 %*% m_chain_2a)

  tik_default <- gmm_covariance_strategy(
    method = "tikhonov", condition_trigger = 1e-10, ridge_relative = 1e-4,
    weighting_requested = "optimal", weighting_canonical = "optimal"
  )
  fit_6_tik <- run_backend_gmm_solve(
    des_6_c1, m_full_6, V_full_6, weighting = "optimal",
    covariance_strategy = tik_default
  )
  expect_equal(unname(fit_6_tik$targets$estimate), 1.1200017499343775, tolerance = 1e-13)

  # 2. Structural class certificates (`certify_ipw_structural_classes` and `certify_ipw_structural_classes_view`)
  #    using a 4-unit balanced panel with zero GLM fits!
  pan_4 <- data.frame(
    id = rep(1:4, each = 4L),
    time = rep(as.double(1:4), 4L),
    y = as.double(c(
      1, 2, 4, 7,
      2, 3, 5, 9,
      1, 1, 2, 3,
      0, 2, 3, 5
    )),
    g = rep(c(4, 4, Inf, Inf), each = 4L),
    stringsAsFactors = FALSE
  )
  Z_4 <- matrix(1.0, nrow = 4L, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  w_4 <- rep(1.0, 4L)

  # (a) `certify_ipw_structural_classes` on sealed outcome-free nuisance + declarations (0 GLM fits)
  pidx_4 <- construct_observed_pairs(pan_4, "id", "time", "y", "g")
  ctx_4 <- construct_ipw_context(pidx_4, id = 1:4, Z = Z_4, weights = w_4)
  ses_4 <- construct_ipw_shared_session(ctx_4, max_nuisance = 2L, max_links = 10L, max_unit_entries = 1000L)
  glm_ctrl_4 <- validate_ipw_glm_control(stats::glm.control(epsilon = 1e-12, maxit = 25L, trace = FALSE))
  coh_4 <- validate_ipw_cohorts(4, Inf)
  sup_4 <- ipw_cohort_support(ctx_4, coh_4)
  mask_4 <- 1:4
  gi_4 <- match(ctx_4$a_g[mask_4], coh_4$cohorts)
  pr_4 <- ipw_pair_rates(gi_4, w_4, sup_4)
  fw_4 <- w_4 / pr_4$rates[gi_4]
  beta_4 <- stats::setNames(0.0, "(Intercept)")
  lv_4 <- ipw_logit_values(Z_4, beta_4)
  ow_4 <- ipw_obs_weights(fw_4, lv_4$odds, gi_4, sup_4$n)
  st_4 <- list(
    n = 4.0, m_cohorts = 2L, p = 1L, Z_use = Z_4, D_use = as.double(gi_4 == 1L),
    risk_gi = gi_4, risk_w = w_4,
    pair_count = pr_4$pair_count, pair_mass = pr_4$pair_mass, rates = pr_4$rates,
    fit_weights = fw_4, beta = beta_4, odds_use = lv_4$odds, prob_use = lv_4$prob,
    is_t = ow_4$is_t, a_obs = ow_4$a_obs, b_obs = ow_4$b_obs,
    sum_a = ow_4$sum_a, mean_a = ow_4$mean_a, sum_b = ow_4$sum_b, mean_b = ow_4$mean_b,
    fit_info = list(converged = TRUE, boundary = FALSE, iter = 1L, rank = 1L)
  )
  hs_4 <- ipw_nuisance_hessian_stage(st_4, ipw_cohort_labels(coh_4))
  nu_4 <- structure(
    list(
      serial = 1L, g = coh_4$g, K = coh_4$K, cohorts = coh_4$cohorts,
      cohort_labels = ipw_cohort_labels(coh_4), control = glm_ctrl_4,
      defining_pair = list(s = 1.0, t = 2.0, s_calendar_index = 1L, t_calendar_index = 2L),
      support = sup_4,
      mask = list(observed_a_index = mask_4, n_risk_a_units = 4L, n_observed_risk_units = 4L, n_pair_nonobservers = 0L),
      stage = st_4[setdiff(names(st_4), c("Z_use", "D_use"))],
      hessian = hs_4, unit_entries = 32.0, context_ref = ctx_4$ref, session_ref = ses_4$ref
    ),
    class = "cdid_ipw_nuisance"
  )
  nu_4 <- ipw_seal(nu_4, "nuisance")
  ses_4$ref$nuisance_refs[[1L]] <- nu_4$ref
  ses_4$ref$unit_entries <- nu_4$unit_entries

  pairs_4 <- list(c(1, 2), c(1, 3), c(2, 3), c(3, 4))
  decls_4 <- lapply(seq_along(pairs_4), function(k) {
    p <- pairs_4[[k]]
    declare_ipw_shared_link(
      ses_4, nu_4, 4, Inf, p[1L], p[2L], glm_ctrl_4,
      link_id = paste0("link_", k), functional_key = "NYT:{Inf}"
    )
  })
  bounds_4 <- list(max_links = 10L, max_classes = 5L, max_class_dates = 10L, max_E_entries = 20L, max_mask_entries = 100L)

  des_full_post <- build_calendar_chained_design(
    links_full, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "post", pt_start = NULL
  )
  des_full_c1 <- build_calendar_chained_design(
    links_full, calendar = c(1, 2, 3, 4), cohorts = c(4, Inf),
    anticipation = 0, pt_assumption = "all", pt_start = 1
  )

  cert_post <- certify_ipw_structural_classes(ses_4, decls_4, des_full_post, 1:4, "optimal", bounds_4)
  cert_c2 <- certify_ipw_structural_classes(ses_4, decls_4, des_full_c2, 1:4, "optimal", bounds_4)
  cert_c1 <- certify_ipw_structural_classes(ses_4, decls_4, des_full_c1, 1:4, "optimal", bounds_4)

  # Sample telescoping classes and E are identical across post, c=2, and c=1!
  expect_identical(cert_c2[c("retained_index", "dropped_index", "E", "classes")], cert_post[c("retained_index", "dropped_index", "E", "classes")])
  expect_identical(cert_c1[c("retained_index", "dropped_index", "E", "classes")], cert_post[c("retained_index", "dropped_index", "E", "classes")])
  expect_identical(cert_c1$retained_index, c(1L, 2L, 4L))
  expect_identical(cert_c1$dropped_index, 3L)

  # `apply_ipw_reduction_map` reconstructs W_all (including (0L, 0L) zero rows) and moments
  W_c2 <- materialize_chained_design(des_full_c2, 100L)$W
  W_c1 <- materialize_chained_design(des_full_c1, 100L)$W
  expect_identical(apply_ipw_reduction_map(cert_c2, W_c2[cert_c2$retained_index, , drop = FALSE]), unname(W_c2 * 1.0))
  expect_identical(apply_ipw_reduction_map(cert_c1, W_c1[cert_c1$retained_index, , drop = FALSE]), unname(W_c1 * 1.0))
  expect_identical(apply_ipw_reduction_map(cert_c2, unname(m_full[cert_c2$retained_index])), unname(m_full))

  # Company check: reconstruct actual complete IFs, independently of the
  # hand-supplied covariance above. Reusing the sealed nuisance fits no GLM.
  company_ref <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "cp03_company_dense_reference.R"), company_ref)
  literal <- company_ref$cp03_company_reference_check()
  expect_equal(unname(fit_full_tik$internal_gauge$theta_gauge_coords), literal$ridge$point, tolerance = 1e-13)
  expect_equal(unname(fit_red_tik$internal_gauge$Lambda), literal$ridge_reduced$Lambda, tolerance = 1e-13)
  expect_equal(unname(fit_red_tik$sandwich_cov), literal$ridge_reduced$covariance, tolerance = 1e-13)
  fits_4 <- lapply(decls_4, function(d) fit_ipw_link_shared(ses_4, d))
  shared_phi <- materialize_ipw_if(ctx_4, fits_4, 16L, 4L)$phi_matrix
  expect_equal(unname(shared_phi), literal$shared_phi, tolerance = 1e-13)
  shared_m <- vapply(fits_4, function(f) f$estimate, numeric(1L))
  expect_equal(unname(shared_m), c(0, 1, 1, 2), tolerance = 1e-13)
  for (cert in list(cert_c2, cert_c1)) {
    expect_equal(t(apply_ipw_reduction_map(cert, t(shared_phi[, cert$retained_index]))),
                 unname(shared_phi), tolerance = 1e-13)
    expect_equal(apply_ipw_reduction_map(cert, shared_m[cert$retained_index]),
                 unname(shared_m), tolerance = 1e-13)
  }

  # (b) `certify_ipw_structural_classes_view` on compact custom IPW links (0 GLM fits)
  o_ctrl <- ipw_owner_control(
    max_chunk_rows = 4L, max_retained_entries = 10000L,
    max_matrix_entries = 10000L, max_work_ops = 100000L, max_output_entries = 10000L
  )
  pt_ctrl <- ipw_point_view_control(
    max_chunk_rows = 4L, max_observer_units = 10L, max_retained_entries = 10000L,
    max_matrix_entries = 10000L, max_work_ops = 100000L, max_output_entries = 10000L
  )
  panel_owner_4 <- construct_ipw_panel_owner(pan_4, "id", "time", "y", "g", o_ctrl)
  frame_owner_4 <- construct_ipw_frame_owner(panel_owner_4, id = 1:4, Z = Z_4, weights = w_4, control = o_ctrl)
  nv_4 <- construct_ipw_nuisance_view(frame_owner_4, g = 4, K = Inf, s = 1, t = 2, control = o_ctrl, nuisance_id = "nv4")
  custom_cb <- function(y1, y0, D, covariates, i.weights, ...) {
    dy <- y1 - y0
    psi <- (2 * D - 1) * dy
    list(ATT = mean(psi), att.inf.func = psi - mean(psi))
  }
  clinks_4 <- lapply(seq_along(pairs_4), function(k) {
    p <- pairs_4[[k]]
    lv <- construct_ipw_link_view(
      frame_owner_4, s = p[1L], t = p[2L], g = 4, K = Inf,
      nuisance_view = nv_4, control = o_ctrl,
      link_id = paste0("link_", k), functional_key = "NYT:{Inf}",
      require_support = TRUE, verify_nuisance_mask = TRUE
    )
    fit_custom_ipw_link_view(frame_owner_4, lv, pt_ctrl, custom_est_method = custom_cb)
  })

  vcert_post <- certify_ipw_structural_classes_view(clinks_4, des_full_post, 1:4, "optimal", bounds_4, o_ctrl)
  vcert_c2 <- certify_ipw_structural_classes_view(clinks_4, des_full_c2, 1:4, "optimal", bounds_4, o_ctrl)
  vcert_c1 <- certify_ipw_structural_classes_view(clinks_4, des_full_c1, 1:4, "optimal", bounds_4, o_ctrl)
  expect_identical(vcert_c2[c("retained_index", "dropped_index", "E", "classes")], vcert_post[c("retained_index", "dropped_index", "E", "classes")])
  expect_identical(vcert_c1[c("retained_index", "dropped_index", "E", "classes")], vcert_post[c("retained_index", "dropped_index", "E", "classes")])

  custom_chunks <- list()
  replay_custom_ipw_if_view(clinks_4, pt_ctrl, function(batch) {
    custom_chunks[[length(custom_chunks) + 1L]] <<- batch$phi
  })
  custom_phi <- do.call(rbind, custom_chunks)
  expect_equal(unname(custom_phi), literal$custom_phi, tolerance = 1e-13)
  custom_m <- vapply(clinks_4, function(f) f$estimate, numeric(1L))
  expect_equal(unname(custom_m), c(0, 1/2, 1/2, 1), tolerance = 1e-13)
  for (cert in list(vcert_c2, vcert_c1)) {
    expect_equal(t(apply_ipw_reduction_map(cert, t(custom_phi[, cert$retained_index]))),
                 unname(custom_phi), tolerance = 1e-13)
    expect_equal(apply_ipw_reduction_map(cert, custom_m[cert$retained_index]),
                 unname(custom_m), tolerance = 1e-13)
  }

  # Fixed dates remain distinct intrinsic vertices even when both load zero.
  merged_fixed_dates <- des_full_c1
  merged_fixed_dates$links$start_vertex_index[1L] <- which(merged_fixed_dates$vertices$role == "reference")
  expect_error(validate_ipw_design_incidence(merged_fixed_dates), class = "cdid_ipw_class_vertex_error")
  expect_error(certify_ipw_structural_classes(ses_4, decls_4, merged_fixed_dates, 1:4, "optimal", bounds_4),
               class = "cdid_ipw_class_vertex_error")
  expect_error(certify_ipw_structural_classes_view(clinks_4, merged_fixed_dates, 1:4, "optimal", bounds_4, o_ctrl),
               class = "cdid_ipw_class_vertex_error")

  # Regression for the company fix: the old new-zero-row branch accepted
  # two forged `pre` vertices at stored dates 4,5 although b_4=3. No incidence,
  # free graph, rank or parameter claim changes, so those later checks cannot
  # repair a missing upper bound in the anchor predicate.
  after_reference <- build_calendar_chained_design(
    links_full, calendar = as.double(1:5), cohorts = c(4, Inf),
    pt_assumption = "all", pt_start = 2
  )
  nv <- nrow(after_reference$vertices)
  forged <- after_reference$vertices[rep(which(after_reference$vertices$time == 2), 2L), ]
  forged$time <- c(4, 5)
  forged$vertex_index <- nv + 1:2
  forged$vertex_id <- c("forged-pre-4", "forged-pre-5")
  after_reference$vertices <- rbind(after_reference$vertices, forged)
  after_reference$links$start_vertex_index[3L] <- nv + 1L
  after_reference$links$end_vertex_index[3L] <- nv + 2L
  after_reference$links$s[3L] <- 4
  after_reference$links$t[3L] <- 5
  expect_error(build_gmm_gauge(after_reference, "lowest"), class = "cdid_gmm_design_error")

  # 3. Invalid / forged anchors, maps, and gauge refusals
  # (a) Forging (0L, 0L) on a post design is refused by both validate_ipw_design_incidence and build_gmm_gauge
  bad_post_zero <- des_full_post
  bad_post_zero$incidence$start[1L] <- 0L
  bad_post_zero$incidence$end[1L] <- 0L
  expect_error(validate_ipw_design_incidence(bad_post_zero), class = "cdid_ipw_class_vertex_error")
  expect_error(build_gmm_gauge(bad_post_zero, "lowest"), class = "cdid_gmm_design_error")

  # (b) Forging an anchored pre-start vertex (u = 1 < c = 2) under PT-all is refused
  bad_pre_start <- des_full_c2
  v_u1 <- which(bad_pre_start$vertices$time == 1)
  bad_pre_start$vertices$param_index[v_u1] <- NA_integer_
  bad_pre_start$vertices$anchored_value[v_u1] <- 0
  expect_error(validate_ipw_design_incidence(bad_pre_start), class = "cdid_ipw_class_vertex_error")
  expect_error(certify_ipw_structural_classes(ses_4, decls_4, bad_pre_start, 1:4, "optimal", bounds_4), class = "cdid_ipw_class_vertex_error")
  expect_error(certify_ipw_structural_classes_view(clinks_4, bad_pre_start, 1:4, "optimal", bounds_4, o_ctrl), class = "cdid_ipw_class_vertex_error")

  # (c) Leaving an in-window pre vertex (u = 2 >= c = 2) free under PT-all is refused
  bad_in_win <- des_full_c2
  v_u2 <- which(bad_in_win$vertices$time == 2)
  bad_in_win$vertices$param_index[v_u2] <- 2L
  bad_in_win$vertices$anchored_value[v_u2] <- NA_real_
  expect_error(validate_ipw_design_incidence(bad_in_win), class = "cdid_ipw_class_vertex_error")
  expect_error(certify_ipw_structural_classes(ses_4, decls_4, bad_in_win, 1:4, "optimal", bounds_4), class = "cdid_ipw_class_vertex_error")

  # (d) Forging (0L, 0L) on link l12 under c = 2 (where s = 1 < c = 2 is NOT in A_4(2)) is refused by build_gmm_gauge
  bad_gauge_zero <- des_full_c2
  bad_gauge_zero$incidence$start[1L] <- 0L
  bad_gauge_zero$incidence$end[1L] <- 0L
  expect_error(build_gmm_gauge(bad_gauge_zero, "lowest"), class = "cdid_gmm_design_error")

  # (e) Forging a (0L, 0L) self-loop (start_vertex_index == end_vertex_index) is refused by build_gmm_gauge
  bad_loop <- des_full_c1
  bad_loop$links$end_vertex_index[1L] <- bad_loop$links$start_vertex_index[1L]
  bad_loop$links$t[1L] <- bad_loop$links$s[1L]
  expect_error(build_gmm_gauge(bad_loop, "lowest"), class = "cdid_gmm_design_error")

  # (f) Forging a free parameter self-loop (inc$start == inc$end > 0L) is refused by build_gmm_gauge
  bad_free_loop <- des_full_c2
  bad_free_loop$incidence$start[1L] <- 1L
  bad_free_loop$incidence$end[1L] <- 1L
  expect_error(build_gmm_gauge(bad_free_loop, "lowest"), class = "cdid_gmm_design_error")

  # (g) Forging `pre_pt_all_zero` coordinate role under post mode or on u < c is refused
  bad_role_post <- des_full_post
  bad_role_post$links$s_coordinate_role <- c("pre_nuisance", "pre_nuisance", "pre_pt_all_zero", "reference")
  bad_role_post$links$t_coordinate_role <- c("pre_pt_all_zero", "reference", "reference", "post_att")
  expect_error(validate_ipw_design_incidence(bad_role_post), class = "cdid_ipw_class_vertex_error")
  expect_error(certify_ipw_structural_classes(ses_4, decls_4, bad_role_post, 1:4, "optimal", bounds_4), class = "cdid_ipw_design_error")

  bad_role_pre_c <- des_full_c2
  bad_role_pre_c$links$s_coordinate_role <- c("pre_pt_all_zero", "pre_pt_all_zero", "pre_pt_all_zero", "reference")
  bad_role_pre_c$links$t_coordinate_role <- c("pre_pt_all_zero", "reference", "reference", "post_att")
  expect_error(validate_ipw_design_incidence(bad_role_pre_c), class = "cdid_ipw_class_vertex_error")
  expect_error(certify_ipw_structural_classes(ses_4, decls_4, bad_role_pre_c, 1:4, "optimal", bounds_4), class = "cdid_ipw_design_error")
})

test_that("company dense reference checks identity, two-step, and the global H=I fit", {
  ref_env <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "cp03_company_dense_reference.R"), ref_env)
  r <- ref_env$cp03_company_reference_check()
  links <- data.frame(g = rep(4, 3L), s = c(1, 2, 3), t = c(2, 3, 4),
                      functional_key = "NT:{Inf}", control_pool = "{Inf}",
                      s_qualified = c(FALSE, FALSE, TRUE), t_qualified = c(FALSE, TRUE, TRUE))
  post <- build_calendar_chained_design(links, 1:4, c(4, Inf))
  all <- restrict_chained_design_pt_all(post, 2)
  expect_equal(unname(materialize_chained_design(all, 100L)$W), r$WH, tolerance = 0)
  for (weighting in c("identity", "two_step")) {
    fit <- run_backend_gmm_solve(all, r$m, r$V, weighting = weighting,
      targets = list(gmm_target_spec("early", 1L, 1), gmm_target_spec("post", 2L, 1)))
    expected <- if (weighting == "identity") r$identity else r$zero
    expect_equal(unname(fit$targets$estimate), expected$point, tolerance = 1e-13)
    expect_equal(unname(fit$internal_gauge$Lambda), expected$Lambda, tolerance = 1e-13)
    expect_equal(unname(fit$sandwich_cov), expected$covariance, tolerance = 1e-13)
  }
  for (weighting in c("identity", "two_step")) {
    post_fit <- run_backend_gmm_solve(post, r$m, r$V, weighting = weighting)
    for (c in c(3, 4)) {
      vacuous <- restrict_chained_design_pt_all(post, c)
      expect_true(vacuous$pt_all_spec$restriction_is_identity)
      expect_true(vacuous$pt_all_spec$no_extra_calendar_window)
      expect_identical(vacuous$calendar_spec$pt_assumption, "all")
      expect_identical(vacuous$calendar_spec$pt_start, c)
      fit <- run_backend_gmm_solve(vacuous, r$m, r$V, weighting = weighting)
      expect_identical(fit$targets$estimate, post_fit$targets$estimate)
      expect_identical(fit$internal_gauge$Lambda, post_fit$internal_gauge$Lambda)
      expect_identical(fit$if_projection$materialized_target_if, post_fit$if_projection$materialized_target_if)
    }
  }
})

test_that("company same-catalog check uses actual NT, original NYT and fixed-horizon operators", {
  cal <- as.double(1:6)
  cohorts <- c(4, 6, Inf)
  pairs <- t(utils::combn(cal, 2L))
  occupancy <- do.call(rbind, lapply(cohorts, function(g) {
    data.frame(cohort = g, s = pairs[, 1], t = pairs[, 2], count = 1)
  }))
  bounds <- list(max_candidates = 1e6, max_rows = 1e6, max_pool_entries = 1e7,
                 max_scratch = 1e8, max_work = 1e8)
  for (cg in c("nevertreated", "notyettreated")) {
    cat <- build_comparison_catalog(cal, cohorts, 4, 0, occupancy, cg,
                                    "company:CP03", NULL, bounds)
    post <- build_comparison_catalog_design(cat, data.frame(g = 4, t = 4:6),
                                            data.frame(g = 4, time = 4:6))
    all <- restrict_chained_design_pt_all(post, 2)
    expect_identical(all$links$catalog_row, cat$rows$row_id)
    for (nm in c("g", "s", "t", "ref", "functional_key", "control_pool",
                 "s_qualified", "t_qualified", "start_vertex_index", "end_vertex_index")) {
      expect_identical(all$links[[nm]], post$links[[nm]])
    }
    for (nm in c("g", "time", "role", "functional_key", "vertex_index", "vertex_id")) {
      expect_identical(all$vertices[[nm]], post$vertices[[nm]])
    }
    Wpost <- materialize_chained_design(post, 100000L)$W
    Wall <- materialize_chained_design(all, 100000L)$W
    # In this literal calendar the only extra anchored date is 2; date 1 stays
    # operator-specific, date 3 was already the reference.
    free <- which(!(post$parameters$role == "pre" & post$parameters$time == 2))
    expect_equal(unname(Wall), unname(Wpost[, free, drop = FALSE]), tolerance = 0)
    expect_identical(all$rank_diagnostics$rank, qr(Wall)$rank)
    for (p in which(all$parameters$level_estimable)) {
      cert <- check_level_estimable(all, p)$witness_links
      w <- numeric(nrow(Wall)); w[cert$link_index] <- cert$sign
      target <- numeric(ncol(Wall)); target[p] <- 1
      expect_equal(drop(crossprod(Wall, w)), stats::setNames(target, colnames(Wall)), tolerance = 0)
    }
    if (cg == "notyettreated") {
      early <- which(post$links$s == 1 & post$links$t == 2)
      expect_identical(sort(post$links$control_pool[early]), c("{6,Inf}", "{Inf}"))
      original <- early[post$links$control_pool[early] == "{6,Inf}"]
      fixed <- early[post$links$control_pool[early] == "{Inf}"]
      expect_true(cat$rows$original[original])
      expect_false(cat$rows$original[fixed])
      expect_gt(cat$rows$n_horizons[fixed], 0L)
      swapped_operator <- all
      swapped_operator$links$end_vertex_index[original] <- all$links$end_vertex_index[fixed]
      expect_error(validate_ipw_design_incidence(swapped_operator), class = "cdid_ipw_class_vertex_error")
      # The pre-start bridge from the accepted two-operator fixture, now drawn
      # from actual catalog rows. Columns are beta_1(K6Inf), beta_1(KInf), ATT4.
      fk1 <- post$links$functional_key[original]
      fk2 <- post$links$functional_key[fixed]
      row <- function(s, t, fk) which(post$links$s == s & post$links$t == t & post$links$functional_key == fk)
      cols <- c(which(all$parameters$time == 1 & all$parameters$functional_key == fk1),
                which(all$parameters$time == 1 & all$parameters$functional_key == fk2),
                which(all$parameters$role == "post" & all$parameters$time == 4))
      rows <- c(row(1, 2, fk1), row(1, 2, fk2), row(1, 4, fk2), row(2, 3, fk1), row(2, 4, fk1))
      expect_length(rows, 5L); expect_length(cols, 3L)
      expect_equal(unname(Wall[rows, cols]), rbind(c(-1, 0, 0), c(0, -1, 0),
        c(0, -1, 1), c(0, 0, 0), c(0, 0, 1)), tolerance = 0)
    }
  }
})

test_that("company literal matrices check irregular anticipation, multiple anchors and newly identified ATT", {
  ref_env <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "cp03_company_dense_reference.R"), ref_env)
  r <- ref_env$cp03_company_reference_check()
  fixtures <- list(
    list(cal = c(2010, 2012, 2015, 2017, 2020), g = 2017, delta = 2, c = 2010,
         s = c(2010, 2012, 2012, 2017), t = c(2012, 2015, 2017, 2020),
         sq = c(FALSE, TRUE, TRUE, TRUE), tq = rep(TRUE, 4L), W = r$irregular),
    list(cal = as.double(1:5), g = 4, delta = 0, c = 2,
         s = c(1, 2, 2, 3), t = c(2, 3, 4, 5),
         sq = c(FALSE, FALSE, FALSE, TRUE), tq = c(FALSE, TRUE, TRUE, TRUE), W = r$multi_anchor),
    list(cal = as.double(1:4), g = 4, delta = 0, c = 2, s = c(1, 2), t = c(3, 4),
         sq = c(FALSE, FALSE), tq = c(TRUE, TRUE), W = r$newly_all)
  )
  for (fx in fixtures) {
    links <- data.frame(g = fx$g, s = fx$s, t = fx$t, s_qualified = fx$sq, t_qualified = fx$tq,
                        functional_key = "NT:{Inf}", control_pool = "{Inf}")
    post <- build_calendar_chained_design(links, fx$cal, c(fx$g, Inf), fx$delta)
    all <- restrict_chained_design_pt_all(post, fx$c)
    W <- materialize_chained_design(all, 1000L)$W
    expect_equal(unname(W), fx$W, tolerance = 0)
    expect_identical(all$rank_diagnostics$rank, qr(fx$W)$rank)
    expect_true(all(all$forest$parent_vertex[is.na(all$vertices$param_index)] == 0L))
    p <- ncol(W)
    for (j in seq_len(p)) {
      cert <- check_level_estimable(all, j)$witness_links
      w <- numeric(nrow(W)); w[cert$link_index] <- cert$sign
      expect_equal(unname(drop(crossprod(W, w))), diag(p)[, j], tolerance = 0)
    }
    for (j in seq_len(p)) for (k in seq_len(p)) {
      cert <- check_contrast_estimable(all, j, k)$witness_links
      w <- numeric(nrow(W)); w[cert$link_index] <- cert$sign
      expect_equal(unname(drop(crossprod(W, w))), diag(p)[, k] - diag(p)[, j], tolerance = 0)
    }
  }
  # This solve uses only the restricted design. The literal unrestricted
  # two-row, three-column matrix cannot identify its third coordinate.
  expect_identical(qr(r$newly_post)$rank, 2L)
  expect_false(post$post_targets$level_estimable)
  expect_true(all$post_targets$level_estimable)
  fit <- run_backend_gmm_solve(all, c(-2, 5), rbind(c(3, 1), c(1, 4)), weighting = "two_step")
  expect_equal(unname(fit$targets$estimate), 5, tolerance = 1e-13)
  expect_equal(unname(fit$targets$se_sandwich^2), 4, tolerance = 1e-13)
})
