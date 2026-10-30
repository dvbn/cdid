require_cdid_namespace <- function() {
  if (!isNamespaceLoaded("cdid") && !requireNamespace("cdid", quietly = TRUE)) {
    stop("Actual installed 'cdid' namespace is required; no fallback is permitted in packaged tests.", call. = FALSE)
  }
  ns <- asNamespace("cdid")
  if (!isNamespace(ns) || !environmentIsLocked(ns)) {
    stop("'cdid' namespace must be a locked R namespace.", call. = FALSE)
  }
  ns
}

# Require installed namespace immediately at file load (never skip)
require_cdid_namespace()

test_that("private IPW link enclosure entrypoints exist in locked cdid namespace and are unexported", {
  ns <- require_cdid_namespace()
  expect_true(isNamespace(ns))
  expect_true(environmentIsLocked(ns))

  priv_fns <- c(
    "ipw_nuisance_problem_control",
    "validate_ipw_nuisance_problem_control",
    "ipw_link_enclosure_control",
    "validate_ipw_link_enclosure_control",
    "preflight_ipw_nuisance_problem",
    "preflight_ipw_link_enclosure",
    "preflight_ipw_link_if_chunk",
    "construct_certified_ipw_nuisance_problem",
    "validate_certified_ipw_nuisance_problem",
    "enclose_ipw_link",
    "validate_ipw_link_enclosure",
    "enclose_ipw_link_if_chunk",
    "stream_ipw_link_if_enclosure",
    "enclose_ipw_shared_xi_chunk",
    "apply_ipw_shared_xi_to_link"
  )
  exports <- getNamespaceExports("cdid")
  for (fn in priv_fns) {
    expect_true(exists(fn, envir = ns, mode = "function", inherits = FALSE), info = fn)
    expect_true(bindingIsLocked(fn, ns), info = fn)
    expect_identical(environment(get(fn, envir = ns, inherits = FALSE)), ns, info = fn)
    expect_false(fn %in% exports, info = fn)
  }
})

test_that("exact rational 3-cohort 2-link fixture encloses all functionals, skew L, target L_1, 4-term IFs, and shared xi actions", {
  ns <- require_cdid_namespace()
  fixture_path <- testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R")
  expect_true(file.exists(fixture_path))
  ref <- dget(fixture_path)
  f1 <- ref$fixture1_exact_rational

  cert_ctrl <- ns$logit_certificate_control(
    taylor_degree = 18L, max_abs_eta = 40, max_squarings = 16L,
    max_chunk_rows = 4L, max_matrix_entries = 2500L, max_work_ops = 500000L
  )
  prob_ctrl <- ns$ipw_nuisance_problem_control(
    max_chunk_rows = 4L, max_retained_entries = 5000L,
    max_matrix_entries = 2500L, max_work_ops = 500000L
  )
  link_ctrl <- ns$ipw_link_enclosure_control(
    max_chunk_rows = 10L, max_retained_entries = 5000L,
    max_matrix_entries = 2500L, max_work_ops = 500000L, max_output_entries = 5000L
  )

  prob <- ns$construct_certified_ipw_nuisance_problem(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G,
    pair_observed = f1$pair_observed, cohorts = f1$cohorts,
    beta0 = f1$beta_star_float, F = f1$F_float, r = 0.05,
    certificate_control = cert_ctrl, problem_control = prob_ctrl
  )
  expect_true(prob$certified)
  expect_identical(prob$status, "CERTIFIED")
  expect_false(prob$separation_certified)
  expect_false(prob$sampling_certified)
  expect_false(prob$public_gate_adopted)
  ns$validate_certified_ipw_nuisance_problem(prob, require_certified = TRUE)

  # Verify all-A normalizer n = 10 vs |R| = 8 and |O| = 5
  expect_identical(prob$dimensions$n_a_units, 10L)
  expect_identical(prob$dimensions$n_risk_units, 8L)
  expect_identical(prob$dimensions$n_observed_risk_units, 5L)
  expect_identical(prob$dimensions$n_outside_risk_units, 2L)
  expect_identical(prob$dimensions$n_observed_outside_risk_units, 1L)

  # Verify cohort masses, rates, shares, Hessian H, and skew nonsymmetric Jacobian L (all 3 cohorts including target j=1)
  nuis <- prob$nuisance_enclosures
  cert_enc <- prob$certificate$enclosures
  for (j in seq_len(3L)) {
    expect_lte(cert_enc$M_cohort[j, 1L], f1$M_cohort_float[j])
    expect_gte(cert_enc$M_cohort[j, 2L], f1$M_cohort_float[j])
    expect_lte(cert_enc$Pi_cohort[j, 1L], f1$Pi_cohort_float[j])
    expect_gte(cert_enc$Pi_cohort[j, 2L], f1$Pi_cohort_float[j])
    expect_lte(cert_enc$r_cohort[j, 1L], f1$r_cohort_float[j])
    expect_gte(cert_enc$r_cohort[j, 2L], f1$r_cohort_float[j])
    expect_lte(cert_enc$q_cohort[j, 1L], f1$q_cohort_float[j])
    expect_gte(cert_enc$q_cohort[j, 2L], f1$q_cohort_float[j])
    for (k in seq_len(2L)) {
      expect_lte(nuis$L_0[k, j, 1L], f1$L_float[k, j])
      expect_gte(nuis$L_0[k, j, 2L], f1$L_float[k, j])
      expect_lte(nuis$L_star[k, j, 1L], f1$L_float[k, j])
      expect_gte(nuis$L_star[k, j, 2L], f1$L_float[k, j])
    }
  }
  # Target column L_1[1, 1] = -9/40 = -0.225 is strictly negative, never zeroed!
  expect_lt(nuis$L_star[1L, 1L, 2L], -0.20)
  # Skew nonsymmetric L: L[1, 2] != L[2, 1]
  expect_gt(abs(f1$L_float[1L, 2L] - f1$L_float[2L, 1L]), 0.05)

  # Reuse the SAME certified nuisance problem across two distinct outcome links
  R_H <- chol(f1$H_float)
  xi_shared <- ns$enclose_ipw_shared_xi_chunk(
    problem = prob, a_start = 1L, a_end = 10L,
    control = link_ctrl, H_chol_candidate = R_H
  )
  expect_true(xi_shared$certified)
  expect_true(xi_shared$has_xi_beta)

  for (lk_nm in c("link1", "link2")) {
    lk_ref <- f1[[lk_nm]]
    yL <- if (lk_nm == "link1") f1$yL_1_float else f1$yL_2_float
    yR <- if (lk_nm == "link1") f1$yR_1_float else f1$yR_2_float
    # Put NA_real_ and Inf outside O (in non-observers u3, u6, u8 and outside-R u9, u10) to prove they are never read
    yL_masked <- yL; yR_masked <- yR
    yL_masked[c(3L, 6L, 8L, 9L, 10L)] <- c(NA_real_, NaN, Inf, -Inf, NA_real_)
    yR_masked[c(3L, 6L, 8L, 9L, 10L)] <- c(Inf, -Inf, NA_real_, NaN, Inf)

    enc <- ns$enclose_ipw_link(
      problem = prob, id = f1$id, y_left = yL_masked, y_right = yR_masked,
      loading_beta_candidate = lk_ref$lambda_beta_float,
      control = link_ctrl, link_id = lk_nm
    )
    expect_true(enc$certified)
    ns$validate_ipw_link_enclosure(enc, require_certified = TRUE)
    expect_identical(enc$problem_ref, prob$ref)

    es <- enc$enclosures_star
    expect_lte(es$T[1L], lk_ref$T_float); expect_gte(es$T[2L], lk_ref$T_float)
    expect_lte(es$C[1L], lk_ref$C_float); expect_gte(es$C[2L], lk_ref$C_float)
    expect_lte(es$mu[1L], lk_ref$mu_float); expect_gte(es$mu[2L], lk_ref$mu_float)
    expect_lte(es$d_plus[1L], lk_ref$d_plus_float); expect_gte(es$d_plus[2L], lk_ref$d_plus_float)
    expect_lte(es$d_minus[1L], lk_ref$d_minus_float); expect_gte(es$d_minus[2L], lk_ref$d_minus_float)
    expect_equal(es$gamma[1L, ], c(lo = 0.0, hi = 0.0))

    for (k in seq_len(2L)) {
      expect_lte(es$v[k, 1L], lk_ref$v_float[k]); expect_gte(es$v[k, 2L], lk_ref$v_float[k])
      expect_lte(es$G_beta[k, 1L], -lk_ref$v_float[k]); expect_gte(es$G_beta[k, 2L], -lk_ref$v_float[k])
      expect_lte(es$lambda_beta[k, 1L], lk_ref$lambda_beta_float[k])
      expect_gte(es$lambda_beta[k, 2L], lk_ref$lambda_beta_float[k])
    }
    for (j in seq_len(3L)) {
      expect_lte(es$gamma[j, 1L], lk_ref$gamma_float[j]); expect_gte(es$gamma[j, 2L], lk_ref$gamma_float[j])
      expect_lte(es$lambda_r[j, 1L], lk_ref$lambda_r_float[j])
      expect_gte(es$lambda_r[j, 2L], lk_ref$lambda_r_float[j])
    }

    ch <- ns$enclose_ipw_link_if_chunk(
      link_enclosure = enc, a_start = 1L, a_end = 10L, control = link_ctrl,
      phi_candidate = lk_ref$phi_float,
      xi_beta_candidate = lk_ref$xi_beta_float,
      include_components = TRUE, include_shared_xi = TRUE
    )
    app <- ns$apply_ipw_shared_xi_to_link(
      link_enclosure = enc, shared_xi_chunk = xi_shared,
      control = link_ctrl, phi_candidate = lk_ref$phi_float
    )

    for (u in seq_len(10L)) {
      expect_lte(ch$phi_star_intersected[u, 1L], lk_ref$phi_float[u])
      expect_gte(ch$phi_star_intersected[u, 2L], lk_ref$phi_float[u])
      expect_lte(ch$components_star$treated[u, 1L], lk_ref$treated_float[u])
      expect_gte(ch$components_star$treated[u, 2L], lk_ref$treated_float[u])
      expect_lte(ch$components_star$control[u, 1L], lk_ref$control_float[u])
      expect_gte(ch$components_star$control[u, 2L], lk_ref$control_float[u])
      expect_lte(ch$components_star$score[u, 1L], lk_ref$score_float[u])
      expect_gte(ch$components_star$score[u, 2L], lk_ref$score_float[u])
      expect_lte(ch$components_star$rate[u, 1L], lk_ref$rate_float[u])
      expect_gte(ch$components_star$rate[u, 2L], lk_ref$rate_float[u])
      expect_lte(ch$shared_xi$phi_shared_star[u, 1L], lk_ref$phi_float[u])
      expect_gte(ch$shared_xi$phi_shared_star[u, 2L], lk_ref$phi_float[u])
      expect_lte(app$shared_xi$phi_shared_star[u, 1L], lk_ref$phi_float[u])
      expect_gte(app$shared_xi$phi_shared_star[u, 2L], lk_ref$phi_float[u])
    }

    # Outside-R units (u9, u10) must be identically [0, 0] across all 4 terms
    for (u_out in c(9L, 10L)) {
      expect_equal(ch$phi_star[u_out, ], c(lo = 0.0, hi = 0.0))
      expect_equal(ch$components_star$treated[u_out, ], c(lo = 0.0, hi = 0.0))
      expect_equal(ch$components_star$control[u_out, ], c(lo = 0.0, hi = 0.0))
      expect_equal(ch$components_star$score[u_out, ], c(lo = 0.0, hi = 0.0))
      expect_equal(ch$components_star$rate[u_out, ], c(lo = 0.0, hi = 0.0))
      expect_identical(ch$E_root[u_out], 0.0)
    }
    # Risk non-observers (u6, u8) have zero treated/control/score but strictly nonzero rate term
    for (u_no in c(6L, 8L)) {
      expect_equal(ch$components_star$treated[u_no, ], c(lo = 0.0, hi = 0.0))
      expect_equal(ch$components_star$control[u_no, ], c(lo = 0.0, hi = 0.0))
      expect_equal(ch$components_star$score[u_no, ], c(lo = 0.0, hi = 0.0))
      expect_gt(abs(lk_ref$rate_float[u_no]), 0.1)
      expect_False <- (ch$components_star$rate[u_no, 1L] <= 0.0 && ch$components_star$rate[u_no, 2L] >= 0.0)
      expect_false(expect_False)
    }

    # Stream in small chunks (max_chunk_rows = 3) and verify bit-identical rows and zero-sum enclosure
    stream_ctrl <- ns$ipw_link_enclosure_control(
      max_chunk_rows = 3L, max_retained_entries = 5000L,
      max_matrix_entries = 2500L, max_work_ops = 500000L, max_output_entries = 5000L
    )
    streamed_phi <- matrix(NA_real_, nrow = 10L, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    st_info <- ns$stream_ipw_link_if_enclosure(
      link_enclosure = enc, control = stream_ctrl,
      phi_candidate = lk_ref$phi_float,
      callback = function(batch) {
        streamed_phi[batch$a_index, ] <<- batch$phi_star_intersected
      }
    )
    expect_identical(st_info$n_a_units, 10L)
    expect_identical(st_info$n_batches, 4L)
    expect_true(st_info$phi_star_sum_contains_zero)
    expect_identical(streamed_phi, ch$phi_star_intersected)
  }
})

test_that("nonzero-beta root-ball fixture validates exponential transfer, active target L_1/lambda_{r,1}, Theorem A/B, and Corollary sep", {
  ns <- require_cdid_namespace()
  fixture_path <- testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R")
  ref <- dget(fixture_path)
  f2 <- ref$fixture2_nonzero_beta_rootball
  s_ref <- f2$at_beta_star
  z_ref <- f2$at_beta0_perturbed

  cert_ctrl <- ns$logit_certificate_control(
    taylor_degree = 18L, max_abs_eta = 40, max_squarings = 16L,
    max_chunk_rows = 4L, max_matrix_entries = 2500L, max_work_ops = 500000L
  )
  prob_ctrl <- ns$ipw_nuisance_problem_control(
    max_chunk_rows = 4L, max_retained_entries = 5000L,
    max_matrix_entries = 2500L, max_work_ops = 500000L
  )
  link_ctrl <- ns$ipw_link_enclosure_control(
    max_chunk_rows = 10L, max_retained_entries = 5000L,
    max_matrix_entries = 2500L, max_work_ops = 500000L, max_output_entries = 5000L
  )

  prob2 <- ns$construct_certified_ipw_nuisance_problem(
    id = f2$id, Z = f2$Z_float, w = f2$w_float, G = f2$G,
    pair_observed = f2$pair_observed, cohorts = f2$cohorts,
    beta0 = f2$beta0_perturbed_float, F = f2$F_float, r = f2$r_ball,
    certificate_control = cert_ctrl, problem_control = prob_ctrl
  )
  expect_true(prob2$certified)
  expect_gt(prob2$nuisance_enclosures$kappa_star_lower, 0.0)

  # Non-diagonal F: Lemma dual(ii) coordinate dual norm ||F^{-T} e_1||_2 > M_{11}^{-1/2}
  M_exact <- crossprod(f2$F_float)
  expect_gt(prob2$nuisance_enclosures$coord_dual_iv[1L, 1L], 1.0 / sqrt(M_exact[1L, 1L]))
  # beta_star lies strictly inside beta_star_box
  for (k in seq_len(2L)) {
    expect_lte(prob2$nuisance_enclosures$beta_star_box[k, 1L], f2$beta_star_float[k])
    expect_gte(prob2$nuisance_enclosures$beta_star_box[k, 2L], f2$beta_star_float[k])
  }

  enc2 <- ns$enclose_ipw_link(
    problem = prob2, id = f2$id, y_left = f2$yL_float, y_right = f2$yR_float,
    loading_beta_candidate = z_ref$lambda_beta,
    candidates = list(
      treated_mean = z_ref$T, control_mean = z_ref$C, mu = z_ref$mu,
      v = z_ref$v, gamma = z_ref$gamma, L = z_ref$L,
      loading_beta = z_ref$lambda_beta, loading_rate = z_ref$lambda_r
    ),
    control = link_ctrl, link_id = "f2_link"
  )
  expect_true(enc2$certified)

  # Distinguish beta0 evaluation enclosures from beta* root-ball enclosures
  expect_lte(enc2$enclosures_0$C[1L], z_ref$C); expect_gte(enc2$enclosures_0$C[2L], z_ref$C)
  expect_lte(enc2$enclosures_star$C[1L], s_ref$C); expect_gte(enc2$enclosures_star$C[2L], s_ref$C)
  expect_lte(enc2$enclosures_0$mu[1L], z_ref$mu); expect_gte(enc2$enclosures_0$mu[2L], z_ref$mu)
  expect_lte(enc2$enclosures_star$mu[1L], s_ref$mu); expect_gte(enc2$enclosures_star$mu[2L], s_ref$mu)
  # s_ref$mu is NOT inside enclosures_0$mu (proving beta0 != beta* separation!)
  expect_false(enc2$enclosures_0$mu[1L] <= s_ref$mu && s_ref$mu <= enc2$enclosures_0$mu[2L])

  # Active target rate loading lambda_{r,1} != 0 (around -0.295) is enclosed at both beta0 and beta*
  expect_lt(enc2$enclosures_star$lambda_r[1L, 2L], -0.25)
  for (j in seq_len(3L)) {
    expect_lte(enc2$enclosures_0$lambda_r[j, 1L], z_ref$lambda_r[j])
    expect_gte(enc2$enclosures_0$lambda_r[j, 2L], z_ref$lambda_r[j])
    expect_lte(enc2$enclosures_star$lambda_r[j, 1L], s_ref$lambda_r[j])
    expect_gte(enc2$enclosures_star$lambda_r[j, 2L], s_ref$lambda_r[j])
  }

  # Verify Theorem B and Corollary sep error bounds
  expect_lte(abs(z_ref$mu - s_ref$mu), enc2$root_perturbation_bounds$point_mu)
  expect_lte(abs(z_ref$mu - s_ref$mu), enc2$separated_errors$mu$total_error_bound)
  expect_lte(abs(z_ref$C - s_ref$C), enc2$separated_errors$control_mean$total_error_bound)
  for (j in seq_len(3L)) {
    expect_lte(abs(z_ref$gamma[j] - s_ref$gamma[j]), enc2$separated_errors$gamma$total_error_bound[j])
    expect_lte(abs(z_ref$lambda_r[j] - s_ref$lambda_r[j]), enc2$separated_errors$loading_rate$total_error_bound[j])
  }

  ch2 <- ns$enclose_ipw_link_if_chunk(
    link_enclosure = enc2, a_start = 1L, a_end = 10L, control = link_ctrl,
    phi_candidate = z_ref$phi,
    H_chol_candidate = chol(z_ref$H),
    include_components = TRUE, include_shared_xi = TRUE
  )
  # Target non-observer u3 (G=2, S=FALSE) has strictly nonzero rate term because lambda_{r,1} != 0!
  expect_gt(ch2$components_star$rate[3L, 1L], 0.5)
  for (u in seq_len(10L)) {
    expect_lte(ch2$phi_0[u, 1L], z_ref$phi[u]); expect_gte(ch2$phi_0[u, 2L], z_ref$phi[u])
    expect_lte(ch2$phi_star[u, 1L], s_ref$phi[u]); expect_gte(ch2$phi_star[u, 2L], s_ref$phi[u])
    expect_lte(ch2$phi_star_via_root[u, 1L], s_ref$phi[u]); expect_gte(ch2$phi_star_via_root[u, 2L], s_ref$phi[u])
    expect_lte(ch2$phi_star_intersected[u, 1L], s_ref$phi[u]); expect_gte(ch2$phi_star_intersected[u, 2L], s_ref$phi[u])
    expect_lte(abs(z_ref$phi[u] - s_ref$phi[u]), ch2$E_root[u])
    expect_lte(abs(z_ref$phi[u] - s_ref$phi[u]), ch2$separated_row_errors$total_error_bound[u])
    expect_lte(ch2$shared_xi$phi_shared_star[u, 1L], s_ref$phi[u])
    expect_gte(ch2$shared_xi$phi_shared_star[u, 2L], s_ref$phi[u])
  }

  # Verify Lemma dual(iii) kappa_* = (1 - a*rho)*(1 - delta) divisor when F has delta ~ 0.25 and candidate is 0
  F_kappa <- chol(z_ref$H) / sqrt(0.75)
  prob_kappa <- ns$construct_certified_ipw_nuisance_problem(
    id = f2$id, Z = f2$Z_float, w = f2$w_float, G = f2$G,
    pair_observed = f2$pair_observed, cohorts = f2$cohorts,
    beta0 = f2$beta0_perturbed_float, F = F_kappa, r = f2$r_ball,
    certificate_control = cert_ctrl, problem_control = prob_ctrl
  )
  expect_true(prob_kappa$certified)
  enc_kappa <- ns$enclose_ipw_link(
    problem = prob_kappa, id = f2$id, y_left = f2$yL_float, y_right = f2$yR_float,
    loading_beta_candidate = c(0.0, 0.0), control = link_ctrl
  )
  for (k in seq_len(2L)) {
    expect_lte(enc_kappa$enclosures_star$lambda_beta[k, 1L], s_ref$lambda_beta[k])
    expect_gte(enc_kappa$enclosures_star$lambda_beta[k, 2L], s_ref$lambda_beta[k])
  }
})

test_that("exact real endpoint subtraction, outcome shift cancellation, and all-A normalizer n invariance hold", {
  ns <- require_cdid_namespace()
  fixture_path <- testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R")
  ref <- dget(fixture_path)
  f1 <- ref$fixture1_exact_rational
  f3 <- ref$fixture3_real_endpoint_subtraction

  cert_ctrl <- ns$logit_certificate_control(18L, 40, 16L, 4L, 2500L, 500000L)
  prob_ctrl <- ns$ipw_nuisance_problem_control(4L, 5000L, 2500L, 500000L)
  link_ctrl <- ns$ipw_link_enclosure_control(12L, 5000L, 2500L, 500000L, 5000L)

  prob <- ns$construct_certified_ipw_nuisance_problem(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G,
    pair_observed = f1$pair_observed, cohorts = f1$cohorts,
    beta0 = f1$beta_star_float, F = f1$F_float, r = 0.05,
    certificate_control = cert_ctrl, problem_control = prob_ctrl
  )

  # 1. Real endpoint subtraction:
  # (a) Target observers u1, u2: y_R = 1.0, y_L = -2^-53 (exact Delta = 1 + 2^-53 > 1.0, while fl(1.0 - (-2^-53)) == 1.0)
  # (b) Control observers u4 (y_R = 2^53, y_L = -1.0 -> exact Delta_4 = 2^53 + 1) vs
  #     u5, u7 (y_R = 2^53, y_L = +1.0 -> exact Delta_5 = Delta_7 = 2^53 - 1).
  #     In IEEE binary64 RNTE, fl(2^53 - (-1)) == 2^53 == fl(2^53 - (+1)), losing the exact control spread omega = 2.0.
  #     Outward real endpoint subtraction must enclose omega_exact = 2.0 in omega_control!
  yL_sub <- rep(0.0, 10L); yR_sub <- rep(0.0, 10L)
  yL_sub[c(1L, 2L)] <- f3$yL_target_float
  yR_sub[c(1L, 2L)] <- f3$yR_target_float
  yR_sub[c(4L, 5L, 7L)] <- f3$yR_control_float
  yL_sub[4L] <- f3$yL_control_u4_float
  yL_sub[c(5L, 7L)] <- f3$yL_control_u5_u7_float
  enc_sub <- ns$enclose_ipw_link(
    problem = prob, id = f1$id, y_left = yL_sub, y_right = yR_sub,
    loading_beta_candidate = c(0.0, 0.0), control = link_ctrl
  )
  expect_equal(f3$fl_rounded_target_diff, 1.0)
  expect_equal(f3$fl_rounded_control_u4, f3$fl_rounded_control_u5)
  expect_gte(enc_sub$enclosures_star$T[2L], f3$expected_target_hi_min)
  expect_gt(enc_sub$enclosures_star$T[2L], 1.0)
  expect_lte(enc_sub$enclosures_star$omega[1L], f3$exact_omega_control_float)
  expect_gte(enc_sub$enclosures_star$omega[2L], f3$exact_omega_control_float)


  # 2. Exact dyadic shift cancellation: shifting y_left and y_right by +1024 on O preserves Delta_i
  enc_base <- ns$enclose_ipw_link(
    problem = prob, id = f1$id, y_left = f1$yL_1_float, y_right = f1$yR_1_float,
    loading_beta_candidate = f1$link1$lambda_beta_float, control = link_ctrl
  )
  enc_shift <- ns$enclose_ipw_link(
    problem = prob, id = f1$id,
    y_left = f1$yL_1_float + 1024.0, y_right = f1$yR_1_float + 1024.0,
    loading_beta_candidate = f1$link1$lambda_beta_float, control = link_ctrl
  )
  expect_lte(enc_shift$enclosures_star$mu[1L], f1$link1$mu_float)
  expect_gte(enc_shift$enclosures_star$mu[2L], f1$link1$mu_float)

  # 3. All-A normalizer n = |A| vs |R|, |O|, sum(w):
  # Append 2 extra outside-R units (n = 12 instead of 10); scale F by sqrt(10/12) so M = F'F scales with H/n
  id12 <- 1:12
  Z12 <- rbind(f1$Z_float, c(1.0, 0.0), c(1.0, 0.0))
  w12 <- c(f1$w_float, 5.0, 7.0)
  G12 <- c(f1$G, 77, 88)
  S12 <- c(f1$pair_observed, TRUE, FALSE)
  F12 <- f1$F_float * sqrt(10.0 / 12.0)
  prob12 <- ns$construct_certified_ipw_nuisance_problem(
    id = id12, Z = Z12, w = w12, G = G12, pair_observed = S12,
    cohorts = f1$cohorts, beta0 = f1$beta_star_float, F = F12, r = 0.05,
    certificate_control = cert_ctrl, problem_control = prob_ctrl
  )
  enc12 <- ns$enclose_ipw_link(
    problem = prob12, id = id12,
    y_left = c(f1$yL_1_float, NA_real_, NA_real_),
    y_right = c(f1$yR_1_float, NA_real_, NA_real_),
    loading_beta_candidate = f1$link1$lambda_beta_float, control = link_ctrl
  )
  # Point estimate mu is invariant to outside-R units, while d_plus = Sigma_plus / n scales by 10/12
  expect_lte(enc12$enclosures_star$mu[1L], f1$link1$mu_float)
  expect_gte(enc12$enclosures_star$mu[2L], f1$link1$mu_float)
  expect_lte(enc12$enclosures_star$d_plus[1L], f1$link1$d_plus_float * (10.0 / 12.0))
  expect_gte(enc12$enclosures_star$d_plus[2L], f1$link1$d_plus_float * (10.0 / 12.0))
  # Normalized IF rows phi_i scale by n_new / n_old = 12 / 10
  ch12 <- ns$enclose_ipw_link_if_chunk(enc12, 1L, 12L, link_ctrl)
  for (u in seq_len(10L)) {
    expected_phi12 <- f1$link1$phi_float[u] * (12.0 / 10.0)
    expect_lte(ch12$phi_star_intersected[u, 1L], expected_phi12)
    expect_gte(ch12$phi_star_intersected[u, 2L], expected_phi12)
  }
  expect_equal(ch12$phi_star[11L, ], c(lo = 0.0, hi = 0.0))
  expect_equal(ch12$phi_star[12L, ], c(lo = 0.0, hi = 0.0))
})

test_that("immutability, structural sealing, memory preflight, refused domains, and no global leaks hold", {
  ns <- require_cdid_namespace()
  fixture_path <- testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R")
  ref <- dget(fixture_path)
  f1 <- ref$fixture1_exact_rational

  global_before <- ls(envir = globalenv(), all.names = TRUE)
  opts_before <- options()

  cert_ctrl <- ns$logit_certificate_control(18L, 40, 16L, 4L, 2500L, 500000L)
  prob_ctrl <- ns$ipw_nuisance_problem_control(4L, 5000L, 2500L, 500000L)
  link_ctrl <- ns$ipw_link_enclosure_control(10L, 5000L, 2500L, 500000L, 5000L)

  id_mut <- f1$id
  Z_mut <- f1$Z_float
  w_mut <- f1$w_float
  yL_mut <- f1$yL_1_float
  yR_mut <- f1$yR_1_float

  prob <- ns$construct_certified_ipw_nuisance_problem(
    id = id_mut, Z = Z_mut, w = w_mut, G = f1$G,
    pair_observed = f1$pair_observed, cohorts = f1$cohorts,
    beta0 = f1$beta_star_float, F = f1$F_float, r = 0.05,
    certificate_control = cert_ctrl, problem_control = prob_ctrl
  )
  enc <- ns$enclose_ipw_link(
    problem = prob, id = id_mut, y_left = yL_mut, y_right = yR_mut,
    loading_beta_candidate = f1$link1$lambda_beta_float, control = link_ctrl
  )

  # Mutating caller variables after construction does not affect sealed problem or link
  id_mut[1L] <- 9999L
  Z_mut[1L, 1L] <- 9999.0
  w_mut[1L] <- 9999.0
  yL_mut[1L] <- 9999.0
  yR_mut[1L] <- -9999.0
  expect_identical(prob$raw$id[1L], 1L)
  expect_identical(prob$raw$Z[1L, 1L], 1.0)
  expect_identical(enc$endpoints$y_left_obs[1L], 0.5)
  ns$validate_certified_ipw_nuisance_problem(prob, require_certified = TRUE)
  ns$validate_ipw_link_enclosure(enc, require_certified = TRUE)

  # Tampering with a copy of prob or enc is refused by structural seal checks
  prob_tampered <- prob
  prob_tampered$raw$w[1L] <- 100.0
  expect_error(
    ns$validate_certified_ipw_nuisance_problem(prob_tampered),
    class = "cdid_ipw_enclosure_integrity_error"
  )
  expect_error(
    ns$enclose_ipw_link(prob_tampered, f1$id, f1$yL_1_float, f1$yR_1_float, c(0, 0), link_ctrl),
    class = "cdid_ipw_enclosure_integrity_error"
  )
  enc_tampered <- enc
  enc_tampered$enclosures_star$mu <- c(-100.0, 100.0)
  expect_error(
    ns$validate_ipw_link_enclosure(enc_tampered),
    class = "cdid_ipw_enclosure_integrity_error"
  )

  # Decorative/foreign certificate attributes on inputs are refused
  Z_bad <- f1$Z_float
  attr(Z_bad, "certified") <- TRUE
  expect_error(
    ns$construct_certified_ipw_nuisance_problem(
      f1$id, Z_bad, f1$w_float, f1$G, f1$pair_observed, f1$cohorts,
      f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, prob_ctrl
    ),
    class = "cdid_ipw_enclosure_foreign_certificate_error"
  )

  # Reordered or duplicated IDs, mismatched masks, or nonfinite endpoints on O are refused
  expect_error(
    ns$enclose_ipw_link(prob, rev(f1$id), f1$yL_1_float, f1$yR_1_float, c(0, 0), link_ctrl),
    class = "cdid_ipw_enclosure_id_error"
  )
  expect_error(
    ns$enclose_ipw_link(prob, f1$id, f1$yL_1_float, f1$yR_1_float, c(0, 0), link_ctrl, risk_mask = rep(TRUE, 10L)),
    class = "cdid_ipw_enclosure_mask_error"
  )
  yL_bad_obs <- f1$yL_1_float
  yL_bad_obs[1L] <- NA_real_
  expect_error(
    ns$enclose_ipw_link(prob, f1$id, yL_bad_obs, f1$yR_1_float, c(0, 0), link_ctrl),
    class = "cdid_ipw_enclosure_endpoint_error"
  )

  # Over-budget requests are refused in preflight BEFORE content scans
  tiny_prob_ctrl <- ns$ipw_nuisance_problem_control(4L, 10L, 2500L, 500000L)
  expect_error(
    ns$construct_certified_ipw_nuisance_problem(
      f1$id, f1$Z_float, f1$w_float, f1$G, f1$pair_observed, f1$cohorts,
      f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, tiny_prob_ctrl
    ),
    class = "cdid_ipw_enclosure_bound_error"
  )
  small_chunk_ctrl <- ns$ipw_link_enclosure_control(3L, 5000L, 2500L, 500000L, 5000L)
  expect_error(
    ns$enclose_ipw_link_if_chunk(enc, 1L, 10L, small_chunk_ctrl),
    class = "cdid_ipw_enclosure_bound_error"
  )

  # Domain/contraction failure returns explicit INCONCLUSIVE, never separation/non-identification
  prob_inc <- ns$construct_certified_ipw_nuisance_problem(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G,
    pair_observed = f1$pair_observed, cohorts = f1$cohorts,
    beta0 = c(3.0, -3.0), F = f1$F_float, r = 0.01,
    certificate_control = cert_ctrl, problem_control = prob_ctrl
  )
  expect_false(prob_inc$certified)
  expect_identical(prob_inc$status, "INCONCLUSIVE")
  expect_false(prob_inc$separation_certified)
  expect_false(prob_inc$sampling_certified)
  expect_false(prob_inc$public_gate_adopted)

  enc_inc <- ns$enclose_ipw_link(prob_inc, f1$id, f1$yL_1_float, f1$yR_1_float, c(0, 0), link_ctrl)
  expect_false(enc_inc$certified)
  expect_identical(enc_inc$status, "INCONCLUSIVE")
  expect_false(enc_inc$separation_certified)

  global_after <- ls(envir = globalenv(), all.names = TRUE)
  expect_identical(global_after, global_before)
  expect_identical(options(), opts_before)
})
