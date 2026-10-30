# Independent testthat suite for private CP02 raw weighted-logit root certificate
# (`R/chained_logit_certificate.R`).

# Test the actual installed private namespace. No source-tree fallback exists.
if (!('cdid' %in% loadedNamespaces()) ||
    !exists('certify_chained_logit_root',asNamespace('cdid'),inherits=FALSE))
  stop('Install/load the candidate cdid namespace before these tests.',call.=FALSE)
cert_ns <- asNamespace('cdid')
cert_private_names <- c('cdid_logit_cert_abort','logit_cert_arith_refuse','logit_cert_require_gmm_core','logit_cert_int_bound','logit_cert_safe_product','logit_cert_count_output','logit_cert_safe_sum','logit_certificate_control','validate_logit_certificate_control','validate_logit_certificate_cohorts','validate_logit_certificate_metric_factor','preflight_logit_certificate','logit_iv_intersect','logit_b64_add','logit_b64_sub','logit_b64_mul','logit_b64_div','logit_iv_exp_scalar','logit_iv_exp','logit_iv_logistic_from_odds','logit_iv_forwardsolve_ft','logit_iv_backsolve_f','logit_iv_whitened_hessian_discrepancy','certify_chained_logit_root','chained_logit_certificate','certify_weighted_logit_root','gmm_iv','gmm_iv_dot','gmm_iv_norm','gmm_next')
for (nm in cert_private_names) assign(nm,get(nm,cert_ns,inherits=FALSE),envir=environment())
fixture_path <- testthat::test_path('fixtures','cp02_logit_certificate_fixtures.R')
if (!identical(digest::digest(file=fixture_path,algo='sha256'),
 '008ac83235b01f77793dc286d6171f8b119eec13776546a97dc1913f1134d759'))
 stop('Packaged raw-logit fixture binding differs.',call.=FALSE)
sys.source(fixture_path,envir=environment())

iv_contains <- function(iv, val) {
  is.numeric(iv) && length(iv) == 2L && all(is.finite(iv)) &&
    iv[1L] <= val && val <= iv[2L]
}

test_that('Actual installed namespace owns raw arithmetic and missing dependencies refuse', {
  expect_identical(environment(certify_chained_logit_root),cert_ns)
  expect_identical(environment(logit_certificate_control),cert_ns)
  expect_identical(environment(logit_iv_exp),cert_ns)
  fx<-make_exact_beta0_zero_fixture()
  ctrl<-logit_certificate_control(18,700,11,4,10000,1000000)
  res<-certify_chained_logit_root(fx$Z,fx$w,fx$G,fx$pair_observed,fx$cohorts,fx$beta0,fx$F,fx$r,ctrl)
  expect_s3_class(res,'cdid_logit_certificate')
  expect_identical(res$provenance$normalizer_n,8.0)
  broken_env<-new.env(parent=baseenv())
  broken_env$cdid_logit_cert_abort<-cdid_logit_cert_abort
  f<-logit_cert_require_gmm_core;environment(f)<-broken_env
  broken_env$logit_cert_require_gmm_core<-f
  expect_error(broken_env$logit_cert_require_gmm_core(),class='cdid_logit_certificate_dependency_error')
})

test_that("Exact rational beta0 = 0 problem encloses raw masses, rates, shares, score, and Hessian", {
  fx <- make_exact_beta0_zero_fixture()
  ctrl <- make_logit_cert_default_control(max_chunk_rows = 3L)
  cert <- certify_chained_logit_root(
    Z = fx$Z, w = fx$w, G = fx$G, pair_observed = fx$pair_observed,
    cohorts = fx$cohorts, beta0 = fx$beta0, F = fx$F, r = fx$r, control = ctrl
  )

  expect_identical(cert$dimensions$n_a_units, 8L)
  expect_identical(cert$dimensions$n_risk_units, 7L)
  expect_identical(cert$dimensions$n_outside_risk_units, 1L)
  expect_identical(cert$dimensions$n_observed_risk_units, 5L)
  expect_identical(cert$dimensions$n_nonobserved_risk_units, 2L)
  expect_identical(cert$dimensions$n_observed_outside_risk_units, 1L)
  expect_identical(as.integer(cert$dimensions$by_cohort_risk_count), c(3L, 3L, 1L))
  expect_identical(as.integer(cert$dimensions$by_cohort_observed_count), c(2L, 2L, 1L))

  # Raw cohort masses M_j = (7, 6, 3) and observer masses Pi_j = (3, 3, 3)
  for (j in 1:3) {
    m_exact <- fx$exact$M_num[j] / fx$exact$M_den[j]
    pi_exact <- fx$exact$Pi_num[j] / fx$exact$Pi_den[j]
    r_exact <- fx$exact$r_num[j] / fx$exact$r_den[j]
    q_exact <- fx$exact$q_num[j] / fx$exact$q_den[j]
    mop_exact <- fx$exact$M_num[j] / fx$exact$Pi_num[j]

    expect_true(iv_contains(cert$enclosures$M_cohort[j, ], m_exact))
    expect_true(iv_contains(cert$enclosures$Pi_cohort[j, ], pi_exact))
    expect_true(iv_contains(cert$enclosures$r_cohort[j, ], r_exact))
    expect_true(iv_contains(cert$enclosures$q_cohort[j, ], q_exact))
    expect_true(iv_contains(cert$enclosures$M_over_Pi_cohort[j, ], mop_exact))
  }

  # Non-dyadic rate r_1 = 3/7 has strictly positive interval width enclosing 3/7:
  # 7 * r_cohort[1, 1] <= 3 <= 7 * r_cohort[1, 2]
  expect_lte(7.0 * cert$enclosures$r_cohort[1L, 1L], 3.0)
  expect_gte(7.0 * cert$enclosures$r_cohort[1L, 2L], 3.0)
  expect_lt(cert$enclosures$r_cohort[1L, 1L], cert$enclosures$r_cohort[1L, 2L])

  # Sum of fitting weights over O is M_1 + M_2 + M_3 = 7 + 6 + 3 = 16
  expect_true(iv_contains(cert$enclosures$c_fitting_weight_sum, 16.0))
  # Total A-unit weight sum is 7 + 6 + 3 + 5 = 21
  expect_true(iv_contains(cert$enclosures$w_a_total, 21.0))

  # At beta0 = (0, 0), odds = 1, prob = 0.5, var = 0.25
  expect_true(iv_contains(cert$enclosures$odds_range, 1.0))
  expect_true(iv_contains(cert$enclosures$prob_range, 0.5))
  expect_true(iv_contains(cert$enclosures$var_range, 0.25))

  # Exact score g(0) = (1/8, 1/48)
  expect_true(iv_contains(cert$enclosures$g[1L, ], 0.125))
  expect_true(48.0 * cert$enclosures$g[2L, 1L] <= 1.0 && 1.0 <= 48.0 * cert$enclosures$g[2L, 2L])

  # Exact Hessian H(0) = [[1/2, 29/96], [29/96, 29/96]]
  expect_true(iv_contains(cert$enclosures$H[1L, 1L, ], 0.5))
  expect_true(96.0 * cert$enclosures$H[1L, 2L, 1L] <= 29.0 && 29.0 <= 96.0 * cert$enclosures$H[1L, 2L, 2L])
  expect_true(96.0 * cert$enclosures$H[2L, 1L, 1L] <= 29.0 && 29.0 <= 96.0 * cert$enclosures$H[2L, 1L, 2L])
  expect_true(96.0 * cert$enclosures$H[2L, 2L, 1L] <= 29.0 && 29.0 <= 96.0 * cert$enclosures$H[2L, 2L, 2L])
})

test_that("Wrong normalizer, outside-risk, non-observer, and target-rate mutations are rejected", {
  fx <- make_exact_beta0_zero_fixture()
  ctrl <- make_logit_cert_default_control()
  cert <- certify_chained_logit_root(
    Z = fx$Z, w = fx$w, G = fx$G, pair_observed = fx$pair_observed,
    cohorts = fx$cohorts, beta0 = fx$beta0, F = fx$F, r = fx$r, control = ctrl
  )

  # True g_1(0) = 1 / 8 = 0.125 (normalized by eligible A-unit count n = 8).
  # Wrong normalizers:
  #   |R| = 7 -> 1/7
  #   |O| = 5 -> 1/5
  #   sum_A w = 21 -> 1/21
  #   sum_R w = 16 -> 1/16
  for (wrong_n in c(7.0, 5.0, 21.0, 16.0)) {
    wrong_g1 <- 1.0 / wrong_n
    expect_false(iv_contains(cert$enclosures$g[1L, ], wrong_g1))
    wrong_H11 <- 4.0 / wrong_n
    expect_false(iv_contains(cert$enclosures$H[1L, 1L, ], wrong_H11))
  }

  # Outside-risk unit (unit 8, G = 9, pair_observed = TRUE, z = (1, 100), w = 5):
  # Changing its covariates z_8 or weight w_8 or pair_observed_8 must NOT change M_j, Pi_j, r_j, g, or H,
  # because outside-risk units enter ONLY through the eligible-unit count n!
  Z_mut_out <- fx$Z
  Z_mut_out[8L, ] <- c(-999.0, 9999.0)
  w_mut_out <- fx$w
  w_mut_out[8L] <- 1234.5
  po_mut_out <- fx$pair_observed
  po_mut_out[8L] <- FALSE
  cert_out <- certify_chained_logit_root(
    Z = Z_mut_out, w = w_mut_out, G = fx$G, pair_observed = po_mut_out,
    cohorts = fx$cohorts, beta0 = fx$beta0, F = fx$F, r = fx$r, control = ctrl
  )
  expect_identical(cert_out$enclosures$M_cohort, cert$enclosures$M_cohort)
  expect_identical(cert_out$enclosures$Pi_cohort, cert$enclosures$Pi_cohort)
  expect_identical(cert_out$enclosures$g, cert$enclosures$g)
  expect_identical(cert_out$enclosures$H, cert$enclosures$H)
  expect_identical(cert_out$delta, cert$delta)
  expect_identical(cert_out$t, cert$t)

  # Dropping outside-risk unit 8 (reducing n from 8 to 7) DOES change q_j, g, and H by 8/7
  cert_drop_out <- certify_chained_logit_root(
    Z = fx$Z[-8L, , drop = FALSE], w = fx$w[-8L], G = fx$G[-8L],
    pair_observed = fx$pair_observed[-8L],
    cohorts = fx$cohorts, beta0 = fx$beta0, F = fx$F, r = fx$r, control = ctrl
  )
  expect_false(iv_contains(cert$enclosures$g[1L, ], mean(cert_drop_out$enclosures$g[1L, ])))

  # Pair non-observer in target cohort (unit 3, G = 2, pair_observed = FALSE, w = 4):
  # Mutating its weight w_3 from 4 to 1 changes M_1 from 7 to 4 (and r_1 from 3/7 to 3/4),
  # which changes target fitting weights c_1, c_2 and changes g_1(0) from 1/8 to (4*(-1/2) + 9*(1/2))/8 = 5/16!
  w_mut_nonobs <- fx$w
  w_mut_nonobs[3L] <- 1.0
  cert_nonobs <- certify_chained_logit_root(
    Z = fx$Z, w = w_mut_nonobs, G = fx$G, pair_observed = fx$pair_observed,
    cohorts = fx$cohorts, beta0 = fx$beta0, F = fx$F, r = fx$r, control = ctrl
  )
  expect_true(iv_contains(cert_nonobs$enclosures$g[1L, ], 5.0 / 16.0))
  expect_false(iv_contains(cert$enclosures$g[1L, ], 5.0 / 16.0))

  # Target-rate omission mutation on unsaturated 3-cell fixture:
  # In make_unsaturated_3cell_fixture(), target cohort has r_1 = 1/2 (M_1 = 59, Pi_1 = 29.5),
  # so true target fitting weights are c_i = w_i / r_1 = 2 * w_i = (16, 18, 25),
  # giving exact root beta* = (0, log(2)) where g(beta*) = (0, 0).
  # If a mutant omits target-rate scaling (setting c_i = w_i for target cohort by dropping
  # the target non-observer unit 4), g(beta*) is NO LONGER zero and the certificate at r = 0.01 refuses!
  unsat <- make_unsaturated_3cell_fixture()
  cert_unsat <- certify_chained_logit_root(
    Z = unsat$Z, w = unsat$w, G = unsat$G, pair_observed = unsat$pair_observed,
    cohorts = unsat$cohorts, beta0 = unsat$beta0, F = unsat$F, r = unsat$r, control = ctrl
  )
  expect_true(cert_unsat$certified)
  expect_identical(cert_unsat$status, "CERTIFIED")
  expect_true(iv_contains(cert_unsat$enclosures$g[1L, ], 0.0))
  expect_true(iv_contains(cert_unsat$enclosures$g[2L, ], 0.0))

  # Mutate target non-observer weight so target rate scaling is omitted (w_4 -> 0.0001 instead of 29.5)
  w_bad_target_rate <- unsat$w
  w_bad_target_rate[4L] <- 1e-4
  cert_bad_tr <- certify_chained_logit_root(
    Z = unsat$Z, w = w_bad_target_rate, G = unsat$G, pair_observed = unsat$pair_observed,
    cohorts = unsat$cohorts, beta0 = unsat$beta0, F = unsat$F, r = unsat$r, control = ctrl
  )
  expect_false(cert_bad_tr$certified)
  expect_identical(cert_bad_tr$status, "INCONCLUSIVE")
  expect_false(iv_contains(cert_bad_tr$enclosures$g[1L, ], 0.0))
})

test_that("Full-rank finite-root balanced and unsaturated fixtures are certified and bound true root", {
  ctrl <- make_logit_cert_default_control()

  # 1. 16-unit propagation fixture (nuisance A)
  fx16 <- make_propagation_16unit_fixture()
  cert16 <- certify_chained_logit_root(
    Z = fx16$Z, w = fx16$w, G = fx16$G, pair_observed = fx16$pair_observed,
    cohorts = fx16$cohorts, beta0 = fx16$beta0, F = fx16$F, r = fx16$r, control = ctrl
  )
  expect_true(cert16$certified)
  expect_identical(cert16$status, "CERTIFIED")
  expect_false(cert16$inconclusive)
  expect_false(cert16$separation_certified)
  expect_false(cert16$sampling_certified)
  expect_false(cert16$public_gate_adopted)
  expect_lt(cert16$delta, 1.0)
  expect_lt(cert16$x, 1.0)
  expect_lt(cert16$q, 1.0)
  expect_lte(cert16$selfmap, fx16$r)
  expect_lte(cert16$rho, fx16$r)

  # Verify ||F (beta* - beta0)||_2 <= rho <= r for 16-unit fixture.
  dist_M_16 <- sqrt(sum((fx16$F %*% (fx16$beta_star_hp - fx16$beta0))^2))
  expect_lte(dist_M_16, cert16$rho)
  expect_lte(cert16$rho, cert16$r)

  # 2. Finite-logit overlap fixture OV at non-root stored coefficient u0 = (11/8, 9/8)
  fx_ov <- make_finite_logit_ov_fixture()
  cert_ov <- certify_chained_logit_root(
    Z = fx_ov$Z, w = fx_ov$w, G = fx_ov$G, pair_observed = fx_ov$pair_observed,
    cohorts = fx_ov$cohorts, beta0 = fx_ov$beta0, F = fx_ov$F, r = fx_ov$r, control = ctrl
  )
  expect_true(cert_ov$certified)
  expect_identical(cert_ov$status, "CERTIFIED")
  dist_M_ov <- sqrt(sum((fx_ov$F %*% (fx_ov$beta_star_hp - fx_ov$beta0))^2))
  expect_lte(dist_M_ov, cert_ov$rho)
  expect_lte(cert_ov$rho, cert_ov$r)

  # 3. Unsaturated 3-cell fixture around slightly perturbed beta0
  fx_un <- make_unsaturated_3cell_fixture()
  beta0_pert <- fx_un$beta0 + c(1e-4, -1e-4)
  cert_un <- certify_chained_logit_root(
    Z = fx_un$Z, w = fx_un$w, G = fx_un$G, pair_observed = fx_un$pair_observed,
    cohorts = fx_un$cohorts, beta0 = beta0_pert, F = fx_un$F, r = 0.01, control = ctrl
  )
  expect_true(cert_un$certified)
  dist_M_un <- sqrt(sum((fx_un$F %*% (fx_un$beta0 - beta0_pert))^2))
  expect_lte(dist_M_un, cert_un$rho)
  expect_lte(cert_un$rho, 0.01)
})

test_that("Raw weights vs cached rounded w/r score discrepancy is detected with no false target substitution", {
  fx <- make_raw_vs_cached_wr_fixture()
  ctrl <- make_logit_cert_default_control()

  # (a) Small-integer normal weights (Pi_1 = 5, M_1 = 29):
  # Verify IEEE 754 binary64 discrepancy between cached w_1 / fl(Pi_1 / M_1) and exact w_1 * M_1 / Pi_1:
  expect_identical(fx$c1_exact - fx$c1_cached, 2^(-48))
  # Cached surrogate score at beta0 = 0 is identically 0.0 in binary64, whereas exact raw score is -2^-49 / 3 < 0
  expect_identical(fx$g0_cached, 0.0)
  expect_lt(fx$g0_exact, 0.0)

  cert <- certify_chained_logit_root(
    Z = fx$Z, w = fx$w, G = fx$G, pair_observed = fx$pair_observed,
    cohorts = fx$cohorts, beta0 = fx$beta0, F = fx$F, r = 1e-13, control = ctrl
  )
  expect_true(cert$certified)
  expect_true(iv_contains(cert$enclosures$c_fitting_weight_range_by_cohort[1L, ], 29.0))
  expect_true(iv_contains(cert$enclosures$g[1L, ], fx$g0_exact))

  # True raw-weight root distance ||F (beta*_raw - 0)||_2 is bounded by certified rho:
  dist_raw_root <- abs(fx$F[1L, 1L] * fx$beta_star_raw)
  expect_lte(dist_raw_root, cert$rho)

  # A tiny radius around 0 that would pass if g(0) were exactly 0 is refused:
  cert_tiny_r <- certify_chained_logit_root(
    Z = fx$Z, w = fx$w, G = fx$G, pair_observed = fx$pair_observed,
    cohorts = fx$cohorts, beta0 = fx$beta0, F = fx$F, r = 1e-18, control = ctrl
  )
  expect_false(cert_tiny_r$certified)
  expect_identical(cert_tiny_r$status, "INCONCLUSIVE")
  expect_true("selfmap_gt_r" %in% cert_tiny_r$inconclusive_reasons)

  # (b) Subnormal cached rate case (Pi_1 = 3 * 2^-1040, M_1 = 3 * 2^-1040 + 7):
  # Cached r_1_fl = fl(Pi_1 / M_1) loses 18 bits in subnormal range, so c1_cached and g0_cached
  # lie strictly outside the certified raw-weight outward enclosures c_iv and g_iv!
  sub_fx <- fx$subnormal_rate
  cert_sub <- certify_chained_logit_root(
    Z = sub_fx$Z, w = sub_fx$w, G = sub_fx$G, pair_observed = sub_fx$pair_observed,
    cohorts = sub_fx$cohorts, beta0 = sub_fx$beta0, F = sub_fx$F, r = 1e-10, control = ctrl
  )
  expect_true(cert_sub$certified)
  expect_true(iv_contains(cert_sub$enclosures$c_fitting_weight_range_by_cohort[1L, ], 7.0))
  expect_false(iv_contains(cert_sub$enclosures$c_fitting_weight_range_by_cohort[1L, ], sub_fx$c1_cached))
  expect_false(iv_contains(cert_sub$enclosures$g[1L, ], sub_fx$g0_cached))
})

test_that("Badly approximate root, quasi-separation, and complete separation never produce false certificates", {
  ctrl <- make_logit_cert_default_control()

  # 1. Badly approximate root on well-behaved overlap data
  fx16 <- make_propagation_16unit_fixture()
  cert_bad_b0 <- certify_chained_logit_root(
    Z = fx16$Z, w = fx16$w, G = fx16$G, pair_observed = fx16$pair_observed,
    cohorts = fx16$cohorts, beta0 = c(3.5, -2.5), F = fx16$F, r = 0.05, control = ctrl
  )
  expect_false(cert_bad_b0$certified)
  expect_identical(cert_bad_b0$status, "INCONCLUSIVE")
  expect_true(cert_bad_b0$inconclusive)
  expect_false(cert_bad_b0$separation_certified)

  # 2. Quasi-separated data (QS: no X = 1 control) across slopes and radius grid
  for (slope in c(log(2^10), log(2^20), log(2^30))) {
    fx_qs <- make_quasi_separated_fixture(slope = slope)
    a_approx <- sqrt(max(colSums(forwardsolve(t(fx_qs$F), t(fx_qs$Z))^2)))
    for (xs in c(0.05, 0.15, 0.29, 0.5, 0.9, 1.5)) {
      r_try <- xs / a_approx
      cert_qs <- certify_chained_logit_root(
        Z = fx_qs$Z, w = fx_qs$w, G = fx_qs$G, pair_observed = fx_qs$pair_observed,
        cohorts = fx_qs$cohorts, beta0 = fx_qs$beta0, F = fx_qs$F, r = r_try, control = ctrl
      )
      expect_false(cert_qs$certified)
      expect_identical(cert_qs$status, "INCONCLUSIVE")
      expect_false(cert_qs$separation_certified)
      if (xs >= 1.0) {
        expect_false(cert_qs$inequalities$x_lt_1)
        expect_true("x_ge_1" %in% cert_qs$inconclusive_reasons)
      }
    }
  }

  # Verify the x < 1 guard mutation explicitly: at xs = 1.5 on QS, unguarded formula
  # (delta + x) / (1 - x) is negative (~ -3 < 1) and t + q*r < 0 <= r, which would falsely pass
  # if x < 1 were not enforced!
  fx_qs20 <- make_quasi_separated_fixture(slope = log(2^20))
  a_qs20 <- sqrt(max(colSums(forwardsolve(t(fx_qs20$F), t(fx_qs20$Z))^2)))
  cert_qs_15 <- certify_chained_logit_root(
    Z = fx_qs20$Z, w = fx_qs20$w, G = fx_qs20$G, pair_observed = fx_qs20$pair_observed,
    cohorts = fx_qs20$cohorts, beta0 = fx_qs20$beta0, F = fx_qs20$F,
    r = 1.5 / a_qs20, control = ctrl
  )
  expect_gt(cert_qs_15$x, 1.0)
  expect_false(cert_qs_15$inequalities$x_lt_1)
  expect_false(cert_qs_15$certified)
  unguarded_q <- (cert_qs_15$delta + cert_qs_15$x) / (1.0 - cert_qs_15$x)
  expect_lt(unguarded_q, 1.0) # confirms why the x < 1 guard is essential!

  # 3. Complete separation (D = 1 iff x > 0)
  Z_cs <- cbind(1.0, c(-2.0, -1.0, 1.0, 2.0))
  w_cs <- c(1.0, 1.0, 1.0, 1.0)
  G_cs <- c(5.0, 5.0, 2.0, 2.0)
  po_cs <- c(TRUE, TRUE, TRUE, TRUE)
  b0_cs <- c(0.0, 8.0)
  p_cs <- 1.0 / (1.0 + exp(-drop(Z_cs %*% b0_cs)))
  H_cs <- crossprod(Z_cs, Z_cs * (p_cs * (1.0 - p_cs))) / 4.0
  F_cs <- chol(H_cs)
  F_cs[2L, 1L] <- 0.0
  cert_cs <- certify_chained_logit_root(
    Z = Z_cs, w = w_cs, G = G_cs, pair_observed = po_cs,
    cohorts = c(2.0, 5.0), beta0 = b0_cs, F = F_cs, r = 0.05, control = ctrl
  )
  expect_false(cert_cs$certified)
  expect_identical(cert_cs$status, "INCONCLUSIVE")
  expect_false(cert_cs$separation_certified)
})

test_that("Rank-deficient Z, all-zero Z, and ill-conditioned metric/Hessian return INCONCLUSIVE", {
  ctrl <- make_logit_cert_default_control()

  # 1. Rank-deficient Z (duplicate columns: col 2 == col 1)
  Z_rd <- cbind(c(1.0, 1.0, 1.0, 1.0), c(1.0, 1.0, 1.0, 1.0))
  w_rd <- c(1.0, 1.0, 1.0, 1.0)
  G_rd <- c(2.0, 2.0, 5.0, 5.0)
  po_rd <- c(TRUE, TRUE, TRUE, TRUE)
  F_rd <- diag(c(0.5, 0.5))
  cert_rd <- certify_chained_logit_root(
    Z = Z_rd, w = w_rd, G = G_rd, pair_observed = po_rd,
    cohorts = c(2.0, 5.0), beta0 = c(0.0, 0.0), F = F_rd, r = 0.1, control = ctrl
  )
  expect_false(cert_rd$certified)
  expect_identical(cert_rd$status, "INCONCLUSIVE")
  expect_gte(cert_rd$delta, 1.0)
  expect_false(cert_rd$inequalities$delta_lt_1)
  expect_true("delta_ge_1" %in% cert_rd$inconclusive_reasons)

  # 2. All-zero Z on observed units
  Z_zero <- matrix(0.0, nrow = 4L, ncol = 2L)
  cert_zero <- certify_chained_logit_root(
    Z = Z_zero, w = w_rd, G = G_rd, pair_observed = po_rd,
    cohorts = c(2.0, 5.0), beta0 = c(0.0, 0.0), F = F_rd, r = 0.1, control = ctrl
  )
  expect_false(cert_zero$certified)
  expect_false(cert_zero$inequalities$a_gt_0)
  expect_true("a_not_positive" %in% cert_zero$inconclusive_reasons)

  # 3. Poorly matched metric factor F = 10 * I (making delta >= 1 even when true root exists)
  fx16 <- make_propagation_16unit_fixture()
  cert_bad_F <- certify_chained_logit_root(
    Z = fx16$Z, w = fx16$w, G = fx16$G, pair_observed = fx16$pair_observed,
    cohorts = fx16$cohorts, beta0 = fx16$beta0, F = diag(c(10.0, 10.0)), r = 0.01, control = ctrl
  )
  expect_false(cert_bad_F$certified)
  expect_identical(cert_bad_F$status, "INCONCLUSIVE")
  expect_true("delta_ge_1" %in% cert_bad_F$inconclusive_reasons)
  expect_false(cert_bad_F$separation_certified)
})

test_that("Non-diagonal F metric uses exact F'F and rejects malformed inputs before allocation", {
  ctrl <- make_logit_cert_default_control()
  fx <- make_exact_beta0_zero_fixture()

  # Construct upper-triangular F whose off-diagonal entry creates an exact F'F[2,2]
  # that differs from rounded fl(F'F)[2,2]:
  # F = [[1, 1], [0, 2^-27]] -> exact (F'F)_{22} = 1 + 2^-54, whereas fl(1 + 2^-54) = 1.0!
  F_exact_test <- matrix(c(1.0, 0.0, 1.0, 2^(-27)), nrow = 2L, ncol = 2L)
  H_dummy <- array(0.0, dim = c(2L, 2L, 2L))
  H_dummy[1L, 1L, ] <- c(1.0, 1.0)
  H_dummy[1L, 2L, ] <- c(1.0, 1.0)
  H_dummy[2L, 1L, ] <- c(1.0, 1.0)
  H_dummy[2L, 2L, ] <- c(1.0, 1.0)
  wh <- logit_iv_whitened_hessian_discrepancy(F_exact_test, H_dummy)
  # Exact M_{22} = 1 + 2^-54 > 1, so M_iv[2, 2, 2] > 1.0 and E_{22} = 1 - (1 + 2^-54) = -2^-54 < 0:
  expect_gt(wh$M_iv[2L, 2L, 2L], 1.0)
  expect_lt(wh$E_iv[2L, 2L, 1L], 0.0)

  # Malformed F: lower-triangular nonzero entry
  F_lower <- matrix(c(1.0, 0.1, 0.2, 1.0), 2L, 2L)
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, fx$cohorts, fx$beta0, F_lower, fx$r, ctrl),
    class = "cdid_logit_certificate_metric_error"
  )

  # Malformed F: zero or negative diagonal
  F_neg <- matrix(c(1.0, 0.0, 0.2, 0.0), 2L, 2L)
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, fx$cohorts, fx$beta0, F_neg, fx$r, ctrl),
    class = "cdid_logit_certificate_metric_error"
  )

  # Decorative metadata on F or beta0
  F_dec <- fx$F
  attr(F_dec, "is_pd") <- TRUE
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, fx$cohorts, fx$beta0, F_dec, fx$r, ctrl),
    class = "cdid_logit_certificate_metric_error"
  )
  b0_dec <- fx$beta0
  attr(b0_dec, "converged") <- TRUE
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, fx$cohorts, b0_dec, fx$F, fx$r, ctrl),
    class = "cdid_logit_certificate_beta_error"
  )

  # Nonpositive or NA weights
  w_neg <- fx$w
  w_neg[2L] <- 0.0
  expect_error(
    certify_chained_logit_root(fx$Z, w_neg, fx$G, fx$pair_observed, fx$cohorts, fx$beta0, fx$F, fx$r, ctrl),
    class = "cdid_logit_certificate_weight_error"
  )
  w_neg[2L] <- -1.0
  expect_error(
    certify_chained_logit_root(fx$Z, w_neg, fx$G, fx$pair_observed, fx$cohorts, fx$beta0, fx$F, fx$r, ctrl),
    class = "cdid_logit_certificate_weight_error"
  )

  # Missing cohort support in R or pair support in O
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, c(2.0, 5.0, 77.0), fx$beta0, fx$F, fx$r, ctrl),
    class = "cdid_logit_certificate_cohort_support_error"
  )
  po_no_target <- fx$pair_observed
  po_no_target[fx$G == 2.0] <- FALSE
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, po_no_target, fx$cohorts, fx$beta0, fx$F, fx$r, ctrl),
    class = "cdid_logit_certificate_pair_support_error"
  )

  # Target cohort cannot be +Inf, and cohorts cannot have duplicates
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, c(Inf, 5.0), fx$beta0, fx$F, fx$r, ctrl),
    class = "cdid_logit_certificate_cohort_error"
  )
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, c(2.0, 2.0), fx$beta0, fx$F, fx$r, ctrl),
    class = "cdid_logit_certificate_cohort_error"
  )

  # Nonpositive radius r
  expect_error(
    certify_chained_logit_root(fx$Z, fx$w, fx$G, fx$pair_observed, fx$cohorts, fx$beta0, fx$F, 0.0, ctrl),
    class = "cdid_logit_certificate_radius_error"
  )

  # Missing control fields or taylor_degree < 2
  expect_error(
    logit_certificate_control(taylor_degree = 1L, max_abs_eta = 700, max_squarings = 11L,
                              max_chunk_rows = 10L, max_matrix_entries = 1000L, max_work_ops = 10000L),
    class = "cdid_logit_certificate_control_error"
  )
  expect_error(
    validate_logit_certificate_control(list(taylor_degree = 18L)),
    class = "cdid_logit_certificate_control_error"
  )
})

test_that("Validated Lemma exp handles signed zero, subnormals, exact dyadic boundaries, large N, and refuses overflow/underflow", {
  # 1. Signed zero +0.0 and -0.0
  e_pos0 <- logit_iv_exp_scalar(0.0, taylor_degree = 18L, max_abs_eta = 700.0, max_squarings = 11L)
  e_neg0 <- logit_iv_exp_scalar(-0.0, taylor_degree = 18L, max_abs_eta = 700.0, max_squarings = 11L)
  expect_identical(e_pos0$iv, e_neg0$iv)
  expect_identical(e_pos0$m, 0L)
  expect_true(iv_contains(e_pos0$iv, 1.0))
  expect_lt(e_pos0$iv[1L], 1.0)
  expect_gt(e_pos0$iv[2L], 1.0)

  # 2. Subnormal inputs +/- 2^-1074 and +/- 2^-1022
  for (sub_y in c(2^(-1074), -2^(-1074), .Machine$double.xmin, -.Machine$double.xmin)) {
    res_sub <- logit_iv_exp_scalar(sub_y, taylor_degree = 18L, max_abs_eta = 700.0, max_squarings = 11L)
    expect_identical(res_sub$m, 0L)
    expect_gt(res_sub$iv[1L], 0.0)
    # Since e^y >= 1 + y, for sub_y > 0 we have e^y > 1 and for sub_y < 0 we have e^y < 1:
    expect_true(iv_contains(res_sub$iv, 1.0))
    if (sub_y > 0) {
      expect_gt(res_sub$iv[2L], 1.0)
    } else {
      expect_lt(res_sub$iv[1L], 1.0)
    }
  }

  # 3. Exact dyadic boundaries and 1-ULP neighbours around +/- 0.5 and +/- 1.0
  half_down <- gmm_next(0.5, FALSE)
  half_up <- gmm_next(0.5, TRUE)
  r_hd <- logit_iv_exp_scalar(half_down, 18L, 700.0, 11L)
  r_h0 <- logit_iv_exp_scalar(0.5, 18L, 700.0, 11L)
  r_hu <- logit_iv_exp_scalar(half_up, 18L, 700.0, 11L)
  expect_identical(r_hd$m, 0L)
  expect_identical(r_h0$m, 0L)
  expect_identical(r_hu$m, 1L)
  expect_gt(abs(r_hu$u), 0.25)
  expect_lte(abs(r_hu$u), 0.5)

  # Interval straddling dyadic boundary [0.5, nextafter(0.5, +Inf)]
  iv_straddle <- logit_iv_exp(c(0.5, half_up), 18L, 700.0, 11L)
  expect_identical(iv_straddle$m_lo, 0L)
  expect_identical(iv_straddle$m_hi, 1L)
  expect_identical(iv_straddle$eval_count, 2L)

  # 4. High Taylor degrees N = 25, 30, 180 (where (N+1)! is not an exact binary64 integer or overflows binary64!)
  for (deg in c(2L, 18L, 25L, 30L, 180L)) {
    res_deg <- logit_iv_exp_scalar(1.25, taylor_degree = deg, max_abs_eta = 700.0, max_squarings = 11L)
    expect_gt(res_deg$iv[1L], 0.0)
    expect_lt(res_deg$iv[1L], res_deg$iv[2L])
  }

  # 5. Domain, squaring, overflow, and underflow refusals map to INCONCLUSIVE in certify_chained_logit_root
  expect_error(
    logit_iv_exp_scalar(705.0, taylor_degree = 18L, max_abs_eta = 700.0, max_squarings = 11L),
    class = "cdid_logit_arithmetic_inconclusive"
  )
  expect_error(
    logit_iv_exp_scalar(710.0, taylor_degree = 18L, max_abs_eta = 800.0, max_squarings = 11L),
    class = "cdid_logit_arithmetic_inconclusive"
  )
  expect_error(
    logit_iv_exp_scalar(-750.0, taylor_degree = 18L, max_abs_eta = 800.0, max_squarings = 11L),
    class = "cdid_logit_arithmetic_inconclusive"
  )

  # Verify that overflow/underflow in certify_chained_logit_root returns INCONCLUSIVE (not an uncaught error)
  fx <- make_exact_beta0_zero_fixture()
  ctrl_wide <- make_logit_cert_default_control(max_abs_eta = 800.0, max_squarings = 11L)
  cert_ovf <- certify_chained_logit_root(
    Z = fx$Z, w = fx$w, G = fx$G, pair_observed = fx$pair_observed,
    cohorts = fx$cohorts, beta0 = c(720.0, 0.0), F = fx$F, r = fx$r, control = ctrl_wide
  )
  expect_false(cert_ovf$certified)
  expect_identical(cert_ovf$status, "INCONCLUSIVE")
  expect_identical(cert_ovf$inconclusive_type, "arithmetic")
  expect_true(cert_ovf$arithmetic_refusal$occurred)
  expect_identical(cert_ovf$arithmetic_refusal$reason, "exp_overflow")

  cert_udf <- certify_chained_logit_root(
    Z = fx$Z, w = fx$w, G = fx$G, pair_observed = fx$pair_observed,
    cohorts = fx$cohorts, beta0 = c(-750.0, 0.0), F = fx$F, r = fx$r, control = ctrl_wide
  )
  expect_false(cert_udf$certified)
  expect_identical(cert_udf$status, "INCONCLUSIVE")
  expect_identical(cert_udf$inconclusive_type, "arithmetic")
  expect_true(cert_udf$arithmetic_refusal$occurred)
  expect_identical(cert_udf$arithmetic_refusal$reason, "exp_lower_underflow")
})

test_that("Cancellation adversaries in eta = z'beta0 and score accumulation are enclosed soundly", {
  # Construct z_1 = (2^53, 1, -2^53), beta0 = (1, 1, 1):
  # Naive float sum: (2^53 + 1) - 2^53 = 0.0, whereas exact dot product is 1.0!
  z_cancel <- c(2^53, 1.0, -2^53)
  b_cancel <- c(1.0, 1.0, 1.0)
  naive_dot <- (z_cancel[1L] * b_cancel[1L] + z_cancel[2L] * b_cancel[2L]) + z_cancel[3L] * b_cancel[3L]
  expect_identical(naive_dot, 0.0)

  iv_dot <- gmm_iv_dot(z_cancel, b_cancel)
  expect_true(iv_contains(iv_dot, 1.0))

  # Consequently, logit_iv_exp(iv_dot, ...) encloses exp(1) = e, NOT only exp(0) = 1!
  exp_cancel <- logit_iv_exp(iv_dot, taylor_degree = 18L, max_abs_eta = 700.0, max_squarings = 11L)
  expect_true(iv_contains(exp_cancel$iv, exp(1.0)))
})

test_that("Streamed chunking invariance, ALTREP inputs, and preflight work/memory bounds hold", {
  fx16 <- make_propagation_16unit_fixture()

  # 1. Streamed chunking invariance across chunk sizes 1, 2, 3, 5, 16
  certs <- lapply(c(1L, 2L, 3L, 5L, 16L), function(ch) {
    ctrl <- make_logit_cert_default_control(max_chunk_rows = ch)
    certify_chained_logit_root(
      Z = fx16$Z, w = fx16$w, G = fx16$G, pair_observed = fx16$pair_observed,
      cohorts = fx16$cohorts, beta0 = fx16$beta0, F = fx16$F, r = fx16$r, control = ctrl
    )
  })
  for (i in 2:length(certs)) {
    expect_identical(certs[[i]]$certified, certs[[1L]]$certified)
    expect_identical(certs[[i]]$delta, certs[[1L]]$delta)
    expect_identical(certs[[i]]$t, certs[[1L]]$t)
    expect_identical(certs[[i]]$a, certs[[1L]]$a)
    expect_identical(certs[[i]]$x, certs[[1L]]$x)
    expect_identical(certs[[i]]$q, certs[[1L]]$q)
    expect_identical(certs[[i]]$selfmap, certs[[1L]]$selfmap)
    expect_identical(certs[[i]]$rho, certs[[1L]]$rho)
    expect_identical(certs[[i]]$enclosures$M_cohort, certs[[1L]]$enclosures$M_cohort)
    expect_identical(certs[[i]]$enclosures$Pi_cohort, certs[[1L]]$enclosures$Pi_cohort)
    expect_identical(certs[[i]]$enclosures$g, certs[[1L]]$enclosures$g)
    expect_identical(certs[[i]]$enclosures$H, certs[[1L]]$enclosures$H)
  }
  expect_identical(certs[[1L]]$work$n_chunks_per_pass, 16L)
  expect_identical(certs[[1L]]$work$n_chunks_processed, 32L)
  expect_lte(certs[[1L]]$work$actual_interval_ops, certs[[1L]]$work$preflight_max_interval_ops)

  # 2. ALTREP integer sequence inputs (1:n) work cleanly
  n_alt <- 12L
  w_alt <- 1:n_alt # ALTREP compact intseq in R
  Z_alt <- cbind(rep(1L, n_alt), rep(c(0L, 1L), length.out = n_alt))
  G_alt <- rep(c(2L, 5L), each = 6L)
  po_alt <- rep(c(TRUE, TRUE, FALSE), length.out = n_alt)
  cert_alt <- certify_chained_logit_root(
    Z = Z_alt, w = w_alt, G = G_alt, pair_observed = po_alt,
    cohorts = c(2L, 5L), beta0 = c(0.0, 0.0), F = diag(c(0.5, 0.5)), r = 0.1,
    control = make_logit_cert_default_control(max_chunk_rows = 4L)
  )
  expect_s3_class(cert_alt, "cdid_logit_certificate")
  expect_identical(cert_alt$dimensions$n_a_units, 12L)

  # 3. Preflight matrix-entry and work-operation bounds refuse before materialization
  ctrl_small_mem <- make_logit_cert_default_control(max_matrix_entries = 10L)
  expect_error(
    certify_chained_logit_root(
      Z = fx16$Z, w = fx16$w, G = fx16$G, pair_observed = fx16$pair_observed,
      cohorts = fx16$cohorts, beta0 = fx16$beta0, F = fx16$F, r = fx16$r, control = ctrl_small_mem
    ),
    class = "cdid_logit_certificate_bound_error"
  )

  ctrl_small_ops <- make_logit_cert_default_control(max_work_ops = 50L)
  expect_error(
    certify_chained_logit_root(
      Z = fx16$Z, w = fx16$w, G = fx16$G, pair_observed = fx16$pair_observed,
      cohorts = fx16$cohorts, beta0 = fx16$beta0, F = fx16$F, r = fx16$r, control = ctrl_small_ops
    ),
    class = "cdid_logit_certificate_bound_error"
  )
})
