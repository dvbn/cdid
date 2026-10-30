test_that("private enclosure work limits precede content and optional arithmetic refuses honestly", {
  ns <- asNamespace("cdid")
  f <- dget(testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R"))$fixture2_nonzero_beta_rootball
  cc <- ns$logit_certificate_control(18,40,16,4,2500,500000)
  pc <- ns$ipw_nuisance_problem_control(4,5000,2500,500000)
  lc <- ns$ipw_link_enclosure_control(10,5000,2500,500000,5000)
  p <- ns$construct_certified_ipw_nuisance_problem(f$id,f$Z_float,f$w_float,f$G,f$pair_observed,f$cohorts,f$beta0_perturbed_float,f$F_float,f$r_ball,cc,pc)
  e <- ns$enclose_ipw_link(p,f$id,f$yL_float,f$yR_float,f$at_beta0_perturbed$lambda_beta,lc)
  expect_true(p$certified);expect_true(e$certified)
  expect_lte(p$work$total_interval_ops,p$preflight$total_work_ops_bound)
  expect_lte(e$work$actual_interval_ops,e$preflight$total_work_ops_bound)
  expect_gte(e$preflight$retained_endpoint_entries,3*nrow(p$raw$Z))
  lo <- ns$ipw_link_enclosure_control(10,5000,2500,736,5000)
  expect_error(ns$enclose_ipw_link(p,f$id,f$yL_float,f$yR_float,f$at_beta0_perturbed$lambda_beta,lo),class="cdid_ipw_enclosure_bound_error")
  lo <- ns$ipw_link_enclosure_control(10,5000,2500,1210,5000)
  expect_error(ns$enclose_ipw_link_if_chunk(e,1,10,lo,include_shared_xi=TRUE,H_chol_candidate=chol(f$at_beta0_perturbed$H)),class="cdid_ipw_enclosure_bound_error")
  tiny <- ns$ipw_link_enclosure_control(1,1,1,1,1)
  expect_error(ns$enclose_ipw_link_if_chunk(e,1,1,tiny,phi_candidate=rep(NaN,10)),class="cdid_ipw_enclosure_bound_error")
  ch <- ns$enclose_ipw_link_if_chunk(e,1,10,lc,phi_candidate=f$at_beta0_perturbed$phi,H_chol_candidate=chol(f$at_beta0_perturbed$H),include_shared_xi=TRUE)
  expect_true(ch$certified);expect_lte(ch$work$actual_interval_ops,ch$preflight$work_ops_bound)
  xmax <- .Machine$double.xmax
  for (z in list(ns$enclose_ipw_shared_xi_chunk(p,1,10,lc,xi_beta_candidate=matrix(xmax,10,2)),ns$enclose_ipw_link_if_chunk(e,1,10,lc,phi_candidate=ifelse(p$masks$in_risk,xmax,0)),ns$enclose_ipw_link(p,f$id,f$yL_float,f$yR_float,f$at_beta0_perturbed$lambda_beta,lc,candidates=list(mu=-xmax)))) {
    expect_false(z$certified);expect_identical(z$status,"INCONCLUSIVE")
    expect_false(z$separation_certified);expect_false(z$sampling_certified)
  }
  # Unqueried candidates have no role in this bounded result.
  ph <- rep(NaN,10);ph[1] <- 0
  xi <- matrix(NaN,10,2);xi[1,] <- 0
  expect_true(ns$enclose_ipw_link_if_chunk(e,1,1,lc,phi_candidate=ph)$certified)
  expect_true(ns$enclose_ipw_shared_xi_chunk(p,1,1,lc,xi_beta_candidate=xi)$certified)
  # Finite differences overflowing off O must remain unread.
  yl <- f$yL_float;yr <- f$yR_float
  yl[!p$masks$in_observed] <- -xmax;yr[!p$masks$in_observed] <- xmax
  ignored <- ns$enclose_ipw_link(p,f$id,yl,yr,f$at_beta0_perturbed$lambda_beta,lc)
  expect_true(ignored$certified);expect_identical(ignored$enclosures_star,e$enclosures_star)
})
