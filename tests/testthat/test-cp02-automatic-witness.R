# SPDX-License-Identifier: GPL-3.0-only
#
# Focused tests for the private automatic nuisance witness (one-radius rule) seam:
# `R/chained_logit_certificate.R` (typed request, selection stage),
# `R/chained_ipw_enclosure.R` (stored selected radius), `R/chained_ipw_point_view.R`
# (consumer identity) and `R/chained_ipw_pipeline.R` (diagnostic retention, charge).
#
# Reference values are EXACT RATIONALS from the accepted proof
# (docs/theory/cp02_automatic_witness.md Section 5 and its cited Theorem root),
# recomputed independently with fractions.Fraction in the retained review evidence
# and hard-coded below. Nothing is read back from the implementation as an expectation.
# Deterministic throughout: no RNG or bootstrap, no 48-unit or N200 development fit.

aw_ulp_above <- function(x, exact, ulps) x >= exact && x <= exact + ulps * .Machine$double.eps * abs(exact)
aw_ulp_below <- function(x, exact, ulps) x <= exact && x >= exact - ulps * .Machine$double.eps * abs(exact)

aw_owner_control <- function() ipw_owner_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
  max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 500000L)
aw_cert_control <- function() logit_certificate_control(taylor_degree = 12L, max_abs_eta = 50.0, max_squarings = 16L,
  max_chunk_rows = 64L, max_matrix_entries = 2000000L, max_work_ops = 500000000L)
aw_problem_control <- function() ipw_nuisance_problem_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
  max_matrix_entries = 2000000L, max_work_ops = 500000000L)

# Four intercept-only eligible units, all pairs observed, n = 4; `g_ids` are the
# target-cohort units (cohort 2, treated at date 2), the rest never treated.
aw_four_unit_view <- function(w, g_ids = 1:2, nuisance_id = "aw") {
  own <- aw_owner_control()
  d <- expand.grid(id = 1:4, t = 1:2)
  d$g <- ifelse(d$id %in% g_ids, 2, Inf); d$y <- as.double(d$id + d$t); d$w <- w[d$id]
  po <- construct_ipw_panel_owner(data = d, idname = "id", tname = "t", yname = "y", gname = "g", control = own)
  fg <- ifelse(1:4 %in% g_ids, 2, Inf)
  fo <- construct_ipw_frame_owner(panel_owner = po, id = 1:4, Z = matrix(1.0, 4, 1, dimnames = list(NULL, "(Intercept)")),
    weights = w, control = own, frame_g = fg, allow_frame_only = FALSE)
  list(nv = construct_ipw_nuisance_view(fo, 2, Inf, 1, 2, own, nuisance_id = nuisance_id), own = own, frame_owner = fo)
}

aw_pipeline_control <- function(require_matching_scale = TRUE) {
  lnk <- ipw_link_enclosure_control(max_chunk_rows = 64L, max_retained_entries = 2000000L, max_matrix_entries = 2000000L,
    max_work_ops = 500000000L, max_output_entries = 500000L)
  ipw_pipeline_control(
    owner_control = aw_owner_control(),
    catalog_bounds = list(max_candidates = 500L, max_rows = 500L, max_pool_entries = 5000L, max_scratch = 100000L, max_work = 500000000L),
    class_bounds = list(max_links = 500L, max_classes = 50L, max_class_dates = 50L, max_E_entries = 5000L, max_mask_entries = 100000L),
    fitted_control = ipw_fitted_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 500000L),
    point_control = ipw_point_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 500000L),
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    certificate_control = aw_cert_control(),
    problem_control = aw_problem_control(),
    link_control = lnk,
    covariance_control = ipw_covariance_transfer_control(max_chunk_rows = 64L, max_operator_block_rows = 64L, max_cluster_rows = 2000L,
      max_matrix_entries = 2000000L, max_work_ops = 500000000L, max_output_entries = 1000000L, link_control = lnk),
    gmm_accuracy = gmm_accuracy_control(method = "certified_eigen_chol_kkt", target_error_unit = "abs", max_target_error_abs = 1e-4,
      max_target_error_se = 1e-2, max_orthogonality_defect = 1e-8, min_certified_lambda_rel = 1e-10, max_factor_residual_rel = 1e-10,
      min_rcond_R = 1e-10, min_rcond_U = 1e-10, min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10, max_loading_residual_rel = 1e-10,
      require_matching_scale = require_matching_scale, max_loading_error_abs = 1e-4),
    gmm_bounds = gmm_bounds_control(max_matrix_entries = 2000000L, max_batch_rows = 64L, max_targets = 100L, max_open_clusters = 2000L,
      allow_materialize_if = TRUE, max_materialized_if_entries = 1000000L),
    multiplier_budget = list(max_units = 2000, max_clusters = 2000, max_targets = 100, max_draws = 2000, max_if_batch_rows = 64,
      max_if_batch_cells = 6400, max_sign_batch_rows = 32, max_sign_batch_cells = 64000, max_workspace_cells = 1000000, max_workspace_bytes = 8000000),
    inference_accuracy = list(max_estimate_error_abs = 1e-3, max_se_error_abs = 1e-2, max_vcov_error_abs = 1e-2,
      max_multiplier_error_abs = 1e-1, max_band_error_abs = 1e-1))
}

# Twelve units, three dates, two treated cohorts (2, 3) and never-treated units, one binary covariate.
aw_twelve_unit_panel <- function() {
  d <- expand.grid(id = 1:12, t = 1:3)
  d$g <- c(rep(2, 4), rep(3, 4), rep(Inf, 4))[d$id]
  d$x <- rep(c(0, 1), 6)[d$id]
  base <- c(0.3, -0.2, 0.5, 0.1, -0.4, 0.2, 0.0, 0.6, -0.1, 0.4, 0.25, -0.35)
  d$y <- base[d$id] + 0.5 * d$t + 0.2 * d$x * d$t + ifelse(d$g == 2 & d$t >= 2, 1 + 0.5 * (d$t - 2), 0) +
    ifelse(d$g == 3 & d$t >= 3, 1.5, 0) + 0.1 * sin(d$id * d$t)
  d
}
aw_auto_witness <- function(fv) list(F = unname(fv$factors$H_chol), r = logit_cert_automatic_radius_request())
aw_run <- function(ctrl, witness, ..., weighting = "optimal", d = aw_twelve_unit_panel()) run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, xformla = ~x,
  control_group = "nevertreated", weighting = weighting, root_witness = witness, cband = FALSE, ...)

# Values that must keep the ORIGINAL typed positive-radius refusal at every compact-view entry.
aw_bad_radii <- function(req) list(NA_real_, NaN, NA_integer_, NA, structure(NA_real_, foo = 1), structure(NaN, certified = TRUE),
  c(0.1, 0.2), 0, -1, Inf, -Inf, "0.5", TRUE, unclass(req), structure(list(request = "automatic_one_radius"), class = class(req)),
  { m <- req; m$request <- "automatic_two_radii"; m }, { m <- req; m$extra <- 1; m }, structure(req, class = c(class(req), "list")))

test_that("explicit typed request versus accidental NA/NaN/decorated values: only the canonical request switches mode", {
  req <- logit_cert_automatic_radius_request()
  expect_s3_class(req, "cdid_logit_automatic_radius_request")
  expect_true(logit_cert_is_automatic_radius_request(req))
  expect_true(logit_cert_is_automatic_radius_request(logit_cert_automatic_radius_request()))
  for (bad in aw_bad_radii(req)) expect_false(logit_cert_is_automatic_radius_request(bad))
  expect_false(logit_cert_is_automatic_radius_request(structure(0.5, certified = TRUE)))
  # The request is a list, never numeric: no numeric coercion path accepts it.
  expect_false(is.numeric(req)); expect_false(gmm_scalar_numeric(req, positive = TRUE))

  v <- aw_four_unit_view(c(1, 1, 1, 9 / 8))
  F1 <- matrix(0.5, 1, 1)
  for (bad in aw_bad_radii(req)) {
    expect_error(certify_chained_logit_root_view(v$nv, 0.0, F1, bad, aw_cert_control(), v$own), class = "cdid_logit_certificate_radius_error")
    expect_error(construct_certified_ipw_nuisance_problem_view(v$nv, 0.0, F1, bad, aw_cert_control(), aw_problem_control(), v$own),
      class = "cdid_ipw_enclosure_radius_error")
  }
  # Dense entry: manual positive radius only; the request itself is refused there.
  expect_error(certify_chained_logit_root(Z = matrix(1, 4, 1), w = rep(1, 4), G = c(2, 2, Inf, Inf), pair_observed = rep(TRUE, 4),
    cohorts = c(2, Inf), beta0 = 0, F = F1, r = req, control = aw_cert_control()), class = "cdid_logit_certificate_radius_error")
  expect_error(certify_chained_logit_root(Z = matrix(1, 4, 1), w = rep(1, 4), G = c(2, 2, Inf, Inf), pair_observed = rep(TRUE, 4),
    cohorts = c(2, Inf), beta0 = 0, F = F1, r = NA_real_, control = aw_cert_control()), class = "cdid_logit_certificate_radius_error")
  # Manual positive radius keeps its interface and records caller_supplied provenance.
  man <- certify_chained_logit_root_view(v$nv, 0.0, F1, 2 / 31, aw_cert_control(), v$own)
  expect_identical(man$radius_selection$mode, "caller_supplied"); expect_false(man$radius_selection$selected)
  expect_identical(man$r, 2 / 31); expect_identical(man$provenance$radius_r, 2 / 31)
  # A decorated POSITIVE numeric is (as before) the manual radius 0.5 at the view entry and never switches mode;
  # the bind entry keeps its existing decorative-attribute refusal.
  dec <- certify_chained_logit_root_view(v$nv, 0.0, F1, structure(0.5, certified = TRUE), aw_cert_control(), v$own)
  expect_identical(dec$radius_selection$mode, "caller_supplied"); expect_identical(dec$r, 0.5)
  pre_manual <- preflight_logit_certificate(4, 1, 2, aw_cert_control())
  expect_false(pre_manual$automatic_radius_reserved); expect_identical(pre_manual$automatic_radius_ops_bound, 0); expect_identical(pre_manual$automatic_radius_entries, 0)
  expect_error(preflight_logit_certificate(4, 1, 2, aw_cert_control(), automatic_radius = NA), class = "cdid_logit_certificate_radius_error")
})

test_that("selection helper on exact interval inputs: exact-zero branch, positive branch, sufficient region and typed refusals", {
  # Accepted proof fixture (H = 1/4, g = 0, F = 1/2, a = 2): r = (1 - 0)/(4 * 2) = 1/8 downward.
  z <- logit_cert_select_automatic_radius(c(0, 0), c(0, 0), c(2, 2), TRUE, TRUE)
  expect_identical(z$branch, "exact_zero_score")
  expect_true(aw_ulp_below(z$r, 1 / 8, 4))
  # Outward helpers widen the exact zero product to at most the smallest subnormal.
  expect_true(z$a_times_t_upper >= 0 && z$a_times_t_upper <= 2^-1074)
  expect_true(aw_ulp_below(z$sufficient_bound_lower, 1 / 6, 8))
  expect_true(z$inside_sufficient_region)
  # Exact rational consequences of r = 1/8 recorded in the note: x = 1/4, q = 1/3, t + q r = 1/24, rho = 0.
  expect_identical(2 * (1 / 8), 1 / 4); expect_equal((0 + 1 / 4) / (1 - 1 / 4), 1 / 3, tolerance = 0); expect_equal((1 / 3) * (1 / 8), 1 / 24, tolerance = 0)

  # Positive branch inside the sufficient region: (delta, t, a) = (0, 1/8, 1) -> r = 1/4 upward; a t = 1/8 <= 1/6.
  p <- logit_cert_select_automatic_radius(c(0, 0), c(1 / 8, 1 / 8), c(1, 1), TRUE, TRUE)
  expect_identical(p$branch, "positive_score")
  expect_true(aw_ulp_above(p$r, 1 / 4, 4))
  expect_true(aw_ulp_above(p$a_times_t_upper, 1 / 8, 4))
  expect_true(p$inside_sufficient_region)
  # Halving that radius breaks the self-map (factor 2 is load bearing): q' = 1/7, t + q' r' = 1/7 > 1/8.
  expect_gt(1 / 8 + (1 / 7) * (1 / 8), 1 / 8)
  # Outside the sufficient region: (1/4, 3/16, 1/2) -> r = 1/2, a t = 3/32 > 9/104; selection still returns one radius.
  o <- logit_cert_select_automatic_radius(c(1 / 4, 1 / 4), c(3 / 16, 3 / 16), c(1 / 2, 1 / 2), TRUE, TRUE)
  expect_true(aw_ulp_above(o$r, 1 / 2, 8)); expect_false(o$inside_sufficient_region)
  expect_true(aw_ulp_below(o$sufficient_bound_lower, 9 / 104, 16))

  # Typed refusals (arithmetic INCONCLUSIVE class), never a numeric fallback.
  e <- tryCatch(logit_cert_select_automatic_radius(c(0, 0), c(0, 0), c(0, 0), FALSE, TRUE), error = function(e) e)
  expect_s3_class(e, "cdid_logit_arithmetic_inconclusive"); expect_identical(e$reason, "auto_radius_a_not_positive")
  e <- tryCatch(logit_cert_select_automatic_radius(c(1, 1.5), c(0, 0), c(2, 2), TRUE, FALSE), error = function(e) e)
  expect_s3_class(e, "cdid_logit_arithmetic_inconclusive"); expect_identical(e$reason, "auto_radius_delta_ge_1")
  e <- tryCatch(logit_cert_select_automatic_radius(c(0, 0), c(NA_real_, Inf), c(2, 2), TRUE, TRUE), error = function(e) e)
  expect_s3_class(e, "cdid_logit_arithmetic_inconclusive"); expect_identical(e$reason, "auto_radius_t_invalid")
  e <- tryCatch(logit_cert_select_automatic_radius(c(0, 0), c(0, 0), c(0, 1e-320), TRUE, TRUE), error = function(e) e)
  expect_s3_class(e, "cdid_logit_arithmetic_inconclusive"); expect_true(e$reason %in% c("auto_radius_4a_nonpositive", "auto_radius_nonpositive"))
})

test_that("compact-view POS fixture: one outward radius encloses the exact rationals, original inequalities certify, manual-at-selected is identical", {
  req <- logit_cert_automatic_radius_request()
  v <- aw_four_unit_view(c(1, 1, 1, 9 / 8))
  cert <- certify_chained_logit_root_view(v$nv, 0.0, matrix(0.5, 1, 1), req, aw_cert_control(), v$own)
  expect_s3_class(cert, "cdid_logit_certificate")
  expect_identical(cert$status, "CERTIFIED"); expect_true(cert$certified); expect_true(all(unlist(cert$inequalities)))
  expect_identical(cert$radius_selection$mode, "automatic_one_radius"); expect_true(cert$radius_selection$selected)
  expect_identical(cert$radius_selection$branch, "positive_score"); expect_true(cert$radius_selection$inside_sufficient_region)
  # Exact rationals: delta = 1/32, t = 1/32, a = 2, r = 2/31, x = 4/31, q = 53/288, t + q r = 385/8928, rho = 9/235.
  tol <- 1e-11  # data-pass widening of delta, t (about 1e-15 absolute) propagates outward; never below the exact value
  expect_true(cert$delta >= 1 / 32 && cert$delta <= (1 / 32) * (1 + tol))
  expect_true(cert$t >= 1 / 32 && cert$t <= (1 / 32) * (1 + tol))
  expect_true(cert$a >= 2 && cert$a <= 2 * (1 + tol))
  expect_true(cert$r >= 2 / 31 && cert$r <= (2 / 31) * (1 + tol))
  expect_true(cert$x >= 4 / 31 && cert$x <= (4 / 31) * (1 + tol))
  expect_true(cert$q >= 53 / 288 && cert$q <= (53 / 288) * (1 + tol))
  expect_true(cert$selfmap >= 385 / 8928 && cert$selfmap <= (385 / 8928) * (1 + tol))
  expect_true(cert$rho >= 9 / 235 && cert$rho <= (9 / 235) * (1 + tol))
  expect_true(cert$radius_selection$a_times_t_upper >= 1 / 16 && cert$radius_selection$sufficient_bound_lower <= 961 / 6208)
  # The one radius formula on the certificate's own upper endpoints (few ulps above the one-rounding ideal).
  ideal <- 2 * cert$t / (1 - cert$delta)
  expect_true(cert$r >= ideal * (1 - 2^-50) && cert$r <= ideal * (1 + 2^-48))
  # Stored consistently.
  expect_identical(cert$r, cert$provenance$radius_r); expect_identical(cert$r, cert$radius_selection$r)
  # Resource charge: reserved 64 ops / 16 entries only in automatic mode; debited exactly 64 relative to the manual run at the selected radius.
  expect_true(cert$work$preflight_max_interval_ops >= cert$work$actual_interval_ops)
  man <- certify_chained_logit_root_view(v$nv, 0.0, matrix(0.5, 1, 1), cert$r, aw_cert_control(), v$own)
  expect_identical(cert$work$actual_interval_ops - man$work$actual_interval_ops, 64)
  expect_identical(cert$work$preflight_max_interval_ops - man$work$preflight_max_interval_ops, 64)
  expect_identical(cert$work$preflight_concurrent_entries - man$work$preflight_concurrent_entries, 16)
  same <- setdiff(names(cert), c("work", "radius_selection"))
  expect_identical(cert[same], man[same])
  pre <- preflight_logit_certificate_view(v$nv, aw_cert_control(), v$own, automatic_radius = TRUE)
  expect_true(pre$automatic_radius_reserved); expect_identical(pre$automatic_radius_ops_bound, 64); expect_identical(pre$automatic_radius_entries, 16)

  # ZERO fixture through the data path: the raw passes widen delta and t to ~1e-15, so the positive branch runs with a tiny radius;
  # the exact-zero branch is reached only on exactly zero certified t (previous test). Certification still holds.
  z <- certify_chained_logit_root_view(aw_four_unit_view(c(1, 1, 1, 1))$nv, 0.0, matrix(0.5, 1, 1), req, aw_cert_control(), v$own)
  expect_identical(z$status, "CERTIFIED"); expect_identical(z$radius_selection$branch, "positive_score")
  expect_true(z$t > 0 && z$t < 1e-13 && z$r < 1e-13 && z$rho <= z$r)
})

test_that("REFUSE fixture: selection returns one radius, the original inequalities refuse, result is INCONCLUSIVE with the stored radius (never nonexistence)", {
  req <- logit_cert_automatic_radius_request()
  v <- aw_four_unit_view(c(1, 1, 1, 1), g_ids = 1:3)
  cert <- certify_chained_logit_root_view(v$nv, 0.0, matrix(0.5, 1, 1), req, aw_cert_control(), v$own)
  # Exact rationals: delta = 0, t = 1/2, a = 2 -> r = 1, x = 2 >= 1; a t = 1 > 1/6.
  tol <- 1e-11
  expect_identical(cert$status, "INCONCLUSIVE"); expect_false(cert$certified); expect_identical(cert$inconclusive_type, "inequality")
  expect_true("x_ge_1" %in% cert$inconclusive_reasons); expect_false(cert$inequalities$x_lt_1)
  expect_true(cert$radius_selection$selected); expect_identical(cert$radius_selection$branch, "positive_score")
  expect_true(cert$r >= 1 && cert$r <= 1 + tol); expect_true(cert$x >= 2 && cert$x <= 2 * (1 + tol))
  expect_false(cert$radius_selection$inside_sufficient_region); expect_true(cert$radius_selection$a_times_t_upper >= 1)
  expect_identical(cert$r, cert$provenance$radius_r)
  expect_false(cert$separation_certified); expect_false(cert$sampling_certified)
  # Through the enclosure: uncertified problem, stored radius identical to the certificate's selected scalar.
  prob <- construct_certified_ipw_nuisance_problem_view(v$nv, 0.0, matrix(0.5, 1, 1), req, aw_cert_control(), aw_problem_control(), v$own)
  expect_false(prob$certified); expect_identical(prob$status, "INCONCLUSIVE")
  expect_identical(prob$raw$r, prob$certificate$r); expect_identical(prob$raw$r, prob$certificate$provenance$radius_r)
})

test_that("selection failure: arithmetic INCONCLUSIVE with r = NA_real_ inside the certificate only; pipeline refuses before any consumer", {
  req <- logit_cert_automatic_radius_request()
  v <- aw_four_unit_view(c(1, 1, 1, 9 / 8))
  F_bad <- matrix(0.5 * 2^-20, 1, 1)  # M = 2^-42 << H: delta >= 1, automatic selection refuses
  cert <- certify_chained_logit_root_view(v$nv, 0.0, F_bad, req, aw_cert_control(), v$own)
  expect_identical(cert$status, "INCONCLUSIVE"); expect_identical(cert$inconclusive_type, "arithmetic")
  expect_identical(cert$inconclusive_reasons, "auto_radius_delta_ge_1")
  expect_true(is.na(cert$r)); expect_true(is.na(cert$provenance$radius_r)); expect_true(is.na(cert$radius_selection$r))
  expect_false(cert$radius_selection$selected); expect_identical(cert$radius_selection$refusal_reason, "auto_radius_delta_ge_1")
  expect_identical(cert$radius_selection$mode, "automatic_one_radius")
  expect_true(cert$raw_score_hessian_enclosures_complete == FALSE)
  prob <- construct_certified_ipw_nuisance_problem_view(v$nv, 0.0, F_bad, req, aw_cert_control(), aw_problem_control(), v$own)
  expect_false(prob$certified); expect_true(is.na(prob$raw$r)); expect_identical(prob$raw$r, prob$certificate$r)
  # Manual radius with the same bad metric: the pre-existing inequality refusal (delta_ge_1), unchanged path.
  man <- certify_chained_logit_root_view(v$nv, 0.0, F_bad, 1e-3, aw_cert_control(), v$own)
  expect_identical(man$status, "INCONCLUSIVE"); expect_true("delta_ge_1" %in% man$inconclusive_reasons); expect_identical(man$r, 1e-3)

  # Pipeline: typed transfer refusal naming the nuisance and the selection reason; no GMM/link consumer ran on a radius.
  e <- tryCatch(aw_run(aw_pipeline_control(), function(fv) list(F = unname(fv$factors$H_chol) * 2^-20, r = req),
    base_period = "varying", aggregations = "simple"), error = function(e) e)
  expect_s3_class(e, "cdid_ipw_pipeline_transfer_inconclusive"); expect_s3_class(e, "cdid_logit_certificate_inconclusive")
  expect_identical(e$diagnostics$reasons, "auto_radius_delta_ge_1"); expect_match(e$diagnostics$nuisance_id, "^nuisance:")
  # Accidental NA in a witness is the typed radius refusal, never automatic mode.
  e2 <- tryCatch(aw_run(aw_pipeline_control(), function(fv) list(F = unname(fv$factors$H_chol), r = NA_real_),
    base_period = "varying", aggregations = "simple"), error = function(e) e)
  expect_s3_class(e2, "cdid_ipw_enclosure_radius_error")
})

test_that("pipeline automatic witness: selected-radius identity in problem/provenance/consumers, unchanged science, retained diagnostics and charge", {
  req <- logit_cert_automatic_radius_request()
  ctrl <- aw_pipeline_control()
  res <- aw_run(ctrl, aw_auto_witness, base_period = "varying", aggregations = c("simple", "group", "dynamic", "calendar"))
  gt <- res$numerical_transfer$gmm_transfer
  expect_true(gt$complete_inference_certified); expect_true(gt$supported_subvector_certified)
  J <- length(gt$total_target_error_cert); expect_gt(J, 0L)
  # Every nuisance problem stores the one selected scalar consistently.
  tab <- res$numerical_transfer$nuisance_root_radii
  expect_s3_class(tab, "data.frame"); expect_identical(nrow(tab), length(res$numerical_transfer$nuisance_problems))
  expect_true(all(tab$mode == "automatic_one_radius")); expect_true(all(tab$branch == "positive_score"))
  expect_true(all(is.finite(tab$r) & tab$r > 0 & tab$rho <= tab$r)); expect_true(all(tab$inside_sufficient_region))
  for (pr in res$numerical_transfer$nuisance_problems) {
    expect_true(pr$certified)
    expect_identical(pr$raw$r, pr$certificate$r); expect_identical(pr$raw$r, pr$certificate$provenance$radius_r)
    expect_identical(pr$raw$r, pr$certificate$radius_selection$r)
    expect_identical(pr$raw$F, pr$certificate$provenance$F_factor)
  }
  # Same science as a manual witness: the radius changes certification only, not estimates, IFs or SEs.
  man <- aw_run(ctrl, function(fv) list(F = unname(fv$factors$H_chol), r = 1e-6), base_period = "varying",
    aggregations = c("simple", "group", "dynamic", "calendar"))
  expect_identical(res$reporting$estimates, man$reporting$estimates); expect_identical(res$reporting$if_full, man$reporting$if_full)
  expect_identical(res$reporting$analytical_se, man$reporting$analytical_se)
  expect_true(all(man$numerical_transfer$nuisance_root_radii$mode == "caller_supplied")); expect_true(all(man$numerical_transfer$nuisance_root_radii$r == 1e-6))
  # Retained components: each outward term is bounded by the unchanged certified total; nothing recomputed on empty arrays.
  pc <- gt$point_error_components; rc <- gt$R_c_components
  for (nm in c("solve", "moment_change", "changing_weights")) { expect_length(pc[[nm]], J); expect_true(all(is.finite(pc[[nm]]) & pc[[nm]] >= 0 & pc[[nm]] <= gt$total_target_error_cert)) }
  for (nm in c("projection", "row_if_change", "numerical_loading", "changing_weights")) { expect_length(rc[[nm]], J); expect_true(all(is.finite(rc[[nm]]) & rc[[nm]] >= 0 & rc[[nm]] <= gt$R_c_cert)) }
  expect_identical(pc$solve, gt$delta_theta_solve); expect_identical(rc$projection, gt$projection_error_cert)
  expect_true(all(pc$changing_weights > 0)); expect_true(all(rc$changing_weights > 0))
  expect_true(gt$certified_lower_eigenvalue$applies); expect_identical(gt$certified_lower_eigenvalue$value, gt$l_lower)
  expect_identical(gt$certified_lower_eigenvalue$value, as.double(res$gmm_fit$propagation_bounds$l_lambda_min_V_cert))
  # Certified bounds depend on the radius through rho (the manual 1e-6 radius is far larger than the selected ~1e-14);
  # both are finite positive outward bounds and the estimates above are identical.
  expect_true(all(is.finite(gt$total_target_error_cert) & gt$total_target_error_cert > 0))
  expect_true(all(is.finite(man$numerical_transfer$gmm_transfer$total_target_error_cert)))
  expect_true(all(tab$r < 1e-6))
  # Aggregation boxes enclose the stored weights/gradients; degree-0 (group) rows have exact-zero gradient boxes.
  for (nm in names(res$aggregations)) {
    z <- res$aggregations[[nm]]
    expect_identical(dim(z$weight_enclosure$lo), dim(z$weights_theta)); expect_identical(dim(z$share_gradient_enclosure$lo), dim(z$gradient_shares))
    expect_true(all(z$weight_enclosure$lo <= z$weights_theta & z$weights_theta <= z$weight_enclosure$hi))
    expect_true(all(z$share_gradient_enclosure$lo <= z$gradient_shares & z$gradient_shares <= z$share_gradient_enclosure$hi))
    av <- z$cell_status$available
    expect_true(all(z$cell_estimate_enclosure$lo[av] <= z$estimates[av] & z$estimates[av] <= z$cell_estimate_enclosure$hi[av]))
    expect_identical(names(z$cell_estimate_enclosure$lo), z$cell_labels)
  }
  expect_true(all(res$aggregations$group$share_gradient_enclosure$lo == 0 & res$aggregations$group$share_gradient_enclosure$hi == 0))
  expect_true(any(res$aggregations$calendar$share_gradient_enclosure$hi != 0))
  # Charge: the retained diagnostic entries are part of the existing simultaneous envelope.
  ra <- res$resource_accounting
  expect_true(ra$d5_arrays$retained_diagnostic_entries > 0)
  expect_true(ra$simultaneous_entries >= ra$initial_envelope_entries + ra$d5_arrays$retained_diagnostic_entries)
  expect_true(ra$reachable_numeric_payload_only); expect_false(ra$peak_heap_or_RSS_claim)

  # Consumer identity (`certify_ipw_point_link_view` with a supplied problem). Its own binding preflight reserves the
  # problem and link caps on top of the point cap, so a larger point cap is needed to reach the identity check.
  pl <- res$point_links_used[[1L]]
  fv <- get("fitted_nuisance_view", pl$fitted_nuisance_ref); lv <- get("link_view", pl$link_view_ref)
  prob <- get("nuisance_problem", res$numerical_transfer$link_enclosures[[1L]]$problem_ref)
  big_pt <- ipw_point_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 40000000,
    max_matrix_entries = 40000000, max_work_ops = 4000000000, max_output_entries = 40000000)
  consume <- function(r, problem = prob) certify_ipw_point_link_view(fv, lv, F = unname(fv$factors$H_chol), r = r, point_control = big_pt,
    certificate_control = ctrl$certificate_control, problem_control = ctrl$problem_control, link_control = ctrl$link_control,
    reader_control = ctrl$owner_control, problem = problem)
  expect_true(consume(req)$certified)
  expect_true(consume(prob$raw$r)$certified)
  expect_error(consume(prob$raw$r * 2), class = "cdid_ipw_point_view_candidate_mismatch_error")
  expect_error(consume(NA_real_), class = "cdid_ipw_point_view_candidate_mismatch_error")
  prob_manual <- bind_fitted_nuisance_to_problem_view(fv, F = unname(fv$factors$H_chol), r = prob$raw$r,
    certificate_control = ctrl$certificate_control, problem_control = ctrl$problem_control, reader_control = ctrl$owner_control)
  expect_identical(prob_manual$certificate$radius_selection$mode, "caller_supplied"); expect_identical(prob_manual$raw$r, prob$raw$r)
  expect_error(consume(req, prob_manual), class = "cdid_ipw_point_view_candidate_mismatch_error")
  # Bind with the request: the stored radius is the certificate's selected scalar, not the request.
  prob_bind <- bind_fitted_nuisance_to_problem_view(fv, F = unname(fv$factors$H_chol), r = req,
    certificate_control = ctrl$certificate_control, problem_control = ctrl$problem_control, reader_control = ctrl$owner_control)
  expect_true(is.numeric(prob_bind$raw$r)); expect_identical(prob_bind$raw$r, prob$raw$r)
  expect_identical(prob_bind$certificate$radius_selection$mode, "automatic_one_radius")
})

test_that("identity weighting and reference-only cells: exact-zero changing-weight terms, no eigenvalue claim, exact-zero reference boxes", {
  ri <- aw_run(aw_pipeline_control(require_matching_scale = FALSE), aw_auto_witness, base_period = "varying",
    aggregations = c("simple", "calendar"), weighting = "identity")
  gi <- ri$numerical_transfer$gmm_transfer
  expect_identical(ri$numerical_transfer$weighting, "identity"); expect_true(gi$supported_subvector_certified)
  J <- length(gi$total_target_error_cert); expect_gt(J, 0L)
  expect_identical(gi$point_error_components$changing_weights, rep(0.0, J)); expect_identical(gi$R_c_components$changing_weights, rep(0.0, J))
  expect_true(all(gi$point_error_components$moment_change <= gi$total_target_error_cert))
  expect_true(all(gi$point_error_components$solve <= gi$total_target_error_cert))
  expect_false(gi$certified_lower_eigenvalue$applies); expect_true(is.na(gi$certified_lower_eigenvalue$value))
  expect_true(all(ri$numerical_transfer$nuisance_root_radii$mode == "automatic_one_radius"))

  ru <- aw_run(aw_pipeline_control(), aw_auto_witness, base_period = "universal", aggregations = c("simple", "calendar"))
  tm <- ru$reporting$target_meta; ref <- which(tm$role == "reference")
  expect_gt(length(ref), 0L)
  expect_identical(unname(ru$reporting$estimates[ref]), rep(0.0, length(ref))); expect_identical(unname(ru$reporting$se[ref]), rep(0.0, length(ref)))
  expect_true(all(ru$reporting$numerical_certificate$se[, ref] == 0))
  expect_true(all(ru$reporting$numerical_certificate$theta$lo[, ref] == 0 & ru$reporting$numerical_certificate$theta$hi[, ref] == 0))
  expect_true(ru$numerical_transfer$gmm_transfer$complete_inference_certified)
})


test_that("automatic requests conflict with a stored manual candidate radius", {
  v <- aw_four_unit_view(c(1, 1, 1, 9 / 8))
  nv <- construct_ipw_nuisance_view(v$frame_owner, 2, Inf, 1, 2, v$own,
    beta0 = 0.0, F = matrix(0.5, 1, 1), r = 2 / 31, nuisance_id = "aw:radius-conflict")
  req <- logit_cert_automatic_radius_request()
  expect_error(certify_chained_logit_root_view(nv, 0.0, matrix(0.5, 1, 1), req, aw_cert_control(), v$own),
    class = "cdid_logit_certificate_candidate_mismatch_error")
  expect_error(construct_certified_ipw_nuisance_problem_view(nv, 0.0, matrix(0.5, 1, 1), req,
    aw_cert_control(), aw_problem_control(), v$own), class = "cdid_ipw_enclosure_candidate_mismatch_error")
})

# Numeric/logical payload only, without traversing shared owner environments.
# Sum repeated appearances conservatively; no allocator/RSS claim follows.
aw_diagnostic_payload <- function(x) {
  if (is.environment(x) || is.function(x)) return(0.0)
  payload <- if (is.list(x)) sum(vapply(x, aw_diagnostic_payload, numeric(1L))) else
    if (is.numeric(x) || is.logical(x)) as.double(length(x)) else 0.0
  # Numeric dimension/row-name metadata are payload; SEXP header bytes are excluded.
  payload + sum(vapply(attributes(x), aw_diagnostic_payload, numeric(1L)))
}

aw_capture_simultaneous_guards <- function(run) {
  env <- environment(run_chained_ipw_pipeline)
  original <- get("ipw_owner_safe_sum", env, inherits = FALSE)
  was_locked <- bindingIsLocked("ipw_owner_safe_sum", env)
  replace <- function(value) {
    if (bindingIsLocked("ipw_owner_safe_sum", env)) unlockBinding("ipw_owner_safe_sum", env)
    assign("ipw_owner_safe_sum", value, env)
    if (was_locked) lockBinding("ipw_owner_safe_sum", env)
  }
  records <- list()
  labels <- c("pipeline simultaneous point producer, interval reader and GMM consumer",
    "simultaneous auxiliary producer, retained post outputs and interval reader",
    "simultaneous auxiliary point producer, compact interval reader, nuisance cache and retained post/report outputs")
  spy <- function(terms, max_bound = 2^53 - 1, label = "sum") {
    if (label %in% labels) records[[length(records) + 1L]] <<- list(label = label, persistent = terms[1L])
    original(terms, max_bound, label)
  }
  replace(spy)
  on.exit(replace(original), add = TRUE)
  value <- run()
  list(value = value, records = records)
}

test_that("retained diagnostics remain reserved in covariance and auxiliary stages and dominate the stored payload", {
  z <- aw_capture_simultaneous_guards(function() aw_run(aw_pipeline_control(), aw_auto_witness,
    base_period = "varying", aggregations = c("simple", "group", "dynamic", "calendar")))
  r <- z$value; nt <- r$numerical_transfer; gt <- nt$gmm_transfer
  charge <- r$resource_accounting$d5_arrays$retained_diagnostic_entries
  expect_gt(length(z$records), 1L)
  for (guard in z$records) expect_gte(guard$persistent, r$resource_accounting$initial_envelope_entries + charge)
  problem_payload <- function(p) aw_diagnostic_payload(p$certificate$radius_selection) +
    aw_diagnostic_payload(p$preflight[c("automatic_radius_reserved", "automatic_radius_ops_bound", "automatic_radius_entries")]) +
    aw_diagnostic_payload(p$preflight$raw_certificate_preflight[c("automatic_radius_reserved", "automatic_radius_ops_bound", "automatic_radius_entries")])
  payload <- aw_diagnostic_payload(gt[c("point_error_components", "R_c_components", "certified_lower_eigenvalue")]) +
    aw_diagnostic_payload(nt$nuisance_root_radii) +
    sum(vapply(nt$nuisance_problems, problem_payload, numeric(1L))) +
    sum(vapply(r$auxiliary_reporting$problems, function(p) problem_payload(p$problem), numeric(1L))) +
    sum(vapply(r$aggregations, function(a) aw_diagnostic_payload(a[c("weight_enclosure", "share_gradient_enclosure", "cell_estimate_enclosure")]), numeric(1L)))
  payload <- payload + 1.0  # stored retained_diagnostic_entries accounting scalar
  expect_lte(payload, charge)
})

test_that("automatic diagnostics preserve zero loading and empty-catalog support branches", {
  d <- expand.grid(id = 1:8, t = 1:3); d$g <- ifelse(d$id <= 4, 3, Inf)
  d$y <- ifelse(d$t == 1, 0, NA_real_)
  d$y[d$t == 3] <- c(1, 1, 3, 3, -1, -1, 1, 1)[d$id[d$t == 3]]
  run <- function(data) run_chained_ipw_pipeline(data, "id", "t", "y", "g", aw_pipeline_control(),
    target_cohorts = 3, base_period = "universal", root_witness = aw_auto_witness, cband = FALSE,
    aggregations = c("simple", "group", "dynamic", "calendar"))
  r <- run(d); gt <- r$numerical_transfer$gmm_transfer
  expect_identical(ncol(r$gmm_fit$if_projection$materialized_target_if), 0L)
  for (x in gt$point_error_components[1:3]) expect_length(x, 0L)
  for (x in gt$R_c_components[1:4]) expect_length(x, 0L)
  expect_false(gt$complete_inference_certified)
  expect_true(all(r$numerical_transfer$nuisance_root_radii$mode == "automatic_one_radius"))
  ref <- which(r$reporting$target_meta$role == "reference")
  expect_identical(unname(r$reporting$se[ref]), rep(0.0, length(ref)))

  # Target units are observed once; controls still form pairs. No target comparison is supported.
  d$y[d$t == 3] <- NA_real_
  d$y[d$t == 2 & d$g == Inf] <- c(-1, -1, 1, 1)
  e <- run(d)
  expect_identical(e$support_status, "empty_catalog"); expect_null(e$gmm_fit)
  expect_length(e$auxiliary_reporting$problems, 0L)
  expect_false(e$reporting$numerical_certificate$whole_request_complete)
  # This universal-base request reports the two dates before g = 3 plus its
  # one post ledger cell. The union of both auxiliary families is work, not
  # six retained reporting columns. Four aggregation types retain lower/upper
  # boxes over J cells, two cohort shares and one point coordinate, plus eight
  # scalar slots per type and the accounting scalar.
  expect_identical(e$resource_accounting$d5_arrays$n_ledger_cells, 1L)
  J <- 1.0 + 2.0
  agg_charge <- 2.0 * 4.0 * J * (J + 2.0 + 1.0)
  expect_equal(e$resource_accounting$d5_arrays$retained_diagnostic_entries,
    agg_charge + 8.0 * 4.0 + 1.0, tolerance = 0)
})
