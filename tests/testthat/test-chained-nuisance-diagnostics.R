# Deterministic checks: source-tree tests use the original direct-source path;
# installed tests resolve the actual locked package namespace. No fit or RNG.
nd_source_files <- c("chained_ipw_owner.R", "chained_ipw_point_view.R", "chained_ipw_pipeline.R",
  "chained_ipw_fitted_view.R", "chained_did_request.R")
nd_r_dir <- testthat::test_path("..", "..", "R")
if (all(file.exists(file.path(nd_r_dir, nd_source_files)))) {
  nd_source <- new.env(parent = baseenv())
  for (nd_file in nd_source_files) sys.source(file.path(nd_r_dir, nd_file), envir = nd_source)
} else {
  if (!requireNamespace("cdid", quietly = TRUE))
    stop("Actual installed cdid required when package R sources are unavailable.", call. = FALSE)
  nd_source <- asNamespace("cdid")
  stopifnot(isNamespace(nd_source), environmentIsLocked(nd_source))
  for (nd_name in c("ipw_owner_seal", "ipw_fitted_diagnostic_summary", "chained_request_H_diagnostics",
    "chained_request_nuisance_diagnostics", "fit_ipw_nuisance_view",
    "preflight_ipw_fitted_nuisance_view", "chained_request_bind_mp")) {
    nd_fun <- get(nd_name, envir = nd_source, inherits = FALSE)
    stopifnot(is.function(nd_fun), identical(environment(nd_fun), nd_source),
      bindingIsLocked(nd_name, nd_source))
  }
}

nd_fixture <- function() {
  summary <- nd_source$ipw_fitted_diagnostic_summary
  make_fit <- function(id) nd_source$ipw_owner_seal(list(
    nuisance_id = id, defining_pair = list(s = 1, t = 2),
    cohorts = list(g = 3, K = c(5, Inf), cohorts = c(3, 5, Inf), cohort_labels = c("3", "5", "Inf")),
    dimensions = list(n_a_units = 12L, n_observed_risk_units = 9L,
      n_covariates = 2L, covariate_names = c("(Intercept)", "x")),
    counts = list(n_a_units = 12L, n_observed_risk_units = 9L, n_treated_obs = 3L,
      n_control_obs = 6L, n_nonobserved_risk_units = 3L,
      by_cohort_a = c(4L, 4L, 4L), by_cohort_pair = c(3L, 3L, 3L)),
    masses = list(cohort_a_weight = c(4, 4, 4), cohort_pair_weight = c(3, 3, 3),
      q = c(1/3, 1/3, 1/3), rates = c(.75, .75, .75), view_chunked_rates = c(.75, .75, .75)),
    means = list(sum_a = 4, sum_b = 4), coefficients = c(0, 1),
    factors = list(H = matrix(c(4, 1, 1, 2), 2), H_chol = diag(2), qr_R = diag(2)),
    diagnostics = list(converged = TRUE, boundary = FALSE, rank = 2L, iter = 5L,
      weight_concentration = list(propensity_fit_weights = summary(c(1, 2, 3), 2, TRUE),
        target_outcome_weights = summary(c(1, 3), 2, TRUE),
        pooled_control_outcome_ipw_weights = summary(c(1, 2, 4), 2, TRUE)),
      fitted_probability_range = summary(c(.25, .5, .75), 2, FALSE),
      fitted_probability_scope = "observed fitted risk set"),
    glm_control = list(epsilon = 1e-8, maxit = 25L, trace = FALSE), control = list(max_chunk_rows = 2L),
    point_summary_disclosure = list(defining_rate_source = "observer_buffer_full_vector_sum_over_frame_owner_cohort_a_mass")),
    "fitted_nuisance_view")
  f1 <- make_fit("same-label"); f2 <- make_fit("same-label"); f3 <- make_fit("aux-only")
  pl <- function(id, f) list(link_id = id, fitted_nuisance_ref = f$ref)
  aux <- list(fitted_views = list(f1, f2, f1, f3),
    point_links = list(pl("reporting:1", f1), pl("reporting:2", f3), NULL),
    descriptors = list(list(g = 3, K = c(5, Inf), s = 1, t = 2),
      list(g = 3, K = c(5, Inf), s = 1, t = 3), list(g = 5, K = Inf, s = 1, t = 2)),
    selections = list(varying = list(list(index = 1L, sign = 1, status = "supported")),
      universal = list(list(index = 2L, sign = -1, status = "supported"),
        list(index = 0L, sign = 0, status = "reference_zero"),
        list(index = 3L, status = "unsupported_pair"))),
    descriptor_status = data.frame(index = 1:3, status = c("supported", "supported", "unsupported_pair"),
      reason = c(NA_character_, NA_character_, "no target observers")),
    unsupported_counts = list(NULL, NULL, list(n_a_units = 12L, n_observed_risk_units = 0L)))
  res <- list(auxiliary_reporting = aux, point_links_all = list(pl("link:r1", f1),
      pl("link:r2", f1), pl("link:r3", f2)),
    nuisance_registry = list(fitted_views = list(f1, f2)),
    frame_owner = list(n_total_units = 20L, n_frame_only_units = 4L, n_a_units = 12L), support_status = "partial_support",
    covariance_handling = list(strategy = "unchanged sentinel"),
    resource_accounting = list(initial_envelope_entries = 1000, simultaneous_entries = 2000, downstream_work_ops = 500))
  list(res = res, control = list(owner_control = list(max_retained_entries = 1e8,
    max_output_entries = 1e8, max_matrix_entries = 1e8, max_work_ops = 1e10)))
}

test_that("Kish uses exact positive scale normalization and distinct fixed vectors", {
  f <- nd_source$ipw_fitted_diagnostic_summary
  # Independently: (1+2+3)^2/(1+4+9)=18/7; (1+3)^2/(1+9)=8/5.
  for (chunk in c(1, 2, 100)) for (scale in c(1, 2^900, 2^-900)) {
    expect_equal(f(c(1, 2, 3) * scale, chunk, TRUE)$value, 18/7)
    expect_equal(f(c(1, 3) * scale, chunk, TRUE)$value, 8/5)
    expect_equal(f(c(1, 2, 4) * scale, chunk, TRUE)$value, 7/3)
    expect_equal(f(rep(scale, 4), chunk, TRUE)$value, 4)
  }
  for (v in list(numeric(), NULL, c(1, Inf), c(1, NA_real_), c(1, 0), c(-1, 1))) {
    z <- f(v, 2, TRUE)
    expect_identical(z$status, "unavailable")
    expect_true(is.character(z$reason) && !is.na(z$reason))
    expect_identical(z$value, NA_real_)
  }
  expect_identical(f(numeric(), 2, TRUE)$reason, "empty_vector")
  expect_identical(f(c(0, 1), 2, TRUE)$reason, "nonpositive_weights")
  expect_identical(f(c(.125, .9, .25), 2, FALSE),
    list(status = "available", reason = NA_character_, min = .125, max = .9))
  expect_identical(f(c(0, 1), 1, FALSE)$min, 0)
  expect_identical(f(numeric(), 2, FALSE)$status, "unavailable")
  expect_identical(f(c(.5, NaN), 2, FALSE)$reason, "nonfinite_vector")
})

test_that("H diagnostics use the stored matrix, 1-norm label and A normalization", {
  f <- nd_source$chained_request_H_diagnostics
  z <- f(matrix(c(4, 1, 1, 2), 2), 12L)
  # ||H||_1=5; H^-1=(1/7)*[2,-1;-1,4], ||H^-1||_1=5/7.
  expect_equal(z$rcond_H_1_estimate$value, 7/25)
  expect_identical(z$method, 'base::rcond(H, norm = "O")')
  expect_identical(z$norm, "O")
  expect_identical(z$norm_label, "1-norm")
  expect_identical(z$n_A, 12L)
  expect_identical(z$diagonal_range, list(status = "available", reason = NA_character_, min = 2, max = 4))
  expect_equal(f(diag(c(2, 8)), 100)$rcond_H_1_estimate$value, 1/4)
  singular <- f(diag(c(1, 0)), 12L)
  expect_identical(singular$rcond_H_1_estimate$status, "available")
  expect_identical(singular$rcond_H_1_estimate$value, 0)
  for (H in list(matrix(NA_real_, 2, 2), matrix(numeric(), 0, 0))) {
    bad <- f(H, 12L)
    expect_identical(bad$rcond_H_1_estimate$status, "unavailable")
    expect_identical(bad$rcond_H_1_estimate$value, NA_real_)
    expect_identical(bad$diagonal_range$status, "unavailable")
  }
})

test_that("true fit identity binds post, auxiliary reuse, descriptors and no references", {
  x <- nd_fixture(); before <- serialize(x$res, NULL)
  z <- nd_source$chained_request_nuisance_diagnostics(x$res, x$control, FALSE)
  expect_identical(serialize(x$res, NULL), before)
  expect_length(z$records, 3L) # repeated ref merges; equal label/g/K/beta does not
  expect_identical(vapply(z$records, `[[`, character(1), "id"), paste0("nuisance:", 1:3))
  expect_identical(vapply(z$post_links, `[[`, character(1), "nuisance_id"), c("nuisance:1", "nuisance:1", "nuisance:2"))
  expect_identical(vapply(z$auxiliary_links, `[[`, character(1), "nuisance_id"), c("nuisance:1", "nuisance:3", NA_character_))
  expect_identical(z$records[[1]]$roles, list(post = TRUE, auxiliary = TRUE))
  expect_identical(z$records[[2]]$roles, list(post = TRUE, auxiliary = FALSE))
  expect_identical(z$records[[3]]$roles, list(post = FALSE, auxiliary = TRUE))
  expect_identical(z$auxiliary_links[[3]]$reason, "no_supported_pair_no_fitted_nuisance")
  for (nm in c("descriptors", "selections", "descriptor_status", "unsupported_counts"))
    expect_identical(z[[paste0("auxiliary_", nm)]], x$res$auxiliary_reporting[[nm]])
  fv <- x$res$auxiliary_reporting$fitted_views[[1]]
  for (nm in c("counts", "masses", "dimensions", "glm_control", "control", "point_summary_disclosure"))
    expect_identical(z$records[[1]][[nm]], fv[[nm]])
  expect_identical(z$records[[1]]$fit, fv$diagnostics)
  expect_identical(z$frame_counts, list(n_full_frame_units = 24L, n_pair_eligible_A_units = 12L))
  expect_identical(unserialize(serialize(z, NULL)), z)
  violations <- character()
  visit <- function(v) {
    if (is.environment(v) || is.function(v) || is.matrix(v) || typeof(v) == "externalptr")
      violations <<- c(violations, "noncompact reference or matrix")
    if (is.list(v)) {
      violations <<- c(violations, intersect(names(v), c("ref", "frame_ref", "panel_ref", "fitted_nuisance_ref", "H", "H_chol", "qr_R",
        "fit_weights", "a_obs", "b_obs", "prob_use", "nuisance_registry", "fitted_views")))
      for (a in v) visit(a)
    }
  }
  visit(z)
  expect_identical(violations, character())
  expect_lt(length(serialize(z, NULL)), 30000L)
  # Owner contents and excluded registry payload never enter the ledger.
  x$res$frame_owner$unrelated <- rep(1, 10000)
  x$res$nuisance_registry$unrelated <- rep(1, 10000)
  expect_identical(nd_source$chained_request_nuisance_diagnostics(x$res, x$control, FALSE), z)
  expect_gt(z$accounting$retained_entries_bound, length(serialize(z, NULL))/8)
  expect_gte(z$accounting$simultaneous_entries_bound, z$accounting$retained_entries_bound + z$accounting$scratch_entries_bound + 2000)
})

test_that("unavailable cases and diagnostic failure leave fit data unchanged", {
  x <- nd_fixture()
  x$res$auxiliary_reporting$fitted_views <- list()
  custom <- nd_source$chained_request_nuisance_diagnostics(x$res, x$control, TRUE)
  expect_identical(custom$status, "unavailable")
  expect_identical(custom$reason, "external_custom_callback_no_builtin_logit_information")
  expect_length(custom$records, 0L)
  expect_identical(custom$post_links[[1]]$reason, custom$reason)
  empty <- nd_source$chained_request_nuisance_diagnostics(x$res, x$control, FALSE)
  expect_identical(empty$reason, "no_supported_pair_no_fitted_nuisance")
  expect_identical(empty$failed_fit_rows, "not returned after failed fit")
  x$res$auxiliary_reporting <- NULL; x$res$point_links_all <- list()
  empty <- nd_source$chained_request_nuisance_diagnostics(x$res, x$control, FALSE)
  expect_identical(empty$status, "unavailable")
  expect_length(empty$post_links, 0L)
  expect_length(empty$auxiliary_links, 0L)
  expect_length(empty$records, 0L)
  x <- nd_fixture()
  x$res$auxiliary_reporting$fitted_views[[4]]$factors$H[,] <- NA_real_
  z <- nd_source$chained_request_nuisance_diagnostics(x$res, x$control, FALSE)
  expect_identical(z$records[[3]]$status, "available")
  expect_true(z$records[[3]]$fit$converged)
  expect_identical(z$records[[3]]$conditioning$rcond_H_1_estimate$status, "unavailable")
  # The output cap applies to the ledger, separately from pre-existing payload.
  x <- nd_fixture()
  x$res$resource_accounting$simultaneous_entries <- 2e6
  x$control$owner_control$max_output_entries <- 1e6
  z <- nd_source$chained_request_nuisance_diagnostics(x$res, x$control, FALSE)
  expect_lt(z$accounting$retained_entries_bound, 1e6)
  for (cap in c("max_retained_entries", "max_matrix_entries", "max_output_entries", "max_work_ops")) {
    x <- nd_fixture(); x$control$owner_control[[cap]] <- 1
    expect_error(nd_source$chained_request_nuisance_diagnostics(x$res, x$control, FALSE), class = "cdid_ipw_owner_bound_error")
  }
})

test_that("capture, release, accounting and actual public binding sites stay explicit", {
  fit <- paste(deparse(body(nd_source$fit_ipw_nuisance_view)), collapse = "\n")
  positions <- vapply(c("st$fit_weights", "st$a_obs", "st$b_obs", "st$prob_use", "st <- NULL"),
    function(s) regexpr(s, fit, fixed = TRUE)[1], integer(1))
  expect_true(all(positions > 0L))
  expect_true(all(positions[1:4] < positions[5]))
  expect_match(fit, "cleanup_temporary_buffers()", fixed = TRUE)
  expect_match(fit, "preflight$nuisance_diagnostics$work_ops_bound", fixed = TRUE)
  pre <- paste(deparse(body(nd_source$preflight_ipw_fitted_nuisance_view)), collapse = "\n")
  expect_match(pre, "D27 three Kish summaries and fitted probability range", fixed = TRUE)
  expect_match(pre, "scratch_entries = 8 * min(B_eff", fixed = TRUE)
  binder <- paste(deparse(body(nd_source$chained_request_bind_mp)), collapse = "\n")
  expect_match(binder, "mp$chained$nuisance_diagnostics <- chained_request_nuisance_diagnostics(res,", fixed = TRUE)
  expect_match(binder, "mp$chained$covariance_handling <- res$covariance_handling", fixed = TRUE)
  condition <- paste(deparse(body(nd_source$chained_request_H_diagnostics)), collapse = "\n")
  expect_match(condition, 'base::rcond(H, norm = "O")', fixed = TRUE)
})
