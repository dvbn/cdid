# SPDX-License-Identifier: GPL-3.0-only
#
# Focused D5 tests for the support-preserving CP02 private pipeline
# (`R/chained_ipw_pipeline.R`): full declared ledger, catalog-relative NA cells,
# estimable post-contrast basis, supported-subvector inference/certificates,
# direct-diagnostic union, and per-aggregate `na.rm` behaviour.
#
# Reference values: `plan/evidence/cp02_d5_design/rational_results.json` (six
# independent Fraction recipes) and `verification/independent/d5_early_pair_reference.json`
# (raw early-pair case, derived before the implementation ran). Exact rationals
# are hard-coded below; the optional workspace files are cross-checked when present.
# Cases 1, 2, 4 and 5 are CONDITIONAL adapter/solver-handoff references with
# declared designs and complete IF inputs: they exercise the D5 composition
# helpers and the accepted core, not raw-data moment validity. Cases 3, 6, the
# early-pair case and the rowless/empty/permutation cases are raw-panel checks.
# Deterministic throughout: no RNG draws; simultaneous probes use fixed signs.

d5_rational <- function(x) {
  vapply(as.character(x), function(s) {
    parts <- strsplit(s, "/", fixed = TRUE)[[1L]]
    if (length(parts) == 1L) as.double(parts[1L]) else as.double(parts[1L]) / as.double(parts[2L])
  }, numeric(1L), USE.NAMES = FALSE)
}

make_d5_pipeline_control <- function(require_matching_scale = TRUE, max_target_error_abs = 1e-4,
                                     max_loading_error_abs = 1e-4, max_work_ops = 500000000L) {
  lnk_ctrl <- ipw_link_enclosure_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
    max_matrix_entries = 2000000L, max_work_ops = max_work_ops, max_output_entries = 500000L)
  ipw_pipeline_control(
    owner_control = ipw_owner_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = max_work_ops, max_output_entries = 500000L),
    catalog_bounds = list(max_candidates = 500L, max_rows = 500L, max_pool_entries = 5000L,
      max_scratch = 100000L, max_work = max_work_ops),
    class_bounds = list(max_links = 500L, max_classes = 50L, max_class_dates = 50L,
      max_E_entries = 5000L, max_mask_entries = 100000L),
    fitted_control = ipw_fitted_view_control(max_chunk_rows = 64L, max_observer_units = 2000L,
      max_retained_entries = 2000000L, max_matrix_entries = 2000000L, max_work_ops = max_work_ops,
      max_output_entries = 500000L),
    point_control = ipw_point_view_control(max_chunk_rows = 64L, max_observer_units = 2000L,
      max_retained_entries = 2000000L, max_matrix_entries = 2000000L, max_work_ops = max_work_ops,
      max_output_entries = 500000L),
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    certificate_control = logit_certificate_control(taylor_degree = 12L, max_abs_eta = 50.0,
      max_squarings = 16L, max_chunk_rows = 64L, max_matrix_entries = 2000000L, max_work_ops = max_work_ops),
    problem_control = ipw_nuisance_problem_control(max_chunk_rows = 64L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = max_work_ops),
    link_control = lnk_ctrl,
    covariance_control = ipw_covariance_transfer_control(max_chunk_rows = 64L, max_operator_block_rows = 64L,
      max_cluster_rows = 2000L, max_matrix_entries = 2000000L, max_work_ops = max_work_ops,
      max_output_entries = 1000000L, link_control = lnk_ctrl),
    gmm_accuracy = gmm_accuracy_control(method = "certified_eigen_chol_kkt", target_error_unit = "abs",
      max_target_error_abs = max_target_error_abs, max_target_error_se = 1e-2, max_orthogonality_defect = 1e-8,
      min_certified_lambda_rel = 1e-10, max_factor_residual_rel = 1e-10, min_rcond_R = 1e-10, min_rcond_U = 1e-10,
      min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10, max_loading_residual_rel = 1e-10,
      require_matching_scale = require_matching_scale, max_loading_error_abs = max_loading_error_abs),
    gmm_bounds = gmm_bounds_control(max_matrix_entries = 2000000L, max_batch_rows = 64L, max_targets = 100L,
      max_open_clusters = 2000L, allow_materialize_if = TRUE, max_materialized_if_entries = 1000000L),
    multiplier_budget = list(max_units = 2000, max_clusters = 2000, max_targets = 100, max_draws = 2000,
      max_if_batch_rows = 64, max_if_batch_cells = 6400, max_sign_batch_rows = 32,
      max_sign_batch_cells = 64000, max_workspace_cells = 1000000, max_workspace_bytes = 8000000),
    inference_accuracy = list(max_estimate_error_abs = 1e-3, max_se_error_abs = 1e-2,
      max_vcov_error_abs = 1e-2, max_multiplier_error_abs = 1e-1, max_band_error_abs = 1e-1))
}

d5_root_witness <- function(fv) list(F = unname(fv$factors$H_chol), r = 1e-6)

d5_locate_fixture <- function(filename) {
  candidates <- c(tryCatch(testthat::test_path("fixtures", filename), error = function(e) character(0L)),
    file.path("fixtures", filename), file.path("tests", "testthat", "fixtures", filename),
    file.path("..", "..", "tests", "testthat", "fixtures", filename))
  for (p in candidates) if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  stop("Fixture file not found: ", filename)
}

d5_locate_optional <- function(rel_path) {
  for (p in c(rel_path, file.path("..", "..", rel_path), file.path("..", "..", "..", rel_path))) {
    if (file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  }
  NULL
}

# Temporary binding replacement for the two injection tests below. Against the
# installed package the enclosing environment is the locked namespace, so a
# plain assign() errors ("cannot change value of locked binding") although it
# works in the unlocked source harness. The binding is unlocked only for the
# swap and re-locked immediately; callers restore the original the same way.
d5_replace_binding <- function(name, value, env) {
  locked <- bindingIsLocked(name, env)
  if (locked) unlockBinding(name, env)
  assign(name, value, envir = env)
  if (locked) lockBinding(name, env)
  invisible(locked)
}

# Adapter-level solver handoff: a declared calendar design, exact moments, an
# explicit IF matrix (units x links) and the D5 ledger/basis target construction.
d5_core_handoff <- function(design, ledger, mhat, Phi, ctrl, weighting = "optimal", target_order = NULL) {
  n <- nrow(Phi); L <- ncol(Phi)
  V <- crossprod(sweep(Phi, 2L, colMeans(Phi), "-")) / n
  dimnames(V) <- list(design$incidence$link_ids, design$incidence$link_ids)
  names(mhat) <- design$incidence$link_ids
  ledger <- ipw_pipeline_ledger_support(ledger, NULL, design)
  basis <- ipw_pipeline_post_basis(design)
  params <- design$parameters
  est <- params[params$level_estimable, , drop = FALSE]
  unavailable_rows <- which(!ledger$available & !is.na(ledger$design_param_index))
  specs <- c(
    lapply(seq_len(nrow(est)), function(j) gmm_target_spec(est$param_id[j], as.integer(est$param_index[j]), 1.0)),
    lapply(unavailable_rows, function(i) gmm_target_spec(ledger$cell_id[i], as.integer(ledger$design_param_index[i]), 1.0)),
    lapply(seq_len(nrow(basis)), function(b) gmm_target_spec(basis$basis_id[b], c(basis$from_param_index[b], basis$to_param_index[b]), c(-1.0, 1.0))))
  if (!is.null(target_order)) specs <- specs[target_order]
  prov <- gmm_covariance_provenance(row_level = "unit", n_units = n, n_clusters = n, n_positive_contributing_rows = n,
    centered = TRUE, weighted = FALSE, weight_type = "unweighted", normalization_factor = 1 / n,
    scale_id = "d5:adapter", nuisance_fit_id = "d5:declared_if", reduction = list(mode = "none"))
  replay <- gmm_if_replay_spec(callback = function(consume) {
    consume(list(unit_id = seq_len(n), phi = structure(Phi, dimnames = list(NULL, design$incidence$link_ids)),
      weight = rep(1, n), nuisance_fit_id = "d5:declared_if")) },
    input_row_level = "unit", aggregate_to_clusters = FALSE, n_input_rows = n, n_output_rows = n,
    n_positive_contributing_rows = n, centered = TRUE, weighted = FALSE, weight_type = "unweighted",
    normalization_factor = 1 / n, scale_id = "d5:adapter", nuisance_fit_id = "d5:declared_if", expected_unit_ids = seq_len(n))
  acc <- ctrl$gmm_accuracy; acc$require_matching_scale <- identical(weighting, "optimal")
  fit <- solve_chained_gmm_core(design = design, mhat = mhat, V = V, provenance = prov, targets = specs,
    weighting = weighting, if_replay = replay, accuracy = acc, bounds = ctrl$gmm_bounds,
    gauge_rule = "lowest", request_param_cov = TRUE, materialize_target_if = TRUE)
  list(fit = fit, ledger = ledger, basis = basis, V = V, specs = specs)
}

test_that("case 1 (conditional adapter): anchored, unanchored and isolated requests share one joint GMM; the basis keeps c-b and its cross-covariance with a", {
  ctrl <- make_d5_pipeline_control()
  # Cohort 2, reference 1, calendar 1..5: a=ATT(2,2), b=ATT(2,3), c=ATT(2,4), z=ATT(2,5).
  # Two distinct valid rows join the reference to a (different operators/pools); one row joins b to c.
  links <- data.frame(g = 2, s = c(1, 1, 3), t = c(2, 2, 4), ref = 1,
    functional_key = c("op:K{Inf}", "op:K{6,Inf}", "op:K{Inf}"), control_pool = c("{Inf}", "{6,Inf}", "{Inf}"),
    s_qualified = c(TRUE, TRUE, TRUE), t_qualified = TRUE, stringsAsFactors = FALSE)
  design <- build_calendar_chained_design(links, calendar = 1:5, cohorts = c(2, 6, Inf), anticipation = 0,
    target_post_times = data.frame(g = 2, t = 2:5))
  ledger <- ipw_pipeline_request_ledger(2, 1:5, c(2, 6, Inf), 0, "nevertreated", target_post_times = data.frame(g = 2, t = 2:5))
  h1 <- c(1, -1, 1, -1, 1, -1, 1, -1); h2 <- c(1, 1, -1, -1, 1, 1, -1, -1); h3 <- c(1, 1, 1, 1, -1, -1, -1, -1); h4 <- h1 * h2
  Phi <- cbind(h1 + h3, h2 + h4, h3)
  z <- d5_core_handoff(design, ledger, c(2, 4, 7), Phi, ctrl)
  expect_equal(unname(z$V), rbind(c(2, 0, 1), c(0, 2, 0), c(1, 0, 1)), tolerance = 0)
  expect_identical(z$ledger$support_status, c("estimable", "disconnected_component", "disconnected_component", "isolated_unsupported"))
  expect_identical(z$ledger$available, c(TRUE, FALSE, FALSE, FALSE))
  expect_false(any(grepl("nonidentif", z$ledger$support_reason)))
  expect_true(all(grepl("catalog-relative", z$ledger$population_claim)))
  expect_equal(nrow(z$basis), 1L)
  expect_identical(z$basis$from_id, "ATT(g=2,t=3)"); expect_identical(z$basis$to_id, "ATT(g=2,t=4)")
  tg <- z$fit$targets
  expect_identical(tg$target_id, c("ATT(g=2,t=2)", "ATT(g=2,t=3)", "ATT(g=2,t=4)", "ATT(g=2,t=5)", "contrast:ATT(g=2,t=4)-ATT(g=2,t=3)"))
  expect_identical(tg$estimable, c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(tg$status[2:4], c("disconnected_component", "disconnected_component", "isolated_unsupported"))
  expect_equal(tg$estimate[c(1L, 5L)], c(3, 15 / 2), tolerance = 1e-12)
  expect_true(all(is.na(tg$estimate[2:4])))
  # Loadings A = [1/2 1/2 0; -1/4 1/4 1]; Omega = [[1,1/2],[1/2,3/4]]; estimator covariance Omega/8.
  expect_equal(unname(t(z$fit$internal_gauge$Lambda)), rbind(c(1 / 2, 1 / 2, 0), c(-1 / 4, 1 / 4, 1)), tolerance = 1e-12)
  expect_equal(unname(z$fit$sandwich_cov), rbind(c(1, 1 / 2), c(1 / 2, 3 / 4)), tolerance = 1e-12)
  Psi <- z$fit$if_projection$materialized_target_if
  expect_equal(unname(crossprod(Psi) / 64), rbind(c(1 / 8, 1 / 16), c(1 / 16, 3 / 32)), tolerance = 1e-12)
  expect_equal(unname(Psi), unname(Phi %*% t(rbind(c(1 / 2, 1 / 2, 0), c(-1 / 4, 1 / 4, 1)))), tolerance = 1e-12)
  # Unanchored levels are never reported from gauge coordinates.
  expect_true(all(is.na(z$fit$parameters$level_estimate[!z$fit$parameters$level_estimable])))
  # Discriminators from the executed probe: m1 + 1/8 moves (a, c-b) by (1/16, -1/32);
  # zeroing the cross-component covariance would give c-b = 7 (blockwise fit is not licensed).
  zp <- d5_core_handoff(design, ledger, c(2 + 1 / 8, 4, 7), Phi, ctrl)
  expect_equal(zp$fit$targets$estimate[c(1L, 5L)] - tg$estimate[c(1L, 5L)], c(1 / 16, -1 / 32), tolerance = 1e-12)
  Vb <- diag(c(2, 2, 1)); dimnames(Vb) <- dimnames(z$V)
  W <- rbind(c(1, 0), c(1, 0), c(0, 1)); m <- c(2, 4, 7)
  blockwise <- solve(t(W) %*% solve(Vb) %*% W, t(W) %*% solve(Vb) %*% m)
  expect_equal(blockwise[2L], 7, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(blockwise[2L], tg$estimate[5L])))
  # Permuting the target order only permutes the core table; loading columns stay aligned by identity.
  perm <- c(5L, 3L, 1L, 4L, 2L)
  zq <- d5_core_handoff(design, ledger, c(2, 4, 7), Phi, ctrl, target_order = perm)
  expect_identical(zq$fit$targets$target_id, tg$target_id[perm])
  expect_equal(zq$fit$targets$estimate, tg$estimate[perm], tolerance = 1e-12)
  expect_equal(unname(zq$fit$if_projection$materialized_target_if[, "ATT(g=2,t=2)"]), unname(Psi[, "ATT(g=2,t=2)"]), tolerance = 1e-12)
  expect_equal(unname(zq$fit$if_projection$materialized_target_if[, 1L]), unname(Psi[, 2L]), tolerance = 1e-12)
  ref <- d5_locate_optional("plan/evidence/cp02_d5_design/rational_results.json")
  if (!is.null(ref)) {
    rr <- jsonlite::fromJSON(ref)$case1
    expect_equal(tg$estimate[c(1L, 5L)], d5_rational(unlist(rr$compact_estimates)), tolerance = 1e-12)
    expect_equal(unname(z$fit$sandwich_cov), matrix(d5_rational(unlist(t(rr$Omega))), 2L, byrow = TRUE), tolerance = 1e-12)
  }
})

test_that("case 2 (conditional adapter): an estimable post contrast survives with zero anchored post levels", {
  ctrl <- make_d5_pipeline_control()
  links <- data.frame(g = 3, s = 1, t = c(3, 4), ref = 2, functional_key = "op", control_pool = "{Inf}",
    s_qualified = FALSE, t_qualified = TRUE, stringsAsFactors = FALSE)
  design <- build_calendar_chained_design(links, calendar = 1:4, cohorts = c(3, Inf), anticipation = 0,
    target_post_times = data.frame(g = 3, t = 3:4))
  ledger <- ipw_pipeline_request_ledger(3, 1:4, c(3, Inf), 0, "nevertreated")
  hp <- rbind(c(1, 1, 1), c(-1, 1, -1), c(1, -1, -1), c(-1, -1, 1))
  Phi <- cbind(hp[, 1L] + hp[, 3L], hp[, 2L] + hp[, 3L])
  z <- d5_core_handoff(design, ledger, c(3, 7), Phi, ctrl)
  expect_equal(unname(z$V), rbind(c(2, 1), c(1, 2)), tolerance = 0)
  expect_identical(z$ledger$support_status, c("disconnected_component", "disconnected_component"))
  expect_equal(design$rank_diagnostics$n_estimable_post_targets, 0L)
  tg <- z$fit$targets
  expect_identical(tg$estimable, c(FALSE, FALSE, TRUE))
  expect_equal(tg$estimate[3L], 4, tolerance = 1e-12)
  expect_equal(unname(z$fit$internal_gauge$Lambda[, 1L]), c(-1, 1), tolerance = 1e-12)
  expect_equal(unname(z$fit$sandwich_cov), matrix(2), tolerance = 1e-12)
  expect_equal(unname(z$fit$if_projection$materialized_target_if[, 1L]), c(0, 2, -2, 0), tolerance = 1e-12)
  expect_equal(sum(z$fit$if_projection$materialized_target_if[, 1L]^2) / 16, 1 / 2, tolerance = 1e-12)
  # Internal gauge coordinates (root contrasts b=3, c=7) never become level estimates.
  expect_equal(unname(z$fit$internal_gauge$theta_gauge_coords), c(3, 7), tolerance = 1e-12)
  expect_true(all(is.na(z$fit$parameters$level_estimate)))
  expect_true(all(is.na(tg$estimate[1:2])))
  zp <- d5_core_handoff(design, ledger, c(3, 7 + 1 / 8), Phi, ctrl)
  expect_equal(zp$fit$targets$estimate[3L], 33 / 8, tolerance = 1e-12)
  zi <- d5_core_handoff(design, ledger, c(3, 7), Phi, ctrl, weighting = "identity")
  expect_equal(zi$fit$targets$estimate[3L], 4, tolerance = 1e-12)
})

test_that("case 5 (conditional adapter): singular weight covariance refuses despite a structurally estimable target, with diagnostics", {
  ctrl <- make_d5_pipeline_control()
  links <- data.frame(g = 2, s = 1, t = 2, ref = 1, functional_key = c("opA", "opB"), control_pool = c("{Inf}", "{6,Inf}"),
    s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE)
  design <- build_calendar_chained_design(links, calendar = 1:2, cohorts = c(2, 6, Inf), anticipation = 0,
    target_post_times = data.frame(g = 2, t = 2))
  expect_true(design$post_targets$level_estimable)
  ledger <- ipw_pipeline_request_ledger(2, 1:2, c(2, 6, Inf), 0, "nevertreated")
  Phi <- rbind(c(1, 0), c(-1, 0), c(1, 0), c(-1, 0))
  err <- tryCatch(d5_core_handoff(design, ledger, c(1, 3), Phi, ctrl), error = function(e) e)
  expect_s3_class(err, "cdid_gmm_pd_error")
  expect_true(!is.null(err$diagnostics))
  expect_false(inherits(err, "cdid_gmm_target_estimability_error"))
  # The same refusal does not depend on an unrelated unsupported request.
  ledger2 <- ipw_pipeline_request_ledger(2, 1:2, c(2, 6, Inf), 0, "nevertreated", target_post_times = data.frame(g = 2, t = c(2, 2)))
  expect_equal(nrow(ledger2), 1L)
})

d5_case4_inputs <- function(frame_g = c(2, 2, 3, Inf)) {
  th <- c("ATT(g=2,t=2)" = 2, "ATT(g=2,t=3)" = NA_real_, "ATT(g=3,t=3)" = 6)
  X <- cbind(c(1, -1, 0, 0), NA_real_, c(0, 1, -1, 0)); colnames(X) <- names(th)
  meta <- data.frame(g = c(2, 2, 3), time = c(2, 3, 3), role = "post", target_id = names(th),
    support_status = c("estimable", "isolated_unsupported", "estimable"), stringsAsFactors = FALSE)
  list(th = th, X = X, meta = meta, g = frame_g, w = rep(1, 4))
}

test_that("case 4 (conditional adapter): per-aggregate partial support under na.rm = FALSE and the explicit restricted target under TRUE", {
  ctrl <- make_d5_pipeline_control(); z4 <- d5_case4_inputs()
  agg <- function(type, na_rm, frame_g = z4$g) aggregate_chained_ipw_targets(z4$th, z4$X, z4$meta, frame_g, z4$w, type = type,
    multiplier_budget = ctrl$multiplier_budget, inference_accuracy = ctrl$inference_accuracy, na.rm = na_rm)
  # na.rm = FALSE: aggregates needing B are NA and name B; unrelated aggregates survive.
  s0 <- agg("simple", FALSE)
  expect_true(is.na(s0$estimates[["ATT_simple"]])); expect_match(s0$cell_status$reason, "ATT\\(g=2,t=3\\)")
  expect_identical(s0$inference$status, "empty_supported_selection"); expect_false(isTRUE(s0$numerical_certificate$certified))
  expect_true(is.na(s0$overall$estimate)); expect_identical(s0$overall$status, "unavailable")
  c0 <- agg("calendar", FALSE)
  expect_equal(unname(c0$estimates), c(2, NA), tolerance = 0); expect_match(c0$cell_status$reason[2L], "ATT\\(g=2,t=3\\)")
  expect_true(is.na(c0$overall$estimate)); expect_match(c0$overall$reason, "calendar_3")
  g0 <- agg("group", FALSE)
  expect_equal(unname(g0$estimates), c(NA, 6), tolerance = 0); expect_true(is.na(g0$overall$estimate))
  d0 <- agg("dynamic", FALSE)
  expect_equal(unname(d0$estimates), c(10 / 3, NA), tolerance = 1e-12)
  expect_identical(d0$cell_status$available, c(TRUE, FALSE)); expect_true(is.na(d0$overall$estimate))
  expect_match(d0$overall$reason, "event_1")
  # na.rm = TRUE: simple restricts to {A, C} with displayed weights, exclusions and full share uncertainty.
  s1 <- agg("simple", TRUE)
  expect_equal(s1$estimates[["ATT_simple"]], 10 / 3, tolerance = 1e-12)
  expect_true(s1$cell_status$restricted); expect_identical(s1$cell_status$excluded_ids, "ATT(g=2,t=3)")
  expect_identical(s1$cell_used_target_ids[[1L]], c("ATT(g=2,t=2)", "ATT(g=3,t=3)"))
  expect_equal(unname(s1$weights_theta[1L, c(1L, 3L)]), c(2 / 3, 1 / 3), tolerance = 1e-12)
  expect_equal(unname(s1$cohort_shares), c(1 / 2, 1 / 4), tolerance = 0)
  expect_equal(unname(s1$gradient_shares[1L, ]), c(-16 / 9, 32 / 9), tolerance = 1e-12)
  expect_equal(unname(s1$if_fixed_full[, 1L]), c(2 / 3, -1 / 3, -1 / 3, 0), tolerance = 1e-12)
  expect_equal(unname(s1$if_share_full[, 1L]), c(-16 / 9, -16 / 9, 32 / 9, 0), tolerance = 1e-12)
  expect_equal(unname(s1$if_full[, 1L]), c(-10 / 9, -19 / 9, 29 / 9, 0), tolerance = 1e-12)
  expect_equal(unname(s1$inference$vcov[1L, 1L]), 217 / 216, tolerance = 1e-12)
  expect_equal(unname(s1$se_fixed_shares_only)^2, 1 / 24, tolerance = 1e-12)
  expect_equal(unname(s1$decomposition$Omega_fixed)[1L, 1L] / 4, 1 / 24, tolerance = 1e-12)
  expect_equal(unname(s1$decomposition$Omega_share)[1L, 1L] / 4, 32 / 27, tolerance = 1e-12)
  expect_equal(unname(s1$decomposition$Omega_cross)[1L, 1L] / 4, -2 / 9, tolerance = 1e-12)
  expect_equal(unname(s1$decomposition$Omega_total)[1L, 1L] / 4, 217 / 216, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(s1$inference$vcov[1L, 1L]), 265 / 216)))
  expect_match(s1$target_definition, "restricted")
  g1 <- agg("group", TRUE)
  expect_equal(unname(g1$estimates), c(2, 6), tolerance = 0)
  expect_equal(g1$overall$estimate, 10 / 3, tolerance = 1e-12); expect_true(g1$overall$share_dependent)
  expect_equal(unname(g1$overall$if_full), c(-10 / 9, -19 / 9, 29 / 9, 0), tolerance = 1e-12)
  c1 <- agg("calendar", TRUE)
  expect_equal(unname(c1$estimates), c(2, 6), tolerance = 0)
  expect_equal(c1$overall$estimate, 4, tolerance = 1e-12)
  expect_equal(unname(c1$overall$if_full), c(1 / 2, 0, -1 / 2, 0), tolerance = 1e-12)
  expect_equal(c1$overall$se^2, 1 / 32, tolerance = 1e-12)
  expect_identical(c1$cell_status$restricted, c(FALSE, TRUE))
  d1 <- agg("dynamic", TRUE)
  expect_equal(unname(d1$estimates), c(10 / 3, NA), tolerance = 1e-12)
  expect_match(d1$cell_status$reason[2L], "empty_restricted_set")
  # did 2.5.1 explicit na.rm = TRUE convention: the empty horizon is omitted and the surviving outer summary renormalized.
  expect_equal(d1$overall$estimate, 10 / 3, tolerance = 1e-12)
  expect_identical(d1$overall$used_bins, "event_0"); expect_identical(d1$overall$excluded_bins, "event_1")
  expect_identical(d1$overall$intended_bins, c("event_0", "event_1"))
  # Metadata variant: unit 4 has an unknown cohort; share-dependent targets explain the missing share.
  s2 <- agg("simple", TRUE, frame_g = c(2, 2, 3, NA))
  expect_true(is.na(s2$estimates[["ATT_simple"]])); expect_match(s2$cell_status$reason, "missing_frame_cohort_share")
  expect_false(s2$shares_known); expect_null(s2$if_shares_full)
  g2 <- agg("group", TRUE, frame_g = c(2, 2, 3, NA))
  expect_equal(unname(g2$estimates), c(2, 6), tolerance = 0)
  expect_identical(g2$overall$status, "unavailable"); expect_match(g2$overall$reason, "missing_frame_cohort_share")
  c2 <- agg("calendar", TRUE, frame_g = c(2, 2, 3, NA))
  expect_equal(unname(c2$estimates), c(2, 6), tolerance = 0)
  expect_equal(c2$overall$estimate, 4, tolerance = 1e-12); expect_equal(c2$overall$se^2, 1 / 32, tolerance = 1e-12)
  d2 <- agg("dynamic", TRUE, frame_g = c(2, 2, 3, NA))
  expect_match(d2$cell_status$reason[1L], "missing_frame_cohort_share")
  expect_error(aggregate_chained_ipw_targets(z4$th, z4$X, z4$meta, z4$g, c(1, 1, 0, 1), type = "simple",
    multiplier_budget = ctrl$multiplier_budget, inference_accuracy = ctrl$inference_accuracy), class = "cdid_ipw_pipeline_aggregation_error")
  ref <- d5_locate_optional("plan/evidence/cp02_d5_design/rational_results.json")
  if (!is.null(ref)) {
    rr <- jsonlite::fromJSON(ref)$case4
    expect_equal(unname(s1$if_full[, 1L]), d5_rational(rr$total_IF), tolerance = 1e-12)
    expect_equal(unname(s1$inference$vcov[1L, 1L]), d5_rational(rr$variance_total), tolerance = 1e-12)
    expect_equal(c1$overall$estimate, d5_rational(rr$restricted_calendar_overall), tolerance = 1e-12)
  }
})

d5_case3_panel <- function() {
  d <- expand.grid(id = 1:8, t = 1:3); d$g <- ifelse(d$id <= 4, 3, Inf)
  d$y <- NA_real_; d$y[d$t == 2] <- 0
  d$y[d$t == 3] <- c(1, 1, 3, 3, -1, -1, 1, 1)[d$id[d$t == 3]]
  d
}

test_that("case 3 (raw panel): an unsupported direct pretrend leaves the supported post output, its IF and its certificate available in both base-period modes", {
  ctrl <- make_d5_pipeline_control(); d <- d5_case3_panel()
  res <- lapply(c("universal", "varying"), function(bp) run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl,
    target_cohorts = 3, control_group = "nevertreated", weighting = "optimal", base_period = bp,
    run_numerical_transfer = TRUE, root_witness = d5_root_witness))
  names(res) <- c("universal", "varying")
  for (bp in names(res)) {
    r <- res[[bp]]; tm <- r$reporting$target_meta
    expect_identical(r$support_status, "partial_support")
    post <- which(tm$role == "post")
    expect_equal(unname(r$reporting$estimates[post]), 2, tolerance = 1e-12)
    expect_equal(unname(r$reporting$if_full[, post]), c(-2, -2, 2, 2, 2, 2, -2, -2), tolerance = 1e-12)
    expect_equal(unname(r$reporting$se[post])^2, 1 / 2, tolerance = 1e-12)
    expect_equal(unname(r$reporting$Omega_a[post, post]), 4, tolerance = 1e-12)
    expect_true(r$numerical_transfer$gmm_transfer$supported_subvector_certified)
    expect_false(r$numerical_transfer$gmm_transfer$whole_request_complete)
    expect_false(r$numerical_transfer$gmm_transfer$complete_inference_certified)
    expect_true("ATT(g=3,t=3)" %in% r$reporting$numerical_certificate$certified_target_ids)
    expect_identical(r$reporting$numerical_certificate$band_family, "pointwise_normal_over_certified_ids")
    pre <- which(tm$role == "pre")
    expect_length(pre, 1L); expect_true(is.na(r$reporting$estimates[pre])); expect_true(all(is.na(r$reporting$if_full[, pre])))
    expect_identical(tm$support_status[pre], "unavailable_pretrend")
    expect_match(tm$support_reason[pre], "zero empirical risk-set observations")
    expect_identical(tm$target_id[pre], if (bp == "universal") "pre(g=3,t=1)" else "pre(g=3,t=2)")
    expect_true(tm$target_id[pre] %in% r$reporting$numerical_certificate$unavailable_target_ids)
    ds <- r$auxiliary_reporting$descriptor_status
    expect_equal(nrow(ds), 1L); expect_identical(ds$status, "unsupported_pair"); expect_equal(c(ds$s, ds$t), c(1, 2)); expect_identical(ds$K, "{Inf}")
    cnt <- r$auxiliary_reporting$unsupported_counts[[1L]]
    expect_true(cnt$counts_available); expect_equal(unname(cnt$by_cohort_pair_count), c(0L, 0L)); expect_equal(cnt$n_observed_pair_all, 0L)
    expect_equal(r$auxiliary_reporting$post_GMM_rows_added, 0L)
    expect_equal(r$aggregations$simple$estimates[["ATT_simple"]], 2, tolerance = 1e-12)
    expect_true(r$aggregations$simple$numerical_certificate$certified)
    expect_equal(r$aggregations$dynamic$cell_status$available, if (bp == "universal") c(FALSE, TRUE, TRUE) else c(FALSE, TRUE))
  }
  expect_true(any(res$universal$reporting$target_meta$role == "reference"))
  expect_equal(unname(res$universal$reporting$estimates[res$universal$reporting$target_meta$role == "reference"]), 0)
  expect_false(any(res$varying$reporting$target_meta$role == "reference"))
  expect_equal(res$universal$gmm_fit$targets, res$varying$gmm_fit$targets, tolerance = 0)
  expect_equal(res$universal$covariance_stream, res$varying$covariance_stream, tolerance = 0)
  expect_equal(length(res$universal$auxiliary_reporting$selections$varying), 1L)
  expect_equal(length(res$universal$auxiliary_reporting$selections$universal), 2L)
  # Genuine accuracy/resource failures are still refusals when a diagnostic is unsupported.
  expect_error(run_chained_ipw_pipeline(d, "id", "t", "y", "g", make_d5_pipeline_control(max_target_error_abs = 1e-18),
    target_cohorts = 3, control_group = "nevertreated", weighting = "optimal", base_period = "universal",
    run_numerical_transfer = TRUE, root_witness = d5_root_witness), class = "cdid_gmm_accuracy_error")
  expect_error(run_chained_ipw_pipeline(d, "id", "t", "y", "g", make_d5_pipeline_control(max_work_ops = 100L),
    target_cohorts = 3, control_group = "nevertreated", weighting = "optimal", base_period = "universal",
    run_numerical_transfer = FALSE), class = "cdid_ipw_owner_bound_error")
})

d5_case6_panel <- function(delta_first = 0) {
  d <- expand.grid(id = 1:8, t = 1:3); d$g <- ifelse(d$id <= 4, 2, Inf)
  i12 <- c(1 + delta_first, 1, 3, 3, -1, -1, 1, 1); i23 <- c(2, 4, 2, 4, -1, 1, -1, 1)
  d$y <- ifelse(d$t == 1, 0, ifelse(d$t == 2, i12[d$id], (i12 + i23)[d$id]))
  d
}

test_that("case 6 (raw panel): complete data is unchanged by the D5 composition, with the proved lossless reduction and whole-request certification", {
  ctrl <- make_d5_pipeline_control(); ctrl_id <- make_d5_pipeline_control(require_matching_scale = FALSE)
  d <- d5_case6_panel()
  ru <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 2, control_group = "nevertreated",
    weighting = "optimal", base_period = "universal", run_numerical_transfer = TRUE, root_witness = d5_root_witness)
  rv <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 2, control_group = "nevertreated",
    weighting = "optimal", base_period = "varying", run_numerical_transfer = TRUE, root_witness = d5_root_witness)
  expect_identical(ru$support_status, "all_requested_cells_supported")
  expect_true(ru$support$whole_request_complete); expect_equal(nrow(ru$support$basis), 0L)
  expect_true(ru$reduction_summary$reduction_active)
  expect_equal(ru$reduction_summary$n_links_all, 3L); expect_equal(ru$reduction_summary$n_links_used, 2L)
  # All original comparisons retain provenance. On the outcome-free route only
  # retained rows are fitted; compare E-expanded moments to the hand truth.
  rows <- ru$metadata$catalog$rows
  expect_equal(nrow(rows), 3L)
  expanded_m <- as.double(ru$reduction_summary$expanded_m)
  key <- paste(rows$s, rows$t)
  expect_equal(unname(expanded_m[match(c("1 2", "2 3", "1 3"), key)]), c(2, 3, 5), tolerance = 1e-12)
  expect_equal(unname(vapply(ru$point_links_used, `[[`, numeric(1L), "estimate")),
    expanded_m[ru$reduction_summary$kept_indices], tolerance = 1e-12)
  expect_identical(ru$reduction_summary$reduction_stage, "outcome_free_before_point_fits")
  expect_equal(ru$reduction_summary$n_rows_fitted, 2L)
  expect_identical(ru$reduction_summary$max_moment_expansion_diff, NA_real_)
  expect_identical(ru$reduction_summary$expansion_self_test,
    "not_evaluated: dropped rows are never fitted (Prop. 5/Cor. 6)")
  expect_length(ru$reduction_summary$dropped_indices, 1L)
  for (m in ru$reduction_summary$dropped_indices) {
    record <- ru$point_links_all[[m]]
    expect_s3_class(record, "cdid_ipw_reduced_link_record")
    expect_identical(names(record), c("link_id", "catalog_row", "link_view_ref", "nuisance_view_ref",
      "fitted_nuisance_ref", "expansion", "point_fit_performed", "outcomes_read", "status", "scope"))
    expect_null(record[["estimate", exact = TRUE]])
    expect_false(record$point_fit_performed); expect_false(record$outcomes_read)
    expect_identical(record$status, "structurally_reduced_outcome_free")
    expect_match(record$scope, "never fitted, replayed or read", fixed = TRUE)
  }
  # The accepted private reference fits EVERY row on this SAME fixture. Keep
  # the original separately-fitted expansion self-test and complete output check.
  rpost <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 2, control_group = "nevertreated",
    weighting = "optimal", base_period = "universal", run_numerical_transfer = TRUE, root_witness = d5_root_witness,
    structural_reduction_stage = "post_point_fit")
  expect_identical(rpost$reduction_summary$reduction_stage, "after_point_fits")
  expect_equal(rpost$reduction_summary$n_rows_fitted, 3L)
  expect_identical(rpost$reduction_summary$expansion_self_test, "evaluated_on_fitted_rows")
  all_m <- vapply(rpost$point_links_all, `[[`, numeric(1L), "estimate")
  expect_equal(unname(all_m[match(c("1 2", "2 3", "1 3"), key)]), c(2, 3, 5), tolerance = 1e-12)
  expect_equal(unname(as.double(rpost$reduction_summary$expanded_m)), unname(all_m), tolerance = 1e-12)
  expect_lt(rpost$reduction_summary$max_moment_expansion_diff, 1e-12)
  expect_equal(expanded_m, as.double(rpost$reduction_summary$expanded_m), tolerance = 1e-12)
  for (nm in c("support", "covariance_stream", "full_frame_transport", "reporting", "basis", "aggregations"))
    expect_equal(ru[[nm]], rpost[[nm]], tolerance = 1e-12, info = nm)
  for (nm in c("targets", "parameters", "internal_gauge", "sandwich_cov", "model_cov", "if_projection", "numerical_guard"))
    expect_equal(ru$gmm_fit[[nm]], rpost$gmm_fit[[nm]], tolerance = 1e-12, info = nm)
  expect_identical(ru$covariance_handling$handler_invoked, rpost$covariance_handling$handler_invoked)
  expect_equal(ru$covariance_handling$lambda, rpost$covariance_handling$lambda, tolerance = 1e-12)
  expect_equal(ru$numerical_transfer$gmm_transfer, rpost$numerical_transfer$gmm_transfer, tolerance = 1e-12)
  post <- which(ru$reporting$target_meta$role == "post")
  expect_equal(unname(ru$reporting$estimates[post]), c(2, 5), tolerance = 1e-12)
  expect_equal(unname(ru$reporting$Omega_a[post, post]), rbind(c(4, 4), c(4, 8)), tolerance = 1e-12)
  expect_equal(unname(ru$reporting$vcov[post, post]), rbind(c(1 / 2, 1 / 2), c(1 / 2, 1)), tolerance = 1e-12)
  f12 <- c(-2, -2, 2, 2, 2, 2, -2, -2); f23 <- c(-2, 2, -2, 2, 2, -2, 2, -2)
  expect_equal(unname(ru$reporting$if_full[, post]), unname(cbind(f12, f12 + f23)), tolerance = 1e-12)
  expect_equal(unname(ru$reporting$estimates[ru$reporting$target_meta$role == "reference"]), 0)
  expect_equal(ru$aggregations$simple$estimates[["ATT_simple"]], 7 / 2, tolerance = 1e-12)
  expect_equal(unname(ru$aggregations$simple$if_full[, 1L]), f12 + f23 / 2, tolerance = 1e-12)
  expect_equal(ru$aggregations$simple$se[["ATT_simple"]]^2, 5 / 8, tolerance = 1e-12)
  expect_equal(unname(ru$aggregations$group$estimates), 7 / 2, tolerance = 1e-12)
  expect_true(ru$numerical_transfer$gmm_transfer$complete_inference_certified)
  expect_true(ru$numerical_transfer$gmm_transfer$whole_request_complete)
  expect_identical(sort(ru$numerical_transfer$gmm_transfer$certified_target_ids), sort(ru$reporting$target_meta$target_id))
  expect_true(all(vapply(ru$aggregations, function(z) isTRUE(z$numerical_certificate$certified) && isTRUE(z$overall$numerical_certificate$certified), logical(1L))))
  expect_equal(rv$reporting$estimates[rv$reporting$target_meta$role == "post"], ru$reporting$estimates[post], tolerance = 1e-12)
  expect_equal(rv$gmm_fit$targets, ru$gmm_fit$targets, tolerance = 0)
  expect_false(any(rv$reporting$target_meta$role == "reference"))
  # Identity weighting keeps every original row and recovers the same exactly telescoping post levels.
  ri <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl_id, target_cohorts = 2, control_group = "nevertreated",
    weighting = "identity", base_period = "universal", run_numerical_transfer = FALSE)
  expect_false(ri$reduction_summary$reduction_active); expect_equal(ri$reduction_summary$n_links_used, 3L)
  expect_equal(unname(ri$reporting$estimates[ri$reporting$target_meta$role == "post"]), c(2, 5), tolerance = 1e-12)
  # Executed perturbation: +1/4 on the first treated increment moves both post targets by 1/16.
  rp <- run_chained_ipw_pipeline(d5_case6_panel(1 / 4), "id", "t", "y", "g", ctrl, target_cohorts = 2, control_group = "nevertreated",
    weighting = "optimal", base_period = "universal", run_numerical_transfer = FALSE)
  expect_equal(unname(rp$reporting$estimates[post] - ru$reporting$estimates[post]), c(1 / 16, 1 / 16), tolerance = 1e-12)
  acc <- ru$resource_accounting$d5_arrays
  expect_equal(acc$n_ledger_cells, 2L); expect_equal(acc$n_basis_contrasts, 0L); expect_equal(acc$ledger_if_entries, 16)
  expect_equal(acc$n_core_targets, 2L); expect_equal(acc$n_estimable_loading_columns, 2L)
  expect_lte(ru$resource_accounting$initial_envelope_entries, ru$resource_accounting$simultaneous_entries)
})

d5_early_pair_panel <- function(delta_unit3 = 0) {
  d <- expand.grid(id = 1:8, t = 1:4); d$g <- ifelse(d$id <= 4, 3, Inf)
  obs13 <- d$id %in% c(1, 2, 5, 6); obs14 <- d$id %in% c(3, 4, 7, 8)
  d13 <- c(1, 3, NA, NA, 0, 2, NA, NA); d14 <- c(NA, NA, 5 + delta_unit3, 7, NA, NA, 1, 3)
  d$y <- NA_real_; d$y[d$t == 1] <- 0
  d$y[d$t == 3 & obs13] <- d13[d$id[d$t == 3 & obs13]]; d$y[d$t == 4 & obs14] <- d14[d$id[d$t == 4 & obs14]]
  d
}

test_that("raw early common-control pairs {1,3}/{1,4} with the reference date unobserved: no anchored post level, one certified estimable post contrast", {
  ctrl <- make_d5_pipeline_control(); ctrl_id <- make_d5_pipeline_control(require_matching_scale = FALSE)
  d <- d5_early_pair_panel()
  res <- list()
  for (bp in c("universal", "varying")) for (wt in c("optimal", "identity")) {
    res[[paste(bp, wt)]] <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", if (wt == "optimal") ctrl else ctrl_id, target_cohorts = 3,
      control_group = "nevertreated", weighting = wt, base_period = bp, run_numerical_transfer = TRUE, root_witness = d5_root_witness)
  }
  for (nm in names(res)) {
    r <- res[[nm]]
    rows <- r$metadata$catalog$rows
    expect_equal(nrow(rows), 2L); expect_setequal(paste(rows$s, rows$t), c("1 3", "1 4"))
    expect_identical(unique(rows$s_role), "pre_private"); expect_equal(length(unique(rows$operator_id)), 1L)
    expect_equal(unname(r$covariance_stream$mhat), c(1, 4), tolerance = 1e-12)
    expect_equal(unname(r$covariance_stream$V), diag(c(8, 8)), tolerance = 1e-12)
    led <- r$support$ledger
    expect_identical(led$cell_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)"))
    expect_identical(led$support_status, c("disconnected_component", "disconnected_component"))
    expect_false(any(led$available)); expect_true(all(is.na(led$loading_column)))
    expect_false(any(grepl("nonidentif", led$support_reason)))
    expect_identical(r$support_status, "partial_support"); expect_equal(r$support$n_available, 0L)
    expect_equal(r$design_all$rank_diagnostics$n_estimable_post_targets, 0L)
    expect_true(all(is.na(r$gmm_fit$parameters$level_estimate)))
    tm <- r$reporting$target_meta; post <- which(tm$role == "post")
    expect_true(all(is.na(r$reporting$estimates[post]))); expect_true(all(is.na(r$reporting$if_full[, post])))
    expect_true(all(is.na(r$reporting$vcov[post, post])))
    expect_true(all(is.na(r$reporting$se[post])))
    b <- r$basis
    expect_identical(b$meta$basis_id, "contrast:ATT(g=3,t=4)-ATT(g=3,t=3)")
    expect_equal(unname(b$estimates), 3, tolerance = 1e-12)
    expect_equal(unname(b$if_full[, 1L]), c(4, -4, -4, 4, -4, 4, 4, -4), tolerance = 1e-12)
    expect_equal(unname(b$inference$Omega[1L, 1L]), 16, tolerance = 1e-12)
    expect_equal(unname(b$vcov[1L, 1L]), 2, tolerance = 1e-12)
    expect_equal(unname(b$se)^2, 2, tolerance = 1e-12)
    expect_true(b$numerical_certificate$certified)
    expect_identical(b$numerical_certificate$certified_target_ids, "contrast:ATT(g=3,t=4)-ATT(g=3,t=3)")
    expect_true(all(b$if_full[, 1L] >= b$numerical_certificate$rows$lo[, 1L] & b$if_full[, 1L] <= b$numerical_certificate$rows$hi[, 1L]))
    expect_equal(nrow(b$witnesses[[1L]]), 2L)
    expect_true(all(is.na(b$cross_vcov_with_reported_cells[post, 1L])))
    gt <- r$numerical_transfer$gmm_transfer
    expect_true(gt$supported_subvector_certified); expect_false(gt$whole_request_complete); expect_false(gt$complete_inference_certified)
    expect_true("contrast:ATT(g=3,t=4)-ATT(g=3,t=3)" %in% gt$certified_target_ids)
    expect_false(any(led$cell_id %in% gt$certified_target_ids))
    expect_setequal(gt$unavailable_target_ids, c(led$cell_id, tm$target_id[tm$role == "pre"]))
    for (a in r$aggregations) { expect_true(all(is.na(a$estimates[a$cell_status$n_used == 0L]))); expect_identical(a$overall$status, "unavailable") }
    expect_match(r$aggregations$simple$cell_status$reason, "ATT\\(g=3,t=3\\) \\[disconnected_component\\]")
    acc <- r$resource_accounting$d5_arrays
    expect_equal(acc$n_ledger_cells, 2L); expect_equal(acc$n_basis_contrasts, 1L); expect_equal(acc$n_core_targets, 3L)
    expect_equal(acc$n_estimable_loading_columns, 1L); expect_equal(acc$basis_if_entries, 8); expect_equal(acc$ledger_if_entries, 16)
    expect_equal(acc$basis_enclosure_entries, 16)
  }
  ru <- res[["universal optimal"]]
  expect_equal(unname(ru$reporting$estimates[ru$reporting$target_meta$role == "reference"]), 0)
  expect_identical(ru$reporting$target_meta$support_status[ru$reporting$target_meta$role == "pre"], "unavailable_pretrend")
  # na.rm = TRUE cannot manufacture an aggregate from an empty supported set.
  rn <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 3, control_group = "nevertreated",
    weighting = "optimal", base_period = "universal", run_numerical_transfer = FALSE, na.rm = TRUE)
  expect_true(is.na(rn$aggregations$simple$estimates[["ATT_simple"]])); expect_match(rn$aggregations$simple$cell_status$reason, "empty_restricted_set")
  expect_identical(rn$aggregations$dynamic$overall$status, "unavailable")
  # Executed perturbation: +1/4 on unit 3's date-4 increment moves m14 and the contrast by 1/8.
  rp <- run_chained_ipw_pipeline(d5_early_pair_panel(1 / 4), "id", "t", "y", "g", ctrl, target_cohorts = 3, control_group = "nevertreated",
    weighting = "optimal", base_period = "universal", run_numerical_transfer = FALSE)
  expect_equal(unname(rp$covariance_stream$mhat), c(1, 4 + 1 / 8), tolerance = 1e-12)
  expect_equal(unname(rp$basis$estimates), 25 / 8, tolerance = 1e-12)
  ref <- d5_locate_optional("verification/independent/d5_early_pair_reference.json")
  if (!is.null(ref)) {
    rr <- jsonlite::fromJSON(ref)
    expect_equal(unname(ru$basis$if_full[, 1L]), d5_rational(rr$contrast_IF), tolerance = 1e-12)
    expect_equal(unname(ru$basis$vcov[1L, 1L]), d5_rational(rr$contrast_variance), tolerance = 1e-12)
  }
})

test_that("isolated requests, rowless cohorts, permutations and explicit na.rm through the private pipeline", {
  ctrl <- make_d5_pipeline_control()
  df80 <- dget(d5_locate_fixture("rotating80_fixture.R"))$panel
  df80b <- df80[!(df80$g == 3 & df80$t == 5), ]
  base <- run_chained_ipw_pipeline(df80b, "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    requested_post = data.frame(g = 3, t = 5), target_post_times = data.frame(g = 3, t = c(3, 4, 5)),
    control_group = "nevertreated", base_period = "universal", weighting = "optimal", run_numerical_transfer = TRUE, root_witness = d5_root_witness)
  led <- base$support$ledger
  expect_identical(led$cell_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)", "ATT(g=3,t=5)"))
  expect_identical(led$support_status, c("estimable", "estimable", "isolated_unsupported"))
  expect_identical(led$requested_by, c("target_post_times", "target_post_times", "requested_post+target_post_times"))
  expect_identical(led$catalog_request_reason[3L], "no_target_pair_occupancy")
  expect_equal(led$n_catalog_rows_touching[3L], 0L)
  expect_true(led$has_clean_control[3L])
  post <- which(base$reporting$target_meta$role == "post")
  expect_equal(unname(base$reporting$estimates[post]), c(1, 1, NA), tolerance = 1e-12)
  expect_true(base$numerical_transfer$gmm_transfer$supported_subvector_certified)
  expect_false(base$numerical_transfer$gmm_transfer$whole_request_complete)
  expect_identical(base$support$unavailable$cell_id, "ATT(g=3,t=5)")
  expect_true(is.na(base$aggregations$simple$estimates[["ATT_simple"]]))
  expect_equal(unname(base$aggregations$dynamic$estimates), c(0, 0, 1, 1, NA), tolerance = 1e-12)
  # Permuting the declared request only reorders the input; the canonical ledger and every output coincide.
  perm <- run_chained_ipw_pipeline(df80b, "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    requested_post = data.frame(g = 3, t = 5), target_post_times = data.frame(g = 3, t = c(5, 3, 4)),
    control_group = "nevertreated", base_period = "universal", weighting = "optimal", run_numerical_transfer = TRUE, root_witness = d5_root_witness)
  expect_identical(perm$support$ledger, led)
  expect_equal(perm$reporting$estimates, base$reporting$estimates, tolerance = 0)
  expect_equal(perm$reporting$if_full, base$reporting$if_full, tolerance = 0)
  expect_equal(perm$gmm_fit$targets, base$gmm_fit$targets, tolerance = 0)
  # Explicit private na.rm = TRUE: restricted targets are displayed, never silently substituted.
  rn <- run_chained_ipw_pipeline(df80b, "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    target_post_times = data.frame(g = 3, t = c(3, 4, 5)), control_group = "nevertreated", base_period = "universal",
    weighting = "optimal", run_numerical_transfer = TRUE, root_witness = d5_root_witness, na.rm = TRUE)
  expect_true(rn$controls$na_rm); expect_true(rn$aggregations$simple$na_rm)
  expect_equal(rn$aggregations$simple$estimates[["ATT_simple"]], 1, tolerance = 1e-12)
  expect_true(rn$aggregations$simple$cell_status$restricted); expect_identical(rn$aggregations$simple$cell_status$excluded_ids, "ATT(g=3,t=5)")
  expect_identical(rn$aggregations$dynamic$overall$excluded_bins, "event_2"); expect_equal(rn$aggregations$dynamic$overall$estimate, 1, tolerance = 1e-12)
  expect_true(rn$aggregations$simple$numerical_certificate$certified)
  expect_false(rn$numerical_transfer$gmm_transfer$whole_request_complete)
  expect_true(rn$numerical_transfer$gmm_transfer$supported_subvector_certified)
  expect_false(rn$numerical_transfer$gmm_transfer$complete_inference_certified)
  # Genuine refusals persist alongside an unsupported cell: singular covariance, accuracy and resource budgets.
  df_zero <- df80b; df_zero$y_const <- ifelse(df_zero$g == 3 & df_zero$t >= 3, 1.0, 0.0)
  expect_error(run_chained_ipw_pipeline(df_zero, "id", "t", "y_const", "g", ctrl, xformla = ~1, target_cohorts = 3,
    target_post_times = data.frame(g = 3, t = c(3, 4, 5)), control_group = "nevertreated", base_period = "universal",
    weighting = "optimal", run_numerical_transfer = FALSE), class = "cdid_gmm_pd_error")
  expect_error(run_chained_ipw_pipeline(df80b, "id", "t", "y_pairperturbed", "g", make_d5_pipeline_control(max_target_error_abs = 1e-18),
    xformla = ~x, weightsname = "weight", target_cohorts = 3, target_post_times = data.frame(g = 3, t = c(3, 4, 5)),
    control_group = "nevertreated", base_period = "universal", weighting = "optimal", run_numerical_transfer = TRUE,
    root_witness = d5_root_witness), class = "cdid_gmm_accuracy_error")
  # Invalid requests remain input errors.
  expect_error(run_chained_ipw_pipeline(df80b, "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    target_post_times = data.frame(g = 3, t = 7), control_group = "nevertreated", run_numerical_transfer = FALSE),
    class = "cdid_ipw_pipeline_target_error")
  expect_error(run_chained_ipw_pipeline(df80b, "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    target_post_times = data.frame(g = 3, t = 2), control_group = "nevertreated", run_numerical_transfer = FALSE),
    class = "cdid_ipw_pipeline_target_error")
  # Rowless requested cohort: two cohort-4 singletons exist in the universe but no catalog row; its cells stay in the ledger.
  extra <- df80b[df80b$t == 1, ][1:2, ]; extra$id <- c(90001L, 90002L); extra$g <- 4
  d4 <- rbind(df80b, extra)
  r4 <- run_chained_ipw_pipeline(d4, "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = c(3, 4),
    control_group = "nevertreated", base_period = "varying", weighting = "optimal", run_numerical_transfer = FALSE)
  led4 <- r4$support$ledger
  expect_identical(led4$cell_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)", "ATT(g=3,t=5)", "ATT(g=4,t=4)", "ATT(g=4,t=5)"))
  expect_identical(led4$support_status[4:5], c("rowless_cohort", "rowless_cohort"))
  expect_false(any(led4$cohort_has_rows[4:5])); expect_true(all(is.na(led4$design_param_index[4:5])))
  expect_equal(r4$metadata$catalog$targets$n_rows[r4$metadata$catalog$targets$g == 4], 0L)
  expect_false(4 %in% r4$design_all$links$g)
  tm4 <- r4$reporting$target_meta
  expect_true(all(tm4$support_status[tm4$g == 4 & tm4$role == "post"] == "rowless_cohort"))
  expect_true(all(tm4$support_status[tm4$g == 4 & tm4$role == "pre"] == "unavailable_pretrend"))
  expect_true(all(is.na(r4$reporting$estimates[tm4$g == 4])))
  expect_equal(unname(r4$reporting$estimates[tm4$g == 3 & tm4$role == "post"]), c(1, 1, NA), tolerance = 1e-12)
  expect_equal(unname(r4$aggregations$group$estimates), c(NA_real_, NA_real_), tolerance = 0)
  expect_match(r4$aggregations$group$cell_status$reason[2L], "rowless_cohort")
  expect_equal(unname(r4$aggregations$calendar$estimates[r4$aggregations$calendar$cell_keys == 3]), 1, tolerance = 1e-12)
  # The reporting ledger holds the rowless cohort's template even without a GMM parameter.
  expect_true(any(tm4$g == 4 & tm4$role == "pre"))
})

test_that("an entirely empty catalog and insufficient pair-observing units return typed support-only results with no accepted GMM, quantile or certificate", {
  ctrl <- make_d5_pipeline_control()
  df80 <- dget(d5_locate_fixture("rotating80_fixture.R"))$panel
  # Every cohort-3 unit observed once: never-treated units still form pairs, but no target pair has occupancy.
  keep <- !(df80$g == 3 & df80$t != 1)
  r0 <- run_chained_ipw_pipeline(df80[keep, ], "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", weighting = "optimal", run_numerical_transfer = TRUE, root_witness = d5_root_witness)
  expect_s3_class(r0, "cdid_ipw_pipeline_result")
  expect_identical(r0$support_status, "empty_catalog")
  expect_null(r0$gmm_fit); expect_null(r0$numerical_transfer); expect_null(r0$basis)
  expect_equal(nrow(r0$metadata$catalog$rows), 0L)
  led <- r0$support$ledger
  expect_identical(led$cell_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)", "ATT(g=3,t=5)"))
  expect_true(all(led$support_status == "empty_catalog")); expect_false(any(led$available))
  tm0 <- r0$reporting$target_meta
  expect_true(all(is.na(r0$reporting$estimates[tm0$role != "reference"])))
  # Only the universal deterministic reference zero is a supported (trivially certified) cell; no ATT cell is certified.
  expect_identical(r0$reporting$numerical_certificate$certified_target_ids, tm0$target_id[tm0$role == "reference"])
  expect_false(any(tm0$target_id[tm0$role == "post"] %in% r0$reporting$numerical_certificate$certified_target_ids))
  expect_false(r0$reporting$numerical_certificate$whole_request_complete)
  expect_true(all(r0$auxiliary_reporting$descriptor_status$status == "unsupported_pair"))
  expect_true(all(vapply(r0$aggregations, function(a) identical(a$overall$status, "unavailable"), logical(1L))))
  for (nm in c("simple", "group", "calendar")) expect_identical(r0$aggregations[[nm]]$numerical_certificate$status, "no_supported_targets")
  expect_identical(r0$aggregations$dynamic$numerical_certificate$certified_target_ids, "event_-1")
  expect_equal(r0$resource_accounting$d5_arrays$n_ledger_cells, 3L)
  rv0 <- run_chained_ipw_pipeline(df80[keep, ], "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    control_group = "nevertreated", base_period = "varying", weighting = "optimal", run_numerical_transfer = FALSE)
  expect_identical(rv0$support_status, "empty_catalog")
  expect_identical(rv0$reporting$numerical_certificate$status, "no_supported_targets")
  expect_false(isTRUE(rv0$reporting$numerical_certificate$certified)); expect_true(is.na(rv0$reporting$cband_crit_val))
  # One pair-observing unit only (calendar 1..5 retained through explicit-NA rows).
  d1 <- df80; d1$y_exactfit <- NA_real_
  first_id <- df80$id[1L]; d1$y_exactfit[d1$id == first_id] <- df80$y_exactfit[df80$id == first_id]
  r1 <- run_chained_ipw_pipeline(d1, "id", "t", "y_exactfit", "g", ctrl, xformla = ~x, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", weighting = "optimal", run_numerical_transfer = FALSE)
  expect_identical(r1$support_status, "insufficient_pair_observing_units")
  expect_equal(r1$panel_owner$n_a_units, 1L); expect_null(r1$gmm_fit); expect_null(r1$metadata)
  expect_true(all(r1$support$ledger$support_status == "insufficient_pair_observing_units"))
  expect_equal(length(r1$aggregations), 0L)
})

test_that("supported-subvector inference scatters NA explicitly and an empty selection is typed", {
  ctrl <- make_d5_pipeline_control()
  X <- cbind(c(-2, -1, 1, 2, 0), NA_real_, c(4, 2, -2, -4, 0)); th <- c(a = 1, b = NA_real_, c = 2)
  z <- ipw_pipeline_supported_inference(th, X, c(TRUE, FALSE, TRUE), NULL, multiplier_budget = ctrl$multiplier_budget,
    inference_accuracy = ctrl$inference_accuracy)
  expect_identical(z$status, "supported_subvector")
  expect_equal(unname(z$estimates), c(1, NA, 2)); expect_true(is.na(z$se[["b"]]))
  expect_true(all(is.na(z$vcov[2L, ])) && all(is.na(z$vcov[, 2L])))
  expect_equal(unname(z$vcov[c(1L, 3L), c(1L, 3L)]), crossprod(X[, c(1L, 3L)]) / 25, tolerance = 1e-14)
  expect_identical(z$support$unavailable_ids, "b"); expect_false(z$support$whole_request_complete)
  e <- ipw_pipeline_supported_inference(th, X, c(FALSE, FALSE, FALSE), NULL, multiplier_budget = ctrl$multiplier_budget,
    inference_accuracy = ctrl$inference_accuracy)
  expect_identical(e$status, "empty_supported_selection"); expect_null(e$multiplier); expect_true(is.na(e$cband_crit_val))
  expect_identical(e$numerical_certificate$band_family, "none"); expect_false(isTRUE(e$numerical_certificate$certified))
  expect_error(ipw_pipeline_supported_inference(th, X, c(TRUE, TRUE, TRUE), NULL, multiplier_budget = ctrl$multiplier_budget,
    inference_accuracy = ctrl$inference_accuracy), class = "cdid_ipw_pipeline_integrity_error")
})

test_that("simultaneous bands under partial support use only certified IDs with caller-fixed signs (deterministic)", {
  ctrl <- make_d5_pipeline_control()
  X2 <- rbind(c(-2, 4), c(-1, 2), c(1, -2), c(2, -4), c(0, 0)); X3 <- cbind(X2[, 1L], NA_real_, X2[, 2L])
  th2 <- c(a = 1, c = 2); th3 <- c(a = 1, b = NA_real_, c = 2); clusters <- c(1L, 1L, 2L, 2L, 3L)
  signs <- unname(as.matrix(expand.grid(rep(list(c(-1.0, 1.0)), 3L)))); storage.mode(signs) <- "double"
  enc2 <- list(theta = list(lo = matrix(th2, 1L), hi = matrix(th2, 1L)), rows = list(lo = X2, hi = X2))
  enc3 <- list(theta = list(lo = matrix(th3, 1L), hi = matrix(th3, 1L)), rows = list(lo = X3, hi = X3))
  direct <- compute_chained_ipw_inference(th2, X2, clusters, cband = TRUE, biters = 8L, multiplier_budget = ctrl$multiplier_budget,
    inference_accuracy = ctrl$inference_accuracy, fixed_signs = signs, enclosure = enc2, eligible_full = c(TRUE, TRUE, TRUE, TRUE, FALSE))
  part <- ipw_pipeline_supported_inference(th3, X3, c(TRUE, FALSE, TRUE), clusters, cband = TRUE, biters = 8L,
    multiplier_budget = ctrl$multiplier_budget, inference_accuracy = ctrl$inference_accuracy, enclosure = enc3,
    fixed_signs = signs, eligible_full = c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(part$cband_crit_val, direct$cband_crit_val, tolerance = 0)
  expect_equal(unname(part$se[c("a", "c")]), unname(direct$se), tolerance = 0); expect_true(is.na(part$se[["b"]]))
  expect_equal(unname(part$cband_lower[c("a", "c")]), unname(direct$cband_lower), tolerance = 0)
  expect_true(is.na(part$cband_lower[["b"]])); expect_true(is.na(part$cband_upper[["b"]]))
  expect_identical(part$numerical_certificate$band_family, "simultaneous_multiplier_sup_t_over_certified_ids")
  expect_identical(part$numerical_certificate$certified_target_ids, c("a", "c"))
  expect_identical(part$numerical_certificate$unavailable_target_ids, "b")
  expect_false(part$numerical_certificate$whole_request_complete); expect_true(part$numerical_certificate$certified)
  expect_true(all(is.na(part$numerical_certificate$band$lo[, 2L])))
  expect_equal(part$numerical_certificate$errors, direct$numerical_certificate$errors, tolerance = 0)
  expect_equal(part$n_zero_eligible_clusters, 1L)
})

test_that("rowless-only requests retain other cohorts' moments and never certify an empty reporting selection", {
  ctrl <- make_d5_pipeline_control()
  d <- d5_case6_panel(); extra <- d[d$t == 1, ][1:2, ]; extra$id <- 9:10; extra$g <- 3
  r <- run_chained_ipw_pipeline(rbind(d, extra), "id", "t", "y", "g", ctrl,
    target_cohorts = c(2, 3), target_post_times = data.frame(g = 3, t = 3),
    control_group = "nevertreated", base_period = "varying", root_witness = d5_root_witness)
  expect_identical(r$support$ledger$cell_id, "ATT(g=3,t=3)")
  expect_identical(r$support$ledger$support_status, "rowless_cohort")
  expect_equal(r$reduction_summary$n_links_all, 3L)
  expect_equal(r$reduction_summary$n_links_used, 2L)
  expect_equal(r$gmm_fit$targets$estimate, c(2, 5), tolerance = 1e-12)
  expect_true(all(is.na(r$reporting$estimates)))
  expect_identical(r$reporting$inference$status, "empty_supported_selection")
  expect_false(r$numerical_transfer$gmm_transfer$supported_subvector_certified)
  expect_false(r$numerical_transfer$gmm_transfer$complete_inference_certified)
  expect_length(r$numerical_transfer$gmm_transfer$certified_target_ids, 0L)
  expect_identical(r$numerical_transfer$gmm_transfer$band_family, "none")
  # A numerical covariance failure in those retained moments still refuses.
  d$y <- ifelse(d$g == 2 & d$t >= 2, 1, 0)
  expect_error(run_chained_ipw_pipeline(rbind(d, extra), "id", "t", "y", "g", ctrl,
    target_cohorts = c(2, 3), target_post_times = data.frame(g = 3, t = 3),
    control_group = "nevertreated", base_period = "varying", run_numerical_transfer = FALSE), class = "cdid_gmm_pd_error")
})

test_that("zero loading columns still pass through the full nonempty GMM and transfer checks", {
  ctrl <- make_d5_pipeline_control()
  d <- d5_case3_panel(); d$y[d$t == 1] <- 0; d$y[d$t == 2] <- NA_real_
  r <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 3,
    control_group = "nevertreated", base_period = "varying", root_witness = d5_root_witness)
  expect_equal(r$reduction_summary$n_links_used, 1L)
  expect_true(all(!r$gmm_fit$targets$estimable))
  expect_equal(ncol(r$gmm_fit$if_projection$materialized_target_if), 0L)
  expect_length(r$numerical_transfer$gmm_transfer$total_target_error_cert, 0L)
  expect_null(r$basis)
  expect_false(r$numerical_transfer$gmm_transfer$supported_subvector_certified)
  expect_identical(r$reporting$inference$status, "empty_supported_selection")
  # The same zero-loading solve is required if the only reported cell is rowless.
  extra <- d[d$t == 1, ][1:2, ]; extra$id <- 9:10; extra$g <- 4; extra$t <- 4
  r <- run_chained_ipw_pipeline(rbind(d, extra), "id", "t", "y", "g", ctrl,
    target_cohorts = c(3, 4), target_post_times = data.frame(g = 4, t = 4),
    control_group = "nevertreated", base_period = "varying", root_witness = d5_root_witness)
  expect_identical(r$support$ledger$support_status, "rowless_cohort")
  expect_identical(r$support$core_targets$kind, "internal_unavailable")
  expect_true(all(is.na(r$gmm_fit$targets$estimate)))
  expect_equal(ncol(r$gmm_fit$if_projection$materialized_target_if), 0L)
  expect_false(r$numerical_transfer$gmm_transfer$supported_subvector_certified)
})

test_that("basis certificates identify their own band family when reported levels have no supported selection", {
  r <- run_chained_ipw_pipeline(d5_early_pair_panel(), "id", "t", "y", "g", make_d5_pipeline_control(),
    target_cohorts = 3, control_group = "nevertreated", base_period = "varying", root_witness = d5_root_witness)
  gt <- r$numerical_transfer$gmm_transfer
  expect_identical(gt$band_family, r$basis$numerical_certificate$band_family)
  expect_identical(gt$band_families$post_basis$certified_target_ids, names(r$basis$estimates))
  expect_null(gt$band_families$reported_cells)
  # Universal adds a distinct exact reference family; it is not the basis band.
  u <- run_chained_ipw_pipeline(d5_early_pair_panel(), "id", "t", "y", "g", make_d5_pipeline_control(),
    target_cohorts = 3, control_group = "nevertreated", base_period = "universal", root_witness = d5_root_witness)
  gt <- u$numerical_transfer$gmm_transfer
  expect_identical(gt$band_family, "separate_families_see_band_families")
  expect_identical(gt$band_families$reported_cells$band_family, "deterministic_exact_zero_over_certified_ids")
  expect_identical(gt$band_families$post_basis$band_family, "pointwise_normal_over_certified_ids")
  expect_setequal(unlist(lapply(gt$band_families, `[[`, "certified_target_ids")), gt$certified_target_ids)
})

test_that("restricted ratio metadata keeps multiplicity and refuses known zero share mass", {
  ctrl <- make_d5_pipeline_control(); z <- d5_case4_inputs()
  agg <- function(type, g = z$g) aggregate_chained_ipw_targets(z$th, z$X, z$meta, g, z$w, type = type,
    na.rm = TRUE, multiplier_budget = ctrl$multiplier_budget, inference_accuracy = ctrl$inference_accuracy)
  expect_equal(agg("simple")$cell_status$share_denominator, 3 / 4, tolerance = 0)
  expect_equal(agg("group")$overall$share_denominator, 3 / 4, tolerance = 0)
  a <- agg("calendar", rep(Inf, 4))
  expect_true(all(is.na(a$estimates)))
  expect_true(all(a$cell_status$reason == "nonpositive_share_mass"))
  expect_false(isTRUE(a$numerical_certificate$certified))
  a <- agg("group", rep(Inf, 4))
  expect_equal(unname(a$estimates), c(2, 6), tolerance = 0)
  expect_identical(a$overall$reason, "nonpositive_share_mass")
  # Restore B as a supported conditional input: two g=2 cells count twice in S.
  z$th[2L] <- 4; z$X[, 2L] <- c(0, 0, 1, -1)
  a <- agg("simple")
  expect_equal(a$estimates[[1L]], 18 / 5, tolerance = 1e-12)
  expect_equal(a$cell_status$share_denominator, 5 / 4, tolerance = 0)
  # Independent Fraction differentiation of (p2*2+p2*4+p3*6)/(2*p2+p3)
  # gives (-24/25,48/25); audit/independent_results.json records the derivation.
  expect_equal(unname(a$gradient_shares[1L, ]), c(-24 / 25, 48 / 25), tolerance = 1e-12)
})

test_that("empty-catalog outputs reserve downstream work and exact references need no quantile", {
  ctrl <- make_d5_pipeline_control(); d <- d5_case3_panel(); d <- d[!(d$g == 3 & d$t != 2), ]
  r <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", run_numerical_transfer = FALSE)
  expect_identical(r$support_status, "empty_catalog")
  expect_null(r$gmm_fit)
  expect_equal(r$reporting$cband_crit_val, 0, tolerance = 0)
  expect_identical(r$reporting$numerical_certificate, NULL)
  expect_gt(r$resource_accounting$downstream_work_ops, 0)
  ctrl$owner_control$max_work_ops <- as.double(r$resource_accounting$downstream_work_ops - 1)
  expect_error(run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", run_numerical_transfer = FALSE),
    "empty-catalog reporting/aggregation/interval/multiplier whole-path work envelope")
})

test_that("support-count readers propagate resource failures after a genuine pair-support absence", {
  env <- environment(build_ipw_auxiliary_reporting)
  original <- get("construct_ipw_nuisance_view", env)
  original_raw <- get("ipw_raw_pair_support_for_cohorts", env)
  replacement <- function(...) {
    args <- list(...)
    if (!is.null(args$nuisance_id) && startsWith(args$nuisance_id, "reporting_probe:"))
      cdid_ipw_pipeline_abort("audit injected support-count reader bound refusal", "cdid_ipw_owner_bound_error")
    do.call(original, args)
  }
  replacement_raw <- function(...) {
    args <- list(...)
    if (isFALSE(args$require_support))
      cdid_ipw_pipeline_abort("audit injected support-count reader bound refusal", "cdid_ipw_owner_bound_error")
    do.call(original_raw, args)
  }
  d5_replace_binding("construct_ipw_nuisance_view", replacement, env)
  d5_replace_binding("ipw_raw_pair_support_for_cohorts", replacement_raw, env)
  on.exit({
    d5_replace_binding("construct_ipw_nuisance_view", original, env)
    d5_replace_binding("ipw_raw_pair_support_for_cohorts", original_raw, env)
  }, add = TRUE)
  expect_error(run_chained_ipw_pipeline(d5_case3_panel(), "id", "t", "y", "g", make_d5_pipeline_control(),
    target_cohorts = 3, control_group = "nevertreated", base_period = "universal", run_numerical_transfer = FALSE),
    class = "cdid_ipw_owner_bound_error")
  expect_error(run_chained_ipw_pipeline(d5_case3_panel(), "id", "t", "y", "g", make_d5_pipeline_control(),
    target_cohorts = 3, control_group = "nevertreated", base_period = "universal", run_numerical_transfer = TRUE,
    root_witness = d5_root_witness),
    class = "cdid_ipw_owner_bound_error")
})

test_that("supported direct comparisons propagate normalization failures and full-frame IDs stay aligned", {
  ctrl <- make_d5_pipeline_control(); d <- d5_case6_panel(); d$g[d$g == 2] <- 3
  local({
    env <- environment(build_ipw_auxiliary_reporting); original <- get("construct_ipw_link_view", env)
    replacement <- function(...) {
      args <- list(...)
      if (!is.null(args$link_id) && startsWith(args$link_id, "reporting:"))
        cdid_ipw_pipeline_abort("audit injected supported auxiliary normalization refusal", "cdid_ipw_normalization_error")
      do.call(original, args)
    }
    d5_replace_binding("construct_ipw_link_view", replacement, env)
    on.exit(d5_replace_binding("construct_ipw_link_view", original, env), add = TRUE)
    expect_error(run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 3,
      control_group = "nevertreated", base_period = "universal", run_numerical_transfer = FALSE), class = "cdid_ipw_normalization_error")
  })
  d <- d5_case6_panel(); fids <- c(9L, 8:1)
  cl <- c("frame_only", "a", "a", "b", "b", "c", "c", "d", "d")
  r <- run_chained_ipw_pipeline(d, "id", "t", "y", "g", ctrl, target_cohorts = 2,
    frame_id = fids, frame_Z = matrix(1, 9, 1), frame_weights = c(3, rep(1, 8)),
    allow_frame_only = TRUE, cluster_vec = cl, control_group = "nevertreated", base_period = "universal",
    root_witness = d5_root_witness)
  post <- which(r$reporting$target_meta$role == "post")
  f12 <- c(-2,-2,2,2,2,2,-2,-2); f23 <- c(-2,2,-2,2,2,-2,2,-2)
  expected <- rbind(c(0, 0), cbind(f12, f12 + f23)[8:1, ]) * 9 / 8
  expect_equal(unname(r$reporting$if_full[, post]), unname(expected), tolerance = 1e-12)
  expect_equal(r$reporting$inference$n_zero_eligible_clusters, 1L)
  expect_false(r$aggregations$simple$shares_known)
  expect_equal(r$aggregations$simple$estimates[[1L]], 7 / 2, tolerance = 1e-12)
  sums <- rowsum(expected, cl, reorder = FALSE); sums <- sweep(sums, 2L, colMeans(sums), "-")
  expect_equal(unname(r$reporting$vcov[post, post]), unname(crossprod(sums) / 81), tolerance = 1e-12)
  expect_true(r$numerical_transfer$gmm_transfer$complete_inference_certified)
})
