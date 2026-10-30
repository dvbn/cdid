# Automatic fourteen-control production for the public CP02 candidate
# (R/chained_public_control.R). Deterministic helper tests: hand-computed
# dimensions and bounds, the exact adopted defaults (D27 memory/worker defaults;
# D78 optional work cap, disabled by default), the accepted automatic witness,
# overflow-safe refusals and a finite-cap refusal of the accepted pipeline
# before any nuisance fit. No DGP draw, RNG, fit or bootstrap. Runs in the
# installed package or in the source harness with the copied native DLL (did is
# not needed here). The former D27 finite work default (1e11) is historical
# evidence only; D78 expectations live here and in test-cp08-optional-work-cap.R.

ac_small_panel <- function() {
  d <- expand.grid(id = 1:6, period = 1:4)
  d$G <- c(0, 0, 3, 3, 4, 4)[d$id]; d$x <- c(0.5, -1, 1, 0, 2, -2)[d$id]; d$y <- d$id + d$period
  d
}

test_that("the NULL control resolves to the exact adopted defaults (D78: no user work budget; D27 memory/worker), recorded as adopted", {
  r <- chained_public_resolve_control(NULL)
  expect_identical(r$profile, list(numeric_workspace_bytes = 1073741824, counted_work_limit = Inf, workers = 1))
  expect_identical(r$profile_status, "ADOPTED_D78_OPTIONAL_WORK_CAP_D27_MEMORY_DEFAULTS")
  expect_identical(r$overrides, character(0L)); expect_false(r$relative_accuracy_override_recorded)
  expect_identical(r$relative_accuracy_source, "deferred_by_D23")
  expect_false("target_cohorts" %in% names(r)); expect_null(r$chunk_rows)
  # carried development settings (not tolerance searches): nuisance solver, certificate precision, GMM guard predicates
  expect_identical(r$glm, list(epsilon = 1e-12, maxit = 50L))
  expect_identical(r$certificate, list(taylor_degree = 12L, max_abs_eta = 50.0, max_squarings = 16L))
  expect_identical(names(r$gmm_guard), c("max_orthogonality_defect", "min_certified_lambda_rel", "max_factor_residual_rel", "min_rcond_R",
                                        "min_rcond_U", "min_sigma_rel_A", "max_kkt_residual_rel", "max_loading_residual_rel"))
  # Historical pilot evidence keeps its old finite D27 work profile; under D78 it
  # is an explicit user cap (an override), not the current default.
  explicit <- chained_public_resolve_control(list(numeric_workspace_bytes = 1073741824,
    counted_work_limit = 100000000000, workers = 1))
  expect_identical(explicit$profile$counted_work_limit, 100000000000)
  expect_identical(explicit$profile[c("numeric_workspace_bytes", "workers")], r$profile[c("numeric_workspace_bytes", "workers")])
  expect_setequal(explicit$overrides, c("numeric_workspace_bytes", "counted_work_limit", "workers"))
  expect_identical(explicit$covariance_strategy, r$covariance_strategy)
})

test_that("explicit advanced overrides are validated and recorded; unknown, parallel or invalid entries refuse", {
  r <- chained_public_resolve_control(list(counted_work_limit = 5e9, numeric_workspace_bytes = 2^30,
                                           chunk_rows = 7, workers = 1))
  # D23 defers both the public guarantee and its override.
  expect_null(r$profile$relative_accuracy); expect_false(r$relative_accuracy_override_recorded)
  expect_identical(r$relative_accuracy_source, "deferred_by_D23")
  expect_error(chained_public_resolve_control(list(relative_accuracy = 1e-5)), class = "cdid_chained_public_accuracy_deferred")
  expect_identical(r$profile$counted_work_limit, 5e9); expect_identical(r$profile$numeric_workspace_bytes, 2^30)
  expect_identical(r$chunk_rows, 7); expect_false("target_cohorts" %in% names(r))
  expect_identical(r$glm, list(epsilon = 1e-12, maxit = 50L)); expect_identical(r$certificate$taylor_degree, 12L)
  expect_identical(r$gmm_guard$min_rcond_R, 1e-10); expect_identical(r$gmm_guard$min_rcond_U, 1e-10)
  expect_setequal(r$overrides, c("counted_work_limit", "numeric_workspace_bytes", "chunk_rows", "workers"))
  e0 <- chained_public_resolve_control(list())
  expect_identical(e0$profile, chained_public_profile_defaults()); expect_false(e0$relative_accuracy_override_recorded)
  refuse <- function(x) expect_error(chained_public_resolve_control(x), class = "cdid_chained_public_control_error")
  # CP08: `workers` declares multiplier-bootstrap worker processes; the default stays 1.
  expect_identical(chained_public_resolve_control(list(workers = 2))$profile$workers, 2)
  refuse(list(bogus = 1)); refuse(list(workers = 0)); refuse(list(workers = 1.5)); refuse(list(relative_accuracy = 0)); refuse(list(relative_accuracy = 1))
  refuse(list(relative_accuracy = c(1e-6, 1e-6))); refuse(list(numeric_workspace_bytes = 4)); refuse(list(counted_work_limit = 0))
  refuse(list(counted_work_limit = 2^53)); refuse(list(chunk_rows = 0.5)); refuse(list(target_cohorts = c(3, 3))); refuse(list(target_cohorts = NA))
  refuse(list(target_cohorts = 3)); refuse(list(glm = list(epsilon = 1e-12))); refuse(list(glm = list(epsilon = -1))); refuse(list(certificate = list(taylor_degree = 1))); refuse(list(gmm_guard = list(max_orthogonality_defect = 0.5)))
  refuse(list(gmm_guard = list(other = 0.1))); refuse(list(1, 2)); refuse("prospective")
  refuse(structure(list(), class = "cdid_ipw_pipeline_control"))      # the private constructor output is not the public control
  refuse(data.frame(relative_accuracy = 1e-6))
})

test_that("dimensions and finite bounds are derived from the one prepared input and match hand arithmetic", {
  prep <- chained_request_prepare_inputs(ac_small_panel(), NULL, "y", "period", "id", "G", ~x)
  dims <- chained_public_dimensions(prep, biters = 50)
  # hand: n = 6 frame units, 24 rows, T = 4, cohorts {3, 4, Inf}, targets {3, 4}, p = 2 (intercept + x)
  expect_identical(dims$n_frame, 6); expect_identical(dims$n_rows, 24); expect_identical(dims$n_dates, 4)
  expect_identical(dims$n_cohorts, 3); expect_identical(dims$n_targets, 2); expect_identical(dims$n_covariates, 2)
  expect_identical(dims$biters, 50); expect_identical(dims$n_clusters, 6)
  expect_identical(dims$pairs_per_target, 6)            # 4 * 3 / 2
  expect_identical(dims$R_full, 3)                      # dates >= 3: {3, 4}; dates >= 4: {4}
  expect_identical(dims$candidates_bound, 60)           # 2 * 6 * (1 + 4)
  expect_identical(dims$links_bound, 36)                # 2 * 6 * 3
  expect_identical(dims$Q_bound, 16)                    # 2 * 2 * 4
  expect_identical(dims$params_bound, 74)               # 2 * 36 + 2
  expect_identical(dims$J_bound, 98)                    # 74 + 16 + 3 + 3 + 2
  expect_identical(dims$mask_entries_bound, 216)        # 36 * 6
  expect_identical(dims$pool_entries_bound, 108)        # 36 * 3
  expect_identical(dims$E_entries_bound, 144)           # 36 * 4
  expect_identical(dims$materialized_if_bound, 588)     # 6 * 98
  # the counts-based entry gives the same numbers from the same counts
  d2 <- chained_public_dimensions_from_counts(6, 6, 24, 4, 3, 2, 2, 50, 6, calendar = 1:4, target_cohorts = c(3, 4))
  expect_identical(d2[setdiff(names(d2), "scope")], dims[setdiff(names(dims), "scope")])
  # declared clusters change n_clusters only
  dc <- ac_small_panel(); dc$cl <- c(1, 1, 2, 2, 3, 3)[dc$id]
  pc <- chained_request_prepare_inputs(dc, NULL, "y", "period", "id", "G", ~x, NULL, "cl")
  expect_identical(chained_public_dimensions(pc, 50, clustered = TRUE)$n_clusters, 3)
})

test_that("all fourteen controls are produced from dimensions and the profile with the expected values", {
  prep <- chained_request_prepare_inputs(ac_small_panel(), NULL, "y", "period", "id", "G", ~x)
  dims <- chained_public_dimensions(prep, biters = 50)
  out <- chained_public_pipeline_control(dims, chained_public_resolve_control(NULL))
  ctrl <- out$control; rec <- out$record
  expect_s3_class(ctrl, "cdid_ipw_pipeline_control")
  expect_setequal(names(ctrl), c("owner_control", "fitted_control", "point_control", "glm_control", "certificate_control", "problem_control",
                                 "link_control", "covariance_control", "class_bounds", "gmm_bounds", "gmm_accuracy", "catalog_bounds",
                                 "multiplier_budget", "inference_accuracy"))
  entries <- 1073741824 / 8
  # D78: the public value stays Inf (no user work budget); the finite private
  # constructors receive the exact-count arithmetic ceiling 2^53 - 1, recorded distinctly.
  ceiling_work <- 2^53 - 1
  expect_identical(rec$numeric_workspace_entries, entries); expect_identical(rec$counted_work_limit, Inf)
  expect_identical(rec$work_cap$status, "disabled"); expect_identical(rec$work_cap$effective_private_max_work_ops, ceiling_work)
  expect_identical(rec$chunk_rows, 6)                                      # min(n_frame, 256, workspace rule)
  expect_identical(ctrl$owner_control$max_chunk_rows, 6); expect_identical(ctrl$owner_control$max_matrix_entries, entries)
  expect_identical(ctrl$owner_control$max_retained_entries, entries); expect_identical(ctrl$owner_control$max_output_entries, entries)
  expect_identical(ctrl$owner_control$max_work_ops, ceiling_work)
  expect_identical(ctrl$fitted_control$max_observer_units, 6); expect_identical(ctrl$point_control$max_observer_units, 6)
  expect_identical(ctrl$glm_control, list(epsilon = 1e-12, maxit = 50L, trace = FALSE))
  expect_identical(ctrl$certificate_control$taylor_degree, 12); expect_identical(ctrl$certificate_control$max_abs_eta, 50)
  expect_identical(ctrl$certificate_control$max_squarings, 16); expect_identical(ctrl$certificate_control$max_work_ops, ceiling_work)
  expect_identical(ctrl$problem_control$max_matrix_entries, entries); expect_identical(ctrl$link_control$max_output_entries, entries)
  expect_identical(ctrl$covariance_control$max_cluster_rows, 6); expect_identical(ctrl$covariance_control$max_operator_block_rows, 6)
  expect_identical(ctrl$class_bounds, list(max_links = 36, max_classes = 36, max_class_dates = 4, max_E_entries = 144, max_mask_entries = 216))
  expect_identical(ctrl$gmm_bounds$max_targets, 98); expect_identical(ctrl$gmm_bounds$max_batch_rows, 6)
  expect_identical(ctrl$gmm_bounds$max_open_clusters, 6); expect_identical(ctrl$gmm_bounds$max_materialized_if_entries, 588)
  expect_true(ctrl$gmm_bounds$allow_materialize_if)
  expect_identical(ctrl$catalog_bounds, list(max_candidates = 60, max_rows = 36, max_pool_entries = 108, max_scratch = entries, max_work = ceiling_work))
  mb <- ctrl$multiplier_budget
  expect_identical(mb$max_units, 6); expect_identical(mb$max_clusters, 6); expect_identical(mb$max_targets, 98); expect_identical(mb$max_draws, 50)
  expect_identical(mb$max_if_batch_rows, 6); expect_identical(mb$max_if_batch_cells, 588); expect_identical(mb$max_sign_batch_rows, 6)
  expect_identical(mb$max_sign_batch_cells, 36); expect_identical(mb$max_workspace_cells, entries); expect_identical(mb$max_workspace_bytes, 1073741824)
  # D23: conventional method; no per-fit interval/accuracy acceptance.
  expect_null(ctrl$inference_accuracy)
  expect_identical(ctrl$gmm_accuracy$method, "conventional_chol_qr")
  expect_identical(ctrl$gmm_accuracy$target_error_unit, "abs"); expect_identical(ctrl$gmm_accuracy$max_target_error_abs, .Machine$double.xmax)
  expect_identical(ctrl$gmm_accuracy$max_loading_error_abs, .Machine$double.xmax); expect_true(is.na(ctrl$gmm_accuracy$max_target_error_se))
  expect_true(is.finite(chained_public_staging_envelope()))
  expect_identical(ctrl$gmm_accuracy$min_certified_lambda_rel, 1e-10); expect_true(ctrl$gmm_accuracy$require_matching_scale)
  # the record states what is and is not claimed
  expect_identical(rec$profile_status, "ADOPTED_D78_OPTIONAL_WORK_CAP_D27_MEMORY_DEFAULTS")
  expect_false(rec$admission$whole_process_RSS_guarantee); expect_match(rec$admission$scope, "not process RSS")
  expect_identical(rec$admission$counted_work_limit, Inf); expect_identical(rec$admission$work_cap_status, "disabled")
  expect_identical(rec$numerical_policy$decision, "D23"); expect_false(rec$numerical_policy$override_recorded)
  expect_null(rec$gmm_stage_absolute_budgets); expect_match(rec$witness, "not requested/performed")
  expect_true(is.data.frame(rec$derived_bounds) && nrow(rec$derived_bounds) == 16L)
  # explicit overrides reach the produced controls and the record
  r2 <- chained_public_resolve_control(list(chunk_rows = 4))
  c2 <- chained_public_pipeline_control(dims, r2)
  expect_null(c2$control$inference_accuracy); expect_identical(c2$record$chunk_rows, 4)
  expect_identical(c2$control$owner_control$max_chunk_rows, 4); expect_identical(c2$control$glm_control$maxit, 50L)
  expect_identical(c2$control$gmm_accuracy$min_rcond_R, 1e-10); expect_false(c2$record$numerical_policy$override_recorded)
  expect_identical(c2$record$chunk_rule, "explicit control$chunk_rows")
})

test_that("chunk rows follow the workspace/link rule, the cap and the frame; multiplier budgets follow the draw count", {
  big <- chained_public_dimensions_from_counts(n_frame = 1e5, n_panel_units = 1e5, n_rows = 1e6, n_dates = 10, n_cohorts = 4, n_targets = 3,
                                               n_covariates = 2, biters = 1000, n_clusters = 1e5)
  expect_identical(big$pairs_per_target, 45); expect_identical(big$links_bound, 540)
  r <- chained_public_resolve_control(NULL)
  # floor(67108864 / (2048 * (540 + 2 + 1))) = 60 < 256
  expect_identical(chained_public_chunk_rows(big, 67108864, r), 60)
  small_links <- chained_public_dimensions_from_counts(1e5, 1e5, 2e5, 2, 2, 1, 1, 1000, 1e5)
  expect_identical(chained_public_chunk_rows(small_links, 67108864, r), 256)             # cap
  tiny <- chained_public_dimensions_from_counts(5, 5, 10, 2, 2, 1, 1, 1000, 5)
  expect_identical(chained_public_chunk_rows(tiny, 67108864, r), 5)                       # frame
  r7 <- chained_public_resolve_control(list(chunk_rows = 7))
  expect_identical(chained_public_chunk_rows(tiny, 67108864, r7), 5)                      # explicit, capped by frame
  expect_identical(chained_public_chunk_rows(big, 67108864, r7), 7)
  mb <- chained_public_multiplier_budget(big, 2000, 60, 67108864, 1073741824)
  expect_identical(mb$max_draws, 2000); expect_identical(mb$max_sign_batch_rows, 60); expect_identical(mb$max_sign_batch_cells, 60 * 1e5)
  mb2 <- chained_public_multiplier_budget(big, 7, 60, 67108864, 1073741824)
  expect_identical(mb2$max_draws, 7); expect_identical(mb2$max_sign_batch_rows, 7)
})

test_that("counts outside exact integer arithmetic or inconsistent counts are refused, never wrapped", {
  expect_error(chained_public_dimensions_from_counts(2^52, 2^52, 2^52, 10, 4, 3, 2, 1000, 2^52), class = "cdid_ipw_owner_bound_error")
  expect_error(chained_public_dimensions_from_counts(10, 10, 20, 2^30, 4, 3, 2, 1000, 10), class = "cdid_ipw_owner_bound_error")
  expect_error(chained_public_dimensions_from_counts(10, 11, 20, 4, 3, 2, 2, 50, 10), class = "cdid_chained_public_dimension_error")   # panel units > frame
  expect_error(chained_public_dimensions_from_counts(10, 10, 20, 4, 3, 4, 2, 50, 10), class = "cdid_chained_public_dimension_error")   # targets > cohorts
  expect_error(chained_public_dimensions_from_counts(10, 10, 20, 4, 3, 2, 2, 50, 11), class = "cdid_chained_public_dimension_error")   # clusters > frame
  expect_error(chained_public_dimensions_from_counts(10.5, 10, 20, 4, 3, 2, 2, 50, 10), class = "cdid_chained_public_dimension_error")
  expect_error(chained_public_dimensions_from_counts(Inf, 10, 20, 4, 3, 2, 2, 50, 10), class = "cdid_chained_public_dimension_error")
})

test_that("the accepted one-radius automatic witness is used verbatim", {
  w <- chained_public_automatic_witness(list(factors = list(H_chol = matrix(c(2, 0, 1, 3), 2L, dimnames = list(c("a", "b"), c("a", "b"))))))
  expect_identical(sort(names(w)), c("F", "r"))
  expect_identical(w$F, matrix(c(2, 0, 1, 3), 2L))                   # names stripped, bytes unchanged
  expect_identical(w$r, logit_cert_automatic_radius_request())
  expect_true(logit_cert_is_automatic_radius_request(w$r))
  expect_false(logit_cert_is_automatic_radius_request(unclass(w$r)))
})

test_that("a finite workspace below the accepted envelope refuses at the owner preflight, before any nuisance fit", {
  prep <- chained_request_prepare_inputs(ac_small_panel(), NULL, "y", "period", "id", "G", ~x)
  dims <- chained_public_dimensions(prep, biters = 50)
  tiny <- chained_public_pipeline_control(dims, chained_public_resolve_control(list(numeric_workspace_bytes = 4096)))
  expect_identical(tiny$record$numeric_workspace_entries, 512)
  e <- tryCatch(run_chained_ipw_pipeline(prep$data_internal, "id", "period", "y", "G", tiny$control, xformla = ~x,
    target_cohorts = prep$target_cohorts, control_group = "nevertreated", run_numerical_transfer = FALSE,
    root_witness = chained_public_automatic_witness, aggregations = character(0L), frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
    frame_weights = prep$frame_w, frame_g = prep$frame_g_arg, bstrap = FALSE, caller_rng = TRUE), error = function(e) e)
  expect_s3_class(e, "cdid_ipw_owner_bound_error")
  expect_match(conditionMessage(e), "exceeds explicit bound")
  expect_identical(e$diagnostics$bound, 512)                           # the explicit finite cap, no relaxation
})


test_that("zero-loading transfer has exact empty inference objects without changing native admission", {
  E <- list(theta=list(lo=matrix(numeric(),1L,0L),hi=matrix(numeric(),1L,0L)),
    rows=list(lo=matrix(numeric(),4L,0L),hi=matrix(numeric(),4L,0L)))
  z <- ipw_pipeline_inference_intervals(E,c(1,1,2,2))
  expect_identical(z$theta,E$theta);expect_identical(z$rows,E$rows)
  expect_identical(dim(z$cluster$lo),c(2L,0L))
  expect_identical(z$vcov,list(lo=matrix(0,0L,0L),hi=matrix(0,0L,0L)))
  expect_identical(z$se,matrix(0,2L,0L))
  expect_error(ipw_pipeline_interval_product(t(E$rows$lo),t(E$rows$hi),E$rows$lo,E$rows$hi),class="cdid_gmm_bound_error")
})
