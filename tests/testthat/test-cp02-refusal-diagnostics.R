# Refusal diagnostic retention (October 17, 2026): the three refusal stages of
# the frozen public pilot (nuisance root certificate INCONCLUSIVE, D20
# scale-aware refusal, Cholesky refusal in the GMM core) keep their classes,
# messages and original diagnostic fields and additionally carry compact
# packets of already computed state. Deterministic unit fixtures only: one
# fixed 24-unit panel written by formula (no RNG, no DGP draw), direct core
# calls on a three-link cycle design, and pure-function checks of the helpers.
# These fixture fits are unit tests, not scientific fits; no did/DRDID is
# needed (the request-path MP construction, which needs did, is covered by the
# root-only host reproduction in tools/cp02_refusal_diagnostic_reproduction.R).

rd_no_refs <- function(x) {
  # TRUE when no environment, function or external pointer occurs anywhere.
  if (is.environment(x) || is.function(x) || typeof(x) == "externalptr") return(FALSE)
  values_ok <- !is.list(x) || all(vapply(x, rd_no_refs, logical(1L)))
  at <- attributes(x)
  values_ok && (is.null(at) || all(vapply(at, rd_no_refs, logical(1L))))
}
rd_max_len <- function(x) {
  if (is.list(x) && !is.data.frame(x)) return(if (length(x)) max(vapply(x, rd_max_len, numeric(1L))) else 0)
  if (is.data.frame(x)) return(as.double(nrow(x)))
  as.double(length(x))
}

rd_control <- function(inference_accuracy = list(max_estimate_error_abs = 1e-3, max_se_error_abs = 1e-2, max_vcov_error_abs = 1e-2,
                                                 max_multiplier_error_abs = 1e-1, max_band_error_abs = 1e-1),
                       max_work_ops = 500000000L) {
  lnk <- ipw_link_enclosure_control(max_chunk_rows = 64L, max_retained_entries = 2000000L, max_matrix_entries = 2000000L,
    max_work_ops = max_work_ops, max_output_entries = 500000L)
  ipw_pipeline_control(
    owner_control = ipw_owner_control(max_chunk_rows = 64L, max_retained_entries = 2000000L, max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops, max_output_entries = 500000L),
    catalog_bounds = list(max_candidates = 500L, max_rows = 500L, max_pool_entries = 5000L, max_scratch = 100000L, max_work = max_work_ops),
    class_bounds = list(max_links = 500L, max_classes = 50L, max_class_dates = 50L, max_E_entries = 5000L, max_mask_entries = 100000L),
    fitted_control = ipw_fitted_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = max_work_ops, max_output_entries = 500000L),
    point_control = ipw_point_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L, max_work_ops = max_work_ops, max_output_entries = 500000L),
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    certificate_control = logit_certificate_control(taylor_degree = 12L, max_abs_eta = 50.0, max_squarings = 16L, max_chunk_rows = 64L,
      max_matrix_entries = 2000000L, max_work_ops = max_work_ops),
    problem_control = ipw_nuisance_problem_control(max_chunk_rows = 64L, max_retained_entries = 2000000L, max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops),
    link_control = lnk,
    covariance_control = ipw_covariance_transfer_control(max_chunk_rows = 64L, max_operator_block_rows = 64L, max_cluster_rows = 2000L,
      max_matrix_entries = 2000000L, max_work_ops = max_work_ops, max_output_entries = 1000000L, link_control = lnk),
    gmm_accuracy = gmm_accuracy_control(method = "certified_eigen_chol_kkt", target_error_unit = "abs", max_target_error_abs = 1e-4,
      max_target_error_se = 1e-2, max_orthogonality_defect = 1e-8, min_certified_lambda_rel = 1e-10, max_factor_residual_rel = 1e-10,
      min_rcond_R = 1e-10, min_rcond_U = 1e-10, min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10, max_loading_residual_rel = 1e-10,
      require_matching_scale = TRUE, max_loading_error_abs = 1e-4),
    gmm_bounds = gmm_bounds_control(max_matrix_entries = 2000000L, max_batch_rows = 64L, max_targets = 100L, max_open_clusters = 2000L,
      allow_materialize_if = TRUE, max_materialized_if_entries = 1000000L),
    multiplier_budget = list(max_units = 2000, max_clusters = 2000, max_targets = 100, max_draws = 2000, max_if_batch_rows = 64,
      max_if_batch_cells = 6400, max_sign_batch_rows = 32, max_sign_batch_cells = 64000, max_workspace_cells = 1000000, max_workspace_bytes = 8000000),
    inference_accuracy = inference_accuracy)
}

# Deterministic 24-unit, 4-period panel (cohorts 3, 4, never): covariate values
# repeat identically across cohorts (no separation) and the outcome is a fixed
# formula with a bounded sine term (no RNG anywhere).
rd_panel <- function() {
  n_u <- 24L; u_g <- rep(c(3, 4, Inf), each = 8L); u_x <- rep(c(-1, -0.5, 0, 0.5, 1, 1.5, -1.5, 0.25), 3L)
  rows <- vector("list", n_u * 4L); k <- 0L
  for (i in seq_len(n_u)) for (tt in 1:4) {
    k <- k + 1L
    eff <- if (u_g[i] == 3 && tt >= 3) 2.0 + 0.5 * (tt - 3) else if (u_g[i] == 4 && tt >= 4) 1.5 else 0.0
    rows[[k]] <- data.frame(id = i, t = tt, g = u_g[i], x = u_x[i], y = 0.25 * i + 0.5 * tt + 0.3 * u_x[i] * tt + eff + 0.1 * sin(i * tt))
  }
  df <- do.call(rbind, rows)
  df$y_const <- ifelse(df$g == 3 & df$t >= 3, 1.0, 0.0)
  df
}
rd_run <- function(r, accuracy = NULL, yname = "y", ...) {
  ctrl <- if (is.null(accuracy)) rd_control() else rd_control(inference_accuracy = accuracy)
  tryCatch(run_chained_ipw_pipeline(data = rd_panel(), idname = "id", tname = "t", yname = yname, gname = "g", xformla = ~x,
    control = ctrl, control_group = "notyettreated", base_period = "varying", weighting = "optimal",
    reduce_structural_redundancy = TRUE, run_numerical_transfer = TRUE,
    root_witness = function(fv) list(F = unname(fv$factors$H_chol), r = r), cband = FALSE, bstrap = FALSE, ...),
    error = function(e) e)
}
rd_reasons_from <- function(rr) c(if (!rr$inequalities$a_gt_0) "a_not_positive", if (!rr$inequalities$delta_lt_1) "delta_ge_1",
  if (!rr$inequalities$x_lt_1) "x_ge_1", if (!rr$inequalities$q_lt_1) "q_ge_1", if (!rr$inequalities$selfmap_le_r) "selfmap_gt_r",
  if (!is.finite(rr$rho)) "rho_nonfinite")

# ------------------------------------------------------------------------------
# A. Helpers: identical re-signal, scalar-only packets, reservation seam
# ------------------------------------------------------------------------------

test_that("the re-signal helper returns the identical condition with one added packet and never swallows anything", {
  mk <- function(cls, diag = list(a = 1)) structure(class = c(cls, "cdid_error", "error", "condition"), list(message = "m", call = NULL, diagnostics = diag))
  e <- tryCatch(ipw_pipeline_resignal_with_state(stop(mk(c("cdid_gmm_pd_error", "cdid_gmm_error"))), "cdid_gmm_error", "extra",
    function(e) list(b = 2, env = globalenv(), f = identity, df = data.frame(k = 1:2))), error = identity)
  expect_identical(class(e), c("cdid_gmm_pd_error", "cdid_gmm_error", "cdid_error", "error", "condition"))
  expect_identical(conditionMessage(e), "m"); expect_null(conditionCall(e))
  expect_identical(e$diagnostics$a, 1); expect_identical(names(e$diagnostics), c("a", "extra"))
  expect_identical(e$diagnostics$extra, list(b = 2, df = data.frame(k = 1:2)))    # environment and function dropped
  # non-matching classes pass through untouched (identical object)
  other <- mk("cdid_ipw_pipeline_error")
  e2 <- tryCatch(ipw_pipeline_resignal_with_state(stop(other), "cdid_gmm_error", "extra", function(e) list(b = 2)), error = identity)
  expect_identical(e2, other)
  e3 <- tryCatch(ipw_pipeline_resignal_with_state(stop("plain"), "cdid_gmm_error", "extra", function(e) list(b = 2)), error = identity)
  expect_s3_class(e3, "simpleError"); expect_identical(conditionMessage(e3), "plain"); expect_null(e3$diagnostics)
  # an existing field is never overwritten; a failing state function still re-signals
  kept <- tryCatch(ipw_pipeline_resignal_with_state(stop(mk("cdid_gmm_error", list(extra = "old"))), "cdid_gmm_error", "extra", function(e) list(b = 2)), error = identity)
  expect_identical(kept$diagnostics$extra, "old")
  bad <- tryCatch(ipw_pipeline_resignal_with_state(stop(mk("cdid_gmm_error")), "cdid_gmm_error", "extra", function(e) stop("state failed")), error = identity)
  expect_identical(conditionMessage(bad), "m"); expect_true(isTRUE(bad$diagnostics$extra$retention_failed))
  expect_match(bad$diagnostics$extra$reason, "state failed")
  # the diagnostics slot may be absent on the original condition
  nodiag <- tryCatch(ipw_pipeline_resignal_with_state(stop(structure(class = c("cdid_gmm_error", "cdid_error", "error", "condition"),
    list(message = "n", call = NULL))), "cdid_gmm_error", "extra", function(e) list(b = 2)), error = identity)
  expect_identical(nodiag$diagnostics, list(extra = list(b = 2)))
  # success values pass through identically and the state function is never evaluated
  touched <- FALSE
  v <- ipw_pipeline_resignal_with_state(list(x = 1:3), "cdid_gmm_error", "extra", function(e) { touched <<- TRUE; list() })
  expect_identical(v, list(x = 1:3)); expect_false(touched)
  # mutation: dropping the re-signal would silently swallow the refusal
  mutant <- ipw_pipeline_resignal_with_state
  txt <- paste(deparse(body(mutant)), collapse = "\n"); changed <- sub("stop(e)", "invisible(NULL)", txt, fixed = TRUE)
  expect_false(identical(changed, txt)); body(mutant) <- parse(text = changed)[[1L]]
  expect_null(mutant(stop(mk("cdid_gmm_error")), "cdid_gmm_error", "extra", function(e) list(b = 2)))
  expect_error(ipw_pipeline_resignal_with_state(stop(mk("cdid_gmm_error")), "cdid_gmm_error", "extra", function(e) list(b = 2)), class = "cdid_gmm_error")
})

test_that("the scalar-state sanitizer keeps atomic/table state and drops every reference", {
  x <- list(a = 1:3, m = matrix(1, 2, 2), e = globalenv(), f = sum, d = data.frame(k = 1:2, z = I(list(sum, 2))),
    nested = list(p = "s", q = new.env(), r = list(s = TRUE, t = identity)), cls = structure(list(u = 1), class = "some_class"))
  z <- ipw_pipeline_scalar_state(x)
  expect_true(rd_no_refs(z))
  expect_identical(names(z), c("a", "m", "d", "nested", "cls"))
  expect_identical(z$a, 1:3); expect_identical(z$m, matrix(1, 2, 2))
  expect_identical(z$nested, list(p = "s", r = list(s = TRUE)))
  expect_null(attr(z$cls, "class")); expect_identical(z$cls$u, 1)
  expect_identical(ncol(z$d), 2L); expect_identical(z$d$z[[2L]], 2); expect_null(z$d$z[[1L]])
  expect_null(ipw_pipeline_scalar_state(globalenv())); expect_identical(ipw_pipeline_scalar_state(NULL), NULL)
})

test_that("the nuisance packet reports the reservation seam and keeps the minimal valid packet when it cannot fit", {
  cert <- list(delta = 1e-14, t = 2e-13, a = 4.1, x = 4e-6, q = 5e-6, selfmap = 3e-12, rho = 2e-13, r = 1e-6,
    radius_selection = list(mode = "caller_supplied", branch = NULL, inside_sufficient_region = NULL),
    inequalities = list(a_gt_0 = TRUE, delta_lt_1 = TRUE, x_lt_1 = TRUE, q_lt_1 = TRUE, selfmap_le_r = FALSE),
    dimensions = list(by_cohort_risk_count = c(5L, 6L, 7L), by_cohort_observed_count = c(4L, 5L, 6L)),
    work = list(actual_interval_ops = 100, actual_reader_work_ops = 50, endpoints_read = 0L))
  prob <- list(certified = FALSE, inconclusive_type = "inequality", raw = list(r = 1e-6), certificate = cert,
    arithmetic_refusal = list(occurred = FALSE, reason = "none"))
  fv <- list(nuisance_id = "nuisance:op9:r9",
    diagnostics = list(converged = TRUE, boundary = FALSE, iter = 7L, rank = 2L, n_covariates = 2L, df_residual = 16L,
      deviance = 1.5, null_deviance = 2.5, aic = 5.5, control = list(epsilon = 1e-12, maxit = 50L)),
    counts = list(n_a_units = 24L, n_risk_a_units = 18L, n_outside_risk_a_units = 6L, n_observed_pair_units = 18L,
      n_observed_risk_units = 15L, n_nonobserved_risk_units = 3L, n_treated_obs = 8L, n_control_obs = 7L),
    dimensions = list(m_cohorts = 3L))
  witness <- list(F = diag(2), r = 1e-6)
  expected_radius <- list(nuisance_id = "nuisance:op9:r9", mode = "caller_supplied", branch = NA_character_, r = 1e-6, delta = 1e-14,
    t = 2e-13, a = 4.1, x = 4e-6, q = 5e-6, selfmap = 3e-12, rho = 2e-13, inside_sufficient_region = NA, certified = FALSE,
    inconclusive_type = "inequality", inequalities = cert$inequalities, arithmetic_refusal = list(occurred = FALSE, reason = "none"))
  full <- ipw_pipeline_nuisance_refusal_state(prob, fv, witness, "link:r9", 1L, 3L)   # one nuisance certified before this one
  expect_identical(full$radius_record, expected_radius)
  expect_identical(full$resource$retained_packet, "full"); expect_identical(full$resource$seam, "none")
  expect_identical(full$resource$full_packet_entries, 54); expect_identical(full$resource$allowance_entries, 72)
  expect_identical(full$fit_diagnostics$iter, 7); expect_identical(full$counts$by_cohort_risk_count, c(5L, 6L, 7L))
  expect_identical(full$witness$r_request, 1e-6); expect_identical(full$witness$r_request_kind, "explicit_numeric_radius")
  expect_identical(full$work$actual_interval_ops, 100); expect_identical(full$link_id, "link:r9"); expect_identical(full$nuisance_index, 2L)
  # one link: the 24-entry reservation cannot hold the 54-entry packet
  minimal <- ipw_pipeline_nuisance_refusal_state(prob, fv, witness, "link:r9", 1L, 1L)
  expect_identical(minimal$radius_record, expected_radius)
  expect_identical(minimal$resource$retained_packet, "minimal"); expect_identical(minimal$resource$packet_entries, 21)
  expect_identical(minimal$resource$allowance_entries, 24); expect_match(minimal$resource$seam, "exceeds the 24 x n_links = 24 reservation")
  expect_false(any(c("fit_diagnostics", "counts", "witness", "work") %in% names(minimal)))
  expect_true(rd_no_refs(full)); expect_true(rd_no_refs(minimal))
  # the automatic request object is recorded by kind, never as a number
  auto <- ipw_pipeline_nuisance_refusal_state(prob, fv, list(F = diag(2), r = logit_cert_automatic_radius_request()), "l", 1L, 3L)
  expect_true(is.na(auto$witness$r_request)); expect_false(identical(auto$witness$r_request_kind, "explicit_numeric_radius"))
})

# ------------------------------------------------------------------------------
# B. GMM core: Cholesky refusal on a three-link cycle design (direct call)
# ------------------------------------------------------------------------------

test_that("the Cholesky refusal keeps its class/message/diagnostics and retains factorization, design-rank and moment state", {
  des <- build_chained_design(data.frame(g = 4, s = c(1, 2, 1), t = c(2, 3, 3), ref = 3, functional_key = "A"))
  g <- build_gmm_gauge(des, "lowest")
  expect_identical(des$incidence$n_parameters, 2L); expect_identical(des$rank_diagnostics$rank, 2L); expect_identical(g$n_links, 3L)
  bounds <- gmm_bounds_control(max_matrix_entries = 50000L, max_batch_rows = 16L, max_targets = 10L, max_open_clusters = 100L,
    allow_materialize_if = TRUE, max_materialized_if_entries = 10000L)
  acc <- gmm_accuracy_control(method = "certified_eigen_chol_kkt", target_error_unit = "abs", max_target_error_abs = 1e-4,
    max_target_error_se = 1e-2, max_orthogonality_defect = 1e-8, min_certified_lambda_rel = 1e-10, max_factor_residual_rel = 1e-10,
    min_rcond_R = 1e-10, min_rcond_U = 1e-10, min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10, max_loading_residual_rel = 1e-10,
    require_matching_scale = TRUE, max_loading_error_abs = 1e-4)
  prov <- gmm_covariance_provenance(row_level = "unit", n_units = 10L, n_clusters = 10L, n_positive_contributing_rows = 10L,
    centered = TRUE, weighted = FALSE, weight_type = "unweighted", normalization_factor = 0.1, scale_id = "scale:test", nuisance_fit_id = "nuisance_fit:test")
  Phi <- cbind(seq(-1, 1, length.out = 10), c(0.3, -0.2, 0.5, -0.4, 0.1, 0.2, -0.3, 0.4, -0.5, -0.1), c(1, -1, 0.5, -0.5, 0.25, -0.25, 0.75, -0.75, 0.1, -0.1))
  Phi <- sweep(Phi, 2L, colMeans(Phi)); called <- 0L
  replay <- gmm_if_replay_spec(callback = function(consume) { called <<- called + 1L
      consume(list(unit_id = 1:10, phi = Phi, weight = rep(1, 10), nuisance_fit_id = "nuisance_fit:test")) },
    input_row_level = "unit", aggregate_to_clusters = FALSE, n_input_rows = 10L, n_output_rows = 10L, n_positive_contributing_rows = 10L,
    centered = TRUE, weighted = FALSE, weight_type = "unweighted", normalization_factor = 0.1, scale_id = "scale:test",
    nuisance_fit_id = "nuisance_fit:test", expected_unit_ids = 1:10)
  targets <- list(gmm_target_spec("level_1", 1L, 1.0), gmm_target_spec("level_2", 2L, 1.0))
  solve_with <- function(V, mhat = c(0.1, 0.2, 0.3)) {
    names(mhat) <- g$link_ids; dimnames(V) <- list(g$link_ids, g$link_ids)
    tryCatch(solve_chained_gmm_core(design = des, mhat = mhat, V = V, provenance = prov, targets = targets, weighting = "optimal",
      if_replay = replay, accuracy = acc, bounds = bounds, gauge_rule = "lowest", request_param_cov = TRUE, materialize_target_if = TRUE),
      error = function(e) e)
  }
  V1 <- matrix(1, 3, 3)                       # rank one: the order-2 leading minor is zero
  e <- solve_with(V1)
  expect_identical(class(e), c("cdid_gmm_pd_error", "cdid_gmm_refusal", "cdid_gmm_error", "cdid_error", "error", "condition"))
  expect_identical(conditionMessage(e), "Unpivoted Cholesky factorization failed on `V`; no fallback or regularization is permitted.")
  expect_identical(called, 0L)                 # the refusal still precedes the IF replay
  expect_identical(names(e$diagnostics), c("m_R", "k", "refusal_state"))
  expect_identical(e$diagnostics$m_R, 3L); expect_identical(e$diagnostics$k, 2L)
  rs <- e$diagnostics$refusal_state
  expect_true(rd_no_refs(rs))
  expect_identical(rs$stage, "optimal_weight_cholesky_whitening")
  expect_identical(rs$factorization$outcome, "factorization_error"); expect_match(rs$factorization$error_message, "leading minor of order 2")
  expect_identical(rs$dimensions[c("m_R", "p", "k", "nullity", "n_targets", "J_est", "M_positive_rows")],
    list(m_R = 3L, p = 2L, k = 2L, nullity = 0L, n_targets = 2L, J_est = 2L, M_positive_rows = 10L))
  expect_match(rs$dimensions$meaning, "None of these is the rank of the m_R x m_R sample covariance")
  expect_identical(rs$design_rank$cycle_rank, 1L); expect_identical(rs$design_rank$n_post_targets, 0L); expect_true(is.data.frame(rs$design_rank$by_cohort))
  expect_identical(rs$covariance$diag_min, 1); expect_identical(rs$covariance$diag_max, 1); expect_true(rs$covariance$exactly_symmetric)
  expect_identical(rs$covariance$link_ids, g$link_ids); expect_identical(rs$covariance$provenance$nuisance_fit_id, "nuisance_fit:test")
  expect_identical(unname(rs$moment_state$mhat), c(0.1, 0.2, 0.3)); expect_identical(names(rs$moment_state$mhat), g$link_ids)
  expect_identical(unname(rs$moment_state$V), V1)
  expect_identical(rs$moment_state$W_G, materialize_gmm_gauge_matrix(g, 50000L))
  expect_identical(rs$moment_state$keep_param_ids, g$param_ids); expect_identical(rs$moment_state$target_ids, c("level_1", "level_2"))
  expect_identical(rs$retained_payload$entries, 3^2 + 3 * 2 + 3)
  expect_match(rs$meaning, "neither exact singularity of V, nor the population covariance rank, nor ATT nonidentification")
  # the packet reflects the actual inputs (a different moment vector is carried verbatim)
  e_b <- solve_with(V1, mhat = c(-0.4, 0.25, 0.05))
  expect_identical(unname(e_b$diagnostics$refusal_state$moment_state$mhat), c(-0.4, 0.25, 0.05))
  # a positive definite covariance passes the unchanged factorization and the solve completes through the replay
  ok <- solve_with(diag(3))
  expect_false(inherits(ok, "error")); expect_true(all(is.finite(ok$targets$estimate))); expect_gte(called, 1L)
  # the neighbouring nonpositive-diagonal refusal is unchanged (no Cholesky packet there)
  V0 <- diag(c(1, 0, 1)); e0 <- solve_with(V0)
  expect_identical(conditionMessage(e0), "Weight covariance `V` has nonpositive diagonal entries.")
  expect_identical(names(e0$diagnostics), "min_diag_V"); expect_identical(e0$diagnostics$min_diag_V, 0)
})

# ------------------------------------------------------------------------------
# C. Pipeline stages on the deterministic fixture (unit fixture fits, no RNG)
# ------------------------------------------------------------------------------

test_that("the nuisance root INCONCLUSIVE refusal carries certificate scalars, fit diagnostics and risk counts consistently", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  tiny <- rd_run(1e-30)
  expect_identical(class(tiny), c("cdid_ipw_pipeline_transfer_inconclusive", "cdid_logit_certificate_inconclusive",
    "cdid_ipw_pipeline_error", "cdid_error", "error", "condition"))
  expect_identical(conditionMessage(tiny), "Nuisance root certificate is INCONCLUSIVE for `nuisance:op1:r1`.")
  d <- tiny$diagnostics
  expect_identical(names(d), c("nuisance_id", "reasons", "refusal_state", "pipeline_state"))
  expect_identical(d$nuisance_id, "nuisance:op1:r1"); expect_identical(d$reasons, "selfmap_gt_r")
  rs <- d$refusal_state; rr <- rs$radius_record
  expect_true(rd_no_refs(d))
  expect_identical(rs$stage, "nuisance_root_certificate"); expect_identical(rs$nuisance_index, 1L); expect_identical(rs$n_links, 8L)
  expect_identical(rr$nuisance_id, "nuisance:op1:r1"); expect_identical(rr$r, 1e-30); expect_false(rr$certified)
  expect_identical(rr$inconclusive_type, "inequality"); expect_identical(rr$mode, "caller_supplied")
  expect_identical(rd_reasons_from(rr), d$reasons)        # reasons are exactly the failed inequalities
  expect_true(all(is.finite(c(rr$delta, rr$t, rr$a, rr$x, rr$q, rr$selfmap, rr$rho))))
  expect_gt(rr$selfmap, rr$r); expect_lt(rr$x, 1); expect_lt(rr$q, 1); expect_gt(rr$a, 0)
  expect_true(rs$fit_diagnostics$converged); expect_identical(rs$fit_diagnostics$n_covariates, 2)
  expect_identical(rs$fit_diagnostics$glm_control, list(epsilon = 1e-12, maxit = 50))
  expect_identical(rs$counts$n_a_units, 24); expect_identical(rs$counts$m_cohorts, 3)
  expect_identical(sum(rs$counts$by_cohort_risk_count), as.integer(rs$counts$n_risk_a_units))
  expect_identical(rs$counts$n_risk_a_units + rs$counts$n_outside_risk_a_units, 24)
  expect_identical(rs$witness$r_request, 1e-30); expect_identical(rs$witness$r_request_kind, "explicit_numeric_radius")
  expect_identical(rs$witness$F_dim, c(2L, 2L))
  expect_identical(rs$resource$retained_packet, "full"); expect_identical(rs$resource$packet_entries, 54)
  expect_identical(rs$resource$allowance_entries, 24 * 8); expect_identical(rs$resource$seam, "none")
  expect_match(rs$meaning, "not evidence against finite-root existence")
  ps <- d$pipeline_state
  expect_identical(ps$stage, "numerical_transfer"); expect_identical(ps$weighting, "optimal"); expect_identical(ps$n_a_units, 24L)
  expect_identical(ps$catalog_class_map$n_links_all, 15L); expect_identical(ps$catalog_class_map$n_links_used, 8L)
  expect_identical(length(ps$catalog_class_map$kept_indices), 8L)
  expect_identical(ps$catalog_class_map$kept_indices, ps$catalog_class_map$class_certificate$retained_index)
  expect_identical(nrow(ps$catalog_class_map$catalog_rows_used), 8L); expect_identical(ps$covariance_stream$n_moments, 8L)
  expect_identical(ps$covariance_stream$link_ids, ps$catalog_class_map$design_link_ids_used)
  expect_identical(ps$design_rank$n_links, 8L); expect_identical(ps$gmm_core$n_estimable, 7L); expect_true(ps$gmm_core$point_certified)
  expect_true(all(c("cell_id", "available", "support_status") %in% names(ps$ledger_support)))
  # a large radius fails the P1-type reason set; the packet tracks every reason from the recorded inequalities
  big <- rd_run(10)
  expect_identical(class(big), class(tiny)); expect_identical(conditionMessage(big), conditionMessage(tiny))
  rb <- big$diagnostics$refusal_state$radius_record
  expect_identical(rb$r, 10); expect_setequal(big$diagnostics$reasons, c("x_ge_1", "q_ge_1", "selfmap_gt_r", "rho_nonfinite"))
  expect_identical(rd_reasons_from(rb), big$diagnostics$reasons); expect_false(is.finite(rb$rho))
  # the same fixture certifies at the accepted radius: the refusal is radius-driven and the record is the same scalar family
  ok <- rd_run(1e-6)
  expect_false(inherits(ok, "error")); expect_true(ok$numerical_transfer$certified)
  radii <- ok$numerical_transfer$nuisance_root_radii
  expect_identical(radii$nuisance_id[1L], "nuisance:op1:r1"); expect_identical(radii$r[1L], 1e-6)
  expect_identical(names(radii)[c(1L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L)], c("nuisance_id", "r", "delta", "t", "a", "x", "q", "selfmap", "rho"))
  expect_identical(c(rr$delta, rr$t, rr$a), c(radii$delta[1L], radii$t[1L], radii$a[1L]))   # radius-free scalars agree across requests
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the D20 refusal carries target labels, the accepted certificate's error bounds and the Theorem 12 attribution", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  ref <- rd_run(1e-6, accuracy = list(max_relative_error_lower_se = 1e-3))
  expect_false(inherits(ref, "error")); d20_ok <- ref$reporting$numerical_certificate$d20; expect_true(d20_ok$certified)
  e <- rd_run(1e-6, accuracy = list(max_relative_error_lower_se = 1e-15))
  expect_identical(class(e), c("cdid_ipw_pipeline_d20_error", "cdid_ipw_pipeline_accuracy_error", "cdid_ipw_pipeline_error", "cdid_error", "error", "condition"))
  expect_match(conditionMessage(e), paste0("^D20 scale-aware numerical acceptance refused: the largest certified error ratio \\(", format(max(d20_ok$worst_ratios), digits = 6),
    " for `", names(d20_ok$worst_ratios)[which.max(d20_ok$worst_ratios)], "`\\) exceeds the level 1e-15; the estimator, sample, weights and inference formulas are unchanged and no relaxation is applied\\.$"))
  d <- e$diagnostics
  expect_identical(names(d)[1:6], c("level", "lower_se_scale", "worst_ratios", "ratios_by_target", "covariance_ratios", "multiplier_ratios"))
  expect_identical(d$level, 1e-15)
  # the refusal packet equals the content the accepted certificate retains on success (same inputs, only the level differs)
  expect_identical(d$worst_ratios, d20_ok$worst_ratios); expect_identical(d$ratios_by_target, d20_ok$ratios_by_target)
  expect_identical(d$covariance_ratios, d20_ok$covariance_ratios); expect_identical(d$lower_se_scale, d20_ok$lower_se_scale)
  expect_identical(d$errors_by_target, d20_ok$errors_by_target); expect_identical(d$covariance_errors, d20_ok$covariance_errors)
  expect_identical(d$multiplier_errors_by_target, d20_ok$multiplier_errors_by_target); expect_identical(d$raw_draw_scale_lower, d20_ok$raw_draw_scale_lower)
  ids <- ref$reporting$numerical_certificate$certified_target_ids
  expect_identical(d$target_ids, ids); expect_identical(colnames(d$ratios_by_target), ids); expect_identical(names(d$lower_se_scale), ids)
  expect_identical(d$structural_zero_targets, character(0L)); expect_identical(d$n_total, 24L); expect_identical(d$n_clusters, 24L)
  expect_identical(rownames(d$errors_by_target), c("estimate", "se", "ci_lower", "ci_upper", "band_lower", "band_upper"))
  # outward ratios bound the plain quotients of the retained errors and scales
  plain <- sweep(d$errors_by_target, 2L, d$lower_se_scale, "/")
  expect_true(all(d$ratios_by_target >= plain)); expect_true(all(d$ratios_by_target <= plain * (1 + 1e-12) + 1e-300))
  expect_identical(d$worst_target_by_family[["ci_lower"]], ids[which.max(d$ratios_by_target["ci_lower", ])])
  expect_true(all(d$worst_covariance_entry %in% ids)); expect_identical(dim(d$se_box), c(2L, length(ids)))
  expect_identical(d$pointwise_constant, stats::qnorm(0.975)); expect_match(d$refusal_meaning, "ATT identification")
  # Theorem 12 attribution of the already certified transfer
  at <- d$error_attribution
  expect_true(rd_no_refs(d)); expect_identical(at$stage, "reported_cell_inference"); expect_true(at$available)
  expect_identical(at$theorem_ref, "thm:gmm"); expect_identical(at$weighting, "optimal")
  expect_true(is.finite(at$epsilon_m) && at$epsilon_m >= 0); expect_true(is.finite(at$epsilon_V) && at$epsilon_V >= 0)
  expect_gt(at$l_lower, 0); expect_gt(at$l_minus_eps_V, 0); expect_lte(at$l_minus_eps_V, at$l_lower)
  expect_identical(at$certified_lower_eigenvalue$value, at$l_lower); expect_true(at$certified_lower_eigenvalue$applies)
  J_est <- 7L
  expect_identical(length(at$loading_target_ids), J_est)
  # post cells are loading columns themselves; pre cells are base-period contrasts of loading columns
  expect_true(all(ids[startsWith(ids, "ATT(")] %in% at$loading_target_ids)); expect_identical(sum(startsWith(ids, "ATT(")), 3L)
  pc <- at$point_error_components; rc <- at$R_c_components
  expect_identical(length(pc$solve), J_est); expect_identical(length(rc$changing_weights), J_est)
  comp_sum <- pc$solve + pc$moment_change + pc$changing_weights
  expect_true(all(at$total_target_error_cert >= pmax(pc$solve, pc$moment_change, pc$changing_weights)))
  expect_true(all(at$total_target_error_cert <= comp_sum * (1 + 1e-12)))
  rc_sum <- rc$projection + rc$row_if_change + rc$numerical_loading + rc$changing_weights
  expect_true(all(at$R_c_cert >= pmax(rc$projection, rc$row_if_change, rc$numerical_loading, rc$changing_weights)))
  expect_true(all(at$R_c_cert <= rc_sum * (1 + 1e-12)))
  expect_true(all(at$total_se_error_cert >= at$R_c_cert))
  split <- at$moment_error_split
  expect_identical(nrow(split), 8L); expect_identical(split$total_error_bound, split$epsilon_m_vec)
  expect_true(all(split$total_error_bound >= pmax(split$eval_error, split$root_error)))
  expect_true(all(split$total_error_bound <= (split$eval_error + split$root_error) * (1 + 1e-12)))
  expect_gte(at$epsilon_m, sqrt(sum(split$total_error_bound^2)) * (1 - 1e-12))
  expect_identical(at$nuisance_root_radii$r, rep(1e-6, 3L))
  expect_identical(at$total_target_error_cert, ref$numerical_transfer$gmm_transfer$total_target_error_cert)
  # no influence-function array travels with the refusal (rows x targets would be 24 x 7)
  expect_lt(rd_max_len(d), 24 * 7)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("a GMM-core refusal inside the pipeline keeps the core's own fields and gains the catalog/class map context", {
  e <- rd_run(1e-6, yname = "y_const")
  expect_identical(class(e), c("cdid_gmm_pd_error", "cdid_gmm_refusal", "cdid_gmm_error", "cdid_error", "error", "condition"))
  expect_identical(conditionMessage(e), "Weight covariance `V` has nonpositive diagonal entries.")
  expect_identical(names(e$diagnostics), c("min_diag_V", "pipeline_state")); expect_identical(e$diagnostics$min_diag_V, 0)
  ps <- e$diagnostics$pipeline_state
  expect_true(rd_no_refs(ps)); expect_identical(ps$stage, "gmm_core_solve"); expect_null(ps$gmm_core)
  expect_identical(ps$control_group, "notyettreated"); expect_identical(ps$base_period, "varying"); expect_identical(ps$n_total, 24L)
  cm <- ps$catalog_class_map
  expect_identical(cm$n_links_all, 15L); expect_identical(cm$n_links_used, 8L); expect_true(cm$reduction_active)
  expect_identical(cm$kept_indices, cm$class_certificate$retained_index); expect_identical(cm$dropped_indices, cm$class_certificate$dropped_index)
  expect_identical(sort(c(cm$kept_indices, cm$dropped_indices)), seq_len(15L))
  expect_identical(nrow(cm$catalog_rows_used), 8L); expect_identical(cm$class_certificate$n_links, 15L)
  expect_identical(ps$covariance_stream$provenance$n_units, 24L); expect_identical(ps$covariance_stream$n_moments, 8L)
  expect_identical(ps$design_rank$n_links, 8L); expect_true(is.data.frame(ps$core_targets)); expect_true(is.data.frame(ps$ledger_support))
  expect_true(all(is.finite(c(ps$resource$simultaneous_entries, ps$resource$retained_with_diagnostics, ps$resource$downstream_work_ops))))
})

test_that("the development request entry re-signals the pipeline refusal unchanged with request-level context", {
  d <- rd_panel(); d$g[is.infinite(d$g)] <- 0
  ctrl <- rd_control()
  e <- tryCatch(chained_att_gt_request(yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~x, data = d,
    control_group = "notyettreated", est_method = "ipw", bstrap = FALSE, cband = FALSE, control = ctrl,
    root_witness = function(fv) list(F = unname(fv$factors$H_chol), r = 1e-30)), error = function(e) e)
  expect_identical(class(e), c("cdid_ipw_pipeline_transfer_inconclusive", "cdid_logit_certificate_inconclusive",
    "cdid_ipw_pipeline_error", "cdid_error", "error", "condition"))
  expect_identical(conditionMessage(e), "Nuisance root certificate is INCONCLUSIVE for `nuisance:op1:r1`.")
  expect_identical(names(e$diagnostics), c("nuisance_id", "reasons", "refusal_state", "pipeline_state", "request_context"))
  rc <- e$diagnostics$request_context
  expect_true(rd_no_refs(rc)); expect_identical(rc$entry, "development")
  expect_identical(rc$input[c("n_frame", "n_frame_only", "n_panel_rows", "n_calendar", "n_units_recoded", "clustered")],
    list(n_frame = 24L, n_frame_only = 0L, n_panel_rows = 96L, n_calendar = 4L, n_units_recoded = 8L, clustered = FALSE))
  expect_identical(rc$input$target_cohorts, c(3, 4)); expect_identical(rc$input$cohort_universe, c(3, 4, Inf))
  expect_identical(rc$settings$control_group, "notyettreated"); expect_identical(rc$settings$xformla, "~x")
  expect_identical(rc$numerical_policy$inference_accuracy, ctrl$inference_accuracy)
  expect_identical(rc$numerical_policy$gmm_accuracy$max_target_error_abs, 1e-4)
  expect_identical(rc$numerical_policy$owner_bounds$max_work_ops, 500000000)
  expect_identical(rc$numerical_policy$source, "explicit fourteen-control development object")
  # the pure attribution helper on a successful transfer certificate is the same object the MP will retain
  ok <- rd_run(1e-6)
  at <- ipw_pipeline_transfer_error_attribution(ok$numerical_transfer, ok$gmm_fit)
  expect_true(at$available); expect_identical(at$epsilon_V, ok$numerical_transfer$gmm_transfer$epsilon_V)
  expect_identical(at$loading_target_ids, ok$gmm_fit$targets$target_id[ok$gmm_fit$targets$estimable])
  expect_identical(nrow(at$moment_error_split), 8L); expect_true(rd_no_refs(ipw_pipeline_scalar_state(at)))
  expect_identical(ipw_pipeline_transfer_error_attribution(NULL)$available, FALSE)
})

# ------------------------------------------------------------------------------
# D. D20 refusal packet on a fixed matrix (independent arithmetic, no fit)
# ------------------------------------------------------------------------------

test_that("ipw_pipeline_enforce_d20 refusals carry labels and bounds with the unchanged message", {
  cert <- list(theta = list(lo = matrix(0.3995, 1L), hi = matrix(0.6005, 1L)),
    rows = list(lo = matrix(c(-1, 1), 2L), hi = matrix(c(-1, 1), 2L)),
    se = matrix(c(1, 1.01), 2L), vcov = list(lo = matrix(1.005^2), hi = matrix(1.005^2)),
    band = list(lo = matrix(0.3995, 1L), hi = matrix(0.6005, 1L)),
    band_upper = list(lo = matrix(0.3995, 1L), hi = matrix(0.6005, 1L)), deterministic = integer(), n_total = 2L, n_clusters = 2L)
  args <- list(cert = cert, estimates = c(a = 0.5), se = c(a = 1.005), vcov = matrix(1.005^2), level = 0.1,
    err_est = 0.1005, err_se = 0.005, err_v = matrix(0), band_lower = 0.5, band_upper = 0.5)
  e <- tryCatch(do.call(ipw_pipeline_enforce_d20, args), error = identity)
  expect_identical(class(e), c("cdid_ipw_pipeline_d20_error", "cdid_ipw_pipeline_accuracy_error", "cdid_ipw_pipeline_error", "cdid_error", "error", "condition"))
  ratio <- gmm_iv_ratio_upper(0.1005, 1)
  d <- e$diagnostics
  # estimate and band endpoints share the 0.1005 error against scale 1; the message names the first maximal family
  worst_family <- names(d$worst_ratios)[which.max(d$worst_ratios)]
  expect_true(worst_family %in% c("estimate", "band_lower", "band_upper"))
  expect_identical(conditionMessage(e), paste0("D20 scale-aware numerical acceptance refused: the largest certified error ratio (", format(max(d$worst_ratios), digits = 6),
    " for `", worst_family, "`) exceeds the level 0.1; the estimator, sample, weights and inference formulas are unchanged and no relaxation is applied."))
  expect_true(abs(max(d$worst_ratios) - ratio) <= 4 * .Machine$double.eps * ratio)
  expect_identical(names(d)[1:6], c("level", "lower_se_scale", "worst_ratios", "ratios_by_target", "covariance_ratios", "multiplier_ratios"))
  expect_identical(d$worst_ratios[["estimate"]], ratio); expect_identical(d$target_ids, "a")
  expect_identical(d$errors_by_target["estimate", "a"], 0.1005); expect_identical(d$errors_by_target["se", "a"], 0.005)
  expect_identical(d$covariance_errors, matrix(0, dimnames = list("a", "a"))); expect_identical(d$structural_zero_targets, character(0L))
  expect_identical(d$worst_target_by_family[["estimate"]], "a"); expect_identical(d$worst_covariance_entry, c("a", "a"))
  expect_identical(d$n_total, 2L); expect_identical(d$se_box, cert$se); expect_true(rd_no_refs(d))
  # the non-positive lower scale branch keeps its message and now names the targets too
  zero <- cert; zero$se <- matrix(c(0, 1.01), 2L)
  e0 <- tryCatch(do.call(ipw_pipeline_enforce_d20, c(list(cert = zero), args[-1L])), error = identity)
  expect_match(conditionMessage(e0), "has no certified positive lower SE scale"); expect_identical(e0$diagnostics$target, "a")
  expect_identical(e0$diagnostics$target_ids, "a")
  # the accepted outcome is unchanged: the same arithmetic passes at a level above the ratio
  ok <- do.call(ipw_pipeline_enforce_d20, c(args[names(args) != "level"], list(level = 0.2)))
  expect_true(ok$d20$certified); expect_identical(ok$d20$errors_by_target["estimate", "a"], 0.1005)
})


test_that("diagnostic state drops hidden attribute references and preserves admitted matrix payload aliases", {
  m <- matrix(as.double(1:4),2L,dimnames=list(c("r1","r2"),c("c1","c2")))
  z <- ipw_pipeline_scalar_state(list(m=m))
  expect_true(gmm_sexp_aliases(m,z$m))
  expect_identical(attributes(z$m),attributes(m))
  hidden <- m;attr(hidden,"owner") <- globalenv()
  packet <- structure(list(m=hidden),owner=globalenv())
  expect_false(rd_no_refs(packet)) # meaningful unsanitized/reference-retaining mutant
  clean <- ipw_pipeline_scalar_state(packet)
  expect_true(rd_no_refs(clean));expect_null(attr(clean,"owner"));expect_null(attr(clean$m,"owner"))
  expect_identical(clean$m,m)
  df <- structure(data.frame(a=1:2),owner=globalenv());attr(df$a,"owner") <- globalenv()
  dd <- ipw_pipeline_scalar_state(df)
  expect_true(rd_no_refs(dd));expect_identical(dd,data.frame(a=1:2))
})

test_that("minimal nuisance admission precedes omitted state and counts actual numeric/logical leaves", {
  # Independent count: radius16 + two integer headers + three resource counters
  # =21 numeric/logical entries; text metadata/containers are not a heap bound.
  count <- function(x) if(is.list(x)) sum(vapply(x,count,0)) else if(is.numeric(x)||is.logical(x)) as.double(length(x)) else 0
  prob <- list(raw=list(r=1e-6),certificate=list(),certified=FALSE)
  fv <- list(nuisance_id="n")
  small <- ipw_pipeline_nuisance_refusal_state(prob,fv,stop("omitted witness forced"),"l",0L,1L)
  expect_identical(small$resource$retained_packet,"minimal")
  expect_identical(count(small),21);expect_lte(count(small),small$resource$packet_entries)
  expect_lte(small$resource$packet_entries,small$resource$allowance_entries)
  fv$diagnostics <- list();fv$counts <- list();fv$dimensions <- list(m_cohorts=3)
  prob$certificate$dimensions <- list(by_cohort_risk_count=c(5L,6L,7L),by_cohort_observed_count=c(4L,5L,6L))
  full <- ipw_pipeline_nuisance_refusal_state(prob,fv,list(F=diag(2),r=1e-6),"l",0L,3L)
  # radius16 + headers/counters5 + fit11 + counts15 + witness3 + work3 =53;
  # the existing conservative54-entry charge covers this exact retained packet.
  expect_identical(count(full),53);expect_identical(full$resource$packet_entries,54)
  expect_lte(count(full),full$resource$packet_entries);expect_lte(full$resource$packet_entries,full$resource$allowance_entries)
  mutant <- ipw_pipeline_nuisance_refusal_state
  txt <- paste(deparse(body(mutant)),collapse="\n")
  changed <- sub("full_fits <- (n_radius + n_extended) <= allowance","full_fits <- TRUE",txt,fixed=TRUE)
  expect_false(identical(changed,txt));body(mutant) <- parse(text=changed)[[1L]]
  expect_error(mutant(prob,fv,stop("omitted witness forced"),"l",0L,1L),"omitted witness forced")
})
