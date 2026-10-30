# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 descriptor context: the single early canonical OR/DR comparison intake
# of `build_cp04_public_registry` (`cp04_descriptor_context`) and its private
# intake in `fit_or_baseline`/`fit_dr_baseline` (`or_baseline_descriptor_intake`).
# On the eligible complete-roster baseline OR/DR route the fits read one sealed
# context of all original pair/cohort/key/ID descriptors instead of M + Q heavy
# link views; every original raw pair still enters the existing common fit.
#
# UNEXECUTED at authoring time (no R in the authoring job). Every expected value
# comes from the unchanged link-view route of the same installed build, forced in
# a disclosed test-local namespace clone whose descriptor predicate declines
# (`dc_ordinary`); counter instrumentation lives only in test-local clones
# (`dc_counted`, the cb_counted pattern). Fixtures and helpers are private to
# this file and follow test-cp08-balanced-or-dr.R / test-cp08-compact-basis-construction.R.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Private fixtures
# ---------------------------------------------------------------------------

# Complete balanced long panel (od_nonadd_panel): cohort 3 against never
# treated (G = 0), unit covariate x, cluster cl, nonadditive outcome; `second`
# adds a second treated cohort 4; `w` is a non-unit unit weight column.
dc_panel <- function(N = 41L, T = 20L, n_g = 18L, seed = 31L, id_perm = FALSE, rows_perm = FALSE, second = FALSE) {
  set.seed(seed)
  g <- rep(0, N)
  g[sample.int(N, n_g)] <- 3
  if (second) g[which(g == 0)[1:6]] <- 4
  ids <- if (id_perm) sample(1000L + seq_len(5L * N), N) else seq_len(N)
  xi <- stats::rbinom(N, 1, 0.5)
  ue <- stats::rnorm(N)
  d <- data.frame(id = rep(ids, each = T), period = rep(seq_len(T), N))
  u <- match(d$id, ids)
  d$G <- g[u]
  d$x <- xi[u]
  d$w <- 1 + (u %% 3L) / 4
  eff <- ifelse(d$G > 0 & d$period >= d$G, 1 + 0.25 * (d$period - d$G), 0)
  d$y <- ue[u] + 0.1 * d$period + 0.35 * sin(0.7 * u * d$period) + 0.2 * d$x * (d$period^2) + eff + stats::rnorm(nrow(d))
  d$cl <- (u - 1L) %/% 4L + 1L
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

# Public OR/DR pipeline run (od_run_pipeline) in a namespace `ns`.
dc_run <- function(d, est = "reg", chunk = 7L, xformla = NULL, weightsname = NULL, pt_assumption = "post",
                   pt_start = NULL, control_group = "nevertreated", reduce = TRUE, weighting = "optimal",
                   aggregations = c("simple", "group", "dynamic", "calendar"), bstrap = FALSE, cband = FALSE,
                   biters = 99L, seed = 101L, ns = api, target_cohorts = NULL) {
  ctrl_list <- if (is.null(chunk)) NULL else list(chunk_rows = chunk)
  resolved <- api$chained_public_resolve_control(ctrl_list, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  prep <- api$chained_request_prepare_inputs(
    data = d, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xformla, weightsname = weightsname, clustervars = "cl", anticipation = 0,
    frame = NULL, outcome_regression = identical(est, "reg"),
    doubly_robust = identical(est, "dr"), or_resource = resolved$profile
  )
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = biters, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = built$control, xformla = prep$xformla, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = if (is.null(target_cohorts)) prep$target_cohorts else target_cohorts,
    anticipation = 0L, control_group = control_group, base_period = "varying",
    weighting = weighting, reduce_structural_redundancy = reduce, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = aggregations, cband = cband, alp = 0.05, biters = biters, seed = seed,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = bstrap, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = pt_assumption, pt_start = pt_start, cp04 = prep$cp04
  )
}

# ---------------------------------------------------------------------------
# Instrumentation (disclosed test-local clones of the namespace; cb_counted)
# ---------------------------------------------------------------------------

dc_counted <- function(wrap = character()) {
  e <- new.env(parent = api)
  counts <- new.env(parent = emptyenv())
  genuine <- "gmm_require_native_dll"
  for (nm in ls(api, all.names = TRUE)) {
    f <- get(nm, envir = api, inherits = FALSE)
    if (is.function(f) && !is.primitive(f) && identical(environment(f), api) && !nm %in% genuine) {
      environment(f) <- e
      assign(nm, f, envir = e)
    }
  }
  for (nm in genuine) assign(nm, get(nm, envir = api, inherits = FALSE), envir = e)
  make_wrap <- function(nm, inner) {
    force(nm); force(inner)
    function(...) {
      assign(nm, get(nm, envir = counts, inherits = FALSE) + 1, envir = counts)
      inner(...)
    }
  }
  for (nm in wrap) {
    if (!exists(nm, envir = e, inherits = FALSE)) stop("No internal function `", nm, "` to count.", call. = FALSE)
    assign(nm, 0, envir = counts)
    assign(nm, make_wrap(nm, get(nm, envir = e, inherits = FALSE)), envir = e)
  }
  seal_inner <- get("ipw_owner_seal", envir = e, inherits = FALSE)
  assign("ipw_owner_seal", function(obj, binding) {
    key <- paste0("seal:", binding)
    assign(key, get0(key, envir = counts, inherits = FALSE, ifnotfound = 0) + 1, envir = counts)
    seal_inner(obj, binding)
  }, envir = e)
  e$.counts <- counts
  e
}
dc_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0)

# The unchanged link-view route of the same build: only the descriptor
# predicate declines (its call is still counted); everything else is installed code.
dc_ordinary <- function(wrap = character()) {
  e <- dc_counted(wrap)
  assign("cp04_descriptor_context", function(...) {
    assign("cp04_descriptor_context", dc_n(e, "cp04_descriptor_context") + 1, envir = e$.counts)
    NULL
  }, envir = e)
  e
}

# Capture the predicate's actual arguments and sealed result inside clone `e`.
dc_capture <- function(e, trace) {
  inner <- get("cp04_descriptor_context", envir = e, inherits = FALSE)
  assign("cp04_descriptor_context", function(metadata, frame, reader, specs, aux, ctrl, ...) {
    trace$args <- c(list(metadata = metadata, frame = frame, reader = reader, specs = specs, aux = aux, ctrl = ctrl),
      list(...))
    out <- inner(metadata, frame, reader, specs, aux, ctrl, ...)
    trace$ctx <- out
    out
  }, envir = e)
  invisible(e)
}

dc_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, dc_strip); attributes(x) <- a }
  x
}
# Route-only accounting differences stay out; every fitted field, record and
# diagnostic (supports, regressions, links, shared fit, treatments) is compared.
dc_backend_view <- function(b) {
  x <- unclass(b)[setdiff(names(b), c("ref", "reader_ref", "accounting"))]
  list(fields = dc_strip(x), accounting = b$accounting[setdiff(names(b$accounting),
    c("work_ops_bound", "fit_scratch_bound", "retained_bound", "owned_entries", "retained_external_entries",
      "external_work_ops", "borrowed_entries"))])
}
# Same strict schema/leaf predicates as od_expect_route_match at the established tolerance.
dc_route_match <- function(a, b, tol = 1e-12, info = NULL) {
  mismatch <- function(x, y, path) {
    if (!identical(class(x), class(y)) || !identical(typeof(x), typeof(y)) ||
        !identical(names(x), names(y)) || !identical(dim(x), dim(y)) ||
        !identical(dimnames(x), dimnames(y))) return(paste(path, "schema"))
    if (is.list(x)) {
      if (length(x) != length(y)) return(paste(path, "length"))
      return(unlist(lapply(seq_along(x), function(i) mismatch(x[[i]], y[[i]],
        paste0(path, "$", if (is.null(names(x))) i else names(x)[i]))), use.names = FALSE))
    }
    same <- if (typeof(x) == "double") length(waldo::compare(x, y, tolerance = tol)) == 0L else identical(x, y)
    if (same) character() else paste(path, "value")
  }
  bad <- mismatch(a, b, "root")
  expect_true(length(bad) == 0L, info = paste(c(info, head(bad, 12L)), collapse = "\n"))
}

dc_wrap <- c("construct_ipw_link_view", "validate_ipw_link_view", "cp04_public_link_views",
  "ipw_owner_link_support_context", "cp04_descriptor_context", "ipw_pipeline_compact_basis_facts",
  "or_baseline_complete_facts", "or_baseline_descriptor_intake", "ipw_owner_make_cohort_descriptor",
  "ipw_owner_make_pair_descriptor", "fit_or_baseline", "fit_dr_baseline", "or_baseline_complete_state",
  "dr_baseline_treatment_batch")

d20 <- dc_panel(41L, 20L)
dc_trace <- new.env(parent = emptyenv())

# ---------------------------------------------------------------------------
# 1. Descriptor/ordinary equality on the 190-comparison panel: all original
#    records, IDs and auxiliary outcomes; the exact global 190-block objective,
#    QR/rank basis, nonconstant PS, global beta, full-F IF/V/SE; the original
#    forest/metric; constructor/validation counts showing the avoided views.
# ---------------------------------------------------------------------------

test_that("descriptor context equals the link-view route on OR and DR with and without X and avoids every heavy view", {
  for (est in c("reg", "dr")) for (xf in list(NULL, ~x)) {
    info <- paste(est, if (is.null(xf)) "noX" else deparse(xf))
    trace <- new.env(parent = emptyenv())
    e <- dc_capture(dc_counted(dc_wrap), trace)
    g <- dc_ordinary(dc_wrap)
    res_d <- dc_run(d20, est = est, xformla = xf, ns = e)
    res_o <- dc_run(d20, est = est, xformla = xf, ns = g)
    M <- nrow(res_d$metadata$catalog$rows)
    Q <- length(res_d$nuisance_registry$cp04_aux$descriptors)
    expect_identical(M, 190L, info = info)
    expect_true(Q >= 1L, info = info)

    # Route records: one descriptor context, zero heavy views; the ordinary route built M + Q.
    cw <- res_d$nuisance_registry$cp04_construction
    co <- res_o$nuisance_registry$cp04_construction
    expect_identical(cw$route, "descriptor_context", info = info)
    expect_identical(cw$heavy_link_views_constructed, 0L, info = info)
    expect_identical(cw$descriptor_rows, as.integer(M + Q), info = info)
    expect_identical(c(cw$catalog_rows, cw$auxiliary_rows), as.integer(c(M, Q)), info = info)
    expect_true(cw$date_facts_checked >= 20, info = info)
    expect_true(cw$context_entries > 0 && cw$context_entries <= cw$context_entries_bound, info = info)
    expect_identical(co$route, "link_views", info = info)
    expect_identical(co$heavy_link_views_constructed, as.integer(M + Q), info = info)

    # Actual constructor/validation/seal counts of both clones.
    expect_identical(dc_n(e, "construct_ipw_link_view"), 0, info = info)
    expect_identical(dc_n(e, "validate_ipw_link_view"), 0, info = info)
    expect_identical(dc_n(e, "cp04_public_link_views"), 0, info = info)
    expect_identical(dc_n(e, "ipw_owner_link_support_context"), 0, info = info)
    expect_identical(dc_n(e, "seal:link_view"), 0, info = info)
    expect_identical(dc_n(e, "cp04_descriptor_context"), 1, info = info)
    expect_identical(dc_n(e, "seal:cp04_descriptor_context"), 1, info = info)
    expect_identical(dc_n(e, "ipw_pipeline_compact_basis_facts"), 1, info = info)
    expect_identical(dc_n(e, "or_baseline_descriptor_intake"), 1, info = info)
    expect_identical(dc_n(e, "or_baseline_complete_facts"), 1, info = info)
    expect_identical(dc_n(e, "or_baseline_complete_state"), 1, info = info)
    # One cohort descriptor in the registry and its re-establishment in the fit;
    # one pair descriptor per original row; the view route built one of each per view.
    expect_identical(dc_n(e, "ipw_owner_make_cohort_descriptor"), 2, info = info)
    MQ <- as.double(M + Q)
    expect_identical(dc_n(e, "ipw_owner_make_pair_descriptor"), MQ, info = info)
    expect_identical(dc_n(g, "ipw_owner_make_cohort_descriptor"), MQ, info = info)
    expect_identical(dc_n(g, "construct_ipw_link_view"), MQ, info = info)
    expect_true(dc_n(g, "validate_ipw_link_view") >= M + Q, info = info)
    expect_identical(dc_n(g, "cp04_public_link_views"), 1, info = info)
    expect_identical(dc_n(g, "seal:link_view"), MQ, info = info)
    expect_identical(dc_n(g, "cp04_descriptor_context"), 1, info = info)
    expect_identical(dc_n(g, "ipw_pipeline_compact_basis_facts"), 0, info = info)
    expect_identical(dc_n(g, "or_baseline_descriptor_intake"), 1, info = info)
    expect_identical(dc_n(g, "seal:cp04_descriptor_context"), 0, info = info)
    expect_identical(dc_n(e, "fit_or_baseline") + dc_n(e, "fit_dr_baseline"), 1, info = info)

    # The sealed context itself: genuine bindings and every original descriptor.
    ctx <- trace$ctx
    expect_true(inherits(ctx, "cdid_cp04_descriptor_context"), info = info)
    expect_true(api$ipw_owner_seal_ok(ctx, "cp04_descriptor_context"), info = info)
    expect_identical(ctx$frame_ref, res_d$frame_owner$ref, info = info)
    expect_identical(ctx$panel_ref, res_d$frame_owner$panel_ref, info = info)
    expect_identical(ctx$reader_ref, res_d$nuisance_registry$cp04_reader$ref, info = info)
    expect_identical(ctx$metadata_ref, res_d$metadata$ref, info = info)
    expect_identical(ctx$method, if (est == "dr") "dr" else "or", info = info)
    expect_identical(length(ctx$descriptors), as.integer(M + Q), info = info)
    expect_identical(vapply(ctx$descriptors, `[[`, character(1), "link_id"),
      c(paste0("cp04:catalog:", seq_len(M)), paste0("cp04:auxiliary:", seq_len(Q))), info = info)
    expect_identical(ctx$K, Inf, info = info)
    expect_identical(ctx$g, 3, info = info)
    expect_identical(unname(ctx$support$by_cohort_pair_count), unname(ctx$cohorts$group_a_count), info = info)
    expect_identical(ctx$support$n_observed_risk_units, sum(ctx$cohorts$group_a_count), info = info)
    expect_identical(ctx$support$n_observed_risk_units, 41L, info = info)
    expect_true(isTRUE(ctx$outcome_recipe$intercept), info = info)
    if (est == "dr") expect_true(inherits(ctx$treatment_recipe, "cdid_dr_treatment_recipe"), info = info) else
      expect_null(ctx$treatment_recipe, info = info)

    # Backends: the same complete-support fit, all original records and IDs,
    # same retained payload (no record difference at all).
    b_d <- res_d$nuisance_registry$cp04
    b_o <- res_o$nuisance_registry$cp04
    expect_false(is.null(b_d$complete_support), info = info)
    expect_false(is.null(b_o$complete_support), info = info)
    dc_route_match(dc_backend_view(b_d), dc_backend_view(b_o), info = info)
    expect_identical(b_d$accounting$owned_entries, b_o$accounting$owned_entries, info = info)
    expect_identical(b_d$accounting$external_work_ops,
      b_o$accounting$external_work_ops + ctx$accounting$work_ops_bound, info = info)
    expect_identical(b_d$accounting$borrowed_entries, b_o$accounting$borrowed_entries, info = info)
    expect_identical(vapply(b_d$links, `[[`, character(1), "link_id"),
      c(paste0("cp04:catalog:", seq_len(M)), paste0("cp04:auxiliary:", seq_len(Q))), info = info)
    expect_identical(vapply(b_d$links, `[[`, character(1), "source_functional_key"),
      vapply(b_o$links, `[[`, character(1), "source_functional_key"), info = info)
    expect_identical(length(b_d$regressions), 190L, info = info)
    expect_identical(b_d$accounting$n_local_intercepts, 190L, info = info)
    expect_identical(length(b_d$supports), 380L, info = info)
    expect_identical(length(b_d$links), as.integer(M + Q), info = info)
    expect_true(all(is.finite(vapply(b_d$links, `[[`, numeric(1), "mu"))), info = info)
    expect_identical(b_d$accounting$complete_support$n_support_row_visits,
      b_o$accounting$complete_support$n_support_row_visits, info = info)
    # Exact global objective and the one stable QR/shared beta basis.
    expect_equal(b_d$shared$objective, b_o$shared$objective, tolerance = 1e-12, info = info)
    expect_equal(b_d$shared$rss, b_o$shared$rss, tolerance = 1e-12, info = info)
    expect_identical(b_d$shared$keep, b_o$shared$keep, info = info)
    expect_identical(b_d$shared$basis$kept_columns, b_o$shared$basis$kept_columns, info = info)
    expect_equal(b_d$shared$beta, b_o$shared$beta, tolerance = 1e-12, info = info)
    if (!is.null(xf)) {
      expect_identical(b_d$shared$p, 1L, info = info)
      expect_true(all(is.finite(b_d$shared$beta)), info = info)
    } else expect_identical(b_d$shared$p, 0L, info = info)
    if (est == "dr") {
      expect_identical(length(b_d$treatments), 190L, info = info)
      expect_identical(b_d$accounting$n_treatment_fits, 190L, info = info)
      acc <- b_d$accounting$complete_support
      expect_identical(c(acc$n_treatment_logit_fits, acc$n_treatment_logit_reuses, acc$n_treatment_models),
        c(1L, 189L, 190L), info = info)
      expect_identical(acc$n_treatment_summary_batches, 1L, info = info)
      expect_identical(dc_n(e, "dr_baseline_treatment_batch"), 1, info = info)
      expect_true(all(vapply(b_d$treatment_specs, `[[`, integer(1), "n_observers") == 41L), info = info)
      expect_identical(vapply(b_d$treatment_specs, `[[`, integer(1), "n_observers"),
        vapply(b_o$treatment_specs, `[[`, integer(1), "n_observers"), info = info)
      coef <- b_d$treatments[[1L]]$coefficients
      expect_identical(length(coef), if (is.null(xf)) 1L else 2L, info = info)
      expect_true(all(is.finite(coef)), info = info)
      # Nonconstant propensity under X: the fitted logit slope is not zero.
      if (!is.null(xf)) expect_true(abs(coef[[2L]]) > 0, info = info)
      dc_route_match(dc_strip(b_d$treatments), dc_strip(b_o$treatments), info = paste(info, "treatments"))
    }

    # Registry records and auxiliary outcomes.
    expect_identical(dc_strip(res_d$nuisance_registry$cp04_points), dc_strip(res_o$nuisance_registry$cp04_points), info = info)
    expect_identical(res_d$nuisance_registry$cp04_aux_index, res_o$nuisance_registry$cp04_aux_index, info = info)
    expect_identical(res_d$nuisance_registry$cp04_aux_index, as.integer(M + seq_len(Q)), info = info)
    expect_identical(res_d$nuisance_registry$cp04_aux_errors, res_o$nuisance_registry$cp04_aux_errors, info = info)
    expect_true(all(vapply(res_d$nuisance_registry$cp04_aux_errors, is.null, logical(1))), info = info)
    expect_identical(dc_strip(res_d$point_links_all), dc_strip(res_o$point_links_all), info = info)
    expect_identical(dc_strip(res_d$auxiliary_reporting), dc_strip(res_o$auxiliary_reporting), info = info)
    expect_identical(res_d$nuisance_registry$nuisance_keys, res_o$nuisance_registry$nuisance_keys, info = info)

    # Public ATT, full original-F IF, V and SE at the established tolerance.
    expect_identical(res_d$reporting$target_meta, res_o$reporting$target_meta, info = info)
    expect_identical(res_d$reporting$available, res_o$reporting$available, info = info)
    expect_true(any(res_d$reporting$available), info = info)
    for (nm in c("estimates", "se", "analytical_se", "if_full", "vcov"))
      expect_equal(res_d$reporting[[nm]], res_o$reporting[[nm]], tolerance = 1e-12, info = paste(info, nm))
    expect_identical(nrow(res_d$reporting$if_full), 41L, info = info)
    av <- res_d$reporting$available
    expect_true(all(is.finite(res_d$reporting$if_full[, av, drop = FALSE])), info = info)
    for (ag in c("simple", "group", "dynamic", "calendar")) {
      expect_equal(res_d$aggregations[[ag]]$overall$estimate, res_o$aggregations[[ag]]$overall$estimate, tolerance = 1e-12, info = info)
      expect_equal(res_d$aggregations[[ag]]$overall$se, res_o$aggregations[[ag]]$overall$se, tolerance = 1e-12, info = info)
      expect_equal(res_d$aggregations[[ag]]$estimates, res_o$aggregations[[ag]]$estimates, tolerance = 1e-12, info = info)
      expect_equal(res_d$aggregations[[ag]]$se, res_o$aggregations[[ag]]$se, tolerance = 1e-12, info = info)
    }

    # The accepted forest governs the post-fit representation only: 190 original
    # links, 19 retained columns, every original record fitted; same E, metric,
    # ridge handling and GMM provenance as the view route.
    cert_d <- res_d$structural_certificate
    cert_o <- res_o$structural_certificate
    expect_false(is.null(cert_d), info = info)
    expect_identical(cert_d$n_links, 190L, info = info)
    expect_identical(length(cert_d$retained_index), 19L, info = info)
    expect_identical(length(cert_d$dropped_index), 171L, info = info)
    for (f in c("n_links", "retained_index", "dropped_index", "E", "link_kind"))
      expect_identical(dc_strip(cert_d[[f]]), dc_strip(cert_o[[f]]), info = paste(info, "cert", f))
    expect_identical(dc_strip(res_d$covariance_handling), dc_strip(res_o$covariance_handling), info = info)
    expect_identical(res_d$gmm_fit$provenance$reduction$kept_indices, res_o$gmm_fit$provenance$reduction$kept_indices, info = info)
    dc_route_match(dc_strip(res_d$gmm_fit$targets), dc_strip(res_o$gmm_fit$targets), info = paste(info, "gmm targets"))
    diag_d <- api$cp04_public_diagnostics(res_d, res_d$controls$pipeline_control)
    diag_o <- api$cp04_public_diagnostics(res_o, res_o$controls$pipeline_control)
    # R01's preparation work is a real additional cumulative charge; every
    # other diagnostic field remains identical to the original view route.
    for (nm in c("external_work_ops", "work_ops_bound"))
      diag_o$accounting[[nm]] <- diag_o$accounting[[nm]] + ctx$accounting$work_ops_bound
    expect_identical(dc_strip(diag_d), dc_strip(diag_o), info = info)

    dc_trace[[paste0(est, if (is.null(xf)) "_noX" else "_x")]] <- list(ctx = ctx, args = trace$args, res = res_d,
      backend = b_d)
  }
})

# ---------------------------------------------------------------------------
# 2. Structural negatives keep the original link-view route (predicate called,
#    nothing sealed, every heavy view built as before).
# ---------------------------------------------------------------------------

test_that("structural negatives keep the ordinary link-view route with its original constructor counts", {
  d5 <- dc_panel(41L, 5L, seed = 7L)
  cases <- list(
    nonunit_weights = list(d = d5, args = list(weightsname = "w")),
    two_targets = list(d = dc_panel(41L, 5L, seed = 7L, second = TRUE), args = list()),
    extra_treated_with_target_subset = list(d = dc_panel(41L, 5L, seed = 7L, second = TRUE),
      args = list(target_cohorts = 3)),
    pt_all = list(d = d5, args = list(pt_assumption = "all", pt_start = 1)),
    missing_obs = list(d = d5[!(d5$id == d5$id[1L] & d5$period == 3L), ], args = list())
  )
  for (nm in names(cases)) for (est in c("reg", "dr")) {
    info <- paste(nm, est)
    cs <- cases[[nm]]
    e <- dc_counted(dc_wrap)
    res <- do.call(dc_run, c(list(d = cs$d, est = est, ns = e), cs$args))
    M <- nrow(res$metadata$catalog$rows)
    Q <- length(res$nuisance_registry$cp04_aux$descriptors)
    expect_true(M >= 1L, info = info)
    cw <- res$nuisance_registry$cp04_construction
    expect_identical(cw$route, "link_views", info = info)
    expect_identical(cw$heavy_link_views_constructed, as.integer(M + sum(!is.na(res$nuisance_registry$cp04_aux_index))), info = info)
    expect_identical(dc_n(e, "cp04_descriptor_context"), 1, info = info)
    expect_identical(dc_n(e, "seal:cp04_descriptor_context"), 0, info = info)
    expect_identical(dc_n(e, "or_baseline_descriptor_intake"), 1, info = info)
    expect_identical(dc_n(e, "construct_ipw_link_view"), as.double(M + Q), info = info)
    expect_identical(dc_n(e, "cp04_public_link_views"), 1, info = info)
    expect_true(dc_n(e, "validate_ipw_link_view") >= M + Q, info = info)
    expect_identical(dc_n(e, "fit_or_baseline") + dc_n(e, "fit_dr_baseline"), 1, info = info)
    expect_true(any(is.finite(res$reporting$estimates)), info = info)
    if (nm == "missing_obs" && est == "reg") dc_trace$foreign <- res
  }
})

# ---------------------------------------------------------------------------
# 3. Direct predicate negatives on the captured actual arguments, and the
#    exact adjacent caps admitted/refused before any posting read or allocation.
# ---------------------------------------------------------------------------

test_that("the descriptor predicate declines unsupported requests and admits exactly at its adjacent caps before any read", {
  base <- dc_trace$reg_noX
  expect_false(is.null(base), info = "section 1 must have captured the OR context")
  args <- base$args
  ctx <- base$ctx
  call <- function(..., ns = api) {
    supplied <- args
    replacement <- list(...)
    supplied[names(replacement)] <- replacement
    do.call(ns$cp04_descriptor_context, supplied)
  }

  # The genuine predicate reproduces the sealed context's facts on the same arguments.
  again <- call()
  expect_true(inherits(again, "cdid_cp04_descriptor_context"))
  drop_acc <- function(x) dc_strip(unclass(x)[setdiff(names(x), c("ref", "accounting"))])
  expect_identical(drop_acc(again), drop_acc(ctx))
  expect_identical(again$accounting, ctx$accounting)

  # Standalone callers without request facts, PT-all, numerical transfer, other
  # methods, a no-intercept recipe, other pools/targets on any row: original route.
  expect_null(call(pt_assumption = NULL))
  expect_null(call(pt_assumption = "all"))
  expect_null(call(run_numerical_transfer = TRUE))
  expect_null(call(or_method = NULL, dr_method = NULL))
  expect_null(call(or_method = "corollary1"))
  cols <- colnames(args$frame$Z_a)
  expect_null(call(outcome_recipe = api$or_pooled_outcome_recipe(cols, FALSE, cols)))
  aux_pool <- args$aux; aux_pool$descriptors[[1L]]$K <- c(4, Inf)
  expect_null(call(aux = aux_pool))
  aux_target <- args$aux; aux_target$descriptors[[1L]]$g <- 4
  expect_null(call(aux = aux_target))
  spec_target <- args$specs; spec_target[[1L]]$g <- 4
  expect_null(call(specs = spec_target))
  spec_key <- args$specs; spec_key[[1L]]$functional_key <- "custom:not-a-catalog-key"
  expect_null(call(specs = spec_key))
  # Owners of another run are not these sealed bindings.
  foreign <- dc_trace$foreign
  expect_false(is.null(foreign))
  expect_null(call(reader = foreign$nuisance_registry$cp04_reader))
  expect_null(call(frame = foreign$frame_owner))

  # Adjacent caps: admitted exactly at the reported totals, refused one entry
  # below before the facts, any descriptor or any posting read.
  acc <- ctx$accounting
  edges <- list(
    list(field = "max_retained_entries", val = acc$retained_total),
    list(field = "max_matrix_entries", val = acc$retained_total + acc$scratch_entries),
    list(field = "max_work_ops", val = acc$work_total),
    list(field = "max_output_entries", val = acc$entries_bound))
  probe <- c("ipw_pipeline_compact_basis_facts", "ipw_owner_make_pair_descriptor", "ipw_owner_make_cohort_descriptor",
    "ipw_owner_date_complete", "ipw_owner_posting_date_complete", "or_public_default_recipe")
  for (edge in edges) {
    tight <- args
    tight$ctrl$owner_control[[edge$field]] <- edge$val
    ok <- do.call(api$cp04_descriptor_context, tight)
    expect_true(inherits(ok, "cdid_cp04_descriptor_context"), info = edge$field)
    expect_identical(drop_acc(ok), drop_acc(ctx), info = edge$field)
    tight$ctrl$owner_control[[edge$field]] <- edge$val - 1
    e <- dc_counted(probe)
    expect_error(do.call(get("cp04_descriptor_context", envir = e), tight), class = "cdid_ipw_owner_bound_error",
      info = edge$field)
    for (nm in probe) expect_identical(dc_n(e, nm), 0, info = paste(edge$field, nm))
    expect_identical(dc_n(e, "seal:cp04_descriptor_context"), 0, info = edge$field)
  }
})

# ---------------------------------------------------------------------------
# 4. Fit-level binding, estimator, recipe and mutation negatives; immediate
#    source-mutation detection of the sealed owners is preserved.
# ---------------------------------------------------------------------------

test_that("the fit intake refuses foreign, mutated, unsealed or mismatched contexts and keeps owner mutation detection", {
  base <- dc_trace$reg_noX
  dr_base <- dc_trace$dr_noX
  expect_false(is.null(base) || is.null(dr_base))
  ctx <- base$ctx
  res <- base$res
  reader <- res$nuisance_registry$cp04_reader
  pc <- res$controls$pipeline_control
  own <- pc$owner_control
  tol <- api$or_public_rank_tol()

  # The sealed context in the view position yields the registry's fitted backend.
  direct <- api$fit_or_baseline(reader, ctx, own, tol, outcome_recipe = ctx$outcome_recipe)
  expect_false(is.null(direct$complete_support))
  dc_route_match(dc_backend_view(direct)$fields, dc_backend_view(base$backend)$fields, info = "direct intake")

  # A modified copy of the sealed context is not the sealed snapshot.
  bad <- ctx; bad$descriptors[[1L]]$link_id <- "cp04:catalog:0"
  expect_error(api$fit_or_baseline(reader, bad, own, tol, outcome_recipe = ctx$outcome_recipe),
    class = "cdid_or_baseline_integrity_error")
  bad2 <- ctx; bad2$support$n_observed_risk_units <- 40L
  expect_error(api$fit_or_baseline(reader, bad2, own, tol, outcome_recipe = ctx$outcome_recipe),
    class = "cdid_or_baseline_integrity_error")
  # An unsealed classed copy, a label without bindings, is refused.
  fake <- unclass(ctx); fake$ref <- NULL
  fake <- structure(fake, class = "cdid_cp04_descriptor_context")
  expect_error(api$fit_or_baseline(reader, fake, own, tol, outcome_recipe = ctx$outcome_recipe),
    class = "cdid_or_baseline_integrity_error")
  # Another sealed reader (other owners) is refused at the intake boundary.
  foreign <- dc_trace$foreign
  expect_error(api$fit_or_baseline(foreign$nuisance_registry$cp04_reader, ctx, own, tol, outcome_recipe = ctx$outcome_recipe),
    class = "cdid_or_baseline_integrity_error")
  # The OR context does not serve the DR fit, nor the DR context the OR fit.
  expect_error(api$fit_dr_baseline(reader, ctx, own, glm_control = pc$glm_control, rank_tol = tol,
    outcome_recipe = ctx$outcome_recipe), class = "cdid_dr_baseline_integrity_error")
  ctx_dr <- dr_base$ctx
  dr_reader <- dr_base$res$nuisance_registry$cp04_reader
  expect_error(api$fit_or_baseline(dr_reader, ctx_dr, own, tol, outcome_recipe = ctx_dr$outcome_recipe),
    class = "cdid_or_baseline_integrity_error")
  # A different valid recipe of the same columns is not the context's checked recipe.
  cols <- ctx$outcome_recipe$columns
  other <- api$or_pooled_outcome_recipe(cols, FALSE, cols)
  expect_error(api$fit_or_baseline(reader, ctx, own, tol, outcome_recipe = other),
    class = "cdid_or_baseline_integrity_error")
  # Only the baseline consumer reads a context.
  expect_error(api$or_public_fit("corollary1", reader, ctx, pc, 0, 0), class = "cdid_cp04_integrity_error")
  expect_error(api$dr_public_fit("corollary1", dr_reader, ctx_dr, pc, 0, 0), class = "cdid_cp04_integrity_error")
  # The DR context in the DR fit reproduces the registry's DR backend.
  direct_dr <- api$fit_dr_baseline(dr_reader, ctx_dr, own, glm_control = pc$glm_control, rank_tol = tol,
    outcome_recipe = ctx_dr$outcome_recipe)
  expect_false(is.null(direct_dr$complete_support))
  dc_route_match(dc_backend_view(direct_dr)$fields, dc_backend_view(dr_base$backend)$fields, info = "direct DR intake")
  # Mutated sealed owners are still refused immediately, before any descriptor is read.
  mutant <- reader; mutant$state$w[1L] <- 2
  expect_error(api$fit_or_baseline(mutant, ctx, own, tol, outcome_recipe = ctx$outcome_recipe),
    class = "cdid_cp04_integrity_error")
})
