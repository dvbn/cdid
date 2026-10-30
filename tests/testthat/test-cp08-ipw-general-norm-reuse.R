# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 missing-panel IPW: one point-link normalization per general-route
# invocation. `run_chained_ipw_pipeline` seals its single
# `preflight_ipw_point_if_view(..., return_norm = TRUE)` with every validated
# object into an invocation-local record (`ipw_point_if_general_plan`); the two
# covariance passes and the two GMM projection passes replay through it
# (`stream_ipw_point_covariance_general`, `replay_ipw_point_if_general`) and
# re-verify it at every boundary and after every callback.
#
# Expected values come from (a) the UNCHANGED standalone wrappers of the same
# installed build (`replay_ipw_point_if_view`, `stream_ipw_point_covariance_view`,
# each normalizing per call exactly as before this change), (b) a forced
# old-route copy of the pipeline whose record consumers are those standalone
# wrappers, (c) the four reported IF components of the standalone replay and
# (d) frozen copies of the previous compact-storage visitor bodies. The record
# never supplies its own truth; floating parity is bitwise where the arithmetic
# path is identical and 1e-12 for whole-structure comparisons; identities,
# orders, shapes, counts and refusal classes compare exactly.
#
# STATUS: UNEXECUTED at authoring time (no R, build, install, fit or test was
# run in the authoring job). The counting namespaces below are disclosed
# test-local clones that re-bind internal functions for counting; they are not
# the installed authority and every numerical value is taken from genuine
# installed functions. The root sequences the independent installed gate.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Test-local helpers (clones of the accepted CP08 suites' fixtures)
# ---------------------------------------------------------------------------

# Re-bound copy of the namespace (as xb_ns / ru_ns): native resolution stays
# genuine; `edits` replace internal bindings; `wrap` names are counted.
gn_ns <- function(edits = list(), wrap = character()) {
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
  for (nm in names(edits)) assign(nm, edits[[nm]], envir = e)
  make_wrap <- function(nm, inner) {
    force(nm); force(inner)
    function(...) {
      assign(nm, get(nm, envir = counts, inherits = FALSE) + 1L, envir = counts)
      inner(...)
    }
  }
  for (nm in wrap) {
    if (!exists(nm, envir = e, inherits = FALSE)) stop("No internal function `", nm, "` to count.", call. = FALSE)
    assign(nm, 0L, envir = counts)
    assign(nm, make_wrap(nm, get(nm, envir = e, inherits = FALSE)), envir = e)
  }
  e$.counts <- counts
  e
}
gn_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0L)
gn_reset <- function(e) for (nm in ls(e$.counts, all.names = TRUE)) assign(nm, 0L, envir = e$.counts)
gn_strip_env <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); y <- lapply(x, gn_strip_env); attributes(y) <- a; return(y) }
  x
}
# Strict schema/leaf comparison, doubles at `tol`, everything else exactly.
gn_expect_match <- function(a, b, tol = 1e-12, info = NULL) {
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
gn_collect <- function(fn, ...) {
  out <- list()
  fn(..., callback = function(b) { force(b); out[[length(out) + 1L]] <<- b })
  out
}
# Rebind one locked reference binding; returns the restorer.
gn_rebind <- function(env, nm, value) {
  old <- get(nm, envir = env, inherits = FALSE)
  unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
  function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
}
gn_outcome <- function(expr) tryCatch({ force(expr); "ok" },
  error = function(e) list(class = class(e), message = conditionMessage(e)))

# ----- Panels: the accepted 41 x 6 design with deterministic missingness ----
# "balanced": every row (the general route is forced by declining the batch
# eligibility). "half": even units keep alternating dates only (half their
# rows). "mixed": scattered one-in-seven missing rows. "rotation": three
# rotation groups cutting across cohorts, each missing two dates. `extra`
# moves unit positions to a second treated cohort 6 (two fitted nuisances).
gn_panel <- function(N = 41L, T = 6L, pattern = "balanced", extra = NULL) {
  u <- seq_len(N)
  x1 <- (u %% 5L) / 4 - 0.5
  x2 <- ((3L * u) %% 7L) / 6 - 0.5
  g <- ifelse(((u %% 3L) == 0L) | (x1 > 0 & (u %% 2L) == 0L), 4, Inf)
  if (!is.null(extra)) g[extra] <- 6
  d <- expand.grid(t = seq_len(T), id = u, KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$x1 <- x1[d$id]
  d$x2 <- x2[d$id]
  d$y <- 0.3 * d$id + 0.2 * d$t + ifelse(d$t >= d$g, 1, 0) + 0.5 * d$x1 * d$t - 0.25 * d$x2 +
    ((7 * d$id + 3 * d$t) %% 5) / 10
  d$w_const <- 1
  keep <- rep(TRUE, nrow(d))
  if (identical(pattern, "half")) keep <- !((d$id %% 2L == 0L) & (((d$t + d$id %/% 2L) %% 2L) == 0L))
  if (identical(pattern, "mixed")) keep <- ((3L * d$id + 5L * d$t) %% 7L) != 0L
  if (identical(pattern, "rotation")) {
    r <- (d$id %/% 3L) %% 3L
    keep <- !((r == 1L & d$t %in% c(1L, 6L)) | (r == 2L & d$t %in% c(2L, 5L)))
  }
  d <- d[keep, ]
  rownames(d) <- NULL
  d
}
gn_formula <- function(p) switch(as.character(p), "1" = ~1, "2" = ~x1, "3" = ~x1 + x2, stop("p must be 1, 2 or 3"))
gn_ctrl <- function(d, p, chunk = 7L, weighting = "optimal") {
  gs <- unique(d$g)
  n_frame <- length(unique(d$id))
  dims <- api$chained_public_dimensions_from_counts(n_frame = n_frame, n_panel_units = n_frame,
    n_rows = nrow(d), n_dates = length(unique(d$t)), n_cohorts = length(gs), n_targets = sum(is.finite(gs)),
    n_covariates = p, biters = 99, n_clusters = n_frame)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(list(chunk_rows = chunk),
    gmm_weighting = weighting), gmm_weighting = weighting)$control
}
gn_run <- function(ns, d, ctrl, p, ...) ns$run_chained_ipw_pipeline(
  data = d, idname = "id", tname = "t", yname = "y", gname = "g", control = ctrl, xformla = gn_formula(p),
  control_group = "nevertreated", base_period = "varying", weighting = "optimal", run_numerical_transfer = FALSE,
  root_witness = NULL, structural_reduction_stage = "outcome_free", ...)

# The general route of the same build: the batch eligibility declines.
gn_general_edit <- list(ipw_pipeline_batch_eligibility = function(...) NULL)
gn_counts <- c("normalize_ipw_point_links_view", "preflight_ipw_point_if_view", "ipw_point_if_general_plan",
  "ipw_point_if_general_verify", "replay_ipw_point_if_general", "stream_ipw_point_covariance_general",
  "replay_ipw_point_if_view", "stream_ipw_point_covariance_view")
gn_new <- function(wrap = gn_counts) gn_ns(gn_general_edit, wrap)
# Forced OLD route of the same build: the record is still built (its
# accounting is part of the pipeline), but both consumers are the unchanged
# standalone wrappers on the record's own links, each normalizing per call as
# every entry did before this change (one covariance-entry normalization, two
# covariance replays, two projection replays).
gn_old <- function(wrap = gn_counts) {
  e <- gn_ns(gn_general_edit, wrap)
  e$stream_ipw_point_covariance_general <- function(plan, control, ...) {
    e$stream_ipw_point_covariance_view(point_links = plan$point_links, control = control)
  }
  e$replay_ipw_point_if_general <- function(plan, control, callback, callback_schema = TRUE, ...) {
    e$replay_ipw_point_if_view(point_links = plan$point_links, control = control, callback = callback,
      callback_schema = callback_schema)
  }
  e
}
gn_fields <- c("support", "basis", "aggregations", "gmm_fit", "structural_certificate", "design_all", "design_used",
  "covariance_stream", "covariance_handling", "full_frame_transport", "auxiliary_reporting")
gn_report <- c("estimates", "se", "analytical_se", "z_stat", "p_value", "ci_lower", "ci_upper", "if_full", "vcov", "available")
# Same fit or the same typed refusal; the whole structure at 1e-12, the
# covariance stream, moment estimates, reported cells and full-F IF bitwise.
gn_expect_same_result <- function(a, b, info = NULL) {
  if (inherits(a, "error") || inherits(b, "error")) {
    expect_true(inherits(a, "error") && inherits(b, "error"), info = paste(info, "one route failed, the other did not"))
    if (inherits(a, "error") && inherits(b, "error")) {
      expect_identical(class(a), class(b), info = info)
      expect_identical(conditionMessage(a), conditionMessage(b), info = info)
    }
    return(invisible(FALSE))
  }
  for (f in gn_fields) gn_expect_match(gn_strip_env(a[[f]]), gn_strip_env(b[[f]]), info = paste(info, f))
  for (nm in c("mhat", "V", "mean_phi", "n_a_units")) {
    expect_identical(a$covariance_stream[[nm]], b$covariance_stream[[nm]], info = paste(info, "covariance_stream", nm))
  }
  for (nm in gn_report) expect_identical(a$reporting[[nm]], b$reporting[[nm]], info = paste(info, "reporting", nm))
  expect_identical(a$reporting$target_meta, b$reporting$target_meta, info = info)
  expect_identical(gn_strip_env(a$resource_accounting$point_producer), gn_strip_env(b$resource_accounting$point_producer), info = info)
  expect_identical(a$resource_accounting$simultaneous_entries, b$resource_accounting$simultaneous_entries, info = info)
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches, info = info)
  expect_identical(a$gmm_fit$if_projection$max_batch_rows_seen, b$gmm_fit$if_projection$max_batch_rows_seen, info = info)
  expect_identical(names(a$point_links_used), names(b$point_links_used), info = info)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 1. Forced general replay on balanced, half, mixed, rotation and mixed-cohort
#    panels: the record route equals the forced old route (same fit or same
#    refusal), with one normalization across the covariance and projection
#    passes and one verification per pass boundary and per chunk callback.
# ---------------------------------------------------------------------------

test_that("general record route equals the forced old route on forced-general, half, mixed, rotation and mixed-cohort panels", {
  cases <- list(
    list(pattern = "balanced", p = 2L, extra = NULL),
    list(pattern = "half", p = 1L, extra = NULL),
    list(pattern = "mixed", p = 2L, extra = NULL),
    list(pattern = "rotation", p = 1L, extra = NULL),
    list(pattern = "balanced", p = 1L, extra = c(5L, 11L, 23L, 29L, 35L))
  )
  n_fitted <- 0L
  for (cs in cases) for (chunk in c(7L, 16L)) {
    d <- gn_panel(pattern = cs$pattern, extra = cs$extra)
    ctrl <- gn_ctrl(d, cs$p, chunk = chunk)
    info <- paste(cs$pattern, "p =", cs$p, "chunk =", chunk, if (is.null(cs$extra)) "" else "two cohorts")
    e_new <- gn_new(); e_old <- gn_old()
    r_new <- tryCatch(gn_run(e_new, d, ctrl, cs$p), error = identity)
    r_old <- tryCatch(gn_run(e_old, d, ctrl, cs$p), error = identity)
    ok <- gn_expect_same_result(r_new, r_old, info = info)
    if (!isTRUE(ok)) next
    n_fitted <- n_fitted + 1L
    expect_null(r_new$resource_accounting$point_batch, info = info)
    pg <- r_new$resource_accounting$point_general
    expect_false(is.null(pg), info = info)
    expect_identical(pg$normalizations_performed, 1L, info = info)
    expect_identical(pg$passes_planned, 4L, info = info)
    expect_identical(pg$n_links, length(r_new$point_links_used), info = info)
    expect_gt(pg$additional_live_entries, 0)
    expect_gt(pg$verification_work_ops_bound, 0)
    expect_false(isTRUE(pg$full_n_by_L_retained))
    n_chunks <- as.integer(r_new$resource_accounting$point_producer$n_chunks)
    expect_identical(n_chunks, as.integer(ceiling(length(unique(d$id)) / min(chunk, length(unique(d$id))))), info = info)
    # Record route: one record, one covariance stream, four record replays, no
    # standalone covariance stream; the standalone replay serves only the
    # auxiliary reporting, exactly as on the old route minus its four passes.
    expect_identical(gn_n(e_new, "ipw_point_if_general_plan"), 1L, info = info)
    expect_identical(gn_n(e_new, "stream_ipw_point_covariance_general"), 1L, info = info)
    expect_identical(gn_n(e_new, "replay_ipw_point_if_general"), 4L, info = info)
    expect_identical(gn_n(e_new, "stream_ipw_point_covariance_view"), 0L, info = info)
    expect_identical(gn_n(e_old, "replay_ipw_point_if_view") - gn_n(e_new, "replay_ipw_point_if_view"), 4L, info = info)
    expect_identical(gn_n(e_old, "stream_ipw_point_covariance_view"), 1L, info = info)
    # Six normalizations (record, covariance entry, two covariance replays, two
    # projection replays) became one: the other five are gone.
    expect_identical(gn_n(e_old, "normalize_ipw_point_links_view") - gn_n(e_new, "normalize_ipw_point_links_view"), 5L, info = info)
    expect_identical(gn_n(e_old, "preflight_ipw_point_if_view") - gn_n(e_new, "preflight_ipw_point_if_view"), 5L, info = info)
    # Raw reservation verifies once; covariance entry verifies once. The fill
    # uses the original producer's entry and C post-callback checks, adds two
    # carrier checks per callback, then verifies the filled state once. Each
    # of the three slice passes verifies at entry and before/after C callbacks.
    # The covariance/fill carrier-only checks with verify_plan = FALSE add no
    # general-plan checks. Without raw admission, reservation and covariance
    # entry plus an optional sparse-state build precede four ordinary replays.
    expected_verifications <- if (isTRUE(pg$retained_if$admitted))
      1L + 1L + (1L + n_chunks) + 2L * n_chunks + 1L + 3L * (1L + 2L * n_chunks) else
      2L + as.integer(isTRUE(pg$replay_state$admitted)) + 4L * (1L + n_chunks)
    expect_identical(gn_n(e_new, "ipw_point_if_general_verify"),
      expected_verifications, info = info)
    # The forced-old clone still verifies at raw reservation (and at any
    # admitted sparse-state build), but replaces all four record consumers.
    expect_identical(gn_n(e_old, "ipw_point_if_general_verify"),
      1L + as.integer(isTRUE(r_old$resource_accounting$point_general$replay_state$admitted)), info = info)
    # The record did not survive into the result.
    found <- FALSE
    walk <- function(x, depth = 0L) {
      if (found || depth > 10L) return(invisible(NULL))
      if (is.list(x)) {
        if (identical(x$kind, "cdid_ipw_point_if_general_plan")) found <<- TRUE
        for (v in x) walk(v, depth + 1L)
      }
      invisible(NULL)
    }
    walk(r_new)
    expect_false(found, info = info)
    expect_true(all(is.finite(r_new$reporting$if_full[, which(r_new$reporting$available), drop = FALSE])), info = info)
  }
  expect_gte(n_fitted, 4L)
  # An unedited legacy-route counting clone takes the general route on a
  # missing panel (batch eligibility declines on incomplete dates) and agrees
  # bitwise with the forced-general counting clone. The installed namespace
  # itself selects the compact standalone raw-IPW backend on incomplete panels
  # (`R/chained_ipw_pipeline.R:4868-4899,5583-5586`) and agrees at 1e-12 on all
  # scientific outputs while recording raw-backend route metadata.
  d <- gn_panel(pattern = "half"); ctrl <- gn_ctrl(d, 1L, chunk = 7L)
  legacy_ns <- gn_ns()
  general_ns <- gn_new()
  r_ns <- tryCatch(gn_run(legacy_ns, d, ctrl, 1L), error = identity)
  r_new <- tryCatch(gn_run(general_ns, d, ctrl, 1L), error = identity)
  gn_expect_same_result(r_ns, r_new, info = "unedited legacy-route clone vs forced-general clone, half")
  if (!inherits(r_ns, "error")) {
    expect_null(r_ns$resource_accounting$point_batch)
    expect_identical(r_ns$resource_accounting$point_general$normalizations_performed, 1L)
  }
  r_api <- tryCatch(gn_run(api, d, ctrl, 1L), error = identity)
  expect_false(inherits(r_api, "error"))
  expect_false(inherits(r_new, "error"))
  expect_s3_class(r_api$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
  expect_null(r_new$nuisance_registry$raw_ipw)
  expect_null(r_api$resource_accounting$point_batch)
  expect_null(r_api$resource_accounting$point_general)
  expect_identical(r_api$resource_accounting$point_producer$chunk_rows, 7)
  expect_identical(r_api$resource_accounting$point_producer$retained_moment_if$status, "streamed")
  expect_identical(r_api$metadata$catalog$rows, r_new$metadata$catalog$rows)
  expect_identical(gn_strip_env(r_api$support), gn_strip_env(r_new$support))
  expect_identical(gn_strip_env(r_api$design_all), gn_strip_env(r_new$design_all))
  expect_identical(gn_strip_env(r_api$design_used), gn_strip_env(r_new$design_used))
  expect_identical(r_api$reporting$target_meta, r_new$reporting$target_meta)
  expect_identical(r_api$reporting$available, r_new$reporting$available)
  expect_identical(r_api$reporting$support, r_new$reporting$support)
  expect_identical(names(r_api$point_links_used), names(r_new$point_links_used))
  expect_identical(r_api$structural_certificate$retained_index, r_new$structural_certificate$retained_index)
  expect_identical(r_api$structural_certificate$dropped_index, r_new$structural_certificate$dropped_index)
  expect_identical(r_api$structural_certificate$classes, r_new$structural_certificate$classes)
  expect_identical(r_api$structural_certificate$E, r_new$structural_certificate$E)
  expect_identical(r_api$covariance_handling$handler_used, r_new$covariance_handling$handler_used)
  penalty_identity <- c("weighting_requested", "weighting_canonical", "handler_selected", "handler_invoked",
                        "condition_trigger", "ridge_relative", "target_availability")
  expect_identical(r_api$covariance_handling[penalty_identity], r_new$covariance_handling[penalty_identity])
  gn_expect_match(r_api$covariance_handling$lambda, r_new$covariance_handling$lambda,
    info = "raw vs legacy lambda, half")
  expect_identical(r_api$gmm_fit$if_projection$n_batches, r_new$gmm_fit$if_projection$n_batches)
  expect_identical(r_api$gmm_fit$if_projection$max_batch_rows_seen, r_new$gmm_fit$if_projection$max_batch_rows_seen)
  expect_identical(r_api$reporting$z_stat, r_api$reporting$estimates / r_api$reporting$se)
  expect_identical(r_new$reporting$z_stat, r_new$reporting$estimates / r_new$reporting$se)
  expect_identical(r_api$covariance_stream$n_a_units, r_new$covariance_stream$n_a_units)
  gn_expect_match(r_api$covariance_stream$mhat, r_new$covariance_stream$mhat, info = "raw vs legacy mhat, half")
  gn_expect_match(r_api$covariance_stream$V, r_new$covariance_stream$V, info = "raw vs legacy V, half")
  gn_expect_match(as.double(r_api$covariance_stream$mean_phi), as.double(r_new$covariance_stream$mean_phi),
    info = "raw vs legacy mean_phi, half")
  for (nm in gn_report) {
    gn_expect_match(r_api$reporting[[nm]], r_new$reporting[[nm]], info = paste("raw vs legacy reporting", nm, "half"))
  }
  gn_expect_match(r_api$reporting$Omega_a, r_new$reporting$Omega_a, info = "raw vs legacy Omega_a, half")
  gn_expect_match(gn_strip_env(r_api$full_frame_transport), gn_strip_env(r_new$full_frame_transport),
    info = "raw vs legacy full_frame_transport, half")
  gn_expect_match(gn_strip_env(r_api$basis), gn_strip_env(r_new$basis), info = "raw vs legacy basis, half")
  gn_expect_match(gn_strip_env(r_api$aggregations), gn_strip_env(r_new$aggregations),
    info = "raw vs legacy aggregations, half")
  gn_expect_match(r_api$gmm_fit$targets, r_new$gmm_fit$targets, info = "raw vs legacy gmm_fit targets, half")
  gn_expect_match(r_api$gmm_fit$if_projection$sandwich_cov, r_new$gmm_fit$if_projection$sandwich_cov,
    info = "raw vs legacy sandwich_cov, half")
  gn_expect_match(r_api$gmm_fit$if_projection$materialized_target_if,
    r_new$gmm_fit$if_projection$materialized_target_if, info = "raw vs legacy materialized_target_if, half")
  gn_expect_match(r_api$auxiliary_reporting$estimates, r_new$auxiliary_reporting$estimates,
    info = "raw vs legacy auxiliary estimates, half")
  expect_identical(r_api$auxiliary_reporting$selections, r_new$auxiliary_reporting$selections)
  expect_identical(r_api$auxiliary_reporting$descriptors, r_new$auxiliary_reporting$descriptors)
  expect_identical(r_api$auxiliary_reporting$descriptor_status, r_new$auxiliary_reporting$descriptor_status)
})

# ---------------------------------------------------------------------------
# 2. Private level: record replays and covariance equal the standalone wrappers
#    bitwise across chunk sizes; the four IF terms and their order are intact.
# ---------------------------------------------------------------------------

gn_fit_links <- function(pattern = "half", p = 1L, chunk = 7L, extra = NULL) {
  d <- gn_panel(pattern = pattern, extra = extra)
  ctrl <- gn_ctrl(d, p, chunk = chunk)
  res <- gn_run(gn_new(character()), d, ctrl, p)
  list(d = d, ctrl = ctrl, res = res, links = res$point_links_used, pc = ctrl$point_control)
}

test_that("record replays and covariance equal the standalone wrappers bitwise; all four IF terms are present in their order", {
  fx <- gn_fit_links("half", p = 1L)
  links <- fx$links
  expect_gte(length(links), 2L)
  for (chunk in c(5L, 7L, 41L)) {
    pc <- fx$pc; pc$max_chunk_rows <- chunk
    info <- paste("chunk", chunk)
    plan <- api$ipw_point_if_general_plan(links, pc)
    expect_identical(plan$kind, "cdid_ipw_point_if_general_plan")
    expect_true(api$ipw_owner_seal_ok(plan, "point_if_general_plan"))
    expect_identical(plan$general$normalizations_performed, 1L)
    expect_identical(plan$n_links, length(links))
    expect_identical(plan$link_ids, names(links))
    # The record's preflight is the general preflight of the same inputs.
    expect_identical(plan$preflight, api$preflight_ipw_point_if_view(links, pc), info = info)
    fitted_refs <- unique(lapply(links, function(lk) lk$fitted_nuisance_ref))
    expect_identical(plan$general$n_distinct_fitted_nuisances, length(fitted_refs), info = info)
    expect_identical(plan$general$n_distinct_fitted_nuisances, as.integer(plan$preflight$n_unique_fitted_nuisances))
    for (ell in seq_along(links)) {
      expect_true(identical(plan$fitted_refs[[plan$link_fitted_index[ell]]], links[[ell]]$fitted_nuisance_ref))
      expect_true(identical(plan$nuisance_refs[[plan$link_nuisance_index[ell]]], links[[ell]]$nuisance_view_ref))
    }
    gb <- gn_collect(api$replay_ipw_point_if_general, plan = plan, control = pc)
    sb <- gn_collect(api$replay_ipw_point_if_view, point_links = links, control = pc, callback_schema = TRUE)
    cb <- gn_collect(api$replay_ipw_point_if_view, point_links = links, control = pc, include_components = TRUE)
    expect_identical(length(gb), length(sb), info = info)
    expect_identical(length(gb), as.integer(plan$preflight$n_chunks), info = info)
    any_nonzero <- c(rate = FALSE, beta_score = FALSE, treated = FALSE, control = FALSE)
    for (i in seq_along(gb)) {
      expect_identical(gb[[i]]$phi, sb[[i]]$phi, info = paste(info, "phi", i))
      expect_identical(gb[[i]]$id, sb[[i]]$id, info = paste(info, "id", i))
      expect_identical(gb[[i]]$a_index, sb[[i]]$a_index, info = paste(info, "a_index", i))
      expect_identical(gb[[i]]$link_ids, sb[[i]]$link_ids, info = paste(info, "link_ids", i))
      expect_identical(colnames(gb[[i]]$phi), plan$link_ids, info = paste(info, "columns", i))
      expect_identical(names(gb[[i]]), c("phi", "id", "link_ids", "a_index"), info = paste(info, "schema", i))
      # The standalone replay forms every column as ((treated + control) +
      # beta_score) + rate; the record route reproduces that sum bitwise.
      comp <- cb[[i]]$components
      expect_identical(gb[[i]]$phi, ((comp$treated + comp$control) + comp$beta_score) + comp$rate,
        info = paste(info, "four-term order", i))
      expect_identical(gb[[i]]$a_index, cb[[i]]$a_index)
      for (nm in names(any_nonzero)) if (any(comp[[nm]] != 0)) any_nonzero[[nm]] <- TRUE
    }
    expect_true(all(any_nonzero), info = paste(info, "every IF term contributes somewhere"))
    cg <- api$stream_ipw_point_covariance_general(plan, pc)
    cv <- api$stream_ipw_point_covariance_view(links, pc)
    for (nm in c("covariance", "mean_phi", "estimates", "n_a_units")) {
      expect_identical(cg[[nm]], cv[[nm]], info = paste(info, "covariance", nm))
    }
    expect_identical(cg$replay$n_batches, cv$replay$n_batches, info = info)
    expect_identical(cg$replay$normalizations_performed, 0L)
    expect_identical(cg$replay$verifications, 1L + as.integer(plan$preflight$n_chunks))
    expect_identical(cg$replay$preflight, cv$replay$preflight, info = info)
    expect_equal(sqrt(diag(cg$covariance)), sqrt(diag(cv$covariance)), tolerance = 1e-12)
  }
  # Two fitted nuisances on one record (mixed cohorts): same equalities.
  fx2 <- tryCatch(gn_fit_links("balanced", p = 1L, extra = c(5L, 11L, 23L, 29L, 35L)), error = identity)
  if (!inherits(fx2, "error")) {
    plan2 <- api$ipw_point_if_general_plan(fx2$links, fx2$pc)
    expect_gt(plan2$general$n_distinct_fitted_nuisances, 1L)
    gb <- gn_collect(api$replay_ipw_point_if_general, plan = plan2, control = fx2$pc)
    sb <- gn_collect(api$replay_ipw_point_if_view, point_links = fx2$links, control = fx2$pc, callback_schema = TRUE)
    expect_identical(lapply(gb, `[[`, "phi"), lapply(sb, `[[`, "phi"))
    expect_identical(api$stream_ipw_point_covariance_general(plan2, fx2$pc)$covariance,
      api$stream_ipw_point_covariance_view(fx2$links, fx2$pc)$covariance)
  } else {
    expect_s3_class(fx2, "cdid_error")
  }
})

# ---------------------------------------------------------------------------
# 3. Refusals: altered controls/chunks/flags, mutated/reordered/foreign records
#    and links, rebound/unlocked references, callback mutation, native identity
#    unavailable. Every refusal happens before any chunk query or callback.
# ---------------------------------------------------------------------------

test_that("record reuse refuses altered controls, flags, forged or foreign records, rebound or unlocked references and callback mutation", {
  fx <- gn_fit_links("mixed", p = 2L, chunk = 7L)
  e <- gn_ns(gn_general_edit, wrap = c("query_ipw_view_chunk", "gmm_sexp_aliases", "ipw_point_if_general_verify"))
  links <- fx$links; pc <- fx$pc
  plan <- e$ipw_point_if_general_plan(links, pc)
  L <- length(links)
  n_chunks <- as.integer(plan$preflight$n_chunks)
  expect_gt(n_chunks, 1L)
  expect_gte(L, 2L)
  cb_calls <- 0L
  cb <- function(b) cb_calls <<- cb_calls + 1L
  refuse <- function(expr, class) {
    gn_reset(e); cb_calls <<- 0L
    expect_error(expr, class = class)
    expect_identical(gn_n(e, "query_ipw_view_chunk"), 0L)
    expect_identical(cb_calls, 0L)
  }
  # The intact record replays; the native alias check ran at every boundary
  # for every bound object (owners, distinct views, each link's two objects).
  gn_reset(e); cb_calls <- 0L
  info <- e$replay_ipw_point_if_general(plan, pc, cb)
  expect_identical(cb_calls, n_chunks)
  expect_identical(info$verifications, 1L + n_chunks)
  expect_identical(gn_n(e, "ipw_point_if_general_verify"), 1L + n_chunks)
  expect_identical(gn_n(e, "query_ipw_view_chunk"), n_chunks * L)
  n_fv <- plan$general$n_distinct_fitted_nuisances; n_nv <- plan$general$n_distinct_nuisance_views
  expect_gte(gn_n(e, "gmm_sexp_aliases"), (1L + n_chunks) * (3L + 2L * L + n_fv + n_nv))

  # Altered controls (a different chunk step, a different cap): refused.
  pc_chunk <- pc; pc_chunk$max_chunk_rows <- pc$max_chunk_rows + 1L
  refuse(e$replay_ipw_point_if_general(plan, pc_chunk, cb), "cdid_ipw_control_mismatch_error")
  refuse(e$stream_ipw_point_covariance_general(plan, pc_chunk), "cdid_ipw_control_mismatch_error")
  pc_work <- pc; pc_work$max_work_ops <- pc$max_work_ops - 1
  refuse(e$replay_ipw_point_if_general(plan, pc_work, cb), "cdid_ipw_control_mismatch_error")
  # The validated form of the same control is the record's own control: accepted.
  gn_reset(e); cb_calls <- 0L
  e$replay_ipw_point_if_general(plan, api$validate_ipw_point_view_control(pc), cb)
  expect_identical(cb_calls, n_chunks)
  # Flags: the record was preflighted without components; malformed scalars refuse.
  refuse(e$ipw_point_if_general_verify(plan, pc, include_components = TRUE), "cdid_ipw_point_view_input_error")
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, callback_schema = NA), "cdid_ipw_point_view_input_error")
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, use_nuisance_cache = "yes"), "cdid_ipw_point_view_input_error")
  refuse(e$replay_ipw_point_if_general(plan, pc, "not a function"), "cdid_ipw_point_view_input_error")
  refuse(e$stream_ipw_point_covariance_general(plan, pc, trace_env = 1), "cdid_ipw_point_view_input_error")
  # Not a record, a mutated copy, a reordered copy: the seal refuses.
  refuse(e$replay_ipw_point_if_general(links, pc, cb), "cdid_ipw_point_view_general_error")
  p_ids <- plan; p_ids$link_ids <- rev(plan$link_ids)
  refuse(e$replay_ipw_point_if_general(p_ids, pc, cb), "cdid_ipw_point_view_general_error")
  p_rev <- plan; p_rev$point_links <- plan$point_links[rev(seq_len(L))]
  refuse(e$replay_ipw_point_if_general(p_rev, pc, cb), "cdid_ipw_point_view_general_error")
  p_ctrl <- plan; p_ctrl$control_in <- pc_chunk
  refuse(e$replay_ipw_point_if_general(p_ctrl, pc_chunk, cb), "cdid_ipw_point_view_general_error")
  # Forged reseals: a consistently reordered link set whose normalization is
  # not reordered, a reordered normalization whose link IDs are not, a foreign
  # link (same data, separate owners) spliced into the record.
  forge <- function(edit) { x <- plan; x$ref <- NULL; x <- edit(x); api$ipw_owner_seal(x, "point_if_general_plan") }
  f1 <- forge(function(x) { o <- rev(seq_len(L)); x$point_links <- x$point_links[o]; x$link_views <- x$link_views[o]
    x$link_fitted_index <- x$link_fitted_index[o]; x$link_nuisance_index <- x$link_nuisance_index[o]; x })
  refuse(e$replay_ipw_point_if_general(f1, pc, cb), "cdid_ipw_point_view_general_error")
  f2 <- forge(function(x) { o <- rev(seq_len(L)); x$point_links <- x$point_links[o]; x$link_views <- x$link_views[o]
    x$link_fitted_index <- x$link_fitted_index[o]; x$link_nuisance_index <- x$link_nuisance_index[o]
    x$norm$links <- x$norm$links[o]; x$norm$validated <- x$norm$validated[o]; x })
  refuse(e$replay_ipw_point_if_general(f2, pc, cb), "cdid_ipw_point_view_general_error")
  fx_f <- gn_fit_links("mixed", p = 2L, chunk = 7L)
  foreign <- fx_f$links[[1L]]
  expect_false(identical(foreign$frame_ref, links[[1L]]$frame_ref))
  f3 <- forge(function(x) { x$point_links[[1L]] <- foreign; x })
  refuse(e$replay_ipw_point_if_general(f3, pc, cb), "cdid_ipw_point_view_general_error")
  f4 <- forge(function(x) { x$point_links[[1L]] <- foreign; x$norm$links[[1L]] <- foreign; x })
  refuse(e$replay_ipw_point_if_general(f4, pc, cb), "cdid_ipw_point_view_general_error")
  f5 <- forge(function(x) { x$link_views[[1L]] <- get("link_view", envir = foreign$link_view_ref, inherits = FALSE); x })
  refuse(e$replay_ipw_point_if_general(f5, pc, cb), "cdid_ipw_point_view_general_error")
  f6 <- forge(function(x) { x$frame_owner <- fx_f$res$frame_owner; x })
  refuse(e$replay_ipw_point_if_general(f6, pc, cb), "cdid_ipw_point_view_general_error")
  # A record built on the foreign links is valid on its own: the refusals above
  # are about mixing, not about the foreign objects.
  plan_f <- e$ipw_point_if_general_plan(fx_f$links, fx_f$pc)
  gn_reset(e); cb_calls <- 0L
  e$replay_ipw_point_if_general(plan_f, fx_f$pc, cb)
  expect_identical(cb_calls, as.integer(plan_f$preflight$n_chunks))

  # Rebound locked references (owner, fitted nuisance, nuisance, link, point
  # link): a content-modified copy in the live binding refuses; restored, the
  # same record replays again (verification is live, never cached).
  pl2 <- plan$point_links[[2L]]
  lv2 <- plan$link_views[[2L]]
  fv1 <- plan$fitted_views[[1L]]
  nv1 <- plan$nuisance_views[[1L]]
  fo <- plan$frame_owner
  bump <- function(x) { x$estimate <- if (is.null(x$estimate)) 1 else x$estimate + 1; x }
  for (case in list(
    list(env = pl2$ref, nm = "point_link_view", value = bump(pl2)),
    list(env = lv2$ref, nm = "link_view", value = bump(lv2)),
    list(env = fv1$ref, nm = "fitted_nuisance_view", value = bump(fv1)),
    list(env = nv1$ref, nm = "nuisance_view", value = bump(nv1)),
    list(env = fo$ref, nm = "frame_owner", value = bump(fo)),
    list(env = fo$panel_ref, nm = "panel_owner", value = bump(plan$panel_owner)))) {
    restore <- gn_rebind(case$env, case$nm, case$value)
    refuse(e$replay_ipw_point_if_general(plan, pc, cb), "cdid_ipw_point_view_general_error")
    refuse(e$stream_ipw_point_covariance_general(plan, pc), "cdid_ipw_point_view_general_error")
    restore()
    gn_reset(e); cb_calls <- 0L
    e$replay_ipw_point_if_general(plan, pc, cb)
    expect_identical(cb_calls, n_chunks, info = case$nm)
  }
  # An identical object rebound into an unlocked binding still refuses: the
  # reference must be locked.
  unlockBinding("point_link_view", pl2$ref)
  refuse(e$replay_ipw_point_if_general(plan, pc, cb), "cdid_ipw_point_view_general_error")
  lockBinding("point_link_view", pl2$ref)
  gn_reset(e); cb_calls <- 0L
  e$replay_ipw_point_if_general(plan, pc, cb)
  expect_identical(cb_calls, n_chunks)

  # Callback mutation: a consumer that rebinds a point link during the first
  # batch is refused right after that callback, before the next chunk's reads.
  restore <- NULL
  mut_cb <- function(b) {
    cb_calls <<- cb_calls + 1L
    if (is.null(restore)) restore <<- gn_rebind(pl2$ref, "point_link_view", bump(pl2))
  }
  gn_reset(e); cb_calls <- 0L
  expect_error(e$replay_ipw_point_if_general(plan, pc, mut_cb), class = "cdid_ipw_point_view_general_error")
  expect_identical(cb_calls, 1L)
  expect_identical(gn_n(e, "query_ipw_view_chunk"), L)
  restore()
  gn_reset(e); cb_calls <- 0L
  e$stream_ipw_point_covariance_general(plan, pc)
  expect_identical(gn_n(e, "query_ipw_view_chunk"), 2L * n_chunks * L)

  # Native identity check unavailable: no record check can pass and nothing
  # falls back to identity-free trust or to a repeated normalization.
  e_dll <- gn_ns(c(gn_general_edit, list(gmm_sexp_aliases = function(a, b) stop(structure(
    class = c("cdid_gmm_native_error", "cdid_gmm_refusal", "cdid_gmm_error", "error", "condition"),
    list(message = "native identity check unavailable", call = NULL))))),
    wrap = c("query_ipw_view_chunk", "normalize_ipw_point_links_view"))
  cb_calls <- 0L
  expect_error(e_dll$replay_ipw_point_if_general(plan, pc, cb), class = "cdid_gmm_native_error")
  expect_identical(cb_calls, 0L)
  expect_identical(gn_n(e_dll, "query_ipw_view_chunk"), 0L)
  expect_identical(gn_n(e_dll, "normalize_ipw_point_links_view"), 0L)
})

# ---------------------------------------------------------------------------
# 4. Explicit finite caps at the adjacent boundary, and the standalone wrappers'
#    unchanged full validation and refusals.
# ---------------------------------------------------------------------------

test_that("record header and verification charges refuse at the adjacent caps before any chunk; standalone wrappers still normalize and refuse", {
  fx <- gn_fit_links("half", p = 1L, chunk = 7L)
  e <- gn_ns(gn_general_edit, wrap = c("query_ipw_view_chunk", "normalize_ipw_point_links_view", "ipw_point_if_general_verify"))
  links <- fx$links; pc <- fx$pc
  plan <- e$ipw_point_if_general_plan(links, pc)
  pre <- plan$preflight
  g <- plan$general
  L <- length(links)
  # Work: the whole-replay bound plus the four passes' verifications.
  need_work <- pre$total_work_ops_bound + g$additional_work_ops_bound
  expect_identical(g$verifications_planned, 1 + 4 * (as.double(pre$n_chunks) + 1))
  expect_identical(g$construction_work_ops_bound, 16 * as.double(L)^2 + 64 * L + 256)
  expect_identical(g$additional_work_ops_bound, g$construction_work_ops_bound + g$verification_work_ops_bound)
  expect_identical(g$verification_work_ops_bound,
    (1 + 4 * (as.double(pre$n_chunks) + 1)) * (64 * L + 128 * pre$n_unique_fitted_nuisances + 256))
  low <- pc; low$max_work_ops <- need_work - 1
  gn_reset(e)
  expect_error(e$ipw_point_if_general_plan(links, low), class = "cdid_ipw_owner_bound_error")
  expect_identical(gn_n(e, "query_ipw_view_chunk"), 0L)
  exact <- pc; exact$max_work_ops <- need_work
  plan_w <- e$ipw_point_if_general_plan(links, exact)
  expect_identical(plan_w$general$verification_work_ops_bound, g$verification_work_ops_bound)
  # Entries: the producer's concurrent entries plus the record header.
  need_me <- pre$peak_concurrent_entries + g$record_header_entries
  expect_identical(g$additional_live_entries, g$record_header_entries)
  low <- pc; low$max_matrix_entries <- need_me - 1
  gn_reset(e)
  expect_error(e$ipw_point_if_general_plan(links, low), class = "cdid_ipw_owner_bound_error")
  expect_identical(gn_n(e, "query_ipw_view_chunk"), 0L)
  exact <- pc; exact$max_matrix_entries <- need_me
  plan_m <- e$ipw_point_if_general_plan(links, exact)
  expect_identical(plan_m$general$record_header_entries, g$record_header_entries)
  # Covariance on the record: the standalone stream's reservation plus the
  # record header, refused at the adjacent cap before any query.
  cov_entries <- as.double(L) * L
  need_cov <- pre$peak_concurrent_entries + 4 * cov_entries + 8 * L + 6 * pre$b_by_L_entries + g$record_header_entries + 128
  low <- pc; low$max_matrix_entries <- need_cov - 1
  plan_low <- e$ipw_point_if_general_plan(links, low)
  gn_reset(e)
  expect_error(e$stream_ipw_point_covariance_general(plan_low, low), class = "cdid_ipw_owner_bound_error")
  expect_identical(gn_n(e, "query_ipw_view_chunk"), 0L)
  exact <- pc; exact$max_matrix_entries <- need_cov
  plan_ex <- e$ipw_point_if_general_plan(links, exact)
  cg <- e$stream_ipw_point_covariance_general(plan_ex, exact)
  expect_identical(cg$covariance, api$stream_ipw_point_covariance_view(links, exact)$covariance)

  # Standalone wrappers: full per-call validation (one normalization per
  # replay call, three per covariance stream), their original refusals.
  gn_reset(e)
  e$replay_ipw_point_if_view(links, pc, callback = function(b) NULL)
  expect_identical(gn_n(e, "normalize_ipw_point_links_view"), 1L)
  gn_reset(e)
  e$stream_ipw_point_covariance_view(links, pc)
  expect_identical(gn_n(e, "normalize_ipw_point_links_view"), 3L)
  expect_identical(gn_n(e, "ipw_point_if_general_verify"), 0L)
  fx_f <- gn_fit_links("half", p = 1L, chunk = 7L)
  mixed_owners <- unname(c(links[1L], fx_f$links[1L]))
  gn_reset(e)
  expect_error(e$replay_ipw_point_if_view(mixed_owners, pc, callback = function(b) NULL),
    class = "cdid_ipw_owner_mismatch_error")
  expect_identical(gn_n(e, "query_ipw_view_chunk"), 0L)
  low_work <- pc; low_work$max_work_ops <- 5
  gn_reset(e)
  expect_error(e$stream_ipw_point_covariance_view(links, low_work), class = "cdid_ipw_owner_bound_error")
  expect_error(e$ipw_point_if_general_plan(links, low_work), class = "cdid_ipw_owner_bound_error")
  expect_identical(gn_n(e, "query_ipw_view_chunk"), 0L)
  # The record builder is as strict as the standalone preflight on its inputs.
  expect_error(e$ipw_point_if_general_plan(mixed_owners, pc), class = "cdid_ipw_owner_mismatch_error")
  expect_error(e$ipw_point_if_general_plan(links), class = "cdid_ipw_point_view_input_error")
  expect_error(e$ipw_point_if_general_plan(list(), pc), class = "cdid_ipw_point_view_input_error")

  # Source checks: the record consumers never normalize or preflight; the
  # visitors no longer call the set helper.
  src <- function(f) paste(deparse(f), collapse = "\n")
  expect_false(grepl("normalize_ipw_point_links_view", src(api$replay_ipw_point_if_general), fixed = TRUE))
  expect_false(grepl("preflight_ipw_point_if_view", src(api$replay_ipw_point_if_general), fixed = TRUE))
  expect_false(grepl("normalize_ipw_point_links_view", src(api$stream_ipw_point_covariance_general), fixed = TRUE))
  expect_false(grepl("preflight_ipw_point_if_view", src(api$stream_ipw_point_covariance_general), fixed = TRUE))
  expect_true(grepl("ipw_point_if_general_verify", src(api$replay_ipw_point_if_general), fixed = TRUE))
  expect_false(grepl("intersect(", src(api$ipw_point_view_verify_compact_storage), fixed = TRUE))
  expect_false(grepl("intersect(", src(api$ipw_fitted_view_verify_compact_storage), fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 5. Compact-storage visitors: identical first offending field, path, classes
#    and dimension checks against frozen copies of the previous bodies.
# ---------------------------------------------------------------------------

# Frozen copy of the previous `ipw_point_view_verify_compact_storage` body.
gn_old_point_visitor <- function(obj, p_max = NULL, m_max = NULL) {
  forbidden_fields <- c(
    "risk_a_idx", "risk_gi", "risk_w", "risk_dy", "dy_c", "dy_t",
    "Z_use", "D_use", "fit_weights", "eta_use", "odds_use", "prob_use",
    "is_t", "is_c", "a_obs", "b_obs", "w_hess", "score_factor",
    "y_left_obs", "y_right_obs", "Delta_iv_obs", "observed_a_index",
    "fitted.values", "linear.predictors", "residuals", "effects",
    "weights", "y", "qr", "phi", "phi_matrix"
  )
  dim_cap <- if (!is.null(p_max) && !is.null(m_max)) max(as.integer(p_max), as.integer(m_max)) else NULL
  check_node <- function(x, path) {
    if (is.environment(x) || is.null(x)) return(invisible(TRUE))
    if (is.list(x)) {
      nms <- names(x)
      if (!is.null(nms)) {
        bad <- intersect(nms, forbidden_fields)
        if (length(bad) > 0L) {
          api$cdid_ipw_point_view_abort(
            paste0("Forbidden non-compact observer/model field `", bad[1L], "` retained at `", path, "`."),
            c("cdid_ipw_point_view_storage_error", "cdid_ipw_owner_integrity_error")
          )
        }
      }
      for (i in seq_along(x)) {
        nm_i <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
        check_node(x[[i]], paste0(path, "$", nm_i))
      }
      return(invisible(TRUE))
    }
    if (!is.null(dim_cap) && is.atomic(x) && (is.numeric(x) || is.integer(x) || is.logical(x))) {
      d <- dim(x)
      if (!is.null(d)) {
        if (any(d > dim_cap)) {
          api$cdid_ipw_point_view_abort(
            paste0("Non-compact matrix dimension (", paste(d, collapse = "x"), ") retained at `", path, "`."),
            c("cdid_ipw_point_view_storage_error", "cdid_ipw_owner_integrity_error")
          )
        }
      } else if (length(x) > dim_cap) {
        api$cdid_ipw_point_view_abort(
          paste0("Non-compact vector length (", length(x), ") retained at `", path, "`."),
          c("cdid_ipw_point_view_storage_error", "cdid_ipw_owner_integrity_error")
        )
      }
    }
    invisible(TRUE)
  }
  check_node(obj, "point_link_view")
}
# Frozen copy of the previous `ipw_fitted_view_verify_compact_storage` body.
gn_old_fitted_visitor <- function(obj) {
  forbidden_fields <- c(
    "risk_a_idx", "risk_gi", "risk_w", "Z_use", "D_use", "fit_weights",
    "eta_use", "odds_use", "prob_use", "is_t", "is_c", "a_obs", "b_obs",
    "w_hess", "score_factor", "y_left_obs", "y_right_obs", "Delta_iv_obs",
    "observed_a_index", "fitted.values", "linear.predictors", "residuals",
    "effects", "weights", "y", "qr"
  )
  check_node <- function(x, path) {
    if (is.environment(x)) return(invisible(TRUE))
    if (is.list(x)) {
      nms <- names(x)
      if (!is.null(nms)) {
        bad <- intersect(nms, forbidden_fields)
        if (length(bad) > 0L) {
          api$cdid_ipw_fitted_view_abort(
            paste0("Forbidden non-compact observer/model field `", bad[1L], "` retained at `", path, "`."),
            c("cdid_ipw_fitted_view_storage_error", "cdid_ipw_owner_integrity_error")
          )
        }
      }
      for (i in seq_along(x)) {
        nm_i <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
        check_node(x[[i]], paste0(path, "$", nm_i))
      }
    }
    invisible(TRUE)
  }
  check_node(obj, "fitted_nuisance_view")
}

test_that("compact-storage visitors report the identical first offending field, path, classes and dimension refusals as their frozen previous bodies", {
  na_named <- list(1, 2); names(na_named) <- c(NA_character_, "effects")
  empty_named <- list(1, 2); names(empty_named) <- c("", "weights")
  dup_names <- list(1, 2, 3); names(dup_names) <- c("zz", "y", "y")
  two_bad <- list(ok = 1, qr = 2, y = 3)
  two_bad_rev <- list(ok = 1, y = 3, qr = 2)
  nested_unnamed <- list(list(a = 1, list(fitted.values = 1:2)))
  deep <- list(ok = list(deeper = list(phi = 1), fine = 2), tail = list(phi_matrix = matrix(0, 2, 2)))
  clean <- list(a = 1, b = list(c = "y", d = c(qr = 1)), e = new.env(), f = NULL, g = c(y = 2, phi = 3))
  mats <- list(m3 = matrix(0, 3, 3), v2 = numeric(2))
  vec5 <- list(a = list(v = numeric(5)))
  lgl <- list(a = matrix(TRUE, 1, 4))
  int_ok <- list(a = 1:3)
  chr_long <- list(a = rep("x", 50))
  attr_list <- structure(list(a = 1, risk_w = 2), class = "foo", extra = "bar")
  cases <- list(na_named, empty_named, dup_names, two_bad, two_bad_rev, nested_unnamed, deep, clean, mats, vec5, lgl,
    int_ok, chr_long, attr_list, list(y = NULL), list(weights = new.env()), list(), list(1, 2, 3), NULL)
  caps <- list(list(NULL, NULL), list(2L, 2L), list(1L, 3L), list(4, 1))
  for (i in seq_along(cases)) for (cp in caps) {
    x <- cases[[i]]
    new <- gn_outcome(api$ipw_point_view_verify_compact_storage(x, p_max = cp[[1L]], m_max = cp[[2L]]))
    old <- gn_outcome(gn_old_point_visitor(x, p_max = cp[[1L]], m_max = cp[[2L]]))
    expect_identical(new, old, info = paste("point visitor case", i, "caps", paste(unlist(cp), collapse = ",")))
    new_f <- gn_outcome(api$ipw_fitted_view_verify_compact_storage(x))
    old_f <- gn_outcome(gn_old_fitted_visitor(x))
    expect_identical(new_f, old_f, info = paste("fitted visitor case", i))
  }
  # The first offending field is the first in `names(x)` order, never the
  # first in the forbidden list.
  expect_match(gn_outcome(api$ipw_point_view_verify_compact_storage(two_bad))$message, "`qr` retained at `point_link_view`", fixed = TRUE)
  expect_match(gn_outcome(api$ipw_point_view_verify_compact_storage(two_bad_rev))$message, "`y` retained at `point_link_view`", fixed = TRUE)
  expect_match(gn_outcome(api$ipw_point_view_verify_compact_storage(nested_unnamed))$message,
    "`fitted.values` retained at `point_link_view$1$2`", fixed = TRUE)
  # `phi` is forbidden in point links only; the fitted-view list is unchanged.
  expect_match(gn_outcome(api$ipw_point_view_verify_compact_storage(deep))$message,
    "`phi` retained at `point_link_view$ok$deeper`", fixed = TRUE)
  expect_identical(gn_outcome(api$ipw_fitted_view_verify_compact_storage(deep)), "ok")
  expect_match(gn_outcome(api$ipw_point_view_verify_compact_storage(mats, 2L, 2L))$message, "(3x3) retained at `point_link_view$m3`", fixed = TRUE)
  expect_match(gn_outcome(api$ipw_point_view_verify_compact_storage(vec5, 2L, 2L))$message, "(5) retained at `point_link_view$a$v`", fixed = TRUE)
  expect_identical(gn_outcome(api$ipw_point_view_verify_compact_storage(vec5)), "ok")
  # Real sealed objects pass both; a retained observer-length field refuses both identically.
  fx <- gn_fit_links("mixed", p = 2L, chunk = 7L)
  pl <- fx$links[[1L]]
  fv <- fx$res$nuisance_registry$fitted_views[[1L]]
  p <- as.integer(fx$res$frame_owner$n_covariates)
  m <- as.integer(get("nuisance_view", envir = fv$nuisance_view_ref, inherits = FALSE)$cohorts$m_cohorts)
  expect_identical(gn_outcome(api$ipw_point_view_verify_compact_storage(pl, p_max = p, m_max = m)), "ok")
  expect_identical(gn_outcome(gn_old_point_visitor(pl, p_max = p, m_max = m)), "ok")
  expect_identical(gn_outcome(api$ipw_fitted_view_verify_compact_storage(fv)), "ok")
  expect_identical(gn_outcome(gn_old_fitted_visitor(fv)), "ok")
  pl_bad <- pl; pl_bad$phi <- numeric(3)
  expect_identical(gn_outcome(api$ipw_point_view_verify_compact_storage(pl_bad, p, m)), gn_outcome(gn_old_point_visitor(pl_bad, p, m)))
  expect_s3_class(tryCatch(api$ipw_point_view_verify_compact_storage(pl_bad, p, m), error = identity), "cdid_ipw_point_view_storage_error")
  pl_wide <- pl; pl_wide$means$extra <- matrix(0, p + m + 1L, 1L)
  expect_identical(gn_outcome(api$ipw_point_view_verify_compact_storage(pl_wide, p, m)), gn_outcome(gn_old_point_visitor(pl_wide, p, m)))
  fv_bad <- fv; fv_bad$factors$weights <- 1
  expect_identical(gn_outcome(api$ipw_fitted_view_verify_compact_storage(fv_bad)), gn_outcome(gn_old_fitted_visitor(fv_bad)))
  expect_s3_class(tryCatch(api$ipw_fitted_view_verify_compact_storage(fv_bad), error = identity), "cdid_ipw_fitted_view_storage_error")
})
