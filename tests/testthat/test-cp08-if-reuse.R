# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 retained raw influence-function reuse through its consumers:
#   (1) IPW batch route: the first covariance replay's raw N_A x L phi is kept
#       in an invocation-local carrier and sliced, in the plan's own partition,
#       for the covariance Gram pass and both GMM projection passes (one
#       producer pass, three slice passes).
#   (2) OR/DR complete-support route: the admitted CP04 raw store stays alive
#       through auxiliary reporting and serves the exact same-recipe auxiliary
#       columns (frame scale on the F = A, lambda = 1 frame, existing sign,
#       original output slot); unmatched columns replay as before.
# Every expected value comes from the unchanged forced-stream route of the same
# installed build (a re-bound namespace copy whose reservation declines), from
# the ordinary plan/producer replays of the same fitted objects, from the four
# reported IF components, or from the existing accounting outputs; the retained
# routes never supply their own truth and no existing gate or tolerance is
# weakened (floating parity at 1e-12, identities/orders/shapes/counts exact).
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
# Shared test-local helpers (clones of the accepted CP08 suites' fixtures)
# ---------------------------------------------------------------------------

# Re-bound copy of the namespace (as xb_ns / cc_counted): native resolution
# stays genuine; `edits` replace internal bindings; `wrap` names are counted.
ru_ns <- function(edits = list(), wrap = character()) {
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
ru_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0L)
ru_strip_env <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); y <- lapply(x, ru_strip_env); attributes(y) <- a; return(y) }
  x
}
# Strict schema/leaf comparison, doubles at `tol`, everything else exactly.
ru_expect_match <- function(a, b, tol = 1e-12, info = NULL) {
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
ru_reseal <- function(x, binding, edit) { x$ref <- NULL; x <- edit(x); api$ipw_owner_seal(x, binding) }
ru_rebind <- function(env, nm, value) {
  old <- get(nm, envir = env)
  unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
  function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
}
ru_collect <- function(fn, ...) {
  out <- list()
  fn(..., callback = function(b) { force(b); out[[length(out) + 1L]] <<- b })
  out
}
# No retained carrier/store survives in a result (environments by kind only).
ru_find_retained <- function(x, depth = 0L) {
  if (depth > 12L) return(FALSE)
  if (is.environment(x)) return(exists("kind", envir = x, inherits = FALSE) &&
    get("kind", envir = x, inherits = FALSE) %in% c("cp04_retained_if", "ipw_point_if_retained"))
  if (inherits(x, "cdid_cp04_retained_if_bound")) return(TRUE)
  if (is.list(x)) for (v in x) if (ru_find_retained(v, depth + 1L)) return(TRUE)
  FALSE
}

# ----- IPW fixtures (clones of test-cp08-ipw-x-batch.R) --------------------
ru_ipw_panel <- function(N = 41L, T = 6L, extra = NULL, drop = NULL) {
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
  d$w_two <- 2
  d$w_vary <- 1 + (d$id %% 4L) / 8
  if (!is.null(drop)) d <- d[!(paste(d$id, d$t) %in% drop), ]
  d
}
ru_ipw_formula <- function(p) switch(as.character(p), "1" = ~1, "2" = ~x1, "3" = ~x1 + x2, stop("p must be 1, 2 or 3"))
ru_ipw_ctrl <- function(d, p, chunk = 7L, weighting = "optimal") {
  gs <- unique(d$g)
  n_frame <- length(unique(d$id))
  dims <- api$chained_public_dimensions_from_counts(n_frame = n_frame, n_panel_units = n_frame,
    n_rows = nrow(d), n_dates = length(unique(d$t)), n_cohorts = length(gs), n_targets = sum(is.finite(gs)),
    n_covariates = p, biters = 99, n_clusters = n_frame)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(list(chunk_rows = chunk),
    gmm_weighting = weighting), gmm_weighting = weighting)$control
}
ru_ipw_run <- function(ns, d, ctrl, p, ...) ns$run_chained_ipw_pipeline(
  data = d, idname = "id", tname = "t", yname = "y", gname = "g", control = ctrl, xformla = ru_ipw_formula(p),
  control_group = "nevertreated", base_period = "varying", weighting = "optimal", run_numerical_transfer = FALSE,
  root_witness = NULL, structural_reduction_stage = "outcome_free", ...)
# Forced plan replay of the same build: the carrier reservation declines.
ru_ipw_stream_edit <- list(ipw_point_if_retained_reserve = function(...) NULL)
ru_ipw_counts <- c("ipw_point_if_plan_chunk", "replay_ipw_point_if_plan", "replay_ipw_point_if_retained",
  "ipw_point_if_retained_fill")
ru_ipw_fields <- c("support", "basis", "aggregations", "gmm_fit", "structural_certificate", "design_all", "design_used",
  "covariance_stream", "covariance_handling", "full_frame_transport", "auxiliary_reporting")
ru_ipw_record <- function(res) {
  if (!is.null(res$resource_accounting$point_batch$retained_if)) res$resource_accounting$point_batch$retained_if else
    res$resource_accounting$point_general$retained_if
}
ru_ipw_expect_same <- function(a, b, info = NULL) {
  for (f in ru_ipw_fields) ru_expect_match(ru_strip_env(a[[f]]), ru_strip_env(b[[f]]), info = paste(info, f))
  for (nm in c("estimates", "se", "analytical_se", "if_full", "vcov", "available", "ci_lower", "ci_upper", "p_value")) {
    ru_expect_match(a$reporting[[nm]], b$reporting[[nm]], info = paste(info, "reporting", nm))
  }
  expect_identical(a$reporting$target_meta, b$reporting$target_meta, info = info)
  rec_a <- ru_ipw_record(a); rec_b <- ru_ipw_record(b)
  add_a <- if (!is.null(rec_a) && identical(rec_a$status, "retained")) rec_a$reserved_entries else 0
  add_b <- if (!is.null(rec_b) && identical(rec_b$status, "retained")) rec_b$reserved_entries else 0
  # Declining the carrier can admit D90's sparse general replay state instead.
  # Check each complete composition before comparing the common base charge.
  base_entries <- function(x, carrier_entries) {
    ra <- x$resource_accounting
    state <- ra$point_general$replay_state
    state_entries <- if (isTRUE(state$admitted)) state$reserved_entries else 0
    expected <- sum(c(ra$initial_envelope_entries,
      ra$point_producer$producer_workspace_entries, ra$gmm_consumer$total_payload_entries,
      ra$d5_arrays$ledger_if_entries, ra$d5_arrays$basis_if_entries,
      ra$d5_arrays$ledger_enclosure_entries, ra$d5_arrays$basis_enclosure_entries,
      ra$d5_arrays$retained_diagnostic_entries, ra$point_batch$if_plan$additional_live_entries,
      ra$point_general$additional_live_entries, ra$point_custom$additional_live_entries))
    expect_identical(ra$simultaneous_entries, expected + carrier_entries + state_entries, info = info)
    ra$simultaneous_entries - carrier_entries - state_entries
  }
  expect_identical(base_entries(a, add_a), base_entries(b, add_b), info = info)
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches, info = info)
  expect_identical(a$gmm_fit$if_projection$max_batch_rows_seen, b$gmm_fit$if_projection$max_batch_rows_seen, info = info)
}

# ----- OR/DR fixtures (clones of test-cp08-retained-moment-if.R) -----------
rm_ctrl <- function(chunk = 7L, work = 1e10, matrix = 2e7) api$ipw_owner_control(max_chunk_rows = chunk,
  max_retained_entries = 1e7, max_matrix_entries = matrix, max_work_ops = work, max_output_entries = 1e7)
rm_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
rm_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))
rm_panel <- function(N = 41L, T = 5L, seed = 11L) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:2] <- c(3, Inf)
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$y <- 0.2 * d$id + 0.4 * d$t + ifelse(d$t >= d$g, 1.25, 0) + stats::rnorm(nrow(d))
  d
}
rm_pairs <- function(T = 5L) {
  pr <- list()
  for (s in seq_len(T - 1L)) for (t in seq.int(s + 1L, T)) pr[[length(pr) + 1L]] <- c(s, t)
  pr
}
# Real sealed owners, reader and link views. `aux_pairs` appends views of the
# same target/pool on repeated pairs under their own labels, as the pipeline's
# auxiliary descriptors are appended at M + j on the same sealed backend.
rm_fixture <- function(d = rm_panel(), chunk = 7L, pairs = rm_pairs(), aux_pairs = list(), K = Inf, target = 3,
                       ctrl = rm_ctrl(chunk)) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  u <- d[!duplicated(d$id), c("id", "g")]
  u <- u[order(u$id), ]
  N <- nrow(u)
  Z <- matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)"))
  V_mat <- matrix(1, N, 1L, dimnames = list(NULL, "reader_intercept"))
  wt <- rep(1, N)
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = wt, control = ctrl, frame_g = u$g,
    allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, u$id, V_mat, u$g, wt, ctrl, allow_unknown_off_A = TRUE)
  views <- lapply(seq_along(pairs), function(k) {
    p <- pairs[[k]]
    api$construct_ipw_link_view(fo, p[1], p[2], g = target, K = K, control = ctrl,
      link_id = paste0("L", p[1], "_", p[2]), functional_key = paste0("cmpcat:rm:", p[1], "-", p[2]))
  })
  for (k in seq_along(aux_pairs)) {
    p <- aux_pairs[[k]]
    views[[length(views) + 1L]] <- api$construct_ipw_link_view(fo, p[1], p[2], g = target, K = K, control = ctrl,
      link_id = paste0("A", k, "_", p[1], "_", p[2]), functional_key = paste0("cmpcat:rm-aux:", k))
  }
  list(d = d, po = po, fo = fo, reader = reader, views = views, ctrl = ctrl, chunk = chunk, u = u, pairs = pairs,
    M = length(pairs), aux_index = length(pairs) + seq_along(aux_pairs))
}
rm_or <- function(fx, ns = api, ctrl = fx$ctrl) ns$fit_or_baseline(fx$reader, fx$views, ctrl, 1e-10, outcome_recipe = rm_icpt)
rm_dr <- function(fx, ns = api, ctrl = fx$ctrl) ns$fit_dr_baseline(fx$reader, fx$views, ctrl, glm_control = rm_glm,
  rank_tol = 1e-10, outcome_recipe = rm_icpt, treatment_columns = NULL)
rm_nonadd_panel <- function(N = 41L, T = 5L, n_g = 18L, seed = 31L, id_perm = FALSE, rows_perm = FALSE, g_date = 3) {
  set.seed(seed)
  g <- rep(0, N)
  g[sample.int(N, n_g)] <- g_date
  ids <- if (id_perm) sample(1000L + seq_len(5L * N), N) else seq_len(N)
  xi <- stats::rbinom(N, 1, 0.5)
  ue <- stats::rnorm(N)
  d <- data.frame(id = rep(ids, each = T), period = rep(seq_len(T), N))
  u <- match(d$id, ids)
  d$G <- g[u]
  d$x <- xi[u]
  eff <- ifelse(d$G > 0 & d$period >= d$G, 1 + 0.25 * (d$period - d$G), 0)
  d$y <- ue[u] + 0.1 * d$period + 0.35 * sin(0.7 * u * d$period) + 0.2 * d$x * (d$period^2) + eff + stats::rnorm(nrow(d))
  d$cl <- (u - 1L) %/% 4L + 1L
  d$w2 <- 2
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}
rm_run_pipeline <- function(d, est = "reg", chunk = 7L, reduce = TRUE, weighting = "optimal",
                            xformla = NULL, weightsname = NULL, control_group = "nevertreated",
                            frame = NULL, aggregations = c("simple", "group", "dynamic", "calendar"),
                            condition_trigger = NULL, ridge_relative = NULL, owner_caps = NULL, clustered = FALSE,
                            base_period = "varying", compute_inffunc = TRUE, ns = api) {
  ctrl_list <- if (is.null(chunk)) list() else list(chunk_rows = chunk)
  if (!is.null(condition_trigger)) ctrl_list$condition_trigger <- condition_trigger
  if (!is.null(ridge_relative)) ctrl_list$ridge_relative <- ridge_relative
  if (!length(ctrl_list)) ctrl_list <- NULL
  resolved <- api$chained_public_resolve_control(ctrl_list, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  prep <- api$chained_request_prepare_inputs(
    data = d, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xformla, weightsname = weightsname, clustervars = if (clustered) "cl" else NULL, anticipation = 0,
    frame = frame, outcome_regression = identical(est, "reg"),
    doubly_robust = identical(est, "dr"), or_resource = resolved$profile
  )
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = 99L, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  if (!is.null(owner_caps)) for (nm in names(owner_caps)) built$control$owner_control[[nm]] <- owner_caps[[nm]]
  tc <- prep$target_cohorts
  ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = built$control, xformla = prep$xformla, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = tc, anticipation = 0L, control_group = control_group, base_period = base_period,
    weighting = weighting, reduce_structural_redundancy = reduce, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = aggregations, cband = FALSE, alp = 0.05, biters = 99L, seed = 101L,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = FALSE, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = "post", pt_start = NULL, cp04 = prep$cp04,
    compute_inffunc = compute_inffunc
  )
}
# Forced-stream reference of the same build: the store reservation declines.
rm_stream_ns <- function() {
  e <- new.env(parent = api)
  f <- get("run_chained_ipw_pipeline", api); environment(f) <- e; e$run_chained_ipw_pipeline <- f
  e$cp04_retained_if_reserve <- function(...) NULL
  e
}
# The store is admitted but its auxiliary extension declines (original lifetime).
rm_no_aux_ns <- function() {
  e <- new.env(parent = api)
  f <- get("run_chained_ipw_pipeline", api); environment(f) <- e; e$run_chained_ipw_pipeline <- f
  e$cp04_retained_if_auxiliary_admit <- function(...) NULL
  e
}
# The auxiliary extension's live envelope is replaced by cap + delta.
rm_aux_live_ns <- function(delta) {
  e <- new.env(parent = api)
  f <- get("run_chained_ipw_pipeline", api); environment(f) <- e; e$run_chained_ipw_pipeline <- f
  e$cp04_retained_if_auxiliary_admit <- function(store, backend, backend_ids, retained, live, work, control)
    api$cp04_retained_if_auxiliary_admit(store, backend, backend_ids, retained, live = store$cap + delta, work, control)
  e
}
# Counting namespace: producer chunk evaluations and replay entries, split
# between the pipeline's main passes and the auxiliary reporting consumer.
rm_counting_ns <- function(force_stream = FALSE) {
  e <- new.env(parent = api)
  for (nm in c("run_chained_ipw_pipeline", "cp04_stream_covariance", "cp04_replay_links", "cp04_retained_if_replay",
               "build_ipw_auxiliary_reporting", "cp04_retained_if_auxiliary_slices")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  n <- new.env(parent = emptyenv())
  n$chunks <- 0L; n$replays <- 0L; n$aux_chunks <- 0L; n$aux_replays <- 0L; n$aux_slices <- 0L; n$in_aux <- FALSE
  n$rows <- list(); n$aux_ids <- list()
  e$chained_full_if_chunk <- function(backend, rows, links, scale, plan = NULL) {
    if (n$in_aux) n$aux_chunks <- n$aux_chunks + 1L else { n$chunks <- n$chunks + 1L; n$rows[[length(n$rows) + 1L]] <- rows }
    api$chained_full_if_chunk(backend, rows, links, scale, plan)
  }
  replay_fn <- get("cp04_replay_links", e)
  e$cp04_replay_links <- function(backend, links, callback, scale = "core", chunk_rows = NULL) {
    if (n$in_aux) { n$aux_replays <- n$aux_replays + 1L; n$aux_ids[[length(n$aux_ids) + 1L]] <- as.integer(links) } else
      n$replays <- n$replays + 1L
    replay_fn(backend, links, callback, scale = scale, chunk_rows = chunk_rows)
  }
  aux_fn <- get("build_ipw_auxiliary_reporting", e)
  e$build_ipw_auxiliary_reporting <- function(...) { n$in_aux <- TRUE; on.exit(n$in_aux <- FALSE); aux_fn(...) }
  slices_fn <- get("cp04_retained_if_auxiliary_slices", e)
  e$cp04_retained_if_auxiliary_slices <- function(...) { n$aux_slices <- n$aux_slices + 1L; slices_fn(...) }
  if (force_stream) e$cp04_retained_if_reserve <- function(...) NULL
  e$n <- n
  e
}
rm_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, rm_strip); attributes(x) <- a }
  x
}
rm_record <- function(res) res$resource_accounting$point_producer$retained_moment_if
rm_collect <- function(producer, links, B, scale = "core") {
  out <- list()
  api$cp04_replay_links(producer, links, function(b) { force(b); out[[length(out) + 1L]] <<- b }, scale = scale, chunk_rows = B)
  out
}
rm_expect_parity <- function(a, b, info = NULL) {
  expect_equal(a$covariance_stream$mhat, b$covariance_stream$mhat, tolerance = 1e-12, info = info)
  expect_equal(a$covariance_stream$V, b$covariance_stream$V, tolerance = 1e-12, info = info)
  expect_equal(a$covariance_stream$mean_phi, b$covariance_stream$mean_phi, tolerance = 1e-12, info = info)
  expect_equal(a$gmm_fit$targets, b$gmm_fit$targets, tolerance = 1e-12, info = info)
  expect_equal(a$gmm_fit$if_projection$sandwich_cov, b$gmm_fit$if_projection$sandwich_cov, tolerance = 1e-12, info = info)
  expect_equal(a$gmm_fit$if_projection$materialized_target_if, b$gmm_fit$if_projection$materialized_target_if,
    tolerance = 1e-12, info = info)
  expect_equal(a$covariance_handling, b$covariance_handling, tolerance = 1e-12, info = info)
  expect_equal(a$full_frame_transport$Phi_full, b$full_frame_transport$Phi_full, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$estimates, b$reporting$estimates, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$se, b$reporting$se, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$analytical_se, b$reporting$analytical_se, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$if_full, b$reporting$if_full, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$vcov, b$reporting$vcov, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$Omega_a, b$reporting$Omega_a, tolerance = 1e-12, info = info)
  expect_identical(a$reporting$target_meta, b$reporting$target_meta, info = info)
  expect_identical(a$reporting$available, b$reporting$available, info = info)
  expect_equal(rm_strip(a$auxiliary_reporting), rm_strip(b$auxiliary_reporting), tolerance = 1e-12, info = info)
  for (ag in names(a$aggregations)) {
    expect_equal(a$aggregations[[ag]]$overall$estimate, b$aggregations[[ag]]$overall$estimate, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$overall$se, b$aggregations[[ag]]$overall$se, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$estimates, b$aggregations[[ag]]$estimates, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$se, b$aggregations[[ag]]$se, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$if_full, b$aggregations[[ag]]$if_full, tolerance = 1e-12, info = info)
  }
  expect_identical(a$nuisance_registry$nuisance_keys, b$nuisance_registry$nuisance_keys, info = info)
  expect_identical(a$nuisance_registry$cp04_aux_index, b$nuisance_registry$cp04_aux_index, info = info)
  expect_identical(a$reduction_summary$kept_indices, b$reduction_summary$kept_indices, info = info)
  expect_identical(a$resource_accounting$simultaneous_entries, b$resource_accounting$simultaneous_entries, info = info)
  expect_identical(a$resource_accounting$point_producer$chunk_rows, b$resource_accounting$point_producer$chunk_rows, info = info)
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches, info = info)
  expect_identical(a$gmm_fit$if_projection$max_batch_rows_seen, b$gmm_fit$if_projection$max_batch_rows_seen, info = info)
}
# Independent structural expectation of the auxiliary map: a requested
# supported pretrend column is servable iff its descriptor's oriented pair of
# calendar dates is the pair of a retained catalog link (same target cohort,
# never-treated pool and complete dates on this route). Returns, in reporting
# order, the requested backend indices and TRUE (servable) / FALSE (must stream).
rm_expected_served <- function(res, base_period) {
  kept <- vapply(res$point_links_used, function(pl) paste(as.double(pl$pair$s), as.double(pl$pair$t)), character(1))
  aux <- res$auxiliary_reporting
  sel <- aux$selections[[base_period]]
  idx <- vapply(sel, function(d) if (!is.na(d$index) && d$index > 0L &&
    identical(aux$descriptor_status$status[d$index], "supported")) d$index else NA_integer_, integer(1))
  idx <- idx[!is.na(idx)]
  list(backend_ids = as.integer(res$nuisance_registry$cp04_aux_index[idx]),
    servable = vapply(idx, function(q) paste(as.double(aux$descriptors[[q]]$s), as.double(aux$descriptors[[q]]$t)) %in% kept,
      logical(1)))
}

rm_d <- rm_nonadd_panel(41L, 5L, 18L, seed = 31L, id_perm = TRUE, rows_perm = TRUE)
# A later cohort gives several pretrend comparisons (varying: 1-2, 2-3, 3-4;
# universal with reference 4: 1-4, 2-4, 3-4), so matched and unmatched
# auxiliary columns can coexist on one retained forest.
rm_d5 <- rm_nonadd_panel(41L, 6L, 18L, seed = 47L, id_perm = TRUE, rows_perm = TRUE, g_date = 5)

# ===========================================================================
# PART 1. IPW batch route: retained point IF carrier
# ===========================================================================

test_that("IPW carrier: pipeline outputs equal the forced plan replay at 1e-12 with and without X; one producer pass, three slices", {
  d <- ru_ipw_panel()
  for (p in c(1L, 2L, 3L)) for (chunk in c(7L, 16L)) {
    info <- paste("p =", p, "chunk =", chunk)
    ctrl <- ru_ipw_ctrl(d, p, chunk = chunk)
    e_r <- ru_ns(wrap = ru_ipw_counts)
    e_s <- ru_ns(ru_ipw_stream_edit, wrap = ru_ipw_counts)
    res_r <- ru_ipw_run(e_r, d, ctrl, p)
    res_s <- ru_ipw_run(e_s, d, ctrl, p)
    n_chunks <- as.integer(ceiling(41 / chunk))
    rec <- ru_ipw_record(res_r)
    expect_false(is.null(res_r$resource_accounting$point_batch), info = info)
    expect_identical(rec$status, "retained", info = info)
    expect_identical(rec$producer_passes, 1L, info = info)
    expect_identical(rec$slice_passes, 3L, info = info)
    expect_identical(rec$n_batches, n_chunks, info = info)
    expect_identical(rec$chunk_rows, min(41L, chunk), info = info)
    expect_identical(rec$n_a_units, 41L, info = info)
    expect_identical(rec$n_links, res_r$reduction_summary$n_links_used, info = info)
    plan <- api$ipw_point_if_batch_plan(res_r$point_links_used, ctrl$point_control, res_r$nuisance_registry$fitted_views[[1L]])
    inv <- api$ipw_point_if_retained_inventory(plan)
    expect_identical(inv$entries, api$ipw_point_if_retained_entries(41, rec$n_links, min(41, chunk)) +
      inv$metadata_entries + inv$covariance_entries, info = info)
    expect_identical(rec$reserved_entries, inv$entries, info = info)
    expect_identical(rec$work_ops_bound, inv$work, info = info)
    expect_identical(ru_ipw_record(res_s)$status, "streamed", info = info)
    # One kernel pass instead of four; the slice entry serves the other three.
    expect_identical(ru_n(e_r, "ipw_point_if_plan_chunk"), n_chunks, info = info)
    expect_identical(ru_n(e_s, "ipw_point_if_plan_chunk"), 4L * n_chunks, info = info)
    expect_identical(ru_n(e_r, "ipw_point_if_retained_fill"), 1L, info = info)
    expect_identical(ru_n(e_r, "replay_ipw_point_if_plan"), 1L, info = info)
    expect_identical(ru_n(e_r, "replay_ipw_point_if_retained"), 3L, info = info)
    expect_identical(ru_n(e_s, "replay_ipw_point_if_plan"), 4L, info = info)
    expect_identical(ru_n(e_s, "replay_ipw_point_if_retained"), 0L, info = info)
    # Original GMM batch boundaries, full-F reporting, covariance, SE, IF.
    expect_identical(res_r$gmm_fit$if_projection$n_batches, n_chunks, info = info)
    ru_ipw_expect_same(res_r, res_s, info = info)
    expect_gt(sum(res_r$reporting$available), 0L)
    expect_equal(sqrt(diag(res_r$covariance_stream$V)), sqrt(diag(res_s$covariance_stream$V)), tolerance = 1e-12)
    # The CP04 record stays absent on the IPW route; nothing retained survives.
    expect_null(res_r$resource_accounting$point_producer$retained_moment_if, info = info)
    expect_false(ru_find_retained(res_r), info = info)
    expect_false(ru_find_retained(res_r$resource_accounting$point_batch), info = info)
  }
})

test_that("IPW carrier: the active original-M Tikhonov ridge and the clustered/non-unit-weight declines keep the forced replay's values", {
  d <- ru_ipw_panel()
  p <- 2L
  gs <- unique(d$g)
  dims <- api$chained_public_dimensions_from_counts(n_frame = 41L, n_panel_units = 41L, n_rows = nrow(d),
    n_dates = 6L, n_cohorts = length(gs), n_targets = sum(is.finite(gs)), n_covariates = p, biters = 99, n_clusters = 41L)
  ridge <- api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(
    list(chunk_rows = 7L, condition_trigger = 1 - 1e-12, ridge_relative = 0.01), gmm_weighting = "optimal"),
    gmm_weighting = "optimal")$control
  res_r <- ru_ipw_run(api, d, ridge, p)
  res_s <- ru_ipw_run(ru_ns(ru_ipw_stream_edit), d, ridge, p)
  expect_identical(ru_ipw_record(res_r)$status, "retained")
  expect_identical(res_r$covariance_handling$handler_used, "tikhonov")
  expect_identical(res_r$covariance_handling$gram_provenance, res_s$covariance_handling$gram_provenance)
  expect_equal(res_r$covariance_handling$lambda, res_s$covariance_handling$lambda, tolerance = 1e-12)
  expect_equal(res_r$covariance_handling$trace_B, res_s$covariance_handling$trace_B, tolerance = 1e-12)
  ru_ipw_expect_same(res_r, res_s, info = "ridge")
  # Clustered frame: the batch plan and retained carrier serve.
  ctrl <- ru_ipw_ctrl(d, p, chunk = 7L)
  cl <- rep(seq_len(11L), length.out = 41L)
  res_c <- ru_ipw_run(api, d, ctrl, p, cluster_vec = cl)
  res_cs <- ru_ipw_run(ru_ns(ru_ipw_stream_edit), d, ctrl, p, cluster_vec = cl)
  expect_false(is.null(res_c$resource_accounting$point_batch))
  expect_identical(ru_ipw_record(res_c)$status, "retained")
  ru_ipw_expect_same(res_c, res_cs, info = "clustered")
  # Constant non-unit weights: eligible plan and retained carrier, same values.
  res_w <- ru_ipw_run(api, d, ctrl, p, weightsname = "w_two")
  res_ws <- ru_ipw_run(ru_ns(ru_ipw_stream_edit), d, ctrl, p, weightsname = "w_two")
  expect_false(is.null(res_w$resource_accounting$point_batch))
  expect_identical(ru_ipw_record(res_w)$status, "retained")
  ru_ipw_expect_same(res_w, res_ws, info = "w_two")
  # Unit weights supplied by name: eligible and retained.
  res_1 <- ru_ipw_run(api, d, ctrl, p, weightsname = "w_const")
  expect_identical(ru_ipw_record(res_1)$status, "retained")
  # Unit-varying weights: the general route, no plan, no batch carrier record.
  res_v <- ru_ipw_run(api, d, ctrl, p, weightsname = "w_vary")
  expect_null(res_v$resource_accounting$point_batch)
  # An incomplete date behaves identically on legacy clones and agrees with production raw-IPW at 1e-12.
  d3 <- ru_ipw_panel(drop = "3 2")
  ctrl3 <- ru_ipw_ctrl(d3, p, chunk = 7L)
  legacy_ns <- ru_ns()
  stream_ns <- ru_ns(ru_ipw_stream_edit)
  r3 <- tryCatch(ru_ipw_run(api, d3, ctrl3, p), error = identity)
  r3_leg <- tryCatch(ru_ipw_run(legacy_ns, d3, ctrl3, p), error = identity)
  s3 <- tryCatch(ru_ipw_run(stream_ns, d3, ctrl3, p), error = identity)
  if (inherits(s3, "error")) {
    expect_identical(class(r3), class(s3))
    expect_identical(conditionMessage(r3), conditionMessage(s3))
    expect_identical(class(r3_leg), class(s3))
    expect_identical(conditionMessage(r3_leg), conditionMessage(s3))
  } else {
    expect_null(r3_leg$resource_accounting$point_batch)
    expect_identical(ru_ipw_record(r3_leg)$status, "retained")
    expect_identical(ru_ipw_record(s3)$status, "streamed")
    expect_null(r3_leg$nuisance_registry$raw_ipw)
    expect_null(s3$nuisance_registry$raw_ipw)
    expect_identical(r3_leg$resource_accounting$point_general$replay_state$status, "superseded_by_retained_if")
    expect_false(r3_leg$resource_accounting$point_general$replay_state$admitted)
    expect_identical(s3$resource_accounting$point_general$replay_state$status, "retained_and_released")
    expect_true(s3$resource_accounting$point_general$replay_state$admitted)
    expect_true(s3$resource_accounting$point_general$replay_state$released)
    ru_ipw_expect_same(r3_leg, s3, info = "incomplete date legacy")
    for (nm in c("mhat", "V", "mean_phi"))
      expect_identical(r3_leg$covariance_stream[[nm]], s3$covariance_stream[[nm]], info = paste("incomplete legacy", nm))
    expect_identical(r3_leg$reporting, s3$reporting)
    expect_null(r3$resource_accounting$point_batch)
    expect_s3_class(r3$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
    for (f in c("support", "basis", "aggregations", "design_all", "design_used", "full_frame_transport")) {
      ru_expect_match(ru_strip_env(r3[[f]]), ru_strip_env(s3[[f]]), info = paste("incomplete date", f))
    }
    for (nm in c("estimates", "se", "analytical_se", "if_full", "vcov", "available", "ci_lower", "ci_upper", "p_value")) {
      ru_expect_match(r3$reporting[[nm]], s3$reporting[[nm]], info = paste("incomplete date reporting", nm))
    }
    expect_identical(r3$reporting$target_meta, s3$reporting$target_meta)
    expect_identical(ru_strip_env(r3$design_used), ru_strip_env(s3$design_used))
    expect_identical(r3$reporting$z_stat, r3$reporting$estimates / r3$reporting$se)
    expect_identical(s3$reporting$z_stat, s3$reporting$estimates / s3$reporting$se)
    ru_expect_match(r3$covariance_stream$mhat, s3$covariance_stream$mhat, info = "incomplete date mhat")
    ru_expect_match(r3$covariance_stream$V, s3$covariance_stream$V, info = "incomplete date V")
    expect_equal(unname(r3$covariance_stream$mean_phi), unname(s3$covariance_stream$mean_phi), tolerance = 1e-12)
    ru_expect_match(r3$gmm_fit$targets, s3$gmm_fit$targets, info = "incomplete date targets")
    for (nm in c("sandwich_cov", "materialized_target_if", "n_batches", "max_batch_rows_seen")) {
      ru_expect_match(r3$gmm_fit$if_projection[[nm]], s3$gmm_fit$if_projection[[nm]], info = paste("incomplete date if_projection", nm))
    }
    for (nm in c("n_links", "retained_index", "dropped_index", "E", "design_link_index", "link_kind", "classes")) {
      ru_expect_match(ru_strip_env(r3$structural_certificate[[nm]]), ru_strip_env(s3$structural_certificate[[nm]]),
        info = paste("incomplete date structural_certificate", nm))
    }
    expect_identical(r3$covariance_handling$handler_used, s3$covariance_handling$handler_used)
    penalty_identity <- c("weighting_requested", "weighting_canonical", "handler_selected", "handler_invoked",
                          "condition_trigger", "ridge_relative", "target_availability")
    expect_identical(r3$covariance_handling[penalty_identity], s3$covariance_handling[penalty_identity])
    expect_equal(r3$covariance_handling$lambda, s3$covariance_handling$lambda, tolerance = 1e-12)
  }
})

test_that("IPW carrier slices equal the plan chunks, all four IF components and the general replay on partial chunks; extra-cohort zero-observer chunks stream", {
  d <- ru_ipw_panel() # complete eligible owner; 41 = 5 * 7 + 6
  p <- 2L
  ctrl <- ru_ipw_ctrl(d, p, chunk = 7L)
  res <- ru_ipw_run(api, d, ctrl, p, target_cohorts = 4)
  expect_identical(ru_ipw_record(res)$status, "retained")
  fv <- res$nuisance_registry$fitted_views[[1L]]
  links <- res$point_links_used
  pc <- ctrl$point_control; oc <- ctrl$owner_control
  po <- res$panel_owner
  plan <- api$ipw_point_if_batch_plan(links, pc, fv)
  L <- length(links)
  carrier <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
  expect_true(is.environment(carrier))
  expect_identical(carrier$kind, "ipw_point_if_retained")
  expect_identical(carrier$state, "reserved")
  expect_false(exists("phi", envir = carrier, inherits = FALSE)) # reserved, not allocated
  expect_identical(carrier$N, 41L); expect_identical(carrier$L, L); expect_identical(carrier$B, 7L)
  expect_identical(carrier$n_chunks, 6L)
  # References of the same build: the component plan replay and the general replay.
  bb <- ru_collect(api$replay_ipw_point_if_plan, plan = plan, control = pc, include_components = TRUE)
  gb <- ru_collect(api$replay_ipw_point_if_view, point_links = links, control = pc, include_components = TRUE)
  cov_p <- api$stream_ipw_point_covariance_plan(plan, pc)
  cov_g <- api$stream_ipw_point_covariance_view(links, pc)
  cov_r <- api$stream_ipw_point_covariance_plan(plan, pc, carrier = carrier)
  for (nm in c("covariance", "mean_phi", "estimates", "n_a_units")) {
    expect_identical(cov_r[[nm]], cov_p[[nm]], info = nm)
    ru_expect_match(cov_r[[nm]], cov_g[[nm]], info = paste("general", nm))
  }
  expect_identical(carrier$state, "filled")
  expect_identical(carrier$producer_passes, 1L)
  expect_identical(carrier$slice_passes, 1L)
  expect_identical(carrier$batches, 6L)
  expect_identical(cov_p$replay$kernel_calls, 6L)
  expect_identical(cov_r$replay$kernel_calls, 0L)
  expect_identical(cov_r$replay$served_from, "retained_point_if")
  expect_identical(cov_r$replay$n_batches, cov_p$replay$n_batches)
  # Slices: bitwise the plan chunks, their four components summed in the kernel's order, the general values.
  sl <- ru_collect(api$replay_ipw_point_if_retained, carrier = carrier, plan = plan, control = pc)
  expect_length(sl, 6L)
  for (i in seq_along(sl)) {
    expect_identical(sl[[i]]$phi, bb[[i]]$phi, info = i)
    expect_identical(sl[[i]]$id, bb[[i]]$id, info = i)
    expect_identical(sl[[i]]$a_index, bb[[i]]$a_index, info = i)
    expect_identical(sl[[i]]$link_ids, bb[[i]]$link_ids, info = i)
    expect_identical(colnames(sl[[i]]$phi), plan$link_ids, info = i)
    expect_identical(sl[[i]]$phi, ((bb[[i]]$components$treated + bb[[i]]$components$control) +
      bb[[i]]$components$beta_score) + bb[[i]]$components$rate, info = paste("four-term order", i))
    ru_expect_match(sl[[i]]$phi, gb[[i]]$phi, info = paste("general phi", i))
    expect_null(sl[[i]]$components)
    expect_identical(names(sl[[i]]), c("phi", "id", "link_ids", "a_index"))
  }
  expect_identical(nrow(sl[[6L]]$phi), 6L)
  expect_identical(do.call(c, lapply(sl, `[[`, "id")), po$a_id)
  expect_identical(unname(carrier$phi), unname(do.call(rbind, lapply(bb, `[[`, "phi"))))
  expect_identical(colnames(carrier$phi), plan$link_ids)
  expect_gt(max(abs(colMeans(carrier$phi))), 0) # raw, uncentered
  expect_identical(carrier$slice_passes, 2L)
  # Every retained term matters: dropping any component from the stored slice breaks the general agreement.
  for (nm in c("rate", "beta_score", "treated", "control")) {
    diffs <- vapply(seq_along(sl), function(i) max(abs((sl[[i]]$phi - bb[[i]]$components[[nm]]) - gb[[i]]$phi)), numeric(1))
    if (nm == "rate" && all(vapply(bb, function(b) all(b$components$rate == 0), logical(1)))) next
    expect_gt(max(diffs), 1e-8, label = nm)
  }
  # The original zero-observer fixture contains another actual treated cohort;
  # its batch plan and retained carrier still equal the general four-term replay.
  dz <- ru_ipw_panel(extra = 15:21)
  cz <- ru_ipw_ctrl(dz, p, chunk = 7L)
  rz <- ru_ipw_run(api, dz, cz, p, target_cohorts = 4)
  expect_identical(ru_ipw_record(rz)$status, "retained")
  pz <- api$ipw_point_if_batch_plan(rz$point_links_used, cz$point_control, rz$nuisance_registry$fitted_views[[1L]])
  expect_true(is.environment(api$ipw_point_if_retained_reserve(pz, cz$point_control, live = 0, work = 0, envelope_control = cz$owner_control)))
  bz <- ru_collect(api$replay_ipw_point_if_plan, plan = pz, control = cz$point_control, include_components = TRUE)
  gz <- ru_collect(api$replay_ipw_point_if_view, point_links = rz$point_links_used, control = cz$point_control, include_components = TRUE)
  expect_true(all(bz[[3L]]$phi == 0))
  for (i in seq_along(bz)) {
    expect_identical(bz[[i]]$phi, ((bz[[i]]$components$treated + bz[[i]]$components$control) +
      bz[[i]]$components$beta_score) + bz[[i]]$components$rate)
    ru_expect_match(bz[[i]]$phi, gz[[i]]$phi, info = paste("zero-observer general", i))
    expect_identical(bz[[i]]$id, gz[[i]]$id)
    expect_identical(bz[[i]]$a_index, gz[[i]]$a_index)
  }
  # A second reservation on the same plan is independent; released carriers serve nothing.
  rec <- api$ipw_point_if_retained_release(carrier)
  expect_identical(rec$status, "retained")
  expect_identical(rec$producer_passes, 1L)
  expect_identical(rec$slice_passes, 2L)
  expect_identical(rec$n_batches, 6L)
  expect_identical(rec$reserved_entries, api$ipw_point_if_retained_inventory(plan)$entries)
  expect_false(exists("phi", envir = carrier, inherits = FALSE))
  expect_false(exists("plan_snapshot", envir = carrier, inherits = FALSE))
  expect_identical(carrier$state, "released")
  expect_error(api$replay_ipw_point_if_retained(carrier, plan, pc, function(b) NULL), class = "cdid_ipw_point_view_batch_error")
  expect_identical(api$ipw_point_if_retained_release(NULL), list(status = "streamed"))
})

test_that("IPW carrier: the GMM core's exact replay, ID and weight checks run unchanged on carrier slices", {
  d <- ru_ipw_panel()
  p <- 2L
  ctrl <- ru_ipw_ctrl(d, p, chunk = 7L)
  res <- ru_ipw_run(api, d, ctrl, p)
  fv <- res$nuisance_registry$fitted_views[[1L]]; links <- res$point_links_used
  pc <- ctrl$point_control; oc <- ctrl$owner_control
  plan <- api$ipw_point_if_batch_plan(links, pc, fv)
  po <- res$panel_owner; n_a <- po$n_a_units
  link_ids <- res$design_used$incidence$link_ids
  nuis_id <- "retained_point_if_test_nuis"
  prov <- api$gmm_covariance_provenance(
    row_level = "unit", n_units = n_a, n_clusters = n_a, n_positive_contributing_rows = n_a,
    centered = TRUE, weighted = FALSE, weight_type = "unweighted", normalization_factor = 1 / as.double(n_a),
    scale_id = "retained_point_if_test_scale", nuisance_fit_id = nuis_id,
    reduction = list(mode = "none", reduction_active = FALSE))
  prov$complete_unit_if_gram <- api$ipw_seal(list(covariance = res$covariance_stream$V, n_units = n_a,
    nuisance_fit_id = nuis_id), "complete_unit_if_gram")
  solve_with <- function(callback, expected_ids = po$a_id) api$solve_chained_gmm_core(
    design = res$design_used, mhat = res$covariance_stream$mhat, V = res$covariance_stream$V, provenance = prov,
    targets = api$ipw_pipeline_core_target_plan(res$support$ledger, res$design_used)$target_specs,
    weighting = "optimal",
    if_replay = api$gmm_if_replay_spec(callback = callback, input_row_level = "unit", aggregate_to_clusters = FALSE,
      n_input_rows = n_a, n_output_rows = n_a, n_positive_contributing_rows = n_a, centered = TRUE,
      weighted = FALSE, weight_type = "unweighted", normalization_factor = 1 / as.double(n_a),
      scale_id = prov$scale_id, nuisance_fit_id = nuis_id, expected_unit_ids = expected_ids),
    accuracy = res$controls$pipeline_control$gmm_accuracy, bounds = res$controls$pipeline_control$gmm_bounds,
    gauge_rule = "lowest", request_param_cov = TRUE, materialize_target_if = TRUE)
  emit_from <- function(consume, replay, tamper = NULL) replay(function(bt) {
    phi <- bt$phi; colnames(phi) <- link_ids
    if (!is.null(tamper)) phi <- tamper(phi)
    consume(list(unit_id = bt$id, phi = phi, weight = rep(1, nrow(phi)), nuisance_fit_id = nuis_id))
  })
  plan_replay <- function(cb) api$replay_ipw_point_if_plan(plan, pc, cb)
  filled_carrier <- function() {
    c0 <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
    api$ipw_point_if_retained_fill(c0, plan, pc, function(b) NULL)
    c0
  }
  g_plan <- solve_with(function(consume) emit_from(consume, plan_replay))
  c1 <- filled_carrier()
  g_ret <- solve_with(function(consume) emit_from(consume, function(cb) api$replay_ipw_point_if_retained(c1, plan, pc, cb)))
  expect_identical(c1$producer_passes, 1L)
  expect_identical(c1$slice_passes, 2L)
  expect_equal(g_ret$targets, g_plan$targets, tolerance = 1e-12)
  expect_equal(g_ret$if_projection$sandwich_cov, g_plan$if_projection$sandwich_cov, tolerance = 1e-12)
  expect_equal(g_ret$if_projection$materialized_target_if, g_plan$if_projection$materialized_target_if, tolerance = 1e-12)
  expect_equal(g_ret$targets$estimate, res$gmm_fit$targets$estimate, tolerance = 1e-12)
  expect_identical(g_ret$if_projection$n_batches, g_plan$if_projection$n_batches)
  expect_identical(g_ret$if_projection$n_batches, as.integer(ceiling(n_a / 7)))
  # A second pass whose projected rows differ is refused by the core's exact replay check.
  c2 <- filled_carrier(); pass <- 0L
  expect_error(solve_with(function(consume) {
    pass <<- pass + 1L
    emit_from(consume, function(cb) api$replay_ipw_point_if_retained(c2, plan, pc, cb),
      tamper = if (pass == 2L) function(phi) { phi[1L, 1L] <- phi[1L, 1L] + 1; phi } else NULL)
  }), class = "cdid_gmm_replay_error")
  expect_identical(pass, 2L)
  # A changed expected ID sequence and a changed weight are refused on slices.
  c3 <- filled_carrier()
  expect_error(solve_with(function(consume) emit_from(consume, function(cb) api$replay_ipw_point_if_retained(c3, plan, pc, cb)),
    expected_ids = rev(po$a_id)), class = "cdid_gmm_replay_error")
  c4 <- filled_carrier()
  expect_error(solve_with(function(consume) api$replay_ipw_point_if_retained(c4, plan, pc, function(bt) {
    phi <- bt$phi; colnames(phi) <- link_ids
    consume(list(unit_id = bt$id, phi = phi, weight = c(2, rep(1, nrow(phi) - 1L)), nuisance_fit_id = nuis_id))
  })), class = "cdid_gmm_replay_error")
})

test_that("IPW carrier admission is exact at the adjacent caps, declines ineligible facts before allocation and refuses every mutant at the boundary", {
  d <- ru_ipw_panel()
  p <- 2L
  ctrl <- ru_ipw_ctrl(d, p, chunk = 7L)
  res <- ru_ipw_run(api, d, ctrl, p)
  fv <- res$nuisance_registry$fitted_views[[1L]]; links <- res$point_links_used
  pc <- ctrl$point_control; oc <- ctrl$owner_control
  plan <- api$ipw_point_if_batch_plan(links, pc, fv)
  N <- 41L; L <- length(links); B <- 7L
  raw_entries <- api$ipw_point_if_retained_entries(N, L, B)
  expect_identical(raw_entries, 2 * N * L + 2 * N + 4 * B * L + 256)
  expect_identical(api$ipw_point_if_retained_work(N, L), 128 * N * L)
  expect_identical(api$ipw_point_if_retained_entries(10000, 19, 171), 2 * 10000 * 19 + 2 * 10000 + 4 * 171 * 19 + 256)
  inv <- api$ipw_point_if_retained_inventory(plan)
  # Source-bound batch inventory (R/chained_ipw_point_view.R:5786-5835).
  # Count the actual compact control payload independently; expected admission
  # caps must not be read back from the inventory under test.
  payload_counts <- function(x, scan_attrs = TRUE) {
    if (is.null(x) || is.environment(x)) return(c(entries = 0, work = 0))
    out <- c(entries = length(x) + 8, work = 64 * (length(x) + 8))
    if (scan_attrs && !is.null(attributes(x))) out <- out + payload_counts(attributes(x), FALSE)
    if (is.list(x)) for (el in x) out <- out + payload_counts(el)
    if (is.character(x)) {
      nb <- nchar(x, type = "bytes", allowNA = TRUE); nb[is.na(nb)] <- 0
      out <- out + c(entries = sum(ceiling(nb / 2)), work = 16 * sum(nb))
    }
    out
  }
  controls <- payload_counts(list(plan$control_in, plan$control)) + c(entries = 0, work = 128)
  label_bytes <- nchar(plan$link_ids, type = "bytes")
  expect_false(identical(plan$panel_owner$id_type, "character")) # no character-ID scratch/work in this fixture
  C <- ceiling(N / B)
  checks <- 16 * (C + 1)
  metadata <- 2 * unname(controls["entries"]) + ceiling(sum(label_bytes) / 2) + 32 + 128 * L + 4096
  covariance <- 4 * L * L + 8 * L + 6 * B * L + 128
  boundary_work <- checks * (64 * L + 512 + 2 * unname(controls["work"]) + 128 * L + 2048)
  entries <- raw_entries + metadata + covariance
  work_bound <- 128 * N * L + C * N * L + boundary_work + unname(controls["work"]) +
    checks * (128 + sum(label_bytes + 8))
  expect_identical(inv$metadata_entries, metadata)
  expect_identical(inv$covariance_entries, covariance)
  expect_identical(inv$boundary_work, boundary_work)
  expect_identical(inv$entries, entries)
  expect_identical(inv$work, work_bound)
  cap <- min(oc$max_retained_entries, oc$max_matrix_entries)
  reserve <- function(..., plan_ = plan, control = pc, live = 0, work = 0, envelope = oc)
    api$ipw_point_if_retained_reserve(plan_, control, live = live, work = work, envelope_control = envelope, ...)
  # Envelope cap: admitted exactly at the cap, declined one above; nothing allocated.
  at <- reserve(live = cap - entries)
  expect_true(is.environment(at))
  expect_identical(at$admission_total, cap)
  expect_identical(at$cap, cap)
  expect_identical(at$entries, entries)
  expect_false(exists("phi", envir = at, inherits = FALSE))
  expect_null(reserve(live = cap - entries + 1))
  # Finite private work cap: ordinary capacity refusal stays below the exact
  # binary64 integer ceiling. Exact-integer overflow remains a separate error.
  work_envelope <- oc
  work_envelope$max_work_ops <- 1e10
  expect_true(is.environment(reserve(work = work_envelope$max_work_ops - work_bound, envelope = work_envelope)))
  expect_null(reserve(work = work_envelope$max_work_ops - work_bound + 1, envelope = work_envelope))
  expect_error(reserve(work = (2^53 - 1) - work_bound + 1),
    "overflows exact IEEE binary64 integer range", class = "cdid_ipw_owner_bound_error")
  # The plan's own covariance-stream charge with the carrier: one below declines, at admits; the
  # plan's original stream at one below still runs (its original charge is unchanged).
  pre <- plan$preflight
  need <- pre$peak_concurrent_entries + plan$batch$additional_live_entries + 4 * L * L + 8 * L + 6 * pre$b_by_L_entries + 128 + entries
  tight_at <- ru_reseal(plan, "point_if_batch_plan", function(x) { x$control$max_matrix_entries <- need; x$control_in <- x$control; x })
  tight_below <- ru_reseal(plan, "point_if_batch_plan", function(x) { x$control$max_matrix_entries <- need - 1; x$control_in <- x$control; x })
  expect_true(is.environment(reserve(plan_ = tight_at, control = tight_at$control)))
  expect_null(reserve(plan_ = tight_below, control = tight_below$control))
  expect_false(inherits(tryCatch(api$stream_ipw_point_covariance_plan(tight_below, tight_below$control), error = identity), "error"))
  # Structural admission, transfer decline, count/input/plan refusals.
  expect_true(is.environment(reserve(clustered = TRUE)))
  expect_null(reserve(numerical_transfer = TRUE))
  expect_true(is.environment(reserve(n_total = N + 1L)))
  expect_error(reserve(n_total = N - 1L), class = "cdid_ipw_count_error")
  expect_true(is.environment(reserve(n_total = N)))
  expect_error(reserve(plan_ = NULL), class = "cdid_ipw_point_view_batch_error")
  expect_error(reserve(plan_ = list(kind = "other")), class = "cdid_ipw_point_view_batch_error")
  expect_error(reserve(live = -1), class = "cdid_ipw_point_view_input_error")
  expect_error(reserve(envelope = NULL), class = "cdid_ipw_point_view_input_error")
  # Integrity refusals of the plan verifier propagate (not declines): an edited unsealed plan, foreign controls.
  bad <- plan; bad$link_par[1L, 1L] <- 0
  expect_error(reserve(plan_ = bad), class = "cdid_ipw_point_view_batch_error")
  expect_error(reserve(control = api$ipw_point_view_control(max_chunk_rows = 8L, max_observer_units = 1e4,
    max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)),
    class = "cdid_ipw_control_mismatch_error")
  # Non-unit constant weights on a genuine plan are admitted under D92.
  fx_w <- {
    oc2 <- api$ipw_owner_control(max_chunk_rows = 7L, max_retained_entries = 1e7, max_matrix_entries = 1e7,
      max_work_ops = 1e10, max_output_entries = 1e7)
    dd <- rm_panel()
    po2 <- api$construct_ipw_panel_owner(dd, "id", "t", "y", "g", oc2)
    uu <- dd[!duplicated(dd$id), ]; uu <- uu[order(uu$id), ]
    fo2 <- api$construct_ipw_frame_owner(po2, uu$id, matrix(1, nrow(uu), 1L, dimnames = list(NULL, "(Intercept)")),
      weights = rep(1.5, nrow(uu)), control = oc2)
    nv2 <- api$construct_ipw_nuisance_view(fo2, 3, Inf, 1, 2, oc2)
    fc2 <- api$ipw_fitted_view_control(max_chunk_rows = 7L, max_observer_units = 1e4, max_retained_entries = 1e7,
      max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
    fv2 <- api$fit_ipw_nuisance_view(nv2, stats::glm.control(epsilon = 1e-12, maxit = 50), fc2)
    lvs2 <- lapply(list(c(1, 2), c(2, 3)), function(pr) api$construct_ipw_link_view(fo2, pr[1], pr[2], nuisance_view = nv2,
      control = oc2, link_id = paste0("W", pr[1], pr[2]), functional_key = "W"))
    pc2 <- api$ipw_point_view_control(max_chunk_rows = 7L, max_observer_units = 1e4, max_retained_entries = 1e7,
      max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
    pls <- lapply(lvs2, function(lv) api$fit_ipw_point_link_view(fv2, lv, pc2))
    list(plan = api$ipw_point_if_batch_plan(pls, pc2, fv2), pc = pc2, oc = oc2)
  }
  expect_true(isTRUE(all(fx_w$plan$frame_owner$w_a == 1.5)))
  expect_true(is.environment(api$ipw_point_if_retained_reserve(fx_w$plan, fx_w$pc, live = 0, work = 0, envelope_control = fx_w$oc)))

  # Lifecycle refusals and mutants on a filled carrier.
  c_unfilled <- reserve()
  expect_error(api$replay_ipw_point_if_retained(c_unfilled, plan, pc, function(b) NULL), class = "cdid_ipw_point_view_batch_error")
  expect_identical(c_unfilled$state, "released")
  carrier <- reserve()
  calls <- 0L
  api$ipw_point_if_retained_fill(carrier, plan, pc, function(b) calls <<- calls + 1L)
  expect_identical(calls, 6L)
  expect_identical(carrier$state, "filled")
  fresh_filled <- function() {
    c0 <- reserve()
    api$ipw_point_if_retained_fill(c0, plan, pc, function(b) NULL)
    c0
  }
  c_dup <- fresh_filled()
  expect_error(api$ipw_point_if_retained_fill(c_dup, plan, pc, function(b) NULL), class = "cdid_ipw_point_view_batch_error")
  expect_identical(c_dup$state, "released")
  expect_error(api$stream_ipw_point_covariance_plan(plan, pc, carrier = fresh_filled()), class = "cdid_ipw_point_view_batch_error")
  refused <- function(plan_ = plan, control = pc, class = "cdid_ipw_point_view_batch_error", carrier_ = fresh_filled()) {
    n <- 0L
    expect_error(api$replay_ipw_point_if_retained(carrier_, plan_, control, function(b) { force(b); n <<- n + 1L }), class = class)
    n
  }
  mutate_phi <- function(fn) {
    cm <- fresh_filled()
    unlockBinding("phi", cm)
    cm$phi <- fn(cm$phi)
    lockBinding("phi", cm)
    cm
  }
  orig <- carrier$phi
  expect_identical(refused(carrier_ = mutate_phi(function(m) { m[3L, 2L] <- m[3L, 2L] + 1e-9; m })), 0L)
  expect_identical(refused(carrier_ = mutate_phi(function(m) { m[c(1L, 2L), ] <- m[c(2L, 1L), ]; m })), 0L)
  expect_identical(refused(carrier_ = mutate_phi(function(m) m[, rev(seq_len(L)), drop = FALSE])), 0L)
  expect_identical(refused(carrier_ = mutate_phi(function(m) m[seq_len(40L), , drop = FALSE])), 0L)
  expect_identical(refused(carrier_ = mutate_phi(function(m) {
    m[1:2, 1:2] <- m[1:2, 1:2] + matrix(c(1, -1, -1, 1), 2L); m
  })), 0L)
  expect_identical(refused(carrier_ = mutate_phi(function(m) { colnames(m) <- rev(colnames(m)); m })), 0L)
  expect_length(ru_collect(api$replay_ipw_point_if_retained, carrier = carrier, plan = plan, control = pc), 6L) # untouched carrier serves
  # Plan identity: a resealed copy, another plan of the same links, or foreign controls refuse.
  resealed <- ru_reseal(plan, "point_if_batch_plan", identity)
  expect_identical(refused(plan_ = resealed, control = resealed$control), 0L)
  plan2 <- api$ipw_point_if_batch_plan(links, pc, fv)
  expect_identical(refused(plan_ = plan2), 0L)
  expect_identical(refused(control = api$ipw_point_view_control(max_chunk_rows = 8L, max_observer_units = 1e4,
    max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7),
    class = "cdid_ipw_control_mismatch_error"), 0L)
  expect_identical(refused(carrier_ = new.env()), 0L)
  # Owner/nuisance/link/point-link rebinding behind the sealed plan refuses before any slice.
  fo <- res$frame_owner; fo2 <- fo; fo2$sum_w_a <- fo2$sum_w_a + 1
  nv <- res$nuisance_registry$nuisance_views[[1L]]
  for (sw in list(list(env = fo$ref, nm = "frame_owner", value = fo2),
                  list(env = fv$ref, nm = "fitted_nuisance_view", value = unclass(fv)),
                  list(env = nv$ref, nm = "nuisance_view", value = unclass(nv)),
                  list(env = links[[1L]]$ref, nm = "point_link_view", value = unclass(links[[1L]])))) {
    restore <- ru_rebind(sw$env, sw$nm, sw$value)
    expect_identical(refused(), 0L, info = sw$nm)
    restore()
  }
  expect_length(ru_collect(api$replay_ipw_point_if_retained, carrier = carrier, plan = plan, control = pc), 6L)
  # Callback rebinding: the next slice is refused before it is emitted.
  c_cb <- fresh_filled()
  restore <- NULL; n <- 0L
  expect_error(api$replay_ipw_point_if_retained(c_cb, plan, pc, function(b) {
    force(b); n <<- n + 1L
    if (is.null(restore)) restore <<- ru_rebind(nv$ref, "nuisance_view", unclass(nv))
  }), class = "cdid_ipw_point_view_batch_error")
  restore()
  expect_identical(n, 1L)
  expect_identical(c_cb$state, "released")
  # The final callback is a boundary too: raw replacement cannot escape when
  # no later slice is due.
  c_last <- fresh_filled()
  n <- 0L
  expect_error(api$replay_ipw_point_if_retained(c_last, plan, pc, function(b) {
    n <<- n + 1L
    if (n == c_last$n_chunks) {
      unlockBinding("phi", c_last)
      c_last$phi <- orig + 1
      lockBinding("phi", c_last)
    }
  }), class = "cdid_ipw_point_view_batch_error")
  expect_identical(n, carrier$n_chunks)
  expect_identical(c_last$state, "released")
  # A producer that fails is never caught: a filling carrier under a producer refusal leaves no filled store
  # and the refusal propagates with its own class.
  e <- ru_ns()
  e$ipw_point_if_plan_chunk <- function(plan, a_start, a_end, kernel, include_components = FALSE) {
    out <- api$ipw_point_if_plan_chunk(plan, a_start, a_end, kernel, include_components)
    if (a_start > 7) api$cdid_abort("forced producer refusal", "cdid_ipw_numeric_error")
    out
  }
  c_fail <- e$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
  err <- tryCatch(e$stream_ipw_point_covariance_plan(plan, pc, carrier = c_fail), error = identity)
  expect_s3_class(err, "cdid_ipw_numeric_error")
  expect_identical(c_fail$state, "released")
  expect_false(exists("raw_ref", envir = c_fail, inherits = FALSE))
  expect_false(exists("phi", envir = c_fail, inherits = FALSE))
  rec_fail <- api$ipw_point_if_retained_release(c_fail)
  expect_identical(rec_fail$status, "incomplete:released")
})

# ===========================================================================
# PART 2. OR/DR complete-support route: auxiliary consumer of the CP04 store
# ===========================================================================

test_that("OR/DR: exact same-recipe auxiliary columns are served from the retained store, unmatched ones stream; outputs equal the forced stream at 1e-12", {
  for (dat in list(list(d = rm_d, nm = "g3"), list(d = rm_d5, nm = "g5"))) for (est in c("reg", "dr")) for (xf in list(NULL, ~x)) for (bp in c("varying", "universal")) {
    info <- paste(dat$nm, est, if (is.null(xf)) "noX" else "~x", bp)
    e_r <- rm_counting_ns()
    res_r <- rm_run_pipeline(dat$d, est = est, chunk = 7L, xformla = xf, base_period = bp, aggregations = "simple", ns = e_r)
    e_s <- rm_counting_ns(force_stream = TRUE)
    res_s <- rm_run_pipeline(dat$d, est = est, chunk = 7L, xformla = xf, base_period = bp, aggregations = "simple", ns = e_s)
    rec <- rm_record(res_r)
    expect_identical(rec$status, "retained", info = info)
    expect_identical(rec$producer_passes, 1L, info = info)
    expect_identical(rec$slice_passes, 3L, info = info)
    expect_identical(rm_record(res_s)$status, "streamed", info = info)
    expect_identical(rm_record(res_s)$auxiliary, NULL, info = info)
    n_chunks <- as.integer(ceiling(41 / 7))
    expect_identical(e_r$n$replays, 4L, info = info)
    expect_identical(e_r$n$chunks, n_chunks, info = info)
    expect_identical(e_s$n$chunks, 4L * n_chunks, info = info)
    # Independent structural expectation: servable iff the descriptor's pair is a retained catalog pair.
    exp_map <- rm_expected_served(res_r, bp)
    servable <- exp_map$servable; backend_ids <- exp_map$backend_ids
    expect_identical(exp_map, rm_expected_served(res_s, bp), info = info)
    expect_gt(length(servable), 0L)
    expect_true(all(backend_ids > res_r$reduction_summary$n_links_all), info = info) # auxiliary IDs are M + j
    expect_identical(e_s$n$aux_ids[[1L]], unique(backend_ids), info = info)
    aux <- rec$auxiliary
    if (any(servable)) {
      expect_identical(aux$status, "extended", info = info)
      expect_identical(aux$columns_requested, length(servable), info = info)
      expect_identical(aux$columns_served, sum(servable), info = info)
      expect_identical(aux$columns_streamed, sum(!servable), info = info)
      expect_identical(aux$slice_passes, 1L, info = info)
      expect_identical(e_r$n$aux_slices, 1L, info = info)
      expect_match(rec$lifetime, "auxiliary columns served", info = info)
      if (any(!servable)) {
        # The unmatched remainder streams once, in its original order, under its own issued partition.
        expect_identical(e_r$n$aux_replays, 1L, info = info)
        expect_identical(e_r$n$aux_chunks, as.integer(ceiling(41 / aux$stream_chunk_rows)), info = info)
        expect_identical(e_r$n$aux_ids[[1L]], unique(backend_ids[!servable]), info = info)
      } else {
        # Every requested column served: zero auxiliary producer pass.
        expect_identical(e_r$n$aux_replays, 0L, info = info)
        expect_identical(e_r$n$aux_chunks, 0L, info = info)
        expect_true(is.na(aux$stream_chunk_rows), info = info)
      }
    } else {
      expect_identical(aux$status, "not_extended", info = info)
      expect_identical(e_r$n$aux_replays, 1L, info = info)
      expect_identical(e_r$n$aux_slices, 0L, info = info)
    }
    expect_identical(e_s$n$aux_replays, 1L, info = info)
    expect_identical(e_s$n$aux_slices, 0L, info = info)
    # Signed values into the original slots: served pretrend columns equal the streamed auxiliary columns
    # (the same recipe on the same sealed backend; lambda = 1 frame identity) at the established 1e-12,
    # and every reported field agrees at 1e-12 with identical labels, roles and auxiliary indices.
    pre_cols <- which(res_r$reporting$target_meta$role == "pre" & res_r$reporting$available)
    expect_gt(length(pre_cols), 0L)
    expect_equal(res_r$reporting$if_full[, pre_cols, drop = FALSE], res_s$reporting$if_full[, pre_cols, drop = FALSE],
      tolerance = 1e-12, info = info)
    expect_identical(colnames(res_r$reporting$if_full), colnames(res_s$reporting$if_full), info = info)
    expect_identical(res_r$reporting$target_meta$auxiliary_index, res_s$reporting$target_meta$auxiliary_index, info = info)
    expect_identical(res_r$reporting$target_meta$role, res_s$reporting$target_meta$role, info = info)
    rm_expect_parity(res_r, res_s, info = info)
    expect_false(ru_find_retained(res_r), info = info)
    expect_identical(res_r$nuisance_registry$cp04_reader$state$lambda, 1, info = info)
  }
})

test_that("OR/DR auxiliary extension: the ridge-active retained metric, the declined extension and the adjacent live cap keep the stream's values", {
  d_ridge <- rm_nonadd_panel(12L, 5L, 6L, seed = 55L)
  for (est in c("reg", "dr")) {
    res_r <- rm_run_pipeline(d_ridge, est = est, chunk = 7L, condition_trigger = 0.99)
    res_s <- rm_run_pipeline(d_ridge, est = est, chunk = 7L, condition_trigger = 0.99, ns = rm_stream_ns())
    expect_identical(rm_record(res_r)$status, "retained")
    expect_identical(res_r$covariance_handling$handler_used, "tikhonov")
    expect_equal(res_r$covariance_handling$lambda, res_s$covariance_handling$lambda, tolerance = 1e-12)
    rm_expect_parity(res_r, res_s, info = paste("ridge", est))
  }
  for (est in c("reg", "dr")) for (bp in c("varying", "universal")) {
    info <- paste(est, bp)
    res <- rm_run_pipeline(rm_d5, est = est, chunk = 7L, base_period = bp, aggregations = "simple")
    ref <- rm_run_pipeline(rm_d5, est = est, chunk = 7L, base_period = bp, aggregations = "simple", ns = rm_stream_ns())
    rec <- rm_record(res)
    expect_identical(rec$status, "retained", info = info)
    # Declined extension: the store keeps its original lifetime and record; values unchanged.
    res_no <- rm_run_pipeline(rm_d5, est = est, chunk = 7L, base_period = bp, aggregations = "simple", ns = rm_no_aux_ns())
    expect_identical(rm_record(res_no)$status, "retained", info = info)
    expect_identical(rm_record(res_no)$auxiliary$status, "not_extended", info = info)
    expect_identical(rm_record(res_no)$slice_passes, 3L, info = info)
    expect_identical(rm_record(res_no)[c("reserved_entries", "admission_total_entries", "cap_entries", "work_ops_bound")],
      rec[c("reserved_entries", "admission_total_entries", "cap_entries", "work_ops_bound")], info = info)
    rm_expect_parity(res_no, ref, info = paste(info, "declined"))
    if (identical(rec$auxiliary$status, "extended")) {
      # Live envelope at the store's cap admits; one above declines the extension only.
      res_at <- rm_run_pipeline(rm_d5, est = est, chunk = 7L, base_period = bp, aggregations = "simple", ns = rm_aux_live_ns(0))
      res_above <- rm_run_pipeline(rm_d5, est = est, chunk = 7L, base_period = bp, aggregations = "simple", ns = rm_aux_live_ns(1))
      expect_identical(rm_record(res_at)$auxiliary$status, "extended", info = info)
      expect_identical(rm_record(res_at)$auxiliary$live_total_entries, rec$cap_entries, info = info)
      expect_identical(rm_record(res_above)$status, "retained", info = info)
      expect_identical(rm_record(res_above)$auxiliary$status, "not_extended", info = info)
      rm_expect_parity(res_at, ref, info = paste(info, "at"))
      rm_expect_parity(res_above, ref, info = paste(info, "above"))
      # The extension's live total exceeds the main admission by the transported target IF charge.
      J <- length(res$gmm_fit$targets$estimable[res$gmm_fit$targets$estimable])
      expect_identical(rec$auxiliary$live_total_entries, rec$admission_total_entries + 2 * 41 * J, info = info)
      expect_identical(rec$auxiliary$work_ops_bound,
        api$cp04_retained_if_auxiliary_work(41, rec$n_links) + (41 + 4) * rec$boundary_check_work_ops, info = info)
    }
    # Point estimates only: no auxiliary IF consumer, no extension, no IF.
    res_pt <- rm_run_pipeline(rm_d5, est = est, chunk = 7L, base_period = bp, aggregations = character(0), compute_inffunc = FALSE)
    expect_identical(rm_record(res_pt)$status, "retained", info = info)
    expect_identical(rm_record(res_pt)$auxiliary$status, "not_extended", info = info)
    expect_null(res_pt$reporting$if_full, info = info)
    expect_equal(res_pt$reporting$estimates, res$reporting$estimates, tolerance = 1e-12, info = info)
  }
  # D92 generalized store: clustered inference, constant non-unit weights, and a frame-only unit retain
  # the main store; auxiliary extension admits when N = n_A and lambda = 1 and declines when N > n_A.
  d0 <- rm_nonadd_panel(41L, 5L, 18L, seed = 31L)
  fr <- d0[!duplicated(d0$id), c("id", "G", "x", "cl")]
  fr <- rbind(fr[order(fr$id), ], data.frame(id = 9999L, G = 0, x = 0, cl = 99L))
  for (case in list(list(nm = "clustered", aux = "extended", run = function(ns) rm_run_pipeline(d0, est = "dr", clustered = TRUE, aggregations = "simple", ns = ns)),
                    list(nm = "constant_weight", aux = "extended", run = function(ns) rm_run_pipeline(d0, est = "reg", weightsname = "w2", aggregations = "simple", ns = ns)),
                    list(nm = "frame_only_unit", aux = "not_extended", run = function(ns) rm_run_pipeline(d0, est = "reg", frame = fr, aggregations = "simple", ns = ns)))) {
    res <- case$run(api); ref <- case$run(rm_stream_ns())
    expect_identical(rm_record(res)$status, "retained", info = case$nm)
    expect_identical(rm_record(res)$auxiliary$status, case$aux, info = case$nm)
    expect_false(ru_find_retained(res), info = case$nm)
    rm_expect_parity(res, ref, info = case$nm)
  }
  # Zero within-cohort variance refuses identically on both routes.
  d_zero <- rm_nonadd_panel(12L, 5L, 6L, seed = 55L)
  d_zero$y <- ifelse(d_zero$G == 3 & d_zero$period >= 3, 2 * (d_zero$period - 2), 0)
  for (est in c("reg", "dr")) {
    e_r <- tryCatch(rm_run_pipeline(d_zero, est = est, chunk = 7L), error = identity)
    e_s <- tryCatch(rm_run_pipeline(d_zero, est = est, chunk = 7L, ns = rm_stream_ns()), error = identity)
    expect_s3_class(e_r, "cdid_ipw_pipeline_scale_error")
    expect_identical(class(e_r), class(e_s))
    expect_identical(conditionMessage(e_r), conditionMessage(e_s))
  }
})

test_that("the auxiliary recipe map is exact: a repeated recipe matches, labels/indices/forged descriptors/references do not; slices are the frame-scale stream bitwise", {
  # Catalog pairs (all ten of five dates) plus appended views that repeat pairs 1-2 and 2-3 under their own labels.
  fx <- rm_fixture(aux_pairs = list(c(1, 2), c(2, 3)))
  st <- fx$reader$state
  expect_identical(st$lambda, 1)
  expect_identical(st$N, st$n_A)
  for (method in c("or", "dr")) {
    fit <- if (method == "or") rm_or(fx) else rm_dr(fx)
    links <- c(1L, 5L, 8L, 10L) # pairs (1,2), (2,3), (3,4), (4,5)
    for (k in seq_along(links)) expect_identical(as.double(fx$pairs[[links[k]]]), c(k, k + 1))
    aux_ids <- as.integer(fx$aux_index) # 11 = (1,2), 12 = (2,3)
    expect_identical(aux_ids, c(11L, 12L))
    expect_false(identical(fit$links[[11L]]$link_id, fit$links[[1L]]$link_id))
    plan <- api$cp04_stream_plan(fit, links, 0)
    B <- plan$chunk_rows
    reserve <- function() api$cp04_retained_if_reserve(fit, links, B, plan, retained = 0, gmm_retained = 0, live = 0, work = 0,
      control = fit$control)
    store <- reserve()
    expect_true(is.environment(store))
    # Recipes: the repeated views carry exactly the catalog recipes; every canonical field is present.
    fields <- c("regression", "pair_id", "cohort", "support_id", if (method == "dr") c("model", "control_support_ids"))
    r1 <- api$cp04_retained_if_recipe(fit, 1L); r11 <- api$cp04_retained_if_recipe(fit, 11L)
    expect_identical(names(r1), fields)
    expect_identical(r1, r11)
    expect_identical(api$cp04_retained_if_recipe(fit, 5L), api$cp04_retained_if_recipe(fit, 12L))
    expect_false(identical(api$cp04_retained_if_recipe(fit, 2L), r1))
    # Map: repeated recipes hit their retained columns; other catalog pairs and the retained indices themselves do not.
    expect_identical(api$cp04_retained_if_auxiliary_map(store, fit, c(11L, 12L, 2L, 3L)), c(1L, 2L, NA_integer_, NA_integer_))
    expect_identical(api$cp04_retained_if_auxiliary_map(store, fit, c(12L, 11L)), c(2L, 1L))
    expect_identical(api$cp04_retained_if_auxiliary_map(store, fit, 1L), NA_integer_)
    expect_identical(api$cp04_retained_if_auxiliary_map(store, fit, c(11L, 11L)), c(1L, 1L))
    expect_identical(api$cp04_retained_if_auxiliary_map(store, fit, integer(0)), integer(0))
    expect_identical(api$cp04_retained_if_auxiliary_map(store, fit, c(11L, 99L)), c(NA_integer_, NA_integer_))
    # Forged descriptors (an unsealed copy; the map is pure): any single canonical field, or a reference, breaks the match.
    for (f in fields) {
      forged <- fit
      forged$links[[11L]][[f]] <- forged$links[[11L]][[f]] + 1L
      expect_identical(api$cp04_retained_if_auxiliary_map(store, forged, 11L), NA_integer_, info = paste(method, f))
      # DR preserves duplicate model/control-support names from the original
      # link record. Remove every exact-name occurrence to witness absence;
      # deleting the active first field alone exposes a valid second field.
      forged$links[[11L]] <- forged$links[[11L]][names(forged$links[[11L]]) != f]
      expect_identical(api$cp04_retained_if_auxiliary_map(store, forged, 11L), NA_integer_, info = paste(method, f, "missing"))
    }
    forged <- fit; forged$links[[11L]]$reference <- 1L
    expect_identical(api$cp04_retained_if_auxiliary_map(store, forged, 11L), NA_integer_)
    forged <- fit; forged$links[[11L]]$cohort <- forged$links[[11L]]$cohort + 0.5
    expect_identical(api$cp04_retained_if_auxiliary_map(store, forged, 11L), NA_integer_)
    # A label change alone never matters either way.
    forged <- fit; forged$links[[11L]]$link_id <- "renamed"
    expect_identical(api$cp04_retained_if_auxiliary_map(store, forged, 11L), 1L)

    # Admission before allocation.
    admit <- function(ids, store_ = store, live = store_$entries, retained = 0, work = 0, backend = fit, control = fit$control)
      api$cp04_retained_if_auxiliary_admit(store_, backend, ids, retained = retained, live = live, work = work, control = control)
    aux_work_bound <- api$cp04_retained_if_auxiliary_work(41, 4) + (41 + 4) * store$boundary_work
    expect_null(admit(c(2L, 3L)))                       # nothing matches
    expect_null(admit(c(11L, 1L)))                      # a retained index is requested
    expect_error(admit(11L, live = store$entries - 1), class = "cdid_cp04_integrity_error") # the envelope must contain the store
    expect_null(admit(11L, live = store$cap + 1))       # one above the cap
    expect_null(admit(11L, work = fit$control$max_work_ops - aux_work_bound + 1))
    expect_true(is.environment(admit(11L, store_ = reserve(), work = fit$control$max_work_ops - aux_work_bound)))
    expect_true(is.environment(admit(11L, store_ = reserve(), live = store$cap)))
    expect_error(admit(11L, live = 0.5), class = "cdid_ipw_owner_bound_error")
    expect_null(store$auxiliary)
    expect_false(exists("phi", envir = store, inherits = FALSE))
    expect_true(is.environment(admit(c(11L, 12L, 2L))))
    aux <- store$auxiliary
    expect_identical(aux$backend_ids, c(11L, 12L, 2L))
    expect_identical(aux$column, c(1L, 2L, NA_integer_))
    expect_identical(aux$ids_stream, 2L)
    expect_identical(aux$columns_served, 2L); expect_identical(aux$columns_streamed, 1L)
    expect_identical(aux$work, aux_work_bound)
    expect_identical(aux$chunk_rows, as.double(api$cp04_preflight_links(fit, c(11L, 12L, 2L), 0)$chunk_rows))
    expect_identical(aux$chunk_rows_stream, aux$chunk_rows)
    expect_identical(aux$chunk_rows_stream, as.double(api$cp04_preflight_links(fit, 2L, store$entries)$chunk_rows))
    expect_null(admit(11L))                              # already admitted
    expect_false(exists("phi", envir = store, inherits = FALSE)) # still unallocated
    # A filled store cannot take the extension afterwards.
    late <- reserve(); bound_late <- api$cp04_retained_if_bind(fit, late)
    invisible(rm_collect(bound_late, links, B))
    expect_identical(late$state, "filled")
    expect_null(admit(11L, store_ = late))

    # Fill through the bound producer, then serve: the stored catalog columns bitwise (they ARE the
    # frame-scale values: lambda = 1), the ordinary frame-scale replay of the repeated views at 1e-12,
    # the live IDs, the frame normalization count, the store's own row partition.
    bound <- api$cp04_retained_if_bind(fit, store)
    invisible(rm_collect(bound, links, B))
    expect_identical(store$producer_passes, 1L)
    frame_batches <- rm_collect(fit, c(11L, 12L), B, scale = "frame")
    frame_ref <- do.call(rbind, lapply(frame_batches, `[[`, "phi"))
    frame_cat <- do.call(rbind, lapply(rm_collect(fit, c(1L, 5L), B, scale = "frame"), `[[`, "phi"))
    served <- list()
    api$cp04_retained_if_auxiliary_slices(store, fit, c(11L, 12L, 2L), function(b) { force(b); served[[length(served) + 1L]] <<- b })
    expect_length(served, length(frame_batches))
    got <- do.call(rbind, lapply(served, `[[`, "phi"))
    expect_identical(unname(got), unname(store$phi[, 1:2]))
    expect_identical(unname(got), unname(frame_cat))              # core = 1 * frame, bitwise
    expect_equal(unname(got), unname(frame_ref), tolerance = 1e-12) # the repeated views' own frame replay
    expect_identical(colnames(got), vapply(fit$links[c(11L, 12L)], `[[`, character(1), "link_id"))
    for (k in seq_along(served)) {
      b <- served[[k]]; fb <- frame_batches[[k]]
      expect_identical(b$rows, fb$rows)
      expect_identical(b$id, fb$id)
      expect_identical(b$id, st$id[b$rows])
      expect_identical(b$links, c(11L, 12L))
      expect_identical(b$scale, "frame")
      expect_identical(b$normalization_count, fb$normalization_count)
      expect_identical(b$normalization_count, st$N)
    }
    expect_identical(store$auxiliary$slice_passes, 1L)
    expect_identical(store$slice_passes, 0L)
    # The signed reporting arithmetic on a served column is the streamed arithmetic on the same values.
    expect_identical(-1 * got[, 1L], -1 * frame_cat[, 1L])
    expect_equal(-1 * got[, 1L], -1 * frame_ref[, 1L], tolerance = 1e-12)
    # Service: the identical request returns the admitted record; any difference is an integrity failure.
    expect_identical(api$cp04_retained_if_auxiliary_service(store, fit, c(11L, 12L, 2L), 0)$column, c(1L, 2L, NA_integer_))
    expect_error(api$cp04_retained_if_auxiliary_service(store, fit, c(12L, 11L, 2L), 0), class = "cdid_cp04_integrity_error")
    expect_error(api$cp04_retained_if_auxiliary_service(store, fit, c(11L, 12L), 0), class = "cdid_cp04_integrity_error")
    expect_error(api$cp04_retained_if_auxiliary_service(store, fit, c(11L, 12L, 2L), 1), class = "cdid_cp04_integrity_error")
    expect_null(api$cp04_retained_if_auxiliary_service(NULL, fit, c(11L, 12L, 2L), 0))
    expect_null(api$cp04_retained_if_auxiliary_service(reserve(), fit, c(11L, 12L, 2L), 0))
    fresh_aux_filled <- function() {
      s0 <- reserve()
      admit(c(11L, 12L, 2L), store_ = s0)
      invisible(rm_collect(api$cp04_retained_if_bind(fit, s0), links, B))
      s0
    }
    expect_error(api$cp04_retained_if_auxiliary_slices(fresh_aux_filled(), fit, c(12L, 11L, 2L), function(b) NULL),
      class = "cdid_cp04_integrity_error")
    # Mutants are refused before any auxiliary slice is emitted.
    refused <- function(backend = fit, ids = c(11L, 12L, 2L), store_ = fresh_aux_filled()) {
      n <- 0L
      expect_error(api$cp04_retained_if_auxiliary_slices(store_, backend, ids, function(b) { force(b); n <<- n + 1L }),
        class = "cdid_cp04_integrity_error")
      n
    }
    orig <- store$phi
    s_phi <- fresh_aux_filled(); s_phi$phi[2L, 1L] <- s_phi$phi[2L, 1L] + 1e-9
    expect_identical(refused(store_ = s_phi), 0L)
    expect_identical(s_phi$state, "released")
    s_lam <- fresh_aux_filled(); unlockBinding("lambda", s_lam); s_lam$lambda <- 2 * s_lam$lambda; lockBinding("lambda", s_lam)
    expect_identical(refused(store_ = s_lam), 0L)
    expect_identical(s_lam$state, "released")
    fit2 <- if (method == "or") rm_or(fx) else rm_dr(fx)
    expect_identical(refused(backend = fit2), 0L)
    served2 <- list()
    api$cp04_retained_if_auxiliary_slices(store, fit, c(11L, 12L, 2L), function(b) served2[[length(served2) + 1L]] <<- b)
    expect_identical(lapply(served2, `[[`, "phi"), lapply(served, `[[`, "phi"))
    # Callback rebinding of the panel owner refuses the next auxiliary slice through the original validator.
    err <- if (method == "or") "cdid_or_baseline_integrity_error" else "cdid_dr_baseline_integrity_error"
    po <- fx$po; shifted <- po; shifted$postings$date_start[2L] <- shifted$postings$date_start[2L] + 1L
    s_cb <- fresh_aux_filled()
    restore <- NULL; n <- 0L
    expect_error(api$cp04_retained_if_auxiliary_slices(s_cb, fit, c(11L, 12L, 2L), function(b) {
      force(b); n <<- n + 1L
      if (is.null(restore)) restore <<- ru_rebind(po$ref, "panel_owner", shifted)
    }), class = err)
    restore()
    expect_identical(n, 1L)
    expect_identical(s_cb$state, "released")
    # A raw replacement at the final auxiliary callback is also refused.
    s_last <- fresh_aux_filled()
    n <- 0L
    expect_error(api$cp04_retained_if_auxiliary_slices(s_last, fit, c(11L, 12L, 2L), function(b) {
      n <<- n + 1L
      if (max(b$rows) == st$N) s_last$phi <- orig + 1
    }), class = "cdid_cp04_integrity_error")
    expect_identical(n, as.integer(ceiling(st$N / store$auxiliary$chunk_rows)))
    expect_identical(s_last$state, "released")
    # Release records the extension and drops everything; nothing serves afterwards.
    rec <- api$cp04_retained_if_release(store)
    expect_identical(rec$status, "retained")
    expect_identical(rec$slice_passes, 0L)
    expect_identical(rec$auxiliary$status, "extended")
    expect_identical(rec$auxiliary$columns_served, 2L)
    expect_identical(rec$auxiliary$columns_streamed, 1L)
    expect_identical(rec$auxiliary$slice_passes, 2L)
    expect_identical(rec$reserved_entries, api$cp04_retained_if_entries(41L, 4L, B) + rec$snapshot_metadata_entries)
    expect_false(exists("phi", envir = store, inherits = FALSE))
    expect_false(exists("auxiliary", envir = store, inherits = FALSE))
    expect_identical(store$state, "released")
    expect_error(api$cp04_retained_if_auxiliary_slices(store, fit, c(11L, 12L, 2L), function(b) NULL), class = "cdid_cp04_integrity_error")
    expect_error(api$cp04_retained_if_auxiliary_service(store, fit, c(11L, 12L, 2L), 0), class = "cdid_cp04_integrity_error")
    expect_identical(api$cp04_retained_if_release(reserve())$auxiliary, list(status = "not_extended"))
  }
})
