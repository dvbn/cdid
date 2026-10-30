# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 retained moment IF: for an eligible fitted sealed baseline
# OR/DR backend and ONE link block, the raw scale = "core" complete moment IF
# is produced once (the covariance mean pass) and sliced, in the original row
# partition, for the covariance Gram pass and both GMM projection passes.
# Every expected value comes from the unchanged forced-stream route of the
# same build (test-local instrumentation declines the admission), from the
# ordinary per-chunk replay of the same fitted backend, or from the existing
# accounting outputs; the retained route never supplies its own truth, and no
# existing gate or tolerance is weakened (public parity stays at 1e-12).
#
# STATUS: Test bodies UNEXECUTED in this adaptation job; only static Base R
# parse/AST checks ran. The root sequences the independent numerical gates.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures (test-local copies; testthat files do not share definitions)
# ---------------------------------------------------------------------------

rm_ctrl <- function(chunk = 7L, work = 1e10, matrix = 2e7) api$ipw_owner_control(max_chunk_rows = chunk,
  max_retained_entries = 1e7, max_matrix_entries = matrix, max_work_ops = work, max_output_entries = 1e7)
rm_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
rm_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))

# Balanced long panel: g = 3 treated, Inf never treated, 41 units, 5 dates.
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

# Real sealed panel/frame owners, the full-roster reader and link views.
rm_fixture <- function(d = rm_panel(), chunk = 7L, pairs = rm_pairs(), K = Inf, target = 3, ctrl = rm_ctrl(chunk)) {
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
  list(d = d, po = po, fo = fo, reader = reader, views = views, ctrl = ctrl, chunk = chunk, u = u, pairs = pairs)
}

# The unchanged general (non-complete-support) route of the same build.
rm_general <- local({
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$or_baseline_complete_facts <- function(...) NULL
  e
})
rm_or <- function(fx, ns = api, ctrl = fx$ctrl) ns$fit_or_baseline(fx$reader, fx$views, ctrl, 1e-10, outcome_recipe = rm_icpt)
rm_dr <- function(fx, ns = api, ctrl = fx$ctrl) ns$fit_dr_baseline(fx$reader, fx$views, ctrl, glm_control = rm_glm,
  rank_tol = 1e-10, outcome_recipe = rm_icpt, treatment_columns = NULL)

rm_nonadd_panel <- function(N = 41L, T = 5L, n_g = 18L, seed = 31L, id_perm = FALSE, rows_perm = FALSE) {
  set.seed(seed)
  g <- rep(0, N)
  g[sample.int(N, n_g)] <- 3
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
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

# `owner_caps` overrides fields of the resolved owner control (adjacent-cap
# cases); every other control is the ordinary public resolution.
rm_run_pipeline <- function(d, est = "reg", chunk = 7L, reduce = TRUE, weighting = "optimal",
                            xformla = NULL, weightsname = NULL, control_group = "nevertreated",
                            frame = NULL, aggregations = c("simple", "group", "dynamic", "calendar"),
                            condition_trigger = NULL, ridge_relative = NULL, owner_caps = NULL, clustered = FALSE, ns = api) {
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
    target_cohorts = tc, anticipation = 0L, control_group = control_group, base_period = "varying",
    weighting = weighting, reduce_structural_redundancy = reduce, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = aggregations, cband = FALSE, alp = 0.05, biters = 99L, seed = 101L,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = FALSE, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = "post", pt_start = NULL, cp04 = prep$cp04
  )
}

# Forced-stream reference of the same build: the admission declines, so every
# pass streams through the unchanged producer loop.
rm_stream_ns <- function() {
  e <- new.env(parent = api)
  for (nm in c("run_chained_ipw_pipeline")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$cp04_retained_if_reserve <- function(...) NULL
  e
}

# Change only the optional reservation's cap in the private adapter. Lowering
# the pipeline's original fit cap instead can correctly refuse before this
# later reservation is reached; those earlier guards are not cache admission.
rm_optional_cap_ns <- function(cap) {
  e <- rm_counting_ns()
  e$cp04_retained_if_reserve <- function(..., control) {
    control$max_retained_entries <- cap
    api$cp04_retained_if_reserve(..., control = control)
  }
  e
}

# Counting namespace: every producer chunk evaluation (the original loop and
# the retained producer pass), every `cp04_replay_links` call and the row
# partition of every evaluated chunk, with an optional forced stream.
rm_counting_ns <- function(force_stream = FALSE, decline_complete = FALSE, legacy_plan = FALSE) {
  e <- new.env(parent = api)
  for (nm in c("run_chained_ipw_pipeline", "cp04_stream_covariance", "cp04_replay_links", "cp04_retained_if_replay",
               "build_cp04_public_registry", "or_public_fit", "dr_public_fit", "fit_or_baseline", "fit_dr_baseline",
               "cp04_stream_plan")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  n <- new.env(parent = emptyenv())
  n$chunks <- 0L; n$replays <- 0L; n$rows <- list()
  e$chained_full_if_chunk <- function(backend, rows, links, scale, plan = NULL) {
    n$chunks <- n$chunks + 1L
    n$rows[[length(n$rows) + 1L]] <- rows
    api$chained_full_if_chunk(backend, rows, links, scale, plan)
  }
  replay_fn <- get("cp04_replay_links", e)
  e$cp04_replay_links <- function(backend, links, callback, scale = "core", chunk_rows = NULL) {
    n$replays <- n$replays + 1L
    replay_fn(backend, links, callback, scale = scale, chunk_rows = chunk_rows)
  }
  if (force_stream) e$cp04_retained_if_reserve <- function(...) NULL
  if (decline_complete) e$or_baseline_complete_facts <- function(...) NULL
  if (legacy_plan) e$cp04_stream_plan_full_block <- function(...) NULL
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
# Closed inventory of the declared storage and boundary work, independent of
# the new admission/cost helpers and their returned accounting records.
rm_retained_inventory <- function(backend, links, B) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  N <- reader$state$N; L <- length(links)
  raw <- 2 * N * L + 2 * N + 4 * B * L + 4 * L * L + 256
  bindings <- length(backend) + length(reader) + length(reader$state) + L +
    length(backend$complete_support$dates)
  label_bytes <- sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  metadata <- 4096 + 64 * bindings + 16 * label_bytes + 4 * B * frame$id_byte_metadata$max
  boundary_work <- 64 * metadata
  list(raw = raw, metadata = metadata, boundary_work = boundary_work,
    entries = raw + metadata, work = 128 * N * L + (40 + 8 * ceiling(N / B)) * boundary_work)
}
rm_expect_retained <- function(res, info = NULL) {
  rec <- rm_record(res)
  expect_identical(rec$status, "retained", info = info)
  expect_identical(rec$producer_passes, 1L, info = info)
  expect_identical(rec$slice_passes, 3L, info = info)
}
rm_rebind <- function(env, nm, value) {
  old <- get(nm, envir = env)
  unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
  function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
}
# Collect every batch a producer hands to its callback.
rm_collect <- function(producer, links, B, scale = "core") {
  out <- list()
  api$cp04_replay_links(producer, links, function(b) {
    force(b); out[[length(out) + 1L]] <<- b
  }, scale = scale, chunk_rows = B)
  out
}
rm_batch_meta <- function(b) b[c("rows", "id", "links", "scale", "normalization_count")]
# Public, covariance, GMM and auxiliary parity against a forced-stream result.
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
  expect_identical(a$reporting$target_meta, b$reporting$target_meta, info = info)
  expect_equal(rm_strip(a$auxiliary_reporting), rm_strip(b$auxiliary_reporting), tolerance = 1e-12, info = info)
  for (ag in names(a$aggregations)) {
    expect_equal(a$aggregations[[ag]]$overall$estimate, b$aggregations[[ag]]$overall$estimate, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$overall$se, b$aggregations[[ag]]$overall$se, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$estimates, b$aggregations[[ag]]$estimates, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$se, b$aggregations[[ag]]$se, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$if_full, b$aggregations[[ag]]$if_full, tolerance = 1e-12, info = info)
  }
  # Fits, nuisance IDs and the structural map are the same objects' values.
  expect_identical(a$nuisance_registry$nuisance_keys, b$nuisance_registry$nuisance_keys, info = info)
  expect_identical(a$nuisance_registry$cp04_aux_index, b$nuisance_registry$cp04_aux_index, info = info)
  expect_identical(a$reduction_summary$kept_indices, b$reduction_summary$kept_indices, info = info)
  expect_identical(a$resource_accounting$simultaneous_entries, b$resource_accounting$simultaneous_entries, info = info)
  expect_identical(a$resource_accounting$point_producer$chunk_rows, b$resource_accounting$point_producer$chunk_rows, info = info)
  # Original GMM batch boundaries.
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches, info = info)
  expect_identical(a$gmm_fit$if_projection$max_batch_rows_seen, b$gmm_fit$if_projection$max_batch_rows_seen, info = info)
  expect_identical(a$gmm_fit$if_projection$n_input_rows, b$gmm_fit$if_projection$n_input_rows, info = info)
}
# No retained store survives in a result (environments are inspected by kind only).
rm_find_store <- function(x, depth = 0L) {
  if (depth > 12L) return(FALSE)
  if (is.environment(x)) return(exists("kind", envir = x, inherits = FALSE) &&
    identical(get("kind", envir = x, inherits = FALSE), "cp04_retained_if"))
  if (inherits(x, "cdid_cp04_retained_if_bound")) return(TRUE)
  if (is.list(x)) for (v in x) if (rm_find_store(v, depth + 1L)) return(TRUE)
  FALSE
}

rm_d <- rm_nonadd_panel(41L, 5L, 18L, seed = 31L, id_perm = TRUE, rows_perm = TRUE)

# ---------------------------------------------------------------------------
# 1. Retained versus forced-stream parity (OR/DR, no-X and ~x, ridge-active):
#    complete IF, V, SE, public results and the active Tikhonov ridge at 1e-12;
#    one producer pass, three slice passes, original batch boundaries.
# ---------------------------------------------------------------------------
test_that("retained moment IF equals the forced stream at 1e-12 on OR and DR with and without X", {
  for (est in c("reg", "dr")) for (xf in list(NULL, ~x)) {
    info <- paste(est, if (is.null(xf)) "noX" else "~x")
    res_r <- rm_run_pipeline(rm_d, est = est, chunk = 7L, xformla = xf)
    res_s <- rm_run_pipeline(rm_d, est = est, chunk = 7L, xformla = xf, ns = rm_stream_ns())
    rec <- rm_record(res_r)
    expect_identical(rec$status, "retained", info = info)
    expect_identical(rec$scale, "core", info = info)
    expect_identical(rec$producer_passes, 1L, info = info)
    expect_identical(rec$slice_passes, 3L, info = info)
    expect_identical(rec$n_links, res_r$reduction_summary$n_links_used, info = info)
    expect_identical(rec$n_rows, res_r$nuisance_registry$cp04_reader$state$N, info = info)
    expect_identical(rec$chunk_rows, res_r$resource_accounting$point_producer$chunk_rows, info = info)
    expect_identical(rec$n_batches, as.integer(ceiling(41 / 7)), info = info)
    expect_identical(res_r$gmm_fit$if_projection$n_batches, as.integer(ceiling(41 / 7)), info = info)
    expect_identical(rm_record(res_s)$status, "streamed", info = info)
    expect_false(is.null(res_r$structural_certificate), info = info)
    expect_false(is.null(res_r$nuisance_registry$cp04$complete_support), info = info)
    rm_expect_parity(res_r, res_s, info = info)
    expect_false(rm_find_store(res_r), info = info)
  }
})

test_that("retained moment IF keeps the active Tikhonov ridge, full-catalog metric and refusals of the forced stream", {
  d_ridge <- rm_nonadd_panel(12L, 5L, 6L, seed = 55L)
  for (est in c("reg", "dr")) {
    res_r <- rm_run_pipeline(d_ridge, est = est, chunk = 7L, condition_trigger = 0.99)
    res_s <- rm_run_pipeline(d_ridge, est = est, chunk = 7L, condition_trigger = 0.99, ns = rm_stream_ns())
    expect_identical(rm_record(res_r)$status, "retained")
    expect_identical(res_r$covariance_handling$handler_used, "tikhonov")
    expect_identical(res_r$covariance_handling$gram_provenance, res_s$covariance_handling$gram_provenance)
    expect_identical(res_r$covariance_handling$dimensions$M_full, 10L)
    expect_identical(res_r$covariance_handling$dimensions$L_retained, 4L)
    expect_equal(res_r$covariance_handling$lambda, res_s$covariance_handling$lambda, tolerance = 1e-12)
    expect_equal(res_r$covariance_handling$trace_B, res_s$covariance_handling$trace_B, tolerance = 1e-12)
    expect_equal(res_r$covariance_handling$scale, res_s$covariance_handling$scale, tolerance = 1e-12)
    rm_expect_parity(res_r, res_s, info = paste("ridge", est))
    big_r <- rm_run_pipeline(d_ridge, est = est, chunk = 7L, condition_trigger = 0.99, ridge_relative = 0.5)
    big_s <- rm_run_pipeline(d_ridge, est = est, chunk = 7L, condition_trigger = 0.99, ridge_relative = 0.5, ns = rm_stream_ns())
    expect_identical(rm_record(big_r)$status, "retained")
    expect_equal(big_r$covariance_handling$lambda, big_s$covariance_handling$lambda, tolerance = 1e-12)
    expect_equal(big_r$reporting$estimates, big_s$reporting$estimates, tolerance = 1e-12)
    expect_equal(big_r$reporting$se, big_s$reporting$se, tolerance = 1e-12)
    expect_equal(big_r$gmm_fit$targets$se_criterion_conditional, big_s$gmm_fit$targets$se_criterion_conditional, tolerance = 1e-12)
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

# ---------------------------------------------------------------------------
# 2. Exactly one retained producer pass and the original batch boundaries:
#    producer chunk evaluations fall from four passes to one while the four
#    `cp04_replay_links` calls and their row partition are unchanged.
# ---------------------------------------------------------------------------
test_that("the retained route evaluates every producer chunk exactly once in the original row partition", {
  for (est in c("reg", "dr")) for (chunk in c(7L, 16L)) {
    info <- paste(est, chunk)
    n_chunks <- as.integer(ceiling(41 / chunk))
    rows_expected <- split(seq_len(41L), ceiling(seq_len(41L) / chunk))
    e_r <- rm_counting_ns()
    res_r <- rm_run_pipeline(rm_d, est = est, chunk = chunk, aggregations = "simple", ns = e_r)
    e_s <- rm_counting_ns(force_stream = TRUE)
    res_s <- rm_run_pipeline(rm_d, est = est, chunk = chunk, aggregations = "simple", ns = e_s)
    expect_identical(rm_record(res_r)$status, "retained", info = info)
    expect_identical(rm_record(res_s)$status, "streamed", info = info)
    # Mean pass, Gram pass and both GMM projection passes still enter through the replay entry point.
    expect_identical(e_r$n$replays, 4L, info = info)
    expect_identical(e_s$n$replays, 4L, info = info)
    # One producer pass instead of four; the stream keeps its four.
    expect_identical(e_r$n$chunks, n_chunks, info = info)
    expect_identical(e_s$n$chunks, 4L * n_chunks, info = info)
    expect_identical(lapply(e_r$n$rows, as.integer), unname(rows_expected), info = info)
    expect_identical(lapply(e_s$n$rows[seq_len(n_chunks)], as.integer), unname(rows_expected), info = info)
    expect_identical(rm_record(res_r)$n_batches, n_chunks, info = info)
    expect_identical(res_r$gmm_fit$if_projection$n_batches, n_chunks, info = info)
    expect_identical(res_r$gmm_fit$if_projection$max_batch_rows_seen, min(41L, chunk), info = info)
    rm_expect_parity(res_r, res_s, info = info)
  }
  # The 190-to-19 reduced T = 20 panel at chunk 50: one chunk, one producer pass, four replay calls.
  d20 <- rm_nonadd_panel(45L, 20L, 20L, seed = 77L)
  e20 <- rm_counting_ns()
  res20 <- rm_run_pipeline(d20, est = "reg", chunk = 50L, reduce = TRUE, aggregations = "simple", ns = e20)
  expect_identical(res20$reduction_summary$n_links_all, 190L)
  expect_identical(res20$reduction_summary$n_links_used, 19L)
  expect_identical(rm_record(res20)$status, "retained")
  expect_identical(rm_record(res20)$n_links, 19L)
  expect_identical(e20$n$replays, 4L)
  expect_identical(e20$n$chunks, 1L)
  res20_s <- rm_run_pipeline(d20, est = "reg", chunk = 50L, reduce = TRUE, aggregations = "simple", ns = rm_stream_ns())
  rm_expect_parity(res20, res20_s, info = "T20")
})

# ---------------------------------------------------------------------------
# 3. Store contents and refusals on the fitted backend: the producer pass hands
#    the original chunks on unchanged, slices equal them bitwise, the store is
#    the raw core-scale complete IF, and any later source/ID/order/value
#    change, boundary rebinding or callback rebinding refuses before emission.
# ---------------------------------------------------------------------------
test_that("retained slices equal the producer chunks and every mutant or rebinding is refused at the boundary", {
  fx <- rm_fixture()
  # A second sealed frame owner on the SAME panel: every existing validator
  # accepts it, so only the reserved frame snapshot can refuse it.
  fo_twin <- api$construct_ipw_frame_owner(fx$po, fx$u$id, matrix(1, nrow(fx$u), 1L, dimnames = list(NULL, "(Intercept)")),
    weights = rep(1, nrow(fx$u)), control = fx$ctrl, frame_g = fx$u$g, allow_frame_only = FALSE, design_on_A_only = TRUE)
  expect_false(identical(fo_twin$ref, fx$fo$ref))
  expect_identical(fo_twin$panel_ref, fx$fo$panel_ref)
  for (method in c("or", "dr")) {
    fit <- if (method == "or") rm_or(fx) else rm_dr(fx)
    fit2 <- if (method == "or") rm_or(fx) else rm_dr(fx)
    err <- if (method == "or") "cdid_or_baseline_integrity_error" else "cdid_dr_baseline_integrity_error"
    st <- fx$reader$state
    links <- 1:4
    plan <- api$cp04_stream_plan(fit, links, 0)
    B <- plan$chunk_rows
    expect_identical(B, 7)
    expect_identical(length(plan$blocks), 1L)
    streamed <- rm_collect(fit, links, B)
    expect_identical(length(streamed), as.integer(ceiling(41 / 7)))
    store <- api$cp04_retained_if_reserve(fit, links, B, plan, retained = 0, gmm_retained = 0, live = 0, work = 0,
      control = fit$control)
    expect_true(is.environment(store))
    expect_identical(store$state, "reserved")
    expect_false(exists("phi", envir = store, inherits = FALSE)) # reserved, not allocated
    expect_identical(store$links, 1:4)
    expect_identical(store$B, B)
    bound <- api$cp04_retained_if_bind(fit, store)
    expect_s3_class(bound, "cdid_cp04_retained_if_bound")
    # Producer pass: the original chunks, handed on unchanged.
    filled <- rm_collect(bound, links, B)
    expect_identical(store$state, "filled")
    expect_identical(store$producer_passes, 1L)
    expect_identical(store$slice_passes, 0L)
    expect_identical(length(filled), length(streamed))
    for (k in seq_along(filled)) {
      expect_identical(rm_batch_meta(filled[[k]]), rm_batch_meta(streamed[[k]]))
      expect_identical(filled[[k]]$diagnostics, streamed[[k]]$diagnostics)
      expect_equal(filled[[k]]$phi, streamed[[k]]$phi, tolerance = 1e-12)
      expect_identical(filled[[k]]$scale, "core")
      expect_identical(filled[[k]]$normalization_count, st$n_A)
    }
    # Slice pass: the stored rows, bitwise equal to the producer's, in the same partition with live IDs.
    sliced <- rm_collect(bound, links, B)
    expect_identical(store$slice_passes, 1L)
    expect_identical(store$producer_passes, 1L)
    for (k in seq_along(sliced)) {
      expect_identical(rm_batch_meta(sliced[[k]]), rm_batch_meta(filled[[k]]))
      expect_identical(sliced[[k]]$phi, filled[[k]]$phi)
      expect_identical(sliced[[k]]$id, st$id[sliced[[k]]$rows])
    }
    # The store is the raw core-scale complete IF: lambda times the frame-scale stream, uncentered.
    H_core <- do.call(rbind, lapply(streamed, `[[`, "phi"))
    H_frame <- do.call(rbind, lapply(rm_collect(fit, links, B, scale = "frame"), `[[`, "phi"))
    expect_identical(dim(store$phi), c(41L, 4L))
    expect_identical(unname(store$phi), unname(do.call(rbind, lapply(filled, `[[`, "phi"))))
    expect_equal(unname(store$phi), unname(H_core), tolerance = 1e-12)
    expect_equal(unname(store$phi), st$lambda * unname(H_frame), tolerance = 1e-12)
    expect_gt(max(abs(colMeans(store$phi))), 0)
    # Covariance through the bound producer equals the streamed covariance.
    cov_s <- api$cp04_stream_covariance(fit, links)
    store2 <- api$cp04_retained_if_reserve(fit, links, B, plan, retained = 0, gmm_retained = 0, live = 0, work = 0,
      control = fit$control)
    cov_r <- api$cp04_stream_covariance(api$cp04_retained_if_bind(fit, store2), links)
    expect_equal(cov_r, cov_s, tolerance = 1e-12)
    expect_identical(store2$state, "filled")
    expect_identical(store2$producer_passes, 1L)
    expect_identical(store2$slice_passes, 1L)

    # Mutants: refused at the first emission boundary, before any callback.
    refused <- function(producer, class = "cdid_cp04_integrity_error", links_ = links, B_ = B, scale_ = "core") {
      calls <- 0L
      expect_error(api$cp04_replay_links(producer, links_, function(b) { force(b); calls <<- calls + 1L },
        scale = scale_, chunk_rows = B_), class = class)
      expect_identical(producer$store$state, "released")
      expect_false(any(c("phi", "raw_ref", "id", "backend_ref", "reader_ref", "frame_ref", "panel_ref",
        "reference_ref", "selection_ref", "design_ref", "backend_snapshot", "reader_snapshot", "state_snapshot",
        "frame_snapshot", "panel_snapshot", "reference_snapshot", "selection_snapshot", "selection_reader_snapshot",
        "design_snapshot", "layout_snapshot", "layout_state_snapshot", "X_F_snapshot", "pair_outcomes_snapshot") %in%
        ls(producer$store, all.names = TRUE)))
      calls
    }
    orig <- store$phi
    # Every independent mutant starts with a valid, freshly filled store.
    fresh <- function() {
      store <<- api$cp04_retained_if_reserve(fit, links, B, plan, retained = 0, gmm_retained = 0,
        live = 0, work = 0, control = fit$control)
      bound <<- api$cp04_retained_if_bind(fit, store)
      expect_identical(length(rm_collect(bound, links, B)), length(filled))
      expect_identical(store$phi, orig)
    }
    store$phi[3L, 2L] <- store$phi[3L, 2L] + 1e-9                     # value
    expect_identical(refused(bound), 0L)
    fresh()
    store$phi[c(1L, 2L), ] <- orig[c(2L, 1L), ]                        # row order
    expect_identical(refused(bound), 0L)
    fresh()
    store$phi <- orig[, c(2L, 1L, 3L, 4L)]                             # column (link) order
    expect_identical(refused(bound), 0L)
    fresh()
    store$phi <- orig[seq_len(40L), , drop = FALSE]                    # shape
    expect_identical(refused(bound), 0L)
    fresh()
    # A balanced 2 x 2 perturbation is invisible to row/column sums on exact
    # dyadic values; retained identity must reject it independently of sums.
    toy <- matrix(c(1, 2, 3, 4), 2L); toy_bad <- toy + matrix(c(1, -1, -1, 1), 2L)
    expect_identical(rowSums(toy_bad), rowSums(toy)); expect_identical(colSums(toy_bad), colSums(toy))
    changed <- orig; changed[1:2, 1:2] <- changed[1:2, 1:2] + matrix(c(1, -1, -1, 1), 2L)
    store$phi <- changed
    expect_identical(refused(bound), 0L)
    fresh()
    colnames(store$phi) <- rev(colnames(orig))
    expect_identical(refused(bound), 0L)
    fresh()
    expect_identical(length(rm_collect(bound, links, B)), length(filled)) # fresh valid store serves
    ids <- store$id
    expect_true(bindingIsLocked("id", store))
    invisible(rm_rebind(store, "id", rev(ids)))                         # locked ID snapshot mutant
    expect_identical(refused(bound), 0L)
    fresh()
    lam <- store$lambda
    expect_true(bindingIsLocked("lambda", store))
    invisible(rm_rebind(store, "lambda", 2 * lam))                      # locked scale snapshot mutant
    expect_identical(refused(bound), 0L)
    # Requests outside the reserved block, scale or partition are refused.
    fresh()
    expect_identical(refused(bound, links_ = c(2L, 1L, 3L, 4L)), 0L)
    fresh()
    expect_identical(refused(bound, links_ = 1:3), 0L)
    fresh()
    expect_identical(refused(bound, B_ = B - 1), 0L)
    fresh()
    expect_identical(refused(bound, scale_ = "frame"), 0L)
    # Source: another fit of the same fixture, or a resealed copy, cannot bind.
    fresh()
    expect_error(api$cp04_retained_if_bind(fit2, store), class = "cdid_cp04_integrity_error")
    resealed <- fit; resealed$ref <- NULL
    resealed <- api$ipw_owner_seal(resealed, if (method == "or") "or_baseline_backend" else "dr_baseline_backend")
    fresh()
    expect_error(api$cp04_retained_if_bind(resealed, store), class = "cdid_cp04_integrity_error")
    # Boundary rebinding: the panel owner (as in the complete-support tests) refuses through the
    # original validator; another sealed frame owner refuses through the frame snapshot.
    po <- fx$po
    shifted <- po; shifted$postings$date_start[2L] <- shifted$postings$date_start[2L] + 1L
    fresh()
    restore <- rm_rebind(po$ref, "panel_owner", shifted)
    expect_identical(refused(bound, class = err), 0L)
    restore()
    # Independently valid replacement objects behind the old binding retain
    # the roster/layout; their own seals cannot establish the reserved source.
    d_twin <- fx$d; d_twin$y <- d_twin$y + d_twin$t
    po_twin <- api$construct_ipw_panel_owner(d_twin, "id", "t", "y", "g", fx$ctrl)
    fresh()
    restore <- rm_rebind(po$ref, "panel_owner", po_twin)
    expect_silent(api$cp04_retained_if_validate_backend(fit))
    expect_identical(refused(bound), 0L)
    restore()
    reader_twin <- api$construct_cp04_frame_reader(fx$fo, fx$u$id,
      matrix(1, nrow(fx$u), 1L, dimnames = list(NULL, "reader_intercept")), fx$u$g,
      rep(1, nrow(fx$u)), fx$ctrl, allow_unknown_off_A = TRUE)
    fresh()
    restore <- rm_rebind(fit$reader_ref, "cp04_reader", reader_twin)
    expect_silent(api$cp04_retained_if_validate_backend(fit))
    expect_identical(refused(bound), 0L)
    restore()
    fresh()
    restore <- rm_rebind(fx$fo$ref, "frame_owner", fo_twin)
    expect_identical(refused(bound), 0L)
    restore()
    fresh()
    expect_identical(length(rm_collect(bound, links, B)), length(filled))
    # Callback rebinding: the next slice is refused before it is emitted.
    restore <- NULL; calls <- 0L
    expect_error(api$cp04_replay_links(bound, links, function(b) {
      force(b); calls <<- calls + 1L
      if (is.null(restore)) restore <<- rm_rebind(po$ref, "panel_owner", shifted)
    }, chunk_rows = B), class = err)
    restore()
    expect_identical(calls, 1L)
    expect_identical(store$state, "released")
    expect_false(exists("phi", envir = store, inherits = FALSE))
    expect_false(exists("backend_ref", envir = store, inherits = FALSE))
    # Explicit successful release is checked on another fresh valid store.
    fresh()
    # Released: nothing remains and the store cannot serve.
    rec <- api$cp04_retained_if_release(store)
    expect_identical(rec$status, "retained")
    expect_identical(rec$producer_passes, 1L)
    expect_identical(rec$n_batches, as.integer(ceiling(41 / 7)))
    inventory <- rm_retained_inventory(fit, links, B)
    expect_identical(rec$reserved_entries, inventory$entries)
    expect_identical(rec$snapshot_metadata_entries, inventory$metadata)
    expect_identical(rec$boundary_check_work_ops, inventory$boundary_work)
    expect_identical(rec$work_ops_bound, inventory$work)
    expect_false(exists("phi", envir = store, inherits = FALSE))
    expect_false(exists("id", envir = store, inherits = FALSE))
    expect_identical(store$state, "released")
    expect_identical(refused(bound), 0L)
    expect_identical(api$cp04_retained_if_release(NULL), list(status = "streamed"))
  }
})

# ---------------------------------------------------------------------------
# 4. Admission: exact entry/work arithmetic, the adjacent caps, partition and
#    block identity, planner probes and ineligible backends, all before any
#    allocation.
# ---------------------------------------------------------------------------
test_that("retained moment IF admission is exact at the adjacent caps and declines before allocation", {
  expect_identical(api$cp04_retained_if_entries(10000, 19, 171), 414696)
  expect_identical(api$cp04_retained_if_work(10000, 19), 24320000)
  fx <- rm_fixture()
  for (method in c("or", "dr")) {
    fit <- if (method == "or") rm_or(fx) else rm_dr(fx)
    fit_g <- if (method == "or") rm_or(fx, rm_general) else rm_dr(fx, rm_general)
    links <- 1:4; L <- 4L; N <- 41L
    plan <- api$cp04_stream_plan(fit, links, 0)
    B <- plan$chunk_rows
    ctrl <- fit$control
    reserve <- function(..., backend = fit, links_ = links, B_ = B, plan_ = plan, retained = 0, gmm_retained = 0,
                        live = 0, work = 0, control = ctrl)
      api$cp04_retained_if_reserve(backend, links_, B_, plan_, retained = retained, gmm_retained = gmm_retained,
        live = live, work = work, control = control)
    inventory <- rm_retained_inventory(fit, links, B)
    raw_entries <- api$cp04_retained_if_entries(N, L, B)
    expect_identical(raw_entries, 2 * N * L + 2 * N + 4 * B * L + 4 * L * L + 256)
    expect_identical(raw_entries, inventory$raw)
    entries <- inventory$entries
    cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
    # Entry cap: admitted exactly at the cap, declined one above.
    at <- reserve(live = cap - entries)
    expect_true(is.environment(at))
    expect_identical(at$admission_total, cap)
    expect_identical(at$cap, cap)
    expect_identical(at$entries, entries)
    expect_identical(at$metadata_entries, inventory$metadata)
    expect_identical(at$boundary_work, inventory$boundary_work)
    expect_identical(at$work, inventory$work)
    expect_false(exists("phi", envir = at, inherits = FALSE))
    expect_null(reserve(live = cap - entries + 1))
    # Work cap: admitted exactly at the cap, declined one above.
    cache_work <- inventory$work
    expect_true(is.environment(reserve(work = ctrl$max_work_ops - cache_work)))
    expect_null(reserve(work = ctrl$max_work_ops - cache_work + 1))
    # Partition identity: a different GMM partition or a two-block covariance plan declines.
    expect_null(reserve(B_ = B - 1))
    plan10 <- api$cp04_stream_plan_legacy(fit, 1:10, 0)
    expect_identical(length(plan10$blocks), 2L)
    expect_null(reserve(links_ = 1:10, plan_ = plan10))
    # Planner probe: the GMM preflight with the cache charged refuses one entry past its bound and
    # that typed bound refusal declines the store; at the bound the store is admitted.
    pre <- api$chained_full_if_budget(fit, B, links)
    own <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
    room <- ctrl$max_matrix_entries - own - pre$output - pre$scratch - entries
    expect_true(is.environment(reserve(gmm_retained = room)))
    expect_null(reserve(gmm_retained = room + 1))
    # General sealed backends are eligible; absence of a backend still declines.
    expect_null(fit_g$complete_support)
    expect_true(is.environment(reserve(backend = fit_g, plan_ = api$cp04_stream_plan(fit_g, links, 0))))
    expect_null(reserve(backend = NULL))
    # Invalid requests preserve input/control failures, before any allocation.
    expect_error(reserve(links_ = c(1L, 1L, 2L, 3L)), class = "cdid_cp04_input_error")
    expect_error(reserve(live = -1), class = "cdid_ipw_owner_control_error")
    expect_error(reserve(live = 0.5), class = "cdid_ipw_owner_control_error")
    # Integrity errors of the existing validators are not swallowed.
    broken <- fit; broken$links <- broken$links[1:3]
    err_class <- if (method == "or") "cdid_or_baseline_integrity_error" else "cdid_dr_baseline_integrity_error"
    expect_error(reserve(backend = broken), class = err_class)
    expect_error(api$cp04_retained_if_bind(fit, new.env()), class = "cdid_cp04_integrity_error")
  }
})

test_that("pipeline admission counts the complete envelope and flips to the ordinary stream one entry below the cap", {
  for (est in c("reg", "dr")) {
    res <- rm_run_pipeline(rm_d, est = est, chunk = 7L, aggregations = "simple")
    rec <- rm_record(res)
    expect_identical(rec$status, "retained")
    oc <- res$controls$pipeline_control$owner_control
    expect_identical(rec$cap_entries, min(oc$max_retained_entries, oc$max_matrix_entries))
    N <- res$nuisance_registry$cp04_reader$state$N
    L <- res$reduction_summary$n_links_used
    B <- res$resource_accounting$point_producer$chunk_rows
    aux_env <- 8 * N * (length(res$nuisance_registry$cp04_aux$selections[["varying"]]) + nrow(res$support$ledger))
    backend <- res$nuisance_registry$cp04
    links <- vapply(res$point_links_used, `[[`, integer(1), "cp04_index")
    inventory <- rm_retained_inventory(backend, links, B)
    expect_identical(rec$reserved_entries, inventory$entries)
    expect_identical(rec$snapshot_metadata_entries, inventory$metadata)
    expect_identical(rec$boundary_check_work_ops, inventory$boundary_work)
    # This reduced catalog uses the legacy single block, whose peak is the
    # original producer budget. Include any excess over the GMM reservation.
    expect_lte(L, backend$control$max_chunk_rows)
    pre <- api$chained_full_if_budget(backend, B, links)
    covariance_peak <- backend$accounting$owned_entries + backend$accounting$borrowed_entries + pre$output + pre$scratch
    covariance_excess <- max(0, covariance_peak - res$resource_accounting$point_producer$producer_workspace_entries)
    total <- res$resource_accounting$simultaneous_entries +
      res$resource_accounting$final_diagnostic_copy_entries + covariance_excess + aux_env + inventory$entries
    expect_identical(rec$admission_total_entries, total)
    expect_identical(rec$work_ops_bound, inventory$work)
    expect_lte(rec$admission_total_entries, rec$cap_entries)
    res_at <- rm_run_pipeline(rm_d, est = est, chunk = 7L, aggregations = "simple",
      ns = rm_optional_cap_ns(total))
    expect_identical(rm_record(res_at)$status, "retained")
    expect_identical(rm_record(res_at)$cap_entries, total)
    e_below <- rm_optional_cap_ns(total - 1)
    res_below <- rm_run_pipeline(rm_d, est = est, chunk = 7L, aggregations = "simple", ns = e_below)
    expect_identical(rm_record(res_below)$status, "streamed")
    expect_identical(e_below$n$replays, 4L)
    expect_identical(e_below$n$chunks, 4L * as.integer(ceiling(N / B)))
    expect_identical(unname(e_below$n$rows), rep(unname(split(seq_len(N), ceiling(seq_len(N) / B))), 4L))
    expect_identical(res_below$resource_accounting$simultaneous_entries, res$resource_accounting$simultaneous_entries)
    rm_expect_parity(res, res_at, info = paste(est, "at"))
    rm_expect_parity(res, res_below, info = paste(est, "below"))
  }
})

# ---------------------------------------------------------------------------
# 5. Multi-block capacity declines keep the original stream. Supported
#    missing/weighted/frame/cluster routes retain when admitted; IPW has no
#    CP04 producer.
# ---------------------------------------------------------------------------
test_that("capacity declines keep stream counts and newly eligible paths retain with full parity", {
  n_chunks <- as.integer(ceiling(41 / 7))
  # Force the legacy branch to retain the original two-block fallback coverage.
  for (est in c("reg", "dr")) {
    e <- rm_counting_ns(legacy_plan = TRUE)
    res <- rm_run_pipeline(rm_d, est = est, chunk = 7L, reduce = FALSE, aggregations = "simple", ns = e)
    expect_identical(res$reduction_summary$n_links_used, 10L)
    expect_identical(rm_record(res)$status, "streamed")
    expect_identical(e$n$replays, 2L + 3L + 2L)
    expect_identical(e$n$chunks, 7L * n_chunks)
    # The general route also declines this forced legacy two-block plan.
    e_g <- rm_counting_ns(decline_complete = TRUE, legacy_plan = TRUE)
    res_g <- rm_run_pipeline(rm_d, est = est, chunk = 7L, aggregations = "simple", ns = e_g)
    expect_null(res_g$nuisance_registry$cp04$complete_support)
    expect_identical(rm_record(res_g)$status, "streamed")
    # Declining complete-support facts also declines structural reduction:
    # ten links give two mean passes, three Gram pairs and two GMM passes.
    expect_identical(res_g$reduction_summary$n_links_used, 10L)
    expect_identical(e_g$n$replays, 7L)
    expect_identical(e_g$n$chunks, 7L * n_chunks)
  }
  d0 <- rm_nonadd_panel(41L, 5L, 18L, seed = 31L)
  dm <- d0[!(d0$id == 5L & d0$period == 3L), ]
  dw <- d0; dw$w <- 1 + (dw$id %% 3) / 4
  fr <- d0[!duplicated(d0$id), c("id", "G", "x", "cl")]
  fr <- rbind(fr[order(fr$id), ], data.frame(id = 9999L, G = 0, x = 0, cl = 99L))
  for (case in list(list(nm = "varying_weights", run = function(ns) rm_run_pipeline(dw, est = "dr", weightsname = "w", ns = ns)),
                    list(nm = "missing_date", run = function(ns) rm_run_pipeline(dm, est = "reg", ns = ns)),
                    list(nm = "frame_only_unit", run = function(ns) rm_run_pipeline(d0, est = "reg", frame = fr, ns = ns)))) {
    e <- rm_counting_ns()
    res <- case$run(e)
    ref <- case$run(rm_stream_ns())
    expect_null(res$nuisance_registry$cp04$complete_support, info = case$nm)
    rm_expect_retained(res, info = case$nm)
    B <- res$resource_accounting$point_producer$chunk_rows
    N <- res$nuisance_registry$cp04_reader$state$N
    expect_identical(e$n$replays, 4L, info = case$nm)
    expect_identical(e$n$chunks, as.integer(ceiling(N / B)), info = case$nm)
    expect_identical(unname(e$n$rows), unname(split(seq_len(N), ceiling(seq_len(N) / B))), info = case$nm)
    rm_expect_parity(res, ref, info = case$nm)
    expect_false(rm_find_store(res), info = case$nm)
  }
  # Clustered inference and constant non-unit weights now retain the same raw
  # moment rows, before the unchanged inference consumers.
  dunit <- d0; dunit$w <- 2
  for (est in c("reg", "dr")) for (kind in c("clustered", "constant_weight")) {
    dat <- if (kind == "clustered") d0 else dunit
    args <- list(d = dat, est = est, aggregations = "simple", clustered = kind == "clustered",
      weightsname = if (kind == "constant_weight") "w" else NULL)
    e <- rm_counting_ns()
    res <- do.call(rm_run_pipeline, c(args, list(ns = e)))
    ref <- do.call(rm_run_pipeline, c(args, list(ns = rm_stream_ns())))
    expect_false(is.null(res$nuisance_registry$cp04$complete_support), info = kind)
    rm_expect_retained(res, info = kind)
    expect_identical(e$n$replays, 4L, info = kind)
    expect_identical(e$n$chunks, n_chunks, info = kind)
    expect_identical(unname(e$n$rows), unname(split(seq_len(41L), ceiling(seq_len(41L) / 7L))), info = kind)
    rm_expect_parity(res, ref, info = paste(est, kind))
    expect_false(rm_find_store(res), info = kind)
  }
  # IPW: no CP04 producer, no record, nothing changed.
  res_ipw <- rm_run_pipeline(d0, est = "ipw", aggregations = "simple")
  expect_null(res_ipw$nuisance_registry$cp04)
  expect_null(rm_record(res_ipw))
  expect_false(rm_find_store(res_ipw))
})

# ---------------------------------------------------------------------------
# 6. Both GMM exact replay checks stay active behind the bound producer: the
#    solve over sliced batches equals the pipeline's, and a changed second
#    pass or a changed ID sequence is refused by the core itself.
# ---------------------------------------------------------------------------
test_that("the GMM core's exact replay, ID, weight and batch checks run unchanged on retained slices", {
  res <- rm_run_pipeline(rm_d, est = "reg", chunk = 7L, aggregations = "simple")
  b <- res$nuisance_registry$cp04
  st <- res$nuisance_registry$cp04_reader$state
  cp04_used <- vapply(res$point_links_used, `[[`, integer(1), "cp04_index")
  link_ids <- res$design_used$incidence$link_ids
  B <- res$resource_accounting$point_producer$chunk_rows
  plan <- api$cp04_stream_plan(b, cp04_used, 0)
  nuis_id <- "retained_if_test_nuis"
  prov <- api$gmm_covariance_provenance(
    row_level = "unit", n_units = st$N, n_clusters = st$N, n_positive_contributing_rows = st$N,
    centered = TRUE, weighted = FALSE, weight_type = "unweighted", normalization_factor = 1 / st$n_A,
    scale_id = "retained_if_test_scale", nuisance_fit_id = nuis_id,
    reduction = list(mode = "none", reduction_active = FALSE))
  prov$complete_unit_if_gram <- api$ipw_seal(list(covariance = res$covariance_stream$V, n_units = st$N,
    nuisance_fit_id = nuis_id), "complete_unit_if_gram")
  solve_with <- function(callback, expected_ids = st$id) api$solve_chained_gmm_core(
    design = res$design_used, mhat = res$covariance_stream$mhat, V = res$covariance_stream$V, provenance = prov,
    targets = api$ipw_pipeline_core_target_plan(res$support$ledger, res$design_used)$target_specs,
    weighting = "optimal",
    if_replay = api$gmm_if_replay_spec(callback = callback, input_row_level = "unit", aggregate_to_clusters = FALSE,
      n_input_rows = st$N, n_output_rows = st$N, n_positive_contributing_rows = st$N, centered = TRUE,
      weighted = FALSE, weight_type = "unweighted", normalization_factor = 1 / st$n_A,
      scale_id = prov$scale_id, nuisance_fit_id = nuis_id, expected_unit_ids = expected_ids),
    accuracy = res$controls$pipeline_control$gmm_accuracy, bounds = res$controls$pipeline_control$gmm_bounds,
    gauge_rule = "lowest", request_param_cov = TRUE, materialize_target_if = TRUE)
  emit <- function(consume, producer, tamper = NULL) {
    api$cp04_replay_links(producer, cp04_used, function(bt) {
      phi <- bt$phi; colnames(phi) <- link_ids
      if (!is.null(tamper)) phi <- tamper(phi)
      consume(list(unit_id = bt$id, phi = phi, weight = rep(1, nrow(phi)), nuisance_fit_id = nuis_id))
    }, chunk_rows = B)
  }
  reserve <- function() api$cp04_retained_if_reserve(b, cp04_used, B, plan, retained = 0, gmm_retained = 0, live = 0,
    work = 0, control = b$control)
  # Streamed and retained solves agree with each other and with the pipeline at 1e-12.
  g_stream <- solve_with(function(consume) emit(consume, b))
  store <- reserve(); bound <- api$cp04_retained_if_bind(b, store)
  g_ret <- solve_with(function(consume) emit(consume, bound))
  expect_identical(store$producer_passes, 1L)
  expect_identical(store$slice_passes, 1L)
  expect_equal(g_ret$targets, g_stream$targets, tolerance = 1e-12)
  expect_equal(g_ret$if_projection$sandwich_cov, g_stream$if_projection$sandwich_cov, tolerance = 1e-12)
  expect_equal(g_ret$if_projection$materialized_target_if, g_stream$if_projection$materialized_target_if, tolerance = 1e-12)
  expect_equal(g_ret$targets$estimate, res$gmm_fit$targets$estimate, tolerance = 1e-12)
  expect_identical(g_ret$if_projection$n_batches, g_stream$if_projection$n_batches)
  expect_identical(g_ret$if_projection$n_batches, as.integer(ceiling(st$N / B)))
  # A second pass whose projected rows differ is refused by the core's exact replay check.
  store2 <- reserve(); bound2 <- api$cp04_retained_if_bind(b, store2)
  pass <- 0L
  expect_error(solve_with(function(consume) {
    pass <<- pass + 1L
    emit(consume, bound2, tamper = if (pass == 2L) function(phi) { phi[1L, 1L] <- phi[1L, 1L] + 1; phi } else NULL)
  }), class = "cdid_gmm_replay_error")
  expect_identical(pass, 2L)
  # A changed expected ID sequence is refused by the core's ID check on sliced batches.
  store3 <- reserve(); bound3 <- api$cp04_retained_if_bind(b, store3)
  expect_error(solve_with(function(consume) emit(consume, bound3), expected_ids = rev(st$id)), class = "cdid_gmm_replay_error")
  # A changed weight is refused on sliced batches.
  store4 <- reserve(); bound4 <- api$cp04_retained_if_bind(b, store4)
  expect_error(solve_with(function(consume) api$cp04_replay_links(bound4, cp04_used, function(bt) {
    phi <- bt$phi; colnames(phi) <- link_ids
    consume(list(unit_id = bt$id, phi = phi, weight = c(2, rep(1, nrow(phi) - 1L)), nuisance_fit_id = nuis_id))
  }, chunk_rows = B)), class = "cdid_gmm_replay_error")
})
