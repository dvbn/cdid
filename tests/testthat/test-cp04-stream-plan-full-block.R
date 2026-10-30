# SPDX-License-Identifier: GPL-3.0-only
#
# Candidate02 / D90 covariance planning (`cp04_stream_plan`): keep the
# original legacy plan and row partition, and prospectively admit ONE
# full-width covariance block at the exact same declared row partition for
# fitted baseline OR/DR (complete-support or general/missing) under the existing
# producer/output budgets and complete simultaneous Gram/scratch/backend/
# retained charges (including full-width consumer scratch, owner metadata
# and boundary validation work) and
# full_work <= legacy_work.
#
# STATUS: Test bodies UNEXECUTED in the bounded review job. Base R parse and
# signature checks only; root sequences the installed numerical gates.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures and namespace clones (from supplied CP08 and CP09 suites)
# ---------------------------------------------------------------------------

rm_ctrl <- function(chunk = 7L, work = 1e10, matrix = 2e7, retained = 1e7, output = 1e7)
  api$ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = retained,
    max_matrix_entries = matrix, max_work_ops = work, max_output_entries = output)
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

rm_general <- local({
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$or_baseline_complete_facts <- function(...) NULL
  e
})
rm_or <- function(fx, ns = api, ctrl = fx$ctrl) ns$fit_or_baseline(fx$reader, fx$views, ctrl, 1e-10, outcome_recipe = rm_icpt)
rm_dr <- function(fx, ns = api, ctrl = fx$ctrl, external_work_ops = 0, automatic_chunk_rows = FALSE)
  ns$fit_dr_baseline(fx$reader, fx$views, ctrl, glm_control = rm_glm, rank_tol = 1e-10,
    outcome_recipe = rm_icpt, treatment_columns = NULL, external_work_ops = external_work_ops,
    automatic_chunk_rows = automatic_chunk_rows)

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

rm_run_pipeline <- function(d, est = "dr", chunk = 7L, reduce = TRUE, weighting = "optimal",
                            xformla = NULL, weightsname = NULL, control_group = "nevertreated",
                            frame = NULL, aggregations = c("simple", "group", "dynamic", "calendar"),
                            condition_trigger = NULL, ridge_relative = NULL, owner_caps = NULL,
                            clustered = FALSE, ns = api) {
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

# Test-local namespace clone supporting forced legacy multi-block planner,
# forced stream (no retained IF store), general nuisance route, optional cap,
# and producer/replay counting.
plan_counting_ns <- function(legacy_plan = FALSE, force_stream = FALSE, decline_complete = FALSE, cap = NULL) {
  e <- new.env(parent = api)
  for (nm in c("att_gt", "chained_public_att_gt", "run_chained_ipw_pipeline",
               "stream_ipw_point_covariance_view", "cp04_stream_covariance",
               "cp04_stream_plan", "cp04_retained_if_reserve", "cp04_replay_links",
               "cp04_retained_if_replay", "build_cp04_public_registry",
               "or_public_fit", "dr_public_fit", "fit_or_baseline", "fit_dr_baseline")) {
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
  if (legacy_plan) e$cp04_stream_plan_full_block <- function(...) NULL
  if (force_stream) e$cp04_retained_if_reserve <- function(...) NULL
  if (decline_complete) e$or_baseline_complete_facts <- function(...) NULL
  if (!is.null(cap)) {
    reserve_fn <- get("cp04_retained_if_reserve", e)
    e$cp04_retained_if_reserve <- function(..., control) {
      control$max_retained_entries <- cap
      reserve_fn(..., control = control)
    }
  }
  e$n <- n
  e
}

rm_record <- function(res) res$resource_accounting$point_producer$retained_moment_if
# Closed storage/owner/header/label/ID-slice inventory and boundary work;
# expected costs never call the new production admission/cost helpers.
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
rm_expect_retained <- function(res) {
  rec <- rm_record(res)
  expect_identical(rec$status, "retained")
  expect_identical(rec$producer_passes, 1L)
  expect_identical(rec$slice_passes, 3L)
}
rm_find_store <- function(x, depth = 0L) {
  if (depth > 12L) return(FALSE)
  if (is.environment(x)) return(exists("kind", envir = x, inherits = FALSE) &&
    identical(get("kind", envir = x, inherits = FALSE), "cp04_retained_if"))
  if (inherits(x, "cdid_cp04_retained_if_bound")) return(TRUE)
  if (is.list(x)) for (v in x) if (rm_find_store(v, depth + 1L)) return(TRUE)
  FALSE
}

cp09_panel <- function() {
  N <- 53L; tt <- c(2, 4, 9, 15, 26, 40)
  u <- seq_len(N)
  gg <- rep(c(9, 26, 0), c(13L, 17L, 23L))[order(sin(u * 0.83))]
  d <- data.frame(id = rep(1000L + rev(u), each = length(tt)), period = rep(tt, N))
  i <- rep(u, each = length(tt)); k <- rep(seq_along(tt), N)
  d$G <- gg[i]
  d$x <- sin(i * 0.71) + i / 100
  d$x2 <- cos(i * 0.37) + sin(i * 0.19)
  d$x3 <- sin(i * 0.27) + cos(i * 0.43)
  d$y <- 0.2 * i + 0.1 * k + 0.12 * d$x * k^2 + 0.09 * d$x2 * k +
    sin(i * k * 0.31) + cos(i * k * 0.17) +
    ifelse(d$G > 0 & d$period >= d$G, 1 + k / 7, 0)
  d$cl <- (i - 1L) %/% 3L + 1L
  d
}
cp09_ctrl <- function(chunk = 9L, matrix = 2e7, work = 1e13, retained = 2e7, output = 2e7)
  api$ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = retained,
    max_matrix_entries = matrix, max_work_ops = work, max_output_entries = output)
cp09_recipe <- function(columns = c("(Intercept)", "x", "x2"), intercept = TRUE)
  api$or_pooled_outcome_recipe(columns, intercept, setdiff(columns, "(Intercept)"))

cp09_fixture <- function(cg = "notyettreated", est = "dr", chunk = 9L,
                         perm = FALSE, wt = 1, ctrl = cp09_ctrl(chunk)) {
  d <- cp09_panel(); d$G[d$G == 0] <- Inf
  if (perm) d <- d[order(cos(seq_len(nrow(d)) * 0.7)), ]
  po <- api$construct_ipw_panel_owner(d, "id", "period", "y", "G", ctrl)
  u <- d[!duplicated(d$id), ]; u <- u[order(sin(u$id * 0.3)), ]
  X <- cbind(`(Intercept)` = 1, x = u$x, x2 = u$x2)
  recipe <- cp09_recipe(colnames(X))
  Z <- if (est == "dr") cbind(X, x3 = u$x3) else X
  ww <- if (length(wt) == 1L) rep(wt, nrow(u)) else wt[match(u$id, sort(u$id))]
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = ww, control = ctrl,
    frame_g = u$G, allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, u$id,
    matrix(1, nrow(u), 1L, dimnames = list(NULL, "reader_intercept")), u$G, ww, ctrl,
    allow_unknown_off_A = TRUE)
  tt <- sort(unique(d$period)); views <- lks <- list()
  for (g in c(9, 26)) for (s in seq_len(length(tt) - 1L)) for (t in seq.int(s + 1L, length(tt))) {
    K <- if (cg == "nevertreated") Inf else c(9, 26, Inf)[c(9, 26, Inf) > max(g, tt[t])]
    j <- length(views) + 1L; id <- paste0("L", j)
    views[[j]] <- api$construct_ipw_link_view(fo, tt[s], tt[t], g = g, K = K, control = ctrl,
      link_id = id, functional_key = paste0("cmpcat:cp09:", j))
    lks[[j]] <- list(s = tt[s], t = tt[t], g = g, K = K, link_id = id)
  }
  if (perm) {
    o <- c(seq.int(length(views), 1L, by = -2L), seq.int(length(views) - 1L, 1L, by = -2L))
    views <- views[o]; lks <- lks[o]
  }
  list(d = d, u = u, po = po, fo = fo, reader = reader, views = views, links = lks,
    X = X, V = if (est == "dr") Z[, c("(Intercept)", "x", "x3"), drop = FALSE] else X,
    calendar = tt, weights = ww, recipe = recipe, ctrl = ctrl, est = est)
}
cp09_fit <- function(fx, ns = api, external_work_ops = 0, automatic_chunk_rows = FALSE) {
  if (fx$est == "or") return(ns$fit_or_baseline(fx$reader, fx$views, fx$ctrl,
    rank_tol = 1e-10, outcome_recipe = fx$recipe, external_work_ops = external_work_ops))
  ns$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = rm_glm, rank_tol = 1e-10,
    outcome_recipe = fx$recipe, treatment_columns = colnames(fx$V),
    external_work_ops = external_work_ops, automatic_chunk_rows = automatic_chunk_rows)
}
cp09_chunk <- function(fit, rows, links, details = FALSE, scale = "frame", plan = NULL) {
  fn <- if (inherits(fit, "cdid_dr_baseline")) api$evaluate_dr_baseline_if_chunk else
    api$evaluate_or_baseline_if_chunk
  fn(fit, rows, links, details = details, scale = scale, plan = plan)
}
cp09_phi <- function(fit, links = seq_along(fit$links), chunk = fit$control$max_chunk_rows, details = FALSE) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  H <- matrix(0, st$N, length(links))
  for (first in seq.int(1L, st$N, by = chunk)) {
    rr <- seq.int(first, min(st$N, first + chunk - 1L))
    H[rr, ] <- cp09_chunk(fit, rr, links, details)$phi
  }
  H
}
cp09_run <- function(d = cp09_panel(), est = "dr", cg = "notyettreated", ns = api,
                     chunk = 9L, xf = ~x + x2, weightsname = NULL, condition_trigger = NULL,
                     ridge_relative = NULL, base_period = "varying", clustered = FALSE) {
  cc <- list(chunk_rows = chunk)
  if (!is.null(condition_trigger)) cc$condition_trigger <- condition_trigger
  if (!is.null(ridge_relative)) cc$ridge_relative <- ridge_relative
  r <- api$chained_public_resolve_control(cc, gmm_weighting = "optimal", gmm_covariance = "tikhonov")
  p <- api$chained_request_prepare_inputs(data = d, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xf, weightsname = weightsname, clustervars = if (clustered) "cl" else NULL, anticipation = 0, frame = NULL,
    outcome_regression = identical(est, "reg"), doubly_robust = identical(est, "dr"), or_resource = r$profile)
  p$cp04$automatic_chunk_rows <- FALSE
  dims <- api$chained_public_dimensions(p, biters = 99L, clustered = clustered)
  b <- api$chained_public_pipeline_control(dims, r, gmm_weighting = "optimal", gmm_covariance = "tikhonov")
  ns$run_chained_ipw_pipeline(data = p$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = b$control, xformla = p$xformla, weightsname = NULL, cluster_vec = if (clustered) p$cluster_full else NULL,
    target_cohorts = p$target_cohorts, anticipation = 0L, control_group = cg, base_period = base_period,
    weighting = "optimal", reduce_structural_redundancy = TRUE, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = c("simple", "group", "dynamic", "calendar"), cband = FALSE, alp = 0.05, biters = 99L,
    seed = 101L, frame_id = p$frame_ids, frame_Z = p$frame_Z, frame_weights = p$frame_w, frame_g = p$frame_g_arg,
    allow_frame_only = p$n_frame_only > 0L, na.rm = FALSE, bstrap = FALSE, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = "post", pt_start = NULL, cp04 = p$cp04)
}
cp09_equal <- function(a, b) {
  expect_true(is.matrix(a$reporting$if_full)); expect_true(is.matrix(a$covariance_stream$V))
  for (nm in c("mhat", "V", "mean_phi"))
    expect_equal(a$covariance_stream[[nm]], b$covariance_stream[[nm]], tolerance = 1e-12)
  expect_equal(a$gmm_fit$targets, b$gmm_fit$targets, tolerance = 1e-12)
  for (nm in c("sandwich_cov", "materialized_target_if"))
    expect_equal(a$gmm_fit$if_projection[[nm]], b$gmm_fit$if_projection[[nm]], tolerance = 1e-12)
  ah <- a$covariance_handling; bh <- b$covariance_handling
  costs <- c("producer_workspace_entries", "peak_concurrent_entries", "work_ops_bound")
  for (z in list(a, b)) {
    conditioning <- z$covariance_handling$conditioning
    count <- conditioning$nonpositive_eigenvalue_count
    expect_true(is.numeric(count) && length(count) == 1L && is.finite(count) && count == floor(count))
    expect_gte(count, 0L); expect_lte(count, length(conditioning$eigenvalues_B))
    expect_identical(count, as.integer(sum(conditioning$eigenvalues_B <= 0)))
    pre <- z$auxiliary_reporting$preflight
    ctrl <- z$controls$pipeline_control$owner_control
    expect_true(all(costs %in% names(pre)))
    for (nm in costs) {
      value <- pre[[nm]]
      expect_true(is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= 0 && value == floor(value) && value <= 2^53 - 1)
      expect_lte(value, if (nm == "work_ops_bound") ctrl$max_work_ops else ctrl$max_matrix_entries)
    }
    expect_lte(pre$producer_workspace_entries, pre$peak_concurrent_entries)
  }
  ah$conditioning$nonpositive_eigenvalue_count <- NULL
  bh$conditioning$nonpositive_eigenvalue_count <- NULL
  expect_equal(ah, bh, tolerance = 1e-12)
  expect_equal(a$full_frame_transport$Phi_full, b$full_frame_transport$Phi_full, tolerance = 1e-12)
  for (nm in c("estimates", "se", "analytical_se", "if_full", "vcov"))
    expect_equal(a$reporting[[nm]], b$reporting[[nm]], tolerance = 1e-12)
  expect_identical(a$reporting$target_meta, b$reporting$target_meta)
  ax <- a$auxiliary_reporting; bx <- b$auxiliary_reporting
  ax$preflight[costs] <- NULL; bx$preflight[costs] <- NULL
  expect_equal(ax, bx, tolerance = 1e-12)
  expect_identical(a$nuisance_registry$cp04_aux_index, b$nuisance_registry$cp04_aux_index)
  expect_identical(a$reduction_summary$kept_indices, b$reduction_summary$kept_indices)
  expect_identical(a$resource_accounting$point_producer$chunk_rows, b$resource_accounting$point_producer$chunk_rows)
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches)
  expect_identical(a$gmm_fit$if_projection$max_batch_rows_seen, b$gmm_fit$if_projection$max_batch_rows_seen)
  # The public pretest consumes these same complete reported cells/covariance.
  pa <- api$chained_pretest_record(a$reporting)
  pb <- api$chained_pretest_record(b$reporting)
  expect_identical(pa$status, pb$status)
  expect_equal(pa, pb, tolerance = 1e-12)
  expect_identical(sort(names(a$aggregations)), sort(names(b$aggregations)))
  for (ag in names(a$aggregations)) for (nm in c("estimates", "se", "if_full", "overall"))
    expect_equal(a$aggregations[[ag]][[nm]], b$aggregations[[ag]][[nm]], tolerance = 1e-12)
}

rm_d <- rm_nonadd_panel(41L, 5L, 18L, seed = 31L, id_perm = TRUE, rows_perm = TRUE)

# ---------------------------------------------------------------------------
# 1. Full-width covariance block admission when links count > row cap:
#    admitted 1-block plan enables the retained store (1
#    producer pass, 3 slice passes) at the exact legacy row partition and
#    matches the forced legacy multi-block planner across the full centered
#    Gram, IF, active Tikhonov ridge, inference and all 4 aggregations.
# ---------------------------------------------------------------------------
test_that("admitted full-width covariance block when L > chunk_rows enables existing complete-support cache and matches forced legacy planner at 1e-12", {
  expect_true(environmentIsLocked(api))
  # Single-cohort unreduced M = 10 > chunk_rows = 7 (with and without X, active ridge).
  n_chunks_rm <- as.integer(ceiling(41 / 7))
  rows_rm <- split(seq_len(41L), ceiling(seq_len(41L) / 7L))
  for (est in c("dr", "reg")) for (xf in list(NULL, ~x)) {
    e_new <- plan_counting_ns()
    e_leg <- plan_counting_ns(legacy_plan = TRUE)
    e_str <- plan_counting_ns(force_stream = TRUE)
    res_new <- rm_run_pipeline(rm_d, est = est, chunk = 7L, reduce = FALSE, xformla = xf,
      condition_trigger = 0.99, ridge_relative = 1e-2, ns = e_new)
    res_leg <- rm_run_pipeline(rm_d, est = est, chunk = 7L, reduce = FALSE, xformla = xf,
      condition_trigger = 0.99, ridge_relative = 1e-2, ns = e_leg)
    res_str <- rm_run_pipeline(rm_d, est = est, chunk = 7L, reduce = FALSE, xformla = xf,
      condition_trigger = 0.99, ridge_relative = 1e-2, ns = e_str)
    b <- res_new$nuisance_registry$cp04
    ids <- vapply(res_new$point_links_used, `[[`, integer(1), "cp04_index")
    expect_identical(length(ids), 10L)
    expect_gt(length(ids), b$control$max_chunk_rows)
    plan_new <- api$cp04_stream_plan(b, ids, 0)
    plan_leg <- api$cp04_stream_plan_legacy(b, ids, 0)
    expect_identical(length(plan_new$blocks), 1L)
    expect_identical(as.integer(plan_new$blocks[[1L]]), seq_along(ids))
    expect_identical(length(plan_leg$blocks), 2L)
    expect_identical(plan_new$chunk_rows, plan_leg$chunk_rows)
    expect_identical(plan_new$chunk_rows, 7)
    expect_lte(plan_new$work, plan_leg$work)
    rec <- rm_record(res_new)
    expect_identical(rec$status, "retained")
    expect_identical(rec$producer_passes, 1L)
    expect_identical(rec$slice_passes, 3L)
    expect_identical(rec$n_links, 10L)
    expect_identical(rec$chunk_rows, 7)
    expect_identical(e_new$n$replays, 4L)
    expect_identical(e_new$n$chunks, 1L * n_chunks_rm)
    expect_identical(unname(e_new$n$rows), unname(rows_rm))
    # Forced legacy planner splits 10 links into 2 blocks (2 mean + 3 Gram + 2 GMM = 7 replays).
    expect_identical(rm_record(res_leg)$status, "streamed")
    expect_identical(e_leg$n$replays, 7L)
    expect_identical(e_leg$n$chunks, 7L * n_chunks_rm)
    # Forced stream with new 1-block planner runs 4 replays (1 mean + 1 Gram + 2 GMM).
    expect_identical(rm_record(res_str)$status, "streamed")
    expect_identical(e_str$n$replays, 4L)
    expect_identical(e_str$n$chunks, 4L * n_chunks_rm)
    expect_identical(res_new$covariance_handling$handler_used, "tikhonov")
    expect_gt(res_new$covariance_handling$lambda, 0)
    cp09_equal(res_new, res_leg)
    cp09_equal(res_new, res_str)
    expect_false(rm_find_store(res_new))
  }

  # Multi-cohort public catalog: NT has 2 * choose(6, 2) = 30 links.
  # NYT also has the (26, Inf) control pool for g = 9 on the choose(4, 2) = 6
  # pairs among dates 2, 4, 9, 15 (all before 26), hence L = 36.
  # Both exceed chunk_rows = 9: legacy has 4 blocks (14 cov replays + 2 GMM = 16);
  # New planner admits 1 block at chunk_rows = 9. All supported sealed
  # baselines retain (1 producer, 3 slices), including general mixed-control NYT OR.
  n_chunks_cp09 <- as.integer(ceiling(53 / 9))
  for (est in c("dr", "reg")) for (cg in c("nevertreated", "notyettreated")) {
    e_new <- plan_counting_ns()
    e_leg <- plan_counting_ns(legacy_plan = TRUE)
    e_gen <- plan_counting_ns(decline_complete = TRUE)
    a <- cp09_run(est = est, cg = cg, chunk = 9L, condition_trigger = 0.99, ridge_relative = 1e-2,
      clustered = FALSE, ns = e_new)
    leg <- cp09_run(est = est, cg = cg, chunk = 9L, condition_trigger = 0.99, ridge_relative = 1e-2,
      clustered = FALSE, ns = e_leg)
    gen <- cp09_run(est = est, cg = cg, chunk = 9L, condition_trigger = 0.99, ridge_relative = 1e-2,
      clustered = FALSE, ns = e_gen)
    b <- a$nuisance_registry$cp04
    st <- get("cp04_reader", b$reader_ref, inherits = FALSE)$state
    ids <- vapply(a$point_links_used, `[[`, integer(1), "cp04_index")
    ra <- a$resource_accounting
    expect_identical(length(ids), if (cg == "nevertreated") 30L else 36L)
    expect_gt(length(ids), b$control$max_chunk_rows)
    retained <- ra$initial_envelope_entries + ra$d5_arrays$retained_diagnostic_entries
    cov_new <- api$cp04_stream_plan(b, ids, retained)
    cov_leg <- api$cp04_stream_plan_legacy(b, ids, retained)
    expect_identical(length(cov_new$blocks), 1L)
    expect_identical(as.integer(cov_new$blocks[[1L]]), seq_along(ids))
    expect_identical(length(cov_leg$blocks), 4L)
    expect_identical(cov_new$chunk_rows, cov_leg$chunk_rows)
    expect_equal(cov_new$chunk_rows, 9, tolerance = 0)
    expect_lte(cov_new$work, cov_leg$work)
    rec <- ra$point_producer$retained_moment_if
    if (!(est == "reg" && cg == "notyettreated")) {
      expect_false(is.null(b$complete_support))
      expect_identical(rec$status, "retained")
      expect_identical(rec$producer_passes, 1L)
      expect_identical(rec$slice_passes, 3L)
      expect_identical(e_new$n$chunks, 1L * n_chunks_cp09)
    } else {
      # OR complete facts require every control pool to be exactly never-treated.
      expect_s3_class(b, "cdid_or_baseline")
      expect_true(any(vapply(b$regression_specs, function(x)
        identical(st$cohorts[x$cohorts], c(26, Inf)), logical(1))))
      expect_null(b$complete_support)
      rm_expect_retained(a)
      expect_true(is.environment(api$cp04_retained_if_reserve(b, ids, cov_new$chunk_rows, cov_new,
        retained = 0, gmm_retained = 0, live = 0, work = 0, control = b$control)))
      expect_identical(e_new$n$chunks, 1L * n_chunks_cp09)
    }
    expect_identical(e_new$n$replays, 4L)
    expect_identical(leg$resource_accounting$point_producer$retained_moment_if$status, "streamed")
    expect_identical(e_leg$n$replays, 4L + 10L + 2L)
    expect_identical(e_leg$n$chunks, 16L * n_chunks_cp09)
    # Charged covariance and GMM probes keep the exact same 1-block and row partition.
    cov_charged <- api$cp04_stream_plan(b, ids, retained + rec$reserved_entries)
    expect_identical(cov_charged$blocks, cov_new$blocks)
    expect_identical(cov_charged$chunk_rows, cov_new$chunk_rows)
    # The pipeline charges these persistent buffers to the DR GMM planner;
    # fitted baseline OR uses the original zero-retained GMM preflight.
    gmm_retained <- if (inherits(b, "cdid_dr_baseline"))
      ra$final_diagnostic_copy_entries + retained + ra$gmm_consumer$total_payload_entries +
        sum(unlist(ra$d5_arrays[c("ledger_if_entries", "basis_if_entries",
          "ledger_enclosure_entries", "basis_enclosure_entries")])) else 0
    gmm_charged <- api$cp04_preflight_links(b, ids, gmm_retained + rec$reserved_entries)
    expect_identical(gmm_charged$chunk_rows, cov_new$chunk_rows)
    # Full centered original-unit Gram, active Tikhonov ridge, IF, inference and all 4 aggregations.
    H <- cp09_phi(b, links = vapply(a$point_links_all, `[[`, integer(1), "cp04_index"))
    Hc <- sweep(H, 2L, colMeans(H), "-")
    expect_equal(unname(a$covariance_stream$V), unname(crossprod(Hc) / nrow(H)), tolerance = 1e-12)
    expect_identical(a$covariance_handling$handler_used, "tikhonov")
    expect_gt(a$covariance_handling$lambda, 0)
    expect_identical(sort(names(a$aggregations)), sort(c("simple", "group", "dynamic", "calendar")))
    cp09_equal(a, leg)
    cp09_equal(a, gen)
    expect_false(rm_find_store(a))
  }
})

# ---------------------------------------------------------------------------
# 2. D92 eligibility: missing/general, clustered, weighted and frame-only
#    sealed OR/DR routes retain at admitted capacity, with the original
#    full-IF/covariance/inference references and legacy partitions.
# ---------------------------------------------------------------------------
test_that("missing and general OR/DR plus weight, cluster and frame routes retain at admitted capacity", {
  # Missing panel observations on multi-cohort NYT OR/DR (L = 36 > chunk = 9):
  # complete_support is NULL; retention still uses 1 producer and 3 slices vs 16 legacy producer replays.
  d <- cp09_panel()
  dm <- d[!((d$id == 1053 & d$period == 4) | (d$id == 1052 & d$period == 15)), ]
  n_chunks_cp09 <- as.integer(ceiling(53 / 9))
  for (est in c("dr", "reg")) {
    e_miss_new <- plan_counting_ns()
    e_miss_leg <- plan_counting_ns(legacy_plan = TRUE)
    res_miss_new <- cp09_run(dm, est = est, cg = "notyettreated", chunk = 9L,
      condition_trigger = 0.99, ridge_relative = 1e-2, clustered = FALSE, ns = e_miss_new)
    res_miss_leg <- cp09_run(dm, est = est, cg = "notyettreated", chunk = 9L,
      condition_trigger = 0.99, ridge_relative = 1e-2, clustered = FALSE, ns = e_miss_leg)
    b_miss <- res_miss_new$nuisance_registry$cp04
    ids_miss <- vapply(res_miss_new$point_links_used, `[[`, integer(1), "cp04_index")
    expect_null(b_miss$complete_support)
    plan_miss <- api$cp04_stream_plan(b_miss, ids_miss, 0)
    plan_miss_leg <- api$cp04_stream_plan_legacy(b_miss, ids_miss, 0)
    expect_identical(length(plan_miss$blocks), 1L)
    expect_identical(length(plan_miss_leg$blocks), 4L)
    expect_identical(plan_miss$chunk_rows, plan_miss_leg$chunk_rows)
    expect_lte(plan_miss$work, plan_miss_leg$work)
    rm_expect_retained(res_miss_new)
    expect_true(is.environment(api$cp04_retained_if_reserve(b_miss, ids_miss, plan_miss$chunk_rows, plan_miss,
      retained = 0, gmm_retained = 0, live = 0, work = 0, control = b_miss$control)))
    expect_identical(e_miss_new$n$replays, 4L)
    expect_identical(e_miss_new$n$chunks, 1L * n_chunks_cp09)
    expect_identical(unname(e_miss_new$n$rows), unname(split(seq_len(53L), ceiling(seq_len(53L) / 9L))))
    expect_identical(e_miss_leg$n$replays, 16L)
    expect_identical(e_miss_leg$n$chunks, 16L * n_chunks_cp09)
    cp09_equal(res_miss_new, res_miss_leg)
    H <- cp09_phi(b_miss, links = vapply(res_miss_new$point_links_all, `[[`, integer(1), "cp04_index"))
    Hc <- sweep(H, 2L, colMeans(H), "-")
    expect_equal(unname(res_miss_new$covariance_stream$V), unname(crossprod(Hc) / nrow(H)), tolerance = 1e-12)
    expect_false(rm_find_store(res_miss_new))

    # General OR/DR route (forced via decline_complete = TRUE) on rm_d (L = 10 > chunk = 7):
    # 1 covariance block (1 producer, 3 slices) vs 2 blocks (7 producer replays) under legacy planner.
    n_chunks_rm <- as.integer(ceiling(41 / 7))
    e_gen_new <- plan_counting_ns(decline_complete = TRUE)
    e_gen_leg <- plan_counting_ns(legacy_plan = TRUE, decline_complete = TRUE)
    res_gen_new <- rm_run_pipeline(rm_d, est = est, chunk = 7L, xformla = ~x, ns = e_gen_new)
    res_gen_leg <- rm_run_pipeline(rm_d, est = est, chunk = 7L, xformla = ~x, ns = e_gen_leg)
    expect_null(res_gen_new$nuisance_registry$cp04$complete_support)
    rm_expect_retained(res_gen_new)
    expect_identical(e_gen_new$n$replays, 4L)
    expect_identical(e_gen_new$n$chunks, 1L * n_chunks_rm)
    expect_identical(unname(e_gen_new$n$rows), unname(split(seq_len(41L), ceiling(seq_len(41L) / 7L))))
    expect_identical(e_gen_leg$n$replays, 7L)
    expect_identical(e_gen_leg$n$chunks, 7L * n_chunks_rm)
    cp09_equal(res_gen_new, res_gen_leg)
  }

  # Newly admitted OR/DR routes keep complete numerical parity with the
  # forced legacy planner at 1e-12.
  for (est in c("dr", "reg")) {
    res_cl_new <- cp09_run(d, est = est, cg = "notyettreated", chunk = 9L, clustered = TRUE)
    res_cl_leg <- cp09_run(d, est = est, cg = "notyettreated", chunk = 9L, clustered = TRUE,
      ns = plan_counting_ns(legacy_plan = TRUE))
    rm_expect_retained(res_cl_new)
    cp09_equal(res_cl_new, res_cl_leg)

    for (w in list(rep(2.5, 53L), 1 + (seq_len(53L) %% 5) / 8)) {
      dw <- d; dw$w <- w[match(dw$id, sort(unique(dw$id)))]
      aa <- cp09_run(dw, est = est, weightsname = "w", chunk = 9L, clustered = FALSE)
      ll <- cp09_run(dw, est = est, weightsname = "w", chunk = 9L, clustered = FALSE,
        ns = plan_counting_ns(legacy_plan = TRUE))
      rm_expect_retained(aa)
      expect_false(rm_find_store(aa))
      cp09_equal(aa, ll)
    }

    fr <- rm_d[!duplicated(rm_d$id), c("id", "G", "x", "cl")]
    fr <- rbind(fr[order(fr$id), ], data.frame(id = 9999L, G = 0, x = 0, cl = 99L))
    res_fr <- rm_run_pipeline(rm_d, est = est, chunk = 7L, frame = fr, aggregations = "simple")
    expect_null(res_fr$nuisance_registry$cp04$complete_support)
    rm_expect_retained(res_fr)
    expect_false(rm_find_store(res_fr))
  }
})

# ---------------------------------------------------------------------------
# 3. Exact adjacent entry/output/work capacity decline, charged planner
#    partition declines, integrity propagation and zero callbacks before
#    preflight.
# ---------------------------------------------------------------------------
test_that("cp04_stream_plan and cp04_retained_if_reserve decline at adjacent entry, output, work and partition boundaries and invoke zero callbacks on preflight refusal", {
  fx <- cp09_fixture(cg = "notyettreated", est = "dr", chunk = 9L)
  fit <- cp09_fit(fx)
  links <- seq_along(fit$links)
  L <- length(links); B <- 9; N <- 53L
  expect_identical(L, 30L)
  legacy0 <- api$cp04_stream_plan_legacy(fit, links, 0)
  expect_identical(length(legacy0$blocks), 4L)
  expect_identical(legacy0$chunk_rows, B)

  pre_full <- api$chained_full_if_budget(fit, B, links)
  producer_peak <- pre_full$output + pre_full$scratch
  # Independent live-buffer and owner/header inventory, in addition to the old four-Gram allowance.
  inventory <- rm_retained_inventory(fit, links, B)
  consumer <- 6 * B * L + 2 * L * L + 16 * L + 256 + inventory$metadata
  peak_full <- producer_peak + consumer
  own <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  gram <- 4 * L * L
  ctrl <- fit$control

  # 3a. Planner adjacent simultaneous entry cap via `retained` argument:
  # admitted at the exact cap, declines to legacy multi-block one entry above.
  ret_cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries - peak_full - gram) - own
  plan_at <- api$cp04_stream_plan(fit, links, ret_cap)
  plan_above <- api$cp04_stream_plan(fit, links, ret_cap + 1)
  expect_identical(length(plan_at$blocks), 1L)
  expect_identical(plan_at$chunk_rows, B)
  expect_identical(plan_at$peak, peak_full)
  expect_identical(plan_at$consumer_workspace_entries, consumer)
  full_work <- 2 * ceiling(N / B) * pre_full$work + 4 * N * L * L + inventory$boundary_work
  expect_identical(plan_at$work, full_work)
  expect_identical(plan_above, api$cp04_stream_plan_legacy(fit, links, ret_cap + 1))
  expect_identical(length(plan_above$blocks), 4L)

  # 3b. Planner adjacent `max_matrix_entries` cap on resealed backend (retained = 0):
  reseal_ctrl <- function(base_fit, ...) {
    b <- base_fit; b$ref <- NULL
    mods <- list(...)
    for (nm in names(mods)) b$control[[nm]] <- mods[[nm]]
    api$ipw_owner_seal(b, "dr_baseline_backend")
  }
  mat_cap <- own + peak_full + gram
  fit_mat_at <- reseal_ctrl(fit, max_matrix_entries = mat_cap, max_retained_entries = mat_cap)
  fit_mat_below <- reseal_ctrl(fit, max_matrix_entries = mat_cap - 1, max_retained_entries = mat_cap - 1)
  expect_identical(length(api$cp04_stream_plan(fit_mat_at, links, 0)$blocks), 1L)
  expect_identical(api$cp04_stream_plan(fit_mat_below, links, 0),
    api$cp04_stream_plan_legacy(fit_mat_below, links, 0))
  expect_identical(length(api$cp04_stream_plan(fit_mat_below, links, 0)$blocks), 4L)

  # 3c. Planner adjacent `max_output_entries` cap on resealed backend:
  out_cap <- max(pre_full$output, gram)
  expect_gt(pre_full$output, gram)
  fit_out_at <- reseal_ctrl(fit, max_output_entries = out_cap)
  fit_out_below <- reseal_ctrl(fit, max_output_entries = out_cap - 1)
  expect_identical(length(api$cp04_stream_plan(fit_out_at, links, 0)$blocks), 1L)
  expect_identical(api$cp04_stream_plan(fit_out_below, links, 0),
    api$cp04_stream_plan_legacy(fit_out_below, links, 0))
  expect_identical(length(api$cp04_stream_plan(fit_out_below, links, 0)$blocks), 4L)

  # 3d. Planner adjacent `max_work_ops` cap:
  # Each legacy tile pair has <= 18 links (< 30), so its chunk work is strictly below `pre_full$work`.
  fit_w <- cp09_fit(fx, external_work_ops = legacy0$work)
  work_cap <- max(fit_w$accounting$work_ops_bound + pre_full$work, full_work)
  expect_gt(work_cap, legacy0$work)
  fit_w_at <- reseal_ctrl(fit_w, max_work_ops = work_cap)
  fit_w_below <- reseal_ctrl(fit_w, max_work_ops = work_cap - 1)
  expect_identical(length(api$cp04_stream_plan(fit_w_at, links, 0)$blocks), 1L)
  expect_identical(api$cp04_stream_plan(fit_w_below, links, 0),
    api$cp04_stream_plan_legacy(fit_w_below, links, 0))
  expect_identical(length(api$cp04_stream_plan(fit_w_below, links, 0)$blocks), 4L)

  # Full proposal work > legacy work is declined even when within max_work_ops.
  e_work_guard <- new.env(parent = api)
  for (nm in c("cp04_stream_plan", "cp04_stream_plan_full_block")) {
    f <- get(nm, api); environment(f) <- e_work_guard; assign(nm, f, e_work_guard)
  }
  e_work_guard$cp04_stream_plan_legacy <- function(backend, links, retained = 0) {
    out <- api$cp04_stream_plan_legacy(backend, links, retained)
    out$work <- plan_at$work - 1
    out
  }
  expect_identical(length(e_work_guard$cp04_stream_plan(fit, links, 0)$blocks), 4L)

  # 3e. Automatic row-partition reduction decline: when the full link width
  # would reduce B below the legacy row partition, the 1-block proposal declines.
  fit_auto <- cp09_fit(fx, automatic_chunk_rows = TRUE)
  own_auto <- fit_auto$accounting$owned_entries + fit_auto$accounting$borrowed_entries
  one_full <- api$dr_baseline_if_envelope(fit_auto, 1, links)
  two_full <- api$dr_baseline_if_envelope(fit_auto, 2, links)
  slope_full <- two_full$output + two_full$scratch - one_full$output - one_full$scratch
  icpt_full <- one_full$output + one_full$scratch - slope_full
  # At `mat_auto_at`, full-width `cp04_plan_rows` gives B = 9; one entry lower gives B_full = 8 < 9
  # while every <=18-link tile pair in legacy still keeps B = 9.
  mat_auto_at <- own_auto + gram + consumer + icpt_full + B * slope_full
  fit_auto_at <- reseal_ctrl(fit_auto, max_matrix_entries = mat_auto_at, max_retained_entries = mat_auto_at)
  fit_auto_shrink <- reseal_ctrl(fit_auto, max_matrix_entries = mat_auto_at - 1, max_retained_entries = mat_auto_at - 1)
  expect_identical(api$cp04_plan_rows(fit_auto_at, list(links), gram + consumer), B)
  expect_identical(api$cp04_plan_rows(fit_auto_shrink, list(links), gram + consumer), B - 1)
  expect_identical(length(api$cp04_stream_plan(fit_auto_at, links, 0)$blocks), 1L)
  plan_shrink <- api$cp04_stream_plan(fit_auto_shrink, links, 0)
  expect_identical(plan_shrink$chunk_rows, B)
  expect_identical(length(plan_shrink$blocks), 4L)

  # The exact issued B is the admission target, even if an automatic planner
  # would suggest another partition. The explicit preflight is still real.
  clone_planner <- function() {
    e <- new.env(parent = api)
    for (nm in c("cp04_stream_plan", "cp04_stream_plan_full_block", "cp04_retained_if_validate_backend")) {
      f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
    }
    e
  }
  e_fixed <- clone_planner()
  auto_calls <- 0L; explicit_calls <- 0L
  e_fixed$cp04_plan_rows <- function(backend, link_sets, retained = 0, chunk_rows = NULL) {
    if (is.null(chunk_rows)) { auto_calls <<- auto_calls + 1L; return(B - 1) }
    explicit_calls <<- explicit_calls + 1L
    api$cp04_plan_rows(backend, link_sets, retained, chunk_rows)
  }
  expect_identical(length(e_fixed$cp04_stream_plan(fit, links, 0)$blocks), 1L)
  expect_identical(auto_calls, 0L)
  expect_identical(explicit_calls, 1L)

  # 3f. Retained store (`cp04_retained_if_reserve`) adjacent entry/work caps and
  # charged covariance/GMM planner probes when L = 30 > B = 9.
  plan0 <- api$cp04_stream_plan(fit, links, 0)
  expect_identical(length(plan0$blocks), 1L)
  entries <- inventory$entries
  cache_work <- inventory$work
  cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  reserve <- function(..., backend = fit, links_ = links, B_ = B, plan_ = plan0,
                      retained = 0, gmm_retained = 0, live = 0, work = 0, control = ctrl)
    api$cp04_retained_if_reserve(backend, links_, B_, plan_, retained = retained,
      gmm_retained = gmm_retained, live = live, work = work, control = control)

  at_live <- reserve(live = cap - entries)
  expect_true(is.environment(at_live))
  expect_identical(at_live$admission_total, cap)
  expect_identical(at_live$entries, inventory$entries)
  expect_identical(at_live$metadata_entries, inventory$metadata)
  expect_identical(at_live$boundary_work, inventory$boundary_work)
  expect_identical(at_live$work, inventory$work)
  expect_false(exists("phi", envir = at_live, inherits = FALSE))
  expect_null(reserve(live = cap - entries + 1))

  expect_true(is.environment(reserve(work = ctrl$max_work_ops - cache_work)))
  expect_null(reserve(work = ctrl$max_work_ops - cache_work + 1))

  # Charged covariance probe: at `room_cov`, `retained + entries` still admits 1 block;
  # at `room_cov + 1`, uncharged `cp04_stream_plan` admits 1 block, but the charged probe
  # (`retained + 1 + entries`) falls back to the 4-block legacy plan and declines the store.
  room_cov <- ret_cap - entries
  plan_cov_at <- api$cp04_stream_plan(fit, links, room_cov)
  plan_cov_above <- api$cp04_stream_plan(fit, links, room_cov + 1)
  expect_identical(length(plan_cov_at$blocks), 1L)
  expect_identical(length(plan_cov_above$blocks), 1L)
  expect_true(is.environment(reserve(plan_ = plan_cov_at, retained = room_cov)))
  expect_null(reserve(plan_ = plan_cov_above, retained = room_cov + 1))

  # Charged GMM probe: admitted at `room_gmm`, declined at `room_gmm + 1`.
  room_gmm <- ctrl$max_matrix_entries - own - producer_peak - entries
  expect_true(is.environment(reserve(gmm_retained = room_gmm)))
  expect_null(reserve(gmm_retained = room_gmm + 1))

  # Automatic charged GMM preflight changes rows rather than throwing; the
  # same store must decline in this case too, without allocation.
  plan_auto <- api$cp04_stream_plan(fit_auto, links, 0)
  expect_identical(length(plan_auto$blocks), 1L)
  pre_auto <- api$chained_full_if_budget(fit_auto, B, links)
  room_auto <- fit_auto$control$max_matrix_entries - own_auto - pre_auto$output - pre_auto$scratch - entries
  gmm_auto_at <- api$cp04_preflight_links(fit_auto, links, room_auto + entries)
  gmm_auto_below <- api$cp04_preflight_links(fit_auto, links, room_auto + entries + 1)
  expect_identical(gmm_auto_at$chunk_rows, B)
  expect_identical(gmm_auto_below$chunk_rows, B - 1)
  expect_true(is.environment(reserve(backend = fit_auto, plan_ = plan_auto, gmm_retained = room_auto)))
  expect_null(reserve(backend = fit_auto, plan_ = plan_auto, gmm_retained = room_auto + 1))

  # Reference inputs are eligible only on a valid sealed backend. This
  # unsealed reference mutation must fail the original validator first.
  ref_fit <- fit; ref_fit$reference_ref <- new.env(parent = emptyenv())
  expect_error(reserve(backend = ref_fit), class = "cdid_dr_baseline_integrity_error")

  # Partition / ordered-links / legacy-plan mismatches decline or refuse at boundary.
  expect_null(reserve(B_ = B - 1))
  expect_null(reserve(plan_ = legacy0))
  expect_null(reserve(links_ = links[-1L]))
  expect_error(reserve(links_ = c(1L, 1L, links[-c(1L, 2L)])), class = "cdid_cp04_input_error")
  st_ord <- reserve()
  bound_ord <- api$cp04_retained_if_bind(fit, st_ord)
  expect_error(api$cp04_replay_links(bound_ord, rev(links), function(...) NULL, chunk_rows = B),
    class = "cdid_cp04_integrity_error")
  expect_identical(st_ord$state, "released")
  st_ord <- reserve()
  bound_ord <- api$cp04_retained_if_bind(fit, st_ord)
  expect_error(api$cp04_replay_links(bound_ord, links, function(...) NULL, chunk_rows = B - 1),
    class = "cdid_cp04_integrity_error")

  expect_identical(st_ord$state, "released")

  # Pipeline optional cap at L = 30 > B = 9 flips from retained to 1-block streamed one entry below cap.
  res_pipe <- cp09_run(cg = "notyettreated", chunk = 9L, clustered = FALSE)
  rec_pipe <- rm_record(res_pipe)
  expect_identical(rec_pipe$status, "retained")
  pipe_backend <- res_pipe$nuisance_registry$cp04
  pipe_links <- vapply(res_pipe$point_links_used, `[[`, integer(1), "cp04_index")
  pipe_inventory <- rm_retained_inventory(pipe_backend, pipe_links, B)
  pipe_pre <- api$chained_full_if_budget(pipe_backend, B, pipe_links)
  pipe_L <- length(pipe_links)
  pipe_consumer <- 6 * B * pipe_L + 2 * pipe_L^2 + 16 * pipe_L + 256 + pipe_inventory$metadata
  pipe_peak <- pipe_backend$accounting$owned_entries + pipe_backend$accounting$borrowed_entries +
    pipe_pre$output + pipe_pre$scratch + pipe_consumer
  pipe_ra <- res_pipe$resource_accounting
  pipe_aux <- 8 * N * (length(res_pipe$nuisance_registry$cp04_aux$selections[["varying"]]) + nrow(res_pipe$support$ledger))
  pipe_total <- pipe_ra$simultaneous_entries + pipe_ra$final_diagnostic_copy_entries +
    max(0, pipe_peak - pipe_ra$point_producer$producer_workspace_entries) + pipe_aux + pipe_inventory$entries
  expect_identical(rec_pipe$reserved_entries, pipe_inventory$entries)
  expect_identical(rec_pipe$snapshot_metadata_entries, pipe_inventory$metadata)
  expect_identical(rec_pipe$boundary_check_work_ops, pipe_inventory$boundary_work)
  expect_identical(rec_pipe$work_ops_bound, pipe_inventory$work)
  expect_identical(rec_pipe$admission_total_entries, pipe_total)
  res_pipe_at <- cp09_run(cg = "notyettreated", chunk = 9L, clustered = FALSE,
    ns = plan_counting_ns(cap = pipe_total))
  e_pipe_below <- plan_counting_ns(cap = pipe_total - 1)
  res_pipe_below <- cp09_run(cg = "notyettreated", chunk = 9L, clustered = FALSE, ns = e_pipe_below)
  expect_identical(rm_record(res_pipe_at)$status, "retained")
  expect_identical(rm_record(res_pipe_below)$status, "streamed")
  expect_identical(e_pipe_below$n$replays, 4L)
  expect_identical(e_pipe_below$n$chunks, 4L * as.integer(ceiling(N / B)))
  expect_identical(unname(e_pipe_below$n$rows), rep(unname(split(seq_len(N), ceiling(seq_len(N) / B))), 4L))
  cp09_equal(res_pipe, res_pipe_at)
  cp09_equal(res_pipe, res_pipe_below)

  # 3g. Integrity errors propagate (never swallowed as capacity decline).
  broken <- fit; broken$links <- broken$links[1:29]
  expect_error(api$cp04_stream_plan(broken, 1:29, 0), class = "cdid_dr_baseline_integrity_error")
  expect_error(reserve(backend = broken), class = "cdid_dr_baseline_integrity_error")

  # A one-block legacy plan keeps its original fast path: the prospective
  # full-width validator is reached only for a multi-block proposal.
  e_fast <- clone_planner()
  e_fast$validate_dr_baseline_backend <- function(...) api$cp04_abort("validation sentinel", "integrity")
  expect_identical(e_fast$cp04_stream_plan(fit, 1:4, 0), api$cp04_stream_plan_legacy(fit, 1:4, 0))
  expect_error(e_fast$cp04_stream_plan(fit, links, 0), class = "cdid_cp04_integrity_error")

  # Real supplied budget refusals may select legacy. Input, numerical,
  # integrity and bound-underestimation defects must propagate unchanged.
  e_error <- clone_planner()
  e_error$chained_full_if_budget <- function(...) {
    api$cp04_budget(ctrl, 0, 0, 0, ctrl$max_output_entries + 1, "forced output cap")
  }
  expect_identical(e_error$cp04_stream_plan(fit, links, 0), legacy0)
  for (kind in c("input", "numeric", "integrity", "bound")) {
    e_error$chained_full_if_budget <- local({
      k <- kind
      function(...) api$cp04_abort("forced non-capacity defect", k)
    })
    expect_error(e_error$cp04_stream_plan(fit, links, 0), class = paste0("cdid_cp04_", kind, "_error"))
  }
  e_error$chained_full_if_budget <- function(...) api$ipw_owner_int_bound(-1, "forced invalid input")
  expect_error(e_error$cp04_stream_plan(fit, links, 0), class = "cdid_ipw_owner_control_error")
  expect_error(api$cp04_stream_plan(fit, replace(links, 1L, links[2L]), 0),
    class = "cdid_dr_baseline_input_error")
  expect_error(api$cp04_stream_plan(fit, links, -1), class = "cdid_ipw_owner_control_error")

  # 3h. Zero callbacks and zero chunk evaluations before preflight refusal.
  e_zero <- plan_counting_ns()
  cb_calls <- 0L
  cb <- function(...) cb_calls <<- cb_calls + 1L
  expect_error(e_zero$cp04_replay_links(fit, links, cb, chunk_rows = B + 1L),
    class = "cdid_ipw_owner_bound_error")
  expect_identical(cb_calls, 0L)
  expect_identical(e_zero$n$chunks, 0L)
  expect_error(e_zero$cp04_stream_covariance(fit, links, retained = ctrl$max_matrix_entries),
    class = "cdid_ipw_owner_bound_error")
  expect_identical(e_zero$n$chunks, 0L)

  changed_store <- reserve()
  changed_bound <- api$cp04_retained_if_bind(fit, changed_store)
  expect_error(e_zero$cp04_stream_covariance(changed_bound, links, retained = ret_cap - entries + 1),
    class = "cdid_ipw_owner_bound_error")
  expect_identical(e_zero$n$chunks, 0L)
  expect_false(exists("phi", envir = changed_store, inherits = FALSE))
  expect_identical(changed_store$state, "released")
  changed_store <- reserve()
  changed_bound <- api$cp04_retained_if_bind(fit, changed_store)
  expect_error(e_zero$cp04_stream_covariance(changed_bound, links, retained = ret_cap + 1),
    class = "cdid_cp04_integrity_error")
  expect_identical(e_zero$n$chunks, 0L)
  expect_false(exists("phi", envir = changed_store, inherits = FALSE))

  expect_identical(changed_store$state, "released")

  # Ordered links, complete centered cross-terms and original row partition
  # also hold for a reversed link request against a genuinely forced legacy run.
  ordered <- rev(links)
  e_order <- plan_counting_ns(force_stream = TRUE)
  e_order_leg <- plan_counting_ns(legacy_plan = TRUE)
  cov_order <- e_order$cp04_stream_covariance(fit, ordered)
  cov_order_leg <- e_order_leg$cp04_stream_covariance(fit, ordered)
  for (nm in c("estimates", "covariance", "mean_phi"))
    expect_equal(cov_order[[nm]], cov_order_leg[[nm]], tolerance = 1e-12)
  expect_identical(cov_order$chunk_rows, cov_order_leg$chunk_rows)
  expect_identical(e_order$n$replays, 2L)
  expect_identical(e_order_leg$n$replays, 14L)
})

test_that("fitted baseline OR full-width planning retains real cap fallbacks, NYT admission and OR errors", {
  fx <- cp09_fixture(cg = "notyettreated", est = "or", chunk = 9L)
  fit <- cp09_fit(fx)
  links <- seq_along(fit$links)
  N <- 53L; B <- 9; L <- length(links); ctrl <- fit$control
  expect_identical(L, 30L)
  expect_true(any(vapply(fx$links, function(x) identical(x$K, c(26, Inf)), logical(1))))
  expect_null(fit$complete_support)
  plan <- api$cp04_stream_plan(fit, links, 0)
  legacy <- api$cp04_stream_plan_legacy(fit, links, 0)
  expect_identical(length(plan$blocks), 1L)
  expect_identical(length(legacy$blocks), 4L)
  expect_identical(plan$chunk_rows, legacy$chunk_rows)
  expect_identical(plan$chunk_rows, B)
  expect_lte(plan$work, legacy$work)

  pre <- api$chained_full_if_budget(fit, B, links)
  producer <- pre$output + pre$scratch
  inventory <- rm_retained_inventory(fit, links, B)
  consumer <- 6 * B * L + 2 * L * L + 16 * L + 256 + inventory$metadata
  gram <- 4 * L * L
  own <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  expect_identical(plan$peak, producer + consumer)
  expect_identical(plan$consumer_workspace_entries, consumer)
  ret_cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries - plan$peak - gram) - own
  expect_identical(length(api$cp04_stream_plan(fit, links, ret_cap)$blocks), 1L)
  expect_identical(api$cp04_stream_plan(fit, links, ret_cap + 1),
    api$cp04_stream_plan_legacy(fit, links, ret_cap + 1))

  reseal <- function(...) {
    b <- fit; b$ref <- NULL
    mods <- list(...)
    for (nm in names(mods)) b$control[[nm]] <- mods[[nm]]
    api$ipw_owner_seal(b, "or_baseline_backend")
  }
  # These are actual adjacent capacity refusals, not class-based exclusions.
  mat_cap <- own + plan$peak + gram
  at <- reseal(max_matrix_entries = mat_cap, max_retained_entries = mat_cap)
  below <- reseal(max_matrix_entries = mat_cap - 1, max_retained_entries = mat_cap - 1)
  expect_identical(length(api$cp04_stream_plan(at, links, 0)$blocks), 1L)
  expect_identical(api$cp04_stream_plan(below, links, 0), api$cp04_stream_plan_legacy(below, links, 0))
  expect_identical(length(api$cp04_stream_plan(below, links, 0)$blocks), 4L)
  expect_gt(pre$output, gram)
  out_at <- reseal(max_output_entries = pre$output)
  out_below <- reseal(max_output_entries = pre$output - 1)
  expect_identical(length(api$cp04_stream_plan(out_at, links, 0)$blocks), 1L)
  expect_identical(api$cp04_stream_plan(out_below, links, 0), api$cp04_stream_plan_legacy(out_below, links, 0))

  # Mixed-control NYT OR is general and now eligible. Exact entry/work caps
  # and charged planner probes still enforce the complete inventory.
  entries <- inventory$entries
  cache_work <- inventory$work
  reserve <- function(retained = 0, gmm_retained = 0, live = 0, work = 0, backend = fit)
    api$cp04_retained_if_reserve(backend, links, B, plan, retained = retained,
      gmm_retained = gmm_retained, live = live, work = work, control = ctrl)
  cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  expect_true(is.environment(reserve()))
  expect_true(is.environment(reserve(live = cap - entries)))
  expect_null(reserve(live = cap - entries + 1))
  expect_true(is.environment(reserve(work = ctrl$max_work_ops - cache_work)))
  expect_null(reserve(work = ctrl$max_work_ops - cache_work + 1))
  room_cov <- ret_cap - entries
  expect_gte(room_cov, 0)
  expect_identical(length(api$cp04_stream_plan(fit, links, room_cov + 1)$blocks), 1L)
  expect_identical(length(api$cp04_stream_plan(fit, links, room_cov + entries)$blocks), 1L)
  expect_identical(api$cp04_stream_plan(fit, links, room_cov + entries + 1),
    api$cp04_stream_plan_legacy(fit, links, room_cov + entries + 1))
  expect_true(is.environment(reserve(retained = room_cov)))
  expect_null(reserve(retained = room_cov + 1))
  room_gmm <- ctrl$max_matrix_entries - own - producer - entries
  expect_true(is.environment(reserve(gmm_retained = room_gmm)))
  expect_null(reserve(gmm_retained = room_gmm + 1))
  ref_fit <- fit; ref_fit$reference_ref <- new.env(parent = emptyenv())
  expect_error(reserve(backend = ref_fit), class = "cdid_or_baseline_integrity_error")

  # Ordered full-width streamed covariance agrees with the forced legacy
  # blocks and the complete detailed original-unit IF, with no cached rows.
  ordered <- rev(links)
  e_new <- plan_counting_ns(force_stream = TRUE)
  e_leg <- plan_counting_ns(legacy_plan = TRUE, force_stream = TRUE)
  a <- e_new$cp04_stream_covariance(fit, ordered)
  b <- e_leg$cp04_stream_covariance(fit, ordered)
  for (nm in c("estimates", "covariance", "mean_phi"))
    expect_equal(a[[nm]], b[[nm]], tolerance = 1e-12)
  H <- cp09_phi(fit, links = ordered)
  Hd <- cp09_phi(fit, links = ordered, details = TRUE)
  expect_equal(H, Hd, tolerance = 1e-12)
  Hc <- sweep(Hd, 2L, colMeans(Hd), "-")
  expect_equal(unname(a$covariance), unname(crossprod(Hc) / N), tolerance = 1e-12)
  expect_equal(unname(a$mean_phi), unname(colMeans(Hd)), tolerance = 1e-12)
  expect_identical(e_new$n$replays, 2L)
  expect_identical(e_leg$n$replays, 14L)

  clone_planner <- function() {
    e <- new.env(parent = api)
    for (nm in c("cp04_stream_plan", "cp04_stream_plan_full_block", "cp04_retained_if_validate_backend")) {
      f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
    }
    e
  }
  e_fast <- clone_planner()
  e_fast$validate_or_baseline_backend <- function(...) api$cp04_abort("OR validation sentinel", "integrity")
  expect_identical(e_fast$cp04_stream_plan(fit, 1:4, 0), api$cp04_stream_plan_legacy(fit, 1:4, 0))
  expect_error(e_fast$cp04_stream_plan(fit, links, 0), class = "cdid_cp04_integrity_error")
  broken <- fit; broken$links <- broken$links[-L]
  expect_error(api$cp04_stream_plan(broken, links[-L], 0), class = "cdid_or_baseline_integrity_error")
  expect_error(api$cp04_stream_plan(fit, replace(links, 1L, links[2L]), 0),
    class = "cdid_or_baseline_input_error")
  expect_error(api$cp04_stream_plan(fit, links, -1), class = "cdid_ipw_owner_control_error")
  e_work <- clone_planner()
  e_work$cp04_stream_plan_legacy <- function(backend, links, retained = 0) {
    p <- api$cp04_stream_plan_legacy(backend, links, retained)
    p$work <- plan$work - 1
    p
  }
  expect_identical(length(e_work$cp04_stream_plan(fit, links, 0)$blocks), 4L)
  e_zero <- plan_counting_ns()
  expect_error(e_zero$cp04_stream_covariance(fit, links, retained = ctrl$max_matrix_entries),
    class = "cdid_ipw_owner_bound_error")
  expect_identical(e_zero$n$chunks, 0L)
  expect_identical(e_zero$n$replays, 0L)
  # These classes are supported now; an empty fake must fail input validation
  # before the planner or any owner lookup, rather than decline by class.
  input_aborts <- c(cdid_or_corollary1 = "or_corollary1_abort",
    cdid_dr_corollary1 = "dr_corollary1_abort", cdid_cp04_custom_backend = "cp04_abort",
    cdid_imputation_backend = "imputation_abort", cdid_augmented_backend = "augmented_abort",
    cdid_reference_ipw_backend = "reference_ipw_abort")
  for (cl in names(input_aborts)) {
    # Each route keeps its original input condition category. Obtain that
    # category from its existing abort, independently of the new planner.
    input_error <- tryCatch(get(input_aborts[[cl]], api)("input-category witness"), error = identity)
    expect_error(api$cp04_stream_plan_full_block(structure(list(), class = cl), links),
      regexp = "Invalid explicit covariance link batch.", fixed = TRUE, class = class(input_error)[1L])
  }
})
