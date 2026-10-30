# SPDX-License-Identifier: GPL-3.0-only
#
# D90 missing-panel IPW: one invocation-local, outcome-free sparse replay state
# on the general route when raw IF retention is declined. The pipeline clones
# below force that decline to isolate this fallback. Admitted against the
# complete pipeline envelope
# (`ipw_point_if_general_state_reserve`), it is built once from one
# `query_ipw_view_chunk(read_outcomes = FALSE, include_design = TRUE)` per chunk
# and link (`ipw_point_if_general_state_build`), sealed, then sliced by the two
# covariance passes and the two GMM projection passes (`state =` on
# `stream_ipw_point_covariance_general`, `replay_ipw_point_if_general` and
# `ipw_point_if_eval_chunk_internal`), which still read every raw endpoint and
# still verify plan and state at every pass boundary and after every callback;
# `ipw_point_if_general_state_release` drops it after the solve.
#
# Expected values come from (a) the dense `materialize_ipw_if` and
# `stream_ipw_covariance` references and the four reported IF components on the
# accepted weighted-96 fixture (test-chained-ipw-point-view.R), (b) the same
# build with `state = NULL`, the standalone `stream_ipw_point_covariance_view`
# and a forced-NULL or typed-capacity-declined pipeline (the old general route)
# on the accepted CP08 panels (test-cp08-ipw-general-norm-reuse.R), (c) the
# fixture's own observation pattern and (d) the admission formulas evaluated
# from pre-existing plan and preflight fields. The state never supplies its own
# truth. Tolerances are the originals: bitwise where the arithmetic path is
# identical, 1e-12 for whole-structure comparisons; identities, orders, counts
# and refusal classes compare exactly.
#
# STATUS: UNEXECUTED at authoring time (base-R parse/AST and static
# signature/arithmetic checks only; no build, install, fit, RNG or test run).
# The re-bound namespace copies are disclosed test-local clones used only to
# count and record calls; every value comes from genuine installed functions.
# The root sequences the installed gate.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api), environmentIsLocked(api))
if (!exists("ipw_point_if_general_state_reserve", envir = api, inherits = FALSE)) {
  stop("The installed cdid lacks the D90 general replay state; install the candidate build.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Accepted helpers: test-cp08-ipw-general-norm-reuse.R (gn_*) and
# test-chained-ipw-point-view.R (cp02_pv_*). Only the simultaneous-entry
# comparison is adapted for the final pipeline's additional live state charge.
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
  state_entries <- function(r) {
    if (isTRUE(r$resource_accounting$point_general$replay_state$admitted)) d90_result_need(r)$entries else 0
  }
  expect_identical(a$resource_accounting$simultaneous_entries - state_entries(a),
    b$resource_accounting$simultaneous_entries - state_entries(b), info = info)
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches, info = info)
  expect_identical(a$gmm_fit$if_projection$max_batch_rows_seen, b$gmm_fit$if_projection$max_batch_rows_seen, info = info)
  expect_identical(names(a$point_links_used), names(b$point_links_used), info = info)
  invisible(TRUE)
}

cp02_pv_fn <- function(name) {
  f <- get(name, envir = api, inherits = FALSE)
  if (!is.function(f) || !identical(environment(f), api)) {
    stop("Foreign function binding: ", name, call. = FALSE)
  }
  f
}

cp02_pv_weighted96_fixture <- function(weighted = TRUE) {
  n <- 96L
  G <- rep(c(4, 7, Inf, 9), each = 24L)
  x <- rep(seq(-1.4, 1.4, length.out = 24L), 4L)
  Z <- cbind(one = 1, x = x)
  w <- if (weighted) 0.7 + ((seq_len(n) * 13L) %% 29L) / 20 else rep(1, n)
  obs <- matrix(FALSE, n, 8L)
  for (i in seq_len(n)) {
    k <- ((i - 1L) %% 24L) + 1L
    if (G[i] != 9 && k <= 18L) {
      obs[i, 1:6] <- TRUE
    } else {
      obs[i, 7:8] <- TRUE
      obs[i, 1L + (k %% 6L)] <- TRUE
    }
  }
  Y <- outer(seq_len(n) / 17, 1:8, function(a, b) sin(a * b) + a + b^2 / 13) +
    outer(x, 1:8, function(a, b) a * b / 3)
  pan <- data.frame(
    id = rep(seq_len(n), each = 8L),
    time = rep(as.double(1:8), n),
    y = as.vector(t(ifelse(obs, Y, NA_real_))),
    g = rep(G, each = 8L)
  )
  idx <- cp02_pv_fn("construct_observed_pairs")(pan, "id", "time", "y", "g")
  ctx <- cp02_pv_fn("construct_ipw_context")(idx, id = seq_len(n), Z = Z, weights = w)
  list(n = n, G = G, Z = Z, w = w, obs = obs, Y = Y, panel = pan, index = idx, context = ctx)
}

# ---------------------------------------------------------------------------
# D90 helpers
# ---------------------------------------------------------------------------

# Final reviewed inventory, evaluated from pre-existing plan/preflight facts
# (never read back from a reservation or state). E = 6 S + (C + 1) L;
# reserve 2 E + 320 L + 2048 + 16 B. Work includes one producer construction,
# (C + 16) E + 7 C L^2 copies/slices, all boundary checks, 4 C chunk guards,
# and setup with two additional full plan verifications.
d90_need <- function(S, C, L, B, pre_work, plan_work) {
  S <- as.double(S); C <- as.double(C); L <- as.double(L); B <- as.double(B)
  E <- 6 * S + (C + 1) * L
  checks <- 1 + 4 * (C + 1)
  one_check <- 64 * L + 256
  copy_work <- (C + 16) * E + 7 * C * L^2
  check_work <- checks * one_check
  guard_work <- 4 * C * one_check
  setup_work <- 2 * as.double(plan_work) + one_check
  list(E = E, entries = 2 * E + 320 * L + 2048 + 16 * B, checks = checks,
    copy_work = copy_work, one_check = one_check, check_work = check_work,
    guard_work = guard_work, setup_work = setup_work,
    work = as.double(pre_work) + copy_work + check_work + guard_work + setup_work)
}
d90_plan_need <- function(plan) d90_need(
  sum(vapply(plan$point_links, function(lk) as.double(lk$counts$n_observed_risk_units), numeric(1L))),
  plan$preflight$n_chunks, plan$n_links, plan$step, plan$preflight$total_work_ops_bound,
  plan$general$one_verification_work_ops_bound)
d90_result_need <- function(r) {
  pp <- r$resource_accounting$point_producer
  d90_need(sum(vapply(r$point_links_used, function(lk) as.double(lk$counts$n_observed_risk_units), numeric(1L))),
    pp$n_chunks, length(r$point_links_used), pp$effective_chunk_rows, pp$total_work_ops_bound,
    r$resource_accounting$point_general$one_verification_work_ops_bound)
}

# Recorders on a re-bound copy: the outcome flag of every chunk query, the
# positions of every raw endpoint read (in order) and the last built state.
d90_spy <- function(e) {
  log <- new.env(parent = emptyenv())
  log$query_outcomes <- logical(0L)
  log$reads <- list()
  log$state <- NULL
  query <- get("query_ipw_view_chunk", envir = e, inherits = FALSE)
  read <- get("ipw_owner_read_y_positions", envir = e, inherits = FALSE)
  build <- get("ipw_point_if_general_state_build", envir = e, inherits = FALSE)
  assign("query_ipw_view_chunk", function(..., read_outcomes = FALSE) {
    log$query_outcomes <- c(log$query_outcomes, isTRUE(read_outcomes))
    query(..., read_outcomes = read_outcomes)
  }, envir = e)
  assign("ipw_owner_read_y_positions", function(panel_owner, pos_s, pos_t, trace_env = NULL) {
    log$reads[[length(log$reads) + 1L]] <- list(pos_s = pos_s, pos_t = pos_t)
    read(panel_owner = panel_owner, pos_s = pos_s, pos_t = pos_t, trace_env = trace_env)
  }, envir = e)
  assign("ipw_point_if_general_state_build", function(...) {
    log$state <- build(...)
    log$state
  }, envir = e)
  log
}
d90_spy_reset <- function(log) {
  log$query_outcomes <- logical(0L)
  log$reads <- list()
}
d90_counts <- c("ipw_point_if_general_plan", "ipw_point_if_general_verify", "ipw_point_if_general_state_reserve",
  "ipw_point_if_general_state_build", "ipw_point_if_general_state_verify", "ipw_point_if_general_state_link_slice",
  "ipw_point_if_general_state_release", "replay_ipw_point_if_general", "stream_ipw_point_covariance_general",
  "stream_ipw_point_covariance_view")
# Isolate the sparse-state fallback in this pipeline comparison only. Both
# clones decline raw IF retention. With that reservation bypassed, its plan
# check has not happened: discard the pipeline's .verified_entry hint so the
# sparse reservation performs its own original check. The state-free clone's
# forced-NULL sparse reservation still performs no check, as before.
d90_new <- function(edits = list(), general_edit = gn_general_edit) {
  e <- gn_ns(c(general_edit, list(ipw_point_if_retained_reserve = function(...) NULL), edits), d90_counts)
  reserve <- e$ipw_point_if_general_state_reserve
  e$ipw_point_if_general_state_reserve <- function(..., .verified_entry = FALSE) reserve(...)
  e$.log <- d90_spy(e)
  e
}
d90_old <- function() d90_new(list(ipw_point_if_general_state_reserve = function(...) NULL))
d90_has_kind <- function(x, kinds = c("cdid_ipw_point_if_general_plan", "cdid_ipw_point_if_general_state",
                                      "cdid_ipw_point_if_general_link_state")) {
  found <- FALSE
  walk <- function(y, depth = 0L) {
    if (found || depth > 10L) return(invisible(NULL))
    if (is.list(y)) {
      k <- y[["kind", exact = TRUE]]
      if (is.character(k) && length(k) == 1L && k %in% kinds) found <<- TRUE
      for (v in y) walk(v, depth + 1L)
    }
    invisible(NULL)
  }
  walk(x)
  found
}

# The accepted weighted-96 point-link set of test-chained-ipw-point-view.R
# (bounded 4-term replay test): off-risk endpoints corrupted to 1e308, three
# links sharing one fitted nuisance and one distinct, the dense shared fits and
# the dense IF reference on the uncorrupted context; plus the fixture's own
# risk set and observed-risk sets per link.
d90_w96 <- function() {
  construct_ipw_panel_owner <- cp02_pv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_pv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_pv_fn("construct_ipw_nuisance_view")
  construct_ipw_link_view <- cp02_pv_fn("construct_ipw_link_view")
  fit_ipw_nuisance_view <- cp02_pv_fn("fit_ipw_nuisance_view")
  fit_ipw_point_link_view_set <- cp02_pv_fn("fit_ipw_point_link_view_set")
  ipw_owner_control <- cp02_pv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_pv_fn("ipw_fitted_view_control")
  ipw_point_view_control <- cp02_pv_fn("ipw_point_view_control")
  construct_ipw_shared_session <- cp02_pv_fn("construct_ipw_shared_session")
  construct_ipw_nuisance <- cp02_pv_fn("construct_ipw_nuisance")
  declare_ipw_shared_link <- cp02_pv_fn("declare_ipw_shared_link")
  fit_ipw_link_shared <- cp02_pv_fn("fit_ipw_link_shared")
  materialize_ipw_if <- cp02_pv_fn("materialize_ipw_if")

  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)
  fx96 <- cp02_pv_weighted96_fixture(weighted = TRUE)

  # Off-risk corruption test: corrupt endpoints for G == 9 (outside risk set g=4, K=c(7, Inf))
  pan_corrupt <- fx96$panel
  pan_corrupt$y[pan_corrupt$g == 9 & pan_corrupt$time %in% c(1, 2, 3)] <- 1e308

  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
  po <- construct_ipw_panel_owner(pan_corrupt, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)

  nv12 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  nv78 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 7, t = 8, control = own_ctrl)
  vfit12 <- fit_ipw_nuisance_view(nv12, glm_control = glm_ctrl, control = fv_ctrl)
  vfit78 <- fit_ipw_nuisance_view(nv78, glm_control = glm_ctrl, control = fv_ctrl)

  lv12 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, nuisance_view = nv12, control = own_ctrl)
  lv23 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 2, t = 3, nuisance_view = nv12, control = own_ctrl)
  lv13 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 1, t = 3, nuisance_view = nv12, control = own_ctrl)
  lv78 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 7, t = 8, nuisance_view = nv78, control = own_ctrl)

  pv_ctrl <- ipw_point_view_control(16, 96, 100000, 5000000, 1000000000000, 100000)
  plinks <- fit_ipw_point_link_view_set(
    list(vfit12, vfit12, vfit12, vfit78),
    list(lv12, lv23, lv13, lv78),
    control = pv_ctrl
  )

  sess <- construct_ipw_shared_session(fx96$context, max_nuisance = 8L, max_links = 16L, max_unit_entries = 10000L)
  nu12 <- construct_ipw_nuisance(sess, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl)
  nu78 <- construct_ipw_nuisance(sess, g = 4, K = c(7, Inf), s = 7, t = 8, control = glm_ctrl)
  d12 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl))
  d23 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 2, t = 3, control = glm_ctrl))
  d13 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 1, t = 3, control = glm_ctrl))
  d78 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu78, g = 4, K = c(7, Inf), s = 7, t = 8, control = glm_ctrl))
  dlinks <- list(d12, d23, d13, d78)
  dense_mat <- materialize_ipw_if(fx96$context, dlinks, max_entries = 10000, chunk_size = 16L, include_components = TRUE)

  pairs <- list(c(1, 2), c(2, 3), c(1, 3), c(7, 8))
  in_risk <- fx96$G %in% c(4, 7, Inf)
  in_obs <- lapply(pairs, function(pr) in_risk & fx96$obs[, pr[1L]] & fx96$obs[, pr[2L]])
  list(fx96 = fx96, po = po, fo = fo, plinks = plinks, dlinks = dlinks, dense_mat = dense_mat,
    pairs = pairs, in_risk = in_risk, in_obs = in_obs)
}
d90_pc <- function(chunk) cp02_pv_fn("ipw_point_view_control")(chunk, 96, 100000, 5000000, 1000000000000, 100000)

# ---------------------------------------------------------------------------
# 1. Pipeline: the admitted state equals the forced-NULL and the typed-
#    capacity-declined general route through covariance, GMM and full-frame
#    inference; one outcome-free setup replaces four query passes, every raw
#    endpoint read and plan check is unchanged, and the state is released.
# ---------------------------------------------------------------------------

test_that("admitted replay state equals the old general route through covariance, GMM and inference; one outcome-free setup, same reads", {
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
    e_new <- d90_new(); e_old <- d90_old()
    r_new <- tryCatch(gn_run(e_new, d, ctrl, cs$p), error = identity)
    r_old <- tryCatch(gn_run(e_old, d, ctrl, cs$p), error = identity)
    ok <- gn_expect_same_result(r_new, r_old, info = info)
    if (!isTRUE(ok)) next
    n_fitted <- n_fitted + 1L
    expect_null(r_new$resource_accounting$point_batch, info = info)
    pp <- r_new$resource_accounting$point_producer
    L <- length(r_new$point_links_used)
    C <- as.integer(pp$n_chunks)
    S <- sum(vapply(r_new$point_links_used, function(lk) as.double(lk$counts$n_observed_risk_units), numeric(1L)))
    expect_identical(2 * S, as.double(pp$total_expected_endpoints_read), info = info)
    need <- d90_result_need(r_new)
    # Admitted: built once, served the four passes, checked at the covariance
    # entry, every pass boundary and after every callback, then released.
    pg_new <- r_new$resource_accounting$point_general
    pg_old <- r_old$resource_accounting$point_general
    rs <- pg_new$replay_state
    expect_true(isTRUE(rs$admitted), info = paste(info, "(state declined under the fixture's public caps)"))
    expect_identical(rs$status, "retained_and_released", info = info)
    expect_identical(c(rs$passes_served, rs$chunks_served, rs$state_verifications),
      c(4L, 4L * C, 1L + 4L * (C + 1L)), info = info)
    expect_identical(c(rs$queries_performed, rs$endpoints_read_during_construction), c(C * L, 0L), info = info)
    expect_identical(rs$total_observed_risk_units, S, info = info)
    expect_identical(c(rs$raw_state_entries, rs$reserved_entries, rs$work_ops_bound, rs$state_checks_planned),
      c(need$E, need$entries, need$work, need$checks), info = info)
    expect_identical(c(rs$copy_slice_work_ops_bound, rs$state_check_work_ops_bound,
      rs$chunk_guard_work_ops_bound, rs$setup_work_ops_bound),
      c(need$copy_work, need$check_work, need$guard_work, need$setup_work), info = info)
    for (nm in c("outcomes_retained", "differences_retained", "full_n_by_L_retained", "per_link_dense_masks_retained",
                 "peak_rss_promised")) expect_false(isTRUE(rs[[nm]]), info = paste(info, nm))
    # Forced NULL: the declined record states the same requirement; the rest of
    # the general record is unchanged; the state's work is charged once.
    ro <- pg_old$replay_state
    expect_identical(ro$status, "declined_capacity", info = info)
    expect_false(isTRUE(ro$admitted), info = info)
    expect_identical(c(ro$reserved_entries, ro$work_ops_bound), c(need$entries, need$work), info = info)
    expect_identical(ro$queries_performed, 0L, info = info)
    keep <- setdiff(names(pg_old), "replay_state")
    expect_identical(pg_new[keep], pg_old[keep], info = info)
    expect_identical(r_new$resource_accounting$downstream_work_ops - r_old$resource_accounting$downstream_work_ops,
      need$work, info = info)
    expect_identical(r_new$resource_accounting$simultaneous_entries - r_old$resource_accounting$simultaneous_entries,
      need$entries, info = info)
    # Spies: C*L outcome-free queries replace 4*C*L outcome queries; every raw
    # endpoint read of the whole fit is the same, in the same order; the plan
    # boundary checks remain; reservation and build add two plan checks.
    qn <- e_new$.log$query_outcomes; qo <- e_old$.log$query_outcomes
    expect_identical(sum(!qn) - sum(!qo), C * L, info = info)
    expect_identical(sum(qo) - sum(qn), 4L * C * L, info = info)
    expect_identical(e_new$.log$reads, e_old$.log$reads, info = paste(info, "endpoint-read trace"))
    for (nm in c("ipw_point_if_general_plan", "ipw_point_if_general_state_reserve",
                 "ipw_point_if_general_state_release", "stream_ipw_point_covariance_general",
                 "replay_ipw_point_if_general", "stream_ipw_point_covariance_view")) {
      expect_identical(gn_n(e_new, nm), gn_n(e_old, nm), info = paste(info, nm))
    }
    expect_identical(gn_n(e_new, "ipw_point_if_general_verify"), 3L + 4L * (1L + C), info = info)
    expect_identical(gn_n(e_old, "ipw_point_if_general_verify"), 1L + 4L * (1L + C), info = info)
    expect_identical(gn_n(e_new, "replay_ipw_point_if_general"), 4L, info = info)
    expect_identical(gn_n(e_new, "ipw_point_if_general_state_build"), 1L, info = info)
    expect_identical(gn_n(e_new, "ipw_point_if_general_state_verify"), 1L + 4L * (1L + C), info = info)
    expect_identical(gn_n(e_new, "ipw_point_if_general_state_link_slice"), 4L * C * L, info = info)
    for (nm in c("ipw_point_if_general_state_build", "ipw_point_if_general_state_verify",
                 "ipw_point_if_general_state_link_slice")) expect_identical(gn_n(e_old, nm), 0L, info = paste(info, nm))
    # Released after the solve; no plan or state reference survives the call.
    st <- e_new$.log$state
    expect_true(is.list(st), info = paste(info, "state built"))
    if (is.list(st)) {
      expect_identical(st$runtime$status, "released", info = info)
      expect_null(get0("point_if_general_state", envir = st$ref, inherits = FALSE), info = info)
      for (lref in st$link_state_refs) expect_null(get0("point_if_general_link_state", envir = lref, inherits = FALSE), info = info)
    }
    expect_false(d90_has_kind(r_new), info = info)
  }
  expect_gte(n_fitted, 4L)

  d <- gn_panel(pattern = "half"); ctrl <- gn_ctrl(d, 1L, chunk = 7L)
  e_ref <- d90_old()
  r_ref <- tryCatch(gn_run(e_ref, d, ctrl, 1L), error = identity)
  # Typed capacity decline through the genuine reservation: one work operation
  # short of the state's requirement declines before any state query.
  e_dec <- d90_new()
  reserve <- get("ipw_point_if_general_state_reserve", envir = e_dec, inherits = FALSE)
  declined <- logical(0L)
  assign("ipw_point_if_general_state_reserve", function(plan, control, live = 0, work = 0, envelope_control = NULL,
                                                       .verified_entry = FALSE) {
    ec <- envelope_control
    ec$max_work_ops <- work + d90_plan_need(plan)$work - 1
    out <- reserve(plan, control, live = live, work = work, envelope_control = ec)
    declined <<- c(declined, is.null(out))
    out
  }, envir = e_dec)
  r_dec <- tryCatch(gn_run(e_dec, d, ctrl, 1L), error = identity)
  gn_expect_same_result(r_dec, r_ref, info = "typed capacity decline vs forced NULL, half")
  expect_identical(declined, TRUE)
  if (!inherits(r_dec, "error")) {
    expect_identical(r_dec$resource_accounting$point_general$replay_state$status, "declined_capacity")
    expect_identical(r_dec$resource_accounting$downstream_work_ops, r_ref$resource_accounting$downstream_work_ops)
  }
  expect_identical(e_dec$.log$query_outcomes, e_ref$.log$query_outcomes)
  expect_identical(e_dec$.log$reads, e_ref$.log$reads)
  for (nm in c("ipw_point_if_general_state_build", "ipw_point_if_general_state_verify",
               "ipw_point_if_general_state_link_slice")) expect_identical(gn_n(e_dec, nm), 0L, info = nm)
  # Installed functions with natural batch eligibility, cloned only to isolate
  # the sparse fallback and record its checks, admit the state on missing data.
  e_api <- d90_new(general_edit = list())
  r_api <- tryCatch(gn_run(e_api, d, ctrl, 1L), error = identity)
  gn_expect_same_result(r_api, r_ref, info = "natural eligibility, raw-NULL clone vs state-free copy, half")
  if (!inherits(r_api, "error")) {
    expect_null(r_api$resource_accounting$point_batch)
    expect_true(isTRUE(r_api$resource_accounting$point_general$replay_state$admitted))
  }
  # A failure after construction still releases the state.
  e_err <- d90_new(list(solve_chained_gmm_core = function(...) stop("D90 injected GMM solve failure", call. = FALSE)))
  r_err <- tryCatch(gn_run(e_err, d, ctrl, 1L), error = identity)
  expect_s3_class(r_err, "error")
  st_err <- e_err$.log$state
  expect_true(is.list(st_err))
  if (is.list(st_err)) {
    expect_identical(st_err$runtime$status, "released")
    expect_null(get0("point_if_general_state", envir = st_err$ref, inherits = FALSE))
  }
})

# ---------------------------------------------------------------------------
# 2. Weighted-96: state slices reproduce the dense IF, the four components and
#    the centered covariance bitwise across chunk sizes (single, tail, empty
#    observer and off-risk chunks), permuted links and cache modes.
# ---------------------------------------------------------------------------

test_that("weighted-96 state replay, component chunks and centered covariance equal the dense references and the state-free replay bitwise", {
  w <- d90_w96()
  dense_mat <- w$dense_mat
  nonobs <- which(w$in_risk & !w$in_obs[[1L]])
  outside <- which(!w$in_risk)
  expect_gt(length(nonobs), 0L)
  dense_replay_num_ok <- function(actual, expected, tol = 1e-12) {
    identical(typeof(actual), "double") &&
      identical(typeof(expected), "double") &&
      identical(is.object(actual), is.object(expected)) &&
      identical(attributes(actual), attributes(expected)) &&
      identical(length(actual), length(expected)) &&
      length(actual) > 0L &&
      all(is.finite(actual)) &&
      all(is.finite(expected)) &&
      all(actual[expected == 0] == 0) &&
      all(abs(actual - expected) <= tol * pmax(1, abs(actual), abs(expected)))
  }
  expect_dense_replay_num <- function(actual, expected, tol = 1e-12, info = NULL) {
    expect_true(dense_replay_num_ok(actual, expected, tol = tol), info = info)
  }
  expect_false(dense_replay_num_ok(numeric(), numeric()))
  expect_false(dense_replay_num_ok(matrix(numeric(), 0L, 4L), matrix(numeric(), 0L, 4L)))
  expect_false(dense_replay_num_ok(NA_real_, NA_real_))
  expect_false(dense_replay_num_ok(5e-13, 0))
  expect_false(dense_replay_num_ok(
    structure(as.integer(dense_mat$phi_matrix), dim = dim(dense_mat$phi_matrix), dimnames = dimnames(dense_mat$phi_matrix)),
    dense_mat$phi_matrix
  ))
  expect_false(dense_replay_num_ok(as.vector(dense_mat$phi_matrix), dense_mat$phi_matrix))
  expect_false(dense_replay_num_ok(
    dense_mat$phi_matrix + 2e-12 * pmax(1, abs(dense_mat$phi_matrix)),
    dense_mat$phi_matrix
  ))
  expect_false(dense_replay_num_ok(-dense_mat$phi_matrix, dense_mat$phi_matrix))
  expect_false(dense_replay_num_ok(
    dense_mat$phi_matrix[c(2L, 1L, 3:nrow(dense_mat$phi_matrix)), , drop = FALSE],
    dense_mat$phi_matrix
  ))
  fields <- c("a_index", "cohort_index", "obs_s_index", "obs_t_index", "odds", "score_base", "chunk_offsets")
  for (ch_sz in c(1L, 5L, 8L, 13L, 96L)) {
    info <- paste("chunk", ch_sz)
    pc <- d90_pc(ch_sz)
    plan <- api$ipw_point_if_general_plan(w$plinks, pc)
    C <- as.integer(plan$preflight$n_chunks); L <- plan$n_links; B <- plan$step
    expect_identical(C, as.integer(ceiling(96 / ch_sz)), info = info)
    st <- api$ipw_point_if_general_state(plan, pc)
    expect_true(is.list(st) && identical(st$kind, "cdid_ipw_point_if_general_state"), info = paste(info, "admitted"))
    if (!is.list(st)) next
    # Sparse layout against the observation pattern: observed risk A-units in A
    # order, per-chunk counts (empty observer chunks are zero steps), treated first.
    for (ell in seq_len(L)) {
      lst <- st$link_states[[ell]]
      a_obs <- which(w$in_obs[[ell]])
      expect_identical(lst$a_index, a_obs, info = paste(info, "link", ell))
      expect_identical(diff(lst$chunk_offsets), tabulate(ceiling(a_obs / B), nbins = C), info = paste(info, "link", ell))
      expect_identical(lst$cohort_index == 1L, w$fx96$G[a_obs] == 4, info = paste(info, "link", ell))
    }
    if (C > 1L) expect_true(any(diff(st$link_states[[4L]]$chunk_offsets) == 0L), info = paste(info, "empty observer chunk"))
    # Replays: the state equals state = NULL (batches and read trace) and the dense IF.
    tr_s <- new.env(parent = emptyenv()); tr_o <- new.env(parent = emptyenv())
    gs <- gn_collect(api$replay_ipw_point_if_general, plan = plan, control = pc, trace_env = tr_s, state = st)
    go <- gn_collect(api$replay_ipw_point_if_general, plan = plan, control = pc, trace_env = tr_o)
    expect_identical(gs, go, info = info)
    expect_identical(length(gs), C, info = info)
    for (i in seq_along(gs)) {
      expect_dense_replay_num(gs[[i]]$phi, dense_mat$phi_matrix[gs[[i]]$a_index, , drop = FALSE], info = paste(info, "batch", i))
    }
    for (nm in c("preflight_completed", "batches_emitted", "units_emitted", "active_query_chunks", "active_nuisance_caches",
                 "max_simultaneous_nuisance_caches", "nuisance_cache_hits", "endpoints_read", "last_chunk_start",
                 "last_chunk_end")) expect_identical(tr_s[[nm]], tr_o[[nm]], info = paste(info, "trace", nm))
    expect_identical(c(tr_s$max_simultaneous_query_chunks, tr_o$max_simultaneous_query_chunks), c(0L, 1L), info = info)
    expect_identical(as.double(tr_s$endpoints_read), as.double(plan$preflight$total_expected_endpoints_read), info = info)
    expect_gt(tr_s$nuisance_cache_hits, 0L)
    # Component-enabled chunks of the sealed partition: state = NULL and dense.
    phi_all <- matrix(0, nrow = 0L, ncol = 4L)
    comp_all <- list(rate = phi_all, beta_score = phi_all, treated = phi_all, control = phi_all)
    for (k in seq_len(C)) {
      a_s <- 1L + (k - 1L) * B
      a_e <- min(96L, k * B)
      bs <- api$ipw_point_if_eval_chunk_internal(plan$norm, a_s, a_e, plan$control, include_components = TRUE,
        batch_index = k, state = st)
      bo <- api$ipw_point_if_eval_chunk_internal(plan$norm, a_s, a_e, plan$control, include_components = TRUE,
        batch_index = k)
      expect_identical(bs[setdiff(names(bs), "work")], bo[setdiff(names(bo), "work")], info = paste(info, "chunk", k))
      expect_identical(bs$work$max_simultaneous_query_chunks, 0L, info = paste(info, "chunk", k))
      idx <- bs$a_index
      expect_dense_replay_num(bs$phi, dense_mat$phi_matrix[idx, , drop = FALSE], info = paste(info, "chunk", k))
      for (nm in names(comp_all)) {
        if (nm %in% c("beta_score", "control")) {
          expect_dense_replay_num(bs$components[[nm]], dense_mat$components[[nm]][idx, , drop = FALSE], info = paste(info, "chunk", k, nm))
        } else {
          expect_identical(bs$components[[nm]], dense_mat$components[[nm]][idx, , drop = FALSE], info = paste(info, "chunk", k, nm))
        }
        comp_all[[nm]] <- rbind(comp_all[[nm]], bs$components[[nm]])
      }
      phi_all <- rbind(phi_all, bs$phi)
    }
    expect_dense_replay_num(phi_all, dense_mat$phi_matrix, info = info)
    expect_identical(comp_all$rate, dense_mat$components$rate, info = info)
    expect_dense_replay_num(comp_all$beta_score, dense_mat$components$beta_score, info = info)
    expect_identical(comp_all$treated, dense_mat$components$treated, info = info)
    expect_dense_replay_num(comp_all$control, dense_mat$components$control, info = info)
    # Nonobserver risk units carry exactly the s0 rate term; off-risk rows are zero;
    # the four signed terms sum in the original order and each contributes.
    expect_true(all(abs(phi_all[nonobs, 1L]) > 1e-6), info = info)
    expect_identical(phi_all[nonobs, 1L], comp_all$rate[nonobs, 1L], info = info)
    for (nm in c("treated", "control", "beta_score")) expect_true(all(comp_all[[nm]][nonobs, 1L] == 0), info = paste(info, nm))
    expect_true(all(phi_all[outside, ] == 0), info = info)
    expect_identical(phi_all, ((comp_all$treated + comp_all$control) + comp_all$beta_score) + comp_all$rate, info = info)
    for (nm in names(comp_all)) expect_true(any(comp_all[[nm]] != 0), info = paste(info, nm, "contributes"))
    # Centered two-pass covariance: state, state = NULL, standalone view, dense.
    cg_s <- api$stream_ipw_point_covariance_general(plan, pc, state = st)
    cg_o <- api$stream_ipw_point_covariance_general(plan, pc)
    cv <- api$stream_ipw_point_covariance_view(w$plinks, control = pc)
    dense_cov <- cp02_pv_fn("stream_ipw_covariance")(w$fx96$context, w$dlinks, chunk_size = ch_sz)
    for (nm in c("covariance", "mean_phi", "estimates", "n_a_units")) {
      expect_identical(cg_s[[nm]], cg_o[[nm]], info = paste(info, nm))
      expect_identical(cg_s[[nm]], cv[[nm]], info = paste(info, nm))
    }
    for (nm in c("covariance", "mean_phi", "estimates")) expect_identical(cg_s[[nm]], dense_cov[[nm]], info = paste(info, "dense", nm))
    expect_identical(c(cg_s$replay$served_from, cg_o$replay$served_from),
      c("general_sparse_replay_state", "general_chunk_query"), info = info)
    expect_identical(cg_s$replay$state_verifications, 1L + C, info = info)
    contrast <- c(1, 1, -1, 0)
    expect_lt(abs(as.double(t(contrast) %*% cg_s$covariance %*% contrast)), 1e-14)
    if (ch_sz %in% c(5L, 13L)) {
      # Permuted links (the distinct nuisance splits the shared ones): the same
      # sealed per-link vectors, the state-free batches, permuted columns.
      perm <- c(2L, 4L, 3L, 1L)
      plan_p <- api$ipw_point_if_general_plan(w$plinks[perm], pc)
      st_p <- api$ipw_point_if_general_state(plan_p, pc)
      for (ell in seq_along(perm)) {
        expect_identical(st_p$link_states[[ell]][fields], st$link_states[[perm[ell]]][fields], info = paste(info, "permuted link", ell))
      }
      gp <- gn_collect(api$replay_ipw_point_if_general, plan = plan_p, control = pc, state = st_p)
      expect_identical(gp, gn_collect(api$replay_ipw_point_if_general, plan = plan_p, control = pc), info = info)
      for (i in seq_along(gp)) {
        expect_identical(unname(gp[[i]]$phi), unname(gs[[i]]$phi[, perm, drop = FALSE]), info = paste(info, "permuted batch", i))
      }
      cg_p <- api$stream_ipw_point_covariance_general(plan_p, pc, state = st_p)
      expect_equal(unname(cg_p$covariance), unname(cg_s$covariance[perm, perm]), tolerance = 1e-12, info = info)
      # Without the one-nuisance chunk cache: the same sealed vectors and batches.
      st_nc <- api$ipw_point_if_general_state(plan, pc, use_nuisance_cache = FALSE)
      expect_identical(lapply(st_nc$link_states, `[`, fields), lapply(st$link_states, `[`, fields), info = info)
      expect_identical(gn_collect(api$replay_ipw_point_if_general, plan = plan, control = pc, use_nuisance_cache = FALSE,
        state = st_nc), gs, info = info)
    }
  }
})

# ---------------------------------------------------------------------------
# 3. Spy: one outcome-free query/score pass serves the four replays; every
#    endpoint read (risk-first), plan check and state check still happens.
# ---------------------------------------------------------------------------

test_that("one outcome-free query and score setup serves four replays while every endpoint read, plan check and state check still happens", {
  w <- d90_w96()
  e <- gn_ns(list(), c("ipw_point_if_general_verify", "ipw_point_if_general_state_verify",
    "ipw_point_if_general_state_link_slice", "ipw_point_if_general_state_build"))
  log <- d90_spy(e)
  pc <- d90_pc(13L)
  plan <- e$ipw_point_if_general_plan(w$plinks, pc)
  C <- as.integer(plan$preflight$n_chunks); L <- plan$n_links
  expect_identical(c(C, L), c(8L, 4L))
  # Construction: one outcome-free query per chunk and link, no endpoint read.
  gn_reset(e); d90_spy_reset(log)
  st <- e$ipw_point_if_general_state(plan, pc)
  expect_identical(gn_n(e, "ipw_point_if_general_state_build"), 1L)
  expect_identical(gn_n(e, "ipw_point_if_general_verify"), 2L)
  expect_identical(log$query_outcomes, rep(FALSE, C * L))
  expect_length(log$reads, 0L)
  four <- function(state) list(
    cov = e$stream_ipw_point_covariance_general(plan, pc, state = state)[c("covariance", "mean_phi", "estimates", "n_a_units")],
    p3 = gn_collect(e$replay_ipw_point_if_general, plan = plan, control = pc, state = state),
    p4 = gn_collect(e$replay_ipw_point_if_general, plan = plan, control = pc, state = state))
  checked <- c("ipw_point_if_general_verify", "ipw_point_if_general_state_verify", "ipw_point_if_general_state_link_slice",
    "ipw_point_if_general_state_build")
  gn_reset(e); d90_spy_reset(log)
  out_s <- four(st)
  q_s <- log$query_outcomes; reads_s <- log$reads
  n_s <- vapply(checked, function(nm) gn_n(e, nm), integer(1L), USE.NAMES = FALSE)
  gn_reset(e); d90_spy_reset(log)
  out_o <- four(NULL)
  n_o <- vapply(checked, function(nm) gn_n(e, nm), integer(1L), USE.NAMES = FALSE)
  expect_identical(out_s, out_o)
  expect_length(q_s, 0L)
  expect_identical(log$query_outcomes, rep(TRUE, 4L * C * L))
  expect_length(reads_s, 4L * C * L)
  expect_identical(reads_s, log$reads)
  expect_identical(n_s, c(1L + 4L * (C + 1L), 1L + 4L * (C + 1L), 4L * C * L, 0L))
  expect_identical(n_o, c(1L + 4L * (C + 1L), 0L, 0L, 0L))
  expect_identical(c(st$runtime$passes_completed, st$runtime$chunks_evaluated, st$runtime$verifications),
    c(4L, 4L * C, 1L + 4L * (C + 1L)))
  # Risk-first: per link, one pass reads exactly the observed risk units at the
  # link's two dates, never the corrupted off-risk endpoints.
  obs <- plan$panel_owner$observations
  pass1 <- reads_s[seq_len(C * L)]
  for (ell in seq_len(L)) {
    calls <- pass1[seq.int(ell, C * L, by = L)]
    ps <- as.integer(unlist(lapply(calls, `[[`, "pos_s")))
    pt <- as.integer(unlist(lapply(calls, `[[`, "pos_t")))
    units <- as.integer(plan$panel_owner$a_unit_index[which(w$in_obs[[ell]])])
    expect_identical(as.integer(obs$unit_index[ps]), units, info = paste("link", ell))
    expect_identical(as.integer(obs$unit_index[pt]), units, info = paste("link", ell))
    expect_true(all(obs$time[ps] == w$pairs[[ell]][1L]) && all(obs$time[pt] == w$pairs[[ell]][2L]), info = paste("link", ell))
    expect_false(any(obs$y[c(ps, pt)] == 1e308), info = paste("link", ell))
  }
  # Release: scalar accounting only, then every sealed binding is dropped.
  rec <- e$ipw_point_if_general_state_release(st, plan)
  expect_identical(rec$status, "retained_and_released")
  expect_true(isTRUE(rec$admitted) && isTRUE(rec$released))
  expect_identical(c(rec$passes_served, rec$chunks_served, rec$state_verifications, rec$queries_performed,
    rec$endpoints_read_during_construction), c(4L, 4L * C, 1L + 4L * (C + 1L), C * L, 0L))
  expect_identical(st$runtime$status, "released")
  expect_null(get0("point_if_general_state", envir = st$ref, inherits = FALSE))
  expect_null(get0("plan_snapshot", envir = st$binding_ref, inherits = FALSE))
  for (lref in st$link_state_refs) expect_null(get0("point_if_general_link_state", envir = lref, inherits = FALSE))
})

# ---------------------------------------------------------------------------
# 4. Refusals: callback mutations (first and last callback), unlocked or
#    rebound bindings, foreign plans/states, a shifted partition, COW-mutated
#    payloads and released states. Each refuses before the
#    next endpoint read or before returning.
# ---------------------------------------------------------------------------

test_that("state checks refuse mutated, foreign, unlocked, rebound, shifted and released states, including a mutation in the last callback", {
  w <- d90_w96()
  e <- gn_ns(list(), character())
  log <- d90_spy(e)
  pc <- d90_pc(13L)
  plan <- e$ipw_point_if_general_plan(w$plinks, pc)
  C <- as.integer(plan$preflight$n_chunks); L <- plan$n_links
  st <- e$ipw_point_if_general_state(plan, pc)
  cb_calls <- 0L
  cb <- function(b) cb_calls <<- cb_calls + 1L
  refuse <- function(expr, class, reads = 0L, calls = 0L, info = NULL) {
    d90_spy_reset(log); cb_calls <<- 0L
    expect_error(expr, class = class, info = info)
    expect_length(log$query_outcomes, 0L)
    expect_identical(length(log$reads), reads, info = info)
    expect_identical(cb_calls, calls, info = info)
  }
  intact <- function(info = NULL) {
    d90_spy_reset(log); cb_calls <<- 0L
    e$replay_ipw_point_if_general(plan, pc, cb, state = st)
    expect_identical(c(cb_calls, length(log$reads), length(log$query_outcomes)), c(C, C * L, 0L), info = info)
  }
  intact("intact state")
  # First callback: the shared runtime marks the state released.
  mut_first <- function(b) {
    cb_calls <<- cb_calls + 1L
    if (cb_calls == 1L) assign("status", "released", envir = st$runtime)
  }
  refuse(e$replay_ipw_point_if_general(plan, pc, mut_first, state = st), "cdid_ipw_point_view_general_error",
    reads = L, calls = 1L, info = "first callback")
  assign("status", "built", envir = st$runtime)
  intact("runtime restored")
  # Last callback: a sealed link state, then a point link, is rebound; the
  # check after the last callback refuses before the replay returns.
  bad_ls <- st$link_states[[L]]
  bad_ls$odds <- rev(bad_ls$odds)
  pl2 <- plan$point_links[[2L]]
  bump <- function(x) { x$estimate <- if (is.null(x$estimate)) 1 else x$estimate + 1; x }
  for (case in list(list(env = st$link_state_refs[[L]], nm = "point_if_general_link_state", value = bad_ls),
                    list(env = pl2$ref, nm = "point_link_view", value = bump(pl2)))) {
    restore <- NULL
    mut_last <- function(b) {
      cb_calls <<- cb_calls + 1L
      if (cb_calls == C) restore <<- gn_rebind(case$env, case$nm, case$value)
    }
    refuse(e$replay_ipw_point_if_general(plan, pc, mut_last, state = st), "cdid_ipw_point_view_general_error",
      reads = C * L, calls = C, info = paste("last callback", case$nm))
    if (is.function(restore)) restore()
    intact(paste("restored", case$nm))
  }
  # An unlocked link-state binding or a rebound plan snapshot refuses before any read.
  unlockBinding("point_if_general_link_state", st$link_state_refs[[1L]])
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, state = st), "cdid_ipw_point_view_general_error", info = "unlocked")
  refuse(e$stream_ipw_point_covariance_general(plan, pc, state = st), "cdid_ipw_point_view_general_error", info = "unlocked, covariance")
  lockBinding("point_if_general_link_state", st$link_state_refs[[1L]])
  plan2 <- e$ipw_point_if_general_plan(w$plinks, pc)
  restore <- gn_rebind(st$binding_ref, "plan_snapshot", plan2)
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, state = st), "cdid_ipw_point_view_general_error", info = "rebound snapshot")
  restore()
  intact("relocked and restored")
  # Foreign plan, foreign state, not a state; a shifted chunk or batch index.
  refuse(e$replay_ipw_point_if_general(plan2, pc, cb, state = st), "cdid_ipw_point_view_general_error", info = "foreign plan")
  st2 <- e$ipw_point_if_general_state(plan2, pc)
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, state = st2), "cdid_ipw_point_view_general_error", info = "foreign state")
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, state = plan), "cdid_ipw_point_view_general_error", info = "not a state")
  refuse(e$ipw_point_if_eval_chunk_internal(plan$norm, 2L, 14L, plan$control, batch_index = 1L, state = st),
    "cdid_ipw_point_view_general_error", info = "shifted chunk")
  refuse(e$ipw_point_if_eval_chunk_internal(plan$norm, 1L, 13L, plan$control, batch_index = 2L, state = st),
    "cdid_ipw_point_view_general_error", info = "wrong batch index")
  # The owner contract is a locked binding plus native SEXP identity
  # (chained_ipw_owner.R:ipw_owner_seal_ok; cdid.c:C_gmm_sexp_aliases).
  # Payload validity is checked during state_build, then protected by COW.
  # Calling the private sealer on invented payloads bypasses that constructor;
  # those fabricated configurations have no replay rejection guarantee.
  # Ordinary aliases work; changing any payload without resealing refuses
  # before a read. No repeated deep validation policy is imposed here.
  same <- st
  expect_identical(gn_collect(e$replay_ipw_point_if_general, plan = plan, control = pc, state = same),
    gn_collect(e$replay_ipw_point_if_general, plan = plan, control = pc))
  edits <- list(
    function(x) { x$odds[1L] <- Inf; x },
    function(x) { x$score_base[1L] <- NaN; x },
    function(x) { x$cohort_index[1L] <- x$cohort_index[1L] + 1L; x },
    function(x) { x$obs_t_index[1L] <- 0L; x },
    function(x) { x$a_index[1:2] <- x$a_index[2:1]; x })
  for (i in seq_along(edits)) {
    changed <- st
    changed$link_states[[1L]] <- edits[[i]](st$link_states[[1L]])
    refuse(e$replay_ipw_point_if_general(plan, pc, cb, state = changed),
      "cdid_ipw_point_view_general_error", info = paste("COW payload mutation", i))
  }
  intact("payload mutations leave the original intact")
  # Construction owns the deep checks. Corrupt a genuine outcome-free query
  # at that boundary to exercise order, endpoint, risk and finite-score
  # refusals, plus failure cleanup, without manufacturing a valid seal.
  bad_queries <- list(
    list(edit = function(x) { x$observed_risk$a_index[1:2] <- x$observed_risk$a_index[2:1]; x },
      class = "cdid_ipw_point_view_general_error"),
    list(edit = function(x) { x$observed_risk$obs_t_index[1L] <- 0L; x },
      class = "cdid_ipw_point_view_general_error"),
    list(edit = function(x) { x$observed_risk$cohort_index[1L] <- x$observed_risk$cohort_index[1L] + 1L; x },
      class = "cdid_ipw_point_view_outcome_leak_error"),
    list(edit = function(x) { x$design$Z_obs[1L, ] <- NaN; x }, class = "cdid_ipw_numeric_error"),
    list(edit = function(x) stop("D90 injected construction allocation failure", call. = FALSE), class = "simpleError"))
  for (i in seq_along(bad_queries)) {
    eb <- gn_ns(list(), character())
    blog <- d90_spy(eb)
    query <- get("query_ipw_view_chunk", envir = eb, inherits = FALSE)
    assign("query_ipw_view_chunk", function(...) bad_queries[[i]]$edit(query(...)), envir = eb)
    reservation <- eb$ipw_point_if_general_state_reserve(plan, pc)
    expect_error(eb$ipw_point_if_general_state_build(reservation, plan, pc),
      class = bad_queries[[i]]$class, info = paste("construction query", i))
    # The fifth injected failure ignores its lazy x argument, so it throws
    # before query(...) is forced; the other corruptions inspect its result.
    expect_identical(blog$query_outcomes, if (i == 5L) logical(0L) else FALSE,
      info = paste("construction query", i))
    expect_length(blog$reads, 0L)
    expect_identical(reservation$state, "failed")
    expect_null(reservation$plan_snapshot)
    expect_true(bindingIsLocked("plan_snapshot", reservation))
    expect_null(blog$state)
  }
  # Released: every later use refuses before any read.
  expect_identical(e$ipw_point_if_general_state_release(st)$status, "retained_and_released")
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, state = st), "cdid_ipw_point_view_general_error", info = "released")
  refuse(e$stream_ipw_point_covariance_general(plan, pc, state = st), "cdid_ipw_point_view_general_error", info = "released, covariance")
  refuse(e$replay_ipw_point_if_general(plan, pc, cb, state = same), "cdid_ipw_point_view_general_error", info = "alias of a released state")
})

# ---------------------------------------------------------------------------
# 5. Admission: the documented composition exactly, adjacent envelope and
#    control caps, exact-integer overflow and malformed
#    requests (refusal), all before any query; the declined record.
# ---------------------------------------------------------------------------

test_that("state admission composes entries and work exactly, declines only at capacity, propagates arithmetic errors and refuses malformed input before querying", {
  w <- d90_w96()
  e <- gn_ns(list(), character())
  log <- d90_spy(e)
  pc <- d90_pc(8L)
  plan <- e$ipw_point_if_general_plan(w$plinks, pc)
  pre <- plan$preflight
  g <- plan$general
  C <- as.double(pre$n_chunks); L <- as.double(plan$n_links); B <- as.double(plan$step); N <- as.double(plan$n_a_units)
  S_ell <- vapply(w$in_obs, sum, integer(1L))
  expect_identical(vapply(plan$point_links, function(lk) as.integer(lk$counts$n_observed_risk_units), integer(1L)), S_ell)
  need <- d90_need(sum(S_ell), C, L, B, pre$total_work_ops_bound, g$one_verification_work_ops_bound)
  expect_identical(d90_plan_need(plan), need)
  # The plan's own caps admit the documented composition.
  d90_spy_reset(log)
  res <- e$ipw_point_if_general_state_reserve(plan, pc)
  expect_true(is.environment(res))
  expect_identical(res$state, "reserved")
  expect_true(environmentIsLocked(res))
  expect_true(all(vapply(setdiff(ls(res, all.names = TRUE), "state"), bindingIsLocked, logical(1L), env = res)))
  expect_false(bindingIsLocked("state", res))
  expect_identical(res$S, sum(as.double(S_ell)))
  expect_identical(c(res$raw_entries, res$entries, res$work, res$state_checks_planned),
    c(need$E, need$entries, need$work, need$checks))
  expect_identical(c(res$construction_work_ops_bound, res$copy_slice_work_ops_bound, res$one_state_check_work_ops_bound,
    res$state_check_work_ops_bound, res$chunk_guard_work_ops_bound, res$setup_work_ops_bound),
    c(as.double(pre$total_work_ops_bound), need$copy_work, need$one_check, need$check_work, need$guard_work, need$setup_work))
  expect_identical(c(res$admission_total, res$work_total), c(need$entries, need$work))
  # Envelope adjacency: exact caps admit; one entry or one operation less declines.
  live0 <- 1000; work0 <- 5000
  ec <- pc
  ec$max_retained_entries <- live0 + need$entries
  ec$max_matrix_entries <- live0 + need$entries
  ec$max_work_ops <- work0 + need$work
  res_x <- e$ipw_point_if_general_state_reserve(plan, pc, live = live0, work = work0, envelope_control = ec)
  expect_true(is.environment(res_x))
  expect_identical(c(res_x$admission_total, res_x$work_total), c(live0 + need$entries, work0 + need$work))
  for (fld in c("max_retained_entries", "max_matrix_entries", "max_work_ops")) {
    low <- ec; low[[fld]] <- low[[fld]] - 1
    expect_null(e$ipw_point_if_general_state_reserve(plan, pc, live = live0, work = work0, envelope_control = low), info = fld)
  }
  expect_length(log$query_outcomes, 0L)
  # The exact-cap reservation builds once and replays as the state-free route.
  st_x <- e$ipw_point_if_general_state_build(res_x, plan, pc)
  expect_identical(res_x$state, "built")
  expect_null(res_x$plan_snapshot)
  expect_true(bindingIsLocked("plan_snapshot", res_x))
  expect_identical(c(st_x$entries, st_x$admission_total), c(need$entries, live0 + need$entries))
  expect_identical(gn_collect(e$replay_ipw_point_if_general, plan = plan, control = pc, state = st_x),
    gn_collect(e$replay_ipw_point_if_general, plan = plan, control = pc))
  expect_error(e$ipw_point_if_general_state_build(res_x, plan, pc), class = "cdid_ipw_point_view_general_error")
  plan2 <- e$ipw_point_if_general_plan(w$plinks, pc)
  expect_error(e$ipw_point_if_general_state_build(e$ipw_point_if_general_state_reserve(plan, pc), plan2, pc),
    class = "cdid_ipw_point_view_general_error")
  # Control adjacency: the covariance stream's concurrent entries and total
  # work with the state; declined, the unchanged covariance runs under the cap.
  v <- as.double(g$one_verification_work_ops_bound)
  need_cov_e <- pre$peak_concurrent_entries + 4 * L * L + 8 * L + 6 * pre$b_by_L_entries + g$additional_live_entries +
    need$entries + 128
  need_cov_w <- 2 * pre$total_work_ops_bound + (1 + 2 * (C + 1)) * v + N * (2 * L * L + 4 * L) + need$work + 128
  need_producer_e <- pre$peak_concurrent_entries + g$additional_live_entries + need$entries
  need_replay_w <- pre$total_work_ops_bound + g$additional_work_ops_bound + need$work
  cov_ref <- e$stream_ipw_point_covariance_general(plan, pc)$covariance
  for (adj in list(list(fld = "max_retained_entries", need = need_producer_e),
                   list(fld = "max_matrix_entries", need = max(need_cov_e, need_producer_e)),
                   list(fld = "max_work_ops", need = max(need_cov_w, need_replay_w)))) {
    pc_x <- pc; pc_x[[adj$fld]] <- adj$need
    pc_1 <- pc; pc_1[[adj$fld]] <- adj$need - 1
    plan_x <- e$ipw_point_if_general_plan(w$plinks, pc_x)
    plan_1 <- e$ipw_point_if_general_plan(w$plinks, pc_1)
    d90_spy_reset(log)
    expect_true(is.environment(e$ipw_point_if_general_state_reserve(plan_x, pc_x)), info = adj$fld)
    expect_null(e$ipw_point_if_general_state_reserve(plan_1, pc_1), info = adj$fld)
    expect_length(log$query_outcomes, 0L)
    expect_identical(e$stream_ipw_point_covariance_general(plan_1, pc_1)$covariance, cov_ref, info = adj$fld)
  }
  # Exact arithmetic failures propagate through reservation, including when a
  # valid capacity cap would also be exceeded. Only requested > bound declines.
  big <- 2^53 - 1
  d90_spy_reset(log)
  expect_error(e$ipw_point_if_general_state_entries(big, 1, 1, 1), class = "cdid_ipw_owner_bound_error")
  expect_error(e$ipw_point_if_general_state_entries(1, big, 1, 1), class = "cdid_ipw_owner_bound_error")
  expect_error(e$ipw_point_if_general_state_entries(1.5, 1, 1, 1), class = "cdid_ipw_owner_bound_error")
  expect_error(e$ipw_point_if_general_state_work(big, 1, 1, 1), class = "cdid_ipw_owner_bound_error")
  expect_identical(e$ipw_point_if_general_state_entries(sum(as.double(S_ell)), C, L, B)$reserved_entries, need$entries)
  expect_identical(e$ipw_point_if_general_state_work(pre$total_work_ops_bound, need$E, C, L,
    plan_check_work = g$one_verification_work_ops_bound)$total_work, need$work)
  expect_error(e$ipw_point_if_general_state_reserve(plan, pc, live = big - 1), class = "cdid_ipw_owner_bound_error")
  expect_error(e$ipw_point_if_general_state_reserve(plan, pc, work = big - 1), class = "cdid_ipw_owner_bound_error")
  for (bad in list(1.5, 2^53)) {
    expect_error(e$ipw_point_if_general_state_reserve(plan, pc, live = bad), class = "cdid_ipw_owner_bound_error")
    expect_error(e$ipw_point_if_general_state_reserve(plan, pc, work = bad), class = "cdid_ipw_owner_bound_error")
  }
  for (bad in list(-1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(e$ipw_point_if_general_state_reserve(plan, pc, live = bad), class = "cdid_ipw_point_view_input_error")
    expect_error(e$ipw_point_if_general_state_reserve(plan, pc, work = bad), class = "cdid_ipw_point_view_input_error")
  }
  expect_error(e$ipw_point_if_general_state_reserve(plan, pc, envelope_control = list(max_retained_entries = 1)),
    class = "cdid_ipw_point_view_input_error")
  for (fld in c("max_retained_entries", "max_matrix_entries", "max_output_entries", "max_work_ops")) {
    bad_ec <- pc; bad_ec[[fld]] <- 1.5
    expect_error(e$ipw_point_if_general_state_reserve(plan, pc, envelope_control = bad_ec),
      class = "cdid_ipw_owner_bound_error", info = fld)
  }
  # An allocation failure inside the capped inventory must propagate, too.
  ef <- gn_ns(list(), character())
  flog <- d90_spy(ef)
  entries <- get("ipw_point_if_general_state_entries", envir = ef, inherits = FALSE)
  assign("ipw_point_if_general_state_entries", function(S, C, L, B, max_lim = NULL) {
    if (!is.null(max_lim)) stop("D90 injected inventory allocation failure", call. = FALSE)
    entries(S, C, L, B, max_lim = max_lim)
  }, envir = ef)
  expect_error(ef$ipw_point_if_general_state_reserve(plan, pc), "D90 injected inventory allocation failure", class = "simpleError")
  expect_length(flog$query_outcomes, 0L)
  expect_length(flog$reads, 0L)
  expect_error(e$ipw_point_if_general_state_reserve(plan), class = "cdid_ipw_point_view_input_error")
  pc_chunk <- pc; pc_chunk$max_chunk_rows <- pc$max_chunk_rows + 1L
  expect_error(e$ipw_point_if_general_state_reserve(plan, pc_chunk), class = "cdid_ipw_control_mismatch_error")
  expect_error(e$ipw_point_if_general_state(plan, pc_chunk), class = "cdid_ipw_control_mismatch_error")
  expect_error(e$ipw_point_if_general_state_reserve(w$plinks, pc), class = "cdid_ipw_point_view_general_error")
  expect_length(log$query_outcomes, 0L)
  # The declined record states the same requirement and no work done.
  dec <- e$ipw_point_if_general_state_release(NULL, plan)
  expect_identical(dec$status, "declined_capacity")
  expect_false(isTRUE(dec$admitted))
  expect_identical(c(dec$raw_state_entries, dec$reserved_entries, dec$work_ops_bound, dec$state_checks_planned),
    c(need$E, need$entries, need$work, need$checks))
  expect_identical(c(dec$queries_performed, dec$passes_served, dec$chunks_served), c(0L, 0L, 0L))
})
