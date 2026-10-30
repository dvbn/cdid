# SPDX-License-Identifier: GPL-3.0-only
#
# D90 missing-panel candidate: a non-detailed general (missing-panel) baseline
# DR IF chunk reads only core$beta, so evaluate_dr_baseline_if_chunk() asks
# or_pooled_if_core() for the beta-only core at p_or > 0 (the same global
# shared-beta score and solve, without block-alpha, direct or standalone-OR
# link phi) and skips the core and its support outcome reads at p_or = 0.
# Checked against the genuine full-core route (the unchanged detailed route,
# and the same caller fed the full core) at the existing 1e-12, with spies
# proving that the discarded core work is not performed. D92 also authorizes
# non-detailed corrected-selection and dated DR consumers. Detailed DR,
# default OR and every other core caller keep the full core.
#
# STATUS: UNEXECUTED at authoring time (base-R parse only; no package load,
# build, install, fit or test body was run in the authoring job). The root
# runs it against the installed candidate.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures (the supplied CP09 panel and owner construction; deterministic)
# ---------------------------------------------------------------------------

d90_ctrl <- function(chunk = 9L)
  api$ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = 2e7,
    max_matrix_entries = 2e7, max_work_ops = 1e13, max_output_entries = 2e7)
d90_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)

# cp09_panel(): 53 units, nonadjacent dates 2, 4, 9, 15, 26, 40 and cohorts
# 9, 26 and never treated (0).
d90_panel <- function() {
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

# 16 missing cells in every cohort and at every date (at most two per unit;
# every (pair, cohort) keeps at least 11 observed units), including both cells
# of the supplied missing-panel case.
d90_drop <- function(d) {
  tt <- sort(unique(d$period))
  d[!((d$id %% 4L == 1L & d$period == tt[1L + (d$id %/% 4L) %% 6L]) |
      (d$id == 1053L & d$period == 4) | (d$id == 1052L & d$period == 15)), , drop = FALSE]
}

# p0: intercept only (p_or = 0); p2: intercept and two shared slopes; p2_noint:
# two shared slopes without intercept (the or_wls_gram_solve branch).
d90_recipe <- function(recipe) switch(recipe,
  p0 = api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0)),
  p2 = api$or_pooled_outcome_recipe(c("(Intercept)", "x", "x2"), TRUE, c("x", "x2")),
  p2_noint = api$or_pooled_outcome_recipe(c("x", "x2"), FALSE, c("x", "x2")))

d90_fixture <- function(cg = "notyettreated", recipe = "p2", missing = TRUE, perm = FALSE,
                        wt = 1, est = "dr", ctrl = d90_ctrl()) {
  d <- d90_panel()
  if (missing) d <- d90_drop(d)
  d$G[d$G == 0] <- Inf
  if (perm) d <- d[order(cos(seq_len(nrow(d)) * 0.7)), ]
  po <- api$construct_ipw_panel_owner(d, "id", "period", "y", "G", ctrl)
  u <- d[!duplicated(d$id), ]; u <- u[order(sin(u$id * 0.3)), ]
  X <- cbind(`(Intercept)` = 1, x = u$x, x2 = u$x2)
  Z <- if (est == "dr") cbind(X, x3 = u$x3) else X
  ww <- if (length(wt) == 1L) rep(wt, nrow(u)) else wt[match(u$id, sort(u$id))]
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = ww, control = ctrl,
    frame_g = u$G, allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, u$id,
    matrix(1, nrow(u), 1L, dimnames = list(NULL, "reader_intercept")), u$G, ww, ctrl,
    allow_unknown_off_A = TRUE)
  tt <- sort(unique(d$period)); views <- list()
  for (g in c(9, 26)) for (s in seq_len(length(tt) - 1L)) for (t in seq.int(s + 1L, length(tt))) {
    K <- if (cg == "nevertreated") Inf else c(9, 26, Inf)[c(9, 26, Inf) > max(g, tt[t])]
    j <- length(views) + 1L
    views[[j]] <- api$construct_ipw_link_view(fo, tt[s], tt[t], g = g, K = K, control = ctrl,
      link_id = paste0("L", j), functional_key = paste0("cmpcat:d90:", j))
  }
  if (perm) views <- views[c(seq.int(length(views), 1L, by = -2L), seq.int(length(views) - 1L, 1L, by = -2L))]
  list(reader = reader, views = views, recipe = d90_recipe(recipe), ctrl = ctrl, est = est,
    treatment_columns = c("(Intercept)", "x", "x3"))
}

# Forced general (missing-panel) route on a complete panel (supplied rm_general).
d90_general <- local({
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$or_baseline_complete_facts <- function(...) NULL
  e
})

d90_fit <- function(fx, ns = api) {
  if (identical(fx$est, "or"))
    return(ns$fit_or_baseline(fx$reader, fx$views, fx$ctrl, rank_tol = 1e-10, outcome_recipe = fx$recipe))
  ns$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = d90_glm, rank_tol = 1e-10,
    outcome_recipe = fx$recipe, treatment_columns = fx$treatment_columns)
}

# ---------------------------------------------------------------------------
# Spies
# ---------------------------------------------------------------------------

# Clone of one installed core caller whose or_pooled_if_core() (and, with
# pair_spy, cp04_pair_rows() and its outcome endpoint reader) resolve to
# recording spies; everything else resolves in the locked installed namespace.
# full_core = TRUE feeds the caller the genuine full core whatever it requests.
d90_spy_ns <- function(caller, full_core = FALSE, pair_spy = FALSE) {
  e <- new.env(parent = api)
  f <- get(caller, api); environment(f) <- e; assign(caller, f, e)
  n <- new.env(parent = emptyenv())
  n$core <- list(); n$pairs <- list(); n$y_reads <- 0
  e$or_pooled_if_core <- function(model, N, w, hits, up, pre, links, details = FALSE, beta_only = FALSE) {
    supplied <- !missing(beta_only)
    out <- if (full_core) api$or_pooled_if_core(model, N, w, hits, up, pre, links, details) else
      api$or_pooled_if_core(model, N, w, hits, up, pre, links, details, beta_only = beta_only)
    n$core[[length(n$core) + 1L]] <- list(supplied = supplied, beta_only = beta_only, details = details,
      names = names(out), out = out, args = list(N = N, w = w, hits = hits, up = up, pre = pre, links = links))
    out
  }
  if (pair_spy) {
    pr <- get("cp04_pair_rows", api); environment(pr) <- e
    e$cp04_pair_rows <- function(reader, rows, pair, read_outcomes = FALSE) {
      out <- pr(reader, rows, pair, read_outcomes)
      n$pairs[[length(n$pairs) + 1L]] <- c(read = as.double(read_outcomes), observed = sum(out$S))
      out
    }
    e$ipw_owner_read_y_positions <- function(panel_owner, pos_s, pos_t, ...) {
      n$y_reads <- n$y_reads + length(pos_s)
      api$ipw_owner_read_y_positions(panel_owner, pos_s, pos_t, ...)
    }
  }
  e$n <- n
  e
}

# Records since the last take (pair reads in call order), then reset.
d90_take <- function(e) {
  n <- e$n
  pairs <- if (length(n$pairs)) do.call(rbind, n$pairs) else
    matrix(numeric(0), 0L, 2L, dimnames = list(NULL, c("read", "observed")))
  out <- list(core = n$core, pairs = pairs, y_reads = n$y_reads)
  n$core <- list(); n$pairs <- list(); n$y_reads <- 0
  out
}

# Access-recording view of a fitted backend: every read of a top-level field
# and of a regression-block field is logged, so a direct core call shows which
# fitted inputs (hence which of its loops) it touched.
d90_recording_model <- function(model, log) {
  bind <- function(x, prefix, swap = list()) {
    env <- new.env(parent = emptyenv())
    for (nm in setdiff(names(x), "")) local({
      key <- nm
      value <- if (is.null(swap[[key]])) x[[key]] else swap[[key]]
      makeActiveBinding(key, function() {
        log$fields <- c(log$fields, paste0(prefix, key))
        value
      }, env)
    })
    env
  }
  bind(model, "", list(regressions = lapply(model$regressions, bind, prefix = "regressions$")))
}

# Every or_pooled_if_core() call expression within a function body.
d90_core_calls <- function(expr) {
  found <- list()
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (identical(e[[1L]], as.name("or_pooled_if_core"))) found[[length(found) + 1L]] <<- e
    for (i in seq_along(e)) if (!identical(e[[i]], quote(expr = ))) walk(e[[i]])
    invisible(NULL)
  }
  walk(expr)
  found
}

# One general-route chunk: the candidate non-detailed route against the
# genuine full-core routes, with the discarded-work spies. Returns the
# candidate frame IF and both routes' outcome endpoint read counts.
d90_check_chunk <- function(fit, rows, links, lab) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  p <- fit$shared$p
  pre <- api$dr_baseline_if_envelope(fit, length(rows), links, FALSE)
  pre_det <- api$dr_baseline_if_envelope(fit, length(rows), links, TRUE)
  e <- d90_spy_ns("evaluate_dr_baseline_if_chunk", pair_spy = TRUE)
  a <- e$evaluate_dr_baseline_if_chunk(fit, rows, links, scale = "frame", details = FALSE)
  act <- d90_take(e)
  dd <- e$evaluate_dr_baseline_if_chunk(fit, rows, links, scale = "frame", details = TRUE)
  det <- d90_take(e)
  # The spies are transparent: the installed routes return the same objects.
  expect_identical(a, api$evaluate_dr_baseline_if_chunk(fit, rows, links, scale = "frame", details = FALSE), info = lab)
  expect_identical(dd, api$evaluate_dr_baseline_if_chunk(fit, rows, links, scale = "frame", details = TRUE), info = lab)
  expect_identical(api$evaluate_dr_baseline_if_chunk(fit, rows, links, scale = "core", details = FALSE)$phi,
    st$lambda * a$phi, info = lab)
  expect_identical(a$diagnostics[c("output_entries_bound", "concurrent_scratch_bound", "work_ops_bound")],
    list(output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work),
    info = lab)

  # Genuine full-core route 1: the unchanged detailed chunk, whose core still
  # forms every output (link phi, direct, alpha, scores, objective).
  expect_equal(a$phi, dd$phi, tolerance = 1e-12, info = lab)
  expect_identical(length(det$core), 1L, info = lab)
  rd <- det$core[[1L]]
  expect_true(identical(rd$beta_only, FALSE) && identical(rd$details, TRUE), info = lab)
  expect_true(all(c("phi", "direct", "beta", "alpha", "beta_score", "beta_block_scores", "alpha_score",
    "denominator", "numerator", "design", "objective") %in% rd$names), info = lab)
  expect_identical(dd$details$beta, rd$out$beta, info = lab)
  expect_identical(dd$details$alpha, rd$out$alpha, info = lab)
  expect_identical(dd$details$direct, rd$out$direct, info = lab)
  expect_identical(dd$details$psi_or_estimator, rd$out$phi, info = lab)
  expect_true(all(det$pairs[, "read"] == 1), info = lab)

  if (p > 0L) {
    # The caller requests and receives only xi_beta, still reading every
    # support's endpoints that xi_beta needs.
    expect_identical(length(act$core), 1L, info = lab)
    rec <- act$core[[1L]]
    expect_true(rec$supplied && identical(rec$beta_only, TRUE) && identical(rec$details, FALSE), info = lab)
    expect_identical(rec$names, "beta", info = lab)
    expect_true(all(act$pairs[, "read"] == 1), info = lab)
    expect_equal(rec$out$beta, dd$details$beta, tolerance = 1e-12, info = lab)
    expect_identical(a$diagnostics$n_global_beta_blocks, length(fit$regression_specs), info = lab)
    # On the captured inputs the beta-only core never reads a link record or a
    # block's intercept forcing/mass (its alpha, direct and standalone-OR phi
    # loops are not entered); the full core does; both give the same xi_beta.
    log_b <- new.env(parent = emptyenv()); log_b$fields <- character()
    log_f <- new.env(parent = emptyenv()); log_f$fields <- character()
    ar <- rec$args
    fast <- api$or_pooled_if_core(d90_recording_model(fit, log_b), ar$N, ar$w, ar$hits, ar$up, ar$pre,
      ar$links, FALSE, beta_only = TRUE)
    full <- api$or_pooled_if_core(d90_recording_model(fit, log_f), ar$N, ar$w, ar$hits, ar$up, ar$pre,
      ar$links, FALSE)
    expect_identical(fast, rec$out, info = lab)
    expect_identical(full$beta, rec$out$beta, info = lab)
    expect_identical(names(full), c("phi", "direct", "beta", "alpha"), info = lab)
    expect_true(all(c("shared", "outcome_recipe", "regression_specs", "supports", "regressions$L",
      "regressions$alpha") %in% log_b$fields), info = lab)
    expect_false("links" %in% log_b$fields, info = lab)
    expect_false(any(c("regressions$forcing_alpha", "regressions$mass") %in% log_b$fields), info = lab)
    expect_true("links" %in% log_f$fields, info = lab)
    if (isTRUE(fit$outcome_recipe$intercept))
      expect_true(all(c("regressions$forcing_alpha", "regressions$mass") %in% log_f$fields), info = lab)
    # Genuine full-core route 2: the same caller consuming the full core.
    ef <- d90_spy_ns("evaluate_dr_baseline_if_chunk", full_core = TRUE)
    af <- ef$evaluate_dr_baseline_if_chunk(fit, rows, links, scale = "frame", details = FALSE)
    fr <- d90_take(ef)$core
    expect_identical(length(fr), 1L, info = lab)
    expect_identical(fr[[1L]]$names, c("phi", "direct", "beta", "alpha"), info = lab)
    expect_equal(a$phi, af$phi, tolerance = 1e-12, info = lab)
  } else {
    # p_or = 0: no core call; every support pair is read for its observation
    # indicator only, and outcome endpoints are read exactly for the treatment
    # models' pairs, as in the full-core detailed route, which also reads every
    # support's endpoints.
    expect_identical(length(act$core), 0L, info = lab)
    ns <- length(pre$sids); nd <- length(pre_det$sids)
    expect_identical(unname(act$pairs[seq_len(ns), "read"]), rep(0, ns), info = lab)
    trt <- act$pairs[-seq_len(ns), , drop = FALSE]
    expect_true(all(trt[, "read"] == 1), info = lab)
    expect_identical(trt, det$pairs[-seq_len(nd), , drop = FALSE], info = lab)
    expect_identical(act$y_reads, sum(trt[, "observed"]), info = lab)
    expect_identical(det$y_reads, sum(det$pairs[, "observed"]), info = lab)
    expect_identical(a$diagnostics$n_regression_score_evaluations, 0L, info = lab)
  }
  list(phi = a$phi, y_act = act$y_reads, y_det = det$y_reads)
}

# ---------------------------------------------------------------------------
# 1. General/missing baseline DR: beta-only (p_or > 0) and skipped (p_or = 0)
#    cores against the genuine full-core route, NT/NYT, permuted panel rows,
#    links and chunk rows, varying weights, nonadjacent link batches.
# ---------------------------------------------------------------------------
test_that("general/missing baseline DR IF: beta-only core (p_or > 0) and skipped core (p_or = 0) match the genuine full-core route at 1e-12 and spies prove the discarded core work is not performed", {
  expect_true(environmentIsLocked(api))
  perm_rows <- c(41L, 3L, 17L, 29L, 8L, 53L, 22L, 11L, 36L)
  sub_links <- c(2L, 7L, 13L, 19L, 30L)
  cases <- list(
    list(missing = TRUE, cg = "nevertreated", perm = FALSE, wt = 1),
    list(missing = TRUE, cg = "notyettreated", perm = TRUE, wt = 1 + (seq_len(53L) %% 5) / 8),
    list(missing = FALSE, cg = "notyettreated", perm = FALSE, wt = 1),
    list(missing = FALSE, cg = "nevertreated", perm = TRUE, wt = 1))
  for (cs in cases) for (recipe in c("p0", "p2", "p2_noint")) {
    lab <- paste(if (cs$missing) "missing" else "forced-general", cs$cg, recipe,
      if (cs$perm) "permuted" else "ordered", if (length(cs$wt) > 1L) "weighted" else "unit-weight")
    fx <- d90_fixture(cg = cs$cg, recipe = recipe, missing = cs$missing, perm = cs$perm, wt = cs$wt)
    fit <- d90_fit(fx, ns = if (cs$missing) api else d90_general)
    expect_true(inherits(fit, "cdid_dr_baseline"), info = lab)
    expect_null(fit$complete_support, info = lab)
    expect_null(fit$reference_ref, info = lab)
    expect_equal(fit$shared$p, if (recipe == "p0") 0 else 2, info = lab)
    expect_identical(length(fit$links), 30L, info = lab)
    N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
    expect_identical(N, 53L, info = lab)
    H <- matrix(NA_real_, N, 30L)
    y <- c(act = 0, det = 0)
    for (rows in split(seq_len(N), ceiling(seq_len(N) / 9L))) {
      chk <- d90_check_chunk(fit, rows, seq_len(30L), lab)
      H[rows, ] <- unname(chk$phi)
      y <- y + c(chk$y_act, chk$y_det)
    }
    expect_false(anyNA(H), info = lab)
    # Permuted chunk rows and a nonadjacent link batch give the same complete
    # unit rows as the ordered full-width chunks.
    chk <- d90_check_chunk(fit, perm_rows, sub_links, lab)
    expect_equal(unname(chk$phi), H[perm_rows, sub_links, drop = FALSE], tolerance = 1e-12, info = lab)
    # At p_or = 0 the skipped support endpoint reads are actual savings.
    if (recipe == "p0") expect_gt(y[["det"]], y[["act"]])
  }
})

# ---------------------------------------------------------------------------
# 2. The flag is private and isolated: the guard refuses detailed or
#    non-logical requests, the default is the full core, and only baseline and
#    corrected-selection DR chunks (including dated references) pass beta_only.
# ---------------------------------------------------------------------------
test_that("beta-only core guard refuses detailed or invalid requests and every other or_pooled_if_core caller keeps the full default core", {
  fx <- d90_fixture(cg = "notyettreated", recipe = "p2", missing = TRUE)
  fit <- d90_fit(fx)
  e <- d90_spy_ns("evaluate_dr_baseline_if_chunk")
  invisible(e$evaluate_dr_baseline_if_chunk(fit, 1:9, 1:30, scale = "frame", details = FALSE))
  ar <- d90_take(e)$core[[1L]]$args
  core <- function(...) api$or_pooled_if_core(fit, ar$N, ar$w, ar$hits, ar$up, ar$pre, ar$links, ...)
  for (bad in list(list(details = TRUE, beta_only = TRUE), list(details = NA, beta_only = TRUE),
                   list(details = FALSE, beta_only = NA), list(details = FALSE, beta_only = 1),
                   list(details = FALSE, beta_only = c(TRUE, TRUE))))
    expect_error(do.call(core, bad), class = "cdid_or_pooled_input_error")
  full <- core(FALSE)
  expect_identical(names(full), c("phi", "direct", "beta", "alpha"))
  expect_identical(core(FALSE, beta_only = FALSE), full)
  expect_identical(core(FALSE, beta_only = TRUE), list(beta = full$beta))
  expect_identical(names(core(TRUE)), c("phi", "direct", "beta", "alpha", "beta_score", "beta_block_scores",
    "alpha_score", "denominator", "numerator", "design", "objective"))

  callers <- character()
  for (nm in ls(api, all.names = TRUE)) {
    f <- get(nm, envir = api, inherits = FALSE)
    if (!is.function(f) || is.primitive(f) || identical(nm, "or_pooled_if_core")) next
    for (cl in d90_core_calls(body(f))) {
      callers <- c(callers, nm)
      tags <- names(cl)
      if (is.null(tags)) tags <- character()
      if (nm %in% c("evaluate_dr_baseline_if_chunk", "evaluate_dr_corollary1_if_chunk")) {
        expect_identical(length(cl), 10L)
        expect_identical(tags[10L], "beta_only")
      } else {
        expect_false("beta_only" %in% tags, info = nm)
        expect_true(length(cl) <= 9L, info = nm)
      }
    }
  }
  expect_true(all(c("evaluate_dr_baseline_if_chunk", "evaluate_or_baseline_if_chunk",
    "evaluate_dr_corollary1_if_chunk") %in% callers))
  expect_identical(sum(callers == "evaluate_dr_baseline_if_chunk"), 1L)
  expect_identical(sum(callers == "evaluate_dr_corollary1_if_chunk"), 1L)
})

# ---------------------------------------------------------------------------
# 3. Default baseline OR on the missing panel consumes the full core's phi.
# ---------------------------------------------------------------------------
test_that("default baseline OR on the missing panel keeps the full default core in non-detailed and detailed chunks", {
  fx <- d90_fixture(cg = "notyettreated", recipe = "p2", missing = TRUE, est = "or")
  fit <- d90_fit(fx)
  expect_true(inherits(fit, "cdid_or_baseline"))
  expect_null(fit$complete_support)
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  links <- seq_along(fit$links)
  e <- d90_spy_ns("evaluate_or_baseline_if_chunk")
  chunks <- c(split(seq_len(N), ceiling(seq_len(N) / 9L)), list(c(41L, 3L, 17L, 29L, 8L, 53L, 22L, 11L, 36L)))
  for (det in c(FALSE, TRUE)) for (rows in chunks) {
    a <- e$evaluate_or_baseline_if_chunk(fit, rows, links, scale = "frame", details = det)
    rec <- d90_take(e)$core
    expect_identical(a, api$evaluate_or_baseline_if_chunk(fit, rows, links, scale = "frame", details = det))
    expect_identical(length(rec), 1L)
    expect_false(rec[[1L]]$supplied)
    expect_false(rec[[1L]]$beta_only)
    expect_true(all(c("phi", "direct", "beta", "alpha") %in% rec[[1L]]$names))
    expect_identical(unname(a$phi), unname(rec[[1L]]$out$phi))
  }
})

# ---------------------------------------------------------------------------
# 4. Corrected-selection (Corollary 1) DR on the supplied F7 fixture.
# ---------------------------------------------------------------------------
test_that("corrected-selection DR on the supplied F7 fixture uses beta-only non-detailed chunks and preserves the full detailed core", {
  fx_dir <- testthat::test_path("fixtures", "dr_finite_v1")
  loader <- new.env(parent = globalenv())
  sys.source(file.path(fx_dir, "load_fixture.R"), envir = loader)
  fx7 <- loader$load_dr_finite_reference(fx_dir, fixture_id = "F7_date_product_F_within_cohort",
    order = NULL, max_chunk_rows = 13L)
  fit7 <- api$fit_dr_backend(reader = fx7$reader, link_views = fx7$views, control = fx7$control,
    glm_control = fx7$glm_control, rank_tol = fx7$rank_tol, outcome_recipe = fx7$outcome_recipe,
    treatment_columns = fx7$treatment_columns, selection = fx7$selection,
    selection_model = fx7$selection_model, selection_pooling = fx7$selection_pooling, sel_fit = NULL)
  expect_true(inherits(fit7, "cdid_dr_corollary1"))
  N7 <- get("cp04_reader", fit7$reader_ref, inherits = FALSE)$state$N
  links7 <- seq_along(fit7$links)
  e7 <- d90_spy_ns("evaluate_dr_corollary1_if_chunk")
  ef7 <- d90_spy_ns("evaluate_dr_corollary1_if_chunk", full_core = TRUE)
  for (det in c(FALSE, TRUE)) for (rows in split(seq_len(N7), ceiling(seq_len(N7) / 13L))) {
    a <- e7$evaluate_dr_corollary1_if_chunk(fit7, rows, links7, scale = "frame", details = det)
    rec <- d90_take(e7)$core
    expect_identical(a, api$evaluate_dr_corollary1_if_chunk(fit7, rows, links7, scale = "frame", details = det))
    expect_identical(length(rec), 1L)
    expect_true(rec[[1L]]$supplied)
    expect_identical(rec[[1L]]$beta_only, !det)
    af <- ef7$evaluate_dr_corollary1_if_chunk(fit7, rows, links7, scale = "frame", details = det)
    ref <- d90_take(ef7)$core
    expect_identical(length(ref), 1L)
    expect_identical(a, af)
    expect_identical(rec[[1L]]$out$beta, ref[[1L]]$out$beta)
    if (det) {
      expect_true(all(c("phi", "direct", "beta", "alpha") %in% rec[[1L]]$names))
      expect_identical(a$details$beta, ref[[1L]]$out$beta)
      expect_identical(a$details$alpha, ref[[1L]]$out$alpha)
      expect_identical(a$details$direct, ref[[1L]]$out$direct)
      expect_identical(a$details$psi_or_estimator, ref[[1L]]$out$phi)
    } else {
      expect_identical(rec[[1L]]$names, "beta")
    }
  }
})
