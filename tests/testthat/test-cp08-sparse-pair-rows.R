# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 general-panel performance: the native sparse two-date pair-row join of
# `cp04_pair_rows` (`C_gmm_cp08_sparse_pair_rows`). Expected values come from a
# frozen verbatim copy of the previous R helper, from an independent keyed join
# on the raw long rows, and from hand-computed search paths; the kernel never
# supplies its own truth. Full OR/DR general-route fits, IF chunks and
# covariances are compared against the same build running the frozen helper.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

sp_ctrl <- function(chunk = 7L, work = 1e10, matrix = 2e7) api$ipw_owner_control(max_chunk_rows = chunk,
  max_retained_entries = 1e7, max_matrix_entries = matrix, max_work_ops = work, max_output_entries = 1e7)
sp_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
sp_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))
sp_g <- function(id) ifelse(id %% 2L == 0L, Inf, 3)

# Frozen verbatim copy of the previous `cp04_pair_rows` (before the native
# join). It is evaluated in the installed namespace (or a clone of it).
sp_old_pair_rows <- function(reader, rows, pair, read_outcomes = FALSE) {
  st <- reader$state
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  n <- length(rows)
  out <- list(S = numeric(n))
  if (read_outcomes) out$dy <- numeric(n)
  ai <- st$a_index[rows]
  keep <- which(ai > 0L)
  ds <- pair$s_calendar_index; dt <- pair$t_calendar_index
  if (is.null(ds)) ds <- match(as.double(pair$s), panel$all_times)
  if (is.null(dt)) dt <- match(as.double(pair$t), panel$all_times)
  post <- panel$postings
  if (length(keep) && !is.na(ds) && !is.na(dt) &&
      post$date_count[ds] > 0L && post$date_count[dt] > 0L) {
    pa <- post$post_a_index
    ms <- ipw_owner_match_sorted_double(ai[keep], pa, post$date_start[ds], post$date_end[ds])
    hs <- which(!is.na(ms))
    if (length(hs)) {
      mt <- ipw_owner_match_sorted_double(ai[keep[hs]], pa, post$date_start[dt], post$date_end[dt])
      ht <- which(!is.na(mt))
      if (length(ht)) {
        obs <- keep[hs[ht]]
        out$S[obs] <- 1
        if (read_outcomes) {
          ps <- post$post_obs_index[ms[hs[ht]]]
          pt <- post$post_obs_index[mt[ht]]
          yy <- ipw_owner_read_y_positions(panel, ps, pt)
          out$dy[obs] <- yy$y_right_obs - yy$y_left_obs
        }
      }
    }
  }
  if (read_outcomes) cp04_finite(out$dy, "observed pair differences")
  out
}
environment(sp_old_pair_rows) <- api

# Independent keyed join on the raw long rows: a frame id is an observed pair
# row of (s, t) exactly when it has a non-NA outcome row at both dates (two
# observed dates imply A membership); dy is that unit's y_t - y_s.
sp_reference <- function(d, ids, s, t) {
  # paste(numeric(0), scalar) produces one key; an empty query has no rows.
  if (!length(ids)) return(list(S = numeric(0), dy = numeric(0)))
  key <- paste(d$id, d$t)
  ys <- d$y[match(paste(ids, s), key)]
  yt <- d$y[match(paste(ids, t), key)]
  S <- as.double(!is.na(ys) & !is.na(yt))
  dy <- numeric(length(ids))
  dy[S == 1] <- yt[S == 1] - ys[S == 1]
  list(S = S, dy = dy)
}

# Long rows of N regular units on T dates under a missingness pattern, plus a
# singleton (one observed row), an all-NA unit (explicit NA rows only), a unit
# observed only on an extra calendar date T + 1 (an A-free date) and, when the
# pattern is general, explicit NA outcomes inside otherwise observed units.
# Unit 3 has exact zero outcomes at dates 1 and 2 (a zero difference).
sp_long <- function(N = 48L, T = 5L, pattern = "half", seed = 3L) {
  set.seed(seed)
  i <- rep(seq_len(N), each = T); t <- rep(seq_len(T), N)
  keep <- switch(pattern,
    balanced = rep(TRUE, N * T),
    half = (i + t) %% 3L != 0L,
    rotating = { start <- (i %% 3L) + 1L; t >= start & t <= start + 2L },
    general = ((i * 7L + t * 3L) %% 5L != 0L) & stats::runif(N * T) > 0.2,
    stop("unknown pattern"))
  keep[i == 3L & t <= 2L] <- TRUE
  d <- data.frame(id = i[keep], t = t[keep])
  d$g <- sp_g(d$id)
  d$y <- 0.3 * d$id + 0.7 * d$t + ifelse(d$t >= d$g, 1.5, 0) + stats::rnorm(nrow(d))
  if (pattern == "general") d$y[sample(which(d$id != 3L), 6L)] <- NA_real_
  d$y[d$id == 3L & d$t %in% c(1L, 2L)] <- 0
  extra <- data.frame(id = c(N + 1L, N + 2L, N + 2L, N + 3L), t = c(2L, 1L, 2L, T + 1L))
  extra$g <- sp_g(extra$id)
  extra$y <- c(1.25, NA_real_, NA_real_, -0.5)
  rbind(d, extra)
}

# Real sealed owners and the full-F reader over every id 1..N + 3 (ids absent
# from the panel under sparse patterns are frame-only rows), optionally in a
# shuffled frame order or with character ids.
sp_fixture <- function(d = sp_long(), chunk = 7L, perm = FALSE, char_ids = FALSE, ctrl = sp_ctrl(chunk), seed = 5L) {
  n_ids <- max(d$id)
  ids <- seq_len(n_ids)
  g <- sp_g(ids)
  if (perm) { set.seed(seed); ord <- sample(n_ids) } else ord <- ids
  ids <- ids[ord]; g <- g[ord]
  if (char_ids) { d$id <- sprintf("u%03d", d$id); ids <- sprintf("u%03d", ids) }
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  w <- rep(1, n_ids)
  Z <- matrix(1, n_ids, 1L, dimnames = list(NULL, "(Intercept)"))
  V <- matrix(1, n_ids, 1L, dimnames = list(NULL, "reader_intercept"))
  fo <- api$construct_ipw_frame_owner(po, ids, Z, weights = w, control = ctrl, frame_g = g, allow_frame_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, ids, V, g, w, ctrl)
  list(d = d, po = po, fo = fo, reader = reader, ids = ids, g = g, ctrl = ctrl, chunk = chunk, N = n_ids)
}
sp_pair <- function(fx, s, t) api$ipw_owner_make_pair_descriptor(fx$po, s, t)
sp_views <- function(fx, pairs) lapply(seq_along(pairs), function(k) {
  p <- pairs[[k]]
  api$construct_ipw_link_view(fx$fo, p[1], p[2], g = 3, K = Inf, control = fx$ctrl,
    link_id = paste0("L", k, ":", p[1], "_", p[2]), functional_key = paste0("cmpcat:sp:", k, ":", p[1], "-", p[2]))
})

# `cp04_pair_rows` re-bound so that its own `.Call` of the sparse join kernel is
# counted (the resolver's and other kernels' calls are not).
sp_counting <- function() {
  e <- new.env(parent = api)
  f <- api$cp04_pair_rows; environment(f) <- e; assign("cp04_pair_rows", f, e)
  n <- new.env(parent = emptyenv()); n$calls <- 0L
  addr <- api$cp04_pair_rows_native_kernel()
  e$.Call <- function(.NAME, ...) {
    if (identical(.NAME, addr)) n$calls <- n$calls + 1L
    base::.Call(.NAME, ...)
  }
  e$n <- n
  e
}

# A clone of the installed namespace in which every closure is re-bound to the
# clone and `cp04_pair_rows` is the frozen helper: the previous general route
# of this same build, for fit/IF/covariance comparisons.
sp_frozen_ns <- function() {
  e <- new.env(parent = api)
  for (nm in ls(api, all.names = TRUE)) {
    f <- get(nm, envir = api, inherits = FALSE)
    if (is.function(f) && !is.primitive(f) && identical(environment(f), api)) {
      environment(f) <- e; assign(nm, f, envir = e)
    }
  }
  for (nm in c("gmm_native_dll_cache", "cp04_payload_native_cache", "or_baseline_native_cache", "cp04_pair_rows_native_cache"))
    assign(nm, new.env(parent = emptyenv()), envir = e)
  old <- sp_old_pair_rows; environment(old) <- e
  assign("cp04_pair_rows", old, envir = e)
  e
}

# A reader-shaped list over an edited copy of a real panel owner: only the
# fields `cp04_pair_rows` reads (state$a_index, panel_ref, control).
sp_fake_reader <- function(fx, edit_panel = identity, edit_a = identity) {
  panel <- edit_panel(fx$po)
  env <- new.env(parent = emptyenv()); assign("panel_owner", panel, envir = env)
  st <- fx$reader$state; st$a_index <- edit_a(st$a_index)
  list(state = st, panel_ref = env, control = fx$reader$control)
}

sp_check <- function(fx, rows, pair, read, info = NULL) {
  new <- api$cp04_pair_rows(fx$reader, rows, pair, read)
  old <- sp_old_pair_rows(fx$reader, rows, pair, read)
  expect_identical(new, old, info = info)
  ref <- sp_reference(fx$d, fx$ids[rows], pair$s, pair$t)
  expect_identical(new$S, ref$S, info = info)
  if (read) expect_identical(new$dy, ref$dy, info = info) else expect_null(new$dy, info = info)
  expect_identical(typeof(new$S), "double"); expect_true(all(new$S %in% c(0, 1)))
  invisible(new)
}

test_that("the native join equals the frozen helper and an independent keyed join on every pattern", {
  for (case in list(list(pattern = "balanced"), list(pattern = "half"), list(pattern = "half", perm = TRUE),
                    list(pattern = "rotating", perm = TRUE, chunk = 16L), list(pattern = "general", char = TRUE),
                    list(pattern = "general", perm = TRUE, chunk = 1L))) {
    fx <- sp_fixture(sp_long(pattern = case$pattern), chunk = if (is.null(case$chunk)) 7L else case$chunk,
      perm = isTRUE(case$perm), char_ids = isTRUE(case$char))
    N <- fx$N; st <- fx$reader$state
    info <- paste(case$pattern, isTRUE(case$perm), isTRUE(case$char))
    set.seed(11L)
    row_sets <- list(seq_len(N), rev(seq_len(N)), sample(N), 1L, N, integer(0), seq_len(fx$chunk),
      seq.int(N - fx$chunk + 1L, N), c(3L, 3L, 5L, 1L, 5L), as.double(seq_len(min(N, 7L))),
      which(st$a_index == 0L), which(fx$g == 3), which(is.infinite(fx$g)), sample(N, 9L))
    pairs <- c(lapply(list(c(1, 2), c(1, 4), c(2, 3), c(2, 5), c(3, 4), c(1, 5), c(4, 5), c(1, 6), c(2, 6)),
      function(p) sp_pair(fx, p[1], p[2])), list(list(s = 2.5, t = 4), list(s = 1, t = 7), list(s = 0, t = 3)))
    for (pair in pairs) for (rows in row_sets) for (read in c(FALSE, TRUE))
      sp_check(fx, rows, pair, read, info = paste(info, pair$s, pair$t, length(rows), read))
    # The A-free extra date 6 and absent dates give all-zero S on every row.
    for (pair in list(sp_pair(fx, 1, 6), list(s = 2.5, t = 4))) {
      hit <- api$cp04_pair_rows(fx$reader, seq_len(N), pair, TRUE)
      expect_identical(hit, list(S = numeric(N), dy = numeric(N)))
    }
    # Singleton, all-NA and frame-only rows are never observed; unit 3 has an
    # observed zero difference at (1, 2), distinct from a non-observed zero.
    off <- which(st$a_index == 0L)
    expect_gte(length(off), 3L)
    hit <- api$cp04_pair_rows(fx$reader, seq_len(N), sp_pair(fx, 1, 2), TRUE)
    expect_identical(hit$S[off], numeric(length(off)))
    r3 <- match(if (isTRUE(case$char)) "u003" else 3L, fx$ids)
    expect_identical(c(hit$S[r3], hit$dy[r3]), c(1, 0))
    # Wrong unit mapping or a reversed difference would be detected: observed
    # differences are nonzero and not symmetric, and the reference is sensitive
    # to the id order it is handed.
    obs <- which(hit$S == 1 & hit$dy != 0)
    expect_gte(length(obs), 2L)
    expect_false(identical(hit$dy[obs], -hit$dy[obs]))
    shuffled <- sp_reference(fx$d, fx$ids[c(obs[2L], obs[1L], obs[-(1:2)])], 1, 2)$dy
    expect_false(identical(shuffled, sp_reference(fx$d, fx$ids[obs], 1, 2)$dy))
    # The native kernel is called once per query with both dates present and a
    # nonempty on-A row set, never for off-A-only/empty rows or absent/A-free dates.
    ec <- sp_counting()
    ec$cp04_pair_rows(fx$reader, seq_len(N), sp_pair(fx, 1, 2), TRUE)
    ec$cp04_pair_rows(fx$reader, off, sp_pair(fx, 1, 2), FALSE)
    expect_identical(ec$n$calls, 1L)
    ec$cp04_pair_rows(fx$reader, integer(0), sp_pair(fx, 1, 2), TRUE)
    ec$cp04_pair_rows(fx$reader, seq_len(N), sp_pair(fx, 1, 6), TRUE)
    ec$cp04_pair_rows(fx$reader, seq_len(N), list(s = 2.5, t = 4), TRUE)
    expect_identical(ec$n$calls, 1L)
  }
})

test_that("kernel positions follow the scalar search path and malformed arguments refuse before any read", {
  k <- api$cp04_pair_rows_native_kernel()
  call <- function(a, pa, po, bounds, read = TRUE, max = 1e7) .Call(k, a, pa, po, bounds, read, max)
  # Date 1 occupies postings 1..4, date 2 postings 5..9 with a duplicated A-unit 6.
  pa <- c(1L, 3L, 4L, 6L, 2L, 3L, 6L, 6L, 7L); po <- 101:109; bd <- c(1L, 4L, 5L, 9L)
  a <- c(6L, 0L, 3L, 5L, 6L, 1L, 7L)
  res <- call(a, pa, po, bd)
  # Hand-computed paths (mid = lo + floor((hi - lo) / 2), first equal stops):
  # a = 6 -> positions 4 and 7 (the first duplicate the path meets), a = 3 -> 2
  # and 6, a = 5 / 1 / 7 absent at one date, a = 0 off A.
  expect_identical(res, list(obs = c(1L, 3L, 5L), pos_s = c(104L, 102L, 104L), pos_t = c(107L, 106L, 107L)))
  # The same positions as the retained scalar matcher chain of the old helper.
  keep <- which(a > 0L)
  ms <- api$ipw_owner_match_sorted_double(a[keep], pa, bd[1L], bd[2L]); hs <- which(!is.na(ms))
  mt <- api$ipw_owner_match_sorted_double(a[keep[hs]], pa, bd[3L], bd[4L]); ht <- which(!is.na(mt))
  expect_identical(res, list(obs = keep[hs[ht]], pos_s = po[ms[hs[ht]]], pos_t = po[mt[ht]]))
  expect_identical(call(a, pa, po, bd, read = FALSE), list(obs = c(1L, 3L, 5L), pos_s = NULL, pos_t = NULL))
  expect_identical(call(integer(0), pa, po, bd), list(obs = integer(0), pos_s = integer(0), pos_t = integer(0)))
  expect_identical(call(c(0L, 0L), pa, po, bd), list(obs = integer(0), pos_s = integer(0), pos_t = integer(0)))
  expect_identical(call(6L, pa, po, c(4L, 4L, 7L, 7L)), list(obs = 1L, pos_s = 104L, pos_t = 107L))
  expect_identical(call(6L, pa, po, c(4L, 4L, 8L, 8L)), list(obs = 1L, pos_s = 104L, pos_t = 108L))
  # Unsorted or duplicated subranges yield the position of the same path as
  # the scalar search (including its path-dependent misses), never a scan.
  pa2 <- c(5L, 2L, 9L, 2L, 7L); q <- c(2L, 7L, 9L, 5L)
  ms2 <- api$ipw_owner_match_sorted_double(q, pa2, 1, 5)
  # Paths: 2 -> 9, 5 -> miss; 7 -> 9, 5, 2 -> miss; 9 found at 3; 5 -> 9, 5 found at 1.
  expect_identical(ms2, c(NA_integer_, NA_integer_, 3L, 1L))
  r2 <- call(q, pa2, 201:205, c(1L, 5L, 1L, 5L))
  expect_identical(r2, list(obs = which(!is.na(ms2)), pos_s = (201:205)[ms2[!is.na(ms2)]], pos_t = (201:205)[ms2[!is.na(ms2)]]))
  # Malformed arguments abort with the typed bound error before any posting is read.
  bad <- list(list(a = as.double(a)), list(a = c(a, NA_integer_)), list(a = c(a, -1L)),
    list(pa = structure(pa, class = "sp_tag")), list(po = po[-1L]), list(po = as.double(po)),
    list(bounds = c(0L, 4L, 5L, 9L)), list(bounds = c(1L, 4L, 5L, 10L)), list(bounds = c(4L, 1L, 5L, 9L)),
    list(bounds = c(1L, 4L, 9L, 5L)), list(bounds = c(1L, 4L, 5L)), list(bounds = as.double(bd)),
    list(bounds = c(1L, NA_integer_, 5L, 9L)), list(read = NA), list(read = "yes"), list(read = c(TRUE, TRUE)),
    list(max = 20), list(max = -1), list(pa = integer(0), po = integer(0)))
  for (ed in bad) {
    args <- list(a = a, pa = pa, po = po, bounds = bd)
    for (nm in names(ed)) args[[nm]] <- ed[[nm]]
    expect_error(do.call(call, args), class = "cdid_gmm_bound_error", info = paste(names(ed), collapse = ","))
  }
  # n=7 scratch21 plus k=3/readTRUE output9 =30 simultaneous entries.
  expect_error(call(a, pa, po, bd, max = 29), class = "cdid_gmm_bound_error")
  expect_identical(call(a, pa, po, bd, max = 30)$obs, c(1L, 3L, 5L))
  # An NA posting value on the visited path is refused rather than ordered.
  pa_na <- pa; pa_na[7L] <- NA_integer_
  expect_error(call(6L, pa_na, po, bd), class = "cdid_gmm_bound_error")
  expect_identical(call(3L, pa_na, po, c(1L, 4L, 5L, 6L)), list(obs = 1L, pos_s = 102L, pos_t = 106L))
  # Registration, forced symbol, arity and the cached live address.
  dll <- getLoadedDLLs()[["cdid"]]
  reg <- getDLLRegisteredRoutines(dll, addNames = TRUE)[[".Call"]]
  expect_identical(reg[["C_gmm_cp08_sparse_pair_rows"]]$numParameters, 6L)
  expect_identical(k, getNativeSymbolInfo("C_gmm_cp08_sparse_pair_rows", PACKAGE = dll)$address)
  expect_identical(api$cp04_pair_rows_native_kernel(), k)
  cloned <- api$cp04_pair_rows_native_kernel; environment(cloned) <- new.env(parent = api)
  expect_identical(cloned(), k)
  expect_error(api$or_baseline_native_kernel("C_gmm_cp08_sparse_pair_rows"), class = "cdid_gmm_native_error")
})

test_that("structurally unsupported inputs keep the retained scalar route before any native call", {
  fx <- sp_fixture(sp_long(pattern = "general"), perm = TRUE)
  N <- fx$N; post <- fx$po$postings; st <- fx$reader$state
  pair <- sp_pair(fx, 1, 4)
  bounds <- c(post$date_start[pair$s_calendar_index], post$date_end[pair$s_calendar_index],
    post$date_start[pair$t_calendar_index], post$date_end[pair$t_calendar_index])
  ai <- st$a_index[seq_len(N)]
  ok <- api$cp04_pair_rows_native_ok
  expect_true(ok(ai, post, bounds, TRUE))
  expect_true(ok(integer(0), post, bounds, FALSE))
  expect_false(ok(as.double(ai), post, bounds, TRUE))
  expect_false(ok(replace(ai, 2L, NA_integer_), post, bounds, TRUE))
  expect_false(ok(structure(ai, class = "sp_tag"), post, bounds, TRUE))
  expect_false(ok(ai, replace(post, "post_a_index", list(as.double(post$post_a_index))), bounds, TRUE))
  expect_false(ok(ai, replace(post, "post_obs_index", list(as.double(post$post_obs_index))), bounds, TRUE))
  expect_false(ok(ai, replace(post, "post_a_index", list(structure(post$post_a_index, class = "sp_tag"))), bounds, TRUE))
  expect_false(ok(ai, post, as.double(bounds), TRUE))
  expect_false(ok(ai, post, replace(bounds, 4L, NA_integer_), TRUE))
  expect_false(ok(ai, post, bounds[c(2L, 1L, 3L, 4L)], TRUE))
  expect_false(ok(ai, post, bounds[1:3], TRUE))
  expect_false(ok(ai, post, bounds, NA))
  expect_false(ok(ai, post, bounds, 1))
  expect_false(ok(ai, post, bounds, c(TRUE, FALSE)))
  # Each unsupported structure takes the scalar route (no native call) and
  # returns exactly the frozen helper's S/dy, including the long-table shape
  # (double observation positions) and NA A-indices treated as off-A rows.
  fakes <- list(
    long_positions = sp_fake_reader(fx, function(p) { p$postings$post_obs_index <- as.double(p$postings$post_obs_index); p }),
    double_postings = sp_fake_reader(fx, function(p) { p$postings$post_a_index <- as.double(p$postings$post_a_index); p }),
    classed_postings = sp_fake_reader(fx, function(p) { p$postings$post_a_index <- structure(p$postings$post_a_index, class = "sp_tag"); p }),
    double_bounds = sp_fake_reader(fx, function(p) { p$postings$date_start <- as.double(p$postings$date_start); p }),
    na_a_index = sp_fake_reader(fx, edit_a = function(a) { a[which(a > 0L)[1:2]] <- NA_integer_; a }))
  supported <- sp_fake_reader(fx)
  ec <- sp_counting()
  for (nm in names(fakes)) for (rows in list(seq_len(N), c(5L, 1L, 9L))) for (read in c(FALSE, TRUE)) {
    before <- ec$n$calls
    new <- ec$cp04_pair_rows(fakes[[nm]], rows, pair, read)
    expect_identical(ec$n$calls, before, info = nm)
    expect_identical(new, sp_old_pair_rows(fakes[[nm]], rows, pair, read), info = nm)
    expect_identical(new, api$cp04_pair_rows(fakes[[nm]], rows, pair, read), info = nm)
    if (nm != "na_a_index") expect_identical(new, api$cp04_pair_rows(supported, rows, pair, read), info = nm)
  }
  # The supported shape of the same panel does call the kernel and agrees.
  before <- ec$n$calls
  expect_identical(ec$cp04_pair_rows(supported, seq_len(N), pair, TRUE), sp_old_pair_rows(supported, seq_len(N), pair, TRUE))
  expect_identical(ec$n$calls, before + 1L)
  # A native failure propagates as the kernel's typed error; nothing is caught.
  # The smallest A-unit of date s is the first posting of its subrange, so its
  # search path always ends on that (now NA) position.
  broken <- sp_fake_reader(fx, function(p) { p$postings$post_a_index[bounds[1L]] <- NA_integer_; p })
  one <- which(st$a_index == post$post_a_index[bounds[1L]])[1L]
  expect_false(is.na(one))
  expect_error(api$cp04_pair_rows(broken, one, pair, TRUE), class = "cdid_gmm_bound_error")
  expect_error(api$cp04_pair_rows(broken, c(one, which(st$a_index > 0L)), pair, FALSE), class = "cdid_gmm_bound_error")
})

test_that("general-route OR and DR fits, full-id IF chunks and covariances equal the frozen-helper build", {
  frozen <- sp_frozen_ns()
  pairs <- list(c(1, 2), c(1, 4), c(2, 3), c(2, 5), c(3, 4))
  for (case in list(list(perm = FALSE, chunk = 7L), list(perm = TRUE, chunk = 16L))) {
    fx <- sp_fixture(sp_long(pattern = "half"), chunk = case$chunk, perm = case$perm)
    N <- fx$N; views <- sp_views(fx, pairs); L <- length(views)
    for (method in c("or", "dr")) {
      fit <- function(ns) if (method == "or") ns$fit_or_baseline(fx$reader, views, fx$ctrl, 1e-10, outcome_recipe = sp_icpt) else
        ns$fit_dr_baseline(fx$reader, views, fx$ctrl, glm_control = sp_glm, rank_tol = 1e-10, outcome_recipe = sp_icpt)
      fn <- fit(api); ff <- fit(frozen)
      # Incomplete dates: both fits take the general route through cp04_pair_rows.
      expect_null(fn$complete_support); expect_null(ff$complete_support)
      expect_identical(fn$supports, ff$supports)
      expect_identical(fn$regressions, ff$regressions)
      expect_identical(fn$shared$beta, ff$shared$beta)
      if (method == "dr") {
        expect_identical(lapply(fn$treatments, `[[`, "coefficients"), lapply(ff$treatments, `[[`, "coefficients"))
        expect_identical(lapply(fn$treatments, `[[`, "means"), lapply(ff$treatments, `[[`, "means"))
      }
      chunk <- function(ns, f, rows, links, scale = "frame", details = FALSE) if (method == "or")
        ns$evaluate_or_baseline_if_chunk(f, rows, links, scale, details) else
        ns$evaluate_dr_baseline_if_chunk(f, rows, links, scale, details)
      set.seed(17L)
      row_sets <- c(lapply(seq.int(1L, N, by = case$chunk), function(s) seq.int(s, min(N, s + case$chunk - 1L))),
        list(sample(N, case$chunk), N, c(2L, 1L)))
      for (rows in row_sets) for (links in list(seq_len(L), c(3L, 1L), L)) for (scale in c("frame", "core")) {
        a <- chunk(api, fn, rows, links, scale); b <- chunk(frozen, ff, rows, links, scale)
        expect_identical(a$phi, b$phi); expect_identical(a$id, b$id); expect_identical(a$rows, b$rows)
        expect_identical(a$links, b$links); expect_identical(a$normalization_count, b$normalization_count)
        expect_identical(a$id, fx$reader$state$id[rows])
      }
      a <- chunk(api, fn, seq_len(min(N, case$chunk)), seq_len(L), "frame", TRUE)
      b <- chunk(frozen, ff, seq_len(min(N, case$chunk)), seq_len(L), "frame", TRUE)
      expect_identical(a$phi, b$phi); expect_identical(a$details, b$details)
      cov <- function(ns, f) if (method == "or") ns$cp04_stream_covariance(f, seq_along(f$links)) else ns$stream_dr_covariance(f)
      cn <- cov(api, fn); cf <- cov(frozen, ff)
      for (nm in c("estimates", "covariance", "mean_phi", "chunk_rows")) expect_identical(cn[[nm]], cf[[nm]], info = nm)
      expect_true(all(is.finite(cn$estimates))); expect_true(all(is.finite(cn$covariance)))
      expect_identical(dim(cn$covariance), c(L, L))
    }
  }
})


test_that("sparse admission preserves no-A short circuit, negative legacy rows and bounded ALTREP access", {
  fx <- sp_fixture(sp_long(pattern = "half"))
  off <- which(fx$reader$state$a_index == 0L)
  expect_gt(length(off), 0L)
  thin <- sp_fake_reader(fx)
  thin$control$max_matrix_entries <- 1
  pair <- sp_pair(fx, 1, 2)
  ec <- sp_counting()
  expect_identical(ec$cp04_pair_rows(thin, off, pair, TRUE), sp_old_pair_rows(thin, off, pair, TRUE))
  expect_identical(ec$n$calls, 0L)
  neg <- sp_fake_reader(fx, edit_a = function(a) { a[1L] <- -1L; a })
  expect_identical(ec$cp04_pair_rows(neg, seq_len(fx$N), pair, TRUE), sp_old_pair_rows(neg, seq_len(fx$N), pair, TRUE))
  expect_identical(ec$n$calls, 0L)
  k <- api$cp04_pair_rows_native_kernel()
  # Both tables are compact sequences; a tiny query reads only search paths.
  a <- seq_len(1000000L)
  got <- .Call(k, 17L, a, a, c(1L, 1000000L, 1L, 1000000L), TRUE, 6)
  expect_identical(got, list(obs = 1L, pos_s = 17L, pos_t = 17L))
  expect_error(.Call(k, 17L, a, a, c(1L, 1000000L, 1L, 1000000L), TRUE, 5), class = "cdid_gmm_bound_error")
  long <- seq_len(2^31)
  expect_false(api$cp04_pair_rows_native_ok(1L, list(post_a_index = long, post_obs_index = long), c(1L, 1L, 1L, 1L), FALSE))
})


test_that("NULL native entry bound is unbounded while finite live-output boundaries stay exact", {
  k <- api$cp04_pair_rows_native_kernel()
  pa <- c(1L, 3L, 4L, 6L, 2L, 3L, 6L, 6L, 7L); po <- 101:109
  bd <- c(1L, 4L, 5L, 9L); a <- c(6L, 0L, 3L, 5L, 6L, 1L, 7L)
  call <- function(read, bound) .Call(k, a, pa, po, bd, read, bound)
  for (read in c(FALSE, TRUE)) {
    expected <- list(obs = c(1L, 3L, 5L),
      pos_s = if (read) c(104L, 102L, 104L) else NULL,
      pos_t = if (read) c(107L, 106L, 107L) else NULL)
    expect_identical(call(read, NULL), expected)
    expect_identical(call(read, if (read) 30 else 24), expected)
  }
  expect_error(call(TRUE, 29), class = "cdid_gmm_bound_error")
  expect_identical(call(TRUE, 30)$obs, c(1L, 3L, 5L))
  expect_error(.Call(k, 6L, pa, po, bd, TRUE, 5), class = "cdid_gmm_bound_error")
  expect_identical(.Call(k, 6L, pa, po, bd, TRUE, 6), list(obs = 1L, pos_s = 104L, pos_t = 107L))
})
