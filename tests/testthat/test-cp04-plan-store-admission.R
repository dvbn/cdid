# SPDX-License-Identifier: GPL-3.0-only
#
# The covariance plan store is optional: it is admitted against the workspace
# limits before it is allocated, holds only the fields its route reads, is not
# live once the fit's own arrays are allocated, and never changes a plan, an
# admission decision, a refusal or a result.
library(testthat)
api <- asNamespace("cdid")

store_panel <- function() {
  set.seed(41)
  n <- 80L
  d <- expand.grid(id = seq_len(n), t = 1:6)
  d$g <- rep(c(0, 3, 5), length.out = n)[d$id]
  d$x <- stats::rnorm(n)[d$id]
  d$y <- d$x + 0.3 * d$t + stats::rnorm(nrow(d)) + (d$g > 0 & d$t >= d$g)
  d$y[stats::runif(nrow(d)) < 0.15] <- NA_real_
  d
}
store_fit <- function(...) {
  base <- list(yname = "y", tname = "t", idname = "id", gname = "g", data = store_panel(), xformla = ~x,
    control_group = "nevertreated", gmm_weighting = "within", bstrap = FALSE, cband = FALSE)
  do.call(cdid::att_gt, utils::modifyList(base, list(...)))
}
# Run f with the store switched off (tests only) or on, from the same seed.
store_run <- function(off, f) {
  sw <- api$cp04_plan_scope
  assign("store_off", off, envir = sw)
  on.exit(assign("store_off", NULL, envir = sw), add = TRUE)
  assign("last", NULL, envir = api$cp04_plan_path_record)
  set.seed(101)
  r <- tryCatch(suppressWarnings(suppressMessages(f())), error = identity)
  list(r = r, paths = api$cp04_plan_path_record$last)
}
store_same <- function(a, b) {
  if (inherits(a, "error")) {
    expect_identical(class(b), class(a))
    expect_identical(conditionMessage(b), conditionMessage(a))
  } else {
    x <- unclass(a); y <- unclass(b); x$call <- y$call <- NULL
    expect_identical(y, x)
  }
}

test_that("the store of a narrow-block OR fit holds only the fields its route reads and is gone before the solve", {
  # This fit plans 30 links in one-link blocks, 465 block pairs. The
  # row planner does not read envelopes on this route (explicit chunk rows),
  # so only the three budget vectors (3 x 465 doubles) may be allocated.
  seen <- new.env()
  seen$slots <- integer(); seen$live <- logical()
  assign(".store_seen", seen, envir = globalenv())
  on.exit(rm(".store_seen", envir = globalenv()), add = TRUE)
  fields <- c("one_output", "one_scratch", "two_output", "two_scratch", "rows", "work", "peak")
  suppressMessages(trace("cp04_stream_plan_legacy", where = api, print = FALSE, exit = bquote({
    mm <- cp04_plan_scope$current$memo
    .store_seen$slots <- c(.store_seen$slots, if (is.environment(mm))
      sum(vapply(.(fields), function(k) length(mm[[k]]), 1L)) else 0L)
  })))
  on.exit(suppressMessages(untrace("cp04_stream_plan_legacy", where = api)), add = TRUE)
  suppressMessages(trace("solve_chained_gmm_core", where = api, print = FALSE, tracer = quote(
    .store_seen$live <- c(.store_seen$live, is.environment(cp04_plan_scope$current) &&
      !is.null(cp04_plan_scope$current$memo)))))
  on.exit(suppressMessages(untrace("solve_chained_gmm_core", where = api)), add = TRUE)
  fit <- store_fit(est_method = "reg", control = list(chunk_rows = 1L))
  expect_true(length(fit$att) > 0L)
  expect_true(length(seen$slots) > 0L)
  expect_true(all(seen$slots == 3L * 465L))
  expect_true(length(seen$live) > 0L)
  expect_false(any(seen$live))
})

test_that("store on and store off give identical fits, refusals and draws at the same limit", {
  cases <- list(
    # narrow blocks (one link per block), explicit chunk rows, bootstrap draws
    or_narrow = function() store_fit(est_method = "reg", bstrap = TRUE, biters = 49L, control = list(chunk_rows = 1L)),
    # envelope store on the row planner's classes, several link blocks
    dr_pooled = function() store_fit(est_method = "dr", outcome_pooling = "pooled", bstrap = TRUE, biters = 99L,
      control = list(numeric_workspace_bytes = 24 * 2^20)),
    # a refusing fit after the store was used
    dr_pooled_refused = function() store_fit(est_method = "dr", outcome_pooling = "pooled", bstrap = TRUE, biters = 99L,
      control = list(numeric_workspace_bytes = 12 * 2^20)),
    # pooled OR, optimal weighting, simultaneous band
    or_pooled_optimal = function() store_fit(est_method = "reg", outcome_pooling = "pooled",
      control_group = "notyettreated", gmm_weighting = "optimal", bstrap = TRUE, biters = 49L, cband = TRUE,
      control = list(chunk_rows = 2L)))
  for (nm in names(cases)) {
    on <- store_run(FALSE, cases[[nm]])
    off <- store_run(TRUE, cases[[nm]])
    store_same(on$r, off$r)
    # The same planning paths; the store was used with it on and declined with it off.
    keep <- c("certified", "complete_single_block", "complete_full_block", "complete_block_pairs")
    expect_identical(on$paths[keep], off$paths[keep])
    expect_gt(on$paths[["store_budget_hits"]], 0)
    expect_identical(on$paths[["store_declined"]], 0)
    expect_gt(off$paths[["store_declined"]], 0)
    expect_identical(off$paths[["store_budget_hits"]] + off$paths[["store_pair_hits"]], 0)
  }
  expect_true(inherits(store_run(FALSE, cases$dr_pooled_refused)$r, "cdid_error"))
})

test_that("a store that does not fit is declined and the plan is the plan built without a store", {
  # Capture the backend and links of the first covariance plan of the narrow-block fit.
  box <- new.env()
  assign(".store_box", box, envir = globalenv())
  on.exit(rm(".store_box", envir = globalenv()), add = TRUE)
  suppressMessages(trace("cp04_stream_plan", where = api, print = FALSE, tracer = quote(
    if (is.null(.store_box$backend)) { .store_box$backend <- backend; .store_box$links <- links })))
  on.exit(suppressMessages(untrace("cp04_stream_plan", where = api)), add = TRUE)
  store_fit(est_method = "reg", control = list(chunk_rows = 1L))
  b <- box$backend; links <- box$links
  L <- length(links); K <- L
  need <- api$cp04_plan_memo_entries(K, L, FALSE)
  expect_identical(need, 3 * K * (K + 1) / 2 + K * (K + 1) / 2 + L + 256)
  cap <- min(b$control$max_matrix_entries, b$control$max_retained_entries)
  own <- b$accounting$owned_entries + b$accounting$borrowed_entries
  # The planner admits this backend's plan at these retained charges; the
  # store fits at the first and not at the second (one entry more).
  for (r in c(cap - own - need, cap - own - need + 1)) {
    plain <- tryCatch(api$cp04_stream_plan_legacy(b, links, r), error = identity)
    prev <- api$cp04_plan_scope_enter()
    scoped <- tryCatch(api$cp04_stream_plan_legacy(b, links, r), error = identity)
    again <- tryCatch(api$cp04_stream_plan_legacy(b, links, r), error = identity)
    memo <- api$cp04_plan_scope$current$memo
    paths <- api$cp04_plan_scope$current$paths
    api$cp04_plan_scope_exit(prev)
    expect_false(inherits(plain, "error"))
    expect_identical(scoped, plain)
    expect_identical(again, plain)
    if (r + own + need <= cap) {
      expect_false(is.null(memo))
      expect_identical(paths[["store_declined"]], 0)
      expect_identical(paths[["store_budget_hits"]], K * (K + 1) / 2)
    } else {
      expect_null(memo)
      expect_identical(paths[["store_declined"]], 2)
      expect_identical(paths[["store_budget_hits"]], 0)
    }
  }
})

test_that("a stored plan is admitted again before it is reused", {
  box <- new.env()
  assign(".store_box", box, envir = globalenv())
  on.exit(rm(".store_box", envir = globalenv()), add = TRUE)
  suppressMessages(trace("cp04_stream_plan", where = api, print = FALSE, tracer = quote(
    if (is.null(.store_box$backend)) { .store_box$backend <- backend; .store_box$links <- links })))
  on.exit(suppressMessages(untrace("cp04_stream_plan", where = api)), add = TRUE)
  store_fit(est_method = "reg", control = list(chunk_rows = 1L))
  b <- box$backend; links <- box$links
  L <- length(links)
  need <- api$cp04_plan_memo_entries(L, L, FALSE)
  cap <- min(b$control$max_matrix_entries, b$control$max_retained_entries)
  own <- b$accounting$owned_entries + b$accounting$borrowed_entries
  r0 <- cap - own - need
  plain <- api$cp04_stream_plan_legacy(b, links, r0 + 1)
  prev <- api$cp04_plan_scope_enter()
  first <- api$cp04_stream_plan_legacy(b, links, r0)
  expect_false(is.null(api$cp04_plan_scope$current$memo))
  # The same key at one entry less of room: the kept store no longer fits.
  second <- api$cp04_stream_plan_legacy(b, links, r0 + 1)
  memo <- api$cp04_plan_scope$current$memo
  paths <- api$cp04_plan_scope$current$paths
  api$cp04_plan_scope_exit(prev)
  expect_identical(second, plain)
  expect_identical(first$blocks, plain$blocks)
  expect_null(memo)
  expect_identical(paths[["store_declined"]], 1)
  expect_identical(paths[["store_budget_hits"]], 0)
})

test_that("the plan store admitted in the covariance stream is charged with the bound retained moment IF store", {
  seen <- new.env()
  seen$pairs <- list(); seen$cov <- NULL
  assign(".store_seen", seen, envir = globalenv())
  on.exit(rm(".store_seen", envir = globalenv()), add = TRUE)
  suppressMessages(trace("cp04_stream_covariance", where = api, print = FALSE, tracer = quote({
    bd <- cp04_retained_if_unbind(backend)
    .store_seen$cov <- list(retained = retained, store = if (is.null(bd$store)) NA_real_ else bd$store$entries)
  })))
  on.exit(suppressMessages(untrace("cp04_stream_covariance", where = api)), add = TRUE)
  suppressMessages(trace("cp04_plan_memo_admit", where = api, print = FALSE, tracer = quote(
    if (!is.null(.store_seen$cov)) {
      .store_seen$pairs[[length(.store_seen$pairs) + 1L]] <- c(.store_seen$cov$retained, .store_seen$cov$store, retained)
      .store_seen$cov <- NULL
    })))
  on.exit(suppressMessages(untrace("cp04_plan_memo_admit", where = api)), add = TRUE)
  fit <- store_fit(est_method = "reg", control = list(chunk_rows = 1L))
  expect_true(length(fit$att) > 0L)
  bound <- Filter(function(x) is.finite(x[2]), seen$pairs)
  # This fit binds a retained moment IF store to the covariance stream, whose
  # plan admits a store.
  expect_true(length(bound) >= 1L)
  for (x in bound) expect_identical(x[3], x[1] + x[2])
  expect_null(api$cp04_plan_scope$current)
})

test_that("an inner pipeline restores the outer plan scope and its store, also when it is refused", {
  box <- new.env()
  assign(".store_box", box, envir = globalenv())
  on.exit(rm(".store_box", envir = globalenv()), add = TRUE)
  suppressMessages(trace("cp04_stream_plan", where = api, print = FALSE, tracer = quote(
    if (is.null(.store_box$backend)) { .store_box$backend <- backend; .store_box$links <- links })))
  store_fit(est_method = "reg", control = list(chunk_rows = 1L))
  suppressMessages(untrace("cp04_stream_plan", where = api))
  prev <- api$cp04_plan_scope_enter()
  outer <- api$cp04_plan_scope$current
  api$cp04_stream_plan_legacy(box$backend, box$links, 0)
  memo <- outer$memo
  expect_false(is.null(memo))
  fit <- store_fit(est_method = "dr", outcome_pooling = "pooled", control = list(numeric_workspace_bytes = 24 * 2^20))
  expect_true(length(fit$att) > 0L)
  expect_identical(api$cp04_plan_scope$current, outer)
  expect_identical(outer$memo, memo)
  refused <- tryCatch(store_fit(est_method = "dr", outcome_pooling = "pooled", bstrap = TRUE, biters = 99L,
    control = list(numeric_workspace_bytes = 12 * 2^20)), error = identity)
  expect_true(inherits(refused, "cdid_error"))
  expect_identical(api$cp04_plan_scope$current, outer)
  expect_identical(outer$memo, memo)
  api$cp04_plan_scope_exit(prev)
  expect_null(api$cp04_plan_scope$current)
})

test_that("the store's size is checked before allocation, also where it would be material", {
  box <- new.env()
  assign(".store_box", box, envir = globalenv())
  on.exit(rm(".store_box", envir = globalenv()), add = TRUE)
  suppressMessages(trace("cp04_stream_plan", where = api, print = FALSE, tracer = quote(
    if (is.null(.store_box$backend)) .store_box$backend <- backend)))
  on.exit(suppressMessages(untrace("cp04_stream_plan", where = api)), add = TRUE)
  store_fit(est_method = "dr", outcome_pooling = "pooled", control = list(numeric_workspace_bytes = 24 * 2^20))
  b <- box$backend
  # K = 45 link blocks: 7 x 1035 doubles = 57,960 bytes of fields.
  expect_identical(api$cp04_plan_memo_entries(45, 975, TRUE), 8 * 1035 + 975 + 256)
  # K = 2251: 7 x 2,534,626 doubles (about 142 MB) of fields, 20.3 million
  # entries with the copy allowance; not admitted under a 64 MiB limit
  # (8,388,608 entries), admitted under 1 GiB with nothing else retained.
  big <- api$cp04_plan_memo_entries(2251, 2251, TRUE)
  expect_identical(big, 8 * 2251 * 2252 / 2 + 2251 + 256)
  at <- function(bytes) {
    bb <- b
    bb$control$max_matrix_entries <- bb$control$max_retained_entries <- floor(bytes / 8)
    bb
  }
  expect_false(api$cp04_plan_memo_admit(at(64 * 2^20), 0, big))
  own <- b$accounting$owned_entries + b$accounting$borrowed_entries
  expect_true(api$cp04_plan_memo_admit(at(2^30), 0, big))
  expect_false(api$cp04_plan_memo_admit(at(2^30), floor(2^30 / 8) - own - big + 1, big))
  # Missing or non-finite quantities decline the store.
  expect_false(api$cp04_plan_memo_admit(at(2^30), NA_real_, big))
  expect_false(api$cp04_plan_memo_admit(at(2^30), 0, Inf))
  bb <- b; bb$accounting$owned_entries <- NULL
  expect_false(api$cp04_plan_memo_admit(bb, 0, 1))
  # No allocation happens for a declined store: the scope's store stays empty.
  prev <- api$cp04_plan_scope_enter()
  expect_null(api$cp04_plan_memo_get(at(64 * 2^20), seq_len(2251), 1L, 2251L, 0))
  expect_null(api$cp04_plan_scope$current$memo)
  api$cp04_plan_scope_exit(prev)
})

test_that("a store for another key is dropped before a new one is allocated and the certified path keeps none", {
  box <- new.env()
  assign(".store_box", box, envir = globalenv())
  on.exit(rm(".store_box", envir = globalenv()), add = TRUE)
  suppressMessages(trace("cp04_stream_plan", where = api, print = FALSE, tracer = quote(
    if (is.null(.store_box$backend)) { .store_box$backend <- backend; .store_box$links <- links })))
  on.exit(suppressMessages(untrace("cp04_stream_plan", where = api)), add = TRUE)
  store_fit(est_method = "reg", control = list(chunk_rows = 1L))
  b <- box$backend; links <- box$links
  prev <- api$cp04_plan_scope_enter()
  first <- api$cp04_plan_memo_get(b, links, 1L, length(links), 0)
  # The previous store is unbound from the scope before the new one is built.
  suppressMessages(trace("cp04_plan_memo_entries", where = api, print = FALSE,
    tracer = quote(assign(".store_prev", cp04_plan_scope$current$memo, envir = globalenv()))))
  second <- api$cp04_plan_memo_get(b, links[-1L], 1L, length(links) - 1L, 0)
  suppressMessages(untrace("cp04_plan_memo_entries", where = api))
  expect_null(get(".store_prev", envir = globalenv()))
  rm(".store_prev", envir = globalenv())
  expect_false(identical(first, second))
  expect_identical(api$cp04_plan_scope$current$memo$key[[2L]], as.integer(links[-1L]))
  api$cp04_plan_scope_exit(prev)
  # A certified plan after a complete one in the same scope leaves no store.
  box2 <- new.env()
  assign(".store_box", box2, envir = globalenv())
  suppressMessages(trace("cp04_stream_plan", where = api, print = FALSE, tracer = quote(
    if (is.null(.store_box$backend)) { .store_box$backend <- backend; .store_box$links <- links; .store_box$retained <- retained })))
  # (This fit is refused later, by the diagnostic copy; its plans are certified.)
  tryCatch(store_fit(est_method = "dr", outcome_pooling = "separate", control = list(numeric_workspace_bytes = 12 * 2^20)),
    cdid_error = function(e) NULL)
  b2 <- box2$backend; l2 <- box2$links; r2 <- box2$retained
  prev <- api$cp04_plan_scope_enter()
  complete <- api$cp04_stream_plan_legacy(b2, l2, r2)
  expect_gt(length(complete$blocks), 1L)
  expect_false(is.null(api$cp04_plan_scope$current$memo))
  certified <- api$cp04_stream_plan(b2, l2, r2)
  paths <- api$cp04_plan_scope$current$paths
  left <- api$cp04_plan_scope$current$memo
  api$cp04_plan_scope_exit(prev)
  expect_identical(paths[["certified"]], 1)
  expect_null(left)
})
