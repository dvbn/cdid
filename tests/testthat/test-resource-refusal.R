# Resource errors at the public boundary (R/chained_resource_refusal.R): a
# refused call is restated once and promptly, from the typed counts of the
# failing check, with no value stated as sufficient; nothing is evaluated
# again (no fit, pretest or callback re-entry), and the caller's time limit and
# random number state are untouched; refusals that a larger workspace limit
# cannot change keep their own reason; the short-pass warning fires only with
# an automatic chunk size. Public calls, no mocks.

rr_panel <- function() {
  set.seed(41)
  n <- 80L
  d <- expand.grid(id = seq_len(n), t = 1:6)
  d$g <- rep(c(0, 3, 5), length.out = n)[d$id]
  d$x <- stats::rnorm(n)[d$id]
  d$y <- d$x + .3 * d$t + stats::rnorm(nrow(d)) + (d$g > 0 & d$t >= d$g)
  m <- d
  m$y[stats::runif(nrow(m)) < .15] <- NA_real_
  m$xv <- round(m$x + .1 * m$t + stats::rnorm(nrow(m), sd = .3), 3)
  list(d = d, m = m)
}

rr_dr <- function(m, mib, biters = 12000L, bstrap = TRUE, control = list(), ...) {
  att_gt("y", "t", "id", "g", data = m, xformla = ~x, est_method = "dr", outcome_pooling = "separate",
    control_group = "nevertreated", bstrap = bstrap, biters = biters, cband = FALSE,
    control = c(list(numeric_workspace_bytes = mib * 2^20), control), ...)
}

rr_expect_no_sufficient_number <- function(e) {
  testthat::expect_true(is.na(e$resource$sufficient_value))
  testthat::expect_false(isTRUE(e$resource$verified))
  m <- conditionMessage(e)
  testthat::expect_false(grepl("numeric_workspace_bytes = [0-9]", m))
  testthat::expect_false(grepl("admits", m, fixed = TRUE))
}

testthat::test_that("a workspace refusal states the failing check's counts, no sufficient number, and ends with the fallback", {
  p <- rr_panel()
  set.seed(101)
  e <- tryCatch(suppressMessages(rr_dr(p$m, 12)), error = identity)
  testthat::expect_s3_class(e, "cdid_workspace_refusal")
  testthat::expect_s3_class(e, "cdid_resource_refusal")
  testthat::expect_s3_class(e, "cdid_ipw_owner_bound_error")
  testthat::expect_s3_class(e, "cdid_error")
  rr_expect_no_sufficient_number(e)
  r <- e$resource
  testthat::expect_identical(r$setting, "numeric_workspace_bytes")
  testthat::expect_identical(r$limit, 12 * 2^20)
  testthat::expect_identical(r$limit_entries, floor(12 * 2^20 / 8))
  testthat::expect_true(r$requested > r$allowed && r$allowed <= r$limit_entries)
  testthat::expect_identical(r$requested_bytes, 8 * r$requested)
  testthat::expect_true(nzchar(r$check) && nzchar(r$stage) && nzchar(r$reason))
  m <- conditionMessage(e)
  testthat::expect_match(m, "Increase the workspace limit.$")
  testthat::expect_match(m, "control = list(numeric_workspace_bytes = ...)", fixed = TRUE)
  testthat::expect_match(m, "not what the whole call needs", fixed = TRUE)
  testthat::expect_match(m, format(r$requested, big.mark = ",", scientific = FALSE), fixed = TRUE)
  testthat::expect_true(nzchar(e$diagnostics$internal_message))
  # The original diagnostic is kept in the message, as its second sentence.
  testthat::expect_match(m, paste0("Check: ", e$diagnostics$internal_message), fixed = TRUE)
  testthat::expect_identical(r$check, r$stage)
  # The 12,000-draw fixture refuses at 25 MiB and completes with 29 MiB.
  set.seed(101)
  e25 <- tryCatch(suppressMessages(rr_dr(p$m, 25)), error = identity)
  testthat::expect_s3_class(e25, "cdid_workspace_refusal")
  rr_expect_no_sufficient_number(e25)
  set.seed(101)
  fit <- suppressMessages(rr_dr(p$m, 29))
  testthat::expect_true(any(is.finite(fit$att)))
})

testthat::test_that("a refused conditional pretest never runs the test: chained_conditional_test() is entered at most once per call", {
  p <- rr_panel()
  ns <- asNamespace("cdid"); hits <- new.env(); hits$n <- 0L
  suppressMessages(trace("chained_conditional_test", where = ns, print = FALSE,
    tracer = bquote(assign("n", get("n", envir = .(hits)) + 1L, envir = .(hits)))))
  on.exit(suppressMessages(untrace("chained_conditional_test", where = ns)), add = TRUE)
  run <- function(mib) tryCatch(conditional_did_pretest("y", "t", "id", "g", xformla = ~x, data = p$m, est_method = "dr",
    bstrap = TRUE, cband = FALSE, biters = 50L, control = list(numeric_workspace_bytes = mib * 2^20)), error = identity)
  set.seed(10); e <- run(2)
  testthat::expect_s3_class(e, "cdid_workspace_refusal")
  rr_expect_no_sufficient_number(e)
  testthat::expect_identical(hits$n, 0L)
  hits$n <- 0L; set.seed(10); ok <- run(1024)
  testthat::expect_lte(hits$n, 1L)
})

testthat::test_that("a caller's elapsed-time limit still fires after a refused call", {
  p <- rr_panel()
  f <- function() {
    setTimeLimit(elapsed = 3, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
    t0 <- proc.time()[["elapsed"]]
    e <- tryCatch(suppressMessages(rr_dr(p$m, 12, biters = 199L)), error = identity)
    testthat::expect_s3_class(e, "cdid_workspace_refusal")
    tryCatch({ while (proc.time()[["elapsed"]] - t0 < 6) NULL; "completed" }, error = function(err) conditionMessage(err))
  }
  testthat::expect_match(f(), "elapsed time limit", fixed = TRUE)
})

testthat::test_that("a refused call leaves the random number state that 811ed24 leaves", {
  p <- rr_panel()
  # Refused before any draw: the state right after set.seed(101).
  set.seed(101); s0 <- get(".Random.seed", envir = globalenv())
  e <- tryCatch(suppressMessages(rr_dr(p$m, 12)), error = identity)
  testthat::expect_s3_class(e, "cdid_workspace_refusal")
  testthat::expect_identical(get(".Random.seed", envir = globalenv()), s0)
  # A late refusal must leave the exact state captured on 811ed24 after the
  # draws. The first element encodes this R version's RNG kinds.
  if (identical(RNGkind()[1], "Mersenne-Twister")) {
    set.seed(101)
    e25 <- tryCatch(suppressMessages(rr_dr(p$m, 25)), error = identity)
    testthat::expect_s3_class(e25, "cdid_workspace_refusal")
    seed <- get(".Random.seed", envir = globalenv())
    ref <- readRDS(testthat::test_path("fixtures", "r4_late_refusal_811ed24_rng.rds"))
    testthat::expect_identical(seed[1L], s0[1L])
    testthat::expect_identical(seed[-1L], ref)
  }
})

testthat::test_that("a custom estimator is evaluated no more often than on 811ed24 when the call is refused", {
  p <- rr_panel()
  k <- 0L
  custom_ipw <- function(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...) {
    k <<- k + 1L
    dy <- y1 - y0; wt <- i.weights * D; wc <- i.weights * (1 - D)
    mt <- sum(wt * dy) / sum(wt); mc <- sum(wc * dy) / sum(wc)
    list(ATT = mt - mc, att.inf.func = length(dy) * (wt * (dy - mt) / sum(wt) - wc * (dy - mc) / sum(wc)))
  }
  # Evaluations on 811ed24 at these limits: 0, 0 and 30 (each call refused).
  for (case in list(c(0.5, 0), c(1, 0), c(2, 30))) {
    k <- 0L
    e <- tryCatch(att_gt("y", "t", "id", "g", data = p$d, est_method = custom_ipw, control_group = "nevertreated",
      bstrap = FALSE, cband = FALSE, control = list(numeric_workspace_bytes = case[1] * 2^20)), error = identity)
    # Refused at every limit; restated only where the failing site records its control.
    testthat::expect_s3_class(e, "cdid_error")
    testthat::expect_identical(k, as.integer(case[2]))
  }
})

testthat::test_that("work-limit refusals and refusals without counts keep their own reason at any workspace limit", {
  p <- rr_panel()
  for (mib in c(1024, 4096)) {
    e <- tryCatch(rr_dr(p$m, mib, bstrap = FALSE, control = list(counted_work_limit = 3e8)), error = identity)
    testthat::expect_s3_class(e, "cdid_work_limit_refusal"); testthat::expect_s3_class(e, "cdid_resource_refusal")
    testthat::expect_false(inherits(e, "cdid_workspace_refusal"))
    testthat::expect_identical(e$resource$setting, "counted_work_limit")
    testthat::expect_match(conditionMessage(e), paste0("Check: ", e$diagnostics$internal_message), fixed = TRUE)
    testthat::expect_match(conditionMessage(e), "The workspace limit does not change it", fixed = TRUE)
    rr_expect_no_sufficient_number(e)
    # The catalog work check records its counts and its control.
    c1 <- tryCatch(rr_dr(p$m, mib, bstrap = FALSE, control = list(counted_work_limit = 1e6)), error = identity)
    testthat::expect_s3_class(c1, "cdid_catalog_bound_error"); testthat::expect_s3_class(c1, "cdid_work_limit_refusal")
    testthat::expect_match(conditionMessage(c1), "Catalog adapter output/string/graph work exceeds its explicit stored controls", fixed = TRUE)
    testthat::expect_match(conditionMessage(c1), "counted_work_limit", fixed = TRUE)
  }
})

testthat::test_that("a check against a remainder of the work limit is not called a workspace refusal (PLAN's B1)", {
  p <- rr_panel()
  ipw <- function(control) tryCatch(att_gt("y", "t", "id", "g", data = p$m, xformla = ~x, est_method = "ipw",
    control_group = "nevertreated", bstrap = FALSE, cband = FALSE, control = control), error = identity)
  a <- ipw(list(counted_work_limit = 30738983))
  b <- ipw(list(counted_work_limit = 30738983, numeric_workspace_bytes = 2^36))
  testthat::expect_s3_class(a, "cdid_error")
  testthat::expect_false(inherits(a, "cdid_workspace_refusal"))
  # The site records no governing control (its bound is a remainder), so the condition passes as on 811ed24,
  # and the same at 64 GiB.
  testthat::expect_identical(conditionMessage(a), "Total compact payload traversal gate (18560) exceeds explicit bound (18267).")
  testthat::expect_identical(class(b), class(a)); testthat::expect_identical(conditionMessage(b), conditionMessage(a))
  testthat::expect_true(length(ipw(list(counted_work_limit = 4e7))$att) > 0L)
})

testthat::test_that("the governing control is read from the bound's expression at the signal site, never from its size", {
  k <- function(x, caller = NULL) resource_bound_kind(x, caller)
  testthat::expect_identical(k(quote(ctrl$owner_control$max_matrix_entries)), "workspace")
  testthat::expect_identical(k(quote(min(ctrl$max_retained_entries, ctrl$max_matrix_entries))), "workspace")
  testthat::expect_identical(k(quote(ctrl$max_work_ops)), "work")
  testthat::expect_true(is.na(k(quote(min(ctrl$max_work_ops, ctrl$max_matrix_entries)))))
  testthat::expect_true(is.na(k(quote(cap_wo - prior_work - work))))
  testthat::expect_true(is.na(k(quote(bounds$max_materialized_if_entries))))
  testthat::expect_true(is.na(k(quote(gmm_bounds$max_targets))))
  testthat::expect_true(is.na(k(134217728)))
  f <- function(ctrl, cap) { max_me <- ctrl$max_matrix_entries; rem <- ctrl$max_work_ops - 1; NULL }
  testthat::expect_identical(k(quote(max_me), f), "workspace")
  testthat::expect_true(is.na(k(quote(rem), f)))
  testthat::expect_true(is.na(k(quote(cap), f)))
  # A condition without a recorded control passes unchanged whatever its numbers.
  e0 <- structure(class = c("cdid_ipw_owner_bound_error", "cdid_error", "error", "condition"),
    list(message = "Total x (300) exceeds explicit bound (200).", call = NULL, diagnostics = list(label = "x", requested = 300, bound = 200)))
  testthat::expect_identical(resource_restate(e0, resource_call_limits(list(numeric_workspace_bytes = 1600))), e0)
})

testthat::test_that("a fixed size limit says that no setting changes it, also at limits above the integer range", {
  e0 <- structure(class = c("cdid_ipw_owner_matrix_dim_error", "cdid_ipw_owner_bound_error", "cdid_error", "error", "condition"),
    list(message = "Real 2D matrix dimension for x (3e+09 x 2) exceeds R 32-bit matrix dimension limit (`INT_MAX` = 2147483647).",
      call = NULL, diagnostics = list(label = "x", nrow = 3e9, ncol = 2, matrix_dim_limit = 2147483647)))
  for (bytes in c(2^30, 2^40)) {
    e <- resource_restate(e0, resource_call_limits(list(numeric_workspace_bytes = bytes)))
    testthat::expect_s3_class(e, "cdid_size_limit_refusal"); testthat::expect_s3_class(e, "cdid_ipw_owner_matrix_dim_error")
    testthat::expect_match(conditionMessage(e), paste0("Check: ", conditionMessage(e0)), fixed = TRUE)
    testthat::expect_match(conditionMessage(e), "a size limit that no setting changes", fixed = TRUE)
    testthat::expect_true(is.na(e$resource$sufficient_value))
  }
  # Typed counts whose limit is neither the workspace limit nor the work limit pass unchanged.
  big <- structure(class = c("cdid_ipw_owner_bound_error", "cdid_error", "error", "condition"),
    list(message = "Total x (300) exceeds explicit bound (200).", call = NULL, diagnostics = list(label = "x", requested = 300, bound = 200)))
  testthat::expect_identical(resource_restate(big, list(bytes = 800, entries = 100, work = Inf)), big)
  plain <- simpleError("boom")
  testthat::expect_identical(resource_restate(plain, resource_call_limits(NULL)), plain)
})

testthat::test_that("the restatement leaves admitted fits and other errors unchanged", {
  p <- rr_panel()
  a <- rr_dr(p$m, 1024, bstrap = FALSE)
  b <- att_gt("y", "t", "id", "g", data = p$m, xformla = ~x, est_method = "dr", outcome_pooling = "separate",
    control_group = "nevertreated", bstrap = FALSE, cband = FALSE)
  testthat::expect_identical(a$att, b$att); testthat::expect_identical(a$se, b$se)
  e <- tryCatch(rr_dr(p$m, -1 / 2^20, bstrap = FALSE), error = identity)
  testthat::expect_s3_class(e, "cdid_chained_public_control_error")
  testthat::expect_false(inherits(e, "cdid_resource_refusal"))
  testthat::expect_null(e$resource)
})

testthat::test_that("an aggregation refused by the fit's workspace limit says that the fit must be rerun", {
  p <- rr_panel()
  fit <- att_gt("y", "t", "id", "g", data = p$m, xformla = ~x, est_method = "dr", outcome_pooling = "separate",
    control_group = "nevertreated", allow_unbalanced_panel = TRUE, bstrap = FALSE, cband = FALSE)
  e <- tryCatch(aggte(fit, type = "simple", bstrap = TRUE, cband = FALSE, biters = 5e7), error = identity)
  testthat::expect_s3_class(e, "cdid_workspace_refusal"); testthat::expect_s3_class(e, "cdid_ipw_pipeline_bound_error")
  rr_expect_no_sufficient_number(e)
  testthat::expect_identical(e$resource$limit, 2^30); testthat::expect_identical(e$resource$limit_from, "fit")
  m <- conditionMessage(e)
  testthat::expect_match(m, "The workspace limit comes from the fit, so the fit must be rerun with a larger limit", fixed = TRUE)
  testthat::expect_match(m, "Increase the workspace limit.$")
  testthat::expect_match(m, "Check: Simultaneous IF/sign/interval/bootstrap producer and consumer payload exceeds caller limits.", fixed = TRUE)
  # A work-limit refusal inside an aggregation (no public call reaches one cheaply): restated from its typed record
  # with the fit's work limit, not as a workspace refusal.
  w0 <- structure(class = c("cdid_ipw_pipeline_bound_error", "cdid_error", "error", "condition"),
    list(message = "Total x work (500) exceeds explicit bound (400).", call = NULL,
      diagnostics = list(label = "x work", requested = 500, bound = 400, limit_kind = "work")))
  mp <- list(chained = list(aggregation_state = list(public = list(bytes = 2^30, resolved = list(profile = list(counted_work_limit = 400))))))
  w <- tryCatch(resource_guard_fitted(stop(w0), mp), error = identity)
  testthat::expect_s3_class(w, "cdid_work_limit_refusal"); testthat::expect_false(inherits(w, "cdid_workspace_refusal"))
  testthat::expect_identical(w$resource$limit, 400); testthat::expect_identical(w$resource$limit_from, "fit")
  testthat::expect_match(conditionMessage(w), "The limit comes from the fit", fixed = TRUE)
})

testthat::test_that("internal reservation failures are internal errors, not model-fit errors", {
  e <- tryCatch(cp04_abort("Retained model bound underestimated.", "bound"), error = identity)
  testthat::expect_s3_class(e, "cdid_cp04_bound_error"); testthat::expect_s3_class(e, "cdid_internal_error")
  testthat::expect_false(inherits(e, "cdid_ipw_fit_error"))
  f <- tryCatch(cp04_abort("Logit encountered boundary probabilities.", "boundary"), error = identity)
  testthat::expect_s3_class(f, "cdid_ipw_fit_error")
})

testthat::test_that("a switching-model refusal (dname) is restated without a sufficient number", {
  set.seed(3); n <- 60L; tt <- 4L
  d <- expand.grid(id = seq_len(n), time = seq_len(tt))
  sw <- rep(c(0L, 3L, 4L), length.out = n)
  d$d <- as.integer(sw[d$id] > 0 & d$time >= sw[d$id]); d$x1 <- stats::rnorm(n)[d$id]
  d$y <- d$x1 + 0.2 * d$time + d$d + stats::rnorm(nrow(d))
  e <- tryCatch(att_gt(yname = "y", tname = "time", idname = "id", dname = "d", xformla = ~x1, data = d,
    bstrap = FALSE, cband = FALSE, control = list(numeric_workspace_bytes = 4096)), error = identity)
  testthat::expect_s3_class(e, "cdid_error")
  if (inherits(e, "cdid_resource_refusal")) rr_expect_no_sufficient_number(e)
})

testthat::test_that("the short-pass warning fires with an automatic chunk size only, once, and states no workspace value", {
  p <- rr_panel()
  run <- function(mib, chunk = NULL) {
    ws <- list()
    ctl <- c(list(numeric_workspace_bytes = mib * 2^20), if (!is.null(chunk)) list(chunk_rows = chunk))
    withCallingHandlers(tryCatch(att_gt("y", "t", "id", "g", data = p$m, xformla = ~xv, est_method = "reg",
      selection = "conditional_mtar", selection_est_method = "dr", selection_xformla = ~1, gmm_weighting = "within",
      bstrap = FALSE, cband = FALSE, control = ctl), error = identity),
      warning = function(w) { if (inherits(w, "cdid_workspace_speed_warning")) ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning") })
    ws
  }
  # An explicit chunk of eight rows at 128 MiB cannot fit.
  testthat::expect_length(run(128, 8L), 0L)
  # Automatic rows on the dated augmented path: the first small limit at which the warning fires (the row count
  # depends on the accounting of the merged branches, so it is not pinned).
  w <- NULL; at <- NA
  for (mib in c(16, 20, 24)) { w <- run(mib); if (length(w)) { at <- mib; break } }
  testthat::expect_length(w, 1L)
  testthat::expect_true(w[[1]]$rows < 16)
  testthat::expect_true(is.na(w[[1]]$suggested_value))
  testthat::expect_match(conditionMessage(w[[1]]), "may take a long time", fixed = TRUE)
  testthat::expect_false(grepl("numeric_workspace_bytes = [0-9]|GiB|MiB", conditionMessage(w[[1]])))
  # The same rows set explicitly at the same limit: no warning.
  testthat::expect_length(run(at, as.integer(w[[1]]$rows)), 0L)
})

testthat::test_that("a work record passes unchanged when the work limit is off", {
  w0 <- structure(class = c("cdid_ipw_owner_bound_error", "cdid_error", "error", "condition"),
    list(message = "Total x work (500) exceeds explicit bound (400).", call = NULL,
      diagnostics = list(label = "x work", requested = 500, bound = 400, limit_kind = "work")))
  testthat::expect_identical(resource_restate(w0, resource_call_limits(NULL)), w0)
  testthat::expect_identical(resource_restate(w0, resource_call_limits(list(counted_work_limit = Inf))), w0)
  testthat::expect_s3_class(resource_restate(w0, resource_call_limits(list(counted_work_limit = 400))), "cdid_work_limit_refusal")
})

testthat::test_that("the alias resolver gives no kind where it cannot trust an assignment (ACCT's probe cases)", {
  k <- function(x, caller) resource_bound_kind(x, caller)
  f1 <- function(ctrl) { a <- ctrl$max_matrix_entries; b <- ctrl$max_retained_entries; cap <- min(a, b); NULL }
  testthat::expect_true(is.na(k(quote(cap), f1)))
  f2 <- function(ctrl, flag) { cap <- if (flag) NULL else min(ctrl$max_retained_entries, ctrl$max_matrix_entries); NULL }
  # An `if` with a NULL branch has the kind of its other branch (the check is not applied on NULL).
  testthat::expect_identical(k(quote(cap), f2), "workspace")
  f3 <- compiler::cmpfun(function(ctrl) { cap <- ctrl$max_matrix_entries; NULL })
  testthat::expect_identical(k(quote(cap), f3), "workspace")
  testthat::expect_true(is.na(k(quote(cap), sum)))
  f5 <- function(ctrl) { cap <- ctrl$max_matrix_entries; g <- function() cap <<- 5; g(); NULL }
  testthat::expect_true(is.na(k(quote(cap), f5)))
  f6 <- function(ctrl, n) { cap <- ctrl$max_matrix_entries; cap[1] <- n; NULL }
  testthat::expect_true(is.na(k(quote(cap), f6)))
  f7 <- function(ctrl) { g <- function() { cap <- ctrl$max_matrix_entries; cap }; NULL }
  testthat::expect_true(is.na(k(quote(cap), f7)))
  testthat::expect_true(is.na(k(quote(x$max_matrix_entries), NULL)))
  testthat::expect_true(is.na(k(quote(dims$max_work_ops), NULL)))
  testthat::expect_identical(k(quote(backend$control$max_matrix_entries), NULL), "workspace")
  testthat::expect_identical(k(quote(own_ctrl$max_work_ops), NULL), "work")
  f9 <- function(ctrl) { cap <- ctrl$max_matrix_entries; assign("cap", 3); NULL }
  testthat::expect_true(is.na(k(quote(cap), f9)))
  # A NULL branch means the check is not applied: `if (cond) NULL else X` has the kind of X.
  testthat::expect_identical(k(quote(if (memory_phase_route) NULL else min(ctrl$owner_control$max_retained_entries,
    ctrl$owner_control$max_matrix_entries)), NULL), "workspace")
  testthat::expect_identical(k(quote(if (flag) ctrl$max_work_ops else NULL), NULL), "work")
  testthat::expect_true(is.na(k(quote(if (flag) ctrl$max_work_ops else ctrl$max_matrix_entries), NULL)))
  testthat::expect_true(is.na(k(quote(if (flag) NULL else cap - used), NULL)))
})

testthat::test_that("every value assigned to a typed field name in the package is one reviewed for that control", {
  # Invariant behind the expression match (ACCT's review): a field named after a workspace or work control
  # holds that control, a min() of such fields, or a remainder of the same control. A new right-hand side
  # fails here until it is reviewed and added.
  allowed <- c(
    "max_matrix_entries <- rem_me",
    "max_matrix_entries <- solve_entries",
    "max_matrix_entries = ",
    "max_matrix_entries = NULL",
    "max_matrix_entries = bnds$max_matrix_entries",
    "max_matrix_entries = bounds$max_matrix_entries",
    "max_matrix_entries = control$max_matrix_entries",
    "max_matrix_entries = control$owner_control$max_matrix_entries",
    "max_matrix_entries = ctrl$max_matrix_entries",
    "max_matrix_entries = entries",
    "max_matrix_entries = gmm_int_bound(max_matrix_entries, \"max_matrix_entries\")",
    "max_matrix_entries = m_e",
    "max_matrix_entries = max_e",
    "max_matrix_entries = max_entries",
    "max_matrix_entries = max_lim",
    "max_matrix_entries = max_matrix_entries",
    "max_matrix_entries = max_me",
    "max_matrix_entries = max_ret",
    "max_matrix_entries = min(ctrl$max_matrix_entries, ctrl$max_output_entries)",
    "max_matrix_entries = min(max_me, max_out)",
    "max_matrix_entries = min(own_ctrl$max_matrix_entries, pt_ctrl$max_matrix_entries)",
    "max_matrix_entries = mme",
    "max_matrix_entries = ws_entries",
    "max_output_entries = ",
    "max_output_entries = control$max_output_entries",
    "max_output_entries = ctrl$max_output_entries",
    "max_output_entries = entries",
    "max_output_entries = max_lim",
    "max_output_entries = max_out",
    "max_output_entries = min(own_ctrl$max_output_entries, pt_ctrl$max_output_entries)",
    "max_retained_entries = ",
    "max_retained_entries = control$max_retained_entries",
    "max_retained_entries = ctrl$max_retained_entries",
    "max_retained_entries = entries",
    "max_retained_entries = max_ret",
    "max_retained_entries = min(own_ctrl$max_retained_entries, pt_ctrl$max_retained_entries)",
    "max_work = ",
    "max_work = max(0, ctrl$max_work_ops - pre$total_work_ops_bound - pre$character_metadata$work_ops_bound -      128)",
    "max_work = work",
    "max_work_ops <- rem_wo",
    "max_work_ops = ",
    "max_work_ops = control$max_work_ops",
    "max_work_ops = ctrl$max_work_ops",
    "max_work_ops = max_wo",
    "max_work_ops = min(own_ctrl$max_work_ops, pt_ctrl$max_work_ops)",
    "max_work_ops = work",
    "max_workspace_cells = as.double(workspace_entries)",
    "max_workspace_cells = entries"
  )
  ns <- asNamespace("cdid")
  fields <- c(resource_workspace_fields, resource_work_fields)
  rhs <- character()
  walk <- function(x) {
    if (is.call(x)) {
      nm <- names(x)
      if (!is.null(nm)) for (i in which(nm %in% fields)) rhs <<- c(rhs, paste0(nm[i], " = ", paste(deparse(x[[i]]), collapse = " ")))
      if ((identical(x[[1L]], as.name("<-")) || identical(x[[1L]], as.name("="))) && length(x) == 3L && is.call(x[[2L]]) &&
          identical(x[[2L]][[1L]], as.name("$")) && as.character(x[[2L]][[3L]]) %in% fields)
        rhs <<- c(rhs, paste0(as.character(x[[2L]][[3L]]), " <- ", paste(deparse(x[[3L]]), collapse = " ")))
      for (el in as.list(x)) if (!missing(el)) walk(el)
    }
  }
  for (f in ls(ns, all.names = TRUE)) {
    g <- get(f, ns)
    if (is.function(g) && !is.primitive(g)) { walk(body(g)); walk(as.call(c(as.name("list"), formals(g)))) }
  }
  testthat::expect_identical(setdiff(unique(rhs), allowed), character(0))
})

testthat::test_that("the pipeline prefix check is typed at its site and restated as a workspace refusal", {
  p <- rr_panel()
  e <- tryCatch(att_gt("y", "t", "id", "g", data = p$m, xformla = ~x, est_method = "dr", control_group = "nevertreated",
    gmm_weighting = "identity", bstrap = FALSE, cband = FALSE, control = list(numeric_workspace_bytes = 2^20)), error = identity)
  testthat::expect_s3_class(e, "cdid_workspace_refusal")
  testthat::expect_identical(e$resource$check, "pipeline simultaneous owners, nuisances, ledger/basis, retained results and inference matrices")
  rr_expect_no_sufficient_number(e)
})

testthat::test_that("cp04_budget() is typed per line: payload lines on the workspace limit, the work charge on the work limit", {
  ctrl <- list(max_retained_entries = 100, max_matrix_entries = 100, max_work_ops = 10, max_output_entries = 100)
  kind <- function(...) tryCatch(cp04_budget(ctrl, ...), error = function(e) e$diagnostics$limit_kind)
  testthat::expect_identical(kind(150, 0, 0, 0, "x"), "workspace")
  testthat::expect_identical(kind(50, 60, 0, 0, "x"), "workspace")
  testthat::expect_identical(kind(1, 1, 50, 0, "x"), "work")
  testthat::expect_identical(kind(1, 1, 1, 150, "x"), "workspace")
})

testthat::test_that("the count helpers signal the original condition when the refusal file is not loaded with them", {
  # The owner helpers copied into an environment that holds only them (as when their file is
  # sourced alone): the optional record is absent and the condition is the original one, as on 811ed24.
  ns <- asNamespace("cdid")
  env <- new.env(parent = baseenv())
  for (nm in ls(ns, pattern = "^(ipw_owner_|cdid_ipw_owner_)")) {
    if (!is.function(get(nm, ns))) next
    f <- get(nm, ns); environment(f) <- env; assign(nm, f, envir = env)
  }
  run <- function(env) list(
    prod = tryCatch(env$ipw_owner_safe_product(10, 10, max_entries = 50, label = "x"), error = identity),
    sum = tryCatch(env$ipw_owner_safe_sum(c(30, 40), max_bound = 50, label = "y"), error = identity))
  cls <- c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error", "cdid_ipw_owner_refusal", "cdid_ipw_owner_error",
    "cdid_ipw_enclosure_error", "cdid_ipw_error", "cdid_error", "error", "condition")
  check <- function(r) {
    testthat::expect_identical(class(r$prod), cls); testthat::expect_identical(class(r$sum), cls)
    testthat::expect_identical(conditionMessage(r$prod), "Count for x (100) exceeds explicit bound (50).")
    testthat::expect_identical(conditionMessage(r$sum), "Total y (70) exceeds explicit bound (50).")
    testthat::expect_true(is.na(r$prod$diagnostics$limit_kind)); testthat::expect_true(is.na(r$sum$diagnostics$limit_kind))
  }
  check(run(env))
  # The same with the file itself sourced alone, when the package sources are present (source-tree runs).
  src <- testthat::test_path("..", "..", "R", "chained_ipw_owner.R")
  if (file.exists(src)) { e2 <- new.env(parent = baseenv()); sys.source(src, envir = e2); check(run(e2)) }
  # The helpers keep their 811ed24 formals.
  testthat::expect_identical(names(formals(ns$ipw_owner_safe_sum)), c("vals", "max_bound", "label"))
  testthat::expect_identical(names(formals(ns$ipw_owner_safe_product)), c("a", "b", "max_entries", "label"))
  testthat::expect_identical(names(formals(ns$gmm_safe_product)), c("a", "b", "max_entries", "label"))
  # With the resolver available, a typed expression still resolves.
  f <- function(ctrl) ns$ipw_owner_safe_sum(c(30, 40), ctrl$max_matrix_entries, "y")
  testthat::expect_identical(tryCatch(f(list(max_matrix_entries = 50)), error = function(e) e$diagnostics$limit_kind), "workspace")
})
