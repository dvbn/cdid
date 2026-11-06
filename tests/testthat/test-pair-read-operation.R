v49_ns <- function(name) getFromNamespace(name, "cdid")
v49_fixture <- local({
  e <- new.env(parent = asNamespace("cdid"))
  e$cp04_reference <- new.env(parent = asNamespace("stats"))
  invisible(capture.output(sys.source(testthat::test_path("reference", "cp04_finite_reference.R"), e$cp04_reference)))
  wanted <- c("cp04_ns", "cp04_control", "cp04_glm_control", "cp04_fixture")
  for (x in parse(file = testthat::test_path("test-cp04-installed.R"))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE)))
  e
})
v49_trace <- function(run, enabled = TRUE) {
  issue <- v49_ns("cp04_pair_read_operation"); operations <- list()
  testthat::local_mocked_bindings(cp04_pair_read_operation = function(...) {
    op <- if (enabled) issue(...) else NULL
    if (!is.null(op)) operations[[length(operations) + 1L]] <<- op
    op
  }, .package = "cdid")
  result <- run()
  list(result = result, operations = operations,
    counts = lapply(operations, function(op) op$counts()))
}
v49_reowner <- function(fx, panel_change = identity, a_change = identity) {
  panel <- fx$panel; panel["ref"] <- list(NULL); panel <- panel_change(panel)
  panel <- v49_ns("ipw_owner_seal")(panel, "panel_owner")
  frame <- fx$frame; frame["ref"] <- list(NULL); frame$panel_ref <- panel$ref
  frame <- v49_ns("ipw_owner_seal")(frame, "frame_owner")
  reader <- fx$reader; reader["ref"] <- list(NULL)
  reader$panel_ref <- panel$ref; reader$frame_ref <- frame$ref
  reader$state$a_index <- a_change(reader$state$a_index)
  reader <- v49_ns("ipw_owner_seal")(reader, "cp04_reader")
  list(panel = panel, frame = frame, reader = reader, views = fx$views)
}
v49_rebind <- function(env, name, value) {
  old <- get(name, env, inherits = FALSE)
  unlockBinding(name, env); assign(name, value, env); lockBinding(name, env)
  function() { unlockBinding(name, env); assign(name, old, env); lockBinding(name, env) }
}
v49_error <- function(run) tryCatch(run(), error = identity)

testthat::test_that("lazy positions match original pair order and avoid repeated sparse searches", {
  fx <- v49_fixture$cp04_fixture(); ref <- v49_fixture$cp04_reference
  rows <- unique(as.integer(c(which(ref$A == 0)[1:5], 1:12, (ref$N - 5L):ref$N)))
  rows <- rev(rows)
  op <- v49_ns("cp04_pair_read_operation")(fx$reader, rows, 1e7)
  testthat::expect_type(op, "list")
  pairs <- lapply(list(c(1, 2), c(2, 4), c(1, 4)), function(d)
    v49_ns("ipw_owner_make_pair_descriptor")(fx$panel, d[1], d[2]))
  for (read in c(FALSE, TRUE)) for (pair in pairs) for (query in list(rows, rev(rows), rows[c(3L, 1L, 3L)]))
    testthat::expect_identical(op$pair(query, pair, read), v49_ns("cp04_pair_rows")(fx$reader, query, pair, read))
  before <- op$counts()
  for (pair in pairs) op$pair(rev(rows), pair, TRUE)
  after <- op$counts()
  testthat::expect_identical(unname(after["searches"]), unname(before["searches"]))
  testthat::expect_gt(after["cache_hits"], before["cache_hits"])
  testthat::expect_lte(after["searches"], length(rows) * length(fx$panel$all_times))
  testthat::expect_equal(after[["fallbacks"]], 0)
  op$finish(); op$abort()
  testthat::expect_identical(op$counts(), after)
  testthat::expect_error(op$pair(rows, pairs[[1L]]), class = "cdid_cp04_integrity_error")
  testthat::expect_silent(op$abort())
  # Complete reservation has an exact capacity-only optional boundary.
  B <- length(rows); T <- length(fx$panel$postings$date_start)
  required <- 4096 + 64 * B + 8 * T + 3 * B * T
  testthat::expect_null(v49_ns("cp04_pair_read_operation")(fx$reader, rows, required - 1))
  exact <- v49_ns("cp04_pair_read_operation")(fx$reader, rows, required)
  testthat::expect_equal(exact$reservation, required)
  testthat::expect_identical(exact$pair(rows, pairs[[1L]], TRUE), v49_ns("cp04_pair_rows")(fx$reader, rows, pairs[[1L]], TRUE))
  exact$finish(); exact$abort()
  # Row indices alone may be package compact sequences; source arrays stay borrowed.
  sequence <- v49_ns("cp04_pair_read_operation")(fx$reader, seq.int(1L, 23L), 1e7)
  testthat::expect_type(sequence, "list")
  testthat::expect_identical(sequence$pair(c(23L, 1L, 5L), pairs[[1L]], TRUE),
    v49_ns("cp04_pair_rows")(fx$reader, c(23L, 1L, 5L), pairs[[1L]], TRUE))
  sequence$finish(); sequence$abort()
})

testthat::test_that("endpoint, difference, visited-posting and queried-A refusals retain original behavior", {
  fx <- v49_fixture$cp04_fixture(); rows <- integer(23L); rows[] <- 1:23; pair <- fx$views[[1L]]$pair
  changes <- list(
    endpoint = function(p) { p$observations$y[] <- Inf; p },
    difference = function(p) {
      p$observations$y[p$observations$time == 1] <- -.Machine$double.xmax
      p$observations$y[p$observations$time == 2] <- .Machine$double.xmax; p
    },
    positions = function(p) { p$postings$post_obs_index[] <- NA_integer_; p },
    posting = function(p) {
      j <- p$postings$date_start[1L] + (p$postings$date_end[1L] - p$postings$date_start[1L]) %/% 2L
      p$postings$post_a_index[j] <- NA_integer_; p
    })
  for (name in names(changes)) {
    bad <- v49_reowner(fx, changes[[name]])
    op <- v49_ns("cp04_pair_read_operation")(bad$reader, rows, 1e7)
    testthat::expect_type(op, "list")
    if (name != "posting") testthat::expect_identical(op$pair(rows, pair, FALSE),
      v49_ns("cp04_pair_rows")(bad$reader, rows, pair, FALSE))
    actual <- v49_error(function() op$pair(rows, pair, TRUE))
    original <- v49_error(function() v49_ns("cp04_pair_rows")(bad$reader, rows, pair, TRUE))
    testthat::expect_s3_class(actual, "error")
    testthat::expect_identical(class(actual), class(original))
    testthat::expect_identical(conditionMessage(actual), conditionMessage(original))
    op$abort()
  }
  # Original native_ok declines these queried A values, leaving its scalar
  # keep=which(ai>0) semantics. The new operation must decline, not throw.
  for (value in list(NA_integer_, -1L)) {
    bad <- v49_reowner(fx, a_change = function(a) { a[2L] <- value; a })
    op <- v49_ns("cp04_pair_read_operation")(bad$reader, rows, 1e7)
    for (read in c(FALSE, TRUE)) testthat::expect_identical(op$pair(rows, pair, read),
      v49_ns("cp04_pair_rows")(bad$reader, rows, pair, read))
    testthat::expect_equal(op$counts()[["fallbacks"]], 2)
    testthat::expect_equal(op$counts()[["searches"]], 0)
    op$finish(); op$abort()
  }
})

testthat::test_that("unsupported sources/helpers and changed operation boundaries keep strict ownership", {
  fx <- v49_fixture$cp04_fixture(); rows <- integer(23L); rows[] <- 1:23; pair <- fx$views[[1L]]$pair
  decorated <- v49_reowner(fx, function(p) { attr(p$observations$y, "v49_metadata") <- TRUE; p })
  testthat::expect_null(v49_ns("cp04_pair_read_operation")(decorated$reader, rows, 1e7))
  compact <- v49_reowner(fx, function(p) { p$postings$post_obs_index <- seq_len(length(p$postings$post_obs_index)); p })
  testthat::expect_null(v49_ns("cp04_pair_read_operation")(compact$reader, rows, 1e7))
  helper_case <- function() {
    original <- v49_ns("ipw_owner_read_y_positions")
    testthat::local_mocked_bindings(ipw_owner_read_y_positions = function(...) original(...), .package = "cdid")
    testthat::expect_null(v49_ns("cp04_pair_read_operation")(fx$reader, rows, 1e7))
  }
  helper_case()
  op <- v49_ns("cp04_pair_read_operation")(fx$reader, rows, 1e7)
  testthat::expect_identical(op$pair(structure(c(2L, 1L), names = c("second", "first")), pair, TRUE),
    v49_ns("cp04_pair_rows")(fx$reader, structure(c(2L, 1L), names = c("second", "first")), pair, TRUE))
  testthat::expect_identical(op$pair(c(1L, 24L), pair, TRUE), v49_ns("cp04_pair_rows")(fx$reader, c(1L, 24L), pair, TRUE))
  testthat::expect_equal(op$counts()[["fallbacks"]], 2)
  testthat::expect_identical(op$pair(c(1L, 2L, 3L), pair, FALSE),
    v49_ns("cp04_pair_rows")(fx$reader, c(1L, 2L, 3L), pair, FALSE))
  testthat::expect_gt(op$counts()[["queries"]], 0)
  changed <- fx$panel; changed$observations$y[1L] <- changed$observations$y[1L] + 1
  restore <- v49_rebind(fx$panel$ref, "panel_owner", changed)
  error <- v49_error(op$finish); restore(); op$abort()
  testthat::expect_s3_class(error, "cdid_cp04_integrity_error")
  dll_failure <- function() {
    op <- v49_ns("cp04_pair_read_operation")(fx$reader, rows, 1e7)
    on.exit(op$abort(), add = TRUE)
    testthat::local_mocked_bindings(gmm_require_native_dll = function()
      stop("v49 live DLL boundary failure"), .package = "cdid")
    testthat::expect_error(op$finish(), "v49 live DLL boundary failure", fixed = TRUE)
  }
  dll_failure()
})

testthat::test_that("legacy CP04 ordinary replay activates cached positions with exact independent full-F IF", {
  fx <- v49_fixture$cp04_fixture(); ref <- v49_fixture$cp04_reference
  fit <- v49_ns("fit_cp04_corollary1")(fx$reader, fx$views, v49_fixture$cp04_control, v49_fixture$cp04_glm_control)
  N <- fx$reader$state$N; actual <- old <- matrix(0, N, 4L)
  run <- function(enabled) v49_trace(function() {
    result <- matrix(0, N, 4L)
    for (start in seq.int(1L, N, by = 23L)) {
      rows <- rev(seq.int(start, min(N, start + 22L)))
      out <- v49_ns("evaluate_cp04_if_chunk")(fit, rows, c(3L, 1L, 4L, 2L))
      result[rows, ] <- out$phi
    }
    result
  }, enabled)
  fast <- run(TRUE); slow <- run(FALSE)
  testthat::expect_gt(length(fast$operations), 0L)
  testthat::expect_identical(fast$result, slow$result)
  testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["queries"]), numeric(1L))), 0)
  testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["cache_hits"]), numeric(1L))), 0)
  for (j in seq_along(fast$operations)) testthat::expect_lte(fast$counts[[j]]["searches"],
    fast$operations[[j]]$chunk_rows * fast$operations[[j]]$n_dates)
  independent <- ref$last[, 22:24] - ref$last[, 25:27]
  testthat::expect_lt(max(abs(fast$result[, c(2L, 4L, 1L)] - independent)), ref$tol)
  testthat::expect_identical(crossprod(fast$result), crossprod(slow$result))
  testthat::expect_identical(sqrt(diag(crossprod(fast$result))) / N, sqrt(diag(crossprod(slow$result))) / N)
  testthat::expect_gt(max(abs(fast$result[ref$A == 0, , drop = FALSE])), .1)
  pre <- v49_ns("cp04_if_budget")(fit, 23L, 1:4)
  capped <- fit; capped["ref"] <- list(NULL)
  capped$control$max_matrix_entries <- sum(capped$accounting$owned_entries, capped$accounting$borrowed_entries) + pre$scratch
  capped <- v49_ns("ipw_owner_seal")(capped, "cp04_backend")
  testthat::expect_identical(v49_trace(function() v49_ns("evaluate_cp04_if_chunk")(capped, 1:23, 1:4))$result,
    v49_trace(function() v49_ns("evaluate_cp04_if_chunk")(capped, 1:23, 1:4), FALSE)$result)
  low <- capped; low["ref"] <- list(NULL); low$control$max_matrix_entries <- low$control$max_matrix_entries - 1
  low <- v49_ns("ipw_owner_seal")(low, "cp04_backend")
  testthat::expect_error(v49_ns("evaluate_cp04_if_chunk")(low, 1:23, 1:4), class = "cdid_ipw_owner_bound_error")
})


testthat::test_that("actual public weighted IPW activates the operation and keeps full outputs and statuses", {
  e <- new.env(parent = environment())
  e$sp_raw <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_baseline", "raw.csv"))
  e$sp_raw$cohort <- as.numeric(replace(e$sp_raw$cohort, e$sp_raw$cohort == "NT", "Inf"))
  e$sp_trans <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_corollary1", "transport.csv"))
  e$sp_trans$cohort <- as.numeric(replace(e$sp_trans$cohort, e$sp_trans$cohort == "NT", "Inf"))
  for (x in parse(file = testthat::test_path("test-selection-public.R"))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name("sp_public_fixture"))) eval(x, e)
  args <- e$sp_public_fixture()$args
  args$selection_est_method <- "ipw"
  args$control <- list(chunk_rows = 23L)
  fast <- v49_trace(function() do.call(cdid::att_gt, args))
  slow <- v49_trace(function() do.call(cdid::att_gt, args), FALSE)
  testthat::expect_gt(length(fast$operations), 0L)
  testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["queries"]), numeric(1L))), 0)
  testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["cache_hits"]), numeric(1L))), 0)
  for (nm in c("att", "inffunc", "V_analytical", "se"))
    testthat::expect_identical(fast$result[[nm]], slow$result[[nm]])
  for (nm in c("available", "target_meta", "support_status", "support_reason", "covariance_handling"))
    testthat::expect_identical(fast$result$chained[[nm]], slow$result$chained[[nm]])
  for (op in fast$operations) testthat::expect_error(op$pair(1L, list()), class = "cdid_cp04_integrity_error")
})


testthat::test_that("source-produced generic and imputation pair forms use cached native positions", {
  fx <- v49_fixture$cp04_fixture(); rows <- integer(23L); rows[] <- 1:23
  selection <- v49_ns("fit_selection_models")(reader = fx$reader, link_views = fx$views,
    control = v49_fixture$cp04_control, glm_control = v49_fixture$cp04_glm_control,
    selection_model = "pair_logit", selection_pooling = "within_cohort")
  # Both descriptors come from their production constructors, not an invented
  # class or the older standalone-context descriptor.
  generic <- selection$pair_cohorts[[1L]]$pair
  st <- fx$reader$state
  records <- lapply(fx$views, function(v) list(pair = v$pair,
    target = as.integer(match(v$cohorts$g, st$cohorts)),
    controls = as.integer(match(sort(v$cohorts$K), st$cohorts)), link_id = v$link_id,
    source_functional_key = v$functional_key, functional_key = v$functional_key))
  catalog <- v49_ns("imputation_build_catalog")(records, st$cohorts, "ipw", "separate")
  reconstructed <- catalog$cell_specs[[1L]]$pair
  testthat::expect_identical(names(generic), c("s", "t", "s_calendar_index", "t_calendar_index"))
  testthat::expect_identical(names(reconstructed), c(names(fx$views[[1L]]$pair), "L", "pair_label", "pair_id"))
  op <- v49_ns("cp04_pair_read_operation")(fx$reader, rows, 1e7)
  for (descriptor in list(generic, reconstructed)) {
    before <- op$counts()
    testthat::expect_identical(op$pair(rows, descriptor, TRUE), v49_ns("cp04_pair_rows")(fx$reader, rows, descriptor, TRUE))
    after <- op$counts()
    testthat::expect_equal(after[["queries"]], before[["queries"]] + 1)
    testthat::expect_identical(after[["fallbacks"]], before[["fallbacks"]])
  }
  # Index-free records retain the old timestamp-match path. Missing or invalid
  # calendar positions likewise fall back before a new native query or error.
  plain <- fx$views[[1L]]$pair[c("s", "t")]
  testthat::expect_identical(op$pair(rows, plain, TRUE), v49_ns("cp04_pair_rows")(fx$reader, rows, plain, TRUE))
  for (bad in list(NA_real_, 0, length(fx$panel$all_times) + 1)) {
    descriptor <- generic; descriptor$s_calendar_index <- bad
    original <- v49_error(function() v49_ns("cp04_pair_rows")(fx$reader, rows, descriptor, TRUE))
    actual <- v49_error(function() op$pair(rows, descriptor, TRUE))
    if (inherits(original, "condition")) {
      testthat::expect_identical(class(actual), class(original))
      testthat::expect_identical(conditionMessage(actual), conditionMessage(original))
    } else testthat::expect_identical(actual, original)
  }
  op$finish(); op$abort()
})


testthat::test_that("warmed payload lookup keeps live native verification without a second resolver call", {
  x <- list(rows = c(3L, 1L), values = c(1.25, -2), label = "\u03bb\u6f22")
  expected <- v49_ns("cp04_payload")(x)
  original <- v49_ns("gmm_require_native_dll"); calls <- 0L
  testthat::local_mocked_bindings(gmm_require_native_dll = function() {
    calls <<- calls + 1L; original()
  }, .package = "cdid")
  testthat::expect_identical(v49_ns("cp04_payload")(x), expected)
  testthat::expect_identical(calls, 1L)
})
