v56_ns <- function(name) getFromNamespace(name, "cdid")
v56_fixture <- function(selection_design = c("plain", "public", "factor")) {
  selection_design <- match.arg(selection_design)
  # Reuse exactly the accepted cp04 deterministic design prefix, stopping
  # before qsolve/finite-refit execution. Its canonical cohorts are 3, 5, Inf.
  design <- readLines(testthat::test_path("reference", "cp04_finite_reference.R"), warn = FALSE)
  first_solver <- grep("^qsolve <-", design)
  stopifnot(length(first_solver) == 1L, first_solver > 1L)
  ref <- new.env(parent = asNamespace("stats"))
  invisible(capture.output(eval(parse(text = design[seq_len(first_solver - 1L)]), ref)))
  g <- c(3, 5, Inf)[ref$dat$k]; Y <- ref$Y; ids <- ref$ids
  long <- data.frame(id = rep(ids, each = 4L), time = rep(1:4, ref$N),
    y = as.vector(t(Y)), g = rep(g, each = 4L))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  long <- long[!long$id %in% ids[zero[-1L]], , drop = FALSE]
  ctrl <- v56_ns("ipw_owner_control")(max_chunk_rows = 23L,
    max_retained_entries = 1e7, max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7)
  panel <- v56_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  Z <- ref$Z
  frame <- v56_ns("construct_ipw_frame_owner")(panel, ids, Z, ref$w, ctrl,
    frame_g = g, allow_frame_only = TRUE)
  SZ <- if (selection_design == "plain") Z else if (selection_design == "public")
    stats::model.matrix(~x, data = ref$dat) else stats::model.matrix(~factor(x), data = ref$dat)
  if (selection_design != "plain") attr(SZ, "intercept") <- TRUE
  reader <- v56_ns("construct_cp04_frame_reader")(frame, ids, SZ, g, ref$w, ctrl)
  make <- function(s, t, K, id) v56_ns("construct_ipw_link_view")(frame, s, t, g = 3,
    K = K, control = ctrl, link_id = id, functional_key = paste0("cmpcat:v56:", id))
  views <- list(make(1, 2, c(5, Inf), "NYT12"), make(2, 4, c(5, Inf), "NYT24"),
    make(1, 2, Inf, "NT12"))
  backend <- v56_ns("fit_cp04_corollary1")(reader, views, ctrl,
    list(epsilon = 1e-13, maxit = 60L, trace = FALSE))
  list(backend = backend, reader = reader, frame = frame, N = length(ids), zero = zero,
    single = which(rowSums(!is.na(Y)) == 1L))
}
v56_oracle <- function(env = new.env(parent = asNamespace("cdid"))) {
  sys.source(testthat::test_path("reference", "legacy_cp04_r_oracle.R"), env)
  env$cp04_v56_oracle
}
v56_run <- function(fx, rows, links = 1:3, native = TRUE) {
  capture <- new.env(parent = emptyenv()); capture$events <- list(); capture$complete <- list(); capture$pair_declines <- list()
  event <- function(x) capture$events[[length(capture$events) + 1L]] <- x
  helpers <- lapply(c("cp04_predict_selection", "cp04_solve", "cp04_finite",
    "ipw_logit_values", "ipw_clip_floor"), v56_ns)
  names(helpers) <- c("cp04_predict_selection", "cp04_solve", "cp04_finite", "ipw_logit_values", "ipw_clip_floor")
  predict <- helpers$cp04_predict_selection; solve <- helpers$cp04_solve; finite <- helpers$cp04_finite
  helpers$cp04_predict_selection <- function(reader, rows, spec, selection) {
    event(list("predict", rows, spec$cohort, spec$pair)); predict(reader, rows, spec, selection)
  }
  helpers$cp04_solve <- function(R, b) { event(list("solve", dim(b))); solve(R, b) }
  helpers$cp04_finite <- function(x, label, positive = FALSE) {
    if (identical(label, "complete frame IF")) capture$complete[[length(capture$complete) + 1L]] <- unname(x)
    finite(x, label, positive)
  }
  issue <- v56_ns("cp04_pair_read_operation")
  logged_issue <- function(reader, rows, funding_entries) {
    op <- issue(reader, rows, funding_entries)
    if (is.null(op)) {
      panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
      B <- length(rows); T <- length(panel$postings$date_start)
      info <- list(rows = rows, links = links, B = B, funding = funding_entries,
        pair_reservation = 4096 + 64 * B + 8 * T + 3 * B * T)
      capture$pair_declines[[length(capture$pair_declines) + 1L]] <- info
      message("V56 actual optional pair decline: rows=", paste(rows, collapse = ","),
        "; links=", paste(links, collapse = ","), "; B=", B, "; funding=", funding_entries,
        "; reservation=", info$pair_reservation)
      return(NULL)
    }
    pair <- op$pair
    op$pair <- function(rows, spec, read_outcomes = FALSE) {
      event(list("pair", rows, spec, read_outcomes)); pair(rows, spec, read_outcomes)
    }
    op
  }
  if (!native) {
    env <- new.env(parent = asNamespace("cdid"))
    for (name in names(helpers)) assign(name, helpers[[name]], env)
    env$cp04_pair_read_operation <- logged_issue
    value <- v56_oracle(env)(fx$backend, rows, links)
    return(list(value = value, complete = capture$complete, events = capture$events,
      pair_declines = capture$pair_declines, admitted = !length(capture$pair_declines)))
  }
  pre <- v56_ns("cp04_if_budget")(fx$backend, length(rows), links)
  st <- fx$reader$state; B <- length(rows); d <- ncol(st$Z); p <- fx$frame$n_covariates
  widths <- length(pre$sids) * (d + 2) + sum(vapply(fx$backend$treatment_specs[pre$tids],
    function(s) 2 * p + 10 + length(s$cohorts) * (d + 1), numeric(1))) +
    3 * length(st$cohorts) + length(links) + 16
  funding <- 8 * B * (widths + p + d + 4)
  op <- logged_issue(fx$reader, rows, funding)
  if (is.null(op)) return(list(value = NULL, complete = capture$complete, events = capture$events,
    pair_declines = capture$pair_declines, admitted = FALSE))
  on.exit(op$abort(), add = TRUE)
  si <- base::getNativeSymbolInfo("C_gmm_legacy_ipw_chunk", PACKAGE = base::getLoadedDLLs()[["cdid"]])
  models <- .Call(si$address, fx$reader, fx$frame, fx$backend, rows, pre$tids, pre$sids,
    op$pair, helpers, pre$scratch - funding)
  testthat::expect_false(is.null(models))
  phi <- matrix(0, B, length(links), dimnames = list(NULL,
    vapply(fx$backend$links[links], `[[`, character(1), "link_id")))
  for (j in seq_along(links)) phi[, j] <- models[[fx$backend$links[[links[j]]]$model]]$phi
  op$finish()
  list(value = phi, complete = capture$complete, events = capture$events,
    pair_declines = capture$pair_declines, admitted = TRUE)
}

testthat::test_that("the legacy worker preserves full-F scores, all coordinates, zeros and query/solve order", {
  fx <- v56_fixture()
  testthat::expect_identical(class(fx$backend), "cdid_cp04_corollary1")
  testthat::expect_identical(class(fx$reader), "cdid_cp04_frame_reader")
  testthat::expect_identical(class(fx$frame), c("cdid_ipw_frame_owner", "cdid_ipw_context"))
  testthat::expect_identical(as.integer(table(fx$reader$state$g)), c(41L, 51L, 59L))
  testthat::expect_identical(as.integer(table(fx$reader$state$g[fx$reader$state$a_index > 0L])), c(28L, 29L, 28L))
  testthat::expect_gt(length(fx$zero), 0L); testthat::expect_gt(length(fx$single), 0L)
  batches <- lapply(seq.int(1L, fx$N, by = 23L), function(a) seq.int(a, min(a + 22L, fx$N)))
  batches <- c(batches, list(c(23L, 1L, 19L, 4L, 13L, 2L), as.double(c(3L, 2L, 1L)),
    as.double(c(fx$zero[1L], fx$single[1L], 1L))))
  native_batches <- fallback_batches <- 0L
  for (rows in batches) for (links in list(1:3, c(3L, 1L))) {
    ref <- v56_run(fx, rows, links, FALSE); got <- v56_run(fx, rows, links, TRUE)
    if (got$admitted) {
      native_batches <- native_batches + 1L
      testthat::expect_equal(got$value, ref$value$phi, tolerance = 1e-12)
      testthat::expect_identical(got$value == 0, ref$value$phi == 0)
      testthat::expect_equal(got$complete, ref$complete, tolerance = 1e-12)
      testthat::expect_identical(lapply(got$complete, function(x) x == 0), lapply(ref$complete, function(x) x == 0))
      testthat::expect_identical(got$events, ref$events)
    } else {
      fallback_batches <- fallback_batches + 1L
      testthat::expect_false(ref$admitted)
      testthat::expect_identical(got$pair_declines, ref$pair_declines)
    }
    # Real public execution, including naturally unsupported optional cache
    # batches, must preserve all output fields and full-F zero/singleton rows.
    actual <- v56_ns("evaluate_cp04_if_chunk")(fx$backend, rows, links)
    testthat::expect_equal(actual, ref$value, tolerance = 1e-12)
    testthat::expect_identical(actual$phi == 0, ref$value$phi == 0)
    if (!got$admitted && any(rows %in% c(fx$zero, fx$single)))
      testthat::expect_true(any(actual$phi[rows %in% c(fx$zero, fx$single), , drop = FALSE] != 0))
  }
  testthat::expect_gt(native_batches, 0L)
  testthat::expect_gt(fallback_batches, 0L)
  message("V56 actual first-test counts: native=", native_batches, "; ordinary fallback=", fallback_batches)
})

testthat::test_that("the canonical frame wrapper follows the original bounded R row subset", {
  fx <- v56_fixture()
  # Validate the canonical matrix header without requiring an R-internal
  # representation; root's actual wrapper probe remains separate evidence.
  testthat::expect_identical(typeof(fx$frame$Z_a), "double")
  testthat::expect_identical(dim(fx$frame$Z_a),
    as.integer(c(fx$frame$n_a_units, fx$frame$n_covariates)))
  testthat::expect_false(is.object(fx$frame$Z_a))
  testthat::expect_true(all(names(attributes(fx$frame$Z_a)) %in% c("dim", "dimnames")))
  rows <- c(23L, 1L, 19L, 4L, 13L, 2L)
  ref <- v56_run(fx, rows, native = FALSE)
  got <- v56_run(fx, rows, native = TRUE)
  testthat::expect_equal(got$value, ref$value$phi, tolerance = 1e-12)
  testthat::expect_equal(got$complete, ref$complete, tolerance = 1e-12)
  testthat::expect_identical(got$events, ref$events)
  testthat::expect_identical(got$value == 0, ref$value$phi == 0)
})

testthat::test_that("ordinary public and factor model-matrix metadata preserve native full coordinates", {
  for (design in c("public", "factor")) {
    fx <- v56_fixture(design); rows <- seq_len(23L)
    Z <- fx$reader$state$Z
    testthat::expect_identical(typeof(Z), "double")
    testthat::expect_false(is.object(Z))
    testthat::expect_true(isTRUE(attr(Z, "intercept")))
    testthat::expect_true(is.integer(attr(Z, "assign")))
    if (design == "factor") testthat::expect_true(is.list(attr(Z, "contrasts")))
    ref <- v56_run(fx, rows, native = FALSE)
    got <- v56_run(fx, rows, native = TRUE)
    testthat::expect_true(got$admitted)
    testthat::expect_equal(got$value, ref$value$phi, tolerance = 1e-12)
    testthat::expect_equal(got$complete, ref$complete, tolerance = 1e-12)
    testthat::expect_identical(got$events, ref$events)
    testthat::expect_identical(got$value == 0, ref$value$phi == 0)
    actual <- v56_ns("evaluate_cp04_if_chunk")(fx$backend, rows, 1:3)
    testthat::expect_equal(actual, ref$value, tolerance = 1e-12)
    pre <- v56_ns("cp04_if_budget")(fx$backend, length(rows), 1:3)
    op <- v56_ns("cp04_pair_read_operation")(fx$reader, rows, 1e6)
    testthat::expect_false(is.null(op))
    run <- v56_ns("cp04_legacy_ipw_native")
    testthat::expect_false(is.null(run(fx$backend, fx$reader, fx$frame, rows, pre$tids, pre$sids, op, pre$scratch)))
    bad <- fx$reader; attr(bad$state$Z, "unrelated") <- TRUE
    testthat::expect_null(run(fx$backend, bad, fx$frame, rows, pre$tids, pre$sids, op, pre$scratch))
    op$abort()
  }
})

testthat::test_that("legacy predictor boundaries and hidden-coordinate overflow retain the first refusal", {
  fx <- v56_fixture(); rows <- seq_len(23L)
  check <- function(fx) {
    failure <- function(native) tryCatch(v56_run(fx, rows, native = native), error = identity)
    ref <- failure(FALSE); got <- failure(TRUE)
    testthat::expect_s3_class(ref, "cdid_cp04_numeric_error")
    testthat::expect_identical(class(got), class(ref))
    testthat::expect_identical(conditionMessage(got), conditionMessage(ref))
  }
  bad <- fx
  bad$backend$selections[[1L]]$coefficients[] <- 1000
  bad$backend <- v56_ns("ipw_owner_seal")(bad$backend, "cp04_backend")
  ref <- tryCatch(v56_run(bad, rows, native = FALSE), error = identity)
  got <- tryCatch(v56_run(bad, rows, native = TRUE), error = identity)
  testthat::expect_s3_class(ref, "cdid_cp04_boundary_error")
  testthat::expect_identical(class(got), class(ref)); testthat::expect_identical(conditionMessage(got), conditionMessage(ref))
  detailed <- v56_oracle()(fx$backend, rows, 1:3, details = TRUE)
  up <- detailed$details$models[[1L]]$upstream
  i <- which.max(rowSums(abs(up)))
  testthat::expect_gt(sum(abs(up[i, ])), 2)
  for (coordinate in c(3L, 4L)) {
    bad <- fx
    bad$backend$treatments[[1L]]$partials$upstream[coordinate, ] <- 1e308 * sign(up[i, ])
    bad$backend <- v56_ns("ipw_owner_seal")(bad$backend, "cp04_backend")
    check(bad)
  }
})

testthat::test_that("legacy helper mutations and unsupported operands decline before native numerics", {
  fx <- v56_fixture(); rows <- seq_len(23L); pre <- v56_ns("cp04_if_budget")(fx$backend, 23L, 1:3)
  issue <- v56_ns("cp04_pair_read_operation")
  op <- issue(fx$reader, rows, 1e6); testthat::expect_false(is.null(op)); on.exit(op$abort(), add = TRUE)
  run <- v56_ns("cp04_legacy_ipw_native")
  testthat::expect_false(is.null(run(fx$backend, fx$reader, fx$frame, rows, pre$tids, pre$sids, op, pre$scratch)))
  for (classes in list("cdid_ipw_frame_owner", c("cdid_ipw_frame_owner", "cdid_ipw_context", "extra"))) {
    bad_frame <- fx$frame; class(bad_frame) <- classes
    testthat::expect_null(run(fx$backend, fx$reader, bad_frame, rows, pre$tids, pre$sids, op, pre$scratch))
  }
  bad_reader <- fx$reader; class(bad_reader) <- c(class(bad_reader), "extra")
  testthat::expect_null(run(fx$backend, bad_reader, fx$frame, rows, pre$tids, pre$sids, op, pre$scratch))
  bad_backend <- fx$backend; class(bad_backend) <- c(class(bad_backend), "extra")
  testthat::expect_null(run(bad_backend, fx$reader, fx$frame, rows, pre$tids, pre$sids, op, pre$scratch))
  bad <- fx$backend; bad$selections[[1L]]$H_chol <- structure(bad$selections[[1L]]$H_chol, unsupported = TRUE)
  testthat::expect_null(run(bad, fx$reader, fx$frame, rows, pre$tids, pre$sids, op, pre$scratch))
  testthat::expect_null(run(fx$backend, fx$reader, fx$frame, rows, pre$tids, pre$sids, op, 0))
  testthat::local_mocked_bindings(cp04_solve = function(...) stop("mutated solve entered"), .package = "cdid")
  testthat::expect_null(run(fx$backend, fx$reader, fx$frame, rows, pre$tids, pre$sids, op, pre$scratch))
})


testthat::test_that("underfunded native scratch retains the unchanged ordinary R fallback", {
  fx <- v56_fixture(); rows <- c(23L, 1L, 19L, 4L, 13L, 2L)
  ref <- v56_oracle()(fx$backend, rows, 1:3)
  run <- v56_ns("cp04_legacy_ipw_native"); calls <- 0L
  testthat::local_mocked_bindings(cp04_legacy_ipw_native = function(backend, reader, frame,
      rows, tids, sids, pair_op, funding_entries) {
    calls <<- calls + 1L
    out <- run(backend, reader, frame, rows, tids, sids, pair_op, funding_entries = 0)
    testthat::expect_null(out)
    out
  }, .package = "cdid")
  got <- v56_ns("evaluate_cp04_if_chunk")(fx$backend, rows, 1:3)
  testthat::expect_identical(calls, 1L)
  testthat::expect_equal(got, ref, tolerance = 1e-12)
  testthat::expect_identical(got$phi == 0, ref$phi == 0)
})
