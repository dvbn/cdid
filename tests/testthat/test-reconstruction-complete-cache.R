v40_ns <- function(name) getFromNamespace(name, "cdid")
v40_old <- local({
  e <- new.env(parent = environment())
  wanted <- c("v39_ns", "v39_helpers", "v39_pure", "v39_aug", "v39_env")
  for (x in parse(file = testthat::test_path("test-reconstruction-active-scalar.R"))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  e
})

testthat::test_that("complete cache exposes actual seen NULLs and final finite centered state only", {
  cells <- list(list(group_id = 1L, mean_d = 0), list(group_id = 1L, mean_d = 0),
    list(group_id = 2L, mean_d = 0))
  groups <- list(list(p = 1L, cell_ids = 1:2), list(p = 1L, cell_ids = 3L))
  state <- list(core = list(cells = cells, groups = groups), cells = cells, groups = groups,
    intercept = TRUE, cids = c(3L, 1L, 2L), N = 4L, B = 4L)
  meta <- list(cells = cells, groups = lapply(groups, function(g) c(g, list(kept_columns = "x"))),
    control = list(max_chunk_rows = 4L), shared_or = NULL)
  X <- matrix(1, 4L, 1L, dimnames = list(NULL, "x"))
  plan <- v40_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, X, 100L)
  calls <- integer(3L)
  raw <- function(cid) {
    calls[cid] <<- calls[cid] + 1L
    if (cid == 2L) NULL else list(rr = 1:4, d_rr = X, s_rr = rep(cid, 4L))
  }
  reader <- v40_ns("reconstruction_cached_cell_reader")(raw, plan, 3L)
  reader(1L)
  testthat::expect_null(v40_ns("reconstruction_complete_cache_witness")(
    v40_ns("reconstruction_score_cache_view")(reader), state))
  reader(2L); reader(3L)
  view <- v40_ns("reconstruction_score_cache_view")(reader)
  testthat::expect_true(all(view$seen)); testthat::expect_null(view$cache[[view$position[2L]]])
  snapshot <- v40_ns("reconstruction_complete_cache_witness")(view, state)
  testthat::expect_type(snapshot, "list")
  testthat::expect_identical(calls, rep(1L, 3L))
  # One loading cell touches every member of its coupled group, in needed order.
  testthat::expect_identical(v40_ns("reconstruction_sparse_adjoint_ids")(snapshot, 1L), c(1L, 2L))
  testthat::expect_identical(v40_ns("reconstruction_sparse_adjoint_ids")(snapshot, c(3L, 1L)), state$cids)
  testthat::expect_null(v40_ns("reconstruction_complete_cache_get")(view, 2L))
  testthat::expect_identical(v40_ns("reconstruction_complete_cache_get")(view, 1L), raw(1L))
  for (field in c("intercept", "cids", "N")) {
    changed <- state
    changed[[field]] <- switch(field, intercept = FALSE, cids = rev(state$cids), N = 5L)
    testthat::expect_false(identical(changed, snapshot$state))
  }
  state_bad <- state; state_bad$cells[[1L]]$mean_d <- Inf
  testthat::expect_null(v40_ns("reconstruction_complete_cache_witness")(view, state_bad))
  state_bad$cells[[1L]]$mean_d <- .Machine$double.xmax
  bad_view <- view; bad_view$cache[[1L]]$d_rr[] <- -.Machine$double.xmax
  testthat::expect_null(v40_ns("reconstruction_complete_cache_witness")(bad_view, state_bad))
  bad_view <- view; class(bad_view$cache[[1L]]$d_rr) <- c("not_plain", "matrix", "array")
  testthat::expect_null(v40_ns("reconstruction_complete_cache_witness")(bad_view, state))
  fake <- function(cid) stop("must not call raw reader")
  attr(fake, "reconstruction_cache_view") <- function() stop("must not call fake getter")
  testthat::expect_null(v40_ns("reconstruction_score_cache_view")(fake))
  # A later actual first-pass read can replace an earlier cached record. The
  # final eligibility pass checks the final cache rather than stale witnesses.
  changing <- NULL
  raw_change <- function(cid) {
    value <- raw(cid)
    if (cid == 3L) {
      f <- environment(changing); cache <- get("cache", f, inherits = FALSE)
      cache[[1L]]$s_rr[] <- Inf; assign("cache", cache, envir = f)
    }
    value
  }
  changing <- v40_ns("reconstruction_cached_cell_reader")(raw_change, plan, 3L)
  changing(1L); changing(2L); changing(3L)
  testthat::expect_null(v40_ns("reconstruction_complete_cache_witness")(
    v40_ns("reconstruction_score_cache_view")(changing), state))
})

v40_run <- function(fit, kind, optimized, cache = "full", mutation = NULL, events = NULL) {
  e <- v40_old$v39_env(FALSE)
  factory <- v40_ns("reconstruction_cached_cell_reader")
  witness <- v40_ns("reconstruction_complete_cache_witness")
  accepted <- logical(); visits <- integer(); scalar_visits <- 0L; raw_calls <- 0L
  first_done <- FALSE; between_done <- FALSE; core_frame <- NULL
  e$reconstruction_complete_cache_witness <- function(view, state) {
    value <- if (optimized) witness(view, state) else NULL
    accepted <<- c(accepted, !is.null(value)); first_done <<- TRUE; value
  }
  e$reconstruction_complete_cache_get <- function(view, cid) {
    visits <<- c(visits, cid)
    if (!is.null(events) && between_done) events$gets <- events$gets + 1L
    v40_ns("reconstruction_complete_cache_get")(view, cid)
  }
  e$reconstruction_cached_cell_reader <- function(read, plan, n_cells) {
    force(read); core_frame <<- environment(read)
    # These plans also match the pre passed below, unlike the corrected V39
    # test's earlier inconsistent no-cache setup.
    counted <- function(cid) {
      raw_calls <<- raw_calls + 1L
      value <- read(cid)
      if (identical(mutation, "earlier_mean") && cid == max(pre$cids)) {
        f <- environment(read); x <- get("core_fit", f, inherits = FALSE)
        x$cells[[fit$links[[1L]]$target_cell_id]]$mean_d[] <- Inf
        assign("core_fit", x, envir = f)
      }
      value
    }
    value <- factory(counted, plan, n_cells)
    if (optimized) value else function(cid) { scalar_visits <<- scalar_visits + 1L; value(cid) }
  }
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  rows <- 1:11; links <- seq_along(fit$links)
  pre <- v40_ns(if (kind == "augmented") "augmented_if_budget" else "imputation_if_budget")(
    fit, length(rows), links, FALSE)
  if (cache == "partial") pre$score_cache$cids <- head(pre$score_cache$cids, 1L)
  if (cache == "none") pre$score_cache <- list(enabled = FALSE, cids = integer(0L), entries = 0)
  pair_reader <- function(r, pair, read_outcomes = TRUE) {
    if (kind == "augmented" && first_done && !between_done && !is.null(mutation) &&
        mutation %in% c("between_mean", "between_center", "between_needed", "between_N")) {
      if (mutation == "between_mean") {
        x <- get("cells", core_frame, inherits = FALSE)
        x[[fit$links[[1L]]$target_cell_id]]$mean_d[] <- Inf
        assign("cells", x, envir = core_frame)
      }
      if (mutation == "between_center") {
        x <- get("recipe", core_frame, inherits = FALSE); x$intercept <- !x$intercept
        assign("recipe", x, envir = core_frame)
      }
      if (mutation == "between_needed") {
        x <- get("sets", core_frame, inherits = FALSE); x$cids <- rev(x$cids)
        assign("sets", x, envir = core_frame)
      }
      if (mutation == "between_N") assign("N", 2 * get("N", core_frame, inherits = FALSE), envir = core_frame)
      between_done <<- TRUE
      if (!is.null(events)) events$mutated <- TRUE
    }
    v40_ns("cp04_pair_rows")(reader, r, pair, read_outcomes)
  }
  out <- if (kind == "augmented") e$augmented_evaluate_if_chunk_core(st = reader$state,
    X_F = v40_ns("augmented_backend_X_F")(fit), pair_reader = pair_reader, backend = fit,
    rows = rows, links = links, scale = "frame", details = FALSE, sets = pre) else
    e$imputation_evaluate_if_chunk_core(st = reader$state, X_F = v40_ns("imputation_backend_X_F")(fit, reader),
      pair_reader = pair_reader, core_fit = fit, recipe = fit$outcome_recipe, or_frame_index = fit$or_frame_index,
      trt_recipe = fit$treatment_recipe, est_method = fit$est_method, rows = rows, links = links,
      scale = "frame", details = FALSE, sets = pre)
  list(out = out, accepted = accepted, visits = visits, scalar_visits = scalar_visits,
    raw_calls = raw_calls, needed = pre$cids)
}

testthat::test_that("full-cache ordinary replay actually visits only ordered active closures", {
  for (kind in c("imputation", "augmented")) {
    a <- if (kind == "augmented") v40_old$v39_aug$rap_case(route_kind = "augmented") else
      v40_old$v39_pure$rcp_case(route_kind = "imputation")
    fit <- if (kind == "augmented") v40_old$v39_aug$rap_fit(a, mode = "pooled") else
      v40_old$v39_pure$rcp_fit(a, method = "ipw", mode = "pooled")
    current <- v40_run(fit, kind, TRUE); old <- v40_run(fit, kind, FALSE)
    testthat::expect_true(any(current$accepted))
    testthat::expect_gt(length(current$visits), 0L)
    testthat::expect_equal(current$out$phi_frame, old$out$phi_frame, tolerance = 1e-10)
    testthat::expect_equal(crossprod(current$out$phi_frame), crossprod(old$out$phi_frame), tolerance = 1e-10)
    testthat::expect_true(all(current$visits %in% current$needed))
    if (kind == "imputation") {
      testthat::expect_equal(current$raw_calls, length(current$needed))
      testthat::expect_lt(current$raw_calls + length(current$visits), old$scalar_visits)
    }
    for (mode in c("partial", "none")) {
      partial <- v40_run(fit, kind, TRUE, mode); plain <- v40_run(fit, kind, FALSE, mode)
      testthat::expect_false(any(partial$accepted))
      testthat::expect_length(partial$visits, 0L)
      testthat::expect_equal(partial$out$phi_frame, plain$out$phi_frame, tolerance = 1e-10)
    }
  }
})

testthat::test_that("a later first-pass callback cannot certify a stale finite mean", {
  a <- v40_old$v39_pure$rcp_case(route_kind = "imputation")
  a$views <- a$views[1L]
  fit <- v40_old$v39_pure$rcp_fit(a, method = "ipw", mode = "pooled")
  current <- tryCatch(v40_run(fit, "imputation", TRUE, mutation = "earlier_mean"), error = identity)
  old <- tryCatch(v40_run(fit, "imputation", FALSE, mutation = "earlier_mean"), error = identity)
  testthat::expect_s3_class(current, "error")
  testthat::expect_identical(class(current), class(old))
  testthat::expect_identical(conditionMessage(current), conditionMessage(old))
})


testthat::test_that("actual augmented cell reads invalidate changed local arithmetic state before sparse arms", {
  a <- v40_old$v39_aug$rap_case(route_kind = "augmented")
  fit <- v40_old$v39_aug$rap_fit(a, mode = "pooled")
  for (what in c("between_mean", "between_center", "between_needed", "between_N")) {
    seen <- new.env(parent = emptyenv()); seen$gets <- 0L; seen$mutated <- FALSE
    plain <- new.env(parent = emptyenv()); plain$gets <- 0L; plain$mutated <- FALSE
    current <- tryCatch(v40_run(fit, "augmented", TRUE, mutation = what, events = seen), error = identity)
    old <- tryCatch(v40_run(fit, "augmented", FALSE, mutation = what, events = plain), error = identity)
    testthat::expect_true(seen$mutated); testthat::expect_true(plain$mutated)
    testthat::expect_identical(seen$gets, 0L)
    if (inherits(old, "error")) {
      testthat::expect_identical(class(current), class(old))
      testthat::expect_identical(conditionMessage(current), conditionMessage(old))
    } else testthat::expect_equal(current$out$phi_frame, old$out$phi_frame, tolerance = 1e-10)
  }
})


testthat::test_that("source annotations do not weaken exact cached-reader executable identity", {
  raw <- function(cid) NULL
  reader <- v40_ns("reconstruction_cached_cell_reader")(raw,
    list(enabled = TRUE, cids = 1L, entries = 10000), 1L)
  norm <- v40_ns("reconstruction_cache_reader_body")
  original <- body(reader)
  testthat::expect_identical(norm(original), utils::removeSource(original))
  testthat::expect_type(v40_ns("reconstruction_score_cache_view")(reader), "list")
  # Restore canonical binding after each private semantic-negative witness so
  # every earlier factory identity guard passes; the body guard must reject.
  check_body <- function(b) {
    changed <- reader; body(changed) <- b
    f <- environment(reader)
    assign("reader", changed, envir = f)
    on.exit(assign("reader", reader, envir = f))
    v40_ns("reconstruction_score_cache_view")(changed)
  }
  literal <- original; literal[[3L]][[2L]][[3L]] <- 1L
  testthat::expect_null(check_body(literal))
  flow <- original; flow[[4L]][[2L]] <- TRUE
  testthat::expect_null(check_body(flow))
  extra <- original; attr(extra, "execution_marker") <- 1L
  testthat::expect_null(check_body(extra))
  extra_inner <- original; attr(extra_inner[[4L]][[3L]], "execution_marker") <- 1L
  testthat::expect_null(check_body(extra_inner))
  malformed <- original; malformed[[4L]][[3L]] <- 1L
  testthat::expect_null(check_body(malformed))
  testthat::expect_null(check_body(quote({ slot <- position[cid]; cache[[slot]] })))
})
