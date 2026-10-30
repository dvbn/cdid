# Standalone adapter checks: immutable source and stats only, no cdid oracle.
reference_start_inspect <- function(code, edit = identity) {
  calls <- list()
  walk <- function(node) {
    if (!is.call(node) && !is.expression(node)) return(node)
    if (is.call(node) && (identical(node[[1L]], quote(glm.fit)) ||
                          identical(node[[1L]], quote(stats::glm.fit)))) {
      calls[[length(calls) + 1L]] <<- node
      return(edit(node))
    }
    for (i in seq_along(node)) {
      if (!identical(node[[i]], quote(expr = ))) node[i] <- list(walk(node[[i]]))
    }
    node
  }
  result <- walk(code)
  list(code = result, calls = calls)
}

test_that("both immutable sources admit exactly ONE reversible matched-start adaptation", {
  hashes <- c(
    chained_ipw_link_original_6550dc.R = "6550dc592451cb7e9b7f45573d083fa14cacc1f3a09421a3330d31997183b0b5",
    weighted_pair_ipw_derivative_reference.R = "c96beff841db195f9958b9e60a8bd32ee8356c4f06b90e4ce0e5fc66bb2dcb45"
  )
  for (file in names(hashes)) {
    path <- testthat::test_path("reference", file)
    expect_identical(digest::digest(file = path, algo = "sha256"), hashes[[file]])
    original <- parse(path, keep.source = FALSE)
    matched <- ipw_reference_match_start(original, file)
    before <- reference_start_inspect(original)$calls
    after <- reference_start_inspect(matched)$calls
    expect_length(before, 1L)
    expect_length(after, 1L)
    expect_false("start" %in% names(before[[1L]]))
    design <- if (file == names(hashes)[1L]) quote(Z_use) else quote(Z[use,,drop=FALSE])
    expect_identical(after[[1L]]$start, substitute(rep(0, ncol(x)), list(x = design)))
    expect_identical(length(after[[1L]]), length(before[[1L]]) + 1L)
    reversed <- reference_start_inspect(matched, function(call) { call$start <- NULL; call })$code
    expect_identical(reversed, original)

    wrong <- reference_start_inspect(original, function(call) { call$weights <- quote(NULL); call })$code
    missing <- reference_start_inspect(original, function(call) quote(NULL))$code
    multiple <- as.expression(c(as.list(original), as.list(original)))
    expect_error(ipw_reference_match_start(wrong, file), "wrong glm.fit call")
    expect_error(ipw_reference_match_start(missing, file), "exactly one glm.fit call")
    expect_error(ipw_reference_match_start(multiple, file), "exactly one glm.fit call")
    expect_error(ipw_reference_match_start(matched, file), "wrong glm.fit call")
  }
  expect_error(ipw_reference_match_start(expression(NULL), "other.R"), "unsupported")
  expect_error(ipw_reference_source("other.R", new.env()), "unsupported")
  expect_error(ipw_reference_source(testthat::test_path("missing", names(hashes)[1L]), new.env()), "missing packaged reference")
})

test_that("the loaded original kernel differs only by the matched-start argument", {
  file <- "chained_ipw_link_original_6550dc.R"
  path <- testthat::test_path("reference", file)
  original <- new.env(parent = baseenv())
  matched <- new.env(parent = baseenv())
  sys.source(path, envir = original, keep.source = FALSE)
  ipw_reference_source(path, envir = matched)
  expect_identical(ls(matched), ls(original))
  for (name in ls(original)) {
    expect_identical(formals(matched[[name]]), formals(original[[name]]))
    reversed <- reference_start_inspect(body(matched[[name]]), function(call) { call$start <- NULL; call })$code
    expect_identical(reversed, body(original[[name]]))
  }
  expect_length(reference_start_inspect(body(matched$fit_ipw_link))$calls, 1L)
  expect_identical(matched$estimate_ipw_link, matched$fit_ipw_link)
})

test_that("the audit loader adapts the original kernel and leaves independent Newton untouched", {
  kernel <- audit_load_reference("chained_ipw_link_original_6550dc.R", new.env(parent = baseenv()))
  calls <- reference_start_inspect(body(kernel$fit_ipw_link))$calls
  expect_length(calls, 1L)
  expect_identical(calls[[1L]]$start, quote(rep(0, ncol(Z_use))))
  ref <- audit_newton_reference()
  original <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "cp02_shared_ipw_audit_newton_reference.R"), envir = original)
  expect_identical(body(ref$independent_fit), body(original$independent_fit))
  expect_identical(formals(ref$independent_fit), formals(original$independent_fit))
})

test_that("both adapted glm calls retain the independent ordinary weighted 2x2 MLE", {
  counts <- c(7L, 3L, 2L, 8L)
  x <- rep(c(0, 0, 1, 1), counts)
  y <- rep(c(0, 1, 0, 1), counts)
  for (cell_weights in list(rep(1, 4L), c(1, 2, 3, 4))) {
    w <- rep(cell_weights, counts)
    cells <- matrix(counts * cell_weights, 2L, byrow = TRUE)
    eta <- stats::qlogis(cells[, 2L] / rowSums(cells))
    beta <- c(eta[1L], eta[2L] - eta[1L])
    for (file in names(ipw_reference_sha)) {
      code <- parse(testthat::test_path("reference", file), keep.source = FALSE)
      call <- reference_start_inspect(ipw_reference_match_start(code, file))$calls[[1L]]
      env <- new.env(parent = asNamespace("stats"))
      env$Z_use <- env$Z <- cbind(1, x)
      env$D_use <- env$D <- y
      env$fit_weights <- env$case <- w
      env$ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)
      env$use <- seq_along(y); env$rr <- rep(1, length(y))
      fit <- eval(call, envir = env)
      expect_true(fit$converged)
      expect_false(fit$boundary)
      expect_lt(max(abs(unname(fit$coefficients) - beta)), 1e-12)
      p <- stats::plogis(drop(env$Z %*% fit$coefficients))
      expect_lt(max(abs(drop(crossprod(env$Z, w * (y - p))))), 1e-12)
    }
  }
})

test_that("the nested root loader changes initialization while retaining original derivative algebra", {
  file <- "weighted_pair_ipw_derivative_reference.R"
  path <- testthat::test_path("reference", file)
  code <- parse(path, keep.source = FALSE)
  # Stop before the existing finite-refit execution to obtain the default
  # reference function/data; the matched loader below runs its unchanged checks.
  last_definition <- which(vapply(as.list(code), function(node)
    is.call(node) && identical(node[[1L]], quote(`<-`)) &&
      identical(node[[2L]], quote(fit_link)), logical(1L)))
  expect_length(last_definition, 1L)
  original <- new.env(parent = asNamespace("stats"))
  eval(code[seq_len(last_definition)], envir = original)
  ref <- new.env(parent = globalenv())
  ref$sys.source <- ipw_reference_source
  ref$root <- new.env(parent = asNamespace("stats"))
  # The exact nested load expression in the immutable numerical reference.
  independent <- parse(testthat::test_path("reference", "cp02_weighted_pair_if_independent_reference.R"), keep.source = FALSE)
  nested <- Filter(function(node) is.call(node) && identical(node[[1L]], quote(`<-`)) &&
    identical(node[[2L]], quote(root_output)), as.list(independent))
  expect_length(nested, 1L)
  ref$root_path <- path
  eval(nested[[1L]], envir = ref)
  root <- ref$root
  expect_identical(root$Z, original$Z)
  expect_identical(root$Y, original$Y)
  expect_identical(root$G, original$G)
  expect_identical(root$observed, original$observed)
  reversed <- reference_start_inspect(body(root$fit_link), function(call) { call$start <- NULL; call })$code
  expect_identical(reversed, body(original$fit_link))
  old_fits <- lapply(root$pairs, original$fit_link, with_if = TRUE)
  # Initialization alone changes rounding on this ordinary fixture. The root
  # equations and strict 1e-14/1e-13 production comparisons are left untouched.
  expect_false(identical(lapply(old_fits, `[[`, "value"), lapply(root$fits, `[[`, "value")))
  expect_false(identical(do.call(cbind, lapply(old_fits, `[[`, "phi")), root$phi))
})
