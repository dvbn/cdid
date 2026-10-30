# SPDX-License-Identifier: GPL-3.0-only
# Installed cdid only; one accepted finite fixture, no simulated data.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE))
  stop("Actual installed cdid required; no source fallback.", call. = FALSE)
gram_api <- asNamespace("cdid")

# Same local planner clone as test-cp04-stream-plan-full-block.R.
gram_namespace <- function(legacy = FALSE, perturb = FALSE) {
  e <- new.env(parent = gram_api)
  for (nm in c("cp04_stream_covariance", "cp04_stream_plan")) {
    f <- get(nm, gram_api); environment(f) <- e; assign(nm, f, e)
  }
  if (legacy) e$cp04_stream_plan_full_block <- function(...) NULL
  if (perturb) e$crossprod <- function(x, y = NULL) {
    if (!is.null(y) && !identical(x, y)) return(base::crossprod(x, y))
    z <- if (is.null(y)) base::crossprod(x) else base::crossprod(x, y)
    # Move only the lower entry one representable double away from zero.
    v <- z[1L, 2L]
    stopifnot(is.finite(v))
    bits <- writeBin(as.double(v), raw(), size = 8L, endian = "little")
    for (k in seq_len(8L)) {
      b <- as.integer(bits[k]); bits[k] <- as.raw((b + 1L) %% 256L)
      if (b < 255L) break
    }
    z[2L, 1L] <- readBin(bits, "double", n = 1L, size = 8L, endian = "little")
    stopifnot(is.finite(z[2L, 1L]), z[2L, 1L] != v)
    z
  }
  e
}

gram_expect <- function(out, expected, means) {
  V <- out$covariance
  expect_true(is.matrix(V) && is.double(V) && all(is.finite(V)))
  expect_identical(V, t(V))
  expect_equal(unname(V), unname(expected), tolerance = 1e-12)
  expect_equal(unname(out$mean_phi), unname(means), tolerance = 1e-12)
}

test_that("finite CP04 Gram is exact symmetric in ordered full and legacy blocks", {
  dir <- testthat::test_path("fixtures", "augmented_finite_v1")
  fixtures <- new.env(parent = gram_api)
  sys.source(file.path(dir, "load_fixture.R"), envir = fixtures)
  fixtures$aug_fixture_verify_hashes(dir)
  fx <- fixtures$load_augmented_finite_reference(dir, "pair_sep_jsep",
    max_chunk_rows = 2L, ns = gram_api)
  # Explicit installed API from test-cp04-installed.R; no guessed fit adapters.
  fit <- gram_api$fit_cp04_corollary1(fx$reader, fx$views, fx$control, fx$glm_control)
  links <- c(4L, 2L, 1L, 3L)
  one_links <- links[1:2]
  one <- gram_namespace(); legacy <- gram_namespace(legacy = TRUE)
  plan <- gram_api$cp04_stream_plan(fit, one_links)
  old_plan <- legacy$cp04_stream_plan(fit, links)
  expect_identical(length(plan$blocks), 1L)
  expect_identical(unname(old_plan$blocks), list(1:2, 3:4))
  expect_identical(plan$chunk_rows, old_plan$chunk_rows)

  # Independent entrywise Gram of complete raw core IF rows, not of V itself.
  # This tests the covariance constructor, not the IF estimator's oracle.
  dom <- gram_api$cp04_backend_domain(fit)
  H <- matrix(NA_real_, dom$N, length(links))
  gram_api$cp04_replay_links(fit, links, function(b) H[b$rows, ] <<- b$phi,
    scale = "core", chunk_rows = plan$chunk_rows)
  expect_true(all(is.finite(H)))
  means <- colMeans(H); X <- sweep(H, 2L, means, "-")
  expected <- vapply(seq_along(links), function(j)
    vapply(seq_along(links), function(i) sum(X[, i] * X[, j]) / dom$n_A,
      numeric(1)), numeric(length(links)))
  full <- one$cp04_stream_covariance(fit, one_links)
  multi <- legacy$cp04_stream_covariance(fit, links)
  gram_expect(full, expected[1:2, 1:2], means[1:2])
  gram_expect(multi, expected, means)
  expect_equal(full$covariance, multi$covariance[1:2, 1:2], tolerance = 1e-12)
  for (force_legacy in c(FALSE, TRUE)) {
    mutated <- gram_namespace(legacy = force_legacy, perturb = TRUE)
    ix <- if (force_legacy) 1:4 else 1:2
    gram_expect(mutated$cp04_stream_covariance(fit, links[ix]), expected[ix, ix], means[ix])
  }

  # Verified reserve/bind/release API from test-cp08-retained-moment-if.R.
  retained <- function() {
    store <- gram_api$cp04_retained_if_reserve(fit, one_links, plan$chunk_rows, plan,
      retained = 0, gmm_retained = 0, live = 0, work = 0, control = fit$control)
    stopifnot(is.environment(store))
    on.exit(gram_api$cp04_retained_if_release(store), add = TRUE)
    out <- gram_api$cp04_stream_covariance(gram_api$cp04_retained_if_bind(fit, store), one_links)
    expect_identical(store$state, "filled")
    expect_identical(store$producer_passes, 1L)
    expect_identical(store$slice_passes, 1L)
    out
  }
  gram_expect(retained(), expected[1:2, 1:2], means[1:2])
})
