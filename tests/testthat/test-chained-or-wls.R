# Source-only verification supplies this environment; installed-package tests
# resolve the same private helpers from the candidate namespace.
if (!exists("or_wls_helpers", inherits = TRUE)) {
  or_wls_helpers <- asNamespace("cdid")
}

or_wls_fixture <- function() {
  # Each cohort alone has deficient rank: x is cohort-constant.
  # Pools A+B and B+C are full rank. The outcome has residual variation.
  k <- rep(c(-1, 0, 2), each = 6L)
  u <- rep(c(-1, 1, -1, 1, -1, 1), 3L)
  Z <- cbind(intercept = 1, x = k, u = u)
  delta <- 2 + 0.4 * k - 0.7 * u +
    rep(c(0.3, -0.2, 0.1, 0.5, -0.4, 0.2), 3L) * (1 + k / 5)
  w <- (1 + (seq_along(k) %% 7) / 4) * rep(c(0.8, 1.6, 1.1), each = 6L)
  list(Z = Z, delta = delta, w = w,
       cohorts = split(seq_along(k), rep(c("A", "B", "C"), each = 6L)),
       pools = list(AB = 1:12, BC = 7:18))
}

or_wls_make <- function(f, rows) {
  or_wls_helpers$or_wls_summary(f$Z[rows, , drop = FALSE], f$delta[rows], f$w[rows])
}

or_wls_objective <- function(summary, beta) {
  p <- length(beta)
  sum((summary$factor[, p + 1L] - summary$factor[, seq_len(p), drop = FALSE] %*% beta)^2)
}

# Fixed BEFORE verification: modest, well-scaled deterministic matrices and
# short QR chains. 1e-10 relative/absolute error allows ordinary double QR
# rounding; it is far below the deliberate 1% weight/column mutations.
or_wls_tol <- 1e-10
or_wls_close <- function(actual, expected) {
  expect_lt(max(c(0, abs(actual - expected))) /
              max(1, abs(expected)), or_wls_tol)
}

test_that("overlapping pools equal direct WLS and a fully interacted reference", {
  f <- or_wls_fixture()
  summaries <- lapply(f$cohorts, function(rows) or_wls_make(f, rows))
  pools <- list(AB = c("A", "B"), BC = c("B", "C"))
  fits <- lapply(pools, function(pool) {
    or_wls_helpers$or_wls_fit(or_wls_helpers$or_wls_merge(summaries[pool]),
                            rank_tol = 1e-12)
  })
  # Reference only: duplicated rows for the overlap, separate coefficient
  # blocks for each actual model. Production never constructs this matrix.
  X <- matrix(0, 24L, 6L)
  X[1:12, 1:3] <- f$Z[f$pools$AB, ]
  X[13:24, 4:6] <- f$Z[f$pools$BC, ]
  stacked_rows <- c(f$pools$AB, f$pools$BC)
  stacked <- stats::lm.wfit(X, f$delta[stacked_rows], f$w[stacked_rows], tol = 1e-12)
  expect_equal(stacked$rank, 6L)
  for (j in seq_along(fits)) {
    rows <- f$pools[[j]]
    direct <- stats::lm.wfit(f$Z[rows, ], f$delta[rows], f$w[rows], tol = 1e-12)
    or_wls_close(fits[[j]]$coefficients, direct$coefficients)
    or_wls_close(fits[[j]]$coefficients, stacked$coefficients[(3L * j - 2L):(3L * j)])
    or_wls_close(fits[[j]]$rss, sum(f$w[rows] * direct$residuals^2))
  }
  or_wls_close(sum(vapply(fits, `[[`, numeric(1), "rss")),
               sum(f$w[stacked_rows] * stacked$residuals^2))
  # Two distinct target loadings consume the SAME already fitted AB model.
  targets <- rbind(c(1, 0.25, -0.2), c(1, 1.5, 0.3))
  direct <- stats::lm.wfit(f$Z[1:12, ], f$delta[1:12], f$w[1:12], tol = 1e-12)
  or_wls_close(drop(targets %*% fits$AB$coefficients), drop(targets %*% direct$coefficients))
  expect_false(identical(fits$AB$coefficients, fits$BC$coefficients))
})

test_that("deficient cohorts, short chunks and empty chunks preserve the whole objective", {
  f <- or_wls_fixture()
  cohorts <- lapply(f$cohorts, function(rows) or_wls_make(f, rows))
  for (s in cohorts) {
    expect_equal(ncol(s$factor), 4L)
    expect_lte(nrow(s$factor), 4L)
    expect_error(or_wls_helpers$or_wls_fit(s, rank_tol = 1e-12),
                 class = "cdid_or_wls_rank_error")
  }
  empty <- or_wls_make(f, integer())
  expect_identical(dim(empty$factor), c(0L, 4L))
  expect_identical(empty$n, 0L)
  expect_identical(or_wls_helpers$or_wls_merge(list(empty, empty)), empty)
  expect_error(or_wls_helpers$or_wls_fit(empty, rank_tol = 1e-12),
               class = "cdid_or_wls_rank_error")
  chunks <- list(integer(), 1L, 2:3, 4:6, 7L, 8:9, 10:12, integer())
  parts <- lapply(chunks, function(rows) or_wls_make(f, rows))
  expect_identical(dim(parts[[2L]]$factor), c(1L, 4L))
  expect_error(or_wls_helpers$or_wls_fit(parts[[2L]], rank_tol = 1e-12),
               class = "cdid_or_wls_rank_error")
  merged <- or_wls_helpers$or_wls_merge(parts)
  nested <- or_wls_helpers$or_wls_merge(list(
    or_wls_helpers$or_wls_merge(parts[1:4]), or_wls_helpers$or_wls_merge(parts[5:8])))
  whole <- or_wls_make(f, 1:12)
  expect_equal(merged$n, 12L)
  expect_lte(nrow(merged$factor), 4L)
  # This coefficient-independent identity verifies every quadratic/cross term
  # and the response constant of the augmented objective.
  A <- cbind(f$Z[1:12, ], f$delta[1:12]) * sqrt(f$w[1:12])
  for (s in list(whole, merged, nested)) {
    or_wls_close(crossprod(s$factor), crossprod(A))
  }
  # Arbitrary coefficients, including beta=0, also test the objective directly.
  probes <- rbind(c(0, 0, 0), c(1, -2, 0.5), c(-0.3, 0.7, 1.9))
  for (beta in split(probes, row(probes))) {
    expected <- sum(f$w[1:12] * (f$delta[1:12] - f$Z[1:12, ] %*% beta)^2)
    for (s in list(whole, merged, nested,
                  or_wls_helpers$or_wls_merge(list(empty, whole, empty)))) {
      or_wls_close(or_wls_objective(s, beta), expected)
    }
  }
  direct <- stats::lm.wfit(f$Z[1:12, ], f$delta[1:12], f$w[1:12], tol = 1e-12)
  or_wls_close(or_wls_helpers$or_wls_fit(merged, rank_tol = 1e-12)$coefficients,
               direct$coefficients)
  # Isolated short/wide chunk objective must also survive its pivot.
  for (beta in split(probes, row(probes))) {
    for (rows in list(1L, 2:3)) {
      s <- or_wls_make(f, rows)
      or_wls_close(or_wls_objective(s, beta),
                   sum(f$w[rows] * (f$delta[rows] - f$Z[rows, , drop = FALSE] %*% beta)^2))
    }
  }
})

test_that("row/cohort/column permutations and nonidentity pivots solve Gram loadings", {
  f <- or_wls_fixture()
  # Move the largest predictor into the middle; both augmented and final
  # LAPACK QR must pivot on this fixture.
  f$Z[, 2L] <- 5 * f$Z[, 2L]
  A <- cbind(f$Z, f$delta) * sqrt(f$w)
  expect_false(identical(qr(A, LAPACK = TRUE)$pivot, 1:4))
  whole <- or_wls_make(f, 1:18)
  direct <- stats::lm.wfit(f$Z, f$delta, f$w, tol = 1e-12)
  perm <- c(18:10, 1:9)
  summaries <- lapply(rev(f$cohorts), function(rows) or_wls_make(f, rev(rows)))
  for (s in list(whole, or_wls_make(f, perm), or_wls_helpers$or_wls_merge(summaries))) {
    fit <- or_wls_helpers$or_wls_fit(s, rank_tol = 1e-12)
    expect_false(identical(fit$pivot, 1:3))
    or_wls_close(fit$coefficients, direct$coefficients)
    or_wls_close(fit$rss, sum(f$w * direct$residuals^2))
    b <- cbind(c(0.7, -0.2, 0.6), c(1.2, 0.3, -0.4))
    # Normalizer deliberately differs from n and sum(weights).
    N <- 37
    gram <- crossprod(f$Z, f$w * f$Z) / N
    expected <- solve(gram, b) # independent small reference only
    actual <- or_wls_helpers$or_wls_gram_solve(fit, b, normalizer = N)
    or_wls_close(actual, expected)
    or_wls_close(gram %*% actual, b)
    or_wls_close(or_wls_helpers$or_wls_gram_solve(fit, b[, 1L], normalizer = N), expected[, 1L])
    or_wls_close(or_wls_helpers$or_wls_gram_solve(fit, b, normalizer = 1), expected / N)
  }
  f2 <- f
  f2$Z <- f$Z[, c(3, 1, 2)]
  reordered <- or_wls_helpers$or_wls_fit(or_wls_make(f2, 1:18), rank_tol = 1e-12)
  or_wls_close(reordered$coefficients, direct$coefficients[c(3, 1, 2)])
  expect_identical(names(reordered$coefficients), colnames(f2$Z))
})

test_that("bounded weight and column mutations fail the independent objective tolerance", {
  f <- or_wls_fixture()
  rows <- 1:12
  reference <- stats::lm.wfit(f$Z[rows, ], f$delta[rows], f$w[rows], tol = 1e-12)
  changed <- f
  changed$w[1:6] <- changed$w[1:6] * 1.01
  bad_weight <- or_wls_helpers$or_wls_fit(or_wls_make(changed, rows), rank_tol = 1e-12)
  expect_gt(max(abs(bad_weight$coefficients - reference$coefficients)) /
              max(1, abs(reference$coefficients)), or_wls_tol)
  summary <- or_wls_make(f, rows)
  summary$factor[, c(1L, 3L)] <- summary$factor[, c(3L, 1L)]
  beta <- c(0.2, -0.3, 0.7)
  expected <- sum(f$w[rows] * (f$delta[rows] - f$Z[rows, ] %*% beta)^2)
  expect_gt(abs(or_wls_objective(summary, beta) - expected) / max(1, abs(expected)),
            or_wls_tol)
})

test_that("invalid inputs and final deficient designs refuse without changing the model", {
  h <- or_wls_helpers
  f <- or_wls_fixture()
  for (bad in list(0, -1, NA_real_, NaN, Inf)) {
    w <- f$w; w[1L] <- bad
    expect_error(h$or_wls_summary(f$Z, f$delta, w), class = "cdid_or_wls_input_error")
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    Z <- f$Z; Z[1L, 1L] <- bad
    y <- f$delta; y[1L] <- bad
    expect_error(h$or_wls_summary(Z, f$delta, f$w), class = "cdid_or_wls_input_error")
    expect_error(h$or_wls_summary(f$Z, y, f$w), class = "cdid_or_wls_input_error")
  }
  expect_error(h$or_wls_summary(f$Z, f$delta[-1L], f$w), class = "cdid_or_wls_input_error")
  expect_error(h$or_wls_summary(f$Z, f$delta, matrix(f$w)), class = "cdid_or_wls_input_error")
  expect_error(h$or_wls_summary(as.data.frame(f$Z), f$delta, f$w), class = "cdid_or_wls_input_error")
  expect_error(h$or_wls_summary(matrix(1e308, 1L, 1L), 1, 1e308),
               class = "cdid_or_wls_numeric_error")
  s <- or_wls_make(f, 1:18)
  expect_error(h$or_wls_fit(s), class = "cdid_or_wls_input_error")
  for (tol in list(0, -1, 1, NA_real_, Inf, c(1e-12, 1e-10))) {
    expect_error(h$or_wls_fit(s, rank_tol = tol), class = "cdid_or_wls_input_error")
  }
  expect_error(h$or_wls_merge(list()), class = "cdid_or_wls_input_error")
  expect_error(h$or_wls_merge(list(s, list())), class = "cdid_or_wls_input_error")
  f2 <- f; f2$Z <- f$Z[, c(2, 1, 3)]
  expect_error(h$or_wls_merge(list(s, or_wls_make(f2, 1:18))), class = "cdid_or_wls_input_error")
  f2$Z <- cbind(f$Z, duplicate = f$Z[, 2L])
  deficient <- or_wls_make(f2, 1:18)
  snapshot <- deficient
  expect_error(h$or_wls_fit(deficient, rank_tol = 1e-12),
               regexp = "model unchanged", class = "cdid_or_wls_rank_error")
  expect_identical(deficient, snapshot)
  expect_equal(ncol(deficient$factor), 5L)
  fit <- h$or_wls_fit(s, rank_tol = 1e-12)
  expect_error(h$or_wls_gram_solve(fit, c(1, 2, 3)), class = "cdid_or_wls_input_error")
  for (N in list(0, -1, Inf, NA_real_)) {
    expect_error(h$or_wls_gram_solve(fit, c(1, 2, 3), N), class = "cdid_or_wls_input_error")
  }
  for (b in list(c(1, 2), c(1, NA, 3), array(1, c(3, 1, 1)))) {
    expect_error(h$or_wls_gram_solve(fit, b, 37), class = "cdid_or_wls_input_error")
  }
})
