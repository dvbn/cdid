# Source-only verification supplies this environment; installed-package tests
# resolve the same private helpers from the candidate namespace.
if (!exists("or_wls_helpers", inherits = TRUE)) {
  or_wls_helpers <- asNamespace("cdid")
}

or_pooled_wls_fixture <- function() {
  # Deterministic 3-cohort panel-like fixture with 3 slope predictors (x, u, v).
  # In cohort A, x is constant (-1); in cohort B, u is constant (0.5);
  # in cohort C, all 3 vary. Overlapping pools and separate blocks are tested.
  D <- rbind(
    # Cohort A (rows 1:6): x is constant (-1)
    cbind(x = rep(-1, 6L),
          u = c(-1, 1, -0.5, 0.5, -1.5, 1.5),
          v = c(0.2, -0.4, 0.8, -0.1, 0.5, -0.7)),
    # Cohort B (rows 7:12): u is constant (0.5)
    cbind(x = c(-1.2, -0.4, 0.3, 0.9, 1.4, -0.6),
          u = rep(0.5, 6L),
          v = c(-0.3, 0.7, -0.5, 0.4, -0.2, 0.9)),
    # Cohort C (rows 13:18): all vary
    cbind(x = c(0.5, 1.5, -0.8, 2.0, -0.2, 1.1),
          u = c(-0.7, 0.4, 1.2, -1.1, 0.8, -0.3),
          v = c(1.0, -0.6, 0.3, -0.9, 0.4, 0.2))
  )
  beta_true <- c(x = 0.6, u = -0.8, v = 1.1)
  block_shift <- rep(c(1.5, -0.7, 2.3), each = 6L)
  noise <- c(0.25, -0.18, 0.11, -0.32, 0.19, -0.04,
             -0.22, 0.31, -0.14, 0.08, -0.27, 0.24,
             0.17, -0.29, 0.36, -0.12, -0.21, 0.09)
  y <- block_shift + drop(D %*% beta_true) + noise
  w <- c(0.9, 1.4, 0.7, 1.8, 1.1, 1.3,
         2.2, 1.7, 2.8, 1.9, 2.5, 1.6,
         0.4, 0.6, 0.5, 0.8, 0.3, 0.7)
  list(
    D = D,
    y = y,
    w = w,
    cohorts = list(A = 1:6, B = 7:12, C = 13:18),
    pools = list(AB = 1:12, BC = 7:18, C_only = 13:18)
  )
}

# Independent reference: direct weighted stacked regression with block dummy
# intercepts and shared slope columns via stats::lm.wfit, plus explicit
# within-block weighted centering for the profiled objective and Gram matrix.
or_pooled_wls_reference <- function(blocks) {
  n_blocks <- length(blocks)
  rows_per_block <- vapply(blocks, function(b) nrow(b$D), integer(1))
  p <- ncol(blocks[[1L]]$D)
  total_n <- sum(rows_per_block)
  E <- matrix(0, nrow = total_n, ncol = n_blocks)
  D_stack <- matrix(0, nrow = total_n, ncol = p)
  D_tilde <- matrix(0, nrow = total_n, ncol = p)
  y_stack <- numeric(total_n)
  y_tilde <- numeric(total_n)
  w_stack <- numeric(total_n)
  offset <- 0L
  for (j in seq_len(n_blocks)) {
    b <- blocks[[j]]
    nb <- nrow(b$D)
    idx <- (offset + 1L):(offset + nb)
    E[idx, j] <- 1
    D_stack[idx, ] <- b$D
    y_stack[idx] <- b$y
    w_stack[idx] <- b$w
    w_sum <- sum(b$w)
    mean_d <- if (p > 0L) colSums(b$D * b$w) / w_sum else numeric(0)
    mean_y <- sum(b$y * b$w) / w_sum
    if (p > 0L) {
      D_tilde[idx, ] <- b$D - rep(mean_d, each = nb)
    }
    y_tilde[idx] <- b$y - mean_y
    offset <- offset + nb
  }
  colnames(D_stack) <- colnames(blocks[[1L]]$D)
  colnames(D_tilde) <- colnames(blocks[[1L]]$D)
  X <- cbind(E, D_stack)
  fit <- stats::lm.wfit(X, y_stack, w_stack, tol = 1e-12)
  intercepts <- fit$coefficients[seq_len(n_blocks)]
  names(intercepts) <- names(blocks)
  slopes <- if (p > 0L) {
    fit$coefficients[n_blocks + seq_len(p)]
  } else {
    numeric(0)
  }
  names(slopes) <- colnames(blocks[[1L]]$D)
  list(
    fit = fit,
    intercepts = intercepts,
    coefficients = slopes,
    rss = sum(w_stack * fit$residuals^2),
    fitted = fit$fitted.values,
    E = E,
    D_stack = D_stack,
    D_tilde = D_tilde,
    y_stack = y_stack,
    y_tilde = y_tilde,
    w_stack = w_stack
  )
}

or_wls_tol <- 1e-10
or_wls_close <- function(actual, expected) {
  expect_lt(max(c(0, abs(actual - expected))) /
              max(1, abs(expected)), or_wls_tol)
}

test_that("shared-slope pooled WLS matches direct stacked WLS with overlapping and unequal-weight blocks", {
  h <- or_wls_helpers
  f <- or_pooled_wls_fixture()
  # Overlapping pools AB, BC, and C_only with distinct block weight scales.
  raw_blocks <- list(
    AB = list(D = f$D[f$pools$AB, , drop = FALSE],
              y = f$y[f$pools$AB],
              w = f$w[f$pools$AB] * 0.75),
    BC = list(D = f$D[f$pools$BC, , drop = FALSE],
              y = f$y[f$pools$BC] + 0.4,
              w = f$w[f$pools$BC] * 2.10),
    C_only = list(D = f$D[f$pools$C_only, , drop = FALSE],
                  y = f$y[f$pools$C_only] - 0.9,
                  w = f$w[f$pools$C_only] * 0.35)
  )
  summaries <- lapply(raw_blocks, function(b) {
    h$or_pooled_wls_block(b$D, b$y, b$w)
  })
  merged <- h$or_pooled_wls_merge(summaries)
  fit <- h$or_pooled_wls_fit(merged, rank_tol = 1e-12)
  ref <- or_pooled_wls_reference(raw_blocks)

  expect_equal(ref$fit$rank, 6L)
  expect_equal(fit$rank, 3L)
  expect_equal(fit$n, 30L)
  expect_equal(fit$n_blocks, 3L)
  or_wls_close(fit$coefficients, ref$coefficients)
  or_wls_close(fit$intercepts, ref$intercepts)
  or_wls_close(fit$rss, ref$rss)

  # Factor-size accounting: residual factor is at most (p + 1) x (p + 1).
  expect_identical(fit$accounting$residual_dim, c(4L, 4L))
  expect_lte(fit$accounting$max_augmented_rows, 5L)
  expect_lte(fit$accounting$max_residual_rows, 4L)
  expect_lte(fit$accounting$max_residual_merge_rows, 8L)

  # Crossproduct of merged residual factor equals exact centered crossproduct.
  centered_aug <- cbind(ref$D_tilde, ref$y_tilde) * sqrt(ref$w_stack)
  or_wls_close(crossprod(merged$residual$factor), crossprod(centered_aug))

  # Profiled and full objectives match at arbitrary slopes and intercepts.
  beta_probes <- rbind(
    c(0, 0, 0),
    c(0.6, -0.8, 1.1),
    c(-0.5, 1.2, -0.4),
    fit$coefficients
  )
  alpha_probes <- rbind(
    c(0, 0, 0),
    c(1.0, -0.5, 2.0),
    fit$intercepts
  )
  for (beta in split(beta_probes, row(beta_probes))) {
    expected_prof <- sum(ref$w_stack * (ref$y_tilde - drop(ref$D_tilde %*% beta))^2)
    or_wls_close(h$or_pooled_wls_objective(merged, beta), expected_prof)
    for (alpha in split(alpha_probes, row(alpha_probes))) {
      pred <- drop(ref$E %*% alpha + ref$D_stack %*% beta)
      expected_joint <- sum(ref$w_stack * (ref$y_stack - pred)^2)
      or_wls_close(h$or_pooled_wls_objective(merged, beta, intercepts = alpha),
                   expected_joint)
    }
  }
  or_wls_close(h$or_pooled_wls_objective(merged, fit$coefficients), fit$rss)
  or_wls_close(h$or_pooled_wls_objective(merged, fit$coefficients,
                                         intercepts = fit$intercepts), fit$rss)

  # In-sample and target predictions match direct stacked WLS.
  block_ids <- rep(names(raw_blocks), vapply(raw_blocks, function(b) nrow(b$D), integer(1)))
  or_wls_close(h$or_pooled_wls_predict(fit, ref$D_stack, block = block_ids),
               ref$fitted)
  target_D <- rbind(c(0.25, -0.4, 0.8), c(-1.1, 0.9, -0.2))
  for (j in seq_along(raw_blocks)) {
    expected_target <- ref$intercepts[j] + drop(target_D %*% ref$coefficients)
    or_wls_close(h$or_pooled_wls_predict(fit, target_D, block = j),
                 expected_target)
    or_wls_close(h$or_pooled_wls_predict(fit, target_D, block = names(raw_blocks)[j]),
                 expected_target)
  }

  # Profiled-Gram loading solve matches centered Gram inverse AND Schur complement
  # of the full block+slope Gram matrix.
  N <- 43
  b_rhs <- cbind(c(0.7, -0.3, 1.1), c(-0.4, 0.9, 0.2))
  rownames(b_rhs) <- colnames(f$D)
  colnames(b_rhs) <- c("g1", "g2")
  gram_prof <- crossprod(ref$D_tilde, ref$w_stack * ref$D_tilde) / N
  expected_load <- solve(gram_prof, b_rhs)
  X_full <- cbind(ref$E, ref$D_stack)
  gram_full <- crossprod(X_full, ref$w_stack * X_full) / N
  expected_schur <- solve(gram_full, rbind(matrix(0, 3L, 2L), b_rhs))[4:6, , drop = FALSE]
  actual_load <- h$or_pooled_wls_gram_solve(fit, b_rhs, normalizer = N)
  or_wls_close(actual_load, expected_load)
  or_wls_close(actual_load, expected_schur)
  or_wls_close(gram_prof %*% actual_load, b_rhs)
  or_wls_close(h$or_pooled_wls_gram_solve(fit, b_rhs[, 1L], normalizer = N),
               expected_load[, 1L])
})

test_that("locally rank-deficient, one-row and constant-X blocks identify shared slopes and local intercepts", {
  h <- or_wls_helpers
  f <- or_pooled_wls_fixture()
  # Block A has constant x; Block B has constant u; Block single has 1 row;
  # Block const_all has 4 rows with all predictors constant but varying y.
  D_const <- matrix(rep(c(2.5, -1.0, 0.4), each = 4L), nrow = 4L,
                    dimnames = list(NULL, colnames(f$D)))
  y_const <- c(1.2, -0.4, 0.8, 2.1)
  w_const <- c(0.6, 1.3, 0.9, 1.7)

  raw_blocks <- list(
    A_x_const = list(D = f$D[f$cohorts$A, , drop = FALSE],
                     y = f$y[f$cohorts$A],
                     w = f$w[f$cohorts$A]),
    B_u_const = list(D = f$D[f$cohorts$B, , drop = FALSE],
                     y = f$y[f$cohorts$B],
                     w = f$w[f$cohorts$B]),
    single_row = list(D = f$D[13L, , drop = FALSE],
                      y = f$y[13L],
                      w = f$w[13L]),
    const_all = list(D = D_const, y = y_const, w = w_const)
  )

  # Every single block is locally rank-deficient for the 3 slopes on its own.
  summaries <- lapply(raw_blocks, function(b) {
    s <- h$or_pooled_wls_block(b$D, b$y, b$w)
    expect_error(h$or_pooled_wls_fit(s, rank_tol = 1e-12),
                 regexp = "model unchanged",
                 class = "cdid_or_wls_rank_error")
    s
  })
  # The 1-row block has 0 residual rows after intercept partialling.
  expect_identical(dim(summaries$single_row$residual$factor), c(0L, 4L))
  # The constant-X block has nonzero y-residual constant in its residual factor.
  expect_gt(summaries$const_all$residual$factor[1L, 4L]^2, 0)

  # Combined across blocks, the centered design has full rank 3.
  merged <- h$or_pooled_wls_merge(summaries)
  fit <- h$or_pooled_wls_fit(merged, rank_tol = 1e-12)
  ref <- or_pooled_wls_reference(raw_blocks)

  expect_equal(fit$rank, 3L)
  or_wls_close(fit$coefficients, ref$coefficients)
  or_wls_close(fit$intercepts, ref$intercepts)
  or_wls_close(fit$rss, ref$rss)

  # 1-row block residual is identically zero: alpha_single + d_13' beta == y_13.
  or_wls_close(h$or_pooled_wls_predict(fit, raw_blocks$single_row$D,
                                       block = "single_row"),
               raw_blocks$single_row$y)

  # Target predictions for all 4 blocks (including 1-row and constant-X blocks)
  # match direct stacked regression.
  target_D <- rbind(c(0.5, -0.5, 0.2), c(-0.8, 1.3, -0.6))
  for (nm in names(raw_blocks)) {
    or_wls_close(h$or_pooled_wls_predict(fit, target_D, block = nm),
                 ref$intercepts[[nm]] + drop(target_D %*% ref$coefficients))
  }
})

test_that("chunked within-block compression, nested merges, permutations and nonidentity pivots preserve all results", {
  h <- or_wls_helpers
  f <- or_pooled_wls_fixture()
  # Scale middle predictor so LAPACK QR pivots non-trivially.
  f$D[, 2L] <- 7 * f$D[, 2L]

  raw_blocks <- list(
    A = list(D = f$D[1:6, , drop = FALSE], y = f$y[1:6], w = f$w[1:6]),
    B = list(D = f$D[7:12, , drop = FALSE], y = f$y[7:12], w = f$w[7:12]),
    C = list(D = f$D[13:18, , drop = FALSE], y = f$y[13:18], w = f$w[13:18])
  )
  ref <- or_pooled_wls_reference(raw_blocks)

  # Compress block A and block B in uneven chunks including empty chunks.
  chunks_A <- lapply(list(integer(), 1L, 2:3, integer(), 4:6), function(idx) {
    h$or_pooled_wls_block_summary(raw_blocks$A$D[idx, , drop = FALSE],
                                  raw_blocks$A$y[idx],
                                  raw_blocks$A$w[idx])
  })
  expect_identical(dim(chunks_A[[1L]]$raw$factor), c(0L, 5L))
  expect_identical(dim(chunks_A[[2L]]$raw$factor), c(1L, 5L))
  block_A_chunked <- h$or_pooled_wls_partial_block(
    h$or_pooled_wls_merge_chunks(chunks_A)
  )

  chunks_B <- lapply(list(c(6L, 4L), c(2L, 5L, 1L), 3L, integer()), function(idx) {
    h$or_pooled_wls_block_summary(raw_blocks$B$D[idx, , drop = FALSE],
                                  raw_blocks$B$y[idx],
                                  raw_blocks$B$w[idx])
  })
  block_B_chunked <- h$or_pooled_wls_partial_block(
    h$or_pooled_wls_merge_chunks(chunks_B)
  )
  block_C_perm <- h$or_pooled_wls_block(raw_blocks$C$D[6:1, , drop = FALSE],
                                        raw_blocks$C$y[6:1],
                                        raw_blocks$C$w[6:1])

  # Nested tree merge in reversed block order C, B, A.
  nested <- h$or_pooled_wls_merge(list(
    h$or_pooled_wls_merge(list(C = block_C_perm, B = block_B_chunked)),
    A = block_A_chunked
  ))

  fit_nested <- h$or_pooled_wls_fit(nested, rank_tol = 1e-12)
  expect_false(identical(fit_nested$pivot, 1:3))
  or_wls_close(fit_nested$coefficients, ref$coefficients)
  or_wls_close(fit_nested$intercepts[c("A", "B", "C")], ref$intercepts)
  or_wls_close(fit_nested$rss, ref$rss)

  b_rhs <- cbind(c(0.5, -0.4, 0.9), c(-1.1, 0.3, 0.7))
  N <- 19
  gram_prof <- crossprod(ref$D_tilde, ref$w_stack * ref$D_tilde) / N
  expected_load <- solve(gram_prof, b_rhs)
  or_wls_close(h$or_pooled_wls_gram_solve(fit_nested, b_rhs, normalizer = N),
               expected_load)
  or_wls_close(h$or_pooled_wls_gram_solve(fit_nested, b_rhs, normalizer = 1),
               expected_load / N)

  # Column permutation equivariance.
  col_perm <- c(3L, 1L, 2L)
  perm_summaries <- lapply(raw_blocks, function(b) {
    h$or_pooled_wls_block(b$D[, col_perm, drop = FALSE], b$y, b$w)
  })
  fit_perm <- h$or_pooled_wls_fit(h$or_pooled_wls_merge(perm_summaries),
                                  rank_tol = 1e-12)
  or_wls_close(fit_perm$coefficients, ref$coefficients[col_perm])
  expect_identical(names(fit_perm$coefficients), colnames(f$D)[col_perm])
  or_wls_close(fit_perm$intercepts, ref$intercepts)
  or_wls_close(fit_perm$rss, ref$rss)
  or_wls_close(
    h$or_pooled_wls_gram_solve(fit_perm, b_rhs[col_perm, , drop = FALSE],
                               normalizer = N),
    expected_load[col_perm, , drop = FALSE]
  )
})

test_that("zero non-intercept predictors work as ordinary intercept-only blocks", {
  h <- or_wls_helpers
  f <- or_pooled_wls_fixture()
  D0_A <- matrix(numeric(0), nrow = 6L, ncol = 0L)
  D0_B <- matrix(numeric(0), nrow = 1L, ncol = 0L)
  D0_C <- matrix(numeric(0), nrow = 6L, ncol = 0L)
  raw_blocks <- list(
    A = list(D = D0_A, y = f$y[1:6], w = f$w[1:6]),
    B_one = list(D = D0_B, y = f$y[7L], w = f$w[7L]),
    C = list(D = D0_C, y = f$y[13:18], w = f$w[13:18])
  )
  # Chunk block A even with p = 0.
  c1 <- h$or_pooled_wls_block_summary(D0_A[1:2, , drop = FALSE], f$y[1:2], f$w[1:2])
  c2 <- h$or_pooled_wls_block_summary(D0_A[3:6, , drop = FALSE], f$y[3:6], f$w[3:6])
  sA <- h$or_pooled_wls_partial_block(h$or_pooled_wls_merge_chunks(list(c1, c2)))
  sB <- h$or_pooled_wls_block(D0_B, f$y[7L], f$w[7L])
  sC <- h$or_pooled_wls_block(D0_C, f$y[13:18], f$w[13:18])

  merged <- h$or_pooled_wls_merge(list(A = sA, B_one = sB, C = sC))
  fit <- h$or_pooled_wls_fit(merged, rank_tol = 1e-12)
  ref <- or_pooled_wls_reference(raw_blocks)

  expect_identical(fit$coefficients, numeric(0))
  expect_identical(fit$rank, 0L)
  expect_identical(dim(fit$R), c(0L, 0L))
  expect_identical(fit$pivot, integer(0))
  or_wls_close(fit$intercepts, ref$intercepts)
  or_wls_close(fit$rss, ref$rss)
  or_wls_close(h$or_pooled_wls_objective(merged, numeric(0)), ref$rss)
  or_wls_close(h$or_pooled_wls_objective(merged, numeric(0),
                                         intercepts = c(0, 0, 0)),
               sum(ref$w_stack * ref$y_stack^2))
  or_wls_close(h$or_pooled_wls_predict(fit, matrix(numeric(0), 3L, 0L),
                                       block = c("A", "B_one", "C")),
               unname(ref$intercepts))
  expect_identical(h$or_pooled_wls_gram_solve(fit, numeric(0), normalizer = 10),
                   numeric(0))
  expect_identical(
    dim(h$or_pooled_wls_gram_solve(fit, matrix(numeric(0), 0L, 2L), normalizer = 10)),
    c(0L, 2L)
  )
})

test_that("global centered deficiency is diagnosed even when raw global X has full rank and block means differ", {
  h <- or_wls_helpers
  f <- or_pooled_wls_fixture()

  # Case 1: Single slope predictor x that is constant within each block (-1, 0, 2)
  # so block means differ and raw global [1, x] has full rank 2, but centered x == 0.
  blocks_1d <- list(
    A = h$or_pooled_wls_block(matrix(-1, 6L, 1L, dimnames = list(NULL, "x")),
                              f$y[1:6], f$w[1:6]),
    B = h$or_pooled_wls_block(matrix(0, 6L, 1L, dimnames = list(NULL, "x")),
                              f$y[7:12], f$w[7:12]),
    C = h$or_pooled_wls_block(matrix(2, 6L, 1L, dimnames = list(NULL, "x")),
                              f$y[13:18], f$w[13:18])
  )
  raw_X_1d <- cbind(1, rep(c(-1, 0, 2), each = 6L))
  expect_equal(qr(raw_X_1d)$rank, 2L)
  merged_1d <- h$or_pooled_wls_merge(blocks_1d)
  snapshot_1d <- merged_1d
  expect_error(h$or_pooled_wls_fit(merged_1d, rank_tol = 1e-12),
               regexp = "model unchanged",
               class = "cdid_or_wls_rank_error")
  expect_identical(merged_1d, snapshot_1d)

  # Case 2: Two slope predictors (d1, d2) where d2 = 2 * d1 + c_b with distinct
  # block shifts c_b = (-3, 1, 4). Raw global [1, d1, d2] has full column rank 3,
  # but within every block (d2 - mean_b(d2)) == 2 * (d1 - mean_b(d1)), so the
  # profiled centered design has rank 1 < 2.
  d1 <- f$D[, "v"]
  c_shift <- rep(c(-3, 1, 4), each = 6L)
  d2 <- 2 * d1 + c_shift
  D_collinear_centered <- cbind(d1 = d1, d2 = d2)
  expect_equal(qr(cbind(1, D_collinear_centered))$rank, 3L)
  blocks_2d <- lapply(f$cohorts, function(idx) {
    h$or_pooled_wls_block(D_collinear_centered[idx, , drop = FALSE],
                          f$y[idx], f$w[idx])
  })
  merged_2d <- h$or_pooled_wls_merge(blocks_2d)
  snapshot_2d <- merged_2d
  expect_error(h$or_pooled_wls_fit(merged_2d, rank_tol = 1e-12),
               regexp = "model unchanged",
               class = "cdid_or_wls_rank_error")
  expect_identical(merged_2d, snapshot_2d)

  # Case 3: Every block has only 1 row (zero total centered degrees of freedom).
  one_row_blocks <- list(
    A = h$or_pooled_wls_block(f$D[1L, , drop = FALSE], f$y[1L], f$w[1L]),
    B = h$or_pooled_wls_block(f$D[7L, , drop = FALSE], f$y[7L], f$w[7L])
  )
  expect_error(h$or_pooled_wls_fit(h$or_pooled_wls_merge(one_row_blocks),
                                  rank_tol = 1e-12),
               regexp = "model unchanged",
               class = "cdid_or_wls_rank_error")
})

test_that("deliberate perturbations of block weight normalization and local slope averaging fail reference checks", {
  h <- or_wls_helpers
  f <- or_pooled_wls_fixture()
  # Construct 3 full-rank blocks with unequal weight scales and heterogeneous
  # covariate dispersions so both block-normalized WLS and local-slope averaging
  # are well-defined and deliberately wrong.
  raw_blocks <- list(
    B1 = list(D = f$D[1:6, , drop = FALSE] +
                cbind(c(-0.5, 0.4, -0.2, 0.6, -0.3, 0.0), 0, 0),
              y = f$y[1:6],
              w = f$w[1:6] * 0.4),
    B2 = list(D = f$D[7:12, , drop = FALSE] +
                cbind(0, c(0.3, -0.5, 0.2, -0.4, 0.6, -0.2), 0),
              y = f$y[7:12],
              w = f$w[7:12] * 3.5),
    B3 = list(D = f$D[13:18, , drop = FALSE],
              y = f$y[13:18],
              w = f$w[13:18] * 1.1)
  )
  ref <- or_pooled_wls_reference(raw_blocks)
  merged <- h$or_pooled_wls_merge(lapply(raw_blocks, function(b) {
    h$or_pooled_wls_block(b$D, b$y, b$w)
  }))
  fit <- h$or_pooled_wls_fit(merged, rank_tol = 1e-12)
  or_wls_close(fit$coefficients, ref$coefficients)

  # Deliberate perturbation 1: normalize each block's weights to sum to 1.
  norm_blocks <- lapply(raw_blocks, function(b) {
    h$or_pooled_wls_block(b$D, b$y, b$w / sum(b$w))
  })
  fit_norm <- h$or_pooled_wls_fit(h$or_pooled_wls_merge(norm_blocks),
                                  rank_tol = 1e-12)
  expect_gt(max(abs(fit_norm$coefficients - ref$coefficients)) /
              max(1, abs(ref$coefficients)), or_wls_tol)
  expect_gt(abs(h$or_pooled_wls_objective(merged, fit_norm$coefficients) - ref$rss) /
              max(1, abs(ref$rss)), or_wls_tol)

  # Deliberate perturbation 2: estimate separate local slopes per block and
  # average them (unweighted or weighted by block total weight).
  local_slopes <- lapply(raw_blocks, function(b) {
    stats::lm.wfit(cbind(1, b$D), b$y, b$w, tol = 1e-12)$coefficients[-1L]
  })
  beta_avg_equal <- Reduce(`+`, local_slopes) / length(local_slopes)
  block_w <- vapply(raw_blocks, function(b) sum(b$w), numeric(1))
  beta_avg_weighted <- drop(do.call(cbind, local_slopes) %*% (block_w / sum(block_w)))
  expect_gt(max(abs(beta_avg_equal - ref$coefficients)) /
              max(1, abs(ref$coefficients)), or_wls_tol)
  expect_gt(max(abs(beta_avg_weighted - ref$coefficients)) /
              max(1, abs(ref$coefficients)), or_wls_tol)
  expect_gt(abs(h$or_pooled_wls_objective(merged, beta_avg_equal) - ref$rss) /
              max(1, abs(ref$rss)), or_wls_tol)
  expect_gt(abs(h$or_pooled_wls_objective(merged, beta_avg_weighted) - ref$rss) /
              max(1, abs(ref$rss)), or_wls_tol)

  # Deliberate perturbation 3: 1% weight perturbation in one block and residual
  # factor column swap.
  perturbed <- raw_blocks
  perturbed$B1$w <- perturbed$B1$w * 1.01
  fit_pert <- h$or_pooled_wls_fit(
    h$or_pooled_wls_merge(lapply(perturbed, function(b) h$or_pooled_wls_block(b$D, b$y, b$w))),
    rank_tol = 1e-12
  )
  expect_gt(max(abs(fit_pert$coefficients - ref$coefficients)) /
              max(1, abs(ref$coefficients)), or_wls_tol)
  swapped <- merged
  swapped$residual$factor[, c(1L, 3L)] <- swapped$residual$factor[, c(3L, 1L)]
  beta_test <- c(0.4, -0.5, 0.8)
  expected_obj <- sum(ref$w_stack * (ref$y_tilde - drop(ref$D_tilde %*% beta_test))^2)
  expect_gt(abs(h$or_pooled_wls_objective(swapped, beta_test) - expected_obj) /
              max(1, abs(expected_obj)), or_wls_tol)
})

test_that("invalid inputs, nonfinite values and empty-intercept blocks abort with typed errors", {
  h <- or_wls_helpers
  f <- or_pooled_wls_fixture()
  D <- f$D[1:6, , drop = FALSE]
  y <- f$y[1:6]
  w <- f$w[1:6]

  for (bad in list(0, -1, NA_real_, NaN, Inf)) {
    w_bad <- w; w_bad[1L] <- bad
    expect_error(h$or_pooled_wls_block(D, y, w_bad), class = "cdid_or_wls_input_error")
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    D_bad <- D; D_bad[1L, 1L] <- bad
    y_bad <- y; y_bad[1L] <- bad
    expect_error(h$or_pooled_wls_block(D_bad, y, w), class = "cdid_or_wls_input_error")
    expect_error(h$or_pooled_wls_block(D, y_bad, w), class = "cdid_or_wls_input_error")
  }
  expect_error(h$or_pooled_wls_block(D, y[-1L], w), class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_block(D, y, matrix(w)), class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_block(as.data.frame(D), y, w), class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_block(matrix(1e308, 1L, 1L), 1, 1e308),
               class = "cdid_or_wls_numeric_error")

  # Empty intercept block diagnostics (direct 0-row block and merged empty chunks).
  expect_error(h$or_pooled_wls_block(D[integer(), , drop = FALSE], numeric(0), numeric(0)),
               regexp = "model unchanged",
               class = "cdid_or_wls_rank_error")
  empty_chunk <- h$or_pooled_wls_block_summary(D[integer(), , drop = FALSE],
                                               numeric(0), numeric(0))
  expect_error(
    h$or_pooled_wls_partial_block(
      h$or_pooled_wls_merge_chunks(list(empty_chunk, empty_chunk))
    ),
    regexp = "model unchanged",
    class = "cdid_or_wls_rank_error"
  )

  # Merge and fit input validation.
  sA <- h$or_pooled_wls_block(f$D[1:6, , drop = FALSE], f$y[1:6], f$w[1:6])
  sB <- h$or_pooled_wls_block(f$D[7:12, , drop = FALSE], f$y[7:12], f$w[7:12])
  sC <- h$or_pooled_wls_block(f$D[13:18, , drop = FALSE], f$y[13:18], f$w[13:18])
  merged <- h$or_pooled_wls_merge(list(A = sA, B = sB, C = sC))

  expect_error(h$or_pooled_wls_merge_chunks(list()), class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_merge(list()), class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_merge(list(sA, list())), class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_merge(list(A = sA, A = sB)), class = "cdid_or_wls_input_error")

  sB_reordered <- h$or_pooled_wls_block(f$D[7:12, c(2, 1, 3), drop = FALSE],
                                        f$y[7:12], f$w[7:12])
  expect_error(h$or_pooled_wls_merge(list(sA, sB_reordered)),
               class = "cdid_or_wls_input_error")

  expect_error(h$or_pooled_wls_fit(merged), class = "cdid_or_wls_input_error")
  for (tol in list(0, -1, 1, NA_real_, Inf, c(1e-12, 1e-10))) {
    expect_error(h$or_pooled_wls_fit(merged, rank_tol = tol),
                 class = "cdid_or_wls_input_error")
  }

  fit <- h$or_pooled_wls_fit(merged, rank_tol = 1e-12)
  expect_error(h$or_pooled_wls_gram_solve(fit, c(1, 2, 3)),
               class = "cdid_or_wls_input_error")
  for (N in list(0, -1, Inf, NA_real_)) {
    expect_error(h$or_pooled_wls_gram_solve(fit, c(1, 2, 3), N),
                 class = "cdid_or_wls_input_error")
  }
  for (b_bad in list(c(1, 2), c(1, NA, 3), array(1, c(3, 1, 1)))) {
    expect_error(h$or_pooled_wls_gram_solve(fit, b_bad, 37),
                 class = "cdid_or_wls_input_error")
  }
  expect_error(h$or_pooled_wls_objective(merged, c(1, 2)),
               class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_objective(merged, c(1, 2, 3), intercepts = c(1, 2)),
               class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_predict(fit, D[, 1:2, drop = FALSE], block = "A"),
               class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_predict(fit, D, block = "UNKNOWN"),
               class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_predict(fit, D, block = 4L),
               class = "cdid_or_wls_input_error")
  expect_error(h$or_pooled_wls_predict(fit, D, block = c("A", "B")),
               class = "cdid_or_wls_input_error")
})

test_that("exact original-row constancy survives chunks without removing shared columns", {
  h <- or_wls_helpers
  D <- matrix(c(1, 1, 3, 3), ncol = 1L, dimnames = list(NULL, "x"))
  y <- c(1, 2, 4, 6)
  w <- c(.7, 1.3, .9, 1.6)
  empty <- h$or_pooled_wls_block_summary(D[integer(), , drop = FALSE], numeric(), numeric())
  a <- h$or_pooled_wls_block_summary(D[1:2, , drop = FALSE], y[1:2], w[1:2])
  b <- h$or_pooled_wls_block_summary(D[3:4, , drop = FALSE], y[3:4], w[3:4])
  expect_identical(a$constant, TRUE)
  expect_identical(b$constant, TRUE)
  expect_identical(empty$constant, TRUE)
  expect_identical(empty$anchor, numeric())
  # Each chunk is constant, but the same block has genuine within variation.
  merged <- h$or_pooled_wls_merge_chunks(list(empty, a, empty, b, empty))
  reversed <- h$or_pooled_wls_merge_chunks(list(empty, b, empty, a, empty))
  expect_identical(merged$constant, FALSE)
  expect_identical(reversed$constant, FALSE)
  expect_identical(merged$anchor, 1)
  expect_identical(reversed$anchor, 3)
  fit <- h$or_pooled_wls_fit(h$or_pooled_wls_partial_block(merged), rank_tol = 1e-12)
  rev_fit <- h$or_pooled_wls_fit(h$or_pooled_wls_partial_block(reversed), rank_tol = 1e-12)
  ref <- stats::lm.wfit(cbind(1, D), y, w, tol = 1e-12)
  or_wls_close(fit$coefficients, ref$coefficients[2L])
  or_wls_close(fit$intercepts, ref$coefficients[1L])
  or_wls_close(rev_fit$coefficients, fit$coefficients)

  # Exact constant support is zero, not a small QR projection residual.
  constant_a <- h$or_pooled_wls_partial_block(
    h$or_pooled_wls_merge_chunks(list(empty, a, empty, a))
  )
  constant_b <- h$or_pooled_wls_partial_block(b)
  expect_equal(constant_a$residual$factor[, 1L], rep(0, nrow(constant_a$residual$factor)))
  expect_equal(constant_b$residual$factor[, 1L], rep(0, nrow(constant_b$residual$factor)))
  no_support <- h$or_pooled_wls_merge(list(A = constant_a, B = constant_b))
  expect_error(h$or_pooled_wls_fit(no_support, rank_tol = 1e-12),
               class = "cdid_or_wls_rank_error", regexp = "model unchanged")
  # A locally constant block retains the column for another informative block.
  pooled <- h$or_pooled_wls_merge(list(informative = h$or_pooled_wls_partial_block(merged),
                                      constant = constant_a))
  pooled_fit <- h$or_pooled_wls_fit(pooled, rank_tol = 1e-12)
  expect_identical(pooled_fit$p, 1L)
  or_wls_close(pooled_fit$coefficients, fit$coefficients)
})

test_that("final centered rank convention is unchanged by blockwise covariate translation", {
  h <- or_wls_helpers
  rank_tol <- 1e-4
  x <- c(-2, -1, 0, 1, 2) * 1e-5
  y <- x + c(.1, -.2, .05, .1, -.05) * 1e-5
  w <- c(.7, 1.2, .9, 1.5, 1.1)
  fits <- lapply(c(0, 1), function(shift) {
    D <- matrix(x + shift, ncol = 1L, dimnames = list(NULL, "x"))
    summary <- h$or_pooled_wls_block(D, y, w)
    f <- h$or_pooled_wls_fit(summary, rank_tol = rank_tol)
    inherited <- h$or_wls_fit(summary$residual, rank_tol = rank_tol)
    expect_identical(f$rank, 1L)
    expect_identical(f$cutoff, inherited$cutoff)
    or_wls_close(f$coefficients, inherited$coefficients)
    f
  })
  or_wls_close(fits[[1L]]$coefficients, fits[[2L]]$coefficients)
  or_wls_close(fits[[1L]]$rss, fits[[2L]]$rss)
})
