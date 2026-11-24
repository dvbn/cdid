# Streamed within and optimal weighting covariances: replayed signed level
# loadings against dense accumulation references, blocked Grams against dense
# Grams, the streamed forest witness, cross-cohort optimal terms, and the
# workspace refusal before covariance allocation.

stream_test_log <- function(seed = 1L, n_a = 23L, n_cal = 5L, L = 9L) {
  set.seed(seed)
  lg <- level_loading_log(n_cal, n_a, L)
  U <- replicate(n_cal, matrix(0.0, n_a, L), simplify = FALSE)
  pair <- function(tp, tn) if (tn > 0L) c(tp, tn) else tp
  for (r in seq_len(40L)) {
    kind <- c("v", "o", "lr")[1L + r %% 3L]
    tp <- sample.int(n_cal, 1L)
    tn <- if (r %% 4L == 0L && kind == "lr") 0L else sample(setdiff(seq_len(n_cal), tp), 1L)
    rows <- sort(sample.int(n_a, sample.int(n_a, 1L)))
    if (r %% 7L == 0L) rows <- c(rows, rows[1L]) # duplicated index keeps last-write semantics
    nr <- length(rows)
    if (kind == "v") {
      col <- sample.int(L, 1L)
      val <- rnorm(nr)
      loading_log_pair_vec(lg, rows, tp, tn, col, val)
      U[[tp]][rows, col] <- U[[tp]][rows, col] + val
      U[[tn]][rows, col] <- U[[tn]][rows, col] - val
    } else {
      cols <- sort(sample.int(L, sample.int(L, 1L)))
      n_src <- nr + 3L
      src <- sample.int(n_src, nr)
      pre <- if (r %% 2L == 0L) runif(nr) else NULL
      if (kind == "o") {
        a <- rnorm(n_src); b <- rnorm(length(cols))
        loading_log_pair_outer(lg, rows, tp, tn, cols, a, b, src = src, pre = pre)
        W <- outer(a, b)[src, , drop = FALSE]
      } else {
        Lf <- matrix(rnorm(n_src * 3L), n_src, 3L); R <- matrix(rnorm(3L * length(cols)), 3L, length(cols))
        loading_log_lowrank(lg, rows, tp, tn, cols, Lf, R, src = src, pre = pre)
        W <- (Lf %*% R)[src, , drop = FALSE]
      }
      if (!is.null(pre)) W <- pre * W
      U[[tp]][rows, cols] <- U[[tp]][rows, cols, drop = FALSE] + W
      if (tn > 0L) U[[tn]][rows, cols] <- U[[tn]][rows, cols, drop = FALSE] - W
    }
  }
  list(log = lg, U = U)
}

test_that("replayed loadings equal the dense accumulation entry for entry", {
  x <- stream_test_log()
  expect_identical(loading_log_dense(x$log), x$U)
  cp <- seq_len(x$log$L)
  for (b in c(1L, 4L, 10L)) for (tau in seq_len(x$log$n_cal)) {
    starts <- seq.int(1L, x$log$n_a, by = b)
    blk <- do.call(rbind, lapply(starts, function(r0)
      loading_log_replay(x$log, tau, r0, min(x$log$n_a, r0 + b - 1L), cp, length(cp))))
    expect_identical(blk, x$U[[tau]])
  }
})

test_that("blocked within Gram equals the dense Gram for several block sizes", {
  x <- stream_test_log(seed = 2L)
  g <- c(3, 3, 4, 4, 4, 5, 5, 3, 4)
  cols <- c(1L, 2L, 4L, 5L, 7L, 9L)
  U_cols <- lapply(x$U, function(U_t) U_t[, cols, drop = FALSE])
  sc_list <- list(`3` = runif(x$log$n_a), `4` = runif(x$log$n_a), `5` = runif(x$log$n_a))
  for (sc in list(NULL, runif(x$log$n_a), sc_list)) {
    ref <- within_cov_from_level_loadings(U_cols, g[cols], sc)
    expect_identical(within_cov_stream(x$log, cols, g[cols], sc, x$log$n_a), ref)
    for (b in c(1L, 3L, 8L)) {
      V_b <- within_cov_stream(x$log, cols, g[cols], sc, b)
      expect_true(isSymmetric(V_b, tol = 0))
      expect_lt(max(abs(V_b - ref)), 1e-12 * max(1, max(abs(ref))))
      # Within weighting has no cross-cohort blocks.
      expect_true(all(V_b[outer(g[cols], g[cols], "!=")] == 0))
    }
  }
  # The same-pair route scales catalog loadings first; both orders agree.
  scaled <- scale_level_loadings_by_g(lapply(x$U, function(U_t) U_t), g, sc_list)
  ref2 <- within_cov_from_level_loadings(lapply(scaled, function(U_t) U_t[, cols, drop = FALSE]), g[cols])
  expect_lt(max(abs(within_cov_stream(x$log, cols, g[cols], sc_list, 5L) - ref2)),
    1e-12 * max(1, max(abs(ref2))))
})

test_that("streamed forest witness reproduces the dense reduction", {
  n_a <- 30L; n_cal <- 4L
  set.seed(3)
  w1 <- rnorm(n_a); w2 <- rnorm(n_a); w3 <- rnorm(n_a)
  # Links: tree (1,2), (2,3), (3,4); chord (1,3) telescopes exactly; chord (1,4)
  # telescopes in moments but not in loadings on one unit; group two repeats.
  s_vec <- c(1, 2, 3, 1, 1, 1, 2, 1)
  t_vec <- c(2, 3, 4, 3, 4, 2, 3, 3)
  grp <- c("a", "a", "a", "a", "a", "b", "b", "b")
  M <- length(s_vec)
  lg <- level_loading_log(n_cal, n_a, M)
  U <- replicate(n_cal, matrix(0.0, n_a, M), simplify = FALSE)
  add <- function(ell, rows, val) {
    loading_log_pair_vec(lg, rows, t_vec[ell], s_vec[ell], ell, val)
    U[[t_vec[ell]]][rows, ell] <<- U[[t_vec[ell]]][rows, ell] + val
    U[[s_vec[ell]]][rows, ell] <<- U[[s_vec[ell]]][rows, ell] - val
  }
  add(1L, seq_len(n_a), w1); add(2L, seq_len(n_a), w1); add(3L, seq_len(n_a), w1)
  add(4L, seq_len(n_a), w1)
  w_bad <- w1; w_bad[7L] <- w_bad[7L] + 1e-3
  add(5L, seq_len(n_a), w_bad)
  add(6L, 1:20, w2[1:20]); add(7L, 1:20, w2[1:20]); add(8L, 1:20, w2[1:20])
  m <- c(1, 2, 3, 3, 6, 0.5, 0.25, 0.75)
  ref <- verify_post_fit_forest_reduction(s_vec, t_vec, grp, m, U_list_all = U)
  for (b in c(1L, 4L, 7L, n_a)) {
    got <- verify_post_fit_forest_reduction(s_vec, t_vec, grp, m, U_log = lg, block_rows = b)
    expect_identical(got, ref)
  }
  expect_true(4L %in% ref$dropped_idx)
  expect_true(5L %in% ref$kept_idx)
  expect_true(8L %in% ref$dropped_idx)
})

test_that("streamed optimal Gram keeps cross-cohort terms of shared controls", {
  set.seed(4)
  N <- 41L; n_a <- 37L
  Phi_raw <- matrix(rnorm(N * 6L), N, 6L)
  Phi_raw[1:10, 4:6] <- Phi_raw[1:10, 1:3] # shared control rows load both cohorts
  mu <- colMeans(Phi_raw)
  sc <- n_a / N
  rows_fn <- make_centered_phi_rows(Phi_raw, mu, sc)
  Phi_core <- sweep(Phi_raw, 2L, mu, "-") * sc
  expect_identical(rows_fn(seq_len(N)), Phi_core)
  expect_identical(rows_fn(5:9), Phi_core[5:9, , drop = FALSE])
  V_ref <- crossprod(Phi_core) / n_a
  V_ref <- 0.5 * (V_ref + t(V_ref))
  src <- list(wt_mode = "optimal", L = 6L, n_a = n_a, N_orig = N, phi_rows = rows_fn, phi_live_entries = length(Phi_raw))
  expect_identical(cross_pair_general_covariance(src, 1e9, block_rows = N), V_ref)
  for (b in c(1L, 5L, 16L)) {
    V_b <- cross_pair_general_covariance(src, 1e9, block_rows = b)
    expect_lt(max(abs(V_b - V_ref)), 1e-13 * max(abs(V_ref)))
    expect_true(isSymmetric(V_b, tol = 0))
    expect_gt(max(abs(V_b[1:3, 4:6])), 0.1)
  }
  # Budget-derived blocks never fall below the minimum block (or all rows,
  # when fewer remain) and are never refused.
  expect_identical(stream_block_rows(100L, 50, 10, 1e6), min(100L, stream_min_block_rows))
  expect_identical(stream_block_rows(100L, 50, 1e9, 0), 100L)
})

stream_fit_data <- function(n = 160L, Tn = 6L, seed = 7L) {
  set.seed(seed)
  g <- sample(c(0, 3, 4, 5), n, replace = TRUE, prob = c(0.4, 0.2, 0.2, 0.2))
  x <- rbinom(n, 1, 0.5)
  d <- expand.grid(time = seq_len(Tn), id = seq_len(n))
  d$g <- g[d$id]; d$x <- x[d$id]
  d$y <- rnorm(n)[d$id] + 0.2 * d$time + 0.3 * d$x * d$time +
    ifelse(d$g > 0 & d$time >= d$g, 1, 0) + rnorm(nrow(d), sd = 0.5)
  keep <- runif(nrow(d)) < 0.85
  d[keep, ]
}

stream_fit <- function(dat, wt, control = NULL) {
  att_gt(yname = "y", tname = "time", idname = "id", gname = "g", xformla = ~x, data = dat,
    panel = TRUE, allow_unbalanced_panel = TRUE, est_method = "dr", bstrap = FALSE, cband = FALSE,
    mode = "chained", cross_pair = TRUE, outcome_pooling = "pooled", gmm_weighting = wt,
    control = control)
}

test_that("cross-pair within and optimal fits do not depend on the replay block size", {
  dat <- stream_fit_data()
  for (wt in c("within", "optimal")) {
    ref <- stream_fit(dat, wt)
    local({
      local_mocked_bindings(stream_block_rows = function(n_rows, ...) min(as.integer(n_rows), 11L), .package = "cdid")
      blk <- stream_fit(dat, wt)
      expect_identical(blk$chained$support_status, ref$chained$support_status)
      expect_identical(blk$chained$available, ref$chained$available)
      expect_identical(blk$chained$covariance_handling$trigger_reason, ref$chained$covariance_handling$trigger_reason)
      expect_identical(blk$chained$covariance_handling$kept_indices, ref$chained$covariance_handling$kept_indices)
      expect_equal(blk$att, ref$att, tolerance = 1e-10)
      expect_equal(blk$se, ref$se, tolerance = 1e-10)
      expect_equal(blk$inffunc, ref$inffunc, tolerance = 1e-10)
    })
  }
})

test_that("replay blocks keep a minimum size within the stated excess and are never refused", {
  m <- stream_min_block_rows
  expect_gt(m, 1L)
  # Live arrays already exceed the workspace: the minimum block, not a refusal.
  expect_identical(stream_block_rows(1000L, 10, 100, 1e6), m)
  # Fewer rows than the minimum remain: all of them.
  expect_identical(stream_block_rows(m - 5L, 10, 100, 1e6), m - 5L)
  # Budget-derived blocks above the minimum are kept; `cap` still applies.
  expect_identical(stream_block_rows(10000L, 10, 4 * 10000 + 10 * 500, 0), 500L)
  expect_identical(stream_block_rows(10000L, 10, 1e12, 0, cap = 700L), 700L)
  expect_identical(stream_block_rows(10000L, 10, 1e12, 0, cap = 3L), 3L)
  set.seed(25L)
  for (i in seq_len(200L)) {
    n <- sample.int(5000L, 1L); per_row <- runif(1, 1, 400)
    max_e <- runif(1, 0, 1e7); live <- runif(1, 0, 1e7)
    b <- stream_block_rows(n, per_row, max_e, live)
    expect_true(b >= min(n, m) && b <= n)
    spare <- max(0, max_e - live - 4 * n)
    expect_lte(b * per_row, spare + m * per_row)
  }
})

test_that("an infeasible retained solve is refused before the covariance is formed", {
  dat <- stream_fit_data()
  ok <- stream_fit(dat, "within")
  L <- ok$chained$covariance_handling$dimensions$L_retained
  expect_gt(L, 20L)
  for (wt in c("within", "optimal")) {
    local({
      local_mocked_bindings(
        within_cov_stream = function(...) stop("covariance formed before the workspace check"),
        cross_pair_general_covariance = function(...) stop("covariance formed before the workspace check"),
        .package = "cdid")
      # The retained solve needs roughly 48 L^2 entries here; 30 L^2 admits every
      # earlier stage and refuses at the solve's established preflight.
      expect_error(stream_fit(dat, wt, control = list(numeric_workspace_bytes = 8 * 30 * L^2)),
        class = "cdid_gmm_bound_error")
    })
  }
})

test_that("coincident loading channels are summed before the Gram, preserving cancellations", {
  # Two derivative channels write the same unit-date-link entries; link 2 is the
  # exact negative of link 1 and the channels of link 5 cancel. Summing each
  # entry before its cross-product keeps both exact null directions, whereas
  # separate per-channel Grams would not.
  set.seed(9918)
  n_cal <- 5L; n_a <- 13L; M <- 6L
  first <- array(rnorm(n_cal * n_a * M), c(n_cal, n_a, M))
  second <- array(rnorm(n_cal * n_a * M), c(n_cal, n_a, M))
  first[, , 2] <- -first[, , 1]; second[, , 2] <- -second[, , 1]
  first[, , 5] <- -second[, , 5]
  sc <- runif(n_a, 0.5, 2)
  g <- c(1, 1, 1, 2, 2, 2)
  lg <- level_loading_log(n_cal, n_a, M)
  for (tau in seq_len(n_cal)) {
    loading_log_lowrank(lg, seq_len(n_a), tau, 0L, seq_len(M), first[tau, , ], diag(M))
    loading_log_lowrank(lg, seq_len(n_a), tau, 0L, seq_len(M), second[tau, , ], diag(M))
  }
  dense <- lapply(seq_len(n_cal), function(tau) (first[tau, , ] + second[tau, , ]) * sc)
  stacked <- do.call(rbind, dense)
  expected <- crossprod(stacked) * outer(g, g, "==")
  for (b in c(1L, 4L, 13L, 50L)) {
    V <- within_cov_stream(lg, seq_len(M), g, sc, min(b, n_a))
    expect_lt(max(abs(V - expected)), 1e-12 * max(abs(expected)))
    expect_true(all(V[5, ] == 0) && all(V[, 5] == 0))
    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    expect_equal(sum(ev > 1e-10 * max(ev)), 4L)
  }
  separate <- crossprod(do.call(rbind, lapply(seq_len(n_cal), function(tau) first[tau, , ] * sc))) +
    crossprod(do.call(rbind, lapply(seq_len(n_cal), function(tau) second[tau, , ] * sc)))
  expect_gt(separate[5, 5], 0)
})

test_that("blocked optimal Gram is centered in two passes on high-offset inputs", {
  set.seed(77)
  N <- 200L; L <- 5L
  base <- matrix(rnorm(N * L, sd = 1e-3), N, L)
  base[, 5] <- base[, 4] + 1e-9 * rnorm(N) # nearly collinear pair
  Phi_raw <- base + 1e6
  mu <- colMeans(Phi_raw)
  explicit <- sweep(Phi_raw, 2L, mu, "-")
  V_ref <- crossprod(explicit) / N
  V_ref <- 0.5 * (V_ref + t(V_ref))
  src <- list(wt_mode = "optimal", L = L, n_a = N, N_orig = N,
    phi_rows = make_centered_phi_rows(Phi_raw, mu, 1), phi_live_entries = length(Phi_raw))
  e_ref <- eigen(V_ref, symmetric = TRUE, only.values = TRUE)$values
  for (b in c(1L, 7L, 64L, N)) {
    V_b <- cross_pair_general_covariance(src, 1e9, block_rows = b)
    expect_lt(max(abs(V_b - V_ref)), 1e-12 * max(abs(V_ref)))
    e_b <- eigen(V_b, symmetric = TRUE, only.values = TRUE)$values
    expect_lt(max(abs(e_b - e_ref)), 1e-12 * max(e_ref))
  }
})
