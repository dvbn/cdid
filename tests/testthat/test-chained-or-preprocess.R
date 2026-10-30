# Portable pure numerical checks; the source harness supplies an environment.
# Installed tests resolve the actual namespace. No fitted expected values use
# production output: objectives/predictions use direct lm.wfit linear algebra.
if (!exists("or_preprocess_helpers", inherits = TRUE)) or_preprocess_helpers <- asNamespace("cdid")

testthat::test_that("common QR aliases preserve objectives, predictions and differentiated fits", {
  h <- or_preprocess_helpers
  x <- rep(c(-2, -1, 0, 1, 2), 4)
  u <- rep(c(-1, 1), 10)
  X <- cbind(intercept = 1, x = x, u = u, duplicate = x, combination = 2 + 3*x - u)
  y <- 1 + x/3 - u/2 + (seq_along(x) %% 7)/9
  w <- 1 + (seq_along(x) %% 3)/4
  b <- h$or_public_basis(X, rep(TRUE, length(x)))
  testthat::expect_identical(b$record$retained_columns, c("intercept", "x", "u"))
  testthat::expect_identical(b$record$aliased_columns, c("duplicate", "combination"))
  testthat::expect_equal(b$Z %*% b$record$alias_map, X, tolerance = 1e-13)
  # Fixed equivalent basis, including a nonsingular change of reduced basis.
  T <- matrix(c(1, .2, 0, 0, 2, .1, .3, 0, 1), 3)
  fit <- function(Z, weights) stats::lm.wfit(Z, y, weights, tol = 1e-7)
  f <- fit(X, w); r <- fit(b$Z, w); z <- fit(b$Z %*% T, w)
  testthat::expect_equal(f$residuals, r$residuals, tolerance = 1e-12)
  testthat::expect_equal(f$residuals, z$residuals, tolerance = 1e-12)
  beta <- c(.3, -.1, .2, -.4, .7)
  testthat::expect_equal(sum(w*(y-X%*%beta)^2), sum(w*(y-b$Z%*%(b$record$alias_map%*%beta))^2), tolerance = 1e-12)
  # A direct analytic WLS case-weight derivative of prediction at every row.
  independent_if <- function(Z, f) {
    G <- crossprod(Z, w*Z)
    nrow(Z) * Z %*% solve(G, t((w*f$residuals)*Z))
  }
  testthat::expect_equal(independent_if(b$Z, r), independent_if(b$Z %*% T, z), tolerance = 1e-11)
  testthat::expect_true(b$record$frozen_during_differentiation)
  testthat::expect_false(b$record$numerical_certificate)
})

testthat::test_that("nontrivial near dependence and treated-only support are not erased", {
  h <- or_preprocess_helpers
  x <- rep(c(-2, -1, 0, 1, 2), 4)
  near <- cbind(intercept = 1, x = x, near = x + 1e-9*(seq_along(x) %% 3 - 1))
  b <- h$or_public_basis(near, rep(TRUE, length(x)))
  testthat::expect_equal(ncol(b$Z), 3L)
  testthat::expect_length(b$record$aliased_columns, 0L)
  testthat::expect_identical(b$record$nontrivial_candidates_retained, "near")
  mixed <- h$or_public_basis(cbind(near, duplicate = x), rep(TRUE, length(x)))
  testthat::expect_identical(mixed$record$aliased_columns, "duplicate")
  testthat::expect_identical(mixed$record$retained_columns, c("intercept", "x", "near"))
  testthat::expect_equal(mixed$Z %*% mixed$record$alias_map, cbind(near, duplicate = x), tolerance = 1e-13)
  testthat::expect_gt(b$record$reconstruction_max_scaled["near"], b$record$reconstruction_roundoff_threshold)
  # Common-row rank is full because treated prediction rows vary in x. The
  # control-only regression must refuse without further column dropping.
  X <- cbind(intercept = 1, x = c(rep(0, 10), x[11:20]))
  b <- h$or_public_basis(X, rep(TRUE, nrow(X)))
  testthat::expect_equal(ncol(b$Z), 2L)
  s <- h$or_wls_summary(b$Z[1:10, ], seq_len(10), rep(1, 10))
  testthat::expect_error(h$or_wls_fit(s, h$or_public_rank_tol()), class = "cdid_or_wls_rank_error")
})

testthat::test_that("preparation preserves extraction and selection designs, environments and off-roster covariates", {
  h <- or_preprocess_helpers
  # Fixed finite input: one repeated treated and control unit, a singleton,
  # and an external zero-outcome identity. Preparation does not fit anything.
  data <- data.frame(id = c("c", "c", "t", "t", "s"), time = c(1, 2, 1, 2, 1),
    y = c(0, 1, 1, 3, 2), g = c(0, 0, 2, 2, 0))
  frame <- data.frame(id = c("z", "t", "c", "s"), g = c(0, 2, 0, 0),
    x = c(NA, 2, 1, NA), duplicate = c(NA, 2, 1, NA), q = c(2, 1, 0, -1))
  f <- ~x + duplicate
  p <- h$chained_request_prepare_inputs(data, frame, "y", "time", "id", "g", f,
    outcome_regression = TRUE)
  testthat::expect_identical(p$xformla, f)
  testthat::expect_identical(environment(p$xformla), environment(f))
  testthat::expect_identical(p$frame_ids, frame$id)
  testthat::expect_equal(p$units$n_obs, c(0, 2, 2, 1))
  testthat::expect_identical(p$or_basis$aliased_columns, "duplicate")
  testthat::expect_equal(ncol(p$frame_Z), 2L)
  testthat::expect_identical(p$cp04$or_method, "baseline")
  testthat::expect_true(all(p$cp04$Z == 1))
  c <- h$chained_request_prepare_inputs(data, frame, "y", "time", "id", "g", f,
    selection = "corollary1", selection_xformla = ~q + I(2*q), outcome_regression = TRUE)
  testthat::expect_identical(colnames(c$cp04$Z), c("(Intercept)", "q", "I(2 * q)"))
  testthat::expect_identical(c$or_basis$aliased_columns, "duplicate")
  testthat::expect_equal(ncol(c$cp04$Z), 3L)
  # Baseline needs neither singleton/zero-row outcome covariates nor an imputed
  # cohort for an external frame-only row whose cohort was not supplied.
  unknown <- frame[, setdiff(names(frame), "g")]
  u <- h$chained_request_prepare_inputs(data, unknown, "y", "time", "id", "g", f,
    outcome_regression = TRUE)
  testthat::expect_true(is.na(u$cp04$frame_g[1]))
  testthat::expect_null(u$frame_g_arg)

  # DR preparation freezes the OR outcome recipe while preserving treatment and selection design identities.
  p_dr <- h$chained_request_prepare_inputs(data, frame, "y", "time", "id", "g", f,
    doubly_robust = TRUE)
  testthat::expect_identical(p_dr$cp04$dr_method, "baseline")
  testthat::expect_null(p_dr$cp04$or_method)
  testthat::expect_identical(p_dr$or_basis$aliased_columns, "duplicate")
  testthat::expect_identical(p_dr$cp04$outcome_recipe$columns, c("(Intercept)", "x"))
  testthat::expect_identical(p_dr$cp04$treatment_columns, c("(Intercept)", "x", "duplicate"))
  testthat::expect_equal(ncol(p_dr$frame_Z), 3L)
  c_dr <- h$chained_request_prepare_inputs(data, frame, "y", "time", "id", "g", f,
    selection = "corollary1", selection_xformla = ~q + I(2*q), doubly_robust = TRUE)
  testthat::expect_identical(c_dr$cp04$dr_method, "corollary1")
  testthat::expect_identical(colnames(c_dr$cp04$Z), c("(Intercept)", "q", "I(2 * q)"))
  testthat::expect_identical(c_dr$cp04$outcome_recipe$columns, c("(Intercept)", "x"))
  testthat::expect_identical(c_dr$cp04$treatment_columns, c("(Intercept)", "x", "duplicate"))
  u_dr <- h$chained_request_prepare_inputs(data, unknown, "y", "time", "id", "g", f,
    doubly_robust = TRUE)
  testthat::expect_true(is.na(u_dr$cp04$frame_g[1]))
  testthat::expect_null(u_dr$frame_g_arg)
})

testthat::test_that("built-in and custom probability routes pass the lazy scope gate while unsupported combinations refuse", {
  h <- or_preprocess_helpers
  for (method in c("dr", "ipw", "reg")) {
    testthat::expect_error(h$att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
      data = stop("prepared data reached"), est_method = method), "prepared data reached")
  }
  testthat::expect_error(h$att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    data = stop("prepared data reached")), "prepared data reached")
  e <- testthat::expect_error(h$att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    data = stop("data forced"), est_method = "stdipw"), class = "cdid_chained_unavailable")
  testthat::expect_no_match(conditionMessage(e), "data forced")
  e <- testthat::expect_error(h$att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    data = stop("data forced"), est_method = function(...) stop("callback forced"),
    selection = "corollary1", selection_est_method = "imputation"),
    class = "cdid_chained_request_unsupported_selection")
  testthat::expect_no_match(conditionMessage(e), "data forced|callback forced")
  testthat::expect_error(h$att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    data = stop("prepared data reached"), est_method = function(...) stop("callback forced"),
    selection = "corollary1"), "prepared data reached")
  invalid <- structure(list(), class = "cdid_ipw_pipeline_control")
  for (method in c("dr", "reg")) {
    testthat::expect_error(h$chained_att_gt_request(yname = "y", tname = "time", idname = "id", gname = "g",
      data = stop("development data reached"), est_method = method, control = invalid,
      root_witness = function(...) NULL), "development data reached")
  }
})

testthat::test_that("fixed residual correction preserves long exact dummy, nested and date aliases", {
  h <- or_preprocess_helpers
  profile <- h$chained_public_profile_defaults()
  for (n in c(125000L, 500000L)) {
    i <- seq_len(n); x <- 50 + 10*sin(i)
    state <- (i - 1L) %% 6L; category <- (i - 1L) %% 3L
    designs <- list(
      nested = cbind(intercept = 1, x = x, s1 = as.double(state == 1L),
        s2 = as.double(state == 2L), s3 = as.double(state == 3L),
        s4 = as.double(state == 4L), s5 = as.double(state == 5L), region2 = as.double(state >= 3L)),
      dummy = cbind(intercept = 1, x = x, d0 = as.double(category == 0L),
        d1 = as.double(category == 1L), d2 = as.double(category == 2L)),
      date = cbind(intercept = 1, year = 1990 + (i - 1L) %% 21L,
        centred = 1990 + (i - 1L) %% 21L - 2000))
    expected_alias <- c(nested = "region2", dummy = "d2", date = "centred")
    for (name in names(designs)) {
      X <- designs[[name]]; p <- ncol(X)
      # Verify the exact algebra independently of either QR implementation.
      expected <- switch(name, nested = X[, "s3"] + X[, "s4"] + X[, "s5"],
        dummy = X[, "intercept"] - X[, "d0"] - X[, "d1"], date = X[, "year"] - 2000)
      testthat::expect_identical(unname(X[, expected_alias[name]]), unname(expected))
      b <- h$or_public_basis(X, rep(TRUE, n), profile)
      testthat::expect_identical(b$record$aliased_columns, unname(expected_alias[name]))
      testthat::expect_length(b$record$nontrivial_candidates_retained, 0L)
      testthat::expect_equal(b$record$reconstruction_roundoff_threshold,
        64 * .Machine$double.eps * p, tolerance = 0)
      testthat::expect_lte(max(b$record$reconstruction_max_scaled),
        b$record$reconstruction_roundoff_threshold)
      testthat::expect_equal(b$record$preprocessing_work_ops_bound, 128 * n * p^2, tolerance = 0)
      y <- sin(i/3) + x/10
      ordinary <- stats::lm.fit(X, y, tol = h$or_public_rank_tol())
      reduced <- stats::lm.fit(b$Z, y, tol = h$or_public_rank_tol())
      testthat::expect_identical(names(which(is.na(ordinary$coefficients))), unname(expected_alias[name]))
      testthat::expect_equal(ordinary$residuals, reduced$residuals, tolerance = 1e-10)
      # The control subset covers every factor level. Compare the actual
      # private compact WLS with lm.wfit on the unreduced original columns.
      control <- which(i %% 5L != 0L); w <- 1 + (i[control] %% 3L)/4
      direct <- stats::lm.wfit(X[control, , drop = FALSE], y[control], w,
        tol = h$or_public_rank_tol())
      fit <- h$or_wls_fit(h$or_wls_summary(b$Z[control, , drop = FALSE], y[control], w),
        h$or_public_rank_tol())
      beta <- direct$coefficients; beta[is.na(beta)] <- 0
      predicted <- as.vector(b$Z %*% fit$coefficients)
      reference <- as.vector(X %*% beta)
      testthat::expect_equal(predicted, reference, tolerance = 1e-10)
      testthat::expect_equal(y[control] - predicted[control], direct$residuals, tolerance = 1e-10)
      testthat::expect_equal(fit$rank, direct$rank, tolerance = 0)
    }
  }
})

testthat::test_that("long duplicate aliases remain distinct from near-dependence", {
  h <- or_preprocess_helpers; n <- 125000L; i <- seq_len(n); x <- 50 + 10*sin(i)
  for (delta in c(1e-9, 1e-10)) {
    X <- cbind(intercept = 1, x = x, near = x + delta*((i %% 3L) - 1L), duplicate = x)
    b <- h$or_public_basis(X, rep(TRUE, n), h$chained_public_profile_defaults())
    testthat::expect_identical(b$record$aliased_columns, "duplicate")
    testthat::expect_identical(b$record$retained_columns, c("intercept", "x", "near"))
    testthat::expect_identical(b$record$nontrivial_candidates_retained, "near")
    testthat::expect_gt(b$record$reconstruction_max_scaled["near"],
      b$record$reconstruction_roundoff_threshold)
    s <- h$or_wls_summary(b$Z, sin(i/3) + x/10, 1 + (i %% 3L)/4)
    testthat::expect_error(h$or_wls_fit(s, h$or_public_rank_tol()), class = "cdid_or_wls_rank_error")
  }
})
