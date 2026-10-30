# Run against the installed candidate only; no sourced-function fallback.
if (!requireNamespace("cdid", quietly = TRUE) ||
    !file.exists(file.path(getNamespaceInfo(asNamespace("cdid"), "path"), "Meta", "package.rds")))
  stop("Install candidate cdid before running test-logit-exact-alias.R.", call. = FALSE)
alias_fn <- function(name) getFromNamespace(name, "cdid")
alias_copy <- function(name, solver = NULL, predicate = NULL) {
  fn <- alias_fn(name); env <- new.env(parent = environment(fn))
  if (!is.null(solver)) env$ipw_glm_fit_internal <- solver
  if (!is.null(predicate)) env$ipw_logit_exact_alias <- predicate
  environment(fn) <- env
  fn
}
alias_flags <- function(converged = TRUE, boundary = FALSE) {
  function(Z_use, D_use, fit_weights, ctrl) list(fit = list(
    rank = ncol(Z_use), converged = converged, boundary = boundary, iter = 2L,
    coefficients = rep(0, ncol(Z_use)), qr = list(qr = diag(ncol(Z_use)))), warnings = character())
}
alias_no_fit <- function(...) stop("unexpected glm.fit", call. = FALSE)
alias_X <- cbind(1, rep(c(-1, 1), each = 4)); alias_y <- rep(c(0, 1, 0, 1), 2)
alias_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)

testthat::test_that("CP04 rejects exact aliases before fitting and preserves full-rank behavior", {
  logit <- alias_fn("cp04_logit"); X <- alias_X; y <- alias_y; w <- rep(1, 8)
  dup <- cbind(X, X[, 2]); z <- rep(c(-1, -1, 1, 1), 2)
  call <- function(fn, Z, yy = y, ww = w, gc = alias_gc) fn(Z, yy, ww, gc, 8, list())
  for (Z in list(dup, cbind(X, 0), cbind(0, X))) {
    e <- testthat::expect_error(call(logit, Z), class = "cdid_cp04_rank_error")
    testthat::expect_s3_class(e, "cdid_ipw_fit_error")
    testthat::expect_identical(conditionMessage(e), "Logit has deficient model rank.")
    testthat::expect_error(call(alias_copy("cp04_logit", alias_no_fit), Z), class = "cdid_cp04_rank_error")
  }
  for (solver in list(alias_flags(), alias_flags(boundary = TRUE), alias_flags(converged = FALSE)))
    for (yy in list(y, as.double(X[, 2] > 0)))
      testthat::expect_error(call(alias_copy("cp04_logit", solver), dup, yy), class = "cdid_cp04_rank_error")
  testthat::expect_error(call(logit, dup, rep(.5, 8)), class = "cdid_cp04_unsupported_error")
  testthat::expect_error(call(logit, dup, ww = c(Inf, w[-1])), class = "cdid_cp04_numeric_error")
  # Nearby full-rank designs must behave identically with the predicate bypassed.
  designs <- list(ordinary = X, orthogonal = cbind(X, z),
    near_collinear = cbind(X, X[, 2] + 2^-20 * z), scaled = cbind(X, 2^-20 * z))
  for (kind in names(designs)) {
    Z <- designs[[kind]]
    fit <- call(logit, Z)
    testthat::expect_identical(fit, call(alias_copy("cp04_logit", predicate = function(...) FALSE), Z))
    if (kind == "near_collinear") {
      # Small score does not imply uniform prediction accuracy near rank loss.
      prob <- plogis(drop(Z %*% fit$coefficients))
      score <- drop(crossprod(Z, w * (y - prob))) / 8
      testthat::expect_lte(max(abs(score)), 1e-12)
    } else
      testthat::expect_equal(alias_fn("ipw_logit_values")(Z, fit$coefficients)$prob, rep(.5, 8), tolerance = 1e-12)
  }
  testthat::expect_equal(call(logit, X)$H_chol, diag(.5, 2L), tolerance = 1e-14)
  for (yy in list(rep(0, 8), rep(1, 8), as.double(X[, 2] > 0)))
    testthat::expect_error(call(logit, X, yy), class = "cdid_cp04_separation_error")
  testthat::expect_error(call(alias_copy("cp04_logit", alias_flags(boundary = TRUE)), X), class = "cdid_cp04_boundary_error")
  testthat::expect_error(call(alias_copy("cp04_logit", alias_flags(converged = FALSE)), X), class = "cdid_cp04_convergence_error")
  testthat::expect_error(call(logit, X, c(0, 0, 1, 0, 1, 0, 1, 1), gc = modifyList(alias_gc, list(maxit = 1L))),
    class = "cdid_cp04_convergence_error")
})

testthat::test_that("original observer worker checks the actual pair rows before glm.fit", {
  # Same panel/context/link API as test-chained-ipw-link.R; last two A units miss date 2.
  X <- alias_X; g <- c(ifelse(alias_y == 1, 3, Inf), 3, Inf)
  panel <- data.frame(id = rep(1:10, each = 3), time = rep(1:3, 10), y = rep(1:3, 10), g = rep(g, each = 3))
  panel$y[panel$id > 8 & panel$time == 2] <- NA_real_
  idx <- alias_fn("construct_observed_pairs")(panel, "id", "time", "y", "g")
  Z_all <- cbind(rbind(X, c(1, 0), c(1, 0)), c(X[, 2], 1, -1))
  ctx <- alias_fn("construct_ipw_context")(idx, 1:10, Z_all)
  coh <- alias_fn("validate_ipw_cohorts")(3, Inf)
  sup <- alias_fn("ipw_cohort_support")(ctx, coh)
  risk <- alias_fn("ipw_pair_risk_set")(ctx, coh$cohorts, alias_fn("query_pair_mask_internal")(ctx, 1, 2))
  rows <- risk$risk_a_idx; Z_use <- ctx$Z_a[rows, , drop = FALSE]
  testthat::expect_identical(dim(Z_use), c(8L, 3L))
  testthat::expect_false(alias_fn("ipw_logit_exact_alias")(Z_all))
  worker <- alias_fn("ipw_nuisance_observer_fit_worker")
  call <- function(fn, Z, d = NULL) fn(sup, rows, risk$risk_gi, ctx$w_a[rows], Z, alias_gc, d)
  for (Z in list(Z_use, cbind(X, 0), cbind(0, X))) {
    d <- new.env(parent = emptyenv())
    e <- testthat::expect_error(call(worker, Z, d), class = "cdid_ipw_rank_error")
    testthat::expect_s3_class(e, "cdid_ipw_fit_error")
    testthat::expect_identical(conditionMessage(e), "Weighted propensity logit has deficient model rank.")
    testthat::expect_identical(d$error_diagnostics$stage, "pre_glm_rank_check")
    testthat::expect_error(call(alias_copy("ipw_nuisance_observer_fit_worker", alias_no_fit), Z), class = "cdid_ipw_rank_error")
  }
  for (solver in list(alias_flags(boundary = TRUE), alias_flags(converged = FALSE)))
    testthat::expect_error(call(alias_copy("ipw_nuisance_observer_fit_worker", solver), Z_use), class = "cdid_ipw_rank_error")
  # Verify equal inverse-rate weights from the actual fixture, not fitted output.
  testthat::expect_identical(ctx$w_a, rep(1, 10))
  testthat::expect_identical(sup$group_a_count, c(5L, 5L))
  testthat::expect_identical(sup$group_mass, c(5, 5))
  pair_count <- tabulate(risk$risk_gi, nbins = sup$m_cohorts)
  testthat::expect_identical(pair_count, c(4L, 4L))
  observer_weights <- ctx$w_a[rows] / (pair_count / sup$group_mass)[risk$risk_gi]
  testthat::expect_identical(observer_weights, rep(5 / 4, 8))
  designs <- list(ordinary = X,
    near_collinear = cbind(X, X[, 2] + 2^-20 * rep(c(-1, -1, 1, 1), 2)))
  for (kind in names(designs)) {
    Z <- designs[[kind]]
    fit <- call(worker, Z)
    testthat::expect_identical(fit, call(alias_copy("ipw_nuisance_observer_fit_worker", predicate = function(...) FALSE), Z))
    testthat::expect_identical(fit$fit_weights, observer_weights)
    if (kind == "near_collinear") {
      prob <- plogis(drop(Z %*% fit$beta))
      D <- as.double(risk$risk_gi == 1L)
      score <- drop(crossprod(Z, observer_weights * (D - prob))) / sup$n
      testthat::expect_lte(max(abs(score)), 1e-12)
    } else
      testthat::expect_equal(fit$prob_use, rep(.5, 8), tolerance = 1e-12)
  }
  testthat::expect_error(call(alias_copy("ipw_nuisance_observer_fit_worker", alias_flags(boundary = TRUE)), X), class = "cdid_ipw_boundary_error")
  testthat::expect_error(call(alias_copy("ipw_nuisance_observer_fit_worker", alias_flags(converged = FALSE)), X), class = "cdid_ipw_convergence_error")
  testthat::expect_error(call(alias_copy("ipw_nuisance_observer_fit_worker", alias_flags()), cbind(1, 2 * alias_y - 1)), class = "cdid_ipw_separation_error")
  testthat::expect_error(alias_fn("fit_ipw_link")(ctx, 3, Inf, 1, 2, alias_gc), class = "cdid_ipw_rank_error")
})
