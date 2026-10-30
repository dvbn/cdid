# Standalone source-only runs supply this environment with the actual helper
# and conditions from R/; package checks use the loaded candidate namespace.
if (!exists("ipw_start_helpers", inherits = TRUE)) {
  ipw_start_helpers <- asNamespace("cdid")
}

ipw_start_control <- stats::glm.control(epsilon = 1e-12, maxit = 50L, trace = FALSE)

# Independent saturated binary likelihood: columns are response 0/1, rows x=0/1.
ipw_start_root <- function(cells, N) {
  total <- rowSums(cells)
  p <- cells[, 2L] / total
  eta <- stats::qlogis(p)
  v <- total * p * (1 - p)
  list(beta = c(eta[1L], eta[2L] - eta[1L]), p = p,
       information = matrix(c(sum(v), v[2L], v[2L], v[2L]), 2L) / N,
       loglik = sum(cells[, 2L] * log(p) + cells[, 1L] * log1p(-p)))
}

ipw_start_metrics <- function(X, y, w, beta, N) {
  eta <- drop(X %*% beta)
  p <- stats::plogis(eta)
  list(p = p, score = drop(crossprod(X, w * (y - p))) / N,
       information = crossprod(X, w * p * (1 - p) * X) / N,
       loglik = sum(w * (y * eta - pmax(eta, 0) - log1p(exp(-abs(eta))))))
}

ipw_start_fixture <- function() {
  # Actual D39 DGP seed 20268107, pair (1,4), cohorts (3,Inf), N=3200.
  # Provenance: .runtime/cp04-selection/logit_diagnosis/primitive_reference.{R,rds,txt}.
  # Raw rows x=0/1, columns S=0/1; A counts are 617/838 and 285/792.
  treated <- matrix(c(399, 382, 27, 30), 2L)
  control <- matrix(c(394, 386, 9, 3), 2L)
  a <- c(617 / sum(treated), 285 / sum(control))
  b_t <- a[1L] / (treated[, 2L] / rowSums(treated))
  b_c <- a[2L] / (control[, 2L] / rowSums(control))
  counts <- c(treated[, 2L], control[, 2L])
  x <- rep(c(0, 1, 0, 1), counts)
  list(X = cbind(intercept = 1, x = x), y = rep(c(1, 1, 0, 0), counts),
       w = rep(c(b_t, b_c), counts), N = 3200L,
       cells = cbind(control[, 2L] * b_c, treated[, 2L] * b_t))
}

ipw_start_worker <- function(X, y, w, ctrl = ipw_start_control) {
  gi <- ifelse(y == 1, 1L, 2L)
  sup <- list(n = nrow(X), m_cohorts = 2L,
              group_mass = c(sum(w[gi == 1L]), sum(w[gi == 2L])))
  ipw_start_helpers$ipw_nuisance_observer_fit_worker(
    sup, seq_along(y), gi, w, X, ctrl
  )
}

test_that("zero start reaches the independent finite D39 root; default initialization does not", {
  f <- ipw_start_fixture()
  root <- ipw_start_root(f$cells, f$N)
  expect_identical(nrow(f$X), 69L)
  expect_true(all(f$cells > 0))
  expect_gt(min(eigen(root$information, symmetric = TRUE)$values), 0.011)

  # Negative control deliberately reproduces the original single-call defaults.
  old_warnings <- character()
  old <- withCallingHandlers(
    stats::glm.fit(f$X, f$y, weights = f$w, family = stats::binomial("logit"),
                   control = ipw_start_control, intercept = FALSE),
    warning = function(w) {
      old_warnings <<- c(old_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  bad <- ipw_start_metrics(f$X, f$y, f$w, old$coefficients, f$N)
  expect_true(old$converged)
  expect_true(any(grepl("fitted probabilities numerically 0 or 1 occurred", old_warnings, fixed = TRUE)))
  expect_gt(max(abs(bad$score)), 0.04)
  expect_lt(bad$loglik, -1e15)

  gf <- ipw_start_helpers$ipw_glm_fit_internal(f$X, f$y, f$w, ipw_start_control)
  fit <- gf$fit
  actual <- ipw_start_metrics(f$X, f$y, f$w, fit$coefficients, f$N)
  expect_true(fit$converged)
  expect_false(fit$boundary)
  expect_identical(gf$warnings, character())
  expect_lt(max(abs(unname(fit$coefficients) - root$beta)), 1e-12)
  expect_lt(max(abs(actual$p - root$p[f$X[, 2L] + 1L])), 1e-12)
  expect_lt(max(abs(actual$score)), 1e-12)
  expect_lt(max(abs(actual$information - root$information)), 1e-12)
  expect_lt(abs(actual$loglik - root$loglik), 1e-9)
})

test_that("ordinary weighted and unweighted roots retain their closed-form likelihoods", {
  counts <- c(7L, 3L, 2L, 8L) # x=0 failures/successes, x=1 failures/successes.
  x <- rep(c(0, 0, 1, 1), counts)
  y <- rep(c(0, 1, 0, 1), counts)
  for (w_cell in list(rep(1, 4L), c(0.7, 1.3, 2.1, 0.8))) {
    w <- rep(w_cell, counts)
    cells <- matrix(counts * w_cell, 2L, byrow = TRUE)
    root <- ipw_start_root(cells, length(y))
    for (cell_design in c(FALSE, TRUE)) {
      X <- if (cell_design) cbind(1 - x, x) else cbind(1, x)
      beta <- if (cell_design) stats::qlogis(root$p) else root$beta
      fit <- ipw_start_helpers$ipw_glm_fit_internal(X, y, w, ipw_start_control)
      actual <- ipw_start_metrics(X, y, w, fit$fit$coefficients, length(y))
      expect_true(fit$fit$converged)
      expect_identical(fit$warnings, character())
      expect_lt(max(abs(unname(fit$fit$coefficients) - beta)), 1e-12)
      expect_lt(max(abs(actual$score)), 1e-12)
      expect_lt(abs(actual$loglik - root$loglik), 1e-12)
    }
    X <- matrix(1, length(y), 1L)
    fit <- ipw_start_helpers$ipw_glm_fit_internal(X, y, w, ipw_start_control)
    expect_lt(abs(fit$fit$coefficients - stats::qlogis(sum(w * y) / sum(w))), 1e-12)
  }
})

test_that("genuine separation, rank deficiency and caller convergence limits still refuse", {
  f <- ipw_start_fixture()
  expect_error(
    ipw_start_worker(cbind(f$X, duplicate = 1 + 2 * f$X[, 2L]), f$y, f$w),
    class = "cdid_ipw_rank_error"
  )
  expect_error(
    ipw_start_worker(f$X, f$y, f$w, stats::glm.control(epsilon = 1e-12, maxit = 1L)),
    class = "cdid_ipw_convergence_error"
  )
  # Quasi separation: mixed outcomes at x=0, pure outcomes at x=-1/+1.
  expect_error(
    ipw_start_worker(cbind(1, c(-1, 0, 0, 1)), c(0, 0, 1, 1), rep(1, 4L)),
    class = "cdid_ipw_separation_error"
  )
  expect_error(
    ipw_start_worker(cbind(1, c(-15, -15, 15, 15)), c(0, 0, 1, 1), rep(1, 4L)),
    class = "cdid_ipw_boundary_error"
  )
})
