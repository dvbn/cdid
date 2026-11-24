# Independent finite-difference oracle for the within working covariance on
# the observation-adjustment routes (observation weighting, imputation and
# augmentation; same-pair and cross-pair). With fixed support, covariates and
# fitted observation and treatment models, the retained link moments are
# affine in the observed outcomes, so a unit step in one observed outcome
# gives that cell's signed level loading exactly. The oracle Gram sums the
# outer products of these loadings over all observed cells within each target
# cohort; without sampling weights the within row scale is the constant n_A.

within_oracle_panel <- function(N = 40L, TT = 4L, seed = 5L) {
  set.seed(seed)
  g <- rep(c(3, 4, 0), length.out = N)
  x <- stats::rnorm(N)
  d <- expand.grid(t = seq_len(TT), id = seq_len(N))
  d$g <- g[d$id]
  d$x <- x[d$id]
  d$y <- stats::rnorm(nrow(d)) + 0.3 * d$x * d$t + ifelse(d$g > 0 & d$t >= d$g, 1, 0)
  d$y[stats::runif(nrow(d)) < 0.18] <- NA_real_
  d
}

# Retained moments, working covariance and link cohorts as passed to the GMM
# solve; the solve itself is not run.
within_oracle_capture <- function(d, args) {
  cap <- NULL
  local_mocked_bindings(solve_chained_gmm_core = function(design, mhat, V, ...) {
    cap <<- list(m = mhat, V = unname(V), g = design$links$g)
    stop(structure(class = c("within_oracle_capture", "error", "condition"),
      list(message = "captured", call = NULL)))
  }, .package = "cdid")
  tryCatch(do.call(cdid::att_gt, c(list(yname = "y", tname = "t", idname = "id", gname = "g",
    xformla = ~x, data = d, bstrap = FALSE, cband = FALSE, gmm_weighting = "within",
    compute_inffunc = FALSE), args)), within_oracle_capture = function(e) NULL)
  cap
}

within_oracle_check <- function(d, args) {
  ref <- within_oracle_capture(d, args)
  expect_false(is.null(ref))
  cells <- which(!is.na(d$y))
  J <- matrix(0.0, length(ref$m), length(cells))
  for (k in seq_along(cells)) {
    d1 <- d
    d1$y[cells[k]] <- d1$y[cells[k]] + 1
    p <- within_oracle_capture(d1, args)
    expect_identical(names(p$m), names(ref$m))
    J[, k] <- p$m - ref$m
  }
  V0 <- tcrossprod(J) * outer(ref$g, ref$g, "==")
  n_scale <- sum(ref$V * V0) / sum(V0 * V0)
  n_a <- sum(tapply(!is.na(d$y), d$id, sum) >= 2L)
  expect_equal(n_scale, n_a, tolerance = 1e-10)
  expect_lt(max(abs(ref$V - n_a * V0)), 1e-10 * max(abs(ref$V)))
  # Cross-cohort blocks stay exactly zero under within weighting.
  expect_true(all(ref$V[!outer(ref$g, ref$g, "==")] == 0))
}

test_that("within covariance matches a finite-difference loading oracle on augmented routes", {
  d <- within_oracle_panel()
  within_oracle_check(d, list(cross_pair = TRUE, est_method = "reg",
    selection = "conditional_mtar", selection_est_method = "dr"))
  within_oracle_check(d, list(cross_pair = TRUE, est_method = "dr",
    selection = "conditional_mtar", selection_est_method = "dr"))
  within_oracle_check(d, list(cross_pair = FALSE, est_method = "reg",
    selection = "conditional_mtar", selection_est_method = "dr"))
})

test_that("within covariance matches a finite-difference loading oracle on imputation and observation weighting", {
  d <- within_oracle_panel()
  within_oracle_check(d, list(cross_pair = TRUE, est_method = "dr", outcome_pooling = "pooled",
    selection = "conditional_mtar", selection_est_method = "imputation"))
  within_oracle_check(d, list(cross_pair = FALSE, est_method = "dr", outcome_pooling = "pooled",
    selection = "conditional_mtar", selection_est_method = "ipw"))
})
