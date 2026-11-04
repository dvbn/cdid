# Unchanged detailed replay and independent raw-row nuisance refits are the
# references. No expected value is produced by the projected helper itself.
dsp <- function(name) getFromNamespace(name, "cdid")
dsp_case <- function(intercept = TRUE, slopes = TRUE, alias = FALSE) {
  j <- rep(1:12, 4); group <- rep(1:4, each = 12); N <- length(j) + 2L
  id <- seq_len(N); g <- c(c(3, 4, 5, Inf)[group], 3, Inf)
  x <- c((j - 6.5) / 5, -.7, .9); z <- c(sin(j) / 2 + j / 30, .2, -.3)
  w <- c((1 + (j %% 4) / 10) * (1 + group * x[1:48] / 30), 1.2, 1.7)
  Y <- vapply(1:5, function(t) c(j / 7 + group * t / 5 + (.2 * t + .03 * t^2) * x[1:48] +
    sin(j + group * t) / 9 + as.double(t >= g[1:48]) * (t + x[1:48] * group / 6), NA_real_, NA_real_), numeric(N))
  Y[c(which(j == 1L), which(j == 4L)), 2L] <- NA_real_
  Y[which(j == 3L), 4L] <- NA_real_
  Y[N, 1L] <- 4.1 # one single-outcome unit; row N-1 has no outcomes
  X <- matrix(numeric(), N, 0L)
  if (intercept) X <- cbind(`(Intercept)` = rep(1, N))
  if (slopes) X <- cbind(X, x = x)
  if (alias) X <- cbind(X, alias = 3 + 2 * x)
  V <- cbind(`(Intercept)` = 1, z = z)
  Xa <- cbind(`(Intercept)` = 1, x = x, z = z)
  if (alias) Xa <- cbind(Xa, alias = X[, "alias"])
  control <- dsp("ipw_owner_control")(7L, 1e8, 2e8, 1e13, 1e8)
  data <- data.frame(id = rep(id, each = 5L), t = rep(1:5, N), g = rep(g, each = 5L), y = as.vector(t(Y)))
  # Exercise both absent rows and stored NA outcomes without changing A.
  data <- data[!(data$id == 2L & data$t == 3L), ]
  Y[2L, 3L] <- NA_real_
  panel <- dsp("construct_ipw_panel_owner")(data, "id", "t", "y", "g", control)
  frame <- dsp("construct_ipw_frame_owner")(panel, id, Xa, w, control)
  reader <- dsp("construct_cp04_frame_reader")(frame, id,
    matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, control)
  defs <- list(list(s = 1, t = 2, g = 3, K = c(4, 5, Inf)),
    list(s = 1, t = 3, g = 3, K = c(4, 5, Inf)), list(s = 2, t = 3, g = 3, K = c(4, 5, Inf)),
    list(s = 1, t = 3, g = 4, K = c(5, Inf)), list(s = 1, t = 3, g = 5, K = Inf),
    list(s = 2, t = 5, g = 3, K = c(4, Inf)), list(s = 1, t = 3, g = 3, K = Inf),
    list(s = 1, t = 3, g = 4, K = Inf), list(s = 1, t = 3, g = 3, K = c(4, 5, Inf)))
  views <- lapply(seq_along(defs), function(i) {
    d <- defs[[i]]
    dsp("construct_ipw_link_view")(frame, d$s, d$t, d$g, d$K, control = control,
      link_id = paste0("projected-", i), functional_key = paste0("cmpcat:projected-", i))
  })
  recipe <- dsp("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  fit <- dsp("fit_dr_baseline")(reader, views, control,
    list(epsilon = 1e-13, maxit = 60L, trace = FALSE), 1e-10, recipe,
    treatment_columns = colnames(V), outcome_pooling = "separate")
  list(fit = fit, reader = reader, Y = Y, X = X, V = V, w = w, g = g, defs = defs,
    intercept = intercept, N = N)
}
dsp_eval <- function(a, rows, links, details = FALSE, scale = "frame")
  dsp("evaluate_dr_separate_if_chunk")(a$fit, rows, links, scale, details)

# Directly refit each requested raw control WLS and pair-risk logit under
# perturbed original-unit weights, recomputing the empirical A-to-pair rates.
# Only base stats fits and the fixture's own data/definitions enter this truth.
dsp_raw <- function(a, frequency = rep(1, a$N)) {
  w <- a$w * frequency; A <- rowSums(!is.na(a$Y)) >= 2L
  vapply(a$defs, function(d) {
    obs <- !is.na(a$Y[, d$s]) & !is.na(a$Y[, d$t]); cohorts <- c(d$g, d$K)
    rates <- vapply(cohorts, function(k) sum(w[a$g == k & obs]) / sum(w[a$g == k & A]), numeric(1))
    rows <- which(obs & a$g %in% cohorts)
    q <- w[rows] / rates[match(a$g[rows], cohorts)]
    dy <- a$Y[rows, d$t] - a$Y[rows, d$s]
    D <- as.double(a$g[rows] == d$g)
    design <- a$X[rows, , drop = FALSE]
    if (a$intercept) design[, -1L] <- design[, -1L, drop = FALSE] * (d$t - d$s) else
      design <- design * (d$t - d$s)
    or_fit <- stats::lm.wfit(design[D == 0, , drop = FALSE], dy[D == 0], q[D == 0], tol = 1e-10)
    stopifnot(all(is.finite(or_fit$coefficients)))
    residual <- dy - as.vector(design %*% or_fit$coefficients)
    trt <- suppressWarnings(stats::glm.fit(a$V[rows, , drop = FALSE], D, weights = q,
      start = rep(0, ncol(a$V)), family = stats::binomial(),
      control = stats::glm.control(epsilon = 1e-13, maxit = 60L), intercept = FALSE))
    stopifnot(trt$converged, all(is.finite(trt$coefficients)))
    prob <- stats::plogis(as.vector(a$V[rows, , drop = FALSE] %*% trt$coefficients))
    stopifnot(all(prob > 1e-6 & prob < 1 - 1e-6))
    weighted.mean(residual, q * D) - weighted.mean(residual, q * (1 - D) * prob / (1 - prob))
  }, numeric(1))
}

testthat::test_that("ordinary incomplete separate DR preserves full coefficient effects and frozen operators", {
  for (options in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE))) {
    a <- do.call(dsp_case, options); fit <- a$fit; L <- length(fit$links)
    testthat::expect_null(fit$separate_complete)
    frozen <- list(groups = fit$outcome_groups, regressions = fit$regressions, treatments = fit$treatments)
    full <- matrix(0, a$N, L); detailed <- full
    for (first in seq.int(1L, a$N, by = 7L)) {
      rows <- seq.int(first, min(a$N, first + 6L))
      ordinary <- dsp_eval(a, rows, seq_len(L))
      old <- dsp_eval(a, rows, seq_len(L), TRUE)
      full[rows, ] <- ordinary$phi; detailed[rows, ] <- old$phi
      testthat::expect_equal(ordinary$phi, old$phi, tolerance = 1e-10)
      testthat::expect_identical(ordinary$id, old$id)
      testthat::expect_identical(ordinary$normalization_count, old$normalization_count)
      # Detailed replay intentionally has larger output/scratch/work bounds.
      # Compare the semantic counts exactly and check each output's own bound.
      bound_fields <- c("output_entries_bound", "concurrent_scratch_bound", "work_ops_bound")
      semantic_fields <- setdiff(names(old$diagnostics), bound_fields)
      testthat::expect_identical(names(ordinary$diagnostics), names(old$diagnostics))
      testthat::expect_identical(ordinary$diagnostics[semantic_fields], old$diagnostics[semantic_fields])
      testthat::expect_lte(unname(dsp("cp04_payload")(ordinary)["entries"]), ordinary$diagnostics$output_entries_bound)
      testthat::expect_lte(unname(dsp("cp04_payload")(old)["entries"]), old$diagnostics$output_entries_bound)
    }
    for (rows in list(7:1, c(50L, 49L, 2L, 1L), 23:19)) for (links in list(c(9L, 3L, 1L), c(8L, 5L, 7L))) {
      ordinary <- dsp_eval(a, rows, links)
      old <- dsp_eval(a, rows, links, TRUE)
      testthat::expect_equal(ordinary$phi, old$phi, tolerance = 1e-10)
      testthat::expect_equal(unname(ordinary$phi), full[rows, links, drop = FALSE], tolerance = 1e-10)
      core <- dsp_eval(a, rows, links, scale = "core")
      testthat::expect_equal(core$phi, a$reader$state$lambda * ordinary$phi, tolerance = 1e-10)
    }
    testthat::expect_equal(crossprod(full) / a$N, crossprod(detailed) / a$N, tolerance = 1e-10)
    testthat::expect_identical(list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments), frozen)
  }
})

testthat::test_that("projected incomplete DR matches independent raw case-weight score derivatives", {
  a <- dsp_case(); fit <- a$fit; L <- length(fit$links)
  testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), dsp_raw(a), tolerance = 1e-9)
  actual <- matrix(0, a$N, L)
  for (first in seq.int(1L, a$N, by = 7L)) {
    rows <- seq.int(first, min(a$N, first + 6L))
    actual[rows, ] <- dsp_eval(a, rows, seq_len(L))$phi
  }
  numerical <- t(vapply(seq_len(a$N), function(i) {
    direction <- rep(-1, a$N); direction[i] <- a$N - 1
    (dsp_raw(a, 1 + 1e-6 * direction) - dsp_raw(a, 1 - 1e-6 * direction)) / 2e-6
  }, numeric(L)))
  testthat::expect_equal(actual, numerical, tolerance = 2e-6)
  testthat::expect_equal(crossprod(actual) / a$N, crossprod(numerical) / a$N, tolerance = 2e-6)
  # Explicitly retain zero- and single-outcome rows in the reference universe.
  testthat::expect_identical(nrow(actual), a$N)
  testthat::expect_equal(actual[49:50, ], numerical[49:50, ], tolerance = 2e-6)
})

testthat::test_that("ordinary baseline projected replay avoids the detailed OR worker", {
  a <- dsp_case(); original <- dsp("or_baseline_separate_if_core")
  n <- 0L
  testthat::local_mocked_bindings(or_baseline_separate_if_core = function(...) {
    n <<- n + 1L; original(...)
  }, .package = "cdid")
  dsp_eval(a, 7:1, c(9L, 3L, 1L))
  testthat::expect_identical(n, 0L)
  dsp_eval(a, 7:1, c(9L, 3L, 1L), TRUE)
  testthat::expect_gt(n, 0L)
})
