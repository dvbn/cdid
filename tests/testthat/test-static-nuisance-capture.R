make_stat_report_fixture <- function() {
  set.seed(917L)
  n <- 240L
  id <- seq_len(n)
  x <- rnorm(n)
  switch <- rep(c(FALSE, TRUE), each = n / 2L)
  cluster <- (id - 1L) %% 20L + 1L
  early <- rnorm(n)
  change <- 0.3 + 0.5 * x + switch + rnorm(n)
  keep_s <- stats::runif(n) < stats::plogis(1.4 + 0.25 * x)
  keep_u <- stats::runif(n) < stats::plogis(1.1 - 0.3 * x)
  complete <- data.frame(id = rep(id, 2L), time = rep(c(2, 5), each = n), x = rep(x, 2L),
    y = c(early, early + change), d = c(rep(0L, n), as.integer(switch)),
    cluster = rep(cluster, 2L))
  missing <- complete
  missing$y[!c(keep_s, keep_u)] <- NA_real_
  list(complete = complete, missing = missing, x = x, switch = switch, cluster = cluster,
    change = change, paired = keep_s & keep_u, any_observed = keep_s | keep_u, n = n)
}

stat_report_fixture_fit <- function(data, ..., xformla = ~ x, control = list(keep_nuisance = TRUE)) {
  att_gt(yname = "y", tname = "time", idname = "id", dname = "d", data = data,
    xformla = xformla, switchers = "entry", clustervars = "cluster", bstrap = FALSE, cband = FALSE,
    control = control, ...)
}

# Independent cluster sandwich from a stacked design, raw scores, and Hessian.
stat_report_direct_coef <- function(X, y, rows, cluster, logit) {
  Xf <- X[rows, , drop = FALSE]
  yf <- y[rows]
  if (logit) {
    fit <- stats::glm.fit(Xf, yf, family = stats::binomial(), control = stats::glm.control(epsilon = 1e-12))
    fitted <- fit$fitted.values
    hessian <- crossprod(Xf, Xf * (fitted * (1 - fitted)))
  } else {
    fit <- stats::lm.fit(Xf, yf)
    fitted <- fit$fitted.values
    hessian <- crossprod(Xf)
  }
  scores <- matrix(0, nrow(X), ncol(X))
  scores[rows, ] <- (yf - fitted) * Xf
  cluster_if <- rowsum(t(solve(hessian, t(scores))), cluster)
  centered <- sweep(cluster_if, 2L, colMeans(cluster_if))
  cbind(estimate = unname(fit$coefficients), std_error = sqrt(colSums(centered^2)))
}

test_that("static complete-data reports use the original risk-set and stayer fits", {
  data <- make_stat_report_fixture()
  pooling <- c(t5_entry = "stayers, custom")
  fit <- stat_report_fixture_fit(data$complete, outcome_pooling = pooling)
  expect_s3_class(nuisance_models(fit), "cdid_nuisance_index")
  expect_identical(nuisance_models(fit, "selection")$code, "no_selection_logit")
  expect_identical(nuisance_models(fit, "imputation")$code, "no_imputation")
  X <- cbind(1, data$x)
  for (role in c("treatment", "outcome")) {
    block <- nuisance_models(fit, role, block = 1L)
    expect_s3_class(block, "cdid_nuisance_block")
    expect_match(block$equation, "X_ref", fixed = TRUE)
    expect_match(block$sample, "stayer", fixed = TRUE)
    logit <- role == "treatment"
    keep <- if (logit) rep(TRUE, data$n) else !data$switch
    y <- if (logit) as.integer(data$switch) else data$change
    direct <- stat_report_direct_coef(X, y, keep, data$cluster, logit)
    expect_equal(unname(as.matrix(block$table[, c("estimate", "std_error")])), unname(direct), tolerance = 2e-8)
    expect_identical(block$covariance$type, "cluster")
    expect_identical(block$covariance$cluster_var, "cluster")
  }

  expect_identical(nuisance_models(fit, "outcome", block = "stayers, custom")$comparisons$comparison,
    fit$target_id)
  point <- stat_report_fixture_fit(data$complete, outcome_pooling = pooling, compute_inffunc = FALSE)
  expect_identical(fit$att, point$att)
  for (role in c("treatment", "outcome")) {
    block <- nuisance_models(point, role, block = 1L)
    expect_identical(block$covariance$status, "unavailable")
    expect_identical(block$covariance$code, "point_only")
    expect_true(all(is.na(block$table$std_error)))
    expect_identical(block$table$estimate, nuisance_models(fit, role, block = 1L)$table$estimate)
  }
})

test_that("static observation and imputation reports match pooled role coefficients and clustered scores", {
  data <- make_stat_report_fixture()
  fit <- stat_report_fixture_fit(data$missing, selection = "conditional_mtar",
    selection_est_method = "dr", selection_model = "pair_logit", selection_pooling = "pooled",
    imputation_pooling = "pooled", imputation_scale = "per_period")
  expect_true(all(nuisance_models(fit)$components$status == "available"))
  for (role in c("selection", "imputation")) {
    comp <- nuisance_models(fit, role)
    expect_identical(comp$blocks$block, "pooled")
    expect_equal(comp$cells$s, c(2, 2))
    expect_equal(comp$cells$t, c(5, 5))
    expect_identical(comp$cells$role, c("switcher", "stayer"))
    expect_identical(comp$coefficients$role, c("local_intercept", "local_intercept", "shared_slope"))
    logit <- role == "selection"
    # Dates 2 and 5 are consecutive observed panel positions, so the
    # imputation design uses index gap 1 even though the dates differ by 3.
    X <- cbind(data$switch, !data$switch, data$x)
    y <- if (logit) data$paired else data$change
    keep <- if (logit) rep(TRUE, data$n) else data$paired
    direct <- stat_report_direct_coef(X, y, keep, data$cluster, logit)
    expect_equal(unname(as.matrix(comp$coefficients[, c("estimate", "std_error")])),
      unname(direct), tolerance = 2e-8)
    expect_equal(comp$blocks$n_fitting_rows, sum(keep))
    expect_output(print(nuisance_models(fit, role, block = "pooled")), "X_ref", fixed = TRUE)
    if (!logit) {
      expect_identical(comp$imputation_scale, "per_period")
      expect_match(comp$blocks$equation, "u_idx - s_idx", fixed = TRUE)
      expect_equal(comp$comparisons$index_gap, c(1, 1))
    }
  }
  outcome <- nuisance_models(fit, "outcome", block = 1L)
  expect_match(outcome$equation, "Reconstructed mean stayer change", fixed = TRUE)
  expect_match(outcome$covariance$accounts_for, "fitted observation or imputation", fixed = TRUE)
  expect_identical(outcome$comparisons$comparison, fit$target_id)
  captured <- outcome$table$estimate
  actual <- fit$fit$nuisance$or
  expect_equal(captured, unname(c(actual$alpha_hat, actual$b_hat[[1L]])))
  xi <- actual$if_callback(seq_len(data$n))
  cluster_if <- rowsum(cbind(xi$chi_alpha, xi$chi_b[[1L]]), data$cluster)
  expected_se <- sqrt(colSums(sweep(cluster_if, 2L, colMeans(cluster_if))^2)) / data$n
  expect_equal(outcome$table$std_error, unname(expected_se), tolerance = 2e-12)
})

test_that("static date logits report endpoint fits and transported pair logits use the requested universe", {
  data <- make_stat_report_fixture()
  fit_d <- stat_report_fixture_fit(data$missing, est_method = "ipw", selection = "conditional_mtar",
    selection_model = "independent_period_logit", selection_pooling = "within_group")
  dates <- nuisance_models(fit_d, "selection")
  expect_equal(sort(dates$cells$t), c(2, 2, 5, 5))
  expect_match(dates$derived, "product of the endpoint probabilities", fixed = TRUE)
  expect_equal(nrow(dates$pair_cohorts), 2L)

  fit_b <- stat_report_fixture_fit(data$missing, est_method = "ipw", selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "pooled", selection_universe = "B",
    selection_justification = "Synthetic witness with a declared B fitting population.")
  pair <- nuisance_models(fit_b, "selection")
  expect_identical(pair$universe, "B")
  expect_equal(pair$blocks$n_fitting_rows, sum(data$any_observed))
  expect_equal(pair$cells$n_observed, c(sum(data$paired & data$switch), sum(data$paired & !data$switch)))
  expect_true(all(is.na(pair$cells$p_min) & is.na(pair$cells$p_max) & is.na(pair$cells$n_floor)))
  X <- cbind(data$switch, !data$switch, data$x)
  direct <- stat_report_direct_coef(X, data$paired, data$any_observed, data$cluster, TRUE)
  expect_equal(unname(as.matrix(pair$coefficients[, c("estimate", "std_error")])),
    unname(direct), tolerance = 2e-8)
})

test_that("tight static score budgets preserve retained metadata and reserve dense reconstruction solves", {
  data <- make_stat_report_fixture()
  data$missing$id <- paste0("unit_", data$missing$id)
  args <- list(est_method = "ipw", selection = "conditional_mtar",
    selection_model = "pair_logit", selection_pooling = "pooled")
  reference <- do.call(stat_report_fixture_fit, c(list(data = data$missing), args))
  count <- reference$request$control$estimated_work_ops
  limited <- do.call(stat_report_fixture_fit, c(list(data = data$missing,
    control = list(keep_nuisance = TRUE, counted_work_limit = count)), args))
  full <- nuisance_models(reference, "selection")
  tight <- nuisance_models(limited, "selection")
  expect_identical(tight$covariance$code, "budget")
  expect_identical(tight$coefficients$estimate, full$coefficients$estimate)
  expect_identical(tight$cells, full$cells)
  expect_equal(tight$cells$n_fitting_rows, c(sum(data$switch), sum(!data$switch)))
  expect_true(all(is.na(tight$cells$p_min) & is.na(tight$cells$p_max) & is.na(tight$cells$n_floor)))

  ample_batch <- do.call(stat_report_fixture_fit, c(list(data = data$missing,
    control = list(keep_nuisance = TRUE, counted_work_limit = count + 100000)), args))
  one_row <- do.call(stat_report_fixture_fit, c(list(data = data$missing,
    control = list(keep_nuisance = TRUE, counted_work_limit = count + 100000, chunk_rows = 1L)), args))
  expect_identical(nuisance_models(ample_batch, "selection")$covariance$code, "ok")
  expect_identical(nuisance_models(one_row, "selection")$covariance$code, "budget")
  expect_equal(nuisance_models(one_row, "selection")$coefficients$estimate, full$coefficients$estimate)

  z <- paste0("z", seq_len(12L))
  data$missing[z] <- lapply(z, function(key) rep(rnorm(data$n), 2L))
  formula <- stats::reformulate(z)
  reference <- stat_report_fixture_fit(data$missing, est_method = "reg", xformla = formula)
  count <- reference$request$control$estimated_work_ops
  limited <- stat_report_fixture_fit(data$missing, est_method = "reg", xformla = formula,
    control = list(keep_nuisance = TRUE, counted_work_limit = count + 500000))
  full <- nuisance_models(reference, "outcome")
  tight <- nuisance_models(limited, "outcome")
  expect_identical(full$covariance$status, "available")
  expect_identical(tight$covariance$code, "budget")
  expect_equal(tight$coefficients$estimate, full$coefficients$estimate, tolerance = 0)
  expect_identical(limited$att, reference$att)
})
