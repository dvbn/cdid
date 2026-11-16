test_that("cohort-time nuisance coefficients and SEs match independently perturbed equations", {
  n <- 36L
  i <- seq_len(n)
  cohort <- rep(c(0L, 2L), each = n / 2L)
  cluster <- (i - 1L) %/% 4L
  x <- cos(i / 5) + cohort / 8
  z <- sin(i) + cos(i / 2)
  pair <- i %% 4L != 1L
  dy <- 0.7 * x + 0.5 * sin(i / 3) + cohort * (0.25 + 0.3 * cos(i / 2))
  Y <- cbind(cos(i), cos(i) + ifelse(pair, dy, NA_real_))
  W <- cbind("(Intercept)" = 1, x = x)
  Z <- cbind("(Intercept)" = 1, z = z)
  covariates <- W[, "x", drop = FALSE]
  at <- which(!is.na(Y), arr.ind = TRUE)
  budget <- 1e6
  report <- list(calendar = c(11L, 12L), clusters = cluster, cluster_var = "cluster",
    compute_inffunc = TRUE, chunk_rows = 5L, workspace_entries = budget)
  catalog <- etwfe_corrected_catalog(list(unit = at[, 1L], period = at[, 2L], y = Y[at]), cohort, 2L, budget)
  control <- list(epsilon = 1e-12, maxit = 100L, trace = FALSE)
  obs <- etwfe_observation_fit(Y, cohort, catalog$cells, W,
    list(Z = Z, U = rep(TRUE, n), rho = rep(1, n), universe = "F", model = "pair_logit",
      pooling = "separate"), "dr", "separate", control, 1e-7, 5L, budget, report = report)
  moments <- etwfe_corrected_treatment(catalog, cohort, etwfe_corrected_features(covariates),
    covariates, obs, "nevertreated", "dr", TRUE, control, 5L, budget, report = report)

  candidates <- rbind(obs$nuisance$selection$coefficients, obs$nuisance$imputation$coefficients,
    moments$nuisance$treatment$coefficients,
    subset(moments$nuisance$outcome$coefficients, role != "fixed_gauge"))
  expect_equal(obs$nuisance$selection$cells$s, rep(11, 2L))
  expect_equal(moments$nuisance$treatment$models$g, 12)
  expect_identical(moments$nuisance$outcome$control_blocks$fixed_gauge, c(TRUE, FALSE))
  expect_true(all(is.na(subset(moments$nuisance$outcome$coefficients, role == "fixed_gauge")$std_error)))

  glm_control <- do.call(stats::glm.control, control)
  logit <- function(X, y, weights) {
    unname(stats::glm.fit(X, y, weights = weights, family = stats::quasibinomial(),
      control = glm_control)$coefficients)
  }
  donor <- cohort == 0L
  treated <- !donor
  W_fit <- cbind(1, (x - mean(x)) / max(abs(x - mean(x))))
  refit <- function(weights) {
    select <- vector("list", 2L)
    impute <- vector("list", 2L)
    reconstruct <- numeric(n)
    for (j in 1:2) {
      in_cohort <- if (j == 1L) donor else treated
      complete <- in_cohort & pair
      select[[j]] <- logit(Z[in_cohort, ], as.double(pair[in_cohort]), weights[in_cohort])
      impute[[j]] <- unname(stats::lm.wfit(W[complete, ], dy[complete], weights[complete])$coefficients)
      predict <- drop(W %*% impute[[j]])
      prob <- plogis(drop(Z %*% select[[j]]))
      reconstruct[in_cohort] <- (predict + pair * (dy - predict) / prob)[in_cohort]
    }
    calendar <- unname(stats::lm.wfit(W_fit[donor, ], reconstruct[donor], weights[donor])$coefficients)
    c(unlist(select), unlist(impute), logit(W, as.double(treated), weights), calendar)
  }
  expect_equal(candidates$estimate, refit(rep(1, n)), tolerance = 1e-9)
  eps <- 1e-4
  derivative <- vapply(unique(cluster), function(k) {
    move <- eps * (cluster == k)
    (refit(1 + move) - refit(1 - move)) / (2 * eps)
  }, numeric(nrow(candidates)))
  centered <- sweep(derivative, 1L, rowMeans(derivative), "-")
  expect_equal(candidates$std_error, sqrt(rowSums(centered^2)), tolerance = 1e-7)
})
