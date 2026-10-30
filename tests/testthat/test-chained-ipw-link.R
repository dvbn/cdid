# Substantive tests for the private CP02 fixed-operator weighted-pair IPW kernel.
# All `glm.control`, Newton, step-size, and comparison thresholds below are
# test-only fixture settings and do not define a production solver policy.

# R CMD check loads the candidate namespace through tests/testthat.R.
# Standalone runs must explicitly load their verified candidate first.
if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}

build_cp02_fixture_panel <- function(n, Y, observed, G, ids = seq_len(n), times = c(1, 2, 3, 4)) {
  rows <- vector("list", n * length(times))
  k <- 0L
  for (i in seq_len(n)) {
    for (t_idx in seq_along(times)) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        id = ids[i],
        time = as.double(times[t_idx]),
        y = if (observed[i, t_idx]) as.double(Y[i, t_idx]) else NA_real_,
        g = as.double(G[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("fit_ipw_link and replay_ipw_if match the immutable two-pair root plus ONE matched-start adaptation and unchanged seven-equation numerical-Jacobian reference", {
  indep_path <- normalizePath(testthat::test_path(
    "reference", "cp02_weighted_pair_if_independent_reference.R"
  ), winslash = "/", mustWork = TRUE)

  ref <- new.env(parent = globalenv())
  ref$commandArgs <- function(...) paste0("--file=", indep_path)
  ref$write.csv <- function(...) invisible(NULL)
  # The unchanged independent source creates root and calls sys.source there.
  # Bind only that nested load: ref$root$fits/phi/d2/wrong come from the
  # immutable derivative reference plus ONE start argument, not Newton output.
  ref$sys.source <- ipw_reference_source
  expect_identical(digest::digest(file = indep_path, algo = "sha256"),
                   "8815f60a292c61c0391e8ee79ccac0df75a3615fb00bc7bbb9f5c9aaa24fe0e7")
  capture.output(sys.source(indep_path, envir = ref))

  panel <- build_cp02_fixture_panel(ref$n, ref$Y, ref$observed, ref$G)
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(ref$n), Z = ref$Z)

  test_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)
  fit12 <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2, control = test_ctrl, link_id = "12")
  fit13 <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 3, control = test_ctrl, link_id = "13")
  fits <- list("12" = fit12, "13" = fit13)

  # Estimates, rates, and coefficients match the matched-start glm comparator
  # and the unchanged independent Newton/Jacobian equations.
  expect_lt(abs(fit12$estimate - ref$root$fits[[1L]]$value), 1e-14)
  expect_lt(abs(fit13$estimate - ref$root$fits[[2L]]$value), 1e-14)
  expect_lt(max(abs(c(fit12$masses$rates, fit12$coefficients, fit12$treated_mean, fit12$control_mean) -
                    ref$results[["12"]]$theta)), 1e-11)
  expect_lt(max(abs(c(fit13$masses$rates, fit13$coefficients, fit13$treated_mean, fit13$control_mean) -
                    ref$results[["13"]]$theta)), 1e-11)

  # Materialized IF via small chunk size (13) matches both independent references.
  mat <- cdid:::materialize_ipw_if(ctx, fits, max_entries = 480, chunk_size = 13, include_components = TRUE)
  expect_identical(dim(mat$phi_matrix), c(240L, 2L))
  expect_lt(max(abs(mat$phi_matrix - ref$root$phi)), 1e-13)
  expect_lt(max(abs(mat$phi_matrix - ref$root$d2)), 1e-6)
  expect_lt(max(abs(mat$phi_matrix - ref$numeric_phi)), 1e-6)
  expect_lt(max(abs(mat$phi_matrix - ref$numeric_finite)), 1e-6)

  # Bridge difference IF (pair 13 minus pair 12) matches independent case-weight refits.
  bridge_if <- mat$phi_matrix[, 2L] - mat$phi_matrix[, 1L]
  bridge_refit <- ref$root$d2[, 2L] - ref$root$d2[, 1L]
  expect_lt(max(abs(bridge_if - bridge_refit)), 1e-6)

  # All four decomposed IF source components (rate, beta_score, treated, control)
  # match the independent numerical-Jacobian forcing components.
  for (j in seq_along(fits)) {
    nm <- names(fits)[j]
    comp_j <- cbind(
      rate = mat$components$rate[, j],
      beta_score = mat$components$beta_score[, j],
      treated = mat$components$treated[, j],
      control = mat$components$control[, j]
    )
    indep_comp <- as.matrix(ref$unit_rows[[nm]][, c("rate", "beta_score", "treated", "control")])
    expect_lt(max(abs(comp_j - indep_comp)), 1e-6)
  }

  # Outside-risk units (cohort 5, 33 units) have exact zeros in every coordinate/component.
  outside <- !ref$risk
  expect_identical(sum(outside), 33L)
  expect_identical(max(abs(mat$phi_matrix[outside, ])), 0)
  expect_identical(max(abs(mat$components$rate[outside, ])), 0)
  expect_identical(max(abs(mat$components$beta_score[outside, ])), 0)
  expect_identical(max(abs(mat$components$treated[outside, ])), 0)
  expect_identical(max(abs(mat$components$control[outside, ])), 0)

  # Risk-set pair non-observers have nonzero rate/propensity nuisance contributions.
  s12 <- ref$unit_rows[["12"]]$S
  s13 <- ref$unit_rows[["13"]]$S
  nonobs12 <- ref$risk & s12 == 0
  nonobs13 <- ref$risk & s13 == 0
  expect_identical(sum(nonobs12), 146L)
  expect_identical(sum(nonobs13), 141L)
  expect_gt(max(abs(mat$phi_matrix[nonobs12, 1L])), 1.8)
  expect_gt(max(abs(mat$phi_matrix[nonobs13, 2L])), 1.5)
  expect_identical(max(abs(mat$components$treated[nonobs12, 1L])), 0)
  expect_identical(max(abs(mat$components$control[nonobs12, 1L])), 0)
  expect_identical(max(abs(mat$components$beta_score[nonobs12, 1L])), 0)

  # Deliberate formula mutations MUST fail against the independent reference:
  # 1. Omitting all score-rate Jacobian columns (A_{beta, r} = 0)
  for (j in 1:2) {
    bad_fit_omit_all <- fits[[j]]
    bad_fit_omit_all$loadings$loading_rate <- bad_fit_omit_all$partials$G_r
    bad_fit_omit_all$loadings$rate_coef_s0 <- -(bad_fit_omit_all$masses$rates / bad_fit_omit_all$masses$q) * bad_fit_omit_all$loadings$loading_rate
    bad_fit_omit_all$loadings$rate_coef_s1 <- ((1 - bad_fit_omit_all$masses$rates) / bad_fit_omit_all$masses$q) * bad_fit_omit_all$loadings$loading_rate
    bad_phi_omit_all <- cdid:::materialize_ipw_if(ctx, bad_fit_omit_all, max_entries = 240, chunk_size = 50)$phi
    expect_lt(max(abs(bad_phi_omit_all - ref$root$fits[[j]]$wrong)), 1e-13)
    expect_gt(max(abs(bad_phi_omit_all - ref$root$d2[, j])), if (j == 1L) 0.10 else 0.89)
  }

  # 2. Omitting only the target rate column A_{beta, r_g} inside the propensity score
  bad_fit_omit_target <- fit12
  A_no_target <- bad_fit_omit_target$partials$A_beta_r
  A_no_target[, 1L] <- 0
  bad_fit_omit_target$loadings$loading_rate <- as.double(
    bad_fit_omit_target$partials$G_r + drop(crossprod(A_no_target, bad_fit_omit_target$loadings$loading_beta))
  )
  bad_fit_omit_target$loadings$rate_coef_s0 <- -(bad_fit_omit_target$masses$rates / bad_fit_omit_target$masses$q) * bad_fit_omit_target$loadings$loading_rate
  bad_fit_omit_target$loadings$rate_coef_s1 <- ((1 - bad_fit_omit_target$masses$rates) / bad_fit_omit_target$masses$q) * bad_fit_omit_target$loadings$loading_rate
  bad_phi_omit_target <- cdid:::materialize_ipw_if(ctx, bad_fit_omit_target, max_entries = 240, chunk_size = 50)$phi
  expect_gt(max(abs(bad_phi_omit_target - ref$root$d2[, 1L])), 0.10)

  # 3. Reversing the sign of the control beta gradient G_beta
  bad_fit_rev_beta <- fit12
  bad_fit_rev_beta$loadings$loading_beta <- -bad_fit_rev_beta$loadings$loading_beta
  bad_fit_rev_beta$loadings$loading_rate <- as.double(
    bad_fit_rev_beta$partials$G_r + drop(crossprod(bad_fit_rev_beta$partials$A_beta_r, bad_fit_rev_beta$loadings$loading_beta))
  )
  bad_fit_rev_beta$loadings$rate_coef_s0 <- -(bad_fit_rev_beta$masses$rates / bad_fit_rev_beta$masses$q) * bad_fit_rev_beta$loadings$loading_rate
  bad_fit_rev_beta$loadings$rate_coef_s1 <- ((1 - bad_fit_rev_beta$masses$rates) / bad_fit_rev_beta$masses$q) * bad_fit_rev_beta$loadings$loading_rate
  bad_phi_rev_beta <- cdid:::materialize_ipw_if(ctx, bad_fit_rev_beta, max_entries = 240, chunk_size = 50)$phi
  expect_gt(max(abs(bad_phi_rev_beta - ref$root$d2[, 1L])), 30)

  # 4. Wrong normalization by observed-pair count instead of n_a_units
  wrong_norm_phi <- mat$phi_matrix[, 1L] * (fit12$counts$n_observed_risk_units / fit12$counts$n_a_units)
  expect_gt(max(abs(wrong_norm_phi - ref$root$d2[, 1L])), 1.0)
})

test_that("singleton and unobserved units are preserved in metadata but excluded from A-unit rate denominators", {
  set.seed(20261014)
  n_a <- 80L
  x_a <- rnorm(n_a)
  Z_a <- cbind(intercept = 1, x = x_a)
  G_a <- sample(c(3, 4, 5, Inf), n_a, replace = TRUE, prob = c(0.35, 0.25, 0.15, 0.25))
  patterns <- list(c(1, 2, 3), c(1, 2), c(1, 3), c(2, 3))
  patt <- sample(seq_along(patterns), n_a, replace = TRUE)
  obs_a <- matrix(FALSE, n_a, 3L)
  for (i in seq_len(n_a)) obs_a[i, patterns[[patt[i]]]] <- TRUE
  Y_a <- matrix(rnorm(n_a * 3L), n_a, 3L)

  panel_a <- build_cp02_fixture_panel(n_a, Y_a, obs_a, G_a, ids = seq_len(n_a), times = c(1, 2, 3))

  # Append 15 singletons (n_obs == 1) and 5 completely unobserved units (n_obs == 0)
  # heavily concentrated in cohort 3 and cohort 4 so wrong inclusion shifts r_j and q_j.
  sing_ids <- n_a + seq_len(15L)
  sing_g <- c(rep(3, 10L), rep(4, 5L))
  sing_rows <- rbind(
    data.frame(id = sing_ids, time = 1, y = rnorm(15L), g = sing_g),
    data.frame(id = sing_ids, time = 2, y = NA_real_, g = sing_g)
  )
  unobs_ids <- n_a + 15L + seq_len(5L)
  unobs_g <- c(rep(3, 3L), rep(Inf, 2L))
  unobs_rows <- data.frame(id = unobs_ids, time = 1, y = NA_real_, g = unobs_g)

  panel_full <- rbind(panel_a, sing_rows, unobs_rows)
  Z_full <- rbind(Z_a, cbind(intercept = 1, x = rnorm(20L)))
  all_ids <- c(seq_len(n_a), sing_ids, unobs_ids)

  idx_a <- cdid:::construct_observed_pairs(panel_a, "id", "time", "y", "g")
  idx_full <- cdid:::construct_observed_pairs(panel_full, "id", "time", "y", "g")

  ctx_a <- cdid:::construct_ipw_context(idx_a, id = seq_len(n_a), Z = Z_a)
  ctx_full <- cdid:::construct_ipw_context(idx_full, id = all_ids, Z = Z_full)

  # Metadata preserves all 15 singletons and 5 unobserved units.
  expect_identical(ctx_full$n_total_units, 100L)
  expect_identical(ctx_full$n_a_units, 80L)
  expect_identical(ctx_full$n_singleton_units, 15L)
  expect_identical(ctx_full$n_unobserved_units, 5L)
  expect_identical(nrow(ctx_full$singletons), 15L)

  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)
  fit_a <- cdid:::fit_ipw_link(ctx_a, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl)
  fit_full <- cdid:::fit_ipw_link(ctx_full, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl)

  # A-unit rates, cohort shares, coefficients, estimate, and IFs are completely unaffected by singletons/unobserved.
  expect_identical(fit_full$masses$rates, fit_a$masses$rates)
  expect_identical(fit_full$masses$q, fit_a$masses$q)
  expect_identical(fit_full$estimate, fit_a$estimate)
  expect_identical(fit_full$coefficients, fit_a$coefficients)

  phi_a <- cdid:::materialize_ipw_if(ctx_a, fit_a, max_entries = 100, chunk_size = 11)$phi
  phi_full <- cdid:::materialize_ipw_if(ctx_full, fit_full, max_entries = 100, chunk_size = 11)$phi
  expect_identical(phi_full, phi_a)

  # Verify that wrongly including singletons in cohort denominators changes r_g and r_4 materially.
  wrong_r_g <- fit_full$counts$by_cohort_pair[["g"]] / (fit_full$counts$by_cohort_a[["g"]] + 10)
  expect_gt(abs(wrong_r_g - fit_full$masses$rates[["g"]]), 0.05)
})

test_that("context, pair queries, fits, and batched IF replay are strictly invariant to panel row and unit-table permutations", {
  set.seed(99102)
  n <- 50L
  ids_chr <- sprintf("unit_%03d", sample(seq_len(n)))
  x <- rnorm(n)
  Z <- cbind(intercept = 1, x = x)
  w <- runif(n, 0.7, 1.8)
  G <- sample(c(10, 20, Inf), n, replace = TRUE, prob = c(0.4, 0.3, 0.3))
  obs <- matrix(TRUE, n, 3L)
  obs[sample(seq_len(n), 12L), 3L] <- FALSE
  Y <- matrix(rnorm(n * 3L), n, 3L)

  panel <- build_cp02_fixture_panel(n, Y, obs, G, ids = ids_chr, times = c(2.5, 5.0, 9.75))
  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)

  pidx_base <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx_base <- cdid:::construct_ipw_context(pidx_base, id = ids_chr, Z = Z, weights = w)
  fit_base <- cdid:::fit_ipw_link(ctx_base, g = 10, K = c(20, Inf), s = 2.5, t = 9.75, control = ctrl)
  mat_base <- cdid:::materialize_ipw_if(ctx_base, fit_base, max_entries = n, chunk_size = 9)

  for (trial in 1:3) {
    row_perm <- sample(seq_len(nrow(panel)))
    unit_perm <- sample(seq_len(n))
    pidx_perm <- cdid:::construct_observed_pairs(panel[row_perm, , drop = FALSE], "id", "time", "y", "g")
    ctx_perm <- cdid:::construct_ipw_context(
      pidx_perm,
      id = ids_chr[unit_perm],
      Z = Z[unit_perm, , drop = FALSE],
      weights = w[unit_perm]
    )
    fit_perm <- cdid:::fit_ipw_link(ctx_perm, g = 10, K = c(20, Inf), s = 2.5, t = 9.75, control = ctrl)
    mat_perm <- cdid:::materialize_ipw_if(ctx_perm, fit_perm, max_entries = n, chunk_size = 9)

    expect_identical(ctx_perm$a_id, ctx_base$a_id)
    expect_identical(ctx_perm$Z_a, ctx_base$Z_a)
    expect_identical(ctx_perm$w_a, ctx_base$w_a)
    expect_identical(fit_perm$estimate, fit_base$estimate)
    expect_identical(fit_perm$coefficients, fit_base$coefficients)
    expect_identical(mat_perm$id, mat_base$id)
    expect_identical(mat_perm$phi, mat_base$phi)
  }
})

test_that("calendar posting queries preserve close float dates/IDs, extreme finite dates, and binary unit bounds", {
  eps <- .Machine$double.eps
  ids <- c(1, 1 + eps, 1 + 2 * eps, 2, 3, 4)
  t1 <- -1e308
  t2 <- 1 + eps
  t3 <- 1 + 2 * eps
  t4 <- 1e308

  panel <- data.frame(
    id = rep(ids, each = 3L),
    time = c(
      t1, t2, t4,
      t1, t3, t4,
      t1, t2, t4,
      t2, t3, t4,
      t1, t2, t4,
      t1, t3, t4
    ),
    y = c(
      1, 2, 5,
      1, 3, 6,
      2, 4, 9,
      0, 1, 2,
      3, 5, 8,
      4, 7, 11
    ),
    g = rep(c(3, 3, 4, 4, Inf, Inf), each = 3L)
  )
  Z <- cbind(intercept = 1, x = c(-1.2, 0.8, -0.5, 1.1, 0.3, -0.7))

  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = rev(ids), Z = Z[rev(seq_along(ids)), , drop = FALSE])

  expect_identical(ctx$n_calendar_dates, 4L)
  expect_identical(ctx$a_id, ids)

  # Distinguish close float dates t2 vs t3 without collision.
  q_t1_t2 <- cdid:::query_observed_pair(ctx, t1, t2)
  q_t1_t3 <- cdid:::query_observed_pair(ctx, t1, t3)
  expect_identical(q_t1_t2$id, c(ids[1], ids[3], ids[5]))
  expect_identical(q_t1_t3$id, c(ids[2], ids[6]))

  # Extreme finite dates (-1e308, 1e308) whose time gap overflows double precision work cleanly.
  q_ext <- cdid:::query_observed_pair(ctx, t1, t4)
  expect_identical(q_ext$id, ids[c(1, 2, 3, 5, 6)])
  expect_identical(q_ext$dy, c(4, 5, 7, 5, 7))

  # Binary bounds on a_index and unit_index match exact filtering.
  for (lo in 1:6) {
    for (hi in lo:6) {
      sub_q <- cdid:::query_observed_pair(ctx, t1, t4, a_bounds = c(lo, hi))
      exp_q <- q_ext[q_ext$a_index >= lo & q_ext$a_index <= hi, , drop = FALSE]
      rownames(sub_q) <- rownames(exp_q) <- NULL
      expect_identical(sub_q, exp_q)
    }
  }

  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)
  fit_ext <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = t1, t = t4, control = ctrl)
  expect_true(is.finite(fit_ext$estimate))
})

test_that("small chunks match full replay, preserve shared-control covariance, and enforce callback/bound contracts", {
  set.seed(7721)
  n <- 90L
  x <- rnorm(n)
  Z <- cbind(intercept = 1, x = x)
  G <- sample(c(3, 4, Inf), n, replace = TRUE, prob = c(0.35, 0.35, 0.30))
  obs <- matrix(TRUE, n, 3L)
  obs[sample(seq_len(n), 20L), 2L] <- FALSE
  obs[sample(seq_len(n), 18L), 3L] <- FALSE
  # Ensure all units have >= 2 observations
  for (i in which(rowSums(obs) < 2L)) obs[i, 1:2] <- TRUE
  Y <- matrix(rnorm(n * 3L), n, 3L)

  panel <- build_cp02_fixture_panel(n, Y, obs, G, times = c(1, 2, 3))
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z)
  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)

  fit12 <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl, link_id = "p12")
  fit13 <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 3, control = ctrl, link_id = "p13")
  fits <- list(p12 = fit12, p13 = fit13)

  full_mat <- cdid:::materialize_ipw_if(ctx, fits, max_entries = n * 2L, chunk_size = n)
  ref_cov <- crossprod(sweep(full_mat$phi_matrix, 2L, colMeans(full_mat$phi_matrix))) / n
  expect_gt(abs(ref_cov[1L, 2L]), 1e-4)

  chunk_replay_num_ok <- function(actual, expected, tol = 1e-12) {
    identical(typeof(actual), "double") &&
      identical(typeof(expected), "double") &&
      identical(class(actual), class(expected)) &&
      identical(is.object(actual), is.object(expected)) &&
      identical(attributes(actual), attributes(expected)) &&
      identical(length(actual), 180L) &&
      identical(length(expected), 180L) &&
      all(is.finite(actual)) &&
      all(is.finite(expected)) &&
      identical(actual == 0, expected == 0) &&
      identical(1 / actual[actual == 0], 1 / expected[expected == 0]) &&
      identical(sign(actual), sign(expected)) &&
      all(abs(actual - expected) <= tol * pmax(1, abs(actual), abs(expected)))
  }
  expect_chunk_replay_num <- function(actual, expected, tol = 1e-12, info = NULL) {
    expect_true(chunk_replay_num_ok(actual, expected, tol = tol), info = info)
  }
  expect_false(chunk_replay_num_ok(
    full_mat$phi_matrix + 2e-12 * pmax(1, abs(full_mat$phi_matrix)),
    full_mat$phi_matrix
  ))
  expect_false(chunk_replay_num_ok(as.vector(full_mat$phi_matrix), full_mat$phi_matrix))
  expect_false(chunk_replay_num_ok(replace(full_mat$phi_matrix, 1L, 0), full_mat$phi_matrix))
  expect_false(chunk_replay_num_ok(-full_mat$phi_matrix, full_mat$phi_matrix))

  for (cs in c(1L, 2L, 7L, 23L, n)) {
    batches <- list()
    info <- cdid:::replay_ipw_if(ctx, fits, chunk_size = cs, callback = function(b) {
      expect_lte(nrow(b$phi), cs)
      batches[[length(batches) + 1L]] <<- b
      "ignored_return_value"
    })
    expect_identical(info$n_a_units, n)
    expect_lte(info$max_batch_units, cs)
    reconstructed <- do.call(rbind, lapply(batches, `[[`, "phi"))
    reconstructed_ids <- unlist(lapply(batches, `[[`, "id"), use.names = FALSE)
    expect_identical(reconstructed_ids, seq_len(n))
    expect_chunk_replay_num(reconstructed, full_mat$phi_matrix)

    cov_cs <- cdid:::stream_ipw_covariance(ctx, fits, chunk_size = cs)
    expect_lt(max(abs(cov_cs$covariance - ref_cov)), 1e-14)
  }

  # No callback is retained in context or fit objects.
  expect_false(any(vapply(ctx, is.function, logical(1L))))
  expect_false(any(vapply(fit12, is.function, logical(1L))))

  # Callback errors propagate unchanged.
  expect_error(
    cdid:::replay_ipw_if(ctx, fits, chunk_size = 5L, callback = function(b) {
      if (b$batch_index == 2L) {
        stop(structure(list(message = "custom_callback_abort", call = NULL),
                       class = c("custom_cb_err", "error", "condition")))
      }
    }),
    class = "custom_cb_err"
  )

  # Caller-explicit materialization bound rejects n * L > max_entries before allocation.
  expect_error(
    cdid:::materialize_ipw_if(ctx, fits, max_entries = n * 2L - 1L, chunk_size = 10L),
    class = "cdid_ipw_bound_error"
  )

  # Mismatched contexts across fits are rejected during joint replay.
  ctx_other <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z)
  fit_other <- cdid:::fit_ipw_link(ctx_other, g = 3, K = c(4, Inf), s = 1, t = 3, control = ctrl)
  expect_error(
    cdid:::replay_ipw_if(ctx, list(fit12, fit_other), chunk_size = 10L, callback = function(b) NULL),
    class = "cdid_ipw_context_mismatch_error"
  )
})

test_that("fit_ipw_link rejects missing control, inferred labels, unsupported groups, and rank/convergence/boundary/numerical/normalization failures", {
  set.seed(31415)
  n <- 60L
  x <- rnorm(n)
  Z <- cbind(intercept = 1, x = x)
  G <- rep(c(3, 4, Inf), each = 20L)
  obs <- matrix(TRUE, n, 3L)
  Y <- matrix(rnorm(n * 3L), n, 3L)
  panel <- build_cp02_fixture_panel(n, Y, obs, G, times = c(1, 2, 3))
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z)
  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)

  # 1. Missing or invalid glm.control
  expect_error(
    cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2),
    class = "cdid_ipw_control_error"
  )
  expect_error(
    cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2, control = list(epsilon = -1, maxit = 10, trace = FALSE)),
    class = "cdid_ipw_control_error"
  )

  # 2. Control label inference ("notyettreated" / "nevertreated") is refused
  expect_error(
    cdid:::fit_ipw_link(ctx, g = 3, K = "notyettreated", s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_cohort_error"
  )
  expect_error(
    cdid:::fit_ipw_link(ctx, g = 3, K = c(3, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_cohort_error"
  )

  # 3. Invalid pair descriptor (s >= t)
  expect_error(
    cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 2, t = 2, control = ctrl),
    class = "cdid_ipw_pair_error"
  )
  expect_error(
    cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 3, t = 1, control = ctrl),
    class = "cdid_ipw_pair_error"
  )

  # 4. Declared cohort absent from A-units (zero A-unit support)
  expect_error(
    cdid:::fit_ipw_link(ctx, g = 3, K = c(4, 99, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_cohort_support_error"
  )

  # 5. Declared cohort with positive A-unit support but zero observed pairs on (1, 3)
  obs_zero_k4 <- obs
  obs_zero_k4[G == 4, 3L] <- FALSE
  panel_zero_k4 <- build_cp02_fixture_panel(n, Y, obs_zero_k4, G, times = c(1, 2, 3))
  ctx_zero_k4 <- cdid:::construct_ipw_context(
    cdid:::construct_observed_pairs(panel_zero_k4, "id", "time", "y", "g"),
    id = seq_len(n), Z = Z
  )
  expect_error(
    cdid:::fit_ipw_link(ctx_zero_k4, g = 3, K = c(4, Inf), s = 1, t = 3, control = ctrl),
    class = "cdid_ipw_pair_support_error"
  )

  # 6. Deficient model rank (collinear columns in Z)
  Z_collinear <- cbind(intercept = 1, x = x, x_dup = 2 * x + 1)
  ctx_rank <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z_collinear)
  expect_error(
    cdid:::fit_ipw_link(ctx_rank, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_rank_error"
  )

  # 7. Nonconvergence under caller-supplied glm.control (maxit = 1, tight epsilon)
  Z_nonconv <- cbind(intercept = 1, x = ifelse(G == 3, 0.8, -0.4) + rnorm(n, sd = 0.5))
  ctx_nonconv <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z_nonconv)
  expect_error(
    cdid:::fit_ipw_link(
      ctx_nonconv, g = 3, K = c(4, Inf), s = 1, t = 2,
      control = stats::glm.control(epsilon = 1e-15, maxit = 1, trace = FALSE)
    ),
    class = "cdid_ipw_convergence_error"
  )

  # 8. Complete separation / boundary fitted probabilities (0 or 1)
  Z_sep <- cbind(intercept = 1, x = ifelse(G == 3, 15, -15))
  ctx_sep <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z_sep)
  expect_error(
    cdid:::fit_ipw_link(ctx_sep, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_boundary_error"
  )

  # 9. Numerical failure (overflow in Hessian / covariate crossprod)
  Z_overflow <- cbind(intercept = 1, x = rep(c(-1e200, 1e200), length.out = n))
  ctx_overflow <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z_overflow)
  expect_error(
    cdid:::fit_ipw_link(ctx_overflow, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_numeric_error"
  )

  # 10. Invalid normalization (overflow in weighted outcome sum)
  Y_norm_fail <- Y
  Y_norm_fail[G == 3, 1L] <- 0
  Y_norm_fail[G == 3, 2L] <- 1e308
  panel_norm_fail <- build_cp02_fixture_panel(n, Y_norm_fail, obs, G, times = c(1, 2, 3))
  ctx_norm_fail <- cdid:::construct_ipw_context(
    cdid:::construct_observed_pairs(panel_norm_fail, "id", "time", "y", "g"),
    id = seq_len(n), Z = Z
  )
  expect_error(
    cdid:::fit_ipw_link(ctx_norm_fail, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_normalization_error"
  )
})

test_that("fitted cache retains no N-row arrays or inverse matrices and shares one locked context", {
  set.seed(8812)
  for (n in c(120L, 600L)) {
    Z <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n))
    G <- sample(c(3, 4, Inf), n, replace = TRUE, prob = c(0.35, 0.35, 0.30))
    obs <- matrix(TRUE, n, 3L)
    Y <- matrix(rnorm(n * 3L), n, 3L)
    panel <- build_cp02_fixture_panel(n, Y, obs, G, times = c(1, 2, 3))
    pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
    ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z)
    expect_true(environmentIsLocked(ctx$ref))

    ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)
    fit1 <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl)
    fit2 <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 3, control = ctrl)
    expect_identical(fit1$context_ref, fit2$context_ref)

    # Every non-environment element inside fit1 has max length <= max(p^2, p * M, M) = 9.
    non_env <- fit1[setdiff(names(fit1), "context_ref")]
    atomic_lengths <- unlist(rapply(non_env, length, how = "unlist"))
    expect_lte(max(atomic_lengths), 9L)
    expect_false("H_inv" %in% names(fit1$factors))
    expect_true(is.matrix(fit1$factors$H_chol))
    expect_null(fit1$pair$functional_key)
  }
})

test_that("non-unit mechanical case weights match independent case-perturbation refits and context validates IDs/weights", {
  set.seed(5123)
  n <- 55L
  x <- rnorm(n)
  Z <- cbind(intercept = 1, x = x)
  G <- sample(c(3, 4, Inf), n, replace = TRUE, prob = c(0.38, 0.32, 0.30))
  w0 <- runif(n, 0.6, 1.9)
  obs <- matrix(TRUE, n, 3L)
  obs[sample(seq_len(n), 14L), 3L] <- FALSE
  Y <- matrix(rnorm(n * 3L), n, 3L)
  panel <- build_cp02_fixture_panel(n, Y, obs, G, times = c(1, 2, 3))
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)

  ctx0 <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z, weights = w0)
  fit0 <- cdid:::fit_ipw_link(ctx0, g = 3, K = c(4, Inf), s = 1, t = 3, control = ctrl, functional_key = "caller_key")
  expect_identical(fit0$pair$functional_key, "caller_key")
  phi0 <- cdid:::materialize_ipw_if(ctx0, fit0, max_entries = n, chunk_size = 12L)$phi

  h <- 5e-5
  refit_if <- vapply(seq_len(n), function(i) {
    wp <- wm <- w0
    wp[i] <- w0[i] * (1 + h)
    wm[i] <- w0[i] * (1 - h)
    fp <- cdid:::fit_ipw_link(cdid:::construct_ipw_context(pidx, seq_len(n), Z, wp), 3, c(4, Inf), 1, 3, ctrl)$estimate
    fm <- cdid:::fit_ipw_link(cdid:::construct_ipw_context(pidx, seq_len(n), Z, wm), 3, c(4, Inf), 1, 3, ctrl)$estimate
    n * (fp - fm) / (2 * h)
  }, numeric(1L))
  expect_lt(max(abs(phi0 - refit_if)), 1e-6)

  # Reject character IDs against numeric pair_index IDs (no string coercion),
  # duplicate IDs, missing eligible A-unit IDs, and nonpositive/nonfinite weights.
  expect_error(
    cdid:::construct_ipw_context(pidx, id = as.character(seq_len(n)), Z = Z),
    class = "cdid_ipw_id_error"
  )
  expect_error(
    cdid:::construct_ipw_context(pidx, id = c(1L, 1L, 3:n), Z = Z),
    class = "cdid_ipw_id_error"
  )
  expect_error(
    cdid:::construct_ipw_context(pidx, id = seq_len(n - 1L), Z = Z[-n, , drop = FALSE]),
    class = "cdid_ipw_id_error"
  )
  expect_error(
    cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z, weights = c(0, w0[-1L])),
    class = "cdid_ipw_weight_error"
  )

  # Cohort present ONLY in singletons has zero A-unit support and is rejected.
  panel_sing_cohort <- rbind(
    panel,
    data.frame(id = c(n + 1L, n + 2L), time = c(1, 1), y = c(1.0, 2.0), g = c(77, 77))
  )
  pidx_sing <- cdid:::construct_observed_pairs(panel_sing_cohort, "id", "time", "y", "g")
  ctx_sing <- cdid:::construct_ipw_context(pidx_sing, id = seq_len(n + 2L), Z = rbind(Z, c(1, 0.1), c(1, -0.2)))
  expect_error(
    cdid:::fit_ipw_link(ctx_sing, g = 3, K = c(4, 77, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_cohort_support_error"
  )
})

test_that("frozen context and explicit bounds prevent inconsistent replay and unbounded allocation", {
  set.seed(632819)
  n <- 60L
  G <- rep(c(3, 4, Inf), each = 20L)
  Z <- cbind(1, rnorm(n))
  Y <- matrix(rnorm(n * 3L), n, 3L)
  panel <- build_cp02_fixture_panel(n, Y, matrix(TRUE, n, 3L), G, times = 1:3)
  idx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(idx, 1:n, Z)
  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100)
  fit <- cdid:::fit_ipw_link(ctx, 3, c(4, Inf), 1, 2, ctrl)
  baseline <- cdid:::materialize_ipw_if(ctx, fit, n, 7)$phi
  changed <- ctx
  changed$observations$y <- changed$observations$y * 2
  expect_identical(changed$ref, ctx$ref)
  expect_error(cdid:::materialize_ipw_if(changed, fit, n, 7), class = "cdid_ipw_context_mismatch_error")
  expect_error(cdid:::fit_ipw_link(changed, 3, c(4, Inf), 1, 2, ctrl), class = "cdid_ipw_context_mismatch_error")
  expect_error(cdid:::query_observed_pair(changed, 1, 2), class = "cdid_ipw_context_mismatch_error")
  expect_error(cdid:::get_calendar_posting(changed, 1), class = "cdid_ipw_context_mismatch_error")
  expect_identical(cdid:::materialize_ipw_if(ctx, fit, n, 7)$phi, baseline)
  expect_error(cdid:::materialize_ipw_if(ctx, fit, NULL, 7), class = "cdid_ipw_bound_error")
  expect_error(cdid:::materialize_ipw_if(ctx, fit, chunk_size = 7), class = "cdid_ipw_bound_error")
  expect_error(cdid:::materialize_ipw_if(ctx, fit, n, 7, NA), class = "cdid_ipw_input_error")
  expect_error(cdid:::materialize_ipw_if(ctx, fit, n, 7, matrix(TRUE)), class = "cdid_ipw_input_error")
  comp <- cdid:::materialize_ipw_if(ctx, fit, n, 7, TRUE)
  expect_identical(comp$preflight$materialized_matrix_count, 5L)
  expect_equal(comp$preflight$materialized_matrix_entries, 5 * n)
  expect_equal(comp$preflight$covariance_matrix_entries, 1)
  for (bad in list(c(ctrl, list(epsilon = 1e-4)), c(ctrl, list(unused = 1)),
                   modifyList(ctrl, list(maxit = 2^31)))) {
    expect_error(cdid:::fit_ipw_link(ctx, 3, c(4, Inf), 1, 2, bad), class = "cdid_ipw_control_error")
  }
  huge_panel <- panel
  huge_panel$y <- huge_panel$y * 1e155
  huge_ctx <- cdid:::construct_ipw_context(
    cdid:::construct_observed_pairs(huge_panel, "id", "time", "y", "g"), 1:n, Z)
  huge_fit <- cdid:::fit_ipw_link(huge_ctx, 3, c(4, Inf), 1, 2, ctrl)
  expect_true(all(is.finite(cdid:::materialize_ipw_if(huge_ctx, huge_fit, n, 7)$phi)))
  expect_error(cdid:::stream_ipw_covariance(huge_ctx, huge_fit, 7), class = "cdid_ipw_numeric_error")
  many <- rep(list(fit), 46341L)
  pre <- cdid:::preflight_ipw_replay(ctx, many, 1)
  expect_gt(pre$covariance_matrix_entries, .Machine$integer.max)
  expect_error(cdid:::stream_ipw_covariance(ctx, many, 1), class = "cdid_ipw_bound_error")
})

