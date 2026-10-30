# SPDX-License-Identifier: GPL-3.0-only
#
# D90/D92 shared full-width covariance planner, generalized retained raw
# original-unit moment IF carrier across all supported CP04 producers, and
# beta-only OR-core evaluation for non-detailed corrected DR and dated baseline DR.
#
# STATUS: UNEXECUTED in the bounded authoring job (base R parse/AST check only;
# root sequences isolated installation and numerical gates).
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# Root runs the existing independent d60 references from the full checkout.
# This focused file uses the installed namespace; no unused source dependency.

# Deterministic route panel with uneven missingness, time-varying covariates/weights,
# known transport ratio rho, and cluster grouping.
d92_route_panel <- function(n = 36L, dates = c(1, 2, 4, 5), cohorts = c(4, 5, 0),
                            vary_x = FALSE, vary_w = FALSE, drop_y = TRUE) {
  id <- seq_len(n)
  g <- rep(cohorts, length.out = n)
  d <- expand.grid(time = dates, id = id, KEEP.OUT.ATTRS = FALSE)[, c("id", "time")]
  ti <- match(d$time, dates)
  d$G <- g[d$id]
  bx <- round(sin(1.7 * id), 2)
  d$x <- if (vary_x) round(bx[d$id] + 0.25 * ti + 0.3 * cos(2.3 * d$id + 1.1 * ti), 2) else bx[d$id]
  bw <- round(0.6 + (id %% 7L) / 7, 3)
  d$w <- if (vary_w) round(bw[d$id] * (1 + 0.15 * ti * ((d$id %% 2L) == 0L)), 3) else bw[d$id]
  te <- ifelse(d$G > 0 & d$time >= d$G, 1 + 0.5 * (d$time - d$G), 0)
  d$y <- round(cos(d$id) + 0.4 * d$time + 0.8 * d$x + te + 0.5 * sin(5.1 * d$id + 2.7 * ti), 3)
  if (drop_y) d$y[((d$id * 37L + ti * 11L) %% 13L) < 3L] <- NA_real_
  d$rho <- round(0.6 + 0.3 * (d$id %% 4L) / 4, 3)
  d$cl <- (d$id - 1L) %/% 3L + 1L
  d
}

d92_custom_reg <- function(y1, y0, D, covariates, i.weights = NULL, inffunc = TRUE, ...) {
  w <- i.weights / mean(i.weights)
  dy <- y1 - y0
  X <- as.matrix(covariates)
  c0 <- (1 - D) * w
  M <- crossprod(X, X * c0) / length(D)
  beta <- solve(crossprod(X, X * c0), crossprod(X, c0 * dy))
  fit <- drop(X %*% beta)
  p1 <- mean(w * D)
  att <- sum(w * D * (dy - fit)) / sum(w * D)
  mx1 <- colMeans(w * D * X) / p1
  list(ATT = att, att.inf.func = w * D * (dy - fit - att) / p1 - drop((c0 * (dy - fit) * X) %*% solve(M, mx1)),
    log_weight_gradient = d92_custom_sel_grad(y1, y0, D, covariates, i.weights))
}

d92_custom_sel_grad <- function(y1, y0, D, covariates, i.weights = NULL, ...) {
  w <- i.weights / mean(i.weights)
  dy <- y1 - y0
  X <- as.matrix(covariates)
  c0 <- (1 - D) * w
  beta <- solve(crossprod(X, X * c0), crossprod(X, c0 * dy))
  r <- dy - drop(X %*% beta)
  p1 <- sum(w * D)
  att <- sum(w * D * r) / p1
  # Exact derivative wrt log(i.weights_i) for the weighted regression callback
  g_t <- (w * D / p1) * (r - att)
  M <- crossprod(X, X * c0)
  mx1 <- colSums(w * D * X) / p1
  g_c <- -drop((c0 * r * X) %*% solve(M, mx1))
  g_t + g_c
}

# Namespace clone that instruments producer chunk calls, replay calls, and captures
# the underlying `run_chained_ipw_pipeline` output while supporting forced legacy
# multi-block planning, forced streaming, or custom retention caps.
d92_counting_ns <- function(legacy_plan = FALSE, force_stream = FALSE, cap = NULL) {
  e <- new.env(parent = api)
  for (nm in c("att_gt", "chained_public_att_gt", "run_chained_ipw_pipeline",
               "cp04_stream_covariance", "cp04_stream_plan", "cp04_stream_plan_full_block",
               "cp04_retained_if_reserve", "cp04_retained_if_phase_reserve",
               "cp04_retained_if_validate_backend",
               "cp04_retained_if_check", "cp04_retained_if_replay", "cp04_retained_if_release",
               "cp04_replay_links", "build_cp04_public_registry")) {
    f <- get(nm, api)
    environment(f) <- e
    assign(nm, f, e)
  }
  n <- new.env(parent = emptyenv())
  n$chunks <- 0L
  n$replays <- 0L
  n$rows <- list()
  n$last_pipeline <- NULL
  e$chained_full_if_chunk <- function(backend, rows, links, scale, plan = NULL) {
    n$chunks <- n$chunks + 1L
    n$rows[[length(n$rows) + 1L]] <- rows
    api$chained_full_if_chunk(backend, rows, links, scale, plan)
  }
  replay_fn <- get("cp04_replay_links", e)
  e$cp04_replay_links <- function(backend, links, callback, scale = "core", chunk_rows = NULL) {
    n$replays <- n$replays + 1L
    replay_fn(backend, links, callback, scale = scale, chunk_rows = chunk_rows)
  }
  pipe_fn <- get("run_chained_ipw_pipeline", e)
  e$run_chained_ipw_pipeline <- function(...) {
    res <- pipe_fn(...)
    n$last_pipeline <- res
    res
  }
  if (legacy_plan) e$cp04_stream_plan_full_block <- function(...) NULL
  if (force_stream) {
    e$cp04_retained_if_reserve <- function(...) NULL
    e$cp04_retained_if_phase_reserve <- function(...) NULL
  }
  if (!is.null(cap)) {
    reserve_fn <- get("cp04_retained_if_reserve", e)
    e$cp04_retained_if_reserve <- function(..., control) {
      control$max_retained_entries <- cap
      reserve_fn(..., control = control)
    }
  }
  e$n <- n
  e
}

d92_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (is.list(x)) {
    a <- attributes(x)
    x <- lapply(x, d92_strip)
    attributes(x) <- a
  }
  x
}

d92_find_store <- function(x, depth = 0L) {
  if (depth > 12L) return(FALSE)
  if (is.environment(x)) return(exists("kind", envir = x, inherits = FALSE) &&
    identical(get("kind", envir = x, inherits = FALSE), "cp04_retained_if"))
  if (inherits(x, "cdid_cp04_retained_if_bound")) return(TRUE)
  if (is.list(x)) for (v in x) if (d92_find_store(v, depth + 1L)) return(TRUE)
  FALSE
}

d92_expect_parity <- function(a, b, info = NULL) {
  expect_equal(a$covariance_stream$mhat, b$covariance_stream$mhat, tolerance = 1e-12, info = info)
  expect_equal(a$covariance_stream$V, b$covariance_stream$V, tolerance = 1e-12, info = info)
  expect_equal(a$covariance_stream$mean_phi, b$covariance_stream$mean_phi, tolerance = 1e-12, info = info)
  expect_equal(a$gmm_fit$targets, b$gmm_fit$targets, tolerance = 1e-12, info = info)
  expect_equal(a$gmm_fit$if_projection$sandwich_cov, b$gmm_fit$if_projection$sandwich_cov, tolerance = 1e-12, info = info)
  expect_equal(a$gmm_fit$if_projection$materialized_target_if, b$gmm_fit$if_projection$materialized_target_if,
    tolerance = 1e-12, info = info)
  expect_equal(a$full_frame_transport$Phi_full, b$full_frame_transport$Phi_full, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$estimates, b$reporting$estimates, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$se, b$reporting$se, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$analytical_se, b$reporting$analytical_se, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$if_full, b$reporting$if_full, tolerance = 1e-12, info = info)
  expect_equal(a$reporting$vcov, b$reporting$vcov, tolerance = 1e-12, info = info)
  expect_equal(d92_strip(a$auxiliary_reporting), d92_strip(b$auxiliary_reporting), tolerance = 1e-12, info = info)
  for (ag in names(a$aggregations)) {
    expect_equal(a$aggregations[[ag]]$overall$estimate, b$aggregations[[ag]]$overall$estimate, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$overall$se, b$aggregations[[ag]]$overall$se, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$estimates, b$aggregations[[ag]]$estimates, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$se, b$aggregations[[ag]]$se, tolerance = 1e-12, info = info)
    expect_equal(a$aggregations[[ag]]$if_full, b$aggregations[[ag]]$if_full, tolerance = 1e-12, info = info)
  }
  expect_identical(a$resource_accounting$simultaneous_entries, b$resource_accounting$simultaneous_entries, info = info)
  expect_identical(a$resource_accounting$point_producer$chunk_rows, b$resource_accounting$point_producer$chunk_rows, info = info)
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches, info = info)
}

d92_run_public <- function(d, chunk = 5L, ns = api, ...) {
  ns$att_gt(
    yname = "y", tname = "time", idname = "id", gname = "G", data = d,
    bstrap = FALSE, cband = FALSE, control = list(chunk_rows = chunk), ...
  )
}

# ---------------------------------------------------------------------------
# 1. Route coverage: full-width covariance planner and retained moment IF across
#    selection IPW/OR/DR/custom, imputation IPW/reg/DR, augmented double-DR,
#    dated reference-input routes, missingness, non-unit weights, F > A, and
#    clustered inference.
# ---------------------------------------------------------------------------
test_that("full-width covariance and retained moment IF cover every supported CP04 producer with 1e-12 parity", {
  d_static <- d92_route_panel(vary_x = FALSE, vary_w = FALSE, drop_y = TRUE)
  d_dated <- d92_route_panel(vary_x = TRUE, vary_w = TRUE, drop_y = TRUE)
  fr_f_gt_a <- unique(d_static[, c("id", "G", "x", "w", "rho", "cl")])
  fr_f_gt_a <- rbind(
    fr_f_gt_a[order(fr_f_gt_a$id), ],
    data.frame(id = c(901L, 902L), G = c(4, 0), x = c(0.15, -0.25), w = c(1.1, 0.9), rho = c(0.75, 0.75), cl = c(13L, 13L))
  )

  cases <- list(
    list(nm = "baseline_or_missing_weights_cluster", d = d_static,
         args = list(xformla = ~x, weightsname = "w", clustervars = "cl", est_method = "reg")),
    list(nm = "baseline_dr_missing_f_gt_a", d = d_static,
         args = list(xformla = ~x, weightsname = "w", frame = fr_f_gt_a[, c("id", "G", "x", "w")], est_method = "dr")),
    list(nm = "selection_ipw_pair_within", d = d_static,
         args = list(xformla = ~x, weightsname = "w", est_method = "ipw",
                     selection = "conditional_mtar", selection_model = "pair_logit", selection_pooling = "within_cohort")),
    list(nm = "selection_or_pair_across_U_B", d = d_static,
         args = list(xformla = ~x, weightsname = "w", est_method = "reg",
                     selection = "conditional_mtar", selection_model = "pair_logit",
                     selection_pooling = "across_cohorts", selection_universe = "B", selection_rho = "rho",
                     selection_justification = "known unit transport ratio rho")),
    # Date factorization is supported on F only; B transport is tested above
    # with the direct-pair model, without inventing a factorization on B.
    list(nm = "selection_or_period_across_F", d = d_static,
         args = list(xformla = ~x, weightsname = "w", est_method = "reg",
                     selection = "conditional_mtar", selection_model = "independent_period_logit",
                     selection_pooling = "across_cohorts")),
    list(nm = "selection_dr_pair_separate_cluster", d = d_static,
         args = list(xformla = ~x, weightsname = "w", clustervars = "cl", est_method = "dr",
                     selection = "conditional_mtar", selection_model = "pair_logit", selection_pooling = "separate")),
    list(nm = "selection_custom_callback", d = d_static,
         args = list(xformla = ~x, weightsname = "w", est_method = d92_custom_reg,
                     selection = "conditional_mtar", selection_model = "pair_logit", selection_pooling = "within_cohort")),
    list(nm = "imputation_ipw_within", d = d_static,
         args = list(xformla = ~x, weightsname = "w", est_method = "ipw",
                     selection = "conditional_mtar", selection_est_method = "imputation", imputation_pooling = "within_cohort")),
    list(nm = "imputation_reg_across", d = d_static,
         args = list(xformla = ~x, weightsname = "w", est_method = "reg",
                     selection = "conditional_mtar", selection_est_method = "imputation", imputation_pooling = "across_cohorts")),
    list(nm = "imputation_dr_within_cluster", d = d_static,
         args = list(xformla = ~x, weightsname = "w", clustervars = "cl", est_method = "dr",
                     selection = "conditional_mtar", selection_est_method = "imputation", imputation_pooling = "within_cohort")),
    list(nm = "augmented_double_dr_within", d = d_static,
         args = list(xformla = ~x, weightsname = "w", est_method = "dr",
                     selection = "conditional_mtar", selection_est_method = "dr",
                     selection_model = "pair_logit", selection_pooling = "within_cohort",
                     imputation_pooling = "within_cohort")),
    list(nm = "dated_reference_ipw", d = d_dated,
         args = list(xformla = ~x, weightsname = "w", est_method = "ipw", fix_weights = "first_period")),
    list(nm = "dated_reference_custom", d = d_dated,
         args = list(xformla = ~x, weightsname = "w", est_method = d92_custom_reg, fix_weights = "base_period")),
    list(nm = "dated_reference_dr", d = d_dated,
         args = list(xformla = ~x, weightsname = "w", est_method = "dr", fix_weights = "first_period")),
    list(nm = "dated_selection_dr", d = d_dated,
         args = list(xformla = ~x, weightsname = "w", est_method = "dr", fix_weights = "first_period",
                     selection = "conditional_mtar", selection_model = "pair_logit", selection_pooling = "within_cohort")),
    list(nm = "dated_imputation_dr", d = d_dated,
         args = list(xformla = ~x, weightsname = "w", est_method = "dr", fix_weights = "base_period",
                     selection = "conditional_mtar", selection_est_method = "imputation", imputation_pooling = "within_cohort")),
    list(nm = "dated_augmented_dr", d = d_dated,
         args = list(xformla = ~x, weightsname = "w", est_method = "dr", fix_weights = "first_period",
                     selection = "conditional_mtar", selection_est_method = "dr",
                     selection_model = "pair_logit", selection_pooling = "within_cohort",
                     imputation_pooling = "within_cohort"))
  )

  for (cs in cases) {
    e_ret <- d92_counting_ns()
    e_str <- d92_counting_ns(force_stream = TRUE)
    e_leg <- d92_counting_ns(legacy_plan = TRUE, force_stream = TRUE)

    fit_ret <- do.call(d92_run_public, c(list(d = cs$d, chunk = 5L, ns = e_ret), cs$args))
    fit_str <- do.call(d92_run_public, c(list(d = cs$d, chunk = 5L, ns = e_str), cs$args))
    fit_leg <- do.call(d92_run_public, c(list(d = cs$d, chunk = 5L, ns = e_leg), cs$args))

    res_ret <- e_ret$n$last_pipeline
    res_str <- e_str$n$last_pipeline
    res_leg <- e_leg$n$last_pipeline
    rec <- res_ret$resource_accounting$point_producer$retained_moment_if

    expect_identical(rec$status, "retained", info = cs$nm)
    expect_identical(rec$producer_passes, 1L, info = cs$nm)
    expect_identical(rec$slice_passes, 3L, info = cs$nm)
    expect_identical(rec$auxiliary$status, "not_extended", info = cs$nm)
    expect_identical(res_str$resource_accounting$point_producer$retained_moment_if$status, "streamed", info = cs$nm)
    expect_identical(res_leg$resource_accounting$point_producer$retained_moment_if$status, "streamed", info = cs$nm)
    expect_false(d92_find_store(res_ret), info = cs$nm)
    expect_false(d92_find_store(fit_ret), info = cs$nm)

    # The clone counts the four main replays (mean, Gram, 2 GMM passes).
    # Auxiliary reporting runs in the installed namespace and is not counted;
    # retained producer chunks here are exactly the one fill pass.
    B <- res_ret$resource_accounting$point_producer$chunk_rows
    N <- res_ret$nuisance_registry$cp04_reader$state$N
    n_batches <- as.integer(ceiling(N / B))
    expect_identical(e_ret$n$replays, 4L, info = cs$nm)
    expect_identical(e_str$n$replays, 4L, info = cs$nm)
    expect_gt(e_leg$n$replays, 4L, label = cs$nm)
    expect_lt(e_ret$n$chunks, e_str$n$chunks, label = cs$nm)
    expect_identical(rec$n_batches, n_batches, info = cs$nm)

    d92_expect_parity(res_ret, res_str, info = paste(cs$nm, "retained vs stream"))
    d92_expect_parity(res_ret, res_leg, info = paste(cs$nm, "retained vs legacy"))
    expect_equal(fit_ret$att, fit_str$att, tolerance = 1e-12, info = cs$nm)
    expect_equal(fit_ret$se, fit_str$se, tolerance = 1e-12, info = cs$nm)
    expect_equal(fit_ret$inffunc, fit_str$inffunc, tolerance = 1e-12, info = cs$nm)
  }
})

# ---------------------------------------------------------------------------
# 2. Beta-only OR-core optimization and details = TRUE vs FALSE exact equality
#    on corrected DR (cdid_dr_corollary1) and dated baseline DR (cdid_dr_baseline).
# ---------------------------------------------------------------------------
test_that("non-detailed corrected DR and dated baseline DR use beta-only OR core and match detailed IFs", {
  d_dated <- d92_route_panel(vary_x = TRUE, vary_w = TRUE, drop_y = TRUE)

  for (has_x in c(TRUE, FALSE)) {
    xf <- if (has_x) ~x else ~1
    # 1. Dated baseline DR
    e_b <- d92_counting_ns()
    d92_run_public(d_dated, chunk = 12L, ns = e_b, xformla = xf, weightsname = "w",
                   est_method = "dr", fix_weights = "first_period")
    b_base <- e_b$n$last_pipeline$nuisance_registry$cp04
    expect_s3_class(b_base, "cdid_dr_baseline")
    links_b <- seq_along(b_base$links)
    env_b_f <- api$dr_baseline_if_envelope(b_base, 12L, links_b, details = FALSE)
    env_b_t <- api$dr_baseline_if_envelope(b_base, 12L, links_b, details = TRUE)
    expect_identical(env_b_f$general$or_core_evaluated, has_x)
    expect_true(env_b_t$general$or_core_evaluated)

    ch_b_f <- api$evaluate_dr_baseline_if_chunk(b_base, 1:12, links_b, scale = "frame", details = FALSE)
    ch_b_t <- api$evaluate_dr_baseline_if_chunk(b_base, 1:12, links_b, scale = "frame", details = TRUE)
    expect_equal(ch_b_f$phi, ch_b_t$phi, tolerance = 1e-14)
    expect_equal(ch_b_t$details$psi_dir + ch_b_t$details$psi_or + ch_b_t$details$psi_ps + ch_b_t$details$psi_sel,
                 ch_b_t$phi, tolerance = 1e-14)
    if (!has_x) {
      expect_identical(ch_b_f$diagnostics$n_regression_score_evaluations, 0L)
      expect_identical(ch_b_f$diagnostics$n_global_beta_blocks, 0L)
      expect_gt(ch_b_t$diagnostics$n_regression_score_evaluations, 0L)
    } else {
      expect_gt(ch_b_f$diagnostics$n_global_beta_blocks, 0L)
    }

    # 2. Corrected DR (both static and dated reference inputs)
    for (dated_cor in c(FALSE, TRUE)) {
      dat <- if (dated_cor) d_dated else d92_route_panel(vary_x = FALSE, vary_w = FALSE, drop_y = TRUE)
      e_c <- d92_counting_ns()
      d92_run_public(dat, chunk = 12L, ns = e_c, xformla = xf, weightsname = "w",
                     est_method = "dr", selection = "conditional_mtar",
                     selection_model = "pair_logit", selection_pooling = "within_cohort",
                     fix_weights = if (dated_cor) "first_period" else NULL)
      b_cor <- e_c$n$last_pipeline$nuisance_registry$cp04
      expect_s3_class(b_cor, "cdid_dr_corollary1")
      links_c <- seq_along(b_cor$links)

      ch_c_f <- api$evaluate_dr_corollary1_if_chunk(b_cor, 1:12, links_c, scale = "frame", details = FALSE)
      ch_c_t <- api$evaluate_dr_corollary1_if_chunk(b_cor, 1:12, links_c, scale = "frame", details = TRUE)
      expect_equal(ch_c_f$phi, ch_c_t$phi, tolerance = 1e-14)
      expect_equal(ch_c_t$details$psi_dir + ch_c_t$details$psi_or + ch_c_t$details$psi_ps + ch_c_t$details$psi_sel,
                   ch_c_t$phi, tolerance = 1e-14)
      expect_true(is.matrix(ch_c_t$details$alpha))
      expect_true(is.matrix(ch_c_t$details$psi_or_estimator))
      if (!has_x) {
        expect_identical(ch_c_f$diagnostics$n_regression_score_evaluations, 0L)
        expect_identical(ch_c_f$diagnostics$n_global_beta_blocks, 0L)
        expect_gt(ch_c_t$diagnostics$n_regression_score_evaluations, 0L)
      } else {
        expect_gt(ch_c_f$diagnostics$n_global_beta_blocks, 0L)
      }
    }
  }
})

# ---------------------------------------------------------------------------
# 3. Tight-cap decline, boundary admission, and LAST callback mutation refusal
#    for generalized producers and route-specific owners.
# ---------------------------------------------------------------------------
test_that("tight-cap decline, error-path release, and LAST callback mutation checks hold", {
  d_static <- d92_route_panel(vary_x = FALSE, vary_w = FALSE, drop_y = TRUE)

  for (rt in list(
    list(nm = "selection_dr", args = list(xformla = ~x, weightsname = "w", est_method = "dr",
                                          selection = "conditional_mtar", selection_model = "pair_logit",
                                          selection_pooling = "within_cohort")),
    list(nm = "imputation_dr", args = list(xformla = ~x, weightsname = "w", est_method = "dr",
                                           selection = "conditional_mtar", selection_est_method = "imputation",
                                           imputation_pooling = "within_cohort")),
    list(nm = "augmented_dr", args = list(xformla = ~x, weightsname = "w", est_method = "dr",
                                          selection = "conditional_mtar", selection_est_method = "dr",
                                          selection_model = "pair_logit", selection_pooling = "within_cohort",
                                          imputation_pooling = "within_cohort"))
  )) {
    e0 <- d92_counting_ns()
    do.call(d92_run_public, c(list(d = d_static, chunk = 5L, ns = e0), rt$args))
    res0 <- e0$n$last_pipeline
    rec0 <- res0$resource_accounting$point_producer$retained_moment_if
    total <- rec0$admission_total_entries
    expect_identical(rec0$status, "retained", info = rt$nm)

    # Exact cap admits; one entry below declines cleanly to ordinary streaming.
    e_at <- d92_counting_ns(cap = total)
    e_lo <- d92_counting_ns(cap = total - 1)
    do.call(d92_run_public, c(list(d = d_static, chunk = 5L, ns = e_at), rt$args))
    do.call(d92_run_public, c(list(d = d_static, chunk = 5L, ns = e_lo), rt$args))
    expect_identical(e_at$n$last_pipeline$resource_accounting$point_producer$retained_moment_if$status, "retained", info = rt$nm)
    expect_identical(e_lo$n$last_pipeline$resource_accounting$point_producer$retained_moment_if$status, "streamed", info = rt$nm)
    d92_expect_parity(res0, e_at$n$last_pipeline, info = paste(rt$nm, "at cap"))
    d92_expect_parity(res0, e_lo$n$last_pipeline, info = paste(rt$nm, "below cap"))

    # LAST callback mutation during producer fill pass and during slice pass is refused.
    b <- res0$nuisance_registry$cp04
    links <- vapply(res0$point_links_used, `[[`, integer(1), "cp04_index")
    B <- res0$resource_accounting$point_producer$chunk_rows
    N <- res0$nuisance_registry$cp04_reader$state$N
    n_batches <- as.integer(ceiling(N / B))
    plan <- api$cp04_stream_plan(b, links, 0)

    # Fill-pass LAST callback mutation of store state or route owner
    store_fill <- api$cp04_retained_if_reserve(b, links, B, plan, 0, 0, 0, 0, b$control)
    bound_fill <- api$cp04_retained_if_bind(b, store_fill)
    seen_fill <- 0L
    expect_error(api$cp04_replay_links(bound_fill, links, function(bt) {
      seen_fill <<- seen_fill + 1L
      if (seen_fill == n_batches) store_fill$state <- "tampered"
    }, chunk_rows = B), class = "cdid_cp04_integrity_error")
    expect_identical(seen_fill, n_batches)
    expect_identical(store_fill$state, "released")

    # Slice-pass LAST callback mutation of locked raw phi or route owner snapshot
    store_slice <- api$cp04_retained_if_reserve(b, links, B, plan, 0, 0, 0, 0, b$control)
    bound_slice <- api$cp04_retained_if_bind(b, store_slice)
    api$cp04_replay_links(bound_slice, links, function(bt) NULL, chunk_rows = B)
    expect_identical(store_slice$state, "filled")
    seen_slice <- 0L
    expect_error(api$cp04_replay_links(bound_slice, links, function(bt) {
      seen_slice <<- seen_slice + 1L
      if (seen_slice == n_batches) store_slice$phi <- store_slice$phi + 1
    }, chunk_rows = B), class = "cdid_cp04_integrity_error")
    expect_identical(seen_slice, n_batches)
    expect_identical(store_slice$state, "released")

    # Release cleans all snapshots
    rel <- api$cp04_retained_if_release(store_slice)
    expect_identical(store_slice$state, "released")
    expect_false(exists("phi", envir = store_slice, inherits = FALSE))
    expect_false(exists("selection_snapshot", envir = store_slice, inherits = FALSE))
    expect_false(exists("design_snapshot", envir = store_slice, inherits = FALSE))
    expect_false(exists("layout_snapshot", envir = store_slice, inherits = FALSE))
  }
})

# ---------------------------------------------------------------------------
# 4. Optional NULL owner: the dated augmented layout carries no state (its
#    sealed layout reader does). The NULL layout-state snapshot fills and
#    slices; a LAST callback that gives it a state is still refused.
# ---------------------------------------------------------------------------
test_that("dated augmented NULL layout state is retained and its LAST callback mutation is refused", {
  d_dated <- d92_route_panel(vary_x = TRUE, vary_w = TRUE, drop_y = TRUE)
  e0 <- d92_counting_ns()
  d92_run_public(d_dated, chunk = 5L, ns = e0, xformla = ~x, weightsname = "w", est_method = "dr",
                 fix_weights = "first_period", selection = "conditional_mtar", selection_est_method = "dr",
                 selection_model = "pair_logit", selection_pooling = "within_cohort",
                 imputation_pooling = "within_cohort")
  res0 <- e0$n$last_pipeline
  expect_identical(res0$resource_accounting$point_producer$retained_moment_if$status, "retained")
  b <- res0$nuisance_registry$cp04
  expect_s3_class(b, "cdid_augmented_backend")
  expect_false(is.null(b$layout))
  expect_null(b$layout$state)

  links <- vapply(res0$point_links_used, `[[`, integer(1), "cp04_index")
  B <- res0$resource_accounting$point_producer$chunk_rows
  plan <- api$cp04_stream_plan(b, links, 0)
  store <- api$cp04_retained_if_reserve(b, links, B, plan, 0, 0, 0, 0, b$control)
  n_batches <- as.integer(ceiling(store$N / B))
  expect_null(store$layout_state_snapshot)
  bound <- api$cp04_retained_if_bind(b, store)
  api$cp04_replay_links(bound, links, function(bt) NULL, chunk_rows = B)
  expect_identical(store$state, "filled")

  seen <- 0L
  expect_error(api$cp04_replay_links(bound, links, function(bt) {
    seen <<- seen + 1L
    if (seen == n_batches) {
      unlockBinding("layout_state_snapshot", store)
      assign("layout_state_snapshot", get("cp04_reader", b$reader_ref, inherits = FALSE)$state, envir = store)
    }
  }, chunk_rows = B), "sealed layout owner", class = "cdid_cp04_integrity_error")
  expect_identical(seen, n_batches)
  expect_identical(store$state, "released")
})
