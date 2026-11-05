# Internal test suite for the D93 standalone OR-free compact raw-IPW backend
# (`R/chained_ipw_raw_backend.R` and private adapters in `R/chained_or_public.R`,
# `R/chained_cp04_public.R`, and `R/chained_ipw_pipeline.R`).

make_d93_compact_ipw_test_panel <- function(
  n_units = 60L,
  times = c(1, 2, 3, 4),
  cohorts = c(2, 3, Inf),
  missing_mode = c("partial", "complete"),
  extreme_z = FALSE,
  unequal_weights = TRUE
) {
  missing_mode <- match.arg(missing_mode)
  unit_ids <- seq_len(n_units)
  g_vec <- cohorts[((unit_ids - 1L) %% length(cohorts)) + 1L]
  z1 <- sin(unit_ids * 0.73) + 0.25 * cos(unit_ids * 0.31)
  if (isTRUE(extreme_z)) {
    # Give treated-cohort and control-cohort units moderate signal plus one
    # outlier in cohort 2 so fitted eta exceeds log(999999) ~ 13.82 (activating
    # D93 odds clipping at 1 - 1e-6) while remaining well below the glm.fit
    # numerical 0/1 boundary (~33.9) and avoiding coordinate separation.
    z1 <- ifelse(g_vec == cohorts[1L], 0.9, -0.6) + 0.8 * sin(unit_ids * 0.73)
    z1[1L] <- 16.0
    z1[2L] <- 12.0
  }
  w_vec <- if (isTRUE(unequal_weights)) {
    0.6 + 0.8 * ((unit_ids %% 7L) + 1L) / 7.0 + 0.15 * abs(cos(unit_ids))
  } else {
    rep(1.0, n_units)
  }
  cl_vec <- ((unit_ids - 1L) %% 12L) + 1L

  rows <- list()
  idx <- 1L
  for (i in seq_len(n_units)) {
    for (t_pos in seq_along(times)) {
      tt <- times[t_pos]
      keep <- TRUE
      if (identical(missing_mode, "partial")) {
        # Deterministic unbalancedness preserving positive support in every (cohort, pair)
        if (tt == times[2L] && (i %% 5L == 0L)) keep <- FALSE
        if (tt == times[3L] && (i %% 7L == 0L)) keep <- FALSE
        if (tt == times[4L] && (i %% 6L == 0L)) keep <- FALSE
      }
      if (keep) {
        y_val <- 1.2 + 0.35 * tt + 0.4 * z1[i] +
          (if (is.finite(g_vec[i]) && tt >= g_vec[i]) 0.85 * (tt - g_vec[i] + 1.0) else 0.0) +
          0.11 * sin(i * 1.1 + tt * 0.7)
        rows[[idx]] <- data.frame(
          id = as.integer(unit_ids[i]),
          time = as.double(tt),
          y = as.double(y_val),
          g = as.double(g_vec[i]),
          x1 = as.double(z1[i]),
          w = as.double(w_vec[i]),
          cluster = as.integer(cl_vec[i]),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  list(
    data = df,
    unit_id = unit_ids,
    unit_g = g_vec,
    unit_x1 = z1,
    unit_w = w_vec,
    unit_cluster = cl_vec
  )
}

make_d93_test_pipeline_control <- function(
  max_chunk_rows = 256L,
  max_observer_units = 10000L,
  max_retained_entries = 5e7,
  max_matrix_entries = 5e7,
  max_work_ops = 5e8,
  max_output_entries = 5e7,
  require_matching_scale = TRUE
) {
  lnk_ctrl <- ipw_link_enclosure_control(
    max_chunk_rows = max_chunk_rows,
    max_retained_entries = max_retained_entries,
    max_matrix_entries = max_matrix_entries,
    max_work_ops = max_work_ops,
    max_output_entries = max_output_entries
  )
  ipw_pipeline_control(
    owner_control = ipw_owner_control(
      max_chunk_rows = max_chunk_rows,
      max_retained_entries = max_retained_entries,
      max_matrix_entries = max_matrix_entries,
      max_work_ops = max_work_ops,
      max_output_entries = max_output_entries
    ),
    catalog_bounds = list(
      max_candidates = 2000L,
      max_rows = 2000L,
      max_pool_entries = 20000L,
      max_scratch = max_matrix_entries,
      max_work = max_work_ops
    ),
    class_bounds = list(
      max_links = 2000L,
      max_classes = 1000L,
      max_class_dates = 200L,
      max_E_entries = 50000L,
      max_mask_entries = 1e7
    ),
    fitted_control = ipw_fitted_view_control(
      max_chunk_rows = max_chunk_rows,
      max_observer_units = max_observer_units,
      max_retained_entries = max_retained_entries,
      max_matrix_entries = max_matrix_entries,
      max_work_ops = max_work_ops,
      max_output_entries = max_output_entries
    ),
    point_control = ipw_point_view_control(
      max_chunk_rows = max_chunk_rows,
      max_observer_units = max_observer_units,
      max_retained_entries = max_retained_entries,
      max_matrix_entries = max_matrix_entries,
      max_work_ops = max_work_ops,
      max_output_entries = max_output_entries
    ),
    glm_control = list(epsilon = 1e-8, maxit = 100L, trace = FALSE),
    certificate_control = logit_certificate_control(
      taylor_degree = 12L,
      max_abs_eta = 50.0,
      max_squarings = 16L,
      max_chunk_rows = max_chunk_rows,
      max_matrix_entries = max_matrix_entries,
      max_work_ops = max_work_ops
    ),
    problem_control = ipw_nuisance_problem_control(
      max_chunk_rows = max_chunk_rows,
      max_retained_entries = max_retained_entries,
      max_matrix_entries = max_matrix_entries,
      max_work_ops = max_work_ops
    ),
    link_control = lnk_ctrl,
    covariance_control = ipw_covariance_transfer_control(
      max_chunk_rows = max_chunk_rows,
      max_operator_block_rows = max_chunk_rows,
      max_cluster_rows = 20000L,
      max_matrix_entries = max_matrix_entries,
      max_work_ops = max_work_ops,
      max_output_entries = max_output_entries,
      link_control = lnk_ctrl
    ),
    gmm_accuracy = gmm_accuracy_control(
      method = "conventional_chol_qr",
      target_error_unit = "abs",
      max_target_error_abs = 1e-4,
      max_target_error_se = 1e-2,
      max_orthogonality_defect = 1e-8,
      min_certified_lambda_rel = 1e-10,
      max_factor_residual_rel = 1e-10,
      min_rcond_R = 1e-10,
      min_rcond_U = 1e-10,
      min_sigma_rel_A = 1e-10,
      max_kkt_residual_rel = 1e-10,
      max_loading_residual_rel = 1e-10,
      require_matching_scale = require_matching_scale,
      max_loading_error_abs = 1e-4,
      covariance_strategy = gmm_covariance_strategy(
        method = "tikhonov",
        weighting_requested = if (require_matching_scale) "optimal" else "identity",
        weighting_canonical = if (require_matching_scale) "optimal" else "identity"
      )
    ),
    gmm_bounds = gmm_bounds_control(
      max_matrix_entries = max_matrix_entries,
      max_batch_rows = max_chunk_rows,
      max_targets = 500L,
      max_open_clusters = 20000L,
      allow_materialize_if = TRUE,
      max_materialized_if_entries = max_output_entries
    ),
    multiplier_budget = list(
      max_units = 20000L,
      max_clusters = 20000L,
      max_targets = 500L,
      max_draws = 2000L,
      max_if_batch_rows = max_chunk_rows,
      max_if_batch_cells = max_matrix_entries,
      max_sign_batch_rows = 64L,
      max_sign_batch_cells = max_matrix_entries,
      max_workspace_cells = max_matrix_entries,
      max_workspace_bytes = 8 * max_matrix_entries
    ),
    inference_accuracy = list(
      max_estimate_error_abs = 1e-3,
      max_se_error_abs = 1e-2,
      max_vcov_error_abs = 1e-2,
      max_multiplier_error_abs = 1e-1,
      max_band_error_abs = 1e-1
    )
  )
}

# Independent scalar mathematical oracle for a single (pair, g, ordered K) link
# operating directly on `panel_owner` and `frame_owner` without calling any
# function in `R/chained_ipw_raw_backend.R`.
scalar_oracle_ipw_link <- function(frame_owner, panel_owner, g, K, s, t, glm_ctrl) {
  n_a <- as.integer(frame_owner$n_a_units)
  w_a <- as.double(frame_owner$w_a)
  Z_a <- frame_owner$Z_a
  a_g <- as.double(panel_owner$a_g)
  cohorts <- c(as.double(g), as.double(K))
  m <- length(cohorts)

  ds <- match(as.double(s), as.double(panel_owner$all_times))
  dt <- match(as.double(t), as.double(panel_owner$all_times))
  post <- panel_owner$postings
  idx_s <- seq.int(post$date_start[ds], length.out = post$date_count[ds])
  idx_t <- seq.int(post$date_start[dt], length.out = post$date_count[dt])
  a_s <- post$post_a_index[idx_s]
  a_t <- post$post_a_index[idx_t]
  obs_s_map <- integer(n_a); obs_s_map[a_s] <- post$post_obs_index[idx_s]
  obs_t_map <- integer(n_a); obs_t_map[a_t] <- post$post_obs_index[idx_t]
  pair_obs <- (obs_s_map > 0L) & (obs_t_map > 0L)

  eps <- 1e-6
  cap <- 1.0 - 1e-6
  cap_odds <- (1.0 - 1e-6) / (1.0 - (1.0 - 1e-6))

  gi <- match(a_g, cohorts)
  q_vec <- numeric(m)
  r_vec <- numeric(m)
  for (j in seq_len(m)) {
    in_c <- which(gi == j)
    obs_c <- in_c[pair_obs[in_c]]
    mass_a_j <- sum(w_a[in_c])
    mass_p_j <- sum(w_a[obs_c])
    # r_j is a ratio of empirical weighted moments. Its unit IF denominator
    # is mean_A(w * 1{G=j}), not a share divided by the total weight mass.
    q_vec[j] <- mass_a_j / n_a
    r_vec[j] <- mass_p_j / mass_a_j
  }
  r_clip <- pmax(r_vec, eps)
  r_act <- as.double(r_vec > eps)

  risk_obs <- which(!is.na(gi) & pair_obs)
  gi_obs <- gi[risk_obs]
  w_obs <- w_a[risk_obs]
  D_obs <- as.double(gi_obs == 1L)
  fit_w <- w_obs / r_clip[gi_obs]
  Z_obs <- Z_a[risk_obs, , drop = FALSE]

  fit <- suppressWarnings(stats::glm.fit(
    x = Z_obs,
    y = D_obs,
    weights = fit_w,
    start = rep(0, ncol(Z_obs)),
    family = stats::binomial(link = "logit"),
    control = list(epsilon = glm_ctrl$epsilon, maxit = glm_ctrl$maxit, trace = FALSE),
    intercept = FALSE
  ))
  beta <- as.double(fit$coefficients)

  # Unclipped raw score/Hessian and D93 capped odds at beta
  eta_raw <- drop(Z_obs %*% beta)
  hi <- is.finite(eta_raw) & eta_raw > 40
  odds_raw <- exp(ifelse(hi, 0, eta_raw))
  p_raw <- odds_raw / (1.0 + odds_raw)
  p_raw[hi] <- 1.0 / (1.0 + exp(-eta_raw[hi]))
  active_odds <- as.double(p_raw < cap)
  odds_use <- ifelse(p_raw < cap, odds_raw, cap_odds)
  prob_use <- p_raw

  # Preserve the scalar raw-information multiplication order.
  w_hess <- fit_w * prob_use * (1.0 - prob_use)
  H <- crossprod(Z_obs, Z_obs * as.numeric(w_hess)) / n_a
  A_beta_r <- matrix(0.0, nrow = ncol(Z_obs), ncol = m)
  for (j in seq_len(m)) {
    if (r_act[j] > 0) {
      sel_j <- which(gi_obs == j)
      s_ij <- w_obs[sel_j] * (D_obs[sel_j] - prob_use[sel_j]) / (r_clip[j]^2)
      A_beta_r[, j] <- -colSums(Z_obs[sel_j, , drop = FALSE] * as.numeric(s_ij)) / n_a
    }
  }

  is_t <- gi_obs == 1L
  is_c <- !is_t
  a_obs <- fit_w[is_t]
  b_obs <- fit_w[is_c] * odds_use[is_c]
  sum_a <- sum(a_obs)
  sum_b <- sum(b_obs)
  mean_a <- sum_a / n_a
  mean_b <- sum_b / n_a

  dy_obs <- as.double(panel_owner$observations$y[obs_t_map[risk_obs]]) -
    as.double(panel_owner$observations$y[obs_s_map[risk_obs]])
  dy_t <- dy_obs[is_t]
  dy_c <- dy_obs[is_c]
  Tmean <- sum(a_obs * dy_t) / sum_a
  Cmean <- sum(b_obs * dy_c) / sum_b
  mu <- Tmean - Cmean

  Z_c <- Z_obs[is_c, , drop = FALSE]
  G_beta <- -colSums(Z_c * as.numeric(active_odds[is_c] * b_obs * (dy_c - Cmean))) / sum_b

  G_r <- numeric(m)
  G_r[1L] <- 0.0
  for (j in seq.int(2L, m)) {
    if (r_act[j] > 0) {
      sel_j <- which(gi_obs == j)
      G_r[j] <- sum(w_obs[sel_j] * odds_use[sel_j] * (dy_obs[sel_j] - Cmean)) / (r_clip[j]^2 * sum_b)
    }
  }

  H_chol <- chol(H)
  loading_beta <- drop(backsolve(H_chol, forwardsolve(t(H_chol), G_beta)))
  loading_rate <- as.double(G_r + drop(crossprod(A_beta_r, loading_beta)))
  rate_coef_s0 <- -(r_vec / q_vec) * loading_rate
  rate_coef_s1 <- ((1.0 - r_vec) / q_vec) * loading_rate

  comp_treated <- numeric(n_a)
  comp_control <- numeric(n_a)
  comp_score <- numeric(n_a)
  comp_rate <- numeric(n_a)

  in_risk_all <- which(!is.na(gi))
  comp_rate[in_risk_all] <- w_a[in_risk_all] * rate_coef_s0[gi[in_risk_all]]
  comp_rate[risk_obs] <- w_obs * rate_coef_s1[gi_obs]
  comp_score[risk_obs] <- fit_w * (D_obs - prob_use) * drop(Z_obs %*% loading_beta)
  comp_treated[risk_obs[is_t]] <- (a_obs / mean_a) * (dy_t - Tmean)
  comp_control[risk_obs[is_c]] <- -(w_obs[is_c] * odds_use[is_c] / r_clip[gi_obs[is_c]]) /
    mean_b * (dy_c - Cmean)
  phi <- comp_treated + comp_control + comp_score + comp_rate

  list(
    estimate = mu,
    treated_mean = Tmean,
    control_mean = Cmean,
    mean_a = mean_a,
    mean_b = mean_b,
    beta = beta,
    H = H,
    A_beta_r = A_beta_r,
    G_beta = G_beta,
    G_r = G_r,
    loading_beta = loading_beta,
    loading_rate = loading_rate,
    rate_coef_s0 = rate_coef_s0,
    rate_coef_s1 = rate_coef_s1,
    components = list(
      treated = comp_treated,
      control = comp_control,
      beta_score = comp_score,
      rate = comp_rate
    ),
    phi = phi
  )
}

# Run legacy scalar/view pipeline route via a cloned environment with a spy hook
run_legacy_scalar_ipw_pipeline <- function(...) {
  fn <- run_chained_ipw_pipeline
  env <- new.env(parent = environment(fn))
  env$ipw_pipeline_batch_eligibility <- ipw_pipeline_batch_eligibility
  environment(fn) <- env
  fn(...)
}

test_that("D93 compact raw-IPW backend matches independent scalar oracle and scalar reference at 1e-12 across weights, rates, NYT/NT, PT-post/all, and base periods", {
  for (miss in c("partial", "complete")) {
    sim <- make_d93_compact_ipw_test_panel(
      n_units = 60L,
      times = c(1, 2, 3, 4),
      cohorts = c(2, 3, Inf),
      missing_mode = miss,
      unequal_weights = TRUE
    )
    ctrl <- make_d93_test_pipeline_control()
    p_owner <- construct_ipw_panel_owner(sim$data, "id", "time", "y", "g", ctrl$owner_control)
    f_owner <- construct_ipw_frame_owner(
      p_owner,
      id = sim$unit_id,
      Z = cbind(`(Intercept)` = 1.0, x1 = sim$unit_x1),
      weights = sim$unit_w,
      control = ctrl$owner_control,
      frame_g = sim$unit_g
    )
    ledger <- ipw_pipeline_request_ledger(c(2, 3), p_owner$all_times, p_owner$cohorts, 0L, "notyettreated")
    aux_plan <- ipw_auxiliary_spec(ledger, p_owner$all_times, p_owner$cohorts, 0L, "notyettreated")
    meta <- build_owner_comparison_catalog(
      panel_owner = p_owner,
      targets = c(2, 3),
      anticipation = 0L,
      control_group = "notyettreated",
      conventions = "D10_D13_admissible_frame",
      control = ctrl$owner_control,
      catalog_bounds = ctrl$catalog_bounds,
      frame_owner = f_owner,
      build_design = TRUE,
      pt_assumption = "all", pt_start = p_owner$all_times[1L]
    )

    reg <- build_ipw_raw_backend_registry(
      metadata_obj = meta,
      frame_owner = f_owner,
      design = meta$design,
      owner_control = ctrl$owner_control,
      fitted_control = ctrl$fitted_control,
      point_control = ctrl$point_control,
      glm_control = ctrl$glm_control,
      class_bounds = ctrl$class_bounds,
      weighting = "identity",
      reduce_structural_redundancy = FALSE,
      aux_spec = aux_plan
    )
    b <- reg$raw_ipw
    expect_s3_class(b, "cdid_ipw_raw_backend")
    expect_identical(names(reg$link_views), paste0("link:r", seq_len(nrow(meta$catalog$rows))))
    expect_identical(vapply(reg$link_views, `[[`, character(1L), "functional_key"),
      stats::setNames(meta$row_mapping$functional_key, names(reg$link_views)))
    expect_null(b$reader_ref)
    expect_false(b$accounting$full_n_masks_retained)
    expect_false(b$accounting$observer_arrays_retained)

    M <- nrow(meta$catalog$rows)
    ch <- evaluate_ipw_raw_if_chunk(
      b,
      rows = seq_len(f_owner$n_a_units),
      links = seq_len(M),
      scale = "core",
      details = TRUE
    )
    reverse_chunk <- evaluate_ipw_raw_if_chunk(b, rev(seq_len(f_owner$n_a_units)), rev(seq_len(M)))
    expect_equal(unname(reverse_chunk$phi), unname(ch$phi[rev(seq_len(f_owner$n_a_units)), rev(seq_len(M)), drop = FALSE]), tolerance = 1e-12)

    ref_links <- vector("list", M)
    for (m in seq_len(M)) {
      sp <- reg$row_specs[[m]]
      ora <- scalar_oracle_ipw_link(f_owner, p_owner, sp$g, sp$K, sp$s, sp$t, ctrl$glm_control)
      nv_ref <- construct_ipw_nuisance_view(f_owner, sp$g, sp$K, sp$s, sp$t, ctrl$owner_control)
      fv_ref <- fit_ipw_nuisance_view(nv_ref, ctrl$glm_control, ctrl$fitted_control)
      lv_ref <- construct_ipw_link_view(f_owner, sp$s, sp$t, sp$g, sp$K, nv_ref, ctrl$owner_control, link_id = sp$link_id)
      pl_ref <- fit_ipw_point_link_view(fv_ref, lv_ref, ctrl$point_control)
      if_ref <- evaluate_ipw_point_if_view_chunk(pl_ref, 1L, f_owner$n_a_units, ctrl$point_control, include_components = TRUE)
      ref_links[[m]] <- pl_ref

      pl_new <- reg$point_links_all[[m]]
      est_info <- ipw_raw_link_estimate(b, m, details = TRUE)

      expect_equal(est_info$estimate, ora$estimate, tolerance = 1e-12)
      expect_equal(est_info$estimate, pl_ref$estimate, tolerance = 1e-12)
      expect_equal(est_info$treated_mean, ora$treated_mean, tolerance = 1e-12)
      expect_equal(est_info$control_mean, ora$control_mean, tolerance = 1e-12)
      expect_equal(est_info$dp, ora$mean_a, tolerance = 1e-12)
      expect_equal(est_info$dm, ora$mean_b, tolerance = 1e-12)
      expect_equal(as.double(pl_new$coefficients), ora$beta, tolerance = 1e-12)
      expect_identical(names(pl_new$coefficients), colnames(f_owner$Z_a))
      expect_equal(unname(pl_new$partials$G_beta), unname(ora$G_beta), tolerance = 1e-12)
      expect_equal(unname(pl_new$partials$G_r), ora$G_r, tolerance = 1e-12)
      expect_equal(unname(pl_new$loadings$loading_beta), ora$loading_beta, tolerance = 1e-12)
      expect_equal(unname(pl_new$loadings$loading_rate), ora$loading_rate, tolerance = 1e-12)

      expect_equal(as.double(ch$components$treated[, m]), ora$components$treated, tolerance = 1e-12)
      expect_equal(as.double(ch$components$control[, m]), ora$components$control, tolerance = 1e-12)
      expect_equal(as.double(ch$components$beta_score[, m]), ora$components$beta_score, tolerance = 1e-12)
      expect_equal(as.double(ch$components$rate[, m]), ora$components$rate, tolerance = 1e-12)
      expect_equal(as.double(ch$phi[, m]), ora$phi, tolerance = 1e-12)
      expect_equal(as.double(ch$phi[, m]), as.double(if_ref$phi), tolerance = 1e-12)
    }

    cov_new <- cp04_stream_covariance(b, seq_len(M))
    cov_ref <- stream_ipw_point_covariance_view(ref_links, ctrl$point_control)
    expect_equal(as.double(cov_new$estimates), as.double(cov_ref$estimates), tolerance = 1e-12)
    expect_equal(unname(cov_new$covariance), unname(cov_ref$covariance), tolerance = 1e-12)
    expect_equal(unname(cov_new$covariance), unname(crossprod(sweep(ch$phi, 2L, colMeans(ch$phi), "-")) / f_owner$n_a_units), tolerance = 1e-12)
    expect_equal(as.double(cov_new$mean_phi), as.double(cov_ref$mean_phi), tolerance = 1e-12)

    for (cg in c("notyettreated", "nevertreated")) {
      for (pt_mode in c("post", "all")) {
        for (bp in c("varying", "universal")) {
          res_new <- run_chained_ipw_pipeline(
            data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
            control = ctrl, xformla = ~x1, weightsname = "w",
            control_group = cg, base_period = bp, pt_assumption = pt_mode,
            pt_start = if (pt_mode == "all") p_owner$all_times[1L] else NULL,
            weighting = "optimal", reduce_structural_redundancy = TRUE,
            run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
          )
          res_ref <- run_legacy_scalar_ipw_pipeline(
            data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
            control = ctrl, xformla = ~x1, weightsname = "w",
            control_group = cg, base_period = bp, pt_assumption = pt_mode,
            pt_start = if (pt_mode == "all") p_owner$all_times[1L] else NULL,
            weighting = "optimal", reduce_structural_redundancy = TRUE,
            run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
          )
          expect_s3_class(res_new$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
          expect_null(res_new$nuisance_registry$cp04)
          expect_identical(res_new$metadata$catalog$rows, res_ref$metadata$catalog$rows)
          expect_identical(res_new$reporting$target_meta, res_ref$reporting$target_meta)
          expect_identical(res_new$auxiliary_reporting$descriptor_status, res_ref$auxiliary_reporting$descriptor_status)
          expect_identical(names(res_new$nuisance_registry$link_views), names(res_ref$nuisance_registry$link_views))
          expect_identical(vapply(res_new$nuisance_registry$nuisance_views, `[[`, character(1L), "nuisance_id"),
            vapply(res_ref$nuisance_registry$nuisance_views, `[[`, character(1L), "nuisance_id"))
          for (j in seq_along(res_new$nuisance_registry$raw_aux_points)) {
            pn <- res_new$nuisance_registry$raw_aux_points[[j]]
            pr <- res_ref$auxiliary_reporting$point_links[[j]]
            if (is.null(pn)) {
              expect_null(pr)
            } else {
              expect_identical(pn$link_id, pr$link_id)
              expect_identical(pn$functional_key, pr$functional_key)
              expect_identical(pn$nuisance_id, pr$nuisance_id)
              expect_identical(res_new$nuisance_registry$raw_ipw$links[[pn$cp04_index]]$link_id, pn$link_id)
            }
          }
          expect_equal(res_new$reporting$estimates, res_ref$reporting$estimates, tolerance = 1e-12)
          expect_equal(res_new$reporting$se, res_ref$reporting$se, tolerance = 1e-12)
          expect_equal(res_new$reporting$vcov, res_ref$reporting$vcov, tolerance = 1e-12)
          expect_equal(res_new$reporting$Omega_a, res_ref$reporting$Omega_a, tolerance = 1e-12)
          expect_equal(res_new$reporting$if_full, res_ref$reporting$if_full, tolerance = 1e-12)
          expect_equal(res_new$covariance_stream$mhat, res_ref$covariance_stream$mhat, tolerance = 1e-12)
          expect_equal(res_new$covariance_stream$V, res_ref$covariance_stream$V, tolerance = 1e-12)
          expect_identical(res_new$covariance_handling$handler_used, res_ref$covariance_handling$handler_used)
          expect_identical(res_new$covariance_handling$handler_invoked, res_ref$covariance_handling$handler_invoked)
          expect_equal(res_new$covariance_handling$lambda, res_ref$covariance_handling$lambda, tolerance = 1e-12)
          for (agg_nm in names(res_ref$aggregations)) {
            expect_equal(res_new$aggregations[[agg_nm]]$estimates, res_ref$aggregations[[agg_nm]]$estimates, tolerance = 1e-12)
            expect_equal(res_new$aggregations[[agg_nm]]$se, res_ref$aggregations[[agg_nm]]$se, tolerance = 1e-12)
            expect_equal(res_new$aggregations[[agg_nm]]$if_full, res_ref$aggregations[[agg_nm]]$if_full, tolerance = 1e-12)
          }
        }
      }
    }
  }
})

test_that("D93 compact raw-IPW backend supports F > A with unknown off-A covariates and cohorts and exact zero transport", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 48L, missing_mode = "partial", unequal_weights = TRUE)
  ctrl <- make_d93_test_pipeline_control()

  n_extra <- 15L
  extra_id <- seq.int(1001L, length.out = n_extra)
  full_id <- c(sim$unit_id, extra_id)
  full_Z <- rbind(
    cbind(`(Intercept)` = 1.0, x1 = sim$unit_x1),
     matrix(NA_real_, nrow = n_extra, ncol = 2L, dimnames = list(NULL, c("(Intercept)", "x1")))
  )
  full_w <- c(sim$unit_w, rep(1.25, n_extra))
  full_g <- c(sim$unit_g, rep(NA_real_, n_extra))
  off_a_rows <- seq.int(length(sim$unit_id) + 1L, length(full_id))

  # The existing A-only constructor accepts unknown off-A X. Omit frame_g:
  # panel cohorts are known, and the assembler records unknown external cohorts
  # as NA. An explicitly supplied frame_g still requires known values.
  po <- construct_ipw_panel_owner(sim$data, "id", "time", "y", "g", ctrl$owner_control)
  fo <- construct_ipw_frame_owner(po, full_id, full_Z, full_w, ctrl$owner_control,
    frame_g = NULL, allow_frame_only = TRUE, design_on_A_only = TRUE)
  expect_identical(fo$frame_only_units$g, full_g[off_a_rows])
  expect_identical(fo$Z_a, full_Z[match(po$a_id, full_id), , drop = FALSE])
  meta <- build_owner_comparison_catalog(po, targets = c(2, 3), anticipation = 0L,
    control_group = "notyettreated", conventions = "D10_D13_admissible_frame",
    control = ctrl$owner_control, catalog_bounds = ctrl$catalog_bounds,
    frame_owner = fo, build_design = TRUE, pt_assumption = "post", pt_start = NULL)
  reg <- build_ipw_raw_backend_registry(meta, fo, owner_control = ctrl$owner_control,
    fitted_control = ctrl$fitted_control, point_control = ctrl$point_control,
    glm_control = ctrl$glm_control, class_bounds = ctrl$class_bounds,
    weighting = "identity", reduce_structural_redundancy = FALSE)
  b <- reg$raw_ipw
  ids <- seq_along(reg$link_views)
  raw <- evaluate_ipw_raw_if_chunk(b, seq_len(fo$n_a_units), ids, details = TRUE)
  scalar <- matrix(0.0, nrow = fo$n_a_units, ncol = length(ids),
    dimnames = dimnames(raw$phi))
  scalar_points <- vector("list", length(ids))
  for (m in ids) {
    sp <- reg$row_specs[[m]]
    ora <- scalar_oracle_ipw_link(fo, po, sp$g, sp$K, sp$s, sp$t, ctrl$glm_control)
    views <- construct_pipeline_catalog_row_views(meta, m, ctrl$owner_control, fo)
    fv <- fit_ipw_nuisance_view(views$nuisance_view, ctrl$glm_control, ctrl$fitted_control)
    pl <- fit_ipw_point_link_view(fv, views$link_view, ctrl$point_control)
    scalar_points[[m]] <- pl
    scalar[, m] <- evaluate_ipw_point_if_view_chunk(pl, 1L, fo$n_a_units,
      ctrl$point_control, include_components = TRUE)$phi
    expect_equal(ipw_raw_link_estimate(b, m), pl$estimate, tolerance = 1e-12)
    expect_equal(ipw_raw_link_estimate(b, m), ora$estimate, tolerance = 1e-12)
    expect_equal(as.double(raw$phi[, m]), ora$phi, tolerance = 1e-12)
  }
  expect_equal(raw$phi, scalar, tolerance = 1e-12)
  expect_equal(unname(cp04_stream_covariance(b, ids)$covariance),
    unname(stream_ipw_point_covariance_view(scalar_points, ctrl$point_control)$covariance), tolerance = 1e-12)
  a_to_full <- match(po$a_id, full_id)
  raw_full <- transport_if_to_full_frame(raw$phi, length(full_id), a_to_full)
  scalar_full <- transport_if_to_full_frame(scalar, length(full_id), a_to_full)
  expect_equal(raw_full, scalar_full, tolerance = 1e-12)
  expect_identical(raw_full$scale_factor, length(full_id) / fo$n_a_units)
  expect_true(all(raw_full$Phi_full[off_a_rows, , drop = FALSE] == 0.0))
  unknown_domain <- cp04_backend_domain(b)
  expect_identical(unknown_domain$N, as.integer(fo$n_a_units))
  expect_identical(unknown_domain$n_A, unknown_domain$N)
  expect_identical(unknown_domain$lambda, 1)
  expect_null(unknown_domain$reader)

  # Preserve the private scalar pipeline's existing full-Z guard. Its valid
  # finite-carrier route is compared below; unknown off-A X was tested above
  # through the existing constructor, backend and original transport.
  expect_error(run_chained_ipw_pipeline(
    data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, xformla = ~x1,
    frame_id = full_id, frame_Z = full_Z, frame_weights = full_w, frame_g = NULL,
    allow_frame_only = TRUE, run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE),
    class = "cdid_ipw_owner_covariate_error")
  finite_full_Z <- full_Z
  finite_full_Z[off_a_rows, ] <- cbind(1.0, rep(0.0, n_extra))

  res_new <- run_chained_ipw_pipeline(
    data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, xformla = ~x1,
    frame_id = full_id, frame_Z = finite_full_Z, frame_weights = full_w, frame_g = NULL,
    allow_frame_only = TRUE,
    control_group = "notyettreated", base_period = "varying", pt_assumption = "post",
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
  )
  res_ref <- run_legacy_scalar_ipw_pipeline(
    data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, xformla = ~x1,
    frame_id = full_id, frame_Z = finite_full_Z, frame_weights = full_w, frame_g = NULL,
    allow_frame_only = TRUE,
    control_group = "notyettreated", base_period = "varying", pt_assumption = "post",
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
  )

  expect_s3_class(res_new$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
  expect_null(res_new$nuisance_registry$cp04)
  domain <- cp04_backend_domain(res_new$nuisance_registry$raw_ipw)
  expect_identical(domain$N, as.integer(length(sim$unit_id)))
  expect_identical(domain$n_A, domain$N)
  expect_identical(domain$lambda, 1)
  expect_null(domain$reader)
  expect_equal(nrow(res_new$reporting$if_full), length(full_id))
  expect_identical(res_new$frame_owner$Z_a, fo$Z_a)
  expect_identical(names(res_new$nuisance_registry$link_views), names(reg$link_views))
  expect_identical(res_new$metadata$catalog$rows, meta$catalog$rows)
  expect_true(all(res_new$reporting$if_full[off_a_rows, res_new$reporting$available] == 0.0))
  expect_equal(res_new$reporting$estimates, res_ref$reporting$estimates, tolerance = 1e-12)
  expect_equal(res_new$reporting$se, res_ref$reporting$se, tolerance = 1e-12)
  expect_equal(res_new$reporting$if_full, res_ref$reporting$if_full, tolerance = 1e-12)
})

test_that("D93 compact raw-IPW backend succeeds when OR control-only WLS rank fails", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 60L, cohorts = c(2, Inf), missing_mode = "partial", unequal_weights = TRUE)
  # Make x2 vary with both signs within cohort 2 and equal zero on cohort Inf,
  # so control-only WLS on K = {Inf} is rank-deficient while pooled {g, K}
  # propensity logit design has full column rank.
  u_match <- match(sim$data$id, sim$unit_id)
  x2_unit <- ifelse(sim$unit_g == 2, cos(sim$unit_id * 0.9) + 0.3, 0.0)
  sim$data$x2 <- x2_unit[u_match]
  ctrl <- make_d93_test_pipeline_control()
  public_data <- sim$data
  public_data$g[public_data$g == Inf] <- 0

  expect_error(
    att_gt(
      yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~x1 + x2, data = public_data,
      est_method = "reg", control_group = "nevertreated",
      bstrap = FALSE, cband = FALSE
    ),
    class = "cdid_or_wls_rank_error"
  )
  expect_error(
    att_gt(
      yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~x1 + x2, data = public_data,
      est_method = "dr", control_group = "nevertreated",
      bstrap = FALSE, cband = FALSE
    ),
    class = "cdid_or_wls_rank_error"
  )

  res_ipw <- run_chained_ipw_pipeline(
    data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, xformla = ~x1 + x2, weightsname = "w",
    control_group = "nevertreated", base_period = "varying", pt_assumption = "post",
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
  )
  res_ref <- run_legacy_scalar_ipw_pipeline(
    data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, xformla = ~x1 + x2, weightsname = "w",
    control_group = "nevertreated", base_period = "varying", pt_assumption = "post",
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
  )
  expect_s3_class(res_ipw$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
  expect_equal(res_ipw$reporting$estimates, res_ref$reporting$estimates, tolerance = 1e-12)
  expect_equal(res_ipw$reporting$se, res_ref$reporting$se, tolerance = 1e-12)
  expect_equal(res_ipw$reporting$if_full, res_ref$reporting$if_full, tolerance = 1e-12)
})

test_that("D93 inactive bounds, capped odds and floored rates preserve raw score/H and all four IF terms", {
  ctrl <- make_d93_test_pipeline_control(require_matching_scale = FALSE)
  for (bound_case in c("inactive", "odds", "rates")) {
    sim <- make_d93_compact_ipw_test_panel(
      n_units = 60L, cohorts = c(2, Inf), missing_mode = "complete", unequal_weights = TRUE
    )
    if (bound_case == "odds") {
      # Two strata each contain both outcomes. Their weighted treatment/control
      # ratios are 1 and exp(14), giving the finite population logit (0, 1)
      # on design (1, x1), x1 in {0, 14}. No coordinate separates the outcomes.
      sim$unit_x1 <- rep(c(0, 0, 14, 14), length.out = length(sim$unit_id))
      sim$unit_w <- ifelse(sim$unit_g == Inf & sim$unit_x1 == 14, exp(-14), 1)
    }
    if (bound_case == "rates") {
      # Match X and weights within target/control pairs, including the heavy
      # pair absent at date 4. Every raw logit then has the finite interior
      # root beta = 0, while date-4 pair rates in both cohorts are below 1e-6.
      # All six date pairs keep positive support and a full-rank design.
      paired <- seq.int(1L, length(sim$unit_id), by = 2L)
      sim$unit_x1 <- rep(sim$unit_x1[paired], each = 2L)
      sim$unit_w <- rep(sim$unit_w[paired], each = 2L)
      sim$unit_w[1:2] <- 1e9
      sim$data <- sim$data[!(sim$data$id %in% 1:2 & sim$data$time == 4), , drop = FALSE]
    }
    ix <- match(sim$data$id, sim$unit_id)
    sim$data$x1 <- sim$unit_x1[ix]
    sim$data$w <- sim$unit_w[ix]
    po <- construct_ipw_panel_owner(sim$data, "id", "time", "y", "g", ctrl$owner_control)
    fo <- construct_ipw_frame_owner(po, sim$unit_id,
      cbind(`(Intercept)` = 1, x1 = sim$unit_x1), sim$unit_w,
      ctrl$owner_control, frame_g = sim$unit_g)
    meta <- build_owner_comparison_catalog(po, targets = 2, anticipation = 0L,
      control_group = "nevertreated", conventions = "D10_D13_admissible_frame",
      control = ctrl$owner_control, catalog_bounds = ctrl$catalog_bounds,
      frame_owner = fo, build_design = TRUE, pt_assumption = "all", pt_start = po$all_times[1L])
    reg <- build_ipw_raw_backend_registry(meta, fo, owner_control = ctrl$owner_control,
      fitted_control = ctrl$fitted_control, point_control = ctrl$point_control,
      glm_control = ctrl$glm_control, class_bounds = ctrl$class_bounds,
      weighting = "identity", reduce_structural_redundancy = FALSE)
    b <- reg$raw_ipw
    ids <- seq_along(b$links)
    ch <- evaluate_ipw_raw_if_chunk(b, seq_len(fo$n_a_units), ids, details = TRUE)
    summaries <- lapply(reg$fitted_views, function(fv) fv$diagnostics$clipping)
    if (bound_case == "odds") expect_true(any(vapply(summaries, function(z) z[["p_capped_ctrl"]] > 0, logical(1))))
    if (bound_case == "rates") expect_true(any(vapply(summaries, function(z) z[["r_floored"]] > 0, logical(1))))
    if (bound_case == "inactive") expect_true(all(vapply(summaries, function(z) z[["p_capped"]] == 0 && z[["r_floored"]] == 0, logical(1))))
    for (m in ids) {
      sp <- reg$row_specs[[m]]
      ora <- scalar_oracle_ipw_link(fo, po, sp$g, sp$K, sp$s, sp$t, ctrl$glm_control)
      nv <- construct_ipw_nuisance_view(fo, sp$g, sp$K, sp$s, sp$t, ctrl$owner_control)
      fv <- fit_ipw_nuisance_view(nv, ctrl$glm_control, ctrl$fitted_control)
      lv <- construct_ipw_link_view(fo, sp$s, sp$t, sp$g, sp$K, nv, ctrl$owner_control)
      pl <- fit_ipw_point_link_view(fv, lv, ctrl$point_control)
      ref <- evaluate_ipw_point_if_view_chunk(pl, 1L, fo$n_a_units, ctrl$point_control, include_components = TRUE)
      fit <- b$treatments[[b$links[[m]]$model]]
      expect_equal(unname(fit$H), unname(ora$H), tolerance = 1e-12)
      expect_equal(unname(fit$A_beta_r), unname(ora$A_beta_r), tolerance = 1e-12)
      expect_equal(unname(fit$partials$G_beta), unname(ora$G_beta), tolerance = 1e-12)
      expect_equal(unname(fit$partials$G_r), unname(ora$G_r), tolerance = 1e-12)
      expect_equal(ipw_raw_link_estimate(b, m), ora$estimate, tolerance = 1e-12)
      for (channel in names(ora$components)) {
        expect_equal(as.double(ch$components[[channel]][, m]), ora$components[[channel]], tolerance = 1e-12)
        expect_equal(as.double(ch$components[[channel]][, m]), as.double(ref$components[[channel]]), tolerance = 1e-12)
      }
      expect_equal(as.double(ch$phi[, m]), ora$phi, tolerance = 1e-12)
      expect_equal(as.double(ch$phi[, m]), as.double(ref$phi), tolerance = 1e-12)
    }
    gram <- crossprod(sweep(ch$phi, 2L, colMeans(ch$phi), "-")) / fo$n_a_units
    expect_equal(unname(cp04_stream_covariance(b, ids)$covariance), unname(gram), tolerance = 1e-12)
  }
})

test_that("D93 compact raw-IPW backend preserves cluster inference, retention decline/streaming, LAST callback integrity, and typed refusals", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 60L, missing_mode = "partial", unequal_weights = TRUE)
  ctrl <- make_d93_test_pipeline_control()

  # 1. Clustered inference parity
  res_cl <- run_chained_ipw_pipeline(
    data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, xformla = ~x1, weightsname = "w", clustername = "cluster",
    control_group = "notyettreated", base_period = "varying", pt_assumption = "post",
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
  )
  res_cl_ref <- run_legacy_scalar_ipw_pipeline(
    data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, xformla = ~x1, weightsname = "w", clustername = "cluster",
    control_group = "notyettreated", base_period = "varying", pt_assumption = "post",
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
  )
  expect_equal(res_cl$reporting$se, res_cl_ref$reporting$se, tolerance = 1e-12)
  expect_equal(res_cl$reporting$vcov, res_cl_ref$reporting$vcov, tolerance = 1e-12)
  expect_equal(res_cl$reporting$inference$cluster_if, res_cl_ref$reporting$inference$cluster_if, tolerance = 1e-12)
  expect_identical(res_cl$resource_accounting$point_producer$retained_moment_if$status, "retained")

  # 2. The original planner admits one block; an already-full explicit live
  # inventory declines optional storage before allocation. Replay stays streamed.
  b <- res_cl$nuisance_registry$raw_ipw
  ids <- vapply(res_cl$point_links_used, `[[`, integer(1L), "cp04_index")
  cov_plan <- cp04_stream_plan(b, ids)
  expect_length(cov_plan$blocks, 1L)
  declined <- cp04_retained_if_reserve(
    b, ids, cov_plan$chunk_rows, cov_plan, retained = 0,
    gmm_retained = 0, live = min(b$control$max_retained_entries, b$control$max_matrix_entries),
    work = 0, control = b$control
  )
  expect_null(declined)
  streamed <- cp04_stream_covariance(b, ids)
  expect_equal(unname(streamed$estimates), unname(res_cl$covariance_stream$mhat), tolerance = 1e-12)
  expect_equal(unname(streamed$covariance), unname(res_cl$covariance_stream$V), tolerance = 1e-12)
  expect_equal(unname(streamed$mean_phi), unname(res_cl$covariance_stream$mean_phi), tolerance = 1e-12)

  # 3. LAST-callback mutation detection on streamed replay. Restore the shared
  # owner binding even when the assertion fails, so later checks are isolated.
  frame_before <- get("frame_owner", envir = b$frame_ref, inherits = FALSE)
  n_clean <- 0L
  cp04_replay_links(b, ids, function(batch) {
    n_clean <<- n_clean + 1L
  }, chunk_rows = 19L)
  expect_equal(n_clean, ceiling(b$n_a_units / 19L))
  expect_identical(frame_before, get("frame_owner", envir = b$frame_ref, inherits = FALSE))
  n_seen <- 0L
  n_mutated <- 0L
  expect_error(
    tryCatch(cp04_replay_links(b, ids, function(batch) {
      n_seen <<- n_seen + 1L
      if (max(batch$rows) == b$n_a_units) {
        # Mutate the frame_owner inside its ref environment on the LAST callback chunk
        n_mutated <<- n_mutated + 1L
        unlockBinding("frame_owner", b$frame_ref)
        fo_bad <- get("frame_owner", envir = b$frame_ref, inherits = FALSE)
        fo_bad$w_a[1L] <- fo_bad$w_a[1L] + 1.0
        assign("frame_owner", fo_bad, envir = b$frame_ref)
        lockBinding("frame_owner", b$frame_ref)
      }
    }, chunk_rows = 19L), finally = {
      unlockBinding("frame_owner", b$frame_ref)
      assign("frame_owner", frame_before, envir = b$frame_ref)
      lockBinding("frame_owner", b$frame_ref)
    }),
    class = "cdid_ipw_owner_integrity_error"
  )

  expect_equal(n_seen, ceiling(b$n_a_units / 19L))
  expect_identical(n_mutated, 1L)
  expect_identical(frame_before, get("frame_owner", envir = b$frame_ref, inherits = FALSE))

  # Invalid seals, active bindings and classed indices are refused before any
  # getter/method can run. These copies never alter the admitted backend.
  bad <- b
  bad$control$max_work_ops <- bad$control$max_work_ops - 1
  expect_error(validate_ipw_raw_backend(bad), class = "cdid_ipw_raw_backend_integrity_error")
  getter_calls <- 0L
  bad_ref <- new.env(parent = emptyenv())
  bad <- b; bad$ref <- bad_ref
  makeActiveBinding("ipw_raw_backend", function() { getter_calls <<- getter_calls + 1L; bad }, bad_ref)
  lockEnvironment(bad_ref, bindings = TRUE)
  expect_error(validate_ipw_raw_backend(bad), class = "cdid_ipw_raw_backend_integrity_error")
  expect_identical(getter_calls, 0L)
  expect_error(ipw_raw_if_budget(b, structure(1, class = "d93_foreign"), ids),
    class = "cdid_ipw_raw_backend_input_error")

  # 4. Original typed refusals (bound error on max_observer_units)
  ctrl_small_obs <- make_d93_test_pipeline_control(max_observer_units = 5L)
  expect_error(
    run_chained_ipw_pipeline(
      data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
      control = ctrl_small_obs, xformla = ~x1, weightsname = "w",
      run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE
    ),
    class = "cdid_ipw_owner_bound_error"
  )
})

# V15: existing D93 complete/scalar tests above retain all formula checks.
# These witnesses check optional reuse admission and producer identity only.
test_that("complete multioperator reuse funds all consumers and declines before fitting under pressure", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 60L, times = 1:4,
    cohorts = c(2, 3, Inf), missing_mode = "complete", unequal_weights = TRUE)
  args <- list(data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = make_d93_test_pipeline_control(), xformla = ~x1, weightsname = "w",
    control_group = "nevertreated", base_period = "varying", pt_assumption = "post",
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE)
  actual <- do.call(run_chained_ipw_pipeline, args)
  legacy <- do.call(run_legacy_scalar_ipw_pipeline, args)
  witness <- actual$resource_accounting$complete_raw_capacity
  expect_identical(witness$admitted, TRUE)
  expect_s3_class(actual$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
  expect_null(legacy$nuisance_registry$raw_ipw)
  expect_lte(actual$resource_accounting$simultaneous_entries, witness$peak_entries)
  expect_lte(actual$nuisance_registry$raw_ipw$accounting$construction_peak_entries, witness$construction_entries)
  # Stores have their own cumulative optional work admission. The route
  # witness bounds mandatory work; independently admitted store/check work
  # remains in the returned total and under the same configured cap.
  store <- actual$resource_accounting$point_producer$retained_moment_if
  optional_work <- sum(c(store$work_ops_bound, store$auxiliary$work_ops_bound))
  expect_gte(optional_work, 0)
  expect_lte(optional_work, actual$resource_accounting$downstream_work_ops)
  expect_lte(actual$resource_accounting$downstream_work_ops - optional_work, witness$mandatory_work_upper)
  expect_lte(actual$resource_accounting$downstream_work_ops, args$control$owner_control$max_work_ops)
  expect_equal(actual$reporting$if_full, legacy$reporting$if_full, tolerance = 1e-12)
  # The old route itself establishes this finite capacity witness. Its actual
  # recorded initial/consumer envelope plus fixed bookkeeping remains allowed;
  # optional raw table/consumer obligations must decline before any fit.
  cap <- ceiling(max(legacy$resource_accounting$initial_envelope_entries,
    legacy$resource_accounting$simultaneous_entries) + 4096)
  args$control <- make_d93_test_pipeline_control(max_retained_entries = cap,
    max_matrix_entries = cap, max_output_entries = cap)
  old_pressure <- do.call(run_legacy_scalar_ipw_pipeline, args)
  pressure <- do.call(run_chained_ipw_pipeline, args)
  expect_identical(pressure$resource_accounting$complete_raw_capacity$admitted, FALSE)
  expect_null(pressure$nuisance_registry$raw_ipw)
  expect_equal(pressure$reporting$estimates, old_pressure$reporting$estimates, tolerance = 1e-12)
  expect_equal(pressure$reporting$if_full, old_pressure$reporting$if_full, tolerance = 1e-12)
  expect_equal(pressure$reporting$vcov, old_pressure$reporting$vcov, tolerance = 1e-12)
  expect_identical(pressure$support$ledger, old_pressure$support$ledger)
  expect_identical(pressure$covariance_handling, old_pressure$covariance_handling)
  # A truly unavailable old request keeps its original typed resource refusal.
  args$control <- make_d93_test_pipeline_control(max_retained_entries = 1000,
    max_matrix_entries = 1000, max_output_entries = 1000)
  expect_error(do.call(run_chained_ipw_pipeline, args), class = "cdid_ipw_owner_bound_error")
})

test_that("complete singleoperator routing and raw callback source checks remain intact", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 60L, times = 1:4,
    cohorts = c(2, Inf), missing_mode = "complete", unequal_weights = TRUE)
  single <- run_chained_ipw_pipeline(data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = make_d93_test_pipeline_control(), xformla = ~x1, weightsname = "w",
    control_group = "nevertreated", pt_assumption = "post", weighting = "optimal",
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE)
  expect_null(single$nuisance_registry$raw_ipw)
  expect_null(single$resource_accounting$complete_raw_capacity)
  sim <- make_d93_compact_ipw_test_panel(n_units = 60L, times = 1:4,
    cohorts = c(2, 3, Inf), missing_mode = "complete", unequal_weights = TRUE)
  multi <- run_chained_ipw_pipeline(data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = make_d93_test_pipeline_control(), xformla = ~x1, weightsname = "w",
    control_group = "nevertreated", pt_assumption = "post", weighting = "optimal",
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE)
  b <- multi$nuisance_registry$raw_ipw
  expect_s3_class(b, "cdid_ipw_raw_backend")
  ids <- vapply(multi$point_links_used, `[[`, integer(1), "cp04_index")
  before <- get("frame_owner", b$frame_ref, inherits = FALSE)
  expect_error(tryCatch(cp04_replay_links(b, ids, function(batch) {
    if (max(batch$rows) == b$n_a_units) {
      unlockBinding("frame_owner", b$frame_ref)
      bad <- get("frame_owner", b$frame_ref, inherits = FALSE); bad$w_a[1L] <- bad$w_a[1L] + 1
      assign("frame_owner", bad, envir = b$frame_ref); lockBinding("frame_owner", b$frame_ref)
    }
  }, chunk_rows = 19L), finally = {
    unlockBinding("frame_owner", b$frame_ref); assign("frame_owner", before, envir = b$frame_ref)
    lockBinding("frame_owner", b$frame_ref)
  }), class = "cdid_ipw_owner_integrity_error")
  expect_identical(before, get("frame_owner", b$frame_ref, inherits = FALSE))
})

# V16: exact per-invocation complete support facts, without coalescing pairs.
test_that("complete pair summary reuse preserves raw chunk-order masses, fingerprints and original pair identities", {
  for (miss in c("complete", "partial", "mixed")) {
    sim <- make_d93_compact_ipw_test_panel(n_units = 60L, times = 1:4,
      cohorts = c(2, 3, Inf), missing_mode = if (miss == "mixed") "complete" else miss, unequal_weights = TRUE)
    if (miss == "mixed") sim$data <- sim$data[!(sim$data$time == 4 & sim$data$id %% 7L == 0L), , drop = FALSE]
    ctrl <- make_d93_test_pipeline_control(max_chunk_rows = 7L)
    po <- construct_ipw_panel_owner(sim$data, "id", "time", "y", "g", ctrl$owner_control)
    fo <- construct_ipw_frame_owner(po, sim$unit_id,
      cbind(`(Intercept)` = 1, x1 = sim$unit_x1), sim$unit_w, ctrl$owner_control, frame_g = sim$unit_g)
    meta <- build_owner_comparison_catalog(po, targets = c(2, 3), anticipation = 0L,
      control_group = "notyettreated", conventions = "D10_D13_admissible_frame", control = ctrl$owner_control,
      catalog_bounds = ctrl$catalog_bounds, frame_owner = fo, build_design = TRUE,
      pt_assumption = "all", pt_start = po$all_times[1L])
    args <- list(metadata_obj = meta, frame_owner = fo, design = meta$design,
      owner_control = ctrl$owner_control, fitted_control = ctrl$fitted_control, point_control = ctrl$point_control,
      glm_control = ctrl$glm_control, class_bounds = ctrl$class_bounds,
      weighting = "identity", reduce_structural_redundancy = FALSE)
    calls <- 0L; summaries <- list()
    cached <- build_ipw_raw_backend_registry
    spy_env <- new.env(parent = environment(cached))
    spy_env$ipw_raw_scan_pair_summary <- function(...) {
      calls <<- calls + 1L
      ipw_raw_scan_pair_summary(...)
    }
    spy_env$ipw_raw_pair_support_for_cohorts <- function(frame_owner, pair_summary, ...) {
      summaries[[as.character(pair_summary$pair_id)]] <<- pair_summary
      ipw_raw_pair_support_for_cohorts(frame_owner, pair_summary, ...)
    }
    environment(cached) <- spy_env
    actual <- do.call(cached, args)
    # Only this builder-local proof is declined. The unchanged standalone
    # scanner still proves/visits every true complete pair independently.
    reference_fn <- build_ipw_raw_backend_registry
    env <- new.env(parent = environment(reference_fn))
    env$ipw_owner_pair_complete <- function(...) FALSE
    environment(reference_fn) <- env
    reference <- do.call(reference_fn, args)
    b <- actual$raw_ipw; old <- reference$raw_ipw
    n_complete <- sum(vapply(b$pairs, `[[`, logical(1), "complete_pair"))
    expect_identical(b$pairs, old$pairs)
    expect_identical(b$support_specs, old$support_specs)
    expect_identical(b$supports, old$supports)
    expect_identical(b$treatment_specs, old$treatment_specs)
    expect_equal(b$treatments, old$treatments, tolerance = 1e-12)
    expect_identical(lapply(actual$nuisance_views, function(v) v$support),
      lapply(reference$nuisance_views, function(v) v$support))
    expect_equal(calls, length(b$pairs) - max(0, n_complete - 1L))
    expect_equal(actual$construction_work$complete_pair_fingerprints_reused, max(0, n_complete - 1L))
    for (pid in seq_along(b$pairs)) {
      raw <- ipw_raw_scan_pair_summary(fo, po, b$pairs[[pid]]$pair, 7L)
      raw$pair_id <- as.integer(pid)
      expect_identical(summaries[[as.character(pid)]], raw)
      expect_identical(raw$complete_pair, b$pairs[[pid]]$complete_pair)
      for (sid in which(vapply(b$support_specs, `[[`, integer(1), "pair_id") == pid)) {
        k <- b$support_specs[[sid]]$cohort
        expect_identical(b$supports[[sid]]$n_pair, raw$cohort_pair_count[k])
        expect_identical(b$supports[[sid]]$pair_mass, raw$cohort_pair_mass[k])
      }
    }
    rows <- 7:1; links <- rev(seq_along(b$links))
    expect_equal(evaluate_ipw_raw_if_chunk(b, rows, links)$phi,
      evaluate_ipw_raw_if_chunk(old, rows, links)$phi, tolerance = 1e-12)
    expect_identical(b$accounting$construction_table_entries_bound, old$accounting$construction_table_entries_bound)
    expect_identical(b$accounting$construction_peak_entries, old$accounting$construction_peak_entries)
    expect_identical(b$accounting$work_ops_bound, old$accounting$work_ops_bound)
  }
})

test_that("raw replay plans extract dispatch metadata once inside the original workspace reservation", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 60L, missing_mode = "partial")
  ctrl <- make_d93_test_pipeline_control(max_chunk_rows = 7L)
  po <- construct_ipw_panel_owner(sim$data, "id", "time", "y", "g", ctrl$owner_control)
  fo <- construct_ipw_frame_owner(po, sim$unit_id,
    cbind(`(Intercept)` = 1, x1 = sim$unit_x1), sim$unit_w, ctrl$owner_control, frame_g = sim$unit_g)
  catalog <- build_owner_comparison_catalog(po, targets = c(2, 3), anticipation = 0L,
    control_group = "notyettreated", conventions = "D10_D13_admissible_frame", control = ctrl$owner_control,
    catalog_bounds = ctrl$catalog_bounds, frame_owner = fo, build_design = TRUE,
    pt_assumption = "all", pt_start = po$all_times[1L])
  registry <- build_ipw_raw_backend_registry(catalog, fo, catalog$design,
    owner_control = ctrl$owner_control, fitted_control = ctrl$fitted_control, point_control = ctrl$point_control,
    glm_control = ctrl$glm_control, class_bounds = ctrl$class_bounds,
    weighting = "identity", reduce_structural_redundancy = FALSE)
  b <- registry$raw_ipw
  links <- rev(seq_along(b$links))[seq.int(1L, length(b$links), by = 2L)]
  calls <- 0L
  spy <- new.env(parent = environment(ipw_raw_if_plan))
  spy$ipw_raw_if_metadata <- function(...) {
    calls <<- calls + 1L
    ipw_raw_if_metadata(...)
  }
  planner <- ipw_raw_if_plan; environment(planner) <- spy
  evaluator <- evaluate_ipw_raw_if_chunk; environment(evaluator) <- spy
  plan <- planner(b, links, 7L)
  expect_identical(calls, 1L)
  expect_false(is.null(plan$replay_index))
  # Independent table invariant: both active pair order and each pair's support
  # order agree with the original backend descriptors, including a link subset.
  expect_identical(plan$replay_index$pair_cohorts, lapply(plan$pids, function(pid) {
    specs <- b$support_specs[plan$sids]
    vapply(Filter(function(s) s$pair_id == pid, specs), `[[`, integer(1L), "cohort")
  }))
  expect_identical(plan$replay_index$treatment_pair_pos,
    vapply(b$treatment_specs[plan$tids], function(s) match(s$pair_id, plan$pids), integer(1L)))
  expect_identical(plan$replay_index$nuisance_ids,
    vapply(b$treatment_specs[plan$tids], `[[`, integer(1L), "nuisance_idx"))
  for (rows in list(7:1, 58:60)) {
    indexed <- evaluator(b, rows, links, plan = plan)
    expect_identical(calls, 1L)
    direct <- evaluate_ipw_raw_if_chunk(b, rows, links, details = TRUE)
    expect_equal(indexed$phi, direct$phi, tolerance = 1e-12)
    expect_identical(indexed$id, direct$id)
    expect_identical(indexed$link_ids, direct$link_ids)
  }
  bad <- plan
  bad$replay_index$pair_cohorts[[1L]][1L] <- bad$replay_index$pair_cohorts[[1L]][1L] + 1L
  expect_error(evaluator(b, 1:7, links, plan = bad), class = "cdid_ipw_raw_backend_integrity_error")
  expect_error(evaluator(b, 1:7, rev(links), plan = plan), class = "cdid_ipw_raw_backend_integrity_error")

  # The exact extra native-counted payload is funded inside the old scratch,
  # including the funded_entries field itself. An external consumer may fill
  # every remaining entry without increasing the original simultaneous peak.
  old <- plan; old$replay_index <- NULL; old$funded_entries <- NULL; old$ref <- NULL
  old$entries <- unname(cp04_payload(old)["entries"])
  expect_identical(plan$entries - plan$funded_entries, old$entries)
  expect_identical(plan$pre, ipw_raw_if_budget(b, 7L, links))
  S <- length(plan$sids); P <- length(plan$pids); T <- length(plan$tids); L <- length(links)
  metadata <- 16 * (S + T + L + b$p + 16)
  expect_lte(plan$funded_entries + 7 * S + 4 * P + 4 * T + 256 + 2 * L + 64, metadata)
  expect_lte(plan$funded_entries + 4 * (S + T + L + b$p) + 128, metadata)
  base_retained <- sum(c(b$accounting$owned_entries, b$accounting$borrowed_entries, old$entries))
  old_peak <- base_retained + old$pre$output + old$pre$scratch
  external <- b$control$max_matrix_entries - old_peak
  expect_equal(external + b$accounting$owned_entries + b$accounting$borrowed_entries +
    plan$entries - plan$funded_entries + plan$pre$output + plan$pre$scratch,
    b$control$max_matrix_entries)
  # Exact matrix/work caps still admit the index because it spends the old
  # reservation; neither cap grows and neither row partition changes.
  for (cap in c("max_matrix_entries", "max_work_ops")) {
    pressure <- b
    pressure$control[[cap]] <- if (cap == "max_matrix_entries") old_peak else old$pre$work
    pressure <- ipw_owner_seal(pressure, "ipw_raw_backend")
    before <- calls
    indexed <- planner(pressure, links, 7L)
    expect_false(is.null(indexed$replay_index))
    expect_identical(calls, before + 1L)
    actual <- evaluator(pressure, 7:1, links, plan = indexed)
    expect_identical(calls, before + 1L)
    expect_equal(actual$phi, evaluate_ipw_raw_if_chunk(b, 7:1, links)$phi, tolerance = 1e-12)
  }
  # The extra retention cap can independently decline, before extracting any
  # metadata. The original plan shape and one-extraction chunk remain valid.
  pressure <- b; pressure$control$max_retained_entries <- base_retained
  pressure <- ipw_owner_seal(pressure, "ipw_raw_backend")
  before <- calls
  plain <- planner(pressure, links, 7L)
  expect_null(plain$replay_index); expect_null(plain$funded_entries)
  expect_identical(calls, before)
  actual <- evaluator(pressure, 7:1, links, plan = plain)
  expect_identical(calls, before + 1L)
  expect_equal(actual$phi, evaluate_ipw_raw_if_chunk(b, 7:1, links)$phi, tolerance = 1e-12)
  # One small request cannot fund this conservative construction inventory.
  # Declining the map preserves the original preflight and direct behavior.
  before <- calls
  small <- planner(b, links[1L], 7L)
  expect_null(small$replay_index); expect_null(small$funded_entries)
  expect_identical(calls, before)
  expect_equal(evaluator(b, 58:60, links[1L], plan = small)$phi,
    evaluate_ipw_raw_if_chunk(b, 58:60, links[1L])$phi, tolerance = 1e-12)
  expect_identical(calls, before + 1L)
})

# V25 frozen collector from 7003ddd, verbatim except its test-local name.
d93_v25_old_observer_collect <- function(frame_owner, panel_owner, pair_summary, coh_desc, n_obs, step) {
  n_a <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  a_cohort_idx <- frame_owner$a_cohort_index
  w_a <- frame_owner$w_a
  Z_a <- frame_owner$Z_a
  pair_desc <- pair_summary$pair
  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  complete_pair <- isTRUE(pair_summary$complete_pair)

  risk_a_idx <- integer(n_obs)
  risk_gi <- integer(n_obs)
  risk_w <- numeric(n_obs)
  Z_use <- matrix(0.0, nrow = n_obs, ncol = p, dimnames = list(NULL, colnames(Z_a)))
  obs_s <- integer(n_obs)
  obs_t <- integer(n_obs)

  pos <- 0L
  a_start <- 1L
  while (a_start <= n_a) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    hit <- if (complete_pair) {
      ipw_owner_query_complete_pair_mask(panel_owner, ds, dt, a_start, a_end)
    } else {
      ipw_owner_query_pair_mask(
        panel_owner, pair_desc$s, pair_desc$t,
        a_start = a_start, a_end = a_end, ds = ds, dt = dt
      )
    }
    a <- hit$a_index
    if (length(a) > 0L) {
      gi <- ipw_owner_view_cohort_index(coh_desc, a_cohort_idx[a])
      in_r <- which(!is.na(gi))
      k <- length(in_r)
      if (k > 0L) {
        if (pos + k > n_obs) {
          ipw_raw_backend_abort("Observed risk support count changed from its pair summary.", "integrity")
        }
        ii <- seq.int(pos + 1L, pos + k)
        a_sub <- a[in_r]
        risk_a_idx[ii] <- a_sub
        risk_gi[ii] <- as.integer(gi[in_r])
        risk_w[ii] <- w_a[a_sub]
        Z_use[ii, ] <- Z_a[a_sub, , drop = FALSE]
        obs_s[ii] <- hit$obs_s[in_r]
        obs_t[ii] <- hit$obs_t[in_r]
        pos <- pos + k
      }
    }
    a_start <- as.double(a_end) + 1
  }
  if (pos != n_obs) {
    ipw_raw_backend_abort("Observed risk support count disagreed with its pair summary.", "integrity")
  }
  list(
    risk_a_idx = risk_a_idx,
    risk_gi = risk_gi,
    risk_w = risk_w,
    Z_use = Z_use,
    obs_s = obs_s,
    obs_t = obs_t
  )
}


# Independent raw-record join. Only the canonical observation-position table
# supplies endpoint row numbers; no posting query or collector supplies truth.
d93_v25_observer_reference <- function(sim, frame, panel, pair, cohorts) {
  observed <- sim$data[!is.na(sim$data$y), c("id", "time"), drop = FALSE]
  ids <- intersect(observed$id[observed$time == pair$s], observed$id[observed$time == pair$t])
  unit_g <- sim$unit_g[match(panel$a_id, sim$unit_id)]
  declared <- match(unit_g, cohorts$cohorts)
  a <- which(panel$a_id %in% ids & !is.na(declared))
  o <- panel$observations
  key <- paste(o$unit_index, o$time)
  ms <- if (length(a)) match(paste(panel$a_unit_index[a], pair$s), key) else integer(0L)
  mt <- if (length(a)) match(paste(panel$a_unit_index[a], pair$t), key) else integer(0L)
  list(risk_a_idx = a, risk_gi = declared[a], risk_w = frame$w_a[a],
    Z_use = frame$Z_a[a, , drop = FALSE], obs_s = ms, obs_t = mt)
}

test_that("coarse raw observers exactly copy the frozen collector and independent ordered risk-set join", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 64L, cohorts = c(2, 3, 4, Inf),
    missing_mode = "partial", unequal_weights = TRUE)
  sim$data$y[sim$data$time == 3 & sim$data$id %% 11L == 0L] <- NA_real_
  empty <- sim$data[1L, , drop = FALSE]; empty$time <- 5; empty$y <- NA_real_
  sim$data <- rbind(sim$data, empty)
  ctrl <- make_d93_test_pipeline_control(max_chunk_rows = 7L)
  fitted <- ctrl$fitted_control; fitted$max_chunk_rows <- 3L
  fitted <- validate_ipw_fitted_view_control(fitted)
  po <- construct_ipw_panel_owner(sim$data, "id", "time", "y", "g", ctrl$owner_control)
  calls <- 0L; reads <- 0L; captured <- NULL
  spy <- new.env(parent = environment(ipw_raw_collect_observer_buffer))
  address <- ipw_raw_observer_collector_kernel()
  spy$.Call <- function(.NAME, ...) {
    if (identical(.NAME, address)) calls <<- calls + 1L
    out <- base::.Call(.NAME, ...)
    if (identical(.NAME, address)) captured <<- out
    out
  }
  spy$ipw_owner_read_y_positions <- function(...) { reads <<- reads + 1L; stop("Collector read an outcome.") }
  collect <- ipw_raw_collect_observer_buffer; environment(collect) <- spy
  designs <- list(matrix(1.0, 64L, 1L, dimnames = list(NULL, "(Intercept)")),
    cbind(`(Intercept)` = 1, x1 = sim$unit_x1, x2 = cos(sim$unit_id * 0.23)))
  for (Z in designs) {
    fo <- construct_ipw_frame_owner(po, sim$unit_id, Z, sim$unit_w, ctrl$owner_control, frame_g = sim$unit_g)
    for (K in list(Inf, c(Inf, 4, 3), c(3, Inf, 4))) {
      cohort <- ipw_owner_make_cohort_descriptor(fo, 2, K)
      for (t in c(4, 5)) {
        pair <- ipw_owner_make_pair_descriptor(po, 1, t)
        summary <- ipw_raw_scan_pair_summary(fo, po, pair, 7L)
        expected <- d93_v25_observer_reference(sim, fo, po, pair, cohort)
        n <- length(expected$risk_a_idx)
        inputs <- serialize(list(po$postings, fo$a_cohort_index, fo$w_a, fo$Z_a, cohort), NULL, version = 3)
        for (step in c(1L, 7L)) {
          before <- calls
          actual <- collect(fo, po, summary, cohort, n, step)
          expect_identical(calls, before + 1L)
          expect_identical(actual, d93_v25_old_observer_collect(fo, po, summary, cohort, n, step))
          expect_identical(actual, expected)
          expect_identical(reads, 0L)
          expect_identical(inputs, serialize(list(po$postings, fo$a_cohort_index, fo$w_a, fo$Z_a, cohort), NULL, version = 3))
          expect_identical(captured$status, 0L)
          for (field in names(actual)) expect_true(gmm_sexp_aliases(actual[[field]], captured[[field]]))
          expect_true(gmm_sexp_aliases(colnames(actual$Z_use), colnames(fo$Z_a)))
          expect_false(gmm_sexp_aliases(actual$Z_use, fo$Z_a))
          E <- (ncol(Z) + 5) * n + 512
          borrowed_names <- unname(cp04_payload(colnames(fo$Z_a))["entries"])
          expect_lte(unname(cp04_payload(captured)["entries"]) - borrowed_names, E)
          B_f <- min(fo$n_a_units, fitted$max_chunk_rows)
          O_old <- (ncol(Z) + 6) * n + 32
          Q_old <- 2 * B_f * ncol(Z) + 88 * B_f + cohort$m_cohorts + 512
          expect_lte(E, O_old + Q_old)
          expect_identical(B_f, 3)
        }
      }
    }
  }
  # Complete support keeps its existing R path, without resolving the new DLL.
  full <- make_d93_compact_ipw_test_panel(n_units = 64L, cohorts = c(2, 3, 4, Inf), missing_mode = "complete")
  po_full <- construct_ipw_panel_owner(full$data, "id", "time", "y", "g", ctrl$owner_control)
  fo_full <- construct_ipw_frame_owner(po_full, full$unit_id, designs[[1L]], full$unit_w,
    ctrl$owner_control, frame_g = full$unit_g)
  cohort <- ipw_owner_make_cohort_descriptor(fo_full, 2, Inf)
  pair <- ipw_owner_make_pair_descriptor(po_full, 1, 4)
  summary <- ipw_raw_scan_pair_summary(fo_full, po_full, pair, 7L)
  n <- length(d93_v25_observer_reference(full, fo_full, po_full, pair, cohort)$risk_a_idx)
  before <- calls
  expect_identical(collect(fo_full, po_full, summary, cohort, n, 7L),
    d93_v25_old_observer_collect(fo_full, po_full, summary, cohort, n, 7L))
  expect_identical(calls, before)
})

test_that("coarse raw observers decline unsupported indices and preserve count failures and native bounds", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 64L, cohorts = c(2, 3, 4, Inf), missing_mode = "partial")
  ctrl <- make_d93_test_pipeline_control(max_chunk_rows = 7L)
  po <- construct_ipw_panel_owner(sim$data, "id", "time", "y", "g", ctrl$owner_control)
  fo <- construct_ipw_frame_owner(po, sim$unit_id, cbind(`(Intercept)` = 1, x1 = sim$unit_x1),
    sim$unit_w, ctrl$owner_control, frame_g = sim$unit_g)
  cohort <- ipw_owner_make_cohort_descriptor(fo, 2, c(Inf, 4, 3))
  pair <- ipw_owner_make_pair_descriptor(po, 1, 4)
  summary <- ipw_raw_scan_pair_summary(fo, po, pair, 7L)
  expected <- d93_v25_observer_reference(sim, fo, po, pair, cohort)
  n <- length(expected$risk_a_idx)
  for (bad_n in c(n - 1L, n + 1L)) {
    old <- tryCatch(d93_v25_old_observer_collect(fo, po, summary, cohort, bad_n, 7L), error = identity)
    actual <- tryCatch(ipw_raw_collect_observer_buffer(fo, po, summary, cohort, bad_n, 7L), error = identity)
    expect_s3_class(actual, "cdid_ipw_raw_backend_integrity_error")
    expect_identical(class(actual), class(old)); expect_identical(conditionMessage(actual), conditionMessage(old))
  }
  spy <- new.env(parent = environment(ipw_raw_collect_observer_buffer))
  spy$ipw_raw_observer_collector_kernel <- function() stop("Ineligible collector resolved native code.")
  fallback <- ipw_raw_collect_observer_buffer; environment(fallback) <- spy
  for (field in c("post_obs_index", "post_a_index")) {
    long <- po; long$postings[[field]] <- as.double(long$postings[[field]])
    expect_identical(fallback(fo, long, summary, cohort, n, 7L),
      d93_v25_old_observer_collect(fo, long, summary, cohort, n, 7L))
  }
  address <- ipw_raw_observer_collector_kernel(); post <- po$postings
  bounds <- c(post$date_start[pair$s_calendar_index], post$date_end[pair$s_calendar_index],
    post$date_start[pair$t_calendar_index], post$date_end[pair$t_calendar_index])
  E <- (fo$n_covariates + 5) * n + 512
  invoke <- function(bd = bounds, owner = cohort$owner_cohort_index, view = cohort$view_cohort_index,
                     observations = post$post_obs_index, ctl = as.integer(c(n, 1L, cohort$m_cohorts)), cap = E) {
    .Call(address, post$post_a_index, observations, bd, fo$a_cohort_index, fo$w_a, fo$Z_a,
      owner, view, ctl, cap)
  }
  expect_identical(invoke()$status, 0L)
  expect_error(invoke(bd = replace(bounds, 1L, 0L)), class = "cdid_gmm_bound_error")
  expect_error(invoke(owner = rev(cohort$owner_cohort_index)), class = "cdid_gmm_bound_error")
  expect_error(invoke(view = replace(cohort$view_cohort_index, 1L, 0L)), class = "cdid_gmm_bound_error")
  expect_error(invoke(cap = E - 1), class = "cdid_gmm_bound_error")
  empty <- invoke(bd = c(1L, 0L, 1L, 0L), ctl = as.integer(c(0, 1, cohort$m_cohorts)), cap = 512)
  expect_identical(empty$status, 0L); expect_length(empty$risk_a_idx, 0L)
  expect_identical(dim(empty$Z_use), c(0L, as.integer(fo$n_covariates)))
  # Deliberate valid-shape mutations must fail the independent buffer truth.
  wrong_view <- cohort$view_cohort_index
  wrong_view[1:2] <- wrong_view[2:1]
  expect_false(identical(invoke(view = wrong_view)$risk_gi, expected$risk_gi))
  expect_false(identical(invoke(observations = post$post_obs_index + 1L)$obs_s, expected$obs_s))
  # Live registered arity and the exact symbol used by the R collector.
  dll <- getLoadedDLLs()[["cdid"]]
  reg <- getDLLRegisteredRoutines(dll, addNames = TRUE)[[".Call"]]
  expect_identical(reg[["C_gmm_ipw_raw_observer_collect"]]$numParameters, 10L)
  expect_identical(address, getNativeSymbolInfo("C_gmm_ipw_raw_observer_collect", PACKAGE = dll)$address)
})

# V29: independent byte/cardinality witnesses for local raw record metadata.
test_that("raw character inventory charges local Unicode/auxiliary roles and scales linearly in D", {
  ctrl <- make_d93_test_pipeline_control()
  long_key <- paste0("catalog:", paste(rep("\u03bb\u6f22", 700L), collapse = ""))
  own_id <- paste0("own:", paste(rep("\u00e9", 180L), collapse = ""))
  keys <- c(rep("same-key", 36L), long_key)
  meta <- list(catalog = list(rows = data.frame(operator_id = rep(1:6, length.out = length(keys)))),
    row_mapping = list(functional_key = keys, link_id = own_id, nuisance_id = "borrowed:nuisance"))
  Z <- matrix(0, 1L, 2L, dimnames = list(NULL, c("(Intercept)", paste0("x_", own_id))))
  frame <- list(Z_a = Z, a_cohorts = c(2, 3, 4, Inf))
  aux <- list(descriptors = list(
    list(g = 2, K = c(3, 4, Inf), link_id = paste0("aux:", own_id)),
    list(g = 3, K = c(4, Inf), nuisance_id = paste0("aux-only:", own_id), functional_key = "custom:aux")))
  inv <- ipw_raw_record_char_inventory(meta, frame, aux, ctrl$fitted_control)
  bytes <- function(x) sum(nchar(x, type = "bytes"))
  common <- bytes(colnames(Z)) + bytes(c("g", "k_1", "k_2", "k_3")) + bytes(as.character(frame$a_cohorts))
  scalar_ids <- c(own_id, "borrowed:nuisance", vapply(aux$descriptors, function(d) {
    x <- c(d$link_id, d$nuisance_id, d$functional_key); if (length(x)) x[which.max(nchar(x, type = "bytes"))] else ""
  }, character(1)))
  id_max <- max(128, nchar(scalar_ids, type = "bytes"))
  # Each actual reporting key includes every K value, not only its auxiliary j.
  aux_key_upper <- vapply(aux$descriptors, function(d) {
    bytes("reporting_operator:") + bytes(as.character(d$g)) + bytes(as.character(d$K)) + max(1, length(d$K))
  }, numeric(1))
  key_max <- max(nchar(keys, type = "bytes"), aux_key_upper, nchar(scalar_ids[-1L], type = "bytes"))
  expected_bytes <- 128 * common + 16 * id_max + 16 * key_max
  expect_equal(inv$records, length(keys) + length(aux$descriptors))
  expect_equal(inv$common_name_bytes, common)
  expect_equal(inv$max_id_bytes, id_max)
  expect_equal(inv$max_key_bytes, key_max)
  expect_equal(inv$per_record_source_bytes, expected_bytes)
  expect_equal(inv$per_record_slots, ceiling(expected_bytes / 2) + 32)
  expect_equal(inv$total_record_slots, inv$records * (ceiling(expected_bytes / 2) + 32))
  expect_gt(bytes(long_key), nchar(long_key, type = "chars"))
  for (d in aux$descriptors) {
    actual_key <- paste0("reporting_operator:", gsub(" ", "", paste(d$g, paste(d$K, collapse = ","), sep = ":")))
    expect_lte(bytes(actual_key), inv$max_key_bytes)
  }
  expect_lte(unname(cp04_payload(list(guard = inv, construction = inv))["entries"]), inv$metadata_entries)
  larger <- meta
  larger$catalog$rows <- rbind(meta$catalog$rows, meta$catalog$rows)
  larger$row_mapping$functional_key <- rep(keys, 2L)
  twice <- ipw_raw_record_char_inventory(larger, frame, aux, ctrl$fitted_control)
  expect_identical(twice$per_record_source_bytes, inv$per_record_source_bytes)
  expect_equal(twice$total_record_slots, (2 * length(keys) + length(aux$descriptors)) * inv$per_record_slots)
  # A source mutation dropping the longest key must fail the same byte truth.
  wrong <- meta; wrong$row_mapping$functional_key <- rep("same-key", length(keys))
  omitted <- ipw_raw_record_char_inventory(wrong, frame, aux, ctrl$fitted_control)
  expect_lt(omitted$per_record_source_bytes, expected_bytes)
  cap <- ctrl$owner_control
  cap$max_matrix_entries <- 4096 + 64 * (ncol(Z) + 4 + length(aux$descriptors) + 1) - 1
  expect_error(ipw_raw_record_char_inventory(meta, frame, aux, ctrl$fitted_control, envelope = cap),
    class = "cdid_ipw_owner_bound_error")
})

test_that("complete NYT guard and constructed raw records share the funded local character inventory", {
  sim <- make_d93_compact_ipw_test_panel(n_units = 64L, times = 1:5,
    cohorts = c(2, 3, 4, Inf), missing_mode = "complete", unequal_weights = TRUE)
  covariate <- paste0("x_", paste(rep("\u03bb\u00e9", 24L), collapse = ""))
  sim$data[[covariate]] <- sim$data$x1
  args <- list(data = sim$data, idname = "id", tname = "time", yname = "y", gname = "g",
    control = make_d93_test_pipeline_control(), xformla = stats::reformulate(covariate), weightsname = "w",
    control_group = "notyettreated", base_period = "varying", pt_assumption = "all", pt_start = 1,
    weighting = "optimal", reduce_structural_redundancy = TRUE,
    run_numerical_transfer = FALSE, cband = FALSE, bstrap = FALSE)
  actual <- do.call(run_chained_ipw_pipeline, args)
  legacy <- do.call(run_legacy_scalar_ipw_pipeline, args)
  reg <- actual$nuisance_registry; b <- reg$raw_ipw
  witness <- actual$resource_accounting$complete_raw_capacity
  expect_identical(witness$admitted, TRUE)
  expect_s3_class(b, "cdid_ipw_raw_backend")
  inv <- b$accounting$construction_character_inventory
  expect_identical(witness$character_inventory, inv)
  expect_identical(witness$construction_table_entries_bound, b$accounting$construction_table_entries_bound)
  compact_entries <- unname(cp04_payload(list(b$pairs, b$support_specs, b$supports,
    b$treatment_specs, b$treatments, b$links))["entries"])
  expect_equal(b$accounting$owned_entries, compact_entries + 256 +
    b$accounting$owner_snapshot_entries + inv$metadata_entries)
  expect_equal(inv$catalog_records, nrow(actual$metadata$catalog$rows))
  expect_gt(inv$auxiliary_records, 0)
  expect_equal(inv$records, inv$catalog_records + inv$auxiliary_records)
  expect_lte(b$accounting$construction_peak_entries, witness$construction_entries)
  expect_lte(actual$resource_accounting$simultaneous_entries, witness$peak_entries)
  # Count every actual retained record, including aliases and auxiliary-only
  # fitted/nuisance records, rather than estimating physical string sharing.
  records <- list(reg$nuisance_views, reg$link_views, reg$fitted_views, reg$point_links_all,
    reg$raw_aux_nuisance_views, reg$raw_aux_fitted_views, reg$raw_aux_points,
    b$pairs, b$support_specs, b$supports, b$treatment_specs, b$treatments, b$links)
  payload <- cp04_payload(records)
  expect_lte(unname(payload["entries"]), witness$construction_table_entries_bound)
  # Static field/disclosure strings remain funded by the unchanged fixed
  # allowance; dynamic own/name/key byte roles are extra and counted fourfold.
  expect_lte(unname(payload["char_bytes"]), 2 * 20480 * inv$records + 4 * inv$records * inv$per_record_source_bytes)
  expect_equal(actual$reporting$estimates, legacy$reporting$estimates, tolerance = 1e-12)
  expect_equal(actual$reporting$if_full, legacy$reporting$if_full, tolerance = 1e-12)
  expect_equal(actual$reporting$vcov, legacy$reporting$vcov, tolerance = 1e-12)
  expect_equal(actual$reporting$se, legacy$reporting$se, tolerance = 1e-12)
  expect_identical(actual$auxiliary_reporting$descriptor_status, legacy$auxiliary_reporting$descriptor_status)
  expect_identical(actual$covariance_handling$handler_used, legacy$covariance_handling$handler_used)
  expect_identical(actual$covariance_handling$handler_invoked, legacy$covariance_handling$handler_invoked)
  expect_equal(actual$covariance_handling$lambda, legacy$covariance_handling$lambda, tolerance = 1e-12)
})
