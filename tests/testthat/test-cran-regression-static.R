# Bounded static (`dname`) regression suite executed unconditionally by
# tests/testthat.R (filter = "cran-regression") alongside test-cran-regression.R
# during standard R CMD check. Exercises the installed exported `att_gt(dname = ...)`
# interface on a small hand-computable static binary fixture with entries, exits,
# re-entry, supported stayers, and zero/single-outcome units across IPW, OR, and DR
# treatment methods and baseline, observation-weighting, imputation, and augmented
# missing-outcome routes, verifying point estimates, full-unit influence functions,
# analytical/clustered covariance, role-graph support, refusal gates, and citations
# against independent base-R references. All 185 test files in tests/testthat/
# run unfiltered (filter = NULL) in CI and via validation/run_full_tests.R.

static_cran_ref_env <- new.env(parent = baseenv())
sys.source(
  testthat::test_path("reference", "static_missing_fit_reference.R"),
  envir = static_cran_ref_env
)
sys.source(
  testthat::test_path("reference", "static_public_inference_reference.R"),
  envir = static_cran_ref_env
)

# Pure base-R panel, transition-cell, and role-graph constructor for the
# independent reference `ref_static_missing_fit()` (uses zero internal cdid symbols).
build_static_cran_ref_inputs <- function(df,
                                         target_ids,
                                         all_ids = target_ids,
                                         idname = "id",
                                         tname = "time",
                                         yname = "y",
                                         dname = "d") {
  dates_num <- sort(unique(as.double(df[[tname]])))
  uids <- sort(unique(as.character(df[[idname]])))
  n <- length(uids)
  T_len <- length(dates_num)

  row_index <- matrix(NA_integer_, nrow = n, ncol = T_len, dimnames = list(uids, as.character(dates_num)))
  Y_mat <- matrix(NA_real_, nrow = n, ncol = T_len, dimnames = list(uids, as.character(dates_num)))
  O_mat <- matrix(0L, nrow = n, ncol = T_len, dimnames = list(uids, as.character(dates_num)))
  D_mat <- matrix(NA_integer_, nrow = n, ncol = T_len, dimnames = list(uids, as.character(dates_num)))

  u_pos <- match(as.character(df[[idname]]), uids)
  t_pos <- match(as.double(df[[tname]]), dates_num)
  for (r in seq_len(nrow(df))) {
    i <- u_pos[r]
    ti <- t_pos[r]
    row_index[i, ti] <- r
    y_val <- as.double(df[[yname]][r])
    if (is.finite(y_val)) {
      Y_mat[i, ti] <- y_val
      O_mat[i, ti] <- 1L
    }
    D_mat[i, ti] <- as.integer(df[[dname]][r])
  }

  panel <- list(
    n = n,
    unit_id = uids,
    dates_num = dates_num,
    row_index = row_index,
    data = df,
    Y = Y_mat,
    O = O_mat,
    D = D_mat
  )

  all_cand <- list()
  for (ti in seq.int(2L, T_len)) {
    t_num <- dates_num[ti]
    t_prev_num <- dates_num[ti - 1L]
    for (d_val in c(0L, 1L)) {
      cid <- sprintf("t%s_%s", as.character(t_num), if (d_val == 0L) "entry" else "exit")
      R_vec <- as.integer(D_mat[, ti - 1L] == d_val)
      G_vec <- as.integer(R_vec == 1L & D_mat[, ti] == (1L - d_val))
      C_vec <- as.integer(R_vec == 1L & D_mat[, ti] == d_val)
      all_cand[[cid]] <- list(
        cell_id = cid,
        t_num = t_num,
        t_prev_num = t_prev_num,
        t_idx = ti,
        t_prev_idx = ti - 1L,
        d = d_val,
        sigma = as.double(1L - 2L * d_val),
        R = R_vec,
        G = G_vec,
        C = C_vec
      )
    }
  }

  build_edges <- function(mask_vec) {
    idx <- which(mask_vec == 1L)
    s_vec <- numeric(0L)
    u_vec <- numeric(0L)
    for (a in seq_len(T_len - 1L)) {
      for (b in seq.int(a + 1L, T_len)) {
        if (any(O_mat[idx, a] == 1L & O_mat[idx, b] == 1L)) {
          s_vec <- c(s_vec, dates_num[a])
          u_vec <- c(u_vec, dates_num[b])
        }
      }
    }
    data.frame(s = s_vec, u = u_vec, stringsAsFactors = FALSE)
  }

  union_ids <- unique(c(target_ids, all_ids))
  role_graphs <- list()
  for (cid in union_ids) {
    cl <- all_cand[[cid]]
    if (cid %in% target_ids) {
      role_graphs[[cid]] <- list(
        G = build_edges(cl$G),
        C = build_edges(cl$C)
      )
    } else {
      role_graphs[[cid]] <- list(
        C = build_edges(cl$C)
      )
    }
  }

  list(
    panel = panel,
    target_cells = all_cand[target_ids],
    all_cells = all_cand[union_ids],
    role_graphs = role_graphs
  )
}

# Deterministic hand-computable static binary fixture with entries, exits,
# re-entry (0 -> 1 -> 0 -> 1), supported stayers (0 -> 0 -> 0 -> 0 and
# 1 -> 1 -> 1 -> 1), and zero- and single-outcome units in the incomplete panel.
make_static_cran_fixture <- function() {
  path_types <- rbind(
    c(0L, 0L, 0L, 0L), # always-0 stayer
    c(0L, 1L, 1L, 1L), # entry at t = 2, then 1-stayer
    c(0L, 0L, 1L, 1L), # entry at t = 3, then 1-stayer
    c(0L, 0L, 0L, 1L), # entry at t = 4
    c(1L, 1L, 1L, 1L), # always-1 stayer
    c(1L, 0L, 0L, 0L), # exit at t = 2, then 0-stayer
    c(1L, 1L, 0L, 0L), # exit at t = 3, then 0-stayer
    c(0L, 1L, 0L, 1L)  # entry at t = 2, exit at t = 3, re-entry at t = 4
  )
  n_rep <- 8L
  D_mat <- path_types[rep(seq_len(nrow(path_types)), each = n_rep), , drop = FALSE]
  n <- nrow(D_mat)
  dates <- c(1, 2, 3, 4)

  # Deterministic covariate sequences and unit fixed effects:
  # replicates 1 and 2 span the extreme min/max of x1 and x2 and are always
  # fully observed (1, 1, 1, 1), guaranteeing strict interior overlap for every
  # pair/date observation logit and treatment logit in every path and role.
  x1_base <- c(-1.20, 1.20, -0.50, 0.50, -0.20, 0.20, -0.10, 0.10)
  x2_base <- c(-1.10, 1.10, 0.60, -0.60, -0.15, 0.25, 0.10, -0.10)
  obs_templates <- rbind(
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 0L, 1L),
    c(1L, 0L, 1L, 1L),
    c(0L, 1L, 1L, 1L),
    c(1L, 0L, 0L, 0L), # single-outcome unit in paths 1-2; (1,1,1,0) in paths 3-8
    c(0L, 0L, 0L, 0L)  # zero-outcome unit in paths 1-2; (1,1,1,1) in paths 3-8
  )

  rows_hand <- vector("list", n)
  rows_comp <- vector("list", n)
  rows_miss <- vector("list", n)

  # True static treatment effect by period: tau_2 = 1.5, tau_3 = 2.0, tau_4 = 2.5
  tau_by_t <- c(0.0, 1.5, 2.0, 2.5)
  trend_by_t <- c(1.0, 1.5, 2.25, 3.25)

  for (i in seq_len(n)) {
    p_idx <- ((i - 1L) %/% n_rep) + 1L
    r_idx <- ((i - 1L) %% n_rep) + 1L
    d_vec <- D_mat[i, ]
    alpha_i <- 2.0 * p_idx + 0.5 * r_idx

    # Hand-computable outcome: Y_{it}(d) = alpha_i + trend_t + tau_t * d_{it} + u_{it}
    # where within each path/replicate pair u_{it} cancels symmetrically across replicates.
    u_sign <- if (r_idx %% 2L == 1L) 1.0 else -1.0
    u_vec <- u_sign * 0.125 * c(1.0, -1.0, 0.5, -0.5)
    y_hand <- alpha_i + trend_by_t + tau_by_t * d_vec + u_vec

    # Covariate-driven outcome for regression/propensity/missingness checks
    x1_vec <- x1_base[r_idx] + 0.03 * (p_idx - 4.5) + 0.02 * dates
    x2_vec <- x2_base[r_idx] - 0.02 * (p_idx - 4.5) + 0.01 * (dates - 2.5)
    y_cov <- alpha_i + trend_by_t + 0.40 * x1_vec - 0.25 * x2_vec +
      tau_by_t * d_vec + u_vec * (1.0 + 0.2 * r_idx)

    o_vec <- obs_templates[r_idx, ]
    # Keep zero-outcome and single-outcome units in the first two path blocks,
    # and observed pairs in the remaining blocks so every role/date has interior support.
    if (p_idx > 2L && r_idx >= 7L) {
      o_vec <- if (r_idx == 7L) c(1L, 1L, 1L, 0L) else c(1L, 1L, 1L, 1L)
    }

    uid <- sprintf("u_%02d", i)
    cl_id <- sprintf("cl_%02d", ((i - 1L) %% 8L) + 1L)

    rows_hand[[i]] <- data.frame(
      id = uid, time = dates, y = y_hand, d = d_vec,
      x1 = x1_vec, x2 = x2_vec, cluster = cl_id, stringsAsFactors = FALSE
    )
    rows_comp[[i]] <- data.frame(
      id = uid, time = dates, y = y_cov, d = d_vec,
      x1 = x1_vec, x2 = x2_vec, cluster = cl_id, stringsAsFactors = FALSE
    )
    rows_miss[[i]] <- data.frame(
      id = uid, time = dates,
      y = ifelse(o_vec == 1L, y_cov, NA_real_),
      d = d_vec, x1 = x1_vec, x2 = x2_vec, cluster = cl_id,
      stringsAsFactors = FALSE
    )
  }

  list(
    hand = do.call(rbind, rows_hand),
    complete = do.call(rbind, rows_comp),
    missing = do.call(rbind, rows_miss)
  )
}

test_that("hand-computable static fixture recovers exact entry, exit, re-entry, and aggregate effects and IFs", {
  fix <- make_static_cran_fixture()
  df_h <- fix$hand

  # Expected transition cells at t = 2, 3, 4:
  # - t2_entry (0 -> 1 at t = 2): paths 2 and 8 (16 switchers) vs paths 1, 3, 4 (24 stayers)
  #   True tau_{2,0} = tau_2 - tau_1 = 1.5; theta_{2,0} = +1.5
  # - t2_exit  (1 -> 0 at t = 2): path 6 (8 switchers) vs paths 5, 7 (16 stayers)
  #   Switchers lose D=1 at t=2 while stayers keep D=1 (effect tau_2=1.5 vs tau_1=0):
  #   theta_{2,1} = -1.5, sigma_{2,1} = -1 => tau_{2,1} = +1.5
  # - t3_entry (0 -> 1 at t = 3): path 3 (8 switchers) vs paths 1, 4, 6 (24 stayers)
  #   theta_{3,0} = +2.0, tau_{3,0} = +2.0
  # - t3_exit  (1 -> 0 at t = 3): paths 7 and 8 (16 switchers) vs paths 2, 5 (16 stayers)
  #   Switchers move from D_2=1 (tau_2=1.5) to D_3=0 (0), stayers move from D_2=1 to D_3=1 (tau_3=2.0):
  #   theta_{3,1} = -2.0, sigma_{3,1} = -1 => tau_{3,1} = +2.0
  # - t4_entry (0 -> 1 at t = 4, including re-entry path 8): paths 4 and 8 (16 switchers) vs paths 1, 6, 7 (24 stayers)
  #   theta_{4,0} = +2.5, tau_{4,0} = +2.5
  expected_targets <- c("t2_entry", "t2_exit", "t3_entry", "t3_exit", "t4_entry")
  expected_tau <- c(t2_entry = 1.5, t2_exit = 1.5, t3_entry = 2.0, t3_exit = 2.0, t4_entry = 2.5)
  expected_sigma <- c(t2_entry = 1.0, t2_exit = -1.0, t3_entry = 1.0, t3_exit = -1.0, t4_entry = 1.0)
  expected_theta <- expected_sigma * expected_tau
  expected_n_sw <- c(t2_entry = 16L, t2_exit = 8L, t3_entry = 8L, t3_exit = 16L, t4_entry = 16L)
  expected_n_st <- c(t2_entry = 24L, t2_exit = 16L, t3_entry = 24L, t3_exit = 16L, t4_entry = 24L)
  expected_agg <- sum(expected_n_sw * expected_tau) / sum(expected_n_sw) # (16*1.5 + 8*1.5 + 8*2.0 + 16*2.0 + 16*2.5) / 64 = 124 / 64 = 1.9375

  # Hand-computed closed-form unit influence functions phi_{tau, c}(i) and aggregate IF phi_delta(i)
  ref_h <- build_static_cran_ref_inputs(df_h, target_ids = expected_targets)
  n_h <- ref_h$panel$n
  hand_phi_tau <- matrix(0.0, nrow = n_h, ncol = length(expected_targets),
                         dimnames = list(ref_h$panel$unit_id, expected_targets))
  G_mat_h <- matrix(0L, nrow = n_h, ncol = length(expected_targets),
                    dimnames = list(ref_h$panel$unit_id, expected_targets))
  for (cid in expected_targets) {
    cl <- ref_h$target_cells[[cid]]
    dy_c <- ref_h$panel$Y[, cl$t_idx] - ref_h$panel$Y[, cl$t_prev_idx]
    g_mean <- mean(dy_c[cl$G == 1L])
    c_mean <- mean(dy_c[cl$C == 1L])
    phi_theta_c <- (n_h / sum(cl$G)) * cl$G * (dy_c - g_mean) -
      (n_h / sum(cl$C)) * cl$C * (dy_c - c_mean)
    hand_phi_tau[, cid] <- cl$sigma * phi_theta_c
    G_mat_h[, cid] <- cl$G
  }
  w_sw <- expected_n_sw / sum(expected_n_sw)
  hand_agg_if <- drop(hand_phi_tau %*% w_sw) +
    drop(sweep(n_h * G_mat_h, 2L, expected_n_sw, FUN = "-") %*% ((expected_tau - expected_agg) / sum(expected_n_sw)))
  hand_phi_cent <- sweep(hand_phi_tau, 2L, colMeans(hand_phi_tau), FUN = "-")
  hand_V_an <- crossprod(hand_phi_cent) / (as.double(n_h)^2)
  hand_se_an <- sqrt(diag(hand_V_an))

  for (meth in c("ipw", "reg", "dr")) {
    res <- cdid::att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ 1, data = df_h, est_method = meth,
      bstrap = FALSE, cband = FALSE
    )
    expect_s3_class(res, "cdid_static_att_gt")
    expect_identical(res$target_id, expected_targets)
    expect_identical(res$n_switchers, expected_n_sw)
    expect_identical(res$n_stayers, expected_n_st)
    expect_equal(res$att, expected_tau, tolerance = 1e-12)
    expect_equal(res$theta_hat, expected_theta, tolerance = 1e-12)
    expect_equal(res$inffunc, hand_phi_tau, tolerance = 1e-11)
    expect_equal(res$V_analytical, hand_V_an, tolerance = 1e-12)
    expect_equal(res$se, hand_se_an, tolerance = 1e-12)
    expect_equal(res$aggregate$estimate, expected_agg, tolerance = 1e-12)
    expect_equal(as.double(res$aggregate$inffunc), as.double(hand_agg_if), tolerance = 1e-11)
    expect_identical(res$aggregate$total_switchers, 64L)
    expect_equal(res$aggregate$weights, w_sw, tolerance = 1e-14)
  }
})

test_that("static dname route matches independent reference across IPW/reg/DR and baseline, weighting, imputation, and augmented missingness", {
  fix <- make_static_cran_fixture()
  df_m <- fix$missing

  # Verify zero- and single-outcome units are retained in the study population
  obs_counts <- tapply(!is.na(df_m$y), df_m$id, sum)
  expect_identical(as.integer(sum(obs_counts == 0L)), 2L)
  expect_identical(as.integer(sum(obs_counts == 1L)), 2L)

  # Representative specifications covering:
  # - baseline (selection = "none") across IPW, reg, DR (separate and pooled OR)
  # - observation weighting (selection_est_method = "ipw") with pair_logit and independent_period_logit
  # - change imputation (selection_est_method = "imputation") with per_period and per_change
  # - augmented correction (selection_est_method = "dr") across all three treatment methods IPW, reg, DR
  #   plus custom outcome_pooling partition with a non-target stayer training cell
  specs <- list(
    list(
      label = "baseline_ipw",
      est_method = "ipw", selection = "none", selection_est_method = "ipw",
      route = "baseline", xformla = ~ x1, selection_xformla = NULL,
      selection_model = "pair_logit", selection_pooling = "separate",
      imputation_pooling = "separate", imputation_scale = "per_period",
      outcome_pooling = "separate", switchers = "both", switch_times = c(2, 3)
    ),
    list(
      label = "baseline_reg_pooled",
      est_method = "reg", selection = "none", selection_est_method = "ipw",
      route = "baseline", xformla = ~ x1 + x2, selection_xformla = NULL,
      selection_model = "pair_logit", selection_pooling = "separate",
      imputation_pooling = "separate", imputation_scale = "per_period",
      outcome_pooling = "pooled", switchers = "both", switch_times = NULL
    ),
    list(
      label = "baseline_dr_separate",
      est_method = "dr", selection = "none", selection_est_method = "ipw",
      route = "baseline", xformla = ~ x1, selection_xformla = NULL,
      selection_model = "pair_logit", selection_pooling = "separate",
      imputation_pooling = "separate", imputation_scale = "per_period",
      outcome_pooling = "separate", switchers = "both", switch_times = NULL
    ),
    list(
      label = "weighting_pair_logit_ipw",
      est_method = "ipw", selection = "conditional_mtar", selection_est_method = "ipw",
      route = "weighting", xformla = ~ x1, selection_xformla = ~ x2,
      selection_model = "pair_logit", selection_pooling = "within_group",
      imputation_pooling = "separate", imputation_scale = "per_period",
      outcome_pooling = "separate", switchers = "both", switch_times = c(2, 3)
    ),
    list(
      label = "weighting_indep_period_dr",
      est_method = "dr", selection = "conditional_mtar", selection_est_method = "ipw",
      route = "weighting", xformla = ~ x1, selection_xformla = ~ x1,
      selection_model = "independent_period_logit", selection_pooling = "pooled",
      imputation_pooling = "separate", imputation_scale = "per_period",
      outcome_pooling = "separate", switchers = "entry", switch_times = c(2, 3, 4)
    ),
    list(
      label = "imputation_per_period_reg",
      est_method = "reg", selection = "conditional_mtar", selection_est_method = "imputation",
      route = "imputation", xformla = ~ x1 + x2, selection_xformla = NULL,
      selection_model = "pair_logit", selection_pooling = "separate",
      imputation_pooling = "within_group", imputation_scale = "per_period",
      outcome_pooling = "separate", switchers = "both", switch_times = c(2, 3)
    ),
    list(
      label = "imputation_per_change_dr",
      est_method = "dr", selection = "conditional_mtar", selection_est_method = "imputation",
      route = "imputation", xformla = ~ x1, selection_xformla = NULL,
      selection_model = "pair_logit", selection_pooling = "separate",
      imputation_pooling = "pooled", imputation_scale = "per_change",
      outcome_pooling = "pooled", switchers = "exit", switch_times = c(2, 3)
    ),
    list(
      label = "augmented_ipw",
      est_method = "ipw", selection = "conditional_mtar", selection_est_method = "dr",
      route = "augmented", xformla = ~ x1, selection_xformla = ~ x2,
      selection_model = "pair_logit", selection_pooling = "within_group",
      imputation_pooling = "within_group", imputation_scale = "per_period",
      outcome_pooling = "separate", switchers = "both", switch_times = c(3, 4),
      # t4_exit is empty so request only t = 3 with switchers = "both"
      override_switch_times = 3
    ),
    list(
      label = "augmented_reg_partition",
      est_method = "reg", selection = "conditional_mtar", selection_est_method = "dr",
      route = "augmented", xformla = ~ x1, selection_xformla = ~ x2,
      selection_model = "pair_logit", selection_pooling = "pooled",
      imputation_pooling = "within_group", imputation_scale = "per_change",
      outcome_pooling = c(t2_entry = "b_en", t3_entry = "b_en", t4_entry = "b_en"),
      switchers = "entry", switch_times = c(2, 3)
    ),
    list(
      label = "augmented_dr",
      est_method = "dr", selection = "conditional_mtar", selection_est_method = "dr",
      route = "augmented", xformla = ~ x1, selection_xformla = ~ x2,
      selection_model = "independent_period_logit", selection_pooling = "within_group",
      imputation_pooling = "pooled", imputation_scale = "per_period",
      outcome_pooling = "separate", switchers = "both", switch_times = c(2, 3)
    )
  )

  for (sp in specs) {
    sw_times <- if (!is.null(sp$override_switch_times)) sp$override_switch_times else sp$switch_times
    res <- cdid::att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = sp$xformla,
      selection_xformla = sp$selection_xformla,
      data = df_m,
      est_method = sp$est_method,
      selection = sp$selection,
      selection_est_method = sp$selection_est_method,
      selection_model = sp$selection_model,
      selection_pooling = sp$selection_pooling,
      imputation_pooling = sp$imputation_pooling,
      imputation_scale = sp$imputation_scale,
      outcome_pooling = sp$outcome_pooling,
      switchers = sp$switchers,
      switch_times = sw_times,
      bstrap = FALSE, cband = FALSE
    )

    expect_s3_class(res, "cdid_static_att_gt")
    expect_identical(res$n, 64L)
    expect_identical(nrow(res$inffunc), 64L)
    expect_identical(res$diagnostics$n_zero_outcome_units, 2L)
    expect_identical(res$diagnostics$n_single_outcome_units, 2L)
    expect_true(all(res$available))
    expect_true(all(is.finite(res$inffunc)))

    t_ids <- res$target_id
    all_ids <- if (sp$est_method %in% c("reg", "dr") && is.character(sp$outcome_pooling) && length(sp$outcome_pooling) > 1L) {
      names(sp$outcome_pooling)
    } else {
      t_ids
    }
    or_blocks_ref <- if (sp$est_method %in% c("reg", "dr")) {
      if (identical(sp$outcome_pooling, "separate")) {
        stats::setNames(t_ids, t_ids)
      } else if (identical(sp$outcome_pooling, "pooled")) {
        stats::setNames(rep("pooled", length(t_ids)), t_ids)
      } else {
        sp$outcome_pooling
      }
    } else {
      NULL
    }

    ref_in <- build_static_cran_ref_inputs(df_m, target_ids = t_ids, all_ids = all_ids)
    ref_fit <- static_cran_ref_env$ref_static_missing_fit(
      panel = ref_in$panel,
      target_cells = ref_in$target_cells,
      all_cells = ref_in$all_cells,
      role_graphs = ref_in$role_graphs,
      or_blocks = or_blocks_ref,
      method = sp$est_method,
      route = sp$route,
      treatment_formula = if (sp$est_method %in% c("ipw", "dr")) sp$xformla else NULL,
      outcome_formula = if (sp$est_method %in% c("reg", "dr")) sp$xformla else NULL,
      selection_formula = if (sp$route %in% c("weighting", "augmented")) {
        if (is.null(sp$selection_xformla)) sp$xformla else sp$selection_xformla
      } else NULL,
      imputation_formula = if (sp$route %in% c("imputation", "augmented")) sp$xformla else NULL,
      selection_model = if (sp$route %in% c("weighting", "augmented")) sp$selection_model else NULL,
      selection_pooling = if (sp$route %in% c("weighting", "augmented")) sp$selection_pooling else NULL,
      imputation_pooling = if (sp$route %in% c("imputation", "augmented")) sp$imputation_pooling else NULL,
      imputation_scale = if (sp$route %in% c("imputation", "augmented")) sp$imputation_scale else NULL,
      weighting = "optimal",
      glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
      rank_tol = 1e-7,
      aggregate_cells = t_ids,
      aggregate_signs = res$sigma
    )

    expect_equal(res$att, ref_fit$tau_hat[t_ids], tolerance = 1e-10, info = sp$label)
    expect_equal(res$theta_hat, ref_fit$theta_hat[t_ids], tolerance = 1e-10, info = sp$label)
    expect_equal(res$inffunc, ref_fit$phi_tau[, t_ids, drop = FALSE], tolerance = 1e-9, info = sp$label)
    expect_equal(res$aggregate$estimate, ref_fit$aggregate$delta_hat, tolerance = 1e-10, info = sp$label)
    expect_equal(as.double(res$aggregate$inffunc), as.double(ref_fit$aggregate$phi), tolerance = 1e-9, info = sp$label)

    ref_inf <- static_cran_ref_env$ref_static_public_cluster_inference(
      phi_mat = res$inffunc,
      cluster_vec = NULL,
      signs_mat = NULL,
      alp = 0.05
    )
    expect_equal(res$V_analytical, ref_inf$V_analytical, tolerance = 1e-12, info = sp$label)
    expect_equal(res$se, ref_inf$analytical_se, tolerance = 1e-12, info = sp$label)
  }
})

test_that("static dname route supports clustered and multiplier inference, support graphs, refusal gates, and citations", {
  fix <- make_static_cran_fixture()
  df_m <- fix$missing

  # 1. Unit-level vs clustered inference and multiplier bootstrap
  set.seed(20261107L)
  res_cl <- cdid::att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m, clustervars = "cluster",
    est_method = "dr", switch_times = c(2, 3), switchers = "both",
    bstrap = TRUE, cband = TRUE, biters = 100
  )
  res_uncl <- cdid::att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m,
    est_method = "dr", switch_times = c(2, 3), switchers = "both",
    bstrap = FALSE, cband = FALSE
  )
  expect_equal(res_cl$att, res_uncl$att, tolerance = 1e-14)
  expect_equal(res_cl$inffunc, res_uncl$inffunc, tolerance = 1e-14)
  expect_equal(res_cl$aggregate$estimate, res_uncl$aggregate$estimate, tolerance = 1e-14)

  # Match clustered analytical covariance against independent cluster reference
  u_order <- match(res_cl$unit_id, df_m$id[df_m$time == 1])
  cl_vec <- df_m$cluster[df_m$time == 1][u_order]
  ref_cl <- static_cran_ref_env$ref_static_public_cluster_inference(
    phi_mat = res_cl$inffunc,
    cluster_vec = cl_vec,
    signs_mat = NULL,
    alp = 0.05
  )
  expect_equal(res_cl$V_analytical, ref_cl$V_analytical, tolerance = 1e-12)
  expect_equal(res_cl$analytical_se, ref_cl$analytical_se, tolerance = 1e-12)
  expect_true(all(is.finite(res_cl$bootstrap_se)))
  expect_true(is.finite(res_cl$cband_crit_val))
  expect_true(is.finite(res_cl$aggregate$se))

  # 2. Explicit empty cell (t4_exit has 0 switchers) is retained as unavailable
  #    and marks the switcher-frequency aggregate unavailable
  res_empty <- cdid::att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m,
    switchers = "exit", switch_times = c(3, 4),
    bstrap = FALSE, cband = FALSE
  )
  expect_identical(res_empty$target_id, c("t3_exit", "t4_exit"))
  expect_identical(unname(res_empty$available), c(TRUE, FALSE))
  expect_identical(unname(res_empty$support_status), c("ok", "empty_switchers"))
  expect_true(is.finite(res_empty$att[["t3_exit"]]))
  expect_true(is.na(res_empty$att[["t4_exit"]]))
  expect_false(isTRUE(res_empty$aggregate$available))
  expect_true(is.na(res_empty$aggregate$estimate))

  # 3. Role-graph summary and tidy/glance/summary/citation contracts
  expect_identical(res_cl$graphs$rule, "all_positive_observed_count_pairs_per_role_no_pruning")
  expect_true(all(res_cl$graphs$summary$endpoints_connected))
  td <- cdid::tidy(res_cl)
  gl <- cdid::glance(res_cl)
  sm <- summary(res_cl)
  out_lines <- paste(utils::capture.output(print(sm)), collapse = "\n")
  expect_s3_class(td, "data.frame")
  expect_identical(nrow(td), length(res_cl$target_id))
  expect_s3_class(gl, "data.frame")
  expect_identical(nrow(gl), 1L)
  expect_identical(names(res_cl$references), c("estimator", "cdid", "static"))
  expect_match(out_lines, "10\\.1016/j\\.jeconom\\.2024\\.105783")
  expect_match(out_lines, "arXiv:2201\\.06898v7")
  expect_no_match(out_lines, "Callaway")

  # 4. Key static refusal gates
  expect_error(
    cdid::att_gt(yname = "y", tname = "time", idname = "id", gname = "d", dname = "d", data = df_m),
    class = "cdid_static_routing_error"
  )
  expect_error(
    cdid::aggte(res_cl, type = "simple"),
    regexp = "MP object"
  )
  df_bad_d <- df_m
  df_bad_d$d[1L] <- NA_integer_
  expect_error(
    cdid::att_gt(yname = "y", tname = "time", idname = "id", dname = "d", data = df_bad_d),
    class = "cdid_static_missing_treatment_error"
  )
})
