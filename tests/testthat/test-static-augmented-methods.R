# Targeted test suite for static observation augmentation (selection_est_method = "dr",
# route = "augmented") combined with treatment est_method = "ipw", "reg", and "dr".
# Verifies against independent references (tests/testthat/reference/static_missing_fit_reference.R
# and tests/testthat/reference/static_public_inference_reference.R) that:
#   1. Augmented IPW does not require or fit counterfactual OR (outcome_pooling lazy,
#      nuisance$or NULL, broken outcome_formula untouched) while fitting treatment logits,
#      selection logits, and imputation WLS with complete original-unit IFs.
#   2. Augmented OR (est_method = "reg") does not estimate treatment propensity
#      (nuisance$logit NULL, Q_hat NA, broken/separated treatment design untouched)
#      while reconstructing counterfactual OR from augmented stayer moments/covariances.
#   3. Both pair_logit and independent_period_logit, active pooling regimes
#      (separate, within_group, pooled, custom OR partition with U \ A cell),
#      imputation scales (per_period, per_change), original-unit IFs (including
#      zero- and single-outcome units and finite-difference unit-weight witnesses),
#      analytical/clustered covariance, and runtime nuisance batch reuse obey contracts.

if (!("cdid" %in% names(getLoadedDLLs())) && file.exists("../../src/cdid.so")) {
  dyn.load("../../src/cdid.so")
} else if (!("cdid" %in% names(getLoadedDLLs())) && file.exists("src/cdid.so")) {
  dyn.load("src/cdid.so")
}

if (!exists("att_gt", mode = "function") || !exists("static_missing_fit", mode = "function")) {
  r_dir <- if (dir.exists("../../R")) "../../R" else "R"
  for (f in sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))) {
    source(f)
  }
}

ref_smf_path <- if (file.exists("reference/static_missing_fit_reference.R")) {
  "reference/static_missing_fit_reference.R"
} else {
  "tests/testthat/reference/static_missing_fit_reference.R"
}
if (!exists("ref_static_missing_fit", mode = "function")) {
  source(ref_smf_path)
}

ref_pub_path <- if (file.exists("reference/static_public_inference_reference.R")) {
  "reference/static_public_inference_reference.R"
} else {
  "tests/testthat/reference/static_public_inference_reference.R"
}
if (!exists("ref_static_public_cluster_inference", mode = "function")) {
  source(ref_pub_path)
}

# Small supported 4-period panel fixture with interior observation support for both
# pair_logit and independent_period_logit across all pooling regimes, while retaining
# zero-outcome and single-outcome units in the original-unit panel.
make_augmented_methods_fixture <- function() {
  set.seed(20261105L)
  path_types <- rbind(
    c(0L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 1L),
    c(0L, 0L, 1L, 1L),
    c(0L, 0L, 0L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 0L, 0L, 0L),
    c(1L, 1L, 0L, 0L),
    c(0L, 1L, 0L, 1L)
  )
  n_rep <- 8L
  D_mat <- path_types[rep(seq_len(nrow(path_types)), each = n_rep), , drop = FALSE]
  n <- nrow(D_mat)
  dates <- c(1, 2, 3, 4)

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    x1 <- rnorm(4, sd = 0.6)
    x2 <- runif(4, -0.8, 0.8)
    d_vec <- D_mat[i, ]
    y_vec <- 0.5 * dates + 0.4 * x1 - 0.3 * x2 + 1.1 * d_vec + rnorm(4, sd = 0.25)
    o_vec <- rbinom(4, 1L, 0.78)
    if (i %in% c(1L, 9L, 33L)) o_vec <- c(0L, 0L, 0L, 0L)       # zero-outcome units
    if (i %in% c(2L, 10L, 41L)) o_vec <- c(1L, 0L, 0L, 0L)      # single-outcome units
    if (i %% 4L == 3L || i %% 4L == 0L) o_vec <- c(1L, 1L, 1L, 1L)
    if (i %% 8L == 5L) o_vec <- c(1L, 0L, 1L, 1L)
    if (i %% 8L == 6L) o_vec <- c(0L, 1L, 1L, 0L)
    y_obs <- ifelse(o_vec == 1L, y_vec, NA_real_)
    cl_id <- sprintf("cl_%02d", ((i - 1L) %% 16L) + 1L)
    rows[[i]] <- data.frame(
      id = sprintf("id_%03d", i),
      cluster = cl_id,
      time = dates,
      y = y_obs,
      d = d_vec,
      x1 = x1,
      x2 = x2,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)

  # Ensure strict interior date observation support in every 8-unit path block
  # so independent_period_logit has interior convex-hull overlap even under "separate".
  for (k in 0:7) {
    full_ids <- sprintf("id_%03d", 8L * k + c(3L, 4L, 7L, 8L))
    x1_t2 <- df$x1[df$id %in% full_ids & df$time == 2]
    ord <- order(abs(x1_t2 - stats::median(x1_t2)))
    df$y[df$id == full_ids[ord[1L]] & df$time == 3] <- NA_real_
    df$y[df$id == full_ids[ord[2L]] & df$time == 2] <- NA_real_
  }

  panel <- static_panel_data(df, idname = "id", tname = "time", yname = "y", dname = "d")
  t_cells <- static_transition_cells(panel, t = c(3, 3), d = c(0L, 1L), cell_id = c("c_3_0", "c_3_1"))
  all_cells <- static_transition_cells(panel, t = c(3, 3, 4), d = c(0L, 1L, 0L), cell_id = c("c_3_0", "c_3_1", "c_4_0"))
  edges_123 <- data.frame(s = c(1, 1, 2), u = c(2, 3, 3))
  edges_234 <- data.frame(s = c(2, 2, 3), u = c(3, 4, 4))
  role_graphs <- list(
    c_3_0 = list(G = edges_123, C = edges_123),
    c_3_1 = list(G = edges_234, C = edges_123),
    c_4_0 = list(C = edges_234)
  )
  or_blocks <- c(c_3_0 = "b_entry", c_3_1 = "b_exit", c_4_0 = "b_entry")
  bounds <- gmm_bounds_control(
    max_matrix_entries = 50000L, max_batch_rows = 15L, max_targets = 20L,
    max_open_clusters = 1L, allow_materialize_if = TRUE, max_materialized_if_entries = 50000L
  )
  acc <- static_graph_gmm_accuracy_control()
  glm_ctrl <- list(epsilon = 1e-12, maxit = 100L, trace = FALSE)

  list(
    df = df,
    panel = panel,
    t_cells = t_cells,
    all_cells = all_cells,
    role_graphs = role_graphs,
    or_blocks = or_blocks,
    bounds = bounds,
    acc = acc,
    glm_ctrl = glm_ctrl
  )
}

# Helper to run independent reference matching public att_gt on the fixture
run_public_reference_for_att_gt <- function(res, df, xformla, selection_xformla = NULL) {
  req <- res$request
  pan <- static_panel_data(df, idname = req$idname, tname = req$tname, yname = req$yname, dname = req$dname)
  t_ids <- res$target_id
  t_cells <- static_transition_cells(
    pan,
    t = vapply(res$fit$cells[t_ids], `[[`, numeric(1L), "t_num"),
    d = vapply(res$fit$cells[t_ids], `[[`, integer(1L), "d"),
    cell_id = t_ids
  )
  all_ids <- names(res$fit$all_cells)
  all_cells <- static_transition_cells(
    pan,
    t = vapply(res$fit$all_cells, `[[`, numeric(1L), "t_num"),
    d = vapply(res$fit$all_cells, `[[`, integer(1L), "d"),
    cell_id = all_ids
  )
  role_graphs <- setNames(vector("list", length(all_ids)), all_ids)
  for (cid in all_ids) {
    g_spec <- if (cid %in% t_ids) res$fit$role_graphs$G[[cid]]$block_spec$edges[, c("s", "u"), drop = FALSE] else NULL
    c_spec <- res$fit$role_graphs$C[[cid]]$block_spec$edges[, c("s", "u"), drop = FALSE]
    role_graphs[[cid]] <- if (!is.null(g_spec)) list(G = g_spec, C = c_spec) else list(C = c_spec)
  }
  sel_f <- if (is.null(selection_xformla)) xformla else selection_xformla
  use_or <- req$est_method %in% c("reg", "dr")
  ref_static_missing_fit(
    panel = pan,
    target_cells = t_cells,
    all_cells = if (use_or) all_cells else NULL,
    or_blocks = if (use_or) res$request$or_blocks else NULL,
    role_graphs = role_graphs,
    treatment_formula = if (req$est_method %in% c("ipw", "dr")) xformla else NULL,
    outcome_formula = if (use_or) xformla else NULL,
    selection_formula = sel_f,
    imputation_formula = xformla,
    method = req$est_method,
    route = req$route,
    selection_model = req$selection_model,
    selection_pooling = req$selection_pooling,
    imputation_pooling = req$imputation_pooling,
    imputation_scale = req$imputation_scale,
    weighting = req$gmm_weighting$canonical,
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    rank_tol = 1e-7,
    centered = FALSE,
    aggregate_cells = t_ids,
    aggregate_signs = res$sigma
  )
}

test_that("independent references use zero production cdid symbols", {
  r_dir <- if (dir.exists("../../R")) "../../R" else "R"
  prod_env <- environment(att_gt)
  if (!isNamespace(prod_env)) {
    prod_env <- new.env(parent = asNamespace("stats"))
    for (f in sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))) {
      sys.source(f, envir = prod_env)
    }
  }
  prod_syms <- ls(prod_env, all.names = TRUE)
  expect_true(all(c("att_gt", "static_public_att_gt", "static_missing_fit") %in% prod_syms))

  for (rpath in c(ref_smf_path, ref_pub_path)) {
    ref_env <- new.env(parent = baseenv())
    sys.source(rpath, envir = ref_env)
    ref_funs <- ls(ref_env)
    used_globals <- unique(unlist(lapply(ref_funs, function(fn) {
      if (is.function(ref_env[[fn]])) {
        codetools::findGlobals(ref_env[[fn]], merge = FALSE)$functions
      } else {
        character(0)
      }
    })))
    expect_identical(intersect(used_globals, prod_syms), character(0))
  }
})

test_that("augmented observation correction with est_method = 'ipw' does not require or fit counterfactual OR and matches reference across pair/date logits and active pooling", {
  fx <- make_augmented_methods_fixture()
  obs_counts <- rowSums(fx$panel$O)
  zero_single_idx <- which(obs_counts <= 1L)
  expect_true(sum(obs_counts == 0L) >= 3L)
  expect_true(sum(obs_counts == 1L) >= 3L)

  # 1. Private kernel check: broken outcome_formula and NULL or_blocks succeed under augmented IPW
  mfit_ipw_lazy <- static_missing_fit(
    panel = fx$panel,
    target_cells = fx$t_cells,
    all_cells = NULL,
    or_blocks = NULL,
    role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    treatment_formula = ~ x1 + x2,
    outcome_formula = ~ nonexistent_or_col,
    selection_formula = ~ x1,
    imputation_formula = ~ x2,
    method = "ipw",
    route = "augmented",
    selection_model = "pair_logit",
    selection_pooling = "within_group",
    imputation_pooling = "within_group",
    imputation_scale = "per_period",
    weighting = "optimal",
    glm_control = fx$glm_ctrl,
    rank_tol = 1e-10,
    bounds = fx$bounds,
    accuracy = fx$acc,
    aggregate_cells = c("c_3_0", "c_3_1"),
    aggregate_signs = c(1, -1)
  )
  expect_null(mfit_ipw_lazy$nuisance$or)
  expect_null(mfit_ipw_lazy$or_blocks)
  expect_identical(mfit_ipw_lazy$non_target_cell_ids, character(0L))
  expect_false(is.null(mfit_ipw_lazy$nuisance$logit))
  expect_true(all(is.finite(mfit_ipw_lazy$Q_hat)))
  expect_false(is.null(mfit_ipw_lazy$nuisance$observation$selection_fit))
  expect_false(is.null(mfit_ipw_lazy$nuisance$observation$imputation_fit))

  # 2. Public att_gt grid across pair_logit / independent_period_logit, active pooling,
  #    imputation scales, and optimal / identity GMM weighting
  ipw_grid <- list(
    list(sel_m = "pair_logit", sel_p = "separate", imp_p = "separate", imp_s = "per_period", gmm_w = "optimal"),
    list(sel_m = "pair_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period", gmm_w = "optimal"),
    list(sel_m = "pair_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change", gmm_w = "identity"),
    list(sel_m = "independent_period_logit", sel_p = "separate", imp_p = "within_group", imp_s = "per_period", gmm_w = "optimal"),
    list(sel_m = "independent_period_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period", gmm_w = "optimal"),
    list(sel_m = "independent_period_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change", gmm_w = "optimal")
  )

  for (cfg in ipw_grid) {
    res <- att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1 + x2, selection_xformla = ~ x1,
      data = fx$df,
      est_method = "ipw",
      selection = "conditional_mtar",
      selection_est_method = "dr",
      selection_model = cfg$sel_m,
      selection_pooling = cfg$sel_p,
      imputation_pooling = cfg$imp_p,
      imputation_scale = cfg$imp_s,
      gmm_weighting = cfg$gmm_w,
      outcome_pooling = stop("outcome_pooling must remain lazy under augmented est_method = 'ipw'"),
      switch_times = 3, switchers = "both",
      bstrap = FALSE, cband = FALSE
    )
    expect_s3_class(res, "cdid_static_att_gt")
    expect_identical(res$request$route, "augmented")
    expect_identical(res$request$est_method, "ipw")
    expect_identical(res$request$selection_est_method, "dr")
    expect_identical(res$request$outcome_pooling_mode, "inactive_ipw")
    expect_null(res$request$or_blocks)
    expect_identical(res$request$non_target_or_cells, character(0L))
    expect_null(res$fit$or_partition)
    expect_null(res$fit$nuisance$or)
    expect_false(is.null(res$fit$nuisance$logit))
    expect_true(all(is.finite(res$fit$Q_hat)))

    rfit <- run_public_reference_for_att_gt(res, fx$df, ~ x1 + x2, selection_xformla = ~ x1)
    rinf <- ref_static_public_cluster_inference(phi_mat = rfit$phi_tau)
    rinf_agg <- ref_static_public_cluster_inference(phi_mat = matrix(rfit$aggregate$phi, ncol = 1L))
    expect_lt(max(abs(res$theta_hat - rfit$theta_hat)), 1e-11)
    expect_lt(max(abs(res$att - rfit$tau_hat)), 1e-11)
    expect_lt(max(abs(res$inffunc - rfit$phi_tau)), 1e-10)
    expect_lt(max(abs(res$V_analytical - rinf$V_analytical)), 1e-11)
    expect_lt(max(abs(res$se - rinf$analytical_se)), 1e-11)
    expect_lt(abs(res$aggregate$estimate - rfit$aggregate$delta_hat), 1e-11)
    expect_lt(max(abs(res$aggregate$inffunc - rfit$aggregate$phi)), 1e-10)
    expect_lt(abs(res$aggregate$se - unname(rinf_agg$analytical_se)), 1e-11)
    expect_true(all(is.finite(res$inffunc[zero_single_idx, , drop = FALSE])))
  }
})

test_that("augmented observation correction with est_method = 'reg' does not estimate treatment propensity and matches reference across pair/date logits and active pooling", {
  fx <- make_augmented_methods_fixture()
  zero_single_idx <- which(rowSums(fx$panel$O) <= 1L)

  # 1. Private kernel check: broken treatment_formula is never evaluated under augmented reg
  mfit_reg_lazy <- static_missing_fit(
    panel = fx$panel,
    target_cells = fx$t_cells,
    all_cells = fx$all_cells,
    or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs,
    treatment_formula = ~ nonexistent_trt_col,
    outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1,
    imputation_formula = ~ x2,
    method = "reg",
    route = "augmented",
    selection_model = "pair_logit",
    selection_pooling = "within_group",
    imputation_pooling = "within_group",
    imputation_scale = "per_period",
    weighting = "optimal",
    glm_control = fx$glm_ctrl,
    rank_tol = 1e-10,
    bounds = fx$bounds,
    accuracy = fx$acc,
    aggregate_cells = c("c_3_0", "c_3_1"),
    aggregate_signs = c(1, -1)
  )
  expect_null(mfit_reg_lazy$nuisance$logit)
  expect_true(all(is.na(mfit_reg_lazy$Q_hat)))
  expect_false(is.null(mfit_reg_lazy$nuisance$or))
  expect_identical(mfit_reg_lazy$non_target_cell_ids, "c_4_0")
  expect_true(all(vapply(mfit_reg_lazy$gmm$block_meta, function(b) identical(b$role, "G"), logical(1L))))

  # 2. Treatment propensity separation invariance: add a covariate x_sep that perfectly
  #    separates switchers G_c from stayers C_c in R_c = 1 while varying smoothly within C_c and G_c.
  #    Augmented IPW and DR must refuse treatment logit separation (cdid_cp04_separation_error),
  #    whereas augmented reg succeeds because it never estimates treatment propensity.
  df_sep <- fx$df
  # In cell t3_entry (baseline t=2, target t=3, d=0): switchers have d(2)=0, d(3)=1.
  sw_ids <- unique(df_sep$id[df_sep$time == 2 & df_sep$d == 0L &
                             df_sep$id %in% df_sep$id[df_sep$time == 3 & df_sep$d == 1L]])
  df_sep$x_sep <- ifelse(df_sep$id %in% sw_ids, 2.0 + 0.2 * df_sep$x1, -2.0 + 0.2 * df_sep$x1)

  for (trt_logit_meth in c("ipw", "dr")) {
    expect_error(
      att_gt(
        yname = "y", tname = "time", idname = "id", dname = "d",
        xformla = ~ x_sep, selection_xformla = ~ x1,
        data = df_sep,
        est_method = trt_logit_meth,
        selection = "conditional_mtar",
        selection_est_method = "dr",
        switch_times = 3, switchers = "entry",
        bstrap = FALSE, cband = FALSE
      ),
      class = "cdid_ipw_fit_error"
    )
  }

  res_sep_reg <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x_sep, selection_xformla = ~ x1,
    data = df_sep,
    est_method = "reg",
    selection = "conditional_mtar",
    selection_est_method = "dr",
    switch_times = 3, switchers = "entry",
    bstrap = FALSE, cband = FALSE
  )
  expect_s3_class(res_sep_reg, "cdid_static_att_gt")
  expect_null(res_sep_reg$fit$nuisance$logit)
  expect_true(all(is.na(res_sep_reg$fit$Q_hat)))
  ref_sep_reg <- run_public_reference_for_att_gt(res_sep_reg, df_sep, ~ x_sep, selection_xformla = ~ x1)
  expect_lt(max(abs(res_sep_reg$att - ref_sep_reg$tau_hat)), 1e-11)
  expect_lt(max(abs(res_sep_reg$inffunc - ref_sep_reg$phi_tau)), 1e-10)

  # 3. Public att_gt grid across pair_logit / independent_period_logit, active pooling
  #    (including custom outcome_pooling with non-target training cell t4_entry in U \ A),
  #    imputation scales, and optimal / identity GMM weighting
  custom_or_pool <- list(
    entry_block = c("t3_entry", "t4_entry"),
    exit_block = "t3_exit"
  )
  reg_grid <- list(
    list(sel_m = "pair_logit", sel_p = "separate", imp_p = "separate", imp_s = "per_period", or_p = "separate", gmm_w = "optimal"),
    list(sel_m = "pair_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period", or_p = "pooled", gmm_w = "optimal"),
    list(sel_m = "pair_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period", or_p = custom_or_pool, gmm_w = "optimal"),
    list(sel_m = "pair_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change", or_p = "separate", gmm_w = "identity"),
    list(sel_m = "independent_period_logit", sel_p = "separate", imp_p = "within_group", imp_s = "per_period", or_p = "separate", gmm_w = "optimal"),
    list(sel_m = "independent_period_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period", or_p = custom_or_pool, gmm_w = "optimal"),
    list(sel_m = "independent_period_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change", or_p = "pooled", gmm_w = "optimal")
  )

  for (cfg in reg_grid) {
    res <- att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1 + x2, selection_xformla = ~ x1,
      data = fx$df,
      est_method = "reg",
      selection = "conditional_mtar",
      selection_est_method = "dr",
      selection_model = cfg$sel_m,
      selection_pooling = cfg$sel_p,
      imputation_pooling = cfg$imp_p,
      imputation_scale = cfg$imp_s,
      outcome_pooling = cfg$or_p,
      gmm_weighting = cfg$gmm_w,
      switch_times = 3, switchers = "both",
      bstrap = FALSE, cband = FALSE
    )
    expect_s3_class(res, "cdid_static_att_gt")
    expect_identical(res$request$route, "augmented")
    expect_identical(res$request$est_method, "reg")
    expect_identical(res$request$selection_est_method, "dr")
    expect_null(res$fit$nuisance$logit)
    expect_true(all(is.na(res$fit$Q_hat)))
    expect_false(is.null(res$fit$nuisance$or))
    if (is.list(cfg$or_p)) {
      expect_identical(res$request$non_target_or_cells, "t4_entry")
    }

    rfit <- run_public_reference_for_att_gt(res, fx$df, ~ x1 + x2, selection_xformla = ~ x1)
    rinf <- ref_static_public_cluster_inference(phi_mat = rfit$phi_tau)
    rinf_agg <- ref_static_public_cluster_inference(phi_mat = matrix(rfit$aggregate$phi, ncol = 1L))
    expect_lt(max(abs(res$theta_hat - rfit$theta_hat)), 1e-11)
    expect_lt(max(abs(res$att - rfit$tau_hat)), 1e-11)
    expect_lt(max(abs(res$inffunc - rfit$phi_tau)), 1e-10)
    expect_lt(max(abs(res$V_analytical - rinf$V_analytical)), 1e-11)
    expect_lt(max(abs(res$se - rinf$analytical_se)), 1e-11)
    expect_lt(abs(res$aggregate$estimate - rfit$aggregate$delta_hat), 1e-11)
    expect_lt(max(abs(res$aggregate$inffunc - rfit$aggregate$phi)), 1e-10)
    expect_lt(abs(res$aggregate$se - unname(rinf_agg$analytical_se)), 1e-11)
    expect_lt(max(abs(res$fit$nuisance$or$alpha_hat[names(rfit$or_summary$alpha_hat)] - rfit$or_summary$alpha_hat)), 1e-11)
    for (bid in names(rfit$or_summary$b_hat)) {
      expect_lt(max(abs(res$fit$nuisance$or$b_hat[[bid]] - rfit$or_summary$b_hat[[bid]])), 1e-11)
    }
    expect_true(all(is.finite(res$inffunc[zero_single_idx, , drop = FALSE])))
  }
})

test_that("finite-difference unit-weight IF witnesses and clustered inference match for augmented ipw and reg", {
  fx <- make_augmented_methods_fixture()
  # Test units 1:10 (includes id_001, id_009 with 0 observed outcomes and id_002, id_010 with 1 observed outcome)
  test_units <- 1:10

  fd_configs <- list(
    list(meth = "ipw", sel_m = "pair_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period"),
    list(meth = "ipw", sel_m = "independent_period_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change"),
    list(meth = "reg", sel_m = "pair_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period"),
    list(meth = "reg", sel_m = "independent_period_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change")
  )

  for (cfg in fd_configs) {
    use_or <- identical(cfg$meth, "reg")
    fd_res <- ref_static_missing_fd_if(
      units = test_units,
      eps = 1e-6,
      hold_W_fixed = TRUE,
      panel = fx$panel,
      target_cells = fx$t_cells,
      all_cells = if (use_or) fx$all_cells else NULL,
      or_blocks = if (use_or) fx$or_blocks else NULL,
      role_graphs = if (use_or) fx$role_graphs else fx$role_graphs[c("c_3_0", "c_3_1")],
      treatment_formula = if (use_or) NULL else ~ x1 + x2,
      outcome_formula = if (use_or) ~ x1 + x2 else NULL,
      selection_formula = ~ x1,
      imputation_formula = ~ x2,
      method = cfg$meth,
      route = "augmented",
      selection_model = cfg$sel_m,
      selection_pooling = cfg$sel_p,
      imputation_pooling = cfg$imp_p,
      imputation_scale = cfg$imp_s,
      weighting = "optimal",
      glm_control = fx$glm_ctrl,
      rank_tol = 1e-10,
      aggregate_cells = c("c_3_0", "c_3_1"),
      aggregate_signs = c(1, -1)
    )
    mfit <- static_missing_fit(
      panel = fx$panel,
      target_cells = fx$t_cells,
      all_cells = if (use_or) fx$all_cells else NULL,
      or_blocks = if (use_or) fx$or_blocks else NULL,
      role_graphs = if (use_or) fx$role_graphs else fx$role_graphs[c("c_3_0", "c_3_1")],
      treatment_formula = if (use_or) NULL else ~ x1 + x2,
      outcome_formula = if (use_or) ~ x1 + x2 else NULL,
      selection_formula = ~ x1,
      imputation_formula = ~ x2,
      method = cfg$meth,
      route = "augmented",
      selection_model = cfg$sel_m,
      selection_pooling = cfg$sel_p,
      imputation_pooling = cfg$imp_p,
      imputation_scale = cfg$imp_s,
      weighting = "optimal",
      glm_control = fx$glm_ctrl,
      rank_tol = 1e-10,
      bounds = fx$bounds,
      accuracy = fx$acc,
      aggregate_cells = c("c_3_0", "c_3_1"),
      aggregate_signs = c(1, -1)
    )
    expect_lt(max(abs(mfit$phi_theta[test_units, ] - fd_res$fd_theta)), 1e-7)
    expect_lt(max(abs(mfit$aggregate$phi[test_units] - fd_res$fd_agg)), 1e-7)
  }

  # Clustered inference check on public att_gt for augmented ipw and reg
  cl_map <- unique(fx$df[, c("id", "cluster")])
  cl_map <- cl_map[match(fx$panel$unit_id, cl_map$id), , drop = FALSE]
  cl_vec <- cl_map$cluster

  for (meth in c("ipw", "reg")) {
    res_uncl <- att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1 + x2, selection_xformla = ~ x1,
      data = fx$df,
      est_method = meth,
      selection = "conditional_mtar",
      selection_est_method = "dr",
      selection_model = "pair_logit",
      selection_pooling = "within_group",
      imputation_pooling = "within_group",
      imputation_scale = "per_period",
      switch_times = 3, switchers = "both",
      bstrap = FALSE, cband = FALSE
    )
    res_cl <- att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1 + x2, selection_xformla = ~ x1,
      data = fx$df, clustervars = "cluster",
      est_method = meth,
      selection = "conditional_mtar",
      selection_est_method = "dr",
      selection_model = "pair_logit",
      selection_pooling = "within_group",
      imputation_pooling = "within_group",
      imputation_scale = "per_period",
      switch_times = 3, switchers = "both",
      bstrap = FALSE, cband = FALSE
    )
    expect_identical(res_cl$att, res_uncl$att)
    expect_identical(res_cl$theta_hat, res_uncl$theta_hat)
    expect_identical(res_cl$inffunc, res_uncl$inffunc)
    expect_identical(res_cl$aggregate$estimate, res_uncl$aggregate$estimate)
    expect_identical(res_cl$aggregate$inffunc, res_uncl$aggregate$inffunc)

    ref_cl <- ref_static_public_cluster_inference(
      phi_mat = res_cl$inffunc,
      cluster_vec = cl_vec,
      signs_mat = NULL,
      alp = 0.05
    )
    ref_cl_agg <- ref_static_public_cluster_inference(
      phi_mat = matrix(res_cl$aggregate$inffunc, ncol = 1L, dimnames = list(NULL, "__agg__")),
      cluster_vec = cl_vec,
      signs_mat = NULL,
      alp = 0.05
    )
    expect_lt(max(abs(res_cl$V_analytical - ref_cl$V_analytical)), 1e-12)
    expect_lt(max(abs(res_cl$se - ref_cl$analytical_se)), 1e-12)
    expect_lt(abs(res_cl$aggregate$se - unname(ref_cl_agg$analytical_se)), 1e-12)
    expect_gt(max(abs(res_cl$V_analytical - res_uncl$V_analytical)), 1e-6)
  }
})

test_that("runtime nuisance batch reuse is preserved for augmented ipw and reg and matches declined-reuse results", {
  fx <- make_augmented_methods_fixture()

  ns <- environment(static_observation_adapter_fit)
  orig_nuisance_if <- get("static_nuisance_if", envir = ns)
  was_locked <- base::bindingIsLocked("static_nuisance_if", ns)
  if (was_locked) base::unlockBinding("static_nuisance_if", ns)
  on.exit({
    if (base::bindingIsLocked("static_nuisance_if", ns)) base::unlockBinding("static_nuisance_if", ns)
    assign("static_nuisance_if", orig_nuisance_if, envir = ns)
    if (was_locked) base::lockBinding("static_nuisance_if", ns)
  }, add = TRUE)

  sel_calls <- 0L
  imp_calls <- 0L
  assign("static_nuisance_if", function(fit, ids) {
    if (identical(fit$kind, "logit")) {
      sel_calls <<- sel_calls + 1L
    } else if (identical(fit$kind, "wls")) {
      imp_calls <<- imp_calls + 1L
    }
    orig_nuisance_if(fit, ids)
  }, envir = ns)
  if (was_locked) base::lockBinding("static_nuisance_if", ns)

  n_batches_per_pass <- as.integer(ceiling(fx$panel$n / fx$bounds$max_batch_rows))
  expect_identical(n_batches_per_pass, 5L)

  for (meth in c("ipw", "reg")) {
    use_or <- identical(meth, "reg")
    run_mfit <- function(bounds_arg) {
      static_missing_fit(
        panel = fx$panel,
        target_cells = fx$t_cells,
        all_cells = if (use_or) fx$all_cells else NULL,
        or_blocks = if (use_or) fx$or_blocks else NULL,
        role_graphs = if (use_or) fx$role_graphs else fx$role_graphs[c("c_3_0", "c_3_1")],
        treatment_formula = if (use_or) NULL else ~ x1 + x2,
        outcome_formula = if (use_or) ~ x1 + x2 else NULL,
        selection_formula = ~ x1,
        imputation_formula = ~ x2,
        method = meth,
        route = "augmented",
        selection_model = "independent_period_logit",
        selection_pooling = "within_group",
        imputation_pooling = "within_group",
        imputation_scale = "per_period",
        weighting = "optimal",
        glm_control = fx$glm_ctrl,
        rank_tol = 1e-10,
        bounds = bounds_arg,
        accuracy = fx$acc,
        materialize_if = TRUE,
        aggregate_cells = c("c_3_0", "c_3_1"),
        aggregate_signs = c(1, -1)
      )
    }

    sel_calls <- 0L; imp_calls <- 0L
    fit_funded <- run_mfit(fx$bounds)
    expect_identical(sel_calls, 4L * n_batches_per_pass)
    expect_identical(imp_calls, 4L * n_batches_per_pass)

    env_funded <- environment(fit_funded$link_if_callback)
    expect_true(isTRUE(env_funded$reuse_nuisance_if))
    budget_info <- env_funded$nuisance_batch_budget

    env_ad <- environment(fit_funded$nuisance$observation$with_if_batch)
    expect_false(env_ad$batch_active)
    expect_null(env_ad$batch_sel_if)
    expect_null(env_ad$batch_imp_if)

    bounds_tight <- fx$bounds
    bounds_tight$max_matrix_entries <- as.integer(floor(budget_info$consumer_entries))
    expect_lt(bounds_tight$max_matrix_entries, budget_info$total_entries)

    sel_calls <- 0L; imp_calls <- 0L
    fit_declined <- run_mfit(bounds_tight)
    env_declined <- environment(fit_declined$link_if_callback)
    expect_false(isTRUE(env_declined$reuse_nuisance_if))
    expect_gt(sel_calls, 4L * n_batches_per_pass)
    expect_gt(imp_calls, 4L * n_batches_per_pass)

    expect_equal(fit_funded$theta_hat, fit_declined$theta_hat, tolerance = 1e-14)
    expect_equal(fit_funded$tau_hat, fit_declined$tau_hat, tolerance = 1e-14)
    expect_equal(fit_funded$cov_theta, fit_declined$cov_theta, tolerance = 1e-14)
    expect_equal(fit_funded$phi_theta, fit_declined$phi_theta, tolerance = 1e-14)
    expect_equal(fit_funded$aggregate$delta_hat, fit_declined$aggregate$delta_hat, tolerance = 1e-14)
    expect_equal(fit_funded$aggregate$phi, fit_declined$aggregate$phi, tolerance = 1e-14)
  }
})
