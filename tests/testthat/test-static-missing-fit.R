# Focused test suite for private end-to-end static missing-outcome assembly
# (R/static_missing_fit.R against independent tests/testthat/reference/static_missing_fit_reference.R).

if (!("cdid" %in% names(getLoadedDLLs())) && file.exists("../../src/cdid.so")) {
  dyn.load("../../src/cdid.so")
} else if (!("cdid" %in% names(getLoadedDLLs())) && file.exists("src/cdid.so")) {
  dyn.load("src/cdid.so")
}

if (!exists("static_missing_fit", mode = "function")) {
  r_dir <- if (dir.exists("../../R")) "../../R" else "R"
  for (f in sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))) {
    source(f)
  }
}

if (!exists("ref_static_missing_fit", mode = "function")) {
  ref_path <- if (file.exists("reference/static_missing_fit_reference.R")) {
    "reference/static_missing_fit_reference.R"
  } else {
    "tests/testthat/reference/static_missing_fit_reference.R"
  }
  source(ref_path)
}

build_saved_recovered_fixture <- function() {
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
    rows[[i]] <- data.frame(
      id = sprintf("id_%03d", i), time = dates, y = y_obs, d = d_vec,
      x1 = x1, x2 = x2, stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
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
    edges_123 = edges_123,
    edges_234 = edges_234,
    role_graphs = role_graphs,
    or_blocks = or_blocks,
    bounds = bounds,
    acc = acc,
    glm_ctrl = glm_ctrl
  )
}

# Companion fixture with interior date observation indicators (0 < sum(O) < n_role)
# in every role population at every date 1..4, so independent_period_logit has
# interior support across all pooling regimes while retaining 0- and 1-outcome units.
build_interior_date_fixture <- function() {
  fx <- build_saved_recovered_fixture()
  df <- fx$df
  # In the saved fixture, lines 52-54 set o_vec[3] = 1 for all i %% 8 in {3,4,5,6,7,0}.
  # For blocks where i %% 8 in {1,2} also drew rbinom = 1 at date 3 (specifically
  # c_3_0__G on id_017..id_024), date 3 is 100% observed (8/8), causing boundary
  # coordinate separation in independent_period_logit.
  # For each 8-unit block k = 0..7, among the 4 fully observed units (i %% 4 in {3, 0}),
  # pick the unit with median x1 at time == 2 and set its outcome at time == 3 to NA,
  # and for the second-closest unit to the median x1, if date 2 or date 4 has <= 1 missing
  # in the block, set one interior-x1 unit to NA so every date cell in every role
  # population has strict interior convex-hull overlap in (1, x1) even under separate pooling.
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
  fx$df <- df
  fx$panel <- panel
  fx$t_cells <- t_cells
  fx$all_cells <- all_cells
  fx
}

test_that("independent reference uses zero production cdid symbols", {
  ref_path <- if (file.exists("reference/static_missing_fit_reference.R")) {
    "reference/static_missing_fit_reference.R"
  } else {
    "tests/testthat/reference/static_missing_fit_reference.R"
  }
  r_dir <- if (dir.exists("../../R")) "../../R" else "R"

  ref_env <- new.env(parent = baseenv())
  sys.source(ref_path, envir = ref_env)
  prod_env <- environment(static_missing_fit)
  if (!isNamespace(prod_env)) {
    prod_env <- new.env(parent = baseenv())
    for (f in sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))) {
      sys.source(f, envir = prod_env)
    }
  }
  prod_syms <- ls(prod_env, all.names = TRUE)
  expect_true("static_missing_fit" %in% prod_syms)
  ref_funs <- ls(ref_env)
  used_globals <- unique(unlist(lapply(ref_funs, function(fn) {
    if (is.function(ref_env[[fn]])) {
      codetools::findGlobals(ref_env[[fn]], merge = FALSE)$functions
    } else {
      character(0)
    }
  })))
  leaked <- intersect(used_globals, prod_syms)
  expect_identical(leaked, character(0))
})

test_that("saved 21 route/method checks and 12-unit FD witness match contract and retain exact separation refusal", {
  fx <- build_saved_recovered_fixture()

  # Verify retained zero-outcome and single-outcome units in panel
  obs_counts <- rowSums(fx$panel$O)
  expect_true(sum(obs_counts == 0L) >= 3L)
  expect_true(sum(obs_counts == 1L) >= 3L)

  # Confirm why independent_period_logit legitimately refuses on the saved fixture:
  # in c_3_0 switcher role G (8 units, id_017..id_024), date 3 is observed for all 8/8 units (no 0s),
  # so the date-cell intercept alpha_{c_3_0__G::d3} has boundary separation.
  g_mask_30 <- fx$t_cells$cells$c_3_0$G == 1L
  expect_identical(sum(fx$panel$O[g_mask_30, 3]), sum(g_mask_30))

  configs <- list(
    list(route = "baseline", sel_m = NULL, sel_p = NULL, imp_p = NULL, imp_s = NULL, expect_sep = FALSE),
    list(route = "weighting", sel_m = "pair_logit", sel_p = "separate", imp_p = NULL, imp_s = NULL, expect_sep = FALSE),
    list(route = "weighting", sel_m = "independent_period_logit", sel_p = "within_group", imp_p = NULL, imp_s = NULL, expect_sep = TRUE),
    list(route = "imputation", sel_m = NULL, sel_p = NULL, imp_p = "within_group", imp_s = "per_period", expect_sep = FALSE),
    list(route = "imputation", sel_m = NULL, sel_p = NULL, imp_p = "pooled", imp_s = "per_change", expect_sep = FALSE),
    list(route = "augmented", sel_m = "pair_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period", expect_sep = FALSE),
    list(route = "augmented", sel_m = "independent_period_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change", expect_sep = TRUE)
  )

  for (cfg in configs) {
    for (meth in c("ipw", "reg", "dr")) {
      if (cfg$expect_sep) {
        expect_error(
          static_missing_fit(
            panel = fx$panel,
            target_cells = fx$t_cells,
            all_cells = fx$all_cells,
            or_blocks = fx$or_blocks,
            role_graphs = fx$role_graphs,
            treatment_formula = ~ x1 + x2,
            outcome_formula = ~ x1 + x2,
            selection_formula = ~ x1,
            imputation_formula = ~ x2,
            method = meth,
            route = cfg$route,
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
          ),
          class = "cdid_cp04_separation_error"
        )
      } else {
        mfit <- static_missing_fit(
          panel = fx$panel,
          target_cells = fx$t_cells,
          all_cells = fx$all_cells,
          or_blocks = fx$or_blocks,
          role_graphs = fx$role_graphs,
          treatment_formula = ~ x1 + x2,
          outcome_formula = ~ x1 + x2,
          selection_formula = ~ x1,
          imputation_formula = ~ x2,
          method = meth,
          route = cfg$route,
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
        rfit <- ref_static_missing_fit(
          panel = fx$panel,
          target_cells = fx$t_cells,
          all_cells = if (meth == "ipw") NULL else fx$all_cells,
          or_blocks = if (meth == "ipw") NULL else fx$or_blocks,
          role_graphs = fx$role_graphs,
          treatment_formula = ~ x1 + x2,
          outcome_formula = ~ x1 + x2,
          selection_formula = ~ x1,
          imputation_formula = ~ x2,
          method = meth,
          route = cfg$route,
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
        expect_lt(max(abs(mfit$theta_hat - rfit$theta_hat)), 1e-11)
        expect_lt(max(abs(mfit$tau_hat - rfit$tau_hat)), 1e-11)
        expect_lt(max(abs(mfit$phi_theta - rfit$phi_theta)), 1e-10)
        expect_lt(max(abs(mfit$phi_tau - rfit$phi_tau)), 1e-10)
        expect_lt(max(abs(mfit$cov_theta - rfit$cov_theta)), 1e-11)
        expect_lt(max(abs(mfit$se_theta - rfit$se_theta)), 1e-11)
        expect_lt(abs(mfit$aggregate$delta_hat - rfit$aggregate$delta_hat), 1e-11)
        expect_lt(max(abs(mfit$aggregate$phi - rfit$aggregate$phi)), 1e-10)
        expect_lt(abs(mfit$aggregate$se - rfit$aggregate$se), 1e-11)
      }
    }
  }

  # 12-unit finite-difference unit-weight perturbation witness (holding W fixed at two-step optimal W)
  fd_res <- ref_static_missing_fd_if(
    units = 1:12,
    eps = 1e-6,
    hold_W_fixed = TRUE,
    panel = fx$panel,
    target_cells = fx$t_cells,
    all_cells = fx$all_cells,
    or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs,
    treatment_formula = ~ x1 + x2,
    outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1,
    imputation_formula = ~ x2,
    method = "dr",
    route = "augmented",
    selection_model = "pair_logit",
    selection_pooling = "within_group",
    imputation_pooling = "within_group",
    imputation_scale = "per_period",
    weighting = "optimal",
    glm_control = fx$glm_ctrl,
    rank_tol = 1e-10,
    aggregate_cells = c("c_3_0", "c_3_1"),
    aggregate_signs = c(1, -1)
  )
  mfit_dr_aug <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs, treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "pair_logit", selection_pooling = "within_group",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    bounds = fx$bounds, accuracy = fx$acc, aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )
  expect_lt(max(abs(mfit_dr_aug$phi_theta[1:12, ] - fd_res$fd_theta)), 1e-7)
  expect_lt(max(abs(mfit_dr_aug$aggregate$phi[1:12] - fd_res$fd_agg)), 1e-7)
})

test_that("selected independent_period_logit and pooling/scale combinations match reference and FD witness on interior-date fixture", {
  fx <- build_interior_date_fixture()

  ind_configs <- list(
    list(route = "weighting", sel_m = "independent_period_logit", sel_p = "separate", imp_p = NULL, imp_s = NULL),
    list(route = "weighting", sel_m = "independent_period_logit", sel_p = "within_group", imp_p = NULL, imp_s = NULL),
    list(route = "weighting", sel_m = "independent_period_logit", sel_p = "pooled", imp_p = NULL, imp_s = NULL),
    list(route = "weighting", sel_m = "pair_logit", sel_p = "pooled", imp_p = NULL, imp_s = NULL),
    list(route = "imputation", sel_m = NULL, sel_p = NULL, imp_p = "separate", imp_s = "per_period"),
    list(route = "augmented", sel_m = "independent_period_logit", sel_p = "within_group", imp_p = "within_group", imp_s = "per_period"),
    list(route = "augmented", sel_m = "independent_period_logit", sel_p = "pooled", imp_p = "pooled", imp_s = "per_change")
  )

  for (cfg in ind_configs) {
    for (meth in c("ipw", "reg", "dr")) {
      mfit <- static_missing_fit(
        panel = fx$panel,
        target_cells = fx$t_cells,
        all_cells = fx$all_cells,
        or_blocks = fx$or_blocks,
        role_graphs = fx$role_graphs,
        treatment_formula = ~ x1 + x2,
        outcome_formula = ~ x1 + x2,
        selection_formula = ~ x1,
        imputation_formula = ~ x2,
        method = meth,
        route = cfg$route,
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
      rfit <- ref_static_missing_fit(
        panel = fx$panel,
        target_cells = fx$t_cells,
        all_cells = fx$all_cells,
        or_blocks = fx$or_blocks,
        role_graphs = fx$role_graphs,
        treatment_formula = ~ x1 + x2,
        outcome_formula = ~ x1 + x2,
        selection_formula = ~ x1,
        imputation_formula = ~ x2,
        method = meth,
        route = cfg$route,
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
      expect_lt(max(abs(mfit$theta_hat - rfit$theta_hat)), 1e-11)
      expect_lt(max(abs(mfit$phi_theta - rfit$phi_theta)), 1e-10)
      expect_lt(max(abs(mfit$cov_theta - rfit$cov_theta)), 1e-11)
      expect_lt(abs(mfit$aggregate$delta_hat - rfit$aggregate$delta_hat), 1e-11)
      expect_lt(max(abs(mfit$aggregate$phi - rfit$aggregate$phi)), 1e-10)
    }
  }

  # FD witness on augmented + independent_period_logit + pooled + per_change (units 1:6 including 0- and 1-outcome units)
  fd_ind <- ref_static_missing_fd_if(
    units = 1:6,
    eps = 1e-6,
    hold_W_fixed = TRUE,
    panel = fx$panel,
    target_cells = fx$t_cells,
    all_cells = fx$all_cells,
    or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs,
    treatment_formula = ~ x1 + x2,
    outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1,
    imputation_formula = ~ x2,
    method = "dr",
    route = "augmented",
    selection_model = "independent_period_logit",
    selection_pooling = "pooled",
    imputation_pooling = "pooled",
    imputation_scale = "per_change",
    weighting = "optimal",
    glm_control = fx$glm_ctrl,
    rank_tol = 1e-10,
    aggregate_cells = c("c_3_0", "c_3_1"),
    aggregate_signs = c(1, -1)
  )
  mfit_ind <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs, treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "independent_period_logit", selection_pooling = "pooled",
    imputation_pooling = "pooled", imputation_scale = "per_change",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    bounds = fx$bounds, accuracy = fx$acc, aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )
  expect_lt(max(abs(mfit_ind$phi_theta[1:6, ] - fd_ind$fd_theta)), 1e-7)
  expect_lt(max(abs(mfit_ind$aggregate$phi[1:6] - fd_ind$fd_agg)), 1e-7)
})

test_that("bounded complete-data reduction matches static_complete_fit and refuses weighting/augmented on complete panels", {
  fx <- build_saved_recovered_fixture()
  df_comp <- fx$df
  # Make panel 100% observed
  set.seed(20261106L)
  df_comp$y <- 0.5 * df_comp$time + 0.4 * df_comp$x1 - 0.3 * df_comp$x2 + 1.1 * df_comp$d + rnorm(nrow(df_comp), sd = 0.25)
  panel_comp <- static_panel_data(df_comp, idname = "id", tname = "time", yname = "y", dname = "d")

  t_cells_comp <- static_transition_cells(panel_comp, t = c(3, 3, 4), d = c(0L, 1L, 0L), cell_id = c("c_3_0", "c_3_1", "c_4_0"))
  role_graphs_comp <- list(
    c_3_0 = list(G = fx$edges_123, C = fx$edges_234),
    c_3_1 = list(G = fx$edges_234, C = fx$edges_123),
    c_4_0 = list(G = fx$edges_234, C = fx$edges_234)
  )
  or_blocks_comp <- c(c_3_0 = "b_entry", c_3_1 = "b_exit", c_4_0 = "b_entry")

  # Build corresponding static_complete_cell objects
  cd_cells <- lapply(names(t_cells_comp$cells), function(cid) {
    cl <- t_cells_comp$cells[[cid]]
    des_t <- static_extract_reference_design(panel_comp, cl, ~ x1 + x2, role = "treatment", mask = rep(TRUE, panel_comp$n))
    des_o <- static_extract_reference_design(panel_comp, cl, ~ x1 + x2, role = "outcome", mask = rep(TRUE, panel_comp$n))
    dy_c <- panel_comp$Y[, cl$t_idx] - panel_comp$Y[, cl$t_prev_idx]
    static_complete_cell(
      cell_id = cid,
      d = cl$d,
      t = cl$t_num,
      unit_id = panel_comp$unit_id,
      R = cl$R,
      G = cl$G,
      C = cl$C,
      dy = dy_c,
      X_tilde = des_t$X_tilde,
      Z_or = des_o$Z,
      or_block = unname(or_blocks_comp[[cid]])
    )
  })
  names(cd_cells) <- names(t_cells_comp$cells)

  for (w_mode in c("identity", "optimal")) {
    for (meth in c("ipw", "reg", "dr")) {
      mfit <- static_missing_fit(
        panel = panel_comp,
        target_cells = t_cells_comp,
        or_blocks = if (meth == "ipw") NULL else or_blocks_comp,
        role_graphs = role_graphs_comp,
        treatment_formula = ~ x1 + x2,
        outcome_formula = ~ x1 + x2,
        method = meth,
        route = "baseline",
        weighting = w_mode,
        glm_control = fx$glm_ctrl,
        rank_tol = 1e-10,
        bounds = fx$bounds,
        accuracy = fx$acc,
        aggregate_cells = c("c_3_0", "c_3_1", "c_4_0"),
        aggregate_signs = c(1, -1, 1)
      )
      cfit <- static_complete_fit(
        cells = cd_cells,
        method = meth,
        or_blocks = if (meth == "ipw") NULL else or_blocks_comp,
        glm_control = fx$glm_ctrl,
        rank_tol = 1e-10
      )
      cagg <- static_complete_aggregate(
        fit = cfit,
        cells = c("c_3_0", "c_3_1", "c_4_0"),
        signs = c(1, -1, 1)
      )

      expect_lt(max(abs(mfit$theta_hat - cfit$theta_hat)), 1e-11)
      expect_lt(max(abs(mfit$tau_hat - cfit$tau_hat)), 1e-11)
      expect_lt(max(abs(mfit$phi_theta - cfit$phi_theta)), 1e-10)
      expect_lt(max(abs(mfit$cov_theta - cfit$cov_theta)), 1e-11)
      expect_lt(max(abs(mfit$se_theta - cfit$se_theta)), 1e-11)
      expect_lt(abs(mfit$aggregate$delta_hat - cagg$delta_hat), 1e-11)
      expect_lt(max(abs(mfit$aggregate$phi - cagg$phi)), 1e-10)
      expect_lt(abs(mfit$aggregate$se - cagg$se), 1e-11)

      if (meth != "ipw") {
        expect_lt(max(abs(mfit$nuisance$or$alpha_hat - cfit$nuisance$or$alpha_hat)), 1e-11)
        for (bid in names(cfit$nuisance$or$b_hat)) {
          expect_lt(max(abs(mfit$nuisance$or$b_hat[[bid]] - cfit$nuisance$or$b_hat[[bid]])), 1e-11)
        }
      }
    }
  }

  # Complete panel + weighting or augmented route refuses with separation error
  for (bad_route in c("weighting", "augmented")) {
    expect_error(
      static_missing_fit(
        panel = panel_comp,
        target_cells = t_cells_comp,
        or_blocks = or_blocks_comp,
        role_graphs = role_graphs_comp,
        formula = ~ x1,
        method = "dr",
        route = bad_route,
        selection_model = "pair_logit",
        selection_pooling = "separate",
        imputation_pooling = if (bad_route == "augmented") "separate" else NULL,
        imputation_scale = if (bad_route == "augmented") "per_period" else NULL,
        weighting = "identity",
        glm_control = fx$glm_ctrl,
        rank_tol = 1e-10,
        bounds = fx$bounds,
        accuracy = fx$acc
      ),
      class = "cdid_cp04_separation_error"
    )
  }

  # Intercept-only reduction: DR == IPW == OR under fixed identity weighting (LS flow)
  fit_ipw_int <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    formula = ~ 1, method = "ipw", route = "baseline", weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  fit_dr_int <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, or_blocks = "separate",
    role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    formula = ~ 1, method = "dr", route = "baseline", weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  fit_or_int <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, or_blocks = "separate",
    role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    formula = ~ 1, method = "reg", route = "baseline", weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  expect_lt(max(abs(fit_ipw_int$theta_hat - fit_dr_int$theta_hat)), 1e-12)
  expect_lt(max(abs(fit_ipw_int$phi_theta - fit_dr_int$phi_theta)), 1e-11)
  expect_lt(max(abs(fit_ipw_int$theta_hat - fit_or_int$theta_hat)), 1e-12)
  expect_lt(max(abs(fit_ipw_int$phi_theta - fit_or_int$phi_theta)), 1e-11)
})

test_that("inference-only clustering leaves point estimates, weights, and unit IFs unchanged while updating sandwich SEs", {
  fx <- build_saved_recovered_fixture()
  # Non-contiguous cluster labels across the 64 units
  cluster_vec <- sprintf("cl_%02d", ((seq_len(fx$panel$n) * 7L) %% 11L) + 1L)

  fit_uncl <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs, treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "pair_logit", selection_pooling = "within_group",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    bounds = fx$bounds, accuracy = fx$acc, cluster = NULL,
    aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )
  fit_cl <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs, treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "pair_logit", selection_pooling = "within_group",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    bounds = fx$bounds, accuracy = fx$acc, cluster = cluster_vec,
    aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )
  ref_cl <- ref_static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs, treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "pair_logit", selection_pooling = "within_group",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    cluster = cluster_vec, aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )

  # Exact invariance of point estimates, moment covariance, GMM loadings, and unit IFs
  expect_identical(fit_cl$theta_hat, fit_uncl$theta_hat)
  expect_identical(fit_cl$tau_hat, fit_uncl$tau_hat)
  expect_identical(fit_cl$mhat, fit_uncl$mhat)
  expect_identical(fit_cl$moment_covariance, fit_uncl$moment_covariance)
  expect_identical(fit_cl$moment_weights$Lambda, fit_uncl$moment_weights$Lambda)
  expect_identical(fit_cl$phi_theta, fit_uncl$phi_theta)
  expect_identical(fit_cl$phi_tau, fit_uncl$phi_tau)
  expect_identical(fit_cl$aggregate$delta_hat, fit_uncl$aggregate$delta_hat)
  expect_identical(fit_cl$aggregate$phi, fit_uncl$aggregate$phi)

  # Clustered covariance and SEs match independent reference and differ from unclustered
  expect_lt(max(abs(fit_cl$cov_theta - ref_cl$cov_theta)), 1e-12)
  expect_lt(max(abs(fit_cl$se_theta - ref_cl$se_theta)), 1e-12)
  expect_lt(abs(fit_cl$aggregate$variance - ref_cl$aggregate$variance), 1e-12)
  expect_lt(abs(fit_cl$aggregate$se - ref_cl$aggregate$se), 1e-12)
  expect_true(max(abs(fit_cl$cov_theta - fit_uncl$cov_theta)) > 1e-6)
})

test_that("ID-vs-position callback permutation, chunking, and unit-order equivariance hold", {
  fx <- build_saved_recovered_fixture()

  bounds_small <- gmm_bounds_control(
    max_matrix_entries = 50000L, max_batch_rows = 5L, max_targets = 20L,
    max_open_clusters = 1L, allow_materialize_if = FALSE, max_materialized_if_entries = 0L
  )
  fit_stream <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs, treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "pair_logit", selection_pooling = "within_group",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    bounds = bounds_small, accuracy = fx$acc, materialize_if = FALSE,
    aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )
  fit_dense <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs, treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "pair_logit", selection_pooling = "within_group",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    bounds = fx$bounds, accuracy = fx$acc, materialize_if = TRUE,
    aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )

  expect_null(fit_stream$phi_theta)
  expect_null(fit_stream$aggregate$phi)
  expect_lt(max(abs(fit_stream$theta_hat - fit_dense$theta_hat)), 1e-14)
  expect_lt(max(abs(fit_stream$cov_theta - fit_dense$cov_theta)), 1e-14)
  expect_lt(abs(fit_stream$aggregate$delta_hat - fit_dense$aggregate$delta_hat), 1e-14)
  expect_lt(abs(fit_stream$aggregate$se - fit_dense$aggregate$se), 1e-14)

  # Query callbacks by permuted integer indices vs character unit_ids vs logical mask
  perm_idx <- c(42L, 1L, 19L, 9L, 64L, 2L, 33L, 8L)
  perm_ids <- fx$panel$unit_id[perm_idx]

  cb_pos <- fit_stream$cell_if_callback(perm_idx)
  cb_ids <- fit_stream$cell_if_callback(perm_ids)
  expect_identical(cb_pos$phi_theta, cb_ids$phi_theta)
  expect_identical(cb_pos$phi_tau, cb_ids$phi_tau)
  expect_lt(max(abs(cb_pos$phi_theta - fit_dense$phi_theta[perm_idx, , drop = FALSE])), 1e-14)

  l_pos <- fit_stream$link_if_callback(perm_idx)
  l_ids <- fit_stream$link_if_callback(perm_ids)
  expect_identical(l_pos, l_ids)

  cc_pos <- fit_stream$centered_links$c_3_0$if_callback(perm_idx)
  cc_ids <- fit_stream$centered_links$c_3_0$if_callback(perm_ids)
  expect_identical(cc_pos$chi_stacked, cc_ids$chi_stacked)

  agg_pos <- fit_stream$aggregate$if_callback(perm_idx)
  agg_ids <- fit_stream$aggregate$if_callback(perm_ids)
  expect_identical(agg_pos$phi, agg_ids$phi)
  expect_lt(max(abs(agg_pos$phi - fit_dense$aggregate$phi[perm_idx])), 1e-14)

  # Unit-order permutation equivariance on panel construction
  set.seed(20261107L)
  uid_perm <- sample(fx$panel$unit_id)
  df_shuf <- fx$df[sample(nrow(fx$df)), , drop = FALSE]
  fit_perm <- static_missing_fit(
    panel = df_shuf, idname = "id", tname = "time", yname = "y", dname = "d", unit_ids = uid_perm,
    target_cells = list(list(cell_id = "c_3_0", t = 3, d = 0L), list(cell_id = "c_3_1", t = 3, d = 1L)),
    all_cells = list(list(cell_id = "c_3_0", t = 3, d = 0L), list(cell_id = "c_3_1", t = 3, d = 1L), list(cell_id = "c_4_0", t = 4, d = 0L)),
    or_blocks = fx$or_blocks, role_graphs = fx$role_graphs,
    treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1, imputation_formula = ~ x2, method = "dr", route = "augmented",
    selection_model = "pair_logit", selection_pooling = "within_group",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "optimal", glm_control = fx$glm_ctrl, rank_tol = 1e-10,
    bounds = fx$bounds, accuracy = fx$acc,
    aggregate_cells = c("c_3_0", "c_3_1"), aggregate_signs = c(1, -1)
  )
  expect_lt(max(abs(fit_perm$theta_hat - fit_dense$theta_hat)), 1e-12)
  expect_lt(max(abs(fit_perm$cov_theta - fit_dense$cov_theta)), 1e-12)
  expect_lt(max(abs(fit_perm$phi_theta[fx$panel$unit_id, ] - fit_dense$phi_theta)), 1e-11)
  expect_lt(max(abs(fit_perm$aggregate$phi[fx$panel$unit_id] - fit_dense$aggregate$phi)), 1e-11)
})

test_that("additional non-target OR training cell in U \\ A contributes to pooled slope and target IF", {
  fx <- build_saved_recovered_fixture()

  fit_with_u <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells,
    or_blocks = c(c_3_0 = "b_entry", c_3_1 = "b_exit", c_4_0 = "b_entry"),
    role_graphs = fx$role_graphs,
    treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    method = "dr", route = "baseline", weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  fit_target_only <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = NULL,
    or_blocks = c(c_3_0 = "b_entry", c_3_1 = "b_exit"),
    role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    method = "dr", route = "baseline", weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )

  # Non-target cell c_4_0 changes the pooled entry slope b_hat[["b_entry"]] and target c_3_0 estimate
  expect_true(max(abs(fit_with_u$nuisance$or$b_hat$b_entry - fit_target_only$nuisance$or$b_hat$b_entry)) > 1e-3)
  expect_true(abs(fit_with_u$theta_hat[["c_3_0"]] - fit_target_only$theta_hat[["c_3_0"]]) > 1e-4)
  expect_identical(fit_with_u$non_target_cell_ids, "c_4_0")

  # FD witness on units that are stayers in c_4_0 (e.g. path type 1: units 1..8)
  fd_u <- ref_static_missing_fd_if(
    units = 1:8, eps = 1e-6, hold_W_fixed = TRUE,
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells,
    or_blocks = c(c_3_0 = "b_entry", c_3_1 = "b_exit", c_4_0 = "b_entry"),
    role_graphs = fx$role_graphs,
    treatment_formula = ~ x1 + x2, outcome_formula = ~ x1 + x2,
    method = "dr", route = "baseline", weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10
  )
  expect_lt(max(abs(fit_with_u$phi_theta[1:8, ] - fd_u$fd_theta)), 1e-7)
})

test_that("refusal propagation and lazy inactive-formula / active-mask covariate checks hold", {
  fx <- build_saved_recovered_fixture()

  # 1. Lazy inactive formulas: broken formulas for unused nuisances are never evaluated
  fit_reg_lazy <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells, all_cells = fx$all_cells, or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs,
    treatment_formula = ~ nonexistent_trt_col,
    outcome_formula = ~ x1 + x2,
    selection_formula = ~ nonexistent_sel_col,
    imputation_formula = ~ nonexistent_imp_col,
    method = "reg", route = "baseline", weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  expect_s3_class(fit_reg_lazy, "cdid_static_missing_fit")

  fit_ipw_lazy <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells,
    role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    treatment_formula = ~ x1 + x2,
    outcome_formula = ~ nonexistent_or_col,
    selection_formula = ~ x1,
    imputation_formula = ~ nonexistent_imp_col,
    method = "ipw", route = "weighting",
    selection_model = "pair_logit", selection_pooling = "separate",
    weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  expect_s3_class(fit_ipw_lazy, "cdid_static_missing_fit")

  fit_imp_lazy <- static_missing_fit(
    panel = fx$panel, target_cells = fx$t_cells,
    role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    treatment_formula = ~ x1 + x2,
    outcome_formula = ~ nonexistent_or_col,
    selection_formula = ~ nonexistent_sel_col,
    imputation_formula = ~ x2,
    method = "ipw", route = "imputation",
    imputation_pooling = "within_group", imputation_scale = "per_period",
    weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  expect_s3_class(fit_imp_lazy, "cdid_static_missing_fit")

  # 2. Covariate NA outside active risk set R_c is ignored; NA inside active R_c is refused
  df_na_outside <- fx$df
  # At t-1 = 1 (time == 1), no cell has baseline t-1 = 1 (cells have t = 3 -> baseline 2, t = 4 -> baseline 3)
  df_na_outside$x1[df_na_outside$time == 1] <- NA_real_
  panel_na_outside <- static_panel_data(df_na_outside, idname = "id", tname = "time", yname = "y", dname = "d")
  fit_na_ok <- static_missing_fit(
    panel = panel_na_outside,
    target_cells = static_transition_cells(panel_na_outside, t = c(3, 3), d = c(0L, 1L), cell_id = c("c_3_0", "c_3_1")),
    role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
    treatment_formula = ~ x1 + x2, selection_formula = ~ x1,
    method = "ipw", route = "weighting",
    selection_model = "pair_logit", selection_pooling = "separate",
    weighting = "identity",
    glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
  )
  expect_equal(fit_na_ok$theta_hat, fit_ipw_lazy$theta_hat, tolerance = 1e-15)
  expect_equal(fit_na_ok$phi_theta, fit_ipw_lazy$phi_theta, tolerance = 1e-15)

  df_na_active <- fx$df
  df_na_active$x1[df_na_active$id == "id_001" & df_na_active$time == 2] <- NA_real_
  panel_na_act <- static_panel_data(df_na_active, idname = "id", tname = "time", yname = "y", dname = "d")
  expect_error(
    static_missing_fit(
      panel = panel_na_act,
      target_cells = static_transition_cells(panel_na_act, t = c(3, 3), d = c(0L, 1L), cell_id = c("c_3_0", "c_3_1")),
      role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
      treatment_formula = ~ x1 + x2, method = "ipw", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
    ),
    class = "cdid_static_missing_covariate_error"
  )

  # 3. Disconnected role graph or missing baseline/target date is refused
  bad_graphs_disc <- fx$role_graphs
  bad_graphs_disc$c_3_0$G <- data.frame(s = c(1, 3), u = c(2, 4)) # {1,2} and {3,4} disconnected
  expect_error(
    static_missing_fit(
      panel = fx$panel, target_cells = fx$t_cells, role_graphs = bad_graphs_disc[c("c_3_0", "c_3_1")],
      formula = ~ x1, method = "ipw", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
    ),
    class = "cdid_static_graph_disconnected_error"
  )

  bad_graphs_nobase <- fx$role_graphs
  bad_graphs_nobase$c_3_0$G <- data.frame(s = c(1, 3), u = c(3, 4)) # missing baseline t-1 = 2
  expect_error(
    static_missing_fit(
      panel = fx$panel, target_cells = fx$t_cells, role_graphs = bad_graphs_nobase[c("c_3_0", "c_3_1")],
      formula = ~ x1, method = "ipw", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
    ),
    class = "cdid_static_graph_missing_baseline_error"
  )

  # 4. Empty switcher role in requested target cell is refused (never silently dropped)
  # At t = 2, d = 1: check if any unit switches 1 -> 0 at t = 2 (path type 6 has 1,0,0,0 -> switches at t=2;
  # let's zero out d at t=1 for path type 6 so no unit switches 1 -> 0 at t=2)
  df_no_exit2 <- fx$df
  df_no_exit2$d[df_no_exit2$id %in% sprintf("id_%03d", 41:48)] <- 1L
  panel_no_exit2 <- static_panel_data(df_no_exit2, idname = "id", tname = "time", yname = "y", dname = "d")
  cells_empty_sw <- static_transition_cells(panel_no_exit2, t = 2, d = 1L, cell_id = "c_2_1")
  expect_error(
    static_missing_fit(
      panel = panel_no_exit2, target_cells = cells_empty_sw,
      role_graphs = list(c_2_1 = list(G = data.frame(s = 1, u = 2), C = data.frame(s = 1, u = 2))),
      formula = ~ x1, method = "ipw", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
    ),
    class = "cdid_static_support_error"
  )

  # 5. Zero observed-pair mass on a retained edge is refused (never pruned)
  df_zero_edge <- fx$df
  # Wipe out all observed outcomes at time == 1 so edge (1, 2) has zero observed pairs
  df_zero_edge$y[df_zero_edge$time == 1] <- NA_real_
  panel_zero_edge <- static_panel_data(df_zero_edge, idname = "id", tname = "time", yname = "y", dname = "d")
  expect_error(
    static_missing_fit(
      panel = panel_zero_edge,
      target_cells = static_transition_cells(panel_zero_edge, t = c(3, 3), d = c(0L, 1L), cell_id = c("c_3_0", "c_3_1")),
      role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
      formula = ~ x1, method = "ipw", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
    ),
    class = "cdid_static_obs_support_error"
  )

  # 6. Missing or_blocks in reg/dr, missing bounds/accuracy, or single-cluster refusal
  expect_error(
    static_missing_fit(
      panel = fx$panel, target_cells = fx$t_cells, or_blocks = NULL,
      role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
      formula = ~ x1, method = "dr", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc
    ),
    class = "cdid_static_missing_partition_error"
  )
  expect_error(
    static_missing_fit(
      panel = fx$panel, target_cells = fx$t_cells, role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
      formula = ~ x1, method = "ipw", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = NULL, accuracy = fx$acc
    ),
    class = "cdid_static_missing_control_error"
  )
  expect_error(
    static_missing_fit(
      panel = fx$panel, target_cells = fx$t_cells, role_graphs = fx$role_graphs[c("c_3_0", "c_3_1")],
      formula = ~ x1, method = "ipw", route = "baseline", weighting = "identity",
      glm_control = fx$glm_ctrl, rank_tol = 1e-10, bounds = fx$bounds, accuracy = fx$acc,
      cluster = rep("only_one_cluster", fx$panel$n)
    ),
    class = "cdid_static_graph_input_error"
  )
})

# Separate deterministic integration fixture; the saved separation failures above
# remain unchanged. No random draw or replacement of the original fixture.
build_component_integration_fixture <- function() {
  uid <- sprintf("review_%02d", 1:32)
  df <- do.call(rbind, lapply(seq_along(uid), function(i) {
    tt <- 1:6
    dd <- if (i <= 12L) as.integer(tt >= 3L) else rep.int(as.integer(i > 24L), 6L)
    data.frame(id = uid[i], time = tt, d = dd,
               y = 0.4 * tt + 0.1 * cos(i * tt) + 0.3 * sin(i) * tt + 0.8 * dd,
               x = rep(sin(i), 6L))
  }))
  panel <- static_panel_data(df, idname = "id", tname = "time", yname = "y", dname = "d")
  cells <- static_transition_cells(panel, t = 3, d = 0L, cell_id = "c")
  bounds <- gmm_bounds_control(
    max_matrix_entries = 100000L, max_batch_rows = 7L, max_targets = 10L,
    max_open_clusters = 1L, allow_materialize_if = TRUE, max_materialized_if_entries = 100000L)
  list(panel = panel, cells = cells, bounds = bounds,
       accuracy = static_graph_gmm_accuracy_control(),
       connected = list(c = list(G = data.frame(s = 2, u = 3), C = data.frame(s = 2, u = 3))),
       extra = list(c = list(
         G = list(dates = 1:6, edges = data.frame(s = c(2, 4, 5, 4), u = c(3, 5, 6, 6))),
         C = list(dates = 1:6, edges = data.frame(s = c(1, 2, 1, 4), u = c(2, 3, 3, 5))))))
}

test_that("data.frame graphs and extra components integrate without pruning or stale private fields", {
  fx <- build_component_integration_fixture()
  graph_df <- data.frame(cell_id = "c", role = c("switcher", "stayer"), s = 2, u = 3)
  for (meth in c("ipw", "reg", "dr")) {
    args <- list(panel = fx$panel, target_cells = fx$cells, method = meth, route = "baseline",
                 or_blocks = if (meth == "ipw") NULL else c(c = "P"), formula = ~ x,
                 weighting = "identity", bounds = fx$bounds, accuracy = fx$accuracy, rank_tol = 1e-10,
                 aggregate_cells = "c")
    connected <- do.call(static_missing_fit, c(args, list(role_graphs = fx$connected)))
    by_df <- do.call(static_missing_fit, c(args, list(role_graphs = graph_df)))
    extra <- do.call(static_missing_fit, c(args, list(role_graphs = fx$extra)))
    ref_args <- args[setdiff(names(args), c("bounds", "accuracy"))]
    reference <- do.call(ref_static_missing_fit, c(ref_args, list(role_graphs = fx$extra)))
    expect_equal(by_df$theta_hat, connected$theta_hat, tolerance = 1e-12)
    expect_equal(by_df$phi_theta, connected$phi_theta, tolerance = 1e-12)
    expect_length(extra$mhat, if (meth == "reg") 4L else 8L)
    expect_identical(nrow(extra$role_graphs$C$c$block_spec$edges), 4L)
    expect_identical(extra$role_graphs$G$c$block_spec$isolated_dates, 1.0)
    expect_identical(extra$role_graphs$C$c$block_spec$isolated_dates, 6.0)
    expect_equal(extra$theta_hat, reference$theta_hat, tolerance = 1e-11)
    expect_equal(extra$phi_theta, reference$phi_theta, tolerance = 1e-10)
    expect_equal(extra$cov_theta, reference$cov_theta, tolerance = 1e-11)
    expect_equal(extra$aggregate$phi, reference$aggregate$phi, tolerance = 1e-10)
    expect_equal(extra$theta_hat, connected$theta_hat, tolerance = 1e-11)
    expect_equal(extra$phi_theta, connected$phi_theta, tolerance = 1e-10)
  }
})

test_that("explicit IF materialization obeys bounds while callback-only inference remains available", {
  fx <- build_component_integration_fixture()
  args <- list(panel = fx$panel, target_cells = fx$cells, role_graphs = fx$connected,
               method = "reg", route = "baseline", or_blocks = c(c = "P"), formula = ~ x,
               weighting = "identity", accuracy = fx$accuracy, rank_tol = 1e-10)
  dense <- do.call(static_missing_fit, c(args, list(bounds = fx$bounds)))
  for (authorization in c(FALSE, TRUE)) {
    limited <- fx$bounds
    limited$allow_materialize_if <- authorization
    limited$max_materialized_if_entries <- 1L
    expect_error(do.call(static_missing_fit, c(args, list(bounds = limited, materialize_if = TRUE))),
                 class = "cdid_gmm_bound_error")
    streamed <- do.call(static_missing_fit, c(args, list(bounds = limited, materialize_if = FALSE)))
    expect_null(streamed$phi_theta)
    expect_equal(streamed$cell_if_callback(c(4L, 1L))$phi_theta,
                 dense$phi_theta[c(4L, 1L), , drop = FALSE], tolerance = 1e-12)
    expect_error(static_missing_aggregate(streamed, cells = "c", materialize_if = TRUE),
                 class = "cdid_gmm_bound_error")
    agg <- static_missing_aggregate(streamed, cells = "c", materialize_if = FALSE)
    dense_agg <- static_missing_aggregate(dense, cells = "c")
    expect_null(agg$phi)
    expect_null(agg$phi_components)
    expect_equal(agg$if_callback(c(4L, 1L))$phi, dense_agg$phi[c(4L, 1L)], tolerance = 1e-12)
    expect_equal(agg$variance, dense_agg$variance, tolerance = 1e-12)
  }
})

test_that("admitted link_if_reader batch reuse computes q/j IFs once per batch and low-budget calls decline reuse with identical results", {
  fx <- build_interior_date_fixture()

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

  run_fit <- function(bounds_arg, mat_if = TRUE) {
    static_missing_fit(
      panel = fx$panel,
      target_cells = fx$t_cells,
      all_cells = fx$all_cells,
      or_blocks = fx$or_blocks,
      role_graphs = fx$role_graphs,
      treatment_formula = ~ x1 + x2,
      outcome_formula = ~ x1 + x2,
      selection_formula = ~ x1,
      imputation_formula = ~ x2,
      method = "dr",
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
      materialize_if = mat_if,
      aggregate_cells = c("c_3_0", "c_3_1"),
      aggregate_signs = c(1, -1)
    )
  }

  # 1. Funded workspace budget admits batch reuse: 64 units / 15 batch_rows = 5 batches per pass
  #    across 4 replay passes (2 covariance + 2 unit IF projection passes) -> 20 calls per fit
  sel_calls <- 0L; imp_calls <- 0L
  fit_funded <- run_fit(fx$bounds, mat_if = TRUE)
  n_batches_per_pass <- as.integer(ceiling(fx$panel$n / fx$bounds$max_batch_rows))
  expect_identical(n_batches_per_pass, 5L)
  expect_identical(sel_calls, 4L * n_batches_per_pass)
  expect_identical(imp_calls, 4L * n_batches_per_pass)

  env_funded <- environment(fit_funded$link_if_callback)
  expect_true(isTRUE(env_funded$reuse_nuisance_if))
  budget_info <- env_funded$nuisance_batch_budget
  expect_identical(
    budget_info$Pq,
    length(fit_funded$nuisance$observation$selection_fit$layout$param_keys)
  )
  expect_identical(
    budget_info$Pj,
    length(fit_funded$nuisance$observation$imputation_fit$layout$param_keys)
  )
  expect_equal(
    budget_info$total_entries,
    budget_info$consumer_entries + budget_info$producer_entries + budget_info$holder_entries,
    tolerance = 1e-15
  )
  expect_lte(budget_info$total_entries, fx$bounds$max_matrix_entries)

  # Returned fit's adapter holder retains no B x P matrices
  env_ad_funded <- environment(fit_funded$nuisance$observation$with_if_batch)
  expect_false(env_ad_funded$batch_active)
  expect_null(env_ad_funded$batch_pos)
  expect_null(env_ad_funded$batch_sel_if)
  expect_null(env_ad_funded$batch_imp_if)

  # Single-batch reader replay calculates once per call and cleans up before returning
  sel_calls <- 0L; imp_calls <- 0L
  b1 <- fit_funded$link_if_callback(1:15)
  expect_identical(sel_calls, 1L)
  expect_identical(imp_calls, 1L)
  expect_false(env_ad_funded$batch_active)
  expect_null(env_ad_funded$batch_sel_if)
  expect_null(env_ad_funded$batch_imp_if)
  b2 <- fit_funded$link_if_callback(1:15)
  expect_identical(sel_calls, 2L)
  expect_identical(imp_calls, 2L)
  expect_identical(b1, b2)

  # 2. Tight workspace budget that funds the GMM solve itself but not the extra
  #    producer + holder reservation privately declines reuse and matches results
  bounds_tight <- fx$bounds
  bounds_tight$max_matrix_entries <- as.integer(floor(budget_info$consumer_entries))
  expect_lt(bounds_tight$max_matrix_entries, budget_info$total_entries)

  sel_calls <- 0L; imp_calls <- 0L
  fit_declined <- run_fit(bounds_tight, mat_if = TRUE)
  env_declined <- environment(fit_declined$link_if_callback)
  expect_false(isTRUE(env_declined$reuse_nuisance_if))
  expect_gt(sel_calls, 4L * n_batches_per_pass)
  expect_gt(imp_calls, 4L * n_batches_per_pass)

  env_ad_declined <- environment(fit_declined$nuisance$observation$with_if_batch)
  expect_false(env_ad_declined$batch_active)
  expect_null(env_ad_declined$batch_sel_if)
  expect_null(env_ad_declined$batch_imp_if)

  # Exact agreement between funded reuse, declined reuse, and independent reference
  rfit <- ref_static_missing_fit(
    panel = fx$panel,
    target_cells = fx$t_cells,
    all_cells = fx$all_cells,
    or_blocks = fx$or_blocks,
    role_graphs = fx$role_graphs,
    treatment_formula = ~ x1 + x2,
    outcome_formula = ~ x1 + x2,
    selection_formula = ~ x1,
    imputation_formula = ~ x2,
    method = "dr",
    route = "augmented",
    selection_model = "independent_period_logit",
    selection_pooling = "within_group",
    imputation_pooling = "within_group",
    imputation_scale = "per_period",
    weighting = "optimal",
    glm_control = fx$glm_ctrl,
    rank_tol = 1e-10,
    aggregate_cells = c("c_3_0", "c_3_1"),
    aggregate_signs = c(1, -1)
  )

  expect_identical(fit_funded$unit_id, fit_declined$unit_id)
  expect_identical(fit_funded$cell_ids, fit_declined$cell_ids)
  expect_equal(fit_funded$theta_hat, fit_declined$theta_hat, tolerance = 1e-14)
  expect_equal(fit_funded$tau_hat, fit_declined$tau_hat, tolerance = 1e-14)
  expect_equal(fit_funded$se_theta, fit_declined$se_theta, tolerance = 1e-14)
  expect_equal(fit_funded$cov_theta, fit_declined$cov_theta, tolerance = 1e-14)
  expect_equal(fit_funded$moment_covariance, fit_declined$moment_covariance, tolerance = 1e-14)
  expect_equal(fit_funded$phi_theta, fit_declined$phi_theta, tolerance = 1e-14)
  expect_equal(fit_funded$phi_tau, fit_declined$phi_tau, tolerance = 1e-14)
  expect_equal(fit_funded$aggregate$delta_hat, fit_declined$aggregate$delta_hat, tolerance = 1e-14)
  expect_equal(fit_funded$aggregate$se, fit_declined$aggregate$se, tolerance = 1e-14)
  expect_equal(fit_funded$aggregate$phi, fit_declined$aggregate$phi, tolerance = 1e-14)

  # Zero-outcome and single-outcome units retain finite, identical IFs
  zero_single_idx <- which(rowSums(fx$panel$O) <= 1L)
  expect_true(length(zero_single_idx) >= 6L)
  expect_true(all(is.finite(fit_funded$phi_theta[zero_single_idx, , drop = FALSE])))
  expect_equal(
    fit_funded$phi_theta[zero_single_idx, , drop = FALSE],
    fit_declined$phi_theta[zero_single_idx, , drop = FALSE],
    tolerance = 1e-14
  )
  expect_lt(max(abs(fit_funded$phi_theta - rfit$phi_theta)), 1e-10)
  expect_lt(max(abs(fit_declined$phi_theta - rfit$phi_theta)), 1e-10)

  # Genuine bound refusal when workspace is below the GMM solver's own requirement
  bounds_unfunded_gmm <- fx$bounds
  bounds_unfunded_gmm$max_matrix_entries <- 100L
  expect_error(run_fit(bounds_unfunded_gmm, mat_if = FALSE), class = "cdid_gmm_bound_error")
})
