# ==============================================================================
# Tests for Public Static Entry / Exit / Re-Entry Frontend (`att_gt(dname = ...)`)
# ==============================================================================

library(testthat)

ref_path <- if (file.exists("tests/testthat/reference/static_public_inference_reference.R")) {
  "tests/testthat/reference/static_public_inference_reference.R"
} else {
  "reference/static_public_inference_reference.R"
}
source(ref_path)

# Fixture builder: static binary entry/exit/re-entry panel with matched complete
# and missing outcomes (including zero- and single-outcome units).
make_static_public_fixture <- function(seed = 20261104L, n_rep = 12L) {
  set.seed(seed)
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
  obs_templates <- rbind(
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 1L),
    c(1L, 1L, 1L, 0L),
    c(1L, 1L, 0L, 1L),
    c(1L, 0L, 1L, 1L),
    c(0L, 1L, 1L, 1L),
    c(1L, 0L, 1L, 0L),
    c(1L, 0L, 0L, 0L),
    c(0L, 0L, 0L, 0L)
  )
  D_mat <- path_types[rep(seq_len(nrow(path_types)), each = n_rep), , drop = FALSE]
  n <- nrow(D_mat)
  dates <- c(1, 2, 3, 4)

  rows_comp <- vector("list", n)
  rows_miss <- vector("list", n)
  for (i in seq_len(n)) {
    rep_k <- ((i - 1L) %% nrow(obs_templates)) + 1L
    o_vec <- obs_templates[rep_k, ]
    x1 <- rnorm(4, sd = 0.35)
    x2 <- runif(4, -0.4, 0.4)
    d_vec <- D_mat[i, ]
    y_vec <- 0.5 * dates + 0.4 * x1 - 0.3 * x2 + 1.2 * d_vec + rnorm(4, sd = 0.25)

    cl_id <- sprintf("cl_%02d", ((i - 1L) %% 8L) + 1L)
    uid <- sprintf("id_%03d", i)

    rows_comp[[i]] <- data.frame(
      id = uid,
      time = dates,
      y = y_vec,
      d = d_vec,
      x1 = x1,
      x2 = x2,
      cluster = cl_id,
      stringsAsFactors = FALSE
    )
    rows_miss[[i]] <- data.frame(
      id = uid,
      time = dates,
      y = ifelse(o_vec == 1L, y_vec, NA_real_),
      d = d_vec,
      x1 = x1,
      x2 = x2,
      cluster = cl_id,
      stringsAsFactors = FALSE
    )
  }

  list(
    complete = do.call(rbind, rows_comp),
    missing = do.call(rbind, rows_miss)
  )
}

# Comparison-only fit against the independently accepted complete-data kernel.
# Production calls never execute this reference fit.
compare_static_public_complete <- function(panel,
                                                    fit_obj,
                                                    avail_target_ids,
                                                    non_target_ids,
                                                    all_cand_cells,
                                                    xformla_use,
                                                    est_method,
                                                    glm_control,
                                                    rank_tol,
                                                    compute_inffunc) {
  if (!all(panel$O == 1L) || !identical(fit_obj$route, "baseline") ||
      length(avail_target_ids) < 1L || length(non_target_ids) > 0L) {
    return(list(
      checked = FALSE,
      reason = if (!all(panel$O == 1L)) {
        "panel_has_missing_outcomes"
      } else if (!identical(fit_obj$route, "baseline")) {
        "active_observation_correction_route"
      } else if (length(non_target_ids) > 0L) {
        "non_target_or_training_cells_present"
      } else {
        "no_available_target_cells"
      }
    ))
  }

  comp_cells <- vector("list", length(avail_target_ids))
  names(comp_cells) <- avail_target_ids
  or_blocks_comp <- if (!is.null(fit_obj$or_blocks)) {
    fit_obj$or_blocks[avail_target_ids]
  } else {
    stats::setNames(avail_target_ids, avail_target_ids)
  }
  for (cid in avail_target_ids) {
    cl <- all_cand_cells[[cid]]
    dy <- as.double(panel$Y[, cl$t_idx] - panel$Y[, cl$t_prev_idx])
    des <- static_extract_reference_design(
      panel = panel,
      cell = cl,
      formula = xformla_use,
      role = "outcome",
      mask = (cl$R == 1L)
    )
    X_safe <- des$X_tilde
    X_safe[cl$R == 0L, ] <- 0.0
    X_safe[, 1L] <- 1.0
    Z_safe <- des$Z
    if (ncol(Z_safe) > 0L) {
      Z_safe[cl$R == 0L, ] <- 0.0
    }
    comp_cells[[cid]] <- static_complete_cell(
      cell_id = cid,
      d = cl$d,
      sigma = cl$sigma,
      t = cl$t_num,
      unit_id = panel$unit_id,
      R = cl$R,
      G = cl$G,
      C = cl$C,
      dy = dy,
      X_tilde = X_safe,
      Z = Z_safe,
      X_or_tilde = X_safe,
      Z_or = Z_safe,
      or_block = or_blocks_comp[[cid]]
    )
  }

  comp_fit <- static_complete_fit(
    cells = comp_cells,
    or_blocks = or_blocks_comp,
    method = est_method,
    glm_control = glm_control,
    rank_tol = rank_tol,
    unit_id = panel$unit_id,
    cluster = NULL,
    aggregate_cells = avail_target_ids,
    aggregate_signs = fit_obj$sigma
  )

  diff_tau <- max(abs(as.double(fit_obj$tau_hat[avail_target_ids]) - as.double(comp_fit$tau_hat[avail_target_ids])))
  diff_theta <- max(abs(as.double(fit_obj$theta_hat[avail_target_ids]) - as.double(comp_fit$theta_hat[avail_target_ids])))
  diff_phi <- if (isTRUE(compute_inffunc) && !is.null(fit_obj$phi_tau)) {
    max(abs(as.matrix(fit_obj$phi_tau[, avail_target_ids, drop = FALSE]) -
              as.matrix(comp_fit$phi_tau[, avail_target_ids, drop = FALSE])))
  } else {
    NA_real_
  }
  diff_agg <- if (!is.null(fit_obj$aggregate) && !is.null(comp_fit$aggregate)) {
    abs(as.double(fit_obj$aggregate$delta_hat) - as.double(comp_fit$aggregate$delta_hat))
  } else {
    NA_real_
  }

  list(
    checked = TRUE,
    max_abs_diff_tau = diff_tau,
    max_abs_diff_theta = diff_theta,
    max_abs_diff_phi_tau = diff_phi,
    max_abs_diff_aggregate = diff_agg
  )
}


complete_reduction_for_result <- function(res, data, formula) {
  panel <- static_panel_data(data, "id", "time", "y", "d")
  compare_static_public_complete(
    panel = panel, fit_obj = res$fit, avail_target_ids = res$target_id,
    non_target_ids = res$request$non_target_or_cells,
    all_cand_cells = res$fit$cells, xformla_use = formula,
    est_method = res$request$est_method,
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    rank_tol = 1e-7, compute_inffunc = TRUE
  )
}

test_that("independent static public inference reference has zero production symbols and att_gt signature is exact", {
  r_dir <- if (dir.exists("R")) "R" else "../../R"
  ref_env <- new.env(parent = baseenv())
  sys.source(ref_path, envir = ref_env)
  prod_env <- new.env(parent = asNamespace("stats"))
  for (f in sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))) {
    sys.source(f, envir = prod_env)
  }
  prod_syms <- ls(prod_env)
  ref_funs <- ls(ref_env)
  used_globals <- unique(unlist(lapply(ref_funs, function(fn) {
    if (is.function(ref_env[[fn]])) {
      codetools::findGlobals(ref_env[[fn]], merge = FALSE)$functions
    } else {
      character(0)
    }
  })))
  expect_identical(intersect(used_globals, prod_syms), character(0))

  f_ours <- as.list(formals(att_gt))
  n_up <- length(upstream_att_gt_formals) + 1L # 24 upstream formals + `...`
  expect_identical(names(f_ours)[seq_along(upstream_att_gt_formals)], upstream_att_gt_formals)
  expect_identical(names(f_ours)[n_up], "...")
  expect_identical(names(f_ours)[-seq_len(n_up)], cdid_chained_only_formals)
  expect_null(f_ours[["allow_unbalanced_panel"]])
  expect_identical(f_ours[["est_method"]], "dr")
  expect_identical(f_ours[["base_period"]], "varying")
  expect_identical(
    tail(names(f_ours), 5L),
    c("outcome_pooling", "dname", "imputation_scale", "switchers", "switch_times")
  )
  expect_identical(f_ours[["outcome_pooling"]], "separate")
  expect_null(f_ours[["dname"]])
  expect_identical(f_ours[["imputation_scale"]], "per_period")
  expect_identical(f_ours[["switchers"]], "both")
  expect_null(f_ours[["switch_times"]])

  if (requireNamespace("did", quietly = TRUE)) {
    f_did <- as.list(formals(did::att_gt))
    expect_identical(names(f_ours)[seq_len(n_up)], names(f_did))
    same_up <- setdiff(names(f_did), "allow_unbalanced_panel")
    expect_identical(f_ours[same_up], f_did[same_up])
  }
})

test_that("routing and lazy promise evaluation preserve ordinary and cohort gname paths and inactive static options", {
  set.seed(101)
  n_u <- 24L
  g_u <- rep(c(0, 2), each = n_u / 2L)
  df_coh <- data.frame(
    id = rep(seq_len(n_u), each = 2L),
    period = rep(c(1, 2), times = n_u),
    G = rep(g_u, each = 2L),
    x = rep(rnorm(n_u), each = 2L),
    y = rnorm(2L * n_u)
  )
  df_coh$d <- as.integer(df_coh$G == 2 & df_coh$period == 2)

  # 1. Ordinary mode never forces dname or any static/chained extension formal
  if (requireNamespace("did", quietly = TRUE)) {
    expect_warning(ord_res <- att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      data = df_coh, mode = "ordinary", bstrap = FALSE,
      dname = stop("dname forced in ordinary mode"),
      outcome_pooling = stop("outcome_pooling forced in ordinary mode"),
      imputation_scale = stop("imputation_scale forced in ordinary mode"),
      switchers = stop("switchers forced in ordinary mode"),
      switch_times = stop("switch_times forced in ordinary mode")
    ), "No pre-treatment periods available for the Wald pre-test")
    expect_s3_class(ord_res, "MP")
  } else {
    # Without `did` installed in the source harness, `forward_call(quote(did::att_gt), ...)`
    # attempts to loadNamespace("did") without forcing any static/chained promise.
    expect_error(
      att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        data = df_coh, mode = "ordinary", bstrap = FALSE,
        dname = stop("dname forced in ordinary mode"),
        outcome_pooling = stop("outcome_pooling forced in ordinary mode"),
        imputation_scale = stop("imputation_scale forced in ordinary mode"),
        switchers = stop("switchers forced in ordinary mode"),
        switch_times = stop("switch_times forced in ordinary mode")
      ),
      class = "packageNotFoundError"
    )
  }

  # 2. Chained cohort gname path (dname = NULL) never forces static-only formals
  #    and routes shared outcome_pooling ("separate" default, explicit "pooled",
  #    lazy under est_method = "ipw", and refusing static-only partitions on gname)
  for (op_val in c("separate", "pooled")) {
    expect_error(
      att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        xformla = ~ x, data = stop("reached chained_public_att_gt without forcing static formals"),
        bstrap = FALSE, cband = FALSE,
        outcome_pooling = op_val,
        imputation_scale = stop("imputation_scale forced on gname route"),
        switchers = stop("switchers forced on gname route"),
        switch_times = stop("switch_times forced on gname route")
      ),
      "reached chained_public_att_gt without forcing static formals"
    )
  }
  expect_error(
    att_gt(
      yname = "y", tname = "period", idname = "id", gname = "G",
      xformla = ~ x, data = stop("reached chained_public_att_gt without forcing static formals"),
      est_method = "ipw", bstrap = FALSE, cband = FALSE,
      outcome_pooling = stop("outcome_pooling forced under gname est_method = ipw"),
      imputation_scale = stop("imputation_scale forced on gname route"),
      switchers = stop("switchers forced on gname route"),
      switch_times = stop("switch_times forced on gname route")
    ),
    "reached chained_public_att_gt without forcing static formals"
  )
  for (bad_gname_pool in list(c(t2_entry = "b1", t3_entry = "b1"), list(b1 = c("t2_entry", "t3_entry")))) {
    e_pool <- expect_error(
      att_gt(
        yname = "y", tname = "period", idname = "id", gname = "G",
        xformla = ~ x, data = stop("data must not be forced on gname outcome_pooling refusal"),
        bstrap = FALSE, cband = FALSE,
        outcome_pooling = bad_gname_pool,
        imputation_scale = stop("imputation_scale forced on gname route"),
        switchers = stop("switchers forced on gname route"),
        switch_times = stop("switch_times forced on gname route")
      ),
      class = "cdid_chained_request_outcome_pooling"
    )
    expect_no_match(conditionMessage(e_pool), "data must not be forced")
  }

  # 3. Conflicting explicit non-NULL gname and dname is diagnosed before data is touched
  expect_error(
    att_gt(
      yname = "y", tname = "period", idname = "id",
      gname = "G", dname = "d",
      data = stop("data must not be forced on gname/dname conflict")
    ),
    class = "cdid_static_routing_error"
  )

  # 4. Missing both gname and dname in chained mode is diagnosed before data is touched
  expect_error(
    att_gt(
      yname = "y", tname = "period", idname = "id",
      data = stop("data must not be forced when both gname and dname are omitted")
    ),
    class = "cdid_static_routing_error"
  )

  # 5. Static dname route with omitted gname (or gname = NULL) never forces inactive nuisance options/formulas
  st_base_ipw <- att_gt(
    yname = "y", tname = "period", idname = "id", dname = "d",
    xformla = ~ x, data = df_coh, est_method = "ipw",
    selection = "none",
    selection_model = stop("selection_model forced under selection = none"),
    selection_pooling = stop("selection_pooling forced under selection = none"),
    selection_xformla = stop("selection_xformla forced under selection = none"),
    imputation_pooling = stop("imputation_pooling forced under selection = none"),
    imputation_scale = stop("imputation_scale forced under selection = none"),
    outcome_pooling = stop("outcome_pooling forced under est_method = ipw"),
    bstrap = FALSE, cband = FALSE
  )
  expect_s3_class(st_base_ipw, "cdid_static_att_gt")

  fix <- make_static_public_fixture()
  # Weighting route never forces imputation options
  st_wt <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = fix$missing,
    est_method = "ipw", selection = "conditional_mtar", selection_est_method = "ipw",
    imputation_pooling = stop("imputation_pooling forced under weighting route"),
    imputation_scale = stop("imputation_scale forced under weighting route"),
    outcome_pooling = stop("outcome_pooling forced under est_method = ipw"),
    bstrap = FALSE, cband = FALSE
  )
  expect_s3_class(st_wt, "cdid_static_att_gt")

  # Imputation route never forces selection model/pooling/formula
  st_imp <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = fix$missing,
    est_method = "reg", selection = "conditional_mtar", selection_est_method = "imputation",
    selection_model = stop("selection_model forced under imputation route"),
    selection_pooling = stop("selection_pooling forced under imputation route"),
    selection_xformla = stop("selection_xformla forced under imputation route"),
    bstrap = FALSE, cband = FALSE
  )
  expect_s3_class(st_imp, "cdid_static_att_gt")

  # Augmented route with est_method = "ipw" never forces outcome_pooling
  st_aug_ipw <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, selection_xformla = ~ x2, data = fix$missing,
    est_method = "ipw", selection = "conditional_mtar", selection_est_method = "dr",
    outcome_pooling = stop("outcome_pooling forced under augmented ipw"),
    bstrap = FALSE, cband = FALSE
  )
  expect_s3_class(st_aug_ipw, "cdid_static_att_gt")
})

test_that("static scope gate refuses inapplicable cohort, weight, frame, anticipation, and PT options before forcing data", {
  bad_data <- function() stop("data promise must not be forced")

  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), panel = FALSE),
    class = "cdid_chained_panel"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), allow_unbalanced_panel = FALSE),
    class = "cdid_chained_balanced"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), est_method = "stdipw"),
    class = "cdid_chained_unavailable"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), est_method = function(...) NULL),
    class = "cdid_static_est_method_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), extra_arg = 1),
    class = "cdid_static_dots_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), weightsname = "w"),
    class = "cdid_static_weights_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), fix_weights = "first_period"),
    class = "cdid_static_weights_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), frame = data.frame(id = 1)),
    class = "cdid_static_frame_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), anticipation = 1),
    class = "cdid_static_anticipation_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), pt_assumption = "all", pt_start = 1),
    class = "cdid_static_pt_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), control_group = "notyettreated"),
    class = "cdid_static_cohort_option_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), base_period = "universal"),
    class = "cdid_static_cohort_option_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), selection = "conditional_mtar", selection_universe = "A"),
    class = "cdid_static_selection_error"
  )
  expect_error(
    att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
           data = bad_data(), selection = "none", selection_est_method = "dr"),
    class = "cdid_static_selection_error"
  )

  # Positive support for augmented observation correction (selection_est_method = "dr")
  # across all three treatment methods: est_method in {"dr", "ipw", "reg"}
  fix_gate <- make_static_public_fixture()
  for (meth in c("dr", "ipw", "reg")) {
    res_aug_m <- att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1, selection_xformla = ~ x2,
      data = fix_gate$missing,
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
    expect_s3_class(res_aug_m, "cdid_static_att_gt")
    expect_identical(res_aug_m$request$route, "augmented")
    expect_identical(res_aug_m$request$est_method, meth)
    expect_identical(res_aug_m$request$selection_est_method, "dr")
    expect_true(all(is.finite(res_aug_m$att)))
    expect_true(is.finite(res_aug_m$aggregate$estimate))
  }

  # The new static interface exposes its approved names, without cohort or
  # draft convenience aliases. Existing gname aliases remain unchanged.
  for (alias in c("within_cohort", "across_cohorts", "by_group")) {
    expect_error(static_public_normalize_pooling(alias, "selection_pooling"),
                 class = "cdid_static_pooling_error")
    expect_error(static_public_normalize_pooling(alias, "imputation_pooling"),
                 class = "cdid_static_pooling_error")
  }
  for (alias in c("or", "OR")) {
    expect_error(att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
                        data = bad_data(), est_method = alias),
                 class = "cdid_static_est_method_error")
  }
  for (alias in c("weighting", "augmented")) {
    expect_error(att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
                        data = bad_data(), selection = "conditional_mtar",
                        selection_est_method = alias),
                 class = "cdid_static_selection_error")
  }
  expect_error(att_gt(yname = "y", tname = "t", idname = "id", dname = "d",
                      data = bad_data(), selection = "corollary1"),
               class = "cdid_static_selection_error")

  # Unknown D in data is refused by static_panel_data
  fix <- make_static_public_fixture()
  df_bad_d <- fix$complete
  df_bad_d$d[1L] <- NA_integer_
  expect_error(
    att_gt(yname = "y", tname = "time", idname = "id", dname = "d",
           data = df_bad_d, bstrap = FALSE, cband = FALSE),
    class = "cdid_static_error"
  )

  # aggte() refuses cdid_static_att_gt objects (cohort-only)
  st_ok <- att_gt(yname = "y", tname = "time", idname = "id", dname = "d",
                  xformla = ~ x1, data = fix$complete, bstrap = FALSE, cband = FALSE)
  expect_error(aggte(st_ok, type = "simple"))
})

test_that("complete-data all-pair static route reduces to static_complete_fit and matches 2-period absorbing entry", {
  fix <- make_static_public_fixture()
  df_c <- fix$complete

  for (meth in c("ipw", "reg", "dr")) {
    for (pool in list("separate", "pooled",
                      c(t2_entry = "entry", t2_exit = "exit", t3_entry = "entry",
                        t3_exit = "exit", t4_entry = "entry"))) {
      res <- att_gt(
        yname = "y", tname = "time", idname = "id", dname = "d",
        xformla = ~ x1 + x2, data = df_c,
        est_method = meth, outcome_pooling = pool,
        bstrap = FALSE, cband = FALSE
      )
      red <- complete_reduction_for_result(res, df_c, ~ x1 + x2)
      expect_true(isTRUE(red$checked))
      expect_lt(red$max_abs_diff_tau, 1e-11)
      expect_lt(red$max_abs_diff_theta, 1e-11)
      expect_lt(red$max_abs_diff_phi_tau, 1e-11)
      expect_lt(red$max_abs_diff_aggregate, 1e-11)
    }
    # Single-cell complete-data reduction on the 4-period panel (all 6 pairs retained)
    res_single <- att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1 + x2, data = df_c,
      est_method = meth, switchers = "entry", switch_times = 2,
      bstrap = FALSE, cband = FALSE
    )
    red_s <- complete_reduction_for_result(res_single, df_c, ~ x1 + x2)
    expect_true(isTRUE(red_s$checked))
    expect_lt(red_s$max_abs_diff_tau, 1e-11)
    expect_lt(red_s$max_abs_diff_theta, 1e-11)
    expect_lt(red_s$max_abs_diff_phi_tau, 1e-11)
    expect_lt(red_s$max_abs_diff_aggregate, 1e-11)
  }

  # 2-period absorbing entry equivalence with chained cohort pipeline (and att_gt(gname=...) when did is present)
  set.seed(303)
  n_u <- 40L
  g_u <- rep(c(0, 2), each = n_u / 2L)
  x1_u <- rnorm(n_u, sd = 0.5)
  df_2p <- data.frame(
    id = rep(seq_len(n_u), each = 2L),
    time = rep(c(1, 2), times = n_u),
    G = rep(g_u, each = 2L),
    x1 = rep(x1_u, each = 2L),
    y = rnorm(2L * n_u)
  )
  df_2p$d <- as.integer(df_2p$G == 2 & df_2p$time == 2)

  # Compare 2-period static IPW/OR/DR against direct 2-period closed-form / complete-core
  for (meth in c("ipw", "reg", "dr")) {
    r_d <- att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1, data = df_2p, est_method = meth,
      bstrap = FALSE, cband = FALSE
    )
    red_2p <- complete_reduction_for_result(r_d, df_2p, ~ x1)
    expect_true(isTRUE(red_2p$checked))
    expect_lt(red_2p$max_abs_diff_tau, 1e-12)
    expect_lt(red_2p$max_abs_diff_phi_tau, 1e-12)

    if (requireNamespace("did", quietly = TRUE)) {
      r_g <- att_gt(
        yname = "y", tname = "time", idname = "id", gname = "G",
        xformla = ~ x1, data = df_2p, est_method = meth,
        bstrap = FALSE, cband = FALSE
      )
      expect_equal(unname(r_d$att), unname(r_g$att), tolerance = 1e-11)
      expect_equal(unname(r_d$se), unname(r_g$se), tolerance = 1e-11)
      expect_equal(as.double(r_d$inffunc[, 1L]), as.double(r_g$inffunc[, 1L]), tolerance = 1e-11)
    }
  }
})

test_that("target discovery, explicit empty cells, disconnected extra graph components, and missing-outcome routes obey contracts", {
  fix <- make_static_public_fixture()
  df_m <- fix$missing

  # 1. Automatic target discovery includes the 5 nonempty cells:
  #    t2_entry, t2_exit, t3_entry, t3_exit, t4_entry.
  res_auto <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1 + x2, data = df_m,
    bstrap = FALSE, cband = FALSE
  )
  expect_identical(
    res_auto$target_id,
    c("t2_entry", "t2_exit", "t3_entry", "t3_exit", "t4_entry")
  )
  expect_true(all(res_auto$available))
  expect_true(isTRUE(res_auto$aggregate$available))

  # Direction filter: switchers = "entry" and switchers = "exit"
  res_entry <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m, switchers = "entry",
    bstrap = FALSE, cband = FALSE
  )
  expect_identical(res_entry$target_id, c("t2_entry", "t3_entry", "t4_entry"))
  expect_true(all(res_entry$d == 0L))
  expect_true(all(res_entry$sigma == 1))

  res_exit <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m, switchers = "exit",
    bstrap = FALSE, cband = FALSE
  )
  expect_identical(res_exit$target_id, c("t2_exit", "t3_exit"))
  expect_true(all(res_exit$d == 1L))
  expect_true(all(res_exit$sigma == -1))

  # 2. Explicit request for empty cell (t = 4, exit): retained as available = FALSE,
  #    band and aggregate marked unavailable (never dropped or renormalized).
  res_explicit_empty <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m,
    switchers = "exit", switch_times = c(3, 4),
    bstrap = TRUE, cband = TRUE, biters = 100
  )
  expect_identical(res_explicit_empty$target_id, c("t3_exit", "t4_exit"))
  expect_identical(unname(res_explicit_empty$available), c(TRUE, FALSE))
  expect_identical(unname(res_explicit_empty$support_status), c("ok", "empty_switchers"))
  expect_true(is.finite(res_explicit_empty$att[["t3_exit"]]))
  expect_true(is.na(res_explicit_empty$att[["t4_exit"]]))
  expect_true(is.na(res_explicit_empty$cband_crit_val))
  expect_false(isTRUE(res_explicit_empty$aggregate$available))
  expect_true(is.na(res_explicit_empty$aggregate$estimate))

  pan <- static_panel_data(df_m, "id", "time", "y", "d")
  cand <- static_transition_cells(pan,
    t = rep(pan$dates_num[-1L], each = 2L), d = rep(0:1, length(pan$dates_num) - 1L),
    cell_id = as.vector(rbind(paste0("t", pan$date_labels[-1L], "_entry"),
                             paste0("t", pan$date_labels[-1L], "_exit"))))
  for (alias in c("by_direction", "within_group", "Within_group")) {
    expect_error(static_public_resolve_outcome_pooling(alias, res_auto$target_id, cand$cells),
                 class = "cdid_static_outcome_pooling_error")
  }
  expect_error(static_public_resolve_targets(pan, cand, "both", "t4_entry"),
               class = "cdid_static_switch_times_error")
  expect_identical(static_public_resolve_targets(pan, cand, "both", "4")$target_ids,
                   c("t4_entry", "t4_exit"))
  res_date_default <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m, switch_times = 4,
    bstrap = FALSE, cband = FALSE
  )
  res_date_both <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m, switch_times = 4, switchers = "both",
    bstrap = FALSE, cband = FALSE
  )
  expect_identical(res_date_default$target_id, c("t4_entry", "t4_exit"))
  expect_identical(res_date_default$target_id, res_date_both$target_id)
  expect_identical(res_date_default$att, res_date_both$att)
  expect_identical(res_date_default$available, res_date_both$available)
  expect_identical(unname(res_date_default$available), c(TRUE, FALSE))
  expect_false(res_date_default$aggregate$available)

  # 3. Disconnected extra component in role graph is retained (K_c^a = 2 > 1)
  df_5 <- rbind(
    fix$complete,
    transform(fix$complete[fix$complete$time == 4, ], time = 5)
  )
  df_5 <- df_5[order(df_5$id, df_5$time), ]
  uid_num <- as.integer(sub("^id_", "", df_5$id))
  in_first_half <- ((uid_num - 1L) %% 12L) < 6L
  keep_y <- (in_first_half & (df_5$time %in% c(1, 2))) |
    (!in_first_half & (df_5$time %in% c(4, 5)))
  df_5$y <- ifelse(keep_y, 1.0 + 0.2 * df_5$x1 + 0.5 * df_5$d, NA_real_)
  res_disc <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_5,
    switchers = "entry", switch_times = 2,
    bstrap = FALSE, cband = FALSE
  )
  expect_true(res_disc$available[["t2_entry"]])
  expect_identical(res_disc$graphs$by_cell_role[["t2_entry::G"]]$n_components, 2L)
  expect_identical(res_disc$graphs$by_cell_role[["t2_entry::C"]]$n_components, 2L)
  expect_true(is.finite(res_disc$att[["t2_entry"]]))

  # 4. Missing-outcome routes, poolings, custom OR partition with non-target training cell,
  #    and distinct selection_xformla match static_missing_fit
  or_part <- c(t3_entry = "b_entry", t3_exit = "b_exit", t4_entry = "b_entry")
  res_aug_custom <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, selection_xformla = ~ x2,
    data = df_m,
    est_method = "dr",
    selection = "conditional_mtar",
    selection_est_method = "dr",
    selection_model = "independent_period_logit",
    selection_pooling = "within_group",
    imputation_pooling = "pooled",
    imputation_scale = "per_change",
    outcome_pooling = or_part,
    switch_times = 3, switchers = "both",
    bstrap = FALSE, cband = FALSE
  )
  expect_identical(res_aug_custom$target_id, c("t3_entry", "t3_exit"))
  expect_identical(res_aug_custom$request$non_target_or_cells, "t4_entry")
  expect_false(res_aug_custom$request$selection_formula_inherited)
  expect_equal(unname(res_aug_custom$att), unname(res_aug_custom$fit$tau_hat), tolerance = 1e-14)
  expect_equal(res_aug_custom$aggregate$estimate, res_aug_custom$fit$aggregate$delta_hat, tolerance = 1e-14)

  # Identity GMM weighting option
  res_id <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m,
    est_method = "ipw", gmm_weighting = "identity",
    switch_times = 3, switchers = "entry",
    bstrap = FALSE, cband = FALSE
  )
  expect_identical(res_id$request$gmm_weighting$canonical, "identity")

  # Opt-in counted_work_limit in control refuses tiny budget before fitting
  expect_error(
    att_gt(
      yname = "y", tname = "time", idname = "id", dname = "d",
      xformla = ~ x1, data = df_m,
      control = list(counted_work_limit = 10),
      bstrap = FALSE, cband = FALSE
    ),
    class = "cdid_static_work_limit_error"
  )
})

test_that("deterministic inference checks match independent reference and obey clustering, degenerate-family, and compute_inffunc contracts", {
  fix <- make_static_public_fixture()
  df_m <- fix$missing

  # Fit unit-level and clustered static models
  res_u <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m,
    est_method = "dr", switch_times = c(2, 3), switchers = "both",
    bstrap = FALSE, cband = FALSE
  )
  res_cl <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m, clustervars = "cluster",
    est_method = "dr", switch_times = c(2, 3), switchers = "both",
    bstrap = FALSE, cband = FALSE
  )

  # Clustering changes inference only: point estimates, IFs, and aggregate point estimate are identical
  expect_identical(res_u$att, res_cl$att)
  expect_identical(res_u$theta_hat, res_cl$theta_hat)
  expect_identical(res_u$inffunc, res_cl$inffunc)
  expect_identical(res_u$aggregate$estimate, res_cl$aggregate$estimate)
  expect_false(isTRUE(all.equal(res_u$V_analytical, res_cl$V_analytical)))

  expect_identical(res_u$n, length(unique(df_m$id)))
  expect_identical(as.character(res_u$unit_id), sort(unique(df_m$id)))
  expect_identical(rownames(res_u$inffunc), as.character(res_u$unit_id))
  expect_identical(names(res_u$aggregate$inffunc), as.character(res_u$unit_id))
  expect_gt(res_u$diagnostics$n_zero_outcome_units, 0L)
  expect_gt(res_u$diagnostics$n_single_outcome_units, 0L)
  expect_true(all(is.finite(res_u$inffunc)))

  # Verify analytical covariance and aggregate IF against independent reference
  G_mat <- do.call(cbind, lapply(res_u$fit$cells[res_u$target_id], `[[`, "G"))
  ref_agg <- ref_static_public_aggregate_if(
    theta_hat = res_u$theta_hat,
    phi_theta = res_u$fit$phi_theta[, res_u$target_id, drop = FALSE],
    G_mat = G_mat,
    signs = res_u$sigma
  )
  expect_equal(res_u$aggregate$estimate, ref_agg$delta_hat, tolerance = 1e-13)
  expect_equal(as.double(res_u$aggregate$inffunc), as.double(ref_agg$phi), tolerance = 1e-13)

  # Fixed-sign Rademacher multiplier check (unit-level and clustered)
  B <- 200L
  set.seed(999L)
  signs_u <- matrix(sample(c(-1.0, 1.0), B * res_u$n, replace = TRUE), nrow = B, ncol = res_u$n)
  inf_u <- static_public_inference(
    estimates = res_u$att,
    if_matrix = res_u$inffunc,
    unit_ids = res_u$unit_id,
    cluster_vec = NULL,
    alp = 0.05,
    bstrap = TRUE,
    cband = TRUE,
    biters = B,
    chunk_rows = 7L,
    fixed_signs = signs_u,
    agg_estimate = res_u$aggregate$estimate,
    agg_if = res_u$aggregate$inffunc,
    all_declared_available = TRUE
  )
  ref_u <- ref_static_public_cluster_inference(
    phi_mat = res_u$inffunc,
    cluster_vec = NULL,
    signs_mat = signs_u,
    alp = 0.05
  )
  ref_u_agg <- ref_static_public_cluster_inference(
    phi_mat = matrix(res_u$aggregate$inffunc, ncol = 1L, dimnames = list(NULL, "__aggregate__")),
    cluster_vec = NULL,
    signs_mat = signs_u,
    alp = 0.05
  )
  expect_equal(inf_u$V_analytical, ref_u$V_analytical, tolerance = 1e-13)
  expect_equal(inf_u$analytical_se, ref_u$analytical_se, tolerance = 1e-13)
  expect_equal(inf_u$bootstrap_iqr_scale, ref_u$iqr_scale, tolerance = 1e-13)
  expect_equal(inf_u$bootstrap_se, ref_u$bootstrap_se, tolerance = 1e-13)
  expect_equal(inf_u$cband_crit_val, ref_u$crit_val, tolerance = 1e-13)
  expect_equal(inf_u$agg_analytical_se, unname(ref_u_agg$analytical_se), tolerance = 1e-13)
  expect_equal(inf_u$agg_bootstrap_se, unname(ref_u_agg$bootstrap_se), tolerance = 1e-13)

  # Clustered fixed-sign check
  cl_info <- static_public_extract_clusters(df_m, static_panel_data(df_m, "id", "time", "y", "d"), "id", "cluster")
  G_cl <- cl_info$n_clusters
  signs_cl <- matrix(sample(c(-1.0, 1.0), B * G_cl, replace = TRUE), nrow = B, ncol = G_cl)
  inf_cl <- static_public_inference(
    estimates = res_cl$att,
    if_matrix = res_cl$inffunc,
    unit_ids = res_cl$unit_id,
    cluster_vec = cl_info$cluster_vec,
    alp = 0.05,
    bstrap = TRUE,
    cband = TRUE,
    biters = B,
    chunk_rows = 7L,
    fixed_signs = signs_cl,
    agg_estimate = res_cl$aggregate$estimate,
    agg_if = res_cl$aggregate$inffunc,
    all_declared_available = TRUE
  )
  ref_cl <- ref_static_public_cluster_inference(
    phi_mat = res_cl$inffunc,
    cluster_vec = cl_info$cluster_vec,
    signs_mat = signs_cl,
    alp = 0.05
  )
  expect_equal(inf_cl$V_analytical, ref_cl$V_analytical, tolerance = 1e-13)
  expect_equal(inf_cl$analytical_se, ref_cl$analytical_se, tolerance = 1e-13)
  expect_equal(inf_cl$bootstrap_se, ref_cl$bootstrap_se, tolerance = 1e-13)
  expect_equal(inf_cl$cband_crit_val, ref_cl$crit_val, tolerance = 1e-13)
  ref_cl_agg <- ref_static_public_cluster_inference(
    phi_mat = matrix(res_cl$aggregate$inffunc, ncol = 1L),
    cluster_vec = cl_info$cluster_vec, signs_mat = signs_cl, alp = 0.05
  )
  expect_equal(inf_cl$agg_analytical_se, unname(ref_cl_agg$analytical_se), tolerance = 1e-13)
  expect_equal(inf_cl$agg_bootstrap_se, unname(ref_cl_agg$bootstrap_se), tolerance = 1e-13)

  # Generated signs use the caller stream once, in bounded batches. Compare
  # both cells and aggregate to independently retained signs for that stream.
  set.seed(117L)
  generated_signs <- matrix(NA_real_, B, G_cl)
  for (start in seq.int(1L, B, by = 7L)) {
    ii <- seq.int(start, min(start + 6L, B))
    generated_signs[ii, ] <- matrix(sample(c(-1.0, 1.0), length(ii) * G_cl,
                                         replace = TRUE), length(ii), G_cl)
  }
  expected_rng <- .Random.seed
  ref_generated <- ref_static_public_cluster_inference(
    res_cl$inffunc, cl_info$cluster_vec, generated_signs)
  ref_generated_agg <- ref_static_public_cluster_inference(
    matrix(res_cl$aggregate$inffunc, ncol = 1L), cl_info$cluster_vec, generated_signs)
  set.seed(117L)
  inf_generated <- static_public_inference(
    res_cl$att, res_cl$inffunc, res_cl$unit_id, cl_info$cluster_vec,
    biters = B, chunk_rows = 7L, agg_if = res_cl$aggregate$inffunc
  )
  expect_identical(.Random.seed, expected_rng)
  expect_equal(inf_generated$bootstrap_se, ref_generated$bootstrap_se, tolerance = 1e-13)
  expect_equal(inf_generated$cband_crit_val, ref_generated$crit_val, tolerance = 1e-13)
  expect_equal(inf_generated$agg_bootstrap_se, unname(ref_generated_agg$bootstrap_se), tolerance = 1e-13)
  expect_error(static_public_inference(
    res_u$att, res_u$inffunc, res_u$unit_id, biters = B,
    workspace_entries = 1, workspace_bytes = 8),
    class = "cdid_static_inference_bound_error")

  # Familiar att_gt FALSE/TRUE behavior: analytic pointwise inference, no draws.
  set.seed(118L)
  before_rng <- .Random.seed
  inf_false_true <- static_public_inference(
    res_u$att, res_u$inffunc, res_u$unit_id, bstrap = FALSE, cband = TRUE,
    biters = B, agg_if = res_u$aggregate$inffunc)
  expect_identical(.Random.seed, before_rng)
  expect_null(inf_false_true$multiplier)
  expect_identical(inf_false_true$se, inf_u$analytical_se)
  expect_identical(inf_false_true$cband_status, "not_requested")

  # Degenerate zero-scale target in declared family makes band unavailable (no survival subfamily)
  phi_deg <- res_u$inffunc
  phi_deg[, 2L] <- 0.0
  inf_deg <- static_public_inference(
    estimates = res_u$att,
    if_matrix = phi_deg,
    unit_ids = res_u$unit_id,
    cluster_vec = NULL,
    alp = 0.05,
    bstrap = TRUE,
    cband = TRUE,
    biters = B,
    fixed_signs = signs_u,
    all_declared_available = TRUE
  )
  expect_true(is.na(inf_deg$cband_crit_val))
  expect_match(inf_deg$cband_status, "^unavailable_degenerate_family")
  expect_true(is.na(inf_deg$bootstrap_se[2L]))
  expect_true(all(is.finite(inf_deg$bootstrap_se[-2L])))

  # Point-only compute_inffunc = FALSE
  before_point_rng <- .Random.seed
  res_pt <- att_gt(
    yname = "y", tname = "time", idname = "id", dname = "d",
    xformla = ~ x1, data = df_m,
    est_method = "dr", switch_times = c(2, 3), switchers = "both",
    compute_inffunc = FALSE
  )
  expect_identical(.Random.seed, before_point_rng)
  expect_equal(res_pt$att, res_u$att, tolerance = 1e-13)
  expect_equal(res_pt$aggregate$estimate, res_u$aggregate$estimate, tolerance = 1e-13)
  expect_null(res_pt$inffunc)
  expect_null(res_pt$V_analytical)
  expect_true(all(is.na(res_pt$se)))
  expect_false(res_pt$bstrap)
  expect_false(res_pt$cband)

  # S3 methods: tidy, glance, summary, print, plot/ggdid
  td <- tidy(res_cl)
  gl <- glance(res_cl)
  sm <- summary(res_cl)
  out_txt <- capture.output(print(res_cl))
  expect_s3_class(td, "data.frame")
  expect_identical(nrow(td), length(res_cl$target_id))
  expect_s3_class(gl, "data.frame")
  expect_identical(nrow(gl), 1L)
  expect_s3_class(sm, "summary.cdid_static_att_gt")
  expect_true(any(grepl("Static No-Carryover Entry/Exit/Re-Entry DiD", out_txt, fixed = TRUE)))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_obj <- ggdid(res_cl)
    expect_s3_class(p_obj, "gg")
  }
})
