# ==============================================================================
# Public Static Entry / Exit / Re-Entry Frontend (`att_gt(dname = ...)`)
# ==============================================================================
#
# Implements the user-facing static binary panel route of `att_gt()` when
# `dname` is supplied, delegating to the accepted private static kernels:
#   - `static_panel_data()`, `static_transition_cells()`,
#     `static_observation_pair_accessor()`, `static_extract_reference_design()`
#     (`R/static_panel_data.R`)
#   - `static_complete_cell()`, `static_complete_fit()`,
#     `static_complete_aggregate()` (`R/static_complete_core.R`)
#   - `static_role_block_spec()`, `static_graph_gmm_accuracy_control()`,
#     `solve_static_graph_gmm()` (`R/static_graph_gmm.R`)
#   - `static_missing_fit()`, `static_missing_aggregate()`
#     (`R/static_missing_fit.R`)
#   - `chained_public_resolve_control()` (`R/chained_public_control.R`)
#   - `chained_cluster_multiplier()` (`R/chained_multiplier.R`)
#
# Scope:
#   - Source-static no-carryover identification benchmark (`Y_{it}(d_{1:t}) = Y_{it}(d_t)`).
#   - Known binary treatment history `D = (D_1, ..., D_T)` across all study dates
#     for every included unit; missing outcomes `Y` allowed, including zero- and
#     single-outcome units.
#   - All-positive-count observation pairs (`n_observed_role >= 1`) per required
#     switcher/stayer role population, retaining all connected components
#     (`K_c^a >= 1`) without model-failure pruning or path selection.
#   - Sign-oriented cell treatment effects `tau_c = sigma_c * theta_c` with
#     `sigma_c = 1 - 2 d`, and switcher-share aggregate `delta_{A, s}` with
#     random switch-share influence function.
#   - Original-unit analytical and common-sign Rademacher multiplier bootstrap
#     inference via `chained_cluster_multiplier()`.
# ==============================================================================

static_public_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    list(message = message, call = NULL, diagnostics = diagnostics),
    class = unique(c(
      class,
      "cdid_static_public_error",
      "cdid_static_error",
      "cdid_chained_request_error",
      "cdid_error",
      "error",
      "condition"
    ))
  ))
}

static_public_need <- function(ok, message, class = character(), diagnostics = list()) {
  if (!isTRUE(ok)) {
    static_public_abort(message = message, class = class, diagnostics = diagnostics)
  }
  invisible(TRUE)
}

static_public_normalize_pooling <- function(val, arg_name) {
  static_public_need(
    chained_request_chr1(val),
    sprintf("`%s` must be a single non-empty character string.", arg_name),
    c("cdid_static_pooling_error", "cdid_chained_request_input_error")
  )
  if (val %in% c("separate", "within_group", "pooled")) {
    return(val)
  }
  static_public_abort(
    sprintf(
      "`%s` must be one of `\"separate\"`, `\"within_group\"`, or `\"pooled\"` (got `%s`).",
      arg_name, val
    ),
    c("cdid_static_pooling_error", "cdid_chained_request_input_error")
  )
}

# Validate all request flags and options before the `data` promise is forced,
# preserving lazy evaluation of inactive nuisance formulas and pooling arguments.
static_public_scope_gate <- function(yname,
                                     tname,
                                     idname,
                                     dname,
                                     panel,
                                     allow_unbalanced_panel,
                                     control_group,
                                     anticipation,
                                     weightsname,
                                     fix_weights,
                                     alp,
                                     bstrap,
                                     cband,
                                     biters,
                                     est_method,
                                     base_period,
                                     faster_mode,
                                     print_details,
                                     pl,
                                     cores,
                                     compute_inffunc,
                                     dots,
                                     mode,
                                     frame,
                                     selection,
                                     selection_model,
                                     selection_pooling,
                                     selection_universe,
                                     selection_rho,
                                     selection_justification,
                                     selection_est_method,
                                     imputation_pooling,
                                     imputation_scale,
                                     gmm_weighting,
                                     gmm_covariance,
                                     pt_assumption,
                                     pt_start,
                                     control,
                                     switchers,
                                     env) {
  mode <- check_mode(mode)
  static_public_need(
    is.null(mode) || identical(mode, "chained"),
    "`dname` static estimation runs in chained panel mode (`mode = NULL` or `mode = \"chained\"`).",
    "cdid_mode_error"
  )
  static_public_need(
    chained_request_lgl1(panel),
    "`panel` must be TRUE or FALSE.",
    "cdid_mode_error"
  )
  static_public_need(
    isTRUE(panel),
    "The static `dname` route requires panel data (`panel = TRUE`); repeated cross sections (`panel = FALSE`) are not supported.",
    c("cdid_static_panel_error", "cdid_chained_panel")
  )
  static_public_need(
    is.null(allow_unbalanced_panel) || chained_request_lgl1(allow_unbalanced_panel),
    "`allow_unbalanced_panel` must be NULL, TRUE or FALSE.",
    "cdid_chained_request_input_error"
  )
  if (identical(allow_unbalanced_panel, FALSE)) {
    static_public_abort(
      paste0(
        "`allow_unbalanced_panel = FALSE` is not supported on the static `dname` route: ",
        "static panel estimation retains all units with known treatment history D = (D_1, ..., D_T), ",
        "including units with missing outcomes. Leave `allow_unbalanced_panel = NULL` or `TRUE`."
      ),
      c("cdid_static_balanced_error", "cdid_chained_balanced")
    )
  }

  static_public_need(
    length(dots) == 0L,
    "Unused `...` arguments are not supported on the static `dname` route.",
    c("cdid_static_dots_error", "cdid_chained_request_dots_error", "cdid_chained_request_input_error")
  )

  # Validate est_method before forcing data
  if (is.function(est_method)) {
    static_public_abort(
      paste0(
        "Custom function `est_method` callbacks are not supported on the static `dname` route. ",
        "Use built-in `est_method = \"dr\"`, `est_method = \"ipw\"`, or `est_method = \"reg\"`."
      ),
      c("cdid_static_est_method_error", "cdid_chained_request_est_method", "cdid_chained_unavailable", "cdid_static_unavailable")
    )
  }
  if (!chained_request_chr1(est_method) || !(est_method %in% c("dr", "ipw", "reg"))) {
    bad_m <- if (chained_request_chr1(est_method)) paste0("\"", est_method, "\"") else paste(class(est_method), collapse = "/")
    static_public_abort(
      paste0(
        "`est_method = ", bad_m, "` is not implemented on the static `dname` route. ",
        "Supported built-in methods are `est_method = \"dr\"` (default), `est_method = \"ipw\"`, and `est_method = \"reg\"`."
      ),
      c("cdid_static_est_method_error", "cdid_chained_request_est_method", "cdid_chained_unavailable", "cdid_static_unavailable")
    )
  }
  est_method_canonical <- est_method

  # Refuse inapplicable cohort / frame / weight / anticipation / PT options with clear diagnostics
  has_frame <- !arg_missing("frame", env) && !identical(substitute(frame, env), NULL) && !is.null(frame)
  if (has_frame) {
    static_public_abort(
      paste0(
        "`frame` is not supported on the static `dname` route: every study unit must have a ",
        "known binary treatment path D = (D_1, ..., D_T) across all panel dates in `data` ",
        "(units with zero or one observed outcome Y are already retained directly from `data`)."
      ),
      c("cdid_static_frame_error", "cdid_chained_request_frame_error", "cdid_static_unavailable")
    )
  }

  has_weights <- (!arg_missing("weightsname", env) && !identical(substitute(weightsname, env), NULL) && !is.null(weightsname)) ||
    (!arg_missing("fix_weights", env) && !identical(substitute(fix_weights, env), NULL) && !is.null(fix_weights))
  if (has_weights) {
    static_public_abort(
      "Sampling weights (`weightsname`, `fix_weights`) are not supported on the static `dname` route; static estimation uses equal unit weights.",
      c("cdid_static_weights_error", "cdid_chained_request_unsupported_weights", "cdid_static_unavailable")
    )
  }

  if (!arg_missing("anticipation", env)) {
    static_public_need(
      chained_request_whole(anticipation, 0),
      "`anticipation` must be a nonnegative whole number.",
      "cdid_chained_request_input_error"
    )
    if (as.integer(anticipation) != 0L) {
      static_public_abort(
        "Non-zero `anticipation` is inapplicable to the static no-carryover `dname` route (each transition at date t is anchored at adjacent reference date t - 1).",
        c("cdid_static_anticipation_error", "cdid_chained_request_input_error", "cdid_static_unavailable")
      )
    }
  }

  if (!arg_missing("pt_assumption", env)) {
    static_public_need(
      chained_request_chr1(pt_assumption) && (pt_assumption %in% c("post", "all")),
      "`pt_assumption` must be \"post\" or \"all\".",
      "cdid_chained_request_pt_start"
    )
    if (identical(pt_assumption, "all")) {
      static_public_abort(
        "`pt_assumption = \"all\"` is inapplicable to the static `dname` route, which identifies switch-time transitions relative to adjacent reference date t - 1 rather than an absorbing pre-adoption window.",
        c("cdid_static_pt_error", "cdid_chained_request_pt_start", "cdid_static_unavailable")
      )
    }
  }
  if (!arg_missing("pt_start", env) && !identical(substitute(pt_start, env), NULL) && !is.null(pt_start)) {
    static_public_abort(
      "`pt_start` is inapplicable to the static `dname` route (`pt_assumption = \"all\"` is a cohort `gname` option).",
      c("cdid_static_pt_error", "cdid_chained_request_pt_start", "cdid_static_unavailable")
    )
  }

  if (!arg_missing("control_group", env)) {
    static_public_abort(
      paste0(
        "`control_group` ('nevertreated' / 'notyettreated') is a cohort (`gname`) option and is inapplicable ",
        "to the static `dname` route, where each transition cell (t, d) uses same-prior-state stayers ",
        "C_c = 1{D_{t-1} = d, D_t = d} as its control group."
      ),
      c("cdid_static_cohort_option_error", "cdid_chained_request_input_error", "cdid_static_unavailable")
    )
  }

  if (!arg_missing("base_period", env)) {
    static_public_need(
      chained_request_chr1(base_period) && (base_period %in% c("varying", "universal")),
      "`base_period` must be \"varying\" or \"universal\".",
      "cdid_chained_request_input_error"
    )
    if (!identical(base_period, "varying")) {
      static_public_abort(
        "`base_period = \"universal\"` is a cohort (`gname`) option and is inapplicable to the static `dname` route (each transition at date t is anchored at t - 1).",
        c("cdid_static_cohort_option_error", "cdid_chained_request_input_error", "cdid_static_unavailable")
      )
    }
  }

  if (!arg_missing("selection_universe", env)) {
    static_public_need(
      chained_request_chr1(selection_universe) && (selection_universe %in% c("F", "B", "A")),
      "`selection_universe` must be \"F\", \"B\", or \"A\".",
      "cdid_chained_request_unsupported_selection"
    )
    if (!identical(selection_universe, "F")) {
      static_public_abort(
        "Sub-universe selection (`selection_universe != \"F\"`) is not supported on the static `dname` route; all units with known treatment histories belong to the study population.",
        c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection", "cdid_static_unavailable")
      )
    }
  }
  if (!arg_missing("selection_rho", env) && !identical(substitute(selection_rho, env), NULL) && !is.null(selection_rho)) {
    static_public_abort(
      "`selection_rho` is not supported on the static `dname` route.",
      c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection", "cdid_static_unavailable")
    )
  }
  if (!arg_missing("selection_justification", env) && !identical(substitute(selection_justification, env), NULL) && !is.null(selection_justification)) {
    static_public_abort(
      "`selection_justification` is not supported on the static `dname` route.",
      c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection", "cdid_static_unavailable")
    )
  }

  # Validate selection and selection_est_method, keeping inactive options lazy
  static_public_need(
    chained_request_chr1(selection) && (selection %in% c("none", "conditional_mtar")),
    "`selection` must be \"none\" or \"conditional_mtar\".",
    c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection")
  )
  static_public_need(
    chained_request_chr1(selection_est_method) &&
      (selection_est_method %in% c("ipw", "imputation", "dr")),
    "`selection_est_method` must be \"ipw\", \"imputation\", or \"dr\".",
    c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection")
  )
  sel_est_canonical <- selection_est_method

  sel_model_use <- NULL
  sel_pool_use <- NULL
  imp_pool_use <- NULL
  imp_scale_use <- NULL

  if (identical(selection, "none")) {
    static_public_need(
      identical(sel_est_canonical, "ipw"),
      "`selection = \"none\"` requires `selection_est_method = \"ipw\"`; set `selection = \"conditional_mtar\"` to use weighting, imputation, or augmented observation correction.",
      c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection")
    )
    route <- "baseline"
  } else {
    if (identical(sel_est_canonical, "ipw")) {
      route <- "weighting"
      static_public_need(
        chained_request_chr1(selection_model) && (selection_model %in% c("pair_logit", "independent_period_logit")),
        "`selection_model` must be \"pair_logit\" or \"independent_period_logit\".",
        c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection")
      )
      sel_model_use <- selection_model
      sel_pool_use <- static_public_normalize_pooling(selection_pooling, "selection_pooling")
    } else if (identical(sel_est_canonical, "imputation")) {
      route <- "imputation"
      imp_pool_use <- static_public_normalize_pooling(imputation_pooling, "imputation_pooling")
      static_public_need(
        chained_request_chr1(imputation_scale) && (imputation_scale %in% c("per_period", "per_change")),
        "`imputation_scale` must be \"per_period\" or \"per_change\".",
        c("cdid_static_imputation_scale_error", "cdid_chained_request_input_error")
      )
      imp_scale_use <- imputation_scale
    } else {
      # sel_est_canonical == "dr" (augmented observation correction; supports est_method in {"dr", "ipw", "reg"})
      route <- "augmented"
      static_public_need(
        chained_request_chr1(selection_model) && (selection_model %in% c("pair_logit", "independent_period_logit")),
        "`selection_model` must be \"pair_logit\" or \"independent_period_logit\".",
        c("cdid_static_selection_error", "cdid_chained_request_unsupported_selection")
      )
      sel_model_use <- selection_model
      sel_pool_use <- static_public_normalize_pooling(selection_pooling, "selection_pooling")
      imp_pool_use <- static_public_normalize_pooling(imputation_pooling, "imputation_pooling")
      static_public_need(
        chained_request_chr1(imputation_scale) && (imputation_scale %in% c("per_period", "per_change")),
        "`imputation_scale` must be \"per_period\" or \"per_change\".",
        c("cdid_static_imputation_scale_error", "cdid_chained_request_input_error")
      )
      imp_scale_use <- imputation_scale
    }
  }

  # Validate inference and execution flags
  static_public_need(
    chained_request_lgl1(bstrap) && chained_request_lgl1(cband) &&
      chained_request_lgl1(faster_mode) && chained_request_lgl1(print_details) &&
      chained_request_lgl1(pl) && chained_request_lgl1(compute_inffunc),
    "`bstrap`, `cband`, `faster_mode`, `print_details`, `pl`, and `compute_inffunc` must be scalar logicals (TRUE or FALSE).",
    "cdid_chained_request_input_error"
  )
  static_public_need(
    chained_request_whole(biters, 1) && biters <= .Machine$integer.max,
    "`biters` must be a positive whole number no greater than .Machine$integer.max.",
    "cdid_chained_request_input_error"
  )
  static_public_need(
    chained_request_whole(cores, 1) && cores <= .Machine$integer.max,
    "`cores` must be a positive whole number no greater than .Machine$integer.max.",
    "cdid_chained_request_input_error"
  )
  static_public_need(
    is.numeric(alp) && length(alp) == 1L && is.finite(alp) && alp > 0 && alp < 1,
    "`alp` must be a finite scalar in (0, 1).",
    "cdid_chained_request_input_error"
  )

  ctrl_resolved <- chained_public_resolve_control(
    control = control,
    gmm_weighting = gmm_weighting,
    gmm_covariance = gmm_covariance
  )
  parallel_plan <- chained_request_parallel_plan(
    pl = pl,
    cores = cores,
    workers_cap = ctrl_resolved$profile$workers
  )

  # Validate column name arguments
  static_public_need(
    chained_request_chr1(yname) && chained_request_chr1(tname) &&
      chained_request_chr1(idname) && chained_request_chr1(dname),
    "`yname`, `tname`, `idname`, and `dname` must be non-empty character scalars.",
    "cdid_chained_request_input_error"
  )
  static_public_need(
    length(unique(c(yname, tname, idname, dname))) == 4L,
    "`yname`, `tname`, `idname`, and `dname` must name distinct columns.",
    "cdid_chained_request_input_error"
  )

  # Validate switchers direction filter
  static_public_need(
    chained_request_chr1(switchers) && (switchers %in% c("both", "entry", "exit")),
    "`switchers` must be one of \"both\", \"entry\", or \"exit\".",
    c("cdid_static_switchers_error", "cdid_chained_request_input_error")
  )

  list(
    mode = "chained",
    est_method_requested = est_method,
    est_method = est_method_canonical,
    selection_requested = selection,
    selection_canonical = selection,
    selection_est_method_requested = selection_est_method,
    selection_est_method = sel_est_canonical,
    route = route,
    selection_model = sel_model_use,
    selection_pooling = sel_pool_use,
    imputation_pooling = imp_pool_use,
    imputation_scale = imp_scale_use,
    ctrl_resolved = ctrl_resolved,
    parallel_plan = parallel_plan,
    switchers = switchers
  )
}

# Extract and validate optional unit-constant cluster vector aligned to `panel$unit_id`.
static_public_extract_clusters <- function(data, panel, idname, clustervars) {
  if (is.null(clustervars)) {
    return(list(
      clustered = FALSE,
      clustername = NULL,
      cluster_vec = NULL,
      n_clusters = panel$n
    ))
  }
  static_public_need(
    is.character(clustervars) && is.null(dim(clustervars)) && !is.object(clustervars) &&
      length(clustervars) >= 1L && !anyNA(clustervars) && all(nzchar(clustervars)),
    "`clustervars` must be NULL or a character vector of column names.",
    c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
  )
  extra <- unique(clustervars[clustervars != idname])
  static_public_need(
    length(extra) <= 1L,
    "At most one cluster variable beyond `idname` is supported in `clustervars`.",
    c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
  )
  if (length(extra) == 0L) {
    return(list(
      clustered = FALSE,
      clustername = NULL,
      cluster_vec = NULL,
      n_clusters = panel$n
    ))
  }
  cname <- extra[1L]
  static_public_need(
    cname %in% names(data),
    sprintf("Cluster variable `%s` was not found in `data`.", cname),
    c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
  )
  raw_col <- data[[cname]]
  static_public_need(
    is.atomic(raw_col) && is.null(dim(raw_col)) &&
      (is.numeric(raw_col) || is.character(raw_col) || is.factor(raw_col) || is.logical(raw_col)) &&
      !anyNA(raw_col),
    sprintf("Cluster variable `%s` must be an atomic vector without missing values.", cname),
    c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
  )
  if (is.numeric(raw_col)) {
    static_public_need(
      all(is.finite(raw_col)),
      sprintf("Numeric cluster variable `%s` must be finite.", cname),
      c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
    )
  }
  c_chr <- as.character(raw_col)
  static_public_need(
    all(nzchar(c_chr)),
    sprintf("Cluster variable `%s` must not contain empty strings.", cname),
    c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
  )

  # Check within-unit constancy across panel rows using panel$row_index (n x T)
  n <- panel$n
  rmat <- panel$row_index
  c_unit <- c_chr[rmat[, 1L]]
  if (ncol(rmat) >= 2L) {
    for (t_j in seq.int(2L, ncol(rmat))) {
      bad_u <- which(c_chr[rmat[, t_j]] != c_unit)
      if (length(bad_u) > 0L) {
        u_bad <- panel$unit_id[bad_u[1L]]
        static_public_abort(
          sprintf(
            "Cluster variable `%s` varies across time within unit `%s`; clusters must be unit-constant.",
            cname, as.character(u_bad)
          ),
          c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
        )
      }
    }
  }
  names(c_unit) <- panel$unit_key
  n_cl <- length(unique(c_unit))
  static_public_need(
    n_cl >= 2L,
    sprintf("Cluster variable `%s` must define at least 2 distinct clusters (got %d).", cname, n_cl),
    c("cdid_static_cluster_error", "cdid_chained_request_cluster_error")
  )
  list(
    clustered = TRUE,
    clustername = cname,
    cluster_vec = c_unit,
    n_clusters = as.integer(n_cl)
  )
}

# Resolve declared target cells A from `switchers` and `switch_times`.
# Auto target discovery (`switch_times = NULL`) retains all candidate transitions
# in the requested direction(s) with `n_switchers >= 1L` across all dates t >= 2.
# Explicit `switch_times` retains requested cells even when `n_switchers == 0L`
# so they are reported as `available = FALSE` (`empty_switchers`), never silently dropped.
static_public_resolve_targets <- function(panel,
                                          all_cand,
                                          switchers,
                                          switch_times) {
  cand_cells <- all_cand$cells
  cand_summary <- all_cand$support_table
  req_d <- if (identical(switchers, "entry")) {
    0L
  } else if (identical(switchers, "exit")) {
    1L
  } else {
    c(0L, 1L)
  }

  post_dates_num <- panel$dates_num[-1L]
  post_dates_lbl <- panel$date_labels[-1L]

  if (is.null(switch_times)) {
    in_dir <- cand_summary$d %in% req_d
    non_empty <- in_dir & (cand_summary$n_switchers >= 1L)
    if (any(non_empty)) {
      sel_ids <- cand_summary$cell_id[non_empty]
    } else {
      # No switchers observed anywhere in the requested direction(s): retain
      # candidate cells in requested direction(s) as unavailable so the result
      # reports them honestly rather than returning an empty vector.
      sel_ids <- cand_summary$cell_id[in_dir]
    }
    return(list(
      target_ids = sel_ids,
      target_cells = cand_cells[sel_ids],
      discovery_mode = "automatic_nonempty_switches"
    ))
  }

  # Explicit `switch_times` supplied by caller
  if (is.character(switch_times) && is.null(dim(switch_times)) && !is.object(switch_times)) {
    static_public_need(
      length(switch_times) >= 1L && !anyNA(switch_times) && all(nzchar(switch_times)),
      "`switch_times` must be NULL, a numeric vector of switch dates, or a character vector of date labels.",
      c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
    )
    t_match <- match(switch_times, post_dates_lbl)
    if (anyNA(t_match)) {
      bad_t <- switch_times[which(is.na(t_match))[1L]]
      if (bad_t %in% panel$date_labels[1L]) {
        static_public_abort(
          sprintf(
            "`switch_times` cannot include the first study date `%s` (transitions require a prior date t - 1).",
            bad_t
          ),
          c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
        )
      }
      static_public_abort(
        sprintf(
          "`switch_times` entry `%s` is not in the panel transition dates [%s].",
          bad_t, paste(post_dates_lbl, collapse = ", ")
        ),
        c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
      )
    }
    req_t_num <- post_dates_num[t_match]
  } else if (is.numeric(switch_times) && is.null(dim(switch_times)) && !is.object(switch_times)) {
    static_public_need(
      length(switch_times) >= 1L && !anyNA(switch_times) && all(is.finite(switch_times)),
      "`switch_times` numeric vector must be non-empty and finite.",
      c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
    )
    t_num_in <- as.double(switch_times)
    if (any(t_num_in == panel$dates_num[1L])) {
      static_public_abort(
        sprintf(
          "`switch_times` cannot include the first study date `%s` (transitions require a prior date t - 1).",
          panel$date_labels[1L]
        ),
        c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
      )
    }
    t_match <- match(t_num_in, post_dates_num)
    if (anyNA(t_match)) {
      bad_t <- switch_times[which(is.na(t_match))[1L]]
      static_public_abort(
        sprintf(
          "`switch_times` date `%s` is not in the panel transition dates [%s].",
          as.character(bad_t), paste(post_dates_lbl, collapse = ", ")
        ),
        c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
      )
    }
    req_t_num <- post_dates_num[t_match]
  } else {
    static_public_abort(
      "`switch_times` must be NULL or a vector of transition dates in `panel$dates[-1]`.",
      c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
    )
  }

  static_public_need(
    !anyDuplicated(req_t_num),
    "`switch_times` must not contain duplicate dates.",
    c("cdid_static_switch_times_error", "cdid_chained_request_input_error")
  )
  req_t_sorted <- sort(req_t_num)

  sel_ids <- character(0L)
  for (tv in req_t_sorted) {
    rows_t <- cand_summary[cand_summary$t == tv & (cand_summary$d %in% req_d), , drop = FALSE]
    # Every direction requested at an explicit date remains declared, including
    # empty cells. Omitted and explicit `switchers = "both"` are equivalent.
    sel_ids <- c(sel_ids, rows_t$cell_id)
  }

  list(
    target_ids = sel_ids,
    target_cells = cand_cells[sel_ids],
    discovery_mode = "explicit_switch_times"
  )
}

# Resolve `outcome_pooling` into an explicit OR partition and identify any
# extra nuisance-only stayer training cells in U \ A.
static_public_resolve_outcome_pooling <- function(outcome_pooling,
                                                  target_ids,
                                                  all_cand_cells) {
  cand_ids <- names(all_cand_cells)
  if (is.character(outcome_pooling) && length(outcome_pooling) == 1L && is.null(names(outcome_pooling))) {
    static_public_need(
      outcome_pooling %in% c("separate", "pooled"),
      "`outcome_pooling` keyword must be \"separate\" or \"pooled\" (or an explicit named partition over cell IDs).",
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    or_map <- if (identical(outcome_pooling, "separate")) {
      stats::setNames(target_ids, target_ids)
    } else {
      stats::setNames(rep("pooled", length(target_ids)), target_ids)
    }
    return(list(
      mode = outcome_pooling,
      or_blocks = or_map,
      all_ids = target_ids,
      non_target_ids = character(0L)
    ))
  }

  if (is.character(outcome_pooling) && !is.null(names(outcome_pooling))) {
    nm <- names(outcome_pooling)
    static_public_need(
      !anyNA(nm) && all(nzchar(nm)) && !anyDuplicated(nm) &&
        !anyNA(outcome_pooling) && all(nzchar(outcome_pooling)),
      "Named character `outcome_pooling` partition must have unique non-empty cell ID names and non-empty block labels.",
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    unknown <- setdiff(nm, cand_ids)
    static_public_need(
      length(unknown) == 0L,
      sprintf(
        "`outcome_pooling` references unknown transition cell ID(s): %s. Valid cell IDs are: %s.",
        paste(unknown, collapse = ", "), paste(cand_ids, collapse = ", ")
      ),
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    miss_t <- setdiff(target_ids, nm)
    static_public_need(
      length(miss_t) == 0L,
      sprintf("`outcome_pooling` omits declared target cell(s): %s.", paste(miss_t, collapse = ", ")),
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    extra_ids <- setdiff(nm, target_ids)
    all_ids <- c(target_ids, extra_ids)
    return(list(
      mode = "custom_named_vector",
      or_blocks = outcome_pooling[all_ids],
      all_ids = all_ids,
      non_target_ids = extra_ids
    ))
  }

  if (is.list(outcome_pooling) && !is.object(outcome_pooling) && length(outcome_pooling) >= 1L) {
    b_names <- names(outcome_pooling)
    if (is.null(b_names)) {
      b_names <- paste0("or_block_", seq_along(outcome_pooling))
    }
    static_public_need(
      !anyNA(b_names) && all(nzchar(b_names)) && !anyDuplicated(b_names),
      "List `outcome_pooling` partition must have unique non-empty block names.",
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    or_map <- character(0L)
    for (k in seq_along(outcome_pooling)) {
      members <- outcome_pooling[[k]]
      static_public_need(
        is.character(members) && is.null(dim(members)) && length(members) >= 1L &&
          !anyNA(members) && all(nzchar(members)),
        sprintf("Block `%s` in `outcome_pooling` must be a non-empty character vector of cell IDs.", b_names[k]),
        c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
      )
      or_map <- c(or_map, stats::setNames(rep(b_names[k], length(members)), members))
    }
    nm <- names(or_map)
    static_public_need(
      !anyDuplicated(nm),
      "`outcome_pooling` blocks must be disjoint (a cell ID cannot belong to multiple OR blocks).",
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    unknown <- setdiff(nm, cand_ids)
    static_public_need(
      length(unknown) == 0L,
      sprintf(
        "`outcome_pooling` references unknown transition cell ID(s): %s. Valid cell IDs are: %s.",
        paste(unknown, collapse = ", "), paste(cand_ids, collapse = ", ")
      ),
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    miss_t <- setdiff(target_ids, nm)
    static_public_need(
      length(miss_t) == 0L,
      sprintf("`outcome_pooling` omits declared target cell(s): %s.", paste(miss_t, collapse = ", ")),
      c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
    )
    extra_ids <- setdiff(nm, target_ids)
    all_ids <- c(target_ids, extra_ids)
    return(list(
      mode = "custom_block_list",
      or_blocks = or_map[all_ids],
      all_ids = all_ids,
      non_target_ids = extra_ids
    ))
  }

  static_public_abort(
    "`outcome_pooling` must be \"separate\", \"pooled\", a named character vector mapping cell_id -> block_id, or a named list of cell_id vectors.",
    c("cdid_static_outcome_pooling_error", "cdid_static_missing_partition_error")
  )
}

# Build the all-positive-observed-count pair graph for one (cell, role) population.
# Retains EVERY date pair (s, u) with s < u and `n_observed_role >= 1L` across all
# connected components (`K_c^a >= 1`); never prunes disconnected extra components
# or multiple paths.
static_public_build_one_role_graph <- function(panel, cell, role, cand_s, cand_u) {
  cid <- cell$cell_id
  role_mask <- if (identical(role, "G")) cell$G else cell$C
  n_role <- as.integer(sum(role_mask))
  base_t <- cell$t_prev_num
  targ_t <- cell$t_num

  acc <- static_observation_pair_accessor(
    panel = panel,
    mask = role_mask,
    s = cand_s,
    u = cand_u
  )
  st <- acc$summary_table
  ret_df <- st[st$n_observed_role >= 1L, , drop = FALSE]
  K_cand <- nrow(st)
  K_ret <- nrow(ret_df)

  if (K_ret == 0L) {
    return(list(
      cell_id = cid,
      role = role,
      baseline_time = base_t,
      target_time = targ_t,
      n_role_units = n_role,
      n_candidate_pairs = as.integer(K_cand),
      n_retained_edges = 0L,
      active_dates = numeric(0L),
      n_active_dates = 0L,
      n_components = 0L,
      cycle_rank = 0L,
      endpoints_incident = FALSE,
      endpoints_connected = FALSE,
      retained_pairs = ret_df,
      candidate_pairs = st,
      block_spec = NULL
    ))
  }

  active_dates <- sort(unique(c(ret_df$s, ret_df$u)))
  from_idx <- match(ret_df$s, active_dates)
  to_idx <- match(ret_df$u, active_dates)
  base_v_idx <- match(base_t, active_dates)
  targ_v_idx <- match(targ_t, active_dates)
  anchor_use <- if (!is.na(base_v_idx)) base_v_idx else 1L

  forest <- incidence_forest(
    nv = length(active_dates),
    from = from_idx,
    to = to_idx,
    anchors = anchor_use,
    edge_order = seq_len(K_ret)
  )
  n_comp <- as.integer(forest$n_components)
  n_vert <- length(active_dates)
  cycle_rank <- as.integer(K_ret - (n_vert - n_comp))
  endpoints_incident <- !is.na(base_v_idx) && !is.na(targ_v_idx) &&
    (forest$degree[base_v_idx] >= 1L) && (forest$degree[targ_v_idx] >= 1L)
  endpoints_connected <- endpoints_incident &&
    (forest$component_id[base_v_idx] == forest$component_id[targ_v_idx])

  b_spec <- NULL
  if (endpoints_connected) {
    edges_clean <- data.frame(
      link_id = sprintf(
        "%s_%s_%s_%s",
        cid, role,
        panel$date_labels[ret_df$s_idx],
        panel$date_labels[ret_df$u_idx]
      ),
      s = ret_df$s,
      u = ret_df$u,
      stringsAsFactors = FALSE
    )
    b_spec <- static_role_block_spec(
      block_id = paste0(cid, "__", role),
      cell_id = cid,
      role = role,
      baseline_time = base_t,
      target_time = targ_t,
      dates = active_dates,
      edges = edges_clean
    )
  }

  list(
    cell_id = cid,
    role = role,
    baseline_time = base_t,
    target_time = targ_t,
    n_role_units = n_role,
    n_candidate_pairs = as.integer(K_cand),
    n_retained_edges = as.integer(K_ret),
    active_dates = active_dates,
    n_active_dates = as.integer(n_vert),
    n_components = n_comp,
    cycle_rank = cycle_rank,
    endpoints_incident = endpoints_incident,
    endpoints_connected = endpoints_connected,
    retained_pairs = ret_df,
    candidate_pairs = st,
    block_spec = b_spec
  )
}

# Build all-pairs role graphs and evaluate cell availability across target and
# non-target OR training cells.
static_public_build_graphs_and_support <- function(panel,
                                                   target_ids,
                                                   non_target_ids,
                                                   all_cand_cells) {
  T_cnt <- length(panel$dates_num)
  idx_mat <- utils::combn(seq_len(T_cnt), 2L)
  cand_s <- panel$dates_num[idx_mat[1L, ]]
  cand_u <- panel$dates_num[idx_mat[2L, ]]

  all_req_ids <- c(target_ids, non_target_ids)
  by_cell_role <- list()
  graph_rows <- list()

  J <- length(target_ids)
  available <- logical(J)
  support_status <- character(J)
  support_reason <- character(J)
  names(available) <- target_ids
  names(support_status) <- target_ids
  names(support_reason) <- target_ids

  for (j in seq_along(target_ids)) {
    cid <- target_ids[j]
    cl <- all_cand_cells[[cid]]
    g_info <- static_public_build_one_role_graph(panel, cl, "G", cand_s, cand_u)
    c_info <- static_public_build_one_role_graph(panel, cl, "C", cand_s, cand_u)
    by_cell_role[[paste0(cid, "::G")]] <- g_info
    by_cell_role[[paste0(cid, "::C")]] <- c_info

    graph_rows[[length(graph_rows) + 1L]] <- data.frame(
      cell_id = cid,
      role = "G",
      is_target_cell = TRUE,
      baseline_time = cl$t_prev_num,
      target_time = cl$t_num,
      n_role_units = g_info$n_role_units,
      n_candidate_pairs = g_info$n_candidate_pairs,
      n_retained_edges = g_info$n_retained_edges,
      n_active_dates = g_info$n_active_dates,
      n_components = g_info$n_components,
      cycle_rank = g_info$cycle_rank,
      endpoints_connected = g_info$endpoints_connected,
      stringsAsFactors = FALSE
    )
    graph_rows[[length(graph_rows) + 1L]] <- data.frame(
      cell_id = cid,
      role = "C",
      is_target_cell = TRUE,
      baseline_time = cl$t_prev_num,
      target_time = cl$t_num,
      n_role_units = c_info$n_role_units,
      n_candidate_pairs = c_info$n_candidate_pairs,
      n_retained_edges = c_info$n_retained_edges,
      n_active_dates = c_info$n_active_dates,
      n_components = c_info$n_components,
      cycle_rank = c_info$cycle_rank,
      endpoints_connected = c_info$endpoints_connected,
      stringsAsFactors = FALSE
    )

    if (cl$n_switchers < 1L) {
      available[j] <- FALSE
      support_status[j] <- "empty_switchers"
      support_reason[j] <- sprintf(
        "No switchers observed at t = %s from prior treatment D_{t-1} = %d.",
        cl$t_label, cl$d
      )
    } else if (cl$n_stayers < 1L) {
      available[j] <- FALSE
      support_status[j] <- "empty_stayers"
      support_reason[j] <- sprintf(
        "No stayers observed at t = %s from prior treatment D_{t-1} = %d.",
        cl$t_label, cl$d
      )
    } else if (!isTRUE(g_info$endpoints_connected)) {
      available[j] <- FALSE
      support_status[j] <- "disconnected_switcher_graph"
      support_reason[j] <- sprintf(
        "Switcher (G) observation-pair graph for cell %s does not connect baseline %s and target %s.",
        cid, cl$t_prev_label, cl$t_label
      )
    } else if (!isTRUE(c_info$endpoints_connected)) {
      available[j] <- FALSE
      support_status[j] <- "disconnected_stayer_graph"
      support_reason[j] <- sprintf(
        "Stayer (C) observation-pair graph for cell %s does not connect baseline %s and target %s.",
        cid, cl$t_prev_label, cl$t_label
      )
    } else {
      available[j] <- TRUE
      support_status[j] <- "ok"
      support_reason[j] <- NA_character_
    }
  }

  for (cid in non_target_ids) {
    cl <- all_cand_cells[[cid]]
    if (cl$n_stayers < 1L) {
      static_public_abort(
        sprintf(
          "Declared OR training cell `%s` (t = %s, d = %d) in `outcome_pooling` has zero stayers (`n_stayers = 0`).",
          cid, cl$t_label, cl$d
        ),
        c("cdid_static_support_error", "cdid_static_missing_support_error")
      )
    }
    c_info <- static_public_build_one_role_graph(panel, cl, "C", cand_s, cand_u)
    if (!isTRUE(c_info$endpoints_connected)) {
      static_public_abort(
        sprintf(
          "Declared OR training cell `%s` (t = %s, d = %d) in `outcome_pooling` has a disconnected stayer (`C`) observation-pair graph between %s and %s.",
          cid, cl$t_label, cl$d, cl$t_prev_label, cl$t_label
        ),
        c("cdid_static_graph_disconnected_error", "cdid_static_support_error")
      )
    }
    by_cell_role[[paste0(cid, "::C")]] <- c_info
    graph_rows[[length(graph_rows) + 1L]] <- data.frame(
      cell_id = cid,
      role = "C",
      is_target_cell = FALSE,
      baseline_time = cl$t_prev_num,
      target_time = cl$t_num,
      n_role_units = c_info$n_role_units,
      n_candidate_pairs = c_info$n_candidate_pairs,
      n_retained_edges = c_info$n_retained_edges,
      n_active_dates = c_info$n_active_dates,
      n_components = c_info$n_components,
      cycle_rank = c_info$cycle_rank,
      endpoints_connected = c_info$endpoints_connected,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, graph_rows)
  rownames(summary_df) <- NULL

  # Assemble `role_graphs` for `static_missing_fit` on available target cells + non-target OR cells
  avail_target_ids <- target_ids[available]
  role_graphs_fit <- list()
  for (cid in avail_target_ids) {
    role_graphs_fit[[cid]] <- list(
      G = by_cell_role[[paste0(cid, "::G")]]$block_spec,
      C = by_cell_role[[paste0(cid, "::C")]]$block_spec
    )
  }
  for (cid in non_target_ids) {
    role_graphs_fit[[cid]] <- list(
      C = by_cell_role[[paste0(cid, "::C")]]$block_spec
    )
  }

  list(
    available = available,
    support_status = support_status,
    support_reason = support_reason,
    avail_target_ids = avail_target_ids,
    role_graphs_fit = role_graphs_fit,
    summary = summary_df,
    by_cell_role = by_cell_role
  )
}

# Analytical and common-sign Rademacher multiplier bootstrap inference on
# original-unit influence functions via `chained_cluster_multiplier()`.
static_public_inference <- function(estimates,
                                    if_matrix,
                                    unit_ids,
                                    cluster_vec = NULL,
                                    alp = 0.05,
                                    bstrap = TRUE,
                                    cband = TRUE,
                                    biters = 1000L,
                                    chunk_rows = 256L,
                                    workspace_entries = 134217728,
                                    workspace_bytes = 1073741824,
                                    workers = 1L,
                                    fixed_signs = NULL,
                                    agg_estimate = NULL,
                                    agg_if = NULL,
                                    all_declared_available = TRUE) {
  Psi <- as.matrix(if_matrix)
  storage.mode(Psi) <- "double"
  n <- nrow(Psi)
  J <- ncol(Psi)
  tids <- colnames(Psi)
  if (is.null(tids)) {
    tids <- paste0("t", seq_len(J))
    colnames(Psi) <- tids
  }
  uid_chr <- as.character(unit_ids)
  rownames(Psi) <- uid_chr

  u_cl <- if (is.null(cluster_vec)) seq_len(n) else unique(as.character(cluster_vec))
  map <- if (is.null(cluster_vec)) seq_len(n) else match(as.character(cluster_vec), u_cl)
  G <- length(u_cl)

  flags <- ipw_pipeline_resolve_bootstrap_flags(bstrap, cband, endpoint = "cells")
  B <- as.integer(biters)
  chunk_u <- as.integer(max(1L, min(n, chunk_rows)))
  sign_r <- as.integer(max(1L, min(B, chunk_u, floor(workspace_entries / max(1, G)))))
  if_r <- as.integer(max(1L, min(n, chunk_u, floor(workspace_entries / max(1, J)))))
  # Explicit numeric payload accounting, including callback-owned inputs. This
  # is an entry/byte bound, not an allocator or peak-RSS guarantee.
  producer_entries <- 4.0 * n * J + 4.0 * G * J + 4.0 * J * J +
    6.0 * n + 4.0 * G +
    (if (is.null(agg_if)) 0.0 else 4.0 * n + 4.0 * G)
  consumer_entries <- 0.0
  if (flags$draws) {
    consumer_entries <- compute_multiplier_workspace_estimate(
      n_units = n, n_clusters = G, n_targets = J, n_draws = B,
      max_if_batch_rows = if_r, max_sign_batch_rows = sign_r,
      workers = workers
    )$total_persistent_cells
    producer_entries <- producer_entries + 4.0 * sign_r * G +
      (if (is.null(agg_if)) 0.0 else 4.0 * B) +
      (if (is.null(fixed_signs)) 0.0 else 4.0 * length(fixed_signs))
  }
  required_entries <- producer_entries + consumer_entries
  static_public_need(
    is.finite(required_entries) && required_entries <= workspace_entries &&
      8.0 * required_entries <= workspace_bytes,
    "Static inference producer and multiplier storage exceeds `control$numeric_workspace_bytes`.",
    c("cdid_static_inference_bound_error", "cdid_multiplier_budget_error")
  )

  # Centered cluster sums tilde{S}_g = S_g - bar{S} and analytical covariance
  S_raw <- matrix(0.0, nrow = G, ncol = J, dimnames = list(as.character(u_cl), tids))
  if (is.null(cluster_vec)) {
    S_raw[,] <- Psi
  } else {
    for (i in seq_len(n)) {
      g_i <- map[i]
      S_raw[g_i, ] <- S_raw[g_i, ] + Psi[i, ]
    }
  }
  S_cent <- sweep(S_raw, 2L, colMeans(S_raw), FUN = "-")
  V_an <- crossprod(S_cent) / (as.double(n)^2)
  dimnames(V_an) <- list(tids, tids)
  an_var <- diag(V_an)
  an_se <- ifelse(is.finite(an_var) & an_var > 0.0, sqrt(pmax(0.0, an_var)), NA_real_)
  names(an_se) <- tids

  # Analytical SE for aggregate if supplied
  agg_an_var <- NA_real_
  agg_an_se <- NA_real_
  if (!is.null(agg_if)) {
    agg_psi <- as.double(agg_if)
    S_agg_raw <- if (is.null(cluster_vec)) {
      agg_psi
    } else {
      vapply(seq_len(G), function(g) sum(agg_psi[map == g]), numeric(1L))
    }
    S_agg_cent <- S_agg_raw - mean(S_agg_raw)
    agg_an_var <- sum(S_agg_cent * S_agg_cent) / (as.double(n)^2)
    agg_an_se <- if (is.finite(agg_an_var) && agg_an_var > 0.0) sqrt(agg_an_var) else NA_real_
  }

  zcrit <- stats::qnorm(1.0 - alp / 2.0)
  run_draws <- flags$draws

  if (!run_draws) {
    return(list(
      V_analytical = V_an,
      analytical_se = an_se,
      bootstrap_se = stats::setNames(rep(NA_real_, J), tids),
      bootstrap_iqr_scale = stats::setNames(rep(NA_real_, J), tids),
      se = an_se,
      se_type = "analytical",
      pointwise_crit_val = zcrit,
      cband_crit_val = NA_real_,
      cband_status = "not_requested",
      multiplier = NULL,
      agg_analytical_var = agg_an_var,
      agg_analytical_se = agg_an_se,
      agg_bootstrap_se = NA_real_,
      agg_se = agg_an_se,
      n_clusters = as.integer(G),
      clustered = !is.null(cluster_vec)
    ))
  }

  mb <- list(
    max_units = as.double(n),
    max_clusters = as.double(G),
    max_targets = as.double(max(1L, J)),
    max_draws = as.double(B),
    max_if_batch_rows = as.double(if_r),
    max_if_batch_cells = as.double(if_r) * as.double(max(1L, J)),
    max_sign_batch_rows = as.double(sign_r),
    max_sign_batch_cells = as.double(sign_r) * as.double(G),
    max_workspace_cells = as.double(workspace_entries),
    max_workspace_bytes = as.double(workspace_bytes)
  )
  if (as.integer(workers) > 1L) {
    mb$workers <- as.double(workers)
  }

  # Draw one bounded sign batch at a time from the caller's stream. The
  # aggregate draw is accumulated from the same batch, so no B x G matrix or
  # RNG replay is needed and the aggregate does not enlarge the band family.
  n_sign_batches <- as.integer(ceiling(B / sign_r))
  if (!is.null(fixed_signs)) {
    static_public_need(
      is.matrix(fixed_signs) && is.numeric(fixed_signs) &&
        identical(dim(fixed_signs), c(B, as.integer(G))) &&
        !anyNA(fixed_signs) && all(fixed_signs == 1 | fixed_signs == -1),
      sprintf("`fixed_signs` must be a %d x %d matrix with entries strictly in {-1, +1}.", B, G),
      "cdid_static_inference_error"
    )
  }
  aggregate_draws <- if (!is.null(agg_if) && isTRUE(all_declared_available)) {
    rep(NA_real_, B)
  } else {
    NULL
  }

  ids_int <- seq_len(n)
  cids_int <- seq_len(G)
  dids_int <- seq_len(B)
  n_if_batches <- as.integer(ceiling(n / if_r))

  mult <- chained_cluster_multiplier(
    target_ids = tids,
    expected_unit_ids = ids_int,
    unit_cluster_map = list(unit_ids = ids_int, cluster_ids = as.integer(map)),
    cluster_ids = cids_int,
    n_denom = as.double(n),
    center = TRUE,
    if_batch_fn = function(k) {
      ii <- seq.int((k - 1L) * if_r + 1L, min(k * if_r, n))
      list(
        unit_ids = ids_int[ii],
        target_ids = tids,
        if_matrix = unname(Psi[ii, , drop = FALSE])
      )
    },
    n_if_batches = n_if_batches,
    sign_batch_fn = function(k) {
      ii <- seq.int((k - 1L) * sign_r + 1L, min(k * sign_r, B))
      signs <- if (is.null(fixed_signs)) {
        matrix(sample(c(-1.0, 1.0), length(ii) * G, replace = TRUE),
               nrow = length(ii), ncol = G)
      } else {
        unname(fixed_signs[ii, , drop = FALSE])
      }
      if (!is.null(aggregate_draws)) {
        aggregate_draws[ii] <<- drop(signs %*% S_agg_cent) / sqrt(as.double(G))
      }
      list(draw_ids = dids_int[ii], cluster_ids = cids_int,
           sign_matrix = signs)
    },
    n_sign_batches = n_sign_batches,
    expected_draw_ids = dids_int,
    n_draws = B,
    alpha = as.double(alp),
    return_raw_draws = FALSE,
    budget = mb
  )

  # Per-target bootstrap SEs (candidate_bootstrap_se is populated for every
  # target whose own analytical and bootstrap scale are positive)
  boot_se <- ifelse(
    mult$target_accepted,
    as.double(mult$candidate_bootstrap_se),
    NA_real_
  )
  names(boot_se) <- tids

  boot_iqr_scale <- ifelse(
    mult$target_accepted,
    as.double(mult$raw_bootstrap_scale),
    NA_real_
  )
  names(boot_iqr_scale) <- tids

  # Simultaneous band over the FULL declared target family (no survival subfamily)
  cband_crit <- NA_real_
  cband_status <- if (!isTRUE(flags$simultaneous)) {
    "not_requested"
  } else if (!isTRUE(all_declared_available)) {
    "unavailable_target_cell"
  } else if (!isTRUE(mult$accepted) || !is.finite(mult$crit_val)) {
    paste0("unavailable_degenerate_family:", paste(mult$refusal_reasons, collapse = ","))
  } else {
    cband_crit <- as.double(mult$crit_val)
    "ok"
  }

  # The scalar aggregate uses the same type-1 IQR and parameter scaling as
  # the accepted worker, using its common sign batches retained above.
  agg_boot_se <- NA_real_
  if (!is.null(aggregate_draws) && is.finite(agg_an_se) && agg_an_se > 0.0 &&
      all(is.finite(aggregate_draws))) {
    quartiles <- stats::quantile(aggregate_draws, probs = c(0.25, 0.75),
                                type = 1L, names = FALSE)
    scale <- (quartiles[2L] - quartiles[1L]) /
      (stats::qnorm(0.75) - stats::qnorm(0.25))
    candidate <- scale * sqrt(as.double(G)) / as.double(n)
    if (is.finite(candidate) && candidate > 0.0) agg_boot_se <- candidate
  }

  rep_se <- if (isTRUE(flags$bootstrap)) boot_se else an_se
  rep_agg_se <- if (isTRUE(flags$bootstrap)) agg_boot_se else agg_an_se

  list(
    V_analytical = V_an,
    analytical_se = an_se,
    bootstrap_se = boot_se,
    bootstrap_iqr_scale = boot_iqr_scale,
    se = rep_se,
    se_type = if (isTRUE(flags$bootstrap)) "multiplier_bootstrap_iqr" else "analytical",
    pointwise_crit_val = zcrit,
    cband_crit_val = cband_crit,
    cband_status = cband_status,
    multiplier = mult,
    agg_analytical_var = agg_an_var,
    agg_analytical_se = agg_an_se,
    agg_bootstrap_se = agg_boot_se,
    agg_se = rep_agg_se,
    n_clusters = as.integer(G),
    clustered = !is.null(cluster_vec)
  )
}

static_public_fit_complete_baseline <- function(panel,
                                                avail_target_ids,
                                                non_target_ids,
                                                all_cand_cells,
                                                or_blocks_fit,
                                                xformla_use,
                                                use_or,
                                                gate,
                                                glm_ctrl,
                                                rank_tol,
                                                ws_entries,
                                                all_declared_available,
                                                avail_sigma,
                                                compute_inffunc) {
  n <- panel$n
  all_fit_ids <- c(avail_target_ids, if (use_or) non_target_ids else character(0L))

  # Known minimum live outputs, before extracting or copying any cell designs:
  # G_matrix and, for OR, both fitted/residual matrices over all training cells.
  gmm_safe_product(
    n,
    as.double(length(avail_target_ids)) + if (use_or) 2.0 * length(all_fit_ids) else 0.0,
    max_entries = ws_entries,
    label = "complete-data baseline live cell outputs"
  )

  build_complete_cell <- function(cid, allow_empty_switchers = FALSE) {
    cl <- all_cand_cells[[cid]]
    dy <- as.double(panel$Y[, cl$t_idx] - panel$Y[, cl$t_prev_idx])
    des <- static_extract_reference_design(
      panel = panel,
      cell = cl,
      formula = xformla_use,
      role = if (!allow_empty_switchers && gate$est_method %in% c("ipw", "dr")) "treatment" else "outcome",
      mask = (cl$R == 1L)
    )
    gmm_safe_product(
      n,
      ncol(des$X_tilde),
      max_entries = ws_entries,
      label = sprintf("complete-data cell `%s` design", cid)
    )
    X_safe <- des$X_tilde
    X_safe[cl$R == 0L, ] <- 0.0
    X_safe[, 1L] <- 1.0
    Z_safe <- des$Z
    if (ncol(Z_safe) > 0L) {
      Z_safe[cl$R == 0L, ] <- 0.0
    }
    static_complete_cell(
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
      allow_empty_switchers = allow_empty_switchers
    )
  }

  comp_cells <- vector("list", length(avail_target_ids))
  names(comp_cells) <- avail_target_ids
  for (cid in avail_target_ids) {
    comp_cells[[cid]] <- build_complete_cell(cid, allow_empty_switchers = FALSE)
  }

  nt_ids_use <- if (use_or) non_target_ids else character(0L)
  nt_comp_cells <- NULL
  if (length(nt_ids_use) > 0L) {
    nt_comp_cells <- vector("list", length(nt_ids_use))
    names(nt_comp_cells) <- nt_ids_use
    for (cid in nt_ids_use) {
      nt_comp_cells[[cid]] <- build_complete_cell(cid, allow_empty_switchers = TRUE)
    }
  }

  or_blocks_use <- if (use_or) {
    or_blocks_fit
  } else {
    stats::setNames(avail_target_ids, avail_target_ids)
  }

  # Admit the complete kernel's numeric payload before entering it. In
  # particular, fitted/residual OR matrices span ALL training cells, even in
  # point-only mode. Count aliased treatment/OR designs once; count masks as
  # double entries conservatively. This bounds numeric payload, not peak RSS.
  admission_cells <- c(comp_cells, nt_comp_cells)
  k <- vapply(admission_cells, function(cl) as.double(ncol(cl$X_tilde)), numeric(1L))
  p <- k - 1.0
  k_target <- k[avail_target_ids]
  max_k <- max(k)
  n_target <- as.double(length(avail_target_ids))
  n_training <- as.double(length(admission_cells))
  use_logit <- gate$est_method %in% c("ipw", "dr")

  # Live cell inputs (dy, masks, X, Z), G_matrix, and cell/aggregate vectors.
  required_entries <- n * sum(4.0 + k + p) + n * n_target + 64.0 * n_target
  if (use_or) {
    block_p <- vapply(
      split(p, unname(or_blocks_use[all_fit_ids])), max, numeric(1L)
    )
    # Fitted/residual matrices; cell and pooled WLS summaries/factors, means,
    # intercepts and slopes, including the current block's merge/solve copies.
    required_entries <- required_entries + 2.0 * n * n_training +
      8.0 * sum((p + 1.0)^2) + 8.0 * sum((block_p + 1.0)^2)
  }
  if (use_logit) {
    # Retained eta/prob/odds/active/score weights and Hessian/QR diagnostics.
    required_entries <- required_entries + 5.0 * n * n_target +
      8.0 * sum(k_target^2) + 16.0 * sum(k_target)
  }
  if (isTRUE(compute_inffunc)) {
    # Two target IFs, four cell components, covariance input copy, four
    # covariance/vhat matrices and their crossproduct scratch.
    required_entries <- required_entries + 7.0 * n * n_target + 5.0 * n_target^2
    if (use_logit) {
      required_entries <- required_entries + 2.0 * n * sum(k_target)
    }
    if (use_or) {
      # Local intercept IFs span all training cells; slope scores/IFs span blocks.
      required_entries <- required_entries + n * n_training + 2.0 * n * sum(block_p)
    }
    if (all_declared_available) {
      # Aggregate components, selected/centered IF and G copies, joint IF,
      # joint covariance/vhat, and covariance conversion/crossproduct scratch.
      required_entries <- required_entries + n * (5.0 * n_target + 5.0) +
        3.0 * (n_target + 1.0)^2
    }
  }
  # One active nuisance fit/solve: risk/weighted/QR designs, centered slopes,
  # score/IF solve transposes and copies, factor scratch and unit work vectors.
  required_entries <- required_entries + n * (16.0 * max_k + 32.0) + 16.0 * max_k^2
  static_public_need(
    is.finite(required_entries) && required_entries <= ws_entries,
    "Complete-data baseline numeric outputs and scratch exceed `control$numeric_workspace_bytes`.",
    "cdid_gmm_bound_error",
    diagnostics = list(required_entries = required_entries, workspace_entries = ws_entries,
                       n_target_cells = n_target, n_training_cells = n_training)
  )

  fit_obj <- static_complete_fit(
    cells = comp_cells,
    or_blocks = or_blocks_use,
    method = gate$est_method,
    glm_control = glm_ctrl,
    rank_tol = rank_tol,
    unit_id = panel$unit_id,
    cluster = NULL,
    aggregate_cells = if (all_declared_available) avail_target_ids else NULL,
    aggregate_signs = if (all_declared_available) avail_sigma else NULL,
    non_target_or_cells = nt_comp_cells,
    materialize_if = isTRUE(compute_inffunc)
  )

  if (!is.null(fit_obj$phi_theta)) {
    rownames(fit_obj$phi_theta) <- panel$unit_key
  }
  if (!is.null(fit_obj$phi_tau)) {
    rownames(fit_obj$phi_tau) <- panel$unit_key
  }
  if (!is.null(fit_obj$aggregate) && !is.null(fit_obj$aggregate$phi)) {
    names(fit_obj$aggregate$phi) <- panel$unit_key
  }

  fit_obj$route <- "baseline"
  fit_obj$target_ids <- avail_target_ids
  fit_obj$non_target_ids <- nt_ids_use
  fit_obj$or_blocks <- if (use_or) or_blocks_fit else NULL
  fit_obj$cells <- all_cand_cells[avail_target_ids]
  fit_obj$all_cells <- all_cand_cells[all_fit_ids]
  fit_obj$diagnostics <- list(
    backend = "static_complete_fit",
    backend_status = "accepted_static_complete_fit_integration",
    method = fit_obj$method,
    requested_method = gate$est_method_requested,
    route = "baseline",
    selection_model = NULL,
    selection_pooling = NULL,
    imputation_pooling = NULL,
    imputation_scale = NULL,
    weighting_requested = gate$ctrl_resolved$gmm_weighting$requested,
    weighting_canonical = gate$ctrl_resolved$gmm_weighting$canonical,
    gmm_executed = FALSE,
    gmm_execution_status = "not_applicable_complete_data_baseline_reduction",
    n_units = n,
    n_zero_outcome_units = 0L,
    n_single_outcome_units = 0L,
    n_target_cells = length(avail_target_ids),
    n_all_cells = length(all_fit_ids),
    n_non_target_cells = length(nt_ids_use),
    n_links = NA_integer_,
    overid_df = NA_integer_,
    gmm_diagnostics = NULL
  )
  fit_obj
}

# Main entry point called by `att_gt()` when `dname` is non-NULL.
static_public_att_gt <- function(yname,
                                 tname,
                                 idname = NULL,
                                 dname,
                                 xformla = NULL,
                                 data,
                                 panel = TRUE,
                                 allow_unbalanced_panel = NULL,
                                 control_group = c("nevertreated", "notyettreated"),
                                 anticipation = 0,
                                 weightsname = NULL,
                                 fix_weights = NULL,
                                 alp = 0.05,
                                 bstrap = TRUE,
                                 cband = TRUE,
                                 biters = 1000,
                                 clustervars = NULL,
                                 est_method = "dr",
                                 base_period = "varying",
                                 faster_mode = TRUE,
                                 print_details = FALSE,
                                 pl = FALSE,
                                 cores = 1,
                                 compute_inffunc = TRUE,
                                 dots = list(),
                                 dots_eval = NULL,
                                 mode = NULL,
                                 frame = NULL,
                                 selection = "none",
                                 selection_model = "pair_logit",
                                 selection_pooling = "separate",
                                 selection_xformla = NULL,
                                 selection_universe = "F",
                                 selection_rho = NULL,
                                 selection_justification = NULL,
                                 selection_est_method = "ipw",
                                 imputation_pooling = "separate",
                                 gmm_weighting = c("optimal", "identity"),
                                 gmm_covariance = "tikhonov",
                                 pt_assumption = c("post", "all"),
                                 pt_start = NULL,
                                 control = NULL,
                                 outcome_pooling = "separate",
                                 imputation_scale = "per_period",
                                 switchers = "both",
                                 switch_times = NULL,
                                 env = parent.frame(),
                                 call = NULL) {
  if (is.null(call)) call <- sys.call()

  # 1. Scope gate before touching `data` (preserving laziness of inactive options)
  gate <- static_public_scope_gate(
    yname = yname,
    tname = tname,
    idname = idname,
    dname = dname,
    panel = panel,
    allow_unbalanced_panel = allow_unbalanced_panel,
    control_group = control_group,
    anticipation = anticipation,
    weightsname = weightsname,
    fix_weights = fix_weights,
    alp = alp,
    bstrap = bstrap,
    cband = cband,
    biters = biters,
    est_method = est_method,
    base_period = base_period,
    faster_mode = faster_mode,
    print_details = print_details,
    pl = pl,
    cores = cores,
    compute_inffunc = compute_inffunc,
    dots = dots,
    mode = mode,
    frame = frame,
    selection = selection,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    selection_universe = selection_universe,
    selection_rho = selection_rho,
    selection_justification = selection_justification,
    selection_est_method = selection_est_method,
    imputation_pooling = imputation_pooling,
    imputation_scale = imputation_scale,
    gmm_weighting = gmm_weighting,
    gmm_covariance = gmm_covariance,
    pt_assumption = pt_assumption,
    pt_start = pt_start,
    control = control,
    switchers = switchers,
    env = env
  )

  # 2. Evaluate active covariate formulas (keeping `selection_xformla` lazy when inactive)
  xformla_use <- if (is.null(xformla)) ~ 1 else xformla
  static_public_need(
    inherits(xformla_use, "formula") && length(xformla_use) == 2L,
    "`xformla` must be NULL or a one-sided formula (e.g. `~ x1 + x2`).",
    c("cdid_static_formula_error", "cdid_chained_request_input_error")
  )

  uses_q <- gate$route %in% c("weighting", "augmented")
  uses_j <- gate$route %in% c("imputation", "augmented")
  use_or <- gate$est_method %in% c("reg", "dr")

  sel_xformla_use <- NULL
  sel_formula_inherited <- NA
  if (uses_q) {
    if (is.null(selection_xformla)) {
      sel_xformla_use <- xformla_use
      sel_formula_inherited <- TRUE
    } else {
      static_public_need(
        inherits(selection_xformla, "formula") && length(selection_xformla) == 2L,
        "`selection_xformla` must be NULL or a one-sided formula.",
        c("cdid_static_formula_error", "cdid_chained_request_input_error")
      )
      sel_xformla_use <- selection_xformla
      sel_formula_inherited <- FALSE
    }
  }

  progress <- chained_progress_emitter(print_details)
  if (is.null(progress)) progress <- function(...) invisible(NULL)
  progress("static panel estimation started (method = ", gate$est_method, ", route = ", gate$route, ")")

  # 3. Build full-unit static panel data (retaining zero- and single-outcome units;
  #    requiring known binary treatment history D across all study dates)
  panel_obj <- static_panel_data(
    data = data,
    idname = idname,
    tname = tname,
    yname = yname,
    dname = dname
  )
  if (!all(panel_obj$d_known == 1L)) {
    static_public_abort(
      sprintf(
        "Static `dname` estimation requires a complete known binary treatment history D = (D_1, ..., D_T) across all %d study dates for every unit; %d unit-date treatment value(s) are missing or absent in `data`.",
        length(panel_obj$dates_num),
        as.integer(sum(panel_obj$d_known == 0L))
      ),
      c("cdid_static_missing_treatment_error", "cdid_static_panel_error", "cdid_static_unavailable")
    )
  }
  n <- panel_obj$n

  # 4. Validate and extract unit-constant clusters if supplied
  cluster_info <- static_public_extract_clusters(
    data = data,
    panel = panel_obj,
    idname = idname,
    clustervars = clustervars
  )

  # 5. Build candidate transition cells across all study dates t >= 2 and directions d in {0, 1}
  cand_grid <- expand.grid(d = c(0L, 1L), t = panel_obj$dates_num[-1L])
  cand_t_lbl <- panel_obj$date_labels[match(cand_grid$t, panel_obj$dates_num)]
  cand_cids <- sprintf("t%s_%s", cand_t_lbl, ifelse(cand_grid$d == 0L, "entry", "exit"))
  all_cand <- static_transition_cells(
    panel = panel_obj,
    t = cand_grid$t,
    d = cand_grid$d,
    cell_id = cand_cids
  )

  # 6. Resolve declared target set A
  targ_res <- static_public_resolve_targets(
    panel = panel_obj,
    all_cand = all_cand,
    switchers = gate$switchers,
    switch_times = switch_times
  )
  target_ids <- targ_res$target_ids
  J <- length(target_ids)

  # 7. Resolve `outcome_pooling` ONLY when `est_method %in% c("reg", "dr")`
  or_res <- if (use_or) {
    static_public_resolve_outcome_pooling(
      outcome_pooling = outcome_pooling,
      target_ids = target_ids,
      all_cand_cells = all_cand$cells
    )
  } else {
    list(
      mode = "inactive_ipw",
      or_blocks = NULL,
      all_ids = target_ids,
      non_target_ids = character(0L)
    )
  }

  # 8. Build all-positive-observed-count role graphs (`K_c^a >= 1`) and check support
  graph_support <- static_public_build_graphs_and_support(
    panel = panel_obj,
    target_ids = target_ids,
    non_target_ids = or_res$non_target_ids,
    all_cand_cells = all_cand$cells
  )
  available <- graph_support$available
  avail_target_ids <- graph_support$avail_target_ids
  all_declared_available <- all(available)

  # 9. Resolve workspace / work budgets from `control`
  prof <- gate$ctrl_resolved$profile
  ws_bytes <- as.double(prof$numeric_workspace_bytes)
  ws_entries <- floor(ws_bytes / 8.0)
  user_work_limit <- prof$counted_work_limit
  m_R_total <- sum(graph_support$summary$n_retained_edges)
  p_est <- length(all.vars(xformla_use)) + (if (uses_q && !is.null(sel_xformla_use)) length(all.vars(sel_xformla_use)) else 0L)
  est_work_ops <- as.double(n) * as.double(max(1L, m_R_total + J)) * as.double(max(1L, p_est + 1L)) * 64.0 +
    (if (isTRUE(compute_inffunc) && (isTRUE(bstrap) || isTRUE(cband))) as.double(biters) * as.double(cluster_info$n_clusters) * as.double(max(1L, J)) * 8.0 else 0.0)
  if (is.finite(user_work_limit) && est_work_ops > user_work_limit) {
    static_public_abort(
      sprintf(
        "Conservative preflight work estimate (%.0f ops) exceeds `control$counted_work_limit` (%.0f ops).",
        est_work_ops, user_work_limit
      ),
      c("cdid_static_work_limit_error", "cdid_chained_public_control_error", "cdid_gmm_bound_error"),
      diagnostics = list(estimated_work_ops = est_work_ops, counted_work_limit = user_work_limit)
    )
  }
  chunk_rows <- if (!is.null(gate$ctrl_resolved$chunk_rows)) {
    as.integer(max(1L, min(n, gate$ctrl_resolved$chunk_rows)))
  } else {
    as.integer(max(1L, min(
      n,
      gate$ctrl_resolved$chunk_rows_cap,
      floor(ws_entries / max(1.0, 2048.0 * (m_R_total + p_est + 1.0)))
    )))
  }

  glm_ctrl <- list(
    epsilon = gate$ctrl_resolved$glm$epsilon,
    maxit = as.integer(gate$ctrl_resolved$glm$maxit),
    trace = FALSE
  )
  rank_tol <- 1e-7

  # 10. Fit available target cells via complete-data baseline kernel when
  #     all outcomes are observed and `route == "baseline"`, or via `static_missing_fit`
  #     when outcomes are missing or a conditional-MTAR correction is requested.
  is_complete_baseline <- all(panel_obj$O == 1L) && identical(gate$route, "baseline")
  fit_obj <- NULL

  if (length(avail_target_ids) >= 1L) {
    all_fit_ids <- c(avail_target_ids, or_res$non_target_ids)
    or_blocks_fit <- if (use_or) or_res$or_blocks[all_fit_ids] else NULL
    avail_sigma <- vapply(
      all_cand$cells[avail_target_ids],
      function(cl) as.double(cl$sigma),
      numeric(1L)
    )
    names(avail_sigma) <- avail_target_ids

    if (is_complete_baseline) {
      fit_obj <- static_public_fit_complete_baseline(
        panel = panel_obj,
        avail_target_ids = avail_target_ids,
        non_target_ids = or_res$non_target_ids,
        all_cand_cells = all_cand$cells,
        or_blocks_fit = or_blocks_fit,
        xformla_use = xformla_use,
        use_or = use_or,
        gate = gate,
        glm_ctrl = glm_ctrl,
        rank_tol = rank_tol,
        ws_entries = ws_entries,
        all_declared_available = all_declared_available,
        avail_sigma = avail_sigma,
        compute_inffunc = isTRUE(compute_inffunc)
      )
    } else {
      gmm_bounds <- gmm_bounds_control(
        max_matrix_entries = ws_entries,
        max_batch_rows = chunk_rows,
        max_targets = max(1L, length(avail_target_ids)),
        max_open_clusters = cluster_info$n_clusters,
        allow_materialize_if = isTRUE(compute_inffunc),
        max_materialized_if_entries = min(as.double(n) * max(1.0, as.double(length(avail_target_ids))), ws_entries)
      )
      gmm_acc <- static_graph_gmm_accuracy_control(
        max_target_error_abs = chained_public_staging_envelope(),
        max_target_error_se = NA_real_,
        max_orthogonality_defect = gate$ctrl_resolved$gmm_guard$max_orthogonality_defect,
        min_certified_lambda_rel = gate$ctrl_resolved$gmm_guard$min_certified_lambda_rel,
        max_factor_residual_rel = gate$ctrl_resolved$gmm_guard$max_factor_residual_rel,
        min_rcond_R = gate$ctrl_resolved$gmm_guard$min_rcond_R,
        min_rcond_U = gate$ctrl_resolved$gmm_guard$min_rcond_U,
        min_sigma_rel_A = gate$ctrl_resolved$gmm_guard$min_sigma_rel_A,
        max_kkt_residual_rel = gate$ctrl_resolved$gmm_guard$max_kkt_residual_rel,
        max_loading_residual_rel = gate$ctrl_resolved$gmm_guard$max_loading_residual_rel,
        require_matching_scale = FALSE,
        max_loading_error_abs = NULL,
        covariance_strategy = gate$ctrl_resolved$covariance_strategy
      )

      fit_obj <- static_missing_fit(
        panel = panel_obj,
        target_cells = avail_target_ids,
        all_cells = all_cand$cells[all_fit_ids],
        or_blocks = or_blocks_fit,
        role_graphs = graph_support$role_graphs_fit,
        treatment_formula = if (gate$est_method %in% c("ipw", "dr")) xformla_use else NULL,
        outcome_formula = if (use_or) xformla_use else NULL,
        selection_formula = if (uses_q) sel_xformla_use else NULL,
        imputation_formula = if (uses_j) xformla_use else NULL,
        method = gate$est_method,
        route = gate$route,
        selection_model = gate$selection_model,
        selection_pooling = gate$selection_pooling,
        imputation_pooling = gate$imputation_pooling,
        imputation_scale = gate$imputation_scale,
        weighting = gate$ctrl_resolved$gmm_weighting$requested,
        glm_control = glm_ctrl,
        rank_tol = rank_tol,
        bounds = gmm_bounds,
        accuracy = gmm_acc,
        cluster = NULL,
        centered = FALSE,
        aggregate_cells = if (all_declared_available) avail_target_ids else NULL,
        aggregate_signs = if (all_declared_available) avail_sigma else NULL,
        materialize_if = isTRUE(compute_inffunc)
      )
    }
  }

  # 11. Scatter fitted values into full declared target vector A (length J)
  decl_cells <- all_cand$cells[target_ids]
  t_vec <- vapply(decl_cells, `[[`, numeric(1L), "t_num")
  t_prev_vec <- vapply(decl_cells, `[[`, numeric(1L), "t_prev_num")
  d_vec <- vapply(decl_cells, `[[`, integer(1L), "d")
  dir_vec <- ifelse(d_vec == 0L, "entry", "exit")
  sigma_vec <- vapply(decl_cells, `[[`, numeric(1L), "sigma")
  n_risk_vec <- vapply(decl_cells, `[[`, integer(1L), "n_risk")
  n_sw_vec <- vapply(decl_cells, `[[`, integer(1L), "n_switchers")
  n_st_vec <- vapply(decl_cells, `[[`, integer(1L), "n_stayers")

  att_vec <- stats::setNames(rep(NA_real_, J), target_ids)
  theta_vec <- stats::setNames(rep(NA_real_, J), target_ids)
  se_vec <- stats::setNames(rep(NA_real_, J), target_ids)
  an_se_vec <- stats::setNames(rep(NA_real_, J), target_ids)
  boot_se_vec <- stats::setNames(rep(NA_real_, J), target_ids)

  if (!is.null(fit_obj)) {
    att_vec[avail_target_ids] <- as.double(fit_obj$tau_hat[avail_target_ids])
    theta_vec[avail_target_ids] <- as.double(fit_obj$theta_hat[avail_target_ids])
  }

  # 12. Inference (`compute_inffunc = TRUE` vs `compute_inffunc = FALSE`)
  point_only <- !isTRUE(compute_inffunc)
  flags <- ipw_pipeline_resolve_bootstrap_flags(bstrap, cband, endpoint = "cells")
  effective_bstrap <- !point_only && flags$bootstrap
  effective_cband <- !point_only && flags$simultaneous
  zcrit <- stats::qnorm(1.0 - alp / 2.0)

  inffunc_full <- NULL
  V_analytical_full <- NULL
  cband_crit_val <- NA_real_
  cband_status <- if (point_only) "unavailable_compute_inffunc_false" else "not_requested"
  inf_res <- NULL

  if (!point_only) {
    static_public_need(
      as.double(n) * as.double(J) <= ws_entries && as.double(J)^2 <= ws_entries,
      "Declared static IF/covariance storage exceeds `control$numeric_workspace_bytes`.",
      c("cdid_static_inference_bound_error", "cdid_multiplier_budget_error")
    )
    inffunc_full <- matrix(NA_real_, nrow = n, ncol = J, dimnames = list(panel_obj$unit_key, target_ids))
    V_analytical_full <- matrix(NA_real_, nrow = J, ncol = J, dimnames = list(target_ids, target_ids))
    if (!is.null(fit_obj) && length(avail_target_ids) >= 1L) {
      phi_avail <- fit_obj$phi_tau[, avail_target_ids, drop = FALSE]
      inffunc_full[, avail_target_ids] <- phi_avail

      agg_est_in <- if (all_declared_available && !is.null(fit_obj$aggregate)) fit_obj$aggregate$delta_hat else NULL
      agg_if_in <- if (all_declared_available && !is.null(fit_obj$aggregate)) fit_obj$aggregate$phi else NULL

      inf_res <- static_public_inference(
        estimates = att_vec[avail_target_ids],
        if_matrix = phi_avail,
        unit_ids = panel_obj$unit_id,
        cluster_vec = cluster_info$cluster_vec,
        alp = alp,
        bstrap = effective_bstrap,
        cband = effective_cband,
        biters = biters,
        chunk_rows = chunk_rows,
        workspace_entries = ws_entries,
        workspace_bytes = ws_bytes,
        workers = gate$parallel_plan$workers,
        fixed_signs = NULL,
        agg_estimate = agg_est_in,
        agg_if = agg_if_in,
        all_declared_available = all_declared_available
      )

      V_analytical_full[avail_target_ids, avail_target_ids] <- inf_res$V_analytical
      an_se_vec[avail_target_ids] <- inf_res$analytical_se
      boot_se_vec[avail_target_ids] <- inf_res$bootstrap_se
      se_vec[avail_target_ids] <- inf_res$se
      cband_crit_val <- inf_res$cband_crit_val
      cband_status <- inf_res$cband_status
    } else if (!all_declared_available && isTRUE(effective_cband)) {
      cband_status <- "unavailable_target_cell"
    }
  }

  ci_lower <- att_vec - zcrit * se_vec
  ci_upper <- att_vec + zcrit * se_vec
  cband_lower <- if (is.finite(cband_crit_val)) att_vec - cband_crit_val * se_vec else stats::setNames(rep(NA_real_, J), target_ids)
  cband_upper <- if (is.finite(cband_crit_val)) att_vec + cband_crit_val * se_vec else stats::setNames(rep(NA_real_, J), target_ids)
  crit_used <- if (isTRUE(effective_cband)) cband_crit_val else zcrit

  # 13. Build switcher-frequency aggregate record (`eq:if-agg`)
  if (all_declared_available && !is.null(fit_obj) && !is.null(fit_obj$aggregate)) {
    agg_obj <- fit_obj$aggregate
    agg_est <- as.double(agg_obj$delta_hat)
    agg_se <- if (!point_only && !is.null(inf_res)) as.double(inf_res$agg_se) else NA_real_
    agg_an_se <- if (!point_only && !is.null(inf_res)) as.double(inf_res$agg_analytical_se) else NA_real_
    agg_bt_se <- if (!point_only && !is.null(inf_res)) as.double(inf_res$agg_bootstrap_se) else NA_real_
    agg_record <- list(
      available = TRUE,
      status = "ok",
      reason = NA_character_,
      cells = target_ids,
      signs = stats::setNames(as.double(sigma_vec), target_ids),
      weights = agg_obj$omega_hat,
      n_switchers_by_cell = stats::setNames(as.integer(n_sw_vec), target_ids),
      total_switchers = as.integer(agg_obj$total_switchers),
      Pi_hat = as.double(agg_obj$Pi_hat),
      estimate = agg_est,
      se = agg_se,
      analytical_se = agg_an_se,
      bootstrap_se = agg_bt_se,
      crit_val = zcrit,
      ci_lower = agg_est - zcrit * agg_se,
      ci_upper = agg_est + zcrit * agg_se,
      inffunc = if (!point_only) agg_obj$phi else NULL,
      scope = "static switcher-frequency weighted average across the declared transition set A with random switch-share influence function (eq:if-agg); not a cohort or dynamic event-time aggregation"
    )
  } else {
    agg_record <- list(
      available = FALSE,
      status = "unavailable_target_cell",
      reason = "One or more declared target cells in A are unavailable; the static aggregate is never silently renormalized over a subset of surviving cells.",
      cells = target_ids,
      signs = stats::setNames(as.double(sigma_vec), target_ids),
      weights = stats::setNames(rep(NA_real_, J), target_ids),
      n_switchers_by_cell = stats::setNames(as.integer(n_sw_vec), target_ids),
      total_switchers = as.integer(sum(n_sw_vec)),
      Pi_hat = as.double(sum(n_sw_vec)) / as.double(n),
      estimate = NA_real_,
      se = NA_real_,
      analytical_se = NA_real_,
      bootstrap_se = NA_real_,
      crit_val = zcrit,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      inffunc = NULL,
      scope = "static switcher-frequency weighted average across the declared transition set A with random switch-share influence function (eq:if-agg); not a cohort or dynamic event-time aggregation"
    )
  }

  # 14. Assemble honest `cdid_static_att_gt` result object
  new_cdid_static_att_gt(
    target_id = target_ids,
    t = t_vec,
    t_prev = t_prev_vec,
    d = d_vec,
    direction = dir_vec,
    sigma = sigma_vec,
    available = available,
    support_status = graph_support$support_status,
    support_reason = graph_support$support_reason,
    n_risk = n_risk_vec,
    n_switchers = n_sw_vec,
    n_stayers = n_st_vec,
    att = att_vec,
    theta_hat = theta_vec,
    se = se_vec,
    analytical_se = an_se_vec,
    bootstrap_se = boot_se_vec,
    V_analytical = V_analytical_full,
    c = crit_used,
    pointwise_crit_val = zcrit,
    cband_crit_val = cband_crit_val,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    cband_lower = cband_lower,
    cband_upper = cband_upper,
    alp = alp,
    bstrap = effective_bstrap,
    cband = effective_cband,
    biters = as.integer(biters),
    n = n,
    unit_id = panel_obj$unit_id,
    dates = panel_obj$dates,
    dates_num = panel_obj$dates_num,
    inffunc = inffunc_full,
    aggregate = agg_record,
    graphs = list(
      rule = "all_positive_observed_count_pairs_per_role_no_pruning",
      summary = graph_support$summary,
      by_cell_role = graph_support$by_cell_role
    ),
    request = list(
      call = call,
      yname = yname,
      tname = tname,
      idname = idname,
      dname = dname,
      xformla = xformla_use,
      est_method_requested = gate$est_method_requested,
      est_method = gate$est_method,
      route = gate$route,
      selection_requested = gate$selection_requested,
      selection = gate$selection_canonical,
      selection_est_method_requested = gate$selection_est_method_requested,
      selection_est_method = gate$selection_est_method,
      selection_model = gate$selection_model,
      selection_pooling = gate$selection_pooling,
      selection_xformla = sel_xformla_use,
      selection_formula_inherited = sel_formula_inherited,
      imputation_pooling = gate$imputation_pooling,
      imputation_scale = gate$imputation_scale,
      outcome_pooling_mode = or_res$mode,
      or_blocks = or_res$or_blocks,
      non_target_or_cells = or_res$non_target_ids,
      switchers = gate$switchers,
      switch_times = switch_times,
      target_discovery = targ_res$discovery_mode,
      gmm_weighting = gate$ctrl_resolved$gmm_weighting,
      gmm_covariance = gate$ctrl_resolved$gmm_covariance,
      compute_inffunc = isTRUE(compute_inffunc),
      requested_bstrap = isTRUE(bstrap),
      requested_cband = isTRUE(cband),
      inference_flag_law = flags$law,
      effective_bstrap = effective_bstrap,
      effective_cband = effective_cband,
      alp = alp,
      biters = as.integer(biters),
      clustervars = clustervars,
      clustername = cluster_info$clustername,
      faster_mode = isTRUE(faster_mode),
      print_details = isTRUE(print_details),
      pl = isTRUE(pl),
      cores = as.integer(cores),
      execution = gate$parallel_plan,
      control = list(
        profile = gate$ctrl_resolved$profile,
        profile_status = gate$ctrl_resolved$profile_status,
        overrides = gate$ctrl_resolved$overrides,
        chunk_rows = chunk_rows,
        estimated_work_ops = est_work_ops
      )
    ),
    inference = list(
      compute_inffunc = isTRUE(compute_inffunc),
      se_reported = if (point_only) "none_compute_inffunc_false" else if (effective_bstrap) "multiplier_bootstrap_iqr" else "analytical",
      band_reported = if (point_only) "none_compute_inffunc_false" else if (effective_cband && is.finite(cband_crit_val)) "multiplier_sup_t" else if (effective_cband) paste0("unavailable (", cband_status, ")") else "pointwise_normal",
      cband_status = cband_status,
      clustered = cluster_info$clustered,
      clustername = cluster_info$clustername,
      n_clusters = cluster_info$n_clusters,
      n_units = n,
      normalization = "1_over_n_centered_cluster_sums_vcov_1_over_n2_bootstrap_se_sqrtG_over_n",
      quantile_type = 1L,
      qualifications = paste(
        "Asymptotic validity requires compatible link moments, regular interior nuisance roots",
        "with inactive probability clipping in population, fixed finite support (T, |A|, m_R),",
        "and many independent units/clusters (G -> Inf); clustering changes inference only."
      ),
      multiplier_diagnostics = if (!is.null(inf_res) && !is.null(inf_res$multiplier)) inf_res$multiplier$diagnostics else NULL
    ),
    diagnostics = list(
      benchmark = "source_static_no_carryover_known_binary_D",
      backend = if (is_complete_baseline) "static_complete_fit" else "static_missing_fit",
      backend_status = if (is_complete_baseline) {
        "accepted_static_complete_fit_integration"
      } else {
        "accepted_static_missing_fit_integration"
      },
      gmm_executed = !is_complete_baseline && !is.null(fit_obj),
      n_units = n,
      n_dates = length(panel_obj$dates_num),
      n_zero_outcome_units = as.integer(sum(rowSums(panel_obj$O) == 0L)),
      n_single_outcome_units = as.integer(sum(rowSums(panel_obj$O) == 1L)),
      n_repeat_outcome_units = as.integer(sum(rowSums(panel_obj$O) >= 2L)),
      n_declared_targets = J,
      n_available_targets = as.integer(sum(available)),
      n_unavailable_targets = as.integer(sum(!available)),
      n_non_target_or_cells = length(or_res$non_target_ids),
      n_total_links = as.integer(m_R_total),
      backend_diagnostics = if (!is.null(fit_obj)) fit_obj$diagnostics else NULL
    ),
    fit = fit_obj
  )
}
