# The regression branch uses within-unit least squares. Treatment IPW and DR
# combine normalized changes on cohort-standardized calendar paths.
chained_cohort_time_fe <- function(data, frame, yname, tname, idname, gname,
                                   xformla, weightsname, clustervars, gate,
                                   est_method, base_period, anticipation, alp,
                                   bstrap, cband, biters, gmm_weighting,
                                   control, call, selection_xformla = NULL,
                                   selection_universe = "F", selection_rho = NULL,
                                   selection_justification = NULL, keep_nuisance = FALSE) {
  need <- function(ok, message, cls = "cdid_etwfe_scope_error") {
    chained_request_need(ok, message, cls)
  }
  need(!is.object(control) && !anyNA(names(control)) && !anyDuplicated(names(control)) &&
         all(names(control) %in% c("preset", "numeric_workspace_bytes", "workers", "chunk_rows")),
       "The cohort_time_fe preset accepts numeric_workspace_bytes, workers and chunk_rows as advanced controls.")
  need(is.character(est_method) && length(est_method) == 1L &&
         est_method %in% c("reg", "ipw", "dr"),
       "Use a built-in treatment estimator with the cohort_time_fe preset.")
  regression <- identical(est_method, "reg")
  corrected <- !identical(gate$selection_canonical, "none")
  if (!is.null(xformla)) {
    need(inherits(xformla, "formula") && length(xformla) == 2L &&
           !length(attr(stats::terms(xformla), "offset")),
         "The cohort_time_fe preset requires a one-sided xformla without offsets.")
  }
  need(is.null(weightsname) && is.null(gate$fix_weights),
       "The linear cohort_time_fe preset requires unit sampling weights.")
  need(anticipation == 0 && identical(base_period, "universal"),
       "Use anticipation = 0 and base_period = \"universal\" for the cohort_time_fe preset.")
  if (regression) {
    need(!gate$cross_pair && identical(gate$gmm_covariance, "tikhonov") &&
           (identical(gmm_weighting, c("optimal", "identity", "within")) ||
              identical(gate$gmm_weighting, "within")),
         "The cohort_time_fe preset uses within-unit least squares. Leave cross_pair and gmm_covariance at their defaults and omit gmm_weighting or set it to \"within\".")
  } else {
    need(!gate$cross_pair && gate$gmm_weighting %in% c("optimal", "identity"),
         "Treatment IPW and DR with cohort_time_fe use optimal or identity moment weighting. Leave cross_pair at its default.")
  }
  nyt <- identical(gate$control_group, "notyettreated")
  need(if (nyt) identical(gate$pt_assumption, "all") else identical(gate$pt_assumption, "post"),
       if (nyt) "Not-yet-treated cohort_time_fe requires pt_assumption = \"all\" and pt_start at the first study date."
       else "Never-treated cohort_time_fe uses pt_assumption = \"post\" and estimates separate nonreference pre-treatment deviations.")

  control$preset <- NULL
  resolved <- if (regression) {
    chained_public_resolve_control(control, gmm_weighting = "within")
  } else {
    chained_public_resolve_control(
      control, gmm_weighting = gate$gmm_weighting_requested,
      gmm_covariance = gate$gmm_covariance
    )
  }
  plan <- chained_request_parallel_plan(gate$pl, gate$cores, resolved$profile$workers)
  progress <- chained_progress_emitter(gate$print_details)
  if (!is.null(progress)) progress("cohort_time_fe: preparing the observed panel")
  prep <- if (corrected) {
    if (!identical(gate$selection_est_method, "imputation") && !is.null(selection_xformla)) {
      need(inherits(selection_xformla, "formula") && length(selection_xformla) == 2L &&
             !length(attr(stats::terms(selection_xformla), "offset")),
           "The cohort_time_fe preset requires a one-sided selection_xformla without offsets.")
    }
    chained_request_prepare_inputs(
      data = data, frame = frame, yname = yname, tname = tname, idname = idname,
      gname = gname, xformla = xformla, clustervars = clustervars,
      selection = gate$selection, selection_model = gate$selection_model,
      selection_pooling = gate$selection_pooling, selection_xformla = selection_xformla,
      selection_universe = selection_universe, selection_rho = selection_rho,
      selection_justification = selection_justification, or_resource = resolved$profile,
      selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
      outcome_pooling = gate$outcome_pooling, imputation_est_method = est_method,
      augmented_est_method = est_method, retain_panel_index = TRUE, full_outcome_frame = TRUE
    )
  } else {
    chained_request_prepare_inputs(
      data = data, frame = frame, yname = yname, tname = tname, idname = idname,
      gname = gname, xformla = xformla, clustervars = clustervars, retain_panel_index = TRUE
    )
  }
  calendar <- prep$calendar
  cohort <- prep$frame_g
  K <- length(calendar)
  treated <- is.finite(cohort)
  need(K >= 2L && all(calendar > 0 & calendar == floor(calendar)) && all(diff(calendar) == 1),
       "The linear cohort_time_fe preset requires consecutive positive integer study dates.")
  need(!nyt || gate$pt_start == calendar[1L], "Set pt_start to the first study date for the not-yet-treated cohort_time_fe preset.")
  need(!anyNA(cohort) && any(treated) && any(!treated) && all(cohort[treated] %in% calendar[-1L]),
       "The cohort_time_fe preset requires known cohorts, a never-treated cohort and treated cohorts with adoption dates after the first study date.")
  need(is.null(prep$reference_source),
       "The cohort_time_fe preset currently requires time-invariant covariates.")
  covariates <- prep$frame_Z[, colnames(prep$frame_Z) != "(Intercept)", drop = FALSE]
  has_covariates <- ncol(covariates) > 0L
  need(!has_covariates || identical(est_method, "ipw") ||
         identical(gate$outcome_pooling, "separate"),
       "Leave outcome_pooling at its default for cohort_time_fe with covariates. The preset specifies its own common calendar outcome model.")
  cohort_index <- match(cohort, calendar, nomatch = 0L)
  n <- prep$n_frame
  index <- prep$panel_index
  units <- prep$units
  eligible <- units$n_obs >= 2L
  need(all(is.finite(index$y)), "Observed outcomes must be finite.", "cdid_etwfe_input_error")
  if (regression && !corrected && gate$compute_inffunc && !is.null(prep$cluster_full)) {
    need(length(unique(prep$cluster_full[eligible])) >= 2L,
         "Cluster inference requires at least two clusters with repeat-observed units.", "cdid_etwfe_cluster_error")
  }
  if (is.null(prep$cluster_full)) units$cluster <- NULL
  units$eligible <- eligible

  # No design matrix with one row per observation and one column per effect.
  # Count only the largest numeric payloads of the slot system and retained IF.
  H <- as.double(length(unique(cohort_index[treated])))
  n_targets <- max(1, H * (K - 1))
  bytes <- resolved$profile$numeric_workspace_bytes
  entries <- floor(bytes / 8)
  need(4 * (K + n_targets)^2 + (if (gate$compute_inffunc) 3 * n * n_targets else 0) <= entries,
       "The cohort_time_fe factorization and requested influence functions exceed numeric_workspace_bytes.",
       "cdid_etwfe_workspace_error")
  chunk <- if (is.null(resolved$chunk_rows)) min(n, 256L) else min(n, resolved$chunk_rows)
  report <- if (keep_nuisance) {
    list(calendar = calendar, clusters = prep$cluster_full, cluster_var = prep$clustervar,
         compute_inffunc = gate$compute_inffunc, chunk_rows = chunk, workspace_entries = entries,
         selection_formula = selection_xformla, outcome_formula = xformla)
  }
  observation <- obs <- NULL
  fit <- if (corrected) {
    catalog <- etwfe_corrected_catalog(index, cohort_index, K, entries)
    geometry <- if (regression) etwfe_corrected_regression_geometry(catalog, gate$control_group, entries)
    support <- if (!regression) etwfe_corrected_treatment_layout(catalog, gate$control_group, est_method, entries)
    fit_cells <- if (regression) geometry$use else support$used_cells
    sel <- if (!identical(gate$selection_est_method, "imputation")) {
      c(prep$cp04[c("Z", "U", "rho", "universe")],
        list(model = gate$selection_model, pooling = gate$selection_pooling))
    }
    glm_control <- list(epsilon = resolved$glm$epsilon, maxit = as.integer(resolved$glm$maxit), trace = FALSE)
    if (!is.null(progress)) progress("cohort_time_fe: fitting observation models on supported cohort pairs")
    obs <- if (length(fit_cells)) {
      fitted <- etwfe_observation_fit(
        catalog$Y, cohort_index, catalog$cells[fit_cells, , drop = FALSE], prep$frame_Z, sel,
        gate$selection_est_method, gate$imputation_pooling, glm_control,
        rank_tol = 1e-7, chunk_rows = chunk, workspace_entries = entries, report = report
      )
      cell_map <- match(seq_len(nrow(catalog$cells)), fit_cells)
      fitted$read <- local({
        read_fit <- fitted$read
        function(cell, rows, gradient = FALSE) read_fit(cell_map[cell], rows, gradient)
      })
      fitted
    } else list(n_parameters = 0L, diagnostics = list(status = "no_required_pairs"))
    obs$local_width <- min(obs$n_parameters, ncol(prep$frame_Z) + 4L +
      if (is.null(sel)) 0L else 2L * ncol(sel$Z))
    observation <- list(
      method = gate$selection_est_method, model = gate$selection_model,
      selection_pooling = gate$selection_pooling, imputation_pooling = gate$imputation_pooling,
      universe = selection_universe, ratio_column = selection_rho,
      justification = selection_justification, support = catalog$cells[, c("cohort", "s", "t", "observed")],
      omitted_pairs = catalog$omitted_pairs, used_pairs = fit_cells,
      model_diagnostics = obs$diagnostics
    )
    features <- etwfe_corrected_features(covariates)
    if (regression) {
      etwfe_corrected_regression(catalog, features, gate$control_group, obs,
                                 gate$compute_inffunc, chunk, entries, geometry = geometry,
                                 report = report, Z = covariates)
    } else {
      scored <- gate$compute_inffunc || identical(gate$gmm_weighting, "optimal")
      moments <- etwfe_corrected_treatment(
        catalog, cohort_index, features, covariates, obs, gate$control_group,
        est_method, scored, glm_control, chunk, entries, report = report, support = support
      )
      etwfe_treatment_core(index, cohort_index, n, calendar, covariates, gate,
                           est_method, resolved, entries, plan$workers, biters,
                           prep$cluster_full, progress, moments = moments, report = report)
    }
  } else if (!regression) {
    etwfe_treatment_core(index, cohort_index, n, calendar, covariates, gate,
                         est_method, resolved, entries, plan$workers, biters,
                         prep$cluster_full, progress, report = report)
  } else if (has_covariates) {
    etwfe_covariate_core(
      index$y, index$unit, index$period, cohort_index, n, K,
      gate$control_group, covariates, compute_inffunc = gate$compute_inffunc,
      workspace_entries = entries, report = report
    )
  } else {
    etwfe_linear_core(
      index$y, index$unit, index$period, cohort_index, n, K,
      gate$control_group, gate$compute_inffunc, report = report
    )
  }

  requested <- fit$targets
  ids <- paste(calendar[requested$cohort], calendar[requested$period], sep = ":")
  available <- is.finite(requested$att)
  tm <- data.frame(
    g = calendar[requested$cohort], time = calendar[requested$period],
    target_id = ids, role = ifelse(requested$is_pre, "pre", "post"),
    ref = calendar[requested$cohort - 1L], base_t = calendar[requested$cohort - 1L],
    available = available, support_status = requested$status,
    support_reason = ifelse(available, "",
      ifelse(requested$status == "no_support", "There is no observed outcome in this cohort-date cell.",
             "The cohort-date and common date levels are not linked by observed differences."))
  )
  estimate <- stats::setNames(requested$att, ids)
  if (!is.null(fit$inffunc)) colnames(fit$inffunc) <- ids
  dims <- list(n_frame = n, n_clusters = if (is.null(prep$cluster_full)) n else length(unique(prep$cluster_full)),
               J_bound = nrow(tm))
  budget <- chained_public_multiplier_budget(dims, biters, chunk, entries, bytes, workers = plan$workers)
  inf <- if (gate$compute_inffunc) {
    ipw_pipeline_supported_inference(
      estimate, fit$inffunc, available, cluster_full = prep$cluster_full, alp = alp,
      cband = cband, biters = biters, multiplier_budget = budget, inference_accuracy = NULL,
      eligible_full = eligible, bstrap = bstrap, caller_rng = TRUE, target_meta = tm
    )
  } else {
    ipw_pipeline_point_only_inference(estimate, available, n, prep$cluster_full, eligible, alp)
  }
  reporting <- c(list(target_meta = tm, if_full = fit$inffunc, inference = inf,
                      base_period = base_period, available = available),
                 inf[c("estimates", "se", "z_stat", "p_value", "ci_lower", "ci_upper",
                       "cband_lower", "cband_upper", "analytical_se", "cband_crit_val", "support")])
  pretest <- if (gate$compute_inffunc && !nyt) {
    chained_pretest_record(reporting, max_matrix_entries = entries)
  } else if (nyt) {
    list(status = "not_applicable", reason = "The not-yet-treated model imposes common pre-treatment changes")
  } else {
    list(status = "not_computed", reason = "compute_inffunc = FALSE returns point estimates only")
  }

  mp <- chained_did_MP(
    reporting, prep$frame_ids, call,
    params = list(yname = yname, tname = tname, idname = idname, gname = gname,
                  alp = alp, xformla = prep$xformla, control_group = gate$control_group,
                  anticipation = anticipation, est_method = est_method, base_period = base_period,
                  tlist = calendar, glist = sort(unique(cohort[treated])),
                  bstrap = gate$bstrap, cband = gate$cband, biters = biters,
                  clustervars = prep$clustervar, compute_inffunc = gate$compute_inffunc,
                  print_details = gate$print_details, pl = gate$pl, cores = gate$cores),
    units = units, pretest = if (identical(pretest$status, "computed")) pretest
  )
  label <- if (regression) "linear cohort-by-date fixed-effects regression" else
    paste("cohort-standardized calendar", if (est_method == "ipw") "IPW" else "doubly robust")
  mp$DIDparams$cdid$estimator <- mp$chained$estimator <- label
  mp$chained$pretest <- pretest
  mp$chained$cohort_time_fe <- list(
    outcome_model = if (identical(est_method, "ipw")) "none" else "linear",
    fixed_effects = c("unit", "date"),
    control_group = gate$control_group, sample_rows = length(index$y),
    covariance = "original-unit score sandwich, no degrees-of-freedom multiplier",
    pt_start = gate$pt_start
  )
  if (has_covariates || !regression) {
    mp$chained$cohort_time_fe$standardization <- "full-frame within-cohort covariate means"
  }
  if (!regression) mp$chained$cohort_time_fe$treatment_moments <- fit$diagnostics
  if (!is.null(observation)) mp$chained$cohort_time_fe$observation <- observation
  mp$chained$aggregation_state <- list(
    calendar = calendar, cohort_grid = sort(unique(cohort)),
    multiplier_budget = budget, inference_accuracy = NULL, enclosure = NULL,
    btype = "rademacher", caller_rng = TRUE, compute_inffunc = gate$compute_inffunc,
    public = list(dims = dims, chunk_rows = chunk, entries = entries, bytes = bytes),
    parallel = plan,
    scope = "Existing cdid full-frame cohort-share aggregation. The fit is not repeated."
  )
  if (keep_nuisance) {
    mp$chained$nuisance_models <- etwfe_nuisance_capture(
      est_method, gate$selection_est_method, gate$selection_canonical, obs$nuisance, fit$nuisance
    )
  }
  mp
}
