# Private CP02 development request path: familiar att_gt()/aggte() requests ->
# accepted chained IPW-GMM pipeline -> reviewed did-compatible result subclasses.
#
#   chained_att_gt_request()   public-style att_gt metadata + EXPLICIT development
#                              controls -> cdid_MP (one bounded fit)
#   chained_aggte_request()    public-style aggte request on that cdid_MP ->
#                              cdid_AGGTEobj (no refit; released selectors)
#
# Public candidate binding (October 17, 2026) adds automatic controls and the
# approved D23/D21 rules. The private development entry retains its explicit
# fourteen-control object and root witness. D27 adopts finite resource defaults;
# estimator, sample and inference laws are unchanged. The
# executor validates the public input ONCE, aligns every unit attribute by
# identity, binds the public never-treated code 0 to the internal Inf
# explicitly, runs the accepted pipeline, and maps the already computed outputs
# through the reviewed adapters (R/chained_did_results.R). Composition note:
# docs/theory/cp02_public_request_path.md.

chained_request_abort <- function(message, class = character()) {
  cdid_abort(message, unique(c(class, "cdid_chained_request_error")))
}

chained_request_need <- function(ok, message, class = "cdid_chained_request_input_error") {
  if (!isTRUE(ok)) chained_request_abort(message, class)
  invisible(TRUE)
}

chained_request_chr1 <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
chained_request_lgl1 <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
chained_request_whole <- function(x, min = 0) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= min && x == floor(x)

# Request-level refusal context (diagnostics only; October 17, 2026). A typed
# pipeline refusal raised inside the one accepted fit is re-signalled unchanged
# (class, message, every existing diagnostic field) with the scalar request
# facts the pipeline cannot see: entry, resolved numerical policy and resource
# profile, the prepared input's counts and the fixed request settings. No owner,
# data, environment or new option is retained; nothing is written to disk.
chained_request_refusal_context <- function(entry, prep, settings, policy) {
  list(entry = entry,
    input = list(n_frame = prep$n_frame, n_frame_only = prep$n_frame_only, n_panel_rows = nrow(prep$data_internal),
      n_calendar = length(prep$calendar), calendar = prep$calendar, cohort_universe = prep$cohort_universe,
      panel_cohorts = prep$panel_cohorts, target_cohorts = prep$target_cohorts,
      cohorts_without_post_dates = prep$cohorts_without_post_dates, n_covariate_columns = ncol(prep$frame_Z),
      clustered = !is.null(prep$cluster_full), n_units_recoded = prep$recoding$n_units_recoded,
      outcome_design = if (is.null(prep$imputation_design)) prep$or_basis else prep$imputation_design),
    settings = settings, numerical_policy = policy,
    meaning = "request facts only; the refusal class, message and stage diagnostics are those of the accepted pipeline")
}

# ------------------------------------------------------------------------------
# CP02 scope gate shared by the public chained path and the development entry
# ------------------------------------------------------------------------------
# Every refusal is explicit and typed, and happens before any data promise is
# forced. The final public default est_method = "dr" is refused for
# chained estimation at CP02 with the familiar `cdid_chained_unavailable` class;
# IPW is never substituted.
chained_request_scope_gate <- function(mode, panel, allow_unbalanced_panel, est_method, selection,
                                       selection_model = "pair_logit", selection_pooling = "separate",
                                       selection_xformla,
                                       selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
                                       selection_est_method = "ipw", imputation_pooling = "separate",
                                       pt_assumption, pt_start = NULL, gmm_weighting, gmm_covariance = "tikhonov",
                                       fix_weights, compute_inffunc, pl, cores,
                                       faster_mode, print_details, control_group, base_period, anticipation, alp,
                                       bstrap, cband, biters, dots, entry = c("public", "development")) {
  entry <- match.arg(entry)
  is_custom_est <- is.function(est_method)
  if (!is_custom_est) {
    chained_request_need(!length(dots),
      paste0("Unknown argument(s) ", paste(names(dots), collapse = ", "), ": chained requests forward nothing through `...`; ",
             "ordinary mode (public att_gt(mode = \"ordinary\")) passes `...` to did unchanged."),
      "cdid_chained_request_input_error")
  }
  # ---- mode / panel gate (same meaning as the public wrapper) ----------------
  mode <- check_mode(mode)
  if (is.null(mode)) mode <- if (check_panel(panel)) "chained" else "ordinary"
  if (!identical(mode, "chained")) {
    chained_request_abort(
      "mode = \"ordinary\" is the public att_gt()'s unchanged delegation to did::att_gt(); this development entry point executes chained requests only.",
      "cdid_chained_request_mode_error")
  }
  chained_request_need(isTRUE(check_panel(panel)),
    "mode = \"chained\" requires panel data (panel = TRUE); repeated cross sections are estimated by did::att_gt() in ordinary mode.",
    "cdid_chained_panel")
  if (!is.null(allow_unbalanced_panel)) {
    chained_request_need(chained_request_lgl1(allow_unbalanced_panel), "`allow_unbalanced_panel` must be NULL, TRUE or FALSE.")
    chained_request_need(isTRUE(allow_unbalanced_panel),
      paste0("allow_unbalanced_panel = FALSE is not available in chained mode: chained estimation keeps every eligible unit ",
             "with at least two observed outcomes instead of balancing the panel; use mode = \"ordinary\" for did's balanced estimator."),
      "cdid_chained_balanced")
  }
  # ---- CP02/CP07 scope: explicit requirements and O12-O15 refusals ------------
  if (!is_custom_est && !(chained_request_chr1(est_method) && est_method %in% c("dr", "ipw", "reg"))) {
    shown <- if (chained_request_chr1(est_method)) paste0("\"", est_method, "\"") else "an invalid value"
    chained_request_abort(
      paste0("Chained estimation (cdid ", cdid_version(), ") supports only est_method = \"dr\", est_method = \"ipw\", est_method = \"reg\" or a custom panel IPW function. ",
             "est_method = ", shown, " is not implemented for chained estimation and is refused here without silent replacement. ",
             "Set est_method = \"dr\", \"ipw\" or \"reg\" (or supply a custom function) for chained GMM, or ",
             "mode = \"ordinary\" for did's Callaway-Sant'Anna estimator unchanged."),
      c("cdid_chained_request_est_method", "cdid_chained_unavailable"))
  }
  # ---- observation-correction axis (D52/D54) -----------------------------------
  # "ipw" keeps the existing probability-weighted normalization and never
  # evaluates `imputation_pooling`; "imputation" never evaluates the inactive
  # selection_model/selection_pooling/selection_xformla promises; "dr"
  # (augmented correction) is implemented with est_method = "dr" only and
  # uses both the imputation and the observation-probability arguments.
  sel_method <- chained_request_selection_est_method(selection_est_method)
  imputation <- identical(sel_method, "imputation")
  augmented <- identical(sel_method, "dr")
  if (augmented) chained_request_augmented_est_method(est_method, is_custom_est = is_custom_est)
  sel_norm <- if (imputation) {
    chained_request_imputation_declaration(selection, selection_universe, selection_rho,
      selection_justification, is_custom_est = is_custom_est)
  } else if (augmented) {
    chained_request_augmented_declaration(selection, selection_model, selection_pooling, selection_universe,
      selection_rho, selection_justification, is_custom_est = is_custom_est)
  } else normalize_selection_declaration(
    selection = selection,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    selection_xformla = selection_xformla,
    selection_universe = selection_universe,
    selection_rho = selection_rho,
    selection_justification = selection_justification,
    is_custom_est = is_custom_est
  )
  imp_pooling <- if (imputation || augmented) chained_request_imputation_pooling(imputation_pooling) else NULL
  pt <- validate_chained_pt_request(pt_assumption, pt_start)
  pt_assumption <- pt$pt_assumption
  pt_start <- pt$pt_start
  w_choice <- resolve_gmm_weighting_choice(
    gmm_weighting,
    arg = "gmm_weighting",
    error_fn = chained_request_abort,
    error_class = "cdid_chained_request_input_error"
  )
  gmm_weighting <- w_choice$canonical
  gmm_weighting_requested <- w_choice$requested
  gmm_covariance <- resolve_gmm_covariance_choice(
    gmm_covariance,
    arg = "gmm_covariance",
    error_fn = chained_request_abort,
    error_class = "cdid_chained_request_input_error"
  )
  # D60: NULL, "first_period" and "base_period" are accepted; static calls are
  # identical under all three. "varying" stays refused before data are forced.
  fix_weights <- chained_reference_fix_weights(fix_weights)
  # Inherited did 2.5.1 computation/output options (CP08). compute_inffunc =
  # FALSE returns point estimates only (bstrap/cband become FALSE; the GMM
  # weights still use the complete internal link IFs); pl/cores run the
  # multiplier draws' products on declared worker processes; print_details
  # reports progress. faster_mode keeps its disclosed single chained engine.
  chained_request_need(chained_request_lgl1(compute_inffunc), "`compute_inffunc` must be TRUE or FALSE.")
  chained_request_need(chained_request_lgl1(pl), "`pl` must be TRUE or FALSE.")
  chained_request_need(chained_request_whole(cores, 1) && cores <= .Machine$integer.max,
    "`cores` must be a positive whole number.")
  chained_request_need(chained_request_lgl1(faster_mode) && chained_request_lgl1(print_details),
    "`faster_mode` and `print_details` must be TRUE or FALSE.")
  control_group <- tryCatch(match.arg(control_group, c("nevertreated", "notyettreated")), error = function(e)
    chained_request_abort("`control_group` must be \"nevertreated\" or \"notyettreated\".", "cdid_chained_request_input_error"))
  chained_request_need(chained_request_chr1(base_period) && base_period %in% c("varying", "universal"),
    "`base_period` must be \"varying\" or \"universal\".")
  chained_request_need(chained_request_whole(anticipation, 0), "`anticipation` must be a nonnegative whole number.")
  chained_request_need(is.numeric(alp) && length(alp) == 1L && is.finite(alp) && alp > 0 && alp < 1, "`alp` must lie in (0, 1).")
  chained_request_need(chained_request_lgl1(bstrap) && chained_request_lgl1(cband), "`bstrap` and `cband` must be TRUE or FALSE.")
  chained_request_need(chained_request_whole(biters, 1), "`biters` must be a positive whole number.")
  # All four Boolean flag pairs are bound with their released did 2.5.1 meaning.
  # bstrap = FALSE, cband = TRUE at this att_gt-style endpoint retains the
  # requested DIDparams flags and computes analytical SEs with the normal
  # pointwise interval without any bootstrap draw (att_gt.R:582-590, 688-695);
  # later aggte-style requests inheriting these flags follow their own released
  # endpoint (R/chained_ipw_pipeline.R, section 6).
  list(mode = mode, pt_assumption = pt_assumption, pt_start = pt_start, gmm_weighting = gmm_weighting,
       gmm_weighting_requested = gmm_weighting_requested, gmm_covariance = gmm_covariance,
       control_group = control_group, is_custom_est = is_custom_est, entry = entry,
       selection = sel_norm$requested, selection_canonical = sel_norm$canonical_assumption,
       selection_model = sel_norm$selection_model, selection_pooling = sel_norm$selection_pooling,
       selection_est_method = sel_method, imputation_pooling = imp_pooling,
       inactive_arguments = if (imputation) c("selection_model", "selection_pooling", "selection_xformla") else
         if (augmented) character(0) else "imputation_pooling",
       fix_weights = fix_weights, compute_inffunc = compute_inffunc, pl = pl, cores = as.double(cores),
       print_details = print_details, faster_mode = faster_mode,
       # did 2.5.1 att_gt(): compute_inffunc = FALSE sets bstrap and cband to FALSE.
       bstrap = if (isTRUE(compute_inffunc)) bstrap else FALSE,
       cband = if (isTRUE(compute_inffunc)) cband else FALSE)
}

# Inherited `pl`/`cores` (did 2.5.1 forwards them to the multiplier bootstrap).
# `pl = TRUE` asks for `cores` worker processes, otherwise one. The request
# must fit the declared resource allowance `workers_cap` (public profile
# `control$workers`, default 1 and never raised here, or the development
# budget's `workers`); otherwise it is refused before any data are forced,
# never run serially in silence. Estimation and every point, covariance and
# GMM computation stay serial; only raw-draw products use the workers.
chained_request_parallel_plan <- function(pl, cores, workers_cap, where = "control$workers") {
  requested <- if (isTRUE(pl)) as.double(cores) else 1
  cap <- if (is.null(workers_cap)) 1 else as.double(workers_cap)
  chained_request_need(requested <= cap,
    paste0("pl = TRUE with cores = ", format(cores), " requests ", format(requested), " multiplier bootstrap worker processes, but only ",
      format(cap), " is declared (", where, "). Declare them explicitly (for example `control = list(workers = ", format(requested),
      ")`); the default is not raised and the request is never run serially in silence."),
    c("cdid_chained_request_parallel_error", "cdid_chained_request_input_error"))
  list(pl = pl, cores = as.double(cores), workers = requested, workers_declared = cap,
    route = if (requested > 1) chained_parallel_route() else "serial",
    applies_to = "multiplier bootstrap raw-draw products of att_gt() and later aggte() calls on this result",
    serial = "nuisance fits, catalog/design, covariance streams, GMM, influence functions and analytical inference",
    rng = "Rademacher signs are drawn from the caller's RNG stream by the calling process in serial batch order; results equal cores = 1")
}

# ------------------------------------------------------------------------------
# Public chained executor (cdid::att_gt, mode = "chained", est_method = "dr", "ipw", "reg", or custom fn)
# ------------------------------------------------------------------------------
# One prepared input, fourteen controls produced from its actual dimensions and
# the resolved profile (R/chained_public_control.R), one conventional pipeline
# call without mandatory interval transfer, the D21 full declared pretrend test, and the reviewed cdid_MP.
chained_public_att_gt <- function(yname, tname, idname, gname, xformla, data, panel, allow_unbalanced_panel, control_group,
                                  anticipation, weightsname, fix_weights, alp, bstrap, cband, biters, clustervars, est_method,
                                  base_period, faster_mode, print_details, pl, cores, compute_inffunc, dots, mode, frame,
                                  selection, selection_model = "pair_logit", selection_pooling = "separate",
                                  selection_xformla, selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
                                  selection_est_method = "ipw", imputation_pooling = "separate", gmm_weighting,
                                  gmm_covariance = "tikhonov",
                                  pt_assumption, pt_start = NULL, control, call,
                                  dots_eval = NULL) {
  gate <- chained_request_scope_gate(mode = mode, panel = panel, allow_unbalanced_panel = allow_unbalanced_panel,
    est_method = est_method, selection = selection, selection_model = selection_model, selection_pooling = selection_pooling,
    selection_xformla = selection_xformla,
    selection_universe = selection_universe, selection_rho = selection_rho, selection_justification = selection_justification,
    selection_est_method = selection_est_method, imputation_pooling = imputation_pooling, pt_assumption = pt_assumption, pt_start = pt_start,
    gmm_weighting = gmm_weighting, gmm_covariance = gmm_covariance, fix_weights = fix_weights,
    compute_inffunc = compute_inffunc, pl = pl, cores = cores,
    faster_mode = faster_mode, print_details = print_details, control_group = control_group, base_period = base_period,
    anticipation = anticipation, alp = alp, bstrap = bstrap, cband = cband, biters = biters, dots = dots, entry = "public")
  resolved <- chained_public_resolve_control(
    control,
    gmm_weighting = gate$gmm_weighting_requested,
    gmm_covariance = gate$gmm_covariance
  )
  # Declared bootstrap workers are checked before the data promise is forced.
  plan <- chained_request_parallel_plan(gate$pl, gate$cores, resolved$profile$workers)
  requested_flags <- list(bstrap = bstrap, cband = cband)
  bstrap <- gate$bstrap; cband <- gate$cband
  progress <- chained_progress_emitter(gate$print_details)
  if (!is.null(progress)) progress("chained att_gt(): validating and aligning the input by unit identity")
  custom_dots <- if (isTRUE(gate$is_custom_est) && is.function(dots_eval)) dots_eval() else list()
  est_label <- if (isTRUE(gate$is_custom_est)) "custom" else est_method
  imputation <- identical(gate$selection_est_method, "imputation")
  augmented <- identical(gate$selection_est_method, "dr")
  # ---- validate and align the public input ONCE ------------------------------
  prep <- chained_request_prepare_inputs(data = data, frame = frame, yname = yname, tname = tname, idname = idname,
    gname = gname, xformla = xformla, weightsname = weightsname, clustervars = clustervars,
    anticipation = anticipation, target_cohorts = NULL, selection = gate$selection,
    selection_model = selection_model, selection_pooling = selection_pooling,
    selection_xformla = selection_xformla, selection_universe = selection_universe,
    selection_rho = selection_rho, selection_justification = selection_justification,
    outcome_regression = !imputation && !augmented && identical(est_label, "reg"), or_resource = resolved$profile,
    doubly_robust = !imputation && !augmented && identical(est_label, "dr"),
    selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
    imputation_est_method = if (imputation) est_label else NULL, fix_weights = gate$fix_weights)
  validate_chained_pt_request(gate$pt_assumption, gate$pt_start, prep$calendar)
  # Every pure-imputation treatment method and the augmented route use the
  # automatic row planner.
  if ((imputation || augmented || est_label %in% c("reg", "dr")) && !is.null(prep$cp04))
    prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  # ---- fourteen controls from the prepared input's actual dimensions ----------
  dims <- chained_public_dimensions(prep, biters = biters, clustered = !is.null(prep$cluster_full))
  produced <- chained_public_pipeline_control(
    dims,
    resolved,
    gmm_weighting = gate$gmm_weighting_requested,
    gmm_covariance = gate$gmm_covariance,
    workers = plan$workers
  )
  # ---- one accepted pipeline fit (explicit flags, caller RNG for draws) --------
  public_context <- function(e) chained_request_refusal_context("public",
    prep, settings = list(control_group = gate$control_group, base_period = base_period, gmm_weighting = gate$gmm_weighting,
      gmm_weighting_requested = gate$gmm_weighting_requested, gmm_covariance = gate$gmm_covariance,
      pt_assumption = gate$pt_assumption, pt_start = gate$pt_start,
      anticipation = as.double(anticipation), alp = alp, bstrap = bstrap, cband = cband, biters = as.double(biters),
      xformla = deparse(prep$xformla), weightsname = weightsname, clustervars = prep$clustervar, est_method = est_label,
      selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
      caller_rng = TRUE, run_numerical_transfer = FALSE, compute_inffunc = gate$compute_inffunc,
      multiplier_workers = plan$workers),
    policy = list(profile = resolved$profile, profile_status = resolved$profile_status, source = resolved$source,
      overrides = resolved$overrides, relative_accuracy_source = resolved$relative_accuracy_source,
      numerical_policy = produced$record$numerical_policy, covariance_strategy = produced$record$covariance_strategy,
      chunk_rows = produced$record$chunk_rows,
      numeric_workspace_entries = produced$record$numeric_workspace_entries, counted_work_limit = produced$record$counted_work_limit))
  res <- ipw_pipeline_resignal_with_state(run_chained_ipw_pipeline(
    data = prep$data_internal, idname = idname, tname = tname, yname = yname, gname = gname,
    control = produced$control, xformla = prep$xformla, weightsname = NULL, clustername = NULL,
    cluster_vec = prep$cluster_full, target_cohorts = prep$target_cohorts,
    anticipation = as.integer(anticipation), control_group = gate$control_group,
    base_period = base_period, weighting = gate$gmm_weighting_requested,
    pt_assumption = gate$pt_assumption, pt_start = gate$pt_start,
    run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = character(0L), cband = cband, alp = alp, biters = as.integer(biters),
    btype = "rademacher", seed = NULL, fixed_signs = NULL,
    cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE,
    bstrap = bstrap, caller_rng = TRUE,
    custom_est_method = if (isTRUE(gate$is_custom_est)) est_method else NULL,
    custom_dots = custom_dots,
    gmm_covariance = gate$gmm_covariance, compute_inffunc = gate$compute_inffunc,
    print_details = gate$print_details, prepared = prep$prepared), "cdid_error", "request_context", public_context)
  if (!is.null(progress)) progress("building the did-compatible result",
    if (isTRUE(gate$compute_inffunc)) " and the full declared nonreference pretrend Wald test" else " (point estimates only)")
  prep$prepared <- NULL
  chained_request_bind_mp(res, prep, call = call, params = list(
    yname = yname, tname = tname, idname = idname, gname = gname, xformla = prep$xformla,
    control_group = gate$control_group, anticipation = as.double(anticipation), weightsname = weightsname,
    alp = alp, bstrap = bstrap, cband = cband, biters = biters, clustervars = prep$clustervar,
    est_method = est_method, extra_args = if (isTRUE(gate$is_custom_est)) custom_dots else NULL,
    base_period = base_period, panel = TRUE, selection_est_method = gate$selection_est_method,
    fix_weights = gate$fix_weights, compute_inffunc = gate$compute_inffunc, print_details = gate$print_details,
    pl = gate$pl, cores = gate$cores),
    request = list(mode = gate$mode, panel = panel, allow_unbalanced_panel = allow_unbalanced_panel, faster_mode = faster_mode,
      print_details = print_details, pl = pl, cores = cores, compute_inffunc = compute_inffunc, selection = gate$selection,
      selection_canonical = gate$selection_canonical, selection_model = gate$selection_model, selection_pooling = gate$selection_pooling,
      selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
      inactive_arguments = gate$inactive_arguments,
      est_method = est_label, gmm_weighting = gate$gmm_weighting, gmm_weighting_requested = gate$gmm_weighting_requested,
      gmm_covariance = gate$gmm_covariance, pt_assumption = gate$pt_assumption, pt_start = gate$pt_start, run_numerical_transfer = FALSE,
      fit_aggregations = character(0L), target_cohorts_supplied = FALSE,
      clustervars_supplied = clustervars, caller_rng = TRUE, seed = NULL, btype = "rademacher", entry = "public",
      requested_bstrap = requested_flags$bstrap, requested_cband = requested_flags$cband, execution = plan,
      note = chained_request_inherited_note(gate)),
    control = produced$control, public = list(record = produced$record, dims = dims, resolved = resolved))
}

# Truthful record of the inherited computation/output flags of one request.
chained_request_inherited_note <- function(gate) {
  paste0("faster_mode only selects did's data branch and keeps the single chained engine; print_details = ",
    gate$print_details, if (isTRUE(gate$print_details)) " reported stage and comparison progress through messages" else
      " kept the fit quiet", "; compute_inffunc = ", gate$compute_inffunc,
    if (!isTRUE(gate$compute_inffunc)) " returned point estimates only (bstrap and cband set to FALSE; no influence functions, standard errors, pre-test or aggregation)" else "",
    "; pl = ", gate$pl, ", cores = ", format(gate$cores),
    " (multiplier bootstrap worker processes used when draws run: ", format(if (isTRUE(gate$pl)) gate$cores else 1), ").")
}

# Formal order of the public att_gt() (did 2.5.1 formals, then cdid's extensions)
# is preserved; development inputs follow after the public extensions.
chained_att_gt_request <- function(yname,
                                   tname,
                                   idname = NULL,
                                   gname,
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
                                   ...,
                                   mode = NULL,
                                   frame = NULL,
                                   selection = "none",
                                   selection_model = "pair_logit",
                                   selection_pooling = "separate",
                                   selection_xformla = NULL,
                   selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
                                   selection_est_method = "ipw",
                                   imputation_pooling = "separate",
                                   gmm_weighting = c("optimal", "identity"),
                                   gmm_covariance = "tikhonov",
                                   pt_assumption = c("post", "all"),
                                   pt_start = NULL,
                                   control,
                                   root_witness = NULL,
                                   run_numerical_transfer = TRUE,
                                   target_cohorts = NULL,
                                   fit_aggregations = character(0L),
                                   call = NULL) {
  if (is.null(call)) call <- sys.call()
  gate <- chained_request_scope_gate(mode = mode, panel = panel, allow_unbalanced_panel = allow_unbalanced_panel,
    est_method = est_method, selection = selection, selection_model = selection_model, selection_pooling = selection_pooling,
    selection_xformla = selection_xformla,
    selection_universe = selection_universe, selection_rho = selection_rho, selection_justification = selection_justification,
    selection_est_method = selection_est_method, imputation_pooling = imputation_pooling, pt_assumption = pt_assumption, pt_start = pt_start,
    gmm_weighting = gmm_weighting, gmm_covariance = gmm_covariance, fix_weights = fix_weights,
    compute_inffunc = compute_inffunc, pl = pl, cores = cores,
    faster_mode = faster_mode, print_details = print_details, control_group = control_group, base_period = base_period,
    anticipation = anticipation, alp = alp, bstrap = bstrap, cband = cband, biters = biters, dots = as.list(substitute(list(...)))[-1L], entry = "development")
  mode <- gate$mode; pt_assumption <- gate$pt_assumption; gmm_weighting <- gate$gmm_weighting
  gmm_weighting_requested <- gate$gmm_weighting_requested; gmm_covariance <- gate$gmm_covariance
  control_group <- gate$control_group
  chained_request_need(chained_request_lgl1(run_numerical_transfer), "`run_numerical_transfer` must be TRUE or FALSE.")
  chained_request_need(!missing(control) && inherits(control, "cdid_ipw_pipeline_control"),
    paste0("Explicit development controls are required: supply `control = ipw_pipeline_control(...)` (all fourteen control ",
           "structures). No public default, environment variable or global option supplies them."),
    "cdid_chained_request_control_error")
  chained_request_need(!isTRUE(run_numerical_transfer) || (!missing(root_witness) && is.function(root_witness)),
    "Supply `root_witness(fitted_view)` explicitly (e.g. the accepted automatic one-radius request); no radius is chosen internally.",
    "cdid_chained_request_control_error")
  chained_request_need(is.character(fit_aggregations) && !anyNA(fit_aggregations) &&
    all(fit_aggregations %in% c("simple", "group", "dynamic", "calendar")),
    "`fit_aggregations` must name fit-time aggregation types (normally none; aggregation is requested later).")
  # The explicit development budget declares its multiplier workers; pl/cores
  # select how many this request uses (refused, never ignored, above it).
  plan <- chained_request_parallel_plan(gate$pl, gate$cores, control$multiplier_budget$workers,
    where = "control$multiplier_budget$workers")
  if (is.list(control$multiplier_budget))
    control$multiplier_budget$workers <- if (plan$workers > 1) plan$workers else NULL
  chained_request_need(isTRUE(gate$compute_inffunc) || (!length(fit_aggregations) && !isTRUE(run_numerical_transfer)),
    "compute_inffunc = FALSE returns point estimates only: fit-time aggregations and numerical transfer need influence functions.")
  requested_flags <- list(bstrap = bstrap, cband = cband)
  bstrap <- gate$bstrap; cband <- gate$cband
  custom_dots <- if (isTRUE(gate$is_custom_est)) list(...) else list()
  est_label <- if (isTRUE(gate$is_custom_est)) "custom" else est_method
  imputation <- identical(gate$selection_est_method, "imputation")
  augmented <- identical(gate$selection_est_method, "dr")

  # ---- validate and align the public input once --------------------------------
  prep <- chained_request_prepare_inputs(data = data, frame = frame, yname = yname, tname = tname, idname = idname,
    gname = gname, xformla = xformla, weightsname = weightsname, clustervars = clustervars,
    anticipation = anticipation, target_cohorts = target_cohorts, selection = gate$selection,
    selection_model = selection_model, selection_pooling = selection_pooling,
    selection_xformla = selection_xformla, selection_universe = selection_universe,
    selection_rho = selection_rho, selection_justification = selection_justification,
    outcome_regression = !imputation && !augmented && identical(est_label, "reg"),
    doubly_robust = !imputation && !augmented && identical(est_label, "dr"),
    selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
    imputation_est_method = if (imputation) est_label else NULL, fix_weights = gate$fix_weights)
  validate_chained_pt_request(gate$pt_assumption, gate$pt_start, prep$calendar)

  # ---- one accepted pipeline fit (explicit flags, caller RNG for draws) --------
  development_context <- function(e) chained_request_refusal_context("development",
    prep, settings = list(control_group = control_group, base_period = base_period, gmm_weighting = gmm_weighting,
      gmm_weighting_requested = gmm_weighting_requested, gmm_covariance = gmm_covariance,
      pt_assumption = gate$pt_assumption, pt_start = gate$pt_start,
      anticipation = as.double(anticipation), alp = alp, bstrap = bstrap, cband = cband, biters = as.double(biters),
      xformla = deparse(prep$xformla), weightsname = weightsname, clustervars = prep$clustervar, est_method = est_label,
      selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
      caller_rng = TRUE, run_numerical_transfer = run_numerical_transfer, fit_aggregations = fit_aggregations),
    policy = list(inference_accuracy = control$inference_accuracy,
      gmm_accuracy = control$gmm_accuracy[intersect(c("target_error_unit", "max_target_error_abs", "max_target_error_se", "max_loading_error_abs"), names(control$gmm_accuracy))],
      owner_bounds = control$owner_control[intersect(c("max_retained_entries", "max_matrix_entries", "max_work_ops", "max_chunk_rows"), names(control$owner_control))],
      source = "explicit fourteen-control development object"))
  res <- ipw_pipeline_resignal_with_state(run_chained_ipw_pipeline(
    data = prep$data_internal, idname = idname, tname = tname, yname = yname, gname = gname,
    control = control, xformla = prep$xformla, weightsname = NULL, clustername = NULL,
    cluster_vec = prep$cluster_full, target_cohorts = prep$target_cohorts,
    anticipation = as.integer(anticipation), control_group = control_group,
    base_period = base_period, weighting = gmm_weighting_requested,
    pt_assumption = gate$pt_assumption, pt_start = gate$pt_start,
    run_numerical_transfer = run_numerical_transfer, root_witness = root_witness,
    aggregations = fit_aggregations, cband = cband, alp = alp, biters = as.integer(biters),
    btype = "rademacher", seed = NULL, fixed_signs = NULL,
    cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE,
    bstrap = bstrap, caller_rng = TRUE,
    custom_est_method = if (isTRUE(gate$is_custom_est)) est_method else NULL,
    custom_dots = custom_dots,
    gmm_covariance = if (!isTRUE(run_numerical_transfer)) gmm_covariance else NULL,
    compute_inffunc = gate$compute_inffunc, print_details = gate$print_details,
    prepared = prep$prepared), "cdid_error", "request_context", development_context)

  prep$prepared <- NULL
  chained_request_bind_mp(res, prep, call = call, params = list(
    yname = yname, tname = tname, idname = idname, gname = gname, xformla = prep$xformla,
    control_group = control_group, anticipation = as.double(anticipation), weightsname = weightsname,
    alp = alp, bstrap = bstrap, cband = cband, biters = biters, clustervars = prep$clustervar,
    est_method = est_method, extra_args = if (isTRUE(gate$is_custom_est)) custom_dots else NULL,
    base_period = base_period, panel = TRUE, selection_est_method = gate$selection_est_method,
    fix_weights = gate$fix_weights, compute_inffunc = gate$compute_inffunc, print_details = gate$print_details,
    pl = gate$pl, cores = gate$cores),
    request = list(mode = mode, panel = panel, allow_unbalanced_panel = allow_unbalanced_panel, faster_mode = faster_mode,
      print_details = print_details, pl = pl, cores = cores, compute_inffunc = compute_inffunc, selection = gate$selection,
      selection_canonical = gate$selection_canonical, selection_model = gate$selection_model, selection_pooling = gate$selection_pooling,
      selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
      inactive_arguments = gate$inactive_arguments,
      est_method = est_label, gmm_weighting = gmm_weighting, gmm_weighting_requested = gmm_weighting_requested,
      gmm_covariance = gmm_covariance, pt_assumption = pt_assumption, pt_start = gate$pt_start, run_numerical_transfer = run_numerical_transfer,
      fit_aggregations = fit_aggregations, target_cohorts_supplied = !is.null(target_cohorts),
      clustervars_supplied = clustervars, caller_rng = TRUE, seed = NULL, btype = "rademacher",
      requested_bstrap = requested_flags$bstrap, requested_cband = requested_flags$cband, execution = plan,
      note = chained_request_inherited_note(gate)),
    control = control)
}

# ------------------------------------------------------------------------------
# Input validation and identity alignment (pure R; no owner is built here)
# ------------------------------------------------------------------------------

# Package-created default covariate formula (`xformla = NULL` means "intercept
# only"). A literal `~1` written inside chained_request_prepare_inputs() would
# carry that preparation frame (data, data_internal, cp04, ...) as its
# environment, and the formula travels into DIDparams/MP through `prep$xformla`.
# The formula is built from character input so that `env` actually applies: the
# AST is the same `~1` call (a one-sided formula, intercept 1, no variables) and
# its environment is the base environment, which model.frame()/terms() accept
# for an intercept-only design. Used only where the caller supplied no formula;
# every caller-supplied formula (including a caller's own ~1) keeps its exact
# identity and environment. One fresh formula per call; nothing is cached.
chained_request_default_formula <- function() stats::as.formula("~1", env = baseenv())

# Returns the internal data copy (0 -> Inf recoded, needed columns only, all
# supplied rows retained), the aligned full frame (ids in frame order or first
# appearance order, covariate matrix, weights, cohorts, clusters), the calendar,
# cohort universe, target cohorts and all recoding/exclusion records.
# D60/D61 (`fix_weights` already resolved by the scope gate): the route's
# existing static gates decide staticness and stay unchanged; only a failure
# caused by a genuinely changing covariate or weight becomes a reference input
# (R/chained_reference_inputs.R), and the first-calendar weight V is retained.
chained_request_prepare_inputs <- function(data, frame = NULL, yname, tname, idname, gname, xformla = NULL,
                                           weightsname = NULL, clustervars = NULL, anticipation = 0,
                                           target_cohorts = NULL, selection = "none",
                                           selection_model = "pair_logit", selection_pooling = "separate",
                                           selection_xformla = NULL,
                                           selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
                                           outcome_regression = FALSE, or_resource = NULL,
                                           doubly_robust = FALSE, selection_est_method = "ipw",
                                           imputation_pooling = NULL, imputation_est_method = NULL,
                                           fix_weights = NULL) {
  chained_request_need(is.data.frame(data) && nrow(data) >= 1L, "`data` must be a nonempty data.frame.")
  chained_request_need(chained_request_chr1(idname),
    "Chained (panel) estimation requires `idname`: units are aligned by identity, never by row position.")
  for (nm in c(yname, tname, gname)) chained_request_need(chained_request_chr1(nm), "`yname`, `tname` and `gname` must be single column names.")
  cols <- c(yname, tname, idname, gname)
  chained_request_need(!anyDuplicated(cols), "`yname`, `tname`, `idname` and `gname` must be distinct columns.")
  missing_cols <- cols[!cols %in% names(data)]
  chained_request_need(!length(missing_cols), paste0("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "), "."))
  id <- data[[idname]]; tt <- data[[tname]]; y <- data[[yname]]; g <- data[[gname]]
  chained_request_need((is.numeric(id) || is.character(id)) && !is.factor(id) && !is.object(id) && !anyNA(id) &&
    (!is.numeric(id) || all(is.finite(id))) && (!is.character(id) || all(nzchar(id))),
    "`idname` must be a plain numeric or character identifier without missing values.")
  chained_request_need(is.numeric(tt) && !is.object(tt) && all(is.finite(tt)), "`tname` must be finite numeric dates.")
  chained_request_need(is.numeric(y) && !is.object(y), "`yname` must be numeric (NA marks an unobserved outcome).")
  chained_request_need(is.numeric(g) && !is.object(g) && !anyNA(g) && all(is.finite(g)) && all(g >= 0),
    paste0("`gname` must be finite nonnegative numeric: the treatment date, or the public never-treated code 0 ",
           "(Inf, NA or negative codes are refused; the internal Inf is bound from 0 here, never guessed)."),
    "cdid_chained_request_cohort_error")
  n_rows <- nrow(data)
  first <- which(!duplicated(id))
  unit_id <- id[first]
  n_units <- length(unit_id)
  pos_unit <- match(id, unit_id)
  # One row per (unit, date), decided on the values themselves: the canonical
  # (sorted unit, date) row permutation of the panel owner is formed once here
  # and adjacent rows are compared as IEEE doubles (0 and -0 are equal), never
  # through formatted strings. The sorted unit order and this permutation are
  # the maps the prepared owner construction reuses (CP08 common boundary).
  unit_key <- if (is.numeric(unit_id)) { k <- as.double(unit_id); k[k == 0] <- 0; k } else as.character(unit_id)
  unit_sort <- order(unit_key, method = "radix")
  unit_rank <- integer(n_units); unit_rank[unit_sort] <- seq_len(n_units)
  tt_d <- as.double(tt)
  ui <- unit_rank[pos_unit]
  row_order <- order(ui, tt_d, method = "radix")
  if (n_rows > 1L) {
    uo <- ui[row_order]; to <- tt_d[row_order]
    chained_request_need(!any(uo[-1L] == uo[-n_rows] & to[-1L] == to[-n_rows]),
      "`data` has duplicated (id, period) rows; one row per unit and date is required.", "cdid_chained_request_input_error")
    uo <- to <- NULL
  }
  # cohorts: time-invariant within unit; public 0 -> internal Inf
  g_unit <- g[first]
  chained_request_need(all(g == g_unit[pos_unit]), "`gname` varies within unit; a unit has exactly one treatment date.",
    "cdid_chained_request_cohort_error")
  n_recoded <- sum(g_unit == 0)
  g_unit_int <- g_unit; g_unit_int[g_unit_int == 0] <- Inf
  calendar <- sort(unique(tt_d))
  # Outcome counts define A/B once (reused by every route and the unit table).
  y_observed <- !is.na(y)
  n_obs_unit <- tabulate(pos_unit[y_observed], nbins = n_units)

  # ---- frame: explicit population frame or the panel's units -----------------
  frame_supplied <- !is.null(frame)
  if (frame_supplied) {
    chained_request_need(is.data.frame(frame) && nrow(frame) >= 1L && idname %in% names(frame),
      "`frame` must be a data.frame with the `idname` column (one row per population unit).", "cdid_chained_request_frame_error")
    fid <- frame[[idname]]
    chained_request_need(identical(class(fid), class(id)) && !anyNA(fid) && !anyDuplicated(fid),
      "`frame` ids must have the same type as `data` ids, be unique and non-missing.", "cdid_chained_request_frame_error")
    absent <- unit_id[!unit_id %in% fid]
    chained_request_need(!length(absent),
      paste0(length(absent), " panel unit(s) are missing from `frame` (e.g. ", paste(utils::head(absent, 5L), collapse = ", "),
             "); the frame must cover every unit in `data`."), "cdid_chained_request_frame_error")
    frame_ids <- fid
  } else {
    frame <- data[first, , drop = FALSE]
    rownames(frame) <- NULL
    frame_ids <- unit_id
  }
  n_frame <- length(frame_ids)
  in_panel <- match(frame_ids, unit_id)
  frame_only <- which(is.na(in_panel))
  # Frame position of every data row (first-appearance unit positions when the
  # frame is the panel's own roster); shared by every alignment below.
  row_frame_pos <- if (frame_supplied) match(id, frame_ids) else pos_unit

  # cohorts on the frame: panel units from data; frame-only units from the frame column (0 -> Inf) or unknown
  frame_g <- rep(NA_real_, n_frame)
  frame_g[!is.na(in_panel)] <- g_unit_int[in_panel[!is.na(in_panel)]]
  frame_has_g <- frame_supplied && gname %in% names(frame)
  if (frame_has_g) {
    fg <- frame[[gname]]
    chained_request_need(is.numeric(fg) && !anyNA(fg) && all(is.finite(fg)) && all(fg >= 0),
      "`frame` cohort column must be complete finite nonnegative numeric (0 = never treated); omit the column if cohorts are unknown.",
      "cdid_chained_request_cohort_error")
    fg_int <- as.double(fg); fg_int[fg_int == 0] <- Inf
    chained_request_need(all(fg_int[!is.na(in_panel)] == frame_g[!is.na(in_panel)]),
      "`frame` cohorts disagree with `data` cohorts for some panel units (aligned by id).", "cdid_chained_request_cohort_error")
    frame_g <- fg_int
    n_recoded <- n_recoded + sum(fg[frame_only] == 0)
  }

  # The accepted owner/catalog keeps the complete population frame. A cohort
  # without a pre-anticipation reference stays declared, with unavailable
  # target cells classified by D5; upstream preprocessing must not change the
  # chained sample, cohort shares or IF row universe (D1/D2/D5).
  rownames(frame) <- NULL

  # Pure imputation never evaluates the inactive probability-model promises;
  # the augmented correction evaluates both designs on every F row.
  imputation <- identical(selection_est_method, "imputation")
  augmented <- identical(selection_est_method, "dr")
  sel_norm <- if (imputation) {
    chained_request_imputation_declaration(selection, selection_universe, selection_rho, selection_justification)
  } else if (augmented) {
    chained_request_augmented_declaration(selection, selection_model, selection_pooling, selection_universe,
      selection_rho, selection_justification)
  } else normalize_selection_declaration(
    selection = selection,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    selection_xformla = selection_xformla,
    selection_universe = selection_universe,
    selection_rho = selection_rho,
    selection_justification = selection_justification
  )
  cp04 <- NULL
  imputation_design <- NULL
  or_baseline <- outcome_regression && !sel_norm$enabled && !imputation && !augmented
  # D60: each route's existing static covariate gate decides staticness and is
  # unchanged (no raw-variable test is added to calls it accepts). Only a gate
  # failure caused by a genuinely changing covariate (or the date variable)
  # becomes a dated recipe of that route; any other failure is re-raised. A
  # formula without a dated variable keeps its static design, re-evaluated by
  # its own unchanged gate (so its own failures keep their original errors).
  prov <- NULL
  covariate_gate <- function(fn, formulas, label) tryCatch(fn(), error = function(e) {
    # Formula constructors such as poly() and model.matrix() can fail before
    # the static validator emits its typed error (for example, one factor
    # level on the first-observed frame). Only ordinary formula errors and
    # the existing covariate refusal may enter the genuinely dated branch.
    if (!inherits(e, c("simpleError", "cdid_chained_request_covariate_error"))) stop(e)
    vars <- unique(unlist(lapply(formulas, function(f) if (inherits(f, "formula")) all.vars(f) else character())))
    cls <- chained_reference_classify(vars, data, frame, frame_supplied, id, tt, tname, frame_ids, pos_unit, label)
    if (!any(cls$kind %in% c("dated", "date"))) stop(e)
    chained_reference_check_classes(cls, label)
    structure(list(prov = cls), class = "cdid_reference_dated_marker")
  })
  is_marker <- function(x) inherits(x, "cdid_reference_dated_marker")
  merge_prov <- function(cls) {
    prov <<- if (is.null(prov)) cls else rbind(prov, cls[!cls$variable %in% prov$variable, , drop = FALSE])
  }
  is_dated <- function(f) chained_reference_formula_dated(f, prov)
  # Dated designs are frozen later on their canonical required rows; owners
  # and readers carry identities/outcome support only (no design is invented).
  carrier <- function(nm) matrix(1, n_frame, 1L, dimnames = list(NULL, nm))
  x_static <- NULL; s_static <- NULL; sf <- NULL
  if (imputation) {
    imp <- covariate_gate(function() chained_request_imputation_inputs(xformla, frame, data, frame_ids, id, frame_g, sel_norm,
      imputation_pooling, imputation_est_method, selection_justification, frame_supplied), list(xformla), "xformla")
    if (is_marker(imp)) {
      merge_prov(imp$prov)
      xf <- if (is.null(xformla)) chained_request_default_formula() else xformla
      chained_request_need(!anyNA(frame_g),
        paste0("Pure imputation (selection_est_method = \"imputation\") requires known cohorts on every F row, including ",
          "frame-only rows: the repeat rates a_k use every roster unit of each cohort."))
      frame_Z <- carrier("reference_carrier")
      cp04 <- chained_request_reference_selection_record("imputation", xf, NULL, sel_norm, selection_justification,
        frame_supplied, n_frame, imputation_pooling = imputation_pooling, imputation_est_method = imputation_est_method)
      imputation_design <- cp04$imputation_design
    } else {
      xf <- imp$xformla; frame_Z <- imp$frame_Z; cp04 <- imp$cp04; imputation_design <- imp$design
      x_static <- frame_Z
    }
  } else if (augmented) {
    aug <- covariate_gate(function() chained_request_augmented_inputs(xformla, selection_xformla, frame, data, frame_ids, id, frame_g,
      sel_norm, imputation_pooling, selection_justification, frame_supplied), list(xformla, selection_xformla), "xformla")
    if (is_marker(aug)) {
      merge_prov(aug$prov)
      xf <- if (is.null(xformla)) chained_request_default_formula() else xformla
      sf <- if (is.null(selection_xformla)) xf else selection_xformla
      chained_request_need(!anyNA(frame_g),
        paste0('The augmented observation correction (selection_est_method = "dr") requires known cohorts on every F row, ',
          "including frame-only rows: the repeat rates a_k and the observation-probability fits use every roster unit of each cohort."))
      U_all <- rep(TRUE, n_frame)
      if (!is_dated(xf)) x_static <- cp04_public_design(xf, frame, data, frame_ids, id, U_all, "xformla",
        row_frame_pos = row_frame_pos)
      if (!is_dated(sf)) s_static <- cp04_public_design(sf, frame, data, frame_ids, id, U_all,
        if (is.null(selection_xformla)) "xformla" else "selection_xformla", allow_zero_col = TRUE,
        row_frame_pos = row_frame_pos)
      frame_Z <- carrier("reference_carrier")
      cp04 <- chained_request_reference_selection_record("augmented", xf, sf, sel_norm, selection_justification,
        frame_supplied, n_frame, imputation_pooling = imputation_pooling,
        selection_formula_inherited = is.null(selection_xformla))
      imputation_design <- cp04$augmented_design
    } else {
      xf <- aug$xformla; frame_Z <- aug$frame_Z; cp04 <- aug$cp04; imputation_design <- aug$design
      sf <- cp04$formula
      x_static <- frame_Z; s_static <- cp04$Z
    }
  } else if (!sel_norm$enabled && !outcome_regression && !doubly_robust) {
    # ---- covariates: unit-level, time-invariant, aligned by identity -----------
    xf <- if (is.null(xformla)) chained_request_default_formula() else xformla
    chained_request_need(inherits(xf, "formula") && length(xf) == 2L, "`xformla` must be NULL or a one-sided formula such as ~ x.")
    xvars <- all.vars(xf)
    frame_Z <- covariate_gate(function() {
      if (length(xvars)) {
        src <- if (frame_supplied) "frame" else "data"
        miss_x <- xvars[!xvars %in% names(frame)]
        chained_request_need(!length(miss_x),
          paste0("Covariate(s) ", paste(miss_x, collapse = ", "), " not found in `", src, "`; chained covariates are unit-level ",
                 "attributes supplied once per unit."), "cdid_chained_request_covariate_error")
      }
      fz <- stats::model.matrix(xf, data = frame)
      chained_request_need(nrow(fz) == n_frame && all(is.finite(fz)),
        "Covariates in `xformla` contain missing or nonfinite values for some frame units; NA covariates are refused.",
        "cdid_chained_request_covariate_error")
      if (length(xvars) && all(xvars %in% names(data))) {
        if (cp04_design_direct_columns(stats::terms(xf), xvars, frame, data)) {
          # Direct plain numeric columns: the panel model matrix is the intercept
          # plus these raw columns, so the original row check is exact on the raw
          # values by unit; no NT-row model matrix is built.
          complete <- TRUE; agree <- TRUE
          for (v in xvars) {
            dv <- data[[v]]
            if (anyNA(dv)) { complete <- FALSE; break }
            if (!all(dv == frame[[v]][row_frame_pos])) agree <- FALSE
          }
          chained_request_need(complete, "Covariates in `xformla` have missing values in panel rows.",
            "cdid_chained_request_covariate_error")
          chained_request_need(agree,
            "Covariates in `xformla` vary within unit across time (or disagree with `frame` by id).",
            "cdid_chained_request_covariate_error")
        } else {
          full_Z <- stats::model.matrix(xf, data = data)
          chained_request_need(nrow(full_Z) == nrow(data), "Covariates in `xformla` have missing values in panel rows.",
            "cdid_chained_request_covariate_error")
          expected <- fz[row_frame_pos, , drop = FALSE]
          chained_request_need(max(abs(full_Z - expected)) == 0,
            "Covariates in `xformla` vary within unit across time (or disagree with `frame` by id).",
            "cdid_chained_request_covariate_error")
        }
      } else if (length(xvars) && any(xvars %in% names(data))) {
        chained_request_abort("Covariates are split between `data` and `frame`; supply every covariate of `xformla` in one place.",
          "cdid_chained_request_covariate_error")
      }
      fz
    }, list(xf), "xformla")
    if (is_marker(frame_Z)) {
      merge_prov(frame_Z$prov)
      frame_Z <- carrier("reference_carrier")
    } else x_static <- frame_Z
    rownames(frame_Z) <- NULL

  } else {
    # Outcome counts define A/B; supplied date rows do not define observation windows.
    frame_n_obs <- integer(n_frame)
    frame_n_obs[!is.na(in_panel)] <- n_obs_unit[in_panel[!is.na(in_panel)]]
    A <- frame_n_obs >= 2L
    U <- if (!sel_norm$enabled || selection_universe == "F") rep(TRUE, n_frame) else
      if (selection_universe == "B") frame_n_obs >= 1L else A
    xf <- if (is.null(xformla)) chained_request_default_formula() else xformla
    frame_Z <- covariate_gate(function() cp04_public_design(xf, frame, data, frame_ids, id,
      A, "xformla", row_frame_pos = row_frame_pos), list(xf, if (sel_norm$enabled) selection_xformla), "xformla")
    if (is_marker(frame_Z)) {
      merge_prov(frame_Z$prov)
      frame_Z <- if (is_dated(xf)) carrier("reference_carrier") else cp04_public_design(xf, frame, data, frame_ids, id, A, "xformla",
        row_frame_pos = row_frame_pos)
    }
    if (!is_dated(xf)) x_static <- frame_Z
    if (sel_norm$enabled) {
      chained_request_need(!anyNA(frame_g), "Corollary1 requires known cohorts on every F row, including frame-only rows.")
      sf <- if (is.null(selection_xformla)) xf else selection_xformla
      s_label <- if (is.null(selection_xformla)) "xformla" else "selection_xformla"
      SZ <- if (identical(sf, xf) && identical(U, A) && !is.null(x_static)) frame_Z else
        covariate_gate(function() cp04_public_design(sf, frame, data, frame_ids, id, U, s_label,
          allow_zero_col = TRUE, row_frame_pos = row_frame_pos), list(xf, sf), s_label)
      if (is_marker(SZ)) {
        merge_prov(SZ$prov)
        SZ <- if (is_dated(sf)) carrier("reader_intercept") else
          cp04_public_design(sf, frame, data, frame_ids, id, U, s_label, allow_zero_col = TRUE, row_frame_pos = row_frame_pos)
      }
      if (!is_dated(sf)) s_static <- SZ
      rho <- rep(NA_real_, n_frame)
      if (is.null(selection_rho)) rho[U] <- 1 else {
        chained_request_need(selection_rho %in% names(frame), "`selection_rho` must name a fixed known normalized ratio column on the frame.")
        rho[U] <- frame[[selection_rho]][U]
        chained_request_need(is.numeric(frame[[selection_rho]]) && all(is.finite(rho[U]) & rho[U] > 0),
          "Known rho must be finite and positive on U; it is never estimated or self-normalized.")
      }
      cp04 <- list(Z = SZ, U = U, rho = rho, universe = selection_universe,
        formula = sf, rho_column = selection_rho, justification = selection_justification,
        selection_requested = sel_norm$requested, canonical_assumption = sel_norm$canonical_assumption,
        legacy_alias = sel_norm$legacy_alias, selection_model = sel_norm$selection_model,
        selection_pooling = sel_norm$selection_pooling,
        population = if (frame_supplied) "supplied representative frame" else "supplied panel unit roster (including singletons)")
    } else {
      # A real full-roster reader supplies identities and actual pair rows only.
      # Its intercept is metadata; baseline OR/DR never fits a selection nuisance.
      cp04 <- list(Z = matrix(1, n_frame, 1L, dimnames = list(NULL, "reader_intercept")),
        U = rep(TRUE, n_frame), rho = rep(1, n_frame), universe = "F",
        frame_g = frame_g)
      if (outcome_regression) cp04$or_method <- "baseline"
      if (doubly_robust) cp04$dr_method <- "baseline"
    }
    if (outcome_regression && sel_norm$enabled) cp04$or_method <- "corollary1"
    if (doubly_robust && sel_norm$enabled) cp04$dr_method <- "corollary1"

  }
  dated_x <- !is.null(prov)
  x_dated <- is_dated(xf)

  # ---- weights: positive, time-invariant, aligned by identity ---------------
  # Unchanged static gate; supplied data weights that genuinely change within a
  # unit become a dated source requiring fix_weights = "first_period"/"base_period".
  wres <- chained_reference_weights(data, frame, weightsname, id, tt, calendar, frame_ids, in_panel,
    frame_only, first, pos_unit, frame_supplied)
  first_w <- if (identical(wres$kind, "varying")) wres$V else wres$frame_w
  chained_request_need(!identical(wres$kind, "varying") || !is.null(fix_weights),
    paste0("Sampling weights change within unit (", wres$example, "). With fix_weights = NULL did uses each cell's ",
      "period-specific weight, which chained estimation does not implement; supply fix_weights = \"first_period\" ",
      "(W = w at the first stored calendar date for every context) or \"base_period\" (W = w at each target cohort's ",
      "pre-exposure date b_g, for controls too)."),
    c("cdid_chained_request_unsupported_weights", "cdid_chained_request_weight_error"))
  weight_reduction <- "static"
  if (!dated_x && identical(wres$kind, "static")) {
    frame_w <- wres$frame_w
  } else if (!dated_x && identical(fix_weights, "first_period") && all(is.finite(first_w))) {
    # Unit-attribute X with first_period and w_t1 on the full roster is the
    # existing static system evaluated with w_t1 (reference TeX section 4).
    frame_w <- first_w; weight_reduction <- "first_period_static"
  } else {
    # Every route: mechanical unit case weights for the identity/support owners
    # only; every scientific weight is read at its anchor through the sealed
    # recipe-indexed reader, and required values are demanded on actual supports.
    frame_w <- rep(1.0, n_frame); weight_reduction <- "reference_recipes"
  }
  # ---- clusters: at most one variable beyond idname, nested within units -----
  clustervar <- NULL; cluster_full <- NULL
  if (!is.null(clustervars)) {
    chained_request_need(is.character(clustervars) && !anyNA(clustervars), "`clustervars` must be a character vector of column names.")
    extra <- unique(clustervars[!clustervars %in% c(idname, "")])
    chained_request_need(length(extra) <= 1L,
      paste0("Only one cluster variable beyond `idname` is supported (did 2.5.1 convention); got ", paste(extra, collapse = ", "), "."),
      "cdid_chained_request_cluster_error")
    if (length(extra) == 1L) {
      clustervar <- extra
      in_data <- clustervar %in% names(data); in_frame <- clustervar %in% names(frame)
      chained_request_need(in_data || in_frame, paste0("Cluster variable ", clustervar, " not found in `data` or `frame`."),
        "cdid_chained_request_cluster_error")
      cl_frame <- rep(NA, n_frame)
      if (in_data) {
        cd <- data[[clustervar]]
        chained_request_need(!anyNA(cd) && is.atomic(cd) && !is.factor(cd), "Cluster labels must be plain non-missing values.",
          "cdid_chained_request_cluster_error")
        cl_unit <- cd[first]
        chained_request_need(all(cd == cl_unit[pos_unit]), "Cluster labels vary within unit across panel rows; units must be nested within clusters.",
          "cdid_chained_request_cluster_error")
        cl_frame <- cl_unit[in_panel]
      }
      if (in_frame) {
        cf <- frame[[clustervar]]
        chained_request_need(!anyNA(cf) && is.atomic(cf) && !is.factor(cf), "Frame cluster labels must be plain non-missing values.",
          "cdid_chained_request_cluster_error")
        if (in_data) chained_request_need(all(as.character(cf[!is.na(in_panel)]) == as.character(cl_frame[!is.na(in_panel)])),
          "`frame` cluster labels disagree with `data` cluster labels for some panel units (aligned by id).", "cdid_chained_request_cluster_error")
        cl_frame <- cf
      } else {
        chained_request_need(!length(frame_only),
          "Frame-only units need cluster labels: supply the cluster variable in `frame` (zero-eligible clusters are retained, never guessed).",
          "cdid_chained_request_cluster_error")
      }
      cluster_full <- cl_frame
    }
  }

  # ---- calendar, cohort universe, targets -------------------------------------
  cohorts_known <- sort(unique(frame_g[!is.na(frame_g)]))
  panel_cohorts <- sort(unique(frame_g[!is.na(in_panel)]))
  finite_panel <- panel_cohorts[is.finite(panel_cohorts)]
  with_post <- finite_panel[vapply(finite_panel, function(gg) any(calendar >= gg), logical(1L))]
  without_post <- setdiff(finite_panel, with_post)
  if (is.null(target_cohorts)) {
    tg <- with_post
  } else {
    chained_request_need(is.numeric(target_cohorts) && length(target_cohorts) >= 1L && all(is.finite(target_cohorts)) &&
      all(target_cohorts %in% with_post),
      paste0("`target_cohorts` must be treated cohorts of `data` with at least one calendar date at or after treatment (",
             paste(with_post, collapse = ", "), ")."), "cdid_chained_request_input_error")
    tg <- sort(unique(as.double(target_cohorts)))
  }
  chained_request_need(length(tg) >= 1L,
    "No treated cohort has a calendar date at or after its treatment date; nothing can be estimated.",
    "cdid_chained_request_input_error")

  # D61: on the reference route the retained unit weight is V = w_t1 by original
  # id (NA where the unit has no supplied weight at the first calendar date).
  units <- data.frame(id = frame_ids, cohort = frame_g,
    weight = if (identical(weight_reduction, "reference_recipes")) first_w else frame_w, stringsAsFactors = FALSE)
  units$cluster <- if (is.null(cluster_full)) NA else cluster_full
  units$n_obs <- ifelse(is.na(in_panel), 0L, n_obs_unit[in_panel])
  units$frame_only <- is.na(in_panel)

  data_internal <- data.frame(id, tt, y, frame_g[row_frame_pos], stringsAsFactors = FALSE)
  names(data_internal) <- c(idname, tname, yname, gname)
  rownames(data_internal) <- NULL

  or_basis <- NULL
  if (!imputation && !augmented && outcome_regression && !x_dated) {
    basis <- or_public_basis(frame_Z, units$n_obs >= 2L, or_resource)
    frame_Z <- basis$Z; or_basis <- basis$record
    cp04$outcome_recipe <- basis$outcome_recipe
    cp04$or_basis <- or_basis
    cp04$or_preprocess <- list(retained_entries = unname(cp04_payload(or_basis)["entries"]),
      work_ops = or_basis$preprocessing_work_ops_bound)
  } else if (!imputation && !augmented && doubly_robust && !x_dated) {
    basis <- or_public_basis(frame_Z, units$n_obs >= 2L, or_resource)
    or_basis <- basis$record
    cp04$outcome_recipe <- basis$outcome_recipe
    cp04$or_basis <- or_basis
    cp04$or_preprocess <- list(retained_entries = unname(cp04_payload(or_basis)["entries"]),
      work_ops = or_basis$preprocessing_work_ops_bound)
    attr(frame_Z, "intercept") <- NULL
    attr(frame_Z, "assign") <- NULL
    cp04$treatment_columns <- colnames(frame_Z)
  } else if (!is.null(frame_Z)) {
    attr(frame_Z, "intercept") <- NULL
    attr(frame_Z, "assign") <- NULL
  }

  # D60 source of every reference route: static formulas keep their existing
  # design and alias map; only genuinely dated formulas are frozen later on
  # their canonical required rows. The selection-roster carrier of the reader is
  # metadata only: every scientific design, weight and rate is read through the
  # sealed recipe-indexed reader.
  reference_source <- NULL; reference_route <- NULL
  if (identical(weight_reduction, "reference_recipes")) {
    reference_route <- if (imputation) "imputation" else if (augmented) "augmented" else
      if (sel_norm$enabled) "corollary1" else "baseline"
    route <- chained_reference_route(reference_route, if (sel_norm$enabled) selection_universe else "F")
    route$or_basis <- !imputation && !augmented && (outcome_regression || doubly_robust)
    if (is.null(prov)) {
      # One provenance row per actual variable, including none (e.g. ~1).
      prov_vars <- unique(c(all.vars(xf), if (!is.null(sf)) all.vars(sf)))
      prov <- data.frame(variable = prov_vars, source = rep("existing static gate", length(prov_vars)),
        kind = rep("unit_attribute", length(prov_vars)), example = rep(NA_character_, length(prov_vars)),
        stringsAsFactors = FALSE)
    }
    formulas <- list(x = list(formula = xf, static = x_static, basis = if (x_dated) NULL else or_basis))
    if (!is.null(sf)) formulas$s <- list(formula = sf, static = s_static)
    reference_source <- chained_reference_source(formulas, prov, data, frame, frame_supplied, id, tt, tname, frame_ids,
      calendar, anticipation, wres, fix_weights, route, resource = or_resource)
    if (is.null(cp04)) cp04 <- list(U = rep(TRUE, n_frame), rho = rep(1, n_frame), universe = "F",
      frame_g = frame_g, ipw_method = "baseline")
    cp04$Z <- matrix(1, n_frame, 1L, dimnames = list(NULL, "reader_intercept"))
    cp04$reference_source <- reference_source
    cp04$reference_route <- reference_route
    if (imputation || augmented) cp04$X_frame <- NULL
    if (x_dated) {
      cp04$outcome_recipe <- NULL
      cp04$or_basis <- NULL
      if (!imputation && !augmented) cp04$treatment_columns <- NULL
    }
  }

  # CP08 common boundary: the exact first eligibility of the prepared owner
  # construction (built-in baseline IPW/OR/DR without selection, imputation or
  # reference inputs; complete full roster F = A with plain finite numeric
  # identifiers; one treated cohort plus never-treated controls; unit weights;
  # static intercept/direct plain numeric-column formula; outcomes without
  # NaN/Inf). Everything else keeps the original constructors unchanged. The
  # record is sealed and bound to exactly these objects; no caller data changes.
  prepared <- NULL
  if (is.numeric(id) && n_frame == n_units && !length(frame_only) && length(calendar) >= 2L &&
      n_rows == as.double(n_units) * length(calendar) && all(n_obs_unit == length(calendar)) &&
      !anyNA(frame_g) && !any(is.nan(y) | is.infinite(y)) &&
      identical(weight_reduction, "static") && identical(wres$kind, "static") && all(frame_w == 1) &&
      !imputation && !augmented && !sel_norm$enabled && is.null(reference_source) && !dated_x && !is.null(x_static) &&
      (is.null(cp04) || (is.null(cp04$reference_source) && is.null(cp04$ipw_method) && is.null(cp04$selection_est_method) &&
        identical(cp04$universe, "F") && (identical(cp04$or_method, "baseline") || identical(cp04$dr_method, "baseline")))) &&
      sum(is.finite(panel_cohorts)) == 1L && any(is.infinite(panel_cohorts)) &&
      cp04_design_direct_columns(stats::terms(xf), all.vars(xf), frame, data)) {
    prepared <- chained_request_prepared_record(data_internal, c(idname = idname, tname = tname, yname = yname, gname = gname),
      n_rows, n_units, typeof(id), unit_key[unit_sort], unit_sort, unit_rank, pos_unit, row_order, n_obs_unit,
      frame_ids, n_frame, in_panel, frame_Z, frame_w, frame_g, cp04)
  }

  # The package default formula no longer carries this call frame
  # (chained_request_default_formula). The invocation-only maps are still
  # released here so that they live solely in the returned sealed record until
  # its owner/reader consumers finish, whatever else may reference this frame.
  on.exit({
    prepared <- unit_key <- unit_sort <- unit_rank <- pos_unit <- row_order <- NULL
    n_obs_unit <- row_frame_pos <- ui <- uo <- to <- tt_d <- NULL
  }, add = TRUE)
  list(
    data_internal = data_internal, xformla = xf, cp04 = cp04, or_basis = or_basis, prepared = prepared,
    fix_weights = fix_weights, weight_reduction = weight_reduction, weight_source = wres$kind,
    reference_source = reference_source, reference_route = reference_route, share_weights = first_w,
    share_law = if (identical(weight_reduction, "reference_recipes")) "first_calendar_masses" else NULL,
    imputation_design = imputation_design,
    frame_ids = frame_ids, frame_Z = frame_Z, frame_w = frame_w, frame_g = frame_g,
    # The frame owner refuses NA cohorts: pass the vector only when complete.
    frame_g_arg = if (anyNA(frame_g)) NULL else frame_g,
    cluster_full = cluster_full, clustervar = clustervar, units = units,
    n_frame = n_frame, n_frame_only = length(frame_only), frame_supplied = frame_supplied,
    calendar = calendar, cohort_universe = cohorts_known, panel_cohorts = panel_cohorts,
    target_cohorts = tg, cohorts_without_post_dates = without_post,
    recoding = list(public_never_treated_code = 0, internal_never_treated_code = Inf, n_units_recoded = n_recoded,
      note = "public gname == 0 bound to internal Inf explicitly on a job-local copy; did-facing fields keep original finite cohorts"),
    dropped_units = list(ids = frame_ids[FALSE], n = 0L,
      rule = "none: the complete frame and supplied panel rows are retained; D5 classifies cohorts without an admissible reference"),
    frame_order = if (frame_supplied) "frame row order" else "first appearance in data (any period)"
  )
}

# ------------------------------------------------------------------------------
# CP08 common boundary: invocation-local, source-bound prepared intake record
# ------------------------------------------------------------------------------
# The record carries, for exactly one prepared input, the maps the request layer
# already derived and validated: first-appearance unit positions of every data
# row, the sorted owner order of the units and its inverse, the canonical
# (sorted unit, date) row permutation, per-unit observed-outcome counts and the
# frame alignment (frame row -> sorted unit, sorted A-unit -> frame row). It
# references, never copies, the internal data copy and the frame arrays it was
# built from. The pipeline accepts it only when every one of those objects is
# the very object it receives (identity through `gmm_sexp_aliases`) and the
# record is the sealed snapshot (`ipw_owner_seal`); any other input, mutation or
# tampering silently keeps the original constructors. Nothing persists across
# calls and no caller object is modified.
chained_request_prepared_record <- function(data_internal, cols, n_rows, n_units, id_storage, unit_key_sorted,
                                            unit_sort, unit_rank, pos_unit, row_order, n_obs_unit,
                                            frame_ids, n_frame, in_panel, frame_Z, frame_w, frame_g, cp04) {
  frame_unit_sorted <- unit_rank[in_panel]
  a_pos_in_frame <- integer(n_units)
  a_pos_in_frame[frame_unit_sorted] <- seq_len(n_frame)
  obj <- structure(list(
    source = data_internal, cols = cols, n_rows = n_rows, n_units = n_units, id_storage = id_storage,
    unit_key = unit_key_sorted, unit_sort = unit_sort, unit_rank = unit_rank,
    pos_unit = pos_unit, row_order = row_order, n_obs_unit = n_obs_unit,
    frame_ids = frame_ids, n_frame = n_frame, in_panel = in_panel,
    frame_unit_sorted = frame_unit_sorted, a_pos_in_frame = a_pos_in_frame,
    frame_Z = frame_Z, frame_w = frame_w, frame_g = frame_g,
    cp04 = if (is.null(cp04)) NULL else list(Z = cp04$Z, U = cp04$U, rho = cp04$rho, frame_g = cp04$frame_g),
    scope = "invocation-local prepared intake of one public request; bound by object identity, never by a flag"),
    class = "cdid_request_prepared_intake")
  ipw_owner_seal(obj, "prepared_intake")
}

# Pipeline-side binding (before any owner is constructed): the sealed record
# must describe exactly the objects the pipeline received, on a request the
# prepared constructors are licensed for (no data weights or cluster column
# lookup, no frame-only units, no custom callback, and either no full-roster
# backend or the baseline OR/DR reader on the same Z/U/rho/cohort arrays).
# NULL means "use the original constructors"; it is never an error.
chained_request_prepared_bind <- function(prepared, data, idname, tname, yname, gname, frame_id, frame_Z, frame_weights,
                                          frame_g, allow_frame_only, weightsname, clustername, cp04, custom_est_method) {
  if (is.null(prepared)) return(NULL)
  ok <- inherits(prepared, "cdid_request_prepared_intake") && is.list(prepared) &&
    ipw_owner_seal_ok(prepared, "prepared_intake") &&
    gmm_sexp_aliases(prepared$source, data) &&
    identical(prepared$cols, c(idname = idname, tname = tname, yname = yname, gname = gname)) &&
    is.null(weightsname) && is.null(clustername) && identical(allow_frame_only, FALSE) &&
    !is.function(custom_est_method) &&
    !is.null(frame_id) && !is.null(frame_Z) && !is.null(frame_weights) && !is.null(frame_g) &&
    gmm_sexp_aliases(prepared$frame_ids, frame_id) && gmm_sexp_aliases(prepared$frame_Z, frame_Z) &&
    gmm_sexp_aliases(prepared$frame_w, frame_weights) && gmm_sexp_aliases(prepared$frame_g, frame_g)
  if (!isTRUE(ok)) return(NULL)
  if (is.null(cp04)) {
    if (!is.null(prepared$cp04)) return(NULL)
  } else {
    pc <- prepared$cp04
    if (is.null(pc) || !is.null(cp04$reference_source) || !is.null(cp04$ipw_method) ||
        !is.null(cp04$selection_est_method) || !identical(cp04$universe, "F") ||
        !(identical(cp04$or_method, "baseline") || identical(cp04$dr_method, "baseline")) ||
        !gmm_sexp_aliases(pc$Z, cp04$Z) || !gmm_sexp_aliases(pc$U, cp04$U) ||
        !gmm_sexp_aliases(pc$rho, cp04$rho) || !gmm_sexp_aliases(pc$frame_g, cp04$frame_g)) return(NULL)
  }
  prepared
}

# ------------------------------------------------------------------------------
# Pipeline outputs -> cdid_MP with retained aggregation state
# ------------------------------------------------------------------------------

# Development entry only: the private explicit-control request leaves W/Wpval
# uncomputed. The three policies it names are now settled (D20 level, D21 full
# declared null with the existing certified covariance rule) and bound on the
# PUBLIC path (`chained_public_att_gt` -> `chained_pretest_record`); this record
# keeps its structure so that the development entry remains an honest
# non-computation rather than a second, differently controlled test.
chained_request_pretest_record <- function(target_meta) {
  pre <- which(target_meta$role == "pre")
  avail <- pre[target_meta$available[pre]]
  list(status = "not_computed",
    reason = paste0("the unconditional pretrend Wald test (W/Wpval) is not computed on this private development entry; the public ",
                    "chained att_gt() path computes the full declared null under conventional computation (see policy_entry_point); the direct ",
                    "pretrend comparisons, their joint analytical covariance and complete IFs are retained in this result"),
    unresolved_policies = c(
      null_set = "full declared set of nonreference pretrend comparisons; exact structural reference zeros are identities; bound on the public path only",
      covariance_rule = "full covariance Cholesky/solve and existing conditioning policy; no pseudoinverse, rank reduction or tiny-SE filter; bound on the public path only",
      acceptance_threshold = "complete per-fit numerical-error certification is not performed; bound on the public path only"),
    policy_entry_point = "chained_request_bind_mp(): pretest binding step, after reported-cell inference and before did::MP construction (R/chained_did_request.R); public path: chained_pretest_record() in R/chained_pretest.R",
    candidate_cells = target_meta$target_id[avail],
    unavailable_pretrend_cells = target_meta$target_id[setdiff(pre, avail)],
    n_declared_pretrend = length(pre), n_available_pretrend = length(avail),
    not_conditional_pretest = "this is never the conditional pretrend test, which is the separate cdid::conditional_did_pretest()")
}

# Reporting only: read the worker's already stored H in its declared coordinates.
# A failed diagnostic never changes the returned nuisance's fit status.
chained_request_H_diagnostics <- function(H, n_A) {
  rc <- tryCatch(base::rcond(H, norm = "O"), error = function(e) NA_real_)
  d <- tryCatch(diag(H), error = function(e) numeric())
  rc_ok <- length(rc) == 1L && is.finite(rc)
  d_ok <- length(d) > 0L && all(is.finite(d))
  list(rcond_H_1_estimate = list(status = if (rc_ok) "available" else "unavailable",
      reason = if (rc_ok) NA_character_ else "nonfinite_or_failed_condition_estimate",
      value = if (rc_ok) as.double(rc) else NA_real_),
    method = 'base::rcond(H, norm = "O")', norm = "O", norm_label = "1-norm",
    diagonal_range = list(status = if (d_ok) "available" else "unavailable",
      reason = if (d_ok) NA_character_ else "nonfinite_or_unavailable_H_diagonal",
      min = if (d_ok) min(d) else NA_real_, max = if (d_ok) max(d) else NA_real_),
    n_A = n_A, normalization = "worker weighted logit Hessian divided by n_A",
    interpretation = "estimated reciprocal condition in declared covariate units; no acceptance threshold")
}

# Compact post-plus-auxiliary ledger. Only existing fitted-view reference identity
# is compared; no coefficients, (g,K) key, mask reconstruction or refit is used.
chained_request_nuisance_diagnostics <- function(res, control, custom) {
  if (!is.null(res$nuisance_registry$cp04)) return(cp04_public_diagnostics(res, control))
  aux <- res$auxiliary_reporting
  fits <- aux$fitted_views
  post <- res$point_links_all
  points <- aux$point_links
  D <- length(fits); M <- length(post); Q <- length(points)
  S <- sum(vapply(aux$selections, length, integer(1L)))
  ctrl <- control$owner_control
  # Reserve all headers, fixed strings, map/copy buffers BEFORE building them.
  # Dynamic source strings/coordinates are gated by the existing payload counter.
  entries <- ipw_owner_safe_product(4096, 1.0 + D + M + Q + S,
    min(ctrl$max_retained_entries, ctrl$max_output_entries, ctrl$max_matrix_entries),
    "compact nuisance diagnostics headers/maps")
  work <- ipw_owner_safe_sum(c(8 * entries, 64.0 * (D * D + (M + Q) * D)),
    ctrl$max_work_ops, "compact nuisance diagnostics header and identity work")
  scratch <- 256.0
  live <- max(c(0, res$resource_accounting$simultaneous_entries,
    res$resource_accounting$initial_envelope_entries))
  prior_work <- res$resource_accounting$downstream_work_ops
  if (is.null(prior_work)) prior_work <- 0.0
  gate <- function() {
    ipw_owner_safe_sum(c(live, entries, scratch),
      min(ctrl$max_retained_entries, ctrl$max_matrix_entries),
      "pipeline plus compact nuisance diagnostics counted workspace")
    ipw_owner_safe_sum(c(entries), ctrl$max_output_entries, "compact nuisance diagnostics output")
    ipw_owner_safe_sum(c(prior_work, work), ctrl$max_work_ops,
      "pipeline plus compact nuisance diagnostics counted work")
  }
  gate()
  charge <- function(x) {
    remaining <- ctrl
    remaining$max_work_ops <- ctrl$max_work_ops - prior_work - work
    remaining$max_matrix_entries <- floor(min(ctrl$max_output_entries - entries,
      ctrl$max_matrix_entries - live - entries - scratch,
      ctrl$max_retained_entries - live - entries - scratch) / 4)
    pb <- ipw_point_view_payload_budget(x, remaining)
    entries <<- ipw_owner_safe_sum(c(entries, 4 * pb$entries))
    work <<- ipw_owner_safe_sum(c(work, pb$work))
    gate()
  }
  # These are compact reporting descriptors, statuses and support counts only.
  charge(list(aux$descriptors, aux$selections, aux$descriptor_status, aux$unsupported_counts))
  for (fv in fits) {
    charge(list(fv$nuisance_id, fv$defining_pair, fv$cohorts, fv$dimensions,
      fv$counts, fv$masses, fv$means, fv$diagnostics, fv$glm_control, fv$control,
      fv$point_summary_disclosure))
    p <- fv$dimensions$n_covariates
    scratch <- max(scratch, 8.0 * p * p + 32.0 * p + 256.0)
    work <- ipw_owner_safe_sum(c(work, 32.0 * p * p * p + 64.0 * p * p + 256.0))
    gate() # LU/rcond/diagonal workspace and arithmetic charged before use
  }
  for (pl in post) charge(pl$link_id)
  for (pl in points) charge(pl$link_id)
  # No matrices are copied; the selected scalar/list payload is stripped of
  # foreign attributes/references using the existing diagnostic scalar copier.
  clean <- ipw_pipeline_scalar_state
  records <- list(); distinct <- integer(D)
  for (j in seq_len(D)) {
    fv <- fits[[j]]; hit <- NA_integer_
    if (j > 1L) for (k in seq_len(j - 1L)) {
      if (identical(fv$ref, fits[[k]]$ref)) { hit <- distinct[k]; break }
    }
    if (is.na(hit)) {
      hit <- length(records) + 1L
      records[[hit]] <- list(id = paste0("nuisance:", hit), status = "available", reason = NA_character_,
        source_nuisance_id = fv$nuisance_id, target_cohort = fv$cohorts$g,
        ordered_control_cohorts = clean(fv$cohorts$K), defining_pair = clean(fv$defining_pair),
        cohorts = clean(fv$cohorts), dimensions = clean(fv$dimensions),
        counts = clean(fv$counts), masses = clean(fv$masses), normalizers = clean(fv$means),
        fit = clean(fv$diagnostics), glm_control = clean(fv$glm_control), control = clean(fv$control),
        point_summary_disclosure = clean(fv$point_summary_disclosure),
        conditioning = chained_request_H_diagnostics(fv$factors$H, fv$dimensions$n_a_units),
        roles = list(post = FALSE, auxiliary = FALSE))
    }
    distinct[j] <- hit
  }
  map_link <- function(pl, index, role) {
    hit <- NA_integer_
    if (!is.null(pl$fitted_nuisance_ref)) for (j in seq_len(D)) {
      if (identical(pl$fitted_nuisance_ref, fits[[j]]$ref)) { hit <- distinct[j]; break }
    }
    if (!is.na(hit)) records[[hit]]$roles[[role]] <<- TRUE
    list(index = index, link_id = clean(pl$link_id),
      nuisance_id = if (is.na(hit)) NA_character_ else records[[hit]]$id,
      status = if (is.na(hit)) "unavailable" else "available",
      reason = if (!is.na(hit)) NA_character_ else if (custom) "external_custom_callback_no_builtin_logit_information" else
        if (is.null(pl)) "no_supported_pair_no_fitted_nuisance" else "no_matching_fitted_nuisance_identity")
  }
  post_map <- lapply(seq_len(M), function(j) map_link(post[[j]], j, "post"))
  aux_map <- lapply(seq_len(Q), function(j) map_link(points[[j]], j, "auxiliary"))
  list(status = if (length(records)) "available" else "unavailable",
    reason = if (length(records)) NA_character_ else if (custom) "external_custom_callback_no_builtin_logit_information" else
      "no_supported_pair_no_fitted_nuisance",
    frame_counts = list(n_full_frame_units = res$frame_owner$n_total_units + res$frame_owner$n_frame_only_units,
      n_pair_eligible_A_units = res$frame_owner$n_a_units),
    records = records, post_links = post_map, auxiliary_links = aux_map,
    auxiliary_descriptors = clean(aux$descriptors), auxiliary_selections = clean(aux$selections),
    auxiliary_descriptor_status = clean(aux$descriptor_status), auxiliary_unsupported_counts = lapply(aux$unsupported_counts, clean),
    support_status = res$support_status,
    failed_fit_rows = "not returned after failed fit",
    interpretation = "Kish is descriptive weight concentration, not inference degrees of freedom, effective clusters or population positivity; no numerical certificate or new acceptance gate",
    accounting = list(retained_entries_bound = entries, scratch_entries_bound = scratch,
      work_ops_bound = work, simultaneous_entries_bound = live + entries + scratch,
      scope = "existing counted workspace/work; not process RSS or CPU time"))
}

chained_request_bind_mp <- function(res, prep, call, params, request, control, public = NULL) {
  chained_request_need(inherits(res, "cdid_ipw_pipeline_result"), "Internal: the pipeline did not return its result class.")
  rep_ <- res$reporting
  chained_request_need(is.data.frame(rep_$target_meta) && !is.null(rep_$inference),
    paste0("The pipeline returned a support-only result (", res$support_status, "): ", res$support_reason),
    "cdid_chained_request_support_only")
  a_ids <- res$panel_owner$a_id
  units <- prep$units
  units$eligible <- units$id %in% a_ids
  chained_request_need(identical(sum(units$eligible), as.integer(res$panel_owner$n_a_units)),
    "Internal: eligibility derived from the panel owner does not match the frame.")
  if (is.null(prep$cluster_full)) units$cluster <- NULL
  calendar <- as.double(res$panel_owner$all_times)
  cohort_universe <- as.double(res$panel_owner$cohorts)
  # did-facing glist: declared treated panel cohorts, including unavailable ones.
  glist <- sort(unique(c(prep$panel_cohorts[is.finite(prep$panel_cohorts)], prep$target_cohorts)))
  params$tlist <- calendar
  params$glist <- glist
  gt <- if (is.null(res$numerical_transfer)) NULL else res$numerical_transfer$gmm_transfer
  # D21 binding (public path): full declared nonreference null on the reported
  # cells and their full joint covariance. The development entry records an
  # honest non-computation (structure retained for its accepted tests).
  point_only <- identical(params$compute_inffunc, FALSE)
  pretest <- if (point_only) list(status = "not_computed",
    reason = paste0("compute_inffunc = FALSE (did 2.5.1 convention): point estimates only, so no influence functions, ",
      "covariance or unconditional pretrend Wald test (W/Wpval are NULL)")) else
    if (is.null(public)) chained_request_pretest_record(rep_$target_meta) else
    chained_pretest_record(rep_,
      max_matrix_entries = control$owner_control$max_matrix_entries)
  mp <- chained_did_MP(rep_, prep$frame_ids, call, params, units = units[, intersect(c("id", "cohort", "weight", "cluster", "eligible"), names(units))],
    support = res$support, basis = res$basis,
    band_families = if (is.null(gt) || !length(gt$band_families)) NULL else gt$band_families,
    pretest = if (identical(pretest$status, "computed")) pretest else NULL,
    covariance_handling = res$covariance_handling, share_law = prep$share_law)
  mp$chained$pt <- chained_pt_assumption_record(res, prep, request)
  if (identical(request$pt_assumption, "all")) {
    pretest$interpretation <- mp$chained$pt$pretest_interpretation
    pretest$window_tags <- rep_$target_meta[, intersect(c("target_id", "pt_window_tag", "pt_operator_covered"), names(rep_$target_meta)), drop = FALSE]
  }
  mp$chained$pretest <- pretest
  mp$chained$request <- c(request, list(call = call))
  if (!is.null(params$fix_weights)) mp$chained$request$fix_weights <- params$fix_weights
  if (!is.null(public)) {
    est_entry <- if (is.function(params$est_method)) "custom IPW callback" else paste0("\"", params$est_method, "\"")
    if (identical(params$selection_est_method, "imputation"))
      est_entry <- paste0(est_entry, ", selection_est_method = \"imputation\"")
    if (identical(params$selection_est_method, "dr"))
      est_entry <- paste0(est_entry, ", selection_est_method = \"dr\"")
    mp$chained$public <- c(list(entry = paste0("cdid::att_gt(mode = \"chained\", est_method = ", est_entry, ", control = ...)"), decisions = c("D23", "D21", "D26")),
      public$record)
    mp$chained$numerical_policy <- public$record$numerical_policy
  }
  mp$chained$selection <- res$selection
  # Pure imputation has no common A-only extraction; its record is the raw
  # full-F design. Probability-weighted routes keep their or_basis record.
  mp$chained$outcome_design <- if (is.null(prep$imputation_design)) prep$or_basis else prep$imputation_design
  # D60/D61 record (absent for static calls without an explicit fix_weights).
  reference_inputs <- chained_request_reference_record(prep, res)
  if (!is.null(reference_inputs)) {
    mp$chained$reference_inputs <- reference_inputs
    if (!is.null(reference_inputs$recipes$basis) && is.null(mp$chained$outcome_design))
      mp$chained$outcome_design <- reference_inputs$recipes$basis
  }
  mp$chained$recoding <- prep$recoding
  mp$chained$dropped_units <- prep$dropped_units
  mp$chained$frame_order <- prep$frame_order
  mp$chained$frame$n_obs <- units$n_obs
  mp$chained$frame$frame_only <- units$frame_only
  mp$chained$calendar <- calendar
  mp$chained$cohort_universe <- cohort_universe
  mp$chained$target_cohorts <- prep$target_cohorts
  mp$chained$cohorts_without_post_dates <- prep$cohorts_without_post_dates
  mp$chained$reduction_summary <- res$reduction_summary[intersect(c("weighting", "weighting_requested", "reduce_structural_redundancy", "reduction_active",
    "n_links_all", "n_links_used", "max_moment_expansion_diff"), names(res$reduction_summary))]
  mp$chained$covariance_handling <- res$covariance_handling
  mp$chained$nuisance_diagnostics <- chained_request_nuisance_diagnostics(res, control, is.function(params$est_method))
  mp$chained$resource_accounting <- res$resource_accounting
  mp$chained$numerical_transfer <- if (is.null(gt)) list(requested = FALSE, performed = FALSE,
    status = "complete_per_fit_error_certification_not_requested_or_performed", certified = FALSE) else list(requested = TRUE,
    certified = isTRUE(res$numerical_transfer$certified),
    complete_inference_certified = isTRUE(gt$complete_inference_certified),
    supported_subvector_certified = isTRUE(gt$supported_subvector_certified),
    whole_request_complete = isTRUE(gt$whole_request_complete),
    certified_target_ids = gt$certified_target_ids, unavailable_target_ids = gt$unavailable_target_ids,
    band_family = gt$band_family, band_families = gt$band_families, certificate_scope = gt$certificate_scope,
    total_target_error_cert = gt$total_target_error_cert, total_se_error_cert = gt$total_se_error_cert,
    total_cov_error_cert = gt$total_cov_error_cert, nuisance_root_radii = res$numerical_transfer$nuisance_root_radii,
    reported_error_bounds = rep_$numerical_certificate$errors,
    # Compact Theorem 12 attribution (section 0b of the request-path note): the
    # certified scalars, the seven component vectors and the per-link
    # root/evaluation split, so that a later aggregation refusal can be
    # attributed without the discarded fit. References to vectors the
    # accepted certificate already retains; no IF array, owner or fit.
    error_attribution = ipw_pipeline_scalar_state(ipw_pipeline_transfer_error_attribution(res$numerical_transfer, res$gmm_fit)))
  mp$chained$fit_aggregations <- if (length(res$aggregations)) lapply(res$aggregations, function(z)
    list(type = z$type, status = z$status, estimates = z$estimates, se = z$se, overall = z$overall[c("estimate", "se", "status")])) else NULL
  # State retained for later aggregation requests (no owner, registry, fitted
  # object, GMM fit or duplicate workspace): the reported enclosure boxes, the
  # recoding grid, and the explicit inference controls. Estimates, IFs, frame
  # cohorts/weights/clusters/eligibility already live in the MP fields.
  mp$chained$aggregation_state <- list(
    enclosure = rep_$enclosure, calendar = calendar, cohort_grid = cohort_universe,
    multiplier_budget = control$multiplier_budget, inference_accuracy = control$inference_accuracy,
    btype = "rademacher", caller_rng = TRUE, numerical_transfer_requested = !is.null(gt),
    # Public requests keep the resolved profile and actual dimensions so that a
    # later aggte() with an explicit `biters` rebuilds the multiplier budget
    # from the same profile instead of being refused by the fit-time draw count.
    public = if (is.null(public)) NULL else list(dims = public$dims, resolved = public$resolved,
      chunk_rows = public$record$chunk_rows, entries = public$record$numeric_workspace_entries,
      bytes = public$record$profile$numeric_workspace_bytes),
    # Inherited pl/cores: later aggte() draws use the same admitted workers.
    parallel = request$execution,
    compute_inffunc = !point_only,
    scope = "retained for chained_aggte_request(); aggregation never refits the estimator")
  # D61: first-calendar masses V by original id (NA where absent) on the
  # reference route; static and first-period static calls keep their weights.
  if (!is.null(prep$share_law)) mp$chained$aggregation_state$share_law <- prep$share_law
  mp
}

# ------------------------------------------------------------------------------
# aggte request on a cdid_MP (released formals and inheritance semantics)
# ------------------------------------------------------------------------------

chained_aggte_request <- function(MP,
                                  type = "group",
                                  balance_e = NULL,
                                  min_e = -Inf,
                                  max_e = Inf,
                                  na.rm = FALSE,
                                  bstrap = NULL,
                                  biters = NULL,
                                  cband = NULL,
                                  alp = NULL,
                                  clustervars = NULL,
                                  call = NULL) {
  if (is.null(call)) call <- sys.call()
  chained_request_need(inherits(MP, "cdid_MP") && is.list(MP$chained$aggregation_state),
    "`MP` must be a cdid_MP produced by chained_att_gt_request() (it carries the retained aggregation state).",
    "cdid_chained_request_input_error")
  chained_request_need(chained_request_chr1(type) && type %in% c("simple", "dynamic", "group", "calendar"),
    "`type` must be one of \"simple\", \"dynamic\", \"group\", \"calendar\".")
  chained_request_need(chained_request_lgl1(na.rm), "`na.rm` must be TRUE or FALSE.")
  # did 2.5.1 compute.aggte(): a compute_inffunc = FALSE result cannot be aggregated.
  chained_request_need(!is.null(MP$inffunc) && !isFALSE(MP$chained$aggregation_state$compute_inffunc),
    paste0("This att_gt() result was produced with compute_inffunc = FALSE (point estimates only), so it has no influence ",
      "functions and cannot be aggregated by aggte(). Re-run att_gt() with compute_inffunc = TRUE (the default) to use aggte()."),
    c("cdid_chained_request_no_inffunc_error", "cdid_chained_request_input_error"))
  dp <- MP$DIDparams
  st <- MP$chained$aggregation_state
  resolved <- chained_request_resolve_aggte_options(dp, bstrap = bstrap, biters = biters, cband = cband, alp = alp,
    clustervars = clustervars, idname = dp$idname, mp_clustered = isTRUE(MP$chained$inference$clustered))
  frame <- MP$chained$frame
  cluster_full <- if (isTRUE(MP$chained$inference$clustered)) frame$cluster else NULL
  ids <- MP$chained$target_id
  theta <- stats::setNames(as.double(MP$att), ids)
  # Public results: the multiplier budget is rebuilt from the retained profile
  # and actual dimensions with the resolved draw count (explicit `biters`
  # overrides are honoured under the same finite workspace); development
  # results keep their explicit caller budget unchanged.
  budget <- if (is.null(st$public)) st$multiplier_budget else
    chained_public_multiplier_budget(st$public$dims, resolved$biters, st$public$chunk_rows, st$public$entries, st$public$bytes,
      workers = if (is.null(st$parallel$workers)) 1 else st$parallel$workers)
  # A post-fit refusal (e.g. D20 on an aggregate family) is re-signalled
  # unchanged with the request facts and the MP's retained numerical-transfer
  # summary, including the compact Theorem 12 attribution when the MP carries it.
  aggregation_context <- function(e) list(entry = "aggte", type = type, na_rm = na.rm,
    balance_e = balance_e, min_e = min_e, max_e = max_e, resolved = resolved, caller_rng = isTRUE(st$caller_rng), refit = FALSE,
    fit_numerical_transfer = MP$chained$numerical_transfer,
    fit_numerical_policy = MP$chained$numerical_policy,
    fit_target_ids = MP$chained$target_id, fit_available = MP$chained$available,
    meaning = "post-fit aggregation of retained state; no refit, draw, radius or tolerance change")
  agg <- ipw_pipeline_resignal_with_state(aggregate_chained_ipw_targets(
    theta_hat = theta, target_if_full = MP$inffunc, target_meta = MP$chained$target_meta,
    full_frame_g = frame$cohort, full_frame_w = frame$weight, type = type, anticipation = as.integer(dp$anticipation),
    cluster_full = cluster_full, alp = resolved$alp, cband = resolved$cband, biters = as.integer(resolved$biters),
    btype = st$btype, seed = NULL, multiplier_budget = budget, inference_accuracy = st$inference_accuracy,
    enclosure = st$enclosure, fixed_signs = NULL, eligible_full = frame$eligible, na.rm = na.rm,
    available = MP$chained$available, bstrap = resolved$bstrap, caller_rng = isTRUE(st$caller_rng),
    balance_e = balance_e, min_e = min_e, max_e = max_e, calendar = st$calendar, cohort_grid = st$cohort_grid,
    share_law = st$share_law, frame_ids = frame$id),
    "cdid_ipw_pipeline_accuracy_error", "aggregation_context", aggregation_context)
  obj <- chained_did_AGGTEobj(agg, MP, call, min_e = min_e, max_e = max_e, balance_e = balance_e, requested = resolved)
  if (!is.null(st$public)) {
    obj$chained$numerical_policy <- MP$chained$numerical_policy
    obj$chained$numerical_transfer <- MP$chained$numerical_transfer
  }
  obj$chained$covariance_handling <- MP$chained$covariance_handling
  obj$chained$request <- list(type = type, na_rm = na.rm, balance_e = balance_e, min_e = min_e, max_e = max_e,
    resolved = resolved, caller_rng = isTRUE(st$caller_rng), refit = FALSE, call = call)
  obj
}

# Inherited versus explicit aggte options with upstream meaning (compute.aggte
# 108-123): NULL inherits the MP's requested DIDparams value; an explicit value
# overrides it. Cluster overrides that cannot be honoured are refused instead of
# upstream's warn-and-fall-back. Every Boolean bstrap/cband pair is bound.
chained_request_resolve_aggte_options <- function(dp, bstrap = NULL, biters = NULL, cband = NULL, alp = NULL,
                                                  clustervars = NULL, idname = NULL, mp_clustered = FALSE) {
  pick <- function(x, inherited, name) {
    if (is.null(x)) list(value = inherited, source = "inherited") else list(value = x, source = "explicit")
  }
  b <- pick(bstrap, dp$bstrap, "bstrap"); cb <- pick(cband, dp$cband, "cband")
  bi <- pick(biters, dp$biters, "biters"); al <- pick(alp, dp$alp, "alp")
  chained_request_need(chained_request_lgl1(b$value) && chained_request_lgl1(cb$value), "`bstrap` and `cband` must be NULL, TRUE or FALSE.")
  chained_request_need(chained_request_whole(bi$value, 1), "`biters` must be NULL or a positive whole number.")
  chained_request_need(is.numeric(al$value) && length(al$value) == 1L && is.finite(al$value) && al$value > 0 && al$value < 1,
    "`alp` must be NULL or lie in (0, 1).")
  # bstrap = FALSE, cband = TRUE (inherited or explicit) follows released
  # compute.aggte: group/dynamic/calendar warn and bootstrap their band with
  # analytical reported SEs (l. 365-370, 509-514, 642-647, getSE 811-835);
  # simple has no bin-band family; overall summaries stay analytical.
  inherited_cv <- dp$clustervars
  cv_source <- "inherited"; cv_value <- inherited_cv
  if (!is.null(clustervars)) {
    chained_request_need(is.character(clustervars) && !anyNA(clustervars), "`clustervars` must be NULL or a character vector.")
    extra <- unique(clustervars[!clustervars %in% c(idname, "")])
    cv_source <- "explicit"
    if (!length(extra)) {
      chained_request_need(!mp_clustered,
        paste0("The chained fit clustered on `", inherited_cv, "`; aggte cannot switch to unit-level inference without changing the ",
               "declared cluster operator (refused rather than silently reporting non-clustered SEs)."),
        "cdid_chained_request_cluster_error")
      cv_value <- NULL
    } else {
      chained_request_need(length(extra) == 1L && mp_clustered && identical(extra, inherited_cv),
        paste0("aggte clustervars `", paste(extra, collapse = ", "), "` cannot be honoured: the chained fit ",
               if (mp_clustered) paste0("clustered on `", inherited_cv, "` and cannot switch variables") else
                 "retained no cluster labels (clustervars was not supplied to the estimation request)",
               "; refused rather than upstream's warn-and-fall-back."),
        "cdid_chained_request_cluster_error")
      cv_value <- inherited_cv
    }
  }
  list(bstrap = b$value, cband = cb$value, biters = bi$value, alp = al$value, clustervars = cv_value,
       sources = c(bstrap = b$source, cband = cb$source, biters = bi$source, alp = al$source, clustervars = cv_source))
}
