# Narrow public adapter for D52 pure actual-change imputation
# (`selection_est_method = "imputation"`; IMPUTATION_SPEC.md Sections 1-4, 7.1).
# The accepted private backend (R/chained_imputation.R) is reused unchanged.
# This file only validates the public request, prepares the raw full-F design,
# calls the backend through its sealed X_frame/design-owner contract, and maps
# its compact fields into the existing selection, nuisance-key and diagnostic
# records. Catalog/auxiliary unions, IF replay, GMM, aggregation, clustering and
# bands stay with the existing CP04 consumers.

# ------------------------------------------------------------------------------
# Request validation (before any data promise is forced)
# ------------------------------------------------------------------------------

# Observation-correction axis. "dr" (augmented observation correction) is
# implemented with treatment est_method = "dr" only; the scope gate refuses its
# other treatment combinations before the data are evaluated
# (chained_request_augmented_est_method, R/chained_augmented_public.R).
chained_request_selection_est_method <- function(selection_est_method) {
  chained_request_need(chained_request_chr1(selection_est_method) &&
      selection_est_method %in% c("ipw", "imputation", "dr"),
    paste0('`selection_est_method` must be "ipw" (default probability-weighted observation correction), ',
      '"imputation" (pure actual-change imputation) or "dr" (augmented observation correction with est_method = "dr").'),
    "cdid_chained_request_unsupported_selection")
  selection_est_method
}

# Pure imputation reads only the assumption name and the representative-F
# declaration. selection_model, selection_pooling and selection_xformla are
# inactive and are deliberately not arguments here, so they are never forced.
chained_request_imputation_declaration <- function(selection, selection_universe, selection_rho,
                                                   selection_justification, is_custom_est = FALSE) {
  if (identical(selection, c("none", "conditional_mtar"))) selection <- "none"
  chained_request_need(chained_request_chr1(selection) && selection %in% c("conditional_mtar", "corollary1"),
    paste0('selection_est_method = "imputation" requires selection = "conditional_mtar" (or the legacy assumption ',
      'alias "corollary1"): pure actual-change imputation relies on conditional mean trend-MAR on the representative ',
      'roster F. selection = "none" is refused, not replaced.'),
    "cdid_chained_request_unsupported_selection")
  chained_request_need(!is_custom_est,
    paste0("Pure imputation with a custom IPW callback is unsupported: corrected callbacks have no specified ",
      "imputation forcing or complete upstream score replay."),
    "cdid_chained_request_unsupported_selection")
  chained_request_need(identical(selection_universe, "F"),
    paste0('selection_est_method = "imputation" is defined only on the representative roster selection_universe = "F" ',
      "(full fixed X and known cohorts on every F row); B/A fitting rosters are unavailable for this method."),
    c("cdid_chained_request_unsupported_selection", "cdid_cp04_unsupported_error"))
  # Existing F convention: rho = 1 (NULL transport), optional justification.
  cp04_validate_declaration("F", selection_rho, selection_justification)
  list(enabled = TRUE, requested = selection, canonical_assumption = "conditional_mtar",
    legacy_alias = identical(selection, "corollary1"), selection_model = NULL, selection_pooling = NULL)
}

# Evaluated only when selection_est_method = "imputation".
chained_request_imputation_pooling <- function(imputation_pooling) {
  chained_request_need(chained_request_chr1(imputation_pooling) &&
      imputation_pooling %in% c("separate", "within_cohort", "across_cohorts"),
    '`imputation_pooling` must be "separate" (default), "within_cohort" or "across_cohorts".',
    "cdid_chained_request_unsupported_selection")
  imputation_pooling
}

# ------------------------------------------------------------------------------
# Input preparation: raw outcome/treatment design on EVERY F row
# ------------------------------------------------------------------------------

# Imputations are predicted, and the treatment/counterfactual designs evaluated,
# on the complete roster (frame-only, singleton and zero-observation units
# included). No common A-only alias reduction is applied: every imputation
# pooling group and the shared counterfactual OR use the backend's group-local
# preprocessing, and the treatment logit keeps the raw columns under the
# existing rank/separation refusals. The reader keeps only its intercept
# metadata; the design reaches the backend through `X_frame`.
chained_request_imputation_inputs <- function(xformla, frame, data, frame_ids, id, frame_g, sel_norm,
                                              imputation_pooling, est_method, justification, frame_supplied) {
  n_frame <- length(frame_ids)
  chained_request_need(!anyNA(frame_g),
    paste0("Pure imputation (selection_est_method = \"imputation\") requires known cohorts on every F row, including ",
      "frame-only rows: the repeat rates a_k use every roster unit of each cohort."))
  xf <- if (is.null(xformla)) ~1 else xformla
  Z <- cp04_public_design(xf, frame, data, frame_ids, id, rep(TRUE, n_frame), "xformla")
  intercept <- isTRUE(attr(Z, "intercept"))
  icpt <- if (intercept) which(as.integer(attr(Z, "assign")) == 0L) else integer(0)
  chained_request_need(!intercept || (length(icpt) == 1L && all(Z[, icpt] == 1)),
    "Declared formula intercept column must be a unit column on every F row.", "cdid_or_design_error")
  attr(Z, "intercept") <- NULL
  attr(Z, "assign") <- NULL
  cols <- colnames(Z)
  slopes <- cols[setdiff(seq_along(cols), icpt)]
  recipe <- or_pooled_outcome_recipe(cols, intercept, slopes)
  trt <- if (est_method %in% c("ipw", "dr")) cols else NULL
  design <- list(route = "pure actual-change imputation (selection_est_method = \"imputation\")",
    formula = xf, columns = cols, intercept = intercept,
    intercept_column = if (intercept) cols[icpt] else NULL, slope_columns = slopes,
    outcome_recipe = recipe, treatment_columns = if (is.null(trt)) character(0) else trt,
    n_rows = n_frame,
    row_universe = "every F roster row, including frame-only, singleton and zero-observation units",
    common_A_alias_reduction = FALSE,
    alias_policy = paste("No global A-only reduction: each imputation pooling group and the shared counterfactual OR",
      "apply the accepted group-local QR alias/support preprocessing on their own fitting rows with prediction-row",
      "alias verification; the treatment logit keeps the raw columns under the existing rank/separation refusals."),
    duration_convention = "stored-calendar transition count L_e = i_t - i_s")
  cp04 <- list(Z = matrix(1, n_frame, 1L, dimnames = list(NULL, "reader_intercept")),
    U = rep(TRUE, n_frame), rho = rep(1, n_frame), universe = "F",
    formula = xf, justification = justification,
    selection_requested = sel_norm$requested, canonical_assumption = "conditional_mtar",
    legacy_alias = sel_norm$legacy_alias,
    population = if (frame_supplied) "supplied representative frame" else "supplied panel unit roster (including singletons)",
    selection_est_method = "imputation", imputation_pooling = imputation_pooling,
    imputation_est_method = est_method, outcome_recipe = recipe, treatment_columns = trt,
    X_frame = Z, imputation_design = design)
  list(xformla = xf, frame_Z = Z, cp04 = cp04, design = design)
}

# ------------------------------------------------------------------------------
# Backend call and registry records
# ------------------------------------------------------------------------------

# `reference` (D60): the sealed recipe-indexed reader; a dated outcome formula's
# columns, intercept and recipe come from its frozen canonical-row design.
imputation_public_fit <- function(cp04, reader, views, ctrl, retained, work, automatic_chunk_rows = FALSE,
                                  reference = NULL) {
  external <- max(0, retained - reader$accounting$borrowed_owner_entries)
  recipe <- cp04$outcome_recipe; trt <- cp04$treatment_columns
  if (!is.null(reference) && is.null(recipe)) {
    rr <- chained_reference_outcome_recipe(reference)
    recipe <- rr$recipe
    trt <- if (cp04$imputation_est_method %in% c("ipw", "dr")) rr$columns else NULL
  }
  tryCatch(fit_imputation_backend(reader, views, ctrl$owner_control, ctrl$glm_control, or_public_rank_tol(),
      outcome_recipe = recipe, treatment_columns = trt,
      est_method = cp04$imputation_est_method, imputation_pooling = cp04$imputation_pooling,
      X_frame = if (is.null(reference)) cp04$X_frame else NULL, retained_external_entries = external, external_work_ops = work,
      automatic_chunk_rows = automatic_chunk_rows, reference = reference),
    cdid_or_wls_rank_error = function(e) {
      # The backend message already states that this is not a non-identification claim.
      e$or_columns <- cp04$outcome_recipe$columns
      e$or_stage <- "group-local actual-change imputation or shared counterfactual OR WLS on the full-F design"
      stop(e)
    })
}

imputation_public_nuisance_keys <- function(backend) {
  law <- backend$law
  c(if (length(backend$cells)) paste0("IMP:", law, ":cell:", seq_along(backend$cells)),
    if (length(backend$groups)) paste0("IMP:", law, ":group:", seq_along(backend$groups)),
    if (!is.null(backend$shared_or)) paste0("IMP:", law, ":shared_counterfactual_or:1"),
    if (length(backend$treatments)) paste0("IMP:", law, ":treatment:", seq_along(backend$treatments)))
}

imputation_public_model_text <- function(est_method, intercept) {
  list(imputation_model = paste("actual-change imputation: observed pair changes on S_e = 1 with fixed unit weights w,",
      "local cohort/pair intercepts when xformla has one, slopes on d_e = L_e z(X) shared by imputation_pooling;",
      "predictions on every F row of each needed cohort, standardized with empirical repeat rates a_k"),
    counterfactual_model = if (est_method %in% c("reg", "dr"))
      paste0("shared one-period counterfactual control outcome-change regression on the imputed synthetic control changes (",
        if (isTRUE(intercept)) "local block intercepts; " else "no intercepts; ",
        "slope coefficients shared across blocks on their retained basis)") else
      "none (est_method = \"ipw\" fits no counterfactual OR)",
    treatment_model = if (est_method %in% c("ipw", "dr"))
      "risk-set treatment propensity logit with synthetic weights w*a_k on F" else
      "none (est_method = \"reg\" fits no treatment logit)")
}

imputation_public_selection_record <- function(cp04, backend = NULL, reader = NULL) {
  m <- cp04$imputation_est_method
  pool <- cp04$imputation_pooling
  acc <- if (is.null(backend)) NULL else backend$accounting
  fields <- c("n_cells", "n_groups", "n_imputation_fits", "n_imputation_local_intercepts",
    "n_imputation_slopes_by_group", "n_or_fits", "n_or_blocks", "n_or_local_intercepts", "n_or_shared_slopes",
    "n_treatment_models", "n_treatment_fits", "n_treatment_columns", "n_selection_fits", "n_selection_cells",
    "n_selection_blocks", "n_selection_parameters", "observation_logit_models")
  trt_text <- switch(m,
    ipw = "The synthetic treatment propensity logit (weights w*a_k on F) must be correct.",
    reg = "The shared counterfactual control-mean regression must be correct across all deduplicated blocks.",
    dr = paste("Either the synthetic treatment propensity logit or the shared counterfactual control-mean regression",
      "across all deduplicated blocks must be correct."))
  c(list(method = "imputation", canonical_assumption = "conditional_mtar",
    requested_selection = cp04$selection_requested, legacy_alias = isTRUE(cp04$legacy_alias),
    selection_est_method = "imputation", est_method = m, imputation_pooling = pool,
    observation_model = "none: no observation-probability model is fitted",
    inactive_arguments = c("selection_model", "selection_pooling", "selection_xformla"),
    selection_model = NA_character_, selection_pooling = NA_character_,
    fitting_universe = "F", population = cp04$population, formula = cp04$formula,
    transport = "rho=1", justification = cp04$justification),
    imputation_public_model_text(m, cp04$outcome_recipe$intercept),
    list(duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    empirical_repeat_rates = if (!is.null(backend) && !is.null(backend$reference_ref))
      chained_reference_record(chained_backend_reference(backend))$recipe_repeat_rates else if (is.null(reader)) NULL else
      stats::setNames(as.double(reader$state$a), as.character(reader$state$cohorts)),
    assumptions = paste0("Each needed actual cohort/pair cell requires actual-pair conditional mean trend-MAR on F, ",
      "E[Delta_e | G=k, X, S_e=1] = E[Delta_e | G=k, X], positive cohort/cell masses and a correct linear actual-change ",
      "mean in d_e = L_e z(X) within its ", pool, " imputation pooling group, with group support/rank on the retained basis",
      if (!identical(pool, "separate")) "; shared slopes are an explicit restriction that can be misspecified, e.g. under heterogeneous treated dynamics" else "",
      ". ", trt_text, " PT/no anticipation and the fixed-weight population law must also hold. ",
      "No observation-probability model, date-independence factorization or observed-level conditioning is used."),
    if_scale = "phi_F on ALL N roster rows; core psi=(n_A/N)*phi_F with 1/n_A covariance normalization",
    unsupported = c("selection_universe B/A", "transport columns (rho)", "corrected custom IPW callbacks",
      "augmented selection_est_method = \"dr\"", "observed-level conditioning"),
    models = if (is.null(acc)) NULL else acc[intersect(fields, names(acc))]))
}

# ------------------------------------------------------------------------------
# Public nuisance diagnostics (compact copy, charged before it is built)
# ------------------------------------------------------------------------------

imputation_public_diagnostics <- function(res, control) {
  b <- res$nuisance_registry$cp04
  bound <- 4 * b$accounting$owned_entries + 4096
  cp04_budget(control$owner_control, c(max(c(0, res$resource_accounting$simultaneous_entries,
    res$resource_accounting$initial_envelope_entries)), bound),
    4 * bound, 8 * bound, bound, "Imputation public diagnostic copy")
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st <- reader$state
  icpt <- isTRUE(b$outcome_recipe$intercept)
  num <- function(x, f) vapply(x, function(r) as.double(r[[f]]), numeric(1))
  int <- function(x, f) vapply(x, function(r) as.integer(r[[f]]), integer(1))
  cells <- data.frame(cell_id = int(b$cells, "cell_id"),
    name = vapply(b$cells, `[[`, character(1), "name"), cohort = num(b$cells, "cohort"),
    s = vapply(b$cells, function(r) as.double(r$pair$s), numeric(1)),
    t = vapply(b$cells, function(r) as.double(r$pair$t), numeric(1)),
    duration_steps = int(b$cells, "L"), group_id = int(b$cells, "group_id"),
    n_observed_pairs = int(b$cells, "n_pair"), mass = num(b$cells, "mass"),
    local_intercept = if (icpt) num(b$cells, "alpha") else NA_real_,
    repeat_rate = num(b$cells, "a_k"), stringsAsFactors = FALSE)
  groups <- lapply(b$groups, function(g) list(group_id = g$group_id, name = g$name, pooling = g$pooling,
    cell_ids = g$cell_ids, cohorts = st$cohorts[g$cohorts], rank = g$p, kept_columns = g$kept_columns,
    slopes = g$beta, aliased_columns = g$basis$aliased_columns,
    nontrivial_candidates_retained = g$basis$nontrivial_candidates_retained, rss = g$rss, objective = g$objective))
  names(groups) <- vapply(b$groups, `[[`, character(1), "name")
  so <- b$shared_or
  block_names <- if (is.null(so)) character(0) else vapply(so$blocks, `[[`, character(1), "name")
  counterfactual_or <- if (is.null(so)) NULL else list(rank = so$p, kept_columns = so$kept_columns,
    shared_slopes = so$beta, has_local_intercepts = isTRUE(so$intercept),
    local_block_intercepts = if (isTRUE(so$intercept)) stats::setNames(as.double(so$alpha), block_names) else numeric(0),
    blocks = data.frame(block_id = int(so$blocks, "block_id"), name = block_names,
      s = vapply(so$blocks, function(r) as.double(r$pair$s), numeric(1)),
      t = vapply(so$blocks, function(r) as.double(r$pair$t), numeric(1)),
      duration_steps = int(so$blocks, "L"),
      control_cohorts = vapply(so$blocks, function(r) paste(r$control_cohort_labels, collapse = ","), character(1)),
      synthetic_mass = num(so$blocks, "mass"), stringsAsFactors = FALSE),
    aliased_columns = so$basis$aliased_columns, rss = so$rss, objective = so$objective)
  treatment <- lapply(b$treatments, function(tr) list(trt_id = tr$trt_id, target_cohort = tr$target_cohort,
    control_cohorts = tr$control_cohorts, s = as.double(tr$pair$s), t = as.double(tr$pair$t),
    coefficients = stats::setNames(as.double(tr$coefficients), tr$column_names), n_risk_F = tr$n_risk_F,
    clipping = tr$clipping))
  coefficients <- list(
    imputation_slopes = lapply(groups, `[[`, "slopes"),
    imputation_local_intercepts = if (icpt) stats::setNames(cells$local_intercept, cells$name) else numeric(0),
    has_local_intercepts = icpt,
    # Counterfactual OR fields keep the existing OR/DR diagnostic names.
    shared_slopes = if (is.null(so)) NULL else so$beta,
    local_intercepts = if (is.null(counterfactual_or)) NULL else counterfactual_or$local_block_intercepts,
    treatment = lapply(treatment, `[[`, "coefficients"))
  link_ids <- vapply(b$links, `[[`, character(1), "link_id")
  fitted_catalog <- data.frame(index = seq_along(b$links), link_id = link_ids,
    role = ifelse(startsWith(link_ids, "cp04:auxiliary:"), "auxiliary", "catalog"),
    cohort = num(b$links, "target_cohort"),
    controls = vapply(b$links, function(l) paste(l$control_cohort_labels, collapse = ","), character(1)),
    s = vapply(b$links, function(l) as.double(l$pair$s), numeric(1)),
    t = vapply(b$links, function(l) as.double(l$pair$t), numeric(1)),
    duration_steps = int(b$links, "L"), target_cell = int(b$links, "target_cell_id"),
    or_block = int(b$links, "or_block_id"), treatment_model = int(b$links, "trt_id"),
    mu = num(b$links, "mu"), mu_ipw = num(b$links, "mu_ipw"), mu_reg = num(b$links, "mu_reg"),
    mu_dr = num(b$links, "mu_dr"), d_plus = num(b$links, "d_plus"), d_minus = num(b$links, "d_minus"),
    stringsAsFactors = FALSE)
  support <- list(n_full_frame_units = st$N, n_pair_eligible_A_units = st$n_A,
    n_frame_only_units = res$frame_owner$n_frame_only_units,
    design_source = b$design_source, design_columns = b$design_columns, design_rows = st$N,
    prediction_rows = "every F row of each needed cohort (prediction standardization, not fill-only)",
    observation_probability_models = 0L, observation_probabilities_fitted = FALSE,
    if_rows = "complete original-F rows, including off-A and frame-only units",
    scale = "core psi = (n_A/N) phi_F with covariance divisor n_A; clusters sum original-unit rows")
  diag_copy <- list(cells = cells, groups = groups, counterfactual_or = counterfactual_or,
    treatment = treatment, coefficients = coefficients, fitted_catalog = fitted_catalog, support = support)
  if (unname(cp04_payload(diag_copy)["entries"]) > bound)
    cp04_abort("Imputation diagnostic copy bound underestimated.", "bound")
  c(list(status = "available", method = "imputation", selection_est_method = "imputation",
    est_method = b$est_method, imputation_pooling = b$imputation_pooling),
    imputation_public_model_text(b$est_method, b$outcome_recipe$intercept),
    list(fit_structure = paste("one WLS fit per imputation pooling group (", length(b$groups), " group(s), ",
      length(b$cells), " cell(s))", if (!is.null(so)) ", one joint shared counterfactual OR fit" else "",
      if (length(b$treatments)) paste0(", ", length(b$treatments), " deduplicated treatment logit(s)") else "", sep = ""),
    duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    selection = "conditional_mtar", selection_route = res$selection,
    inactive_selection = b$selection_record,
    outcome_recipe = b$outcome_recipe, treatment_recipe = b$treatment_recipe),
    diag_copy,
    list(accounting = b$accounting,
    frame_counts = list(n_full_frame_units = res$reporting$inference$n_total,
      n_pair_eligible_A_units = res$frame_owner$n_a_units),
    auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
    interpretation = paste("Pure actual-change imputation: the imputation fits use observed pair changes only, and the",
      "counterfactual OR (reg/dr) is a separate model fitted to the imputed synthetic control changes. No observation-probability",
      "model is fitted; selection_model/selection_pooling/selection_xformla are inactive. Complete F-row influence functions",
      paste0("include repeat-rate and imputation paths",
        if (length(b$treatments)) ", treatment propensity paths" else "",
        if (!is.null(so)) ", counterfactual OR paths" else "", ". No numerical certificate."))))
}
