# Narrow public adapter for the D52 augmented observation correction
# (`selection_est_method = "dr"` with treatment `est_method = "dr"`, `"ipw"` or
# `"reg"`; AUGMENTED_SELECTION_SPEC.md sections 1-3 and IMPLEMENTATION_HANDOFF.md
# section 2). The private backend (R/chained_augmented.R) composes the accepted
# observation-probability, actual-change imputation, roster treatment and
# counterfactual OR kernels. This file only validates the public request,
# prepares the two full-F designs (xformla for the imputation, counterfactual OR
# and treatment logits; selection_xformla, or the inherited xformla, for the
# observation probabilities), calls the backend through its sealed X_frame /
# design-owner contract and maps its compact fields into the existing
# selection, nuisance-key and diagnostic records. Catalog/auxiliary unions, IF
# replay, GMM, aggregation, clustering, bands and the pretest stay with the
# existing CP04 consumers.

# ------------------------------------------------------------------------------
# Request validation (before any data promise is forced)
# ------------------------------------------------------------------------------

# Built-in treatment DR, IPW and OR are implemented. Custom callbacks and
# unimplemented est_method strings are refused explicitly; nothing is substituted.
chained_request_augmented_est_method <- function(est_method, is_custom_est = FALSE) {
  chained_request_need(!is_custom_est,
    paste0('selection_est_method = "dr" (augmented observation correction) with a custom IPW callback is ',
      "unsupported: callbacks have no specified augmented forcing or complete upstream score replay. ",
      "Nothing is substituted."),
    "cdid_chained_request_unsupported_selection")
  if (!(chained_request_chr1(est_method) && est_method %in% c("dr", "ipw", "reg"))) {
    chained_request_abort(paste0(
      'selection_est_method = "dr" (augmented observation correction) supports only est_method = "dr", "ipw" or "reg" in cdid ',
      cdid_version(), '; est_method = "', est_method, '" with the augmented correction is not implemented. ',
      "It is refused before the data are evaluated and no other correction or treatment method is substituted."),
      c("cdid_chained_request_unsupported_selection", "cdid_chained_unavailable"))
  }
  invisible(TRUE)
}

# The assumption name, the representative-F declaration and the two active
# probability-model axes. selection_model/selection_pooling keep their meanings
# (and the legacy alias rule) through normalize_selection_declaration().
chained_request_augmented_declaration <- function(selection, selection_model, selection_pooling,
                                                  selection_universe, selection_rho,
                                                  selection_justification, is_custom_est = FALSE) {
  if (identical(selection, c("none", "conditional_mtar"))) selection <- "none"
  chained_request_need(chained_request_chr1(selection) && selection %in% c("conditional_mtar", "corollary1"),
    paste0('selection_est_method = "dr" requires selection = "conditional_mtar" (or the legacy assumption alias ',
      '"corollary1"): the augmented observation correction relies on conditional mean trend-MAR on the ',
      'representative roster F. selection = "none" is refused, not replaced.'),
    "cdid_chained_request_unsupported_selection")
  chained_request_need(identical(selection_universe, "F"),
    paste0('selection_est_method = "dr" is defined only on the representative roster selection_universe = "F" ',
      "(fixed X and known cohorts on every F row); B/A fitting rosters are unavailable for this method."),
    c("cdid_chained_request_unsupported_selection", "cdid_cp04_unsupported_error"))
  normalize_selection_declaration(
    selection = selection,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    selection_universe = "F",
    selection_rho = selection_rho,
    selection_justification = selection_justification,
    is_custom_est = is_custom_est
  )
}

# ------------------------------------------------------------------------------
# Input preparation: both designs on EVERY F row
# ------------------------------------------------------------------------------

# xformla defines the imputation j, the counterfactual OR (reg/dr) and the
# treatment logit (ipw/dr); it is evaluated on every F row exactly as for pure
# imputation (no common A-only alias reduction). selection_xformla (NULL
# inherits xformla) defines the observation logits and is the reader design.
# The selection design never doubles as the outcome design: the latter reaches
# the backend only through X_frame.
chained_request_augmented_inputs <- function(xformla, selection_xformla, frame, data, frame_ids, id, frame_g,
                                             sel_norm, imputation_pooling, justification, frame_supplied,
                                             est_method = "dr") {
  n_frame <- length(frame_ids)
  chained_request_need(!anyNA(frame_g),
    paste0('The augmented observation correction (selection_est_method = "dr") requires known cohorts on every F row, ',
      "including frame-only rows: the repeat rates a_k and the observation-probability fits use every roster unit of each cohort."))
  U_all <- rep(TRUE, n_frame)
  xf <- if (is.null(xformla)) ~1 else xformla
  Z <- cp04_public_design(xf, frame, data, frame_ids, id, U_all, "xformla")
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
  inherited <- is.null(selection_xformla)
  sf <- if (inherited) xf else selection_xformla
  SZ <- cp04_public_design(sf, frame, data, frame_ids, id, U_all,
    if (inherited) "xformla" else "selection_xformla", allow_zero_col = TRUE)
  design <- list(route = paste0('augmented observation correction (selection_est_method = "dr", est_method = "', est_method, '")'),
    formula = xf, columns = cols, intercept = intercept,
    intercept_column = if (intercept) cols[icpt] else NULL, slope_columns = slopes,
    outcome_recipe = recipe, treatment_columns = if (is.null(trt)) character(0) else trt,
    selection_formula = sf, selection_formula_inherited = inherited,
    selection_columns = colnames(SZ), selection_intercept = isTRUE(attr(SZ, "intercept")),
    n_rows = n_frame,
    row_universe = "every F roster row, including frame-only, singleton and zero-observation units",
    common_A_alias_reduction = FALSE,
    alias_policy = paste("No global A-only reduction: each imputation pooling group and the shared counterfactual OR",
      "apply the accepted group-local QR alias/support preprocessing on their own fitting rows with prediction-row",
      "alias verification; the treatment and observation-probability logits keep their raw columns under the existing",
      "rank/separation/boundary refusals."),
    duration_convention = "stored-calendar transition count L_e = i_t - i_s")
  cp04 <- list(Z = SZ, U = U_all, rho = rep(1, n_frame), universe = "F",
    formula = sf, justification = justification,
    selection_requested = sel_norm$requested, canonical_assumption = "conditional_mtar",
    legacy_alias = sel_norm$legacy_alias, selection_model = sel_norm$selection_model,
    selection_pooling = sel_norm$selection_pooling,
    population = if (frame_supplied) "supplied representative frame" else "supplied panel unit roster (including singletons)",
    selection_est_method = "dr", augmented_est_method = est_method, imputation_pooling = imputation_pooling,
    outcome_formula = xf, selection_formula_inherited = inherited,
    outcome_recipe = recipe, treatment_columns = trt, X_frame = Z, augmented_design = design)
  list(xformla = xf, frame_Z = Z, cp04 = cp04, design = design)
}

# ------------------------------------------------------------------------------
# Backend call and registry records
# ------------------------------------------------------------------------------

# `reference` (D60): the sealed recipe-indexed reader; a dated outcome formula's
# columns, intercept and recipe come from its frozen canonical-row design.
augmented_public_fit <- function(cp04, reader, views, ctrl, retained, work, automatic_chunk_rows = FALSE,
                                 reference = NULL) {
  external <- max(0, retained - reader$accounting$borrowed_owner_entries)
  est_m <- if (is.null(cp04$augmented_est_method)) "dr" else cp04$augmented_est_method
  recipe <- cp04$outcome_recipe; trt <- if (est_m %in% c("ipw", "dr")) cp04$treatment_columns else NULL
  if (!is.null(reference) && is.null(recipe)) {
    rr <- chained_reference_outcome_recipe(reference)
    recipe <- rr$recipe; trt <- if (est_m %in% c("ipw", "dr")) rr$columns else NULL
  }
  tryCatch(fit_augmented_backend(reader, views, ctrl$owner_control, ctrl$glm_control, or_public_rank_tol(),
      outcome_recipe = recipe, treatment_columns = trt,
      est_method = est_m, imputation_pooling = cp04$imputation_pooling,
      selection_model = cp04$selection_model, selection_pooling = cp04$selection_pooling,
      X_frame = if (is.null(reference)) cp04$X_frame else NULL, retained_external_entries = external,
      external_work_ops = work, automatic_chunk_rows = automatic_chunk_rows, reference = reference,
      outcome_pooling = if (is.null(cp04$outcome_pooling)) "pooled" else cp04$outcome_pooling),
    cdid_or_wls_rank_error = function(e) {
      # The backend message already states that this is not a non-identification claim.
      e$or_columns <- recipe$columns
      e$or_stage <- paste("group-local actual-change imputation or shared counterfactual OR WLS on the full-F design",
        "(augmented observation correction)")
      stop(e)
    })
}

augmented_public_nuisance_keys <- function(backend) {
  law <- backend$law
  sb <- backend$selection_backend
  c(if (length(backend$cells)) paste0("AUG:", law, ":cell:", seq_along(backend$cells)),
    if (length(backend$groups)) paste0("AUG:", law, ":group:", seq_along(backend$groups)),
    if (length(sb$blocks)) paste0("AUG:", law, ":selection_block:", seq_along(sb$blocks)),
    if (!is.null(backend$shared_or)) paste0("AUG:", law, if (reconstruction_separate(backend$shared_or)) ":counterfactual_block:" else ":shared_counterfactual_or:",
      if (reconstruction_separate(backend$shared_or)) seq_len(backend$shared_or$n_blocks) else 1L),
    if (length(backend$treatments)) paste0("AUG:", law, ":treatment:", seq_along(backend$treatments)))
}

augmented_public_model_text <- function(intercept, selection_model, selection_pooling, imputation_pooling,
                                        outcome_pooling = "pooled", est_method = "dr") {
  list(
    imputation_model = paste0("actual-change imputation j: observed pair changes on S_e = 1 with fixed unit weights w, ",
      "local cohort/pair intercepts when xformla has one, slopes on d_e = L_e z(X) shared by imputation_pooling = \"",
      imputation_pooling, "\"; predictions on every F row of each needed cohort"),
    observation_model = paste0("selection_model = \"", selection_model, "\" observation-probability logit q fitted on every F row ",
      "of each cohort with selection_pooling = \"", selection_pooling, "\" from selection_xformla (NULL inherits xformla)",
      if (identical(selection_model, "independent_period_logit"))
        "; the pair probability is the product of the two endpoint date probabilities" else
        "; one direct logit for observing both endpoints of each actual pair"),
    augmented_response = paste("Z = j + S_e (Delta_e - j) / max(q, 1e-6) on observed pairs and Z = j elsewhere;",
      "Delta and the floored inverse are evaluated only where S_e = 1; the observed residual term is not renormalized"),
    counterfactual_model = if (est_method %in% c("reg", "dr"))
      paste0(if (identical(outcome_pooling, "separate")) "comparison-specific one-period counterfactual control outcome-change regressions on the reconstructed " else "shared one-period counterfactual control outcome-change regression on the reconstructed ",
        "synthetic control changes Z (", if (isTRUE(intercept)) "local block intercepts; " else "no intercepts; ",
        if (identical(outcome_pooling, "separate")) "unrestricted block slopes on their local retained bases)" else "slope coefficients shared across every deduplicated block on their retained basis)") else
      "none (est_method = \"ipw\" fits no counterfactual OR)",
    treatment_model = if (est_method %in% c("ipw", "dr"))
      "risk-set treatment propensity logit with roster synthetic weights w*a_k on F (no observation-probability weight)" else
      "none (est_method = \"reg\" fits no treatment logit)")
}

augmented_public_selection_record <- function(cp04, backend = NULL, reader = NULL) {
  m <- if (!is.null(backend) && !is.null(backend$est_method)) backend$est_method else
    if (is.null(cp04$augmented_est_method)) "dr" else cp04$augmented_est_method
  acc <- if (is.null(backend)) NULL else backend$accounting
  fields <- c("n_cells", "n_groups", "n_imputation_fits", "n_imputation_local_intercepts",
    "n_imputation_slopes_by_group", "n_or_fits", "n_or_blocks", "n_or_local_intercepts", "n_or_shared_slopes",
    "n_treatment_models", "n_treatment_fits", "n_treatment_columns", "n_selection_fits", "n_selection_cells",
    "n_selection_blocks", "n_selection_parameters", "n_selection_pair_cohorts", "observation_logit_models")
  sb_sum <- selection_backend_summary(backend)
  date <- identical(cp04$selection_model, "independent_period_logit")
  pool <- cp04$imputation_pooling
  icpt <- if (!is.null(backend) && !is.null(backend$outcome_recipe)) backend$outcome_recipe$intercept else cp04$outcome_recipe$intercept
  trt_assump <- switch(m,
    ipw = "In addition the synthetic treatment propensity logit (weights w*a_k on F) must be correct, ",
    reg = paste0("In addition the counterfactual control-mean regressions under the chosen outcome pooling (",
      if (identical(cp04$outcome_pooling, "separate")) "comparison-specific" else "shared across deduplicated blocks",
      ", each used cell needing one reconstruction arm) must be correct, "),
    dr = paste0("In addition EITHER the synthetic treatment propensity logit OR the ",
      "counterfactual control-mean regressions under the chosen outcome pooling (each used cell needs one reconstruction arm) must be correct, "))
  inf_paths <- switch(m,
    ipw = "probabilities, imputation, repeat rates, treatment logits, both normalizers); core",
    reg = "probabilities, imputation, repeat rates, counterfactual OR, treated normalizer); core",
    dr = "probabilities, imputation, repeat rates, treatment logits, shared counterfactual OR, both normalizers); core")
  inf_clip <- switch(m,
    reg = "Clipping can change the target and the observation robustness branches; ordinary IF inference",
    "Clipping can change the target and the observation/treatment robustness branches; ordinary IF inference")
  c(list(method = "augmented", canonical_assumption = "conditional_mtar",
    requested_selection = cp04$selection_requested, legacy_alias = isTRUE(cp04$legacy_alias),
    selection_est_method = "dr", est_method = m, imputation_pooling = pool, outcome_pooling = cp04$outcome_pooling,
    selection_model = cp04$selection_model, selection_pooling = cp04$selection_pooling,
    inactive_arguments = if (m %in% c("reg", "dr")) character(0) else "outcome_pooling",
    fitting_universe = "F", population = cp04$population, formula = cp04$formula,
    outcome_formula = cp04$outcome_formula,
    selection_formula_inherited = isTRUE(cp04$selection_formula_inherited),
    transport = "rho=1", justification = cp04$justification),
    augmented_public_model_text(icpt, cp04$selection_model, cp04$selection_pooling, pool, cp04$outcome_pooling, m),
    list(duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    active_cells = if (is.null(sb_sum)) data.frame(cell_id = integer(), kind = character(), cohort = numeric(),
      label = character(), block_id = integer(), stringsAsFactors = FALSE) else sb_sum$active_cells,
    n_active_cells = if (is.null(sb_sum)) 0L else sb_sum$n_active_cells,
    n_selection_blocks = if (is.null(sb_sum)) 0L else sb_sum$n_selection_blocks,
    n_selection_parameters = if (is.null(sb_sum)) 0L else sb_sum$n_selection_parameters,
    empirical_repeat_rates = if (!is.null(backend) && !is.null(backend$reference_ref))
      chained_reference_record(chained_backend_reference(backend))$recipe_repeat_rates else if (is.null(reader)) NULL else
      stats::setNames(as.double(reader$state$a), as.character(reader$state$cohorts)),
    assumptions = paste0(
      "Each needed actual cohort/pair cell requires conditional mean trend-MAR given the common conditioning set ",
      "W = (xformla covariates, variables appearing only in selection_xformla): E[Delta_e | G=k, W, S_e=1] = E[Delta_e | G=k, W], ",
      "with positive observation probability. Each needed cell then needs EITHER the imputation arm (a correct linear actual-change ",
      "mean given W, represented by the xformla design, for every cell of its ", pool, " imputation pooling group) ",
      "OR the observation arm (a correct ", cp04$selection_model, " model given W for every cell of its ",
      cp04$selection_pooling, " selection block",
      if (date) ", correct date logits at both endpoint blocks and used-pair conditional date independence R_s _||_ R_t | G, W on F" else "",
      "); the arms may differ across cells. ", trt_assump,
      "with PT/no anticipation and the fixed-weight population law. Original-target branch claims also require ",
      "clipping to be inactive in population or the clipped moments to preserve that target. ",
      "These conditions are sufficient, not necessary. ",
      "Correctness of an xformla-only imputer alone does not establish reconstruction when selection-only variables shift ",
      "the outcome-change mean; it requires the stated compatibility (for example q given W correct, or its error depending ",
      "on xformla covariates only). No observed-level conditioning is used."),
    inference = paste("Complete original-unit sandwich on all F rows with every fitted-model path retained (observation",
      inf_paths,
      "psi = (n_A/N) phi_F with 1/n_A covariance normalization; clusters sum original-unit rows. Recovered joint ATT/GMM",
      "coordinates additionally need compatible catalog moments preserving the target and the stated identification",
      "premises.", inf_clip,
      "also requires a negligible threshold-crossing remainder."),
    if_scale = "phi_F on ALL N roster rows; core psi=(n_A/N)*phi_F with 1/n_A covariance normalization",
    unsupported = c("selection_universe B/A", "transport columns (rho)", "corrected custom IPW callbacks",
      "observed-level conditioning", "B/A augmentation",
      "SMAR / Corollary 2", "singleton outcome levels", "cross-fitting / ML nuisances"),
    models = if (is.null(acc)) NULL else acc[intersect(fields, names(acc))]))
}

# ------------------------------------------------------------------------------
# Public nuisance diagnostics (compact copy, charged before it is built)
# ------------------------------------------------------------------------------

augmented_public_diagnostics <- function(res, control, public_object = NULL) {
  b <- res$nuisance_registry$cp04
  bound <- if (cp04_field_diagnostics(b))
    cp04_public_diagnostic_admit(res, control, public_object, "Augmented public diagnostic copy") else {
    bound <- 4 * b$accounting$owned_entries + 4096
    cp04_budget(control$owner_control, c(max(c(0, res$resource_accounting$simultaneous_entries,
      res$resource_accounting$initial_envelope_entries)), bound),
      4 * bound, 8 * bound, bound, "Augmented public diagnostic copy")
    bound
  }
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st <- reader$state
  sb <- b$selection_backend
  icpt <- isTRUE(b$outcome_recipe$intercept)
  num <- function(x, f) vapply(x, function(r) as.double(r[[f]]), numeric(1))
  int <- function(x, f) vapply(x, function(r) as.integer(r[[f]]), integer(1))
  cells <- data.frame(cell_id = int(b$cells, "cell_id"),
    name = vapply(b$cells, `[[`, character(1), "name"), cohort = num(b$cells, "cohort"),
    s = vapply(b$cells, function(r) as.double(r$pair$s), numeric(1)),
    t = vapply(b$cells, function(r) as.double(r$pair$t), numeric(1)),
    duration_steps = int(b$cells, "L"), group_id = int(b$cells, "group_id"),
    selection_pair_cohort = int(b$cells, "q_id"),
    n_observed_pairs = int(b$cells, "n_pair"), mass = num(b$cells, "mass"),
    local_intercept = if (icpt) num(b$cells, "alpha") else NA_real_,
    repeat_rate = num(b$cells, "a_k"),
    roster_mass = num(b$cells, "M_one"),
    probability_weighted_mass = num(b$cells, "MW_one"),
    mass_discrepancy = num(b$cells, "MD_one"),
    imputation_total = num(b$cells, "Y_one"),
    augmented_total = num(b$cells, "Yaug_one"),
    stringsAsFactors = FALSE)
  groups <- lapply(b$groups, function(g) list(group_id = g$group_id, name = g$name, pooling = g$pooling,
    cell_ids = g$cell_ids, cohorts = st$cohorts[g$cohorts], rank = g$p, kept_columns = g$kept_columns,
    slopes = g$beta, aliased_columns = g$basis$aliased_columns,
    nontrivial_candidates_retained = g$basis$nontrivial_candidates_retained, rss = g$rss, objective = g$objective))
  names(groups) <- vapply(b$groups, `[[`, character(1), "name")
  sb_sum <- selection_backend_summary(b)
  selection <- list(selection_model = sb$selection_model, selection_pooling = sb$selection_pooling,
    formula_columns = b$selection_columns, has_intercept = isTRUE(sb$design_spec$has_intercept),
    n_slope_columns = sb$design_spec$d_s, design_branch = sb$design_spec$branch,
    n_blocks = length(sb$blocks), n_cells = length(sb$cells), n_parameters = length(sb$coefficients),
    active_cells = sb_sum$active_cells, blocks = selection_blocks_diagnostic_list(sb),
    coefficients = sb$coefficients)
  so <- b$shared_or
  block_names <- if (is.null(so)) character(0) else vapply(so$blocks, `[[`, character(1), "name")
  counterfactual_or <- if (is.null(so)) NULL else list(rank = so$p, kept_columns = so$kept_columns,
    shared_slopes = so$beta, has_local_intercepts = isTRUE(so$intercept),
    local_block_intercepts = if (isTRUE(so$intercept)) stats::setNames(as.double(so$alpha), block_names) else numeric(0),
    response = "reconstructed augmented response Z",
    blocks = data.frame(block_id = int(so$blocks, "block_id"), name = block_names,
      s = vapply(so$blocks, function(r) as.double(r$pair$s), numeric(1)),
      t = vapply(so$blocks, function(r) as.double(r$pair$t), numeric(1)),
      duration_steps = int(so$blocks, "L"),
      control_cohorts = vapply(so$blocks, function(r) paste(r$control_cohort_labels, collapse = ","), character(1)),
      synthetic_mass = num(so$blocks, "mass"), stringsAsFactors = FALSE),
    aliased_columns = so$basis$aliased_columns, rss = so$rss, objective = so$objective)
  if (reconstruction_separate(so)) counterfactual_or <- list(outcome_pooling = "separate",
    response = "reconstructed control outcome changes", blocks = outcome_public_blocks(b),
    n_fits = so$n_blocks, has_local_intercepts = isTRUE(b$outcome_recipe$intercept),
    local_block_intercepts = if (isTRUE(b$outcome_recipe$intercept)) stats::setNames(
      vapply(so$blocks, `[[`, numeric(1), "alpha"), block_names) else numeric(),
    rss = so$rss, objective = so$objective)
  treatment <- lapply(b$treatments, function(tr) list(trt_id = tr$trt_id, target_cohort = tr$target_cohort,
    control_cohorts = tr$control_cohorts, s = as.double(tr$pair$s), t = as.double(tr$pair$t),
    coefficients = stats::setNames(as.double(tr$coefficients), tr$column_names), n_risk_F = tr$n_risk_F,
    clipping = tr$clipping))
  coefficients <- list(
    imputation_slopes = lapply(groups, `[[`, "slopes"),
    imputation_local_intercepts = if (icpt) stats::setNames(cells$local_intercept, cells$name) else numeric(0),
    has_local_intercepts = icpt,
    selection = sb$coefficients,
    # Counterfactual OR fields keep the existing OR/DR diagnostic names.
    shared_slopes = if (is.null(so)) NULL else so$beta,
    block_slopes = if (reconstruction_separate(so)) lapply(counterfactual_or$blocks, `[[`, "beta") else NULL,
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
    mu = num(b$links, "mu"), T = num(b$links, "T_gl"), C = num(b$links, "C_ol"),
    Bg = num(b$links, "B_gl"), Bo = num(b$links, "B_ol"),
    d_plus = num(b$links, "d_plus"), d_minus = num(b$links, "d_minus"),
    stringsAsFactors = FALSE)
  support <- list(n_full_frame_units = st$N, n_pair_eligible_A_units = st$n_A,
    n_frame_only_units = res$frame_owner$n_frame_only_units,
    design_source = b$design_source, design_columns = b$design_columns, design_rows = st$N,
    selection_design_columns = b$selection_columns, selection_design_rows = st$N,
    prediction_rows = "every F row of each needed cohort (prediction standardization, not fill-only)",
    observation_probability_models = length(sb$blocks), observation_probabilities_fitted = TRUE,
    observation_probability_rows = "every F row of each needed cohort; 1/max(q, 1e-6) evaluated only on observed pairs",
    if_rows = "complete original-F rows, including off-A and frame-only units",
    scale = "core psi = (n_A/N) phi_F with covariance divisor n_A; clusters sum original-unit rows")
  diag_copy <- list(cells = cells,
    clipping = stats::setNames(lapply(b$cells, `[[`, "clipping"), cells$name),
    groups = groups, records = list(selection = selection), counterfactual_or = counterfactual_or,
    treatment = treatment, coefficients = coefficients, fitted_catalog = fitted_catalog, support = support)
  if (unname(cp04_payload(diag_copy)["entries"]) > bound)
    cp04_abort("Augmented diagnostic copy bound underestimated.", "bound")
  out <- c(list(status = "available", method = "augmented", selection_est_method = "dr",
    est_method = b$est_method, imputation_pooling = b$imputation_pooling, outcome_pooling = b$outcome_pooling,
    selection_model = b$selection_model, selection_pooling = b$selection_pooling),
    augmented_public_model_text(b$outcome_recipe$intercept, b$selection_model, b$selection_pooling, b$imputation_pooling, b$outcome_pooling, b$est_method),
    list(fit_structure = paste0("one WLS fit per imputation pooling group (", length(b$groups), " group(s), ",
      length(b$cells), " cell(s)), ", length(sb$blocks), " observation-probability logit block fit(s) over ",
      length(sb$cells), " selection cell(s)",
      if (reconstruction_separate(so)) paste0(", ", so$n_blocks, " comparison-specific counterfactual OR fits on Z") else
        if (!is.null(so)) ", one joint shared counterfactual OR fit on Z" else "",
      if (length(b$treatments)) paste0(", ", length(b$treatments), " deduplicated treatment logit(s)") else ""),
    duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    selection = "conditional_mtar", selection_route = res$selection,
    outcome_recipe = b$outcome_recipe, treatment_recipe = b$treatment_recipe),
    diag_copy,
    list(accounting = b$accounting,
    frame_counts = list(n_full_frame_units = res$reporting$inference$n_total,
      n_pair_eligible_A_units = res$frame_owner$n_a_units),
    auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
    interpretation = if (identical(b$est_method, "dr"))
      paste("Augmented observation correction with treatment DR: each needed cell uses",
        "Z = j + S(Delta - j)/q, combining the actual-change imputation fit (observed pairs, fixed weights) with the",
        "observation-probability fit (every F row). Complete F-row influence functions include the observation-probability,",
        "imputation, repeat-rate, treatment-propensity and shared counterfactual OR paths and both normalizers.",
        "T - C and T - B_g are internal algebra, not reported routes. No numerical certificate.") else
      if (identical(b$est_method, "ipw"))
      paste("Augmented observation correction with treatment IPW: each needed target and control cell uses",
        "Z = j + S(Delta - j)/q, combining the actual-change imputation fit (observed pairs, fixed weights) with the",
        "observation-probability fit (every F row); no counterfactual OR is fitted (outcome_pooling is inactive).",
        "Complete F-row influence functions include the observation-probability, imputation, repeat-rate and",
        "treatment-propensity paths and both normalizers. No numerical certificate.") else
      paste("Augmented observation correction with treatment OR (est_method = \"reg\"): each needed cell uses",
        "Z = j + S(Delta - j)/q, combining the actual-change imputation fit (observed pairs, fixed weights) with the",
        "observation-probability fit (every F row); counterfactual control outcome-change regressions are fitted on Z,",
        "and no treatment logit is fitted. Complete F-row influence functions include the observation-probability,",
        "imputation, repeat-rate and counterfactual OR paths and the treated normalizer. No numerical certificate.")))
  cp04_public_diagnostic_check(out, b, bound)
}
