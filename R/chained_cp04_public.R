# Narrow public Corollary1 adapter. The supplied roster declares the population;
# application assumptions are recorded, never estimated from fit diagnostics.
cp04_validate_declaration <- function(universe, rho, justification) {
  chained_request_need(chained_request_chr1(universe) && universe %in% c("F", "B", "A"),
    '`selection_universe` must be "F", "B" (at least once), or "A" (repeat).')
  chained_request_need(is.null(rho) || chained_request_chr1(rho),
    "`selection_rho` must be NULL or a column name for a fixed KNOWN normalized population-law density ratio.")
  chained_request_need(is.null(justification) || chained_request_chr1(justification),
    "`selection_justification` must be NULL or a nonempty explanation of representativeness/known transport.")
  if (universe == "F") chained_request_need(is.null(rho), "F uses rho = 1; a transport column is only supported for A/B.")
  if (universe != "F") chained_request_need(chained_request_chr1(justification),
    paste0("U=", universe, " requires explicit `selection_justification`: NULL rho asserts U independent of X within cohort; ",
      "otherwise supply a fixed known normalized rho=f_P(X|G,U)/f_P(X|G). U-conditional logit correctness is a separate assumption."))
  invisible(TRUE)
}

# CP08 common-boundary performance: a one-sided formula whose every term is one
# plain variable name (no function call, interaction, offset or backquoted
# label) over plain numeric, non-object, dimensionless columns on both the frame
# and the data has a model matrix made of the intercept column and those raw
# columns coerced to double. On that route the within-unit/agreement check of
# the panel rows is exact on the raw values themselves, so the N-row frame
# model matrix is built once and no NT-row model matrix or row copy is made.
# Everything else (context-dependent calls such as scale()/poly(), factors,
# logical columns, matrix columns, backquoted names) keeps the original route.
cp04_design_direct_columns <- function(tf, vars, frame, data) {
  labels <- attr(tf, "term.labels")
  if (length(attr(tf, "offset")) || !identical(as.integer(attr(tf, "response")), 0L)) return(FALSE)
  if (length(labels) != length(vars) || anyDuplicated(labels) || !setequal(labels, vars)) return(FALSE)
  variables <- as.list(attr(tf, "variables"))[-1L]
  if (length(variables) != length(vars) || !all(vapply(variables, is.symbol, logical(1L)))) return(FALSE)
  plain <- function(col) is.numeric(col) && !is.object(col) && is.null(dim(col))
  for (v in vars) {
    if (!plain(frame[[v]])) return(FALSE)
    if (v %in% names(data) && !plain(data[[v]])) return(FALSE)
  }
  TRUE
}

# Evaluate and validate formulas ONLY on the declared fitting roster. No missing
# X/rho outside U is evaluated and subsequently multiplied by zero.
# `row_frame_pos` may carry the request layer's already computed frame position
# of every data row (match(id, frame_ids)); standalone callers leave it NULL.
cp04_public_design <- function(formula, frame, data, frame_ids, id, U, label, allow_zero_col = FALSE,
                               row_frame_pos = NULL) {
  chained_request_need(inherits(formula, "formula") && length(formula) == 2L,
    paste0("`", label, "` must be a one-sided formula."))
  vars <- all.vars(formula)
  chained_request_need(all(vars %in% names(frame)), paste0("Covariates in `", label, "` must all be supplied on the unit frame."),
    "cdid_chained_request_covariate_error")
  tf <- stats::terms(formula)
  has_intercept <- identical(as.integer(attr(tf, "intercept")), 1L)
  term_labels <- attr(tf, "term.labels")
  active <- if (all(U)) frame else frame[U, , drop = FALSE]
  mm <- function(d) {
    m <- stats::model.matrix(tf, stats::model.frame(tf, d, na.action = stats::na.pass))
    asg <- as.integer(attr(m, "assign"))
    if (ncol(m) > 0L && length(term_labels) > 0L && length(asg) == ncol(m)) {
      cn <- colnames(m)
      for (k in seq_along(term_labels)) {
        idx <- which(asg == k)
        lbl <- term_labels[k]
        if (length(idx) == 1L && nchar(lbl) >= 3L && startsWith(lbl, "`") && endsWith(lbl, "`")) {
          raw_nm <- substr(lbl, 2L, nchar(lbl) - 1L)
          if (raw_nm %in% names(d) && !(raw_nm %in% cn[-idx])) cn[idx] <- raw_nm
        }
      }
      colnames(m) <- cn
    }
    attr(m, "intercept") <- has_intercept
    attr(m, "assign") <- asg
    m
  }
  z <- mm(active)
  chained_request_need(nrow(z) == sum(U) && (if (isTRUE(allow_zero_col)) ncol(z) >= 0L else ncol(z) > 0L) && all(is.finite(z)),
    paste0("Covariates in `", label, "` contain missing or nonfinite values on the fitting roster; NA covariates are refused."),
    "cdid_chained_request_covariate_error")
  Z <- matrix(NA_real_, nrow(frame), ncol(z), dimnames = list(NULL, colnames(z)))
  Z[U, ] <- z
  present <- vars %in% names(data)
  if (length(vars) && all(present)) {
    if (is.null(row_frame_pos)) row_frame_pos <- match(id, frame_ids)
    ii <- which(U[row_frame_pos])
    if (cp04_design_direct_columns(tf, vars, frame, data)) {
      # Raw-value check by unit: the model matrix of a direct formula is the
      # intercept plus these columns as doubles, so equality of the raw values
      # on every panel row is the original comparison without the NT-row matrix.
      agree <- TRUE
      fpos <- row_frame_pos[ii]
      for (v in vars) {
        dv <- data[[v]][ii]
        if (!all(is.finite(dv)) || !all(dv == frame[[v]][fpos])) { agree <- FALSE; break }
      }
      chained_request_need(agree,
        paste0("Covariates in `", label, "` vary within unit across time, have missing values, or disagree with `frame` by id."),
        "cdid_chained_request_covariate_error")
    } else {
      dz <- mm(data[ii, , drop = FALSE])
      expected <- Z[row_frame_pos[ii], , drop = FALSE]
      chained_request_need(identical(dim(dz), dim(expected)) && identical(colnames(dz), colnames(expected)) &&
        all(is.finite(dz)) && all(dz == expected),
        paste0("Covariates in `", label, "` vary within unit across time, have missing values, or disagree with `frame` by id."),
        "cdid_chained_request_covariate_error")
    }
  } else if (any(present)) chained_request_abort(paste0("Covariates in `", label, "` are split between data and frame."),
    "cdid_chained_request_covariate_error")
  attr(Z, "intercept") <- has_intercept
  attr(Z, "assign") <- as.integer(attr(z, "assign"))
  Z
}

cp04_selection_record <- function(cp04, backend = NULL, reader = NULL) {
  if (is.null(cp04)) return(list(method = "none", canonical_assumption = "none"))
  if (identical(cp04$ipw_method, "baseline")) return(list(method = "none", canonical_assumption = "none",
    treatment_model = "risk-set treatment propensity logit (or the custom callback) per (pair, g, K, recipe)",
    reference_inputs = "D60 recipe-indexed designs z(C^r) and weights W_g; baseline pair rates per (pair, cohort, recipe)",
    assumptions = paste("Correct synthetic risk-set treatment propensity logits (built-in) under valid baseline observation",
      "reconstruction, chained parallel trends/no anticipation and the fixed reference law P^g of each context."),
    models = if (is.null(backend)) NULL else backend$accounting[intersect(c("n_rate_fits", "n_treatment_fits",
      "n_custom_callbacks", "n_reference_recipes"), names(backend$accounting))]))
  if (identical(cp04$selection_est_method, "imputation"))
    return(imputation_public_selection_record(cp04, backend, reader))
  if (identical(cp04$selection_est_method, "dr"))
    return(augmented_public_selection_record(cp04, backend, reader))
  model_fields <- c("n_selection_fits", "n_selection_cells", "n_selection_blocks", "n_selection_parameters",
    "n_regression_fits", "n_regression_blocks",
    "n_local_intercepts", "n_shared_slopes", "n_supports", "n_treatment_fits")
  if (identical(cp04$dr_method, "baseline")) return(list(method = "none", canonical_assumption = "none",
    outcome_model = "shared one-period outcome-change control-mean regression",
    treatment_model = "risk-set treatment propensity logit",
    duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    assumptions = paste("Either correct shared one-period outcome-change control-mean regression with local block intercepts",
      "or correct synthetic risk-set treatment propensity logits under valid observation reconstruction,",
      "chained parallel trends/no anticipation and fixed-weight population law."),
    models = if (is.null(backend)) NULL else backend$accounting[intersect(model_fields, names(backend$accounting))]))
  if (identical(cp04$or_method, "baseline")) return(list(method = "none", canonical_assumption = "none",
    outcome_model = "shared one-period outcome-change control-mean regression",
    duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    assumptions = "Shared one-period outcome-change control-mean regression with local block intercepts, chained parallel trends/no anticipation and fixed-weight population law.",
    models = if (is.null(backend)) NULL else backend$accounting[intersect(model_fields, names(backend$accounting))]))
  sel_model <- if (is.null(cp04$selection_model)) "pair_logit" else cp04$selection_model
  sel_pool <- if (is.null(cp04$selection_pooling)) "separate" else cp04$selection_pooling
  sel_req <- if (is.null(cp04$selection_requested)) "corollary1" else cp04$selection_requested
  sb_sum <- selection_backend_summary(backend)
  a_rates <- if (!is.null(sb_sum)) sb_sum$empirical_repeat_rates else if (!is.null(reader))
    stats::setNames(as.double(reader$state$a), as.character(reader$state$cohorts)) else NULL
  # D60: repeat rates a_k^{g,U} belong to each recipe's reference law.
  if (!is.null(backend) && !is.null(backend$reference_ref))
    a_rates <- chained_reference_record(chained_backend_reference(backend))$recipe_repeat_rates
  assump_prefix <- if (identical(sel_model, "independent_period_logit")) {
    paste0("Each active pair/cohort requires population conditional mean trend-MAR on F, positivity, used-pair conditional date independence P(S_e=1|X,G) = P(R_s=1|X,G)P(R_t=1|X,G) on F, and correct marginal date logits in F under ", sel_pool, " slope pooling.")
  } else {
    paste0("Each actual pair/cohort requires population mean trend-MAR, positivity, and correct logit conditional on all X in the declared U",
           if (!identical(sel_pool, "separate")) paste0(" under ", sel_pool, " slope pooling.") else ".")
  }
  out <- list(method = "corollary1", canonical_assumption = "conditional_mtar",
    requested_selection = sel_req, selection_model = sel_model, selection_pooling = sel_pool,
    link = "logit", fitting_universe = cp04$universe,
    population = cp04$population, formula = cp04$formula,
    transport = if (is.null(cp04$rho_column)) "rho=1" else paste0("fixed known normalized rho column: ", cp04$rho_column),
    justification = cp04$justification,
    active_cells = if (is.null(sb_sum)) data.frame(cell_id = integer(), kind = character(), cohort = numeric(), label = character(), block_id = integer(), stringsAsFactors = FALSE) else sb_sum$active_cells,
    n_active_cells = if (is.null(sb_sum)) 0L else sb_sum$n_active_cells,
    n_selection_blocks = if (is.null(sb_sum)) 0L else sb_sum$n_selection_blocks,
    n_selection_parameters = if (is.null(sb_sum)) 0L else sb_sum$n_selection_parameters,
    empirical_repeat_rates = a_rates,
    assumptions = paste(assump_prefix,
      if (!is.null(cp04$dr_method)) "Either correct shared one-period outcome-change control-mean regression with local block intercepts or correct synthetic risk-set treatment propensity logits under valid observation reconstruction, PT/no anticipation and the fixed-weight population law must also hold; F logit does not imply U logit." else
      if (!is.null(cp04$or_method)) "Correct shared one-period outcome-change control-mean regression with local block intercepts, PT/no anticipation and the fixed-weight population law must also hold; F logit does not imply U logit." else
      "Synthetic treatment odds, PT/no anticipation and the fixed-weight population law must also hold; F logit does not imply U logit."),
    if_scale = "phi_F on ALL N roster rows; core psi=(n_A/N)*phi_F with 1/n_A covariance normalization",
    unsupported = c("known pair probabilities (including design p=1)", "custom IPW callbacks without an explicit log_weight_gradient", "estimated or unnormalized transport", "probit"),
    models = if (is.null(backend)) NULL else backend$accounting[intersect(model_fields, names(backend$accounting))])
  if (!is.null(cp04$or_method) || !is.null(cp04$dr_method)) {
    out$outcome_model <- "shared one-period outcome-change control-mean regression"
    out$duration_convention <- "stored-calendar transition count L_e = i_t - i_s"
  }
  if (!is.null(cp04$dr_method)) {
    out$treatment_model <- "risk-set treatment propensity logit"
  }
  out
}

# Same auxiliary descriptors as selection=none. Post and auxiliary links are
# fitted jointly so an actual pair/cohort score is shared across all uses.
# `imputation` is the prepared pure-imputation cp04 record (or NULL) and
# `augmented` the prepared augmented-correction record (or NULL); either
# backend receives the same deduplicated union of post and auxiliary views.
# `reference_source` (D60) turns the fixed union of post and both auxiliary
# views into recipes and one sealed reference reader of every route's designs,
# anchored weights and recipe cohort rates; `ipw_method = "baseline"` is the
# dated baseline IPW/custom-callback route (R/chained_reference_routes.R).
# `pt_assumption`/`run_numerical_transfer` (CP08 descriptor context) are the
# request facts the early canonical OR/DR intake predicate needs; NULL (a
# standalone caller) keeps the link-view route.
build_cp04_public_registry <- function(metadata, frame, reader, ledger, ctrl, control_group, anticipation,
                                       retained_pipeline_entries = 0, external_work_ops = 0, or_method = NULL,
                                       automatic_chunk_rows = FALSE, outcome_recipe = NULL,
                                       selection_model = "pair_logit", selection_pooling = "separate",
                                       dr_method = NULL, treatment_columns = NULL, imputation = NULL,
                                       augmented = NULL, reference_source = NULL, ipw_method = NULL,
                                       custom_est_method = NULL, custom_dots = list(), aux_spec = NULL,
                                       pt_assumption = NULL, run_numerical_transfer = NULL) {
  M <- if (is.null(metadata)) 0L else nrow(metadata$catalog$rows)
  # CP08/D78 work ownership: the caller's cumulative prefix as received. The
  # descriptor-batch and descriptor-context work below is added to
  # `external_work_ops` (the prefix the backend fit receives); the composed
  # total returned in `cp04_work_accounting` charges `work_prefix` exactly once.
  work_prefix <- external_work_ops
  # The pipeline's one auxiliary union plan (CP08 scaling action 3) is reused;
  # standalone callers rebuild the same plan here.
  aux <- if (!is.null(aux_spec)) aux_spec else ipw_auxiliary_spec(ledger, get("panel_owner", frame$panel_ref)$all_times,
    get("panel_owner", frame$panel_ref)$cohorts,
    anticipation, control_group)
  # CP08 construction: the actual catalog is validated once before its
  # descriptor batch; every row keeps the single extractor's bound check and
  # preflight, and reads its fields through one invocation-local source.
  descriptors <- if (M > 0L) extract_pipeline_catalog_row_specs(metadata, seq_len(M),
    retained_entries = retained_pipeline_entries, external_work_ops = external_work_ops,
    control = ctrl$owner_control, return_accounting = TRUE) else list(specs = list(), accounting = NULL)
  specs <- descriptors$specs
  if (!is.null(descriptors$accounting)) external_work_ops <- ipw_owner_safe_sum(
    c(external_work_ops, descriptors$accounting$work_ops_bound), ctrl$owner_control$max_work_ops,
    "CP04 registry plus descriptor preparation work")
  descriptors <- NULL
  # CP08 descriptor context: on the eligible complete-roster baseline OR/DR
  # route (decided once from the sealed owners, validated catalog, request
  # facts and the existing recipe checks, before any link view) the fits read
  # one sealed context of all original pair/cohort/key/ID descriptors instead
  # of M + Q heavy link views; NULL keeps the link-view route unchanged.
  context <- if (M > 0L && is.null(imputation) && is.null(augmented) && is.null(reference_source) &&
      is.null(ipw_method) && !is.function(custom_est_method))
    cp04_descriptor_context(metadata, frame, reader, specs, aux, ctrl, or_method = or_method,
      dr_method = dr_method, outcome_recipe = outcome_recipe, treatment_columns = treatment_columns,
      pt_assumption = pt_assumption, run_numerical_transfer = run_numerical_transfer,
      retained_pipeline_entries = retained_pipeline_entries, external_work_ops = external_work_ops) else NULL
  if (!is.null(context)) {
    # The context's fact/descriptor work precedes the fit and must be part of
    # its cumulative work envelope, not an independently reusable full cap.
    external_work_ops <- ipw_owner_safe_sum(c(external_work_ops, context$accounting$work_ops_bound),
      ctrl$owner_control$max_work_ops, "CP04 context preparation plus fitting work")
    # The already-live pipeline reservation includes 2048 slots per original
    # row plus its 32*(p^2+p*C+64) payload allowance; it covers the descriptor
    # context in place of the omitted heavy views throughout the fit. The
    # existing fit metadata header also remains live through construction.
    views <- list(); aux_index <- as.integer(M + seq_len(context$n_auxiliary))
    aux_errors <- vector("list", context$n_auxiliary)
    construction <- list(route = "descriptor_context", heavy_link_views_constructed = 0L,
      descriptor_rows = context$accounting$descriptor_rows, catalog_rows = as.integer(M),
      auxiliary_rows = context$n_auxiliary, date_facts_checked = context$accounting$date_facts_checked,
      context_entries = context$accounting$entries, context_entries_bound = context$accounting$entries_bound,
      context_work_ops_bound = context$accounting$work_ops_bound, link_support_context = NULL)
  } else {
    built <- cp04_public_link_views(frame, specs, aux, ctrl$owner_control, retained_pipeline_entries)
    views <- built$views; aux_index <- built$aux_index; aux_errors <- built$aux_errors
    construction <- list(route = "link_views", heavy_link_views_constructed = length(views),
      descriptor_rows = 0L, catalog_rows = as.integer(M), auxiliary_rows = length(aux$descriptors),
      date_facts_checked = 0L, context_entries = 0, context_entries_bound = 0, context_work_ops_bound = 0,
      link_support_context = built$construction_work)
    built <- NULL
  }
  reference <- NULL
  if (!is.null(reference_source) && length(views)) {
    reference <- build_chained_reference_reader(reference_source, reader, views, M, ctrl$owner_control,
      retained = retained_pipeline_entries, work = external_work_ops)
    if (!is.null(reference$outcome_recipe)) outcome_recipe <- reference$outcome_recipe
    if (is.null(treatment_columns) && !is.null(dr_method)) treatment_columns <- reference$x_columns
  }
  if (!is.null(ipw_method) && is.null(reference))
    cp04_abort("The full-roster baseline IPW route exists only for reference-date recipes.", "unsupported")
  backend <- if (!is.null(context) && !is.null(dr_method))
    dr_public_fit(dr_method, reader, context, ctrl, retained_pipeline_entries, external_work_ops,
      automatic_chunk_rows, outcome_recipe = context$outcome_recipe, treatment_columns = treatment_columns,
      selection_model = selection_model, selection_pooling = selection_pooling, reference = NULL) else
    if (!is.null(context))
    or_public_fit(or_method, reader, context, ctrl, retained_pipeline_entries, external_work_ops,
      automatic_chunk_rows, outcome_recipe = context$outcome_recipe,
      selection_model = selection_model, selection_pooling = selection_pooling, reference = NULL) else
    if (length(views) && !is.null(augmented))
    augmented_public_fit(augmented, reader, views, ctrl, retained_pipeline_entries, external_work_ops,
      automatic_chunk_rows, reference = reference) else if (length(views) && !is.null(imputation))
    imputation_public_fit(imputation, reader, views, ctrl, retained_pipeline_entries, external_work_ops,
      automatic_chunk_rows, reference = reference) else if (length(views) && !is.null(dr_method))
    dr_public_fit(dr_method, reader, views, ctrl, retained_pipeline_entries, external_work_ops,
      automatic_chunk_rows, outcome_recipe = outcome_recipe, treatment_columns = treatment_columns,
      selection_model = selection_model, selection_pooling = selection_pooling, reference = reference) else
    if (length(views) && !is.null(or_method))
    or_public_fit(or_method, reader, views, ctrl, retained_pipeline_entries, external_work_ops,
      automatic_chunk_rows, outcome_recipe = outcome_recipe,
      selection_model = selection_model, selection_pooling = selection_pooling, reference = reference) else
    if (length(views) && !is.null(ipw_method))
    reference_ipw_public_fit(reader, views, ctrl, retained_pipeline_entries, external_work_ops, reference,
      custom_est_method = custom_est_method, custom_dots = custom_dots) else if (length(views) && is.function(custom_est_method))
    # O15 (D64): conditional_mtar custom callbacks with the explicit log-weight gradient.
    fit_cp04_custom_backend(reader, views, ctrl$owner_control, ctrl$glm_control,
      retained_external_entries = max(0, retained_pipeline_entries - reader$accounting$borrowed_owner_entries),
      external_work_ops = external_work_ops, selection_model = selection_model, selection_pooling = selection_pooling,
      reference = reference, custom_est_method = custom_est_method, custom_dots = custom_dots) else if (length(views))
    fit_cp04_corollary1(reader, views, ctrl$owner_control, ctrl$glm_control,
      retained_external_entries = max(0, retained_pipeline_entries - reader$accounting$borrowed_owner_entries),
      external_work_ops = external_work_ops,
      selection_model = selection_model, selection_pooling = selection_pooling, reference = reference) else NULL
  points <- if (is.null(backend)) list() else if (!is.null(context)) lapply(seq_along(context$descriptors), function(j) {
    # The same point record from the row's own pair descriptor and the one
    # cohort descriptor every row shares (the view route's per-view copies).
    d <- context$descriptors[[j]]
    list(link_id = d$link_id, estimate = chained_full_link_estimate(backend, j),
      cohorts = context$cohorts, pair = d$pair, cp04_index = j)
  }) else lapply(seq_along(views), function(j) {
    v <- views[[j]]
    list(link_id = v$link_id, estimate = chained_full_link_estimate(backend, j),
      cohorts = v$cohorts, pair = v$pair, cp04_index = j)
  })
  context <- NULL
  work_accounting <- cp04_registry_work_accounting(work_prefix, external_work_ops, backend, ctrl$owner_control$max_work_ops)
  list(cp04 = backend, cp04_reader = reader, cp04_aux = aux, cp04_aux_index = aux_index,
    cp04_aux_errors = aux_errors, point_links = points[seq_len(M)], cp04_points = points,
    cp04_construction = construction, cp04_work_accounting = work_accounting,
    nuisance_keys = if (is.null(backend)) character() else
    if (inherits(backend, "cdid_augmented_backend")) augmented_public_nuisance_keys(backend) else
    if (inherits(backend, "cdid_imputation_backend")) imputation_public_nuisance_keys(backend) else
    if (inherits(backend, "cdid_cp04_custom_backend"))
      c(paste0("cp04:", backend$universe, ":selection:", seq_along(backend$selections)),
        paste0("cp04:", backend$universe, ":callback:", seq_along(backend$treatments))) else
    if (inherits(backend, "cdid_reference_ipw_backend"))
      c(paste0("IPW:", backend$law, ":rate:", seq_along(backend$supports)),
        paste0("IPW:", backend$law, if (isTRUE(backend$custom)) ":callback:" else ":treatment:",
          seq_along(backend$treatments))) else if (!is.null(dr_method))
      c(paste0("DR:", backend$law, ":shared_regression:1"),
        paste0("DR:", backend$law, ":block:", seq_along(backend$regressions)),
        paste0("DR:", backend$law, ":treatment:", seq_along(backend$treatments)),
        if (length(backend$selections)) paste0("cp04:", backend$universe, ":selection:", seq_along(backend$selections))) else
    if (!is.null(or_method))
      c(paste0("OR:", backend$law, ":shared_regression:1"),
        paste0("OR:", backend$law, ":block:", seq_along(backend$regressions)),
        if (length(backend$selections)) paste0("cp04:", backend$universe, ":selection:", seq_along(backend$selections))) else
      paste0("cp04:", backend$universe, ":selection:", seq_along(backend$selections)),
    nuisance_views = list(), fitted_views = list())
}

# ---------------------------------------------------------------------------
# CP08/D78 work ownership at the registry seam (one proved accounting fix).
#
#   prefix   : the pipeline's cumulative work handed to build_cp04_public_registry
#   prepared : prefix + descriptor-batch work + descriptor-context work, i.e. the
#              `external_work_ops` the backend fit actually received
#   backend  : the fitted backend. Baseline OR/DR (fit_or_baseline,
#              fit_dr_baseline) declare `accounting$work_ops_convention =
#              "cumulative"`: work_ops_bound = prepared + own_work_ops, and every
#              downstream consumer of those backends reads that cumulative bound.
#              Corollary 1 / selection / imputation / augmented / custom / reference
#              backends keep their own conventions (Corollary 1 records own work
#              only) and are NOT rewritten here.
#
# Composition rule (charges the prefix exactly once, never subtracts a prefix
# from a route whose convention is not declared, never bypasses a finite cap):
#   cumulative backend : total = backend work_ops_bound
#                        (verified: recorded external_work_ops == prepared and
#                         prepared + own_work_ops == work_ops_bound; a mismatch
#                         is an integrity refusal, not a guess)
#   proved incremental : total = prepared + backend work_ops_bound (including
#                        descriptor/context preparation exactly once)
#   unresolved wrapper : unchanged legacy prefix + backend work_ops_bound
#   no backend         : total = prepared
# Before this fix the pipeline added `prefix` to the cumulative baseline bound,
# charging the whole upstream prefix twice.
# ---------------------------------------------------------------------------
cp04_registry_work_accounting <- function(prefix, prepared, backend, max_work_ops) {
  prefix <- ipw_owner_int_bound(prefix, "external_work_ops", min_val = 0)
  prepared <- ipw_owner_int_bound(prepared, "prepared external_work_ops", min_val = 0)
  if (prepared < prefix)
    cp04_abort("CP04 registry preparation work cannot be below the pipeline prefix it extends.", "integrity")
  registry_preparation <- prepared - prefix
  if (is.null(backend)) {
    return(list(prefix = prefix, registry_preparation = registry_preparation, prepared = prepared,
      backend_convention = "none", backend_work_ops_bound = 0, backend_own_work_ops = 0,
      composed_total = ipw_owner_safe_sum(prepared, max_work_ops, "CP04 registry preparation work"),
      rule = "no backend: total = prefix + registry preparation"))
  }
  acc <- backend$accounting
  if (!is.list(acc)) acc <- list()
  # A legacy backend without a recorded bound composed as `c(prefix, NULL)`
  # before; that unchanged behaviour is kept and labelled, not silently guessed.
  bound <- if (is.null(acc$work_ops_bound)) 0 else
    ipw_owner_int_bound(acc$work_ops_bound, "backend work_ops_bound", min_val = 0)
  convention <- acc$work_ops_convention
  if (is.null(acc$work_ops_bound)) convention <- "legacy_no_recorded_work_ops_bound"
  # Most-specific source-proved selection subclasses must precede their
  # Corollary1 parents: their stored bounds exclude external work. Sources:
  # chained_selection_consumers.R:623-664,1230-1284;
  # chained_cp04_corollary1.R:538-549; chained_dr_backend.R:2012-2092.
  known_incremental <- inherits(backend, c("cdid_cp04_corollary1",
    "cdid_or_selection_backend", "cdid_dr_corollary1"))
  # The older direct-pair OR parent is different: its bound includes external
  # work (chained_or_corollary1.R:279-301), unlike its selection subclass.
  legacy_cumulative_or <- !known_incremental && inherits(backend, "cdid_or_corollary1")
  allowed_cumulative <- inherits(backend, c("cdid_or_baseline", "cdid_dr_baseline")) || legacy_cumulative_or
  if (identical(convention, "cumulative") ||
      (is.null(convention) && legacy_cumulative_or && !is.null(acc$work_ops_bound))) {
    if (!allowed_cumulative)
      cp04_abort("The backend's source-proved work convention is not cumulative.", "integrity")
    # Keep the established explicit own-work validator before integrity checks.
    own <- if (is.null(acc$own_work_ops) && legacy_cumulative_or) NULL else
      ipw_owner_int_bound(acc$own_work_ops, "backend own_work_ops", min_val = 0)
    recorded_external <- ipw_owner_int_bound(acc$external_work_ops, "backend external_work_ops", min_val = 0)
    if (is.null(own)) own <- ipw_owner_int_bound(bound - recorded_external, "backend own_work_ops", min_val = 0)
    if (recorded_external != prepared || ipw_owner_safe_sum(c(recorded_external, own)) != bound)
      cp04_abort("Cumulative work accounting of the backend disagrees with the registry's prepared prefix.", "integrity")
    convention <- "cumulative"
    total <- ipw_owner_safe_sum(bound, max_work_ops, "CP04 fitting plus downstream work (cumulative backend)")
    rule <- "cumulative backend: total = backend work_ops_bound = prefix + registry preparation + own work (prefix charged once)"
  } else if (known_incremental && !is.null(acc$work_ops_bound) &&
             (is.null(convention) || identical(convention, "incremental"))) {
    if (!is.null(acc$external_work_ops) &&
        ipw_owner_int_bound(acc$external_work_ops, "backend external_work_ops", min_val = 0) != prepared)
      cp04_abort("Incremental backend was fitted with a different prepared prefix.", "integrity")
    convention <- "incremental"
    own <- bound
    total <- ipw_owner_safe_sum(c(prepared, own), max_work_ops, "CP04 fitting plus downstream work")
    rule <- "source-proved incremental backend: total = prefix + registry preparation + own work (each charged once)"
  } else {
    if (is.null(convention)) convention <- "legacy_unlabelled"
    own <- bound
    total <- ipw_owner_safe_sum(c(prefix, bound), max_work_ops, "CP04 fitting plus downstream work")
    rule <- "unresolved wrapper: unchanged legacy composition = prefix + backend work_ops_bound; no new convention inferred"
  }
  list(prefix = prefix, registry_preparation = registry_preparation, prepared = prepared,
    backend_convention = convention, backend_work_ops_bound = bound, backend_own_work_ops = own,
    composed_total = total, rule = rule)
}

# CP08 stage 5: the registry's catalog and auxiliary link views, in the original
# order, with the original typed auxiliary support-failure capture. One bounded
# construction context, admitted before its creation against the registry's
# envelope (ipw_owner_link_support_context), lets a view borrow the exact
# support summary of an equal-mask view scanned earlier in this invocation; it
# is dropped before any fit. Not admitted, every view scans as before. The
# actual scans and borrowed records are counted apart, as evidence only.
cp04_public_link_views <- function(frame, specs, aux, control, retained_pipeline_entries = 0) {
  M <- length(specs)
  # Only descriptor lengths are read; no pair/date/population scan. Include
  # auxiliary classes before admission so an earlier cache cannot coexist
  # with an unadmitted later constructor. Original failure order stays below.
  max_cohorts <- 2
  for (s in specs) max_cohorts <- max(max_cohorts, 1 + length(s$K))
  for (s in aux$descriptors) max_cohorts <- max(max_cohorts, 1 + length(s$K))
  ctx <- ipw_owner_link_support_context(frame, M + length(aux$descriptors), retained_pipeline_entries, control,
    max_cohorts = max_cohorts)
  views <- lapply(seq_len(M), function(j) {
    s <- specs[[j]]
    construct_ipw_link_view(frame, s$s, s$t, g = s$g, K = s$K,
      control = control, link_id = paste0("cp04:catalog:", j), functional_key = s$functional_key,
      construction_context = ctx)
  })
  aux_index <- rep(NA_integer_, length(aux$descriptors)); aux_errors <- vector("list", length(aux_index))
  for (j in seq_along(aux_index)) {
    s <- aux$descriptors[[j]]
    v <- tryCatch(construct_ipw_link_view(frame, s$s, s$t, g = s$g, K = s$K,
      control = control, link_id = paste0("cp04:auxiliary:", j),
      functional_key = paste0("cmpcat:CP04-aux:g", s$g, ":K", paste(s$K, collapse = ",")),
      construction_context = ctx),
      cdid_ipw_support_error = identity)
    if (inherits(v, "error")) aux_errors[[j]] <- v else {
      views[[length(views) + 1L]] <- v; aux_index[j] <- length(views)
    }
  }
  work <- if (is.null(ctx)) NULL else list(support_scans = ctx$link_scans, support_reuses = ctx$link_reuses,
    support_records = ctx$link_records, reserved_entries = ctx$link_entries)
  ctx <- NULL
  list(views = views, aux_index = aux_index, aux_errors = aux_errors, construction_work = work)
}

# ---------------------------------------------------------------------------
# CP08 descriptor context: one early canonical OR/DR comparison intake.
#
# On the complete-roster baseline OR/DR route every one of the M catalog rows
# and Q auxiliary comparisons is a link view of the same target cohort against
# the never-treated cohort alone on two complete dates, and the fits read only
# its pair descriptor, cohorts, link ID, functional key and (DR) complete-date
# support counts before feeding every original raw pair into the existing
# common fit. This predicate decides that route ONCE, before any link view,
# from actual facts and the existing checks: the estimator-independent sealed
# owner/catalog facts of `ipw_pipeline_compact_basis_facts` (one operator whose
# target is one finite cohort and whose pool is the never-treated cohort, the
# supplied roster exactly the A roster, one fixed n_A x p double design with
# unit weights, PT-post without numerical transfer, every catalog date posting
# the whole A roster), the reader's F = A roster and one never-treated cohort,
# the public default/explicit outcome recipe and the DR treatment recipe through
# their existing checks (an intercept recipe), every auxiliary descriptor on the
# same target/pool with both dates complete, and both cohorts with positive A
# count and mass. It then keeps ONE immutable context: the sealed owner/
# metadata/reader bindings, the one cohort descriptor and every row's pair
# descriptor from the owner's own constructors, every original link ID and
# functional key in registry order, the complete-date support facts (every
# risk unit observed; under unit weights the scanned pair masses are exact
# integer sums equal to the cohort masses) and the checked recipes. Nothing is
# estimated or labelled; NULL keeps the link-view route (a declined request
# raises its original refusal there in its original order), and a refusal
# inside the consuming fit is never caught to retry the other route. The whole
# reservation (context payload, fact-read scratch, descriptor and fact work) is
# admitted against the unchanged caps beside the caller's live envelope before
# any posting is read or any descriptor allocated.
# ---------------------------------------------------------------------------
cp04_descriptor_context <- function(metadata, frame, reader, specs, aux, ctrl, or_method = NULL, dr_method = NULL,
                                    outcome_recipe = NULL, treatment_columns = NULL, pt_assumption = NULL,
                                    run_numerical_transfer = NULL, retained_pipeline_entries = 0,
                                    external_work_ops = 0) {
  # 1. Route shape and sealed bindings (no owner read).
  method <- if (identical(dr_method, "baseline")) "dr" else if (identical(or_method, "baseline")) "or" else return(NULL)
  if (!is.list(specs) || !length(specs) || !is.list(aux) || !is.list(aux$descriptors) ||
      !is.list(metadata) || !inherits(metadata, "cdid_pipeline_metadata") || !is.environment(metadata$ref) ||
      !is.list(frame) || !inherits(frame, "cdid_ipw_frame_owner") || !ipw_owner_seal_ok(frame, "frame_owner") ||
      !is.list(reader) || !inherits(reader, "cdid_cp04_frame_reader") || !ipw_owner_seal_ok(reader, "cp04_reader") ||
      !identical(reader$frame_ref, frame$ref) || !identical(reader$panel_ref, frame$panel_ref) ||
      !is.list(ctrl) || !is.list(ctrl$owner_control) ||
      !gmm_scalar_numeric(retained_pipeline_entries, nonneg = TRUE) || !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    return(NULL)
  own <- tryCatch(validate_ipw_owner_control(ctrl$owner_control), cdid_error = function(e) NULL)
  if (is.null(own)) return(NULL)
  st <- reader$state
  if (isTRUE(st$reference_carrier) || !gmm_scalar_numeric(st$N, positive = TRUE) || !gmm_scalar_numeric(st$n_A) ||
      st$N != st$n_A || !is.double(st$w) || !is.numeric(st$cohorts) || length(st$cohorts) != 2L || anyNA(st$cohorts) ||
      sum(is.infinite(st$cohorts) & st$cohorts > 0) != 1L) return(NULL)
  po <- get0("panel_owner", envir = frame$panel_ref, inherits = FALSE)
  if (!is.list(po) || !inherits(po, "cdid_ipw_panel_owner") || !ipw_owner_seal_ok(po, "panel_owner") ||
      !gmm_scalar_numeric(po$n_calendar_dates, positive = TRUE) || !gmm_scalar_numeric(po$n_a_units, positive = TRUE) ||
      !identical(as.double(po$n_a_units), as.double(st$n_A)) || !is.list(po$postings) ||
      !is.matrix(frame$Z_a) || is.null(colnames(frame$Z_a))) return(NULL)
  M <- length(specs); Q <- length(aux$descriptors)
  for (s in specs) if (!is.list(s) || !gmm_scalar_character(s$functional_key) || !startsWith(s$functional_key, "cmpcat:") ||
      !gmm_scalar_numeric(s$s) || !gmm_scalar_numeric(s$t) || !(s$s < s$t) || !gmm_scalar_numeric(s$g) ||
      !is.numeric(s$K) || is.object(s$K) || !is.null(dim(s$K)) || length(s$K) != 1L) return(NULL)
  for (d in aux$descriptors) if (!is.list(d) || !gmm_scalar_numeric(d$s) || !gmm_scalar_numeric(d$t) || !(d$s < d$t) ||
      !gmm_scalar_numeric(d$g) || !is.numeric(d$K) || is.object(d$K) || !is.null(dim(d$K)) || length(d$K) != 1L) return(NULL)
  # 2. Whole reservation before any posting read or allocation. Payload: one
  # descriptor record per row (five slots, six pair scalars, two bounded ID/key
  # strings, two row indices: at most 320 slots including names/attributes, plus half the key bytes), the
  # one cohort descriptor and support facts (m = 2 cohorts), four date-fact
  # vectors, both recipes and the context headers. Scratch: one posting chunk
  # slice of a date fact, the descriptor temporaries and the completeness memo.
  # Work: every catalog/auxiliary date's complete-posting fact, the facts'
  # row/operator/weight/design scans, two descriptor searches per row, the
  # cohort/recipe checks and the payload verification walk.
  n_a <- as.double(st$n_A); T_cal <- as.double(po$n_calendar_dates); p <- as.double(ncol(frame$Z_a))
  m <- 2; n_desc <- as.double(M + Q); step <- min(n_a, own$max_chunk_rows)
  C_a <- as.double(length(frame$a_cohorts)); J <- as.double(length(st$cohorts))
  key_bytes <- sum(nchar(vapply(specs, `[[`, character(1), "functional_key"), type = "bytes")) + 128 * Q
  col_bytes <- sum(nchar(colnames(frame$Z_a), type = "bytes"))
  D_bound <- min(T_cal, 2 * n_desc)
  cap <- min(own$max_retained_entries, own$max_matrix_entries)
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    max(0, retained_pipeline_entries - reader$accounting$borrowed_owner_entries)
  entries <- ipw_owner_safe_sum(c(4096, 320 * n_desc, ceiling(key_bytes / 2), 64 * m + 512, 8 * D_bound + 64,
    512 + 4 * col_bytes + 32 * p), cap, "CP04 descriptor context payload")
  scratch <- ipw_owner_safe_sum(c(16 * step, 64 * m, 16 * C_a, 1024), own$max_matrix_entries,
    "CP04 descriptor context scratch")
  fact_work <- ipw_owner_safe_product(D_bound, 16 * n_a + 64 * ceiling(n_a / step) + 64,
    label = "CP04 descriptor complete-date posting facts")
  work <- ipw_owner_safe_sum(c(512, fact_work, 64 * n_desc * (2 + ceiling(log2(T_cal + 1))),
    64 * M * (1 + ceiling(log2(M + 1))), 32 * n_a + 16 * p, 128 * (m + C_a + J + p + 16), 8 * entries),
    own$max_work_ops, "CP04 descriptor context work")
  cp04_budget(own, c(base, entries), scratch, c(external_work_ops, work), entries, "CP04 descriptor context")
  # 3. Estimator-independent sealed-owner/catalog facts (reads the catalog dates' postings).
  facts <- ipw_pipeline_compact_basis_facts(metadata, frame, own, pt_assumption, run_numerical_transfer, st$N)
  if (is.null(facts) || !identical(as.integer(facts$n_rows), as.integer(M)) || !identical(facts$frame_ref, frame$ref) ||
      !identical(facts$panel_ref, frame$panel_ref) || !identical(facts$metadata_ref, metadata$ref) ||
      !gmm_scalar_numeric(facts$g) || !identical(facts$K, Inf)) return(NULL)
  g <- as.double(facts$g); K <- Inf
  target <- match(g, st$cohorts); control <- match(Inf, st$cohorts)
  if (is.na(target) || is.na(control) || target == control) return(NULL)
  # 4. The existing recipe checks on the static Z_a columns; an intercept recipe.
  z_cols <- colnames(frame$Z_a)
  recipe <- tryCatch({
    r <- if (is.null(outcome_recipe)) or_public_default_recipe(frame) else outcome_recipe
    if (method == "or") or_pooled_check_recipe(r, z_cols, "baseline") else dr_check_outcome_recipe(r, z_cols, "baseline")$recipe
  }, cdid_error = function(e) NULL)
  if (is.null(recipe) || !isTRUE(recipe$intercept)) return(NULL)
  treatment <- NULL
  if (method == "dr") {
    treatment <- tryCatch(dr_check_treatment_recipe(treatment_columns, z_cols, "baseline"), cdid_error = function(e) NULL)
    if (is.null(treatment)) return(NULL)
  }
  # 5. Every row and auxiliary descriptor: this target against the never-treated
  # cohort alone; both cohorts with positive A count and mass, then the one
  # cohort descriptor by the owner's constructor.
  for (s in specs) if (!identical(as.double(s$g), g) || !identical(as.double(s$K), K)) return(NULL)
  for (d in aux$descriptors) if (!identical(as.double(d$g), g) || !identical(as.double(d$K), K)) return(NULL)
  ci <- ipw_owner_match_sorted_double(c(g, K), frame$a_cohorts)
  if (anyNA(ci) || any(frame$cohort_a_count[ci] < 1L) ||
      any(!is.finite(frame$cohort_a_weight[ci]) | frame$cohort_a_weight[ci] <= 0)) return(NULL)
  cohorts <- ipw_owner_make_cohort_descriptor(frame, g, K, require_support = TRUE)
  # 6. Pair descriptors by the owner's constructor (the link constructor's own
  # statement) and complete-date facts of every auxiliary date not among the
  # catalog dates, on this sealed panel's postings.
  post <- po$postings
  memo <- list2env(list(panel_owner = po), parent = emptyenv())
  dates <- as.integer(facts$date_indices)
  if (anyNA(dates) || !length(dates)) return(NULL)
  desc <- vector("list", M + Q)
  for (j in seq_len(M)) {
    s <- specs[[j]]
    pr <- ipw_owner_make_pair_descriptor(po, s$s, s$t)
    if (is.na(pr$s_calendar_index) || is.na(pr$t_calendar_index) ||
        !all(c(pr$s_calendar_index, pr$t_calendar_index) %in% dates)) return(NULL)
    desc[[j]] <- list(pair = pr, link_id = paste0("cp04:catalog:", j), functional_key = s$functional_key,
      catalog_row = as.integer(j), auxiliary_row = NA_integer_)
  }
  for (j in seq_len(Q)) {
    d <- aux$descriptors[[j]]
    pr <- ipw_owner_make_pair_descriptor(po, d$s, d$t)
    if (is.na(pr$s_calendar_index) || is.na(pr$t_calendar_index)) return(NULL)
    for (dd in c(pr$s_calendar_index, pr$t_calendar_index)) {
      if (dd %in% dates) next
      if (!isTRUE(ipw_owner_date_complete(po, dd, step, memo))) return(NULL)
      dates <- c(dates, as.integer(dd))
    }
    desc[[M + j]] <- list(pair = pr, link_id = paste0("cp04:auxiliary:", j),
      functional_key = paste0("cmpcat:CP04-aux:g", d$g, ":K", paste(d$K, collapse = ",")),
      catalog_row = NA_integer_, auxiliary_row = as.integer(j))
  }
  memo <- NULL
  dates <- sort(unique(dates))
  if (!isTRUE(all(post$date_count[dates] == po$n_a_units))) return(NULL)
  n_risk <- sum(cohorts$group_a_count)
  support <- list(n_observed_pair_all = as.integer(po$n_a_units), n_risk_a_units = n_risk,
    n_outside_risk_a_units = as.integer(po$n_a_units) - n_risk, n_observed_risk_units = n_risk,
    n_nonobserved_risk_units = 0L, by_cohort_pair_count = cohorts$group_a_count,
    by_cohort_pair_mass = cohorts$group_a_mass,
    rates = stats::setNames(rep(1, cohorts$m_cohorts), cohorts$cohort_labels))
  context <- structure(list(
    route = "descriptor_context", method = method,
    frame_ref = frame$ref, panel_ref = frame$panel_ref, metadata_ref = metadata$ref, reader_ref = reader$ref,
    n_a = as.integer(po$n_a_units), N = st$N, step = step, g = g, K = K,
    target = as.integer(target), control = as.integer(control), cohorts = cohorts,
    dates = dates, date_start = post$date_start[dates], date_end = post$date_end[dates],
    date_count = post$date_count[dates],
    n_catalog = as.integer(M), n_auxiliary = as.integer(Q), descriptors = desc, support = support,
    outcome_recipe = recipe, treatment_recipe = treatment,
    accounting = list(entries = 0, entries_bound = entries, scratch_entries = scratch, work_ops_bound = work,
      retained_total = base + entries, work_total = external_work_ops + work,
      date_facts_checked = length(dates), descriptor_rows = as.integer(M + Q), heavy_link_views_constructed = 0L),
    scope = paste("sealed-owner/catalog/request facts and the existing recipe checks only; every original pair,",
      "cohort, link ID and functional key kept; no estimate, label or outcome read")),
    class = "cdid_cp04_descriptor_context")
  context["ref"] <- list(NULL)
  actual <- unname(cp04_payload(context)["entries"])
  if (actual > entries) cp04_abort("CP04 descriptor context bound underestimated.", "bound")
  context$accounting$entries <- actual
  ipw_owner_seal(context, "cp04_descriptor_context")
}

cp04_backend_domain <- function(backend) {
  if (inherits(backend, "cdid_ipw_raw_backend")) {
    validate_ipw_raw_backend(backend)
    frame <- get("frame_owner", backend$frame_ref, inherits = FALSE)
    panel <- get("panel_owner", backend$panel_ref, inherits = FALSE)
    n_a <- as.integer(frame$n_a_units)
    return(list(
      is_raw = TRUE,
      reader = NULL,
      state = NULL,
      frame = frame,
      panel = panel,
      frame_ref = backend$frame_ref,
      panel_ref = backend$panel_ref,
      N = n_a,
      n_A = n_a,
      lambda = 1,
      id = panel$a_id
    ))
  }
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  list(
    is_raw = FALSE,
    reader = reader,
    state = st,
    frame = frame,
    panel = panel,
    frame_ref = frame$ref,
    panel_ref = reader$panel_ref,
    N = st$N,
    n_A = st$n_A,
    lambda = st$lambda,
    id = cp04_retained_if_ids(backend, reader)
  )
}

# Only the public automatic Corollary 1 OR route may reduce its row batch.
# The existing affine charges determine admission before any IF is evaluated;
# explicit public/private chunk sizes and the link blocks remain fixed.
cp04_plan_rows <- function(backend, link_sets, retained = 0, chunk_rows = NULL) {
  st <- cp04_backend_domain(backend); ctrl <- backend$control
  if (!is.null(chunk_rows)) {
    chunk_rows <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
    if (chunk_rows > ctrl$max_chunk_rows)
      cdid_ipw_owner_abort("Replay rows exceed the explicit chunk bound.", "cdid_ipw_owner_bound_error")
  }
  B <- min(st$N, if (is.null(chunk_rows)) ctrl$max_chunk_rows else chunk_rows)
  if (!inherits(backend, c("cdid_or_corollary1", "cdid_dr_corollary1", "cdid_dr_baseline", "cdid_imputation_backend",
                           "cdid_augmented_backend")) ||
      !isTRUE(backend$accounting$automatic_chunk_rows) || !is.null(chunk_rows)) return(B)
  env_fn <- if (inherits(backend, "cdid_augmented_backend")) augmented_if_envelope else
    if (inherits(backend, "cdid_imputation_backend")) imputation_if_envelope else
    if (inherits(backend, "cdid_dr_corollary1")) dr_corollary1_if_envelope else
    if (inherits(backend, "cdid_dr_baseline")) dr_baseline_if_envelope else or_corollary1_if_envelope
  own <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  for (links in link_sets) {
    one <- env_fn(backend, 1, links)
    two <- env_fn(backend, 2, links)
    slope <- two$output + two$scratch - one$output - one$scratch
    intercept <- one$output + one$scratch - slope
    output_slope <- two$output - one$output
    output_intercept <- one$output - output_slope
    B <- min(B, floor((ctrl$max_matrix_entries - retained - own - intercept) / slope),
      floor((ctrl$max_output_entries - output_intercept) / output_slope))
  }
  if (B < 1) {
    # Preserve the ordinary typed refusal and its complete one-row charge.
    pre <- env_fn(backend, 1, link_sets[[1L]])
    ipw_owner_safe_sum(c(retained, own, pre$output, pre$scratch), ctrl$max_matrix_entries,
      "CP04 producer and downstream simultaneous entries")
    ipw_owner_safe_sum(pre$output, ctrl$max_output_entries, "Corollary 1 OR IF chunk output payload")
    cdid_ipw_owner_abort("No admissible Corollary 1 OR row batch.", "cdid_ipw_owner_bound_error")
  }
  B
}

cp04_preflight_links <- function(backend, links, retained = 0, chunk_rows = NULL) {
  dom <- cp04_backend_domain(backend); ctrl <- backend$control
  B <- cp04_plan_rows(backend, list(links), retained, chunk_rows)
  pre <- chained_full_if_budget(backend, B, links)
  own <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  ipw_owner_safe_sum(c(retained, own, pre$output, pre$scratch), ctrl$max_matrix_entries,
    "CP04 producer and downstream simultaneous entries")
  work <- ipw_owner_safe_sum(ceiling(dom$N / B) * pre$work, ctrl$max_work_ops,
    "CP04 requested replay cumulative work")
  list(producer_workspace_entries = own + pre$output + pre$scratch,
    peak_concurrent_entries = own + pre$output + pre$scratch, work_ops_bound = work,
    chunk_rows = B)
}

# ---------------------------------------------------------------------------
# CP08 retained moment IF (invocation-local). For ONE supported sealed
# CP04 backend and ONE link block, the ordinary
# producer's raw scale = "core" complete link IF (every N roster row, the
# requested links in their order, the original global beta and propensity
# contributions) is written once into an N x L matrix during the covariance
# mean pass and sliced, in the original row partition, for the covariance Gram
# pass and both GMM projection passes. Nothing projected or centered is ever
# stored; nothing is attached to the backend, reader, frame, panel or any
# global. The store reaches a producer only through the classed carrier that
# `cp04_retained_if_bind` returns in the backend position, which the pipeline
# creates after admission and drops right after the GMM solve. Admission is
# decided from exact facts and arithmetic BEFORE any allocation; a declined
# store leaves every original charge, partition and stream untouched.
# ---------------------------------------------------------------------------

cp04_supported_backend_classes <- c(
  "cdid_or_baseline", "cdid_dr_baseline",
  "cdid_cp04_corollary1", "cdid_or_corollary1", "cdid_dr_corollary1",
  "cdid_imputation_backend", "cdid_augmented_backend",
  "cdid_cp04_custom_backend", "cdid_reference_ipw_backend",
  "cdid_ipw_raw_backend"
)

# Charged entries: the N x L store plus one possible duplicate, 2N ID/row
# slots, 4B x L slice scratch, 4 L^2 Gram-side temporaries and 256 headers.
cp04_retained_if_entries <- function(N, L, B) {
  ipw_owner_safe_sum(c(ipw_owner_safe_product(2, ipw_owner_safe_product(N, L)),
    ipw_owner_safe_product(2, N), ipw_owner_safe_product(4, ipw_owner_safe_product(B, L)),
    ipw_owner_safe_product(4, ipw_owner_safe_product(L, L)), 256))
}
# Copy/slice work; route-specific metadata and boundary validation are added below.
cp04_retained_if_work <- function(N, L) ipw_owner_safe_product(128, ipw_owner_safe_product(N, L))

# Only a finite exact requested count above an explicit cap proves capacity.
# Source03's bound class also covers invalid counts, exact overflow and real
# matrix dimension errors; those must retain their original failure category.
cp04_capacity_decline <- function(e) {
  d <- e$diagnostics
  count <- function(x) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) == 1L && is.finite(x) && x >= 0 && x == floor(x) && x <= 2^53 - 1
  if (inherits(e, "cdid_ipw_owner_bound_error") &&
      !inherits(e, c("cdid_ipw_owner_control_error", "cdid_ipw_owner_matrix_dim_error")) &&
      is.list(d) && count(d$requested) && count(d$bound) && d$requested > d$bound) return(NULL)
  stop(e)
}

# Augmented dated fits own a sealed recipe-layout reader. Its repeated IDs
# belong to selection evaluation; the producer emits lay$frame_id in original
# unit order. Pure imputation keeps the original reader and a separate layout.
cp04_retained_if_ids <- function(backend, reader) {
  if (inherits(backend, "cdid_ipw_raw_backend"))
    get("panel_owner", backend$panel_ref, inherits = FALSE)$a_id else
  if (inherits(backend, "cdid_augmented_backend") && !is.null(backend$layout))
    backend$layout$frame_id else reader$state$id
}

# Snapshots borrow sealed objects already charged by the producer. Charge their
# bindings/headers, link labels, ID slices and validation temporaries separately.
# Complete-support posting checks are linear in dates; every other owner check
# uses the sealed SEXP identity, without walking the owner's arrays.
cp04_retained_if_boundary_cost <- function(backend, links, B) {
  dom <- cp04_backend_domain(backend)
  frame <- dom$frame
  n <- ipw_owner_safe_sum(c(length(backend), length(dom$reader), length(dom$state),
    length(links), length(backend$complete_support$dates)))
  bytes <- ipw_owner_safe_sum(vapply(backend$links[links], function(l)
    ipw_owner_safe_sum(nchar(l$link_id, type = "bytes")), numeric(1)))
  metadata <- ipw_owner_safe_sum(c(4096, ipw_owner_safe_product(64, n),
    ipw_owner_safe_product(16, bytes),
    ipw_owner_safe_product(4, ipw_owner_safe_product(B, frame$id_byte_metadata$max))))
  list(entries = metadata, work = ipw_owner_safe_product(64, metadata))
}

cp04_retained_if_validate_backend <- function(backend) {
  if (!is.list(backend) || !inherits(backend, cp04_supported_backend_classes))
    cp04_abort("Retained moment IF storage serves only supported sealed CP04 backends.", "integrity")
  if (inherits(backend, "cdid_ipw_raw_backend")) validate_ipw_raw_backend(backend) else
  if (inherits(backend, "cdid_dr_baseline")) validate_dr_baseline_backend(backend) else
  if (inherits(backend, "cdid_or_baseline")) validate_or_baseline_backend(backend) else
  if (inherits(backend, "cdid_cp04_corollary1")) validate_cp04_backend(backend) else
  if (inherits(backend, "cdid_or_corollary1")) validate_or_corollary1_backend(backend) else
  if (inherits(backend, "cdid_dr_corollary1")) validate_dr_corollary1_backend(backend) else
  if (inherits(backend, "cdid_imputation_backend")) validate_imputation_backend(backend) else
  if (inherits(backend, "cdid_augmented_backend")) validate_augmented_backend(backend) else
  if (inherits(backend, "cdid_cp04_custom_backend")) validate_cp04_custom_backend(backend) else
  if (inherits(backend, "cdid_reference_ipw_backend")) validate_reference_ipw_backend(backend) else
    cp04_abort("Retained moment IF storage serves only supported sealed CP04 backends.", "integrity")
  dom <- cp04_backend_domain(backend)
  ids <- dom$id
  if (length(ids) != dom$N)
    cp04_abort("Retained moment IF storage requires the producer's original-unit roster.", "integrity")
  if (!is.null(backend$selection_backend)) {
    validate_selection_models(backend$selection_backend)
    # For dated augmented fits this is the layout reader, not the original
    # reader: validate the actual selection producer's owner and frame/panel.
    sr <- get("cp04_reader", backend$selection_backend$reader_ref, inherits = FALSE)
    if (!identical(backend$selection_backend$reader_ref, backend$reader_ref) ||
        !identical(sr$frame_ref, dom$reader$frame_ref) || !identical(sr$panel_ref, dom$reader$panel_ref))
      cp04_abort("Retained moment IF storage selection backend does not match its producer reader.", "integrity")
  }
  invisible(TRUE)
}

# Reserve a store or return NULL. Structural facts first (a supported sealed
# CP04 backend, one covariance link block, identical covariance and GMM row
# partitions), then the exact arithmetic admission of the cache inside the
# caller's complete simultaneous live envelope and work cap, then both original
# planners re-run with the cache charged: a planner bound refusal or any
# partition change declines. Only that typed bound refusal of the two probes is
# a decline; integrity and numerical errors propagate. Nothing is allocated here.
cp04_retained_if_reserve <- function(backend, links, chunk_rows, cov_plan, retained, gmm_retained, live, work, control) {
  if (!is.list(backend) || !inherits(backend, cp04_supported_backend_classes)) return(NULL)
  if (!is.list(cov_plan) || !is.list(cov_plan$blocks) || !length(cov_plan$blocks))
    cp04_abort("Invalid covariance plan for retained moment IF storage.", "integrity")
  if (length(cov_plan$blocks) != 1L) return(NULL)
  cp04_retained_if_validate_backend(backend)
  dom <- cp04_backend_domain(backend)
  reader <- dom$reader; st <- dom$state; frame <- dom$frame
  N <- dom$N; L <- length(links)
  control <- validate_ipw_owner_control(control)
  for (nm in c("retained", "gmm_retained", "live", "work"))
    ipw_owner_int_bound(get(nm), nm, min_val = 0)
  ipw_owner_int_bound(chunk_rows, "chunk_rows")
  if (!cp04_reduction_count(N, 1L) || !cp04_reduction_count(L, 1L) ||
      !is.numeric(links) || is.object(links) || !is.null(dim(links)) || any(!is.finite(links)) ||
      any(links != floor(links)) || any(links < 1L | links > length(backend$links)) || anyDuplicated(links) > 0L)
    cp04_abort("Invalid explicit retained moment IF link batch.", "input")
  if (!identical(as.integer(cov_plan$blocks[[1L]]), seq_len(L)) ||
      !isTRUE(cov_plan$chunk_rows == chunk_rows)) return(NULL)
  if (chunk_rows > N || chunk_rows > control$max_chunk_rows)
    cp04_abort("Invalid retained moment IF row partition.", "input")
  ipw_owner_check_matrix_dims(N, L, label = "CP04 retained moment IF")
  boundary <- cp04_retained_if_boundary_cost(backend, links, chunk_rows)
  entries <- ipw_owner_safe_sum(c(cp04_retained_if_entries(N, L, chunk_rows), boundary$entries))
  # Four passes, including both sides of each callback and admission/planner
  # checks. Keep the old producer work in full: saved replays never fund checks.
  checks <- ipw_owner_safe_sum(c(32, ipw_owner_safe_product(8, ceiling(N / chunk_rows) + 1)))
  cache_work <- ipw_owner_safe_sum(c(cp04_retained_if_work(N, L),
    ipw_owner_safe_product(checks, boundary$work)))
  cap <- min(control$max_retained_entries, control$max_matrix_entries)
  total <- ipw_owner_safe_sum(c(live, entries))
  work_total <- ipw_owner_safe_sum(c(work, cache_work))
  if (total > cap || work_total > control$max_work_ops) return(NULL)
  retained_probe <- ipw_owner_safe_sum(c(retained, entries))
  gmm_retained_probe <- ipw_owner_safe_sum(c(gmm_retained, entries))
  probe <- function(expr) tryCatch(expr, cdid_ipw_owner_bound_error = cp04_capacity_decline)
  cov_probe <- probe(cp04_stream_plan(backend, links, retained_probe))
  if (is.null(cov_probe) || length(cov_probe$blocks) != 1L || !isTRUE(cov_probe$chunk_rows == chunk_rows)) return(NULL)
  gmm_probe <- probe(cp04_preflight_links(backend, links, gmm_retained_probe))
  if (is.null(gmm_probe) || !isTRUE(gmm_probe$chunk_rows == chunk_rows)) return(NULL)
  store <- new.env(parent = emptyenv())
  store$kind <- "cp04_retained_if"; store$state <- "reserved"
  store$backend_snapshot <- backend
  store$backend_ref <- backend$ref; store$reader_ref <- backend$reader_ref
  store$frame_ref <- dom$frame_ref; store$panel_ref <- dom$panel_ref
  store$reader_snapshot <- reader; store$state_snapshot <- st; store$frame_snapshot <- frame
  store$panel_snapshot <- dom$panel
  store$reference_ref <- backend$reference_ref
  store$reference_snapshot <- if (is.null(backend$reference_ref)) NULL else chained_backend_reference(backend)
  store$selection_ref <- if (is.null(backend$selection_backend)) NULL else backend$selection_backend$ref
  store$selection_snapshot <- backend$selection_backend
  store$selection_reader_snapshot <- if (is.null(backend$selection_backend)) NULL else
    get("cp04_reader", backend$selection_backend$reader_ref, inherits = FALSE)
  store$design_ref <- backend$design_ref
  store$design_snapshot <- if (is.null(backend$design_ref)) NULL else get("imputation_design_owner", backend$design_ref, inherits = FALSE)
  store$layout_snapshot <- backend$layout
  store$layout_state_snapshot <- if (is.null(backend$layout)) NULL else backend$layout$state
  store$X_F_snapshot <- backend$X_F
  store$pair_outcomes_snapshot <- backend$pair_outcomes
  for (nm in c("reader_snapshot", "state_snapshot", "frame_snapshot", "panel_snapshot",
               "reference_snapshot", "selection_snapshot", "design_snapshot",
               "layout_snapshot", "layout_state_snapshot", "X_F_snapshot", "pair_outcomes_snapshot"))
    lockBinding(nm, store)
  store$id <- dom$id; store$N <- N; store$n_A <- dom$n_A; store$lambda <- dom$lambda
  store$links <- as.integer(links); store$L <- L; store$B <- chunk_rows; store$scale <- "core"
  store$entries <- entries; store$work <- cache_work; store$admission_total <- total; store$cap <- cap
  store$producer_passes <- 0L; store$slice_passes <- 0L; store$batches <- 0L
  store$metadata_entries <- boundary$entries; store$boundary_work <- boundary$work
  for (nm in setdiff(ls(store, all.names = TRUE), c("state", "producer_passes", "slice_passes", "batches")))
    if (!bindingIsLocked(nm, store)) lockBinding(nm, store)
  store
}

# Source-derived simultaneous phase inventories for the two CP04 replay phases:
# - Covariance phase (`cp04_stream_covariance`): persistent pipeline envelope
#   and retained diagnostics (`retained_with_diagnostics`), persistent D5
#   ledger/basis/enclosure arrays (`d5_arrays`), covariance consumer buffers
#   (full-width sweep centering and crossprod slices `6 * B_cov * L`, Gram
#   product/copy `2 * L^2`, mean/estimate/header slots `16 * L + 256` whenever
#   not already in `cov_plan$peak`), complete Gram reserve (`4 * L^2`), backend
#   owned + borrowed entries, and covariance producer output + scratch (`cov_plan$peak`).
# - GMM phase (`solve_chained_gmm_core`, after covariance store release):
#   `cp04_planner_retained` (which includes `combined_initial`,
#   `d5_arrays$retained_diagnostic_entries`, persistent D5 arrays,
#   `gmm_preflight$total_payload_entries` and `final_diagnostic_copy_entries`),
#   retained covariance outputs (`cov_stream_used`, `mhat_used`, `V_used`,
#   `expanded_m`, `prov_used` and its sealed `complete_unit_if_gram`:
#   `4 * L^2 + 16 * L + 2 * M_all + 512`), and GMM producer workspace at the
#   issued GMM partition (`point_preflight$producer_workspace_entries`).
cp04_retained_if_phase_inventory <- function(phase = c("covariance", "gmm"), backend, cov_plan, point_preflight,
                                             retained_with_diagnostics, cp04_planner_retained, d5_arrays,
                                             L_used, M_all) {
  if (!is.character(phase) || length(phase) != 1L || is.na(phase) || !(phase %in% c("covariance", "gmm")))
    cp04_abort("Invalid retained moment IF phase.", "input")
  d5_persistent <- ipw_owner_safe_sum(c(d5_arrays$ledger_if_entries, d5_arrays$basis_if_entries,
    d5_arrays$ledger_enclosure_entries, d5_arrays$basis_enclosure_entries))
  if (identical(phase, "covariance")) {
    B_cov <- cov_plan$chunk_rows
    consumer_extra <- if (is.null(cov_plan$consumer_workspace_entries))
      ipw_owner_safe_sum(c(
        ipw_owner_safe_product(6, ipw_owner_safe_product(B_cov, L_used)),
        ipw_owner_safe_product(2, ipw_owner_safe_product(L_used, L_used)),
        ipw_owner_safe_product(16, L_used), 256)) else 0
    gram_entries <- ipw_owner_safe_product(4, ipw_owner_safe_product(L_used, L_used))
    own <- ipw_owner_safe_sum(c(backend$accounting$owned_entries, backend$accounting$borrowed_entries))
    retained <- ipw_owner_safe_sum(c(retained_with_diagnostics, d5_persistent, consumer_extra))
    live <- ipw_owner_safe_sum(c(retained, own, cov_plan$peak, gram_entries))
    return(list(retained = retained, live = live, consumer_extra = consumer_extra,
      gram_entries = gram_entries, d5_persistent_entries = d5_persistent))
  }
  cov_outputs <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(4, ipw_owner_safe_product(L_used, L_used)),
    ipw_owner_safe_product(16, L_used),
    ipw_owner_safe_product(2, M_all), 512))
  retained <- ipw_owner_safe_sum(c(cp04_planner_retained, cov_outputs))
  live <- ipw_owner_safe_sum(c(retained, point_preflight$producer_workspace_entries))
  list(retained = retained, live = live, cov_retained_outputs = cov_outputs,
    d5_persistent_entries = d5_persistent)
}

# D92 phase-local retained moment IF (invocation-local, single-phase lifetime).
# Considered only when the accepted whole-lifetime store declines for a
# one-block, non-transfer, dated augmented backend:
#   - phase = "covariance": covariance-only store at `cov_plan$chunk_rows`,
#     filled during the covariance mean pass, sliced during the covariance Gram
#     pass, and released before allocating GMM consumers.
#   - phase = "gmm": GMM-only store at `point_preflight$chunk_rows`, filled
#     during the first GMM projection pass, sliced during the second GMM
#     projection pass, and released before target IF transport, auxiliary
#     reporting and final diagnostic copying.
# Never keeps two stores live simultaneously (`active_store` must be NULL or
# already released). Only proved finite capacity refusals decline (NULL);
# integrity, input, malformed, overflow and numerical failures propagate.
cp04_retained_if_phase_reserve <- function(backend, links, phase = c("covariance", "gmm"), chunk_rows, cov_plan,
                                           retained, live, work, control, numerical_transfer = FALSE,
                                           active_store = NULL) {
  if (!is.character(phase) || length(phase) != 1L || is.na(phase) || !(phase %in% c("covariance", "gmm")))
    cp04_abort("Invalid retained moment IF phase.", "input")
  if (!is.logical(numerical_transfer) || length(numerical_transfer) != 1L || is.na(numerical_transfer))
    cp04_abort("Invalid numerical_transfer flag for retained moment IF storage.", "input")
  if (!is.null(active_store) && (!is.environment(active_store) || !identical(active_store$kind, "cp04_retained_if") ||
      !identical(active_store$state, "released")))
    cp04_abort("Retained moment IF phase stores cannot be simultaneously live.", "integrity")
  if (inherits(backend, "cdid_cp04_retained_if_bound"))
    cp04_abort("Retained moment IF phase storage cannot be reserved on an already bound producer.", "integrity")
  if (isTRUE(numerical_transfer)) return(NULL)
  if (!is.list(backend) || !inherits(backend, c("cdid_augmented_backend", "cdid_ipw_raw_backend"))) return(NULL)
  if (!is.list(cov_plan) || !is.list(cov_plan$blocks) || !length(cov_plan$blocks))
    cp04_abort("Invalid covariance plan for retained moment IF storage.", "integrity")
  if (length(cov_plan$blocks) != 1L) return(NULL)
  cp04_retained_if_validate_backend(backend)
  if (inherits(backend, "cdid_augmented_backend") && (is.null(backend$reference_ref) || is.null(backend$layout))) return(NULL)
  dom <- cp04_backend_domain(backend)
  reader <- dom$reader; st <- dom$state; frame <- dom$frame
  N <- dom$N; L <- length(links)
  control <- validate_ipw_owner_control(control)
  for (nm in c("retained", "live", "work"))
    ipw_owner_int_bound(get(nm), nm, min_val = 0)
  ipw_owner_int_bound(chunk_rows, "chunk_rows")
  if (!cp04_reduction_count(N, 1L) || !cp04_reduction_count(L, 1L) ||
      !is.numeric(links) || is.object(links) || !is.null(dim(links)) || any(!is.finite(links)) ||
      any(links != floor(links)) || any(links < 1L | links > length(backend$links)) || anyDuplicated(links) > 0L)
    cp04_abort("Invalid explicit retained moment IF link batch.", "input")
  if (!identical(as.integer(cov_plan$blocks[[1L]]), seq_len(L))) return(NULL)
  if (identical(phase, "covariance") && !isTRUE(cov_plan$chunk_rows == chunk_rows)) return(NULL)
  if (chunk_rows > N || chunk_rows > control$max_chunk_rows)
    cp04_abort("Invalid retained moment IF row partition.", "input")
  ipw_owner_check_matrix_dims(N, L, label = "CP04 retained moment IF")
  boundary <- cp04_retained_if_boundary_cost(backend, links, chunk_rows)
  entries <- ipw_owner_safe_sum(c(cp04_retained_if_entries(N, L, chunk_rows), boundary$entries))
  checks <- ipw_owner_safe_sum(c(32, ipw_owner_safe_product(8, ceiling(N / chunk_rows) + 1)))
  cache_work <- ipw_owner_safe_sum(c(cp04_retained_if_work(N, L),
    ipw_owner_safe_product(checks, boundary$work)))
  cap <- min(control$max_retained_entries, control$max_matrix_entries)
  total <- ipw_owner_safe_sum(c(live, entries))
  work_total <- ipw_owner_safe_sum(c(work, cache_work))
  retained_probe <- ipw_owner_safe_sum(c(retained, entries))
  if (total > cap || retained_probe > cap || work_total > control$max_work_ops) return(NULL)
  probe <- function(expr) tryCatch(expr, cdid_ipw_owner_bound_error = cp04_capacity_decline)
  if (identical(phase, "covariance")) {
    cov_probe <- probe({
      ipw_owner_check_matrix_dims(N, L, backend$control$max_matrix_entries, "CP04 retained moment IF")
      LL <- ipw_owner_check_matrix_dims(L, L, backend$control$max_matrix_entries, "CP04 complete Gram")
      gram_entries <- ipw_owner_safe_product(4, LL, backend$control$max_matrix_entries, "CP04 complete Gram buffers")
      pl <- cp04_stream_plan(backend, links, retained_probe)
      if (length(pl$blocks) == 1L && identical(as.integer(pl$blocks[[1L]]), seq_len(L)) &&
          isTRUE(pl$chunk_rows == chunk_rows)) {
        cp04_budget(backend$control,
          c(retained_probe, backend$accounting$owned_entries, backend$accounting$borrowed_entries),
          c(pl$peak, gram_entries), pl$work, gram_entries, "CP04 complete Gram")
        pre_cov <- cp04_preflight_links(backend, links, retained_probe, chunk_rows = chunk_rows)
        ipw_owner_safe_sum(c(pre_cov$work_ops_bound, cache_work), backend$control$max_work_ops,
          "CP04 retained moment IF fill work")
      }
      pl
    })
    if (is.null(cov_probe) || length(cov_probe$blocks) != 1L ||
        !identical(as.integer(cov_probe$blocks[[1L]]), seq_len(L)) ||
        !isTRUE(cov_probe$chunk_rows == chunk_rows)) return(NULL)
  } else {
    gmm_probe <- probe({
      ipw_owner_check_matrix_dims(N, L, backend$control$max_matrix_entries, "CP04 retained moment IF")
      pre_gmm <- cp04_preflight_links(backend, links, retained_probe)
      if (isTRUE(pre_gmm$chunk_rows == chunk_rows)) {
        ipw_owner_safe_sum(c(retained_probe, pre_gmm$peak_concurrent_entries),
          min(backend$control$max_retained_entries, backend$control$max_matrix_entries),
          "CP04 GMM producer and retained moment IF simultaneous entries")
        ipw_owner_safe_sum(c(pre_gmm$work_ops_bound, cache_work), backend$control$max_work_ops,
          "CP04 retained moment IF fill work")
      }
      pre_gmm
    })
    if (is.null(gmm_probe) || !isTRUE(gmm_probe$chunk_rows == chunk_rows)) return(NULL)
  }
  store <- new.env(parent = emptyenv())
  store$kind <- "cp04_retained_if"; store$state <- "reserved"
  store$phase <- phase; store$max_slice_passes <- 1L
  store$backend_snapshot <- backend
  store$backend_ref <- backend$ref; store$reader_ref <- backend$reader_ref
  store$frame_ref <- dom$frame_ref; store$panel_ref <- dom$panel_ref
  store$reader_snapshot <- reader; store$state_snapshot <- st; store$frame_snapshot <- frame
  store$panel_snapshot <- dom$panel
  store$reference_ref <- backend$reference_ref
  store$reference_snapshot <- if (is.null(backend$reference_ref)) NULL else chained_backend_reference(backend)
  store$selection_ref <- if (is.null(backend$selection_backend)) NULL else backend$selection_backend$ref
  store$selection_snapshot <- backend$selection_backend
  store$selection_reader_snapshot <- if (is.null(backend$selection_backend)) NULL else
    get("cp04_reader", backend$selection_backend$reader_ref, inherits = FALSE)
  store$design_ref <- backend$design_ref
  store$design_snapshot <- if (is.null(backend$design_ref)) NULL else get("imputation_design_owner", backend$design_ref, inherits = FALSE)
  store$layout_snapshot <- backend$layout
  store$layout_state_snapshot <- if (is.null(backend$layout)) NULL else backend$layout$state
  store$X_F_snapshot <- backend$X_F
  store$pair_outcomes_snapshot <- backend$pair_outcomes
  for (nm in c("reader_snapshot", "state_snapshot", "frame_snapshot", "panel_snapshot",
               "reference_snapshot", "selection_snapshot", "design_snapshot",
               "layout_snapshot", "layout_state_snapshot", "X_F_snapshot", "pair_outcomes_snapshot"))
    lockBinding(nm, store)
  store$id <- dom$id; store$N <- N; store$n_A <- dom$n_A; store$lambda <- dom$lambda
  store$links <- as.integer(links); store$L <- L; store$B <- chunk_rows; store$scale <- "core"
  store$entries <- entries; store$work <- cache_work; store$admission_total <- total; store$cap <- cap
  store$producer_passes <- 0L; store$slice_passes <- 0L; store$batches <- 0L
  store$metadata_entries <- boundary$entries; store$boundary_work <- boundary$work
  for (nm in setdiff(ls(store, all.names = TRUE), c("state", "producer_passes", "slice_passes", "batches")))
    if (!bindingIsLocked(nm, store)) lockBinding(nm, store)
  store
}

cp04_retained_if_phase_record <- function(cov_record = NULL, gmm_record = NULL) {
  if (is.null(cov_record)) cov_record <- list(status = "streamed")
  if (is.null(gmm_record)) gmm_record <- list(status = "streamed")
  cov_ok <- !identical(cov_record$status, "streamed")
  gmm_ok <- !identical(gmm_record$status, "streamed")
  if (!cov_ok && !gmm_ok) return(list(status = "streamed"))
  active <- if (cov_ok) cov_record else gmm_record
  status <- if ((!cov_ok || identical(cov_record$status, "retained")) &&
      (!gmm_ok || identical(gmm_record$status, "retained"))) "retained" else
    paste0("incomplete:", if (cov_ok && !identical(cov_record$status, "retained")) cov_record$status else gmm_record$status)
  cov_B <- if (cov_ok) as.double(cov_record$chunk_rows) else NA_real_
  gmm_B <- if (gmm_ok) as.double(gmm_record$chunk_rows) else NA_real_
  cov_batches <- if (cov_ok) as.integer(cov_record$n_batches) else NA_integer_
  gmm_batches <- if (gmm_ok) as.integer(gmm_record$n_batches) else NA_integer_
  prod_passes <- as.integer((if (cov_ok) cov_record$producer_passes else 0L) +
    (if (gmm_ok) gmm_record$producer_passes else 0L))
  slice_passes <- as.integer((if (cov_ok) cov_record$slice_passes else 0L) +
    (if (gmm_ok) gmm_record$slice_passes else 0L))
  streamed_passes <- as.integer((if (cov_ok) 0L else 2L) + (if (gmm_ok) 0L else 2L))
  lifetime <- if (cov_ok && gmm_ok) {
    paste("two independent phase-local stores (never simultaneously live): covariance store reserved before allocation,",
      "filled by the covariance mean pass, sliced for the covariance Gram pass at the original covariance partition,",
      "and released before allocating GMM consumers; GMM store reserved before allocation, filled by the first GMM",
      "projection pass, sliced for the second GMM projection pass at the original GMM partition, and released before",
      "target IF transport, auxiliary reporting and final diagnostic copying")
  } else if (cov_ok) {
    paste("covariance-only phase-local store: reserved before allocation, filled by the covariance mean pass,",
      "sliced for the covariance Gram pass at the original covariance partition, and released before allocating GMM",
      "consumers; GMM phase streamed at its original partition")
  } else {
    paste("GMM-only phase-local store: covariance phase streamed at its original partition; GMM store reserved before",
      "allocation, filled by the first GMM projection pass, sliced for the second GMM projection pass at the original",
      "GMM partition, and released before target IF transport, auxiliary reporting and final diagnostic copying")
  }
  list(
    status = status,
    mode = "phase_local",
    scale = active$scale,
    n_rows = active$n_rows,
    n_links = active$n_links,
    chunk_rows = c(covariance = cov_B, gmm = gmm_B),
    covariance_chunk_rows = cov_B,
    gmm_chunk_rows = gmm_B,
    n_batches = c(covariance = cov_batches, gmm = gmm_batches),
    producer_passes = prod_passes,
    slice_passes = slice_passes,
    streamed_passes = streamed_passes,
    simultaneous_stores = 1L,
    reserved_entries = max(c(if (cov_ok) cov_record$reserved_entries else 0,
      if (gmm_ok) gmm_record$reserved_entries else 0)),
    admission_total_entries = max(c(if (cov_ok) cov_record$admission_total_entries else 0,
      if (gmm_ok) gmm_record$admission_total_entries else 0)),
    cap_entries = min(c(if (cov_ok) cov_record$cap_entries else Inf,
      if (gmm_ok) gmm_record$cap_entries else Inf)),
    snapshot_metadata_entries = max(c(if (cov_ok) cov_record$snapshot_metadata_entries else 0,
      if (gmm_ok) gmm_record$snapshot_metadata_entries else 0)),
    boundary_check_work_ops = ipw_owner_safe_sum(c(
      if (cov_ok) cov_record$boundary_check_work_ops else 0,
      if (gmm_ok) gmm_record$boundary_check_work_ops else 0)),
    work_ops_bound = ipw_owner_safe_sum(c(
      if (cov_ok) cov_record$work_ops_bound else 0,
      if (gmm_ok) gmm_record$work_ops_bound else 0)),
    covariance = cov_record,
    gmm = gmm_record,
    auxiliary = list(status = "not_extended"),
    lifetime = lifetime
  )
}

cp04_retained_if_bind <- function(backend, store) {
  if (!is.environment(store) || !identical(store$kind, "cp04_retained_if") || !is.list(backend) ||
      !is.environment(backend$ref) || !identical(store$backend_ref, backend$ref))
    cp04_abort("Retained moment IF storage can only be bound to the backend it was reserved for.", "integrity")
  structure(list(backend = backend, store = store), class = "cdid_cp04_retained_if_bound")
}

cp04_retained_if_unbind <- function(backend) {
  if (inherits(backend, "cdid_cp04_retained_if_bound")) list(backend = backend$backend, store = backend$store) else
    list(backend = backend, store = NULL)
}

# Every emission boundary: the original validators (seal, reader, frame owner,
# route owners) and the reserved snapshot of the reader's roster, scale
# factors, frame/panel/reference/selection/design/layout references, link block
# and row partition.
cp04_retained_if_check <- function(store, backend, links, B, scale) {
  if (!is.environment(store) || !identical(store$kind, "cp04_retained_if") || identical(store$state, "released"))
    cp04_abort("Retained moment IF storage is not available.", "integrity")
  cp04_retained_if_validate_backend(backend)
  dom <- cp04_backend_domain(backend)
  reader <- dom$reader; st <- dom$state
  frame <- dom$frame; panel <- dom$panel
  ids <- dom$id
  reader_ok <- if (dom$is_raw) {
    is.null(reader) && is.null(store$reader_snapshot) && is.null(st) && is.null(store$state_snapshot) &&
      is.null(backend$reader_ref) && is.null(store$reader_ref)
  } else {
    gmm_sexp_aliases(reader, store$reader_snapshot) && gmm_sexp_aliases(st, store$state_snapshot) &&
      identical(backend$reader_ref, store$reader_ref)
  }
  if (!gmm_sexp_aliases(backend, store$backend_snapshot) || !reader_ok ||
      !gmm_sexp_aliases(frame, store$frame_snapshot) || !gmm_sexp_aliases(panel, store$panel_snapshot) ||
      !identical(backend$ref, store$backend_ref) ||
      !identical(dom$frame_ref, store$frame_ref) || !identical(dom$panel_ref, store$panel_ref) ||
      !identical(ids, store$id) || !identical(dom$N, store$N) || !identical(dom$n_A, store$n_A) ||
      !identical(dom$lambda, store$lambda))
    cp04_abort("Retained moment IF storage does not describe the backend's sealed reader, frame and panel.", "integrity")
  if (!identical(backend$reference_ref, store$reference_ref) ||
      (!is.null(store$reference_ref) && !gmm_sexp_aliases(chained_backend_reference(backend), store$reference_snapshot)))
    cp04_abort("Retained moment IF storage does not describe the backend's sealed reference owner.", "integrity")
  if (!identical(if (is.null(backend$selection_backend)) NULL else backend$selection_backend$ref, store$selection_ref) ||
      (!is.null(store$selection_ref) && (!gmm_sexp_aliases(backend$selection_backend, store$selection_snapshot) ||
        !gmm_sexp_aliases(get("cp04_reader", backend$selection_backend$reader_ref, inherits = FALSE),
          store$selection_reader_snapshot))))
    cp04_abort("Retained moment IF storage does not describe the backend's sealed selection owner.", "integrity")
  if (!identical(backend$design_ref, store$design_ref) ||
      (!is.null(store$design_ref) && !gmm_sexp_aliases(get("imputation_design_owner", backend$design_ref, inherits = FALSE), store$design_snapshot)))
    cp04_abort("Retained moment IF storage does not describe the backend's sealed imputation design owner.", "integrity")
  # The layout state is optional: pure imputation layouts carry it; the dated
  # augmented layout leaves it to its sealed layout reader (checked above). The
  # native alias test is FALSE for NULL, so presence and identity are separate.
  if (!identical(is.null(backend$layout), is.null(store$layout_snapshot)) ||
      (!is.null(store$layout_snapshot) && (!gmm_sexp_aliases(backend$layout, store$layout_snapshot) ||
        !identical(is.null(backend$layout$state), is.null(store$layout_state_snapshot)) ||
        (!is.null(store$layout_state_snapshot) &&
          !gmm_sexp_aliases(backend$layout$state, store$layout_state_snapshot)))))
    cp04_abort("Retained moment IF storage does not describe the backend's sealed layout owner.", "integrity")
  if (!identical(is.null(backend$X_F), is.null(store$X_F_snapshot)) ||
      (!is.null(store$X_F_snapshot) && !gmm_sexp_aliases(backend$X_F, store$X_F_snapshot)) ||
      !identical(is.null(backend$pair_outcomes), is.null(store$pair_outcomes_snapshot)) ||
      (!is.null(store$pair_outcomes_snapshot) && !gmm_sexp_aliases(backend$pair_outcomes, store$pair_outcomes_snapshot)))
    cp04_abort("Retained moment IF storage does not describe the backend's sealed custom covariates or outcomes.", "integrity")
  if (!identical(as.integer(links), store$links) || !identical(scale, store$scale) || !isTRUE(B == store$B))
    cp04_abort("Retained moment IF storage serves one link block, one scale and one row partition.", "integrity")
  list(id = ids, N = dom$N, n_A = dom$n_A, lambda = dom$lambda)
}

# One producer pass (the original plan and chunk loop; each chunk is stored raw
# and handed on unchanged), then slice passes over the same row partition. A
# locked raw-matrix binding retains the exact completed object; copy-on-modify
# value, label or order changes cannot alias it. IDs come from the live reader,
# checked against the compact source snapshots. The caller's preflight admitted the original replay
# work; the fill/slice work is charged in addition.
cp04_retained_if_replay <- function(store, backend, links, callback, scale, B, pre) {
  completed <- FALSE
  on.exit(if (!completed && is.environment(store) && identical(store$kind, "cp04_retained_if"))
    cp04_retained_if_release(store), add = TRUE)
  st <- cp04_retained_if_check(store, backend, links, B, scale)
  ctrl <- backend$control; N <- store$N; L <- store$L
  starts <- seq.int(1L, N, by = B)
  if (identical(store$state, "reserved")) {
    ipw_owner_safe_sum(c(pre$work_ops_bound, store$work), ctrl$max_work_ops, "CP04 retained moment IF fill work")
    ipw_owner_check_matrix_dims(N, L, ctrl$max_matrix_entries, "CP04 retained moment IF")
    plan <- chained_full_if_plan(backend, links, B)
    store$state <- "filling"
    phi <- matrix(0, N, L); filled <- 0; cn <- NULL
    for (s in starts) {
      rows <- seq.int(s, min(N, s + B - 1L))
      b <- chained_full_if_chunk(backend, rows, links, scale, plan)
      if (!is.list(b) || !identical(b$rows, rows) || !identical(as.integer(b$links), store$links) ||
          !identical(b$scale, scale) || !is.matrix(b$phi) || !is.double(b$phi) ||
          nrow(b$phi) != length(rows) || ncol(b$phi) != L || any(!is.finite(b$phi)) ||
          !identical(b$id, st$id[rows]) || !identical(b$normalization_count, st$n_A))
        cp04_abort("The producer chunk does not match the reserved retained moment IF layout.", "integrity")
      if (filled == 0) cn <- colnames(b$phi) else if (!identical(colnames(b$phi), cn))
        cp04_abort("The producer chunk link labels changed within one pass.", "integrity")
      phi[rows, ] <- b$phi
      filled <- filled + length(rows)
      callback(b)
      st <- cp04_retained_if_check(store, backend, links, B, scale)
      if (!identical(store$state, "filling") || exists("phi", store, inherits = FALSE) ||
          exists("raw_ref", store, inherits = FALSE))
        cp04_abort("Retained moment IF storage was modified during the producer pass.", "integrity")
      # A phase fill must not retain the previous source output while the next
      # producer allocates its output and scratch. Keep the accepted common-
      # partition path unchanged.
      if (exists("phase", store, inherits = FALSE)) rm(b)
    }
    if (filled != N) cp04_abort("The producer pass did not cover every roster row exactly once.", "integrity")
    colnames(phi) <- cn
    raw_ref <- new.env(parent = emptyenv()); raw_ref$phi <- phi
    lockEnvironment(raw_ref, bindings = TRUE)
    store$phi <- phi; store$raw_ref <- raw_ref; lockBinding("raw_ref", store)
    store$state <- "filled"; store$producer_passes <- store$producer_passes + 1L; store$batches <- length(starts)
    completed <- TRUE
    return(invisible(NULL))
  }
  if (!identical(store$state, "filled")) cp04_abort("Retained moment IF storage was not completely filled.", "integrity")
  if (exists("max_slice_passes", store, inherits = FALSE) && !is.null(store$max_slice_passes) &&
      store$slice_passes >= store$max_slice_passes)
    cp04_abort("Retained moment IF storage exceeded its admitted slice passes.", "integrity")
  ipw_owner_safe_sum(c(pre$work_ops_bound, store$work), ctrl$max_work_ops, "CP04 retained moment IF slice work")
  check_slice <- function() {
    st <- cp04_retained_if_check(store, backend, links, B, scale)
    if (!identical(store$state, "filled"))
      cp04_abort("Retained moment IF storage was not completely filled.", "integrity")
    full <- store$phi
    if (!exists("raw_ref", store, inherits = FALSE) || !bindingIsLocked("raw_ref", store) ||
        !is.environment(store$raw_ref) || !identical(parent.env(store$raw_ref), emptyenv()) ||
        !environmentIsLocked(store$raw_ref) || !exists("phi", store$raw_ref, inherits = FALSE) ||
        !bindingIsLocked("phi", store$raw_ref) || !gmm_sexp_aliases(full, get("phi", store$raw_ref, inherits = FALSE)) ||
        !is.matrix(full) || !is.double(full) || nrow(full) != N || ncol(full) != L)
      cp04_abort("Retained moment IF storage lost its shape after the producer pass.", "integrity")
    st
  }
  st <- check_slice()
  for (s in starts) {
    rows <- seq.int(s, min(N, s + B - 1L))
    phi <- store$phi[rows, , drop = FALSE]
    callback(list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
      normalization_count = st$n_A))
    st <- check_slice()
  }
  store$slice_passes <- store$slice_passes + 1L
  completed <- TRUE
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# CP08 retained moment IF, auxiliary consumer. An auxiliary descriptor is
# appended to the same sealed backend at index M + j with its own label; its
# complete link IF is nevertheless the SAME evaluation recipe as a catalog row
# whenever the backend evaluates both from identical canonical descriptors:
# OR the regression pool, the oriented pair, the target cohort and the support
# (rates); DR additionally the treatment model and its ordered control
# supports. Within one sealed backend those descriptors fix the global beta,
# the local intercept, the propensity fit and every rate the column reads, so
# an exact descriptor match (never a label, a numeric index membership or a
# value closeness) identifies the retained column. With F = A (N = n_A) and
# lambda = 1 the stored scale = "core" value `lambda * phi` is the frame-scale
# value itself and the frame normalization count n_A is N, so the matched
# column is served unchanged into the original output slot with the existing
# reporting sign. Unmatched columns replay through the ordinary producer in
# their original order. The extension is admitted BEFORE the store is
# allocated, against the complete simultaneous envelope through the auxiliary
# consumer (transported target IF, auxiliary outputs, diagnostic copies) and
# both original auxiliary planners re-run with the store charged; declined,
# the store keeps its original lifetime (released after the GMM solve).
# ---------------------------------------------------------------------------

# Extra charged work for the one auxiliary slice pass (the store's 128 N L
# covers the fill and three slice passes at 32 N L each).
cp04_retained_if_auxiliary_work <- function(N, L) ipw_owner_safe_product(32, ipw_owner_safe_product(N, L))

# The canonical evaluation recipe of backend link `j`, or NULL when the record
# is not a plain complete-support recipe (a dated/reference recipe, a missing
# or non-integral descriptor). Whole-number descriptors are compared exactly as
# integers; nothing is rounded or matched by proximity.
cp04_retained_if_recipe <- function(backend, j) {
  if (!is.list(backend) || !is.list(backend$links) || !is.numeric(j) || length(j) != 1L || is.na(j) ||
      j < 1 || j > length(backend$links)) return(NULL)
  l <- backend$links[[j]]
  if (!is.list(l) || !is.null(l$reference)) return(NULL)
  whole <- function(x) is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) >= 1L && all(is.finite(x)) &&
    all(x == floor(x)) && all(abs(x) <= .Machine$integer.max)
  fields <- c("regression", "pair_id", "cohort", "support_id")
  dr <- inherits(backend, c("cdid_dr_baseline", "cdid_ipw_raw_backend"))
  if (dr) fields <- c(fields, "model")
  rec <- vector("list", length(fields)); names(rec) <- fields
  for (f in fields) {
    x <- l[[f]]
    if (!whole(x) || length(x) != 1L) return(NULL)
    rec[[f]] <- as.integer(x)
  }
  if (dr) {
    cs <- l$control_support_ids
    if (!whole(cs)) return(NULL)
    rec$control_support_ids <- as.integer(cs)
  }
  rec
}

# For every requested auxiliary backend index the position of the retained
# store column with the identical canonical recipe, or NA. A requested index
# that is itself a retained link is never served (membership is not a recipe);
# a recipe retained more than once is never served either.
cp04_retained_if_auxiliary_map <- function(store, backend, backend_ids) {
  col <- rep(NA_integer_, length(backend_ids))
  if (!is.environment(store) || !identical(store$kind, "cp04_retained_if") || !is.list(backend) ||
      !is.numeric(backend_ids) || is.object(backend_ids) || !is.null(dim(backend_ids)) || !length(backend_ids) ||
      anyNA(backend_ids) || any(backend_ids != floor(backend_ids)) || any(backend_ids < 1) ||
      any(backend_ids > length(backend$links))) return(col)
  retained <- lapply(store$links, function(k) cp04_retained_if_recipe(backend, k))
  for (q in seq_along(backend_ids)) {
    id <- as.integer(backend_ids[q])
    if (id %in% store$links) next
    rq <- cp04_retained_if_recipe(backend, id)
    if (is.null(rq)) next
    hit <- which(vapply(retained, function(r) !is.null(r) && identical(r, rq), logical(1)))
    if (length(hit) == 1L) col[q] <- hit
  }
  col
}

# Admit the auxiliary extension of a RESERVED (unallocated) store or leave it
# unextended (NULL). `backend_ids` are the requested reporting columns' backend
# indices in their reporting order; `retained` is the auxiliary replay's
# retained argument (pipeline envelope plus the reporting output envelope);
# `live` is the complete simultaneous envelope through the auxiliary consumer
# INCLUDING the store's entries; `work` the cumulative work before this charge.
# Structural facts (no reference inputs, N = n_A, lambda = 1, at least one
# exact recipe match), then the arithmetic admission, then the original
# auxiliary planner probes: the whole-request preflight and its simultaneous
# check as the consumer runs them, and for the unmatched remainder the planner
# with and without the store charged, which must issue the same row partition.
# Only typed bound refusals of the probes decline; integrity errors propagate.
cp04_retained_if_auxiliary_admit <- function(store, backend, backend_ids, retained, live, work, control) {
  if (!is.environment(store) || !identical(store$kind, "cp04_retained_if") || !identical(store$state, "reserved") ||
      (exists("phase", store, inherits = FALSE) && !is.null(store$phase)) ||
      !is.null(store$auxiliary)) return(NULL)
  if (!is.list(backend) || !inherits(backend, c("cdid_or_baseline", "cdid_dr_baseline", "cdid_ipw_raw_backend")) ||
      !identical(store$backend_ref, backend$ref) || !is.null(backend$reference_ref) ||
      is.null(backend$complete_support)) return(NULL)
  if (!is.numeric(backend_ids) || is.object(backend_ids) || !is.null(dim(backend_ids)) || !length(backend_ids) ||
      anyNA(backend_ids) || any(backend_ids != floor(backend_ids)) || any(backend_ids < 1) ||
      any(backend_ids > length(backend$links)) || any(backend_ids %in% store$links)) return(NULL)
  if (!isTRUE(store$lambda == 1) || !identical(as.double(store$N), as.double(store$n_A))) return(NULL)
  control <- validate_ipw_owner_control(control)
  for (nm in c("retained", "live", "work")) ipw_owner_int_bound(get(nm), nm, min_val = 0)
  col <- cp04_retained_if_auxiliary_map(store, backend, backend_ids)
  if (!any(!is.na(col))) return(NULL)
  N <- store$N; L <- store$L
  aux_work <- ipw_owner_safe_sum(c(cp04_retained_if_auxiliary_work(N, L),
    ipw_owner_safe_product(N + 4, store$boundary_work)))
  work_total <- ipw_owner_safe_sum(c(work, aux_work))
  if (live < store$entries) cp04_abort("Auxiliary lifetime omits the retained moment IF charge.", "integrity")
  if (live > store$cap || work_total > control$max_work_ops) return(NULL)
  probe <- function(expr) tryCatch(expr, cdid_ipw_owner_bound_error = cp04_capacity_decline)
  ids <- unique(as.integer(backend_ids))
  whole <- probe({
    pre <- cp04_preflight_links(backend, ids, retained)
    ipw_owner_safe_sum(c(retained, pre$peak_concurrent_entries), control$max_matrix_entries,
      "CP04 auxiliary producer and requested reporting outputs")
    pre
  })
  if (is.null(whole)) return(NULL)
  ids_stream <- unique(as.integer(backend_ids[is.na(col)]))
  chunk_rows_stream <- NA_real_
  if (length(ids_stream)) {
    without <- probe(cp04_preflight_links(backend, ids_stream, retained, chunk_rows = whole$chunk_rows))
    with_store <- probe({
      pre <- cp04_preflight_links(backend, ids_stream, ipw_owner_safe_sum(c(retained, store$entries)),
        chunk_rows = whole$chunk_rows)
      ipw_owner_safe_sum(c(retained, store$entries, pre$peak_concurrent_entries), control$max_matrix_entries,
        "CP04 auxiliary producer, retained moment IF and requested reporting outputs")
      pre
    })
    if (is.null(without) || is.null(with_store) || !isTRUE(without$chunk_rows == whole$chunk_rows) ||
        !isTRUE(with_store$chunk_rows == whole$chunk_rows)) return(NULL)
    chunk_rows_stream <- as.double(with_store$chunk_rows)
  } else if (is.null(probe(ipw_owner_safe_sum(c(retained, store$entries), control$max_matrix_entries,
    "CP04 retained moment IF and requested reporting outputs")))) return(NULL)
  store$auxiliary <- list(backend_ids = as.integer(backend_ids), ids = ids, column = col, ids_stream = ids_stream,
    chunk_rows = as.double(whole$chunk_rows), chunk_rows_stream = chunk_rows_stream,
    retained = as.double(retained), live_total = live, work = aux_work,
    columns_requested = length(backend_ids), columns_served = sum(!is.na(col)), columns_streamed = sum(is.na(col)),
    slice_passes = 0L)
  invisible(store)
}

# The consumer's view of the admitted extension: NULL when the store carries
# none (the ordinary replay runs), otherwise the admitted record after the
# actual request (indices, order, retained argument, recomputed recipe map) is
# found identical to the admitted one; any difference is an integrity failure
# of the sealed state, never a silent re-admission.
cp04_retained_if_auxiliary_service <- function(store, backend, backend_ids, retained) {
  if (is.null(store)) return(NULL)
  if (!is.environment(store) || !identical(store$kind, "cp04_retained_if") || identical(store$state, "released"))
    cp04_abort("Retained moment IF storage is not available.", "integrity")
  if (exists("phase", store, inherits = FALSE) && !is.null(store$phase))
    cp04_abort("Phase-local retained moment IF storage does not serve auxiliary reporting.", "integrity")
  aux <- store$auxiliary
  if (is.null(aux)) return(NULL)
  cp04_retained_if_check(store, backend, store$links, store$B, "core")
  if (!identical(store$state, "filled")) cp04_abort("Retained moment IF storage was not completely filled.", "integrity")
  if (!identical(aux$backend_ids, as.integer(backend_ids)) || !identical(aux$retained, as.double(retained)) ||
      !identical(aux$column, cp04_retained_if_auxiliary_map(store, backend, backend_ids)))
    cp04_abort("The auxiliary reporting request differs from the admitted retained moment IF extension.", "integrity")
  aux
}

# Serve the matched columns: one slice pass over the store's own row partition,
# the frame-scale batch schema (`rows`, live `id`, served backend `links`,
# `phi`, scale = "frame", normalization count N) with the raw identity, the
# sealed source snapshots and the F = A, lambda = 1 identity re-checked at every
# boundary from the live reader state.
cp04_retained_if_auxiliary_slices <- function(store, backend, backend_ids, callback) {
  completed <- FALSE
  on.exit(if (!completed && is.environment(store) && identical(store$kind, "cp04_retained_if"))
    cp04_retained_if_release(store), add = TRUE)
  if (!is.function(callback)) cp04_abort("Auxiliary slice callback must be a function.", "input")
  st <- cp04_retained_if_check(store, backend, store$links, store$B, "core")
  if (exists("phase", store, inherits = FALSE) && !is.null(store$phase))
    cp04_abort("Phase-local retained moment IF storage does not serve auxiliary reporting.", "integrity")
  aux <- store$auxiliary
  if (is.null(aux) || !identical(aux$backend_ids, as.integer(backend_ids)))
    cp04_abort("Retained moment IF storage carries no admitted auxiliary extension for this request.", "integrity")
  if (!identical(store$state, "filled")) cp04_abort("Retained moment IF storage was not completely filled.", "integrity")
  ctrl <- backend$control; N <- store$N; L <- store$L; B <- aux$chunk_rows
  ipw_owner_safe_sum(aux$work, ctrl$max_work_ops, "CP04 retained moment IF auxiliary slice work")
  served <- which(!is.na(aux$column)); cols <- aux$column[served]
  links_out <- as.integer(backend_ids[served])
  labels <- vapply(backend$links[links_out], `[[`, character(1), "link_id")
  # The auxiliary consumer keeps its original whole-request partition, which
  # need not equal the main producer's stored partition. Its original output
  # envelope covers these requested-column slices at this issued B.
  check <- function() {
    st <- cp04_retained_if_check(store, backend, store$links, store$B, "core")
    if (!identical(store$state, "filled")) cp04_abort("Retained moment IF storage was not completely filled.", "integrity")
    if (!isTRUE(st$lambda == 1) || !identical(as.double(st$N), as.double(st$n_A)) ||
        !isTRUE(store$lambda == 1) || !identical(as.double(store$N), as.double(store$n_A)))
      cp04_abort("Retained moment IF storage serves auxiliary columns only on the F = A, lambda = 1 frame.", "integrity")
    full <- store$phi
    if (!exists("raw_ref", store, inherits = FALSE) || !bindingIsLocked("raw_ref", store) ||
        !is.environment(store$raw_ref) || !identical(parent.env(store$raw_ref), emptyenv()) ||
        !environmentIsLocked(store$raw_ref) || !exists("phi", store$raw_ref, inherits = FALSE) ||
        !bindingIsLocked("phi", store$raw_ref) || !gmm_sexp_aliases(full, get("phi", store$raw_ref, inherits = FALSE)) ||
        !is.matrix(full) || !is.double(full) || nrow(full) != N || ncol(full) != L)
      cp04_abort("Retained moment IF storage lost its shape after the producer pass.", "integrity")
    st
  }
  st <- check()
  for (s in seq.int(1L, N, by = B)) {
    rows <- seq.int(s, min(N, s + B - 1L))
    phi <- store$phi[rows, cols, drop = FALSE]
    dimnames(phi) <- list(NULL, labels)
    callback(list(rows = rows, id = st$id[rows], links = links_out, phi = phi, scale = "frame",
      normalization_count = st$N))
    # This is also the next emission's boundary check, with a final check after
    # the last callback rather than leaving its source/raw mutation unobserved.
    st <- check()
  }
  store$auxiliary$slice_passes <- store$auxiliary$slice_passes + 1L
  completed <- TRUE
  invisible(NULL)
}

# Release after the last consumer (the GMM solve, or the admitted auxiliary
# reporting): the matrix and borrowed snapshots are dropped; the scalar record
# of what happened is returned (NULL store: the ordinary stream served every
# pass). `slice_passes` counts the covariance Gram and GMM passes only; the
# auxiliary service is recorded apart.
cp04_retained_if_release <- function(store) {
  if (is.null(store)) return(list(status = "streamed"))
  aux <- store$auxiliary
  phase <- if (exists("phase", store, inherits = FALSE)) store$phase else NULL
  lifetime <- if (!is.null(phase) && identical(phase, "covariance")) {
    paste("reserved before allocation inside the covariance phase simultaneous envelope; filled by the covariance",
      "mean pass; sliced for the covariance Gram pass at the original covariance partition; released before",
      "allocating GMM consumers")
  } else if (!is.null(phase)) {
    paste("reserved before allocation inside the GMM phase simultaneous envelope (after covariance store release);",
      "filled by the first GMM projection pass; sliced for the second GMM projection pass at the original GMM",
      "partition; released before target IF transport, auxiliary reporting and final diagnostic copying")
  } else if (is.null(aux)) {
    paste("reserved before allocation inside the complete simultaneous envelope; filled by the covariance",
      "mean pass; sliced for the covariance Gram pass and both GMM projection passes; released before",
      "auxiliary replay and final output")
  } else {
    paste("reserved before allocation inside the complete simultaneous",
      "envelope through the auxiliary consumer; filled by the covariance mean pass; sliced for the covariance Gram",
      "pass and both GMM projection passes; exact same-recipe auxiliary columns served from the store in frame",
      "scale (F = A, lambda = 1); released after auxiliary reporting and before final output")
  }
  rec <- c(
    list(status = if (identical(store$state, "filled")) "retained" else paste0("incomplete:", store$state)),
    if (is.null(phase)) list() else list(phase = phase),
    list(scale = store$scale, n_rows = store$N, n_links = store$L, chunk_rows = store$B, n_batches = store$batches,
      producer_passes = store$producer_passes, slice_passes = store$slice_passes,
      reserved_entries = store$entries, admission_total_entries = store$admission_total, cap_entries = store$cap,
      snapshot_metadata_entries = store$metadata_entries, boundary_check_work_ops = store$boundary_work,
      work_ops_bound = store$work,
      auxiliary = if (is.null(aux)) list(status = "not_extended") else list(status = "extended",
        columns_requested = aux$columns_requested, columns_served = aux$columns_served,
        columns_streamed = aux$columns_streamed, slice_passes = aux$slice_passes,
        chunk_rows = aux$chunk_rows, stream_chunk_rows = aux$chunk_rows_stream,
        live_total_entries = aux$live_total, work_ops_bound = aux$work),
      lifetime = lifetime)
  )
  drop <- intersect(c("phi", "raw_ref", "id", "backend_snapshot", "backend_ref", "reader_ref", "frame_ref", "panel_ref",
    "reference_ref", "selection_ref", "design_ref", "reader_snapshot", "state_snapshot", "frame_snapshot", "panel_snapshot",
    "reference_snapshot", "selection_snapshot", "selection_reader_snapshot", "design_snapshot", "layout_snapshot", "layout_state_snapshot",
    "X_F_snapshot", "pair_outcomes_snapshot", "auxiliary"), ls(store, all.names = TRUE))
  for (nm in drop) if (bindingIsLocked(nm, store)) unlockBinding(nm, store)
  rm(list = drop, envir = store)
  store$state <- "released"
  rec
}

cp04_replay_links <- function(backend, links, callback, scale = "core", chunk_rows = NULL) {
  bound <- cp04_retained_if_unbind(backend); backend <- bound$backend
  completed <- FALSE
  on.exit(if (!completed && is.environment(bound$store) && identical(bound$store$kind, "cp04_retained_if"))
    cp04_retained_if_release(bound$store), add = TRUE)
  pre <- cp04_preflight_links(backend, links, chunk_rows = chunk_rows)
  dom <- cp04_backend_domain(backend)
  N <- dom$N; B <- pre$chunk_rows
  # CP08 retained moment IF: a bound store fills on its first pass and slices afterwards.
  if (!is.null(bound$store)) {
    cp04_retained_if_replay(bound$store, backend, links, callback, scale, B, pre)
    completed <- TRUE
    return(invisible(NULL))
  }
  # CP08 stage 4: one compact native plan per eligible complete-support replay.
  plan <- chained_full_if_plan(backend, links, B)
  for (s in seq.int(1L, N, by = B)) {
    callback(chained_full_if_chunk(backend, seq.int(s, min(N, s + B - 1L)), links, scale, plan))
    if (dom$is_raw) validate_ipw_raw_backend(backend)
  }
  completed <- TRUE
  invisible(NULL)
}

# Complete Gram in bounded row/link blocks. No N x moments array, masks per
# model, inverses, or dense whole-outcome cache is retained.
cp04_stream_plan_legacy <- function(backend, links, retained = 0) {
  st <- cp04_backend_domain(backend); ctrl <- backend$control
  L <- length(links); LB <- min(L, ctrl$max_chunk_rows)
  blocks <- split(seq_len(L), ceiling(seq_len(L) / LB))
  B <- min(st$N, ctrl$max_chunk_rows)
  for (a in seq_along(blocks)) for (b in a:length(blocks)) {
    ids <- unique(c(links[blocks[[a]]], links[blocks[[b]]]))
    B <- min(B, cp04_plan_rows(backend, list(ids), retained + 4 * L * L))
  }
  work <- 0; peak <- 0
  # Complete preflight before replay/allocations, including both Gram buffers.
  for (a in seq_along(blocks)) for (b in a:length(blocks)) {
    ids <- unique(c(links[blocks[[a]]], links[blocks[[b]]]))
    pre <- chained_full_if_budget(backend, B, ids)
    work <- ipw_owner_safe_sum(c(work, 2 * ceiling(st$N / B) * pre$work + 4 * st$N * length(ids)^2),
      ctrl$max_work_ops, "CP04 covariance cumulative work")
    peak <- max(peak, pre$output + pre$scratch)
  }
  list(work = work, peak = peak, blocks = blocks, L = L, chunk_rows = B)
}

# Prospective single full-width covariance block for supported sealed CP04
# backends at the exact legacy row partition, under the actual producer/output
# budgets and simultaneous Gram/scratch/backend/retained caps.
cp04_stream_plan_full_block <- function(backend, links, retained = 0, legacy = NULL) {
  if (!is.list(backend) || !inherits(backend, cp04_supported_backend_classes)) return(NULL)
  abort <- if (inherits(backend, "cdid_ipw_raw_backend")) ipw_raw_backend_abort else
    if (inherits(backend, "cdid_dr_baseline")) dr_baseline_abort else
    if (inherits(backend, "cdid_or_baseline")) or_baseline_abort else
    if (inherits(backend, "cdid_or_corollary1")) or_corollary1_abort else
    if (inherits(backend, "cdid_dr_corollary1")) dr_corollary1_abort else
    if (inherits(backend, "cdid_imputation_backend")) imputation_abort else
    if (inherits(backend, "cdid_augmented_backend")) augmented_abort else
    if (inherits(backend, "cdid_reference_ipw_backend")) reference_ipw_abort else cp04_abort
  L <- length(links)
  if (!cp04_reduction_count(L, 1L) || !is.numeric(links) || is.object(links) || !is.null(dim(links)) ||
      any(!is.finite(links)) || any(links != floor(links)) ||
      any(links < 1L | links > length(backend$links)) || anyDuplicated(links) > 0L)
    abort("Invalid explicit covariance link batch.")
  ipw_owner_int_bound(retained, "retained", min_val = 0)
  if (is.null(legacy)) legacy <- cp04_stream_plan_legacy(backend, links, retained)
  if (!is.list(legacy) || !is.list(legacy$blocks) || !length(legacy$blocks))
    cp04_abort("Invalid legacy covariance plan.", "integrity")
  if (length(legacy$blocks) <= 1L) return(NULL)
  cp04_retained_if_validate_backend(backend)
  st <- cp04_backend_domain(backend); ctrl <- backend$control
  B <- legacy$chunk_rows
  if (!cp04_reduction_count(B, 1L) || B > st$N || B > ctrl$max_chunk_rows)
    cp04_abort("Invalid legacy covariance row partition.", "integrity")
  ipw_owner_int_bound(legacy$work, "legacy covariance work", min_val = 0)
  boundary <- cp04_retained_if_boundary_cost(backend, links, B)
  tryCatch({
    LL <- ipw_owner_check_matrix_dims(L, L, ctrl$max_matrix_entries, "CP04 complete Gram")
    BL <- ipw_owner_check_matrix_dims(B, L, ctrl$max_matrix_entries, "CP04 full-width covariance chunk")
    gram_entries <- ipw_owner_safe_product(4, LL, ctrl$max_matrix_entries, "CP04 complete Gram buffers")
    # In addition to the legacy four-Gram charge, reserve the full-width
    # sweep expansion/centering and both crossprod slices (6 B L), two extra
    # Gram-sized product/copy buffers, and linear mean/index/output scratch.
    # Producer output AND scratch stay charged concurrently with all of it.
    consumer <- ipw_owner_safe_sum(c(
      ipw_owner_safe_product(6, BL), ipw_owner_safe_product(2, LL),
      ipw_owner_safe_product(16, L), 256, boundary$entries), ctrl$max_matrix_entries,
      "CP04 full-width covariance consumer workspace")
    ret_gram <- ipw_owner_safe_sum(c(retained, gram_entries, consumer), ctrl$max_matrix_entries,
      "CP04 covariance retained, Gram and consumer entries")
    # Test the issued legacy partition explicitly. Automatic replanning of a
    # different width is unnecessary when this exact B is safely admitted.
    B_exp <- cp04_plan_rows(backend, list(links), ret_gram, chunk_rows = B)
    if (!isTRUE(B_exp == B)) return(NULL)
    pre <- chained_full_if_budget(backend, B, links)
    peak <- ipw_owner_safe_sum(c(pre$output, pre$scratch, consumer), ctrl$max_matrix_entries,
      "CP04 covariance producer and consumer workspace")
    ipw_owner_safe_sum(c(retained, gram_entries, backend$accounting$owned_entries,
      backend$accounting$borrowed_entries, peak), ctrl$max_matrix_entries,
      "CP04 producer and downstream simultaneous entries")
    chunk_work <- ipw_owner_safe_sum(pre$work, ctrl$max_work_ops, "CP04 full-width producer chunk work")
    replay_work <- ipw_owner_safe_product(ceiling(st$N / B), chunk_work, ctrl$max_work_ops,
      "CP04 requested replay cumulative work")
    work <- ipw_owner_safe_sum(c(ipw_owner_safe_product(2, replay_work),
      ipw_owner_safe_product(st$N, gram_entries), boundary$work),
      ctrl$max_work_ops, "CP04 covariance cumulative work")
    if (work > legacy$work) return(NULL)
    cp04_budget(ctrl, c(retained, backend$accounting$owned_entries, backend$accounting$borrowed_entries),
      c(peak, gram_entries), work, gram_entries, "CP04 complete Gram")
    list(work = work, peak = peak, blocks = split(seq_len(L), rep.int(1L, L)), L = L, chunk_rows = B,
      consumer_workspace_entries = consumer)
  }, cdid_ipw_owner_bound_error = cp04_capacity_decline)
}

cp04_stream_plan <- function(backend, links, retained = 0) {
  legacy <- cp04_stream_plan_legacy(backend, links, retained)
  full <- cp04_stream_plan_full_block(backend, links, retained, legacy)
  if (is.null(full)) legacy else full
}

cp04_stream_covariance <- function(backend, links, retained = 0) {
  bound <- cp04_retained_if_unbind(backend); backend <- bound$backend; store <- bound$store
  completed <- FALSE
  on.exit(if (!completed && is.environment(store) && identical(store$kind, "cp04_retained_if"))
    cp04_retained_if_release(store), add = TRUE)
  st <- cp04_backend_domain(backend); ctrl <- backend$control
  plan <- cp04_stream_plan(backend, links, retained)
  work <- plan$work; peak <- plan$peak; blocks <- plan$blocks; L <- plan$L
  # CP08 retained moment IF: a bound store must describe this one link block and
  # the original row partition; its reserved entries join the concurrent check.
  store_entries <- 0
  if (!is.null(store)) {
    if (exists("phase", store, inherits = FALSE) && !is.null(store$phase) && !identical(store$phase, "covariance"))
      cp04_abort("Retained moment IF phase storage serves only its admitted phase.", "integrity")
    if (length(blocks) != 1L)
      cp04_abort("Retained moment IF storage serves exactly one covariance link block.", "integrity")
    cp04_retained_if_check(store, backend, links, plan$chunk_rows, "core")
    store_entries <- store$entries
  }
  cp04_budget(ctrl, c(retained, backend$accounting$owned_entries, backend$accounting$borrowed_entries, store_entries),
    c(peak, 4 * L * L), work, 4 * L * L, "CP04 complete Gram")
  producer <- if (is.null(store)) backend else cp04_retained_if_bind(backend, store)
  means <- numeric(L)
  for (ix in blocks) cp04_replay_links(producer, links[ix], function(b) {
    means[ix] <<- means[ix] + colSums(b$phi) / st$N
  }, chunk_rows = plan$chunk_rows)
  V <- matrix(0, L, L)
  for (a in seq_along(blocks)) for (b in a:length(blocks)) {
    ia <- blocks[[a]]; ib <- blocks[[b]]; ids <- unique(c(ia, ib)); gram <- matrix(0, length(ia), length(ib))
    cp04_replay_links(producer, links[ids], function(z) {
      x <- sweep(z$phi, 2L, means[ids], "-")
      gram <<- gram + (if (a == b) crossprod(x) else
        crossprod(x[, match(ia, ids), drop = FALSE], x[, match(ib, ids), drop = FALSE])) / st$n_A
    }, chunk_rows = plan$chunk_rows)
    if (a == b) {
      for (j in seq_along(ia)) if (j > 1L) gram[j, seq_len(j - 1L)] <- gram[seq_len(j - 1L), j]
      V[ia, ia] <- gram
    } else {
      V[ia, ib] <- gram; V[ib, ia] <- t(gram)
    }
  }
  estimates <- vapply(links, function(j) chained_full_link_estimate(backend, j), numeric(1))
  out <- list(estimates = estimates, covariance = V, mean_phi = means,
    work_ops_bound = work, producer_workspace_entries = peak, chunk_rows = plan$chunk_rows)
  completed <- TRUE
  out
}

# D60: every reference-input backend's diagnostics carry its recipe table,
# extraction conventions, provenance, frozen bases and recipe repeat rates.
cp04_public_diagnostics <- function(res, control) {
  out <- cp04_public_diagnostics_route(res, control)
  b <- res$nuisance_registry$cp04
  if (!is.null(b$reference_ref) && is.null(out$reference_inputs))
    out$reference_inputs <- chained_reference_record(chained_backend_reference(b))
  out
}

cp04_public_diagnostics_route <- function(res, control) {
  b <- res$nuisance_registry$cp04
  if (inherits(b, "cdid_reference_ipw_backend")) return(reference_ipw_public_diagnostics(res, control))
  if (inherits(b, "cdid_augmented_backend")) return(augmented_public_diagnostics(res, control))
  if (inherits(b, "cdid_imputation_backend")) return(imputation_public_diagnostics(res, control))
  if (inherits(b, c("cdid_dr_baseline", "cdid_dr_corollary1"))) return(dr_public_diagnostics(res, control))
  if (inherits(b, c("cdid_or_baseline", "cdid_or_corollary1"))) return(or_public_diagnostics(res, control))
  # Reserve the compact copy before constructing lists; source factors are
  # borrowed, without owner/design/IF copies. Conservative header allowance.
  bound <- 4 * b$accounting$owned_entries + 4096
  cp04_budget(control$owner_control, c(max(c(0, res$resource_accounting$simultaneous_entries,
    res$resource_accounting$initial_envelope_entries)), bound),
    4 * bound, 8 * bound, bound, "CP04 public diagnostic copy")
  records <- list(selection = lapply(seq_along(b$selections), function(j)
    c(list(id = j, spec = b$selection_specs[[j]]), b$selections[[j]])),
    treatment = lapply(seq_along(b$treatments), function(j) c(list(id = j, spec = b$treatment_specs[[j]]), b$treatments[[j]])))
  bytes <- unname(cp04_payload(records)["entries"])
  if (bytes > bound) cp04_abort("Diagnostic copy bound underestimated.", "bound")
  list(status = "available", method = "corollary1",
    selection_model = if (is.null(b$selection_model)) "pair_logit" else b$selection_model,
    selection_pooling = if (is.null(b$selection_pooling)) "separate" else b$selection_pooling,
    records = records,
    frame_counts = list(n_full_frame_units = res$reporting$inference$n_total, n_pair_eligible_A_units = res$frame_owner$n_a_units),
    links = b$links, accounting = b$accounting,
    auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
    interpretation = "Actual U/cohort/pair selection scores shared across post/auxiliary/NT/NYT uses; complete F rows, N-normalized information factors. No numerical certificate.")
}

# Exact nonnegative integer-valued scalar count (no tolerance, no `abs`).
cp04_reduction_count <- function(x, min_val = 0L) {
  is.numeric(x) && !is.object(x) && length(x) == 1L && is.finite(x) &&
    x == floor(x) && x >= min_val
}

# CP08 stage 6 (extended): post-fit complete-support structural reduction map
# for baseline OR/DR under optimal GMM weighting and PT-post, for the
# intercept-only (p = 0) recipe and for structurally common fixed-X outcome
# recipes (p >= 1) with ONE global shared slope, including a nonconstant shared
# DR treatment propensity. The licensing identity is the fitted formulas of
# or_pooled_fit_common / or_pooled_link_moments and dr_baseline_treatment_fit:
# on the common complete roster with constant weights every link's OR/DR
# moment is a difference of common weighted date means minus
# L(s,t) * (common covariate contrast)' beta with the ORIGINAL global beta and,
# for DR, the ONE shared odds functional; L(s,t) = i_t - i_s on the stored
# modeled calendar, so cycles telescope. Its derivative (the full original-unit
# IF, including the global beta IF and the propensity score IF) is linear in
# the same date-indexed terms and telescopes along the same forest paths.
# Eligibility is decided from exact sealed facts and bitwise identities only;
# no numerical near-equality, no pairwise/ordinary-DiD refit and no caught
# failure admits a map. Returns NULL when the fitted backend or design is
# structurally ineligible; never fabricates an IPW certificate or alters
# backend fits, diagnostics or auxiliary link indices.
cp04_complete_reduction <- function(backend, point_links, design, design_link_index,
                                    frame_owner, bounds, owner_control,
                                    retained_pipeline_entries = 0, external_work_ops = 0) {
  if (is.null(backend) || !inherits(backend, c("cdid_or_baseline", "cdid_dr_baseline"))) {
    return(NULL)
  }
  is_dr <- inherits(backend, "cdid_dr_baseline")
  if (is_dr) validate_dr_baseline_backend(backend) else validate_or_baseline_backend(backend)
  cs <- backend$complete_support
  recipe <- backend$outcome_recipe
  shared <- backend$shared
  if (is.null(cs) || !is.null(backend$reference_ref) ||
      length(cs$target) != 1L || length(cs$control) != 1L ||
      !is.list(recipe) || !inherits(recipe, "cdid_or_pooled_recipe") || !isTRUE(recipe$intercept) ||
      !is.list(shared)) {
    return(NULL)
  }
  # Recipe identity: the explicit private recipe roles must be self-consistent
  # (one intercept column, slope columns = recipe$columns[slope_index]).
  p_rec <- recipe$p
  if (!cp04_reduction_count(p_rec, 0L) ||
      !is.character(recipe$columns) || !length(recipe$columns) || anyNA(recipe$columns) ||
      anyDuplicated(recipe$columns) > 0L ||
      !is.numeric(recipe$slope_index) || length(recipe$slope_index) != p_rec ||
      !identical(as.integer(recipe$slope_index), which(recipe$columns %in% recipe$slope_columns)) ||
      !identical(recipe$slope_columns, recipe$columns[recipe$slope_index]) ||
      !identical(as.integer(recipe$intercept_index), setdiff(seq_along(recipe$columns), as.integer(recipe$slope_index))) ||
      length(recipe$intercept_index) != 1L) {
    return(NULL)
  }
  # Shared-fit identity: ONE global kept slope basis of dimension p <= p_rec,
  # frozen during differentiation, with the fitted global beta of that length.
  p <- shared$p
  keep <- shared$keep
  if (!cp04_reduction_count(p, 0L) || p > p_rec ||
      !is.numeric(keep) || length(keep) != p || anyNA(keep) ||
      (p > 0L && (any(keep != floor(keep)) || any(keep < 1L | keep > p_rec) ||
                  is.unsorted(keep, strictly = TRUE))) ||
      !is.numeric(shared$beta) || length(shared$beta) != p || any(!is.finite(shared$beta)) ||
      !is.list(shared$basis) || !isTRUE(shared$basis$frozen_during_differentiation) ||
      !identical(as.integer(shared$basis$kept_index), as.integer(keep)) ||
      !identical(shared$basis$slope_columns, recipe$slope_columns) ||
      !identical(shared$basis$kept_columns, recipe$slope_columns[keep]) ||
      !identical(cs$route, "complete_support") ||
      !cp04_reduction_count(cs$p, 0L) || !identical(as.integer(cs$p), as.integer(p_rec)) ||
      !cp04_reduction_count(cs$n_z, 1L)) {
    return(NULL)
  }
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  frame_ref_obj <- if (is.environment(reader$frame_ref) &&
                       exists("frame_owner", reader$frame_ref, inherits = FALSE)) {
    get("frame_owner", reader$frame_ref, inherits = FALSE)
  } else NULL
  if (!inherits(frame_owner, "cdid_ipw_frame_owner") || !is.list(frame_ref_obj) ||
      !identical(frame_ref_obj$ref, frame_owner$ref) ||
      !identical(reader$panel_ref, frame_owner$panel_ref)) {
    return(NULL)
  }
  # Fixed-X design binding: the recipe reads exactly the sealed frame design
  # Z_a columns (fixed per unit by construction); the complete-support facts
  # counted exactly the used design columns.
  Z_a <- frame_ref_obj$Z_a
  z_names <- colnames(Z_a)
  if (!is.matrix(Z_a) || !is.double(Z_a) || ncol(Z_a) < 1L || is.null(z_names) ||
      anyNA(z_names) || anyDuplicated(z_names) > 0L) {
    return(NULL)
  }
  p_trt <- 1L
  trt <- NULL
  if (is_dr) {
    trt <- backend$treatment_recipe
    ofi <- backend$or_frame_index
    if (!is.list(trt) || !inherits(trt, "cdid_dr_treatment_recipe") ||
        !cp04_reduction_count(trt$p, 1L) || !is.integer(trt$index) ||
        length(trt$index) != trt$p || anyNA(trt$index) || anyDuplicated(trt$index) > 0L ||
        any(trt$index < 1L | trt$index > ncol(Z_a)) ||
        !identical(trt$columns, z_names[trt$index]) ||
        !is.integer(ofi) || length(ofi) != length(recipe$columns) || anyNA(ofi) ||
        any(ofi < 1L | ofi > ncol(Z_a)) || anyDuplicated(ofi) > 0L ||
        !identical(recipe$columns, z_names[ofi]) ||
        !identical(as.integer(cs$n_z), length(unique(c(ofi, trt$index))))) {
      return(NULL)
    }
    p_trt <- as.integer(trt$p)
  } else if (!identical(recipe$columns, z_names) ||
             !identical(as.integer(cs$n_z), length(recipe$columns))) {
    return(NULL)
  }
  # Global slope fit facts (p >= 1): full rank on the kept basis, one block per
  # regression, finite factor, and the forcing/control-support bookkeeping of
  # the global beta IF spans every control block on the never-treated roster.
  n_sup <- length(backend$supports)
  n_reg <- length(backend$regression_specs)
  if (n_sup < 1L || n_reg < 1L || length(backend$support_specs) != n_sup ||
      length(backend$regressions) != n_reg) {
    return(NULL)
  }
  # Admit the new all-original G/S binding scans and their simultaneous
  # temporary vectors before lapply/unique or finite-matrix copies. The
  # intercept-only route performs no new global-slope scan and keeps its
  # established descriptor/work charges exactly.
  binding_entries <- binding_work <- 0
  if (p_rec > 0L || p_trt > 1L) {
    binding_ctrl <- validate_ipw_owner_control(owner_control)
    binding_cap <- min(binding_ctrl$max_retained_entries, binding_ctrl$max_matrix_entries)
    binding_entries <- ipw_owner_safe_sum(c(16 * as.double(n_reg), 16 * as.double(n_sup),
      8 * as.double(p_rec), 8 * as.double(p) * as.double(n_sup),
      8 * as.double(p)^2, 8 * as.double(p_trt)^2, 128), binding_cap,
      "CP04 shared-X original G/S binding temporary entries")
    binding_work <- ipw_owner_safe_sum(c(256 * as.double(n_reg), 256 * as.double(n_sup),
      256 * as.double(p_rec), 256 * as.double(p) * as.double(n_sup),
      256 * as.double(p)^2, 256 * as.double(p_trt)^2, 128), binding_ctrl$max_work_ops,
      "CP04 shared-X original G/S binding scan work")
    binding_retained <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
    binding_external <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
    ipw_owner_safe_sum(c(max(binding_retained,
      backend$accounting$owned_entries + backend$accounting$borrowed_entries), binding_entries), binding_cap,
      "pipeline/backend plus shared-X binding temporaries before scans")
    ipw_owner_safe_sum(c(max(binding_external, backend$accounting$work_ops_bound), binding_work),
      binding_ctrl$max_work_ops, "pipeline/backend plus shared-X binding work before scans")
  }
  if (p > 0L) {
    reg_sids <- lapply(backend$regression_specs, `[[`, "support_ids")
    if (any(lengths(reg_sids) != 1L)) {
      return(NULL)
    }
    ctrl_sids <- sort(unique(as.integer(unlist(reg_sids))))
    if (!or_baseline_native_index(ctrl_sids, n_sup) ||
        !all(vapply(backend$support_specs[ctrl_sids], function(s)
          identical(as.integer(s$cohort), as.integer(cs$control)), logical(1)))) {
      return(NULL)
    }
    fit <- shared$fit
    if (!is.list(fit) || !inherits(fit, "cdid_or_pooled_wls_fit") ||
        !identical(as.integer(fit$p), as.integer(p)) ||
        !identical(as.integer(fit$rank), as.integer(p)) ||
        !identical(as.integer(fit$n_blocks), as.integer(n_reg)) ||
        !is.matrix(fit$R) || !is.numeric(fit$R) || ncol(fit$R) != p || nrow(fit$R) < p ||
        any(!is.finite(fit$R)) ||
        !is.numeric(fit$pivot) || length(fit$pivot) != p || anyNA(fit$pivot) ||
        !identical(sort(as.integer(fit$pivot)), seq_len(p)) ||
        !is.matrix(shared$forcing_beta) || !identical(dim(shared$forcing_beta), c(as.integer(p), as.integer(n_sup))) ||
        any(!is.finite(shared$forcing_beta)) ||
        !identical(as.integer(shared$control_support_ids), ctrl_sids)) {
      return(NULL)
    }
  }
  if (!is.list(point_links) || is.object(point_links) || length(point_links) < 1L) {
    return(NULL)
  }
  M <- length(point_links)
  b <- validate_ipw_certificate_bounds(bounds)
  ctrl <- validate_ipw_owner_control(owner_control)
  retained_pipeline_entries <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
  external_work_ops <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
  if (M > b$max_links || 1L > b$max_classes || !inherits(design, "cdid_chained_design") ||
      !is.data.frame(design$vertices)) {
    cdid_abort("CP04 structural reduction exceeds certificate link or class bounds.",
               c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  cspec <- design$calendar_spec
  if (is.null(cspec) || !inherits(cspec, "cdid_calendar_spec") ||
      !identical(if (is.null(cspec$pt_assumption)) "post" else cspec$pt_assumption, "post")) {
    return(NULL)
  }

  # Reserve linear descriptor/forest buffers before allocation, composed with live charges.
  cap <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  V_n <- as.double(nrow(design$vertices))
  own_entries <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  base_entries <- max(retained_pipeline_entries, own_entries)
  base_work <- max(external_work_ops, backend$accounting$work_ops_bound)
  linear <- ipw_owner_safe_sum(c(256 * as.double(M), 64 * V_n, 256, binding_entries), cap,
    "CP04 complete-support reduction linear descriptor and forest plan")
  # Issued complete dates bound every current endpoint. Reserve both path
  # passes, sorting and incidence checks before any descriptor/forest work.
  date_bound <- min(as.double(length(cs$dates)), 2 * as.double(M))
  path_bound <- min(as.double(M), max(0, date_bound - 1))
  E_bound <- ipw_owner_safe_product(M, path_bound,
    label = "CP04 complete-support reduction path entry upper bound")
  # Shared fixed-X descriptor comparisons per link (roster design sums on the
  # full recipe slope set, kept-basis vectors and, for DR, the shared logit's
  # coefficient, factor and loading arrays). Zero for the intercept-only
  # route, whose existing charges are unchanged. No fit, replay or beta-IF
  # work is charged or removed here: the retained links' covariance/replay
  # charges, including the all-block global beta IF cost, stay in the
  # pipeline's own envelopes.
  slope_check_work <- if (p_rec > 0L || p_trt > 1L) {
    64 * as.double(M) * (2 * as.double(p_rec) + 2 * as.double(p) +
      as.double(p_trt) * (as.double(p_trt) + 8))
  } else 0
  preflight_work <- ipw_owner_safe_sum(
    c(256 * as.double(M) * (date_bound + ceiling(log2(M + 1))), 64 * E_bound, slope_check_work, binding_work),
    ctrl$max_work_ops, "CP04 complete-support reduction full preflight work")
  ipw_owner_safe_sum(c(base_entries, linear), cap,
    "pipeline/backend plus CP04 complete-support reduction linear plan")
  ipw_owner_safe_sum(linear, ctrl$max_output_entries,
    "CP04 complete-support reduction linear output bound")
  ipw_owner_safe_sum(c(base_work, preflight_work), ctrl$max_work_ops,
    "pipeline/backend plus CP04 complete-support reduction work before construction")

  dli <- validate_ipw_design_for_classes(design, M, design_link_index)
  validate_ipw_design_incidence(design)
  DL <- design$links
  DV <- design$vertices

  target_cohort_val <- as.double(reader$state$cohorts[cs$target])
  control_cohort_val <- as.double(reader$state$cohorts[cs$control])
  if (length(target_cohort_val) != 1L || length(control_cohort_val) != 1L ||
      !is.finite(target_cohort_val) || !is.infinite(control_cohort_val)) {
    return(NULL)
  }

  cp04_idx <- integer(M)
  s_vec <- t_vec <- numeric(M)
  L_vec <- numeric(M)
  vs <- vt <- integer(M)
  sids_t <- sids_c <- rids <- tids <- integer(M)
  src_fk_first <- NULL
  dl_pools <- if ("control_pool" %in% names(DL)) design_pool_cache(DL$control_pool) else NULL

  for (m in seq_len(M)) {
    pt <- point_links[[m]]
    if (!is.list(pt) || !or_baseline_native_index(pt$cp04_index, length(backend$links))) {
      return(NULL)
    }
    idx <- as.integer(pt$cp04_index)
    cp04_idx[m] <- idx
    bl <- backend$links[[idx]]
    if (!identical(pt$link_id, bl$link_id) ||
        !is.character(bl$source_functional_key) || length(bl$source_functional_key) != 1L ||
        is.na(bl$source_functional_key) || !nzchar(bl$source_functional_key)) {
      return(NULL)
    }
    if (m == 1L) {
      src_fk_first <- bl$source_functional_key
    } else if (!identical(bl$source_functional_key, src_fk_first)) {
      return(NULL)
    }

    sid <- bl$support_id
    rid <- bl$regression
    if (!or_baseline_native_index(sid, length(backend$support_specs)) ||
        !or_baseline_native_index(rid, length(backend$regression_specs))) {
      return(NULL)
    }
    ss <- backend$support_specs[[sid]]
    rs <- backend$regression_specs[[rid]]
    if (!identical(as.integer(bl$pair_id), as.integer(ss$pair_id)) ||
        !identical(as.integer(bl$pair_id), as.integer(rs$pair_id)) ||
        !identical(as.integer(bl$cohort), as.integer(ss$cohort)) ||
        !identical(as.integer(ss$cohort), as.integer(cs$target)) ||
        !identical(as.integer(rs$cohorts), as.integer(cs$control)) ||
        length(rs$support_ids) != 1L ||
        !or_baseline_native_index(rs$support_ids[1L], length(backend$support_specs))) {
      return(NULL)
    }
    csid <- as.integer(rs$support_ids[1L])
    css <- backend$support_specs[[csid]]
    if (!identical(as.integer(css$pair_id), as.integer(bl$pair_id)) ||
        !identical(as.integer(css$cohort), as.integer(cs$control)) ||
        !identical(ss$pair, css$pair) ||
        !identical(pt$pair, ss$pair) ||
        !is.list(pt$cohorts) ||
        !identical(as.double(pt$cohorts$g), target_cohort_val) ||
        !identical(as.double(pt$cohorts$K), control_cohort_val)) {
      return(NULL)
    }
    ds <- as.integer(ss$pair$s_calendar_index)
    dt <- as.integer(ss$pair$t_calendar_index)
    if (!all(c(ds, dt) %in% cs$dates) ||
        !identical(as.double(cspec$calendar[ds]), as.double(ss$pair$s)) ||
        !identical(as.double(cspec$calendar[dt]), as.double(ss$pair$t))) {
      return(NULL)
    }
    # Duration identity: the fitted L of this link's supports and block is the
    # stored-calendar transition count i_t - i_s (or_pooled_transitions), the
    # quantity whose linearity makes the slope terms telescope along cycles.
    L_m <- as.double(dt) - as.double(ds)
    if (!(L_m >= 1) || !identical(as.double(ss$transitions), L_m) ||
        !identical(as.double(css$transitions), L_m) ||
        !identical(as.double(rs$transitions), L_m)) {
      return(NULL)
    }
    # Slope descriptor dimensions of the fitted identity: the link's treated
    # covariate contrast c = L * zbar_T(keep) and its block's mean design
    # L * zbar_C(keep) both live on the global kept basis.
    if (!is.numeric(bl$c) || length(bl$c) != p || any(!is.finite(bl$c))) {
      return(NULL)
    }

    j <- dli[m]
    if (!identical(as.double(DL$g[j]), target_cohort_val) ||
        !identical(as.double(DL$s[j]), as.double(ss$pair$s)) ||
        !identical(as.double(DL$t[j]), as.double(ss$pair$t)) ||
        !identical(DL$functional_key[j], src_fk_first)) {
      return(NULL)
    }
    if (!is.null(dl_pools)) {
      cp_vals <- dl_pools$values(j)
      if (is.null(cp_vals) || !identical(as.double(cp_vals), control_cohort_val)) {
        return(NULL)
      }
    }
    v_start <- as.integer(DL$start_vertex_index[j])
    v_end <- as.integer(DL$end_vertex_index[j])
    if (is.na(v_start) || is.na(v_end) || v_start == v_end ||
        !identical(as.double(DV$g[v_start]), target_cohort_val) ||
        !identical(as.double(DV$g[v_end]), target_cohort_val) ||
        !identical(as.double(DV$time[v_start]), as.double(ss$pair$s)) ||
        !identical(as.double(DV$time[v_end]), as.double(ss$pair$t))) {
      return(NULL)
    }

    s_vec[m] <- as.double(ss$pair$s)
    t_vec[m] <- as.double(ss$pair$t)
    L_vec[m] <- L_m
    vs[m] <- v_start
    vt[m] <- v_end
    sids_t[m] <- as.integer(sid)
    sids_c[m] <- csid
    rids[m] <- as.integer(rid)
    if (is_dr) {
      if (!or_baseline_native_index(bl$model, length(backend$treatments)) ||
          !or_baseline_native_index(bl$model, length(backend$treatment_specs))) {
        return(NULL)
      }
      tids[m] <- as.integer(bl$model)
    }
  }
  if (anyDuplicated(cp04_idx) > 0L) {
    return(NULL)
  }

  # Verify exact same cohort rates/masses/counts across pairs, with rates == 1 and positive denominators.
  ut1 <- backend$supports[[sids_t[1L]]]
  uc1 <- backend$supports[[sids_c[1L]]]
  rg1 <- backend$regressions[[rids[1L]]]
  if (!identical(ut1$rate, 1) || !identical(uc1$rate, 1) ||
      !identical(ut1$v, 1) || !identical(uc1$v, 1) ||
      !is.finite(ut1$massA) || ut1$massA <= 0 ||
      !is.finite(uc1$massA) || uc1$massA <= 0 ||
      !is.finite(ut1$pair_mass) || ut1$pair_mass <= 0 ||
      !is.finite(uc1$pair_mass) || uc1$pair_mass <= 0 ||
      !is.finite(ut1$repeat_mass) || ut1$repeat_mass <= 0 ||
      !is.finite(uc1$repeat_mass) || uc1$repeat_mass <= 0 ||
      !identical(ut1$pair_mass, ut1$repeat_mass) ||
      !identical(uc1$pair_mass, uc1$repeat_mass) ||
      !is.numeric(ut1$n_pair) || length(ut1$n_pair) != 1L || !is.finite(ut1$n_pair) || ut1$n_pair <= 0 ||
      !is.numeric(uc1$n_pair) || length(uc1$n_pair) != 1L || !is.finite(uc1$n_pair) || uc1$n_pair <= 0 ||
      !identical(as.integer(ut1$n_pair), as.integer(ut1$n_repeat)) ||
      !identical(as.integer(uc1$n_pair), as.integer(uc1$n_repeat)) ||
      !is.finite(rg1$mass) || rg1$mass <= 0) {
    return(NULL)
  }
  # Common covariate contrast: both cohort rosters are complete and fixed-X, so
  # every support's raw design sum over its roster (length p_rec, the full
  # recipe slope set) is one object per cohort; the complete fit stores it per
  # support, and the common complete roster plus bitwise equality bind
  # zbar_T(keep) and zbar_C(keep) across links without any tolerance.
  if (!is.numeric(ut1$sumZ) || length(ut1$sumZ) != p_rec || any(!is.finite(ut1$sumZ)) ||
      !is.numeric(uc1$sumZ) || length(uc1$sumZ) != p_rec || any(!is.finite(uc1$sumZ)) ||
      as.double(ut1$n_pair) + as.double(uc1$n_pair) > as.double(reader$state$N)) {
    return(NULL)
  }
  # The sealed complete roster and exact support/recipe bindings establish
  # the common fitting operator. QR-derived regression masses can round
  # differently from support sums; no near-equality admission rule is needed.
  for (m in seq_len(M)) {
    ut <- backend$supports[[sids_t[m]]]
    uc <- backend$supports[[sids_c[m]]]
    rg <- backend$regressions[[rids[m]]]
    if (!identical(ut$rate, ut1$rate) || !identical(ut$v, ut1$v) ||
        !identical(ut$massA, ut1$massA) || !identical(ut$pair_mass, ut1$pair_mass) ||
        !identical(ut$repeat_mass, ut1$repeat_mass) ||
        !identical(ut$n_pair, ut1$n_pair) || !identical(ut$n_repeat, ut1$n_repeat) ||
        !identical(uc$rate, uc1$rate) || !identical(uc$v, uc1$v) ||
        !identical(uc$massA, uc1$massA) || !identical(uc$pair_mass, uc1$pair_mass) ||
        !identical(uc$repeat_mass, uc1$repeat_mass) ||
        !identical(uc$n_pair, uc1$n_pair) || !identical(uc$n_repeat, uc1$n_repeat) ||
        !is.finite(rg$mass) || rg$mass <= 0) {
      return(NULL)
    }
    # Shared fixed-X identity per link: same roster design sums (bitwise), the
    # support/block L equal to this link's transition count, and the block's
    # kept-basis mean design of the global fit's dimension.
    if (!identical(ut$sumZ, ut1$sumZ) || !identical(uc$sumZ, uc1$sumZ) ||
        !identical(as.double(ut$L), L_vec[m]) || !identical(as.double(uc$L), L_vec[m]) ||
        !identical(as.double(rg$L), L_vec[m]) ||
        !is.numeric(rg$mean_d) || length(rg$mean_d) != p || any(!is.finite(rg$mean_d)) ||
        !is.numeric(rg$alpha) || length(rg$alpha) != 1L || !is.finite(rg$alpha)) {
      return(NULL)
    }
  }

  if (is_dr) {
    # ONE shared odds functional: every link's treatment model is the same
    # fitted logit on the same complete two-cohort risk set (bitwise identical
    # coefficients, information factor and support factors; the identities
    # proved at fit and re-checked by dr_baseline_native_plan), of the recipe's
    # dimension p_trt >= 1 on any fixed-X columns. The PS-IF loadings
    # (partials$beta, 4 x p_trt) and the OR-slope contrast (d_bar_diff, length
    # p) are the date-linear terms whose derivative telescopes with the map.
    s1 <- backend$treatment_specs[[tids[1L]]]
    f1 <- backend$treatments[[tids[1L]]]
    v_of <- function(s) vapply(backend$supports[s$support_ids], `[[`, numeric(1), "v")
    v1 <- v_of(s1)
    m_coh <- length(s1$support_ids)
    if (m_coh != 2L || !identical(as.integer(s1$cohorts), as.integer(c(cs$target, cs$control))) ||
        !is.numeric(f1$coefficients) || length(f1$coefficients) != p_trt ||
        any(!is.finite(f1$coefficients)) ||
        !is.numeric(f1$H_chol) || length(f1$H_chol) != p_trt * p_trt || any(!is.finite(f1$H_chol))) {
      return(NULL)
    }
    for (m in seq_len(M)) {
      sk <- backend$treatment_specs[[tids[m]]]
      fk_fit <- backend$treatments[[tids[m]]]
      if (!identical(sk$cohorts, s1$cohorts) ||
          !identical(as.integer(sk$support_ids), c(sids_t[m], sids_c[m])) ||
          !identical(as.integer(sk$regression), rids[m]) ||
          !identical(fk_fit$coefficients, f1$coefficients) ||
          !identical(fk_fit$H_chol, f1$H_chol) ||
          !identical(v_of(sk), v1) ||
          !identical(dim(fk_fit$cross), c(p_trt, m_coh)) ||
          !identical(dim(fk_fit$partials$beta), c(4L, p_trt)) ||
          !identical(dim(fk_fit$partials$upstream), c(4L, m_coh)) ||
          !is.numeric(fk_fit$d_bar_diff) || length(fk_fit$d_bar_diff) != p ||
          any(!is.finite(fk_fit$d_bar_diff)) ||
          !is.numeric(fk_fit$means) || anyNA(fk_fit$means) ||
          !all(c("T_dr", "M_dr", "dp", "dm", "mu_dr") %in% names(fk_fit$means))) {
        return(NULL)
      }
    }
  }

  members <- order(s_vec, t_vec, seq_len(M), method = "radix")
  dates <- sort(unique(c(s_vec[members], t_vec[members])))
  if (length(dates) > b$max_class_dates) {
    cdid_abort("CP04 class date count exceeds `max_class_dates`.",
               c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  a_idx <- match(s_vec[members], dates)
  b_idx <- match(t_vec[members], dates)
  if (anyDuplicated(paste(a_idx, b_idx)) > 0L) {
    return(NULL)
  }
  vmap <- rep(NA_integer_, length(dates))
  for (side in 1:2) {
    di <- if (side == 1L) a_idx else b_idx
    vv <- if (side == 1L) vs[members] else vt[members]
    for (h in seq_along(di)) {
      if (is.na(vmap[di[h]])) {
        vmap[di[h]] <- vv[h]
      } else if (vmap[di[h]] != vv[h]) {
        return(NULL)
      }
    }
  }
  if (anyDuplicated(vmap) > 0L) {
    return(NULL)
  }

  forest <- incidence_forest(length(dates), a_idx, b_idx, integer(0L))
  tree <- sort(forest$parent_link[forest$parent_link > 0L])
  dropped <- setdiff(seq_along(members), tree)
  depth <- ipw_forest_depth(forest)
  total_E <- 0L
  for (h in dropped) {
    total_E <- total_E + as.integer(ipw_forest_path(a_idx[h], b_idx[h], forest, depth, TRUE))
  }
  if (total_E > b$max_E_entries) {
    cdid_abort("CP04 reduction map entries exceed `max_E_entries`.",
               c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  additional_entries <- ipw_owner_safe_sum(c(linear, 16 * as.double(total_E)), cap,
    "CP04 complete-support reduction construction/sort/copy buffers and map")
  work_ops_bound <- ipw_owner_safe_sum(
    c(256 * as.double(M) * (length(dates) + ceiling(log2(M + 1))), 64 * as.double(total_E),
      slope_check_work, binding_work),
    ctrl$max_work_ops,
    "CP04 complete-support reduction forest/path/sort work"
  )
  ipw_owner_safe_sum(c(base_entries, additional_entries), cap,
    "pipeline/backend plus CP04 complete-support reduction map")
  ipw_owner_safe_sum(additional_entries, ctrl$max_output_entries,
    "CP04 complete-support reduction map output")
  ipw_owner_safe_sum(c(base_work, work_ops_bound), ctrl$max_work_ops,
    "pipeline/backend plus CP04 complete-support reduction work")

  E_i <- E_j <- integer(total_E)
  E_x <- double(total_E)
  pos <- 0L
  for (h in dropped) {
    path <- ipw_forest_path(a_idx[h], b_idx[h], forest, depth, FALSE)
    o <- order(members[path$edges])
    idx <- pos + seq_along(o)
    E_i[idx] <- as.integer(members[h])
    E_j[idx] <- as.integer(members[path$edges[o]])
    E_x[idx] <- as.double(path$signs[o])
    pos <- pos + length(o)
  }
  if (total_E > 0L) {
    ord_e <- order(E_i, E_j, method = "radix")
    E_i <- E_i[ord_e]
    E_j <- E_j[ord_e]
    E_x <- E_x[ord_e]
  }

  retained_index <- as.integer(sort(members[tree]))
  dropped_index <- as.integer(sort(members[dropped]))
  kind_str <- if (is_dr) "complete_support_DR" else "complete_support_OR"
  class_rec <- list(
    source_functional_key = src_fk_first,
    method = if (is_dr) "dr" else "or",
    g = target_cohort_val,
    K = control_cohort_val,
    dates = dates,
    vertex_index = vmap,
    link_index = as.integer(members),
    forest_link_index = as.integer(members[tree]),
    dropped_link_index = as.integer(members[dropped]),
    roots = dates[forest$roots],
    n_components = forest$n_components,
    cycle_rank = length(members) - (length(dates) - forest$n_components),
    n_target_units = as.integer(ut1$n_pair),
    n_control_units = as.integer(uc1$n_pair)
  )

  out <- list(
    n_links = as.integer(M),
    retained_index = retained_index,
    dropped_index = dropped_index,
    E = list(i = E_i, j = E_j, x = E_x),
    design_link_index = as.integer(dli),
    link_kind = rep(kind_str, M),
    classes = list(class_rec),
    backend_ref = backend$ref,
    context_ref = frame_owner$ref,
    accounting = list(
      additional_entries = additional_entries,
      work_ops_bound = work_ops_bound
    )
  )
  ipw_owner_seal(out, "cp04_complete_reduction")
}

cp04_expand_reduction_map <- function(map, x_retained) {
  if (!is.list(map) || is.object(map) || !ipw_owner_seal_ok(map, "cp04_complete_reduction")) {
    cdid_abort("Expected an unchanged sealed `cp04_complete_reduction` map.",
               c("cdid_ipw_certificate_error", "cdid_ipw_input_error"))
  }
  vec <- is.null(dim(x_retained))
  X <- if (vec) matrix(x_retained, ncol = 1L) else x_retained
  n_r <- length(map$retained_index)
  if (!is.numeric(X) || !is.matrix(X) || nrow(X) != n_r || any(!is.finite(X))) {
    cdid_abort("`x_retained` must be finite with one row per retained link.",
               c("cdid_ipw_input_error"))
  }
  out <- matrix(0, map$n_links, ncol(X))
  out[map$retained_index, ] <- X
  rpos <- match(map$E$j, map$retained_index)
  for (k in seq_along(map$E$i)) {
    i <- map$E$i[k]
    out[i, ] <- out[i, ] + map$E$x[k] * X[rpos[k], ]
  }
  if (vec) drop(out) else out
}
