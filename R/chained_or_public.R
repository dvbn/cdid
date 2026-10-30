# Public OR adapter: one frozen common extraction, shared one-period outcome-change
# control regression with local block intercepts, and the existing full-roster
# replay/GMM/reporting consumers.

or_public_rank_tol <- function() eval(formals(stats::lm.fit)$tol)

# R's lm.fit / dqrdc2 convention proposes aliases (tol=1e-7). A proposal is
# accepted only when every common A row reconstructs to arithmetic roundoff.
# This second check is a numerical discrepancy diagnostic, NOT an exact-span
# or effect-error certificate. Near dependencies are retained for the separate
# control-only fit/conditioning check. No selection-design columns are touched.
or_public_basis <- function(Z, active, resource = NULL, intercept = NULL, intercept_index = NULL) {
  n <- sum(active); p <- ncol(Z)
  if (!is.null(resource)) {
    ipw_owner_safe_sum(c(8 * length(Z), 64 * n * p, 64 * p * p),
      resource$numeric_workspace_bytes, "OR common design QR preprocessing bytes")
    ipw_owner_safe_sum(128 * n * p * p, resource$counted_work_limit,
      "OR common design QR preprocessing work")
  }
  X <- Z[active, , drop = FALSE]
  chained_request_need(n > 0L && p > 0L && all(is.finite(X)),
    "Outcome regression requires a finite common outcome design on repeat-observer rows.",
    "cdid_or_design_error")
  if (!is.null(intercept)) {
    has_intercept <- isTRUE(intercept)
    icpt_idx <- if (has_intercept) as.integer(if (is.null(intercept_index)) 1L else intercept_index) else integer(0)
  } else if (!is.null(attr(Z, "intercept")) || !is.null(attr(Z, "assign"))) {
    asg <- attr(Z, "assign")
    has_intercept <- if (!is.null(attr(Z, "intercept"))) isTRUE(attr(Z, "intercept")) else (0L %in% as.integer(asg))
    icpt_idx <- if (has_intercept) {
      if (!is.null(asg) && (0L %in% as.integer(asg))) which(as.integer(asg) == 0L)[1L] else 1L
    } else integer(0)
  } else {
    has_intercept <- p >= 1L && !is.null(colnames(X)) &&
      (colnames(X)[1L] %in% c("(Intercept)", "intercept")) && all(X[, 1L] == 1)
    icpt_idx <- if (has_intercept) 1L else integer(0)
  }
  if (has_intercept) {
    chained_request_need(length(icpt_idx) == 1L && !is.na(icpt_idx) &&
      icpt_idx >= 1L && icpt_idx <= p && all(X[, icpt_idx] == 1),
      "Declared formula intercept column must be a unit column on repeat-observer rows.",
      "cdid_or_design_error")
  }
  ord <- if (has_intercept && icpt_idx != 1L) c(icpt_idx, setdiff(seq_len(p), icpt_idx)) else seq_len(p)
  q <- qr(X[, ord, drop = FALSE], tol = or_public_rank_tol(), LAPACK = FALSE)
  retained <- sort(ord[q$pivot[seq_len(q$rank)]])
  proposed <- setdiff(seq_len(p), retained)
  map <- matrix(0, length(retained), p,
    dimnames = list(colnames(X)[retained], colnames(X)))
  discrepancy <- numeric(p); scaled <- numeric(p)
  # A zero-column span cannot fit the private WLS recipe; refuse explicitly.
  chained_request_need(length(retained) > 0L,
    paste0("Common outcome design has no nonzero retained columns: ", paste(colnames(X), collapse = ", "),
      ". No reduced control regression can be computed."), "cdid_or_design_error")
  map[, retained] <- diag(length(retained))
  threshold <- 64 * .Machine$double.eps * max(1L, p)
  if (length(proposed)) {
    kept_qr <- qr(X[, retained, drop = FALSE], tol = or_public_rank_tol(), LAPACK = FALSE)
    map[, proposed] <- qr.coef(kept_qr, X[, proposed, drop = FALSE])
    # Bound replay memory by columns; no second common n-by-p fitted matrix.
    for (j in proposed) {
      predicted <- as.vector(X[, retained, drop = FALSE] %*% map[, j])
      # One fixed least-squares residual correction counters accumulated QR
      # roundoff on long exact aliases. The discrepancy guard is unchanged.
      map[, j] <- map[, j] + qr.coef(kept_qr, X[, j] - predicted)
      predicted <- as.vector(X[, retained, drop = FALSE] %*% map[, j])
      difference <- abs(X[, j] - predicted)
      scale <- abs(X[, j]) + as.vector(abs(X[, retained, drop = FALSE]) %*% abs(map[, j]))
      discrepancy[j] <- max(difference)
      scaled[j] <- max(difference / pmax(scale, max(abs(X[, j])), .Machine$double.xmin))
    }
  }
  rejected <- proposed[!is.finite(scaled[proposed]) | scaled[proposed] > threshold]
  if (length(rejected)) {
    # Keep every nontrivial candidate, while still removing verified perfect
    # aliases. Their old reconstruction uses the same retained columns and
    # remains valid with zero loadings on the newly retained candidates.
    previous <- retained; retained <- sort(c(retained, rejected))
    expanded <- matrix(0, length(retained), p,
      dimnames = list(colnames(X)[retained], colnames(X)))
    expanded[match(previous, retained), ] <- map
    expanded[, retained] <- diag(length(retained))
    map <- expanded
  }
  aliased <- setdiff(seq_len(p), retained)
  Z_out <- Z[, retained, drop = FALSE]
  attr(Z_out, "intercept") <- NULL
  attr(Z_out, "assign") <- NULL
  slope_cols <- if (!is.null(colnames(X))) {
    if (has_intercept) colnames(X)[setdiff(retained, icpt_idx)] else colnames(X)[retained]
  } else character(0)
  outcome_recipe <- if (!is.null(colnames(X))) {
    or_pooled_outcome_recipe(colnames(X)[retained], has_intercept, slope_cols)
  } else NULL
  list(Z = Z_out, outcome_recipe = outcome_recipe, record = list(
    original_columns = colnames(X), retained_columns = colnames(X)[retained],
    aliased_columns = colnames(X)[aliased], retained_index = retained,
    intercept = has_intercept,
    intercept_column = if (has_intercept && !is.null(colnames(X))) colnames(X)[icpt_idx] else NULL,
    slope_columns = slope_cols,
    outcome_recipe = outcome_recipe,
    alias_map = map, proposed_aliases = colnames(X)[proposed],
    nontrivial_candidates_retained = colnames(X)[rejected],
    reconstruction_max_abs = stats::setNames(discrepancy, colnames(X)),
    reconstruction_max_scaled = stats::setNames(scaled, colnames(X)),
    reconstruction_roundoff_threshold = threshold,
    rank_convention = "R stats::lm.fit / dqrdc2 pivoted QR", rank_tol = or_public_rank_tol(),
    control_fit_rank_check = "private compact LAPACK QR relative-diagonal criterion at the recorded lm.fit tolerance; separate from alias removal",
    preprocessing_work_ops_bound = 128 * n * p * p,
    checked_rows = n, row_universe = "all common repeat-observer rows, including controls and treated prediction rows",
    frozen_during_differentiation = TRUE, numerical_certificate = FALSE))
}

# The public default outcome recipe of a frame's static Z_a columns: an
# intercept recipe when the first column is the model-matrix intercept and
# equals one on every A row, otherwise every column is a slope. The same
# statement serves both public fits and the CP08 descriptor-context predicate.
or_public_default_recipe <- function(frame) {
  cols <- colnames(frame$Z_a)
  has_icpt <- length(cols) >= 1L && (cols[1L] %in% c("(Intercept)", "intercept")) &&
    all(frame$Z_a[, 1L] == 1)
  or_pooled_outcome_recipe(cols, has_icpt, if (has_icpt) cols[-1L] else cols)
}

# CP08 descriptor context: only the baseline OR/DR fits read a
# `cdid_cp04_descriptor_context` in the view position, and never with
# reference inputs; every other consumer keeps its link views.
or_public_views_or_context <- function(method, views, reference) {
  if (inherits(views, "cdid_cp04_descriptor_context") && (!identical(method, "baseline") || !is.null(reference)))
    cp04_abort("The CP04 descriptor context serves only the baseline OR/DR route without reference inputs.", "integrity")
  invisible(views)
}

# `reference` (D60) is the sealed recipe-indexed reader; its frozen union
# design replaces the static Z_a columns for every block and link. A dated
# conditional-selection request always uses the general selection-model
# consumer (recipe-keyed cells, pooling inside each recipe).
or_public_fit <- function(method, reader, views, ctrl, retained, work,
                          automatic_chunk_rows = FALSE, outcome_recipe = NULL,
                          selection_model = "pair_logit", selection_pooling = "separate",
                          reference = NULL) {
  external <- max(0, retained - reader$accounting$borrowed_owner_entries)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(outcome_recipe) && !is.null(reference))
    cp04_abort("A reference-input outcome regression needs its frozen outcome recipe.", "integrity")
  if (is.null(outcome_recipe)) outcome_recipe <- or_public_default_recipe(frame)
  or_public_views_or_context(method, views, reference)
  tryCatch(if (identical(method, "baseline"))
    fit_or_baseline(reader, views, ctrl$owner_control, or_public_rank_tol(),
      outcome_recipe = outcome_recipe,
      retained_external_entries = external, external_work_ops = work, reference = reference) else
    if (!is.null(reference)) fit_or_selection_models(reader, views, ctrl$owner_control, ctrl$glm_control,
      or_public_rank_tol(), outcome_recipe = outcome_recipe,
      retained_external_entries = external, external_work_ops = work,
      automatic_chunk_rows = automatic_chunk_rows,
      selection_model = selection_model, selection_pooling = selection_pooling, reference = reference) else
    fit_or_corollary1(reader, views, ctrl$owner_control, ctrl$glm_control, or_public_rank_tol(),
      outcome_recipe = outcome_recipe,
      retained_external_entries = external, external_work_ops = work,
      automatic_chunk_rows = automatic_chunk_rows,
      selection_model = selection_model, selection_pooling = selection_pooling),
    cdid_or_wls_rank_error = function(e) {
      e$message <- paste0(conditionMessage(e), " Control-only support/rank or residual conditioning failure after common alias preprocessing; ",
        "no further terms, samples or moments were removed. This computation failure does not establish population ATT nonidentification.")
      e$or_columns <- if (is.null(reference)) colnames(frame$Z_a) else reference$columns
      e$or_stage <- "control-only weighted regression after frozen common extraction"
      stop(e)
    })
}

dr_public_fit <- function(method, reader, views, ctrl, retained, work,
                          automatic_chunk_rows = FALSE, outcome_recipe = NULL,
                          treatment_columns = NULL,
                          selection_model = "pair_logit", selection_pooling = "separate",
                          reference = NULL) {
  external <- max(0, retained - reader$accounting$borrowed_owner_entries)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(outcome_recipe) && !is.null(reference))
    cp04_abort("A reference-input DR needs its frozen outcome recipe.", "integrity")
  if (is.null(outcome_recipe)) outcome_recipe <- or_public_default_recipe(frame)
  or_public_views_or_context(method, views, reference)
  tryCatch(if (identical(method, "baseline"))
    fit_dr_baseline(reader, views, ctrl$owner_control, ctrl$glm_control, or_public_rank_tol(),
      outcome_recipe = outcome_recipe, treatment_columns = treatment_columns,
      retained_external_entries = external, external_work_ops = work,
      automatic_chunk_rows = automatic_chunk_rows, reference = reference) else
    fit_dr_corollary1(reader, views, ctrl$owner_control, ctrl$glm_control, or_public_rank_tol(),
      outcome_recipe = outcome_recipe, treatment_columns = treatment_columns,
      retained_external_entries = external, external_work_ops = work,
      automatic_chunk_rows = automatic_chunk_rows,
      selection_model = selection_model, selection_pooling = selection_pooling, reference = reference),
    cdid_or_wls_rank_error = function(e) {
      e$message <- paste0(conditionMessage(e), " Control-only support/rank or residual conditioning failure after common alias preprocessing; ",
        "no further terms, samples or moments were removed. This computation failure does not establish population ATT nonidentification.")
      e$or_columns <- outcome_recipe$columns
      e$or_stage <- "control-only weighted regression after frozen common extraction"
      stop(e)
    })
}

chained_full_if_budget <- function(backend, B, links) {
  if (inherits(backend, "cdid_ipw_raw_backend")) return(ipw_raw_if_budget(backend, B, links))
  if (inherits(backend, "cdid_reference_ipw_backend")) return(reference_ipw_if_budget(backend, B, links))
  if (inherits(backend, "cdid_cp04_custom_backend")) return(cp04_custom_if_budget(backend, B, links))
  if (inherits(backend, "cdid_augmented_backend")) return(augmented_if_budget(backend, B, links))
  if (inherits(backend, "cdid_imputation_backend")) return(imputation_if_budget(backend, B, links))
  if (inherits(backend, "cdid_dr_baseline")) return(dr_baseline_if_budget(backend, B, links))
  if (inherits(backend, "cdid_dr_corollary1")) return(dr_corollary1_if_budget(backend, B, links))
  if (inherits(backend, "cdid_or_baseline")) return(or_baseline_if_budget(backend, B, links))
  if (inherits(backend, "cdid_or_corollary1")) return(or_corollary1_if_budget(backend, B, links))
  cp04_if_budget(backend, B, links)
}
# CP08 stage 4: the compact native plan of an eligible complete-support OR/DR
# replay, built once per replay from the validated sealed backend (else NULL).
chained_full_if_plan <- function(backend, links, B) {
  if (inherits(backend, "cdid_ipw_raw_backend")) return(ipw_raw_if_plan(backend, links, B))
  if (inherits(backend, c("cdid_reference_ipw_backend", "cdid_cp04_custom_backend", "cdid_augmented_backend",
      "cdid_imputation_backend")) || !is.list(backend) || is.null(backend$complete_support)) return(NULL)
  if (inherits(backend, "cdid_dr_baseline")) return(dr_baseline_native_plan(backend, links, B))
  if (inherits(backend, "cdid_or_baseline")) return(or_baseline_native_plan(backend, links, B))
  NULL
}
chained_full_if_chunk <- function(backend, rows, links, scale, plan = NULL) {
  if (inherits(backend, "cdid_ipw_raw_backend")) return(evaluate_ipw_raw_if_chunk(backend, rows, links, scale = scale, plan = plan))
  if (inherits(backend, "cdid_reference_ipw_backend")) return(evaluate_reference_ipw_if_chunk(backend, rows, links, scale))
  if (inherits(backend, "cdid_cp04_custom_backend")) return(evaluate_cp04_custom_if_chunk(backend, rows, links, scale))
  if (inherits(backend, "cdid_augmented_backend")) return(evaluate_augmented_if_chunk(backend, rows, links, scale))
  if (inherits(backend, "cdid_imputation_backend")) return(evaluate_imputation_if_chunk(backend, rows, links, scale))
  if (inherits(backend, "cdid_dr_baseline")) return(evaluate_dr_baseline_if_chunk(backend, rows, links, scale, plan = plan))
  if (inherits(backend, "cdid_dr_corollary1")) return(evaluate_dr_corollary1_if_chunk(backend, rows, links, scale))
  if (inherits(backend, "cdid_or_baseline")) return(evaluate_or_baseline_if_chunk(backend, rows, links, scale, plan = plan))
  if (inherits(backend, "cdid_or_corollary1")) return(evaluate_or_corollary1_if_chunk(backend, rows, links, scale))
  evaluate_cp04_if_chunk(backend, rows, links, scale)
}
chained_full_link_estimate <- function(backend, j) {
  if (inherits(backend, "cdid_ipw_raw_backend")) return(ipw_raw_link_estimate(backend, j))
  l <- backend$links[[j]]
  if (inherits(backend, c("cdid_augmented_backend", "cdid_imputation_backend", "cdid_dr_baseline",
                          "cdid_dr_corollary1", "cdid_or_baseline", "cdid_or_corollary1",
                          "cdid_reference_ipw_backend", "cdid_cp04_custom_backend"))) return(l$mu)
  f <- backend$treatments[[l$model]]; unname(f$means["T"] - f$means["M"])
}

or_public_diagnostics <- function(res, control) {
  b <- res$nuisance_registry$cp04
  bound <- 4 * b$accounting$owned_entries + 4096
  cp04_budget(control$owner_control, c(max(c(0, res$resource_accounting$simultaneous_entries,
    res$resource_accounting$initial_envelope_entries)), bound),
    4 * bound, 8 * bound, bound, "OR public diagnostic copy")
  shared <- b$shared[intersect(c("p", "keep", "beta", "coefficients", "fit", "summary",
    "basis", "control_support_ids", "rss", "objective", "verification_scan"), names(b$shared))]
  records <- list(
    selection = lapply(seq_along(b$selections), function(j)
      c(list(id = j, spec = b$selection_specs[[j]]), b$selections[[j]])),
    regression = lapply(seq_along(b$regressions), function(j) {
      reg <- b$regressions[[j]]
      reg_compact <- reg[setdiff(names(reg), "qr")]
      c(list(id = j, role = "local_intercept_block_of_single_joint_fit",
        spec = b$regression_specs[[j]]), reg_compact)
    })
  )
  block_keys <- vapply(b$regression_specs, function(s) {
    if (!is.null(s$block_key)) s$block_key else if (!is.null(s$key)) s$key else ""
  }, character(1))
  local_icpts <- stats::setNames(
    vapply(b$regressions, function(r) as.double(r$alpha), numeric(1)),
    block_keys
  )
  coefficients <- list(
    shared_slopes = b$shared$coefficients,
    local_intercepts = local_icpts,
    has_local_intercepts = isTRUE(b$outcome_recipe$intercept)
  )
  supports_diag <- lapply(seq_along(b$supports), function(j)
    c(list(id = j, spec = b$support_specs[[j]]), b$supports[[j]]))
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st_cohorts <- reader$state$cohorts
  link_ids <- vapply(b$links, `[[`, character(1), "link_id")
  fitted_catalog <- data.frame(
    index = seq_along(b$links),
    link_id = link_ids,
    role = ifelse(startsWith(link_ids, "cp04:auxiliary:"), "auxiliary", "catalog"),
    cohort = st_cohorts[vapply(b$links, `[[`, integer(1), "cohort")],
    s = vapply(b$links, function(l) as.double(b$support_specs[[l$support_id]]$pair$s), numeric(1)),
    t = vapply(b$links, function(l) as.double(b$support_specs[[l$support_id]]$pair$t), numeric(1)),
    duration_steps = vapply(b$links, function(l) as.double(b$regressions[[l$regression]]$L), numeric(1)),
    regression_block = vapply(b$links, `[[`, integer(1), "regression"),
    support_id = vapply(b$links, `[[`, integer(1), "support_id"),
    mu = vapply(b$links, `[[`, numeric(1), "mu"),
    stringsAsFactors = FALSE
  )
  diag_copy <- list(records = records, shared = shared, coefficients = coefficients,
    fitted_catalog = fitted_catalog, supports = supports_diag)
  if (unname(cp04_payload(diag_copy)["entries"]) > bound)
    cp04_abort("OR diagnostic copy bound underestimated.", "bound")
  out <- list(status = "available", method = "reg",
    outcome_model = "shared one-period outcome-change control-mean regression",
    fit_structure = "single_joint_fit_with_block_local_intercepts",
    duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    selection = if (inherits(b, "cdid_or_baseline")) "none" else "corollary1",
    selection_model = if (inherits(b, "cdid_or_baseline")) NULL else if (is.null(b$selection_model)) "pair_logit" else b$selection_model,
    selection_pooling = if (inherits(b, "cdid_or_baseline")) NULL else if (is.null(b$selection_pooling)) "separate" else b$selection_pooling,
    selection_route = res$selection,
    outcome_recipe = b$outcome_recipe,
    fitted_catalog = fitted_catalog,
    coefficients = coefficients,
    basis = list(common_design = if (is.null(b$reference_ref)) res$controls$or_basis else chained_backend_reference(b)$basis,
      pooled_control = b$shared$basis),
    shared = shared,
    records = records,
    supports = supports_diag,
    links = b$links,
    accounting = b$accounting,
    frame_counts = list(n_full_frame_units = res$reporting$inference$n_total,
      n_pair_eligible_A_units = res$frame_owner$n_a_units),
    auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
    interpretation = paste("Single joint weighted least-squares fit of shared one-period outcome-change slopes",
      "across all catalog and auxiliary control blocks with local block intercepts (when the outcome formula includes an intercept)",
      "and stored-calendar transition durations L_e = i_t - i_s; block records describe local intercept blocks of that single joint fit,",
      "not independent regression fits. Complete F-row influence functions with global slope dependencies. No treatment fits or numerical certificate."))
  # D60: recipe table, extraction conventions and provenance of the dated route.
  if (!is.null(b$reference_ref)) out$reference_inputs <- chained_reference_record(chained_backend_reference(b))
  out
}

dr_public_diagnostics <- function(res, control) {
  b <- res$nuisance_registry$cp04
  bound <- 4 * b$accounting$owned_entries + 4096
  cp04_budget(control$owner_control, c(max(c(0, res$resource_accounting$simultaneous_entries,
    res$resource_accounting$initial_envelope_entries)), bound),
    4 * bound, 8 * bound, bound, "DR public diagnostic copy")
  shared <- b$shared[intersect(c("p", "keep", "beta", "coefficients", "fit", "summary",
    "basis", "control_support_ids", "rss", "objective", "verification_scan"), names(b$shared))]
  records <- list(
    selection = lapply(seq_along(b$selections), function(j)
      c(list(id = j, spec = b$selection_specs[[j]]), b$selections[[j]])),
    regression = lapply(seq_along(b$regressions), function(j) {
      reg <- b$regressions[[j]]
      reg_compact <- reg[setdiff(names(reg), "qr")]
      c(list(id = j, role = "local_intercept_block_of_single_joint_fit",
        spec = b$regression_specs[[j]]), reg_compact)
    }),
    treatment = lapply(seq_along(b$treatments), function(j)
      c(list(id = j, spec = b$treatment_specs[[j]]), b$treatments[[j]]))
  )
  block_keys <- vapply(b$regression_specs, function(s) {
    if (!is.null(s$block_key)) s$block_key else if (!is.null(s$key)) s$key else ""
  }, character(1))
  local_icpts <- stats::setNames(
    vapply(b$regressions, function(r) as.double(r$alpha), numeric(1)),
    block_keys
  )
  coefficients <- list(
    shared_slopes = b$shared$coefficients,
    local_intercepts = local_icpts,
    has_local_intercepts = isTRUE(b$outcome_recipe$intercept),
    treatment = lapply(b$treatments, `[[`, "coefficients")
  )
  supports_diag <- lapply(seq_along(b$supports), function(j)
    c(list(id = j, spec = b$support_specs[[j]]), b$supports[[j]]))
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st_cohorts <- reader$state$cohorts
  link_ids <- vapply(b$links, `[[`, character(1), "link_id")
  fitted_catalog <- data.frame(
    index = seq_along(b$links),
    link_id = link_ids,
    role = ifelse(startsWith(link_ids, "cp04:auxiliary:"), "auxiliary", "catalog"),
    cohort = st_cohorts[vapply(b$links, `[[`, integer(1), "cohort")],
    s = vapply(b$links, function(l) as.double(b$support_specs[[l$support_id]]$pair$s), numeric(1)),
    t = vapply(b$links, function(l) as.double(b$support_specs[[l$support_id]]$pair$t), numeric(1)),
    duration_steps = vapply(b$links, function(l) as.double(b$regressions[[l$regression]]$L), numeric(1)),
    regression_block = vapply(b$links, `[[`, integer(1), "regression"),
    treatment_model = vapply(b$links, `[[`, integer(1), "model"),
    support_id = vapply(b$links, `[[`, integer(1), "support_id"),
    mu = vapply(b$links, `[[`, numeric(1), "mu"),
    mu_dr = vapply(b$links, `[[`, numeric(1), "mu_dr"),
    mu_or = vapply(b$links, `[[`, numeric(1), "mu_or"),
    mu_ipw = vapply(b$links, `[[`, numeric(1), "mu_ipw"),
    stringsAsFactors = FALSE
  )
  diag_copy <- list(records = records, shared = shared, coefficients = coefficients,
    fitted_catalog = fitted_catalog, supports = supports_diag)
  if (unname(cp04_payload(diag_copy)["entries"]) > bound)
    cp04_abort("DR diagnostic copy bound underestimated.", "bound")
  out <- list(status = "available", method = "dr",
    outcome_model = "shared one-period outcome-change control-mean regression",
    treatment_model = "risk-set treatment propensity logit",
    fit_structure = "single_joint_or_fit_with_block_local_intercepts_plus_deduplicated_treatment_logits",
    duration_convention = "stored-calendar transition count L_e = i_t - i_s",
    selection = if (inherits(b, "cdid_dr_baseline")) "none" else "corollary1",
    selection_model = if (inherits(b, "cdid_dr_baseline")) NULL else if (is.null(b$selection_model)) "pair_logit" else b$selection_model,
    selection_pooling = if (inherits(b, "cdid_dr_baseline")) NULL else if (is.null(b$selection_pooling)) "separate" else b$selection_pooling,
    selection_route = res$selection,
    outcome_recipe = b$outcome_recipe,
    treatment_recipe = b$treatment_recipe,
    fitted_catalog = fitted_catalog,
    coefficients = coefficients,
    basis = list(common_design = if (is.null(b$reference_ref)) res$controls$or_basis else chained_backend_reference(b)$basis,
      pooled_control = b$shared$basis),
    shared = shared,
    records = records,
    supports = supports_diag,
    links = b$links,
    accounting = b$accounting,
    frame_counts = list(n_full_frame_units = res$reporting$inference$n_total,
      n_pair_eligible_A_units = res$frame_owner$n_a_units),
    auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
    interpretation = paste("Treatment-DR extension of the shared one-period outcome-change control-mean regression",
      "and residualized Hajek IPW over deduplicated risk-set treatment propensity logits under the declared observation route;",
      "complete F-row influence functions with global OR slope, treatment propensity, and observation/selection paths.",
      "No numerical certificate."))
  if (!is.null(b$reference_ref)) out$reference_inputs <- chained_reference_record(chained_backend_reference(b))
  out
}
