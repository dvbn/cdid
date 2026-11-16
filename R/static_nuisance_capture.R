# The static backends already hold the fitted coefficients and their unit scores.
# The optional public report keeps only tables and consumes scores in row chunks.

static_report_equation <- function(lhs, slopes, scale = "X_ref", local = "alpha[cell]") {
  design <- if (length(slopes)) paste0(" + sum_j beta[j] * ", scale, "[j]") else ""
  paste0(lhs, " = ", local, design)
}

static_report_block <- function(id, keys, alpha, beta, slopes, lhs, population,
                                scale = "X_ref", converged = TRUE, iterations = NA_integer_) {
  list(id = id, keys = keys, estimate = c(alpha, beta), slopes = slopes,
    equation = static_report_equation(lhs, slopes, scale),
    population = population, converged = converged, iterations = iterations)
}

static_report_covariance <- function(report, evaluate, group, cluster, n, chunk_rows,
                                     compute_if, workspace_entries, work_left, source,
                                     accounts_for, producer_width, work_per_row, factor_entries,
                                     frame_matches) {
  width <- sum(vapply(report, function(b) as.double(length(b$estimate)), numeric(1)))
  groups <- cluster$n_clusters
  available <- compute_if
  code <- if (available) "ok" else "point_only"
  reason <- if (available) NULL else "Coefficient inference is disabled by compute_inffunc = FALSE."
  # Bound additional scratch and work. Retain coefficients if the optional
  # influence pass does not fit the requested resource limits.
  accumulator <- if (cluster$clustered) 4 * as.double(groups) * width + n else 8 * width
  match_entries <- if (frame_matches > 0 || cluster$clustered) 8 * as.double(n) else 0
  room <- floor((workspace_entries - accumulator - factor_entries - match_entries) / producer_width)
  chunk <- min(chunk_rows, max(1, room), n)
  # match() can rescan or hash the full frame on every call, even in a
  # one-row batch. The budget uses the final admitted chunk size.
  work <- as.double(n) * (work_per_row + 16 * frame_matches * ceiling(n / chunk) +
    if (cluster$clustered) 16 else 0)
  if (available && (room < 1 || work > work_left)) {
    available <- FALSE
    code <- "budget"
    reason <- "The optional coefficient influence pass exceeds the declared workspace or work limit."
  }
  variance <- rep(NA_real_, width)
  if (available) {
    acc <- chained_nuisance_moments_new(width, if (is.null(group)) NULL else groups)
    for (first in seq.int(1, n, by = chunk)) {
      rows <- seq.int(first, min(n, first + chunk - 1))
      acc <- chained_nuisance_moments_add(acc, evaluate(rows), group[rows])
    }
    variance <- chained_nuisance_moments_variance(acc, n, groups)
  }
  list(variance = variance, work = if (available) work else 0, record = list(
    status = if (available) "available" else "unavailable", code = code, reason = reason,
    type = if (cluster$clustered) "cluster" else "unit", cluster_var = cluster$clustername,
    n_groups = groups, N = n, diagonal_only = TRUE,
    estimator = paste(
      "Centered sandwich with one coefficient influence contribution per original unit.",
      "Sum within clusters when requested, center across all frame units or clusters,",
      "and divide the sum of squares by N^2, with no small-sample factor."
    ),
    source = source, accounts_for = accounts_for))
}

# The selection and other role views use the same coefficient schema. Use the
# fitted cell order for intercepts and the retained design-column order for slopes.
static_report_tables <- function(report, detail, comparisons, variance) {
  ids <- vapply(report, `[[`, "", "id")
  sizes <- vapply(report, function(b) length(b$estimate), integer(1L))
  var_by_block <- split(variance, rep(seq_along(report), sizes))
  coefs <- lapply(seq_along(report), function(m) {
    b <- report[[m]]
    n_cell <- length(b$keys)
    cells <- match(b$keys, detail$key)
    intercept <- if (n_cell == 1L) "(Intercept)" else paste0("(Intercept)[", detail$label[cells], "]")
    chained_nuisance_coef_frame(b$id, c(intercept, b$slopes),
      c(rep("local_intercept", n_cell),
        rep(if (n_cell == 1L) "slope" else "shared_slope", length(b$slopes))),
      c(detail$cell[cells], rep(NA_integer_, length(b$slopes))), b$estimate, var_by_block[[m]])
  })
  blocks <- do.call(rbind, lapply(report, function(b) {
    cells <- match(b$keys, detail$key)
    rows <- sum(detail$n_fitting_rows[cells])
    uses <- if (!is.null(comparisons$block)) sum(comparisons$block == b$id) else
      sum(vapply(strsplit(comparisons$fits, ",", fixed = TRUE),
        function(x) b$id %in% x, logical(1L)))
    data.frame(block = b$id, label = paste0(b$id, " [", length(cells), " cell(s)]"),
      population = b$population, n_cells = length(cells), n_parameters = length(b$estimate),
      n_fitting_rows = rows, n_units = b$n_units, weight_sum = rows,
      kept_slopes = paste(b$slopes, collapse = ", "), aliased_slopes = "",
      n_comparisons = uses, equation = b$equation,
      branch = if (length(b$slopes)) "intercept_slope" else "intercept_only",
      converged = b$converged, iterations = b$iterations, rank = length(b$estimate),
      n_fitting_clusters = NA_integer_, n_contributing_clusters = NA_integer_)
  }))
  detail$block <- rep(ids, lengths(lapply(report, `[[`, "keys")))[
    match(detail$key, unlist(lapply(report, `[[`, "keys"), use.names = FALSE))]
  list(blocks = blocks, coefficients = do.call(rbind, coefs), detail = detail)
}

static_report_transition_detail <- function(cells, keys, field) {
  cells <- cells[keys]
  dates <- vapply(cells, `[[`, numeric(1L), "t_num")
  prior <- vapply(cells, `[[`, integer(1L), "d")
  rows <- vapply(cells, `[[`, integer(1L), field)
  data.frame(cell = seq_along(keys), key = keys, label = paste0("t=", dates, ", prior D=", prior),
    s = vapply(cells, `[[`, numeric(1L), "t_prev_num"), t = dates,
    prior_state = prior, n_fitting_rows = rows)
}

static_report_role <- function(role, report, detail, comparisons, variance, request, n,
                               response, template, sample, use, pooling) {
  comparisons$block <- comparisons$fits
  tables <- static_report_tables(report, detail, comparisons, variance$variance)
  slopes <- unique(unlist(lapply(report, `[[`, "slopes"), use.names = FALSE))
  label <- switch(role, treatment = "switcher propensity", outcome = "counterfactual stayer outcome change",
    imputation = "actual-change imputation")
  out <- list(status = "available", role = label, route = request$route, pooling = pooling,
    formula = paste(deparse(request$xformla), collapse = " "),
    template = template, response = response, sample = sample, use = use,
    design_columns = c("(Intercept)", slopes), intercept_term = "(Intercept)", slope_terms = slopes,
    weights = list(text = "Unit weights"), n_frame = n, blocks = tables$blocks,
    coefficients = tables$coefficients, comparisons = comparisons, covariance = variance$record)
  out[[switch(role, treatment = "models", outcome = "control_blocks", imputation = "cells")]] <- tables$detail
  out
}

static_report_treatment <- function(fit, result, variance_for) {
  models <- fit$nuisance$logit
  keys <- names(models)
  detail <- static_report_transition_detail(fit$all_cells, keys, "n_risk")
  report <- lapply(seq_along(models), function(j) {
    m <- models[[j]]
    b <- static_report_block(keys[j], keys[j], m$gamma_hat[1L], m$gamma_hat[-1L],
      names(m$gamma_hat)[-1L], "logit Pr(switcher | at risk, X_ref)",
      paste("Transition", detail$label[j], "among switchers and same-prior-state stayers"),
      converged = m$diagnostics$converged, iterations = m$diagnostics$iter)
    b$n_units <- detail$n_fitting_rows[j]
    b
  })
  width <- sum(lengths(lapply(models, `[[`, "gamma_hat")))
  v <- variance_for(report, function(rows) do.call(cbind, lapply(models, function(m) {
    m$ell_gamma[rows, , drop = FALSE]
  })), "Static treatment logit's retained ell_gamma on the original units.",
    "Accounts for the logit score on the full transition risk set.", 10 * width, 8 * width)
  comparisons <- data.frame(comparison = keys, s = detail$s, t = detail$t,
    prior_state = detail$prior_state, fits = keys)
  static_report_role("treatment", report, detail, comparisons, v, result$request, result$n,
    "Indicator of switching treatment state at the transition date",
    "logit Pr(switcher | at risk, X_ref) = alpha[transition] + X_ref' beta[transition]",
    "Switchers and same-prior-state stayers", "Treatment odds for reweighting the stayer comparisons.",
    "transition risk set")
}

static_report_outcome <- function(fit, result, variance_for) {
  model <- fit$nuisance$or
  partition <- model$partition
  keys <- names(model$alpha_hat)
  detail <- static_report_transition_detail(fit$all_cells, keys, "n_stayers")
  reconstructed <- inherits(model, "cdid_static_reconstructed_or_fit")
  slopes <- if (reconstructed) model$b_hat else model$slope_hat
  lhs <- if (reconstructed) "Reconstructed mean stayer change given X_ref" else
    "Mean adjacent stayer change given X_ref"
  report <- lapply(partition$block_ids, function(id) {
    cids <- partition$blocks[[id]]
    b <- static_report_block(id, cids, model$alpha_hat[cids], slopes[[id]],
      partition$block_colnames[[id]], lhs,
      "Same-prior-state stayers in the listed transition cells")
    b$n_units <- NA_integer_
    b
  })
  evaluate <- function(rows) {
    if (reconstructed) {
      chunk <- fit$nuisance$observation$with_if_batch(fit$unit_id[rows], function() {
        model$if_callback(rows, include_chi_A_tensor = FALSE)
      })
      alpha <- chunk$chi_alpha
      beta <- chunk$chi_b
    } else {
      alpha <- model$ell_alpha[rows, , drop = FALSE]
      beta <- lapply(model$block_ell_slope, function(m) m[rows, , drop = FALSE])
    }
    do.call(cbind, lapply(report, function(b) cbind(alpha[, b$keys, drop = FALSE], beta[[b$id]])))
  }
  width <- length(keys) + sum(lengths(slopes))
  p <- as.double(lengths(slopes))
  cell_widths <- vapply(unname(partition$cell_block[keys]),
    function(id) as.double(length(slopes[[id]])), numeric(1))
  p_cells <- sum(cell_widths)
  producer_width <- 16 * (width + p_cells)
  work <- 32 * (width + p_cells)
  factor_entries <- 0
  frame_matches <- 0
  if (reconstructed) {
    observation <- fit$nuisance$observation
    selection_cost <- static_report_observation_score_cost(observation$selection_fit)
    imputation_cost <- static_report_observation_score_cost(observation$imputation_fit)
    obs_width <- length(observation$selection_fit$coefficients) + length(observation$imputation_fit$coefficients)
    fits_with_scores <- as.integer(!is.null(observation$selection_fit)) +
      as.integer(!is.null(observation$imputation_fit))
    n_edges <- sum(lengths(lapply(model$cells, `[[`, "varpi_c")))
    producer_width <- producer_width + 16 * obs_width
    work <- work + 32 * (n_edges * (1 + obs_width) * (1 + max(p)) + 2 * sum(p^2) + sum(cell_widths^2)) +
      selection_cost$work + imputation_cost$work
    factor_entries <- 4 * max(p^2) + selection_cost$factor_entries + imputation_cost$factor_entries
    # One batch setup, at most one score lookup per active observation fit,
    # and the mean and covariance evaluators for each edge. Each evaluator
    # may match positions and IDs in its wrapper and its observation inputs.
    frame_matches <- 1 + fits_with_scores + 2 * n_edges * (4 + fits_with_scores)
  }
  v <- variance_for(report, evaluate,
    if (reconstructed) "Static reconstructed OR callback: chi_alpha and chi_b." else
      "Static complete-data OR's retained ell_alpha and block_ell_slope.",
    if (reconstructed) paste(
      "Accounts for original-unit stayer shares, covariate means, Gram matrices and reconstructed",
      "edge moments, including their fitted observation or imputation channels when active."
    ) else "Accounts for the original-unit stayer score across all training cells of the block.",
    producer_width, work, factor_entries, frame_matches)
  target <- detail[match(fit$cell_ids, detail$key), , drop = FALSE]
  comparisons <- data.frame(comparison = target$key, s = target$s, t = target$t,
    prior_state = target$prior_state, fits = unname(partition$cell_block[target$key]))
  component <- static_report_role("outcome", report, detail, comparisons, v, result$request, result$n,
    if (reconstructed) "Adjacent stayer change reconstructed from observation-graph edge moments" else
      "Observed adjacent outcome change of stayers",
    static_report_equation(lhs, unique(unlist(partition$block_colnames)), local = "alpha[transition]"),
    "Same-prior-state stayers in the block's transition cells (including training-only cells, when used)",
    "Counterfactual stayer change evaluated at the comparison's covariates.",
    result$request$outcome_pooling_mode)
  if (reconstructed) {
    component$weights$text <- "Unit weights in the stayer covariate Gram, with the fitted route's observation correction in the edge moments"
  }
  component
}

# Each score replay can transpose the largest slope factor and perform two
# triangular solves per active block. During reconstruction, reserve space for
# observation-factor scratch alongside the enclosing OR callback.
static_report_observation_score_cost <- function(model) {
  blocks <- model$layout$blocks
  p <- vapply(blocks, function(bl) as.double(bl$d_s), 0)
  cells <- lengths(lapply(blocks, `[[`, "cells"))
  list(work = 32 * sum(cells * (1 + p) + 2 * p^2), factor_entries = 4 * max(0, p^2))
}

static_report_observation_detail <- function(model, adapter, result, is_date, is_logit) {
  date_value <- function(keys) result$dates_num[as.integer(substring(keys, 2L))]
  groups <- adapter$groups
  map <- lapply(groups, function(group) {
    pairs <- group$pairs
    cells <- if (is_date) unique(unlist(lapply(pairs, `[[`, "date_cell_keys"), use.names = FALSE)) else
      vapply(pairs, `[[`, "", "cell_key")
    pair_start <- if (is_date) rep(NA_real_, length(cells)) else date_value(vapply(pairs, `[[`, "", "start"))
    end <- if (is_date) date_value(sub(".*::", "", cells)) else date_value(vapply(pairs, `[[`, "", "end"))
    role <- if (group$role == "G") "switcher" else "stayer"
    label <- paste("transition t =", group$target, ", prior D =", group$prior_state, ",", role,
      if (is_date) paste("date", end) else paste0("pair (", pair_start, ", ", end, ")"))
    data.frame(key = cells, label = label, group = group$group, role = role, transition = group$target,
      prior_state = group$prior_state, s = pair_start, t = end)
  })
  detail <- do.call(rbind, map)
  rownames(detail) <- NULL
  detail$cell <- seq_len(nrow(detail))
  # Reuse the fitted observation rates and WLS row counts. Probability
  # extrema would require replaying the fitting designs and are not retained.
  if (is_logit) {
    detail$n_fitting_rows <- vapply(groups[detail$group], function(g) {
      if (is.null(g$n_fit)) g$n_role else g$n_fit
    }, 0L)
    rate <- vapply(detail$key, function(key) {
      model$block_fits[[model$layout$cell_block[[key]]]]$pi[[model$layout$cell_local[[key]]]]
    }, 0)
    detail$n_observed <- round(detail$n_fitting_rows * rate)
    detail$p_min <- detail$p_max <- NA_real_
    detail$n_floor <- NA_integer_
  } else {
    detail$n_fitting_rows <- vapply(detail$key, function(key) {
      model$block_fits[[model$layout$cell_block[[key]]]]$blocks[[key]]$n
    }, 0)
    detail$n_observed <- detail$n_fitting_rows
  }
  detail
}

static_report_observation <- function(role, adapter, result, variance_for) {
  logit <- identical(role, "selection")
  model <- if (logit) adapter$selection_fit else adapter$imputation_fit
  dates <- logit && identical(adapter$selection_model, "independent_period_logit")
  transport <- adapter$selection_transport
  universe <- if (is.null(transport)) "F" else transport$universe
  sample <- if (logit && universe != "F") paste("Transition-role units in U =", universe) else
    if (logit) "Full switcher or stayer transition-role populations" else
      "Observed outcome pairs in the switcher or stayer transition role"
  lhs <- if (dates) "logit Pr(Y_date observed | transition role, X_ref)" else
    if (logit) "logit Pr(both outcomes observed | transition role, U, X_ref)" else
      "Mean actual outcome change given transition role and X_ref"
  scale <- if (!logit && identical(adapter$imputation_scale, "per_period"))
    "(u_idx - s_idx) * X_ref" else "X_ref"
  detail <- static_report_observation_detail(model, adapter, result, dates, logit)
  report <- lapply(model$layout$blocks, function(bl) {
    fitted <- model$block_fits[[bl$block]]
    b <- static_report_block(bl$block, bl$cells, fitted$alpha, fitted$beta, bl$slope_names, lhs, sample,
      scale = scale, converged = if (logit) fitted$diagnostics$converged else TRUE,
      iterations = if (logit) fitted$diagnostics$iter else NA_integer_)
    cells <- detail[match(bl$cells, detail$key), , drop = FALSE]
    # At a single transition date, distinct roles have disjoint populations.
    # Across dates or across imputation pair masks, the union needs a row scan.
    b$n_units <- if (logit && length(unique(cells$transition)) == 1L) {
      sum(cells$n_fitting_rows[!duplicated(cells$group)])
    } else if (nrow(cells) == 1L) cells$n_fitting_rows else NA_integer_
    b
  })
  pair_df <- do.call(rbind, lapply(adapter$groups, function(group) {
    do.call(rbind, lapply(group$pairs, function(pair) {
      keys <- if (dates) unname(pair$date_cell_keys) else pair$cell_key
      blocks <- unique(unname(model$layout$cell_block[keys]))
      endpoints <- result$dates_num[as.integer(substring(c(pair$start, pair$end), 2L))]
      data.frame(comparison = pair$cell_key, role = if (group$role == "G") "switcher" else "stayer",
        transition = group$target, prior_state = group$prior_state, s = endpoints[1L], t = endpoints[2L],
        index_gap = pair$gap, fits = paste(blocks, collapse = ","), blocks = paste(blocks, collapse = ","))
    }))
  }))
  rownames(pair_df) <- NULL
  width <- length(model$coefficients)
  cost <- static_report_observation_score_cost(model)
  v <- variance_for(report, function(rows) static_nuisance_if(model, model$unit_id[rows]),
    "Static observation adapter's own static_nuisance_if.",
    if (logit) "Accounts for the original-unit observation logit score across cells of the block." else
      "Accounts for the original-unit observed-pair regression score across cells of the block.",
    12 * width, cost$work, cost$factor_entries, frame_matches = 1)
  if (!logit) {
    slopes <- unique(unlist(lapply(report, `[[`, "slopes"), use.names = FALSE))
    component <- static_report_role("imputation", report, detail, pair_df, v, result$request, result$n,
      "Observed actual outcome change Y_u - Y_s", static_report_equation(lhs, slopes, scale),
      sample, "Prediction of actual outcome change for observation correction.", adapter$imputation_pooling)
    component$imputation_scale <- adapter$imputation_scale
    return(component)
  }
  tables <- static_report_tables(report, detail, pair_df, v$variance)
  slopes <- unique(unlist(lapply(report, `[[`, "slopes"), use.names = FALSE))
  list(status = "available", role = "observation-selection propensity", route = result$request$route,
    selection_model = adapter$selection_model, selection_pooling = adapter$selection_pooling,
    link = "logit", response = if (dates) "Outcome observed at the date" else "Both outcomes of the pair observed",
    template = static_report_equation(lhs, slopes),
    derived = if (dates) "The pair probability is the product of the endpoint probabilities under conditional date independence.",
    formula = paste(deparse(result$request$selection_xformla), collapse = " "),
    design_columns = c("(Intercept)", slopes), intercept_term = "(Intercept)", slope_terms = slopes,
    universe = universe, universe_text = sample, population = sample, transport = if (is.null(transport)) NULL else
      list(universe = universe, ratio_column = result$request$selection_rho_column, justification = transport$justification),
    weights = list(text = "Unit weights in the logit"), n_frame = result$n,
    use = if (!is.null(transport)) "Full-role moments use inverse observation weights based on q * rho. The role fitting rate cancels." else
      "The fitted observation probability enters the observation correction for its transition role.",
    glm_control = result$fit$glm_control, probability_floor = ipw_clip_eps,
    blocks = tables$blocks, cells = tables$detail, coefficients = tables$coefficients,
    pair_cohorts = pair_df, comparisons = pair_df, covariance = v$record)
}

static_nuisance_capture <- function(result, cluster, chunk_rows, workspace_entries, work_left) {
  fit <- result$fit
  request <- result$request
  group <- NULL
  variance_for <- function(report, evaluate, source, accounts_for, producer_width, work_per_row,
                           factor_entries = 0, frame_matches = 0) {
    # Make the common cluster map only if a coefficient pass is admitted.
    v <- static_report_covariance(report, evaluate, {
      if (cluster$clustered && is.null(group)) {
        group <<- match(cluster$cluster_vec, unique(cluster$cluster_vec))
      }
      group
    }, cluster, result$n, chunk_rows,
      request$compute_inffunc, workspace_entries, work_left, source, accounts_for,
      producer_width, work_per_row, factor_entries, frame_matches)
    work_left <<- work_left - v$work
    v
  }
  unused <- function(used, code, reason) {
    if (used) chained_nuisance_status("unavailable", "no_fitted_targets",
      "The declared transition cells have no available fitted target.") else
      chained_nuisance_status("not_applicable", code, reason)
  }
  obs <- fit$nuisance$observation
  has_q <- request$route %in% c("weighting", "augmented")
  has_j <- request$route %in% c("imputation", "augmented")
  selection <- if (!is.null(obs$selection_fit)) static_report_observation("selection", obs, result, variance_for) else
    unused(has_q, "no_selection_logit", "This route fits no observation-selection logit.")
  outcome <- if (!is.null(fit$nuisance$or)) static_report_outcome(fit, result, variance_for) else
    unused(request$est_method != "ipw", "no_outcome_regression", "Treatment IPW fits no counterfactual outcome regression.")
  treatment <- if (length(fit$nuisance$logit)) static_report_treatment(fit, result, variance_for) else
    unused(request$est_method != "reg", "no_treatment_propensity", "Treatment regression fits no treatment propensity.")
  imputation <- if (!is.null(obs$imputation_fit)) static_report_observation("imputation", obs, result, variance_for) else
    unused(has_j, "no_imputation", "This route fits no actual-change imputation model.")
  chained_nuisance_capture_assemble(selection, request, request$est_method,
    roles = list(outcome = outcome, treatment = treatment, imputation = imputation))
}
