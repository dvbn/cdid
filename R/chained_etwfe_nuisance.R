# Optional coefficient reporting for the cohort-time preset. Report builders
# receive the fitted roots and factors while they are still live. Only tables
# return to the public fit; row scores, data and model closures stay private.
etwfe_nuisance_capture <- function(est_method, selection_est_method, selection, observed, fitted) {
  components <- list(
    selection = if (selection == "none") {
      chained_nuisance_status("not_applicable", "no_selection_model",
        'selection = "none" fits no observation model.')
    } else if (selection_est_method == "imputation") {
      chained_nuisance_status("not_applicable", "pure_imputation",
        'selection_est_method = "imputation" fits no observation logit.')
    } else etwfe_nuisance_supported(observed$selection, "observation logit"),
    outcome = if (est_method == "ipw") {
      chained_nuisance_status("not_applicable", "no_outcome_regression",
        'est_method = "ipw" fits no outcome regression.')
    } else etwfe_nuisance_supported(fitted$outcome, "outcome regression"),
    treatment = if (est_method == "reg") {
      chained_nuisance_status("not_applicable", "no_treatment_propensity",
        'est_method = "reg" fits no treatment propensity.')
    } else etwfe_nuisance_supported(fitted$treatment, "treatment propensity"),
    imputation = if (selection == "none" || selection_est_method == "ipw") {
      chained_nuisance_status("not_applicable", "no_imputation",
        "This observation-correction choice fits no actual-change imputation model.")
    } else etwfe_nuisance_supported(observed$imputation, "actual-change imputation model")
  )
  structure(list(schema = chained_nuisance_schema, status = "captured", entry = "att_gt",
    route = "cohort_time_fe", est_method = est_method,
    selection_est_method = selection_est_method, selection = selection, components = components,
    source = "Coefficients and score operators of the fitted cohort-time preset."),
    class = "cdid_nuisance_capture")
}

etwfe_nuisance_supported <- function(component, label) {
  if (!is.null(component)) return(component)
  chained_nuisance_status("unavailable", "no_supported_model",
    paste0("No supported equation required a fitted ", label, " in this call."))
}

etwfe_nuisance_stats <- function(width, n, report, score, spare, row_entries, source, blocks = rep.int(1L, width)) {
  clusters <- report$clusters
  group <- if (!is.null(clusters)) match(clusters, unique(clusters))
  G <- if (is.null(group)) n else max(group)
  models <- sort(unique(blocks))
  n_models <- length(models)
  hit <- NULL
  size <- min(n, report$chunk_rows, floor((spare - G * (3 * width + n_models) * !is.null(group) - 8 * width) /
    max(1, row_entries + 4 * width)))
  status <- "ok"
  reason <- NULL
  if (!width) {
    status <- "fixed_probability"
    reason <- "The specified model has no fitted coefficients."
  } else if (!report$compute_inffunc) {
    status <- "point_only"
    reason <- "The fit was requested with compute_inffunc = FALSE."
  } else if (G < 2L) {
    status <- "insufficient_clusters"
    reason <- "At least two sampling units or requested clusters are required for a sandwich SE."
  } else if (size < 1L) {
    status <- "workspace"
    reason <- "The coefficient score pass exceeds the available numeric workspace."
  }

  variance <- rep(NA_real_, width)
  if (identical(status, "ok")) {
    acc <- chained_nuisance_moments_new(width, if (!is.null(group)) G)
    if (!is.null(group)) hit <- matrix(FALSE, G, n_models)
    for (first in seq.int(1L, n, by = size)) {
      rows <- seq.int(first, min(n, first + size - 1L))
      xi <- score(rows)
      group_rows <- if (!is.null(group)) group[rows]
      acc <- chained_nuisance_moments_add(acc, xi, group_rows)
      if (!is.null(group)) {
        for (j in seq_along(models)) {
          local <- rowSums(xi[, blocks == models[j], drop = FALSE] != 0) > 0
          hit[group_rows[local], j] <- TRUE
        }
      }
    }
    variance <- chained_nuisance_moments_variance(acc, n, G)
  }
  count <- if (is.null(hit)) rep(NA_integer_, n_models) else as.integer(colSums(hit))
  names(count) <- models
  covariance <- chained_nuisance_covariance_record(identical(status, "ok"), status, reason,
    !is.null(group), report$cluster_var, G, n, source)
  covariance$estimator <- "Centered original-unit sandwich, summed within requested clusters, divided by N^2, with no small-sample factor."
  covariance$accounts_for <- source
  list(variance = variance, covariance = covariance, clusters = count)
}

etwfe_nuisance_block_covariance <- function(stat, blocks) {
  count <- stat$clusters[match(as.character(unique(blocks)), names(stat$clusters))]
  covariance <- stat$covariance
  if (!length(blocks)) {
    covariance$status <- "not_applicable"
    covariance$code <- "no_free_coefficients"
    covariance$reason <- "The specified model has no fitted coefficients."
  } else if (any(!is.na(count) & count < 2L)) {
    covariance$status <- if (all(count < 2L)) "unavailable" else "partial"
    covariance$code <- "insufficient_clusters"
    covariance$reason <- "At least two clusters must contain nonzero original-unit coefficient contributions."
    n <- stat$clusters[as.character(blocks)]
    stat$variance[!is.na(n) & n < 2L] <- NA_real_
  }
  stat$covariance <- covariance
  stat
}

etwfe_nuisance_cohort <- function(g, report) c(0, report$calendar)[g + 1L]

etwfe_nuisance_formula <- function(formula, columns) {
  if (is.null(formula)) paste(c("Design columns", columns), collapse = " ") else
    paste(deparse(formula), collapse = " ")
}

etwfe_nuisance_role <- function(role, template, response, sample, use, pooling, blocks, coefs,
                                detail, comparisons, covariance, columns, formula) {
  result <- list(status = "available", route = "cohort_time_fe", role = role, template = template,
    response = response, sample = sample, weights = list(text = "unit weights"), use = use,
    pooling = pooling, blocks = blocks, coefficients = coefs, comparisons = comparisons,
    covariance = covariance, design_columns = columns, formula = formula, n_frame = covariance$N)
  result[[names(detail)]] <- detail[[1L]]
  result
}

etwfe_nuisance_comparisons <- function(g, s, t, donor, fit, report) {
  g <- etwfe_nuisance_cohort(g, report)
  donor <- etwfe_nuisance_cohort(donor, report)
  s <- report$calendar[s]
  t <- report$calendar[t]
  data.frame(comparison = paste(g, s, t, donor, sep = ":"), g = g, s = s, t = t,
    controls = ifelse(is.na(donor), "", chained_nuisance_num(donor)), fits = fit)
}

etwfe_nuisance_treatment <- function(fits, g, k, cohort, Z, links, report, spare) {
  if (!length(fits)) return(NULL)
  n <- length(cohort)
  terms <- c("(Intercept)", colnames(Z))
  width <- ncol(Z) + 1L
  if (length(terms) != width) terms <- c("(Intercept)", paste0("z", seq_len(width - 1L)))
  ids <- paste0("treatment_", seq_along(fits))
  cols <- ifelse(is.na(links$donor), 0L, match(
    paste(links$g, links$donor), paste(g, k)))
  use <- which(cols > 0L)
  comparisons <- etwfe_nuisance_comparisons(links$g[use], links$s[use], links$t[use],
    links$donor[use], ids[cols[use]], report)
  tg <- etwfe_nuisance_cohort(g, report)
  dk <- etwfe_nuisance_cohort(k, report)
  cohorts <- paste(tg, dk, sep = ",")
  fitted_n <- tabulate(match(cohort, sort(unique(cohort))))
  n_by <- stats::setNames(fitted_n, sort(unique(cohort)))
  rows <- unname(as.double(n_by[as.character(g)] + n_by[as.character(k)]))
  source <- "Uses the full two-cohort frame logit score and its fitted information, on original units."
  statistics <- lapply(seq_along(fits), function(m) {
    fit <- fits[[m]]
    score <- function(rows) etwfe_tm_odds_score(fit, Z, cohort, g[m], k[m], rows)
    stat <- etwfe_nuisance_stats(width, n, report, score, spare, 8 * width + 10, source)
    cv <- etwfe_nuisance_block_covariance(stat, rep.int(1L, width))
    stat$variance <- cv$variance
    stat$covariance <- cv$covariance
    stat
  })
  coefs <- do.call(rbind, lapply(seq_along(fits), function(m) {
    chained_nuisance_coef_frame(ids[m], terms, c("intercept", rep("slope", width - 1L)),
      rep(NA_integer_, width), fits[[m]]$coefficients, statistics[[m]]$variance)
  }))
  covariance <- statistics[[1L]]$covariance
  status <- vapply(statistics, function(s) s$covariance$status, "")
  if (any(status == "available") && !all(status == "available")) {
    covariance <- statistics[[which(status == "available")[1L]]]$covariance
    covariance$status <- "partial"
    covariance$code <- "insufficient_clusters"
    covariance$reason <- "Some treatment logits were fitted in fewer than two clusters."
  }
  blocks <- data.frame(block = ids, label = paste0("target ", tg, " vs donor ", dk), cohorts = cohorts,
    n_parameters = width, n_fitting_rows = rows, weight_sum = rows,
    n_comparisons = tabulate(cols[use], length(ids)), kept_slopes = paste(terms[-1L], collapse = ","),
    aliased_slopes = "", n_contributing_clusters = vapply(statistics, function(x) unname(x$clusters[1L]), integer(1)),
    equation = paste0("logit Pr(G_i = ", tg, " | Z_i, G_i in {", cohorts, "}) = ",
      paste(c("a", paste0("b[", terms[-1L], "] * ", terms[-1L])[seq_len(width - 1L)]), collapse = " + ")))
  etwfe_nuisance_role("two-cohort treatment propensity", "logit Pr(G_i = g | Z_i, G_i in {g,k}) = a[g,k] + Z_i'b[g,k]",
    "D_i = 1{G_i = g}", "All frame units in each target-donor risk set",
    "Supplies normalized donor odds to the treatment moments.", "two-cohort risk set", blocks, coefs,
    list(models = data.frame(block = ids, g = tg, donor = dk)), comparisons, covariance, terms,
    etwfe_nuisance_formula(report$outcome_formula, terms))
}

etwfe_tm_odds_score <- function(fit, Z, cohort, g, donor, rows) {
  result <- matrix(0, length(rows), ncol(Z) + 1L)
  loc <- which(cohort[rows] %in% c(g, donor))
  if (length(loc)) {
    X <- etwfe_tm_design(Z, rows[loc])
    probability <- ipw_logit_values(X, fit$coefficients)$prob
    result[loc, ] <- t(cp04_solve(fit$H_chol, t(X * (as.double(cohort[rows[loc]] == g) - probability))))
  }
  result
}

etwfe_nuisance_slot_codes <- function(slot, n_t, groups) {
  data.frame(slot = seq_along(slot), date = (slot - 1L) %% n_t + 1L,
    cohort = groups[(slot - 1L) %/% n_t + 1L])
}

etwfe_nuisance_outcome <- function(beta, keep, coef_index, slots, Z, center, amplitude, n,
                                  report, score, spare, row_entries, fitting_rows, model,
                                  sample, source, comparisons = NULL) {
  if (!length(keep)) return(NULL)
  all_slots <- unique(c(slots$slot[slots$gauge], match(keep, coef_index[, 1L], nomatch = 0L)))
  all_slots <- sort(all_slots[all_slots > 0L])
  coord <- as.vector(t(coef_index[all_slots, , drop = FALSE]))
  dates <- report$calendar[slots$date[all_slots]]
  cohorts <- etwfe_nuisance_cohort(slots$cohort[all_slots], report)
  nm <- colnames(Z)
  if (is.null(nm) && length(center)) nm <- paste0("z", seq_along(center))
  features <- c("(Intercept)", nm)
  gauge <- slots$gauge[all_slots]
  term <- paste0(rep(ifelse(cohorts == 0, paste0("date[", dates, "]"),
    paste0("cohort_date[", cohorts, ",", dates, "]")), each = length(features)), ":",
    rep(features, length(all_slots)))
  if (length(center)) source <- paste(source, "The reported centering and scaling define fixed coordinates for these coefficient SEs.")
  stat <- etwfe_nuisance_stats(length(keep), n, report, score, spare, row_entries, source)
  stat <- etwfe_nuisance_block_covariance(stat, rep.int(1L, length(keep)))
  variance <- stat$variance[match(coord, keep)]
  coef_role <- rep(c("level", rep("scaled_slope", length(nm))), length(all_slots))
  coef_role[rep(gauge, each = length(features))] <- "fixed_gauge"
  coefs <- chained_nuisance_coef_frame("calendar", term, coef_role, rep(seq_along(all_slots),
    each = length(features)), beta[coord], variance)
  transform <- data.frame(term = nm, center = as.double(center), amplitude = as.double(amplitude))
  eq <- switch(model,
    levels = "Y_it = alpha_i + W_i' b[slot(G_i,t)]",
    calendar_changes = "Y_it - Y_is = W_i' (b[t] - b[s]) + e_i(s,t)",
    corrected_levels = "Z_i(k;s,t) = W_i' (b[slot(k,t)] - b[slot(k,s)]) + e_i(k;s,t)",
    corrected_calendar = "Z_i(k;s,t) = W_i' (b[t] - b[s]) + e_i(k;s,t)")
  eq <- paste0(eq, if (length(nm)) " with W_i = (1, (Z_i - center) / amplitude)" else " with W_i = 1",
    ". The lowest slot of each fitted connected component is fixed at zero.")
  if (identical(model, "levels")) eq <- paste(eq, "Unit effects are profiled out.")
  if (is.null(comparisons)) comparisons <- data.frame(comparison = character(), fits = character())
  blocks <- data.frame(block = "calendar", label = "Anchored calendar coefficients",
    cohorts = chained_nuisance_gtxt(unique(cohorts)), n_parameters = length(keep),
    n_fitting_rows = fitting_rows, weight_sum = fitting_rows, n_comparisons = nrow(comparisons),
    kept_slopes = paste(nm, collapse = ","), aliased_slopes = "",
    n_contributing_clusters = unname(stat$clusters[1L]), equation = eq)
  fields <- data.frame(block = "calendar", cell = seq_along(all_slots), cohort = cohorts,
    date = dates, slot_type = ifelse(cohorts == 0, "common_date", "cohort_date"), fixed_gauge = gauge)
  res <- etwfe_nuisance_role("calendar outcome regression", eq,
    switch(model, levels = "Observed outcome levels Y_it", calendar_changes = "Observed untreated pair changes",
      "Observation-corrected pair response Z_i(k;s,t)"),
    sample, "Calendar levels enter fitted changes and cohort-standardized contrasts.",
    "one linked slot system", blocks, coefs, list(control_blocks = fields), comparisons,
    stat$covariance, features, etwfe_nuisance_formula(report$outcome_formula, features))
  res$basis_transform <- transform
  res
}

# Scores of within-unit least squares in the fitted feature coordinates. A unit's
# level residuals already sum to zero, so its raw slot design gives the
# same score as the within-transformed design.
etwfe_nuisance_within_score <- function(rows, y, slot, features, beta, n_obs, first,
                                       keep, solve_factor) {
  nb <- length(rows)
  n <- length(n_obs)
  widths <- n_obs[rows]
  origin <- rep.int(seq_len(nb), widths)
  ids <- rep.int(first[rows] - 1, widths) + sequence(widths)
  use <- slot[ids]
  W <- if (is.null(features)) matrix(1, length(ids), 1L) else features[rows[origin], , drop = FALSE]
  residual <- y[ids] - rowSums(W * beta[use, , drop = FALSE])
  centered <- residual
  if (length(ids)) centered <- residual - rep.int(drop(rowsum(residual, origin)) / widths[widths > 0L],
                                                   widths[widths > 0L])
  A <- nrow(beta)
  q <- ncol(beta)
  u <- matrix(0, nb * A, q)
  if (length(ids)) u[(use - 1L) * nb + origin, ] <- centered * W
  score <- matrix(u, nb, A * q)[, keep, drop = FALSE]
  n * t(solve_factor(t(score)))
}

etwfe_nuisance_observation <- function(st, scat, diagnostics, score, report) {
  blocks <- c(lapply(st$sel_blocks, `[[`, "index"), as.list(st$pi_index),
    lapply(st$imp_groups, `[[`, "index"))
  owner <- integer(st$P)
  for (m in seq_along(blocks)) owner[blocks[[m]]] <- m
  budget <- diagnostics$budget
  source <- "Uses the observation provider's fitted selection and actual-change imputation scores, summed by original unit."
  # The borrowed outcome, feature and roster arrays remain live during replay.
  spare <- budget$workspace_entries - budget$retained_entries - budget$borrowed_entries - budget$scratch_entries
  stat <- etwfe_nuisance_stats(st$P, st$N, report, score, spare, 2 * st$P, source, owner)
  part <- function(index) {
    ix <- as.integer(unlist(index, use.names = FALSE))
    local <- stat
    local$variance <- stat$variance[ix]
    etwfe_nuisance_block_covariance(local, owner[ix])
  }
  list(
    selection = if (st$need_q) etwfe_nuisance_observation_selection(st, scat, diagnostics,
      part(blocks[seq_along(st$sel_blocks)]), part(as.list(st$pi_index)), report),
    imputation = if (st$need_j) etwfe_nuisance_observation_imputation(st, diagnostics,
      part(lapply(st$imp_groups, `[[`, "index")), report))
}

etwfe_nuisance_observation_selection <- function(st, scat, diagnostics, stat, roster, report) {
  ds <- st$design
  N <- st$N
  is_date <- identical(scat$selection_model, "independent_period_logit")
  nc <- length(st$sel_cells)
  fits <- st$sel_blocks
  ids <- paste0("selection_", seq_along(fits))
  se_cells <- lapply(st$sel_cells, function(cell) {
    cell$cohort_value <- etwfe_nuisance_cohort(cell$cohort_value, report)
    if (is_date) {
      cell$t <- report$calendar[cell$calendar_index]
    } else {
      cell$pair$s <- report$calendar[cell$pair$s_calendar_index]
      cell$pair$t <- report$calendar[cell$pair$t_calendar_index]
    }
    cell
  })
  ci <- vapply(se_cells, `[[`, integer(1), "cohort")
  cblock <- vapply(se_cells, `[[`, integer(1), "block_id")
  keys <- vapply(se_cells, chained_nuisance_cell_key, "")
  var <- numeric(st$P)
  sel_pos <- as.integer(unlist(lapply(fits, `[[`, "index"), use.names = FALSE))
  var[sel_pos] <- stat$variance
  coefs <- do.call(rbind, lapply(seq_along(fits), function(m) {
    bf <- fits[[m]]
    intercepts <- if (bf$has_intercept) bf$cell_ids else integer()
    term <- c(paste0("a[", keys[intercepts], "]")[seq_along(intercepts)], ds$slope_names)
    role <- c(rep("local_intercept", length(intercepts)),
      rep(if (bf$B_m > 1L) "shared_slope" else "slope", bf$d_s))
    chained_nuisance_coef_frame(ids[m], term, role, c(intercepts, rep(NA_integer_, bf$d_s)),
      bf$coefficients, var[bf$index])
  }))
  pcs <- scat$pair_cohorts
  pc_cells <- lapply(pcs, function(pc) if (is_date) c(pc$s_cell_id, pc$t_cell_id) else pc$pair_cell_id)
  pc_blocks <- lapply(pc_cells, function(index) unique(cblock[index]))
  pair_cohorts <- data.frame(pair_cohort = seq_along(pcs),
    cohort = etwfe_nuisance_cohort(vapply(pcs, `[[`, 0, "cohort_value"), report),
    s = report$calendar[vapply(pcs, function(pc) pc$pair$s_calendar_index, integer(1))],
    t = report$calendar[vapply(pcs, function(pc) pc$pair$t_calendar_index, integer(1))],
    probability = if (is_date) "Product of two date logits under conditional date independence" else "Direct pair logit",
    cells = vapply(pc_cells, paste, "", collapse = ","),
    blocks = vapply(pc_blocks, function(b) paste(ids[b], collapse = ","), ""))
  diag <- lapply(fits, `[[`, "diagnostics")
  blocks <- data.frame(block = ids,
    label = vapply(fits, function(x) chained_nuisance_block_label(scat$selection_pooling, se_cells, x$cell_ids), ""),
    cohorts = vapply(fits, function(x)
      chained_nuisance_gtxt(unique(etwfe_nuisance_cohort(st$coh_vals[ci[x$cell_ids]], report))), ""),
    n_cells = lengths(lapply(fits, `[[`, "cell_ids")),
    n_parameters = lengths(lapply(fits, `[[`, "index")),
    n_fitting_rows = vapply(fits, function(x) sum(st$nU[ci[x$cell_ids]]), 0),
    n_units = vapply(fits, function(x) sum(st$nU[unique(ci[x$cell_ids])]), 0),
    n_fitting_clusters = unname(stat$clusters[as.character(seq_along(fits))]),
    branch = vapply(fits, `[[`, "", "branch"),
    converged = vapply(diag, function(x) isTRUE(x$converged), FALSE),
    iterations = vapply(diag, function(x) etwfe_tm_scalar(x$iter), 0),
    rank = vapply(diag, function(x) etwfe_tm_scalar(x$rank), 0),
    n_comparisons = tabulate(unlist(pc_blocks, use.names = FALSE), length(fits)),
    equation = vapply(fits, function(x) chained_nuisance_block_equation(scat$selection_model,
      se_cells, x$cell_ids, x$has_intercept, ds$slope_names), ""))
  rate <- vapply(seq_len(nc), function(i) fits[[cblock[i]]]$pi[se_cells[[i]]$local_b], 0)
  cells <- data.frame(cell = seq_len(nc), block = ids[cblock], type = if (is_date) "date" else "pair",
    cohort = etwfe_nuisance_cohort(st$coh_vals[ci], report), key = keys,
    s = vapply(se_cells, function(x) if (is_date) NA_real_ else x$pair$s, 0),
    t = vapply(se_cells, function(x) if (is_date) x$t else x$pair$t, 0),
    n_rows = st$nU[ci], n_observed = round(st$nU[ci] * rate), observed_share_weighted = rate,
    p_min = NA_real_, p_max = NA_real_, n_floor = NA_integer_)
  roster_table <- if (length(st$pi_index)) {
    groups <- etwfe_nuisance_cohort(st$coh_vals, report)
    chained_nuisance_coef_frame("roster", paste0("Pr(U | G=", groups, ")"),
      rep("roster_probability", length(groups)), seq_along(groups), st$pi_val, roster$variance)
  }
  list(status = "available", route = "cohort_time_fe", role = "observation propensity",
    selection_model = scat$selection_model, selection_pooling = scat$selection_pooling, link = "logit",
    response = if (is_date) "R_it = 1{Y_it observed}" else "S_i(s,t) = 1{Y_is and Y_it both observed}",
    template = chained_nuisance_selection_template(scat$selection_model, scat$selection_pooling, ds$branch),
    derived = if (is_date) "Each pair probability is the product of its two fitted date probabilities.",
    formula = etwfe_nuisance_formula(report$selection_formula, colnames(st$Z)),
    design_columns = colnames(st$Z), slope_terms = ds$slope_names, design_branch = ds$branch,
    universe = st$universe, universe_text = chained_nuisance_universe_text(st$universe),
    weights = list(text = "unit weights on the fitting roster"), n_frame = N,
    use = diagnostics$selection$law, probability_floor = ipw_clip_eps,
    blocks = blocks, coefficients = coefs, cells = cells, pair_cohorts = pair_cohorts,
    covariance = stat$covariance, roster = roster_table,
    roster_covariance = if (!is.null(roster_table)) roster$covariance)
}

etwfe_nuisance_observation_imputation <- function(st, diagnostics, stat, report) {
  fits <- st$imp_groups
  groups <- etwfe_nuisance_cohort(st$cell_k_val, report)
  s <- report$calendar[st$cell_s]
  t <- report$calendar[st$cell_t]
  keys <- paste0("g=", groups, ",", s, ":", t)
  ids <- paste0("imputation_", seq_along(fits))
  var <- numeric(st$P)
  pos <- as.integer(unlist(lapply(fits, `[[`, "index"), use.names = FALSE))
  var[pos] <- stat$variance
  coefs <- do.call(rbind, lapply(seq_along(fits), function(m) {
    f <- fits[[m]]
    intercepts <- if (st$icpt) f$cells else integer()
    term <- c(paste0("a[", keys[intercepts], "]")[seq_along(intercepts)], f$kept_columns)
    role <- c(rep("local_intercept", length(intercepts)), rep("slope", f$p))
    chained_nuisance_coef_frame(ids[m], term, role, c(intercepts, rep(NA_integer_, f$p)),
      c(f$alpha, f$beta), var[f$index])
  }))
  fit_count <- vapply(diagnostics$cells, `[[`, 0, "n_pair")
  rows <- vapply(fits, function(f) sum(fit_count[f$cells]), 0)
  formula <- function(cols) {
    rhs <- c(if (st$icpt) "a[k,s,t]", if (length(cols)) paste0("(t-s) * (",
      paste(paste0("b[", cols, "] * ", cols), collapse = " + "), ")"))
    paste("E[Y_it - Y_is | G_i = k, X_i, S_i(s,t) = 1] =", if (length(rhs)) paste(rhs, collapse = " + ") else "0")
  }
  offset <- length(st$sel_blocks) + length(st$pi_index)
  blocks <- data.frame(block = ids,
    label = vapply(fits, function(f) paste(keys[f$cells], collapse = ", "), ""),
    cohorts = vapply(fits, function(f) chained_nuisance_gtxt(unique(groups[f$cells])), ""),
    n_parameters = lengths(lapply(fits, `[[`, "index")), n_fitting_rows = rows, weight_sum = rows,
    n_comparisons = lengths(lapply(fits, `[[`, "cells")),
    n_contributing_clusters = unname(stat$clusters[as.character(offset + seq_along(fits))]),
    kept_slopes = vapply(fits, function(f) paste(f$kept_columns, collapse = ","), ""),
    aliased_slopes = vapply(fits, function(f) paste(f$aliased_columns, collapse = ","), ""),
    equation = vapply(fits, function(f) formula(f$kept_columns), ""))
  cells <- data.frame(block = ids[st$cell_group], cell = seq_along(keys), cohort = groups,
    s = s, t = t, L = st$cell_L, n_observed_pairs = fit_count)
  comparisons <- etwfe_nuisance_comparisons(st$cell_k_val, st$cell_s, st$cell_t,
    rep(NA_integer_, length(groups)), cells$block, report)
  etwfe_nuisance_role("actual-change imputation", formula(diagnostics$imputation$slope_columns),
    "Observed outcome change on a supported pair", "Observed pairs within each cohort",
    "Predicts the actual change for each frame unit of the cohort in the reconstructed pair response.",
    diagnostics$imputation$pooling, blocks, coefs, list(cells = cells), comparisons,
    stat$covariance, colnames(st$X), etwfe_nuisance_formula(report$outcome_formula, colnames(st$X)))
}

etwfe_nuisance_projection <- function(catalog, geometry, fit, obs, Z, report, entries, model, comparisons) {
  keep_slots <- geometry$keep_slots
  if (!length(keep_slots)) return(NULL)
  n <- nrow(catalog$Y)
  K <- ncol(catalog$Y)
  p <- length(fit$beta)
  coord <- matrix(seq_len(p), nrow(fit$beta))
  keep <- as.vector(coord[keep_slots, , drop = FALSE])
  groups <- if (identical(model, "corrected_calendar")) 0L else catalog$groups
  slots <- etwfe_nuisance_slot_codes(seq_len(nrow(fit$beta)), K, groups)
  slots$gauge <- slots$slot %in% setdiff(unique(c(geometry$left, geometry$right)), keep_slots)
  center <- if (is.null(Z)) numeric() else colMeans(Z)
  amplitude <- vapply(seq_along(center), function(j) max(abs(Z[, j] - center[j])), 0)
  amplitude[amplitude == 0] <- 1
  cells <- catalog$cells[geometry$use, , drop = FALSE]
  fitting_rows <- sum(lengths(catalog$units)[cells$group])
  obs_budget <- obs$diagnostics$budget
  spare <- entries - catalog$entries - obs_budget$retained_entries -
    obs_budget$scratch_entries - fit$live_entries
  if (!is.null(fit$block_rows)) report$chunk_rows <- min(report$chunk_rows, fit$block_rows)
  score <- function(rows) fit$score(rows, obs$score_rows(rows))[, keep, drop = FALSE]
  source <- "Includes the generated pair response contributions from the fitted observation and imputation scores, alongside the projection's direct original-unit score."
  etwfe_nuisance_outcome(fit$beta, keep, coord, slots, Z, center, amplitude, n, report, score, spare,
    6 * p + 2 * obs$n_parameters + nrow(catalog$cells) + 3 * ncol(fit$beta) + 8,
    fitting_rows, model, "Full cohort-frame units on supported pair equations", source, comparisons)
}
