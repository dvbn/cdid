# The treatment moments have the same incidence map as the existing cross-pair
# solver. Each "calendar" vertex is standardized to the declared target cohort.
# At its reference date the target and calendar paths have the same origin.
etwfe_treatment_core <- function(index, cohort, n, calendar, covariates, gate,
                                 est_method, resolved, entries, workers, biters,
                                 cluster_full, progress, moments = NULL, report = NULL) {
  K <- length(calendar)
  treated <- sort(unique(cohort[cohort > 0L]))
  nyt <- identical(gate$control_group, "notyettreated")
  post <- expand.grid(period = seq_len(K), cohort = treated)
  pre <- post$period < post$cohort
  targets <- post[if (nyt) !pre else post$period != post$cohort - 1L, ]
  targets$is_pre <- targets$period < targets$cohort
  targets$att <- NA_real_
  observed <- tabulate((match(cohort[index$unit], treated, nomatch = 0L) - 1L) * K + index$period,
                       nbins = K * length(treated))
  targets$status <- ifelse(observed[(match(targets$cohort, treated) - 1L) * K + targets$period] > 0L,
                           "unidentified", "no_support")

  weighting <- gate$gmm_weighting_requested
  optimal <- identical(gate$gmm_weighting, "optimal")
  scored <- gate$compute_inffunc || optimal
  result <- list(targets = targets, inffunc = NULL)
  if (gate$compute_inffunc) result$inffunc <- matrix(NA_real_, n, nrow(targets))

  if (!is.null(progress)) progress("cohort_time_fe: fitting treatment moments on cohort-specific averages")
  if (is.null(moments)) {
    moments <- etwfe_treatment_moments(
      index$y, index$unit, index$period, cohort, n, K, covariates,
      gate$control_group, est_method, scored, entries,
      list(epsilon = resolved$glm$epsilon, maxit = as.integer(resolved$glm$maxit), trace = FALSE),
      report = report
    )
  }
  result$diagnostics <- moments$diagnostics
  if (!is.null(report)) result$nuisance <- moments$nuisance
  links <- moments$links
  L <- nrow(links)
  result$diagnostics$gmm_weighting <- weighting
  result$diagnostics$n_moments <- L
  if (!L) return(result)
  ledger <- data.frame(g = targets$cohort[!targets$is_pre], time = targets$period[!targets$is_pre])
  design <- moments$design
  if (is.null(design)) {
    design <- build_cross_pair_chained_design(
      links, seq_len(K), c(treated, Inf), anticipation = 0L, pt_assumption = gate$pt_assumption,
      pt_start = if (nyt) 1L else NULL, ledger = ledger
    )
  }

  ids <- paste(targets$cohort, targets$period, sep = ":")
  specs <- lapply(seq_len(nrow(targets)), function(j) {
    g <- targets$cohort[j]
    t <- targets$period[j]
    if (targets$status[j] == "no_support") return(NULL)
    roles <- if (targets$is_pre[j]) c("pre_target", "pre") else c("post_target", "trend_post")
    a <- cross_pair_find_param(design$parameters, g, t, roles[1L])
    b <- cross_pair_find_param(design$parameters, g, t, roles[2L])
    if (is.na(a) || is.na(b)) return(NULL)
    gmm_target_spec(ids[j], c(a, b), c(1, -1))
  })
  keep <- which(!vapply(specs, is.null, logical(1)))
  if (!length(keep)) return(result)

  dims <- chained_public_dimensions_from_counts(
    n_frame = n, n_panel_units = length(unique(index$unit)), n_rows = length(index$y),
    n_dates = K, n_cohorts = length(treated) + 1L, n_targets = length(treated),
    n_covariates = ncol(covariates) + 1L, biters = biters,
    n_clusters = if (is.null(cluster_full)) n else length(unique(cluster_full)),
    calendar = seq_len(K), target_cohorts = treated
  )
  dims$links_bound <- L
  ctrl <- chained_public_pipeline_control(dims, resolved, workers = workers)$control
  chunk <- min(n, ctrl$gmm_bounds$max_batch_rows)
  if (!is.null(moments$solver_reserve)) {
    reserve <- moments$solver_reserve
    if (scored) chunk <- min(chunk, reserve$max_batch_rows)
    solve_entries <- entries - reserve$fixed - scored * chunk * reserve$row_scratch -
      gate$compute_inffunc * as.double(n) * nrow(targets)
    if (solve_entries < 1) etwfe_tm_workspace_error()
    ctrl$gmm_bounds$max_matrix_entries <- solve_entries
    ctrl$gmm_bounds$max_materialized_if_entries <- min(
      ctrl$gmm_bounds$max_materialized_if_entries, solve_entries)
    ctrl$gmm_bounds$max_batch_rows <- chunk
  }
  if (6 * as.double(L)^2 + 2 * as.double(chunk) * L * scored +
      3 * as.double(n) * nrow(targets) * gate$compute_inffunc > entries) {
    cdid_abort("The cohort-by-date treatment moment solve exceeds numeric_workspace_bytes.",
               "cdid_etwfe_workspace_error")
  }
  if (!is.null(progress)) progress("cohort_time_fe: combining treatment moments and original-unit scores")
  unit_ids <- as.character(seq_len(n))
  row_first <- seq.int(1L, n, by = chunk)
  provenance <- list(
    n_positive_contributing_rows = n, centered = TRUE, weighted = FALSE,
    weight_type = "unweighted", normalization_factor = 1 / n,
    scale_id = "etwfe_target_standardized_full_units",
    nuisance_fit_id = paste0("cohort_time_fe:", est_method)
  )
  prov <- do.call(gmm_covariance_provenance,
                  c(list(row_level = "unit", n_units = n, n_clusters = n), provenance))
  cov <- matrix(0, L, L)
  replay <- function(consume) {
    for (s in row_first) {
      ids <- seq.int(s, min(n, s + chunk - 1L))
      values <- moments$score_rows(ids)
      colnames(values) <- design$incidence$link_ids
      consume(list(unit_id = unit_ids[ids], phi = values,
                   weight = rep.int(1, length(ids)), nuisance_fit_id = provenance$nuisance_fit_id))
    }
  }
  holder <- NULL
  if (gate$compute_inffunc) {
    fit_entries <- preflight_gmm_allocations(
      L, design$rank_diagnostics$rank, length(keep), gate$gmm_weighting, ctrl$gmm_bounds,
      n_out_rows = n, n_parameters = design$incidence$n_parameters, materialize_target_if = TRUE,
      certify_error = FALSE, covariance_strategy = ctrl$gmm_accuracy$covariance_strategy
    )$total_payload_entries
    score_budget <- moments$diagnostics$workspace
    outstanding <- as.double(n) * nrow(targets) + as.double(L)^2 + fit_entries +
      score_budget$fixed + chunk * (score_budget$score_row + L)
    if (outstanding + static_missing_link_if_replay_entries(n, chunk, L) <= entries) {
      holder <- static_missing_link_if_replay_holder(replay, unit_ids, chunk, L)
      replay <- holder$callback
    }
  }
  result$diagnostics$retained_score_batches <- !is.null(holder)
  if (optimal) {
    replay(function(batch) cov <<- cov + crossprod(batch$phi))
    cov <- cov / n
    prov$complete_unit_if_gram <- ipw_seal(
      list(covariance = cov, n_units = n, nuisance_fit_id = provenance$nuisance_fit_id),
      "complete_unit_if_gram"
    )
  }
  replay_spec <- do.call(gmm_if_replay_spec, c(list(
    callback = replay, input_row_level = "unit", aggregate_to_clusters = FALSE,
    n_input_rows = n, n_output_rows = n, expected_unit_ids = unit_ids), provenance))
  fit <- solve_chained_gmm_core(
    design, moments$mhat, cov, prov, specs[keep], weighting, replay_spec,
    ctrl$gmm_accuracy, ctrl$gmm_bounds, gauge_rule = "lowest",
    materialize_target_if = gate$compute_inffunc,
    project_influence = gate$compute_inffunc
  )
  if (!is.null(holder)) holder$release()
  present <- which(fit$targets$estimable)
  mapped <- match(fit$targets$target_id[present], ids)
  result$targets$att[mapped] <- fit$targets$estimate[present]
  result$targets$status[mapped] <- "identified"
  if (gate$compute_inffunc && length(mapped)) {
    result$inffunc[, mapped] <- fit$if_projection$materialized_target_if
  }
  result$diagnostics$covariance_handling <- fit$covariance_handling[
    c("handler_used", "trigger_reason", "lambda", "covariance_rank")
  ]
  result
}
