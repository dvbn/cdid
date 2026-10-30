# CP08 O14 (D63/D66): chained conditional pretrend test executor and result.
#
# H0: delta_c = 0 target-almost-everywhere for every declared cell c of the
# selected base-period family (D21 cells, reference zeros excluded), i.e.
# E_g[delta_c | z_c] = 0 with z_c the model.matrix(xformla) row on the cell's
# own recipe. The process J_c(u) = sigma_c E_g[1{z_c <= u} delta_c] is the
# accepted auxiliary link with h_u inserted into its numerator integrands; the
# union-catalog nuisance fits are reused unchanged (no refit by grid point, no
# masking). CvM = (N/M) sum_c sum_m J_c(u_m)^2 with one grid row per original
# unit (M = N); multiplier draws share one sign stream with the existing
# pipeline conventions. This is not the D21 Wald test (att_gt W/Wpval).

cpt_meta_base <- function(status = "unavailable", reason = NA_character_) {
  list(test = "chained conditional pretrend test (Cramer-von Mises)",
    status = status, reason = reason,
    null = paste("For every declared pre-treatment cell c of the selected base-period family (reference zeros excluded):",
      "E_g[delta_c | z_c] = 0, with z_c the model.matrix(xformla) row of the unit at the cell's own reference date;",
      "it can be weaker than conditioning on all scientific variables."),
    not_wald = "not the unconditional Wald test: att_gt()'s unconditional W/Wpval is separate and unchanged",
    ks = "KS fields are NULL (not computed), as in did")
}

cpt_result <- function(meta, clustervars, xformla, CvM = NA_real_, CvMb = NULL, CvMcval = NA_real_, CvMpval = NA_real_) {
  structure(list(CvM = CvM, CvMb = CvMb, CvMcval = CvMcval, CvMpval = CvMpval, KS = NULL, KSb = NULL, KScval = NULL,
    KSpval = NULL, clustervars = clustervars, xformla = xformla, chained = meta), class = c("cdid_MP.TEST", "MP.TEST"))
}

cpt_cell_table <- function(cells) {
  if (!length(cells)) return(data.frame(cell_id = character(), g = numeric(), t = numeric(), stringsAsFactors = FALSE))
  data.frame(cell_id = vapply(cells, `[[`, character(1L), "id"),
    g = vapply(cells, `[[`, numeric(1L), "g"), t = vapply(cells, `[[`, numeric(1L), "t"),
    from = vapply(cells, `[[`, numeric(1L), "from"), sign = vapply(cells, `[[`, numeric(1L), "sign"),
    pair = vapply(cells, function(x) if (is.na(x$s)) NA_character_ else
      paste0(fmt_design_num(x$s), "->", fmt_design_num(x$t2)), character(1L)),
    controls = vapply(cells, function(x) paste0("{", paste(fmt_design_num(x$K), collapse = ","), "}"), character(1L)),
    status = vapply(cells, `[[`, character(1L), "status"),
    reason = vapply(cells, function(x) if (is.null(x$reason)) NA_character_ else x$reason, character(1L)),
    stringsAsFactors = FALSE)
}

# Declared family of the selected base period, with the auxiliary descriptor
# status of the fitted registry.
cpt_declared_cells <- function(state) {
  sels <- state$aux_plan$selections[[state$base_period]]
  descs <- state$aux_plan$descriptors
  ds <- state$auxiliary$descriptor_status
  out <- list()
  for (d in sels) {
    if (identical(d$status, "reference_zero")) next
    id <- paste0("pre(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")")
    cell <- list(id = id, g = as.double(d$g), t = as.double(d$t), from = as.double(d$from), sign = as.double(d$sign),
      K = as.double(d$K), q = as.integer(d$index), s = NA_real_, t2 = NA_real_, status = "pending", reason = NULL)
    if (identical(d$status, "unavailable")) {
      cell$status <- "unavailable"; cell$reason <- d$reason
    } else {
      cell$s <- descs[[d$index]]$s; cell$t2 <- descs[[d$index]]$t
      st <- if (is.null(ds)) NA_character_ else ds$status[d$index]
      if (!identical(st, "supported")) {
        cell$status <- "unavailable"
        cell$reason <- paste0("unsupported auxiliary comparison: ", if (is.null(ds)) "no fitted registry" else ds$reason[d$index])
      }
    }
    out[[length(out) + 1L]] <- cell
  }
  out
}

# Prospective numeric payload and work of the test's own buffers, charged on
# top of the fitted registry/owners with the existing caps and multiplier budget.
cpt_preflight <- function(dims, N, C, B, br, gchunk, echunk, ctrl, mb, retained, work, workers) {
  oc <- ctrl$owner_control
  cap <- min(oc$max_retained_entries, oc$max_matrix_entries)
  S <- function(x, ...) ipw_owner_safe_sum(as.double(x), ...)
  # Cache retains the full N x p grid as well as its distinct design/map.
  grids <- S(vapply(dims, function(d) N * (d$p + 1) + d$Md * (d$p + 3), numeric(1L)),
    label = "conditional test grids")
  cells <- S(vapply(dims, function(d) d$Md * (d$R + d$Dl + d$p + 3) + C * (d$R + d$Dl) + 3 * N + d$Dl * d$R + 64 +
    if (is.null(d$rows_retained)) 0 else d$rows_retained, numeric(1L)), label = "conditional test cell states")
  # Extraction holds nu/rho, their products, all term loadings, nuisance
  # details, accumulator copies and the rowsum/matrix-product temporaries.
  build <- max(vapply(dims, function(d) 4 * d$Md * (d$R + d$Dl) +
    64 * d$chunk * (2 * d$R + 3 * d$Dl + 8) + 4 * d$Md * gchunk +
    4 * gchunk * (d$R + d$Dl) + 8 * N, numeric(1L)))
  batch <- max(vapply(dims, function(d) br * (2 * C + 2 * d$Dl + 2 * d$R + 3 * d$Md +
    6 * gchunk + 6 * echunk) + 4 * d$Md * gchunk, numeric(1L)))
  # Parent signs/results, complete worker cell copies and their serialization
  # buffers coexist on PSOCK. This also conservatively covers the fork route.
  parallel <- if (workers > 1) workers * (2 * cells + 2 * br * C + 2 * br + 256) else 0
  draws <- batch * max(1, workers) + parallel + 2 * B
  degeneracy <- 6 * C * gchunk + 4 * echunk * gchunk + 4 * gchunk *
    max(vapply(dims, function(d) d$p + d$R + d$Dl, numeric(1L)))
  total <- ipw_owner_safe_sum(c(retained, grids, cells, max(build, draws, degeneracy)), cap,
    "conditional pretest grids, cell states and simultaneous build, multiplier or degeneracy buffers")
  wk <- ipw_owner_safe_sum(c(work,
    S(vapply(dims, function(d) 64 * N * (d$R + d$Dl) + 16 * d$Md * d$Md * (d$p + d$R + d$Dl), numeric(1L)), label = "build"),
    S(vapply(dims, function(d) 8 * B * (C * (d$R + d$Dl) + d$Md * d$Md * (d$p + 1) + d$Md * (d$R + d$Dl) + N), numeric(1L)),
      label = "draws"),
    S(vapply(dims, function(d) 8 * (C * d$Md * (d$R + d$Dl) + N * d$Md * (d$p + 1)), numeric(1L)),
      label = "cluster-process degeneracy check")), oc$max_work_ops, "conditional pretest build and multiplier work")
  mult <- batch * max(1, workers) + parallel + 2 * B + C
  if (B > mb$max_draws || br < 1 || mult > mb$max_workspace_cells || 8 * mult > mb$max_workspace_bytes)
    cpt_abort("Conditional pretest multiplier draws exceed the caller multiplier budget.",
      c("cdid_conditional_pretest_bound_error", "cdid_ipw_owner_bound_error"))
  list(retained_entries = total, work_ops = wk, multiplier_entries = mult)
}

# The test on a conditional pretest state of the accepted pipeline.
chained_conditional_test <- function(state, prep, src, xformla, alp, biters, mb, workers, route_text,
                                     fixed_signs = NULL, keep_rows = FALSE, progress = NULL) {
  meta <- cpt_meta_base()
  meta$base_period <- state$base_period
  meta$route <- route_text
  N <- state$n_total
  if (!identical(state$status, "available")) {
    meta$reason <- state$reason
    return(list(meta = meta))
  }
  cells <- cpt_declared_cells(state)
  meta$cells <- cpt_cell_table(cells)
  if (!length(cells)) {
    meta$reason <- paste0("No declared pre-treatment cell exists in the base_period = \"", state$base_period,
      "\" family (reference zeros excluded): there is no conditional pretrend restriction to test.")
    return(list(meta = meta))
  }
  bad <- which(vapply(cells, function(x) identical(x$status, "unavailable"), logical(1L)))
  if (length(bad)) {
    meta$reason <- paste0("The complete requested test is unavailable: ", length(bad), " declared cell(s) are unavailable (",
      paste(vapply(cells[bad], function(x) paste0(x$id, ": ", x$reason), character(1L)), collapse = "; "),
      "). No available-only subset is tested.")
    return(list(meta = meta))
  }
  reg <- state$registry
  backend <- reg$cp04
  ref <- if (is.null(backend)) NULL else chained_backend_reference(backend)
  # Grids on each cell's own recipe.
  src$grid_budget <- list(retained = state$resource$retained_entries, cached = 0,
    cap = min(state$control$owner_control$max_retained_entries, state$control$owner_control$max_matrix_entries))
  cache <- list(); keys <- character(length(cells)); gbad <- list()
  for (i in seq_along(cells)) {
    recipe <- if (is.null(ref)) NA_integer_ else ref$view_recipe[reg$cp04_aux_index[cells[[i]]$q]]
    g <- cpt_grid_for_cell(cache, src, prep, xformla, ref, recipe, cells[[i]]$s)
    cache <- g$cache; keys[i] <- g$key
    if (!identical(g$grid$status, "available")) gbad[[cells[[i]]$id]] <- g$grid
  }
  if (length(gbad)) {
    g1 <- gbad[[1L]]
    meta$reason <- paste0("The complete requested test is unavailable: required grid inputs are missing or unsupported (",
      paste(paste0(names(gbad), ": ", vapply(gbad, `[[`, character(1L), "reason")), collapse = "; "),
      "). Estimated ATTs are unaffected; no unit is deleted.")
    meta$grid <- list(missing = g1$missing, unsupported = g1$unsupported)
    return(list(meta = meta))
  }
  # Clusters: full original-unit alignment, zero-IF clusters included.
  cl <- state$cluster_full
  u_cl <- if (is.null(cl)) seq_len(N) else unique(cl)
  cmap <- if (is.null(cl)) seq_len(N) else match(cl, u_cl)
  C <- length(u_cl)
  if (!is.null(cl) && C < 2L) {
    meta$reason <- "Clustered multiplier inference needs at least two clusters; the conditional test is unavailable."
    return(list(meta = meta))
  }
  oc <- state$control$owner_control
  chunk_rows <- oc$max_chunk_rows
  state$retained_now <- state$resource$retained_entries
  if (!is.null(progress)) progress("conditional pretest: building ", length(cells), " cell process(es) on ", N, " grid rows")
  ex <- lapply(cells, function(x) cpt_cell_extractor(state, x$q, x, chunk_rows))
  B <- as.integer(biters)
  br <- as.integer(min(mb$max_sign_batch_rows, floor(mb$max_sign_batch_cells / C), B))
  gchunk <- as.integer(chunk_rows); echunk <- as.integer(chunk_rows)
  dims <- lapply(seq_along(cells), function(i) {
    gd <- cache[[keys[i]]]
    list(Md = gd$distinct$n_distinct, p = ncol(gd$Z), R = length(ex[[i]]$spec$s), Dl = length(ex[[i]]$spec$coords),
      chunk = chunk_rows, rows_retained = if (keep_rows) N * (1 + length(ex[[i]]$spec$s) + length(ex[[i]]$spec$coords)) else 0)
  })
  pf <- cpt_preflight(dims, N, C, B, br, gchunk, echunk, state$control, mb, state$resource$retained_entries,
    state$resource$work_ops, workers)
  built <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    grid <- cache[[keys[i]]]$distinct
    e <- ex[[i]]
    acc <- cpt_cell_acc(e$spec, grid, C, N, keep_rows)
    for (k in seq_len(e$n_chunks)) for (pc in e$pieces(k)) acc <- cpt_cell_add(acc, e$spec, pc, grid, cmap)
    built[[i]] <- cpt_cell_finalize(acc, e$spec, grid, N, gchunk)
    built[[i]]$estimate_scalar <- cells[[i]]$sign * e$estimate
    built[[i]]$route <- e$route
    built[[i]]$recipe_key <- keys[i]
  }
  CvM <- cpt_statistic(built, N)
  if (!is.null(progress)) progress("conditional pretest: ", B, " multiplier draws on ", C, if (is.null(cl)) " units" else " clusters")
  mult <- cpt_multiplier(built, N, C, B, br, gchunk, echunk, fixed_signs = fixed_signs, workers = workers)
  CvMb <- mult$CvMb
  scal <- data.frame(cell_id = vapply(built, `[[`, character(1L), "id"),
    estimate = vapply(built, `[[`, numeric(1L), "estimate_scalar"),
    scalar_from_process = vapply(built, function(x) x$sigma * x$mu1, numeric(1L)),
    n_grid_rows = N, n_distinct_grid_points = vapply(built, function(x) nrow(x$zd), integer(1L)),
    recipe = keys, n_coordinates = vapply(built, `[[`, integer(1L), "n_coordinates"),
    stringsAsFactors = FALSE)
  scal$cell_id <- NULL
  meta$cells <- cbind(meta$cells, scal[match(meta$cells$cell_id, vapply(built, `[[`, character(1L), "id")), , drop = FALSE])
  meta$cells$status <- "tested"
  rownames(meta$cells) <- NULL
  meta$grid <- list(rule = paste("one row per supplied original unit (frame order; singletons, zero-outcome and frame-only",
      "units included; duplicate rows keep their mass), model.matrix(xformla) columns at each cell's reference date",
      "through the route's frozen terms, levels, contrasts and transform state"),
    columns = lapply(cache, `[[`, "columns"), n_units = N, M = N,
    indicator = "1{z_i <= u} componentwise (did's indicator weight function)")
  meta$statistic <- list(CvM = "(N/M) sum_c sum_m J_c(u_m)^2 with M = N and unstandardized J",
    multiplier = "J*_b(u) = N^-1 sum_k eta_bk S_k(u); CvM_b analogous; p = mean(CvM_b > CvM); critical value type-1 quantile",
    standardized = FALSE)
  meta$inference <- list(B = B, n_clusters = C, clustered = !is.null(cl), cluster_variable = prep$clustervar,
    rng = if (is.null(fixed_signs)) list(mode = "caller_rng_stream_Rademacher", kind = RNGkind(), batch_rows = br,
      ordering = "batches in increasing draw order; column-major within each batch; one sign matrix shared by every cell and grid point",
      caller_rng_restored = FALSE) else list(mode = "caller_fixed_signs", batch_rows = br),
    workers = workers, route = mult$route,
    qualifications = paste("Calibration requires the uniform process and multiplier-bootstrap regularity conditions of the",
      "technical note, including a growing number of independent clusters for clustered inference; scalar ATT",
      "asymptotic normality alone does not establish them, and finite numerical checks establish algebra only,",
      "not nominal size or power."))
  meta$resource <- pf
  meta$denominators <- lapply(built, `[[`, "denominator_check")
  if (!cpt_process_nonzero(built, C, gchunk, echunk)) {
    meta$reason <- paste("Globally degenerate process: every cluster process row is exactly zero (no estimated process",
      "varies across units), so the conditional test is unavailable.")
    meta$CvM <- CvM
    return(list(meta = meta, cells = if (keep_rows) built else NULL))
  }
  meta$status <- "available"; meta$reason <- NA_character_
  list(meta = meta, CvM = CvM, CvMb = CvMb, CvMpval = mean(CvMb > CvM),
    CvMcval = unname(stats::quantile(CvMb, probs = 1 - alp, type = 1, names = FALSE)),
    cells = if (keep_rows) built else NULL)
}

# Chained-mode executor of the public conditional_did_pretest().
chained_conditional_did_pretest <- function(yname, tname, idname, gname, xformla, data, panel, allow_unbalanced_panel,
                                            control_group, weightsname, alp, bstrap, cband, biters, clustervars, est_method,
                                            print_details, pl, cores, anticipation, base_period, fix_weights, mode, frame,
                                            selection, selection_model, selection_pooling, selection_xformla,
                                            selection_universe, selection_rho, selection_justification,
                                            selection_est_method, imputation_pooling, control, call = NULL,
                                            fixed_signs = NULL, keep_rows = FALSE) {
  gate <- chained_request_scope_gate(mode = mode, panel = panel, allow_unbalanced_panel = allow_unbalanced_panel,
    est_method = est_method, selection = selection, selection_model = selection_model, selection_pooling = selection_pooling,
    selection_xformla = selection_xformla, selection_universe = selection_universe, selection_rho = selection_rho,
    selection_justification = selection_justification, selection_est_method = selection_est_method,
    imputation_pooling = imputation_pooling, pt_assumption = "post", pt_start = NULL,
    gmm_weighting = c("optimal", "identity"), gmm_covariance = "tikhonov", fix_weights = fix_weights,
    compute_inffunc = TRUE, pl = pl, cores = cores, faster_mode = TRUE, print_details = print_details,
    control_group = control_group, base_period = base_period, anticipation = anticipation, alp = alp,
    bstrap = bstrap, cband = cband, biters = biters, dots = list(), entry = "public")
  meta <- cpt_meta_base()
  meta$base_period <- base_period
  meta$request <- list(est_method = if (isTRUE(gate$is_custom_est)) "custom" else est_method,
    control_group = gate$control_group, anticipation = as.double(anticipation), selection = gate$selection,
    selection_model = gate$selection_model, selection_pooling = gate$selection_pooling,
    selection_est_method = gate$selection_est_method, imputation_pooling = gate$imputation_pooling,
    fix_weights = gate$fix_weights, biters = as.double(biters), alp = alp,
    bstrap_cband_note = "as in did, the conditional test always uses multiplier draws; bstrap and cband are validated only",
    call = call)
  if (isTRUE(gate$is_custom_est)) {
    meta$status <- "not_applicable"
    meta$reason <- paste("A custom function-valued est_method has no claimed conditional process: the conditional",
      "test is defined for the built-in IPW, OR and DR routes only. Nothing was fitted.")
    return(cpt_result(meta, clustervars, xformla))
  }
  resolved <- chained_public_resolve_control(control, gmm_weighting = gate$gmm_weighting_requested,
    gmm_covariance = gate$gmm_covariance)
  plan <- chained_request_parallel_plan(gate$pl, gate$cores, resolved$profile$workers)
  progress <- chained_progress_emitter(gate$print_details)
  if (!is.null(progress)) progress("chained conditional_did_pretest(): validating and aligning the input by unit identity")
  est_label <- est_method
  imputation <- identical(gate$selection_est_method, "imputation")
  augmented <- identical(gate$selection_est_method, "dr")
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
  if ((imputation || augmented || est_label %in% c("reg", "dr")) && !is.null(prep$cp04))
    prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- chained_public_dimensions(prep, biters = biters, clustered = !is.null(prep$cluster_full))
  produced <- chained_public_pipeline_control(dims, resolved, gmm_weighting = gate$gmm_weighting_requested,
    gmm_covariance = gate$gmm_covariance, workers = plan$workers)
  context <- function(e) chained_request_refusal_context("public", prep,
    settings = list(control_group = gate$control_group, base_period = base_period, anticipation = as.double(anticipation),
      alp = alp, biters = as.double(biters), xformla = deparse(prep$xformla), weightsname = weightsname,
      clustervars = prep$clustervar, est_method = est_label, selection_est_method = gate$selection_est_method,
      imputation_pooling = gate$imputation_pooling, entry = "conditional_did_pretest", multiplier_workers = plan$workers),
    policy = list(profile = resolved$profile, profile_status = resolved$profile_status, source = resolved$source,
      overrides = resolved$overrides, chunk_rows = produced$record$chunk_rows,
      numeric_workspace_entries = produced$record$numeric_workspace_entries,
      counted_work_limit = produced$record$counted_work_limit))
  state <- ipw_pipeline_resignal_with_state(run_chained_ipw_pipeline(
    data = prep$data_internal, idname = idname, tname = tname, yname = yname, gname = gname,
    control = produced$control, xformla = prep$xformla, weightsname = NULL, clustername = NULL,
    cluster_vec = prep$cluster_full, target_cohorts = prep$target_cohorts,
    anticipation = as.integer(anticipation), control_group = gate$control_group,
    base_period = base_period, weighting = gate$gmm_weighting_requested,
    pt_assumption = "post", pt_start = NULL, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = character(0L), cband = FALSE, alp = alp, biters = as.integer(biters),
    btype = "rademacher", seed = NULL, fixed_signs = NULL,
    cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w,
    frame_g = prep$frame_g_arg, allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE,
    bstrap = FALSE, caller_rng = TRUE, custom_est_method = NULL, custom_dots = list(),
    gmm_covariance = gate$gmm_covariance, compute_inffunc = TRUE, print_details = gate$print_details,
    conditional_pretest = TRUE), "cdid_error", "request_context", context)
  src <- cpt_grid_source(data, frame, idname, tname, prep)
  route_text <- paste0("est_method = \"", est_label, "\", selection = \"", gate$selection, "\", selection_est_method = \"",
    gate$selection_est_method, "\"", if (!is.null(state$registry$cp04)) paste0("; backend law ", state$registry$cp04$law) else
      "; risk-set IPW links", if (!is.null(prep$reference_source)) "; reference-date recipes" else "")
  test <- ipw_pipeline_resignal_with_state(chained_conditional_test(state, prep, src, prep$xformla, alp, biters,
    produced$control$multiplier_budget, plan$workers, route_text, fixed_signs = fixed_signs, keep_rows = keep_rows,
    progress = progress), "cdid_error", "request_context", context)
  m <- test$meta
  m$request <- meta$request
  m$execution <- plan
  m$controls <- list(chunk_rows = produced$record$chunk_rows, numeric_workspace_entries = produced$record$numeric_workspace_entries,
    counted_work_limit = produced$record$counted_work_limit, profile_source = resolved$source)
  out <- if (identical(m$status, "available")) cpt_result(m, clustervars, xformla, CvM = test$CvM, CvMb = test$CvMb,
    CvMcval = test$CvMcval, CvMpval = test$CvMpval) else cpt_result(m, clustervars, xformla)
  if (isTRUE(keep_rows)) {
    # Private diagnostics only (small inputs): original-unit rows and the state.
    attr(out, "cells") <- test$cells
    attr(out, "state") <- state
  }
  out
}

print.cdid_MP.TEST <- function(x, ...) {
  summary.cdid_MP.TEST(x, ...)
  invisible(x)
}

summary.cdid_MP.TEST <- function(object, ...) {
  m <- object$chained
  cat("Chained conditional pretrend test (cdid, mode = \"chained\")\n")
  if (!identical(m$status, "available")) {
    cat("Status          : ", m$status, "\n", sep = "")
    cat("Reason          : ", m$reason, "\n", sep = "")
    return(invisible(object))
  }
  cat("Cramer von Mises: \n")
  cat("  Test Statistic: ", object$CvM, "\n")
  cat("  Critical Value: ", object$CvMcval, "\n")
  cat("  P-value       : ", object$CvMpval, "\n \n")
  cat("Cells (base_period = \"", m$base_period, "\"): ", paste(m$cells$cell_id, collapse = ", "), "\n", sep = "")
  cat("Multiplier draws: ", m$inference$B, " on ", m$inference$n_clusters,
    if (isTRUE(m$inference$clustered)) " clusters" else " units", "\n", sep = "")
  cat("Clustering on   : ", paste0(object$clustervars, sep = ","), "\n")
  cat("X formula       : ", paste(deparse(object$xformla), collapse = " "), "\n")
  cat("Note            : ", m$not_wald, "\n", sep = "")
  invisible(object)
}
