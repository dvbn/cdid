# ==============================================================================
# Private Chained Augmented Observation-Correction Backend (D52 seventh route)
# ==============================================================================
# selection_est_method = "dr" combined with treatment est_method = "dr" on the
# representative roster F (AUGMENTED_SELECTION_SPEC.md A1-A9 and E1-E7;
# IMPLEMENTATION_HANDOFF.md sections 1-9).
#
# Both observation-correction arms are active. For every needed actual
# cohort/pair cell c = (k, e) the reconstructed response is
#   Z_c = j_c + pi_c (Delta_e - j_c),   pi_c = S_e / max(q_c, eps),
# with Delta_e and 1/max(q_c, eps) evaluated only on observed pairs (S_e = 1);
# the selection path d log max(q_c, eps) vanishes where the floor binds. The mass is
# the pure-imputation roster mass M_c(h) = H_N[w I_k a_k h]; the change total is
# Y_c(h) = a_k H_N[w I_k Z_c h] (A1-A2). The roster treatment logits are the
# pure-imputation roots (no j or q path); the shared counterfactual OR fits Z on
# all F rows with the unchanged synthetic weights over every deduplicated block;
# the link keeps both normalizers (A3).
#
# Composition (every accepted kernel is called unchanged):
#   S0  reader/design/recipe/catalog checks. reader$state$Z is the
#       observation-probability design and is never used as the outcome design.
#   S1  fit_selection_models() on the same reader (U = F, rho = 1); the
#       selection pair-cohort catalog must equal the imputation cells C_imp.
#   S2-S3 imputation_fit_core(..., est_method = "ipw"): the accepted group-local
#       WLS imputation fits and roster treatment logits. Its pure-imputation link
#       pass is computed by that kernel and discarded (charged below).
#   S4  per-cell augmentation scan: Y_c(1), Y_c(z), M^W_c, M^Delta_c and R^q_c(1).
#   S5  shared counterfactual OR on Z with rate, j (M^Delta) and q (R^q) forcing.
#   S6  link ratios, M^Delta adjoints and the merged selection loadings.
# Replay forms the complete original-unit IF (A6-A9) in bounded row chunks:
# dir + imputation + selection + OR + PS + rate. The selection rows xi_theta are
# the accepted original-unit block rows on theta coordinates only; the repeat
# rate a_k is counted once, in the rate channel (spec 5.3).
# ==============================================================================

augmented_abort <- function(message, class = "input", extra = list()) {
  stop(structure(
    c(list(message = message, call = NULL), extra),
    class = unique(c(
      paste0("cdid_augmented_", class, "_error"),
      "cdid_augmented_error",
      paste0("cdid_cp04_", class, "_error"),
      "cdid_ipw_fit_error", "cdid_error", "error", "condition"
    ))
  ))
}

# One actual cell on rows of its cohort: fitted imputation j on every row and,
# on observed pairs only, the raw q, pi = 1/max(q, eps), e = Delta - j and
# (optionally) the selection log-derivative channel rows
# T_logit = -d log max(q, eps) / d theta in the cell's theta coordinates (zero
# where the floor binds). Missing outcomes are never read; Z = j off pairs.
augmented_cell_rows <- function(sel_fit, up_map, X_F, pair_reader, cf, gf, rows,
                                need_logderiv = FALSE) {
  n <- length(rows)
  jhat <- if (gf$p > 0L) {
    cf$alpha + as.vector((cf$L * X_F[rows, gf$kept_cols_F, drop = FALSE]) %*% gf$beta)
  } else {
    rep(cf$alpha, n)
  }
  hit <- pair_reader(rows, cf$pair, TRUE)
  obs <- which(hit$S == 1)
  out <- list(jhat = jhat, Z = jhat, obs = obs, q = numeric(0L), pi = numeric(0L),
    e = numeric(0L), T_logit = NULL)
  if (length(obs)) {
    if (need_logderiv) {
      ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, cf$q_id, rows[obs])
      q <- ld$pred$prob
      out$T_logit <- ld$T_logit
    } else {
      q <- predict_pair_selection(sel_fit, rows[obs], cf$q_id)$prob
    }
    cp04_finite(q, "fitted observation probabilities on observed pairs", TRUE)
    pi <- 1 / ipw_clip_floor(q)
    e <- hit$dy[obs] - jhat[obs]
    out$Z[obs] <- jhat[obs] + pi * e
    out$q <- q
    out$pi <- pi
    out$e <- e
  }
  out
}

# Logical payload bound of the compact fitted fields (cp04_payload convention,
# not RSS). Groups, treatments and the pure-imputation shared-OR fields reuse
# imputation_retained_core_bound(); augmented cell, block, global q-forcing and
# link fields are charged here with the unreduced p as the slope bound.
augmented_retained_core_bound <- function(cat, recipe, trt_recipe, theta_ids) {
  E <- function(x) unname(cp04_payload(x)["entries"])
  A <- function(n, nm = NULL) n + 8 + if (is.null(nm)) 0 else 9 + E(nm)
  M <- function(nr, nc, rn = NULL, cn = NULL, dn = FALSE) {
    nr * nc + 27 + if (dn) 1 + E(list(rn, cn)) else 0
  }
  R <- function(fields, cls = NULL, nm = names(fields)) {
    length(fields) + 8 + sum(fields) +
      if (is.null(nm) && is.null(cls)) 0 else
        8 + as.integer(!is.null(nm)) + as.integer(!is.null(cls)) + E(nm) + E(cls)
  }
  V <- function(x) length(x) + 8 + sum(x)
  p <- as.double(recipe$p)
  z <- recipe$slope_columns
  q <- as.double(trt_recipe$p)
  nc <- length(cat$cell_specs)
  ng <- length(cat$group_specs)
  dq <- vapply(theta_ids, length, integer(1L))
  pure <- imputation_retained_core_bound(cat, recipe, trt_recipe)

  cells <- V(vapply(seq_len(nc), function(cid) {
    cs <- cat$cell_specs[[cid]]
    E(cs) + R(c(group_id = A(1), local_id = A(1), n_pair = A(1), weight_sum = A(1), mass = A(1),
      alpha = A(1), mean_d = A(p, z), mean_y = A(1), a_k = A(1), M_one = A(1),
      Y_one = A(1), M_z = A(p, z), Y_z = A(p, z), M_d_kept = A(p, z))) +
      R(c(q_id = A(1), theta_ids = A(dq[cid]), Yaug_one = A(1), Yaug_z = A(p, z),
        MW_one = A(1), MW_d = A(p, z), MD_one = A(1), MD_d = A(p, z), Rq_one = A(dq[cid]),
        clipping = E(ipw_clip_summary(q = 0.5, weights = 1))))
  }, numeric(1L)))

  block_q <- vapply(cat$or_blocks, function(b) {
    d <- sum(dq[b$control_cell_ids])
    2 * A(d) + 64
  }, numeric(1L))
  d_global <- sum(dq[sort(unique(unlist(lapply(cat$or_blocks, `[[`, "control_cell_ids"))))])
  shared <- unname(pure["shared_or"]) + sum(block_q) + A(d_global) + M(p, d_global) + 64

  adjoint <- R(c(tilde_lambda = A(nc), tilde_nu = ng + 8 + ng * A(p)))
  links <- V(vapply(cat$links, function(lk) {
    k <- length(lk$control_cell_ids)
    dg <- dq[lk$target_cell_id]
    dc <- sum(dq[lk$control_cell_ids])
    E(lk) + R(c(mu = A(1), mu_dr = A(1),
      d_plus = A(1), num_T = A(1), T_gl = A(1), num_Bg = A(1), B_gl = A(1),
      d_minus = A(1), num_C = A(1), C_ol = A(1), num_Bo = A(1), B_ol = A(1),
      d_bar_plus_or = A(p, z), d_bar_minus_or = A(p, z), d_bar_diff_or = A(p, z),
      a_coeff_T_g = A(1), a_coeff_Bg_g = A(1), a_coeff_C_k = A(k), a_coeff_Bo_k = A(k),
      m_minus_k = A(k), lam_T = A(1), nu_T = A(p, z), lam_C_k = A(k),
      nu_C_k = k + 8 + k * A(p, z), denminus_gamma = A(q),
      g_C_gamma = A(q, trt_recipe$columns), g_Bo_gamma = A(q, trt_recipe$columns),
      g_dr_gamma = A(q, trt_recipe$columns),
      sel_T_ids = A(dg), sel_T = A(dg), sel_C_ids = A(dc), sel_C = A(dc),
      adjoint_T = adjoint, adjoint_C = adjoint))
  }, numeric(1L)))

  # The pure-imputation convention charges the shared-OR alias field twice.
  c(cells = cells, groups = unname(pure["groups"]), shared_or = shared, shared = shared,
    treatments = unname(pure["treatments"]), links = links)
}

fit_augmented_backend <- function(reader, link_views, control, glm_control, rank_tol,
                                  outcome_recipe, treatment_columns = NULL,
                                  est_method = "dr",
                                  imputation_pooling = "separate",
                                  selection_model = "pair_logit",
                                  selection_pooling = "separate",
                                  design_owner = NULL,
                                  X_frame = NULL,
                                  retained_external_entries = 0,
                                  external_work_ops = 0,
                                  automatic_chunk_rows = FALSE,
                                  pair_outcomes = NULL,
                                  reference = NULL) {
  validate_cp04_reader(reader)
  # D60: every stage runs on the recipe-cohort layout and its sealed layout
  # reader (cohort k under recipe r is one layout cohort with W_g, a_k^{g,F},
  # z(C^r) and the selection design of the same extraction date); influence
  # rows are summed back to one row per original unit.
  lay <- NULL
  if (!is.null(reference)) {
    if (!is.null(design_owner) || !is.null(X_frame) || !is.null(pair_outcomes))
      augmented_abort("A reference-input augmented fit reads its designs only through the sealed reference reader.", "identity")
    lay <- chained_reference_layout(reader, reference, link_views, validate_ipw_owner_control(control),
      retained = reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
        retained_external_entries + reference$accounting$owned_entries, work = external_work_ops, selection = TRUE)
    reader <- lay$reader
    validate_cp04_reader(reader)
  }
  if (!is.character(est_method) || is.object(est_method) || length(est_method) != 1L ||
      is.na(est_method) || !identical(est_method, "dr")) {
    augmented_abort(paste0(
      'The augmented observation correction (selection_est_method = "dr") is implemented only with ',
      'treatment est_method = "dr"; no other treatment method or observation correction is substituted.'
    ), "method")
  }
  imputation_pooling <- imputation_pooling_choice(imputation_pooling)
  validate_selection_mode_args(selection_model, selection_pooling)
  if (missing(rank_tol)) augmented_abort("Supply explicit private `rank_tol`.")
  rank_tol <- or_baseline_rank_tol(rank_tol)
  ctrl <- validate_ipw_owner_control(control)
  gc <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE)) {
    augmented_abort("Invalid external pipeline reservation.")
  }
  if (!identical(automatic_chunk_rows, isTRUE(automatic_chunk_rows))) {
    augmented_abort("Invalid private automatic row-planning flag.")
  }

  st <- reader$state
  N <- st$N
  if (!identical(st$universe, "F") || !all(st$U) || any(st$rho != 1)) {
    augmented_abort(paste0(
      'The augmented observation correction requires the full representative roster `universe = "F"` ',
      "with U = TRUE and rho = 1 on every F row."
    ), "universe")
  }
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (missing(outcome_recipe)) {
    augmented_abort("Supply the explicit private outcome recipe (intercept flag and covariate roles).")
  }

  # S0. Two designs, one reader: the outcome/treatment design reaches this
  # backend only through X_frame or a sealed design owner (or, D60, the
  # sealed reference reader's recipe-cohort layout).
  if (is.null(lay) && is.null(design_owner) && is.null(X_frame)) {
    augmented_abort(paste0(
      "The augmented observation correction never reads the outcome/treatment design from `reader$state$Z`, ",
      "which holds the observation-probability (selection_xformla) design. Supply the full-F xformla design ",
      "through `X_frame` or `design_owner`."
    ), "identity")
  }
  des <- if (!is.null(lay)) list(X = lay$X, columns = reference$x_columns, source = "reference_recipe_layout",
    design_ref = NULL, owned_entries = 0, borrowed_entries = reference$accounting$owned_entries, work_ops = 0) else
    imputation_resolve_full_design(
    reader = reader,
    frame_owner = frame,
    design_owner = design_owner,
    X_frame = X_frame,
    control = ctrl
  )
  if (!des$source %in% c("design_owner", "X_frame", "reference_recipe_layout")) {
    augmented_abort("The outcome/treatment design must come from `X_frame` or `design_owner`.", "identity")
  }
  X_F <- des$X
  f_cols <- des$columns
  or_chk <- dr_check_outcome_recipe(outcome_recipe, f_cols, "baseline")
  recipe <- or_chk$recipe
  or_frame_index <- or_chk$frame_index
  trt_recipe <- dr_check_treatment_recipe(treatment_columns, f_cols, "baseline")

  if (!is.list(link_views) || !length(link_views)) augmented_abort("Supply actual private link views.")
  L <- length(link_views)
  G <- length(st$cohorts)
  p_or_raw <- recipe$p
  p_trt <- trt_recipe$p
  slope_cols_F <- if (p_or_raw > 0L) or_frame_index[recipe$slope_index] else integer(0L)

  base_entries <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    des$owned_entries + des$borrowed_entries + retained_external_entries
  header <- 4096 + 1024 * L * (G + 1)
  cp04_budget(
    ctrl, c(base_entries, header), header,
    c(external_work_ops, des$work_ops, 128 * L * L * (G + 1)),
    header, "augmented observation-correction metadata"
  )

  records <- vector("list", L)
  for (ell in seq_len(L)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:")) {
      augmented_abort("The augmented backend requires the supplied comparison-catalog functional key.")
    }
    key_bytes <- nchar(source, type = "bytes")
    cp04_budget(
      ctrl, c(base_entries, header), header + 16 * key_bytes,
      c(external_work_ops, des$work_ops, 128 * L * L * (G + 1) + 64 * key_bytes),
      header, "augmented link functional key"
    )
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate)) {
      augmented_abort("Legacy nuisance candidates are unsupported.")
    }
    if (!is.null(v$nuisance_ref)) {
      nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair)) {
        augmented_abort("Cross-pair or unverified nuisance sharing is unsupported.")
      }
    }
    if (is.null(lay)) {
      g_pos <- match(v$cohorts$g, st$cohorts)
      k_pos <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    } else {
      g_pos <- lay$view_cohorts[[ell]][1L]; k_pos <- lay$view_cohorts[[ell]][-1L]
    }
    if (is.na(g_pos) || anyNA(k_pos) || !isTRUE(all(v$support$by_cohort_pair_count >= 1))) {
      augmented_abort("Positive actual-pair cohort support required.", "support")
    }
    or_pooled_transitions(v$pair, "baseline")
    lid <- if (!is.null(v$link_id) && nzchar(as.character(v$link_id))) as.character(v$link_id) else paste0("link_", ell)
    records[[ell]] <- list(
      pair = v$pair,
      target = as.integer(g_pos),
      controls = as.integer(k_pos),
      link_id = lid,
      source_functional_key = source,
      functional_key = paste0("cmpcat:AUG-DR:", imputation_pooling, ":", selection_model, ":",
        selection_pooling, ":", source,
        if (is.null(lay)) "" else paste0("|reference:", reference$recipes$key[reference$view_recipe[ell]]))
    )
  }
  lids <- vapply(records, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1L))
  if (any(!nzchar(lids)) || anyDuplicated(lids)) augmented_abort("Distinct nonempty link IDs required.")

  cat <- imputation_build_catalog(records, st$cohorts, "dr", imputation_pooling,
    families = if (is.null(lay)) NULL else lay$families)
  n_cells <- length(cat$cell_specs)
  n_groups <- length(cat$group_specs)
  n_blocks <- length(cat$or_blocks)
  n_trt <- length(cat$treatment_specs)

  pair_entries <- if (is.null(pair_outcomes)) 0 else unname(cp04_payload(pair_outcomes)["entries"])
  if (!is.null(pair_outcomes) && (!is.list(pair_outcomes) || is.object(pair_outcomes))) {
    augmented_abort("`pair_outcomes` must be a named private list of full-F differences.")
  }
  base_entries <- base_entries + pair_entries

  # S1. Observation-probability models on the SAME reader. Refusals propagate
  # unchanged; nothing falls back to pure imputation or to weighting.
  sel_fit <- fit_selection_models(
    reader = reader,
    link_views = link_views,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    control = ctrl,
    glm_control = gc,
    retained_external_entries = des$owned_entries + des$borrowed_entries +
      retained_external_entries + pair_entries + header,
    external_work_ops = external_work_ops + des$work_ops,
    view_cohorts = if (is.null(lay)) NULL else lay$view_cohorts,
    families = if (is.null(lay)) NULL else lay$families
  )

  # Catalog identity: the selection pair-cohort catalog is C_imp itself.
  q_of_cell <- vapply(cat$cell_specs, function(cs) {
    tryCatch(
      as.integer(resolve_pair_cohort_id(sel_fit, list(pair = cs$pair, cohort = cs$cohort_index))),
      cdid_cp04_input_error = function(e) NA_integer_
    )
  }, integer(1L))
  if (anyNA(q_of_cell) || anyDuplicated(q_of_cell) ||
      length(sel_fit$pair_cohorts) != n_cells || !setequal(q_of_cell, seq_len(n_cells))) {
    augmented_abort(paste0(
      "The actual imputation cells and the fitted observation-probability pair-cohort catalog are not ",
      "one-to-one; no cell is reconstructed without its own fitted probability."
    ), "identity")
  }
  up_map <- build_selection_upstream_map(sel_fit)
  theta_ids <- lapply(q_of_cell, function(q) as.integer(up_map$pair_channels[[q]][-1L]))
  d_q <- vapply(theta_ids, length, integer(1L))

  # Measure the metadata skeleton, then charge every future compact field.
  obj <- structure(
    list(
      selection_est_method = "dr",
      est_method = "dr",
      dr_method = "dr",
      imputation_pooling = imputation_pooling,
      selection = "conditional_mtar",
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      uses_selection_models = TRUE,
      universe = "F",
      law = paste0(
        "AUG-DR:", imputation_pooling, ":", selection_model, ":", selection_pooling, ":",
        if (recipe$intercept) "local-intercepts" else "no-intercept",
        if (is.null(lay)) ":actual-pair:fixed-w:full-F" else ":actual-pair:W_g:full-F:recipe-cohort-layout"
      ),
      outcome_recipe = recipe,
      or_frame_index = as.integer(or_frame_index),
      treatment_recipe = trt_recipe,
      selection_design = sel_fit$design_spec,
      pair_outcomes = pair_outcomes,
      design_source = des$source,
      design_columns = f_cols,
      selection_columns = colnames(st$Z),
      X_F = if (is.null(lay)) NULL else des$X, # otherwise retained by the sealed design owner
      layout = if (is.null(lay)) NULL else lay[c("orig_row", "by_orig", "cohort", "recipe", "families",
        "N", "N_layout", "frame_id")],
      pairs = cat$pairs,
      cell_specs = cat$cell_specs,
      cell_lookup = cat$cell_lookup,
      active_cohorts = cat$active_cohorts,
      group_specs = cat$group_specs,
      or_blocks = cat$or_blocks,
      treatment_specs = cat$treatment_specs,
      selection_backend = sel_fit,
      up_map = up_map,
      q_of_cell = as.integer(q_of_cell),
      theta_ids = theta_ids,
      cells = NULL,
      groups = NULL,
      shared_or = NULL,
      shared = NULL,
      treatments = NULL,
      links = NULL,
      cohorts = st$cohorts,
      a = st$a,
      tau = st$tau,
      eta = st$eta,
      lambda = st$lambda,
      rank_tol = rank_tol,
      glm_control = gc,
      automatic_chunk_rows = automatic_chunk_rows,
      reader_ref = reader$ref,
      design_ref = des$design_ref,
      frame_ref = reader$frame_ref,
      panel_ref = reader$panel_ref,
      control = ctrl,
      accounting = list(
        owned_entries = 0,
        character_bytes = 0,
        borrowed_entries = base_entries - pair_entries,
        retained_external_entries = retained_external_entries,
        design_owned_entries = des$owned_entries,
        design_borrowed_entries = des$borrowed_entries,
        external_work_ops = external_work_ops,
        retained_bound = 0,
        fit_scratch_bound = 0,
        work_ops_bound = 0,
        selection_work_ops_bound = sel_fit$accounting$work_ops_bound,
        automatic_chunk_rows = automatic_chunk_rows,
        n_cells = n_cells,
        n_groups = n_groups,
        n_imputation_fits = n_groups,
        n_imputation_local_intercepts = if (recipe$intercept) n_cells else 0L,
        n_imputation_slopes_by_group = rep.int(as.integer(p_or_raw), n_groups),
        n_or_fits = 1L,
        n_or_blocks = n_blocks,
        n_or_local_intercepts = if (recipe$intercept) n_blocks else 0L,
        n_or_shared_slopes = as.integer(p_or_raw),
        n_treatment_models = n_trt,
        n_treatment_fits = n_trt,
        n_treatment_columns = p_trt,
        n_selection_fits = length(sel_fit$blocks),
        n_selection_cells = length(sel_fit$cells),
        n_selection_blocks = length(sel_fit$blocks),
        n_selection_parameters = length(sel_fit$coefficients),
        n_selection_pair_cohorts = length(sel_fit$pair_cohorts),
        n_selection_design_columns = ncol(st$Z),
        observation_logit_models = length(sel_fit$blocks),
        selection_owned_entries = sel_fit$accounting$owned_entries,
        augmentation_passes = 3L,
        full_F_masks_per_model = 0L,
        full_F_prediction_matrices = 0L,
        full_F_by_link_IF_retained = FALSE,
        inverse_matrices_retained = FALSE,
        dense_dummy_regression = FALSE,
        peak_rss_measured = FALSE
      )
    ),
    class = c("cdid_augmented_backend", "cdid_chained_augmented")
  )
  if (!is.null(lay)) {
    obj$reference_ref <- reference$ref
    obj$accounting$n_reference_recipes <- nrow(reference$recipes)
    obj$accounting$n_layout_rows <- lay$N_layout
  }
  obj["ref"] <- list(NULL)

  # Reference scans visit every recipe-cohort layout row; N remains the
  # original-unit normalization count in every estimating equation.
  N_scan <- if (is.null(lay)) N else lay$N_layout
  B_chunk <- min(N_scan, ctrl$max_chunk_rows)
  retained_bound <- unname(cp04_payload(obj)["entries"]) - pair_entries +
    sum(augmented_retained_core_bound(cat, recipe, trt_recipe, theta_ids))
  d_X <- ncol(X_F)
  d_q_max <- max(c(0L, d_q))
  fit_scratch <- max(
    64 * B_chunk * (p_or_raw + p_trt + G + 8) + 64 * N_scan * (p_trt + 4) +
      128 * ((p_or_raw + 2)^2 * max(1L, n_cells) + (p_trt + 2)^2) + 8192 +
      64 * B_chunk * (d_X + ncol(st$Z) + 2 * d_q_max + 2 * p_or_raw + 16) +
      128 * (p_or_raw + 2) * (d_q_max + 2) * max(1L, n_cells),
    sel_fit$accounting$fit_scratch_bound
  )
  log_n <- 1 + ceiling(log2(N_scan + 1))
  n_k <- function(ci) st$cohort_end[ci] - st$cohort_start[ci] + 1
  scan_cells <- sum(vapply(cat$cell_specs, function(cs) n_k(cs$cohort_index), numeric(1L)))
  scan_blocks <- sum(vapply(cat$or_blocks, function(b)
    sum(vapply(b$control_cohort_indices, n_k, numeric(1L))), numeric(1L)))
  scan_links <- sum(vapply(cat$links, function(lk) n_k(lk$target_cohort_index) +
    sum(vapply(lk$control_cohort_indices, n_k, numeric(1L))), numeric(1L)))
  # Stage 1-5 of imputation_fit_core (including its discarded link pass), the
  # three augmentation scans, the augmented link pass and the M^Delta adjoints.
  pure_work <- 256 * N_scan * (
    (n_cells + n_blocks + 1) * (p_or_raw + 2)^2 * log_n +
      n_trt * (p_trt + 2)^2 * gc$maxit +
      L * (p_or_raw + p_trt + G + 4)
  )
  aug_work <- 256 * (2 * scan_cells + 2 * scan_blocks + scan_links) *
    (d_X + ncol(st$Z) + d_q_max + log_n + p_trt + (p_or_raw + 2) * (p_or_raw + d_q_max + 2)) +
    256 * L * (n_groups + 1) * (p_or_raw + 2)^3
  fit_work <- pure_work + aug_work + 128 * retained_bound

  obj$accounting$retained_bound <- base_entries + retained_bound
  obj$accounting$fit_scratch_bound <- fit_scratch
  obj$accounting$work_ops_bound <- external_work_ops + des$work_ops +
    sel_fit$accounting$work_ops_bound + fit_work
  cp04_budget(
    ctrl, c(base_entries, retained_bound), fit_scratch,
    c(external_work_ops, des$work_ops, sel_fit$accounting$work_ops_bound, fit_work),
    retained_bound + pair_entries, "augmented observation-correction fitting preflight"
  )

  pair_reader <- function(rows, pair, read_outcomes = TRUE) {
    pr <- cp04_pair_rows(reader, rows, pair, read_outcomes)
    if (read_outcomes && !is.null(pair_outcomes)) {
      pk <- paste0(as.character(pair$s), ":", as.character(pair$t))
      dy_vec <- pair_outcomes[[pk]]
      if (is.null(dy_vec) || !is.numeric(dy_vec) || length(dy_vec) != N) {
        augmented_abort(paste0("Invalid supplied `pair_outcomes` entry for pair ", pk, "."))
      }
      dy_sub <- as.double(dy_vec[rows])
      obs_sub <- !is.na(dy_sub)
      if (!identical(obs_sub, as.logical(pr$S))) {
        augmented_abort(paste0("Supplied `pair_outcomes` mask disagrees with panel observed pair indicator for ", pk, "."),
          "support")
      }
      if (any(!is.finite(dy_sub[obs_sub]))) {
        augmented_abort(paste0("Non-finite observed outcome difference in `pair_outcomes` for ", pk, "."))
      }
      pr$dy <- numeric(length(rows))
      pr$dy[obs_sub] <- dy_sub[obs_sub]
    }
    pr
  }
  # Bounded cohort row chunks: one index vector of at most B_chunk rows is live.
  chunk_starts <- function(ci) seq.int(st$cohort_start[ci], st$cohort_end[ci], by = B_chunk)
  chunk_rows <- function(ci, start) {
    st$cohort_order[seq.int(start, min(st$cohort_end[ci], start + B_chunk - 1L))]
  }

  # S2-S3. Accepted group-local imputation WLS and roster treatment logits.
  core <- imputation_fit_core(
    st = st,
    X_F = X_F,
    pair_reader = pair_reader,
    cat = cat,
    recipe = recipe,
    or_frame_index = or_frame_index,
    trt_recipe = trt_recipe,
    est_method = "ipw",
    imputation_pooling = imputation_pooling,
    glm_control = gc,
    rank_tol = rank_tol,
    ctrl = ctrl
  )
  groups <- core$groups
  treatments <- core$treatments

  # Full-F cohort column maxima for the shared-OR alias verification (same rule
  # as the pure-imputation OR stage).
  cohort_colmax_F <- vector("list", G)
  for (ci in cat$active_cohorts) {
    cmax <- numeric(p_or_raw)
    if (p_or_raw > 0L) for (start in chunk_starts(ci)) {
      Z_sl <- abs(X_F[chunk_rows(ci, start), slope_cols_F, drop = FALSE])
      if (nrow(Z_sl) > 0L) cmax <- pmax(cmax, apply(Z_sl, 2L, max))
    }
    cohort_colmax_F[[ci]] <- cmax
  }

  # S4. Augmentation scan per cell (link independent).
  cells <- vector("list", n_cells)
  for (cid in seq_len(n_cells)) {
    cf <- core$cells[[cid]]
    gf <- groups[[cf$group_id]]
    cf$q_id <- q_of_cell[cid]
    cf$theta_ids <- theta_ids[[cid]]
    ci <- cf$cohort_index
    a_k <- st$a[ci]
    Y_aug <- 0
    Y_aug_z <- numeric(p_or_raw)
    MW_one <- 0
    MW_d <- numeric(gf$p)
    Rq_one <- numeric(d_q[cid])
    n_obs <- 0
    clip_c <- ipw_clip_summary(q = numeric(0), weights = numeric(0))
    for (start in chunk_starts(ci)) {
      rows <- chunk_rows(ci, start)
      ev <- augmented_cell_rows(sel_fit, up_map, X_F, pair_reader, cf, gf, rows, TRUE)
      w_syn <- st$w[rows] * a_k
      Y_aug <- Y_aug + sum(w_syn * ev$Z)
      if (p_or_raw > 0L) {
        Y_aug_z <- Y_aug_z + as.vector(crossprod(w_syn * ev$Z, X_F[rows, slope_cols_F, drop = FALSE]))
      }
      ob <- ev$obs
      if (length(ob)) {
        w_pi <- w_syn[ob] * ev$pi
        clip_c <- ipw_clip_summary_merge(clip_c, ipw_clip_summary(q = ev$q, weights = w_pi))
        MW_one <- MW_one + sum(w_pi)
        if (gf$p > 0L) {
          MW_d <- MW_d + as.vector(crossprod(w_pi, cf$L * X_F[rows[ob], gf$kept_cols_F, drop = FALSE]))
        }
        if (d_q[cid] > 0L) Rq_one <- Rq_one - as.vector(crossprod(w_pi * ev$e, ev$T_logit))
        n_obs <- n_obs + length(ob)
      }
    }
    if (n_obs != cf$n_pair) {
      augmented_abort(paste0("Observed-pair support of cell `", cf$name,
        "` changed between the imputation fit and the augmentation scan."), "integrity")
    }
    Y_aug <- Y_aug / N
    Y_aug_z <- Y_aug_z / N
    MW_one <- MW_one / N
    MW_d <- MW_d / N
    Rq_one <- Rq_one / N
    names(Y_aug_z) <- recipe$slope_columns
    names(MW_d) <- gf$kept_columns
    MD_one <- cf$M_one - MW_one
    MD_d <- as.double(cf$M_d_kept) - MW_d
    names(MD_d) <- gf$kept_columns
    cp04_finite(c(Y_aug, Y_aug_z, MW_one, MW_d, MD_one, MD_d, Rq_one), "augmented cell moments")
    cells[[cid]] <- c(cf, list(
      Yaug_one = Y_aug,
      Yaug_z = Y_aug_z,
      MW_one = MW_one,
      MW_d = MW_d,
      MD_one = MD_one,
      MD_d = MD_d,
      Rq_one = Rq_one,
      clipping = clip_c
    ))
  }

  # S5. Shared counterfactual OR on Z over every deduplicated block (roster
  # synthetic weights w a_k on all F rows of each control cohort).
  used_ctrl <- sort(unique(unlist(lapply(cat$or_blocks, `[[`, "control_cell_ids"))))
  ctrl_acc <- vector("list", n_cells)
  for (cid in used_ctrl) {
    cf <- cells[[cid]]
    gf <- groups[[cf$group_id]]
    a_k <- st$a[cf$cohort_index]
    acc <- or_pooled_support_start(recipe, cf$L)
    for (start in chunk_starts(cf$cohort_index)) {
      rows <- chunk_rows(cf$cohort_index, start)
      ev <- augmented_cell_rows(sel_fit, up_map, X_F, pair_reader, cf, gf, rows, FALSE)
      acc <- or_pooled_support_add(acc, recipe, X_F[rows, or_frame_index, drop = FALSE], ev$Z,
        st$w[rows] * a_k, "augmented_or")
    }
    ctrl_acc[[cid]] <- acc
  }
  supports_or <- lapply(seq_len(n_cells), function(cid) {
    acc <- ctrl_acc[[cid]]
    list(summary = if (!is.null(acc)) acc$summary else NULL, scale = 1.0,
      colmax = cohort_colmax_F[[cells[[cid]]$cohort_index]])
  })
  blocks_or_spec <- lapply(cat$or_blocks, function(b) {
    list(name = b$name, support_ids = as.integer(b$control_cell_ids), L = as.integer(b$L))
  })
  links_or_spec <- lapply(cat$links, function(lk) {
    list(regression = as.integer(lk$or_block_id), support_id = as.integer(lk$target_cell_id))
  })
  visit_or_rows <- function(sid, visit) {
    ci_s <- cells[[sid]]$cohort_index
    for (start in chunk_starts(ci_s)) visit(X_F[chunk_rows(ci_s, start), or_frame_index, drop = FALSE])
    invisible(NULL)
  }
  fit_or <- tryCatch(or_pooled_fit_blocks(
    recipe = recipe,
    supports = supports_or,
    blocks = blocks_or_spec,
    links = links_or_spec,
    rank_tol = rank_tol,
    N = N,
    visit_rows = visit_or_rows,
    route = "augmented_or"
  ), cdid_or_wls_rank_error = function(e) {
    e$message <- paste0(conditionMessage(e), " This computation failure does not establish population ATT nonidentification.")
    e$is_identifiability_claim <- FALSE
    stop(e)
  })
  keep_or <- fit_or$keep
  p_or <- length(keep_or)
  kept_or_cols_F <- if (p_or > 0L) slope_cols_F[keep_or] else integer(0L)
  kept_or_names <- if (p_or > 0L) recipe$slope_columns[keep_or] else character(0L)
  beta_or <- as.double(fit_or$beta)
  names(beta_or) <- kept_or_names
  alpha_or <- if (recipe$intercept) {
    vapply(fit_or$regressions, `[[`, numeric(1L), "alpha")
  } else {
    numeric(n_blocks)
  }
  cp04_finite(c(alpha_or, beta_or), "augmented shared counterfactual OR coefficients")

  # Forcing of every block constituent: rate (with r = Z - b), j (M^Delta) and
  # q (-R^q). The global shared-slope q forcing collects all blocks once.
  q_beta_ids <- as.integer(sort(unique(c(integer(0L), unlist(theta_ids[used_ctrl])))))
  q_beta <- matrix(0.0, p_or, length(q_beta_ids))
  or_block_fits <- vector("list", n_blocks)
  for (bid in seq_len(n_blocks)) {
    b_spec <- cat$or_blocks[[bid]]
    reg_b <- fit_or$regressions[[bid]]
    cids_b <- b_spec$control_cell_ids
    w_sum_b <- sum(vapply(cids_b, function(cid) ctrl_acc[[cid]]$weight_mass, numeric(1L)))
    mass_b <- w_sum_b / N
    if (!is.finite(mass_b) || mass_b <= 0) {
      augmented_abort(paste0("Shared OR block `", b_spec$name, "` has nonpositive synthetic control mass."), "support")
    }
    mean_d_b <- if (recipe$intercept && p_or > 0L) as.double(reg_b$mean_d) else numeric(p_or)
    names(mean_d_b) <- kept_or_names
    alpha_b <- if (recipe$intercept) as.double(reg_b$alpha) else 0.0
    L_b <- b_spec$L
    n_kb <- length(cids_b)
    a_force_alpha <- numeric(n_kb)
    a_force_beta <- matrix(0.0, p_or, n_kb)
    j_force_alpha_alpha <- numeric(n_kb)
    j_force_alpha_beta <- vector("list", n_kb)
    j_force_beta_alpha <- matrix(0.0, p_or, n_kb)
    j_force_beta_beta <- vector("list", n_kb)
    q_alpha_ids <- as.integer(sort(unique(c(integer(0L), unlist(theta_ids[cids_b])))))
    q_alpha <- numeric(length(q_alpha_ids))

    for (idx_k in seq_len(n_kb)) {
      cid <- cids_b[idx_k]
      cf <- cells[[cid]]
      gf <- groups[[cf$group_id]]
      a_k <- st$a[cf$cohort_index]
      sum_w_r <- 0
      sum_w_r_dtilde <- numeric(p_or)
      sum_wD_dtilde <- numeric(p_or)
      sum_wD_dtilde_dimp <- matrix(0.0, p_or, gf$p)
      Rq_dtilde <- matrix(0.0, p_or, d_q[cid])
      for (start in chunk_starts(cf$cohort_index)) {
        rows <- chunk_rows(cf$cohort_index, start)
        n_r <- length(rows)
        ev <- augmented_cell_rows(sel_fit, up_map, X_F, pair_reader, cf, gf, rows, p_or > 0L)
        w_raw <- st$w[rows]
        w_syn <- w_raw * a_k
        d_or <- if (p_or > 0L) L_b * X_F[rows, kept_or_cols_F, drop = FALSE] else matrix(0.0, n_r, 0L)
        mhat <- if (p_or > 0L) alpha_b + as.vector(d_or %*% beta_or) else rep(alpha_b, n_r)
        r_or <- ev$Z - mhat
        sum_w_r <- sum_w_r + sum(w_raw * r_or)
        if (p_or > 0L) {
          one_minus_pi <- rep(1, n_r)
          if (length(ev$obs)) one_minus_pi[ev$obs] <- 1 - ev$pi
          w_D <- w_syn * one_minus_pi
          d_tilde <- if (recipe$intercept) sweep(d_or, 2L, mean_d_b, FUN = "-") else d_or
          sum_w_r_dtilde <- sum_w_r_dtilde + as.vector(crossprod(w_raw * r_or, d_tilde))
          sum_wD_dtilde <- sum_wD_dtilde + as.vector(crossprod(w_D, d_tilde))
          if (gf$p > 0L) {
            d_imp <- L_b * X_F[rows, gf$kept_cols_F, drop = FALSE]
            sum_wD_dtilde_dimp <- sum_wD_dtilde_dimp + crossprod(w_D * d_tilde, d_imp)
          }
          if (length(ev$obs) && d_q[cid] > 0L) {
            Rq_dtilde <- Rq_dtilde - crossprod(
              d_tilde[ev$obs, , drop = FALSE] * (w_syn[ev$obs] * ev$pi * ev$e), ev$T_logit)
          }
        }
      }
      a_force_alpha[idx_k] <- sum_w_r / N
      a_force_beta[, idx_k] <- sum_w_r_dtilde / N
      j_force_alpha_alpha[idx_k] <- if (recipe$intercept) cf$MD_one else 0.0
      j_force_alpha_beta[[idx_k]] <- as.double(cf$MD_d)
      j_force_beta_alpha[, idx_k] <- if (recipe$intercept) (sum_wD_dtilde / N) else numeric(p_or)
      j_force_beta_beta[[idx_k]] <- sum_wD_dtilde_dimp / N
      if (d_q[cid] > 0L) {
        pos <- match(cf$theta_ids, q_alpha_ids)
        q_alpha[pos] <- q_alpha[pos] - cf$Rq_one
        if (p_or > 0L) {
          gpos <- match(cf$theta_ids, q_beta_ids)
          q_beta[, gpos] <- q_beta[, gpos, drop = FALSE] - Rq_dtilde / N
        }
      }
    }
    cp04_finite(c(a_force_alpha, a_force_beta, j_force_alpha_alpha, j_force_beta_alpha,
      unlist(j_force_alpha_beta), unlist(j_force_beta_beta), q_alpha), "augmented OR forcing")
    or_block_fits[[bid]] <- c(
      b_spec,
      list(
        alpha = alpha_b,
        mass = mass_b,
        weight_sum = w_sum_b,
        mean_d = mean_d_b,
        a_force_alpha = a_force_alpha,
        a_force_beta = a_force_beta,
        j_force_alpha_alpha = j_force_alpha_alpha,
        j_force_alpha_beta = j_force_alpha_beta,
        j_force_beta_alpha = j_force_beta_alpha,
        j_force_beta_beta = j_force_beta_beta,
        q_alpha_ids = q_alpha_ids,
        q_alpha = q_alpha
      )
    )
  }
  cp04_finite(q_beta, "augmented global shared-slope selection forcing")
  theta_or <- if (recipe$intercept) c(alpha_or, beta_or) else beta_or
  names(theta_or) <- paste0("OR.theta", seq.int(0L, length(theta_or) - 1L))
  shared_or <- list(
    keep = as.integer(keep_or),
    p = as.integer(p_or),
    kept_columns = kept_or_names,
    kept_cols_F = as.integer(kept_or_cols_F),
    intercept = isTRUE(recipe$intercept),
    n_blocks = as.integer(n_blocks),
    alpha = alpha_or,
    beta = beta_or,
    theta = theta_or,
    blocks = or_block_fits,
    q_beta_ids = q_beta_ids,
    q_beta = q_beta,
    fit = fit_or$fit,
    summary = fit_or$summary,
    basis = fit_or$basis,
    rss = fit_or$rss,
    objective = fit_or$objective
  )

  # S6. Links: augmented T/C numerators, unchanged normalizers, M^Delta adjoint
  # weights and merged selection loadings over the theta coordinates.
  solve_adjoint <- function(lambda_vec, nu_list) {
    tilde_lambda <- numeric(n_cells)
    tilde_nu <- vector("list", n_groups)
    for (gid in seq_len(n_groups)) {
      gf <- groups[[gid]]
      rhs_r <- numeric(gf$p)
      for (cid in gf$cell_ids) {
        cf <- cells[[cid]]
        lam_c <- if (recipe$intercept) lambda_vec[cid] else 0.0
        if (recipe$intercept) tilde_lambda[cid] <- lam_c / cf$mass
        if (gf$p > 0L) rhs_r <- rhs_r + nu_list[[cid]] - lam_c * cf$mean_d
      }
      if (gf$p > 0L) {
        rhs_mat <- matrix(rhs_r, nrow = gf$p, ncol = 1L)
        sol <- if (recipe$intercept) {
          or_pooled_wls_gram_solve(gf$fit, rhs_mat, normalizer = N)
        } else {
          or_wls_gram_solve(gf$fit, rhs_mat, normalizer = N)
        }
        tilde_nu[[gid]] <- as.double(sol[, 1L])
      } else {
        tilde_nu[[gid]] <- numeric(0L)
      }
    }
    list(tilde_lambda = tilde_lambda, tilde_nu = tilde_nu)
  }
  zero_loadings <- function() lapply(cells, function(cf) numeric(groups[[cf$group_id]]$p))

  links <- vector("list", L)
  trt_cols_F <- trt_recipe$index
  for (ell in seq_len(L)) {
    lk <- cat$links[[ell]]
    cf_g <- cells[[lk$target_cell_id]]
    a_g <- st$a[lk$target_cohort_index]
    d_plus <- cf_g$M_one
    if (!is.finite(d_plus) || d_plus <= 0) {
      augmented_abort(paste0("Link `", lk$link_id, "` has nonpositive treated denominator d_{l,+}."), "support")
    }
    num_T <- cf_g$Yaug_one
    T_gl <- num_T / d_plus
    # Identically zero treated-rate loading (E4), computed explicitly.
    a_coeff_T_g <- ((cf_g$Yaug_one / a_g) - T_gl * (cf_g$M_one / a_g)) / d_plus
    b_fit <- or_block_fits[[lk$or_block_id]]
    d_bar_plus_or <- if (p_or > 0L) (lk$L * cf_g$M_z[keep_or]) / d_plus else numeric(0L)
    names(d_bar_plus_or) <- kept_or_names
    B_gl <- (if (recipe$intercept) b_fit$alpha else 0.0) +
      (if (p_or > 0L) sum(d_bar_plus_or * beta_or) else 0.0)
    num_Bg <- d_plus * B_gl
    a_coeff_Bg_g <- ((num_Bg / a_g) - B_gl * (d_plus / a_g)) / d_plus
    lam_T <- cf_g$MD_one / d_plus
    nu_T <- cf_g$MD_d / d_plus
    sel_T_ids <- cf_g$theta_ids
    sel_T <- -cf_g$Rq_one / d_plus
    cp04_finite(c(d_plus, num_T, T_gl, a_coeff_T_g, d_bar_plus_or, B_gl, lam_T, nu_T, sel_T),
      "augmented link treated summary")

    tfit <- treatments[[lk$trt_id]]
    n_ctrl <- length(lk$control_cohort_indices)
    sum_wm <- 0
    sum_wm_Z <- 0
    sum_wm_m <- 0
    sum_wm_dor <- numeric(p_or)
    sum_wm_by_k <- numeric(n_ctrl)
    sum_w_odds_Z_by_k <- numeric(n_ctrl)
    sum_w_odds_by_k <- numeric(n_ctrl)
    sum_w_odds_m_by_k <- numeric(n_ctrl)
    sum_wdot_Z_X <- numeric(p_trt)
    sum_wdot_X <- numeric(p_trt)
    sum_wdot_m_X <- numeric(p_trt)
    lam_C_k <- numeric(n_ctrl)
    nu_C_k <- vector("list", n_ctrl)
    sel_C_ids <- as.integer(sort(unique(c(integer(0L), unlist(lapply(lk$control_cell_ids,
      function(cid) cells[[cid]]$theta_ids))))))
    sel_C <- numeric(length(sel_C_ids))
    for (idx_k in seq_len(n_ctrl)) {
      ci_k <- lk$control_cohort_indices[idx_k]
      cid_k <- lk$control_cell_ids[idx_k]
      cf_k <- cells[[cid_k]]
      gf_k <- groups[[cf_k$group_id]]
      a_k <- st$a[ci_k]
      sum_wmD_k <- 0
      sum_dimpD_k <- numeric(gf_k$p)
      Rq_o_k <- numeric(d_q[cid_k])
      for (start in chunk_starts(ci_k)) {
        rows <- chunk_rows(ci_k, start)
        n_r <- length(rows)
        w_raw <- st$w[rows]
        w_syn <- w_raw * a_k
        X_trt_rows <- X_F[rows, trt_cols_F, drop = FALSE]
        lv_t <- ipw_logit_values(X_trt_rows, tfit$coefficients)
        odds <- lv_t$odds
        wm <- w_syn * odds
        # Odds derivative rows: wm on active rows, zero where the cap binds.
        wdot <- wm * lv_t$active
        ev <- augmented_cell_rows(sel_fit, up_map, X_F, pair_reader, cf_k, gf_k, rows, TRUE)
        Z_rows <- ev$Z
        one_minus_pi <- rep(1, n_r)
        if (length(ev$obs)) one_minus_pi[ev$obs] <- 1 - ev$pi
        wmD <- wm * one_minus_pi
        d_or <- if (p_or > 0L) lk$L * X_F[rows, kept_or_cols_F, drop = FALSE] else matrix(0.0, n_r, 0L)
        mhat <- if (p_or > 0L) {
          (if (recipe$intercept) b_fit$alpha else 0.0) + as.vector(d_or %*% beta_or)
        } else {
          rep(if (recipe$intercept) b_fit$alpha else 0.0, n_r)
        }
        sum_wm <- sum_wm + sum(wm)
        sum_wm_Z <- sum_wm_Z + sum(wm * Z_rows)
        sum_wm_by_k[idx_k] <- sum_wm_by_k[idx_k] + sum(wm)
        sum_wmD_k <- sum_wmD_k + sum(wmD)
        if (gf_k$p > 0L) {
          sum_dimpD_k <- sum_dimpD_k +
            as.vector(crossprod(wmD, lk$L * X_F[rows, gf_k$kept_cols_F, drop = FALSE]))
        }
        sum_w_odds_Z_by_k[idx_k] <- sum_w_odds_Z_by_k[idx_k] + sum(w_raw * odds * Z_rows)
        sum_w_odds_by_k[idx_k] <- sum_w_odds_by_k[idx_k] + sum(w_raw * odds)
        sum_wdot_Z_X <- sum_wdot_Z_X + as.vector(crossprod(wdot * Z_rows, X_trt_rows))
        sum_wdot_X <- sum_wdot_X + as.vector(crossprod(wdot, X_trt_rows))
        sum_wm_m <- sum_wm_m + sum(wm * mhat)
        if (p_or > 0L) sum_wm_dor <- sum_wm_dor + as.vector(crossprod(wm, d_or))
        sum_w_odds_m_by_k[idx_k] <- sum_w_odds_m_by_k[idx_k] + sum(w_raw * odds * mhat)
        sum_wdot_m_X <- sum_wdot_m_X + as.vector(crossprod(wdot * mhat, X_trt_rows))
        if (length(ev$obs) && d_q[cid_k] > 0L) {
          Rq_o_k <- Rq_o_k - as.vector(crossprod(wm[ev$obs] * ev$pi * ev$e, ev$T_logit))
        }
      }
      lam_C_k[idx_k] <- sum_wmD_k
      nu_C_k[[idx_k]] <- sum_dimpD_k
      if (d_q[cid_k] > 0L) {
        pos <- match(cf_k$theta_ids, sel_C_ids)
        sel_C[pos] <- sel_C[pos] - Rq_o_k
      }
    }
    d_minus <- sum_wm / N
    if (!is.finite(d_minus) || d_minus <= 0) {
      augmented_abort(paste0("Link `", lk$link_id, "` has nonpositive control denominator d_{l,-}."), "support")
    }
    num_C <- sum_wm_Z / N
    C_ol <- num_C / d_minus
    m_minus_k <- (sum_wm_by_k / N) / d_minus
    lam_C_k <- (lam_C_k / N) / d_minus
    nu_C_k <- lapply(nu_C_k, function(v) (v / N) / d_minus)
    sel_C <- (sel_C / N) / d_minus
    a_coeff_C_k <- ((sum_w_odds_Z_by_k / N) - C_ol * (sum_w_odds_by_k / N)) / d_minus
    denminus_gamma <- sum_wdot_X / N
    g_C_gamma <- ((sum_wdot_Z_X / N) - C_ol * (sum_wdot_X / N)) / d_minus
    names(g_C_gamma) <- trt_recipe$columns
    d_bar_minus_or <- if (p_or > 0L) (sum_wm_dor / N) / d_minus else numeric(0L)
    names(d_bar_minus_or) <- kept_or_names
    d_bar_diff_or <- d_bar_plus_or - d_bar_minus_or
    num_Bo <- sum_wm_m / N
    B_ol <- (if (recipe$intercept) b_fit$alpha else 0.0) +
      (if (p_or > 0L) sum(d_bar_minus_or * beta_or) else 0.0)
    mu_dr <- (T_gl - B_gl) - (C_ol - B_ol)
    a_coeff_Bo_k <- ((sum_w_odds_m_by_k / N) - B_ol * (sum_w_odds_by_k / N)) / d_minus
    g_Bo_gamma <- ((sum_wdot_m_X / N) - B_ol * (sum_wdot_X / N)) / d_minus
    g_dr_gamma <- g_C_gamma - g_Bo_gamma
    names(g_Bo_gamma) <- names(g_dr_gamma) <- trt_recipe$columns
    cp04_finite(c(d_minus, num_C, C_ol, m_minus_k, lam_C_k, unlist(nu_C_k), sel_C, a_coeff_C_k,
      g_C_gamma, num_Bo, B_ol, mu_dr, d_bar_minus_or, a_coeff_Bo_k, g_Bo_gamma, g_dr_gamma),
      "augmented link control and DR summary")

    lam_T_vec <- numeric(n_cells)
    nu_T_list <- zero_loadings()
    if (recipe$intercept) lam_T_vec[lk$target_cell_id] <- lam_T
    nu_T_list[[lk$target_cell_id]] <- as.double(nu_T)
    lam_C_vec <- numeric(n_cells)
    nu_C_list <- zero_loadings()
    for (idx_k in seq_len(n_ctrl)) {
      cid_k <- lk$control_cell_ids[idx_k]
      if (recipe$intercept) lam_C_vec[cid_k] <- lam_C_vec[cid_k] + lam_C_k[idx_k]
      nu_C_list[[cid_k]] <- nu_C_list[[cid_k]] + nu_C_k[[idx_k]]
    }

    links[[ell]] <- c(
      lk,
      list(
        mu = mu_dr,
        mu_dr = mu_dr,
        d_plus = d_plus,
        num_T = num_T,
        T_gl = T_gl,
        num_Bg = num_Bg,
        B_gl = B_gl,
        d_minus = d_minus,
        num_C = num_C,
        C_ol = C_ol,
        num_Bo = num_Bo,
        B_ol = B_ol,
        d_bar_plus_or = d_bar_plus_or,
        d_bar_minus_or = d_bar_minus_or,
        d_bar_diff_or = d_bar_diff_or,
        a_coeff_T_g = a_coeff_T_g,
        a_coeff_Bg_g = a_coeff_Bg_g,
        a_coeff_C_k = a_coeff_C_k,
        a_coeff_Bo_k = a_coeff_Bo_k,
        m_minus_k = m_minus_k,
        lam_T = lam_T,
        nu_T = nu_T,
        lam_C_k = lam_C_k,
        nu_C_k = nu_C_k,
        denminus_gamma = denminus_gamma,
        g_C_gamma = g_C_gamma,
        g_Bo_gamma = g_Bo_gamma,
        g_dr_gamma = g_dr_gamma,
        sel_T_ids = as.integer(sel_T_ids),
        sel_T = sel_T,
        sel_C_ids = sel_C_ids,
        sel_C = sel_C,
        adjoint_T = solve_adjoint(lam_T_vec, nu_T_list),
        adjoint_C = solve_adjoint(lam_C_vec, nu_C_list)
      )
    )
  }

  obj[c("cells", "groups", "shared_or", "shared", "treatments", "links")] <- list(
    cells, groups, shared_or, shared_or, treatments, links
  )
  obj$accounting$n_imputation_slopes_by_group <- vapply(groups, `[[`, integer(1L), "p")
  obj$accounting$n_or_shared_slopes <- as.integer(p_or)
  payload <- cp04_payload(obj)
  actual <- unname(payload["entries"])
  if (actual > retained_bound + pair_entries) {
    augmented_abort("Retained augmented backend bound underestimated.", "bound")
  }
  cp04_budget(
    ctrl, c(obj$accounting$borrowed_entries, actual), fit_scratch,
    obj$accounting$work_ops_bound, actual, "retained augmented observation-correction backend"
  )
  obj$accounting$owned_entries <- actual
  obj$accounting$character_bytes <- unname(payload["char_bytes"])
  ipw_owner_seal(obj, "augmented_backend")
}

fit_chained_augmented <- fit_augmented_backend

validate_augmented_backend <- function(backend) {
  if (!inherits(backend, "cdid_augmented_backend") ||
      !ipw_owner_seal_ok(backend, "augmented_backend")) {
    augmented_abort("Expected the immutable sealed augmented observation-correction backend (`cdid_augmented_backend`).",
      "integrity")
  }
  reader <- if (is.environment(backend$reader_ref) &&
                exists("cp04_reader", backend$reader_ref, inherits = FALSE)) {
    get("cp04_reader", backend$reader_ref, inherits = FALSE)
  } else {
    NULL
  }
  validate_cp04_reader(reader)
  validate_selection_models(backend$selection_backend)
  if (!identical(backend$selection_backend$reader_ref, backend$reader_ref)) {
    augmented_abort("The retained observation-probability models do not belong to this reader.", "integrity")
  }
  d_owner <- if (is.environment(backend$design_ref) &&
                 exists("imputation_design_owner", backend$design_ref, inherits = FALSE)) {
    get("imputation_design_owner", backend$design_ref, inherits = FALSE)
  } else {
    NULL
  }
  if (is.null(backend$layout)) validate_imputation_design_owner(d_owner, reader = reader) else
    chained_backend_reference(backend)
  invisible(backend)
}

augmented_backend_X_F <- function(backend) {
  if (!is.null(backend$X_F)) return(backend$X_F)
  get("imputation_design_owner", backend$design_ref, inherits = FALSE)$X
}

# Rows/coordinates needed by one link batch. Every cell whose Z enters a
# requested link or a contributing OR block is reconstructed; every block
# contributes whenever a shared OR slope exists (Trap 2); imputation groups
# with slopes need all of their cells' scores; selection rows are needed on
# the theta coordinates of the reconstructed cells only (never rate rows).
augmented_if_needed_sets <- function(backend, links, details = FALSE) {
  cells <- backend$cells
  groups <- backend$groups
  sh <- backend$shared_or
  if (details) {
    zids <- seq_along(cells)
    cids <- zids
    gids <- seq_along(groups)
    bids <- seq_len(sh$n_blocks)
    tids <- seq_along(backend$treatments)
    uids <- sort(unique(c(integer(0L), unlist(backend$up_map$theta_coords))))
  } else {
    l_sub <- backend$links[links]
    cids_link <- sort(unique(c(
      vapply(l_sub, `[[`, integer(1L), "target_cell_id"),
      unlist(lapply(l_sub, `[[`, "control_cell_ids"))
    )))
    bids <- if (sh$p > 0L) {
      seq_len(sh$n_blocks)
    } else {
      sort(unique(vapply(l_sub, `[[`, integer(1L), "or_block_id")))
    }
    cids_or <- sort(unique(unlist(lapply(sh$blocks[bids], `[[`, "control_cell_ids"))))
    zids <- sort(unique(c(cids_link, cids_or)))
    gids <- sort(unique(vapply(cells[zids], `[[`, integer(1L), "group_id")))
    cids <- sort(unique(c(
      zids,
      unlist(lapply(groups[gids], function(gf) if (gf$p > 0L) gf$cell_ids else integer(0L)))
    )))
    tids <- sort(unique(vapply(l_sub, `[[`, integer(1L), "trt_id")))
    uids <- sort(unique(c(integer(0L), unlist(lapply(cells[zids], `[[`, "theta_ids")))))
  }
  list(
    zids = as.integer(zids),
    cids = as.integer(cids),
    gids = as.integer(gids),
    bids = as.integer(bids),
    tids = as.integer(tids),
    uids = as.integer(uids),
    q_ids = as.integer(unique(vapply(cells[zids], `[[`, integer(1L), "q_id")))
  )
}

# Affine-in-B output/scratch/work envelope of one row/link batch: the
# pure-imputation envelope terms, the B x |uids| selection rows, the per-cell
# q evaluation on observed rows and the accepted selection IF envelope
# (dr_corollary1_if_envelope template).
augmented_if_envelope <- function(backend, B, links = seq_along(backend$links), details = FALSE) {
  validate_augmented_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  # The retained layout matrix holds all recipe copies of the original rows,
  # even though the evaluator visits those copies in bounded sub-chunks.
  if (!is.null(backend$layout)) B <- B * max(1, backend$layout$by_orig$count)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  sets <- augmented_if_needed_sets(backend, links, details)
  nl <- length(links)
  nc <- length(sets$cids)
  nz <- length(sets$zids)
  ng <- length(sets$gids)
  nb <- length(sets$bids)
  nt <- length(sets$tids)
  nu <- length(sets$uids)
  G <- length(st$cohorts)
  p_or <- backend$outcome_recipe$p
  p_trt <- backend$treatment_recipe$p
  n_cells <- length(backend$cells)
  d_sel <- ncol(st$Z)

  group_width <- if (details) sum(vapply(backend$groups, `[[`, integer(1L), "p")) else 0
  n_alpha <- if (details) {
    backend$shared_or$n_blocks
  } else {
    length(unique(vapply(backend$links[links], `[[`, integer(1L), "or_block_id")))
  }
  width <- 3 * G + group_width + 24 * n_alpha + 40 * p_or + nu + nt * (6 * p_trt + 8) + 48 * nl +
    if (details) 4 * n_cells + 32 * nl + nu + 32 else 8
  metadata <- 16384 + 1024 * (nc + length(backend$groups) + nb + nt + nl + G + nu) +
    64 * n_cells * (p_or + 2) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1L), "link_id"), type = "bytes"))
  output <- B * (width + 16) + metadata
  scratch <- 64 * B * (width + length(backend$design_columns) + d_sel + p_or + p_trt + nz + 8) +
    128 * ((p_or + 2)^2 + (p_trt + 2)^2) + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- 128 * B * (
    (nc + 1) * (p_or + 2)^2 * log_a +
      (ng + nb + nl + 1) * (p_or + 1) * (p_or + nc + 1) +
      nt * ((p_trt + 2)^2 + p_trt * (G + 4)) +
      G + nl * (nc + 1) * (p_or + p_trt + G + 8) +
      (nz + nb + 2 * nl) * (d_sel + nu + log_a + 4) +
      (nb + p_or + 2 * nl) * (nu + 1)
  ) + 64 * metadata
  selection_pre <- selection_if_budget(backend$selection_backend, B, sets$q_ids, details = details)
  scratch <- scratch + selection_pre$output + selection_pre$scratch
  work <- work + selection_pre$work
  c(sets, list(output = output, scratch = scratch, work = work))
}

augmented_if_budget <- function(backend, B, links = seq_along(backend$links), details = FALSE) {
  pre <- augmented_if_envelope(backend, B, links, details)
  cp04_budget(
    backend$control,
    c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    pre$scratch,
    c(backend$accounting$work_ops_bound, pre$work),
    pre$output,
    "augmented observation-correction IF chunk"
  )
  pre
}

# Complete original-unit IF rows of one bounded row chunk (A6-A9). Channels:
# dir (augmented residuals through Z), imputation (M^Delta adjoint), selection
# (merged -R^q loadings on theta rows), OR (direct, rate, j, q forcing over all
# blocks), PS (unchanged roster treatment rows) and rate (single count).
augmented_evaluate_if_chunk_core <- function(st, X_F, pair_reader, backend, rows, links,
                                             scale = "frame", details = FALSE, sets) {
  B <- length(rows)
  N <- st$N
  G <- length(st$cohorts)
  recipe <- backend$outcome_recipe
  trt_recipe <- backend$treatment_recipe
  sel_fit <- backend$selection_backend
  up_map <- backend$up_map
  cells <- backend$cells
  groups <- backend$groups
  sh <- backend$shared_or
  n_cells <- length(cells)
  n_groups <- length(groups)

  ci <- st$cohort_index[rows]
  w <- st$w[rows]
  A <- as.double(st$a_index[rows] > 0L)
  X_chunk <- X_F[rows, , drop = FALSE]

  # 1. Repeat-rate and cohort-share rows (unchanged pure-imputation rows).
  xi_a <- matrix(0.0, B, G)
  xi_tau <- matrix(0.0, B, G)
  xi_eta <- matrix(0.0, B, G)
  coh_labels <- imputation_cohort_label(st$cohorts)
  colnames(xi_a) <- paste0("a.k", coh_labels)
  colnames(xi_tau) <- paste0("tau.k", coh_labels)
  colnames(xi_eta) <- paste0("eta.k", coh_labels)
  for (k in seq_len(G)) {
    I_k <- as.double(ci == k)
    if (st$mass_F[k] > 0) xi_a[, k] <- w * I_k * (A - st$a[k]) / st$mass_F[k]
    xi_tau[, k] <- w * (I_k - st$tau[k]) / st$mean_w
    xi_eta[, k] <- w * A * (I_k - st$eta[k]) / st$mean_wA
  }

  # 2. Accepted original-unit selection rows on the needed theta coordinates.
  uids <- sets$uids
  up <- if (length(uids)) {
    evaluate_selection_upstream_chunk(sel_fit, up_map, rows, uids)$up
  } else {
    matrix(0.0, B, 0L)
  }
  load_theta <- function(ids, coef) {
    if (!length(ids)) return(numeric(B))
    pos <- match(ids, uids)
    if (anyNA(pos)) augmented_abort("Selection loading outside the evaluated theta coordinates.", "integrity")
    as.vector(up[, pos, drop = FALSE] %*% coef)
  }

  # 3. Reconstructed responses of the needed cells on this chunk (once each).
  zcache <- vector("list", n_cells)
  cell_Z <- function(cid) {
    if (!is.null(zcache[[cid]])) return(zcache[[cid]])
    cf <- cells[[cid]]
    ix <- which(ci == cf$cohort_index)
    out <- if (length(ix)) {
      ev <- augmented_cell_rows(sel_fit, up_map, X_F, pair_reader, cf, groups[[cf$group_id]], rows[ix], FALSE)
      list(ix = ix, Z = ev$Z)
    } else {
      list(ix = integer(0L), Z = numeric(0L))
    }
    zcache[[cid]] <<- out
    out
  }

  # 4. Imputation scores and the sparse adjoint pullback (pure-imputation algebra).
  read_cell_score <- function(cid) {
    cf <- cells[[cid]]
    gf <- groups[[cf$group_id]]
    ix <- which(ci == cf$cohort_index)
    if (!length(ix)) return(NULL)
    hit <- pair_reader(rows[ix], cf$pair, TRUE)
    obs_local <- which(hit$S == 1)
    if (!length(obs_local)) return(NULL)
    rr <- ix[obs_local]
    dy_rr <- hit$dy[obs_local]
    Z_rr <- if (gf$p > 0L) X_chunk[rr, gf$kept_cols_F, drop = FALSE] else matrix(0.0, length(rr), 0L)
    d_rr <- cf$L * Z_rr
    jhat_rr <- if (gf$p > 0L) cf$alpha + as.vector(d_rr %*% gf$beta) else rep(cf$alpha, length(rr))
    list(rr = rr, d_rr = d_rr, s_rr = w[rr] * (dy_rr - jhat_rr))
  }

  xi_beta_imp <- if (details) vector("list", n_groups) else NULL
  if (details) for (gid in seq_len(n_groups)) {
    gf <- groups[[gid]]
    h_beta <- matrix(0.0, B, gf$p, dimnames = list(NULL, gf$kept_columns))
    if (gf$p > 0L) {
      for (cid in gf$cell_ids) {
        obs_rec <- read_cell_score(cid)
        if (is.null(obs_rec)) next
        cf <- cells[[cid]]
        d_tilde <- if (recipe$intercept) sweep(obs_rec$d_rr, 2L, cf$mean_d, FUN = "-") else obs_rec$d_rr
        h_beta[obs_rec$rr, ] <- h_beta[obs_rec$rr, , drop = FALSE] + obs_rec$s_rr * d_tilde
      }
      sol_b <- if (recipe$intercept) {
        t(or_pooled_wls_gram_solve(gf$fit, t(h_beta), normalizer = N))
      } else {
        t(or_wls_gram_solve(gf$fit, t(h_beta), normalizer = N))
      }
      colnames(sol_b) <- gf$kept_columns
      xi_beta_imp[[gid]] <- sol_b
    } else {
      xi_beta_imp[[gid]] <- h_beta
    }
  }
  alpha_if <- function(cid, obs_rec = NULL) {
    cf <- cells[[cid]]
    gf <- groups[[cf$group_id]]
    out <- numeric(B)
    if (!recipe$intercept) return(out)
    if (missing(obs_rec)) obs_rec <- read_cell_score(cid)
    if (!is.null(obs_rec)) out[obs_rec$rr] <- obs_rec$s_rr / cf$mass
    if (gf$p > 0L) out <- out - as.vector(xi_beta_imp[[cf$group_id]] %*% cf$mean_d)
    out
  }
  eval_adjoint_chunk <- function(adj) {
    out_vec <- numeric(B)
    for (cid in sets$cids) {
      obs_rec <- read_cell_score(cid)
      if (is.null(obs_rec)) next
      cf <- cells[[cid]]
      gf <- groups[[cf$group_id]]
      t_lam <- if (recipe$intercept) adj$tilde_lambda[cid] else 0.0
      t_nu <- adj$tilde_nu[[cf$group_id]]
      val <- rep(t_lam, length(obs_rec$rr))
      if (gf$p > 0L && length(t_nu) > 0L) {
        d_tilde <- if (recipe$intercept) sweep(obs_rec$d_rr, 2L, cf$mean_d, FUN = "-") else obs_rec$d_rr
        val <- val + as.vector(d_tilde %*% t_nu)
      }
      out_vec[obs_rec$rr] <- out_vec[obs_rec$rr] + obs_rec$s_rr * val
    }
    out_vec
  }
  eval_imp_linear <- function(lambda_vec, nu_list) {
    tilde_lambda <- numeric(n_cells)
    tilde_nu <- vector("list", n_groups)
    for (gid in seq_len(n_groups)) {
      gf <- groups[[gid]]
      rhs <- numeric(gf$p)
      for (cid in gf$cell_ids) {
        cf <- cells[[cid]]
        lam <- if (recipe$intercept) lambda_vec[cid] else 0
        if (recipe$intercept) tilde_lambda[cid] <- lam / cf$mass
        if (gf$p) rhs <- rhs + nu_list[[cid]] - lam * cf$mean_d
      }
      tilde_nu[[gid]] <- if (!gf$p) numeric(0L) else as.double(
        if (recipe$intercept) or_pooled_wls_gram_solve(gf$fit, rhs, normalizer = N)
        else or_wls_gram_solve(gf$fit, rhs, normalizer = N))
    }
    eval_adjoint_chunk(list(tilde_lambda = tilde_lambda, tilde_nu = tilde_nu))
  }
  zero_imp_loadings <- function() lapply(cells, function(cf) numeric(groups[[cf$group_id]]$p))

  # 5. Shared counterfactual OR rows: direct, rate, j and q channels.
  or_if <- NULL
  alpha_bids <- integer(0L)
  if (length(sets$bids) > 0L) {
    p_or <- sh$p
    n_b <- sh$n_blocks
    alpha_bids <- if (details) seq_len(n_b) else
      sort(unique(vapply(backend$links[links], `[[`, integer(1L), "or_block_id")))
    n_alpha <- length(alpha_bids)
    rhs_beta_dir <- rhs_beta_a <- rhs_beta_j <- rhs_beta_q <- matrix(0.0, B, p_or)
    rhs_alpha_dir <- rhs_alpha_a <- rhs_alpha_j <- rhs_alpha_q <- matrix(0.0, B, n_alpha)

    for (bid in sets$bids) {
      apos <- match(bid, alpha_bids)
      bf <- sh$blocks[[bid]]
      for (idx_k in seq_along(bf$control_cell_ids)) {
        cid <- bf$control_cell_ids[idx_k]
        cf <- cells[[cid]]
        gf <- groups[[cf$group_id]]
        ci_k <- cf$cohort_index
        a_k <- st$a[ci_k]
        zc <- cell_Z(cid)
        ix <- zc$ix
        if (length(ix)) {
          w_syn_ix <- w[ix] * a_k
          d_or_ix <- if (p_or > 0L) bf$L * X_chunk[ix, sh$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix), 0L)
          mhat_ix <- if (p_or > 0L) {
            (if (recipe$intercept) bf$alpha else 0.0) + as.vector(d_or_ix %*% sh$beta)
          } else {
            rep(if (recipe$intercept) bf$alpha else 0.0, length(ix))
          }
          r_ix <- zc$Z - mhat_ix
          if (recipe$intercept && !is.na(apos)) {
            rhs_alpha_dir[ix, apos] <- rhs_alpha_dir[ix, apos] + w_syn_ix * r_ix
          }
          if (p_or > 0L) {
            d_tilde_ix <- if (recipe$intercept) sweep(d_or_ix, 2L, bf$mean_d, FUN = "-") else d_or_ix
            rhs_beta_dir[ix, ] <- rhs_beta_dir[ix, , drop = FALSE] + (w_syn_ix * r_ix) * d_tilde_ix
          }
        }
        xi_ak <- xi_a[, ci_k]
        if (recipe$intercept && !is.na(apos)) {
          rhs_alpha_a[, apos] <- rhs_alpha_a[, apos] + bf$a_force_alpha[idx_k] * xi_ak
        }
        if (p_or > 0L) rhs_beta_a <- rhs_beta_a + outer(xi_ak, bf$a_force_beta[, idx_k])
        if (details) {
          xi_ac <- alpha_if(cid)
          xi_br <- xi_beta_imp[[cf$group_id]]
          if (recipe$intercept && !is.na(apos)) {
            rhs_alpha_j[, apos] <- rhs_alpha_j[, apos] + bf$j_force_alpha_alpha[idx_k] * xi_ac +
              (if (gf$p > 0L) as.vector(xi_br %*% bf$j_force_alpha_beta[[idx_k]]) else 0.0)
          }
          if (p_or > 0L) {
            if (recipe$intercept) rhs_beta_j <- rhs_beta_j + outer(xi_ac, bf$j_force_beta_alpha[, idx_k])
            if (gf$p > 0L) rhs_beta_j <- rhs_beta_j + xi_br %*% t(bf$j_force_beta_beta[[idx_k]])
          }
        }
      }
      if (recipe$intercept && !is.na(apos)) rhs_alpha_q[, apos] <- load_theta(bf$q_alpha_ids, bf$q_alpha)
    }
    if (p_or > 0L && length(sh$q_beta_ids)) {
      pos <- match(sh$q_beta_ids, uids)
      if (anyNA(pos)) augmented_abort("Shared OR selection forcing outside the evaluated theta coordinates.", "integrity")
      rhs_beta_q <- up[, pos, drop = FALSE] %*% t(sh$q_beta)
    }

    if (!details) {
      # Global shared-slope j forcing collects every contributing block before
      # each scalar pullback (row width p_OR, independent of the cell count).
      if (p_or) for (qq in seq_len(p_or)) {
        lam <- numeric(n_cells)
        nu <- zero_imp_loadings()
        for (bid in sets$bids) {
          bf <- sh$blocks[[bid]]
          for (idx_k in seq_along(bf$control_cell_ids)) {
            cid <- bf$control_cell_ids[idx_k]
            if (recipe$intercept) lam[cid] <- lam[cid] + bf$j_force_beta_alpha[qq, idx_k]
            nu[[cid]] <- nu[[cid]] + bf$j_force_beta_beta[[idx_k]][qq, ]
          }
        }
        rhs_beta_j[, qq] <- eval_imp_linear(lam, nu)
      }
      if (recipe$intercept) for (apos in seq_along(alpha_bids)) {
        bf <- sh$blocks[[alpha_bids[apos]]]
        lam <- numeric(n_cells)
        nu <- zero_imp_loadings()
        for (idx_k in seq_along(bf$control_cell_ids)) {
          cid <- bf$control_cell_ids[idx_k]
          lam[cid] <- lam[cid] + bf$j_force_alpha_alpha[idx_k]
          nu[[cid]] <- nu[[cid]] + bf$j_force_alpha_beta[[idx_k]]
        }
        rhs_alpha_j[, apos] <- eval_imp_linear(lam, nu)
      }
    }

    solve_or_channel <- function(rhs_a_mat, rhs_b_mat) {
      xi_b <- if (p_or > 0L) {
        if (recipe$intercept) {
          t(or_pooled_wls_gram_solve(sh$fit, t(rhs_b_mat), normalizer = N))
        } else {
          t(or_wls_gram_solve(sh$fit, t(rhs_b_mat), normalizer = N))
        }
      } else {
        matrix(0.0, B, 0L)
      }
      colnames(xi_b) <- sh$kept_columns
      xi_al <- matrix(0.0, B, n_alpha)
      if (recipe$intercept) {
        for (apos in seq_along(alpha_bids)) {
          bf <- sh$blocks[[alpha_bids[apos]]]
          xi_al[, apos] <- (rhs_a_mat[, apos] / bf$mass) -
            (if (p_or > 0L) as.vector(xi_b %*% bf$mean_d) else 0.0)
        }
      }
      list(alpha = xi_al, beta = xi_b)
    }
    ch_dir <- solve_or_channel(rhs_alpha_dir, rhs_beta_dir)
    ch_a <- solve_or_channel(rhs_alpha_a, rhs_beta_a)
    ch_j <- solve_or_channel(rhs_alpha_j, rhs_beta_j)
    ch_q <- solve_or_channel(rhs_alpha_q, rhs_beta_q)
    alpha_tot <- ch_dir$alpha + ch_a$alpha + ch_j$alpha + ch_q$alpha
    beta_tot <- ch_dir$beta + ch_a$beta + ch_j$beta + ch_q$beta
    theta_tot <- if (recipe$intercept) cbind(alpha_tot, beta_tot) else beta_tot
    if (ncol(theta_tot) > 0L) colnames(theta_tot) <- paste0("OR.theta", seq.int(0L, ncol(theta_tot) - 1L))
    or_if <- list(
      direct = ch_dir, a = ch_a, j = ch_j, q = ch_q,
      alpha = alpha_tot, beta = beta_tot, theta = theta_tot,
      alpha_block_ids = alpha_bids,
      rhs_beta_score = rhs_beta_dir + rhs_beta_a + rhs_beta_j + rhs_beta_q,
      rhs_alpha_score = rhs_alpha_dir + rhs_alpha_a + rhs_alpha_j + rhs_alpha_q
    )
  }
  or_term <- function(ch, bid, d_bar) {
    p_or <- sh$p
    (if (recipe$intercept) ch$alpha[, match(bid, alpha_bids)] else 0.0) +
      (if (p_or > 0L) as.vector(ch$beta %*% d_bar) else 0.0)
  }

  # 6. Roster treatment logit rows (unchanged pure-imputation rows; no j/q path).
  trt_if <- NULL
  if (length(sets$tids) > 0L) {
    trt_if <- vector("list", length(backend$treatments))
    p_trt <- trt_recipe$p
    trt_cols_F <- trt_recipe$index
    for (tid in sets$tids) {
      tf <- backend$treatments[[tid]]
      score_dir <- matrix(0.0, B, p_trt)
      ix <- which(ci %in% tf$all_cohort_indices)
      if (length(ix)) {
        ci_ix <- ci[ix]
        y_ix <- as.double(ci_ix == tf$target_cohort_index)
        w_syn_ix <- w[ix] * st$a[ci_ix]
        X_trt_ix <- X_chunk[ix, trt_cols_F, drop = FALSE]
        lv_ix <- ipw_logit_values(X_trt_ix, tf$coefficients)
        score_dir[ix, ] <- (w_syn_ix * (y_ix - lv_ix$prob)) * X_trt_ix
      }
      score_a <- xi_a[, tf$all_cohort_indices, drop = FALSE] %*% t(tf$a_forcing)
      xi_g_dir <- t(cp04_solve(tf$H_chol, t(score_dir)))
      xi_g_a <- t(cp04_solve(tf$H_chol, t(score_a)))
      xi_g_tot <- xi_g_dir + xi_g_a
      colnames(xi_g_dir) <- colnames(xi_g_a) <- colnames(xi_g_tot) <- paste0("gamma", seq.int(0L, p_trt - 1L))
      trt_if[[tid]] <- list(direct = xi_g_dir, a = xi_g_a, gamma = xi_g_tot,
        score_direct = score_dir, score_a = score_a)
    }
  }

  # 7. Link rows: the six channels of (A9).
  nl <- length(links)
  link_names <- vapply(backend$links[links], `[[`, character(1L), "link_id")
  phi <- matrix(0.0, B, nl, dimnames = list(NULL, link_names))
  link_details <- if (details) vector("list", nl) else NULL
  for (j in seq_along(links)) {
    lf <- backend$links[[links[j]]]
    cf_g <- cells[[lf$target_cell_id]]
    ci_g <- lf$target_cohort_index
    a_g <- st$a[ci_g]
    bid <- lf$or_block_id
    bf <- sh$blocks[[bid]]
    p_or <- sh$p

    # Treated ratios T and B_g.
    zg <- cell_Z(lf$target_cell_id)
    ix_g <- zg$ix
    w_plus <- numeric(B)
    Z_g_vec <- numeric(B)
    m_g_vec <- numeric(B)
    if (length(ix_g)) {
      w_plus[ix_g] <- w[ix_g] * a_g
      Z_g_vec[ix_g] <- zg$Z
      Z_or_g <- if (p_or > 0L) X_chunk[ix_g, sh$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix_g), 0L)
      m_g_vec[ix_g] <- if (p_or > 0L) {
        (if (recipe$intercept) bf$alpha else 0.0) + as.vector((lf$L * Z_or_g) %*% sh$beta)
      } else {
        rep(if (recipe$intercept) bf$alpha else 0.0, length(ix_g))
      }
    }
    psi_denplus <- w_plus + st$mass_F[ci_g] * xi_a[, ci_g]
    direct_T <- w_plus * (Z_g_vec - lf$T_gl) / lf$d_plus
    a_T <- lf$a_coeff_T_g * xi_a[, ci_g]
    j_T_adj <- eval_adjoint_chunk(lf$adjoint_T)
    j_T <- if (details) {
      lf$lam_T * alpha_if(lf$target_cell_id) +
        (if (length(lf$nu_T)) as.vector(xi_beta_imp[[cf_g$group_id]] %*% lf$nu_T) else 0.0)
    } else {
      j_T_adj
    }
    sel_T <- load_theta(lf$sel_T_ids, lf$sel_T)
    psi_T <- direct_T + a_T + j_T + sel_T

    direct_Bg <- w_plus * (m_g_vec - lf$B_gl) / lf$d_plus
    a_Bg <- lf$a_coeff_Bg_g * xi_a[, ci_g]
    or_Bg_dir <- or_term(or_if$direct, bid, lf$d_bar_plus_or)
    or_Bg_a <- or_term(or_if$a, bid, lf$d_bar_plus_or)
    or_Bg_j <- or_term(or_if$j, bid, lf$d_bar_plus_or)
    or_Bg_q <- or_term(or_if$q, bid, lf$d_bar_plus_or)
    psi_Bg <- direct_Bg + a_Bg + or_Bg_dir + or_Bg_a + or_Bg_j + or_Bg_q

    # Control ratios C and B_o.
    tf <- backend$treatments[[lf$trt_id]]
    tif <- trt_if[[lf$trt_id]]
    w_minus <- numeric(B)
    Z_k_vec <- numeric(B)
    m_k_vec <- numeric(B)
    a_denminus <- numeric(B)
    a_C <- numeric(B)
    a_Bo <- numeric(B)
    j_C_fwd <- numeric(B)
    for (idx_k in seq_along(lf$control_cohort_indices)) {
      ci_k <- lf$control_cohort_indices[idx_k]
      cid_k <- lf$control_cell_ids[idx_k]
      cf_k <- cells[[cid_k]]
      a_k <- st$a[ci_k]
      zk <- cell_Z(cid_k)
      ix_k <- zk$ix
      if (length(ix_k)) {
        lv_k <- ipw_logit_values(X_chunk[ix_k, trt_recipe$index, drop = FALSE], tf$coefficients)
        w_minus[ix_k] <- (w[ix_k] * a_k) * lv_k$odds
        Z_k_vec[ix_k] <- zk$Z
        Z_or_k <- if (p_or > 0L) X_chunk[ix_k, sh$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix_k), 0L)
        m_k_vec[ix_k] <- if (p_or > 0L) {
          (if (recipe$intercept) bf$alpha else 0.0) + as.vector((lf$L * Z_or_k) %*% sh$beta)
        } else {
          rep(if (recipe$intercept) bf$alpha else 0.0, length(ix_k))
        }
      }
      xi_ak <- xi_a[, ci_k]
      a_C <- a_C + lf$a_coeff_C_k[idx_k] * xi_ak
      a_Bo <- a_Bo + lf$a_coeff_Bo_k[idx_k] * xi_ak
      a_denminus <- a_denminus + ((lf$d_minus * lf$m_minus_k[idx_k]) / a_k) * xi_ak
      if (details) {
        gf_k <- groups[[cf_k$group_id]]
        j_C_fwd <- j_C_fwd + lf$lam_C_k[idx_k] * alpha_if(cid_k) +
          (if (gf_k$p > 0L) as.vector(xi_beta_imp[[cf_k$group_id]] %*% lf$nu_C_k[[idx_k]]) else 0.0)
      }
    }
    direct_C <- w_minus * (Z_k_vec - lf$C_ol) / lf$d_minus
    j_C_adj <- eval_adjoint_chunk(lf$adjoint_C)
    j_C <- if (details) j_C_fwd else j_C_adj
    sel_C <- load_theta(lf$sel_C_ids, lf$sel_C)
    trt_C_dir <- as.vector(tif$direct %*% lf$g_C_gamma)
    trt_C_a <- as.vector(tif$a %*% lf$g_C_gamma)
    psi_C <- direct_C + a_C + j_C + sel_C + trt_C_dir + trt_C_a

    direct_Bo <- w_minus * (m_k_vec - lf$B_ol) / lf$d_minus
    trt_Bo_dir <- as.vector(tif$direct %*% lf$g_Bo_gamma)
    trt_Bo_a <- as.vector(tif$a %*% lf$g_Bo_gamma)
    or_Bo_dir <- or_term(or_if$direct, bid, lf$d_bar_minus_or)
    or_Bo_a <- or_term(or_if$a, bid, lf$d_bar_minus_or)
    or_Bo_j <- or_term(or_if$j, bid, lf$d_bar_minus_or)
    or_Bo_q <- or_term(or_if$q, bid, lf$d_bar_minus_or)
    psi_Bo <- direct_Bo + a_Bo + trt_Bo_dir + trt_Bo_a + or_Bo_dir + or_Bo_a + or_Bo_j + or_Bo_q

    comp <- cbind(
      direct = (direct_T - direct_Bg) - (direct_C - direct_Bo),
      rate = (a_T - a_Bg) - (a_C - a_Bo),
      imputation = j_T - j_C,
      selection = sel_T - sel_C,
      treatment_direct = -(trt_C_dir - trt_Bo_dir),
      treatment_rate = -(trt_C_a - trt_Bo_a),
      OR_direct = -or_Bg_dir + or_Bo_dir,
      OR_rate = -or_Bg_a + or_Bo_a,
      OR_imputation = -or_Bg_j + or_Bo_j,
      OR_selection = -or_Bg_q + or_Bo_q
    )
    phi[, j] <- rowSums(comp)

    if (details) {
      psi_denminus <- w_minus + a_denminus + as.vector(tif$gamma %*% lf$denminus_gamma)
      link_details[[j]] <- list(
        link_id = lf$link_id,
        denplus = psi_denplus,
        numT = lf$d_plus * psi_T + lf$T_gl * psi_denplus,
        T = psi_T,
        numBg = lf$d_plus * psi_Bg + lf$B_gl * psi_denplus,
        Bg = psi_Bg,
        denminus = psi_denminus,
        numC = lf$d_minus * psi_C + lf$C_ol * psi_denminus,
        C = psi_C,
        numBo = lf$d_minus * psi_Bo + lf$B_ol * psi_denminus,
        Bo = psi_Bo,
        dr = phi[, j],
        components = comp,
        reference_components = cbind(
          j_score_pullback = comp[, "imputation"] + comp[, "OR_imputation"],
          q_score_pullback = comp[, "selection"] + comp[, "OR_selection"],
          all_rate_paths = comp[, "rate"] + comp[, "OR_rate"] + comp[, "treatment_rate"],
          remaining_direct_and_downstream_scores = comp[, "direct"] + comp[, "OR_direct"] +
            comp[, "treatment_direct"]
        ),
        j_T_forward = j_T,
        j_T_adjoint = j_T_adj,
        j_C_forward = j_C_fwd,
        j_C_adjoint = j_C_adj
      )
    }
  }

  cp04_finite(phi, "complete augmented link IF")
  phi_out <- if (identical(scale, "core")) st$lambda * phi else phi

  det_out <- NULL
  if (details) {
    xi_alpha_imp <- matrix(0.0, B, n_cells,
      dimnames = list(NULL, vapply(cells, `[[`, character(1L), "name")))
    cell_scores <- matrix(0.0, B, n_cells)
    for (cid in seq_len(n_cells)) {
      obs_rec <- read_cell_score(cid)
      xi_alpha_imp[, cid] <- alpha_if(cid, obs_rec)
      if (!is.null(obs_rec)) cell_scores[obs_rec$rr, cid] <- obs_rec$s_rr
    }
    det_out <- list(
      a = xi_a,
      tau = xi_tau,
      eta = xi_eta,
      lambda = A - st$lambda,
      imputation_alpha = xi_alpha_imp,
      imputation_beta = xi_beta_imp,
      imputation_cell_scores = cell_scores,
      selection_theta = up,
      selection_theta_ids = uids,
      or_if = or_if,
      trt_if = trt_if,
      links = link_details
    )
  }

  list(
    rows = rows,
    id = st$id[rows],
    links = links,
    phi = phi_out,
    phi_frame = phi,
    scale = scale,
    normalization_count = if (identical(scale, "core")) st$n_A else st$N,
    details = det_out
  )
}

evaluate_augmented_if_chunk <- function(backend, rows, links = seq_along(backend$links),
                                        scale = "frame", details = FALSE) {
  validate_augmented_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  valid_index <- function(x, n) {
    is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
      length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) &&
      all(x >= 1 & x <= n) && !anyDuplicated(x)
  }
  if (length(rows) > backend$control$max_chunk_rows) {
    cdid_ipw_owner_abort("Replay rows exceed the explicit chunk bound.", "cdid_ipw_owner_bound_error")
  }
  if (!valid_index(rows, st$N) ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) ||
      !scale %in% c("frame", "core")) {
    augmented_abort("Invalid explicit IF row/link batch, scale or details flag.")
  }

  if (!is.null(backend$layout)) return(evaluate_augmented_layout_if_chunk(backend, reader, rows, links, scale, details))
  pre <- augmented_if_budget(backend, length(rows), links, details)
  X_F <- augmented_backend_X_F(backend)
  pair_reader <- function(r_sub, pair, read_outcomes = TRUE) {
    pr <- cp04_pair_rows(reader, r_sub, pair, read_outcomes)
    if (read_outcomes && !is.null(backend$pair_outcomes)) {
      pk <- paste0(as.character(pair$s), ":", as.character(pair$t))
      dy <- backend$pair_outcomes[[pk]][r_sub]
      observed <- pr$S == 1
      if (!identical(!is.na(dy), observed) || any(!is.finite(dy[observed]))) {
        augmented_abort("Retained private pair differences disagree with the sealed panel.", "support")
      }
      pr$dy <- numeric(length(r_sub))
      pr$dy[observed] <- dy[observed]
    }
    pr
  }

  res <- augmented_evaluate_if_chunk_core(
    st = st,
    X_F = X_F,
    pair_reader = pair_reader,
    backend = backend,
    rows = as.integer(rows),
    links = as.integer(links),
    scale = scale,
    details = details,
    sets = pre
  )

  out <- list(
    rows = res$rows,
    id = res$id,
    links = res$links,
    phi = res$phi,
    scale = res$scale,
    normalization_count = res$normalization_count,
    diagnostics = list(
      n_imputation_cell_evaluations = length(pre$cids),
      n_imputation_group_solves = length(pre$gids),
      n_reconstructed_cells = length(pre$zids),
      n_or_block_evaluations = length(pre$bids),
      n_treatment_score_evaluations = length(pre$tids),
      n_selection_theta_coordinates = length(pre$uids),
      n_observation_logit_evaluations = length(pre$q_ids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  if (details) out$details <- res$details
  if (unname(cp04_payload(out)["entries"]) > pre$output) {
    augmented_abort("Augmented IF output bound underestimated.", "bound")
  }
  out
}

# D60 recipe-cohort layout: layout rows of these original rows are evaluated by
# the unchanged core in bounded sub-chunks and summed into one row per unit.
evaluate_augmented_layout_if_chunk <- function(backend, reader, rows, links, scale, details) {
  if (details) augmented_abort("Diagnostic channel details are not retained on the recipe-cohort layout.", "unsupported")
  lay <- backend$layout; st <- reader$state
  pre <- augmented_if_budget(backend, length(rows), links, FALSE)
  map <- chained_reference_layout_rows(lay, as.integer(rows))
  BL <- backend$control$max_chunk_rows
  pair_reader <- function(r_sub, pair, read_outcomes = TRUE) cp04_pair_rows(reader, r_sub, pair, read_outcomes)
  phi_l <- matrix(0, length(map$rows), length(links))
  if (length(map$rows)) for (s0 in seq.int(1L, length(map$rows), by = BL)) {
    ii <- seq.int(s0, min(length(map$rows), s0 + BL - 1L))
    res <- augmented_evaluate_if_chunk_core(st = st, X_F = backend$X_F, pair_reader = pair_reader, backend = backend,
      rows = as.integer(map$rows[ii]), links = as.integer(links), scale = "frame", details = FALSE, sets = pre)
    phi_l[ii, ] <- res$phi
  }
  colnames(phi_l) <- vapply(backend$links[links], `[[`, character(1L), "link_id")
  B <- length(rows)
  phi <- chained_reference_layout_sum(phi_l, map$owner, B)
  cp04_finite(phi, "complete augmented link IF (recipe-cohort layout)")
  if (identical(scale, "core")) phi <- st$lambda * phi
  out <- list(rows = rows, id = lay$frame_id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (identical(scale, "core")) st$n_A else st$N,
    diagnostics = list(n_imputation_cell_evaluations = length(pre$cids), n_imputation_group_solves = length(pre$gids),
      n_reconstructed_cells = length(pre$zids), n_or_block_evaluations = length(pre$bids),
      n_treatment_score_evaluations = length(pre$tids), n_selection_theta_coordinates = length(pre$uids),
      n_observation_logit_evaluations = length(pre$q_ids), n_layout_rows = length(map$rows),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (unname(cp04_payload(out)["entries"]) > pre$output) augmented_abort("Augmented IF output bound underestimated.", "bound")
  out
}

replay_augmented_if <- function(backend, chunk_rows, link_batch_size, callback,
                                scale = "frame", details = FALSE) {
  validate_augmented_backend(backend)
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  LB <- ipw_owner_int_bound(link_batch_size, "link_batch_size")
  if (B > backend$control$max_chunk_rows || !is.function(callback) ||
      !identical(details, isTRUE(details)) || !is.character(scale) ||
      length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core")) {
    augmented_abort("Invalid explicit replay controls/callback.")
  }
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N
  L <- length(backend$links)
  total <- backend$accounting$work_ops_bound
  # Preflight the cumulative work of every batch before the first callback.
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    pre <- augmented_if_budget(backend, min(B, N), ll, details)
    total <- ipw_owner_safe_sum(
      c(total, ceiling(N / B) * pre$work),
      backend$control$max_work_ops,
      "total augmented observation-correction replay work"
    )
  }
  batches <- 0L
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1L))
      batch <- evaluate_augmented_if_chunk(backend, rr, ll, scale, details)
      callback(batch)
      rm(batch)
      batches <- batches + 1L
    }
  }
  invisible(list(
    n_frame_rows = N,
    n_links = L,
    batches = batches,
    scale = scale,
    normalization_count = if (identical(scale, "core")) reader$state$n_A else N,
    work_ops_bound = total,
    all_F_identities = TRUE
  ))
}

stream_augmented_covariance <- function(backend, links = seq_along(backend$links), retained = 0) {
  validate_augmented_backend(backend)
  cp04_stream_covariance(backend, links, retained = retained)
}
