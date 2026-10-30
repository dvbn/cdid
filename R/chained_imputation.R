# ==============================================================================
# Private Chained Pure-Imputation Backend (D52; IMPUTATION_SPEC.md)
# ==============================================================================
# Scope and architecture:
#   1. Supports `selection_est_method = "imputation"` across all three treatment
#      comparison methods (`est_method %in% c("ipw", "reg", "dr")`) and all three
#      imputation slope-pooling modes (`imputation_pooling %in% c("separate",
#      "within_cohort", "across_cohorts")`, default `"separate"`).
#   2. Fits NO observation-probability model (`pi` or `q`); `selection_model`,
#      `selection_pooling`, and `selection_xformla` are inactive and never forced.
#   3. Fitting dependencies:
#      - `est_method = "ipw"`: fits repeat rates `a_k`, actual-change imputation
#        `j_hat_{k,e}` on observed pairs with fixed unit weights `w_i`, and
#        synthetic-weight treatment propensity logits `gamma_hat_l` on F; fits
#        NO counterfactual OR model `theta_hat_OR`.
#      - `est_method = "reg"`: fits repeat rates `a_k`, actual-change imputation
#        `j_hat_{k,e}`, and shared counterfactual control OR `theta_hat_OR` on F;
#        fits NO treatment propensity logits `gamma_hat_l`.
#      - `est_method = "dr"`: fits `a_k`, `j_hat_{k,e}`, `gamma_hat_l`, and
#        `theta_hat_OR`.
#   4. Uses index-based stored-calendar transition count `L_e = i_t - i_s` via
#      `or_pooled_transitions`, partialled QR block WLS (`or_pooled_wls_*` /
#      `or_wls_*`) without dense cell-dummy or C^2 Gram matrices, group-local
#      retained slope dimensions `p_{imp,r}` with full-roster F alias verification,
#      sparse adjoint IF pullback, and complete original-F IF rows normalized by
#      `N = |F|` (and `(n_A / N)` on `"core"` scale with divisor `n_A`).
# ==============================================================================

imputation_abort <- function(message, class = "input", extra = list()) {
  stop(structure(
    c(list(message = message, call = NULL), extra),
    class = unique(c(
      paste0("cdid_imputation_", class, "_error"),
      "cdid_imputation_error",
      paste0("cdid_dr_", class, "_error"),
      paste0("cdid_dr_baseline_", class, "_error"),
      paste0("cdid_or_", class, "_error"),
      paste0("cdid_or_baseline_", class, "_error"),
      paste0("cdid_cp04_", class, "_error"),
      "cdid_ipw_fit_error", "cdid_error", "error", "condition"
    ))
  ))
}

imputation_pooling_choice <- function(imputation_pooling = "separate") {
  if (!is.character(imputation_pooling) || is.object(imputation_pooling) ||
      !is.null(dim(imputation_pooling)) || length(imputation_pooling) != 1L ||
      is.na(imputation_pooling) ||
      !imputation_pooling %in% c("separate", "within_cohort", "across_cohorts")) {
    imputation_abort(
      '`imputation_pooling` must be one of "separate", "within_cohort", or "across_cohorts".',
      "pooling"
    )
  }
  imputation_pooling
}

imputation_method_choice <- function(est_method = "dr", dr_method = NULL) {
  m <- if (!is.null(dr_method) && missing(est_method)) dr_method else est_method
  if (!is.character(m) || is.object(m) || !is.null(dim(m)) ||
      length(m) != 1L || is.na(m) || !m %in% c("dr", "ipw", "reg")) {
    imputation_abort(
      '`est_method` must be one of "dr", "ipw", or "reg".',
      "method"
    )
  }
  m
}

imputation_inactive_selection_record <- function() {
  list(
    selection_est_method = "imputation",
    active = FALSE,
    observation_model_fitted = FALSE,
    selection_formula_forced = FALSE,
    selection_model = NA_character_,
    selection_pooling = NA_character_,
    selection_xformla = NULL,
    universe = "F",
    note = "Pure imputation fits no observation-probability model; selection_model, selection_pooling, and selection_xformla are inactive and unforced."
  )
}

imputation_cohort_sort_key <- function(k) {
  k <- as.double(k)
  ifelse(is.infinite(k) | k == 0.0, -Inf, k)
}

imputation_cohort_label <- function(k) {
  vapply(as.double(k), function(x) {
    if (is.infinite(x) || x == 0.0) "0" else format(x, scientific = FALSE, trim = TRUE)
  }, character(1L))
}

# ------------------------------------------------------------------------------
# 1. Explicit Sealed Full-F Outcome/Treatment Design Owner & Resolution
# ------------------------------------------------------------------------------
# `frame_owner$Z_a` stores covariates on A only (`n_a_units` rows), while
# `reader$state$Z` may be either the full-F outcome/treatment design or a
# selection-probability design. `construct_imputation_design_owner` seals an
# explicit verified full-F outcome/treatment design matrix `X` (`N x d`) aligned
# to `reader$state$id` and verified against `frame_owner$Z_a` on all A rows.

construct_imputation_design_owner <- function(frame_owner, reader, id, X,
                                              control = frame_owner$control) {
  validate_ipw_frame_owner(frame_owner)
  validate_cp04_reader(reader)
  ctrl <- validate_ipw_owner_control(control)
  if (!identical(reader$frame_ref, frame_owner$ref) ||
      !identical(reader$panel_ref, frame_owner$panel_ref)) {
    imputation_abort("`reader` and `frame_owner` must share the same sealed owner references.", "identity")
  }
  st <- reader$state
  if (!identical(st$universe, "F") || !all(st$U) || any(st$rho != 1)) {
    imputation_abort(
      'Pure imputation (`selection_est_method = "imputation"`) requires the full representative roster `universe = "F"` with `U = TRUE` and `rho = 1` on all F rows.',
      "universe"
    )
  }
  N <- st$N
  if (!is.atomic(id) || is.object(id) || !is.null(dim(id)) || length(id) != N ||
      !identical(id, st$id) ||
      !is.matrix(X) || is.object(X) || !is.numeric(X) || nrow(X) != N || ncol(X) < 1L) {
    imputation_abort("Supply full-F `id` matching `reader$state$id` and a numeric `N x d` matrix `X`.", "design")
  }
  d <- ncol(X)
  cnames <- colnames(X)
  if (!is.character(cnames) || length(cnames) != d || anyNA(cnames) ||
      any(!nzchar(cnames)) || anyDuplicated(cnames)) {
    imputation_abort("Full-F outcome/treatment design `X` requires distinct nonempty column names.", "design")
  }
  ipw_owner_refuse_decorative_attrs(id, "full-F design id")
  ipw_owner_refuse_decorative_attrs(X, "full-F design X")
  ipw_owner_check_matrix_dims(N, d, ctrl$max_matrix_entries, "full-F imputation design")

  borrowed <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries
  bound <- N * (d + 2) + 4096 + unname(cp04_payload(list(dimnames(X), attributes(id)))["entries"])
  scratch <- 8 * min(N, ctrl$max_chunk_rows) * (d + 1) + 256
  cp04_budget(ctrl, c(borrowed, bound), scratch, 64 * N * (d + 1), bound, "full-F imputation design owner")

  for (start in seq.int(1L, N, by = min(N, ctrl$max_chunk_rows))) {
    rr <- seq.int(start, min(N, start + ctrl$max_chunk_rows - 1L))
    cp04_finite(X[rr, , drop = FALSE], "full-F imputation design")
  }

  # Verify exact agreement with frame_owner$Z_a on all A rows for shared columns
  a_pos <- which(st$a_index > 0L)
  a_idx <- st$a_index[a_pos]
  z_cols <- colnames(frame_owner$Z_a)
  if (!is.character(z_cols) || !all(z_cols %in% cnames)) {
    imputation_abort("Full-F design `X` must include all columns of `frame_owner$Z_a`.", "design")
  }
  x_match <- match(z_cols, cnames)
  if (length(a_pos) > 0L &&
      !identical(unname(storage.mode(X)), "double")) {
    storage.mode(X) <- "double"
  }
  if (length(a_pos) > 0L) for (start in seq.int(1L, length(a_pos), by = ctrl$max_chunk_rows)) {
    ix <- seq.int(start, min(length(a_pos), start + ctrl$max_chunk_rows - 1L))
    if (max(abs(X[a_pos[ix], x_match, drop = FALSE] - frame_owner$Z_a[a_idx[ix], , drop = FALSE])) > 0) {
      imputation_abort(
        "Full-F design `X` disagrees with `frame_owner$Z_a` on repeat-observed A units.",
        "identity"
      )
    }
  }

  obj <- structure(
    list(
      id = id,
      X = X,
      columns = cnames,
      N = N,
      d = d,
      reader_ref = reader$ref,
      frame_ref = frame_owner$ref,
      panel_ref = frame_owner$panel_ref,
      control = ctrl,
      accounting = list(
        owned_entries = 0,
        borrowed_entries = borrowed,
        scratch_bound = scratch,
        work_ops_bound = 64 * N * (d + 1)
      )
    ),
    class = "cdid_imputation_design_owner"
  )
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  if (actual > bound) imputation_abort("Imputation design owner bound underestimated.", "bound")
  cp04_budget(ctrl, c(borrowed, actual), scratch, 0, actual, "retained imputation design owner")
  obj$accounting$owned_entries <- actual
  ipw_owner_seal(obj, "imputation_design_owner")
}

validate_imputation_design_owner <- function(design_owner, reader = NULL, frame_owner = NULL) {
  if (!inherits(design_owner, "cdid_imputation_design_owner") ||
      !ipw_owner_seal_ok(design_owner, "imputation_design_owner")) {
    imputation_abort("Expected a sealed `cdid_imputation_design_owner`.", "integrity")
  }
  if (!is.null(reader) && !identical(design_owner$reader_ref, reader$ref)) {
    imputation_abort("`design_owner` does not belong to the supplied `reader`.", "identity")
  }
  if (!is.null(frame_owner) && !identical(design_owner$frame_ref, frame_owner$ref)) {
    imputation_abort("`design_owner` does not belong to the supplied `frame_owner`.", "identity")
  }
  invisible(design_owner)
}

imputation_resolve_full_design <- function(reader, frame_owner,
                                           design_owner = NULL,
                                           X_frame = NULL,
                                           control = reader$control) {
  st <- reader$state
  if (!identical(st$universe, "F") || !all(st$U) || any(st$rho != 1)) {
    imputation_abort(
      'Pure imputation (`selection_est_method = "imputation"`) requires the full representative roster `universe = "F"` with `U = TRUE` and `rho = 1` on all F rows.',
      "universe"
    )
  }
  if (!is.null(design_owner)) {
    validate_imputation_design_owner(design_owner, reader = reader, frame_owner = frame_owner)
    return(list(
      X = design_owner$X,
      columns = design_owner$columns,
      source = "design_owner",
      design_ref = design_owner$ref,
      owned_entries = 0,
      borrowed_entries = design_owner$accounting$owned_entries,
      work_ops = 0
    ))
  }
  if (!is.null(X_frame)) {
    d_owner <- construct_imputation_design_owner(
      frame_owner = frame_owner,
      reader = reader,
      id = st$id,
      X = X_frame,
      control = control
    )
    return(list(
      X = d_owner$X,
      columns = d_owner$columns,
      source = "X_frame",
      design_ref = d_owner$ref,
      owned_entries = d_owner$accounting$owned_entries,
      borrowed_entries = 0,
      work_ops = d_owner$accounting$work_ops_bound
    ))
  }
  # Fallback: check whether `reader$state$Z` is already the full-F extension of `frame_owner$Z_a`
  Z_f <- st$Z
  z_cols <- colnames(frame_owner$Z_a)
  f_cols <- colnames(Z_f)
  if (!is.matrix(Z_f) || nrow(Z_f) != st$N || ncol(Z_f) < ncol(frame_owner$Z_a) ||
      !is.character(z_cols) || !is.character(f_cols) ||
      anyNA(f_cols) || any(!nzchar(f_cols)) || anyDuplicated(f_cols) ||
      !all(z_cols %in% f_cols)) {
    imputation_abort(
      paste0(
        "`reader$state$Z` does not contain the columns of `frame_owner$Z_a` (it may hold a selection-probability recipe while `frame_owner$Z_a` is A-only). ",
        "Supply full-F outcome/treatment covariates via `design_owner` / `X_frame` or construct `reader` with full-F outcome/treatment covariates."
      ),
      "design"
    )
  }
  scan_retained <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries
  scan_scratch <- 16 * st$N + 16 * min(st$N, control$max_chunk_rows) * (ncol(Z_f) + 1) + 4096
  cp04_budget(control, scan_retained, scan_scratch, 64 * st$N * (ncol(Z_f) + 1), 4096,
    "full-F imputation design identity scan")
  a_pos <- which(st$a_index > 0L)
  a_idx <- st$a_index[a_pos]
  m_cols <- match(z_cols, f_cols)
  if (length(a_pos) > 0L) for (start in seq.int(1L, length(a_pos), by = control$max_chunk_rows)) {
    ix <- seq.int(start, min(length(a_pos), start + control$max_chunk_rows - 1L))
    if (max(abs(Z_f[a_pos[ix], m_cols, drop = FALSE] - frame_owner$Z_a[a_idx[ix], , drop = FALSE])) > 0) {
      imputation_abort(
        paste0(
          "`reader$state$Z` disagrees with `frame_owner$Z_a` on repeat-observed A units (do not confuse a selection-probability matrix with full-F outcome/treatment covariates). ",
          "Supply full-F outcome/treatment covariates via `design_owner` / `X_frame`."
        ),
        "design"
      )
    }
  }
  for (start in seq.int(1L, st$N, by = min(st$N, control$max_chunk_rows))) {
    rr <- seq.int(start, min(st$N, start + control$max_chunk_rows - 1L))
    cp04_finite(Z_f[rr, , drop = FALSE], "full-F outcome/treatment design on reader")
  }
  list(
    X = Z_f,
    columns = f_cols,
    source = "reader$state$Z",
    design_ref = NULL,
    owned_entries = 0,
    borrowed_entries = 0,
    work_ops = 16 * st$N * ncol(Z_f)
  )
}

# ------------------------------------------------------------------------------
# 2. Canonical Registries: Pairs, Active Cells C_imp, Pooling Partition R_imp,
#    Shared OR Blocks B (reg/dr only), and Treatment Propensity Models (ipw/dr only)
# ------------------------------------------------------------------------------

# `families` (D60 recipe-cohort layout): the recipe of each layout cohort.
# across_cohorts pooling then forms one group per recipe (pooling partitions
# cells inside one recipe); NULL keeps the single global group unchanged.
imputation_build_catalog <- function(records, cohorts, est_method, imputation_pooling, families = NULL) {
  L_links <- length(records)
  G <- length(cohorts)

  # 1. Deduplicate pairs and order canonically by (L_e, s_calendar_index, t_calendar_index)
  pair_keys <- character()
  raw_pairs <- list()
  link_raw_pair <- integer(L_links)
  for (ell in seq_len(L_links)) {
    pr <- records[[ell]]$pair
    L_e <- or_pooled_transitions(pr)
    pr$L <- as.integer(L_e)
    pr$pair_label <- paste0("e", pr$s_calendar_index, pr$t_calendar_index)
    pk <- paste0(pr$s_calendar_index, ":", pr$t_calendar_index)
    idx <- match(pk, pair_keys)
    if (is.na(idx)) {
      idx <- length(raw_pairs) + 1L
      pair_keys[idx] <- pk
      raw_pairs[[idx]] <- pr
    }
    link_raw_pair[ell] <- idx
  }
  pair_ord <- order(
    vapply(raw_pairs, `[[`, integer(1L), "L"),
    vapply(raw_pairs, `[[`, integer(1L), "s_calendar_index"),
    vapply(raw_pairs, `[[`, integer(1L), "t_calendar_index")
  )
  pairs <- raw_pairs[pair_ord]
  for (p_id in seq_along(pairs)) pairs[[p_id]]$pair_id <- as.integer(p_id)
  raw_to_canon_pair <- match(seq_along(raw_pairs), pair_ord)
  n_pairs <- length(pairs)

  # 2. Active actual-change imputation cells C_imp = union_l {(g_l, e_l)} U {(k, e_l) : k in K_l}
  cell_mat <- matrix(FALSE, G, n_pairs)
  for (ell in seq_len(L_links)) {
    pid <- raw_to_canon_pair[link_raw_pair[ell]]
    rec <- records[[ell]]
    cell_mat[rec$target, pid] <- TRUE
    cell_mat[rec$controls, pid] <- TRUE
  }
  # Canonical cohort sort order: never-treated (0 / Inf) first, then finite cohorts ascending
  coh_order <- order(imputation_cohort_sort_key(cohorts))
  coh_labels <- imputation_cohort_label(cohorts)

  cell_specs <- list()
  cell_lookup <- matrix(0L, G, n_pairs)
  for (ci in coh_order) {
    for (pid in seq_len(n_pairs)) {
      if (!cell_mat[ci, pid]) next
      cid <- length(cell_specs) + 1L
      pr <- pairs[[pid]]
      c_lab <- coh_labels[ci]
      cell_specs[[cid]] <- list(
        cell_id = as.integer(cid),
        cohort_index = as.integer(ci),
        cohort = as.double(cohorts[ci]),
        cohort_label = c_lab,
        pair_id = as.integer(pid),
        pair = pr,
        pair_label = pr$pair_label,
        transitions = as.integer(pr$L),
        L = as.integer(pr$L),
        name = paste0("k", c_lab, ".", pr$pair_label, if (is.null(families)) "" else paste0(".r", families[ci]))
      )
      cell_lookup[ci, pid] <- as.integer(cid)
    }
  }
  n_cells <- length(cell_specs)
  active_cohorts <- coh_order[vapply(coh_order, function(ci) any(cell_mat[ci, ]), logical(1L))]

  # 3. Imputation slope-pooling partition R_imp of C_imp
  group_specs <- list()
  if (identical(imputation_pooling, "separate")) {
    for (cid in seq_len(n_cells)) {
      cs <- cell_specs[[cid]]
      gid <- cid
      cell_specs[[cid]]$group_id <- as.integer(gid)
      cell_specs[[cid]]$local_id <- 1L
      group_specs[[gid]] <- list(
        group_id = as.integer(gid),
        group_index0 = as.integer(gid - 1L),
        name = paste0("group", gid - 1L),
        pooling = "separate",
        cell_ids = as.integer(cid),
        cohorts = as.integer(cs$cohort_index)
      )
    }
  } else if (identical(imputation_pooling, "within_cohort")) {
    for (g_idx in seq_along(active_cohorts)) {
      ci <- active_cohorts[g_idx]
      cids <- which(vapply(cell_specs, `[[`, integer(1L), "cohort_index") == ci)
      for (loc in seq_along(cids)) {
        cell_specs[[cids[loc]]]$group_id <- as.integer(g_idx)
        cell_specs[[cids[loc]]]$local_id <- as.integer(loc)
      }
      group_specs[[g_idx]] <- list(
        group_id = as.integer(g_idx),
        group_index0 = as.integer(g_idx - 1L),
        name = paste0("group", g_idx - 1L),
        pooling = "within_cohort",
        cell_ids = as.integer(cids),
        cohorts = as.integer(ci)
      )
    }
  } else if (is.null(families)) {
    # across_cohorts: single global slope-pooling group across all active cells
    cids <- seq_len(n_cells)
    for (loc in cids) {
      cell_specs[[loc]]$group_id <- 1L
      cell_specs[[loc]]$local_id <- as.integer(loc)
    }
    group_specs[[1L]] <- list(
      group_id = 1L,
      group_index0 = 0L,
      name = "group0",
      pooling = "across_cohorts",
      cell_ids = as.integer(cids),
      cohorts = as.integer(active_cohorts)
    )
  } else {
    # across_cohorts inside each recipe family (D60): one group per recipe.
    cell_fam <- families[vapply(cell_specs, `[[`, integer(1L), "cohort_index")]
    fams <- sort(unique(cell_fam))
    for (g_idx in seq_along(fams)) {
      cids <- which(cell_fam == fams[g_idx])
      for (loc in seq_along(cids)) {
        cell_specs[[cids[loc]]]$group_id <- as.integer(g_idx)
        cell_specs[[cids[loc]]]$local_id <- as.integer(loc)
      }
      group_specs[[g_idx]] <- list(
        group_id = as.integer(g_idx),
        group_index0 = as.integer(g_idx - 1L),
        name = paste0("group", g_idx - 1L),
        pooling = "across_cohorts",
        recipe = as.integer(fams[g_idx]),
        cell_ids = as.integer(cids),
        cohorts = as.integer(active_cohorts[families[active_cohorts] == fams[g_idx]])
      )
    }
  }

  # 4. Link records with target & control cell mappings
  links <- vector("list", L_links)
  for (ell in seq_len(L_links)) {
    rec <- records[[ell]]
    pid <- raw_to_canon_pair[link_raw_pair[ell]]
    pr <- pairs[[pid]]
    k_ord <- rec$controls[order(imputation_cohort_sort_key(cohorts[rec$controls]))]
    c_plus <- cell_lookup[rec$target, pid]
    c_minus <- cell_lookup[k_ord, pid]
    links[[ell]] <- list(
      link_index = as.integer(ell),
      link_id = rec$link_id,
      source_functional_key = rec$source_functional_key,
      functional_key = rec$functional_key,
      pair_id = as.integer(pid),
      pair = pr,
      L = as.integer(pr$L),
      target_cohort_index = as.integer(rec$target),
      target_cohort = as.double(cohorts[rec$target]),
      target_cohort_label = coh_labels[rec$target],
      control_cohort_indices = as.integer(k_ord),
      control_cohorts = as.double(cohorts[k_ord]),
      control_cohort_labels = coh_labels[k_ord],
      target_cell_id = as.integer(c_plus),
      control_cell_ids = as.integer(c_minus),
      or_block_id = NA_integer_,
      trt_id = NA_integer_
    )
  }

  # 5. Deduplicated shared counterfactual OR blocks B (built ONLY when est_method %in% c("reg", "dr"))
  or_blocks <- list()
  if (est_method %in% c("reg", "dr")) {
    bkeys <- character()
    raw_blocks <- list()
    link_raw_block <- integer(L_links)
    for (ell in seq_len(L_links)) {
      lk <- links[[ell]]
      bk <- paste0(lk$pair_id, "|", paste(lk$control_cohort_indices, collapse = ","))
      bid <- match(bk, bkeys)
      if (is.na(bid)) {
        bid <- length(raw_blocks) + 1L
        bkeys[bid] <- bk
        raw_blocks[[bid]] <- list(
          pair_id = lk$pair_id,
          pair = lk$pair,
          L = lk$L,
          control_cohort_indices = lk$control_cohort_indices,
          control_cohorts = lk$control_cohorts,
          control_cohort_labels = lk$control_cohort_labels,
          control_cell_ids = lk$control_cell_ids,
          key = bk
        )
      }
      link_raw_block[ell] <- bid
    }
    b_ord <- order(
      vapply(raw_blocks, `[[`, integer(1L), "pair_id"),
      vapply(raw_blocks, `[[`, character(1L), "key")
    )
    or_blocks <- raw_blocks[b_ord]
    raw_to_canon_block <- match(seq_along(raw_blocks), b_ord)
    for (bid in seq_along(or_blocks)) {
      b <- or_blocks[[bid]]
      or_blocks[[bid]]$block_id <- as.integer(bid)
      or_blocks[[bid]]$block_index0 <- as.integer(bid - 1L)
      or_blocks[[bid]]$name <- paste0(
        b$pair$pair_label, ":K={", paste(b$control_cohort_labels, collapse = ","), "}",
        if (is.null(families)) "" else paste0(":r", families[b$control_cohort_indices[1L]])
      )
      or_blocks[[bid]]$links <- which(raw_to_canon_block[link_raw_block] == bid)
    }
    for (ell in seq_len(L_links)) {
      links[[ell]]$or_block_id <- as.integer(raw_to_canon_block[link_raw_block[ell]])
    }
  }

  # 6. Deduplicated treatment propensity models (built ONLY when est_method %in% c("ipw", "dr"))
  trt_specs <- list()
  if (est_method %in% c("ipw", "dr")) {
    tkeys <- character()
    for (ell in seq_len(L_links)) {
      lk <- links[[ell]]
      tk <- paste0(
        lk$pair_id, "|", lk$target_cohort_index, "|",
        paste(lk$control_cohort_indices, collapse = ",")
      )
      tid <- match(tk, tkeys)
      if (is.na(tid)) {
        tid <- length(trt_specs) + 1L
        tkeys[tid] <- tk
        trt_specs[[tid]] <- list(
          trt_id = as.integer(tid),
          pair_id = lk$pair_id,
          pair = lk$pair,
          target_cohort_index = lk$target_cohort_index,
          target_cohort = lk$target_cohort,
          target_cohort_label = lk$target_cohort_label,
          control_cohort_indices = lk$control_cohort_indices,
          control_cohorts = lk$control_cohorts,
          control_cohort_labels = lk$control_cohort_labels,
          all_cohort_indices = as.integer(c(lk$target_cohort_index, lk$control_cohort_indices)),
          target_cell_id = lk$target_cell_id,
          control_cell_ids = lk$control_cell_ids,
          links = integer(0L)
        )
      }
      trt_specs[[tid]]$links <- c(trt_specs[[tid]]$links, ell)
      links[[ell]]$trt_id <- as.integer(tid)
    }
  }

  list(
    pairs = pairs,
    cell_specs = cell_specs,
    cell_lookup = cell_lookup,
    active_cohorts = as.integer(active_cohorts),
    group_specs = group_specs,
    or_blocks = or_blocks,
    treatment_specs = trt_specs,
    links = links
  )
}

# ------------------------------------------------------------------------------
# 3. Core Pure-Imputation Fitter (Stage 1 WLS, Stage 2 Shared OR, Stage 3 Trt Logit,
#    Stage 4 Link Ratios & Sparse Adjoint Precomputation)
# ------------------------------------------------------------------------------

imputation_fit_core <- function(st, X_F, pair_reader, cat, recipe, or_frame_index,
                                trt_recipe, est_method, imputation_pooling,
                                glm_control, rank_tol, ctrl) {
  N <- st$N
  G <- length(st$cohorts)
  B_chunk <- min(N, ctrl$max_chunk_rows)
  p_or_raw <- recipe$p
  slope_cols_F <- if (p_or_raw > 0L) or_frame_index[recipe$slope_index] else integer(0L)

  # 0. Verify positive roster mass and repeat-observation rate a_k > 0 on every active cohort
  for (ci in cat$active_cohorts) {
    if (!is.finite(st$mass_F[ci]) || st$mass_F[ci] <= 0 ||
        !is.finite(st$a[ci]) || st$a[ci] <= 0) {
      imputation_abort(
        paste0(
          "Cohort ", imputation_cohort_label(st$cohorts[ci]),
          " has zero full-roster or repeat-observation mass on F (a_k <= 0); synthetic imputation weights require positive repeat mass."
        ),
        "support"
      )
    }
  }

  # Precompute full-F cohort colmax for alias verification
  cohort_colmax_F <- vector("list", G)
  for (ci in cat$active_cohorts) {
    first <- st$cohort_start[ci]
    last <- st$cohort_end[ci]
    cmax <- numeric(p_or_raw)
    if (p_or_raw > 0L) {
      for (start in seq.int(first, last, by = B_chunk)) {
        rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
        Z_sl <- abs(X_F[rows, slope_cols_F, drop = FALSE])
        if (nrow(Z_sl) > 0L) cmax <- pmax(cmax, apply(Z_sl, 2L, max))
      }
    }
    cohort_colmax_F[[ci]] <- cmax
  }

  # ----------------------------------------------------------------------------
  # Stage 1: Actual-change imputation WLS on observed pairs with FIXED unit weights w_i
  # ----------------------------------------------------------------------------
  n_cells <- length(cat$cell_specs)
  cell_supports <- vector("list", n_cells)
  for (cid in seq_len(n_cells)) {
    cs <- cat$cell_specs[[cid]]
    ci <- cs$cohort_index
    first <- st$cohort_start[ci]
    last <- st$cohort_end[ci]
    acc <- or_pooled_support_start(recipe, cs$L)
    for (start in seq.int(first, last, by = B_chunk)) {
      rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
      hit <- pair_reader(rows, cs$pair, TRUE)
      obs <- which(hit$S == 1)
      if (!length(obs)) next
      r_obs <- rows[obs]
      w_obs <- st$w[r_obs]
      dy_obs <- hit$dy[obs]
      Z_obs <- X_F[r_obs, or_frame_index, drop = FALSE]
      acc <- or_pooled_support_add(acc, recipe, Z_obs, dy_obs, w_obs, "imputation")
    }
    if (acc$n < 1L || !is.finite(acc$weight_mass) || acc$weight_mass <= 0 ||
        acc$weight_mass / N <= 0) {
      imputation_abort(
        paste0("Actual-change imputation cell `", cs$name, "` has zero observed-pair support."),
        "support"
      )
    }
    cell_supports[[cid]] <- acc
  }

  # Fit each pooling group r in R_imp with full-F prediction alias verification
  n_groups <- length(cat$group_specs)
  group_fits <- vector("list", n_groups)
  cells_fitted <- vector("list", n_cells)

  for (gid in seq_len(n_groups)) {
    gs <- cat$group_specs[[gid]]
    cids <- gs$cell_ids
    B_r <- length(cids)

    # Build supports for or_pooled_fit_blocks:
    # 1..B_r: observed-pair training supports for cells in group r
    # (B_r+1)..(2*B_r): full-F cohort prediction supports for cells in group r
    supports_r <- vector("list", 2L * B_r)
    blocks_r <- vector("list", B_r)
    links_r <- vector("list", B_r)
    for (loc in seq_len(B_r)) {
      cid <- cids[loc]
      cs <- cat$cell_specs[[cid]]
      acc_c <- cell_supports[[cid]]
      supports_r[[loc]] <- list(
        summary = acc_c$summary,
        scale = 1.0,
        colmax = acc_c$colmax
      )
      supports_r[[B_r + loc]] <- list(
        summary = NULL,
        scale = 1.0,
        colmax = cohort_colmax_F[[cs$cohort_index]]
      )
      blocks_r[[loc]] <- list(
        name = cs$name,
        support_ids = as.integer(loc),
        L = as.integer(cs$L)
      )
      links_r[[loc]] <- list(
        regression = as.integer(loc),
        support_id = as.integer(B_r + loc)
      )
    }

    visit_group_rows <- function(sid, visit) {
      if (sid <= B_r) {
        cid <- cids[sid]
        cs <- cat$cell_specs[[cid]]
        ci <- cs$cohort_index
        first <- st$cohort_start[ci]
        last <- st$cohort_end[ci]
        for (start in seq.int(first, last, by = B_chunk)) {
          rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
          hit <- pair_reader(rows, cs$pair, FALSE)
          obs <- which(hit$S == 1)
          if (length(obs)) visit(X_F[rows[obs], or_frame_index, drop = FALSE])
        }
      } else {
        cid <- cids[sid - B_r]
        cs <- cat$cell_specs[[cid]]
        ci <- cs$cohort_index
        first <- st$cohort_start[ci]
        last <- st$cohort_end[ci]
        for (start in seq.int(first, last, by = B_chunk)) {
          rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
          visit(X_F[rows, or_frame_index, drop = FALSE])
        }
      }
      invisible(NULL)
    }

    fit_r <- tryCatch(or_pooled_fit_blocks(
      recipe = recipe,
      supports = supports_r,
      blocks = blocks_r,
      links = links_r,
      rank_tol = rank_tol,
      N = N,
      visit_rows = visit_group_rows,
      route = "imputation"
    ), cdid_or_wls_rank_error = function(e) {
      e$message <- paste0(conditionMessage(e), " This computation failure does not establish population ATT nonidentification.")
      e$is_identifiability_claim <- FALSE
      stop(e)
    })
    keep_r <- fit_r$keep
    p_r <- length(keep_r)
    kept_cols_F <- if (p_r > 0L) slope_cols_F[keep_r] else integer(0L)
    kept_names <- if (p_r > 0L) recipe$slope_columns[keep_r] else character(0L)
    beta_r <- as.double(fit_r$beta)
    names(beta_r) <- kept_names
    cp04_finite(beta_r, "imputation group slope coefficients")

    alpha_vec <- numeric(B_r)
    for (loc in seq_len(B_r)) {
      cid <- cids[loc]
      cs <- cat$cell_specs[[cid]]
      acc_c <- cell_supports[[cid]]
      reg_loc <- fit_r$regressions[[loc]]
      alpha_c <- if (recipe$intercept) as.double(reg_loc$alpha) else 0.0
      cp04_finite(alpha_c, "imputation cell intercept")
      alpha_vec[loc] <- alpha_c
      w_sum_c <- acc_c$weight_mass
      mass_c <- w_sum_c / N
      mean_d_c <- if (recipe$intercept && p_r > 0L) {
        as.double(reg_loc$mean_d)
      } else {
        numeric(p_r)
      }
      names(mean_d_c) <- kept_names
      mean_y_c <- acc_c$sumY / w_sum_c

      # Compute full-roster F cell moments M_{k,e}(1), Y_{k,e}(1), M_{k,e}(z), Y_{k,e}(z)
      ci <- cs$cohort_index
      a_k <- st$a[ci]
      first <- st$cohort_start[ci]
      last <- st$cohort_end[ci]
      M_one <- 0.0
      Y_one <- 0.0
      M_z <- numeric(p_or_raw)
      Y_z <- numeric(p_or_raw)
      for (start in seq.int(first, last, by = B_chunk)) {
        rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
        w_syn <- st$w[rows] * a_k
        Z_r <- if (p_r > 0L) X_F[rows, kept_cols_F, drop = FALSE] else matrix(0.0, length(rows), 0L)
        jhat_rows <- if (p_r > 0L) {
          alpha_c + as.vector((cs$L * Z_r) %*% beta_r)
        } else {
          rep(alpha_c, length(rows))
        }
        M_one <- M_one + sum(w_syn)
        Y_one <- Y_one + sum(w_syn * jhat_rows)
        if (p_or_raw > 0L) {
          Z_all_rows <- X_F[rows, slope_cols_F, drop = FALSE]
          M_z <- M_z + as.vector(crossprod(w_syn, Z_all_rows))
          Y_z <- Y_z + as.vector(crossprod(w_syn * jhat_rows, Z_all_rows))
        }
      }
      M_one <- M_one / N
      Y_one <- Y_one / N
      M_z <- M_z / N
      Y_z <- Y_z / N
      names(M_z) <- names(Y_z) <- recipe$slope_columns
      M_d_kept <- if (p_r > 0L) cs$L * M_z[keep_r] else numeric(0L)
      cp04_finite(c(M_one, Y_one, M_z, Y_z), "imputation full-F cell moments")

      cells_fitted[[cid]] <- c(
        cs,
        list(
          group_id = as.integer(gid),
          local_id = as.integer(loc),
          n_pair = as.integer(acc_c$n),
          weight_sum = w_sum_c,
          mass = mass_c,
          alpha = alpha_c,
          mean_d = mean_d_c,
          mean_y = mean_y_c,
          a_k = a_k,
          M_one = M_one,
          Y_one = Y_one,
          M_z = M_z,
          Y_z = Y_z,
          M_d_kept = M_d_kept
        )
      )
    }
    names(alpha_vec) <- vapply(cat$cell_specs[cids], `[[`, character(1L), "name")

    group_fits[[gid]] <- c(
      gs,
      list(
        keep = as.integer(keep_r),
        p = as.integer(p_r),
        kept_columns = kept_names,
        kept_cols_F = as.integer(kept_cols_F),
        beta = beta_r,
        alpha = if (recipe$intercept) alpha_vec else numeric(0L),
        fit = fit_r$fit,
        summary = fit_r$summary,
        basis = fit_r$basis,
        rss = fit_r$rss,
        objective = fit_r$objective
      )
    )
  }

  # ----------------------------------------------------------------------------
  # Stage 2: Shared counterfactual control OR across B (ONLY when est_method %in% c("reg", "dr"))
  # ----------------------------------------------------------------------------
  shared_or <- NULL
  if (est_method %in% c("reg", "dr")) {
    n_blocks <- length(cat$or_blocks)
    # Build per-cell full-F synthetic-weight OR support summaries for all control cells in B
    used_ctrl_cells <- sort(unique(unlist(lapply(cat$or_blocks, `[[`, "control_cell_ids"))))
    ctrl_cell_or_acc <- vector("list", n_cells)
    for (cid in used_ctrl_cells) {
      cf <- cells_fitted[[cid]]
      gf <- group_fits[[cf$group_id]]
      ci <- cf$cohort_index
      a_k <- st$a[ci]
      first <- st$cohort_start[ci]
      last <- st$cohort_end[ci]
      acc_or <- or_pooled_support_start(recipe, cf$L)
      for (start in seq.int(first, last, by = B_chunk)) {
        rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
        w_syn <- st$w[rows] * a_k
        Z_r <- if (gf$p > 0L) X_F[rows, gf$kept_cols_F, drop = FALSE] else matrix(0.0, length(rows), 0L)
        jhat_rows <- if (gf$p > 0L) {
          cf$alpha + as.vector((cf$L * Z_r) %*% gf$beta)
        } else {
          rep(cf$alpha, length(rows))
        }
        acc_or <- or_pooled_support_add(
          acc_or, recipe, X_F[rows, or_frame_index, drop = FALSE], jhat_rows, w_syn, "imputation_or"
        )
      }
      ctrl_cell_or_acc[[cid]] <- acc_or
    }

    # Map supports for or_pooled_fit_blocks: one support per active cell cid in C_imp
    supports_or <- lapply(seq_len(n_cells), function(cid) {
      cs <- cells_fitted[[cid]]
      acc_or <- ctrl_cell_or_acc[[cid]]
      list(
        summary = if (!is.null(acc_or)) acc_or$summary else NULL,
        scale = 1.0,
        colmax = cohort_colmax_F[[cs$cohort_index]]
      )
    })
    blocks_or_spec <- lapply(cat$or_blocks, function(b) {
      list(
        name = b$name,
        support_ids = as.integer(b$control_cell_ids),
        L = as.integer(b$L)
      )
    })
    links_or_spec <- lapply(cat$links, function(lk) {
      list(
        regression = as.integer(lk$or_block_id),
        support_id = as.integer(lk$target_cell_id)
      )
    })
    visit_or_rows <- function(sid, visit) {
      cs <- cells_fitted[[sid]]
      ci <- cs$cohort_index
      first <- st$cohort_start[ci]
      last <- st$cohort_end[ci]
      for (start in seq.int(first, last, by = B_chunk)) {
        rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
        visit(X_F[rows, or_frame_index, drop = FALSE])
      }
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
      route = "imputation_or"
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
    cp04_finite(c(alpha_or, beta_or), "shared counterfactual OR coefficients")

    # Precompute per-block masses, mean_d, and upstream forcing from a_k and j_hat_{k,e_b}
    or_block_fits <- vector("list", n_blocks)
    for (bid in seq_len(n_blocks)) {
      b_spec <- cat$or_blocks[[bid]]
      reg_b <- fit_or$regressions[[bid]]
      cids_b <- b_spec$control_cell_ids
      w_sum_b <- sum(vapply(cids_b, function(cid) ctrl_cell_or_acc[[cid]]$weight_mass, numeric(1L)))
      mass_b <- w_sum_b / N
      if (!is.finite(mass_b) || mass_b <= 0) {
        imputation_abort(paste0("Shared OR block `", b_spec$name, "` has nonpositive synthetic control mass."), "support")
      }
      mean_d_b <- if (recipe$intercept && p_or > 0L) as.double(reg_b$mean_d) else numeric(p_or)
      names(mean_d_b) <- kept_or_names
      alpha_b <- if (recipe$intercept) as.double(reg_b$alpha) else 0.0

      # Precompute forcing for each control cohort k in K_b (corresponding to cell cid in cids_b)
      n_kb <- length(cids_b)
      a_force_alpha <- numeric(n_kb)
      a_force_beta <- matrix(0.0, p_or, n_kb)
      j_force_alpha_alpha <- numeric(n_kb)
      j_force_alpha_beta <- vector("list", n_kb)
      j_force_beta_alpha <- matrix(0.0, p_or, n_kb)
      j_force_beta_beta <- vector("list", n_kb)

      for (idx_k in seq_len(n_kb)) {
        cid <- cids_b[idx_k]
        cf <- cells_fitted[[cid]]
        gf <- group_fits[[cf$group_id]]
        ci <- cf$cohort_index
        a_k <- st$a[ci]
        L_b <- b_spec$L
        first <- st$cohort_start[ci]
        last <- st$cohort_end[ci]

        sum_w_r <- 0.0
        sum_w_r_dtilde <- numeric(p_or)
        sum_wsyn_dtilde <- numeric(p_or)
        sum_wsyn_dtilde_dimp <- matrix(0.0, p_or, gf$p)

        for (start in seq.int(first, last, by = B_chunk)) {
          rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
          w_raw <- st$w[rows]
          w_syn <- w_raw * a_k
          Z_imp <- if (gf$p > 0L) X_F[rows, gf$kept_cols_F, drop = FALSE] else matrix(0.0, length(rows), 0L)
          d_imp <- L_b * Z_imp
          jhat_rows <- if (gf$p > 0L) {
            cf$alpha + as.vector(d_imp %*% gf$beta)
          } else {
            rep(cf$alpha, length(rows))
          }

          Z_or <- if (p_or > 0L) X_F[rows, kept_or_cols_F, drop = FALSE] else matrix(0.0, length(rows), 0L)
          d_or <- L_b * Z_or
          mhat_rows <- if (p_or > 0L) {
            alpha_b + as.vector(d_or %*% beta_or)
          } else {
            rep(alpha_b, length(rows))
          }
          r_or <- jhat_rows - mhat_rows

          sum_w_r <- sum_w_r + sum(w_raw * r_or)
          if (p_or > 0L) {
            d_tilde <- if (recipe$intercept) {
               sweep(d_or, 2L, mean_d_b, FUN = "-")
            } else {
              d_or
            }
            sum_w_r_dtilde <- sum_w_r_dtilde + as.vector(crossprod(w_raw * r_or, d_tilde))
            sum_wsyn_dtilde <- sum_wsyn_dtilde + as.vector(crossprod(w_syn, d_tilde))
            if (gf$p > 0L) {
              sum_wsyn_dtilde_dimp <- sum_wsyn_dtilde_dimp + crossprod(w_syn * d_tilde, d_imp)
            }
          }
        }
        a_force_alpha[idx_k] <- sum_w_r / N
        a_force_beta[, idx_k] <- sum_w_r_dtilde / N
        j_force_alpha_alpha[idx_k] <- if (recipe$intercept) cf$M_one else 0.0
        j_force_alpha_beta[[idx_k]] <- as.double(cf$M_d_kept)
        j_force_beta_alpha[, idx_k] <- if (recipe$intercept) (sum_wsyn_dtilde / N) else numeric(p_or)
        j_force_beta_beta[[idx_k]] <- sum_wsyn_dtilde_dimp / N
      }

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
          j_force_beta_beta = j_force_beta_beta
        )
      )
    }

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
      fit = fit_or$fit,
      summary = fit_or$summary,
      basis = fit_or$basis,
      rss = fit_or$rss,
      objective = fit_or$objective
    )
  }

  # ----------------------------------------------------------------------------
  # Stage 3: Synthetic-weight treatment propensity logits on F (ONLY when est_method %in% c("ipw", "dr"))
  # ----------------------------------------------------------------------------
  trt_fits <- NULL
  if (est_method %in% c("ipw", "dr")) {
    n_trt <- length(cat$treatment_specs)
    trt_fits <- vector("list", n_trt)
    p_trt <- trt_recipe$p
    trt_cols_F <- trt_recipe$index

    for (tid in seq_len(n_trt)) {
      ts <- cat$treatment_specs[[tid]]
      # Collect all F rows in {g} U K in canonical supplied row order
      risk_rows <- sort(unlist(lapply(ts$all_cohort_indices, function(ci) {
        st$cohort_order[seq.int(st$cohort_start[ci], st$cohort_end[ci])]
      }), use.names = FALSE))
      n_risk <- length(risk_rows)
      ci_risk <- st$cohort_index[risk_rows]
      y_risk <- as.double(ci_risk == ts$target_cohort_index)
      w_raw <- st$w[risk_rows]
      w_syn <- w_raw * st$a[ci_risk]
      X_trt <- X_F[risk_rows, trt_cols_F, drop = FALSE]

      lfit <- cp04_logit(X_trt, y_risk, w_syn, glm_control, N, ts)
      gamma_hat <- as.double(lfit$coefficients)
      names(gamma_hat) <- paste0("gamma", seq.int(0L, p_trt - 1L))
      lv <- ipw_logit_values(X_trt, gamma_hat)
      cp04_finite(c(gamma_hat, lv$prob, lv$eta, lv$odds), "imputation treatment propensity fit")

      # Repeat-rate forcing C_{gamma, k} = H_N[w I_k (I_{k=g} - h_hat) X^trt]
      m_coh <- length(ts$all_cohort_indices)
      C_a <- matrix(0.0, p_trt, m_coh)
      raw_score <- w_raw * (y_risk - lv$prob) * X_trt
      for (idx_k in seq_len(m_coh)) {
        ci_k <- ts$all_cohort_indices[idx_k]
        sub_k <- which(ci_risk == ci_k)
        if (length(sub_k)) {
          C_a[, idx_k] <- colSums(raw_score[sub_k, , drop = FALSE]) / N
        }
      }
      cp04_finite(C_a, "imputation treatment repeat-rate forcing")

      trt_fits[[tid]] <- c(
        ts,
        list(
          coefficients = gamma_hat,
          column_names = trt_recipe$columns,
          H_chol = lfit$H_chol,
          H = lfit$H,
          n_risk_F = as.integer(n_risk),
          a_forcing = C_a,
          clipping = ipw_clip_summary(prob = lv$prob, controls = y_risk == 0, weights = w_syn)
        )
      )
    }
  }

  # ----------------------------------------------------------------------------
  # Stage 4: Link Point Estimates, Partial Derivatives, and Sparse Adjoint Weights
  # ----------------------------------------------------------------------------
  L_links <- length(cat$links)
  links_fitted <- vector("list", L_links)

  # Precompute OR sensitivity vectors to imputation cells (for OR_j channel)
  # When shared_or is present, let's collect the total RHS forcing matrix from
  # each imputation cell c on (xi_alpha_b^OR, xi_beta^OR).
  for (ell in seq_len(L_links)) {
    lk <- cat$links[[ell]]
    cf_g <- cells_fitted[[lk$target_cell_id]]
    gf_g <- group_fits[[cf_g$group_id]]

    # Treated arm T_{gl} (available for all 3 est_methods)
    d_plus <- cf_g$M_one
    num_T <- cf_g$Y_one
    if (!is.finite(d_plus) || d_plus <= 0) {
      imputation_abort(paste0("Link `", lk$link_id, "` has nonpositive treated denominator d_{l,+}."), "support")
    }
    T_gl <- num_T / d_plus
    d_bar_plus_j <- cf_g$M_d_kept / d_plus
    # a_g partial of T_gl (identically 0, computed explicitly for completeness)
    a_coeff_T_g <- ((cf_g$Y_one / cf_g$a_k) - T_gl * (cf_g$M_one / cf_g$a_k)) / d_plus
    cp04_finite(c(d_plus, num_T, T_gl, d_bar_plus_j, a_coeff_T_g), "link treated imputation summary")

    # Treated OR baseline B_{gl} (when shared_or is fitted)
    num_Bg <- NA_real_
    B_gl <- NA_real_
    mu_reg <- NA_real_
    d_bar_plus_or <- numeric(0L)
    a_coeff_Bg_g <- 0.0
    if (!is.null(shared_or)) {
      b_fit <- shared_or$blocks[[lk$or_block_id]]
      p_or <- shared_or$p
      d_bar_plus_or <- if (p_or > 0L) {
        (lk$L * cf_g$M_z[shared_or$keep]) / d_plus
      } else {
        numeric(0L)
      }
      names(d_bar_plus_or) <- shared_or$kept_columns
      B_gl <- (if (recipe$intercept) b_fit$alpha else 0.0) +
        (if (p_or > 0L) sum(d_bar_plus_or * shared_or$beta) else 0.0)
      num_Bg <- d_plus * B_gl
      mu_reg <- T_gl - B_gl
      a_coeff_Bg_g <- ((num_Bg / cf_g$a_k) - B_gl * (d_plus / cf_g$a_k)) / d_plus
      cp04_finite(c(num_Bg, B_gl, mu_reg, d_bar_plus_or), "link treated OR baseline summary")
    }

    # Control IPW arm C_{ol} and (when DR) Control OR baseline B_{ol}
    d_minus <- NA_real_
    num_C <- NA_real_
    C_ol <- NA_real_
    mu_ipw <- NA_real_
    num_Bo <- NA_real_
    B_ol <- NA_real_
    mu_dr <- NA_real_
    d_bar_minus_or <- numeric(0L)
    d_bar_diff_or <- numeric(0L)
    a_coeff_C_k <- numeric(0L)
    a_coeff_Bo_k <- numeric(0L)
    m_minus_k <- numeric(0L)
    d_bar_minus_j_k <- list()
    denminus_gamma <- numeric(0L)
    g_C_gamma <- numeric(0L)
    g_Bo_gamma <- numeric(0L)
    g_dr_gamma <- numeric(0L)

    if (!is.null(trt_fits)) {
      tfit <- trt_fits[[lk$trt_id]]
      p_trt <- trt_recipe$p
      trt_cols_F <- trt_recipe$index
      n_ctrl_coh <- length(lk$control_cohort_indices)

      sum_wm <- 0.0
      sum_wm_j <- 0.0
      sum_wm_m <- 0.0
      sum_wm_dor <- if (!is.null(shared_or)) numeric(shared_or$p) else numeric(0L)
      sum_wm_by_k <- numeric(n_ctrl_coh)
      sum_wm_dimp_by_k <- vector("list", n_ctrl_coh)
      sum_w_odds_j_by_k <- numeric(n_ctrl_coh)
      sum_w_odds_by_k <- numeric(n_ctrl_coh)
      sum_w_odds_m_by_k <- numeric(n_ctrl_coh)
      sum_wdot_j_X <- numeric(p_trt)
      sum_wdot_X <- numeric(p_trt)
      sum_wdot_m_X <- numeric(p_trt)

      b_fit <- if (!is.null(shared_or)) shared_or$blocks[[lk$or_block_id]] else NULL
      p_or <- if (!is.null(shared_or)) shared_or$p else 0L

      for (idx_k in seq_len(n_ctrl_coh)) {
        ci_k <- lk$control_cohort_indices[idx_k]
        cid_k <- lk$control_cell_ids[idx_k]
        cf_k <- cells_fitted[[cid_k]]
        gf_k <- group_fits[[cf_k$group_id]]
        a_k <- st$a[ci_k]
        first <- st$cohort_start[ci_k]
        last <- st$cohort_end[ci_k]
        sum_dimp_k <- numeric(gf_k$p)

        for (start in seq.int(first, last, by = B_chunk)) {
          rows <- st$cohort_order[seq.int(start, min(last, start + B_chunk - 1L))]
          w_raw <- st$w[rows]
          w_syn <- w_raw * a_k
          X_trt_rows <- X_F[rows, trt_cols_F, drop = FALSE]
          lv <- ipw_logit_values(X_trt_rows, tfit$coefficients)
          odds <- lv$odds
          wm <- w_syn * odds
          # d wm / d eta = wm on active rows; zero where the odds cap binds.
          wdot <- wm * lv$active

          Z_imp <- if (gf_k$p > 0L) X_F[rows, gf_k$kept_cols_F, drop = FALSE] else matrix(0.0, length(rows), 0L)
          d_imp <- lk$L * Z_imp
          jhat_rows <- if (gf_k$p > 0L) {
            cf_k$alpha + as.vector(d_imp %*% gf_k$beta)
          } else {
            rep(cf_k$alpha, length(rows))
          }

          sum_wm <- sum_wm + sum(wm)
          sum_wm_j <- sum_wm_j + sum(wm * jhat_rows)
          sum_wm_by_k[idx_k] <- sum_wm_by_k[idx_k] + sum(wm)
          if (gf_k$p > 0L) {
            sum_dimp_k <- sum_dimp_k + as.vector(crossprod(wm, d_imp))
          }
          sum_w_odds_j_by_k[idx_k] <- sum_w_odds_j_by_k[idx_k] + sum(w_raw * odds * jhat_rows)
          sum_w_odds_by_k[idx_k] <- sum_w_odds_by_k[idx_k] + sum(w_raw * odds)
          sum_wdot_j_X <- sum_wdot_j_X + as.vector(crossprod(wdot * jhat_rows, X_trt_rows))
          sum_wdot_X <- sum_wdot_X + as.vector(crossprod(wdot, X_trt_rows))

          if (!is.null(shared_or)) {
            Z_or <- if (p_or > 0L) X_F[rows, shared_or$kept_cols_F, drop = FALSE] else matrix(0.0, length(rows), 0L)
            d_or <- lk$L * Z_or
            mhat_rows <- if (p_or > 0L) {
              (if (recipe$intercept) b_fit$alpha else 0.0) + as.vector(d_or %*% shared_or$beta)
            } else {
              rep(if (recipe$intercept) b_fit$alpha else 0.0, length(rows))
            }
            sum_wm_m <- sum_wm_m + sum(wm * mhat_rows)
            if (p_or > 0L) {
              sum_wm_dor <- sum_wm_dor + as.vector(crossprod(wm, d_or))
            }
            sum_w_odds_m_by_k[idx_k] <- sum_w_odds_m_by_k[idx_k] + sum(w_raw * odds * mhat_rows)
            sum_wdot_m_X <- sum_wdot_m_X + as.vector(crossprod(wdot * mhat_rows, X_trt_rows))
          }
        }
        sum_wm_dimp_by_k[[idx_k]] <- sum_dimp_k
      }

      d_minus <- sum_wm / N
      if (!is.finite(d_minus) || d_minus <= 0) {
        imputation_abort(paste0("Link `", lk$link_id, "` has nonpositive control IPW denominator d_{l,-}."), "support")
      }
      num_C <- sum_wm_j / N
      C_ol <- num_C / d_minus
      mu_ipw <- T_gl - C_ol

      m_minus_k <- (sum_wm_by_k / N) / d_minus
      d_bar_minus_j_k <- lapply(sum_wm_dimp_by_k, function(v) (v / N) / d_minus)
      a_coeff_C_k <- ((sum_w_odds_j_by_k / N) - C_ol * (sum_w_odds_by_k / N)) / d_minus
      denminus_gamma <- sum_wdot_X / N
      g_C_gamma <- ((sum_wdot_j_X / N) - C_ol * (sum_wdot_X / N)) / d_minus
      names(g_C_gamma) <- trt_recipe$columns
      cp04_finite(c(d_minus, num_C, C_ol, mu_ipw, m_minus_k, a_coeff_C_k, g_C_gamma), "link control IPW summary")

      if (!is.null(shared_or)) {
        d_bar_minus_or <- if (p_or > 0L) (sum_wm_dor / N) / d_minus else numeric(0L)
        names(d_bar_minus_or) <- shared_or$kept_columns
        d_bar_diff_or <- d_bar_plus_or - d_bar_minus_or
        num_Bo <- sum_wm_m / N
        B_ol <- (if (recipe$intercept) b_fit$alpha else 0.0) +
          (if (p_or > 0L) sum(d_bar_minus_or * shared_or$beta) else 0.0)
        mu_dr <- (T_gl - B_gl) - (C_ol - B_ol)
        a_coeff_Bo_k <- ((sum_w_odds_m_by_k / N) - B_ol * (sum_w_odds_by_k / N)) / d_minus
        g_Bo_gamma <- ((sum_wdot_m_X / N) - B_ol * (sum_wdot_X / N)) / d_minus
        g_dr_gamma <- g_C_gamma - g_Bo_gamma
        names(g_Bo_gamma) <- names(g_dr_gamma) <- trt_recipe$columns
        cp04_finite(c(num_Bo, B_ol, mu_dr, d_bar_minus_or, a_coeff_Bo_k, g_Bo_gamma, g_dr_gamma), "link DR summary")
      }
    }

    mu_primary <- switch(
      est_method,
      ipw = mu_ipw,
      reg = mu_reg,
      dr = mu_dr
    )
    cp04_finite(mu_primary, "link primary pure-imputation estimate")

    links_fitted[[ell]] <- c(
      lk,
      list(
        mu = mu_primary,
        mu_ipw = mu_ipw,
        mu_or = mu_reg,
        mu_reg = mu_reg,
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
        d_bar_plus_j = d_bar_plus_j,
        d_bar_plus_or = d_bar_plus_or,
        d_bar_minus_or = d_bar_minus_or,
        d_bar_diff_or = d_bar_diff_or,
        a_coeff_T_g = a_coeff_T_g,
        a_coeff_Bg_g = a_coeff_Bg_g,
        a_coeff_C_k = a_coeff_C_k,
        a_coeff_Bo_k = a_coeff_Bo_k,
        m_minus_k = m_minus_k,
        d_bar_minus_j_k = d_bar_minus_j_k,
        denminus_gamma = denminus_gamma,
        g_C_gamma = g_C_gamma,
        g_Bo_gamma = g_Bo_gamma,
        g_dr_gamma = g_dr_gamma
      )
    )
  }

  # ----------------------------------------------------------------------------
  # Stage 5: Precompute Sparse Adjoint Pullback Weights on R_imp
  # ----------------------------------------------------------------------------
  # Any linear functional sum_{c in C_imp} (lambda_c * xi_{alpha_c}^imp + nu_c' * xi_{beta_{r(c)}}^imp)
  # pulls back to observed-pair scores s_{i,c} * (tilde_lambda_c + tilde_d_{i,c}' tilde_nu_{r(c)})
  # where tilde_nu_r = (S_r^imp)^{-1} sum_{c in C_r} (nu_c - lambda_c * mean_d_c)
  # and tilde_lambda_c = lambda_c / mass_c (if intercept = TRUE, else 0).
  solve_imp_adjoint <- function(lambda_vec, nu_list) {
    tilde_lambda <- numeric(n_cells)
    tilde_nu <- vector("list", n_groups)
    for (gid in seq_len(n_groups)) {
      gf <- group_fits[[gid]]
      cids <- gf$cell_ids
      rhs_r <- numeric(gf$p)
      for (cid in cids) {
        cf <- cells_fitted[[cid]]
        lam_c <- if (recipe$intercept) lambda_vec[cid] else 0.0
        if (recipe$intercept) tilde_lambda[cid] <- lam_c / cf$mass
        if (gf$p > 0L) {
          rhs_r <- rhs_r + nu_list[[cid]] - lam_c * cf$mean_d
        }
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

  for (ell in seq_len(L_links)) {
    lf <- links_fitted[[ell]]
    # 1. Treated arm T_gl j_direct weights
    lam_T <- numeric(n_cells)
    nu_T <- lapply(cells_fitted, function(cf) numeric(group_fits[[cf$group_id]]$p))
    if (recipe$intercept) lam_T[lf$target_cell_id] <- 1.0
    nu_T[[lf$target_cell_id]] <- lf$d_bar_plus_j
    links_fitted[[ell]]$adjoint_T <- solve_imp_adjoint(lam_T, nu_T)

    # 2. Control IPW arm C_ol j_direct weights (if trt_fits present)
    if (!is.null(trt_fits)) {
      lam_C <- numeric(n_cells)
      nu_C <- lapply(cells_fitted, function(cf) numeric(group_fits[[cf$group_id]]$p))
      for (idx_k in seq_along(lf$control_cell_ids)) {
        cid_k <- lf$control_cell_ids[idx_k]
        if (recipe$intercept) lam_C[cid_k] <- lam_C[cid_k] + lf$m_minus_k[idx_k]
        nu_C[[cid_k]] <- nu_C[[cid_k]] + lf$d_bar_minus_j_k[[idx_k]]
      }
      links_fitted[[ell]]$adjoint_C <- solve_imp_adjoint(lam_C, nu_C)
    }
  }

  list(
    cells = cells_fitted,
    groups = group_fits,
    shared_or = shared_or,
    treatments = trt_fits,
    links = links_fitted
  )
}

# ------------------------------------------------------------------------------
# 4. Chunked Complete Original-Unit IF Evaluation Kernel
# ------------------------------------------------------------------------------

imputation_if_needed_sets <- function(core_fit, links, est_method, details = FALSE) {
  if (details) {
    cids <- seq_along(core_fit$cells)
    gids <- seq_along(core_fit$groups)
    bids <- if (!is.null(core_fit$shared_or)) seq_len(core_fit$shared_or$n_blocks) else integer(0L)
    tids <- if (!is.null(core_fit$treatments)) seq_along(core_fit$treatments) else integer(0L)
    return(list(cids = cids, gids = gids, bids = bids, tids = tids))
  }
  l_sub <- core_fit$links[links]
  cids_link <- sort(unique(c(
    vapply(l_sub, `[[`, integer(1L), "target_cell_id"),
    if (est_method %in% c("ipw", "dr")) unlist(lapply(l_sub, `[[`, "control_cell_ids")) else integer(0L)
  )))
  bids <- if (est_method %in% c("reg", "dr") && !is.null(core_fit$shared_or)) {
    if (core_fit$shared_or$p > 0L) {
      seq_len(core_fit$shared_or$n_blocks)
    } else {
      sort(unique(vapply(l_sub, `[[`, integer(1L), "or_block_id")))
    }
  } else {
    integer(0L)
  }
  cids_or <- if (length(bids)) {
    sort(unique(unlist(lapply(core_fit$shared_or$blocks[bids], `[[`, "control_cell_ids"))))
  } else {
    integer(0L)
  }
  cids_need <- sort(unique(c(cids_link, cids_or)))
  gids <- sort(unique(vapply(core_fit$cells[cids_need], `[[`, integer(1L), "group_id")))
  # Any group with p_r > 0 needs all of its member cells to form h_{i,beta_r}^imp
  cids_all <- sort(unique(c(
    cids_need,
    unlist(lapply(core_fit$groups[gids], function(gf) if (gf$p > 0L) gf$cell_ids else integer(0L)))
  )))
  tids <- if (est_method %in% c("ipw", "dr") && !is.null(core_fit$treatments)) {
    sort(unique(vapply(l_sub, `[[`, integer(1L), "trt_id")))
  } else {
    integer(0L)
  }
  list(cids = cids_all, gids = gids, bids = bids, tids = tids)
}

imputation_evaluate_if_chunk_core <- function(st, X_F, pair_reader, core_fit, recipe,
                                              or_frame_index, trt_recipe, est_method,
                                              rows, links, scale = "frame",
                                              details = FALSE, sets = NULL) {
  B <- length(rows)
  N <- st$N
  G <- length(st$cohorts)
  p_or_raw <- recipe$p
  slope_cols_F <- if (p_or_raw > 0L) or_frame_index[recipe$slope_index] else integer(0L)
  if (is.null(sets)) sets <- imputation_if_needed_sets(core_fit, links, est_method, details)

  ci <- st$cohort_index[rows]
  w <- st$w[rows]
  A <- as.double(st$a_index[rows] > 0L)
  X_chunk <- X_F[rows, , drop = FALSE]

  # 1. Repeat-observation rate IFs xi_{i,a_k} and cohort share IFs xi_{i,tau_k}, xi_{i,eta_k}
  xi_a <- matrix(0.0, B, G)
  xi_tau <- matrix(0.0, B, G)
  xi_eta <- matrix(0.0, B, G)
  coh_labels <- imputation_cohort_label(st$cohorts)
  colnames(xi_a) <- paste0("a.k", coh_labels)
  colnames(xi_tau) <- paste0("tau.k", coh_labels)
  colnames(xi_eta) <- paste0("eta.k", coh_labels)
  for (k in seq_len(G)) {
    I_k <- as.double(ci == k)
    if (st$mass_F[k] > 0) {
      xi_a[, k] <- w * I_k * (A - st$a[k]) / st$mass_F[k]
    }
    xi_tau[, k] <- w * (I_k - st$tau[k]) / st$mean_w
    xi_eta[, k] <- w * A * (I_k - st$eta[k]) / st$mean_wA
  }

  # 2. Actual-change imputation cell scores s_{i,c} = w_i I_{ik} S_{ie} (Delta_{ie} - j_hat_{k,e}(X_i))
  n_cells <- length(core_fit$cells)
  n_groups <- length(core_fit$groups)
  read_cell_score <- function(cid) {
    cf <- core_fit$cells[[cid]]
    gf <- core_fit$groups[[cf$group_id]]
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

  # Solve per-group slope IFs xi_{i,beta_r}^imp and per-cell intercept IFs xi_{i,alpha_c}^imp
  xi_beta_imp <- if (details) vector("list", n_groups) else NULL

  if (details) for (gid in seq_len(n_groups)) {
    gf <- core_fit$groups[[gid]]
    if (!gid %in% sets$gids) {
      xi_beta_imp[[gid]] <- matrix(0.0, B, gf$p, dimnames = list(NULL, gf$kept_columns))
      next
    }
    h_beta <- matrix(0.0, B, gf$p, dimnames = list(NULL, gf$kept_columns))
    if (gf$p > 0L) {
      for (cid in gf$cell_ids) {
        obs_rec <- read_cell_score(cid)
        if (is.null(obs_rec)) next
        cf <- core_fit$cells[[cid]]
        d_tilde <- if (recipe$intercept) {
          sweep(obs_rec$d_rr, 2L, cf$mean_d, FUN = "-")
        } else {
          obs_rec$d_rr
        }
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

  # Local intercepts are pulled back only when requested; group slope rows are
  # the only imputation parameter matrices retained during ordinary replay.
  alpha_if <- function(cid, obs_rec = NULL) {
    cf <- core_fit$cells[[cid]]
    gf <- core_fit$groups[[cf$group_id]]
    out <- numeric(B)
    if (!recipe$intercept) return(out)
    if (missing(obs_rec)) obs_rec <- read_cell_score(cid)
    if (!is.null(obs_rec)) out[obs_rec$rr] <- obs_rec$s_rr / cf$mass
    if (gf$p > 0L) out <- out - as.vector(xi_beta_imp[[cf$group_id]] %*% cf$mean_d)
    out
  }

  # Sparse adjoint pullback evaluator on the chunk
  eval_adjoint_chunk <- function(adj) {
    out_vec <- numeric(B)
    for (cid in sets$cids) {
      obs_rec <- read_cell_score(cid)
      if (is.null(obs_rec)) next
      cf <- core_fit$cells[[cid]]
      gf <- core_fit$groups[[cf$group_id]]
      t_lam <- if (recipe$intercept) adj$tilde_lambda[cid] else 0.0
      t_nu <- adj$tilde_nu[[cf$group_id]]
      val <- rep(t_lam, length(obs_rec$rr))
      if (gf$p > 0L && length(t_nu) > 0L) {
        d_tilde <- if (recipe$intercept) {
          sweep(obs_rec$d_rr, 2L, cf$mean_d, FUN = "-")
        } else {
          obs_rec$d_rr
        }
        val <- val + as.vector(d_tilde %*% t_nu)
      }
      out_vec[obs_rec$rr] <- out_vec[obs_rec$rr] + obs_rec$s_rr * val
    }
    out_vec
  }

  # Pull back one downstream scalar loading without retaining any row-by-cell
  # or row-by-group parameter IF arrays. This uses the same coupled Gram factors
  # and normalized inverse-transpose identity as the precomputed link adjoints.
  eval_imp_linear <- function(lambda_vec, nu_list) {
    tilde_lambda <- numeric(n_cells)
    tilde_nu <- vector("list", n_groups)
    for (gid in seq_len(n_groups)) {
      gf <- core_fit$groups[[gid]]
      rhs <- numeric(gf$p)
      for (cid in gf$cell_ids) {
        cf <- core_fit$cells[[cid]]
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
  zero_imp_loadings <- function() {
    lapply(core_fit$cells, function(cf) numeric(core_fit$groups[[cf$group_id]]$p))
  }

  # 3. Shared counterfactual OR IFs (decomposed into direct, a, j channels)
  or_if <- NULL
  if (!is.null(core_fit$shared_or) && length(sets$bids) > 0L) {
    sh <- core_fit$shared_or
    p_or <- sh$p
    n_b <- sh$n_blocks
    rhs_beta_dir <- matrix(0.0, B, p_or)
    rhs_beta_a <- matrix(0.0, B, p_or)
    rhs_beta_j <- matrix(0.0, B, p_or)
    alpha_bids <- if (details) seq_len(n_b) else sort(unique(vapply(core_fit$links[links], `[[`, integer(1L), "or_block_id")))
    n_alpha <- length(alpha_bids)
    rhs_alpha_dir <- matrix(0.0, B, n_alpha)
    rhs_alpha_a <- matrix(0.0, B, n_alpha)
    rhs_alpha_j <- matrix(0.0, B, n_alpha)

    for (bid in sets$bids) {
      apos <- match(bid, alpha_bids)
      bf <- sh$blocks[[bid]]
      n_kb <- length(bf$control_cell_ids)
      for (idx_k in seq_len(n_kb)) {
        cid <- bf$control_cell_ids[idx_k]
        cf <- core_fit$cells[[cid]]
        gf <- core_fit$groups[[cf$group_id]]
        ci_k <- cf$cohort_index
        a_k <- st$a[ci_k]

        # Direct score on F rows of control cohort k
        ix <- which(ci == ci_k)
        if (length(ix)) {
          w_syn_ix <- w[ix] * a_k
          Z_imp_ix <- if (gf$p > 0L) X_chunk[ix, gf$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix), 0L)
          d_imp_ix <- bf$L * Z_imp_ix
          jhat_ix <- if (gf$p > 0L) {
            cf$alpha + as.vector(d_imp_ix %*% gf$beta)
          } else {
            rep(cf$alpha, length(ix))
          }

          Z_or_ix <- if (p_or > 0L) X_chunk[ix, sh$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix), 0L)
          d_or_ix <- bf$L * Z_or_ix
          mhat_ix <- if (p_or > 0L) {
            (if (recipe$intercept) bf$alpha else 0.0) + as.vector(d_or_ix %*% sh$beta)
          } else {
            rep(if (recipe$intercept) bf$alpha else 0.0, length(ix))
          }
          r_ix <- jhat_ix - mhat_ix
          if (recipe$intercept && !is.na(apos)) {
            rhs_alpha_dir[ix, apos] <- rhs_alpha_dir[ix, apos] + w_syn_ix * r_ix
          }
          if (p_or > 0L) {
            d_tilde_ix <- if (recipe$intercept) {
              sweep(d_or_ix, 2L, bf$mean_d, FUN = "-")
            } else {
              d_or_ix
            }
            rhs_beta_dir[ix, ] <- rhs_beta_dir[ix, , drop = FALSE] + (w_syn_ix * r_ix) * d_tilde_ix
          }
        }

        # Repeat-rate forcing on xi_{i,a_k}
        xi_ak <- xi_a[, ci_k]
        if (recipe$intercept && !is.na(apos)) {
          rhs_alpha_a[, apos] <- rhs_alpha_a[, apos] + bf$a_force_alpha[idx_k] * xi_ak
        }
        if (p_or > 0L) {
          rhs_beta_a <- rhs_beta_a + outer(xi_ak, bf$a_force_beta[, idx_k])
        }

        if (details) {
          # Imputation forcing on (xi_{alpha_c}^imp, xi_{beta_r}^imp)
          xi_ac <- alpha_if(cid)
          xi_br <- xi_beta_imp[[cf$group_id]]
          if (recipe$intercept && !is.na(apos)) {
            rhs_alpha_j[, apos] <- rhs_alpha_j[, apos] + bf$j_force_alpha_alpha[idx_k] * xi_ac +
              (if (gf$p > 0L) as.vector(xi_br %*% bf$j_force_alpha_beta[[idx_k]]) else 0.0)
          }
          if (p_or > 0L) {
            if (recipe$intercept) {
              rhs_beta_j <- rhs_beta_j + outer(xi_ac, bf$j_force_beta_alpha[, idx_k])
            }
            if (gf$p > 0L) {
              rhs_beta_j <- rhs_beta_j + xi_br %*% t(bf$j_force_beta_beta[[idx_k]])
            }
          }
        }
      }
    }

    if (!details) {
      # Global shared-slope forcing collects every contributing block before
      # each scalar pullback. The row width is p_OR, independent of cell count.
      if (p_or) for (q in seq_len(p_or)) {
        lam <- numeric(n_cells)
        nu <- zero_imp_loadings()
        for (bid in sets$bids) {
          bf <- sh$blocks[[bid]]
          for (idx_k in seq_along(bf$control_cell_ids)) {
            cid <- bf$control_cell_ids[idx_k]
            if (recipe$intercept) lam[cid] <- lam[cid] + bf$j_force_beta_alpha[q, idx_k]
            nu[[cid]] <- nu[[cid]] + bf$j_force_beta_beta[[idx_k]][q, ]
          }
        }
        rhs_beta_j[, q] <- eval_imp_linear(lam, nu)
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
          bid <- alpha_bids[apos]
          bf <- sh$blocks[[bid]]
          xi_al[, apos] <- (rhs_a_mat[, apos] / bf$mass) -
            (if (p_or > 0L) as.vector(xi_b %*% bf$mean_d) else 0.0)
        }
      }
      list(alpha = xi_al, beta = xi_b)
    }

    ch_dir <- solve_or_channel(rhs_alpha_dir, rhs_beta_dir)
    ch_a <- solve_or_channel(rhs_alpha_a, rhs_beta_a)
    ch_j <- solve_or_channel(rhs_alpha_j, rhs_beta_j)
    alpha_tot <- ch_dir$alpha + ch_a$alpha + ch_j$alpha
    beta_tot <- ch_dir$beta + ch_a$beta + ch_j$beta
    theta_tot <- if (recipe$intercept) cbind(alpha_tot, beta_tot) else beta_tot
    if (ncol(theta_tot) > 0L) {
      colnames(theta_tot) <- paste0("OR.theta", seq.int(0L, ncol(theta_tot) - 1L))
    }
    or_if <- list(
      direct = ch_dir,
      a = ch_a,
      j = ch_j,
      alpha = alpha_tot,
      beta = beta_tot,
      theta = theta_tot,
      rhs_beta_score = rhs_beta_dir + rhs_beta_a + rhs_beta_j,
      rhs_alpha_score = rhs_alpha_dir + rhs_alpha_a + rhs_alpha_j
    )
  }

  # 4. Synthetic-weight treatment propensity logit IFs (decomposed into direct and a channels)
  trt_if <- NULL
  if (!is.null(core_fit$treatments) && length(sets$tids) > 0L) {
    n_trt <- length(core_fit$treatments)
    trt_if <- vector("list", n_trt)
    p_trt <- trt_recipe$p
    trt_cols_F <- trt_recipe$index

    for (tid in sets$tids) {
      tf <- core_fit$treatments[[tid]]
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
      colnames(xi_g_dir) <- colnames(xi_g_a) <- colnames(xi_g_tot) <-
        paste0("gamma", seq.int(0L, p_trt - 1L))
      trt_if[[tid]] <- list(
        direct = xi_g_dir,
        a = xi_g_a,
        gamma = xi_g_tot,
        score_direct = score_dir,
        score_a = score_a
      )
    }
  }

  # 5. Link IFs and 8-channel linear derivative decomposition
  nl <- length(links)
  link_names <- vapply(core_fit$links[links], `[[`, character(1L), "link_id")
  phi <- matrix(0.0, B, nl, dimnames = list(NULL, link_names))
  psi_ipw_mat <- matrix(NA_real_, B, nl, dimnames = list(NULL, link_names))
  psi_reg_mat <- matrix(NA_real_, B, nl, dimnames = list(NULL, link_names))
  psi_dr_mat <- matrix(NA_real_, B, nl, dimnames = list(NULL, link_names))
  link_details <- if (details) vector("list", nl) else NULL

  for (j in seq_along(links)) {
    ell <- links[j]
    lf <- core_fit$links[[ell]]
    cf_g <- core_fit$cells[[lf$target_cell_id]]
    gf_g <- core_fit$groups[[cf_g$group_id]]
    ci_g <- lf$target_cohort_index
    a_g <- st$a[ci_g]

    # --- Treated arm T_gl ---
    ix_g <- which(ci == ci_g)
    w_plus <- numeric(B)
    j_g_vec <- numeric(B)
    if (length(ix_g)) {
      w_plus[ix_g] <- w[ix_g] * a_g
      Z_imp_g <- if (gf_g$p > 0L) X_chunk[ix_g, gf_g$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix_g), 0L)
      j_g_vec[ix_g] <- if (gf_g$p > 0L) {
        cf_g$alpha + as.vector((lf$L * Z_imp_g) %*% gf_g$beta)
      } else {
        rep(cf_g$alpha, length(ix_g))
      }
    }
    psi_denplus <- w_plus + st$mass_F[ci_g] * xi_a[, ci_g]
    psi_numT <- if (details) w_plus * j_g_vec +
      (cf_g$Y_one / a_g) * xi_a[, ci_g] +
      (if (recipe$intercept) cf_g$M_one * alpha_if(lf$target_cell_id) else 0.0) +
      (if (gf_g$p > 0L) as.vector(xi_beta_imp[[cf_g$group_id]] %*% cf_g$M_d_kept) else 0.0) else NULL

    direct_T <- w_plus * (j_g_vec - lf$T_gl) / lf$d_plus
    a_T <- lf$a_coeff_T_g * xi_a[, ci_g]
    j_T_adj <- eval_adjoint_chunk(lf$adjoint_T)
    j_T <- if (details) (if (recipe$intercept) alpha_if(lf$target_cell_id) else 0.0) +
      (if (gf_g$p > 0L) as.vector(xi_beta_imp[[cf_g$group_id]] %*% lf$d_bar_plus_j) else 0.0) else j_T_adj
    psi_T <- direct_T + a_T + j_T

    # --- Treated OR baseline B_gl (if shared_or fitted) ---
    direct_Bg <- a_Bg <- or_Bg_dir <- or_Bg_a <- or_Bg_j <- psi_Bg <- psi_numBg <- rep(NA_real_, B)
    if (!is.null(or_if)) {
      sh <- core_fit$shared_or
      bid <- lf$or_block_id
      bf <- sh$blocks[[bid]]
      p_or <- sh$p
      m_g_vec <- numeric(B)
      if (length(ix_g)) {
        Z_or_g <- if (p_or > 0L) X_chunk[ix_g, sh$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix_g), 0L)
        m_g_vec[ix_g] <- if (p_or > 0L) {
          (if (recipe$intercept) bf$alpha else 0.0) + as.vector((lf$L * Z_or_g) %*% sh$beta)
        } else {
          rep(if (recipe$intercept) bf$alpha else 0.0, length(ix_g))
        }
      }
      direct_Bg <- w_plus * (m_g_vec - lf$B_gl) / lf$d_plus
      a_Bg <- lf$a_coeff_Bg_g * xi_a[, ci_g]
      or_Bg_dir <- (if (recipe$intercept) or_if$direct$alpha[, match(bid, alpha_bids)] else 0.0) +
        (if (p_or > 0L) as.vector(or_if$direct$beta %*% lf$d_bar_plus_or) else 0.0)
      or_Bg_a <- (if (recipe$intercept) or_if$a$alpha[, match(bid, alpha_bids)] else 0.0) +
        (if (p_or > 0L) as.vector(or_if$a$beta %*% lf$d_bar_plus_or) else 0.0)
      or_Bg_j <- (if (recipe$intercept) or_if$j$alpha[, match(bid, alpha_bids)] else 0.0) +
        (if (p_or > 0L) as.vector(or_if$j$beta %*% lf$d_bar_plus_or) else 0.0)
      psi_Bg <- direct_Bg + a_Bg + or_Bg_dir + or_Bg_a + or_Bg_j
      psi_numBg <- lf$d_plus * psi_Bg + lf$B_gl * psi_denplus
    }

    # --- Control IPW arm C_ol and Control OR baseline B_ol (if trt_if fitted) ---
    direct_C <- a_C <- j_C <- j_C_adj <- trt_C_dir <- trt_C_a <- psi_C <- psi_denminus <- psi_numC <- rep(NA_real_, B)
    direct_Bo <- a_Bo <- trt_Bo_dir <- trt_Bo_a <- or_Bo_dir <- or_Bo_a <- or_Bo_j <- psi_Bo <- psi_numBo <- rep(NA_real_, B)

    if (!is.null(trt_if)) {
      tf <- core_fit$treatments[[lf$trt_id]]
      tif <- trt_if[[lf$trt_id]]
      trt_cols_F <- trt_recipe$index
      n_ctrl_coh <- length(lf$control_cohort_indices)

      w_minus <- numeric(B)
      j_k_vec <- numeric(B)
      m_k_vec <- numeric(B)
      a_denminus <- numeric(B)
      a_numC <- numeric(B)
      a_C_vec <- numeric(B)
      j_C_vec <- numeric(B)
      a_Bo_vec <- numeric(B)

      sh <- core_fit$shared_or
      p_or <- if (!is.null(sh)) sh$p else 0L
      bid <- lf$or_block_id
      bf <- if (!is.null(sh)) sh$blocks[[bid]] else NULL

      # Precompute control cohort sums for psi_denminus and psi_numC
      for (idx_k in seq_len(n_ctrl_coh)) {
        ci_k <- lf$control_cohort_indices[idx_k]
        cid_k <- lf$control_cell_ids[idx_k]
        cf_k <- core_fit$cells[[cid_k]]
        gf_k <- core_fit$groups[[cf_k$group_id]]
        a_k <- st$a[ci_k]

        ix_k <- which(ci == ci_k)
        if (length(ix_k)) {
          X_trt_k <- X_chunk[ix_k, trt_cols_F, drop = FALSE]
          lv_k <- ipw_logit_values(X_trt_k, tf$coefficients)
          wm_k <- (w[ix_k] * a_k) * lv_k$odds
          w_minus[ix_k] <- wm_k

          Z_imp_k <- if (gf_k$p > 0L) X_chunk[ix_k, gf_k$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix_k), 0L)
          j_k_vec[ix_k] <- if (gf_k$p > 0L) {
            cf_k$alpha + as.vector((lf$L * Z_imp_k) %*% gf_k$beta)
          } else {
            rep(cf_k$alpha, length(ix_k))
          }

          if (!is.null(sh)) {
            Z_or_k <- if (p_or > 0L) X_chunk[ix_k, sh$kept_cols_F, drop = FALSE] else matrix(0.0, length(ix_k), 0L)
            m_k_vec[ix_k] <- if (p_or > 0L) {
              (if (recipe$intercept) bf$alpha else 0.0) + as.vector((lf$L * Z_or_k) %*% sh$beta)
            } else {
              rep(if (recipe$intercept) bf$alpha else 0.0, length(ix_k))
            }
          }
        }

        xi_ak <- xi_a[, ci_k]
        a_C_vec <- a_C_vec + lf$a_coeff_C_k[idx_k] * xi_ak
        if (details) j_C_vec <- j_C_vec +
          (if (recipe$intercept) lf$m_minus_k[idx_k] * alpha_if(cid_k) else 0.0) +
          (if (gf_k$p > 0L) as.vector(xi_beta_imp[[cf_k$group_id]] %*% lf$d_bar_minus_j_k[[idx_k]]) else 0.0)
        if (!is.null(sh)) {
          a_Bo_vec <- a_Bo_vec + lf$a_coeff_Bo_k[idx_k] * xi_ak
        }
        # For denminus: (H_N[w I_k o_l]) * xi_{a_k} = (d_minus * m_minus_k / a_k) * xi_{a_k}
        a_denminus <- a_denminus + ((lf$d_minus * lf$m_minus_k[idx_k]) / a_k) * xi_ak
      }

      direct_C <- w_minus * (j_k_vec - lf$C_ol) / lf$d_minus
      a_C <- a_C_vec
      j_C_adj <- eval_adjoint_chunk(lf$adjoint_C)
      j_C <- if (details) j_C_vec else j_C_adj
      trt_C_dir <- as.vector(tif$direct %*% lf$g_C_gamma)
      trt_C_a <- as.vector(tif$a %*% lf$g_C_gamma)
      psi_C <- direct_C + a_C + j_C + trt_C_dir + trt_C_a

      # Denminus and numC IFs
      # Note: H_N[wdot_- X^trt] is recovered from the fitted control slice or from
      # psi_numC = d_minus * psi_C + C_ol * psi_denminus
      if (details) {
        mean_wdot_X <- lf$denminus_gamma
        psi_denminus <- w_minus + a_denminus + as.vector(tif$gamma %*% mean_wdot_X)
        psi_numC <- lf$d_minus * psi_C + lf$C_ol * psi_denminus
      }

      if (!is.null(or_if)) {
        direct_Bo <- w_minus * (m_k_vec - lf$B_ol) / lf$d_minus
        a_Bo <- a_Bo_vec
        trt_Bo_dir <- as.vector(tif$direct %*% lf$g_Bo_gamma)
        trt_Bo_a <- as.vector(tif$a %*% lf$g_Bo_gamma)
        or_Bo_dir <- (if (recipe$intercept) or_if$direct$alpha[, match(bid, alpha_bids)] else 0.0) +
          (if (p_or > 0L) as.vector(or_if$direct$beta %*% lf$d_bar_minus_or) else 0.0)
        or_Bo_a <- (if (recipe$intercept) or_if$a$alpha[, match(bid, alpha_bids)] else 0.0) +
          (if (p_or > 0L) as.vector(or_if$a$beta %*% lf$d_bar_minus_or) else 0.0)
        or_Bo_j <- (if (recipe$intercept) or_if$j$alpha[, match(bid, alpha_bids)] else 0.0) +
          (if (p_or > 0L) as.vector(or_if$j$beta %*% lf$d_bar_minus_or) else 0.0)
        psi_Bo <- direct_Bo + a_Bo + trt_Bo_dir + trt_Bo_a + or_Bo_dir + or_Bo_a + or_Bo_j
        if (details) {
          psi_numBo <- lf$d_minus * psi_Bo + lf$B_ol * psi_denminus
        }
      }
    }

    # Assemble method IFs and 8-channel decompositions
    if (!is.null(trt_if)) {
      comp_ipw <- cbind(
        direct = direct_T - direct_C,
        a_direct = a_T - a_C,
        j_direct = j_T - j_C,
        treatment_direct = -trt_C_dir,
        treatment_a = -trt_C_a,
        OR_direct = numeric(B),
        OR_a = numeric(B),
        OR_j = numeric(B)
      )
      psi_ipw_mat[, j] <- rowSums(comp_ipw)
    } else {
      comp_ipw <- NULL
    }

    if (!is.null(or_if)) {
      comp_reg <- cbind(
        direct = direct_T - direct_Bg,
        a_direct = a_T - a_Bg,
        j_direct = j_T,
        treatment_direct = numeric(B),
        treatment_a = numeric(B),
        OR_direct = -or_Bg_dir,
        OR_a = -or_Bg_a,
        OR_j = -or_Bg_j
      )
      psi_reg_mat[, j] <- rowSums(comp_reg)
    } else {
      comp_reg <- NULL
    }

    if (!is.null(trt_if) && !is.null(or_if)) {
      comp_dr <- cbind(
        direct = (direct_T - direct_Bg) - (direct_C - direct_Bo),
        a_direct = (a_T - a_Bg) - (a_C - a_Bo),
        j_direct = j_T - j_C,
        treatment_direct = -(trt_C_dir - trt_Bo_dir),
        treatment_a = -(trt_C_a - trt_Bo_a),
        OR_direct = -or_Bg_dir + or_Bo_dir,
        OR_a = -or_Bg_a + or_Bo_a,
        OR_j = -or_Bg_j + or_Bo_j
      )
      psi_dr_mat[, j] <- rowSums(comp_dr)
    } else {
      comp_dr <- NULL
    }

    phi[, j] <- switch(
      est_method,
      ipw = psi_ipw_mat[, j],
      reg = psi_reg_mat[, j],
      dr = psi_dr_mat[, j]
    )

    if (details) {
      link_details[[j]] <- list(
        link_id = lf$link_id,
        denplus = psi_denplus,
        numT = psi_numT,
        T = psi_T,
        numBg = psi_numBg,
        Bg = psi_Bg,
        denminus = psi_denminus,
        numC = psi_numC,
        C = psi_C,
        numBo = psi_numBo,
        Bo = psi_Bo,
        ipw = psi_ipw_mat[, j],
        reg = psi_reg_mat[, j],
        dr = psi_dr_mat[, j],
        components_ipw = comp_ipw,
        components_reg = comp_reg,
        components_dr = comp_dr,
        j_T_forward = j_T,
        j_T_adjoint = j_T_adj,
        j_C_forward = j_C,
        j_C_adjoint = j_C_adj
      )
    }
  }

  cp04_finite(phi, "complete pure-imputation link IF")
  phi_out <- if (identical(scale, "core")) st$lambda * phi else phi

  det_out <- NULL
  if (details) {
    xi_alpha_imp <- matrix(0.0, B, n_cells, dimnames = list(NULL, vapply(core_fit$cells, `[[`, character(1L), "name")))
    cell_scores <- matrix(0.0, B, n_cells)
    for (cid in sets$cids) {
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
      or_if = or_if,
      trt_if = trt_if,
      links = link_details,
      psi_ipw = psi_ipw_mat,
      psi_reg = psi_reg_mat,
      psi_dr = psi_dr_mat
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

# ------------------------------------------------------------------------------
# 5. Sealed-Owner Pure-Imputation Backend (`fit_imputation_backend`)
# ------------------------------------------------------------------------------

validate_imputation_backend <- function(backend) {
  if (!inherits(backend, "cdid_imputation_backend") ||
      !ipw_owner_seal_ok(backend, "imputation_backend")) {
    imputation_abort("Expected the immutable sealed pure-imputation backend (`cdid_imputation_backend`).", "integrity")
  }
  reader <- if (is.environment(backend$reader_ref) &&
                exists("cp04_reader", backend$reader_ref, inherits = FALSE)) {
    get("cp04_reader", backend$reader_ref, inherits = FALSE)
  } else {
    NULL
  }
  validate_cp04_reader(reader)
  if (!is.null(backend$design_ref)) {
    d_owner <- if (is.environment(backend$design_ref) &&
                   exists("imputation_design_owner", backend$design_ref, inherits = FALSE)) {
      get("imputation_design_owner", backend$design_ref, inherits = FALSE)
    } else {
      NULL
    }
    validate_imputation_design_owner(d_owner, reader = reader)
  }
  if (!is.null(backend$reference_ref)) chained_backend_reference(backend)
  invisible(backend)
}

# Logical payload, not RSS. These are the six fields populated by the core fit.
# Use the unreduced p as an upper bound for every retained/aliased slope subset.
# List/vector headers, attributes, names and string bytes follow cp04_payload;
# matrix charges use dimensions only, without allocating p-by-p placeholders.
imputation_retained_core_bound <- function(cat, recipe, trt_recipe) {
  E <- function(x) unname(cp04_payload(x)["entries"])
  # A names attribute costs 9 + E(names). A matrix costs its entries + 8,
  # the one-slot attributes list (9), and its two-element dim vector (10).
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
  has_or <- length(cat$or_blocks) > 0L
  has_trt <- length(cat$treatment_specs) > 0L

  # A pooled summary and fit each own block records and accounting. The fit
  # additionally owns its slope_fit (including a second R and coefficients).
  wls <- function(bn) {
    block <- R(c(r11 = A(1), r1d = A(p, z), r1y = A(1), n = A(1),
      weight_sum = A(1), mean_d = A(p, z), mean_y = A(1),
      augmented_rows = A(1), residual_rows = A(1)))
    blocks <- R(rep(block, length(bn)), nm = bn)
    accounting <- R(c(n = A(1), p = A(1), n_blocks = A(1), residual_dim = A(2),
      max_augmented_rows = A(1), max_residual_rows = A(1),
      max_chunk_merge_rows = A(1), max_residual_merge_rows = A(1), n_compressions = A(1)))
    simple_summary <- R(c(factor = M(p + 1, p + 1), n = A(1), columns = E(z)),
      "cdid_or_wls_summary")
    simple_fit <- R(c(coefficients = A(p, z), rss = A(1), R = M(p, p),
      pivot = A(p), n = A(1), rank = A(1), rank_tol = A(1), cutoff = A(1)), "cdid_or_wls_fit")
    if (recipe$intercept) {
      summary <- R(c(residual = simple_summary, blocks = blocks, n = A(1), p = A(1),
        columns = E(z), accounting = accounting), "cdid_or_pooled_wls_summary")
      fit <- R(c(coefficients = A(p, z), intercepts = A(length(bn), bn), rss = A(1),
        R = M(p, p), pivot = A(p), n = A(1), p = A(1), n_blocks = A(1),
        rank = A(1), rank_tol = A(1), cutoff = A(1), blocks = blocks,
        slope_fit = if (p > 0) simple_fit else 0, accounting = accounting), "cdid_or_pooled_wls_fit")
    } else {
      summary <- simple_summary
      fit <- simple_fit
    }
    basis <- R(c(slope_columns = E(z), kept_columns = E(z), kept_index = A(p),
      proposed_columns = E(z), aliased_columns = E(z), nontrivial_candidates_retained = E(z),
      alias_map = M(p, p, z, z, TRUE), block_constants = M(length(bn), p, bn, z, TRUE),
      reconstruction_control_scaled = A(p), reconstruction_treated_scaled = A(p),
      reconstruction_roundoff_threshold = A(1), rows_checked = A(1),
      proposal = E("R stats::lm.fit / dqrdc2 pivoted QR at rank_tol on the pooled (partialled) factor"),
      rank_check = E("unchanged relative-diagonal rule of the accepted kernel on the reduced factor"),
      frozen_during_differentiation = A(1), numerical_certificate = A(1)))
    c(fit = fit, summary = summary, basis = basis, rss = A(1), objective = A(1))
  }

  cells <- V(vapply(cat$cell_specs, function(cs) E(cs) + R(c(
    group_id = A(1), local_id = A(1), n_pair = A(1), weight_sum = A(1), mass = A(1),
    alpha = A(1), mean_d = A(p, z), mean_y = A(1), a_k = A(1), M_one = A(1),
    Y_one = A(1), M_z = A(p, z), Y_z = A(p, z), M_d_kept = A(p, z))), numeric(1)))
  groups <- V(vapply(cat$group_specs, function(gs) {
    bn <- vapply(cat$cell_specs[gs$cell_ids], `[[`, character(1), "name")
    E(gs) + R(c(keep = A(p), p = A(1), kept_columns = E(z), kept_cols_F = A(p),
      beta = A(p, z), alpha = A(length(bn), bn), wls(bn)))
  }, numeric(1)))
  shared <- if (has_or) {
    bn <- vapply(cat$or_blocks, `[[`, character(1), "name")
    blocks <- V(vapply(cat$or_blocks, function(b) {
      k <- length(b$control_cell_ids)
      E(b) + R(c(alpha = A(1), mass = A(1), weight_sum = A(1), mean_d = A(p, z),
        a_force_alpha = A(k), a_force_beta = M(p, k), j_force_alpha_alpha = A(k),
        j_force_alpha_beta = k + 8 + k * A(p), j_force_beta_alpha = M(p, k),
        j_force_beta_beta = k + 8 + k * M(p, p, z, z, TRUE)))
    }, numeric(1)))
    nt <- length(bn) + p
    R(c(keep = A(p), p = A(1), kept_columns = E(z), kept_cols_F = A(p),
      intercept = A(1), n_blocks = A(1), alpha = A(length(bn), bn), beta = A(p, z),
      theta = A(nt, paste0("OR.theta", seq_len(nt) - 1L)), blocks = blocks, wls(bn)))
  } else 0
  treatments <- if (has_trt) V(vapply(cat$treatment_specs, function(ts) {
    E(ts) + R(c(coefficients = A(q, paste0("gamma", seq_len(q) - 1L)),
      column_names = E(trt_recipe$columns),
      H_chol = M(q, q, trt_recipe$columns, trt_recipe$columns, TRUE),
      # The current lfit$H partially matches H_chol and retains that alias too.
      H = M(q, q, trt_recipe$columns, trt_recipe$columns, TRUE),
      n_risk_F = A(1), a_forcing = M(q, length(ts$all_cohort_indices)),
      clipping = E(ipw_clip_summary(prob = 0.5, controls = FALSE, weights = 1))))
  }, numeric(1))) else 0
  # Every adjoint retains one cell vector and one slope vector per group,
  # including the headers of zero-length vectors in intercept-only fits.
  adjoint <- R(c(tilde_lambda = A(nc), tilde_nu = ng + 8 + ng * A(p)))
  links <- V(vapply(cat$links, function(lk) {
    k <- length(lk$control_cell_ids)
    E(lk) + R(c(mu = A(1), mu_ipw = A(1), mu_or = A(1), mu_reg = A(1), mu_dr = A(1),
      d_plus = A(1), num_T = A(1), T_gl = A(1), num_Bg = A(1), B_gl = A(1),
      d_minus = A(1), num_C = A(1), C_ol = A(1), num_Bo = A(1), B_ol = A(1),
      d_bar_plus_j = A(p, z), d_bar_plus_or = A(p, z), d_bar_minus_or = A(p, z),
      d_bar_diff_or = A(p, z), a_coeff_T_g = A(1), a_coeff_Bg_g = A(1),
      a_coeff_C_k = A(k), a_coeff_Bo_k = A(k), m_minus_k = A(k),
      d_bar_minus_j_k = k + 8 + k * A(p), denminus_gamma = A(q),
      g_C_gamma = A(q, trt_recipe$columns), g_Bo_gamma = A(q, trt_recipe$columns),
      g_dr_gamma = A(q, trt_recipe$columns), adjoint_T = adjoint,
      if (has_trt) c(adjoint_C = adjoint)))
  }, numeric(1)))
  # cp04_payload charges logical aliases separately, even if R shares storage.
  c(cells = cells, groups = groups, shared_or = shared,
    shared = if (has_or) shared else E(list(p = 0L, keep = integer(), beta = numeric())),
    treatments = treatments, links = links)
}

fit_imputation_backend <- function(reader, link_views, control = ipw_owner_control(),
                                   glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
                                   rank_tol = or_public_rank_tol(),
                                   outcome_recipe,
                                   treatment_columns = NULL,
                                   est_method = "dr",
                                   imputation_pooling = "separate",
                                   design_owner = NULL,
                                   X_frame = NULL,
                                   retained_external_entries = 0,
                                   external_work_ops = 0,
                                   automatic_chunk_rows = FALSE,
                                   selection_model = NULL,
                                   selection_pooling = NULL,
                                   selection_xformla = NULL,
                                   dr_method = NULL,
                                   pair_outcomes = NULL,
                                   reference = NULL) {
  # Note: selection_model, selection_pooling, and selection_xformla are NEVER
  # evaluated or forced under pure imputation.
  sel_rec <- imputation_inactive_selection_record()

  validate_cp04_reader(reader)
  est_method <- imputation_method_choice(est_method = est_method, dr_method = dr_method)
  imputation_pooling <- imputation_pooling_choice(imputation_pooling)
  if (missing(rank_tol)) imputation_abort("Supply explicit private `rank_tol`.")
  rank_tol <- or_baseline_rank_tol(rank_tol)
  ctrl <- validate_ipw_owner_control(control)
  gc <- if (est_method %in% c("ipw", "dr")) {
    validate_ipw_glm_control(glm_control)
  } else {
    list(epsilon = 1e-12, maxit = 50L, trace = FALSE)
  }
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE)) {
    imputation_abort("Invalid external pipeline reservation.")
  }
  if (!identical(automatic_chunk_rows, isTRUE(automatic_chunk_rows))) {
    imputation_abort("Invalid private automatic row-planning flag.")
  }

  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (missing(outcome_recipe)) {
    imputation_abort("Supply the explicit private outcome recipe (intercept flag and covariate roles).")
  }

  # D60: the accepted core runs on the recipe-cohort layout (each cohort used
  # under recipe r is one layout cohort with W_g, a_k^{g,F} and z(C^r)); its
  # rows are summed back to one influence row per original unit.
  lay <- NULL
  if (!is.null(reference)) {
    if (!is.null(design_owner) || !is.null(X_frame) || !is.null(pair_outcomes))
      imputation_abort("A reference-input imputation reads its design only through the sealed reference reader.", "identity")
    lay <- chained_reference_layout(reader, reference, link_views, ctrl,
      retained = reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
        retained_external_entries + reference$accounting$owned_entries, work = external_work_ops)
    st <- lay$state
    des <- list(X = lay$X, columns = reference$x_columns, source = "reference_recipe_layout", design_ref = NULL,
      owned_entries = 0, borrowed_entries = reference$accounting$owned_entries, work_ops = 0)
  } else {
  # Resolve full-F outcome/treatment design matrix X_F (never confusing probability Q with X_F)
  des <- imputation_resolve_full_design(
    reader = reader,
    frame_owner = frame,
    design_owner = design_owner,
    X_frame = X_frame,
    control = ctrl
  )
  }
  X_F <- des$X
  f_cols <- des$columns

  or_chk <- dr_check_outcome_recipe(outcome_recipe, f_cols, "baseline")
  recipe <- or_chk$recipe
  or_frame_index <- or_chk$frame_index
  trt_recipe <- if (est_method %in% c("ipw", "dr")) {
    dr_check_treatment_recipe(treatment_columns, f_cols, "baseline")
  } else {
    structure(list(columns = character(0L), index = integer(0L), p = 0L), class = "cdid_dr_treatment_recipe")
  }

  if (!is.list(link_views) || !length(link_views)) {
    imputation_abort("Supply actual private link views.")
  }
  L <- length(link_views)
  G <- length(st$cohorts)
  p_or <- recipe$p
  p_trt <- trt_recipe$p

  base_entries <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    des$owned_entries + des$borrowed_entries + retained_external_entries
  header <- 4096 + 1024 * L * (G + 1)
  cp04_budget(
    ctrl, c(base_entries, header), header,
    c(external_work_ops, des$work_ops, 128 * L * L * (G + 1)),
    header, "pure-imputation metadata"
  )

  records <- vector("list", L)
  for (ell in seq_len(L)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:")) {
      imputation_abort("Pure-imputation backend requires the supplied comparison-catalog functional key.")
    }
    key_bytes <- nchar(source, type = "bytes")
    cp04_budget(
      ctrl, c(base_entries, header), header + 16 * key_bytes,
      c(external_work_ops, des$work_ops, 128 * L * L * (G + 1) + 64 * key_bytes),
      header, "pure-imputation link functional key"
    )
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate)) {
      imputation_abort("Legacy nuisance candidates are unsupported.")
    }
    if (!is.null(v$nuisance_ref)) {
      nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair)) {
        imputation_abort("Cross-pair or unverified nuisance sharing is unsupported.")
      }
    }
    if (is.null(lay)) {
      g_pos <- match(v$cohorts$g, st$cohorts)
      k_pos <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    } else {
      g_pos <- lay$view_cohorts[[ell]][1L]; k_pos <- lay$view_cohorts[[ell]][-1L]
    }
    if (is.na(g_pos) || anyNA(k_pos) || !isTRUE(all(v$support$by_cohort_pair_count >= 1))) {
      imputation_abort("Positive actual-pair cohort support required.", "support")
    }
    or_pooled_transitions(v$pair, "baseline")
    lid <- if (!is.null(v$link_id) && nzchar(as.character(v$link_id))) as.character(v$link_id) else paste0("link_", ell)
    records[[ell]] <- list(
      pair = v$pair,
      target = as.integer(g_pos),
      controls = as.integer(k_pos),
      link_id = lid,
      source_functional_key = source,
      functional_key = paste0("cmpcat:IMP-", toupper(est_method), ":", imputation_pooling, ":", source,
        if (is.null(lay)) "" else paste0("|reference:", reference$recipes$key[reference$view_recipe[ell]]))
    )
  }
  lids <- vapply(records, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(lids)) || anyDuplicated(lids)) {
    imputation_abort("Distinct nonempty link IDs required.")
  }

  cat <- imputation_build_catalog(records, st$cohorts, est_method, imputation_pooling,
    families = if (is.null(lay)) NULL else lay$families)
  n_cells <- length(cat$cell_specs)
  n_groups <- length(cat$group_specs)
  n_blocks <- length(cat$or_blocks)
  n_trt <- length(cat$treatment_specs)

  # Exact private fixture differences are retained for the same fit/IF law.
  # Charge their payload before scanning or retaining them.
  pair_entries <- if (is.null(pair_outcomes)) 0 else unname(cp04_payload(pair_outcomes)["entries"])
  base_entries <- base_entries + pair_entries
  if (!is.null(pair_outcomes) && (!is.list(pair_outcomes) || is.object(pair_outcomes))) {
    imputation_abort("`pair_outcomes` must be a named private list of full-F differences.")
  }

  # Measure already-built metadata, then charge every future compact field.
  # The placeholder accounting vectors have their maximum retained lengths.
  obj <- structure(
    list(
      selection_est_method = "imputation",
      est_method = est_method,
      dr_method = est_method,
      imputation_pooling = imputation_pooling,
      uses_selection_models = FALSE,
      selection = "none",
      selection_model = "none",
      selection_pooling = "none",
      selection_record = sel_rec,
      universe = "F",
      law = paste0(
        "IMPUTATION:", est_method, ":", imputation_pooling, ":",
        if (recipe$intercept) "local-intercepts" else "no-intercept",
        if (is.null(lay)) ":actual-pair:fixed-w:full-F" else ":actual-pair:W_g:full-F:recipe-cohort-layout"
      ),
      outcome_recipe = recipe,
      or_frame_index = as.integer(or_frame_index),
      treatment_recipe = trt_recipe,
      pair_outcomes = pair_outcomes,
      design_source = des$source,
      design_columns = f_cols,
      X_F = if (is.null(lay)) NULL else X_F, # A supplied X_frame is retained by its sealed design owner.
      layout = if (is.null(lay)) NULL else lay[c("state", "orig_row", "by_orig", "cohort", "recipe", "families",
        "N", "N_layout")],
      pairs = cat$pairs,
      cell_specs = cat$cell_specs,
      cell_lookup = cat$cell_lookup,
      active_cohorts = cat$active_cohorts,
      group_specs = cat$group_specs,
      or_blocks = cat$or_blocks,
      treatment_specs = cat$treatment_specs,
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
        automatic_chunk_rows = automatic_chunk_rows,
        n_cells = n_cells,
        n_groups = n_groups,
        n_imputation_fits = n_groups,
        n_imputation_local_intercepts = if (recipe$intercept) n_cells else 0L,
        n_imputation_slopes_by_group = rep.int(as.integer(p_or), n_groups),
        n_or_fits = if (n_blocks > 0L) 1L else 0L,
        n_or_blocks = n_blocks,
        n_or_local_intercepts = if (n_blocks > 0L && recipe$intercept) n_blocks else 0L,
        n_or_shared_slopes = if (n_blocks > 0L) as.integer(p_or) else 0L,
        n_treatment_models = n_trt,
        n_treatment_fits = n_trt,
        n_treatment_columns = p_trt,
        n_selection_fits = 0L,
        n_selection_cells = 0L,
        n_selection_blocks = 0L,
        n_selection_parameters = 0L,
        observation_logit_models = 0L,
        full_F_masks_per_model = 0L,
        full_F_prediction_matrices = 0L,
        full_F_by_link_IF_retained = FALSE,
        inverse_matrices_retained = FALSE,
        dense_dummy_regression = FALSE,
        peak_rss_measured = FALSE
      )
    ),
    class = c("cdid_imputation_backend", "cdid_chained_imputation")
  )
  if (!is.null(lay)) {
    obj$reference_ref <- reference$ref
    obj$accounting$n_reference_recipes <- nrow(reference$recipes)
    obj$accounting$n_layout_rows <- lay$N_layout
  }
  obj["ref"] <- list(NULL)
  # Scans run over the layout rows on the D60 recipe-cohort layout.
  N_scan <- if (is.null(lay)) st$N else lay$N_layout
  B_chunk <- min(N_scan, ctrl$max_chunk_rows)
  retained_bound <- unname(cp04_payload(obj)["entries"]) - pair_entries +
    sum(imputation_retained_core_bound(cat, recipe, trt_recipe))
  fit_scratch <- 64 * B_chunk * (p_or + p_trt + G + 8) +
    64 * N_scan * (if (n_trt > 0L) (p_trt + 4) else 2) +
    128 * ((p_or + 2)^2 * max(1L, n_cells) + (p_trt + 2)^2) + 8192
  log_n <- 1 + ceiling(log2(N_scan + 1))
  fit_work <- 256 * N_scan * (
    (n_cells + n_blocks + 1) * (p_or + 2)^2 * log_n +
      n_trt * (p_trt + 2)^2 * gc$maxit +
      L * (p_or + p_trt + G + 4)
  ) + 128 * retained_bound

  obj$accounting$retained_bound <- base_entries + retained_bound
  obj$accounting$fit_scratch_bound <- fit_scratch
  obj$accounting$work_ops_bound <- external_work_ops + des$work_ops + fit_work

  cp04_budget(
    ctrl, c(base_entries, retained_bound), fit_scratch,
    c(external_work_ops, des$work_ops, fit_work),
    retained_bound + pair_entries, "pure-imputation fitting preflight"
  )

  pair_reader <- function(rows, pair, read_outcomes = TRUE) {
    if (!is.null(lay)) return(cp04_pair_rows(reader, lay$orig_row[rows], pair, read_outcomes))
    pr <- cp04_pair_rows(reader, rows, pair, read_outcomes)
    if (read_outcomes && !is.null(pair_outcomes)) {
      pk <- paste0(as.character(pair$s), ":", as.character(pair$t))
      dy_vec <- pair_outcomes[[pk]]
      if (is.null(dy_vec) || !is.numeric(dy_vec) || length(dy_vec) != st$N) {
        imputation_abort(paste0("Invalid supplied `pair_outcomes` entry for pair ", pk, "."))
      }
      dy_sub <- as.double(dy_vec[rows])
      obs_sub <- !is.na(dy_sub)
      if (!identical(obs_sub, as.logical(pr$S))) {
        imputation_abort(paste0("Supplied `pair_outcomes` mask disagrees with panel observed pair indicator for ", pk, "."), "support")
      }
      if (any(!is.finite(dy_sub[obs_sub]))) {
        imputation_abort(paste0("Non-finite observed outcome difference in `pair_outcomes` for ", pk, "."))
      }
      pr$dy <- numeric(length(rows))
      pr$dy[obs_sub] <- dy_sub[obs_sub]
    }
    pr
  }

  core_fit <- imputation_fit_core(
    st = st,
    X_F = X_F,
    pair_reader = pair_reader,
    cat = cat,
    recipe = recipe,
    or_frame_index = or_frame_index,
    trt_recipe = trt_recipe,
    est_method = est_method,
    imputation_pooling = imputation_pooling,
    glm_control = gc,
    rank_tol = rank_tol,
    ctrl = ctrl
  )

  # Fill only the preflighted compact fields; preserve NULL list slots.
  obj[c("cells", "groups", "shared_or", "shared", "treatments", "links")] <- list(
    core_fit$cells, core_fit$groups, core_fit$shared_or,
    if (!is.null(core_fit$shared_or)) core_fit$shared_or else list(p = 0L, keep = integer(0L), beta = numeric(0L)),
    core_fit$treatments, core_fit$links
  )
  obj$accounting$n_imputation_slopes_by_group <- vapply(core_fit$groups, `[[`, integer(1L), "p")
  obj$accounting$n_or_shared_slopes <- if (!is.null(core_fit$shared_or)) core_fit$shared_or$p else 0L
  payload <- cp04_payload(obj)
  actual <- unname(payload["entries"])
  if (actual > retained_bound + pair_entries) {
    imputation_abort("Retained pure-imputation backend bound underestimated.", "bound")
  }
  cp04_budget(
    ctrl, c(obj$accounting$borrowed_entries, actual), fit_scratch,
    obj$accounting$work_ops_bound, actual, "retained pure-imputation backend"
  )
  obj$accounting$owned_entries <- actual
  obj$accounting$character_bytes <- unname(payload["char_bytes"])
  ipw_owner_seal(obj, "imputation_backend")
}

fit_chained_imputation <- fit_imputation_backend

imputation_backend_X_F <- function(backend, reader) {
  if (!is.null(backend$design_ref)) {
    d_owner <- get("imputation_design_owner", backend$design_ref, inherits = FALSE)
    return(d_owner$X)
  }
  if (!is.null(backend$X_F)) return(backend$X_F)
  reader$state$Z
}

imputation_if_envelope <- function(backend, B, links = seq_along(backend$links), details = FALSE) {
  validate_imputation_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- if (is.null(backend$layout)) reader$state else backend$layout$state
  # A recipe-cohort layout evaluates up to (layout rows / original rows) rows
  # per original row; charge the largest such multiplicity.
  if (!is.null(backend$layout)) B <- B * max(1, backend$layout$by_orig$count)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  sets <- imputation_if_needed_sets(backend, links, backend$est_method, details)
  nl <- length(links)
  nc <- length(sets$cids)
  ng <- length(sets$gids)
  nb <- length(sets$bids)
  nt <- length(sets$tids)
  G <- length(st$cohorts)
  p_or <- backend$outcome_recipe$p
  p_trt <- backend$treatment_recipe$p

  # Parameter IF matrices exist only in explicit fixture diagnostics. Ordinary
  # replay uses scalar adjoints and keeps just the shared OR and link row widths.
  group_width <- if (details) sum(vapply(backend$groups, `[[`, integer(1L), "p")) else 0
  n_alpha <- if (details) nb else length(unique(vapply(backend$links[links], `[[`, integer(1L), "or_block_id")))
  if (is.null(backend$shared_or)) n_alpha <- 0L
  width <- 3 * G + group_width + 16 * n_alpha + 24 * p_or +
    nt * (6 * p_trt + 8) + 32 * nl +
    if (details) 4 * length(backend$cells) + 64 * nl + 32 else 8
  metadata <- 16384 + 1024 * (nc + length(backend$groups) + nb + nt + nl + G) +
    64 * length(backend$cells) * (p_or + 2) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1L), "link_id"), type = "bytes"))
  output <- B * (width + 16) + metadata
  scratch <- 64 * B * (width + length(backend$design_columns) + p_or + p_trt + 8) + 128 * ((p_or + 2)^2 + (p_trt + 2)^2) + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- 128 * B * (
    (nc + 1) * (p_or + 2)^2 * log_a +
      (ng + nb + nl + 1) * (p_or + 1) * (p_or + nc + 1) +
      nt * ((p_trt + 2)^2 + p_trt * (G + 4)) +
      G + nl * (nc + 1) * (p_or + p_trt + G + 8)
  ) + 64 * metadata
  list(
    cids = sets$cids,
    gids = sets$gids,
    bids = sets$bids,
    tids = sets$tids,
    output = output,
    scratch = scratch,
    work = work
  )
}

imputation_if_budget <- function(backend, B, links = seq_along(backend$links), details = FALSE) {
  pre <- imputation_if_envelope(backend, B, links, details)
  cp04_budget(
    backend$control,
    c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    pre$scratch,
    c(backend$accounting$work_ops_bound, pre$work),
    pre$output,
    "pure-imputation IF chunk"
  )
  pre
}

evaluate_imputation_if_chunk <- function(backend, rows, links = seq_along(backend$links),
                                         scale = "frame", details = FALSE) {
  validate_imputation_backend(backend)
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
    imputation_abort("Invalid explicit IF row/link batch, scale or details flag.")
  }

  if (!is.null(backend$layout)) return(evaluate_imputation_layout_if_chunk(backend, reader, rows, links, scale, details))
  pre <- imputation_if_budget(backend, length(rows), links, details)
  X_F <- imputation_backend_X_F(backend, reader)
  pair_reader <- function(r_sub, pair, read_outcomes = TRUE) {
    pr <- cp04_pair_rows(reader, r_sub, pair, read_outcomes)
    if (read_outcomes && !is.null(backend$pair_outcomes)) {
      pk <- paste0(as.character(pair$s), ":", as.character(pair$t))
      dy <- backend$pair_outcomes[[pk]][r_sub]
      observed <- pr$S == 1
      if (!identical(!is.na(dy), observed) || any(!is.finite(dy[observed]))) {
        imputation_abort("Retained private pair differences disagree with the sealed panel.", "support")
      }
      pr$dy <- numeric(length(r_sub))
      pr$dy[observed] <- dy[observed]
    }
    pr
  }

  res <- imputation_evaluate_if_chunk_core(
    st = st,
    X_F = X_F,
    pair_reader = pair_reader,
    core_fit = backend,
    recipe = backend$outcome_recipe,
    or_frame_index = backend$or_frame_index,
    trt_recipe = backend$treatment_recipe,
    est_method = backend$est_method,
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
      n_or_block_evaluations = length(pre$bids),
      n_treatment_score_evaluations = length(pre$tids),
      n_observation_logit_evaluations = 0L,
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  if (details) out$details <- res$details
  if (unname(cp04_payload(out)["entries"]) > pre$output) {
    imputation_abort("Pure-imputation IF output bound underestimated.", "bound")
  }
  out
}

# D60 recipe-cohort layout: the layout rows of these original rows are
# evaluated by the unchanged core (in bounded sub-chunks) and summed into one
# row per original unit, so every channel of a unit shares one influence row.
evaluate_imputation_layout_if_chunk <- function(backend, reader, rows, links, scale, details) {
  if (details) imputation_abort("Diagnostic channel details are not retained on the recipe-cohort layout.", "unsupported")
  lay <- backend$layout; st <- lay$state
  pre <- imputation_if_budget(backend, length(rows), links, FALSE)
  map <- chained_reference_layout_rows(lay, as.integer(rows))
  pair_reader <- function(r_sub, pair, read_outcomes = TRUE) cp04_pair_rows(reader, lay$orig_row[r_sub], pair, read_outcomes)
  B <- length(rows); BL <- backend$control$max_chunk_rows
  phi_l <- matrix(0, length(map$rows), length(links))
  if (length(map$rows)) for (s0 in seq.int(1L, length(map$rows), by = BL)) {
    ii <- seq.int(s0, min(length(map$rows), s0 + BL - 1L))
    res <- imputation_evaluate_if_chunk_core(st = st, X_F = backend$X_F, pair_reader = pair_reader,
      core_fit = backend, recipe = backend$outcome_recipe, or_frame_index = backend$or_frame_index,
      trt_recipe = backend$treatment_recipe, est_method = backend$est_method, rows = as.integer(map$rows[ii]),
      links = as.integer(links), scale = "frame", details = FALSE, sets = pre)
    phi_l[ii, ] <- res$phi
  }
  colnames(phi_l) <- vapply(backend$links[links], `[[`, character(1L), "link_id")
  phi <- chained_reference_layout_sum(phi_l, map$owner, B)
  cp04_finite(phi, "complete pure-imputation link IF (recipe-cohort layout)")
  ost <- reader$state
  if (identical(scale, "core")) phi <- ost$lambda * phi
  out <- list(rows = rows, id = ost$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (identical(scale, "core")) ost$n_A else ost$N,
    diagnostics = list(n_imputation_cell_evaluations = length(pre$cids), n_imputation_group_solves = length(pre$gids),
      n_or_block_evaluations = length(pre$bids), n_treatment_score_evaluations = length(pre$tids),
      n_observation_logit_evaluations = 0L, n_layout_rows = length(map$rows), output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (unname(cp04_payload(out)["entries"]) > pre$output) imputation_abort("Pure-imputation IF output bound underestimated.", "bound")
  out
}

replay_imputation_if <- function(backend, chunk_rows, link_batch_size, callback,
                                 scale = "frame", details = FALSE) {
  validate_imputation_backend(backend)
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  LB <- ipw_owner_int_bound(link_batch_size, "link_batch_size")
  if (B > backend$control$max_chunk_rows || !is.function(callback) ||
      !identical(details, isTRUE(details)) || !is.character(scale) ||
      length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core")) {
    imputation_abort("Invalid explicit replay controls/callback.")
  }
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N
  L <- length(backend$links)
  total <- backend$accounting$work_ops_bound
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    pre <- imputation_if_budget(backend, min(B, N), ll, details)
    total <- ipw_owner_safe_sum(
      c(total, ceiling(N / B) * pre$work),
      backend$control$max_work_ops,
      "total pure-imputation replay work"
    )
  }
  batches <- 0L
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1L))
      batch <- evaluate_imputation_if_chunk(backend, rr, ll, scale, details)
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

stream_imputation_covariance <- function(backend, links = seq_along(backend$links), retained = 0) {
  validate_imputation_backend(backend)
  cp04_stream_covariance(backend, links, retained = retained)
}

