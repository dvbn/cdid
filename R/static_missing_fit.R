# Private end-to-end static binary missing-outcome assembly.
#
# Wires the accepted private component APIs:
#   1. R/static_panel_data.R         (panel normalization, cell masks, lazy pair accessor, reference designs)
#   2. R/static_complete_core.R      (full-R_c treatment propensity logit and covariance helper)
#   3. R/static_observation_moments.R (route-specific normalized edge means and stayer covariances)
#   4. R/static_nuisance_models.R    (role-mask selection logit and actual-change imputation WLS)
#   5. R/static_observation_adapter.R (single-fit observation nuisance and moment adapter)
#   6. R/static_reconstructed_or.R   (deterministic LS stayer flow, pooled known-X OR, switch-step centering)
#   7. R/static_graph_gmm.R          (joint role-graph GMM solver and bounded original-unit IF streaming)
#
# Scope:
#   - Private assembly only. No public att_gt wrapper, no automatic target/graph
#     selection, and no silent pruning of unsupported cells or disconnected graphs.
#   - Treatment propensity logits are fitted on full R_c = 1{D_{t-1} = d} using
#     static_complete_fit_cell_logit (only when method %in% c("ipw", "dr")).
#   - Observation nuisances (q and/or j) are fitted once across all active role
#     populations via static_observation_adapter_fit using lazy pair callbacks.
#   - Reconstructed stayer OR is fitted across all declared cells in U >= A
#     using deterministic graph-LS flows varpi^or_c = D_c (D_c' D_c)^(-1) ell_c
#     (only when method %in% c("reg", "or", "dr")).
#   - All requested target cells in A and their role graphs are assembled into a
#     single solve_static_graph_gmm call with chunked original-unit IF replay.

static_missing_abort <- function(message, subclass = character(), diagnostics = list()) {
  stop(structure(
    list(message = message, call = NULL, diagnostics = diagnostics),
    class = unique(c(
      subclass,
      "cdid_static_missing_fit_error",
      "cdid_static_error",
      "cdid_error",
      "error",
      "condition"
    ))
  ))
}

static_missing_validate_method <- function(method) {
  if (missing(method) || !is.character(method) || !is.null(dim(method)) ||
      is.object(method) || length(method) != 1L || is.na(method) ||
      !(method %in% c("ipw", "reg", "or", "dr"))) {
    static_missing_abort(
      "Supply an explicit `method` in `c(\"ipw\", \"reg\", \"or\", \"dr\")`.",
      "cdid_static_missing_method_error"
    )
  }
  canonical <- if (identical(method, "or")) "reg" else method
  gmm_method <- if (identical(method, "reg")) "or" else method
  list(
    requested = method,
    canonical = canonical,
    gmm_method = gmm_method,
    use_logit = canonical %in% c("ipw", "dr"),
    use_or = canonical %in% c("reg", "dr")
  )
}

static_missing_validate_route <- function(route) {
  if (missing(route) || !is.character(route) || !is.null(dim(route)) ||
      is.object(route) || length(route) != 1L || is.na(route) ||
      !(route %in% c("baseline", "weighting", "imputation", "augmented"))) {
    static_missing_abort(
      "Supply an explicit observation `route` in `c(\"baseline\", \"weighting\", \"imputation\", \"augmented\")`.",
      "cdid_static_missing_route_error"
    )
  }
  list(
    route = route,
    uses_q = route %in% c("weighting", "augmented"),
    uses_j = route %in% c("imputation", "augmented")
  )
}

static_missing_resolve_panel <- function(panel,
                                         idname = NULL,
                                         tname = NULL,
                                         yname = NULL,
                                         dname = NULL,
                                         dates = NULL,
                                         unit_ids = NULL) {
  if (missing(panel) || is.null(panel)) {
    static_missing_abort(
      "`panel` must be supplied as a `cdid_static_panel_data` object (or a `data.frame` with explicit `idname`, `tname`, `yname`, `dname`).",
      "cdid_static_missing_input_error"
    )
  }
  if (inherits(panel, "cdid_static_panel_data")) {
    return(panel)
  }
  if (is.data.frame(panel)) {
    if (is.null(idname) || is.null(tname) || is.null(yname) || is.null(dname)) {
      static_missing_abort(
        "When `panel` is a `data.frame`, `idname`, `tname`, `yname`, and `dname` must all be supplied explicitly.",
        "cdid_static_missing_input_error"
      )
    }
    return(static_panel_data(
      data = panel,
      idname = idname,
      tname = tname,
      yname = yname,
      dname = dname,
      dates = dates,
      unit_ids = unit_ids
    ))
  }
  static_missing_abort(
    "`panel` must inherit from `cdid_static_panel_data` or be a `data.frame`.",
    "cdid_static_missing_input_error"
  )
}

static_missing_coerce_cell_map <- function(panel, spec, what) {
  uid_chr <- as.character(panel$unit_id)
  if (inherits(spec, "cdid_static_transition_cells")) {
    if (!identical(spec$n, panel$n) || !identical(as.character(spec$unit_id), uid_chr)) {
      static_missing_abort(
        sprintf("`%s` original-unit ordering does not match `panel`.", what),
        c("cdid_static_missing_alignment_error", "cdid_static_order_error")
      )
    }
    return(spec$cells)
  }
  if (inherits(spec, "cdid_static_transition_cell")) {
    if (!identical(spec$n, panel$n) || !identical(as.character(spec$unit_id), uid_chr)) {
      static_missing_abort(
        sprintf("`%s` original-unit ordering does not match `panel`.", what),
        c("cdid_static_missing_alignment_error", "cdid_static_order_error")
      )
    }
    out <- list(spec)
    names(out) <- spec$cell_id
    return(out)
  }
  if (is.list(spec) && !is.data.frame(spec) && length(spec) >= 1L &&
      all(vapply(spec, inherits, logical(1L), "cdid_static_transition_cell"))) {
    cids <- vapply(spec, `[[`, character(1L), "cell_id")
    if (anyDuplicated(cids)) {
      static_missing_abort(
        sprintf("`%s` contains duplicate `cell_id` values.", what),
        "cdid_static_duplicate_cell_error"
      )
    }
    for (cl in spec) {
      if (!identical(cl$n, panel$n) || !identical(as.character(cl$unit_id), uid_chr)) {
        static_missing_abort(
          sprintf("Cell `%s` in `%s` does not match `panel` original-unit ordering.", cl$cell_id, what),
          c("cdid_static_missing_alignment_error", "cdid_static_order_error")
        )
      }
    }
    names(spec) <- cids
    return(spec)
  }
  resolved <- static_transition_cells(panel = panel, cells = spec)
  resolved$cells
}

static_missing_resolve_cells_and_partition <- function(panel,
                                                       target_cells,
                                                       all_cells = NULL,
                                                       or_blocks = NULL,
                                                       use_or) {
  if (missing(target_cells) || is.null(target_cells)) {
    static_missing_abort(
      "Explicit `target_cells` must be supplied.",
      "cdid_static_missing_input_error"
    )
  }

  # Allow target_cells to be a character vector of cell_ids when all_cells is supplied
  if (is.character(target_cells) && is.null(dim(target_cells)) && !is.object(target_cells)) {
    if (is.null(all_cells)) {
      static_missing_abort(
        "When `target_cells` is a character vector of cell IDs, `all_cells` must supply the cell specifications.",
        "cdid_static_missing_input_error"
      )
    }
    if (length(target_cells) < 1L || anyNA(target_cells) || any(!nzchar(target_cells)) ||
        anyDuplicated(target_cells)) {
      static_missing_abort(
        "`target_cells` character vector must be non-empty, non-NA, and unique.",
        "cdid_static_missing_input_error"
      )
    }
    all_map_raw <- static_missing_coerce_cell_map(panel, all_cells, "all_cells")
    miss_t <- setdiff(target_cells, names(all_map_raw))
    if (length(miss_t) > 0L) {
      static_missing_abort(
        sprintf("`target_cells` not found in `all_cells`: %s.", paste(miss_t, collapse = ", ")),
        "cdid_static_missing_alignment_error"
      )
    }
    target_map <- all_map_raw[target_cells]
    extra_map <- all_map_raw[setdiff(names(all_map_raw), target_cells)]
  } else {
    target_map <- static_missing_coerce_cell_map(panel, target_cells, "target_cells")
    if (!is.null(all_cells)) {
      extra_raw <- static_missing_coerce_cell_map(panel, all_cells, "all_cells")
      # Verify that any shared cell_id in all_cells matches target_map
      shared_ids <- intersect(names(target_map), names(extra_raw))
      for (cid in shared_ids) {
        if (!identical(target_map[[cid]]$t_idx, extra_raw[[cid]]$t_idx) ||
            !identical(target_map[[cid]]$d, extra_raw[[cid]]$d)) {
          static_missing_abort(
            sprintf("Cell `%s` has conflicting `(t, d)` definitions between `target_cells` and `all_cells`.", cid),
            "cdid_static_missing_alignment_error"
          )
        }
      }
      extra_map <- extra_raw[setdiff(names(extra_raw), names(target_map))]
    } else {
      extra_map <- list()
    }
  }

  target_ids <- names(target_map)
  if (length(target_ids) < 1L) {
    static_missing_abort(
      "At least one target transition cell must be supplied.",
      "cdid_static_missing_input_error"
    )
  }

  # Validate support on every requested target cell in A (never prune empty cells)
  for (cid in target_ids) {
    cl <- target_map[[cid]]
    if (cl$n_risk < 1L || cl$n_switchers < 1L || cl$n_stayers < 1L ||
        !isTRUE(cl$support$is_supported)) {
      static_missing_abort(
        sprintf(
          "Requested target cell `%s` (t = %s, d = %d) has insufficient support (n_risk = %d, n_switchers = %d, n_stayers = %d); empty switcher or stayer roles are refused and never dropped.",
          cid, cl$t_label, cl$d, cl$n_risk, cl$n_switchers, cl$n_stayers
        ),
        c("cdid_static_support_error", "cdid_static_missing_support_error"),
        diagnostics = list(
          cell_id = cid,
          t = cl$t_num,
          d = cl$d,
          n_risk = cl$n_risk,
          n_switchers = cl$n_switchers,
          n_stayers = cl$n_stayers
        )
      )
    }
  }

  if (!use_or) {
    return(list(
      target_ids = target_ids,
      target_cells = target_map,
      all_ids = target_ids,
      all_cells = target_map,
      non_target_ids = character(0L),
      or_blocks = NULL,
      or_partition = NULL
    ))
  }

  # Outcome regression is active (method %in% c("reg", "or", "dr")): resolve U and P_or
  if (is.null(or_blocks)) {
    static_missing_abort(
      "Explicit `or_blocks` partition must be supplied when `method` is `\"reg\"`, `\"or\"`, or `\"dr\"`.",
      "cdid_static_missing_partition_error"
    )
  }

  combined_map <- c(target_map, extra_map)
  combined_ids <- names(combined_map)

  # Resolve keyword shorthands if explicitly supplied
  if (is.character(or_blocks) && length(or_blocks) == 1L && is.null(names(or_blocks)) &&
      or_blocks %in% c("separate", "pooled", "by_direction")) {
    if (identical(or_blocks, "separate")) {
      or_blocks_use <- stats::setNames(combined_ids, combined_ids)
    } else if (identical(or_blocks, "pooled")) {
      or_blocks_use <- stats::setNames(rep("pooled", length(combined_ids)), combined_ids)
    } else {
      dirs <- vapply(combined_map, function(cl) if (cl$d == 0L) "entry_block" else "exit_block", character(1L))
      or_blocks_use <- stats::setNames(dirs, combined_ids)
    }
    all_ids <- combined_ids
  } else if (is.list(or_blocks) && !is.object(or_blocks)) {
    mem_raw <- unlist(or_blocks, use.names = FALSE)
    if (is.character(mem_raw)) {
      miss_c <- setdiff(mem_raw, combined_ids)
      if (length(miss_c) > 0L) {
        static_missing_abort(
          sprintf("`or_blocks` references unknown cells not in `target_cells` or `all_cells`: %s.", paste(miss_c, collapse = ", ")),
          "cdid_static_missing_partition_error"
        )
      }
      miss_t <- setdiff(target_ids, mem_raw)
      if (length(miss_t) > 0L) {
        static_missing_abort(
          sprintf("`or_blocks` omits requested target cells: %s.", paste(miss_t, collapse = ", ")),
          "cdid_static_missing_partition_error"
        )
      }
      all_ids <- c(target_ids, setdiff(unique(mem_raw), target_ids))
    } else {
      all_ids <- combined_ids
    }
    or_blocks_use <- or_blocks
  } else if (is.character(or_blocks) && !is.null(names(or_blocks))) {
    nm_ob <- names(or_blocks)
    if (anyDuplicated(nm_ob)) {
      static_missing_abort(
        "`names(or_blocks)` must not contain duplicate cell IDs.",
        "cdid_static_missing_partition_error"
      )
    }
    miss_c <- setdiff(nm_ob, combined_ids)
    if (length(miss_c) > 0L) {
      static_missing_abort(
        sprintf("`names(or_blocks)` references unknown cells not in `target_cells` or `all_cells`: %s.", paste(miss_c, collapse = ", ")),
        "cdid_static_missing_partition_error"
      )
    }
    miss_t <- setdiff(target_ids, nm_ob)
    if (length(miss_t) > 0L) {
      static_missing_abort(
        sprintf("`or_blocks` omits requested target cells: %s.", paste(miss_t, collapse = ", ")),
        "cdid_static_missing_partition_error"
      )
    }
    all_ids <- c(target_ids, setdiff(nm_ob, target_ids))
    or_blocks_use <- or_blocks[all_ids]
  } else {
    all_ids <- combined_ids
    or_blocks_use <- or_blocks
  }

  all_map <- combined_map[all_ids]
  stub_cells <- lapply(all_map, function(cl) {
    list(p = 0L, z_names = character(0L), or_block = cl$or_block)
  })
  names(stub_cells) <- all_ids
  or_part <- static_reconstructed_or_validate_partition(stub_cells, or_blocks_use)

  # Only keep non-target cells in U \ A that actually belong to a block containing
  # at least one target cell, OR validate all declared cells in U.
  # Per Pipeline Step 4: "including declared non-target cells in pooled blocks".
  non_target_ids <- setdiff(all_ids, target_ids)
  for (cid in non_target_ids) {
    cl <- all_map[[cid]]
    if (cl$n_stayers < 1L) {
      static_missing_abort(
        sprintf(
          "Declared OR training cell `%s` (t = %s, d = %d) has zero stayers (`n_stayers = 0`).",
          cid, cl$t_label, cl$d
        ),
        c("cdid_static_support_error", "cdid_static_missing_support_error")
      )
    }
  }

  list(
    target_ids = target_ids,
    target_cells = target_map,
    all_ids = all_ids,
    all_cells = all_map,
    non_target_ids = non_target_ids,
    or_blocks = or_part$cell_block,
    or_partition = or_part
  )
}

static_missing_match_panel_dates <- function(panel, date_vec, what) {
  if (is.character(date_vec)) {
    idx <- match(date_vec, panel$date_labels)
  } else if (is.numeric(date_vec) && all(is.finite(date_vec))) {
    idx <- match(as.double(date_vec), panel$dates_num)
  } else {
    static_missing_abort(
      sprintf("`%s` must contain finite numeric or character study dates.", what),
      c("cdid_static_pair_endpoint_error", "cdid_static_missing_graph_error")
    )
  }
  if (anyNA(idx)) {
    bad <- date_vec[which(is.na(idx))[1L]]
    static_missing_abort(
      sprintf("`%s` contains date `%s` outside the panel date grid [%s].",
              what, as.character(bad), paste(panel$date_labels, collapse = ", ")),
      c("cdid_static_pair_endpoint_error", "cdid_static_missing_graph_error")
    )
  }
  as.integer(idx)
}

static_missing_coerce_one_role_block <- function(panel,
                                                 cell,
                                                 role,
                                                 spec,
                                                 default_block_id = paste0(cell$cell_id, "__", role)) {
  cid <- cell$cell_id
  base_t <- cell$t_prev_num
  targ_t <- cell$t_num

  if (inherits(spec, "cdid_static_role_block_spec")) {
    b_spec <- static_role_block_spec(
      block_id = spec$block_id,
      cell_id = cid,
      role = role,
      baseline_time = spec$baseline_time,
      target_time = spec$target_time,
      dates = spec$dates,
      edges = spec$edges,
      mhat = spec$mhat
    )
    if (b_spec$baseline_time != base_t || b_spec$target_time != targ_t) {
      static_missing_abort(
        sprintf(
          "Role block `%s` (cell `%s`, role `%s`) has `(baseline_time = %g, target_time = %g)` which disagrees with cell transition dates `(%g, %g)`.",
          b_spec$block_id, cid, role, b_spec$baseline_time, b_spec$target_time, base_t, targ_t
        ),
        "cdid_static_missing_alignment_error"
      )
    }
  } else if (is.data.frame(spec) || is.matrix(spec)) {
    edges_df <- as.data.frame(spec, stringsAsFactors = FALSE)
    if (!all(c("s", "u") %in% names(edges_df))) {
      if (ncol(edges_df) == 2L) {
        names(edges_df) <- c("s", "u")
      } else {
        static_missing_abort(
          sprintf("Cell `%s` role `%s` edge specification must have columns `s` and `u`.", cid, role),
          "cdid_static_missing_graph_error"
        )
      }
    }
    s_num <- static_missing_match_panel_dates(panel, edges_df$s, sprintf("role_graphs[[%s]]$%s$s", cid, role))
    u_num <- static_missing_match_panel_dates(panel, edges_df$u, sprintf("role_graphs[[%s]]$%s$u", cid, role))
    s_val <- panel$dates_num[s_num]
    u_val <- panel$dates_num[u_num]
    link_ids <- if ("link_id" %in% names(edges_df)) {
      as.character(edges_df$link_id)
    } else {
      sprintf("%s_%s_%s_%s", cid, role, panel$date_labels[s_num], panel$date_labels[u_num])
    }
    edges_clean <- data.frame(
      link_id = link_ids,
      s = s_val,
      u = u_val,
      stringsAsFactors = FALSE
    )
    dates_use <- sort(unique(c(edges_clean$s, edges_clean$u)))
    b_spec <- static_role_block_spec(
      block_id = default_block_id,
      cell_id = cid,
      role = role,
      baseline_time = base_t,
      target_time = targ_t,
      dates = dates_use,
      edges = edges_clean
    )
  } else if (is.list(spec) && !is.object(spec)) {
    edges_raw <- if (!is.null(spec$edges)) spec$edges else spec$pairs
    if (is.null(edges_raw)) {
      static_missing_abort(
        sprintf("Cell `%s` role `%s` graph specification must provide `edges`.", cid, role),
        "cdid_static_missing_graph_error"
      )
    }
    if (is.list(edges_raw) && !is.data.frame(edges_raw)) {
      pairs_mat <- do.call(rbind, lapply(edges_raw, function(p) {
        if (length(p) != 2L) {
          static_missing_abort(
            sprintf("Cell `%s` role `%s` edge pair must have length 2.", cid, role),
            "cdid_static_missing_graph_error"
          )
        }
        p
      }))
      edges_df <- data.frame(s = pairs_mat[, 1L], u = pairs_mat[, 2L], stringsAsFactors = FALSE)
    } else {
      edges_df <- as.data.frame(edges_raw, stringsAsFactors = FALSE)
      if (!all(c("s", "u") %in% names(edges_df)) && ncol(edges_df) == 2L) {
        names(edges_df) <- c("s", "u")
      }
    }
    if (!all(c("s", "u") %in% names(edges_df))) {
      static_missing_abort(
        sprintf("Cell `%s` role `%s` `edges` must provide columns `s` and `u`.", cid, role),
        "cdid_static_missing_graph_error"
      )
    }
    s_idx <- static_missing_match_panel_dates(panel, edges_df$s, sprintf("role_graphs[[%s]]$%s$edges$s", cid, role))
    u_idx <- static_missing_match_panel_dates(panel, edges_df$u, sprintf("role_graphs[[%s]]$%s$edges$u", cid, role))
    s_val <- panel$dates_num[s_idx]
    u_val <- panel$dates_num[u_idx]
    link_ids <- if ("link_id" %in% names(edges_df)) {
      as.character(edges_df$link_id)
    } else if (!is.null(spec$link_ids)) {
      as.character(spec$link_ids)
    } else {
      sprintf("%s_%s_%s_%s", cid, role, panel$date_labels[s_idx], panel$date_labels[u_idx])
    }
    edges_clean <- data.frame(
      link_id = link_ids,
      s = s_val,
      u = u_val,
      stringsAsFactors = FALSE
    )
    if (!is.null(spec$dates)) {
      d_idx <- static_missing_match_panel_dates(panel, spec$dates, sprintf("role_graphs[[%s]]$%s$dates", cid, role))
      dates_use <- panel$dates_num[d_idx]
    } else {
      dates_use <- sort(unique(c(edges_clean$s, edges_clean$u)))
    }
    bid <- if (!is.null(spec$block_id)) as.character(spec$block_id) else default_block_id
    base_use <- if (!is.null(spec$baseline_time)) {
      b_idx <- static_missing_match_panel_dates(panel, spec$baseline_time, sprintf("role_graphs[[%s]]$%s$baseline_time", cid, role))
      panel$dates_num[b_idx]
    } else {
      base_t
    }
    targ_use <- if (!is.null(spec$target_time)) {
      t_idx <- static_missing_match_panel_dates(panel, spec$target_time, sprintf("role_graphs[[%s]]$%s$target_time", cid, role))
      panel$dates_num[t_idx]
    } else {
      targ_t
    }
    if (base_use != base_t || targ_use != targ_t) {
      static_missing_abort(
        sprintf(
          "Role block `%s` (cell `%s`, role `%s`) has `(baseline_time = %g, target_time = %g)` which disagrees with cell transition dates `(%g, %g)`.",
          bid, cid, role, base_use, targ_use, base_t, targ_t
        ),
        "cdid_static_missing_alignment_error"
      )
    }
    b_spec <- static_role_block_spec(
      block_id = bid,
      cell_id = cid,
      role = role,
      baseline_time = base_use,
      target_time = targ_use,
      dates = dates_use,
      edges = edges_clean
    )
  } else {
    static_missing_abort(
      sprintf("Invalid role graph specification for cell `%s`, role `%s`.", cid, role),
      "cdid_static_missing_graph_error"
    )
  }

  # Validate all declared dates and edges against panel dates and create lazy pair accessor
  d_idx <- static_missing_match_panel_dates(panel, b_spec$dates, sprintf("role_block[[%s]]$dates", b_spec$block_id))
  K <- nrow(b_spec$edges)
  adapter_pair_keys <- paste0("e", seq_len(K))
  role_mask <- if (identical(role, "G")) cell$G else cell$C
  accessor <- static_observation_pair_accessor(
    panel = panel,
    mask = role_mask,
    s = b_spec$edges$s,
    u = b_spec$edges$u,
    pair_id = adapter_pair_keys
  )

  list(
    block_spec = b_spec,
    accessor = accessor,
    adapter_pair_keys = adapter_pair_keys,
    date_indices = d_idx
  )
}

static_missing_resolve_role_graphs <- function(panel,
                                               target_ids,
                                               all_ids,
                                               all_cells,
                                               role_graphs) {
  if (missing(role_graphs) || is.null(role_graphs)) {
    static_missing_abort(
      "Explicit `role_graphs` must be supplied; automatic graph selection is not permitted.",
      "cdid_static_missing_graph_error"
    )
  }

  raw_by_cell_role <- list()

  if (is.list(role_graphs) && length(role_graphs) >= 1L &&
      all(vapply(role_graphs, inherits, logical(1L), "cdid_static_role_block_spec"))) {
    for (b in role_graphs) {
      key <- paste0(b$cell_id, "::", b$role)
      if (!is.null(raw_by_cell_role[[key]])) {
        static_missing_abort(
          sprintf("Duplicate role block specification for cell `%s`, role `%s`.", b$cell_id, b$role),
          "cdid_static_graph_duplicate_block_error"
        )
      }
      raw_by_cell_role[[key]] <- b
    }
  } else if (is.data.frame(role_graphs)) {
    if (!all(c("cell_id", "role", "s", "u") %in% names(role_graphs))) {
      static_missing_abort(
        "A `data.frame` `role_graphs` specification must have columns `cell_id`, `role`, `s`, and `u`.",
        "cdid_static_missing_graph_error"
      )
    }
    roles_norm <- vapply(as.character(role_graphs$role), function(role) static_graph_normalize_role(role)$code, character(1L))
    split_keys <- paste0(as.character(role_graphs$cell_id), "::", roles_norm)
    for (sk in unique(split_keys)) {
      sub_df <- role_graphs[split_keys == sk, , drop = FALSE]
      raw_by_cell_role[[sk]] <- sub_df
    }
  } else if (is.list(role_graphs) && !is.object(role_graphs)) {
    c_names <- names(role_graphs)
    if (is.null(c_names) || anyNA(c_names) || any(!nzchar(c_names)) || anyDuplicated(c_names)) {
      static_missing_abort(
        "Named list `role_graphs` must be keyed by unique `cell_id` names.",
        "cdid_static_missing_graph_error"
      )
    }
    for (cid in c_names) {
      entry <- role_graphs[[cid]]
      if (inherits(entry, "cdid_static_role_block_spec")) {
        raw_by_cell_role[[paste0(cid, "::", entry$role)]] <- entry
      } else if (is.list(entry) && !is.data.frame(entry)) {
        g_spec <- if (!is.null(entry$G)) entry$G else entry$switcher
        c_spec <- if (!is.null(entry$C)) entry$C else entry$stayer
        if (!is.null(g_spec)) raw_by_cell_role[[paste0(cid, "::G")]] <- g_spec
        if (!is.null(c_spec)) raw_by_cell_role[[paste0(cid, "::C")]] <- c_spec
        if (is.null(g_spec) && is.null(c_spec) && (!is.null(entry$edges) || !is.null(entry$pairs))) {
          # For non-target OR training cells only C is required; if a single edge spec is given for a non-target cell, assign to C
          if (!(cid %in% target_ids) && (cid %in% all_ids)) {
            raw_by_cell_role[[paste0(cid, "::C")]] <- entry
          } else {
            static_missing_abort(
              sprintf("Cell `%s` in `role_graphs` must explicitly specify `$G` and `$C` role graphs.", cid),
              "cdid_static_missing_graph_error"
            )
          }
        }
      } else if ((is.data.frame(entry) || is.matrix(entry)) && !(cid %in% target_ids) && (cid %in% all_ids)) {
        raw_by_cell_role[[paste0(cid, "::C")]] <- entry
      } else {
        static_missing_abort(
          sprintf("Cell `%s` in `role_graphs` must be a list with `$G` and `$C` role graph specifications.", cid),
          "cdid_static_missing_graph_error"
        )
      }
    }
  } else {
    static_missing_abort(
      "`role_graphs` must be a list of `static_role_block_spec` objects, a named list keyed by `cell_id`, or a `data.frame(cell_id, role, s, u)`.",
      "cdid_static_missing_graph_error"
    )
  }

  # Validate switcher graphs for all target cells in A and stayer graphs for all cells in U
  G_resolved <- vector("list", length(target_ids))
  names(G_resolved) <- target_ids
  for (cid in target_ids) {
    key_g <- paste0(cid, "::G")
    if (is.null(raw_by_cell_role[[key_g]])) {
      static_missing_abort(
        sprintf("Target cell `%s` is missing an explicit switcher (`G`) role graph.", cid),
        c("cdid_static_graph_role_mismatch_error", "cdid_static_missing_graph_error")
      )
    }
    G_resolved[[cid]] <- static_missing_coerce_one_role_block(
      panel = panel,
      cell = all_cells[[cid]],
      role = "G",
      spec = raw_by_cell_role[[key_g]],
      default_block_id = paste0(cid, "__G")
    )
  }

  C_resolved <- vector("list", length(all_ids))
  names(C_resolved) <- all_ids
  for (cid in all_ids) {
    key_c <- paste0(cid, "::C")
    if (is.null(raw_by_cell_role[[key_c]])) {
      static_missing_abort(
        sprintf("Cell `%s` is missing an explicit stayer (`C`) role graph.", cid),
        c("cdid_static_graph_role_mismatch_error", "cdid_static_missing_graph_error")
      )
    }
    C_resolved[[cid]] <- static_missing_coerce_one_role_block(
      panel = panel,
      cell = all_cells[[cid]],
      role = "C",
      spec = raw_by_cell_role[[key_c]],
      default_block_id = paste0(cid, "__C")
    )
  }

  list(G = G_resolved, C = C_resolved)
}

static_missing_extract_one_cell_design <- function(panel,
                                                   cell,
                                                   role,
                                                   mask,
                                                   formula_arg = NULL,
                                                   design_arg = NULL,
                                                   require_intercept = TRUE) {
  cid <- cell$cell_id
  n <- panel$n
  uid_chr <- as.character(panel$unit_id)
  active_idx <- which(as.logical(mask))

  if (!is.null(design_arg)) {
    des_c <- if (inherits(design_arg, "cdid_static_reference_design")) {
      design_arg
    } else if (is.list(design_arg) && !is.data.frame(design_arg) && (cid %in% names(design_arg))) {
      design_arg[[cid]]
    } else {
      design_arg
    }

    if (inherits(des_c, "cdid_static_reference_design")) {
      if (!identical(des_c$n, n) || !identical(as.character(des_c$unit_id), uid_chr)) {
        static_missing_abort(
          sprintf("Cell `%s` `%s` design `unit_id` does not match `panel$unit_id`.", cid, role),
          c("cdid_static_missing_alignment_error", "cdid_static_order_error")
        )
      }
      X_tilde <- des_c$X_tilde
      Z <- des_c$Z
    } else if (is.list(des_c) && (!is.null(des_c$X_tilde) || !is.null(des_c$Z))) {
      X_tilde <- des_c$X_tilde
      Z <- des_c$Z
      if (is.null(Z) && !is.null(X_tilde)) Z <- X_tilde[, -1L, drop = FALSE]
      if (is.null(X_tilde) && !is.null(Z)) {
        X_tilde <- matrix(NA_real_, nrow = n, ncol = 1L + ncol(Z),
                          dimnames = list(uid_chr, c("(Intercept)", colnames(Z))))
        X_tilde[active_idx, 1L] <- 1.0
        if (ncol(Z) > 0L) X_tilde[active_idx, -1L] <- Z[active_idx, , drop = FALSE]
      }
    } else if (is.matrix(des_c) && is.numeric(des_c)) {
      if (nrow(des_c) != n) {
        static_missing_abort(
          sprintf("Cell `%s` `%s` design matrix must have `nrow == %d`.", cid, role, n),
          "cdid_static_design_error"
        )
      }
      if (!is.null(rownames(des_c)) && !identical(as.character(rownames(des_c)), uid_chr)) {
        static_missing_abort(
          sprintf("Cell `%s` `%s` design matrix rownames do not match `panel$unit_id`.", cid, role),
          "cdid_static_order_error"
        )
      }
      cn <- colnames(des_c)
      has_int <- (!is.null(cn) && length(cn) >= 1L && cn[1L] == "(Intercept)") ||
        (require_intercept && ncol(des_c) >= 1L && is.null(cn) &&
         length(active_idx) > 0L && all(!is.na(des_c[active_idx, 1L])) && all(des_c[active_idx, 1L] == 1))
      if (has_int) {
        X_tilde <- des_c
        if (is.null(colnames(X_tilde))) {
          colnames(X_tilde) <- c("(Intercept)", if (ncol(X_tilde) > 1L) paste0("z", seq_len(ncol(X_tilde) - 1L)) else character(0L))
        }
        Z <- X_tilde[, -1L, drop = FALSE]
      } else {
        Z <- des_c
        if (ncol(Z) > 0L && is.null(colnames(Z))) {
          colnames(Z) <- paste0("z", seq_len(ncol(Z)))
        }
        X_tilde <- matrix(NA_real_, nrow = n, ncol = 1L + ncol(Z),
                          dimnames = list(uid_chr, c("(Intercept)", colnames(Z))))
        X_tilde[active_idx, 1L] <- 1.0
        if (ncol(Z) > 0L) X_tilde[active_idx, -1L] <- Z[active_idx, , drop = FALSE]
      }
    } else {
      static_missing_abort(
        sprintf("Invalid `%s` design supplied for cell `%s`.", role, cid),
        "cdid_static_design_error"
      )
    }

    if (length(active_idx) > 0L) {
      X_act <- X_tilde[active_idx, , drop = FALSE]
      Z_act <- Z[active_idx, , drop = FALSE]
      if (anyNA(X_act) || any(!is.finite(X_act)) || anyNA(Z_act) || any(!is.finite(Z_act))) {
        static_missing_abort(
          sprintf("Cell `%s` (`%s` design) has missing or non-finite reference covariate values on active mask rows.", cid, role),
          "cdid_static_missing_covariate_error"
        )
      }
      if (!all(X_act[, 1L] == 1)) {
        static_missing_abort(
          sprintf("Cell `%s` (`%s` design) must include a leading unit `(Intercept)` column on active rows.", cid, role),
          "cdid_static_no_intercept_error"
        )
      }
    }
    rownames(X_tilde) <- uid_chr
    rownames(Z) <- uid_chr
    return(list(X_tilde = X_tilde, Z = Z, p = ncol(Z), z_names = colnames(Z)))
  }

  if (is.null(formula_arg)) {
    static_missing_abort(
      sprintf("Cell `%s` requires an explicit `%s` formula or pre-extracted design.", cid, role),
      "cdid_static_formula_error"
    )
  }
  f_use <- if (is.list(formula_arg) && !inherits(formula_arg, "formula")) {
    if (!(cid %in% names(formula_arg))) {
      static_missing_abort(
        sprintf("Formula list for `%s` is missing entry for cell `%s`.", role, cid),
        "cdid_static_formula_error"
      )
    }
    formula_arg[[cid]]
  } else {
    formula_arg
  }

  des <- static_extract_reference_design(
    panel = panel,
    cell = cell,
    formula = f_use,
    role = role,
    mask = mask
  )
  list(X_tilde = des$X_tilde, Z = des$Z, p = des$p_slopes, z_names = des$slope_names)
}

static_missing_resolve_rows <- function(rows, unit_id) {
  n <- length(unit_id)
  if (is.null(rows)) return(seq_len(n))
  if (!is.atomic(rows) || !is.null(dim(rows)) || is.object(rows) || length(rows) < 1L || anyNA(rows)) {
    static_missing_abort(
      "`rows` must be a non-empty vector of integer positions in `1..n`, logical mask of length `n`, or original `unit_id` values.",
      "cdid_static_missing_input_error"
    )
  }
  if (is.logical(rows)) {
    if (length(rows) != n) {
      static_missing_abort(
        sprintf("Logical `rows` mask must have length `n = %d`.", n),
        "cdid_static_missing_input_error"
      )
    }
    return(which(rows))
  }
  if (is.character(rows)) {
    uid_chr <- as.character(unit_id)
    idx <- match(rows, uid_chr)
    if (anyNA(idx) || anyDuplicated(idx)) {
      static_missing_abort(
        "Character `rows` must be unique original `unit_id` values present in `panel$unit_id`.",
        "cdid_static_missing_input_error"
      )
    }
    return(as.integer(idx))
  }
  if (is.numeric(rows)) {
    if (any(!is.finite(rows))) {
      static_missing_abort("Numeric `rows` must be finite.", "cdid_static_missing_input_error")
    }
    if (is.numeric(unit_id) && all(rows %in% unit_id) && !all(rows %in% seq_len(n))) {
      idx <- match(as.double(rows), as.double(unit_id))
      if (anyNA(idx) || anyDuplicated(idx)) {
        static_missing_abort("Numeric `rows` must be unique `unit_id` values or integer positions in `1..n`.", "cdid_static_missing_input_error")
      }
      return(as.integer(idx))
    }
    if (any(rows != floor(rows) | rows < 1 | rows > n) || anyDuplicated(as.integer(rows))) {
      static_missing_abort(
        sprintf("Integer `rows` positions must be unique values in `1..%d`.", n),
        "cdid_static_missing_input_error"
      )
    }
    return(as.integer(rows))
  }
  static_missing_abort("Unsupported `rows` specification.", "cdid_static_missing_input_error")
}

# Wrap a `cdid_static_observation_moment` object (whose `$eval_if(rows)` expects
# original `unit_id` values) with `$if_fn(rows)` accepting integer indices `1..n`
# or original `unit_id` values as required by `static_reconstructed_or`.
static_missing_wrap_adapter_moment <- function(moment_obj, unit_id) {
  n <- length(unit_id)
  raw_eval_if <- moment_obj$eval_if
  moment_obj$if_fn <- function(rows = seq_len(n)) {
    idx <- static_missing_resolve_rows(rows, unit_id)
    raw_eval_if(rows = unit_id[idx])
  }
  moment_obj
}

# Uncentered IPW cell link container matching the interface of
# `static_reconstructed_or_center_cell_links` for uniform downstream consumption.
static_missing_ipw_cell_links <- function(cell_id,
                                          unit_id,
                                          t,
                                          G_edges,
                                          G_links,
                                          C_edges,
                                          C_links) {
  uid_chr <- as.character(unit_id)
  n <- length(uid_chr)
  K_G <- nrow(G_edges)
  K_C <- nrow(C_edges)
  g_edge_ids <- paste0(G_edges[, 1L], "->", G_edges[, 2L])
  c_edge_ids <- paste0(C_edges[, 1L], "->", C_edges[, 2L])

  raw_G <- stats::setNames(vapply(G_links, `[[`, numeric(1L), "estimate"), g_edge_ids)
  raw_C <- stats::setNames(vapply(C_links, `[[`, numeric(1L), "estimate"), c_edge_ids)
  stacked_links <- c(
    stats::setNames(raw_G, paste0("G:", g_edge_ids)),
    stats::setNames(raw_C, paste0("C:", c_edge_ids))
  )

  if_callback <- function(rows = seq_len(n)) {
    idx <- static_missing_resolve_rows(rows, unit_id)
    m <- length(idx)
    row_ids <- uid_chr[idx]
    chi_G <- matrix(0.0, nrow = m, ncol = K_G, dimnames = list(row_ids, g_edge_ids))
    for (k in seq_len(K_G)) {
      chi_G[, k] <- G_links[[k]]$if_fn(idx)
    }
    chi_C <- matrix(0.0, nrow = m, ncol = K_C, dimnames = list(row_ids, c_edge_ids))
    for (k in seq_len(K_C)) {
      chi_C[, k] <- C_links[[k]]$if_fn(idx)
    }
    chi_stacked <- cbind(chi_G, chi_C)
    colnames(chi_stacked) <- names(stacked_links)
    list(
      rows = idx,
      unit_id = unit_id[idx],
      chi_raw_G = chi_G,
      chi_raw_C = chi_C,
      chi_center_G = numeric(m),
      chi_center_C = numeric(m),
      chi_adj_G = chi_G,
      chi_adj_C = chi_C,
      chi_stacked = chi_stacked
    )
  }

  structure(
    list(
      cell_id = cell_id,
      block_id = NULL,
      method = "ipw",
      paper_method = "ipw",
      t = as.double(t),
      alpha_c = 0.0,
      b_P = numeric(0L),
      zG = NULL,
      zCr = NULL,
      center_G = 0.0,
      center_C = 0.0,
      step_G = stats::setNames(numeric(K_G), g_edge_ids),
      step_C = stats::setNames(numeric(K_C), c_edge_ids),
      raw_G = raw_G,
      raw_C = raw_C,
      adjusted_G = raw_G,
      adjusted_C = raw_C,
      stacked_links = stacked_links,
      if_callback = if_callback,
      influence = NULL
    ),
    class = "cdid_static_centered_cell_links"
  )
}

static_missing_wrap_centered_cell_links <- function(cc_obj, unit_id) {
  force(cc_obj)
  force(unit_id)
  n <- length(unit_id)
  raw_cc_cb <- cc_obj$if_callback
  cc_obj$if_callback <- function(rows = seq_len(n)) {
    idx <- static_missing_resolve_rows(rows, unit_id)
    raw_cc_cb(idx)
  }
  cc_obj
}

# Simultaneous workspace admission for the optional adapter-owned batch holder
# of built-in selection (`q`) and imputation (`j`) coefficient IF matrices and
# the fit-owned compact batch holder of reconstructed OR IF (`chi_alpha`, `chi_b`)
# arrays. Composes the live consumer reservation (max of moment-covariance
# streaming and joint-GMM replay phases), the static producer scratch reservation
# inside `link_if_reader`, the adapter q/j holder reservation (`B * (Pq + Pj)`
# plus owned row key, dims, and dimnames metadata), and the compact OR holder
# reservation (`B * (U + sum(p_block))` plus owned row key, dims, dimnames, and
# list metadata) against `bounds$max_matrix_entries`.
# If an additional holder is not funded, only that optimization is privately
# declined without altering chunk sizes, public bounds, or scientific behavior.
static_missing_admit_nuisance_batch_reuse <- function(
  n,
  batch_rows,
  bounds,
  obs_adapter,
  gmm_blocks,
  m_R,
  target_ids,
  all_ids,
  m_info,
  or_fit,
  or_designs,
  trt_X_safe_by_cell,
  graph_res,
  w_choice,
  accuracy,
  cluster,
  materialize_target_if,
  aggregate_cells
) {
  sel_fit <- obs_adapter$selection_fit
  imp_fit <- obs_adapter$imputation_fit
  Pq <- if (!is.null(sel_fit)) length(sel_fit$layout$param_keys) else 0L
  Pj <- if (!is.null(imp_fit)) length(imp_fit$layout$param_keys) else 0L
  P_tot <- as.integer(Pq + Pj)
  use_or_fit <- isTRUE(m_info$use_or) && !is.null(or_fit)
  n_all_cells <- if (use_or_fit) length(all_ids) else 0L
  block_ps <- if (use_or_fit) unlist(or_fit$partition$block_p, use.names = FALSE) else integer(0L)
  n_or_blocks <- length(block_ps)
  sum_p_or_blocks <- if (n_or_blocks > 0L) sum(as.integer(block_ps)) else 0L
  max_e <- as.double(bounds$max_matrix_entries)
  B <- as.integer(batch_rows)

  if ((P_tot <= 0L && !use_or_fit) || B <= 0L) {
    return(list(
      admitted = FALSE,
      qj_admitted = FALSE,
      or_admitted = FALSE,
      Pq = as.integer(Pq),
      Pj = as.integer(Pj),
      U_or = as.integer(n_all_cells),
      n_or_blocks = as.integer(n_or_blocks),
      sum_p_or_blocks = as.integer(sum_p_or_blocks),
      batch_rows = B,
      qj_holder_entries = 0.0,
      or_holder_entries = 0.0,
      holder_entries = 0.0,
      producer_entries = 0.0,
      consumer_entries = 0.0,
      cov_consumer_entries = 0.0,
      gmm_consumer_entries = 0.0,
      qj_total_entries = 0.0,
      or_total_entries = 0.0,
      total_entries = 0.0,
      max_matrix_entries = max_e
    ))
  }

  tryCatch({
    max_safe <- 2^53 - 1
    P <- function(a, b, label) gmm_safe_product(a, b, max_safe, label)
    S <- function(vals, label) {
      tot <- sum(as.double(vals))
      if (!is.finite(tot) || tot > max_safe) {
        cdid_gmm_abort(paste0(label, " exceeds `bounds$max_matrix_entries`."), "cdid_gmm_bound_error")
      }
      tot
    }

    # 1a. Adapter-owned q/j coefficient IF batch holder reservation:
    #     - B * (Pq + Pj) numeric entries for at most one selection and one imputation IF matrix
    #     - B + 4 entries for the owned integer row-position key and two matrix dim vectors
    #     - owned row/column-name vectors and small holder metadata per active fit
    if (P_tot > 0L) {
      n_active_fits <- as.integer(Pq > 0L) + as.integer(Pj > 0L)
      holder_numeric <- P(B, P_tot, "nuisance IF batch holder matrices")
      holder_key_dim <- as.double(B) + 4.0
      holder_meta <- P(n_active_fits, B, "nuisance IF batch holder rownames") + as.double(P_tot) + 8.0
      qj_holder_entries <- S(c(holder_numeric, holder_key_dim, holder_meta), "new holder reservation")
    } else {
      qj_holder_entries <- 0.0
    }

    # 1b. Fit-owned compact reconstructed-OR batch holder reservation (when OR/DR active):
    #     - B * (U + sum(p_block)) numeric entries for `chi_alpha` (B x U) and `chi_b` blocks (B x p_b)
    #     - B + 2 * (1 + n_or_blocks) entries for the owned integer row-position key and matrix dim vectors
    #     - owned row/column-name vectors and list metadata for `chi_alpha` and `chi_b`
    if (use_or_fit) {
      or_holder_numeric <- P(B, n_all_cells + sum_p_or_blocks, "reconstructed OR batch holder matrices")
      or_holder_key_dim <- as.double(B) + 2.0 * (1.0 + as.double(n_or_blocks))
      or_holder_meta <- P(1L + n_or_blocks, B, "reconstructed OR batch holder rownames") +
        as.double(n_all_cells + sum_p_or_blocks) + 2.0 * as.double(n_or_blocks) + 8.0
      or_holder_entries <- S(
        c(or_holder_numeric, or_holder_key_dim, or_holder_meta),
        "reconstructed OR batch holder reservation"
      )
    } else {
      or_holder_entries <- 0.0
    }

    # 2. Static producer scratch reservation coexisting during `link_if_reader(idx)`:
    #    - Assembled link IF matrix `out` (B x m_R) allocated inside the `with_if_batch`
    #      callback `fn`, plus row index/name vectors
    #    - On later row batches (`n > B`), the unchanged streaming wrapper
    #      `make_static_if_batch_callback` (`R/static_graph_gmm.R:970-978`) still
    #      retains the previous batch's shared `phi_b` / `batch` (`B x m_R` plus
    #      `unit_id` and row names) while evaluating `reader_fn(idx)`
    link_out_entries <- P(B, m_R, "link_if_reader batch output") + 2.0 * as.double(B)
    prev_batch_entries <- if (n > B) {
      P(B, m_R, "previous streamed IF batch retained during reader_fn") + 2.0 * as.double(B)
    } else {
      0.0
    }

    #    - Active target-cell centering callback (`chi_raw_G`, `chi_raw_C`, `chi_adj_G`,
    #      `chi_adj_C`, `chi_stacked`, `chi_center_G`, `chi_center_C`); `slice_c` is
    #      released before evaluating the next cell's callback
    K_cell_max <- max(vapply(target_ids, function(cid) {
      k_g <- nrow(graph_res$G[[cid]]$block_spec$edges)
      k_c <- if (m_info$gmm_method %in% c("ipw", "dr")) nrow(graph_res$C[[cid]]$block_spec$edges) else 0L
      as.integer(k_g + k_c)
    }, integer(1L)))
    cell_center_entries <- P(B, 3L * K_cell_max + 2L, "target-cell link centering scratch")

    #    - Reconstructed OR batch callback (`or_fit$if_callback`) when active
    if (use_or_fit) {
      sum_p_or_cells <- sum(vapply(all_ids, function(cid) ncol(or_designs[[cid]]$Z), integer(1L)))
      max_p_or <- if (n_or_blocks > 0L) max(0L, as.integer(block_ps)) else 0L
      or_entries <- P(
        B,
        3L * n_all_cells + 4L * sum_p_or_cells + 3L * sum_p_or_blocks + 5L * max_p_or,
        "reconstructed OR batch IF scratch"
      )
    } else {
      max_p_or <- 0L
      or_entries <- 0.0
    }

    #    - Active moment `eval_if` direct output, requested channel slices, contraction,
    #      and `static_nuisance_if` current-block scratch during second matrix construction
    d_out_max <- max(1L, max_p_or)
    p_trt_max <- if (isTRUE(m_info$use_logit) && length(trt_X_safe_by_cell)) {
      max(vapply(trt_X_safe_by_cell, ncol, integer(1L)))
    } else 0L
    moment_direct_entries <- P(B, 6L * d_out_max + 12L, "moment direct and contraction IF scratch")
    channel_slice_entries <- P(B, P_tot + p_trt_max, "moment nuisance channel slices")
    sel_block_cols <- if (!is.null(sel_fit)) {
      vapply(sel_fit$layout$blocks, function(bl) length(bl$index), integer(1L))
    } else integer(0L)
    imp_block_cols <- if (!is.null(imp_fit)) {
      vapply(imp_fit$layout$blocks, function(bl) length(bl$index), integer(1L))
    } else integer(0L)
    max_nuis_block_cols <- max(0L, sel_block_cols, imp_block_cols)
    kernel_block_scratch <- P(B, 4L * max_nuis_block_cols + 4L, "static_nuisance_if block scratch")

    producer_entries <- S(
      c(
        link_out_entries,
        prev_batch_entries,
        cell_center_entries,
        or_entries,
        moment_direct_entries,
        channel_slice_entries,
        kernel_block_scratch
      ),
      "static producer scratch reservation"
    )

    # 3. Live consumer reservation: maximum of moment-covariance streaming and
    #    joint-GMM solve / IF replay phases using their actual dimensions/options.
    cov_consumer_entries <- S(
      c(
        P(m_R, m_R, "moment covariance accumulator S"),
        P(2L * B, m_R, "covariance stream IF batch and centered batch"),
        P(7L, n, "covariance and cross-stage replay unit vectors"),
        2.0 * as.double(m_R) + 6.0 * as.double(B)
      ),
      "covariance consumer reservation"
    )

    p_gmm <- sum(vapply(gmm_blocks, function(blk) length(blk$dates) - 1L, integer(1L)))
    k_gmm <- sum(vapply(gmm_blocks, function(blk) length(blk$dates) - as.integer(blk$n_components), integer(1L)))
    J_est <- length(target_ids)

    acc_use <- accuracy
    if (is.list(acc_use) && inherits(acc_use$covariance_strategy, "cdid_gmm_covariance_strategy")) {
      acc_use$covariance_strategy$weighting_requested <- w_choice$requested
      acc_use$covariance_strategy$weighting_canonical <- w_choice$canonical
      acc_use$covariance_strategy$alias_used <- w_choice$alias_used
    }
    acc_use <- validate_gmm_accuracy_control(acc_use, w_choice$canonical)
    if (w_choice$canonical == "optimal" && is.null(acc_use$covariance_strategy)) {
      acc_use$covariance_strategy <- gmm_covariance_strategy(
        method = "tikhonov",
        condition_trigger = 1e-10,
        ridge_relative = 1e-4,
        weighting_requested = w_choice$requested,
        weighting_canonical = w_choice$canonical
      )
    }

    gmm_pre <- preflight_gmm_allocations(
      m_R = m_R,
      k = k_gmm,
      J_est = J_est,
      weighting = w_choice$canonical,
      bounds = bounds,
      n_out_rows = n,
      request_param_cov = FALSE,
      materialize_target_if = materialize_target_if,
      n_input_rows = n,
      n_parameters = p_gmm,
      certify_error = !identical(acc_use$method, "conventional_chol_qr"),
      covariance_strategy = acc_use$covariance_strategy,
      n_dropped_triplets = 0L
    )

    outer_gmm_entries <- S(
      c(
        P(m_R, m_R, "outer moment covariance V_use"),
        P(3L, n, "cross-stage replay verification vectors"),
        3.0 * as.double(B)
      ),
      "outer GMM wrapper reservation"
    )

    G_cl <- if (!is.null(cluster)) length(unique(as.character(cluster))) else 0L
    cluster_entries <- if (G_cl >= 2L) {
      S(
        c(
          P(2L * G_cl, J_est, "cluster IF accumulators"),
          P(J_est, J_est, "cluster sandwich covariance"),
          P(B, J_est, "cluster projected IF batch"),
          as.double(n) + as.double(G_cl)
        ),
        "cluster sandwich reservation"
      )
    } else 0.0

    agg_entries <- if (!is.null(aggregate_cells)) {
      S(
        c(
          P(3L, n, "aggregate IF vectors"),
          P(B, J_est, "aggregate cell_if_callback projected batch")
        ),
        "aggregate replay reservation"
      )
    } else 0.0

    gmm_consumer_entries <- S(
      c(gmm_pre$total_payload_entries, outer_gmm_entries, max(cluster_entries, agg_entries)),
      "GMM consumer reservation"
    )

    consumer_entries <- max(cov_consumer_entries, gmm_consumer_entries)
    base_entries <- S(
      c(consumer_entries, producer_entries),
      "simultaneous consumer + producer reservation"
    )
    qj_total_entries <- S(
      c(base_entries, qj_holder_entries),
      "simultaneous consumer + producer + q/j holder reservation"
    )
    or_total_entries <- S(
      c(base_entries, qj_holder_entries, or_holder_entries),
      "simultaneous consumer + producer + q/j holder + OR holder reservation"
    )
    holder_entries <- S(
      c(qj_holder_entries, or_holder_entries),
      "combined nuisance and OR holder reservation"
    )
    total_entries <- or_total_entries

    qj_admitted <- (P_tot > 0L) && is.finite(qj_total_entries) && (qj_total_entries <= max_e)
    or_admitted <- use_or_fit && is.finite(or_total_entries) && (or_total_entries <= max_e)

    list(
      admitted = qj_admitted,
      qj_admitted = qj_admitted,
      or_admitted = or_admitted,
      Pq = as.integer(Pq),
      Pj = as.integer(Pj),
      U_or = as.integer(n_all_cells),
      n_or_blocks = as.integer(n_or_blocks),
      sum_p_or_blocks = as.integer(sum_p_or_blocks),
      batch_rows = B,
      qj_holder_entries = qj_holder_entries,
      or_holder_entries = or_holder_entries,
      holder_entries = holder_entries,
      producer_entries = producer_entries,
      consumer_entries = consumer_entries,
      cov_consumer_entries = cov_consumer_entries,
      gmm_consumer_entries = gmm_consumer_entries,
      qj_total_entries = qj_total_entries,
      or_total_entries = or_total_entries,
      total_entries = total_entries,
      max_matrix_entries = max_e
    )
  }, error = function(e) {
    list(
      admitted = FALSE,
      qj_admitted = FALSE,
      or_admitted = FALSE,
      Pq = as.integer(Pq),
      Pj = as.integer(Pj),
      U_or = as.integer(n_all_cells),
      n_or_blocks = as.integer(n_or_blocks),
      sum_p_or_blocks = as.integer(sum_p_or_blocks),
      batch_rows = B,
      qj_holder_entries = NA_real_,
      or_holder_entries = NA_real_,
      holder_entries = NA_real_,
      producer_entries = NA_real_,
      consumer_entries = NA_real_,
      cov_consumer_entries = NA_real_,
      gmm_consumer_entries = NA_real_,
      qj_total_entries = NA_real_,
      or_total_entries = NA_real_,
      total_entries = NA_real_,
      max_matrix_entries = max_e
    )
  })
}

# Main private end-to-end static missing-outcome assembly entrypoint
static_missing_fit <- function(panel,
                               target_cells,
                               role_graphs,
                               method,
                               route,
                               or_blocks = NULL,
                               all_cells = NULL,
                               formula = NULL,
                               treatment_formula = NULL,
                               outcome_formula = NULL,
                               selection_formula = NULL,
                               imputation_formula = NULL,
                               treatment_design = NULL,
                               outcome_design = NULL,
                               selection_design = NULL,
                               imputation_design = NULL,
                               selection_model = NULL,
                               selection_pooling = NULL,
                               imputation_pooling = NULL,
                               imputation_scale = NULL,
                               weighting = "optimal",
                               glm_control = list(epsilon = 1e-12, maxit = 25L, trace = FALSE),
                               rank_tol = 1e-7,
                               bounds = NULL,
                               accuracy = NULL,
                               cluster = NULL,
                               centered = FALSE,
                               aggregate_cells = NULL,
                               aggregate_signs = NULL,
                               materialize_if = TRUE,
                               idname = NULL,
                               tname = NULL,
                               yname = NULL,
                               dname = NULL,
                               dates = NULL,
                               unit_ids = NULL) {
  # 1. Validate method, route, panel, transition cells, and role graphs
  m_info <- static_missing_validate_method(method)
  r_info <- static_missing_validate_route(route)
  panel_obj <- static_missing_resolve_panel(
    panel = panel,
    idname = idname,
    tname = tname,
    yname = yname,
    dname = dname,
    dates = dates,
    unit_ids = unit_ids
  )
  n <- panel_obj$n
  n_num <- as.double(n)
  unit_id <- panel_obj$unit_id
  unit_id_chr <- as.character(unit_id)

  cell_res <- static_missing_resolve_cells_and_partition(
    panel = panel_obj,
    target_cells = target_cells,
    all_cells = all_cells,
    or_blocks = or_blocks,
    use_or = m_info$use_or
  )
  target_ids <- cell_res$target_ids
  target_map <- cell_res$target_cells
  all_ids <- cell_res$all_ids
  all_map <- cell_res$all_cells

  graph_res <- static_missing_resolve_role_graphs(
    panel = panel_obj,
    target_ids = target_ids,
    all_ids = all_ids,
    all_cells = all_map,
    role_graphs = role_graphs
  )

  # Validate bounds and accuracy controls for joint GMM
  w_choice <- resolve_gmm_weighting_choice(weighting, arg = "weighting")
  if (is.null(bounds)) {
    static_missing_abort(
      "Explicit `bounds` (`gmm_bounds_control(...)`) must be supplied.",
      c("cdid_gmm_bound_error", "cdid_static_missing_control_error")
    )
  }
  bounds <- validate_gmm_bounds_control(bounds)
  if (is.null(accuracy)) {
    static_missing_abort(
      "Explicit `accuracy` (`static_graph_gmm_accuracy_control(...)` or `gmm_accuracy_control(...)`) must be supplied.",
      c("cdid_gmm_control_error", "cdid_static_missing_control_error")
    )
  }
  batch_rows <- as.integer(min(n, bounds$max_batch_rows))

  # 2. Extract and validate reference covariate designs ONLY for active nuisances
  trt_formula_use <- if (!is.null(treatment_formula)) treatment_formula else formula
  or_formula_use <- if (!is.null(outcome_formula)) outcome_formula else formula
  sel_formula_use <- if (!is.null(selection_formula)) selection_formula else formula
  imp_formula_use <- if (!is.null(imputation_formula)) imputation_formula else formula

  # Even when route == "baseline", reference covariates at t - 1 must be validated
  # on R_c for every cell in U (via treatment design for IPW, or outcome design for OR/DR).
  trt_designs <- list()
  if (m_info$use_logit) {
    for (cid in target_ids) {
      cl <- target_map[[cid]]
      trt_designs[[cid]] <- static_missing_extract_one_cell_design(
        panel = panel_obj,
        cell = cl,
        role = "treatment",
        mask = (cl$R == 1L),
        formula_arg = trt_formula_use,
        design_arg = treatment_design,
        require_intercept = TRUE
      )
    }
  }

  or_designs <- list()
  if (m_info$use_or) {
    for (cid in all_ids) {
      cl <- all_map[[cid]]
      # Target cells need Z_or on R_c = G_c + C_c (stayer OR + switcher/stayer centering);
      # non-target training cells in U \ A also validate reference covariates on R_c (S0: X_c on R_c).
      or_designs[[cid]] <- static_missing_extract_one_cell_design(
        panel = panel_obj,
        cell = cl,
        role = "outcome",
        mask = (cl$R == 1L),
        formula_arg = or_formula_use,
        design_arg = outcome_design,
        require_intercept = TRUE
      )
    }
  }

  # 3. Treatment propensity logit fits on FULL R_c (only when method %in% c("ipw", "dr"))
  logit_fits <- list()
  trt_if_map <- list()
  trt_h_by_cell <- list()
  trt_jac_by_cell <- list()
  trt_coef_if_by_cell <- list()
  trt_X_safe_by_cell <- list()
  trt_active_by_cell <- list()

  if (m_info$use_logit) {
    for (cid in target_ids) {
      cl <- target_map[[cid]]
      td <- trt_designs[[cid]]
      # Outside-risk entries set to neutral 0.0 solely for full-N score matrix multiplication
      X_safe <- td$X_tilde
      X_safe[cl$R == 0L, ] <- 0.0
      colnames(X_safe) <- colnames(td$X_tilde)
      rownames(X_safe) <- unit_id_chr

      thin_cell <- list(
        cell_id = cid,
        n = n,
        R = cl$R,
        G = cl$G,
        C = cl$C,
        n_switchers = cl$n_switchers,
        n_stayers = cl$n_stayers,
        X_tilde = X_safe
      )
      lfit <- static_complete_fit_cell_logit(
        cell = thin_cell,
        glm_control = glm_control,
        rank_tol = rank_tol
      )
      ell_g <- lfit$ell_gamma
      dimnames(ell_g) <- list(unit_id_chr, colnames(X_safe))
      lfit$ell_gamma <- ell_g
      logit_fits[[cid]] <- lfit

      trt_key <- paste0("treatment:", cid)
      trt_if_map[[trt_key]] <- ell_g
      h_c <- stats::setNames(as.double(lfit$odds), unit_id_chr)
      act_c <- stats::setNames(as.double(lfit$active), unit_id_chr)
      rdot_c <- X_safe * (as.double(lfit$active) * as.double(lfit$odds))
      dimnames(rdot_c) <- list(unit_id_chr, colnames(X_safe))

      trt_h_by_cell[[cid]] <- h_c
      trt_active_by_cell[[cid]] <- act_c
      trt_X_safe_by_cell[[cid]] <- X_safe
      trt_jac_by_cell[[cid]] <- stats::setNames(list(rdot_c), trt_key)
      trt_coef_if_by_cell[[cid]] <- stats::setNames(list(ell_g), trt_key)
    }
  }

  # 4. Build observation adapter groups and fit q / j ONCE across all declared role populations
  # Active role groups:
  #   - For every target cell c in A: switcher group (cid, "G") and stayer group (cid, "C").
  #   - For every non-target cell c in U \ A (when use_or): stayer group (cid, "C").
  obs_groups <- list()
  group_accessors <- list()

  add_obs_group <- function(cid, role, g_resolved) {
    cl <- all_map[[cid]]
    mask_int <- if (role == "G") cl$G else cl$C
    mask_log <- (mask_int == 1L)
    g_key <- paste0(cid, "__", role)
    b_spec <- g_resolved$block_spec
    acc <- g_resolved$accessor
    p_keys <- g_resolved$adapter_pair_keys
    K <- length(p_keys)

    Z_sel_g <- NULL
    if (r_info$uses_q) {
      des_sel <- static_missing_extract_one_cell_design(
        panel = panel_obj,
        cell = cl,
        role = "selection",
        mask = mask_log,
        formula_arg = sel_formula_use,
        design_arg = selection_design,
        require_intercept = TRUE
      )
      Z_sel_g <- des_sel$Z
    }

    Z_imp_g <- NULL
    if (r_info$uses_j) {
      des_imp <- static_missing_extract_one_cell_design(
        panel = panel_obj,
        cell = cl,
        role = "imputation",
        mask = mask_log,
        formula_arg = imp_formula_use,
        design_arg = imputation_design,
        require_intercept = TRUE
      )
      Z_imp_g <- des_imp$Z
    }

    pairs_list <- vector("list", K)
    names(pairs_list) <- p_keys
    for (k in seq_len(K)) {
      s_idx_k <- acc$summary_table$s_idx[k]
      u_idx_k <- acc$summary_table$u_idx[k]
      pairs_list[[k]] <- list(
        pair_key = p_keys[k],
        start = paste0("d", s_idx_k),
        end = paste0("d", u_idx_k),
        gap = acc$summary_table$gap[k]
      )
    }

    dates_list <- NULL
    if (r_info$uses_q && identical(selection_model, "independent_period_logit")) {
      # Preserve unique edge dates in order of appearance in declared pairs
      edge_date_indices <- unique(as.integer(rbind(acc$summary_table$s_idx, acc$summary_table$u_idx)))
      edge_date_keys <- paste0("d", edge_date_indices)
      dates_list <- lapply(seq_along(edge_date_indices), function(j_d) {
        d_idx <- edge_date_indices[j_d]
        list(
          date_key = edge_date_keys[j_d],
          observed = (panel_obj$O[, d_idx] == 1L)
        )
      })
      names(dates_list) <- edge_date_keys
    }

    group_accessors[[g_key]] <<- acc
    obs_groups[[g_key]] <<- list(
      group_key = g_key,
      role = role,
      cell_id = cid,
      target = cl$t_num,
      prior_state = cl$d,
      mask = mask_log,
      A = mask_log,
      selection_design = Z_sel_g,
      imputation_design = Z_imp_g,
      pairs = pairs_list,
      dates = dates_list
    )
  }

  for (cid in target_ids) {
    add_obs_group(cid, "G", graph_res$G[[cid]])
    add_obs_group(cid, "C", graph_res$C[[cid]])
  }
  for (cid in cell_res$non_target_ids) {
    add_obs_group(cid, "C", graph_res$C[[cid]])
  }

  pair_cb <- function(group_key, pair_key) {
    acc <- group_accessors[[group_key]]
    p_dat <- acc$get_pair(pair_key)
    list(
      observed = (p_dat$O == 1L),
      dy = p_dat$dy,
      gap = p_dat$gap,
      start = paste0("d", p_dat$s_idx),
      end = paste0("d", p_dat$u_idx)
    )
  }

  obs_adapter <- static_observation_adapter_fit(
    unit_id = unit_id,
    groups = obs_groups,
    route = r_info$route,
    selection_model = if (r_info$uses_q) selection_model else NULL,
    selection_pooling = if (r_info$uses_q) selection_pooling else NULL,
    imputation_pooling = if (r_info$uses_j) imputation_pooling else NULL,
    imputation_scale = if (r_info$uses_j) imputation_scale else NULL,
    glm_control = glm_control,
    rank_tol = rank_tol,
    max_chunk_rows = batch_rows,
    pair_callback = pair_cb
  )

  ones_N <- stats::setNames(rep(1.0, n), unit_id_chr)

  # 5. Reconstruct stayer OR across U and center cell links on A
  or_fit <- NULL
  raw_or_if_callback <- NULL
  or_fit_centering <- NULL
  or_batch_active <- FALSE
  or_batch_pos <- NULL
  or_batch_slice <- NULL
  centered_by_cell <- vector("list", length(target_ids))
  names(centered_by_cell) <- target_ids
  raw_links_by_cell <- vector("list", length(all_ids))
  names(raw_links_by_cell) <- all_ids

  if (m_info$use_or) {
    recon_cells <- vector("list", length(all_ids))
    names(recon_cells) <- all_ids

    for (cid in all_ids) {
      cl <- all_map[[cid]]
      g_res_C <- graph_res$C[[cid]]
      b_spec_C <- g_res_C$block_spec
      p_keys_C <- g_res_C$adapter_pair_keys
      K_C <- length(p_keys_C)
      g_key_C <- paste0(cid, "__C")
      Z_or_c <- or_designs[[cid]]$Z
      p_or_c <- ncol(Z_or_c)

      edge_means_c <- vector("list", K_C)
      edge_covs_c <- if (p_or_c > 0L) vector("list", K_C) else NULL

      for (k in seq_len(K_C)) {
        pk <- p_keys_C[k]
        m_k <- obs_adapter$mean(
          group = g_key_C,
          pair = pk,
          h = ones_N,
          h_jacobian = list(),
          h_coef_if = NULL
        )
        edge_means_c[[k]] <- static_missing_wrap_adapter_moment(m_k, unit_id)
        if (p_or_c > 0L) {
          c_k <- obs_adapter$cov(
            group = g_key_C,
            pair = pk,
            Z = Z_or_c
          )
          edge_covs_c[[k]] <- static_missing_wrap_adapter_moment(c_k, unit_id)
        }
      }

      raw_links_by_cell[[cid]] <- list(C_unweighted_means = edge_means_c, C_covs = edge_covs_c)
      recon_cells[[cid]] <- static_reconstructed_or_cell(
        cell_id = cid,
        unit_id = unit_id,
        C = cl$C,
        Z = Z_or_c,
        edges = b_spec_C$edges[, c("s", "u"), drop = FALSE],
        t = cl$t_num,
        baseline_date = cl$t_prev_num,
        dates = b_spec_C$dates,
        d = cl$d,
        edge_means = edge_means_c,
        edge_covs = edge_covs_c
      )
    }

    or_fit <- static_reconstructed_or_fit(
      cells = recon_cells,
      or_blocks = cell_res$or_blocks,
      rank_tol = rank_tol,
      unit_id = unit_id,
      target_cells = target_ids,
      chunk_size = batch_rows,
      materialize_if = FALSE
    )

    # Private internal view of `or_fit` for target-cell centering inside
    # `link_if_reader`. Retains the original `or_fit` and its complete callback
    # contract in `nuisance$or`, while sharing only `chi_alpha` and `chi_b`
    # across target cells within one admitted ordered row batch.
    raw_or_if_callback <- or_fit$if_callback
    or_fit_centering <- or_fit
    or_fit_centering$if_callback <- function(rows = seq_len(n), include_chi_A_tensor = FALSE) {
      if (!or_batch_active || !identical(include_chi_A_tensor, FALSE)) {
        return(raw_or_if_callback(rows = rows, include_chi_A_tensor = include_chi_A_tensor))
      }
      rows_int <- static_reconstructed_or_validate_rows(rows, n)
      if (!identical(rows_int, or_batch_pos)) {
        return(raw_or_if_callback(rows = rows_int, include_chi_A_tensor = FALSE))
      }
      if (is.null(or_batch_slice)) {
        full_slice <- raw_or_if_callback(rows = rows_int, include_chi_A_tensor = FALSE)
        or_batch_slice <<- list(
          chi_alpha = full_slice$chi_alpha,
          chi_b = full_slice$chi_b
        )
        full_slice <- NULL
      }
      or_batch_slice
    }

    for (cid in target_ids) {
      cl <- target_map[[cid]]
      Z_or_c <- or_designs[[cid]]$Z
      g_res_G <- graph_res$G[[cid]]
      b_spec_G <- g_res_G$block_spec
      p_keys_G <- g_res_G$adapter_pair_keys
      K_G <- length(p_keys_G)
      g_key_G <- paste0(cid, "__G")

      G_links_c <- vector("list", K_G)
      for (k in seq_len(K_G)) {
        m_g <- obs_adapter$mean(
          group = g_key_G,
          pair = p_keys_G[k],
          h = ones_N,
          h_jacobian = list(),
          h_coef_if = NULL
        )
        G_links_c[[k]] <- static_missing_wrap_adapter_moment(m_g, unit_id)
      }

      zG_c <- static_reconstructed_or_role_mean_z(
        unit_id = unit_id,
        A = cl$G,
        Z = Z_or_c,
        materialize_if = FALSE
      )

      C_edges_arg <- NULL
      C_links_c <- NULL
      zCr_c <- NULL
      if (identical(m_info$canonical, "dr")) {
        g_res_C <- graph_res$C[[cid]]
        b_spec_C <- g_res_C$block_spec
        p_keys_C <- g_res_C$adapter_pair_keys
        K_C <- length(p_keys_C)
        g_key_C <- paste0(cid, "__C")
        C_edges_arg <- b_spec_C$edges[, c("s", "u"), drop = FALSE]
        C_links_c <- vector("list", K_C)
        for (k in seq_len(K_C)) {
          m_c <- obs_adapter$mean(
            group = g_key_C,
            pair = p_keys_C[k],
            h = trt_h_by_cell[[cid]],
            h_jacobian = trt_jac_by_cell[[cid]],
            h_coef_if = trt_coef_if_by_cell[[cid]]
          )
          C_links_c[[k]] <- static_missing_wrap_adapter_moment(m_c, unit_id)
        }
        zCr_c <- static_reconstructed_or_role_mean_z(
          unit_id = unit_id,
          A = cl$C,
          Z = Z_or_c,
          r = trt_h_by_cell[[cid]],
          X_tilde = trt_X_safe_by_cell[[cid]],
          active = trt_active_by_cell[[cid]],
          ell_gamma = logit_fits[[cid]]$ell_gamma,
          materialize_if = FALSE
        )
      }

      raw_links_by_cell[[cid]]$G_means <- G_links_c
      raw_links_by_cell[[cid]]$C_tilted_means <- C_links_c

      cc_obj <- static_reconstructed_or_center_cell_links(
        or_fit = or_fit_centering,
        cell_id = cid,
        method = m_info$gmm_method,
        t = cl$t_num,
        G_edges = b_spec_G$edges[, c("s", "u"), drop = FALSE],
        G_links = G_links_c,
        C_edges = C_edges_arg,
        C_links = C_links_c,
        zG = zG_c,
        zCr = zCr_c,
        materialize_if = FALSE
      )
      centered_by_cell[[cid]] <- static_missing_wrap_centered_cell_links(cc_obj, unit_id)
    }
  } else {
    # method == "ipw": uncentered G and tilted C links for every target cell in A
    for (cid in target_ids) {
      cl <- target_map[[cid]]
      g_res_G <- graph_res$G[[cid]]
      b_spec_G <- g_res_G$block_spec
      p_keys_G <- g_res_G$adapter_pair_keys
      K_G <- length(p_keys_G)
      g_key_G <- paste0(cid, "__G")

      G_links_c <- vector("list", K_G)
      for (k in seq_len(K_G)) {
        m_g <- obs_adapter$mean(
          group = g_key_G,
          pair = p_keys_G[k],
          h = ones_N,
          h_jacobian = list(),
          h_coef_if = NULL
        )
        G_links_c[[k]] <- static_missing_wrap_adapter_moment(m_g, unit_id)
      }

      g_res_C <- graph_res$C[[cid]]
      b_spec_C <- g_res_C$block_spec
      p_keys_C <- g_res_C$adapter_pair_keys
      K_C <- length(p_keys_C)
      g_key_C <- paste0(cid, "__C")

      C_links_c <- vector("list", K_C)
      for (k in seq_len(K_C)) {
        m_c <- obs_adapter$mean(
          group = g_key_C,
          pair = p_keys_C[k],
          h = trt_h_by_cell[[cid]],
          h_jacobian = trt_jac_by_cell[[cid]],
          h_coef_if = trt_coef_if_by_cell[[cid]]
        )
        C_links_c[[k]] <- static_missing_wrap_adapter_moment(m_c, unit_id)
      }

      raw_links_by_cell[[cid]] <- list(G_means = G_links_c, C_tilted_means = C_links_c)
      centered_by_cell[[cid]] <- static_missing_ipw_cell_links(
        cell_id = cid,
        unit_id = unit_id,
        t = cl$t_num,
        G_edges = b_spec_G$edges[, c("s", "u"), drop = FALSE],
        G_links = G_links_c,
        C_edges = b_spec_C$edges[, c("s", "u"), drop = FALSE],
        C_links = C_links_c
      )
    }
  }

  # 6. Assemble ALL target cells and role blocks into ONE joint static_graph_gmm solve
  gmm_blocks <- list()
  mhat_parts <- list()
  cell_col_ranges <- vector("list", length(target_ids))
  names(cell_col_ranges) <- target_ids
  col_cursor <- 0L

  for (cid in target_ids) {
    cc <- centered_by_cell[[cid]]
    b_spec_G <- graph_res$G[[cid]]$block_spec
    links_G <- stats::setNames(as.double(cc$adjusted_G), b_spec_G$edges$link_id)
    gmm_blocks[[length(gmm_blocks) + 1L]] <- b_spec_G
    mhat_parts[[length(mhat_parts) + 1L]] <- links_G

    if (m_info$gmm_method %in% c("ipw", "dr")) {
      b_spec_C <- graph_res$C[[cid]]$block_spec
      links_C <- stats::setNames(as.double(cc$adjusted_C), b_spec_C$edges$link_id)
      gmm_blocks[[length(gmm_blocks) + 1L]] <- b_spec_C
      mhat_parts[[length(mhat_parts) + 1L]] <- links_C
      n_c_links <- length(links_G) + length(links_C)
      c_link_ids <- c(names(links_G), names(links_C))
    } else {
      n_c_links <- length(links_G)
      c_link_ids <- names(links_G)
    }

    cell_col_ranges[[cid]] <- list(
      cols = col_cursor + seq_len(n_c_links),
      link_ids = c_link_ids
    )
    col_cursor <- col_cursor + n_c_links
  }

  mhat_vec <- unlist(mhat_parts, use.names = TRUE)
  all_link_ids <- names(mhat_vec)
  m_R <- length(all_link_ids)

  # The accepted GMM API validates explicit materialization requests against
  # both caller authorization and allocation limits, with typed refusal.
  mat_gmm_if <- materialize_if

  # Admit operation-local batch reuse of built-in selection/imputation IFs and
  # compact reconstructed-OR centering IFs only when the simultaneous
  # consumer + producer + holder reservations fit within `bounds$max_matrix_entries`.
  nuisance_batch_budget <- static_missing_admit_nuisance_batch_reuse(
    n = n,
    batch_rows = batch_rows,
    bounds = bounds,
    obs_adapter = obs_adapter,
    gmm_blocks = gmm_blocks,
    m_R = m_R,
    target_ids = target_ids,
    all_ids = all_ids,
    m_info = m_info,
    or_fit = or_fit,
    or_designs = or_designs,
    trt_X_safe_by_cell = trt_X_safe_by_cell,
    graph_res = graph_res,
    w_choice = w_choice,
    accuracy = accuracy,
    cluster = cluster,
    materialize_target_if = mat_gmm_if,
    aggregate_cells = aggregate_cells
  )
  reuse_nuisance_if <- isTRUE(nuisance_batch_budget$admitted)
  reuse_or_if <- isTRUE(nuisance_batch_budget$or_admitted)

  with_or_if_batch <- function(rows_pos, fn, enabled = TRUE) {
    if (!is.function(fn)) {
      static_missing_abort("`fn` must be a function.", "cdid_static_missing_input_error")
    }
    if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
      static_missing_abort("`enabled` must be a scalar logical.", "cdid_static_missing_control_error")
    }
    if (or_batch_active || !isTRUE(enabled) || is.null(or_fit) || length(rows_pos) > batch_rows) {
      return(fn())
    }
    or_batch_active <<- TRUE
    or_batch_pos <<- as.integer(rows_pos)
    or_batch_slice <<- NULL
    on.exit({
      or_batch_active <<- FALSE
      or_batch_pos <<- NULL
      or_batch_slice <<- NULL
    }, add = TRUE)
    fn()
  }

  # Row-batch IF reader for the stacked centered link vector psi^m
  link_if_reader <- function(rows = seq_len(n)) {
    idx <- static_missing_resolve_rows(rows, unit_id)
    m <- length(idx)
    obs_adapter$with_if_batch(
      rows = unit_id[idx],
      enabled = reuse_nuisance_if && (m <= batch_rows),
      fn = function() {
        with_or_if_batch(
          rows_pos = idx,
          enabled = reuse_or_if && (m <= batch_rows),
          fn = function() {
            out <- matrix(0.0, nrow = m, ncol = m_R, dimnames = list(unit_id_chr[idx], all_link_ids))
            for (cid in target_ids) {
              slice_c <- centered_by_cell[[cid]]$if_callback(idx)
              rng <- cell_col_ranges[[cid]]
              out[, rng$cols] <- slice_c$chi_stacked
              slice_c <- NULL
            }
            out
          }
        )
      }
    )
  }

  gmm_if_cb <- make_static_if_batch_callback(
    phi_source = link_if_reader,
    unit_ids = unit_id,
    batch_rows = batch_rows,
    link_ids = all_link_ids
  )

  gmm_res <- solve_static_graph_gmm(
    blocks = gmm_blocks,
    method = m_info$gmm_method,
    mhat = mhat_vec,
    if_callback = gmm_if_cb,
    unit_ids = unit_id,
    cluster_ids = cluster,
    weighting = w_choice$requested,
    bounds = bounds,
    accuracy = accuracy,
    cell_ids = target_ids,
    centered = centered,
    materialize_target_if = mat_gmm_if
  )

  # 7. Extract cell estimates, sign-oriented effects, covariances, and original-unit IF callbacks
  Lambda <- gmm_res$moment_weights$Lambda
  theta_hat <- gmm_res$estimates[target_ids]
  sigma_vec <- stats::setNames(
    vapply(target_map, function(cl) as.double(cl$sigma), numeric(1L)),
    target_ids
  )
  tau_hat <- stats::setNames(sigma_vec * theta_hat, target_ids)
  pi_hat <- stats::setNames(
    vapply(target_map, function(cl) as.double(cl$pi_hat), numeric(1L)),
    target_ids
  )
  q0_hat <- stats::setNames(
    vapply(target_map, function(cl) as.double(cl$q0_hat), numeric(1L)),
    target_ids
  )
  n_switchers <- stats::setNames(
    vapply(target_map, function(cl) as.integer(cl$n_switchers), integer(1L)),
    target_ids
  )
  n_stayers <- stats::setNames(
    vapply(target_map, function(cl) as.integer(cl$n_stayers), integer(1L)),
    target_ids
  )
  n_risk <- stats::setNames(
    vapply(target_map, function(cl) as.integer(cl$n_risk), integer(1L)),
    target_ids
  )
  Q_hat <- stats::setNames(rep(NA_real_, length(target_ids)), target_ids)
  if (m_info$use_logit) {
    for (cid in target_ids) Q_hat[[cid]] <- logit_fits[[cid]]$Q_hat
  }

  cov_theta <- gmm_res$vcov_n2[target_ids, target_ids, drop = FALSE]
  vhat_theta <- gmm_res$asymptotic_cov[target_ids, target_ids, drop = FALSE]
  se_theta <- stats::setNames(sqrt(pmax(0.0, diag(cov_theta))), target_ids)

  sign_diag <- diag(as.double(sigma_vec), nrow = length(target_ids), ncol = length(target_ids))
  cov_tau <- sign_diag %*% cov_theta %*% sign_diag
  vhat_tau <- sign_diag %*% vhat_theta %*% sign_diag
  dimnames(cov_tau) <- list(target_ids, target_ids)
  dimnames(vhat_tau) <- list(target_ids, target_ids)
  se_tau <- se_theta

  cell_if_callback <- function(rows = seq_len(n)) {
    idx <- static_missing_resolve_rows(rows, unit_id)
    chi_sub <- link_if_reader(idx)
    phi_th <- chi_sub %*% Lambda[, target_ids, drop = FALSE]
    rownames(phi_th) <- unit_id_chr[idx]
    colnames(phi_th) <- target_ids
    phi_ta <- sweep(phi_th, 2L, sigma_vec, FUN = "*")
    list(
      rows = idx,
      unit_id = unit_id[idx],
      phi_theta = phi_th,
      phi_tau = phi_ta
    )
  }

  phi_theta <- NULL
  phi_tau <- NULL
  if (mat_gmm_if && !is.null(gmm_res$target_if)) {
    phi_theta <- gmm_res$target_if[, target_ids, drop = FALSE]
    phi_tau <- sweep(phi_theta, 2L, sigma_vec, FUN = "*")
  }

  cell_summary <- data.frame(
    cell_id = target_ids,
    t = vapply(target_map, `[[`, numeric(1L), "t_num"),
    d = vapply(target_map, `[[`, integer(1L), "d"),
    sigma = as.double(sigma_vec),
    n_risk = as.integer(n_risk),
    n_switchers = as.integer(n_switchers),
    n_stayers = as.integer(n_stayers),
    pi_hat = as.double(pi_hat),
    q0_hat = as.double(q0_hat),
    Q_hat = as.double(Q_hat),
    theta_hat = as.double(theta_hat),
    se_theta = as.double(se_theta),
    tau_hat = as.double(tau_hat),
    se_tau = as.double(se_tau),
    row.names = target_ids,
    stringsAsFactors = FALSE
  )

  # Diagnostics using author-approved vocabulary
  zero_outcome_units <- sum(rowSums(panel_obj$O) == 0L)
  single_outcome_units <- sum(rowSums(panel_obj$O) == 1L)

  fit_obj <- structure(
    list(
      method = m_info$canonical,
      requested_method = m_info$requested,
      gmm_method = m_info$gmm_method,
      route = r_info$route,
      weighting = w_choice$canonical,
      requested_weighting = w_choice$requested,
      n = n,
      unit_id = unit_id,
      dates = panel_obj$dates,
      dates_num = panel_obj$dates_num,
      cell_ids = target_ids,
      target_cells = target_ids,
      all_cell_ids = all_ids,
      non_target_cell_ids = cell_res$non_target_ids,
      or_blocks = cell_res$or_blocks,
      or_partition = cell_res$or_partition,
      theta_hat = theta_hat,
      tau_hat = tau_hat,
      sigma = sigma_vec,
      pi_hat = pi_hat,
      q0_hat = q0_hat,
      Q_hat = Q_hat,
      n_risk = n_risk,
      n_switchers = n_switchers,
      n_stayers = n_stayers,
      cov_theta = cov_theta,
      vhat_theta = vhat_theta,
      se_theta = se_theta,
      cov_tau = cov_tau,
      vhat_tau = vhat_tau,
      se_tau = se_tau,
      phi_theta = phi_theta,
      phi_tau = phi_tau,
      cell_if_callback = cell_if_callback,
      link_if_callback = link_if_reader,
      mhat = mhat_vec,
      moment_weights = gmm_res$moment_weights,
      moment_covariance = gmm_res$moment_covariance,
      cell_summary = cell_summary,
      cells = target_map,
      all_cells = all_map,
      role_graphs = graph_res,
      centered_links = centered_by_cell,
      raw_links = raw_links_by_cell,
      nuisance = list(
        logit = if (m_info$use_logit) logit_fits else NULL,
        observation = obs_adapter,
        or = or_fit
      ),
      gmm = gmm_res,
      cluster = cluster,
      centered = centered,
      bounds = bounds,
      accuracy = accuracy,
      rank_tol = rank_tol,
      glm_control = glm_control,
      aggregate = NULL,
      diagnostics = list(
        method = m_info$canonical,
        requested_method = m_info$requested,
        route = r_info$route,
        selection_model = obs_adapter$selection_model,
        selection_pooling = obs_adapter$selection_pooling,
        imputation_pooling = obs_adapter$imputation_pooling,
        imputation_scale = obs_adapter$imputation_scale,
        weighting_requested = w_choice$requested,
        weighting_canonical = w_choice$canonical,
        n_units = n,
        n_zero_outcome_units = as.integer(zero_outcome_units),
        n_single_outcome_units = as.integer(single_outcome_units),
        n_target_cells = length(target_ids),
        n_all_cells = length(all_ids),
        n_non_target_cells = length(cell_res$non_target_ids),
        n_links = m_R,
        overid_df = gmm_res$diagnostics$overid_df,
        gmm_diagnostics = gmm_res$diagnostics
      )
    ),
    class = "cdid_static_missing_fit"
  )

  if (!is.null(aggregate_cells)) {
    fit_obj$aggregate <- static_missing_aggregate(
      fit = fit_obj,
      cells = aggregate_cells,
      signs = aggregate_signs,
      cluster = cluster,
      materialize_if = materialize_if
    )
  }

  fit_obj
}

# Switcher-weighted aggregate effect delta_{A_agg, s} and random switch-share IF
# for a fitted `cdid_static_missing_fit` object. Does NOT forge `cdid_static_complete_fit`.
static_missing_aggregate <- function(fit,
                                     cells = fit$cell_ids,
                                     signs = NULL,
                                     cluster = fit$cluster,
                                     materialize_if = !is.null(fit$phi_theta)) {
  if (!inherits(fit, "cdid_static_missing_fit")) {
    static_missing_abort(
      "`fit` must be a `cdid_static_missing_fit` object from `static_missing_fit()`.",
      "cdid_static_missing_input_error"
    )
  }

  cell_ids <- fit$cell_ids
  if (missing(cells) || is.null(cells) || length(cells) < 1L) {
    static_missing_abort(
      "`cells` must be a non-empty subset of `fit$cell_ids`.",
      "cdid_static_missing_aggregate_error"
    )
  }
  if (is.numeric(cells)) {
    if (any(!is.finite(cells)) || any(cells != as.integer(cells)) ||
        any(cells < 1L | cells > length(cell_ids))) {
      static_missing_abort(
        "Integer `cells` indices are out of range for `fit$cell_ids`.",
        "cdid_static_missing_aggregate_error"
      )
    }
    agg_ids <- cell_ids[as.integer(cells)]
  } else if (is.character(cells)) {
    if (anyNA(cells) || any(!nzchar(cells))) {
      static_missing_abort(
        "`cells` must not contain NA or empty strings.",
        "cdid_static_missing_aggregate_error"
      )
    }
    unknown <- setdiff(cells, cell_ids)
    if (length(unknown) > 0L) {
      static_missing_abort(
        sprintf("`cells` contains unknown target cell IDs: %s.", paste(unknown, collapse = ", ")),
        "cdid_static_missing_aggregate_error"
      )
    }
    agg_ids <- cells
  } else {
    static_missing_abort(
      "`cells` must be a character vector of target cell IDs or integer indices.",
      "cdid_static_missing_aggregate_error"
    )
  }
  if (anyDuplicated(agg_ids)) {
    static_missing_abort(
      "`cells` must not contain duplicate cell IDs.",
      "cdid_static_missing_aggregate_error"
    )
  }

  k_agg <- length(agg_ids)
  if (is.null(signs)) {
    s_vec <- stats::setNames(rep(1.0, k_agg), agg_ids)
  } else {
    if (!is.numeric(signs) || is.object(signs) || !is.null(dim(signs))) {
      static_missing_abort(
        "`signs` must be a numeric vector with entries in `{-1, +1}`.",
        "cdid_static_missing_aggregate_error"
      )
    }
    if (!is.null(names(signs))) {
      if (anyDuplicated(names(signs)) || !all(agg_ids %in% names(signs))) {
        static_missing_abort(
          "Named `signs` must provide a unique entry for every cell in `cells`.",
          "cdid_static_missing_aggregate_error"
        )
      }
      s_vec <- as.double(signs[agg_ids])
    } else {
      if (length(signs) != k_agg) {
        static_missing_abort(
          sprintf("`signs` must have length %d matching `cells`.", k_agg),
          "cdid_static_missing_aggregate_error"
        )
      }
      s_vec <- as.double(signs)
    }
    if (any(!is.finite(s_vec)) || any(!(s_vec %in% c(-1, 1)))) {
      static_missing_abort(
        "`signs` must take values strictly in `{-1, +1}`.",
        "cdid_static_missing_aggregate_error"
      )
    }
    names(s_vec) <- agg_ids
  }

  n <- fit$n
  uid_chr <- as.character(fit$unit_id)
  if (!gmm_scalar_logical(materialize_if)) {
    static_missing_abort("`materialize_if` must be a scalar logical.", "cdid_gmm_control_error")
  }
  if (materialize_if) {
    if (!isTRUE(fit$bounds$allow_materialize_if)) {
      static_missing_abort("Aggregate IF materialization is not authorized by bounds.", "cdid_gmm_bound_error")
    }
    gmm_safe_product(n, 1L, fit$bounds$max_materialized_if_entries, "Materialized aggregate IF")
    gmm_safe_product(n, 1L, fit$bounds$max_matrix_entries, "Materialized aggregate IF")
  }
  n_c <- fit$n_switchers[agg_ids]
  pi_c <- fit$pi_hat[agg_ids]
  total_switchers <- sum(n_c)
  Pi_hat <- sum(pi_c)
  if (total_switchers < 1L || !is.finite(Pi_hat) || Pi_hat <= 0) {
    static_missing_abort(
      "Total switcher mass `Pi_hat` across `cells` is zero.",
      c("cdid_static_support_error", "cdid_static_missing_aggregate_error")
    )
  }

  omega_hat <- stats::setNames(as.double(n_c) / as.double(total_switchers), agg_ids)
  theta_sub <- fit$theta_hat[agg_ids]
  signed_theta <- stats::setNames(s_vec * theta_sub, agg_ids)
  delta_hat <- sum(omega_hat * signed_theta)
  weight_phi <- omega_hat * s_vec
  dev_theta <- signed_theta - delta_hat

  if_callback <- function(rows = seq_len(n)) {
    idx <- static_missing_resolve_rows(rows, fit$unit_id)
    phi_th_sub <- if (!is.null(fit$phi_theta)) {
      fit$phi_theta[idx, agg_ids, drop = FALSE]
    } else {
      fit$cell_if_callback(idx)$phi_theta[, agg_ids, drop = FALSE]
    }
    cell_term <- drop(phi_th_sub %*% weight_phi)
    share_term <- numeric(length(idx))
    for (j in seq_len(k_agg)) {
      share_term <- share_term + (as.double(fit$cells[[agg_ids[j]]]$G[idx]) - pi_c[j]) * dev_theta[j]
    }
    share_term <- share_term / Pi_hat
    phi_r <- stats::setNames(cell_term + share_term, uid_chr[idx])
    list(
      rows = idx,
      unit_id = fit$unit_id[idx],
      phi = phi_r,
      cell_if_term = stats::setNames(cell_term, uid_chr[idx]),
      share_if_term = stats::setNames(share_term, uid_chr[idx])
    )
  }

  # Evaluate full or chunked aggregate IF for covariance
  batch_rows <- as.integer(min(n, fit$bounds$max_batch_rows))
  phi_full <- numeric(n)
  cell_term_full <- if (materialize_if) numeric(n) else NULL
  share_term_full <- if (materialize_if) numeric(n) else NULL
  starts <- seq.int(1L, n, by = batch_rows)
  for (s_i in starts) {
    idx <- s_i:min(n, s_i + batch_rows - 1L)
    part <- if_callback(idx)
    phi_full[idx] <- part$phi
    if (materialize_if) {
      cell_term_full[idx] <- part$cell_if_term
      share_term_full[idx] <- part$share_if_term
    }
  }
  names(phi_full) <- uid_chr
  if (materialize_if) {
    names(cell_term_full) <- uid_chr
    names(share_term_full) <- uid_chr
  }

  cov_agg <- static_complete_covariance(phi_full, cluster = cluster)

  structure(
    list(
      cells = agg_ids,
      signs = s_vec,
      delta_hat = delta_hat,
      omega_hat = omega_hat,
      pi_hat = pi_c,
      Pi_hat = Pi_hat,
      total_switchers = as.integer(total_switchers),
      signed_theta = signed_theta,
      phi = if (isTRUE(materialize_if)) phi_full else NULL,
      phi_components = if (isTRUE(materialize_if)) {
        list(cell_if_term = cell_term_full, share_if_term = share_term_full)
      } else {
        NULL
      },
      if_callback = if_callback,
      variance = as.double(cov_agg$cov),
      vhat = as.double(cov_agg$vhat),
      se = as.double(cov_agg$se),
      clustered = cov_agg$clustered,
      n_clusters = cov_agg$n_clusters
    ),
    class = "cdid_static_missing_aggregate"
  )
}
