# ==============================================================================
# Private Static Role-Graph / Joint-GMM Numerical Adapter
# ==============================================================================
#
# Layer: Private numerical adapter ONLY for the binary static missing-outcome
# extension adapting the static binary specialization of de Chaisemartin,
# D'Haultfoeuille, Pasquier, Sow and Vazquez-Bare (2026, arXiv:2201.06898v7)
# with the chained GMM framework of Bellego, Benatia and Dortet-Bernadet (2025).
#
# Scope:
#   - Accepts caller-explicit fixed role graphs with unique block/cell/role IDs,
#     strictly ordered dates, baseline t-1 and target t, explicit retained edges
#     s < u, aligned centered point moments `mhat`, and a complete original-unit
#     influence-function (IF) bounded-batch replay callback over all N units
#     (including zero/single-outcome units).
#   - Does NOT fit nuisances, define working models, select supported edges or
#     cells from data, prune graphs, drop observations/moments, or introduce
#     switch-step outcome centering.
#   - Builds a single joint `cdid_chained_design` incidence representation
#     retaining all valid edge rows across K_c^a >= 1 connected components per
#     role block, grounded at each role block's baseline t-1 in the endpoint
#     component T_{c,1}^a and one deterministic coordinate gauge per additional
#     component (rank = |vertices| - K); all other date levels in the block
#     (including earlier dates s < t-1) remain unrestricted parameters.
#   - Constructs one target per cell: G(t) - C(t) for IPW/DR, and G(t) for OR.
#   - Streams the full cross-cell / cross-role m_R x m_R moment covariance V
#     from the complete original-unit IF stream in bounded batches across two
#     passes, and seals `provenance$complete_unit_if_gram` via `ipw_seal()` ONLY
#     when V was genuinely formed here from the complete unit IF stream.
#   - Solves the joint GMM via `solve_chained_gmm_core()` with existing identity
#     and two-step + Tikhonov covariance handling.
#   - Clusters affect inference only (empirical cluster-sum sandwich covariance);
#     point estimates, unit moment covariance V, and unit-GMM weights/loadings
#     remain unchanged.
# ==============================================================================

cdid_static_graph_abort <- function(message, subclass = NULL, diagnostics = list()) {
  cdid_gmm_abort(
    message = message,
    class = unique(c(subclass, "cdid_static_graph_error")),
    diagnostics = if (is.null(diagnostics)) list() else diagnostics
  )
}

static_graph_normalize_role <- function(role, where = "role") {
  if (!gmm_scalar_character(role)) {
    cdid_static_graph_abort(
      paste0("`", where, "` must be a non-empty character scalar."),
      "cdid_static_graph_input_error"
    )
  }
  if (role %in% c("G", "switcher")) {
    return(list(code = "G", name = "switcher"))
  }
  if (role %in% c("C", "stayer")) {
    return(list(code = "C", name = "stayer"))
  }
  cdid_static_graph_abort(
    paste0("Invalid role `", role, "` in `", where, "`; must be \"G\"/\"switcher\" or \"C\"/\"stayer\"."),
    c("cdid_static_graph_role_mismatch_error", "cdid_static_graph_input_error")
  )
}

static_graph_validate_method <- function(method) {
  if (missing(method) || !gmm_scalar_character(method)) {
    cdid_static_graph_abort(
      "`method` must be supplied explicitly as one of \"ipw\", \"dr\", or \"or\".",
      "cdid_static_graph_input_error"
    )
  }
  m <- tolower(method)
  if (!m %in% c("ipw", "dr", "or")) {
    cdid_static_graph_abort(
      paste0("Unsupported `method = \"", method, "\"`; must be \"ipw\", \"dr\", or \"or\"."),
      "cdid_static_graph_input_error"
    )
  }
  m
}

static_graph_validate_ids <- function(ids, label, allow_duplicates = FALSE) {
  if (is.null(ids) || !is.null(dim(ids)) || is.object(ids) || length(ids) < 1L ||
      anyNA(ids) || !(is.character(ids) || is.numeric(ids))) {
    cdid_static_graph_abort(
      paste0("`", label, "` must be a non-empty atomic character or numeric vector without NAs."),
      "cdid_static_graph_input_error"
    )
  }
  if (is.numeric(ids)) {
    if (any(!is.finite(ids))) {
      cdid_static_graph_abort(
        paste0("Numeric `", label, "` must be finite."),
        "cdid_static_graph_input_error"
      )
    }
    chr <- sprintf("%.17g", as.double(ids))
  } else {
    if (any(!nzchar(ids))) {
      cdid_static_graph_abort(
        paste0("Character `", label, "` must be non-empty strings."),
        "cdid_static_graph_input_error"
      )
    }
    chr <- as.character(ids)
  }
  if (!allow_duplicates && anyDuplicated(chr)) {
    cdid_static_graph_abort(
      paste0("Duplicate identifiers in `", label, "`."),
      "cdid_static_graph_input_error"
    )
  }
  chr
}

# Convenience constructor for the `conventional_chol_qr` + Tikhonov
# GMM accuracy control on the static missing-outcome route.
static_graph_gmm_accuracy_control <- function(
  method = "conventional_chol_qr",
  target_error_unit = "abs",
  max_target_error_abs = 1e-6,
  max_target_error_se = 1e-3,
  max_orthogonality_defect = 1e-10,
  min_certified_lambda_rel = 1e-12,
  max_factor_residual_rel = 1e-10,
  min_rcond_R = 1e-12,
  min_rcond_U = 1e-12,
  min_sigma_rel_A = 1e-10,
  max_kkt_residual_rel = 1e-8,
  max_loading_residual_rel = 1e-8,
  require_matching_scale = FALSE,
  max_loading_error_abs = NULL,
  covariance_strategy = gmm_covariance_strategy(
    method = "tikhonov",
    condition_trigger = 1e-10,
    ridge_relative = 1e-4,
    weighting_requested = "optimal",
    weighting_canonical = "optimal"
  )
) {
  gmm_accuracy_control(
    method = method,
    target_error_unit = target_error_unit,
    max_target_error_abs = max_target_error_abs,
    max_target_error_se = max_target_error_se,
    max_orthogonality_defect = max_orthogonality_defect,
    min_certified_lambda_rel = min_certified_lambda_rel,
    max_factor_residual_rel = max_factor_residual_rel,
    min_rcond_R = min_rcond_R,
    min_rcond_U = min_rcond_U,
    min_sigma_rel_A = min_sigma_rel_A,
    max_kkt_residual_rel = max_kkt_residual_rel,
    max_loading_residual_rel = max_loading_residual_rel,
    require_matching_scale = require_matching_scale,
    max_loading_error_abs = max_loading_error_abs,
    covariance_strategy = covariance_strategy
  )
}

# ------------------------------------------------------------------------------
# 1. Role-Block Specification & Validation
# ------------------------------------------------------------------------------

static_role_block_spec <- function(
  block_id,
  cell_id,
  role,
  baseline_time,
  target_time,
  dates,
  edges,
  mhat = NULL
) {
  if (missing(block_id) || !gmm_scalar_character(block_id)) {
    cdid_static_graph_abort("`block_id` must be a non-empty character scalar.", "cdid_static_graph_input_error")
  }
  if (missing(cell_id) || !gmm_scalar_character(cell_id)) {
    cdid_static_graph_abort("`cell_id` must be a non-empty character scalar.", "cdid_static_graph_input_error")
  }
  if (missing(role)) {
    cdid_static_graph_abort("`role` must be supplied explicitly.", "cdid_static_graph_input_error")
  }
  role_info <- static_graph_normalize_role(role, where = paste0("block `", block_id, "` role"))

  if (missing(baseline_time) || !gmm_scalar_numeric(baseline_time, positive = FALSE)) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `baseline_time` must be a finite numeric scalar."),
      "cdid_static_graph_input_error"
    )
  }
  if (missing(target_time) || !gmm_scalar_numeric(target_time, positive = FALSE)) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `target_time` must be a finite numeric scalar."),
      "cdid_static_graph_input_error"
    )
  }
  t_base <- as.double(baseline_time)
  t_targ <- as.double(target_time)
  if (!(t_base < t_targ)) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `baseline_time` (", fmt_design_num(t_base),
             ") must be strictly less than `target_time` (", fmt_design_num(t_targ), ")."),
      c("cdid_static_graph_date_order_error", "cdid_static_graph_input_error")
    )
  }

  if (missing(dates) || !is.atomic(dates) || !is.null(dim(dates)) || is.object(dates) ||
      !is.numeric(dates) || length(dates) < 2L || anyNA(dates) || any(!is.finite(dates))) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `dates` must be a finite numeric vector of length >= 2."),
      "cdid_static_graph_input_error"
    )
  }
  dates_num <- as.double(dates)
  if (anyDuplicated(dates_num) || is.unsorted(dates_num, strictly = TRUE)) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `dates` must be strictly increasing without duplicates."),
      c("cdid_static_graph_date_order_error", "cdid_static_graph_input_error")
    )
  }

  if (!t_base %in% dates_num) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `baseline_time` (", fmt_design_num(t_base), ") is missing from `dates`."),
      "cdid_static_graph_missing_baseline_error"
    )
  }
  if (!t_targ %in% dates_num) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `target_time` (", fmt_design_num(t_targ), ") is missing from `dates`."),
      "cdid_static_graph_missing_target_error"
    )
  }

  if (missing(edges) || !is.data.frame(edges) || nrow(edges) < 1L || anyDuplicated(names(edges))) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `edges` must be a non-empty data.frame with unique column names."),
      "cdid_static_graph_input_error"
    )
  }
  s_col <- if ("s" %in% names(edges)) "s" else if ("start_time" %in% names(edges)) "start_time" else NULL
  u_col <- if ("u" %in% names(edges)) "u" else if ("end_time" %in% names(edges)) "end_time" else if ("t" %in% names(edges)) "t" else NULL
  if (is.null(s_col) || is.null(u_col)) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: `edges` must contain endpoint columns `s` and `u`."),
      "cdid_static_graph_input_error"
    )
  }
  s_raw <- edges[[s_col]]
  u_raw <- edges[[u_col]]
  if (!is.atomic(s_raw) || !is.null(dim(s_raw)) || is.object(s_raw) || !is.numeric(s_raw) ||
      !is.atomic(u_raw) || !is.null(dim(u_raw)) || is.object(u_raw) || !is.numeric(u_raw) ||
      anyNA(s_raw) || anyNA(u_raw) || any(!is.finite(s_raw)) || any(!is.finite(u_raw))) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: edge endpoints `s` and `u` must be finite numeric vectors."),
      "cdid_static_graph_input_error"
    )
  }
  s_num <- as.double(s_raw)
  u_num <- as.double(u_raw)
  n_edges <- length(s_num)
  if (any(!(s_num < u_num))) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: every retained edge must satisfy s < u."),
      c("cdid_static_graph_edge_order_error", "cdid_static_graph_input_error")
    )
  }

  from_idx <- match(s_num, dates_num)
  to_idx <- match(u_num, dates_num)
  if (anyNA(from_idx) || anyNA(to_idx)) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: all edge endpoints (s, u) must belong to declared `dates`."),
      c("cdid_static_graph_missing_vertex_error", "cdid_static_graph_input_error")
    )
  }
  pair_keys <- paste0(from_idx, ":", to_idx)
  if (anyDuplicated(pair_keys)) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: duplicate (s, u) edges within a single role block are forbidden."),
      c("cdid_static_graph_duplicate_edge_error", "cdid_static_graph_input_error")
    )
  }

  id_col <- if ("link_id" %in% names(edges)) "link_id" else if ("edge_id" %in% names(edges)) "edge_id" else NULL
  if (!is.null(id_col)) {
    link_ids <- static_graph_validate_ids(edges[[id_col]], paste0("block `", block_id, "` edge IDs"), allow_duplicates = FALSE)
  } else {
    link_ids <- paste0(
      block_id, ":",
      vapply(s_num, fmt_design_num, character(1L)), "->",
      vapply(u_num, fmt_design_num, character(1L))
    )
  }

  # Validate endpoint connectedness and incident support of required dates t-1 and t.
  # Additional connected components (K_c^a > 1) and harmless isolated non-target dates
  # s not in {t-1, t} are retained and receive one deterministic coordinate gauge per component.
  base_v_idx <- match(t_base, dates_num)
  targ_v_idx <- match(t_targ, dates_num)
  forest <- incidence_forest(
    nv = length(dates_num),
    from = from_idx,
    to = to_idx,
    anchors = base_v_idx,
    edge_order = seq_len(n_edges)
  )
  if (forest$degree[base_v_idx] == 0L) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: baseline t-1 (", fmt_design_num(t_base), ") has no incident edges."),
      c("cdid_static_graph_missing_baseline_error", "cdid_static_graph_disconnected_error")
    )
  }
  if (forest$degree[targ_v_idx] == 0L) {
    cdid_static_graph_abort(
      paste0("Block `", block_id, "`: target t (", fmt_design_num(t_targ), ") has no incident edges."),
      c("cdid_static_graph_missing_target_error", "cdid_static_graph_disconnected_error")
    )
  }
  if (forest$component_id[targ_v_idx] != forest$component_id[base_v_idx]) {
    cdid_static_graph_abort(
      paste0(
        "Block `", block_id, "`: baseline t-1 (", fmt_design_num(t_base),
        ") and target t (", fmt_design_num(t_targ),
        ") lie in distinct connected components of role graph (T_c^a, K_c^a) (`as:graph` requires t-1 and t to belong to the same connected component T_{c,1}^a)."
      ),
      "cdid_static_graph_disconnected_error",
      diagnostics = list(
        block_id = block_id,
        cell_id = cell_id,
        role = role_info$code,
        n_components = forest$n_components,
        baseline_component_id = forest$component_id[base_v_idx],
        target_component_id = forest$component_id[targ_v_idx]
      )
    )
  }

  # Optional block-level mhat alignment
  edge_mhat <- if ("mhat" %in% names(edges)) edges$mhat else NULL
  if (!is.null(edge_mhat) && !is.null(mhat)) {
    if (!is.numeric(edge_mhat) || !is.numeric(mhat) || length(edge_mhat) != n_edges ||
        length(mhat) != n_edges || !identical(as.double(edge_mhat), as.double(mhat))) {
      cdid_static_graph_abort(
        paste0("Block `", block_id, "`: `edges$mhat` and block `mhat` are both supplied and disagree."),
        "cdid_static_graph_alignment_error"
      )
    }
  }
  mhat_use <- if (!is.null(mhat)) mhat else edge_mhat
  if (!is.null(mhat_use)) {
    if (!is.atomic(mhat_use) || !is.null(dim(mhat_use)) || is.object(mhat_use) ||
        !is.numeric(mhat_use) || length(mhat_use) != n_edges || anyNA(mhat_use) || any(!is.finite(mhat_use))) {
      cdid_static_graph_abort(
        paste0("Block `", block_id, "`: block `mhat` must be a finite numeric vector of length `nrow(edges)` (", n_edges, ")."),
        "cdid_static_graph_alignment_error"
      )
    }
    if (!is.null(names(mhat_use)) && !identical(as.character(names(mhat_use)), link_ids)) {
      cdid_static_graph_abort(
        paste0("Block `", block_id, "`: named block `mhat` does not match edge `link_id` order."),
        "cdid_static_graph_alignment_error"
      )
    }
    mhat_use <- as.double(mhat_use)
  }

  edges_clean <- data.frame(
    link_id = link_ids,
    s = s_num,
    u = u_num,
    from_date_index = as.integer(from_idx),
    to_date_index = as.integer(to_idx),
    stringsAsFactors = FALSE
  )
  if (!is.null(mhat_use)) {
    edges_clean$mhat <- mhat_use
  }

  structure(
    list(
      block_id = as.character(block_id),
      cell_id = as.character(cell_id),
      role = role_info$code,
      role_name = role_info$name,
      baseline_time = t_base,
      target_time = t_targ,
      dates = dates_num,
      baseline_date_index = as.integer(base_v_idx),
      target_date_index = as.integer(targ_v_idx),
      edges = edges_clean,
      mhat = mhat_use,
      forest = forest,
      n_components = as.integer(forest$n_components),
      endpoint_component_id = as.integer(forest$component_id[base_v_idx]),
      reference_date_indices = as.integer(forest$roots),
      reference_dates = as.double(dates_num[forest$roots]),
      isolated_dates = as.double(dates_num[forest$degree == 0L])
    ),
    class = "cdid_static_role_block_spec"
  )
}

validate_static_role_block_spec <- function(block) {
  if (inherits(block, "cdid_static_role_block_spec") && is.list(block)) {
    return(static_role_block_spec(
      block_id = block$block_id,
      cell_id = block$cell_id,
      role = block$role,
      baseline_time = block$baseline_time,
      target_time = block$target_time,
      dates = block$dates,
      edges = block$edges,
      mhat = block$mhat
    ))
  }
  if (!is.list(block)) {
    cdid_static_graph_abort(
      "Each element of `blocks` must be a `cdid_static_role_block_spec` or a named list.",
      "cdid_static_graph_input_error"
    )
  }
  req <- c("block_id", "cell_id", "role", "baseline_time", "target_time", "dates", "edges")
  if (!all(req %in% names(block))) {
    cdid_static_graph_abort(
      paste0("Role block list is missing required fields: ", paste(setdiff(req, names(block)), collapse = ", "), "."),
      "cdid_static_graph_input_error"
    )
  }
  static_role_block_spec(
    block_id = block$block_id,
    cell_id = block$cell_id,
    role = block$role,
    baseline_time = block$baseline_time,
    target_time = block$target_time,
    dates = block$dates,
    edges = block$edges,
    mhat = block$mhat
  )
}

# ------------------------------------------------------------------------------
# 2. Joint `cdid_chained_design` Incidence & Target Construction
# ------------------------------------------------------------------------------

build_static_graph_gmm_design <- function(
  blocks,
  method,
  mhat = NULL,
  cell_ids = NULL
) {
  method <- static_graph_validate_method(method)
  if (missing(blocks) || !is.list(blocks) || is.object(blocks) || length(blocks) < 1L) {
    cdid_static_graph_abort(
      "`blocks` must be a non-empty plain list of role-block specifications.",
      "cdid_static_graph_input_error"
    )
  }

  n_blocks <- length(blocks)
  v_blocks <- lapply(blocks, validate_static_role_block_spec)

  block_ids <- vapply(v_blocks, `[[`, character(1L), "block_id")
  if (anyDuplicated(block_ids)) {
    cdid_static_graph_abort(
      "Duplicate `block_id` across role blocks.",
      c("cdid_static_graph_duplicate_block_error", "cdid_static_graph_input_error")
    )
  }

  block_cells <- vapply(v_blocks, `[[`, character(1L), "cell_id")
  block_roles <- vapply(v_blocks, `[[`, character(1L), "role")
  cell_role_keys <- paste0(block_cells, "::", block_roles)
  if (anyDuplicated(cell_role_keys)) {
    cdid_static_graph_abort(
      "Duplicate `(cell_id, role)` block supplied.",
      c("cdid_static_graph_duplicate_role_error", "cdid_static_graph_input_error")
    )
  }

  unique_cells <- unique(block_cells)
  if (!is.null(cell_ids)) {
    cell_ids <- static_graph_validate_ids(cell_ids, "cell_ids", allow_duplicates = FALSE)
    if (!setequal(cell_ids, unique_cells)) {
      cdid_static_graph_abort(
        "`cell_ids` must contain the exact set of `cell_id` values present in `blocks`.",
        "cdid_static_graph_input_error"
      )
    }
  } else {
    cell_ids <- unique_cells
  }

  n_cells <- length(cell_ids)
  expected_roles <- if (method == "or") "G" else c("G", "C")
  cell_base <- numeric(n_cells)
  cell_targ <- numeric(n_cells)
  cell_G_block <- integer(n_cells)
  cell_C_block <- rep(NA_integer_, n_cells)

  for (ci in seq_len(n_cells)) {
    cid <- cell_ids[ci]
    idx <- which(block_cells == cid)
    roles_c <- block_roles[idx]
    if (!setequal(roles_c, expected_roles) || length(roles_c) != length(expected_roles)) {
      cdid_static_graph_abort(
        paste0(
          "Cell `", cid, "` has roles {", paste(roles_c, collapse = ", "),
          "}, but method `", method, "` requires exactly {", paste(expected_roles, collapse = ", "), "}."
        ),
        "cdid_static_graph_role_mismatch_error",
        diagnostics = list(cell_id = cid, method = method, supplied_roles = roles_c, expected_roles = expected_roles)
      )
    }
    bases_c <- vapply(v_blocks[idx], `[[`, numeric(1L), "baseline_time")
    targs_c <- vapply(v_blocks[idx], `[[`, numeric(1L), "target_time")
    if (length(unique(bases_c)) != 1L || length(unique(targs_c)) != 1L) {
      cdid_static_graph_abort(
        paste0("Cell `", cid, "` role blocks disagree on `baseline_time` or `target_time`."),
        "cdid_static_graph_date_order_error"
      )
    }
    cell_base[ci] <- bases_c[1L]
    cell_targ[ci] <- targs_c[1L]
    cell_G_block[ci] <- idx[roles_c == "G"]
    if (method != "or") {
      cell_C_block[ci] <- idx[roles_c == "C"]
    }
  }

  # Assemble joint vertices, parameters, links, and incidence grounded at each
  # block's baseline t-1.
  vert_rows <- vector("list", n_blocks)
  param_rows <- vector("list", n_blocks)
  link_rows <- vector("list", n_blocks)
  block_meta_rows <- vector("list", n_blocks)
  anchor_vertices <- integer(n_blocks)
  block_target_param <- integer(n_blocks)

  v_offset <- 0L
  p_offset <- 0L
  l_offset <- 0L

  for (b in seq_len(n_blocks)) {
    blk <- v_blocks[[b]]
    nv_b <- length(blk$dates)
    np_b <- nv_b - 1L
    nl_b <- nrow(blk$edges)

    v_indices <- v_offset + seq_len(nv_b)
    is_base <- seq_len(nv_b) == blk$baseline_date_index
    is_targ <- seq_len(nv_b) == blk$target_date_index

    p_local <- rep(NA_integer_, nv_b)
    p_local[!is_base] <- p_offset + seq_len(np_b)

    v_ids <- paste0(
      "v(block=", blk$block_id,
      ",cell=", blk$cell_id,
      ",role=", blk$role,
      ",t=", vapply(blk$dates, fmt_design_num, character(1L)),
      ")"
    )
    p_ids <- paste0(
      "nu(block=", blk$block_id,
      ",cell=", blk$cell_id,
      ",role=", blk$role,
      ",t=", vapply(blk$dates[!is_base], fmt_design_num, character(1L)),
      ")"
    )
    v_role_label <- ifelse(
      is_base,
      "reference",
      ifelse(is_targ, "target", ifelse(blk$dates < blk$baseline_time, "earlier_free", "intermediate_free"))
    )

    vert_rows[[b]] <- data.frame(
      vertex_index = as.integer(v_indices),
      vertex_id = v_ids,
      block_index = rep(as.integer(b), nv_b),
      block_id = rep(blk$block_id, nv_b),
      cell_id = rep(blk$cell_id, nv_b),
      g = rep(as.double(match(blk$cell_id, cell_ids)), nv_b),
      time = blk$dates,
      baseline_time = rep(blk$baseline_time, nv_b),
      target_time = rep(blk$target_time, nv_b),
      role = v_role_label,
      role_code = rep(blk$role, nv_b),
      role_name = rep(blk$role_name, nv_b),
      functional_key = rep(paste0("static_block:", blk$block_id), nv_b),
      component_id = rep(as.integer(b), nv_b),
      is_anchor = is_base,
      is_target_date = is_targ,
      anchored_value = ifelse(is_base, 0.0, NA_real_),
      is_isolated = rep(FALSE, nv_b),
      connected_to_ref = rep(TRUE, nv_b),
      param_index = as.integer(p_local),
      stringsAsFactors = FALSE
    )

    non_base_idx <- which(!is_base)
    param_rows[[b]] <- data.frame(
      param_index = as.integer(p_local[non_base_idx]),
      param_id = p_ids,
      vertex_index = as.integer(v_indices[non_base_idx]),
      block_index = rep(as.integer(b), np_b),
      block_id = rep(blk$block_id, np_b),
      cell_id = rep(blk$cell_id, np_b),
      g = rep(as.double(match(blk$cell_id, cell_ids)), np_b),
      time = blk$dates[non_base_idx],
      baseline_time = rep(blk$baseline_time, np_b),
      target_time = rep(blk$target_time, np_b),
      role = v_role_label[non_base_idx],
      role_code = rep(blk$role, np_b),
      role_name = rep(blk$role_name, np_b),
      functional_key = rep(paste0("static_block:", blk$block_id), np_b),
      component_id = rep(as.integer(b), np_b),
      is_target_date = is_targ[non_base_idx],
      is_isolated = rep(FALSE, np_b),
      level_estimable = rep(TRUE, np_b),
      stringsAsFactors = FALSE
    )

    sv_global <- v_indices[blk$edges$from_date_index]
    ev_global <- v_indices[blk$edges$to_date_index]
    sp_global <- p_local[blk$edges$from_date_index]
    ep_global <- p_local[blk$edges$to_date_index]

    link_rows[[b]] <- data.frame(
      link_index = as.integer(l_offset + seq_len(nl_b)),
      link_id = blk$edges$link_id,
      block_index = rep(as.integer(b), nl_b),
      block_id = rep(blk$block_id, nl_b),
      cell_id = rep(blk$cell_id, nl_b),
      role_code = rep(blk$role, nl_b),
      role_name = rep(blk$role_name, nl_b),
      g = rep(as.double(match(blk$cell_id, cell_ids)), nl_b),
      s = blk$edges$s,
      u = blk$edges$u,
      t = blk$edges$u,
      ref = rep(blk$baseline_time, nl_b),
      baseline_time = rep(blk$baseline_time, nl_b),
      target_time = rep(blk$target_time, nl_b),
      functional_key = rep(paste0("static_block:", blk$block_id), nl_b),
      start_vertex_index = as.integer(sv_global),
      end_vertex_index = as.integer(ev_global),
      start_param_index = as.integer(sp_global),
      end_param_index = as.integer(ep_global),
      component_id = rep(as.integer(b), nl_b),
      stringsAsFactors = FALSE
    )

    anchor_v <- v_indices[blk$baseline_date_index]
    target_v <- v_indices[blk$target_date_index]
    target_p <- p_local[blk$target_date_index]
    anchor_vertices[b] <- as.integer(anchor_v)
    block_target_param[b] <- as.integer(target_p)

    block_meta_rows[[b]] <- data.frame(
      block_index = as.integer(b),
      block_id = blk$block_id,
      cell_id = blk$cell_id,
      role = blk$role,
      role_name = blk$role_name,
      baseline_time = blk$baseline_time,
      target_time = blk$target_time,
      n_dates = as.integer(nv_b),
      n_parameters = as.integer(np_b),
      n_edges = as.integer(nl_b),
      component_id = as.integer(b),
      anchor_vertex_index = as.integer(anchor_v),
      target_vertex_index = as.integer(target_v),
      target_param_index = as.integer(target_p),
      stringsAsFactors = FALSE
    )

    v_offset <- v_offset + nv_b
    p_offset <- p_offset + np_b
    l_offset <- l_offset + nl_b
  }

  vertices_df <- do.call(rbind, vert_rows)
  parameters_df <- do.call(rbind, param_rows)
  links_df <- do.call(rbind, link_rows)
  blocks_df <- do.call(rbind, block_meta_rows)
  rownames(vertices_df) <- rownames(parameters_df) <- rownames(links_df) <- rownames(blocks_df) <- NULL

  if (anyDuplicated(links_df$link_id)) {
    cdid_static_graph_abort(
      "Duplicate `link_id` across the joint role-graph system.",
      c("cdid_static_graph_duplicate_link_error", "cdid_static_graph_input_error")
    )
  }

  n_links_total <- nrow(links_df)
  n_vertices_total <- nrow(vertices_df)
  n_params_total <- nrow(parameters_df)

  joint_forest <- incidence_forest(
    nv = n_vertices_total,
    from = links_df$start_vertex_index,
    to = links_df$end_vertex_index,
    anchors = anchor_vertices,
    edge_order = seq_len(n_links_total)
  )
  anchored_comps <- as.integer(joint_forest$component_id[anchor_vertices])
  target_comps <- as.integer(joint_forest$component_id[blocks_df$target_vertex_index])
  if (length(unique(anchored_comps)) != n_blocks || !identical(anchored_comps, target_comps)) {
    cdid_static_graph_abort(
      "Joint role-graph forest failed component/anchor verification.",
      "cdid_static_graph_disconnected_error"
    )
  }
  unanchored_comps <- as.integer(setdiff(seq_len(joint_forest$n_components), anchored_comps))

  # Align joint forest component_id, isolation, and level estimability across tables
  vertices_df$component_id <- as.integer(joint_forest$component_id)
  vertices_df$is_isolated <- as.logical(joint_forest$degree == 0L)
  vertices_df$connected_to_ref <- as.logical(vertices_df$component_id %in% anchored_comps)

  parameters_df$component_id <- as.integer(joint_forest$component_id[parameters_df$vertex_index])
  parameters_df$is_isolated <- as.logical(joint_forest$degree[parameters_df$vertex_index] == 0L)
  parameters_df$level_estimable <- as.logical(parameters_df$component_id %in% anchored_comps)

  links_df$component_id <- as.integer(joint_forest$component_id[links_df$start_vertex_index])
  blocks_df$component_id <- as.integer(joint_forest$component_id[blocks_df$anchor_vertex_index])

  inc_start <- ifelse(is.na(links_df$start_param_index), 0L, links_df$start_param_index)
  inc_end <- ifelse(is.na(links_df$end_param_index), 0L, links_df$end_param_index)

  incidence <- list(
    start = as.integer(inc_start),
    end = as.integer(inc_end),
    n_links = as.integer(n_links_total),
    n_parameters = as.integer(n_params_total),
    link_ids = links_df$link_id,
    param_ids = parameters_df$param_id
  )

  n_unanchored <- length(unanchored_comps)
  rank_val <- as.integer(n_params_total - n_unanchored)
  rank_diagnostics <- list(
    n_links = as.integer(n_links_total),
    n_vertices = as.integer(n_vertices_total),
    n_parameters = as.integer(n_params_total),
    n_components = as.integer(joint_forest$n_components),
    rank = rank_val,
    nullity = as.integer(n_unanchored),
    cycle_rank = as.integer(n_links_total - rank_val),
    has_full_column_rank = n_unanchored == 0L,
    unanchored_components = unanchored_comps
  )

  # Construct one target per cell: G(target) - C(target) for IPW/DR, G(target) for OR
  target_specs <- vector("list", n_cells)
  cell_rows <- vector("list", n_cells)
  for (ci in seq_len(n_cells)) {
    cid <- cell_ids[ci]
    b_g <- cell_G_block[ci]
    p_g <- block_target_param[b_g]
    blk_g <- v_blocks[[b_g]]
    tid <- cid

    if (method == "or") {
      target_specs[[ci]] <- gmm_target_spec(
        target_id = tid,
        param_index = p_g,
        coeff = 1.0
      )
      cell_rows[[ci]] <- data.frame(
        cell_index = as.integer(ci),
        cell_id = cid,
        target_id = tid,
        method = method,
        baseline_time = cell_base[ci],
        target_time = cell_targ[ci],
        switcher_block_id = blk_g$block_id,
        stayer_block_id = NA_character_,
        switcher_target_param_index = as.integer(p_g),
        stayer_target_param_index = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      b_c <- cell_C_block[ci]
      p_c <- block_target_param[b_c]
      blk_c <- v_blocks[[b_c]]
      target_specs[[ci]] <- gmm_target_spec(
        target_id = tid,
        param_index = c(p_g, p_c),
        coeff = c(1.0, -1.0)
      )
      cell_rows[[ci]] <- data.frame(
        cell_index = as.integer(ci),
        cell_id = cid,
        target_id = tid,
        method = method,
        baseline_time = cell_base[ci],
        target_time = cell_targ[ci],
        switcher_block_id = blk_g$block_id,
        stayer_block_id = blk_c$block_id,
        switcher_target_param_index = as.integer(p_g),
        stayer_target_param_index = as.integer(p_c),
        stringsAsFactors = FALSE
      )
    }
  }
  cells_df <- do.call(rbind, cell_rows)
  rownames(cells_df) <- NULL

  # Resolve and validate aligned point moments mhat
  has_block_mhat <- vapply(v_blocks, function(x) !is.null(x$mhat), logical(1L))
  if (any(has_block_mhat) && !all(has_block_mhat) && is.null(mhat)) {
    cdid_static_graph_abort(
      "Block-level `mhat` was supplied for only a subset of blocks while joint `mhat` is NULL.",
      "cdid_static_graph_alignment_error"
    )
  }
  block_mhat_vec <- if (all(has_block_mhat)) {
    unlist(lapply(v_blocks, `[[`, "mhat"), use.names = FALSE)
  } else {
    NULL
  }
  if (!is.null(mhat)) {
    if (!is.atomic(mhat) || !is.null(dim(mhat)) || is.object(mhat) || !is.numeric(mhat) ||
        length(mhat) != n_links_total || anyNA(mhat) || any(!is.finite(mhat))) {
      cdid_static_graph_abort(
        paste0("`mhat` must be a finite numeric vector of length `m_R` (", n_links_total, ")."),
        "cdid_static_graph_alignment_error"
      )
    }
    if (!is.null(names(mhat)) && !identical(as.character(names(mhat)), links_df$link_id)) {
      cdid_static_graph_abort(
        "Named `mhat` does not match the joint `link_id` order.",
        "cdid_static_graph_alignment_error"
      )
    }
    mhat_vec <- as.double(mhat)
    if (!is.null(block_mhat_vec) && !identical(mhat_vec, as.double(block_mhat_vec))) {
      cdid_static_graph_abort(
        "Joint `mhat` and block-level `mhat` are both supplied and disagree.",
        "cdid_static_graph_alignment_error"
      )
    }
  } else if (!is.null(block_mhat_vec)) {
    mhat_vec <- as.double(block_mhat_vec)
  } else {
    mhat_vec <- NULL
  }
  if (!is.null(mhat_vec)) {
    names(mhat_vec) <- links_df$link_id
    links_df$mhat <- as.double(mhat_vec)
  }

  design <- structure(
    list(
      links = links_df,
      vertices = vertices_df,
      parameters = parameters_df,
      incidence = incidence,
      forest = joint_forest,
      rank_diagnostics = rank_diagnostics,
      static_graph_spec = list(
        method = method,
        blocks = blocks_df,
        cells = cells_df,
        target_specs = target_specs,
        role_blocks = v_blocks,
        mhat = mhat_vec
      )
    ),
    class = c("cdid_static_graph_design", "cdid_chained_design")
  )
  design
}

materialize_static_graph_incidence_matrix <- function(
  design,
  max_matrix_entries = NULL,
  drop_component_gauges = FALSE,
  gauge_rule = "lowest",
  custom_gauge_roots = NULL
) {
  if (!inherits(design, "cdid_chained_design") || !is.list(design$incidence)) {
    cdid_static_graph_abort("`design` must be a valid `cdid_chained_design` object.", "cdid_static_graph_input_error")
  }
  if (!gmm_scalar_logical(drop_component_gauges)) {
    cdid_static_graph_abort("`drop_component_gauges` must be a scalar logical.", "cdid_static_graph_input_error")
  }
  if (isTRUE(drop_component_gauges)) {
    gauge <- build_gmm_gauge(design, gauge_rule = gauge_rule, custom_roots = custom_gauge_roots)
    m_e <- if (is.null(max_matrix_entries)) max(1L, gauge$n_links * max(1L, gauge$rank_k)) else max_matrix_entries
    return(materialize_gmm_gauge_matrix(gauge, max_matrix_entries = m_e))
  }
  inc <- design$incidence
  m <- inc$n_links
  p <- inc$n_parameters
  if (!is.null(max_matrix_entries)) {
    gmm_safe_product(m, p, gmm_int_bound(max_matrix_entries, "max_matrix_entries"), "dense incidence matrix D")
  }
  D <- matrix(0.0, nrow = m, ncol = p, dimnames = list(inc$link_ids, inc$param_ids))
  for (e in seq_len(m)) {
    s_idx <- inc$start[e]
    u_idx <- inc$end[e]
    if (s_idx > 0L) D[e, s_idx] <- D[e, s_idx] - 1.0
    if (u_idx > 0L) D[e, u_idx] <- D[e, u_idx] + 1.0
  }
  D
}

# ------------------------------------------------------------------------------
# 3. Bounded-Batch Complete-Unit IF Callback Helper & Moment Covariance Stream
# ------------------------------------------------------------------------------

make_static_if_batch_callback <- function(
  phi_source,
  unit_ids,
  batch_rows,
  link_ids = NULL,
  nuisance_fit_id = NULL
) {
  u_chr <- static_graph_validate_ids(unit_ids, "unit_ids", allow_duplicates = FALSE)
  n_units <- length(u_chr)
  b_rows <- gmm_int_bound(batch_rows, "batch_rows")
  if (b_rows > .Machine$integer.max) {
    cdid_static_graph_abort("`batch_rows` exceeds integer index range.", "cdid_static_graph_input_error")
  }
  b_rows <- as.integer(b_rows)

  if (is.matrix(phi_source)) {
    if (!is.numeric(phi_source) || is.object(phi_source) || nrow(phi_source) != n_units ||
        ncol(phi_source) < 1L || anyNA(phi_source) || any(!is.finite(phi_source))) {
      cdid_static_graph_abort(
        "`phi_source` matrix must be a finite numeric matrix with `nrow(phi_source) == length(unit_ids)`.",
        "cdid_static_graph_alignment_error"
      )
    }
    if (!is.null(rownames(phi_source)) && !identical(as.character(rownames(phi_source)), u_chr)) {
      cdid_static_graph_abort(
        "`rownames(phi_source)` do not match `unit_ids`.",
        "cdid_static_graph_alignment_error"
      )
    }
    if (!is.null(link_ids) && !is.null(colnames(phi_source)) &&
        !identical(as.character(colnames(phi_source)), as.character(link_ids))) {
      cdid_static_graph_abort(
        "`colnames(phi_source)` do not match `link_ids`.",
        "cdid_static_graph_alignment_error"
      )
    }
    reader_fn <- function(idx) {
      sub <- phi_source[idx, , drop = FALSE]
      if (!is.null(link_ids) && is.null(colnames(sub))) {
        colnames(sub) <- link_ids
      }
      sub
    }
  } else if (is.function(phi_source)) {
    reader_fn <- phi_source
  } else {
    cdid_static_graph_abort(
      "`phi_source` must be a numeric matrix or a batch reader function `function(idx)`.",
      "cdid_static_graph_input_error"
    )
  }

  function(consume) {
    starts <- seq.int(1L, n_units, by = b_rows)
    for (st in starts) {
      en <- min(n_units, st + b_rows - 1L)
      idx <- seq.int(st, en)
      phi_b <- reader_fn(idx)
      batch <- list(
        unit_id = unit_ids[idx],
        phi = phi_b
      )
      if (!is.null(nuisance_fit_id)) {
        batch$nuisance_fit_id <- nuisance_fit_id
      }
      consume(batch)
    }
    invisible(NULL)
  }
}

# One deterministic replay contract spans covariance, unit inference and
# cluster inference. Reuse the existing bounded row-sum/squared-row checks;
# these diagnose inconsistent callbacks, not equality against arbitrary
# colliding streams. The caller must still provide deterministic replay.
static_graph_cross_stage_replay <- function(if_callback, unit_ids, link_ids, bounds) {
  if (!is.function(if_callback)) {
    cdid_static_graph_abort("Require a deterministic IF replay callback.", "cdid_gmm_replay_error")
  }
  N <- length(unit_ids)
  M <- length(link_ids)
  gmm_safe_product(3, N, bounds$max_matrix_entries, "cross-stage replay verification vectors")
  batch_sizes <- integer(N)
  row_sum <- row_square <- numeric(N)
  passes <- 0L
  n_batches <- 0L
  function(consume) {
    first <- passes == 0L
    nin <- nb <- 0L
    if_callback(function(batch) {
      # Let the stage's existing validator preserve its precise diagnostics.
      # A later consistency error aborts the whole solve before any result.
      consume(batch)
      if (!is.list(batch) || anyDuplicated(names(batch)) || !("phi" %in% names(batch))) {
        cdid_static_graph_abort("Invalid cross-stage IF batch.", "cdid_gmm_replay_error")
      }
      phi <- batch$phi
      if (!is.matrix(phi) || !is.numeric(phi) || is.object(phi) ||
          nrow(phi) < 1L || ncol(phi) != M || nrow(phi) > bounds$max_batch_rows ||
          any(!is.finite(phi))) {
        cdid_static_graph_abort("Invalid cross-stage IF batch values or dimensions.", "cdid_gmm_replay_error")
      }
      B <- nrow(phi)
      if (nin + B > N) {
        cdid_static_graph_abort("Cross-stage IF replay exceeds its original-unit count.", "cdid_gmm_replay_error")
      }
      raw_ids <- if ("unit_id" %in% names(batch)) batch$unit_id else batch$id
      ids <- static_graph_validate_ids(raw_ids, "batch unit_id", allow_duplicates = TRUE)
      idx <- seq.int(nin + 1L, nin + B)
      rs <- unname(rowSums(phi))
      rq <- unname(rowSums(phi * phi))
      if (!identical(ids, unit_ids[idx]) || any(!is.finite(rs)) || any(!is.finite(rq))) {
        cdid_static_graph_abort("Cross-stage IF replay unit order or row checks are invalid.", "cdid_gmm_replay_error")
      }
      nb <<- nb + 1L
      if (first) {
        batch_sizes[nb] <<- B
        row_sum[idx] <<- rs
        row_square[idx] <<- rq
      } else if (nb > n_batches || B != batch_sizes[nb] ||
                 !identical(rs, row_sum[idx]) || !identical(rq, row_square[idx])) {
        cdid_static_graph_abort("IF replay changed between covariance and inference stages.", "cdid_gmm_replay_error")
      }
      nin <<- nin + B
    })
    if (nin != N || (!first && nb != n_batches)) {
      cdid_static_graph_abort("Cross-stage IF replay has inconsistent row or batch counts.", "cdid_gmm_replay_error")
    }
    if (first) n_batches <<- nb
    passes <<- passes + 1L
    invisible(NULL)
  }
}

stream_static_graph_moment_covariance <- function(
  if_callback,
  unit_ids,
  link_ids,
  bounds,
  centered = TRUE,
  normalization_factor = NULL,
  nuisance_fit_id = "nuisance_fit:static_caller_supplied"
) {
  if (!is.function(if_callback)) {
    cdid_static_graph_abort("`if_callback` must be a replay function `function(consume)`.", "cdid_static_graph_input_error")
  }
  bounds <- validate_gmm_bounds_control(bounds)
  u_expected <- static_graph_validate_ids(unit_ids, "unit_ids", allow_duplicates = FALSE)
  l_expected <- static_graph_validate_ids(link_ids, "link_ids", allow_duplicates = FALSE)
  if (!gmm_scalar_logical(centered)) {
    cdid_static_graph_abort("`centered` must be a scalar logical.", "cdid_static_graph_input_error")
  }
  if (!gmm_scalar_character(nuisance_fit_id)) {
    cdid_static_graph_abort("`nuisance_fit_id` must be a non-empty character scalar.", "cdid_static_graph_input_error")
  }

  N <- length(u_expected)
  m_R <- length(l_expected)
  if (N > .Machine$integer.max || m_R > .Machine$integer.max) {
    cdid_static_graph_abort("Unit or moment dimension exceeds R integer range.", "cdid_gmm_bound_error")
  }
  if (is.null(normalization_factor)) {
    normalization_factor <- 1.0 / as.double(N)
  }
  if (!gmm_scalar_numeric(normalization_factor, positive = TRUE)) {
    cdid_static_graph_abort("`normalization_factor` must be a finite positive numeric scalar.", "cdid_static_graph_input_error")
  }
  norm_fac <- as.double(normalization_factor)

  max_e <- bounds$max_matrix_entries
  B_max <- min(N, bounds$max_batch_rows)
  gmm_safe_product(m_R, m_R, max_e, "static moment covariance accumulator S (m_R x m_R)")
  gmm_safe_product(2, gmm_safe_product(B_max, m_R, max_e, "static IF batch (B x m_R)"), max_e, "static IF batch and centered batch")
  gmm_safe_product(4, N, max_e, "static IF stream unit verification vectors")

  batch_sizes <- integer(N)
  row_sum1 <- numeric(N)
  row_sq1 <- numeric(N)
  raw_sum <- numeric(m_R)
  nb1 <- 0L
  max_B_seen <- 0L

  extract_batch <- function(batch) {
    if (!is.list(batch) || anyDuplicated(names(batch)) || !("phi" %in% names(batch))) {
      cdid_static_graph_abort("Invalid IF batch: must be a list with field `phi` and `unit_id` (or `id`).", "cdid_gmm_replay_error")
    }
    uid_raw <- if ("unit_id" %in% names(batch)) batch$unit_id else if ("id" %in% names(batch)) batch$id else NULL
    if (is.null(uid_raw)) {
      cdid_static_graph_abort("IF batch is missing `unit_id`.", "cdid_gmm_replay_error")
    }
    phi <- batch$phi
    if (!is.matrix(phi) || !is.numeric(phi) || is.object(phi) || nrow(phi) < 1L) {
      cdid_static_graph_abort("IF batch `phi` must be a non-empty numeric matrix.", "cdid_gmm_replay_error")
    }
    if (ncol(phi) != m_R) {
      cdid_static_graph_abort(
        paste0("IF batch `phi` column count (", ncol(phi), ") does not match joint moment count m_R (", m_R, ")."),
        c("cdid_static_graph_alignment_error", "cdid_gmm_replay_error")
      )
    }
    if (!is.null(colnames(phi)) && !identical(as.character(colnames(phi)), l_expected)) {
      cdid_static_graph_abort(
        "IF batch `colnames(phi)` do not match the joint `link_id` order.",
        c("cdid_static_graph_alignment_error", "cdid_gmm_replay_error")
      )
    }
    B <- nrow(phi)
    if (B > bounds$max_batch_rows) {
      cdid_static_graph_abort(
        paste0("IF batch row count (", B, ") exceeds `bounds$max_batch_rows` (", bounds$max_batch_rows, ")."),
        c("cdid_gmm_bound_error", "cdid_gmm_replay_error")
      )
    }
    gmm_safe_product(B, m_R, max_e, "IF batch validation")
    if (anyNA(phi) || any(!is.finite(phi))) {
      cdid_static_graph_abort("IF batch `phi` contains NA or non-finite values.", "cdid_gmm_replay_error")
    }
    if (!is.null(batch$nuisance_fit_id) && !identical(batch$nuisance_fit_id, nuisance_fit_id)) {
      cdid_static_graph_abort("IF batch `nuisance_fit_id` disagrees with declared nuisance identity.", "cdid_gmm_replay_error")
    }
    if (!is.null(batch$weight)) {
      w <- batch$weight
      if (!is.numeric(w) || !is.null(dim(w)) || is.object(w) || length(w) != B || anyNA(w) || any(w != 1.0)) {
        cdid_static_graph_abort("Survey weights are not supported in this slice; batch weights must be omitted or all 1.", "cdid_gmm_replay_error")
      }
    }
    u <- static_graph_validate_ids(uid_raw, "batch unit_id", allow_duplicates = TRUE)
    if (length(u) != B) {
      cdid_static_graph_abort("IF batch `unit_id` length does not match `nrow(phi)`.", "cdid_gmm_replay_error")
    }
    phi_num <- as.double(phi)
    dim(phi_num) <- c(B, m_R)
    list(u = u, phi = phi_num, B = B)
  }

  # Pass 1: validate batches, unit ID alignment, and accumulate sample mean
  nin1 <- 0L
  consume_pass1 <- function(batch) {
    eb <- extract_batch(batch)
    B <- eb$B
    if (nin1 + B > N) {
      cdid_static_graph_abort("IF stream emitted more rows than `length(unit_ids)`.", "cdid_gmm_replay_error")
    }
    idx <- seq.int(nin1 + 1L, nin1 + B)
    if (!identical(eb$u, u_expected[idx])) {
      cdid_static_graph_abort(
        "Streamed `unit_id` sequence does not match declared `unit_ids`.",
        c("cdid_static_graph_alignment_error", "cdid_gmm_replay_error")
      )
    }
    nb1 <<- nb1 + 1L
    batch_sizes[nb1] <<- B
    max_B_seen <<- max(max_B_seen, B)
    rs <- rowSums(eb$phi)
    rq <- rowSums(eb$phi * eb$phi)
    cs <- colSums(eb$phi)
    if (any(!is.finite(rs)) || any(!is.finite(rq)) || any(!is.finite(cs))) {
      cdid_static_graph_abort("IF stream row/column sum overflowed in Pass 1.", "cdid_gmm_replay_error")
    }
    row_sum1[idx] <<- rs
    row_sq1[idx] <<- rq
    raw_sum <<- raw_sum + cs
    if (any(!is.finite(raw_sum))) {
      cdid_static_graph_abort("IF stream mean accumulation overflowed in Pass 1.", "cdid_gmm_replay_error")
    }
    nin1 <<- nin1 + B
  }

  tryCatch(
    withCallingHandlers(
      if_callback(consume_pass1),
      warning = function(w) cdid_static_graph_abort(paste0("IF stream warning in Pass 1: ", conditionMessage(w)), "cdid_gmm_replay_error")
    ),
    error = function(e) {
      if (inherits(e, "cdid_error")) stop(e)
      cdid_static_graph_abort(paste0("IF stream failed in Pass 1: ", conditionMessage(e)), "cdid_gmm_replay_error")
    }
  )
  if (nin1 != N || nb1 < 1L) {
    cdid_static_graph_abort(
      paste0("IF stream Pass 1 emitted ", nin1, " rows across ", nb1, " batches; expected ", N, " rows."),
      "cdid_gmm_replay_error"
    )
  }

  mean_phi <- if (centered) raw_sum / as.double(N) else numeric(m_R)
  if (any(!is.finite(mean_phi))) {
    cdid_static_graph_abort("IF stream sample mean is non-finite.", "cdid_gmm_replay_error")
  }

  # Pass 2: verify identical batch structure and accumulate outer product S
  S <- matrix(0.0, nrow = m_R, ncol = m_R)
  nin2 <- 0L
  nb2 <- 0L
  consume_pass2 <- function(batch) {
    eb <- extract_batch(batch)
    B <- eb$B
    nb2 <<- nb2 + 1L
    if (nb2 > nb1 || B != batch_sizes[nb2] || nin2 + B > N) {
      cdid_static_graph_abort("IF stream batch structure changed between Pass 1 and Pass 2.", "cdid_gmm_replay_error")
    }
    idx <- seq.int(nin2 + 1L, nin2 + B)
    rs2 <- rowSums(eb$phi)
    rq2 <- rowSums(eb$phi * eb$phi)
    if (!identical(eb$u, u_expected[idx]) ||
        !identical(rs2, row_sum1[idx]) ||
        !identical(rq2, row_sq1[idx])) {
      cdid_static_graph_abort("IF stream unit IDs or row values changed between Pass 1 and Pass 2.", "cdid_gmm_replay_error")
    }
    phi_c <- if (centered) sweep(eb$phi, 2L, mean_phi, "-") else eb$phi
    if (any(!is.finite(phi_c))) {
      cdid_static_graph_abort("Centered IF batch overflowed in Pass 2.", "cdid_gmm_replay_error")
    }
    S <<- S + crossprod(phi_c)
    if (any(!is.finite(S))) {
      cdid_static_graph_abort("Moment covariance accumulation overflowed in Pass 2.", "cdid_gmm_replay_error")
    }
    nin2 <<- nin2 + B
  }

  tryCatch(
    withCallingHandlers(
      if_callback(consume_pass2),
      warning = function(w) cdid_static_graph_abort(paste0("IF stream warning in Pass 2: ", conditionMessage(w)), "cdid_gmm_replay_error")
    ),
    error = function(e) {
      if (inherits(e, "cdid_error")) stop(e)
      cdid_static_graph_abort(paste0("IF stream failed in Pass 2: ", conditionMessage(e)), "cdid_gmm_replay_error")
    }
  )
  if (nin2 != N || nb2 != nb1) {
    cdid_static_graph_abort("IF stream Pass 2 row or batch count disagreed with Pass 1.", "cdid_gmm_replay_error")
  }

  V <- norm_fac * S
  if (anyNA(V) || any(!is.finite(V)) || any(diag(V) < 0.0)) {
    cdid_static_graph_abort("Streamed moment covariance `V` contains non-finite values or negative diagonal entries.", "cdid_gmm_covariance_error")
  }
  for (j in seq_len(m_R)) {
    for (i in seq_len(j)) {
      V[j, i] <- V[i, j]
    }
  }
  dim(V) <- c(m_R, m_R)
  dimnames(V) <- NULL

  gram_seal <- ipw_seal(
    list(
      covariance = V,
      n_units = as.integer(N),
      nuisance_fit_id = nuisance_fit_id
    ),
    "complete_unit_if_gram"
  )

  list(
    V = V,
    mean_phi = stats::setNames(mean_phi, l_expected),
    raw_sum_phi = stats::setNames(raw_sum, l_expected),
    n_units = as.integer(N),
    n_links = as.integer(m_R),
    n_batches = as.integer(nb1),
    max_batch_rows_seen = as.integer(max_B_seen),
    centered = centered,
    normalization_factor = norm_fac,
    nuisance_fit_id = nuisance_fit_id,
    complete_unit_if_gram = gram_seal
  )
}

# Bounded two-pass cluster-sum sandwich accumulation on the replay stream.
# Does not require units to be pre-sorted by cluster and preserves unit-level
# point estimates, moment covariance V, and unit-GMM weights Lambda.
stream_static_graph_cluster_sandwich <- function(
  if_callback,
  Lambda,
  unit_ids,
  cluster_ids,
  link_ids,
  bounds,
  centered = TRUE,
  normalization_factor,
  nuisance_fit_id,
  materialize_cluster_if = FALSE
) {
  u_expected <- static_graph_validate_ids(unit_ids, "unit_ids", allow_duplicates = FALSE)
  l_expected <- static_graph_validate_ids(link_ids, "link_ids", allow_duplicates = FALSE)
  c_chr <- static_graph_validate_ids(cluster_ids, "cluster_ids", allow_duplicates = TRUE)
  N <- length(u_expected)
  if (length(c_chr) != N) {
    cdid_static_graph_abort(
      paste0("`cluster_ids` length (", length(c_chr), ") must equal `length(unit_ids)` (", N, ")."),
      "cdid_static_graph_input_error"
    )
  }
  u_clusters <- unique(c_chr)
  G <- length(u_clusters)
  if (G < 2L) {
    cdid_static_graph_abort(
      "At least 2 distinct clusters are required when `cluster_ids` is supplied.",
      "cdid_static_graph_input_error"
    )
  }
  c_map <- match(c_chr, u_clusters)
  m_R <- nrow(Lambda)
  J <- ncol(Lambda)
  max_e <- bounds$max_matrix_entries

  gmm_safe_product(G, J, max_e, "cluster IF sum accumulator (G x J)")
  gmm_safe_product(J, J, max_e, "cluster sandwich covariance (J x J)")
  if (isTRUE(materialize_cluster_if)) {
    if (!isTRUE(bounds$allow_materialize_if)) {
      cdid_static_graph_abort("IF materialization is not authorized by `bounds$allow_materialize_if`.", "cdid_gmm_bound_error")
    }
    gmm_safe_product(G, J, bounds$max_materialized_if_entries, "materialized cluster target IF (G x J)")
  }

  S_cluster1 <- matrix(0.0, nrow = G, ncol = J)
  batch_sizes <- integer(N)
  nb1 <- 0L
  nin1 <- 0L

  consume_c1 <- function(batch) {
    uid_raw <- if ("unit_id" %in% names(batch)) batch$unit_id else batch$id
    u <- static_graph_validate_ids(uid_raw, "batch unit_id", allow_duplicates = TRUE)
    phi <- batch$phi
    B <- nrow(phi)
    if (B > bounds$max_batch_rows || nin1 + B > N || ncol(phi) != m_R) {
      cdid_static_graph_abort("Invalid batch in cluster Pass 1.", "cdid_gmm_replay_error")
    }
    if (!is.null(colnames(phi)) && !identical(as.character(colnames(phi)), l_expected)) {
      cdid_static_graph_abort("IF batch `colnames(phi)` do not match `link_ids` in cluster Pass 1.", "cdid_static_graph_alignment_error")
    }
    idx <- seq.int(nin1 + 1L, nin1 + B)
    if (!identical(u, u_expected[idx])) {
      cdid_static_graph_abort("Batch `unit_id` mismatch in cluster Pass 1.", "cdid_gmm_replay_error")
    }
    if (!is.null(batch$cluster_id)) {
      b_cid <- static_graph_validate_ids(batch$cluster_id, "batch cluster_id", allow_duplicates = TRUE)
      if (!identical(b_cid, c_chr[idx])) {
        cdid_static_graph_abort("Batch `cluster_id` disagrees with declared `cluster_ids`.", "cdid_gmm_replay_error")
      }
    }
    psi_b <- phi %*% Lambda
    if (any(!is.finite(psi_b))) {
      cdid_static_graph_abort("Projected IF batch overflowed in cluster Pass 1.", "cdid_gmm_replay_error")
    }
    for (i in seq_len(B)) {
      g_i <- c_map[nin1 + i]
      S_cluster1[g_i, ] <<- S_cluster1[g_i, ] + psi_b[i, ]
    }
    if (any(!is.finite(S_cluster1))) {
      cdid_static_graph_abort("Cluster IF sum overflowed in Pass 1.", "cdid_gmm_replay_error")
    }
    nb1 <<- nb1 + 1L
    batch_sizes[nb1] <<- B
    nin1 <<- nin1 + B
  }

  if_callback(consume_c1)
  if (nin1 != N || nb1 < 1L) {
    cdid_static_graph_abort("Cluster Pass 1 row count mismatch.", "cdid_gmm_replay_error")
  }

  # Pass 2: verify identical cluster sums across replay
  S_cluster2 <- matrix(0.0, nrow = G, ncol = J)
  nb2 <- 0L
  nin2 <- 0L
  consume_c2 <- function(batch) {
    uid_raw <- if ("unit_id" %in% names(batch)) batch$unit_id else batch$id
    u <- static_graph_validate_ids(uid_raw, "batch unit_id", allow_duplicates = TRUE)
    phi <- batch$phi
    B <- nrow(phi)
    nb2 <<- nb2 + 1L
    if (nb2 > nb1 || B != batch_sizes[nb2] || nin2 + B > N) {
      cdid_static_graph_abort("Cluster Pass 2 batch structure changed.", "cdid_gmm_replay_error")
    }
    idx <- seq.int(nin2 + 1L, nin2 + B)
    if (!identical(u, u_expected[idx])) {
      cdid_static_graph_abort("Cluster Pass 2 unit IDs changed.", "cdid_gmm_replay_error")
    }
    psi_b <- phi %*% Lambda
    for (i in seq_len(B)) {
      g_i <- c_map[nin2 + i]
      S_cluster2[g_i, ] <<- S_cluster2[g_i, ] + psi_b[i, ]
    }
    nin2 <<- nin2 + B
  }
  if_callback(consume_c2)
  if (nin2 != N || nb2 != nb1 || !identical(S_cluster1, S_cluster2)) {
    cdid_static_graph_abort("Cluster IF sums changed between Pass 1 and Pass 2.", "cdid_gmm_replay_error")
  }

  cluster_mean <- if (centered) colMeans(S_cluster1) else numeric(J)
  S_centered <- if (centered) sweep(S_cluster1, 2L, cluster_mean, "-") else S_cluster1
  cov_raw <- crossprod(S_centered)
  sandwich <- normalization_factor * cov_raw
  if (anyNA(sandwich) || any(!is.finite(sandwich)) || any(diag(sandwich) < 0.0)) {
    cdid_static_graph_abort("Cluster sandwich covariance is non-finite or has negative diagonal.", "cdid_gmm_replay_error")
  }
  for (j in seq_len(J)) {
    for (i in seq_len(j)) {
      sandwich[j, i] <- sandwich[i, j]
    }
  }
  dimnames(sandwich) <- list(colnames(Lambda), colnames(Lambda))
  rownames(S_cluster1) <- u_clusters
  colnames(S_cluster1) <- colnames(Lambda)

  list(
    sandwich_cov = sandwich,
    se_sandwich = sqrt(diag(sandwich)),
    cluster_mean = stats::setNames(cluster_mean, colnames(Lambda)),
    n_clusters = as.integer(G),
    cluster_ids_unique = u_clusters,
    cluster_target_if = if (isTRUE(materialize_cluster_if)) S_cluster1 else NULL
  )
}

# ------------------------------------------------------------------------------
# 4. Main Private Static Role-Graph / Joint-GMM Adapter Entrypoint
# ------------------------------------------------------------------------------

solve_static_graph_gmm <- function(
  blocks,
  method,
  mhat = NULL,
  if_callback,
  unit_ids,
  cluster_ids,
  weighting,
  bounds,
  accuracy,
  cell_ids = NULL,
  V = NULL,
  provenance = NULL,
  centered = TRUE,
  normalization_factor = NULL,
  scale_id = "scale:static_complete_unit_if",
  nuisance_fit_id = "nuisance_fit:static_caller_supplied",
  request_param_cov = FALSE,
  materialize_target_if = FALSE,
  gauge_rule = "lowest",
  custom_gauge_roots = NULL
) {
  if (missing(blocks) || missing(method) || missing(if_callback) || missing(unit_ids) ||
      missing(cluster_ids) || missing(weighting) || missing(bounds) || missing(accuracy)) {
    cdid_static_graph_abort(
      "All adapter arguments (`blocks`, `method`, `if_callback`, `unit_ids`, `cluster_ids`, `weighting`, `bounds`, `accuracy`) must be supplied explicitly.",
      "cdid_static_graph_input_error"
    )
  }
  if (!is.function(if_callback)) {
    cdid_static_graph_abort("`if_callback` must be a replay function `function(consume)`.", "cdid_static_graph_input_error")
  }
  if (!gmm_scalar_logical(centered) || !gmm_scalar_logical(request_param_cov) || !gmm_scalar_logical(materialize_target_if)) {
    cdid_static_graph_abort("`centered`, `request_param_cov`, and `materialize_target_if` must be scalar logicals.", "cdid_static_graph_input_error")
  }
  if (!gmm_scalar_character(scale_id) || !gmm_scalar_character(nuisance_fit_id)) {
    cdid_static_graph_abort("`scale_id` and `nuisance_fit_id` must be non-empty character scalars.", "cdid_static_graph_input_error")
  }

  # 1. Validate weighting, bounds, accuracy, and unit/cluster IDs
  w_choice <- resolve_gmm_weighting_choice(weighting, arg = "weighting")
  bounds <- validate_gmm_bounds_control(bounds)
  if (is.list(accuracy) && inherits(accuracy$covariance_strategy, "cdid_gmm_covariance_strategy")) {
    accuracy$covariance_strategy$weighting_requested <- w_choice$requested
    accuracy$covariance_strategy$weighting_canonical <- w_choice$canonical
    accuracy$covariance_strategy$alias_used <- w_choice$alias_used
  }
  accuracy <- validate_gmm_accuracy_control(accuracy, w_choice$canonical)
  if (w_choice$canonical == "optimal" && is.null(accuracy$covariance_strategy)) {
    accuracy$covariance_strategy <- gmm_covariance_strategy(
      method = "tikhonov",
      condition_trigger = 1e-10,
      ridge_relative = 1e-4,
      weighting_requested = w_choice$requested,
      weighting_canonical = w_choice$canonical
    )
  }

  u_chr <- static_graph_validate_ids(unit_ids, "unit_ids", allow_duplicates = FALSE)
  N <- length(u_chr)
  clustered <- !is.null(cluster_ids)
  if (clustered) {
    c_chr <- static_graph_validate_ids(cluster_ids, "cluster_ids", allow_duplicates = TRUE)
    if (length(c_chr) != N) {
      cdid_static_graph_abort(
        paste0("`cluster_ids` length (", length(c_chr), ") must equal `length(unit_ids)` (", N, ")."),
        "cdid_static_graph_input_error"
      )
    }
    n_clusters <- length(unique(c_chr))
    if (n_clusters < 2L) {
      cdid_static_graph_abort(
        "At least 2 distinct clusters are required when `cluster_ids` is non-NULL.",
        "cdid_static_graph_input_error"
      )
    }
  } else {
    n_clusters <- N
  }

  if (is.null(normalization_factor)) {
    normalization_factor <- 1.0 / as.double(N)
  }
  if (!gmm_scalar_numeric(normalization_factor, positive = TRUE)) {
    cdid_static_graph_abort("`normalization_factor` must be a finite positive numeric scalar.", "cdid_static_graph_input_error")
  }
  norm_fac <- as.double(normalization_factor)

  # 2. Build joint `cdid_chained_design` and cell targets
  design <- build_static_graph_gmm_design(
    blocks = blocks,
    method = method,
    mhat = mhat,
    cell_ids = cell_ids
  )
  spec <- design$static_graph_spec
  mhat_vec <- spec$mhat
  if (is.null(mhat_vec)) {
    cdid_static_graph_abort(
      "Point moments `mhat` must be supplied either via `mhat` or on every role block.",
      "cdid_static_graph_alignment_error"
    )
  }
  link_ids <- design$incidence$link_ids
  m_R <- design$incidence$n_links
  p <- design$incidence$n_parameters
  J <- nrow(spec$cells)

  # All later passes use the same original-unit replay contract.
  if_callback <- static_graph_cross_stage_replay(if_callback, u_chr, link_ids, bounds)

  # 3. Stream complete-unit moment covariance V in bounded batches (or validate
  #    externally supplied V without fabricating `complete_unit_if_gram`)
  if (is.null(V)) {
    cov_stream <- stream_static_graph_moment_covariance(
      if_callback = if_callback,
      unit_ids = unit_ids,
      link_ids = link_ids,
      bounds = bounds,
      centered = centered,
      normalization_factor = norm_fac,
      nuisance_fit_id = nuisance_fit_id
    )
    V_use <- cov_stream$V
    if (is.null(provenance)) {
      prov_use <- gmm_covariance_provenance(
        row_level = "unit",
        n_units = N,
        n_clusters = N,
        n_positive_contributing_rows = N,
        centered = centered,
        weighted = FALSE,
        weight_type = "none",
        normalization_factor = norm_fac,
        scale_id = scale_id,
        nuisance_fit_id = nuisance_fit_id,
        reduction = list(mode = "none", reduction_active = FALSE)
      )
    } else {
      prov_use <- provenance
    }
    if (w_choice$canonical == "optimal" && !is.null(accuracy$covariance_strategy)) {
      prov_use$complete_unit_if_gram <- cov_stream$complete_unit_if_gram
    }
    gram_sealed <- isTRUE(w_choice$canonical == "optimal" && !is.null(accuracy$covariance_strategy))
  } else {
    # Externally supplied V: NEVER fabricate `complete_unit_if_gram` provenance
    if (is.matrix(V) && (!is.null(rownames(V)) || !is.null(colnames(V)))) {
      if ((!is.null(rownames(V)) && !identical(as.character(rownames(V)), link_ids)) ||
          (!is.null(colnames(V)) && !identical(as.character(colnames(V)), link_ids))) {
        cdid_static_graph_abort(
          "Externally supplied `V` row/column names do not match `design$incidence$link_ids`.",
          c("cdid_static_graph_alignment_error", "cdid_gmm_covariance_error")
        )
      }
    }
    V_use <- V
    if (is.null(provenance)) {
      prov_use <- gmm_covariance_provenance(
        row_level = "unit",
        n_units = N,
        n_clusters = N,
        n_positive_contributing_rows = N,
        centered = centered,
        weighted = FALSE,
        weight_type = "none",
        normalization_factor = norm_fac,
        scale_id = scale_id,
        nuisance_fit_id = nuisance_fit_id,
        reduction = list(mode = "none", reduction_active = FALSE)
      )
    } else {
      prov_use <- provenance
    }
    cov_stream <- list(
      V = V_use,
      mean_phi = NULL,
      raw_sum_phi = NULL,
      n_units = as.integer(N),
      n_links = as.integer(m_R),
      n_batches = NA_integer_,
      max_batch_rows_seen = NA_integer_,
      centered = prov_use$centered,
      normalization_factor = prov_use$normalization_factor,
      nuisance_fit_id = prov_use$nuisance_fit_id,
      complete_unit_if_gram = prov_use$complete_unit_if_gram
    )
    gram_sealed <- FALSE
  }

  # 4. Wrap callback for unit-level IF replay in `solve_chained_gmm_core`
  core_if_callback <- function(consume) {
    if_callback(function(batch) {
      if (!is.list(batch) || !("phi" %in% names(batch))) {
        cdid_static_graph_abort("Invalid IF batch in replay.", "cdid_gmm_replay_error")
      }
      phi_b <- batch$phi
      if (is.matrix(phi_b) && !is.null(colnames(phi_b)) &&
          !identical(as.character(colnames(phi_b)), link_ids)) {
        cdid_static_graph_abort(
          "IF batch `colnames(phi)` do not match `design$incidence$link_ids`.",
          c("cdid_static_graph_alignment_error", "cdid_gmm_replay_error")
        )
      }
      uid_b <- if ("unit_id" %in% names(batch)) batch$unit_id else batch$id
      out_b <- list(
        unit_id = uid_b,
        phi = phi_b,
        weight = if (!is.null(batch$weight)) batch$weight else rep(1.0, NROW(phi_b)),
        nuisance_fit_id = if (!is.null(batch$nuisance_fit_id)) batch$nuisance_fit_id else prov_use$nuisance_fit_id
      )
      consume(out_b)
    })
  }

  if_replay_unit <- gmm_if_replay_spec(
    callback = core_if_callback,
    input_row_level = "unit",
    aggregate_to_clusters = FALSE,
    n_input_rows = N,
    n_output_rows = N,
    n_positive_contributing_rows = N,
    centered = prov_use$centered,
    weighted = FALSE,
    weight_type = "none",
    normalization_factor = prov_use$normalization_factor,
    scale_id = prov_use$scale_id,
    nuisance_fit_id = prov_use$nuisance_fit_id,
    expected_unit_ids = unit_ids
  )

  # 5. Solve joint GMM system via `solve_chained_gmm_core`
  gmm_fit <- solve_chained_gmm_core(
    design = design,
    mhat = mhat_vec,
    V = V_use,
    provenance = prov_use,
    targets = spec$target_specs,
    weighting = weighting,
    if_replay = if_replay_unit,
    accuracy = accuracy,
    bounds = bounds,
    gauge_rule = gauge_rule,
    custom_gauge_roots = custom_gauge_roots,
    request_param_cov = request_param_cov,
    materialize_target_if = materialize_target_if
  )

  Lambda <- gmm_fit$internal_gauge$Lambda
  unit_sandwich_cov <- gmm_fit$sandwich_cov
  unit_se_sandwich <- gmm_fit$targets$se_sandwich

  # 6. Cluster inference-only step (points and unit-GMM weights unchanged)
  cluster_res <- NULL
  if (clustered) {
    cluster_res <- stream_static_graph_cluster_sandwich(
      if_callback = if_callback,
      Lambda = Lambda,
      unit_ids = unit_ids,
      cluster_ids = cluster_ids,
      link_ids = link_ids,
      bounds = bounds,
      centered = prov_use$centered,
      normalization_factor = prov_use$normalization_factor,
      nuisance_fit_id = prov_use$nuisance_fit_id,
      materialize_cluster_if = materialize_target_if
    )
    active_sandwich_cov <- cluster_res$sandwich_cov
    active_se <- cluster_res$se_sandwich
  } else {
    active_sandwich_cov <- unit_sandwich_cov
    active_se <- unit_se_sandwich
  }

  # Asymptotic covariance (1/N scale) and finite-sample estimator vcov (1/N^2 scale)
  raw_S_cov <- active_sandwich_cov / prov_use$normalization_factor
  asymptotic_cov <- raw_S_cov / as.double(N)
  vcov_n2 <- raw_S_cov / (as.double(N)^2)
  se_n2 <- sqrt(pmax(0.0, diag(vcov_n2)))

  target_if_mat <- if (isTRUE(materialize_target_if)) {
    mat <- gmm_fit$if_projection$materialized_target_if
    rownames(mat) <- as.character(unit_ids)
    colnames(mat) <- spec$cells$cell_id
    mat
  } else {
    NULL
  }

  cell_effects <- spec$cells
  cell_effects$estimable <- gmm_fit$targets$estimable
  cell_effects$status <- gmm_fit$targets$status
  cell_effects$estimate <- gmm_fit$targets$estimate
  cell_effects$se <- as.double(se_n2)
  cell_effects$se_unit <- sqrt(pmax(0.0, diag(unit_sandwich_cov) /
                                    prov_use$normalization_factor / as.double(N)^2))
  cell_effects$se_n2 <- as.double(se_n2)

  V_named <- V_use
  if (is.matrix(V_named) && nrow(V_named) == m_R && ncol(V_named) == m_R) {
    dimnames(V_named) <- list(link_ids, link_ids)
  }

  param_out <- design$parameters
  param_out$level_estimate <- gmm_fit$parameters$level_estimate
  param_out$status <- gmm_fit$parameters$status

  structure(
    list(
      cell_effects = cell_effects,
      estimates = stats::setNames(cell_effects$estimate, cell_effects$cell_id),
      se = stats::setNames(cell_effects$se, cell_effects$cell_id),
      vcov = vcov_n2,
      sandwich_cov = vcov_n2,
      core_se = stats::setNames(as.double(active_se), cell_effects$cell_id),
      core_sandwich_cov = active_sandwich_cov,
      unit_sandwich_cov = unit_sandwich_cov,
      asymptotic_cov = asymptotic_cov,
      vcov_n2 = vcov_n2,
      target_if = target_if_mat,
      cluster_target_if = if (!is.null(cluster_res)) cluster_res$cluster_target_if else NULL,
      moment_covariance = V_named,
      moment_weights = list(
        Lambda = Lambda,
        lambda = Lambda
      ),
      mapping = list(
        blocks = spec$blocks,
        cells = spec$cells,
        parameters = param_out,
        links = design$links,
        targets = spec$target_specs,
        unit_ids = unit_ids,
        cluster_ids = cluster_ids
      ),
      diagnostics = list(
        method = spec$method,
        weighting_requested = w_choice$requested,
        weighting_canonical = w_choice$canonical,
        n_units = as.integer(N),
        n_clusters = as.integer(n_clusters),
        clustered = clustered,
        n_cells = as.integer(J),
        n_blocks = nrow(spec$blocks),
        n_links = as.integer(m_R),
        n_parameters = as.integer(p),
        rank = design$rank_diagnostics$rank,
        overid_df = gmm_fit$overid_df,
        raw_quadratic_residual = gmm_fit$raw_quadratic_residual,
        covariance_handling = gmm_fit$covariance_handling,
        rank_diagnostics = design$rank_diagnostics,
        numerical_guard = gmm_fit$numerical_guard,
        allocations = gmm_fit$allocations,
        stream_summary = list(
          n_units = cov_stream$n_units,
          n_links = cov_stream$n_links,
          n_batches = cov_stream$n_batches,
          max_batch_rows_seen = cov_stream$max_batch_rows_seen,
          mean_phi = cov_stream$mean_phi,
          centered = prov_use$centered,
          normalization_factor = prov_use$normalization_factor,
          scale_id = prov_use$scale_id,
          nuisance_fit_id = prov_use$nuisance_fit_id,
          gram_provenance_sealed = gram_sealed
        ),
        normalization_convention = list(
          centered = prov_use$centered,
          normalization_factor = prov_use$normalization_factor,
          moment_covariance_formula = if (prov_use$centered) {
            "V = normalization_factor * sum_{i=1}^N (chi_i - bar_chi)(chi_i - bar_chi)'"
          } else {
            "V = normalization_factor * sum_{i=1}^N chi_i chi_i'"
          },
          sandwich_covariance_formula = if (clustered) {
            if (prov_use$centered) {
              "sandwich_cov = normalization_factor * sum_{g=1}^G (tilde_phi_g - bar_tilde_phi)(tilde_phi_g - bar_tilde_phi)', tilde_phi_g = sum_{i in g} chi_i %*% Lambda"
            } else {
              "sandwich_cov = normalization_factor * sum_{g=1}^G tilde_phi_g tilde_phi_g', tilde_phi_g = sum_{i in g} chi_i %*% Lambda"
            }
          } else {
            if (prov_use$centered) {
              "sandwich_cov = normalization_factor * sum_{i=1}^N (phi_i - bar_phi)(phi_i - bar_phi)', phi_i = chi_i %*% Lambda"
            } else {
              "sandwich_cov = normalization_factor * sum_{i=1}^N phi_i phi_i', phi_i = chi_i %*% Lambda"
            }
          },
          asymptotic_cov_scale = "1 / N",
          vcov_n2_scale = "1 / N^2"
        )
      ),
      design = design,
      gmm_fit = gmm_fit
    ),
    class = "cdid_static_graph_gmm_result"
  )
}
