# Private bounded adapter between the static observation-nuisance kernel
# (R/static_nuisance_models.R: observation logits, change imputation WLS, and
# first-step original-unit influence functions) and the static observation-moment
# kernel (R/static_observation_moments.R: reconstructed role-pair means and
# covariances, final pair-probability lower floor, and complete influence functions).
#
# Responsibilities:
#   * validate the compact fixed-group / pair / date input contract on a common
#     original unit_id vector of length N (zero- and single-outcome units stay
#     in N and in role masks);
#   * construct at most ONE selection logit fit (pair_logit or
#     independent_period_logit) and at most ONE actual-change imputation WLS
#     fit across all declared fixed role groups and edges when the route uses
#     them, under explicit pooling (separate / within_group / pooled) and
#     imputation_scale (per_period / per_change);
#   * expose fitted selection/imputation objects and an on-demand edge moment
#     interface (mean with caller-supplied tilt h + h_jacobian + h_coef_if;
#     covariance with stayer Z) that calls static_obs_mean / static_obs_cov
#     with full-N normalization, discards full-edge prediction/Jacobian/direct
#     IF temporaries, and provides eval_if(rows) by replaying the direct IF from
#     the frozen full-N summary (N, S, M and C, Zbar) plus contracting
#     loadings against static_nuisance_if(fit, rows) and h_coef_if(rows).
#
# Not implemented here: public att_gt/formula parsing, graph selection,
# counterfactual treatment OR reconstruction, or final GMM assembly.

static_obs_adapter_routes <- c("baseline", "weighting", "imputation", "augmented")
static_obs_adapter_sel_models <- c("pair_logit", "independent_period_logit")
static_obs_adapter_poolings <- c("separate", "within_group", "pooled")
static_obs_adapter_imp_scales <- c("per_period", "per_change")
static_obs_adapter_roles <- c("G", "C")

static_obs_adapter_abort <- function(message, subclass = "cdid_static_obs_adapter_input_error") {
  cdid_abort(message, c(subclass, "cdid_static_obs_adapter_error"))
}

static_obs_adapter_key_ok <- function(x) {
  is.character(x) && length(x) >= 1L && !anyNA(x) && all(nzchar(x)) &&
    !any(grepl("|", x, fixed = TRUE)) && !any(grepl("::", x, fixed = TRUE)) &&
    !anyDuplicated(x)
}

static_obs_adapter_scalar_key_ok <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x) &&
    !grepl("|", x, fixed = TRUE) && !grepl("::", x, fixed = TRUE)
}

static_obs_adapter_choice <- function(x, allowed, what) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !(x %in% allowed)) {
    static_obs_adapter_abort(
      sprintf("`%s` must be explicitly supplied as one of: %s.",
              what, paste(allowed, collapse = ", "))
    )
  }
  x
}

static_obs_adapter_mask <- function(x, N, ids, what) {
  if (is.null(x) || !is.atomic(x) || !is.null(dim(x)) || is.object(x) ||
      !(is.logical(x) || is.numeric(x)) || length(x) != N || anyNA(x)) {
    static_obs_adapter_abort(
      sprintf("`%s` must be a non-NA logical (or 0/1) vector of length N.", what)
    )
  }
  if (!is.null(names(x)) && !identical(as.character(names(x)), ids)) {
    static_obs_adapter_abort(
      sprintf("Names of `%s` do not match `unit_id` in order.", what),
      "cdid_static_obs_order_error"
    )
  }
  if (is.numeric(x)) {
    v <- as.double(x)
    if (any(v != 0 & v != 1)) {
      static_obs_adapter_abort(sprintf("`%s` must take only 0/1 values.", what))
    }
    return(v == 1)
  }
  as.logical(x)
}

static_obs_adapter_role_obs <- function(x, N, ids, inA, what) {
  if (is.null(x) || !is.atomic(x) || !is.null(dim(x)) || is.object(x) ||
      !(is.logical(x) || is.numeric(x)) || length(x) != N) {
    static_obs_adapter_abort(
      sprintf("`%s` must be a logical or 0/1 vector of length N.", what)
    )
  }
  if (!is.null(names(x)) && !identical(as.character(names(x)), ids)) {
    static_obs_adapter_abort(
      sprintf("Names of `%s` do not match `unit_id` in order.", what),
      "cdid_static_obs_order_error"
    )
  }
  v <- as.double(x[inA])
  if (anyNA(v) || any(v != 0 & v != 1)) {
    static_obs_adapter_abort(
      sprintf("`%s` must be 0/1 (or TRUE/FALSE) without NA on the group role mask.", what)
    )
  }
  out_log <- logical(N)
  out_log[inA] <- (v == 1)
  out_num <- numeric(N)
  out_num[inA] <- v
  list(logical = out_log, numeric = out_num)
}

static_obs_adapter_dy <- function(x, N, ids, R_o, what) {
  if (is.null(x) || !is.atomic(x) || !is.null(dim(x)) || is.object(x) ||
      !(is.numeric(x) || (is.logical(x) && all(is.na(x)))) || length(x) != N) {
    static_obs_adapter_abort(
      sprintf("`%s` must be a numeric vector of length N.", what)
    )
  }
  if (!is.null(names(x)) && !identical(as.character(names(x)), ids)) {
    static_obs_adapter_abort(
      sprintf("Names of `%s` do not match `unit_id` in order.", what),
      "cdid_static_obs_order_error"
    )
  }
  v_o <- as.double(x[R_o])
  if (length(v_o) && any(!is.finite(v_o))) {
    static_obs_adapter_abort(
      sprintf("`%s` must be finite on observed role units (A & O).", what)
    )
  }
  out <- rep(NA_real_, N)
  out[R_o] <- v_o
  out
}

# Extracts and validates a design matrix, reading ONLY rows in `read_rows`.
# Rows outside `read_rows` are filled with NA_real_ so nonfinite entries outside
# the relevant population are never touched or multiplied.
static_obs_adapter_design <- function(X, N, ids, read_rows, scale_factor = 1, what = "design") {
  if (is.null(X)) return(NULL)
  if (!is.matrix(X) || !is.numeric(X) || is.object(X) || nrow(X) != N) {
    static_obs_adapter_abort(
      sprintf("`%s` must be NULL or a numeric matrix with N rows.", what)
    )
  }
  if (!is.null(rownames(X)) && !identical(as.character(rownames(X)), ids)) {
    static_obs_adapter_abort(
      sprintf("Rownames of `%s` do not match `unit_id` in order.", what),
      "cdid_static_obs_order_error"
    )
  }
  d_s <- ncol(X)
  if (d_s == 0L) {
    return(matrix(numeric(0), nrow = N, ncol = 0L))
  }
  cn <- colnames(X)
  if (is.null(cn) || !static_nuisance_key_ok(cn)) {
    static_obs_adapter_abort(
      sprintf("`%s` must have unique non-empty column names without '|'.", what)
    )
  }
  if (any(cn %in% c("(Intercept)", "Intercept"))) {
    static_obs_adapter_abort(
      sprintf("`%s` must not contain an explicit intercept column.", what)
    )
  }
  sub <- unname(X[read_rows, , drop = FALSE])
  if (length(sub) && any(!is.finite(sub))) {
    static_obs_adapter_abort(
      sprintf("`%s` must be finite on its required population.", what)
    )
  }
  if (scale_factor != 1) {
    sub <- scale_factor * sub
  }
  out <- matrix(NA_real_, nrow = N, ncol = d_s, dimnames = list(NULL, cn))
  out[read_rows, ] <- sub
  out
}

static_obs_adapter_gap <- function(gap, what) {
  if (is.null(gap) || !is.numeric(gap) || length(gap) != 1L || !is.finite(gap) ||
      gap < 1 || gap != floor(gap)) {
    static_obs_adapter_abort(
      sprintf("`%s` must be an explicit positive integer panel-index gap (never inferred from calendar dates).", what)
    )
  }
  as.integer(gap)
}

static_obs_adapter_partition <- function(pooling, group_cell_map) {
  all_cells <- unlist(group_cell_map, use.names = FALSE)
  if (pooling == "separate") {
    part <- as.list(all_cells)
    names(part) <- all_cells
    return(part)
  }
  if (pooling == "within_group") {
    return(group_cell_map)
  }
  list(pooled = all_cells)
}

# Resolves requested original-unit IDs (or NULL -> all unit_id in order) into
# integer positions in 1..N while enforcing type/uniqueness/membership checks.
static_obs_adapter_rows <- function(unit_id, rows) {
  if (is.null(rows)) {
    return(list(ids = unit_id, pos = seq_along(unit_id)))
  }
  if (!is.atomic(rows) || is.object(rows) || !is.null(dim(rows)) || length(rows) < 1L ||
      anyNA(rows) || !(is.character(rows) || is.numeric(rows)) ||
      anyDuplicated(rows) || is.character(rows) != is.character(unit_id)) {
    static_obs_adapter_abort(
      "Requested `rows` must be unique non-NA original unit IDs of `unit_id`'s type.",
      "cdid_static_nuisance_input_error"
    )
  }
  pos <- match(rows, unit_id)
  if (anyNA(pos)) {
    static_obs_adapter_abort(
      "Requested `rows` contain unknown original unit IDs.",
      "cdid_static_nuisance_input_error"
    )
  }
  list(ids = rows, pos = pos)
}

# Embeds a compact prediction derivative (`pred$deriv`, `pred$param_index`) on
# role rows `inA` into full-N block-coordinate Jacobian matrices keyed by
# "<prefix>|<block>". Each block channel always carries the full parameter
# coordinate vector of that fitted block in declared fit order.
static_obs_adapter_block_jacobians <- function(fit, prefix, inA, pred, N, ids_chr) {
  out <- list()
  for (b in names(fit$layout$blocks)) {
    bl <- fit$layout$blocks[[b]]
    hit <- which(pred$param_index %in% bl$index)
    if (!length(hit)) next
    p_b <- length(bl$index)
    pkeys <- fit$layout$param_keys[bl$index]
    M_b <- matrix(0, nrow = N, ncol = p_b, dimnames = list(ids_chr, pkeys))
    col_pos <- match(pred$param_index[hit], bl$index)
    M_b[inA, col_pos] <- pred$deriv[, hit, drop = FALSE]
    out[[paste0(prefix, "|", b)]] <- M_b
  }
  out
}

# Evaluates and slices block-keyed original-unit coefficient IFs for the
# required channel keys on `ids` (in chunks of at most `max_chunk_rows`), or
# slices a caller-supplied operation-local `full_if` matrix directly.
static_obs_adapter_fit_if_channels <- function(fit, prefix, channel_keys, ids, max_chunk_rows,
                                               full_if = NULL) {
  if (!length(channel_keys)) return(list())
  if (is.null(fit)) {
    static_obs_adapter_abort(
      sprintf("Internal error: missing `%s` fit for active loading channels.", prefix),
      "cdid_static_obs_channel_error"
    )
  }
  if (is.null(full_if)) {
    n_req <- length(ids)
    P <- length(fit$coefficients)
    full_if <- matrix(0, nrow = n_req, ncol = P,
                      dimnames = list(as.character(ids), fit$layout$param_keys))
    for (start in seq.int(1L, n_req, by = max_chunk_rows)) {
      idx <- seq.int(start, min(n_req, start + max_chunk_rows - 1L))
      full_if[idx, ] <- static_nuisance_if(fit, ids[idx])
    }
  }
  out <- list()
  for (ch in channel_keys) {
    b <- sub(paste0("^", prefix, "\\|"), "", ch)
    bl <- fit$layout$blocks[[b]]
    if (is.null(bl)) {
      static_obs_adapter_abort(
        sprintf("Unknown %s block channel '%s'.", prefix, ch),
        "cdid_static_obs_channel_error"
      )
    }
    out[[ch]] <- full_if[, bl$index, drop = FALSE]
  }
  out
}

# Resolves caller-supplied `h_coef_if` (provider function or keyed list of
# matrices/functions) on requested `ids`.
static_obs_adapter_eval_h_if <- function(h_coef_if, h_keys, ids, pos, unit_id_chr) {
  if (!length(h_keys)) return(list())
  ids_chr <- as.character(ids)
  n_req <- length(ids)
  raw <- if (is.function(h_coef_if)) h_coef_if(ids) else h_coef_if
  if (is.null(raw) || !is.list(raw) || is.object(raw)) {
    static_obs_adapter_abort(
      "`h_coef_if` must provide a keyed list of coefficient-IF matrices for all `h_jacobian` channels.",
      "cdid_static_obs_channel_error"
    )
  }
  rk <- names(raw)
  if (is.null(rk) || anyNA(rk) || any(!nzchar(rk)) || anyDuplicated(rk) ||
      !setequal(rk, h_keys)) {
    static_obs_adapter_abort(
      "`h_coef_if` keys must match `h_jacobian` channel keys exactly.",
      "cdid_static_obs_channel_error"
    )
  }
  out <- list()
  for (k in h_keys) {
    entry <- raw[[k]]
    if (is.function(entry)) {
      mat <- entry(ids)
    } else if (is.matrix(entry) && is.numeric(entry)) {
      if (nrow(entry) == length(unit_id_chr) &&
          !is.null(rownames(entry)) && identical(rownames(entry), unit_id_chr)) {
        mat <- entry[pos, , drop = FALSE]
      } else if (nrow(entry) == n_req &&
                 !is.null(rownames(entry)) && identical(rownames(entry), ids_chr)) {
        mat <- entry
      } else {
        static_obs_adapter_abort(
          sprintf("`h_coef_if[['%s']]` must have labeled rows matching `unit_id` or requested `rows` in order.", k),
          "cdid_static_obs_order_error"
        )
      }
    } else {
      static_obs_adapter_abort(
        sprintf("`h_coef_if[['%s']]` must be a numeric matrix or a provider function(ids).", k),
        "cdid_static_obs_channel_error"
      )
    }
    out[[k]] <- mat
  }
  out
}

static_observation_adapter_fit <- function(unit_id,
                                           groups,
                                           route,
                                           selection_model = NULL,
                                           selection_pooling = NULL,
                                           imputation_pooling = NULL,
                                           imputation_scale = NULL,
                                           glm_control = NULL,
                                           rank_tol = NULL,
                                           max_chunk_rows,
                                           selection_design = NULL,
                                           imputation_design = NULL,
                                           dates = NULL,
                                           pair_callback = NULL) {
  static_nuisance_validate_ids(unit_id)
  ids_chr <- static_obs_ids(unit_id)
  N <- length(unit_id)

  route <- static_obs_adapter_choice(route, static_obs_adapter_routes, "route")
  uses_O <- route != "imputation"
  uses_q <- route %in% c("weighting", "augmented")
  uses_j <- route %in% c("imputation", "augmented")

  if (uses_q || !is.null(selection_model)) {
    selection_model <- static_obs_adapter_choice(
      selection_model, static_obs_adapter_sel_models, "selection_model"
    )
  }
  if (uses_q || !is.null(selection_pooling)) {
    selection_pooling <- static_obs_adapter_choice(
      selection_pooling, static_obs_adapter_poolings, "selection_pooling"
    )
  }
  if (uses_j || !is.null(imputation_pooling)) {
    imputation_pooling <- static_obs_adapter_choice(
      imputation_pooling, static_obs_adapter_poolings, "imputation_pooling"
    )
  }
  if (uses_j || !is.null(imputation_scale)) {
    imputation_scale <- static_obs_adapter_choice(
      imputation_scale, static_obs_adapter_imp_scales, "imputation_scale"
    )
  }
  if (uses_q && is.null(glm_control)) {
    static_obs_adapter_abort("`glm_control` must be supplied explicitly when route uses selection.")
  }
  if (!is.null(glm_control)) {
    glm_control <- validate_ipw_glm_control(glm_control)
  }
  if (uses_j && (is.null(rank_tol) || !is.numeric(rank_tol) || length(rank_tol) != 1L ||
                 !is.finite(rank_tol) || rank_tol <= 0)) {
    static_obs_adapter_abort("`rank_tol` must be supplied explicitly as a positive finite scalar when route uses imputation.")
  }
  if (missing(max_chunk_rows) || is.null(max_chunk_rows)) {
    static_obs_adapter_abort(
      "`max_chunk_rows` must be supplied explicitly as a positive integer.",
      c("cdid_static_nuisance_input_error", "cdid_static_obs_adapter_input_error")
    )
  }
  chunk <- static_nuisance_chunk(max_chunk_rows, N)

  if (!is.null(pair_callback) && !is.function(pair_callback)) {
    static_obs_adapter_abort("`pair_callback` must be NULL or a function(group_key, pair_key).")
  }
  if (!is.list(groups) || length(groups) < 1L || !static_obs_adapter_key_ok(names(groups))) {
    static_obs_adapter_abort(
      "`groups` must be a non-empty named list with unique keys without '|' or '::'."
    )
  }
  if (!is.null(dates) && (!is.list(dates) || length(dates) < 1L || !static_obs_adapter_key_ok(names(dates)))) {
    static_obs_adapter_abort(
      "`dates` must be NULL or a non-empty named list with unique keys without '|' or '::'."
    )
  }

  sel_cells <- list()
  sel_group_cells <- list()
  imp_cells <- list()
  imp_group_cells <- list()
  parsed_groups <- vector("list", length(groups))
  names(parsed_groups) <- names(groups)

  for (g_key in names(groups)) {
    grp <- groups[[g_key]]
    if (!is.list(grp)) {
      static_obs_adapter_abort(sprintf("Group '%s' must be a list.", g_key))
    }
    role <- if (!is.null(grp$role)) grp$role else NULL
    if (is.null(role) || !is.character(role) || length(role) != 1L ||
        is.na(role) || !(role %in% static_obs_adapter_roles)) {
      static_obs_adapter_abort(
        sprintf("Group '%s' must declare `role` as 'G' or 'C'.", g_key)
      )
    }
    A <- static_obs_adapter_mask(grp$mask, N, ids_chr, sprintf("groups[['%s']]$mask", g_key))
    inA <- which(A)
    if (!length(inA)) {
      static_obs_adapter_abort(
        sprintf("Group '%s' has an empty role population.", g_key),
        "cdid_static_obs_support_error"
      )
    }

    prs <- grp$pairs
    if (!is.list(prs) || length(prs) < 1L || !static_obs_adapter_key_ok(names(prs))) {
      static_obs_adapter_abort(
        sprintf("Group '%s': `pairs` must be a non-empty named list with unique keys without '|' or '::'.", g_key)
      )
    }

    g_sel_design_raw <- if (!is.null(grp$selection_design)) grp$selection_design else selection_design
    g_imp_design_raw <- if (!is.null(grp$imputation_design)) grp$imputation_design else imputation_design

    parsed_pairs <- vector("list", length(prs))
    names(parsed_pairs) <- names(prs)
    g_pair_cells <- character(0)
    needed_dates <- character(0)

    for (p_key in names(prs)) {
      pr <- prs[[p_key]]
      if (!is.list(pr)) {
        static_obs_adapter_abort(sprintf("Group '%s', pair '%s' must be a list.", g_key, p_key))
      }
      if (!is.null(pair_callback) && (is.null(pr$observed) && is.null(pr$O) && is.null(pr$dy))) {
        cb <- pair_callback(g_key, p_key)
        if (!is.list(cb)) {
          static_obs_adapter_abort(
            sprintf("`pair_callback('%s', '%s')` must return a list.", g_key, p_key)
          )
        }
        for (nm in names(cb)) {
          if (is.null(pr[[nm]])) pr[[nm]] <- cb[[nm]]
        }
      }
      start_d <- if (!is.null(pr$start)) as.character(pr$start) else NULL
      end_d <- if (!is.null(pr$end)) as.character(pr$end) else NULL
      if (!static_obs_adapter_scalar_key_ok(start_d) || !static_obs_adapter_scalar_key_ok(end_d) ||
          identical(start_d, end_d)) {
        static_obs_adapter_abort(
          sprintf("Group '%s', pair '%s' must supply distinct scalar `start` and `end` date keys without '|' or '::'.",
                  g_key, p_key)
        )
      }
      gap <- static_obs_adapter_gap(pr$gap, sprintf("groups[['%s']]$pairs[['%s']]$gap", g_key, p_key))

      raw_obs <- if (!is.null(pr$observed)) pr$observed else pr$O
      need_pair_obs <- uses_O || uses_j || (uses_q && identical(selection_model, "pair_logit"))
      obs_info <- NULL
      R_o <- integer(0)
      dy_clean <- NULL
      if (need_pair_obs) {
        obs_info <- static_obs_adapter_role_obs(
          raw_obs, N, ids_chr, inA,
          sprintf("groups[['%s']]$pairs[['%s']]$observed", g_key, p_key)
        )
        R_o <- which(A & obs_info$logical)
      }
      if (uses_O || uses_j) {
        dy_clean <- static_obs_adapter_dy(
          pr$dy, N, ids_chr, R_o,
          sprintf("groups[['%s']]$pairs[['%s']]$dy", g_key, p_key)
        )
      }

      cell_key <- paste0(g_key, "::", p_key)
      g_pair_cells <- c(g_pair_cells, cell_key)
      needed_dates <- c(needed_dates, start_d, end_d)

      if (uses_q && identical(selection_model, "pair_logit")) {
        p_sel_raw <- if (!is.null(pr$selection_design)) pr$selection_design else g_sel_design_raw
        p_sel_des <- static_obs_adapter_design(
          p_sel_raw, N, ids_chr, inA, scale_factor = 1,
          what = sprintf("selection_design (%s, %s)", g_key, p_key)
        )
        sel_cells[[cell_key]] <- list(
          mask = A,
          response = obs_info$numeric,
          design = p_sel_des
        )
      }

      if (uses_j) {
        p_imp_raw <- if (!is.null(pr$imputation_design)) pr$imputation_design else g_imp_design_raw
        scale_fac <- if (identical(imputation_scale, "per_period")) gap else 1L
        # Read imputation design on inA (needed for j_predict over A and wls_fit on A & O)
        p_imp_des <- static_obs_adapter_design(
          p_imp_raw, N, ids_chr, inA, scale_factor = scale_fac,
          what = sprintf("imputation_design (%s, %s)", g_key, p_key)
        )
        imp_cells[[cell_key]] <- list(
          mask = A,
          observed = obs_info$logical,
          dy = dy_clean,
          design = p_imp_des
        )
      }

      parsed_pairs[[p_key]] <- list(
        pair = p_key,
        cell_key = cell_key,
        start = start_d,
        end = end_d,
        gap = gap,
        date_cell_keys = c(start = paste0(g_key, "::", start_d),
                           end = paste0(g_key, "::", end_d)),
        observed = if (!is.null(obs_info)) obs_info$logical else NULL,
        dy = dy_clean
      )
    }

    if (uses_q && identical(selection_model, "pair_logit")) {
      sel_group_cells[[g_key]] <- g_pair_cells
    }
    if (uses_j) {
      imp_group_cells[[g_key]] <- g_pair_cells
    }

    if (uses_q && identical(selection_model, "independent_period_logit")) {
      g_dates_src <- if (!is.null(grp$dates)) grp$dates else dates
      if (is.null(g_dates_src) || !is.list(g_dates_src) || length(g_dates_src) < 1L ||
          !static_obs_adapter_key_ok(names(g_dates_src))) {
        static_obs_adapter_abort(
          sprintf("Group '%s': `dates` records must be supplied explicitly for `independent_period_logit`.", g_key)
        )
      }
      req_dates <- unique(needed_dates)
      if (!all(req_dates %in% names(g_dates_src))) {
        miss_d <- setdiff(req_dates, names(g_dates_src))
        static_obs_adapter_abort(
          sprintf("Group '%s': missing per-date observation records for date(s): %s.",
                  g_key, paste(miss_d, collapse = ", "))
        )
      }
      # Preserve declared date order in g_dates_src for dates used by this group
      # (or all declared dates in grp$dates if group-specific).
      date_keys_use <- if (!is.null(grp$dates)) {
        names(grp$dates)
      } else {
        names(g_dates_src)[names(g_dates_src) %in% req_dates]
      }
      g_date_cells <- character(0)
      for (d_key in date_keys_use) {
        dr <- g_dates_src[[d_key]]
        if (!is.list(dr)) {
          static_obs_adapter_abort(
            sprintf("Group '%s', date '%s' record must be a list.", g_key, d_key)
          )
        }
        raw_d_obs <- if (!is.null(dr$observed)) dr$observed else dr$response
        if (is.null(raw_d_obs)) raw_d_obs <- dr$O
        d_obs <- static_obs_adapter_role_obs(
          raw_d_obs, N, ids_chr, inA,
          sprintf("dates[['%s']]$observed (group '%s')", d_key, g_key)
        )
        d_sel_raw <- if (!is.null(dr$selection_design)) dr$selection_design else
          if (!is.null(dr$design)) dr$design else g_sel_design_raw
        d_sel_des <- static_obs_adapter_design(
          d_sel_raw, N, ids_chr, inA, scale_factor = 1,
          what = sprintf("selection_design (group '%s', date '%s')", g_key, d_key)
        )
        d_cell_key <- paste0(g_key, "::", d_key)
        g_date_cells <- c(g_date_cells, d_cell_key)
        sel_cells[[d_cell_key]] <- list(
          mask = A,
          response = d_obs$numeric,
          design = d_sel_des
        )
      }
      sel_group_cells[[g_key]] <- g_date_cells
    }

    parsed_groups[[g_key]] <- list(
      group = g_key,
      role = role,
      target = grp$target,
      prior_state = grp$prior_state,
      mask = A,
      inA = inA,
      Z = grp$Z,
      pairs = parsed_pairs
    )
  }

  selection_fit <- NULL
  selection_spec <- NULL
  if (uses_q) {
    sel_partition <- static_obs_adapter_partition(selection_pooling, sel_group_cells)
    selection_fit <- static_nuisance_logit_fit(
      unit_id = unit_id,
      cells = sel_cells,
      partition = sel_partition,
      glm_control = glm_control,
      max_chunk_rows = chunk
    )
    selection_spec <- list(
      model = selection_model,
      pooling = selection_pooling,
      partition = sel_partition,
      group_cells = sel_group_cells
    )
  }

  imputation_fit <- NULL
  imputation_spec <- NULL
  if (uses_j) {
    imp_partition <- static_obs_adapter_partition(imputation_pooling, imp_group_cells)
    imputation_fit <- static_nuisance_wls_fit(
      unit_id = unit_id,
      cells = imp_cells,
      partition = imp_partition,
      rank_tol = rank_tol,
      max_chunk_rows = chunk
    )
    imputation_spec <- list(
      pooling = imputation_pooling,
      scale = imputation_scale,
      partition = imp_partition,
      group_cells = imp_group_cells
    )
  }

  # Private operation-local holder for one ordered original-unit batch of
  # built-in selection and imputation coefficient IF matrices. Activated only
  # within an admitted `with_if_batch()` scope (at most `chunk` rows) and
  # cleared on normal return or error via `on.exit`.
  batch_active <- FALSE
  batch_pos <- NULL
  batch_sel_if <- NULL
  batch_imp_if <- NULL

  with_if_batch <- function(rows, fn, enabled = TRUE) {
    if (!is.function(fn)) {
      static_obs_adapter_abort("`fn` must be a function.", "cdid_static_obs_adapter_error")
    }
    if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
      static_obs_adapter_abort("`enabled` must be a scalar logical.", "cdid_static_obs_adapter_error")
    }
    rk <- static_obs_adapter_rows(unit_id, rows)
    if (batch_active || !isTRUE(enabled) || length(rk$pos) > chunk || !(uses_q || uses_j)) {
      return(fn())
    }
    batch_active <<- TRUE
    batch_pos <<- as.integer(rk$pos)
    batch_sel_if <<- NULL
    batch_imp_if <<- NULL
    on.exit({
      batch_active <<- FALSE
      batch_pos <<- NULL
      batch_sel_if <<- NULL
      batch_imp_if <<- NULL
    }, add = TRUE)
    fn()
  }

  get_batch_selection_if <- function(req_ids, req_pos) {
    if (!batch_active || is.null(selection_fit) || !identical(as.integer(req_pos), batch_pos)) {
      return(NULL)
    }
    if (is.null(batch_sel_if)) {
      batch_sel_if <<- static_nuisance_if(selection_fit, req_ids)
    }
    batch_sel_if
  }

  get_batch_imputation_if <- function(req_ids, req_pos) {
    if (!batch_active || is.null(imputation_fit) || !identical(as.integer(req_pos), batch_pos)) {
      return(NULL)
    }
    if (is.null(batch_imp_if)) {
      batch_imp_if <<- static_nuisance_if(imputation_fit, req_ids)
    }
    batch_imp_if
  }

  resolve_edge <- function(group, pair) {
    if (!is.character(group) || length(group) != 1L || is.na(group) ||
        !(group %in% names(parsed_groups))) {
      static_obs_adapter_abort(sprintf("Unknown group key '%s'.", as.character(group)))
    }
    grp <- parsed_groups[[group]]
    if (!is.character(pair) || length(pair) != 1L || is.na(pair) ||
        !(pair %in% names(grp$pairs))) {
      static_obs_adapter_abort(
        sprintf("Unknown pair key '%s' in group '%s'.", as.character(pair), group)
      )
    }
    list(grp = grp, pr = grp$pairs[[pair]])
  }

  predict_edge_q_full <- function(grp, pr) {
    role_ids <- unit_id[grp$inA]
    pred <- if (identical(selection_model, "pair_logit")) {
      static_nuisance_q_predict(selection_fit, ids = role_ids, pair = pr$cell_key)
    } else {
      static_nuisance_q_predict(selection_fit, ids = role_ids, dates = unname(pr$date_cell_keys))
    }
    q_full <- rep(NA_real_, N)
    names(q_full) <- ids_chr
    q_full[grp$inA] <- pred$value
    q_jac <- static_obs_adapter_block_jacobians(
      selection_fit, "selection", grp$inA, pred, N, ids_chr
    )
    list(q = q_full, q_jacobian = q_jac, assumption = pred$assumption)
  }

  predict_edge_j_full <- function(grp, pr) {
    role_ids <- unit_id[grp$inA]
    pred <- static_nuisance_j_predict(imputation_fit, ids = role_ids, cell = pr$cell_key)
    j_full <- rep(NA_real_, N)
    names(j_full) <- ids_chr
    j_full[grp$inA] <- pred$value
    j_jac <- static_obs_adapter_block_jacobians(
      imputation_fit, "imputation", grp$inA, pred, N, ids_chr
    )
    list(j = j_full, j_jacobian = j_jac)
  }

  predict_edge_q_values <- function(pr, req_ids) {
    if (!length(req_ids)) return(numeric(0))
    out <- numeric(length(req_ids))
    for (start in seq.int(1L, length(req_ids), by = chunk)) {
      idx <- seq.int(start, min(length(req_ids), start + chunk - 1L))
      sub_ids <- req_ids[idx]
      pred <- if (identical(selection_model, "pair_logit")) {
        static_nuisance_q_predict(selection_fit, ids = sub_ids, pair = pr$cell_key)
      } else {
        static_nuisance_q_predict(selection_fit, ids = sub_ids, dates = unname(pr$date_cell_keys))
      }
      out[idx] <- pred$value
    }
    out
  }

  predict_edge_j_values <- function(pr, req_ids) {
    if (!length(req_ids)) return(numeric(0))
    out <- numeric(length(req_ids))
    for (start in seq.int(1L, length(req_ids), by = chunk)) {
      idx <- seq.int(start, min(length(req_ids), start + chunk - 1L))
      sub_ids <- req_ids[idx]
      out[idx] <- static_nuisance_j_predict(imputation_fit, ids = sub_ids, cell = pr$cell_key)$value
    }
    out
  }

  build_moment <- function(kind, group, pair, h = NULL, h_jacobian = list(),
                           h_coef_if = NULL, Z = NULL) {
    ed <- resolve_edge(group, pair)
    grp <- ed$grp
    pr <- ed$pr

    if (kind == "mean") {
      if (is.null(h_jacobian) || !is.list(h_jacobian) || is.object(h_jacobian)) {
        static_obs_abort(
          "Supply `h_jacobian` explicitly as a keyed list of Jacobian matrices (use list() for a known, non-estimated prediction).",
          "cdid_static_obs_channel_error"
        )
      }
      h_keys <- names(h_jacobian)
      if (length(h_jacobian)) {
        if (is.null(h_keys) || anyNA(h_keys) || any(!nzchar(h_keys)) || anyDuplicated(h_keys)) {
          static_obs_abort(
            "Every channel of `h_jacobian` needs a unique non-empty key.",
            "cdid_static_obs_channel_error"
          )
        }
        if (any(grepl("^(selection|imputation)\\|", h_keys))) {
          static_obs_abort(
            "`h_jacobian` channel keys cannot use the reserved 'selection|' or 'imputation|' namespaces.",
            "cdid_static_obs_channel_error"
          )
        }
        if (is.null(h_coef_if) || (is.list(h_coef_if) && !length(h_coef_if))) {
          static_obs_abort(
            "Non-empty `h_jacobian` requires `h_coef_if` (matrix list or provider function) so treatment nuisance IF is never omitted.",
            "cdid_static_obs_channel_error"
          )
        }
      } else {
        h_keys <- character(0)
        if (!is.null(h_coef_if) && !(is.list(h_coef_if) && !length(h_coef_if))) {
          static_obs_abort(
            "`h_coef_if` must be NULL or list() when `h_jacobian` is empty.",
            "cdid_static_obs_channel_error"
          )
        }
      }
    } else {
      h_keys <- character(0)
      if (is.null(Z)) Z <- grp$Z
      if (is.null(Z)) {
        static_obs_adapter_abort(
          sprintf("Covariance for group '%s', pair '%s' requires stayer `Z`.", group, pair)
        )
      }
    }

    q_tmp <- if (uses_q) predict_edge_q_full(grp, pr) else list(q = NULL, q_jacobian = NULL, assumption = NULL)
    j_tmp <- if (uses_j) predict_edge_j_full(grp, pr) else list(j = NULL, j_jacobian = NULL)
    O_arg <- if (uses_O) pr$observed else NULL
    dy_arg <- if (uses_O) pr$dy else NULL

    res <- if (kind == "mean") {
      static_obs_mean(
        id = unit_id,
        A = grp$mask,
        h = h,
        route = route,
        O = O_arg,
        dy = dy_arg,
        q = q_tmp$q,
        j = j_tmp$j,
        h_jacobian = h_jacobian,
        q_jacobian = q_tmp$q_jacobian,
        j_jacobian = j_tmp$j_jacobian,
        coef_if = NULL
      )
    } else {
      static_obs_cov(
        id = unit_id,
        A = grp$mask,
        Z = Z,
        route = route,
        O = O_arg,
        dy = dy_arg,
        q = q_tmp$q,
        j = j_tmp$j,
        q_jacobian = q_tmp$q_jacobian,
        j_jacobian = j_tmp$j_jacobian,
        coef_if = NULL
      )
    }

    # Freeze compact full-N scalar/vector summaries and loadings, then discard
    # full-edge temporary predictions, Jacobians and dense direct IF.
    est <- res$estimate
    S <- res$denominator_sum
    S_mean <- res$denominator_mean
    M <- res$pseudo_mean
    Zbar <- res$weighted_design_mean
    loadings <- res$loadings
    n_weight_rows <- res$n_weight_rows
    n_observed_rows <- res$n_observed_rows
    components <- res$components
    kink <- res$kink
    assumption <- q_tmp$assumption

    # Retain only weight-row slices of external inputs (h for mean, Z for cov)
    # so the closure never holds unnecessary full-N temporary arrays.
    is_w_mask <- if (route %in% c("baseline", "weighting")) (grp$mask & pr$observed) else grp$mask
    R_w <- which(is_w_mask)
    h_on_w <- if (kind == "mean") as.double(h[R_w]) else rep.int(1, length(R_w))
    Z_on_w <- if (kind == "cov") {
      Zm <- static_obs_matrix(Z, N, ids_chr, "Z")
      Zm[R_w, , drop = FALSE]
    } else NULL
    out_cols <- if (kind == "cov") colnames(Z_on_w) else NULL

    # Validate h_coef_if upfront on a 1-row probe if it is a static list of matrices.
    if (kind == "mean" && length(h_keys) && is.list(h_coef_if)) {
      probe <- static_obs_adapter_eval_h_if(h_coef_if, h_keys, unit_id[1L], 1L, ids_chr)
      static_obs_contract(loadings[h_keys], probe, unit_id[1L])
      rm(probe)
    }

    rm(res, q_tmp, j_tmp, O_arg, dy_arg, h, h_jacobian, Z)

    all_loading_keys <- names(loadings)
    if (is.null(all_loading_keys)) all_loading_keys <- character(0)
    sel_loading_keys <- all_loading_keys[grepl("^selection\\|", all_loading_keys)]
    imp_loading_keys <- all_loading_keys[grepl("^imputation\\|", all_loading_keys)]

    eval_if <- function(rows = unit_id, include_components = FALSE) {
      rk <- static_obs_adapter_rows(unit_id, rows)
      req_ids <- rk$ids
      req_pos <- rk$pos
      req_chr <- as.character(req_ids)
      n_req <- length(req_pos)

      # 1. Replay direct IF on requested rows from frozen full-N summary (N, S, M, est, Zbar).
      w_match <- match(req_pos, R_w)
      sel_w <- which(!is.na(w_match))
      direct <- if (kind == "mean") {
        stats::setNames(numeric(n_req), req_chr)
      } else {
        matrix(0, nrow = n_req, ncol = ncol(Z_on_w), dimnames = list(req_chr, out_cols))
      }

      if (length(sel_w)) {
        idx_w <- w_match[sel_w]
        pos_w <- req_pos[sel_w]
        ids_w <- unit_id[pos_w]
        nw_sub <- length(sel_w)
        h_sub <- h_on_w[idx_w]

        obs_sub <- if (uses_O) pr$observed[pos_w] else rep(FALSE, nw_sub)
        o_local <- which(obs_sub)
        pos_o <- pos_w[o_local]
        ids_o <- unit_id[pos_o]
        dy_o <- if (uses_O && length(o_local)) pr$dy[pos_o] else numeric(0)

        phi_sub <- numeric(nw_sub)
        if (uses_q && length(o_local)) {
          q_o <- predict_edge_q_values(pr, ids_o)
          phi_sub[o_local] <- static_obs_floor_inverse(q_o)$phi
        }
        j_sub <- if (uses_j) predict_edge_j_values(pr, ids_w) else NULL

        if (route == "baseline") {
          w_sub <- h_sub
          Y_sub <- dy_o
        } else if (route == "weighting") {
          w_sub <- h_sub * phi_sub
          Y_sub <- dy_o
        } else if (route == "imputation") {
          w_sub <- h_sub
          Y_sub <- j_sub
        } else {
          resid_sub <- numeric(nw_sub)
          if (length(o_local)) resid_sub[o_local] <- dy_o - j_sub[o_local]
          w_sub <- h_sub
          Y_sub <- j_sub + phi_sub * resid_sub
        }

        if (kind == "mean") {
          direct[sel_w] <- N * w_sub * (Y_sub - M) / S
        } else {
          Zw_sub <- Z_on_w[idx_w, , drop = FALSE]
          Zc_sub <- sweep(Zw_sub, 2L, Zbar, "-")
          Fw_sub <- Zc_sub * (Y_sub - M) - matrix(est, nw_sub, length(est), byrow = TRUE)
          direct[sel_w, ] <- N * w_sub * Fw_sub / S
        }
      }

      # 2. Contract fixed moment loadings against original-unit nuisance IF rows.
      nuis <- if (kind == "mean") {
        stats::setNames(numeric(n_req), req_chr)
      } else {
        matrix(0, nrow = n_req, ncol = ncol(Z_on_w), dimnames = list(req_chr, out_cols))
      }

      if (length(all_loading_keys)) {
        sel_full_if <- if (length(sel_loading_keys)) {
          get_batch_selection_if(req_ids, req_pos)
        } else NULL
        imp_full_if <- if (length(imp_loading_keys)) {
          get_batch_imputation_if(req_ids, req_pos)
        } else NULL
        sel_if <- static_obs_adapter_fit_if_channels(
          selection_fit, "selection", sel_loading_keys, req_ids, chunk,
          full_if = sel_full_if
        )
        imp_if <- static_obs_adapter_fit_if_channels(
          imputation_fit, "imputation", imp_loading_keys, req_ids, chunk,
          full_if = imp_full_if
        )
        h_if <- if (kind == "mean") {
          static_obs_adapter_eval_h_if(h_coef_if, h_keys, req_ids, req_pos, ids_chr)
        } else list()
        coef_if_req <- c(h_if, sel_if, imp_if)[all_loading_keys]
        contrib <- static_obs_contract(loadings, coef_if_req, req_ids)
        if (kind == "mean") {
          nuis <- stats::setNames(drop(contrib), req_chr)
        } else {
          nuis <- contrib
          dimnames(nuis) <- list(req_chr, out_cols)
        }
      }

      full_if <- direct + nuis
      if (isTRUE(include_components)) {
        return(list(influence = full_if, direct_if = direct, nuisance_if = nuis))
      }
      full_if
    }

    diagnostics <- list(
      n_units = N,
      n_role_rows = length(grp$inA),
      n_weight_rows = n_weight_rows,
      n_observed_rows = n_observed_rows,
      denominator_sum = S,
      denominator_mean = S_mean,
      pseudo_mean = M,
      weighted_design_mean = Zbar,
      components = components,
      kink = kink,
      selection_model = if (uses_q) selection_model else NULL,
      selection_pooling = if (uses_q) selection_pooling else NULL,
      imputation_pooling = if (uses_j) imputation_pooling else NULL,
      imputation_scale = if (uses_j) imputation_scale else NULL,
      assumption = assumption
    )

    structure(
      list(
        kind = kind,
        group = group,
        role = grp$role,
        pair = pair,
        endpoints = c(start = pr$start, end = pr$end),
        gap = pr$gap,
        route = route,
        estimate = est,
        denominator_sum = S,
        denominator_mean = S_mean,
        pseudo_mean = M,
        weighted_design_mean = Zbar,
        loadings = loadings,
        diagnostics = diagnostics,
        kink = kink,
        components = components,
        eval_if = eval_if
      ),
      class = "cdid_static_observation_moment"
    )
  }

  edge_mean <- function(group, pair, h, h_jacobian = list(), h_coef_if = NULL) {
    build_moment(
      "mean", group = group, pair = pair, h = h,
      h_jacobian = h_jacobian, h_coef_if = h_coef_if, Z = NULL
    )
  }

  edge_cov <- function(group, pair, Z = NULL) {
    build_moment("cov", group = group, pair = pair, Z = Z)
  }

  nuisance_if <- function(rows = unit_id) {
    rk <- static_obs_adapter_rows(unit_id, rows)
    req_ids <- rk$ids
    sel_ch <- if (!is.null(selection_fit)) {
      paste0("selection|", names(selection_fit$layout$blocks))
    } else character(0)
    imp_ch <- if (!is.null(imputation_fit)) {
      paste0("imputation|", names(imputation_fit$layout$blocks))
    } else character(0)
    c(
      static_obs_adapter_fit_if_channels(selection_fit, "selection", sel_ch, req_ids, chunk),
      static_obs_adapter_fit_if_channels(imputation_fit, "imputation", imp_ch, req_ids, chunk)
    )
  }

  group_summary <- lapply(parsed_groups, function(g) {
    list(
      group = g$group,
      role = g$role,
      target = g$target,
      prior_state = g$prior_state,
      n_role = length(g$inA),
      pairs = lapply(g$pairs, function(p) {
        list(pair = p$pair, cell_key = p$cell_key, start = p$start, end = p$end,
             gap = p$gap, date_cell_keys = p$date_cell_keys)
      })
    )
  })

  structure(
    list(
      unit_id = unit_id,
      N = N,
      route = route,
      selection_model = if (uses_q) selection_model else NULL,
      selection_pooling = if (uses_q) selection_pooling else NULL,
      imputation_pooling = if (uses_j) imputation_pooling else NULL,
      imputation_scale = if (uses_j) imputation_scale else NULL,
      max_chunk_rows = chunk,
      selection_fit = selection_fit,
      imputation_fit = imputation_fit,
      selection_spec = selection_spec,
      imputation_spec = imputation_spec,
      groups = group_summary,
      mean = edge_mean,
      cov = edge_cov,
      nuisance_if = nuisance_if,
      with_if_batch = with_if_batch
    ),
    class = "cdid_static_observation_adapter"
  )
}

static_obs_adapter_fit <- static_observation_adapter_fit

static_observation_adapter_check <- function(adapter) {
  if (!inherits(adapter, "cdid_static_observation_adapter")) {
    static_obs_adapter_abort("Expected a `cdid_static_observation_adapter` object.")
  }
  invisible(adapter)
}

static_observation_adapter_mean <- function(adapter, group, pair, h,
                                            h_jacobian = list(), h_coef_if = NULL) {
  static_observation_adapter_check(adapter)
  adapter$mean(group = group, pair = pair, h = h,
               h_jacobian = h_jacobian, h_coef_if = h_coef_if)
}

static_obs_adapter_mean <- static_observation_adapter_mean

static_observation_adapter_cov <- function(adapter, group, pair, Z = NULL) {
  static_observation_adapter_check(adapter)
  adapter$cov(group = group, pair = pair, Z = Z)
}

static_obs_adapter_cov <- static_observation_adapter_cov
