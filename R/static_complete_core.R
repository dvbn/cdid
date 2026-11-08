# Private complete-data static binary switcher/stayer mathematical kernel.
#
# Adapts the static no-carryover binary specialization of de Chaisemartin,
# D'Haultfoeuille, Pasquier, Sow and Vazquez-Bare (2026, arXiv:2201.06898v7)
# with the chained estimation and inference framework of Bellego, Benatia and
# Dortet-Bernadet (2025, doi:10.1016/j.jeconom.2024.105783):
#   - Setup and targets: transition cells c = (t, d), masks R_c = G_c + C_c,
#     natural cell contrast theta_c, sign-oriented cell contrast
#     tau_c = sigma_c * theta_c with sigma_c = 1 - 2 d(c),
#     full switcher-frequency aggregate delta and subaggregates delta_{A,s}.
#   - Estimators: cell-specific risk-set logit working model, clipped odds
#     r_tilde(u) = exp(min(u, kbar)) with ceiling 1 - 1e-6, stayer-only
#     least-squares outcome regression with cell-specific intercepts alpha_c
#     and block-pooled slopes b_{P(c)} across a declared partition P of C,
#     separately normalized cell estimators theta_hat^{ipw,or,dr}_c, and
#     switch-share weighted aggregate delta_hat^m_{A,s}.
#   - Asymptotics and inference: first-step influence functions ell^gamma_c
#     and ell^beta, complete original-unit cell influence functions
#     phi^{dr,ipw,or}_c, aggregate influence function phi^m_{A,s} including
#     random switch-share derivatives, original-unit variance estimation,
#     and inference-only cluster aggregation.
#
# Internal contract:
#   - Private kernel functions only; no public wrapper, no formula parser, no
#     automatic cell discovery from observed switch counts, no silent dropping
#     of unsupported cells, no default outcome-regression pooling partition,
#     no default method, no default glm_control or rank_tol, and no user
#     survey-weight option.
#   - All declared cells must share the exact same n original units in the
#     exact same unit_id order.
#   - Output influence functions have n rows corresponding to original units,
#     and output covariances aggregate over original units (or supplied
#     inference-only clusters of original units), never over switch rows.

static_complete_abort <- function(message, class = character(), extra = list()) {
  cond <- c(
    list(message = message, call = NULL),
    extra
  )
  stop(structure(
    cond,
    class = unique(c(class, "cdid_static_complete_error", "cdid_error", "error", "condition"))
  ))
}

static_complete_scalar_logical <- function(x) {
  is.logical(x) && is.null(dim(x)) && !is.object(x) && length(x) == 1L && !is.na(x)
}

static_complete_validate_glm_control <- function(glm_control) {
  if (missing(glm_control) || is.null(glm_control) || !is.list(glm_control) || is.object(glm_control)) {
    static_complete_abort(
      "Supply an explicit `glm_control` list (e.g. `list(epsilon = 1e-12, maxit = 50L, trace = FALSE)`); no default solver policy is adopted.",
      c("cdid_static_control_error", "cdid_ipw_control_error", "cdid_ipw_input_error", "cdid_static_input_error")
    )
  }
  ctrl <- glm_control
  nm <- names(ctrl)
  if (!is.null(nm) && length(ctrl) == 2L && !anyDuplicated(nm) && setequal(nm, c("epsilon", "maxit"))) {
    ctrl$trace <- FALSE
  }
  validate_ipw_glm_control(ctrl)
}

static_complete_validate_rank_tol <- function(rank_tol) {
  if (missing(rank_tol) || !ipw_scalar_numeric(rank_tol, finite = TRUE) ||
      rank_tol <= 0 || rank_tol >= 1) {
    static_complete_abort(
      "Supply an explicit finite `rank_tol` in `(0, 1)` (e.g. `1e-10`); no default rank tolerance is adopted.",
      c("cdid_static_rank_tol_error", "cdid_or_wls_input_error", "cdid_static_input_error")
    )
  }
  as.double(rank_tol)
}

static_complete_validate_method <- function(method) {
  if (missing(method) || !is.character(method) || !is.null(dim(method)) ||
      is.object(method) || length(method) != 1L || is.na(method) ||
      !(method %in% c("ipw", "reg", "or", "dr"))) {
    static_complete_abort(
      "Supply an explicit `method` in `c(\"ipw\", \"reg\", \"dr\")` (or paper alias `\"or\"`); no default estimator is adopted.",
      c("cdid_static_method_error", "cdid_static_input_error")
    )
  }
  canonical <- if (identical(method, "or")) "reg" else method
  paper <- if (identical(canonical, "reg")) "or" else canonical
  list(requested = method, method = canonical, paper_method = paper)
}

static_complete_validate_unit_id <- function(unit_id, label = "unit_id") {
  if (missing(unit_id) || is.null(unit_id)) {
    static_complete_abort(
      sprintf("Explicit `%s` vector is required to verify original-unit identity and order.", label),
      c("cdid_static_unit_id_error", "cdid_ipw_id_error", "cdid_static_input_error")
    )
  }
  if (is.factor(unit_id)) unit_id <- as.character(unit_id)
  if (!is.atomic(unit_id) || !is.null(dim(unit_id)) || is.object(unit_id) ||
      length(unit_id) < 1L || anyNA(unit_id)) {
    static_complete_abort(
      sprintf("`%s` must be a nonmissing atomic vector of original-unit identifiers.", label),
      c("cdid_static_unit_id_error", "cdid_ipw_id_error", "cdid_static_input_error")
    )
  }
  is_num <- is.numeric(unit_id)
  is_chr <- is.character(unit_id)
  if (!(is_num || is_chr) ||
      (is_num && any(!is.finite(unit_id))) ||
      (is_chr && any(!nzchar(unit_id)))) {
    static_complete_abort(
      sprintf("`%s` must be finite numeric or nonempty character identifiers.", label),
      c("cdid_static_unit_id_error", "cdid_ipw_id_error", "cdid_static_input_error")
    )
  }
  key <- if (is_num) as.double(unit_id) else unit_id
  if (anyDuplicated(key)) {
    static_complete_abort(
      sprintf("`%s` must not contain duplicate original-unit identifiers.", label),
      c("cdid_static_unit_id_error", "cdid_ipw_id_error", "cdid_static_input_error")
    )
  }
  list(id = unit_id, key = key, is_num = is_num, n = length(unit_id))
}

static_complete_same_unit_id <- function(u1, u2) {
  identical(u1$is_num, u2$is_num) &&
    identical(u1$n, u2$n) &&
    identical(u1$key, u2$key)
}

static_complete_validate_binary_mask <- function(mask, n, name, cell_id) {
  if (missing(mask) || is.null(mask) || !is.atomic(mask) || !is.null(dim(mask)) ||
      is.object(mask) || length(mask) != n || anyNA(mask)) {
    static_complete_abort(
      sprintf("Cell `%s` mask `%s` must be a nonmissing binary vector of length `n = %d`.", cell_id, name, n),
      c("cdid_static_mask_error", "cdid_static_input_error")
    )
  }
  if (is.logical(mask)) {
    return(as.integer(mask))
  }
  if (!is.numeric(mask) || any(!is.finite(mask) | !(mask %in% c(0, 1)))) {
    static_complete_abort(
      sprintf("Cell `%s` mask `%s` must contain only binary values in `{0, 1}`.", cell_id, name),
      c("cdid_static_mask_error", "cdid_static_input_error")
    )
  }
  as.integer(mask)
}

static_complete_validate_design_pair <- function(X_tilde = NULL, Z = NULL, n, cell_id, role = "treatment") {
  has_X <- !is.null(X_tilde)
  has_Z <- !is.null(Z)
  if (!has_X && !has_Z) {
    static_complete_abort(
      sprintf("Cell `%s` requires an explicit %s covariate design (`X_tilde` with constant column 1 and/or slope matrix `Z`).", cell_id, role),
      c("cdid_static_design_error", "cdid_ipw_covariate_error", "cdid_static_input_error")
    )
  }
  if (has_X) {
    if (!is.matrix(X_tilde) || !is.numeric(X_tilde) || is.object(X_tilde) ||
        nrow(X_tilde) != n || ncol(X_tilde) < 1L || anyNA(X_tilde) || any(!is.finite(X_tilde))) {
      static_complete_abort(
        sprintf("Cell `%s` %s `X_tilde` must be a finite numeric matrix with `nrow == %d` and `ncol >= 1`.", cell_id, role, n),
        c("cdid_static_design_error", "cdid_ipw_covariate_error", "cdid_static_input_error")
      )
    }
    if (!all(X_tilde[, 1L] == 1)) {
      static_complete_abort(
        sprintf("Cell `%s` %s `X_tilde` must include a leading constant column of ones (`X_tilde[, 1] == 1`).", cell_id, role),
        c("cdid_static_design_error", "cdid_ipw_covariate_error", "cdid_static_input_error")
      )
    }
    storage.mode(X_tilde) <- "double"
    rownames(X_tilde) <- NULL
  }
  if (has_Z) {
    if (!is.matrix(Z) || !is.numeric(Z) || is.object(Z) ||
        nrow(Z) != n || anyNA(Z) || any(!is.finite(Z))) {
      static_complete_abort(
        sprintf("Cell `%s` %s `Z` must be a finite numeric matrix with `nrow == %d`.", cell_id, role, n),
        c("cdid_static_design_error", "cdid_ipw_covariate_error", "cdid_static_input_error")
      )
    }
    storage.mode(Z) <- "double"
    rownames(Z) <- NULL
  }
  if (has_X && has_Z) {
    if (ncol(X_tilde) != ncol(Z) + 1L ||
        !all(X_tilde[, -1L, drop = FALSE] == Z)) {
      static_complete_abort(
        sprintf("Cell `%s` %s `X_tilde` and `Z` are inconsistent; require `X_tilde == cbind(1, Z)`.", cell_id, role),
        c("cdid_static_design_error", "cdid_ipw_covariate_error", "cdid_static_input_error")
      )
    }
  } else if (has_X) {
    Z <- X_tilde[, -1L, drop = FALSE]
  } else {
    X_tilde <- cbind("(Intercept)" = rep(1.0, n), Z)
  }

  p <- ncol(Z)
  z_names <- colnames(Z)
  if (is.null(z_names)) {
    if (p > 0L) {
      x_names <- colnames(X_tilde)
      if (!is.null(x_names) && length(x_names) == p + 1L && all(!is.na(x_names[-1L])) && all(nzchar(x_names[-1L]))) {
        z_names <- x_names[-1L]
      } else {
        z_names <- paste0("z", seq_len(p))
      }
    } else {
      z_names <- character(0L)
    }
  } else {
    if (anyNA(z_names) || any(!nzchar(z_names)) || anyDuplicated(z_names)) {
      static_complete_abort(
        sprintf("Cell `%s` %s slope column names must be non-empty and unique.", cell_id, role),
        c("cdid_static_design_error", "cdid_static_input_error")
      )
    }
  }
  if ( "(Intercept)" %in% z_names ) {
    static_complete_abort(
      sprintf("Cell `%s` %s slope matrix `Z` must not contain an `(Intercept)` column name.", cell_id, role),
      c("cdid_static_design_error", "cdid_static_input_error")
    )
  }
  colnames(Z) <- z_names
  colnames(X_tilde) <- c("(Intercept)", z_names)
  list(X_tilde = X_tilde, Z = Z, k = ncol(X_tilde), p = p, z_names = z_names)
}

# Construct a single declared cell specification for the static complete-data kernel.
static_complete_cell <- function(cell_id,
                                 d = NULL,
                                 sigma = NULL,
                                 t = NULL,
                                 unit_id,
                                 R,
                                 G,
                                 C,
                                 dy,
                                 X_tilde = NULL,
                                 Z = NULL,
                                 X_or_tilde = NULL,
                                 Z_or = NULL,
                                 or_block = NULL,
                                 allow_empty_switchers = FALSE) {
  if (!static_complete_scalar_logical(allow_empty_switchers)) {
    static_complete_abort(
      "`allow_empty_switchers` must be a non-missing logical scalar.",
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }
  if (missing(cell_id) || !is.character(cell_id) || !is.null(dim(cell_id)) ||
      is.object(cell_id) || length(cell_id) != 1L || is.na(cell_id) || !nzchar(cell_id)) {
    static_complete_abort(
      "`cell_id` must be a non-empty character scalar.",
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }
  uid <- static_complete_validate_unit_id(unit_id, label = sprintf("cells[[%s]]$unit_id", cell_id))
  n <- uid$n

  has_d <- !is.null(d)
  has_sigma <- !is.null(sigma)
  if (!has_d && !has_sigma) {
    static_complete_abort(
      sprintf("Cell `%s` must explicitly supply initial treatment status `d` in `{0, 1}` and/or switcher sign `sigma` in `{-1, 1}`.", cell_id),
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }
  if (has_d) {
    if (!ipw_scalar_numeric(d, finite = TRUE) || !(d %in% c(0, 1))) {
      static_complete_abort(
        sprintf("Cell `%s` field `d` must be a binary scalar in `{0, 1}`.", cell_id),
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
    d_val <- as.integer(d)
  }
  if (has_sigma) {
    if (!ipw_scalar_numeric(sigma, finite = TRUE) || !(sigma %in% c(-1, 1))) {
      static_complete_abort(
        sprintf("Cell `%s` field `sigma` must be a sign scalar in `{-1, 1}`.", cell_id),
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
    sigma_val <- as.double(sigma)
  }
  if (has_d && has_sigma) {
    if (sigma_val != 1 - 2 * d_val) {
      static_complete_abort(
        sprintf("Cell `%s` has inconsistent `d = %d` and `sigma = %g`; require `sigma == 1 - 2 * d`.", cell_id, d_val, sigma_val),
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
  } else if (has_d) {
    sigma_val <- as.double(1L - 2L * d_val)
  } else {
    d_val <- as.integer((1 - sigma_val) / 2)
  }

  t_val <- NULL
  if (!is.null(t)) {
    if (!ipw_scalar_numeric(t, finite = TRUE)) {
      static_complete_abort(
        sprintf("Cell `%s` field `t` must be a finite numeric scalar when supplied.", cell_id),
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
    t_val <- as.double(t)
  }

  R_int <- static_complete_validate_binary_mask(R, n, "R", cell_id)
  G_int <- static_complete_validate_binary_mask(G, n, "G", cell_id)
  C_int <- static_complete_validate_binary_mask(C, n, "C", cell_id)
  if (any(G_int * C_int != 0L) || any(R_int != (G_int + C_int))) {
    static_complete_abort(
      sprintf("Cell `%s` masks violate `R_c = G_c + C_c` with disjoint switcher `G_c` and stayer `C_c`.", cell_id),
      c("cdid_static_mask_error", "cdid_static_input_error")
    )
  }

  n_risk <- sum(R_int)
  n_switchers <- sum(G_int)
  n_stayers <- sum(C_int)
  if (n_risk < 1L || (!isTRUE(allow_empty_switchers) && n_switchers < 1L) || n_stayers < 1L) {
    static_complete_abort(
      sprintf(
        "Cell `%s` has insufficient empirical sample support (`n_risk = %d`, `n_switchers = %d`, `n_stayers = %d`): availability event E_n fails in this sample (sample-support obstruction, not population non-identification); unsupported cells are never silently dropped.",
        cell_id, n_risk, n_switchers, n_stayers
      ),
      c("cdid_static_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error"),
      extra = list(
        cell_id = cell_id,
        n_risk = n_risk,
        n_switchers = n_switchers,
        n_stayers = n_stayers,
        sample_support_only = TRUE,
        identification_claim = "none: empty sample support in a requested cell violates sample event E_n and does not alter or re-define the population target"
      )
    )
  }

  if (missing(dy) || is.null(dy) || !is.atomic(dy) || !is.null(dim(dy)) ||
      is.object(dy) || !is.numeric(dy) || length(dy) != n ||
      anyNA(dy) || any(!is.finite(dy))) {
    static_complete_abort(
      sprintf("Cell `%s` outcome change `dy` must be a finite numeric vector of length `n = %d`.", cell_id, n),
      c("cdid_static_outcome_error", "cdid_ipw_outcome_error", "cdid_static_input_error")
    )
  }
  dy_num <- as.double(dy)

  ipw_des <- static_complete_validate_design_pair(
    X_tilde = X_tilde, Z = Z, n = n, cell_id = cell_id, role = "treatment"
  )
  if (!is.null(X_or_tilde) || !is.null(Z_or)) {
    or_des <- static_complete_validate_design_pair(
      X_tilde = X_or_tilde, Z = Z_or, n = n, cell_id = cell_id, role = "outcome"
    )
    separate_or_design <- TRUE
  } else {
    or_des <- ipw_des
    separate_or_design <- FALSE
  }

  block_val <- NULL
  if (!is.null(or_block)) {
    if (is.factor(or_block)) or_block <- as.character(or_block)
    if (!is.atomic(or_block) || !is.null(dim(or_block)) || is.object(or_block) ||
        length(or_block) != 1L || is.na(or_block) ||
        (is.numeric(or_block) && !is.finite(or_block)) ||
        (is.character(or_block) && !nzchar(or_block))) {
      static_complete_abort(
        sprintf("Cell `%s` field `or_block` must be a non-empty character or finite numeric scalar.", cell_id),
        c("cdid_static_partition_error", "cdid_static_input_error")
      )
    }
    block_val <- as.character(or_block)
  }

  structure(
    list(
      cell_id = cell_id,
      t = t_val,
      d = d_val,
      sigma = sigma_val,
      n = n,
      unit_id = uid$id,
      unit_key = uid$key,
      unit_id_is_num = uid$is_num,
      R = R_int,
      G = G_int,
      C = C_int,
      n_risk = n_risk,
      n_switchers = n_switchers,
      n_stayers = n_stayers,
      dy = dy_num,
      X_tilde = ipw_des$X_tilde,
      Z = ipw_des$Z,
      k_ipw = ipw_des$k,
      p_ipw = ipw_des$p,
      z_names_ipw = ipw_des$z_names,
      X_or_tilde = or_des$X_tilde,
      Z_or = or_des$Z,
      p_or = or_des$p,
      z_names_or = or_des$z_names,
      separate_or_design = separate_or_design,
      or_block = block_val
    ),
    class = "cdid_static_complete_cell"
  )
}

static_complete_normalize_cell <- function(raw_cell,
                                           fallback_id = NULL,
                                           top_unit_id = NULL,
                                           allow_empty_switchers = FALSE) {
  if (inherits(raw_cell, "cdid_static_complete_cell")) {
    if (!is.null(fallback_id) && nzchar(fallback_id) && !identical(raw_cell$cell_id, fallback_id)) {
      static_complete_abort(
        sprintf("List name `%s` does not match `cell$cell_id` (`%s`).", fallback_id, raw_cell$cell_id),
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
    if (!isTRUE(allow_empty_switchers) && raw_cell$n_switchers < 1L) {
      static_complete_abort(
        sprintf(
          "Cell `%s` has insufficient empirical sample support (`n_risk = %d`, `n_switchers = %d`, `n_stayers = %d`): availability event E_n fails in this sample (sample-support obstruction, not population non-identification); unsupported cells are never silently dropped.",
          raw_cell$cell_id, raw_cell$n_risk, raw_cell$n_switchers, raw_cell$n_stayers
        ),
        c("cdid_static_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error")
      )
    }
    return(raw_cell)
  }
  if (!is.list(raw_cell) || is.object(raw_cell)) {
    static_complete_abort(
      "Every element of `cells` must be a `cdid_static_complete_cell` or a plain named list.",
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }
  cid <- raw_cell$cell_id
  if (is.null(cid)) {
    cid <- fallback_id
  } else if (!is.null(fallback_id) && nzchar(fallback_id) && !identical(cid, fallback_id)) {
    static_complete_abort(
      sprintf("List name `%s` does not match element `cell_id` (`%s`).", fallback_id, cid),
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }
  uid <- if (!is.null(raw_cell$unit_id)) raw_cell$unit_id else top_unit_id
  if (is.null(uid)) {
    static_complete_abort(
      sprintf("Cell `%s` is missing `unit_id`; original-unit identity and order must be explicitly supplied.", if (is.null(cid)) "<unnamed>" else cid),
      c("cdid_static_unit_id_error", "cdid_static_input_error")
    )
  }
  static_complete_cell(
    cell_id = cid,
    d = raw_cell$d,
    sigma = raw_cell$sigma,
    t = raw_cell$t,
    unit_id = uid,
    R = raw_cell$R,
    G = raw_cell$G,
    C = raw_cell$C,
    dy = raw_cell$dy,
    X_tilde = raw_cell$X_tilde,
    Z = raw_cell$Z,
    X_or_tilde = raw_cell$X_or_tilde,
    Z_or = raw_cell$Z_or,
    or_block = raw_cell$or_block,
    allow_empty_switchers = allow_empty_switchers
  )
}

static_complete_validate_partition <- function(cells, or_blocks = NULL) {
  cell_ids <- names(cells)
  n_cells <- length(cell_ids)
  per_cell_blocks <- vapply(cells, function(cl) {
    if (is.null(cl$or_block)) NA_character_ else as.character(cl$or_block)
  }, character(1L))
  has_per_cell <- !all(is.na(per_cell_blocks))
  if (has_per_cell && anyNA(per_cell_blocks)) {
    static_complete_abort(
      "When `or_block` is supplied on cells, every declared cell must specify `or_block`.",
      c("cdid_static_partition_error", "cdid_static_input_error")
    )
  }

  has_arg <- !is.null(or_blocks)
  if (!has_per_cell && !has_arg) {
    static_complete_abort(
      "Supply an explicit outcome-regression block partition (`or_blocks` or per-cell `or_block`); no pooling default is adopted.",
      c("cdid_static_partition_error", "cdid_static_input_error")
    )
  }

  parse_arg_blocks <- function(ob) {
    if (is.factor(ob)) ob <- as.character(ob)
    if (is.atomic(ob) && is.null(dim(ob)) && !is.object(ob)) {
      if (length(ob) != n_cells || anyNA(ob) ||
          (is.numeric(ob) && any(!is.finite(ob))) ||
          (is.character(ob) && any(!nzchar(ob)))) {
        static_complete_abort(
          sprintf("`or_blocks` vector must have length `%d` with non-empty finite block labels for every declared cell.", n_cells),
          c("cdid_static_partition_error", "cdid_static_input_error")
        )
      }
      ob_chr <- as.character(ob)
      nm <- names(ob)
      if (!is.null(nm)) {
        if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) || !setequal(nm, cell_ids)) {
          static_complete_abort(
            "Named `or_blocks` vector must have unique names matching the declared `cell_ids`.",
            c("cdid_static_partition_error", "cdid_static_input_error")
          )
        }
        ob_chr <- ob_chr[match(cell_ids, nm)]
      }
      names(ob_chr) <- cell_ids
      return(ob_chr)
    }
    if (is.list(ob) && !is.object(ob) && length(ob) >= 1L) {
      list_names <- names(ob)
      if (!is.null(list_names) && (anyNA(list_names) || any(!nzchar(list_names)) || anyDuplicated(list_names))) {
        static_complete_abort(
          "Named `or_blocks` partition list must have non-empty unique block names.",
          c("cdid_static_partition_error", "cdid_static_input_error")
        )
      }
      seen <- character(0L)
      cell_to_block <- stats::setNames(rep(NA_character_, n_cells), cell_ids)
      for (b_idx in seq_along(ob)) {
        members <- ob[[b_idx]]
        if (is.factor(members)) members <- as.character(members)
        if (!is.atomic(members) || !is.null(dim(members)) || is.object(members) ||
            length(members) < 1L || anyNA(members)) {
          static_complete_abort(
            "`or_blocks` partition list elements must be non-empty vectors of declared cell IDs or indices.",
            c("cdid_static_partition_error", "cdid_static_input_error")
          )
        }
        if (is.numeric(members)) {
          if (any(!is.finite(members) | members != floor(members) | members < 1 | members > n_cells)) {
            static_complete_abort(
              "`or_blocks` numeric cell indices must be integers in `1..length(cells)`.",
              c("cdid_static_partition_error", "cdid_static_input_error")
            )
          }
          member_ids <- cell_ids[as.integer(members)]
        } else if (is.character(members)) {
          if (any(!nzchar(members)) || !all(members %in% cell_ids)) {
            static_complete_abort(
              "`or_blocks` character members must all be declared `cell_ids`.",
              c("cdid_static_partition_error", "cdid_static_input_error")
            )
          }
          member_ids <- members
        } else {
          static_complete_abort(
            "`or_blocks` partition list elements must be character cell IDs or integer cell indices.",
            c("cdid_static_partition_error", "cdid_static_input_error")
          )
        }
        if (anyDuplicated(member_ids) || any(member_ids %in% seen)) {
          static_complete_abort(
            "`or_blocks` partition list must be disjoint across cells.",
            c("cdid_static_partition_error", "cdid_static_input_error")
          )
        }
        b_label <- if (!is.null(list_names)) list_names[b_idx] else paste0("block_", b_idx)
        seen <- c(seen, member_ids)
        cell_to_block[member_ids] <- b_label
      }
      if (!setequal(seen, cell_ids)) {
        static_complete_abort(
          "`or_blocks` partition list must cover every declared cell in `cells`.",
          c("cdid_static_partition_error", "cdid_static_input_error")
        )
      }
      return(cell_to_block)
    }
    static_complete_abort(
      "`or_blocks` must be an explicit per-cell block vector or a disjoint list partitioning the declared cells.",
      c("cdid_static_partition_error", "cdid_static_input_error")
    )
  }

  cell_block <- if (has_arg) parse_arg_blocks(or_blocks) else stats::setNames(per_cell_blocks, cell_ids)
  if (has_arg && has_per_cell) {
    # Require exact equivalence of the partition induced by `or_blocks` and `cell$or_block`.
    if (!identical(unname(cell_block), unname(per_cell_blocks))) {
      static_complete_abort(
        "Supplied `or_blocks` argument conflicts with per-cell `or_block` fields.",
        c("cdid_static_partition_error", "cdid_static_input_error")
      )
    }
  }

  block_ids <- unique(unname(cell_block))
  blocks <- lapply(block_ids, function(bid) {
    cell_ids[cell_block == bid]
  })
  names(blocks) <- block_ids

  # Verify outcome-regression covariate dimension and column names agree within each block P.
  block_p <- integer(length(block_ids))
  names(block_p) <- block_ids
  block_colnames <- vector("list", length(block_ids))
  names(block_colnames) <- block_ids
  for (bid in block_ids) {
    cids <- blocks[[bid]]
    c1 <- cells[[cids[1L]]]
    p_ref <- c1$p_or
    names_ref <- c1$z_names_or
    for (cid in cids[-1L]) {
      ck <- cells[[cid]]
      if (ck$p_or != p_ref || !identical(ck$z_names_or, names_ref)) {
        static_complete_abort(
          sprintf(
            "Outcome-regression block `%s` has incompatible slope designs between cells `%s` and `%s`.",
            bid, cids[1L], cid
          ),
          c("cdid_static_partition_error", "cdid_static_design_error", "cdid_or_wls_input_error")
        )
      }
    }
    block_p[[bid]] <- p_ref
    block_colnames[[bid]] <- names_ref
  }

  list(
    cell_block = cell_block,
    block_ids = block_ids,
    blocks = blocks,
    n_blocks = length(block_ids),
    block_p = block_p,
    block_colnames = block_colnames
  )
}

static_complete_validate_cluster <- function(cluster, n) {
  if (is.null(cluster)) return(NULL)
  if (is.factor(cluster)) cluster <- as.character(cluster)
  if (!is.atomic(cluster) || !is.null(dim(cluster)) || is.object(cluster) ||
      length(cluster) != n || anyNA(cluster)) {
    static_complete_abort(
      sprintf("`cluster` must be NULL or a nonmissing atomic vector of length `n = %d`.", n),
      c("cdid_static_cluster_error", "cdid_static_input_error")
    )
  }
  is_num <- is.numeric(cluster)
  is_chr <- is.character(cluster)
  if (!(is_num || is_chr) ||
      (is_num && any(!is.finite(cluster))) ||
      (is_chr && any(!nzchar(cluster)))) {
    static_complete_abort(
      "`cluster` identifiers must be finite numeric or non-empty character values.",
      c("cdid_static_cluster_error", "cdid_static_input_error")
    )
  }
  key <- if (is_num) as.double(cluster) else cluster
  unique_keys <- unique(key)
  cluster_index <- match(key, unique_keys)
  n_clusters <- length(unique_keys)
  cluster_sizes <- tabulate(cluster_index, nbins = n_clusters)
  list(
    cluster = cluster,
    cluster_index = cluster_index,
    unique_clusters = unique(cluster),
    n_clusters = n_clusters,
    cluster_sizes = cluster_sizes
  )
}

static_complete_validate_inputs <- function(cells,
                                            or_blocks = NULL,
                                            method,
                                            glm_control,
                                            rank_tol,
                                            unit_id = NULL,
                                            cluster = NULL,
                                            non_target_or_cells = NULL) {
  m_info <- static_complete_validate_method(method)
  ctrl <- static_complete_validate_glm_control(glm_control)
  rtol <- static_complete_validate_rank_tol(rank_tol)

  if (missing(cells) || !is.list(cells) || is.object(cells) || length(cells) < 1L) {
    static_complete_abort(
      "Supply a non-empty list of declared `cells`; population cells are never inferred from observed switch counts.",
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }

  top_uid <- if (!is.null(unit_id)) static_complete_validate_unit_id(unit_id, label = "unit_id") else NULL
  raw_names <- names(cells)
  if (!is.null(raw_names) && (anyNA(raw_names) || anyDuplicated(raw_names[nzchar(raw_names)]))) {
    static_complete_abort(
      "Names of `cells`, when supplied, must not be NA or duplicated.",
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }

  norm_cells <- vector("list", length(cells))
  for (j in seq_along(cells)) {
    fb <- if (!is.null(raw_names) && nzchar(raw_names[j])) raw_names[j] else NULL
    if (is.null(top_uid) && !inherits(cells[[j]], "cdid_static_complete_cell") &&
        is.list(cells[[j]]) && is.null(cells[[j]]$unit_id)) {
      static_complete_abort(
        sprintf("Cell %d must supply an explicit `unit_id` vector (or supply matching `unit_id` on every cell).", j),
        c("cdid_static_unit_id_error", "cdid_static_input_error")
      )
    }
    norm_cells[[j]] <- static_complete_normalize_cell(
      cells[[j]],
      fallback_id = fb,
      top_unit_id = if (!is.null(top_uid)) top_uid$id else NULL,
      allow_empty_switchers = FALSE
    )
  }

  cell_ids <- vapply(norm_cells, `[[`, character(1L), "cell_id")
  if (anyDuplicated(cell_ids)) {
    static_complete_abort(
      "Declared `cells` must have unique `cell_id` identifiers.",
      c("cdid_static_cell_error", "cdid_static_input_error")
    )
  }
  names(norm_cells) <- cell_ids

  ref_uid <- list(
    id = norm_cells[[1L]]$unit_id,
    key = norm_cells[[1L]]$unit_key,
    is_num = norm_cells[[1L]]$unit_id_is_num,
    n = norm_cells[[1L]]$n
  )
  if (!is.null(top_uid) && !static_complete_same_unit_id(ref_uid, top_uid)) {
    static_complete_abort(
      sprintf("Cell `%s` `unit_id` does not match top-level `unit_id` in identity, type, or order.", cell_ids[1L]),
      c("cdid_static_unit_id_error", "cdid_ipw_id_error", "cdid_static_input_error")
    )
  }
  for (j in seq_along(norm_cells)[-1L]) {
    cj <- norm_cells[[j]]
    uid_j <- list(id = cj$unit_id, key = cj$unit_key, is_num = cj$unit_id_is_num, n = cj$n)
    if (!static_complete_same_unit_id(ref_uid, uid_j)) {
      static_complete_abort(
        sprintf(
          "Cell `%s` `unit_id` does not match cell `%s` `unit_id`; all declared cells must share the exact same original units in the exact same order.",
          cj$cell_id, cell_ids[1L]
        ),
        c("cdid_static_unit_id_error", "cdid_ipw_id_error", "cdid_static_input_error")
      )
    }
  }

  norm_non_target <- list()
  non_target_ids <- character(0L)
  if (!is.null(non_target_or_cells) && length(non_target_or_cells) > 0L) {
    if (!is.list(non_target_or_cells) || is.object(non_target_or_cells)) {
      static_complete_abort(
        "`non_target_or_cells` must be NULL or a list of stayer training cell specifications.",
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
    nt_names <- names(non_target_or_cells)
    if (!is.null(nt_names) && (anyNA(nt_names) || anyDuplicated(nt_names[nzchar(nt_names)]))) {
      static_complete_abort(
        "Names of `non_target_or_cells`, when supplied, must not be NA or duplicated.",
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
    norm_non_target <- vector("list", length(non_target_or_cells))
    for (j in seq_along(non_target_or_cells)) {
      fb <- if (!is.null(nt_names) && nzchar(nt_names[j])) nt_names[j] else NULL
      norm_non_target[[j]] <- static_complete_normalize_cell(
        non_target_or_cells[[j]],
        fallback_id = fb,
        top_unit_id = if (!is.null(top_uid)) top_uid$id else ref_uid$id,
        allow_empty_switchers = TRUE
      )
      cj <- norm_non_target[[j]]
      uid_j <- list(id = cj$unit_id, key = cj$unit_key, is_num = cj$unit_id_is_num, n = cj$n)
      if (!static_complete_same_unit_id(ref_uid, uid_j)) {
        static_complete_abort(
          sprintf(
            "Non-target OR cell `%s` `unit_id` does not match cell `%s` `unit_id`.",
            cj$cell_id, cell_ids[1L]
          ),
          c("cdid_static_unit_id_error", "cdid_ipw_id_error", "cdid_static_input_error")
        )
      }
    }
    non_target_ids <- vapply(norm_non_target, `[[`, character(1L), "cell_id")
    if (anyDuplicated(non_target_ids) || any(non_target_ids %in% cell_ids)) {
      static_complete_abort(
        "Declared `cells` and `non_target_or_cells` must have disjoint unique `cell_id` identifiers.",
        c("cdid_static_cell_error", "cdid_static_input_error")
      )
    }
    names(norm_non_target) <- non_target_ids
  }

  or_cells <- if (length(norm_non_target) > 0L) c(norm_cells, norm_non_target) else norm_cells
  partition <- static_complete_validate_partition(or_cells, or_blocks = or_blocks)
  cluster_info <- static_complete_validate_cluster(cluster, n = ref_uid$n)

  list(
    cells = norm_cells,
    cell_ids = cell_ids,
    non_target_cells = norm_non_target,
    non_target_ids = non_target_ids,
    or_cells = or_cells,
    n = ref_uid$n,
    unit_id = ref_uid$id,
    partition = partition,
    method_info = m_info,
    glm_control = ctrl,
    rank_tol = rtol,
    cluster_info = cluster_info
  )
}

# Fit cell-specific switching logit model on the full known risk set R_c = 1
# with clipped odds and complete original-unit first-step influence function.
# Reuses `ipw_nuisance_observer_fit_worker`, `ipw_nuisance_hessian_stage`,
# `ipw_clip_logit`, `ipw_clip_odds_active`, and `ipw_clip_summary` from
# R/chained_ipw_link.R with complete-data observation rates equal to 1.
static_complete_fit_cell_logit <- function(cell, glm_control, rank_tol, materialize_if = TRUE) {
  if (!static_complete_scalar_logical(materialize_if)) {
    static_complete_abort(
      "`materialize_if` must be a non-missing logical scalar.",
      c("cdid_static_input_error")
    )
  }
  n <- cell$n
  n_num <- as.double(n)
  risk_idx <- which(cell$R == 1L)
  risk_gi <- ifelse(cell$G[risk_idx] == 1L, 1L, 2L)
  risk_w <- rep(1.0, length(risk_idx))
  Z_use <- cell$X_tilde[risk_idx, , drop = FALSE]
  p <- ncol(Z_use)

  group_a_count <- c(as.integer(cell$n_switchers), as.integer(cell$n_stayers))
  group_mass <- as.double(group_a_count)
  sup <- list(
    cohorts = c(1, 0),
    m_cohorts = 2L,
    ci = c(1L, 2L),
    n = n_num,
    group_a_count = group_a_count,
    group_mass = group_mass,
    q = group_mass / n_num
  )

  diag_env <- new.env(parent = emptyenv())
  st <- ipw_nuisance_observer_fit_worker(
    sup = sup,
    risk_a_idx = risk_idx,
    risk_gi = risk_gi,
    risk_w = risk_w,
    Z_use = Z_use,
    ctrl = glm_control,
    diag_env = diag_env
  )

  # Apply caller's explicit `rank_tol` relative-diagonal criterion on both the
  # risk-set design QR and the final IRLS QR factor so near/exact collinearity
  # above dqrls's internal `epsilon / 1000` threshold is deterministically refused.
  q_design <- qr(Z_use, LAPACK = TRUE)
  R_design <- qr.R(q_design)
  k_diag <- min(nrow(R_design), ncol(R_design))
  diag_design <- abs(R_design[cbind(seq_len(k_diag), seq_len(k_diag))])
  cutoff_design <- rank_tol * max(diag_design)
  rank_design <- if (cutoff_design > 0) sum(diag_design > cutoff_design) else 0L
  if (rank_design < p) {
    static_complete_abort(
      sprintf(
        "Cell `%s` propensity logit design has deficient numerical rank %d < %d on the risk set (n_risk=%d, rank_tol=%g, cutoff=%g); sample rank obstruction, not population non-identification.",
        cell$cell_id, rank_design, p, nrow(Z_use), rank_tol, cutoff_design
      ),
      c("cdid_static_logit_rank_error", "cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }

  qr_diag <- abs(diag(st$fit_info$qr_R))
  cutoff_qr <- rank_tol * max(qr_diag)
  rank_qr <- if (cutoff_qr > 0) sum(qr_diag > cutoff_qr) else 0L
  if (rank_qr < p) {
    static_complete_abort(
      sprintf(
        "Cell `%s` propensity logit IRLS QR factor has deficient numerical rank %d < %d (n_risk=%d, rank_tol=%g, cutoff=%g); sample rank obstruction, not population non-identification.",
        cell$cell_id, rank_qr, p, nrow(Z_use), rank_tol, cutoff_qr
      ),
      c("cdid_static_logit_rank_error", "cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }

  hs <- ipw_nuisance_hessian_stage(st, cohort_labels = c("switcher", "stayer"))
  chol_diag <- abs(diag(hs$H_chol))
  cutoff_chol <- sqrt(rank_tol) * max(chol_diag)
  if (sum(chol_diag > cutoff_chol) < p) {
    static_complete_abort(
      sprintf(
        "Cell `%s` propensity logit Hessian factor is ill-conditioned under `rank_tol = %g`.",
        cell$cell_id, rank_tol
      ),
      c("cdid_static_logit_rank_error", "cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }

  lv_risk <- ipw_logit_values(Z_use, st$beta)
  eta_full <- rep(NA_real_, n)
  prob_full <- numeric(n)
  odds_full <- numeric(n)
  active_full <- numeric(n)
  eta_full[risk_idx] <- lv_risk$eta
  prob_full[risk_idx] <- lv_risk$prob
  odds_full[risk_idx] <- lv_risk$odds
  active_full[risk_idx] <- lv_risk$active

  # Raw score on all n units: R_{ic} (G_{ic} - Lambda(xtilde_{ic}' gamma_hat_c)) xtilde_{ic}
  score_weight <- numeric(n)
  score_weight[risk_idx] <- st$D_use - st$prob_use
  if (isTRUE(materialize_if)) {
    score_matrix <- cell$X_tilde * score_weight
    # Unit-level logit IF ell^gamma_c(W_i) = H_hat_c^{-1} score_i via triangular solves on H_chol.
    ell_gamma <- t(backsolve(hs$H_chol, forwardsolve(t(hs$H_chol), t(score_matrix))))
    colnames(ell_gamma) <- colnames(Z_use)
  } else {
    score_matrix <- NULL
    ell_gamma <- NULL
  }

  list(
    cell_id = cell$cell_id,
    gamma_hat = st$beta,
    H_hat = hs$H,
    H_chol = hs$H_chol,
    pi_hat = st$mean_a,
    Q_hat = st$mean_b,
    sum_switcher_weights = st$sum_a,
    sum_stayer_odds = st$sum_b,
    risk_idx = risk_idx,
    eta = eta_full,
    prob = prob_full,
    odds = odds_full,
    active = active_full,
    score_weight = score_weight,
    score_matrix = score_matrix,
    ell_gamma = ell_gamma,
    diagnostics = list(
      converged = st$fit_info$converged,
      boundary = st$fit_info$boundary,
      iter = st$fit_info$iter,
      rank = st$fit_info$rank,
      rank_design = rank_design,
      rank_qr = rank_qr,
      rank_tol = rank_tol,
      cutoff_design = cutoff_design,
      cutoff_qr = cutoff_qr,
      qr_R = st$fit_info$qr_R,
      qr_pivot = st$fit_info$qr_pivot,
      qr_tol = st$fit_info$qr_tol,
      df_residual = st$fit_info$df_residual,
      deviance = st$fit_info$deviance,
      null_deviance = st$fit_info$null_deviance,
      aic = st$fit_info$aic,
      rates = st$rates,
      clipping = st$fit_info$clipping
    )
  )
}

# Solve H_hat_c^{-1} rhs using the Cholesky factor H_chol of H_hat_c without matrix inversion.
static_complete_logit_hessian_solve <- function(logit_fit, rhs) {
  H_chol <- logit_fit$H_chol
  sol <- drop(backsolve(H_chol, forwardsolve(t(H_chol), rhs)))
  names(sol) <- names(logit_fit$gamma_hat)
  if (any(!is.finite(sol))) {
    static_complete_abort(
      sprintf("Nonfinite logit Hessian triangular solve in cell `%s`.", logit_fit$cell_id),
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  sol
}

# Fit stayer-only least-squares outcome regression across the declared block
# partition P of C with cell-specific intercepts and block-shared slopes,
# reusing `or_pooled_wls_block`, `or_pooled_wls_merge`,
# `or_pooled_wls_fit`, `or_pooled_wls_gram_solve`, and `or_pooled_wls_predict`
# from R/chained_or_pooled_wls.R without forming an n * |C| block design or
# inverting the normal matrix M.
static_complete_fit_or_partition <- function(cells, partition, rank_tol, n, materialize_if = TRUE) {
  if (!static_complete_scalar_logical(materialize_if)) {
    static_complete_abort(
      "`materialize_if` must be a non-missing logical scalar.",
      c("cdid_static_input_error")
    )
  }
  n_num <- as.double(n)
  cell_ids <- names(cells)
  block_ids <- partition$block_ids

  block_fits <- vector("list", length(block_ids))
  names(block_fits) <- block_ids
  block_slope_scores <- if (isTRUE(materialize_if)) {
    out_sc <- vector("list", length(block_ids))
    names(out_sc) <- block_ids
    out_sc
  } else {
    NULL
  }
  block_ell_slope <- if (isTRUE(materialize_if)) {
    out_el <- vector("list", length(block_ids))
    names(out_el) <- block_ids
    out_el
  } else {
    NULL
  }

  alpha_hat <- stats::setNames(numeric(length(cell_ids)), cell_ids)
  slope_by_cell <- vector("list", length(cell_ids))
  names(slope_by_cell) <- cell_ids
  mean_z0_by_cell <- vector("list", length(cell_ids))
  names(mean_z0_by_cell) <- cell_ids
  mean_dy0_by_cell <- stats::setNames(numeric(length(cell_ids)), cell_ids)
  fitted_by_cell <- matrix(0.0, nrow = n, ncol = length(cell_ids), dimnames = list(NULL, cell_ids))
  residuals_by_cell <- matrix(0.0, nrow = n, ncol = length(cell_ids), dimnames = list(NULL, cell_ids))
  ell_alpha_by_cell <- if (isTRUE(materialize_if)) {
    matrix(0.0, nrow = n, ncol = length(cell_ids), dimnames = list(NULL, cell_ids))
  } else {
    NULL
  }

  for (bid in block_ids) {
    cids <- partition$blocks[[bid]]
    p_b <- partition$block_p[[bid]]
    z_names_b <- partition$block_colnames[[bid]]

    cell_summaries <- lapply(cids, function(cid) {
      cell <- cells[[cid]]
      stayer_idx <- which(cell$C == 1L)
      D_stayer <- cell$Z_or[stayer_idx, , drop = FALSE]
      y_stayer <- cell$dy[stayer_idx]
      w_stayer <- rep(1.0, length(stayer_idx))
      or_pooled_wls_block(D = D_stayer, y = y_stayer, weights = w_stayer)
    })
    names(cell_summaries) <- cids

    pooled_summary <- if (length(cell_summaries) == 1L) {
      s1 <- cell_summaries[[1L]]
      names(s1$blocks) <- cids[1L]
      s1
    } else {
      or_pooled_wls_merge(cell_summaries)
    }

    bfit <- or_pooled_wls_fit(pooled_summary, rank_tol = rank_tol)
    block_fits[[bid]] <- bfit

    # Accumulate each original unit's summed stayer centered slope score across
    # all cells c' in block P:
    #   s^{b_P}_i = sum_{c' in P} C_{ic'} e_hat_{ic'} (Z_{or, ic'} - bar_Z_{0c'})
    score_b <- if (isTRUE(materialize_if)) {
      matrix(0.0, nrow = n, ncol = p_b, dimnames = list(NULL, z_names_b))
    } else {
      NULL
    }
    for (cid in cids) {
      cell <- cells[[cid]]
      b_rec <- bfit$blocks[[cid]]
      mean_z0 <- as.double(b_rec$mean_d)
      names(mean_z0) <- z_names_b
      mean_z0_by_cell[[cid]] <- mean_z0
      mean_dy0_by_cell[[cid]] <- as.double(b_rec$mean_y)
      alpha_hat[[cid]] <- as.double(bfit$intercepts[[cid]])
      slope_by_cell[[cid]] <- bfit$coefficients

      pred_c <- or_pooled_wls_predict(bfit, cell$Z_or, block = cid)
      res_c <- cell$dy - pred_c
      fitted_by_cell[, cid] <- pred_c
      residuals_by_cell[, cid] <- res_c

      if (isTRUE(materialize_if) && p_b > 0L) {
        Z_centered_c <- sweep(cell$Z_or, 2L, mean_z0, FUN = "-")
        score_b <- score_b + (cell$C * res_c) * Z_centered_c
      }
    }
    if (isTRUE(materialize_if)) {
      block_slope_scores[[bid]] <- score_b

      # Unit-level slope IF ell^{b_P}(W_i) = S_hat_P^{-1} s^{b_P}_i via triangular solve on R_P.
      ell_b <- if (p_b > 0L) {
        sol <- t(or_pooled_wls_gram_solve(bfit, t(score_b), normalizer = n_num))
        colnames(sol) <- z_names_b
        sol
      } else {
        matrix(0.0, nrow = n, ncol = 0L)
      }
      block_ell_slope[[bid]] <- ell_b

      # Unit-level local intercept IF:
      #   ell^{alpha_c}(W_i) = C_{ic} e_hat_{ic} / q0_hat_c - bar_Z_{0c}' ell^{b_P}(W_i)
      for (cid in cids) {
        cell <- cells[[cid]]
        q0_c <- cell$n_stayers / n_num
        res_c <- residuals_by_cell[, cid]
        mean_z0 <- mean_z0_by_cell[[cid]]
        slope_corr <- if (p_b > 0L) drop(ell_b %*% mean_z0) else numeric(n)
        ell_alpha_by_cell[, cid] <- (cell$C * res_c) / q0_c - slope_corr
      }
    }
  }

  slope_hat <- lapply(block_fits, `[[`, "coefficients")

  list(
    partition = partition,
    block_fits = block_fits,
    alpha_hat = alpha_hat,
    slope_hat = slope_hat,
    slope_by_cell = slope_by_cell,
    mean_z0_by_cell = mean_z0_by_cell,
    mean_dy0_by_cell = mean_dy0_by_cell,
    fitted = fitted_by_cell,
    residuals = residuals_by_cell,
    block_slope_scores = block_slope_scores,
    block_ell_slope = block_ell_slope,
    ell_alpha = ell_alpha_by_cell,
    rank_tol = rank_tol
  )
}

# Original-unit covariance and inference-only cluster covariance.
# Independent-unit covariance: crossprod(phi) / n^2.
# Cluster covariance: sum original-unit IF within each cluster before crossprod / n^2.
# No degrees-of-freedom or few-cluster corrections are applied.
static_complete_covariance <- function(phi, cluster = NULL) {
  is_vec <- is.null(dim(phi))
  if (!is.numeric(phi) || is.object(phi) || !(is_vec || is.matrix(phi)) || anyNA(phi) || any(!is.finite(phi))) {
    static_complete_abort(
      "`phi` must be a finite numeric vector or matrix of original-unit influence functions.",
      c("cdid_static_inference_error", "cdid_static_input_error")
    )
  }
  phi_mat <- if (is_vec) matrix(as.double(phi), ncol = 1L) else {
    storage.mode(phi) <- "double"
    phi
  }
  n <- nrow(phi_mat)
  q <- ncol(phi_mat)
  if (n < 1L || q < 1L) {
    static_complete_abort(
      "`phi` must have at least one original-unit row and one column.",
      c("cdid_static_inference_error", "cdid_static_input_error")
    )
  }

  c_info <- static_complete_validate_cluster(cluster, n = n)
  n_num <- as.double(n)

  if (is.null(c_info)) {
    cov_mat <- crossprod(phi_mat) / (n_num^2)
    n_clusters <- n
    clustered <- FALSE
    cluster_phi <- NULL
  } else {
    J <- c_info$n_clusters
    c_idx <- c_info$cluster_index
    if (q == 1L) {
      c_sums <- vapply(seq_len(J), function(j) sum(phi_mat[c_idx == j, 1L]), numeric(1L))
      cluster_phi <- matrix(c_sums, nrow = J, ncol = 1L, dimnames = list(NULL, colnames(phi_mat)))
    } else {
      c_sums <- vapply(
        seq_len(J),
        function(j) colSums(phi_mat[c_idx == j, , drop = FALSE]),
        numeric(q)
      )
      cluster_phi <- t(c_sums)
      colnames(cluster_phi) <- colnames(phi_mat)
    }
    cov_mat <- crossprod(cluster_phi) / (n_num^2)
    n_clusters <- J
    clustered <- TRUE
  }

  if (any(!is.finite(cov_mat))) {
    static_complete_abort(
      "Nonfinite covariance matrix encountered.",
      c("cdid_static_inference_error", "cdid_ipw_numeric_error")
    )
  }
  vhat_mat <- cov_mat * n_num
  se_vec <- sqrt(pmax(0, diag(cov_mat)))
  names(se_vec) <- colnames(phi_mat)

  list(
    cov = if (is_vec) as.double(cov_mat[1L, 1L]) else cov_mat,
    cov_matrix = cov_mat,
    vhat = if (is_vec) as.double(vhat_mat[1L, 1L]) else vhat_mat,
    vhat_matrix = vhat_mat,
    se = if (is_vec) as.double(se_vec[1L]) else se_vec,
    n = n,
    n_clusters = n_clusters,
    clustered = clustered,
    cluster_phi = cluster_phi
  )
}

# Fit the complete-data static binary switcher/stayer kernel across declared cells.
static_complete_fit <- function(cells,
                                or_blocks = NULL,
                                method,
                                glm_control,
                                rank_tol,
                                unit_id = NULL,
                                cluster = NULL,
                                aggregate_cells = NULL,
                                aggregate_signs = NULL,
                                non_target_or_cells = NULL,
                                materialize_if = TRUE) {
  if (!static_complete_scalar_logical(materialize_if)) {
    static_complete_abort(
      "`materialize_if` must be a non-missing logical scalar.",
      c("cdid_static_input_error")
    )
  }
  val <- static_complete_validate_inputs(
    cells = cells,
    or_blocks = or_blocks,
    method = method,
    glm_control = glm_control,
    rank_tol = rank_tol,
    unit_id = unit_id,
    cluster = cluster,
    non_target_or_cells = non_target_or_cells
  )

  norm_cells <- val$cells
  cell_ids <- val$cell_ids
  n_cells <- length(cell_ids)
  n <- val$n
  n_num <- as.double(n)
  partition <- val$partition
  m <- val$method_info$method
  paper_m <- val$method_info$paper_method
  ctrl <- val$glm_control
  rtol <- val$rank_tol

  use_logit <- m %in% c("ipw", "dr")
  use_or <- m %in% c("reg", "dr")

  logit_fits <- if (use_logit) {
    fits <- lapply(
      norm_cells,
      static_complete_fit_cell_logit,
      glm_control = ctrl,
      rank_tol = rtol,
      materialize_if = materialize_if
    )
    names(fits) <- cell_ids
    fits
  } else {
    NULL
  }

  or_fits <- if (use_or) {
    static_complete_fit_or_partition(
      val$or_cells,
      partition = partition,
      rank_tol = rtol,
      n = n,
      materialize_if = materialize_if
    )
  } else {
    NULL
  }

  d_vec <- stats::setNames(vapply(norm_cells, `[[`, integer(1L), "d"), cell_ids)
  sigma_vec <- stats::setNames(vapply(norm_cells, `[[`, numeric(1L), "sigma"), cell_ids)
  n_risk_vec <- stats::setNames(vapply(norm_cells, `[[`, integer(1L), "n_risk"), cell_ids)
  n_switchers_vec <- stats::setNames(vapply(norm_cells, `[[`, integer(1L), "n_switchers"), cell_ids)
  n_stayers_vec <- stats::setNames(vapply(norm_cells, `[[`, integer(1L), "n_stayers"), cell_ids)
  pi_hat_vec <- n_switchers_vec / n_num
  q0_hat_vec <- n_stayers_vec / n_num
  Q_hat_vec <- stats::setNames(rep(NA_real_, n_cells), cell_ids)

  theta_hat <- stats::setNames(numeric(n_cells), cell_ids)
  tau_hat <- stats::setNames(numeric(n_cells), cell_ids)
  mu_1_vec <- stats::setNames(numeric(n_cells), cell_ids)
  mu_0_vec <- stats::setNames(rep(NA_real_, n_cells), cell_ids)

  phi_theta <- if (isTRUE(materialize_if)) {
    matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  } else {
    NULL
  }
  phi_tau <- if (isTRUE(materialize_if)) {
    matrix(0.0, nrow = n, ncol = n_cells, dimnames = list(NULL, cell_ids))
  } else {
    NULL
  }
  cell_components <- vector("list", n_cells)
  names(cell_components) <- cell_ids

  for (cid in cell_ids) {
    cell <- norm_cells[[cid]]
    G_c <- cell$G
    C_c <- cell$C
    R_c <- cell$R
    dy_c <- cell$dy
    pi_c <- pi_hat_vec[[cid]]
    q0_c <- q0_hat_vec[[cid]]
    sigma_c <- sigma_vec[[cid]]

    if (identical(m, "ipw")) {
      lfit <- logit_fits[[cid]]
      r_c <- lfit$odds
      act_c <- lfit$active
      Q_c <- lfit$Q_hat
      Q_hat_vec[[cid]] <- Q_c

      # Normalized IPW cell contrast
      mu_1c <- sum(G_c * dy_c) / cell$n_switchers
      mu_0c <- sum(C_c * r_c * dy_c) / lfit$sum_stayer_odds
      theta_c <- mu_1c - mu_0c

      k_gamma_c <- colSums(cell$X_tilde * (C_c * act_c * r_c * (dy_c - mu_0c))) / lfit$sum_stayer_odds
      lambda_gamma_c <- static_complete_logit_hessian_solve(lfit, k_gamma_c)
      k_beta_slope_c <- NULL
      lambda_beta_slope_c <- NULL

      if (isTRUE(materialize_if)) {
        # Complete original-unit IPW influence function
        direct_switcher <- (G_c / pi_c) * (dy_c - mu_1c)
        direct_stayer <- -(C_c * r_c / Q_c) * (dy_c - mu_0c)
        logit_score_term <- -(R_c * lfit$score_weight) * drop(cell$X_tilde %*% lambda_gamma_c)
        or_score_term <- numeric(n)
        phi_c <- direct_switcher + direct_stayer + logit_score_term
      } else {
        direct_switcher <- NULL
        direct_stayer <- NULL
        logit_score_term <- NULL
        or_score_term <- NULL
        phi_c <- NULL
      }
    } else if (identical(m, "reg")) {
      bid <- partition$cell_block[[cid]]
      bfit <- or_fits$block_fits[[bid]]
      p_b <- partition$block_p[[bid]]
      e_c <- or_fits$residuals[, cid]

      # Outcome-regression cell contrast
      mu_1c <- sum(G_c * e_c) / cell$n_switchers
      mu_0c <- 0.0
      theta_c <- mu_1c

      if (p_b > 0L) {
        bar_Z_1c <- colSums(cell$Z_or * G_c) / cell$n_switchers
        bar_Z_0c <- or_fits$mean_z0_by_cell[[cid]]
        k_beta_slope_c <- bar_Z_1c - bar_Z_0c
        lambda_beta_slope_c <- or_pooled_wls_gram_solve(bfit, k_beta_slope_c, normalizer = n_num)
      } else {
        k_beta_slope_c <- numeric(0L)
        lambda_beta_slope_c <- numeric(0L)
      }
      k_gamma_c <- NULL
      lambda_gamma_c <- NULL

      if (isTRUE(materialize_if)) {
        # Complete original-unit OR influence function:
        # k^{or}_{beta,c}' ell^beta(W_i) = C_{ic} e_hat_{ic} / q0_c + (bar_Z_{1c} - bar_Z_{0c})' ell^{b_P}(W_i)
        direct_switcher <- (G_c / pi_c) * (e_c - mu_1c)
        direct_stayer <- numeric(n)
        logit_score_term <- numeric(n)
        slope_if_term <- if (p_b > 0L) {
          drop(or_fits$block_slope_scores[[bid]] %*% lambda_beta_slope_c)
        } else {
          numeric(n)
        }
        intercept_if_term <- (C_c * e_c) / q0_c
        or_score_term <- -(intercept_if_term + slope_if_term)
        phi_c <- direct_switcher + or_score_term
      } else {
        direct_switcher <- NULL
        direct_stayer <- NULL
        logit_score_term <- NULL
        or_score_term <- NULL
        phi_c <- NULL
      }
    } else {
      # method == "dr"
      lfit <- logit_fits[[cid]]
      r_c <- lfit$odds
      act_c <- lfit$active
      Q_c <- lfit$Q_hat
      Q_hat_vec[[cid]] <- Q_c

      bid <- partition$cell_block[[cid]]
      bfit <- or_fits$block_fits[[bid]]
      p_b <- partition$block_p[[bid]]
      e_c <- or_fits$residuals[, cid]

      # Doubly robust cell contrast
      mu_1c <- sum(G_c * e_c) / cell$n_switchers
      mu_0c <- sum(C_c * r_c * e_c) / lfit$sum_stayer_odds
      theta_c <- mu_1c - mu_0c

      k_gamma_c <- colSums(cell$X_tilde * (C_c * act_c * r_c * (e_c - mu_0c))) / lfit$sum_stayer_odds
      lambda_gamma_c <- static_complete_logit_hessian_solve(lfit, k_gamma_c)

      if (p_b > 0L) {
        bar_Z_1c <- colSums(cell$Z_or * G_c) / cell$n_switchers
        bar_Z_0c_r <- colSums(cell$Z_or * (C_c * r_c)) / lfit$sum_stayer_odds
        k_beta_slope_c <- bar_Z_1c - bar_Z_0c_r
        lambda_beta_slope_c <- or_pooled_wls_gram_solve(bfit, k_beta_slope_c, normalizer = n_num)
      } else {
        k_beta_slope_c <- numeric(0L)
        lambda_beta_slope_c <- numeric(0L)
      }

      if (isTRUE(materialize_if)) {
        # Complete original-unit DR influence function
        direct_switcher <- (G_c / pi_c) * (e_c - mu_1c)
        direct_stayer <- -(C_c * r_c / Q_c) * (e_c - mu_0c)
        logit_score_term <- -(R_c * lfit$score_weight) * drop(cell$X_tilde %*% lambda_gamma_c)
        or_score_term <- if (p_b > 0L) {
          -drop(or_fits$block_slope_scores[[bid]] %*% lambda_beta_slope_c)
        } else {
          numeric(n)
        }
        phi_c <- direct_switcher + direct_stayer + logit_score_term + or_score_term
      } else {
        direct_switcher <- NULL
        direct_stayer <- NULL
        logit_score_term <- NULL
        or_score_term <- NULL
        phi_c <- NULL
      }
    }

    if (!is.finite(theta_c) || (isTRUE(materialize_if) && any(!is.finite(phi_c)))) {
      static_complete_abort(
        sprintf("Nonfinite cell estimate or influence function in cell `%s`.", cid),
        c("cdid_static_numeric_error", "cdid_ipw_numeric_error")
      )
    }

    theta_hat[[cid]] <- theta_c
    tau_hat[[cid]] <- sigma_c * theta_c
    mu_1_vec[[cid]] <- mu_1c
    mu_0_vec[[cid]] <- mu_0c
    if (isTRUE(materialize_if)) {
      phi_theta[, cid] <- phi_c
      phi_tau[, cid] <- sigma_c * phi_c
    }

    cell_components[[cid]] <- list(
      direct_switcher = direct_switcher,
      direct_stayer = direct_stayer,
      logit_score = logit_score_term,
      or_score = or_score_term,
      mu_1 = mu_1c,
      mu_0 = mu_0c,
      k_gamma = k_gamma_c,
      lambda_gamma = lambda_gamma_c,
      k_beta_slope = k_beta_slope_c,
      lambda_beta_slope = lambda_beta_slope_c
    )
  }

  if (isTRUE(materialize_if)) {
    inf_theta <- static_complete_covariance(phi_theta, cluster = cluster)
    inf_tau <- static_complete_covariance(phi_tau, cluster = cluster)
    cov_theta_out <- inf_theta$cov_matrix
    vhat_theta_out <- inf_theta$vhat_matrix
    se_theta_out <- inf_theta$se
    cov_tau_out <- inf_tau$cov_matrix
    vhat_tau_out <- inf_tau$vhat_matrix
    se_tau_out <- inf_tau$se
    clustered_out <- inf_theta$clustered
    n_clusters_out <- inf_theta$n_clusters
  } else {
    cov_theta_out <- NULL
    vhat_theta_out <- NULL
    se_theta_out <- stats::setNames(rep(NA_real_, n_cells), cell_ids)
    cov_tau_out <- NULL
    vhat_tau_out <- NULL
    se_tau_out <- stats::setNames(rep(NA_real_, n_cells), cell_ids)
    clustered_out <- !is.null(val$cluster_info)
    n_clusters_out <- if (!is.null(val$cluster_info)) val$cluster_info$n_clusters else n
  }

  G_mat <- vapply(norm_cells, function(cl) as.double(cl$G), numeric(n))
  if (n == 1L) G_mat <- matrix(G_mat, nrow = 1L, ncol = n_cells)
  colnames(G_mat) <- cell_ids

  out <- structure(
    list(
      n = n,
      unit_id = val$unit_id,
      method = m,
      requested_method = val$method_info$requested,
      paper_method = paper_m,
      cell_ids = cell_ids,
      non_target_ids = val$non_target_ids,
      d = d_vec,
      sigma = sigma_vec,
      n_risk = n_risk_vec,
      n_switchers = n_switchers_vec,
      n_stayers = n_stayers_vec,
      pi_hat = pi_hat_vec,
      q0_hat = q0_hat_vec,
      Q_hat = Q_hat_vec,
      theta_hat = theta_hat,
      tau_hat = tau_hat,
      mu_1 = mu_1_vec,
      mu_0 = mu_0_vec,
      phi_theta = phi_theta,
      phi_tau = phi_tau,
      cell_components = cell_components,
      cov_theta = cov_theta_out,
      vhat_theta = vhat_theta_out,
      se_theta = se_theta_out,
      cov_tau = cov_tau_out,
      vhat_tau = vhat_tau_out,
      se_tau = se_tau_out,
      clustered = clustered_out,
      n_clusters = n_clusters_out,
      cluster = cluster,
      G_matrix = G_mat,
      partition = partition,
      glm_control = ctrl,
      rank_tol = rtol,
      nuisance = list(
        logit = logit_fits,
        or = or_fits
      ),
      aggregate = NULL
    ),
    class = "cdid_static_complete_fit"
  )

  if (!is.null(aggregate_cells) || !is.null(aggregate_signs)) {
    if (is.null(aggregate_cells) || is.null(aggregate_signs)) {
      static_complete_abort(
        "Supply both `aggregate_cells` and `aggregate_signs` when requesting aggregation inside `static_complete_fit()`.",
        c("cdid_static_aggregate_error", "cdid_static_input_error")
      )
    }
    out$aggregate <- static_complete_aggregate(
      fit = out,
      cells = aggregate_cells,
      signs = aggregate_signs,
      cluster = cluster,
      materialize_if = materialize_if
    )
  }

  out
}

# Aggregate a subset of fitted cells A <= C with explicit signs s in {-1, +1}^A
# and random switch-share weights omega_hat_{A,c} = n_c / sum_{c' in A} n_{c'}.
static_complete_aggregate <- function(fit,
                                      cells,
                                      signs,
                                      cluster = fit$cluster,
                                      materialize_if = !is.null(fit$phi_theta)) {
  if (!inherits(fit, "cdid_static_complete_fit") || !is.list(fit)) {
    static_complete_abort(
      "Expected a `cdid_static_complete_fit` object from `static_complete_fit()`.",
      c("cdid_static_aggregate_error", "cdid_static_input_error")
    )
  }
  if (!static_complete_scalar_logical(materialize_if)) {
    static_complete_abort(
      "`materialize_if` must be a non-missing logical scalar.",
      c("cdid_static_aggregate_error", "cdid_static_input_error")
    )
  }
  if (isTRUE(materialize_if) && is.null(fit$phi_theta)) {
    static_complete_abort(
      "Cannot materialize aggregate influence function from a point-only `cdid_static_complete_fit` object (`phi_theta` is NULL).",
      c("cdid_static_inference_error", "cdid_static_aggregate_error")
    )
  }
  if (missing(cells) || is.null(cells)) {
    static_complete_abort(
      "Supply explicit requested `cells` (`A <= C`) for aggregation; cell selection is never inferred.",
      c("cdid_static_aggregate_error", "cdid_static_input_error")
    )
  }
  if (is.factor(cells)) cells <- as.character(cells)
  if (!is.atomic(cells) || !is.null(dim(cells)) || is.object(cells) ||
      length(cells) < 1L || anyNA(cells)) {
    static_complete_abort(
      "`cells` must be a non-empty vector of declared cell IDs or indices.",
      c("cdid_static_aggregate_error", "cdid_static_input_error")
    )
  }
  all_cids <- fit$cell_ids
  if (is.numeric(cells)) {
    if (any(!is.finite(cells) | cells != floor(cells) | cells < 1 | cells > length(all_cids))) {
      static_complete_abort(
        "`cells` numeric indices must be valid integers in `1..length(fit$cell_ids)`.",
        c("cdid_static_aggregate_error", "cdid_static_input_error")
      )
    }
    sel_cids <- all_cids[as.integer(cells)]
  } else if (is.character(cells)) {
    if (any(!nzchar(cells)) || !all(cells %in% all_cids)) {
      static_complete_abort(
        "Requested `cells` contains unknown cell identifiers absent from `fit$cell_ids`.",
        c("cdid_static_aggregate_error", "cdid_static_input_error")
      )
    }
    sel_cids <- cells
  } else {
    static_complete_abort(
      "`cells` must be character cell IDs or integer indices.",
      c("cdid_static_aggregate_error", "cdid_static_input_error")
    )
  }
  if (anyDuplicated(sel_cids)) {
    static_complete_abort(
      "Requested `cells` for aggregation must not contain duplicate cells.",
      c("cdid_static_aggregate_error", "cdid_static_input_error")
    )
  }

  n_sel <- length(sel_cids)
  if (missing(signs) || is.null(signs) || !is.atomic(signs) || !is.null(dim(signs)) ||
      is.object(signs) || !is.numeric(signs) || length(signs) != n_sel ||
      anyNA(signs) || any(!is.finite(signs) | !(signs %in% c(-1, 1)))) {
    static_complete_abort(
      sprintf("Supply explicit `signs` in `{-1, 1}` of length `length(cells) = %d`.", n_sel),
      c("cdid_static_aggregate_error", "cdid_static_input_error")
    )
  }
  s_vec <- as.double(signs)
  s_names <- names(signs)
  if (!is.null(s_names)) {
    if (anyNA(s_names) || any(!nzchar(s_names)) || anyDuplicated(s_names) || !setequal(s_names, sel_cids)) {
      static_complete_abort(
        "Named `signs` must match the requested `cells`.",
        c("cdid_static_aggregate_error", "cdid_static_input_error")
      )
    }
    s_vec <- unname(s_vec[match(sel_cids, s_names)])
  }
  names(s_vec) <- sel_cids

  n <- fit$n
  n_num <- as.double(n)
  n_sw_A <- fit$n_switchers[sel_cids]
  total_sw_A <- sum(n_sw_A)
  if (!is.finite(total_sw_A) || total_sw_A < 1) {
    static_complete_abort(
      "Requested cell set `A` has zero total switchers; aggregate estimand is unavailable.",
      c("cdid_static_support_error", "cdid_ipw_support_error")
    )
  }

  pi_A <- fit$pi_hat[sel_cids]
  Pi_A <- sum(pi_A)
  omega_A <- n_sw_A / as.double(total_sw_A)
  names(omega_A) <- sel_cids

  theta_A <- fit$theta_hat[sel_cids]
  signed_theta_A <- s_vec * theta_A
  delta_hat <- sum(omega_A * signed_theta_A)

  if (isTRUE(materialize_if)) {
    # Random switch-share aggregate influence function:
    #   phi^m_{A,s} = sum_{c in A} omega_{A,c} s_c phi^m_c
    #               + (1 / Pi_A) sum_{c in A} (s_c theta^m_c - delta^m_{A,s}) (G_c - pi_c)
    phi_A <- fit$phi_theta[, sel_cids, drop = FALSE]
    G_A <- fit$G_matrix[, sel_cids, drop = FALSE]

    cell_if_term <- drop(phi_A %*% (omega_A * s_vec))
    G_centered_A <- sweep(G_A, 2L, pi_A, FUN = "-")
    share_weights_A <- (signed_theta_A - delta_hat) / Pi_A
    share_if_term <- drop(G_centered_A %*% share_weights_A)
    phi_agg <- cell_if_term + share_if_term

    if (!is.finite(delta_hat) || any(!is.finite(phi_agg))) {
      static_complete_abort(
        "Nonfinite aggregate estimate or influence function.",
        c("cdid_static_numeric_error", "cdid_ipw_numeric_error")
      )
    }

    inf_agg <- static_complete_covariance(phi_agg, cluster = cluster)
    joint_if <- cbind(phi_A, aggregate = phi_agg)
    inf_joint <- static_complete_covariance(joint_if, cluster = cluster)
    phi_comp_out <- list(
      cell_if_term = cell_if_term,
      share_if_term = share_if_term
    )
    var_out <- inf_agg$cov
    vhat_out <- inf_agg$vhat
    se_out <- inf_agg$se
    jcov_out <- inf_joint$cov_matrix
    jvhat_out <- inf_joint$vhat_matrix
    cl_out <- inf_agg$clustered
    ncl_out <- inf_agg$n_clusters
  } else {
    if (!is.finite(delta_hat)) {
      static_complete_abort(
        "Nonfinite aggregate estimate or influence function.",
        c("cdid_static_numeric_error", "cdid_ipw_numeric_error")
      )
    }
    c_info <- static_complete_validate_cluster(cluster, n = n)
    phi_agg <- NULL
    phi_comp_out <- NULL
    var_out <- NA_real_
    vhat_out <- NA_real_
    se_out <- NA_real_
    joint_if <- NULL
    jcov_out <- NULL
    jvhat_out <- NULL
    cl_out <- !is.null(c_info)
    ncl_out <- if (!is.null(c_info)) c_info$n_clusters else n
  }

  structure(
    list(
      delta_hat = delta_hat,
      method = fit$method,
      paper_method = fit$paper_method,
      cells = sel_cids,
      signs = s_vec,
      n = n,
      n_switchers = n_sw_A,
      total_switchers = total_sw_A,
      pi_hat = pi_A,
      Pi_hat = Pi_A,
      omega_hat = omega_A,
      theta_hat = theta_A,
      signed_theta_hat = signed_theta_A,
      phi = phi_agg,
      phi_components = phi_comp_out,
      variance = var_out,
      vhat = vhat_out,
      se = se_out,
      joint_if = joint_if,
      joint_cov = jcov_out,
      joint_vhat = jvhat_out,
      clustered = cl_out,
      n_clusters = ncl_out,
      cluster = cluster
    ),
    class = "cdid_static_complete_aggregate"
  )
}

# Re-evaluate original-unit or cluster covariance on a fitted static complete
# object without altering point estimates or unit-level influence functions.
static_complete_inference <- function(object, cluster = NULL) {
  if (inherits(object, "cdid_static_complete_fit")) {
    inf_theta <- static_complete_covariance(object$phi_theta, cluster = cluster)
    inf_tau <- static_complete_covariance(object$phi_tau, cluster = cluster)
    agg_inf <- if (!is.null(object$aggregate)) {
      static_complete_inference(object$aggregate, cluster = cluster)
    } else {
      NULL
    }
    return(list(
      theta_hat = object$theta_hat,
      tau_hat = object$tau_hat,
      cov_theta = inf_theta$cov_matrix,
      vhat_theta = inf_theta$vhat_matrix,
      se_theta = inf_theta$se,
      cov_tau = inf_tau$cov_matrix,
      vhat_tau = inf_tau$vhat_matrix,
      se_tau = inf_tau$se,
      clustered = inf_theta$clustered,
      n_clusters = inf_theta$n_clusters,
      aggregate = agg_inf
    ))
  }
  if (inherits(object, "cdid_static_complete_aggregate")) {
    inf_agg <- static_complete_covariance(object$phi, cluster = cluster)
    inf_joint <- static_complete_covariance(object$joint_if, cluster = cluster)
    return(list(
      delta_hat = object$delta_hat,
      variance = inf_agg$cov,
      vhat = inf_agg$vhat,
      se = inf_agg$se,
      joint_cov = inf_joint$cov_matrix,
      joint_vhat = inf_joint$vhat_matrix,
      clustered = inf_agg$clustered,
      n_clusters = inf_agg$n_clusters
    ))
  }
  static_complete_abort(
    "Expected a `cdid_static_complete_fit` or `cdid_static_complete_aggregate` object.",
    c("cdid_static_inference_error", "cdid_static_input_error")
  )
}
