# Private deterministic stayer graph reconstruction of source-adjacent outcome
# regression (OR), its complete original-unit influence function, and fixed
# switch-step link centering for the static binary missing-outcome extension.
#
# Mathematical components implemented here (Bellego, Benatia and Dortet-Bernadet
# 2025; de Chaisemartin, D'Haultfoeuille, Pasquier, Sow and Vazquez-Bare 2026,
# arXiv:2201.06898v7):
#   - Componentwise-grounded incidence D^a_c (rank |T^a_c| - K^a_c),
#     target ell^a_c = e_t (zero outside endpoint component T^a_{c,1}),
#     and affine flow space {varpi : (D^a_c)' varpi = ell^a_c}
#   - Deterministic LS flow varpi^or_c = D_c (D_c' D_c)^(-1) ell_c
#     (identically zero outside endpoint component T^C_{c,1}),
#     m_hat^C_c = sum_k varpi_{c,k} psi_hat^{C,0}_{c,k},
#     Cov_hat_c = sum_k varpi_{c,k} Cov_hat_{c,k},
#     pi_hat^C_c = E_n[C_c], zbar^C_c = E_n[C_c Z_c] / pi_hat^C_c,
#     S_hat_c = E_n[C_c (Z_c - zbar^C_c)(Z_c - zbar^C_c)'] / pi_hat^C_c
#   - Reconstructed stayer OR slopes and intercepts:
#     b_hat_P = (sum_{c in P} pi_hat^C_c S_hat_c)^(-1) sum_{c in P} pi_hat^C_c Cov_hat_c,
#     alpha_hat_c = m_hat^C_c - (zbar^C_c)' b_hat_{P(c)}
#   - Complete original-unit IFs chi^b_P and chi^alpha_c, retaining random stayer
#     shares pi_hat^C_c, means zbar^C_c, and known-X stayer Gram matrices S_hat_c
#   - Fixed switch-step link centering at [1{u >= t} - 1{s >= t}] using full-role
#     design means zG = E_G[Z_c] and zCr = E_{C,r}[Z_c], invariant to date-level
#     centering at non-switch dates
#
# Depends only on `cdid_abort` (R/forward.R) and compact QR helpers
# `or_wls_*` / `or_pooled_wls_*` (R/chained_or_wls.R, R/chained_or_pooled_wls.R).

static_reconstructed_or_abort <- function(message, class = character(), extra = list()) {
  cond <- c(list(message = message, call = NULL), extra)
  stop(structure(
    cond,
    class = unique(c(
      class,
      "cdid_static_reconstructed_or_error",
      "cdid_error",
      "error",
      "condition"
    ))
  ))
}

static_reconstructed_or_scalar_num <- function(x, finite = TRUE) {
  is.atomic(x) && is.null(dim(x)) && !is.object(x) &&
    is.numeric(x) && length(x) == 1L && !is.na(x) &&
    (!finite || is.finite(x))
}

static_reconstructed_or_validate_rank_tol <- function(rank_tol) {
  if (missing(rank_tol) || !static_reconstructed_or_scalar_num(rank_tol, finite = TRUE) ||
      rank_tol <= 0 || rank_tol >= 1) {
    static_reconstructed_or_abort(
      "Supply an explicit finite `rank_tol` in `(0, 1)`; no default rank tolerance is adopted.",
      c("cdid_static_reconstructed_or_rank_tol_error", "cdid_or_wls_input_error")
    )
  }
  as.double(rank_tol)
}

static_reconstructed_or_validate_chunk_size <- function(chunk_size, n_default = NULL) {
  if (is.null(chunk_size)) {
    if (is.null(n_default)) {
      static_reconstructed_or_abort(
        "Supply an explicit positive integer `chunk_size`.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
    return(as.integer(max(1L, n_default)))
  }
  if (!static_reconstructed_or_scalar_num(chunk_size, finite = TRUE) ||
      chunk_size < 1 || chunk_size != floor(chunk_size)) {
    static_reconstructed_or_abort(
      "`chunk_size` must be a finite positive integer scalar.",
      "cdid_static_reconstructed_or_input_error"
    )
  }
  as.integer(chunk_size)
}

static_reconstructed_or_validate_unit_id <- function(unit_id, what = "unit_id") {
  if (missing(unit_id) || is.null(unit_id)) {
    static_reconstructed_or_abort(
      sprintf("Explicit `%s` vector is required to verify original-unit identity and order.", what),
      "cdid_static_reconstructed_or_unit_id_error"
    )
  }
  if (is.factor(unit_id)) unit_id <- as.character(unit_id)
  if (!is.atomic(unit_id) || !is.null(dim(unit_id)) || is.object(unit_id) ||
      length(unit_id) < 1L || anyNA(unit_id)) {
    static_reconstructed_or_abort(
      sprintf("`%s` must be a non-empty atomic vector of original-unit identifiers without NA.", what),
      "cdid_static_reconstructed_or_unit_id_error"
    )
  }
  is_num <- is.numeric(unit_id)
  is_chr <- is.character(unit_id)
  if (!(is_num || is_chr) ||
      (is_num && any(!is.finite(unit_id))) ||
      (is_chr && any(!nzchar(unit_id)))) {
    static_reconstructed_or_abort(
      sprintf("`%s` must contain finite numeric or non-empty character identifiers.", what),
      "cdid_static_reconstructed_or_unit_id_error"
    )
  }
  key <- if (is_num) as.double(unit_id) else unit_id
  ids_chr <- as.character(unit_id)
  if (anyDuplicated(key) || anyDuplicated(ids_chr)) {
    static_reconstructed_or_abort(
      sprintf("`%s` must not contain duplicate original-unit identifiers.", what),
      "cdid_static_reconstructed_or_unit_id_error"
    )
  }
  list(id = unit_id, key = key, ids_chr = ids_chr, is_num = is_num, n = length(unit_id))
}

static_reconstructed_or_same_unit_id <- function(u1, u2) {
  identical(u1$is_num, u2$is_num) &&
    identical(u1$n, u2$n) &&
    identical(u1$key, u2$key)
}

static_reconstructed_or_validate_rows <- function(rows, n) {
  if (is.null(rows)) return(seq_len(n))
  if (!is.atomic(rows) || !is.null(dim(rows)) || is.object(rows) || anyNA(rows)) {
    static_reconstructed_or_abort(
      "`rows` must be NULL, a logical mask of length `n`, or integer row indices in `1..n`.",
      "cdid_static_reconstructed_or_input_error"
    )
  }
  if (is.logical(rows)) {
    if (length(rows) != n) {
      static_reconstructed_or_abort(
        sprintf("Logical `rows` mask must have length `n = %d`.", n),
        "cdid_static_reconstructed_or_input_error"
      )
    }
    return(which(rows))
  }
  if (!is.numeric(rows) || any(!is.finite(rows) | rows != floor(rows) | rows < 1 | rows > n)) {
    static_reconstructed_or_abort(
      sprintf("Integer `rows` indices must lie in `1..%d`.", n),
      "cdid_static_reconstructed_or_input_error"
    )
  }
  as.integer(rows)
}

static_reconstructed_or_validate_mask <- function(mask, n, ids_chr, what, cell_id = NULL) {
  prefix <- if (is.null(cell_id)) "" else sprintf("Cell `%s` ", cell_id)
  if (missing(mask) || is.null(mask) || !is.atomic(mask) || !is.null(dim(mask)) ||
      is.object(mask) || !(is.numeric(mask) || is.logical(mask)) ||
      length(mask) != n || anyNA(mask)) {
    static_reconstructed_or_abort(
      sprintf("%smask `%s` must be a nonmissing 0/1 or logical vector of length `n = %d`.", prefix, what, n),
      "cdid_static_reconstructed_or_mask_error"
    )
  }
  nm <- names(mask)
  if (!is.null(nm) && !identical(as.character(nm), ids_chr)) {
    static_reconstructed_or_abort(
      sprintf("%smask `%s` names do not match `unit_id` in the same order.", prefix, what),
      "cdid_static_reconstructed_or_order_error"
    )
  }
  if (is.logical(mask)) return(as.integer(mask))
  if (any(!is.finite(mask) | !(mask %in% c(0, 1)))) {
    static_reconstructed_or_abort(
      sprintf("%smask `%s` must contain only binary values in `{0, 1}`.", prefix, what),
      "cdid_static_reconstructed_or_mask_error"
    )
  }
  as.integer(mask)
}

static_reconstructed_or_validate_Z <- function(Z = NULL, X_tilde = NULL, n, ids_chr,
                                               active_rows, cell_id) {
  has_Z <- !is.null(Z)
  has_X <- !is.null(X_tilde)
  if (!has_Z && !has_X) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` requires explicit slope design `Z` (without intercept, `ncol >= 0`) and/or `X_tilde`.", cell_id),
      "cdid_static_reconstructed_or_design_error"
    )
  }
  if (has_X) {
    if (!is.matrix(X_tilde) || !is.numeric(X_tilde) || is.object(X_tilde) ||
        nrow(X_tilde) != n || ncol(X_tilde) < 1L) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `X_tilde` must be a numeric matrix with `nrow == %d` and `ncol >= 1`.", cell_id, n),
        "cdid_static_reconstructed_or_design_error"
      )
    }
    if (!is.null(rownames(X_tilde)) && !identical(as.character(rownames(X_tilde)), ids_chr)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `X_tilde` rownames do not match `unit_id` in order.", cell_id),
        "cdid_static_reconstructed_or_order_error"
      )
    }
    X_act <- X_tilde[active_rows, , drop = FALSE]
    if (anyNA(X_act) || any(!is.finite(X_act))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `X_tilde` must be finite on all active role units.", cell_id),
        "cdid_static_reconstructed_or_design_error"
      )
    }
    if (!all(X_act[, 1L] == 1)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `X_tilde` must have leading intercept column equal to 1.", cell_id),
        "cdid_static_reconstructed_or_design_error"
      )
    }
    storage.mode(X_tilde) <- "double"
  }
  if (has_Z) {
    if (is.atomic(Z) && is.null(dim(Z)) && !is.object(Z) && is.numeric(Z) && length(Z) == n) {
      nm <- names(Z)
      if (!is.null(nm) && !identical(as.character(nm), ids_chr)) {
        static_reconstructed_or_abort(
          sprintf("Cell `%s` `Z` names do not match `unit_id` in order.", cell_id),
          "cdid_static_reconstructed_or_order_error"
        )
      }
      Z <- matrix(as.double(Z), nrow = n, ncol = 1L)
    }
    if (!is.matrix(Z) || !is.numeric(Z) || is.object(Z) || nrow(Z) != n) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `Z` must be a numeric matrix with `nrow == %d` (0 or more columns).", cell_id, n),
        "cdid_static_reconstructed_or_design_error"
      )
    }
    if (!is.null(rownames(Z)) && !identical(as.character(rownames(Z)), ids_chr)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `Z` rownames do not match `unit_id` in order.", cell_id),
        "cdid_static_reconstructed_or_order_error"
      )
    }
    Z_act <- Z[active_rows, , drop = FALSE]
    if (anyNA(Z_act) || any(!is.finite(Z_act))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `Z` must be finite on all active stayer units.", cell_id),
        "cdid_static_reconstructed_or_design_error"
      )
    }
    storage.mode(Z) <- "double"
  }
  if (has_X && has_Z) {
    if (ncol(X_tilde) != ncol(Z) + 1L ||
        !identical(unname(X_tilde[active_rows, -1L, drop = FALSE]),
                   unname(Z[active_rows, , drop = FALSE]))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `X_tilde` and `Z` disagree on active stayer units.", cell_id),
        "cdid_static_reconstructed_or_design_error"
      )
    }
  } else if (has_X) {
    Z <- X_tilde[, -1L, drop = FALSE]
  }

  p <- ncol(Z)
  z_names <- colnames(Z)
  if (is.null(z_names)) {
    if (p > 0L) {
      x_names <- if (has_X) colnames(X_tilde) else NULL
      if (!is.null(x_names) && length(x_names) == p + 1L &&
          all(!is.na(x_names[-1L])) && all(nzchar(x_names[-1L]))) {
        z_names <- x_names[-1L]
      } else {
        z_names <- paste0("z", seq_len(p))
      }
    } else {
      z_names <- character(0L)
    }
  } else {
    if (anyNA(z_names) || any(!nzchar(z_names)) || anyDuplicated(z_names)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` slope matrix `Z` column names must be non-empty and unique.", cell_id),
        "cdid_static_reconstructed_or_design_error"
      )
    }
  }
  if ("(Intercept)" %in% z_names) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` slope matrix `Z` must not contain an `(Intercept)` column.", cell_id),
      "cdid_static_reconstructed_or_design_error"
    )
  }
  colnames(Z) <- z_names
  list(Z = Z, p = p, z_names = z_names)
}

# Parse and validate a two-column edge specification (s_k, u_k) with s_k < u_k.
static_reconstructed_or_parse_edges <- function(edges, cell_id = "<cell>") {
  if (is.data.frame(edges)) {
    if (ncol(edges) != 2L) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `edges` data.frame must have exactly 2 columns `(s, u)`.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    edges <- cbind(as.double(edges[[1L]]), as.double(edges[[2L]]))
  } else if (is.list(edges) && !is.object(edges)) {
    if (!length(edges)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `edges` must be non-empty.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
     pairs <- lapply(edges, function(e) {
      if (!is.atomic(e) || !is.numeric(e) || length(e) != 2L || any(!is.finite(e))) {
        static_reconstructed_or_abort(
          sprintf("Cell `%s` `edges` list elements must be finite numeric pairs `c(s, u)`.", cell_id),
          "cdid_static_reconstructed_or_graph_error"
        )
      }
      as.double(e)
    })
    edges <- do.call(rbind, pairs)
  } else if (is.matrix(edges) && is.numeric(edges) && !is.object(edges)) {
    if (ncol(edges) != 2L || nrow(edges) < 1L || any(!is.finite(edges))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `edges` must be a finite numeric matrix with 2 columns and >= 1 row.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    storage.mode(edges) <- "double"
  } else {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` `edges` must be a 2-column matrix, data.frame, or list of `(s, u)` pairs.", cell_id),
      "cdid_static_reconstructed_or_graph_error"
    )
  }
  s <- edges[, 1L]
  u <- edges[, 2L]
  if (any(s >= u)) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` `edges` must satisfy `s < u` on every edge.", cell_id),
      "cdid_static_reconstructed_or_graph_error"
    )
  }
  edge_keys <- paste0(s, "->", u)
  if (anyDuplicated(edge_keys)) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` `edges` contains duplicate `(s, u)` pairs.", cell_id),
      "cdid_static_reconstructed_or_graph_error"
    )
  }
  colnames(edges) <- c("s", "u")
  rownames(edges) <- edge_keys
  list(edges = edges, s = s, u = u, edge_ids = edge_keys, n_edges = nrow(edges))
}

# Build and/or validate componentwise-grounded incidence matrix D_c, target
# vector ell_c = e_t, endpoint connectivity between baseline (t - 1) and target
# (t), and one deterministic coordinate gauge per connected component.
static_reconstructed_or_graph <- function(D = NULL,
                                          ell = NULL,
                                          edges = NULL,
                                          t = NULL,
                                          baseline_date = NULL,
                                          dates = NULL,
                                          reference_dates = NULL,
                                          cell_id = "<cell>") {
  has_edges <- !is.null(edges)
  has_D <- !is.null(D) || !is.null(ell)
  if (!has_edges && !has_D) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` requires explicit grounded incidence `(D, ell)` and/or `(edges, t)`.", cell_id),
      "cdid_static_reconstructed_or_graph_error"
    )
  }

  t_val <- NULL
  base_val <- NULL
  if (!is.null(t)) {
    if (!static_reconstructed_or_scalar_num(t, finite = TRUE)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` transition date `t` must be a finite numeric scalar.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    t_val <- as.double(t)
  }
  if (!is.null(baseline_date)) {
    if (!static_reconstructed_or_scalar_num(baseline_date, finite = TRUE)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `baseline_date` must be a finite numeric scalar.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    base_val <- as.double(baseline_date)
  } else if (!is.null(t_val)) {
    base_val <- t_val - 1.0
  }
  if (!is.null(t_val) && !is.null(base_val) && !(base_val < t_val)) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` requires `baseline_date < t`.", cell_id),
      "cdid_static_reconstructed_or_graph_error"
    )
  }

  edge_info <- NULL
  vertices <- NULL
  n_comp <- NA_integer_
  comp_id <- NULL
  edge_comp_id <- NULL
  col_comp_id <- NULL
  ref_dates <- NULL
  non_ground <- NULL
  iso_dates <- numeric(0L)

  if (has_edges) {
    if (is.null(t_val)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` with `edges` requires transition date `t`.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    edge_info <- static_reconstructed_or_parse_edges(edges, cell_id = cell_id)
    edge_vertices <- sort(unique(c(edge_info$s, edge_info$u)))
    if (!is.null(dates)) {
      if (!is.atomic(dates) || !is.null(dim(dates)) || is.object(dates) ||
          !is.numeric(dates) || length(dates) < 2L || anyNA(dates) || any(!is.finite(dates)) ||
          anyDuplicated(as.double(dates))) {
        static_reconstructed_or_abort(
          sprintf("Cell `%s` `dates` must be a finite numeric vector of unique dates.", cell_id),
          "cdid_static_reconstructed_or_graph_error"
        )
      }
      dates_num <- sort(as.double(dates))
      if (!all(edge_vertices %in% dates_num)) {
        static_reconstructed_or_abort(
          sprintf("Cell `%s` `edges` contain endpoints not present in declared `dates`.", cell_id),
          "cdid_static_reconstructed_or_graph_error"
        )
      }
      vertices <- dates_num
    } else {
      vertices <- edge_vertices
    }

    if (!(base_val %in% vertices) || !(t_val %in% vertices) ||
        !(base_val %in% edge_vertices) || !(t_val %in% edge_vertices)) {
      static_reconstructed_or_abort(
        sprintf(
          "Cell `%s` role graph does not contain non-isolated baseline date `%g` and target date `%g` (Assumption as:graph).",
          cell_id, base_val, t_val
        ),
        "cdid_static_reconstructed_or_graph_error"
      )
    }

    # Partition vertices into connected components T_{c,1}, ..., T_{c,K_c},
    # placing the endpoint component containing baseline_date = t - 1 first (m = 1).
    n_v <- length(vertices)
    s_idx <- match(edge_info$s, vertices)
    u_idx <- match(edge_info$u, vertices)
    comp_id <- integer(n_v)
    deg <- integer(n_v)
    for (k in seq_len(edge_info$n_edges)) {
      deg[s_idx[k]] <- deg[s_idx[k]] + 1L
      deg[u_idx[k]] <- deg[u_idx[k]] + 1L
    }
    iso_dates <- vertices[deg == 0L]

    bfs_component <- function(start_v, cid_val) {
      comp_id[start_v] <<- cid_val
      queue <- start_v
      while (length(queue) > 0L) {
        curr <- queue[1L]
        queue <- queue[-1L]
        inc <- which(s_idx == curr | u_idx == curr)
        if (length(inc) > 0L) {
          nbrs <- unique(c(s_idx[inc], u_idx[inc]))
          new_nbrs <- nbrs[comp_id[nbrs] == 0L]
          if (length(new_nbrs) > 0L) {
            comp_id[new_nbrs] <<- cid_val
            queue <- c(queue, new_nbrs)
          }
        }
      }
    }

    root_idx <- match(base_val, vertices)
    targ_idx <- match(t_val, vertices)
    bfs_component(root_idx, 1L)
    if (comp_id[targ_idx] != 1L) {
      static_reconstructed_or_abort(
        sprintf(
          "Cell `%s` role graph is disconnected between baseline date `%g` and target date `%g` (Assumption as:graph).",
          cell_id, base_val, t_val
        ),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    n_comp <- 1L
    for (v in seq_len(n_v)) {
      if (comp_id[v] == 0L) {
        n_comp <- n_comp + 1L
        bfs_component(v, n_comp)
      }
    }

    # Select one deterministic reference date r_{c,m} per connected component,
    # with r_{c,1} = t - 1 (base_val) in the required endpoint component.
    ref_arg <- NULL
    if (!is.null(reference_dates)) {
      if (!is.atomic(reference_dates) || !is.null(dim(reference_dates)) || is.object(reference_dates) ||
          !is.numeric(reference_dates) || length(reference_dates) < 1L ||
          anyNA(reference_dates) || any(!is.finite(reference_dates)) ||
          anyDuplicated(as.double(reference_dates)) || !all(as.double(reference_dates) %in% vertices)) {
        static_reconstructed_or_abort(
          sprintf("Cell `%s` `reference_dates` must be unique finite dates belonging to the cell graph.", cell_id),
          "cdid_static_reconstructed_or_graph_error"
        )
      }
      ref_arg <- as.double(reference_dates)
    }

    ref_dates <- numeric(n_comp)
    ref_dates[1L] <- base_val
    if (!is.null(ref_arg)) {
      c1_v <- vertices[comp_id == 1L]
      c1_cand <- c1_v[c1_v %in% ref_arg]
      if (length(c1_cand) > 0L && !identical(c1_cand, base_val)) {
        static_reconstructed_or_abort(
          sprintf("Cell `%s` endpoint component reference date must be `baseline_date = %g`.", cell_id, base_val),
          "cdid_static_reconstructed_or_graph_error"
        )
      }
    }
    if (n_comp >= 2L) {
      for (m_c in 2L:n_comp) {
        cm_v <- vertices[comp_id == m_c]
        if (length(cm_v) == 1L) {
          ref_dates[m_c] <- cm_v[1L]
        } else if (is.null(ref_arg)) {
          ref_dates[m_c] <- cm_v[1L]
        } else {
          cm_cand <- cm_v[cm_v %in% ref_arg]
          if (length(cm_cand) != 1L) {
            static_reconstructed_or_abort(
              sprintf("Cell `%s` `reference_dates` must specify exactly one reference date in component %d.", cell_id, m_c),
              "cdid_static_reconstructed_or_graph_error"
            )
          }
          ref_dates[m_c] <- cm_cand[1L]
        }
      }
    }

    non_ground <- vertices[!(vertices %in% ref_dates)]
    K <- edge_info$n_edges
    V <- length(non_ground)
    col_comp_id <- comp_id[match(non_ground, vertices)]
    edge_comp_id <- comp_id[s_idx]

    D_from_edges <- matrix(0.0, nrow = K, ncol = V,
                           dimnames = list(edge_info$edge_ids, as.character(non_ground)))
    for (k in seq_len(K)) {
      if (!(edge_info$s[k] %in% ref_dates)) {
        D_from_edges[k, match(edge_info$s[k], non_ground)] <- -1.0
      }
      if (!(edge_info$u[k] %in% ref_dates)) {
        D_from_edges[k, match(edge_info$u[k], non_ground)] <- 1.0
      }
    }
    ell_from_edges <- numeric(V)
    names(ell_from_edges) <- as.character(non_ground)
    ell_from_edges[match(t_val, non_ground)] <- 1.0

    if (!has_D) {
      D <- D_from_edges
      ell <- ell_from_edges
    } else {
      if (is.null(D) || is.null(ell) || !is.matrix(D) || !is.numeric(D) ||
          !identical(dim(D), dim(D_from_edges)) || length(ell) != V ||
          max(abs(unname(D) - unname(D_from_edges))) > 0 ||
          max(abs(unname(as.double(ell)) - unname(ell_from_edges))) > 0) {
        static_reconstructed_or_abort(
          sprintf("Cell `%s` supplied `(D, ell)` does not match grounded incidence implied by `(edges, t)`.", cell_id),
          "cdid_static_reconstructed_or_graph_error"
        )
      }
    }
  } else {
    if (is.null(D) || is.null(ell)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` must supply both `D` and `ell` when `edges` is omitted.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    if (!is.matrix(D) || !is.numeric(D) || is.object(D) ||
        nrow(D) < 1L || ncol(D) < 1L || anyNA(D) || any(!is.finite(D))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `D` must be a finite numeric matrix with `nrow >= ncol >= 1`.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    if (!is.atomic(ell) || !is.null(dim(ell)) || is.object(ell) ||
        !is.numeric(ell) || length(ell) != ncol(D) || anyNA(ell) || any(!is.finite(ell))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `ell` must be a finite numeric vector of length `ncol(D) = %d`.", cell_id, ncol(D)),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    storage.mode(D) <- "double"
    ell <- as.double(ell)
    K <- nrow(D)
    V <- ncol(D)
    if (K < V) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` incidence `D` has fewer edges (%d) than non-ground dates (%d); graph is disconnected.", cell_id, K, V),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    # Validate grounded incidence structure: each row has entries in {-1, 0, 1},
    # at most one +1, at most one -1, and at least one nonzero entry.
    if (any(!(D %in% c(-1, 0, 1)))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` grounded incidence `D` must contain only entries in `{-1, 0, 1}`.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    pos_counts <- rowSums(D == 1)
    neg_counts <- rowSums(D == -1)
    if (any(pos_counts > 1L | neg_counts > 1L | (pos_counts + neg_counts) == 0L)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` grounded incidence `D` contains invalid edge rows.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    # Ground reference vertices are represented by rows where pos_counts + neg_counts == 1.
    if (!any((pos_counts + neg_counts) == 1L)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` grounded incidence `D` has no edge incident to baseline ground `t - 1`.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    # Target vector ell must be e_t (exactly one 1 and V - 1 zeros).
    if (!(sum(ell == 1) == 1L && sum(ell == 0) == (V - 1L))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` target `ell` must be a standard basis vector `e_t` in `{0, 1}^ncol(D)`.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    # Verify that every connected component of columns 1..V is connected to a
    # component reference vertex (represented by vertex 1L on rows with |nz| == 1).
    visited <- logical(V + 1L)
    visited[1L] <- TRUE
    queue <- 1L
    while (length(queue) > 0L) {
      curr <- queue[1L]
      queue <- queue[-1L]
      for (k in seq_len(K)) {
        nz <- which(D[k, ] != 0) + 1L
        endpts <- if (length(nz) == 1L) c(1L, nz) else nz
        if (curr %in% endpts) {
          other <- endpts[endpts != curr]
          new_v <- other[!visited[other]]
          if (length(new_v) > 0L) {
            visited[new_v] <- TRUE
            queue <- c(queue, new_v)
          }
        }
      }
    }
    if (!all(visited)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` grounded incidence `D` describes a disconnected role graph.", cell_id),
        "cdid_static_reconstructed_or_graph_error"
      )
    }
    # Identify independent column/row factorization blocks of D, target block first.
    # Grounding can split one original component into several such blocks, so
    # these blocks do not identify the original graph topology from D alone.
    col_comp_id <- integer(V)
    bfs_col_comp <- function(start_c, cid_val) {
      col_comp_id[start_c] <<- cid_val
      q_c <- start_c
      while (length(q_c) > 0L) {
        curr_c <- q_c[1L]
        q_c <- q_c[-1L]
        inc_rows <- which(D[, curr_c] != 0)
        if (length(inc_rows) > 0L) {
          nbr_cols <- which(colSums(D[inc_rows, , drop = FALSE] != 0) > 0L)
          new_cols <- nbr_cols[col_comp_id[nbr_cols] == 0L]
          if (length(new_cols) > 0L) {
            col_comp_id[new_cols] <<- cid_val
            q_c <- c(q_c, new_cols)
          }
        }
      }
    }
    targ_col <- which(ell == 1)[1L]
    bfs_col_comp(targ_col, 1L)
    n_comp <- 1L
    for (cc in seq_len(V)) {
      if (col_comp_id[cc] == 0L) {
        n_comp <- n_comp + 1L
        bfs_col_comp(cc, n_comp)
      }
    }
    edge_comp_id <- integer(K)
    for (k in seq_len(K)) {
      edge_comp_id[k] <- col_comp_id[which(D[k, ] != 0)[1L]]
    }
  }

  list(
    D = D,
    ell = ell,
    edges = if (!is.null(edge_info)) edge_info$edges else NULL,
    edge_ids = if (!is.null(edge_info)) edge_info$edge_ids else rownames(D),
    t = t_val,
    baseline_date = base_val,
    n_edges = nrow(D),
    n_non_ground_dates = ncol(D),
    dates = vertices,
    n_components = if (!is.null(edge_info)) as.integer(n_comp) else NA_integer_,
    component_id = if (!is.null(comp_id)) stats::setNames(as.integer(comp_id), as.character(vertices)) else NULL,
    factor_col_block_id = as.integer(col_comp_id),
    factor_edge_block_id = as.integer(edge_comp_id),
    edge_component_id = if (!is.null(edge_info)) {
      stats::setNames(as.integer(edge_comp_id), edge_info$edge_ids)
    } else {
      NULL
    },
    reference_dates = ref_dates,
    non_reference_dates = non_ground,
    isolated_dates = iso_dates
  )
}

# Deterministic least-squares flow varpi = D (D' D)^(-1) ell via column-pivoted
# QR and triangular solves (forwardsolve / backsolve), never matrix inversion.
# Validates D' varpi = ell. Outside the endpoint component T_{c,1}, ell = 0 and
# the minimum-norm flow varpi is identically zero.
static_reconstructed_or_ls_flow <- function(D = NULL,
                                            ell = NULL,
                                            edges = NULL,
                                            t = NULL,
                                            baseline_date = NULL,
                                            dates = NULL,
                                            reference_dates = NULL,
                                            rank_tol,
                                            flow = NULL,
                                            cell_id = "<cell>") {
  rtol <- static_reconstructed_or_validate_rank_tol(rank_tol)
  g <- static_reconstructed_or_graph(
    D = D, ell = ell, edges = edges, t = t,
    baseline_date = baseline_date, dates = dates,
    reference_dates = reference_dates, cell_id = cell_id
  )
  D_mat <- g$D
  ell_vec <- g$ell
  K <- nrow(D_mat)
  V <- ncol(D_mat)

  active_comps <- unique(g$factor_col_block_id)
  if (length(active_comps) <= 1L) {
    q_D <- qr(D_mat, LAPACK = TRUE)
    R_D <- qr.R(q_D)
    piv_D <- q_D$pivot
  } else {
    R_D <- matrix(0.0, nrow = V, ncol = V)
    piv_D <- integer(V)
    offset <- 0L
    for (m_comp in active_comps) {
      cols_m <- which(g$factor_col_block_id == m_comp)
      rows_m <- which(g$factor_edge_block_id == m_comp)
      k_m <- length(cols_m)
      if (length(rows_m) < k_m) {
        static_reconstructed_or_abort(
          sprintf(
            "Cell `%s` grounded incidence `D` component %d has fewer edges (%d) than non-reference dates (%d).",
            cell_id, m_comp, length(rows_m), k_m
          ),
          c("cdid_static_reconstructed_or_graph_error", "cdid_or_wls_rank_error")
        )
      }
      q_m <- qr(D_mat[rows_m, cols_m, drop = FALSE], LAPACK = TRUE)
      R_m <- qr.R(q_m)
      idx_m <- offset + seq_len(k_m)
      R_D[idx_m, idx_m] <- R_m
      piv_D[idx_m] <- cols_m[q_m$pivot]
      offset <- offset + k_m
    }
  }
  diag_R <- abs(diag(R_D))
  cutoff_R <- rtol * max(diag_R)
  rank_R <- if (cutoff_R > 0) sum(diag_R > cutoff_R) else 0L
  if (rank_R < V) {
    static_reconstructed_or_abort(
      sprintf(
        "Cell `%s` grounded incidence `D` has deficient rank %d < %d under `rank_tol = %g`.",
        cell_id, rank_R, V, rtol
      ),
      c("cdid_static_reconstructed_or_graph_error", "cdid_or_wls_rank_error")
    )
  }

  if (is.null(flow)) {
    # Solve (D' D) u = ell in pivoted order R_D' R_D u_piv = ell[pivot]
    ell_piv <- ell_vec[piv_D]
    u_piv <- backsolve(R_D, forwardsolve(t(R_D), ell_piv))
    u <- numeric(V)
    u[piv_D] <- u_piv
    if (!is.null(colnames(D_mat))) {
      names(u) <- colnames(D_mat)
    }
    varpi <- drop(D_mat %*% u)
    flow_rule <- "least_squares"
  } else {
    if (!is.atomic(flow) || !is.null(dim(flow)) || is.object(flow) ||
        !is.numeric(flow) || length(flow) != K || anyNA(flow) || any(!is.finite(flow))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` explicit `flow` must be a finite numeric vector of length `nrow(D) = %d`.", cell_id, K),
        "cdid_static_reconstructed_or_flow_error"
      )
    }
    varpi <- as.double(flow)
    u <- NULL
    flow_rule <- "explicit"
  }

  if (!is.null(g$edge_ids)) {
    names(varpi) <- g$edge_ids
  }
  lhs <- drop(crossprod(D_mat, varpi))
  residual <- lhs - ell_vec
  max_abs_resid <- max(abs(residual))
  flow_tol <- max(rtol, sqrt(.Machine$double.eps)) * (1.0 + max(abs(ell_vec)))
  if (!is.finite(max_abs_resid) || max_abs_resid > flow_tol) {
    static_reconstructed_or_abort(
      sprintf(
        "Cell `%s` flow violates `D' varpi = ell` (max|D' varpi - ell| = %g > %g).",
        cell_id, max_abs_resid, flow_tol
      ),
      "cdid_static_reconstructed_or_flow_error"
    )
  }

  list(
    varpi = varpi,
    flow_rule = flow_rule,
    D = D_mat,
    ell = ell_vec,
    edges = g$edges,
    edge_ids = g$edge_ids,
    t = g$t,
    baseline_date = g$baseline_date,
    dates = g$dates,
    n_components = g$n_components,
    component_id = g$component_id,
    edge_component_id = g$edge_component_id,
    reference_dates = g$reference_dates,
    non_reference_dates = g$non_reference_dates,
    isolated_dates = g$isolated_dates,
    qr_R = unname(R_D),
    qr_pivot = piv_D,
    rank = rank_R,
    rank_tol = rtol,
    cutoff = cutoff_R,
    potential_u = u,
    flow_residual = residual,
    max_abs_flow_residual = max_abs_resid
  )
}

# Normalize a single edge mean or edge covariance primitive into a uniform
# representation with:
#   - $estimate (scalar for mean, length-p vector for cov)
#   - $eval_if(rows) -> length(rows) vector (mean) or length(rows) x p matrix (cov)
#   - $loadings (optional keyed list of parameter loadings)
static_reconstructed_or_normalize_edge_moment <- function(obj, kind, n, ids_chr,
                                                          p = NULL, z_names = NULL,
                                                          coef_if = NULL,
                                                          cell_id, edge_idx) {
  what <- sprintf("Cell `%s` edge %d %s", cell_id, edge_idx, kind)
  if (is.null(obj) || !is.list(obj)) {
    static_reconstructed_or_abort(
      sprintf("%s must be a list (e.g. `static_obs_%s()` output or `list(estimate = ..., influence = ...)`).", what, kind),
      "cdid_static_reconstructed_or_edge_error"
    )
  }
  est <- obj[["estimate"]]
  if (kind == "mean") {
    if (!static_reconstructed_or_scalar_num(est, finite = TRUE)) {
      static_reconstructed_or_abort(
        sprintf("%s `$estimate` must be a finite numeric scalar.", what),
        "cdid_static_reconstructed_or_edge_error"
      )
    }
    est <- as.double(est)
  } else {
    if (is.null(est) || !is.atomic(est) || !is.null(dim(est)) || is.object(est) ||
        !is.numeric(est) || length(est) != p || anyNA(est) || any(!is.finite(est))) {
      static_reconstructed_or_abort(
        sprintf("%s `$estimate` must be a finite numeric vector of length `p = %d`.", what, p),
        "cdid_static_reconstructed_or_edge_error"
      )
    }
    est_names <- names(est)
    if (!is.null(est_names) && !is.null(z_names) && !identical(est_names, z_names)) {
      static_reconstructed_or_abort(
        sprintf("%s `$estimate` names do not match cell `Z` column names.", what),
        "cdid_static_reconstructed_or_edge_error"
      )
    }
    est <- stats::setNames(as.double(est), z_names)
  }

  loadings <- if (!is.null(obj[["loadings"]])) obj[["loadings"]] else list()
  if (!is.list(loadings) || is.object(loadings)) {
    static_reconstructed_or_abort(
      sprintf("%s `$loadings` must be a list when present.", what),
      "cdid_static_reconstructed_or_edge_error"
    )
  }

  check_if_slice <- function(val, rows) {
    m <- length(rows)
    row_ids <- ids_chr[rows]
    if (kind == "mean") {
      if (is.matrix(val) && ncol(val) == 1L) {
        rn <- rownames(val)
        val <- drop(val[, 1L])
        names(val) <- rn
      }
      if (!is.atomic(val) || !is.null(dim(val)) || is.object(val) ||
          !is.numeric(val) || length(val) != m || anyNA(val) || any(!is.finite(val))) {
        static_reconstructed_or_abort(
          sprintf("%s influence function must be a finite numeric vector of length %d.", what, m),
          "cdid_static_reconstructed_or_edge_error"
        )
      }
      if (!is.null(names(val)) && !identical(as.character(names(val)), row_ids)) {
        static_reconstructed_or_abort(
          sprintf("%s influence function names do not match `unit_id[rows]` in order.", what),
          "cdid_static_reconstructed_or_order_error"
        )
      }
      return(as.double(val))
    } else {
      if (p == 1L && is.atomic(val) && is.null(dim(val)) && is.numeric(val) && length(val) == m) {
        nm <- names(val)
        val <- matrix(as.double(val), nrow = m, ncol = 1L, dimnames = list(nm, z_names))
      }
      if (!is.matrix(val) || !is.numeric(val) || is.object(val) ||
          nrow(val) != m || ncol(val) != p || anyNA(val) || any(!is.finite(val))) {
        static_reconstructed_or_abort(
          sprintf("%s influence function must be a finite numeric `%d x %d` matrix.", what, m, p),
          "cdid_static_reconstructed_or_edge_error"
        )
      }
      if (!is.null(rownames(val)) && !identical(as.character(rownames(val)), row_ids)) {
        static_reconstructed_or_abort(
          sprintf("%s influence function rownames do not match `unit_id[rows]` in order.", what),
          "cdid_static_reconstructed_or_order_error"
        )
      }
      if (!is.null(colnames(val)) && !is.null(z_names) && !identical(colnames(val), z_names)) {
        static_reconstructed_or_abort(
          sprintf("%s influence function colnames do not match cell `Z` column names.", what),
          "cdid_static_reconstructed_or_edge_error"
        )
      }
      storage.mode(val) <- "double"
      colnames(val) <- z_names
      return(val)
    }
  }

  fn_cb <- if (!is.null(obj[["if_fn"]])) obj[["if_fn"]] else obj[["eval_if"]]
  inf_stored <- NULL
  if (!is.null(fn_cb) && is.function(fn_cb)) {
    eval_if <- function(rows = seq_len(n)) {
      rows <- static_reconstructed_or_validate_rows(rows, n)
      check_if_slice(fn_cb(rows), rows)
    }
    eval_direct <- if (!is.null(obj[["direct_if"]])) {
      dir_full <- check_if_slice(obj[["direct_if"]], seq_len(n))
      function(rows = seq_len(n)) {
        rows <- static_reconstructed_or_validate_rows(rows, n)
        if (kind == "mean") dir_full[rows] else dir_full[rows, , drop = FALSE]
      }
    } else {
      eval_if
    }
    if (!is.null(obj[["influence"]])) {
      inf_stored <- check_if_slice(obj[["influence"]], seq_len(n))
    }
  } else if (!is.null(obj[["influence"]])) {
    inf_full <- check_if_slice(obj[["influence"]], seq_len(n))
    inf_stored <- inf_full
    dir_full <- if (!is.null(obj[["direct_if"]])) check_if_slice(obj[["direct_if"]], seq_len(n)) else inf_full
    eval_if <- function(rows = seq_len(n)) {
      rows <- static_reconstructed_or_validate_rows(rows, n)
      if (kind == "mean") inf_full[rows] else inf_full[rows, , drop = FALSE]
    }
    eval_direct <- function(rows = seq_len(n)) {
      rows <- static_reconstructed_or_validate_rows(rows, n)
      if (kind == "mean") dir_full[rows] else dir_full[rows, , drop = FALSE]
    }
  } else if (!is.null(obj[["direct_if"]])) {
    dir_full <- check_if_slice(obj[["direct_if"]], seq_len(n))
    if (length(loadings) > 0L) {
      if (is.null(coef_if)) {
        static_reconstructed_or_abort(
          sprintf("%s has nuisance `$loadings` but neither `$influence` nor `coef_if` was supplied.", what),
          "cdid_static_reconstructed_or_edge_error"
        )
      }
      static_reconstructed_or_contract_loadings(
        loadings = loadings, coef_if = coef_if, rows = 1L, n = n, ids_chr = ids_chr
      )
    }
    eval_direct <- function(rows = seq_len(n)) {
      rows <- static_reconstructed_or_validate_rows(rows, n)
      if (kind == "mean") dir_full[rows] else dir_full[rows, , drop = FALSE]
    }
    eval_if <- function(rows = seq_len(n)) {
      rows <- static_reconstructed_or_validate_rows(rows, n)
      base_val <- if (kind == "mean") dir_full[rows] else dir_full[rows, , drop = FALSE]
      if (!length(loadings)) return(base_val)
      contrib <- static_reconstructed_or_contract_loadings(
        loadings = loadings, coef_if = coef_if, rows = rows, n = n, ids_chr = ids_chr
      )
      if (kind == "mean") base_val + drop(contrib) else base_val + contrib
    }
  } else {
    static_reconstructed_or_abort(
      sprintf("%s must provide `$influence`, `$if_fn`, or (`$direct_if` + `coef_if`).", what),
      "cdid_static_reconstructed_or_edge_error"
    )
  }

  list(
    kind = kind,
    estimate = est,
    loadings = loadings,
    influence = inf_stored,
    eval_if = eval_if,
    eval_direct = eval_direct
  )
}

# Contract keyed loadings with keyed coefficient IF matrices (or row callbacks)
# over requested original-unit `rows` without materializing N x total_params.
static_reconstructed_or_contract_loadings <- function(loadings, coef_if, rows, n, ids_chr) {
  rows <- static_reconstructed_or_validate_rows(rows, n)
  m <- length(rows)
  lk <- names(loadings)
  if (is.null(lk) || !length(lk)) {
    return(matrix(0.0, nrow = m, ncol = 1L))
  }
  if (is.null(coef_if) || !is.list(coef_if) || is.object(coef_if)) {
    static_reconstructed_or_abort(
      "`coef_if` must be a keyed list of `n x p_k` matrices or `function(rows)` callbacks.",
      "cdid_static_reconstructed_or_channel_error"
    )
  }
  ck <- names(coef_if)
  if (is.null(ck) || anyNA(ck) || any(!nzchar(ck)) || anyDuplicated(ck)) {
    static_reconstructed_or_abort(
      "`coef_if` keys must be non-empty and unique.",
      "cdid_static_reconstructed_or_channel_error"
    )
  }
  miss <- setdiff(lk, ck)
  if (length(miss)) {
    static_reconstructed_or_abort(
      sprintf("Missing `coef_if` entries for loading keys: %s.", paste(miss, collapse = ", ")),
      "cdid_static_reconstructed_or_channel_error"
    )
  }
  out <- NULL
  row_ids <- ids_chr[rows]
  for (k in lk) {
    J <- loadings[[k]]
    Jm <- if (is.matrix(J)) J else matrix(J, nrow = 1L, dimnames = list(NULL, names(J)))
    L_src <- coef_if[[k]]
    L_sub <- if (is.function(L_src)) {
      L_src(rows)
    } else if (is.matrix(L_src) && is.numeric(L_src) && nrow(L_src) == n) {
      if (!is.null(rownames(L_src)) && !identical(as.character(rownames(L_src)), ids_chr)) {
        static_reconstructed_or_abort(
          sprintf("`coef_if[[\"%s\"]]` rownames do not match `unit_id` in order.", k),
          "cdid_static_reconstructed_or_order_error"
        )
      }
      L_src[rows, , drop = FALSE]
    } else {
      static_reconstructed_or_abort(
        sprintf("`coef_if[[\"%s\"]]` must be an `n x p` numeric matrix or `function(rows)`.", k),
        "cdid_static_reconstructed_or_channel_error"
      )
    }
    if (!is.matrix(L_sub) || !is.numeric(L_sub) || nrow(L_sub) != m ||
        ncol(L_sub) != ncol(Jm) || anyNA(L_sub) || any(!is.finite(L_sub))) {
      static_reconstructed_or_abort(
        sprintf("`coef_if[[\"%s\"]]` slice must be a finite `%d x %d` numeric matrix.", k, m, ncol(Jm)),
        "cdid_static_reconstructed_or_channel_error"
      )
    }
    if (!is.null(rownames(L_sub)) && !identical(as.character(rownames(L_sub)), row_ids)) {
      static_reconstructed_or_abort(
        sprintf("`coef_if[[\"%s\"]]` slice rownames do not match `unit_id[rows]`.", k),
        "cdid_static_reconstructed_or_order_error"
      )
    }
    if (!is.null(colnames(L_sub)) && !is.null(colnames(Jm)) &&
        !identical(colnames(L_sub), colnames(Jm))) {
      static_reconstructed_or_abort(
        sprintf("`coef_if[[\"%s\"]]` column names differ from its loading parameter names.", k),
        "cdid_static_reconstructed_or_channel_error"
      )
    }
    term <- L_sub %*% t(Jm)
    out <- if (is.null(out)) term else out + term
  }
  out
}

static_reconstructed_or_scale_add_loadings <- function(acc, loadings, weight) {
  if (!length(loadings) || weight == 0) return(acc)
  for (k in names(loadings)) {
    val <- weight * loadings[[k]]
    if (is.null(acc[[k]])) {
      acc[[k]] <- val
    } else {
      if (!identical(dim(acc[[k]]), dim(val)) || length(acc[[k]]) != length(val)) {
        static_reconstructed_or_abort(
          sprintf("Incompatible loading dimensions for parameter key `%s`.", k),
          "cdid_static_reconstructed_or_channel_error"
        )
      }
      c_acc <- if (is.matrix(acc[[k]])) colnames(acc[[k]]) else names(acc[[k]])
      c_val <- if (is.matrix(val)) colnames(val) else names(val)
      if (!is.null(c_acc) && !is.null(c_val) && !identical(c_acc, c_val)) {
        static_reconstructed_or_abort(
          sprintf("Incompatible loading parameter names for shared key `%s`.", k),
          "cdid_static_reconstructed_or_channel_error"
        )
      }
      acc[[k]] <- acc[[k]] + val
      if (is.null(c_acc) && !is.null(c_val)) {
        if (is.matrix(acc[[k]])) colnames(acc[[k]]) <- c_val else names(acc[[k]]) <- c_val
      }
    }
  }
  acc
}

# Construct and validate one stayer cell specification for reconstructed OR.
static_reconstructed_or_cell <- function(cell_id,
                                         unit_id,
                                         C,
                                         Z = NULL,
                                         X_tilde = NULL,
                                         D = NULL,
                                         ell = NULL,
                                         edges = NULL,
                                         t = NULL,
                                         baseline_date = NULL,
                                         dates = NULL,
                                         reference_dates = NULL,
                                         d = NULL,
                                         edge_means,
                                         edge_covs = NULL,
                                         coef_if = NULL,
                                         flow = NULL,
                                         or_block = NULL) {
  if (missing(cell_id) || !is.character(cell_id) || !is.null(dim(cell_id)) ||
      is.object(cell_id) || length(cell_id) != 1L || is.na(cell_id) || !nzchar(cell_id)) {
    static_reconstructed_or_abort(
      "`cell_id` must be a non-empty character scalar.",
      "cdid_static_reconstructed_or_cell_error"
    )
  }
  uid <- static_reconstructed_or_validate_unit_id(
    unit_id,
    what = sprintf("cells[[%s]]$unit_id", cell_id)
  )
  n <- uid$n
  C_int <- static_reconstructed_or_validate_mask(C, n, uid$ids_chr, "C", cell_id = cell_id)
  stayer_idx <- which(C_int == 1L)
  n_stayers <- length(stayer_idx)
  if (n_stayers < 1L) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` has zero stayers (`sum(C) == 0`); availability event E_n fails.", cell_id),
      "cdid_static_reconstructed_or_support_error"
    )
  }

  d_val <- NULL
  if (!is.null(d)) {
    if (!static_reconstructed_or_scalar_num(d, finite = TRUE) || !(d %in% c(0, 1))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` field `d` must lie in `{0, 1}` when supplied.", cell_id),
        "cdid_static_reconstructed_or_cell_error"
      )
    }
    d_val <- as.integer(d)
  }

  z_info <- static_reconstructed_or_validate_Z(
    Z = Z, X_tilde = X_tilde, n = n, ids_chr = uid$ids_chr,
    active_rows = stayer_idx, cell_id = cell_id
  )
  p <- z_info$p
  z_names <- z_info$z_names

  g <- static_reconstructed_or_graph(
    D = D, ell = ell, edges = edges, t = t,
    baseline_date = baseline_date, dates = dates,
    reference_dates = reference_dates, cell_id = cell_id
  )
  K <- g$n_edges

  if (missing(edge_means) || !is.list(edge_means) || length(edge_means) != K) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` `edge_means` must be a list of length `nrow(D) = %d`.", cell_id, K),
      "cdid_static_reconstructed_or_edge_error"
    )
  }
  norm_means <- lapply(seq_len(K), function(k) {
    static_reconstructed_or_normalize_edge_moment(
      edge_means[[k]], kind = "mean", n = n, ids_chr = uid$ids_chr,
      coef_if = coef_if, cell_id = cell_id, edge_idx = k
    )
  })

  norm_covs <- vector("list", K)
  if (p > 0L) {
    if (is.null(edge_covs) || !is.list(edge_covs) || length(edge_covs) != K) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` has `p = %d > 0` covariates and requires `edge_covs` of length `nrow(D) = %d`.", cell_id, p, K),
        "cdid_static_reconstructed_or_edge_error"
      )
    }
    norm_covs <- lapply(seq_len(K), function(k) {
      static_reconstructed_or_normalize_edge_moment(
        edge_covs[[k]], kind = "cov", n = n, ids_chr = uid$ids_chr,
        p = p, z_names = z_names, coef_if = coef_if, cell_id = cell_id, edge_idx = k
      )
    })
  } else {
    if (!is.null(edge_covs) && (!is.list(edge_covs) || !(length(edge_covs) %in% c(0L, K)))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` with `p = 0` must omit `edge_covs` or supply an empty/0-column list.", cell_id),
        "cdid_static_reconstructed_or_edge_error"
      )
    }
    norm_covs <- lapply(seq_len(K), function(k) {
      list(
        kind = "cov",
        estimate = stats::setNames(numeric(0L), character(0L)),
        loadings = list(),
        eval_if = function(rows = seq_len(n)) {
          rows <- static_reconstructed_or_validate_rows(rows, n)
          matrix(0.0, nrow = length(rows), ncol = 0L)
        },
        eval_direct = function(rows = seq_len(n)) {
          rows <- static_reconstructed_or_validate_rows(rows, n)
          matrix(0.0, nrow = length(rows), ncol = 0L)
        }
      )
    })
  }

  block_val <- NULL
  if (!is.null(or_block)) {
    if (is.factor(or_block)) or_block <- as.character(or_block)
    if (!is.atomic(or_block) || !is.null(dim(or_block)) || is.object(or_block) ||
        length(or_block) != 1L || is.na(or_block) ||
        (is.numeric(or_block) && !is.finite(or_block)) ||
        (is.character(or_block) && !nzchar(or_block))) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `or_block` must be a non-empty character or finite numeric scalar.", cell_id),
        "cdid_static_reconstructed_or_partition_error"
      )
    }
    block_val <- as.character(or_block)
  }

  structure(
    list(
      cell_id = cell_id,
      t = g$t,
      baseline_date = g$baseline_date,
      d = d_val,
      n = n,
      unit_id = uid$id,
      unit_key = uid$key,
      unit_ids_chr = uid$ids_chr,
      unit_id_is_num = uid$is_num,
      C = C_int,
      stayer_idx = stayer_idx,
      n_stayers = n_stayers,
      Z = z_info$Z,
      p = p,
      z_names = z_names,
      D = g$D,
      ell = g$ell,
      edges = g$edges,
      edge_ids = g$edge_ids,
      n_edges = K,
      dates = g$dates,
      n_components = g$n_components,
      component_id = g$component_id,
      edge_component_id = g$edge_component_id,
      reference_dates = g$reference_dates,
      non_reference_dates = g$non_reference_dates,
      isolated_dates = g$isolated_dates,
      edge_means = norm_means,
      edge_covs = norm_covs,
      coef_if = coef_if,
      flow = flow,
      or_block = block_val
    ),
    class = "cdid_static_reconstructed_or_cell"
  )
}

static_reconstructed_or_normalize_cell <- function(raw_cell, fallback_id = NULL,
                                                   top_unit_id = NULL, top_coef_if = NULL) {
  if (inherits(raw_cell, "cdid_static_reconstructed_or_cell")) {
    if (!is.null(fallback_id) && nzchar(fallback_id) && !identical(raw_cell$cell_id, fallback_id)) {
      static_reconstructed_or_abort(
        sprintf("List name `%s` does not match `cell$cell_id` (`%s`).", fallback_id, raw_cell$cell_id),
        "cdid_static_reconstructed_or_cell_error"
      )
    }
    return(raw_cell)
  }
  if (!is.list(raw_cell) || is.object(raw_cell)) {
    static_reconstructed_or_abort(
      "Every element of `cells` must be a `cdid_static_reconstructed_or_cell` or a plain named list.",
      "cdid_static_reconstructed_or_cell_error"
    )
  }
  cid <- raw_cell$cell_id
  if (is.null(cid)) {
    cid <- fallback_id
  } else if (!is.null(fallback_id) && nzchar(fallback_id) && !identical(cid, fallback_id)) {
    static_reconstructed_or_abort(
      sprintf("List name `%s` does not match element `cell_id` (`%s`).", fallback_id, cid),
      "cdid_static_reconstructed_or_cell_error"
    )
  }
  uid <- if (!is.null(raw_cell$unit_id)) raw_cell$unit_id else top_unit_id
  cif <- if (!is.null(raw_cell$coef_if)) raw_cell$coef_if else top_coef_if
  Z_arg <- if (!is.null(raw_cell$Z_or)) raw_cell$Z_or else raw_cell$Z
  X_arg <- if (!is.null(raw_cell$X_or_tilde)) raw_cell$X_or_tilde else raw_cell$X_tilde
  static_reconstructed_or_cell(
    cell_id = cid,
    unit_id = uid,
    C = raw_cell$C,
    Z = Z_arg,
    X_tilde = X_arg,
    D = raw_cell$D,
    ell = raw_cell$ell,
    edges = raw_cell$edges,
    t = raw_cell$t,
    baseline_date = raw_cell$baseline_date,
    dates = raw_cell$dates,
    reference_dates = raw_cell$reference_dates,
    d = raw_cell$d,
    edge_means = raw_cell$edge_means,
    edge_covs = raw_cell$edge_covs,
    coef_if = cif,
    flow = raw_cell$flow,
    or_block = raw_cell$or_block
  )
}

static_reconstructed_or_validate_partition <- function(cells, or_blocks = NULL) {
  cell_ids <- names(cells)
  n_cells <- length(cell_ids)
  per_cell_blocks <- vapply(cells, function(cl) {
    if (is.null(cl$or_block)) NA_character_ else as.character(cl$or_block)
  }, character(1L))
  has_per_cell <- !all(is.na(per_cell_blocks))
  if (has_per_cell && anyNA(per_cell_blocks)) {
    static_reconstructed_or_abort(
      "When `or_block` is supplied on cells, every declared cell must specify `or_block`.",
      "cdid_static_reconstructed_or_partition_error"
    )
  }
  has_arg <- !is.null(or_blocks)
  if (!has_per_cell && !has_arg) {
    static_reconstructed_or_abort(
      "Supply an explicit OR slope partition (`or_blocks` or per-cell `or_block`); no pooling default is adopted.",
      "cdid_static_reconstructed_or_partition_error"
    )
  }

  parse_arg_blocks <- function(ob) {
    if (is.factor(ob)) ob <- as.character(ob)
    if (is.atomic(ob) && is.null(dim(ob)) && !is.object(ob)) {
      if (length(ob) != n_cells || anyNA(ob) ||
          (is.numeric(ob) && any(!is.finite(ob))) ||
          (is.character(ob) && any(!nzchar(ob)))) {
        static_reconstructed_or_abort(
          sprintf("`or_blocks` vector must have length `%d` with non-empty finite block labels.", n_cells),
          "cdid_static_reconstructed_or_partition_error"
        )
      }
      ob_chr <- as.character(ob)
      nm <- names(ob)
      if (!is.null(nm)) {
        if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) || !setequal(nm, cell_ids)) {
          static_reconstructed_or_abort(
            "Named `or_blocks` vector must have unique names matching the declared `cell_ids`.",
            "cdid_static_reconstructed_or_partition_error"
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
        static_reconstructed_or_abort(
          "Named `or_blocks` partition list must have non-empty unique block names.",
          "cdid_static_reconstructed_or_partition_error"
        )
      }
      seen <- character(0L)
      cell_to_block <- stats::setNames(rep(NA_character_, n_cells), cell_ids)
      for (b_idx in seq_along(ob)) {
        members <- ob[[b_idx]]
        if (is.factor(members)) members <- as.character(members)
        if (!is.atomic(members) || !is.null(dim(members)) || is.object(members) ||
            length(members) < 1L || anyNA(members)) {
          static_reconstructed_or_abort(
            "`or_blocks` partition list elements must be non-empty vectors of declared cell IDs or indices.",
            "cdid_static_reconstructed_or_partition_error"
          )
        }
        if (is.numeric(members)) {
          if (any(!is.finite(members) | members != floor(members) | members < 1 | members > n_cells)) {
            static_reconstructed_or_abort(
              "`or_blocks` numeric cell indices must be integers in `1..length(cells)`.",
              "cdid_static_reconstructed_or_partition_error"
            )
          }
          member_ids <- cell_ids[as.integer(members)]
        } else if (is.character(members)) {
          if (any(!nzchar(members)) || !all(members %in% cell_ids)) {
            static_reconstructed_or_abort(
              "`or_blocks` character members must all be declared `cell_ids`.",
              "cdid_static_reconstructed_or_partition_error"
            )
          }
          member_ids <- members
        } else {
          static_reconstructed_or_abort(
            "`or_blocks` partition list elements must be character cell IDs or integer indices.",
            "cdid_static_reconstructed_or_partition_error"
          )
        }
        if (anyDuplicated(member_ids) || any(member_ids %in% seen)) {
          static_reconstructed_or_abort(
            "`or_blocks` partition list must be disjoint across cells.",
            "cdid_static_reconstructed_or_partition_error"
          )
        }
        b_label <- if (!is.null(list_names)) list_names[b_idx] else paste0("block_", b_idx)
        seen <- c(seen, member_ids)
        cell_to_block[member_ids] <- b_label
      }
      if (!setequal(seen, cell_ids)) {
        static_reconstructed_or_abort(
          "`or_blocks` partition list must cover every declared cell in `cells`.",
          "cdid_static_reconstructed_or_partition_error"
        )
      }
      return(cell_to_block)
    }
    static_reconstructed_or_abort(
      "`or_blocks` must be an explicit per-cell block vector or a disjoint list partitioning the declared cells.",
      "cdid_static_reconstructed_or_partition_error"
    )
  }

  cell_block <- if (has_arg) parse_arg_blocks(or_blocks) else stats::setNames(per_cell_blocks, cell_ids)
  if (has_arg && has_per_cell && !identical(unname(cell_block), unname(per_cell_blocks))) {
    static_reconstructed_or_abort(
      "Supplied `or_blocks` argument conflicts with per-cell `or_block` fields.",
      "cdid_static_reconstructed_or_partition_error"
    )
  }

  block_ids <- unique(unname(cell_block))
  blocks <- lapply(block_ids, function(bid) cell_ids[cell_block == bid])
  names(blocks) <- block_ids

  block_p <- integer(length(block_ids))
  names(block_p) <- block_ids
  block_colnames <- vector("list", length(block_ids))
  names(block_colnames) <- block_ids
  for (bid in block_ids) {
    cids <- blocks[[bid]]
    c1 <- cells[[cids[1L]]]
    p_ref <- c1$p
    names_ref <- c1$z_names
    for (cid in cids[-1L]) {
      ck <- cells[[cid]]
      if (ck$p != p_ref || !identical(ck$z_names, names_ref)) {
        static_reconstructed_or_abort(
          sprintf(
            "OR block `%s` has incompatible slope designs between cells `%s` and `%s`.",
            bid, cids[1L], cid
          ),
          c("cdid_static_reconstructed_or_partition_error", "cdid_or_wls_input_error")
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

# Factor a single cell's known-X stayer covariate design in bounded row chunks
# using `or_pooled_wls_block_summary`, `or_pooled_wls_merge_chunks`, and
# `or_pooled_wls_partial_block`.
# NOTE: A zero dummy vector is passed to `or_pooled_wls_block_summary` SOLELY to
# compute the partialled QR design factor of `(1, Z_c)` on `{C_c = 1}`. It is NOT
# a regression fit to zero outcomes; the actual right-hand side comes from the
# graph-reconstructed stayer moments `m_hat^C_c` and `Cov_hat_c`.
static_reconstructed_or_factor_cell_design <- function(cell, chunk_size) {
  stayer_idx <- cell$stayer_idx
  n_stayers <- cell$n_stayers
  Z_stayer <- cell$Z[stayer_idx, , drop = FALSE]
  csize <- static_reconstructed_or_validate_chunk_size(chunk_size, n_default = n_stayers)

  n_chunks <- as.integer(ceiling(n_stayers / csize))
  chunks <- vector("list", n_chunks)
  for (j in seq_len(n_chunks)) {
    lo <- (j - 1L) * csize + 1L
    hi <- min(n_stayers, j * csize)
    idx_j <- seq.int(lo, hi)
    D_j <- Z_stayer[idx_j, , drop = FALSE]
    # Zero dummy vector used strictly for QR design factorization of known X.
    y_zero_factor_only <- numeric(length(idx_j))
    w_j <- rep(1.0, length(idx_j))
    chunks[[j]] <- or_pooled_wls_block_summary(D = D_j, y = y_zero_factor_only, weights = w_j)
  }
  merged_chunk <- if (n_chunks == 1L) chunks[[1L]] else or_pooled_wls_merge_chunks(chunks)
  partialled <- or_pooled_wls_partial_block(merged_chunk)
  names(partialled$blocks) <- cell$cell_id
  partialled
}

# Reconstruct stayer outcome regression coefficients (alpha_hat_c, b_hat_P) and
# their complete original-unit influence functions across declared cells and
# OR slope partition P.
# All cells in each block P in P_or contribute to b_hat_P and chi^b_P, even when
# a subset `target_cells` is designated for downstream reporting.
static_reconstructed_or_fit <- function(cells,
                                        or_blocks = NULL,
                                        rank_tol,
                                        unit_id = NULL,
                                        coef_if = NULL,
                                        target_cells = NULL,
                                        chunk_size = NULL,
                                        materialize_if = TRUE) {
  rtol <- static_reconstructed_or_validate_rank_tol(rank_tol)
  if (missing(cells) || !is.list(cells) || is.object(cells) || length(cells) < 1L) {
    static_reconstructed_or_abort(
      "Supply a non-empty list of declared `cells`.",
      "cdid_static_reconstructed_or_cell_error"
    )
  }
  if (!is.logical(materialize_if) || length(materialize_if) != 1L || is.na(materialize_if)) {
    static_reconstructed_or_abort(
      "`materialize_if` must be TRUE or FALSE.",
      "cdid_static_reconstructed_or_input_error"
    )
  }

  top_uid <- if (!is.null(unit_id)) static_reconstructed_or_validate_unit_id(unit_id, "unit_id") else NULL
  raw_names <- names(cells)
  if (!is.null(raw_names) && (anyNA(raw_names) || anyDuplicated(raw_names[nzchar(raw_names)]))) {
    static_reconstructed_or_abort(
      "Names of `cells`, when supplied, must not be NA or duplicated.",
      "cdid_static_reconstructed_or_cell_error"
    )
  }

  norm_cells <- vector("list", length(cells))
  for (j in seq_along(cells)) {
    fb <- if (!is.null(raw_names) && nzchar(raw_names[j])) raw_names[j] else NULL
    norm_cells[[j]] <- static_reconstructed_or_normalize_cell(
      cells[[j]],
      fallback_id = fb,
      top_unit_id = if (!is.null(top_uid)) top_uid$id else NULL,
      top_coef_if = coef_if
    )
  }
  cell_ids <- vapply(norm_cells, `[[`, character(1L), "cell_id")
  if (anyDuplicated(cell_ids)) {
    static_reconstructed_or_abort(
      "Declared `cells` must have unique `cell_id` identifiers.",
      "cdid_static_reconstructed_or_cell_error"
    )
  }
  names(norm_cells) <- cell_ids

  ref_uid <- list(
    id = norm_cells[[1L]]$unit_id,
    key = norm_cells[[1L]]$unit_key,
    ids_chr = norm_cells[[1L]]$unit_ids_chr,
    is_num = norm_cells[[1L]]$unit_id_is_num,
    n = norm_cells[[1L]]$n
  )
  if (!is.null(top_uid) && !static_reconstructed_or_same_unit_id(ref_uid, top_uid)) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` `unit_id` does not match top-level `unit_id`.", cell_ids[1L]),
      "cdid_static_reconstructed_or_unit_id_error"
    )
  }
  for (j in seq_along(norm_cells)[-1L]) {
    cj <- norm_cells[[j]]
    uid_j <- list(id = cj$unit_id, key = cj$unit_key, ids_chr = cj$unit_ids_chr,
                  is_num = cj$unit_id_is_num, n = cj$n)
    if (!static_reconstructed_or_same_unit_id(ref_uid, uid_j)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `unit_id` does not match cell `%s` `unit_id`.", cj$cell_id, cell_ids[1L]),
        "cdid_static_reconstructed_or_unit_id_error"
      )
    }
  }

  n <- ref_uid$n
  n_num <- as.double(n)
  ids_chr <- ref_uid$ids_chr
  csize <- static_reconstructed_or_validate_chunk_size(chunk_size, n_default = n)
  partition <- static_reconstructed_or_validate_partition(norm_cells, or_blocks = or_blocks)

  if (is.null(target_cells)) {
    target_cids <- cell_ids
  } else {
    if (is.factor(target_cells)) target_cells <- as.character(target_cells)
    if (!is.atomic(target_cells) || !is.null(dim(target_cells)) || is.object(target_cells) ||
        length(target_cells) < 1L || anyNA(target_cells) ||
        !is.character(target_cells) || any(!nzchar(target_cells)) ||
        anyDuplicated(target_cells) || !all(target_cells %in% cell_ids)) {
      static_reconstructed_or_abort(
        "`target_cells` must be a non-empty character vector of unique declared `cell_ids`.",
        "cdid_static_reconstructed_or_cell_error"
      )
    }
    target_cids <- target_cells
  }

  # Step 1: Per-cell deterministic LS flow varpi_c and reconstructed stayer moments
  # m_hat^C_c and Cov_hat_c, plus known-X stayer quantities pi_c, zbar_c, S_c.
  cell_flows <- vector("list", length(cell_ids))
  names(cell_flows) <- cell_ids
  cell_design_summaries <- vector("list", length(cell_ids))
  names(cell_design_summaries) <- cell_ids
  cell_pieces <- vector("list", length(cell_ids))
  names(cell_pieces) <- cell_ids

  pi_hat <- stats::setNames(numeric(length(cell_ids)), cell_ids)
  mC_hat <- stats::setNames(numeric(length(cell_ids)), cell_ids)
  Cov_hat <- vector("list", length(cell_ids))
  names(Cov_hat) <- cell_ids
  zbar_hat <- vector("list", length(cell_ids))
  names(zbar_hat) <- cell_ids
  S_hat <- vector("list", length(cell_ids))
  names(S_hat) <- cell_ids

  for (cid in cell_ids) {
    cell <- norm_cells[[cid]]
    fl <- static_reconstructed_or_ls_flow(
      D = cell$D, ell = cell$ell, edges = cell$edges, t = cell$t,
      baseline_date = cell$baseline_date, dates = cell$dates,
      reference_dates = cell$reference_dates, rank_tol = rtol,
      flow = cell$flow, cell_id = cid
    )
    cell_flows[[cid]] <- fl
    varpi_c <- fl$varpi
    K_c <- cell$n_edges
    p_c <- cell$p
    z_names_c <- cell$z_names

    edge_m_vec <- vapply(cell$edge_means, `[[`, numeric(1L), "estimate")
    names(edge_m_vec) <- cell$edge_ids
    mC_c <- sum(varpi_c * edge_m_vec)

    loadings_mC_c <- list()
    for (k in seq_len(K_c)) {
      loadings_mC_c <- static_reconstructed_or_scale_add_loadings(
        loadings_mC_c, cell$edge_means[[k]]$loadings, varpi_c[k]
      )
    }

    if (p_c > 0L) {
      edge_cov_mat <- do.call(rbind, lapply(cell$edge_covs, `[[`, "estimate"))
      rownames(edge_cov_mat) <- cell$edge_ids
      colnames(edge_cov_mat) <- z_names_c
      Cov_c <- drop(crossprod(varpi_c, edge_cov_mat))
      names(Cov_c) <- z_names_c
      loadings_Cov_c <- list()
      for (k in seq_len(K_c)) {
        loadings_Cov_c <- static_reconstructed_or_scale_add_loadings(
          loadings_Cov_c, cell$edge_covs[[k]]$loadings, varpi_c[k]
        )
      }
    } else {
      edge_cov_mat <- matrix(0.0, nrow = K_c, ncol = 0L, dimnames = list(cell$edge_ids, character(0L)))
      Cov_c <- stats::setNames(numeric(0L), character(0L))
      loadings_Cov_c <- list()
    }

    # Factor known-X stayer covariates in bounded chunks.
    des_sum <- static_reconstructed_or_factor_cell_design(cell, chunk_size = csize)
    cell_design_summaries[[cid]] <- des_sum
    b_rec <- des_sum$blocks[[cid]]
    pi_c <- cell$n_stayers / n_num
    zbar_c <- stats::setNames(as.double(b_rec$mean_d), z_names_c)
    if (p_c > 0L) {
      R_cell <- des_sum$residual$factor[, seq_len(p_c), drop = FALSE]
      S_c <- crossprod(R_cell) / as.double(cell$n_stayers)
      dimnames(S_c) <- list(z_names_c, z_names_c)
    } else {
      S_c <- matrix(0.0, nrow = 0L, ncol = 0L, dimnames = list(character(0L), character(0L)))
    }

    pi_hat[[cid]] <- pi_c
    mC_hat[[cid]] <- mC_c
    Cov_hat[[cid]] <- Cov_c
    zbar_hat[[cid]] <- zbar_c
    S_hat[[cid]] <- S_c

    # Row-chunk evaluators for cell c:
    eval_chi_mC <- local({
      cl <- cell
      vp <- varpi_c
      Kc <- K_c
      function(rows = seq_len(n)) {
        rows <- static_reconstructed_or_validate_rows(rows, n)
        acc <- numeric(length(rows))
        for (k in seq_len(Kc)) {
          if (vp[k] != 0) {
            acc <- acc + vp[k] * cl$edge_means[[k]]$eval_if(rows)
          }
        }
        acc
      }
    })

    eval_chi_Cov <- local({
      cl <- cell
      vp <- varpi_c
      Kc <- K_c
      pc <- p_c
      zn <- z_names_c
      function(rows = seq_len(n)) {
        rows <- static_reconstructed_or_validate_rows(rows, n)
        m <- length(rows)
        if (pc == 0L) {
          return(matrix(0.0, nrow = m, ncol = 0L, dimnames = list(NULL, character(0L))))
        }
        acc <- matrix(0.0, nrow = m, ncol = pc, dimnames = list(NULL, zn))
        for (k in seq_len(Kc)) {
          if (vp[k] != 0) {
            acc <- acc + vp[k] * cl$edge_covs[[k]]$eval_if(rows)
          }
        }
        acc
      }
    })

    cell_pieces[[cid]] <- list(
      cell_id = cid,
      pi_c = pi_c,
      n_stayers = cell$n_stayers,
      zbar_c = zbar_c,
      S_c = S_c,
      pi_S_c = pi_c * S_c,
      mC_c = mC_c,
      Cov_c = Cov_c,
      pi_Cov_c = pi_c * Cov_c,
      varpi_c = varpi_c,
      edge_means = edge_m_vec,
      edge_covs = edge_cov_mat,
      loadings_mC = loadings_mC_c,
      loadings_Cov = loadings_Cov_c,
      eval_chi_mC = eval_chi_mC,
      eval_chi_Cov = eval_chi_Cov
    )
  }

  # Step 2: Per-block pooled known-X Gram factorization and reconstructed slope solve
  # b_hat_P = A_P^(-1) B_P, followed by cell intercepts alpha_hat_c = mC_c - zbar_c' b_hat_P.
  block_ids <- partition$block_ids
  block_fits <- vector("list", length(block_ids))
  names(block_fits) <- block_fits_names <- block_ids
  A_hat <- vector("list", length(block_ids))
  names(A_hat) <- block_ids
  B_hat <- vector("list", length(block_ids))
  names(B_hat) <- block_ids
  b_hat <- vector("list", length(block_ids))
  names(b_hat) <- block_ids
  alpha_hat <- stats::setNames(numeric(length(cell_ids)), cell_ids)
  slope_by_cell <- vector("list", length(cell_ids))
  names(slope_by_cell) <- cell_ids

  loadings_B <- vector("list", length(block_ids))
  names(loadings_B) <- block_ids
  loadings_b <- vector("list", length(block_ids))
  names(loadings_b) <- block_ids
  loadings_alpha <- vector("list", length(cell_ids))
  names(loadings_alpha) <- cell_ids

  for (bid in block_ids) {
    cids <- partition$blocks[[bid]]
    p_b <- partition$block_p[[bid]]
    z_names_b <- partition$block_colnames[[bid]]

    pooled_des <- if (length(cids) == 1L) {
      cell_design_summaries[[cids[1L]]]
    } else {
      or_pooled_wls_merge(cell_design_summaries[cids])
    }
    # Factor pooled known-X stayer Gram matrix sum_{c in P} N * pi_c * S_c via QR.
    gram_factor <- or_pooled_wls_fit(pooled_des, rank_tol = rtol)

    A_P <- matrix(0.0, nrow = p_b, ncol = p_b, dimnames = list(z_names_b, z_names_b))
    B_P <- stats::setNames(numeric(p_b), z_names_b)
    load_B_P <- list()
    for (cid in cids) {
      cp <- cell_pieces[[cid]]
      if (p_b > 0L) {
        A_P <- A_P + cp$pi_S_c
        B_P <- B_P + cp$pi_Cov_c
        load_B_P <- static_reconstructed_or_scale_add_loadings(
          load_B_P, cp$loadings_Cov, cp$pi_c
        )
      }
    }

    b_P <- if (p_b > 0L) {
      sol <- or_pooled_wls_gram_solve(gram_factor, B_P, normalizer = n_num)
      names(sol) <- z_names_b
      sol
    } else {
      stats::setNames(numeric(0L), character(0L))
    }

    load_b_P <- list()
    if (p_b > 0L && length(load_B_P) > 0L) {
      for (k in names(load_B_P)) {
        Jk <- load_B_P[[k]]
        sol_J <- or_pooled_wls_gram_solve(gram_factor, Jk, normalizer = n_num)
        dimnames(sol_J) <- dimnames(Jk)
        load_b_P[[k]] <- sol_J
      }
    }

    A_hat[[bid]] <- A_P
    B_hat[[bid]] <- B_P
    b_hat[[bid]] <- b_P
    loadings_B[[bid]] <- load_B_P
    loadings_b[[bid]] <- load_b_P

    for (cid in cids) {
      cp <- cell_pieces[[cid]]
      zbar_c <- cp$zbar_c
      alpha_c <- if (p_b > 0L) cp$mC_c - sum(zbar_c * b_P) else cp$mC_c
      alpha_hat[[cid]] <- alpha_c
      slope_by_cell[[cid]] <- b_P

      load_alpha_c <- cp$loadings_mC
      if (p_b > 0L && length(load_b_P) > 0L) {
        for (k in names(load_b_P)) {
          term_k <- stats::setNames(drop(crossprod(zbar_c, load_b_P[[k]])), colnames(load_b_P[[k]]))
          load_alpha_c <- static_reconstructed_or_scale_add_loadings(
            load_alpha_c, stats::setNames(list(term_k), k), -1.0
          )
        }
      }
      loadings_alpha[[cid]] <- load_alpha_c
    }

    block_fits[[bid]] <- list(
      block_id = bid,
      cells = cids,
      p = p_b,
      z_names = z_names_b,
      A_P = A_P,
      B_P = B_P,
      b_P = b_P,
      gram_factor = gram_factor,
      R = gram_factor$R,
      pivot = gram_factor$pivot,
      rank = gram_factor$rank,
      rank_tol = gram_factor$rank_tol,
      cutoff = gram_factor$cutoff,
      loadings_B = load_B_P,
      loadings_b = load_b_P
    )
  }

  # Step 3: Row-chunk IF callback implementing exact original-unit IF identities:
  #   chi_pi = C - pi
  #   chi_z = C * (Z - zbar) / pi
  #   chi_A = sum_{c in P} [C_i (Z_i - zbar_c)(Z_i - zbar_c)' - pi_c S_c]
  #   chi_B = sum_{c in P} [(C_i - pi_c) Cov_c + pi_c chi_Cov_c]
  #   chi_b = A_P^(-1) (chi_B - chi_A b_P)
  #   chi_alpha = chi_mC - chi_z' b_P - zbar_c' chi_b
  if_callback <- function(rows = seq_len(n), include_chi_A_tensor = FALSE) {
    rows <- static_reconstructed_or_validate_rows(rows, n)
    m <- length(rows)
    row_ids <- ids_chr[rows]

    chi_pi_mat <- matrix(0.0, nrow = m, ncol = length(cell_ids), dimnames = list(row_ids, cell_ids))
    chi_mC_mat <- matrix(0.0, nrow = m, ncol = length(cell_ids), dimnames = list(row_ids, cell_ids))
    chi_alpha_mat <- matrix(0.0, nrow = m, ncol = length(cell_ids), dimnames = list(row_ids, cell_ids))

    chi_z_by_cell <- vector("list", length(cell_ids))
    names(chi_z_by_cell) <- cell_ids
    chi_Cov_by_cell <- vector("list", length(cell_ids))
    names(chi_Cov_by_cell) <- cell_ids
    chi_A_b_by_cell <- vector("list", length(cell_ids))
    names(chi_A_b_by_cell) <- cell_ids
    chi_B_by_cell <- vector("list", length(cell_ids))
    names(chi_B_by_cell) <- cell_ids
    chi_A_tensor_by_cell <- if (include_chi_A_tensor) {
      out_l <- vector("list", length(cell_ids))
      names(out_l) <- cell_ids
      out_l
    } else {
      NULL
    }

    chi_A_b_by_block <- vector("list", length(block_ids))
    names(chi_A_b_by_block) <- block_ids
    chi_B_by_block <- vector("list", length(block_ids))
    names(chi_B_by_block) <- block_ids
    chi_b_by_block <- vector("list", length(block_ids))
    names(chi_b_by_block) <- block_ids
    chi_A_tensor_by_block <- if (include_chi_A_tensor) {
      out_b <- vector("list", length(block_ids))
      names(out_b) <- block_ids
      out_b
    } else {
      NULL
    }

    for (bid in block_ids) {
      bf <- block_fits[[bid]]
      cids <- bf$cells
      p_b <- bf$p
      zn_b <- bf$z_names
      b_P <- bf$b_P

      chi_A_b_P <- matrix(0.0, nrow = m, ncol = p_b, dimnames = list(row_ids, zn_b))
      chi_B_P <- matrix(0.0, nrow = m, ncol = p_b, dimnames = list(row_ids, zn_b))
      if (include_chi_A_tensor) {
        chi_A_P_arr <- array(0.0, dim = c(m, p_b, p_b), dimnames = list(row_ids, zn_b, zn_b))
      }

      for (cid in cids) {
        cell <- norm_cells[[cid]]
        cp <- cell_pieces[[cid]]
        pi_c <- cp$pi_c
        zbar_c <- cp$zbar_c
        S_c <- cp$S_c
        Cov_c <- cp$Cov_c

        C_sub <- as.double(cell$C[rows])
        chi_pi_c <- C_sub - pi_c
        chi_pi_mat[, cid] <- chi_pi_c

        chi_mC_c <- cp$eval_chi_mC(rows)
        chi_mC_mat[, cid] <- chi_mC_c

        chi_Cov_c <- cp$eval_chi_Cov(rows)
        rownames(chi_Cov_c) <- row_ids
        chi_Cov_by_cell[[cid]] <- chi_Cov_c

        if (p_b > 0L) {
          # Only read Z on rows where C_sub == 1 so non-stayer rows never propagate NA.
          act_local <- which(C_sub == 1)
          Zc_sub <- matrix(0.0, nrow = m, ncol = p_b, dimnames = list(row_ids, zn_b))
          if (length(act_local) > 0L) {
            Z_act <- cell$Z[rows[act_local], , drop = FALSE]
            Zc_sub[act_local, ] <- sweep(Z_act, 2L, zbar_c, FUN = "-")
          }
          chi_z_c <- Zc_sub / pi_c
          chi_z_by_cell[[cid]] <- chi_z_c

          # chi_A_c * b_P = C_i (Z_i - zbar_c)((Z_i - zbar_c)' b_P) - pi_c S_c b_P
          proj_b <- drop(Zc_sub %*% b_P)
          pi_S_b <- drop(cp$pi_S_c %*% b_P)
          chi_A_b_c <- sweep(Zc_sub * proj_b, 2L, pi_S_b, FUN = "-")
          chi_A_b_by_cell[[cid]] <- chi_A_b_c
          chi_A_b_P <- chi_A_b_P + chi_A_b_c

          # chi_B_c = (C_i - pi_c) Cov_c + pi_c chi_Cov_c
          chi_B_c <- outer(chi_pi_c, Cov_c) + pi_c * chi_Cov_c
          dimnames(chi_B_c) <- list(row_ids, zn_b)
          chi_B_by_cell[[cid]] <- chi_B_c
          chi_B_P <- chi_B_P + chi_B_c

          if (include_chi_A_tensor) {
            chi_A_c_arr <- array(0.0, dim = c(m, p_b, p_b), dimnames = list(row_ids, zn_b, zn_b))
            for (ii in seq_len(m)) {
              if (C_sub[ii] == 1) {
                z_ii <- Zc_sub[ii, ]
                chi_A_c_arr[ii, , ] <- tcrossprod(z_ii) - cp$pi_S_c
              } else {
                chi_A_c_arr[ii, , ] <- -cp$pi_S_c
              }
            }
            chi_A_tensor_by_cell[[cid]] <- chi_A_c_arr
            chi_A_P_arr <- chi_A_P_arr + chi_A_c_arr
          }
        } else {
          empty_mat <- matrix(0.0, nrow = m, ncol = 0L, dimnames = list(row_ids, character(0L)))
          chi_z_by_cell[[cid]] <- empty_mat
          chi_A_b_by_cell[[cid]] <- empty_mat
          chi_B_by_cell[[cid]] <- empty_mat
          if (include_chi_A_tensor) {
            chi_A_tensor_by_cell[[cid]] <- array(0.0, dim = c(m, 0L, 0L))
          }
        }
      }

      # Solve A_P^(-1) (chi_B_P - chi_A_P b_P) via compact QR factor without matrix inversion.
      if (p_b > 0L && m > 0L) {
        rhs_P <- chi_B_P - chi_A_b_P
        sol_P <- t(or_pooled_wls_gram_solve(bf$gram_factor, t(rhs_P), normalizer = n_num))
        if (p_b == 1L) sol_P <- matrix(sol_P, nrow = m, ncol = 1L)
        dimnames(sol_P) <- list(row_ids, zn_b)
        chi_b_P <- sol_P
      } else {
        chi_b_P <- matrix(0.0, nrow = m, ncol = p_b, dimnames = list(row_ids, zn_b))
      }

      chi_A_b_by_block[[bid]] <- chi_A_b_P
      chi_B_by_block[[bid]] <- chi_B_P
      chi_b_by_block[[bid]] <- chi_b_P
      if (include_chi_A_tensor) {
        chi_A_tensor_by_block[[bid]] <- chi_A_P_arr
      }

      for (cid in cids) {
        cp <- cell_pieces[[cid]]
        zbar_c <- cp$zbar_c
        chi_mC_c <- chi_mC_mat[, cid]
        if (p_b > 0L) {
          chi_z_c <- chi_z_by_cell[[cid]]
          term_z <- drop(chi_z_c %*% b_P)
          term_b <- drop(chi_b_P %*% zbar_c)
          chi_alpha_mat[, cid] <- chi_mC_c - term_z - term_b
        } else {
          chi_alpha_mat[, cid] <- chi_mC_c
        }
      }
    }

    chi_b_by_cell <- lapply(cell_ids, function(cid) {
      chi_b_by_block[[partition$cell_block[[cid]]]]
    })
    names(chi_b_by_cell) <- cell_ids

    list(
      rows = rows,
      unit_id = ref_uid$id[rows],
      chi_pi = chi_pi_mat,
      chi_z = chi_z_by_cell,
      chi_mC = chi_mC_mat,
      chi_Cov = chi_Cov_by_cell,
      chi_A_b_cell = chi_A_b_by_cell,
      chi_B_cell = chi_B_by_cell,
      chi_A_b = chi_A_b_by_block,
      chi_B = chi_B_by_block,
      chi_b = chi_b_by_block,
      chi_b_by_cell = chi_b_by_cell,
      chi_alpha = chi_alpha_mat,
      chi_A_cell = chi_A_tensor_by_cell,
      chi_A = chi_A_tensor_by_block
    )
  }

  full_if <- if (materialize_if) if_callback(seq_len(n), include_chi_A_tensor = FALSE) else NULL

  structure(
    list(
      n = n,
      unit_id = ref_uid$id,
      cell_ids = cell_ids,
      target_cells = target_cids,
      partition = partition,
      rank_tol = rtol,
      chunk_size = csize,
      pi_hat = pi_hat,
      zbar_hat = zbar_hat,
      S_hat = S_hat,
      mC_hat = mC_hat,
      Cov_hat = Cov_hat,
      A_hat = A_hat,
      B_hat = B_hat,
      b_hat = b_hat,
      slope_by_cell = slope_by_cell,
      alpha_hat = alpha_hat,
      flows = cell_flows,
      cells = cell_pieces,
      blocks = block_fits,
      loadings = list(
        B = loadings_B,
        b = loadings_b,
        alpha = loadings_alpha
      ),
      if_callback = if_callback,
      influence = full_if,
      chi_pi = if (materialize_if) full_if$chi_pi else NULL,
      chi_z = if (materialize_if) full_if$chi_z else NULL,
      chi_mC = if (materialize_if) full_if$chi_mC else NULL,
      chi_Cov = if (materialize_if) full_if$chi_Cov else NULL,
      chi_A_b = if (materialize_if) full_if$chi_A_b else NULL,
      chi_B = if (materialize_if) full_if$chi_B else NULL,
      chi_b = if (materialize_if) full_if$chi_b else NULL,
      chi_b_by_cell = if (materialize_if) full_if$chi_b_by_cell else NULL,
      chi_alpha = if (materialize_if) full_if$chi_alpha else NULL
    ),
    class = "cdid_static_reconstructed_or_fit"
  )
}

# Construct an explicit full-role design mean object {estimate, eval_if, influence}
# for switcher role G (`zG = E_G[Z_c]`) or odds-tilted stayer role C (`zCr = E_{C,r}[Z_c]`)
# when caller wishes to form the mean-Z primitive from masks/odds/ell_gamma
# (with influence functions chi^{zG}_c and chi^{zr}_c).
static_reconstructed_or_role_mean_z <- function(unit_id,
                                                A,
                                                Z,
                                                r = NULL,
                                                X_tilde = NULL,
                                                active = NULL,
                                                ell_gamma = NULL,
                                                materialize_if = TRUE) {
  uid <- static_reconstructed_or_validate_unit_id(unit_id, "unit_id")
  n <- uid$n
  n_num <- as.double(n)
  ids_chr <- uid$ids_chr
  A_int <- static_reconstructed_or_validate_mask(A, n, ids_chr, "A")
  act_rows <- which(A_int == 1L)
  if (!length(act_rows)) {
    static_reconstructed_or_abort(
      "Role mask `A` is empty (`sum(A) == 0`).",
      "cdid_static_reconstructed_or_support_error"
    )
  }
  z_info <- static_reconstructed_or_validate_Z(
    Z = Z, X_tilde = NULL, n = n, ids_chr = ids_chr,
    active_rows = act_rows, cell_id = "<role_mean_z>"
  )
  Zm <- z_info$Z
  p <- z_info$p
  z_names <- z_info$z_names

  if (is.null(r)) {
    w_act <- rep(1.0, length(act_rows))
  } else {
    if (!is.atomic(r) || !is.null(dim(r)) || is.object(r) || !is.numeric(r) || length(r) != n) {
      static_reconstructed_or_abort(
        "`r` must be a numeric vector of length `length(unit_id)`.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
    if (!is.null(names(r)) && !identical(as.character(names(r)), ids_chr)) {
      static_reconstructed_or_abort(
        "Names of `r` do not match `unit_id` in order.",
        "cdid_static_reconstructed_or_order_error"
      )
    }
    w_act <- as.double(r[act_rows])
    if (anyNA(w_act) || any(!is.finite(w_act) | w_act < 0)) {
      static_reconstructed_or_abort(
        "`r` must be finite and nonnegative on role units.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
  }
  S_w <- sum(w_act)
  if (!is.finite(S_w) || S_w <= 0) {
    static_reconstructed_or_abort(
      "Role weight denominator `sum(A * r)` must be strictly positive.",
      "cdid_static_reconstructed_or_support_error"
    )
  }
  denom_mean <- S_w / n_num

  if (p == 0L) {
    est <- stats::setNames(numeric(0L), character(0L))
    eval_if <- function(rows = seq_len(n)) {
      rows <- static_reconstructed_or_validate_rows(rows, n)
      matrix(0.0, nrow = length(rows), ncol = 0L, dimnames = list(ids_chr[rows], character(0L)))
    }
    return(list(
      estimate = est,
      denominator_mean = denom_mean,
      loadings_gamma = NULL,
      eval_if = eval_if,
      influence = if (materialize_if) eval_if(seq_len(n)) else NULL
    ))
  }

  Z_act <- Zm[act_rows, , drop = FALSE]
  z_bar <- stats::setNames(colSums(Z_act * w_act) / S_w, z_names)
  Zc_act <- sweep(Z_act, 2L, z_bar, FUN = "-")

  # Optional treatment-logit loading J^gamma_zr = E_n[C_c (Z_c - zbar^r_c) rdot_c'] / Q_c
  J_gamma <- NULL
  has_logit <- !is.null(ell_gamma) || !is.null(X_tilde) || !is.null(active)
  if (has_logit) {
    if (is.null(r) || is.null(X_tilde) || is.null(active) || is.null(ell_gamma)) {
      static_reconstructed_or_abort(
        "Supply all of `r`, `X_tilde`, `active`, and `ell_gamma` when including treatment-logit IF in `zCr`.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
    if (!is.matrix(X_tilde) || !is.numeric(X_tilde) || nrow(X_tilde) != n ||
        !is.matrix(ell_gamma) || !is.numeric(ell_gamma) || nrow(ell_gamma) != n ||
        ncol(X_tilde) != ncol(ell_gamma)) {
      static_reconstructed_or_abort(
        "`X_tilde` and `ell_gamma` must be numeric `n x k` matrices with matching dimensions.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
    if (!is.null(rownames(X_tilde)) && !identical(as.character(rownames(X_tilde)), ids_chr)) {
      static_reconstructed_or_abort(
        "`X_tilde` rownames do not match `unit_id` in order.",
        "cdid_static_reconstructed_or_order_error"
      )
    }
    if (!is.null(rownames(ell_gamma)) && !identical(as.character(rownames(ell_gamma)), ids_chr)) {
      static_reconstructed_or_abort(
        "`ell_gamma` rownames do not match `unit_id` in order.",
        "cdid_static_reconstructed_or_order_error"
      )
    }
    if (!is.null(colnames(X_tilde)) && !is.null(colnames(ell_gamma)) &&
        !identical(colnames(ell_gamma), colnames(X_tilde))) {
      static_reconstructed_or_abort(
        "`ell_gamma` column names differ from `X_tilde` column names.",
        "cdid_static_reconstructed_or_channel_error"
      )
    }
    if (!is.atomic(active) || !is.null(dim(active)) || is.object(active) || length(active) != n) {
      static_reconstructed_or_abort(
        "`active` must be a vector of length `length(unit_id)`.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
    if (!is.null(names(active)) && !identical(as.character(names(active)), ids_chr)) {
      static_reconstructed_or_abort(
        "Names of `active` do not match `unit_id` in order.",
        "cdid_static_reconstructed_or_order_error"
      )
    }
    act_mask <- as.double(active[act_rows])
    if (anyNA(act_mask) || any(!(act_mask %in% c(0, 1)))) {
      static_reconstructed_or_abort(
        "`active` must be binary `{0, 1}` on active role units.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
    X_act <- X_tilde[act_rows, , drop = FALSE]
    if (anyNA(X_act) || any(!is.finite(X_act)) || anyNA(ell_gamma) || any(!is.finite(ell_gamma))) {
      static_reconstructed_or_abort(
        "`X_tilde` (on active rows) and `ell_gamma` must be finite.",
        "cdid_static_reconstructed_or_input_error"
      )
    }
    # J_gamma is p x k: sum_{i in act} (Z_ic - zbar^r_c) (act_ic * r_ic * xtilde_ic') / S_w
    J_gamma <- crossprod(Zc_act, X_act * (act_mask * w_act)) / S_w
    dimnames(J_gamma) <- list(z_names, colnames(X_tilde))
  }

  direct_full <- matrix(0.0, nrow = n, ncol = p, dimnames = list(ids_chr, z_names))
  direct_full[act_rows, ] <- (Zc_act * w_act) / denom_mean

  eval_if <- function(rows = seq_len(n)) {
    rows <- static_reconstructed_or_validate_rows(rows, n)
    out <- direct_full[rows, , drop = FALSE]
    if (!is.null(J_gamma)) {
      out <- out + ell_gamma[rows, , drop = FALSE] %*% t(J_gamma)
    }
    out
  }

  list(
    estimate = z_bar,
    denominator_mean = denom_mean,
    loadings_gamma = J_gamma,
    direct_if = direct_full,
    eval_if = eval_if,
    influence = if (materialize_if) eval_if(seq_len(n)) else NULL
  )
}

static_reconstructed_or_normalize_mean_z <- function(z_obj, p, z_names, n, ids_chr, what) {
  if (p == 0L && is.null(z_obj)) {
    return(list(
      estimate = stats::setNames(numeric(0L), character(0L)),
      eval_if = function(rows = seq_len(n)) {
        rows <- static_reconstructed_or_validate_rows(rows, n)
        matrix(0.0, nrow = length(rows), ncol = 0L, dimnames = list(ids_chr[rows], character(0L)))
      }
    ))
  }
  if (is.null(z_obj) || !is.list(z_obj)) {
    static_reconstructed_or_abort(
      sprintf("`%s` must be an explicit mean-Z list containing `$estimate` and `$influence` (or `$eval_if` / `$if_fn`).", what),
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  norm <- static_reconstructed_or_normalize_edge_moment(
    z_obj, kind = "cov", n = n, ids_chr = ids_chr,
    p = p, z_names = z_names, cell_id = what, edge_idx = 1L
  )
  list(estimate = norm$estimate, eval_if = norm$eval_if)
}

# Fixed switch-step link centering.
# Consumes a fitted `cdid_static_reconstructed_or_fit` (or explicit alpha/b + IFs),
# explicit full-role design mean objects `zG` and (for DR) `zCr`, and role edge
# link specifications for one cell `c = (t, d)`.
# Edge `(s, u)` with `s < u` is adjusted by `[1{u >= t} - 1{s >= t}]` times the
# role center `alpha_c + z_role' b_{P(c)}` with the exact IF chain rule:
#   chi^{c, role}_c = chi^alpha_c + z_role' chi^b_{P(c)} + b_{P(c)}' chi^{z_role}_c.
# Methods:
#   - "ipw": both G and C links retained, uncentered (adjustment == 0)
#   - "or" / "reg": only centered G links enter final solver
#   - "dr": both centered G and centered C links enter final solver
static_reconstructed_or_center_cell_links <- function(or_fit,
                                                      cell_id,
                                                      method,
                                                      t = NULL,
                                                      G_edges = NULL,
                                                      G_links = NULL,
                                                      C_edges = NULL,
                                                      C_links = NULL,
                                                      zG = NULL,
                                                      zCr = NULL,
                                                      date_shift_G = NULL,
                                                      date_shift_C = NULL,
                                                      materialize_if = TRUE) {
  if (!inherits(or_fit, "cdid_static_reconstructed_or_fit")) {
    static_reconstructed_or_abort(
      "`or_fit` must be a `cdid_static_reconstructed_or_fit` object from `static_reconstructed_or_fit()`.",
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  if (missing(cell_id) || !is.character(cell_id) || length(cell_id) != 1L ||
      !(cell_id %in% or_fit$cell_ids)) {
    static_reconstructed_or_abort(
      "`cell_id` must name a cell present in `or_fit$cell_ids`.",
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  if (missing(method) || !is.character(method) || length(method) != 1L ||
      is.na(method) || !(method %in% c("ipw", "or", "reg", "dr"))) {
    static_reconstructed_or_abort(
      "Supply an explicit `method` in `c(\"ipw\", \"or\", \"reg\", \"dr\")`.",
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  paper_method <- if (identical(method, "reg")) "or" else method
  canonical_method <- if (identical(method, "or")) "reg" else method

  n <- or_fit$n
  ids_chr <- as.character(or_fit$unit_id)
  bid <- or_fit$partition$cell_block[[cell_id]]
  p <- or_fit$partition$block_p[[bid]]
  z_names <- or_fit$partition$block_colnames[[bid]]
  alpha_c <- or_fit$alpha_hat[[cell_id]]
  b_P <- or_fit$b_hat[[bid]]

  t_cell <- or_fit$flows[[cell_id]]$t
  if (is.null(t)) t <- t_cell
  if (!static_reconstructed_or_scalar_num(t, finite = TRUE)) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` centering requires a finite numeric transition date `t`.", cell_id),
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  t_val <- as.double(t)
  if (!is.null(t_cell) && t_val != t_cell) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` supplied `t = %g` conflicts with cell transition date `%g`.", cell_id, t_val, t_cell),
      "cdid_static_reconstructed_or_centering_error"
    )
  }

  if (is.null(G_edges) || is.null(G_links)) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` requires switcher `G_edges` and `G_links`.", cell_id),
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  g_edge_info <- static_reconstructed_or_parse_edges(G_edges, cell_id = paste0(cell_id, ":G"))
  K_G <- g_edge_info$n_edges
  if (!is.list(G_links) || length(G_links) != K_G) {
    static_reconstructed_or_abort(
      sprintf("Cell `%s` `G_links` must be a list of length `nrow(G_edges) = %d`.", cell_id, K_G),
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  norm_G_links <- lapply(seq_len(K_G), function(k) {
    static_reconstructed_or_normalize_edge_moment(
      G_links[[k]], kind = "mean", n = n, ids_chr = ids_chr,
      cell_id = paste0(cell_id, ":G"), edge_idx = k
    )
  })
  # Fixed switch-step multiplier [1{u >= t} - 1{s >= t}] (never scaled by gap length u - s).
  step_G <- as.double(g_edge_info$u >= t_val) - as.double(g_edge_info$s >= t_val)
  names(step_G) <- g_edge_info$edge_ids

  uses_C <- paper_method %in% c("ipw", "dr")
  if (!uses_C && (!is.null(C_edges) || !is.null(C_links) || !is.null(zCr))) {
    static_reconstructed_or_abort(
      sprintf("Method `%s` uses only switcher centered links; do not pass `C_edges`, `C_links`, or `zCr`.", method),
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  if (uses_C) {
    if (is.null(C_edges) || is.null(C_links)) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` with method `%s` requires stayer `C_edges` and `C_links`.", cell_id, method),
        "cdid_static_reconstructed_or_centering_error"
      )
    }
    c_edge_info <- static_reconstructed_or_parse_edges(C_edges, cell_id = paste0(cell_id, ":C"))
    K_C <- c_edge_info$n_edges
    if (!is.list(C_links) || length(C_links) != K_C) {
      static_reconstructed_or_abort(
        sprintf("Cell `%s` `C_links` must be a list of length `nrow(C_edges) = %d`.", cell_id, K_C),
        "cdid_static_reconstructed_or_centering_error"
      )
    }
    norm_C_links <- lapply(seq_len(K_C), function(k) {
      static_reconstructed_or_normalize_edge_moment(
        C_links[[k]], kind = "mean", n = n, ids_chr = ids_chr,
        cell_id = paste0(cell_id, ":C"), edge_idx = k
      )
    })
    step_C <- as.double(c_edge_info$u >= t_val) - as.double(c_edge_info$s >= t_val)
    names(step_C) <- c_edge_info$edge_ids
  } else {
    c_edge_info <- NULL
    K_C <- 0L
    norm_C_links <- NULL
    step_C <- NULL
  }

  # Optional internal non-switch date shift for testing non-switch centering invariance:
  # h_{c,s} = 1{s >= t} center + shift(s) with shift(t) - shift(t - 1) == 0.
  extra_edge_shift <- function(edge_info, date_shift, role_label) {
    if (is.null(date_shift)) return(list(est = numeric(edge_info$n_edges), eval_if = NULL))
    if (!is.list(date_shift) || is.null(date_shift$dates) || is.null(date_shift$values)) {
      static_reconstructed_or_abort(
        sprintf("`date_shift_%s` must be a list with `$dates`, `$values`, and optional `$if_mat`.", role_label),
        "cdid_static_reconstructed_or_centering_error"
      )
    }
    ds_dates <- as.double(date_shift$dates)
    ds_vals <- as.double(date_shift$values)
    idx_t <- match(t_val, ds_dates, nomatch = 0L)
    idx_tm1 <- match(t_val - 1.0, ds_dates, nomatch = 0L)
    val_t <- if (idx_t > 0L) ds_vals[idx_t] else 0.0
    val_tm1 <- if (idx_tm1 > 0L) ds_vals[idx_tm1] else 0.0
    if (abs(val_t - val_tm1) > 0) {
      static_reconstructed_or_abort(
        sprintf("`date_shift_%s` must have zero increment across `(t - 1, t)` (prop:invariance).", role_label),
        "cdid_static_reconstructed_or_centering_error"
      )
    }
    u_pos <- match(edge_info$u, ds_dates, nomatch = 0L)
    s_pos <- match(edge_info$s, ds_dates, nomatch = 0L)
    u_val <- ifelse(u_pos > 0L, ds_vals[u_pos], 0.0)
    s_val <- ifelse(s_pos > 0L, ds_vals[s_pos], 0.0)
    est_diff <- u_val - s_val
    if_mat <- date_shift$if_mat
    eval_shift_if <- if (!is.null(if_mat)) {
      function(rows) {
        sub <- if_mat[rows, , drop = FALSE]
        out <- matrix(0.0, nrow = length(rows), ncol = edge_info$n_edges)
        for (k in seq_len(edge_info$n_edges)) {
          if (u_pos[k] > 0L) out[, k] <- out[, k] + sub[, u_pos[k]]
          if (s_pos[k] > 0L) out[, k] <- out[, k] - sub[, s_pos[k]]
        }
        out
      }
    } else {
      NULL
    }
    list(est = est_diff, eval_if = eval_shift_if)
  }

  shift_G <- extra_edge_shift(g_edge_info, date_shift_G, "G")
  shift_C <- if (uses_C) extra_edge_shift(c_edge_info, date_shift_C, "C") else NULL

  raw_G <- vapply(norm_G_links, `[[`, numeric(1L), "estimate")
  names(raw_G) <- g_edge_info$edge_ids
  raw_C <- if (uses_C) {
    rc <- vapply(norm_C_links, `[[`, numeric(1L), "estimate")
    names(rc) <- c_edge_info$edge_ids
    rc
  } else {
    NULL
  }

  if (paper_method == "ipw") {
    if (!is.null(zG) || !is.null(zCr)) {
      static_reconstructed_or_abort(
        "Method `ipw` is uncentered; do not pass `zG` or `zCr`.",
        "cdid_static_reconstructed_or_centering_error"
      )
    }
    center_G <- 0.0
    center_C <- 0.0
    adj_G <- raw_G
    adj_C <- raw_C
    zG_norm <- NULL
    zCr_norm <- NULL
  } else if (paper_method == "or") {
    zG_norm <- static_reconstructed_or_normalize_mean_z(zG, p, z_names, n, ids_chr, "zG")
    zCr_norm <- NULL
    center_G <- if (p > 0L) alpha_c + sum(zG_norm$estimate * b_P) else alpha_c
    center_C <- NA_real_
    adj_G <- raw_G - step_G * center_G - shift_G$est
    adj_C <- NULL
  } else {
    # paper_method == "dr"
    zG_norm <- static_reconstructed_or_normalize_mean_z(zG, p, z_names, n, ids_chr, "zG")
    zCr_norm <- static_reconstructed_or_normalize_mean_z(zCr, p, z_names, n, ids_chr, "zCr")
    center_G <- if (p > 0L) alpha_c + sum(zG_norm$estimate * b_P) else alpha_c
    center_C <- if (p > 0L) alpha_c + sum(zCr_norm$estimate * b_P) else alpha_c
    adj_G <- raw_G - step_G * center_G - shift_G$est
    adj_C <- raw_C - step_C * center_C - shift_C$est
  }

  stacked_links <- if (uses_C) {
    c(stats::setNames(adj_G, paste0("G:", names(adj_G))),
      stats::setNames(adj_C, paste0("C:", names(adj_C))))
  } else {
    stats::setNames(adj_G, paste0("G:", names(adj_G)))
  }

  if_callback <- function(rows = seq_len(n)) {
    rows <- static_reconstructed_or_validate_rows(rows, n)
    m <- length(rows)
    row_ids <- ids_chr[rows]

    chi_raw_G <- matrix(0.0, nrow = m, ncol = K_G, dimnames = list(row_ids, g_edge_info$edge_ids))
    for (k in seq_len(K_G)) {
      chi_raw_G[, k] <- norm_G_links[[k]]$eval_if(rows)
    }
    chi_raw_C <- if (uses_C) {
      mat_c <- matrix(0.0, nrow = m, ncol = K_C, dimnames = list(row_ids, c_edge_info$edge_ids))
      for (k in seq_len(K_C)) {
        mat_c[, k] <- norm_C_links[[k]]$eval_if(rows)
      }
      mat_c
    } else {
      NULL
    }

    if (paper_method == "ipw") {
      chi_center_G <- numeric(m)
      chi_center_C <- numeric(m)
      chi_adj_G <- chi_raw_G
      chi_adj_C <- chi_raw_C
    } else {
      or_slice <- if (!is.null(or_fit$influence) && identical(rows, seq_len(n))) {
        or_fit$influence
      } else {
        or_fit$if_callback(rows)
      }
      chi_alpha_c <- or_slice$chi_alpha[, cell_id]
      chi_b_P <- or_slice$chi_b[[bid]]

      if (p > 0L) {
        chi_zG <- zG_norm$eval_if(rows)
        chi_center_G <- chi_alpha_c + drop(chi_b_P %*% zG_norm$estimate) + drop(chi_zG %*% b_P)
      } else {
        chi_center_G <- chi_alpha_c
      }
      chi_adj_G <- chi_raw_G - outer(chi_center_G, step_G)
      if (!is.null(shift_G$eval_if)) {
        chi_adj_G <- chi_adj_G - shift_G$eval_if(rows)
      }

      if (paper_method == "dr") {
        if (p > 0L) {
          chi_zCr <- zCr_norm$eval_if(rows)
          chi_center_C <- chi_alpha_c + drop(chi_b_P %*% zCr_norm$estimate) + drop(chi_zCr %*% b_P)
        } else {
          chi_center_C <- chi_alpha_c
        }
        chi_adj_C <- chi_raw_C - outer(chi_center_C, step_C)
        if (!is.null(shift_C$eval_if)) {
          chi_adj_C <- chi_adj_C - shift_C$eval_if(rows)
        }
      } else {
        chi_center_C <- NULL
        chi_adj_C <- NULL
      }
    }

    chi_stacked <- if (uses_C) {
      mat <- cbind(chi_adj_G, chi_adj_C)
      colnames(mat) <- names(stacked_links)
      mat
    } else {
      mat <- chi_adj_G
      colnames(mat) <- names(stacked_links)
      mat
    }

    list(
      rows = rows,
      unit_id = or_fit$unit_id[rows],
      chi_raw_G = chi_raw_G,
      chi_raw_C = chi_raw_C,
      chi_center_G = chi_center_G,
      chi_center_C = chi_center_C,
      chi_adj_G = chi_adj_G,
      chi_adj_C = chi_adj_C,
      chi_stacked = chi_stacked
    )
  }

  full_if <- if (materialize_if) if_callback(seq_len(n)) else NULL

  structure(
    list(
      cell_id = cell_id,
      block_id = bid,
      method = canonical_method,
      paper_method = paper_method,
      t = t_val,
      alpha_c = alpha_c,
      b_P = b_P,
      zG = if (!is.null(zG_norm)) zG_norm$estimate else NULL,
      zCr = if (!is.null(zCr_norm)) zCr_norm$estimate else NULL,
      center_G = center_G,
      center_C = center_C,
      G_edges = g_edge_info$edges,
      C_edges = if (uses_C) c_edge_info$edges else NULL,
      step_G = step_G,
      step_C = step_C,
      raw_G = raw_G,
      raw_C = raw_C,
      adjusted_G = adj_G,
      adjusted_C = adj_C,
      stacked_links = stacked_links,
      if_callback = if_callback,
      influence = full_if,
      chi_center_G = if (materialize_if) full_if$chi_center_G else NULL,
      chi_center_C = if (materialize_if) full_if$chi_center_C else NULL,
      chi_adjusted_G = if (materialize_if) full_if$chi_adj_G else NULL,
      chi_adjusted_C = if (materialize_if) full_if$chi_adj_C else NULL,
      chi_stacked = if (materialize_if) full_if$chi_stacked else NULL
    ),
    class = "cdid_static_reconstructed_or_centered_links"
  )
}

# Multi-cell fixed switch-step centering wrapper across a requested cell set
# A <= or_fit$cell_ids. Stacks adjusted links in cell order (G links, and for
# IPW/DR also C links) and exposes a row-chunk IF callback `if_callback(rows)`.
static_reconstructed_or_center_links <- function(or_fit,
                                                 cell_links,
                                                 method,
                                                 materialize_if = TRUE) {
  if (!inherits(or_fit, "cdid_static_reconstructed_or_fit")) {
    static_reconstructed_or_abort(
      "`or_fit` must be a `cdid_static_reconstructed_or_fit` object.",
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  if (missing(cell_links) || !is.list(cell_links) || is.object(cell_links) || !length(cell_links)) {
    static_reconstructed_or_abort(
      "Supply a non-empty named list `cell_links` of per-cell link specifications.",
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  cids <- names(cell_links)
  if (is.null(cids) || anyNA(cids) || any(!nzchar(cids)) || anyDuplicated(cids) ||
      !all(cids %in% or_fit$cell_ids)) {
    static_reconstructed_or_abort(
      "`cell_links` must have unique names belonging to `or_fit$cell_ids`.",
      "cdid_static_reconstructed_or_centering_error"
    )
  }

  per_cell <- lapply(cids, function(cid) {
    spec <- cell_links[[cid]]
    if (inherits(spec, "cdid_static_reconstructed_or_centered_links")) {
      if (!identical(spec$cell_id, cid)) {
        static_reconstructed_or_abort(
          sprintf("Element `%s` has mismatched `cell_id = %s`.", cid, spec$cell_id),
          "cdid_static_reconstructed_or_centering_error"
        )
      }
      return(spec)
    }
    if (!is.list(spec)) {
      static_reconstructed_or_abort(
        sprintf("`cell_links[[%s]]` must be a list.", cid),
        "cdid_static_reconstructed_or_centering_error"
      )
    }
    static_reconstructed_or_center_cell_links(
      or_fit = or_fit,
      cell_id = cid,
      method = method,
      t = spec$t,
      G_edges = spec$G_edges,
      G_links = spec$G_links,
      C_edges = spec$C_edges,
      C_links = spec$C_links,
      zG = spec$zG,
      zCr = spec$zCr,
      date_shift_G = spec$date_shift_G,
      date_shift_C = spec$date_shift_C,
      materialize_if = materialize_if
    )
  })
  names(per_cell) <- cids

  stacked_list <- lapply(cids, function(cid) {
    v <- per_cell[[cid]]$stacked_links
    stats::setNames(v, paste0(cid, "|", names(v)))
  })
  stacked_links <- unlist(stacked_list, use.names = TRUE)

  n <- or_fit$n
  if_callback <- function(rows = seq_len(n)) {
    rows <- static_reconstructed_or_validate_rows(rows, n)
    slices <- lapply(cids, function(cid) per_cell[[cid]]$if_callback(rows))
    names(slices) <- cids
    chi_mat <- do.call(cbind, lapply(cids, function(cid) {
      m_c <- slices[[cid]]$chi_stacked
      colnames(m_c) <- paste0(cid, "|", colnames(m_c))
      m_c
    }))
    list(
      rows = rows,
      unit_id = or_fit$unit_id[rows],
      cells = slices,
      chi_stacked = chi_mat
    )
  }

  full_if <- if (materialize_if) if_callback(seq_len(n)) else NULL

  structure(
    list(
      n = n,
      unit_id = or_fit$unit_id,
      cell_ids = cids,
      method = per_cell[[1L]]$method,
      paper_method = per_cell[[1L]]$paper_method,
      cells = per_cell,
      stacked_links = stacked_links,
      if_callback = if_callback,
      influence = full_if,
      chi_stacked = if (materialize_if) full_if$chi_stacked else NULL
    ),
    class = "cdid_static_reconstructed_or_centered_catalog"
  )
}

# Evaluate flow-contracted cell contrast theta_hat^m_c and its original-unit IF
# phi^m_c from centered links and valid role flows varpi^G_c, varpi^C_c.
static_reconstructed_or_cell_contrast <- function(centered_cell,
                                                  rank_tol,
                                                  flow_G = NULL,
                                                  flow_C = NULL) {
  if (!inherits(centered_cell, "cdid_static_reconstructed_or_centered_links")) {
    static_reconstructed_or_abort(
      "`centered_cell` must be a `cdid_static_reconstructed_or_centered_links` object.",
      "cdid_static_reconstructed_or_centering_error"
    )
  }
  rtol <- static_reconstructed_or_validate_rank_tol(rank_tol)
  cid <- centered_cell$cell_id
  t_val <- centered_cell$t
  pm <- centered_cell$paper_method

  fl_G <- static_reconstructed_or_ls_flow(
    edges = centered_cell$G_edges, t = t_val,
    rank_tol = rtol, flow = flow_G, cell_id = paste0(cid, ":G")
  )
  varpi_G <- fl_G$varpi

  if (pm %in% c("ipw", "dr")) {
    fl_C <- static_reconstructed_or_ls_flow(
      edges = centered_cell$C_edges, t = t_val,
      rank_tol = rtol, flow = flow_C, cell_id = paste0(cid, ":C")
    )
    varpi_C <- fl_C$varpi
    theta_hat <- sum(varpi_G * centered_cell$adjusted_G) - sum(varpi_C * centered_cell$adjusted_C)
  } else {
    fl_C <- NULL
    varpi_C <- NULL
    theta_hat <- sum(varpi_G * centered_cell$adjusted_G)
  }

  if_slice <- if (!is.null(centered_cell$influence)) {
    centered_cell$influence
  } else {
    centered_cell$if_callback()
  }
  phi_G <- drop(if_slice$chi_adj_G %*% varpi_G)
  phi_c <- if (pm %in% c("ipw", "dr")) {
    phi_G - drop(if_slice$chi_adj_C %*% varpi_C)
  } else {
    phi_G
  }
  names(phi_c) <- rownames(if_slice$chi_adj_G)

  list(
    cell_id = cid,
    method = centered_cell$method,
    paper_method = pm,
    theta_hat = theta_hat,
    phi = phi_c,
    flow_G = varpi_G,
    flow_C = varpi_C
  )
}
