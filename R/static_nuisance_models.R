# Private role-mask observation-nuisance kernels for the static binary
# entry/exit/reentry missing-outcome extension adapting de Chaisemartin,
# D'Haultfoeuille, Pasquier, Sow and Vazquez-Bare (2026, arXiv:2201.06898v7)
# with Bellego, Benatia and Dortet-Bernadet (2025, doi:10.1016/j.jeconom.2024.105783).
# Generic fit/predict/IF only:
#   * observation logits (pair-observed or date-observed responses) on explicit
#     role masks A with local cell intercepts and slopes shared only inside a
#     declared cell-to-block partition, solved by the existing compact
#     fit_selection_block_provider IRLS (or_pooled_wls_* factors);
#   * raw pair probabilities q = Lambda(u) or the declared date product
#     q_s * q_u with keyed parameter derivatives. NO floor or cap is applied
#     here; the lower probability floor is applied to the final product by the
#     later observation-moment primitive. Conditional date independence is a
#     declared working assumption of the product route, never verified by fitting;
#   * actual-change imputation regressions fitted on A * O rows with local
#     intercepts and declared slope blocks (or_pooled_wls_* QR/partialling),
#     predicted over A;
#   * coefficient influence functions for requested original-unit rows with
#     full-N normalization: all score terms of one original unit (overlapping
#     roles, several cells of one block) are SUMMED before the Gram solve, and
#     each cell's own design is centred at its own Fisher/Gram mean.
# Callers supply retained full-rank designs WITHOUT a constant column (local
# intercepts always exist), unit weights are identically one, and no public
# routing, pooling default, formula handling or gap scaling is chosen here.
#
# Parameter order of a fit: blocks in partition order; inside a block the local
# intercepts in the block's declared cell order, then the block slopes in
# design-column order. Keys are "<block>|alpha|<cell>" and "<block>|beta|<col>".

static_nuisance_abort <- function(message) {
  cdid_abort(message, "cdid_static_nuisance_input_error")
}

static_nuisance_key_ok <- function(x) {
  is.character(x) && !anyNA(x) && all(nzchar(x)) && !any(grepl("|", x, fixed = TRUE)) &&
    !anyDuplicated(x)
}

static_nuisance_validate_ids <- function(unit_id) {
  if (!is.atomic(unit_id) || is.object(unit_id) || !is.null(dim(unit_id)) ||
      !(is.character(unit_id) || is.numeric(unit_id)) || length(unit_id) < 1L ||
      anyNA(unit_id) || anyDuplicated(unit_id))
    static_nuisance_abort("`unit_id` must be a nonempty unique non-NA character/numeric vector (full N).")
  invisible(unit_id)
}

# Requested original-unit IDs -> row positions (order preserved, no duplicates).
static_nuisance_rows <- function(fit, ids) {
  if (!is.atomic(ids) || is.object(ids) || !is.null(dim(ids)) || anyNA(ids) ||
      !(is.character(ids) || is.numeric(ids)) ||
      anyDuplicated(ids) || is.character(ids) != is.character(fit$unit_id))
    static_nuisance_abort("Requested `ids` must be unique non-NA original unit IDs of the fit's ID type.")
  pos <- match(ids, fit$unit_id)
  if (anyNA(pos)) static_nuisance_abort("Requested `ids` contain unknown original unit IDs.")
  pos
}

# Validates the explicit cell list, the cell-to-block partition and per-block
# design dimensions/names. Only structure is checked here; values are read
# solely on the needed rows by the fitting helpers (which keep their own
# finite/binary/rank refusals).
static_nuisance_layout <- function(N, cells, partition, kind) {
  if (!is.list(cells) || length(cells) < 1L || !static_nuisance_key_ok(names(cells)))
    static_nuisance_abort("`cells` must be a nonempty list with unique nonempty keys without '|'.")
  if (!is.list(partition) || length(partition) < 1L || !static_nuisance_key_ok(names(partition)))
    static_nuisance_abort("`partition` must be a nonempty named list of blocks with unique keys without '|'.")
  members <- unlist(partition, use.names = FALSE)
  if (!all(vapply(partition, function(p) is.character(p) && length(p) >= 1L && !anyNA(p), logical(1))) ||
      anyDuplicated(members) || !setequal(members, names(cells)) || length(members) != length(cells))
    static_nuisance_abort("`partition` must assign every cell key to exactly one nonempty block.")
  need <- if (kind == "logit") c("mask", "response") else c("mask", "observed", "dy")
  for (k in names(cells)) {
    cl <- cells[[k]]
    if (!is.list(cl) || !all(need %in% names(cl)))
      static_nuisance_abort(sprintf("Cell '%s' must supply %s and `design`.", k, paste0("`", need, "`", collapse = ", ")))
    if (!is.logical(cl$mask) || !is.null(dim(cl$mask)) || length(cl$mask) != N || anyNA(cl$mask))
      static_nuisance_abort(sprintf("Cell '%s': `mask` must be a non-NA logical vector of length N.", k))
    vec_n <- if (kind == "logit") "response" else c("observed", "dy")
    for (v in vec_n) if (!is.atomic(cl[[v]]) || !is.null(dim(cl[[v]])) || length(cl[[v]]) != N)
      static_nuisance_abort(sprintf("Cell '%s': `%s` must be a vector of length N.", k, v))
    if (kind == "wls") {
      o <- cl$observed[cl$mask]
      if (!is.logical(cl$observed) || anyNA(o))
        static_nuisance_abort(sprintf("Cell '%s': `observed` must be logical and non-NA on the role mask.", k))
      if (!is.numeric(cl$dy)) static_nuisance_abort(sprintf("Cell '%s': `dy` must be numeric.", k))
    }
    if (!is.null(cl$design) && (!is.matrix(cl$design) || !is.numeric(cl$design) ||
                                is.object(cl$design) || nrow(cl$design) != N))
      static_nuisance_abort(sprintf("Cell '%s': `design` must be NULL or a numeric N-row matrix.", k))
  }
  blocks <- vector("list", length(partition))
  names(blocks) <- names(partition)
  off <- 0L
  keys <- character(0)
  for (m in names(partition)) {
    ks <- partition[[m]]
    dims <- vapply(ks, function(k) if (is.null(cells[[k]]$design)) 0L else ncol(cells[[k]]$design), integer(1))
    if (length(unique(dims)) != 1L)
      static_nuisance_abort(sprintf("Block '%s': cell designs have different column counts.", m))
    d_s <- dims[[1L]]
    cn <- if (d_s > 0L) colnames(cells[[ks[1L]]]$design) else character(0)
    if (d_s > 0L && (!static_nuisance_key_ok(cn) ||
        !all(vapply(ks, function(k) identical(colnames(cells[[k]]$design), cn), logical(1)))))
      static_nuisance_abort(sprintf("Block '%s': designs need identical unique nonempty column names without '|'.", m))
    p_m <- length(ks) + d_s
    blocks[[m]] <- list(block = m, cells = ks, d_s = d_s, slope_names = cn,
                        index = off + seq_len(p_m),
                        alpha_index = off + seq_along(ks),
                        beta_index = off + length(ks) + seq_len(d_s))
    keys <- c(keys, paste(m, "alpha", ks, sep = "|"), if (d_s) paste(m, "beta", cn, sep = "|"))
    off <- off + p_m
  }
  cell_block <- stats::setNames(rep(names(partition), lengths(partition)), members)
  cell_local <- stats::setNames(unlist(lapply(partition, seq_along), use.names = FALSE), members)
  list(blocks = blocks, param_keys = keys, cell_block = cell_block, cell_local = cell_local)
}

static_nuisance_chunk <- function(max_chunk_rows, N) {
  if (is.null(max_chunk_rows)) return(N)
  if (!is.numeric(max_chunk_rows) || length(max_chunk_rows) != 1L || !is.finite(max_chunk_rows) ||
      max_chunk_rows < 1 || max_chunk_rows != floor(max_chunk_rows))
    static_nuisance_abort("`max_chunk_rows` must be a positive integer.")
  as.integer(max_chunk_rows)
}

static_nuisance_design_rows <- function(cl, rr, d_s) {
  if (d_s > 0L) unname(cl$design[rr, , drop = FALSE]) else matrix(numeric(0), length(rr), 0L)
}

# ---------------------------------------------------------------------------
# Observation-probability logit fit: cells[[key]] = list(mask, response, design);
# response is the pair- or date-observed indicator (read on mask only; units
# with no recorded outcome are zeros there and stay in N and in every mask).
# Unit weights are one; the repeat-rate channel is off (a_k = 1); mass_U_k = n_b / N
# only feeds the provider's diagnostic observed rate `pi`.
static_nuisance_logit_fit <- function(unit_id, cells, partition, glm_control,
                                      max_chunk_rows = NULL) {
  static_nuisance_validate_ids(unit_id)
  N <- length(unit_id)
  lay <- static_nuisance_layout(N, cells, partition, "logit")
  gc <- validate_ipw_glm_control(glm_control)
  chunk <- static_nuisance_chunk(max_chunk_rows, N)
  fits <- lapply(lay$blocks, function(bl) {
    ks <- bl$cells
    rows <- lapply(ks, function(k) which(cells[[k]]$mask))
    provider <- function(b, consume) {
      cl <- cells[[ks[b]]]
      rb <- rows[[b]]
      n_b <- length(rb)
      if (n_b == 0L) return(invisible(NULL))
      for (start in seq.int(1L, n_b, by = chunk)) {
        rr <- rb[seq.int(start, min(n_b, start + chunk - 1L))]
        consume(list(D = static_nuisance_design_rows(cl, rr, bl$d_s),
                     y = as.double(cl$response[rr]), w = rep.int(1, length(rr)),
                     a_k = 1, mass_U_k = n_b / N))
      }
    }
    spec <- list(N = N, d = 1L + bl$d_s, has_intercept = TRUE, int_col = 1L,
                 slope_cols = 1L + seq_len(bl$d_s), d_s = bl$d_s,
                 branch = if (bl$d_s > 0L) "intercept_slope" else "intercept_only",
                 slope_names = bl$slope_names)
    f <- fit_selection_block_provider(provider, length(ks), spec, gc, N,
                                      spec_id = list(model = "static_nuisance_logit",
                                                     block = bl$block, cells = ks))
    f$n_fit <- sum(lengths(rows))
    f
  })
  coef <- unlist(lapply(fits, function(f) c(f$alpha, f$beta)), use.names = FALSE)
  names(coef) <- lay$param_keys
  structure(list(kind = "logit", unit_id = unit_id, N = N, cells = cells,
                 partition = partition, layout = lay, block_fits = fits,
                 coefficients = coef, max_chunk_rows = chunk),
            class = "cdid_static_nuisance_fit")
}

# ---------------------------------------------------------------------------
# Actual-change imputation WLS fit: cells[[key]] = list(mask = A, observed = O,
# dy, design). The fit rows are A & O; dy and design are read only there
# (prediction designs only on A). Designs are used exactly as supplied
# (no gap scaling is chosen here).
static_nuisance_wls_fit <- function(unit_id, cells, partition, rank_tol,
                                    max_chunk_rows = NULL) {
  static_nuisance_validate_ids(unit_id)
  N <- length(unit_id)
  lay <- static_nuisance_layout(N, cells, partition, "wls")
  chunk <- static_nuisance_chunk(max_chunk_rows, N)
  fits <- lapply(lay$blocks, function(bl) {
    partials <- lapply(bl$cells, function(k) {
      cl <- cells[[k]]
      rb <- which(cl$mask & cl$observed %in% TRUE)
      n_b <- length(rb)
      if (n_b == 0L)  # existing empty-intercept-block rank refusal
        return(or_pooled_wls_partial_block(or_pooled_wls_block_summary(
          matrix(numeric(0), 0L, bl$d_s), numeric(0), numeric(0))))
      chunks <- lapply(seq.int(1L, n_b, by = chunk), function(start) {
        rr <- rb[seq.int(start, min(n_b, start + chunk - 1L))]
        or_pooled_wls_block_summary(static_nuisance_design_rows(cl, rr, bl$d_s),
                                    as.double(cl$dy[rr]), rep.int(1, length(rr)))
      })
      or_pooled_wls_partial_block(if (length(chunks) == 1L) chunks[[1L]] else
        or_pooled_wls_merge_chunks(chunks))
    })
    names(partials) <- bl$cells
    f <- or_pooled_wls_fit(or_pooled_wls_merge(partials), rank_tol = rank_tol)
    f$alpha <- unname(f$intercepts)
    f$beta <- unname(f$coefficients)
    f$m_b <- vapply(f$blocks, function(r) r$weight_sum / N, numeric(1))
    f$z_bar <- if (bl$d_s > 0L) do.call(rbind, lapply(f$blocks, function(r) unname(r$mean_d))) else
      matrix(numeric(0), length(bl$cells), 0L)
    f
  })
  coef <- unlist(lapply(fits, function(f) c(f$alpha, f$beta)), use.names = FALSE)
  names(coef) <- lay$param_keys
  structure(list(kind = "wls", unit_id = unit_id, N = N, cells = cells,
                 partition = partition, layout = lay, block_fits = fits,
                 coefficients = coef, max_chunk_rows = chunk),
            class = "cdid_static_nuisance_fit")
}

static_nuisance_check_fit <- function(fit, kind = NULL) {
  if (!inherits(fit, "cdid_static_nuisance_fit") || (!is.null(kind) && !identical(fit$kind, kind)))
    static_nuisance_abort(sprintf("Expected a static nuisance %s fit.", if (is.null(kind)) "" else kind))
  invisible(fit)
}

static_nuisance_theta <- function(fit, coefficients) {
  if (is.null(coefficients)) return(fit$coefficients)
  if (!is.numeric(coefficients) || length(coefficients) != length(fit$coefficients) ||
      any(!is.finite(coefficients)) ||
      (!is.null(names(coefficients)) && !identical(names(coefficients), fit$layout$param_keys)))
    static_nuisance_abort("`coefficients` must be a finite full parameter vector in fit key order.")
  coefficients
}

# Linear index of one cell on rows `pos` (all inside its mask), with the
# global parameter indices and the local regressor matrix [1, X].
static_nuisance_index <- function(fit, cell, pos, theta) {
  if (!is.character(cell) || length(cell) != 1L || !cell %in% names(fit$cells))
    static_nuisance_abort("Unknown cell key.")
  cl <- fit$cells[[cell]]
  if (!all(cl$mask[pos]))
    static_nuisance_abort(sprintf("Requested rows lie outside the role mask of cell '%s'.", cell))
  bl <- fit$layout$blocks[[fit$layout$cell_block[[cell]]]]
  X <- static_nuisance_design_rows(cl, pos, bl$d_s)
  cp04_finite(X, "static nuisance prediction design")
  a_idx <- bl$alpha_index[fit$layout$cell_local[[cell]]]
  eta <- theta[a_idx] + (if (bl$d_s) drop(X %*% theta[bl$beta_index]) else 0)
  list(eta = rep_len(eta, length(pos)), idx = c(a_idx, bl$beta_index),
       B = cbind(1, X))
}

static_nuisance_pred_out <- function(fit, ids, value, deriv, idx, extra = list()) {
  rownames(deriv) <- as.character(ids)
  colnames(deriv) <- fit$layout$param_keys[idx]
  c(list(ids = ids, value = stats::setNames(value, as.character(ids)), deriv = deriv,
         param_index = idx, param_keys = fit$layout$param_keys[idx]), extra)
}

# Raw observation probability for requested rows. Exactly one of
#   pair  = one fitted pair-observation cell key       -> q = Lambda(u), or
#   dates = c(start-date cell key, end-date cell key)  -> q = q_s * q_u.
# `deriv` is dq/dtheta on the union of touched coordinates (`param_index`).
# No floor/cap: the raw product is returned for the later lower-floor primitive.
static_nuisance_q_predict <- function(fit, ids, pair = NULL, dates = NULL, coefficients = NULL) {
  static_nuisance_check_fit(fit, "logit")
  if (is.null(pair) == is.null(dates))
    static_nuisance_abort("Supply exactly one of `pair` or `dates`.")
  pos <- static_nuisance_rows(fit, ids)
  theta <- static_nuisance_theta(fit, coefficients)
  comp <- function(cell) {
    ix <- static_nuisance_index(fit, cell, pos, theta)
    q <- ipw_clip_logit(ix$eta)$prob  # raw probability; the cap acts on odds only
    cp04_finite(q, "static observation probabilities")
    list(q = q, dq = (q * (1 - q)) * ix$B, idx = ix$idx)
  }
  if (!is.null(pair)) {
    cq <- comp(pair)
    return(static_nuisance_pred_out(fit, ids, cq$q, cq$dq, cq$idx,
                                    list(route = "pair", floored = FALSE)))
  }
  if (!is.character(dates) || length(dates) != 2L || anyNA(dates) || dates[1L] == dates[2L])
    static_nuisance_abort("`dates` must be two distinct declared date-cell keys.")
  cs <- comp(dates[1L])
  cu <- comp(dates[2L])
  idx <- unique(c(cs$idx, cu$idx))
  D <- matrix(0, length(pos), length(idx))
  js <- match(cs$idx, idx)
  ju <- match(cu$idx, idx)
  D[, js] <- D[, js] + cu$q * cs$dq
  D[, ju] <- D[, ju] + cs$q * cu$dq
  static_nuisance_pred_out(fit, ids, cs$q * cu$q, D, idx,
    list(route = "date_product", floored = FALSE,
         q_dates = cbind(start = cs$q, end = cu$q),
         assumption = "conditional date independence: declared working assumption, not tested by the fit"))
}

# Imputed actual change j = alpha_cell + d' beta on requested rows of A, with
# dj/dtheta = [1, d] on the cell's intercept and block slopes.
static_nuisance_j_predict <- function(fit, ids, cell, coefficients = NULL) {
  static_nuisance_check_fit(fit, "wls")
  pos <- static_nuisance_rows(fit, ids)
  ix <- static_nuisance_index(fit, cell, pos, static_nuisance_theta(fit, coefficients))
  static_nuisance_pred_out(fit, ids, ix$eta, ix$B, ix$idx, list(route = "imputation"))
}

# ---------------------------------------------------------------------------
# First-step nuisance influence function on requested original-unit rows (any
# rows, including units with no recorded outcome or outside every mask).
# Returns length(ids) x P with theta_hat - theta* ~ N^{-1} sum_i IF_i (full N).
# Per block, with cell score s_ib (logit: A_ib (y_ib - q_ib); WLS: A_ib O_ib (dy_ib - j_ib))
# and the block's Fisher/Gram mass m_b and mean zbar_b of cell b:
#   xi_beta  = (Gram_slope / N)^{-1} sum_b s_ib (X_ib - zbar_b)
#   xi_alpha_b = s_ib / m_b - zbar_b' xi_beta,
# using the retained compact slope factor (or_wls_gram_solve). Each cell's own
# X_ib is used, so time/cell-varying designs are exact (no common-Z shortcut).
static_nuisance_if <- function(fit, ids) {
  static_nuisance_check_fit(fit)
  pos <- static_nuisance_rows(fit, ids)
  n <- length(pos)
  N <- fit$N
  out <- matrix(0, n, length(fit$coefficients),
                dimnames = list(as.character(ids), fit$layout$param_keys))
  for (m in names(fit$layout$blocks)) {
    bl <- fit$layout$blocks[[m]]
    bf <- fit$block_fits[[m]]
    B_m <- length(bl$cells)
    cs <- matrix(0, n, B_m)
    st <- matrix(0, n, bl$d_s)
    for (b in seq_len(B_m)) {
      cl <- fit$cells[[bl$cells[b]]]
      sel <- which(cl$mask[pos])
      if (!length(sel)) next
      r <- pos[sel]
      if (identical(fit$kind, "logit")) {
        X <- static_nuisance_design_rows(cl, r, bl$d_s)
        eta <- bf$alpha[b] + (if (bl$d_s) drop(X %*% bf$beta) else 0)
        s <- as.double(cl$response[r]) - ipw_clip_logit(eta)$prob
      } else {
        o <- which(cl$observed[r])  # A & O; dy/design never read elsewhere
        sel <- sel[o]
        r <- r[o]
        if (!length(sel)) next
        X <- static_nuisance_design_rows(cl, r, bl$d_s)
        s <- as.double(cl$dy[r]) - (bf$alpha[b] + (if (bl$d_s) drop(X %*% bf$beta) else 0))
      }
      cp04_finite(s, "static nuisance cell scores")
      cs[sel, b] <- s
      if (bl$d_s) st[sel, ] <- st[sel, ] + s * sweep(X, 2L, bf$z_bar[b, ], "-")
    }
    xi_beta <- if (bl$d_s) {
      t(if (identical(fit$kind, "logit")) or_wls_gram_solve(bf$slope_fit, t(st), normalizer = N) else
          or_pooled_wls_gram_solve(bf, t(st), normalizer = N))
    } else matrix(0, n, 0L)
    xi_alpha <- sweep(cs, 2L, bf$m_b, "/") - xi_beta %*% t(bf$z_bar)
    out[, bl$index] <- cbind(xi_alpha, xi_beta)
  }
  out
}
