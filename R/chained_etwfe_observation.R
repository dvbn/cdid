# ==============================================================================
# Private bounded observation-correction provider for linear cohort_time_fe
# ==============================================================================
# etwfe_observation_fit(Y, cohort, cells, X, selection, method, imputation_pooling,
#                       glm_control, rank_tol, chunk_rows, workspace_entries)
#
# Fits ONLY the observation nuisances of the full-frame, observation-corrected
# cohort-by-date linear outcome regression and serves, per supplied cell
# c = (cohort k, s, t), the corrected pair change of every unit of cohort k:
#   ipw         Zc = S Delta / q_safe
#   imputation  Zc = j
#   dr          Zc = j + S (Delta - j) / q_safe
# with S = Y[, s] and Y[, t] both non-NA, Delta = Y[, t] - Y[, s] read ONLY where
# S is TRUE (NA is never multiplied), and Zc defined on the full unit cohort frame
# (units without any outcome stay in N, in their cohort and in every F roster).
# The pair catalog, public calls, REG/IPW/DR treatment geometry, inference and
# integration stay with the caller.
#
# Inputs. Y: N x T double, NA = missing. cohort: length-N integer labels (0 =
# never treated, positive = adoption index). cells: data.frame(cohort, s, t);
# s < t are 1-based column (stored calendar) indices of Y; distinct rows; the
# caller lists only required, supported edges. X: full-frame outcome formula
# design (named "(Intercept)" marks the intercept; no offset; may hold zero
# slope columns); read only when j is fitted and only on active cohort rows.
# selection: NULL for pure imputation, else list(Z, U, rho, universe, model,
# pooling); Z and rho are read only on U.
#
# Observation probability q (weighting arms only):
#   * fitted on the declared roster U of each cell's cohort with the separate
#     selection design Z: extract_selection_design_spec() fixes the local
#     intercept / no-intercept / zero-column (~0: eta = 0, q = 1/2 per pair or
#     per date) branch; partition_selection_cells() gives the canonical pair or
#     date cells and pooling blocks; fit_selection_block_provider() runs the
#     accepted block IRLS (unit weights 1; repeat-rate channel off, a_k = 1)
#     with its rank / separation / boundary / convergence refusals.
#   * pair_logit: q = Lambda(u). independent_period_logit: q = q_s q_t, the
#     declared conditional date-independence working law (U = F, rho = 1 only).
#   * universe F: q_safe = max(q, 1e-6), floor on the final (product) value.
#   * universe A / B (direct pair-logit IPW only): q_safe = pi_k rho max(q_U, 1e-6),
#     pi_k = Pr(U | cohort k) = n_{k,U} / n_k is ONE estimated coordinate per
#     active cohort (score N / n_k I(g = k) (U - pi_k), d Zc / d pi_k = -Zc / pi_k).
#     The previous repeat-rate parameter a_k is not carried.
# Actual-change imputation j (imputation / dr): j_c = alpha_c + (t - s) z' beta_r
# with local cell intercepts (when X has "(Intercept)") and slopes pooled by
# imputation_pooling, fitted by unit-weight LS on observed pairs through the
# accepted or_pooled_support_* / or_pooled_fit_blocks kernels; the alias
# convention keeps exactly the estimable columns, verified on the full cohort
# prediction frame (a dropped column is reported, never silently replaced).
#
# Gradients (local coordinates, global order):
#   q coefficients   -S (Delta - j) 1{q > eps} dq/dtheta / (rho pi_U max(q, eps)^2)
#                    (j = 0 for weighting; rho pi_U = 1 on F; product rule on q_s q_t)
#   pi_k (A / B)     -Zc / pi_k
#   j coefficients   (1 - S / q_safe) dj/db (dr), dj/db (imputation)
# Global parameter order: selection blocks in canonical block order (local
# intercepts in block cell order, then block slopes), then pi_k for active
# cohorts in ascending cohort order (A / B only), then imputation groups in
# canonical order (local intercepts in group cell order, then kept slopes).
# score_rows() returns N-scaled coefficient influence rows (theta_hat - theta ~
# N^{-1} sum_i IF_i); every raw pair/date row of one original unit inside one
# pooled fit is summed BEFORE the Gram / Fisher solve.
#
# Memory. Y, X, Z, U, rho and cohort are borrowed (never copied or modified);
# units are located through their cohort labels. Retained: fitted coefficients,
# compact block / group factors and per-cell coordinate metadata (the per-cohort
# row lists are fitting scratch only). No N x pairs, N x cells or N x parameters
# array is kept; every scan and callback works in row batches of at most
# chunk_rows units. The retained charge is the actual payload of the state and
# diagnostics (cp04_payload convention, computed in R); a pre-fit bound is
# charged first.
#
# Package helpers used unchanged (no other file is modified): cp04_abort,
# cp04_finite, imputation_abort, imputation_pooling_choice,
# validate_ipw_glm_control, validate_selection_mode_args,
# extract_selection_design_spec, partition_selection_cells,
# fit_selection_block_provider, evaluate_selection_block_if,
# or_pooled_outcome_recipe, or_pooled_support_start, or_pooled_support_add,
# or_pooled_fit_blocks, or_pooled_wls_gram_solve, or_wls_gram_solve,
# ipw_clip_logit, ipw_clip_floor, ipw_clip_floor_active, ipw_clip_summary,
# ipw_clip_summary_merge, ipw_owner_safe_sum.
# ==============================================================================

etwfe_obs_count <- function(x) {
  is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x == floor(x) && x >= 1
}

etwfe_obs_rows <- function(rows, N, what) {
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || anyNA(rows) ||
      any(!is.finite(rows)) || any(rows != floor(rows)) || any(rows < 1 | rows > N) ||
      anyDuplicated(rows))
    cp04_abort(paste0("Invalid `rows` for ", what, ": supply unique 1-based original unit positions."),
               "input")
  as.integer(rows)
}

# Row batches are generated lazily: only the batch starts (n / B integers) exist
# up front; each batch index range is built inside its loop iteration.
etwfe_obs_starts <- function(n, B) if (n < 1L) integer(0) else seq.int(1L, n, by = B)
etwfe_obs_ix <- function(a, n, B) seq.int(a, min(n, a + B - 1L))

# Observed pair indicator and Delta on `rows`; outcomes are read only where both
# dates are observed. NaN is not a missing-outcome code.
etwfe_obs_pair <- function(Y, rows, s, t, read = TRUE) {
  ys <- Y[rows, s]
  yt <- Y[rows, t]
  if (any(is.nan(ys)) || any(is.nan(yt)))
    cp04_abort("Outcome matrix holds NaN; only NA marks a missing outcome.", "numeric")
  S <- !is.na(ys) & !is.na(yt)
  out <- list(S = S)
  if (read) {
    dy <- numeric(length(rows))
    o <- which(S)
    if (length(o)) dy[o] <- yt[o] - ys[o]
    cp04_finite(dy, "observed pair differences")
    out$dy <- dy
  }
  out
}

etwfe_obs_date <- function(Y, rows, d) {
  y <- Y[rows, d]
  if (any(is.nan(y)))
    cp04_abort("Outcome matrix holds NaN; only NA marks a missing outcome.", "numeric")
  as.double(!is.na(y))
}

etwfe_obs_sel_D <- function(st, rr) {
  if (st$design$d_s < 1L) return(matrix(0, length(rr), 0L))
  D <- st$Z[rr, st$design$slope_cols, drop = FALSE]
  dimnames(D) <- NULL
  D
}

etwfe_obs_sel_y <- function(st, cell, rr) {
  if (identical(cell$type, "date")) return(etwfe_obs_date(st$Y, rr, cell$calendar_index))
  as.double(etwfe_obs_pair(st$Y, rr, cell$pair$s_calendar_index, cell$pair$t_calendar_index,
                           FALSE)$S)
}

# Raw fitted probability of one selection cell (block m, local b) on rows with
# slope design D. Zero-column blocks give eta = 0, q = 1/2. No floor / cap here.
etwfe_obs_sel_prob <- function(st, m, b, D) {
  bf <- st$sel_blocks[[m]]
  n <- nrow(D)
  eta <- (if (bf$has_intercept) bf$alpha[b] else 0) +
    (if (bf$d_s > 0L) drop(D %*% bf$beta) else 0)
  ipw_clip_logit(rep_len(as.double(eta), n))$prob
}

# Raw pair probability of data cell i on observed-pair rows `ro` (all in U),
# with dq/dtheta on the cell's selection coordinates when requested.
etwfe_obs_pair_q <- function(st, i, ro, deriv) {
  cs <- st$cell_sel[[i]]
  D <- etwfe_obs_sel_D(st, ro)
  qs <- lapply(cs$comps, function(cp) etwfe_obs_sel_prob(st, cp$block, cp$local, D))
  q <- if (length(qs) == 1L) qs[[1L]] else qs[[1L]] * qs[[2L]]
  cp04_finite(q, "fitted observation probabilities on observed pairs", TRUE)
  dq <- NULL
  if (deriv) {
    dq <- matrix(0, length(ro), length(cs$sel_coords))
    if (length(cs$sel_coords)) {
      Breg <- if (st$design$has_intercept) cbind(1, D) else D
      for (a in seq_along(cs$comps)) {
        pos <- cs$comp_pos[[a]]
        if (!length(pos)) next
        other <- if (length(qs) == 1L) 1 else qs[[3L - a]]
        qa <- qs[[a]]
        dq[, pos] <- dq[, pos] + (other * qa * (1 - qa)) * Breg
      }
    }
  }
  list(q = q, dq = dq)
}

# Imputed actual change j = alpha_c + L z_keep' beta_r on rows of the cell's
# cohort (full cohort frame), with dj/db on the cell's j coordinates.
etwfe_obs_j <- function(st, i, rr, deriv) {
  g <- st$imp_groups[[st$cell_group[i]]]
  loc <- st$cell_local[i]
  n <- length(rr)
  Dz <- if (g$p > 0L) st$cell_L[i] * unname(st$X[rr, g$kept_X, drop = FALSE]) else matrix(0, n, 0L)
  j <- (if (st$icpt) g$alpha[loc] else 0) + (if (g$p > 0L) drop(Dz %*% g$beta) else 0)
  j <- rep_len(as.double(j), n)
  cp04_finite(j, "imputed pair changes")
  list(j = j, B = if (deriv) (if (st$icpt) cbind(1, Dz) else Dz) else NULL)
}

etwfe_obs_cell_chunk <- function(st, i, rr, gradient) {
  n <- length(rr)
  lay <- st$cell_layout[[i]]
  pr <- etwfe_obs_pair(st$Y, rr, st$cell_s[i], st$cell_t[i], TRUE)
  obs <- which(pr$S)
  G <- if (gradient) matrix(0, n, lay$n_coords) else NULL
  j <- numeric(n)
  jj <- NULL
  if (st$need_j) {
    jj <- etwfe_obs_j(st, i, rr, gradient)
    j <- jj$j
  }
  value <- j
  den <- NULL
  if (st$need_q && length(obs)) {
    ro <- rr[obs]
    if (!all(st$U[ro])) cp04_abort("U must contain all A structurally.", "support")
    qq <- etwfe_obs_pair_q(st, i, ro, gradient)
    qf <- ipw_clip_floor(qq$q)
    scl <- if (identical(st$universe, "F")) rep.int(1, length(obs)) else
      st$pi_val[st$cell_kk[i]] * st$rho[ro]
    den <- scl * qf
    e <- pr$dy[obs] - j[obs]
    value[obs] <- j[obs] + e / den
    if (gradient && length(lay$sel_pos))
      G[obs, lay$sel_pos] <- (-e * ipw_clip_floor_active(qq$q) / (scl * qf^2)) * qq$dq
  }
  cp04_finite(value, "reconstructed pair changes")
  if (gradient) {
    if (length(lay$pi_pos)) G[, lay$pi_pos] <- -value / st$pi_val[st$cell_kk[i]]
    if (st$need_j && length(lay$j_pos)) {
      fac <- rep.int(1, n)
      if (!is.null(den)) fac[obs] <- 1 - 1 / den
      G[, lay$j_pos] <- fac * jj$B
    }
    cp04_finite(G, "reconstructed pair-change gradients")
  }
  list(value = value, gradient = G)
}

etwfe_obs_read <- function(st, cell, rows, gradient = FALSE) {
  if (!is.numeric(cell) || is.object(cell) || !is.null(dim(cell)) || length(cell) != 1L ||
      is.na(cell) || !is.finite(cell) || cell != floor(cell) || cell < 1 || cell > st$C)
    cp04_abort("`cell` must be a 1-based index into the supplied cells.", "input")
  if (!identical(gradient, TRUE) && !identical(gradient, FALSE))
    cp04_abort("`gradient` must be TRUE or FALSE.", "input")
  i <- as.integer(cell)
  rows <- etwfe_obs_rows(rows, st$N, "read")
  if (length(rows) && any(st$cohort[rows] != st$cell_k_val[i]))
    cp04_abort("`read` rows must be units of the cell's cohort.", "input")
  n <- length(rows)
  co <- st$coords[[i]]
  value <- numeric(n)
  G <- if (gradient) matrix(0, n, length(co)) else NULL
  for (a in etwfe_obs_starts(n, st$B)) {
    ix <- etwfe_obs_ix(a, n, st$B)
    r <- etwfe_obs_cell_chunk(st, i, rows[ix], gradient)
    value[ix] <- r$value
    if (gradient) G[ix, ] <- r$gradient
  }
  list(value = value, coordinates = co, gradient = G)
}

etwfe_obs_score_chunk <- function(st, rc) {
  nb <- length(rc)
  N <- st$N
  res <- matrix(0, nb, st$P)
  kc <- st$cohort[rc]
  if (st$need_q) {
    Uc <- st$U[rc]
    inU <- which(Uc)
    d_s <- st$design$d_s
    Zs <- matrix(0, nb, d_s)
    if (d_s > 0L && length(inU)) Zs[inU, ] <- etwfe_obs_sel_D(st, rc[inU])
    for (m in seq_along(st$sel_blocks)) {
      bf <- st$sel_blocks[[m]]
      if (!length(bf$index)) next
      cs <- matrix(0, nb, bf$B_m)
      for (b in seq_len(bf$B_m)) {
        cell <- st$sel_cells[[bf$cell_ids[b]]]
        ii <- which(kc == st$coh_vals[cell$cohort] & Uc)
        if (!length(ii)) next
        y <- etwfe_obs_sel_y(st, cell, rc[ii])
        cs[ii, b] <- y - etwfe_obs_sel_prob(st, m, b, Zs[ii, , drop = FALSE])
      }
      xi <- evaluate_selection_block_if(bf, Zs, cs, N)$xi_theta
      res[, bf$index] <- matrix(xi, nb, length(bf$index))
    }
    if (length(st$pi_index)) for (kk in seq_len(st$K)) {
      ii <- which(kc == st$coh_vals[kk])
      if (length(ii))
        res[ii, st$pi_index[kk]] <- (N / st$n_k[kk]) * (as.double(Uc[ii]) - st$pi_val[kk])
    }
  }
  if (st$need_j) for (g in seq_along(st$imp_groups)) {
    gr <- st$imp_groups[[g]]
    B_r <- length(gr$cells)
    p <- gr$p
    h_a <- matrix(0, nb, B_r)
    h_b <- matrix(0, nb, p)
    for (loc in seq_len(B_r)) {
      i <- gr$cells[loc]
      ii <- which(kc == st$cell_k_val[i])
      if (!length(ii)) next
      pr <- etwfe_obs_pair(st$Y, rc[ii], st$cell_s[i], st$cell_t[i], TRUE)
      o <- which(pr$S)
      if (!length(o)) next
      ro <- rc[ii[o]]
      r <- pr$dy[o] - etwfe_obs_j(st, i, ro, FALSE)$j
      h_a[ii[o], loc] <- h_a[ii[o], loc] + r
      if (p > 0L) {
        d <- st$cell_L[i] * unname(st$X[ro, gr$kept_X, drop = FALSE])
        if (st$icpt) d <- sweep(d, 2L, gr$mean_d[loc, ], "-")
        h_b[ii[o], ] <- h_b[ii[o], ] + r * d
      }
    }
    xi_b <- if (p > 0L) {
      sol <- if (st$icpt) or_pooled_wls_gram_solve(gr$fit, t(h_b), normalizer = N) else
        or_wls_gram_solve(gr$fit, t(h_b), normalizer = N)
      t(matrix(sol, p, nb))
    } else matrix(0, nb, 0L)
    xi_a <- if (st$icpt) {
      sweep(h_a, 2L, gr$mass, "/") - (if (p > 0L) xi_b %*% t(gr$mean_d) else 0)
    } else matrix(0, nb, 0L)
    res[, gr$index] <- cbind(xi_a, xi_b)
  }
  cp04_finite(res, "nuisance coefficient influence rows")
  res
}

etwfe_obs_score_rows <- function(st, rows) {
  rows <- etwfe_obs_rows(rows, st$N, "score_rows")
  n <- length(rows)
  out <- matrix(0, n, st$P)
  if (!n || !st$P) return(out)
  for (a in etwfe_obs_starts(n, st$B)) {
    ix <- etwfe_obs_ix(a, n, st$B)
    out[ix, ] <- etwfe_obs_score_chunk(st, rows[ix])
  }
  out
}

# The returned callbacks close over ONE compact state list only.
etwfe_obs_callbacks <- function(st) {
  force(st)
  list(
    read = function(cell, rows, gradient = FALSE) etwfe_obs_read(st, cell, rows, gradient),
    score_rows = function(rows) etwfe_obs_score_rows(st, rows)
  )
}

# Logical payload count with the cp04_payload convention (length + 8 per node,
# recursive lists, ceil(bytes / 2) per string, attributes once); environments,
# closures and NULL charge nothing. Plain R, used on compact objects only.
etwfe_obs_payload <- function(x) {
  if (is.null(x) || is.environment(x) || is.function(x)) return(0)
  n <- length(x) + 8
  if (is.character(x)) n <- n + sum(ceiling(nchar(x, type = "bytes", allowNA = TRUE, keepNA = FALSE) / 2))
  if (is.list(x)) n <- n + sum(vapply(x, etwfe_obs_payload, numeric(1)))
  a <- attributes(x)
  if (length(a)) n <- n + length(a) + 8 + sum(vapply(a, etwfe_obs_payload, numeric(1)))
  n
}

# Pre-fit charge (entries; not RSS). Retained: coefficients and parameter names,
# selection blocks (local intercepts, slope factor, Fisher means), imputation
# group fits (pooled factor and per-cell records), per-cell coordinate metadata
# and diagnostics. Scratch: one row batch of the widest scan / callback plus
# block QR partials. The post-fit retained charge is the actual payload.
etwfe_obs_budget <- function(N, TT, B, K, C, P, sel_blocks_dims, imp_groups_dims, d_sel,
                             p_X, coord_max) {
  sel_ret <- sum(vapply(sel_blocks_dims, function(x) x[1L] * (x[2L] + 16) + 2 * (x[2L] + 2)^2 + 512,
                        numeric(1)))
  imp_ret <- sum(vapply(imp_groups_dims, function(x) x[1L] * (3 * x[2L] + 160) + 3 * (x[2L] + 2)^2 + 768,
                        numeric(1)))
  Bm_max <- max(c(1, vapply(sel_blocks_dims, `[`, numeric(1), 1L)))
  Br_max <- max(c(1, vapply(imp_groups_dims, `[`, numeric(1), 1L)))
  retained <- 40 * P + sel_ret + imp_ret + C * (6 * coord_max + 640) + 64 * K + 8192
  scratch <- B * (P + Bm_max + Br_max + d_sel + p_X + TT + coord_max + 16) + N +
    (Bm_max + 1) * (d_sel + 2)^2 + (Br_max + 1) * (p_X + 2)^2 + 1024
  c(retained = retained, scratch = scratch)
}

etwfe_observation_fit <- function(Y, cohort, cells, X, selection, method, imputation_pooling,
                                  glm_control, rank_tol, chunk_rows, workspace_entries, report = NULL) {
  # ---- 0. Arguments, frame and cells -----------------------------------------
  if (!is.character(method) || is.object(method) || !is.null(dim(method)) ||
      length(method) != 1L || is.na(method) || !method %in% c("ipw", "imputation", "dr"))
    cp04_abort('`method` must be "ipw", "imputation" or "dr" (observation correction).',
               "unsupported")
  need_q <- method %in% c("ipw", "dr")
  need_j <- method %in% c("imputation", "dr")
  if (!is.matrix(Y) || is.object(Y) || !is.numeric(Y) || nrow(Y) < 1L || ncol(Y) < 2L)
    cp04_abort("Supply `Y` as a numeric N x T outcome matrix (T >= 2) with NA for missing outcomes.",
               "input")
  N <- nrow(Y)
  TT <- ncol(Y)
  if (!is.numeric(cohort) || is.object(cohort) || !is.null(dim(cohort)) || length(cohort) != N ||
      anyNA(cohort) || any(!is.finite(cohort)) || any(cohort != floor(cohort)) || any(cohort < 0))
    cp04_abort("`cohort` must be a length-N vector of nonnegative integer labels (0 = never treated).",
               "input")
  if (!etwfe_obs_count(chunk_rows)) cp04_abort("`chunk_rows` must be a positive integer.", "input")
  B <- as.integer(min(chunk_rows, N))
  if (!is.numeric(workspace_entries) || is.object(workspace_entries) ||
      !is.null(dim(workspace_entries)) || length(workspace_entries) != 1L ||
      is.na(workspace_entries) || !is.finite(workspace_entries) || workspace_entries <= 0)
    cp04_abort("`workspace_entries` must be a finite positive numeric bound.", "input")
  if (!is.data.frame(cells) || nrow(cells) < 1L || !all(c("cohort", "s", "t") %in% names(cells)))
    cp04_abort("`cells` must be a nonempty data.frame with columns cohort, s, t.", "input")
  int_col <- function(x) is.numeric(x) && !is.object(x) && !anyNA(x) && all(is.finite(x)) &&
    all(x == floor(x))
  if (!int_col(cells$cohort) || !int_col(cells$s) || !int_col(cells$t) || any(cells$cohort < 0) ||
      any(cells$s < 1 | cells$s > TT) || any(cells$t < 1 | cells$t > TT) || any(cells$s >= cells$t))
    cp04_abort("Cells need integer cohort >= 0 and Y column indices 1 <= s < t <= T.", "input")
  C <- nrow(cells)
  cell_k_val <- as.double(cells$cohort)
  cell_s <- as.integer(cells$s)
  cell_t <- as.integer(cells$t)
  cell_L <- cell_t - cell_s
  labels <- paste0("k", cell_k_val, "_s", cell_s, "_t", cell_t)
  if (anyDuplicated(labels)) cp04_abort("`cells` must hold distinct (cohort, s, t) rows.", "input")
  coh_vals <- sort(unique(cell_k_val))
  K <- length(coh_vals)
  cell_kk <- match(cell_k_val, coh_vals)
  coh_rows <- unname(split(seq_len(N), factor(match(cohort, coh_vals), levels = seq_len(K))))
  n_k <- lengths(coh_rows)
  if (any(n_k < 1L)) cp04_abort("A cell cohort has no frame units.", "support")
  # Observed pair support of every listed cell (outcome-free scan, refused before any fit).
  cell_npair <- numeric(C)
  for (i in seq_len(C)) {
    rk <- coh_rows[[cell_kk[i]]]
    for (a in etwfe_obs_starts(length(rk), B)) {
      ix <- etwfe_obs_ix(a, length(rk), B)
      cell_npair[i] <- cell_npair[i] + sum(etwfe_obs_pair(Y, rk[ix], cell_s[i], cell_t[i], FALSE)$S)
    }
    if (cell_npair[i] < 1) {
      msg <- paste0("Actual-change cell `", labels[i], "` has zero observed-pair support.")
      if (need_j) imputation_abort(msg, "support") else cp04_abort(msg, "support")
    }
  }

  # ---- 1. Selection arguments, design spec and canonical catalog ------------
  universe <- "F"
  U <- rho <- Zsel <- NULL
  design <- NULL
  scat <- NULL
  if (need_q) {
    if (!is.list(selection) || is.object(selection) || is.null(names(selection)) ||
        !all(c("Z", "U", "rho", "universe", "model", "pooling") %in% names(selection)))
      cp04_abort("Weighted correction needs `selection` = list(Z, U, rho, universe, model, pooling).",
                 "input")
    U <- selection$U
    rho <- selection$rho
    Zsel <- selection$Z
    universe <- selection$universe
    if (!is.logical(U) || is.object(U) || !is.null(dim(U)) || length(U) != N || anyNA(U) ||
        !any(U) || !is.numeric(rho) || is.object(rho) || !is.null(dim(rho)) || length(rho) != N)
      cp04_abort("Invalid fitting roster U or known rho.")
    if (!is.character(universe) || is.object(universe) || length(universe) != 1L ||
        is.na(universe) || !universe %in% c("F", "A", "B"))
      cp04_abort('Selection `universe` must be "F", "A" or "B".', "unsupported")
    if (identical(universe, "F") && !all(U))
      cp04_abort("Universe F fits on the full frame: U must be TRUE and rho = 1 on every unit.",
                 "unsupported")
    for (a in etwfe_obs_starts(N, B)) {
      ix <- etwfe_obs_ix(a, N, B)
      ru <- ix[U[ix]]
      if (!length(ru)) next
      cp04_finite(rho[ru], "known normalized rho on U", TRUE)
      if (identical(universe, "F") && any(rho[ru] != 1))
        cp04_abort("Universe F fits on the full frame: U must be TRUE and rho = 1 on every unit.",
                   "unsupported")
    }
    validate_selection_mode_args(selection$model, selection$pooling, universe, U, rho)
    if (!identical(universe, "F") && !identical(method, "ipw"))
      cp04_abort("Fitting rosters A/B support only the direct pair-logit IPW observation correction.",
                 "unsupported")
    if (!is.matrix(Zsel) || nrow(Zsel) != N)
      cp04_abort("Supply a numeric selection design matrix Z with one row per frame unit.", "input")
    design <- extract_selection_design_spec(Zsel, U, chunk_rows = B)
    if (!identical(universe, "F")) for (a in etwfe_obs_starts(N, B)) {
      ix <- etwfe_obs_ix(a, N, B)
      nobs <- rowSums(!is.na(Y[ix, , drop = FALSE]))
      if (any(nobs >= 2L & !U[ix])) cp04_abort("U must contain all A structurally.", "support")
    }
    pair_specs <- lapply(seq_len(C), function(i) list(cohort = cell_kk[i], pair = list(
      s = cell_s[i], t = cell_t[i], s_calendar_index = cell_s[i], t_calendar_index = cell_t[i])))
    scat <- partition_selection_cells(pair_specs, coh_vals, all_times = seq_len(TT),
                                      selection_model = selection$model,
                                      selection_pooling = selection$pooling)
    rm(pair_specs)
  } else if (!is.null(selection)) {
    imputation_abort(paste("Pure imputation fits no observation-probability model;",
                           "`selection` must be NULL."), "input")
  }

  # ---- 2. Imputation arguments, recipe and canonical groups -----------------
  recipe <- NULL
  groups <- list()
  group_labels <- character(0)
  if (need_j) {
    imputation_pooling <- imputation_pooling_choice(imputation_pooling)
    if (!is.numeric(rank_tol) || is.object(rank_tol) || !is.null(dim(rank_tol)) ||
        length(rank_tol) != 1L || is.na(rank_tol) || !is.finite(rank_tol) || rank_tol <= 0)
      imputation_abort("`rank_tol` must be a finite positive scalar.", "input")
    if (!is.matrix(X) || is.object(X) || !is.numeric(X) || nrow(X) != N)
      imputation_abort("Supply the full-frame outcome design `X` as a numeric N-row matrix.", "design")
    cn <- colnames(X)
    if (is.null(cn)) cn <- character(0)
    if (length(cn) != ncol(X))
      imputation_abort("Full-F outcome/treatment design `X` requires distinct nonempty column names.",
                       "design")
    recipe <- or_pooled_outcome_recipe(cn, "(Intercept)" %in% cn, setdiff(cn, "(Intercept)"))
    for (kk in seq_len(K)) for (a in etwfe_obs_starts(n_k[kk], B))
      cp04_finite(X[coh_rows[[kk]][etwfe_obs_ix(a, n_k[kk], B)], , drop = FALSE], "full-F imputation design")
    ord <- order(cell_k_val, cell_L, cell_s, cell_t, method = "radix")
    if (identical(imputation_pooling, "separate")) {
      groups <- lapply(ord, function(i) i)
      group_labels <- paste0("separate_", labels[ord])
    } else if (identical(imputation_pooling, "within_cohort")) {
      groups <- lapply(coh_vals, function(k) ord[cell_k_val[ord] == k])
      group_labels <- paste0("within_cohort_k", coh_vals)
    } else {
      groups <- list(ord)
      group_labels <- "across_cohorts_all"
    }
  } else {
    imputation_pooling <- NA_character_
  }

  # ---- 3. Preflight workspace charge (upper bounds before any fit) ---------
  d_sel <- if (need_q) design$d else 0L
  d_s <- if (need_q) design$d_s else 0L
  has_int <- need_q && design$has_intercept
  p_X <- if (need_j) recipe$p else 0L
  icpt <- need_j && recipe$intercept
  sel_dims <- if (need_q) lapply(scat$block_specs, function(bs) c(bs$B_m, d_s)) else list()
  imp_dims <- lapply(groups, function(g) c(length(g), p_X))
  P_ub <- (if (need_q) sum(vapply(scat$block_specs, function(bs)
    (if (has_int) bs$B_m else 0) + d_s, numeric(1))) else 0) +
    (if (need_q && !identical(universe, "F")) K else 0) +
    sum(vapply(groups, function(g) (if (icpt) length(g) else 0) + p_X, numeric(1)))
  coord_ub <- 2 * (1 + d_s) + 1 + 1 + p_X
  pre <- etwfe_obs_budget(N, TT, B, K, C, P_ub, sel_dims, imp_dims, d_sel, p_X, coord_ub)
  ipw_owner_safe_sum(pre, workspace_entries, "ETWFE observation provider workspace")

  # ---- 4. Selection fits -----------------------------------------------------
  off <- 0L
  sel_blocks <- list()
  nU <- numeric(K)
  pi_val <- numeric(0)
  pi_index <- integer(0)
  par_names <- character(0)
  coef_all <- numeric(0)
  if (need_q) {
    gc <- validate_ipw_glm_control(glm_control)
    nU <- vapply(coh_rows, function(r) as.double(sum(U[r])), numeric(1))
    stq <- list(Y = Y, Z = Zsel, design = design)
    sel_cell_label <- function(cell) {
      if (identical(cell$type, "date")) paste0("k", coh_vals[cell$cohort], "_d", cell$calendar_index) else
        paste0("k", coh_vals[cell$cohort], "_s", cell$pair$s_calendar_index, "_t",
               cell$pair$t_calendar_index)
    }
    make_provider <- function(cids) {
      force(cids)
      function(b, consume) {
        cell <- scat$cells[[cids[b]]]
        kk <- cell$cohort
        rk <- coh_rows[[kk]]
        mass <- nU[kk] / N
        for (a in etwfe_obs_starts(length(rk), B)) {
          rr <- rk[etwfe_obs_ix(a, length(rk), B)]
          rr <- rr[U[rr]]
          if (!length(rr)) next
          consume(list(D = etwfe_obs_sel_D(stq, rr), y = etwfe_obs_sel_y(stq, cell, rr),
                       w = rep.int(1, length(rr)), a_k = 1, mass_U_k = mass))
        }
        invisible(NULL)
      }
    }
    for (m in seq_along(scat$block_specs)) {
      bs <- scat$block_specs[[m]]
      cids <- bs$cell_ids
      cl <- vapply(scat$cells[cids], sel_cell_label, character(1))
      # Labels carry cohort VALUES (the catalog's own labels carry cohort indices).
      blabel <- switch(scat$selection_pooling,
                       separate = paste0("separate_", cl),
                       within_cohort = paste0("within_cohort_k", coh_vals[bs$cohorts]),
                       "across_cohorts_all")
      bf <- fit_selection_block_provider(make_provider(cids), length(cids), design, gc, N,
        spec_id = list(model = "etwfe_observation", block = blabel, cells = cl))
      B_m <- length(cids)
      p_m <- length(bf$coefficients)
      bf$label <- blabel
      bf$cell_ids <- as.integer(cids)
      bf$alpha_index <- if (bf$has_intercept) off + seq_len(B_m) else integer(0)
      bf$beta_index <- off + (if (bf$has_intercept) B_m else 0L) + seq_len(bf$d_s)
      bf$index <- off + seq_len(p_m)
      off <- off + p_m
      if (p_m) {
        par_names <- c(par_names,
                       if (bf$has_intercept) paste0("q|", blabel, "|alpha|", cl),
                       if (bf$d_s > 0L) paste0("q|", blabel, "|beta|", design$slope_names))
        coef_all <- c(coef_all, bf$coefficients)
      }
      sel_blocks[[m]] <- bf
    }
    if (!identical(universe, "F")) {
      pi_val <- nU / n_k
      cp04_finite(pi_val, "roster probability Pr(U | cohort)", TRUE)
      pi_index <- off + seq_len(K)
      off <- off + K
      par_names <- c(par_names, paste0("piU|k", coh_vals))
      coef_all <- c(coef_all, pi_val)
    }
  }

  # ---- 5. Imputation fits (or_pooled_fit_blocks, one call per group) -------
  imp_groups <- list()
  cell_group <- cell_local <- integer(C)
  if (need_j) {
    slope_X <- recipe$slope_index
    colmax_k <- lapply(seq_len(K), function(kk) {
      cm <- numeric(p_X)
      if (p_X > 0L) for (a in etwfe_obs_starts(n_k[kk], B))
        cm <- pmax(cm, apply(abs(X[coh_rows[[kk]][etwfe_obs_ix(a, n_k[kk], B)], slope_X, drop = FALSE]),
                             2L, max))
      cm
    })
    cell_support <- function(i) {
      acc <- or_pooled_support_start(recipe, cell_L[i])
      rk <- coh_rows[[cell_kk[i]]]
      for (a in etwfe_obs_starts(length(rk), B)) {
        rr <- rk[etwfe_obs_ix(a, length(rk), B)]
        pr <- etwfe_obs_pair(Y, rr, cell_s[i], cell_t[i], TRUE)
        o <- which(pr$S)
        if (!length(o)) next
        acc <- or_pooled_support_add(acc, recipe, X[rr[o], , drop = FALSE], pr$dy[o],
                                     rep.int(1, length(o)), "imputation")
      }
      if (acc$n < 1L || !is.finite(acc$weight_mass) || acc$weight_mass <= 0)
        imputation_abort(paste0("Actual-change imputation cell `", labels[i],
                                "` has zero observed-pair support."), "support")
      acc
    }
    for (g in seq_along(groups)) {
      ids <- groups[[g]]
      B_r <- length(ids)
      supports <- vector("list", 2L * B_r)
      blocks_r <- links_r <- vector("list", B_r)
      for (loc in seq_len(B_r)) {
        i <- ids[loc]
        acc <- cell_support(i)
        supports[[loc]] <- list(summary = acc$summary, scale = 1.0, colmax = acc$colmax)
        supports[[B_r + loc]] <- list(summary = NULL, scale = 1.0, colmax = colmax_k[[cell_kk[i]]])
        blocks_r[[loc]] <- list(name = labels[i], support_ids = as.integer(loc),
                                L = as.integer(cell_L[i]))
        links_r[[loc]] <- list(regression = as.integer(loc), support_id = as.integer(B_r + loc))
      }
      rm(acc)
      visit <- function(sid, visit_fn) {
        tr <- sid <= B_r
        i <- ids[if (tr) sid else sid - B_r]
        rk <- coh_rows[[cell_kk[i]]]
        for (a in etwfe_obs_starts(length(rk), B)) {
          rr <- rk[etwfe_obs_ix(a, length(rk), B)]
          if (tr) {
            o <- which(etwfe_obs_pair(Y, rr, cell_s[i], cell_t[i], FALSE)$S)
            if (length(o)) visit_fn(X[rr[o], , drop = FALSE])
          } else {
            visit_fn(X[rr, , drop = FALSE])
          }
        }
        invisible(NULL)
      }
      fit_r <- tryCatch(or_pooled_fit_blocks(
        recipe = recipe, supports = supports, blocks = blocks_r, links = links_r,
        rank_tol = rank_tol, N = N, visit_rows = visit, route = "imputation"
      ), cdid_or_wls_rank_error = function(e) {
        e$message <- paste0(conditionMessage(e),
          " This computation failure does not establish population ATT nonidentification.")
        e$is_identifiability_claim <- FALSE
        stop(e)
      })
      rm(supports)
      keep <- as.integer(fit_r$keep)
      p_g <- length(keep)
      regs <- fit_r$regressions
      beta <- as.double(fit_r$beta)
      alpha <- if (icpt) vapply(regs, function(r) as.double(r$alpha), numeric(1)) else numeric(0)
      cp04_finite(c(alpha, beta), "imputation coefficients")
      mass <- if (icpt) vapply(regs, function(r) as.double(r$mass), numeric(1)) else numeric(0)
      mean_d <- if (icpt && p_g > 0L) do.call(rbind, lapply(regs, function(r) as.double(r$mean_d))) else
        matrix(0, B_r, p_g)
      if (icpt) cp04_finite(mass, "imputation cell mass", TRUE)
      kept_names <- recipe$slope_columns[keep]
      alpha_index <- if (icpt) off + seq_len(B_r) else integer(0)
      beta_index <- off + (if (icpt) B_r else 0L) + seq_len(p_g)
      imp_groups[[g]] <- list(
        label = group_labels[g], cells = as.integer(ids), keep = keep, kept_X = slope_X[keep],
        kept_columns = kept_names, aliased_columns = fit_r$basis$aliased_columns, p = p_g,
        alpha = alpha, beta = beta, mass = mass, mean_d = mean_d, fit = fit_r$fit,
        rss = fit_r$rss, alpha_index = alpha_index, beta_index = beta_index,
        index = c(alpha_index, beta_index))
      off <- off + length(alpha_index) + p_g
      cell_group[ids] <- g
      cell_local[ids] <- seq_len(B_r)
      par_names <- c(par_names,
                     if (icpt) paste0("j|", group_labels[g], "|alpha|", labels[ids]),
                     if (p_g) paste0("j|", group_labels[g], "|beta|", kept_names))
      coef_all <- c(coef_all, alpha, beta)
      rm(fit_r, regs)
    }
  }
  P <- off
  if (length(par_names) != P || length(coef_all) != P)
    cp04_abort("Internal parameter layout mismatch.", "integrity")
  names(coef_all) <- par_names

  # ---- 6. Per-cell coordinates (sorted global indices, fixed per cell) ------
  cell_sel <- vector("list", C)
  coords <- vector("list", C)
  cell_layout <- vector("list", C)
  for (i in seq_len(C)) {
    sel_coords <- integer(0)
    if (need_q) {
      pc <- scat$pair_cohorts[[scat$request_to_pc[i]]]
      cids <- if (identical(scat$selection_model, "independent_period_logit"))
        c(pc$s_cell_id, pc$t_cell_id) else pc$pair_cell_id
      comps <- lapply(cids, function(cid) {
        cell <- scat$cells[[cid]]
        bf <- sel_blocks[[cell$block_id]]
        list(cid = cid, block = cell$block_id, local = cell$local_b,
             idx = c(if (bf$has_intercept) bf$alpha_index[cell$local_b], bf$beta_index))
      })
      sel_coords <- sort(unique(as.integer(unlist(lapply(comps, `[[`, "idx")))))
      comp_pos <- lapply(comps, function(cp) match(cp$idx, sel_coords))
      cell_sel[[i]] <- list(comps = comps, sel_coords = sel_coords, comp_pos = comp_pos)
    }
    pi_c <- if (length(pi_index)) pi_index[cell_kk[i]] else integer(0)
    j_c <- integer(0)
    if (need_j) {
      gr <- imp_groups[[cell_group[i]]]
      j_c <- c(if (icpt) gr$alpha_index[cell_local[i]], gr$beta_index)
    }
    co <- as.integer(c(sel_coords, pi_c, j_c))
    if (is.unsorted(co, strictly = TRUE)) cp04_abort("Internal coordinate order mismatch.", "integrity")
    coords[[i]] <- co
    ns <- length(sel_coords)
    cell_layout[[i]] <- list(n_coords = length(co), sel_pos = seq_len(ns),
                             pi_pos = if (length(pi_c)) ns + 1L else integer(0),
                             j_pos = ns + length(pi_c) + seq_along(j_c))
  }

  st <- list(
    Y = Y, cohort = cohort, X = if (need_j) X else NULL, Z = Zsel, U = U, rho = rho,
    N = N, TT = TT, B = B, C = C, K = K, P = P, method = method,
    need_q = need_q, need_j = need_j, universe = universe, design = design, icpt = icpt,
    cell_k_val = cell_k_val, cell_kk = cell_kk, cell_s = cell_s, cell_t = cell_t, cell_L = cell_L,
    coh_vals = coh_vals, n_k = n_k, nU = nU, pi_val = pi_val, pi_index = pi_index,
    sel_blocks = sel_blocks, sel_cells = if (need_q) scat$cells else NULL,
    cell_sel = cell_sel, imp_groups = imp_groups, cell_group = cell_group,
    cell_local = cell_local, coords = coords, cell_layout = cell_layout)

  # ---- 7. Compact diagnostics: per-cell support and active floor counts ----
  cell_diag <- vector("list", C)
  for (i in seq_len(C)) {
    clip <- NULL
    if (need_q) {
      rk <- coh_rows[[cell_kk[i]]]
      for (a in etwfe_obs_starts(length(rk), B)) {
        rr <- rk[etwfe_obs_ix(a, length(rk), B)]
        o <- which(etwfe_obs_pair(Y, rr, cell_s[i], cell_t[i], FALSE)$S)
        if (!length(o)) next
        ro <- rr[o]
        if (!all(U[ro])) cp04_abort("U must contain all A structurally.", "support")
        q <- etwfe_obs_pair_q(st, i, ro, FALSE)$q
        den <- (if (identical(universe, "F")) 1 else pi_val[cell_kk[i]] * rho[ro]) * ipw_clip_floor(q)
        clip <- ipw_clip_summary_merge(clip, ipw_clip_summary(q = q, weights = 1 / den))
      }
    }
    cell_diag[[i]] <- list(label = labels[i], cohort = cell_k_val[i], s = cell_s[i], t = cell_t[i],
                           L = cell_L[i], n_cohort = n_k[cell_kk[i]], n_pair = cell_npair[i],
                           coordinates = coords[[i]], clipping = clip)
  }

  sel_desc <- if (need_q) list(
    model = scat$selection_model, pooling = scat$selection_pooling, universe = universe,
    design = list(d = design$d, d_s = design$d_s, has_intercept = design$has_intercept,
                  branch = design$branch, slope_names = design$slope_names),
    blocks = lapply(sel_blocks, function(bf) list(
      label = bf$label, branch = bf$branch, B_m = bf$B_m, d_s = bf$d_s, index = bf$index,
      observed_rate = bf$pi, converged = bf$diagnostics$converged, iter = bf$diagnostics$iter,
      warnings = bf$diagnostics$warnings)),
    roster = if (!identical(universe, "F")) list(cohorts = coh_vals, n = n_k, n_U = nU,
                                                 pi = pi_val, index = pi_index) else NULL,
    law = if (identical(universe, "F")) {
      if (identical(scat$selection_model, "independent_period_logit"))
        "q = q_s q_t (declared conditional date independence); q_safe = max(q, 1e-6) on the final product" else
        "q = Lambda(u); q_safe = max(q, 1e-6)"
    } else "q_safe = Pr(U | cohort) rho max(q_U, 1e-6); Pr(U | cohort) estimated, a_k not carried"
  ) else NULL
  imp_desc <- if (need_j) list(
    pooling = imputation_pooling, intercept = icpt, slope_columns = recipe$slope_columns,
    normalization = "d = (t - s) z(X) on the stored calendar index",
    groups = lapply(imp_groups, function(gr) list(
      label = gr$label, cells = labels[gr$cells], index = gr$index, p = gr$p,
      kept_columns = gr$kept_columns, aliased_columns = gr$aliased_columns, rss = gr$rss))
  ) else NULL

  diagnostics <- list(
    method = method, N = N, T = TT, n_cells = C, n_parameters = P,
    parameter_names = par_names, coefficients = coef_all,
    parameter_order = paste("selection blocks (local intercepts, slopes)",
                            "-> Pr(U | cohort) (A/B) -> imputation groups (local intercepts, kept slopes)"),
    selection = sel_desc, imputation = imp_desc, cells = cell_diag
  )
  # Post-fit charge: ACTUAL payload of the retained state (borrowed inputs
  # excluded) and of the returned diagnostics, plus the callback batch scratch.
  coord_max <- max(c(0L, lengths(coords)))
  bud <- etwfe_obs_budget(N, TT, B, K, C, P,
                          lapply(sel_blocks, function(bf) c(bf$B_m, bf$d_s)),
                          lapply(imp_groups, function(gr) c(length(gr$cells), gr$p)),
                          d_sel, p_X, coord_max)
  borrowed <- c("Y", "cohort", "X", "Z", "U", "rho")
  retained <- etwfe_obs_payload(st[setdiff(names(st), borrowed)]) + etwfe_obs_payload(diagnostics) + 512
  ipw_owner_safe_sum(c(retained, bud[["scratch"]]), workspace_entries,
                     "ETWFE observation provider workspace")
  diagnostics$budget <- list(
    retained_entries = retained, scratch_entries = unname(bud["scratch"]),
    preflight_retained_entries = unname(pre["retained"]),
    preflight_scratch_entries = unname(pre["scratch"]),
    borrowed_entries = as.double(N) * (TT + 1 + (if (need_q) d_sel + 2 else 0) +
                                        (if (need_j) ncol(X) else 0)),
    workspace_entries = as.double(workspace_entries), chunk_rows = B,
    note = paste("logical payload entries (cp04_payload convention), not RSS; Y/X/Z/U/rho/cohort",
                 "are borrowed; callback outputs (rows x coordinates or rows x n_parameters) are the caller's"))
  cb <- etwfe_obs_callbacks(st)
  result <- list(read = cb$read, score_rows = cb$score_rows, n_parameters = P, diagnostics = diagnostics)
  if (!is.null(report)) result$nuisance <- etwfe_nuisance_observation(st, scat, diagnostics, cb$score_rows, report)
  result
}
