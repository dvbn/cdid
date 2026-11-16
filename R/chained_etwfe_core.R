# Private ETWFE linear computational core (initial integration step).
#
# Scope: unweighted linear levels with unit fixed effects, common date fixed
# effects and saturated treated cohort-by-date deviations (etwfe with ivar and
# cgroup; jwdid with ivar and optional `never`, on aligned samples).
# Public validation, sample construction, clustered VCE, aggregation and the
# public entry point are owned by the root orchestrator. Nothing here is
# exported and no user-facing default is introduced. Base R only.
#
# Contract summary (details in notes.md):
#   etwfe_linear_core(y, unit, period, cohort, n_units, n_periods,
#                     control_group, compute_inffunc = TRUE)
#   -> list(targets = data.frame(cohort, period, is_pre, att, status),
#           inffunc = n_units x nrow(targets) matrix (N-scaled) or NULL)
#
# Algebra: every observed row (i, t) is assigned to exactly one slot k: the
# common date slot t (never-treated rows, NYT pre-treatment rows, the NT base
# row g - 1) or the free cohort-date slot (g, t). The fitted model is
# y_it = alpha_i + mu_k(i,t); the target is delta_gt = mu_(g,t) - mu_t.
# Unit intercepts are eliminated (M = within-unit transform):
#   G   = Z'MZ = diag(slot counts) - sum_i c_i c_i' / n_i
#   rhs = Z'My = rowsum(y - ybar_i, slot)
# c_i depends only on (cohort, observed-period mask), so G is accumulated
# once per distinct mask pattern. Estimability is decided structurally: slots
# are linked when some unit observes both; one gauge slot (lowest index, a
# date slot whenever the component has one) is fixed at zero per component,
# which makes the reduced G positive definite in exact arithmetic. A contrast
# is identified iff its cell slot is observed and lies in the same component as
# its date slot. One Cholesky factor serves the coefficients and the
# influence-function transfer H' = G^-1 C'.

# Working-memory budget (doubles) per block of the mask/residual passes
# (2^22 doubles = 32 MiB). Performance setting only; results do not depend on it.
etwfe_core_block_doubles <- 4194304

# A cohort's influence-function transfer uses dense zero-filled BLAS blocks when
# its observed share of unit-by-period cells is at least this value, else
# row-gather sums. Performance setting only (both paths are exact); unbenchmarked.
etwfe_core_dense_min_density <- 0.125

#' Private ETWFE linear core (see notes.md for the exact contract)
#' @noRd
etwfe_linear_core <- function(y, unit, period, cohort, n_units, n_periods,
                              control_group, compute_inffunc = TRUE, report = NULL) {
  etwfe_core_fit(y, unit, period, cohort, n_units, n_periods, control_group,
                 compute_inffunc,
                 block_doubles = etwfe_core_block_doubles,
                 dense_min_density = etwfe_core_dense_min_density, report = report)
}

#' Implementation with private performance hooks (used by developer checks)
#' @noRd
etwfe_core_fit <- function(y, unit, period, cohort, n_units, n_periods,
                           control_group, compute_inffunc,
                           block_doubles, dense_min_density, report = NULL) {
  never <- etwfe_core_contract(y, unit, period, cohort, n_units, n_periods,
                               control_group, compute_inffunc)
  n_units <- as.integer(n_units)
  n_t <- as.integer(n_periods)
  n_obs <- length(y)
  tt <- seq_len(n_t)

  ## Slot structure ---------------------------------------------------------
  tc <- sort.int(unique(cohort[cohort > 0L]))  # treated cohorts
  n_q <- length(tc)
  q_unit <- match(cohort, tc, nomatch = 0L)    # 0 = never treated
  # full_map[t, q + 1]: full slot of cohort index q at period t; date slot t
  # is full index t, cell (q, t) is n_t * q + t.
  full_map <- matrix(tt, n_t, n_q + 1L)
  if (n_q) {
    free_tq <- if (never) outer(tt, tc - 1L, "!=") else outer(tt, tc, ">=")
    cell_full <- matrix(n_t * rep(seq_len(n_q), each = n_t) + tt, n_t, n_q)
    full_map[, -1L] <- ifelse(free_tq, cell_full, tt)
  }
  slot_full <- full_map[cbind(period, q_unit[unit] + 1L)]
  n_full <- n_t * (n_q + 1L)
  slot_n <- tabulate(slot_full, n_full)
  obs_slots <- which(slot_n > 0L)
  p <- length(obs_slots)
  cmap <- integer(n_full)
  cmap[obs_slots] <- seq_len(p)
  slot <- cmap[slot_full]                       # compressed slot per row
  smap <- matrix(cmap[full_map], n_t, n_q + 1L) # 0 where the slot is unobserved

  ## Unit records (rows are sorted by unit, period) ---------------------------
  n_i <- tabulate(unit, n_units)
  present <- which(n_i > 0L)
  n_pres <- n_i[present]
  first_row <- cumsum(c(1L, n_pres))[seq_along(present)]
  q_pres <- q_unit[present]

  # One grouped pass: unit outcome sums plus exact observed-mask keys
  # (distinct powers of two below 2^52 per 52-period chunk sum exactly).
  n_chunk <- (n_t + 51L) %/% 52L
  bits <- 2^((period - 1L) %% 52L)
  if (n_chunk == 1L) {
    agg <- cbind(y, bits)
  } else {
    agg <- matrix(0, n_obs, n_chunk + 1L)
    agg[, 1L] <- y
    agg[cbind(seq_len(n_obs), (period - 1L) %/% 52L + 2L)] <- bits
  }
  if (n_obs) agg <- rowsum(agg, unit, reorder = TRUE) # rows follow `present`
  ybar <- numeric(n_units)
  ybar[present] <- agg[, 1L] / n_pres
  ytil <- y - ybar[unit]
  rhs <- if (n_obs) as.vector(rowsum(ytil, slot, reorder = TRUE)) else numeric()

  ## Distinct (cohort, mask) patterns ----------------------------------------
  keys <- lapply(seq_len(n_chunk) + 1L, function(j) agg[, j])
  o <- do.call(order, c(list(q_pres), keys, list(method = "radix")))
  n_p <- length(present)
  qo <- q_pres[o]
  brk <- qo[-1L] != qo[-n_p]
  for (k in keys) {
    ko <- k[o]
    brk <- brk | ko[-1L] != ko[-n_p]
  }
  pstart <- if (n_p) which(c(TRUE, brk)) else integer()
  pat_rep <- o[pstart]                          # representative (index into present)
  pat_m <- diff(c(pstart, n_p + 1L))            # units sharing the pattern
  pat_q <- q_pres[pat_rep]
  pat_n <- n_pres[pat_rep]
  pat_row <- first_row[pat_rep]

  ## G = Z'MZ and the slot link graph ----------------------------------------
  G <- matrix(0, p, p)
  diag(G) <- slot_n[obs_slots]
  link <- matrix(FALSE, p, p)
  for (qq in unique(pat_q)) {
    ip <- which(pat_q == qq)
    A <- matrix(0, n_t, n_t)
    cnt <- matrix(0, n_t, n_t)
    bl <- etwfe_core_blocks(rep.int(2 * n_t, length(ip)), block_doubles)
    for (b in seq_along(bl$start)) {
      jj <- ip[bl$start[b]:bl$end[b]]
      nb <- length(jj)
      idx <- rep.int(pat_row[jj] - 1L, pat_n[jj]) + sequence(pat_n[jj])
      mk <- matrix(0, nb, n_t)
      mk[cbind(rep.int(seq_len(nb), pat_n[jj]), period[idx])] <- 1
      A <- A + crossprod(mk * (pat_m[jj] / pat_n[jj]), mk)
      cnt <- cnt + crossprod(mk)
    }
    use <- which(smap[, qq + 1L] > 0L)
    s <- smap[use, qq + 1L]
    G[s, s] <- G[s, s] - A[use, use]
    link[s, s] <- link[s, s] | cnt[use, use] > 0
  }
  comp <- etwfe_core_components(link)
  keep <- which(duplicated(comp))               # drop one gauge slot per component

  ## Targets and structural status -------------------------------------------
  tq <- rep(seq_len(n_q), each = n_t)
  tper <- rep.int(tt, n_q)
  tg <- tc[tq]
  sel <- if (never) tper != tg - 1L else tper >= tg
  tq <- tq[sel]
  tper <- tper[sel]
  tg <- tg[sel]
  cell_c <- smap[cbind(tper, tq + 1L)]
  date_c <- smap[cbind(tper, rep.int(1L, length(tper)))]
  has <- cell_c > 0L
  ident <- has & date_c > 0L
  ident[ident] <- comp[cell_c[ident]] == comp[date_c[ident]]
  status <- rep.int("no_support", length(tper))
  status[has] <- "unidentified"
  status[ident] <- "identified"

  ## One factorization -------------------------------------------------------
  R <- etwfe_core_factor(G[keep, keep, drop = FALSE])
  mu <- numeric(p)
  if (length(keep)) {
    mu[keep] <- backsolve(R, backsolve(R, rhs[keep], transpose = TRUE))
  }
  att <- rep.int(NA_real_, length(tper))
  att[ident] <- mu[cell_c[ident]] - mu[date_c[ident]]
  targets <- data.frame(cohort = as.integer(tg), period = as.integer(tper),
                        is_pre = tper < tg, att = att, status = status,
                        stringsAsFactors = FALSE)
  result <- list(targets = targets, inffunc = NULL)
  if (!is.null(report)) {
    slots <- etwfe_nuisance_slot_codes(obs_slots, n_t, c(0L, tc))
    slots$gauge <- !duplicated(comp)
    first <- cumsum(c(1, as.double(n_i)))[seq_len(n_units)]
    beta <- matrix(mu, p, 1L)
    solve_factor <- function(z) cp04_solve(R, z)
    score <- function(rows) etwfe_nuisance_within_score(rows, y, slot, NULL,
      beta, n_i, first, keep, solve_factor)
    spare <- report$workspace_entries - 5 * as.double(p)^2 - 6 * length(y) - 8 * as.double(n_units)
    result$nuisance <- list(outcome = etwfe_nuisance_outcome(mu, keep, matrix(seq_len(p), p, 1L),
      slots, NULL, numeric(), numeric(), n_units, report, score, spare, 4 * p + 12 * n_t,
      sum(n_i[n_i > 1L]), "levels", "Observed rows on repeat-observed units",
      "Within-OLS score and fitted slot factor on original units."))
  }
  if (!compute_inffunc) return(result)

  ## Influence functions ------------------------------------------------------
  id <- which(ident)
  n_id <- length(id)
  Ht <- matrix(0, p + 1L, n_id)                 # row 1: zero row for slot 0
  if (n_id) {
    Ct <- matrix(0, p, n_id)
    Ct[cbind(cell_c[id], seq_len(n_id))] <- 1
    Ct[cbind(date_c[id], seq_len(n_id))] <- -1
    Ht[keep + 1L, ] <- backsolve(R, backsolve(R, Ct[keep, , drop = FALSE],
                                              transpose = TRUE))
  }
  # Residuals of the FE fit; they sum to zero within unit, so Z'e = (MZ)'e is
  # the within-transformed score of each original unit.
  fit <- mu[slot]
  fbar <- numeric(n_units)
  if (n_obs) fbar[present] <- as.vector(rowsum(fit, unit, reorder = TRUE)) / n_pres
  ehat <- ytil - fit + fbar[unit]
  inffunc <- etwfe_core_inffunc(ehat, slot, unit, period, Ht, smap, q_pres,
                                present, n_pres, first_row, n_units, n_t, id,
                                length(tper), block_doubles, dense_min_density)
  result$inffunc <- inffunc
  result
}

#' inffunc[i, j] = n_units * (C G^-1 Z_i' e_i)_j for identified j, 0 for units
#' without rows, NA columns for non-identified targets.
#' @noRd
etwfe_core_inffunc <- function(ehat, slot, unit, period, Ht, smap, q_pres,
                               present, n_pres, first_row, n_units, n_t, id,
                               n_targets, block_doubles, dense_min_density) {
  n_id <- length(id)
  inf <- matrix(0, n_units, n_targets)
  if (n_id < n_targets) {
    na_col <- rep.int(TRUE, n_targets)
    na_col[id] <- FALSE                         # (-id would select nothing if id is empty)
    inf[, na_col] <- NA_real_
  }
  if (!n_id) return(inf)
  for (qq in sort.int(unique(q_pres))) {
    js <- which(q_pres == qq)
    nj <- n_pres[js]
    if (sum(nj) >= dense_min_density * length(js) * n_t) {
      # Dense: cohort residual block (units x periods) times the cohort's
      # period-to-target transfer (periods x targets).
      Hq <- Ht[smap[, qq + 1L] + 1L, , drop = FALSE]
      bl <- etwfe_core_blocks(rep.int(n_t + 2 * n_id, length(js)), block_doubles)
      for (b in seq_along(bl$start)) {
        jb <- js[bl$start[b]:bl$end[b]]
        nb <- length(jb)
        njb <- n_pres[jb]
        idx <- rep.int(first_row[jb] - 1L, njb) + sequence(njb)
        if (all(njb == n_t)) {
          E <- matrix(ehat[idx], nb, n_t, byrow = TRUE)
        } else {
          E <- matrix(0, nb, n_t)
          E[cbind(rep.int(seq_len(nb), njb), period[idx])] <- ehat[idx]
        }
        inf[present[jb], id] <- n_units * (E %*% Hq)
      }
    } else {
      # Sparse: gather transfer rows per observed row, sum within unit.
      bl <- etwfe_core_blocks(nj * (2 * n_id + 1), block_doubles)
      for (b in seq_along(bl$start)) {
        jb <- js[bl$start[b]:bl$end[b]]
        njb <- n_pres[jb]
        idx <- rep.int(first_row[jb] - 1L, njb) + sequence(njb)
        inf[present[jb], id] <- n_units *
          rowsum(Ht[slot[idx] + 1L, , drop = FALSE] * ehat[idx], unit[idx],
                 reorder = TRUE)
      }
    }
  }
  inf
}

#' Contiguous blocks whose summed cost stays near `budget`
#' @noRd
etwfe_core_blocks <- function(cost, budget) {
  n <- length(cost)
  if (!n) return(list(start = integer(), end = integer()))
  b <- (cumsum(as.double(cost)) - 1) %/% budget
  end <- c(which(b[-1L] != b[-n]), n)
  list(start = c(1L, end[-length(end)] + 1L), end = end)
}

#' Connected components of a symmetric logical adjacency (breadth-first)
#' @noRd
etwfe_core_components <- function(link) {
  p <- nrow(link)
  comp <- integer(p)
  nc <- 0L
  for (s in seq_len(p)) {
    if (comp[s]) next
    nc <- nc + 1L
    comp[s] <- nc
    front <- s
    while (length(front)) {
      reach <- if (length(front) == 1L) link[front, ] else
        colSums(link[front, , drop = FALSE]) > 0
      front <- which(reach & comp == 0L)
      comp[front] <- nc
    }
  }
  comp
}

#' Cholesky factor of the gauge-fixed G. Positive definiteness is guaranteed
#' structurally by the graph gauge, so failure here is numerical, not a rank
#' decision.
#' @noRd
etwfe_core_factor <- function(G) {
  if (!nrow(G)) return(G)
  R <- tryCatch(chol(G), error = function(e) NULL)
  rc <- if (is.null(R)) NA_real_ else rcond(R, triangular = TRUE)
  if (!is.finite(rc) || rc * rc <= .Machine$double.eps) {
    cdid_abort(paste0(
      "ETWFE core: the structurally identified fixed-effect system is ",
      "numerically singular (Cholesky failed or reciprocal condition <= ",
      "machine epsilon); no estimates are returned."),
      "cdid_error_etwfe_numeric")
  }
  R
}

#' Cheap private-contract guard; public validation is owned by the caller.
#' Returns TRUE for "nevertreated", FALSE for "notyettreated".
#' @noRd
etwfe_core_contract <- function(y, unit, period, cohort, n_units, n_periods,
                                control_group, compute_inffunc) {
  never <- if (identical(control_group, "nevertreated")) TRUE else
    if (identical(control_group, "notyettreated")) FALSE else NA
  n <- length(y)
  ok <- !is.na(never) &&
    (isTRUE(compute_inffunc) || isFALSE(compute_inffunc)) &&
    is.numeric(y) && is.integer(unit) && is.integer(period) &&
    is.integer(cohort) &&
    length(n_units) == 1L && length(n_periods) == 1L &&
    isTRUE(n_units >= 1 && n_units == trunc(n_units)) &&
    isTRUE(n_periods >= 1 && n_periods == trunc(n_periods)) &&
    length(unit) == n && length(period) == n &&
    length(cohort) == n_units &&
    all(is.finite(y)) &&
    isTRUE(all(unit >= 1L & unit <= n_units)) &&
    isTRUE(all(period >= 1L & period <= n_periods)) &&
    isTRUE(all(cohort == 0L | (cohort >= 2L & cohort <= n_periods))) &&
    etwfe_core_sorted(unit, period)
  if (!ok) {
    cdid_abort(paste0(
      "ETWFE core: inputs violate the private contract (integer unit/period ",
      "rows sorted by unit then period without duplicates, integer cohort ",
      "per unit in {0, 2..n_periods}, finite numeric outcomes, ",
      "control_group 'nevertreated' or 'notyettreated')."),
      "cdid_error_etwfe_contract")
  }
  never
}

#' @noRd
etwfe_core_sorted <- function(unit, period) {
  n <- length(unit)
  if (n < 2L) return(TRUE)
  du <- unit[-1L] - unit[-n]
  isTRUE(all(du > 0L | (du == 0L & period[-1L] > period[-n])))
}
