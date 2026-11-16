# Cohort-time treatment moments under selection = "none", unit weights, zero
# anticipation and time-invariant covariates.
#
# A link is one date pair s < t of target cohort g (calendar indices).
#   target : raw mean of y_t - y_s over target units observing both dates.
#   control: donor cohort k observing both dates (NT: k = 0; NYT: k != g with
#            k = 0 or k > t). Value = full-target-frame mean of the fitted
#            common increment m(W; s, t) plus the donor odds-weighted mean of
#            dy - m, odds from the target-vs-k logit on the full two-cohort
#            frame. IPW sets m = 0 and fits no outcome model.
# m(W; s, t) = W'(b_t - b_s), W = (1, Z), is fitted once by unweighted least
# squares over all observed clean pairs (NT: never units; NYT: every cohort
# before adoption). For one unit the pair loss equals m_i times its
# within-unit sum of squares on its m_i clean dates, so the Gram is built from
# dense unit-by-date blocks, never from pair rows.
# Unobserved cells are zero-filled with a zero mask (inputs are finite), so no
# arithmetic ever multiplies 0 by NA. score_rows(rows) returns N-scaled
# original-unit scores, length(rows) x nrow(links), streamed in row blocks.

#' @noRd
etwfe_treatment_moments <- function(y, unit, period, cohort, n_units, n_periods, Z,
                                    control_group, est_method, compute_scores,
                                    workspace_entries, glm_control, report = NULL) {
  never <- etwfe_core_contract(y, unit, period, cohort, n_units, n_periods,
                               control_group, compute_scores)
  if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != n_units || any(!is.finite(Z)) ||
      !(identical(est_method, "ipw") || identical(est_method, "dr"))) {
    etwfe_tm_fail("Z must be a finite numeric unit-frame matrix and est_method 'ipw' or 'dr'.",
                  "cdid_error_etwfe_contract")
  }
  if (!is.numeric(workspace_entries) || length(workspace_entries) != 1L ||
      is.na(workspace_entries) || workspace_entries < 1) etwfe_tm_workspace_error()
  glm_control <- validate_ipw_glm_control(glm_control)
  dr <- identical(est_method, "dr")
  N <- as.integer(n_units)
  n_t <- as.integer(n_periods)
  q <- ncol(Z) + 1L
  level <- sort(unique(cohort))                 # never (0) first when present
  roster <- split(seq_len(N), factor(cohort, levels = level))
  n_roster <- lengths(roster, use.names = FALSE)
  n_obs <- tabulate(unit, N)
  pn <- list(y = y, period = period, cohort = cohort, never = never, n_t = n_t,
             n_obs = n_obs, first = cumsum(c(1, as.double(n_obs)))[seq_len(N)])

  # Workspace in doubles: retained/fixed arrays plus one block row of scratch.
  per_unit <- n_t * (3 * q + 8) + 4 * q
  fixed <- 4 * length(y) + N * (q + 10) + 3 * length(level) * n_t^2 + 1024
  if (fixed + per_unit > workspace_entries) etwfe_tm_workspace_error()
  chunks <- function(ids) {
    size <- max(1, floor(min(workspace_entries - fixed, 2^22) / per_unit))
    split(ids, (seq_along(ids) - 1L) %/% size)
  }

  # Centred, amplitude-scaled calendar coordinates (predictions invariant).
  # The treatment logit uses the supplied design cbind(1, Z).
  center <- colMeans(Z)
  amplitude <- rep(1, q - 1L)
  W <- matrix(1, N, q)
  for (a in seq_len(q - 1L)) {
    z <- Z[, a] - center[a]
    if (max(abs(z)) > 0) amplitude[a] <- max(abs(z))
    W[, a + 1L] <- z / amplitude[a]
  }

  ## Pair support and target change sums per cohort ---------------------------
  pair_n <- pair_dy <- vector("list", length(level))
  for (c in seq_along(level)) {
    pair_n[[c]] <- pair_dy[[c]] <- matrix(0, n_t, n_t)
    members <- roster[[c]]
    for (ids in chunks(members[n_obs[members] > 1L])) {
      b <- etwfe_tm_read(pn, ids)
      pair_n[[c]] <- pair_n[[c]] + crossprod(b$O)
      pair_dy[[c]] <- pair_dy[[c]] + crossprod(b$O, b$Y)
    }
    pair_dy[[c]] <- pair_dy[[c]] - t(pair_dy[[c]])  # [s, t]: sum of y_t - y_s
  }

  ## Links: every supported target edge and eligible donor edge ---------------
  upper <- which(upper.tri(diag(n_t)), arr.ind = TRUE)
  edges <- function(g, k) {                     # k = NA: the target's own edges
    if (is.na(k)) return(pair_n[[match(g, level)]][upper] > 0)
    pair_n[[match(k, level)]][upper] > 0 & (k == 0L | upper[, 2L] < k)
  }
  role_g <- role_k <- integer()
  for (g in level[level > 0L]) {
    k <- c(NA_integer_, if (never) level[level == 0L] else level[level != g])
    role_g <- c(role_g, rep.int(g, length(k)))
    role_k <- c(role_k, k)
  }
  # Count first; positions are materialized only after the workspace check.
  n_hit <- vapply(seq_along(role_g), function(r) sum(edges(role_g[r], role_k[r])), 0)
  n_link <- sum(n_hit)
  is_model <- n_hit > 0 & !is.na(role_k)
  model_g <- role_g[is_model]
  model_k <- role_k[is_model]
  n_model <- length(model_g)
  n_coef <- n_t * q
  frame_max <- max(c(0, n_roster[match(model_g, level)] + n_roster[match(model_k, level)]))
  fixed <- fixed + n_link * (16 + 2 * q) + n_model * (2 * q^2 + 8) +
    (3 * q + 2) * n_t^2 + 12 * frame_max * (q + 2) +
    (if (dr && n_model) 4 * n_coef^2 + 4 * n_coef else 0)
  score_row <- per_unit + 2 * n_coef + 5 * n_link
  if (fixed + per_unit > workspace_entries ||
      (compute_scores && fixed + n_link + score_row > workspace_entries)) {
    etwfe_tm_workspace_error()
  }
  j <- as.integer(unlist(lapply(which(n_hit > 0), function(r) which(edges(role_g[r], role_k[r])))))
  g <- rep(role_g, n_hit)
  donor <- rep(role_k, n_hit)
  s <- upper[j, 1L]
  t <- upper[j, 2L]
  o <- order(g, s, t, !is.na(donor), donor)
  g <- g[o]; s <- s[o]; t <- t[o]; donor <- donor[o]
  target <- is.na(donor)
  links <- data.frame(g = g, s = s, t = t, ref = g - 1L,
                      link_role = c("control", "target")[target + 1L],
                      functional_key = c("calendar", "target")[target + 1L],
                      donor = donor, stringsAsFactors = FALSE)
  model <- match(g * (n_t + 1L) + donor, model_g * (n_t + 1L) + model_k)  # NA: target
  target_cols <- split(which(target), factor(g[target], levels = level))

  mhat <- support <- mass <- resid_mean <- numeric(n_link)
  for (c in which(level > 0L)) {
    cols <- target_cols[[c]]
    st <- cbind(s[cols], t[cols])
    support[cols] <- mass[cols] <- pair_n[[c]][st]
    mhat[cols] <- pair_dy[[c]][st] / mass[cols]
  }

  ## Shared untreated calendar increment (DR with a donor edge only) ----------
  B <- matrix(0, n_t, q)                        # date levels, gauge rows zero
  cal <- NULL
  if (dr && n_model) {
    clean_ids <- which(n_obs > 1L & (cohort == 0L | (!never & cohort > 2L)))
    cal <- etwfe_tm_calendar(pn, W, chunks(clean_ids), c(s[!target], t[!target]), count_rows = !is.null(report))
    B <- cal$B
  }

  ## Treatment odds per used (target, donor) and donor moments ----------------
  Wbar <- rowsum(W, cohort, reorder = TRUE) / n_roster  # rows follow `level`
  gamma <- imbalance <- if (compute_scores) matrix(0, q, n_link) else NULL
  rep_t <- rep(seq_len(n_t), q)                 # feature-major (a - 1) n_t + d
  rep_f <- rep(seq_len(q), each = n_t)
  fits <- vector("list", n_model)
  for (m in seq_len(n_model)) {
    cg <- match(model_g[m], level)
    ck <- match(model_k[m], level)
    ids <- c(roster[[cg]], roster[[ck]])
    fit <- cp04_logit(etwfe_tm_design(Z, ids), as.double(cohort[ids] == model_g[m]),
                      rep(1, length(ids)), glm_control, N,
                      list(model = "etwfe_treatment_odds", g = model_g[m], donor = model_k[m]))
    cols <- which(model == m)
    st <- cbind(s[cols], t[cols])
    M0 <- Den <- matrix(0, n_t, n_t)
    M1 <- M2 <- M3 <- NULL
    if (compute_scores) M1 <- M2 <- M3 <- matrix(0, n_t * q, n_t)
    n_capped <- 0
    members <- roster[[ck]]
    for (ids in chunks(members[n_obs[members] > 1L])) {
      b <- etwfe_tm_read(pn, ids)
      X <- etwfe_tm_design(Z, ids)
      lv <- ipw_logit_values(X, fit$coefficients)
      n_capped <- n_capped + sum(lv$active == 0)
      E <- if (dr) b$O * (b$Y - tcrossprod(W[ids, , drop = FALSE], B)) else b$Y
      Ow <- b$O * lv$odds
      M0 <- M0 + crossprod(Ow, E)
      Den <- Den + crossprod(Ow, b$O)
      if (compute_scores) {
        OX <- Ow[, rep_t, drop = FALSE] * (X * lv$active)[, rep_f, drop = FALSE]
        M1 <- M1 + crossprod(OX, E)
        M2 <- M2 + crossprod(OX, b$O)
        if (dr) {
          M3 <- M3 + crossprod(Ow[, rep_t, drop = FALSE] * W[ids, rep_f, drop = FALSE], b$O)
        }
      }
    }
    support[cols] <- pair_n[[ck]][st]
    mass[cols] <- Den[st]
    resid_mean[cols] <- (M0[st] - M0[st[, 2:1, drop = FALSE]]) / mass[cols]
    increment <- B[t[cols], , drop = FALSE] - B[s[cols], , drop = FALSE]
    mhat[cols] <- drop(increment %*% Wbar[cg, ]) + resid_mean[cols]
    if (compute_scores) {
      for (a in seq_len(q)) {
        at <- cbind((a - 1L) * n_t + s[cols], t[cols])
        ta <- cbind((a - 1L) * n_t + t[cols], s[cols])
        # Active odds derivative of the normalized donor residual mean.
        gamma[a, cols] <- (M1[at] - M1[ta] - resid_mean[cols] * M2[at]) / mass[cols]
        if (dr) imbalance[a, cols] <- Wbar[cg, a] - M3[at] / mass[cols]
      }
    }
    fd <- fit$diagnostics
    fits[[m]] <- list(g = model_g[m], k = model_k[m], cg = cg, cols = cols,
                      coefficients = fit$coefficients, H_chol = fit$H_chol,
                      n_fit = n_roster[cg] + n_roster[ck], n_capped = n_capped,
                      iter = etwfe_tm_scalar(fd$iter), converged = etwfe_tm_scalar(fd$converged))
  }
  if (any(!is.finite(mhat))) {
    etwfe_tm_fail("nonfinite moment values.", "cdid_error_etwfe_numeric")
  }

  diagnostics <- list(
    est_method = est_method, control_group = control_group,
    n_units = N, n_periods = n_t, n_features = q,
    link_support = data.frame(n_units = support, mass = mass, residual_mean = resid_mean),
    odds = data.frame(g = model_g, donor = model_k,
                      n_fit = vapply(fits, `[[`, 0, "n_fit"),
                      iter = vapply(fits, `[[`, 0, "iter"),
                      converged = vapply(fits, `[[`, 0, "converged"),
                      n_capped = vapply(fits, `[[`, 0, "n_capped")),
    calendar = if (!is.null(cal)) list(
      levels = cal$B, center = center, amplitude = amplitude,
      component = cal$component, gauge = cal$gauge, n_coef = length(cal$keep),
      rcond = cal$rcond, n_units = cal$n_units),
    workspace = list(entries = workspace_entries, fixed = fixed,
                     block_row = per_unit, score_row = score_row))
  score_rows <- NULL
  if (compute_scores) {
    score_rows <- etwfe_tm_scorer(list(
      pn = pn, Z = Z, W = W, B = B, Wbar = Wbar, n_roster = n_roster, level = level,
      s = s, t = t, target = target, target_cols = target_cols, mhat = mhat,
      resid_mean = resid_mean, mass = mass, gamma = gamma, imbalance = imbalance,
      fits = lapply(fits, function(f) f[c("g", "k", "cg", "cols", "coefficients", "H_chol")]),
      cal = cal[c("keep", "R", "scale")], dr = dr, N = N, q = q,
      retained = fixed, score_row = score_row, entries = workspace_entries))
  }
  result <- list(links = links, mhat = mhat, score_rows = score_rows, diagnostics = diagnostics)
  if (!is.null(report)) {
    spare <- workspace_entries - fixed
    result$nuisance <- list(treatment = etwfe_nuisance_treatment(fits, model_g, model_k,
      cohort, Z, links, report, spare))
    if (!is.null(cal)) {
      slots <- etwfe_nuisance_slot_codes(seq_len(n_t), n_t, 0L)
      slots$gauge <- slots$slot %in% cal$gauge
      score <- function(rows) {
        b <- etwfe_tm_read(pn, rows)
        w <- W[rows, , drop = FALSE]
        predict <- tcrossprod(w, B)
        etwfe_tm_calendar_score(pn, rows, N, q, cal, b, w, predict)
      }
      j <- which(!target)
      comparisons <- etwfe_nuisance_comparisons(g[j], s[j], t[j], donor[j],
        rep("calendar", length(j)), report)
      result$nuisance$outcome <- etwfe_nuisance_outcome(as.vector(t(B)), cal$keep,
        t(matrix(seq_len(n_t * q), q)), slots, Z, center, amplitude, N, report, score,
        spare, 4 * n_t * q + 8 * n_t, cal$n_pairs, "calendar_changes", "Clean observed pairs",
        "Clean-pair OLS score and fitted calendar factor on original units.", comparisons)
    }
  }
  result
}

#' Shared calendar regression: Gram and right side from clean unit blocks,
#' one full-vector gauge per connected calendar component, covariate-core rank
#' policy (Jacobi-scaled Cholesky, rcond^2 > eps), no ridge fallback.
#' @noRd
etwfe_tm_calendar <- function(pn, W, blocks, used_dates, count_rows = FALSE) {
  n_t <- pn$n_t
  q <- ncol(W)
  n_coef <- n_t * q
  rep_d <- rep(seq_len(n_t), each = q)          # coefficient (d - 1) q + a
  rep_a <- rep(seq_len(q), n_t)
  gram <- matrix(0, n_coef, n_coef)
  rhs <- numeric(n_coef)
  link <- matrix(FALSE, n_t, n_t)
  n_fit <- 0
  pair_count <- if (count_rows) matrix(0, n_t, n_t)
  for (ids in blocks) {
    b <- etwfe_tm_read(pn, ids)
    Cm <- b$O * etwfe_tm_clean(pn, ids)
    m <- rowSums(Cm)
    use <- m > 1
    if (!any(use)) next
    Cm <- Cm[use, , drop = FALSE]
    m <- m[use]
    w <- W[ids[use], , drop = FALSE]
    yc <- b$Y[use, , drop = FALSE] * Cm
    n_fit <- n_fit + length(m)
    connected <- crossprod(Cm)
    link <- link | connected > 0
    if (count_rows) pair_count <- pair_count + connected
    # sum over pairs of (u_t - u_s)(u_t - u_s)' = m D - c c', kron W W'.
    K <- Cm[, rep_d, drop = FALSE] * w[, rep_a, drop = FALSE]
    gram <- gram - crossprod(K)
    for (d in seq_len(n_t)) {
      jj <- (d - 1L) * q + seq_len(q)
      gram[jj, jj] <- gram[jj, jj] + crossprod(w, w * (m * Cm[, d]))
    }
    rhs <- rhs + colSums(K * (m * yc - rowSums(yc))[, rep_d, drop = FALSE])
  }
  # The system is block diagonal across date components; only components that
  # carry a donor edge enter it.
  comp <- etwfe_core_components(link)
  used <- comp %in% comp[used_dates]
  keep <- which((duplicated(comp) & used)[rep_d])
  G <- gram[keep, keep, drop = FALSE]
  rank_error <- function() {
    etwfe_tm_fail(paste("the shared untreated calendar regression is rank deficient",
                        "after fixing one date gauge per calendar component."),
                  "cdid_error_etwfe_numeric")
  }
  if (!length(keep) || any(!is.finite(G)) || any(diag(G) <= 0)) rank_error()
  scale <- sqrt(diag(G))
  R <- tryCatch(chol(t(t(G / scale) / scale)), error = function(e) NULL)
  rc <- if (is.null(R)) NA_real_ else rcond(R, triangular = TRUE)
  if (!is.finite(rc) || rc^2 <= .Machine$double.eps) rank_error()
  beta <- numeric(n_coef)
  beta[keep] <- cp04_solve(R, rhs[keep] / scale) / scale
  result <- list(B = matrix(beta, n_t, q, byrow = TRUE), keep = keep, R = R, scale = scale,
    component = comp, gauge = which(!duplicated(comp) & used), rcond = rc, n_units = n_fit)
  if (count_rows) result$n_pairs <- sum(pair_count[upper.tri(pair_count) & outer(used, used, "&")])
  result
}

#' Score closure over a compact state; no N-by-links matrix is retained.
#' @noRd
etwfe_tm_scorer <- function(st) {
  force(st)
  function(rows) {
    if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || anyNA(rows) ||
        any(rows != trunc(rows) | rows < 1 | rows > st$N) || anyDuplicated(rows)) {
      etwfe_tm_fail("score rows must be unique original 1-based unit indices.",
                    "cdid_error_etwfe_contract")
    }
    used <- st$retained + as.double(length(rows)) * length(st$mhat)
    if (used + st$score_row > st$entries) etwfe_tm_workspace_error()
    size <- max(1, floor(min(st$entries - used, 2^22) / st$score_row))
    out <- matrix(0, length(rows), length(st$mhat))
    for (pos in split(seq_along(rows), (seq_along(rows) - 1L) %/% size)) {
      out[pos, ] <- etwfe_tm_score_block(st, as.integer(rows[pos]))
    }
    out
  }
}

#' N-scaled scores of one row block, in link order.
#' @noRd
etwfe_tm_score_block <- function(st, ids) {
  N <- st$N
  s <- st$s
  t <- st$t
  nb <- length(ids)
  out <- matrix(0, nb, length(st$mhat))
  b <- etwfe_tm_read(st$pn, ids)
  k_of <- st$pn$cohort[ids]
  w <- st$W[ids, , drop = FALSE]
  P <- tcrossprod(w, st$B)                      # zero for IPW
  E <- if (st$dr) b$O * (b$Y - P) else b$Y
  # Normalized observed-pair deviation: N mult 1{s,t} (V_t - V_s - centre) / mass.
  direct <- function(r, cols, V, centre, mult) {
    hit <- b$O[r, s[cols], drop = FALSE] * b$O[r, t[cols], drop = FALSE]
    dev <- V[r, t[cols], drop = FALSE] - V[r, s[cols], drop = FALSE] -
      rep(centre[cols], each = length(r))
    N * mult * hit * dev / rep(st$mass[cols], each = length(r))
  }
  for (c in which(st$level > 0L)) {
    r <- which(k_of == st$level[c])
    cols <- st$target_cols[[c]]
    if (length(r) && length(cols)) out[r, cols] <- direct(r, cols, b$Y, st$mhat, 1)
  }
  for (f in st$fits) {
    r <- which(k_of == f$g | k_of == f$k)
    if (!length(r)) next
    cols <- f$cols
    X <- etwfe_tm_design(st$Z, ids[r])
    lv <- ipw_logit_values(X, f$coefficients)
    treated <- k_of[r] == f$g
    # Logit coefficient IF (H = information / N) times the active odds gradient.
    u <- X * (treated - lv$prob)
    out[r, cols] <- out[r, cols] +
      t(cp04_solve(f$H_chol, t(u))) %*% st$gamma[, cols, drop = FALSE]
    rk <- r[!treated]                           # donor residual, odds-normalized
    if (length(rk)) {
      out[rk, cols] <- out[rk, cols] + direct(rk, cols, E, st$resid_mean, lv$odds[!treated])
    }
    rg <- r[treated]                            # full-target-frame increment mean
    if (st$dr && length(rg)) {
      centred <- w[rg, , drop = FALSE] - rep(st$Wbar[f$cg, ], each = length(rg))
      increment <- st$B[t[cols], , drop = FALSE] - st$B[s[cols], , drop = FALSE]
      out[rg, cols] <- out[rg, cols] + (N / st$n_roster[f$cg]) * tcrossprod(centred, increment)
    }
  }
  if (!is.null(st$cal$keep)) {
    # N G^-1 times each unit's clean-pair score W_i (x) (m_i e_id - sum_d e_id).
    q <- st$q
    phi <- matrix(0, nb, st$pn$n_t * q)
    phi[, st$cal$keep] <- etwfe_tm_calendar_score(st$pn, ids, N, q, st$cal, b, w, P)
    cols <- which(!st$target)
    for (a in seq_len(q)) {
      out[, cols] <- out[, cols] +
        (phi[, (t[cols] - 1L) * q + a, drop = FALSE] - phi[, (s[cols] - 1L) * q + a, drop = FALSE]) *
        rep(st$imbalance[a, cols], each = nb)
    }
  }
  if (any(!is.finite(out))) etwfe_tm_fail("nonfinite scores.", "cdid_error_etwfe_numeric")
  out
}

# The fitted date-major calendar coordinates, without gauge columns.
etwfe_tm_calendar_score <- function(pn, ids, N, q, cal, b, w, P) {
  rep_d <- rep(seq_len(pn$n_t), each = q)
  rep_a <- rep(seq_len(q), pn$n_t)
  Cm <- b$O * etwfe_tm_clean(pn, ids)
  e <- Cm * (b$Y - P)
  u <- Cm * (rowSums(Cm) * e - rowSums(e))
  keep <- cal$keep
  score <- u[, rep_d[keep], drop = FALSE] * w[, rep_a[keep], drop = FALSE]
  N * t(cp04_solve(cal$R, t(score) / cal$scale) / cal$scale)
}

#' Dense zero-filled outcome block and observation mask (units x dates).
#' @noRd
etwfe_tm_read <- function(pn, ids) {
  n <- pn$n_obs[ids]
  rows <- rep(pn$first[ids] - 1, n) + sequence(n)
  at <- cbind(rep(seq_along(ids), n), pn$period[rows])
  O <- Y <- matrix(0, length(ids), pn$n_t)
  O[at] <- 1
  Y[at] <- pn$y[rows]
  list(O = O, Y = Y)
}

#' Untreated (clean) dates: never units always, NYT cohorts before adoption.
#' @noRd
etwfe_tm_clean <- function(pn, ids) {
  outer(pn$cohort[ids], seq_len(pn$n_t), function(k, d) k == 0L | (!pn$never & d < k))
}

#' @noRd
etwfe_tm_design <- function(Z, ids) cbind(1, unname(Z[ids, , drop = FALSE]))

#' @noRd
etwfe_tm_scalar <- function(x) if (length(x) == 1L) as.double(x) else NA_real_

#' @noRd
etwfe_tm_fail <- function(message, class) {
  cdid_abort(paste("ETWFE treatment moments:", message), class)
}

#' @noRd
etwfe_tm_workspace_error <- function() {
  etwfe_tm_fail("the request exceeds the available numeric workspace.",
                "cdid_etwfe_workspace_error")
}
