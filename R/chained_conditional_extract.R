# CP08 O14: per-backend process extractors of the chained conditional pretest.
#
# Each declared cell is one auxiliary link of the fitted union registry. Its
# extractor returns the (S4) pieces of R/chained_conditional_engine.R for
# bounded chunks of original rows: term weights nu_r and responses rho_r (the
# link-level outcome change, residual or imputed/augmented response that h_u
# multiplies), the N-normalized nuisance IF rows xi of every coordinate the
# link reads (rates, selection, propensity, imputation and the global OR
# slope with all their upstream forcing, taken from the backend's own
# details-mode IF evaluator), and the loading integrands of lambda = d log nu
# and pi = d rho. No nuisance is refit and no row is masked: the shared fits
# of the union catalog keep their assumed models; only loadings depend on u.
# Under probability clipping, nu carries the floored observation denominators
# and capped odds, and lambda the active masks of the same derivatives.

cpt_piece <- function(rows, nu, rho, lam, pim, sD, xi) {
  n <- length(rows); R <- ncol(nu); Dl <- ncol(xi)
  P <- matrix(0, n, Dl); blam <- matrix(0, Dl, R)
  for (r in seq_len(R)) {
    if (!is.null(lam[[r]])) {
      P <- P + (sD[r] * nu[, r] * rho[, r]) * lam[[r]]
      blam[, r] <- colSums(nu[, r] * lam[[r]])
    }
    if (!is.null(pim[[r]])) P <- P + (sD[r] * nu[, r]) * pim[[r]]
  }
  list(rows = as.integer(rows), nu = nu, rho = rho, P = P, xi = xi, blam = blam)
}

# Treated/control weighted ratio pair (IPW and DR families). Coordinates are
# [upstream (n_up), propensity (p_trt), OR intercept (DR), OR slopes (DR)].
# `active` is the odds mask 1{p < 1 - eps}: d log o_c = active X' d beta.
cpt_two_arm_piece <- function(rows, rr, D, fw, X, odds, active, rho, lam_up, xi, D_T, D_M, n_up, p_trt,
                              d_or = NULL, or_icpt = FALSE) {
  n <- length(rows); Dl <- ncol(xi)
  nu <- rh <- matrix(0, n, 2L)
  lamT <- lamM <- matrix(0, n, Dl)
  if (length(rr)) {
    nu[rr, 1L] <- fw * D
    nu[rr, 2L] <- fw * (1 - D) * odds
    rh[rr, 1L] <- rho
    rh[rr, 2L] <- rho
    if (n_up) {
      lamT[rr, seq_len(n_up)] <- lam_up
      lamM[rr, seq_len(n_up)] <- lam_up
    }
    if (p_trt) lamM[rr, n_up + seq_len(p_trt)] <- active * X
  }
  pim <- NULL
  if (!is.null(d_or)) {
    pim <- matrix(0, n, Dl); off <- n_up + p_trt
    if (or_icpt) {
      if (length(rr)) pim[rr, off + 1L] <- -1
      off <- off + 1L
    }
    if (ncol(d_or) && length(rr)) pim[rr, off + seq_len(ncol(d_or))] <- -d_or
  }
  cpt_piece(rows, nu, rh, list(lamT, lamM), list(pim, pim), c(1 / D_T, -1 / D_M), xi)
}

cpt_backend_kind <- function(backend) {
  if (inherits(backend, "cdid_reference_ipw_backend")) return(if (isTRUE(backend$custom)) "custom" else "reference_ipw")
  if (inherits(backend, "cdid_cp04_custom_backend")) return("custom")
  if (inherits(backend, "cdid_augmented_backend")) return("augmented")
  if (inherits(backend, "cdid_imputation_backend")) return("imputation")
  if (inherits(backend, "cdid_dr_baseline")) return("dr_baseline")
  if (inherits(backend, "cdid_dr_corollary1")) return("dr_selection")
  if (inherits(backend, "cdid_or_baseline")) return("or_baseline")
  if (inherits(backend, "cdid_or_corollary1"))
    return(if (isTRUE(backend$uses_selection_models)) "or_selection" else "or_legacy")
  if (inherits(backend, "cdid_cp04_corollary1"))
    return(if (isTRUE(backend$uses_selection_models)) "ipw_selection" else "ipw_legacy")
  "unsupported"
}

# Details-mode rows of the backend's own complete IF evaluator (frame scale).
cpt_backend_details <- function(backend, kind, rows, j) {
  rows <- as.integer(rows); j <- as.integer(j)
  out <- switch(kind,
    reference_ipw = evaluate_reference_ipw_if_chunk(backend, rows, j, "frame", details = TRUE),
    ipw_legacy = , ipw_selection = evaluate_cp04_if_chunk(backend, rows, j, "frame", details = TRUE),
    dr_baseline = evaluate_dr_baseline_if_chunk(backend, rows, j, "frame", details = TRUE),
    dr_selection = evaluate_dr_corollary1_if_chunk(backend, rows, j, "frame", details = TRUE),
    or_baseline = evaluate_or_baseline_if_chunk(backend, rows, j, "frame", details = TRUE),
    or_legacy = , or_selection = evaluate_or_corollary1_if_chunk(backend, rows, j, "frame", details = TRUE),
    imputation = evaluate_imputation_if_chunk(backend, rows, j, "frame", details = TRUE),
    augmented = evaluate_augmented_if_chunk(backend, rows, j, "frame", details = TRUE),
    cpt_abort("No details-mode evaluator for this backend.", "cdid_conditional_pretest_integrity_error"))
  out$details
}

# Row chunk of a backend whose IF envelope is affine in the row count: the
# largest batch whose details-mode output and scratch fit the existing caps.
cpt_details_rows <- function(backend, kind, j, retained) {
  st <- get("cp04_reader", backend$reader_ref, inherits = FALSE)$state; ctrl <- backend$control
  B <- min(st$N, ctrl$max_chunk_rows)
  env_fn <- switch(kind, augmented = augmented_if_envelope, imputation = imputation_if_envelope,
    dr_selection = dr_corollary1_if_envelope, dr_baseline = dr_baseline_if_envelope,
    or_legacy = , or_selection = or_corollary1_if_envelope, NULL)
  if (is.null(env_fn)) return(B)
  one <- env_fn(backend, 1, j, TRUE); two <- env_fn(backend, 2, j, TRUE)
  slope <- two$output + two$scratch - one$output - one$scratch
  icpt <- one$output + one$scratch - slope
  oslope <- two$output - one$output; oicpt <- one$output - oslope
  own <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  if (slope > 0) B <- min(B, floor((ctrl$max_matrix_entries - retained - own - icpt) / slope))
  if (oslope > 0) B <- min(B, floor((ctrl$max_output_entries - oicpt) / oslope))
  if (B < 1) cdid_ipw_owner_abort("No admissible conditional-pretest details row batch within the existing caps.",
    "cdid_ipw_owner_bound_error")
  as.integer(B)
}

cpt_design_rows <- function(frame, st, ref, recipe, rows, cols = NULL) {
  X <- if (is.null(ref)) frame$Z_a[st$a_index[rows], , drop = FALSE] else chained_reference_x(ref, recipe, rows)
  if (is.null(cols)) X else X[, cols, drop = FALSE]
}

# ------------------------------------------------------------------------------
# CP02 risk-set IPW (A frame; rates and propensity of the auxiliary link)
# ------------------------------------------------------------------------------
cpt_extract_legacy_ipw <- function(state, pl, sel, chunk_rows) {
  frame <- state$frame_owner
  lv <- get("link_view", envir = pl$link_view_ref, inherits = FALSE)
  owner_ctrl <- state$control$owner_control
  nA <- as.integer(frame$n_a_units); N <- state$n_total; scale <- N / nA
  r <- as.double(pl$masses$rates); q <- as.double(pl$masses$q); m <- length(r)
  beta <- as.double(pl$coefficients); p <- length(beta)
  Abr <- matrix(as.double(pl$partials$A_beta_r), p, m)
  D_T <- as.double(pl$means$sum_a) / N; D_C <- as.double(pl$means$sum_b) / N
  spec <- list(id = sel$id, sigma = sel$sign, s = c(1, -1), D = c(D_T, D_C),
    coords = c(paste0("rate:", seq_len(m)), paste0("ps:", seq_len(p))))
  B <- as.integer(min(nA, chunk_rows))
  starts <- seq.int(1L, nA, by = B)
  pieces <- function(k) {
    a_s <- starts[k]; a_e <- min(nA, a_s + B - 1L)
    ch <- query_ipw_view_chunk(lv, a_s, a_e, owner_ctrl, read_outcomes = TRUE, include_design = TRUE)
    n <- a_e - a_s + 1L; rows <- state$a_to_full[seq.int(a_s, a_e)]
    w <- as.double(ch$weight); gi <- as.integer(ch$cohort_index); risk <- which(ch$in_risk)
    pos <- as.integer(ch$observed_risk$local_pos); gj <- as.integer(ch$observed_risk$cohort_index)
    S <- numeric(n); S[pos] <- 1
    xr <- matrix(0, n, m)
    if (length(risk)) xr[cbind(risk, gi[risk])] <- w[risk] * (S[risk] - r[gi[risk]]) / q[gi[risk]]
    score <- matrix(0, n, p)
    Z <- matrix(0, length(pos), p); odds <- act <- fw <- dy <- numeric(length(pos))
    lam_up <- matrix(0, length(pos), m)
    if (length(pos)) {
      Z <- ch$design$Z_obs
      lv_ <- ipw_logit_values(Z, beta)
      fw <- w[pos] / ipw_clip_floor(r)[gj]
      score[pos, ] <- (fw * (as.double(gj == 1L) - lv_$prob)) * Z
      odds <- lv_$odds
      act <- lv_$active
      dy <- as.double(ch$observed_risk$y_right_obs) - as.double(ch$observed_risk$y_left_obs)
      lam_up[cbind(seq_along(pos), gj)] <- ifelse(r[gj] > ipw_clip_eps, -1 / r[gj], 0)
    }
    xb <- t(cp04_solve(pl$factors$H_chol, t(score + xr %*% t(Abr))))
    xi <- scale * cbind(xr, xb)
    list(cpt_two_arm_piece(rows, pos, as.double(gj == 1L), fw, Z, odds, act, dy, lam_up, xi, D_T, D_C, m, p))
  }
  list(spec = spec, n_chunks = length(starts), pieces = pieces, estimate = as.double(pl$estimate),
    route = "CP02 risk-set IPW (A frame, observed-pair rates and one propensity logit per auxiliary link)")
}

# ------------------------------------------------------------------------------
# Full-roster IPW and DR families (reference baseline IPW, CP04 conditional_mtar
# IPW, baseline and conditional_mtar DR)
# ------------------------------------------------------------------------------
cpt_two_arm_rows <- function(kind, backend, tj, rows, reader, frame, ref) {
  st <- reader$state
  spec <- backend$treatment_specs[[tj]]; f <- backend$treatments[[tj]]
  ci <- st$cohort_index[rows]
  ci_in <- which(ci %in% spec$cohorts)
  if (identical(kind, "reference_ipw")) ci_in <- ci_in[st$a_index[rows[ci_in]] > 0L]
  if (identical(kind, "ipw_legacy")) ci_in <- ci_in[st$U[rows[ci_in]]]
  hit <- if (length(ci_in)) cp04_pair_rows(reader, rows[ci_in], spec$pair, TRUE) else list(S = numeric(0), dy = numeric(0))
  ob <- which(hit$S == 1)
  rr <- ci_in[ob]; dy <- as.double(hit$dy[ob]); n <- length(rr)
  kk <- match(ci[rr], spec$cohorts)
  rows_rr <- rows[rr]
  D <- as.double(kk == 1L)
  dr <- kind %in% c("dr_baseline", "dr_selection")
  trt_cols <- if (dr) backend$treatment_recipe$index else NULL
  fw <- numeric(n)
  if (kind %in% c("reference_ipw", "dr_baseline")) {
    up_n <- length(spec$cohorts)
    rate <- vapply(spec$support_ids, function(s) backend$supports[[s]]$rate, numeric(1L))
    lam_up <- matrix(0, n, up_n)
    if (n) lam_up[cbind(seq_len(n), kk)] <- ifelse(rate[kk] > ipw_clip_eps, -1 / rate[kk], 0)
    if (identical(kind, "reference_ipw")) {
      W <- if (n) chained_reference_weight(ref, spec$reference, rows_rr) else numeric(0)
      fw <- W / ipw_clip_floor(rate)[kk]
    } else {
      wr <- if (!n) numeric(0) else if (is.null(ref)) st$w[rows_rr] else chained_reference_weight(ref, spec$reference, rows_rr)
      v <- vapply(spec$support_ids, function(s) backend$supports[[s]]$v, numeric(1L))
      fw <- wr * v[kk]
    }
  } else if (identical(kind, "ipw_legacy")) {
    d <- ncol(st$Z); up_n <- length(spec$cohorts) * (d + 1L)
    lam_up <- matrix(0, n, up_n)
    for (k in seq_along(spec$cohorts)) {
      sub <- which(kk == k)
      if (!length(sub)) next
      sj <- spec$selection_ids[k]
      pp <- cp04_predict_selection(reader, rows_rr[sub], backend$selection_specs[[sj]], backend$selections[[sj]])
      fw[sub] <- st$w[rows_rr[sub]] * st$a[spec$cohorts[k]] / (st$rho[rows_rr[sub]] * ipw_clip_floor(pp))
      lam_up[sub, (k - 1L) * (d + 1L) + seq_len(d + 1L)] <- cbind(1 / st$a[spec$cohorts[k]],
        -((1 - pp) * ipw_clip_floor_active(pp)) * st$Z[rows_rr[sub], , drop = FALSE])
    }
  } else {
    sel_fit <- backend$selection_backend; up_map <- backend$up_map
    ids <- f$upstream_ids; up_n <- length(ids)
    lam_up <- matrix(0, n, up_n)
    for (k in seq_along(spec$cohorts)) {
      sub <- which(kk == k)
      if (!length(sub)) next
      qk <- spec$selection_ids[k]
      ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, qk, rows_rr[sub])
      fw[sub] <- ld$pred$weight
      lam_up[sub, match(up_map$pair_channels[[qk]], ids)] <- ld$log_deriv
    }
  }
  X <- if (n) cpt_design_rows(frame, st, ref, spec$reference, rows_rr, trt_cols) else
    matrix(0, 0L, length(f$coefficients))
  lv <- if (n) ipw_logit_values(X, f$coefficients) else list(odds = numeric(0), active = numeric(0))
  out <- list(rr = rr, D = D, fw = fw, X = X, odds = lv$odds, active = lv$active, dy = dy,
    lam_up = lam_up, up_n = up_n, p_trt = length(f$coefficients))
  if (dr) {
    reg <- backend$regressions[[spec$regression]]; p_or <- backend$shared$p
    kept <- if (p_or > 0L) backend$or_frame_index[backend$outcome_recipe$slope_index[backend$shared$keep]] else integer(0)
    Z_or <- if (n && p_or > 0L) cpt_design_rows(frame, st, ref, spec$reference, rows_rr, kept) else matrix(0, n, p_or)
    d_or <- reg$L * Z_or
    out$d_or <- d_or
    out$rho <- dy - (if (backend$outcome_recipe$intercept) reg$alpha else 0) -
      (if (p_or > 0L) as.vector(d_or %*% backend$shared$beta) else 0)
  } else out$rho <- dy
  out
}

cpt_extract_two_arm <- function(state, backend, kind, j, sel, chunk_rows) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  ref <- chained_backend_reference(backend)
  tj <- backend$links[[j]]$model
  spec_t <- backend$treatment_specs[[tj]]; f <- backend$treatments[[tj]]
  dr <- kind %in% c("dr_baseline", "dr_selection")
  up_n <- switch(kind, reference_ipw = , dr_baseline = length(spec_t$cohorts),
    ipw_legacy = length(spec_t$cohorts) * (ncol(st$Z) + 1L), length(f$upstream_ids))
  p_trt <- length(f$coefficients)
  icpt <- dr && isTRUE(backend$outcome_recipe$intercept)
  p_or <- if (dr) backend$shared$p else 0L
  coords <- c(paste0("up:", seq_len(up_n)), paste0("ps:", seq_len(p_trt)),
    if (icpt) paste0("or_alpha:", spec_t$regression), if (p_or) paste0("or_beta:", seq_len(p_or)))
  D_T <- unname(f$means["dp"]); D_M <- unname(f$means["dm"])
  spec <- list(id = sel$id, sigma = sel$sign, s = c(1, -1), D = c(D_T, D_M), coords = coords)
  B <- as.integer(min(st$N, chunk_rows, cpt_details_rows(backend, kind, j, state$retained_now)))
  starts <- seq.int(1L, st$N, by = B)
  pieces <- function(k) {
    rows <- seq.int(starts[k], min(st$N, starts[k] + B - 1L))
    det <- cpt_backend_details(backend, kind, rows, j)
    md <- det$models[[tj]]
    xi_ps <- if (!is.null(md$gamma)) md$gamma else md$beta
    xi <- cbind(md$upstream, xi_ps)
    if (dr) {
      if (icpt) xi <- cbind(xi, det$alpha[, match(spec_t$regression, det$regression_ids)])
      if (p_or) xi <- cbind(xi, det$beta)
    }
    xi <- unname(xi)
    rd <- cpt_two_arm_rows(kind, backend, tj, rows, reader, frame, ref)
    list(cpt_two_arm_piece(rows, rd$rr, rd$D, rd$fw, rd$X, rd$odds, rd$active, rd$rho, rd$lam_up, xi, D_T, D_M,
      up_n, p_trt, d_or = if (dr) rd$d_or else NULL, or_icpt = icpt))
  }
  list(spec = spec, n_chunks = length(starts), pieces = pieces,
    estimate = chained_full_link_estimate(backend, j), route = backend$law)
}

# ------------------------------------------------------------------------------
# Control-only shared-slope OR (baseline, conditional_mtar legacy and models):
# one treated residual ratio with fitted weights q.
# ------------------------------------------------------------------------------
cpt_extract_or <- function(state, backend, kind, j, sel, chunk_rows) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  ref <- chained_backend_reference(backend)
  l <- backend$links[[j]]; s <- l$support_id; b <- l$regression
  sp <- backend$support_specs[[s]]; u <- backend$supports[[s]]
  blk <- backend$regressions[[b]]
  icpt <- isTRUE(backend$outcome_recipe$intercept); p <- backend$shared$p
  kept <- backend$outcome_recipe$slope_index[backend$shared$keep]
  channels <- if (identical(kind, "or_baseline")) s else u$channels
  coords <- c(paste0("up:", channels), if (icpt) paste0("or_alpha:", b), if (p) paste0("or_beta:", seq_len(p)))
  spec <- list(id = sel$id, sigma = sel$sign, s = 1, D = l$d, coords = coords)
  Dl <- length(coords); n_ch <- length(channels)
  B <- as.integer(min(st$N, chunk_rows, cpt_details_rows(backend, kind, j, state$retained_now)))
  starts <- seq.int(1L, st$N, by = B)
  pieces <- function(k) {
    rows <- seq.int(starts[k], min(st$N, starts[k] + B - 1L)); n <- length(rows)
    det <- cpt_backend_details(backend, kind, rows, j)
    xi_ch <- if (identical(kind, "or_baseline")) det$rates[, match(s, det$support_ids), drop = FALSE] else
      det$upstream[, match(channels, det$upstream_ids), drop = FALSE]
    xi <- xi_ch
    if (icpt) xi <- cbind(xi, det$alpha[, match(b, det$regression_ids)])
    if (p) xi <- cbind(xi, det$beta)
    xi <- unname(xi)
    ci <- st$cohort_index[rows]
    ix <- which(ci == sp$cohort)
    if (!identical(kind, "or_baseline")) ix <- ix[st$U[rows[ix]]]
    hit <- if (length(ix)) cp04_pair_rows(reader, rows[ix], sp$pair, TRUE) else list(S = numeric(0), dy = numeric(0))
    ob <- which(hit$S == 1); rr <- ix[ob]; rows_rr <- rows[rr]; dy <- as.double(hit$dy[ob]); m <- length(rr)
    nu <- rh <- matrix(0, n, 1L); lam <- pim <- matrix(0, n, Dl)
    if (m) {
      z <- if (is.null(ref)) frame$Z_a[st$a_index[rows_rr], kept, drop = FALSE] else
        chained_reference_design(ref, sp$reference, rows_rr)[, kept, drop = FALSE]
      if (identical(kind, "or_baseline")) {
        w <- if (is.null(ref)) st$w[rows_rr] else chained_reference_weight(ref, sp$reference, rows_rr)
        q <- w * u$v
        lam[rr, 1L] <- if (u$rate > ipw_clip_eps) -1 / u$rate else 0
      } else if (identical(kind, "or_legacy")) {
        pr <- cp04_predict_selection(reader, rows_rr, backend$selection_specs[[s]], backend$selections[[s]])
        q <- st$w[rows_rr] * or_corollary1_factor(st$a[sp$cohort], st$rho[rows_rr], pr)
        lam[rr, seq_len(n_ch)] <- cbind(1 / st$a[sp$cohort],
          -((1 - pr) * ipw_clip_floor_active(pr)) * st$Z[rows_rr, , drop = FALSE])
      } else {
        ld <- selection_pair_cohort_log_deriv(backend$selection_backend, backend$up_map, backend$q_of_s[s], rows_rr)
        a <- if (is.null(ref)) st$a[sp$cohort] else chained_reference_stat(ref, sp$reference, sp$cohort, "a")
        w <- if (is.null(ref)) st$w[rows_rr] else chained_reference_weight(ref, sp$reference, rows_rr)
        q <- w * or_corollary1_factor(a, st$rho[rows_rr], ld$pred$prob)
        lam[rr, seq_len(n_ch)] <- ld$log_deriv
      }
      d <- blk$L * z
      nu[rr, 1L] <- q
      rh[rr, 1L] <- dy - blk$alpha - (if (p) as.vector(d %*% backend$shared$beta) else 0)
      off <- n_ch
      if (icpt) {
        pim[rr, off + 1L] <- -1
        off <- off + 1L
      }
      if (p) pim[rr, off + seq_len(p)] <- -d
    }
    list(cpt_piece(rows, nu, rh, list(lam), list(pim), 1 / l$d, xi))
  }
  list(spec = spec, n_chunks = length(starts), pieces = pieces, estimate = l$mu, route = backend$law)
}

# ------------------------------------------------------------------------------
# Pure imputation (ipw / reg / dr) and the augmented correction (dr): treated
# T and B_g, control C and B_o ratios with imputed (or augmented) responses.
# ------------------------------------------------------------------------------
cpt_extract_imputation <- function(state, backend, kind, j, sel, chunk_rows) {
  augmented <- identical(kind, "augmented")
  if (!is.null(backend$pair_outcomes))
    cpt_abort("Retained private pair outcomes are not a public conditional-pretest input.", "cdid_conditional_pretest_integrity_error")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  lay <- backend$layout
  # The augmented layout fit owns a sealed layout reader; pure imputation keeps
  # the original reader and its layout state.
  st <- if (is.null(lay) || augmented) reader$state else lay$state
  N <- st$N
  X_F <- if (!is.null(lay)) backend$X_F else if (augmented) augmented_backend_X_F(backend) else imputation_backend_X_F(backend, reader)
  pair_reader <- if (!is.null(lay) && !augmented) function(r_sub, pair, read_outcomes = TRUE)
    cp04_pair_rows(reader, lay$orig_row[r_sub], pair, read_outcomes) else
    function(r_sub, pair, read_outcomes = TRUE) cp04_pair_rows(reader, r_sub, pair, read_outcomes)
  recipe <- backend$outcome_recipe; trt_recipe <- backend$treatment_recipe
  est <- if (augmented) "dr" else backend$est_method
  lf <- backend$links[[j]]
  cells <- backend$cells; groups <- backend$groups; sh <- backend$shared_or
  ci_g <- lf$target_cohort_index; cid_g <- lf$target_cell_id
  k_ctrl <- if (est %in% c("ipw", "dr")) lf$control_cohort_indices else integer(0)
  c_ctrl <- if (est %in% c("ipw", "dr")) lf$control_cell_ids else integer(0)
  use_or <- est %in% c("reg", "dr")
  icpt <- isTRUE(recipe$intercept)
  # Coordinates.
  a_k <- unique(c(ci_g, if (est %in% c("ipw", "dr")) lf$control_cohort_indices))
  cids <- unique(c(cid_g, c_ctrl))
  gids <- unique(vapply(cells[cids], `[[`, integer(1L), "group_id"))
  gids <- gids[vapply(groups[gids], function(gf) gf$p > 0L, logical(1L))]
  tid <- if (est %in% c("ipw", "dr")) lf$trt_id else NA_integer_
  p_trt <- if (is.na(tid)) 0L else trt_recipe$p
  bid <- if (use_or) lf$or_block_id else NA_integer_
  p_or <- if (use_or) sh$p else 0L
  th <- if (augmented) sort(unique(unlist(lapply(cells[cids], `[[`, "theta_ids")))) else integer(0)
  coords <- c(paste0("a:", a_k), if (icpt) paste0("ia:", cids),
    unlist(lapply(gids, function(g) paste0("ib:", g, ":", seq_len(groups[[g]]$p)))),
    if (p_trt) paste0("ps:", seq_len(p_trt)), if (use_or && icpt) paste0("oa:", bid),
    if (p_or) paste0("ob:", seq_len(p_or)), if (length(th)) paste0("th:", th))
  Dl <- length(coords)
  col_a <- function(k) match(paste0("a:", k), coords)
  col_ia <- function(cid) match(paste0("ia:", cid), coords)
  col_ib <- function(g) match(paste0("ib:", g, ":", seq_len(groups[[g]]$p)), coords)
  col_ps <- match(paste0("ps:", seq_len(p_trt)), coords)
  col_oa <- if (use_or && icpt) match(paste0("oa:", bid), coords) else integer(0)
  col_ob <- match(paste0("ob:", seq_len(p_or)), coords)
  col_th <- function(ids) match(paste0("th:", ids), coords)
  terms <- c("T", if (use_or) "Bg", if (est %in% c("ipw", "dr")) "C", if (est == "dr") "Bo")
  s <- c(T = 1, Bg = -1, C = -1, Bo = 1)[terms]
  Dv <- c(T = lf$d_plus, Bg = lf$d_plus, C = lf$d_minus, Bo = lf$d_minus)[terms]
  spec <- list(id = sel$id, sigma = sel$sign, s = unname(s), D = unname(Dv), coords = coords)
  R <- length(terms); sD <- unname(s / Dv)
  sel_fit <- backend$selection_backend; up_map <- backend$up_map

  cell_resp <- function(cid, lr) {
    cf <- cells[[cid]]; gf <- groups[[cf$group_id]]; n <- length(lr)
    d <- if (gf$p > 0L) lf$L * X_F[lr, gf$kept_cols_F, drop = FALSE] else matrix(0, n, 0L)
    if (!augmented) {
      jh <- cf$alpha + (if (gf$p > 0L) as.vector(d %*% gf$beta) else 0)
      return(list(rho = jh, fac = rep(1, n), d = d, pt = NULL))
    }
    ev <- augmented_cell_rows(sel_fit, up_map, X_F, pair_reader, cf, gf, lr, need_logderiv = TRUE)
    fac <- rep(1, n)
    pt <- matrix(0, n, length(cf$theta_ids))
    if (length(ev$obs)) {
      fac[ev$obs] <- 1 - ev$pi
      if (length(cf$theta_ids)) pt[ev$obs, ] <- (ev$e * ev$pi) * ev$T_logit
    }
    list(rho = ev$Z, fac = fac, d = d, pt = pt)
  }
  or_resp <- function(lr) {
    z <- if (p_or > 0L) lf$L * X_F[lr, sh$kept_cols_F, drop = FALSE] else matrix(0, length(lr), 0L)
    bf <- sh$blocks[[bid]]
    list(rho = (if (icpt) bf$alpha else 0) + (if (p_or > 0L) as.vector(z %*% sh$beta) else 0), d = z)
  }
  put_resp <- function(pim_r, ii, cid, rs) {
    cf <- cells[[cid]]
    if (icpt) pim_r[ii, col_ia(cid)] <- rs$fac
    gid <- cf$group_id
    if (gid %in% gids) pim_r[ii, col_ib(gid)] <- rs$fac * rs$d
    if (augmented && length(cf$theta_ids)) pim_r[ii, col_th(cf$theta_ids)] <- rs$pt
    pim_r
  }
  rows_piece <- function(lr, det, orig) {
    n <- length(lr); ci <- st$cohort_index[lr]; w <- st$w[lr]
    nu <- rho <- matrix(0, n, R, dimnames = list(NULL, terms))
    lam <- lapply(seq_len(R), function(r) matrix(0, n, Dl)); pim <- lam
    names(lam) <- names(pim) <- terms
    ig <- which(ci == ci_g)
    if (length(ig)) {
      a_g <- st$a[ci_g]; nug <- w[ig] * a_g
      rs <- cell_resp(cid_g, lr[ig])
      nu[ig, "T"] <- nug; rho[ig, "T"] <- rs$rho
      lam$T[ig, col_a(ci_g)] <- 1 / a_g
      pim$T <- put_resp(pim$T, ig, cid_g, rs)
      if (use_or) {
        orr <- or_resp(lr[ig])
        nu[ig, "Bg"] <- nug; rho[ig, "Bg"] <- orr$rho
        lam$Bg[ig, col_a(ci_g)] <- 1 / a_g
        if (icpt) pim$Bg[ig, col_oa] <- 1
        if (p_or) pim$Bg[ig, col_ob] <- orr$d
      }
    }
    if (est %in% c("ipw", "dr")) {
      tf <- backend$treatments[[tid]]
      for (idx in seq_along(lf$control_cohort_indices)) {
        kc <- lf$control_cohort_indices[idx]; cidk <- lf$control_cell_ids[idx]
        ik <- which(ci == kc)
        if (!length(ik)) next
        Xt <- X_F[lr[ik], trt_recipe$index, drop = FALSE]
        lvt <- ipw_logit_values(Xt, tf$coefficients)
        odds <- lvt$odds
        a_kc <- st$a[kc]; nuk <- w[ik] * a_kc * odds
        rs <- cell_resp(cidk, lr[ik])
        nu[ik, "C"] <- nuk; rho[ik, "C"] <- rs$rho
        lam$C[ik, col_a(kc)] <- 1 / a_kc
        lam$C[ik, col_ps] <- lvt$active * Xt
        pim$C <- put_resp(pim$C, ik, cidk, rs)
        if (est == "dr") {
          orr <- or_resp(lr[ik])
          nu[ik, "Bo"] <- nuk; rho[ik, "Bo"] <- orr$rho
          lam$Bo[ik, col_a(kc)] <- 1 / a_kc
          lam$Bo[ik, col_ps] <- lvt$active * Xt
          if (icpt) pim$Bo[ik, col_oa] <- 1
          if (p_or) pim$Bo[ik, col_ob] <- orr$d
        }
      }
    }
    xi <- matrix(0, n, Dl)
    xi[, col_a(a_k)] <- det$a[, a_k, drop = FALSE]
    if (icpt) xi[, col_ia(cids)] <- det$imputation_alpha[, cids, drop = FALSE]
    for (g in gids) xi[, col_ib(g)] <- det$imputation_beta[[g]]
    if (p_trt) xi[, col_ps] <- det$trt_if[[tid]]$gamma
    if (use_or && icpt) {
      abids <- if (!is.null(det$or_if$alpha_block_ids)) det$or_if$alpha_block_ids else seq_len(sh$n_blocks)
      xi[, col_oa] <- det$or_if$alpha[, match(bid, abids)]
    }
    if (p_or) xi[, col_ob] <- det$or_if$beta
    if (length(th)) xi[, col_th(th)] <- det$selection_theta[, match(th, det$selection_theta_ids), drop = FALSE]
    cpt_piece(orig, nu, rho, lam, pim, sD, xi)
  }
  B <- as.integer(min(N, chunk_rows, cpt_details_rows(backend, kind, j, state$retained_now)))
  starts <- seq.int(1L, N, by = B)
  pieces <- function(k) {
    rows <- seq.int(starts[k], min(N, starts[k] + B - 1L))
    if (is.null(lay)) return(list(rows_piece(rows, cpt_backend_details(backend, kind, rows, j), rows)))
    # D60 recipe-cohort layout: the unchanged core runs on the layout rows of
    # these original rows in details mode; rows keep their original identity.
    map <- chained_reference_layout_rows(lay, rows)
    if (!length(map$rows)) return(list())
    pre <- if (augmented) augmented_if_budget(backend, length(rows), j, TRUE) else imputation_if_budget(backend, length(rows), j, TRUE)
    out <- list(); BL <- backend$control$max_chunk_rows
    for (s0 in seq.int(1L, length(map$rows), by = BL)) {
      lr <- as.integer(map$rows[seq.int(s0, min(length(map$rows), s0 + BL - 1L))])
      det <- if (augmented) augmented_evaluate_if_chunk_core(st = st, X_F = X_F, pair_reader = pair_reader,
          backend = backend, rows = lr, links = as.integer(j), scale = "frame", details = TRUE, sets = pre)$details else
        imputation_evaluate_if_chunk_core(st = st, X_F = X_F, pair_reader = pair_reader, core_fit = backend,
          recipe = recipe, or_frame_index = backend$or_frame_index, trt_recipe = trt_recipe, est_method = backend$est_method,
          rows = lr, links = as.integer(j), scale = "frame", details = TRUE, sets = pre)$details
      out[[length(out) + 1L]] <- rows_piece(lr, det, lay$orig_row[lr])
    }
    out
  }
  list(spec = spec, n_chunks = length(starts), pieces = pieces, estimate = lf$mu, route = backend$law)
}

# Extractor of the auxiliary descriptor q of one declared cell.
cpt_cell_extractor <- function(state, q, sel, chunk_rows) {
  if (is.null(state$cp04_reader)) {
    pl <- state$auxiliary$point_links[[q]]
    if (is.null(pl)) cpt_abort("Missing auxiliary point link.", "cdid_conditional_pretest_integrity_error")
    return(cpt_extract_legacy_ipw(state, pl, sel, chunk_rows))
  }
  reg <- state$registry
  j <- reg$cp04_aux_index[q]
  backend <- reg$cp04
  kind <- cpt_backend_kind(backend)
  switch(kind,
    reference_ipw = , ipw_legacy = , ipw_selection = , dr_baseline = , dr_selection =
      cpt_extract_two_arm(state, backend, kind, j, sel, chunk_rows),
    or_baseline = , or_legacy = , or_selection = cpt_extract_or(state, backend, kind, j, sel, chunk_rows),
    imputation = , augmented = cpt_extract_imputation(state, backend, kind, j, sel, chunk_rows),
    cpt_abort(paste0("The conditional pretest has no process evaluator for this backend (", kind, ")."),
      "cdid_conditional_pretest_unsupported_error"))
}
