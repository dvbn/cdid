# Private and public selection-model consumers for chained IPW and shared-slope OR
# (D48 Step 2 Section B; CONNECTION_BRIEF.md and DATE_SELECTION_SPEC.md).
#
# Connects the Section A observation-probability backend (R/chained_selection_models.R)
# to:
#   1. Chained IPW (fit_cp04_selection_models, evaluate_cp04_selection_if_chunk)
#   2. Shared one-period outcome-change OR (fit_or_selection_models, evaluate_or_selection_if_chunk)
#   3. Public request normalization and selection/nuisance diagnostics

normalize_selection_declaration <- function(selection = "none",
                                            selection_model = c("pair_logit", "independent_period_logit"),
                                            selection_pooling = c("separate", "within_cohort", "across_cohorts"),
                                            selection_xformla = NULL,
                                            selection_universe = "F",
                                            selection_rho = NULL,
                                            selection_justification = NULL,
                                            is_custom_est = FALSE) {
  if (identical(selection, c("none", "conditional_mtar"))) {
    selection <- "none"
  }
  model_default <- c("pair_logit", "independent_period_logit")
  pooling_default <- c("separate", "within_cohort", "across_cohorts")
  if (identical(selection_model, model_default) || is.null(selection_model)) {
    selection_model <- "pair_logit"
  }
  if (identical(selection_pooling, pooling_default) || is.null(selection_pooling)) {
    selection_pooling <- "separate"
  }

  chained_request_need(
    chained_request_chr1(selection) && selection %in% c("none", "conditional_mtar", "corollary1"),
    '`selection` must be "none", "conditional_mtar", or legacy alias "corollary1".',
    "cdid_chained_request_unsupported_selection"
  )
  chained_request_need(
    chained_request_chr1(selection_model) && selection_model %in% model_default,
    '`selection_model` must be "pair_logit" or "independent_period_logit".',
    "cdid_chained_request_unsupported_selection"
  )
  chained_request_need(
    chained_request_chr1(selection_pooling) && selection_pooling %in% pooling_default,
    '`selection_pooling` must be "separate", "within_cohort", or "across_cohorts".',
    "cdid_chained_request_unsupported_selection"
  )

  if (identical(selection, "none")) {
    chained_request_need(
      is.null(selection_xformla) &&
        identical(selection_universe, "F") &&
        is.null(selection_rho) &&
        is.null(selection_justification) &&
        identical(selection_model, "pair_logit") &&
        identical(selection_pooling, "separate"),
      'Selection formula/model/pooling/universe/transport declarations require `selection = "conditional_mtar"` (or legacy `"corollary1"`).',
      "cdid_chained_request_unsupported_selection"
    )
    return(list(
      enabled = FALSE,
      requested = "none",
      canonical_assumption = "none",
      legacy_alias = FALSE,
      selection_model = NULL,
      selection_pooling = NULL
    ))
  }

  # O15 (D64): custom IPW callbacks are accepted under selection_est_method =
  # "ipw" (the only caller passing is_custom_est here besides the augmented
  # route, which refuses callbacks first); each callback must return its
  # explicit `log_weight_gradient`, checked when the callback is called.

  if (identical(selection, "corollary1")) {
    chained_request_need(
      identical(selection_model, "pair_logit"),
      '`selection = "corollary1"` is the legacy alias for `selection = "conditional_mtar", selection_model = "pair_logit"` and cannot be combined with `selection_model = "independent_period_logit"`.',
      c("cdid_chained_request_unsupported_selection", "cdid_cp04_unsupported_error")
    )
  }

  if (identical(selection_model, "independent_period_logit")) {
    chained_request_need(
      chained_request_chr1(selection_universe) && selection_universe %in% c("F", "B", "A"),
      '`selection_universe` must be "F", "B" (at least once), or "A" (repeat).'
    )
    chained_request_need(
      identical(selection_universe, "F") && is.null(selection_rho) && is.null(selection_justification),
      '`selection_model = "independent_period_logit"` is defined only on the representative frame `selection_universe = "F"` with `selection_rho = NULL` and `selection_justification = NULL`; conditioning on A or B breaks date factorization.',
      c("cdid_chained_request_unsupported_selection", "cdid_cp04_unsupported_error")
    )
  } else {
    cp04_validate_declaration(selection_universe, selection_rho, selection_justification)
  }

  list(
    enabled = TRUE,
    requested = selection,
    canonical_assumption = "conditional_mtar",
    legacy_alias = identical(selection, "corollary1"),
    selection_model = selection_model,
    selection_pooling = selection_pooling
  )
}

selection_use_legacy_consumer <- function(reader, selection_model = "pair_logit", selection_pooling = "separate") {
  identical(selection_model, "pair_logit") &&
    identical(selection_pooling, "separate") &&
    is.matrix(reader$state$Z) &&
    ncol(reader$state$Z) >= 1L
}

# Build the global upstream coordinate layout u = ((a_k^U)_{k in K_sel}, (theta_m)_{m in M})
# and per-pair-cohort channel maps for a fitted cdid_selection_models object.
# D60: a recipe-indexed fit has one repeat-rate coordinate a_k^{g,U} per
# (cohort, recipe) of its pair-cohorts (`rates_ref` gives each one's recipe).
build_selection_upstream_map <- function(sel_fit, cohorts_a = NULL) {
  if (!inherits(sel_fit, "cdid_selection_models_table")) validate_selection_models(sel_fit)
  pc_ref <- if (is.null(sel_fit$reference_ref)) NULL else vapply(sel_fit$pair_cohorts, `[[`, integer(1), "reference")
  rates_ref <- NULL
  if (!is.null(pc_ref)) {
    kr <- unique(data.frame(k = vapply(sel_fit$pair_cohorts, `[[`, integer(1), "cohort"), r = pc_ref))
    kr <- kr[order(kr$k, kr$r, method = "radix"), , drop = FALSE]
    cohorts_a <- kr$k; rates_ref <- as.integer(kr$r)
  } else if (is.null(cohorts_a)) {
    cohorts_a <- sort(unique(vapply(sel_fit$pair_cohorts, `[[`, integer(1), "cohort")))
  }
  cohorts_a <- as.integer(cohorts_a)
  n_a <- length(cohorts_a)
  M <- length(sel_fit$blocks)
  alpha_coords <- vector("list", M)
  beta_coords <- vector("list", M)
  theta_coords <- vector("list", M)

  pos <- n_a
  coh_vec <- cohorts_a
  blk_vec <- rep.int(0L, n_a)
  cell_vec <- rep.int(0L, n_a)
  col_vec <- rep.int(0L, n_a)
  role_vec <- rep.int("rate", n_a)

  for (m in seq_len(M)) {
    bl <- sel_fit$blocks[[m]]
    if (identical(bl$branch, "zero_col")) {
      alpha_coords[[m]] <- integer(0)
      beta_coords[[m]] <- integer(0)
      theta_coords[[m]] <- integer(0)
      next
    }
    n_al <- if (bl$has_intercept) bl$B_m else 0L
    n_be <- bl$d_s
    al_c <- if (n_al > 0L) seq.int(pos + 1L, pos + n_al) else integer(0)
    pos <- pos + n_al
    be_c <- if (n_be > 0L) seq.int(pos + 1L, pos + n_be) else integer(0)
    pos <- pos + n_be
    alpha_coords[[m]] <- al_c
    beta_coords[[m]] <- be_c
    theta_coords[[m]] <- c(al_c, be_c)

    if (n_al > 0L) {
      cell_cohs <- vapply(sel_fit$cells[bl$cell_ids], `[[`, integer(1), "cohort")
      coh_vec <- c(coh_vec, cell_cohs)
      blk_vec <- c(blk_vec, rep.int(m, n_al))
      cell_vec <- c(cell_vec, as.integer(bl$cell_ids))
      col_vec <- c(col_vec, seq_len(n_al))
      role_vec <- c(role_vec, rep.int("alpha", n_al))
    }
    if (n_be > 0L) {
      b_coh <- if (identical(sel_fit$selection_pooling, "across_cohorts")) 0L else
        as.integer(sel_fit$block_specs[[m]]$cohort)
      coh_vec <- c(coh_vec, rep.int(b_coh, n_be))
      blk_vec <- c(blk_vec, rep.int(m, n_be))
      cell_vec <- c(cell_vec, rep.int(0L, n_be))
      col_vec <- c(col_vec, seq_len(n_be))
      role_vec <- c(role_vec, rep.int("beta", n_be))
    }
  }

  D <- pos
  n_pc <- length(sel_fit$pair_cohorts)
  pair_channels <- vector("list", n_pc)
  for (q in seq_len(n_pc)) {
    pc <- sel_fit$pair_cohorts[[q]]
    ka <- if (is.null(rates_ref)) match(pc$cohort, cohorts_a) else which(cohorts_a == pc$cohort & rates_ref == pc$reference)
    if (identical(sel_fit$selection_model, "pair_logit")) {
      m <- pc$block_id
      b <- pc$local_b
      bl <- sel_fit$blocks[[m]]
      kl <- c(
        if (!identical(bl$branch, "zero_col") && bl$has_intercept) alpha_coords[[m]][b] else integer(0),
        beta_coords[[m]]
      )
    } else {
      ms <- pc$s_block_id; bs <- pc$s_local_b; bl_s <- sel_fit$blocks[[ms]]
      mt <- pc$t_block_id; bt <- pc$t_local_b; bl_t <- sel_fit$blocks[[mt]]
      if (ms != mt) {
        kl <- c(
          if (!identical(bl_s$branch, "zero_col") && bl_s$has_intercept) alpha_coords[[ms]][bs] else integer(0),
          beta_coords[[ms]],
          if (!identical(bl_t$branch, "zero_col") && bl_t$has_intercept) alpha_coords[[mt]][bt] else integer(0),
          beta_coords[[mt]]
        )
      } else {
        m <- ms; bl <- bl_s
        kl <- c(
          if (!identical(bl$branch, "zero_col") && bl$has_intercept) c(alpha_coords[[m]][bs], alpha_coords[[m]][bt]) else integer(0),
          beta_coords[[m]]
        )
      }
    }
    pair_channels[[q]] <- as.integer(c(ka, kl))
  }

  max_dq <- max(c(0L, vapply(pair_channels, function(ch) length(ch) - 1L, integer(1))))
  list(
    n_a = n_a,
    d = max_dq,
    D = D,
    cohorts_a = cohorts_a,
    rates_ref = rates_ref,
    alpha_coords = alpha_coords,
    beta_coords = beta_coords,
    theta_coords = theta_coords,
    pair_channels = pair_channels,
    cohort = as.integer(coh_vec),
    support = as.integer(blk_vec),
    cell = as.integer(cell_vec),
    column = as.integer(col_vec),
    role = role_vec
  )
}

# Evaluates predicted pair probabilities, weights w_i * a_k^U / (rho_i * max(p_{e,k}(X_i), eps)),
# and the log-weight derivative matrix with respect to pair_channels[[q]] on `rows`.
# The selection columns carry the final pair mask 1{p_{e,k} > eps} (for date
# products, on both endpoint contributions); the a_k^U column is unmasked.
selection_pair_cohort_log_deriv <- function(sel_fit, up_map, q, rows, pred = NULL, Z = NULL, a_vec = NULL) {
  pc <- sel_fit$pair_cohorts[[q]]
  ref <- NULL
  if (is.null(Z) || is.null(a_vec)) {
    ref <- chained_backend_reference(sel_fit)
    if (is.null(ref)) {
      reader <- get("cp04_reader", sel_fit$reader_ref, inherits = FALSE)
      st <- reader$state
      Z <- st$Z
      a_vec <- st$a
    }
  }
  a_pc <- if (is.null(ref)) a_vec[pc$cohort] else chained_reference_stat(ref, pc$reference, pc$cohort, "a")
  ch <- up_map$pair_channels[[q]]
  dq <- length(ch) - 1L
  n_r <- length(rows)

  if (!n_r) {
    return(list(
      pred = list(prob = numeric(0), weight = numeric(0), a = a_pc),
      T_logit = matrix(0, 0L, dq),
      log_deriv = matrix(0, 0L, length(ch))
    ))
  }

  if (is.null(pred)) {
    pred <- predict_pair_selection(sel_fit, rows, q)
  }
  act <- ipw_clip_floor_active(pred$prob)
  d_s <- sel_fit$design_spec$d_s
  Z_slope <- if (d_s > 0L) {
    if (is.null(ref)) Z[rows, sel_fit$design_spec$slope_cols, drop = FALSE] else
      chained_reference_s(ref, pc$reference, rows)[, sel_fit$design_spec$slope_cols, drop = FALSE]
  } else matrix(0, n_r, 0L)

  if (identical(sel_fit$selection_model, "pair_logit")) {
    m <- pc$block_id
    bl <- sel_fit$blocks[[m]]
    if (identical(bl$branch, "zero_col")) {
      T_logit <- matrix(0, n_r, 0L)
    } else {
      w_p <- -((1 - pred$prob) * act)
      T_logit <- cbind(
        if (bl$has_intercept) matrix(w_p, n_r, 1L) else matrix(0, n_r, 0L),
        if (d_s > 0L) w_p * Z_slope else matrix(0, n_r, 0L)
      )
    }
  } else {
    ms <- pc$s_block_id; bl_s <- sel_fit$blocks[[ms]]
    mt <- pc$t_block_id; bl_t <- sel_fit$blocks[[mt]]
    if (identical(bl_s$branch, "zero_col")) {
      T_logit <- matrix(0, n_r, 0L)
    } else {
      w_s <- -((1 - pred$prob_s) * act)
      w_t <- -((1 - pred$prob_t) * act)
      if (ms != mt) {
        T_s <- cbind(
          if (bl_s$has_intercept) matrix(w_s, n_r, 1L) else matrix(0, n_r, 0L),
          if (d_s > 0L) w_s * Z_slope else matrix(0, n_r, 0L)
        )
        T_t <- cbind(
          if (bl_t$has_intercept) matrix(w_t, n_r, 1L) else matrix(0, n_r, 0L),
          if (d_s > 0L) w_t * Z_slope else matrix(0, n_r, 0L)
        )
        T_logit <- cbind(T_s, T_t)
      } else {
        T_logit <- cbind(
          if (bl_s$has_intercept) cbind(w_s, w_t) else matrix(0, n_r, 0L),
          if (d_s > 0L) (w_s + w_t) * Z_slope else matrix(0, n_r, 0L)
        )
      }
    }
  }

  log_deriv <- cbind(rep.int(1 / a_pc, n_r), T_logit)
  list(pred = pred, T_logit = T_logit, log_deriv = log_deriv)
}

# Evaluates upstream IF columns `uids` on original-unit row chunk `rows`.
# Reused date/pair cell scores are summed by original unit i within each block
# inside evaluate_selection_block_if before Fisher-information inversion.
evaluate_selection_upstream_chunk <- function(sel_fit, up_map, rows, uids, details = FALSE) {
  reader <- get("cp04_reader", sel_fit$reader_ref, inherits = FALSE)
  st <- reader$state
  ref <- chained_backend_reference(sel_fit)
  B <- length(rows)
  ci <- st$cohort_index[rows]
  # D60: no common weight; each rate coordinate and cell reads its recipe's W_g.
  w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
  U_rows <- st$U[rows]
  A <- as.double(st$a_index[rows] > 0L)
  w_of <- function(r, ix) if (is.null(ref)) w[ix] else chained_reference_weight(ref, r, rows[ix])

  up <- matrix(0, B, length(uids))
  for (j in which(uids <= up_map$n_a)) {
    k <- up_map$cohorts_a[uids[j]]
    r <- up_map$rates_ref[uids[j]]
    ix <- which(ci == k & U_rows)
    if (length(ix)) {
      a_k <- if (is.null(ref)) st$a[k] else chained_reference_stat(ref, r, k, "a")
      m_k <- if (is.null(ref)) st$mass_U[k] else chained_reference_stat(ref, r, k, "mass_U")
      up[ix, j] <- w_of(r, ix) * (A[ix] - a_k) / m_k
    }
  }

  active_blocks <- if (details) {
    seq_along(sel_fit$blocks)
  } else {
    which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids), logical(1)))
  }

  d_s <- sel_fit$design_spec$d_s
  # Selection covariates are defined only on the declared fitting universe.
  # Outside U the nuisance IF is structurally zero; do not read an NA covariate
  # there and then multiply it by a zero score.
  Z_slope <- if (is.null(ref) && d_s > 0L) st$Z[rows[U_rows], sel_fit$design_spec$slope_cols, drop = FALSE] else
    matrix(0, sum(U_rows), d_s)
  block_ifs <- vector("list", length(sel_fit$blocks))
  active_cells <- integer(0)
  u_pos <- cumsum(U_rows)

  for (m in active_blocks) {
    bl <- sel_fit$blocks[[m]]
    cids <- bl$cell_ids
    B_m <- bl$B_m
    active_cells <- c(active_cells, cids)
    cell_scores <- matrix(0, B, B_m)
    # A D60 block lives inside one recipe: its slope rows are that recipe's
    # selection design on the block cohorts' U rows (zero scores elsewhere).
    Zb <- Z_slope
    if (!is.null(ref) && d_s > 0L) {
      rb <- sel_fit$block_specs[[m]]$reference
      on <- which(U_rows & ci %in% vapply(sel_fit$cells[cids], `[[`, integer(1), "cohort"))
      if (length(on)) Zb[u_pos[on], ] <- chained_reference_s(ref, rb, rows[on])[, sel_fit$design_spec$slope_cols, drop = FALSE]
    }
    for (b in seq_len(B_m)) {
      cid <- cids[b]
      cell <- sel_fit$cells[[cid]]
      rr <- which(ci == cell$cohort & U_rows)
      if (!length(rr)) next
      pr <- predict_selection_cell(sel_fit, rows[rr], cid)
      Y <- if (identical(sel_fit$selection_model, "independent_period_logit")) {
        cp04_date_rows(reader, rows[rr], cell$date_spec)$R
      } else {
        cp04_pair_rows(reader, rows[rr], cell$pair)$S
      }
      cell_scores[rr, b] <- w_of(cell$reference, rr) * (Y - pr)
    }
    bif <- list(
      xi_alpha = matrix(0, B, if (bl$has_intercept) B_m else 0L),
      xi_beta = matrix(0, B, bl$d_s),
      xi_theta = matrix(0, B, length(up_map$theta_coords[[m]]))
    )
    if (any(U_rows)) {
      on_U <- evaluate_selection_block_if(bl, Zb, cell_scores[U_rows, , drop = FALSE], st$N)
      for (name in names(bif)) bif[[name]][U_rows, ] <- on_U[[name]]
    }
    block_ifs[[m]] <- bif
    tc <- up_map$theta_coords[[m]]
    keep_tc <- which(tc %in% uids)
    if (length(keep_tc)) {
      up[, match(tc[keep_tc], uids)] <- bif$xi_theta[, keep_tc, drop = FALSE]
    }
  }

  list(
    up = up,
    block_ifs = block_ifs,
    active_blocks = active_blocks,
    active_cells = unique(active_cells)
  )
}

selection_blocks_diagnostic_list <- function(sel_fit) {
  lapply(seq_along(sel_fit$blocks), function(m) {
    bl <- sel_fit$blocks[[m]]
    list(
      block_id = m,
      branch = bl$branch,
      cell_ids = bl$cell_ids,
      B_m = bl$B_m,
      d_s = bl$d_s,
      has_intercept = bl$has_intercept,
      alpha = bl$alpha,
      beta = bl$beta,
      coefficients = c(bl$alpha, bl$beta),
      H_chol = bl$H_chol,
      m_b = bl$m_b,
      z_bar = bl$z_bar,
      pi = bl$pi,
      r = bl$r,
      diagnostics = bl$diagnostics
    )
  })
}

# ------------------------------------------------------------------------------
# 1. Chained IPW with General Selection Models
# ------------------------------------------------------------------------------

cp04_selection_treatment_fit <- function(reader, spec, sel_fit, up_map, glm_control, ctrl, reference = NULL) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p <- if (is.null(reference)) frame$n_covariates else length(reference$x_columns)
  m <- length(spec$cohorts)
  n <- spec$n_observers

  X <- matrix(0, n, p)
  y <- w <- pp <- dy <- numeric(n)
  kk <- obs_rows <- integer(n)
  pos <- 0L

  for (start in seq.int(1L, st$N, by = ctrl$max_chunk_rows)) {
    rows <- seq.int(start, min(st$N, start + ctrl$max_chunk_rows - 1L))
    ci <- match(st$cohort_index[rows], spec$cohorts)
    rows <- rows[!is.na(ci)]
    ci <- ci[!is.na(ci)]
    if (!length(rows)) next
    hit <- cp04_pair_rows(reader, rows, spec$pair, TRUE)
    keep <- which(hit$S == 1)
    if (!length(keep)) next
    rows <- rows[keep]
    ci <- ci[keep]
    ii <- seq.int(pos + 1L, pos + length(rows))
    pos <- pos + length(rows)
    X[ii, ] <- if (is.null(reference)) frame$Z_a[st$a_index[rows], , drop = FALSE] else
      chained_reference_x(reference, spec$reference, rows)
    y[ii] <- as.double(ci == 1L)
    kk[ii] <- ci
    dy[ii] <- hit$dy[keep]
    obs_rows[ii] <- rows
    for (k in unique(ci)) {
      jj <- which(ci == k)
      q <- spec$selection_ids[k]
      pred_k <- predict_pair_selection(sel_fit, rows[jj], q)
      pp[ii[jj]] <- pred_k$prob
      w[ii[jj]] <- pred_k$weight
    }
  }
  if (pos != n) cp04_abort("Observed risk support changed from its owner view.", "identity")

  fit <- cp04_logit(X, y, w, glm_control, st$N, spec)
  lv <- ipw_logit_values(X, fit$coefficients)
  wp <- w * y
  wm <- w * (1 - y) * lv$odds
  dp <- sum(wp) / st$N
  dm <- sum(wm) / st$N
  cp04_finite(c(dp, dm), "treated/control normalizers", TRUE)
  T_val <- sum(wp * dy) / sum(wp)
  M_val <- sum(wm * dy) / sum(wm)
  cp04_finite(c(T_val, M_val), "treated/control means")

  u_ids <- sort(unique(unlist(up_map$pair_channels[spec$selection_ids])))
  u_len <- length(u_ids)
  C <- matrix(0, p, u_len)
  U <- matrix(0, 4L, u_len, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  score <- w * (y - lv$prob) * X
  residual <- cbind(wp * (dy - T_val) / dp, wm * (dy - M_val) / dm, wp, wm)

  for (k in seq_len(m)) {
    rr <- which(kk == k)
    if (!length(rr)) next
    q <- spec$selection_ids[k]
    ch_k <- up_map$pair_channels[[q]]
    pos_k <- match(ch_k, u_ids)
    # Chunked log-derivative evaluation respecting max_chunk_rows
    ld_mat <- matrix(0, length(rr), length(ch_k))
    for (cs in seq.int(1L, length(rr), by = ctrl$max_chunk_rows)) {
      ce <- min(length(rr), cs + ctrl$max_chunk_rows - 1L)
      sub <- cs:ce
      ld_mat[sub, ] <- selection_pair_cohort_log_deriv(sel_fit, up_map, q, obs_rows[rr[sub]])$log_deriv
    }
    C[, pos_k] <- C[, pos_k] + crossprod(score[rr, , drop = FALSE], ld_mat) / st$N
    U[, pos_k] <- U[, pos_k] + crossprod(residual[rr, , drop = FALSE], ld_mat) / st$N
  }

  V <- matrix(0, 4L, p, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  wm_act <- wm * lv$active
  V[2L, ] <- colSums(wm_act * (dy - M_val) * X) / (st$N * dm)
  V[4L, ] <- colSums(wm_act * X) / st$N
  cp04_finite(c(C, U, V), "analytic cross derivatives")
  fit$upstream_ids <- as.integer(u_ids)
  fit$cross <- C
  fit$partials <- list(upstream = U, beta = V)
  fit$means <- c(T = T_val, M = M_val, dp = dp, dm = dm)
  fit$clipping <- ipw_clip_summary(prob = lv$prob, q = pp, controls = y == 0, weights = w)
  fit
}

# `reference` (D60): selection cells, rate coordinates and treatment logits are
# keyed by recipe; designs and weights are read through the sealed reader.
fit_cp04_selection_models <- function(reader, link_views, control = ipw_owner_control(),
                                      glm_control,
                                      retained_external_entries = 0, external_work_ops = 0,
                                      selection_model = "pair_logit",
                                      selection_pooling = "separate",
                                      sel_fit = NULL, reference = NULL) {
  validate_cp04_reader(reader)
  ctrl <- validate_ipw_owner_control(control)
  gc <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    cp04_abort("Invalid external pipeline reservation.")
  st <- reader$state
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p <- if (is.null(reference)) frame$n_covariates else length(reference$x_columns)
  d <- if (is.null(reference)) ncol(st$Z) else length(reference$designs[[reference$role_slot["s"]]]$columns)
  if (!is.list(link_views) || !length(link_views))
    cp04_abort("Supply actual private link views.")

  if (is.null(sel_fit)) {
    sel_fit <- fit_selection_models(
      reader = reader,
      link_views = link_views,
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      control = ctrl,
      glm_control = gc,
      retained_external_entries = retained_external_entries,
      external_work_ops = external_work_ops,
      reference = reference
    )
  } else {
    validate_selection_models(sel_fit)
    if (!identical(sel_fit$reader_ref, reader$ref) ||
        !identical(sel_fit$selection_model, selection_model) ||
        !identical(sel_fit$selection_pooling, selection_pooling) ||
        !identical(is.null(sel_fit$reference_ref), is.null(reference)))
      cp04_abort("Supplied `sel_fit` does not match reader or selection specification.", "integrity")
  }
  up_map <- build_selection_upstream_map(sel_fit)

  L <- length(link_views)
  tkeys <- character()
  ts <- list()
  links <- vector("list", L)

  for (i in seq_len(L)) {
    v <- link_views[[i]]
    pair <- list(s = v$pair$s, t = v$pair$t,
                 s_calendar_index = v$pair$s_calendar_index,
                 t_calendar_index = v$pair$t_calendar_index)
    target <- match(v$cohorts$g, st$cohorts)
    controls <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    cohorts <- c(target, controls)
    r <- if (is.null(reference)) NULL else reference$view_recipe[i]
    si <- vapply(cohorts, function(k) {
      resolve_pair_cohort_id(sel_fit, if (is.null(r)) list(pair = pair, cohort = k) else
        list(pair = pair, cohort = k, reference = r))
    }, integer(1))
    tkey <- paste(pair$s_calendar_index, pair$t_calendar_index, target,
                  paste(controls, collapse = ","), sep = ":")
    if (!is.null(r)) tkey <- paste0(tkey, "|recipe:", r)
    tid <- match(tkey, tkeys)
    if (is.na(tid)) {
      tid <- length(ts) + 1L
      tkeys[tid] <- tkey
      ts[[tid]] <- list(
        pair = pair,
        cohorts = cohorts,
        selection_ids = si,
        n_observers = as.integer(v$support$n_observed_risk_units)
      )
      if (!is.null(r)) ts[[tid]]$reference <- as.integer(r)
    }
    source <- v$functional_key
    links[[i]] <- list(
      link_id = v$link_id,
      model = tid,
      pair = c(s = pair$s, t = pair$t),
      g = v$cohorts$g,
      K = sort(v$cohorts$K, method = "radix"),
      source_functional_key = source,
      functional_key = paste0(
        "cmpcat:CP04-", st$universe, "-weighted-",
        selection_model, "-", selection_pooling, "-",
        if (all(st$rho[st$U] == 1)) "a/p" else "a/(rho*p)", "-etaQG|",
        substring(source, 8L), if (is.null(r)) "" else paste0("|reference:", reference$recipes$key[r])
      )
    )
  }

  T_count <- length(ts)
  max_m <- max(vapply(ts, function(x) length(x$cohorts), integer(1)))
  max_trt_n <- max(vapply(ts, `[[`, integer(1), "n_observers"))
  u_max <- up_map$D
  trt_models_bound <- T_count * (p * p + 6 * p + 4 * u_max + p * u_max + 2048) + 8192
  trt_scratch <- max_trt_n * (p + d + 16) + 32 * (p + u_max + 4)^2 + 8192
  trt_scan_n <- sum(vapply(ts, `[[`, integer(1), "n_observers"))
  trt_work <- 128 * (gc$maxit + 2) * trt_scan_n * (p + 2)^2 +
    128 * T_count * st$N * max_m + 128 * trt_models_bound
  base_retained <- sel_fit$accounting$owned_entries + sel_fit$accounting$borrowed_entries
  cum_work <- ipw_owner_safe_sum(
    c(sel_fit$accounting$work_ops_bound, trt_work),
    ctrl$max_work_ops,
    "CP04 selection-models plus treatment pipeline fit work"
  )
  cp04_budget(ctrl, c(base_retained, trt_models_bound), trt_scratch, c(cum_work, external_work_ops),
              sel_fit$accounting$owned_entries + trt_models_bound, "CP04 selection-models fits")

  treatments <- lapply(ts, function(s) cp04_selection_treatment_fit(reader, s, sel_fit, up_map, gc, ctrl, reference))
  selections <- selection_blocks_diagnostic_list(sel_fit)
  selection_specs <- sel_fit$block_specs

  state <- list(
    uses_selection_models = TRUE,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    selection_backend = sel_fit,
    up_map = up_map,
    selection_specs = selection_specs,
    selections = selections,
    treatment_specs = ts,
    treatments = treatments,
    links = links,
    reader_ref = reader$ref,
    control = ctrl,
    glm_control = gc,
    universe = st$universe,
    justification = st$justification,
    law = paste0(
      "CP04:", st$universe, ":", selection_model, ":", selection_pooling,
      if (is.null(reference)) ":w*a_U/(rho_U*p_U):eta_k=P(G=k|A=1):common-frame-Z_a" else
        ":W_g*a_U/(rho_U*p_U):recipe-indexed-reference-design"
    ),
    accounting = list(
      owned_entries = 0,
      borrowed_entries = reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
        retained_external_entries,
      retained_external_entries = retained_external_entries,
      external_work_ops = external_work_ops,
      retained_bound = base_retained + trt_models_bound,
      fit_scratch_bound = max(sel_fit$accounting$fit_scratch_bound, trt_scratch),
      work_ops_bound = cum_work,
      n_selection_fits = length(sel_fit$blocks),
      n_selection_cells = length(sel_fit$cells),
      n_selection_blocks = length(sel_fit$blocks),
      n_selection_parameters = length(sel_fit$coefficients),
      n_treatment_fits = T_count,
      n_rate_fits = up_map$n_a,
      rates_from_reader_shares = TRUE,
      full_F_masks_per_model = 0L,
      full_F_by_link_IF_retained = FALSE,
      inverse_matrices_retained = FALSE,
      peak_rss_measured = FALSE
    )
  )
  class(state) <- c("cdid_cp04_selection_backend", "cdid_cp04_corollary1")
  if (!is.null(reference)) {
    state$reference_ref <- reference$ref
    state$accounting$n_reference_recipes <- nrow(reference$recipes)
  }
  state["ref"] <- list(NULL)
  actual <- unname(cp04_payload(state)["entries"])
  state$accounting$owned_entries <- actual
  state$accounting$retained_bound <- max(state$accounting$retained_bound, state$accounting$borrowed_entries + actual)
  cp04_budget(ctrl, c(state$accounting$borrowed_entries, actual), state$accounting$fit_scratch_bound,
              c(cum_work, external_work_ops), actual, "retained CP04 selection-models backend")
  ipw_owner_seal(state, "cp04_backend")
}

cp04_selection_if_budget <- function(backend, B, links, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p <- if (is.null(backend$reference_ref)) frame$n_covariates else length(chained_backend_reference(backend)$x_columns)
  sel_fit <- backend$selection_backend
  up_map <- backend$up_map

  tids <- unique(vapply(backend$links[links], `[[`, integer(1), "model"))
  sids <- unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "selection_ids")))
  uids <- if (details) {
    seq_len(up_map$D)
  } else {
    sort(unique(unlist(lapply(backend$treatments[tids], `[[`, "upstream_ids"))))
  }
  active_blocks <- if (details) {
    seq_along(sel_fit$blocks)
  } else {
    which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids), logical(1)))
  }
  active_cells <- unique(unlist(lapply(sel_fit$block_specs[active_blocks], `[[`, "cell_ids")))
  nu <- length(uids)
  nc <- length(active_cells)

  widths <- length(links) + length(tids) * (2 * p + nu + 12) +
    if (details) 3 * length(st$cohorts) + nu + nc else 0
  max_bytes <- frame$id_byte_metadata$max
  link_bytes <- sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  metadata <- 8192 + 1024 * (length(sel_fit$blocks) + length(backend$treatments) +
    length(links) + length(st$cohorts)) + B * max_bytes + 8 * link_bytes
  output <- B * (widths + 16) + metadata
  scratch <- 32 * B * (widths + p + sel_fit$design_spec$d + 8) +
    64 * ((p + 2)^2 + (sel_fit$design_spec$d_s + 2)^2) + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- 128 * B * (
    (nc + 1) * (sel_fit$design_spec$d_s + 2)^2 * log_a +
      length(tids) * ((p + 2)^2 * log_a + p * (nu + 4)) +
      length(links)
  ) + 64 * metadata
  # Shared selection slopes consume every score cell in their active blocks,
  # including cells outside the requested treatment supports. Reuse the
  # accepted selection envelope for those transient score and IF matrices.
  selection_pre <- selection_if_budget(sel_fit, B, sids, details = details)
  scratch <- scratch + selection_pre$output + selection_pre$scratch
  work <- work + selection_pre$work
  cp04_budget(
    backend$control,
    c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    scratch,
    c(backend$accounting$work_ops_bound, work),
    output,
    "CP04 selection-models IF chunk"
  )
  list(
    output = output,
    scratch = scratch,
    work = work,
    tids = tids,
    sids = sids,
    uids = uids,
    active_blocks = active_blocks,
    active_cells = active_cells
  )
}

evaluate_cp04_selection_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_cp04_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) ||
      anyNA(rows) || any(rows != floor(rows)) || any(rows < 1 | rows > st$N) ||
      anyDuplicated(rows) || length(rows) > backend$control$max_chunk_rows ||
      !is.numeric(links) || is.object(links) || !is.null(dim(links)) || !length(links) ||
      anyNA(links) || any(links != floor(links)) || any(links < 1 | links > length(backend$links)) ||
      anyDuplicated(links) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    cp04_abort("Invalid explicit IF row/link batch, scale or details flag.")

  pre <- cp04_selection_if_budget(backend, length(rows), links, details)
  tids <- pre$tids
  sids <- pre$sids
  uids <- pre$uids
  B <- length(rows)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p <- if (is.null(backend$reference_ref)) frame$n_covariates else length(chained_backend_reference(backend)$x_columns)
  sel_fit <- backend$selection_backend
  up_map <- backend$up_map

  ref <- chained_backend_reference(backend)
  ci <- st$cohort_index[rows]
  w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
  A <- as.double(st$a_index[rows] > 0L)
  U_rows <- st$U[rows]

  up_res <- evaluate_selection_upstream_chunk(sel_fit, up_map, rows, uids, details = details)
  up <- up_res$up

  models <- vector("list", length(backend$treatments))
  for (tj in tids) {
    spec <- backend$treatment_specs[[tj]]
    f <- backend$treatments[[tj]]
    u_pos <- match(f$upstream_ids, uids)
    upstream <- up[, u_pos, drop = FALSE]

    ci_in <- which(ci %in% spec$cohorts)
    hit_all <- if (length(ci_in)) cp04_pair_rows(reader, rows[ci_in], spec$pair, TRUE) else list(S = numeric(0), dy = numeric(0))
    obs_local <- which(hit_all$S == 1)
    rr <- ci_in[obs_local]

    beta_score <- matrix(0, B, p)
    direct <- matrix(0, B, 4L)
    direct[, 3L] <- -f$means["dp"]
    direct[, 4L] <- -f$means["dm"]
    if (length(rr)) {
      X <- if (is.null(ref)) frame$Z_a[st$a_index[rows[rr]], , drop = FALSE] else
        chained_reference_x(ref, spec$reference, rows[rr])
      dy_rr <- hit_all$dy[obs_local]
      lv <- ipw_logit_values(X, f$coefficients)
      cp04_finite(c(lv$prob, lv$eta, lv$odds), "treatment replay")
      D_ind <- as.double(ci[rr] == spec$cohorts[1L])
      fw <- numeric(length(rr))
      for (k in seq_along(spec$cohorts)) {
        sub_k <- which(ci[rr] == spec$cohorts[k])
        if (length(sub_k)) {
          fw[sub_k] <- predict_pair_selection(sel_fit, rows[rr[sub_k]], spec$selection_ids[k])$weight
        }
      }
      wp <- fw * D_ind
      wm <- fw * (1 - D_ind) * lv$odds
      beta_score[rr, ] <- fw * (D_ind - lv$prob) * X
      direct[rr, 1L] <- wp * (dy_rr - f$means["T"]) / f$means["dp"]
      direct[rr, 2L] <- wm * (dy_rr - f$means["M"]) / f$means["dm"]
      direct[rr, 3L] <- direct[rr, 3L] + wp
      direct[rr, 4L] <- direct[rr, 4L] + wm
    }
    xb <- t(cp04_solve(f$H_chol, t(beta_score + upstream %*% t(f$cross))))
    coord <- direct + upstream %*% t(f$partials$upstream) + xb %*% t(f$partials$beta)
    colnames(coord) <- c("T", "M", "dp", "dm")
    cp04_finite(c(xb, coord), "complete frame IF")
    models[[tj]] <- if (details) {
      list(
        beta = xb,
        coordinates = coord,
        direct = direct,
        beta_score = beta_score,
        upstream = upstream,
        normalizers_A = coord[, 3:4, drop = FALSE] / st$lambda -
          outer(A - st$lambda, f$means[c("dp", "dm")] / st$lambda^2)
      )
    } else {
      list(phi = coord[, 1L] - coord[, 2L])
    }
  }

  phi <- matrix(0, B, length(links), dimnames = list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id")))
  for (j in seq_along(links)) {
    m <- models[[backend$links[[links[j]]]$model]]
    phi[, j] <- if (details) m$coordinates[, 1L] - m$coordinates[, 2L] else m$phi
  }
  if (scale == "core") phi <- st$lambda * phi

  out <- list(
    rows = rows,
    id = st$id[rows],
    links = links,
    phi = phi,
    scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(
      n_selection_score_evaluations = length(pre$active_cells),
      n_selection_block_evaluations = length(pre$active_blocks),
      n_treatment_score_evaluations = length(tids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  if (details) {
    # D60 rates/shares belong to each recipe (read them from `upstream`);
    # frame-level a/tau/eta summaries exist only on the static path.
    a_if <- tau <- eta <- NULL
    if (is.null(ref)) {
      a_if <- matrix(0, B, length(st$cohorts))
      for (k in unique(ci)) {
        rr_k <- which(ci == k & U_rows)
        if (length(rr_k)) a_if[rr_k, k] <- w[rr_k] * (A[rr_k] - st$a[k]) / st$mass_U[k]
      }
      tau <- eta <- matrix(0, B, length(st$cohorts))
      for (k in seq_along(st$cohorts)) {
        tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
        eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
      }
    }
    alpha_list <- lapply(up_res$block_ifs, function(bif) if (is.null(bif)) NULL else bif$xi_theta)
    out$details <- list(
      a = a_if,
      alpha = alpha_list,
      upstream = up,
      upstream_ids = uids,
      models = models,
      tau = tau,
      eta = eta,
      lambda = A - st$lambda
    )
  }
  actual_output <- unname(cp04_payload(out)["entries"])
  if (actual_output > pre$output) cp04_abort("IF output bound underestimated.", "bound")
  out
}

# ------------------------------------------------------------------------------
# 2. Shared-Slope OR with General Selection Models
# ------------------------------------------------------------------------------

or_selection_moments_add <- function(M, z, dy, q, T_logit) {
  if (!ncol(T_logit)) return(M)
  x <- cbind(1, z, dy)
  for (j in seq_len(ncol(T_logit))) {
    M[, , j] <- M[, , j] + crossprod(x, (q * T_logit[, j]) * x)
  }
  M
}

or_selection_support <- function(reader, spec, sel_fit, up_map, q_id, recipe, ctrl, reference = NULL) {
  st <- reader$state
  a <- if (is.null(reference)) st$a[spec$cohort] else chained_reference_stat(reference, spec$reference, spec$cohort, "a")
  m <- recipe$p + 2L
  dq <- length(up_map$pair_channels[[q_id]]) - 1L
  acc <- or_pooled_support_start(recipe, spec$transitions)
  moments <- array(0, c(m, m, dq))
  clip <- ipw_clip_summary(q = numeric(0), weights = numeric(0))

  or_baseline_visit(reader, spec, ctrl, function(rows, S, obs, dy, Z) {
    if (!length(obs)) return(invisible(NULL))
    rr <- rows[obs]
    ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, q_id, rr)
    q <- ld$pred$weight
    cp04_finite(q, "Corollary 1 fitted pair weights", TRUE)
    clip <<- ipw_clip_summary_merge(clip, ipw_clip_summary(q = ld$pred$prob, weights = q))
    acc <<- or_pooled_support_add(acc, recipe, Z, dy, q, "corollary1")
    moments <<- or_selection_moments_add(moments, Z[, recipe$slope_index, drop = FALSE], dy, q, ld$T_logit)
  }, reference)
  if (!is.finite(acc$weight_mass) || acc$weight_mass <= 0 || acc$weight_mass / st$N <= 0)
    or_corollary1_abort("Positive actual-pair cohort support required.", "support")
  cp04_finite(c(acc$sumZ, acc$sumY, moments), "Corollary 1 compact support moments")
  list(design = acc, moments = moments, a = a, pair_mass = acc$weight_mass, n_pair = acc$n,
    clipping = clip)
}

or_selection_forcing <- function(recipe, fitted, supports, pools, channels, D, N) {
  keep <- fitted$keep
  p <- length(keep)
  m <- recipe$p + 2L
  beta_total <- matrix(0, p, D)
  energy_total <- numeric(D)
  regressions <- fitted$regressions

  for (b in seq_along(pools)) {
    sids <- pools[[b]]$support_ids
    rec <- regressions[[b]]
    fa <- fe <- numeric(D)
    fb <- matrix(0, p, D)
    cr <- numeric(m)
    cr[1L] <- -rec$alpha
    cr[1L + keep] <- -rec$L * fitted$beta
    cr[m] <- 1
    for (s in sids) {
      u <- supports[[s]]
      ka <- channels[[s]][1L]
      kl <- channels[[s]][-1L]
      g <- or_pooled_factor_moments(recipe, u$design$summary, keep, rec$alpha, fitted$beta)
      if (recipe$intercept) {
        fa[ka] <- fa[ka] + g$g[1L] / (N * u$a)
        fb[, ka] <- fb[, ka] + (g$g[-1L] - rec$mean_d * g$g[1L]) / (N * u$a)
      } else {
        fb[, ka] <- fb[, ka] + g$g / (N * u$a)
      }
      fe[ka] <- fe[ka] + g$e2 / (N * u$a)

      if (length(kl) > 0L) {
        Mr <- vapply(seq_along(kl), function(j) as.vector(u$moments[, , j] %*% cr), numeric(m))
        if (!is.matrix(Mr)) Mr <- matrix(Mr, nrow = m, ncol = length(kl))
        if (recipe$intercept) {
          fa[kl] <- fa[kl] + Mr[1L, ] / N
          fb[, kl] <- fb[, kl] + (rec$L * Mr[1L + keep, , drop = FALSE] -
            outer(rec$mean_d, Mr[1L, ])) / N
        } else {
          fb[, kl] <- fb[, kl] + rec$L * Mr[1L + keep, , drop = FALSE] / N
        }
        fe[kl] <- fe[kl] + colSums(cr * Mr) / N
      }
    }
    ids <- sort(unique(unlist(channels[sids])))
    cp04_finite(c(fa, fb, fe), "Corollary 1 regression forcing")
    regressions[[b]]$channel_ids <- ids
    regressions[[b]]$forcing_alpha <- fa[ids]
    regressions[[b]]$forcing_beta <- fb[, ids, drop = FALSE]
    regressions[[b]]$forcing_energy <- fe[ids]
    beta_total <- beta_total + fb
    energy_total <- energy_total + fe
  }
  list(regressions = regressions, beta = beta_total, energy = energy_total)
}

or_selection_link <- function(recipe, fitted, u, rec, N) {
  keep <- fitted$keep
  p <- length(keep)
  m <- recipe$p + 2L
  dq <- dim(u$moments)[3L]
  lm <- or_pooled_link_moments(recipe, keep, u$pair_mass, u$design$sumZ, u$design$sumY,
    rec$L, rec$alpha, fitted$beta)
  d <- u$pair_mass / N
  num <- lm$mu * d
  if (dq > 0L) {
    cr <- numeric(m)
    cr[1L] <- -rec$alpha
    cr[1L + keep] <- -rec$L * fitted$beta
    cr[m] <- 1
    Mr <- vapply(seq_len(dq), function(j) as.vector(u$moments[, , j] %*% cr), numeric(m))
    if (!is.matrix(Mr)) Mr <- matrix(Mr, nrow = m, ncol = dq)
    ans <- list(
      mu = lm$mu, c = lm$c, d = d, numerator = num,
      rate_loading = c(lm$residual / (N * u$a * d),
        (Mr[1L, ] - lm$mu * u$moments[1L, 1L, ]) / (N * d)),
      weight_loading = c(d / u$a, u$moments[1L, 1L, ] / N),
      numerator_loading = c(num / u$a, Mr[1L, ] / N),
      design_loading = matrix(c(lm$c * d / u$a, rec$L * u$moments[1L + keep, 1L, ] / N),
        p, 1L + dq)
    )
  } else {
    ans <- list(
      mu = lm$mu, c = lm$c, d = d, numerator = num,
      rate_loading = c(lm$residual / (N * u$a * d)),
      weight_loading = c(d / u$a),
      numerator_loading = c(num / u$a),
      design_loading = matrix(lm$c * d / u$a, p, 1L)
    )
  }
  cp04_finite(unlist(ans), "Corollary 1 link moments")
  ans
}

or_selection_assemble <- function(recipe, supports, pools, links, channels, D, rank_tol, N, visit_rows) {
  blocks <- lapply(seq_along(pools), function(b) list(
    name = paste0("b", b),
    support_ids = pools[[b]]$support_ids,
    L = pools[[b]]$transitions
  ))
  fitted <- or_pooled_fit_blocks(
    recipe,
    lapply(supports, function(u) list(summary = u$design$summary, scale = 1, colmax = u$design$colmax)),
    blocks,
    lapply(links, function(l) list(regression = l$regression, support_id = l$support_id)),
    rank_tol, N, visit_rows, "corollary1"
  )
  forcing <- or_selection_forcing(recipe, fitted, supports, pools, channels, D, N)
  for (i in seq_along(links)) {
    l <- links[[i]]
    links[[i]] <- c(l, or_selection_link(recipe, fitted, supports[[l$support_id]],
      forcing$regressions[[l$regression]], N))
  }
  beta <- fitted$beta
  names(beta) <- fitted$basis$kept_columns
  list(
    supports = lapply(seq_along(supports), function(s) {
      u <- supports[[s]]
      list(
        a = u$a, channels = channels[[s]], pair_mass = u$pair_mass, n_pair = u$n_pair,
        sumZ = u$design$sumZ, sumY = u$design$sumY, colmax = u$design$colmax, L = u$design$L,
        clipping = u$clipping
      )
    }),
    regressions = forcing$regressions,
    links = links,
    shared = list(
      p = length(fitted$keep), keep = fitted$keep, beta = fitted$beta,
      coefficients = beta, fit = fitted$fit, summary = fitted$summary, basis = fitted$basis,
      forcing_beta = forcing$beta, forcing_energy = forcing$energy,
      control_support_ids = sort(unique(unlist(lapply(pools, `[[`, "support_ids")))),
      rss = fitted$rss, objective = fitted$objective,
      verification_scan = fitted$verification_scan
    )
  )
}

# `reference` (D60): supports (pair, cohort, recipe), blocks (pair, pool,
# recipe) and selection cells keyed by recipe; ONE global shared slope.
fit_or_selection_models <- function(reader, link_views, control, glm_control, rank_tol, outcome_recipe,
                                    retained_external_entries = 0, external_work_ops = 0,
                                    automatic_chunk_rows = FALSE,
                                    selection_model = "pair_logit",
                                    selection_pooling = "separate",
                                    sel_fit = NULL, reference = NULL) {
  validate_cp04_reader(reader)
  if (missing(rank_tol) || !gmm_scalar_numeric(rank_tol) || rank_tol <= 0 || rank_tol >= 1)
    or_corollary1_abort("Supply explicit finite private rank_tol in (0, 1).")
  ctrl <- validate_ipw_owner_control(control)
  gc <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    or_corollary1_abort("Invalid external pipeline reservation.")
  if (!identical(automatic_chunk_rows, isTRUE(automatic_chunk_rows)))
    or_corollary1_abort("Invalid private automatic row-planning flag.")
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (missing(outcome_recipe))
    or_corollary1_abort("Supply the explicit private outcome recipe (intercept flag and Z_a roles).")
  if (!is.null(reference)) {
    validate_chained_reference_reader(reference)
    if (!is.list(link_views) || length(reference$view_recipe) != length(link_views))
      or_corollary1_abort("The reference reader must index exactly the supplied link views.", "integrity")
  }
  recipe <- or_pooled_check_recipe(outcome_recipe, if (is.null(reference)) colnames(frame$Z_a) else reference$columns,
    "corollary1")
  if (!is.list(link_views) || !length(link_views))
    or_corollary1_abort("Supply actual private link views.")

  L <- length(link_views)
  J <- length(st$cohorts)
  p <- recipe$p
  d <- if (is.null(reference)) ncol(st$Z) else length(reference$designs[[reference$role_slot["s"]]]$columns)
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    retained_external_entries + if (is.null(reference)) 0 else reference$accounting$owned_entries
  header <- 4096 + 1024 * L * (J + 1)
  cp04_budget(ctrl, c(base, header), header, c(external_work_ops, 128 * L * L * (J + 1)), header,
    "Corollary 1 OR metadata")

  records <- list()
  for (ell in seq_along(link_views)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:"))
      or_corollary1_abort("Corollary 1 OR requires the supplied comparison-catalog functional key.")
    key_bytes <- nchar(source, type = "bytes")
    cp04_budget(ctrl, c(base, header), header + 16 * key_bytes,
      c(external_work_ops, 128 * L * L * (J + 1) + 64 * key_bytes), header + 16 * key_bytes,
      "Corollary 1 catalog key")
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
      or_corollary1_abort("Legacy nuisance candidates are unsupported.")
    if (!is.null(v$nuisance_ref)) {
      nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair))
        or_corollary1_abort("Cross-pair or unverified nuisance sharing is unsupported.")
    }
    target <- match(v$cohorts$g, st$cohorts)
    control_cohorts <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    if (anyNA(c(target, control_cohorts)) || !isTRUE(all(v$support$by_cohort_pair_count >= 1)))
      or_corollary1_abort("Positive actual-pair cohort support required.", "support")
    or_pooled_transitions(v$pair, "corollary1")
    records[[ell]] <- list(
      pair = v$pair, target = target, controls = control_cohorts,
      link_id = v$link_id, source_functional_key = source,
      functional_key = paste0(
        "cmpcat:OR-Cor1-", st$universe, "-", selection_model, "-", selection_pooling,
        "-shared-slope-WLS-a/(rho*p)|", substring(source, 8L)
      )
    )
    if (!is.null(reference)) {
      records[[ell]]$reference <- reference$view_recipe[ell]
      records[[ell]]$functional_key <- paste0(records[[ell]]$functional_key, "|reference:",
        reference$recipes$key[reference$view_recipe[ell]])
    }
    cp04_budget(ctrl, c(base, header), header + cp04_payload(records)["entries"],
      c(external_work_ops, 128 * L * L * (J + 1) + 64 * cp04_payload(records)["char_bytes"]),
      header + cp04_payload(records)["entries"], "Corollary 1 catalog metadata")
  }
  ids <- vapply(records, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids))
    or_corollary1_abort("Distinct nonempty link IDs required.")

  if (is.null(sel_fit)) {
    sel_fit <- fit_selection_models(
      reader = reader,
      link_views = link_views,
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      control = ctrl,
      glm_control = gc,
      retained_external_entries = retained_external_entries,
      external_work_ops = external_work_ops,
      reference = reference
    )
  } else {
    validate_selection_models(sel_fit)
    if (!identical(sel_fit$reader_ref, reader$ref) ||
        !identical(sel_fit$selection_model, selection_model) ||
        !identical(sel_fit$selection_pooling, selection_pooling) ||
        !identical(is.null(sel_fit$reference_ref), is.null(reference)))
      or_corollary1_abort("Supplied `sel_fit` does not match reader or selection specification.", "integrity")
  }

  catalog <- or_baseline_catalog(records, st$cohorts, recipe, reference)
  specs <- catalog$specs
  pools <- catalog$pools
  links <- catalog$links
  for (i in seq_along(pools)) {
    pr <- specs[[pools[[i]]$support_ids[1L]]]$pair
    pools[[i]]$block_key <- paste0(
      "pair:", pr$s, "-", pr$t, "|pool:",
      paste(st$cohorts[pools[[i]]$cohorts], collapse = ","), "|",
      if (is.null(reference)) recipe$key else paste0("reference-design:L_e*z(C^r)|",
        reference$recipes$key[pools[[i]]$reference]),
      if (is.null(reference)) "|weights:w*a_U/(rho_U*p_U)|U:" else "|weights:W_g*a_U/(rho_U*p_U)|U:", st$universe,
      "|model:", selection_model, "|pooling:", selection_pooling
    )
  }

  S <- length(specs)
  G <- length(pools)
  cohorts_a <- sort(unique(vapply(specs, `[[`, integer(1), "cohort")))
  up_map <- build_selection_upstream_map(sel_fit, cohorts_a = cohorts_a)
  n_a <- up_map$n_a
  D <- up_map$D
  q_of_s <- vapply(specs, function(s) {
    resolve_pair_cohort_id(sel_fit, if (is.null(s$reference)) list(pair = s$pair, cohort = s$cohort) else
      list(pair = s$pair, cohort = s$cohort, reference = s$reference))
  }, integer(1))
  channels <- lapply(seq_len(S), function(s) up_map$pair_channels[[q_of_s[s]]])

  upstream <- list(
    n_a = n_a,
    d = up_map$d,
    D = D,
    cohort = up_map$cohort,
    support = up_map$support,
    cell = up_map$cell,
    column = up_map$column,
    role = up_map$role
  )
  selection_specs <- sel_fit$block_specs
  selections <- selection_blocks_diagnostic_list(sel_fit)

  metadata <- unname(cp04_payload(list(specs, pools, links, selection_specs, upstream, sel_fit))["entries"])
  cohort_n <- vapply(specs, function(s) st$cohort_end[s$cohort] - st$cohort_start[s$cohort] + 1, numeric(1))
  charges <- or_corollary1_fit_charges(
    S, G, L, p, max(1L, up_map$d), D, st$n_A, min(st$N, ctrl$max_chunk_rows),
    max(cohort_n), sum(cohort_n), gc$maxit, metadata, header,
    sum(nchar(c(recipe$columns, if (is.null(reference)) colnames(st$Z) else
      reference$designs[[reference$role_slot["s"]]]$columns), type = "bytes"))
  )
  scratch <- max(sel_fit$accounting$fit_scratch_bound, charges$scratch)
  work <- ipw_owner_safe_sum(
    c(sel_fit$accounting$work_ops_bound, charges$work),
    ctrl$max_work_ops,
    "Corollary 1 OR selection-models pipeline fit work"
  )
  bound <- metadata + charges$models
  cp04_budget(ctrl, c(base, bound), scratch, c(work, external_work_ops), bound, "Corollary 1 OR selection-models fits")

  supports <- lapply(seq_len(S), function(j) {
    or_selection_support(reader, specs[[j]], sel_fit, up_map, q_of_s[j], recipe, ctrl, reference)
  })
  state <- or_selection_assemble(
    recipe, supports, pools, links, channels, D, rank_tol, st$N,
    function(sid, visit) or_baseline_visit(reader, specs[[sid]], ctrl,
      function(rows, S_vec, obs, dy, Z) visit(Z), reference)
  )
  rm(supports)

  obj <- structure(
    list(
      uses_selection_models = TRUE,
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      selection_backend = sel_fit,
      up_map = up_map,
      q_of_s = q_of_s,
      reader_ref = reader$ref,
      control = ctrl,
      glm_control = gc,
      rank_tol = rank_tol,
      outcome_recipe = recipe,
      support_specs = specs,
      selection_specs = selection_specs,
      selections = selections,
      upstream = upstream,
      supports = state$supports,
      regression_specs = pools,
      regressions = state$regressions,
      shared = state$shared,
      links = state$links,
      law = paste0(
        "OR:Cor1:", st$universe, ":", selection_model, ":", selection_pooling,
        ":control-only:shared-one-period-slope:",
        if (recipe$intercept) "block-intercepts" else "no-intercept",
        if (is.null(reference)) ":w*a_U/(rho_U*p_U):common-frame-Z_a" else
          ":W_g*a_U/(rho_U*p_U):recipe-indexed-reference-design"
      ),
      accounting = list(
        owned_entries = 0,
        borrowed_entries = base,
        retained_external_entries = retained_external_entries,
        external_work_ops = external_work_ops,
        retained_bound = base + bound,
        fit_scratch_bound = scratch,
        work_ops_bound = work,
        support_cohort_row_visits = sum(cohort_n),
        automatic_chunk_rows = automatic_chunk_rows,
        n_cohort_summaries = S,
        n_rate_fits = n_a,
        rates_from_reader_shares = TRUE,
        n_selection_fits = length(sel_fit$blocks),
        n_selection_cells = length(sel_fit$cells),
        n_selection_blocks = length(sel_fit$blocks),
        n_selection_parameters = length(sel_fit$coefficients),
        n_upstream_coordinates = D,
        n_regression_fits = 1L,
        n_regression_blocks = G,
        n_local_intercepts = if (recipe$intercept) G else 0L,
        n_shared_slopes = state$shared$p,
        n_slope_columns = p,
        n_aliased_slopes = length(state$shared$basis$aliased_columns),
        n_support_scans = S * (1L + state$shared$verification_scan),
        n_treatment_fits = 0L,
        full_F_masks_per_model = 0L,
        full_F_by_link_IF_retained = FALSE,
        inverse_matrices_retained = FALSE,
        dense_dummy_regression = FALSE,
        peak_rss_measured = FALSE
      )
    ),
    class = c("cdid_or_selection_backend", "cdid_or_corollary1")
  )
  if (!is.null(reference)) {
    obj$reference_ref <- reference$ref
    obj$accounting$n_reference_recipes <- nrow(reference$recipes)
  }
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  obj$accounting$owned_entries <- actual
  obj$accounting$retained_bound <- max(obj$accounting$retained_bound, base + actual)
  cp04_budget(ctrl, c(base, actual), scratch, c(work, external_work_ops), actual, "retained Corollary 1 OR selection-models")
  ipw_owner_seal(obj, "or_corollary1_backend")
}

evaluate_or_selection_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_or_corollary1_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  valid_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) &&
    all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid_index(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    or_corollary1_abort("Invalid explicit IF row/link batch, scale or details flag.")

  pre <- or_corollary1_if_budget(backend, length(rows), links, details)
  B <- length(rows)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  kept <- backend$outcome_recipe$slope_index[backend$shared$keep]
  ref <- chained_backend_reference(backend)
  ci <- st$cohort_index[rows]
  # D60: no common weight; supports read their recipe's W_g on observed rows.
  w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
  U <- st$U[rows]
  A <- as.double(st$a_index[rows] > 0L)

  sel_fit <- backend$selection_backend
  up_map <- backend$up_map
  up_res <- evaluate_selection_upstream_chunk(sel_fit, up_map, rows, pre$uids, details = details)
  up <- up_res$up

  hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]
    sp <- backend$support_specs[[sid]]
    q <- backend$q_of_s[sid]
    ix <- which(ci == sp$cohort & U)
    hit <- cp04_pair_rows(reader, rows[ix], sp$pair, TRUE)
    obs <- which(hit$S == 1)
    rr <- rows[ix[obs]]
    pr_obs <- if (length(rr)) predict_pair_selection(sel_fit, rr, q)$prob else numeric(0)
    hits[[j]] <- if (is.null(ref)) list(
      obs = ix[obs],
      dy = hit$dy[obs],
      z = frame$Z_a[st$a_index[rr], kept, drop = FALSE],
      v = or_corollary1_factor(st$a[sp$cohort], st$rho[rr], pr_obs)
    ) else list(
      obs = ix[obs],
      dy = hit$dy[obs],
      z = chained_reference_design(ref, sp$reference, rr)[, kept, drop = FALSE],
      v = or_corollary1_factor(chained_reference_stat(ref, sp$reference, sp$cohort, "a"), st$rho[rr], pr_obs),
      w = chained_reference_weight(ref, sp$reference, rr)
    )
  }

  core <- or_pooled_if_core(backend, st$N, w, hits, up, pre, links, details)
  phi <- core$phi
  dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id"))
  cp04_finite(phi, "complete Corollary 1 OR link IF")
  if (scale == "core") phi <- st$lambda * phi

  out <- list(
    rows = rows,
    id = st$id[rows],
    links = links,
    phi = phi,
    scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(
      n_rate_score_evaluations = sum(pre$uids <= backend$upstream$n_a),
      n_selection_score_evaluations = length(up_res$active_cells),
      n_selection_block_evaluations = length(up_res$active_blocks),
      n_upstream_coordinates = length(pre$uids),
      n_regression_score_evaluations = length(union(pre$gids, pre$beta_gids)),
      n_global_beta_blocks = length(pre$beta_gids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  if (details) {
    tau <- eta <- if (is.null(ref)) matrix(0, B, length(st$cohorts)) else NULL
    if (is.null(ref)) for (k in seq_along(st$cohorts)) {
      tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
      eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
    }
    out$details <- list(
      support_ids = pre$sids,
      upstream_ids = pre$uids,
      upstream = up,
      beta_regression_ids = pre$beta_gids,
      beta = core$beta,
      beta_score = core$beta_score,
      beta_block_scores = core$beta_block_scores,
      regression_ids = pre$gids,
      alpha = core$alpha,
      alpha_score = core$alpha_score,
      direct = core$direct,
      denominator = core$denominator,
      numerator = core$numerator,
      design = core$design,
      objective = core$objective,
      tau = tau,
      eta = eta,
      lambda = A - st$lambda
    )
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    or_corollary1_abort("IF output bound underestimated.", "bound")
  out
}

# Summarize active selection cells and parameter counts for public records
selection_backend_summary <- function(backend) {
  if (is.null(backend)) return(NULL)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  if (isTRUE(backend$uses_selection_models) && !is.null(backend$selection_backend)) {
    sb <- backend$selection_backend
    cells_df <- if (!length(sb$cells)) {
      data.frame(cell_id = integer(), kind = character(), cohort = numeric(),
                 label = character(), block_id = integer(), stringsAsFactors = FALSE)
    } else {
      data.frame(
        cell_id = seq_along(sb$cells),
        kind = vapply(sb$cells, `[[`, character(1), "type"),
        cohort = st$cohorts[vapply(sb$cells, `[[`, integer(1), "cohort")],
        label = vapply(sb$cells, `[[`, character(1), "label"),
        block_id = vapply(sb$cells, `[[`, integer(1), "block_id"),
        stringsAsFactors = FALSE
      )
    }
    return(list(
      selection_model = sb$selection_model,
      selection_pooling = sb$selection_pooling,
      n_active_cells = length(sb$cells),
      n_selection_cells = length(sb$cells),
      n_selection_blocks = length(sb$blocks),
      n_selection_parameters = length(sb$coefficients),
      has_intercept = sb$design_spec$has_intercept,
      n_slope_columns = sb$design_spec$d_s,
      active_cells = cells_df,
      empirical_repeat_rates = if (is.null(sb$reference_ref)) stats::setNames(as.double(st$a), as.character(st$cohorts)) else
        chained_reference_record(chained_backend_reference(sb))$recipe_repeat_rates,
      repeat_rates = if (is.null(sb$reference_ref)) stats::setNames(as.double(st$a), as.character(st$cohorts)) else
        chained_reference_record(chained_backend_reference(sb))$recipe_repeat_rates
    ))
  }
  # Legacy separate pair_logit path
  n_cells <- length(backend$selection_specs)
  d_col <- ncol(st$Z)
  cells_df <- if (!n_cells) {
    data.frame(cell_id = integer(), kind = character(), cohort = numeric(),
               label = character(), block_id = integer(), stringsAsFactors = FALSE)
  } else {
    data.frame(
      cell_id = seq_len(n_cells),
      kind = rep.int("pair", n_cells),
      cohort = st$cohorts[vapply(backend$selection_specs, `[[`, integer(1), "cohort")],
      label = vapply(seq_len(n_cells), function(j) {
        sp <- backend$selection_specs[[j]]
        paste0("pair_", sp$pair$s, "_", sp$pair$t, "_g", st$cohorts[sp$cohort])
      }, character(1)),
      block_id = seq_len(n_cells),
      stringsAsFactors = FALSE
    )
  }
  list(
    selection_model = if (!is.null(backend$selection_model)) backend$selection_model else "pair_logit",
    selection_pooling = if (!is.null(backend$selection_pooling)) backend$selection_pooling else "separate",
    n_active_cells = n_cells,
    n_selection_cells = n_cells,
    n_selection_blocks = n_cells,
    n_selection_parameters = n_cells * d_col,
    has_intercept = !identical(attr(st$Z, "has_intercept"), FALSE) && !identical(attr(st$Z, "intercept"), FALSE) && d_col >= 1L,
    n_slope_columns = max(0L, d_col - 1L),
    active_cells = cells_df,
    empirical_repeat_rates = stats::setNames(as.double(st$a), as.character(st$cohorts)),
    repeat_rates = stats::setNames(as.double(st$a), as.character(st$cohorts))
  )
}

# Pure-helper table entry points for selection-model IPW and OR consumers.
# Exercise the exact same upstream map, log-weight derivative, forcing, link,
# WLS and IF core helpers as fit_cp04_selection_models / fit_or_selection_models
# without fabricating a fake C owner or fake namespace when cdid is not installed.
predict_table_pair_selection <- function(table_fit, Z, w, rho, rows, q) {
  pc <- table_fit$pair_cohorts[[q]]
  d_s <- table_fit$design_spec$d_s
  slope_cols <- table_fit$design_spec$slope_cols
  Z_slope <- if (d_s > 0L) Z[rows, slope_cols, drop = FALSE] else matrix(0, length(rows), 0L)
  eval_c <- function(cid) {
    cell <- table_fit$cells[[cid]]
    bl <- table_fit$blocks[[cell$block_id]]
    b <- cell$local_b
    if (identical(bl$branch, "zero_col")) return(rep.int(0.5, length(rows)))
    eta <- (if (bl$has_intercept) bl$alpha[b] else 0) +
      (if (bl$d_s > 0L) drop(Z_slope %*% bl$beta) else rep.int(0, length(rows)))
    ipw_clip_logit(eta)$prob
  }
  if (identical(table_fit$selection_model, "independent_period_logit")) {
    ps <- eval_c(pc$s_cell_id)
    pt <- eval_c(pc$t_cell_id)
    pp <- ps * pt
    w_sel <- w[rows] * table_fit$a[pc$cohort] / (rho[rows] * ipw_clip_floor(pp))
    list(prob = pp, prob_s = ps, prob_t = pt, prob_weight = ipw_clip_floor(pp),
         active = ipw_clip_floor_active(pp), weight = w_sel, a = table_fit$a[pc$cohort], pair_cohort_id = q)
  } else {
    pp <- eval_c(pc$pair_cell_id)
    w_sel <- w[rows] * table_fit$a[pc$cohort] / (rho[rows] * ipw_clip_floor(pp))
    list(prob = pp, prob_weight = ipw_clip_floor(pp), active = ipw_clip_floor_active(pp),
         weight = w_sel, a = table_fit$a[pc$cohort], pair_cohort_id = q)
  }
}

selection_table_upstream_if <- function(sel_fit, up_map) {
  N <- nrow(sel_fit$a_if)
  up_all <- matrix(0, N, up_map$D)
  if (up_map$n_a > 0L) {
    up_all[, seq_len(up_map$n_a)] <- sel_fit$a_if[, up_map$cohorts_a, drop = FALSE]
  }
  for (m in seq_along(sel_fit$blocks)) {
    bl <- sel_fit$blocks[[m]]
    tc <- up_map$theta_coords[[m]]
    if (!length(tc)) next
    xi_al <- if (bl$has_intercept) sel_fit$xi_alpha[, bl$cell_ids, drop = FALSE] else matrix(0, N, 0L)
    xi_be <- sel_fit$xi_beta_by_block[[m]]
    up_all[, tc] <- cbind(xi_al, xi_be)
  }
  up_all
}

fit_cp04_selection_table <- function(Z, V, Y, obs, g, w, U = rep(TRUE, nrow(Z)), rho = rep(1, nrow(Z)),
                                     universe = NULL, calendar, links, selection_model = "pair_logit",
                                     selection_pooling = "separate", has_intercept = NULL,
                                     glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
                                     chunk_rows = nrow(Z)) {
  N <- nrow(Z)
  cohorts <- sort(unique(as.double(g)), method = "radix")
  ci <- match(as.double(g), cohorts)
  A <- as.double(rowSums(obs != 0) >= 2L)
  if (is.null(universe)) {
    universe <- if (all(U)) "F" else if (identical(as.logical(U), as.logical(A > 0))) "A" else "B"
  }
  pair_specs <- list()
  ts <- list()
  for (ell in seq_along(links)) {
    lk <- links[[ell]]
    ks <- sort(as.double(lk$K), method = "radix")
    coh_idx <- match(c(as.double(lk$g), ks), cohorts)
    s_cal <- match(as.double(lk$s), calendar)
    t_cal <- match(as.double(lk$t), calendar)
    pair <- list(s = as.double(lk$s), t = as.double(lk$t), s_calendar_index = s_cal, t_calendar_index = t_cal)
    si <- integer(length(coh_idx))
    for (k in seq_along(coh_idx)) {
      sp <- list(pair = pair, cohort = coh_idx[k])
      same <- which(vapply(pair_specs, identical, logical(1), sp))
      if (length(same)) si[k] <- same[1L] else {
        pair_specs[[length(pair_specs) + 1L]] <- sp
        si[k] <- length(pair_specs)
      }
    }
    ts[[ell]] <- list(pair = pair, cohorts = coh_idx, selection_ids = si, link_id = lk$link_id)
  }
  pair_responses <- list()
  for (q in seq_along(pair_specs)) {
    pr <- pair_specs[[q]]$pair
    key <- paste0(pr$s_calendar_index, ":", pr$t_calendar_index)
    if (is.null(pair_responses[[key]])) {
      pair_responses[[key]] <- as.double(obs[, pr$s_calendar_index] != 0 & obs[, pr$t_calendar_index] != 0)
    }
  }
  sel_fit <- evaluate_selection_table_if(
    Z = Z, cohort_index = as.double(g), w = w, A = A, U = U, rho = rho,
    universe = universe, cohorts = cohorts, all_times = calendar,
    date_responses = obs, pair_responses = pair_responses,
    pair_specs = pair_specs, selection_model = selection_model,
    selection_pooling = selection_pooling, glm_control = glm_control,
    has_intercept = has_intercept, chunk_rows = chunk_rows,
    scale = "frame", details = TRUE
  )
  class(sel_fit) <- c("cdid_selection_models_table", class(sel_fit))
  for (ell in seq_along(ts)) {
    ts[[ell]]$selection_ids <- vapply(ts[[ell]]$cohorts, function(k) {
      resolve_pair_cohort_id(sel_fit, list(pair = ts[[ell]]$pair, cohort = k))
    }, integer(1))
  }
  up_map <- build_selection_upstream_map(sel_fit)
  up_all <- selection_table_upstream_if(sel_fit, up_map)
  p <- ncol(V)
  L <- length(ts)
  phi <- matrix(0, N, L, dimnames = list(NULL, vapply(ts, `[[`, character(1), "link_id")))
  treatments <- vector("list", L)
  for (ell in seq_len(L)) {
    spec <- ts[[ell]]
    s_idx <- spec$pair$s_calendar_index
    t_idx <- spec$pair$t_calendar_index
    S_pair <- as.double(obs[, s_idx] != 0 & obs[, t_idx] != 0)
    obs_rows <- which(ci %in% spec$cohorts & S_pair == 1)
    ci_obs <- match(ci[obs_rows], spec$cohorts)
    X <- V[obs_rows, , drop = FALSE]
    y <- as.double(ci_obs == 1L)
    dy <- Y[obs_rows, t_idx] - Y[obs_rows, s_idx]
    fw <- numeric(length(obs_rows))
    m <- length(spec$cohorts)
    ch_by_k <- vector("list", m)
    ld_by_k <- vector("list", m)
    for (k in seq_len(m)) {
      jj <- which(ci_obs == k)
      q <- spec$selection_ids[k]
      pred <- predict_table_pair_selection(sel_fit, Z, w, rho, obs_rows[jj], q)
      ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, q, obs_rows[jj], pred = pred, Z = Z, a_vec = sel_fit$a)
      fw[jj] <- pred$weight
      ch_by_k[[k]] <- up_map$pair_channels[[q]]
      ld_by_k[[k]] <- ld$log_deriv
    }
    uids <- sort(unique(unlist(ch_by_k)))
    u <- length(uids)
    fit <- cp04_logit(X, y, fw, validate_ipw_glm_control(glm_control), N, spec)
    lv <- ipw_logit_values(X, fit$coefficients)
    wp <- fw * y
    wm <- fw * (1 - y) * lv$odds
    dp <- sum(wp) / N
    dm <- sum(wm) / N
    T_mean <- sum(wp * dy) / sum(wp)
    M_mean <- sum(wm * dy) / sum(wm)
    C <- matrix(0, p, u)
    U_part <- matrix(0, 4L, u, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    score <- fw * (y - lv$prob) * X
    residual <- cbind(wp * (dy - T_mean) / dp, wm * (dy - M_mean) / dm, wp, wm)
    for (k in seq_len(m)) {
      jj <- which(ci_obs == k)
      pos_u <- match(ch_by_k[[k]], uids)
      ld <- ld_by_k[[k]]
      C[, pos_u] <- C[, pos_u] + crossprod(score[jj, , drop = FALSE], ld) / N
      U_part[, pos_u] <- U_part[, pos_u] + crossprod(residual[jj, , drop = FALSE], ld) / N
    }
    V_part <- matrix(0, 4L, p, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    wm_act <- wm * lv$active
    V_part[2L, ] <- colSums(wm_act * (dy - M_mean) * X) / (N * dm)
    V_part[4L, ] <- colSums(wm_act * X) / N
    beta_score <- matrix(0, N, p)
    direct <- matrix(0, N, 4L)
    direct[, 3L] <- -dp
    direct[, 4L] <- -dm
    beta_score[obs_rows, ] <- score
    direct[obs_rows, 1L] <- residual[, 1L]
    direct[obs_rows, 2L] <- residual[, 2L]
    direct[obs_rows, 3L] <- direct[obs_rows, 3L] + wp
    direct[obs_rows, 4L] <- direct[obs_rows, 4L] + wm
    up_sub <- up_all[, uids, drop = FALSE]
    xb <- t(cp04_solve(fit$H_chol, t(beta_score + up_sub %*% t(C))))
    coord <- direct + up_sub %*% t(U_part) + xb %*% t(V_part)
    colnames(coord) <- c("T", "M", "dp", "dm")
    phi[, ell] <- coord[, 1L] - coord[, 2L]
    fit$means <- c(T = T_mean, M = M_mean, dp = dp, dm = dm, mu = T_mean - M_mean)
    fit$upstream_ids <- uids
    fit$cross <- C
    fit$partials <- list(upstream = U_part, beta = V_part)
    fit$beta_if <- xb
    fit$coordinates_if <- coord
    treatments[[ell]] <- fit
  }
  list(
    selection = sel_fit,
    up_map = up_map,
    upstream_if = up_all,
    treatments = treatments,
    estimates = vapply(treatments, function(f) unname(f$means["T"] - f$means["M"]), numeric(1)),
    phi = phi
  )
}

fit_or_selection_table <- function(Z, X_or, Y, obs, g, w, U = rep(TRUE, nrow(Z)), rho = rep(1, nrow(Z)),
                                   universe = NULL, calendar, links, outcome_recipe, selection_model = "pair_logit",
                                   selection_pooling = "separate", has_intercept = NULL,
                                   glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
                                   rank_tol = 1e-10, chunk_rows = nrow(Z)) {
  N <- nrow(Z)
  cohorts <- sort(unique(as.double(g)), method = "radix")
  ci <- match(as.double(g), cohorts)
  A <- as.double(rowSums(obs != 0) >= 2L)
  if (is.null(universe)) {
    universe <- if (all(U)) "F" else if (identical(as.logical(U), as.logical(A > 0))) "A" else "B"
  }
  recipe <- or_pooled_check_recipe(outcome_recipe, colnames(X_or), "corollary1")
  records <- vector("list", length(links))
  for (ell in seq_along(links)) {
    lk <- links[[ell]]
    s_cal <- match(as.double(lk$s), calendar)
    t_cal <- match(as.double(lk$t), calendar)
    pair <- list(s = as.double(lk$s), t = as.double(lk$t), s_calendar_index = s_cal, t_calendar_index = t_cal)
    records[[ell]] <- list(
      pair = pair,
      target = match(as.double(lk$g), cohorts),
      controls = match(sort(as.double(lk$K), method = "radix"), cohorts),
      link_id = lk$link_id,
      source_functional_key = paste0("cmpcat:", lk$link_id),
      functional_key = paste0("cmpcat:OR-Cor1:", lk$link_id)
    )
  }
  catalog <- or_baseline_catalog(records, cohorts, recipe)
  specs <- catalog$specs
  pools <- catalog$pools
  cat_links <- catalog$links
  selection_specs <- lapply(specs, function(s) list(pair = s$pair, cohort = s$cohort))
  S <- length(specs)
  S_mat <- matrix(0, N, S)
  pair_responses <- list()
  for (s in seq_len(S)) {
    pr <- specs[[s]]$pair
    S_mat[, s] <- as.double(obs[, pr$s_calendar_index] != 0 & obs[, pr$t_calendar_index] != 0)
    key <- paste0(pr$s_calendar_index, ":", pr$t_calendar_index)
    if (is.null(pair_responses[[key]])) pair_responses[[key]] <- S_mat[, s]
  }
  sel_fit <- evaluate_selection_table_if(
    Z = Z, cohort_index = as.double(g), w = w, A = A, U = U, rho = rho,
    universe = universe, cohorts = cohorts, all_times = calendar,
    date_responses = obs, pair_responses = pair_responses,
    pair_specs = selection_specs, selection_model = selection_model,
    selection_pooling = selection_pooling, glm_control = glm_control,
    has_intercept = has_intercept, chunk_rows = chunk_rows,
    scale = "frame", details = TRUE
  )
  class(sel_fit) <- c("cdid_selection_models_table", class(sel_fit))
  cohorts_a <- sort(unique(vapply(specs, `[[`, integer(1), "cohort")))
  up_map <- build_selection_upstream_map(sel_fit, cohorts_a = cohorts_a)
  D <- up_map$D
  q_of_s <- vapply(specs, function(s) {
    resolve_pair_cohort_id(sel_fit, list(pair = s$pair, cohort = s$cohort))
  }, integer(1))
  channels <- lapply(seq_len(S), function(s) up_map$pair_channels[[q_of_s[s]]])
  up_all <- selection_table_upstream_if(sel_fit, up_map)

  supports <- vector("list", S)
  for (s in seq_len(S)) {
    sp <- specs[[s]]
    q <- q_of_s[s]
    s_idx <- sp$pair$s_calendar_index
    t_idx <- sp$pair$t_calendar_index
    obs_rows <- which(ci == sp$cohort & S_mat[, s] == 1)
    dy <- Y[obs_rows, t_idx] - Y[obs_rows, s_idx]
    Z_obs <- X_or[obs_rows, , drop = FALSE]
    pred <- predict_table_pair_selection(sel_fit, Z, w, rho, obs_rows, q)
    ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, q, obs_rows, pred = pred, Z = Z, a_vec = sel_fit$a)
    acc <- or_pooled_support_start(recipe, sp$transitions)
    acc <- or_pooled_support_add(acc, recipe, Z_obs, dy, pred$weight, "corollary1")
    m_dim <- recipe$p + 2L
    moments <- array(0, c(m_dim, m_dim, ncol(ld$T_logit)))
    moments <- or_selection_moments_add(
      moments,
      Z_obs[, recipe$slope_index, drop = FALSE],
      dy,
      pred$weight,
      ld$T_logit
    )
    supports[[s]] <- list(
      design = acc,
      moments = moments,
      a = sel_fit$a[sp$cohort],
      pair_mass = acc$weight_mass,
      n_pair = acc$n
    )
  }
  state <- or_selection_assemble(
    recipe, supports, pools, cat_links, channels, D, rank_tol, N,
    function(sid, visit) {
      sp <- specs[[sid]]
      obs_rows <- which(ci == sp$cohort & S_mat[, sid] == 1)
      if (length(obs_rows)) visit(X_or[obs_rows, , drop = FALSE])
    }
  )
  fake_backend <- list(
    outcome_recipe = recipe,
    support_specs = specs,
    supports = state$supports,
    regression_specs = pools,
    regressions = state$regressions,
    shared = state$shared,
    links = state$links,
    upstream = list(n_a = up_map$n_a, d = up_map$d, D = D)
  )
  L <- length(state$links)
  sets <- or_pooled_if_sets(fake_backend, seq_len(L), details = TRUE)
  uids <- sort(unique(unlist(lapply(state$supports[sets$sids], `[[`, "channels"))))
  pre <- list(gids = sets$gids, beta_gids = sets$beta_gids, sids = sets$sids, uids = uids)
  kept <- recipe$slope_index[state$shared$keep]
  hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]
    sp <- specs[[sid]]
    q <- q_of_s[sid]
    s_idx <- sp$pair$s_calendar_index
    t_idx <- sp$pair$t_calendar_index
    obs_rows <- which(ci == sp$cohort & S_mat[, sid] == 1)
    pred <- predict_table_pair_selection(sel_fit, Z, w, rho, obs_rows, q)
    hits[[j]] <- list(
      obs = obs_rows,
      dy = Y[obs_rows, t_idx] - Y[obs_rows, s_idx],
      z = X_or[obs_rows, kept, drop = FALSE],
      v = or_corollary1_factor(sel_fit$a[sp$cohort], rho[obs_rows], pred$prob)
    )
  }
  core <- or_pooled_if_core(fake_backend, N, w, hits, up_all[, uids, drop = FALSE], pre, seq_len(L), details = TRUE)
  phi <- core$phi
  dimnames(phi) <- list(NULL, vapply(state$links, `[[`, character(1), "link_id"))
  list(
    selection = sel_fit,
    up_map = up_map,
    upstream_if = up_all,
    shared = state$shared,
    regressions = state$regressions,
    links = state$links,
    estimates = vapply(state$links, `[[`, numeric(1), "mu"),
    phi = phi,
    beta_if = core$beta,
    alpha_if = core$alpha
  )
}
