# Private Chained Doubly Robust (DR) backend (CP07; DR_SPEC.md and BACKEND_REUSE_NOTES.md).
# Connects:
#   1. Baseline MTAR empirical pair rates or conditional_mtar selection models (fitted ONCE)
#   2. Shared one-period outcome-change control regression over deduplicated blocks (fitted ONCE)
#   3. Distinct risk-set treatment propensity logits (fitted ONCE per unique (pair, g, K))
#   4. Complete original-unit four-term influence function:
#      psi_DR = psi_dir + psi_or + psi_ps + psi_sel
#      with raw core = (n_A / N) * psi_DR scaling and bounded row/link chunking.

dr_baseline_abort <- function(message, class = "input", extra = list()) {
  stop(structure(
    c(list(message = message, call = NULL), extra),
    class = unique(c(
      paste0("cdid_dr_baseline_", class, "_error"),
      paste0("cdid_or_baseline_", class, "_error"),
      paste0("cdid_cp04_", class, "_error"),
      "cdid_ipw_fit_error", "cdid_error", "error", "condition"
    ))
  ))
}

dr_corollary1_abort <- function(message, class = "input", extra = list()) {
  stop(structure(
    c(list(message = message, call = NULL), extra),
    class = unique(c(
      paste0("cdid_dr_corollary1_", class, "_error"),
      paste0("cdid_or_corollary1_", class, "_error"),
      paste0("cdid_cp04_", class, "_error"),
      "cdid_ipw_fit_error", "cdid_error", "error", "condition"
    ))
  ))
}

dr_route_abort <- function(route, message, class = "input", extra = list()) {
  if (identical(route, "baseline")) {
    dr_baseline_abort(message, class = class, extra = extra)
  } else {
    dr_corollary1_abort(message, class = class, extra = extra)
  }
}

# Validate the private OR recipe against the common frame Z_a columns.
# Allows frame$Z_a to hold additional columns used by the treatment design while
# requiring recipe$columns to be an exact ordered subset of colnames(frame$Z_a).
dr_check_outcome_recipe <- function(recipe, frame_columns, route = "baseline") {
  ok <- inherits(recipe, "cdid_or_pooled_recipe") && is.list(recipe) &&
    identical(
      tryCatch(
        or_pooled_outcome_recipe(recipe$columns, recipe$intercept, recipe$slope_columns),
        error = function(e) NULL
      ),
      recipe
    ) &&
    is.character(frame_columns) && length(frame_columns) >= 1L &&
    all(recipe$columns %in% frame_columns) && !anyDuplicated(recipe$columns)
  if (!isTRUE(ok)) {
    dr_route_abort(
      route,
      "The outcome recipe must be a valid `cdid_or_pooled_recipe` whose columns belong to `frame$Z_a`.",
      "input"
    )
  }
  frame_index <- match(recipe$columns, frame_columns)
  list(recipe = recipe, frame_index = as.integer(frame_index))
}

# Validate the treatment propensity design column selection on frame$Z_a.
# Defaults to all columns of frame$Z_a when treatment_columns is NULL.
dr_check_treatment_recipe <- function(treatment_columns, frame_columns, route = "baseline") {
  if (!is.character(frame_columns) || !length(frame_columns) || anyNA(frame_columns) ||
      any(!nzchar(frame_columns)) || anyDuplicated(frame_columns)) {
    dr_route_abort(route, "Common frame `Z_a` requires distinct nonempty column names.", "input")
  }
  if (is.null(treatment_columns)) {
    idx <- seq_along(frame_columns)
  } else if (is.character(treatment_columns) && !is.object(treatment_columns) && is.null(dim(treatment_columns))) {
    if (!length(treatment_columns) || anyNA(treatment_columns) || any(!nzchar(treatment_columns)) ||
        anyDuplicated(treatment_columns) || !all(treatment_columns %in% frame_columns)) {
      dr_route_abort(route, "`treatment_columns` must be a nonempty subset of distinct `frame$Z_a` column names.", "input")
    }
    idx <- match(treatment_columns, frame_columns)
  } else if (is.numeric(treatment_columns) && !is.object(treatment_columns) && is.null(dim(treatment_columns))) {
    if (!length(treatment_columns) || anyNA(treatment_columns) || any(!is.finite(treatment_columns)) ||
        any(treatment_columns != floor(treatment_columns)) ||
        any(treatment_columns < 1L | treatment_columns > length(frame_columns)) ||
        anyDuplicated(as.integer(treatment_columns))) {
      dr_route_abort(route, "`treatment_columns` integer indices must be distinct positions within `1..ncol(frame$Z_a)`.", "input")
    }
    idx <- as.integer(treatment_columns)
  } else {
    dr_route_abort(route, "`treatment_columns` must be NULL, character column names, or integer column indices.", "input")
  }
  structure(
    list(
      columns = frame_columns[idx],
      index = as.integer(idx),
      p = length(idx)
    ),
    class = "cdid_dr_treatment_recipe"
  )
}

# Bounded cohort chunks of ONE actual (pair, cohort) support, projecting Z_a onto
# the OR recipe columns `or_frame_index`. A D60 support reads its recipe's
# frozen design (x role, same column layout) instead of Z_a. `complete` (CP08
# stage 3) reads the same pair rows by complete-date position.
dr_support_visit <- function(reader, spec, ctrl, or_frame_index, visit, reference = NULL, complete = NULL) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- if (is.null(complete)) NULL else or_baseline_complete_panel(reader, complete, dr_baseline_abort)
  first <- st$cohort_start[spec$cohort]
  last <- st$cohort_end[spec$cohort]
  for (start in seq.int(first, last, by = ctrl$max_chunk_rows)) {
    rows <- st$cohort_order[seq.int(start, min(last, start + ctrl$max_chunk_rows - 1L))]
    hit <- if (is.null(panel)) cp04_pair_rows(reader, rows, spec$pair, TRUE) else
      or_baseline_complete_pair_rows(reader, panel, rows, spec$pair, TRUE, dr_baseline_abort)
    obs <- which(hit$S == 1)
    Z_sub <- if (is.null(reference)) frame$Z_a[st$a_index[rows[obs]], or_frame_index, drop = FALSE] else
      chained_reference_x(reference, spec$reference, rows[obs])[, or_frame_index, drop = FALSE]
    visit(rows, hit$S, obs, hit$dy[obs], Z_sub)
  }
  invisible(NULL)
}

# A D60 support reads its context weight W_g on A rows only.
dr_baseline_support <- function(reader, spec, recipe, or_frame_index, ctrl, reference = NULL, complete = NULL) {
  st <- reader$state
  acc <- or_pooled_support_start(recipe, spec$transitions)
  repeat_mass <- 0
  n_repeat <- 0
  dr_support_visit(reader, spec, ctrl, or_frame_index, function(rows, S, obs, dy, Z) {
    A <- st$a_index[rows] > 0L
    w <- if (is.null(reference)) st$w[rows] else {
      x <- rep(NA_real_, length(rows))
      x[A] <- chained_reference_weight(reference, spec$reference, rows[A])
      x
    }
    repeat_mass <<- repeat_mass + sum(w[A])
    n_repeat <<- n_repeat + sum(A)
    acc <<- or_pooled_support_add(acc, recipe, Z, dy, w[obs], "baseline")
  }, reference, complete)
  or_baseline_support_rate(acc, repeat_mass, n_repeat, st$N)
}

dr_selection_support <- function(reader, spec, sel_fit, up_map, q_id, recipe, or_frame_index, ctrl, reference = NULL) {
  st <- reader$state
  a <- if (is.null(reference)) st$a[spec$cohort] else chained_reference_stat(reference, spec$reference, spec$cohort, "a")
  m <- recipe$p + 2L
  dq <- length(up_map$pair_channels[[q_id]]) - 1L
  acc <- or_pooled_support_start(recipe, spec$transitions)
  moments <- array(0, c(m, m, dq))

  dr_support_visit(reader, spec, ctrl, or_frame_index, function(rows, S, obs, dy, Z) {
    if (!length(obs)) return(invisible(NULL))
    rr <- rows[obs]
    ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, q_id, rr)
    q <- ld$pred$weight
    cp04_finite(q, "Corollary 1 fitted pair weights", TRUE)
    acc <<- or_pooled_support_add(acc, recipe, Z, dy, q, "corollary1")
    moments <<- or_selection_moments_add(moments, Z[, recipe$slope_index, drop = FALSE], dy, q, ld$T_logit)
  }, reference)
  if (!is.finite(acc$weight_mass) || acc$weight_mass <= 0 || acc$weight_mass / st$N <= 0)
    dr_corollary1_abort("Positive actual-pair cohort support required.", "support")
  cp04_finite(c(acc$sumZ, acc$sumY, moments), "Corollary 1 compact support moments")
  list(design = acc, moments = moments, a = a, pair_mass = acc$weight_mass, n_pair = acc$n)
}

# Build deduplicated treatment model specifications across links.
# Each treatment model is uniquely determined by (pair, target g, sorted controls K)
# and, on a D60 reference route, its recipe. `n_observers` (CP08 descriptor
# intake) supplies each link's observed risk-unit count when no link view
# exists; NULL reads it from the link views as before.
dr_build_treatment_catalog <- function(links, pools, specs, link_views, q_of_s = NULL, n_observers = NULL) {
  L <- length(links)
  if (!is.null(n_observers) && (!is.numeric(n_observers) || length(n_observers) != L || anyNA(n_observers)))
    dr_baseline_abort("Descriptor observer counts must cover every link.", "integrity")
  tkeys <- character()
  ts <- list()
  for (i in seq_len(L)) {
    l <- links[[i]]
    s_plus <- l$support_id
    s_minus <- pools[[l$regression]]$support_ids
    s_all <- c(s_plus, s_minus)
    cohorts <- c(l$cohort, pools[[l$regression]]$cohorts)
    pair <- specs[[s_plus]]$pair
    tkey <- paste(pair$s_calendar_index, pair$t_calendar_index, l$cohort,
                  paste(pools[[l$regression]]$cohorts, collapse = ","), sep = ":")
    if (!is.null(l$reference)) tkey <- paste0(tkey, "|recipe:", l$reference)
    tid <- match(tkey, tkeys)
    if (is.na(tid)) {
      tid <- length(ts) + 1L
      tkeys[tid] <- tkey
      ts[[tid]] <- list(
        pair = pair,
        pair_id = l$pair_id,
        cohorts = as.integer(cohorts),
        support_ids = as.integer(s_all),
        selection_ids = if (is.null(q_of_s)) NULL else as.integer(q_of_s[s_all]),
        regression = as.integer(l$regression),
        n_observers = if (is.null(n_observers)) as.integer(link_views[[i]]$support$n_observed_risk_units) else
          as.integer(n_observers[i])
      )
      if (!is.null(l$reference)) ts[[tid]]$reference <- as.integer(l$reference)
    }
    links[[i]]$model <- as.integer(tid)
    links[[i]]$control_support_ids <- as.integer(s_minus)
  }
  list(treatment_specs = ts, links = links)
}

# CP08 stage 3: the eligible route's one treatment logit. The first model's
# actual ordered inputs are fitted; every later model must present exactly the
# same ordered observer rows, design, binary response, fitting weights, support
# factors v, cohorts, GLM controls and normalization, and then reuses that fit
# and its fitted values (a difference refuses; nothing is refitted elsewhere).
dr_baseline_shared_logit <- function(share, spec, rows, X, y, w, v, glm_control, N) {
  if (is.null(share$fit)) {
    share$fit <- cp04_logit(X, y, w, glm_control, N, spec)
    share$inputs <- list(rows = rows, X = X, y = y, w = w, v = v, cohorts = spec$cohorts,
      glm_control = glm_control, N = N)
    share$calls <- 1L
    share$reuses <- 0L
  } else {
    inp <- share$inputs
    if (!identical(rows, inp$rows) || !identical(X, inp$X) || !identical(y, inp$y) ||
        !identical(w, inp$w) || !identical(v, inp$v) || !identical(spec$cohorts, inp$cohorts) ||
        !identical(glm_control, inp$glm_control) || !identical(N, inp$N))
      dr_baseline_abort("Complete-support treatment inputs differ from the shared fit's ordered inputs.", "identity")
    share$reuses <- share$reuses + 1L
  }
  share$fit
}

# CP08 construction: the outcome-independent invariants of an eligible fit's
# treatment models, kept in the existing `share` environment beside the one
# shared logit and released with it. The fit's gate admits the complete route
# only when every treatment model has identical cohorts and observer counts,
# so every model's ordered observer rows, treatment design, binary response
# and cohort positions are the same deterministic reads of the same sealed
# reader and frame, and its fitting weights are the same whenever its support
# factors v are the same; the first model's arrays ARE these objects (its
# `share$inputs` plus the derived `share$prep`). A later model borrows them
# only when its compact bindings are identical: cohorts, observer count,
# support factors v, GLM controls and normalization (the shared logit's own
# inputs), the treatment design columns, kept OR columns and frame reference.
# NULL means construct: the original per-model path runs and the shared
# logit's unchanged exact-input check then decides (a difference refuses).
# No numerical near-equality admits a borrow; nothing is refitted.
dr_baseline_shared_prepared <- function(share, spec, supports, reader, trt_recipe, kept_or_cols, p_or,
                                        glm_control, N) {
  if (!is.environment(share)) return(NULL)
  prep <- share$prep; inp <- share$inputs
  if (!is.list(prep) || !is.list(inp) || is.null(share$fit) || !is.list(spec) ||
      !is.list(supports) || !is.list(reader) || !is.list(trt_recipe)) return(NULL)
  n <- spec$n_observers
  sids <- spec$support_ids
  if (!is.integer(sids) || length(sids) != length(spec$cohorts) || anyNA(sids) ||
      any(sids < 1L | sids > length(supports))) return(NULL)
  v <- vapply(supports[sids], function(u) {
    if (is.list(u) && is.double(u$v) && length(u$v) == 1L && is.finite(u$v)) u$v else NA_real_
  }, numeric(1))
  if (anyNA(v)) return(NULL)
  frame <- if (is.environment(reader$frame_ref))
    get0("frame_owner", envir = reader$frame_ref, inherits = FALSE) else NULL
  if (!identical(spec$cohorts, inp$cohorts) || !identical(as.integer(n), length(inp$rows)) || !isTRUE(n > 0L) ||
      !identical(v, inp$v) || !identical(glm_control, inp$glm_control) || !identical(N, inp$N) ||
      !identical(prep$reader_ref, reader$ref) || !identical(prep$reader_state, reader$state) ||
      !identical(prep$frame_ref, reader$frame_ref) || !identical(prep$frame, frame) ||
      !identical(prep$trt_recipe, trt_recipe) || !identical(prep$trt_index, trt_recipe$index) ||
      !identical(prep$kept_or_cols, kept_or_cols) || !identical(prep$p_or, p_or) || !identical(prep$n, n))
    return(NULL)
  prep
}

# Fit ONE distinct baseline treatment propensity model (pair, g, K) and compute
# both residualized DR and raw IPW moments/partials without refitting the logit.
# `complete`/`share` (CP08 stage 3): complete-date positions and the one shared
# logit; every pair-specific moment, rate derivative and partial is its own.
dr_baseline_treatment_fit <- function(reader, spec, supports, reg, shared,
                                      recipe, or_frame_index, trt_recipe,
                                      glm_control, ctrl, reference = NULL,
                                      complete = NULL, share = NULL) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- if (is.null(complete)) NULL else or_baseline_complete_panel(reader, complete, dr_baseline_abort)
  p_trt <- trt_recipe$p
  p_or <- shared$p
  kept_or_cols <- if (p_or > 0L) or_frame_index[recipe$slope_index[shared$keep]] else integer(0)
  kept_columns <- if (p_or > 0L) shared$basis$kept_columns else NULL
  m <- length(spec$cohorts)
  n <- spec$n_observers

  # CP08 construction: a later model of an eligible fit borrows the first
  # model's prepared invariants under identical compact bindings
  # (dr_baseline_shared_prepared); otherwise the original arrays are built.
  prep <- if (is.null(panel) || is.null(share) || !is.null(reference)) NULL else
    dr_baseline_shared_prepared(share, spec, supports, reader, trt_recipe, kept_or_cols, p_or, glm_control, st$N)
  if (!is.null(prep) && !identical(prep$kept_columns, kept_columns)) prep <- NULL
  if (!is.null(prep)) {
    share$prep_reuses <- (if (is.null(share$prep_reuses)) 0L else share$prep_reuses) + 1L
    rows <- obs_rows <- share$inputs$rows
    X_trt <- share$inputs$X; y <- share$inputs$y; w <- share$inputs$w
    Z_or <- prep$Z_or; kk <- prep$kk
    dy <- numeric(n)
  } else {
    X_trt <- matrix(0, n, p_trt, dimnames = list(NULL, trt_recipe$columns))
    Z_or <- matrix(0, n, p_or, dimnames = if (p_or > 0L) list(NULL, kept_columns) else NULL)
    y <- w <- dy <- numeric(n)
    kk <- integer(n)
    obs_rows <- if (is.null(share)) NULL else integer(n)
  }
  pos <- 0L

  if (!is.null(panel)) {
    # CP08 stage 4: the eligible fit's observer rows in the scan's order (every
    # frame row of the model's cohorts; complete dates observe all), with one
    # native call reading their raw pair differences under the checks and
    # errors of the scan's per-chunk or_baseline_complete_pair_rows calls; the
    # design, response, cohort and weight vectors take the scan's values. The
    # general scan below runs only without complete-date facts.
    if (p_or < 0L || !is.null(reference) || is.null(obs_rows))
      dr_baseline_abort("Complete-support treatment reads require the facts of their fit.", "integrity")
    # The frame rows of the model's cohorts in increasing order, from the
    # reader's cohort ranges (no frame-length temporary); a borrowing model
    # reads its own pair differences on the shared rows.
    if (is.null(prep))
      rows <- sort(unlist(lapply(spec$cohorts, function(k) st$cohort_order[seq_len(st$cohort_end[k] -
        st$cohort_start[k] + 1) + st$cohort_start[k] - 1])), method = "radix")
    post <- panel$postings
    res <- .Call(or_baseline_native_kernel("C_gmm_cp08_pair_rows"), panel$observations$y, post$date_start,
      post$post_a_index, post$post_obs_index, rows, st$a_index,
      as.integer(c(spec$pair$s_calendar_index, spec$pair$t_calendar_index)), as.integer(ctrl$max_chunk_rows),
      ctrl$max_matrix_entries)
    if (!identical(res$status[1L], 0L)) or_baseline_native_refuse(res$status, dr_baseline_abort)
    pos <- length(rows)
    if (pos == n && n > 0L) {
      if (is.null(prep)) {
        ci <- match(st$cohort_index[rows], spec$cohorts)
        a_idx <- st$a_index[rows]
        obs_rows[] <- rows
        X_trt[] <- frame$Z_a[a_idx, trt_recipe$index, drop = FALSE]
        if (p_or > 0L) Z_or[] <- frame$Z_a[a_idx, kept_or_cols, drop = FALSE]
        y[] <- as.double(ci == 1L)
        kk[] <- ci
        w[] <- st$w[rows] * vapply(supports[spec$support_ids], `[[`, numeric(1), "v")[ci]
      }
      dy[] <- res$dy
    }
  }
  for (start in if (is.null(panel)) seq.int(1L, st$N, by = ctrl$max_chunk_rows)) {
    rows <- seq.int(start, min(st$N, start + ctrl$max_chunk_rows - 1L))
    ci <- match(st$cohort_index[rows], spec$cohorts)
    rows <- rows[!is.na(ci)]
    ci <- ci[!is.na(ci)]
    if (!length(rows)) next
    hit <- if (is.null(panel)) cp04_pair_rows(reader, rows, spec$pair, TRUE) else
      or_baseline_complete_pair_rows(reader, panel, rows, spec$pair, TRUE, dr_baseline_abort)
    keep <- which(hit$S == 1)
    if (!length(keep)) next
    rows <- rows[keep]
    ci <- ci[keep]
    ii <- seq.int(pos + 1L, pos + length(rows))
    pos <- pos + length(rows)
    if (!is.null(obs_rows)) obs_rows[ii] <- rows
    a_idx <- st$a_index[rows]
    if (is.null(reference)) {
      X_trt[ii, ] <- frame$Z_a[a_idx, trt_recipe$index, drop = FALSE]
      if (p_or > 0L) Z_or[ii, ] <- frame$Z_a[a_idx, kept_or_cols, drop = FALSE]
    } else {
      Xr <- chained_reference_x(reference, spec$reference, rows)
      X_trt[ii, ] <- Xr[, trt_recipe$index, drop = FALSE]
      if (p_or > 0L) Z_or[ii, ] <- Xr[, kept_or_cols, drop = FALSE]
    }
    y[ii] <- as.double(ci == 1L)
    kk[ii] <- ci
    dy[ii] <- hit$dy[keep]
    w_rows <- if (is.null(reference)) st$w[rows] else chained_reference_weight(reference, spec$reference, rows)
    for (k in unique(ci)) {
      jj <- which(ci == k)
      sid <- spec$support_ids[k]
      w[ii[jj]] <- w_rows[jj] * supports[[sid]]$v
    }
  }
  if (pos != n) dr_baseline_abort("Observed risk support changed from its owner view.", "identity")

  fit <- if (is.null(share)) cp04_logit(X_trt, y, w, glm_control, st$N, spec) else
    dr_baseline_shared_logit(share, spec, obs_rows, X_trt, y, w,
      vapply(supports[spec$support_ids], `[[`, numeric(1), "v"), glm_control, st$N)
  names(fit$coefficients) <- trt_recipe$columns
  lv <- if (is.null(share)) ipw_logit_values(X_trt, fit$coefficients) else {
    if (is.null(share$lv)) share$lv <- ipw_logit_values(X_trt, fit$coefficients)
    share$lv
  }
  # Outcome-independent weights, their sums and normalizers: the borrowed
  # objects of the first model, or this model's own (the original statements).
  if (!is.null(prep)) {
    wp <- prep$wp; wm <- prep$wm; sum_wp <- prep$sum_wp; sum_wm <- prep$sum_wm; dp <- prep$dp; dm <- prep$dm
  } else {
    wp <- w * y
    wm <- w * (1 - y) * lv$odds
    sum_wp <- sum(wp)
    sum_wm <- sum(wm)
    dp <- sum_wp / st$N
    dm <- sum_wm / st$N
  }
  cp04_finite(c(dp, dm), "treated/control normalizers", TRUE)

  # Raw outcome IPW quantities
  T_ipw <- sum(wp * dy) / sum_wp
  M_ipw <- sum(wm * dy) / sum_wm
  mu_ipw <- T_ipw - M_ipw

  # Residualized DR quantities reusing fitted OR parameters (no treatment refit)
  d_obs <- reg$L * Z_or
  m_hat <- (if (recipe$intercept) reg$alpha else 0) +
    (if (p_or > 0L) as.vector(d_obs %*% shared$beta) else 0)
  r_obs <- dy - m_hat
  T_dr <- sum(wp * r_obs) / sum_wp
  M_dr <- sum(wm * r_obs) / sum_wm
  mu_dr <- T_dr - M_dr
  d_bar_plus <- if (p_or > 0L) as.vector(colSums(wp * d_obs) / sum_wp) else numeric(0)
  d_bar_minus <- if (p_or > 0L) as.vector(colSums(wm * d_obs) / sum_wm) else numeric(0)
  d_bar_diff <- d_bar_plus - d_bar_minus
  if (p_or > 0L) {
    names(d_bar_plus) <- names(d_bar_minus) <- names(d_bar_diff) <- shared$basis$kept_columns
  }
  mu_dr_decomp <- mu_ipw - (if (p_or > 0L) sum(d_bar_diff * shared$beta) else 0)
  cp04_finite(c(T_ipw, M_ipw, mu_ipw, T_dr, M_dr, mu_dr, d_bar_plus, d_bar_minus, mu_dr_decomp),
              "treated/control DR and IPW means")

  # Upstream rate cross-derivatives on spec$support_ids (log_deriv = -1{rate_k > eps} / rate_k)
  C <- matrix(0, p_trt, m)
  U_dr <- matrix(0, 4L, m, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  U_ipw <- matrix(0, 4L, m, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  # The score's per-cohort column sums and the wm-weighted design sums are
  # outcome independent: borrowed, or computed once here by the original
  # statements; each model's residual records and rate scales stay its own.
  score <- if (is.null(prep)) w * (y - lv$prob) * X_trt else NULL
  cohort_rows <- if (is.null(prep)) lapply(seq_len(m), function(k) which(kk == k)) else prep$cohort_rows
  score_sums <- if (is.null(prep)) matrix(0, p_trt, m) else prep$score_sums
  residual_dr <- cbind(wp * (r_obs - T_dr) / dp, wm * (r_obs - M_dr) / dm, wp, wm)
  residual_ipw <- cbind(wp * (dy - T_ipw) / dp, wm * (dy - M_ipw) / dm, wp, wm)

  for (k in seq_len(m)) {
    rr <- cohort_rows[[k]]
    if (!length(rr)) next
    sid <- spec$support_ids[k]
    r_k <- supports[[sid]]$rate
    scale_k <- if (r_k > ipw_clip_eps) -1 / (st$N * r_k) else 0
    if (is.null(prep)) score_sums[, k] <- colSums(score[rr, , drop = FALSE])
    C[, k] <- score_sums[, k] * scale_k
    U_dr[, k] <- colSums(residual_dr[rr, , drop = FALSE]) * scale_k
    U_ipw[, k] <- colSums(residual_ipw[rr, , drop = FALSE]) * scale_k
  }

  # Odds derivatives carry the active-odds mask; the raw score does not.
  wm_act <- wm * lv$active
  wm_design <- if (is.null(prep)) colSums(wm_act * X_trt) else prep$wm_design
  V_dr <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), trt_recipe$columns))
  V_dr[2L, ] <- colSums(wm_act * (r_obs - M_dr) * X_trt) / (st$N * dm)
  V_dr[4L, ] <- wm_design / st$N

  V_ipw <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), trt_recipe$columns))
  V_ipw[2L, ] <- colSums(wm_act * (dy - M_ipw) * X_trt) / (st$N * dm)
  V_ipw[4L, ] <- V_dr[4L, ]

  cp04_finite(c(C, U_dr, U_ipw, V_dr, V_ipw), "baseline DR analytic cross derivatives")

  # The first model of an eligible fit (the one that fitted the shared logit)
  # keeps its outcome-independent invariants for the later models; the
  # arrays are its own objects, not copies. Released with `share` by the fit.
  if (!is.null(share) && !is.null(panel) && is.null(prep) && is.null(share$prep) &&
      identical(share$calls, 1L) && identical(share$reuses, 0L)) {
    share$prep <- list(reader_ref = reader$ref, reader_state = st, frame_ref = reader$frame_ref, frame = frame,
      trt_recipe = trt_recipe, trt_index = trt_recipe$index, kept_or_cols = kept_or_cols,
      kept_columns = kept_columns, p_or = p_or, n = n, kk = kk, Z_or = Z_or, wp = wp, wm = wm,
      sum_wp = sum_wp, sum_wm = sum_wm, dp = dp, dm = dm, cohort_rows = cohort_rows,
      score_sums = score_sums, wm_design = wm_design)
    share$prep_calls <- 1L
    share$prep_reuses <- 0L
  }

  fit$support_ids <- as.integer(spec$support_ids)
  fit$cross <- C
  fit$partials <- list(upstream = U_dr, beta = V_dr)
  fit$partials_ipw <- list(upstream = U_ipw, beta = V_ipw)
  fit$means <- c(
    T = T_dr, M = M_dr, dp = dp, dm = dm,
    T_dr = T_dr, M_dr = M_dr, mu_dr = mu_dr,
    T_ipw = T_ipw, M_ipw = M_ipw, mu_ipw = mu_ipw,
    mu_or = T_dr, mu_dr_decomp = mu_dr_decomp
  )
  fit$d_bar_plus <- d_bar_plus
  fit$d_bar_minus <- d_bar_minus
  fit$d_bar_diff <- d_bar_diff
  fit$clipping <- ipw_clip_summary(prob = lv$prob,
    rates = vapply(supports[spec$support_ids], `[[`, numeric(1), "rate"), controls = y == 0, weights = w)
  fit
}

# CP08 stage 7: complete-roster DR treatment-summary batch for the remaining
# models after the first model fits the shared logit and records share$prep.
# Requires exact immutable bindings (dr_baseline_shared_prepared) for every
# model in `specs`, reserves simultaneous scratch/output/work under `ctrl`
# before allocating descriptors or calling C, and evaluates every model's raw
# differences (subtracted before summation in original observer/chunk order),
# IPW/DR means, covariate contrasts and analytic cross derivatives in one
# bounded native pass without materializing N-length per-model temporaries.
# Returns NULL (declining to the ordinary per-model path) if any structural
# or binding invariant is not satisfied.
dr_baseline_treatment_batch <- function(reader, specs, supports, regressions, shared,
                                        recipe, or_frame_index, trt_recipe,
                                        glm_control, ctrl, reference = NULL,
                                        complete = NULL, share = NULL,
                                        scratch = NULL, call_fn = .Call) {
  if (is.null(complete) || !is.null(reference) || !is.environment(share) ||
      !is.list(specs) || !length(specs) || !is.list(supports) ||
      !is.list(regressions) || !is.list(shared) || !is.list(reader))
    return(NULL)
  if (!identical(share$calls, 1L) || !identical(share$reuses, 0L) ||
      !identical(share$prep_calls, 1L) || !identical(share$prep_reuses, 0L))
    return(NULL)
  panel <- or_baseline_complete_panel(reader, complete, dr_baseline_abort)
  st <- reader$state
  frame <- if (is.environment(reader$frame_ref))
    get0("frame_owner", envir = reader$frame_ref, inherits = FALSE) else NULL
  if (is.null(frame) || !is.list(st)) return(NULL)
  if (!inherits(recipe, "cdid_or_pooled_recipe") || !is.logical(recipe$intercept) ||
      length(recipe$intercept) != 1L || is.na(recipe$intercept))
    return(NULL)
  if (!inherits(trt_recipe, "cdid_dr_treatment_recipe") || !is.numeric(trt_recipe$p) ||
      length(trt_recipe$p) != 1L || !is.finite(trt_recipe$p) || trt_recipe$p < 1L ||
      trt_recipe$p != floor(trt_recipe$p))
    return(NULL)
  p_trt <- as.integer(trt_recipe$p)
  p_or <- shared$p
  if (!is.numeric(p_or) || length(p_or) != 1L || !is.finite(p_or) ||
      p_or < 0L || p_or != floor(p_or))
    return(NULL)
  p_or <- as.integer(p_or)
  kept_or_cols <- if (p_or > 0L) or_frame_index[recipe$slope_index[shared$keep]] else integer(0)
  kept_columns <- if (p_or > 0L) shared$basis$kept_columns else NULL
  beta <- if (p_or > 0L) shared$beta else numeric(0)
  if (!is.double(beta) || length(beta) != p_or || any(!is.finite(beta))) return(NULL)
  if (p_or > 0L && (!is.character(kept_columns) || length(kept_columns) != p_or)) return(NULL)

  for (t in seq_along(specs)) {
    spec <- specs[[t]]
    if (!is.list(spec) || !is.null(spec$reference)) return(NULL)
    prep_t <- dr_baseline_shared_prepared(share, spec, supports, reader, trt_recipe,
      kept_or_cols, p_or, glm_control, st$N)
    if (is.null(prep_t) || !identical(prep_t$kept_columns, kept_columns)) return(NULL)
  }

  prep <- share$prep; inp <- share$inputs
  # The native summaries use unmasked odds derivatives: a capped control odds
  # declines to the per-model path, which applies the active-odds mask.
  lv_s <- share$lv
  if (!is.list(lv_s) || !is.double(lv_s$active) || length(lv_s$active) != length(inp$y) ||
      any(lv_s$active[inp$y == 0] != 1))
    return(NULL)
  n <- as.integer(prep$n)
  m <- length(inp$cohorts)
  M_batch <- length(specs)
  rows <- inp$rows
  X_trt <- inp$X
  if (!is.integer(rows) || length(rows) != n || n < 1L || n > st$N || anyNA(rows) ||
      any(rows < 1L | rows > st$N) || is.unsorted(rows, strictly = TRUE))
    return(NULL)
  if (!is.matrix(X_trt) || !is.double(X_trt) || nrow(X_trt) != n || ncol(X_trt) != p_trt ||
      !identical(colnames(X_trt), trt_recipe$columns))
    return(NULL)
  if (!is.matrix(prep$Z_or) || !is.double(prep$Z_or) || nrow(prep$Z_or) != n || ncol(prep$Z_or) != p_or)
    return(NULL)
  if (!is.double(prep$wp) || length(prep$wp) != n || !is.double(prep$wm) || length(prep$wm) != n)
    return(NULL)
  if (!is.double(prep$sum_wp) || length(prep$sum_wp) != 1L || !is.finite(prep$sum_wp) || prep$sum_wp <= 0 ||
      !is.double(prep$sum_wm) || length(prep$sum_wm) != 1L || !is.finite(prep$sum_wm) || prep$sum_wm <= 0 ||
      !is.double(prep$dp) || length(prep$dp) != 1L || !is.double(prep$dm) || length(prep$dm) != 1L)
    return(NULL)
  if (!is.list(prep$cohort_rows) || length(prep$cohort_rows) != m ||
      sum(lengths(prep$cohort_rows)) != n)
    return(NULL)
  for (k in seq_len(m)) {
    rr_k <- prep$cohort_rows[[k]]
    if (!is.integer(rr_k) || anyNA(rr_k) || any(rr_k < 1L | rr_k > n) ||
        is.unsorted(rr_k, strictly = TRUE))
      return(NULL)
  }
  if (!is.matrix(prep$score_sums) || !is.double(prep$score_sums) ||
      nrow(prep$score_sums) != p_trt || ncol(prep$score_sums) != m)
    return(NULL)
  if (!is.double(prep$wm_design) || length(prep$wm_design) != p_trt) return(NULL)

  post <- panel$postings
  n_dates <- length(post$date_start)
  for (t in seq_len(M_batch)) {
    spec <- specs[[t]]
    pair <- spec$pair
    if (!is.list(pair)) return(NULL)
    ds <- pair$s_calendar_index; dt <- pair$t_calendar_index
    if (!is.numeric(ds) || length(ds) != 1L || !is.finite(ds) || ds != floor(ds) ||
        ds < 1L || ds > n_dates || !(as.integer(ds) %in% complete$dates))
      return(NULL)
    if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt != floor(dt) ||
        dt < 1L || dt > n_dates || !(as.integer(dt) %in% complete$dates))
      return(NULL)
    rid <- spec$regression
    if (!is.integer(rid) || length(rid) != 1L || is.na(rid) || rid < 1L || rid > length(regressions))
      return(NULL)
    reg <- regressions[[rid]]
    if (!is.list(reg) || !is.double(reg$L) || length(reg$L) != 1L || !is.finite(reg$L))
      return(NULL)
    alpha_t <- if (recipe$intercept) reg$alpha else 0
    if (!is.double(alpha_t) || length(alpha_t) != 1L || !is.finite(alpha_t))
      return(NULL)
    sids <- spec$support_ids
    rates_t <- vapply(supports[sids], function(u) {
      if (is.list(u) && is.double(u$rate) && length(u$rate) == 1L && is.finite(u$rate) && u$rate != 0)
        u$rate else NA_real_
    }, numeric(1))
    if (anyNA(rates_t)) return(NULL)
  }

  share_entries <- n * (p_trt + p_or + 14) + p_trt * p_trt + p_trt * (m + 1) + 2 * m + 128
  per_model_summary <- 12 + 3 * p_or + p_trt * m + 8 * (m + p_trt) + 16
  batch_native_entries <- 2 * n + M_batch * (per_model_summary + 4 + m) + 128
  batch_work <- 64 * M_batch * n * (p_or + p_trt + m + 8) + 128 * M_batch * per_model_summary
  if (!is.null(scratch)) {
    ipw_owner_safe_sum(c(share_entries, batch_native_entries), scratch,
      "complete-support DR treatment summary batch")
  }
  retained_res <- if (!is.null(complete$fit_retained_reserved)) complete$fit_retained_reserved else
    c(reader$accounting$owned_entries, reader$accounting$borrowed_owner_entries)
  work_res <- if (!is.null(complete$fit_work_reserved)) complete$fit_work_reserved else batch_work
  out_res <- M_batch * (p_trt * p_trt + 12 * (p_trt + p_or) + 8 * m + p_trt * m + 4096)
  cp04_budget(ctrl, retained_res, c(share_entries, batch_native_entries), work_res, out_res,
    "baseline DR treatment summary batch")
  cp04_finite(c(prep$dp, prep$dm), "treated/control normalizers", TRUE)

  pairs_mat <- matrix(0L, 2L, M_batch)
  model_par <- matrix(0, 2L + m, M_batch)
  for (t in seq_len(M_batch)) {
    spec <- specs[[t]]
    reg <- regressions[[spec$regression]]
    pairs_mat[1L, t] <- as.integer(spec$pair$s_calendar_index)
    pairs_mat[2L, t] <- as.integer(spec$pair$t_calendar_index)
    model_par[1L, t] <- reg$L
    model_par[2L, t] <- if (recipe$intercept) reg$alpha else 0
    model_par[seq.int(3L, 2L + m), t] <- vapply(supports[spec$support_ids], `[[`, numeric(1), "rate")
  }
  ld <- isTRUE(unname(capabilities("long.double")))
  norm_vec <- as.double(c(prep$sum_wp, prep$sum_wm, prep$dp, prep$dm, st$N))
  ctl_vec <- as.integer(c(ctrl$max_chunk_rows, ld))
  kernel <- or_baseline_native_kernel("C_gmm_cp08_dr_treatment_batch")
  res <- call_fn(kernel, panel$observations$y, post$date_start, post$post_a_index, post$post_obs_index,
    rows, st$a_index, X_trt, prep$Z_or, beta, prep$wp, prep$wm, prep$cohort_rows,
    prep$score_sums, prep$wm_design, norm_vec, pairs_mat, model_par, ctl_vec, ctrl$max_matrix_entries)
  if (!identical(res$status[1L], 0L)) {
    if (identical(res$status[1L], 14L)) cp04_finite(NA_real_, "treated/control DR and IPW means")
    if (identical(res$status[1L], 15L)) cp04_finite(NA_real_, "baseline DR analytic cross derivatives")
    or_baseline_native_refuse(res$status, dr_baseline_abort)
  }

  share$reuses <- share$reuses + M_batch
  share$prep_reuses <- share$prep_reuses + M_batch
  share$summary_batches <- (if (is.null(share$summary_batches)) 0L else share$summary_batches) + 1L
  share$summary_batched_models <- (if (is.null(share$summary_batched_models)) 0L else share$summary_batched_models) + M_batch

  base_fit <- share$fit
  names(base_fit$coefficients) <- trt_recipe$columns
  u_dimnames <- list(c("T", "M", "dp", "dm"), NULL)
  v_dimnames <- list(c("T", "M", "dp", "dm"), trt_recipe$columns)
  mean_names <- c("T", "M", "dp", "dm", "T_dr", "M_dr", "mu_dr", "T_ipw", "M_ipw", "mu_ipw", "mu_or", "mu_dr_decomp")
  empty_dbar <- numeric(0)
  out <- vector("list", M_batch)
  for (t in seq_len(M_batch)) {
    fit <- base_fit
    fit$support_ids <- as.integer(specs[[t]]$support_ids)
    fit$cross <- matrix(res$cross[, t], nrow = p_trt, ncol = m)
    fit$partials <- list(
      upstream = matrix(res$U_dr[, t], nrow = 4L, ncol = m, dimnames = u_dimnames),
      beta = matrix(res$V_dr[, t], nrow = 4L, ncol = p_trt, dimnames = v_dimnames)
    )
    fit$partials_ipw <- list(
      upstream = matrix(res$U_ipw[, t], nrow = 4L, ncol = m, dimnames = u_dimnames),
      beta = matrix(res$V_ipw[, t], nrow = 4L, ncol = p_trt, dimnames = v_dimnames)
    )
    mn <- res$means[, t]
    names(mn) <- mean_names
    fit$means <- mn
    if (p_or > 0L) {
      dbp <- as.vector(res$d_bar_plus[, t])
      dbm <- as.vector(res$d_bar_minus[, t])
      dbd <- as.vector(res$d_bar_diff[, t])
      names(dbp) <- names(dbm) <- names(dbd) <- kept_columns
      fit$d_bar_plus <- dbp
      fit$d_bar_minus <- dbm
      fit$d_bar_diff <- dbd
    } else {
      fit$d_bar_plus <- empty_dbar
      fit$d_bar_minus <- empty_dbar
      fit$d_bar_diff <- empty_dbar
    }
    fit$clipping <- ipw_clip_summary(prob = lv_s$prob,
      rates = vapply(supports[specs[[t]]$support_ids], `[[`, numeric(1), "rate"),
      controls = inp$y == 0, weights = inp$w)
    out[[t]] <- fit
  }
  out
}

# Fit ONE distinct conditional_mtar treatment propensity model (pair, g, K) and
# compute both residualized DR and raw IPW moments/partials without refitting.
dr_selection_treatment_fit <- function(reader, spec, sel_fit, up_map, reg, shared,
                                       recipe, or_frame_index, trt_recipe,
                                       glm_control, ctrl, reference = NULL) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p_trt <- trt_recipe$p
  p_or <- shared$p
  kept_or_cols <- if (p_or > 0L) or_frame_index[recipe$slope_index[shared$keep]] else integer(0)
  m <- length(spec$cohorts)
  n <- spec$n_observers

  X_trt <- matrix(0, n, p_trt, dimnames = list(NULL, trt_recipe$columns))
  Z_or <- matrix(0, n, p_or, dimnames = if (p_or > 0L) list(NULL, shared$basis$kept_columns) else NULL)
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
    a_idx <- st$a_index[rows]
    if (is.null(reference)) {
      X_trt[ii, ] <- frame$Z_a[a_idx, trt_recipe$index, drop = FALSE]
      if (p_or > 0L) Z_or[ii, ] <- frame$Z_a[a_idx, kept_or_cols, drop = FALSE]
    } else {
      Xr <- chained_reference_x(reference, spec$reference, rows)
      X_trt[ii, ] <- Xr[, trt_recipe$index, drop = FALSE]
      if (p_or > 0L) Z_or[ii, ] <- Xr[, kept_or_cols, drop = FALSE]
    }
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
  if (pos != n) dr_corollary1_abort("Observed risk support changed from its owner view.", "identity")

  fit <- cp04_logit(X_trt, y, w, glm_control, st$N, spec)
  names(fit$coefficients) <- trt_recipe$columns
  lv <- ipw_logit_values(X_trt, fit$coefficients)
  wp <- w * y
  wm <- w * (1 - y) * lv$odds
  dp <- sum(wp) / st$N
  dm <- sum(wm) / st$N
  cp04_finite(c(dp, dm), "treated/control normalizers", TRUE)

  # Raw outcome IPW quantities
  T_ipw <- sum(wp * dy) / sum(wp)
  M_ipw <- sum(wm * dy) / sum(wm)
  mu_ipw <- T_ipw - M_ipw

  # Residualized DR quantities reusing fitted OR parameters (no treatment refit)
  d_obs <- reg$L * Z_or
  m_hat <- (if (recipe$intercept) reg$alpha else 0) +
    (if (p_or > 0L) as.vector(d_obs %*% shared$beta) else 0)
  r_obs <- dy - m_hat
  T_dr <- sum(wp * r_obs) / sum(wp)
  M_dr <- sum(wm * r_obs) / sum(wm)
  mu_dr <- T_dr - M_dr
  d_bar_plus <- if (p_or > 0L) as.vector(colSums(wp * d_obs) / sum(wp)) else numeric(0)
  d_bar_minus <- if (p_or > 0L) as.vector(colSums(wm * d_obs) / sum(wm)) else numeric(0)
  d_bar_diff <- d_bar_plus - d_bar_minus
  if (p_or > 0L) {
    names(d_bar_plus) <- names(d_bar_minus) <- names(d_bar_diff) <- shared$basis$kept_columns
  }
  mu_dr_decomp <- mu_ipw - (if (p_or > 0L) sum(d_bar_diff * shared$beta) else 0)
  cp04_finite(c(T_ipw, M_ipw, mu_ipw, T_dr, M_dr, mu_dr, d_bar_plus, d_bar_minus, mu_dr_decomp),
              "treated/control DR and IPW means")

  u_ids <- sort(unique(unlist(up_map$pair_channels[spec$selection_ids])))
  u_len <- length(u_ids)
  C <- matrix(0, p_trt, u_len)
  U_dr <- matrix(0, 4L, u_len, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  U_ipw <- matrix(0, 4L, u_len, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  score <- w * (y - lv$prob) * X_trt
  residual_dr <- cbind(wp * (r_obs - T_dr) / dp, wm * (r_obs - M_dr) / dm, wp, wm)
  residual_ipw <- cbind(wp * (dy - T_ipw) / dp, wm * (dy - M_ipw) / dm, wp, wm)

  for (k in seq_len(m)) {
    rr <- which(kk == k)
    if (!length(rr)) next
    q <- spec$selection_ids[k]
    ch_k <- up_map$pair_channels[[q]]
    pos_k <- match(ch_k, u_ids)
    ld_mat <- matrix(0, length(rr), length(ch_k))
    for (cs in seq.int(1L, length(rr), by = ctrl$max_chunk_rows)) {
      ce <- min(length(rr), cs + ctrl$max_chunk_rows - 1L)
      sub <- cs:ce
      ld_mat[sub, ] <- selection_pair_cohort_log_deriv(sel_fit, up_map, q, obs_rows[rr[sub]])$log_deriv
    }
    C[, pos_k] <- C[, pos_k] + crossprod(score[rr, , drop = FALSE], ld_mat) / st$N
    U_dr[, pos_k] <- U_dr[, pos_k] + crossprod(residual_dr[rr, , drop = FALSE], ld_mat) / st$N
    U_ipw[, pos_k] <- U_ipw[, pos_k] + crossprod(residual_ipw[rr, , drop = FALSE], ld_mat) / st$N
  }

  wm_act <- wm * lv$active
  V_dr <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), trt_recipe$columns))
  V_dr[2L, ] <- colSums(wm_act * (r_obs - M_dr) * X_trt) / (st$N * dm)
  V_dr[4L, ] <- colSums(wm_act * X_trt) / st$N

  V_ipw <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), trt_recipe$columns))
  V_ipw[2L, ] <- colSums(wm_act * (dy - M_ipw) * X_trt) / (st$N * dm)
  V_ipw[4L, ] <- V_dr[4L, ]

  cp04_finite(c(C, U_dr, U_ipw, V_dr, V_ipw), "Corollary 1 DR analytic cross derivatives")

  fit$upstream_ids <- as.integer(u_ids)
  fit$cross <- C
  fit$partials <- list(upstream = U_dr, beta = V_dr)
  fit$partials_ipw <- list(upstream = U_ipw, beta = V_ipw)
  fit$means <- c(
    T = T_dr, M = M_dr, dp = dp, dm = dm,
    T_dr = T_dr, M_dr = M_dr, mu_dr = mu_dr,
    T_ipw = T_ipw, M_ipw = M_ipw, mu_ipw = mu_ipw,
    mu_or = T_dr, mu_dr_decomp = mu_dr_decomp
  )
  fit$d_bar_plus <- d_bar_plus
  fit$d_bar_minus <- d_bar_minus
  fit$d_bar_diff <- d_bar_diff
  fit$clipping <- ipw_clip_summary(prob = lv$prob, q = pp, controls = y == 0, weights = w)
  fit
}

dr_populate_link_records <- function(links, state_links, treatments, pools) {
  for (i in seq_along(links)) {
    l_or <- state_links[[i]]
    tid <- links[[i]]$model
    trt <- treatments[[tid]]
    links[[i]] <- c(
      l_or,
      list(
        model = as.integer(tid),
        control_support_ids = as.integer(pools[[l_or$regression]]$support_ids),
        mu = unname(trt$means["mu_dr"]),
        mu_dr = unname(trt$means["mu_dr"]),
        mu_or = unname(l_or$mu),
        mu_ipw = unname(trt$means["mu_ipw"]),
        mu_dr_decomp = unname(trt$means["mu_dr_decomp"]),
        T_dr = unname(trt$means["T_dr"]),
        M_dr = unname(trt$means["M_dr"]),
        T_ipw = unname(trt$means["T_ipw"]),
        M_ipw = unname(trt$means["M_ipw"]),
        d_plus = unname(trt$means["dp"]),
        d_minus = unname(trt$means["dm"]),
        d_bar_plus = trt$d_bar_plus,
        d_bar_minus = trt$d_bar_minus,
        d_bar_diff = trt$d_bar_diff,
        v_ps = trt$partials$beta[2L, ],
        v_ipw = trt$partials_ipw$beta[2L, ]
      )
    )
    # Ensure l$mu is the DR point estimate while preserving l$mu_or for OR
    links[[i]]$mu <- unname(trt$means["mu_dr"])
  }
  links
}

validate_dr_baseline_backend <- function(backend) {
  if (!inherits(backend, "cdid_dr_baseline") ||
      !ipw_owner_seal_ok(backend, "dr_baseline_backend"))
    dr_baseline_abort("Expected the immutable private baseline DR backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) &&
                exists("cp04_reader", backend$reader_ref, inherits = FALSE))
    get("cp04_reader", backend$reader_ref, inherits = FALSE) else NULL
  validate_cp04_reader(reader)
  if (!is.null(backend$reference_ref)) chained_backend_reference(backend)
  if (!is.null(backend$complete_support))
    or_baseline_complete_panel(reader, backend$complete_support, dr_baseline_abort)
  invisible(backend)
}

validate_dr_corollary1_backend <- function(backend) {
  if (!inherits(backend, "cdid_dr_corollary1") ||
      !ipw_owner_seal_ok(backend, "dr_corollary1_backend"))
    dr_corollary1_abort("Expected the immutable private Corollary 1 DR backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) &&
                exists("cp04_reader", backend$reader_ref, inherits = FALSE))
    get("cp04_reader", backend$reader_ref, inherits = FALSE) else NULL
  validate_cp04_reader(reader)
  if (!is.null(backend$reference_ref)) chained_backend_reference(backend)
  invisible(backend)
}

validate_dr_backend <- function(backend) {
  if (inherits(backend, "cdid_dr_baseline")) return(validate_dr_baseline_backend(backend))
  if (inherits(backend, "cdid_dr_corollary1")) return(validate_dr_corollary1_backend(backend))
  dr_baseline_abort("Expected an immutable private DR backend (`cdid_dr_baseline` or `cdid_dr_corollary1`).", "integrity")
}

# ------------------------------------------------------------------------------
# Route 1: Baseline MTAR Chained DR (`selection = "none"`)
# ------------------------------------------------------------------------------

# `reference` (D60) is the sealed recipe-indexed reader of exactly these views
# (post catalog rows, then auxiliary descriptors): supports (pair, cohort,
# recipe), blocks (pair, pool, recipe) and treatment logits (pair, g, K, recipe)
# each occur once; ONE global shared slope spans every block.
fit_dr_baseline <- function(reader, link_views, control = ipw_owner_control(),
                            glm_control, rank_tol = or_public_rank_tol(),
                            outcome_recipe, treatment_columns = NULL,
                            retained_external_entries = 0, external_work_ops = 0,
                            automatic_chunk_rows = FALSE, reference = NULL) {
  validate_cp04_reader(reader)
  if (missing(rank_tol)) dr_baseline_abort("Supply explicit private rank_tol.")
  rank_tol <- or_baseline_rank_tol(rank_tol)
  ctrl <- validate_ipw_owner_control(control)
  gc <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    dr_baseline_abort("Invalid external pipeline reservation.")
  if (!identical(automatic_chunk_rows, isTRUE(automatic_chunk_rows)))
    dr_baseline_abort("Invalid private automatic row-planning flag.")

  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (missing(outcome_recipe))
    dr_baseline_abort("Supply the explicit private outcome recipe (intercept flag and Z_a roles).")
  if (!is.null(reference)) {
    validate_chained_reference_reader(reference)
    if (!is.list(link_views) || length(reference$view_recipe) != length(link_views))
      dr_baseline_abort("The reference reader must index exactly the supplied link views.", "integrity")
  }
  x_cols <- if (is.null(reference)) colnames(frame$Z_a) else reference$x_columns
  or_chk <- dr_check_outcome_recipe(outcome_recipe, x_cols, "baseline")
  recipe <- or_chk$recipe
  or_frame_index <- or_chk$frame_index
  trt_recipe <- dr_check_treatment_recipe(treatment_columns, x_cols, "baseline")
  if (!is.list(link_views) || !length(link_views))
    dr_baseline_abort("Supply actual private link views.")
  # CP08 descriptor intake: the registry's eligible-route context in the view
  # position yields the same records in the same order (or_baseline_descriptor_intake).
  desc <- or_baseline_descriptor_intake(link_views, reader, frame, recipe, treatment = trt_recipe,
    reference = reference, abort = dr_baseline_abort)

  L <- if (is.null(desc)) length(link_views) else length(desc)
  J <- length(st$cohorts)
  p_or <- recipe$p
  p_trt <- trt_recipe$p
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    retained_external_entries + if (is.null(reference)) 0 else reference$accounting$owned_entries
  header <- 4096 + 1024 * L * (J + 1)
  cp04_budget(ctrl, c(base, header), header, c(external_work_ops, 128 * L * L * (J + 1)), header,
              "baseline DR metadata")

  records <- list(); charge <- list(n = 0L, records = c(entries = 0, char_bytes = 0))
  for (ell in seq_len(L)) {
    if (is.null(desc)) {
      v <- link_views[[ell]]
      validate_ipw_link_view(v, frame_owner = frame)
      source <- v$functional_key
    } else {
      v <- desc[[ell]]
      source <- v$functional_key
    }
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:"))
      dr_baseline_abort("Baseline DR requires the supplied comparison-catalog functional key.")
    key_bytes <- nchar(source, type = "bytes")
    cp04_budget(ctrl, c(base, header), header + 16 * key_bytes,
                c(external_work_ops, 128 * L * L * (J + 1) + 64 * key_bytes), header + 16 * key_bytes,
                "baseline DR catalog key")
    if (is.null(desc)) {
      if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
        dr_baseline_abort("Legacy nuisance candidates are unsupported.")
      if (!is.null(v$nuisance_ref)) {
        nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
        if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair))
          dr_baseline_abort("Cross-pair or unverified nuisance sharing is unsupported.")
      }
      g <- v$cohorts$g; K <- v$cohorts$K; pair_count <- v$support$by_cohort_pair_count
    } else {
      g <- link_views$g; K <- link_views$K; pair_count <- link_views$support$by_cohort_pair_count
    }
    target <- match(g, st$cohorts)
    control_cohorts <- match(sort(K, method = "radix"), st$cohorts)
    if (anyNA(c(target, control_cohorts)) || !isTRUE(all(pair_count >= 1)))
      dr_baseline_abort("Positive actual-pair cohort support required.", "support")
    or_pooled_transitions(v$pair, "baseline")
    records[[ell]] <- list(
      pair = v$pair, target = target, controls = control_cohorts,
      link_id = v$link_id, source_functional_key = source,
      functional_key = paste0("cmpcat:DR-baseline-shared-slope-WLS-IPW|", substring(source, 8L))
    )
    if (!is.null(reference)) {
      records[[ell]]$reference <- reference$view_recipe[ell]
      records[[ell]]$functional_key <- paste0(records[[ell]]$functional_key, "|reference:",
        reference$recipes$key[reference$view_recipe[ell]])
    }
    charge <- or_baseline_payload_append(charge, records, dr_baseline_abort)
    cp04_budget(ctrl, c(base, header), header + charge$total["entries"],
                c(external_work_ops, 128 * L * L * (J + 1) + 64 * charge$total["char_bytes"]),
                header + charge$total["entries"], "baseline DR catalog metadata")
  }
  or_baseline_payload_verify(charge, records, dr_baseline_abort)
  ids <- vapply(records, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids))
    dr_baseline_abort("Distinct nonempty link IDs required.")

  catalog <- or_baseline_catalog(records, st$cohorts, recipe, reference)
  specs <- catalog$specs
  pools <- catalog$pools
  links <- catalog$links
  # Every descriptor row's observer count is the context's complete-date fact
  # (all risk units observed), the value its heavy view's support would hold.
  tcat <- dr_build_treatment_catalog(links, pools, specs, link_views,
    n_observers = if (is.null(desc)) NULL else rep(link_views$support$n_observed_risk_units, L))
  ts <- tcat$treatment_specs
  links <- tcat$links

  S <- length(specs)
  G <- length(pools)
  T_count <- length(ts)
  max_m <- max(vapply(ts, function(x) length(x$cohorts), integer(1)))
  max_trt_n <- max(vapply(ts, `[[`, integer(1), "n_observers"))
  trt_scan_n <- sum(vapply(ts, `[[`, integer(1), "n_observers"))

  metadata <- unname(cp04_payload(list(specs, pools, links, ts, trt_recipe))["entries"])
  # CP08 stage 3: the route is chosen here from actual facts, before any scan:
  # the OR facts plus each model's exact complete risk-set observer count.
  # Treatment roots are shared only within identical ordered (g, K) groups.
  complete <- or_baseline_complete_facts(reader, recipe, records, specs, ctrl, reference, trt_recipe,
    retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
    metadata_entries = header + metadata)
  if (!is.null(complete) && !all(vapply(ts, function(s) is.integer(s$cohorts) && length(s$cohorts) >= 2L &&
      s$cohorts[1L] %in% complete$target && all(s$cohorts[-1L] %in% complete$control) &&
      identical(as.integer(s$n_observers), as.integer(sum(st$cohort_end[s$cohorts] - st$cohort_start[s$cohorts] + 1L))),
      logical(1)))) complete <- NULL
  trt_groups <- if (is.null(complete)) NULL else {
    trt_keys <- vapply(ts, function(s) paste(s$cohorts, collapse = ":"), character(1))
    lapply(unique(trt_keys), function(k) which(trt_keys == k))
  }
  name_bytes <- sum(nchar(c(recipe$columns, trt_recipe$columns), type = "bytes"))
  or_charges <- if (is.null(complete)) or_baseline_fit_charges(
    S, G, L, p_or, st$N, st$n_A, min(st$N, ctrl$max_chunk_rows),
    metadata, header, name_bytes
  ) else or_baseline_complete_fit_charges(st, specs, S, G, L, p_or, min(st$N, ctrl$max_chunk_rows),
    metadata, header, name_bytes, complete)
  trt_models_bound <- T_count * (p_trt * p_trt + 12 * (p_trt + p_or) + 8 * max_m + p_trt * max_m + 4096) +
    L * (8 * p_or + 8 * p_trt + 1024) + 8192
  trt_scratch <- max_trt_n * (p_trt + p_or + 16) + 32 * (p_trt + max_m + 4)^2 + 8192
  trt_work <- 128 * (gc$maxit + 2) * trt_scan_n * (p_trt + 2)^2 +
    128 * T_count * st$N * max_m + 128 * trt_models_bound
  if (!is.null(complete)) {
    # One logit on the first model's observers plus its prepared invariants
    # (share$inputs and share$prep) per (g, K) cohort group, followed by one
    # bounded native summary batch over the remaining models of that group.
    m_rem <- max(0L, T_count - length(trt_groups))
    share_entries <- max_trt_n * (p_trt + p_or + 14) + p_trt * p_trt + p_trt * (max_m + 1) + 2 * max_m + 128
    per_model_summary <- 12 + 3 * p_or + p_trt * max_m + 8 * (max_m + p_trt) + 16
    batch_native_entries <- 2 * max_trt_n + m_rem * (per_model_summary + 4 + max_m) + 128
    trt_work <- sum(vapply(trt_groups, function(g_idx) {
      ng <- ts[[g_idx[1L]]]$n_observers
      rem_g <- max(0L, length(g_idx) - 1L)
      128 * (gc$maxit + 2) * ng * (p_trt + 2)^2 +
        128 * ng * ((p_trt + 2)^2 + p_or + max_m) +
        64 * rem_g * ng * (p_or + p_trt + max_m + 8)
    }, numeric(1))) + 128 * trt_models_bound
    trt_scratch <- max(trt_scratch + share_entries, share_entries + batch_native_entries)
  }

  bound <- metadata + or_charges$models + trt_models_bound
  scratch <- max(or_charges$scratch, trt_scratch)
  # CP08/D78 work ownership (see fit_or_baseline): `own_work` is this fit's
  # INCREMENTAL charge (OR charges plus every treatment logit term); the
  # returned `accounting$work_ops_bound` stays CUMULATIVE (caller prefix plus
  # own work) for every downstream consumer.
  own_work <- ipw_owner_safe_sum(c(or_charges$work, trt_work), label = "baseline DR own fit work")
  work <- ipw_owner_safe_sum(
    c(own_work, external_work_ops),
    ctrl$max_work_ops,
    "baseline DR pipeline fit work"
  )
  cp04_budget(ctrl, c(base, bound), scratch, work, bound, "baseline DR fits")

  state <- if (is.null(complete)) {
    supports <- lapply(specs, function(s) dr_baseline_support(reader, s, recipe, or_frame_index,
      ctrl, reference, complete))
    ans <- or_baseline_assemble(
      recipe, supports, pools, links, rank_tol, st$N,
      function(sid, visit) dr_support_visit(reader, specs[[sid]], ctrl, or_frame_index,
                                            function(rows, S_vec, obs, dy, Z) visit(Z), reference, complete)
    )
    rm(supports)
    ans
  } else {
    complete$fit_work_reserved <- work
    complete$fit_retained_reserved <- c(base, bound)
    or_baseline_complete_state(reader, specs, pools, links, recipe, ctrl, complete, scratch, rank_tol,
      or_frame_index, dr_baseline_abort)
  }

  share <- if (is.null(complete)) NULL else new.env(parent = emptyenv())
  treatments <- if (is.null(complete)) {
    lapply(ts, function(s) {
      dr_baseline_treatment_fit(
        reader, s, state$supports, state$regressions[[s$regression]],
        state$shared, recipe, or_frame_index, trt_recipe, gc, ctrl, reference, complete, share
      )
    })
  } else if (length(trt_groups) <= 1L) {
    complete$fit_work_reserved <- ipw_owner_safe_sum(
      c(work, state$original_factor_validation_work),
      ctrl$max_work_ops,
      "fit plus original-factor validation"
    )
    first_trt <- dr_baseline_treatment_fit(
      reader, ts[[1L]], state$supports, state$regressions[[ts[[1L]]$regression]],
      state$shared, recipe, or_frame_index, trt_recipe, gc, ctrl, reference, complete, share
    )
    if (T_count <= 1L) {
      share$summary_batches <- 0L
      share$summary_batched_models <- 0L
      list(first_trt)
    } else {
      rest_trt <- dr_baseline_treatment_batch(
        reader, ts[-1L], state$supports, state$regressions, state$shared,
        recipe, or_frame_index, trt_recipe, gc, ctrl, reference, complete,
        share, scratch = scratch, call_fn = .Call
      )
      if (is.null(rest_trt)) {
        # Structural decline precedes any native summary execution. Reserve the
        # original per-model route plus the attempted batch validation before
        # its remaining models run; a reduced batch charge cannot fund it.
        n1 <- ts[[1L]]$n_observers
        validation_work <- 64 * m_rem * n1 * (p_or + p_trt + max_m + 8) +
          128 * m_rem * per_model_summary
        original_trt_work <- 128 * (gc$maxit + 2) * n1 * (p_trt + 2)^2 +
          128 * trt_scan_n * ((p_trt + 2)^2 + p_or + max_m) +
          128 * T_count * st$N * max_m + 128 * trt_models_bound
        own_work <- ipw_owner_safe_sum(c(or_charges$work, original_trt_work, validation_work),
          label = "baseline DR structural-decline original-loop own work")
        work <- ipw_owner_safe_sum(c(own_work, external_work_ops), ctrl$max_work_ops,
          "baseline DR structural-decline original-loop work")
        complete$fit_work_reserved <- ipw_owner_safe_sum(
          c(work, state$original_factor_validation_work), ctrl$max_work_ops,
          "original-loop fit plus original-factor validation")
        cp04_budget(ctrl, c(base, bound), scratch, complete$fit_work_reserved,
          bound, "baseline DR structural-decline original-loop fits")
        share$summary_batches <- 0L
        share$summary_batched_models <- 0L
        rest_trt <- lapply(ts[-1L], function(s) {
          dr_baseline_treatment_fit(
            reader, s, state$supports, state$regressions[[s$regression]],
            state$shared, recipe, or_frame_index, trt_recipe, gc, ctrl, reference, complete, share
          )
        })
      }
      c(list(first_trt), rest_trt)
    }
  } else {
    complete$fit_work_reserved <- ipw_owner_safe_sum(
      c(work, state$original_factor_validation_work),
      ctrl$max_work_ops,
      "fit plus original-factor validation"
    )
    share$calls <- share$reuses <- share$summary_batches <- share$summary_batched_models <- 0L
    out_trt <- vector("list", T_count)
    for (g_idx in trt_groups) {
      share_g <- new.env(parent = emptyenv())
      first_g <- dr_baseline_treatment_fit(
        reader, ts[[g_idx[1L]]], state$supports, state$regressions[[ts[[g_idx[1L]]]$regression]],
        state$shared, recipe, or_frame_index, trt_recipe, gc, ctrl, reference, complete, share_g
      )
      if (length(g_idx) <= 1L) {
        share_g$summary_batches <- 0L
        share_g$summary_batched_models <- 0L
        out_trt[[g_idx[1L]]] <- first_g
      } else {
        rest_g <- dr_baseline_treatment_batch(
          reader, ts[g_idx[-1L]], state$supports, state$regressions, state$shared,
          recipe, or_frame_index, trt_recipe, gc, ctrl, reference, complete,
          share_g, scratch = scratch, call_fn = .Call
        )
        if (is.null(rest_g)) {
          ng <- ts[[g_idx[1L]]]$n_observers
          rem_g <- length(g_idx) - 1L
          validation_work <- 64 * rem_g * ng * (p_or + p_trt + max_m + 8) +
            128 * rem_g * per_model_summary
          own_work <- ipw_owner_safe_sum(c(own_work,
            128 * rem_g * ng * ((p_trt + 2)^2 + p_or + max_m) + 128 * rem_g * st$N * max_m,
            validation_work), label = "baseline DR structural-decline original-loop own work")
          work <- ipw_owner_safe_sum(c(own_work, external_work_ops), ctrl$max_work_ops,
            "baseline DR structural-decline original-loop work")
          complete$fit_work_reserved <- ipw_owner_safe_sum(
            c(work, state$original_factor_validation_work), ctrl$max_work_ops,
            "original-loop fit plus original-factor validation")
          cp04_budget(ctrl, c(base, bound), scratch, complete$fit_work_reserved,
            bound, "baseline DR structural-decline original-loop fits")
          share_g$summary_batches <- 0L
          share_g$summary_batched_models <- 0L
          rest_g <- lapply(ts[g_idx[-1L]], function(s) {
            dr_baseline_treatment_fit(
              reader, s, state$supports, state$regressions[[s$regression]],
              state$shared, recipe, or_frame_index, trt_recipe, gc, ctrl, reference, complete, share_g
            )
          })
        }
        out_trt[g_idx] <- c(list(first_g), rest_g)
      }
      share$calls <- share$calls + share_g$calls
      share$reuses <- share$reuses + share_g$reuses
      share$summary_batches <- share$summary_batches +
        (if (is.null(share_g$summary_batches)) 0L else share_g$summary_batches)
      share$summary_batched_models <- share$summary_batched_models +
        (if (is.null(share_g$summary_batched_models)) 0L else share_g$summary_batched_models)
      # Only the fitted model records survive a group. Its prepared n-row
      # arrays must die before the next group's first fit allocates.
      rm(share_g, first_g)
      if (length(g_idx) > 1L) rm(rest_g)
    }
    out_trt
  }
  logit_counts <- if (is.null(share)) NULL else c(
    fits = share$calls,
    reuses = share$reuses,
    batches = if (is.null(share$summary_batches)) 0L else share$summary_batches,
    batched_models = if (is.null(share$summary_batched_models)) 0L else share$summary_batched_models
  )
  rm(share)
  links_out <- dr_populate_link_records(links, state$links, treatments, pools)

  obj <- structure(
    list(
      uses_selection_models = FALSE,
      selection_model = "none",
      selection_pooling = "none",
      reader_ref = reader$ref,
      control = ctrl,
      glm_control = gc,
      rank_tol = rank_tol,
      outcome_recipe = recipe,
      or_frame_index = or_frame_index,
      treatment_recipe = trt_recipe,
      support_specs = specs,
      supports = state$supports,
      regression_specs = pools,
      regressions = state$regressions,
      shared = state$shared,
      treatment_specs = ts,
      treatments = treatments,
      links = links_out,
      universe = "A",
      law = paste0(
        "DR:baseline:shared-one-period-slope:",
        if (recipe$intercept) "block-intercepts" else "no-intercept",
        if (is.null(reference)) ":actual-pair:w/r:common-frame-Z_a" else
          ":actual-pair:W_g/r:recipe-indexed-reference-design"
      ),
      accounting = list(
        owned_entries = 0,
        borrowed_entries = base,
        retained_external_entries = retained_external_entries,
        external_work_ops = external_work_ops,
        retained_bound = base + bound,
        fit_scratch_bound = scratch,
        work_ops_bound = work,
        own_work_ops = own_work,
        work_ops_convention = "cumulative",
        automatic_chunk_rows = automatic_chunk_rows,
        n_cohort_summaries = S,
        n_supports = S,
        n_rate_fits = S,
        n_selection_fits = 0L,
        n_selection_cells = 0L,
        n_selection_blocks = 0L,
        n_selection_parameters = 0L,
        n_regression_fits = 1L,
        n_regression_blocks = G,
        n_local_intercepts = if (recipe$intercept) G else 0L,
        n_shared_slopes = state$shared$p,
        n_slope_columns = p_or,
        n_aliased_slopes = length(state$shared$basis$aliased_columns),
        n_support_scans = S * (1L + state$shared$verification_scan),
        n_treatment_fits = T_count,
        n_treatment_columns = p_trt,
        full_F_masks_per_model = 0L,
        full_F_by_link_IF_retained = FALSE,
        inverse_matrices_retained = FALSE,
        dense_dummy_regression = FALSE,
        peak_rss_measured = FALSE
      )
    ),
    class = c("cdid_dr_baseline", "cdid_dr_backend")
  )
  if (!is.null(reference)) {
    obj$reference_ref <- reference$ref
    obj$accounting$n_reference_recipes <- nrow(reference$recipes)
  }
  if (!is.null(complete)) {
    # n_treatment_fits keeps counting the T_count retained treatment models.
    complete$fit_work_reserved <- complete$fit_retained_reserved <- NULL
    obj$complete_support <- complete
    obj$accounting$work_ops_bound <- ipw_owner_safe_sum(c(obj$accounting$work_ops_bound,
      state$original_factor_validation_work), ctrl$max_work_ops, "fit plus original-factor validation")
    obj$accounting$own_work_ops <- ipw_owner_safe_sum(c(obj$accounting$own_work_ops,
      state$original_factor_validation_work), ctrl$max_work_ops, "own fit plus original-factor validation")
    obj$accounting$complete_support <- c(or_baseline_complete_record(or_charges, complete,
        state$shared$verification_scan, state$shared$basis$rows_checked, state$common_qr_compressions),
      list(n_treatment_models = T_count, n_treatment_logit_fits = unname(logit_counts["fits"]),
        n_treatment_logit_reuses = unname(logit_counts["reuses"]),
        n_treatment_summary_batches = unname(logit_counts["batches"]),
        n_treatment_summary_batched_models = unname(logit_counts["batched_models"])))
    obj$accounting$complete_support$original_factor_validations <- state$original_factor_validations
    obj$accounting$complete_support$original_factor_qr_compressions <- state$original_factor_qr_compressions
    obj$accounting$complete_support$original_factor_validation_work <- state$original_factor_validation_work
  }
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  if (actual > bound) dr_baseline_abort("Retained baseline DR bound underestimated.", "bound")
  cp04_budget(ctrl, c(base, actual), scratch, work, actual, "retained baseline DR")
  obj$accounting$owned_entries <- actual
  ipw_owner_seal(obj, "dr_baseline_backend")
}

dr_baseline_if_envelope <- function(backend, B, links, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p_or <- backend$shared$p
  p_trt <- backend$treatment_recipe$p
  G <- length(backend$regression_specs)
  J <- length(st$cohorts)
  sets <- or_pooled_if_sets(backend, links, details)
  gids <- sets$gids
  beta_gids <- sets$beta_gids
  tids <- sort(unique(vapply(backend$links[links], `[[`, integer(1), "model")))
  trt_sids <- sort(unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids"))))
  sids <- sort(unique(c(sets$sids, trt_sids)))
  nl <- length(links)
  ns <- length(sids)
  nt <- length(tids)
  if (!is.null(backend$complete_support)) {
    sets$sids <- sids
    return(dr_baseline_complete_if_envelope(backend, B, links, details, sets, tids, J, frame))
  }

  # General (missing-panel) route: the same simultaneous-buffer convention as
  # the complete route (twice the live per-row buffers, plus the small-matrix
  # and metadata terms), from the chunk's ACTUAL incidences instead of dense
  # model-by-support products. Nothing of the complete route's eligibility is
  # assumed: every support and every model runs its own endpoint search and
  # logit, the full OR core and its hits are always charged (D90: a
  # non-detailed non-reference chunk evaluates only xi_beta, or nothing at
  # p_or = 0, within these unchanged bounds) and D60 reference reads are
  # charged where they occur.
  # Per-row inventory of evaluate_dr_baseline_if_chunk (one support, then one
  # treatment model, is active at a time; earlier locals stay bound until
  # overwritten or function exit, so the stages are SUMMED, never maximized):
  #  support stage: rows/id/ci/w/A (5), the dense rate matrix up (ns), one
  #    support's read temporaries (ix, rows[ix], S, dy, ai, keep, the native
  #    observed rows and two posting positions with the kernel's three
  #    borrowed ints per row, both outcome endpoints and their difference, the
  #    hit mask, the gathered a-indices: 16); with a reference, on_a/ia/wa,
  #    the weight and design key/hit lookups and the FULL p_ref-column design
  #    read before its kept slope columns (p_ref + 9).
  #  OR core: hits (obs, dy, z and the D60 weight): each support holds at most
  #    B rows and a row lies in exactly one cohort, so it is an observed row of
  #    at most h_max supports (the largest number of requested supports that
  #    share one cohort); the resid cache (ix, q, r, d) once per block-support
  #    incidence I_beta + I_g and once per link; h_beta, xi_beta, hb, rhs,
  #    their transposes and solve temporaries (8 p_or) with the
  #    up[, control supports] product temporary (n_ctrl); h_alpha, alpha, the
  #    up[, block supports] temporary (c_max) and the alpha-row temporaries;
  #    phi, direct and the link-row temporaries (4 (p_or + 4) in all); with
  #    details, block scores, den, num, design, dm with its products and
  #    tw/tr/energy/objective.
  #  treatment models: retained records (five B-vectors; details
  #    3 p_trt + m_t + 16) plus ONE model's transients: upstream (m_max), the
  #    pair read, ci_in/obs_local/rr, beta_score, direct, X_trt, Z_or, d_rr,
  #    w/dy/m_hat/r, prob/eta/odds and their finite check, D_ind, fw, wp, wm,
  #    the score product, xi_gamma with its cross product, transposes, solve
  #    temporaries and name copy, the five psi vectors, coord with its two
  #    products and finite check, and the previous model's still-bound
  #    xi_gamma and coord (13 p_trt + 2 p_or + 52); details add direct_ipw,
  #    coord_ipw with its products, psi_ipw and normalizers_A (32); a
  #    reference adds its full design read and lookups (p_ref + 6).
  #  link assembly: phi, its finite check and the core scaling copy (3 nl);
  #    details: five psi matrices, tau, eta, lambda and loop temporaries
  #    (5 nl + 2 J + 6).
  # Output keeps the generic per-row width with each model at its own support
  # count: 3 p_trt + m_t + 16 covers beta, gamma, beta_score, upstream,
  # coordinates, direct, six psi vectors and normalizers_A. Work keeps 128
  # operations per row unit: endpoint searches and logits per support and per
  # model (log_a), block-support residuals at their incidence, each model's
  # own upstream columns and the reference reads; the per-chunk lookups add
  # the resid-cache key scans (a list grown per incidence) and the support
  # matches. Everything is affine in B, so cp04_plan_rows plans unchanged.
  nb <- length(beta_gids); ng <- length(gids)
  m_t <- vapply(backend$treatment_specs[tids], function(s) length(s$support_ids), numeric(1))
  m_max <- max(m_t)
  inc <- function(g) vapply(backend$regression_specs[g], function(s) length(s$support_ids), numeric(1))
  I_beta <- sum(inc(beta_gids)); I_g <- sum(inc(gids)); c_max <- max(c(0, inc(gids)))
  n_ctrl <- length(backend$shared$control_support_ids)
  coh <- vapply(backend$support_specs[sids], function(s)
    if (length(s$cohort) == 1L && !is.na(s$cohort)) as.integer(s$cohort) else NA_integer_, integer(1))
  h_max <- if (anyNA(coh) || any(coh < 1L)) ns else max(tabulate(coh))
  ref <- chained_backend_reference(backend)
  has_ref <- !is.null(ref)
  p_ref <- if (has_ref) length(ref$x_columns) else 0L
  # D60 hits also carry their recipe weights (one slot per support row).
  width <- 4 * ns + 4 * (p_or + 1) * (nb + 1) + 4 * ng +
    sum(3 * p_trt + m_t + 16) + 8 * nl * (p_or + 4) +
    (if (details) 3 * J + 8 * nl + 16 else 8) + (if (has_ref) 4 * ns else 0)
  metadata <- 16384 + 1024 * (ns + G + nt + nl + J) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  support_live <- 5 + ns + 16 + (if (has_ref) p_ref + 9 else 0)
  core_live <- (2 + p_or + (if (has_ref) 1 else 0)) * h_max + (3 + p_or) * (I_beta + I_g + nl) +
    8 * p_or + n_ctrl + 2 * ng + c_max + 2 * nl + 4 * (p_or + 4) +
    (if (details) p_or * nb + (2 + p_or) * nl + 5 * p_or + 8 else 0)
  trt_live <- (if (details) sum(3 * p_trt + m_t + 16) else 5 * nt) +
    m_max + 13 * p_trt + 2 * p_or + 52 + (if (details) 32 else 0) + (if (has_ref) p_ref + 6 else 0)
  tail_live <- 3 * nl + (if (details) 5 * nl + 2 * J + 6 else 0)
  live <- support_live + core_live + trt_live + tail_live
  output <- B * (width + 16) + metadata
  scratch <- 2 * B * live + 128 * ((p_or + 2)^2 + (p_trt + 2)^2) + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  units <- (ns + 1) * (p_or + 2)^2 * log_a +
    (p_or + 1) * (p_or * (nb + ng) + I_beta + I_g) +
    nt * ((p_trt + 2)^2 * log_a + p_or + 8) + sum((p_trt + 10) * (m_t + 4)) +
    J + nl * (p_or + p_trt + 4) + (if (has_ref) (ns + nt) * (p_ref + 2) else 0)
  I_all <- I_beta + I_g + nl
  lookups <- (I_all + ng + nl + 1) * (ns + ng + 1) + 2 * I_all^2 + nt * (ns + m_max + 1)
  work <- 128 * B * units + 64 * lookups + 64 * metadata
  list(
    gids = gids, beta_gids = beta_gids, sids = sids, tids = tids,
    output = output, scratch = scratch, work = work,
    general = list(route = "general", live_entries_per_row = live, row_units = units,
      lookup_elements = lookups, cohort_support_max = h_max,
      block_support_incidence = c(beta = I_beta, local = I_g),
      treatment_support_incidence = sum(m_t), or_core_evaluated = details || p_or > 0L,
      reference_columns = p_ref)
  )
}

# Prospective multi-group descriptor reservation. This uses only fitted counts
# and link IDs, before support arrays, expanded slope descriptors or group
# matrices are allocated. Count repeated global descriptors even where R can
# share their storage. The single-group CP08 reservation is unchanged.
dr_baseline_group_plan_bound <- function(backend, links, sets, tids) {
  specs <- backend$treatment_specs[tids]
  keys <- vapply(specs, function(s) paste(s$cohorts, collapse = ":"), character(1))
  groups <- unique(keys)
  if (length(groups) <= 1L) return(0)
  p_or <- backend$shared$p; p_trt <- backend$treatment_recipe$p
  nb <- if (p_or > 0L) sum(vapply(backend$regression_specs[sets$beta_gids],
    function(s) length(s$support_ids), integer(1))) else 0L
  lk <- backend$links[links]
  lm <- match(vapply(lk, `[[`, integer(1), "model"), tids)
  ids <- vapply(lk, `[[`, character(1), "link_id")
  bound <- 512
  for (key in groups) {
    ti <- which(keys == key); li <- which(lm %in% ti)
    m <- length(specs[[ti[1L]]]$support_ids)
    bound <- bound + or_baseline_native_plan_bound(length(sets$sids),
      length(ti) * (5 + (4 + m) * p_trt + 6 * m + p_or + 1),
      length(li), ids[li], p_or, nb) + p_trt * p_trt + p_trt + 128
  }
  ipw_owner_safe_sum(bound, label = "multi-group native DR descriptor bound")
}

# Eligible DR IF chunks: support preparation, one treatment replay per (g, K),
# and each model's own m_t support columns. The detailed R route evaluates the
# global OR core once; native group batches repeat support/global-slope work.
# At p_or = 0, non-detailed IFs omit the zero OR term. Concurrent payload counts
# twice the live buffers, the small matrices, metadata and group descriptors.
dr_baseline_complete_if_envelope <- function(backend, B, links, details, sets, tids, J, frame) {
  p_or <- backend$shared$p; p_trt <- backend$treatment_recipe$p
  G <- length(backend$regression_specs)
  ns <- length(sets$sids); nl <- length(links); nt <- length(tids)
  m_t <- vapply(backend$treatment_specs[tids], function(s) length(s$support_ids), numeric(1))
  n_replays <- length(unique(lapply(backend$treatment_specs[tids], `[[`, "cohorts")))
  cs <- or_baseline_complete_core_charges(backend, links, details, sets)
  core <- details || p_or > 0L
  width <- 4 * ns + 4 * (p_or + 1) * (length(sets$beta_gids) + 1) + 4 * length(sets$gids) +
    sum(3 * p_trt + m_t + 16) + 8 * nl * (p_or + 4) + (if (details) 3 * J + 8 * nl + 16 else 8)
  metadata <- 16384 + 4096 * (n_replays - 1L) + 1024 * (ns + G + nt + nl + J) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  # Plus: model records, the shared replay (rows, X, eta/odds/prob, D, fw, wp,
  # wm, score), one model's transients (upstream, reads, residuals, direct,
  # xi_gamma, psi terms, coordinates, checks), link assembly (phi, copies,
  # checks; details: five psi matrices) and the details tail (tau, eta, lambda).
  live <- cs$support_live + (if (core) cs$core_live else 0) +
    (if (details) sum(3 * p_trt + m_t + 16) else 5 * nt) + 2 * p_trt + 12 +
    max(m_t) + 2 * p_or + 6 * p_trt + 47 + 3 * nl + (if (details) 20 + 5 * nl + 2 * J + 5 else 0)
  # Detailed replay retains one immutable (rr, wp, wm, score) record per
  # group. Native replay processes groups sequentially, but repeats the full
  # support/global-slope preparation and keeps the assembled B x nl output.
  native_groups <- if (details) 1L else n_replays
  if (n_replays > 1L) live <- live +
    (if (details) (n_replays - 1L) * (p_trt + 3) else nl)
  plan_bound <- if (details || n_replays <= 1L) 0 else
    dr_baseline_group_plan_bound(backend, links, sets, tids)
  units <- native_groups * (cs$support_units + (if (core) cs$core_units else 0)) + n_replays * (p_trt + 2)^2 + sum(p_trt * (m_t + 4)) +
    J + nl * (p_or + p_trt + 4)
  lookups <- native_groups * (if (core) cs$core_lookups else 0) + nt * (ns + max(m_t) + 1)
  list(gids = sets$gids, beta_gids = sets$beta_gids, sids = sets$sids, tids = tids,
    output = B * (width + 16) + metadata,
    scratch = 2 * B * live + 128 * ((p_or + 2)^2 + (p_trt + 2)^2) + metadata + 2 * plan_bound,
    work = 128 * B * units + 64 * lookups + 64 * (metadata + plan_bound),
    complete = list(route = "complete_support", live_entries_per_row = live, row_units = units,
      lookup_elements = lookups, distinct_pairs = cs$pairs, block_support_incidence = cs$incidence,
      treatment_support_incidence = sum(m_t), treatment_replays_per_chunk = n_replays, or_core_evaluated = core))
}

dr_baseline_if_budget <- function(backend, B, links, details = FALSE) {
  pre <- dr_baseline_if_envelope(backend, B, links, details)
  cp04_budget(
    backend$control,
    c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    pre$scratch,
    c(backend$accounting$work_ops_bound, pre$work),
    pre$output,
    "baseline DR IF chunk"
  )
  pre
}

# CP08 stage 3: one treatment replay per (g, K) group per eligible chunk. Every
# model in a (g, K) group shares that group's ordered rows, cohorts,
# coefficients and support factors v (checked here; the inputs were proved equal
# at fit), so lv, fw, wp, wm and the treatment score are computed once by the
# general expressions; every model keeps its own outcomes, residuals, partials,
# rate loadings and terms.
dr_baseline_complete_replay <- function(backend, replay, spec, f, rr, ci, w_rr, X_trt) {
  v <- vapply(backend$supports[spec$support_ids], `[[`, numeric(1), "v")
  cs <- backend$complete_support
  multi_grp <- is.list(cs) && (length(cs$target) > 1L || length(cs$control) > 1L)
  key <- if (multi_grp) paste(spec$cohorts, collapse = ":") else NULL
  if (!is.null(replay)) {
    if (!multi_grp || identical(spec$cohorts, replay$cohorts)) {
      if (!identical(rr, replay$rr) || !identical(spec$cohorts, replay$cohorts) ||
          !identical(f$coefficients, replay$coefficients) || !identical(v, replay$v))
        dr_baseline_abort("Complete-support treatment replay inputs differ between models.", "identity")
      return(replay)
    }
    if (!is.integer(spec$cohorts) || length(spec$cohorts) < 2L ||
        !(spec$cohorts[1L] %in% cs$target) || !all(spec$cohorts[-1L] %in% cs$control))
      dr_baseline_abort("Complete-support treatment replay inputs differ between models.", "identity")
    hit <- replay$by_group[[key]]
    if (!is.null(hit)) {
      if (!identical(rr, hit$rr) || !identical(spec$cohorts, hit$cohorts) ||
          !identical(f$coefficients, hit$coefficients) || !identical(v, hit$v))
        dr_baseline_abort("Complete-support treatment replay inputs differ between models.", "identity")
      hit$by_group <- replay$by_group
      return(hit)
    }
  }
  lv <- ipw_logit_values(X_trt, f$coefficients)
  cp04_finite(c(lv$prob, lv$eta, lv$odds), "baseline DR treatment replay")
  D_ind <- as.double(ci[rr] == spec$cohorts[1L])
  fw <- numeric(length(rr))
  for (k in seq_along(spec$cohorts)) {
    sub_k <- which(ci[rr] == spec$cohorts[k])
    if (length(sub_k)) fw[sub_k] <- w_rr[sub_k] * backend$supports[[spec$support_ids[k]]]$v
  }
  out <- list(rr = rr, cohorts = spec$cohorts, coefficients = f$coefficients, v = v,
    wp = fw * D_ind, wm = fw * (1 - D_ind) * lv$odds, beta_score = fw * (D_ind - lv$prob) * X_trt)
  if (multi_grp) {
    bg <- if (!is.null(replay) && is.list(replay$by_group)) replay$by_group else list()
    bg[[key]] <- out
    out$by_group <- bg
  }
  out
}

# CP08 stage 4, replay: compact descriptors of one requested link batch of an
# eligible complete-support DR fit (pre$sids, pre$tids and links in their
# original orders), derived once per replay from the validated sealed backend
# and sealed. Every model within a (g, K) group must present that group's
# cohorts, logit coefficients, information factor and support factors (the
# identities proved at fit), so the chunk's per-group treatment replay and solve
# serve them all.
dr_baseline_native_live <- function(B, nl, nt, nd, p_or = 0L, p_trt = 1L)
  B * ((p_trt + 1) * nt + nl + nd + 2 * p_or + 24) + ceiling(B * nd / 8) + 512 +
    if (p_or > 0L || p_trt > 1L) 4 * p_or + 2 * p_trt * p_trt else 0

dr_baseline_native_plan <- function(backend, links, B, pre = NULL) {
  validate_dr_baseline_backend(backend)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) || B < 1 ||
      B != floor(B) || B > backend$control$max_chunk_rows)
    dr_baseline_abort("Invalid native IF link batch or chunk rows.")
  if (is.null(pre)) pre <- dr_baseline_if_budget(backend, B, links, FALSE)
  trt <- backend$treatment_recipe
  p_or <- backend$shared$p; p_trt <- trt$p
  cs <- backend$complete_support
  if (is.null(cs) || !is.numeric(p_or) || length(p_or) != 1L || !is.finite(p_or) ||
      p_or < 0L || p_or != floor(p_or) || !isTRUE(backend$outcome_recipe$intercept) ||
      !is.null(backend$reference_ref) || !is.numeric(p_trt) || length(p_trt) != 1L || !is.finite(p_trt) ||
      p_trt < 1L || p_trt != floor(p_trt))
    dr_baseline_abort("The native IF requires the facts of a complete-support fit.", "integrity")
  specs <- backend$treatment_specs[pre$tids]; fits <- backend$treatments[pre$tids]
  s1 <- specs[[1L]]; f1 <- fits[[1L]]
  v_of <- function(s) vapply(backend$supports[s$support_ids], `[[`, numeric(1), "v")
  multi_grp <- length(cs$target) > 1L || length(cs$control) > 1L
  grp_keys <- vapply(specs, function(s) paste(s$cohorts, collapse = ":"), character(1))
  grp_u <- if (multi_grp) unique(grp_keys) else grp_keys[1L]
  I_beta <- if (p_or > 0L) sum(vapply(backend$regression_specs[pre$beta_gids], function(s)
    length(s$support_ids), integer(1))) else 0L
  if (length(grp_u) <= 1L) {
    v1 <- v_of(s1); m <- length(s1$support_ids)
    if (multi_grp && (!is.integer(s1$cohorts) || length(s1$cohorts) != m || m < 2L ||
        !(s1$cohorts[1L] %in% cs$target) || !all(s1$cohorts[-1L] %in% cs$control)))
      dr_baseline_abort("Complete-support treatment replay inputs differ between models.", "identity")
    for (k in seq_along(specs))
      if (!identical(specs[[k]]$cohorts, s1$cohorts) || !identical(fits[[k]]$coefficients, f1$coefficients) ||
          !identical(fits[[k]]$H_chol, f1$H_chol) || !identical(v_of(specs[[k]]), v1) ||
          !identical(dim(fits[[k]]$cross), c(as.integer(p_trt), m)) ||
          !identical(dim(fits[[k]]$partials$beta), c(4L, as.integer(p_trt))) ||
          !identical(dim(fits[[k]]$partials$upstream), c(4L, m)))
        dr_baseline_abort("Complete-support treatment replay inputs differ between models.", "identity")
    md <- vapply(specs, function(s) as.double(c(s$pair$s_calendar_index, s$pair$t_calendar_index)), numeric(2))
    sa <- or_baseline_native_support_arrays(backend, pre$sids, as.vector(md))
    mdl <- rbind(matrix(match(md, sa$dates), 2L), vapply(specs, function(s)
      as.double(c(match(s$support_ids, pre$sids), s$cohorts)), numeric(2 * m)))
    mdl_cross <- matrix(vapply(fits, function(f) as.double(t(f$cross)), numeric(p_trt * m)), p_trt * m)
    mdl_par <- vapply(seq_along(specs), function(k) {
      f <- fits[[k]]
      as.double(c(backend$regressions[[specs[[k]]$regression]]$alpha, f$means[c("T_dr", "M_dr", "dp", "dm")],
        as.double(f$partials$beta), as.double(f$partials$upstream)))
    }, numeric(5 + 4 * p_trt + 4 * m))
    lk <- backend$links[links]
    lnk <- match(vapply(lk, `[[`, integer(1), "model"), pre$tids)
    ids <- vapply(lk, `[[`, character(1), "link_id")
    dates <- or_baseline_native_int(sa$dates); sup <- or_baseline_native_int(sa$sup)
    mdl <- or_baseline_native_int(mdl); lnk <- or_baseline_native_int(lnk)
    if (is.null(dates) || is.null(sup) || is.null(mdl) || is.null(lnk) || anyNA(sa$sup_par) || anyNA(mdl_cross) ||
        anyNA(mdl_par))
      dr_baseline_abort("Native IF descriptors do not match the fitted supports, models and links.", "integrity")
    slope <- NULL
    if (p_or > 0L) {
      reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
      sa$sup_par[3L, ] <- -sa$sup_par[3L, ]
      mdl_or_par <- vapply(seq_along(specs), function(k) {
        as.double(c(backend$regressions[[specs[[k]]$regression]]$L, fits[[k]]$d_bar_diff))
      }, numeric(1 + p_or))
      if (anyNA(mdl_or_par))
        dr_baseline_abort("Native IF descriptors do not match the fitted model OR slopes.", "integrity")
      kept_cols <- backend$or_frame_index[backend$outcome_recipe$slope_index[backend$shared$keep]]
      slope <- or_baseline_native_slope_desc(backend, pre, mdl_or_par, kept_cols,
        reader$state$N, dr_baseline_abort)
    }
    bound <- or_baseline_native_plan_bound(length(pre$sids),
      length(pre$tids) * (5 + (4 + m) * p_trt + 6 * m + p_or + 1), length(links), ids,
      p_or, I_beta) + p_trt * p_trt + p_trt + 32
    ipw_owner_safe_sum(c(bound, dr_baseline_native_live(B, length(links), length(pre$tids),
      length(dates), p_or, p_trt)), pre$scratch, "complete-support native DR IF reservation")
    plan <- list(kind = "dr_baseline_native_plan", backend_ref = backend$ref, links = as.integer(links), B = B,
      sids = pre$sids, gids = pre$gids, tids = pre$tids, dates = dates, sup = sup, sup_par = sa$sup_par,
      mdl = mdl, mdl_cross = mdl_cross, mdl_par = mdl_par, lnk = lnk, dimnames = list(slope, ids),
      cohorts = s1$cohorts, replay_spec = list(support_ids = s1$support_ids, cohorts = s1$cohorts),
      replay_fit = list(coefficients = f1$coefficients), H = f1$H_chol, trt_index = trt$index, entries = 0)
    plan$entries <- unname(cp04_payload(plan)["entries"])
    if (plan$entries > bound) dr_baseline_abort("Native IF plan bound underestimated.", "bound")
    return(ipw_owner_seal(plan, "dr_baseline_native_plan"))
  }
  bound <- dr_baseline_group_plan_bound(backend, links, pre, pre$tids)
  nd_bound <- length(cs$dates)
  ipw_owner_safe_sum(c(2 * bound, dr_baseline_native_live(B, length(links), length(pre$tids),
    nd_bound, p_or, p_trt), B * length(links)), pre$scratch,
    "complete-support multi-group native DR IF reservation")
  md <- vapply(specs, function(s) as.double(c(s$pair$s_calendar_index, s$pair$t_calendar_index)), numeric(2))
  sa <- or_baseline_native_support_arrays(backend, pre$sids, as.vector(md))
  lk <- backend$links[links]
  lnk_all <- match(vapply(lk, `[[`, integer(1), "model"), pre$tids)
  ids <- vapply(lk, `[[`, character(1), "link_id")
  dates <- or_baseline_native_int(sa$dates); sup <- or_baseline_native_int(sa$sup)
  lnk_all <- or_baseline_native_int(lnk_all)
  if (is.null(dates) || is.null(sup) || is.null(lnk_all) || anyNA(sa$sup_par))
    dr_baseline_abort("Native IF descriptors do not match the fitted supports, models and links.", "integrity")
  slope <- NULL
  mdl_or_par <- NULL
  if (p_or > 0L) {
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    sa$sup_par[3L, ] <- -sa$sup_par[3L, ]
    mdl_or_par <- vapply(seq_along(specs), function(k) {
      as.double(c(backend$regressions[[specs[[k]]$regression]]$L, fits[[k]]$d_bar_diff))
    }, numeric(1 + p_or))
    if (anyNA(mdl_or_par))
      dr_baseline_abort("Native IF descriptors do not match the fitted model OR slopes.", "integrity")
    kept_cols <- backend$or_frame_index[backend$outcome_recipe$slope_index[backend$shared$keep]]
    slope <- or_baseline_native_slope_desc(backend, pre, mdl_or_par, kept_cols,
      reader$state$N, dr_baseline_abort)
  }
  groups <- vector("list", length(grp_u))
  ng_bound <- 0
  for (g in seq_along(grp_u)) {
    t_idx <- which(grp_keys == grp_u[g])
    l_idx <- which(lnk_all %in% t_idx)
    if (!length(t_idx) || !length(l_idx))
      dr_baseline_abort("Native IF descriptors do not match the fitted supports, models and links.", "integrity")
    specs_g <- specs[t_idx]; fits_g <- fits[t_idx]
    s_g <- specs_g[[1L]]; f_g <- fits_g[[1L]]
    v_g <- v_of(s_g); m_g <- length(s_g$support_ids)
    if (!is.integer(s_g$cohorts) || length(s_g$cohorts) != m_g || m_g < 2L ||
        !(s_g$cohorts[1L] %in% cs$target) || !all(s_g$cohorts[-1L] %in% cs$control))
      dr_baseline_abort("Complete-support treatment replay inputs differ between models.", "identity")
    for (k in seq_along(specs_g))
      if (!identical(specs_g[[k]]$cohorts, s_g$cohorts) || !identical(fits_g[[k]]$coefficients, f_g$coefficients) ||
          !identical(fits_g[[k]]$H_chol, f_g$H_chol) || !identical(v_of(specs_g[[k]]), v_g) ||
          !identical(dim(fits_g[[k]]$cross), c(as.integer(p_trt), m_g)) ||
          !identical(dim(fits_g[[k]]$partials$beta), c(4L, as.integer(p_trt))) ||
          !identical(dim(fits_g[[k]]$partials$upstream), c(4L, m_g)))
        dr_baseline_abort("Complete-support treatment replay inputs differ between models.", "identity")
    md_g <- md[, t_idx, drop = FALSE]
    mdl_g <- or_baseline_native_int(rbind(matrix(match(md_g, sa$dates), 2L), vapply(specs_g, function(s)
      as.double(c(match(s$support_ids, pre$sids), s$cohorts)), numeric(2 * m_g))))
    mdl_cross_g <- matrix(vapply(fits_g, function(f) as.double(t(f$cross)), numeric(p_trt * m_g)), p_trt * m_g)
    mdl_par_g <- vapply(seq_along(specs_g), function(k) {
      f <- fits_g[[k]]
      as.double(c(backend$regressions[[specs_g[[k]]$regression]]$alpha, f$means[c("T_dr", "M_dr", "dp", "dm")],
        as.double(f$partials$beta), as.double(f$partials$upstream)))
    }, numeric(5 + 4 * p_trt + 4 * m_g))
    lnk_g <- or_baseline_native_int(match(lnk_all[l_idx], t_idx))
    if (is.null(mdl_g) || is.null(lnk_g) || anyNA(mdl_cross_g) || anyNA(mdl_par_g))
      dr_baseline_abort("Native IF descriptors do not match the fitted supports, models and links.", "integrity")
    slope_g <- if (p_or > 0L) {
      sl <- slope
      sl[[10L]] <- mdl_or_par[, t_idx, drop = FALSE]
      sl
    } else NULL
    groups[[g]] <- list(t_idx = as.integer(t_idx), l_idx = as.integer(l_idx), mdl = mdl_g,
      mdl_cross = mdl_cross_g, mdl_par = mdl_par_g, lnk = lnk_g, dimnames = list(slope_g, ids[l_idx]),
      cohorts = s_g$cohorts, replay_spec = list(support_ids = s_g$support_ids, cohorts = s_g$cohorts),
      replay_fit = list(coefficients = f_g$coefficients), H = f_g$H_chol)
    ng_bound <- ng_bound + or_baseline_native_plan_bound(length(pre$sids),
      length(t_idx) * (5 + (4 + m_g) * p_trt + 6 * m_g + p_or + 1), length(l_idx), ids[l_idx],
      p_or, I_beta) + p_trt * p_trt + p_trt + 128
  }
  if (ng_bound + 512 > bound)
    dr_baseline_abort("Native IF group descriptor bound underestimated.", "bound")
  plan <- list(kind = "dr_baseline_native_plan", backend_ref = backend$ref, links = as.integer(links), B = B,
    sids = pre$sids, gids = pre$gids, tids = pre$tids, dates = dates, sup = sup, sup_par = sa$sup_par,
    ids = ids, groups = groups, trt_index = trt$index, entries = 0)
  plan$entries <- unname(cp04_payload(plan)["entries"])
  if (plan$entries > bound) dr_baseline_abort("Native IF plan bound underestimated.", "bound")
  ipw_owner_seal(plan, "dr_baseline_native_plan")
}

# One eligible non-detailed DR chunk in the original order: the rate stage's
# slots and every model's reads (kernel part 1); the chunk's per-group treatment
# replay (dr_baseline_complete_replay) after the group's reads; the
# unchanged cp04_solve once when every model's right-hand side is bitwise the
# same (else once per model); the four terms, frame-coordinate checks and link
# columns (kernel part 2). A later model's read failure is raised after the
# earlier models' frame-IF checks, as in the original model loop.
dr_baseline_native_if_chunk <- function(backend, reader, rows, links, scale, pre, plan = NULL) {
  st <- reader$state; B <- length(rows)
  p_or <- backend$shared$p; p_trt <- backend$treatment_recipe$p
  if (is.null(plan)) plan <- dr_baseline_native_plan(backend, links, B, pre) else
    or_baseline_native_check(plan, backend, links, B, pre, "dr_baseline_native_plan", dr_baseline_abort)
  ipw_owner_safe_sum(c(plan$entries, dr_baseline_native_live(B, length(links), length(plan$tids),
    length(plan$dates), p_or, p_trt)), pre$scratch, "complete-support native DR IF chunk")
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  po <- or_baseline_complete_panel(reader, backend$complete_support, dr_baseline_abort)
  post <- po$postings; rows_i <- as.integer(rows); cap <- backend$control$max_matrix_entries
  if (is.null(plan$groups)) {
    a <- .Call(or_baseline_native_kernel("C_gmm_cp08_dr_if_prepare"), po$observations$y, post$date_start,
      post$post_a_index, post$post_obs_index, rows_i, st$cohort_index, st$a_index, st$w, plan$dates,
      plan$sup, plan$sup_par, plan$mdl, plan$mdl_cross, cap)
    code <- a$status
    if (!identical(code[1L], 0L) && (identical(code[2L], 1L) || identical(code[3L], 1L)))
      or_baseline_native_refuse(code, dr_baseline_abort)
    n_ok <- if (identical(code[1L], 0L)) length(plan$tids) else code[3L] - 1L
    ci <- st$cohort_index[rows]; w <- st$w[rows]
    rr <- which(ci %in% plan$cohorts)
    beta_score <- matrix(0, B, p_trt); wp <- wm <- numeric()
    if (length(rr)) {
      X_trt <- frame$Z_a[st$a_index[rows[rr]], plan$trt_index, drop = FALSE]
      rp <- dr_baseline_complete_replay(backend, NULL, plan$replay_spec, plan$replay_fit, rr, ci, w[rr], X_trt)
      wp <- rp$wp; wm <- rp$wm; beta_score[rr, ] <- rp$beta_score
    }
    P <- a$upstream
    xi <- if (isTRUE(a$identical)) t(cp04_solve(plan$H, t(beta_score + P))) else {
      n_m <- ncol(P) %/% p_trt
      do.call(cbind, lapply(seq_len(n_m), function(k) {
        cols <- (k - 1L) * p_trt + seq_len(p_trt)
        t(cp04_solve(plan$H, t(beta_score + P[, cols, drop = FALSE])))
      }))
    }
    dn <- plan$dimnames
    if (!is.null(dn) && !is.null(dn[[1L]])) {
      dn[[1L]][[1L]][[1L]] <- frame$Z_a
    }
    b <- .Call(or_baseline_native_kernel("C_gmm_cp08_dr_if_assemble"), rows_i, st$cohort_index, st$a_index, st$w,
      a$gather, plan$sup, plan$sup_par, plan$mdl, plan$mdl_par, xi, as.integer(rr), wp, wm, plan$lnk,
      as.integer(n_ok), dn, cap)
    if (!identical(b$status[1L], 0L)) or_baseline_native_refuse(b$status, dr_baseline_abort)
    if (!identical(code[1L], 0L)) or_baseline_native_refuse(code, dr_baseline_abort)
    phi <- b$phi
  } else {
    phi <- matrix(0, B, length(links), dimnames = list(NULL, plan$ids))
    ci <- st$cohort_index[rows]; w <- st$w[rows]
    prep_kernel <- or_baseline_native_kernel("C_gmm_cp08_dr_if_prepare")
    asm_kernel <- or_baseline_native_kernel("C_gmm_cp08_dr_if_assemble")
    for (grp in plan$groups) {
      a <- .Call(prep_kernel, po$observations$y, post$date_start,
        post$post_a_index, post$post_obs_index, rows_i, st$cohort_index, st$a_index, st$w, plan$dates,
        plan$sup, plan$sup_par, grp$mdl, grp$mdl_cross, cap)
      code <- a$status
      if (!identical(code[1L], 0L) && (identical(code[2L], 1L) || identical(code[3L], 1L)))
        or_baseline_native_refuse(code, dr_baseline_abort)
      n_ok <- if (identical(code[1L], 0L)) length(grp$t_idx) else code[3L] - 1L
      rr <- which(ci %in% grp$cohorts)
      beta_score <- matrix(0, B, p_trt); wp <- wm <- numeric()
      if (length(rr)) {
        X_trt <- frame$Z_a[st$a_index[rows[rr]], plan$trt_index, drop = FALSE]
        rp <- dr_baseline_complete_replay(backend, NULL, grp$replay_spec, grp$replay_fit, rr, ci, w[rr], X_trt)
        wp <- rp$wp; wm <- rp$wm; beta_score[rr, ] <- rp$beta_score
      }
      P <- a$upstream
      xi <- if (isTRUE(a$identical)) t(cp04_solve(grp$H, t(beta_score + P))) else {
        n_m <- ncol(P) %/% p_trt
        do.call(cbind, lapply(seq_len(n_m), function(k) {
          cols <- (k - 1L) * p_trt + seq_len(p_trt)
          t(cp04_solve(grp$H, t(beta_score + P[, cols, drop = FALSE])))
        }))
      }
      dn <- grp$dimnames
      if (!is.null(dn) && !is.null(dn[[1L]])) {
        dn[[1L]][[1L]][[1L]] <- frame$Z_a
      }
      b <- .Call(asm_kernel, rows_i, st$cohort_index, st$a_index, st$w,
        a$gather, plan$sup, plan$sup_par, grp$mdl, grp$mdl_par, xi, as.integer(rr), wp, wm, grp$lnk,
        as.integer(n_ok), dn, cap)
      if (!identical(b$status[1L], 0L)) or_baseline_native_refuse(b$status, dr_baseline_abort)
      if (!identical(code[1L], 0L)) or_baseline_native_refuse(code, dr_baseline_abort)
      phi[, grp$l_idx] <- b$phi
      # Release this group's returned arrays and R replay/solve temporaries
      # before the next .Call. Only phi and the sealed plan span groups.
      rm(a, b, P, xi, beta_score, wp, wm, dn)
      if (length(rr)) rm(rp, X_trt)
    }
  }
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_rate_score_evaluations = length(pre$sids),
      n_regression_score_evaluations = if (p_or > 0L) length(union(pre$gids, pre$beta_gids)) else 0L,
      n_global_beta_blocks = if (p_or > 0L) length(pre$beta_gids) else 0L,
      n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    dr_baseline_abort("IF output bound underestimated.", "bound")
  out
}

evaluate_dr_baseline_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE, plan = NULL) {
  validate_dr_baseline_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  valid_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) &&
    all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid_index(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    dr_baseline_abort("Invalid explicit IF row/link batch, scale or details flag.")

  pre <- dr_baseline_if_budget(backend, length(rows), links, details)
  # CP08 stage 4: an eligible non-detailed chunk is one native batch; detailed
  # requests keep the route below.
  if (!details && !is.null(backend$complete_support))
    return(dr_baseline_native_if_chunk(backend, reader, rows, links, scale, pre, plan))
  B <- length(rows)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p_or <- backend$shared$p
  p_trt <- backend$treatment_recipe$p
  kept_or_cols <- if (p_or > 0L) backend$or_frame_index[backend$outcome_recipe$slope_index[backend$shared$keep]] else integer(0)
  trt_cols <- backend$treatment_recipe$index
  ref <- chained_backend_reference(backend)
  panel <- if (is.null(backend$complete_support)) NULL else
    or_baseline_complete_panel(reader, backend$complete_support, dr_baseline_abort)
  pair_rows <- function(rr, pair, read = TRUE) if (is.null(panel)) cp04_pair_rows(reader, rr, pair, read) else
    or_baseline_complete_pair_rows(reader, panel, rr, pair, read, dr_baseline_abort)
  # CP08 stage 3 / D90 / D92: the non-detailed DR IF reads only core$beta
  # (psi_or, the existing zero vector at p_or = 0). A non-detailed chunk
  # (static or D60 reference) therefore requests the beta-only core at p_or > 0
  # and at p_or = 0 skips the core and its support outcome reads (the fit read
  # and checked these sealed endpoints); every rate column is still formed from
  # its support's observed pair rows and each model reads its own pair endpoints
  # and reference recipes below. Detailed chunks keep the full core.
  core_used <- details || p_or > 0L
  beta_only <- !details

  ci <- st$cohort_index[rows]
  # A D60 chunk has no common weight: every support and treatment model reads
  # its recipe's W_g (NA elsewhere, so any stray common-weight use is refused).
  w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
  A <- as.double(st$a_index[rows] > 0L)

  up <- matrix(0, B, length(pre$sids))
  hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]
    s <- backend$support_specs[[sid]]
    u <- backend$supports[[sid]]
    ix <- which(ci == s$cohort)
    hit <- pair_rows(rows[ix], s$pair, core_used)
    obs <- ix[hit$S == 1]
    if (is.null(ref)) {
      up[ix, j] <- w[ix] * (hit$S - A[ix] * u$rate) / u$massA
      if (core_used) hits[[j]] <- list(
        obs = obs,
        dy = hit$dy[hit$S == 1],
        z = if (p_or > 0L) frame$Z_a[st$a_index[rows[obs]], kept_or_cols, drop = FALSE] else matrix(0, length(obs), 0L)
      )
    } else {
      on_a <- A[ix] > 0; ia <- ix[on_a]
      wa <- chained_reference_weight(ref, s$reference, rows[ia])
      up[ia, j] <- wa * (hit$S[on_a] - u$rate) / u$massA
      if (core_used) hits[[j]] <- list(
        obs = obs,
        dy = hit$dy[hit$S == 1],
        z = if (p_or > 0L) chained_reference_x(ref, s$reference, rows[obs])[, kept_or_cols, drop = FALSE] else
          matrix(0, length(obs), 0L),
        w = wa[match(obs, ia)]
      )
    }
  }

  # 1. Shared OR IF core (global xi_beta; full core when detailed)
  core <- if (core_used) or_pooled_if_core(backend, st$N, w, hits, up, pre, links, details, beta_only = beta_only) else NULL

  # 2. Treatment propensity roots and residualized link components
  models <- vector("list", length(backend$treatments))
  replay <- NULL # Invocation-local immutable treatment records keyed by (g, K).
  for (tj in pre$tids) {
    spec <- backend$treatment_specs[[tj]]
    f <- backend$treatments[[tj]]
    reg <- backend$regressions[[spec$regression]]
    u_pos <- match(spec$support_ids, pre$sids)
    upstream <- up[, u_pos, drop = FALSE]

    ci_in <- which(ci %in% spec$cohorts)
    hit_all <- if (length(ci_in)) pair_rows(rows[ci_in], spec$pair) else list(S = numeric(0), dy = numeric(0))
    obs_local <- which(hit_all$S == 1)
    rr <- ci_in[obs_local]

    beta_score <- matrix(0, B, p_trt)
    direct <- matrix(0, B, 4L)
    direct[, 3L] <- -f$means["dp"]
    direct[, 4L] <- -f$means["dm"]
    direct_ipw <- if (details) {
      d_ipw <- matrix(0, B, 4L)
      d_ipw[, 3L] <- -f$means["dp"]
      d_ipw[, 4L] <- -f$means["dm"]
      d_ipw
    } else NULL

    if (length(rr)) {
      a_idx <- st$a_index[rows[rr]]
      if (is.null(ref)) {
        X_trt <- frame$Z_a[a_idx, trt_cols, drop = FALSE]
        Z_or <- if (p_or > 0L) frame$Z_a[a_idx, kept_or_cols, drop = FALSE] else matrix(0, length(rr), 0L)
        w_rr <- w[rr]
      } else {
        Xr <- chained_reference_x(ref, spec$reference, rows[rr])
        X_trt <- Xr[, trt_cols, drop = FALSE]
        Z_or <- if (p_or > 0L) Xr[, kept_or_cols, drop = FALSE] else matrix(0, length(rr), 0L)
        w_rr <- chained_reference_weight(ref, spec$reference, rows[rr])
      }
      dy_rr <- hit_all$dy[obs_local]
      d_rr <- reg$L * Z_or
      m_hat_rr <- (if (backend$outcome_recipe$intercept) reg$alpha else 0) +
        (if (p_or > 0L) as.vector(d_rr %*% backend$shared$beta) else 0)
      r_rr <- dy_rr - m_hat_rr

      if (is.null(panel)) {
        lv <- ipw_logit_values(X_trt, f$coefficients)
        cp04_finite(c(lv$prob, lv$eta, lv$odds), "baseline DR treatment replay")
        D_ind <- as.double(ci[rr] == spec$cohorts[1L])
        fw <- numeric(length(rr))
        for (k in seq_along(spec$cohorts)) {
          sub_k <- which(ci[rr] == spec$cohorts[k])
          if (length(sub_k)) {
            sid <- spec$support_ids[k]
            fw[sub_k] <- w_rr[sub_k] * backend$supports[[sid]]$v
          }
        }
        wp <- fw * D_ind
        wm <- fw * (1 - D_ind) * lv$odds
        beta_score[rr, ] <- fw * (D_ind - lv$prob) * X_trt
      } else {
        replay <- dr_baseline_complete_replay(backend, replay, spec, f, rr, ci, w_rr, X_trt)
        wp <- replay$wp
        wm <- replay$wm
        beta_score[rr, ] <- replay$beta_score
      }
      direct[rr, 1L] <- wp * (r_rr - f$means["T_dr"]) / f$means["dp"]
      direct[rr, 2L] <- wm * (r_rr - f$means["M_dr"]) / f$means["dm"]
      direct[rr, 3L] <- direct[rr, 3L] + wp
      direct[rr, 4L] <- direct[rr, 4L] + wm
      if (details) {
        direct_ipw[rr, 1L] <- wp * (dy_rr - f$means["T_ipw"]) / f$means["dp"]
        direct_ipw[rr, 2L] <- wm * (dy_rr - f$means["M_ipw"]) / f$means["dm"]
        direct_ipw[rr, 3L] <- direct_ipw[rr, 3L] + wp
        direct_ipw[rr, 4L] <- direct_ipw[rr, 4L] + wm
      }
    }

    xi_gamma <- t(cp04_solve(f$H_chol, t(beta_score + upstream %*% t(f$cross))))
    colnames(xi_gamma) <- backend$treatment_recipe$columns
    psi_dir_m <- direct[, 1L] - direct[, 2L]
    psi_or_m <- if (p_or > 0L) -as.vector(core$beta %*% f$d_bar_diff) else numeric(B)
    psi_ps_m <- -as.vector(xi_gamma %*% f$partials$beta[2L, ])
    psi_sel_m <- as.vector(upstream %*% (f$partials$upstream[1L, ] - f$partials$upstream[2L, ]))
    psi_dr_m <- psi_dir_m + psi_or_m + psi_ps_m + psi_sel_m
    coord <- direct + upstream %*% t(f$partials$upstream) + xi_gamma %*% t(f$partials$beta)
    colnames(coord) <- c("T", "M", "dp", "dm")
    cp04_finite(c(xi_gamma, coord, psi_dr_m), "complete baseline DR frame IF")

    if (details) {
      coord_ipw <- direct_ipw + upstream %*% t(f$partials_ipw$upstream) + xi_gamma %*% t(f$partials_ipw$beta)
      colnames(coord_ipw) <- c("T", "M", "dp", "dm")
      models[[tj]] <- list(
        beta = xi_gamma,
        gamma = xi_gamma,
        coordinates = coord,
        direct = direct,
        beta_score = beta_score,
        upstream = upstream,
        psi_dr = psi_dr_m,
        psi_dir = psi_dir_m,
        psi_or = psi_or_m,
        psi_ps = psi_ps_m,
        psi_sel = psi_sel_m,
        psi_ipw = coord_ipw[, 1L] - coord_ipw[, 2L],
        normalizers_A = coord[, 3:4, drop = FALSE] / st$lambda -
          outer(A - st$lambda, f$means[c("dp", "dm")] / st$lambda^2)
      )
    } else {
      models[[tj]] <- list(
        phi = psi_dr_m,
        psi_dir = psi_dir_m,
        psi_or = psi_or_m,
        psi_ps = psi_ps_m,
        psi_sel = psi_sel_m
      )
    }
  }

  link_names <- vapply(backend$links[links], `[[`, character(1), "link_id")
  phi <- matrix(0, B, length(links), dimnames = list(NULL, link_names))
  if (details) {
    psi_dir_mat <- psi_or_mat <- psi_ps_mat <- psi_sel_mat <- psi_ipw_mat <-
      matrix(0, B, length(links), dimnames = list(NULL, link_names))
  }
  for (j in seq_along(links)) {
    m <- models[[backend$links[[links[j]]]$model]]
    phi[, j] <- if (details) m$psi_dr else m$phi
    if (details) {
      psi_dir_mat[, j] <- m$psi_dir
      psi_or_mat[, j] <- m$psi_or
      psi_ps_mat[, j] <- m$psi_ps
      psi_sel_mat[, j] <- m$psi_sel
      psi_ipw_mat[, j] <- m$psi_ipw
    }
  }
  cp04_finite(phi, "complete baseline DR link IF")
  if (scale == "core") phi <- st$lambda * phi

  out <- list(
    rows = rows,
    id = st$id[rows],
    links = links,
    phi = phi,
    scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(
      n_rate_score_evaluations = length(pre$sids),
      n_regression_score_evaluations = if (core_used) length(union(pre$gids, pre$beta_gids)) else 0L,
      n_global_beta_blocks = if (core_used) length(pre$beta_gids) else 0L,
      n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  if (details) {
    # Frame tau/eta summaries need a common fixed weight; a D60 context forms
    # them only where an existing consumer reads them (none here), so NULL.
    tau <- eta <- if (is.null(ref)) matrix(0, B, length(st$cohorts)) else NULL
    if (is.null(ref)) for (k in seq_along(st$cohorts)) {
      tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
      eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
    }
    out$details <- list(
      support_ids = pre$sids,
      rates = up,
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
      treatment_ids = pre$tids,
      models = models,
      psi_dir = psi_dir_mat,
      psi_or = psi_or_mat,
      psi_ps = psi_ps_mat,
      psi_sel = psi_sel_mat,
      psi_ipw = psi_ipw_mat,
      psi_or_estimator = core$phi,
      tau = tau,
      eta = eta,
      lambda = A - st$lambda
    )
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    dr_baseline_abort("IF output bound underestimated.", "bound")
  out
}

replay_dr_baseline_if <- function(backend, chunk_rows, link_batch_size, callback,
                                  scale = "frame", details = FALSE) {
  validate_dr_baseline_backend(backend)
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  LB <- ipw_owner_int_bound(link_batch_size, "link_batch_size")
  if (B > backend$control$max_chunk_rows || !is.function(callback) ||
      !identical(details, isTRUE(details)) || !is.character(scale) || length(scale) != 1L ||
      is.na(scale) || !scale %in% c("frame", "core"))
    dr_baseline_abort("Invalid explicit replay controls/callback.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N
  L <- length(backend$links)
  total <- backend$accounting$work_ops_bound
  native <- !details && !is.null(backend$complete_support)
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    pre <- dr_baseline_if_budget(backend, min(B, N), ll, details)
    total <- ipw_owner_safe_sum(
      c(total, ceiling(N / B) * pre$work),
      backend$control$max_work_ops,
      "total baseline DR replay work"
    )
    if (native) dr_baseline_native_plan(backend, ll, min(B, N), pre)
  }
  batches <- 0L
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    plan <- if (native) dr_baseline_native_plan(backend, ll, min(B, N)) else NULL
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1L))
      batch <- evaluate_dr_baseline_if_chunk(backend, rr, ll, scale, details, plan)
      callback(batch)
      rm(batch)
      batches <- batches + 1L
    }
  }
  invisible(list(
    n_frame_rows = N, n_links = L, batches = batches, scale = scale,
    normalization_count = if (scale == "core") reader$state$n_A else N,
    work_ops_bound = total, all_F_identities = TRUE
  ))
}

# ------------------------------------------------------------------------------
# Route 2 & Route 3: Conditional MTAR Chained DR (`pair_logit` & `independent_period_logit`)
# ------------------------------------------------------------------------------

# `reference` (D60): selection cells, supports, blocks and treatment logits
# are keyed by recipe; ONE global shared slope spans every block.
fit_dr_corollary1 <- function(reader, link_views, control = ipw_owner_control(),
                              glm_control, rank_tol = or_public_rank_tol(),
                              outcome_recipe, treatment_columns = NULL,
                              retained_external_entries = 0, external_work_ops = 0,
                              automatic_chunk_rows = FALSE,
                              selection_model = "pair_logit",
                              selection_pooling = "separate",
                              sel_fit = NULL, reference = NULL) {
  validate_cp04_reader(reader)
  if (missing(rank_tol) || !gmm_scalar_numeric(rank_tol) || rank_tol <= 0 || rank_tol >= 1)
    dr_corollary1_abort("Supply explicit finite private rank_tol in (0, 1).")
  ctrl <- validate_ipw_owner_control(control)
  gc <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    dr_corollary1_abort("Invalid external pipeline reservation.")
  if (!identical(automatic_chunk_rows, isTRUE(automatic_chunk_rows)))
    dr_corollary1_abort("Invalid private automatic row-planning flag.")

  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (missing(outcome_recipe))
    dr_corollary1_abort("Supply the explicit private outcome recipe (intercept flag and Z_a roles).")
  if (!is.null(reference)) {
    validate_chained_reference_reader(reference)
    if (!is.list(link_views) || length(reference$view_recipe) != length(link_views))
      dr_corollary1_abort("The reference reader must index exactly the supplied link views.", "integrity")
  }
  x_cols <- if (is.null(reference)) colnames(frame$Z_a) else reference$x_columns
  or_chk <- dr_check_outcome_recipe(outcome_recipe, x_cols, "corollary1")
  recipe <- or_chk$recipe
  or_frame_index <- or_chk$frame_index
  trt_recipe <- dr_check_treatment_recipe(treatment_columns, x_cols, "corollary1")
  if (!is.list(link_views) || !length(link_views))
    dr_corollary1_abort("Supply actual private link views.")

  L <- length(link_views)
  J <- length(st$cohorts)
  p_or <- recipe$p
  p_trt <- trt_recipe$p
  d <- if (is.null(reference)) ncol(st$Z) else length(reference$designs[[reference$role_slot["s"]]]$columns)
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    retained_external_entries + if (is.null(reference)) 0 else reference$accounting$owned_entries
  header <- 4096 + 1024 * L * (J + 1)
  cp04_budget(ctrl, c(base, header), header, c(external_work_ops, 128 * L * L * (J + 1)), header,
              "Corollary 1 DR metadata")

  records <- list()
  for (ell in seq_along(link_views)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:"))
      dr_corollary1_abort("Corollary 1 DR requires the supplied comparison-catalog functional key.")
    key_bytes <- nchar(source, type = "bytes")
    cp04_budget(ctrl, c(base, header), header + 16 * key_bytes,
                c(external_work_ops, 128 * L * L * (J + 1) + 64 * key_bytes), header + 16 * key_bytes,
                "Corollary 1 DR catalog key")
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
      dr_corollary1_abort("Legacy nuisance candidates are unsupported.")
    if (!is.null(v$nuisance_ref)) {
      nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair))
        dr_corollary1_abort("Cross-pair or unverified nuisance sharing is unsupported.")
    }
    target <- match(v$cohorts$g, st$cohorts)
    control_cohorts <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    if (anyNA(c(target, control_cohorts)) || !isTRUE(all(v$support$by_cohort_pair_count >= 1)))
      dr_corollary1_abort("Positive actual-pair cohort support required.", "support")
    or_pooled_transitions(v$pair, "corollary1")
    records[[ell]] <- list(
      pair = v$pair, target = target, controls = control_cohorts,
      link_id = v$link_id, source_functional_key = source,
      functional_key = paste0(
        "cmpcat:DR-Cor1-", st$universe, "-", selection_model, "-", selection_pooling,
        "-shared-slope-WLS-IPW-a/(rho*p)|", substring(source, 8L)
      )
    )
    if (!is.null(reference)) {
      records[[ell]]$reference <- reference$view_recipe[ell]
      records[[ell]]$functional_key <- paste0(records[[ell]]$functional_key, "|reference:",
        reference$recipes$key[reference$view_recipe[ell]])
    }
    cp04_budget(ctrl, c(base, header), header + cp04_payload(records)["entries"],
                c(external_work_ops, 128 * L * L * (J + 1) + 64 * cp04_payload(records)["char_bytes"]),
                header + cp04_payload(records)["entries"], "Corollary 1 DR catalog metadata")
  }
  ids <- vapply(records, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids))
    dr_corollary1_abort("Distinct nonempty link IDs required.")

  # Fit selection nuisance roots ONCE per request (or reuse supplied sel_fit)
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
      dr_corollary1_abort("Supplied `sel_fit` does not match reader or selection specification.", "integrity")
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
  # Build ONE canonical upstream map covering all active cohorts and selection blocks
  cohorts_a <- sort(unique(vapply(specs, `[[`, integer(1), "cohort")))
  up_map <- build_selection_upstream_map(sel_fit, cohorts_a = cohorts_a)
  n_a <- up_map$n_a
  D <- up_map$D
  q_of_s <- vapply(specs, function(s) {
    resolve_pair_cohort_id(sel_fit, if (is.null(s$reference)) list(pair = s$pair, cohort = s$cohort) else
      list(pair = s$pair, cohort = s$cohort, reference = s$reference))
  }, integer(1))
  channels <- lapply(seq_len(S), function(s) up_map$pair_channels[[q_of_s[s]]])

  tcat <- dr_build_treatment_catalog(links, pools, specs, link_views, q_of_s = q_of_s)
  ts <- tcat$treatment_specs
  links <- tcat$links
  T_count <- length(ts)

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

  metadata <- unname(cp04_payload(list(specs, pools, links, ts, selection_specs, upstream, sel_fit, trt_recipe))["entries"])
  cohort_n <- vapply(specs, function(s) st$cohort_end[s$cohort] - st$cohort_start[s$cohort] + 1, numeric(1))
  or_charges <- or_corollary1_fit_charges(
    S, G, L, p_or, max(1L, up_map$d), D, st$n_A, min(st$N, ctrl$max_chunk_rows),
    max(cohort_n), sum(cohort_n), gc$maxit, metadata, header,
    sum(nchar(c(recipe$columns, trt_recipe$columns, if (is.null(reference)) colnames(st$Z) else
      reference$designs[[reference$role_slot["s"]]]$columns), type = "bytes"))
  )
  max_m <- max(vapply(ts, function(x) length(x$cohorts), integer(1)))
  max_trt_n <- max(vapply(ts, `[[`, integer(1), "n_observers"))
  trt_scan_n <- sum(vapply(ts, `[[`, integer(1), "n_observers"))
  trt_models_bound <- T_count * (p_trt * p_trt + 12 * (p_trt + p_or) + 8 * D + p_trt * D + 4096) +
    L * (8 * p_or + 8 * p_trt + 1024) + 8192
  trt_scratch <- max_trt_n * (p_trt + p_or + d + 16) + 32 * (p_trt + D + 4)^2 + 8192
  trt_work <- 128 * (gc$maxit + 2) * trt_scan_n * (p_trt + 2)^2 +
    128 * T_count * st$N * max_m + 128 * trt_models_bound

  scratch <- max(sel_fit$accounting$fit_scratch_bound, or_charges$scratch, trt_scratch)
  work <- ipw_owner_safe_sum(
    c(sel_fit$accounting$work_ops_bound, or_charges$work, trt_work),
    ctrl$max_work_ops,
    "Corollary 1 DR selection-models pipeline fit work"
  )
  bound <- metadata + or_charges$models + trt_models_bound
  cp04_budget(ctrl, c(base, bound), scratch, c(work, external_work_ops), bound, "Corollary 1 DR selection-models fits")

  supports <- lapply(seq_len(S), function(j) {
    dr_selection_support(reader, specs[[j]], sel_fit, up_map, q_of_s[j], recipe, or_frame_index, ctrl, reference)
  })
  state <- or_selection_assemble(
    recipe, supports, pools, links, channels, D, rank_tol, st$N,
    function(sid, visit) dr_support_visit(reader, specs[[sid]], ctrl, or_frame_index,
                                          function(rows, S_vec, obs, dy, Z) visit(Z), reference)
  )
  rm(supports)

  treatments <- lapply(ts, function(s) {
    dr_selection_treatment_fit(
      reader, s, sel_fit, up_map, state$regressions[[s$regression]],
      state$shared, recipe, or_frame_index, trt_recipe, gc, ctrl, reference
    )
  })
  links_out <- dr_populate_link_records(links, state$links, treatments, pools)

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
      or_frame_index = or_frame_index,
      treatment_recipe = trt_recipe,
      support_specs = specs,
      selection_specs = selection_specs,
      selections = selections,
      upstream = upstream,
      supports = state$supports,
      regression_specs = pools,
      regressions = state$regressions,
      shared = state$shared,
      treatment_specs = ts,
      treatments = treatments,
      links = links_out,
      universe = st$universe,
      justification = st$justification,
      law = paste0(
        "DR:Cor1:", st$universe, ":", selection_model, ":", selection_pooling,
        ":shared-one-period-slope:",
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
        n_supports = S,
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
        n_slope_columns = p_or,
        n_aliased_slopes = length(state$shared$basis$aliased_columns),
        n_support_scans = S * (1L + state$shared$verification_scan),
        n_treatment_fits = T_count,
        n_treatment_columns = p_trt,
        full_F_masks_per_model = 0L,
        full_F_by_link_IF_retained = FALSE,
        inverse_matrices_retained = FALSE,
        dense_dummy_regression = FALSE,
        peak_rss_measured = FALSE
      )
    ),
    class = c("cdid_dr_selection_backend", "cdid_dr_corollary1", "cdid_dr_backend")
  )
  if (!is.null(reference)) {
    obj$reference_ref <- reference$ref
    obj$accounting$n_reference_recipes <- nrow(reference$recipes)
  }
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  obj$accounting$owned_entries <- actual
  obj$accounting$retained_bound <- max(obj$accounting$retained_bound, base + actual)
  cp04_budget(ctrl, c(base, actual), scratch, c(work, external_work_ops), actual,
              "retained Corollary 1 DR selection-models")
  ipw_owner_seal(obj, "dr_corollary1_backend")
}

fit_dr_selection_models <- fit_dr_corollary1
fit_dr_conditional_mtar <- fit_dr_corollary1

dr_corollary1_if_envelope <- function(backend, B, links, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p_or <- backend$shared$p
  p_trt <- backend$treatment_recipe$p
  d <- max(1L, backend$upstream$d)
  G <- length(backend$regression_specs)
  J <- length(st$cohorts)
  sel_fit <- backend$selection_backend
  up_map <- backend$up_map

  sets <- or_pooled_if_sets(backend, links, details)
  gids <- sets$gids
  beta_gids <- sets$beta_gids
  tids <- sort(unique(vapply(backend$links[links], `[[`, integer(1), "model")))
  trt_sids <- sort(unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids"))))
  sids <- sort(unique(c(sets$sids, trt_sids)))
  uids <- if (details) {
    seq_len(up_map$D)
  } else {
    sort(unique(c(
      unlist(lapply(backend$supports[sids], `[[`, "channels")),
      unlist(lapply(backend$treatments[tids], `[[`, "upstream_ids"))
    )))
  }
  active_blocks <- if (details) {
    seq_along(sel_fit$blocks)
  } else {
    which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids), logical(1)))
  }
  active_cells <- unique(unlist(lapply(sel_fit$block_specs[active_blocks], `[[`, "cell_ids")))

  nl <- length(links)
  ns <- length(sids)
  nt <- length(tids)
  nu <- length(uids)
  nc <- length(active_cells)

  width <- 4 * nu + 4 * (p_or + 1) * (length(beta_gids) + 1) + 4 * length(gids) +
    nt * (2 * p_trt + nu + 16) + 8 * nl * (p_or + 4) +
    if (details) 3 * J + nu + nc + 8 * nl + 16 else 8
  metadata <- 16384 + 1024 * (ns + nu + G + nt + nl + J + length(sel_fit$blocks)) +
    B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  output <- B * (width + 16) + metadata
  scratch <- 64 * B * (width + p_or + p_trt + d + 8) +
    128 * ((p_or + 2)^2 + (p_trt + 2)^2 + (d + 2)^2) + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- 128 * B * (
    (ns + 1) * ((p_or + 2)^2 + (d + 2)^2) * log_a +
      (length(beta_gids) + length(gids)) * (p_or + 1) * (p_or + nu) +
      nt * ((p_trt + 2)^2 * log_a + p_trt * (nu + 4)) +
      J + nl * (p_or + p_trt + d + 4)
  ) + 64 * metadata

  sel_q_ids <- unique(backend$q_of_s[sids])
  selection_pre <- selection_if_budget(sel_fit, B, sel_q_ids, details = details)
  scratch <- scratch + selection_pre$output + selection_pre$scratch
  work <- work + selection_pre$work

  list(
    gids = gids, beta_gids = beta_gids, sids = sids, tids = tids,
    uids = uids, active_blocks = active_blocks, active_cells = active_cells,
    output = output, scratch = scratch, work = work
  )
}

dr_corollary1_if_budget <- function(backend, B, links, details = FALSE) {
  pre <- dr_corollary1_if_envelope(backend, B, links, details)
  cp04_budget(
    backend$control,
    c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    pre$scratch,
    c(backend$accounting$work_ops_bound, pre$work),
    pre$output,
    "Corollary 1 DR IF chunk"
  )
  pre
}

evaluate_dr_corollary1_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_dr_corollary1_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  valid_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) &&
    all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid_index(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    dr_corollary1_abort("Invalid explicit IF row/link batch, scale or details flag.")

  pre <- dr_corollary1_if_budget(backend, length(rows), links, details)
  B <- length(rows)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p_or <- backend$shared$p
  p_trt <- backend$treatment_recipe$p
  kept_or_cols <- if (p_or > 0L) backend$or_frame_index[backend$outcome_recipe$slope_index[backend$shared$keep]] else integer(0)
  trt_cols <- backend$treatment_recipe$index

  ref <- chained_backend_reference(backend)
  ci <- st$cohort_index[rows]
  # D60: no common weight; every support and treatment model reads its recipe.
  w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
  U_rows <- st$U[rows]
  A <- as.double(st$a_index[rows] > 0L)

  sel_fit <- backend$selection_backend
  up_map <- backend$up_map
  # 1. Evaluate upstream selection & rate IFs ONCE for the chunk
  up_res <- evaluate_selection_upstream_chunk(sel_fit, up_map, rows, pre$uids, details = details)
  up <- up_res$up

  # D92: the non-detailed Corollary 1 DR IF reads only core$beta (psi_or, the
  # existing zero vector at p_or = 0). A non-detailed chunk therefore requests
  # the beta-only core at p_or > 0 and skips the core and its support outcome
  # reads at p_or = 0; each treatment model still reads its own pair endpoints
  # and selection weights below. Detailed chunks keep the full core.
  core_used <- details || p_or > 0L
  beta_only <- !details

  hits <- vector("list", length(pre$sids))
  if (core_used) {
    for (j in seq_along(pre$sids)) {
      sid <- pre$sids[j]
      sp <- backend$support_specs[[sid]]
      q <- backend$q_of_s[sid]
      ix <- which(ci == sp$cohort & U_rows)
      hit <- cp04_pair_rows(reader, rows[ix], sp$pair, TRUE)
      obs <- which(hit$S == 1)
      rr <- rows[ix[obs]]
      pr_obs <- if (length(rr)) predict_pair_selection(sel_fit, rr, q)$prob else numeric(0)
      hits[[j]] <- if (is.null(ref)) list(
        obs = ix[obs],
        dy = hit$dy[obs],
        z = if (p_or > 0L) frame$Z_a[st$a_index[rr], kept_or_cols, drop = FALSE] else matrix(0, length(obs), 0L),
        v = or_corollary1_factor(st$a[sp$cohort], st$rho[rr], pr_obs)
      ) else list(
        obs = ix[obs],
        dy = hit$dy[obs],
        z = if (p_or > 0L) chained_reference_x(ref, sp$reference, rr)[, kept_or_cols, drop = FALSE] else
          matrix(0, length(obs), 0L),
        v = or_corollary1_factor(chained_reference_stat(ref, sp$reference, sp$cohort, "a"), st$rho[rr], pr_obs),
        w = chained_reference_weight(ref, sp$reference, rr)
      )
    }
  }

  # 2. Shared OR IF core (global xi_beta; full core when detailed)
  core <- if (core_used) or_pooled_if_core(backend, st$N, w, hits, up, pre, links, details, beta_only = beta_only) else NULL

  # 3. Treatment propensity roots and residualized link components
  models <- vector("list", length(backend$treatments))
  for (tj in pre$tids) {
    spec <- backend$treatment_specs[[tj]]
    f <- backend$treatments[[tj]]
    reg <- backend$regressions[[spec$regression]]
    u_pos <- match(f$upstream_ids, pre$uids)
    upstream <- up[, u_pos, drop = FALSE]

    ci_in <- which(ci %in% spec$cohorts)
    hit_all <- if (length(ci_in)) cp04_pair_rows(reader, rows[ci_in], spec$pair, TRUE) else list(S = numeric(0), dy = numeric(0))
    obs_local <- which(hit_all$S == 1)
    rr <- ci_in[obs_local]

    beta_score <- matrix(0, B, p_trt)
    direct <- matrix(0, B, 4L)
    direct[, 3L] <- -f$means["dp"]
    direct[, 4L] <- -f$means["dm"]
    direct_ipw <- if (details) {
      d_ipw <- matrix(0, B, 4L)
      d_ipw[, 3L] <- -f$means["dp"]
      d_ipw[, 4L] <- -f$means["dm"]
      d_ipw
    } else NULL

    if (length(rr)) {
      a_idx <- st$a_index[rows[rr]]
      if (is.null(ref)) {
        X_trt <- frame$Z_a[a_idx, trt_cols, drop = FALSE]
        Z_or <- if (p_or > 0L) frame$Z_a[a_idx, kept_or_cols, drop = FALSE] else matrix(0, length(rr), 0L)
      } else {
        Xr <- chained_reference_x(ref, spec$reference, rows[rr])
        X_trt <- Xr[, trt_cols, drop = FALSE]
        Z_or <- if (p_or > 0L) Xr[, kept_or_cols, drop = FALSE] else matrix(0, length(rr), 0L)
      }
      dy_rr <- hit_all$dy[obs_local]
      d_rr <- reg$L * Z_or
      m_hat_rr <- (if (backend$outcome_recipe$intercept) reg$alpha else 0) +
        (if (p_or > 0L) as.vector(d_rr %*% backend$shared$beta) else 0)
      r_rr <- dy_rr - m_hat_rr

      lv <- ipw_logit_values(X_trt, f$coefficients)
      cp04_finite(c(lv$prob, lv$eta, lv$odds), "Corollary 1 DR treatment replay")
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
      beta_score[rr, ] <- fw * (D_ind - lv$prob) * X_trt
      direct[rr, 1L] <- wp * (r_rr - f$means["T_dr"]) / f$means["dp"]
      direct[rr, 2L] <- wm * (r_rr - f$means["M_dr"]) / f$means["dm"]
      direct[rr, 3L] <- direct[rr, 3L] + wp
      direct[rr, 4L] <- direct[rr, 4L] + wm
      if (details) {
        direct_ipw[rr, 1L] <- wp * (dy_rr - f$means["T_ipw"]) / f$means["dp"]
        direct_ipw[rr, 2L] <- wm * (dy_rr - f$means["M_ipw"]) / f$means["dm"]
        direct_ipw[rr, 3L] <- direct_ipw[rr, 3L] + wp
        direct_ipw[rr, 4L] <- direct_ipw[rr, 4L] + wm
      }
    }

    xi_gamma <- t(cp04_solve(f$H_chol, t(beta_score + upstream %*% t(f$cross))))
    colnames(xi_gamma) <- backend$treatment_recipe$columns
    psi_dir_m <- direct[, 1L] - direct[, 2L]
    psi_or_m <- if (p_or > 0L) -as.vector(core$beta %*% f$d_bar_diff) else numeric(B)
    psi_ps_m <- -as.vector(xi_gamma %*% f$partials$beta[2L, ])
    psi_sel_m <- as.vector(upstream %*% (f$partials$upstream[1L, ] - f$partials$upstream[2L, ]))
    psi_dr_m <- psi_dir_m + psi_or_m + psi_ps_m + psi_sel_m
    coord <- direct + upstream %*% t(f$partials$upstream) + xi_gamma %*% t(f$partials$beta)
    colnames(coord) <- c("T", "M", "dp", "dm")
    cp04_finite(c(xi_gamma, coord, psi_dr_m), "complete Corollary 1 DR frame IF")

    if (details) {
      coord_ipw <- direct_ipw + upstream %*% t(f$partials_ipw$upstream) + xi_gamma %*% t(f$partials_ipw$beta)
      colnames(coord_ipw) <- c("T", "M", "dp", "dm")
      models[[tj]] <- list(
        beta = xi_gamma,
        gamma = xi_gamma,
        coordinates = coord,
        direct = direct,
        beta_score = beta_score,
        upstream = upstream,
        psi_dr = psi_dr_m,
        psi_dir = psi_dir_m,
        psi_or = psi_or_m,
        psi_ps = psi_ps_m,
        psi_sel = psi_sel_m,
        psi_ipw = coord_ipw[, 1L] - coord_ipw[, 2L],
        normalizers_A = coord[, 3:4, drop = FALSE] / st$lambda -
          outer(A - st$lambda, f$means[c("dp", "dm")] / st$lambda^2)
      )
    } else {
      models[[tj]] <- list(
        phi = psi_dr_m,
        psi_dir = psi_dir_m,
        psi_or = psi_or_m,
        psi_ps = psi_ps_m,
        psi_sel = psi_sel_m
      )
    }
  }

  link_names <- vapply(backend$links[links], `[[`, character(1), "link_id")
  phi <- matrix(0, B, length(links), dimnames = list(NULL, link_names))
  if (details) {
    psi_dir_mat <- psi_or_mat <- psi_ps_mat <- psi_sel_mat <- psi_ipw_mat <-
      matrix(0, B, length(links), dimnames = list(NULL, link_names))
  }
  for (j in seq_along(links)) {
    m <- models[[backend$links[[links[j]]]$model]]
    phi[, j] <- if (details) m$psi_dr else m$phi
    if (details) {
      psi_dir_mat[, j] <- m$psi_dir
      psi_or_mat[, j] <- m$psi_or
      psi_ps_mat[, j] <- m$psi_ps
      psi_sel_mat[, j] <- m$psi_sel
      psi_ipw_mat[, j] <- m$psi_ipw
    }
  }
  cp04_finite(phi, "complete Corollary 1 DR link IF")
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
      n_regression_score_evaluations = if (core_used) length(union(pre$gids, pre$beta_gids)) else 0L,
      n_global_beta_blocks = if (core_used) length(pre$beta_gids) else 0L,
      n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  if (details) {
    # D60 rates/shares belong to each recipe (in `upstream`); frame-level
    # a/tau/eta summaries exist only on the static path.
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
    alpha_sel_list <- lapply(up_res$block_ifs, function(bif) if (is.null(bif)) NULL else bif$xi_theta)
    out$details <- list(
      support_ids = pre$sids,
      upstream_ids = pre$uids,
      upstream = up,
      a = a_if,
      alpha_selection = alpha_sel_list,
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
      treatment_ids = pre$tids,
      models = models,
      psi_dir = psi_dir_mat,
      psi_or = psi_or_mat,
      psi_ps = psi_ps_mat,
      psi_sel = psi_sel_mat,
      psi_ipw = psi_ipw_mat,
      psi_or_estimator = core$phi,
      tau = tau,
      eta = eta,
      lambda = A - st$lambda
    )
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    dr_corollary1_abort("IF output bound underestimated.", "bound")
  out
}

replay_dr_corollary1_if <- function(backend, chunk_rows, link_batch_size, callback,
                                    scale = "frame", details = FALSE) {
  validate_dr_corollary1_backend(backend)
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  LB <- ipw_owner_int_bound(link_batch_size, "link_batch_size")
  if (B > backend$control$max_chunk_rows || !is.function(callback) ||
      !identical(details, isTRUE(details)) || !is.character(scale) || length(scale) != 1L ||
      is.na(scale) || !scale %in% c("frame", "core"))
    dr_corollary1_abort("Invalid explicit replay controls/callback.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N
  L <- length(backend$links)
  total <- backend$accounting$work_ops_bound
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    pre <- dr_corollary1_if_budget(backend, min(B, N), ll, details)
    total <- ipw_owner_safe_sum(
      c(total, ceiling(N / B) * pre$work),
      backend$control$max_work_ops,
      "total Corollary 1 DR pipeline replay work"
    )
  }
  batches <- 0L
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1L))
      batch <- evaluate_dr_corollary1_if_chunk(backend, rr, ll, scale, details)
      callback(batch)
      rm(batch)
      batches <- batches + 1L
    }
  }
  invisible(list(
    n_frame_rows = N, n_links = L, batches = batches, scale = scale,
    normalization_count = if (scale == "core") reader$state$n_A else N,
    work_ops_bound = total, all_F_identities = TRUE
  ))
}

# ------------------------------------------------------------------------------
# Unified Private DR Dispatchers & Covariance Streaming
# ------------------------------------------------------------------------------

fit_dr_backend <- function(reader, link_views, control = ipw_owner_control(),
                           glm_control, rank_tol = or_public_rank_tol(),
                           outcome_recipe, treatment_columns = NULL,
                           selection = "none",
                           selection_model = "pair_logit",
                           selection_pooling = "separate",
                           retained_external_entries = 0, external_work_ops = 0,
                           automatic_chunk_rows = FALSE,
                           sel_fit = NULL) {
  if (!is.character(selection) || length(selection) != 1L || is.na(selection) ||
      !selection %in% c("none", "baseline", "conditional_mtar", "corollary1")) {
    dr_baseline_abort('`selection` must be "none" ("baseline") or "conditional_mtar" ("corollary1").', "input")
  }
  if (selection %in% c("none", "baseline")) {
    fit_dr_baseline(
      reader = reader, link_views = link_views, control = control,
      glm_control = glm_control, rank_tol = rank_tol,
      outcome_recipe = outcome_recipe, treatment_columns = treatment_columns,
      retained_external_entries = retained_external_entries,
      external_work_ops = external_work_ops,
      automatic_chunk_rows = automatic_chunk_rows
    )
  } else {
    fit_dr_corollary1(
      reader = reader, link_views = link_views, control = control,
      glm_control = glm_control, rank_tol = rank_tol,
      outcome_recipe = outcome_recipe, treatment_columns = treatment_columns,
      retained_external_entries = retained_external_entries,
      external_work_ops = external_work_ops,
      automatic_chunk_rows = automatic_chunk_rows,
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      sel_fit = sel_fit
    )
  }
}

dr_if_envelope <- function(backend, B, links, details = FALSE) {
  if (inherits(backend, "cdid_dr_baseline")) return(dr_baseline_if_envelope(backend, B, links, details))
  if (inherits(backend, "cdid_dr_corollary1")) return(dr_corollary1_if_envelope(backend, B, links, details))
  dr_baseline_abort("Expected a private DR backend.", "integrity")
}

dr_if_budget <- function(backend, B, links, details = FALSE) {
  if (inherits(backend, "cdid_dr_baseline")) return(dr_baseline_if_budget(backend, B, links, details))
  if (inherits(backend, "cdid_dr_corollary1")) return(dr_corollary1_if_budget(backend, B, links, details))
  dr_baseline_abort("Expected a private DR backend.", "integrity")
}

evaluate_dr_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  if (inherits(backend, "cdid_dr_baseline")) return(evaluate_dr_baseline_if_chunk(backend, rows, links, scale, details))
  if (inherits(backend, "cdid_dr_corollary1")) return(evaluate_dr_corollary1_if_chunk(backend, rows, links, scale, details))
  dr_baseline_abort("Expected a private DR backend.", "integrity")
}

replay_dr_if <- function(backend, chunk_rows, link_batch_size, callback,
                         scale = "frame", details = FALSE) {
  if (inherits(backend, "cdid_dr_baseline"))
    return(replay_dr_baseline_if(backend, chunk_rows, link_batch_size, callback, scale, details))
  if (inherits(backend, "cdid_dr_corollary1"))
    return(replay_dr_corollary1_if(backend, chunk_rows, link_batch_size, callback, scale, details))
  dr_baseline_abort("Expected a private DR backend.", "integrity")
}

stream_dr_covariance <- function(backend, links = seq_along(backend$links), retained = 0) {
  validate_dr_backend(backend)
  cp04_stream_covariance(backend, links, retained = retained)
}

# ------------------------------------------------------------------------------
# Pure-Helper Table Entry Points for Chained DR
# (Exercise the exact same selection/rate, OR WLS, treatment logit, and 4-term IF
#  helpers without requiring a C owner when testing pure table inputs.)
# ------------------------------------------------------------------------------

dr_table_pair_dy <- function(Y, pair_outcomes, pair, obs_rows) {
  if (!is.null(pair_outcomes)) {
    k1 <- paste0(pair$s, ":", pair$t)
    k2 <- paste0(pair$s_calendar_index, ":", pair$t_calendar_index)
    vec <- if (!is.null(pair_outcomes[[k1]])) pair_outcomes[[k1]] else pair_outcomes[[k2]]
    if (!is.null(vec)) return(as.double(vec[obs_rows]))
  }
  as.double(Y[obs_rows, pair$t_calendar_index] - Y[obs_rows, pair$s_calendar_index])
}

fit_dr_baseline_table <- function(X_or, V, Y = NULL, obs, g, w,
                                  calendar, links, outcome_recipe,
                                  pair_outcomes = NULL,
                                  glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
                                  rank_tol = 1e-10) {
  N <- nrow(X_or)
  cohorts <- sort(unique(as.double(g)), method = "radix")
  ci <- match(as.double(g), cohorts)
  A <- as.double(rowSums(obs != 0) >= 2L)
  n_A <- sum(A > 0)
  mean_w <- sum(w) / N
  mean_wA <- sum(w * A) / N
  tau_hat <- eta_hat <- numeric(length(cohorts))
  tau_if <- eta_if <- matrix(0, N, length(cohorts))
  for (j in seq_along(cohorts)) {
    I_j <- as.double(ci == j)
    tau_hat[j] <- (sum(w * I_j) / N) / mean_w
    eta_hat[j] <- (sum(w * A * I_j) / N) / mean_wA
    tau_if[, j] <- w * (I_j - tau_hat[j]) / mean_w
    eta_if[, j] <- w * A * (I_j - eta_hat[j]) / mean_wA
  }

  recipe <- or_pooled_check_recipe(outcome_recipe, colnames(X_or), "baseline")
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
      functional_key = paste0("cmpcat:DR-Base:", lk$link_id)
    )
  }
  catalog <- or_baseline_catalog(records, cohorts, recipe)
  specs <- catalog$specs
  pools <- catalog$pools
  cat_links <- catalog$links
  S <- length(specs)
  S_mat <- matrix(0, N, S)
  supports_raw <- vector("list", S)
  rate_if <- matrix(0, N, S)
  for (s in seq_len(S)) {
    sp <- specs[[s]]
    pr <- sp$pair
    S_mat[, s] <- as.double(obs[, pr$s_calendar_index] != 0 & obs[, pr$t_calendar_index] != 0)
    obs_rows <- which(ci == sp$cohort & S_mat[, s] == 1)
    dy <- dr_table_pair_dy(Y, pair_outcomes, pr, obs_rows)
    Z_obs <- X_or[obs_rows, , drop = FALSE]
    A_k <- ci == sp$cohort & A > 0
    repeat_mass <- sum(w[A_k])
    n_repeat <- sum(A_k)
    acc <- or_pooled_support_start(recipe, sp$transitions)
    acc <- or_pooled_support_add(acc, recipe, Z_obs, dy, w[obs_rows], "baseline")
    supports_raw[[s]] <- or_baseline_support_rate(acc, repeat_mass, n_repeat, N)
    I_k <- as.double(ci == sp$cohort)
    rate_if[, s] <- (w * I_k * (S_mat[, s] - supports_raw[[s]]$rate * A)) / (repeat_mass / N)
  }

  state <- or_baseline_assemble(
    recipe, supports_raw, pools, cat_links, rank_tol, N,
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
    links = state$links
  )
  L <- length(state$links)
  sets <- or_pooled_if_sets(fake_backend, seq_len(L), details = TRUE)
  pre <- list(gids = sets$gids, beta_gids = sets$beta_gids, sids = sets$sids)
  p_or <- state$shared$p
  kept <- if (p_or > 0L) recipe$slope_index[state$shared$keep] else integer(0)
  hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]
    sp <- specs[[sid]]
    obs_rows <- which(ci == sp$cohort & S_mat[, sid] == 1)
    dy <- dr_table_pair_dy(Y, pair_outcomes, sp$pair, obs_rows)
    hits[[j]] <- list(
      obs = obs_rows,
      dy = dy,
      z = if (p_or > 0L) X_or[obs_rows, kept, drop = FALSE] else matrix(0, length(obs_rows), 0L),
      v = rep(state$supports[[sid]]$v, length(obs_rows))
    )
  }
  core <- or_pooled_if_core(fake_backend, N, w, hits, rate_if[, pre$sids, drop = FALSE], pre, seq_len(L), details = TRUE)

  p_trt <- ncol(V)
  link_names <- vapply(state$links, `[[`, character(1), "link_id")
  phi <- psi_dir <- psi_or <- psi_ps <- psi_sel <- psi_ipw <-
    matrix(0, N, L, dimnames = list(NULL, link_names))
  treatments <- vector("list", L)
  gc <- validate_ipw_glm_control(glm_control)

  for (ell in seq_len(L)) {
    l_or <- state$links[[ell]]
    reg <- state$regressions[[l_or$regression]]
    s_plus <- l_or$support_id
    s_minus <- pools[[l_or$regression]]$support_ids
    s_all <- c(s_plus, s_minus)
    coh_idx <- c(l_or$cohort, pools[[l_or$regression]]$cohorts)
    pair <- specs[[s_plus]]$pair
    S_pair <- as.double(obs[, pair$s_calendar_index] != 0 & obs[, pair$t_calendar_index] != 0)
    obs_rows <- which(ci %in% coh_idx & S_pair == 1)
    ci_obs <- match(ci[obs_rows], coh_idx)

    X_t <- V[obs_rows, , drop = FALSE]
    Z_o <- if (p_or > 0L) X_or[obs_rows, kept, drop = FALSE] else matrix(0, length(obs_rows), 0L)
    y_t <- as.double(ci_obs == 1L)
    dy_obs <- dr_table_pair_dy(Y, pair_outcomes, pair, obs_rows)
    d_obs <- reg$L * Z_o
    m_hat <- (if (recipe$intercept) reg$alpha else 0) +
      (if (p_or > 0L) as.vector(d_obs %*% state$shared$beta) else 0)
    r_obs <- dy_obs - m_hat

    m <- length(coh_idx)
    rates <- vapply(s_all, function(sid) state$supports[[sid]]$rate, numeric(1))
    fw <- w[obs_rows] / ipw_clip_floor(rates)[ci_obs]

    fit <- cp04_logit(X_t, y_t, fw, gc, N, list(pair = pair, cohorts = coh_idx))
    lv <- ipw_logit_values(X_t, fit$coefficients)
    wp <- fw * y_t
    wm <- fw * (1 - y_t) * lv$odds
    dp <- sum(wp) / N
    dm <- sum(wm) / N
    T_dr <- sum(wp * r_obs) / sum(wp)
    M_dr <- sum(wm * r_obs) / sum(wm)
    T_ipw <- sum(wp * dy_obs) / sum(wp)
    M_ipw <- sum(wm * dy_obs) / sum(wm)
    d_bar_plus <- if (p_or > 0L) as.vector(colSums(wp * d_obs) / sum(wp)) else numeric(0)
    d_bar_minus <- if (p_or > 0L) as.vector(colSums(wm * d_obs) / sum(wm)) else numeric(0)
    d_bar_diff <- d_bar_plus - d_bar_minus

    C <- matrix(0, p_trt, m)
    U_dr <- matrix(0, 4L, m, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    U_ipw <- matrix(0, 4L, m, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    score <- fw * (y_t - lv$prob) * X_t
    residual_dr <- cbind(wp * (r_obs - T_dr) / dp, wm * (r_obs - M_dr) / dm, wp, wm)
    residual_ipw <- cbind(wp * (dy_obs - T_ipw) / dp, wm * (dy_obs - M_ipw) / dm, wp, wm)

    for (k in seq_len(m)) {
      jj <- which(ci_obs == k)
      scale_k <- if (rates[k] > ipw_clip_eps) -1 / (N * rates[k]) else 0
      C[, k] <- colSums(score[jj, , drop = FALSE]) * scale_k
      U_dr[, k] <- colSums(residual_dr[jj, , drop = FALSE]) * scale_k
      U_ipw[, k] <- colSums(residual_ipw[jj, , drop = FALSE]) * scale_k
    }
    wm_act <- wm * lv$active
    V_dr <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    V_dr[2L, ] <- colSums(wm_act * (r_obs - M_dr) * X_t) / (N * dm)
    V_dr[4L, ] <- colSums(wm_act * X_t) / N
    V_ipw <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    V_ipw[2L, ] <- colSums(wm_act * (dy_obs - M_ipw) * X_t) / (N * dm)
    V_ipw[4L, ] <- V_dr[4L, ]

    beta_score <- matrix(0, N, p_trt)
    direct <- direct_ipw <- matrix(0, N, 4L)
    direct[, 3L] <- direct_ipw[, 3L] <- -dp
    direct[, 4L] <- direct_ipw[, 4L] <- -dm
    beta_score[obs_rows, ] <- score
    direct[obs_rows, 1L] <- residual_dr[, 1L]
    direct[obs_rows, 2L] <- residual_dr[, 2L]
    direct[obs_rows, 3L] <- direct[obs_rows, 3L] + wp
    direct[obs_rows, 4L] <- direct[obs_rows, 4L] + wm
    direct_ipw[obs_rows, 1L] <- residual_ipw[, 1L]
    direct_ipw[obs_rows, 2L] <- residual_ipw[, 2L]
    direct_ipw[obs_rows, 3L] <- direct_ipw[obs_rows, 3L] + wp
    direct_ipw[obs_rows, 4L] <- direct_ipw[obs_rows, 4L] + wm

    up_sub <- rate_if[, s_all, drop = FALSE]
    xi_gamma <- t(cp04_solve(fit$H_chol, t(beta_score + up_sub %*% t(C))))
    psi_dir[, ell] <- direct[, 1L] - direct[, 2L]
    psi_or[, ell] <- if (p_or > 0L) -as.vector(core$beta %*% d_bar_diff) else numeric(N)
    psi_ps[, ell] <- -as.vector(xi_gamma %*% V_dr[2L, ])
    psi_sel[, ell] <- as.vector(up_sub %*% (U_dr[1L, ] - U_dr[2L, ]))
    phi[, ell] <- psi_dir[, ell] + psi_or[, ell] + psi_ps[, ell] + psi_sel[, ell]
    coord_ipw <- direct_ipw + up_sub %*% t(U_ipw) + xi_gamma %*% t(V_ipw)
    psi_ipw[, ell] <- coord_ipw[, 1L] - coord_ipw[, 2L]

    fit$means <- c(
      T = T_dr, M = M_dr, dp = dp, dm = dm, mu = T_dr - M_dr,
      T_dr = T_dr, M_dr = M_dr, mu_dr = T_dr - M_dr,
      T_ipw = T_ipw, M_ipw = M_ipw, mu_ipw = T_ipw - M_ipw,
      mu_or = l_or$mu
    )
    fit$d_bar_plus <- d_bar_plus
    fit$d_bar_minus <- d_bar_minus
    fit$d_bar_diff <- d_bar_diff
    fit$support_ids <- s_all
    fit$cross <- C
    fit$partials <- list(upstream = U_dr, beta = V_dr)
    fit$beta_if <- xi_gamma
    treatments[[ell]] <- fit
  }

  list(
    supports = state$supports,
    support_specs = specs,
    rate_if = rate_if,
    tau_hat = tau_hat,
    eta_hat = eta_hat,
    tau_if = tau_if,
    eta_if = eta_if,
    shared = state$shared,
    regressions = state$regressions,
    treatments = treatments,
    links = state$links,
    estimates = vapply(treatments, function(f) unname(f$means["mu_dr"]), numeric(1)),
    estimates_or = vapply(state$links, `[[`, numeric(1), "mu"),
    estimates_ipw = vapply(treatments, function(f) unname(f$means["mu_ipw"]), numeric(1)),
    phi = phi,
    phi_core = phi * (n_A / N),
    psi_dir = psi_dir,
    psi_or = psi_or,
    psi_ps = psi_ps,
    psi_sel = psi_sel,
    psi_ipw = psi_ipw,
    psi_or_estimator = core$phi,
    beta_if = core$beta,
    alpha_if = core$alpha
  )
}

fit_dr_selection_table <- function(Z, X_or, V, Y = NULL, obs, g, w,
                                   U = rep(TRUE, nrow(Z)), rho = rep(1, nrow(Z)),
                                   universe = NULL, calendar, links, outcome_recipe,
                                   pair_outcomes = NULL,
                                   selection_model = "pair_logit",
                                   selection_pooling = "separate",
                                   has_intercept = NULL,
                                   glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
                                   rank_tol = 1e-10, chunk_rows = nrow(Z),
                                   sel_fit = NULL) {
  N <- nrow(X_or)
  cohorts <- sort(unique(as.double(g)), method = "radix")
  ci <- match(as.double(g), cohorts)
  A <- as.double(rowSums(obs != 0) >= 2L)
  n_A <- sum(A > 0)
  mean_w <- sum(w) / N
  mean_wA <- sum(w * A) / N
  tau_hat <- eta_hat <- numeric(length(cohorts))
  tau_if <- eta_if <- matrix(0, N, length(cohorts))
  for (j in seq_along(cohorts)) {
    I_j <- as.double(ci == j)
    tau_hat[j] <- (sum(w * I_j) / N) / mean_w
    eta_hat[j] <- (sum(w * A * I_j) / N) / mean_wA
    tau_if[, j] <- w * (I_j - tau_hat[j]) / mean_w
    eta_if[, j] <- w * A * (I_j - eta_hat[j]) / mean_wA
  }
  if (is.null(universe)) {
    universe <- if (all(U)) "F" else if (identical(as.logical(U), as.logical(A > 0))) "A" else "B"
  }
  cp04_finite(Z[U, , drop = FALSE], "selection design on U")
  cp04_finite(rho[U], "known normalized rho on U", TRUE)
  if (any(!U)) {
    Z[!U, ] <- 0
    rho[!U] <- 1
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
      functional_key = paste0("cmpcat:DR-Cor1:", lk$link_id)
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

  if (is.null(sel_fit)) {
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
  }

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
    obs_rows <- which(ci == sp$cohort & S_mat[, s] == 1)
    dy <- dr_table_pair_dy(Y, pair_outcomes, sp$pair, obs_rows)
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
  uids_or <- sort(unique(unlist(lapply(state$supports[sets$sids], `[[`, "channels"))))
  pre <- list(gids = sets$gids, beta_gids = sets$beta_gids, sids = sets$sids, uids = uids_or)
  p_or <- state$shared$p
  kept <- if (p_or > 0L) recipe$slope_index[state$shared$keep] else integer(0)
  hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]
    sp <- specs[[sid]]
    q <- q_of_s[sid]
    obs_rows <- which(ci == sp$cohort & S_mat[, sid] == 1)
    pred <- predict_table_pair_selection(sel_fit, Z, w, rho, obs_rows, q)
    hits[[j]] <- list(
      obs = obs_rows,
      dy = dr_table_pair_dy(Y, pair_outcomes, sp$pair, obs_rows),
      z = if (p_or > 0L) X_or[obs_rows, kept, drop = FALSE] else matrix(0, length(obs_rows), 0L),
      v = or_corollary1_factor(sel_fit$a[sp$cohort], rho[obs_rows], pred$prob)
    )
  }
  core <- or_pooled_if_core(fake_backend, N, w, hits, up_all[, uids_or, drop = FALSE], pre, seq_len(L), details = TRUE)

  p_trt <- ncol(V)
  link_names <- vapply(state$links, `[[`, character(1), "link_id")
  phi <- psi_dir <- psi_or <- psi_ps <- psi_sel <- psi_ipw <-
    matrix(0, N, L, dimnames = list(NULL, link_names))
  treatments <- vector("list", L)
  gc <- validate_ipw_glm_control(glm_control)

  for (ell in seq_len(L)) {
    l_or <- state$links[[ell]]
    reg <- state$regressions[[l_or$regression]]
    s_plus <- l_or$support_id
    s_minus <- pools[[l_or$regression]]$support_ids
    s_all <- c(s_plus, s_minus)
    coh_idx <- c(l_or$cohort, pools[[l_or$regression]]$cohorts)
    q_ids <- q_of_s[s_all]
    pair <- specs[[s_plus]]$pair
    s_idx <- pair$s_calendar_index
    t_idx <- pair$t_calendar_index
    S_pair <- as.double(obs[, s_idx] != 0 & obs[, t_idx] != 0)
    obs_rows <- which(ci %in% coh_idx & S_pair == 1)
    ci_obs <- match(ci[obs_rows], coh_idx)

    X_t <- V[obs_rows, , drop = FALSE]
    Z_o <- if (p_or > 0L) X_or[obs_rows, kept, drop = FALSE] else matrix(0, length(obs_rows), 0L)
    y_t <- as.double(ci_obs == 1L)
    dy_obs <- dr_table_pair_dy(Y, pair_outcomes, pair, obs_rows)
    d_obs <- reg$L * Z_o
    m_hat <- (if (recipe$intercept) reg$alpha else 0) +
      (if (p_or > 0L) as.vector(d_obs %*% state$shared$beta) else 0)
    r_obs <- dy_obs - m_hat

    fw <- numeric(length(obs_rows))
    m <- length(coh_idx)
    ch_by_k <- vector("list", m)
    ld_by_k <- vector("list", m)
    for (k in seq_len(m)) {
      jj <- which(ci_obs == k)
      q <- q_ids[k]
      pred <- predict_table_pair_selection(sel_fit, Z, w, rho, obs_rows[jj], q)
      ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, q, obs_rows[jj], pred = pred, Z = Z, a_vec = sel_fit$a)
      fw[jj] <- pred$weight
      ch_by_k[[k]] <- up_map$pair_channels[[q]]
      ld_by_k[[k]] <- ld$log_deriv
    }
    uids_t <- sort(unique(unlist(ch_by_k)))
    u_len <- length(uids_t)
    fit <- cp04_logit(X_t, y_t, fw, gc, N, list(pair = pair, cohorts = coh_idx))
    lv <- ipw_logit_values(X_t, fit$coefficients)
    wp <- fw * y_t
    wm <- fw * (1 - y_t) * lv$odds
    dp <- sum(wp) / N
    dm <- sum(wm) / N
    T_dr <- sum(wp * r_obs) / sum(wp)
    M_dr <- sum(wm * r_obs) / sum(wm)
    T_ipw <- sum(wp * dy_obs) / sum(wp)
    M_ipw <- sum(wm * dy_obs) / sum(wm)
    d_bar_plus <- if (p_or > 0L) as.vector(colSums(wp * d_obs) / sum(wp)) else numeric(0)
    d_bar_minus <- if (p_or > 0L) as.vector(colSums(wm * d_obs) / sum(wm)) else numeric(0)
    d_bar_diff <- d_bar_plus - d_bar_minus

    C <- matrix(0, p_trt, u_len)
    U_dr <- matrix(0, 4L, u_len, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    U_ipw <- matrix(0, 4L, u_len, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    score <- fw * (y_t - lv$prob) * X_t
    residual_dr <- cbind(wp * (r_obs - T_dr) / dp, wm * (r_obs - M_dr) / dm, wp, wm)
    residual_ipw <- cbind(wp * (dy_obs - T_ipw) / dp, wm * (dy_obs - M_ipw) / dm, wp, wm)

    for (k in seq_len(m)) {
      jj <- which(ci_obs == k)
      pos_u <- match(ch_by_k[[k]], uids_t)
      ld <- ld_by_k[[k]]
      C[, pos_u] <- C[, pos_u] + crossprod(score[jj, , drop = FALSE], ld) / N
      U_dr[, pos_u] <- U_dr[, pos_u] + crossprod(residual_dr[jj, , drop = FALSE], ld) / N
      U_ipw[, pos_u] <- U_ipw[, pos_u] + crossprod(residual_ipw[jj, , drop = FALSE], ld) / N
    }
    wm_act <- wm * lv$active
    V_dr <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    V_dr[2L, ] <- colSums(wm_act * (r_obs - M_dr) * X_t) / (N * dm)
    V_dr[4L, ] <- colSums(wm_act * X_t) / N
    V_ipw <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), NULL))
    V_ipw[2L, ] <- colSums(wm_act * (dy_obs - M_ipw) * X_t) / (N * dm)
    V_ipw[4L, ] <- V_dr[4L, ]

    beta_score <- matrix(0, N, p_trt)
    direct <- direct_ipw <- matrix(0, N, 4L)
    direct[, 3L] <- direct_ipw[, 3L] <- -dp
    direct[, 4L] <- direct_ipw[, 4L] <- -dm
    beta_score[obs_rows, ] <- score
    direct[obs_rows, 1L] <- residual_dr[, 1L]
    direct[obs_rows, 2L] <- residual_dr[, 2L]
    direct[obs_rows, 3L] <- direct[obs_rows, 3L] + wp
    direct[obs_rows, 4L] <- direct[obs_rows, 4L] + wm
    direct_ipw[obs_rows, 1L] <- residual_ipw[, 1L]
    direct_ipw[obs_rows, 2L] <- residual_ipw[, 2L]
    direct_ipw[obs_rows, 3L] <- direct_ipw[obs_rows, 3L] + wp
    direct_ipw[obs_rows, 4L] <- direct_ipw[obs_rows, 4L] + wm

    up_sub <- up_all[, uids_t, drop = FALSE]
    xi_gamma <- t(cp04_solve(fit$H_chol, t(beta_score + up_sub %*% t(C))))
    psi_dir[, ell] <- direct[, 1L] - direct[, 2L]
    psi_or[, ell] <- if (p_or > 0L) -as.vector(core$beta %*% d_bar_diff) else numeric(N)
    psi_ps[, ell] <- -as.vector(xi_gamma %*% V_dr[2L, ])
    psi_sel[, ell] <- as.vector(up_sub %*% (U_dr[1L, ] - U_dr[2L, ]))
    phi[, ell] <- psi_dir[, ell] + psi_or[, ell] + psi_ps[, ell] + psi_sel[, ell]
    coord_ipw <- direct_ipw + up_sub %*% t(U_ipw) + xi_gamma %*% t(V_ipw)
    psi_ipw[, ell] <- coord_ipw[, 1L] - coord_ipw[, 2L]

    fit$means <- c(
      T = T_dr, M = M_dr, dp = dp, dm = dm, mu = T_dr - M_dr,
      T_dr = T_dr, M_dr = M_dr, mu_dr = T_dr - M_dr,
      T_ipw = T_ipw, M_ipw = M_ipw, mu_ipw = T_ipw - M_ipw,
      mu_or = l_or$mu
    )
    fit$d_bar_plus <- d_bar_plus
    fit$d_bar_minus <- d_bar_minus
    fit$d_bar_diff <- d_bar_diff
    fit$upstream_ids <- uids_t
    fit$cross <- C
    fit$partials <- list(upstream = U_dr, beta = V_dr)
    fit$beta_if <- xi_gamma
    treatments[[ell]] <- fit
  }

  list(
    selection = sel_fit,
    up_map = up_map,
    upstream_if = up_all,
    tau_hat = tau_hat,
    eta_hat = eta_hat,
    tau_if = tau_if,
    eta_if = eta_if,
    shared = state$shared,
    regressions = state$regressions,
    treatments = treatments,
    links = state$links,
    estimates = vapply(treatments, function(f) unname(f$means["mu_dr"]), numeric(1)),
    estimates_or = vapply(state$links, `[[`, numeric(1), "mu"),
    estimates_ipw = vapply(treatments, function(f) unname(f$means["mu_ipw"]), numeric(1)),
    phi = phi,
    phi_core = phi * (n_A / N),
    psi_dir = psi_dir,
    psi_or = psi_or,
    psi_ps = psi_ps,
    psi_sel = psi_sel,
    psi_ipw = psi_ipw,
    psi_or_estimator = core$phi,
    beta_if = core$beta,
    alpha_if = core$alpha
  )
}
