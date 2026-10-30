# Independent dense-small point and central finite-difference IF reference for
# D48 Step 2 Section B (IPW and shared outcome-change OR consumers across all 3
# independent_period_logit modes on F and all 9 pair_logit modes on F/B/A).
#
# Uses base R (stats::glm.fit and stats::lm.wfit) only; zero RNG, zero package code.

sp_ref_tol <- 3e-7
sp_ref_root_tol <- 2e-11
sp_ref_steps <- c(1e-4, 5e-5)

# 9 base observation patterns per half-cohort (6 repeat observers in A, 2 singletons
# in B \ A, 1 zero-observer in F \ B) across 4 calendar dates t in (1, 2, 3, 4),
# duplicated across j in (-1, +1) for 18 rows per cohort (N = 54 total).
sp_ref_patterns9 <- rbind(
  c(1L, 1L, 1L, 1L),
  c(1L, 1L, 1L, 0L),
  c(1L, 1L, 0L, 1L),
  c(1L, 0L, 1L, 1L),
  c(0L, 1L, 1L, 1L),
  c(1L, 1L, 1L, 1L),
  c(1L, 0L, 0L, 0L),
  c(0L, 1L, 0L, 0L),
  c(0L, 0L, 0L, 0L)
)
sp_ref_patterns <- rbind(sp_ref_patterns9, sp_ref_patterns9)

sp_ref_cohorts <- c(3, 5, Inf)
sp_ref_calendar <- c(1, 2, 3, 4)

sp_ref_build_fixture <- function() {
  rows <- list()
  idx <- 0L
  for (k in seq_along(sp_ref_cohorts)) {
    for (j in c(-1L, 1L)) {
      for (h in 1:9) {
        idx <- idx + 1L
        hh <- if (j == -1L) h else h + 9L
        z1 <- 0.7 * j + 0.25 * ((h %% 3L) - 1.0) + 0.08 * (k - 2.0)
        z2 <- 0.5 * ((h %% 5L) - 2.0) - 0.2 * j * ((h %% 2L) - 0.5) - 0.06 * (k - 2.0)
        x1 <- 0.65 * j - 0.2 * ((h %% 3L) - 1.0) + 0.10 * (k - 2.0)
        x2 <- 0.45 * ((h %% 4L) - 1.5) + 0.25 * j * ((h %% 2L) - 0.5) - 0.08 * (k - 2.0)
        w_i <- 0.7 + ((7L * idx) %% 11L) / 10.0
        cl_i <- 1L + ((idx - 1L) %% 9L)
        rows[[idx]] <- data.frame(
          id = idx,
          k = k,
          g = sp_ref_cohorts[k],
          h = hh,
          z1 = z1,
          z2 = z2,
          x1 = x1,
          x2 = x2,
          w = w_i,
          cluster = cl_i,
          frame_only = (hh == 18L),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dat <- do.call(rbind, rows)
  N <- nrow(dat)
  obs <- sp_ref_patterns[dat$h, , drop = FALSE]
  n_obs <- rowSums(obs)
  A <- as.double(n_obs >= 2L)
  B <- as.double(n_obs >= 1L)

  # Potential outcomes with heterogeneous level and trend; masked outside observed dates
  Y_full <- sapply(seq_along(sp_ref_calendar), function(t_idx) {
    tt <- sp_ref_calendar[t_idx]
    0.35 * dat$k +
      tt * (0.45 * dat$x1 - 0.30 * dat$x2 + 0.22 * dat$k) +
      (tt^2) * (0.08 * dat$x1^2 + 0.06 * dat$x2^2 + 0.05 * dat$k) +
      ((tt + dat$h) %% 2L) * 0.14 * dat$z1
  })
  Y <- Y_full
  Y[obs == 0L] <- NA_real_

  # Fixed known normalized transport ratios rho_U on B and A within each cohort
  rho_F <- rep.int(1.0, N)
  rho_B <- rep.int(NA_real_, N)
  rho_A <- rep.int(NA_real_, N)
  raw_tilt <- exp(0.18 * dat$z1 - 0.12 * dat$z2 + 0.07 * dat$k)
  for (k in seq_along(sp_ref_cohorts)) {
    ib <- which(dat$k == k & B == 1)
    ia <- which(dat$k == k & A == 1)
    # Normalize so sum_{i in U, G=k} w_i / rho_i = sum_{i in U, G=k} w_i
    c_b <- sum(dat$w[ib] / raw_tilt[ib]) / sum(dat$w[ib])
    c_a <- sum(dat$w[ia] / raw_tilt[ia]) / sum(dat$w[ia])
    rho_B[ib] <- raw_tilt[ib] * c_b
    rho_A[ia] <- raw_tilt[ia] * c_a
  }

  Z <- cbind(`(Intercept)` = 1.0, z1 = dat$z1, z2 = dat$z2)
  attr(Z, "intercept") <- TRUE
  attr(Z, "assign") <- c(0L, 1L, 2L)

  V <- cbind(`(Intercept)` = 1.0, x1 = dat$x1, x2 = dat$x2)
  attr(V, "intercept") <- TRUE
  attr(V, "assign") <- c(0L, 1L, 2L)

  links <- list(
    list(link_id = "NYT12", s = 1, t = 2, g = 3, K = c(5, Inf)),
    list(link_id = "NYT24", s = 2, t = 4, g = 3, K = c(5, Inf)),
    list(link_id = "NT12",  s = 1, t = 2, g = 3, K = c(Inf)),
    list(link_id = "NT13",  s = 1, t = 3, g = 5, K = c(Inf))
  )

  # Long panel data (excluding frame_only rows; h == 16 is retained in panel as all-NA unit)
  panel_ids <- which(!dat$frame_only)
  long <- data.frame(
    id = rep(dat$id[panel_ids], each = length(sp_ref_calendar)),
    time = rep(sp_ref_calendar, length(panel_ids)),
    y = as.vector(t(Y[panel_ids, , drop = FALSE])),
    g = rep(ifelse(is.infinite(dat$g[panel_ids]), 0, dat$g[panel_ids]), each = length(sp_ref_calendar)),
    x1 = rep(dat$x1[panel_ids], each = length(sp_ref_calendar)),
    x2 = rep(dat$x2[panel_ids], each = length(sp_ref_calendar)),
    z1 = rep(dat$z1[panel_ids], each = length(sp_ref_calendar)),
    z2 = rep(dat$z2[panel_ids], each = length(sp_ref_calendar)),
    w = rep(dat$w[panel_ids], each = length(sp_ref_calendar)),
    cluster = rep(dat$cluster[panel_ids], each = length(sp_ref_calendar)),
    stringsAsFactors = FALSE
  )
  frame_df <- data.frame(
    id = dat$id,
    g = ifelse(is.infinite(dat$g), 0, dat$g),
    x1 = dat$x1,
    x2 = dat$x2,
    z1 = dat$z1,
    z2 = dat$z2,
    w = dat$w,
    cluster = dat$cluster,
    rho_B = rho_B,
    rho_A = rho_A,
    stringsAsFactors = FALSE
  )

  list(
    N = N,
    dat = dat,
    obs = obs,
    A = A,
    B = B,
    Y = Y,
    Z = Z,
    V = V,
    rho = list(F = rho_F, B = rho_B, A = rho_A),
    U = list(F = rep.int(TRUE, N), B = (B == 1), A = (A == 1)),
    links = links,
    long = long,
    frame_df = frame_df
  )
}

sp_ref_fit_logit <- function(X, y, ww) {
  f <- withCallingHandlers(
    stats::glm.fit(
      x = X,
      y = y,
      weights = ww,
      family = stats::binomial(),
      control = stats::glm.control(epsilon = 1e-13, maxit = 60L)
    ),
    warning = function(w) {
      if (identical(conditionMessage(w), "non-integer #successes in a binomial glm!")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  stopifnot(f$converged, f$rank == ncol(X), all(is.finite(f$coefficients)))
  pr <- stats::plogis(drop(X %*% f$coefficients))
  stopifnot(all(pr > 0 & pr < 1), max(abs(crossprod(X, ww * (y - pr)))) < sp_ref_root_tol)
  unname(as.double(f$coefficients))
}

# Independent dense selection fitter using explicit cell-indicator dummy matrices
# and base stats::glm.fit.
sp_ref_fit_selection <- function(fx, cw = rep.int(1.0, fx$N),
                                 selection_model = "pair_logit",
                                 selection_pooling = "separate",
                                 universe = "F",
                                 has_intercept = TRUE,
                                 zero_col = FALSE,
                                 slope_cols = c("z1", "z2")) {
  N <- fx$N
  ww <- fx$dat$w * cw
  U <- fx$U[[universe]]
  rho <- fx$rho[[universe]]
  A <- fx$A
  K <- length(sp_ref_cohorts)

  mass_U <- numeric(K)
  repeat_mass <- numeric(K)
  for (k in seq_len(K)) {
    ik_u <- which(fx$dat$k == k & U)
    ik_a <- which(fx$dat$k == k & A == 1)
    mass_U[k] <- sum(ww[ik_u])
    repeat_mass[k] <- sum(ww[ik_a])
  }
  a <- repeat_mass / mass_U

  # Collect active pair-cohorts from fx$links in link traversal order
  pair_cohorts <- list()
  for (ell in seq_along(fx$links)) {
    lk <- fx$links[[ell]]
    ks <- sort(lk$K)
    cohs <- match(c(lk$g, ks), sp_ref_cohorts)
    s_idx <- match(lk$s, sp_ref_calendar)
    t_idx <- match(lk$t, sp_ref_calendar)
    for (k in cohs) {
      pc <- list(cohort = k, s = lk$s, t = lk$t, s_idx = s_idx, t_idx = t_idx)
      same <- which(vapply(pair_cohorts, identical, logical(1), pc))
      if (!length(same)) pair_cohorts[[length(pair_cohorts) + 1L]] <- pc
    }
  }

  # Build active cells in canonical order
  cells <- list()
  if (identical(selection_model, "independent_period_logit")) {
    for (k in seq_len(K)) {
      for (t_idx in seq_along(sp_ref_calendar)) {
        used <- any(vapply(pair_cohorts, function(pc) {
          pc$cohort == k && (pc$s_idx == t_idx || pc$t_idx == t_idx)
        }, logical(1)))
        if (used) {
          cells[[length(cells) + 1L]] <- list(kind = "date", cohort = k, t_idx = t_idx)
        }
      }
    }
  } else {
    u_pairs <- list()
    for (pc in pair_cohorts) {
      pr <- list(s_idx = pc$s_idx, t_idx = pc$t_idx)
      if (!any(vapply(u_pairs, identical, logical(1), pr))) {
        u_pairs[[length(u_pairs) + 1L]] <- pr
      }
    }
    for (pr in u_pairs) {
      for (k in seq_len(K)) {
        used <- any(vapply(pair_cohorts, function(pc) {
          pc$cohort == k && pc$s_idx == pr$s_idx && pc$t_idx == pr$t_idx
        }, logical(1)))
        if (used) {
          cells[[length(cells) + 1L]] <- list(kind = "pair", cohort = k, s_idx = pr$s_idx, t_idx = pr$t_idx)
        }
      }
    }
  }

  C <- length(cells)
  block_groups <- if (identical(selection_pooling, "separate")) {
    lapply(seq_len(C), function(cid) cid)
  } else if (identical(selection_pooling, "within_cohort")) {
    cohs_active <- sort(unique(vapply(cells, `[[`, integer(1), "cohort")))
    lapply(cohs_active, function(k) which(vapply(cells, `[[`, integer(1), "cohort") == k))
  } else {
    list(seq_len(C))
  }

  Z_slope <- if (zero_col || !length(slope_cols)) matrix(0, N, 0L) else fx$Z[, slope_cols, drop = FALSE]
  d_s <- ncol(Z_slope)

  cell_prob <- matrix(0, N, C)
  alpha_list <- list()
  beta_list <- list()

  for (m in seq_along(block_groups)) {
    cids <- block_groups[[m]]
    B_m <- length(cids)
    if (zero_col) {
      for (cid in cids) cell_prob[, cid] <- 0.5
      next
    }
    X_rows <- list()
    y_vec <- numeric(0)
    w_vec <- numeric(0)
    for (b in seq_len(B_m)) {
      cid <- cids[b]
      cell <- cells[[cid]]
      ik <- which(fx$dat$k == cell$cohort & U)
      y_cell <- if (identical(cell$kind, "date")) {
        fx$obs[ik, cell$t_idx]
      } else {
        as.double(fx$obs[ik, cell$s_idx] == 1L & fx$obs[ik, cell$t_idx] == 1L)
      }
      if (has_intercept) {
        E_b <- matrix(0, length(ik), B_m)
        E_b[, b] <- 1.0
        X_rows[[b]] <- cbind(E_b, Z_slope[ik, , drop = FALSE])
      } else {
        X_rows[[b]] <- Z_slope[ik, , drop = FALSE]
      }
      y_vec <- c(y_vec, y_cell)
      w_vec <- c(w_vec, ww[ik])
    }
    X_blk <- do.call(rbind, X_rows)
    coef_blk <- sp_ref_fit_logit(X_blk, y_vec, w_vec)
    al_m <- if (has_intercept) coef_blk[seq_len(B_m)] else numeric(0)
    be_m <- if (d_s > 0L) coef_blk[(if (has_intercept) B_m else 0L) + seq_len(d_s)] else numeric(0)
    alpha_list[[m]] <- al_m
    beta_list[[m]] <- be_m
    for (b in seq_len(B_m)) {
      cid <- cids[b]
      eta <- (if (has_intercept) al_m[b] else 0.0) +
        (if (d_s > 0L) drop(Z_slope %*% be_m) else 0.0)
      cell_prob[, cid] <- stats::plogis(eta)
    }
  }

  # Pair probabilities and fitted selection weights for each pair_cohort
  n_pc <- length(pair_cohorts)
  pair_prob <- matrix(0, N, n_pc)
  for (q in seq_len(n_pc)) {
    pc <- pair_cohorts[[q]]
    if (identical(selection_model, "independent_period_logit")) {
      cid_s <- which(vapply(cells, function(c) c$cohort == pc$cohort && c$t_idx == pc$s_idx, logical(1)))
      cid_t <- which(vapply(cells, function(c) c$cohort == pc$cohort && c$t_idx == pc$t_idx, logical(1)))
      pair_prob[, q] <- cell_prob[, cid_s] * cell_prob[, cid_t]
    } else {
      cid_p <- which(vapply(cells, function(c) {
        c$cohort == pc$cohort && c$s_idx == pc$s_idx && c$t_idx == pc$t_idx
      }, logical(1)))
      pair_prob[, q] <- cell_prob[, cid_p]
    }
  }

  list(
    a = a,
    cells = cells,
    block_groups = block_groups,
    alpha = unlist(alpha_list, use.names = FALSE),
    beta = unlist(beta_list, use.names = FALSE),
    pair_cohorts = pair_cohorts,
    pair_prob = pair_prob,
    ww = ww,
    rho = rho
  )
}

sp_ref_ipw_point <- function(fx, cw = rep.int(1.0, fx$N),
                             selection_model = "pair_logit",
                             selection_pooling = "separate",
                             universe = "F",
                             has_intercept = TRUE,
                             zero_col = FALSE,
                             slope_cols = c("z1", "z2")) {
  sel <- sp_ref_fit_selection(
    fx, cw = cw, selection_model = selection_model,
    selection_pooling = selection_pooling, universe = universe,
    has_intercept = has_intercept, zero_col = zero_col,
    slope_cols = slope_cols
  )
  nc <- sum(cw)
  L <- length(fx$links)
  trt_beta <- vector("list", L)
  T_vec <- M_vec <- dp_vec <- dm_vec <- mu_vec <- numeric(L)

  for (ell in seq_len(L)) {
    lk <- fx$links[[ell]]
    ks <- sort(lk$K)
    cohs <- match(c(lk$g, ks), sp_ref_cohorts)
    s_idx <- match(lk$s, sp_ref_calendar)
    t_idx <- match(lk$t, sp_ref_calendar)
    S_pair <- as.double(fx$obs[, s_idx] == 1L & fx$obs[, t_idx] == 1L)
    obs_rows <- which(fx$dat$k %in% cohs & S_pair == 1)
    ci_obs <- match(fx$dat$k[obs_rows], cohs)
    D <- as.double(ci_obs == 1L)
    dy <- fx$Y[obs_rows, t_idx] - fx$Y[obs_rows, s_idx]
    fw <- numeric(length(obs_rows))
    for (k_pos in seq_along(cohs)) {
      jj <- which(ci_obs == k_pos)
      k_val <- cohs[k_pos]
      q <- which(vapply(sel$pair_cohorts, function(pc) {
        pc$cohort == k_val && pc$s_idx == s_idx && pc$t_idx == t_idx
      }, logical(1)))
      rr <- obs_rows[jj]
      fw[jj] <- sel$ww[rr] * sel$a[k_val] / (sel$rho[rr] * sel$pair_prob[rr, q])
    }
    X <- fx$V[obs_rows, , drop = FALSE]
    b_coef <- sp_ref_fit_logit(X, D, fw)
    odds <- exp(drop(X %*% b_coef))
    wp <- fw * D
    wm <- fw * (1 - D) * odds
    dp <- sum(wp) / nc
    dm <- sum(wm) / nc
    T_m <- sum(wp * dy) / sum(wp)
    M_m <- sum(wm * dy) / sum(wm)
    trt_beta[[ell]] <- b_coef
    dp_vec[ell] <- dp
    dm_vec[ell] <- dm
    T_vec[ell] <- T_m
    M_vec[ell] <- M_m
    mu_vec[ell] <- T_m - M_m
  }

  vec <- c(
    a = sel$a,
    sel_alpha = sel$alpha,
    sel_beta = sel$beta,
    trt_beta = unlist(trt_beta, use.names = FALSE),
    T = T_vec,
    M = M_vec,
    dp = dp_vec,
    dm = dm_vec,
    mu = mu_vec
  )
  list(
    a = sel$a,
    sel_cells = sel$cells,
    sel_block_groups = sel$block_groups,
    sel_alpha = sel$alpha,
    sel_beta = sel$beta,
    trt_beta = trt_beta,
    T = T_vec,
    M = M_vec,
    dp = dp_vec,
    dm = dm_vec,
    mu = mu_vec,
    vector = vec
  )
}

sp_ref_or_point <- function(fx, cw = rep.int(1.0, fx$N),
                            selection_model = "pair_logit",
                            selection_pooling = "separate",
                            universe = "F",
                            has_intercept = TRUE,
                            zero_col = FALSE,
                            slope_cols = c("z1", "z2"),
                            or_intercept = TRUE) {
  sel <- sp_ref_fit_selection(
    fx, cw = cw, selection_model = selection_model,
    selection_pooling = selection_pooling, universe = universe,
    has_intercept = has_intercept, zero_col = zero_col,
    slope_cols = slope_cols
  )
  L <- length(fx$links)
  # Unique control blocks in link traversal order: (pair, sorted K)
  blocks <- list()
  link_block <- integer(L)
  for (ell in seq_len(L)) {
    lk <- fx$links[[ell]]
    s_idx <- match(lk$s, sp_ref_calendar)
    t_idx <- match(lk$t, sp_ref_calendar)
    ctrl_cohs <- match(sort(lk$K), sp_ref_cohorts)
    bspec <- list(s_idx = s_idx, t_idx = t_idx, controls = ctrl_cohs, L_step = t_idx - s_idx)
    same <- which(vapply(blocks, identical, logical(1), bspec))
    if (length(same)) {
      link_block[ell] <- same[1L]
    } else {
      blocks[[length(blocks) + 1L]] <- bspec
      link_block[ell] <- length(blocks)
    }
  }

  G <- length(blocks)
  X_slope <- fx$V[, c("x1", "x2"), drop = FALSE]
  p_or <- ncol(X_slope)
  XX_list <- list()
  YY_vec <- numeric(0)
  WW_vec <- numeric(0)

  for (b in seq_len(G)) {
    bl <- blocks[[b]]
    s_idx <- bl$s_idx
    t_idx <- bl$t_idx
    S_pair <- as.double(fx$obs[, s_idx] == 1L & fx$obs[, t_idx] == 1L)
    for (k_val in bl$controls) {
      rows <- which(fx$dat$k == k_val & S_pair == 1)
      q <- which(vapply(sel$pair_cohorts, function(pc) {
        pc$cohort == k_val && pc$s_idx == s_idx && pc$t_idx == t_idx
      }, logical(1)))
      q_w <- sel$ww[rows] * sel$a[k_val] / (sel$rho[rows] * sel$pair_prob[rows, q])
      dy <- fx$Y[rows, t_idx] - fx$Y[rows, s_idx]
      d_mat <- bl$L_step * X_slope[rows, , drop = FALSE]
      if (or_intercept) {
        E_b <- matrix(0, length(rows), G)
        E_b[, b] <- 1.0
        XX_list[[length(XX_list) + 1L]] <- cbind(E_b, d_mat)
      } else {
        XX_list[[length(XX_list) + 1L]] <- d_mat
      }
      YY_vec <- c(YY_vec, dy)
      WW_vec <- c(WW_vec, q_w)
    }
  }

  XX <- do.call(rbind, XX_list)
  wfit <- stats::lm.wfit(XX, YY_vec, WW_vec, tol = 1e-12)
  coef <- unname(as.double(wfit$coefficients))
  or_alpha <- if (or_intercept) coef[seq_len(G)] else numeric(G)
  or_beta <- coef[(if (or_intercept) G else 0L) + seq_len(p_or)]

  mu_vec <- numeric(L)
  for (ell in seq_len(L)) {
    lk <- fx$links[[ell]]
    b <- link_block[ell]
    bl <- blocks[[b]]
    s_idx <- bl$s_idx
    t_idx <- bl$t_idx
    g_idx <- match(lk$g, sp_ref_cohorts)
    S_pair <- as.double(fx$obs[, s_idx] == 1L & fx$obs[, t_idx] == 1L)
    rows <- which(fx$dat$k == g_idx & S_pair == 1)
    q <- which(vapply(sel$pair_cohorts, function(pc) {
      pc$cohort == g_idx && pc$s_idx == s_idx && pc$t_idx == t_idx
    }, logical(1)))
    T_w <- sel$ww[rows] * sel$a[g_idx] / (sel$rho[rows] * sel$pair_prob[rows, q])
    dy <- fx$Y[rows, t_idx] - fx$Y[rows, s_idx]
    res <- dy - or_alpha[b] - drop(bl$L_step * X_slope[rows, , drop = FALSE] %*% or_beta)
    mu_vec[ell] <- sum(T_w * res) / sum(T_w)
  }

  vec <- c(
    a = sel$a,
    sel_alpha = sel$alpha,
    sel_beta = sel$beta,
    or_alpha = if (or_intercept) or_alpha else numeric(0),
    or_beta = or_beta,
    mu = mu_vec
  )
  list(
    a = sel$a,
    sel_cells = sel$cells,
    sel_block_groups = sel$block_groups,
    sel_alpha = sel$alpha,
    sel_beta = sel$beta,
    or_block_first_link = match(seq_len(G), link_block),
    or_alpha = or_alpha,
    or_beta = or_beta,
    mu = mu_vec,
    vector = vec
  )
}

sp_ref_central_fd <- function(N, eval_fn, steps = sp_ref_steps) {
  base <- eval_fn(rep.int(1.0, N))
  p <- length(base$vector)
  L <- length(base$mu)
  fds <- vector("list", length(steps))
  for (s_i in seq_along(steps)) {
    h <- steps[s_i]
    fd <- matrix(0, N, p)
    for (i in seq_len(N)) {
      cp <- cm <- rep.int(1.0, N)
      cp[i] <- 1.0 + h
      cm[i] <- 1.0 - h
      fd[i, ] <- N * (eval_fn(cp)$vector - eval_fn(cm)$vector) / (2.0 * h)
    }
    fds[[s_i]] <- fd
  }
  step_diff <- if (length(steps) >= 2L) max(abs(fds[[1L]] - fds[[2L]])) else 0.0
  stopifnot(step_diff < sp_ref_tol)
  fd_last <- fds[[length(steps)]]
  mu_cols <- seq.int(p - L + 1L, p)
  list(
    point = base,
    fd_stack = fd_last,
    fd_phi = fd_last[, mu_cols, drop = FALSE],
    step_diff = step_diff
  )
}

# Public catalog and aggregation reference for both est_method = "ipw" and
# est_method = "reg" across all 3 independent_period_logit modes on F and all 9
# pair_logit modes on F/B/A. Reuses cp06_pooled_public_catalog and
# cp06_pooled_public_cells_from_mu from cp06_pooled_public_or_reference.R.
sp_ref_public_eval <- function(raw, w = raw$w, est_method = "ipw",
                               selection_model = "pair_logit",
                               selection_pooling = "separate",
                               route = "F",
                               sel_intercept = TRUE,
                               sel_zero_col = FALSE,
                               Q_slope = cbind(x = raw$x),
                               or_intercept = TRUE,
                               X_slope = cbind(x = raw$x),
                               V_treat = cbind(`(Intercept)` = 1, x = raw$x),
                               dates = 1:4,
                               control_group = "nevertreated",
                               anticipation = 0L,
                               pt_assumption = "post",
                               pt_start = NULL,
                               base_period = "varying",
                               cat_info = NULL) {
  if (is.null(cat_info)) {
    cat_info <- cp06_pooled_public_catalog(
      raw, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption,
      pt_start = pt_start, base_period = base_period
    )
  }
  cal <- cat_info$calendar
  links_df <- cat_info$links_df
  N <- nrow(raw)
  idx_time <- function(tt) match(tt, cal)
  R_mat <- !is.na(as.matrix(raw[paste0("Y", seq_along(cal))]))
  A_vec <- as.numeric(rowSums(R_mat) >= 2L)
  B_vec <- as.numeric(rowSums(R_mat) >= 1L)
  U_vec <- if (route == "F") rep(1, N) else if (route == "B") B_vec else A_vec
  rho_vec <- if (route %in% c("B", "A")) raw$rho else rep(1, N)

  cohorts_all <- cat_info$cohorts_all
  K_len <- length(cohorts_all)
  mass_F <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k]), numeric(1))
  mass_U <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & U_vec == 1]), numeric(1))
  mass_A <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & A_vec == 1]), numeric(1))
  tau <- mass_F / sum(mass_F)
  eta <- mass_A / sum(mass_A)
  a_U <- ifelse(mass_U > 0, mass_A / mass_U, 0)
  names(tau) <- names(eta) <- names(a_U) <- as.character(cohorts_all)

  pc_list <- list()
  for (r in seq_len(nrow(links_df))) {
    is <- links_df$is[r]; it <- links_df$it[r]
    ks <- c(links_df$g[r], sort(as.numeric(strsplit(links_df$controls[r], "|", fixed = TRUE)[[1L]])))
    for (kv in ks) {
      k_idx <- match(kv, cohorts_all)
      pc <- list(k_idx = k_idx, cohort = kv, is = is, it = it)
      if (!any(vapply(pc_list, identical, logical(1), pc))) {
        pc_list[[length(pc_list) + 1L]] <- pc
      }
    }
  }
  ord_pc <- order(
    vapply(pc_list, `[[`, integer(1), "k_idx"),
    vapply(pc_list, `[[`, integer(1), "is"),
    vapply(pc_list, `[[`, integer(1), "it")
  )
  pc_list <- pc_list[ord_pc]

  cells <- list()
  if (identical(selection_model, "independent_period_logit")) {
    for (k_idx in seq_len(K_len)) {
      for (t_idx in seq_along(cal)) {
        used <- any(vapply(pc_list, function(pc) {
          pc$k_idx == k_idx && (pc$is == t_idx || pc$it == t_idx)
        }, logical(1)))
        if (used) {
          cells[[length(cells) + 1L]] <- list(
            kind = "date", k_idx = k_idx, cohort = cohorts_all[k_idx], t_idx = t_idx
          )
        }
      }
    }
  } else {
    up_pairs <- unique(data.frame(
      is = vapply(pc_list, `[[`, integer(1), "is"),
      it = vapply(pc_list, `[[`, integer(1), "it")
    ))
    up_pairs <- up_pairs[order(up_pairs$is, up_pairs$it), , drop = FALSE]
    for (p_i in seq_len(nrow(up_pairs))) {
      is <- up_pairs$is[p_i]; it <- up_pairs$it[p_i]
      for (k_idx in seq_len(K_len)) {
        used <- any(vapply(pc_list, function(pc) {
          pc$k_idx == k_idx && pc$is == is && pc$it == it
        }, logical(1)))
        if (used) {
          cells[[length(cells) + 1L]] <- list(
            kind = "pair", k_idx = k_idx, cohort = cohorts_all[k_idx], is = is, it = it
          )
        }
      }
    }
  }
  C_cells <- length(cells)
  block_groups <- if (identical(selection_pooling, "separate")) {
    lapply(seq_len(C_cells), function(cid) cid)
  } else if (identical(selection_pooling, "within_cohort")) {
    active_k <- sort(unique(vapply(cells, `[[`, integer(1), "k_idx")))
    lapply(active_k, function(k_idx) which(vapply(cells, `[[`, integer(1), "k_idx") == k_idx))
  } else {
    list(seq_len(C_cells))
  }

  Z_sl <- if (sel_zero_col || is.null(Q_slope)) matrix(0, N, 0L) else as.matrix(Q_slope)
  d_s <- ncol(Z_sl)
  cell_prob <- matrix(0, N, C_cells)
  sel_alpha <- list()
  sel_beta <- list()
  for (m in seq_along(block_groups)) {
    cids <- block_groups[[m]]
    B_m <- length(cids)
    if (sel_zero_col && !sel_intercept) {
      for (cid in cids) cell_prob[, cid] <- 0.5
      next
    }
    X_rows <- list(); y_vec <- numeric(0); w_vec <- numeric(0)
    for (b in seq_len(B_m)) {
      cid <- cids[b]; cell <- cells[[cid]]
      ik <- which(raw$cohort == cell$cohort & U_vec == 1)
      y_c <- if (identical(cell$kind, "date")) {
        as.double(R_mat[ik, cell$t_idx])
      } else {
        as.double(R_mat[ik, cell$is] & R_mat[ik, cell$it])
      }
      if (sel_intercept) {
        E_b <- matrix(0, length(ik), B_m); E_b[, b] <- 1.0
        X_rows[[b]] <- cbind(E_b, Z_sl[ik, , drop = FALSE])
      } else {
        X_rows[[b]] <- Z_sl[ik, , drop = FALSE]
      }
      y_vec <- c(y_vec, y_c); w_vec <- c(w_vec, w[ik])
    }
    coef_blk <- sp_ref_fit_logit(do.call(rbind, X_rows), y_vec, w_vec)
    al_m <- if (sel_intercept) coef_blk[seq_len(B_m)] else numeric(0)
    be_m <- if (d_s > 0L) coef_blk[(if (sel_intercept) B_m else 0L) + seq_len(d_s)] else numeric(0)
    sel_alpha[[m]] <- al_m
    sel_beta[[m]] <- be_m
    for (b in seq_len(B_m)) {
      cid <- cids[b]
      eta_lp <- (if (sel_intercept) al_m[b] else 0.0) + (if (d_s > 0L) drop(Z_sl %*% be_m) else 0.0)
      cell_prob[, cid] <- stats::plogis(eta_lp)
    }
  }

  get_pair_prob <- function(is, it, kv) {
    k_idx <- match(kv, cohorts_all)
    if (identical(selection_model, "independent_period_logit")) {
      cid_s <- which(vapply(cells, function(c) c$k_idx == k_idx && c$t_idx == is, logical(1)))
      cid_t <- which(vapply(cells, function(c) c$k_idx == k_idx && c$t_idx == it, logical(1)))
      cell_prob[, cid_s] * cell_prob[, cid_t]
    } else {
      cid_p <- which(vapply(cells, function(c) c$k_idx == k_idx && c$is == is && c$it == it, logical(1)))
      cell_prob[, cid_p]
    }
  }

  alpha <- numeric(0)
  beta <- numeric(0)
  if (identical(est_method, "ipw")) {
    mu <- numeric(nrow(links_df))
    for (l in seq_len(nrow(links_df))) {
      is <- links_df$is[l]; it <- links_df$it[l]
      g_val <- links_df$g[l]
      ctrl_cohs <- sort(as.numeric(strsplit(links_df$controls[l], "|", fixed = TRUE)[[1L]]))
      cohs <- c(g_val, ctrl_cohs)
      S_pair <- as.numeric(R_mat[, is] & R_mat[, it])
      obs_rows <- which(raw$cohort %in% cohs & S_pair == 1)
      D <- as.double(raw$cohort[obs_rows] == g_val)
      dy <- raw[[paste0("Y", it)]][obs_rows] - raw[[paste0("Y", is)]][obs_rows]
      fw <- numeric(length(obs_rows))
      for (kv in cohs) {
        jj <- which(raw$cohort[obs_rows] == kv)
        rr <- obs_rows[jj]
        pp <- get_pair_prob(is, it, kv)[rr]
        fw[jj] <- w[rr] * a_U[as.character(kv)] / (rho_vec[rr] * pp)
      }
      X_tr <- V_treat[obs_rows, , drop = FALSE]
      b_tr <- sp_ref_fit_logit(X_tr, D, fw)
      odds <- exp(drop(X_tr %*% b_tr))
      wp <- fw * D; wm <- fw * (1 - D) * odds
      mu[l] <- sum(wp * dy) / sum(wp) - sum(wm * dy) / sum(wm)
    }
  } else {
    upairs <- unique(data.frame(
      s = links_df$s, t = links_df$t,
      is = idx_time(links_df$s), it = idx_time(links_df$t),
      stringsAsFactors = FALSE
    ))
    upairs <- upairs[order(upairs$is, upairs$it), , drop = FALSE]
    upairs$pair_id <- seq_len(nrow(upairs)); upairs$L <- upairs$it - upairs$is
    link_pair_id <- match(paste0(links_df$s, "_", links_df$t), paste0(upairs$s, "_", upairs$t))
    block_df <- unique(data.frame(pair_id = link_pair_id, controls = links_df$controls, stringsAsFactors = FALSE))
    pool_sort_key <- vapply(block_df$controls, function(cs) {
      ks <- match(as.numeric(strsplit(cs, "|", fixed = TRUE)[[1L]]), cohorts_all)
      paste(sprintf("%012d", ks), collapse = "|")
    }, character(1))
    block_df <- block_df[order(sprintf("%012d|%s", block_df$pair_id, pool_sort_key), method = "radix"), , drop = FALSE]
    G <- nrow(block_df)
    link_block_id <- match(
      paste0(link_pair_id, "|", links_df$controls),
      paste0(block_df$pair_id, "|", block_df$controls)
    )
    p_or <- if (is.null(X_slope)) 0L else ncol(X_slope)
    XX <- NULL; YY <- numeric(0); WW <- numeric(0)
    for (b in seq_len(G)) {
      pid <- block_df$pair_id[b]
      is <- upairs$is[pid]; it <- upairs$it[pid]; L_b <- upairs$L[pid]
      ctrl_cohs <- as.numeric(strsplit(block_df$controls[b], "|", fixed = TRUE)[[1L]])
      for (kv in ctrl_cohs) {
        rows <- which(raw$cohort == kv & R_mat[, is] & R_mat[, it])
        pp <- get_pair_prob(is, it, kv)[rows]
        qw <- w[rows] * a_U[as.character(kv)] / (rho_vec[rows] * pp)
        E <- if (or_intercept) {
          m <- matrix(0, length(rows), G); m[, b] <- 1; m
        } else matrix(0, length(rows), 0L)
        D_m <- if (p_or > 0L) L_b * X_slope[rows, , drop = FALSE] else matrix(0, length(rows), 0L)
        XX <- rbind(XX, cbind(E, D_m))
        YY <- c(YY, raw[[paste0("Y", it)]][rows] - raw[[paste0("Y", is)]][rows])
        WW <- c(WW, qw)
      }
    }
    wfit <- stats::lm.wfit(x = XX, y = YY, w = WW, tol = 1e-12)
    coef <- unname(wfit$coefficients)
    alpha <- if (or_intercept) coef[seq_len(G)] else numeric(G)
    beta <- if (p_or > 0L) coef[(if (or_intercept) G else 0L) + seq_len(p_or)] else numeric(0)
    if (p_or > 0L && !is.null(colnames(X_slope))) names(beta) <- colnames(X_slope)
    mu <- vapply(seq_len(nrow(links_df)), function(l) {
      b <- link_block_id[l]; pid <- link_pair_id[l]
      is <- upairs$is[pid]; it <- upairs$it[pid]; L_l <- upairs$L[pid]
      kv <- links_df$g[l]
      rows <- which(raw$cohort == kv & R_mat[, is] & R_mat[, it])
      pp <- get_pair_prob(is, it, kv)[rows]
      qw <- w[rows] * a_U[as.character(kv)] / (rho_vec[rows] * pp)
      dy <- raw[[paste0("Y", it)]][rows] - raw[[paste0("Y", is)]][rows]
      pred <- alpha[b] + if (p_or > 0L) as.vector(L_l * X_slope[rows, , drop = FALSE] %*% beta) else 0
      sum(qw * (dy - pred)) / sum(qw)
    }, numeric(1))
  }

  cells_out <- cp06_pooled_public_cells_from_mu(
    cat_info, mu, anticipation = anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start, base_period = base_period
  )
  post <- which(cells_out$time >= cells_out$g)
  agg_overall <- c(simple = NA_real_, group = NA_real_, calendar = NA_real_, dynamic = NA_real_)
  if (length(post) && all(!is.na(cells_out$att[post]))) {
    g_post <- as.character(cells_out$g[post])
    agg_overall["simple"] <- sum(tau[g_post] * cells_out$att[post]) / sum(tau[g_post])
    ug <- sort(unique(cells_out$g[post]))
    est_g <- vapply(ug, function(gv) mean(cells_out$att[post[cells_out$g[post] == gv]]), numeric(1))
    agg_overall["group"] <- sum(tau[as.character(ug)] * est_g) / sum(tau[as.character(ug)])
    ut <- sort(unique(cells_out$time[post]))
    est_t <- vapply(ut, function(tv) {
      rt <- post[cells_out$time[post] == tv]
      gt <- as.character(cells_out$g[rt])
      sum(tau[gt] * cells_out$att[rt]) / sum(tau[gt])
    }, numeric(1))
    agg_overall["calendar"] <- mean(est_t)
    ue <- sort(unique(cells_out$time[post] - cells_out$g[post]))
    est_e <- vapply(ue, function(ev) {
      re <- post[(cells_out$time[post] - cells_out$g[post]) == ev]
      ge <- as.character(cells_out$g[re])
      sum(tau[ge] * cells_out$att[re]) / sum(tau[ge])
    }, numeric(1))
    agg_overall["dynamic"] <- mean(est_e)
  }

  list(
    alpha = alpha, beta = beta,
    sel_alpha = unlist(sel_alpha, use.names = FALSE),
    sel_beta = unlist(sel_beta, use.names = FALSE),
    mu = mu, cells = cells_out, att = cells_out$att,
    agg_overall = agg_overall,
    tau = tau, eta = eta, a_U = a_U,
    n_active_cells = C_cells,
    n_selection_blocks = length(block_groups),
    cat_info = cat_info
  )
}

sp_ref_public_reference <- function(raw, est_method = "ipw",
                                    selection_model = "pair_logit",
                                    selection_pooling = "separate",
                                    route = "F",
                                    sel_intercept = TRUE,
                                    sel_zero_col = FALSE,
                                    Q_slope = cbind(x = raw$x),
                                    or_intercept = TRUE,
                                    X_slope = cbind(x = raw$x),
                                    V_treat = cbind(`(Intercept)` = 1, x = raw$x),
                                    dates = 1:4,
                                    control_group = "nevertreated",
                                    anticipation = 0L,
                                    pt_assumption = "post",
                                    pt_start = NULL,
                                    base_period = "varying",
                                    eps = 5e-5) {
  cat_info <- cp06_pooled_public_catalog(
    raw, dates = dates, control_group = control_group,
    anticipation = anticipation, pt_assumption = pt_assumption,
    pt_start = pt_start, base_period = base_period
  )
  ev0 <- sp_ref_public_eval(
    raw, w = raw$w, est_method = est_method,
    selection_model = selection_model, selection_pooling = selection_pooling,
    route = route, sel_intercept = sel_intercept, sel_zero_col = sel_zero_col,
    Q_slope = Q_slope, or_intercept = or_intercept, X_slope = X_slope,
    V_treat = V_treat, dates = dates, control_group = control_group,
    anticipation = anticipation, pt_assumption = pt_assumption,
    pt_start = pt_start, base_period = base_period, cat_info = cat_info
  )
  N <- nrow(raw)
  L <- length(ev0$mu)
  R <- length(ev0$att)
  avail <- !is.na(ev0$att)
  avail_agg <- !is.na(ev0$agg_overall)
  phi_mu <- matrix(0, N, L)
  inffunc <- matrix(NA_real_, N, R)
  inffunc[, avail] <- 0
  phi_agg_overall <- matrix(NA_real_, N, length(ev0$agg_overall), dimnames = list(NULL, names(ev0$agg_overall)))
  phi_agg_overall[, avail_agg] <- 0

  for (i in seq_len(N)) {
    wp <- raw$w; wp[i] <- raw$w[i] * (1 + eps)
    wm <- raw$w; wm[i] <- raw$w[i] * (1 - eps)
    up <- sp_ref_public_eval(
      raw, w = wp, est_method = est_method,
      selection_model = selection_model, selection_pooling = selection_pooling,
      route = route, sel_intercept = sel_intercept, sel_zero_col = sel_zero_col,
      Q_slope = Q_slope, or_intercept = or_intercept, X_slope = X_slope,
      V_treat = V_treat, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption,
      pt_start = pt_start, base_period = base_period, cat_info = cat_info
    )
    um <- sp_ref_public_eval(
      raw, w = wm, est_method = est_method,
      selection_model = selection_model, selection_pooling = selection_pooling,
      route = route, sel_intercept = sel_intercept, sel_zero_col = sel_zero_col,
      Q_slope = Q_slope, or_intercept = or_intercept, X_slope = X_slope,
      V_treat = V_treat, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption,
      pt_start = pt_start, base_period = base_period, cat_info = cat_info
    )
    phi_mu[i, ] <- N * (up$mu - um$mu) / (2 * eps)
    inffunc[i, avail] <- N * (up$att[avail] - um$att[avail]) / (2 * eps)
    if (any(avail_agg)) {
      phi_agg_overall[i, avail_agg] <- N * (up$agg_overall[avail_agg] - um$agg_overall[avail_agg]) / (2 * eps)
    }
  }
  V_analytical <- matrix(NA_real_, R, R)
  if (any(avail)) {
    centered <- scale(inffunc[, avail, drop = FALSE], center = TRUE, scale = FALSE)
    V_analytical[avail, avail] <- crossprod(centered) / N
  }
  c(ev0, list(phi_mu = phi_mu, inffunc = inffunc, phi_agg_overall = phi_agg_overall, V_analytical = V_analytical))
}

