# Independent dense-small public reference for D48 CP07 chained doubly robust (DR)
# estimation across baseline MTAR (selection = "none") and conditional MTAR
# (pair_logit on F/B/A and independent_period_logit on F with separate,
# within_cohort, and across_cohorts pooling).
#
# Uses base R (stats::glm.fit and stats::lm.wfit) and reuses the independent
# catalog/reporting mapping from cp06_pooled_public_or_reference.R.

if (!exists("cp06_pooled_public_catalog", mode = "function")) {
  ref_path <- if (requireNamespace("testthat", quietly = TRUE)) {
    testthat::test_path("reference", "cp06_pooled_public_or_reference.R")
  } else {
    file.path("tests", "testthat", "reference", "cp06_pooled_public_or_reference.R")
  }
  sys.source(ref_path, envir = environment())
}

cp07_public_dr_fit_logit <- function(X, y, ww, root_tol = 2e-11) {
  if (ncol(X) == 0L) {
    return(numeric(0))
  }
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
  stopifnot(all(pr > 0 & pr < 1), max(abs(crossprod(X, ww * (y - pr)))) < root_tol)
  unname(as.double(f$coefficients))
}

cp07_public_dr_dense_eval <- function(raw, w = raw$w, route = "none",
                                      selection_model = "pair_logit",
                                      selection_pooling = "separate",
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

  cohorts_all <- cat_info$cohorts_all
  K_len <- length(cohorts_all)
  mass_F <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k]), numeric(1))
  mass_A <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & A_vec == 1]), numeric(1))
  tau <- mass_F / sum(mass_F)
  eta <- mass_A / sum(mass_A)
  names(tau) <- names(eta) <- as.character(cohorts_all)

  sel_alpha <- list()
  sel_beta <- list()
  C_cells <- 0L
  n_sel_blocks <- 0L

  if (identical(route, "none")) {
    a_U <- eta
    get_q_weights <- function(rows, is, it, kv) {
      # DR_SPEC: r_e,k = H_N(w I_k S_e) / H_N(w I_k A), v_e,k = 1/r_e,k.
      # A uses every stored date. Recompute both weighted masses at each
      # case-weight perturbation so the baseline rate path is differentiated.
      pair_rows <- raw$cohort == kv & R_mat[, is] & R_mat[, it]
      rate <- sum(w[pair_rows]) / mass_A[match(kv, cohorts_all)]
      w[rows] / rate
    }
  } else {
    U_vec <- if (route == "F") rep(1, N) else if (route == "B") B_vec else A_vec
    rho_vec <- if (route %in% c("B", "A")) raw$rho else rep(1, N)
    mass_U <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & U_vec == 1]), numeric(1))
    a_U <- ifelse(mass_U > 0, mass_A / mass_U, 0)
    names(a_U) <- as.character(cohorts_all)

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
    n_sel_blocks <- length(block_groups)

    Z_sl <- if (sel_zero_col || is.null(Q_slope)) matrix(0, N, 0L) else as.matrix(Q_slope)
    d_s <- ncol(Z_sl)
    cell_prob <- matrix(0, N, C_cells)
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
      coef_blk <- cp07_public_dr_fit_logit(do.call(rbind, X_rows), y_vec, w_vec)
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
    get_q_weights <- function(rows, is, it, kv) {
      pp <- get_pair_prob(is, it, kv)[rows]
      w[rows] * a_U[as.character(kv)] / (rho_vec[rows] * pp)
    }
  }

  # 1. Shared pooled one-period outcome-change OR across unique control blocks
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
  if (!or_intercept && p_or == 0L) {
    alpha <- numeric(G)
    beta <- numeric(0)
  } else {
    XX <- NULL; YY <- numeric(0); WW <- numeric(0)
    for (b in seq_len(G)) {
      pid <- block_df$pair_id[b]
      is <- upairs$is[pid]; it <- upairs$it[pid]; L_b <- upairs$L[pid]
      ctrl_cohs <- as.numeric(strsplit(block_df$controls[b], "|", fixed = TRUE)[[1L]])
      for (kv in ctrl_cohs) {
        rows <- which(raw$cohort == kv & R_mat[, is] & R_mat[, it])
        qw <- get_q_weights(rows, is, it, kv)
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
  }

  # 2. Treatment propensity logit and residualized Hájek DR link functional
  L_links <- nrow(links_df)
  mu <- numeric(L_links)
  trt_beta <- vector("list", L_links)
  V_tr_mat <- if (is.null(V_treat)) matrix(0, N, 0L) else as.matrix(V_treat)

  for (l in seq_len(L_links)) {
    b <- link_block_id[l]
    pid <- link_pair_id[l]
    is <- upairs$is[pid]; it <- upairs$it[pid]; L_l <- upairs$L[pid]
    g_val <- links_df$g[l]
    ctrl_cohs <- sort(as.numeric(strsplit(links_df$controls[l], "|", fixed = TRUE)[[1L]]))
    cohs <- c(g_val, ctrl_cohs)
    S_pair <- as.numeric(R_mat[, is] & R_mat[, it])
    obs_rows <- which(raw$cohort %in% cohs & S_pair == 1)
    D <- as.double(raw$cohort[obs_rows] == g_val)
    dy <- raw[[paste0("Y", it)]][obs_rows] - raw[[paste0("Y", is)]][obs_rows]
    pred <- alpha[b] + if (p_or > 0L) as.vector(L_l * X_slope[obs_rows, , drop = FALSE] %*% beta) else 0
    res <- dy - pred

    fw <- numeric(length(obs_rows))
    for (kv in cohs) {
      jj <- which(raw$cohort[obs_rows] == kv)
      rr <- obs_rows[jj]
      fw[jj] <- get_q_weights(rr, is, it, kv)
    }
    X_tr <- V_tr_mat[obs_rows, , drop = FALSE]
    b_tr <- cp07_public_dr_fit_logit(X_tr, D, fw)
    trt_beta[[l]] <- b_tr
    odds <- if (ncol(X_tr) > 0L) exp(drop(X_tr %*% b_tr)) else rep(1.0, length(obs_rows))
    wp <- fw * D
    wm <- fw * (1 - D) * odds
    mu[l] <- sum(wp * res) / sum(wp) - sum(wm * res) / sum(wm)
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
    alpha = alpha,
    beta = beta,
    trt_beta = trt_beta,
    sel_alpha = unlist(sel_alpha, use.names = FALSE),
    sel_beta = unlist(sel_beta, use.names = FALSE),
    mu = mu,
    cells = cells_out,
    att = cells_out$att,
    agg_overall = agg_overall,
    tau = tau,
    eta = eta,
    a_U = a_U,
    n_active_cells = C_cells,
    n_selection_blocks = n_sel_blocks,
    blocks = block_df,
    pairs = upairs,
    cat_info = cat_info
  )
}

cp07_public_dr_reference <- function(raw, route = "none",
                                     selection_model = "pair_logit",
                                     selection_pooling = "separate",
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
  ev0 <- cp07_public_dr_dense_eval(
    raw, w = raw$w, route = route,
    selection_model = selection_model, selection_pooling = selection_pooling,
    sel_intercept = sel_intercept, sel_zero_col = sel_zero_col,
    Q_slope = Q_slope, or_intercept = or_intercept, X_slope = X_slope,
    V_treat = V_treat, dates = dates, control_group = control_group,
    anticipation = anticipation, pt_assumption = pt_assumption,
    pt_start = pt_start, base_period = base_period, cat_info = cat_info
  )

  N <- nrow(raw)
  p <- length(ev0$beta)
  G <- length(ev0$alpha)
  L <- length(ev0$mu)
  R <- length(ev0$att)
  C <- length(ev0$tau)

  phi_beta <- matrix(0, N, p, dimnames = list(NULL, names(ev0$beta)))
  phi_alpha <- matrix(0, N, G)
  phi_mu <- matrix(0, N, L)
  inffunc <- matrix(NA_real_, N, R)
  avail <- !is.na(ev0$att)
  inffunc[, avail] <- 0
  phi_tau <- matrix(0, N, C, dimnames = list(NULL, names(ev0$tau)))
  phi_eta <- matrix(0, N, C, dimnames = list(NULL, names(ev0$eta)))
  phi_a_U <- matrix(0, N, C, dimnames = list(NULL, names(ev0$a_U)))
  phi_agg_overall <- matrix(NA_real_, N, length(ev0$agg_overall),
    dimnames = list(NULL, names(ev0$agg_overall)))
  avail_agg <- !is.na(ev0$agg_overall)
  phi_agg_overall[, avail_agg] <- 0

  for (i in seq_len(N)) {
    wp <- raw$w; wp[i] <- raw$w[i] * (1 + eps)
    wm <- raw$w; wm[i] <- raw$w[i] * (1 - eps)
    up <- cp07_public_dr_dense_eval(
      raw, w = wp, route = route,
      selection_model = selection_model, selection_pooling = selection_pooling,
      sel_intercept = sel_intercept, sel_zero_col = sel_zero_col,
      Q_slope = Q_slope, or_intercept = or_intercept, X_slope = X_slope,
      V_treat = V_treat, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption,
      pt_start = pt_start, base_period = base_period, cat_info = cat_info
    )
    um <- cp07_public_dr_dense_eval(
      raw, w = wm, route = route,
      selection_model = selection_model, selection_pooling = selection_pooling,
      sel_intercept = sel_intercept, sel_zero_col = sel_zero_col,
      Q_slope = Q_slope, or_intercept = or_intercept, X_slope = X_slope,
      V_treat = V_treat, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption,
      pt_start = pt_start, base_period = base_period, cat_info = cat_info
    )
    if (p > 0L) phi_beta[i, ] <- N * (up$beta - um$beta) / (2 * eps)
    if (G > 0L) phi_alpha[i, ] <- N * (up$alpha - um$alpha) / (2 * eps)
    phi_mu[i, ] <- N * (up$mu - um$mu) / (2 * eps)
    inffunc[i, avail] <- N * (up$att[avail] - um$att[avail]) / (2 * eps)
    if (any(avail_agg)) {
      phi_agg_overall[i, avail_agg] <- N * (up$agg_overall[avail_agg] - um$agg_overall[avail_agg]) / (2 * eps)
    }
    phi_tau[i, ] <- N * (up$tau - um$tau) / (2 * eps)
    phi_eta[i, ] <- N * (up$eta - um$eta) / (2 * eps)
    phi_a_U[i, ] <- N * (up$a_U - um$a_U) / (2 * eps)
  }

  V_analytical <- matrix(NA_real_, R, R)
  if (any(avail)) {
    centered <- scale(inffunc[, avail, drop = FALSE], center = TRUE, scale = FALSE)
    V_analytical[avail, avail] <- crossprod(centered) / N
  }

  c(
    ev0,
    list(
      phi_beta = phi_beta,
      phi_alpha = phi_alpha,
      phi_mu = phi_mu,
      inffunc = inffunc,
      phi_agg_overall = phi_agg_overall,
      phi_tau = phi_tau,
      phi_eta = phi_eta,
      phi_a_U = phi_a_U,
      V_analytical = V_analytical
    )
  )
}
