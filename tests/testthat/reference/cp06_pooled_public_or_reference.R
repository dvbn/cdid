# Independent pooled dense-functional reference for public control-mean OR (D47).
# Uses only direct stacked weighted block-dummy + duration-scaled slope WLS
# (stats::lm.wfit / QR normal equations), an independent Newton score root for
# Corollary 1 selection logits, incidence least squares for identity GMM, and
# central finite-difference multiplicative case-weight derivatives
#   phi_i(psi) = N * (psi(w * (1 + eps * e_i)) - psi(w * (1 - eps * e_i))) / (2 * eps).
# No production cdid fit, summary, replay, IF, catalog or GMM helper is called.
# This fixed-fixture oracle assumes supported positive weights, interior finite
# selection roots and identified dense designs; it is not an input/rank validator.
# Finite-difference rows and their column means are numerical approximations,
# not exact analytic-centering certificates. Aggregate maps cover the four-date,
# anticipation-zero fixture with every requested post cell available.

cp06_pooled_public_catalog <- function(raw, dates = 1:4, control_group = "nevertreated",
                                       anticipation = 0L, pt_assumption = "post", pt_start = NULL,
                                       base_period = "varying") {
  cal <- sort(unique(as.numeric(dates)))
  T_len <- length(cal)
  idx_time <- function(tt) match(tt, cal)
  Y_mat <- as.matrix(raw[paste0("Y", seq_len(T_len))])
  R_mat <- !is.na(Y_mat)
  A_vec <- rowSums(R_mat) >= 2L

  cohorts_all <- sort(unique(as.numeric(raw$cohort[A_vec])))
  finite_cohorts <- cohorts_all[is.finite(cohorts_all)]
  target_cohorts <- finite_cohorts[vapply(finite_cohorts, function(g) any(cal >= g), logical(1))]
  max_coh <- max(cohorts_all)
  target_cohorts <- target_cohorts[target_cohorts < max_coh]

  has_pair <- function(k, is, it) {
    any(A_vec & raw$cohort == k & R_mat[, is] & R_mat[, it])
  }
  ctrl_pool <- function(g, t_val) {
    if (identical(control_group, "nevertreated")) {
      cohorts_all[is.infinite(cohorts_all)]
    } else {
      cohorts_all[cohorts_all != g & cohorts_all > t_val + anticipation]
    }
  }

  cat_list <- list()
  for (g in target_cohorts) {
    if (!any(cal + anticipation < g)) next
    h0_g <- which(cal >= g)[1L]
    for (is in seq_len(T_len - 1L)) {
      for (it in (is + 1L):T_len) {
        if (!has_pair(g, is, it)) next
        s_val <- cal[is]; t_val <- cal[it]
        if (identical(control_group, "nevertreated")) {
          K_sup <- cohorts_all[is.infinite(cohorts_all) & vapply(cohorts_all, has_pair, logical(1), is = is, it = it)]
          if (!length(K_sup)) next
          cat_list[[length(cat_list) + 1L]] <- data.frame(
            role = "catalog", g = g, s = s_val, t = t_val, is = is, it = it,
            controls = paste(K_sup, collapse = "|"), stringsAsFactors = FALSE
          )
        } else {
          cutoffs <- c(it, seq.int(max(h0_g, it), T_len))
          pools_seen <- character(0)
          cand_pools <- lapply(cutoffs, function(hc) {
            cut_t <- cal[hc]
            cohorts_all[(is.infinite(cohorts_all) | (cohorts_all > g & cohorts_all > cut_t + anticipation)) &
                        vapply(cohorts_all, has_pair, logical(1), is = is, it = it)]
          })
          ord_c <- order(-lengths(cand_pools))
          for (idx_c in ord_c) {
            K_sup <- cand_pools[[idx_c]]
            if (!length(K_sup)) next
            p_str <- paste(K_sup, collapse = "|")
            if (p_str %in% pools_seen) next
            pools_seen <- c(pools_seen, p_str)
            cat_list[[length(cat_list) + 1L]] <- data.frame(
              role = "catalog", g = g, s = s_val, t = t_val, is = is, it = it,
              controls = p_str, stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  cat_df <- do.call(rbind, cat_list)

  aux_list <- list()
  seen_aux <- character(0)
  add_aux <- function(g, s_val, t_val) {
    is <- idx_time(s_val); it <- idx_time(t_val)
    if (is.na(is) || is.na(it) || is >= it) return(NULL)
    K <- ctrl_pool(g, t_val)
    key <- paste0(g, ":", s_val, ":", t_val, ":", paste(K, collapse = ","))
    if (key %in% seen_aux) return(NULL)
    seen_aux <<- c(seen_aux, key)
    if (!has_pair(g, is, it) || !length(K) || !all(vapply(K, has_pair, logical(1), is = is, it = it))) return(NULL)
    aux_list[[length(aux_list) + 1L]] <<- data.frame(
      role = "auxiliary", g = g, s = s_val, t = t_val, is = is, it = it,
      controls = paste(K, collapse = "|"), stringsAsFactors = FALSE
    )
  }
  for (bp in c("varying", "universal")) {
    for (g in target_cohorts) {
      ref_val <- max(cal[cal + anticipation < g])
      times <- if (bp == "varying") cal[-1L] else cal
      for (u in times[times < g]) {
        from_val <- if (bp == "varying") cal[idx_time(u) - 1L] else ref_val
        if (bp == "universal" && u == ref_val) next
        add_aux(g, min(u, from_val), max(u, from_val))
      }
    }
  }
  aux_df <- if (length(aux_list)) do.call(rbind, aux_list) else cat_df[0, , drop = FALSE]
  links_df <- rbind(cat_df, aux_df)
  rownames(links_df) <- NULL
  list(links_df = links_df, cat_df = cat_df, aux_df = aux_df,
       calendar = cal, target_cohorts = target_cohorts, cohorts_all = cohorts_all)
}

cp06_pooled_public_cells_from_mu <- function(cat_info, mu_all, anticipation = 0L,
                                             pt_assumption = "post", pt_start = NULL,
                                             base_period = "varying") {
  cal <- cat_info$calendar
  target_cohorts <- cat_info$target_cohorts
  cat_df <- cat_info$cat_df
  aux_df <- cat_info$aux_df
  M_cat <- nrow(cat_df)
  mu_cat <- mu_all[seq_len(M_cat)]
  mu_aux <- if (nrow(aux_df)) mu_all[M_cat + seq_len(nrow(aux_df))] else numeric(0)

  ps <- if (identical(pt_assumption, "all")) {
    if (is.null(pt_start)) min(cal) else as.numeric(pt_start)
  } else NA_real_

  coord_name <- function(g, u, controls_str) {
    ref_g <- max(cal[cal + anticipation < g])
    K_num <- as.numeric(strsplit(controls_str, "|", fixed = TRUE)[[1L]])
    kmin <- min(K_num)
    elig <- kmin > u + anticipation
    if (u == ref_g) return("ZERO")
    if (identical(pt_assumption, "all") && u >= ps && u < ref_g && elig) return("ZERO")
    if (u > ref_g && elig) return(paste0("shared:", g, ":", u))
    paste0("priv:", g, ":", controls_str, ":", u)
  }

  s_names <- vapply(seq_len(M_cat), function(r) coord_name(cat_df$g[r], cat_df$s[r], cat_df$controls[r]), character(1))
  t_names <- vapply(seq_len(M_cat), function(r) coord_name(cat_df$g[r], cat_df$t[r], cat_df$controls[r]), character(1))
  all_coords <- unique(c(s_names, t_names))
  all_coords <- all_coords[all_coords != "ZERO"]
  A_gmm <- matrix(0, M_cat, length(all_coords), dimnames = list(NULL, all_coords))
  for (r in seq_len(M_cat)) {
    if (s_names[r] != "ZERO") A_gmm[r, s_names[r]] <- A_gmm[r, s_names[r]] - 1
    if (t_names[r] != "ZERO") A_gmm[r, t_names[r]] <- A_gmm[r, t_names[r]] + 1
  }
  psi_hat <- as.vector(qr.solve(crossprod(A_gmm), crossprod(A_gmm, mu_cat)))
  names(psi_hat) <- all_coords

  rows_out <- list()
  for (g in target_cohorts) {
    ref_g <- max(cal[cal + anticipation < g])
    for (u in cal) {
      if (u < g) {
        if (identical(base_period, "varying")) {
          iu <- match(u, cal)
          if (iu == 1L) next
          s_prev <- cal[iu - 1L]
          hit <- which(aux_df$g == g & aux_df$s == s_prev & aux_df$t == u)
          val <- if (length(hit)) mu_aux[hit[1L]] else NA_real_
          rows_out[[length(rows_out) + 1L]] <- data.frame(g = g, time = u, role = "pre", att = val, stringsAsFactors = FALSE)
        } else {
          if (u == ref_g) {
            rows_out[[length(rows_out) + 1L]] <- data.frame(g = g, time = u, role = "reference", att = 0.0, stringsAsFactors = FALSE)
          } else if (u < ref_g) {
            hit <- which(aux_df$g == g & aux_df$s == u & aux_df$t == ref_g)
            val <- if (length(hit)) -mu_aux[hit[1L]] else NA_real_
            rows_out[[length(rows_out) + 1L]] <- data.frame(g = g, time = u, role = "pre", att = val, stringsAsFactors = FALSE)
          } else {
            cn <- paste0("shared:", g, ":", u)
            val <- if (cn %in% names(psi_hat)) unname(psi_hat[cn]) else NA_real_
            rows_out[[length(rows_out) + 1L]] <- data.frame(g = g, time = u, role = "pre", att = val, stringsAsFactors = FALSE)
          }
        }
      } else {
        cn <- paste0("shared:", g, ":", u)
        val <- if (cn %in% names(psi_hat)) unname(psi_hat[cn]) else NA_real_
        rows_out[[length(rows_out) + 1L]] <- data.frame(g = g, time = u, role = "post", att = val, stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows_out)
  rownames(out) <- NULL
  out
}

cp06_pooled_public_dense_eval <- function(raw, w = raw$w, route = "none", intercept = TRUE,
                                          slope_mat = cbind(x = raw$x),
                                          Q_mat = cbind(`(Intercept)` = 1, x = raw$x),
                                          dates = 1:4, control_group = "nevertreated",
                                          anticipation = 0L, pt_assumption = "post",
                                          pt_start = NULL, base_period = "varying",
                                          cat_info = NULL) {
  if (is.null(cat_info)) {
    cat_info <- cp06_pooled_public_catalog(raw, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption, pt_start = pt_start,
      base_period = base_period)
  }
  cal <- cat_info$calendar
  links_df <- cat_info$links_df
  N <- nrow(raw)
  idx_time <- function(tt) match(tt, cal)
  R_mat <- !is.na(as.matrix(raw[paste0("Y", seq_along(cal))]))
  A_vec <- as.numeric(rowSums(R_mat) >= 2L)
  B_vec <- as.numeric(rowSums(R_mat) >= 1L)
  U_vec <- if (route == "none" || route == "F") rep(1, N) else if (route == "B") B_vec else A_vec
  rho_vec <- if (route %in% c("B", "A")) raw$rho else rep(1, N)

  cohorts_all <- cat_info$cohorts_all
  mass_F <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k]), numeric(1))
  mass_U <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & U_vec == 1]), numeric(1))
  mass_A <- vapply(cohorts_all, function(k) sum(w[raw$cohort == k & A_vec == 1]), numeric(1))
  tau <- mass_F / sum(mass_F)
  eta <- mass_A / sum(mass_A)
  a_U <- ifelse(mass_U > 0, mass_A / mass_U, 0)
  names(tau) <- names(eta) <- names(a_U) <- as.character(cohorts_all)

  upairs <- unique(data.frame(s = links_df$s, t = links_df$t,
    is = idx_time(links_df$s), it = idx_time(links_df$t), stringsAsFactors = FALSE))
  upairs <- upairs[order(upairs$is, upairs$it), , drop = FALSE]
  rownames(upairs) <- NULL
  upairs$pair_id <- seq_len(nrow(upairs))
  upairs$L <- upairs$it - upairs$is

  link_pair_id <- match(paste0(links_df$s, "_", links_df$t), paste0(upairs$s, "_", upairs$t))
  block_df <- unique(data.frame(pair_id = link_pair_id, controls = links_df$controls, stringsAsFactors = FALSE))
  pool_sort_key <- vapply(block_df$controls, function(cs) {
    ks <- match(as.numeric(strsplit(cs, "|", fixed = TRUE)[[1L]]), cohorts_all)
    paste(sprintf("%012d", ks), collapse = "|")
  }, character(1))
  block_df <- block_df[order(sprintf("%012d|%s", block_df$pair_id, pool_sort_key), method = "radix"), , drop = FALSE]
  rownames(block_df) <- NULL
  block_df$block_id <- seq_len(nrow(block_df))
  G <- nrow(block_df)

  link_block_id <- match(paste0(link_pair_id, "|", links_df$controls),
                         paste0(block_df$pair_id, "|", block_df$controls))

  pair_cohort_weights <- function(s_idx, t_idx, k) {
    S <- as.numeric(R_mat[, s_idx] & R_mat[, t_idx])
    ik <- which(raw$cohort == k)
    obs <- which(raw$cohort == k & S == 1)
    if (route == "none") {
      rate <- sum(w[ik] * S[ik]) / sum(w[ik] * A_vec[ik])
      list(rows = obs, q = w[obs] / rate)
    } else {
      a_k <- sum(w[ik] * A_vec[ik]) / sum(w[ik] * U_vec[ik])
      iu <- which(raw$cohort == k & U_vec == 1)
      Xk <- Q_mat[iu, , drop = FALSE]
      yk <- S[iu]
      wk <- w[iu]
      gamma <- numeric(ncol(Xk))
      for (iter in seq_len(25L)) {
        pk <- stats::plogis(as.vector(Xk %*% gamma))
        grad <- as.vector(crossprod(Xk, wk * (yk - pk)))
        if (max(abs(grad)) <= 1e-15) break
        hess <- crossprod(Xk, (wk * pk * (1 - pk)) * Xk)
        step <- as.vector(qr.solve(hess, grad, tol = 1e-15))
        gamma <- gamma + step
        if (max(abs(step)) <= 1e-15) break
      }
      eta_lp <- as.vector(Q_mat[obs, , drop = FALSE] %*% gamma)
      p_hat <- stats::plogis(eta_lp)
      list(rows = obs, q = w[obs] * a_k / (rho_vec[obs] * p_hat))
    }
  }

  p <- if (is.null(slope_mat)) 0L else ncol(slope_mat)
  XX <- NULL; YY <- numeric(0); WW <- numeric(0)
  for (b in seq_len(G)) {
    pid <- block_df$pair_id[b]
    s_idx <- upairs$is[pid]; t_idx <- upairs$it[pid]; L_b <- upairs$L[pid]
    ctrl_cohorts <- as.numeric(strsplit(block_df$controls[b], "|", fixed = TRUE)[[1L]])
    for (k in ctrl_cohorts) {
      pw <- pair_cohort_weights(s_idx, t_idx, k)
      rows <- pw$rows
      E <- if (intercept) {
        m <- matrix(0, length(rows), G); m[, b] <- 1; m
      } else matrix(0, length(rows), 0L)
      D <- if (p > 0L) L_b * slope_mat[rows, , drop = FALSE] else matrix(0, length(rows), 0L)
      XX <- rbind(XX, cbind(E, D))
      YY <- c(YY, raw[[paste0("Y", t_idx)]][rows] - raw[[paste0("Y", s_idx)]][rows])
      WW <- c(WW, pw$q)
    }
  }
  wfit <- stats::lm.wfit(x = XX, y = YY, w = WW, tol = 1e-12)
  coef <- unname(wfit$coefficients)
  alpha <- if (intercept) coef[seq_len(G)] else numeric(G)
  beta <- if (p > 0L) coef[(if (intercept) G else 0L) + seq_len(p)] else numeric(0)
  if (p > 0L && !is.null(colnames(slope_mat))) names(beta) <- colnames(slope_mat)

  mu <- vapply(seq_len(nrow(links_df)), function(l) {
    b <- link_block_id[l]
    pid <- link_pair_id[l]
    s_idx <- upairs$is[pid]; t_idx <- upairs$it[pid]; L_l <- upairs$L[pid]
    pw <- pair_cohort_weights(s_idx, t_idx, links_df$g[l])
    rows <- pw$rows
    dy <- raw[[paste0("Y", t_idx)]][rows] - raw[[paste0("Y", s_idx)]][rows]
    pred <- alpha[b] + if (p > 0L) as.vector(L_l * slope_mat[rows, , drop = FALSE] %*% beta) else 0
    sum(pw$q * (dy - pred)) / sum(pw$q)
  }, numeric(1))

  cells <- cp06_pooled_public_cells_from_mu(cat_info, mu, anticipation = anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start, base_period = base_period)

  post <- which(cells$time >= cells$g)
  agg_overall <- c(simple = NA_real_, group = NA_real_, calendar = NA_real_, dynamic = NA_real_)
  if (length(post) && all(!is.na(cells$att[post]))) {
    g_post <- as.character(cells$g[post])
    agg_overall["simple"] <- sum(tau[g_post] * cells$att[post]) / sum(tau[g_post])
    ug <- sort(unique(cells$g[post]))
    est_g <- vapply(ug, function(gv) mean(cells$att[post[cells$g[post] == gv]]), numeric(1))
    agg_overall["group"] <- sum(tau[as.character(ug)] * est_g) / sum(tau[as.character(ug)])
    ut <- sort(unique(cells$time[post]))
    est_t <- vapply(ut, function(tv) {
      rt <- post[cells$time[post] == tv]
      gt <- as.character(cells$g[rt])
      sum(tau[gt] * cells$att[rt]) / sum(tau[gt])
    }, numeric(1))
    agg_overall["calendar"] <- mean(est_t)
    ue <- sort(unique(cells$time[post] - cells$g[post]))
    est_e <- vapply(ue, function(ev) {
      re <- post[(cells$time[post] - cells$g[post]) == ev]
      ge <- as.character(cells$g[re])
      sum(tau[ge] * cells$att[re]) / sum(tau[ge])
    }, numeric(1))
    agg_overall["dynamic"] <- mean(est_e)
  }

  list(alpha = alpha, beta = beta, mu = mu, cells = cells, att = cells$att,
       agg_overall = agg_overall,
       tau = tau, eta = eta, a_U = a_U, blocks = block_df, pairs = upairs,
       cat_info = cat_info)
}

cp06_pooled_public_or_reference <- function(raw, route = "none", intercept = TRUE,
                                            slope_mat = cbind(x = raw$x),
                                            Q_mat = cbind(`(Intercept)` = 1, x = raw$x),
                                            dates = 1:4, control_group = "nevertreated",
                                            anticipation = 0L, pt_assumption = "post",
                                            pt_start = NULL, base_period = "varying",
                                            eps = 1e-5) {
  cat_info <- cp06_pooled_public_catalog(raw, dates = dates, control_group = control_group,
    anticipation = anticipation, pt_assumption = pt_assumption, pt_start = pt_start,
    base_period = base_period)
  ev0 <- cp06_pooled_public_dense_eval(raw, w = raw$w, route = route, intercept = intercept,
    slope_mat = slope_mat, Q_mat = Q_mat, dates = dates, control_group = control_group,
    anticipation = anticipation, pt_assumption = pt_assumption, pt_start = pt_start,
    base_period = base_period, cat_info = cat_info)

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
    up <- cp06_pooled_public_dense_eval(raw, w = wp, route = route, intercept = intercept,
      slope_mat = slope_mat, Q_mat = Q_mat, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption, pt_start = pt_start,
      base_period = base_period, cat_info = cat_info)
    um <- cp06_pooled_public_dense_eval(raw, w = wm, route = route, intercept = intercept,
      slope_mat = slope_mat, Q_mat = Q_mat, dates = dates, control_group = control_group,
      anticipation = anticipation, pt_assumption = pt_assumption, pt_start = pt_start,
      base_period = base_period, cat_info = cat_info)
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

  list(beta = ev0$beta, alpha = ev0$alpha, mu = ev0$mu, att = ev0$att,
       agg_overall = ev0$agg_overall,
       cells = ev0$cells, tau = ev0$tau, eta = ev0$eta, a_U = ev0$a_U,
       phi_beta = phi_beta, phi_alpha = phi_alpha, phi_mu = phi_mu,
       inffunc = inffunc, phi_agg_overall = phi_agg_overall,
       phi_tau = phi_tau, phi_eta = phi_eta, phi_a_U = phi_a_U,
       V_analytical = V_analytical, blocks = ev0$blocks, pairs = ev0$pairs,
       cat_info = cat_info)
}
