# ------------------------------------------------------------------------------
# Cross-Pair Chained DiD Backend
# ------------------------------------------------------------------------------
# Separates per-cohort within-unit outcome change moments of target cohort g and
# clean later-treated / never-treated control cohorts k > g into single directed
# incidence edges in a `cdid_chained_design`, solved jointly by
# `solve_chained_gmm_core()`.
#
# Key properties:
# - Identifies target contrasts ATT(g, t) across disjoint or staggered
#   cohort-by-period observation windows where target cohort g and clean control
#   cohorts need not be observed on the same date pair, under PT-post, PT-ctrl,
#   and PT-all.
# - Under `pt_assumption = "ctrl"`, `gmm_weighting = "within"`, and
#   `base_period = "universal"`, solves the homoskedastic within-unit GLS
#   system with exact mask-group spanning-forest reduction and exact full-catalog
#   Gram reconstruction G = E'E, while extending seamlessly to optimal two-step
#   GMM (`gmm_weighting = "optimal"`), static IPW / OR / DR covariates, survey
#   weights, cluster-robust inference, and simultaneous multiplier-bootstrap
#   confidence bands.

cross_pair_extract_panel_matrices <- function(panel_owner, n_a_units) {
  all_times <- as.double(panel_owner$all_times)
  n_cal <- length(all_times)
  Y_mat <- matrix(NA_real_, nrow = n_a_units, ncol = n_cal)
  O_mat <- matrix(FALSE, nrow = n_a_units, ncol = n_cal)
  post <- panel_owner$postings
  y_obs <- panel_owner$observations$y
  for (c_idx in seq_len(n_cal)) {
    cnt <- post$date_count[c_idx]
    if (cnt > 0L) {
      rng <- post$date_start[c_idx]:post$date_end[c_idx]
      ai <- post$post_a_index[rng]
      oi <- post$post_obs_index[rng]
      O_mat[ai, c_idx] <- TRUE
      Y_mat[ai, c_idx] <- y_obs[oi]
    }
  }
  list(Y_mat = Y_mat, O_mat = O_mat, all_times = all_times, n_cal = n_cal)
}

# Exact spanning-forest reduction + exact integer Gram matrix G = E'E
# across disjoint mask groups. Within each group, shorter date spans (t - s)
# are prioritized so retained witness chains follow consecutive date steps
# whenever adjacent pairs exist, and every dropped link's tree-path expansion
# row in E_drop is accumulated exactly into G_grp = I + E_drop' E_drop.
build_mask_forest_and_gram <- function(s_num, t_num, group_keys) {
  M <- length(s_num)
  if (M == 0L) {
    return(list(kept_idx = integer(0L), dropped_idx = integer(0L), G = matrix(0.0, 0L, 0L)))
  }
  if (M == 1L) {
    return(list(kept_idx = 1L, dropped_idx = integer(0L), G = matrix(1.0, 1L, 1L)))
  }

  grps <- split(seq_len(M), group_keys)
  n_grp <- length(grps)
  kept_chunks <- vector("list", n_grp)
  G_blocks <- vector("list", n_grp)

  for (g_i in seq_len(n_grp)) {
    idx <- grps[[g_i]]
    if (length(idx) == 1L) {
      kept_chunks[[g_i]] <- idx
      G_blocks[[g_i]] <- matrix(1.0, 1L, 1L)
      next
    }
    s_g <- as.double(s_num[idx])
    t_g <- as.double(t_num[idx])
    ord_span <- order(t_g - s_g, s_g, t_g)
    idx_ord <- idx[ord_span]
    s_ord <- s_g[ord_span]
    t_ord <- t_g[ord_span]

    dates_u <- sort(unique(c(s_ord, t_ord)))
    nv <- length(dates_u)
    us_vec <- match(s_ord, dates_u)
    ut_vec <- match(t_ord, dates_u)

    uf_parent <- seq_len(nv)
    find_root <- function(v) {
      while (uf_parent[v] != v) {
        uf_parent[v] <<- uf_parent[uf_parent[v]]
        v <- uf_parent[v]
      }
      v
    }

    m_grp <- length(idx_ord)
    is_tree <- logical(m_grp)
    for (e in seq_len(m_grp)) {
      rs <- find_root(us_vec[e])
      rt <- find_root(ut_vec[e])
      if (rs != rt) {
        uf_parent[rt] <- rs
        is_tree[e] <- TRUE
      }
    }

    tree_pos <- which(is_tree)
    drop_pos <- which(!is_tree)
    k_grp <- length(tree_pos)
    d_grp <- length(drop_pos)
    kept_chunks[[g_i]] <- idx_ord[tree_pos]

    if (d_grp == 0L) {
      G_blocks[[g_i]] <- diag(1.0, k_grp)
      next
    }

    # Build oriented tree adjacency over the k_grp retained edges
    adj_to <- vector("list", nv)
    adj_le <- vector("list", nv)
    adj_sgn <- vector("list", nv)
    for (le in seq_len(k_grp)) {
      u <- us_vec[tree_pos[le]]
      v <- ut_vec[tree_pos[le]]
      # Traversing u -> v adds +edge(le); traversing v -> u adds -edge(le)
      adj_to[[u]] <- c(adj_to[[u]], v)
      adj_le[[u]] <- c(adj_le[[u]], le)
      adj_sgn[[u]] <- c(adj_sgn[[u]], 1.0)
      adj_to[[v]] <- c(adj_to[[v]], u)
      adj_le[[v]] <- c(adj_le[[v]], le)
      adj_sgn[[v]] <- c(adj_sgn[[v]], -1.0)
    }

    visited <- logical(nv)
    depth <- integer(nv)
    par_v <- integer(nv)
    par_le <- integer(nv)
    par_sgn <- numeric(nv)
    queue <- integer(nv)

    for (r0 in seq_len(nv)) {
      if (visited[r0] || is.null(adj_to[[r0]])) next
      visited[r0] <- TRUE
      depth[r0] <- 0L
      head_q <- 1L
      tail_q <- 1L
      queue[1L] <- r0
      while (head_q <= tail_q) {
        curr <- queue[head_q]
        head_q <- head_q + 1L
        nb_v <- adj_to[[curr]]
        nb_le <- adj_le[[curr]]
        nb_s <- adj_sgn[[curr]]
        for (jj in seq_along(nb_v)) {
          nxt <- nb_v[jj]
          if (!visited[nxt]) {
            visited[nxt] <- TRUE
            depth[nxt] <- depth[curr] + 1L
            par_v[nxt] <- curr
            par_le[nxt] <- nb_le[jj]
            # Moving from par_v[nxt] (= curr) to nxt multiplies edge(par_le) by nb_s[jj]
            par_sgn[nxt] <- nb_s[jj]
            tail_q <- tail_q + 1L
            queue[tail_q] <- nxt
          }
        }
      }
    }

    E_drop <- matrix(0.0, nrow = d_grp, ncol = k_grp)
    for (dd in seq_len(d_grp)) {
      e_idx <- drop_pos[dd]
      u <- us_vec[e_idx]
      v <- ut_vec[e_idx]
      while (depth[u] > depth[v]) {
        E_drop[dd, par_le[u]] <- E_drop[dd, par_le[u]] - par_sgn[u]
        u <- par_v[u]
      }
      while (depth[v] > depth[u]) {
        E_drop[dd, par_le[v]] <- E_drop[dd, par_le[v]] + par_sgn[v]
        v <- par_v[v]
      }
      while (u != v) {
        E_drop[dd, par_le[u]] <- E_drop[dd, par_le[u]] - par_sgn[u]
        u <- par_v[u]
        E_drop[dd, par_le[v]] <- E_drop[dd, par_le[v]] + par_sgn[v]
        v <- par_v[v]
      }
    }

    G_blocks[[g_i]] <- diag(1.0, k_grp) + crossprod(E_drop)
  }

  kept_raw <- unlist(kept_chunks, use.names = FALSE)
  L <- length(kept_raw)
  G_mat <- matrix(0.0, nrow = L, ncol = L)
  pos <- 0L
  for (g_i in seq_len(n_grp)) {
    B_g <- G_blocks[[g_i]]
    kg <- nrow(B_g)
    rng <- (pos + 1L):(pos + kg)
    G_mat[rng, rng] <- B_g
    pos <- pos + kg
  }

  ord_k <- order(kept_raw)
  kept <- as.integer(kept_raw[ord_k])
  G_mat <- G_mat[ord_k, ord_k, drop = FALSE]
  # Guarantee exact bit-for-bit IEEE symmetry
  if (L > 1L) {
    G_mat <- 0.5 * (G_mat + t(G_mat))
  }

  list(kept_idx = kept, dropped_idx = setdiff(seq_len(M), kept), G = G_mat)
}

cross_pair_spanning_forest <- function(raw_links) {
  build_mask_forest_and_gram(raw_links$s, raw_links$t, raw_links$mask_key)
}

# Exact, collision-free pair-observation mask partitioning for a single cohort k.
# Uses two BLAS `crossprod` calls on `n_k x n_cal` to identify positive pairs and
# bucket candidate masks, then verifies every multi-pair bucket element-by-element
# against the exact boolean representative mask while caching `obs_idx`.
compute_cohort_pair_masks <- function(O_k, idx_k, n_cal) {
  n_k <- length(idx_k)
  if (n_k == 0L || n_cal < 2L) {
    return(list(cs = integer(0L), ct = integer(0L), mask_id = integer(0L), obs_list = list()))
  }
  if (n_k == 1L) {
    obs_cols <- which(O_k[1L, ])
    if (length(obs_cols) < 2L) {
      return(list(cs = integer(0L), ct = integer(0L), mask_id = integer(0L), obs_list = list()))
    }
    comb <- utils::combn(obs_cols, 2L)
    n_p <- ncol(comb)
    return(list(
      cs = as.integer(comb[1L, ]),
      ct = as.integer(comb[2L, ]),
      mask_id = rep(1L, n_p),
      obs_list = list(as.integer(idx_k))
    ))
  }

  O_dbl <- matrix(as.double(O_k), nrow = n_k, ncol = n_cal)
  cnt_mat <- crossprod(O_dbl)
  ut <- which(upper.tri(cnt_mat) & cnt_mat > 0.5, arr.ind = TRUE)
  if (nrow(ut) == 0L) {
    return(list(cs = integer(0L), ct = integer(0L), mask_id = integer(0L), obs_list = list()))
  }
  ord_ut <- order(ut[, 1L], ut[, 2L])
  cs_v <- as.integer(ut[ord_ut, 1L])
  ct_v <- as.integer(ut[ord_ut, 2L])
  n_p <- length(cs_v)
  pair_lin <- cbind(cs_v, ct_v)
  cnt_v <- as.integer(round(cnt_mat[pair_lin]))

  # Fast shortcut when every positive pair in cohort k is completely balanced (cnt == n_k)
  if (all(cnt_v == n_k)) {
    return(list(
      cs = cs_v,
      ct = ct_v,
      mask_id = rep(1L, n_p),
      obs_list = list(as.integer(idx_k))
    ))
  }

  u_wt <- sin(seq_len(n_k) * 1.618033988749895) + 2.0
  fp_mat <- crossprod(O_dbl, O_dbl * u_wt)
  fp_key <- paste(cnt_v, formatC(fp_mat[pair_lin], digits = 12L, format = "g"), sep = ":")
  buckets <- split(seq_len(n_p), fp_key)

  mask_id <- integer(n_p)
  obs_list <- list()
  cur_id <- 0L

  for (b_idx in buckets) {
    if (length(b_idx) == 1L) {
      p1 <- b_idx[1L]
      cur_id <- cur_id + 1L
      mask_id[p1] <- cur_id
      obs_list[[cur_id]] <- idx_k[which(O_k[, cs_v[p1]] & O_k[, ct_v[p1]])]
      next
    }
    reps_mask <- list()
    reps_id <- integer(0L)
    for (p_i in b_idx) {
      m_bool <- O_k[, cs_v[p_i]] & O_k[, ct_v[p_i]]
      matched_id <- 0L
      for (r_j in seq_along(reps_id)) {
        if (identical(m_bool, reps_mask[[r_j]])) {
          matched_id <- reps_id[r_j]
          break
        }
      }
      if (matched_id == 0L) {
        cur_id <- cur_id + 1L
        matched_id <- cur_id
        reps_mask[[length(reps_mask) + 1L]] <- m_bool
        reps_id <- c(reps_id, matched_id)
        obs_list[[matched_id]] <- idx_k[which(m_bool)]
      }
      mask_id[p_i] <- matched_id
    }
  }

  list(cs = cs_v, ct = ct_v, mask_id = mask_id, obs_list = obs_list)
}

build_cross_pair_link_catalog <- function(O_mat, all_times, a_g, w_a, tgt_cohorts,
                                          cohort_universe, anticipation, control_group,
                                          eligibility = NULL, a_to_full = NULL) {
  n_cal <- length(all_times)
  coh <- sort(unique(as.double(cohort_universe)))
  delta <- as.double(anticipation)

  empty_df <- data.frame(
    g = numeric(0L), s = numeric(0L), t = numeric(0L),
    cs = integer(0L), ct = integer(0L), ref = numeric(0L),
    link_role = character(0L), cohort = numeric(0L),
    coh_idx = integer(0L), coh_mask_id = integer(0L),
    control_pool = character(0L), functional_key = character(0L),
    mask_key = character(0L), stringsAsFactors = FALSE
  )
  if (n_cal < 2L) {
    return(list(links = empty_df, obs_by_coh = list(), units_by_coh = list(), coh_vals = coh))
  }

  n_coh <- length(coh)
  coh_lbls <- vapply(coh, fmt_design_num, character(1L))
  pos_w <- which(w_a > 0)
  coh_match <- match(a_g[pos_w], coh)
  units_by_coh <- vector("list", n_coh)
  for (ci in seq_len(n_coh)) units_by_coh[[ci]] <- integer(0L)
  if (length(pos_w) > 0L) {
    sp_u <- split(pos_w, coh_match)
    for (nm in names(sp_u)) {
      ci <- as.integer(nm)
      if (!is.na(ci) && ci >= 1L && ci <= n_coh) units_by_coh[[ci]] <- as.integer(sp_u[[nm]])
    }
  }

  pair_info <- vector("list", n_coh)
  obs_by_coh <- vector("list", n_coh)
  clean_ct_by_coh <- vector("list", n_coh)
  for (ci in seq_len(n_coh)) {
    idx_k <- units_by_coh[[ci]]
    pm_k <- compute_cohort_pair_masks(O_mat[idx_k, , drop = FALSE], idx_k, n_cal)
    pair_info[[ci]] <- pm_k
    obs_by_coh[[ci]] <- pm_k$obs_list
    k_val <- coh[ci]
    if (is.infinite(k_val)) {
      clean_ct_by_coh[[ci]] <- rep(TRUE, n_cal)
    } else {
      clean_ct_by_coh[[ci]] <- vapply(all_times, function(tv) exact_cmp_sum(tv, delta, k_val) < 0L, logical(1L))
    }
  }

  tgt_sorted <- sort(unique(as.double(tgt_cohorts)))
  chunks <- list()

  for (g in tgt_sorted) {
    ref <- tryCatch(calendar_reference_date(g, all_times, delta), cdid_error = function(e) NA_real_)
    if (is.na(ref)) next
    g_ci <- match(g, coh)
    if (is.na(g_ci) || !length(units_by_coh[[g_ci]])) next
    g_lbl <- coh_lbls[g_ci]
    e <- chained_reference_eligibility_of(eligibility, match(ref, all_times))
    ex <- if (e == 0L) integer() else eligibility$excluded[[match(e, eligibility$index)]]
    pair_info_for <- function(k) {
      idx <- units_by_coh[[k]]
      keep <- !(a_to_full[idx] %in% ex)
      if (all(keep)) return(pair_info[[k]])
      idx <- idx[keep]
      compute_cohort_pair_masks(O_mat[idx, , drop = FALSE], idx, n_cal)
    }

    # 1. Target cohort within-unit pairs (s < t)
    pi_t <- if (length(ex)) pair_info_for(g_ci) else pair_info[[g_ci]]
    if (!length(pi_t$cs) && length(ex)) next
    n_pt <- length(pi_t$cs)
    if (n_pt > 0L) {
      chunks[[length(chunks) + 1L]] <- data.frame(
        g = rep(g, n_pt),
        s = all_times[pi_t$cs],
        t = all_times[pi_t$ct],
        cs = pi_t$cs,
        ct = pi_t$ct,
        ref = rep(ref, n_pt),
        link_role = rep("target", n_pt),
        cohort = rep(g, n_pt),
        coh_idx = rep(g_ci, n_pt),
        coh_mask_id = pi_t$mask_id,
        control_pool = rep(paste0("target:", g_lbl), n_pt),
        functional_key = rep("target", n_pt),
        mask_key = paste0(g_ci, "|target|", pi_t$mask_id),
        stringsAsFactors = FALSE
      )
      if (!is.null(eligibility))
        chunks[[length(chunks)]]$n_observed <- lengths(pi_t$obs_list)[pi_t$mask_id]
    }

    # 2. Clean control cohort within-unit pairs (s < t)
    ctrl_cis <- if (identical(control_group, "nevertreated")) {
      which(is.infinite(coh))
    } else {
      which(coh > g)
    }
    for (k_ci in ctrl_cis) {
      pi_c <- if (length(ex)) pair_info_for(k_ci) else pair_info[[k_ci]]
      if (!length(pi_c$cs)) next
      keep_c <- clean_ct_by_coh[[k_ci]][pi_c$ct]
      if (!any(keep_c)) next
      cs_sub <- pi_c$cs[keep_c]
      ct_sub <- pi_c$ct[keep_c]
      mid_sub <- pi_c$mask_id[keep_c]
      n_pc <- length(cs_sub)
      k_val <- coh[k_ci]
      k_lbl <- coh_lbls[k_ci]
      fk_str <- paste0("ctrl:", k_lbl)
      chunks[[length(chunks) + 1L]] <- data.frame(
        g = rep(g, n_pc),
        s = all_times[cs_sub],
        t = all_times[ct_sub],
        cs = cs_sub,
        ct = ct_sub,
        ref = rep(ref, n_pc),
        link_role = rep("control", n_pc),
        cohort = rep(k_val, n_pc),
        coh_idx = rep(k_ci, n_pc),
        coh_mask_id = mid_sub,
        control_pool = rep(paste0("{", k_lbl, "}"), n_pc),
        functional_key = rep(fk_str, n_pc),
        mask_key = paste0(g_ci, "|ctrl:", k_ci, "|", mid_sub),
        stringsAsFactors = FALSE
      )
      if (!is.null(eligibility))
        chunks[[length(chunks)]]$n_observed <- lengths(pi_c$obs_list)[mid_sub]
    }
  }

  if (!length(chunks)) {
    return(list(links = empty_df, obs_by_coh = obs_by_coh, units_by_coh = units_by_coh, coh_vals = coh))
  }
  out <- do.call(rbind, chunks)
  rownames(out) <- NULL
  list(links = out, obs_by_coh = obs_by_coh, units_by_coh = units_by_coh, coh_vals = coh)
}

build_cross_pair_chained_design <- function(links_df, calendar, cohorts, anticipation,
                                            pt_assumption, pt_start, ledger = NULL) {
  cal <- sort(unique(as.double(calendar)))
  coh <- sort(unique(as.double(cohorts)))
  delta <- as.double(anticipation)
  eff_pt_start <- if (identical(pt_assumption, "ctrl") && is.null(pt_start)) min(cal) else pt_start
  cal_spec <- validate_calendar_spec(cal, coh, delta, pt_assumption = pt_assumption, pt_start = eff_pt_start)

  m <- nrow(links_df)
  row_cohorts <- sort(unique(as.double(links_df$g)))
  nh <- length(row_cohorts)
  refs <- vapply(row_cohorts, function(g) calendar_reference_date(g, cal, delta), numeric(1L))
  n_cal <- length(cal)
  cmp_gt_mat <- matrix(0L, nrow = nh, ncol = n_cal)
  for (r_i in seq_len(nh)) {
    g_r <- row_cohorts[r_i]
    for (c_j in seq_len(n_cal)) {
      cmp_gt_mat[r_i, c_j] <- exact_cmp_sum(cal[c_j], delta, g_r)
    }
  }

  classify_endpoint_vec <- function(g_vec, u_vec, ref_vec, l_role_vec, fk_vec) {
    n_e <- length(g_vec)
    role_out <- character(n_e)
    fk_out <- rep(NA_character_, n_e)
    rk_out <- integer(n_e)

    is_ref <- (u_vec == ref_vec)
    role_out[is_ref] <- "reference"
    rk_out[is_ref] <- 1L

    rem <- which(!is_ref)
    if (length(rem) > 0L) {
      cmp_ant <- cmp_gt_mat[cbind(match(g_vec[rem], row_cohorts), match(u_vec[rem], cal))]
      pre_mask <- cmp_ant < 0L
      post_mask <- !pre_mask

      if (any(pre_mask)) {
        idx_pre <- rem[pre_mask]
        u_pre <- u_vec[idx_pre]
        lr_pre <- l_role_vec[idx_pre]
        fk_pre <- fk_vec[idx_pre]

        if (identical(pt_assumption, "all")) {
          sh <- u_pre >= eff_pt_start
          if (any(sh)) {
            role_out[idx_pre[sh]] <- "pre_shared"
            rk_out[idx_pre[sh]] <- 3L
          }
          nsh <- !sh
          if (any(nsh)) {
            is_t <- lr_pre[nsh] == "target"
            idx_nsh <- idx_pre[nsh]
            role_out[idx_nsh[is_t]] <- "pre_target"
            fk_out[idx_nsh[is_t]] <- "target"
            rk_out[idx_nsh[is_t]] <- 2L
            role_out[idx_nsh[!is_t]] <- "pre"
            fk_out[idx_nsh[!is_t]] <- fk_pre[nsh][!is_t]
            rk_out[idx_nsh[!is_t]] <- 4L
          }
        } else {
          is_t <- lr_pre == "target"
          if (any(is_t)) {
            role_out[idx_pre[is_t]] <- "pre_target"
            fk_out[idx_pre[is_t]] <- "target"
            rk_out[idx_pre[is_t]] <- 2L
          }
          if (any(!is_t)) {
            idx_c <- idx_pre[!is_t]
            u_c <- u_pre[!is_t]
            fk_c <- fk_pre[!is_t]
            if (identical(pt_assumption, "ctrl")) {
              sh_c <- u_c >= eff_pt_start
              role_out[idx_c[sh_c]] <- "pre_shared"
              rk_out[idx_c[sh_c]] <- 3L
              role_out[idx_c[!sh_c]] <- "pre"
              fk_out[idx_c[!sh_c]] <- fk_c[!sh_c]
              rk_out[idx_c[!sh_c]] <- 4L
            } else {
              role_out[idx_c] <- "pre"
              fk_out[idx_c] <- fk_c
              rk_out[idx_c] <- 4L
            }
          }
        }
      }

      if (any(post_mask)) {
        idx_po <- rem[post_mask]
        is_t <- l_role_vec[idx_po] == "target"
        if (any(is_t)) {
          idx_pt <- idx_po[is_t]
          is_ant <- u_vec[idx_pt] < g_vec[idx_pt]
          role_out[idx_pt[is_ant]] <- "anticipation_target"
          fk_out[idx_pt[is_ant]] <- "target"
          rk_out[idx_pt[is_ant]] <- 5L
          role_out[idx_pt[!is_ant]] <- "post_target"
          fk_out[idx_pt[!is_ant]] <- "target"
          rk_out[idx_pt[!is_ant]] <- 7L
        }
        if (any(!is_t)) {
          idx_pc <- idx_po[!is_t]
          role_out[idx_pc] <- "trend_post"
          rk_out[idx_pc] <- 6L
        }
      }
    }

    list(role = role_out, fk = fk_out, kind_rank = rk_out)
  }

  s_ep <- classify_endpoint_vec(links_df$g, links_df$s, links_df$ref, links_df$link_role, links_df$functional_key)
  t_ep <- classify_endpoint_vec(links_df$g, links_df$t, links_df$ref, links_df$link_role, links_df$functional_key)

  sub_led <- if (is.data.frame(ledger) && nrow(ledger) > 0L) {
    ledger[ledger$g %in% row_cohorts, , drop = FALSE]
  } else {
    NULL
  }
  n_led <- if (!is.null(sub_led)) nrow(sub_led) else 0L

  v_g <- c(row_cohorts, links_df$g, links_df$g, if (n_led > 0L) c(sub_led$g, sub_led$g) else numeric(0L))
  v_t <- c(refs, links_df$s, links_df$t, if (n_led > 0L) c(sub_led$time, sub_led$time) else numeric(0L))
  v_role <- c(
    rep("reference", nh), s_ep$role, t_ep$role,
    if (n_led > 0L) c(rep("post_target", n_led), rep("trend_post", n_led)) else character(0L)
  )
  v_fk <- c(
    rep(NA_character_, nh), s_ep$fk, t_ep$fk,
    if (n_led > 0L) c(rep("target", n_led), rep(NA_character_, n_led)) else character(0L)
  )
  v_rank <- c(
    rep(1L, nh), s_ep$kind_rank, t_ep$kind_rank,
    if (n_led > 0L) c(rep(7L, n_led), rep(6L, n_led)) else integer(0L)
  )

  v_key <- paste(
    match(v_g, row_cohorts),
    match(v_t, cal),
    v_role,
    ifelse(is.na(v_fk), "", v_fk),
    sep = "|"
  )
  keep_v <- !duplicated(v_key)
  fk_sort <- ifelse(is.na(v_fk[keep_v]), "", v_fk[keep_v])
  ord_v <- which(keep_v)[order(match(v_g[keep_v], row_cohorts), v_rank[keep_v], v_t[keep_v], fk_sort)]

  vertices <- data.frame(
    g = v_g[ord_v],
    time = v_t[ord_v],
    functional_key = v_fk[ord_v],
    role = v_role[ord_v],
    vertex_index = seq_along(ord_v),
    stringsAsFactors = FALSE
  )
  vk <- v_key[ord_v]
  nv <- nrow(vertices)
  is_anchor_v <- vertices$role == "reference"
  pv <- which(!is_anchor_v)
  p <- length(pv)
  vertices$param_index <- NA_integer_
  if (p > 0L) vertices$param_index[pv] <- seq_len(p)

  g_fmt_v <- fmt_design_num(vertices$g)
  t_fmt_v <- fmt_design_num(vertices$time)
  label <- ifelse(
    vertices$role == "reference",
    paste0("ref(g=", g_fmt_v, ",t=", t_fmt_v, ")"),
    ifelse(
      vertices$role == "post_target",
      paste0("gamma_post(g=", g_fmt_v, ",t=", t_fmt_v, ")"),
      ifelse(
        vertices$role == "anticipation_target",
        paste0("gamma_ant(g=", g_fmt_v, ",t=", t_fmt_v, ")"),
        ifelse(
          vertices$role == "pre_target",
          paste0("gamma_pre(g=", g_fmt_v, ",t=", t_fmt_v, ")"),
          ifelse(
            vertices$role == "trend_post",
            paste0("lambda_post(g=", g_fmt_v, ",t=", t_fmt_v, ")"),
            ifelse(
              vertices$role == "pre_shared",
              paste0("lambda_pre(g=", g_fmt_v, ",t=", t_fmt_v, ")"),
              paste0("lambda_ctrl(g=", g_fmt_v, ",t=", t_fmt_v, ",fn=", vertices$functional_key, ")")
            )
          )
        )
      )
    )
  )
  vertices$vertex_id <- make.unique(label)

  gi <- match(links_df$g, row_cohorts)
  sk <- paste(gi, match(links_df$s, cal), s_ep$role, ifelse(is.na(s_ep$fk), "", s_ep$fk), sep = "|")
  tk <- paste(gi, match(links_df$t, cal), t_ep$role, ifelse(is.na(t_ep$fk), "", t_ep$fk), sep = "|")
  from <- match(sk, vk)
  to <- match(tk, vk)

  edge_order <- seq_len(m)
  anchors <- which(is_anchor_v)
  forest <- incidence_forest(nv, from, to, anchors, edge_order)
  anchored <- unique(forest$component_id[anchors])
  vertices$degree <- forest$degree
  vertices$is_isolated <- forest$degree == 0L
  vertices$component_id <- forest$component_id
  vertices$connected_to_ref <- vertices$component_id %in% anchored
  vertices$anchored_value <- ifelse(is_anchor_v, 0, NA_real_)

  parameters <- vertices[pv, , drop = FALSE]
  rownames(parameters) <- NULL
  parameters$param_id <- parameters$vertex_id
  parameters$level_estimable <- parameters$connected_to_ref

  links <- links_df
  links$link_index <- seq_len(m)
  links$link_id <- paste0("cp_link_", seq_len(m))
  links$start_vertex_index <- from
  links$end_vertex_index <- to
  links$start_param_index <- vertices$param_index[from]
  links$end_param_index <- vertices$param_index[to]
  s_cmp_link <- cmp_gt_mat[cbind(gi, match(links$s, cal))]
  t_cmp_link <- cmp_gt_mat[cbind(gi, match(links$t, cal))]
  links$link_type <- ifelse(
    t_cmp_link < 0L, "pre_only",
    ifelse(s_cmp_link < 0L & links$s != links$ref, "straddling", "post_only")
  )
  links$component_id <- forest$component_id[from]
  rownames(links) <- NULL

  unanchored <- setdiff(seq_len(forest$n_components), anchored)
  rank <- as.integer(p - length(unanchored))
  all_ci <- match(vertices$g, row_cohorts)
  pci <- match(parameters$g, row_cohorts)
  vcount <- tabulate(all_ci, nh)
  ccount <- tabulate(all_ci[forest$roots], nh)
  pcount <- tabulate(pci, nh)
  mcount <- tabulate(gi, nh)
  nullity_coh <- tabulate(all_ci[forest$roots[unanchored]], nh)
  rank_coh <- pcount - nullity_coh

  by_cohort <- data.frame(
    g = row_cohorts, ref = refs, n_links = mcount,
    n_vertices = vcount, n_parameters = pcount,
    n_components = ccount, rank = rank_coh, nullity = nullity_coh,
    cycle_rank = mcount - rank_coh, stringsAsFactors = FALSE
  )

  if (identical(pt_assumption, "ctrl")) {
    ctrl_pre <- links$link_role == "control"
    s_ps <- ctrl_pre & links$s != links$ref & s_cmp_link < 0L & links$s >= eff_pt_start
    t_ps <- ctrl_pre & links$t != links$ref & t_cmp_link < 0L & links$t >= eff_pt_start
    priv_keys <- paste(c(gi[s_ps], gi[t_ps]), c(links$s[s_ps], links$t[t_ps]),
                       c(links$functional_key[s_ps], links$functional_key[t_ps]), sep = ":")
    unmerged_coh <- if (length(priv_keys)) tabulate(c(gi[s_ps], gi[t_ps])[!duplicated(priv_keys)], nh) else integer(nh)
    shared_coh <- tabulate(all_ci[vertices$role == "pre_shared"], nh)
    merged_coh <- pmax(0L, unmerged_coh - shared_coh)
    cal_extra_coh <- vapply(seq_len(nh), function(h) sum(cal >= eff_pt_start & cal < refs[h]), integer(1L))
    by_cohort$pt_start <- rep(eff_pt_start, nh)
    by_cohort$n_extra_calendar_dates <- cal_extra_coh
    by_cohort$n_shared_pre_dates <- shared_coh
    by_cohort$n_merged_pre_parameters <- merged_coh
    by_cohort$n_pre_parameters <- tabulate(pci[parameters$role %in% c("pre", "pre_shared", "pre_target")], nh)
    by_cohort$no_extra_calendar_window <- cal_extra_coh == 0L
    by_cohort$no_extra_restrictions <- merged_coh == 0L
  } else if (identical(pt_assumption, "all")) {
    s_pa <- links$s != links$ref & s_cmp_link < 0L & links$s >= eff_pt_start
    t_pa <- links$t != links$ref & t_cmp_link < 0L & links$t >= eff_pt_start
    priv_keys <- paste(c(gi[s_pa], gi[t_pa]), c(links$s[s_pa], links$t[t_pa]),
                       c(links$functional_key[s_pa], links$functional_key[t_pa]), sep = ":")
    unmerged_coh <- if (length(priv_keys)) tabulate(c(gi[s_pa], gi[t_pa])[!duplicated(priv_keys)], nh) else integer(nh)
    shared_coh <- tabulate(all_ci[vertices$role == "pre_shared"], nh)
    restr_coh <- pmax(0L, unmerged_coh - shared_coh)
    cal_extra_coh <- vapply(seq_len(nh), function(h) sum(cal >= eff_pt_start & cal < refs[h]), integer(1L))
    by_cohort$pt_start <- rep(eff_pt_start, nh)
    by_cohort$n_extra_calendar_dates <- cal_extra_coh
    by_cohort$n_extra_represented_dates <- shared_coh
    by_cohort$n_restricted_pre_parameters <- restr_coh
    by_cohort$n_pre_parameters <- tabulate(pci[parameters$role %in% c("pre", "pre_shared", "pre_target")], nh)
    by_cohort$n_zero_loading_links <- integer(nh)
    by_cohort$no_extra_calendar_window <- cal_extra_coh == 0L
    by_cohort$no_extra_represented_restriction <- restr_coh == 0L
    by_cohort$no_extra_restrictions <- restr_coh == 0L
  }

  rank_diag <- list(
    n_links = m, n_vertices = nv, n_parameters = p,
    n_cohorts = nh, n_components = forest$n_components, rank = rank,
    nullity = length(unanchored), cycle_rank = m - rank,
    has_full_column_rank = length(unanchored) == 0L,
    unanchored_components = unanchored,
    by_cohort = by_cohort
  )

  incidence <- list(
    start = ifelse(is.na(links$start_param_index), 0L, links$start_param_index),
    end = ifelse(is.na(links$end_param_index), 0L, links$end_param_index),
    n_links = m,
    n_parameters = p,
    link_ids = links$link_id,
    param_ids = parameters$param_id
  )

  out <- list(
    incidence = incidence,
    forest = forest,
    parameters = parameters,
    vertices = vertices,
    links = links,
    calendar_spec = cal_spec,
    rank_diagnostics = rank_diag,
    cross_pair = TRUE
  )
  if (identical(pt_assumption, "ctrl")) {
    out$pt_ctrl_spec <- list(
      pt_assumption = "ctrl",
      pt_start = eff_pt_start,
      unrestricted_n_parameters = p + sum(by_cohort$n_merged_pre_parameters),
      restricted_n_parameters = p,
      n_extra_calendar_dates = sum(by_cohort$n_extra_calendar_dates),
      n_shared_pre_dates = sum(by_cohort$n_shared_pre_dates),
      n_merged_pre_parameters = sum(by_cohort$n_merged_pre_parameters),
      no_extra_calendar_window = all(by_cohort$no_extra_calendar_window),
      no_extra_restrictions = all(by_cohort$no_extra_restrictions),
      restriction_is_identity = sum(by_cohort$n_merged_pre_parameters) == 0L,
      by_cohort = by_cohort[, c("g", "ref", "pt_start", "n_extra_calendar_dates",
                                "n_shared_pre_dates", "n_merged_pre_parameters",
                                "n_pre_parameters", "no_extra_calendar_window",
                                "no_extra_restrictions"), drop = FALSE]
    )
  } else if (identical(pt_assumption, "all")) {
    out$pt_all_spec <- list(
      pt_assumption = "all",
      pt_start = eff_pt_start,
      unrestricted_n_parameters = p + sum(by_cohort$n_restricted_pre_parameters),
      restricted_n_parameters = p,
      free_post_param_indices = seq_len(p),
      restricted_post_param_indices = integer(0L),
      zero_loading_link_indices = integer(0L),
      n_extra_calendar_dates = sum(by_cohort$n_extra_calendar_dates),
      n_extra_represented_dates = sum(by_cohort$n_extra_represented_dates),
      n_restricted_pre_parameters = sum(by_cohort$n_restricted_pre_parameters),
      n_zero_loading_links = 0L,
      no_extra_calendar_window = all(by_cohort$no_extra_calendar_window),
      no_extra_represented_restriction = all(by_cohort$no_extra_represented_restriction),
      no_extra_restrictions = all(by_cohort$no_extra_restrictions),
      restriction_is_identity = sum(by_cohort$n_restricted_pre_parameters) == 0L,
      by_cohort = by_cohort[, c("g", "ref", "pt_start", "n_extra_calendar_dates",
                                "n_extra_represented_dates", "n_restricted_pre_parameters",
                                "n_pre_parameters", "n_zero_loading_links",
                                "no_extra_calendar_window", "no_extra_represented_restriction",
                                "no_extra_restrictions"), drop = FALSE]
    )
  }
  structure(out, class = "cdid_chained_design")
}

# Checked propensity-logit solver aligned with `R/chained_ipw_link.R`
fit_checked_propensity_logit <- function(Z_u, D_u, w_u, c_pos, glm_ctrl, n_a, keep_coef = FALSE) {
  p <- ncol(Z_u)
  if (nrow(Z_u) < p) {
    cdid_abort(
      "Observed risk-set pair count is smaller than covariate dimension `ncol(Z)`.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }
  if (any(!is.finite(w_u) | w_u <= 0)) {
    cdid_abort(
      "Nonfinite or nonpositive inverse-rate propensity weights `w_i / r_G`.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  if (ipw_logit_exact_alias(Z_u)) {
    cdid_abort(
      "Weighted propensity logit has deficient model rank.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }
  gf <- ipw_glm_fit_internal(Z_u, D_u, w_u, glm_ctrl)
  fit <- gf$fit
  glm_warnings <- gf$warnings

  has_boundary_warn <- any(grepl("fitted probabilities numerically 0 or 1 occurred", glm_warnings, fixed = TRUE))
  has_conv_warn <- any(grepl("did not converge", glm_warnings, fixed = TRUE))

  if (is.null(fit$rank) || fit$rank < p || is.null(fit$qr$qr) || anyNA(fit$coefficients)) {
    cdid_abort(
      "Weighted propensity logit has deficient model rank.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }
  qr_R <- fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  qr_diag <- abs(diag(qr_R))
  if (any(!is.finite(qr_diag) | qr_diag <= 0)) {
    cdid_abort(
      "Weighted propensity logit QR factor is singular or nonfinite.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }
  if (isTRUE(fit$boundary) || has_boundary_warn) {
    cdid_abort(
      "Weighted propensity logit encountered boundary probabilities (0 or 1).",
      c("cdid_ipw_boundary_error", "cdid_ipw_fit_error")
    )
  }
  if (!isTRUE(fit$converged) || has_conv_warn) {
    cdid_abort(
      "Weighted propensity logit did not converge under caller-supplied `glm.control`.",
      c("cdid_ipw_convergence_error", "cdid_ipw_fit_error")
    )
  }

  beta <- as.double(fit$coefficients)
  if (length(beta) != p || any(!is.finite(beta))) {
    cdid_abort(
      "Weighted propensity logit produced nonfinite coefficients.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }

  separation <- ipw_logit_coordinate_separation(Z_u, D_u, w_u)
  if (!is.null(separation)) {
    sep_diag <- list(
      separation_witness = separation,
      observed_nuisance_only = TRUE,
      identification_claim = "none: an invalid observed nuisance is not ATT nonidentification"
    )
    stop(structure(
      list(
        message = "Weighted propensity logit has no finite root: exact signed-coordinate separation witness.",
        call = NULL,
        diagnostics = sep_diag
      ),
      class = c("cdid_ipw_separation_error", "cdid_ipw_fit_error", "cdid_error", "error", "condition")
    ))
  }

  lv_u <- ipw_logit_values(Z_u, beta)
  p_u <- lv_u$prob
  if (any(!is.finite(lv_u$eta)) || any(!is.finite(lv_u$odds)) ||
      any(lv_u$odds <= 0) || any(!is.finite(p_u))) {
    cdid_abort(
      "Weighted propensity logit produced nonfinite linear predictors, odds, or probabilities.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  prob_eps <- 10 * .Machine$double.eps
  if (any(p_u <= prob_eps | p_u >= 1 - prob_eps | p_u * (1 - p_u) <= 0)) {
    cdid_abort(
      "Weighted propensity logit fitted probabilities reached numerical boundary (0 or 1).",
      c("cdid_ipw_boundary_error", "cdid_ipw_fit_error")
    )
  }

  o_c <- lv_u$odds[c_pos]
  act_c <- lv_u$active[c_pos]
  w_c <- w_u[c_pos]
  b_c <- w_c * o_c
  B_c <- sum(b_c)
  if (!is.finite(B_c) || B_c <= 0) {
    cdid_abort(
      "Invalid control inverse-probability normalization denominator.",
      c("cdid_ipw_normalization_error", "cdid_ipw_fit_error")
    )
  }

  H_mat <- crossprod(Z_u, Z_u * (w_u * p_u * (1.0 - p_u))) / as.double(n_a)
  if (any(!is.finite(H_mat))) {
    cdid_abort("Nonfinite weighted propensity Hessian `H`.", c("cdid_ipw_numeric_error", "cdid_ipw_fit_error"))
  }
  H_chol <- tryCatch(base::chol(H_mat), error = function(e) NULL)
  if (is.null(H_chol) || any(!is.finite(H_chol)) || any(diag(H_chol) <= 0)) {
    cdid_abort("Weighted propensity Hessian `H` is not positive definite.", c("cdid_ipw_rank_error", "cdid_ipw_fit_error"))
  }

  out <- list(
    Z_u = Z_u, D_u = D_u, w_u = w_u, p_u = p_u,
    o_c = o_c, act_c = act_c, b_c = b_c, B_c = B_c,
    H_chol = H_chol
  )
  # Opt-in nuisance capture only: the fitted coefficients of this exact solve.
  if (isTRUE(keep_coef)) {
    out$beta <- beta
    out$diagnostics <- list(converged = isTRUE(fit$converged), iter = fit$iter, rank = fit$rank)
  }
  out
}

# Checked WLS QR factorizer + target-support alias verification aligned with `R/chained_or_wls.R` and `R/chained_or_pooled.R`
fit_checked_wls_qr <- function(X_c, w_c, rank_tol = 1e-7) {
  p <- ncol(X_c)
  if (nrow(X_c) == 0L || p == 0L) {
    cdid_abort("Outcome-regression control design matrix is empty on observed pair support.", "cdid_or_wls_rank_error")
  }
  sq_wc <- sqrt(w_c)
  Xw <- X_c * sq_wc
  qr_c <- qr(Xw, LAPACK = TRUE)
  R_full <- qr.R(qr_c)
  d_len <- min(nrow(R_full), ncol(R_full))
  diag_abs <- if (d_len > 0L) abs(R_full[cbind(seq_len(d_len), seq_len(d_len))]) else numeric(0L)
  max_d <- if (length(diag_abs)) max(diag_abs) else 0.0
  rk <- if (is.finite(max_d) && max_d > 0.0) sum(diag_abs > rank_tol * max_d) else 0L
  if (rk < 1L) {
    cdid_abort("Outcome-regression control design matrix has zero numerical rank on observed pair support.", "cdid_or_wls_rank_error")
  }
  piv <- qr_c$pivot[seq_len(rk)]
  R_rk <- R_full[seq_len(rk), seq_len(rk), drop = FALSE]
  drop_cols <- if (rk < p) qr_c$pivot[(rk + 1L):p] else integer(0L)
  A_alias <- NULL
  err_c <- 0.0
  if (length(drop_cols) > 0L) {
    A_alias <- backsolve(R_rk, qr.qty(qr_c, Xw[, drop_cols, drop = FALSE])[seq_len(rk), , drop = FALSE])
    err_c <- max(abs(X_c[, drop_cols, drop = FALSE] - X_c[, piv, drop = FALSE] %*% A_alias))
  }
  list(
    X_c = X_c, w_c = w_c, sq_wc = sq_wc, qr_c = qr_c,
    rk = rk, p = p, piv = piv, R_rk = R_rk,
    drop_cols = drop_cols, A_alias = A_alias, err_c = err_c, rank_tol = rank_tol
  )
}

verify_wls_target_support <- function(qr_rec, X_t, g_val = NA_real_) {
  if (!length(qr_rec$drop_cols)) return(invisible(TRUE))
  scale_x <- max(1.0, max(abs(qr_rec$X_c)), max(abs(X_t)))
  alias_tol <- 64 * .Machine$double.eps * max(1L, qr_rec$p)
  err_t <- max(abs(X_t[, qr_rec$drop_cols, drop = FALSE] - X_t[, qr_rec$piv, drop = FALSE] %*% qr_rec$A_alias))
  if (!is.finite(qr_rec$err_c) || !is.finite(err_t) || max(qr_rec$err_c, err_t) > alias_tol * scale_x) {
    g_msg <- if (is.finite(g_val)) paste0(" for target cohort g=", fmt_design_num(g_val)) else ""
    cdid_abort(
      paste0(
        "Outcome-regression control design matrix is rank-deficient (", qr_rec$rk, " < ", qr_rec$p,
        ") without an exact linear alias across control and target support", g_msg,
        "; counterfactual mean is not identified without extrapolation."
      ),
      "cdid_or_wls_rank_error"
    )
  }
  invisible(TRUE)
}

evaluate_cross_pair_moments <- function(links_df, Y_mat, a_g, w_a, Z_a,
                                        obs_by_coh, units_by_coh, coh_vals,
                                        est_method, cp04, glm_control, wt_mode,
                                        compute_if = TRUE, nuisance_sink = NULL) {
  L <- nrow(links_df)
  n_a <- nrow(Y_mat)
  mhat <- numeric(L)
  sum_w_a <- sum(w_a)
  inv_sqrt_w <- ifelse(w_a > 0, 1.0 / sqrt(w_a), 0.0)

  is_or <- identical(est_method, "reg")
  is_dr <- identical(est_method, "dr")
  Z_trt <- if (is_dr && !is.null(cp04) && !is.null(cp04$treatment_columns)) {
    Z_a[, cp04$treatment_columns, drop = FALSE]
  } else Z_a
  Z_or <- if ((is_or || is_dr) && !is.null(cp04) && !is.null(cp04$outcome_recipe)) {
    Z_a[, cp04$outcome_recipe$columns, drop = FALSE]
  } else Z_a
  uncond_trt <- ncol(Z_trt) == 1L && all(Z_trt[, 1L] == Z_trt[1L, 1L]) && Z_trt[1L, 1L] != 0
  uncond_or <- ncol(Z_or) == 1L && all(Z_or[, 1L] == Z_or[1L, 1L]) && Z_or[1L, 1L] != 0
  unconditional <- (is_or && uncond_or) || (!is_or && !is_dr && uncond_trt) || (is_dr && uncond_trt && uncond_or)

  make_fast_backend <- function(n_sup, n_reg, n_trt, n_rt) {
    structure(
      list(
        mode = "baseline",
        est_method = est_method,
        outcome_pooling = "separate",
        imputation_pooling = NULL,
        outcome_recipe = if (!is.null(cp04)) cp04$outcome_recipe else NULL,
        reader_ref = NULL,
        reference_ref = NULL,
        uses_selection_models = FALSE,
        selection_backend = NULL,
        accounting = list(
          automatic_chunk_rows = if (!is.null(cp04)) isTRUE(cp04$automatic_chunk_rows) else TRUE,
          n_selection_fits = 0L,
          n_selection_cells = 0L,
          n_selection_blocks = 0L,
          n_selection_parameters = 0L,
          n_regression_fits = as.integer(n_reg),
          n_regression_blocks = as.integer(n_reg),
          n_supports = as.integer(n_sup),
          n_treatment_fits = as.integer(n_trt),
          n_rate_fits = as.integer(n_rt),
          n_reference_recipes = 0L
        )
      ),
      class = "cdid_cross_pair_backend"
    )
  }

  if (unconditional) {
    # Deduplicate (coh_idx, cs, ct) evaluations across target slices g
    u_key <- paste(links_df$coh_idx, links_df$cs, links_df$ct, sep = ":")
    first_pos <- which(!duplicated(u_key))
    u_map <- match(u_key, u_key[first_pos])
    n_u <- length(first_pos)
    u_mhat <- numeric(n_u)
    u_Phi <- if (isTRUE(compute_if)) matrix(0.0, nrow = n_a, ncol = n_u) else NULL
    u_W <- numeric(n_u)

    for (uu in seq_len(n_u)) {
      ell0 <- first_pos[uu]
      ci <- links_df$coh_idx[ell0]
      mid <- links_df$coh_mask_id[ell0]
      cs <- links_df$cs[ell0]
      ct <- links_df$ct[ell0]
      idx <- obs_by_coh[[ci]][[mid]]
      dy <- Y_mat[idx, ct] - Y_mat[idx, cs]
      w_sub <- w_a[idx]
      W_sum <- sum(w_sub)
      omega <- w_sub / W_sum
      m_val <- sum(omega * dy)
      u_mhat[uu] <- m_val
      u_W[uu] <- W_sum
      if (isTRUE(compute_if)) {
        u_Phi[idx, uu] <- as.double(n_a) * omega * (dy - m_val)
      }
    }

    mhat <- u_mhat[u_map]
    if (isTRUE(compute_if)) {
      Phi <- u_Phi[, u_map, drop = FALSE]
      u_Phi <- NULL
      mean_phi <- colMeans(Phi)
      Phi_centered <- sweep(Phi, 2L, mean_phi, "-")
      Phi <- NULL
    } else {
      Phi_centered <- NULL
      mean_phi <- numeric(L)
    }

    if (identical(wt_mode, "within")) {
      # Exact analytical within-unit covariance via cohort-mask weight overlaps
      # (zero L x n_a dense Omega_scaled allocation)
      g_vec <- as.double(links_df$g)
      ci_vec <- as.integer(links_df$coh_idx)
      mid_vec <- as.integer(links_df$coh_mask_id)
      s_vec <- as.double(links_df$s)
      t_vec <- as.double(links_df$t)
      W_vec <- u_W[u_map]

      same_g_coh <- outer(g_vec, g_vec, "==") & outer(ci_vec, ci_vec, "==")
      time_sign <- (outer(t_vec, t_vec, "==") + outer(s_vec, s_vec, "==") -
                      outer(t_vec, s_vec, "==") - outer(s_vec, t_vec, "=="))
      active_mat <- same_g_coh & (time_sign != 0)
      Vw <- matrix(0.0, nrow = L, ncol = L)
      act_pairs <- which(upper.tri(active_mat, diag = TRUE) & active_mat, arr.ind = TRUE)
      if (nrow(act_pairs) > 0L) {
        ov_cache <- new.env(hash = TRUE, parent = emptyenv())
        in_mask <- logical(n_a)
        for (r in seq_len(nrow(act_pairs))) {
          i1 <- act_pairs[r, 1L]
          i2 <- act_pairs[r, 2L]
          ci <- ci_vec[i1]
          m1 <- mid_vec[i1]
          m2 <- mid_vec[i2]
          ov_w <- if (m1 == m2) {
            W_vec[i1]
          } else {
            lo_m <- min(m1, m2); hi_m <- max(m1, m2)
            ck <- paste0(ci, ":", lo_m, ":", hi_m)
            val_cached <- ov_cache[[ck]]
            if (is.null(val_cached)) {
              idx1 <- obs_by_coh[[ci]][[lo_m]]
              idx2 <- obs_by_coh[[ci]][[hi_m]]
              in_mask[idx1] <- TRUE
              val_cached <- sum(w_a[idx2[in_mask[idx2]]])
              in_mask[idx1] <- FALSE
              ov_cache[[ck]] <- val_cached
            }
            val_cached
          }
          v_ij <- time_sign[i1, i2] * sum_w_a * ov_w / (W_vec[i1] * W_vec[i2])
          Vw[i1, i2] <- v_ij
          Vw[i2, i1] <- v_ij
        }
      }
      V_used <- 0.5 * (Vw + t(Vw))
    } else if (identical(wt_mode, "optimal")) {
      V_opt <- crossprod(Phi_centered) / as.double(n_a)
      V_used <- 0.5 * (V_opt + t(V_opt))
    } else {
      V_used <- matrix(0.0, nrow = L, ncol = L)
    }

    if (!is.null(nuisance_sink)) cross_pair_nuisance_baseline_collect(nuisance_sink, environment(), unconditional = TRUE)
    return(list(
      mhat = mhat, Phi = Phi_centered, mean_phi = mean_phi, V = V_used,
      backend = make_fast_backend(n_u, 0L, 0L, 0L)
    ))
  }

  # Covariate-adjusted path with per-mask-group WLS QR and propensity-logit caching
  Phi <- if (isTRUE(compute_if)) matrix(0.0, nrow = n_a, ncol = L) else NULL
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  Omega_scaled <- if (identical(wt_mode, "within")) matrix(0.0, nrow = L, ncol = n_a) else NULL
  coh_mass <- vapply(units_by_coh, function(ix) sum(w_a[ix]), numeric(1L))

  # Cache target cohort summaries once per g
  u_g <- sort(unique(as.double(links_df$g)))
  tgt_cache <- vector("list", length(u_g))
  for (gi in seq_along(u_g)) {
    gv <- u_g[gi]
    ci_t <- match(gv, coh_vals)
    idx_t <- units_by_coh[[ci_t]]
    w_t <- w_a[idx_t]
    W_t <- sum(w_t)
    tgt_cache[[gi]] <- list(
      idx_t = idx_t,
      w_t = w_t,
      W_t = W_t,
      Z_trt_t = Z_trt[idx_t, , drop = FALSE],
      X_t = Z_or[idx_t, , drop = FALSE]
    )
  }

  qr_cache <- new.env(hash = TRUE, parent = emptyenv())
  logit_cache <- new.env(hash = TRUE, parent = emptyenv())
  alias_checked <- new.env(hash = TRUE, parent = emptyenv())

  for (ell in seq_len(L)) {
    cs <- links_df$cs[ell]
    ct <- links_df$ct[ell]
    ci <- links_df$coh_idx[ell]
    mid <- links_df$coh_mask_id[ell]

    if (identical(links_df$link_role[ell], "target")) {
      idx <- obs_by_coh[[ci]][[mid]]
      dy <- Y_mat[idx, ct] - Y_mat[idx, cs]
      w_sub <- w_a[idx]
      W_sum <- sum(w_sub)
      omega <- w_sub / W_sum
      m_val <- sum(omega * dy)
      mhat[ell] <- m_val
      if (isTRUE(compute_if)) {
        Phi[idx, ell] <- as.double(n_a) * omega * (dy - m_val)
      }
      if (!is.null(Omega_scaled)) Omega_scaled[ell, idx] <- omega * inv_sqrt_w[idx]
      next
    }

    g_ell <- links_df$g[ell]
    gi <- match(g_ell, u_g)
    tc <- tgt_cache[[gi]]
    idx_t <- tc$idx_t
    w_t <- tc$w_t
    W_t <- tc$W_t
    idx_c <- obs_by_coh[[ci]][[mid]]
    w_c <- w_a[idx_c]
    dy_c <- Y_mat[idx_c, ct] - Y_mat[idx_c, cs]

    lg_rec <- NULL
    if (!is_or) {
      lg_key <- paste0(gi, ":", ci, ":", mid)
      lg_rec <- logit_cache[[lg_key]]
      if (is.null(lg_rec)) {
        idx_k_all <- units_by_coh[[ci]]
        M_k <- coh_mass[ci]
        is_sub_k <- length(idx_c) < length(idx_k_all)
        r_raw <- if (is_sub_k) sum(w_c) / M_k else 1.0
        r_c <- max(r_raw, ipw_clip_eps)
        r_act <- is_sub_k && (r_raw > ipw_clip_eps)
        v_c <- w_c / r_c
        Z_c_mat <- Z_trt[idx_c, , drop = FALSE]
        if (uncond_trt) {
          b_c_unc <- v_c * (W_t / sum(v_c))
          lg_rec <- list(
            uncond = TRUE,
            b_c = b_c_unc,
            B_c = sum(b_c_unc),
            A_gamma_r = NULL
          )
        } else {
          idx_u <- c(idx_t, idx_c)
          D_u <- c(rep(1.0, length(idx_t)), rep(0.0, length(idx_c)))
          w_u <- c(w_t, v_c)
          Z_u <- rbind(tc$Z_trt_t, Z_c_mat)
          c_pos <- (length(idx_t) + 1L):length(idx_u)
          lg_fit <- fit_checked_propensity_logit(Z_u, D_u, w_u, c_pos, glm_ctrl, n_a,
                                                 keep_coef = !is.null(nuisance_sink))
          if (isTRUE(compute_if) || !is.null(nuisance_sink)) {
            lg_rec <- lg_fit
            lg_rec$uncond <- FALSE
            lg_rec$idx_u <- idx_u
            lg_rec$Z_c <- Z_c_mat
            lg_rec$idx_k_all <- idx_k_all
            lg_rec$M_k <- M_k
            lg_rec$r_raw <- r_raw
            lg_rec$A_gamma_r <- if (r_act) {
              colSums(Z_c_mat * (w_c * lg_fit$p_u[c_pos])) / (as.double(n_a) * r_c * r_c)
            } else NULL
          } else {
            lg_rec <- list(uncond = FALSE, b_c = lg_fit$b_c, B_c = lg_fit$B_c)
          }
        }
        logit_cache[[lg_key]] <- lg_rec
      }
    }

    if (!is_or && !is_dr) {
      b_c <- lg_rec$b_c
      B_c <- lg_rec$B_c
      omega_c <- b_c / B_c
      m_val <- sum(omega_c * dy_c)
      mhat[ell] <- m_val
      if (isTRUE(compute_if)) {
        Phi[idx_c, ell] <- Phi[idx_c, ell] + as.double(n_a) * omega_c * (dy_c - m_val)
        if (!isTRUE(lg_rec$uncond)) {
          G_beta <- colSums(lg_rec$Z_c * (b_c * (dy_c - m_val) * lg_rec$act_c)) / B_c
          h_inv_g <- backsolve(lg_rec$H_chol, backsolve(lg_rec$H_chol, G_beta, transpose = TRUE))
          Phi[lg_rec$idx_u, ell] <- Phi[lg_rec$idx_u, ell] + lg_rec$w_u * (lg_rec$D_u - lg_rec$p_u) * as.double(lg_rec$Z_u %*% h_inv_g)
          if (!is.null(lg_rec$A_gamma_r)) {
            lam_r <- sum(lg_rec$A_gamma_r * h_inv_g)
            ik_all <- lg_rec$idx_k_all
            inv_Mk <- as.double(n_a) / lg_rec$M_k
            Phi[ik_all, ell] <- Phi[ik_all, ell] - (w_a[ik_all] * inv_Mk) * (lg_rec$r_raw * lam_r)
            Phi[idx_c, ell] <- Phi[idx_c, ell] + (w_c * inv_Mk) * lam_r
          }
        }
      }
      if (!is.null(Omega_scaled)) Omega_scaled[ell, idx_c] <- omega_c * inv_sqrt_w[idx_c]
    } else {
      qr_key <- paste0(ci, ":", mid)
      qr_rec <- qr_cache[[qr_key]]
      if (is.null(qr_rec)) {
        qr_rec <- fit_checked_wls_qr(Z_or[idx_c, , drop = FALSE], w_c)
        qr_cache[[qr_key]] <- qr_rec
      }
      al_key <- paste0(gi, ":", qr_key)
      if (is.null(alias_checked[[al_key]])) {
        verify_wls_target_support(qr_rec, tc$X_t, g_ell)
        alias_checked[[al_key]] <- TRUE
      }

      rk <- qr_rec$rk
      piv <- qr_rec$piv
      R_rk <- qr_rec$R_rk
      X_c_piv <- qr_rec$X_c[, piv, drop = FALSE]
      X_t_piv <- tc$X_t[, piv, drop = FALSE]
      gamma_piv <- backsolve(R_rk, qr.qty(qr_rec$qr_c, dy_c * qr_rec$sq_wc)[seq_len(rk)])
      mu_t <- as.double(X_t_piv %*% gamma_piv)
      mu_c <- as.double(X_c_piv %*% gamma_piv)
      e_c <- dy_c - mu_c
      if (!is.null(nuisance_sink)) cross_pair_nuisance_baseline_or(nuisance_sink, environment())
      xbar_t <- colSums(X_t_piv * w_t) / W_t
      mu_bar_t <- sum(w_t * mu_t) / W_t

      if (is_or) {
        mhat[ell] <- mu_bar_t
        if (isTRUE(compute_if) || !is.null(Omega_scaled)) {
          coef_x <- backsolve(R_rk, forwardsolve(t(R_rk), xbar_t))
          omega_c <- w_c * as.double(X_c_piv %*% coef_x)
          if (isTRUE(compute_if)) {
            Phi[idx_t, ell] <- Phi[idx_t, ell] + (as.double(n_a) * w_t / W_t) * (mu_t - mu_bar_t)
            Phi[idx_c, ell] <- Phi[idx_c, ell] + as.double(n_a) * omega_c * e_c
          }
          if (!is.null(Omega_scaled)) Omega_scaled[ell, idx_c] <- omega_c * inv_sqrt_w[idx_c]
        }
      } else {
        b_c <- lg_rec$b_c
        B_c <- lg_rec$B_c
        e_bar_att <- sum(b_c * e_c) / B_c
        mhat[ell] <- mu_bar_t + e_bar_att
        if (isTRUE(compute_if) || !is.null(Omega_scaled)) {
          xbar_att_c <- colSums(X_c_piv * b_c) / B_c
          coef_dr <- backsolve(R_rk, forwardsolve(t(R_rk), xbar_t - xbar_att_c))
          w_or_adj <- w_c * as.double(X_c_piv %*% coef_dr)
          omega_c <- (b_c / B_c) + w_or_adj
          if (isTRUE(compute_if)) {
            Phi[idx_t, ell] <- Phi[idx_t, ell] + (as.double(n_a) * w_t / W_t) * (mu_t - mu_bar_t)
            Phi[idx_c, ell] <- Phi[idx_c, ell] + (as.double(n_a) * b_c / B_c) * (e_c - e_bar_att) +
              as.double(n_a) * w_or_adj * e_c
            if (!isTRUE(lg_rec$uncond)) {
              G_beta <- colSums(lg_rec$Z_c * (b_c * (e_c - e_bar_att) * lg_rec$act_c)) / B_c
              h_inv_g <- backsolve(lg_rec$H_chol, backsolve(lg_rec$H_chol, G_beta, transpose = TRUE))
              Phi[lg_rec$idx_u, ell] <- Phi[lg_rec$idx_u, ell] + lg_rec$w_u * (lg_rec$D_u - lg_rec$p_u) * as.double(lg_rec$Z_u %*% h_inv_g)
              if (!is.null(lg_rec$A_gamma_r)) {
                lam_r <- sum(lg_rec$A_gamma_r * h_inv_g)
                ik_all <- lg_rec$idx_k_all
                inv_Mk <- as.double(n_a) / lg_rec$M_k
                Phi[ik_all, ell] <- Phi[ik_all, ell] - (w_a[ik_all] * inv_Mk) * (lg_rec$r_raw * lam_r)
                Phi[idx_c, ell] <- Phi[idx_c, ell] + (w_c * inv_Mk) * lam_r
              }
            }
          }
          if (!is.null(Omega_scaled)) Omega_scaled[ell, idx_c] <- omega_c * inv_sqrt_w[idx_c]
        }
      }
    }
  }

  if (isTRUE(compute_if)) {
    mean_phi <- colMeans(Phi)
    Phi_centered <- sweep(Phi, 2L, mean_phi, "-")
    Phi <- NULL
  } else {
    Phi_centered <- NULL
    mean_phi <- numeric(L)
  }

  V_used <- if (identical(wt_mode, "within")) {
    g_vec <- as.double(links_df$g)
    s_vec <- as.double(links_df$s)
    t_vec <- as.double(links_df$t)
    cols_by_g <- split(seq_len(L), match(g_vec, unique(g_vec)))
    Vw <- matrix(0.0, nrow = L, ncol = L)
    for (cols_g in cols_by_g) {
      sg <- s_vec[cols_g]
      tg <- t_vec[cols_g]
      ts_g <- (outer(tg, tg, "==") + outer(sg, sg, "==") -
                 outer(tg, sg, "==") - outer(sg, tg, "=="))
      Om_g <- Omega_scaled[cols_g, , drop = FALSE]
      nz_c <- which(colSums(Om_g != 0) > 0L)
      G_g <- if (length(nz_c)) tcrossprod(Om_g[, nz_c, drop = FALSE]) * sum_w_a else
        matrix(0.0, length(cols_g), length(cols_g))
      V_blk <- ts_g * G_g
      Vw[cols_g, cols_g] <- 0.5 * (V_blk + t(V_blk))
    }
    Vw
  } else if (identical(wt_mode, "optimal")) {
    V_opt <- crossprod(Phi_centered) / as.double(n_a)
    0.5 * (V_opt + t(V_opt))
  } else {
    matrix(0.0, nrow = L, ncol = L)
  }

  n_reg_fits <- length(ls(qr_cache, all.names = TRUE))
  n_trt_fits <- length(ls(logit_cache, all.names = TRUE))
  if (!is.null(nuisance_sink)) cross_pair_nuisance_baseline_collect(nuisance_sink, environment())
  list(
    mhat = mhat, Phi = Phi_centered, mean_phi = mean_phi, V = V_used,
    backend = make_fast_backend(L, n_reg_fits, n_trt_fits, n_trt_fits)
  )
}

cross_pair_find_param <- function(params_df, g_val, t_val, roles) {
  hit <- which(params_df$g == g_val & params_df$time == t_val & params_df$role %in% roles)
  if (length(hit) == 1L) as.integer(params_df$param_index[hit]) else NA_integer_
}

build_cross_pair_target_specs <- function(design, ledger, pre_selections = NULL, include_pre = FALSE) {
  params <- design$parameters
  R_post <- nrow(ledger)
  specs <- list()
  meta_rows <- list()

  for (h in seq_len(R_post)) {
    gv <- ledger$g[h]
    tv <- ledger$time[h]
    cid <- ledger$cell_id[h]
    idx_g <- cross_pair_find_param(params, gv, tv, "post_target")
    idx_l <- cross_pair_find_param(params, gv, tv, "trend_post")
    if (!is.na(idx_g) && !is.na(idx_l)) {
      specs[[length(specs) + 1L]] <- gmm_target_spec(cid, c(idx_g, idx_l), c(1.0, -1.0))
      meta_rows[[length(meta_rows) + 1L]] <- data.frame(
        target_id = cid, kind = "post", ledger_index = h, pre_index = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (isTRUE(include_pre) && length(pre_selections) > 0L) {
    for (j in seq_along(pre_selections)) {
      d <- pre_selections[[j]]
      if (identical(d$status, "reference_zero") || identical(d$status, "unavailable")) next
      tid <- paste0("pre(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")")
      p_idx <- integer(0L); coeff <- numeric(0L)
      if (d$t != d$ref) {
        ig_t <- cross_pair_find_param(params, d$g, d$t, c("pre_target", "anticipation_target"))
        il_t <- cross_pair_find_param(params, d$g, d$t, c("pre_shared", "trend_post"))
        if (is.na(ig_t) || is.na(il_t)) next
        p_idx <- c(p_idx, ig_t, il_t)
        coeff <- c(coeff, 1.0, -1.0)
      }
      if (d$from != d$ref) {
        ig_f <- cross_pair_find_param(params, d$g, d$from, c("pre_target", "anticipation_target"))
        il_f <- cross_pair_find_param(params, d$g, d$from, c("pre_shared", "trend_post"))
        if (is.na(ig_f) || is.na(il_f)) next
        p_idx <- c(p_idx, ig_f, il_f)
        coeff <- c(coeff, -1.0, 1.0)
      }
      if (length(p_idx) > 0L) {
        specs[[length(specs) + 1L]] <- gmm_target_spec(tid, p_idx, coeff)
        meta_rows[[length(meta_rows) + 1L]] <- data.frame(
          target_id = tid, kind = "pre", ledger_index = NA_integer_, pre_index = j,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(specs)) {
    p1 <- params$param_index[1L]
    specs[[1L]] <- gmm_target_spec(params$param_id[1L], p1, 1.0)
    meta_rows[[1L]] <- data.frame(
      target_id = params$param_id[1L], kind = "internal_unavailable",
      ledger_index = NA_integer_, pre_index = NA_integer_, stringsAsFactors = FALSE
    )
  }

  list(specs = specs, table = do.call(rbind, meta_rows))
}

run_chained_cross_pair_fit <- function(
  data, idname, tname, yname, gname,
  control, xformla = ~1, cluster_vec = NULL,
  target_cohorts = NULL, anticipation = 0L,
  control_group = "notyettreated", base_period = "universal",
  weighting = "optimal", pt_assumption = "post", pt_start = NULL,
  est_method = "dr", cband = TRUE, alp = 0.05, biters = 1000L,
  btype = "rademacher", seed = NULL, fixed_signs = NULL,
  cp04 = NULL, frame_id = NULL, frame_Z = NULL,
  frame_weights = NULL, frame_g = NULL,
  allow_frame_only = FALSE, na.rm = FALSE,
  bstrap = TRUE, caller_rng = TRUE, gmm_covariance = "tikhonov",
  compute_inffunc = TRUE, print_details = FALSE,
  aggregations = character(0L), prepared = NULL,
  reduce_structural_redundancy = TRUE,
  automatic_chunk_rows = if (!is.null(cp04)) isTRUE(cp04$automatic_chunk_rows) else TRUE,
  keep_nuisance = FALSE
) {
  if (!identical(keep_nuisance, FALSE) && !identical(keep_nuisance, TRUE)) {
    cdid_ipw_pipeline_abort("`keep_nuisance` must be TRUE or FALSE.", "cdid_ipw_pipeline_input_error")
  }
  wt_info <- resolve_gmm_weighting_choice(
    weighting, arg = "weighting",
    error_fn = cdid_ipw_pipeline_abort,
    error_class = "cdid_ipw_pipeline_input_error"
  )
  wt_mode <- wt_info$canonical
  wt_requested <- wt_info$requested
  cg_mode <- match.arg(control_group, c("nevertreated", "notyettreated"))
  bp_mode <- match.arg(base_period, c("varying", "universal"))
  boot_flags <- ipw_pipeline_resolve_bootstrap_flags(bstrap, cband, where = "cross-pair pipeline inference")
  ctrl <- validate_ipw_pipeline_control(control, weighting = wt_mode)

  if (!is.null(gmm_covariance)) {
    cov_method <- resolve_gmm_covariance_choice(
      gmm_covariance, arg = "gmm_covariance",
      error_fn = cdid_ipw_pipeline_abort,
      error_class = "cdid_ipw_pipeline_input_error"
    )
    if (!is.null(ctrl$gmm_accuracy$covariance_strategy)) {
      ctrl$gmm_accuracy$covariance_strategy$method <- cov_method
      ctrl$gmm_accuracy$covariance_strategy$weighting_requested <- wt_requested
      ctrl$gmm_accuracy$covariance_strategy$weighting_canonical <- wt_mode
    } else if (identical(ctrl$gmm_accuracy$method, "conventional_chol_qr")) {
      ctrl$gmm_accuracy$covariance_strategy <- gmm_covariance_strategy(
        method = cov_method,
        weighting_requested = wt_requested,
        weighting_canonical = wt_mode
      )
    }
  } else if (!is.null(ctrl$gmm_accuracy$covariance_strategy)) {
    ctrl$gmm_accuracy$covariance_strategy$weighting_requested <- wt_requested
    ctrl$gmm_accuracy$covariance_strategy$weighting_canonical <- wt_mode
  }

  bound <- chained_request_prepared_bind(
    prepared, data, idname, tname, yname, gname,
    frame_id, frame_Z, frame_weights, frame_g,
    allow_frame_only, NULL, NULL, cp04, NULL
  )
  panel_owner <- if (is.null(bound)) {
    construct_ipw_panel_owner(
      data = data, idname = idname, tname = tname, yname = yname, gname = gname,
      control = ctrl$owner_control
    )
  } else {
    construct_ipw_panel_owner_prepared(
      data = data, idname = idname, tname = tname, yname = yname, gname = gname,
      control = ctrl$owner_control, prepared = bound
    )
  }

  f_id <- frame_id
  f_Z <- as.matrix(frame_Z)
  f_w <- if (!is.null(frame_weights)) as.double(frame_weights) else rep(1.0, length(f_id))
  f_g <- frame_g
  frame_owner <- if (is.null(bound)) {
    construct_ipw_frame_owner(
      panel_owner = panel_owner, id = f_id, Z = f_Z, weights = f_w,
      control = ctrl$owner_control, frame_g = f_g,
      allow_frame_only = allow_frame_only, design_on_A_only = !is.null(cp04)
    )
  } else {
    construct_ipw_frame_owner_prepared(
      panel_owner = panel_owner, id = frame_id, Z = frame_Z, weights = frame_weights,
      control = ctrl$owner_control, frame_g = frame_g,
      prepared = bound, design_on_A_only = !is.null(cp04)
    )
  }

  use_general_cp04 <- !is.null(cp04) && (
    !is.null(cp04$reference_source) ||
    identical(cp04$outcome_pooling, "pooled") ||
    !is.null(cp04$selection_est_method) ||
    identical(cp04$canonical_assumption, "conditional_mtar")
  )
  cp04_baseline <- !is.null(cp04) && (
    identical(cp04$or_method, "baseline") ||
    identical(cp04$dr_method, "baseline") ||
    identical(cp04$ipw_method, "baseline")
  )
  ctrl_cp04_owner <- ctrl$owner_control
  if (use_general_cp04 && isTRUE(automatic_chunk_rows)) {
    d_z_auto <- max(
      1L,
      if (!is.null(cp04$Z)) ncol(cp04$Z) else 1L,
      ncol(f_Z),
      if (!is.null(cp04$reference_source$p_bound)) as.integer(cp04$reference_source$p_bound) else 1L
    )
    n_rec_auto <- max(1L, length(unique(if (!is.null(target_cohorts)) target_cohorts else panel_owner$cohorts)))
    n_upper_auto <- as.integer(min(
      .Machine$integer.max,
      as.double(length(f_id)) * (if (!is.null(cp04$reference_source)) as.double(n_rec_auto) else 1.0)
    ))
    entries_lim_auto <- min(
      ctrl_cp04_owner$max_retained_entries,
      ctrl_cp04_owner$max_matrix_entries,
      ctrl_cp04_owner$max_output_entries
    )
    b_auto <- as.integer(min(
      n_upper_auto,
      max(ctrl_cp04_owner$max_chunk_rows, floor((entries_lim_auto / 4.0) / (128.0 * (d_z_auto + 8.0))))
    ))
    if (is.finite(b_auto) && b_auto > ctrl_cp04_owner$max_chunk_rows) {
      ctrl_cp04_owner$max_chunk_rows <- b_auto
    }
  }
  cp04_reader <- if (!use_general_cp04) NULL else if (is.null(bound)) {
    construct_cp04_frame_reader(
      frame_owner, f_id, cp04$Z, if (cp04_baseline) cp04$frame_g else f_g, f_w,
      ctrl_cp04_owner, U = cp04$U, rho = cp04$rho,
      universe = cp04$universe, justification = cp04$justification,
      allow_unknown_off_A = cp04_baseline, reference_carrier = !is.null(cp04$reference_source)
    )
  } else {
    construct_cp04_frame_reader_prepared(
      frame_owner, frame_id, cp04$Z, cp04$frame_g, frame_weights,
      ctrl_cp04_owner, prepared = bound, U = cp04$U, rho = cp04$rho,
      universe = cp04$universe, justification = cp04$justification,
      allow_unknown_off_A = TRUE, reference_carrier = FALSE
    )
  }
  prepared <- bound <- NULL

  n_total <- as.integer(length(f_id))
  n_a_units <- as.integer(panel_owner$n_a_units)
  a_to_full <- match(panel_owner$a_id, f_id)
  eligible_full <- seq_len(n_total) %in% a_to_full
  full_frame_g <- if (!is.null(f_g)) as.double(f_g) else {
    out_g <- rep(NA_real_, n_total)
    known <- match(f_id, panel_owner$units$id)
    out_g[!is.na(known)] <- as.double(panel_owner$units$g[known[!is.na(known)]])
    out_g
  }
  full_frame_w <- as.double(f_w)
  # Opt-in fitted-nuisance capture (R/chained_cross_pair_nuisance.R); NULL by default.
  nz_sink <- if (keep_nuisance) cross_pair_nuisance_sink(cluster_vec, f_id, a_to_full, ctrl$owner_control) else NULL
  share_law <- NULL
  if (!is.null(cp04) && !is.null(cp04$reference_source)) {
    full_frame_w <- as.double(cp04$reference_source$V)
    share_law <- "first_calendar_masses"
  }
  eligibility <- cp04$reference_source$eligibility
  share_eligibility <- chained_reference_share_eligibility(eligibility, full_frame_g)
  cluster_full <- cluster_vec

  max_coh <- max(panel_owner$cohorts)
  tgt_cohorts <- if (is.null(target_cohorts)) {
    cand_g <- panel_owner$cohorts[is.finite(panel_owner$cohorts) & panel_owner$cohorts < max_coh]
    if (!length(cand_g)) {
      cdid_ipw_pipeline_abort(
        "No treated target cohorts strictly precede the latest control cohort.",
        c("cdid_ipw_pipeline_target_error", "cdid_catalog_target_error")
      )
    }
    sort(unique(as.double(cand_g)))
  } else {
    sort(unique(as.double(target_cohorts)))
  }

  pt_valid <- validate_chained_pt_request(pt_assumption, pt_start, panel_owner$all_times)
  eff_pt_start <- pt_valid$pt_start
  ledger <- ipw_pipeline_request_ledger(tgt_cohorts, panel_owner$all_times, panel_owner$cohorts, anticipation, cg_mode)
  R_full <- nrow(ledger)
  aux_plan <- ipw_auxiliary_spec(ledger, panel_owner$all_times, panel_owner$cohorts, anticipation, cg_mode)
  pre_specs <- aux_plan$selections[[bp_mode]]
  Q_rep <- length(pre_specs)

  agg_types <- if (!is.null(aggregations) && length(aggregations) > 0L) unique(as.character(aggregations)) else character(0L)
  reserve_draws <- boot_flags$draws ||
    (isTRUE(boot_flags$released_false_true) && any(agg_types %in% c("group", "dynamic", "calendar")))
  reduce_structural_redundancy <- isTRUE(reduce_structural_redundancy) && (wt_mode %in% c("optimal", "within"))

  controls_out <- list(
    anticipation = as.integer(anticipation), control_group = cg_mode, base_period = bp_mode,
    pt_assumption = pt_assumption, pt_start = pt_start, cross_pair = TRUE,
    weighting = wt_mode, weighting_requested = wt_requested,
    gmm_covariance = if (!is.null(ctrl$gmm_accuracy$covariance_strategy)) ctrl$gmm_accuracy$covariance_strategy$method else NULL,
    covariance_strategy = ctrl$gmm_accuracy$covariance_strategy,
    reduce_structural_redundancy = reduce_structural_redundancy, run_numerical_transfer = FALSE,
    root_witness_supplied = FALSE, custom_est_method_supplied = FALSE,
    cband = cband, alp = alp, biters = as.integer(biters), btype = btype,
    seed = if (is.null(seed)) NULL else as.integer(seed), fixed_signs_supplied = !is.null(fixed_signs), na_rm = na.rm,
    bstrap = bstrap, caller_rng = caller_rng, bootstrap_requested = boot_flags$bootstrap, flag_law = boot_flags$law,
    draws_reserved = reserve_draws, compute_inffunc = compute_inffunc, print_details = print_details,
    or_basis = cp04$or_basis, outcome_recipe = cp04$outcome_recipe, treatment_columns = cp04$treatment_columns,
    pipeline_control = ctrl
  )

  run_cell_inference <- function(rm) {
    if (!isTRUE(compute_inffunc))
      return(ipw_pipeline_point_only_inference(rm$theta_hat, rm$available, n_total, cluster_full, eligible_full, alp))
    ipw_pipeline_supported_inference(rm$theta_hat, rm$if_full, rm$available, cluster_full, alp, cband, biters, btype, seed,
      ctrl$multiplier_budget, ctrl$inference_accuracy, NULL, fixed_signs, eligible_full,
      reason = "no reported cell is supported by the cross-pair comparison catalog",
      bstrap = bstrap, caller_rng = caller_rng, target_meta = rm$target_meta)
  }
  run_aggregations <- function(rm) {
    out <- list()
    for (agg_t in agg_types) {
      out[[agg_t]] <- aggregate_chained_ipw_targets(
        theta_hat = rm$theta_hat, target_if_full = rm$if_full, target_meta = rm$target_meta,
        full_frame_g = full_frame_g, full_frame_w = full_frame_w, type = agg_t, anticipation = anticipation,
        cluster_full = cluster_full, alp = alp, cband = cband, biters = biters, btype = btype, seed = seed,
        multiplier_budget = ctrl$multiplier_budget, inference_accuracy = ctrl$inference_accuracy,
        enclosure = NULL, fixed_signs = fixed_signs, eligible_full = eligible_full, na.rm = na.rm, available = rm$available,
        bstrap = bstrap, caller_rng = caller_rng, calendar = panel_owner$all_times, cohort_grid = panel_owner$cohorts,
        share_law = share_law, frame_ids = f_id, share_eligibility = share_eligibility)
    }
    out
  }

  pm <- cross_pair_extract_panel_matrices(panel_owner, n_a_units)
  Z_full_a <- frame_owner$Z_a
  cat_res <- build_cross_pair_link_catalog(
    O_mat = pm$O_mat, all_times = pm$all_times, a_g = as.double(frame_owner$a_g),
    w_a = as.double(frame_owner$w_a), tgt_cohorts = tgt_cohorts,
    cohort_universe = panel_owner$cohorts, anticipation = anticipation, control_group = cg_mode,
    eligibility = eligibility, a_to_full = a_to_full
  )
  links_all <- cat_res$links
  M_all <- nrow(links_all)

  if (n_a_units < 2L || M_all < 1L) {
    ledger <- ipw_pipeline_ledger_support(ledger, NULL, NULL, global_status = list(
      status = "empty_catalog",
      reason = "the cross-pair comparison catalog retained zero links for the requested target cohorts and control group"
    ))
    empty_rm <- build_ipw_auxiliary_reporting(ledger, panel_owner$all_times, panel_owner$cohorts, anticipation, cg_mode,
      NULL, frame_owner, NULL, ctrl, NULL, NULL, FALSE, bp_mode, n_total, a_to_full, 1024.0,
      pt_assumption = pt_assumption, pt_start = pt_start, aux_spec = aux_plan, compute_if = compute_inffunc)
    empty_out <- ipw_pipeline_support_only_result("empty_catalog",
      "Cross-pair comparison catalog retained zero links for the requested target cohorts and control group.",
      ledger, panel_owner, frame_owner, NULL, n_total, empty_rm, run_cell_inference(empty_rm), run_aggregations(empty_rm),
      bp_mode, wt_mode, cg_mode, controls_out,
      list(initial_envelope_entries = 1024.0, simultaneous_entries = 1024.0, downstream_work_ops = 1024.0,
           reachable_numeric_payload_only = TRUE, peak_heap_or_RSS_claim = FALSE),
      nuisance_registry = NULL, selection = cp04_selection_record(cp04, NULL, cp04_reader))
    if (!is.null(nz_sink)) {
      nz_sink$route <- list(core = if (use_general_cp04) "general" else "baseline", est_method = est_method,
        n_links = as.integer(M_all), n_a_units = n_a_units)
      empty_out$nuisance_capture_input <- cross_pair_nuisance_input(nz_sink, "empty_catalog", paste0(
        "The cross-pair catalog retained ", M_all, " link(s) on ", n_a_units,
        " pair-eligible unit(s); no nuisance model was fitted."))
    }
    return(empty_out)
  }

  progress <- chained_progress_emitter(print_details)
  if (!is.null(progress)) {
    progress("cross-pair catalog: ", M_all, " candidate moment link(s) across ", length(tgt_cohorts), " target cohort(s)")
  }
  need_if <- isTRUE(compute_inffunc) || identical(wt_mode, "optimal")
  gen_mom <- NULL
  active_xp_backend <- NULL
  if (use_general_cp04) {
    gen_mom <- evaluate_cross_pair_general_moments(
      links_all = links_all, panel_owner = panel_owner, frame_owner = frame_owner,
      cp04 = cp04, cp04_reader = cp04_reader, f_Z = f_Z, pm = pm,
      est_method = est_method, ctrl = ctrl, wt_mode = wt_mode,
      reduce_structural_redundancy = reduce_structural_redundancy,
      compute_if = need_if, nuisance_sink = nz_sink
    )
    active_xp_backend <- gen_mom$backend
    kept_idx <- gen_mom$kept_idx
    dropped_idx <- gen_mom$dropped_idx
    sf_gram <- gen_mom$sf_gram
    links_used <- links_all[kept_idx, , drop = FALSE]
    rownames(links_used) <- NULL
    L_used <- nrow(links_used)

    design_all <- build_cross_pair_chained_design(
      links_all, pm$all_times, panel_owner$cohorts, anticipation, pt_assumption, eff_pt_start, ledger
    )
    design_used <- if (L_used == M_all) {
      design_all
    } else {
      build_cross_pair_chained_design(
        links_used, pm$all_times, panel_owner$cohorts, anticipation, pt_assumption, eff_pt_start, ledger
      )
    }

    mhat_raw_used <- gen_mom$mhat_used
    mhat_used <- stats::setNames(as.double(mhat_raw_used), design_used$incidence$link_ids)
    V_used <- gen_mom$V_used
    dimnames(V_used) <- list(design_used$incidence$link_ids, design_used$incidence$link_ids)
    phi_replay_mat <- if (isTRUE(compute_inffunc)) gen_mom$Phi_core else NULL
    gen_mom$Phi_core <- NULL
    mean_phi_out <- gen_mom$mean_phi
    replay_n <- n_total
    replay_ids <- f_id
    scale_id_str <- "scale:CP04_complete_F_psi_centered_1_over_nA"
    max_exp_diff <- gen_mom$max_moment_diff
    red_stage_str <- if (reduce_structural_redundancy) "post_fit_verified_forest" else "none"
  } else {
    sf_gram <- NULL
    if (reduce_structural_redundancy) {
      sf <- cross_pair_spanning_forest(links_all)
      kept_idx <- sf$kept_idx
      dropped_idx <- sf$dropped_idx
      if (length(kept_idx) < M_all) sf_gram <- sf$G
    } else {
      kept_idx <- seq_len(M_all)
      dropped_idx <- integer(0L)
    }
    links_used <- links_all[kept_idx, , drop = FALSE]
    rownames(links_used) <- NULL
    L_used <- nrow(links_used)

    design_all <- build_cross_pair_chained_design(
      links_all, pm$all_times, panel_owner$cohorts, anticipation, pt_assumption, eff_pt_start, ledger
    )
    design_used <- if (L_used == M_all) {
      design_all
    } else {
      build_cross_pair_chained_design(
        links_used, pm$all_times, panel_owner$cohorts, anticipation, pt_assumption, eff_pt_start, ledger
      )
    }

    mom <- evaluate_cross_pair_moments(
      links_df = links_used, Y_mat = pm$Y_mat,
      a_g = as.double(frame_owner$a_g), w_a = as.double(frame_owner$w_a), Z_a = Z_full_a,
      obs_by_coh = cat_res$obs_by_coh, units_by_coh = cat_res$units_by_coh, coh_vals = cat_res$coh_vals,
      est_method = est_method, cp04 = cp04, glm_control = ctrl$glm_control, wt_mode = wt_mode,
      compute_if = need_if, nuisance_sink = nz_sink
    )
    active_xp_backend <- mom$backend
    mhat_raw_used <- mom$mhat
    mhat_used <- stats::setNames(as.double(mhat_raw_used), design_used$incidence$link_ids)
    V_used <- mom$V
    dimnames(V_used) <- list(design_used$incidence$link_ids, design_used$incidence$link_ids)
    phi_replay_mat <- if (isTRUE(compute_inffunc)) mom$Phi else NULL
    mom$Phi <- NULL
    mean_phi_out <- mom$mean_phi
    replay_n <- n_a_units
    replay_ids <- panel_owner$a_id
    scale_id_str <- "scale:admissible_unit_centered_1_over_nA"
    max_exp_diff <- 0.0
    red_stage_str <- if (reduce_structural_redundancy) "outcome_free_before_point_fits" else "none"
  }

  R_rep <- Q_rep + R_full
  if (isTRUE(automatic_chunk_rows)) {
    b_gmm_auto <- as.integer(min(
      replay_n, 8192L,
      max(16L, floor((ctrl$gmm_bounds$max_matrix_entries / 32.0) / max(1.0, as.double(L_used))))
    ))
    if (is.finite(b_gmm_auto) && b_gmm_auto > ctrl$gmm_bounds$max_batch_rows) {
      ctrl$gmm_bounds$max_batch_rows <- b_gmm_auto
    }
    b_mult_auto <- as.integer(min(
      n_total, 8192L,
      max(16L, floor((ctrl$multiplier_budget$max_if_batch_cells / 2.0) / max(1.0, as.double(R_rep))))
    ))
    if (is.finite(b_mult_auto) && b_mult_auto > ctrl$multiplier_budget$max_if_batch_rows) {
      ctrl$multiplier_budget$max_if_batch_rows <- b_mult_auto
    }
  }

  nuis_id_str <- paste0("nuisance_fit:cross_pair:", est_method)
  red_spec <- if (!is.null(sf_gram)) {
    list(
      mode = "sparse_class_certificate",
      G = sf_gram,
      class_certificate = NULL,
      kept_indices = kept_idx,
      n_links_all = M_all,
      reduction_active = TRUE,
      gram_provenance = "mask_spanning_forest_exact_gram"
    )
  } else {
    list(
      mode = "none",
      class_certificate = NULL,
      kept_indices = kept_idx,
      n_links_all = M_all,
      reduction_active = FALSE
    )
  }

  prov_used <- gmm_covariance_provenance(
    row_level = "unit", n_units = replay_n, n_clusters = replay_n,
    n_positive_contributing_rows = replay_n, centered = TRUE,
    weighted = FALSE, weight_type = "unweighted",
    normalization_factor = 1.0 / as.double(n_a_units),
    scale_id = scale_id_str, nuisance_fit_id = nuis_id_str,
    reduction = red_spec
  )
  if (wt_mode %in% c("optimal", "within") && !is.null(ctrl$gmm_accuracy$covariance_strategy)) {
    prov_used$complete_unit_if_gram <- ipw_seal(
      list(covariance = V_used, n_units = prov_used$n_units, nuisance_fit_id = prov_used$nuisance_fit_id),
      "complete_unit_if_gram"
    )
  }

  batch_cap <- max(1L, min(replay_n, as.integer(ctrl$gmm_bounds$max_batch_rows)))
  make_replay_cb <- function(link_ids_cb) {
    function(consume) {
      phi_mat <- phi_replay_mat
      if (!is.null(phi_mat)) colnames(phi_mat) <- link_ids_cb
      pos <- 1L
      while (pos <= replay_n) {
        end_pos <- min(replay_n, pos + batch_cap - 1L)
        rr <- pos:end_pos
        phi_blk <- if (!is.null(phi_mat)) {
          phi_mat[rr, , drop = FALSE]
        } else {
          matrix(0.0, nrow = length(rr), ncol = length(link_ids_cb), dimnames = list(NULL, link_ids_cb))
        }
        consume(list(
          unit_id = replay_ids[rr],
          phi = phi_blk,
          weight = rep(1.0, length(rr)),
          nuisance_fit_id = nuis_id_str
        ))
        pos <- end_pos + 1L
      }
    }
  }

  if_replay_spec <- gmm_if_replay_spec(
    callback = make_replay_cb(design_used$incidence$link_ids),
    input_row_level = "unit", aggregate_to_clusters = FALSE,
    n_input_rows = replay_n, n_output_rows = replay_n,
    n_positive_contributing_rows = replay_n, centered = TRUE,
    weighted = FALSE, weight_type = "unweighted",
    normalization_factor = 1.0 / as.double(n_a_units),
    scale_id = scale_id_str, nuisance_fit_id = nuis_id_str,
    expected_unit_ids = replay_ids
  )

  # All causal post-treatment and causal anticipation cells (plus any PT-ctrl
  # pre-treatment cells on `>= eff_pt_start`) live directly in `design_used` and
  # are solved in `gmm_fit` under the caller's exact `pt_assumption`.
  tgt_plan <- build_cross_pair_target_specs(design_used, ledger, pre_specs, include_pre = TRUE)

  if (!is.null(progress)) {
    progress("cross-pair GMM: solving ", L_used, " retained moment(s) with weighting = \"", wt_mode, "\"")
  }
  gmm_fit <- solve_chained_gmm_core(
    design = design_used, mhat = mhat_used, V = V_used,
    provenance = prov_used, targets = tgt_plan$specs,
    weighting = wt_requested, if_replay = if_replay_spec,
    accuracy = ctrl$gmm_accuracy, bounds = ctrl$gmm_bounds,
    gauge_rule = "lowest", request_param_cov = FALSE,
    materialize_target_if = isTRUE(compute_inffunc)
  )

  # Evaluate strictly pre-anticipation (`t < g - delta`) diagnostic pretrend
  # contrasts that were not shared in `design_used` (under PT-post, PT-all, or
  # partial-window PT-ctrl) on the same moment vector for D4/D21 reporting.
  aux_gmm_fit <- NULL
  aux_tgt_plan <- NULL
  if (Q_rep > 0L) {
    rem_pre_specs <- Filter(function(d) {
      if (identical(d$status, "reference_zero") || identical(d$status, "unavailable")) return(FALSE)
      if (exact_cmp_sum(d$t, as.double(anticipation), d$g) >= 0L) return(FALSE)
      tid <- paste0("pre(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")")
      !(tid %in% gmm_fit$targets$target_id)
    }, pre_specs)
    if (length(rem_pre_specs) > 0L) {
      design_pre <- build_cross_pair_chained_design(
        links_used, pm$all_times, panel_owner$cohorts, anticipation, "ctrl", min(pm$all_times), ledger = NULL
      )
      mhat_pre <- stats::setNames(as.double(mhat_raw_used), design_pre$incidence$link_ids)
      V_pre <- V_used
      dimnames(V_pre) <- list(design_pre$incidence$link_ids, design_pre$incidence$link_ids)
      prov_pre <- prov_used
      if (wt_mode %in% c("optimal", "within") && !is.null(ctrl$gmm_accuracy$covariance_strategy)) {
        prov_pre$complete_unit_if_gram <- ipw_seal(
          list(covariance = V_pre, n_units = prov_pre$n_units, nuisance_fit_id = prov_pre$nuisance_fit_id),
          "complete_unit_if_gram"
        )
      }
      if_replay_pre <- if_replay_spec
      if_replay_pre$callback <- make_replay_cb(design_pre$incidence$link_ids)
      aux_tgt_plan <- build_cross_pair_target_specs(design_pre, ledger[0L, , drop = FALSE], rem_pre_specs, include_pre = TRUE)
      if (!identical(aux_tgt_plan$table$kind, "internal_unavailable")) {
        aux_gmm_fit <- solve_chained_gmm_core(
          design = design_pre, mhat = mhat_pre, V = V_pre,
          provenance = prov_pre, targets = aux_tgt_plan$specs,
          weighting = wt_requested, if_replay = if_replay_pre,
          accuracy = ctrl$gmm_accuracy, bounds = ctrl$gmm_bounds,
          gauge_rule = "lowest", request_param_cov = FALSE,
          materialize_target_if = isTRUE(compute_inffunc)
        )
      }
    }
  }
  phi_replay_mat <- NULL

  # Map post ledger support & IFs
  est_idx <- which(gmm_fit$targets$estimable)
  Phi_full_main <- if (isTRUE(compute_inffunc)) {
    if (length(est_idx) > 0L) {
      target_if_A <- as.matrix(gmm_fit$if_projection$materialized_target_if)
      colnames(target_if_A) <- gmm_fit$targets$target_id[est_idx]
      if (use_general_cp04) {
        target_if_A * (as.double(n_total) / as.double(n_a_units))
      } else {
        transport_if_to_full_frame(Phi_A = target_if_A, n_total = n_total, a_to_full = a_to_full)$Phi_full
      }
    } else {
      matrix(0.0, nrow = n_total, ncol = 0L)
    }
  } else NULL

  Phi_full_aux <- if (isTRUE(compute_inffunc) && !is.null(aux_gmm_fit)) {
    aux_est_idx <- which(aux_gmm_fit$targets$estimable)
    if (length(aux_est_idx) > 0L) {
      aux_if_A <- as.matrix(aux_gmm_fit$if_projection$materialized_target_if)
      colnames(aux_if_A) <- aux_gmm_fit$targets$target_id[aux_est_idx]
      if (use_general_cp04) {
        aux_if_A * (as.double(n_total) / as.double(n_a_units))
      } else {
        transport_if_to_full_frame(Phi_A = aux_if_A, n_total = n_total, a_to_full = a_to_full)$Phi_full
      }
    } else {
      matrix(0.0, nrow = n_total, ncol = 0L)
    }
  } else NULL

  ledger$cohort_has_rows <- ledger$g %in% links_used$g
  ledger$core_target_row <- match(ledger$cell_id, gmm_fit$targets$target_id)
  ledger$loading_column <- match(ledger$core_target_row, est_idx)
  ledger$available <- !is.na(ledger$loading_column)
  ledger$support_status <- ifelse(
    is.na(ledger$ref), "no_admissible_reference",
    ifelse(!ledger$has_clean_control, "no_clean_control_cohort",
      ifelse(!ledger$cohort_has_rows, "cohort_absent_from_catalog",
        ifelse(ledger$available, "estimable",
          as.character(gmm_fit$targets$status[ledger$core_target_row]))))
  )
  ledger$support_reason <- ifelse(
    ledger$available,
    "cross-pair contrast gamma_g(t) - lambda_g(t) identified in the cross-pair incidence graph",
    ifelse(is.na(ledger$core_target_row),
      paste0("target cohort g=", fmt_design_num(ledger$g), " has no supported cross-pair moments"),
      paste0("cross-pair contrast not estimable (", gmm_fit$targets$support_reason[ledger$core_target_row], ")"))
  )
  ledger$population_claim <- "none: catalog-relative support only; not a population non-identification statement"

  core_targets <- data.frame(
    target_id = gmm_fit$targets$target_id,
    kind = ifelse(gmm_fit$targets$estimable, "level", "level_unavailable"),
    ledger_index = tgt_plan$table$ledger_index,
    core_target_row = seq_len(nrow(gmm_fit$targets)),
    loading_column = match(seq_len(nrow(gmm_fit$targets)), est_idx),
    status = gmm_fit$targets$status,
    support_reason = gmm_fit$targets$support_reason,
    stringsAsFactors = FALSE
  )
  basis <- data.frame(
    basis_id = character(), component_id = integer(), from_param_index = integer(),
    to_param_index = integer(), from_id = character(), to_id = character(),
    g = numeric(), time_from = numeric(), time_to = numeric(),
    core_target_row = integer(), loading_column = integer(),
    available = logical(), support_status = character(), support_reason = character(),
    stringsAsFactors = FALSE
  )

  # Assemble unified reported cells (pre + post) in canonical (g, time) order
  out_if <- if (isTRUE(compute_inffunc)) matrix(NA_real_, nrow = n_total, ncol = R_rep) else NULL
  th <- rep(NA_real_, R_rep)
  avail_rep <- logical(R_rep)
  metas <- vector("list", R_rep)
  meta_row <- function(g, time, role, target_id, ref, base_t, is_post, ledger_index, auxiliary_index, status, reason, loading_column, functional_key) {
    data.frame(
      g = g, time = time, role = role, target_id = target_id, param_id = target_id,
      param_index = NA_integer_, functional_key = functional_key, ref = ref,
      base_t = base_t, base_period_mode = bp_mode, is_post = is_post,
      ledger_index = ledger_index, auxiliary_index = auxiliary_index,
      support_status = status, support_reason = reason,
      loading_column = loading_column, stringsAsFactors = FALSE
    )
  }

  aux_est_idx <- if (!is.null(aux_gmm_fit)) which(aux_gmm_fit$targets$estimable) else integer(0L)

  for (j in seq_len(Q_rep)) {
    d <- pre_specs[[j]]
    if (identical(d$status, "reference_zero")) {
      metas[[j]] <- meta_row(
        d$g, d$t, "reference",
        paste0("reference(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")"),
        d$ref, d$from, FALSE, NA_integer_, NA_integer_,
        "deterministic_reference_zero",
        "universal base-period reference cell: exact zero with zero influence function",
        NA_integer_, NA_character_
      )
      th[j] <- 0.0
      if (isTRUE(compute_inffunc)) out_if[, j] <- 0.0
      avail_rep[j] <- TRUE
      next
    }
    tid <- paste0("pre(g=", fmt_design_num(d$g), ",t=", fmt_design_num(d$t), ")")
    row_m <- match(tid, gmm_fit$targets$target_id)
    if (!is.na(row_m)) {
      col_m <- match(row_m, est_idx)
      if (!is.na(col_m)) {
        metas[[j]] <- meta_row(
          d$g, d$t, "pre", tid, d$ref, d$from, FALSE, NA_integer_, as.integer(d$index),
          "supported_cross_pair_pretrend",
          paste0("cross-pair pretrend contrast ", fmt_design_num(d$from), "->", fmt_design_num(d$t)),
          col_m, "cross_pair_pre"
        )
        th[j] <- gmm_fit$targets$estimate[row_m]
        if (isTRUE(compute_inffunc)) out_if[, j] <- Phi_full_main[, col_m]
        avail_rep[j] <- TRUE
      } else {
        why <- paste0("cross-pair pretrend contrast not estimable (", gmm_fit$targets$support_reason[row_m], ")")
        metas[[j]] <- meta_row(
          d$g, d$t, "pre", tid, d$ref, d$from, FALSE, NA_integer_, as.integer(d$index),
          "unavailable_pretrend", why, NA_integer_, NA_character_
        )
      }
    } else {
      row_a <- if (!is.null(aux_gmm_fit)) match(tid, aux_gmm_fit$targets$target_id) else NA_integer_
      col_a <- if (!is.na(row_a)) match(row_a, aux_est_idx) else NA_integer_
      if (!is.na(col_a)) {
        metas[[j]] <- meta_row(
          d$g, d$t, "pre", tid, d$ref, d$from, FALSE, NA_integer_, as.integer(d$index),
          "supported_cross_pair_pretrend",
          paste0("cross-pair pretrend contrast ", fmt_design_num(d$from), "->", fmt_design_num(d$t)),
          NA_integer_, "cross_pair_pre"
        )
        th[j] <- aux_gmm_fit$targets$estimate[row_a]
        if (isTRUE(compute_inffunc)) out_if[, j] <- Phi_full_aux[, col_a]
        avail_rep[j] <- TRUE
      } else {
        why <- if (identical(d$status, "unavailable")) d$reason else "unsupported_cross_pair_pretrend: target or clean control cohort has no connected path on this pre-treatment date"
        metas[[j]] <- meta_row(
          d$g, d$t, "pre", tid, d$ref, d$from, FALSE, NA_integer_, as.integer(d$index),
          "unavailable_pretrend", why, NA_integer_, NA_character_
        )
      }
    }
  }

  for (h in seq_len(R_full)) {
    j <- Q_rep + h
    col <- ledger$loading_column[h]
    metas[[j]] <- meta_row(
      ledger$g[h], ledger$time[h], "post", ledger$cell_id[h], ledger$ref[h], ledger$ref[h],
      TRUE, as.integer(h), NA_integer_, ledger$support_status[h], ledger$support_reason[h],
      col, NA_character_
    )
    if (!is.na(col)) {
      th[j] <- gmm_fit$targets$estimate[ledger$core_target_row[h]]
      if (isTRUE(compute_inffunc)) out_if[, j] <- Phi_full_main[, col]
      avail_rep[j] <- TRUE
    }
  }

  meta <- do.call(rbind, metas)
  o_rep <- order(meta$g, meta$time)
  meta <- meta[o_rep, , drop = FALSE]
  rownames(meta) <- NULL
  th <- th[o_rep]
  avail_rep <- avail_rep[o_rep]
  names(th) <- meta$target_id
  if (isTRUE(compute_inffunc)) {
    out_if <- out_if[, o_rep, drop = FALSE]
    colnames(out_if) <- meta$target_id
  }
  meta <- chained_pt_reporting_tags(meta, pt_assumption, eff_pt_start, as.double(anticipation), aux_plan$descriptors, design_used$links)
  meta$available <- avail_rep
  meta$reported_column <- seq_len(R_rep)

  S_cov <- if (isTRUE(compute_inffunc)) matrix(NA_real_, R_rep, R_rep, dimnames = list(meta$target_id, meta$target_id)) else NULL
  av_idx <- which(avail_rep)
  if (length(av_idx) && isTRUE(compute_inffunc)) {
    centered_if <- sweep(out_if[, av_idx, drop = FALSE], 2L, colMeans(out_if[, av_idx, drop = FALSE]), "-")
    S_cov[av_idx, av_idx] <- crossprod(centered_if) * as.double(n_a_units) / (as.double(n_total)^2)
  }

  empty_desc_tbl <- data.frame(
    index = integer(), g = numeric(), s = numeric(), t = numeric(),
    K = character(), status = character(), reason = character(), stringsAsFactors = FALSE
  )
  reporting_map <- list(
    base_period = bp_mode, target_meta = meta, theta_hat = th, if_full = out_if,
    available = avail_rep, Sigma_a = S_cov, Sigma_a_enclosure = NULL, enclosure = NULL, L_map = NULL,
    auxiliary = list(
      descriptors = list(), selections = list(varying = list(), universal = list()),
      point_links = list(), descriptor_status = empty_desc_tbl,
      unsupported_counts = list(), fitted_views = list(), nuisance_views = list(),
      new_root_problems = list(), link_enclosures = list(), preflight = NULL,
      post_GMM_rows_added = 0L, n_descriptors = 0L, n_supported = 0L, n_unsupported = 0L,
      scope = "cross-pair pre-trend contrasts evaluated directly on the cross-pair incidence graph"
    )
  )

  cell_inference <- run_cell_inference(reporting_map)
  agg_results <- run_aggregations(reporting_map)
  n_available <- sum(ledger$available)
  whole_request_complete <- n_available == R_full && all(avail_rep)
  point_links_stub <- lapply(seq_len(M_all), function(i) list(link_id = design_all$links$link_id[i], fitted_nuisance_ref = NULL))

  approx_entries <- as.double(L_used) * as.double(L_used) +
    as.double(n_a_units) * as.double(L_used) +
    as.double(n_total) * as.double(R_rep) + 1024.0
  approx_work <- 32.0 * as.double(n_a_units) * as.double(L_used) * as.double(L_used) + 1024.0

  out <- structure(
    list(
      support_status = if (whole_request_complete) "all_requested_cells_supported" else "partial_support",
      support = list(
        ledger = ledger, basis = basis, core_targets = core_targets,
        n_requested = R_full, n_available = n_available, n_unavailable = R_full - n_available,
        whole_request_complete = whole_request_complete,
        unavailable = ledger[!ledger$available, c("cell_id", "g", "time", "support_status", "support_reason"), drop = FALSE],
        population_claim = "none: catalog-relative support only; unavailable cells are not population non-identification claims",
        maps = list(
          ledger_to_core_target_row = ledger$core_target_row,
          ledger_to_loading_column = ledger$loading_column,
          basis_to_loading_column = integer(0L),
          core_target_to_loading_column = core_targets$loading_column
        )
      ),
      resource_accounting = list(
        initial_envelope_entries = approx_entries,
        simultaneous_entries = approx_entries,
        historical_replay_peak_entries = approx_entries,
        downstream_work_ops = approx_work,
        reachable_numeric_payload_only = TRUE,
        peak_heap_or_RSS_claim = FALSE,
        numerical_transfer_requested = FALSE
      ),
      panel_owner = panel_owner,
      frame_owner = frame_owner,
      metadata = NULL,
      design_all = design_all,
      design_used = design_used,
      selection = cp04_selection_record(cp04, active_xp_backend, cp04_reader),
      nuisance_registry = if (!is.null(active_xp_backend)) list(cp04 = active_xp_backend) else NULL,
      point_links_all = point_links_stub,
      point_links_used = point_links_stub[kept_idx],
      structural_certificate = NULL,
      reduction_summary = list(
        weighting = wt_mode,
        weighting_requested = wt_requested,
        reduce_structural_redundancy = reduce_structural_redundancy,
        reduction_active = (L_used < M_all),
        n_links_all = M_all,
        n_links_used = L_used,
        kept_indices = kept_idx,
        dropped_indices = dropped_idx,
        max_moment_expansion_diff = max_exp_diff,
        reduction_stage = red_stage_str
      ),
      covariance_stream = list(
        mhat = mhat_used,
        V = V_used,
        mean_phi = mean_phi_out,
        n_a_units = n_a_units
      ),
      gmm_fit = gmm_fit,
      aux_gmm_fit = aux_gmm_fit,
      covariance_handling = gmm_fit$covariance_handling,
      numerical_transfer = NULL,
      auxiliary_reporting = reporting_map$auxiliary,
      full_frame_transport = if (isTRUE(compute_inffunc)) list(Phi_full = Phi_full_main, n_total = n_total, n_a_units = n_a_units) else NULL,
      reporting = list(
        base_period = bp_mode,
        target_meta = reporting_map$target_meta,
        L_map = NULL,
        estimates = cell_inference$estimates,
        se = cell_inference$se,
        analytical_se = cell_inference$analytical_se,
        numerical_certificate = cell_inference$numerical_certificate,
        z_stat = cell_inference$z_stat,
        p_value = cell_inference$p_value,
        ci_lower = cell_inference$ci_lower,
        ci_upper = cell_inference$ci_upper,
        cband_crit_val = cell_inference$cband_crit_val,
        cband_lower = cell_inference$cband_lower,
        cband_upper = cell_inference$cband_upper,
        vcov = cell_inference$vcov,
        Omega_a = reporting_map$Sigma_a,
        Omega_a_enclosure = NULL,
        if_full = reporting_map$if_full,
        available = reporting_map$available,
        support = cell_inference$support,
        inference = cell_inference,
        enclosure = NULL
      ),
      basis = NULL,
      aggregations = agg_results,
      controls = controls_out
    ),
    class = "cdid_ipw_pipeline_result"
  )
  # Opt-in only: compact capture input for chained_cross_pair_nuisance_capture().
  if (!is.null(nz_sink)) out$nuisance_capture_input <- cross_pair_nuisance_input(nz_sink)
  out
}
