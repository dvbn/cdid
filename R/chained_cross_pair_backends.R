# General estimator, selection-model, reference-input, and within-covariance
# backends for cross-pair and same-pair chained estimation.

cross_pair_gram_solve <- function(fit, intercept, rhs, N) {
  is_vec <- is.null(dim(rhs))
  rhs_mat <- if (is_vec) matrix(as.double(rhs), ncol = 1L) else as.matrix(rhs)
  if (nrow(rhs_mat) == 0L) {
    out <- matrix(0.0, 0L, ncol(rhs_mat))
    return(if (is_vec) numeric(0L) else out)
  }
  sol <- if (isTRUE(intercept)) {
    or_pooled_wls_gram_solve(fit, rhs_mat, normalizer = N)
  } else {
    or_wls_gram_solve(fit, rhs_mat, normalizer = N)
  }
  if (is_vec) as.double(sol) else as.matrix(sol)
}

within_cov_from_level_loadings <- function(U_list, g_vec, scale_by_g = NULL) {
  L <- length(g_vec)
  if (L == 0L) return(matrix(0.0, 0L, 0L))
  g_num <- as.double(g_vec)
  u_g <- unique(g_num)
  cols_by_g <- split(seq_len(L), match(g_num, u_g))
  V <- matrix(0.0, L, L)
  for (gi in seq_along(u_g)) {
    cols_g <- cols_by_g[[gi]]
    sc_g <- if (is.null(scale_by_g)) NULL else if (is.list(scale_by_g)) {
      scale_by_g[[as.character(u_g[gi])]]
    } else {
      scale_by_g
    }
    k_g <- length(cols_g)
    S_g <- matrix(0.0, k_g, k_g)
    for (tau in seq_along(U_list)) {
      U_tau <- U_list[[tau]]
      if (!is.null(U_tau) && nrow(U_tau) > 0L) {
        U_sub <- U_tau[, cols_g, drop = FALSE]
        nz_r <- which(rowSums(U_sub != 0) > 0L)
        if (length(nz_r) > 0L) {
          U_nz <- if (is.null(sc_g)) U_sub[nz_r, , drop = FALSE] else U_sub[nz_r, , drop = FALSE] * sc_g[nz_r]
          S_g <- S_g + crossprod(U_nz)
        }
      }
    }
    V[cols_g, cols_g] <- 0.5 * (S_g + t(S_g))
  }
  V
}

scale_level_loadings_by_g <- function(U_list, g_vec, scale_by_g) {
  L <- length(g_vec)
  if (L == 0L || is.null(scale_by_g)) return(U_list)
  g_num <- as.double(g_vec)
  u_g <- unique(g_num)
  cols_by_g <- split(seq_len(L), match(g_num, u_g))
  for (tau in seq_along(U_list)) {
    U_tau <- U_list[[tau]]
    if (is.null(U_tau) || nrow(U_tau) == 0L) next
    if (!is.list(scale_by_g)) {
      nz_r <- which(rowSums(U_tau != 0) > 0L)
      if (length(nz_r) > 0L) {
        U_tau[nz_r, ] <- U_tau[nz_r, , drop = FALSE] * scale_by_g[nz_r]
        U_list[[tau]] <- U_tau
      }
    } else {
      for (gi in seq_along(u_g)) {
        cols_g <- cols_by_g[[gi]]
        sc_g <- scale_by_g[[as.character(u_g[gi])]]
        U_sub <- U_tau[, cols_g, drop = FALSE]
        nz_r <- which(rowSums(U_sub != 0) > 0L)
        if (length(nz_r) > 0L) {
          U_tau[nz_r, cols_g] <- U_sub[nz_r, , drop = FALSE] * sc_g[nz_r]
        }
      }
      U_list[[tau]] <- U_tau
    }
  }
  U_list
}

build_within_scale_by_g <- function(g_vec, view_rec, reader, frame_owner, cp04, reference) {
  n_a <- as.integer(frame_owner$n_a_units)
  if (n_a == 0L || !length(g_vec)) return(numeric(0L))
  if (is.null(reference)) {
    w_a <- as.double(frame_owner$w_a)
    sw <- sum(w_a)
    return(ifelse(w_a > 0, sqrt(sw / w_a), 0.0))
  }
  st <- reader$state
  orig_of_a <- integer(n_a)
  is_a <- st$a_index > 0L
  orig_of_a[st$a_index[is_a]] <- which(is_a)

  if (chained_reference_restricted(reference)) {
    V_a <- as.double(cp04$reference_source$V[orig_of_a])
    shared_weights <- !identical(reference$weight_kind, "varying") ||
      (identical(reference$fix_weights, "first_period") && all(is.finite(V_a) & V_a > 0))
    groups <- unique(g_vec)
    recipe <- view_rec[match(groups, g_vec)]
    out <- lapply(recipe, function(r) {
      w <- if (shared_weights) V_a else {
        keys <- reference$recipes$weight_index[r] * (reference$N + 1) + orig_of_a
        reference$w[chained_reference_lookup(keys, reference$weight_keys)]
      }
      supported <- which(is.finite(w) & w > 0 & chained_reference_eligible(reference, r, orig_of_a))
      scale <- numeric(n_a)
      # Use the eligible repeat observers' weight mean with the original
      # repeat-observer normalization shared by all target groups.
      if (length(supported))
        scale[supported] <- sqrt(n_a * mean(w[supported]) / w[supported])
      scale
    })
    return(stats::setNames(out, as.character(groups)))
  }

  if (!identical(reference$weight_kind, "varying")) {
    w_a <- as.double(cp04$reference_source$V[orig_of_a])
    sw <- sum(w_a)
    return(ifelse(w_a > 0, sqrt(sw / w_a), 0.0))
  }

  if (identical(reference$fix_weights, "first_period") && !is.null(cp04$reference_source$V)) {
    V_a <- as.double(cp04$reference_source$V[orig_of_a])
    if (all(is.finite(V_a) & V_a > 0)) {
      sw <- sum(V_a)
      return(sqrt(sw / V_a))
    }
  }

  w_indices <- unique(reference$recipes$weight_index[view_rec])
  rel_by_widx <- list()
  rel_fallback <- numeric(n_a)
  N1 <- reference$N + 1
  for (w_idx in w_indices) {
    keys_a <- w_idx * N1 + orig_of_a
    hit <- chained_reference_lookup(keys_a, reference$weight_keys)
    ok <- !is.na(hit)
    w_sup <- numeric(n_a)
    if (any(ok)) w_sup[ok] <- reference$w[hit[ok]]
    mean_w <- if (any(ok)) mean(w_sup[ok]) else 1.0
    rel_w <- ifelse(w_sup > 0, w_sup / mean_w, 0.0)
    rel_by_widx[[as.character(w_idx)]] <- rel_w
    fill_fb <- (rel_fallback <= 0) & (rel_w > 0)
    if (any(fill_fb)) rel_fallback[fill_fb] <- rel_w[fill_fb]
  }

  g_num <- as.double(g_vec)
  u_g <- unique(g_num)
  out_list <- list()
  for (gv in u_g) {
    ell1 <- which(g_num == gv)[1L]
    w_idx_g <- reference$recipes$weight_index[view_rec[ell1]]
    rel_g <- rel_by_widx[[as.character(w_idx_g)]]
    miss <- (rel_g <= 0) & (rel_fallback > 0)
    if (any(miss)) rel_g[miss] <- rel_fallback[miss]
    out_list[[as.character(gv)]] <- ifelse(rel_g > 0, sqrt(as.double(n_a) / rel_g), 0.0)
  }
  out_list
}

# Post-fit exact-mask spanning-forest cycle reducer.
# Nuisance parameters are already fitted on the FULL unreduced catalog; a
# candidate chord edge (s_m, t_m) is pruned if and only if its post-fit moment,
# full-roster IF (when supplied), and per-date level loadings (when supplied)
# telescope along the tree path within `tol`.
# `Phi_block(r0, r1, cols)` and `U_log` supply the same IF columns and
# per-date loadings in bounded row blocks (`block_rows`) instead of dense
# matrices; the checked maxima, and hence the reduction, are identical.
verify_post_fit_forest_reduction <- function(s_vec, t_vec, group_key, mhat_all,
                                             Phi_all = NULL, U_list_all = NULL,
                                             tol = 1e-9, Phi_block = NULL, n_phi_rows = 0L,
                                             U_log = NULL, block_rows = NULL) {
  M <- length(s_vec)
  if (M == 0L) {
    return(list(kept_idx = integer(0L), dropped_idx = integer(0L),
                G = matrix(0.0, 0L, 0L), max_moment_diff = 0.0))
  }
  if (M == 1L) {
    return(list(kept_idx = 1L, dropped_idx = integer(0L),
                G = matrix(1.0, 1L, 1L), max_moment_diff = 0.0))
  }

  grp_int <- match(group_key, unique(group_key))
  n_grp <- max(grp_int)
  grp_indices <- split(seq_len(M), grp_int)

  kept_chunks <- vector("list", n_grp)
  G_blocks <- vector("list", n_grp)
  max_m_diff <- 0.0
  pending <- vector("list", n_grp)
  streamed <- !is.null(Phi_block) || !is.null(U_log)
  b_rows <- if (is.null(block_rows)) .Machine$integer.max else max(1L, as.integer(block_rows))

  for (g_i in seq_len(n_grp)) {
    idx <- grp_indices[[g_i]]
    m_g <- length(idx)
    if (m_g == 1L) {
      kept_chunks[[g_i]] <- idx
      G_blocks[[g_i]] <- matrix(1.0, 1L, 1L)
      next
    }
    sv <- as.double(s_vec[idx])
    tv <- as.double(t_vec[idx])
    ord_local <- order(abs(tv - sv), sv, tv, idx)
    idx_ord <- idx[ord_local]
    sv_ord <- sv[ord_local]
    tv_ord <- tv[ord_local]

    u_nodes <- sort(unique(c(sv_ord, tv_ord)))
    n_nodes <- length(u_nodes)
    u_i <- match(sv_ord, u_nodes)
    v_i <- match(tv_ord, u_nodes)

    uf_parent <- seq_len(n_nodes)
    uf_find <- function(x) {
      while (uf_parent[x] != x) {
        uf_parent[x] <<- uf_parent[uf_parent[x]]
        x <- uf_parent[x]
      }
      x
    }

    is_tree <- logical(m_g)
    for (e in seq_len(m_g)) {
      ra <- uf_find(u_i[e])
      rb <- uf_find(v_i[e])
      if (ra != rb) {
        uf_parent[rb] <- ra
        is_tree[e] <- TRUE
      }
    }

    tree_pos <- which(is_tree)
    cand_pos <- which(!is_tree)
    k_g <- length(tree_pos)
    d_g <- length(cand_pos)
    tree_ids <- idx_ord[tree_pos]

    if (d_g == 0L) {
      kept_chunks[[g_i]] <- tree_ids
      G_blocks[[g_i]] <- diag(1.0, k_g)
      next
    }

    adj_to <- vector("list", n_nodes)
    adj_le <- vector("list", n_nodes)
    adj_sgn <- vector("list", n_nodes)
    for (le in seq_len(k_g)) {
      a <- u_i[tree_pos[le]]
      b <- v_i[tree_pos[le]]
      adj_to[[a]] <- c(adj_to[[a]], b)
      adj_le[[a]] <- c(adj_le[[a]], le)
      adj_sgn[[a]] <- c(adj_sgn[[a]], 1.0)
      adj_to[[b]] <- c(adj_to[[b]], a)
      adj_le[[b]] <- c(adj_le[[b]], le)
      adj_sgn[[b]] <- c(adj_sgn[[b]], -1.0)
    }

    visited <- logical(n_nodes)
    depth <- integer(n_nodes)
    par_v <- integer(n_nodes)
    par_le <- integer(n_nodes)
    par_sgn <- numeric(n_nodes)
    queue <- integer(n_nodes)

    for (r0 in seq_len(n_nodes)) {
      if (visited[r0] || is.null(adj_to[[r0]])) next
      visited[r0] <- TRUE
      depth[r0] <- 0L
      qh <- 1L
      qt <- 1L
      queue[1L] <- r0
      while (qh <= qt) {
        cur <- queue[qh]
        qh <- qh + 1L
        nbs <- adj_to[[cur]]
        les <- adj_le[[cur]]
        sgs <- adj_sgn[[cur]]
        for (ii in seq_along(nbs)) {
          nxt <- nbs[ii]
          if (!visited[nxt]) {
            visited[nxt] <- TRUE
            depth[nxt] <- depth[cur] + 1L
            par_v[nxt] <- cur
            par_le[nxt] <- les[ii]
            par_sgn[nxt] <- sgs[ii]
            qt <- qt + 1L
            queue[qt] <- nxt
          }
        }
      }
    }

    E_cand <- matrix(0.0, nrow = k_g, ncol = d_g)
    for (dd in seq_len(d_g)) {
      ep <- cand_pos[dd]
      a <- u_i[ep]
      b <- v_i[ep]
      while (depth[a] > depth[b]) {
        E_cand[par_le[a], dd] <- E_cand[par_le[a], dd] - par_sgn[a]
        a <- par_v[a]
      }
      while (depth[b] > depth[a]) {
        E_cand[par_le[b], dd] <- E_cand[par_le[b], dd] + par_sgn[b]
        b <- par_v[b]
      }
      while (a != b) {
        E_cand[par_le[a], dd] <- E_cand[par_le[a], dd] - par_sgn[a]
        a <- par_v[a]
        E_cand[par_le[b], dd] <- E_cand[par_le[b], dd] + par_sgn[b]
        b <- par_v[b]
      }
    }

    cand_ids <- idx_ord[cand_pos]
    m_tree <- mhat_all[tree_ids]
    m_cand <- mhat_all[cand_ids]
    m_pred <- as.vector(crossprod(E_cand, m_tree))
    max_tree_abs <- apply(abs(E_cand) * abs(m_tree), 2L, max)
    m_scale <- pmax(1.0, abs(m_cand), max_tree_abs)
    m_err <- abs(m_cand - m_pred)
    ok_vec <- is.finite(m_err) & (m_err <= tol * m_scale)

    if (any(ok_vec) && !is.null(Phi_block)) {
      surv <- which(ok_vec)
      E_s <- E_cand[, surv, drop = FALSE]
      d_s <- length(surv)
      mx_ca <- numeric(d_s)
      mx_pr <- numeric(d_s)
      mx_err <- numeric(d_s)
      any_nz <- FALSE
      n_ph <- as.integer(n_phi_rows)
      for (r0 in seq.int(1L, max(1L, n_ph), by = min(b_rows, max(1L, n_ph)))) {
        r1 <- min(n_ph, r0 + b_rows - 1L)
        if (r1 < r0) next
        # Release the previous block's buffers before forming the next block.
        Phi_tr <- Phi_ca <- Phi_pred <- Phi_ca_nz <- NULL
        Phi_tr <- Phi_block(r0, r1, tree_ids)
        Phi_ca <- Phi_block(r0, r1, cand_ids[surv])
        nz_phi <- which(rowSums(Phi_tr != 0) > 0L | rowSums(Phi_ca != 0) > 0L)
        if (!length(nz_phi)) next
        any_nz <- TRUE
        Phi_pred <- Phi_tr[nz_phi, , drop = FALSE] %*% E_s
        Phi_ca_nz <- Phi_ca[nz_phi, , drop = FALSE]
        mx_ca <- pmax(mx_ca, apply(abs(Phi_ca_nz), 2L, max))
        mx_pr <- pmax(mx_pr, apply(abs(Phi_pred), 2L, max))
        mx_err <- pmax(mx_err, apply(abs(Phi_ca_nz - Phi_pred), 2L, max))
      }
      if (any_nz) {
        phi_scale <- pmax(1.0, mx_ca, mx_pr)
        ok_vec[surv] <- is.finite(mx_err) & (mx_err <= tol * phi_scale)
      }
    } else if (any(ok_vec) && !is.null(Phi_all)) {
      surv <- which(ok_vec)
      Phi_tr <- Phi_all[, tree_ids, drop = FALSE]
      Phi_ca <- Phi_all[, cand_ids[surv], drop = FALSE]
      nz_phi <- which(rowSums(Phi_tr != 0) > 0L | rowSums(Phi_ca != 0) > 0L)
      if (length(nz_phi) > 0L) {
        Phi_pred <- Phi_tr[nz_phi, , drop = FALSE] %*% E_cand[, surv, drop = FALSE]
        Phi_ca_nz <- Phi_ca[nz_phi, , drop = FALSE]
        phi_scale <- pmax(1.0, apply(abs(Phi_ca_nz), 2L, max), apply(abs(Phi_pred), 2L, max))
        phi_err <- apply(abs(Phi_ca_nz - Phi_pred), 2L, max)
        ok_vec[surv] <- is.finite(phi_err) & (phi_err <= tol * phi_scale)
      }
    }

    if (any(ok_vec) && !is.null(U_log)) {
      surv <- which(ok_vec)
      pending[[g_i]] <- list(tree = tree_ids, cand = cand_ids[surv], surv = surv,
        E = E_cand[, surv, drop = FALSE])
    } else if (any(ok_vec) && !is.null(U_list_all)) {
      for (tau in seq_along(U_list_all)) {
        surv <- which(ok_vec)
        if (!length(surv)) break
        U_tau <- U_list_all[[tau]]
        U_tr <- U_tau[, tree_ids, drop = FALSE]
        U_ca <- U_tau[, cand_ids[surv], drop = FALSE]
        nz_u <- which(rowSums(U_tr != 0) > 0L | rowSums(U_ca != 0) > 0L)
        if (!length(nz_u)) next
        U_pred <- U_tr[nz_u, , drop = FALSE] %*% E_cand[, surv, drop = FALSE]
        U_ca_nz <- U_ca[nz_u, , drop = FALSE]
        has_any <- (colSums(U_ca_nz != 0) > 0L) | (colSums(U_pred != 0) > 0L)
        if (any(has_any)) {
          u_scale <- pmax(1.0, apply(abs(U_ca_nz), 2L, max), apply(abs(U_pred), 2L, max))
          u_err <- apply(abs(U_ca_nz - U_pred), 2L, max)
          bad_u <- has_any & (!is.finite(u_err) | u_err > tol * u_scale)
          if (any(bad_u)) ok_vec[surv[bad_u]] <- FALSE
        }
      }
    }

    if (streamed) {
      if (is.null(pending[[g_i]])) pending[[g_i]] <- list()
      pending[[g_i]]$final <- list(ok_vec = ok_vec, m_err = m_err, cand_ids = cand_ids,
        tree_ids = tree_ids, k_g = k_g, E_cand = E_cand)
      next
    }
    fin <- forest_group_finalize(ok_vec, m_err, cand_ids, tree_ids, k_g, E_cand)
    if (fin$max_ok_err > max_m_diff) max_m_diff <- fin$max_ok_err
    kept_chunks[[g_i]] <- fin$kept_g
    G_blocks[[g_i]] <- fin$B_g
  }

  if (streamed) {
    chk_ids <- which(vapply(pending, function(pg) !is.null(pg$cand), logical(1L)))
    if (length(chk_ids)) {
      bad <- loading_log_witness(U_log, lapply(pending[chk_ids], function(pg)
        list(tree = pg$tree, cand = pg$cand, E = pg$E)), tol,
        min(b_rows, max(1L, U_log$n_a)))
      for (jj in seq_along(chk_ids)) {
        pg <- pending[[chk_ids[jj]]]
        if (any(bad[[jj]])) pending[[chk_ids[jj]]]$final$ok_vec[pg$surv[bad[[jj]]]] <- FALSE
      }
    }
    for (g_i in seq_len(n_grp)) {
      fn <- pending[[g_i]]$final
      if (is.null(fn)) next
      fin <- forest_group_finalize(fn$ok_vec, fn$m_err, fn$cand_ids, fn$tree_ids, fn$k_g, fn$E_cand)
      if (fin$max_ok_err > max_m_diff) max_m_diff <- fin$max_ok_err
      kept_chunks[[g_i]] <- fin$kept_g
      G_blocks[[g_i]] <- fin$B_g
    }
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
  kept_idx <- as.integer(kept_raw[ord_k])
  dropped_idx <- as.integer(setdiff(seq_len(M), kept_idx))
  G_mat <- G_mat[ord_k, ord_k, drop = FALSE]
  G_mat <- 0.5 * (G_mat + t(G_mat))

  list(kept_idx = kept_idx, dropped_idx = dropped_idx,
       G = G_mat, max_moment_diff = max_m_diff)
}

forest_group_finalize <- function(ok_vec, m_err, cand_ids, tree_ids, k_g, E_cand) {
  ok_pos <- which(ok_vec)
  max_ok_err <- if (length(ok_pos) > 0L) max(m_err[ok_pos]) else 0.0
  extra_keep <- cand_ids[!ok_vec]
  kept_g <- c(tree_ids, extra_keep)
  k_tot <- length(kept_g)
  B_g <- diag(1.0, k_tot)
  if (length(ok_pos) > 0L) {
    B_tree <- diag(1.0, k_g) + tcrossprod(E_cand[, ok_pos, drop = FALSE])
    B_g[seq_len(k_g), seq_len(k_g)] <- 0.5 * (B_tree + t(B_tree))
  }
  list(kept_g = kept_g, B_g = B_g, max_ok_err = max_ok_err)
}

# Build canonical link descriptors and lightweight cdid_cross_pair_cell_view
# objects from either cross-pair `links_df` or same-pair `design_links`.
build_general_link_descriptors <- function(links_df, all_times, reader, O_mat,
                                           sel_mode, is_cross_pair) {
  st <- reader$state
  L <- nrow(links_df)
  n_cal <- length(all_times)
  G_coh <- length(st$cohorts)
  pair_cnt_by_ci <- vector("list", G_coh)
  if (is.null(links_df$n_observed)) for (ci_s in seq_len(G_coh)) {
    r_f <- st$cohort_order[seq.int(st$cohort_start[ci_s], st$cohort_end[ci_s])]
    a_ix <- st$a_index[r_f]
    a_ix <- a_ix[a_ix > 0L]
    if (length(a_ix) > 0L) {
      O_ci <- matrix(as.double(O_mat[a_ix, , drop = FALSE]), nrow = length(a_ix), ncol = n_cal)
      pair_cnt_by_ci[[ci_s]] <- crossprod(O_ci)
    } else {
      pair_cnt_by_ci[[ci_s]] <- matrix(0.0, n_cal, n_cal)
    }
  }

  descs <- vector("list", L)
  views <- vector("list", L)
  dl_pools <- if (!is_cross_pair) design_pool_cache(links_df$control_pool) else NULL
  cs_all <- match(as.double(links_df$s), all_times)
  ct_all <- match(as.double(links_df$t), all_times)

  for (ell in seq_len(L)) {
    g_val <- as.double(links_df$g[ell])
    s_val <- as.double(links_df$s[ell])
    t_val <- as.double(links_df$t[ell])
    cs <- cs_all[ell]
    ct <- ct_all[ell]
    pair_obj <- list(s = s_val, t = t_val, s_calendar_index = cs, t_calendar_index = ct)

    if (is_cross_pair) {
      is_tgt <- identical(links_df$link_role[ell], "target")
      has_tgt_out <- is_tgt
      has_ctrl <- !is_tgt
      tgt_on_pair <- is_tgt && sel_mode %in% c("none", "ipw")
      K_vals <- if (is_tgt) numeric(0L) else as.double(links_df$cohort[ell])
      ctrl_sign <- if (is_tgt) 0.0 else 1.0
      sel_cohs <- if (is_tgt) g_val else K_vals
      fk <- as.character(links_df$functional_key[ell])
      lid <- paste0("xplink_", ell)
    } else {
      has_tgt_out <- TRUE
      has_ctrl <- TRUE
      tgt_on_pair <- sel_mode %in% c("none", "ipw")
      K_vals <- sort(as.double(dl_pools$values(ell)), method = "radix")
      ctrl_sign <- -1.0
      sel_cohs <- c(g_val, K_vals)
      fk <- as.character(links_df$functional_key[ell])
      lid <- if (!is.null(links_df$link_id)) as.character(links_df$link_id[ell]) else paste0("splink_", ell)
    }
    if (!startsWith(fk, "cmpcat:")) fk <- paste0("cmpcat:", fk)

    sel_cis <- match(sel_cohs, st$cohorts)
    n_obs_risk <- if (is.null(links_df$n_observed)) {
      sum(vapply(sel_cis, function(ci_s) pair_cnt_by_ci[[ci_s]][cs, ct], numeric(1)))
    } else links_df$n_observed[ell]

    # Per-cohort requirement modes for `chained_reference_requirements`
    n_k <- length(K_vals)
    if (sel_mode == "none") {
      coh_x <- if (!has_ctrl) "pair" else c(if (tgt_on_pair) "pair" else "A", rep("pair", n_k))
      coh_s <- rep("none", 1L + n_k)
      coh_w <- rep("A", 1L + n_k)
    } else if (sel_mode == "ipw") {
      coh_x <- if (!has_ctrl) "pair" else c(if (tgt_on_pair) "pair" else "U", rep("pair", n_k))
      coh_s <- if (!has_ctrl) "U" else c(if (tgt_on_pair) "U" else "none", rep("U", n_k))
      coh_w <- rep("U", 1L + n_k)
    } else {
      coh_x <- rep("F", 1L + n_k)
      coh_s <- if (!has_ctrl) "U" else c(if (tgt_on_pair) "U" else "none", rep("U", n_k))
      coh_w <- rep("F", 1L + n_k)
    }

    views[[ell]] <- structure(
      list(
        link_id = lid,
        functional_key = fk,
        pair = pair_obj,
        cohorts = list(g = g_val, K = K_vals, cohorts = sel_cohs),
        support = list(
          n_observed_risk_units = as.integer(max(1L, n_obs_risk)),
          by_cohort_pair_count = rep(1L, length(sel_cohs))
        ),
        cohort_x_rows = coh_x,
        cohort_s_rows = coh_s,
        cohort_w_rows = coh_w
      ),
      class = "cdid_cross_pair_cell_view"
    )

    descs[[ell]] <- list(
      ell = ell,
      link_id = lid,
      functional_key = fk,
      g = g_val,
      ci_g = match(g_val, st$cohorts),
      K = K_vals,
      ci_K = match(K_vals, st$cohorts),
      s = s_val,
      t = t_val,
      cs = cs,
      ct = ct,
      L_steps = as.integer(ct - cs),
      pair = pair_obj,
      has_target_outcome = has_tgt_out,
      has_control = has_ctrl,
      target_on_pair = tgt_on_pair,
      ctrl_sign = ctrl_sign
    )
  }

  list(descs = descs, views = views)
}

# Accumulate unscaled pair-difference Jacobian loadings into per-date A-unit matrices U_list.
add_pair_diff_to_U <- function(U_list, a_idx, cs, ct, w_mat, cols = NULL) {
  if (!length(a_idx)) return(U_list)
  if (is.matrix(w_mat)) {
    if (ncol(w_mat) == 0L) return(U_list)
    if (is.null(cols)) {
      U_list[[ct]][a_idx, ] <- U_list[[ct]][a_idx, , drop = FALSE] + w_mat
      U_list[[cs]][a_idx, ] <- U_list[[cs]][a_idx, , drop = FALSE] - w_mat
    } else {
      U_list[[ct]][a_idx, cols] <- U_list[[ct]][a_idx, cols, drop = FALSE] + w_mat
      U_list[[cs]][a_idx, cols] <- U_list[[cs]][a_idx, cols, drop = FALSE] - w_mat
    }
  } else {
    col_j <- if (is.null(cols)) 1L else as.integer(cols[1L])
    w_vec <- as.double(w_mat)
    U_list[[ct]][a_idx, col_j] <- U_list[[ct]][a_idx, col_j] + w_vec
    U_list[[cs]][a_idx, col_j] <- U_list[[cs]][a_idx, col_j] - w_vec
  }
  U_list
}

# Core unified evaluator for Regime A (`selection = "none"`) and Regime B
# (`selection = "conditional_mtar", selection_est_method = "ipw"`).
evaluate_weighted_regimes_core <- function(descs, views, reader, frame_owner, f_Z,
                                           Y_mat, O_mat, all_times, cp04, reference,
                                           est_method, outcome_pooling, sel_mode,
                                           glm_control, rank_tol, ctrl,
                                           compute_if = TRUE, compute_within = FALSE,
                                           existing_sel_fit = NULL, nuisance_sink = NULL,
                                           within_log_limit = Inf, center_if = FALSE) {
  st <- reader$state
  N <- st$N
  n_a <- nrow(Y_mat)
  n_cal <- length(all_times)
  L <- length(descs)
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  is_or <- identical(est_method, "reg")
  is_dr <- identical(est_method, "dr")
  need_or <- is_or || is_dr
  need_trt <- !is_or

  view_rec <- if (is.null(reference)) rep.int(1L, L) else as.integer(reference$view_recipe)

  # Precompute cohort row subsets on F, U, A
  G_coh <- length(st$cohorts)
  rows_F_by_ci <- vector("list", G_coh)
  rows_U_by_ci <- vector("list", G_coh)
  rows_A_by_ci <- vector("list", G_coh)
  for (ci in seq_len(G_coh)) {
    rf <- st$cohort_order[seq.int(st$cohort_start[ci], st$cohort_end[ci])]
    rows_F_by_ci[[ci]] <- rf
    rows_U_by_ci[[ci]] <- rf[st$U[rf]]
    rows_A_by_ci[[ci]] <- rf[st$a_index[rf] > 0L]
  }

  eligible_rows <- function(rows, recipe) {
    if (!chained_reference_restricted(reference)) return(rows)
    rows[chained_reference_eligible(reference, recipe, rows)]
  }

  # Helper to read unit weights and designs under static or reference route
  get_weights <- function(r_id, rows) {
    if (is.null(reference)) st$w[rows] else chained_reference_weight(reference, r_id, rows)
  }
  X_src <- if (!is.null(cp04$X_U)) cp04$X_U else f_Z
  recipe_or <- if (need_or) {
    if (!is.null(reference) && !is.null(reference$outcome_recipe)) reference$outcome_recipe else cp04$outcome_recipe
  } else NULL

  get_X_trt <- function(r_id, rows) {
    if (!is.null(reference)) return(chained_reference_x(reference, r_id, rows, "x"))
    if (is_dr && !is.null(cp04$treatment_columns)) X_src[rows, cp04$treatment_columns, drop = FALSE] else X_src[rows, , drop = FALSE]
  }
  get_X_or <- function(r_id, rows) {
    if (!is.null(reference)) {
      if (!is.null(reference$or_index)) return(chained_reference_design(reference, r_id, rows))
      return(chained_reference_x(reference, r_id, rows, "x")[, recipe_or$columns, drop = FALSE])
    }
    X_src[rows, recipe_or$columns, drop = FALSE]
  }

  # Upstream selection / rate setup
  sel_fit <- NULL
  up_map <- NULL
  up_if <- NULL
  sel_op <- NULL

  if (sel_mode == "ipw") {
    sel_fit <- if (!is.null(existing_sel_fit)) existing_sel_fit else {
      fit_selection_models(
        reader = reader, link_views = views,
        selection_model = if (is.null(cp04$selection_model)) "pair_logit" else cp04$selection_model,
        selection_pooling = if (is.null(cp04$selection_pooling)) "separate" else cp04$selection_pooling,
        control = ctrl$owner_control, glm_control = glm_ctrl,
        reference = reference
      )
    }
    all_ci <- sort(unique(unlist(lapply(descs, function(d) c(d$ci_g, d$ci_K)))))
    extra_kr <- if (!is.null(reference)) {
      do.call(rbind, lapply(seq_len(L), function(ell) {
        d <- descs[[ell]]
        data.frame(k = c(d$ci_g, d$ci_K), r = view_rec[ell])
      }))
    } else NULL
    up_map <- build_selection_upstream_map(sel_fit, cohorts_a = all_ci, extra_rates_kr = extra_kr)
    if (compute_if && up_map$D > 0L) {
      up_if <- matrix(0.0, N, up_map$D)
      uids_all <- seq_len(up_map$D)
      B_ch <- min(N, ctrl$owner_control$max_chunk_rows)
      op_up <- selection_prediction_operation(sel_fit, source_reader = reader, source_reference = reference)
      on.exit(if (!is.null(op_up)) op_up$abort(), add = TRUE)
      for (s0 in seq.int(1L, N, by = B_ch)) {
        rr <- seq.int(s0, min(N, s0 + B_ch - 1L))
        up_if[rr, ] <- op_up$upstream(up_map, rr, uids_all, FALSE)$up
      }
      op_up$finish(); op_up <- NULL
    }
    sel_op <- selection_prediction_operation(sel_fit, source_reader = reader, source_reference = reference)
    on.exit(if (!is.null(sel_op)) sel_op$abort(), add = TRUE)
  }

  rate_coord_of <- function(ci, r_id) {
    if (is.null(reference)) match(ci, up_map$cohorts_a) else {
      which(up_map$cohorts_a == ci & up_map$rates_ref == r_id)[1L]
    }
  }
  a_val_of <- function(ci, r_id) {
    if (is.null(reference)) st$a[ci] else chained_reference_stat(reference, r_id, ci, "a")
  }

  # In Regime A (`selection = "none"`), deduplicate baseline pair-rate cells (ci, cs, ct, r_id)
  rate_env <- new.env(hash = TRUE, parent = emptyenv())
  rate_list <- list()
  get_baseline_cell <- function(ci, cs, ct, r_id) {
    key <- paste(ci, cs, ct, r_id, sep = ":")
    rec <- rate_env[[key]]
    if (!is.null(rec)) return(rec)
    rA <- eligible_rows(rows_A_by_ci[[ci]], r_id)
    a_ix <- st$a_index[rA]
    is_obs <- O_mat[a_ix, cs] & O_mat[a_ix, ct]
    rows_obs <- rA[is_obs]
    a_obs <- a_ix[is_obs]
    wA <- get_weights(r_id, rA)
    w_obs <- wA[is_obs]
    M_k <- sum(wA)
    if (!length(rows_obs) || !is.finite(M_k) || M_k <= 0 || sum(w_obs) <= 0) {
      cp04_abort("Unsupported baseline cohort or pair support.", "support")
    }
    is_sub <- length(rows_obs) < length(rA)
    r_raw <- if (is_sub) sum(w_obs) / M_k else 1.0
    r_c <- max(r_raw, ipw_clip_eps)
    r_act <- is_sub && (r_raw > ipw_clip_eps)
    id_r <- length(rate_list) + 1L
    rec <- list(
      rate_id = id_r, ci = ci, cs = cs, ct = ct, r_id = r_id,
      rows_A = rA, rows_obs = rows_obs, a_obs = a_obs,
      w_A = wA, w_obs = w_obs, fw = w_obs / r_c,
      M_k = M_k, is_sub = is_sub, r_raw = r_raw, r_c = r_c, r_act = r_act,
      dy = Y_mat[a_obs, ct] - Y_mat[a_obs, cs]
    )
    rate_list[[id_r]] <<- rec
    rate_env[[key]] <- rec
    rec
  }

  # In Regime B (`selection_est_method = "ipw"`), deduplicate selection pair cells (q_id)
  sel_cell_cache <- if (sel_mode == "ipw") vector("list", length(sel_fit$pair_cohorts)) else NULL
  get_sel_cell <- function(ell, pos_in_view, cs, ct, r_id) {
    q_id <- if (!is.null(existing_sel_fit)) {
      ci <- if (pos_in_view == 1L && descs[[ell]]$target_on_pair) descs[[ell]]$ci_g else
        descs[[ell]]$ci_K[pos_in_view - as.integer(descs[[ell]]$target_on_pair)]
      sp <- if (is.null(reference)) list(cohort = ci, pair = descs[[ell]]$pair) else
        list(cohort = ci, pair = descs[[ell]]$pair, reference = r_id)
      resolve_pair_cohort_id(sel_fit, sp)
    } else {
      sel_fit$links[[ell]]$pair_cohort_ids[pos_in_view]
    }
    rec <- sel_cell_cache[[q_id]]
    if (!is.null(rec)) return(rec)
    pc <- sel_fit$pair_cohorts[[q_id]]
    ci <- pc$cohort
    rA <- eligible_rows(rows_A_by_ci[[ci]], r_id)
    a_ix <- st$a_index[rA]
    is_obs <- O_mat[a_ix, cs] & O_mat[a_ix, ct]
    rows_obs <- rA[is_obs]
    a_obs <- a_ix[is_obs]
    pred <- sel_op$pair(rows_obs, q_id)
    ld <- if (compute_if) sel_op$deriv(up_map, q_id, rows_obs, pred = pred)$log_deriv else NULL
    rec <- list(
      q_id = q_id, ci = ci, cs = cs, ct = ct, r_id = r_id,
      rows_obs = rows_obs, a_obs = a_obs,
      prob = pred$prob, fw = pred$weight,
      ch = up_map$pair_channels[[q_id]], ld = ld,
      dy = Y_mat[a_obs, ct] - Y_mat[a_obs, cs]
    )
    sel_cell_cache[[q_id]] <<- rec
    rec
  }

  # Target domain cache for off-pair target evaluation (cross-pair control links)
  tgt_dom_env <- new.env(hash = TRUE, parent = emptyenv())
  get_off_pair_target <- function(ci_g, r_id) {
    key <- paste(ci_g, r_id, sep = ":")
    rec <- tgt_dom_env[[key]]
    if (!is.null(rec)) return(rec)
    if (sel_mode == "none") {
      rows_g <- eligible_rows(rows_A_by_ci[[ci_g]], r_id)
      if (!length(rows_g)) cp04_abort("Target cohort has zero positive weight mass.", "support")
      w_g <- get_weights(r_id, rows_g)
      fw_g <- w_g
      ch_g <- integer(0L)
      ld_g <- NULL
    } else {
      rows_g <- eligible_rows(rows_U_by_ci[[ci_g]], r_id)
      if (!length(rows_g)) cp04_abort("Target cohort has zero positive weight mass.", "support")
      w_g <- get_weights(r_id, rows_g)
      a_g_val <- a_val_of(ci_g, r_id)
      if (!is.finite(a_g_val) || a_g_val <= 0) {
        cp04_abort("Target cohort has zero positive repeat-observer mass.", "support")
      }
      fw_g <- w_g * a_g_val / st$rho[rows_g]
      rc <- rate_coord_of(ci_g, r_id)
      ch_g <- as.integer(rc)
      ld_g <- if (compute_if) matrix(1.0 / a_g_val, nrow = length(rows_g), ncol = 1L) else NULL
    }
    if (sum(fw_g) <= 0) cp04_abort("Target cohort has zero positive weight mass.", "support")
    X_or_g <- if (need_or) get_X_or(r_id, rows_g) else NULL
    cmax_g <- if (need_or && recipe_or$p > 0L) {
      apply(abs(X_or_g[, recipe_or$slope_index, drop = FALSE]), 2L, max)
    } else numeric(0L)
    rec <- list(
      ci = ci_g, r_id = r_id, rows = rows_g, fw = fw_g, d_plus = sum(fw_g) / N,
      ch = ch_g, ld = ld_g,
      X_trt = if (need_trt) get_X_trt(r_id, rows_g) else NULL,
      X_or = X_or_g, colmax = cmax_g
    )
    tgt_dom_env[[key]] <- rec
    rec
  }

  # ----------------------------------------------------------------------------
  # Stage 1 (if need_or): Fit counterfactual OR (`or_outcome_fit_blocks`) across
  # all control links using exact `cdid` alias & rank checks!
  # ----------------------------------------------------------------------------
  ctrl_links <- which(vapply(descs, `[[`, logical(1L), "has_control"))
  has_ctrl_any <- length(ctrl_links) > 0L
  or_groups <- NULL
  sup_list <- list()
  sup_env <- new.env(hash = TRUE, parent = emptyenv())
  blk_list <- list()
  blk_env <- new.env(hash = TRUE, parent = emptyenv())
  link_blk_id <- integer(L)
  link_tgt_sup_id <- integer(L)

  if (need_or && has_ctrl_any) {
    add_ctrl_support <- function(c_rec, L_steps) {
      key <- paste("ctrl", c_rec$ci, c_rec$cs, c_rec$ct, c_rec$r_id, sep = ":")
      sid <- sup_env[[key]]
      if (!is.null(sid)) return(sid)
      Z_obs <- get_X_or(c_rec$r_id, c_rec$rows_obs)
      acc <- or_pooled_support_start(recipe_or, L_steps)
      acc <- or_pooled_support_add(acc, recipe_or, Z_obs, c_rec$dy, c_rec$fw,
                                   if (sel_mode == "none") "baseline" else "corollary1")
      sid <- length(sup_list) + 1L
      sup_list[[sid]] <<- list(
        kind = "ctrl", cell = c_rec, Z_obs = Z_obs, acc = acc,
        spec = list(summary = acc$summary, scale = 1.0, colmax = acc$colmax)
      )
      sup_env[[key]] <- sid
      sid
    }
    add_tgt_support <- function(ell) {
      d <- descs[[ell]]
      r_id <- view_rec[ell]
      if (d$target_on_pair) {
        c_rec <- if (sel_mode == "none") get_baseline_cell(d$ci_g, d$cs, d$ct, r_id) else
          get_sel_cell(ell, 1L, d$cs, d$ct, r_id)
        key <- paste("tgtonpair", c_rec$ci, c_rec$cs, c_rec$ct, r_id, sep = ":")
        sid <- sup_env[[key]]
        if (!is.null(sid)) return(sid)
        Z_obs <- get_X_or(r_id, c_rec$rows_obs)
        cmax <- if (recipe_or$p > 0L) apply(abs(Z_obs[, recipe_or$slope_index, drop = FALSE]), 2L, max) else numeric(0L)
        sid <- length(sup_list) + 1L
        sup_list[[sid]] <<- list(
          kind = "tgt_pair", cell = c_rec, Z_obs = Z_obs,
          spec = list(summary = NULL, scale = 1.0, colmax = cmax)
        )
        sup_env[[key]] <- sid
        sid
      } else {
        t_dom <- get_off_pair_target(d$ci_g, r_id)
        key <- paste("tgtoffpair", d$ci_g, r_id, sep = ":")
        sid <- sup_env[[key]]
        if (!is.null(sid)) return(sid)
        sid <- length(sup_list) + 1L
        sup_list[[sid]] <<- list(
          kind = "tgt_dom", dom = t_dom, Z_obs = t_dom$X_or,
          spec = list(summary = NULL, scale = 1.0, colmax = t_dom$colmax)
        )
        sup_env[[key]] <- sid
        sid
      }
    }

    or_links_spec <- vector("list", length(ctrl_links))
    for (idx_c in seq_along(ctrl_links)) {
      ell <- ctrl_links[idx_c]
      d <- descs[[ell]]
      r_id <- view_rec[ell]
      m_k <- length(d$ci_K)
      ctrl_sids <- integer(m_k)
      for (kk in seq_len(m_k)) {
        c_rec <- if (sel_mode == "none") get_baseline_cell(d$ci_K[kk], d$cs, d$ct, r_id) else
          get_sel_cell(ell, kk + as.integer(d$target_on_pair), d$cs, d$ct, r_id)
        ctrl_sids[kk] <- add_ctrl_support(c_rec, d$L_steps)
      }
      bkey <- paste(paste(ctrl_sids, collapse = ","), d$cs, d$ct, r_id, sep = "|")
      bid <- blk_env[[bkey]]
      if (is.null(bid)) {
        bid <- length(blk_list) + 1L
        blk_list[[bid]] <- list(
          name = paste0("or_blk_", bid),
          support_ids = as.integer(ctrl_sids),
          L = as.integer(d$L_steps)
        )
        blk_env[[bkey]] <- bid
      }
      tsid <- add_tgt_support(ell)
      link_blk_id[ell] <- bid
      link_tgt_sup_id[ell] <- tsid
      or_links_spec[[idx_c]] <- list(regression = as.integer(bid), support_id = as.integer(tsid))
    }

    visit_or <- function(sid, visit) {
      Z_s <- sup_list[[sid]]$Z_obs
      if (nrow(Z_s) > 0L) visit(Z_s)
      invisible(NULL)
    }
    or_groups <- or_outcome_fit_blocks(
      recipe = recipe_or,
      supports = lapply(sup_list, `[[`, "spec"),
      blocks = blk_list,
      links = or_links_spec,
      rank_tol = rank_tol,
      N = N,
      visit_rows = visit_or,
      mode = outcome_pooling,
      route = if (sel_mode == "none") "baseline" else "corollary1"
    )
  }

  # Precompute OR block fitted values & residuals on control supports
  n_or_blks <- length(blk_list)
  blk_fits <- vector("list", n_or_blks)
  if (n_or_blks > 0L) {
    for (bid in seq_len(n_or_blks)) {
      bf <- or_outcome_block_fit(or_groups, bid)
      sids <- blk_list[[bid]]$support_ids
      w_tot <- sum(vapply(sids, function(sid) sup_list[[sid]]$acc$weight_mass, numeric(1L)))
      d_c_b <- NULL
      mhat_c_b <- NULL
      if (is_dr) {
        kept_sl <- recipe_or$slope_index[bf$keep]
        X_or_b <- if (length(sids) == 1L) sup_list[[sids[1L]]]$Z_obs else
          do.call(rbind, lapply(sids, function(sid) sup_list[[sid]]$Z_obs))
        d_c_b <- if (bf$p > 0L) blk_list[[bid]]$L * X_or_b[, kept_sl, drop = FALSE] else NULL
        mhat_c_b <- (if (recipe_or$intercept) bf$alpha else 0.0) +
          (if (bf$p > 0L) as.vector(d_c_b %*% bf$beta) else numeric(nrow(X_or_b)))
      }
      blk_fits[[bid]] <- c(bf, list(
        support_ids = sids, weight_sum = w_tot, mass = w_tot / N,
        d_c_all = d_c_b, mhat_c_all = mhat_c_b
      ))
    }
  }

  # ----------------------------------------------------------------------------
  # Stage 2: Link point estimates, direct scores, propensity logits, and
  # sensitivities to OR blocks and upstream channels
  # ----------------------------------------------------------------------------
  mhat <- numeric(L)
  Phi_F <- if (compute_if) matrix(0.0, N, L) else NULL
  U_log <- if (compute_within) level_loading_log(n_cal, n_a, L, within_log_limit) else NULL

  # Accumulators for Stage-1 OR adjoints across all L links
  Lambda_or_alpha <- if (n_or_blks > 0L) matrix(0.0, n_or_blks, L) else NULL
  n_or_grps <- if (!is.null(or_groups)) length(or_groups$groups) else 0L
  Lambda_or_beta <- if (n_or_grps > 0L) {
    lapply(or_groups$groups, function(g) matrix(0.0, length(g$keep), L))
  } else NULL

  # Accumulators for upstream rate / selection channels
  D_up <- if (sel_mode == "ipw") up_map$D else 0L
  Lambda_up <- if (compute_if && D_up > 0L) matrix(0.0, D_up, L) else NULL
  # Baseline rate accumulator (initialized after all baseline cells are registered)
  baseline_rate_loads <- vector("list", L)

  trt_cache <- new.env(hash = TRUE, parent = emptyenv())
  n_trt_fits <- 0L

  for (ell in seq_len(L)) {
    d <- descs[[ell]]
    r_id <- view_rec[ell]

    # ---- Target outcome arm (+T_gl) when has_target_outcome == TRUE ----------
    T_gl <- 0.0
    if (d$has_target_outcome) {
      t_cell <- if (sel_mode == "none") get_baseline_cell(d$ci_g, d$cs, d$ct, r_id) else
        get_sel_cell(ell, 1L, d$cs, d$ct, r_id)
      fw_t <- t_cell$fw
      sum_fwt <- sum(fw_t)
      dp_t <- sum_fwt / N
      omega_t <- fw_t / sum_fwt
      T_gl <- sum(omega_t * t_cell$dy)
      if (compute_if) {
        Phi_F[t_cell$rows_obs, ell] <- Phi_F[t_cell$rows_obs, ell] +
          (fw_t / dp_t) * (t_cell$dy - T_gl)
        if (sel_mode == "ipw" && length(t_cell$ch)) {
          u_vec <- as.vector(crossprod((fw_t * (t_cell$dy - T_gl)) / (N * dp_t), t_cell$ld))
          Lambda_up[t_cell$ch, ell] <- Lambda_up[t_cell$ch, ell] + u_vec
        }
      }
      if (compute_within) {
        loading_log_pair_vec(U_log, t_cell$a_obs, d$ct, d$cs, ell, omega_t)
      }
    }

    if (!d$has_control) {
      mhat[ell] <- T_gl
      next
    }

    # ---- Control link / control arm (ctrl_sign * (...)) ----------------------
    sgn <- d$ctrl_sign
    m_k <- length(d$ci_K)
    c_cells <- vector("list", m_k)
    for (kk in seq_len(m_k)) {
      c_cells[[kk]] <- if (sel_mode == "none") get_baseline_cell(d$ci_K[kk], d$cs, d$ct, r_id) else
        get_sel_cell(ell, kk + as.integer(d$target_on_pair), d$cs, d$ct, r_id)
    }

    # Target covariate domain (on pair for same-pair, or on A_g / U_g for cross-pair)
    if (d$target_on_pair) {
      t_dom_cell <- if (sel_mode == "none") get_baseline_cell(d$ci_g, d$cs, d$ct, r_id) else
        get_sel_cell(ell, 1L, d$cs, d$ct, r_id)
      rows_g <- t_dom_cell$rows_obs
      fw_g <- t_dom_cell$fw
      dp_g <- sum(fw_g) / N
      X_trt_g <- if (need_trt) get_X_trt(r_id, rows_g) else NULL
      X_or_g <- if (need_or) sup_list[[link_tgt_sup_id[ell]]]$Z_obs else NULL
      ch_g <- if (sel_mode == "ipw") t_dom_cell$ch else integer(0L)
      ld_g <- if (sel_mode == "ipw") t_dom_cell$ld else NULL
      rate_g <- if (sel_mode == "none") t_dom_cell else NULL
    } else {
      t_dom <- get_off_pair_target(d$ci_g, r_id)
      rows_g <- t_dom$rows
      fw_g <- t_dom$fw
      dp_g <- t_dom$d_plus
      X_trt_g <- t_dom$X_trt
      X_or_g <- t_dom$X_or
      ch_g <- t_dom$ch
      ld_g <- t_dom$ld
      rate_g <- NULL
    }

    # Counterfactual OR prediction on target (`B_gl`) and controls (`mhat_c`)
    B_gl <- 0.0
    d_bar_plus <- numeric(0L)
    bf <- NULL
    if (need_or) {
      bid <- link_blk_id[ell]
      bf <- blk_fits[[bid]]
      grp_id <- bf$group_id
      kept_sl <- recipe_or$slope_index[bf$keep]
      d_g <- d$L_steps * X_or_g[, kept_sl, drop = FALSE]
      mhat_g <- (if (recipe_or$intercept) bf$alpha else 0.0) +
        (if (bf$p > 0L) as.vector(d_g %*% bf$beta) else numeric(nrow(d_g)))
      B_gl <- sum(fw_g * mhat_g) / sum(fw_g)
      d_bar_plus <- if (bf$p > 0L) colSums(fw_g * d_g) / sum(fw_g) else numeric(0L)
      if (compute_if) {
        Phi_F[rows_g, ell] <- Phi_F[rows_g, ell] + sgn * (fw_g / dp_g) * (mhat_g - B_gl)
        if (sel_mode == "ipw" && d$target_on_pair && length(ch_g)) {
          u_bg <- as.vector(crossprod((fw_g * (mhat_g - B_gl)) / (N * dp_g), ld_g))
          Lambda_up[ch_g, ell] <- Lambda_up[ch_g, ell] + sgn * u_bg
        }
      }
    }

    # Treatment propensity logit (`C_ol` and `B_ol`) when `need_trt`
    C_ol <- 0.0
    B_ol <- 0.0
    d_bar_minus <- numeric(0L)
    if (need_trt) {
      pair_tag <- if (!d$target_on_pair && sel_mode == "none" &&
                      all(!vapply(c_cells, `[[`, logical(1L), "is_sub"))) {
        "complete_A"
      } else {
        paste(d$cs, d$ct, sep = ":")
      }
      tkey <- paste(
        d$ci_g, d$target_on_pair, paste(d$ci_K, collapse = ","),
        pair_tag, r_id, sep = "|"
      )
      tfit <- trt_cache[[tkey]]
      rows_c_all <- unlist(lapply(c_cells, `[[`, "rows_obs"), use.names = FALSE)
      a_c_all <- unlist(lapply(c_cells, `[[`, "a_obs"), use.names = FALSE)
      fw_c_all <- unlist(lapply(c_cells, `[[`, "fw"), use.names = FALSE)
      dy_c_all <- unlist(lapply(c_cells, `[[`, "dy"), use.names = FALSE)
      len_c_by_k <- vapply(c_cells, function(cc) length(cc$rows_obs), integer(1L))
      ends_k <- cumsum(len_c_by_k)
      starts_k <- ends_k - len_c_by_k + 1L

      if (is.null(tfit)) {
        X_trt_c <- get_X_trt(r_id, rows_c_all)
        X_all <- rbind(X_trt_g, X_trt_c)
        y_all <- c(rep(1.0, length(rows_g)), rep(0.0, length(rows_c_all)))
        w_all <- c(fw_g, fw_c_all)
        lfit <- cp04_logit(X_all, y_all, w_all, glm_ctrl, N, d)
        lv <- ipw_logit_values(X_all, lfit$coefficients)
        c_pos <- seq.int(length(rows_g) + 1L, nrow(X_all))
        n_trt_fits <- n_trt_fits + 1L
        tfit <- if (compute_if || !is.null(nuisance_sink)) {
          list(
            coefficients = lfit$coefficients, H_chol = lfit$H_chol, H = lfit$H,
            X_all = X_all, y_all = y_all, w_all = w_all, rows_all = c(rows_g, rows_c_all),
            X_c = X_trt_c, prob_all = lv$prob, odds_c = lv$odds[c_pos], act_c = lv$active[c_pos],
            c_pos = c_pos, g_pos = seq_len(length(rows_g))
          )
        } else {
          list(odds_c = lv$odds[c_pos])
        }
        trt_cache[[tkey]] <- tfit
      }
      if (!is.null(nuisance_sink)) cross_pair_nuisance_weighted_trt(nuisance_sink, environment(), tkey)

      wm <- fw_c_all * tfit$odds_c
      sum_wm <- sum(wm)
      dm <- sum_wm / N
      cp04_finite(dm, "control odds normalizer", TRUE)
      omega_c <- wm / sum_wm
      C_ol <- sum(omega_c * dy_c_all)

      mhat_c_all <- NULL
      d_c_all <- NULL
      if (is_dr) {
        d_c_all <- bf$d_c_all
        mhat_c_all <- bf$mhat_c_all
        B_ol <- sum(omega_c * mhat_c_all)
        d_bar_minus <- if (bf$p > 0L) colSums(omega_c * d_c_all) else numeric(0L)
      }

      eff_c_res <- if (is_dr) (dy_c_all - mhat_c_all) - (C_ol - B_ol) else (dy_c_all - C_ol)

      if (compute_if) {
        # Direct control residual on observed control pairs
        Phi_F[rows_c_all, ell] <- Phi_F[rows_c_all, ell] + sgn * (wm / dm) * eff_c_res

        # Propensity logit beta partial & score backprop
        wm_act <- wm * tfit$act_c
        V_gamma <- colSums(wm_act * eff_c_res * tfit$X_c) / (N * dm)
        h_inv_v <- as.double(cp04_solve(tfit$H_chol, matrix(V_gamma, ncol = 1L)))
        score_all <- tfit$w_all * (tfit$y_all - tfit$prob_all) * as.vector(tfit$X_all %*% h_inv_v)
        Phi_F[tfit$rows_all, ell] <- Phi_F[tfit$rows_all, ell] + sgn * score_all

        if (sel_mode == "ipw") {
          # Selection / repeat-rate forcing through propensity score AND control moment
          if (length(ch_g)) {
            sc_g <- tfit$w_all[tfit$g_pos] * (1.0 - tfit$prob_all[tfit$g_pos]) *
              as.vector(tfit$X_all[tfit$g_pos, , drop = FALSE] %*% h_inv_v)
            Lambda_up[ch_g, ell] <- Lambda_up[ch_g, ell] + sgn * as.vector(crossprod(sc_g, ld_g)) / N
          }
          for (kk in seq_len(m_k)) {
            cc <- c_cells[[kk]]
            if (!length(cc$ch)) next
            ix_k <- seq.int(starts_k[kk], ends_k[kk])
            sc_k <- tfit$w_all[tfit$c_pos[ix_k]] * (-tfit$prob_all[tfit$c_pos[ix_k]]) *
              as.vector(tfit$X_c[ix_k, , drop = FALSE] %*% h_inv_v)
            mom_k <- (wm[ix_k] * eff_c_res[ix_k]) / dm
            Lambda_up[cc$ch, ell] <- Lambda_up[cc$ch, ell] +
              sgn * as.vector(crossprod(sc_k + mom_k, cc$ld)) / N
          }
        } else {
          # Baseline pair-rate forcing (target rate when target_on_pair, plus control rates)
          r_loads <- list()
          if (!is.null(rate_g) && rate_g$r_act) {
            sc_g <- tfit$w_all[tfit$g_pos] * (1.0 - tfit$prob_all[tfit$g_pos]) *
              as.vector(tfit$X_all[tfit$g_pos, , drop = FALSE] %*% h_inv_v)
            r_loads[[length(r_loads) + 1L]] <- list(
              rate_id = rate_g$rate_id,
              val = sgn * (-sum(sc_g) / (N * rate_g$r_c))
            )
          }
          for (kk in seq_len(m_k)) {
            cc <- c_cells[[kk]]
            if (!cc$r_act) next
            ix_k <- seq.int(starts_k[kk], ends_k[kk])
            sc_k <- tfit$w_all[tfit$c_pos[ix_k]] * (-tfit$prob_all[tfit$c_pos[ix_k]]) *
              as.vector(tfit$X_c[ix_k, , drop = FALSE] %*% h_inv_v)
            mom_k <- (wm[ix_k] * eff_c_res[ix_k]) / dm
            r_loads[[length(r_loads) + 1L]] <- list(
              rate_id = cc$rate_id,
              val = sgn * (-(sum(sc_k) + sum(mom_k)) / (N * cc$r_c))
            )
          }
          baseline_rate_loads[[ell]] <- r_loads
        }
      }

      if (compute_within) {
        loading_log_pair_vec(U_log, a_c_all, d$ct, d$cs, ell, sgn * omega_c)
      }
    }

    # Record link estimate & OR sensitivities
    if (is_or) {
      mhat[ell] <- T_gl + sgn * B_gl
      if (recipe_or$intercept) Lambda_or_alpha[link_blk_id[ell], ell] <- Lambda_or_alpha[link_blk_id[ell], ell] + sgn
      if (bf$p > 0L) Lambda_or_beta[[bf$group_id]][, ell] <- Lambda_or_beta[[bf$group_id]][, ell] + sgn * d_bar_plus
    } else if (is_dr) {
      mhat[ell] <- T_gl + sgn * (B_gl + C_ol - B_ol)
      if (!recipe_or$intercept) {
        # No intercept to cancel; B_gl - B_ol only has slope term
      }
      if (bf$p > 0L) {
        Lambda_or_beta[[bf$group_id]][, ell] <- Lambda_or_beta[[bf$group_id]][, ell] + sgn * (d_bar_plus - d_bar_minus)
      }
    } else {
      mhat[ell] <- T_gl + sgn * C_ol
    }
  }

  # ----------------------------------------------------------------------------
  # Stage 3: Backpropagate counterfactual OR (`Lambda_or_alpha`, `Lambda_or_beta`)
  # to observed control pairs (for both `Phi_F` and `U_list`) AND upstream rates/logits
  # ----------------------------------------------------------------------------
  n_rates <- length(rate_list)
  Lam_r <- if (compute_if && sel_mode == "none" && n_rates > 0L) {
    mat_r <- matrix(0.0, n_rates, L)
    for (ell in seq_len(L)) {
      rls <- baseline_rate_loads[[ell]]
      for (rl in rls) mat_r[rl$rate_id, ell] <- mat_r[rl$rate_id, ell] + rl$val
    }
    mat_r
  } else NULL

  if (n_or_blks > 0L) {
    # Fold alpha_b into group slope rhs: alpha_b = ybar_b - mean_d_b' beta_grp
    for (bid in seq_len(n_or_blks)) {
      bf <- blk_fits[[bid]]
      if (recipe_or$intercept && bf$p > 0L) {
        lam_a <- Lambda_or_alpha[bid, ]
        if (any(lam_a != 0)) {
          Lambda_or_beta[[bf$group_id]] <- Lambda_or_beta[[bf$group_id]] -
            outer(as.double(bf$mean_d), lam_a)
        }
      }
    }

    nu_or_beta <- vector("list", n_or_grps)
    cols_or_beta <- vector("list", n_or_grps)
    for (gid in seq_len(n_or_grps)) {
      g_obj <- or_groups$groups[[gid]]
      Lam_b <- Lambda_or_beta[[gid]]
      p_g <- nrow(Lam_b)
      cb <- if (p_g > 0L) which(colSums(Lam_b != 0) > 0L) else integer(0L)
      if (length(cb) > 0L) {
        nu_raw <- cross_pair_gram_solve(g_obj$fit, recipe_or$intercept, Lam_b[, cb, drop = FALSE], N)
        keep_c <- which(colSums(nu_raw != 0) > 0L)
        cols_or_beta[[gid]] <- cb[keep_c]
        nu_or_beta[[gid]] <- nu_raw[, keep_c, drop = FALSE]
      } else {
        cols_or_beta[[gid]] <- integer(0L)
        nu_or_beta[[gid]] <- matrix(0.0, nrow = p_g, ncol = 0L)
      }
    }

    blks_by_gid <- split(seq_len(n_or_blks), vapply(blk_fits, `[[`, integer(1L), "group_id"))
    for (gid_chr in names(blks_by_gid)) {
      gid <- as.integer(gid_chr)
      bids_g <- blks_by_gid[[gid_chr]]
      nu_sub <- nu_or_beta[[gid]]
      cols_b <- cols_or_beta[[gid]]
      p_g <- nrow(nu_sub)
      has_slope <- length(cols_b) > 0L
      all_sids_g <- unique(unlist(lapply(bids_g, function(b) blk_fits[[b]]$support_ids), use.names = FALSE))

      if (length(bids_g) == 1L && length(all_sids_g) == 1L) {
        bid <- bids_g[1L]
        bf <- blk_fits[[bid]]
        lam_a <- if (recipe_or$intercept) Lambda_or_alpha[bid, ] else numeric(L)
        cols_a <- which(lam_a != 0)
        has_int <- length(cols_a) > 0L
        if (!has_int && !has_slope) next
        sup <- sup_list[[bf$support_ids[1L]]]
        cc <- sup$cell
        fw_s <- cc$fw
        Z_s <- sup$Z_obs
        kept_sl <- recipe_or$slope_index[bf$keep]
        d_s <- if (p_g > 0L) bf$L * Z_s[, kept_sl, drop = FALSE] else NULL
        mhat_s <- (if (recipe_or$intercept) bf$alpha else 0.0) +
          (if (p_g > 0L) as.vector(d_s %*% bf$beta) else numeric(length(fw_s)))
        e_s <- cc$dy - mhat_s

        if (has_int) {
          W_int <- outer(fw_s / bf$weight_sum, lam_a[cols_a])
          if (compute_if) {
            Phi_F[cc$rows_obs, cols_a] <- Phi_F[cc$rows_obs, cols_a, drop = FALSE] + (N * e_s) * W_int
            if (sel_mode == "ipw" && length(cc$ch)) {
              Lambda_up[cc$ch, cols_a] <- Lambda_up[cc$ch, cols_a, drop = FALSE] + crossprod(cc$ld, e_s * W_int)
            } else if (sel_mode == "none" && cc$r_act) {
              Lam_r[cc$rate_id, cols_a] <- Lam_r[cc$rate_id, cols_a] - colSums(e_s * W_int) / cc$r_c
            }
          }
          if (compute_within) {
            loading_log_pair_outer(U_log, cc$a_obs, cc$ct, cc$cs, cols_a, fw_s / bf$weight_sum, lam_a[cols_a])
          }
        }

        if (has_slope) {
          d_tilde <- if (recipe_or$intercept) sweep(d_s, 2L, bf$mean_d, "-") else d_s
          H_nu <- ((fw_s / N) * d_tilde) %*% nu_sub
          if (compute_if) {
            Phi_F[cc$rows_obs, cols_b] <- Phi_F[cc$rows_obs, cols_b, drop = FALSE] + (N * e_s) * H_nu
            if (sel_mode == "ipw" && D_up > 0L && length(cc$ch)) {
              Lambda_up[cc$ch, cols_b] <- Lambda_up[cc$ch, cols_b, drop = FALSE] + crossprod(cc$ld, e_s * H_nu)
            } else if (sel_mode == "none" && n_rates > 0L && cc$r_act) {
              Lam_r[cc$rate_id, cols_b] <- Lam_r[cc$rate_id, cols_b] - colSums(e_s * H_nu) / cc$r_c
            }
          }
          if (compute_within) {
            loading_log_lowrank(U_log, cc$a_obs, cc$ct, cc$cs, cols_b, (fw_s / N) * d_tilde, nu_sub)
          }
        }
        next
      }

      u_rows <- if (has_slope && compute_if) {
        sort(unique(unlist(lapply(all_sids_g, function(s) sup_list[[s]]$cell$rows_obs), use.names = FALSE)))
      } else integer(0L)
      u_ch <- if (has_slope && compute_if && sel_mode == "ipw" && D_up > 0L) {
        sort(unique(unlist(lapply(all_sids_g, function(s) sup_list[[s]]$cell$ch), use.names = FALSE)))
      } else integer(0L)
      u_rates <- if (has_slope && compute_if && sel_mode == "none" && n_rates > 0L) {
        sort(unique(unlist(lapply(all_sids_g, function(s) {
          if (sup_list[[s]]$cell$r_act) sup_list[[s]]$cell$rate_id else integer(0L)
        }), use.names = FALSE)))
      } else integer(0L)
      u_a_rows <- if (has_slope && compute_within) {
        sort(unique(unlist(lapply(all_sids_g, function(s) sup_list[[s]]$cell$a_obs), use.names = FALSE)))
      } else integer(0L)
      active_taus <- if (has_slope && compute_within) {
        sort(unique(unlist(lapply(all_sids_g, function(s) c(sup_list[[s]]$cell$cs, sup_list[[s]]$cell$ct)), use.names = FALSE)))
      } else integer(0L)

      S_phi_g <- if (length(u_rows) > 0L) matrix(0.0, length(u_rows), p_g) else NULL
      L_up_g <- if (length(u_ch) > 0L) matrix(0.0, length(u_ch), p_g) else NULL
      R_rate_g <- if (length(u_rates) > 0L) matrix(0.0, length(u_rates), p_g) else NULL
      U_grp <- if (length(active_taus) > 0L && length(u_a_rows) > 0L) {
        replicate(length(active_taus), matrix(0.0, length(u_a_rows), p_g), simplify = FALSE)
      } else NULL

      for (bid in bids_g) {
        bf <- blk_fits[[bid]]
        lam_a <- if (recipe_or$intercept) Lambda_or_alpha[bid, ] else numeric(L)
        cols_a <- which(lam_a != 0)
        has_int <- length(cols_a) > 0L
        if (!has_int && !has_slope) next
        lam_a_sub <- if (has_int) lam_a[cols_a] else NULL
        kept_sl <- recipe_or$slope_index[bf$keep]

        for (sid in bf$support_ids) {
          sup <- sup_list[[sid]]
          cc <- sup$cell
          fw_s <- cc$fw
          Z_s <- sup$Z_obs
          d_s <- if (p_g > 0L) bf$L * Z_s[, kept_sl, drop = FALSE] else NULL
          mhat_s <- (if (recipe_or$intercept) bf$alpha else 0.0) +
            (if (p_g > 0L) as.vector(d_s %*% bf$beta) else numeric(length(fw_s)))
          e_s <- cc$dy - mhat_s

          if (has_int) {
            W_int <- outer(fw_s / bf$weight_sum, lam_a_sub)
            if (compute_if) {
              Phi_F[cc$rows_obs, cols_a] <- Phi_F[cc$rows_obs, cols_a, drop = FALSE] + (N * e_s) * W_int
              if (sel_mode == "ipw" && length(cc$ch)) {
                Lambda_up[cc$ch, cols_a] <- Lambda_up[cc$ch, cols_a, drop = FALSE] + crossprod(cc$ld, e_s * W_int)
              } else if (sel_mode == "none" && cc$r_act) {
                Lam_r[cc$rate_id, cols_a] <- Lam_r[cc$rate_id, cols_a] - colSums(e_s * W_int) / cc$r_c
              }
            }
            if (compute_within) {
              loading_log_pair_outer(U_log, cc$a_obs, cc$ct, cc$cs, cols_a, fw_s / bf$weight_sum, lam_a_sub)
            }
          }

          if (has_slope) {
            d_tilde <- if (recipe_or$intercept) sweep(d_s, 2L, bf$mean_d, "-") else d_s
            H_s <- (fw_s / N) * d_tilde
            if (compute_if) {
              if (!is.null(S_phi_g) && length(cc$rows_obs)) {
                r_loc <- match(cc$rows_obs, u_rows)
                S_phi_g[r_loc, ] <- S_phi_g[r_loc, , drop = FALSE] + (N * e_s) * H_s
              }
              if (!is.null(L_up_g) && length(cc$ch)) {
                ch_loc <- match(cc$ch, u_ch)
                L_up_g[ch_loc, ] <- L_up_g[ch_loc, , drop = FALSE] + crossprod(cc$ld, e_s * H_s)
              } else if (!is.null(R_rate_g) && cc$r_act) {
                rt_loc <- match(cc$rate_id, u_rates)
                R_rate_g[rt_loc, ] <- R_rate_g[rt_loc, ] - colSums(e_s * H_s) / cc$r_c
              }
            }
            if (compute_within && !is.null(U_grp) && length(cc$a_obs)) {
              a_loc <- match(cc$a_obs, u_a_rows)
              ts_loc <- match(cc$cs, active_taus)
              tt_loc <- match(cc$ct, active_taus)
              U_grp <- add_pair_diff_to_U(U_grp, a_loc, ts_loc, tt_loc, H_s)
            }
          }
        }
      }

      if (has_slope) {
        if (compute_if) {
          if (!is.null(S_phi_g)) {
            Phi_F[u_rows, cols_b] <- Phi_F[u_rows, cols_b, drop = FALSE] + S_phi_g %*% nu_sub
          }
          if (!is.null(L_up_g) && any(L_up_g != 0)) {
            Lambda_up[u_ch, cols_b] <- Lambda_up[u_ch, cols_b, drop = FALSE] + L_up_g %*% nu_sub
          }
          if (!is.null(R_rate_g) && any(R_rate_g != 0)) {
            Lam_r[u_rates, cols_b] <- Lam_r[u_rates, cols_b, drop = FALSE] + R_rate_g %*% nu_sub
          }
        }
        if (compute_within && !is.null(U_grp)) {
          for (ti in seq_along(active_taus)) {
            if (any(U_grp[[ti]] != 0)) {
              tau <- active_taus[ti]
              loading_log_lowrank(U_log, u_a_rows, tau, 0L, cols_b, U_grp[[ti]], nu_sub)
            }
          }
        }
      }
    }
  }

  # ----------------------------------------------------------------------------
  # Stage 4: Contract upstream rate / selection IFs into Phi_F
  # ----------------------------------------------------------------------------
  if (compute_if) {
    if (sel_mode == "ipw" && D_up > 0L && any(Lambda_up != 0)) {
      cols_nz <- which(colSums(Lambda_up != 0) > 0L)
      if (length(cols_nz)) {
        Phi_F[, cols_nz] <- Phi_F[, cols_nz, drop = FALSE] + up_if %*% Lambda_up[, cols_nz, drop = FALSE]
      }
    } else if (sel_mode == "none" && n_rates > 0L) {
      for (rid in seq_len(n_rates)) {
        lr <- Lam_r[rid, ]
        cols_nz <- which(lr != 0)
        if (!length(cols_nz)) next
        rc <- rate_list[[rid]]
        inv_Mk <- as.double(N) / rc$M_k
        lr_sub <- lr[cols_nz]
        Phi_F[rc$rows_A, cols_nz] <- Phi_F[rc$rows_A, cols_nz, drop = FALSE] -
          outer((rc$w_A * inv_Mk) * rc$r_raw, lr_sub)
        Phi_F[rc$rows_obs, cols_nz] <- Phi_F[rc$rows_obs, cols_nz, drop = FALSE] +
          outer(rc$w_obs * inv_Mk, lr_sub)
      }
    }
  }

  if (!is.null(sel_op)) {
    sel_op$finish(); sel_op <- NULL
  }
  if (!is.null(nuisance_sink)) cross_pair_nuisance_weighted_collect(nuisance_sink, environment())

  ref_env <- if (!is.null(reference)) {
    e <- new.env(parent = emptyenv()); e$reference_reader <- reference; e
  } else NULL

  acc <- list(
    automatic_chunk_rows = isTRUE(cp04$automatic_chunk_rows),
    n_selection_fits = if (!is.null(sel_fit)) length(sel_fit$blocks) else 0L,
    n_selection_cells = if (!is.null(sel_fit)) length(sel_fit$cells) else 0L,
    n_selection_blocks = if (!is.null(sel_fit)) length(sel_fit$blocks) else 0L,
    n_selection_parameters = if (!is.null(sel_fit)) length(sel_fit$coefficients) else 0L,
    n_regression_fits = if (!is.null(or_groups)) length(or_groups$groups) else 0L,
    n_regression_blocks = as.integer(n_or_blks),
    n_supports = length(sup_list),
    n_treatment_fits = as.integer(n_trt_fits),
    n_rate_fits = length(rate_list),
    n_reference_recipes = if (!is.null(reference)) nrow(reference$recipes) else 0L
  )

  backend <- structure(
    list(
      mode = if (sel_mode == "none") "baseline" else "corollary1",
      est_method = est_method,
      outcome_pooling = outcome_pooling,
      outcome_recipe = recipe_or,
      reader_ref = reader$ref,
      reference_ref = ref_env,
      uses_selection_models = !is.null(sel_fit),
      selection_backend = sel_fit,
      accounting = acc
    ),
    class = "cdid_cross_pair_backend"
  )

  # Centered in place, one column at a time, when requested: entrywise
  # `Phi_F - colMeans(Phi_F)`, with transients of one column.
  Phi_mean <- NULL
  if (isTRUE(center_if) && !is.null(Phi_F)) {
    Phi_mean <- colMeans(Phi_F)
    for (j in seq_len(ncol(Phi_F))) Phi_F[, j] <- Phi_F[, j] - Phi_mean[j]
  }
  list(mhat = mhat, Phi_F = Phi_F, Phi_mean = Phi_mean, U_log = U_log, backend = backend)
}

# Core unified evaluator for Regime C (`selection_est_method = "imputation"`)
# and Regime D (`selection_est_method = "dr"` / augmented).
evaluate_reconstruction_regimes_core <- function(descs, views, reader, frame_owner,
                                                 Y_mat, O_mat, all_times, cp04, reference,
                                                 est_method, outcome_pooling, imputation_pooling,
                                                 augmented, glm_control, rank_tol, ctrl,
                                                 compute_if = TRUE, compute_within = FALSE,
                                                 existing_sel_fit = NULL, nuisance_sink = NULL,
                                                 within_log_limit = Inf, center_if = FALSE) {
  L <- length(descs)
  n_a <- nrow(Y_mat)
  n_cal <- length(all_times)
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  is_or <- identical(est_method, "reg")
  is_dr <- identical(est_method, "dr")
  need_or <- is_or || is_dr
  need_trt <- !is_or

  # Resolve static full-F design or recipe-cohort layout
  lay <- NULL
  if (!is.null(reference)) {
    lay <- chained_reference_layout(reader, reference, views, ctrl$owner_control, selection = augmented)
    st <- lay$state
    X_F <- lay$X
    f_cols <- reference$x_columns
    families <- lay$families
    sel_reader <- if (augmented) lay$reader else NULL
    pair_reader <- function(rows, pair, need_dy = TRUE) cp04_pair_rows(reader, lay$orig_row[rows], pair, need_dy)
  } else {
    des <- imputation_resolve_full_design(reader = reader, frame_owner = frame_owner,
                                          X_frame = cp04$X_frame, control = ctrl$owner_control)
    st <- reader$state
    X_F <- des$X
    f_cols <- des$columns
    families <- NULL
    sel_reader <- reader
    pair_reader <- function(rows, pair, need_dy = TRUE) cp04_pair_rows(reader, rows, pair, need_dy)
  }

  N_orig <- st$N
  N_lay <- length(st$w)
  G_coh <- length(st$cohorts)

  rec_raw <- cp04$outcome_recipe
  trt_cols_raw <- cp04$treatment_columns
  if (!is.null(reference) && is.null(rec_raw)) {
    rr <- chained_reference_outcome_recipe(reference)
    rec_raw <- rr$recipe
    trt_cols_raw <- if (need_trt) rr$columns else NULL
  }
  or_chk <- dr_check_outcome_recipe(rec_raw, f_cols, "baseline")
  recipe <- or_chk$recipe
  or_frame_index <- or_chk$frame_index
  slope_cols_F <- if (recipe$p > 0L) or_frame_index[recipe$slope_index] else integer(0L)
  trt_recipe <- if (need_trt) {
    dr_check_treatment_recipe(trt_cols_raw, f_cols, "baseline")
  } else {
    structure(list(columns = character(0L), index = integer(0L), p = 0L), class = "cdid_dr_treatment_recipe")
  }

  # Map each link's target cohort & active cell cohorts in `st` (static or layout)
  link_ci_g <- integer(L)
  link_ci_K <- vector("list", L)
  sel_view_cohorts <- if (!is.null(lay) && augmented) vector("list", L) else NULL
  for (ell in seq_len(L)) {
    d <- descs[[ell]]
    if (is.null(lay)) {
      link_ci_g[ell] <- d$ci_g
      link_ci_K[[ell]] <- d$ci_K
    } else {
      vc <- lay$view_cohorts[[ell]]
      link_ci_g[ell] <- vc[1L]
      link_ci_K[[ell]] <- if (length(vc) > 1L) vc[-1L] else integer(0L)
      if (augmented) {
        sel_view_cohorts[[ell]] <- if (d$has_target_outcome && !d$has_control) vc[1L] else
          if (!d$has_target_outcome) vc[-1L] else vc
      }
    }
  }

  # Deduplicate Stage-1 actual-change cells (ci, cs, ct)
  cell_env <- new.env(hash = TRUE, parent = emptyenv())
  cell_specs <- list()
  add_s1_cell <- function(ci, d) {
    key <- paste(ci, d$cs, d$ct, sep = ":")
    cid <- cell_env[[key]]
    if (!is.null(cid)) return(cid)
    cid <- length(cell_specs) + 1L
    cell_specs[[cid]] <<- list(
      cell_id = as.integer(cid),
      name = paste0("cell_", cid, "_c", ci, "_", d$cs, "_", d$ct),
      cohort_index = as.integer(ci),
      cohort = st$cohorts[ci],
      pair = d$pair,
      L = as.integer(d$L_steps)
    )
    cell_env[[key]] <- cid
    cid
  }

  link_tgt_cid <- integer(L)
  link_ctrl_cids <- vector("list", L)
  for (ell in seq_len(L)) {
    d <- descs[[ell]]
    if (d$has_target_outcome) link_tgt_cid[ell] <- add_s1_cell(link_ci_g[ell], d)
    if (d$has_control) {
      cids_k <- integer(length(link_ci_K[[ell]]))
      for (kk in seq_along(cids_k)) cids_k[kk] <- add_s1_cell(link_ci_K[[ell]][kk], d)
      link_ctrl_cids[[ell]] <- cids_k
    }
  }

  n_cells <- length(cell_specs)
  grp_keys <- vapply(cell_specs, function(cs) {
    fam <- if (is.null(families)) 1L else families[cs$cohort_index]
    if (identical(imputation_pooling, "separate")) paste0("sep:", cs$cell_id) else
      if (identical(imputation_pooling, "within_cohort")) paste0("wc:", fam, ":", cs$cohort_index) else
        paste0("ac:", fam)
  }, character(1L))
  u_gkeys <- unique(grp_keys)
  group_specs <- lapply(seq_along(u_gkeys), function(gid) {
    cids <- which(grp_keys == u_gkeys[gid])
    list(group_id = as.integer(gid), name = paste0("imp_grp_", gid), cell_ids = as.integer(cids))
  })
  active_cohorts <- sort(unique(c(link_ci_g, unlist(link_ci_K))))

  cat_s1 <- list(
    active_cohorts = as.integer(active_cohorts),
    cell_specs = cell_specs,
    group_specs = group_specs
  )

  # Stage 1: Actual-change WLS on observed pairs
  s1 <- imputation_fit_core(
    st = st, X_F = X_F, pair_reader = pair_reader, cat = cat_s1,
    recipe = recipe, or_frame_index = or_frame_index, trt_recipe = trt_recipe,
    est_method = est_method, imputation_pooling = imputation_pooling,
    glm_control = glm_ctrl, rank_tol = rank_tol, ctrl = ctrl$owner_control,
    outcome_pooling = outcome_pooling, stage1_only = TRUE
  )
  cells <- s1$cells
  imp_groups <- s1$groups
  n_imp_grps <- length(imp_groups)

  # Optional observation-probability model fit for `augmented`
  sel_fit <- NULL
  up_map <- NULL
  Xi_theta <- NULL
  D_theta <- 0L
  if (augmented) {
    sel_fit <- if (!is.null(existing_sel_fit) && identical(existing_sel_fit$reader_ref, sel_reader$ref)) {
      existing_sel_fit
    } else {
      fit_selection_models(
        reader = sel_reader, link_views = views,
        selection_model = cp04$selection_model, selection_pooling = cp04$selection_pooling,
        control = ctrl$owner_control, glm_control = glm_ctrl,
        view_cohorts = sel_view_cohorts, families = families
      )
    }
    up_map <- build_selection_upstream_map(sel_fit, cohorts_a = active_cohorts)
    theta_uids <- if (up_map$D > up_map$n_a) seq.int(up_map$n_a + 1L, up_map$D) else integer(0L)
    D_theta <- length(theta_uids)
    if (compute_if && D_theta > 0L) {
      Xi_theta <- matrix(0.0, N_lay, D_theta)
      B_ch <- min(N_lay, ctrl$owner_control$max_chunk_rows)
      op_up <- selection_prediction_operation(sel_fit, source_reader = sel_reader)
      on.exit(if (!is.null(op_up)) op_up$abort(), add = TRUE)
      for (s0 in seq.int(1L, N_lay, by = B_ch)) {
        rr <- seq.int(s0, min(N_lay, s0 + B_ch - 1L))
        Xi_theta[rr, ] <- op_up$upstream(up_map, rr, theta_uids, FALSE)$up
      }
      op_up$finish(); op_up <- NULL
    }
  }

  # Cohort rows on `st` (`1..N_lay`), synthetic weights, and repeat-rate IFs
  rows_by_ci <- vector("list", G_coh)
  colmax_F <- vector("list", G_coh)
  M_one_coh <- numeric(G_coh)
  M_z_coh <- vector("list", G_coh)
  for (ci in active_cohorts) {
    rf <- st$cohort_order[seq.int(st$cohort_start[ci], st$cohort_end[ci])]
    rows_by_ci[[ci]] <- rf
    w_syn <- st$w[rf] * st$a[ci]
    M_one_coh[ci] <- sum(w_syn) / N_orig
    if (recipe$p > 0L) {
      Z_sl <- X_F[rf, slope_cols_F, drop = FALSE]
      colmax_F[[ci]] <- apply(abs(Z_sl), 2L, max)
      M_z_coh[[ci]] <- as.vector(crossprod(w_syn, Z_sl)) / N_orig
    } else {
      colmax_F[[ci]] <- numeric(0L)
      M_z_coh[[ci]] <- numeric(0L)
    }
  }

  # Materialize per-cell pseudo-outcomes Z_c on F_{ci} (jhat for imputation;
  # jhat + S * pi * (dy - jhat) for augmented)
  cell_data <- vector("list", n_cells)
  sel_op <- if (augmented) selection_prediction_operation(sel_fit, source_reader = sel_reader) else NULL
  on.exit(if (!is.null(sel_op)) sel_op$abort(), add = TRUE)

  for (cid in seq_len(n_cells)) {
    cf <- cells[[cid]]
    gf <- imp_groups[[cf$group_id]]
    ci <- cf$cohort_index
    rf <- rows_by_ci[[ci]]
    w_raw <- st$w[rf]
    w_syn <- w_raw * st$a[ci]
    d_imp <- if (gf$p > 0L) cf$L * X_F[rf, gf$kept_cols_F, drop = FALSE] else matrix(0.0, length(rf), 0L)

    hit <- pair_reader(rf, cf$pair, TRUE)
    obs_local <- which(hit$S == 1)
    rows_obs <- rf[obs_local]
    orig_obs <- if (is.null(lay)) rows_obs else lay$orig_row[rows_obs]
    a_obs <- reader$state$a_index[orig_obs]
    dy_obs <- hit$dy[obs_local]

    if (!augmented) {
      a0_imp <- if (recipe$intercept) cf$alpha else 0.0
      jhat <- if (gf$p > 0L) a0_imp + as.vector(d_imp %*% gf$beta) else rep(a0_imp, length(rf))
      Z_val <- jhat
      one_minus_spi <- rep(1.0, length(rf))
      pi_obs <- numeric(0L)
      T_logit <- NULL
      theta_ch <- integer(0L)
    } else {
      q_id <- resolve_pair_cohort_id(sel_fit, list(cohort = ci, pair = cf$pair))
      cf_q <- cf; cf_q$q_id <- q_id
      aug_r <- augmented_cell_rows(
        sel_fit = sel_fit, up_map = up_map, X_F = X_F, pair_reader = pair_reader,
        cf = cf_q, gf = gf, rows = rf, need_logderiv = compute_if,
        pair_predictor = function(sf, rws, q) sel_op$pair(rws, q),
        log_derivative = function(sf, um, q, rws) {
          pr <- sel_op$pair(rws, q)
          dv <- sel_op$deriv(um, q, rws, pred = pr)
          list(pred = pr, T_logit = dv$T_logit)
        }
      )
      jhat <- aug_r$jhat
      Z_val <- aug_r$Z
      pi_obs <- aug_r$pi
      one_minus_spi <- rep(1.0, length(rf))
      if (length(obs_local)) one_minus_spi[obs_local] <- 1.0 - pi_obs
      T_logit <- aug_r$T_logit
      raw_ch <- up_map$pair_channels[[q_id]]
      theta_ch <- raw_ch[raw_ch > up_map$n_a] - up_map$n_a
    }

    e_imp_obs <- dy_obs - jhat[obs_local]
    d_imp_obs <- d_imp[obs_local, , drop = FALSE]
    d_imp_tilde <- if (recipe$intercept && gf$p > 0L) sweep(d_imp_obs, 2L, cf$mean_d, "-") else d_imp_obs

    cell_data[[cid]] <- list(
      cid = cid, cf = cf, gf = gf, ci = ci, rf = rf,
      w_raw = w_raw, w_syn = w_syn, d_imp = d_imp, jhat = jhat, Z = Z_val,
      obs_local = obs_local, rows_obs = rows_obs, a_obs = a_obs,
      e_imp_obs = e_imp_obs, d_imp_tilde = d_imp_tilde,
      one_minus_spi = one_minus_spi, pi_obs = pi_obs,
      T_logit = T_logit, theta_ch = as.integer(theta_ch)
    )
  }
  if (!is.null(sel_op)) { sel_op$finish(); sel_op <- NULL }

  # ----------------------------------------------------------------------------
  # Stage 2 (when need_or & has_ctrl_any): Counterfactual OR on pseudo-outcomes Z
  # ----------------------------------------------------------------------------
  ctrl_links <- which(vapply(descs, `[[`, logical(1L), "has_control"))
  has_ctrl_any <- length(ctrl_links) > 0L
  or_groups <- NULL
  sup_or <- list()
  sup_or_env <- new.env(hash = TRUE, parent = emptyenv())
  blk_or <- list()
  blk_or_env <- new.env(hash = TRUE, parent = emptyenv())
  link_or_bid <- integer(L)
  link_or_tsid <- integer(L)

  if (need_or && has_ctrl_any) {
    add_or_ctrl_sup <- function(cid) {
      key <- paste0("ctrl:", cid)
      sid <- sup_or_env[[key]]
      if (!is.null(sid)) return(sid)
      cd <- cell_data[[cid]]
      acc <- or_pooled_support_start(recipe, cd$cf$L)
      acc <- or_pooled_support_add(acc, recipe, X_F[cd$rf, or_frame_index, drop = FALSE],
                                   cd$Z, cd$w_syn, "imputation_or")
      sid <- length(sup_or) + 1L
      sup_or[[sid]] <<- list(kind = "ctrl", cid = cid, ci = cd$ci, acc = acc,
                             spec = list(summary = acc$summary, scale = 1.0, colmax = colmax_F[[cd$ci]]))
      sup_or_env[[key]] <- sid
      sid
    }
    add_or_tgt_sup <- function(ci_g) {
      key <- paste0("tgt:", ci_g)
      sid <- sup_or_env[[key]]
      if (!is.null(sid)) return(sid)
      sid <- length(sup_or) + 1L
      sup_or[[sid]] <<- list(kind = "tgt", cid = 0L, ci = ci_g, acc = NULL,
                             spec = list(summary = NULL, scale = 1.0, colmax = colmax_F[[ci_g]]))
      sup_or_env[[key]] <- sid
      sid
    }

    or_links_spec <- vector("list", length(ctrl_links))
    for (idx_c in seq_along(ctrl_links)) {
      ell <- ctrl_links[idx_c]
      d <- descs[[ell]]
      cids_k <- link_ctrl_cids[[ell]]
      sids_k <- vapply(cids_k, add_or_ctrl_sup, integer(1L))
      bkey <- paste(paste(sids_k, collapse = ","), d$L_steps, sep = "|")
      bid <- blk_or_env[[bkey]]
      if (is.null(bid)) {
        bid <- length(blk_or) + 1L
        blk_or[[bid]] <- list(name = paste0("imp_or_blk_", bid),
                              support_ids = as.integer(sids_k), L = as.integer(d$L_steps))
        blk_or_env[[bkey]] <- bid
      }
      tsid <- add_or_tgt_sup(link_ci_g[ell])
      link_or_bid[ell] <- bid
      link_or_tsid[ell] <- tsid
      or_links_spec[[idx_c]] <- list(regression = as.integer(bid), support_id = as.integer(tsid))
    }

    visit_or_F <- function(sid, visit) {
      ci <- sup_or[[sid]]$ci
      rf <- rows_by_ci[[ci]]
      if (length(rf)) visit(X_F[rf, or_frame_index, drop = FALSE])
      invisible(NULL)
    }

    or_groups <- tryCatch(or_outcome_fit_blocks(
      recipe = recipe, supports = lapply(sup_or, `[[`, "spec"),
      blocks = blk_or, links = or_links_spec, rank_tol = rank_tol,
      N = N_orig, visit_rows = visit_or_F, mode = outcome_pooling, route = "imputation_or"
    ), cdid_or_wls_rank_error = function(e) {
      e$message <- paste0(conditionMessage(e), " This computation failure does not establish population ATT nonidentification.")
      e$is_identifiability_claim <- FALSE
      stop(e)
    })
  }

  n_or_blks <- length(blk_or)
  blk_fits <- vector("list", n_or_blks)
  for (bid in seq_len(n_or_blks)) {
    bf <- or_outcome_block_fit(or_groups, bid)
    sids <- blk_or[[bid]]$support_ids
    w_tot <- sum(vapply(sids, function(sid) sup_or[[sid]]$acc$weight_mass, numeric(1L)))
    blk_fits[[bid]] <- c(bf, list(support_ids = sids, weight_sum = w_tot, mass = w_tot / N_orig))
  }

  # ----------------------------------------------------------------------------
  # Stage 3 (when need_trt & has_ctrl_any): Synthetic propensity logits on F
  # ----------------------------------------------------------------------------
  trt_env <- new.env(hash = TRUE, parent = emptyenv())
  trt_list <- list()
  get_syn_trt <- function(ci_g, ci_K_vec, spec_d) {
    all_ci <- c(ci_g, sort(unique(ci_K_vec)))
    key <- paste(ci_g, paste(all_ci[-1L], collapse = ","), sep = "|")
    tid <- trt_env[[key]]
    if (!is.null(tid)) return(trt_list[[tid]])
    risk_rows <- sort(unlist(rows_by_ci[all_ci], use.names = FALSE))
    ci_r <- st$cohort_index[risk_rows]
    y_r <- as.double(ci_r == ci_g)
    w_raw <- st$w[risk_rows]
    w_syn <- w_raw * st$a[ci_r]
    X_trt <- X_F[risk_rows, trt_recipe$index, drop = FALSE]
    lfit <- cp04_logit(X_trt, y_r, w_syn, glm_ctrl, N_orig, spec_d)
    lv <- ipw_logit_values(X_trt, lfit$coefficients)
    pos_by_ci <- lapply(all_ci, function(ci_x) match(rows_by_ci[[ci_x]], risk_rows))

    tid <- length(trt_list) + 1L
    rec <- if (compute_if || !is.null(nuisance_sink)) {
      C_a <- matrix(0.0, trt_recipe$p, length(all_ci))
      raw_sc <- w_raw * (y_r - lv$prob) * X_trt
      for (kk in seq_along(all_ci)) {
        sub_k <- pos_by_ci[[kk]]
        if (length(sub_k)) C_a[, kk] <- colSums(raw_sc[sub_k, , drop = FALSE]) / N_orig
      }
      list(
        tid = tid, all_ci = all_ci, pos_by_ci = pos_by_ci, risk_rows = risk_rows,
        ci_r = ci_r, y_r = y_r, w_syn = w_syn, X_trt = X_trt,
        coefficients = lfit$coefficients, H_chol = lfit$H_chol, H = lfit$H,
        prob = lv$prob, odds = lv$odds, active = lv$active, C_a = C_a
      )
    } else {
      list(tid = tid, all_ci = all_ci, pos_by_ci = pos_by_ci, odds = lv$odds)
    }
    trt_list[[tid]] <<- rec
    trt_env[[key]] <- tid
    rec
  }

  # ----------------------------------------------------------------------------
  # Stage 4: Link estimates, direct full-F scores, and reverse-mode sensitivities
  # ----------------------------------------------------------------------------
  mhat <- numeric(L)
  Phi_lay <- if (compute_if) matrix(0.0, N_lay, L) else NULL
  U_log <- if (compute_within) level_loading_log(n_cal, n_a, L, within_log_limit) else NULL

  Lambda_imp_alpha <- matrix(0.0, n_cells, L)
  Lambda_imp_beta <- lapply(imp_groups, function(g) matrix(0.0, g$p, L))
  Lambda_or_alpha <- if (n_or_blks > 0L) matrix(0.0, n_or_blks, L) else NULL
  n_or_grps <- if (!is.null(or_groups)) length(or_groups$groups) else 0L
  Lambda_or_beta <- if (n_or_grps > 0L) lapply(or_groups$groups, function(g) matrix(0.0, length(g$keep), L)) else NULL
  Lambda_a <- if (compute_if) matrix(0.0, G_coh, L) else NULL
  Lambda_theta <- if (compute_if && D_theta > 0L) matrix(0.0, D_theta, L) else NULL
  Lambda_trt <- list()

  for (ell in seq_len(L)) {
    d <- descs[[ell]]
    ci_g <- link_ci_g[ell]

    # ---- Target outcome arm (+T_gl) ------------------------------------------
    T_gl <- 0.0
    if (d$has_target_outcome) {
      cid_g <- link_tgt_cid[ell]
      cd_g <- cell_data[[cid_g]]
      dp <- M_one_coh[ci_g]
      T_gl <- sum(cd_g$w_syn * cd_g$Z) / (N_orig * dp)

      if (compute_if) {
        Phi_lay[cd_g$rf, ell] <- Phi_lay[cd_g$rf, ell] + (cd_g$w_syn / dp) * (cd_g$Z - T_gl)
      }
      # Sensitivity to Stage-1 cell `cid_g`
      w_j <- cd_g$w_syn * cd_g$one_minus_spi / (N_orig * dp)
      if (recipe$intercept) Lambda_imp_alpha[cid_g, ell] <- Lambda_imp_alpha[cid_g, ell] + sum(w_j)
      if (cd_g$gf$p > 0L) {
        Lambda_imp_beta[[cd_g$gf$group_id]][, ell] <- Lambda_imp_beta[[cd_g$gf$group_id]][, ell] +
          as.vector(crossprod(w_j, cd_g$d_imp))
      }
      if (augmented) {
        if (compute_if && length(cd_g$theta_ch) && length(cd_g$obs_local)) {
          w_th <- (cd_g$w_syn[cd_g$obs_local] / (N_orig * dp)) * cd_g$pi_obs * cd_g$e_imp_obs
          Lambda_theta[cd_g$theta_ch, ell] <- Lambda_theta[cd_g$theta_ch, ell] +
            as.vector(crossprod(w_th, cd_g$T_logit))
        }
        if (compute_within && length(cd_g$obs_local)) {
          w_dy <- (cd_g$w_syn[cd_g$obs_local] / (N_orig * dp)) * cd_g$pi_obs
          loading_log_pair_vec(U_log, cd_g$a_obs, d$ct, d$cs, ell, w_dy)
        }
      }
    }

    if (!d$has_control) {
      mhat[ell] <- T_gl
      next
    }

    # ---- Control arm (sgn * (...)) -------------------------------------------
    sgn <- d$ctrl_sign
    dp_g <- M_one_coh[ci_g]
    rf_g <- rows_by_ci[[ci_g]]
    w_syn_g <- st$w[rf_g] * st$a[ci_g]

    B_gl <- 0.0
    d_bar_plus <- numeric(0L)
    bf <- NULL
    if (need_or) {
      bid <- link_or_bid[ell]
      bf <- blk_fits[[bid]]
      kept_F <- slope_cols_F[bf$keep]
      d_bar_plus <- if (bf$p > 0L) (d$L_steps * M_z_coh[[ci_g]][bf$keep]) / dp_g else numeric(0L)
      B_gl <- (if (recipe$intercept) bf$alpha else 0.0) +
        (if (bf$p > 0L) sum(d_bar_plus * bf$beta) else 0.0)
      if (compute_if) {
        d_g_or <- if (bf$p > 0L) d$L_steps * X_F[rf_g, kept_F, drop = FALSE] else matrix(0.0, length(rf_g), 0L)
        mhat_g <- (if (recipe$intercept) bf$alpha else 0.0) +
          (if (bf$p > 0L) as.vector(d_g_or %*% bf$beta) else numeric(length(rf_g)))
        Phi_lay[rf_g, ell] <- Phi_lay[rf_g, ell] + sgn * (w_syn_g / dp_g) * (mhat_g - B_gl)
      }
    }

    C_ol <- 0.0
    B_ol <- 0.0
    d_bar_minus <- numeric(0L)
    if (need_trt) {
      ci_K_vec <- link_ci_K[[ell]]
      cids_k <- link_ctrl_cids[[ell]]
      tfit <- get_syn_trt(ci_g, ci_K_vec, d)
      tid <- tfit$tid
      if (length(Lambda_trt) < tid || is.null(Lambda_trt[[tid]])) {
        Lambda_trt[[tid]] <- matrix(0.0, trt_recipe$p, L)
      }

      # Gather control rows across k in K in their cell order
      sum_wm <- 0.0
      sum_wm_Z <- 0.0
      sum_wm_m <- 0.0
      sum_wm_dor <- if (is_dr) numeric(bf$p) else numeric(0L)
      per_k <- vector("list", length(cids_k))

      for (kk in seq_along(cids_k)) {
        cd_k <- cell_data[[cids_k[kk]]]
        pos_in_trt <- tfit$pos_by_ci[[match(cd_k$ci, tfit$all_ci)]]
        odds_k <- tfit$odds[pos_in_trt]
        act_k <- if (compute_if) tfit$active[pos_in_trt] else NULL
        X_trt_k <- if (compute_if) tfit$X_trt[pos_in_trt, , drop = FALSE] else NULL
        wm_k <- cd_k$w_syn * odds_k
        sum_wm <- sum_wm + sum(wm_k)
        sum_wm_Z <- sum_wm_Z + sum(wm_k * cd_k$Z)

        mhat_k <- NULL
        d_or_k <- NULL
        if (is_dr) {
          kept_F <- slope_cols_F[bf$keep]
          d_or_k <- if (bf$p > 0L) d$L_steps * X_F[cd_k$rf, kept_F, drop = FALSE] else matrix(0.0, length(cd_k$rf), 0L)
          mhat_k <- (if (recipe$intercept) bf$alpha else 0.0) +
            (if (bf$p > 0L) as.vector(d_or_k %*% bf$beta) else numeric(length(cd_k$rf)))
          sum_wm_m <- sum_wm_m + sum(wm_k * mhat_k)
          if (bf$p > 0L) sum_wm_dor <- sum_wm_dor + colSums(wm_k * d_or_k)
        }
        per_k[[kk]] <- list(cd = cd_k, wm = wm_k, act = act_k, X_trt = X_trt_k, mhat = mhat_k)
      }

      dm <- sum_wm / N_orig
      cp04_finite(dm, "synthetic control normalizer", TRUE)
      C_ol <- sum_wm_Z / sum_wm
      if (is_dr) {
        B_ol <- sum_wm_m / sum_wm
        d_bar_minus <- if (bf$p > 0L) sum_wm_dor / sum_wm else numeric(0L)
      }

      for (kk in seq_along(cids_k)) {
        pk <- per_k[[kk]]
        cd_k <- pk$cd
        cid_k <- cd_k$cid
        wm_k <- pk$wm
        eff_res <- if (is_dr) (cd_k$Z - pk$mhat) - (C_ol - B_ol) else (cd_k$Z - C_ol)

        if (compute_if) {
          Phi_lay[cd_k$rf, ell] <- Phi_lay[cd_k$rf, ell] + sgn * (wm_k / dm) * eff_res
          # Multi-cohort a_k partial on C_ol - B_ol (auto-zero when |K| == 1)
          Lambda_a[cd_k$ci, ell] <- Lambda_a[cd_k$ci, ell] +
            sgn * sum(wm_k * eff_res) / (N_orig * dm * st$a[cd_k$ci])
          # Propensity gamma sensitivity
          Lambda_trt[[tid]][, ell] <- Lambda_trt[[tid]][, ell] +
            sgn * colSums((wm_k * pk$act * eff_res) * pk$X_trt) / (N_orig * dm)
        }

        # Stage-1 sensitivity through C_ol
        w_j_k <- sgn * wm_k * cd_k$one_minus_spi / (N_orig * dm)
        if (recipe$intercept) Lambda_imp_alpha[cid_k, ell] <- Lambda_imp_alpha[cid_k, ell] + sum(w_j_k)
        if (cd_k$gf$p > 0L) {
          Lambda_imp_beta[[cd_k$gf$group_id]][, ell] <- Lambda_imp_beta[[cd_k$gf$group_id]][, ell] +
            as.vector(crossprod(w_j_k, cd_k$d_imp))
        }
        if (augmented) {
          if (compute_if && length(cd_k$theta_ch) && length(cd_k$obs_local)) {
            w_th_k <- sgn * (wm_k[cd_k$obs_local] / (N_orig * dm)) * cd_k$pi_obs * cd_k$e_imp_obs
            Lambda_theta[cd_k$theta_ch, ell] <- Lambda_theta[cd_k$theta_ch, ell] +
              as.vector(crossprod(w_th_k, cd_k$T_logit))
          }
          if (compute_within && length(cd_k$obs_local)) {
            w_dy_k <- sgn * (wm_k[cd_k$obs_local] / (N_orig * dm)) * cd_k$pi_obs
            loading_log_pair_vec(U_log, cd_k$a_obs, d$ct, d$cs, ell, w_dy_k)
          }
        }
      }
    }

    if (is_or) {
      mhat[ell] <- T_gl + sgn * B_gl
      if (recipe$intercept) Lambda_or_alpha[link_or_bid[ell], ell] <- Lambda_or_alpha[link_or_bid[ell], ell] + sgn
      if (bf$p > 0L) Lambda_or_beta[[bf$group_id]][, ell] <- Lambda_or_beta[[bf$group_id]][, ell] + sgn * d_bar_plus
    } else if (is_dr) {
      mhat[ell] <- T_gl + sgn * (B_gl + C_ol - B_ol)
      if (bf$p > 0L) {
        Lambda_or_beta[[bf$group_id]][, ell] <- Lambda_or_beta[[bf$group_id]][, ell] + sgn * (d_bar_plus - d_bar_minus)
      }
    } else {
      mhat[ell] <- T_gl + sgn * C_ol
    }
  }

  # Backpropagate Stage-3 propensity logits -> Phi_lay + Lambda_a
  if (compute_if && length(trt_list) > 0L) {
    for (tid in seq_along(trt_list)) {
      if (tid > length(Lambda_trt) || is.null(Lambda_trt[[tid]])) next
      L_gam <- Lambda_trt[[tid]]
      cols_nz <- which(colSums(L_gam != 0) > 0L)
      if (!length(cols_nz)) next
      tfit <- trt_list[[tid]]
      nu_gam <- cp04_solve(tfit$H_chol, L_gam[, cols_nz, drop = FALSE])
      sc_rows <- (tfit$w_syn * (tfit$y_r - tfit$prob)) * (tfit$X_trt %*% nu_gam)
      Phi_lay[tfit$risk_rows, cols_nz] <- Phi_lay[tfit$risk_rows, cols_nz, drop = FALSE] + sc_rows
      Lambda_a[tfit$all_ci, cols_nz] <- Lambda_a[tfit$all_ci, cols_nz, drop = FALSE] + crossprod(tfit$C_a, nu_gam)
    }
  }

  # Backpropagate Stage-2 OR -> Phi_lay + Stage-1 (Lambda_imp_*) + Stage-0 (Lambda_a, Lambda_theta) + U_log
  if (n_or_blks > 0L) {
    for (bid in seq_len(n_or_blks)) {
      bf <- blk_fits[[bid]]
      if (recipe$intercept && bf$p > 0L) {
        lam_a <- Lambda_or_alpha[bid, ]
        if (any(lam_a != 0)) {
          Lambda_or_beta[[bf$group_id]] <- Lambda_or_beta[[bf$group_id]] - outer(as.double(bf$mean_d), lam_a)
        }
      }
    }
    nu_or_beta <- vector("list", n_or_grps)
    cols_or_beta <- vector("list", n_or_grps)
    for (gid in seq_len(n_or_grps)) {
      g_obj <- or_groups$groups[[gid]]
      Lam_b <- Lambda_or_beta[[gid]]
      p_g <- nrow(Lam_b)
      cb <- if (p_g > 0L) which(colSums(Lam_b != 0) > 0L) else integer(0L)
      if (length(cb) > 0L) {
        nu_raw <- cross_pair_gram_solve(g_obj$fit, recipe$intercept, Lam_b[, cb, drop = FALSE], N_orig)
        keep_c <- which(colSums(nu_raw != 0) > 0L)
        cols_or_beta[[gid]] <- cb[keep_c]
        nu_or_beta[[gid]] <- nu_raw[, keep_c, drop = FALSE]
      } else {
        cols_or_beta[[gid]] <- integer(0L)
        nu_or_beta[[gid]] <- matrix(0.0, nrow = p_g, ncol = 0L)
      }
    }

    blks_by_gid <- split(seq_len(n_or_blks), vapply(blk_fits, `[[`, integer(1L), "group_id"))
    for (gid_chr in names(blks_by_gid)) {
      gid <- as.integer(gid_chr)
      bids_g <- blks_by_gid[[gid_chr]]
      nu_sub <- nu_or_beta[[gid]]
      cols_b <- cols_or_beta[[gid]]
      p_g <- nrow(nu_sub)
      has_slope <- length(cols_b) > 0L
      all_sids_g <- unique(unlist(lapply(bids_g, function(b) blk_fits[[b]]$support_ids), use.names = FALSE))
      all_cids_g <- vapply(all_sids_g, function(s) sup_or[[s]]$cid, integer(1L))

      if (length(bids_g) == 1L && length(all_sids_g) == 1L) {
        bid <- bids_g[1L]
        bf <- blk_fits[[bid]]
        lam_a <- if (recipe$intercept) Lambda_or_alpha[bid, ] else numeric(L)
        cols_a <- which(lam_a != 0)
        has_int <- length(cols_a) > 0L
        if (!has_int && !has_slope) next
        cid <- sup_or[[bf$support_ids[1L]]]$cid
        cd <- cell_data[[cid]]
        kept_F <- slope_cols_F[bf$keep]
        d_or <- if (p_g > 0L) bf$L * X_F[cd$rf, kept_F, drop = FALSE] else NULL
        mhat_or <- (if (recipe$intercept) bf$alpha else 0.0) +
          (if (p_g > 0L) as.vector(d_or %*% bf$beta) else numeric(length(cd$rf)))
        r_or <- cd$Z - mhat_or

        if (has_int) {
          Q_int <- outer(cd$w_syn / bf$weight_sum, lam_a[cols_a])
          if (compute_if) {
            Phi_lay[cd$rf, cols_a] <- Phi_lay[cd$rf, cols_a, drop = FALSE] + (N_orig * r_or) * Q_int
            Lambda_a[cd$ci, cols_a] <- Lambda_a[cd$ci, cols_a] + colSums((r_or / st$a[cd$ci]) * Q_int)
          }
          Q_j_int <- cd$one_minus_spi * Q_int
          if (recipe$intercept) Lambda_imp_alpha[cid, cols_a] <- Lambda_imp_alpha[cid, cols_a] + colSums(Q_j_int)
          if (cd$gf$p > 0L) {
            Lambda_imp_beta[[cd$gf$group_id]][, cols_a] <-
              Lambda_imp_beta[[cd$gf$group_id]][, cols_a, drop = FALSE] + crossprod(cd$d_imp, Q_j_int)
          }
          if (augmented && length(cd$obs_local)) {
            Q_obs_int <- Q_int[cd$obs_local, , drop = FALSE]
            if (compute_if && length(cd$theta_ch)) {
              Lambda_theta[cd$theta_ch, cols_a] <- Lambda_theta[cd$theta_ch, cols_a, drop = FALSE] +
                crossprod(cd$T_logit, (cd$pi_obs * cd$e_imp_obs) * Q_obs_int)
            }
            if (compute_within) {
              loading_log_pair_outer(U_log, cd$a_obs, cd$cf$pair$t_calendar_index,
                cd$cf$pair$s_calendar_index, cols_a, cd$w_syn / bf$weight_sum, lam_a[cols_a],
                src = cd$obs_local, pre = cd$pi_obs)
            }
          }
        }

        if (has_slope) {
          d_tilde <- if (recipe$intercept) sweep(d_or, 2L, bf$mean_d, "-") else d_or
          H_nu <- ((cd$w_syn / N_orig) * d_tilde) %*% nu_sub
          if (compute_if) {
            Phi_lay[cd$rf, cols_b] <- Phi_lay[cd$rf, cols_b, drop = FALSE] + (N_orig * r_or) * H_nu
            Lambda_a[cd$ci, cols_b] <- Lambda_a[cd$ci, cols_b] + colSums((r_or / st$a[cd$ci]) * H_nu)
          }
          H_j_nu <- cd$one_minus_spi * H_nu
          if (recipe$intercept) Lambda_imp_alpha[cid, cols_b] <- Lambda_imp_alpha[cid, cols_b] + colSums(H_j_nu)
          if (cd$gf$p > 0L) {
            Lambda_imp_beta[[cd$gf$group_id]][, cols_b] <-
              Lambda_imp_beta[[cd$gf$group_id]][, cols_b, drop = FALSE] + crossprod(cd$d_imp, H_j_nu)
          }
          if (augmented && length(cd$obs_local)) {
            H_obs_nu <- H_nu[cd$obs_local, , drop = FALSE]
            if (compute_if && D_theta > 0L && length(cd$theta_ch)) {
              Lambda_theta[cd$theta_ch, cols_b] <- Lambda_theta[cd$theta_ch, cols_b, drop = FALSE] +
                crossprod(cd$T_logit, (cd$pi_obs * cd$e_imp_obs) * H_obs_nu)
            }
            if (compute_within) {
              loading_log_lowrank(U_log, cd$a_obs, cd$cf$pair$t_calendar_index,
                cd$cf$pair$s_calendar_index, cols_b, (cd$w_syn / N_orig) * d_tilde, nu_sub,
                src = cd$obs_local, pre = cd$pi_obs)
            }
          }
        }
        next
      }

      u_rf <- if (has_slope && compute_if) {
        sort(unique(unlist(lapply(all_cids_g, function(cid) cell_data[[cid]]$rf), use.names = FALSE)))
      } else integer(0L)
      u_ci <- if (has_slope && compute_if) {
        sort(unique(vapply(all_cids_g, function(cid) cell_data[[cid]]$ci, integer(1L))))
      } else integer(0L)
      u_cids <- if (has_slope && recipe$intercept) sort(unique(all_cids_g)) else integer(0L)
      touched_ig <- if (has_slope) {
        sort(unique(vapply(all_cids_g, function(cid) cell_data[[cid]]$gf$group_id, integer(1L))))
      } else integer(0L)
      u_th <- if (has_slope && augmented && compute_if && D_theta > 0L) {
        sort(unique(unlist(lapply(all_cids_g, function(cid) cell_data[[cid]]$theta_ch), use.names = FALSE)))
      } else integer(0L)
      u_a_obs <- if (has_slope && augmented && compute_within) {
        sort(unique(unlist(lapply(all_cids_g, function(cid) cell_data[[cid]]$a_obs), use.names = FALSE)))
      } else integer(0L)
      active_taus <- if (has_slope && augmented && compute_within) {
        sort(unique(unlist(lapply(all_cids_g, function(cid) {
          c(cell_data[[cid]]$cf$pair$s_calendar_index, cell_data[[cid]]$cf$pair$t_calendar_index)
        }), use.names = FALSE)))
      } else integer(0L)

      S_phi_g <- if (length(u_rf) > 0L) matrix(0.0, length(u_rf), p_g) else NULL
      L_a_g <- if (length(u_ci) > 0L) matrix(0.0, length(u_ci), p_g) else NULL
      L_imp_a_g <- if (length(u_cids) > 0L) matrix(0.0, length(u_cids), p_g) else NULL
      L_imp_b_g <- if (length(touched_ig) > 0L) {
        lapply(touched_ig, function(ig_id) matrix(0.0, imp_groups[[ig_id]]$p, p_g))
      } else NULL
      L_th_g <- if (length(u_th) > 0L) matrix(0.0, length(u_th), p_g) else NULL
      U_grp <- if (length(active_taus) > 0L && length(u_a_obs) > 0L) {
        replicate(length(active_taus), matrix(0.0, length(u_a_obs), p_g), simplify = FALSE)
      } else NULL

      for (bid in bids_g) {
        bf <- blk_fits[[bid]]
        lam_a <- if (recipe$intercept) Lambda_or_alpha[bid, ] else numeric(L)
        cols_a <- which(lam_a != 0)
        has_int <- length(cols_a) > 0L
        if (!has_int && !has_slope) next
        lam_a_sub <- if (has_int) lam_a[cols_a] else NULL
        kept_F <- slope_cols_F[bf$keep]

        for (sid in bf$support_ids) {
          cid <- sup_or[[sid]]$cid
          cd <- cell_data[[cid]]
          d_or <- if (p_g > 0L) bf$L * X_F[cd$rf, kept_F, drop = FALSE] else NULL
          mhat_or <- (if (recipe$intercept) bf$alpha else 0.0) +
            (if (p_g > 0L) as.vector(d_or %*% bf$beta) else numeric(length(cd$rf)))
          r_or <- cd$Z - mhat_or

          if (has_int) {
            Q_int <- outer(cd$w_syn / bf$weight_sum, lam_a_sub)
            if (compute_if) {
              Phi_lay[cd$rf, cols_a] <- Phi_lay[cd$rf, cols_a, drop = FALSE] + (N_orig * r_or) * Q_int
              Lambda_a[cd$ci, cols_a] <- Lambda_a[cd$ci, cols_a] + colSums((r_or / st$a[cd$ci]) * Q_int)
            }
            Q_j_int <- cd$one_minus_spi * Q_int
            if (recipe$intercept) Lambda_imp_alpha[cid, cols_a] <- Lambda_imp_alpha[cid, cols_a] + colSums(Q_j_int)
            if (cd$gf$p > 0L) {
              Lambda_imp_beta[[cd$gf$group_id]][, cols_a] <-
                Lambda_imp_beta[[cd$gf$group_id]][, cols_a, drop = FALSE] + crossprod(cd$d_imp, Q_j_int)
            }
            if (augmented && length(cd$obs_local)) {
              Q_obs_int <- Q_int[cd$obs_local, , drop = FALSE]
              if (compute_if && length(cd$theta_ch)) {
                Lambda_theta[cd$theta_ch, cols_a] <- Lambda_theta[cd$theta_ch, cols_a, drop = FALSE] +
                  crossprod(cd$T_logit, (cd$pi_obs * cd$e_imp_obs) * Q_obs_int)
              }
              if (compute_within) {
                loading_log_pair_outer(U_log, cd$a_obs, cd$cf$pair$t_calendar_index,
                  cd$cf$pair$s_calendar_index, cols_a, cd$w_syn / bf$weight_sum, lam_a[cols_a],
                  src = cd$obs_local, pre = cd$pi_obs)
              }
            }
          }

          if (has_slope) {
            d_tilde <- if (recipe$intercept) sweep(d_or, 2L, bf$mean_d, "-") else d_or
            H_or <- (cd$w_syn / N_orig) * d_tilde
            if (compute_if) {
              if (!is.null(S_phi_g) && length(cd$rf)) {
                rf_loc <- match(cd$rf, u_rf)
                S_phi_g[rf_loc, ] <- S_phi_g[rf_loc, , drop = FALSE] + (N_orig * r_or) * H_or
              }
              if (!is.null(L_a_g)) {
                ci_loc <- match(cd$ci, u_ci)
                L_a_g[ci_loc, ] <- L_a_g[ci_loc, ] + colSums((r_or / st$a[cd$ci]) * H_or)
              }
            }
            H_j <- cd$one_minus_spi * H_or
            if (recipe$intercept && !is.null(L_imp_a_g)) {
              cid_loc <- match(cid, u_cids)
              L_imp_a_g[cid_loc, ] <- L_imp_a_g[cid_loc, ] + colSums(H_j)
            }
            if (cd$gf$p > 0L && !is.null(L_imp_b_g)) {
              ig_loc <- match(cd$gf$group_id, touched_ig)
              L_imp_b_g[[ig_loc]] <- L_imp_b_g[[ig_loc]] + crossprod(cd$d_imp, H_j)
            }
            if (augmented && length(cd$obs_local)) {
              H_obs <- H_or[cd$obs_local, , drop = FALSE]
              if (!is.null(L_th_g) && length(cd$theta_ch)) {
                th_loc <- match(cd$theta_ch, u_th)
                L_th_g[th_loc, ] <- L_th_g[th_loc, , drop = FALSE] +
                  crossprod(cd$T_logit, (cd$pi_obs * cd$e_imp_obs) * H_obs)
              }
              if (compute_within && !is.null(U_grp) && length(cd$a_obs)) {
                a_loc <- match(cd$a_obs, u_a_obs)
                ts_loc <- match(cd$cf$pair$s_calendar_index, active_taus)
                tt_loc <- match(cd$cf$pair$t_calendar_index, active_taus)
                U_grp <- add_pair_diff_to_U(U_grp, a_loc, ts_loc, tt_loc, cd$pi_obs * H_obs)
              }
            }
          }
        }
      }

      if (has_slope) {
        if (compute_if) {
          if (!is.null(S_phi_g)) {
            Phi_lay[u_rf, cols_b] <- Phi_lay[u_rf, cols_b, drop = FALSE] + S_phi_g %*% nu_sub
          }
          if (!is.null(L_a_g)) {
            Lambda_a[u_ci, cols_b] <- Lambda_a[u_ci, cols_b, drop = FALSE] + L_a_g %*% nu_sub
          }
          if (!is.null(L_th_g) && any(L_th_g != 0)) {
            Lambda_theta[u_th, cols_b] <- Lambda_theta[u_th, cols_b, drop = FALSE] + L_th_g %*% nu_sub
          }
        }
        if (recipe$intercept && !is.null(L_imp_a_g) && any(L_imp_a_g != 0)) {
          Lambda_imp_alpha[u_cids, cols_b] <- Lambda_imp_alpha[u_cids, cols_b, drop = FALSE] + L_imp_a_g %*% nu_sub
        }
        if (!is.null(L_imp_b_g)) {
          for (ig_i in seq_along(touched_ig)) {
            ig_id <- touched_ig[ig_i]
            M_ig <- L_imp_b_g[[ig_i]]
            if (nrow(M_ig) > 0L && any(M_ig != 0)) {
              Lambda_imp_beta[[ig_id]][, cols_b] <- Lambda_imp_beta[[ig_id]][, cols_b, drop = FALSE] +
                M_ig %*% nu_sub
            }
          }
        }
        if (augmented && compute_within && !is.null(U_grp)) {
          for (ti in seq_along(active_taus)) {
            if (any(U_grp[[ti]] != 0)) {
              tau <- active_taus[ti]
              loading_log_lowrank(U_log, u_a_obs, tau, 0L, cols_b, U_grp[[ti]], nu_sub)
            }
          }
        }
      }
    }
  }

  # Backpropagate Stage-1 (`Lambda_imp_alpha`, `Lambda_imp_beta`) -> observed pairs
  for (cid in seq_len(n_cells)) {
    cd <- cell_data[[cid]]
    if (recipe$intercept && cd$gf$p > 0L) {
      lam_a <- Lambda_imp_alpha[cid, ]
      if (any(lam_a != 0)) {
        Lambda_imp_beta[[cd$gf$group_id]] <- Lambda_imp_beta[[cd$gf$group_id]] -
          outer(as.double(cd$cf$mean_d), lam_a)
      }
    }
  }
  nu_imp_beta <- vector("list", n_imp_grps)
  cols_imp_beta <- vector("list", n_imp_grps)
  for (gid in seq_len(n_imp_grps)) {
    gf <- imp_groups[[gid]]
    Lam_b <- Lambda_imp_beta[[gid]]
    p_g <- nrow(Lam_b)
    cb <- if (p_g > 0L) which(colSums(Lam_b != 0) > 0L) else integer(0L)
    if (length(cb) > 0L) {
      nu_raw <- cross_pair_gram_solve(gf$fit, recipe$intercept, Lam_b[, cb, drop = FALSE], N_orig)
      keep_c <- which(colSums(nu_raw != 0) > 0L)
      cols_imp_beta[[gid]] <- cb[keep_c]
      nu_imp_beta[[gid]] <- nu_raw[, keep_c, drop = FALSE]
    } else {
      cols_imp_beta[[gid]] <- integer(0L)
      nu_imp_beta[[gid]] <- matrix(0.0, nrow = p_g, ncol = 0L)
    }
  }

  cells_by_gid <- split(seq_len(n_cells), vapply(cell_data, function(cd) cd$gf$group_id, integer(1L)))
  for (gid_chr in names(cells_by_gid)) {
    gid <- as.integer(gid_chr)
    cids_g <- cells_by_gid[[gid_chr]]
    nu_sub <- nu_imp_beta[[gid]]
    cols_b <- cols_imp_beta[[gid]]
    p_g <- nrow(nu_sub)
    has_slope <- length(cols_b) > 0L

    if (length(cids_g) == 1L) {
      cid <- cids_g[1L]
      cd <- cell_data[[cid]]
      lam_a <- if (recipe$intercept) Lambda_imp_alpha[cid, ] else numeric(L)
      cols_a <- which(lam_a != 0)
      has_int <- length(cols_a) > 0L
      if (!has_int && !has_slope) next
      w_obs <- cd$w_raw[cd$obs_local]

      if (has_int) {
        W_int <- outer(w_obs / cd$cf$weight_sum, lam_a[cols_a])
        if (compute_if) {
          Phi_lay[cd$rows_obs, cols_a] <- Phi_lay[cd$rows_obs, cols_a, drop = FALSE] + (N_orig * cd$e_imp_obs) * W_int
        }
        if (compute_within) {
          loading_log_pair_outer(U_log, cd$a_obs, cd$cf$pair$t_calendar_index,
            cd$cf$pair$s_calendar_index, cols_a, w_obs / cd$cf$weight_sum, lam_a[cols_a])
        }
      }

      if (has_slope) {
        H_nu <- ((w_obs / N_orig) * cd$d_imp_tilde) %*% nu_sub
        if (compute_if) {
          Phi_lay[cd$rows_obs, cols_b] <- Phi_lay[cd$rows_obs, cols_b, drop = FALSE] + (N_orig * cd$e_imp_obs) * H_nu
        }
        if (compute_within) {
          loading_log_lowrank(U_log, cd$a_obs, cd$cf$pair$t_calendar_index,
            cd$cf$pair$s_calendar_index, cols_b, (w_obs / N_orig) * cd$d_imp_tilde, nu_sub)
        }
      }
      next
    }

    u_rows_obs <- if (has_slope && compute_if) {
      sort(unique(unlist(lapply(cids_g, function(cid) cell_data[[cid]]$rows_obs), use.names = FALSE)))
    } else integer(0L)
    u_a_obs <- if (has_slope && compute_within) {
      sort(unique(unlist(lapply(cids_g, function(cid) cell_data[[cid]]$a_obs), use.names = FALSE)))
    } else integer(0L)
    active_taus <- if (has_slope && compute_within) {
      sort(unique(unlist(lapply(cids_g, function(cid) {
        c(cell_data[[cid]]$cf$pair$s_calendar_index, cell_data[[cid]]$cf$pair$t_calendar_index)
      }), use.names = FALSE)))
    } else integer(0L)

    S_phi_imp <- if (length(u_rows_obs) > 0L) matrix(0.0, length(u_rows_obs), p_g) else NULL
    U_grp_imp <- if (length(active_taus) > 0L && length(u_a_obs) > 0L) {
      replicate(length(active_taus), matrix(0.0, length(u_a_obs), p_g), simplify = FALSE)
    } else NULL

    for (cid in cids_g) {
      cd <- cell_data[[cid]]
      lam_a <- if (recipe$intercept) Lambda_imp_alpha[cid, ] else numeric(L)
      cols_a <- which(lam_a != 0)
      has_int <- length(cols_a) > 0L
      if (!has_int && !has_slope) next
      w_obs <- cd$w_raw[cd$obs_local]

      if (has_int) {
        W_int <- outer(w_obs / cd$cf$weight_sum, lam_a[cols_a])
        if (compute_if) {
          Phi_lay[cd$rows_obs, cols_a] <- Phi_lay[cd$rows_obs, cols_a, drop = FALSE] + (N_orig * cd$e_imp_obs) * W_int
        }
        if (compute_within) {
          loading_log_pair_outer(U_log, cd$a_obs, cd$cf$pair$t_calendar_index,
            cd$cf$pair$s_calendar_index, cols_a, w_obs / cd$cf$weight_sum, lam_a[cols_a])
        }
      }

      if (has_slope) {
        H_imp <- (w_obs / N_orig) * cd$d_imp_tilde
        if (compute_if && !is.null(S_phi_imp) && length(cd$rows_obs)) {
          r_loc <- match(cd$rows_obs, u_rows_obs)
          S_phi_imp[r_loc, ] <- S_phi_imp[r_loc, , drop = FALSE] + (N_orig * cd$e_imp_obs) * H_imp
        }
        if (compute_within && !is.null(U_grp_imp) && length(cd$a_obs)) {
          a_loc <- match(cd$a_obs, u_a_obs)
          ts_loc <- match(cd$cf$pair$s_calendar_index, active_taus)
          tt_loc <- match(cd$cf$pair$t_calendar_index, active_taus)
          U_grp_imp <- add_pair_diff_to_U(U_grp_imp, a_loc, ts_loc, tt_loc, H_imp)
        }
      }
    }

    if (has_slope) {
      if (compute_if && !is.null(S_phi_imp)) {
        Phi_lay[u_rows_obs, cols_b] <- Phi_lay[u_rows_obs, cols_b, drop = FALSE] + S_phi_imp %*% nu_sub
      }
      if (compute_within && !is.null(U_grp_imp)) {
        for (ti in seq_along(active_taus)) {
          if (any(U_grp_imp[[ti]] != 0)) {
            tau <- active_taus[ti]
            loading_log_lowrank(U_log, u_a_obs, tau, 0L, cols_b, U_grp_imp[[ti]], nu_sub)
          }
        }
      }
    }
  }

  # Backpropagate Stage-0 (`Lambda_a` and `Lambda_theta`) -> Phi_lay, then fold layout to N_orig
  Phi_F <- NULL
  if (compute_if) {
    for (ci in active_cohorts) {
      la <- Lambda_a[ci, ]
      cols_nz <- which(la != 0)
      if (!length(cols_nz)) next
      rf <- rows_by_ci[[ci]]
      A_rf <- as.double(st$a_index[rf] > 0L)
      xi_a <- (st$w[rf] / st$mass_F[ci]) * (A_rf - st$a[ci])
      Phi_lay[rf, cols_nz] <- Phi_lay[rf, cols_nz, drop = FALSE] + outer(xi_a, la[cols_nz])
    }
    if (D_theta > 0L && any(Lambda_theta != 0)) {
      cols_nz <- which(colSums(Lambda_theta != 0) > 0L)
      if (length(cols_nz)) {
        Phi_lay[, cols_nz] <- Phi_lay[, cols_nz, drop = FALSE] + Xi_theta %*% Lambda_theta[, cols_nz, drop = FALSE]
      }
    }
    Phi_F <- if (is.null(lay)) Phi_lay else chained_reference_layout_sum(Phi_lay, lay$orig_row, N_orig)
  }
  if (!is.null(nuisance_sink)) cross_pair_nuisance_reconstruction_collect(nuisance_sink, environment())

  ref_env <- if (!is.null(reference)) {
    e <- new.env(parent = emptyenv()); e$reference_reader <- reference; e
  } else NULL

  acc <- list(
    automatic_chunk_rows = isTRUE(cp04$automatic_chunk_rows),
    n_cells = as.integer(n_cells),
    n_groups = as.integer(n_imp_grps),
    n_imputation_fits = as.integer(n_imp_grps),
    n_imputation_local_intercepts = if (recipe$intercept) as.integer(n_cells) else 0L,
    n_or_fits = if (!is.null(or_groups)) length(or_groups$groups) else 0L,
    n_or_blocks = as.integer(n_or_blks),
    n_treatment_models = length(trt_list),
    n_treatment_fits = length(trt_list),
    n_treatment_columns = as.integer(trt_recipe$p),
    n_selection_fits = if (!is.null(sel_fit)) length(sel_fit$blocks) else 0L,
    n_selection_cells = if (!is.null(sel_fit)) length(sel_fit$cells) else 0L,
    n_selection_blocks = if (!is.null(sel_fit)) length(sel_fit$blocks) else 0L,
    n_selection_parameters = if (!is.null(sel_fit)) length(sel_fit$coefficients) else 0L
  )

  backend <- structure(
    list(
      mode = if (augmented) "augmented" else "imputation",
      est_method = est_method,
      outcome_pooling = outcome_pooling,
      imputation_pooling = imputation_pooling,
      outcome_recipe = recipe,
      reader_ref = reader$ref,
      reference_ref = ref_env,
      uses_selection_models = !is.null(sel_fit),
      selection_backend = sel_fit,
      accounting = acc
    ),
    class = "cdid_cross_pair_backend"
  )

  # Centered in place, one column at a time, when requested: entrywise
  # `Phi_F - colMeans(Phi_F)`, with transients of one column.
  Phi_mean <- NULL
  if (isTRUE(center_if) && !is.null(Phi_F)) {
    Phi_mean <- colMeans(Phi_F)
    for (j in seq_len(ncol(Phi_F))) Phi_F[, j] <- Phi_F[, j] - Phi_mean[j]
  }
  list(mhat = mhat, Phi_F = Phi_F, Phi_mean = Phi_mean, U_log = U_log, backend = backend)
}

# Public entry for cross-pair general evaluation (all selection modes, pooled OR,
# reference inputs, and within covariance).
evaluate_cross_pair_general_moments <- function(links_all, panel_owner, frame_owner,
                                                cp04, cp04_reader, f_Z, pm,
                                                est_method, ctrl, wt_mode,
                                                reduce_structural_redundancy = TRUE,
                                                compute_if = TRUE, nuisance_sink = NULL,
                                                defer_covariance = FALSE) {
  sel_mode <- if (is.null(cp04)) "none" else if (!is.null(cp04$selection_est_method)) {
    cp04$selection_est_method
  } else if (identical(cp04$canonical_assumption, "conditional_mtar")) {
    "ipw"
  } else {
    "none"
  }
  outcome_pooling <- if (is.null(cp04) || is.null(cp04$outcome_pooling)) "separate" else cp04$outcome_pooling
  imputation_pooling <- if (!is.null(cp04) && !is.null(cp04$imputation_pooling)) cp04$imputation_pooling else "separate"

  built <- build_general_link_descriptors(
    links_df = links_all, all_times = pm$all_times, reader = cp04_reader,
    O_mat = pm$O_mat, sel_mode = sel_mode, is_cross_pair = TRUE
  )
  descs <- built$descs
  views <- built$views
  M_all <- length(descs)

  ctrl_eff <- ctrl
  ctrl_eff$owner_control <- cp04_reader$control

  reference <- if (!is.null(cp04) && !is.null(cp04$reference_source)) {
    src_ref <- cp04$reference_source
    recover <- identical(cp04$missing_reference, "recover")
    can_dedup_views <- isTRUE(cp04$automatic_chunk_rows) &&
      (!isTRUE(src_ref$formulas$x$dated) || !(sel_mode %in% c("none", "ipw"))) &&
      M_all > 64L
    if (can_dedup_views) {
      st_r <- cp04_reader$state
      xf_r <- src_ref$formulas$x
      sf_r <- src_ref$formulas$s
      n_roles <- 1 + isTRUE(xf_r$dated) + (!is.null(sf_r) && isTRUE(sf_r$dated))
      n_req_full <- 0
      v_keys <- character(M_all)
      for (ell in seq_len(M_all)) {
        v <- views[[ell]]
        ks <- match(c(v$cohorts$g, v$cohorts$K), st_r$cohorts)
        for (kk in ks) n_req_full <- n_req_full + n_roles * (st_r$cohort_end[kk] - st_r$cohort_start[kk] + 1)
        v_keys[ell] <- paste(
          fmt_design_num(v$cohorts$g),
          paste(fmt_design_num(v$cohorts$K), collapse = ","),
          paste(v$cohort_x_rows, collapse = ","),
          paste(v$cohort_s_rows, collapse = ","),
          paste(v$cohort_w_rows, collapse = ","),
          sep = "|"
        )
      }
      u_idx <- which(!duplicated(v_keys))
      u_map <- match(v_keys, v_keys[u_idx])
      tryCatch({
        ref_u <- build_chained_reference_reader(src_ref, cp04_reader, views[u_idx], length(u_idx),
          ctrl_eff$owner_control, recover_missing_reference = recover)
        ref_obj <- unclass(ref_u)
        attr(ref_obj, ".cdid_owner_seal") <- NULL
        ref_obj$view_recipe <- as.integer(ref_u$view_recipe[u_map])
        ref_obj$counts$n_views <- as.integer(M_all)
        ref_obj$counts$n_requirement_slots <- n_req_full
        ipw_owner_seal(structure(ref_obj, class = "cdid_reference_reader"), "reference_reader")
      }, cdid_reference_missing_error = function(e) {
        build_chained_reference_reader(src_ref, cp04_reader, views, M_all, ctrl_eff$owner_control,
        recover_missing_reference = recover)
      })
    } else {
      build_chained_reference_reader(src_ref, cp04_reader, views, M_all, ctrl_eff$owner_control,
        recover_missing_reference = recover)
    }
  } else NULL

  want_within <- identical(wt_mode, "within")
  res <- if (sel_mode %in% c("none", "ipw")) {
    evaluate_weighted_regimes_core(
      descs = descs, views = views, reader = cp04_reader, frame_owner = frame_owner,
      f_Z = f_Z, Y_mat = pm$Y_mat, O_mat = pm$O_mat, all_times = pm$all_times,
      cp04 = cp04, reference = reference, est_method = est_method,
      outcome_pooling = outcome_pooling, sel_mode = sel_mode,
      glm_control = ctrl_eff$glm_control, rank_tol = or_public_rank_tol(), ctrl = ctrl_eff,
      compute_if = compute_if, compute_within = want_within, nuisance_sink = nuisance_sink,
      center_if = TRUE
    )
  } else {
    evaluate_reconstruction_regimes_core(
      descs = descs, views = views, reader = cp04_reader, frame_owner = frame_owner,
      Y_mat = pm$Y_mat, O_mat = pm$O_mat, all_times = pm$all_times,
      cp04 = cp04, reference = reference, est_method = est_method,
      outcome_pooling = outcome_pooling, imputation_pooling = imputation_pooling,
      augmented = identical(sel_mode, "dr"),
      glm_control = ctrl_eff$glm_control, rank_tol = or_public_rank_tol(), ctrl = ctrl_eff,
      compute_if = compute_if, compute_within = want_within, nuisance_sink = nuisance_sink,
      center_if = TRUE
    )
  }

  mhat_all <- res$mhat
  N_orig <- cp04_reader$state$N
  n_a <- nrow(pm$Y_mat)
  max_e <- ctrl$gmm_bounds$max_matrix_entries
  Phi_F_all <- NULL
  mean_phi_all <- NULL
  if (isTRUE(compute_if)) {
    # One complete original-unit IF matrix, centered by the producer in place
    # (`x - colMeans(x)`, entrywise as `sweep()` computes it).
    Phi_F_all <- res$Phi_F
    mean_phi_all <- res$Phi_mean
    res$Phi_F <- NULL
  }
  U_log <- res$U_log
  res$U_log <- NULL

  # Post-fit exact-mask spanning-forest reduction (preserves 100% of pooled fits)
  rec_tag <- if (is.null(reference)) rep.int(1L, M_all) else as.integer(reference$view_recipe)
  grp_key <- paste(links_all$mask_key, rec_tag, sep = "|")
  sf_gram <- NULL
  max_m_diff <- 0.0
  if (isTRUE(reduce_structural_redundancy) && wt_mode %in% c("optimal", "within")) {
    use_phi <- identical(wt_mode, "optimal") && !is.null(Phi_F_all)
    live_w <- (if (use_phi) as.double(length(Phi_F_all)) else 0) +
      (if (want_within && !is.null(U_log)) U_log$entries else 0)
    b_w <- stream_block_rows(if (use_phi) N_orig else n_a, 7.0 * M_all, max_e, live_w)
    sf <- verify_post_fit_forest_reduction(
      s_vec = links_all$s, t_vec = links_all$t, group_key = grp_key,
      mhat_all = mhat_all,
      Phi_block = if (use_phi) function(r0, r1, cols) Phi_F_all[r0:r1, cols, drop = FALSE] else NULL,
      n_phi_rows = if (use_phi) nrow(Phi_F_all) else 0L,
      U_log = if (want_within) U_log else NULL,
      block_rows = b_w
    )
    kept_idx <- sf$kept_idx
    dropped_idx <- sf$dropped_idx
    max_m_diff <- sf$max_moment_diff
    if (length(kept_idx) < M_all) sf_gram <- sf$G
  } else {
    kept_idx <- seq_len(M_all)
    dropped_idx <- integer(0L)
  }

  L_kept <- length(kept_idx)
  mhat_used <- mhat_all[kept_idx]
  phi_scale <- as.double(n_a) / as.double(N_orig)
  if (isTRUE(compute_if)) {
    Phi_kept <- if (L_kept == M_all) Phi_F_all else Phi_F_all[, kept_idx, drop = FALSE]
    Phi_F_all <- NULL
    mean_kept <- mean_phi_all[kept_idx]
    mean_phi_used <- mean_kept * phi_scale
    # Complete centered, scaled original-unit rows with every cohort's column,
    # entrywise equal to `(Phi_F - colMeans(Phi_F))[, kept] * n_A / N`: every
    # replay scales the centered rows it reads.
    phi_rows <- phi_row_reader(Phi_kept, phi_scale)
  } else {
    Phi_kept <- NULL
    phi_rows <- NULL
    mean_phi_used <- numeric(L_kept)
  }

  cov_source <- list(
    wt_mode = wt_mode, L = L_kept, n_a = n_a, N_orig = N_orig,
    U_log = if (want_within) U_log else NULL,
    within_cols = kept_idx,
    within_g = links_all$g[kept_idx],
    within_scale = if (want_within) build_within_scale_by_g(
      g_vec = links_all$g[kept_idx],
      view_rec = rec_tag[kept_idx],
      reader = cp04_reader,
      frame_owner = frame_owner,
      cp04 = cp04,
      reference = reference
    ) else NULL,
    phi_rows = phi_rows,
    phi_live_entries = if (is.null(Phi_kept)) 0 else as.double(length(Phi_kept))
  )
  U_log <- NULL

  out <- list(
    mhat_used = mhat_used,
    V_used = NULL,
    phi_rows = phi_rows,
    cov_source = cov_source,
    mean_phi = mean_phi_used,
    kept_idx = kept_idx,
    dropped_idx = dropped_idx,
    sf_gram = sf_gram,
    max_moment_diff = max_m_diff,
    backend = res$backend
  )
  if (!isTRUE(defer_covariance)) {
    out$V_used <- cross_pair_general_covariance(cov_source, max_e)
  }
  out
}

make_centered_phi_rows <- function(Phi_kept, mean_kept, phi_scale) {
  for (j in seq_len(ncol(Phi_kept))) Phi_kept[, j] <- Phi_kept[, j] - mean_kept[j]
  phi_row_reader(Phi_kept, phi_scale)
}

phi_row_reader <- function(Phi_centered, phi_scale) {
  force(Phi_centered)
  force(phi_scale)
  function(rr) Phi_centered[rr, , drop = FALSE] * phi_scale
}

# Retained weighting covariance of the cross-pair general route from streamed
# sources. `within`: per-cohort Gram of scaled signed level loadings replayed by
# date and A-row block. `optimal`: `crossprod` of complete centered, scaled
# original-unit IF rows (all cohorts in each row, so shared controls carry
# cross-cohort covariance) divided by `n_A`. `identity`: zero matrix. The row
# block is chosen from the workspace left after the arrays live during the
# accumulation; at least one row is always admitted, so this pass never
# refuses a fit that the retained solve admits.
cross_pair_general_covariance <- function(cov_source, max_entries, block_rows = NULL) {
  L <- cov_source$L
  if (identical(cov_source$wt_mode, "within")) {
    lg <- cov_source$U_log
    if (is.null(block_rows)) {
      block_rows <- stream_block_rows(lg$n_a, 4.0 * L, max_entries,
        lg$entries + cov_source$phi_live_entries + 5.0 * as.double(L) * L)
    }
    return(within_cov_stream(lg, cov_source$within_cols, cov_source$within_g,
      cov_source$within_scale, block_rows))
  }
  if (identical(cov_source$wt_mode, "optimal")) {
    n_rows <- cov_source$N_orig
    if (is.null(block_rows)) {
      block_rows <- stream_block_rows(n_rows, 3.0 * L, max_entries,
        cov_source$phi_live_entries + 5.0 * as.double(L) * L)
    }
    S <- matrix(0.0, L, L)
    for (r0 in seq.int(1L, max(1L, n_rows), by = block_rows)) {
      r1 <- min(n_rows, r0 + block_rows - 1L)
      if (r1 < r0) next
      S <- S + crossprod(cov_source$phi_rows(r0:r1))
    }
    V_opt <- S / as.double(cov_source$n_a)
    return(0.5 * (V_opt + t(V_opt)))
  }
  matrix(0.0, L, L)
}

# Same-pair `gmm_weighting = "within"` helper for all CP04 routes (pooled OR,
# selection IPW, imputation, augmented, and reference inputs).
compute_same_pair_within_bundle <- function(design_links, point_links_all, panel_owner,
                                            frame_owner, cp04, cp04_reader, nuis_res,
                                            f_Z, ctrl, reduce_structural_redundancy = TRUE,
                                            within_log_limit = Inf) {
  M_all <- nrow(design_links)
  n_a <- as.integer(frame_owner$n_a_units)
  pm <- cross_pair_extract_panel_matrices(panel_owner, n_a)

  sel_mode <- if (!is.null(cp04$selection_est_method)) {
    cp04$selection_est_method
  } else if (identical(cp04$canonical_assumption, "conditional_mtar")) {
    "ipw"
  } else {
    "none"
  }
  est_method <- if (!is.null(cp04$imputation_est_method)) cp04$imputation_est_method else
    if (!is.null(cp04$augmented_est_method)) cp04$augmented_est_method else
      if (!is.null(cp04$dr_method)) "dr" else if (!is.null(cp04$or_method)) "reg" else "ipw"
  outcome_pooling <- if (is.null(cp04$outcome_pooling)) "pooled" else cp04$outcome_pooling
  imputation_pooling <- if (is.null(cp04$imputation_pooling)) "separate" else cp04$imputation_pooling

  # Include any supported auxiliary pretrend links in `nuis_res$cp04_points`
  # so pooled OR / pooled imputation slopes match `nuis_res$cp04` identically.
  cp_pts <- nuis_res$cp04_points
  M_union <- if (!is.null(cp_pts) && length(cp_pts) >= M_all) length(cp_pts) else M_all
  if (M_union > M_all) {
    aux_idx <- seq.int(M_all + 1L, M_union)
    aux_df <- data.frame(
      g = vapply(cp_pts[aux_idx], function(p) as.double(p$cohorts$g), numeric(1L)),
      s = vapply(cp_pts[aux_idx], function(p) as.double(p$pair$s), numeric(1L)),
      t = vapply(cp_pts[aux_idx], function(p) as.double(p$pair$t), numeric(1L)),
      control_pool = vapply(cp_pts[aux_idx], function(p) {
        paste0("{", paste(sprintf("%.17g", as.double(p$cohorts$K)), collapse = ","), "}")
      }, character(1L)),
      link_id = vapply(cp_pts[aux_idx], function(p) as.character(p$link_id), character(1L)),
      stringsAsFactors = FALSE
    )
    aux_df$functional_key <- paste0("cmpcat:CP04-aux:g", aux_df$g, ":K", aux_df$control_pool)
    cat_df <- data.frame(
      g = as.double(design_links$g),
      s = as.double(design_links$s),
      t = as.double(design_links$t),
      control_pool = as.character(design_links$control_pool),
      link_id = if (!is.null(design_links$link_id)) as.character(design_links$link_id) else paste0("cp04:catalog:", seq_len(M_all)),
      functional_key = as.character(design_links$functional_key),
      stringsAsFactors = FALSE
    )
    union_df <- rbind(cat_df, aux_df)
  } else {
    union_df <- design_links
  }

  built <- build_general_link_descriptors(
    links_df = union_df, all_times = pm$all_times, reader = cp04_reader,
    O_mat = pm$O_mat, sel_mode = sel_mode, is_cross_pair = FALSE
  )
  descs <- built$descs
  views <- built$views

  ctrl_eff <- ctrl
  ctrl_eff$owner_control <- cp04_reader$control

  active_b <- nuis_res$cp04
  reference <- if (!is.null(active_b) && !is.null(active_b$reference_ref)) {
    chained_backend_reference(active_b)
  } else if (!is.null(cp04$reference_source)) {
    build_chained_reference_reader(cp04$reference_source, cp04_reader, views, M_all, ctrl_eff$owner_control)
  } else NULL

  existing_sb <- if (!is.null(active_b) && isTRUE(active_b$uses_selection_models) &&
                     !is.null(active_b$selection_backend) &&
                     length(active_b$selection_backend$links) == length(views)) {
    active_b$selection_backend
  } else NULL

  res <- if (sel_mode %in% c("none", "ipw")) {
    evaluate_weighted_regimes_core(
      descs = descs, views = views, reader = cp04_reader, frame_owner = frame_owner,
      f_Z = f_Z, Y_mat = pm$Y_mat, O_mat = pm$O_mat, all_times = pm$all_times,
      cp04 = cp04, reference = reference, est_method = est_method,
      outcome_pooling = outcome_pooling, sel_mode = sel_mode,
      glm_control = ctrl_eff$glm_control, rank_tol = or_public_rank_tol(), ctrl = ctrl_eff,
      compute_if = FALSE, compute_within = TRUE, existing_sel_fit = existing_sb,
      within_log_limit = within_log_limit
    )
  } else {
    evaluate_reconstruction_regimes_core(
      descs = descs, views = views, reader = cp04_reader, frame_owner = frame_owner,
      Y_mat = pm$Y_mat, O_mat = pm$O_mat, all_times = pm$all_times,
      cp04 = cp04, reference = reference, est_method = est_method,
      outcome_pooling = outcome_pooling, imputation_pooling = imputation_pooling,
      augmented = identical(sel_mode, "dr"),
      glm_control = ctrl_eff$glm_control, rank_tol = or_public_rank_tol(), ctrl = ctrl_eff,
      compute_if = FALSE, compute_within = TRUE, existing_sel_fit = existing_sb,
      within_log_limit = within_log_limit
    )
  }

  # Signed loadings stay in the replay log; only catalog columns are read.
  U_log <- res$U_log
  res$U_log <- NULL

  # Build exact mask key on A per link via per-cohort BLAS pair-mask IDs
  a_g <- as.double(frame_owner$a_g)
  w_a <- as.double(frame_owner$w_a)
  u_ag <- sort(unique(a_g))
  mid_mat_by_g <- list()
  for (gv in u_ag) {
    idx_k <- which(a_g == gv & w_a > 0)
    mid_mat <- matrix(0L, nrow = pm$n_cal, ncol = pm$n_cal)
    if (length(idx_k) > 0L) {
      pm_k <- compute_cohort_pair_masks(pm$O_mat[idx_k, , drop = FALSE], idx_k, pm$n_cal)
      if (length(pm_k$cs) > 0L) {
        mid_mat[cbind(pm_k$cs, pm_k$ct)] <- pm_k$mask_id
        mid_mat[cbind(pm_k$ct, pm_k$cs)] <- pm_k$mask_id
      }
    }
    mid_mat_by_g[[as.character(gv)]] <- mid_mat
  }

  mask_keys <- character(M_all)
  for (m in seq_len(M_all)) {
    d <- descs[[m]]
    cohs_m <- c(d$g, d$K)
    mids_m <- integer(length(cohs_m))
    for (hh in seq_along(cohs_m)) {
      mm <- mid_mat_by_g[[as.character(cohs_m[hh])]]
      mids_m[hh] <- if (is.null(mm)) 0L else mm[d$cs, d$ct]
    }
    rec_m <- if (is.null(reference)) 1L else reference$view_recipe[m]
    mask_keys[m] <- paste(
      fmt_design_num(d$g),
      design_links$control_pool[m],
      paste(mids_m, collapse = ","),
      rec_m,
      sep = "|"
    )
  }

  all_est <- vapply(point_links_all, `[[`, numeric(1L), "estimate")
  sf_gram <- NULL
  if (isTRUE(reduce_structural_redundancy)) {
    max_e <- ctrl$gmm_bounds$max_matrix_entries
    sf <- verify_post_fit_forest_reduction(
      s_vec = design_links$s, t_vec = design_links$t, group_key = mask_keys,
      mhat_all = all_est, Phi_all = NULL, U_log = U_log,
      block_rows = stream_block_rows(n_a, 7.0 * M_all, if (is.null(max_e)) Inf else max_e, U_log$entries)
    )
    kept_idx <- sf$kept_idx
    dropped_idx <- sf$dropped_idx
    if (length(kept_idx) < M_all) sf_gram <- sf$G
  } else {
    kept_idx <- seq_len(M_all)
    dropped_idx <- integer(0L)
  }

  rec_cat <- if (is.null(reference)) rep.int(1L, M_all) else as.integer(reference$view_recipe[seq_len(M_all)])
  sc_by_g <- build_within_scale_by_g(
    g_vec = design_links$g,
    view_rec = rec_cat,
    reader = cp04_reader,
    frame_owner = frame_owner,
    cp04 = cp04,
    reference = reference
  )

  # The scale is applied to replayed rows when the retained Gram is formed.
  list(
    U_log = U_log,
    scale_by_g = sc_by_g,
    kept_idx = kept_idx,
    dropped_idx = dropped_idx,
    G = sf_gram
  )
}

cross_pair_public_diagnostics <- function(res, control, public_object = NULL) {
  b <- res$nuisance_registry$cp04
  sb_sum <- if (!is.null(b$reader_ref)) selection_backend_summary(b) else NULL
  list(
    status = "available",
    method = b$mode,
    est_method = b$est_method,
    outcome_pooling = b$outcome_pooling,
    imputation_pooling = b$imputation_pooling,
    selection = sb_sum,
    frame_counts = list(
      n_full_frame_units = res$reporting$inference$n_total,
      n_pair_eligible_A_units = res$frame_owner$n_a_units
    ),
    accounting = b$accounting,
    auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
    interpretation = "Cross-pair single-cohort nuisance fits and exact full-roster influence-function propagation."
  )
}
