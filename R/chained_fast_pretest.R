# Private native fast-path wrappers for the chained conditional pretrend test
# (`conditional_did_pretest`). All routines dispatch through `cdid_fast_kernel()`
# and fall back cleanly to the pure-R reference implementations whenever a mock
# or probe DLL omits optional `C_cdid_cpt_*` symbols or internal functions are
# traced/replaced.

cdid_fast_cpt_distinct_grid <- function(Z) {
  fk <- cdid_fast_kernel("C_cdid_cpt_distinct_grid")
  if (is.null(fk)) return(NULL)
  fk(Z)
}

cdid_fast_cpt_acc_init <- function(spec, grid, n_clusters, n_total, keep_rows = FALSE) {
  fk_init <- cdid_fast_kernel("C_cdid_cpt_acc_init")
  fk_add <- cdid_fast_kernel("C_cdid_cpt_acc_add_piece")
  fk_fin <- cdid_fast_kernel("C_cdid_cpt_acc_finalize")
  if (is.null(fk_init) || is.null(fk_add) || is.null(fk_fin)) return(NULL)
  dims <- as.integer(c(
    grid$n_distinct,
    length(spec$s),
    length(spec$coords),
    n_clusters,
    n_total,
    if (isTRUE(keep_rows)) 1L else 0L
  ))
  ptr <- fk_init(dims)
  if (is.null(ptr)) return(NULL)
  list(.cpt_c_acc = ptr, cmap = NULL, keep_rows = isTRUE(keep_rows))
}

cdid_fast_cpt_cell_add <- function(acc, spec, piece, grid, cmap) {
  fk <- cdid_fast_kernel("C_cdid_cpt_acc_add_piece")
  if (is.null(fk)) return(NULL)
  rows <- piece$rows
  if (!length(rows)) return(acc)
  rows_int <- as.integer(rows)
  if (is.null(acc$cmap)) acc$cmap <- as.integer(cmap)
  ok <- fk(
    acc$.cpt_c_acc,
    rows_int,
    piece$nu,
    piece$rho,
    piece$P,
    piece$xi,
    piece$blam,
    as.double(spec$s / spec$D),
    as.integer(grid$map),
    acc$cmap
  )
  if (identical(ok, FALSE)) {
    cpt_abort(
      paste0("Process extractor returned malformed or nonfinite rows for cell ", spec$id, "."),
      "cdid_conditional_pretest_integrity_error"
    )
  }
  if (!isTRUE(ok)) return(NULL)
  acc
}

cdid_fast_cpt_finalize_ptr <- function(fk_fin, ptr, spec, grid, cmap_int, N) {
  spec_num <- list(
    as.double(spec$sigma),
    as.double(spec$s),
    as.double(spec$D),
    as.integer(grid$M),
    as.double(N)
  )
  res <- fk_fin(
    ptr,
    grid$z,
    as.integer(grid$mult),
    as.integer(grid$map),
    cmap_int,
    spec_num
  )
  if (is.null(res)) return(NULL)
  cell <- list(
    id = spec$id,
    sigma = spec$sigma,
    J = res$J,
    mult = res$mult,
    M = res$M,
    zd = grid$z,
    Rmat = res$Rmat,
    Lmat = res$Lmat,
    AggX = res$AggX,
    AggK = res$AggK,
    A_k = res$A_k,
    A_d = res$A_d,
    A_v = res$A_v,
    mu1 = res$mu1,
    L1 = res$L1,
    R1 = res$R1,
    denominator_check = list(supplied = spec$D, accumulated = res$accumulated),
    coords = spec$coords,
    n_coordinates = length(spec$coords),
    n_terms = length(spec$s)
  )
  if (!is.null(res$rows)) cell$rows <- res$rows
  cell
}

cdid_fast_cpt_cell_finalize <- function(acc, spec, grid, N) {
  fk_fin <- cdid_fast_kernel("C_cdid_cpt_acc_finalize")
  if (is.null(fk_fin)) return(NULL)
  cmap_int <- if (is.null(acc$cmap)) rep.int(1L, as.integer(N)) else acc$cmap
  cdid_fast_cpt_finalize_ptr(fk_fin, acc$.cpt_c_acc, spec, grid, cmap_int, N)
}

cdid_fast_cpt_batch_cvmb <- function(E, cells, N) {
  fk <- cdid_fast_kernel("C_cdid_cpt_batch_cvmb")
  if (is.null(fk)) return(NULL)
  fk(E, cells, as.double(N))
}

cdid_fast_cpt_process_nonzero <- function(cells, n_clusters) {
  fk <- cdid_fast_kernel("C_cdid_cpt_process_nonzero")
  if (is.null(fk)) return(NULL)
  fk(cells, as.integer(n_clusters))
}

cdid_fast_cpt_pair_key <- function(pair) {
  paste0(sprintf("%a", as.double(pair$s) + 0), "|", sprintf("%a", as.double(pair$t) + 0))
}

cdid_fast_cpt_build_cells <- function(state, cells, ex, cache, keys, cmap, C, N, chunk_rows, keep_rows) {
  fk_init <- cdid_fast_kernel("C_cdid_cpt_acc_init")
  fk_fast <- cdid_fast_kernel("C_cdid_cpt_acc_add_fast")
  fk_fin <- cdid_fast_kernel("C_cdid_cpt_acc_finalize")
  if (is.null(fk_init) || is.null(fk_fast) || is.null(fk_fin)) return(NULL)

  n_cells <- length(cells)
  if (n_cells < 1L) return(list())
  cmap_int <- as.integer(cmap)
  keep_int <- if (isTRUE(keep_rows)) 1L else 0L

  # 1. Legacy CP02 risk-set IPW (A frame)
  if (is.null(state$cp04_reader)) {
    frame <- state$frame_owner
    owner_ctrl <- state$control$owner_control
    nA <- as.integer(frame$n_a_units)
    scale <- N / nA
    B <- as.integer(min(nA, chunk_rows))
    starts <- seq.int(1L, nA, by = B)
    built <- vector("list", n_cells)

    for (i in seq_len(n_cells)) {
      pl <- state$auxiliary$point_links[[cells[[i]]$q]]
      if (is.null(pl)) return(NULL)
      lv <- get("link_view", envir = pl$link_view_ref, inherits = FALSE)
      r <- as.double(pl$masses$rates)
      q <- as.double(pl$masses$q)
      m <- length(r)
      beta <- as.double(pl$coefficients)
      p <- length(beta)
      t_Abr <- t(matrix(as.double(pl$partials$A_beta_r), p, m))
      D_T <- as.double(pl$means$sum_a) / N
      D_C <- as.double(pl$means$sum_b) / N
      r_floor <- ipw_clip_floor(r)
      r_log_d <- ifelse(r > ipw_clip_eps, -1 / r, 0)

      grid <- cache[[keys[i]]]$distinct
      gmap_i <- as.integer(grid$map)
      spec_i <- ex[[i]]$spec
      ptr <- fk_init(as.integer(c(grid$n_distinct, 2L, m + p, C, N, keep_int)))
      if (is.null(ptr)) return(NULL)
      prm <- as.double(c(D_T, D_C, m, p, 0))

      for (k in seq_along(starts)) {
        a_s <- starts[k]
        a_e <- min(nA, a_s + B - 1L)
        ch <- query_ipw_view_chunk(lv, a_s, a_e, owner_ctrl, read_outcomes = TRUE, include_design = TRUE)
        n <- a_e - a_s + 1L
        rows <- as.integer(state$a_to_full[seq.int(a_s, a_e)])
        w <- as.double(ch$weight)
        gi <- as.integer(ch$cohort_index)
        risk <- which(ch$in_risk)
        pos <- as.integer(ch$observed_risk$local_pos)
        gj <- as.integer(ch$observed_risk$cohort_index)
        S <- numeric(n)
        S[pos] <- 1
        xr <- matrix(0, n, m)
        if (length(risk)) xr[cbind(risk, gi[risk])] <- w[risk] * (S[risk] - r[gi[risk]]) / q[gi[risk]]
        score <- matrix(0, n, p)
        m_pos <- length(pos)
        if (m_pos) {
          Z <- ch$design$Z_obs
          lv_ <- ipw_logit_values(Z, beta)
          fw <- w[pos] / r_floor[gj]
          D_obs <- as.double(gj == 1L)
          score[pos, ] <- (fw * (D_obs - lv_$prob)) * Z
          odds <- lv_$odds
          act <- lv_$active
          dy <- as.double(ch$observed_risk$y_right_obs) - as.double(ch$observed_risk$y_left_obs)
          lam_up <- matrix(0, m_pos, m)
          lam_up[cbind(seq_len(m_pos), gj)] <- r_log_d[gj]
        } else {
          Z <- matrix(0, 0L, p)
          D_obs <- fw <- odds <- act <- dy <- numeric(0L)
          lam_up <- matrix(0, 0L, m)
        }
        xb <- t(cp04_solve(pl$factors$H_chol, t(score + xr %*% t_Abr)))
        ok <- fk_fast(
          ptr, 1L, rows, pos,
          list(D_obs, fw, odds, act, dy),
          list(lam_up, Z, NULL),
          list(scale * xr, scale * xb),
          prm, gmap_i, cmap_int
        )
        if (identical(ok, FALSE)) {
          cpt_abort(
            paste0("Process extractor returned malformed or nonfinite rows for cell ", spec_i$id, "."),
            "cdid_conditional_pretest_integrity_error"
          )
        }
        if (!isTRUE(ok)) return(NULL)
      }

      cl_out <- cdid_fast_cpt_finalize_ptr(fk_fin, ptr, spec_i, grid, cmap_int, N)
      if (is.null(cl_out)) return(NULL)
      cl_out$estimate_scalar <- cells[[i]]$sign * ex[[i]]$estimate
      cl_out$route <- ex[[i]]$route
      cl_out$recipe_key <- keys[i]
      built[[i]] <- cl_out
    }
    return(built)
  }

  # 2. Full-roster CP04 / OR / DR backends
  reg <- state$registry
  backend <- reg$cp04
  kind <- cpt_backend_kind(backend)
  if (!kind %in% c("reference_ipw", "ipw_legacy", "ipw_selection", "dr_baseline", "dr_selection",
                   "or_baseline", "or_legacy", "or_selection")) {
    return(NULL)
  }

  js <- as.integer(reg$cp04_aux_index[vapply(cells, `[[`, integer(1L), "q")])
  u_js <- unique(js)
  B_det <- tryCatch(cpt_details_rows(backend, kind, u_js, state$retained_now), error = function(e) 0L)
  if (B_det < 1L) return(NULL)

  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  ref <- chained_backend_reference(backend)
  B <- as.integer(min(st$N, chunk_rows, B_det))
  starts <- seq.int(1L, st$N, by = B)

  ptrs <- vector("list", n_cells)
  gmaps <- vector("list", n_cells)
  for (i in seq_len(n_cells)) {
    grid <- cache[[keys[i]]]$distinct
    gmaps[[i]] <- as.integer(grid$map)
    spec_i <- ex[[i]]$spec
    ptrs[[i]] <- fk_init(as.integer(c(
      grid$n_distinct,
      length(spec_i$s),
      length(spec_i$coords),
      C,
      N,
      keep_int
    )))
    if (is.null(ptrs[[i]])) return(NULL)
  }

  # 2a. Two-arm IPW & DR backends
  if (kind %in% c("reference_ipw", "ipw_legacy", "ipw_selection", "dr_baseline", "dr_selection")) {
    dr <- kind %in% c("dr_baseline", "dr_selection")
    icpt <- dr && isTRUE(backend$outcome_recipe$intercept)
    separate <- dr && identical(backend$outcome_pooling, "separate")
    trt_cols <- if (dr) backend$treatment_recipe$index else NULL
    st_Z <- st$Z
    st_a <- st$a
    d_Z <- if (identical(kind, "ipw_legacy")) ncol(st_Z) else 0L
    sel_fit <- if (kind %in% c("ipw_selection", "dr_selection")) backend$selection_backend else NULL
    up_map <- if (kind %in% c("ipw_selection", "dr_selection")) backend$up_map else NULL

    tjs <- vapply(js, function(j) as.integer(backend$links[[j]]$model), integer(1L))
    u_tjs <- unique(tjs)

    # Precompute unique pair descriptors and per-tj metadata
    pair_keys <- character(0L)
    pair_objs <- list()
    pair_cohs <- list()
    tj_meta <- vector("list", length(backend$treatments))

    # Register helper to add a pair+cohorts to the chunk-level pair table
    add_pair_entry <- function(pair_obj, cohs) {
      pk <- cdid_fast_cpt_pair_key(pair_obj)
      p_idx <- match(pk, pair_keys)
      if (is.na(p_idx)) {
        p_idx <- length(pair_keys) + 1L
        pair_keys[p_idx] <<- pk
        pair_objs[[p_idx]] <<- pair_obj
        pair_cohs[[p_idx]] <<- as.integer(cohs)
      } else {
        pair_cohs[[p_idx]] <<- unique(c(pair_cohs[[p_idx]], as.integer(cohs)))
      }
      p_idx
    }

    for (tj in u_tjs) {
      spec_t <- backend$treatment_specs[[tj]]
      f <- backend$treatments[[tj]]
      p_idx <- add_pair_entry(spec_t$pair, spec_t$cohorts)

      up_n <- switch(
        kind,
        reference_ipw = ,
        dr_baseline = length(spec_t$cohorts),
        ipw_legacy = length(spec_t$cohorts) * (d_Z + 1L),
        length(f$upstream_ids)
      )
      p_trt <- length(f$coefficients)
      outcome <- if (dr) {
        if (separate) or_outcome_block_fit(backend$outcome_groups, spec_t$regression) else backend$shared
      } else NULL
      p_or <- if (dr) outcome$p else 0L
      kept <- if (dr && p_or > 0L) backend$or_frame_index[backend$outcome_recipe$slope_index[outcome$keep]] else integer(0L)
      reg_blk <- if (dr) backend$regressions[[spec_t$regression]] else NULL
      rate_vec <- if (kind %in% c("reference_ipw", "dr_baseline")) {
        vapply(spec_t$support_ids, function(s) backend$supports[[s]]$rate, numeric(1L))
      } else NULL
      v_vec <- if (identical(kind, "dr_baseline")) {
        vapply(spec_t$support_ids, function(s) backend$supports[[s]]$v, numeric(1L))
      } else NULL
      chan_pos <- if (kind %in% c("ipw_selection", "dr_selection")) {
        lapply(spec_t$selection_ids, function(qk) match(up_map$pair_channels[[qk]], f$upstream_ids))
      } else NULL

      tj_meta[[tj]] <- list(
        spec = spec_t,
        f = f,
        pair_idx = p_idx,
        up_n = up_n,
        p_trt = p_trt,
        outcome = outcome,
        p_or = p_or,
        kept = kept,
        reg_blk = reg_blk,
        reg_id = if (dr) spec_t$regression else NA_integer_,
        reg_chr = if (dr && separate) as.character(spec_t$regression) else NULL,
        rate = rate_vec,
        rate_floor = if (!is.null(rate_vec)) ipw_clip_floor(rate_vec) else NULL,
        rate_log_d = if (!is.null(rate_vec)) ifelse(rate_vec > ipw_clip_eps, -1 / rate_vec, 0) else NULL,
        v = v_vec,
        chan_pos = chan_pos,
        t_cross = t(f$cross),
        prm = as.double(c(
          unname(f$means["dp"]),
          unname(f$means["dm"]),
          up_n,
          p_trt,
          if (icpt) 1 else 0
        ))
      )
    }

    sel_ok <- !is.null(sel_fit) && cdid_are_default_fns(c("selection_prediction_operation", "predict_pair_selection"))
    fuse_det <- is.null(ref) && switch(
      kind,
      ipw_selection = sel_ok && cdid_are_default_fns(c(
        "cpt_backend_details", "evaluate_cp04_if_chunk", "evaluate_cp04_selection_if_chunk",
        "cp04_selection_if_budget", "validate_cp04_backend"
      )),
      ipw_legacy = cdid_are_default_fns(c(
        "cpt_backend_details", "evaluate_cp04_if_chunk", "cp04_if_budget",
        "validate_cp04_backend", "cp04_predict_selection"
      )),
      dr_selection = separate && sel_ok && cdid_are_default_fns(c(
        "cpt_backend_details", "evaluate_dr_corollary1_if_chunk", "evaluate_dr_separate_if_chunk",
        "dr_separate_if_envelope", "or_selection_separate_if_core", "validate_dr_backend", "cp04_budget"
      )),
      dr_baseline = separate && cdid_are_default_fns(c(
        "cpt_backend_details", "evaluate_dr_baseline_if_chunk", "evaluate_dr_separate_if_chunk",
        "dr_separate_if_envelope", "or_baseline_separate_if_core", "validate_dr_backend", "cp04_budget"
      )),
      FALSE
    )

    gids_need <- sids_need <- uids_need <- NULL
    sup_pair_idx <- sel_pair_idx <- NULL
    blk_meta <- NULL

    if (fuse_det) {
      if (kind %in% c("ipw_selection", "ipw_legacy")) {
        validate_cp04_backend(backend)
      } else {
        validate_dr_backend(backend)
      }
      if (identical(kind, "ipw_selection")) {
        uids_need <- sort(unique(unlist(lapply(backend$treatments[u_tjs], `[[`, "upstream_ids"))))
        for (tj in u_tjs) tj_meta[[tj]]$up_pos <- match(tj_meta[[tj]]$f$upstream_ids, uids_need)
      } else if (identical(kind, "ipw_legacy")) {
        sids_need <- sort(unique(unlist(lapply(backend$treatment_specs[u_tjs], `[[`, "selection_ids"))))
        sel_pair_idx <- integer(length(backend$selections))
        for (sj in sids_need) {
          s_sp <- backend$selection_specs[[sj]]
          sel_pair_idx[sj] <- add_pair_entry(s_sp$pair, s_sp$cohort)
        }
      } else {
        gids_need <- sort(unique(c(
          vapply(backend$links[u_js], `[[`, integer(1L), "regression"),
          vapply(backend$treatment_specs[u_tjs], `[[`, integer(1L), "regression")
        )))
        sids_need <- sort(unique(c(
          vapply(backend$links[u_js], `[[`, integer(1L), "support_id"),
          unlist(lapply(backend$regression_specs[gids_need], `[[`, "support_ids")),
          unlist(lapply(backend$treatment_specs[u_tjs], `[[`, "support_ids"))
        )))
        uids_need <- if (identical(kind, "dr_selection")) {
          sort(unique(c(
            unlist(lapply(backend$supports[sids_need], `[[`, "channels")),
            unlist(lapply(backend$treatments[u_tjs], `[[`, "upstream_ids"))
          )))
        } else sids_need
        sup_pair_idx <- integer(length(backend$supports))
        for (sid in sids_need) {
          sp_s <- backend$support_specs[[sid]]
          sup_pair_idx[sid] <- add_pair_entry(sp_s$pair, sp_s$cohort)
        }
        blk_meta <- vector("list", length(gids_need))
        for (kb in seq_along(gids_need)) {
          b <- gids_need[kb]
          of <- or_outcome_block_fit(backend$outcome_groups, b)
          rec <- backend$regressions[[b]]
          spos <- match(backend$regression_specs[[b]]$support_ids, sids_need)
          blk_meta[[kb]] <- list(
            of = of,
            rec = rec,
            spos = spos,
            upos = if (identical(kind, "dr_selection")) match(rec$channel_ids, uids_need) else spos,
            kept_or = if (of$p > 0L) backend$or_frame_index[backend$outcome_recipe$slope_index[of$keep]] else integer(0L),
            t_forcing_beta = if (of$p > 0L) t(rec$forcing_beta) else NULL
          )
        }
        for (tj in u_tjs) {
          tm <- tj_meta[[tj]]
          tj_meta[[tj]]$up_pos <- match(
            if (identical(kind, "dr_selection")) tm$f$upstream_ids else tm$spec$support_ids,
            uids_need
          )
          tj_meta[[tj]]$reg_pos <- match(tm$reg_id, gids_need)
        }
      }
      if (!is.null(sel_fit) && identical(sel_fit$selection_model, "pair_logit")) {
        act_blks <- which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids_need), logical(1L)))
        for (m_b in act_blks) {
          for (cid in sel_fit$blocks[[m_b]]$cell_ids) {
            cl_s <- sel_fit$cells[[cid]]
            add_pair_entry(cl_s$pair, cl_s$cohort)
          }
        }
      }
    }

    n_pairs <- length(pair_keys)
    n_sel_max <- if (identical(kind, "ipw_legacy")) {
      length(backend$selections)
    } else if (kind %in% c("ipw_selection", "dr_selection")) {
      length(sel_fit$pair_cohorts)
    } else 0L

    pair_env <- new.env(parent = emptyenv())
    pair_env$offset <- 0L
    pair_env$pair_S <- NULL
    chunk_pair_reader <- if (fuse_det && !is.null(sel_fit) && identical(sel_fit$selection_model, "pair_logit")) {
      function(r_sub, pair_obj, read_outcomes = FALSE) {
        if (!isTRUE(read_outcomes) && !is.null(pair_env$pair_S)) {
          p_idx <- match(cdid_fast_cpt_pair_key(pair_obj), pair_keys)
          if (!is.na(p_idx)) {
            return(list(S = as.double(pair_env$pair_S[[p_idx]][r_sub - pair_env$offset])))
          }
        }
        cp04_pair_rows(reader, r_sub, pair_obj, isTRUE(read_outcomes))
      }
    } else NULL

    sel_op <- if (sel_ok) {
      selection_prediction_operation(
        sel_fit,
        source_reader = reader,
        source_reference = ref,
        pair_response_reader = chunk_pair_reader
      )
    } else NULL
    on.exit(if (!is.null(sel_op)) sel_op$abort(), add = TRUE)

    for (k_ch in seq_along(starts)) {
      rows <- seq.int(starts[k_ch], min(st$N, starts[k_ch] + B - 1L))
      n_r <- length(rows)
      ci <- st$cohort_index[rows]

      if (fuse_det) {
        if (identical(kind, "ipw_selection")) {
          cp04_selection_if_budget(backend, n_r, u_js, TRUE)
        } else if (identical(kind, "ipw_legacy")) {
          cp04_if_budget(backend, n_r, u_js)
        } else {
          pre_chk <- dr_separate_if_envelope(backend, n_r, u_js, TRUE)
          cp04_budget(
            backend$control,
            c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
            pre_chk$scratch,
            c(backend$accounting$work_ops_bound, pre_chk$work),
            pre_chk$output,
            "separate DR IF chunk"
          )
        }
      }

      # Evaluate each unique pair at most once on this chunk
      pair_S <- vector("list", n_pairs)
      pair_dy <- vector("list", n_pairs)
      for (p_idx in seq_len(n_pairs)) {
        cand <- which(ci %in% pair_cohs[[p_idx]])
        if (identical(kind, "reference_ipw")) cand <- cand[st$a_index[rows[cand]] > 0L]
        if (identical(kind, "ipw_legacy")) cand <- cand[st$U[rows[cand]]]
        S_vec <- logical(n_r)
        dy_vec <- numeric(n_r)
        if (length(cand)) {
          hit_p <- cp04_pair_rows(reader, rows[cand], pair_objs[[p_idx]], TRUE)
          ob_p <- which(hit_p$S == 1)
          if (length(ob_p)) {
            obs_c <- cand[ob_p]
            S_vec[obs_c] <- TRUE
            dy_vec[obs_c] <- as.double(hit_p$dy[ob_p])
          }
        }
        pair_S[[p_idx]] <- S_vec
        pair_dy[[p_idx]] <- dy_vec
      }
      pair_env$offset <- starts[k_ch] - 1L
      pair_env$pair_S <- pair_S

      sel_cache <- if (n_sel_max > 0L) vector("list", n_sel_max) else NULL
      det <- NULL
      up_chunk <- NULL
      a_if_chunk <- NULL
      alpha_xi <- NULL
      prob_vec <- NULL
      blk_ifs <- NULL
      w_ch <- if (fuse_det) st$w[rows] else NULL

      if (!fuse_det) {
        det <- cpt_backend_details(backend, kind, rows, u_js)
      } else if (identical(kind, "ipw_selection")) {
        up_chunk <- sel_op$upstream(up_map, rows, uids_need, FALSE)$up
      } else if (identical(kind, "ipw_legacy")) {
        A_ch <- as.double(st$a_index[rows] > 0L)
        U_ch <- st$U[rows]
        a_if_chunk <- matrix(0, n_r, length(st$cohorts))
        for (k_c in unique(ci)) {
          rr_k <- which(ci == k_c & U_ch)
          if (length(rr_k)) a_if_chunk[rr_k, k_c] <- w_ch[rr_k] * (A_ch[rr_k] - st_a[k_c]) / st$mass_U[k_c]
        }
        alpha_xi <- prob_vec <- vector("list", length(backend$selections))
        for (sj in sids_need) {
          s_sp <- backend$selection_specs[[sj]]
          f_s <- backend$selections[[sj]]
          rr_s <- which(ci == s_sp$cohort & U_ch)
          xi_s <- matrix(0, n_r, d_Z)
          pr_s <- numeric(n_r)
          if (length(rr_s)) {
            pr_s[rr_s] <- cp04_predict_selection(reader, rows[rr_s], s_sp, f_s)
            S_s <- as.double(pair_S[[sel_pair_idx[sj]]][rr_s])
            sc_s <- w_ch[rr_s] * (S_s - pr_s[rr_s]) * st_Z[rows[rr_s], , drop = FALSE]
            xi_s[rr_s, ] <- t(cp04_solve(f_s$H_chol, t(sc_s)))
          }
          alpha_xi[[sj]] <- xi_s
          prob_vec[[sj]] <- pr_s
        }
      } else {
        # Separate DR (`dr_selection` or `dr_baseline`)
        n_s <- length(sids_need)
        sup_hits <- vector("list", n_s)
        if (identical(kind, "dr_selection")) {
          up_chunk <- sel_op$upstream(up_map, rows, uids_need, FALSE)$up
          for (j_s in seq_len(n_s)) {
            sid <- sids_need[j_s]
            sp_s <- backend$support_specs[[sid]]
            p_idx <- sup_pair_idx[sid]
            obs_s <- which(pair_S[[p_idx]] & (ci == sp_s$cohort))
            if (length(obs_s)) {
              r_obs <- rows[obs_s]
              qk <- backend$q_of_s[sid]
              sc_q <- sel_cache[[qk]]
              if (is.null(sc_q)) {
                pr_qk <- sel_op$pair(r_obs, qk)
                ld_qk <- sel_op$deriv(up_map, qk, r_obs, pred = pr_qk, Z = st_Z, a_vec = st_a)
                sc_q <- list(
                  v = or_corollary1_factor(st_a[sp_s$cohort], st$rho[r_obs], pr_qk$prob),
                  fw = ld_qk$pred$weight,
                  lam = ld_qk$log_deriv
                )
                sel_cache[[qk]] <- sc_q
              }
              sup_hits[[j_s]] <- list(obs = obs_s, dy = pair_dy[[p_idx]][obs_s], v = sc_q$v)
            } else {
              sup_hits[[j_s]] <- list(obs = integer(0L), dy = numeric(0L), v = numeric(0L))
            }
          }
        } else {
          A_ch <- as.double(st$a_index[rows] > 0L)
          up_chunk <- matrix(0, n_r, n_s)
          for (j_s in seq_len(n_s)) {
            sid <- sids_need[j_s]
            sp_s <- backend$support_specs[[sid]]
            u_s <- backend$supports[[sid]]
            p_idx <- sup_pair_idx[sid]
            ix_s <- which(ci == sp_s$cohort)
            if (length(ix_s)) {
              S_ix <- as.double(pair_S[[p_idx]][ix_s])
              up_chunk[ix_s, j_s] <- w_ch[ix_s] * (S_ix - A_ch[ix_s] * u_s$rate) / u_s$massA
              obs_s <- ix_s[S_ix == 1]
            } else {
              obs_s <- integer(0L)
            }
            sup_hits[[j_s]] <- list(
              obs = obs_s,
              dy = pair_dy[[p_idx]][obs_s],
              v = rep.int(u_s$v, length(obs_s))
            )
          }
        }

        blk_ifs <- vector("list", length(gids_need))
        N_grp <- backend$outcome_groups$N
        for (kb in seq_along(gids_need)) {
          bm <- blk_meta[[kb]]
          f_or <- bm$of
          rec <- bm$rec
          ha <- numeric(n_r)
          hb <- matrix(0, n_r, f_or$p)
          for (s_pos in bm$spos) {
            h <- sup_hits[[s_pos]]
            ix <- h$obs
            if (length(ix)) {
              d_b <- if (f_or$p > 0L) {
                f_or$L * frame$Z_a[st$a_index[rows[ix]], bm$kept_or, drop = FALSE]
              } else {
                matrix(0, length(ix), 0L)
              }
              q_b <- w_ch[ix] * h$v
              r_b <- h$dy - f_or$alpha - (if (f_or$p > 0L) as.vector(d_b %*% f_or$beta) else 0)
              qr_b <- q_b * r_b
              ha[ix] <- ha[ix] + qr_b
              if (f_or$p > 0L) {
                hb[ix, ] <- hb[ix, , drop = FALSE] + qr_b * (if (icpt) sweep(d_b, 2L, f_or$mean_d) else d_b)
              }
            }
          }
          us_b <- up_chunk[, bm$upos, drop = FALSE]
          xb_b <- hb
          if (f_or$p > 0L) {
            rhs <- sweep(hb, 2L, rec$score_mean_beta) + us_b %*% bm$t_forcing_beta
            xb_b <- t(if (icpt) {
              or_pooled_wls_gram_solve(f_or$fit, t(rhs), normalizer = N_grp)
            } else {
              or_wls_gram_solve(f_or$fit, t(rhs), normalizer = N_grp)
            })
          }
          xa_b <- if (icpt) {
            (ha - rec$score_mean_alpha + as.vector(us_b %*% rec$forcing_alpha)) / rec$mass -
              (if (f_or$p > 0L) as.vector(xb_b %*% rec$mean_d) else 0)
          } else NULL
          blk_ifs[[kb]] <- list(alpha = xa_b, beta = xb_b)
        }
      }

      rd_cache <- vector("list", length(backend$treatments))

      for (tj in u_tjs) {
        tm <- tj_meta[[tj]]
        spec_t <- tm$spec
        f <- tm$f
        S_vec <- pair_S[[tm$pair_idx]]
        rr <- which(S_vec & (ci %in% spec_t$cohorts))
        m_obs <- length(rr)
        dy <- pair_dy[[tm$pair_idx]][rr]
        kk <- match(ci[rr], spec_t$cohorts)
        rows_rr <- rows[rr]
        D_obs <- as.double(kk == 1L)
        fw <- numeric(m_obs)
        fw_score <- if (fuse_det && identical(kind, "dr_selection")) numeric(m_obs) else NULL
        lam_up <- matrix(0, m_obs, tm$up_n)

        if (kind %in% c("reference_ipw", "dr_baseline")) {
          if (m_obs) {
            lam_up[cbind(seq_len(m_obs), kk)] <- tm$rate_log_d[kk]
            if (identical(kind, "reference_ipw")) {
              W <- chained_reference_weight(ref, spec_t$reference, rows_rr)
              fw <- W / tm$rate_floor[kk]
            } else {
              wr <- if (is.null(ref)) st$w[rows_rr] else chained_reference_weight(ref, spec_t$reference, rows_rr)
              fw <- wr * tm$v[kk]
            }
          }
        } else if (identical(kind, "ipw_legacy")) {
          for (kc in seq_along(spec_t$cohorts)) {
            sub <- which(kk == kc)
            if (!length(sub)) next
            sj <- spec_t$selection_ids[kc]
            sc <- sel_cache[[sj]]
            if (is.null(sc)) {
              r_sub <- rows_rr[sub]
              pp <- if (fuse_det) {
                prob_vec[[sj]][rr[sub]]
              } else {
                cp04_predict_selection(reader, r_sub, backend$selection_specs[[sj]], backend$selections[[sj]])
              }
              coh_k <- spec_t$cohorts[kc]
              sc <- list(
                fw = st$w[r_sub] * st_a[coh_k] / (st$rho[r_sub] * ipw_clip_floor(pp)),
                lam = cbind(
                  1 / st_a[coh_k],
                  -((1 - pp) * ipw_clip_floor_active(pp)) * st_Z[r_sub, , drop = FALSE]
                )
              )
              sel_cache[[sj]] <- sc
            }
            fw[sub] <- sc$fw
            lam_up[sub, (kc - 1L) * (d_Z + 1L) + seq_len(d_Z + 1L)] <- sc$lam
          }
        } else {
          for (kc in seq_along(spec_t$cohorts)) {
            sub <- which(kk == kc)
            if (!length(sub)) next
            qk <- spec_t$selection_ids[kc]
            sc <- sel_cache[[qk]]
            if (is.null(sc)) {
              r_sub <- rows_rr[sub]
              pr_qk <- if (!is.null(sel_op)) sel_op$pair(r_sub, qk) else NULL
              ld <- if (!is.null(sel_op)) {
                if (is.null(ref)) {
                  sel_op$deriv(up_map, qk, r_sub, pred = pr_qk, Z = st_Z, a_vec = st_a)
                } else {
                  sel_op$deriv(up_map, qk, r_sub, pred = pr_qk)
                }
              } else if (is.null(ref)) {
                selection_pair_cohort_log_deriv(sel_fit, up_map, qk, r_sub, pred = pr_qk, Z = st_Z, a_vec = st_a)
              } else {
                selection_pair_cohort_log_deriv(sel_fit, up_map, qk, r_sub, pred = pr_qk)
              }
              sc <- list(
                v = if (fuse_det && identical(kind, "dr_selection")) {
                  or_corollary1_factor(st_a[spec_t$cohorts[kc]], st$rho[r_sub], ld$pred$prob)
                } else NULL,
                fw = ld$pred$weight,
                lam = ld$log_deriv
              )
              sel_cache[[qk]] <- sc
            }
            fw[sub] <- sc$fw
            if (!is.null(fw_score)) fw_score[sub] <- w_ch[rr[sub]] * sc$v
            lam_up[sub, tm$chan_pos[[kc]]] <- sc$lam
          }
        }

        X <- if (m_obs) cpt_design_rows(frame, st, ref, spec_t$reference, rows_rr, trt_cols) else
          matrix(0, 0L, tm$p_trt)
        lv <- if (m_obs) ipw_logit_values(X, f$coefficients) else
          list(prob = numeric(0L), eta = numeric(0L), odds = numeric(0L), active = numeric(0L))
        if (dr) {
          p_or <- tm$p_or
          Z_or <- if (m_obs && p_or > 0L) cpt_design_rows(frame, st, ref, spec_t$reference, rows_rr, tm$kept) else
            matrix(0, m_obs, p_or)
          d_or <- tm$reg_blk$L * Z_or
          rho <- dy - (if (icpt) tm$reg_blk$alpha else 0) -
            (if (p_or > 0L) as.vector(d_or %*% tm$outcome$beta) else 0)
          d_or_arg <- if (p_or > 0L) d_or else NULL
        } else {
          rho <- dy
          d_or_arg <- NULL
        }

        if (!fuse_det) {
          md <- det$models[[tj]]
          xi_up <- md$upstream
          xi_ps <- if (!is.null(md$gamma)) md$gamma else md$beta
          local <- if (separate) det$blocks[[tm$reg_chr]] else NULL
          xi_alpha <- if (icpt) {
            as.double(if (separate) local$alpha else det$alpha[, match(tm$reg_id, det$regression_ids)])
          } else NULL
          xi_beta <- if (dr && tm$p_or > 0L) {
            if (separate) local$beta else det$beta
          } else NULL
        } else {
          if (identical(kind, "ipw_legacy")) {
            xi_up <- matrix(0, n_r, tm$up_n)
            for (kc in seq_along(spec_t$cohorts)) {
              sj <- spec_t$selection_ids[kc]
              ix_u <- (kc - 1L) * (d_Z + 1L) + seq_len(d_Z + 1L)
              xi_up[, ix_u] <- cbind(a_if_chunk[, spec_t$cohorts[kc]], alpha_xi[[sj]])
            }
          } else {
            xi_up <- up_chunk[, tm$up_pos, drop = FALSE]
          }
          score <- matrix(0, n_r, tm$p_trt)
          if (m_obs) {
            if (kind %in% c("ipw_selection", "ipw_legacy")) {
              cp04_finite(c(lv$prob, lv$eta, lv$odds), "treatment replay")
            }
            w_sc <- if (!is.null(fw_score)) fw_score else fw
            score[rr, ] <- (w_sc * (D_obs - lv$prob)) * X
          }
          xi_ps <- t(cp04_solve(f$H_chol, t(score + xi_up %*% tm$t_cross)))
          if (kind %in% c("ipw_selection", "ipw_legacy")) {
            cp04_finite(xi_ps, "complete frame IF")
          }
          if (dr) {
            local <- blk_ifs[[tm$reg_pos]]
            xi_alpha <- if (icpt) as.double(local$alpha) else NULL
            xi_beta <- if (tm$p_or > 0L) local$beta else NULL
          } else {
            xi_alpha <- NULL
            xi_beta <- NULL
          }
        }

        rd_cache[[tj]] <- list(
          rr = rr,
          obs_vecs = list(D_obs, fw, lv$odds, lv$active, rho),
          obs_mats = list(lam_up, X, d_or_arg),
          xi_blocks = list(xi_up, xi_ps, xi_alpha, xi_beta)
        )
      }

      for (i in seq_len(n_cells)) {
        tj <- tjs[i]
        tm <- tj_meta[[tj]]
        rd <- rd_cache[[tj]]
        ok <- fk_fast(
          ptrs[[i]], 1L, rows, rd$rr,
          rd$obs_vecs, rd$obs_mats, rd$xi_blocks,
          tm$prm, gmaps[[i]], cmap_int
        )
        if (identical(ok, FALSE)) {
          cpt_abort(
            paste0("Process extractor returned malformed or nonfinite rows for cell ", ex[[i]]$spec$id, "."),
            "cdid_conditional_pretest_integrity_error"
          )
        }
        if (!isTRUE(ok)) return(NULL)
      }
    }
    if (!is.null(sel_op)) {
      sel_op$finish()
      sel_op <- NULL
    }
  } else {
    # 2b. Control-only OR backends (`or_baseline`, `or_legacy`, `or_selection`)
    separate <- identical(backend$outcome_pooling, "separate")
    icpt <- isTRUE(backend$outcome_recipe$intercept)
    st_Z <- st$Z
    st_a <- st$a
    sel_fit <- if (identical(kind, "or_selection")) backend$selection_backend else NULL
    up_map <- if (identical(kind, "or_selection")) backend$up_map else NULL

    pair_keys <- character(0L)
    pair_objs <- list()
    pair_cohs <- list()
    cell_meta <- vector("list", n_cells)

    add_pair_entry_or <- function(pair_obj, cohs) {
      pk <- cdid_fast_cpt_pair_key(pair_obj)
      p_idx <- match(pk, pair_keys)
      if (is.na(p_idx)) {
        p_idx <- length(pair_keys) + 1L
        pair_keys[p_idx] <<- pk
        pair_objs[[p_idx]] <<- pair_obj
        pair_cohs[[p_idx]] <<- as.integer(cohs)
      } else {
        pair_cohs[[p_idx]] <<- unique(c(pair_cohs[[p_idx]], as.integer(cohs)))
      }
      p_idx
    }

    for (i in seq_len(n_cells)) {
      j <- js[i]
      l <- backend$links[[j]]
      s <- l$support_id
      b <- l$regression
      sp <- backend$support_specs[[s]]
      u <- backend$supports[[s]]
      blk <- backend$regressions[[b]]
      outcome <- if (separate) or_outcome_block_fit(backend$outcome_groups, b) else backend$shared
      p <- outcome$p
      kept <- backend$outcome_recipe$slope_index[outcome$keep]
      channels <- if (identical(kind, "or_baseline")) s else u$channels
      n_ch <- length(channels)
      p_idx <- add_pair_entry_or(sp$pair, sp$cohort)

      cell_meta[[i]] <- list(
        j = j,
        l = l,
        s = s,
        b = b,
        b_chr = if (separate) as.character(b) else NULL,
        sp = sp,
        u = u,
        blk = blk,
        outcome = outcome,
        p = p,
        kept = kept,
        channels = channels,
        n_ch = n_ch,
        pair_idx = p_idx,
        rate_log_d = if (identical(kind, "or_baseline")) {
          if (u$rate > ipw_clip_eps) -1 / u$rate else 0
        } else 0,
        prm = as.double(c(l$d, n_ch, if (icpt) 1 else 0))
      )
    }

    sel_ok <- !is.null(sel_fit) && cdid_are_default_fns(c("selection_prediction_operation", "predict_pair_selection"))
    fuse_det <- separate && is.null(ref) && switch(
      kind,
      or_selection = sel_ok && cdid_are_default_fns(c(
        "cpt_backend_details", "evaluate_or_corollary1_if_chunk", "evaluate_or_selection_if_chunk",
        "or_corollary1_if_budget", "or_selection_separate_if_core", "validate_or_corollary1_backend"
      )),
      or_legacy = cdid_are_default_fns(c(
        "cpt_backend_details", "evaluate_or_corollary1_if_chunk", "or_corollary1_if_budget",
        "or_selection_separate_if_core", "validate_or_corollary1_backend", "cp04_predict_selection"
      )),
      or_baseline = cdid_are_default_fns(c(
        "cpt_backend_details", "evaluate_or_baseline_if_chunk", "or_baseline_if_budget",
        "or_baseline_separate_if_core", "validate_or_baseline_backend"
      )),
      FALSE
    )

    gids_need <- sids_need <- uids_need <- NULL
    sup_pair_idx <- NULL
    blk_meta <- NULL

    if (fuse_det) {
      if (identical(kind, "or_baseline")) {
        validate_or_baseline_backend(backend)
      } else {
        validate_or_corollary1_backend(backend)
      }
      gids_need <- sort(unique(vapply(backend$links[u_js], `[[`, integer(1L), "regression")))
      sids_need <- sort(unique(c(
        vapply(backend$links[u_js], `[[`, integer(1L), "support_id"),
        unlist(lapply(backend$regression_specs[gids_need], `[[`, "support_ids"))
      )))
      uids_need <- if (identical(kind, "or_baseline")) {
        sids_need
      } else {
        sort(unique(unlist(lapply(backend$supports[sids_need], `[[`, "channels"))))
      }
      sup_pair_idx <- integer(length(backend$supports))
      for (sid in sids_need) {
        sp_s <- backend$support_specs[[sid]]
        sup_pair_idx[sid] <- add_pair_entry_or(sp_s$pair, sp_s$cohort)
      }
      blk_meta <- vector("list", length(gids_need))
      for (kb in seq_along(gids_need)) {
        b <- gids_need[kb]
        of <- or_outcome_block_fit(backend$outcome_groups, b)
        rec <- backend$regressions[[b]]
        spos <- match(backend$regression_specs[[b]]$support_ids, sids_need)
        blk_meta[[kb]] <- list(
          of = of,
          rec = rec,
          spos = spos,
          upos = if (identical(kind, "or_baseline")) spos else match(rec$channel_ids, uids_need),
          kept_or = backend$outcome_recipe$slope_index[of$keep],
          t_forcing_beta = if (of$p > 0L) t(rec$forcing_beta) else NULL
        )
      }
      for (i in seq_len(n_cells)) {
        cm <- cell_meta[[i]]
        cell_meta[[i]]$s_pos <- match(cm$s, sids_need)
        cell_meta[[i]]$b_pos <- match(cm$b, gids_need)
        cell_meta[[i]]$ch_upos <- match(cm$channels, uids_need)
      }
      if (!is.null(sel_fit) && identical(sel_fit$selection_model, "pair_logit")) {
        act_blks <- which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids_need), logical(1L)))
        for (m_b in act_blks) {
          for (cid in sel_fit$blocks[[m_b]]$cell_ids) {
            cl_s <- sel_fit$cells[[cid]]
            add_pair_entry_or(cl_s$pair, cl_s$cohort)
          }
        }
      }
    }

    n_pairs <- length(pair_keys)
    n_sup_max <- length(backend$supports)
    pair_env <- new.env(parent = emptyenv())
    pair_env$offset <- 0L
    pair_env$pair_S <- NULL
    chunk_pair_reader <- if (fuse_det && !is.null(sel_fit) && identical(sel_fit$selection_model, "pair_logit")) {
      function(r_sub, pair_obj, read_outcomes = FALSE) {
        if (!isTRUE(read_outcomes) && !is.null(pair_env$pair_S)) {
          p_idx <- match(cdid_fast_cpt_pair_key(pair_obj), pair_keys)
          if (!is.na(p_idx)) {
            return(list(S = as.double(pair_env$pair_S[[p_idx]][r_sub - pair_env$offset])))
          }
        }
        cp04_pair_rows(reader, r_sub, pair_obj, isTRUE(read_outcomes))
      }
    } else NULL

    sel_op <- if (sel_ok) {
      selection_prediction_operation(
        sel_fit,
        source_reader = reader,
        source_reference = ref,
        pair_response_reader = chunk_pair_reader
      )
    } else NULL
    on.exit(if (!is.null(sel_op)) sel_op$abort(), add = TRUE)

    for (k_ch in seq_along(starts)) {
      rows <- seq.int(starts[k_ch], min(st$N, starts[k_ch] + B - 1L))
      n_r <- length(rows)
      ci <- st$cohort_index[rows]

      if (fuse_det) {
        if (identical(kind, "or_baseline")) {
          or_baseline_if_budget(backend, n_r, u_js, TRUE)
        } else {
          or_corollary1_if_budget(backend, n_r, u_js, TRUE)
        }
      }

      pair_S <- vector("list", n_pairs)
      pair_dy <- vector("list", n_pairs)
      for (p_idx in seq_len(n_pairs)) {
        cand <- which(ci %in% pair_cohs[[p_idx]])
        if (!identical(kind, "or_baseline")) cand <- cand[st$U[rows[cand]]]
        S_vec <- logical(n_r)
        dy_vec <- numeric(n_r)
        if (length(cand)) {
          hit_p <- cp04_pair_rows(reader, rows[cand], pair_objs[[p_idx]], TRUE)
          ob_p <- which(hit_p$S == 1)
          if (length(ob_p)) {
            obs_c <- cand[ob_p]
            S_vec[obs_c] <- TRUE
            dy_vec[obs_c] <- as.double(hit_p$dy[ob_p])
          }
        }
        pair_S[[p_idx]] <- S_vec
        pair_dy[[p_idx]] <- dy_vec
      }
      pair_env$offset <- starts[k_ch] - 1L
      pair_env$pair_S <- pair_S

      det <- NULL
      up_chunk <- NULL
      sup_hits <- NULL
      blk_ifs <- NULL
      w_ch <- if (fuse_det) st$w[rows] else NULL

      if (!fuse_det) {
        det <- cpt_backend_details(backend, kind, rows, u_js)
      } else {
        n_s <- length(sids_need)
        sup_hits <- vector("list", n_s)
        if (identical(kind, "or_selection")) {
          up_chunk <- sel_op$upstream(up_map, rows, uids_need, FALSE)$up
          for (j_s in seq_len(n_s)) {
            sid <- sids_need[j_s]
            sp_s <- backend$support_specs[[sid]]
            p_idx <- sup_pair_idx[sid]
            obs_s <- which(pair_S[[p_idx]] & (ci == sp_s$cohort))
            r_obs <- rows[obs_s]
            if (length(obs_s)) {
              pr_qk <- sel_op$pair(r_obs, backend$q_of_s[sid])
              v_obs <- or_corollary1_factor(st_a[sp_s$cohort], st$rho[r_obs], pr_qk$prob)
            } else {
              pr_qk <- NULL
              v_obs <- numeric(0L)
            }
            sup_hits[[j_s]] <- list(
              obs = obs_s,
              r_obs = r_obs,
              dy = pair_dy[[p_idx]][obs_s],
              v = v_obs,
              pr = pr_qk
            )
          }
        } else if (identical(kind, "or_legacy")) {
          U_ch <- st$U[rows]
          A_ch <- as.double(st$a_index[rows] > 0L)
          up_chunk <- matrix(0, n_r, length(uids_need))
          for (ju in which(uids_need <= backend$upstream$n_a)) {
            k_c <- backend$upstream$cohort[uids_need[ju]]
            ix_c <- which(ci == k_c & U_ch)
            if (length(ix_c)) up_chunk[ix_c, ju] <- w_ch[ix_c] * (A_ch[ix_c] - st_a[k_c]) / st$mass_U[k_c]
          }
          for (j_s in seq_len(n_s)) {
            sid <- sids_need[j_s]
            sp_s <- backend$support_specs[[sid]]
            f_s <- backend$selections[[sid]]
            p_idx <- sup_pair_idx[sid]
            ix_s <- which(ci == sp_s$cohort & U_ch)
            if (length(ix_s)) {
              pr_ix <- cp04_predict_selection(reader, rows[ix_s], backend$selection_specs[[sid]], f_s)
              S_ix <- as.double(pair_S[[p_idx]][ix_s])
              up_chunk[ix_s, match(backend$supports[[sid]]$channels[-1L], uids_need)] <-
                t(cp04_solve(f_s$H_chol, t(w_ch[ix_s] * (S_ix - pr_ix) * st_Z[rows[ix_s], , drop = FALSE])))
              ob_l <- which(S_ix == 1)
              obs_s <- ix_s[ob_l]
              r_obs <- rows[obs_s]
              pr_obs <- pr_ix[ob_l]
              v_obs <- if (length(obs_s)) or_corollary1_factor(st_a[sp_s$cohort], st$rho[r_obs], pr_obs) else numeric(0L)
            } else {
              obs_s <- integer(0L)
              r_obs <- integer(0L)
              pr_obs <- numeric(0L)
              v_obs <- numeric(0L)
            }
            sup_hits[[j_s]] <- list(
              obs = obs_s,
              r_obs = r_obs,
              dy = pair_dy[[p_idx]][obs_s],
              v = v_obs,
              pr_obs = pr_obs
            )
          }
        } else {
          A_ch <- as.double(st$a_index[rows] > 0L)
          up_chunk <- matrix(0, n_r, n_s)
          for (j_s in seq_len(n_s)) {
            sid <- sids_need[j_s]
            sp_s <- backend$support_specs[[sid]]
            u_s <- backend$supports[[sid]]
            p_idx <- sup_pair_idx[sid]
            ix_s <- which(ci == sp_s$cohort)
            if (length(ix_s)) {
              S_ix <- as.double(pair_S[[p_idx]][ix_s])
              up_chunk[ix_s, j_s] <- w_ch[ix_s] * (S_ix - A_ch[ix_s] * u_s$rate) / u_s$massA
              obs_s <- ix_s[S_ix == 1]
            } else {
              obs_s <- integer(0L)
            }
            sup_hits[[j_s]] <- list(
              obs = obs_s,
              r_obs = rows[obs_s],
              dy = pair_dy[[p_idx]][obs_s],
              v = rep.int(u_s$v, length(obs_s))
            )
          }
        }

        blk_ifs <- vector("list", length(gids_need))
        N_grp <- backend$outcome_groups$N
        for (kb in seq_along(gids_need)) {
          bm <- blk_meta[[kb]]
          f_or <- bm$of
          rec <- bm$rec
          ha <- numeric(n_r)
          hb <- matrix(0, n_r, f_or$p)
          for (s_pos in bm$spos) {
            h <- sup_hits[[s_pos]]
            ix <- h$obs
            if (length(ix)) {
              d_b <- if (f_or$p > 0L) {
                f_or$L * frame$Z_a[st$a_index[h$r_obs], bm$kept_or, drop = FALSE]
              } else {
                matrix(0, length(ix), 0L)
              }
              q_b <- w_ch[ix] * h$v
              r_b <- h$dy - f_or$alpha - (if (f_or$p > 0L) as.vector(d_b %*% f_or$beta) else 0)
              qr_b <- q_b * r_b
              ha[ix] <- ha[ix] + qr_b
              if (f_or$p > 0L) {
                hb[ix, ] <- hb[ix, , drop = FALSE] + qr_b * (if (icpt) sweep(d_b, 2L, f_or$mean_d) else d_b)
              }
            }
          }
          us_b <- up_chunk[, bm$upos, drop = FALSE]
          xb_b <- hb
          if (f_or$p > 0L) {
            rhs <- sweep(hb, 2L, rec$score_mean_beta) + us_b %*% bm$t_forcing_beta
            xb_b <- t(if (icpt) {
              or_pooled_wls_gram_solve(f_or$fit, t(rhs), normalizer = N_grp)
            } else {
              or_wls_gram_solve(f_or$fit, t(rhs), normalizer = N_grp)
            })
          }
          xa_b <- if (icpt) {
            (ha - rec$score_mean_alpha + as.vector(us_b %*% rec$forcing_alpha)) / rec$mass -
              (if (f_or$p > 0L) as.vector(xb_b %*% rec$mean_d) else 0)
          } else NULL
          blk_ifs[[kb]] <- list(alpha = xa_b, beta = xb_b)
        }
      }

      sup_cache <- vector("list", n_sup_max)

      for (i in seq_len(n_cells)) {
        cm <- cell_meta[[i]]
        s <- cm$s
        sp <- cm$sp
        u <- cm$u
        sc <- sup_cache[[s]]
        if (is.null(sc)) {
          if (fuse_det) {
            h_s <- sup_hits[[cm$s_pos]]
            rr <- h_s$obs
            rows_rr <- h_s$r_obs
            dy <- h_s$dy
            m_obs <- length(rr)
            if (m_obs) {
              q_obs <- w_ch[rr] * h_s$v
              if (identical(kind, "or_baseline")) {
                lam_obs <- matrix(cm$rate_log_d, m_obs, 1L)
              } else if (identical(kind, "or_legacy")) {
                lam_obs <- cbind(
                  1 / st_a[sp$cohort],
                  -((1 - h_s$pr_obs) * ipw_clip_floor_active(h_s$pr_obs)) * st_Z[rows_rr, , drop = FALSE]
                )
              } else {
                ld <- sel_op$deriv(up_map, backend$q_of_s[s], rows_rr, pred = h_s$pr, Z = st_Z, a_vec = st_a)
                lam_obs <- ld$log_deriv
              }
            } else {
              q_obs <- numeric(0L)
              lam_obs <- matrix(0, 0L, cm$n_ch)
            }
          } else {
            S_vec <- pair_S[[cm$pair_idx]]
            rr <- which(S_vec & (ci == sp$cohort))
            rows_rr <- rows[rr]
            dy <- pair_dy[[cm$pair_idx]][rr]
            m_obs <- length(rr)
            if (m_obs) {
              if (identical(kind, "or_baseline")) {
                w_obs <- if (is.null(ref)) st$w[rows_rr] else chained_reference_weight(ref, sp$reference, rows_rr)
                q_obs <- w_obs * u$v
                lam_obs <- matrix(cm$rate_log_d, m_obs, 1L)
              } else if (identical(kind, "or_legacy")) {
                pr <- cp04_predict_selection(reader, rows_rr, backend$selection_specs[[s]], backend$selections[[s]])
                q_obs <- st$w[rows_rr] * or_corollary1_factor(st_a[sp$cohort], st$rho[rows_rr], pr)
                lam_obs <- cbind(
                  1 / st_a[sp$cohort],
                  -((1 - pr) * ipw_clip_floor_active(pr)) * st_Z[rows_rr, , drop = FALSE]
                )
              } else {
                qk <- backend$q_of_s[s]
                pr_qk <- if (!is.null(sel_op)) sel_op$pair(rows_rr, qk) else NULL
                ld <- if (!is.null(sel_op)) {
                  if (is.null(ref)) {
                    sel_op$deriv(up_map, qk, rows_rr, pred = pr_qk, Z = st_Z, a_vec = st_a)
                  } else {
                    sel_op$deriv(up_map, qk, rows_rr, pred = pr_qk)
                  }
                } else if (is.null(ref)) {
                  selection_pair_cohort_log_deriv(sel_fit, up_map, qk, rows_rr, pred = pr_qk, Z = st_Z, a_vec = st_a)
                } else {
                  selection_pair_cohort_log_deriv(sel_fit, up_map, qk, rows_rr, pred = pr_qk)
                }
                a_val <- if (is.null(ref)) st_a[sp$cohort] else chained_reference_stat(ref, sp$reference, sp$cohort, "a")
                w_obs <- if (is.null(ref)) st$w[rows_rr] else chained_reference_weight(ref, sp$reference, rows_rr)
                q_obs <- w_obs * or_corollary1_factor(a_val, st$rho[rows_rr], ld$pred$prob)
                lam_obs <- ld$log_deriv
              }
            } else {
              q_obs <- numeric(0L)
              lam_obs <- matrix(0, 0L, cm$n_ch)
            }
          }
          sc <- list(rr = rr, rows_rr = rows_rr, dy = dy, m_obs = m_obs, q = q_obs, lam = lam_obs)
          sup_cache[[s]] <- sc
        }

        m_obs <- sc$m_obs
        p <- cm$p
        if (m_obs) {
          z <- if (is.null(ref)) frame$Z_a[st$a_index[sc$rows_rr], cm$kept, drop = FALSE] else
            chained_reference_design(ref, sp$reference, sc$rows_rr)[, cm$kept, drop = FALSE]
          d_mat <- cm$blk$L * z
          rho <- sc$dy - cm$blk$alpha - (if (p > 0L) as.vector(d_mat %*% cm$outcome$beta) else 0)
        } else {
          d_mat <- matrix(0, 0L, p)
          rho <- numeric(0L)
        }
        d_or_arg <- if (p > 0L) d_mat else NULL

        if (!fuse_det) {
          xi_ch <- if (identical(kind, "or_baseline")) {
            det$rates[, match(s, det$support_ids), drop = FALSE]
          } else {
            det$upstream[, match(cm$channels, det$upstream_ids), drop = FALSE]
          }
          local <- if (separate) det$blocks[[cm$b_chr]] else NULL
          xi_alpha <- if (icpt) {
            as.double(if (separate) local$alpha else det$alpha[, match(cm$b, det$regression_ids)])
          } else NULL
          xi_beta <- if (p > 0L) {
            if (separate) local$beta else det$beta
          } else NULL
        } else {
          xi_ch <- up_chunk[, cm$ch_upos, drop = FALSE]
          local <- blk_ifs[[cm$b_pos]]
          xi_alpha <- if (icpt) as.double(local$alpha) else NULL
          xi_beta <- if (p > 0L) local$beta else NULL
        }

        ok <- fk_fast(
          ptrs[[i]], 2L, rows, sc$rr,
          list(sc$q, rho),
          list(sc$lam, d_or_arg),
          list(xi_ch, xi_alpha, xi_beta),
          cm$prm, gmaps[[i]], cmap_int
        )
        if (identical(ok, FALSE)) {
          cpt_abort(
            paste0("Process extractor returned malformed or nonfinite rows for cell ", ex[[i]]$spec$id, "."),
            "cdid_conditional_pretest_integrity_error"
          )
        }
        if (!isTRUE(ok)) return(NULL)
      }
    }
    if (!is.null(sel_op)) {
      sel_op$finish()
      sel_op <- NULL
    }
  }

  built <- vector("list", n_cells)
  for (i in seq_len(n_cells)) {
    grid <- cache[[keys[i]]]$distinct
    cl_out <- cdid_fast_cpt_finalize_ptr(fk_fin, ptrs[[i]], ex[[i]]$spec, grid, cmap_int, N)
    if (is.null(cl_out)) return(NULL)
    cl_out$estimate_scalar <- cells[[i]]$sign * ex[[i]]$estimate
    cl_out$route <- ex[[i]]$route
    cl_out$recipe_key <- keys[i]
    built[[i]] <- cl_out
  }
  built
}
