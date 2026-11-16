# Fast C-backed wrappers for Imputation, Augmented, Reconstruction, and Multiplier paths.

cdid_ia_fast_imp_if_chunk_core <- function(st, X_F, backend, rows, links, scale, pre) {
  fk <- cdid_fast_kernel("C_cdid_ia_imp_if_chunk")
  if (is.null(fk)) return(NULL)
  rows_i <- as.integer(rows)
  links_i <- as.integer(links)
  c_env <- cdid_ia_get_env_cache(backend, augmented = FALSE)
  if (!is.null(c_env) && st$N <= 10000L && is.null(rownames(X_F))) {
    phi_frame <- c_env$phi_all_frame
    if (is.null(phi_frame) && !isTRUE(c_env$phi_failed)) {
      reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
      panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
      all_links <- seq_along(backend$links)
      all_pre <- if (identical(links_i, all_links)) {
        pre
      } else {
        imputation_if_needed_sets(backend, all_links, backend$est_method, FALSE)
      }
      phi_frame <- fk(backend, st, panel, X_F, seq_len(st$N), all_links, all_pre, "frame")
      if (is.null(phi_frame)) {
        c_env$phi_failed <- TRUE
      } else {
        c_env$phi_all_frame <- phi_frame
        c_env$phi_all_core <- if (st$lambda == 1.0) phi_frame else st$lambda * phi_frame
        cdid_fast_kernel_cache$last_retained_if <- list(
          backend_ref = backend$ref,
          links = all_links,
          phi_core = phi_frame,
          lambda = 1.0,
          N = as.integer(st$N)
        )
      }
    }
    if (!is.null(phi_frame)) {
      src_mat <- if (identical(scale, "core")) c_env$phi_all_core else phi_frame
      n_all_links <- ncol(src_mat)
      phi_out <- if (length(rows_i) == st$N && rows_i[1L] == 1L && rows_i[st$N] == st$N &&
                     identical(rows_i, seq_len(st$N)) &&
                     length(links_i) == n_all_links && identical(links_i, seq_len(n_all_links))) {
        src_mat
      } else {
        src_mat[rows_i, links_i, drop = FALSE]
      }
      return(list(
        rows = rows_i,
        id = st$id[rows_i],
        links = links_i,
        phi = phi_out,
        scale = scale,
        normalization_count = if (identical(scale, "core")) st$n_A else st$N,
        details = NULL
      ))
    }
  }
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  phi_out <- fk(backend, st, panel, X_F, rows_i, links_i, pre, as.character(scale))
  if (is.null(phi_out)) return(NULL)
  list(
    rows = rows_i,
    id = st$id[rows_i],
    links = links_i,
    phi = phi_out,
    scale = scale,
    normalization_count = if (identical(scale, "core")) st$n_A else st$N,
    details = NULL
  )
}

cdid_ia_fast_aug_if_chunk_core <- function(st, X_F, backend, rows, links, scale, pre, pair_op = NULL) {
  fk <- cdid_fast_kernel("C_cdid_ia_aug_if_chunk")
  if (is.null(fk)) return(NULL)
  rows_i <- as.integer(rows)
  links_i <- as.integer(links)
  c_env <- cdid_ia_get_env_cache(backend, augmented = TRUE)
  if (!is.null(c_env) && st$N <= 10000L && is.null(rownames(X_F))) {
    phi_frame <- c_env$phi_all_frame
    if (is.null(phi_frame) && !isTRUE(c_env$phi_failed)) {
      reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
      panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
      all_links <- seq_along(backend$links)
      all_pre <- if (identical(links_i, all_links)) {
        pre
      } else {
        augmented_if_needed_sets(backend, all_links, FALSE)
      }
      phi_frame <- fk(backend, st, panel, X_F, seq_len(st$N), all_links, all_pre, "frame", NULL)
      if (is.null(phi_frame)) {
        c_env$phi_failed <- TRUE
      } else {
        c_env$phi_all_frame <- phi_frame
        c_env$phi_all_core <- if (st$lambda == 1.0) phi_frame else st$lambda * phi_frame
        cdid_fast_kernel_cache$last_retained_if <- list(
          backend_ref = backend$ref,
          links = all_links,
          phi_core = phi_frame,
          lambda = 1.0,
          N = as.integer(st$N)
        )
      }
    }
    if (!is.null(phi_frame)) {
      src_mat <- if (identical(scale, "core")) c_env$phi_all_core else phi_frame
      n_all_links <- ncol(src_mat)
      phi_out <- if (length(rows_i) == st$N && rows_i[1L] == 1L && rows_i[st$N] == st$N &&
                     identical(rows_i, seq_len(st$N)) &&
                     length(links_i) == n_all_links && identical(links_i, seq_len(n_all_links))) {
        src_mat
      } else {
        src_mat[rows_i, links_i, drop = FALSE]
      }
      return(list(
        rows = rows_i,
        id = st$id[rows_i],
        links = links_i,
        phi = phi_out,
        scale = scale,
        normalization_count = if (identical(scale, "core")) st$n_A else st$N,
        details = NULL
      ))
    }
  }
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  phi_out <- fk(backend, st, panel, X_F, rows_i, links_i, pre, as.character(scale), NULL)
  if (is.null(phi_out)) return(NULL)
  list(
    rows = rows_i,
    id = st$id[rows_i],
    links = links_i,
    phi = phi_out,
    scale = scale,
    normalization_count = if (identical(scale, "core")) st$n_A else st$N,
    details = NULL
  )
}

cdid_ia_fast_cell_replay_build <- function(backend, st, reader, X, ids, cids, zids, B_build) {
  fk <- cdid_fast_kernel("C_cdid_ia_cell_replay_build")
  if (is.null(fk)) return(NULL)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  fk(backend, st, panel, X, as.integer(ids), as.integer(cids), as.integer(zids), as.integer(B_build))
}

cdid_ia_fast_imp_stage1 <- function(st, panel, X_F, cat, recipe, or_frame_index, B_chunk) {
  if (is.null(panel)) return(NULL)
  fk <- cdid_fast_kernel("C_cdid_ia_imp_stage1")
  if (is.null(fk)) return(NULL)
  fk(st, panel, X_F, cat, recipe, as.integer(or_frame_index), as.integer(B_chunk))
}

cdid_ia_fast_imp_stage2_cells <- function(st, X_F, cat, recipe, or_frame_index, B_chunk, cell_meta, used_ctrl) {
  fk <- cdid_fast_kernel("C_cdid_ia_imp_stage2_cells")
  if (is.null(fk)) return(NULL)
  fk(st, X_F, cat, recipe, as.integer(or_frame_index), as.integer(B_chunk), cell_meta,
    if (is.null(used_ctrl)) NULL else as.integer(used_ctrl))
}

cdid_ia_fast_imp_stage24_forcing_links <- function(st, X_F, cat, recipe, trt_recipe, B_chunk,
                                                   cells_fitted, group_fits, block_specs = NULL,
                                                   trt_fits = NULL, do_block_forcing = TRUE) {
  fk <- cdid_fast_kernel("C_cdid_ia_imp_stage24_forcing_links")
  if (is.null(fk)) return(NULL)
  fk(st, X_F, cat, recipe, trt_recipe, as.integer(B_chunk), cells_fitted, group_fits,
    list(block_specs = block_specs, trt_fits = trt_fits, do_block_forcing = isTRUE(do_block_forcing)))
}

cdid_ia_fast_aug_stage456 <- function(mode, st, panel, X_F, sel_fit, up_map, cat, recipe,
                                      or_frame_index, trt_recipe, B_chunk, cells, groups,
                                      q_of_cell, theta_ids, d_q,
                                      cell_ev_cache = NULL, block_specs = NULL,
                                      trt_fits = NULL, separate_or = FALSE) {
  fk <- cdid_fast_kernel("C_cdid_ia_aug_stage456")
  if (is.null(fk)) return(NULL)
  mode_i <- as.integer(mode)
  if (mode_i == 1L) {
    if (is.null(panel)) return(NULL)
    used_ctrl <- if (length(cat$or_blocks) > 0L) {
      as.integer(sort(unique(unlist(lapply(cat$or_blocks, `[[`, "control_cell_ids")))))
    } else {
      NULL
    }
    res <- fk(
      1L, st, panel, X_F, cat, recipe, as.integer(or_frame_index), trt_recipe,
      as.integer(B_chunk), cells,
      list(
        groups = groups,
        sel_fit = sel_fit,
        q_of_cell = as.integer(q_of_cell),
        theta_ids = theta_ids,
        used_ctrl = used_ctrl
      )
    )
    if (is.null(res)) return(NULL)
    n_cells <- length(cells)
    cells_out <- vector("list", n_cells)
    for (cid in seq_len(n_cells)) {
      cf <- cells[[cid]]
      cf$q_id <- q_of_cell[cid]
      cf$theta_ids <- theta_ids[[cid]]
      cells_out[[cid]] <- c(cf, res$cell_aug[[cid]])
    }
    return(list(
      cells = cells_out,
      ctrl_acc = res$ctrl_acc,
      cohort_colmax_F = res$cohort_colmax_F,
      cell_ev_cache = res$cell_ev_cache
    ))
  }
  if (mode_i == 2L) {
    res <- fk(
      2L, st, panel, X_F, cat, recipe, as.integer(or_frame_index), trt_recipe,
      as.integer(B_chunk), cells,
      list(
        groups = groups,
        cell_ev_cache = cell_ev_cache,
        theta_ids = theta_ids,
        block_specs = block_specs,
        treatments = trt_fits
      )
    )
    if (is.null(res)) return(NULL)
    N <- st$N
    bf_out <- res$block_forcing
    q_beta_global <- NULL
    if (!is.null(bf_out) && !is.null(block_specs)) {
      n_blocks <- length(block_specs)
      if (!separate_or && n_blocks > 0L) {
        bs1 <- block_specs[[1L]]
        p_or_g <- length(bs1$beta_or)
        q_beta_ids_g <- bs1$q_beta_ids
        q_beta_global <- matrix(0.0, p_or_g, length(q_beta_ids_g))
      }
      for (bid in seq_len(n_blocks)) {
        b_spec <- cat$or_blocks[[bid]]
        bs_f <- block_specs[[bid]]
        cids_b <- b_spec$control_cell_ids
        p_or <- length(bs_f$beta_or)
        q_alpha_ids <- bs_f$q_alpha_ids
        q_beta_ids <- bs_f$q_beta_ids
        q_alpha <- numeric(length(q_alpha_ids))
        if (separate_or) q_beta_loc <- matrix(0.0, p_or, length(q_beta_ids))
        rq_list <- bf_out[[bid]]$Rq_dtilde_list
        for (idx_k in seq_along(cids_b)) {
          cid <- cids_b[idx_k]
          if (d_q[cid] > 0L) {
            cf <- cells[[cid]]
            pos <- match(cf$theta_ids, q_alpha_ids)
            q_alpha[pos] <- q_alpha[pos] - cf$Rq_one
            if (p_or > 0L) {
              gpos <- match(cf$theta_ids, q_beta_ids)
              if (separate_or) {
                q_beta_loc[, gpos] <- q_beta_loc[, gpos, drop = FALSE] - rq_list[[idx_k]] / N
              } else {
                q_beta_global[, gpos] <- q_beta_global[, gpos, drop = FALSE] - rq_list[[idx_k]] / N
              }
            }
          }
        }
        bf_out[[bid]]$q_alpha <- q_alpha
        if (separate_or) bf_out[[bid]]$q_beta <- q_beta_loc
      }
    }
    lc_out <- res$link_ctrl
    if (!is.null(lc_out)) {
      for (ell in seq_along(lc_out)) {
        lk <- cat$links[[ell]]
        cids_l <- lk$control_cell_ids
        sel_C_ids <- as.integer(sort(unique(c(integer(0L), unlist(theta_ids[cids_l])))))
        sel_C <- numeric(length(sel_C_ids))
        rq_o_k <- lc_out[[ell]]$Rq_o_k
        for (idx_k in seq_along(cids_l)) {
          cid_k <- cids_l[idx_k]
          if (d_q[cid_k] > 0L) {
            pos <- match(theta_ids[[cid_k]], sel_C_ids)
            sel_C[pos] <- sel_C[pos] - rq_o_k[[idx_k]]
          }
        }
        lc_out[[ell]]$sel_C_ids <- sel_C_ids
        lc_out[[ell]]$sel_C <- sel_C
      }
    }
    return(list(
      block_forcing = bf_out,
      q_beta_global = q_beta_global,
      link_ctrl = lc_out
    ))
  }
  NULL
}

cdid_ia_fast_check_design <- function(X, Z_a, a_index, x_match) {
  fk <- cdid_fast_kernel("C_cdid_ia_check_design_matrix")
  if (is.null(fk)) return(FALSE)
  res <- fk(X, Z_a, as.integer(a_index), as.integer(x_match))
  identical(res, 0L)
}

cdid_ia_sexp_aliases <- function(a, b) {
  fk <- cdid_fast_kernel("C_gmm_sexp_aliases")
  if (!is.null(fk)) isTRUE(fk(a, b)) else gmm_sexp_aliases(a, b)
}

cdid_ia_get_env_cache <- function(backend, augmented = FALSE) {
  ext <- get0(if (augmented) "val_aug_extra" else "val_imp_extra",
              envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (!is.list(ext) || !is.environment(ext$cache) ||
      !identical(ext$backend_ref, .subset2(backend, "ref")) ||
      !cdid_ia_sexp_aliases(ext$backend, backend)) {
    return(NULL)
  }
  if (augmented) {
    if (!cdid_ia_check_validated_aug_backend(backend)) return(NULL)
  } else {
    if (!cdid_ia_check_validated_imp_backend(backend)) return(NULL)
  }
  ext$cache
}

cdid_ia_record_validated_imp_backend <- function(backend, reader, frame, panel, d_owner) {
  cdid_record_validated_backend(backend, "imputation_backend", reader, frame, panel, reader$frame_ref)
  assign("val_imp_extra", list(
    backend = backend,
    backend_ref = .subset2(backend, "ref"),
    design_ref = .subset2(backend, "design_ref"),
    d_owner = d_owner,
    cache = new.env(parent = emptyenv())
  ), envir = cdid_fast_kernel_cache)
  invisible(NULL)
}

cdid_ia_check_validated_imp_backend <- function(backend) {
  if (!is.list(backend) || !inherits(backend, "cdid_imputation_backend") ||
      !is.null(.subset2(backend, "reference_ref")) ||
      !cdid_are_default_fns(c("validate_imputation_design_owner", "chained_backend_reference"))) {
    return(FALSE)
  }
  if (!cdid_check_validated_backend(backend, "imputation_backend")) return(FALSE)
  ext <- get0("val_imp_extra", envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (!is.list(ext) || !cdid_ia_sexp_aliases(ext$backend, backend) ||
      !identical(ext$backend_ref, .subset2(backend, "ref"))) {
    return(FALSE)
  }
  dref <- .subset2(backend, "design_ref")
  if (!is.null(dref)) {
    if (!is.environment(dref) || !identical(dref, ext$design_ref) ||
        !identical(parent.env(dref), emptyenv()) || !environmentIsLocked(dref) ||
        !exists("imputation_design_owner", envir = dref, inherits = FALSE) ||
        !bindingIsLocked("imputation_design_owner", dref) ||
        bindingIsActive("imputation_design_owner", dref) ||
        !cdid_ia_sexp_aliases(get("imputation_design_owner", envir = dref, inherits = FALSE), ext$d_owner)) {
      return(FALSE)
    }
  }
  TRUE
}

cdid_ia_record_validated_aug_backend <- function(backend, reader, frame, panel, sel_b, d_owner) {
  cdid_record_validated_backend(backend, "augmented_backend", reader, frame, panel, reader$frame_ref)
  assign("val_aug_extra", list(
    backend = backend,
    backend_ref = .subset2(backend, "ref"),
    sel_b = sel_b,
    sel_ref = .subset2(sel_b, "ref"),
    design_ref = .subset2(backend, "design_ref"),
    d_owner = d_owner,
    cache = new.env(parent = emptyenv())
  ), envir = cdid_fast_kernel_cache)
  invisible(NULL)
}

cdid_ia_check_validated_aug_backend <- function(backend) {
  if (!is.list(backend) || !inherits(backend, "cdid_augmented_backend") ||
      !is.null(.subset2(backend, "layout")) || !is.null(.subset2(backend, "reference_ref")) ||
      !cdid_are_default_fns(c(
        "validate_selection_models", "validate_imputation_design_owner", "chained_backend_reference"
      ))) {
    return(FALSE)
  }
  if (!cdid_check_validated_backend(backend, "augmented_backend")) return(FALSE)
  ext <- get0("val_aug_extra", envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (!is.list(ext) || !cdid_ia_sexp_aliases(ext$backend, backend) ||
      !identical(ext$backend_ref, .subset2(backend, "ref"))) {
    return(FALSE)
  }
  sref <- ext$sel_ref
  if (!is.environment(sref) || !identical(parent.env(sref), emptyenv()) ||
      !environmentIsLocked(sref) ||
      !exists("selection_models", envir = sref, inherits = FALSE) ||
      !bindingIsLocked("selection_models", sref) ||
      bindingIsActive("selection_models", sref) ||
      !cdid_ia_sexp_aliases(get("selection_models", envir = sref, inherits = FALSE), ext$sel_b)) {
    return(FALSE)
  }
  dref <- .subset2(backend, "design_ref")
  if (!is.environment(dref) || !identical(dref, ext$design_ref) ||
      !identical(parent.env(dref), emptyenv()) || !environmentIsLocked(dref) ||
      !exists("imputation_design_owner", envir = dref, inherits = FALSE) ||
      !bindingIsLocked("imputation_design_owner", dref) ||
      bindingIsActive("imputation_design_owner", dref) ||
      !cdid_ia_sexp_aliases(get("imputation_design_owner", envir = dref, inherits = FALSE), ext$d_owner)) {
    return(FALSE)
  }
  TRUE
}

cdid_ia_solve_adjoint_fast <- function(groups, cells, n_cells, N, intercept,
                                       zero_tilde_nu, structural_cids, lam_vals, nu_vals) {
  tilde_lambda <- numeric(n_cells)
  tilde_nu <- zero_tilde_nu
  n_sc <- length(structural_cids)
  if (n_sc == 1L) {
    cid <- structural_cids[1L]
    cf <- cells[[cid]]
    gid <- cf$group_id
    gf <- groups[[gid]]
    lam_c <- if (intercept) lam_vals[1L] else 0.0
    if (intercept) tilde_lambda[cid] <- lam_c / cf$mass
    if (gf$p > 0L) {
      rhs_r <- nu_vals[[1L]] - lam_c * cf$mean_d
      sfit <- if (intercept) gf$fit$slope_fit else gf$fit
      R_g <- sfit$R
      piv_g <- sfit$pivot
      sol <- backsolve(R_g, forwardsolve(t(R_g), rhs_r[piv_g])) * N
      ans <- numeric(gf$p)
      ans[piv_g] <- sol
      if (any(!is.finite(ans))) {
        rhs_mat <- matrix(rhs_r, nrow = gf$p, ncol = 1L)
        sol_m <- if (intercept) {
          or_pooled_wls_gram_solve(gf$fit, rhs_mat, normalizer = N)
        } else {
          or_wls_gram_solve(gf$fit, rhs_mat, normalizer = N)
        }
        ans <- as.double(sol_m[, 1L])
      }
      tilde_nu[[gid]] <- ans
    }
    return(reconstruction_pack_plain_adjoint(tilde_lambda, tilde_nu, structural_cids, cells))
  }
  active_groups <- unique(vapply(cells[structural_cids], `[[`, integer(1L), "group_id"))
  for (gid in active_groups) {
    gf <- groups[[gid]]
    rhs_r <- numeric(gf$p)
    for (idx in seq_len(n_sc)) {
      cid <- structural_cids[idx]
      cf <- cells[[cid]]
      if (cf$group_id == gid) {
        lam_c <- if (intercept) lam_vals[idx] else 0.0
        if (intercept) tilde_lambda[cid] <- lam_c / cf$mass
        if (gf$p > 0L) rhs_r <- rhs_r + nu_vals[[idx]] - lam_c * cf$mean_d
      }
    }
    if (gf$p > 0L) {
      sfit <- if (intercept) gf$fit$slope_fit else gf$fit
      R_g <- sfit$R
      piv_g <- sfit$pivot
      sol <- backsolve(R_g, forwardsolve(t(R_g), rhs_r[piv_g])) * N
      ans <- numeric(gf$p)
      ans[piv_g] <- sol
      if (any(!is.finite(ans))) {
        rhs_mat <- matrix(rhs_r, nrow = gf$p, ncol = 1L)
        sol_m <- if (intercept) {
          or_pooled_wls_gram_solve(gf$fit, rhs_mat, normalizer = N)
        } else {
          or_wls_gram_solve(gf$fit, rhs_mat, normalizer = N)
        }
        ans <- as.double(sol_m[, 1L])
      }
      tilde_nu[[gid]] <- ans
    }
  }
  reconstruction_pack_plain_adjoint(tilde_lambda, tilde_nu, structural_cids, cells)
}

cdid_ia_fast_or_pooled_fit <- function(recipe, supports, blocks, rank_tol, N, fk_comp = NULL) {
  if (!isTRUE(recipe$intercept)) return(NULL)
  if (is.null(fk_comp)) fk_comp <- cdid_fast_kernel("C_cdid_wls_compress")
  if (is.null(fk_comp)) return(NULL)
  p <- as.integer(recipe$p)
  n_b <- length(blocks)
  if (n_b < 1L) return(NULL)
  block_names <- vapply(blocks, `[[`, character(1L), "name")
  if (anyNA(block_names) || any(!nzchar(block_names)) || anyDuplicated(block_names)) return(NULL)
  block_recs <- vector("list", n_b)
  acc_fac <- NULL
  acc_n <- 0L
  max_aug <- 0L
  max_res <- 0L
  max_chunk_merge <- 0L
  max_res_merge <- 0L
  n_compressions <- 0L
  cols <- if (p > 0L) recipe$slope_columns else NULL
  keep <- seq_len(p)
  for (b in seq_len(n_b)) {
    sids <- blocks[[b]]$support_ids
    if (length(sids) < 1L) return(NULL)
    sup1 <- supports[[sids[1L]]]
    if (!identical(sup1$scale, 1.0)) return(NULL)
    chunk1 <- sup1$summary
    if (!inherits(chunk1, "cdid_or_pooled_wls_chunk")) return(NULL)
    fac <- chunk1$raw$factor
    c_n <- as.integer(chunk1$n)
    c_const <- chunk1$constant
    c_anchor <- chunk1$anchor
    c_merge <- chunk1$accounting$max_chunk_merge_rows
    c_comp <- chunk1$accounting$n_compressions
    if (length(sids) > 1L) {
      for (sid in sids[-1L]) {
        sup_k <- supports[[sid]]
        if (!identical(sup_k$scale, 1.0)) return(NULL)
        nxt <- sup_k$summary
        if (!inherits(nxt, "cdid_or_pooled_wls_chunk")) return(NULL)
        st_rows <- nrow(fac) + nrow(nxt$raw$factor)
        if (st_rows > c_merge) c_merge <- st_rows
        if (nxt$n > 0L) {
          if (c_n == 0L) {
            c_anchor <- nxt$anchor
            c_const <- nxt$constant
          } else {
            c_const <- c_const & nxt$constant & (c_anchor == nxt$anchor)
          }
        }
        st_mat <- rbind(fac, nxt$raw$factor)
        fac <- if (nrow(st_mat) > 0L) fk_comp(st_mat) else st_mat
        if (is.null(fac)) return(NULL)
        c_n <- c_n + as.integer(nxt$n)
        c_comp <- c_comp + nxt$accounting$n_compressions + if (st_rows > 0L) 1L else 0L
      }
    }
    m_fac <- nrow(fac)
    if (c_n < 1L || m_fac < 1L || ncol(fac) != p + 2L || any(!is.finite(fac))) return(NULL)
    is_upper_int <- (m_fac == 1L || all(fac[-1L, 1L] == 0.0))
    if (is_upper_int) {
      r11 <- fac[1L, 1L]
      if (!is.finite(r11) || abs(r11) <= 0.0) return(NULL)
      qty <- fac[, -1L, drop = FALSE]
    } else {
      q_int <- qr.default(fac[, 1L, drop = FALSE], LAPACK = TRUE)
      r11 <- q_int$qr[1L, 1L]
      if (!is.finite(r11) || abs(r11) <= 0.0) return(NULL)
      qty <- qr.qty(q_int, fac[, -1L, drop = FALSE])
      if (any(!is.finite(qty))) return(NULL)
    }
    r1d <- if (p > 0L) qty[1L, keep] else numeric(0L)
    names(r1d) <- cols
    r1y <- qty[1L, p + 1L]
    res_raw <- qty[-1L, , drop = FALSE]
    has_const <- (p > 0L && any(c_const))
    if (has_const) {
      res_raw[, which(c_const)] <- 0.0
    }
    res_fac <- if (nrow(res_raw) == 0L) {
      res_raw
    } else if (is_upper_int && !has_const && length(sids) == 1L && nrow(res_raw) == p + 1L &&
               all(res_raw[lower.tri(res_raw)] == 0.0) && all(diag(res_raw) >= 0.0)) {
      res_raw
    } else {
      fk_comp(res_raw)
    }
    if (is.null(res_fac)) return(NULL)
    m_res <- nrow(res_fac)
    block_recs[[b]] <- list(
      r11 = r11,
      r1d = r1d,
      r1y = r1y,
      n = c_n,
      weight_sum = r11 * r11,
      mean_d = r1d / r11,
      mean_y = r1y / r11,
      augmented_rows = m_fac,
      residual_rows = m_res
    )
    c_comp <- c_comp + if (nrow(res_raw) > 0L) 1L else 0L
    if (b == 1L) {
      acc_fac <- res_fac
      acc_n <- c_n
      max_aug <- m_fac
      max_res <- m_res
      max_chunk_merge <- c_merge
      max_res_merge <- m_res
      n_compressions <- c_comp
    } else {
      stacked <- nrow(acc_fac) + m_res
      if (stacked > max_res_merge) max_res_merge <- stacked
      stacked_mat <- rbind(acc_fac, res_fac)
      acc_fac <- if (nrow(stacked_mat) > 0L) fk_comp(stacked_mat) else stacked_mat
      if (is.null(acc_fac)) return(NULL)
      acc_n <- acc_n + c_n
      if (m_fac > max_aug) max_aug <- m_fac
      if (m_res > max_res) max_res <- m_res
      if (c_merge > max_chunk_merge) max_chunk_merge <- c_merge
      n_compressions <- n_compressions + c_comp + if (stacked > 0L) 1L else 0L
    }
  }
  names(block_recs) <- block_names
  if (p > 0L) {
    if (nrow(acc_fac) < p) return(NULL)
    X <- acc_fac[, keep, drop = FALSE]
    q <- qr.default(X, LAPACK = TRUE)
    R <- q$qr[keep, , drop = FALSE]
    if (p > 1L) R[lower.tri(R)] <- 0.0
    diagonal <- abs(diag(R))
    max_diag <- max(diagonal)
    cutoff <- rank_tol * max_diag
    rank <- if (cutoff > 0.0) sum(diagonal > cutoff) else 0L
    if (rank < p) return(NULL)
    if (min(diagonal) <= 100.0 * p * cutoff) {
      proposal <- qr.default(X, tol = rank_tol, LAPACK = FALSE)
      if (proposal$rank < p) return(NULL)
    }
    qty <- qr.qty(q, acc_fac[, p + 1L, drop = FALSE])[, 1L]
    coefficients <- numeric(p)
    coefficients[q$pivot] <- backsolve(R, qty[keep])
    names(coefficients) <- cols
    rss <- sum(qty[-keep]^2)
    R_out <- unname(R)
    pivot <- q$pivot
    if (any(!is.finite(c(coefficients, rss, R_out)))) return(NULL)
    slope_fit <- structure(
      list(
        coefficients = coefficients,
        rss = rss,
        R = R_out,
        pivot = pivot,
        n = acc_n,
        rank = rank,
        rank_tol = rank_tol,
        cutoff = cutoff
      ),
      class = "cdid_or_wls_fit"
    )
  } else {
    slope_fit <- NULL
    coefficients <- numeric(0L)
    if (!is.null(cols)) names(coefficients) <- cols
    rss <- sum(acc_fac[, 1L]^2)
    R_out <- matrix(numeric(0L), 0L, 0L)
    pivot <- integer(0L)
    rank <- 0L
    cutoff <- 0.0
    if (!is.finite(rss)) return(NULL)
  }
  intercepts <- numeric(n_b)
  regressions <- vector("list", n_b)
  beta <- coefficients
  for (b in seq_len(n_b)) {
    rec <- block_recs[[b]]
    rhs <- rec$r1y - sum(rec$r1d * coefficients)
    alpha_b <- rhs / rec$r11
    if (!is.finite(alpha_b)) return(NULL)
    intercepts[b] <- alpha_b
    regressions[[b]] <- list(
      name = block_names[b],
      L = blocks[[b]]$L,
      n = rec$n,
      alpha = alpha_b,
      weight_sum = rec$weight_sum,
      mass = rec$weight_sum / N,
      mean_d = unname(rec$mean_d),
      mean_y = rec$mean_y,
      residual_rows = rec$residual_rows
    )
  }
  names(intercepts) <- block_names
  acct <- list(
    n = acc_n,
    p = p,
    n_blocks = n_b,
    residual_dim = dim(acc_fac),
    max_augmented_rows = max_aug,
    max_residual_rows = max_res,
    max_chunk_merge_rows = max_chunk_merge,
    max_residual_merge_rows = max_res_merge,
    n_compressions = n_compressions
  )
  res_obj <- if (p > 0L) {
    structure(list(factor = acc_fac, n = acc_n, columns = cols), class = "cdid_or_wls_summary")
  } else {
    list(factor = acc_fac, n = acc_n, columns = cols)
  }
  pooled <- structure(
    list(
      residual = res_obj,
      blocks = block_recs,
      n = acc_n,
      p = p,
      columns = cols,
      accounting = acct
    ),
    class = "cdid_or_pooled_wls_summary"
  )
  fit <- structure(
    list(
      coefficients = coefficients,
      intercepts = intercepts,
      rss = rss,
      R = R_out,
      pivot = pivot,
      n = acc_n,
      p = p,
      n_blocks = n_b,
      rank = rank,
      rank_tol = rank_tol,
      cutoff = cutoff,
      blocks = block_recs,
      slope_fit = slope_fit,
      accounting = acct
    ),
    class = "cdid_or_pooled_wls_fit"
  )
  map <- matrix(0.0, p, 0L, dimnames = list(cols, character(0L)))
  constants <- matrix(0.0, n_b, 0L, dimnames = list(block_names, character(0L)))
  basis <- list(
    slope_columns = recipe$slope_columns,
    kept_columns = recipe$slope_columns[keep],
    kept_index = keep,
    proposed_columns = character(0L),
    aliased_columns = character(0L),
    nontrivial_candidates_retained = character(0L),
    alias_map = map,
    block_constants = constants,
    reconstruction_control_scaled = numeric(0L),
    reconstruction_treated_scaled = numeric(0L),
    reconstruction_roundoff_threshold = NA_real_,
    rows_checked = 0,
    proposal = "R stats::lm.fit / dqrdc2 pivoted QR at rank_tol on the pooled (partialled) factor",
    rank_check = "unchanged relative-diagonal rule of the accepted kernel on the reduced factor",
    frozen_during_differentiation = TRUE,
    numerical_certificate = FALSE
  )
  list(
    fit = fit,
    summary = pooled,
    keep = keep,
    beta = beta,
    regressions = regressions,
    basis = basis,
    rss = rss,
    objective = rss / N,
    verification_scan = FALSE
  )
}

cdid_ia_fast_or_outcome_fit_blocks <- function(recipe, supports, blocks, links, rank_tol, N,
                                               mode, route = "pooled") {
  if (!isTRUE(recipe$intercept) || !mode %in% c("pooled", "separate")) return(NULL)
  fk_comp <- cdid_fast_kernel("C_cdid_wls_compress")
  if (is.null(fk_comp)) return(NULL)
  B <- length(blocks)
  if (B < 1L) return(NULL)
  block_names <- vapply(blocks, `[[`, character(1L), "name")
  if (anyNA(block_names) || any(!nzchar(block_names)) || anyDuplicated(block_names)) return(NULL)
  if (mode == "pooled") {
    g1 <- cdid_ia_fast_or_pooled_fit(recipe, supports, blocks, rank_tol, N, fk_comp = fk_comp)
    if (is.null(g1)) return(NULL)
    groups <- list(g1)
    group_id <- rep.int(1L, B)
    local_block_id <- seq_len(B)
    group_map <- list(list(
      block_ids = seq_len(B),
      support_ids = seq_along(supports),
      link_ids = which(vapply(links, function(l) l$regression %in% seq_len(B), logical(1L))),
      slope_dimension_bound = recipe$p,
      coefficient_dimension_bound = recipe$p + B,
      operator_dimension_bound = c(recipe$p, recipe$p),
      summary_dimension_bound = c(recipe$p + 1L, recipe$p + 1L)
    ))
  } else {
    bids <- vapply(links, `[[`, integer(1L), "regression")
    lsup <- vapply(links, `[[`, integer(1L), "support_id")
    groups <- vector("list", B)
    group_map <- vector("list", B)
    l_names <- names(links)
    for (g in seq_len(B)) {
      lids <- which(bids == g)
      if (!is.null(l_names)) names(lids) <- l_names[lids]
      b_g <- blocks[[g]]
      sids <- sort.int(unique(c(as.integer(b_g$support_ids), lsup[lids])))
      local_b <- b_g
      local_b$support_ids <- match(b_g$support_ids, sids)
      gg <- cdid_ia_fast_or_pooled_fit(recipe, supports[sids], list(local_b), rank_tol, N, fk_comp = fk_comp)
      if (is.null(gg)) return(NULL)
      groups[[g]] <- gg
      group_map[[g]] <- list(
        block_ids = g,
        support_ids = as.integer(sids),
        link_ids = lids,
        slope_dimension_bound = recipe$p,
        coefficient_dimension_bound = recipe$p + 1L,
        operator_dimension_bound = c(recipe$p, recipe$p),
        summary_dimension_bound = c(recipe$p + 1L, recipe$p + 1L)
      )
    }
    group_id <- seq_len(B)
    local_block_id <- rep.int(1L, B)
  }
  block_map <- structure(
    list(
      block_id = seq_len(B),
      name = block_names,
      group_id = group_id,
      local_block_id = local_block_id
    ),
    class = "data.frame",
    row.names = .set_row_names(B)
  )
  rss <- sum(vapply(groups, `[[`, numeric(1L), "rss"))
  retained <- lapply(groups, function(g) list(
    slope_dimension = length(g$keep),
    operator_dimension = dim(g$fit$R),
    summary_dimension = dim(g$summary$residual$factor),
    alias_map_dimension = dim(g$basis$alias_map),
    block_constants_dimension = dim(g$basis$block_constants),
    verification_scan = g$verification_scan,
    alias_rows_checked = g$basis$rows_checked
  ))
  structure(
    list(
      mode = mode,
      recipe = recipe,
      N = N,
      route = route,
      block_map = block_map,
      group_map = group_map,
      groups = groups,
      rss = rss,
      objective = rss / N,
      retained = retained
    ),
    class = "cdid_or_outcome_groups"
  )
}

cdid_ia_fast_fit_selection_models <- function(reader, panel, link_views, selection_model,
                                              selection_pooling, ctrl, gc,
                                              retained_external_entries, external_work_ops) {
  if (is.null(panel) || isTRUE(gc$trace)) return(NULL)
  fk_sel <- cdid_fast_kernel("C_cdid_ia_fit_sel_blocks")
  if (is.null(fk_sel)) return(NULL)
  st <- reader$state
  validate_selection_mode_args(selection_model, selection_pooling, st$universe, st$U, st$rho)
  L <- length(link_views)
  if (L < 1L) return(NULL)
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries + retained_external_entries
  d_cols <- ncol(st$Z)
  design_B <- min(ctrl$max_chunk_rows, st$N)
  cp04_budget(ctrl, base, 32 * design_B * (d_cols + 4),
              c(32 * st$N * (d_cols + 4), external_work_ops), 0,
              "selection design scan")
  design_spec <- extract_selection_design_spec(st$Z, st$U, has_intercept = NULL,
                                               chunk_rows = ctrl$max_chunk_rows)
  if (!design_spec$has_intercept && design_spec$d_s == 0L) return(NULL)
  header_bound <- 2048 + 1024 * L + sum(vapply(link_views, function(v) {
    if (is.list(v) && is.list(v$cohorts)) 128 * length(v$cohorts$cohorts) else 0
  }, numeric(1L)))
  cp04_budget(ctrl, c(base, header_bound), header_bound, 64 * L * L,
              header_bound, "selection model catalog metadata")

  cat_obj <- build_selection_block_catalog(
    reader = reader,
    link_views = link_views,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    control = ctrl,
    reference = NULL,
    view_cohorts = NULL,
    families = NULL
  )
  cells <- cat_obj$cells
  block_specs <- cat_obj$block_specs
  pair_cohorts <- cat_obj$pair_cohorts
  links <- cat_obj$links
  C <- length(cells)
  M <- length(block_specs)
  d_s <- design_spec$d_s
  d_full <- design_spec$d
  has_int <- design_spec$has_intercept

  metadata_payload <- cp04_payload(list(cells, block_specs, pair_cohorts, links, design_spec))
  metadata_entries <- unname(metadata_payload["entries"])
  model_bound <- sum(vapply(block_specs, function(bs) {
    B_m <- bs$B_m
    B_m * (d_s + 8) + (d_s + 2)^2 + 4096
  }, numeric(1L))) + C * 64 + 4096
  resources <- selection_fit_resources(st, block_specs, cells, design_spec, gc, ctrl,
                                       selection_model, selection_pooling)
  scratch <- resources$scratch + header_bound
  sel_work <- resources$work + 32 * st$N * (d_full + 4)
  work <- sel_work + 64 * (L * L + unname(metadata_payload["char_bytes"]))
  retained_bound <- base + metadata_entries + model_bound + 4096
  cp04_budget(ctrl, retained_bound, scratch, c(work, external_work_ops),
              metadata_entries + model_bound + 4096, "selection models fit")

  is_date <- identical(selection_model, "independent_period_logit")
  is_sep_pair <- identical(selection_model, "pair_logit") &&
    identical(selection_pooling, "separate") && d_full > 0L
  raw_blocks <- fk_sel(
    st, panel, cells, block_specs, design_spec,
    c(if (is_date) 1.0 else 0.0, if (is_sep_pair) 1.0 else 0.0, gc$epsilon, as.double(gc$maxit))
  )
  if (is.null(raw_blocks)) return(NULL)

  cn <- colnames(st$Z)
  dn_chol <- if (!is.null(cn)) list(cn, cn) else NULL
  coord_ord <- c(design_spec$int_col, design_spec$slope_cols)
  empty_warn <- character(0L)
  blocks <- vector("list", M)
  for (m in seq_len(M)) {
    bs <- block_specs[[m]]
    cids <- as.integer(bs$cell_ids)
    B_m <- length(cids)
    rb <- raw_blocks[[m]]
    if (is_sep_pair && B_m == 1L) {
      coef <- rb[[1L]]
      H_chol <- rb[[2L]]
      if (!is.null(dn_chol)) dimnames(H_chol) <- dn_chol
      blocks[[m]] <- list(
        branch = "cp04_separate_pair",
        has_intercept = has_int,
        B_m = 1L,
        d_s = d_s,
        alpha = if (has_int) coef[design_spec$int_col] else numeric(0L),
        beta = coef[design_spec$slope_cols],
        coefficients = coef,
        H_chol = H_chol,
        coordinate_order = coord_ord,
        m_b = rb[[6L]],
        z_bar = rb[[7L]],
        slope_fit = NULL,
        pi = rb[[4L]],
        r = rb[[5L]],
        diagnostics = list(
          converged = TRUE,
          boundary = FALSE,
          rank = length(coef),
          iter = rb[[3L]],
          warnings = empty_warn
        ),
        block_id = bs$block_id,
        cell_ids = cids
      )
    } else {
      alpha <- rb[[1L]]
      beta <- rb[[2L]]
      p_m <- (if (has_int) B_m else 0L) + d_s
      slope_fit <- if (d_s > 0L) {
        structure(
          list(
            coefficients = numeric(d_s),
            rss = 0,
            R = rb[[5L]],
            pivot = seq_len(d_s),
            n = rb[[9L]],
            p = d_s,
            rank = d_s,
            rank_tol = rb[[10L]],
            cutoff = rb[[11L]]
          ),
          class = "cdid_or_wls_fit"
        )
      } else {
        NULL
      }
      blocks[[m]] <- list(
        branch = design_spec$branch,
        has_intercept = has_int,
        B_m = B_m,
        d_s = d_s,
        alpha = alpha,
        beta = beta,
        coefficients = c(alpha, beta),
        H_chol = NULL,
        m_b = rb[[3L]],
        z_bar = rb[[4L]],
        slope_fit = slope_fit,
        pi = rb[[6L]],
        r = rb[[7L]],
        diagnostics = list(
          converged = TRUE,
          boundary = FALSE,
          rank = p_m,
          iter = rb[[8L]],
          warnings = empty_warn
        ),
        block_id = bs$block_id,
        cell_ids = cids
      )
    }
  }

  alpha_all <- if (has_int) unlist(lapply(blocks, `[[`, "alpha"), use.names = FALSE) else numeric(0L)
  beta_all <- unlist(lapply(blocks, `[[`, "beta"), use.names = FALSE)
  theta_all <- c(alpha_all, beta_all)
  law_str <- paste0(
    "CP04:", st$universe, ":weighted-", selection_model, ":", selection_pooling, ":a_U/(rho*p_U)"
  )
  obj <- structure(
    list(
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      design_spec = design_spec,
      cells = cells,
      block_specs = block_specs,
      blocks = blocks,
      pair_cohorts = pair_cohorts,
      links = links,
      alpha = alpha_all,
      beta = beta_all,
      coefficients = theta_all,
      glm_control = gc,
      law = law_str,
      universe = st$universe,
      transport_justification = st$justification,
      reader_ref = reader$ref,
      control = ctrl,
      accounting = list(
        owned_entries = 0,
        character_bytes = 0,
        borrowed_entries = base,
        retained_bound = retained_bound,
        fit_scratch_bound = scratch,
        work_ops_bound = work,
        n_active_cells = C,
        n_selection_blocks = M,
        n_pair_cohorts = length(pair_cohorts),
        n_parameters = length(theta_all),
        full_F_masks_per_model = 0L,
        full_F_by_cell_IF_retained = FALSE,
        global_dummy_matrix_constructed = FALSE,
        inverse_matrices_retained = FALSE,
        peak_rss_measured = FALSE
      )
    ),
    class = "cdid_selection_models"
  )
  obj["ref"] <- list(NULL)
  payload <- cp04_payload(obj)
  actual <- unname(payload["entries"])
  if (actual > metadata_entries + model_bound + 4096) {
    cp04_abort("Retained selection model bound underestimated.", "bound")
  }
  obj$accounting$owned_entries <- actual
  obj$accounting$character_bytes <- unname(payload["char_bytes"])
  cdid_fast_seal(obj, "selection_models")
}


