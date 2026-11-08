# Independent test-only reference for test-static-observation-adapter.R.
# Uses stats::glm.fit / stats::lm.wfit and literal ratio/pseudo-change formulas
# for independent verification of coefficients, moment estimates, and
# unit-weight-perturbation finite-difference original-unit influence functions.

adapter_ref_partition <- function(pooling, group_cell_map) {
  all_cells <- unlist(group_cell_map, use.names = FALSE)
  if (pooling == "separate") {
    out <- as.list(all_cells)
    names(out) <- all_cells
    return(out)
  }
  if (pooling == "within_group") {
    return(group_cell_map)
  }
  list(pooled = all_cells)
}

adapter_ref_build_specs <- function(N, groups, route, selection_model,
                                    selection_pooling, imputation_pooling,
                                    imputation_scale, selection_design = NULL,
                                    imputation_design = NULL, dates = NULL) {
  uses_q <- route %in% c("weighting", "augmented")
  uses_j <- route %in% c("imputation", "augmented")
  sel_cells <- list()
  sel_group_cells <- list()
  imp_cells <- list()
  imp_group_cells <- list()

  for (g_key in names(groups)) {
    grp <- groups[[g_key]]
    A <- as.logical(grp$mask)
    inA <- which(A)
    g_sel <- if (!is.null(grp$selection_design)) grp$selection_design else selection_design
    g_imp <- if (!is.null(grp$imputation_design)) grp$imputation_design else imputation_design
    g_p_cells <- character(0)
    req_dates <- character(0)

    for (p_key in names(grp$pairs)) {
      pr <- grp$pairs[[p_key]]
      ck <- paste0(g_key, "::", p_key)
      g_p_cells <- c(g_p_cells, ck)
      req_dates <- c(req_dates, as.character(pr$start), as.character(pr$end))
      raw_obs <- if (!is.null(pr$observed)) pr$observed else pr$O
      obs_log <- logical(N)
      if (!is.null(raw_obs)) obs_log[inA] <- as.double(raw_obs[inA]) == 1
      obs_num <- as.double(obs_log)

      if (uses_q && identical(selection_model, "pair_logit")) {
        X_s <- if (!is.null(pr$selection_design)) pr$selection_design else g_sel
        sel_cells[[ck]] <- list(mask = A, response = obs_num, design = X_s)
      }
      if (uses_j) {
        X_m <- if (!is.null(pr$imputation_design)) pr$imputation_design else g_imp
        if (!is.null(X_m) && ncol(X_m) > 0L && identical(imputation_scale, "per_period")) {
          X_scaled <- matrix(NA_real_, N, ncol(X_m), dimnames = dimnames(X_m))
          X_scaled[inA, ] <- as.integer(pr$gap) * X_m[inA, , drop = FALSE]
          X_m <- X_scaled
        }
        dy_vec <- rep(NA_real_, N)
        R_o <- which(A & obs_log)
        dy_vec[R_o] <- as.double(pr$dy[R_o])
        imp_cells[[ck]] <- list(mask = A, observed = obs_log, dy = dy_vec, design = X_m)
      }
    }
    if (uses_q && identical(selection_model, "pair_logit")) {
      sel_group_cells[[g_key]] <- g_p_cells
    }
    if (uses_j) {
      imp_group_cells[[g_key]] <- g_p_cells
    }
    if (uses_q && identical(selection_model, "independent_period_logit")) {
      d_src <- if (!is.null(grp$dates)) grp$dates else dates
      d_keys <- if (!is.null(grp$dates)) names(grp$dates) else names(d_src)[names(d_src) %in% unique(req_dates)]
      g_d_cells <- character(0)
      for (d_key in d_keys) {
        dr <- d_src[[d_key]]
        raw_d <- if (!is.null(dr$observed)) dr$observed else dr$response
        d_num <- numeric(N)
        d_num[inA] <- as.double(raw_d[inA])
        X_d <- if (!is.null(dr$selection_design)) dr$selection_design else
          if (!is.null(dr$design)) dr$design else g_sel
        d_ck <- paste0(g_key, "::", d_key)
        g_d_cells <- c(g_d_cells, d_ck)
        sel_cells[[d_ck]] <- list(mask = A, response = d_num, design = X_d)
      }
      sel_group_cells[[g_key]] <- g_d_cells
    }
  }

  list(
    sel_cells = sel_cells,
    sel_partition = if (uses_q) adapter_ref_partition(selection_pooling, sel_group_cells) else NULL,
    imp_cells = imp_cells,
    imp_partition = if (uses_j) adapter_ref_partition(imputation_pooling, imp_group_cells) else NULL
  )
}

adapter_ref_fit_block <- function(cells, keys, kind, w_unit) {
  parts <- lapply(seq_along(keys), function(b) {
    cl <- cells[[keys[b]]]
    rows <- if (kind == "logit") which(cl$mask) else which(cl$mask & cl$observed)
    X <- if (is.null(cl$design) || ncol(cl$design) == 0L) {
      matrix(numeric(0), length(rows), 0L)
    } else {
      unname(cl$design[rows, , drop = FALSE])
    }
    y <- if (kind == "logit") as.double(cl$response[rows]) else as.double(cl$dy[rows])
    list(unit = rows, cell = rep.int(b, length(rows)), y = y, X = X)
  })
  unit <- unlist(lapply(parts, `[[`, "unit"))
  cell <- unlist(lapply(parts, `[[`, "cell"))
  y <- unlist(lapply(parts, `[[`, "y"))
  X <- do.call(rbind, lapply(parts, `[[`, "X"))
  B <- length(keys)
  Xf <- cbind(outer(cell, seq_len(B), "==") * 1, X)
  ww <- w_unit[unit]
  coef <- if (kind == "logit") {
    f <- suppressWarnings(stats::glm.fit(
      Xf, y, weights = ww, family = stats::binomial(),
      control = stats::glm.control(epsilon = 1e-14, maxit = 200)
    ))
    stopifnot(f$converged)
    unname(f$coefficients)
  } else {
    unname(stats::lm.wfit(Xf, y, ww)$coefficients)
  }
  list(alpha = stats::setNames(coef[seq_len(B)], keys),
       beta = if (ncol(X) > 0L) coef[B + seq_len(ncol(X))] else numeric(0))
}

adapter_ref_fit_all <- function(cells, partition, kind, w_unit) {
  blocks <- lapply(partition, function(keys) adapter_ref_fit_block(cells, keys, kind, w_unit))
  coef <- unlist(lapply(blocks, function(b) c(unname(b$alpha), b$beta)), use.names = FALSE)
  list(blocks = blocks, coefficients = coef)
}

adapter_ref_predict_cell <- function(cells, partition, fit_obj, cell_key, rows, kind) {
  b_name <- names(partition)[vapply(partition, function(ks) cell_key %in% ks, logical(1))][1L]
  bf <- fit_obj$blocks[[b_name]]
  cl <- cells[[cell_key]]
  alpha <- bf$alpha[[cell_key]]
  eta <- if (is.null(cl$design) || ncol(cl$design) == 0L) {
    rep(alpha, length(rows))
  } else {
    alpha + drop(unname(cl$design[rows, , drop = FALSE]) %*% bf$beta)
  }
  if (kind == "logit") stats::plogis(eta) else eta
}

# Literal ratio / pseudo-change functional for one edge under unit weights w_unit.
adapter_ref_edge_functional <- function(N, groups, group_key, pair_key, route,
                                        selection_model = NULL, selection_pooling = NULL,
                                        imputation_pooling = NULL, imputation_scale = NULL,
                                        selection_design = NULL, imputation_design = NULL,
                                        dates = NULL, kind = "mean", h_vec = NULL,
                                        h_fitter = NULL, Z = NULL, w_unit = rep.int(1, N),
                                        floor_eps = 1e-6) {
  specs <- adapter_ref_build_specs(
    N = N, groups = groups, route = route, selection_model = selection_model,
    selection_pooling = selection_pooling, imputation_pooling = imputation_pooling,
    imputation_scale = imputation_scale, selection_design = selection_design,
    imputation_design = imputation_design, dates = dates
  )
  uses_q <- route %in% c("weighting", "augmented")
  uses_j <- route %in% c("imputation", "augmented")
  sel_fit <- if (uses_q) adapter_ref_fit_all(specs$sel_cells, specs$sel_partition, "logit", w_unit) else NULL
  imp_fit <- if (uses_j) adapter_ref_fit_all(specs$imp_cells, specs$imp_partition, "wls", w_unit) else NULL

  grp <- groups[[group_key]]
  pr <- grp$pairs[[pair_key]]
  A <- as.logical(grp$mask)
  inA <- which(A)
  ck <- paste0(group_key, "::", pair_key)

  raw_obs <- if (!is.null(pr$observed)) pr$observed else pr$O
  O_inA <- if (route != "imputation") (as.double(raw_obs[inA]) == 1) else rep(FALSE, length(inA))
  R_o <- inA[O_inA]
  R_w <- if (route %in% c("baseline", "weighting")) R_o else inA

  h_use <- if (kind == "mean") {
    if (!is.null(h_fitter)) h_fitter(w_unit)[R_w] else as.double(h_vec[R_w])
  } else {
    rep.int(1, length(R_w))
  }
  w_emp <- w_unit[R_w]

  phi_w <- numeric(length(R_w))
  if (uses_q && length(R_o)) {
    q_o <- if (identical(selection_model, "pair_logit")) {
      adapter_ref_predict_cell(specs$sel_cells, specs$sel_partition, sel_fit, ck, R_o, "logit")
    } else {
      ck_s <- paste0(group_key, "::", as.character(pr$start))
      ck_u <- paste0(group_key, "::", as.character(pr$end))
      q_s <- adapter_ref_predict_cell(specs$sel_cells, specs$sel_partition, sel_fit, ck_s, R_o, "logit")
      q_u <- adapter_ref_predict_cell(specs$sel_cells, specs$sel_partition, sel_fit, ck_u, R_o, "logit")
      q_s * q_u
    }
    phi_o <- 1 / pmax(q_o, floor_eps)
    phi_w[match(R_o, R_w)] <- phi_o
  }

  j_w <- if (uses_j) {
    adapter_ref_predict_cell(specs$imp_cells, specs$imp_partition, imp_fit, ck, R_w, "wls")
  } else NULL

  dy_o <- if (route != "imputation" && length(R_o)) as.double(pr$dy[R_o]) else numeric(0)

  if (route == "baseline") {
    w <- w_emp * h_use
    Y <- dy_o
  } else if (route == "weighting") {
    w <- w_emp * h_use * phi_w
    Y <- dy_o
  } else if (route == "imputation") {
    w <- w_emp * h_use
    Y <- j_w
  } else {
    resid_w <- numeric(length(R_w))
    if (length(R_o)) {
      o_pos <- match(R_o, R_w)
      resid_w[o_pos] <- dy_o - j_w[o_pos]
    }
    w <- w_emp * h_use
    Y <- j_w + phi_w * resid_w
  }

  S <- sum(w)
  M <- sum(w * Y) / S
  if (kind == "mean") return(M)

  Z_use <- if (!is.null(Z)) Z else grp$Z
  Zm <- if (is.matrix(Z_use)) Z_use else matrix(as.double(Z_use), ncol = 1L)
  Zw <- Zm[R_w, , drop = FALSE]
  Zbar <- colSums(w * Zw) / S
  Zc <- sweep(Zw, 2L, Zbar, "-")
  drop(crossprod(Zc, w * Y)) / S
}

adapter_ref_fd_if <- function(units, eval_fn, N, eps = 1e-5) {
  f0 <- eval_fn(rep.int(1, N))
  m <- length(f0)
  out <- matrix(0, nrow = length(units), ncol = m)
  for (idx in seq_along(units)) {
    i <- units[idx]
    wp <- wm <- rep.int(1, N)
    wp[i] <- 1 + eps
    wm[i] <- 1 - eps
    out[idx, ] <- N * (eval_fn(wp) - eval_fn(wm)) / (2 * eps)
  }
  if (m == 1L) drop(out) else out
}
