# Independent test-only reference implementation for R/static_missing_fit.R.
# Uses only base R and stats (stats::glm.fit, stats::lm.wfit, stats::model.matrix)
# and calls ZERO production cdid symbols.

ref_smf_incidence <- function(edges_df, baseline_time, target_time, dates = NULL) {
  s <- as.double(edges_df$s)
  u <- as.double(edges_df$u)
  verts <- sort(unique(c(s, u, as.double(dates))))
  # Independent Boolean transitive closure, not production forest/gauge helpers.
  adj <- diag(TRUE, length(verts))
  for (k in seq_along(s)) {
    si <- match(s[k], verts)
    ui <- match(u[k], verts)
    adj[si, ui] <- adj[ui, si] <- TRUE
  }
  repeat {
    next_adj <- (adj %*% adj) > 0
    if (identical(next_adj, adj)) break
    adj <- next_adj
  }
  endpoint_vertices <- which(adj[match(as.double(baseline_time), verts), ])
  assigned <- logical(length(verts))
  assigned[endpoint_vertices] <- TRUE
  references <- as.double(baseline_time)
  for (v in seq_along(verts)) {
    if (!assigned[v]) {
      members <- which(adj[v, ])
      assigned[members] <- TRUE
      references <- c(references, min(verts[members]))
    }
  }
  non_ground <- verts[!(verts %in% references)]
  K <- length(s)
  V <- length(non_ground)
  D <- matrix(0.0, nrow = K, ncol = V,
              dimnames = list(paste0(s, "->", u), as.character(non_ground)))
  for (k in seq_len(K)) {
    if (!(s[k] %in% references)) D[k, match(s[k], non_ground)] <- -1.0
    if (!(u[k] %in% references)) D[k, match(u[k], non_ground)] <- 1.0
  }
  ell <- numeric(V)
  names(ell) <- as.character(non_ground)
  ell[match(as.double(target_time), non_ground)] <- 1.0
  potential <- drop(solve(crossprod(D), ell))
  potential[!(non_ground %in% verts[endpoint_vertices])] <- 0.0
  varpi <- drop(D %*% potential)
  names(varpi) <- rownames(D)
  list(D = D, ell = ell, varpi = varpi, vertices = verts, non_ground = non_ground)
}

ref_smf_extract_design <- function(panel, t_prev_idx, formula_or_design, mask, cell_id = NULL) {
  n <- panel$n
  uid_chr <- as.character(panel$unit_id)
  active_idx <- which(as.logical(mask))

  if (inherits(formula_or_design, "formula")) {
    row_idx <- panel$row_index[, t_prev_idx]
    act_rows <- row_idx[active_idx]
    sub_df <- panel$data[act_rows, , drop = FALSE]
    mm <- stats::model.matrix(formula_or_design, data = sub_df)
    cn <- colnames(mm)
    X_tilde <- matrix(NA_real_, nrow = n, ncol = ncol(mm), dimnames = list(uid_chr, cn))
    X_tilde[active_idx, ] <- mm
    Z <- X_tilde[, -1L, drop = FALSE]
    return(list(X_tilde = X_tilde, Z = Z))
  }

  des_c <- if (is.list(formula_or_design) && !is.data.frame(formula_or_design) &&
               !is.null(cell_id) && (cell_id %in% names(formula_or_design))) {
    formula_or_design[[cell_id]]
  } else {
    formula_or_design
  }

  if (is.list(des_c) && (!is.null(des_c$X_tilde) || !is.null(des_c$Z))) {
    Z <- if (!is.null(des_c$Z)) des_c$Z else des_c$X_tilde[, -1L, drop = FALSE]
    X_tilde <- if (!is.null(des_c$X_tilde)) des_c$X_tilde else {
      out_X <- matrix(NA_real_, nrow = n, ncol = 1L + ncol(Z),
                      dimnames = list(uid_chr, c("(Intercept)", colnames(Z))))
      out_X[active_idx, 1L] <- 1.0
      if (ncol(Z) > 0L) out_X[active_idx, -1L] <- Z[active_idx, , drop = FALSE]
      out_X
    }
    return(list(X_tilde = X_tilde, Z = Z))
  }

  mat <- as.matrix(des_c)
  storage.mode(mat) <- "double"
  cn <- colnames(mat)
  has_int <- !is.null(cn) && length(cn) >= 1L && cn[1L] == "(Intercept)"
  if (has_int) {
    X_tilde <- mat
    Z <- mat[, -1L, drop = FALSE]
  } else {
    Z <- mat
    if (ncol(Z) > 0L && is.null(colnames(Z))) colnames(Z) <- paste0("z", seq_len(ncol(Z)))
    X_tilde <- matrix(NA_real_, nrow = n, ncol = 1L + ncol(Z),
                      dimnames = list(uid_chr, c("(Intercept)", colnames(Z))))
    X_tilde[active_idx, 1L] <- 1.0
    if (ncol(Z) > 0L) X_tilde[active_idx, -1L] <- Z[active_idx, , drop = FALSE]
  }
  list(X_tilde = X_tilde, Z = Z)
}

ref_smf_fit_treatment_logit <- function(R, G, X_tilde, unit_weights, glm_control) {
  n <- length(R)
  n_num <- as.double(n)
  w_u <- as.double(unit_weights)
  idx <- which(R == 1L)
  X_r <- X_tilde[idx, , drop = FALSE]
  y_r <- as.double(G[idx])
  w_r <- w_u[idx]

  fit <- suppressWarnings(stats::glm.fit(
    x = X_r,
    y = y_r,
    weights = w_r,
    family = stats::binomial(link = "logit"),
    control = stats::glm.control(
      epsilon = glm_control$epsilon,
      maxit = glm_control$maxit,
      trace = FALSE
    )
  ))
  gamma_hat <- as.double(fit$coefficients)
  names(gamma_hat) <- colnames(X_tilde)

  eta_r <- drop(X_r %*% gamma_hat)
  p_r <- stats::plogis(eta_r)
  W_r <- w_r * p_r * (1.0 - p_r)
  # Normalized by sum(w_u) for unit-weight perturbation compatibility (equals n when w_u = 1)
  W_sum <- sum(w_u)
  H_hat <- crossprod(X_r, X_r * W_r) / W_sum
  H_inv <- solve(H_hat)

  score_mat <- matrix(0.0, nrow = n, ncol = ncol(X_tilde),
                      dimnames = list(rownames(X_tilde), colnames(X_tilde)))
  score_mat[idx, ] <- X_r * (w_r * (y_r - p_r))
  ell_gamma <- (n_num / W_sum) * (score_mat %*% H_inv)
  dimnames(ell_gamma) <- list(rownames(X_tilde), colnames(X_tilde))

  p_hat <- rep(0.0, n)
  p_hat[idx] <- p_r
  p_bar <- 1.0 - 1e-6
  p_tilde <- pmin(p_hat, p_bar)
  r_hat <- rep(0.0, n)
  r_hat[idx] <- p_tilde[idx] / (1.0 - p_tilde[idx])
  active <- (R == 1L) & (p_hat < p_bar)

  # Derivative dr/dgamma on R == 1: where active, r_i * X_tilde_i
  dr_dgamma <- matrix(0.0, nrow = n, ncol = ncol(X_tilde),
                      dimnames = list(rownames(X_tilde), colnames(X_tilde)))
  act_idx <- which(active)
  if (length(act_idx) > 0L) {
    dr_dgamma[act_idx, ] <- X_tilde[act_idx, , drop = FALSE] * r_hat[act_idx]
  }

  list(
    gamma_hat = gamma_hat,
    p_hat = p_hat,
    r_hat = r_hat,
    active = active,
    H_hat = H_hat,
    ell_gamma = ell_gamma,
    dr_dgamma = dr_dgamma
  )
}

ref_smf_partition_blocks <- function(pooling, group_cell_map) {
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

ref_smf_fit_nuisance_block <- function(cells, keys, kind, w_unit, glm_control) {
  n <- length(w_unit)
  W_sum <- sum(w_unit)
  n_num <- as.double(n)
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
  Xf <- cbind(outer(cell, seq_len(B), "==") * 1.0, X)
  ww <- w_unit[unit]
  p_tot <- ncol(Xf)

  if (kind == "logit") {
    f <- suppressWarnings(stats::glm.fit(
      x = Xf, y = y, weights = ww, family = stats::binomial(link = "logit"),
      control = stats::glm.control(
        epsilon = glm_control$epsilon,
        maxit = glm_control$maxit,
        trace = FALSE
      )
    ))
    coef <- as.double(f$coefficients)
    eta <- drop(Xf %*% coef)
    q_fit <- stats::plogis(eta)
    wt_hess <- ww * q_fit * (1.0 - q_fit)
    G_mat <- crossprod(Xf, Xf * wt_hess) / W_sum
    resid <- ww * (y - q_fit)
  } else {
    f <- stats::lm.wfit(x = Xf, y = y, w = ww)
    coef <- as.double(f$coefficients)
    fit_val <- drop(Xf %*% coef)
    G_mat <- crossprod(Xf, Xf * ww) / W_sum
    resid <- ww * (y - fit_val)
  }

  G_inv <- solve(G_mat)
  # Aggregate unit scores across stacked rows belonging to each original unit i in 1..n
  U_score <- matrix(0.0, nrow = n, ncol = p_tot)
  row_scores <- Xf * resid
  for (r in seq_along(unit)) {
    u_i <- unit[r]
    U_score[u_i, ] <- U_score[u_i, ] + row_scores[r, ]
  }
  coef_if <- (n_num / W_sum) * (U_score %*% G_inv)

  list(
    keys = keys,
    alpha = stats::setNames(coef[seq_len(B)], keys),
    beta = if (ncol(X) > 0L) coef[B + seq_len(ncol(X))] else numeric(0L),
    coef = coef,
    coef_if = coef_if
  )
}

ref_smf_predict_nuisance <- function(cells, partition, block_fits, cell_key, mask_rows, n, kind) {
  b_name <- names(partition)[vapply(partition, function(ks) cell_key %in% ks, logical(1L))][1L]
  bf <- block_fits[[b_name]]
  cl <- cells[[cell_key]]
  B <- length(bf$keys)
  pos_c <- match(cell_key, bf$keys)
  p_tot <- length(bf$coef)
  X_sub <- if (is.null(cl$design) || ncol(cl$design) == 0L) {
    matrix(numeric(0), length(mask_rows), 0L)
  } else {
    unname(cl$design[mask_rows, , drop = FALSE])
  }
  eta <- bf$alpha[[cell_key]] + if (ncol(X_sub) > 0L) drop(X_sub %*% bf$beta) else 0.0

  X_full_block <- matrix(0.0, nrow = length(mask_rows), ncol = p_tot)
  X_full_block[, pos_c] <- 1.0
  if (ncol(X_sub) > 0L) {
    X_full_block[, B + seq_len(ncol(X_sub))] <- X_sub
  }

  val_vec <- rep(NA_real_, n)
  jac_mat <- matrix(0.0, nrow = n, ncol = p_tot)
  if (kind == "logit") {
    q_sub <- stats::plogis(eta)
    val_vec[mask_rows] <- q_sub
    jac_mat[mask_rows, ] <- X_full_block * (q_sub * (1.0 - q_sub))
  } else {
    val_vec[mask_rows] <- eta
    jac_mat[mask_rows, ] <- X_full_block
  }
  list(block = b_name, value = val_vec, jac = jac_mat)
}

# Full end-to-end independent reference for static_missing_fit
ref_static_missing_fit <- function(panel,
                                   target_cells,
                                   role_graphs,
                                   method,
                                   route,
                                   or_blocks = NULL,
                                   all_cells = NULL,
                                   formula = NULL,
                                   treatment_formula = NULL,
                                   outcome_formula = NULL,
                                   selection_formula = NULL,
                                   imputation_formula = NULL,
                                   treatment_design = NULL,
                                   outcome_design = NULL,
                                   selection_design = NULL,
                                   imputation_design = NULL,
                                   selection_model = NULL,
                                   selection_pooling = NULL,
                                   imputation_pooling = NULL,
                                   imputation_scale = NULL,
                                   weighting = "optimal",
                                   glm_control = list(epsilon = 1e-12, maxit = 100L, trace = FALSE),
                                   rank_tol = 1e-10,
                                   condition_trigger = 1e-10,
                                   ridge_relative = 1e-4,
                                   centered = FALSE,
                                   cluster = NULL,
                                   aggregate_cells = NULL,
                                   aggregate_signs = NULL,
                                   unit_weights = NULL,
                                   fixed_W = NULL) {
  n <- panel$n
  n_num <- as.double(n)
  uid <- panel$unit_id
  uid_chr <- as.character(uid)
  w_u <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  W_sum <- sum(w_u)

  if (identical(method, "or")) method <- "reg"
  use_logit <- method %in% c("ipw", "dr")
  use_or <- method %in% c("reg", "dr")
  uses_q <- route %in% c("weighting", "augmented")
  uses_j <- route %in% c("imputation", "augmented")

  # 1. Resolve cells in A and U
  t_map <- if (!is.null(target_cells$cells)) target_cells$cells else target_cells
  target_ids <- names(t_map)
  c_map <- t_map
  if (use_or && !is.null(all_cells)) {
    a_map <- if (!is.null(all_cells$cells)) all_cells$cells else all_cells
    for (cid in names(a_map)) c_map[[cid]] <- a_map[[cid]]
  }
  all_ids <- names(c_map)
  non_target_ids <- setdiff(all_ids, target_ids)

  trt_form <- if (!is.null(treatment_formula)) treatment_formula else formula
  or_form <- if (!is.null(outcome_formula)) outcome_formula else formula
  sel_form <- if (!is.null(selection_formula)) selection_formula else formula
  imp_form <- if (!is.null(imputation_formula)) imputation_formula else formula

  # 2. Extract designs and fit treatment logits on R_c == 1 for c in A
  trt_fits <- list()
  for (cid in target_ids) {
    cl <- c_map[[cid]]
    if (use_logit) {
      des_t <- ref_smf_extract_design(
        panel, cl$t_prev_idx,
        if (!is.null(treatment_design)) treatment_design else trt_form,
        cl$R == 1L, cid
      )
      X_safe <- des_t$X_tilde
      X_safe[cl$R == 0L, ] <- 0.0
      trt_fits[[cid]] <- ref_smf_fit_treatment_logit(cl$R, cl$G, X_safe, w_u, glm_control)
    }
  }

  or_Z_by_cell <- list()
  if (use_or) {
    for (cid in all_ids) {
      cl <- c_map[[cid]]
      des_o <- ref_smf_extract_design(
        panel, cl$t_prev_idx,
        if (!is.null(outcome_design)) outcome_design else or_form,
        cl$R == 1L, cid
      )
      or_Z_by_cell[[cid]] <- des_o$Z
    }
  }

  # 3. Build observation groups and fit selection / imputation models once across all declared groups
  sel_cells <- list()
  sel_group_cells <- list()
  imp_cells <- list()
  imp_group_cells <- list()
  group_edges <- list()

  add_group <- function(cid, role) {
    cl <- c_map[[cid]]
    g_key <- paste0(cid, "__", role)
    mask_log <- if (role == "G") (cl$G == 1L) else (cl$C == 1L)
    inA <- which(mask_log)

    rg_entry <- role_graphs[[cid]]
    spec_rc <- if (inherits(rg_entry, "cdid_static_role_block_spec")) {
      rg_entry
    } else if (is.list(rg_entry) && !is.data.frame(rg_entry) && (!is.null(rg_entry[[role]]))) {
      rg_entry[[role]]
    } else {
      rg_entry
    }
    edges_df <- if (!is.null(spec_rc$edges)) spec_rc$edges else if (!is.null(spec_rc$pairs)) spec_rc$pairs else as.data.frame(spec_rc)
    edges_df <- data.frame(
      link_id = if ("link_id" %in% names(edges_df)) as.character(edges_df$link_id) else paste0(g_key, "_e", seq_len(nrow(edges_df))),
      s = as.double(edges_df$s),
      u = as.double(edges_df$u),
      stringsAsFactors = FALSE
    )
    K <- nrow(edges_df)
    p_keys <- paste0("e", seq_len(K))
    group_edges[[g_key]] <<- list(cid = cid, role = role, mask = mask_log, edges = edges_df, p_keys = p_keys,
                                 dates = if (!is.null(spec_rc$dates)) as.double(spec_rc$dates) else NULL)

    Z_sel_g <- if (uses_q) {
      ref_smf_extract_design(
        panel, cl$t_prev_idx,
        if (!is.null(selection_design)) selection_design else sel_form,
        mask_log, cid
      )$Z
    } else NULL

    Z_imp_g <- if (uses_j) {
      ref_smf_extract_design(
        panel, cl$t_prev_idx,
        if (!is.null(imputation_design)) imputation_design else imp_form,
        mask_log, cid
      )$Z
    } else NULL

    g_p_cells <- character(0L)
    s_indices <- match(edges_df$s, panel$dates_num)
    u_indices <- match(edges_df$u, panel$dates_num)

    for (k in seq_len(K)) {
      pk <- p_keys[k]
      ck <- paste0(g_key, "::", pk)
      g_p_cells <- c(g_p_cells, ck)
      si <- s_indices[k]
      ui <- u_indices[k]
      gap_k <- as.integer(ui - si)
      O_pair <- (panel$O[, si] == 1L) & (panel$O[, ui] == 1L)
      dy_pair <- rep(NA_real_, n)
      obs_i <- which(mask_log & O_pair)
      dy_pair[obs_i] <- panel$Y[obs_i, ui] - panel$Y[obs_i, si]

      if (uses_q && identical(selection_model, "pair_logit")) {
        sel_cells[[ck]] <<- list(mask = mask_log, response = as.double(O_pair), design = Z_sel_g)
      }
      if (uses_j) {
        X_m <- Z_imp_g
        if (!is.null(X_m) && ncol(X_m) > 0L && identical(imputation_scale, "per_period")) {
          X_sc <- matrix(NA_real_, nrow = n, ncol = ncol(X_m), dimnames = dimnames(X_m))
          X_sc[inA, ] <- gap_k * X_m[inA, , drop = FALSE]
          X_m <- X_sc
        }
        imp_cells[[ck]] <<- list(mask = mask_log, observed = O_pair, dy = dy_pair, design = X_m)
      }
    }

    if (uses_q && identical(selection_model, "pair_logit")) {
      sel_group_cells[[g_key]] <<- g_p_cells
    }
    if (uses_j) {
      imp_group_cells[[g_key]] <<- g_p_cells
    }
    if (uses_q && identical(selection_model, "independent_period_logit")) {
      edge_date_indices <- unique(as.integer(rbind(s_indices, u_indices)))
      g_d_cells <- character(0L)
      for (di in edge_date_indices) {
        dk <- paste0("d", di)
        d_ck <- paste0(g_key, "::", dk)
        g_d_cells <- c(g_d_cells, d_ck)
        sel_cells[[d_ck]] <<- list(
          mask = mask_log,
          response = as.double(panel$O[, di] == 1L),
          design = Z_sel_g
        )
      }
      sel_group_cells[[g_key]] <<- g_d_cells
    }
  }

  for (cid in target_ids) {
    add_group(cid, "G")
    add_group(cid, "C")
  }
  for (cid in non_target_ids) {
    add_group(cid, "C")
  }

  sel_part <- if (uses_q) ref_smf_partition_blocks(selection_pooling, sel_group_cells) else NULL
  imp_part <- if (uses_j) ref_smf_partition_blocks(imputation_pooling, imp_group_cells) else NULL

  sel_fits <- if (uses_q) {
    lapply(sel_part, function(ks) ref_smf_fit_nuisance_block(sel_cells, ks, "logit", w_u, glm_control))
  } else NULL
  imp_fits <- if (uses_j) {
    lapply(imp_part, function(ks) ref_smf_fit_nuisance_block(imp_cells, ks, "wls", w_u, glm_control))
  } else NULL

  coef_if_all <- list()
  if (use_logit) {
    for (cid in target_ids) {
      coef_if_all[[paste0("trt|", cid)]] <- trt_fits[[cid]]$ell_gamma
    }
  }
  if (uses_q) {
    for (bn in names(sel_fits)) {
      coef_if_all[[paste0("sel|", bn)]] <- sel_fits[[bn]]$coef_if
    }
  }
  if (uses_j) {
    for (bn in names(imp_fits)) {
      coef_if_all[[paste0("imp|", bn)]] <- imp_fits[[bn]]$coef_if
    }
  }

  # Helper to evaluate an edge mean or covariance in group g_key, pair index k
  eval_edge <- function(g_key, k, kind = c("mean", "cov"), h_vec = NULL, h_jac = list(), Z_cov = NULL) {
    kind <- match.arg(kind)
    ge <- group_edges[[g_key]]
    mask_log <- ge$mask
    inA <- which(mask_log)
    pk <- ge$p_keys[k]
    ck <- paste0(g_key, "::", pk)
    si <- match(ge$edges$s[k], panel$dates_num)
    ui <- match(ge$edges$u[k], panel$dates_num)
    O_pair <- (panel$O[, si] == 1L) & (panel$O[, ui] == 1L)
    dy_pair <- rep(NA_real_, n)
    obs_i <- which(mask_log & O_pair)
    dy_pair[obs_i] <- panel$Y[obs_i, ui] - panel$Y[obs_i, si]

    q_vec <- NULL
    q_jac <- list()
    if (uses_q) {
      if (identical(selection_model, "pair_logit")) {
        pr_q <- ref_smf_predict_nuisance(sel_cells, sel_part, sel_fits, ck, obs_i, n, "logit")
        q_vec <- pr_q$value
        q_jac[[paste0("sel|", pr_q$block)]] <- pr_q$jac
      } else {
        ck_s <- paste0(g_key, "::d", si)
        ck_u <- paste0(g_key, "::d", ui)
        pr_s <- ref_smf_predict_nuisance(sel_cells, sel_part, sel_fits, ck_s, obs_i, n, "logit")
        pr_u <- ref_smf_predict_nuisance(sel_cells, sel_part, sel_fits, ck_u, obs_i, n, "logit")
        q_vec <- rep(NA_real_, n)
        q_vec[obs_i] <- pr_s$value[obs_i] * pr_u$value[obs_i]
        # Product rule: d(q_s * q_u) = q_u * dq_s + q_s * dq_u
        ks_ch <- paste0("sel|", pr_s$block)
        ku_ch <- paste0("sel|", pr_u$block)
        J_s <- matrix(0.0, nrow = n, ncol = ncol(pr_s$jac))
        J_u <- matrix(0.0, nrow = n, ncol = ncol(pr_u$jac))
        J_s[obs_i, ] <- pr_s$jac[obs_i, , drop = FALSE] * pr_u$value[obs_i]
        J_u[obs_i, ] <- pr_u$jac[obs_i, , drop = FALSE] * pr_s$value[obs_i]
        if (identical(ks_ch, ku_ch)) {
          q_jac[[ks_ch]] <- J_s + J_u
        } else {
          q_jac[[ks_ch]] <- J_s
          q_jac[[ku_ch]] <- J_u
        }
      }
    }

    j_vec <- NULL
    j_jac <- list()
    if (uses_j) {
      pr_j <- ref_smf_predict_nuisance(imp_cells, imp_part, imp_fits, ck, inA, n, "wls")
      j_vec <- pr_j$value
      j_jac[[paste0("imp|", pr_j$block)]] <- pr_j$jac
    }

    # Evaluate moment & complete IF via ref_recon_or_edge_moment algebra
    A_num <- as.double(mask_log)
    O_num <- as.double(O_pair)
    eps_floor <- 1e-6
    h_use <- if (kind == "mean") as.double(h_vec) else rep(1.0, n)

    phi_vec <- numeric(n)
    dphi_dq <- numeric(n)
    if (uses_q) {
      q_o <- as.double(q_vec[obs_i])
      act_o <- as.double(q_o > eps_floor)
      qf_o <- pmax(q_o, eps_floor)
      phi_vec[obs_i] <- 1.0 / qf_o
      dphi_dq[obs_i] <- -act_o / (qf_o * qf_o)
    }

    w_route <- numeric(n)
    Y_pseudo <- numeric(n)
    if (route == "baseline") {
      idx <- obs_i
      w_route[idx] <- w_u[idx] * h_use[idx]
      Y_pseudo[idx] <- dy_pair[idx]
    } else if (route == "weighting") {
      idx <- obs_i
      w_route[idx] <- w_u[idx] * h_use[idx] * phi_vec[idx]
      Y_pseudo[idx] <- dy_pair[idx]
    } else if (route == "imputation") {
      idx <- inA
      w_route[idx] <- w_u[idx] * h_use[idx]
      Y_pseudo[idx] <- j_vec[idx]
    } else {
      idx <- inA
      w_route[idx] <- w_u[idx] * h_use[idx]
      j_idx <- j_vec[idx]
      corr <- numeric(length(idx))
      obs_local <- which(O_num[idx] == 1)
      if (length(obs_local) > 0L) {
        obs_rows <- idx[obs_local]
        corr[obs_local] <- phi_vec[obs_rows] * (dy_pair[obs_rows] - j_idx[obs_local])
      }
      Y_pseudo[idx] <- j_idx + corr
    }

    S <- sum(w_route)
    M <- sum(w_route * Y_pseudo) / S

    if (kind == "mean") {
      est <- M
      Fw <- matrix(Y_pseudo - M, ncol = 1L)
      Gy <- matrix(1.0, nrow = n, ncol = 1L)
      direct <- (n_num / W_sum) * w_route * (Y_pseudo - M) / (S / W_sum)
    } else {
      Zm <- as.matrix(Z_cov)
      storage.mode(Zm) <- "double"
      p <- ncol(Zm)
      Zbar <- colSums(Zm[idx, , drop = FALSE] * w_route[idx]) / S
      Zc <- matrix(0.0, nrow = n, ncol = p, dimnames = dimnames(Zm))
      Zc[inA, ] <- sweep(Zm[inA, , drop = FALSE], 2L, Zbar, FUN = "-")
      est <- drop(crossprod(Zc[idx, , drop = FALSE], w_route[idx] * Y_pseudo[idx])) / S
      names(est) <- colnames(Zm)
      Fw <- Zc * (Y_pseudo - M) - matrix(est, nrow = n, ncol = p, byrow = TRUE)
      Gy <- Zc
      direct <- matrix(0.0, nrow = n, ncol = p, dimnames = list(uid_chr, colnames(Zm)))
      direct[idx, ] <- (n_num * w_route[idx] / S) * Fw[idx, , drop = FALSE]
    }

    loadings <- list()
    add_load <- function(lst, key, val) {
      if (is.null(lst[[key]])) lst[[key]] <- val else lst[[key]] <- lst[[key]] + val
      lst
    }
    if (kind == "mean" && length(h_jac) > 0L) {
      dw_dh <- if (route == "weighting") A_num * O_num * phi_vec else if (route == "baseline") A_num * O_num else A_num
      coef_h <- (w_u * dw_dh) * Fw
      for (k_ch in names(h_jac)) {
        Jk <- drop(crossprod(h_jac[[k_ch]][idx, , drop = FALSE], coef_h[idx, , drop = FALSE])) / S
        loadings <- add_load(loadings, k_ch, Jk)
      }
    }
    if (uses_q && length(q_jac) > 0L) {
      if (route == "weighting") {
        coef_q <- (w_u * A_num * O_num * h_use * dphi_dq) * Fw
      } else {
        resid_o <- numeric(n)
        resid_o[obs_i] <- dy_pair[obs_i] - j_vec[obs_i]
        coef_q <- (w_u * A_num * O_num * h_use * resid_o * dphi_dq) * Gy
      }
      for (k_ch in names(q_jac)) {
        val <- crossprod(q_jac[[k_ch]][obs_i, , drop = FALSE], coef_q[obs_i, , drop = FALSE]) / S
        Jk <- if (kind == "mean") drop(val[, 1L]) else t(val)
        loadings <- add_load(loadings, k_ch, Jk)
      }
    }
    if (uses_j && length(j_jac) > 0L) {
      mult_j <- if (route == "imputation") rep(1.0, n) else (1.0 - O_num * phi_vec)
      coef_j <- (w_u * A_num * h_use * mult_j) * Gy
      for (k_ch in names(j_jac)) {
        val <- crossprod(j_jac[[k_ch]][idx, , drop = FALSE], coef_j[idx, , drop = FALSE]) / S
        Jk <- if (kind == "mean") drop(val[, 1L]) else t(val)
        loadings <- add_load(loadings, k_ch, Jk)
      }
    }

    inf <- direct
    for (k_ch in names(loadings)) {
      Lk <- coef_if_all[[k_ch]]
      Jk <- loadings[[k_ch]]
      if (kind == "mean") {
        inf <- inf + drop(Lk %*% Jk)
      } else {
        inf <- inf + Lk %*% t(Jk)
      }
    }
    list(estimate = est, influence = inf)
  }

  # 4. Reconstruct stayer OR across U (if method %in% c("reg", "dr")) and build adjusted link moments on A
  mhat_list <- list()
  chi_list <- list()
  blocks_spec <- list()
  or_summary <- NULL

  if (use_or) {
    # Parse OR partition over all_ids
    if (is.list(or_blocks)) {
      b_ids <- names(or_blocks)
      if (is.null(b_ids)) b_ids <- paste0("block_", seq_along(or_blocks))
      cell_block <- stats::setNames(character(length(all_ids)), all_ids)
      for (bi in seq_along(or_blocks)) {
        mems <- or_blocks[[bi]]
        if (is.numeric(mems)) mems <- all_ids[as.integer(mems)]
        cell_block[mems] <- b_ids[bi]
      }
    } else {
      ob_chr <- as.character(or_blocks)
      if (!is.null(names(or_blocks))) ob_chr <- ob_chr[match(all_ids, names(or_blocks))]
      cell_block <- stats::setNames(ob_chr, all_ids)
      b_ids <- unique(ob_chr)
    }
    or_part_blocks <- lapply(b_ids, function(bid) all_ids[cell_block == bid])
    names(or_part_blocks) <- b_ids

    pi_C <- stats::setNames(numeric(length(all_ids)), all_ids)
    zbar_C <- list()
    S_C_mat <- list()
    varpi_C <- list()
    mC_hat <- stats::setNames(numeric(length(all_ids)), all_ids)
    Cov_hat <- list()
    psi_mC <- matrix(0.0, nrow = n, ncol = length(all_ids), dimnames = list(uid_chr, all_ids))
    psi_Cov <- list()
    psi_zbar_C <- list()
    psi_S_C <- list()

    for (cid in all_ids) {
      cl <- c_map[[cid]]
      g_key_C <- paste0(cid, "__C")
      ge_C <- group_edges[[g_key_C]]
      K_C <- nrow(ge_C$edges)
      C_num <- as.double(cl$C)
      w_C <- w_u * C_num
      S_c_sum <- sum(w_C)
      pi_C[[cid]] <- S_c_sum / W_sum

      Zm <- or_Z_by_cell[[cid]]
      p_c <- ncol(Zm)
      zn <- colnames(Zm)
      st_idx <- which(cl$C == 1L)

      if (p_c > 0L) {
        Z_st <- Zm[st_idx, , drop = FALSE]
        w_st <- w_C[st_idx]
        zb <- stats::setNames(colSums(Z_st * w_st) / S_c_sum, zn)
        Zc_st <- sweep(Z_st, 2L, zb, FUN = "-")
        Sc <- crossprod(Zc_st, Zc_st * w_st) / S_c_sum
        dimnames(Sc) <- list(zn, zn)
        pz <- matrix(0.0, nrow = n, ncol = p_c, dimnames = list(uid_chr, zn))
        pz[st_idx, ] <- (n_num * w_st / S_c_sum) * Zc_st
      } else {
        zb <- numeric(0L)
        Sc <- matrix(0.0, 0L, 0L)
        pz <- matrix(0.0, nrow = n, ncol = 0L)
      }
      zbar_C[[cid]] <- zb
      S_C_mat[[cid]] <- Sc
      psi_zbar_C[[cid]] <- pz

      inc_C <- ref_smf_incidence(ge_C$edges, cl$t_prev_num, cl$t_num, ge_C$dates)
      vp <- inc_C$varpi
      varpi_C[[cid]] <- vp

      em_C <- lapply(seq_len(K_C), function(k) eval_edge(g_key_C, k, "mean", h_vec = rep(1.0, n)))
      mC_hat[[cid]] <- sum(vp * vapply(em_C, `[[`, numeric(1L), "estimate"))
      for (k in seq_len(K_C)) {
        psi_mC[, cid] <- psi_mC[, cid] + vp[k] * em_C[[k]]$influence
      }

      if (p_c > 0L) {
        ec_C <- lapply(seq_len(K_C), function(k) eval_edge(g_key_C, k, "cov", Z_cov = Zm))
        ec_mat <- do.call(rbind, lapply(ec_C, `[[`, "estimate"))
        Cov_hat[[cid]] <- stats::setNames(drop(crossprod(vp, ec_mat)), zn)
        pCov <- matrix(0.0, nrow = n, ncol = p_c, dimnames = list(uid_chr, zn))
        for (k in seq_len(K_C)) {
          pCov <- pCov + vp[k] * ec_C[[k]]$influence
        }
        psi_Cov[[cid]] <- pCov
      } else {
        Cov_hat[[cid]] <- numeric(0L)
        psi_Cov[[cid]] <- matrix(0.0, nrow = n, ncol = 0L)
      }
    }

    b_hat <- list()
    psi_b <- list()
    alpha_hat <- stats::setNames(numeric(length(all_ids)), all_ids)
    psi_alpha <- matrix(0.0, nrow = n, ncol = length(all_ids), dimnames = list(uid_chr, all_ids))

    for (bid in b_ids) {
      cids <- or_part_blocks[[bid]]
      p_b <- ncol(or_Z_by_cell[[cids[1L]]])
      zn_b <- colnames(or_Z_by_cell[[cids[1L]]])
      if (p_b > 0L) {
        A_P <- matrix(0.0, p_b, p_b, dimnames = list(zn_b, zn_b))
        B_P <- stats::setNames(numeric(p_b), zn_b)
        for (cid in cids) {
          A_P <- A_P + pi_C[[cid]] * S_C_mat[[cid]]
          B_P <- B_P + pi_C[[cid]] * Cov_hat[[cid]]
        }
        b_P <- solve(A_P, B_P)
        names(b_P) <- zn_b
        b_hat[[bid]] <- b_P

        rhs_P <- matrix(0.0, nrow = n, ncol = p_b, dimnames = list(uid_chr, zn_b))
        for (cid in cids) {
          cl <- c_map[[cid]]
          st_idx <- which(cl$C == 1L)
          w_st <- w_u[st_idx]
          Zc_st <- sweep(or_Z_by_cell[[cid]][st_idx, , drop = FALSE], 2L, zbar_C[[cid]], FUN = "-")
          Zc_b <- drop(Zc_st %*% b_P)
          Sb <- drop(S_C_mat[[cid]] %*% b_P)
          # pi_c * (psi_Cov_c - psi_S_c b_P) + (C_c - pi_c)*(Cov_c - S_c b_P)
          # Note: pi_c * psi_S_c b_P on stayer i is (n/W_sum)*w_i*(Zc_i (Zc_i' b_P) - S_c b_P)
          term_S <- matrix(0.0, nrow = n, ncol = p_b)
          term_S[st_idx, ] <- (n_num * w_st / W_sum) * (Zc_st * Zc_b - matrix(Sb, nrow = length(st_idx), ncol = p_b, byrow = TRUE))
          C_cent <- (n_num * w_u * as.double(cl$C) / W_sum) - pi_C[[cid]]
          diff_vec <- Cov_hat[[cid]] - Sb
          rhs_P <- rhs_P + pi_C[[cid]] * psi_Cov[[cid]] - term_S + outer(C_cent, diff_vec)
        }
        psi_b[[bid]] <- rhs_P %*% solve(A_P)
      } else {
        b_hat[[bid]] <- numeric(0L)
        psi_b[[bid]] <- matrix(0.0, nrow = n, ncol = 0L)
      }

      for (cid in cids) {
        b_P <- b_hat[[bid]]
        alpha_hat[[cid]] <- mC_hat[[cid]] - if (length(b_P) > 0L) sum(zbar_C[[cid]] * b_P) else 0.0
        pa <- psi_mC[, cid]
        if (length(b_P) > 0L) {
          pa <- pa - drop(psi_zbar_C[[cid]] %*% b_P) - drop(psi_b[[bid]] %*% zbar_C[[cid]])
        }
        psi_alpha[, cid] <- pa
      }
    }

    or_summary <- list(
      alpha_hat = alpha_hat,
      b_hat = b_hat,
      mC_hat = mC_hat,
      Cov_hat = Cov_hat,
      zbar_hat = zbar_C,
      S_hat = S_C_mat,
      psi_alpha = psi_alpha,
      psi_b = psi_b
    )

    # Center links for every target cell c in A
    for (cid in target_ids) {
      cl <- c_map[[cid]]
      bid <- cell_block[[cid]]
      b_P <- b_hat[[bid]]
      alpha_c <- alpha_hat[[cid]]
      pa_c <- psi_alpha[, cid]
      pb_P <- psi_b[[bid]]
      Zm <- or_Z_by_cell[[cid]]
      p_c <- ncol(Zm)

      # Switcher role block
      g_key_G <- paste0(cid, "__G")
      ge_G <- group_edges[[g_key_G]]
      K_G <- nrow(ge_G$edges)
      g_idx <- which(cl$G == 1L)
      w_g <- w_u[g_idx]
      S_g_sum <- sum(w_g)
      if (p_c > 0L) {
        Z_g <- Zm[g_idx, , drop = FALSE]
        zG <- colSums(Z_g * w_g) / S_g_sum
        pzG <- matrix(0.0, nrow = n, ncol = p_c)
        pzG[g_idx, ] <- (n_num * w_g / S_g_sum) * sweep(Z_g, 2L, zG, FUN = "-")
      } else {
        zG <- numeric(0L)
        pzG <- matrix(0.0, nrow = n, ncol = 0L)
      }
      center_G <- alpha_c + if (p_c > 0L) sum(zG * b_P) else 0.0
      chi_center_G <- pa_c + if (p_c > 0L) (drop(pzG %*% b_P) + drop(pb_P %*% zG)) else 0.0

      step_G <- as.double(ge_G$edges$s < cl$t_num & cl$t_num <= ge_G$edges$u)
      for (k in seq_len(K_G)) {
        lid <- ge_G$edges$link_id[k]
        em_k <- eval_edge(g_key_G, k, "mean", h_vec = rep(1.0, n))
        mhat_list[[lid]] <- em_k$estimate - step_G[k] * center_G
        chi_list[[lid]] <- em_k$influence - step_G[k] * chi_center_G
      }
      blocks_spec[[length(blocks_spec) + 1L]] <- list(
        block_id = g_key_G, cell_id = cid, role = "G",
        baseline_time = cl$t_prev_num, target_time = cl$t_num,
        edges = ge_G$edges, dates = ge_G$dates
      )

      if (identical(method, "dr")) {
        g_key_C <- paste0(cid, "__C")
        ge_C <- group_edges[[g_key_C]]
        K_C <- nrow(ge_C$edges)
        tf <- trt_fits[[cid]]
        st_idx <- which(cl$C == 1L)
        w_st_r <- w_u[st_idx] * tf$r_hat[st_idx]
        S_cr <- sum(w_st_r)
        if (p_c > 0L) {
          Z_st <- Zm[st_idx, , drop = FALSE]
          zCr <- colSums(Z_st * w_st_r) / S_cr
          Zc_r <- sweep(Z_st, 2L, zCr, FUN = "-")
          pzCr <- matrix(0.0, nrow = n, ncol = p_c)
          pzCr[st_idx, ] <- (n_num * w_st_r / S_cr) * Zc_r
          # Add propensity-score derivative contribution to zCr
          J_gamma_z <- crossprod(tf$dr_dgamma[st_idx, , drop = FALSE], Zc_r * w_u[st_idx]) / S_cr
          pzCr <- pzCr + tf$ell_gamma %*% J_gamma_z
        } else {
          zCr <- numeric(0L)
          pzCr <- matrix(0.0, nrow = n, ncol = 0L)
        }
        center_C <- alpha_c + if (p_c > 0L) sum(zCr * b_P) else 0.0
        chi_center_C <- pa_c + if (p_c > 0L) (drop(pzCr %*% b_P) + drop(pb_P %*% zCr)) else 0.0

        step_C <- as.double(ge_C$edges$s < cl$t_num & cl$t_num <= ge_C$edges$u)
        h_jac_c <- list()
        h_jac_c[[paste0("trt|", cid)]] <- tf$dr_dgamma
        for (k in seq_len(K_C)) {
          lid <- ge_C$edges$link_id[k]
          em_k <- eval_edge(g_key_C, k, "mean", h_vec = tf$r_hat, h_jac = h_jac_c)
          mhat_list[[lid]] <- em_k$estimate - step_C[k] * center_C
          chi_list[[lid]] <- em_k$influence - step_C[k] * chi_center_C
        }
        blocks_spec[[length(blocks_spec) + 1L]] <- list(
          block_id = g_key_C, cell_id = cid, role = "C",
          baseline_time = cl$t_prev_num, target_time = cl$t_num,
          edges = ge_C$edges, dates = ge_C$dates
        )
      }
    }
  } else {
    # method == "ipw"
    for (cid in target_ids) {
      cl <- c_map[[cid]]
      g_key_G <- paste0(cid, "__G")
      ge_G <- group_edges[[g_key_G]]
      for (k in seq_len(nrow(ge_G$edges))) {
        lid <- ge_G$edges$link_id[k]
        em_k <- eval_edge(g_key_G, k, "mean", h_vec = rep(1.0, n))
        mhat_list[[lid]] <- em_k$estimate
        chi_list[[lid]] <- em_k$influence
      }
      blocks_spec[[length(blocks_spec) + 1L]] <- list(
        block_id = g_key_G, cell_id = cid, role = "G",
        baseline_time = cl$t_prev_num, target_time = cl$t_num,
        edges = ge_G$edges, dates = ge_G$dates
      )

      g_key_C <- paste0(cid, "__C")
      ge_C <- group_edges[[g_key_C]]
      tf <- trt_fits[[cid]]
      h_jac_c <- list()
      h_jac_c[[paste0("trt|", cid)]] <- tf$dr_dgamma
      for (k in seq_len(nrow(ge_C$edges))) {
        lid <- ge_C$edges$link_id[k]
        em_k <- eval_edge(g_key_C, k, "mean", h_vec = tf$r_hat, h_jac = h_jac_c)
        mhat_list[[lid]] <- em_k$estimate
        chi_list[[lid]] <- em_k$influence
      }
      blocks_spec[[length(blocks_spec) + 1L]] <- list(
        block_id = g_key_C, cell_id = cid, role = "C",
        baseline_time = cl$t_prev_num, target_time = cl$t_num,
        edges = ge_C$edges, dates = ge_C$dates
      )
    }
  }

  # 5. Joint static graph GMM assembly across all target cells in A
  mhat_vec <- unlist(mhat_list)
  chi_mat <- do.call(cbind, chi_list)
  rownames(chi_mat) <- uid_chr
  K_total <- length(mhat_vec)

  # Build global L matrix in target parameter parameterization
  # For ipw/dr: each cell c has C(t) + theta_c for G(t), plus non-target free dates in G and C.
  # Regardless of parameterization, the target projection Lambda_A = W L (L' W L)^{-1} E_A is invariant!
  param_cols <- list()
  target_col_map <- integer(length(target_ids))
  names(target_col_map) <- target_ids

  if (identical(method, "reg")) {
    # OR mode: only switcher blocks G; target is G(t) for each cell c
    for (b in blocks_spec) {
      cid <- b$cell_id
      inc <- ref_smf_incidence(b$edges, b$baseline_time, b$target_time, b$dates)
      t_pos <- match(as.double(b$target_time), inc$non_ground)
      # Put target column first, then any other non-ground vertices
      col_mat <- matrix(0.0, nrow = K_total, ncol = ncol(inc$D))
      rownames(col_mat) <- names(mhat_vec)
      colnames(col_mat) <- paste0(b$block_id, "::", colnames(inc$D))
      col_mat[b$edges$link_id, ] <- inc$D
      param_cols[[b$block_id]] <- col_mat
    }
    L_mat <- do.call(cbind, param_cols)
    for (cid in target_ids) {
      cl <- c_map[[cid]]
      target_col_map[[cid]] <- match(paste0(cid, "__G::", cl$t_num), colnames(L_mat))
    }
  } else {
    # IPW/DR mode: for each cell c, parameterize G(t) = C(t) + theta_c
    for (cid in target_ids) {
      cl <- c_map[[cid]]
      bG <- blocks_spec[[which(vapply(blocks_spec, function(b) b$cell_id == cid && b$role == "G", logical(1L)))]]
      bC <- blocks_spec[[which(vapply(blocks_spec, function(b) b$cell_id == cid && b$role == "C", logical(1L)))]]
      incG <- ref_smf_incidence(bG$edges, bG$baseline_time, bG$target_time, bG$dates)
      incC <- ref_smf_incidence(bC$edges, bC$baseline_time, bC$target_time, bC$dates)
      tG_pos <- match(as.double(bG$target_time), incG$non_ground)
      tC_pos <- match(as.double(bC$target_time), incC$non_ground)

      otherG <- setdiff(seq_len(ncol(incG$D)), tG_pos)
      otherC <- setdiff(seq_len(ncol(incC$D)), tC_pos)
      n_p_c <- 2L + length(otherG) + length(otherC)
      col_names_c <- c(
        paste0(cid, "::theta"),
        paste0(cid, "::C_t"),
        if (length(otherG) > 0L) paste0(cid, "::G_", incG$non_ground[otherG]) else character(0L),
        if (length(otherC) > 0L) paste0(cid, "::C_", incC$non_ground[otherC]) else character(0L)
      )
      Lc <- matrix(0.0, nrow = K_total, ncol = n_p_c, dimnames = list(names(mhat_vec), col_names_c))
      # G(t) = C(t) + theta_c => column for G(t) in incG$D loads onto both theta and C_t
      Lc[bG$edges$link_id, 1L] <- incG$D[, tG_pos]
      Lc[bG$edges$link_id, 2L] <- incG$D[, tG_pos]
      if (length(otherG) > 0L) {
        Lc[bG$edges$link_id, 2L + seq_along(otherG)] <- incG$D[, otherG, drop = FALSE]
      }
      Lc[bC$edges$link_id, 2L] <- incC$D[, tC_pos]
      if (length(otherC) > 0L) {
        Lc[bC$edges$link_id, 2L + length(otherG) + seq_along(otherC)] <- incC$D[, otherC, drop = FALSE]
      }
      param_cols[[cid]] <- Lc
    }
    L_mat <- do.call(cbind, param_cols)
    for (cid in target_ids) {
      target_col_map[[cid]] <- match(paste0(cid, "::theta"), colnames(L_mat))
    }
  }

  # Weighting matrix W
  chi_mean <- colSums(chi_mat * w_u) / W_sum
  chi_use <- if (isTRUE(centered)) sweep(chi_mat, 2L, chi_mean, FUN = "-") else chi_mat

  if (is.null(cluster)) {
    Omega_hat <- crossprod(chi_use, chi_use * (w_u^2)) / (W_sum^2) * n_num
  } else {
    cl_chr <- as.character(cluster)
    cl_levels <- unique(cl_chr)
    cl_sums <- matrix(0.0, nrow = length(cl_levels), ncol = K_total)
    for (j in seq_along(cl_levels)) {
      idx_j <- which(cl_chr == cl_levels[j])
      cl_sums[j, ] <- colSums(chi_use[idx_j, , drop = FALSE] * w_u[idx_j])
    }
    if (isTRUE(centered)) {
      cl_sums <- sweep(cl_sums, 2L, colMeans(cl_sums), FUN = "-")
    }
    Omega_hat <- crossprod(cl_sums) / (W_sum^2) * n_num
  }

  if (!is.null(fixed_W)) {
    W_mat <- fixed_W
  } else if (identical(weighting, "identity")) {
    W_mat <- diag(1.0, K_total)
  } else {
    # Optimal weighting with conditional Tikhonov ridge matching D26/D27 covariance strategy
    # Note: point estimates use unclustered Omega_1 for W_hat so point estimates are invariant to cluster!
    Omega_uncl <- crossprod(chi_use, chi_use * (w_u^2)) / (W_sum^2) * n_num
    Omega_uncl <- 0.5 * (Omega_uncl + t(Omega_uncl))
    eig_vals <- eigen(Omega_uncl, symmetric = TRUE, only.values = TRUE)$values
    d_min <- min(eig_vals)
    d_max <- max(eig_vals)
    rho <- if (d_max <= 0.0) 0.0 else max(d_min, 0.0) / d_max
    if (d_min <= 0.0 || rho < condition_trigger) {
      s_bar <- mean(diag(Omega_uncl))
      ridge <- if (s_bar > 0) ridge_relative * s_bar else ridge_relative
      Omega_reg <- Omega_uncl + diag(ridge, K_total)
    } else {
      Omega_reg <- Omega_uncl
    }
    W_mat <- solve(Omega_reg)
  }
  dimnames(W_mat) <- list(names(mhat_vec), names(mhat_vec))

  # Retain the original singular-information failure gate. The covariance-based
  # point and loading solves below do not use this inverse.
  bread_inv <- solve( t(L_mat) %*% W_mat %*% L_mat )

  E_A <- matrix(0.0, nrow = ncol(L_mat), ncol = length(target_ids),
                dimnames = list(colnames(L_mat), target_ids))
  for (cid in target_ids) {
    E_A[target_col_map[[cid]], cid] <- 1.0
  }

  if (!is.null(fixed_W) || identical(weighting, "identity")) {
    # Preserve the fixed-weight perturbation and identity branches.
    beta_hat <- drop(bread_inv %*% t(L_mat) %*% W_mat %*% mhat_vec)
    Lambda_A <- W_mat %*% L_mat %*% bread_inv %*% E_A
  } else {
    # Whiten by Omega_reg = C C', then factor A P = Q R with column pivots.
    # Neither W_mat nor bread_inv enters these point/loading solves.
    covariance_chol <- chol(Omega_reg)
    A <- forwardsolve(t(covariance_chol), L_mat)
    design_qr <- qr(A, LAPACK = TRUE)
    design_R <- qr.R(design_qr)
    pivot <- design_qr$pivot
    p <- ncol(L_mat)
    whitened_mhat <- forwardsolve(t(covariance_chol), as.double(mhat_vec))
    qty <- qr.qty(design_qr, whitened_mhat)
    beta_hat <- numeric(p)
    beta_hat[pivot] <- as.double(backsolve(design_R, qty[seq_len(p)]))

    # Lambda_A = C^{-T} Q R^{-T} P' E_A, without normal equations.
    target_coords <- forwardsolve(t(design_R), E_A[pivot, , drop = FALSE])
    Lambda_A <- backsolve(covariance_chol, qr.qy(
      design_qr, rbind(target_coords, matrix(0.0, K_total - p, length(target_ids)))
    ))
    rownames(Lambda_A) <- names(mhat_vec)
  }
  names(beta_hat) <- colnames(L_mat)
  colnames(Lambda_A) <- target_ids

  theta_hat <- stats::setNames(beta_hat[target_col_map], target_ids)
  sigmas <- vapply(c_map[target_ids], function(cl) as.double(cl$sigma), numeric(1L))
  names(sigmas) <- target_ids
  tau_hat <- sigmas * theta_hat

  # First-order IF of target contrasts: phi_theta = chi_use %*% Lambda_A
  phi_theta <- chi_use %*% Lambda_A
  dimnames(phi_theta) <- list(uid_chr, target_ids)
  phi_tau <- sweep(phi_theta, 2L, sigmas, FUN = "*")

  if (is.null(cluster)) {
    phi_c <- if (isTRUE(centered)) sweep(phi_theta, 2L, colMeans(phi_theta), FUN = "-") else phi_theta
    cov_theta <- crossprod(phi_c) / (n_num^2)
  } else {
    cl_chr <- as.character(cluster)
    cl_levels <- unique(cl_chr)
    cl_phi <- matrix(0.0, nrow = length(cl_levels), ncol = length(target_ids),
                     dimnames = list(cl_levels, target_ids))
    for (j in seq_along(cl_levels)) {
      idx_j <- which(cl_chr == cl_levels[j])
      cl_phi[j, ] <- colSums(phi_theta[idx_j, , drop = FALSE])
    }
    if (isTRUE(centered)) {
      cl_phi <- sweep(cl_phi, 2L, colMeans(cl_phi), FUN = "-")
    }
    cov_theta <- crossprod(cl_phi) / (n_num^2)
  }
  vhat_theta <- cov_theta * n_num
  se_theta <- sqrt(pmax(0.0, diag(cov_theta)))
  cov_tau <- (sigmas %o% sigmas) * cov_theta
  vhat_tau <- cov_tau * n_num

  agg_res <- NULL
  if (!is.null(aggregate_cells)) {
    agg_res <- ref_static_missing_aggregate(
      c_map = c_map,
      theta_hat = theta_hat,
      phi_theta = phi_theta,
      cells = aggregate_cells,
      signs = aggregate_signs,
      unit_weights = w_u,
      cluster = cluster
    )
  }

  list(
    theta_hat = theta_hat,
    tau_hat = tau_hat,
    se_theta = se_theta,
    cov_theta = cov_theta,
    vhat_theta = vhat_theta,
    cov_tau = cov_tau,
    vhat_tau = vhat_tau,
    phi_theta = phi_theta,
    phi_tau = phi_tau,
    mhat = mhat_vec,
    chi = chi_use,
    L = L_mat,
    W = W_mat,
    Lambda = Lambda_A,
    beta_hat = beta_hat,
    overid_residual = drop(mhat_vec - L_mat %*% beta_hat),
    trt_fits = trt_fits,
    or_summary = or_summary,
    aggregate = agg_res
  )
}

ref_static_missing_aggregate <- function(c_map,
                                         theta_hat,
                                         phi_theta,
                                         cells,
                                         signs,
                                         unit_weights = NULL,
                                         cluster = NULL) {
  n <- nrow(phi_theta)
  n_num <- as.double(n)
  w_u <- if (is.null(unit_weights)) rep(1.0, n) else as.double(unit_weights)
  W_sum <- sum(w_u)

  agg_ids <- as.character(cells)
  s_vec <- if (is.null(signs)) rep(1.0, length(agg_ids)) else as.double(signs)
  names(s_vec) <- agg_ids

  G_mat <- vapply(agg_ids, function(cid) as.double(c_map[[cid]]$G), numeric(n))
  if (n == 1L) G_mat <- matrix(G_mat, nrow = 1L, dimnames = list(rownames(phi_theta), agg_ids))
  colnames(G_mat) <- agg_ids

  pi_c <- colSums(G_mat * w_u) / W_sum
  Pi_hat <- sum(pi_c)
  omega_hat <- pi_c / Pi_hat
  names(omega_hat) <- agg_ids

  th_sub <- theta_hat[agg_ids]
  signed_theta <- s_vec * th_sub
  delta_hat <- sum(omega_hat * signed_theta)

  phi_sub <- phi_theta[, agg_ids, drop = FALSE]
  cell_if_term <- drop(phi_sub %*% (omega_hat * s_vec))
  G_cent <- sweep(G_mat * (n_num * w_u / W_sum), 2L, pi_c, FUN = "-")
  share_weights <- (signed_theta - delta_hat) / Pi_hat
  share_if_term <- drop(G_cent %*% share_weights)
  phi_agg <- cell_if_term + share_if_term
  names(phi_agg) <- rownames(phi_theta)
  names(cell_if_term) <- rownames(phi_theta)
  names(share_if_term) <- rownames(phi_theta)

  if (is.null(cluster)) {
    var_agg <- sum(phi_agg^2) / (n_num^2)
  } else {
    cl_chr <- as.character(cluster)
    cl_levels <- unique(cl_chr)
    cl_sums <- vapply(cl_levels, function(lv) sum(phi_agg[cl_chr == lv]), numeric(1L))
    var_agg <- sum(cl_sums^2) / (n_num^2)
  }

  list(
    delta_hat = delta_hat,
    omega_hat = omega_hat,
    pi_hat = pi_c,
    Pi_hat = Pi_hat,
    phi = phi_agg,
    phi_components = list(
      cell_if_term = cell_if_term,
      share_if_term = share_if_term
    ),
    variance = var_agg,
    vhat = var_agg * n_num,
    se = sqrt(max(0.0, var_agg))
  )
}

# Central finite-difference original-unit perturbation evaluator.
# When hold_W_fixed = TRUE, holds the final GMM weighting matrix W fixed at the
# unperturbed estimate while refitting all nuisances, normalizations, reconstructed OR,
# centering, and linear GMM projection (exact finite-sample witness of the first-order IF).
ref_static_missing_fd_if <- function(units = NULL,
                                     eps = 1e-6,
                                     hold_W_fixed = TRUE,
                                     ...) {
  base_fit <- ref_static_missing_fit(...)
  n <- nrow(base_fit$phi_theta)
  if (is.null(units)) units <- seq_len(n)
  W_fix <- if (isTRUE(hold_W_fixed)) base_fit$W else NULL

  q_cells <- length(base_fit$theta_hat)
  fd_theta <- matrix(0.0, nrow = length(units), ncol = q_cells,
                     dimnames = list(rownames(base_fit$phi_theta)[units], names(base_fit$theta_hat)))
  fd_agg <- if (!is.null(base_fit$aggregate)) numeric(length(units)) else NULL
  if (!is.null(fd_agg)) names(fd_agg) <- rownames(base_fit$phi_theta)[units]

  for (idx in seq_along(units)) {
    i <- units[idx]
    wp <- wm <- rep(1.0, n)
    wp[i] <- 1.0 + eps
    wm[i] <- 1.0 - eps
    fit_p <- ref_static_missing_fit(..., unit_weights = wp, fixed_W = W_fix)
    fit_m <- ref_static_missing_fit(..., unit_weights = wm, fixed_W = W_fix)
    fd_theta[idx, ] <- n * (fit_p$theta_hat - fit_m$theta_hat) / (2.0 * eps)
    if (!is.null(fd_agg)) {
      fd_agg[idx] <- n * (fit_p$aggregate$delta_hat - fit_m$aggregate$delta_hat) / (2.0 * eps)
    }
  }

  list(
    base_fit = base_fit,
    fd_theta = fd_theta,
    fd_agg = fd_agg
  )
}
