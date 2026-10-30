# Private tiny-fixture diagnostics; intentionally outside the package R namespace.
# The test environment inherits the actual locked installed cdid namespace.
# ------------------------------------------------------------------------------
# 6. Pure-Helper Table Entry Point (`fit_imputation_table`)
# ------------------------------------------------------------------------------
# Exercises the exact same `imputation_build_catalog`, `imputation_fit_core`,
# and `imputation_evaluate_if_chunk_core` kernels on table inputs.

fit_imputation_table <- function(X, Y = NULL, obs, g, w,
                                 calendar, links, outcome_recipe,
                                 treatment_columns = NULL,
                                 est_method = "dr",
                                 imputation_pooling = "separate",
                                 pair_outcomes = NULL,
                                 unit_id = NULL,
                                 cluster_id = NULL,
                                 max_chunk_rows = NULL,
                                 glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
                                 rank_tol = 1e-10,
                                 selection_model = NULL,
                                 selection_pooling = NULL,
                                 selection_xformla = NULL) {
  sel_rec <- imputation_inactive_selection_record()
  est_method <- imputation_method_choice(est_method = est_method)
  imputation_pooling <- imputation_pooling_choice(imputation_pooling)
  rank_tol <- or_baseline_rank_tol(rank_tol)

  N <- nrow(X)
  if (is.null(unit_id)) unit_id <- as.double(seq_len(N)) else unit_id <- as.double(unit_id)
  if (is.null(cluster_id)) cluster_id <- unit_id else cluster_id <- as.double(cluster_id)
  if (is.null(max_chunk_rows)) max_chunk_rows <- N
  ctrl <- list(
    max_chunk_rows = as.integer(max(1L, min(N, max_chunk_rows))),
    max_retained_entries = 1e9,
    max_matrix_entries = 1e9,
    max_work_ops = 1e12,
    max_output_entries = 1e9
  )

  g_num <- as.double(g)
  w_num <- as.double(w)
  A_ind <- as.integer(rowSums(obs != 0) >= 2L)
  n_A <- sum(A_ind > 0L)
  a_index <- integer(N)
  a_index[A_ind > 0L] <- seq_len(n_A)

  cohorts <- sort(unique(g_num), method = "radix")
  ci <- match(g_num, cohorts)
  ord <- order(ci, method = "radix")
  counts <- tabulate(ci, length(cohorts))
  ends <- cumsum(as.double(counts))
  starts <- ends - counts + 1
  mass <- repeat_mass <- numeric(length(cohorts))
  for (j in seq_along(cohorts)) {
    rr <- ord[seq.int(starts[j], ends[j])]
    mass[j] <- sum(w_num[rr])
    repeat_mass[j] <- sum(w_num[rr[a_index[rr] > 0L]])
  }

  st <- list(
    id = unit_id,
    Z = X,
    g = g_num,
    w = w_num,
    U = rep(TRUE, N),
    rho = rep(1.0, N),
    universe = "F",
    mass_U = mass / N,
    a_index = a_index,
    cohort_index = ci,
    cohort_order = ord,
    cohort_start = starts,
    cohort_end = ends,
    cohorts = cohorts,
    mass_F = mass / N,
    a = ifelse(mass > 0, repeat_mass / mass, 0),
    eta = repeat_mass / sum(repeat_mass),
    tau = mass / sum(mass),
    mean_w = sum(mass) / N,
    mean_wA = sum(repeat_mass) / N,
    N = N,
    n_A = n_A,
    lambda = n_A / N
  )

  f_cols <- colnames(X)
  or_chk <- dr_check_outcome_recipe(outcome_recipe, f_cols, "baseline")
  recipe <- or_chk$recipe
  or_frame_index <- or_chk$frame_index
  trt_recipe <- if (est_method %in% c("ipw", "dr")) {
    dr_check_treatment_recipe(treatment_columns, f_cols, "baseline")
  } else {
    structure(list(columns = character(0L), index = integer(0L), p = 0L), class = "cdid_dr_treatment_recipe")
  }

  records <- vector("list", length(links))
  for (ell in seq_along(links)) {
    lk <- links[[ell]]
    s_cal <- match(as.double(lk$s), calendar)
    t_cal <- match(as.double(lk$t), calendar)
    pair <- list(
      s = as.double(lk$s),
      t = as.double(lk$t),
      s_calendar_index = as.integer(s_cal),
      t_calendar_index = as.integer(t_cal)
    )
    records[[ell]] <- list(
      pair = pair,
      target = match(as.double(lk$g), cohorts),
      controls = match(sort(as.double(lk$K), method = "radix"), cohorts),
      link_id = lk$link_id,
      source_functional_key = paste0("cmpcat:", lk$link_id),
      functional_key = paste0("cmpcat:IMP-", toupper(est_method), ":", imputation_pooling, ":", lk$link_id)
    )
  }

  cat <- imputation_build_catalog(records, cohorts, est_method, imputation_pooling)

  pair_reader <- function(rows, pair, read_outcomes = TRUE) {
    s_col <- pair$s_calendar_index
    t_col <- pair$t_calendar_index
    S_vec <- as.double(obs[rows, s_col] != 0 & obs[rows, t_col] != 0)
    out <- list(S = S_vec)
    if (read_outcomes) {
      dy_vec <- numeric(length(rows))
      obs_idx <- which(S_vec == 1)
      if (length(obs_idx)) {
        r_obs <- rows[obs_idx]
        dy_vec[obs_idx] <- dr_table_pair_dy(Y, pair_outcomes, pair, r_obs)
      }
      out$dy <- dy_vec
    }
    out
  }

  core_fit <- imputation_fit_core(
    st = st,
    X_F = X,
    pair_reader = pair_reader,
    cat = cat,
    recipe = recipe,
    or_frame_index = or_frame_index,
    trt_recipe = trt_recipe,
    est_method = est_method,
    imputation_pooling = imputation_pooling,
    glm_control = glm_control,
    rank_tol = rank_tol,
    ctrl = ctrl
  )

  # Evaluate complete original-unit IF in chunks of `ctrl$max_chunk_rows`
  L_links <- length(links)
  all_links <- seq_len(L_links)
  chunks <- list()
  for (start in seq.int(1L, N, by = ctrl$max_chunk_rows)) {
    rr <- seq.int(start, min(N, start + ctrl$max_chunk_rows - 1L))
    chunks[[length(chunks) + 1L]] <- imputation_evaluate_if_chunk_core(
      st = st,
      X_F = X,
      pair_reader = pair_reader,
      core_fit = core_fit,
      recipe = recipe,
      or_frame_index = or_frame_index,
      trt_recipe = trt_recipe,
      est_method = est_method,
      rows = rr,
      links = all_links,
      scale = "frame",
      details = TRUE
    )
  }

  # Reassemble full-N IF matrices from chunks
  phi <- do.call(rbind, lapply(chunks, `[[`, "phi"))
  psi_ipw <- do.call(rbind, lapply(chunks, function(ch) ch$details$psi_ipw))
  psi_reg <- do.call(rbind, lapply(chunks, function(ch) ch$details$psi_reg))
  psi_dr <- do.call(rbind, lapply(chunks, function(ch) ch$details$psi_dr))
  xi_a <- do.call(rbind, lapply(chunks, function(ch) ch$details$a))
  xi_tau <- do.call(rbind, lapply(chunks, function(ch) ch$details$tau))
  xi_eta <- do.call(rbind, lapply(chunks, function(ch) ch$details$eta))
  xi_alpha_imp <- do.call(rbind, lapply(chunks, function(ch) ch$details$imputation_alpha))
  xi_beta_imp <- lapply(seq_along(core_fit$groups), function(gid) {
    do.call(rbind, lapply(chunks, function(ch) ch$details$imputation_beta[[gid]]))
  })
  or_theta_if <- if (!is.null(core_fit$shared_or)) {
    do.call(rbind, lapply(chunks, function(ch) ch$details$or_if$theta))
  } else {
    NULL
  }
  trt_gamma_if <- if (!is.null(core_fit$treatments)) {
    lapply(seq_along(core_fit$treatments), function(tid) {
      do.call(rbind, lapply(chunks, function(ch) ch$details$trt_if[[tid]]$gamma))
    })
  } else {
    NULL
  }

  list(
    selection_est_method = "imputation",
    est_method = est_method,
    imputation_pooling = imputation_pooling,
    selection_record = sel_rec,
    state = st,
    catalog = cat,
    core = core_fit,
    estimates = vapply(core_fit$links, `[[`, numeric(1L), "mu"),
    phi = phi,
    phi_core = st$lambda * phi,
    psi_ipw = psi_ipw,
    psi_reg = psi_reg,
    psi_dr = psi_dr,
    xi_a = xi_a,
    xi_tau = xi_tau,
    xi_eta = xi_eta,
    xi_alpha_imp = xi_alpha_imp,
    xi_beta_imp = xi_beta_imp,
    or_theta_if = or_theta_if,
    trt_gamma_if = trt_gamma_if,
    chunks = chunks
  )
}

# ------------------------------------------------------------------------------
# 7. Complete Coordinate, Component, Prediction, and Joint Covariance Extractor
# ------------------------------------------------------------------------------
# Works on both a sealed `cdid_imputation_backend` (replaying bounded chunks via
# `evaluate_imputation_if_chunk(..., details = TRUE)`) and a `fit_imputation_table`
# result, exposing all nuisance scalars, complete N-row IFs, 8-channel component
# decompositions, cell predictions, and joint unit/cluster/core covariances.

imputation_reference_summary <- function(object, cluster_id = NULL, chunk_rows = NULL, fixed_test_functions = NULL) {
  is_sealed <- inherits(object, "cdid_imputation_backend")
  if (is_sealed) {
    validate_imputation_backend(object)
    reader <- get("cp04_reader", object$reader_ref, inherits = FALSE)
    st <- reader$state
    X_F <- imputation_backend_X_F(object, reader)
    core <- object
    recipe <- object$outcome_recipe
    B <- if (is.null(chunk_rows)) object$control$max_chunk_rows else min(st$N, as.integer(chunk_rows))
    L_links <- length(object$links)
    all_links <- seq_len(L_links)
    chunks <- list()
    for (start in seq.int(1L, st$N, by = B)) {
      rr <- seq.int(start, min(st$N, start + B - 1L))
      chunks[[length(chunks) + 1L]] <- evaluate_imputation_if_chunk(
        object, rr, all_links, scale = "frame", details = TRUE
      )
    }
  } else if (is.list(object) && !is.null(object$core) && !is.null(object$state)) {
    st <- object$state
    X_F <- st$Z
    core <- object$core
    chunks <- object$chunks
    recipe <- NULL
  } else {
    imputation_abort("Expected a `cdid_imputation_backend` or `fit_imputation_table` result.")
  }

  N <- st$N
  G <- length(st$cohorts)
  coh_ord <- order(imputation_cohort_sort_key(st$cohorts))
  coh_labels <- imputation_cohort_label(st$cohorts)
  w <- st$w
  A <- as.double(st$a_index > 0L)
  ci <- st$cohort_index

  # Reassemble chunked detail arrays
  xi_a <- do.call(rbind, lapply(chunks, function(ch) ch$details$a))
  xi_tau <- do.call(rbind, lapply(chunks, function(ch) ch$details$tau))
  xi_eta <- do.call(rbind, lapply(chunks, function(ch) ch$details$eta))
  xi_alpha_imp <- do.call(rbind, lapply(chunks, function(ch) ch$details$imputation_alpha))
  xi_beta_imp <- lapply(seq_along(core$groups), function(gid) {
    do.call(rbind, lapply(chunks, function(ch) ch$details$imputation_beta[[gid]]))
  })

  nuis_vals <- numeric(0L)
  if_cols <- list()
  add_coord <- function(nm, val, psi_vec) {
    nuis_vals[[nm]] <<- as.double(val)
    if_cols[[nm]] <<- as.double(psi_vec)
  }

  # 1. Repeat rates a, shares tau, eta, and degree-one counts
  for (k in coh_ord) {
    clab <- coh_labels[k]
    add_coord(paste0("a.k", clab), st$a[k], xi_a[, k])
  }
  for (k in coh_ord) {
    clab <- coh_labels[k]
    add_coord(paste0("tau.k", clab), st$tau[k], xi_tau[, k])
  }
  for (k in coh_ord) {
    clab <- coh_labels[k]
    add_coord(paste0("eta.k", clab), st$eta[k], xi_eta[, k])
  }
  for (k in coh_ord) {
    clab <- coh_labels[k]
    I_k <- as.double(ci == k)
    cF_k <- st$mass_F[k]
    cA_k <- sum(w * I_k * A) / N
    add_coord(paste0("countF.k", clab), cF_k, w * I_k)
    add_coord(paste0("countA.k", clab), cA_k, w * I_k * A)
  }
  add_coord("countF.total", st$mean_w, w)
  add_coord("countA.total", st$mean_wA, w * A)

  # 2. Imputation group alphas and betas
  has_icpt <- length(core$groups[[1L]]$alpha) > 0L
  for (gid in seq_along(core$groups)) {
    gf <- core$groups[[gid]]
    g_prefix <- paste0("imputation.group", gf$group_index0)
    if (has_icpt) {
      for (loc in seq_along(gf$cell_ids)) {
        cid <- gf$cell_ids[loc]
        cf <- core$cells[[cid]]
        add_coord(
          paste0(g_prefix, ".alpha.", cf$name),
          cf$alpha,
          xi_alpha_imp[, cid]
        )
      }
    }
    if (gf$p > 0L) {
      for (j in seq_len(gf$p)) {
        col_nm <- gf$kept_columns[j]
        add_coord(
          paste0(g_prefix, ".beta.", col_nm),
          gf$beta[j],
          xi_beta_imp[[gid]][, j]
        )
      }
    }
  }

  # 3. Shared counterfactual OR parameters (when fitted)
  if (!is.null(core$shared_or)) {
    sh <- core$shared_or
    or_theta_if <- do.call(rbind, lapply(chunks, function(ch) ch$details$or_if$theta))
    for (j in seq_along(sh$theta)) {
      nm <- names(sh$theta)[j]
      add_coord(nm, sh$theta[j], or_theta_if[, j])
    }
  }

  # 4. Treatment propensity parameters per link (when fitted)
  if (!is.null(core$treatments)) {
    trt_if_all <- lapply(seq_along(core$treatments), function(tid) {
      do.call(rbind, lapply(chunks, function(ch) ch$details$trt_if[[tid]]$gamma))
    })
    for (ell in seq_along(core$links)) {
      lf <- core$links[[ell]]
      tid <- lf$trt_id
      tf <- core$treatments[[tid]]
      tif <- trt_if_all[[tid]]
      lid_short <- sub("^link_", "l", lf$link_id)
      for (j in seq_along(tf$coefficients)) {
        add_coord(
          paste0("treatment.", lid_short, ".gamma", j - 1L),
          tf$coefficients[j],
          tif[, j]
        )
      }
    }
  }

  # 5. Roster moments: tiny fixture diagnostics only. These raw degree-one
  # rows intentionally retain the case-weight normalization convention.
  h_cols <- which(colnames(X_F) != "(Intercept)")
  for (cid in seq_along(core$cells)) {
    cf <- core$cells[[cid]]
    gf <- core$groups[[cf$group_id]]
    k <- cf$cohort_index
    I_k <- as.double(ci == k)
    a_k <- st$a[k]
    d <- if (gf$p) cf$L * X_F[, gf$kept_cols_F, drop = FALSE] else matrix(0, N, 0L)
    jhat <- cf$alpha + if (gf$p) as.vector(d %*% gf$beta) else rep(0, N)
    if (is.null(fixed_test_functions)) {
      hs <- c(list(one = rep(1, N)), lapply(h_cols, function(j) X_F[, j]))
      names(hs) <- c("one", colnames(X_F)[h_cols])
    } else {
      hs <- fixed_test_functions
      stopifnot(is.list(hs), !is.null(names(hs)), !anyDuplicated(names(hs)),
        all(vapply(hs, function(h) is.numeric(h) && length(h) == N && all(is.finite(h)), TRUE)))
    }
    for (hn in names(hs)) {
      hh <- hs[[hn]]
      raw <- w * I_k * hh
      M <- sum(raw * a_k) / N
      Y <- sum(raw * a_k * jhat) / N
      psi_M <- raw * a_k + sum(raw) / N * xi_a[, k]
      psi_Y <- raw * a_k * jhat + sum(raw * jhat) / N * xi_a[, k] +
        (if (length(gf$alpha)) M * xi_alpha_imp[, cid] else 0) +
        (if (gf$p) as.vector(xi_beta_imp[[cf$group_id]] %*% (colSums(raw * a_k * d) / N)) else 0)
      add_coord(paste0("M.", cf$name, ".", hn), M, psi_M)
      add_coord(paste0("Y.", cf$name, ".", hn), Y, psi_Y)
    }
  }

  # 6. Link ratios, points, and 8-channel component decompositions
  components <- list()
  for (ell in seq_along(core$links)) {
    lf <- core$links[[ell]]
    lid_short <- sub("^link_", "l", lf$link_id)
    get_l_psi <- function(field) {
      unlist(lapply(chunks, function(ch) ch$details$links[[ell]][[field]]), use.names = FALSE)
    }
    add_coord(paste0("denplus.", lid_short), lf$d_plus, get_l_psi("denplus"))
    add_coord(paste0("numT.", lid_short), lf$num_T, get_l_psi("numT"))
    add_coord(paste0("T.", lid_short), lf$T_gl, get_l_psi("T"))
    if (!is.null(core$shared_or)) {
      add_coord(paste0("numBg.", lid_short), lf$num_Bg, get_l_psi("numBg"))
      add_coord(paste0("Bg.", lid_short), lf$B_gl, get_l_psi("Bg"))
    }
    if (!is.null(core$treatments)) {
      add_coord(paste0("denminus.", lid_short), lf$d_minus, get_l_psi("denminus"))
      add_coord(paste0("numC.", lid_short), lf$num_C, get_l_psi("numC"))
      add_coord(paste0("C.", lid_short), lf$C_ol, get_l_psi("C"))
    }
    if (!is.null(core$shared_or) && !is.null(core$treatments)) {
      add_coord(paste0("numBo.", lid_short), lf$num_Bo, get_l_psi("numBo"))
      add_coord(paste0("Bo.", lid_short), lf$B_ol, get_l_psi("Bo"))
    }
    if (!is.null(core$treatments)) {
      add_coord(paste0("ipw.", lid_short), lf$mu_ipw, get_l_psi("ipw"))
      components[[paste0("ipw.", lid_short)]] <- do.call(
        rbind, lapply(chunks, function(ch) ch$details$links[[ell]]$components_ipw)
      )
    }
    if (!is.null(core$shared_or)) {
      add_coord(paste0("reg.", lid_short), lf$mu_reg, get_l_psi("reg"))
      components[[paste0("reg.", lid_short)]] <- do.call(
        rbind, lapply(chunks, function(ch) ch$details$links[[ell]]$components_reg)
      )
    }
    if (!is.null(core$shared_or) && !is.null(core$treatments)) {
      add_coord(paste0("dr.", lid_short), lf$mu_dr, get_l_psi("dr"))
      components[[paste0("dr.", lid_short)]] <- do.call(
        rbind, lapply(chunks, function(ch) ch$details$links[[ell]]$components_dr)
      )
    }
  }

  unit_if_mat <- do.call(cbind, if_cols)

  # 7. Cell predictions j_hat_{k,e}(X_i) for all i in 1..N
  pred_mat <- matrix(0.0, N, length(core$cells))
  colnames(pred_mat) <- vapply(core$cells, `[[`, character(1L), "name")
  for (cid in seq_along(core$cells)) {
    cf <- core$cells[[cid]]
    gf <- core$groups[[cf$group_id]]
    Z_imp <- if (gf$p > 0L) X_F[, gf$kept_cols_F, drop = FALSE] else matrix(0.0, N, 0L)
    pred_mat[, cid] <- if (gf$p > 0L) {
      cf$alpha + as.vector((cf$L * Z_imp) %*% gf$beta)
    } else {
      rep(cf$alpha, N)
    }
  }

  # 8. Joint covariances when all 3 methods are present
  covariances <- NULL
  l_shorts <- vapply(core$links, function(lf) sub("^link_", "l", lf$link_id), character(1L))
  joint_coords <- c(
    paste0("ipw.", l_shorts),
    paste0("reg.", l_shorts),
    paste0("dr.", l_shorts),
    paste0("tau.k", coh_labels[coh_ord]),
    paste0("eta.k", coh_labels[coh_ord])
  )
  if (all(joint_coords %in% colnames(unit_if_mat))) {
    Psi <- unit_if_mat[, joint_coords, drop = FALSE]
    Psi_c <- sweep(Psi, 2L, colMeans(Psi), FUN = "-")
    V_unit <- crossprod(Psi_c) / (N^2)
    if (is.null(cluster_id)) cluster_id <- st$id
    u_cl <- sort(unique(cluster_id))
    C_mat <- matrix(0.0, length(u_cl), ncol(Psi_c), dimnames = list(NULL, joint_coords))
    for (g_i in seq_along(u_cl)) {
      ix <- which(cluster_id == u_cl[g_i])
      C_mat[g_i, ] <- colSums(Psi_c[ix, , drop = FALSE])
    }
    V_cluster <- crossprod(C_mat) / (N^2)
    Omega_cluster <- crossprod(C_mat) / N
    # Core conversion: Psi_core = (n_A / N) * Psi_c, normalized by n_A^2
    Psi_core <- st$lambda * Psi_c
    C_core <- st$lambda * C_mat
    V_core_cluster <- crossprod(C_core) / (st$n_A^2)
    covariances <- list(
      unit_coefficient = V_unit,
      cluster_coefficient = V_cluster,
      cluster_sqrtN = Omega_cluster,
      core_cluster_coefficient = V_core_cluster
    )
  }

  list(
    nuisance = nuis_vals,
    unit_if = unit_if_mat,
    components = components,
    predictions = pred_mat,
    covariances = covariances
  )
}
