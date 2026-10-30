# Deterministic installed-namespace and pure-helper tests for the D48 Step 2
# Section A probability-model backend (R/chained_selection_models.R).
#
# When cdid is installed, every binding is resolved and verified from the locked
# cdid namespace via selection_ns(), exercising full sealed panel/frame/reader/link
# owners as well as pure-helper entry points. When sourced in pure-helper mode
# (with `selection_models_helpers` supplied in the calling environment), all
# pure-helper and reference tests run directly without fabricating a fake namespace
# or fake owner.

selection_has_installed_cdid <- function() {
  if (exists("selection_models_helpers", inherits = TRUE)) return(FALSE)
  isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE)
}

selection_ns <- function(name) {
  if (exists("selection_models_helpers", inherits = TRUE)) {
    env <- get("selection_models_helpers", inherits = TRUE)
    stopifnot(is.environment(env), exists(name, envir = env, inherits = FALSE))
    fn <- get(name, envir = env, inherits = FALSE)
    stopifnot(is.function(fn))
    return(fn)
  }
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}

selection_glm_control <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
selection_glm_control_tight <- list(epsilon = 1e-15, maxit = 60L, trace = FALSE)

# Deterministic 32-atom generator reproducing the exact frozen reference laws
# (validation/reference/selection_models/independent_reference/reference.py,
# pair_pooling_reference/reference.py, and separate_date_reference.py) with zero RNG.
selection_build_32atom_laws <- function() {
  expit <- function(v) 1 / (1 + exp(-v))
  grid <- expand.grid(
    R3 = c(0, 1), R2 = c(0, 1), R1 = c(0, 1),
    x = c(-1, 1), cohort = c(0L, 1L)
  )
  K <- as.integer(grid$cohort)
  X <- as.double(grid$x)
  R <- cbind(R1 = as.double(grid$R1), R2 = as.double(grid$R2), R3 = as.double(grid$R3))
  N <- length(K)
  A <- as.double(rowSums(R) >= 2)
  B <- as.double(rowSums(R) >= 1)
  S <- cbind(S12 = R[, 1] * R[, 2], S13 = R[, 1] * R[, 3], S23 = R[, 2] * R[, 3])
  cell_mass <- matrix(c(0.21, 0.29, 0.31, 0.19), nrow = 2L, byrow = TRUE)
  base_mass <- cell_mass[cbind(K + 1L, ifelse(X > 0, 2L, 1L))]

  make_date_law <- function(theta, mode) {
    d <- length(theta)
    p <- matrix(0, N, 3L)
    for (i in seq_len(N)) {
      for (t in 1:3) {
        b <- 3L * K[i] + t
        slope <- if (mode == "date_only") 0 else if (mode == "within_cohort") theta[6L + K[i] + 1L] * X[i] else if (mode == "across_cohorts") theta[7L] * X[i] else theta[6L + b] * X[i]
        p[i, t] <- expit(theta[b] + slope)
      }
    }
    atom_prob <- apply(p^R * (1 - p)^(1 - R), 1L, prod)
    w <- N * base_mass * atom_prob
    list(theta = theta, w = w, p = p)
  }

  date_laws <- list(
    within_cohort = make_date_law(c(-0.8, 0.2, 0.7, 0.4, -0.5, 1.1, 0.55, -0.35), "within_cohort"),
    across_cohorts = make_date_law(c(-0.8, 0.2, 0.7, 0.4, -0.5, 1.1, 0.4), "across_cohorts"),
    date_only = make_date_law(c(-0.8, 0.2, 0.7, 0.4, -0.5, 1.1), "date_only"),
    separate_date_F = make_date_law(c(-0.8, 0.2, 0.7, 0.4, -0.5, 1.1, 0.55, -0.2, 0.35, -0.35, 0.45, 0.15), "separate")
  )

  # Direct-pair 32-atom law
  pop_pair <- expit(c(-2.2, -1.9)[K + 1L] + 0.3 * X)
  hist_mass <- numeric(N)
  sing_w <- c(0.2, 0.15, 0.15)
  for (i in seq_len(N)) {
    ones <- sum(R[i, ])
    p_i <- pop_pair[i]
    hist_mass[i] <- if (ones >= 2) 0.5 * p_i else if (ones == 0) (1 - 2.5 * p_i) else p_i * sing_w[which.max(R[i, ])]
  }
  W_pair <- N * base_mass * hist_mass
  masks <- list(F = rep(TRUE, N), B = (B == 1), A = (A == 1))
  rho_list <- list()
  for (u in c("F", "B", "A")) {
    m_u <- as.double(masks[[u]])
    up <- if (u == "F") rep(1, N) else (if (u == "B") 2.5 else 2.0) * pop_pair
    rho_u <- rep(1, N)
    for (k in 0:1) {
      ik <- (K == k)
      ubar <- sum(W_pair * ik * m_u) / sum(W_pair * ik)
      rho_u[ik] <- up[ik] / ubar
    }
    rho_list[[u]] <- rho_u
  }

  pair_specs <- list(
    list(cohort = 1L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 1L, pair = list(s = 1, t = 3, s_calendar_index = 1L, t_calendar_index = 3L)),
    list(cohort = 1L, pair = list(s = 2, t = 3, s_calendar_index = 2L, t_calendar_index = 3L)),
    list(cohort = 2L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 2L, pair = list(s = 1, t = 3, s_calendar_index = 1L, t_calendar_index = 3L)),
    list(cohort = 2L, pair = list(s = 2, t = 3, s_calendar_index = 2L, t_calendar_index = 3L))
  )
  probes_df <- data.frame(cohort = c(0L, 0L, 1L, 1L), x = c(-1, 1, -1, 1))

  list(
    N = N, K = K, X = X, R = R, S = S, A = A, B = B,
    date_laws = date_laws,
    W_pair = W_pair, masks = masks, rho_list = rho_list,
    pair_specs = pair_specs, probes_df = probes_df
  )
}

# Independent dense dummy-matrix GLM + dense Fisher information reference
# for any selection model / pooling / intercept specification.
selection_dense_glm_reference <- function(Z, cohort_index, w, A,
                                          U = rep(TRUE, nrow(Z)),
                                          rho = rep(1, nrow(Z)),
                                          cohorts = sort(unique(cohort_index)),
                                          date_responses = NULL,
                                          pair_responses = NULL,
                                          pair_specs,
                                          selection_model = "pair_logit",
                                          selection_pooling = "separate",
                                          glm_control = selection_glm_control_tight,
                                          has_intercept = NULL,
                                          probes = NULL,
                                          case_weights = rep(1, nrow(Z))) {
  N <- nrow(Z)
  cw <- as.double(case_weights)
  w_eff <- cw * as.double(w)
  d_spec <- selection_ns("extract_selection_design_spec")(Z, U, has_intercept = has_intercept)
  cat_obj <- selection_ns("partition_selection_cells")(
    pair_specs = pair_specs, cohorts = cohorts,
    selection_model = selection_model, selection_pooling = selection_pooling
  )
  cells <- cat_obj$cells
  block_specs <- cat_obj$block_specs
  pair_cohorts <- cat_obj$pair_cohorts
  C <- length(cells)
  M <- length(block_specs)
  K <- length(cohorts)
  ci <- match(cohort_index, cohorts)

  a_vec <- numeric(K)
  mass_U <- numeric(K)
  a_if <- matrix(0, N, K)
  for (k in seq_len(K)) {
    ik_u <- which(ci == k & U)
    mass_U[k] <- sum(w_eff[ik_u]) / sum(cw)
    a_vec[k] <- sum(w_eff[ik_u] * (A[ik_u] > 0)) / sum(w_eff[ik_u])
    a_if[ik_u, k] <- w_eff[ik_u] * (as.double(A[ik_u] > 0) - a_vec[k]) / (sum(w_eff[ik_u]) / N)
  }

  get_y <- function(cid, rr) {
    cell <- cells[[cid]]
    if (identical(selection_model, "independent_period_logit")) {
      as.double(date_responses[rr, cell$calendar_index])
    } else {
      key <- paste0(cell$pair$s_calendar_index, ":", cell$pair$t_calendar_index)
      as.double(pair_responses[[key]][rr])
    }
  }

  d_s <- d_spec$d_s
  slope_cols <- d_spec$slope_cols
  alpha_all <- if (d_spec$has_intercept) numeric(C) else numeric(0)
  beta_list <- vector("list", M)
  xi_alpha <- if (d_spec$has_intercept) matrix(0, N, C) else matrix(0, N, 0L)
  xi_beta_list <- vector("list", M)

  for (m in seq_len(M)) {
    bs <- block_specs[[m]]
    cids <- bs$cell_ids
    B_m <- bs$B_m
    p_m <- (if (d_spec$has_intercept) B_m else 0L) + d_s

    if (p_m == 0L) {
      beta_list[[m]] <- numeric(0)
      xi_beta_list[[m]] <- matrix(0, N, 0L)
      next
    }

    # Build explicit dense stacked design X_stack across cells in block m
    rows_list <- lapply(cids, function(cid) which(ci == cells[[cid]]$cohort & U))
    n_b_vec <- vapply(rows_list, length, integer(1))
    n_stack <- sum(n_b_vec)
    X_stack <- matrix(0, n_stack, p_m)
    y_stack <- numeric(n_stack)
    w_stack <- numeric(n_stack)
    orig_unit <- integer(n_stack)

    pos <- 0L
    for (b in seq_len(B_m)) {
      cid <- cids[b]
      rr <- rows_list[[b]]
      nb <- length(rr)
      idx <- seq.int(pos + 1L, pos + nb)
      if (d_spec$has_intercept) {
        X_stack[idx, b] <- 1
        if (d_s > 0L) X_stack[idx, B_m + seq_len(d_s)] <- Z[rr, slope_cols, drop = FALSE]
      } else {
        X_stack[idx, seq_len(d_s)] <- Z[rr, slope_cols, drop = FALSE]
      }
      y_stack[idx] <- get_y(cid, rr)
      w_stack[idx] <- w_eff[rr]
      orig_unit[idx] <- rr
      pos <- pos + nb
    }

    gf <- suppressWarnings(stats::glm.fit(
      x = X_stack, y = y_stack, weights = w_stack,
      start = rep(0, p_m), family = stats::binomial(link = "logit"),
      control = glm_control, intercept = FALSE
    ))
    theta_m <- as.double(gf$coefficients)
    eta_stack <- drop(X_stack %*% theta_m)
    p_stack <- 1 / (1 + exp(-eta_stack))

    # Dense Fisher information using ONE fitting weight w_stack
    H_dense <- crossprod(X_stack, w_stack * p_stack * (1 - p_stack) * X_stack) / sum(cw)
    raw_score_stack <- w_stack * (y_stack - p_stack) * X_stack
    # Sum repeated cell scores by ORIGINAL unit i in 1..N BEFORE solving H_dense
    unit_score_m <- matrix(0, N, p_m)
    for (r_idx in seq_len(n_stack)) {
      u_i <- orig_unit[r_idx]
      unit_score_m[u_i, ] <- unit_score_m[u_i, ] + raw_score_stack[r_idx, ]
    }
    xi_m <- t(solve(H_dense, t(unit_score_m)))

    if (d_spec$has_intercept) {
      alpha_all[cids] <- theta_m[seq_len(B_m)]
      xi_alpha[, cids] <- xi_m[, seq_len(B_m), drop = FALSE]
      beta_list[[m]] <- if (d_s > 0L) theta_m[B_m + seq_len(d_s)] else numeric(0)
      xi_beta_list[[m]] <- if (d_s > 0L) xi_m[, B_m + seq_len(d_s), drop = FALSE] else matrix(0, N, 0L)
    } else {
      beta_list[[m]] <- theta_m
      xi_beta_list[[m]] <- xi_m
    }
  }

  beta_all <- unlist(beta_list, use.names = FALSE)
  theta_all <- c(alpha_all, beta_all)
  xi_beta_all <- if (d_s > 0L) do.call(cbind, xi_beta_list) else matrix(0, N, 0L)
  theta_if <- cbind(xi_alpha, xi_beta_all)

  # Probe predictions
  targets <- list()
  if (is.null(probes)) {
    for (q in seq_along(pair_cohorts)) {
      pc <- pair_cohorts[[q]]
      rk <- which(ci == pc$cohort & U)
      z_eval <- if (d_s > 0L) colSums(w_eff[rk] * Z[rk, slope_cols, drop = FALSE]) / sum(w_eff[rk]) else numeric(0)
      rho_eval <- sum(w_eff[rk] * rho[rk]) / sum(w_eff[rk])
      targets[[length(targets) + 1L]] <- list(pair_cohort_id = q, cohort = pc$cohort, z_slope = z_eval, rho = rho_eval)
    }
  } else {
    Z_p <- probes$Z
    coh_p <- if (all(probes$cohort %in% cohorts)) match(probes$cohort, cohorts) else as.integer(probes$cohort)
    rho_p <- if (is.null(probes$rho)) rep(1, nrow(Z_p)) else as.double(probes$rho)
    Z_p_s <- if (ncol(Z_p) == d_s) Z_p else Z_p[, slope_cols, drop = FALSE]
    for (p_idx in seq_len(nrow(Z_p))) {
      k <- coh_p[p_idx]
      z_row <- as.double(Z_p_s[p_idx, ])
      for (q in which(vapply(pair_cohorts, function(pc) pc$cohort == k, logical(1)))) {
        targets[[length(targets) + 1L]] <- list(pair_cohort_id = q, cohort = k, z_slope = z_row, rho = rho_p[p_idx])
      }
    }
  }

  eval_c <- function(cid, z_s) {
    cell <- cells[[cid]]
    m <- cell$block_id
    if (!d_spec$has_intercept && d_s == 0L) return(0.5)
    eta <- (if (d_spec$has_intercept) alpha_all[cid] else 0) + (if (d_s > 0L) sum(z_s * beta_list[[m]]) else 0)
    1 / (1 + exp(-eta))
  }
  eta_if_c <- function(cid, z_s) {
    cell <- cells[[cid]]
    m <- cell$block_id
    if (!d_spec$has_intercept && d_s == 0L) return(numeric(N))
    ans <- if (d_spec$has_intercept) xi_alpha[, cid] else numeric(N)
    if (d_s > 0L) ans <- ans + drop(xi_beta_list[[m]] %*% z_s)
    ans
  }

  n_t <- length(targets)
  pair_prob <- numeric(n_t)
  mult <- numeric(n_t)
  pair_prob_if <- matrix(0, N, n_t)
  mult_if <- matrix(0, N, n_t)
  for (j in seq_len(n_t)) {
    tg <- targets[[j]]
    pc <- pair_cohorts[[tg$pair_cohort_id]]
    if (identical(selection_model, "independent_period_logit")) {
      ps <- eval_c(pc$s_cell_id, tg$z_slope)
      pt <- eval_c(pc$t_cell_id, tg$z_slope)
      pp <- ps * pt
      lp_if <- (1 - ps) * eta_if_c(pc$s_cell_id, tg$z_slope) + (1 - pt) * eta_if_c(pc$t_cell_id, tg$z_slope)
    } else {
      pp <- eval_c(pc$pair_cell_id, tg$z_slope)
      lp_if <- (1 - pp) * eta_if_c(pc$pair_cell_id, tg$z_slope)
    }
    xp <- pp * lp_if
    xa <- a_if[, pc$cohort]
    xm <- xa / (tg$rho * pp) - (a_vec[pc$cohort] / (tg$rho * pp^2)) * xp
    pair_prob[j] <- pp
    mult[j] <- a_vec[pc$cohort] / (tg$rho * pp)
    pair_prob_if[, j] <- xp
    mult_if[, j] <- xm
  }

  list(
    coefficients = theta_all,
    a = a_vec,
    pair_probabilities = pair_prob,
    multipliers = mult,
    vector = c(theta_all, a_vec, pair_prob, mult),
    analytic_stack = cbind(theta_if, a_if, pair_prob_if, mult_if),
    theta_if = theta_if,
    a_if = a_if,
    pair_prob_if = pair_prob_if,
    multiplier_if = mult_if
  )
}

testthat::test_that("frozen 32-atom independent_period_logit models match reference, dense GLM, FD, and mutation witnesses", {
  laws <- selection_build_32atom_laws()
  df32 <- data.frame(id = seq_len(laws$N), cohort = laws$K, x = laws$X)

  # Self-contained packaged reference cross-check
  ref_dir <- testthat::test_path("fixtures", "selection_models", "independent_reference")

  for (mode in c("within_cohort", "across_cohorts", "date_only", "separate_date_F")) {
    law <- laws$date_laws[[mode]]
    Z_use <- if (mode == "date_only") {
      selection_ns("selection_model_design")(~ 1, df32)
    } else {
      selection_ns("selection_model_design")(~ 1 + x, df32)
    }
    pooling_use <- if (mode %in% c("date_only", "separate_date_F")) "separate" else mode
    probes_Z <- if (mode == "date_only") {
      matrix(numeric(0), nrow(laws$probes_df), 0L)
    } else {
      matrix(laws$probes_df$x, ncol = 1L)
    }
    probes_arg <- list(cohort = laws$probes_df$cohort, Z = probes_Z)

    out <- selection_ns("evaluate_selection_table_if")(
      Z = Z_use, cohort_index = laws$K, w = law$w, A = laws$A,
      cohorts = c(0L, 1L), all_times = c(1, 2, 3), date_responses = laws$R,
      pair_specs = laws$pair_specs, selection_model = "independent_period_logit",
      selection_pooling = pooling_use, glm_control = selection_glm_control_tight,
      probes = probes_arg, chunk_rows = 7L, details = TRUE
    )
    dense <- selection_dense_glm_reference(
      Z = Z_use, cohort_index = laws$K, w = law$w, A = laws$A,
      cohorts = c(0L, 1L), date_responses = laws$R,
      pair_specs = laws$pair_specs, selection_model = "independent_period_logit",
      selection_pooling = pooling_use, glm_control = selection_glm_control_tight,
      probes = probes_arg
    )

    # Exact population score root recovery and agreement with dense dummy GLM
    testthat::expect_lt(max(abs(out$coefficients - law$theta)), 2e-11)
    testthat::expect_lt(max(abs(out$coefficients - dense$coefficients)), 1e-13)
    testthat::expect_lt(max(abs(out$pair_probabilities - dense$pair_probabilities)), 1e-14)
    testthat::expect_lt(max(abs(out$multipliers - dense$multipliers)), 1e-13)
    testthat::expect_lt(max(abs(out$analytic_stack - dense$analytic_stack)), 1e-12)
    testthat::expect_lt(max(abs(colMeans(out$analytic_stack))), 2e-11)

    # If on-disk CSV exists in workspace, verify exact match to frozen CSV
    if (dir.exists(ref_dir) && mode != "separate_date_F") {
      csv_path <- file.path(ref_dir, paste0(mode, "_if.csv"))
      if (file.exists(csv_path)) {
        ref_csv <- as.matrix(utils::read.csv(csv_path))
        testthat::expect_lt(max(abs(out$analytic_stack - ref_csv)), 1e-12)
      }
    }

    # Central finite-difference case-weight derivative check at h = 1e-4 and 5e-5
    for (h in c(1e-4, 5e-5)) {
      fd <- matrix(0, laws$N, ncol(out$analytic_stack))
      for (i in seq_len(laws$N)) {
        cp <- cm <- rep(1, laws$N)
        cp[i] <- 1 + h; cm[i] <- 1 - h
        vp <- selection_dense_glm_reference(
          Z = Z_use, cohort_index = laws$K, w = law$w, A = laws$A,
          cohorts = c(0L, 1L), date_responses = laws$R,
          pair_specs = laws$pair_specs, selection_model = "independent_period_logit",
          selection_pooling = pooling_use, glm_control = selection_glm_control_tight,
          probes = probes_arg, case_weights = cp
        )$vector
        vm <- selection_dense_glm_reference(
          Z = Z_use, cohort_index = laws$K, w = law$w, A = laws$A,
          cohorts = c(0L, 1L), date_responses = laws$R,
          pair_specs = laws$pair_specs, selection_model = "independent_period_logit",
          selection_pooling = pooling_use, glm_control = selection_glm_control_tight,
          probes = probes_arg, case_weights = cm
        )$vector
        fd[i, ] <- laws$N * (vp - vm) / (2 * h)
      }
      testthat::expect_lt(max(abs(fd - out$analytic_stack)), 3e-7)
    }

    # Off-A rows (singletons B \ A and zero-observation F \ B) must be nonzero on U = F
    off_a <- which(laws$A == 0)
    testthat::expect_gt(max(abs(out$theta_if[off_a, ])), 1e-3)
    testthat::expect_gt(max(abs(out$a_if[off_a, ])), 1e-3)
    testthat::expect_gt(max(abs(out$multiplier_if[off_a, ])), 1e-3)

    # Reused endpoint covariance: pair (1,2) and pair (1,3) share date 1
    pair_cov <- crossprod(out$pair_prob_if) / (laws$N^2)
    testthat::expect_gt(abs(pair_cov[1L, 2L]), 1e-4)

    # Mutation Witness 1: omitting endpoint t from pair (s, t) log-probability derivative
    pc1 <- out$pair_cohorts[[1L]] # k=1, pair (1,2)
    ps1 <- 1 / (1 + exp(-(out$alpha[pc1$s_cell_id] + (if (ncol(probes_Z) > 0L) out$blocks[[pc1$s_block_id]]$beta[1L] * (-1) else 0))))
    pt1 <- 1 / (1 + exp(-(out$alpha[pc1$t_cell_id] + (if (ncol(probes_Z) > 0L) out$blocks[[pc1$t_block_id]]$beta[1L] * (-1) else 0))))
    eta_s_if <- out$xi_alpha[, pc1$s_cell_id] + (if (ncol(probes_Z) > 0L) out$xi_beta_by_block[[pc1$s_block_id]][, 1L] * (-1) else 0)
    omitted_ep_if <- (ps1 * pt1) * (1 - ps1) * eta_s_if
    testthat::expect_gt(max(abs(out$pair_prob_if[, 1L] - omitted_ep_if)), 1e-3)

    # Mutation Witness 2: treating repeated date scores as independent unit copies
    # instead of summing by original unit i before Fisher solve / covariance
    cov_true <- crossprod(out$theta_if) / (laws$N^2)
    cov_copied <- matrix(0, ncol(out$theta_if), ncol(out$theta_if))
    for (t_idx in 1:3) {
      # Isolate only date t_idx cells (cells t_idx and 3 + t_idx)
      cs_t <- out$cell_scores
      keep_c <- c(t_idx, 3L + t_idx)
      cs_t[, setdiff(1:6, keep_c)] <- 0
      xi_a_t <- matrix(0, laws$N, 6L)
      xi_b_t <- lapply(out$blocks, function(bl) matrix(0, laws$N, bl$d_s))
      Z_s <- if (out$design_spec$d_s > 0L) Z_use[, out$design_spec$slope_cols, drop = FALSE] else matrix(0, laws$N, 0L)
      for (m in seq_along(out$blocks)) {
        bl <- out$blocks[[m]]
        bif_t <- selection_ns("evaluate_selection_block_if")(bl, Z_s, cs_t[, bl$cell_ids, drop = FALSE], laws$N)
        xi_a_t[, bl$cell_ids] <- bif_t$xi_alpha
        xi_b_t[[m]] <- bif_t$xi_beta
      }
      th_t <- cbind(xi_a_t, if (out$design_spec$d_s > 0L) do.call(cbind, xi_b_t) else matrix(0, laws$N, 0L))
      cov_copied <- cov_copied + crossprod(th_t) / (laws$N^2)
    }
    testthat::expect_gt(max(abs(cov_true - cov_copied)), 1e-4)

    # Mutation Witness 3: omitting non-endpoint date 3 score from pair (1,2) IF
    # Must fail (> 1e-3) when slope is shared across dates (within_cohort, across_cohorts)
    # and be 0 (< 1e-13) when parameters are separate across dates (date_only, separate_date_F).
    cs_no3 <- out$cell_scores
    cs_no3[, c(3L, 6L)] <- 0
    xi_a_no3 <- matrix(0, laws$N, 6L)
    xi_b_no3 <- lapply(out$blocks, function(bl) matrix(0, laws$N, bl$d_s))
    Z_s <- if (out$design_spec$d_s > 0L) Z_use[, out$design_spec$slope_cols, drop = FALSE] else matrix(0, laws$N, 0L)
    for (m in seq_along(out$blocks)) {
      bl <- out$blocks[[m]]
      bif_no3 <- selection_ns("evaluate_selection_block_if")(bl, Z_s, cs_no3[, bl$cell_ids, drop = FALSE], laws$N)
      xi_a_no3[, bl$cell_ids] <- bif_no3$xi_alpha
      xi_b_no3[[m]] <- bif_no3$xi_beta
    }
    eta1_no3 <- xi_a_no3[, pc1$s_cell_id] + (if (out$design_spec$d_s > 0L) xi_b_no3[[pc1$s_block_id]][, 1L] * (-1) else 0)
    eta2_no3 <- xi_a_no3[, pc1$t_cell_id] + (if (out$design_spec$d_s > 0L) xi_b_no3[[pc1$t_block_id]][, 1L] * (-1) else 0)
    p12_no3_if <- (ps1 * pt1) * ((1 - ps1) * eta1_no3 + (1 - pt1) * eta2_no3)
    no3_gap <- max(abs(out$pair_prob_if[, 1L] - p12_no3_if))
    if (mode %in% c("within_cohort", "across_cohorts")) {
      testthat::expect_gt(no3_gap, 1e-3)
    } else {
      testthat::expect_lt(no3_gap, 1e-13)
    }
  }
})

testthat::test_that("frozen 32-atom pair_logit models match reference across separate/within_cohort/across_cohorts and F/B/A", {
  laws <- selection_build_32atom_laws()
  df32 <- data.frame(id = seq_len(laws$N), cohort = laws$K, x = laws$X)
  Z_use <- selection_ns("selection_model_design")(~ 1 + x, df32)
  S_list <- list("1:2" = laws$S[, 1L], "1:3" = laws$S[, 2L], "2:3" = laws$S[, 3L])
  ref_dir <- testthat::test_path("fixtures", "selection_models", "pair_pooling_reference")

  for (mode in c("separate", "within_cohort", "across_cohorts")) {
    for (u in c("F", "B", "A")) {
      U_mask <- laws$masks[[u]]
      rho_vec <- laws$rho_list[[u]]
      rho_probes <- vapply(seq_len(nrow(laws$probes_df)), function(j) {
        idx <- which(laws$K == laws$probes_df$cohort[j] & laws$X == laws$probes_df$x[j])[1L]
        rho_vec[idx]
      }, numeric(1))
      probes_arg <- list(
        cohort = laws$probes_df$cohort,
        Z = matrix(laws$probes_df$x, ncol = 1L),
        rho = rho_probes
      )

      out <- selection_ns("evaluate_selection_table_if")(
        Z = Z_use, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
        U = U_mask, rho = rho_vec, universe = u,
        cohorts = c(0L, 1L), all_times = c(1, 2, 3), pair_responses = S_list,
        pair_specs = laws$pair_specs, selection_model = "pair_logit",
        selection_pooling = mode, glm_control = selection_glm_control_tight,
        probes = probes_arg, chunk_rows = 6L, details = TRUE
      )
      dense <- selection_dense_glm_reference(
        Z = Z_use, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
        U = U_mask, rho = rho_vec, cohorts = c(0L, 1L), pair_responses = S_list,
        pair_specs = laws$pair_specs, selection_model = "pair_logit",
        selection_pooling = mode, glm_control = selection_glm_control_tight,
        probes = probes_arg
      )

      testthat::expect_lt(max(abs(out$coefficients - dense$coefficients)), 1e-13)
      testthat::expect_lt(max(abs(out$pair_probabilities - dense$pair_probabilities)), 1e-14)
      testthat::expect_lt(max(abs(out$multipliers - dense$multipliers)), 1e-13)
      testthat::expect_lt(max(abs(out$analytic_stack - dense$analytic_stack)), 1e-12)
      testthat::expect_lt(max(abs(colMeans(out$analytic_stack))), 2e-11)

      if (dir.exists(ref_dir)) {
        csv_path <- file.path(ref_dir, paste0(mode, "_", u, "_if.csv"))
        if (file.exists(csv_path)) {
          ref_csv <- as.matrix(utils::read.csv(csv_path, header = FALSE))
          testthat::expect_lt(max(abs(out$analytic_stack - ref_csv)), 1e-12)
        }
      }

      # Finite-difference case-weight derivative check at h = 1e-4 and 5e-5
      steps_use <- if ((mode == "separate" && u == "F") || (mode == "within_cohort" && u == "B") || (mode == "across_cohorts" && u == "A")) c(1e-4, 5e-5) else 5e-5
      for (h in steps_use) {
        fd <- matrix(0, laws$N, ncol(out$analytic_stack))
        for (i in seq_len(laws$N)) {
          cp <- cm <- rep(1, laws$N)
          cp[i] <- 1 + h; cm[i] <- 1 - h
          vp <- selection_dense_glm_reference(
            Z = Z_use, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
            U = U_mask, rho = rho_vec, cohorts = c(0L, 1L), pair_responses = S_list,
            pair_specs = laws$pair_specs, selection_model = "pair_logit",
            selection_pooling = mode, glm_control = selection_glm_control_tight,
            probes = probes_arg, case_weights = cp
          )$vector
          vm <- selection_dense_glm_reference(
            Z = Z_use, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
            U = U_mask, rho = rho_vec, cohorts = c(0L, 1L), pair_responses = S_list,
            pair_specs = laws$pair_specs, selection_model = "pair_logit",
            selection_pooling = mode, glm_control = selection_glm_control_tight,
            probes = probes_arg, case_weights = cm
          )$vector
          fd[i, ] <- laws$N * (vp - vm) / (2 * h)
        }
        testthat::expect_lt(max(abs(fd - out$analytic_stack)), 3e-7)
      }

      # Cross-cohort slope contribution to cohort-0 pair probability IF from cohort-1 units:
      # must be nonzero ONLY in across_cohorts (and 0 on U = A where slope is 0 & S_e independent).
      cross_coh_max <- max(abs(out$pair_prob_if[laws$K == 1L & U_mask, 1L]))
      if (mode == "across_cohorts") {
        testthat::expect_gt(cross_coh_max, 1e-3)
      } else {
        testthat::expect_lt(cross_coh_max, 1e-14)
      }

      # Off-A rows: identically 0 on U = A, strictly nonzero on U in {F, B}
      off_a_max <- max(abs(out$analytic_stack[laws$A == 0, ]))
      if (u == "A") {
        testthat::expect_lt(off_a_max, 1e-14)
        testthat::expect_lt(max(abs(out$a_if)), 1e-14)
      } else {
        testthat::expect_gt(off_a_max, 1e-3)
      }
    }
  }
})

testthat::test_that("3-Bernoulli conditioning witness proves non-factorization on A and B and refuses date model off F", {
  # Exact enumeration of {0, 1}^3 under independent Bernoulli(1/2) marginal dates
  h <- as.matrix(expand.grid(R1 = c(0, 1), R2 = c(0, 1), R3 = c(0, 1)))
  in_A <- rowSums(h) >= 2
  in_B <- rowSums(h) >= 1
  # On F: marginal = 1/2, joint = 1/4 = (1/2)^2
  testthat::expect_equal(mean(h[, 1]), 0.5)
  testthat::expect_equal(mean(h[, 1] * h[, 2]), 0.25)
  # On A: marginal = 3/4, joint = 1/2 != (3/4)^2 = 9/16
  marg_A <- mean(h[in_A, 1])
  joint_A <- mean(h[in_A, 1] * h[in_A, 2])
  testthat::expect_equal(marg_A, 0.75)
  testthat::expect_equal(joint_A, 0.5)
  testthat::expect_gt(abs(joint_A - marg_A^2), 0.06)
  # On B: marginal = 4/7, joint = 2/7 != (4/7)^2 = 16/49
  marg_B <- mean(h[in_B, 1])
  joint_B <- mean(h[in_B, 1] * h[in_B, 2])
  testthat::expect_equal(marg_B, 4 / 7)
  testthat::expect_equal(joint_B, 2 / 7)
  testthat::expect_gt(abs(joint_B - marg_B^2), 0.04)

  # Backend must refuse independent_period_logit whenever universe != "F", !all(U), or rho != 1
  laws <- selection_build_32atom_laws()
  Z_use <- cbind(1, laws$X)
  for (u in c("B", "A")) {
    testthat::expect_error(
      selection_ns("evaluate_selection_table_if")(
        Z = Z_use, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
        U = laws$masks[[u]], rho = laws$rho_list[[u]], universe = u,
        cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
        selection_model = "independent_period_logit", selection_pooling = "separate",
        glm_control = selection_glm_control
      ),
      class = "cdid_cp04_unsupported_error"
    )
  }
  bad_rho <- rep(1, laws$N); bad_rho[1L] <- 1.1
  testthat::expect_error(
    selection_ns("evaluate_selection_table_if")(
      Z = Z_use, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
      U = rep(TRUE, laws$N), rho = bad_rho, universe = "F",
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "separate",
      glm_control = selection_glm_control
    ),
    class = "cdid_cp04_unsupported_error"
  )
})

testthat::test_that("formula intercept metadata, homogeneous no-intercept (~ 0 + ...), and zero-column (~ 0) branches work", {
  laws <- selection_build_32atom_laws()
  df32 <- data.frame(
    id = seq_len(laws$N),
    cohort = laws$K,
    x = laws$X,
    z2 = ifelse(laws$K == 1L, 0.5, -0.5) + 0.25 * laws$X,
    const1 = rep(1, laws$N)
  )
  law_w <- laws$date_laws$within_cohort$w

  # 1. Zero-column (~ 0 and ~ -1): eta = 0, p_marginal = 0.5, p_pair = 0.25 (date) or 0.5 (pair)
  for (fml in list(~ 0, ~ -1)) {
    Z0 <- selection_ns("selection_model_design")(fml, df32)
    spec0 <- selection_ns("extract_selection_design_spec")(Z0)
    testthat::expect_false(spec0$has_intercept)
    testthat::expect_identical(spec0$d_s, 0L)
    testthat::expect_identical(spec0$branch, "zero_col")

    out_date0 <- selection_ns("evaluate_selection_table_if")(
      Z = Z0, cohort_index = laws$K, w = law_w, A = laws$A,
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "across_cohorts",
      glm_control = selection_glm_control, details = TRUE
    )
    testthat::expect_length(out_date0$coefficients, 0L)
    testthat::expect_equal(ncol(out_date0$theta_if), 0L)
    testthat::expect_equal(unname(out_date0$pair_probabilities), rep(0.25, 6L))
    testthat::expect_equal(unname(out_date0$multipliers), rep(out_date0$a, each = 3L) / 0.25)
    testthat::expect_equal(max(abs(out_date0$pair_prob_if)), 0)
    # Empirical repeat-rate IF remains active and drives multiplier_if = a_if / 0.25
    testthat::expect_gt(max(abs(out_date0$multiplier_if)), 0.1)
    for (q in 1:6) {
      k_q <- out_date0$pair_cohorts[[q]]$cohort
      testthat::expect_equal(out_date0$multiplier_if[, q], out_date0$a_if[, k_q] / 0.25, tolerance = 1e-14)
    }

    S_list <- list("1:2" = laws$S[, 1L], "1:3" = laws$S[, 2L], "2:3" = laws$S[, 3L])
    out_pair0 <- selection_ns("evaluate_selection_table_if")(
      Z = Z0, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
      cohorts = c(0L, 1L), pair_responses = S_list, pair_specs = laws$pair_specs,
      selection_model = "pair_logit", selection_pooling = "within_cohort",
      glm_control = selection_glm_control
    )
    testthat::expect_equal(unname(out_pair0$pair_probabilities), rep(0.5, 6L))
    testthat::expect_equal(unname(out_pair0$multipliers), rep(out_pair0$a, each = 3L) / 0.5)
  }

  # 2. Homogeneous no-intercept (~ 0 + x + z2 and ~ -1 + x + z2): uncentered, no cell intercepts
  Z_no_int <- selection_ns("selection_model_design")(~ 0 + x + z2, df32)
  spec_no_int <- selection_ns("extract_selection_design_spec")(Z_no_int)
  testthat::expect_false(spec_no_int$has_intercept)
  testthat::expect_identical(spec_no_int$d_s, 2L)
  testthat::expect_identical(spec_no_int$branch, "no_intercept")

  for (pooling in c("separate", "within_cohort", "across_cohorts")) {
    out_ni <- selection_ns("evaluate_selection_table_if")(
      Z = Z_no_int, cohort_index = laws$K, w = law_w, A = laws$A,
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = pooling,
      glm_control = selection_glm_control_tight, chunk_rows = 9L, details = TRUE
    )
    dense_ni <- selection_dense_glm_reference(
      Z = Z_no_int, cohort_index = laws$K, w = law_w, A = laws$A,
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = pooling,
      glm_control = selection_glm_control_tight
    )
    expected_p <- switch(pooling, separate = 12L, within_cohort = 4L, across_cohorts = 2L)
    testthat::expect_length(out_ni$coefficients, expected_p)
    testthat::expect_length(out_ni$alpha, 0L)
    testthat::expect_lt(max(abs(out_ni$coefficients - dense_ni$coefficients)), 1e-13)
    testthat::expect_lt(max(abs(out_ni$analytic_stack - dense_ni$analytic_stack)), 1e-12)
  }

  # 3. Critical Formula Intercept Metadata Witness: ~ 0 + const1 vs ~ 1
  # Under ~ 0 + const1 (has_intercept = FALSE), const1 = 1 is a single shared
  # homogeneous slope across all 3 dates in each cohort (2 parameters total in within_cohort).
  # Under ~ 1 (has_intercept = TRUE), each (cohort, date) cell gets its own free
  # local intercept (6 parameters total in within_cohort).
  Z_homog_const <- selection_ns("selection_model_design")(~ 0 + const1, df32)
  Z_local_int <- selection_ns("selection_model_design")(~ 1, df32)
  out_homog <- selection_ns("evaluate_selection_table_if")(
    Z = Z_homog_const, cohort_index = laws$K, w = law_w, A = laws$A,
    cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
    selection_model = "independent_period_logit", selection_pooling = "within_cohort",
    glm_control = selection_glm_control_tight
  )
  out_local <- selection_ns("evaluate_selection_table_if")(
    Z = Z_local_int, cohort_index = laws$K, w = law_w, A = laws$A,
    cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
    selection_model = "independent_period_logit", selection_pooling = "within_cohort",
    glm_control = selection_glm_control_tight
  )
  testthat::expect_length(out_homog$coefficients, 2L)
  testthat::expect_length(out_local$coefficients, 6L)
  testthat::expect_gt(max(abs(out_homog$pair_probabilities - out_local$pair_probabilities)), 0.05)
})

testthat::test_that("144-unit reference fixture verifies deduplication, permutation invariance, moment contractions, and refusals", {
  cp04_ref <- new.env(parent = asNamespace("stats"))
  invisible(capture.output(sys.source(testthat::test_path("reference", "cp04_finite_reference.R"), cp04_ref)))

  N <- cp04_ref$N
  g_vals <- c(3, 5, Inf)[cp04_ref$dat$k]
  df144 <- data.frame(id = cp04_ref$ids, g = g_vals, x = cp04_ref$Z[, 2L], w = cp04_ref$w)
  Z144 <- selection_ns("selection_model_design")(~ 1 + x, df144)

  # Active pairs from the 4 links in cp04_fixture: (1,2) and (2,4) across cohorts {3, 5, Inf}
  # Note: date 3 is never used in any link pair; only dates {1, 2, 4} are active!
  raw_specs <- list(
    # Link 1 (NYT12): pair (1,2), cohorts 1, 2, 3
    list(cohort = 1L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 2L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 3L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    # Link 2 (NYT24): pair (2,4), cohorts 1, 2, 3
    list(cohort = 1L, pair = list(s = 2, t = 4, s_calendar_index = 2L, t_calendar_index = 4L)),
    list(cohort = 2L, pair = list(s = 2, t = 4, s_calendar_index = 2L, t_calendar_index = 4L)),
    list(cohort = 3L, pair = list(s = 2, t = 4, s_calendar_index = 2L, t_calendar_index = 4L)),
    # Link 3 (NT12): duplicate pair (1,2) on cohorts 1, 3
    list(cohort = 1L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 3L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    # Link 4 (horizon12): duplicate pair (1,2) on cohorts 1, 3, 2
    list(cohort = 1L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 3L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 2L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L))
  )
  S_list144 <- list("1:2" = cp04_ref$S[, 1L], "2:4" = cp04_ref$S[, 2L])

  for (sel_model in c("independent_period_logit", "pair_logit")) {
    for (pooling in c("separate", "within_cohort", "across_cohorts")) {
      out <- selection_ns("evaluate_selection_table_if")(
        Z = Z144, cohort_index = g_vals, w = cp04_ref$w, A = cp04_ref$A,
        cohorts = c(3, 5, Inf), all_times = 1:4,
        date_responses = cp04_ref$obs, pair_responses = S_list144,
        pair_specs = raw_specs, selection_model = sel_model,
        selection_pooling = pooling, glm_control = selection_glm_control,
        chunk_rows = 19L, details = TRUE
      )
      dense <- selection_dense_glm_reference(
        Z = Z144, cohort_index = g_vals, w = cp04_ref$w, A = cp04_ref$A,
        cohorts = c(3, 5, Inf), date_responses = cp04_ref$obs, pair_responses = S_list144,
        pair_specs = raw_specs, selection_model = sel_model,
        selection_pooling = pooling, glm_control = selection_glm_control
      )

      # Active cell deduplication check:
      # 9 unique date cells (3 cohorts x dates {1, 2, 4}; date 2 is shared by (1,2) and (2,4), date 3 is unused)
      # 6 unique pair cells (3 cohorts x pairs {(1,2), (2,4)})
      expected_cells <- if (sel_model == "independent_period_logit") 9L else 6L
      testthat::expect_length(out$cells, expected_cells)
      testthat::expect_length(out$pair_cohorts, 6L)
      if (sel_model == "independent_period_logit") {
        used_dates <- sort(unique(vapply(out$cells, `[[`, integer(1), "calendar_index")))
        testthat::expect_identical(used_dates, c(1L, 2L, 4L))
      }

      testthat::expect_lt(max(abs(out$coefficients - dense$coefficients)), 1e-12)
      testthat::expect_lt(max(abs(out$analytic_stack - dense$analytic_stack)), 1e-11)

      # Permutation and duplication invariance of pair_specs (bit-identical at same chunk_rows)
      out_rev <- selection_ns("evaluate_selection_table_if")(
        Z = Z144, cohort_index = g_vals, w = cp04_ref$w, A = cp04_ref$A,
        cohorts = c(3, 5, Inf), all_times = 1:4,
        date_responses = cp04_ref$obs, pair_responses = S_list144,
        pair_specs = rev(raw_specs[1:6]), selection_model = sel_model,
        selection_pooling = pooling, glm_control = selection_glm_control,
        chunk_rows = 19L, details = TRUE
      )
      testthat::expect_identical(out$coefficients, out_rev$coefficients)
      testthat::expect_identical(out$analytic_stack, out_rev$analytic_stack)

      # Chunk-size invariance across different chunk_rows
      out_ch23 <- selection_ns("evaluate_selection_table_if")(
        Z = Z144, cohort_index = g_vals, w = cp04_ref$w, A = cp04_ref$A,
        cohorts = c(3, 5, Inf), all_times = 1:4,
        date_responses = cp04_ref$obs, pair_responses = S_list144,
        pair_specs = raw_specs[1:6], selection_model = sel_model,
        selection_pooling = pooling, glm_control = selection_glm_control,
        chunk_rows = 23L, details = TRUE
      )
      testthat::expect_equal(out$coefficients, out_ch23$coefficients, tolerance = 1e-13)
      testthat::expect_equal(out$analytic_stack, out_ch23$analytic_stack, tolerance = 1e-13)

      # Verify Eq. (7) profiled moment contraction matches direct (xi_alpha, xi_beta) contraction
      acc <- NULL
      direct_loading <- numeric(N)
      for (q in c(1L, 4L)) {
        pc <- out$pair_cohorts[[q]]
        obs_q <- which(g_vals == c(3, 5, Inf)[pc$cohort] & S_list144[[if (q <= 2L) "1:2" else "2:4"]] == 1)
        ow <- cp04_ref$w[obs_q] * 0.75
        acc <- selection_ns("accumulate_selection_table_moments")(out, Z144, q, obs_q, ow, acc = acc, sign = -1)
      }
      prof_loading <- selection_ns("evaluate_selection_table_loading")(out, Z144, acc, include_repeat_rate = TRUE)
      t_full <- c(
        unlist(lapply(acc$blocks, `[[`, "t_alpha"), use.names = FALSE),
        unlist(lapply(acc$blocks, `[[`, "t_beta"), use.names = FALSE)
      )
      direct_loading <- drop(out$a_if %*% acc$t_a) + drop(out$theta_if %*% t_full)
      testthat::expect_lt(max(abs(prof_loading - direct_loading)), 1e-13)
    }
  }

  # Mutation witness on 144-unit fixture: dropping singleton date observations (B \ A)
  # corrupts independent_period_logit fits and breaks agreement with full-F dense GLM
  obs_A_only <- cp04_ref$obs
  obs_A_only[cp04_ref$A == 0, ] <- 0
  out_full_B <- selection_ns("evaluate_selection_table_if")(
    Z = Z144, cohort_index = g_vals, w = cp04_ref$w, A = cp04_ref$A,
    cohorts = c(3, 5, Inf), all_times = 1:4, date_responses = cp04_ref$obs,
    pair_specs = raw_specs[1:6], selection_model = "independent_period_logit",
    selection_pooling = "within_cohort", glm_control = selection_glm_control
  )
  out_zero_sing <- selection_ns("evaluate_selection_table_if")(
    Z = Z144, cohort_index = g_vals, w = cp04_ref$w, A = cp04_ref$A,
    cohorts = c(3, 5, Inf), all_times = 1:4, date_responses = obs_A_only,
    pair_specs = raw_specs[1:6], selection_model = "independent_period_logit",
    selection_pooling = "within_cohort", glm_control = selection_glm_control
  )
  testthat::expect_gt(max(abs(out_full_B$coefficients - out_zero_sing$coefficients)), 1e-2)
  testthat::expect_gt(max(abs(out_full_B$pair_probabilities - out_zero_sing$pair_probabilities)), 1e-2)
})

testthat::test_that("rank, separation, boundary, convergence, and support refusals honor actual fitted model", {
  laws <- selection_build_32atom_laws()
  df32 <- data.frame(id = seq_len(laws$N), cohort = laws$K, x = laws$X)
  Z_int <- selection_ns("selection_model_design")(~ 1 + x, df32)

  # 1. Rank deficiency: duplicate slope column
  Z_collinear <- cbind(1, laws$X, laws$X)
  attr(Z_collinear, "has_intercept") <- TRUE
  testthat::expect_error(
    selection_ns("evaluate_selection_table_if")(
      Z = Z_collinear, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "across_cohorts",
      glm_control = selection_glm_control
    ),
    class = "cdid_cp04_rank_error"
  )

  # 2. Rank deficiency: cohort-constant covariate in within_cohort with cell intercepts
  Z_coh_const <- cbind(1, as.double(laws$K))
  attr(Z_coh_const, "has_intercept") <- TRUE
  testthat::expect_error(
    selection_ns("evaluate_selection_table_if")(
      Z = Z_coh_const, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "within_cohort",
      glm_control = selection_glm_control
    ),
    class = "cdid_cp04_rank_error"
  )

  # 3. Coordinate separation on cell indicator vs homogeneous no-intercept model:
  # Make date 1 in cohort 0 all-failure (R[, 1] == 0 on cohort 0).
  R_cell0_zero <- laws$R
  R_cell0_zero[laws$K == 0L, 1L] <- 0
  # Under ~ 1 + x (has_intercept = TRUE), cell (k=0, t=1) has its own indicator column -> separation!
  testthat::expect_error(
    selection_ns("evaluate_selection_table_if")(
      Z = Z_int, cohort_index = laws$K, w = laws$date_laws$within_cohort$w, A = laws$A,
      cohorts = c(0L, 1L), date_responses = R_cell0_zero, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "within_cohort",
      glm_control = selection_glm_control
    ),
    class = "cdid_cp04_separation_error"
  )
  # Under ~ 0 + x (has_intercept = FALSE) in within_cohort, there is NO cell indicator column,
  # and x takes both -1 and +1 on successes in dates 2 and 3 -> finite root exists and succeeds!
  Z_homog <- selection_ns("selection_model_design")(~ 0 + x, df32)
  fit_homog_ok <- selection_ns("evaluate_selection_table_if")(
    Z = Z_homog, cohort_index = laws$K, w = laws$date_laws$within_cohort$w, A = laws$A,
    cohorts = c(0L, 1L), date_responses = R_cell0_zero, pair_specs = laws$pair_specs,
    selection_model = "independent_period_logit", selection_pooling = "within_cohort",
    glm_control = selection_glm_control
  )
  testthat::expect_true(all(is.finite(fit_homog_ok$coefficients)))

  # 4. Signed slope separation across block: Y == (x > 0)
  R_slope_sep <- matrix(as.double(laws$X > 0), laws$N, 3L)
  testthat::expect_error(
    selection_ns("evaluate_selection_table_if")(
      Z = Z_homog, cohort_index = laws$K, w = laws$W_pair, A = laws$A,
      cohorts = c(0L, 1L), date_responses = R_slope_sep, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "across_cohorts",
      glm_control = selection_glm_control
    ),
    class = "cdid_cp04_separation_error"
  )

  # 5. Convergence refusal under maxit = 1
  bad_gc <- selection_glm_control
  bad_gc$maxit <- 1L
  testthat::expect_error(
    selection_ns("evaluate_selection_table_if")(
      Z = Z_int, cohort_index = laws$K, w = laws$date_laws$within_cohort$w, A = laws$A,
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "within_cohort",
      glm_control = bad_gc
    ),
    class = "cdid_cp04_convergence_error"
  )

  # 6. Support refusal when cohort repeat-observation mass a_k == 0
  A_no_k0 <- laws$A
  A_no_k0[laws$K == 0L] <- 0
  testthat::expect_error(
    selection_ns("evaluate_selection_table_if")(
      Z = Z_int, cohort_index = laws$K, w = laws$date_laws$within_cohort$w, A = A_no_k0,
      cohorts = c(0L, 1L), date_responses = laws$R, pair_specs = laws$pair_specs,
      selection_model = "independent_period_logit", selection_pooling = "within_cohort",
      glm_control = selection_glm_control
    ),
    class = "cdid_cp04_support_error"
  )
})

testthat::test_that("installed-namespace sealed CP04 reader, cp04_date_rows, fit_selection_models, replay, and covariance stream work", {
  testthat::skip_if_not(
    selection_has_installed_cdid(),
    "Installed cdid namespace required for C-backed panel/frame/link owner integration."
  )

  cp04_ref <- new.env(parent = asNamespace("stats"))
  invisible(capture.output(sys.source(testthat::test_path("reference", "cp04_finite_reference.R"), cp04_ref)))
  ctrl <- selection_ns("ipw_owner_control")(
    max_chunk_rows = 23L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7
  )

  build_owner_fixture <- function(order = seq_len(cp04_ref$N), Z = cp04_ref$Z,
                                  U = rep(TRUE, cp04_ref$N), rho = rep(1, cp04_ref$N),
                                  universe = "F", justification = NULL, has_intercept = NULL) {
    g <- c(3, 5, Inf)[cp04_ref$dat$k]
    long <- data.frame(
      id = rep(cp04_ref$ids, each = 4L),
      time = rep(1:4, cp04_ref$N),
      y = as.vector(t(cp04_ref$Y)),
      g = rep(g, each = 4L)
    )
    zero <- which(cp04_ref$A == 0 & rowSums(cp04_ref$obs) == 0)
    long <- long[!long$id %in% zero[-1L], ]
    panel <- selection_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
    frame <- selection_ns("construct_ipw_frame_owner")(
      panel, cp04_ref$ids[order], cp04_ref$V[order, , drop = FALSE],
      weights = cp04_ref$w[order], control = ctrl, frame_g = g[order], allow_frame_only = TRUE
    )
    Z_ord <- Z[order, , drop = FALSE]
    attr(Z_ord, "has_intercept") <- attr(Z, "has_intercept")
    attr(Z_ord, "assign") <- attr(Z, "assign")
    attr(Z_ord, "intercept") <- attr(Z, "intercept")
    reader <- selection_ns("construct_cp04_frame_reader")(
      frame, cp04_ref$ids[order], Z_ord, g[order], cp04_ref$w[order], ctrl,
      U = U[order], rho = rho[order], universe = universe,
      justification = justification, has_intercept = has_intercept
    )
    make_lv <- function(s, t, K, id) {
      selection_ns("construct_ipw_link_view")(
        frame, s, t, g = 3, K = K, control = ctrl, link_id = id,
        functional_key = paste0("cmpcat:reference-F-law:g1:K", if (length(K) == 1L) 2L else 1L)
      )
    }
    views <- list(
      make_lv(1, 2, c(5, Inf), "NYT12"),
      make_lv(2, 4, c(5, Inf), "NYT24"),
      make_lv(1, 2, Inf, "NT12"),
      make_lv(1, 2, c(Inf, 5), "horizon12")
    )
    list(panel = panel, frame = frame, reader = reader, views = views, g = g)
  }

  fx <- build_owner_fixture()

  # 1. Verify cp04_date_rows recovers exact full-B marginal date indicators
  # including singletons (B \ A, a_index == 0) while keeping zero-observers and F \ B at 0
  for (t_cal in 1:4) {
    got_R <- numeric(cp04_ref$N)
    for (start in seq.int(1L, cp04_ref$N, by = 17L)) {
      rr <- seq.int(start, min(cp04_ref$N, start + 16L))
      got_R[rr] <- selection_ns("cp04_date_rows")(fx$reader, rr, list(t = t_cal, calendar_index = t_cal))$R
    }
    testthat::expect_identical(got_R, as.double(cp04_ref$obs[, t_cal]))
  }
  singletons <- which(cp04_ref$A == 0 & rowSums(cp04_ref$obs) == 1)
  zero_obs <- which(cp04_ref$A == 0 & rowSums(cp04_ref$obs) == 0)
  testthat::expect_gt(length(singletons), 0L)
  testthat::expect_gt(length(zero_obs), 0L)

  # 2. Preserve existing separate positive-column pair_logit CP04 fits bit-for-bit
  cp04_fit <- selection_ns("fit_cp04_corollary1")(fx$reader, fx$views, ctrl, selection_glm_control)
  sel_sep_pair <- selection_ns("fit_selection_models")(
    fx$reader, fx$views, selection_model = "pair_logit", selection_pooling = "separate",
    control = ctrl, glm_control = selection_glm_control
  )
  testthat::expect_s3_class(sel_sep_pair, "cdid_selection_models")
  testthat::expect_true(selection_ns("ipw_owner_seal_ok")(sel_sep_pair, "selection_models"))
  testthat::expect_equal(unname(selection_ns("cp04_payload")(sel_sep_pair)["entries"]),
                         sel_sep_pair$accounting$owned_entries)

  # Match every selection cell against fit_cp04_corollary1$selections
  for (j in seq_along(cp04_fit$selection_specs)) {
    sp <- cp04_fit$selection_specs[[j]]
    qid <- selection_ns("resolve_pair_cohort_id")(sel_sep_pair, sp)
    bl <- sel_sep_pair$blocks[[sel_sep_pair$pair_cohorts[[qid]]$block_id]]
    testthat::expect_identical(bl$branch, "cp04_separate_pair")
    testthat::expect_identical(bl$coefficients, cp04_fit$selections[[j]]$coefficients)
    testthat::expect_identical(bl$H_chol, cp04_fit$selections[[j]]$H_chol)
    testthat::expect_identical(bl$pi, cp04_fit$selections[[j]]$pi)
    testthat::expect_identical(bl$r, cp04_fit$selections[[j]]$r)
  }

  # Existing CP04 positive-column semantics include homogeneous designs and a
  # supplied intercept outside column 1, with logical public metadata.
  Z_hom <- matrix(cp04_ref$Z[, 2L], ncol = 1L)
  attr(Z_hom, "intercept") <- FALSE
  attr(Z_hom, "assign") <- 1L
  Z_perm <- cp04_ref$Z[, c(2L, 1L), drop = FALSE]
  attr(Z_perm, "intercept") <- TRUE
  attr(Z_perm, "assign") <- c(1L, 0L)
  for (Z_legacy in list(Z_hom, Z_perm)) {
    fx_legacy <- build_owner_fixture(Z = Z_legacy)
    old <- selection_ns("fit_cp04_corollary1")(
      fx_legacy$reader, fx_legacy$views, ctrl, selection_glm_control
    )
    new <- selection_ns("fit_selection_models")(
      fx_legacy$reader, fx_legacy$views, "pair_logit", "separate", ctrl, selection_glm_control
    )
    for (j in seq_along(old$selection_specs)) {
      q <- selection_ns("resolve_pair_cohort_id")(new, old$selection_specs[[j]])
      bl <- new$blocks[[new$pair_cohorts[[q]]$block_id]]
      testthat::expect_identical(bl$branch, "cp04_separate_pair")
      testthat::expect_identical(bl$coefficients, old$selections[[j]]$coefficients)
      testthat::expect_identical(bl$H_chol, old$selections[[j]]$H_chol)
      testthat::expect_identical(bl$pi, old$selections[[j]]$pi)
      testthat::expect_identical(bl$r, old$selections[[j]]$r)
    }
    dense_legacy <- selection_dense_glm_reference(
      Z_legacy, fx_legacy$g, cp04_ref$w, cp04_ref$A,
      cohorts = c(3, 5, Inf), pair_responses = list("1:2" = cp04_ref$S[, 1L], "2:4" = cp04_ref$S[, 2L]),
      pair_specs = old$selection_specs, glm_control = selection_glm_control
    )
    actual_if <- matrix(NA_real_, cp04_ref$N, ncol(dense_legacy$analytic_stack))
    selection_ns("replay_selection_if")(new, 23L, callback = function(b) {
      actual_if[b$rows, ] <<- cbind(b$theta_if, b$a_if, b$pair_prob_if, b$multiplier_if)
    })
    testthat::expect_lt(max(abs(actual_if - dense_legacy$analytic_stack)), 1e-12)
    acc_legacy <- selection_ns("add_selection_moments")(new, 1L, 1:20, rep(0.5, 20L))
    tt_legacy <- c(unlist(lapply(acc_legacy$blocks, `[[`, "t_alpha"), use.names = FALSE),
                   unlist(lapply(acc_legacy$blocks, `[[`, "t_beta"), use.names = FALSE))
    expected_load <- drop(dense_legacy$theta_if %*% tt_legacy) + drop(dense_legacy$a_if %*% acc_legacy$t_a)
    actual_load <- numeric(cp04_ref$N)
    for (start in seq.int(1L, cp04_ref$N, by = 23L)) {
      rr <- seq.int(start, min(cp04_ref$N, start + 22L))
      actual_load[rr] <- selection_ns("evaluate_pair_selection_loading")(new, rr, acc_legacy)
    }
    testthat::expect_lt(max(abs(actual_load - expected_load)), 1e-12)
  }

  # 3. Verify all 6 (selection_model x selection_pooling) combinations on the sealed reader
  raw_specs <- list(
    list(cohort = 1L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 2L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 3L, pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L)),
    list(cohort = 1L, pair = list(s = 2, t = 4, s_calendar_index = 2L, t_calendar_index = 4L)),
    list(cohort = 2L, pair = list(s = 2, t = 4, s_calendar_index = 2L, t_calendar_index = 4L)),
    list(cohort = 3L, pair = list(s = 2, t = 4, s_calendar_index = 2L, t_calendar_index = 4L))
  )
  S_list144 <- list("1:2" = cp04_ref$S[, 1L], "2:4" = cp04_ref$S[, 2L])

  for (sel_model in c("independent_period_logit", "pair_logit")) {
    for (pooling in c("separate", "within_cohort", "across_cohorts")) {
      fit_m <- selection_ns("fit_selection_models")(
        fx$reader, fx$views, selection_model = sel_model, selection_pooling = pooling,
        control = ctrl, glm_control = selection_glm_control
      )
      tbl_m <- selection_ns("evaluate_selection_table_if")(
        Z = cp04_ref$Z, cohort_index = fx$g, w = cp04_ref$w, A = cp04_ref$A,
        cohorts = c(3, 5, Inf), all_times = 1:4,
        date_responses = cp04_ref$obs, pair_responses = S_list144,
        pair_specs = raw_specs, selection_model = sel_model,
        selection_pooling = pooling, glm_control = selection_glm_control,
        details = TRUE
      )

      testthat::expect_equal(fit_m$coefficients, tbl_m$coefficients, tolerance = 1e-13)
      testthat::expect_false(fit_m$accounting$full_F_by_cell_IF_retained)
      testthat::expect_false(fit_m$accounting$global_dummy_matrix_constructed)
      testthat::expect_false(fit_m$accounting$inverse_matrices_retained)

      # Replay across different chunk sizes and compare against table IF & streaming covariance
      for (ch_sz in c(7L, 23L)) {
        replayed_mult <- matrix(NA_real_, cp04_ref$N, length(fit_m$pair_cohorts))
        replayed_theta <- matrix(NA_real_, cp04_ref$N, length(fit_m$coefficients))
        info_r <- selection_ns("replay_selection_if")(
          fit_m, chunk_rows = ch_sz, pair_batch_size = 2L,
          callback = function(b) {
            replayed_mult[b$rows, b$pair_cohort_ids] <<- b$multiplier_if
            if ( all(b$pair_cohort_ids == c(1L, 2L)) ) replayed_theta[b$rows, ] <<- b$theta_if
          }
        )
        testthat::expect_equal(unname(replayed_mult), unname(tbl_m$multiplier_if), tolerance = 1e-13)
      }

      full_stack <- matrix(NA_real_, cp04_ref$N, ncol(tbl_m$analytic_stack))
      selection_ns("replay_selection_if")(
        fit_m, chunk_rows = 23L, callback = function(b) {
          full_stack[b$rows, ] <<- cbind(b$theta_if, b$a_if, b$pair_prob_if, b$multiplier_if)
        }
      )
      testthat::expect_equal(unname(full_stack), unname(tbl_m$analytic_stack), tolerance = 1e-13)
      cov_stream <- selection_ns("stream_selection_covariance")(fit_m, chunk_rows = 13L)
      expected_cov <- crossprod(tbl_m$analytic_stack) / (cp04_ref$N^2)
      testthat::expect_equal(unname(cov_stream$covariance), unname(expected_cov), tolerance = 1e-13)

      # Core scale check
      cov_core <- selection_ns("stream_selection_covariance")(fit_m, chunk_rows = 13L, scale = "core")
      testthat::expect_equal(cov_core$normalization_count, sum(cp04_ref$A))
      testthat::expect_equal(unname(cov_core$covariance), unname(expected_cov), tolerance = 1e-13)

      # Moment accumulation & profiled loading evaluation on chunks
      acc <- selection_ns("add_selection_moments")(fit_m, 1L, 1:20, rep(0.5, 20L), sign = -1)
      acc <- selection_ns("add_selection_moments")(fit_m, 4L, 21:40, rep(0.25, 20L), acc = acc, sign = -1)
      load_full <- numeric(cp04_ref$N)
      for (start in seq.int(1L, cp04_ref$N, by = 19L)) {
        rr <- seq.int(start, min(cp04_ref$N, start + 18L))
        load_full[rr] <- selection_ns("evaluate_pair_selection_loading")(fit_m, rr, acc)
      }
      t_vec <- c(
        unlist(lapply(acc$blocks, `[[`, "t_alpha"), use.names = FALSE),
        unlist(lapply(acc$blocks, `[[`, "t_beta"), use.names = FALSE)
      )
      expected_load <- drop(tbl_m$a_if %*% acc$t_a) + drop(tbl_m$theta_if %*% t_vec)
      testthat::expect_equal(load_full, expected_load, tolerance = 1e-13)
    }
  }

  # 4. Zero-column (~ 0) and no-intercept (~ 0 + x) through construct_cp04_frame_reader
  df144 <- data.frame(x = cp04_ref$Z[, 2L])
  Z0_144 <- selection_ns("selection_model_design")(~ 0, df144)
  fx0 <- build_owner_fixture(Z = Z0_144, has_intercept = FALSE)
  testthat::expect_equal(ncol(fx0$reader$state$Z), 0L)
  fit0 <- selection_ns("fit_selection_models")(
    fx0$reader, fx0$views, selection_model = "independent_period_logit",
    selection_pooling = "within_cohort", control = ctrl, glm_control = selection_glm_control
  )
  testthat::expect_length(fit0$coefficients, 0L)
  pred0 <- selection_ns("predict_period_pair_selection")(fit0, 1:10, 1L)
  testthat::expect_equal(pred0$prob, rep(0.25, 10L))

  Z_ni_144 <- selection_ns("selection_model_design")(~ 0 + x, df144)
  fx_ni <- build_owner_fixture(Z = Z_ni_144, has_intercept = FALSE)
  fit_ni <- selection_ns("fit_selection_models")(
    fx_ni$reader, fx_ni$views, selection_model = "independent_period_logit",
    selection_pooling = "across_cohorts", control = ctrl, glm_control = selection_glm_control
  )
  testthat::expect_length(fit_ni$coefficients, 1L)
  testthat::expect_length(fit_ni$alpha, 0L)

  # 5. Integrity & resource budget preflight refusals before callback
  fit_ok <- selection_ns("fit_selection_models")(
    fx$reader, fx$views, selection_model = "independent_period_logit",
    selection_pooling = "within_cohort", control = ctrl, glm_control = selection_glm_control
  )
  tampered <- fit_ok
  tampered$alpha[1L] <- tampered$alpha[1L] + 1
  testthat::expect_error(
    selection_ns("evaluate_selection_if_chunk")(tampered, 1:5),
    class = "cdid_cp04_integrity_error"
  )

  limited <- function(field) {
    x <- fit_ok
    x$ref <- NULL
    x$control[[field]] <- 1
    selection_ns("ipw_owner_seal")(x, "selection_models")
  }
  for (field in c("max_output_entries", "max_work_ops", "max_matrix_entries", "max_retained_entries")) {
    called <- FALSE
    testthat::expect_error(
      selection_ns("replay_selection_if")(limited(field), 7L, 2L, function(b) called <<- TRUE),
      class = "cdid_ipw_owner_bound_error"
    )
    testthat::expect_false(called)
    testthat::expect_error(
      selection_ns("stream_selection_covariance")(limited(field), 7L),
      class = "cdid_ipw_owner_bound_error"
    )
  }
})

# These are pure numerical/metadata probes, not sealed-owner substitutes.
testthat::test_that("logical public intercept metadata and original CP04 coordinates survive separate pair fits", {
  laws <- selection_build_32atom_laws()
  base <- cbind(laws$X, 1)
  attr(base, "assign") <- c(1L, 0L) # actual intercept is column 2
  designs <- list(base)
  for (flag in list(TRUE, 1L, 1)) {
    z <- base
    attr(z, "intercept") <- flag
    designs[[length(designs) + 1L]] <- z
  }
  for (flag in list(FALSE, 0L, 0)) {
    z <- matrix(laws$X, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    attr(z, "intercept") <- flag
    designs[[length(designs) + 1L]] <- z
  }
  for (Z in designs) {
    out <- selection_ns("evaluate_selection_table_if")(
      Z, laws$K, laws$W_pair, laws$A, cohorts = c(0L, 1L),
      pair_responses = list("1:2" = laws$S[, 1L], "1:3" = laws$S[, 2L], "2:3" = laws$S[, 3L]),
      pair_specs = laws$pair_specs, selection_model = "pair_logit",
      selection_pooling = "separate", glm_control = selection_glm_control,
      chunk_rows = 7L, details = TRUE
    )
    dense <- selection_dense_glm_reference(
      Z, laws$K, laws$W_pair, laws$A, cohorts = c(0L, 1L),
      pair_responses = list("1:2" = laws$S[, 1L], "1:3" = laws$S[, 2L], "2:3" = laws$S[, 3L]),
      pair_specs = laws$pair_specs, glm_control = selection_glm_control
    )
    testthat::expect_lt(max(abs(out$analytic_stack - dense$analytic_stack)), 1e-12)
    testthat::expect_identical(out$design_spec$has_intercept, ncol(Z) == 2L)
    for (b in seq_along(out$blocks)) {
      cell <- out$cells[[b]]
      rr <- which(laws$K == cell$cohort - 1L)
      key <- paste0(cell$pair$s_calendar_index, ":", cell$pair$t_calendar_index)
      yy <- list("1:2" = laws$S[, 1L], "1:3" = laws$S[, 2L], "2:3" = laws$S[, 3L])[[key]][rr]
      X <- Z[rr, , drop = FALSE]
      colnames(X) <- NULL
      old <- selection_ns("cp04_logit")(X, yy, laws$W_pair[rr], selection_glm_control, laws$N, NULL)
      testthat::expect_identical(out$blocks[[b]]$branch, "cp04_separate_pair")
      testthat::expect_identical(out$blocks[[b]]$coefficients, old$coefficients)
      testthat::expect_identical(out$blocks[[b]]$H_chol, old$H_chol)
    }
  }
})

testthat::test_that("provider IRLS consumes bounded repeated scans and fitting scratch scales linearly in cells", {
  laws <- selection_build_32atom_laws()
  Z <- selection_ns("selection_model_design")(~ 1 + x, data.frame(x = laws$X))
  ds <- selection_ns("extract_selection_design_spec")(Z)
  cat <- selection_ns("partition_selection_cells")(
    laws$pair_specs, c(0L, 1L), selection_model = "independent_period_logit",
    selection_pooling = "across_cohorts"
  )
  visits <- matrix(0L, 6L, laws$N)
  largest <- 0L
  provider <- function(b, consume) {
    cell <- cat$cells[[b]]
    rr <- which(laws$K == cell$cohort - 1L)
    for (start in seq.int(1L, length(rr), by = 7L)) {
      rows <- rr[seq.int(start, min(length(rr), start + 6L))]
      largest <<- max(largest, length(rows))
      visits[b, rows] <<- visits[b, rows] + 1L
      w <- laws$date_laws$across_cohorts$w
      mass <- sum(w[rr]) / laws$N
      consume(list(D = matrix(laws$X[rows], ncol = 1L),
                   y = laws$R[rows, cell$calendar_index], w = w[rows],
                   a_k = sum(w[rr] * laws$A[rr]) / sum(w[rr]), mass_U_k = mass))
    }
  }
  fit <- selection_ns("fit_selection_block_provider")(
    provider, 6L, ds, selection_glm_control_tight, laws$N
  )
  testthat::expect_lt(max(abs(fit$coefficients - laws$date_laws$across_cohorts$theta)), 2e-11)
  testthat::expect_lte(largest, 7L)
  passes <- 2L * fit$diagnostics$iter + 3L
  for (b in 1:6) {
    testthat::expect_true(all(visits[b, laws$K == cat$cells[[b]]$cohort - 1L] == passes))
    testthat::expect_true(all(visits[b, laws$K != cat$cells[[b]]$cohort - 1L] == 0L))
  }
  testthat::expect_lt(unname(selection_ns("cp04_payload")(fit)["entries"]), 4096)
  # Resource algebra only: repeat catalog metadata, with the SAME 32 units and
  # fixed chunk size. No observations or scientific reference law is changed.
  st <- list(N = laws$N, cohort_start = c(1, 17), cohort_end = c(16, 32))
  ctrl <- list(max_chunk_rows = 7L) # scalar inputs to the pure bounds helper only
  bounds <- lapply(c(1L, 2L, 4L), function(mult) {
    cells <- rep(cat$cells, mult)
    specs <- list(list(B_m = length(cells), cell_ids = seq_along(cells)))
    selection_ns("selection_fit_resources")(
      st, specs, cells, ds, selection_glm_control, ctrl,
      "independent_period_logit", "across_cohorts"
    )
  })
  testthat::expect_equal(bounds[[3L]]$scratch - bounds[[2L]]$scratch,
                         2 * (bounds[[2L]]$scratch - bounds[[1L]]$scratch))
  testthat::expect_true(all(vapply(bounds, function(x) x$max_chunk_rows == 7 && x$legacy_rows == 0, logical(1))))
})

testthat::test_that("installed owners retain endpoint fields and match frozen full-F IF rows in every date and pair universe mode", {
  testthat::skip_if_not(selection_has_installed_cdid(), "Real installed cdid owners required.")
  laws <- selection_build_32atom_laws()
  root <- testthat::test_path("fixtures", "selection_models")
  testthat::expect_true(dir.exists(root))
  date_ref <- jsonlite::fromJSON(file.path(root, "independent_reference", "results.json"))$cases
  pair_ref <- jsonlite::fromJSON(file.path(root, "pair_pooling_reference", "results.json"))$cases
  sep_ref <- jsonlite::fromJSON(file.path(root, "pair_pooling_reference", "separate_date_results.json"))
  ctrl <- selection_ns("ipw_owner_control")(
    max_chunk_rows = 23L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7
  )
  ids <- seq_len(laws$N)
  # Embed exactly the frozen observation histories in real panel owners. The
  # calendar-valued outcome is an inert constructor field, never a fit input.
  g <- c(3, Inf)[laws$K + 1L]
  long <- data.frame(id = rep(ids, each = 3L), time = rep(1:3, laws$N),
                     y = ifelse(as.vector(t(laws$R)) == 1, rep(1:3, laws$N), NA_real_),
                     g = rep(g, each = 3L))
  panel <- selection_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  verify_case <- function(model, mode, u, w, coef_ref, frozen, Z, probes, U, rho) {
    frame <- selection_ns("construct_ipw_frame_owner")(
      panel, ids, cbind(1, laws$X), weights = w, control = ctrl, frame_g = g, allow_frame_only = TRUE
    )
    reader <- selection_ns("construct_cp04_frame_reader")(
      frame, ids, Z, g, w, ctrl, U = U, rho = rho, universe = u,
      justification = if (u == "F") NULL else "Known normalized ratio from the frozen 32-atom law"
    )
    views <- lapply(list(c(1, 2), c(1, 3), c(2, 3)), function(e) {
      selection_ns("construct_ipw_link_view")(
        frame, e[1L], e[2L], g = 3, K = Inf, control = ctrl,
        link_id = paste0("frozen", e[1L], e[2L]), functional_key = "cmpcat:frozen32"
      )
    })
    fit <- selection_ns("fit_selection_models")(
      reader, views, model, mode, ctrl, selection_glm_control_tight
    )
    testthat::expect_lt(max(abs(fit$coefficients - coef_ref)), 2e-11)
    pp <- selection_ns("predict_selection_probes")(fit, probes)
    if (model == "independent_period_logit") {
      # This is the actual field-flow regression for the host length-zero error.
      testthat::expect_true(all(vapply(pp$targets, function(tg) {
        length(tg$prob_s) == 1L && length(tg$prob_t) == 1L &&
          is.finite(tg$prob_s) && is.finite(tg$prob_t) &&
          abs(tg$prob_s * tg$prob_t - tg$pair_prob) < 1e-14
      }, logical(1))))
    }
    P <- length(coef_ref)
    Q <- length(pp$targets)
    testthat::expect_equal(ncol(frozen), P + 2L + 2L * Q)
    for (sz in c(7L, 23L)) {
      got <- matrix(NA_real_, laws$N, ncol(frozen))
      selection_ns("replay_selection_if")(
        fit, sz, callback = function(b) {
          got[b$rows, ] <<- cbind(b$theta_if, b$a_if, b$pair_prob_if, b$multiplier_if)
          testthat::expect_identical(b$id, ids[b$rows])
        }, probes = probes
      )
      testthat::expect_lt(max(abs(got - frozen)), 1e-12)
      # Log-weight IF is multiplier IF divided by the fixed evaluation multiplier.
      testthat::expect_lt(max(abs(sweep(got[, P + 2L + Q + seq_len(Q), drop = FALSE], 2L, pp$multipliers, "/") -
                                   sweep(frozen[, P + 2L + Q + seq_len(Q), drop = FALSE], 2L, pp$multipliers, "/"))), 1e-12)
    }
    cov <- selection_ns("stream_selection_covariance")(fit, 13L, probes = probes)
    testthat::expect_lt(max(abs(cov$covariance - crossprod(frozen) / laws$N^2)), 1e-12)

    # Independent forcing moments from FROZEN coefficients, retaining every pair
    # and both endpoints, contracted against the frozen original-unit IF arrays.
    C <- length(fit$cells)
    t_theta <- numeric(P)
    t_a <- numeric(2L)
    a <- vapply(0:1, function(k) sum(w[laws$K == k & U] * laws$A[laws$K == k & U]) /
                  sum(w[laws$K == k & U]), numeric(1))
    acc <- NULL
    for (q in seq_along(fit$pair_cohorts)) {
      pc <- fit$pair_cohorts[[q]]
      obs <- which(laws$K == pc$cohort - 1L &
                     laws$R[, pc$pair$s_calendar_index] * laws$R[, pc$pair$t_calendar_index] == 1)
      ow <- 0.75 * w[obs]
      acc <- selection_ns("add_selection_moments")(fit, q, obs, ow, acc = acc)
      t_a[pc$cohort] <- t_a[pc$cohort] + sum(ow / a[pc$cohort]) / laws$N
      cids <- if (model == "independent_period_logit") c(pc$s_cell_id, pc$t_cell_id) else pc$pair_cell_id
      for (cid in cids) {
        m <- fit$cells[[cid]]$block_id
        slope_idx <- if (ncol(Z) > 1L) C + m else integer(0)
        pr <- stats::plogis(coef_ref[cid] + if (length(slope_idx)) laws$X[obs] * coef_ref[slope_idx] else 0)
        force <- -ow * (1 - pr)
        t_theta[cid] <- t_theta[cid] + sum(force) / laws$N
        if (length(slope_idx)) t_theta[slope_idx] <- t_theta[slope_idx] + sum(force * laws$X[obs]) / laws$N
      }
    }
    actual_moments <- c(unlist(lapply(acc$blocks, `[[`, "t_alpha"), use.names = FALSE),
                        unlist(lapply(acc$blocks, `[[`, "t_beta"), use.names = FALSE))
    testthat::expect_lt(max(abs(actual_moments - t_theta)), 1e-13)
    testthat::expect_lt(max(abs(acc$t_a - t_a)), 1e-13)
    load <- numeric(laws$N)
    for (start in seq.int(1L, laws$N, by = 7L)) {
      rr <- seq.int(start, min(laws$N, start + 6L))
      load[rr] <- selection_ns("evaluate_pair_selection_loading")(fit, rr, acc)
    }
    expected <- drop(frozen[, seq_len(P), drop = FALSE] %*% t_theta) +
      drop(frozen[, P + 1:2, drop = FALSE] %*% t_a)
    testthat::expect_lt(max(abs(load - expected)), 1e-12)
  }
  for (key in c("within_cohort", "across_cohorts", "date_only", "separate_date_F")) {
    law <- laws$date_laws[[key]]
    mode <- if (key %in% c("date_only", "separate_date_F")) "separate" else key
    Z <- selection_ns("selection_model_design")(if (key == "date_only") ~ 1 else ~ 1 + x, data.frame(x = laws$X))
    probes <- list(cohort = c(3, 3, Inf, Inf),
                   Z = if (key == "date_only") matrix(numeric(0), 4L, 0L) else matrix(laws$probes_df$x, ncol = 1L))
    path <- if (key == "separate_date_F") file.path(root, "pair_pooling_reference", "separate_date_F_if.csv") else
      file.path(root, "independent_reference", paste0(key, "_if.csv"))
    frozen <- as.matrix(utils::read.csv(path, header = key != "separate_date_F"))
    coefs <- if (key == "separate_date_F") sep_ref$coefficients else date_ref[[key]]$coefficients
    verify_case("independent_period_logit", mode, "F", law$w, coefs, frozen, Z, probes,
                rep(TRUE, laws$N), rep(1, laws$N))
  }
  for (mode in c("separate", "within_cohort", "across_cohorts")) for (u in c("F", "B", "A")) {
    key <- paste0(mode, "_", u)
    Z <- selection_ns("selection_model_design")(~ 1 + x, data.frame(x = laws$X))
    rho <- laws$rho_list[[u]]
    rho_p <- vapply(1:4, function(j) rho[which(laws$K == laws$probes_df$cohort[j] &
                                              laws$X == laws$probes_df$x[j])[1L]], numeric(1))
    probes <- list(cohort = c(3, 3, Inf, Inf), Z = matrix(laws$probes_df$x, ncol = 1L), rho = rho_p)
    frozen <- as.matrix(utils::read.csv(file.path(root, "pair_pooling_reference", paste0(key, "_if.csv")), header = FALSE))
    verify_case("pair_logit", mode, u, laws$W_pair, pair_ref[[key]]$coefficients,
                frozen, Z, probes, laws$masks[[u]], rho)
  }
})
