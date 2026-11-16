# Comprehensive tests for R/static_complete_core.R against independent
# full-regression / estimating-equation and finite-difference references.

build_static_panel_fixture <- function(seed = 20261104L) {
  set.seed(seed)
  # Construct 60 original units across T = 4 dates with all 12 transition path
  # types represented (5 units per path type), covering:
  #   - initially untreated stayers (0,0,0,0)
  #   - single entrants (0,1,1,1), (0,0,1,1), (0,0,0,1)
  #   - entry + exit + reentry (0,1,0,1), (0,1,0,0), (0,0,1,0)
  #   - initially treated stayers (1,1,1,1)
  #   - single exiters (1,0,0,0), (1,1,0,0), (1,1,1,0)
  #   - exit + reentry + exit (1,0,1,0), (1,0,1,1), (1,1,0,1)
  path_types <- rbind(
    c(0L, 0L, 0L, 0L),
    c(0L, 1L, 1L, 1L),
    c(0L, 0L, 1L, 1L),
    c(0L, 0L, 0L, 1L),
    c(0L, 1L, 0L, 1L),
    c(0L, 1L, 0L, 0L),
    c(0L, 0L, 1L, 0L),
    c(1L, 1L, 1L, 1L),
    c(1L, 0L, 0L, 0L),
    c(1L, 1L, 0L, 0L),
    c(1L, 0L, 1L, 0L),
    c(1L, 0L, 1L, 1L)
  )
  n_rep <- 5L
  D_mat <- path_types[rep(seq_len(nrow(path_types)), each = n_rep), , drop = FALSE]
  n <- nrow(D_mat)
  unit_id <- sprintf("u_%03d", seq_len(n))

  # Reference-date covariates X_{i, t-1} and outcome levels Y_{i, t} with
  # nonlinear/misspecified components so all first-step IF terms are active.
  Z_list <- vector("list", 3L)
  Y_mat <- matrix(0.0, nrow = n, ncol = 4L)
  u_effect <- rnorm(n, sd = 0.4)
  Y_mat[, 1L] <- 1.0 + u_effect + 0.5 * D_mat[, 1L] + rnorm(n, sd = 0.3)

  for (step in 1:3) {
    t_cur <- step + 1L
    z1 <- rnorm(n, mean = 0.1 * step, sd = 0.8)
    z2 <- rnorm(n, mean = -0.1 * step, sd = 0.7)
    Z_step <- cbind(z1 = z1, z2 = z2)
    Z_list[[step]] <- Z_step

    # Outcome change depends on (t, D_{t-1}, D_t, Z_{t-1}) plus nonlinearity
    dy_step <- (0.3 * step) +
      (0.6 + 0.15 * step) * z1 -
      (0.4 - 0.10 * step) * z2 +
      0.25 * (z1^2 - 0.64) * (1 - 2 * D_mat[, step]) +
      (1.2 + 0.3 * step + 0.4 * z1 - 0.2 * z2) * (D_mat[, t_cur] - D_mat[, step]) +
      0.2 * u_effect +
      rnorm(n, sd = 0.25)
    Y_mat[, t_cur] <- Y_mat[, step] + dy_step
  }

  # Build all 6 declared cells c = (t, d) for t in {2, 3, 4}, d in {0, 1}
  cells <- list()
  for (step in 1:3) {
    t_cur <- step + 1L
    dy_t <- Y_mat[, t_cur] - Y_mat[, step]
    Z_t <- Z_list[[step]]
    X_t <- cbind("(Intercept)" = rep(1.0, n), Z_t)

    for (d_val in c(0L, 1L)) {
      cid <- sprintf("c_t%d_d%d", t_cur, d_val)
      R_c <- as.integer(D_mat[, step] == d_val)
      G_c <- as.integer(D_mat[, step] == d_val & D_mat[, t_cur] == (1L - d_val))
      C_c <- as.integer(D_mat[, step] == d_val & D_mat[, t_cur] == d_val)

      cells[[cid]] <- static_complete_cell(
        cell_id = cid,
        t = t_cur,
        d = d_val,
        sigma = 1 - 2 * d_val,
        unit_id = unit_id,
        R = R_c,
        G = G_c,
        C = C_c,
        dy = dy_t,
        X_tilde = X_t,
        Z = Z_t
      )
    }
  }

  cell_ids <- names(cells)
  partitions <- list(
    separate = stats::setNames(cell_ids, cell_ids),
    pooled = stats::setNames(rep("all_cells", length(cell_ids)), cell_ids),
    by_direction = stats::setNames(
      ifelse(vapply(cells, `[[`, integer(1L), "d") == 0L, "entry_block", "exit_block"),
      cell_ids
    )
  )

  cluster <- rep(sprintf("cl_%02d", seq_len(15L)), length.out = n)

  list(
    n = n,
    unit_id = unit_id,
    D_mat = D_mat,
    Y_mat = Y_mat,
    cells = cells,
    cell_ids = cell_ids,
    partitions = partitions,
    cluster = cluster,
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    rank_tol = 1e-10
  )
}

test_that("panel fixture covers entry, exit, repeated switches, and initially treated units", {
  fx <- build_static_panel_fixture()
  G_mat <- vapply(fx$cells, function(cl) cl$G, integer(fx$n))
  C_mat <- vapply(fx$cells, function(cl) cl$C, integer(fx$n))

  # Repeated switchers (units with >= 2 switches across the panel)
  switches_per_unit <- rowSums(G_mat)
  expect_true(any(switches_per_unit == 0L))
  expect_true(any(switches_per_unit == 1L))
  expect_true(any(switches_per_unit == 2L))
  expect_true(any(switches_per_unit == 3L))

  # Units that act as a switcher in one cell and a stayer in another cell
  expect_true(any(rowSums(G_mat) >= 1L & rowSums(C_mat) >= 1L))

  # Initially treated units (D_1 == 1) serve as switchers and stayers in exit cells
  # and as switchers/stayers in later entry cells after exiting
  init_treated <- which(fx$D_mat[, 1L] == 1L)
  expect_true(length(init_treated) > 0L)
  expect_true(sum(fx$cells[["c_t2_d1"]]$G[init_treated]) > 0L)
  expect_true(sum(fx$cells[["c_t2_d1"]]$C[init_treated]) > 0L)
  expect_true(sum(fx$cells[["c_t3_d0"]]$G[init_treated]) > 0L)
})

test_that("production kernel matches independent full-design reference across ipw/reg/dr and all OR partitions", {
  fx <- build_static_panel_fixture()
  sigma_signs <- vapply(fx$cells, `[[`, numeric(1L), "sigma")
  entry_cells <- c("c_t2_d0", "c_t3_d0", "c_t4_d0")
  exit_cells <- c("c_t2_d1", "c_t3_d1", "c_t4_d1")

  for (part_name in names(fx$partitions)) {
    or_blocks <- fx$partitions[[part_name]]
    for (m in c("ipw", "reg", "dr")) {
      prod_fit <- static_complete_fit(
        cells = fx$cells,
        or_blocks = or_blocks,
        method = m,
        glm_control = fx$glm_control,
        rank_tol = fx$rank_tol,
        aggregate_cells = fx$cell_ids,
        aggregate_signs = sigma_signs
      )
      ref_fit <- ref_static_complete_eval(
        cells = fx$cells,
        or_blocks = or_blocks,
        method = m,
        aggregate_cells = fx$cell_ids,
        aggregate_signs = sigma_signs,
        glm_control = fx$glm_control
      )

      expect_equal(prod_fit$theta_hat, ref_fit$theta_hat, tolerance = 1e-11)
      expect_equal(prod_fit$tau_hat, ref_fit$tau_hat, tolerance = 1e-11)
      expect_equal(prod_fit$pi_hat, ref_fit$pi_hat, tolerance = 1e-12)
      expect_equal(prod_fit$mu_1, ref_fit$mu_1, tolerance = 1e-11)
      if (m %in% c("ipw", "dr")) {
        expect_equal(prod_fit$Q_hat, ref_fit$Q_hat, tolerance = 1e-11)
        expect_equal(prod_fit$mu_0, ref_fit$mu_0, tolerance = 1e-11)
        for (cid in fx$cell_ids) {
          expect_equal(prod_fit$nuisance$logit[[cid]]$gamma_hat, ref_fit$gamma_hat[[cid]], tolerance = 1e-11)
          expect_equal(prod_fit$nuisance$logit[[cid]]$H_hat, ref_fit$H_hat[[cid]], tolerance = 1e-11)
          expect_equal(prod_fit$nuisance$logit[[cid]]$ell_gamma, ref_fit$ell_gamma[[cid]], tolerance = 1e-10)
        }
      }
      if (m %in% c("reg", "dr")) {
        expect_equal(prod_fit$nuisance$or$alpha_hat, ref_fit$alpha_hat, tolerance = 1e-11)
        for (bid in names(ref_fit$slope_hat)) {
          expect_equal(unname(prod_fit$nuisance$or$slope_hat[[bid]]), unname(ref_fit$slope_hat[[bid]]), tolerance = 1e-11)
        }
      }

      # Original-unit cell influence functions and cross-cell covariances
      expect_equal(prod_fit$phi_theta, ref_fit$phi_theta, tolerance = 1e-10)
      expect_equal(prod_fit$phi_tau, ref_fit$phi_tau, tolerance = 1e-10)
      expect_equal(prod_fit$cov_theta, ref_fit$cov_theta, tolerance = 1e-10)
      expect_equal(prod_fit$cov_tau, ref_fit$cov_tau, tolerance = 1e-10)
      expect_equal(unname(colMeans(prod_fit$phi_theta)), rep(0.0, length(fx$cell_ids)), tolerance = 1e-11)

      # Full source aggregate delta = delta_{C, sigma}
      expect_equal(prod_fit$aggregate$delta_hat, ref_fit$aggregate$delta_hat, tolerance = 1e-11)
      expect_equal(prod_fit$aggregate$omega_hat, ref_fit$aggregate$omega_hat, tolerance = 1e-12)
      expect_equal(prod_fit$aggregate$phi, ref_fit$aggregate$phi, tolerance = 1e-10)
      expect_equal(
        prod_fit$aggregate$phi_components$cell_if_term,
        ref_fit$aggregate$cell_if_term,
        tolerance = 1e-10
      )
      expect_equal(
        prod_fit$aggregate$phi_components$share_if_term,
        ref_fit$aggregate$share_if_term,
        tolerance = 1e-10
      )
      expect_equal(prod_fit$aggregate$variance, ref_fit$aggregate$variance, tolerance = 1e-10)
      expect_equal(prod_fit$aggregate$vhat, ref_fit$aggregate$vhat, tolerance = 1e-10)
      expect_equal(prod_fit$aggregate$se, ref_fit$aggregate$se, tolerance = 1e-10)
      expect_equal(mean(prod_fit$aggregate$phi), 0.0, tolerance = 1e-11)

      # Subaggregates: entry-only, exit natural (s = +1), exit sign-oriented (s = -1)
      agg_entry <- static_complete_aggregate(prod_fit, cells = entry_cells, signs = c(1, 1, 1))
      ref_entry <- ref_static_complete_eval(
        fx$cells, or_blocks, m, aggregate_cells = entry_cells, aggregate_signs = c(1, 1, 1)
      )$aggregate
      expect_equal(agg_entry$delta_hat, ref_entry$delta_hat, tolerance = 1e-11)
      expect_equal(agg_entry$phi, ref_entry$phi, tolerance = 1e-10)
      expect_equal(agg_entry$variance, ref_entry$variance, tolerance = 1e-10)

      agg_exit_nat <- static_complete_aggregate(prod_fit, cells = exit_cells, signs = c(1, 1, 1))
      ref_exit_nat <- ref_static_complete_eval(
        fx$cells, or_blocks, m, aggregate_cells = exit_cells, aggregate_signs = c(1, 1, 1)
      )$aggregate
      expect_equal(agg_exit_nat$delta_hat, ref_exit_nat$delta_hat, tolerance = 1e-11)
      expect_equal(agg_exit_nat$phi, ref_exit_nat$phi, tolerance = 1e-10)

      agg_exit_sig <- static_complete_aggregate(prod_fit, cells = exit_cells, signs = c(-1, -1, -1))
      expect_equal(agg_exit_sig$delta_hat, -agg_exit_nat$delta_hat, tolerance = 1e-12)
      expect_equal(agg_exit_sig$phi, -agg_exit_nat$phi, tolerance = 1e-12)
      expect_equal(agg_exit_sig$variance, agg_exit_nat$variance, tolerance = 1e-12)
    }
  }

  # Paper method alias "or" and list-partition representation produce identical results
  fit_reg <- static_complete_fit(
    fx$cells,
    or_blocks = list(entry = entry_cells, exit = exit_cells),
    method = "or",
    glm_control = list(epsilon = 1e-12, maxit = 50L),
    rank_tol = fx$rank_tol
  )
  fit_vec <- static_complete_fit(
    fx$cells,
    or_blocks = fx$partitions$by_direction,
    method = "reg",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol
  )
  expect_equal(fit_reg$theta_hat, fit_vec$theta_hat, tolerance = 1e-14)
  expect_equal(fit_reg$phi_theta, fit_vec$phi_theta, tolerance = 1e-14)
})

test_that("analytic IF matches n times centered unit-weight finite-difference directional derivative", {
  fx <- build_static_panel_fixture()
  sigma_signs <- vapply(fx$cells, `[[`, numeric(1L), "sigma")

  for (part_name in c("separate", "by_direction", "pooled")) {
    or_blocks <- fx$partitions[[part_name]]
    for (m in c("ipw", "reg", "dr")) {
      prod_fit <- static_complete_fit(
        cells = fx$cells,
        or_blocks = or_blocks,
        method = m,
        glm_control = fx$glm_control,
        rank_tol = fx$rank_tol,
        aggregate_cells = fx$cell_ids,
        aggregate_signs = sigma_signs
      )
      fd <- ref_static_complete_fd_if(
        cells = fx$cells,
        or_blocks = or_blocks,
        method = m,
        aggregate_cells = fx$cell_ids,
        aggregate_signs = sigma_signs,
        eps = 1e-6,
        glm_control = fx$glm_control
      )

      expect_equal(prod_fit$phi_theta, fd$phi_theta, tolerance = 1e-7)
      expect_equal(prod_fit$phi_tau, fd$phi_tau, tolerance = 1e-7)
      expect_equal(prod_fit$aggregate$phi, fd$phi_agg, tolerance = 1e-7)
    }
  }

  # Also test differing treatment vs outcome design matrices (`X_tilde` vs `Z_or`)
  cells_diff_design <- lapply(fx$cells, function(cl) {
    z_ipw <- cbind(cl$Z, z_inter = cl$Z[, 1L] * cl$Z[, 2L])
    z_or <- cl$Z[, 1L, drop = FALSE]
    static_complete_cell(
      cell_id = cl$cell_id,
      t = cl$t,
      d = cl$d,
      sigma = cl$sigma,
      unit_id = cl$unit_id,
      R = cl$R,
      G = cl$G,
      C = cl$C,
      dy = cl$dy,
      Z = z_ipw,
      Z_or = z_or,
      or_block = fx$partitions$by_direction[[cl$cell_id]]
    )
  })
  fit_diff <- static_complete_fit(
    cells = cells_diff_design,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    aggregate_cells = c("c_t2_d0", "c_t2_d1", "c_t3_d0"),
    aggregate_signs = c(1, -1, 1)
  )
  ref_diff <- ref_static_complete_eval(
    cells = cells_diff_design,
    or_blocks = fx$partitions$by_direction,
    method = "dr",
    aggregate_cells = c("c_t2_d0", "c_t2_d1", "c_t3_d0"),
    aggregate_signs = c(1, -1, 1),
    glm_control = fx$glm_control
  )
  fd_diff <- ref_static_complete_fd_if(
    cells = cells_diff_design,
    or_blocks = fx$partitions$by_direction,
    method = "dr",
    aggregate_cells = c("c_t2_d0", "c_t2_d1", "c_t3_d0"),
    aggregate_signs = c(1, -1, 1),
    eps = 1e-6,
    glm_control = fx$glm_control
  )
  expect_equal(fit_diff$theta_hat, ref_diff$theta_hat, tolerance = 1e-11)
  expect_equal(fit_diff$phi_theta, ref_diff$phi_theta, tolerance = 1e-10)
  expect_equal(fit_diff$phi_theta, fd_diff$phi_theta, tolerance = 1e-7)
  expect_equal(fit_diff$aggregate$phi, fd_diff$phi_agg, tolerance = 1e-7)

  # Test subaggregate A < C under fully pooled OR partition P = {C}, where
  # stayers in cells c' in C \ A contribute to phi_{A,s} through shared slopes,
  # plus mixed intercept-only / covariate designs and multi-unit directional FD.
  entry_cells <- c("c_t2_d0", "c_t3_d0", "c_t4_d0")
  fit_sub_pool <- static_complete_fit(
    cells = fx$cells,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    aggregate_cells = entry_cells,
    aggregate_signs = c(1, 1, 1)
  )
  fd_sub_pool <- ref_static_complete_fd_if(
    cells = fx$cells,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    aggregate_cells = entry_cells,
    aggregate_signs = c(1, 1, 1),
    eps = 1e-6,
    glm_control = fx$glm_control
  )
  expect_equal(fit_sub_pool$aggregate$phi, fd_sub_pool$phi_agg, tolerance = 1e-7)

  # Multi-unit centered directional derivative <d, phi> == n * (T(1 + eps*d) - T(1 - eps*d)) / (2*eps)
  set.seed(20261106L)
  dir_vec <- rnorm(fx$n)
  dir_vec <- dir_vec - mean(dir_vec)
  eps_dir <- 1e-6
  plus_dir <- ref_static_complete_eval(
    fx$cells, fx$partitions$pooled, "dr",
    unit_weights = 1.0 + eps_dir * dir_vec,
    aggregate_cells = entry_cells, aggregate_signs = c(1, 1, 1),
    glm_control = fx$glm_control, compute_if = FALSE, polish = TRUE
  )
  minus_dir <- ref_static_complete_eval(
    fx$cells, fx$partitions$pooled, "dr",
    unit_weights = 1.0 - eps_dir * dir_vec,
    aggregate_cells = entry_cells, aggregate_signs = c(1, 1, 1),
    glm_control = fx$glm_control, compute_if = FALSE, polish = TRUE
  )
  fd_dir_agg <- fx$n * (plus_dir$aggregate$delta_hat - minus_dir$aggregate$delta_hat) / (2.0 * eps_dir)
  expect_equal(sum(dir_vec * fit_sub_pool$aggregate$phi), fd_dir_agg, tolerance = 1e-7)
})

test_that("formula-perturbation witnesses and Proposition 5(b) normalization gap hold", {
  fx <- build_static_panel_fixture()
  sigma_signs <- vapply(fx$cells, `[[`, numeric(1L), "sigma")

  fit_dr <- static_complete_fit(
    cells = fx$cells,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    aggregate_cells = fx$cell_ids,
    aggregate_signs = sigma_signs
  )
  fd_dr <- ref_static_complete_fd_if(
    cells = fx$cells,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    aggregate_cells = fx$cell_ids,
    aggregate_signs = sigma_signs,
    eps = 1e-6,
    glm_control = fx$glm_control
  )

  # Witness 1: Omitting logit first-step score term `-k_{gamma,c}' ell^gamma_c`
  # under misspecified OR breaks derivative agreement by the exact logit_score component.
  for (cid in fx$cell_ids) {
    comp <- fit_dr$cell_components[[cid]]
    phi_no_logit <- comp$direct_switcher + comp$direct_stayer + comp$or_score
    expect_gt(max(abs(comp$logit_score)), 1e-2)
    expect_gt(max(abs(phi_no_logit - fd_dr$phi_theta[, cid])), 1e-2)
    expect_equal(
      (fit_dr$phi_theta[, cid] - phi_no_logit),
      comp$logit_score,
      tolerance = 1e-14
    )
  }

  # Witness 2: Omitting OR first-step score term `-k_{beta,c}' ell^beta` under
  # misspecified logit breaks derivative agreement, and stayers from other cells
  # in the pooled block contribute non-zero IF values to cell c.
  for (cid in fx$cell_ids) {
    comp <- fit_dr$cell_components[[cid]]
    phi_no_or <- comp$direct_switcher + comp$direct_stayer + comp$logit_score
    expect_gt(max(abs(comp$or_score)), 1e-2)
    expect_gt(max(abs(phi_no_or - fd_dr$phi_theta[, cid])), 1e-2)

    # Units not at risk in cell `cid` (`R_{ic} == 0`) that are stayers in other
    # cells of the pooled block have non-zero `phi_theta[i, cid]` solely from `or_score`!
    outside_risk <- which(fx$cells[[cid]]$R == 0L)
    expect_true(length(outside_risk) > 0L)
    expect_equal(fit_dr$phi_theta[outside_risk, cid], comp$or_score[outside_risk], tolerance = 1e-14)
    expect_gt(max(abs(fit_dr$phi_theta[outside_risk, cid])), 1e-3)
  }

  # Witness 3: Omitting the random switch-share term in the aggregate IF breaks
  # agreement with the finite-difference derivative of delta_hat_{A,s}.
  cell_only_if <- fit_dr$aggregate$phi_components$cell_if_term
  share_if <- fit_dr$aggregate$phi_components$share_if_term
  expect_gt(max(abs(share_if)), 1e-2)
  expect_gt(max(abs(cell_only_if - fd_dr$phi_agg)), 1e-2)
  expect_equal(cell_only_if + share_if, fd_dr$phi_agg, tolerance = 1e-7)

  # Witness 4: Exact difference between Hajek-normalized `delta_hat^dr` and
  # Horvitz-Thompson-normalized `delta_hat^src` (de Chaisemartin et al. 2026).
  n_sw <- fit_dr$n_switchers
  total_sw <- sum(n_sw)
  delta_src_terms <- numeric(length(fx$cell_ids))
  prop5b_gap_terms <- numeric(length(fx$cell_ids))
  names(delta_src_terms) <- names(prop5b_gap_terms) <- fx$cell_ids

  for (cid in fx$cell_ids) {
    cl <- fx$cells[[cid]]
    r_c <- fit_dr$nuisance$logit[[cid]]$odds
    e_c <- fit_dr$nuisance$or$residuals[, cid]
    n_c <- as.double(cl$n_switchers)
    sum_Cr <- sum(cl$C * r_c)
    sum_Cre <- sum(cl$C * r_c * e_c)
    sum_Ge <- sum(cl$G * e_c)

    delta_src_terms[[cid]] <- (n_c / total_sw) * cl$sigma * ((sum_Ge - sum_Cre) / n_c)
    prop5b_gap_terms[[cid]] <- (n_c / total_sw) * cl$sigma * sum_Cre * (1.0 / n_c - 1.0 / sum_Cr)
  }
  delta_src <- sum(delta_src_terms)
  prop5b_gap <- sum(prop5b_gap_terms)

  expect_gt(abs(prop5b_gap), 1e-4)
  # Preserve the scalar absolute-or-relative bound across testthat editions.
  prop5b_error <- abs((fit_dr$aggregate$delta_hat - delta_src) - prop5b_gap) /
    max(1, abs(prop5b_gap))
  expect_lt(prop5b_error, 1e-14)
})

test_that("Proposition 5(c) intercept-only equality holds across ipw/reg/dr and all OR partitions", {
  fx <- build_static_panel_fixture()
  sigma_signs <- vapply(fx$cells, `[[`, numeric(1L), "sigma")

  cells_int <- lapply(fx$cells, function(cl) {
    static_complete_cell(
      cell_id = cl$cell_id,
      t = cl$t,
      d = cl$d,
      sigma = cl$sigma,
      unit_id = cl$unit_id,
      R = cl$R,
      G = cl$G,
      C = cl$C,
      dy = cl$dy,
      X_tilde = matrix(1.0, nrow = fx$n, ncol = 1L),
      Z = matrix(numeric(0L), nrow = fx$n, ncol = 0L)
    )
  })

  # Direct difference in mean changes between switchers and stayers in each cell
  diff_means <- vapply(cells_int, function(cl) {
    mean(cl$dy[cl$G == 1L]) - mean(cl$dy[cl$C == 1L])
  }, numeric(1L))
  n_sw <- vapply(cells_int, `[[`, integer(1L), "n_switchers")
  delta_src_int <- sum((n_sw / sum(n_sw)) * sigma_signs * diff_means)

  ref_fit_int <- NULL
  for (part_name in names(fx$partitions)) {
    or_blocks <- fx$partitions[[part_name]]
    for (m in c("ipw", "reg", "dr")) {
      fit_m <- static_complete_fit(
        cells = cells_int,
        or_blocks = or_blocks,
        method = m,
        glm_control = fx$glm_control,
        rank_tol = fx$rank_tol,
        aggregate_cells = fx$cell_ids,
        aggregate_signs = sigma_signs
      )
      expect_equal(fit_m$theta_hat, diff_means, tolerance = 1e-12)
      expect_equal(fit_m$tau_hat, sigma_signs * diff_means, tolerance = 1e-12)
      expect_equal(fit_m$aggregate$delta_hat, delta_src_int, tolerance = 1e-12)

      if (is.null(ref_fit_int)) {
        ref_fit_int <- fit_m
      } else {
        expect_equal(fit_m$phi_theta, ref_fit_int$phi_theta, tolerance = 1e-12)
        expect_equal(fit_m$phi_tau, ref_fit_int$phi_tau, tolerance = 1e-12)
        expect_equal(fit_m$aggregate$phi, ref_fit_int$aggregate$phi, tolerance = 1e-12)
        expect_equal(fit_m$cov_theta, ref_fit_int$cov_theta, tolerance = 1e-12)
        expect_equal(fit_m$aggregate$variance, ref_fit_int$aggregate$variance, tolerance = 1e-12)
      }
    }
  }
})

test_that("probability ceiling 1 - 1e-6 caps odds and zeroes clipped-odds derivative above cap", {
  # Symmetric two-covariate risk-set design around (z1 = 0, z2 = 0, G = 0.5) so
  # gamma_0 = 0, gamma_1 = gamma_2 ~ 1.901, and the stayer at (z1, z2) = (4.2, 4.2)
  # has linear predictor u_hat ~ 15.97 > kbar (~ 13.8155), giving
  # p_hat ~ 1 - 1.16e-7 in (1 - 1e-6, 1 - 10 * .Machine$double.eps).
  Z_cap <- rbind(
    matrix(rep(c(1.0, 0.0), each = 40L), ncol = 2L),
    matrix(rep(c(0.0, 1.0), each = 40L), ncol = 2L),
    matrix(rep(c(-1.0, 0.0), each = 40L), ncol = 2L),
    matrix(rep(c(0.0, -1.0), each = 40L), ncol = 2L),
    c(4.2, 4.2),
    c(-4.2, -4.2)
  )
  colnames(Z_cap) <- c("z1", "z2")
  G <- c(
    rep(1L, 39L), 0L,
    rep(1L, 39L), 0L,
    rep(0L, 39L), 1L,
    rep(0L, 39L), 1L,
    0L,
    1L
  )
  n <- nrow(Z_cap)
  unit_id <- sprintf("u_%03d", seq_len(n))
  C <- 1L - G
  R <- rep(1L, n)
  set.seed(20261105L)
  dy <- 1.5 + 0.8 * Z_cap[, 1L] - 0.5 * Z_cap[, 2L] +
    0.4 * (Z_cap[, 1L]^2 + Z_cap[, 2L]^2) + 2.0 * G + rnorm(n, sd = 0.2)

  cell_cap <- static_complete_cell(
    cell_id = "c_cap",
    t = 2,
    d = 0L,
    unit_id = unit_id,
    R = R,
    G = G,
    C = C,
    dy = dy,
    Z = Z_cap
  )

  ctrl <- list(epsilon = 1e-12, maxit = 50L, trace = FALSE)
  for (m in c("ipw", "dr")) {
    fit_cap <- static_complete_fit(
      cells = list(c_cap = cell_cap),
      or_blocks = c(c_cap = "b1"),
      method = m,
      glm_control = ctrl,
      rank_tol = 1e-10,
      aggregate_cells = "c_cap",
      aggregate_signs = 1
    )
    lfit <- fit_cap$nuisance$logit[["c_cap"]]
    clip_diag <- lfit$diagnostics$clipping

    # Verify ceiling is active on the (4.2, 4.2) stayer unit
    expect_equal(clip_diag$p_capped, 1)
    expect_equal(clip_diag$p_capped_ctrl, 1)
    capped_idx <- which(Z_cap[, 1L] == 4.2)
    uncapped_idx <- which(Z_cap[, 1L] < 4.2)
    expect_true(all(lfit$prob[capped_idx] > ipw_clip_cap))
    expect_true(all(lfit$prob[uncapped_idx] < ipw_clip_cap))
    expect_equal(lfit$active[capped_idx], rep(0.0, length(capped_idx)))
    expect_equal(lfit$active[uncapped_idx], rep(1.0, length(uncapped_idx)))
    expect_equal(lfit$odds[capped_idx], rep(ipw_clip_cap_odds, length(capped_idx)), tolerance = 1e-12)

    # Compare analytic IF against reference and centered unit-weight finite-difference derivative
    ref_cap <- ref_static_complete_eval(
      cells = list(c_cap = cell_cap),
      or_blocks = c(c_cap = "b1"),
      method = m,
      aggregate_cells = "c_cap",
      aggregate_signs = 1,
      glm_control = ctrl
    )
    ref_cap_polished <- ref_static_complete_eval(
      cells = list(c_cap = cell_cap),
      or_blocks = c(c_cap = "b1"),
      method = m,
      aggregate_cells = "c_cap",
      aggregate_signs = 1,
      glm_control = ctrl,
      polish = TRUE
    )
    fd_cap <- ref_static_complete_fd_if(
      cells = list(c_cap = cell_cap),
      or_blocks = c(c_cap = "b1"),
      method = m,
      aggregate_cells = "c_cap",
      aggregate_signs = 1,
      eps = 1e-5,
      glm_control = ctrl
    )
    expect_equal(fit_cap$phi_theta[, "c_cap"], ref_cap$phi_theta[, "c_cap"], tolerance = 1e-12)
    expect_equal(ref_cap_polished$phi_theta[, "c_cap"], fd_cap$phi_theta[, "c_cap"], tolerance = 1e-7)
    expect_equal(fit_cap$phi_theta[, "c_cap"], fd_cap$phi_theta[, "c_cap"], tolerance = 1e-7)

    # Witness: If the zero-derivative indicator `1{u < kbar}` were omitted (treating
    # capped units as having derivative `r_hat * xtilde`), the IF would disagree with FD!
    res_c <- if (m == "ipw") {
      dy - fit_cap$mu_0[["c_cap"]]
    } else {
      fit_cap$nuisance$or$residuals[, "c_cap"] - fit_cap$mu_0[["c_cap"]]
    }
    k_gamma_wrong <- colSums(cell_cap$X_tilde * (C * lfit$odds * res_c)) / lfit$sum_stayer_odds
    lambda_wrong <- static_complete_logit_hessian_solve(lfit, k_gamma_wrong)
    logit_score_wrong <- -(R * lfit$score_weight) * drop(cell_cap$X_tilde %*% lambda_wrong)
    comp <- fit_cap$cell_components[["c_cap"]]
    phi_wrong <- comp$direct_switcher + comp$direct_stayer + logit_score_wrong + comp$or_score
    expect_gt(max(abs(phi_wrong - fd_cap$phi_theta[, "c_cap"])), 1e-2)
  }
})

test_that("Corollary 1 branches and pooled OR block identification with cell-constant covariate hold", {
  # Part A: Exact linear outcome model for block P(c) -> mu_0c = 0, k_{gamma,c} = 0,
  # and both OR and DR recover exact theta_c even under misspecified logit or active ceiling.
  fx <- build_static_panel_fixture()
  true_theta <- c(c_t2_d0 = 2.5, c_t2_d1 = -1.8, c_t3_d0 = 3.0, c_t3_d1 = -2.2, c_t4_d0 = 1.7, c_t4_d1 = -1.1)
  true_alpha <- c(c_t2_d0 = 0.4, c_t2_d1 = -0.3, c_t3_d0 = 0.7, c_t3_d1 = -0.5, c_t4_d0 = 0.2, c_t4_d1 = 0.1)
  true_b_pool <- c(z1 = 1.25, z2 = -0.75)

  cells_exact_or <- lapply(fx$cell_ids, function(cid) {
    cl <- fx$cells[[cid]]
    dy_exact <- true_alpha[[cid]] + drop(cl$Z %*% true_b_pool) + true_theta[[cid]] * cl$G
    static_complete_cell(
      cell_id = cid,
      t = cl$t,
      d = cl$d,
      sigma = cl$sigma,
      unit_id = cl$unit_id,
      R = cl$R,
      G = cl$G,
      C = cl$C,
      dy = dy_exact,
      Z = cl$Z
    )
  })
  names(cells_exact_or) <- fx$cell_ids

  fit_dr_exact <- static_complete_fit(
    cells = cells_exact_or,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol
  )
  fit_or_exact <- static_complete_fit(
    cells = cells_exact_or,
    or_blocks = fx$partitions$pooled,
    method = "reg",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol
  )
  expect_equal(fit_dr_exact$theta_hat, true_theta, tolerance = 1e-12)
  expect_equal(fit_or_exact$theta_hat, true_theta, tolerance = 1e-12)
  expect_equal(unname(fit_dr_exact$mu_0), rep(0.0, length(fx$cell_ids)), tolerance = 1e-12)
  for (cid in fx$cell_ids) {
    expect_equal(unname(fit_dr_exact$cell_components[[cid]]$k_gamma), rep(0.0, 3L), tolerance = 1e-12)
    expect_equal(fit_dr_exact$cell_components[[cid]]$logit_score, rep(0.0, fx$n), tolerance = 1e-12)
  }

  # Part B: Assumption 4(c) pooled block identification when z1 is constant among
  # stayers of cell 1 (`c_t2_d0`) but varies among stayers of cell 2 (`c_t3_d0`).
  c1 <- fx$cells[["c_t2_d0"]]
  c2 <- fx$cells[["c_t3_d0"]]
  Z_c1_const <- c1$Z
  Z_c1_const[c1$C == 1L, 1L] <- 3.0 # Constant on stayers of c1, varies on switchers of c1 and stayers of c2

  cell_1_const <- static_complete_cell(
    cell_id = "c_t2_d0", t = c1$t, d = c1$d, unit_id = c1$unit_id,
    R = c1$R, G = c1$G, C = c1$C, dy = c1$dy,
    Z = c1$Z, Z_or = Z_c1_const
  )
  cell_2_var <- static_complete_cell(
    cell_id = "c_t3_d0", t = c2$t, d = c2$d, unit_id = c2$unit_id,
    R = c2$R, G = c2$G, C = c2$C, dy = c2$dy,
    Z = c2$Z, Z_or = c2$Z
  )
  two_cells <- list(c_t2_d0 = cell_1_const, c_t3_d0 = cell_2_var)

  # Separate OR partition must fail rank check on cell 1
  expect_error(
    static_complete_fit(
      cells = two_cells,
      or_blocks = c(c_t2_d0 = "b1", c_t3_d0 = "b2"),
      method = "dr",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_or_wls_rank_error"
  )

  # Pooled OR partition identifies the common slope from cell 2 and matches full-design reference & FD!
  fit_pooled_id <- static_complete_fit(
    cells = two_cells,
    or_blocks = c(c_t2_d0 = "shared", c_t3_d0 = "shared"),
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    aggregate_cells = c("c_t2_d0", "c_t3_d0"),
    aggregate_signs = c(1, 1)
  )
  ref_pooled_id <- ref_static_complete_eval(
    cells = two_cells,
    or_blocks = c(c_t2_d0 = "shared", c_t3_d0 = "shared"),
    method = "dr",
    aggregate_cells = c("c_t2_d0", "c_t3_d0"),
    aggregate_signs = c(1, 1),
    glm_control = fx$glm_control
  )
  fd_pooled_id <- ref_static_complete_fd_if(
    cells = two_cells,
    or_blocks = c(c_t2_d0 = "shared", c_t3_d0 = "shared"),
    method = "dr",
    aggregate_cells = c("c_t2_d0", "c_t3_d0"),
    aggregate_signs = c(1, 1),
    eps = 1e-6,
    glm_control = fx$glm_control
  )
  expect_equal(fit_pooled_id$theta_hat, ref_pooled_id$theta_hat, tolerance = 1e-11)
  expect_equal(fit_pooled_id$phi_theta, ref_pooled_id$phi_theta, tolerance = 1e-10)
  expect_equal(fit_pooled_id$phi_theta, fd_pooled_id$phi_theta, tolerance = 1e-7)
  expect_equal(fit_pooled_id$aggregate$phi, fd_pooled_id$phi_agg, tolerance = 1e-7)
})

test_that("original-unit covariance, cross-cell covariance, and inference-only clusters work properly", {
  fx <- build_static_panel_fixture()
  sigma_signs <- vapply(fx$cells, `[[`, numeric(1L), "sigma")

  fit_uncl <- static_complete_fit(
    cells = fx$cells,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    cluster = NULL,
    aggregate_cells = fx$cell_ids,
    aggregate_signs = sigma_signs
  )
  fit_cl <- static_complete_fit(
    cells = fx$cells,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    cluster = fx$cluster,
    aggregate_cells = fx$cell_ids,
    aggregate_signs = sigma_signs
  )

  # Cluster choice must leave all point estimates, nuisance fits, and unit-level IFs unchanged
  expect_identical(fit_cl$theta_hat, fit_uncl$theta_hat)
  expect_identical(fit_cl$tau_hat, fit_uncl$tau_hat)
  expect_identical(fit_cl$phi_theta, fit_uncl$phi_theta)
  expect_identical(fit_cl$phi_tau, fit_uncl$phi_tau)
  expect_identical(fit_cl$aggregate$delta_hat, fit_uncl$aggregate$delta_hat)
  expect_identical(fit_cl$aggregate$phi, fit_uncl$aggregate$phi)

  # Cross-cell covariance has non-zero off-diagonal entries
  off_diag <- fit_uncl$cov_theta[row(fit_uncl$cov_theta) != col(fit_uncl$cov_theta)]
  expect_gt(max(abs(off_diag)), 1e-4)

  # Original-unit variance of cell-IF component includes cross-cell covariance,
  # so it differs from the naive sum of cell variances ignoring shared original units:
  omega <- fit_uncl$aggregate$omega_hat
  w_signed <- omega * sigma_signs
  var_cell_term_true <- as.double(t(w_signed) %*% fit_uncl$cov_theta %*% w_signed)
  var_cell_term_naive <- sum((w_signed^2) * diag(fit_uncl$cov_theta))
  expect_equal(
    sum(fit_uncl$aggregate$phi_components$cell_if_term^2) / (fx$n^2),
    var_cell_term_true,
    tolerance = 1e-14
  )
  expect_gt(abs(var_cell_term_true - var_cell_term_naive), 1e-5)

  # Cluster covariance matches Corollary 2 (`crossprod(rowsum_by_cluster(IF)) / n^2`)
  # with no silent degrees-of-freedom or few-cluster multiplier
  ref_cl <- ref_static_complete_eval(
    cells = fx$cells,
    or_blocks = fx$partitions$pooled,
    method = "dr",
    aggregate_cells = fx$cell_ids,
    aggregate_signs = sigma_signs,
    cluster = fx$cluster,
    glm_control = fx$glm_control
  )
  expect_equal(fit_cl$cov_theta, ref_cl$cov_theta, tolerance = 1e-12)
  expect_equal(fit_cl$cov_tau, ref_cl$cov_tau, tolerance = 1e-12)
  expect_equal(fit_cl$aggregate$variance, ref_cl$aggregate$variance, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(fit_cl$cov_theta, fit_uncl$cov_theta)))

  # Singleton clusters (each unit in its own cluster) match unclustered covariance identically
  inf_singleton <- static_complete_inference(fit_uncl, cluster = fx$unit_id)
  expect_equal(inf_singleton$cov_theta, fit_uncl$cov_theta, tolerance = 1e-15)
  expect_equal(inf_singleton$aggregate$variance, fit_uncl$aggregate$variance, tolerance = 1e-15)
})

test_that("exact alias, rank, support, convergence, separation, boundary, and input refusals are enforced", {
  fx <- build_static_panel_fixture()
  base_cell <- fx$cells[["c_t2_d0"]]

  # 1. No-switcher refusal (n_c = 0) and no-stayer refusal (n_{0c} = 0)
  expect_error(
    static_complete_cell(
      cell_id = "no_sw", t = 2, d = 0L, unit_id = fx$unit_id,
      R = base_cell$C, G = rep(0L, fx$n), C = base_cell$C,
      dy = base_cell$dy, Z = base_cell$Z
    ),
    class = "cdid_static_support_error"
  )
  err_no_stayer <- tryCatch(
    static_complete_cell(
      cell_id = "no_st", t = 2, d = 0L, unit_id = fx$unit_id,
      R = base_cell$G, G = base_cell$G, C = rep(0L, fx$n),
      dy = base_cell$dy, Z = base_cell$Z
    ),
    error = function(e) e
  )
  expect_s3_class(err_no_stayer, "cdid_static_support_error")
  expect_s3_class(err_no_stayer, "cdid_ipw_support_error")
  expect_true(isTRUE(err_no_stayer$sample_support_only))

  # 2. Mismatched original-unit identity or order across cells
  permuted_cell <- static_complete_cell(
    cell_id = "c_perm", t = 3, d = 0L,
    unit_id = rev(fx$unit_id),
    R = rev(base_cell$R), G = rev(base_cell$G), C = rev(base_cell$C),
    dy = rev(base_cell$dy), Z = base_cell$Z[rev(seq_len(fx$n)), , drop = FALSE]
  )
  expect_error(
    static_complete_fit(
      cells = list(c_t2_d0 = base_cell, c_perm = permuted_cell),
      or_blocks = c(c_t2_d0 = "b1", c_perm = "b2"),
      method = "dr",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_static_unit_id_error"
  )

  # 3. Mask inconsistency (R_c != G_c + C_c or G_c * C_c != 0)
  expect_error(
    static_complete_cell(
      cell_id = "bad_mask", t = 2, d = 0L, unit_id = fx$unit_id,
      R = rep(1L, fx$n), G = base_cell$G, C = base_cell$C,
      dy = base_cell$dy, Z = base_cell$Z
    ),
    class = "cdid_static_mask_error"
  )

  # 4. Missing required explicit controls / partition / method (no silent defaults)
  expect_error(
    static_complete_fit(
      cells = fx$cells,
      or_blocks = NULL,
      method = "dr",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_static_partition_error"
  )
  expect_error(
    static_complete_fit(
      cells = fx$cells,
      or_blocks = fx$partitions$separate,
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_static_method_error"
  )
  expect_error(
    static_complete_fit(
      cells = fx$cells,
      or_blocks = fx$partitions$separate,
      method = "dr",
      rank_tol = fx$rank_tol
    ),
    class = "cdid_static_control_error"
  )
  expect_error(
    static_complete_fit(
      cells = fx$cells,
      or_blocks = fx$partitions$separate,
      method = "dr",
      glm_control = fx$glm_control
    ),
    class = "cdid_static_rank_tol_error"
  )

  # 5. Logit exact alias (zero column or duplicate column) and collinearity (z3 = z1 + 2 * z2)
  cell_zero_col <- static_complete_cell(
    cell_id = "c_zero", t = 2, d = 0L, unit_id = fx$unit_id,
    R = base_cell$R, G = base_cell$G, C = base_cell$C, dy = base_cell$dy,
    Z = cbind(z1 = base_cell$Z[, 1L], z_zero = rep(0.0, fx$n))
  )
  expect_error(
    static_complete_fit(
      cells = list(c_zero = cell_zero_col),
      or_blocks = c(c_zero = "b1"),
      method = "ipw",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_ipw_rank_error"
  )

  cell_collinear <- static_complete_cell(
    cell_id = "c_col", t = 2, d = 0L, unit_id = fx$unit_id,
    R = base_cell$R, G = base_cell$G, C = base_cell$C, dy = base_cell$dy,
    Z = cbind(base_cell$Z, z3 = base_cell$Z[, 1L] + 2.0 * base_cell$Z[, 2L])
  )
  expect_error(
    static_complete_fit(
      cells = list(c_col = cell_collinear),
      or_blocks = c(c_col = "b1"),
      method = "ipw",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_ipw_rank_error"
  )
  expect_error(
    static_complete_fit(
      cells = list(c_col = cell_collinear),
      or_blocks = c(c_col = "b1"),
      method = "reg",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_or_wls_rank_error"
  )

  # 6. Logit non-convergence under caller-supplied glm_control (maxit = 1L)
  expect_error(
    static_complete_fit(
      cells = list(c_t2_d0 = base_cell),
      or_blocks = c(c_t2_d0 = "b1"),
      method = "ipw",
      glm_control = list(epsilon = 1e-12, maxit = 1L, trace = FALSE),
      rank_tol = fx$rank_tol
    ),
    class = "cdid_ipw_convergence_error"
  )

  # 7. Logit coordinate separation (accepted separation lemma witness:
  # z_sep == 0 has both switchers and stayers, while z_sep > 0 has only switchers;
  # glm.fit stops before hitting 1 - 10*eps boundary, and ipw_logit_coordinate_separation
  # catches the exact coordinate separation with observed_nuisance_only = TRUE).
  z_sep <- rep(0.0, fx$n)
  sw_idx <- which(base_cell$G == 1L)
  z_sep[sw_idx[1:3]] <- c(1.0, 1.0, 1.0)
  cell_sep <- static_complete_cell(
    cell_id = "c_sep", t = 2, d = 0L, unit_id = fx$unit_id,
    R = base_cell$R, G = base_cell$G, C = base_cell$C, dy = base_cell$dy,
    Z = cbind(z_sep = z_sep)
  )
  err_sep <- tryCatch(
    static_complete_fit(
      cells = list(c_sep = cell_sep),
      or_blocks = c(c_sep = "b1"),
      method = "dr",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    error = function(e) e
  )
  expect_s3_class(err_sep, "cdid_ipw_separation_error")
  expect_true(isTRUE(err_sep$diagnostics$observed_nuisance_only))
  expect_match(err_sep$diagnostics$identification_claim, "not ATT nonidentification", fixed = TRUE)

  # 8. Logit boundary probability refusal (extreme separation hitting 0/1 boundary)
  z_bnd <- ifelse(base_cell$G == 1L, 50.0, -50.0)
  cell_bnd <- static_complete_cell(
    cell_id = "c_bnd", t = 2, d = 0L, unit_id = fx$unit_id,
    R = base_cell$R, G = base_cell$G, C = base_cell$C, dy = base_cell$dy,
    Z = cbind(z_bnd = z_bnd)
  )
  expect_error(
    static_complete_fit(
      cells = list(c_bnd = cell_bnd),
      or_blocks = c(c_bnd = "b1"),
      method = "ipw",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    ),
    class = "cdid_ipw_boundary_error"
  )
})

test_that("non_target_or_cells and materialize_if = FALSE preserve estimates, rank checks, and skip IF/cov allocation", {
  fx <- build_static_panel_fixture()

  # Subset of target cells A = {c_t2_d0, c_t2_d1} sharing OR slopes with non-target U \ A = {c_t3_d0, c_t3_d1}
  target_ids <- c("c_t2_d0", "c_t2_d1")
  nt_ids <- c("c_t3_d0", "c_t3_d1")
  part_u <- c(c_t2_d0 = "b_d0", c_t2_d1 = "b_d1", c_t3_d0 = "b_d0", c_t3_d1 = "b_d1")

  # Reference full-U fit evaluated on A
  full_fit <- static_complete_fit(
    cells = fx$cells[c(target_ids, nt_ids)],
    or_blocks = part_u,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    aggregate_cells = target_ids,
    aggregate_signs = c(c_t2_d0 = 1, c_t2_d1 = -1)
  )

  # Direct fit on A with non_target_or_cells = U \ A
  nt_fit <- static_complete_fit(
    cells = fx$cells[target_ids],
    or_blocks = part_u,
    method = "dr",
    glm_control = fx$glm_control,
    rank_tol = fx$rank_tol,
    aggregate_cells = target_ids,
    aggregate_signs = c(c_t2_d0 = 1, c_t2_d1 = -1),
    non_target_or_cells = fx$cells[nt_ids],
    materialize_if = TRUE
  )
  expect_equal(nt_fit$tau_hat, full_fit$tau_hat[target_ids], tolerance = 1e-14)
  expect_equal(nt_fit$theta_hat, full_fit$theta_hat[target_ids], tolerance = 1e-14)
  expect_equal(nt_fit$phi_tau, full_fit$phi_tau[, target_ids, drop = FALSE], tolerance = 1e-14)
  expect_equal(nt_fit$cov_tau, full_fit$cov_tau[target_ids, target_ids, drop = FALSE], tolerance = 1e-14)
  expect_equal(nt_fit$aggregate$delta_hat, full_fit$aggregate$delta_hat, tolerance = 1e-14)
  expect_equal(nt_fit$aggregate$phi, full_fit$aggregate$phi, tolerance = 1e-14)

  # Non-target stayer-only cell (0 switchers) allowed when allow_empty_switchers = TRUE
  base_c <- fx$cells[[1L]]
  R_st_only <- base_c$C
  G_st_only <- integer(fx$n)
  C_st_only <- base_c$C
  expect_error(
    static_complete_cell(
      cell_id = "c_st_only", t = 4, d = 0L, unit_id = fx$unit_id,
      R = R_st_only, G = G_st_only, C = C_st_only, dy = base_c$dy,
      X_tilde = base_c$X_tilde, Z = base_c$Z,
      allow_empty_switchers = FALSE
    ),
    class = "cdid_static_support_error"
  )
  c_st_only <- static_complete_cell(
    cell_id = "c_st_only", t = 4, d = 0L, unit_id = fx$unit_id,
    R = R_st_only, G = G_st_only, C = C_st_only, dy = base_c$dy,
    X_tilde = base_c$X_tilde, Z = base_c$Z,
    allow_empty_switchers = TRUE
  )
  expect_identical(c_st_only$n_switchers, 0L)
  expect_gt(c_st_only$n_stayers, 0L)

  # Point-only materialize_if = FALSE across ipw, reg, dr
  for (m in c("ipw", "reg", "dr")) {
    fit_if_true <- static_complete_fit(
      cells = fx$cells[target_ids],
      or_blocks = part_u,
      method = m,
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol,
      aggregate_cells = target_ids,
      aggregate_signs = c(c_t2_d0 = 1, c_t2_d1 = -1),
      non_target_or_cells = fx$cells[nt_ids],
      materialize_if = TRUE
    )
    fit_if_false <- static_complete_fit(
      cells = fx$cells[target_ids],
      or_blocks = part_u,
      method = m,
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol,
      aggregate_cells = target_ids,
      aggregate_signs = c(c_t2_d0 = 1, c_t2_d1 = -1),
      non_target_or_cells = fx$cells[nt_ids],
      materialize_if = FALSE
    )
    expect_equal(fit_if_false$tau_hat, fit_if_true$tau_hat, tolerance = 1e-14)
    expect_equal(fit_if_false$theta_hat, fit_if_true$theta_hat, tolerance = 1e-14)
    expect_equal(fit_if_false$aggregate$delta_hat, fit_if_true$aggregate$delta_hat, tolerance = 1e-14)
    expect_null(fit_if_false$phi_theta)
    expect_null(fit_if_false$phi_tau)
    expect_null(fit_if_false$cov_theta)
    expect_null(fit_if_false$cov_tau)
    expect_null(fit_if_false$vhat_theta)
    expect_null(fit_if_false$vhat_tau)
    expect_null(fit_if_false$aggregate$phi)
    expect_null(fit_if_false$aggregate$joint_if)
    expect_null(fit_if_false$aggregate$joint_cov)
  }
})
