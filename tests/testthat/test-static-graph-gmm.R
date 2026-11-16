library(testthat)

source(testthat::test_path("reference", "static_graph_gmm_dense_reference.R"), local = TRUE)

make_test_bounds <- function(
  max_matrix_entries = 50000L,
  max_batch_rows = 16L,
  max_targets = 20L,
  max_open_clusters = 1L,
  allow_materialize_if = TRUE,
  max_materialized_if_entries = 50000L
) {
  gmm_bounds_control(
    max_matrix_entries = max_matrix_entries,
    max_batch_rows = max_batch_rows,
    max_targets = max_targets,
    max_open_clusters = max_open_clusters,
    allow_materialize_if = allow_materialize_if,
    max_materialized_if_entries = max_materialized_if_entries
  )
}

test_that("1. Single-edge exact identification matches closed-form contrast for IPW/DR and OR (including zero/single-outcome units)", {
  bG <- static_role_block_spec(
    block_id = "cellA_G",
    cell_id = "cellA",
    role = "G",
    baseline_time = 2,
    target_time = 3,
    dates = c(2, 3),
    edges = data.frame(link_id = "e_G_23", s = 2, u = 3, stringsAsFactors = FALSE)
  )
  bC <- static_role_block_spec(
    block_id = "cellA_C",
    cell_id = "cellA",
    role = "C",
    baseline_time = 2,
    target_time = 3,
    dates = c(2, 3),
    edges = data.frame(link_id = "e_C_23", s = 2, u = 3, stringsAsFactors = FALSE)
  )

  N <- 30L
  unit_ids <- sprintf("unit_%02d", seq_len(N))
  set.seed(101)
  phi_raw <- matrix(rnorm(N * 2L), nrow = N, ncol = 2L)
  # Units 1..5 have zero/single observed outcomes so their complete-unit IF rows are identically 0
  phi_raw[1:5, ] <- 0.0
  phi <- sweep(phi_raw, 2L, colMeans(phi_raw), "-")
  dimnames(phi) <- list(unit_ids, c("e_G_23", "e_C_23"))
  mhat <- c(e_G_23 = 2.75, e_C_23 = 0.65)

  cb <- make_static_if_batch_callback(phi, unit_ids, batch_rows = 8L, link_ids = colnames(phi))
  bounds <- make_test_bounds(max_batch_rows = 8L)
  acc <- static_graph_gmm_accuracy_control()

  for (wt in c("identity", "optimal", "two_step")) {
    fit_dr <- solve_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = mhat,
      if_callback = cb,
      unit_ids = unit_ids,
      cluster_ids = NULL,
      weighting = wt,
      bounds = bounds,
      accuracy = acc,
      materialize_target_if = TRUE
    )
    ref_dr <- dense_reference_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = mhat,
      phi = phi,
      weighting = wt
    )
    expect_equal(unname(fit_dr$estimates), 2.75 - 0.65, tolerance = 1e-12)
    expect_equal(fit_dr$estimates, ref_dr$estimates, tolerance = 1e-12)
    expect_equal(fit_dr$core_se, ref_dr$se, tolerance = 1e-12)
    expect_equal(fit_dr$core_sandwich_cov, ref_dr$sandwich_cov, tolerance = 1e-12)
    expect_equal(fit_dr$asymptotic_cov, ref_dr$asymptotic_cov, tolerance = 1e-12)
    expect_equal(fit_dr$vcov_n2, ref_dr$vcov_n2, tolerance = 1e-12)
    expect_equal(as.double(fit_dr$moment_weights$Lambda), c(1.0, -1.0), tolerance = 1e-12)
    expect_equal(fit_dr$target_if, ref_dr$target_if, tolerance = 1e-12)
    expect_identical(rownames(fit_dr$target_if), unit_ids)
    expect_equal(as.double(fit_dr$target_if[1:5, 1L] - mean(phi_raw[, 1L] - phi_raw[, 2L]) * 0),
                 as.double(phi[1:5, 1L] - phi[1:5, 2L]), tolerance = 1e-12)
  }

  # OR mode: switcher block only, target is G(t)
  phi_or <- phi[, "e_G_23", drop = FALSE]
  mhat_or <- mhat["e_G_23"]
  cb_or <- make_static_if_batch_callback(phi_or, unit_ids, batch_rows = 8L, link_ids = "e_G_23")
  fit_or <- solve_static_graph_gmm(
    blocks = list(bG),
    method = "or",
    mhat = mhat_or,
    if_callback = cb_or,
    unit_ids = unit_ids,
    cluster_ids = NULL,
    weighting = "optimal",
    bounds = bounds,
    accuracy = acc,
    materialize_target_if = TRUE
  )
  ref_or <- dense_reference_static_graph_gmm(
    blocks = list(bG),
    method = "or",
    mhat = mhat_or,
    phi = phi_or,
    weighting = "optimal"
  )
  expect_equal(unname(fit_or$estimates), 2.75, tolerance = 1e-12)
  expect_equal(fit_or$core_se, ref_or$se, tolerance = 1e-12)
  expect_equal(as.double(fit_or$moment_weights$Lambda), 1.0, tolerance = 1e-12)
  expect_equal(fit_or$target_if, ref_or$target_if, tolerance = 1e-12)
})

test_that("2. Valid earlier-date bridge (s < t-1) identifies target without imposing early PT (`ex:bridge`)", {
  # Baseline t-1 = 2, target t = 3, earlier date 1.
  # Switcher block has ONLY bridge edges (1,2) and (1,3), missing direct edge (2,3).
  # Grounded at baseline 2 (nu(2) = 0), date 1 is an unrestricted free parameter nu(1):
  #   m_(1,2) = nu(2) - nu(1) = -nu(1)
  #   m_(1,3) = nu(3) - nu(1)
  # Therefore nu(3) = m_(1,3) - m_(1,2) and nu(1) = -m_(1,2) != 0.
  bG_bridge <- static_role_block_spec(
    block_id = "bridge_G",
    cell_id = "cell_bridge",
    role = "G",
    baseline_time = 2,
    target_time = 3,
    dates = c(1, 2, 3),
    edges = data.frame(
      link_id = c("G_12", "G_13"),
      s = c(1, 1),
      u = c(2, 3),
      stringsAsFactors = FALSE
    )
  )
  bC_direct <- static_role_block_spec(
    block_id = "bridge_C",
    cell_id = "cell_bridge",
    role = "C",
    baseline_time = 2,
    target_time = 3,
    dates = c(2, 3),
    edges = data.frame(
      link_id = "C_23",
      s = 2,
      u = 3,
      stringsAsFactors = FALSE
    )
  )

  N <- 36L
  unit_ids <- sprintf("u%02d", seq_len(N))
  set.seed(202)
  phi <- matrix(rnorm(N * 3L), nrow = N, ncol = 3L, dimnames = list(unit_ids, c("G_12", "G_13", "C_23")))
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  # Choose m_(1,2) = -4.2 (so unrestricted early level nu^G(1) = +4.2 != 0), m_(1,3) = -1.7, m^C_(2,3) = 0.8
  mhat <- c(G_12 = -4.2, G_13 = -1.7, C_23 = 0.8)

  cb <- make_static_if_batch_callback(phi, unit_ids, batch_rows = 10L, link_ids = names(mhat))
  bounds <- make_test_bounds(max_batch_rows = 10L)
  acc <- static_graph_gmm_accuracy_control()

  fit <- solve_static_graph_gmm(
    blocks = list(bG_bridge, bC_direct),
    method = "ipw",
    mhat = mhat,
    if_callback = cb,
    unit_ids = unit_ids,
    cluster_ids = NULL,
    weighting = "optimal",
    bounds = bounds,
    accuracy = acc,
    materialize_target_if = TRUE
  )
  ref <- dense_reference_static_graph_gmm(
    blocks = list(bG_bridge, bC_direct),
    method = "ipw",
    mhat = mhat,
    phi = phi,
    weighting = "optimal"
  )

  # Exact bridge algebra: G(3) - C(3) = (m_(1,3) - m_(1,2)) - m^C_(2,3) = (-1.7 - (-4.2)) - 0.8 = 1.7
  expect_equal(unname(fit$estimates), (-1.7 - (-4.2)) - 0.8, tolerance = 1e-12)
  expect_equal(fit$estimates, ref$estimates, tolerance = 1e-12)
  expect_equal(fit$core_se, ref$se, tolerance = 1e-12)
  expect_equal(as.double(fit$moment_weights$Lambda), c(-1.0, 1.0, -1.0), tolerance = 1e-12)
  # Early date t=1 parameter level is unrestricted and equals -m_(1,2) = +4.2 (not anchored to 0)
  p_tab <- fit$mapping$parameters
  nu_G_1 <- p_tab$level_estimate[p_tab$block_id == "bridge_G" & p_tab$time == 1]
  nu_G_3 <- p_tab$level_estimate[p_tab$block_id == "bridge_G" & p_tab$time == 3]
  expect_equal(nu_G_1, 4.2, tolerance = 1e-12)
  expect_equal(nu_G_3, 2.5, tolerance = 1e-12)
})

test_that("3. Unequal role graphs across switcher and stayer blocks match dense reference", {
  # Switcher block: dates {1, 2, 3, 4}, baseline = 2, target = 4, 5 edges (overid df = 2)
  bG <- static_role_block_spec(
    block_id = "uneq_G",
    cell_id = "cell_uneq",
    role = "switcher",
    baseline_time = 2,
    target_time = 4,
    dates = c(1, 2, 3, 4),
    edges = data.frame(
      link_id = c("G_12", "G_14", "G_23", "G_24", "G_34"),
      s = c(1, 1, 2, 2, 3),
      u = c(2, 4, 3, 4, 4),
      stringsAsFactors = FALSE
    )
  )
  # Stayer block: dates {2, 3, 4}, baseline = 2, target = 4, 2 edges (just-identified tree, df = 0)
  bC <- static_role_block_spec(
    block_id = "uneq_C",
    cell_id = "cell_uneq",
    role = "stayer",
    baseline_time = 2,
    target_time = 4,
    dates = c(2, 3, 4),
    edges = data.frame(
      link_id = c("C_23", "C_34"),
      s = c(2, 3),
      u = c(3, 4),
      stringsAsFactors = FALSE
    )
  )

  N <- 45L
  unit_ids <- sprintf("id_%02d", seq_len(N))
  link_ids <- c("G_12", "G_14", "G_23", "G_24", "G_34", "C_23", "C_34")
  set.seed(303)
  L_factor <- matrix(rnorm(7L * 7L), 7L, 7L)
  phi <- matrix(rnorm(N * 7L), N, 7L) %*% L_factor
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  dimnames(phi) <- list(unit_ids, link_ids)
  mhat <- c(G_12 = 0.9, G_14 = 3.1, G_23 = 1.1, G_24 = 2.3, G_34 = 1.0, C_23 = 0.4, C_34 = 0.3)

  cb <- make_static_if_batch_callback(phi, unit_ids, batch_rows = 12L, link_ids = link_ids)
  bounds <- make_test_bounds(max_batch_rows = 12L)
  acc <- static_graph_gmm_accuracy_control()

  for (wt in c("identity", "two_step")) {
    fit <- solve_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = mhat,
      if_callback = cb,
      unit_ids = unit_ids,
      cluster_ids = NULL,
      weighting = wt,
      bounds = bounds,
      accuracy = acc,
      materialize_target_if = TRUE
    )
    ref <- dense_reference_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = mhat,
      phi = phi,
      weighting = wt
    )
    expect_equal(fit$diagnostics$overid_df, 2L)
    expect_equal(fit$estimates, ref$estimates, tolerance = 1e-11)
    expect_equal(fit$core_se, ref$se, tolerance = 1e-11)
    expect_equal(fit$core_sandwich_cov, ref$sandwich_cov, tolerance = 1e-11)
    expect_equal(fit$moment_weights$Lambda, ref$Lambda, tolerance = 1e-11)
    expect_equal(fit$target_if, ref$target_if, tolerance = 1e-11)
  }
})

test_that("4. Two cells with nonzero cross-covariance: joint two-step differs from blockwise/separate GMM while identity matches", {
  # Cell 1: overidentified triangle on {1, 2, 3} for G, single edge (2,3) for C
  b1G <- static_role_block_spec(
    block_id = "c1_G", cell_id = "cell_1", role = "G",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(link_id = c("c1_G_12", "c1_G_13", "c1_G_23"), s = c(1, 1, 2), u = c(2, 3, 3))
  )
  b1C <- static_role_block_spec(
    block_id = "c1_C", cell_id = "cell_1", role = "C",
    baseline_time = 2, target_time = 3, dates = c(2, 3),
    edges = data.frame(link_id = "c1_C_23", s = 2, u = 3)
  )
  # Cell 2: overidentified triangle on {2, 3, 4} for G, single edge (3,4) for C
  b2G <- static_role_block_spec(
    block_id = "c2_G", cell_id = "cell_2", role = "G",
    baseline_time = 3, target_time = 4, dates = c(2, 3, 4),
    edges = data.frame(link_id = c("c2_G_23", "c2_G_24", "c2_G_34"), s = c(2, 2, 3), u = c(3, 4, 4))
  )
  b2C <- static_role_block_spec(
    block_id = "c2_C", cell_id = "cell_2", role = "C",
    baseline_time = 3, target_time = 4, dates = c(3, 4),
    edges = data.frame(link_id = "c2_C_34", s = 3, u = 4)
  )

  all_blocks <- list(b1G, b1C, b2G, b2C)
  link_ids <- c("c1_G_12", "c1_G_13", "c1_G_23", "c1_C_23", "c2_G_23", "c2_G_24", "c2_G_34", "c2_C_34")
  N <- 60L
  unit_ids <- sprintf("u_%02d", seq_len(N))
  set.seed(404)
  Z <- matrix(rnorm(N * 8L), N, 8L)
  # Induce strong cross-cell and cross-role correlation between overidentification residuals and target edges
  C_cross <- diag(8L)
  C_cross[3L, 5L] <- 0.85
  C_cross[1L, 7L] <- -0.75
  C_cross[4L, 6L] <- 0.80
  C_cross[2L, 8L] <- -0.70
  phi <- Z %*% C_cross
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  dimnames(phi) <- list(unit_ids, link_ids)

  # Point moments with nonzero overidentification residuals in both cells
  mhat <- c(
    c1_G_12 = 1.0, c1_G_13 = 3.4, c1_G_23 = 1.9, c1_C_23 = 0.5,
    c2_G_23 = 0.8, c2_G_24 = 2.9, c2_G_34 = 1.4, c2_C_34 = 0.3
  )

  bounds <- make_test_bounds(max_batch_rows = 15L)
  acc <- static_graph_gmm_accuracy_control()
  cb_all <- make_static_if_batch_callback(phi, unit_ids, batch_rows = 15L, link_ids = link_ids)

  fit_joint_opt <- solve_static_graph_gmm(
    blocks = all_blocks, method = "dr", mhat = mhat,
    if_callback = cb_all, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "two_step", bounds = bounds, accuracy = acc, materialize_target_if = TRUE
  )
  fit_joint_id <- solve_static_graph_gmm(
    blocks = all_blocks, method = "dr", mhat = mhat,
    if_callback = cb_all, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "identity", bounds = bounds, accuracy = acc, materialize_target_if = TRUE
  )

  ref_joint_opt <- dense_reference_static_graph_gmm(all_blocks, "dr", mhat, phi, weighting = "two_step")
  ref_joint_id <- dense_reference_static_graph_gmm(all_blocks, "dr", mhat, phi, weighting = "identity")
  ref_blockdiag_opt <- dense_reference_static_graph_gmm(
    all_blocks, "dr", mhat, phi, weighting = "two_step", mutate_block_diagonal_V = TRUE
  )

  # Separate cell-by-cell fits
  idx_c1 <- 1:4
  idx_c2 <- 5:8
  fit_c1_opt <- solve_static_graph_gmm(
    blocks = list(b1G, b1C), method = "dr", mhat = mhat[idx_c1],
    if_callback = make_static_if_batch_callback(phi[, idx_c1, drop = FALSE], unit_ids, 15L, link_ids[idx_c1]),
    unit_ids = unit_ids, cluster_ids = NULL, weighting = "two_step", bounds = bounds, accuracy = acc
  )
  fit_c2_opt <- solve_static_graph_gmm(
    blocks = list(b2G, b2C), method = "dr", mhat = mhat[idx_c2],
    if_callback = make_static_if_batch_callback(phi[, idx_c2, drop = FALSE], unit_ids, 15L, link_ids[idx_c2]),
    unit_ids = unit_ids, cluster_ids = NULL, weighting = "two_step", bounds = bounds, accuracy = acc
  )
  fit_c1_id <- solve_static_graph_gmm(
    blocks = list(b1G, b1C), method = "dr", mhat = mhat[idx_c1],
    if_callback = make_static_if_batch_callback(phi[, idx_c1, drop = FALSE], unit_ids, 15L, link_ids[idx_c1]),
    unit_ids = unit_ids, cluster_ids = NULL, weighting = "identity", bounds = bounds, accuracy = acc
  )
  fit_c2_id <- solve_static_graph_gmm(
    blocks = list(b2G, b2C), method = "dr", mhat = mhat[idx_c2],
    if_callback = make_static_if_batch_callback(phi[, idx_c2, drop = FALSE], unit_ids, 15L, link_ids[idx_c2]),
    unit_ids = unit_ids, cluster_ids = NULL, weighting = "identity", bounds = bounds, accuracy = acc
  )

  # Joint optimal matches dense joint reference
  expect_equal(fit_joint_opt$estimates, ref_joint_opt$estimates, tolerance = 1e-11)
  expect_equal(fit_joint_opt$core_sandwich_cov, ref_joint_opt$sandwich_cov, tolerance = 1e-11)
  expect_equal(fit_joint_opt$moment_weights$Lambda, ref_joint_opt$Lambda, tolerance = 1e-11)

  # Identity weighting has block-diagonal D' I D, so joint identity == separate identity
  expect_equal(fit_joint_id$estimates, ref_joint_id$estimates, tolerance = 1e-12)
  expect_equal(fit_joint_id$estimates, c(fit_c1_id$estimates, fit_c2_id$estimates), tolerance = 1e-12)
  # Yet joint identity still preserves the non-zero off-diagonal cross-cell sandwich covariance!
  expect_gt(abs(fit_joint_id$core_sandwich_cov["cell_1", "cell_2"]), 1e-3)

  # Joint two-step strictly differs from separate per-cell two-step and block-diagonal V two-step!
  est_separate_opt <- c(fit_c1_opt$estimates, fit_c2_opt$estimates)
  expect_gt(max(abs(fit_joint_opt$estimates - est_separate_opt)), 1e-2)
  expect_gt(max(abs(fit_joint_opt$estimates - ref_blockdiag_opt$estimates)), 1e-2)
  # Cross-role / cross-cell loading entries are nonzero in joint two-step
  expect_gt(max(abs(fit_joint_opt$moment_weights$Lambda[idx_c2, "cell_1"])), 1e-2)
  expect_gt(max(abs(fit_joint_opt$moment_weights$Lambda[idx_c1, "cell_2"])), 1e-2)
})

test_that("5. Singular complete-data-style IF Gram triggers existing Tikhonov (`complete_unit_if_gram` seal), while unsealed external negative-eigenvalue V is refused", {
  bG <- static_role_block_spec(
    block_id = "sing_G", cell_id = "cell_sing", role = "G",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(link_id = c("G_12", "G_13", "G_23"), s = c(1, 1, 2), u = c(2, 3, 3))
  )
  bC <- static_role_block_spec(
    block_id = "sing_C", cell_id = "cell_sing", role = "C",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(link_id = c("C_12", "C_13", "C_23"), s = c(1, 1, 2), u = c(2, 3, 3))
  )
  link_ids <- c("G_12", "G_13", "G_23", "C_12", "C_13", "C_23")

  # Complete-data-style IF redundancy: chi_(1,3) == chi_(1,2) + chi_(2,3) for every unit i,
  # plus N = 5 < m_R = 6 so the 6x6 Gram matrix has rank at most 4 and floating roundoff zero/negative eigenvalues.
  N <- 5L
  unit_ids <- paste0("u", seq_len(N))
  set.seed(505)
  g12 <- rnorm(N) * 1e3
  g23 <- rnorm(N) * 1e3
  c12 <- rnorm(N) * 1e3
  c23 <- rnorm(N) * 1e3
  phi <- cbind(
    G_12 = g12,
    G_13 = g12 + g23,
    G_23 = g23,
    C_12 = c12,
    C_13 = c12 + c23,
    C_23 = c23
  )
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  rownames(phi) <- unit_ids

  mhat <- c(G_12 = 1.2, G_13 = 3.5, G_23 = 2.1, C_12 = 0.4, C_13 = 1.1, C_23 = 0.6)
  cb <- make_static_if_batch_callback(phi, unit_ids, batch_rows = 3L, link_ids = link_ids)
  bounds <- make_test_bounds(max_batch_rows = 3L)
  acc <- static_graph_gmm_accuracy_control()

  fit <- solve_static_graph_gmm(
    blocks = list(bG, bC),
    method = "dr",
    mhat = mhat,
    if_callback = cb,
    unit_ids = unit_ids,
    cluster_ids = NULL,
    weighting = "optimal",
    bounds = bounds,
    accuracy = acc,
    materialize_target_if = TRUE
  )
  ref <- dense_reference_static_graph_gmm(
    blocks = list(bG, bC),
    method = "dr",
    mhat = mhat,
    phi = phi,
    weighting = "optimal"
  )

  expect_identical(fit$diagnostics$covariance_handling$handler_used, "tikhonov")
  expect_true(fit$diagnostics$stream_summary$gram_provenance_sealed)
  expect_equal(fit$diagnostics$covariance_handling$lambda, ref$ridge_lambda, tolerance = 1e-12)
  expect_equal(fit$estimates, ref$estimates, tolerance = 1e-10)
  expect_equal(fit$moment_weights$Lambda, ref$Lambda, tolerance = 1e-10)
  # Crucial: vcov is the actual empirical sandwich Lambda' V Lambda, NOT regularized inverse-information!
  expect_equal(fit$core_sandwich_cov, ref$sandwich_cov, tolerance = 1e-10)
  expect_equal(fit$core_se, ref$se, tolerance = 1e-10)

  # Externally supplied matrix with a genuine negative eigenvalue does NOT get fabricated provenance and is refused
  V_bad <- ref$V
  V_bad[1L, 1L] <- -1.0
  dimnames(V_bad) <- NULL
  expect_error(
    solve_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = mhat,
      if_callback = cb,
      unit_ids = unit_ids,
      cluster_ids = NULL,
      weighting = "optimal",
      bounds = bounds,
      accuracy = acc,
      V = V_bad
    ),
    class = "cdid_gmm_error"
  )
})

test_that("6. Column/row (edge, block, cell, and unit) reordering invariance with consistent metadata", {
  b1G <- static_role_block_spec(
    block_id = "c1_G", cell_id = "cell_1", role = "G",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(link_id = c("c1_G_12", "c1_G_13", "c1_G_23"), s = c(1, 1, 2), u = c(2, 3, 3))
  )
  b1C <- static_role_block_spec(
    block_id = "c1_C", cell_id = "cell_1", role = "C",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(link_id = c("c1_C_12", "c1_C_23"), s = c(1, 2), u = c(2, 3))
  )
  b2G <- static_role_block_spec(
    block_id = "c2_G", cell_id = "cell_2", role = "G",
    baseline_time = 2, target_time = 4, dates = c(2, 3, 4),
    edges = data.frame(link_id = c("c2_G_23", "c2_G_24", "c2_G_34"), s = c(2, 2, 3), u = c(3, 4, 4))
  )
  b2C <- static_role_block_spec(
    block_id = "c2_C", cell_id = "cell_2", role = "C",
    baseline_time = 2, target_time = 4, dates = c(2, 4),
    edges = data.frame(link_id = "c2_C_24", s = 2, u = 4)
  )

  blocks_orig <- list(b1G, b1C, b2G, b2C)
  link_ids_orig <- c("c1_G_12", "c1_G_13", "c1_G_23", "c1_C_12", "c1_C_23", "c2_G_23", "c2_G_24", "c2_G_34", "c2_C_24")
  N <- 40L
  unit_ids_orig <- sprintf("u_%02d", seq_len(N))
  cluster_ids_orig <- sprintf("cl_%02d", rep(1:10, length.out = N))

  set.seed(606)
  phi_orig <- matrix(rnorm(N * 9L), N, 9L, dimnames = list(unit_ids_orig, link_ids_orig))
  phi_orig <- sweep(phi_orig, 2L, colMeans(phi_orig), "-")
  mhat_orig <- stats::setNames(c(0.5, 2.2, 1.4, 0.2, 0.4, 1.1, 2.6, 1.3, 0.7), link_ids_orig)

  bounds <- make_test_bounds(max_batch_rows = 11L)
  acc <- static_graph_gmm_accuracy_control()

  fit_orig <- solve_static_graph_gmm(
    blocks = blocks_orig,
    method = "dr",
    mhat = mhat_orig,
    if_callback = make_static_if_batch_callback(phi_orig, unit_ids_orig, 11L, link_ids_orig),
    unit_ids = unit_ids_orig,
    cluster_ids = cluster_ids_orig,
    weighting = "optimal",
    bounds = bounds,
    accuracy = acc,
    cell_ids = c("cell_1", "cell_2"),
    materialize_target_if = TRUE
  )

  # Permute edges within blocks, permute block order, permute cell order, and permute unit row order
  b1G_perm <- static_role_block_spec(
    block_id = "c1_G", cell_id = "cell_1", role = "G",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = b1G$edges[c(3L, 1L, 2L), c("link_id", "s", "u")]
  )
  b1C_perm <- static_role_block_spec(
    block_id = "c1_C", cell_id = "cell_1", role = "C",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = b1C$edges[c(2L, 1L), c("link_id", "s", "u")]
  )
  b2G_perm <- static_role_block_spec(
    block_id = "c2_G", cell_id = "cell_2", role = "G",
    baseline_time = 2, target_time = 4, dates = c(2, 3, 4),
    edges = b2G$edges[c(2L, 3L, 1L), c("link_id", "s", "u")]
  )
  blocks_perm <- list(b2C, b1C_perm, b2G_perm, b1G_perm)
  link_ids_perm <- unlist(lapply(blocks_perm, function(b) b$edges$link_id), use.names = FALSE)

  set.seed(607)
  row_perm <- sample.int(N)
  unit_ids_perm <- unit_ids_orig[row_perm]
  cluster_ids_perm <- cluster_ids_orig[row_perm]
  phi_perm <- phi_orig[row_perm, link_ids_perm, drop = FALSE]
  mhat_perm <- mhat_orig[link_ids_perm]
  cells_perm <- c("cell_2", "cell_1")

  fit_perm <- solve_static_graph_gmm(
    blocks = blocks_perm,
    method = "dr",
    mhat = mhat_perm,
    if_callback = make_static_if_batch_callback(phi_perm, unit_ids_perm, 7L, link_ids_perm),
    unit_ids = unit_ids_perm,
    cluster_ids = cluster_ids_perm,
    weighting = "optimal",
    bounds = bounds,
    accuracy = acc,
    cell_ids = cells_perm,
    materialize_target_if = TRUE
  )

  expect_equal(fit_perm$estimates[c("cell_1", "cell_2")], fit_orig$estimates, tolerance = 1e-11)
  expect_equal(fit_perm$core_se[c("cell_1", "cell_2")], fit_orig$core_se, tolerance = 1e-11)
  expect_equal(fit_perm$core_sandwich_cov[c("cell_1", "cell_2"), c("cell_1", "cell_2")], fit_orig$core_sandwich_cov, tolerance = 1e-11)
  expect_equal(fit_perm$asymptotic_cov[c("cell_1", "cell_2"), c("cell_1", "cell_2")], fit_orig$asymptotic_cov, tolerance = 1e-11)
  expect_equal(fit_perm$vcov_n2[c("cell_1", "cell_2"), c("cell_1", "cell_2")], fit_orig$vcov_n2, tolerance = 1e-11)
  expect_equal(
    fit_perm$moment_weights$Lambda[link_ids_orig, c("cell_1", "cell_2")],
    fit_orig$moment_weights$Lambda,
    tolerance = 1e-11
  )
  expect_equal(
    fit_perm$target_if[unit_ids_orig, c("cell_1", "cell_2")],
    fit_orig$target_if,
    tolerance = 1e-11
  )
  expect_identical(fit_perm$mapping$unit_ids, unit_ids_perm)
  expect_identical(fit_perm$mapping$cluster_ids, cluster_ids_perm)
})

test_that("7. Cluster inference-only: point estimates, moment covariance V, Lambda, and unit IFs are unchanged", {
  bG <- static_role_block_spec(
    block_id = "cl_G", cell_id = "cell_cl", role = "G",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(link_id = c("G_12", "G_13", "G_23"), s = c(1, 1, 2), u = c(2, 3, 3))
  )
  bC <- static_role_block_spec(
    block_id = "cl_C", cell_id = "cell_cl", role = "C",
    baseline_time = 2, target_time = 3, dates = c(2, 3),
    edges = data.frame(link_id = "C_23", s = 2, u = 3)
  )

  N <- 40L
  unit_ids <- sprintf("u%02d", seq_len(N))
  # Interleaved (non-contiguous) cluster IDs with positive within-cluster correlation
  cluster_ids <- sprintf("cluster_%02d", rep(1:8, length.out = N))
  set.seed(707)
  cl_effects <- matrix(rnorm(8L * 4L, sd = 1.5), nrow = 8L, ncol = 4L)
  cl_idx <- match(cluster_ids, sprintf("cluster_%02d", 1:8))
  phi <- cl_effects[cl_idx, ] + matrix(rnorm(N * 4L, sd = 0.5), nrow = N, ncol = 4L)
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  dimnames(phi) <- list(unit_ids, c("G_12", "G_13", "G_23", "C_23"))
  mhat <- c(G_12 = 0.8, G_13 = 2.5, G_23 = 1.4, C_23 = 0.3)

  cb <- make_static_if_batch_callback(phi, unit_ids, batch_rows = 9L, link_ids = colnames(phi))
  bounds <- make_test_bounds(max_batch_rows = 9L)
  acc <- static_graph_gmm_accuracy_control()

  fit_uncl <- solve_static_graph_gmm(
    blocks = list(bG, bC), method = "dr", mhat = mhat,
    if_callback = cb, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "optimal", bounds = bounds, accuracy = acc, materialize_target_if = TRUE
  )
  fit_cl <- solve_static_graph_gmm(
    blocks = list(bG, bC), method = "dr", mhat = mhat,
    if_callback = cb, unit_ids = unit_ids, cluster_ids = cluster_ids,
    weighting = "optimal", bounds = bounds, accuracy = acc, materialize_target_if = TRUE
  )
  ref_cl <- dense_reference_static_graph_gmm(
    blocks = list(bG, bC), method = "dr", mhat = mhat,
    phi = phi, cluster_ids = cluster_ids, weighting = "optimal"
  )

  # Point estimates, moment covariance V, unit-GMM weights Lambda, and full-unit target_if are IDENTICAL
  expect_identical(fit_cl$estimates, fit_uncl$estimates)
  expect_identical(fit_cl$moment_covariance, fit_uncl$moment_covariance)
  expect_identical(fit_cl$moment_weights$Lambda, fit_uncl$moment_weights$Lambda)
  expect_identical(fit_cl$target_if, fit_uncl$target_if)
  expect_identical(fit_cl$unit_sandwich_cov, fit_uncl$core_sandwich_cov)

  # Clustered sandwich covariance and SE match dense cluster-sum reference and differ from unclustered SE
  expect_equal(fit_cl$core_sandwich_cov, ref_cl$sandwich_cov, tolerance = 1e-12)
  expect_equal(fit_cl$core_se, ref_cl$se, tolerance = 1e-12)
  expect_equal(fit_cl$asymptotic_cov, ref_cl$asymptotic_cov, tolerance = 1e-12)
  expect_equal(fit_cl$vcov_n2, ref_cl$vcov_n2, tolerance = 1e-12)
  expect_gt(abs(fit_cl$core_se - fit_uncl$core_se), 1e-2)
  expect_identical(nrow(fit_cl$cluster_target_if), 8L)

  # Contiguous cluster ordering yields identical estimates, Lambda, V, and cluster sandwich vcov/SE
  ord_contig <- order(cluster_ids)
  phi_contig <- phi[ord_contig, , drop = FALSE]
  u_contig <- unit_ids[ord_contig]
  cl_contig <- cluster_ids[ord_contig]
  fit_cl_contig <- solve_static_graph_gmm(
    blocks = list(bG, bC), method = "dr", mhat = mhat,
    if_callback = make_static_if_batch_callback(phi_contig, u_contig, 9L, colnames(phi)),
    unit_ids = u_contig, cluster_ids = cl_contig,
    weighting = "optimal", bounds = bounds, accuracy = acc, materialize_target_if = TRUE
  )
  expect_equal(fit_cl_contig$estimates, fit_cl$estimates, tolerance = 1e-12)
  expect_equal(fit_cl_contig$core_sandwich_cov, fit_cl$core_sandwich_cov, tolerance = 1e-12)
  expect_equal(fit_cl_contig$core_se, fit_cl$core_se, tolerance = 1e-12)
})

test_that("8. Bounded-batch streaming succeeds when dense N x m_R exceeds max_matrix_entries", {
  b1G <- static_role_block_spec(
    block_id = "stream_1G", cell_id = "cell_s1", role = "G",
    baseline_time = 1, target_time = 4, dates = c(1, 2, 3, 4),
    edges = data.frame(
      link_id = paste0("s1_G_", c("12", "13", "14", "23", "24", "34")),
      s = c(1, 1, 1, 2, 2, 3),
      u = c(2, 3, 4, 3, 4, 4)
    )
  )
  b1C <- static_role_block_spec(
    block_id = "stream_1C", cell_id = "cell_s1", role = "C",
    baseline_time = 1, target_time = 4, dates = c(1, 2, 3, 4),
    edges = data.frame(
      link_id = paste0("s1_C_", c("12", "13", "14", "23", "24", "34")),
      s = c(1, 1, 1, 2, 2, 3),
      u = c(2, 3, 4, 3, 4, 4)
    )
  )
  b2G <- static_role_block_spec(
    block_id = "stream_2G", cell_id = "cell_s2", role = "G",
    baseline_time = 1, target_time = 4, dates = c(1, 2, 3, 4),
    edges = data.frame(
      link_id = paste0("s2_G_", c("12", "13", "14", "23", "24", "34")),
      s = c(1, 1, 1, 2, 2, 3),
      u = c(2, 3, 4, 3, 4, 4)
    )
  )
  b2C <- static_role_block_spec(
    block_id = "stream_2C", cell_id = "cell_s2", role = "C",
    baseline_time = 1, target_time = 4, dates = c(1, 2, 3, 4),
    edges = data.frame(
      link_id = paste0("s2_C_", c("12", "13", "14", "23", "24", "34")),
      s = c(1, 1, 1, 2, 2, 3),
      u = c(2, 3, 4, 3, 4, 4)
    )
  )
  blocks_s <- list(b1G, b1C, b2G, b2C)
  link_ids <- unlist(lapply(blocks_s, function(b) b$edges$link_id), use.names = FALSE)
  N <- 7000L
  m_R <- length(link_ids) # 24; dense N * m_R = 168000
  unit_ids <- sprintf("u_%04d", seq_len(N))
  set.seed(808)
  phi_backing <- matrix(rnorm(N * m_R), nrow = N, ncol = m_R, dimnames = list(unit_ids, link_ids))
  phi_backing <- sweep(phi_backing, 2L, colMeans(phi_backing), "-")
  mhat <- stats::setNames(seq(0.5, by = 0.15, length.out = m_R), link_ids)

  # Reader function generates slices on demand; max_matrix_entries = 160000 < dense N * m_R = 168000
  cb_stream <- make_static_if_batch_callback(
    phi_source = function(idx) phi_backing[idx, , drop = FALSE],
    unit_ids = unit_ids,
    batch_rows = 50L,
    link_ids = link_ids
  )
  tight_bounds <- gmm_bounds_control(
    max_matrix_entries = 160000L,
    max_batch_rows = 50L,
    max_targets = 5L,
    max_open_clusters = 1L,
    allow_materialize_if = FALSE,
    max_materialized_if_entries = 100L
  )
  acc <- static_graph_gmm_accuracy_control()

  fit_stream <- solve_static_graph_gmm(
    blocks = blocks_s,
    method = "dr",
    mhat = mhat,
    if_callback = cb_stream,
    unit_ids = unit_ids,
    cluster_ids = NULL,
    weighting = "optimal",
    bounds = tight_bounds,
    accuracy = acc,
    materialize_target_if = FALSE
  )
  ref_stream <- dense_reference_static_graph_gmm(blocks_s, "dr", mhat, phi_backing, weighting = "optimal")

  expect_null(fit_stream$target_if)
  expect_equal(fit_stream$diagnostics$stream_summary$n_batches, 140L)
  expect_equal(fit_stream$diagnostics$stream_summary$max_batch_rows_seen, 50L)
  expect_equal(fit_stream$estimates, ref_stream$estimates, tolerance = 1e-11)
  expect_equal(fit_stream$core_se, ref_stream$se, tolerance = 1e-11)
})

test_that("9. Disconnected graphs, missing baseline/target, role mismatches, and alignment errors are refused", {
  # Disconnected graph: {1, 2} and {3, 4} with no connecting edge
  expect_error(
    static_role_block_spec(
      block_id = "disc_G", cell_id = "c1", role = "G",
      baseline_time = 1, target_time = 4, dates = c(1, 2, 3, 4),
      edges = data.frame(s = c(1, 3), u = c(2, 4))
    ),
    class = "cdid_static_graph_disconnected_error"
  )

  # Isolated required endpoint (baseline t-1 or target t has degree 0) is refused,
  # whereas a harmless isolated non-target vertex is accepted as gauge-only.
  expect_error(
    static_role_block_spec(
      block_id = "iso_base_G", cell_id = "c1", role = "G",
      baseline_time = 1, target_time = 3, dates = c(1, 2, 3),
      edges = data.frame(s = 2, u = 3)
    ),
    class = "cdid_static_graph_disconnected_error"
  )
  expect_error(
    static_role_block_spec(
      block_id = "iso_targ_G", cell_id = "c1", role = "G",
      baseline_time = 1, target_time = 3, dates = c(1, 2, 3),
      edges = data.frame(s = 1, u = 2)
    ),
    class = "cdid_static_graph_disconnected_error"
  )
  iso_harmless <- static_role_block_spec(
    block_id = "iso_G", cell_id = "c1", role = "G",
    baseline_time = 1, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(s = 1, u = 3)
  )
  expect_identical(iso_harmless$n_components, 2L)
  expect_identical(iso_harmless$isolated_dates, 2.0)

  # Baseline missing from dates
  expect_error(
    static_role_block_spec(
      block_id = "nobase_G", cell_id = "c1", role = "G",
      baseline_time = 2, target_time = 3, dates = c(1, 3),
      edges = data.frame(s = 1, u = 3)
    ),
    class = "cdid_static_graph_missing_baseline_error"
  )

  # Target missing from dates
  expect_error(
    static_role_block_spec(
      block_id = "notarg_G", cell_id = "c1", role = "G",
      baseline_time = 1, target_time = 3, dates = c(1, 2),
      edges = data.frame(s = 1, u = 2)
    ),
    class = "cdid_static_graph_missing_target_error"
  )

  # Invalid edge order s >= u
  expect_error(
    static_role_block_spec(
      block_id = "badedge_G", cell_id = "c1", role = "G",
      baseline_time = 1, target_time = 2, dates = c(1, 2),
      edges = data.frame(s = 2, u = 1)
    ),
    class = "cdid_static_graph_edge_order_error"
  )

  # Duplicate edge (s, u) within a role block
  expect_error(
    static_role_block_spec(
      block_id = "dupedge_G", cell_id = "c1", role = "G",
      baseline_time = 1, target_time = 2, dates = c(1, 2),
      edges = data.frame(s = c(1, 1), u = c(2, 2))
    ),
    class = "cdid_static_graph_duplicate_edge_error"
  )

  bG <- static_role_block_spec(
    block_id = "ok_G", cell_id = "c1", role = "G",
    baseline_time = 1, target_time = 2, dates = c(1, 2),
    edges = data.frame(link_id = "G_12", s = 1, u = 2)
  )
  bC <- static_role_block_spec(
    block_id = "ok_C", cell_id = "c1", role = "C",
    baseline_time = 1, target_time = 2, dates = c(1, 2),
    edges = data.frame(link_id = "C_12", s = 1, u = 2)
  )

  # IPW/DR missing stayer block -> role mismatch
  expect_error(
    build_static_graph_gmm_design(list(bG), method = "ipw"),
    class = "cdid_static_graph_role_mismatch_error"
  )
  # OR with stayer block -> role mismatch
  expect_error(
    build_static_graph_gmm_design(list(bG, bC), method = "or"),
    class = "cdid_static_graph_role_mismatch_error"
  )

  # Misaligned mhat names
  expect_error(
    build_static_graph_gmm_design(list(bG, bC), method = "dr", mhat = c(C_12 = 1.0, G_12 = 2.0)),
    class = "cdid_static_graph_alignment_error"
  )

  # Misaligned IF batch column names
  phi_bad <- matrix(rnorm(20L), nrow = 10L, ncol = 2L, dimnames = list(paste0("u", 1:10), c("C_12", "G_12")))
  expect_error(
    solve_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = c(G_12 = 1.0, C_12 = 0.5),
      if_callback = function(consume) consume(list(unit_id = paste0("u", 1:10), phi = phi_bad)),
      unit_ids = paste0("u", 1:10),
      cluster_ids = NULL,
      weighting = "identity",
      bounds = make_test_bounds(),
      accuracy = static_graph_gmm_accuracy_control()
    ),
    class = "cdid_static_graph_alignment_error"
  )

  # Non-unit survey weights in IF batch are refused in this slice
  phi_ok <- matrix(rnorm(20L), nrow = 10L, ncol = 2L, dimnames = list(paste0("u", 1:10), c("G_12", "C_12")))
  expect_error(
    solve_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = c(G_12 = 1.0, C_12 = 0.5),
      if_callback = function(consume) consume(list(unit_id = paste0("u", 1:10), phi = phi_ok, weight = rep(2.0, 10L))),
      unit_ids = paste0("u", 1:10),
      cluster_ids = NULL,
      weighting = "identity",
      bounds = make_test_bounds(),
      accuracy = static_graph_gmm_accuracy_control()
    ),
    class = "cdid_gmm_replay_error"
  )

  # Missing caller-explicit controls are refused
  expect_error(
    solve_static_graph_gmm(
      blocks = list(bG, bC),
      method = "dr",
      mhat = c(G_12 = 1.0, C_12 = 0.5),
      if_callback = function(consume) consume(list(unit_id = paste0("u", 1:10), phi = phi_ok)),
      unit_ids = paste0("u", 1:10),
      cluster_ids = NULL,
      weighting = "identity"
    ),
    class = "cdid_static_graph_input_error"
  )
})

test_that("10. Formula mutation sensitivity: mutating target sign, cross-covariance, baseline grounding, or sandwich vcov fails checks", {
  b1G <- static_role_block_spec(
    block_id = "m_c1_G", cell_id = "cell_1", role = "G",
    baseline_time = 2, target_time = 3, dates = c(1, 2, 3),
    edges = data.frame(link_id = c("c1_G_12", "c1_G_13", "c1_G_23"), s = c(1, 1, 2), u = c(2, 3, 3))
  )
  b1C <- static_role_block_spec(
    block_id = "m_c1_C", cell_id = "cell_1", role = "C",
    baseline_time = 2, target_time = 3, dates = c(2, 3),
    edges = data.frame(link_id = "c1_C_23", s = 2, u = 3)
  )
  blocks <- list(b1G, b1C)
  link_ids <- c("c1_G_12", "c1_G_13", "c1_G_23", "c1_C_23")
  N <- 40L
  unit_ids <- sprintf("u%02d", seq_len(N))
  set.seed(1010)
  Z <- matrix(rnorm(N * 4L), N, 4L)
  Z[, 4L] <- 0.8 * Z[, 1L] + 0.6 * Z[, 4L]
  phi <- sweep(Z, 2L, colMeans(Z), "-")
  dimnames(phi) <- list(unit_ids, link_ids)
  mhat <- c(c1_G_12 = -2.5, c1_G_13 = -0.3, c1_G_23 = 1.7, c1_C_23 = 0.9)

  fit <- solve_static_graph_gmm(
    blocks = blocks, method = "dr", mhat = mhat,
    if_callback = make_static_if_batch_callback(phi, unit_ids, 10L, link_ids),
    unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "optimal", bounds = make_test_bounds(max_batch_rows = 10L),
    accuracy = static_graph_gmm_accuracy_control()
  )
  fit_id <- solve_static_graph_gmm(
    blocks = blocks, method = "dr", mhat = mhat,
    if_callback = make_static_if_batch_callback(phi, unit_ids, 10L, link_ids),
    unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "identity", bounds = make_test_bounds(max_batch_rows = 10L),
    accuracy = static_graph_gmm_accuracy_control()
  )

  # Mutation A: +1/+1 target sign instead of G(t) - C(t)
  mut_sign <- dense_reference_static_graph_gmm(blocks, "dr", mhat, phi, weighting = "optimal", mutate_target_sign = TRUE)
  expect_gt(abs(fit$estimates - mut_sign$estimates), 1.0)

  # Mutation B: block-diagonal covariance instead of full joint covariance
  mut_bdiag <- dense_reference_static_graph_gmm(blocks, "dr", mhat, phi, weighting = "optimal", mutate_block_diagonal_V = TRUE)
  expect_gt(abs(fit$estimates - mut_bdiag$estimates), 1e-2)

  # Mutation C: grounding at earliest date (t=1) instead of baseline t-1=2
  mut_ground <- dense_reference_static_graph_gmm(blocks, "dr", mhat, phi, weighting = "optimal", mutate_ground_first_date = TRUE)
  expect_gt(abs(fit$estimates - mut_ground$estimates), 1.0)

  # Mutation D: inverse-information (D' W D)^{-1} instead of empirical sandwich Lambda' V Lambda under identity weighting
  mut_invinfo <- dense_reference_static_graph_gmm(blocks, "dr", mhat, phi, weighting = "identity", mutate_inverse_info_vcov = TRUE)
  expect_gt(abs(fit_id$core_sandwich_cov - mut_invinfo$sandwich_cov), 1e-1)
})


test_that("principal static uncertainty is estimator scale and replay stays consistent", {
  b <- static_role_block_spec(
    block_id = "g", cell_id = "c", role = "G", baseline_time = 1,
    target_time = 2, dates = c(1, 2),
    edges = data.frame(link_id = "e", s = 1, u = 2)
  )
  ids <- paste0("u", 1:4)
  phi <- matrix(c(-2, 2, -2, 2), ncol = 1,
                dimnames = list(ids, "e"))
  bounds <- make_test_bounds(max_batch_rows = 2L)
  accuracy <- static_graph_gmm_accuracy_control()
  fit_one <- function(callback, cluster = NULL) solve_static_graph_gmm(
    blocks = list(b), method = "or", mhat = c(e = 3),
    if_callback = callback, unit_ids = ids, cluster_ids = cluster,
    weighting = "identity", bounds = bounds, accuracy = accuracy,
    materialize_target_if = TRUE
  )
  callback <- make_static_if_batch_callback(phi, ids, batch_rows = 2L, link_ids = "e")
  fit <- fit_one(callback)
  # Hand calculation: sum(psi^2)/N^2 = 16/16 = 1.
  expect_equal(unname(fit$vcov), matrix(1), tolerance = 1e-12)
  expect_equal(unname(fit$se), 1, tolerance = 1e-12)
  expect_equal(fit$cell_effects$se, 1, tolerance = 1e-12)
  expect_equal(fit$cell_effects$se_unit, 1, tolerance = 1e-12)
  expect_identical(fit$vcov, fit$vcov_n2)
  expect_identical(fit$sandwich_cov, fit$vcov)
  expect_equal(unname(fit$core_sandwich_cov), matrix(4), tolerance = 1e-12)
  # Two clusters each have sum -4 or +4: 32/16 = 2.
  cluster <- c("a", "b", "a", "b")
  clustered <- fit_one(callback, cluster)
  expect_equal(unname(clustered$vcov), matrix(2), tolerance = 1e-12)
  expect_equal(unname(clustered$se), sqrt(2), tolerance = 1e-12)
  expect_identical(clustered$estimates, fit$estimates)
  changing <- function(change_at) {
    pass <- 0L
    function(consume) {
      pass <<- pass + 1L
      callback(function(batch) {
        if (pass >= change_at) batch$phi <- 2 * batch$phi
        consume(batch)
      })
    }
  }
  # Passes1-2 make covariance,3-4 unit inference,5-6 cluster inference.
  expect_error(fit_one(changing(3L)), class = "cdid_gmm_replay_error")
  expect_error(fit_one(changing(5L), cluster), class = "cdid_gmm_replay_error")
})

test_that("12. Multi-component role graphs: component gauges, disconnected trees/cycles, isolated dates, gauge/relabeling invariance, and full-catalog Tikhonov trace", {
  N <- 40L
  unit_ids <- sprintf("u%02d", seq_len(N))
  cluster_ids <- sprintf("cl%02d", ((seq_len(N) - 1L) %% 8L) + 1L)
  bounds <- make_test_bounds(max_batch_rows = 9L)
  acc_id <- static_graph_gmm_accuracy_control()
  acc_opt <- static_graph_gmm_accuracy_control()

  # (A) Target component + disconnected tree + harmless isolated non-target vertex
  #     with unequal role-specific nonconsecutive graphs across G and C.
  # G dates: {10, 17, 25} (endpoint cycle, base=17, targ=25), {40, 52, 65} (disconnected tree), {80} (isolated vertex).
  # C dates: {12, 17, 25} (endpoint tree, base=17, targ=25), {90, 105} (disconnected tree edge).
  blk_G_tree <- static_role_block_spec(
    block_id = "cA_G",
    cell_id = "cell_A",
    role = "G",
    baseline_time = 17,
    target_time = 25,
    dates = c(10, 17, 25, 40, 52, 65, 80),
    edges = data.frame(
      link_id = c("g_10_17", "g_17_25", "g_10_25", "g_40_52", "g_52_65"),
      s = c(10, 17, 10, 40, 52),
      u = c(17, 25, 25, 52, 65),
      stringsAsFactors = FALSE
    )
  )
  blk_C_tree <- static_role_block_spec(
    block_id = "cA_C",
    cell_id = "cell_A",
    role = "C",
    baseline_time = 17,
    target_time = 25,
    dates = c(12, 17, 25, 90, 105),
    edges = data.frame(
      link_id = c("c_12_17", "c_17_25", "c_90_105"),
      s = c(12, 17, 90),
      u = c(17, 25, 105),
      stringsAsFactors = FALSE
    )
  )
  blocks_A <- list(blk_G_tree, blk_C_tree)
  mhat_A <- c(
    g_10_17 = -0.35, g_17_25 = 0.82, g_10_25 = 0.51, g_40_52 = 1.40, g_52_65 = -0.65,
    c_12_17 = 0.15, c_17_25 = 0.22, c_90_105 = -0.95
  )

  # Deterministic full-rank complete-unit IF matrix with cross-component correlation
  idx_u <- seq_len(N)
  phi_A <- cbind(
    g_10_17  = 0.50 * sin(0.31 * idx_u) + 0.15 * cos(0.73 * idx_u),
    g_17_25  = 0.65 * cos(0.29 * idx_u) - 0.20 * sin(0.53 * idx_u),
    g_10_25  = 0.45 * sin(0.47 * idx_u) + 0.25 * cos(0.19 * idx_u),
    g_40_52  = 0.70 * cos(0.41 * idx_u) + 0.30 * sin(0.31 * idx_u),
    g_52_65  = 0.55 * sin(0.61 * idx_u) - 0.25 * cos(0.29 * idx_u),
    c_12_17  = 0.40 * cos(0.37 * idx_u) + 0.18 * sin(0.67 * idx_u),
    c_17_25  = 0.60 * sin(0.23 * idx_u) - 0.22 * cos(0.43 * idx_u),
    c_90_105 = 0.75 * cos(0.59 * idx_u) + 0.35 * sin(0.23 * idx_u)
  )
  phi_A <- sweep(phi_A, 2L, colMeans(phi_A), "-")
  rownames(phi_A) <- unit_ids

  design_A <- build_static_graph_gmm_design(blocks = blocks_A, method = "dr", mhat = mhat_A)
  # All 8 valid edge rows retained in exact declared order
  expect_identical(design_A$links$link_id, names(mhat_A))
  expect_identical(design_A$rank_diagnostics$n_links, 8L)
  expect_identical(design_A$rank_diagnostics$n_vertices, 12L)
  expect_identical(design_A$rank_diagnostics$n_parameters, 10L)
  expect_identical(design_A$rank_diagnostics$n_components, 5L)
  expect_identical(design_A$rank_diagnostics$nullity, 3L)
  # rank = |vertices| - number_of_components = 12 - 5 = 7
  expect_identical(design_A$rank_diagnostics$rank, 7L)
  expect_identical(design_A$rank_diagnostics$cycle_rank, 1L)
  expect_false(design_A$rank_diagnostics$has_full_column_rank)

  D_gauge_A <- materialize_static_graph_incidence_matrix(design_A, drop_component_gauges = TRUE)
  expect_identical(dim(D_gauge_A), c(8L, 7L))
  expect_identical(qr(D_gauge_A)$rank, 7L)

  cb_A <- make_static_if_batch_callback(phi_A, unit_ids = unit_ids, batch_rows = 7L, link_ids = names(mhat_A))
  fit_A_id <- solve_static_graph_gmm(
    blocks = blocks_A, method = "dr", mhat = mhat_A,
    if_callback = cb_A, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "identity", bounds = bounds, accuracy = acc_id,
    materialize_target_if = TRUE
  )
  fit_A_opt <- solve_static_graph_gmm(
    blocks = blocks_A, method = "dr", mhat = mhat_A,
    if_callback = cb_A, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "two_step", bounds = bounds, accuracy = acc_opt,
    materialize_target_if = TRUE
  )
  fit_A_cl <- solve_static_graph_gmm(
    blocks = blocks_A, method = "dr", mhat = mhat_A,
    if_callback = cb_A, unit_ids = unit_ids, cluster_ids = cluster_ids,
    weighting = "two_step", bounds = bounds, accuracy = acc_opt,
    materialize_target_if = TRUE
  )

  ref_A_id <- dense_reference_static_graph_gmm(
    blocks = blocks_A, method = "dr", mhat = mhat_A, phi = phi_A, weighting = "identity"
  )
  ref_A_opt <- dense_reference_static_graph_gmm(
    blocks = blocks_A, method = "dr", mhat = mhat_A, phi = phi_A, weighting = "two_step"
  )
  ref_A_cl <- dense_reference_static_graph_gmm(
    blocks = blocks_A, method = "dr", mhat = mhat_A, phi = phi_A,
    cluster_ids = cluster_ids, weighting = "two_step"
  )

  expect_equal(fit_A_id$estimates, ref_A_id$estimates, tolerance = 1e-12)
  expect_equal(fit_A_id$moment_weights$Lambda, ref_A_id$Lambda, tolerance = 1e-12)
  expect_equal(fit_A_id$target_if, ref_A_id$target_if, tolerance = 1e-12)
  expect_equal(fit_A_id$vcov, ref_A_id$vcov_n2, tolerance = 1e-12)
  expect_equal(fit_A_id$se, ref_A_id$se_n2, tolerance = 1e-12)

  expect_equal(fit_A_opt$estimates, ref_A_opt$estimates, tolerance = 1e-12)
  expect_equal(fit_A_opt$moment_weights$Lambda, ref_A_opt$Lambda, tolerance = 1e-12)
  expect_equal(fit_A_opt$target_if, ref_A_opt$target_if, tolerance = 1e-12)
  expect_equal(fit_A_opt$vcov, ref_A_opt$vcov_n2, tolerance = 1e-12)
  expect_equal(fit_A_opt$se, ref_A_opt$se_n2, tolerance = 1e-12)

  # Disconnected tree edges have zero cycle rank, so their direct GMM loadings are 0 under any weight matrix
  tree_edges_A <- c("g_40_52", "g_52_65", "c_90_105")
  expect_equal(unname(fit_A_id$moment_weights$Lambda[tree_edges_A, 1L]), rep(0.0, 3L), tolerance = 1e-13)
  expect_equal(unname(fit_A_opt$moment_weights$Lambda[tree_edges_A, 1L]), rep(0.0, 3L), tolerance = 1e-13)

  # Clustering is inference-only: point estimates and Lambda are identical to unclustered two-step
  expect_equal(fit_A_cl$estimates, fit_A_opt$estimates, tolerance = 1e-14)
  expect_equal(fit_A_cl$moment_weights$Lambda, fit_A_opt$moment_weights$Lambda, tolerance = 1e-14)
  expect_equal(fit_A_cl$vcov, ref_A_cl$vcov_n2, tolerance = 1e-12)
  expect_equal(fit_A_cl$se, ref_A_cl$se_n2, tolerance = 1e-12)

  # Parameter table statuses distinguish endpoint-connected vs disconnected vs isolated non-target dates
  p_tab <- fit_A_opt$gmm_fit$parameters
  iso_row <- p_tab[p_tab$time == 80, ]
  expect_true(iso_row$is_isolated)
  expect_false(iso_row$level_estimable)
  expect_true(iso_row$is_gauge_root)
  expect_identical(iso_row$status, "isolated_unsupported")
  expect_true(is.na(iso_row$level_estimate))

  disc_rows <- p_tab[p_tab$time %in% c(40, 52, 65, 90, 105), ]
  expect_true(all(!disc_rows$is_isolated))
  expect_true(all(!disc_rows$level_estimable))
  expect_true(all(disc_rows$status == "disconnected_component"))
  expect_true(all(is.na(disc_rows$level_estimate)))

  # (B) Target component + disconnected cycle with cross-component covariance:
  #     two-step GMM loads the disconnected cycle's zero-mean circulation.
  blk_G_cyc <- static_role_block_spec(
    block_id = "cB_G",
    cell_id = "cell_B",
    role = "G",
    baseline_time = 10,
    target_time = 20,
    dates = c(10, 20, 40, 50, 60),
    edges = data.frame(
      link_id = c("g_10_20", "g_40_50", "g_50_60", "g_40_60"),
      s = c(10, 40, 50, 40),
      u = c(20, 50, 60, 60),
      stringsAsFactors = FALSE
    )
  )
  blocks_B <- list(blk_G_cyc)
  # Nonzero cycle discrepancy on (40,50)+(50,60)-(40,60) = 0.30 + 0.40 - 0.55 = 0.15 != 0
  mhat_B <- c(g_10_20 = 1.25, g_40_50 = 0.30, g_50_60 = 0.40, g_40_60 = 0.55)
  phi_B <- cbind(
    g_10_20 = 0.80 * sin(0.25 * idx_u) + 0.35 * cos(0.45 * idx_u),
    g_40_50 = 0.60 * sin(0.25 * idx_u) + 0.50 * sin(0.71 * idx_u),
    g_50_60 = 0.45 * cos(0.45 * idx_u) - 0.40 * cos(0.83 * idx_u),
    g_40_60 = -0.35 * sin(0.25 * idx_u) + 0.55 * sin(0.39 * idx_u)
  )
  phi_B <- sweep(phi_B, 2L, colMeans(phi_B), "-")
  rownames(phi_B) <- unit_ids
  cb_B <- make_static_if_batch_callback(phi_B, unit_ids = unit_ids, batch_rows = 8L, link_ids = names(mhat_B))

  fit_B_id <- solve_static_graph_gmm(
    blocks = blocks_B, method = "or", mhat = mhat_B,
    if_callback = cb_B, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "identity", bounds = bounds, accuracy = acc_id,
    materialize_target_if = TRUE
  )
  fit_B_opt <- solve_static_graph_gmm(
    blocks = blocks_B, method = "or", mhat = mhat_B,
    if_callback = cb_B, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "two_step", bounds = bounds, accuracy = acc_opt,
    materialize_target_if = TRUE
  )
  ref_B_id <- dense_reference_static_graph_gmm(
    blocks = blocks_B, method = "or", mhat = mhat_B, phi = phi_B, weighting = "identity"
  )
  ref_B_opt <- dense_reference_static_graph_gmm(
    blocks = blocks_B, method = "or", mhat = mhat_B, phi = phi_B, weighting = "two_step"
  )

  expect_equal(fit_B_id$estimates, ref_B_id$estimates, tolerance = 1e-12)
  expect_equal(fit_B_id$moment_weights$Lambda, ref_B_id$Lambda, tolerance = 1e-12)
  # Under identity weighting, disconnected cycle loadings are 0
  expect_equal(unname(fit_B_id$moment_weights$Lambda[2:4, 1L]), rep(0.0, 3L), tolerance = 1e-13)
  expect_equal(unname(fit_B_id$estimates), 1.25, tolerance = 1e-13)

  # Under two-step weighting with cross-component covariance, disconnected cycle has nonzero circulation loading
  expect_equal(fit_B_opt$estimates, ref_B_opt$estimates, tolerance = 1e-12)
  expect_equal(fit_B_opt$moment_weights$Lambda, ref_B_opt$Lambda, tolerance = 1e-12)
  expect_equal(fit_B_opt$target_if, ref_B_opt$target_if, tolerance = 1e-12)
  expect_equal(fit_B_opt$vcov, ref_B_opt$vcov_n2, tolerance = 1e-12)
  expect_equal(fit_B_opt$se, ref_B_opt$se_n2, tolerance = 1e-12)
  expect_true(all(abs(fit_B_opt$moment_weights$Lambda[2:4, 1L]) > 1e-3))
  expect_true(abs(fit_B_opt$estimates[["cell_B"]] - 1.25) > 1e-3)
  # Cycle loading is a valid zero-boundary circulation on component 2:
  # lambda_(40,50) == lambda_(50,60) == -lambda_(40,60)
  lam_cyc <- unname(fit_B_opt$moment_weights$Lambda[2:4, 1L])
  expect_equal(lam_cyc[1L], lam_cyc[2L], tolerance = 1e-12)
  expect_equal(lam_cyc[1L], -lam_cyc[3L], tolerance = 1e-12)

  # (C) Invariance to coordinate gauge choice in non-endpoint components and order-preserving date relabeling
  fit_B_high <- solve_static_graph_gmm(
    blocks = blocks_B, method = "or", mhat = mhat_B,
    if_callback = cb_B, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "two_step", bounds = bounds, accuracy = acc_opt,
    materialize_target_if = TRUE, gauge_rule = "highest"
  )
  # Custom gauge root: date 50 is parameter index 3 in blk_G_cyc (non-baseline dates are 20, 40, 50, 60)
  fit_B_cust <- solve_static_graph_gmm(
    blocks = blocks_B, method = "or", mhat = mhat_B,
    if_callback = cb_B, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "two_step", bounds = bounds, accuracy = acc_opt,
    materialize_target_if = TRUE, gauge_rule = "custom", custom_gauge_roots = 3L
  )
  expect_equal(fit_B_high$estimates, fit_B_opt$estimates, tolerance = 1e-12)
  expect_equal(fit_B_high$moment_weights$Lambda, fit_B_opt$moment_weights$Lambda, tolerance = 1e-12)
  expect_equal(fit_B_high$target_if, fit_B_opt$target_if, tolerance = 1e-12)
  expect_equal(fit_B_high$vcov, fit_B_opt$vcov, tolerance = 1e-12)
  expect_equal(fit_B_cust$estimates, fit_B_opt$estimates, tolerance = 1e-12)
  expect_equal(fit_B_cust$moment_weights$Lambda, fit_B_opt$moment_weights$Lambda, tolerance = 1e-12)
  expect_equal(fit_B_cust$target_if, fit_B_opt$target_if, tolerance = 1e-12)
  expect_equal(fit_B_cust$vcov, fit_B_opt$vcov, tolerance = 1e-12)

  # Order-preserving date relabeling s -> 3 * s + 7
  relabel_fn <- function(d) 3.0 * d + 7.0
  blk_G_relab <- static_role_block_spec(
    block_id = "cB_G",
    cell_id = "cell_B",
    role = "G",
    baseline_time = relabel_fn(10),
    target_time = relabel_fn(20),
    dates = relabel_fn(c(10, 20, 40, 50, 60)),
    edges = data.frame(
      link_id = c("g_10_20", "g_40_50", "g_50_60", "g_40_60"),
      s = relabel_fn(c(10, 40, 50, 40)),
      u = relabel_fn(c(20, 50, 60, 60)),
      stringsAsFactors = FALSE
    )
  )
  fit_B_relab <- solve_static_graph_gmm(
    blocks = list(blk_G_relab), method = "or", mhat = mhat_B,
    if_callback = cb_B, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "two_step", bounds = bounds, accuracy = acc_opt,
    materialize_target_if = TRUE
  )
  expect_equal(fit_B_relab$estimates, fit_B_opt$estimates, tolerance = 1e-13)
  expect_equal(fit_B_relab$moment_weights$Lambda, fit_B_opt$moment_weights$Lambda, tolerance = 1e-13)
  expect_equal(fit_B_relab$target_if, fit_B_opt$target_if, tolerance = 1e-13)
  expect_equal(fit_B_relab$vcov, fit_B_opt$vcov, tolerance = 1e-13)

  # (D) Singular / active Tikhonov covariance case where pruning the disconnected
  #     component would alter the full-catalog trace/scale tr(V) / M_full.
  blk_G_tik <- static_role_block_spec(
    block_id = "cT_G",
    cell_id = "cell_T",
    role = "G",
    baseline_time = 10,
    target_time = 20,
    dates = c(5, 10, 20, 40, 50),
    edges = data.frame(
      link_id = c("e_5_10", "e_10_20", "e_5_20", "e_40_50"),
      s = c(5, 10, 5, 40),
      u = c(10, 20, 20, 50),
      stringsAsFactors = FALSE
    )
  )
  mhat_T <- c(e_5_10 = 0.40, e_10_20 = 0.90, e_5_20 = 1.55, e_40_50 = 2.10)
  z1 <- sin(0.33 * idx_u) - mean(sin(0.33 * idx_u))
  z2 <- cos(0.51 * idx_u) - mean(cos(0.51 * idx_u))
  z3 <- sin(0.79 * idx_u) - mean(sin(0.79 * idx_u))
  # Collinear moment 3 (non-parallel to cycle vector) triggers Tikhonov ridge,
  # while disconnected edge e_40_50 has large variance 100 * z3^2 that enters tr(V_full) / M_full.
  phi_T <- cbind(
    e_5_10  = 1.0 * z1 + 0.2 * z2,
    e_10_20 = 0.3 * z1 + 1.4 * z2,
    e_5_20  = 2.0 * z1 - 0.5 * z2,
    e_40_50 = 10.0 * z3
  )
  rownames(phi_T) <- unit_ids
  cb_T <- make_static_if_batch_callback(phi_T, unit_ids = unit_ids, batch_rows = 8L, link_ids = names(mhat_T))

  fit_T <- solve_static_graph_gmm(
    blocks = list(blk_G_tik), method = "or", mhat = mhat_T,
    if_callback = cb_T, unit_ids = unit_ids, cluster_ids = NULL,
    weighting = "two_step", bounds = bounds, accuracy = acc_opt,
    materialize_target_if = TRUE
  )
  ref_T_full <- dense_reference_static_graph_gmm(
    blocks = list(blk_G_tik), method = "or", mhat = mhat_T, phi = phi_T, weighting = "two_step"
  )

  expect_identical(fit_T$diagnostics$covariance_handling$handler_used, "tikhonov")
  expect_equal(fit_T$diagnostics$covariance_handling$lambda, ref_T_full$ridge_lambda, tolerance = 1e-12)
  expect_equal(fit_T$estimates, ref_T_full$estimates, tolerance = 1e-11)
  expect_equal(fit_T$moment_weights$Lambda, ref_T_full$Lambda, tolerance = 1e-11)
  expect_equal(fit_T$target_if, ref_T_full$target_if, tolerance = 1e-11)
  expect_equal(fit_T$vcov, ref_T_full$vcov_n2, tolerance = 1e-11)

  # Compare against a hypothetical pruned catalog omitting e_40_50:
  # pruning alters tr(V)/M and therefore alters ridge_lambda and the target estimate.
  blk_G_pruned <- static_role_block_spec(
    block_id = "cT_G",
    cell_id = "cell_T",
    role = "G",
    baseline_time = 10,
    target_time = 20,
    dates = c(5, 10, 20),
    edges = data.frame(
      link_id = c("e_5_10", "e_10_20", "e_5_20"),
      s = c(5, 10, 5),
      u = c(10, 20, 20),
      stringsAsFactors = FALSE
    )
  )
  ref_T_pruned <- dense_reference_static_graph_gmm(
    blocks = list(blk_G_pruned), method = "or",
    mhat = mhat_T[1:3], phi = phi_T[, 1:3, drop = FALSE], weighting = "two_step"
  )
  expect_true(abs(fit_T$diagnostics$covariance_handling$lambda - ref_T_pruned$ridge_lambda) > 1e-5)
  expect_true(abs(fit_T$estimates[["cell_T"]] - ref_T_pruned$estimates[["cell_T"]]) > 1e-5)
})
