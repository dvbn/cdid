# Focused packaged tests for the private CP02 exact-GMM computation core
# (R/chained_gmm_core.R).

find_gmm_ref_file <- function(filename) {
  candidates <- c(
    file.path("reference", filename),
    file.path("tests", "testthat", "reference", filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Cannot locate packaged reference file: ", filename)
  }
  hit[1L]
}

hex16_to_double <- function(hex_vec) {
  out <- numeric(length(hex_vec))
  for (i in seq_along(hex_vec)) {
    h <- trimws(hex_vec[i])
    if (grepl("^-?0x", h, ignore.case = TRUE)) {
      out[i] <- as.numeric(h)
    } else {
      bytes <- vapply(seq(1L, 15L, by = 2L), function(pos) {
        as.raw(strtoi(substr(h, pos, pos + 1L), base = 16L))
      }, raw(1L))
      if (.Platform$endian == "little") {
        bytes <- rev(bytes)
      }
      out[i] <- readBin(bytes, what = "numeric", n = 1L, size = 8L, endian = .Platform$endian)
    }
  }
  out
}

read_stored_3x3_matrix <- function(matrix_name) {
  lines <- readLines(find_gmm_ref_file("stored_6c_and_singular_hex.csv"), warn = FALSE)
  prefix <- paste0(matrix_name, ",")
  hit <- lines[startsWith(lines, prefix)]
  if (length(hit) != 3L) {
    stop("Expected 3 rows for ", matrix_name)
  }
  M <- matrix(0, nrow = 3L, ncol = 3L)
  for (i in 1:3) {
    hex_entries <- strsplit(sub(prefix, "", hit[i], fixed = TRUE), ";", fixed = TRUE)[[1L]]
    M[i, ] <- hex16_to_double(hex_entries)
  }
  M
}

default_test_bounds <- function(
  max_matrix_entries = 50000L,
  max_batch_rows = 64L,
  max_targets = 32L,
  max_open_clusters = 1L,
  allow_materialize_if = FALSE,
  max_materialized_if_entries = 0L
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

# D12 control plumbing: allow callers to pass optional max_loading_error_abs (default NULL = not supplied)
# without selecting any hidden default or modifying existing control values.
default_test_accuracy <- function(
  weighting = "optimal",
  target_error_unit = if (weighting == "optimal") "both" else "abs",
  require_matching_scale = (weighting == "optimal"),
  max_loading_error_abs = NULL
) {
  gmm_accuracy_control(
    method = if (weighting == "optimal") "certified_eigen_chol_kkt" else "certified_qr_kkt",
    target_error_unit = target_error_unit,
    max_target_error_abs = 1e-10,
    max_target_error_se = if (weighting == "optimal") 1e-10 else NA_real_,
    max_orthogonality_defect = 1e-11,
    min_certified_lambda_rel = if (weighting == "optimal") 1e-8 else NA_real_,
    max_factor_residual_rel = 1e-12,
    min_rcond_R = if (weighting == "optimal") 1e-6 else NA_real_,
    min_rcond_U = 1e-6,
    min_sigma_rel_A = 1e-6,
    max_kkt_residual_rel = 1e-12,
    max_loading_residual_rel = 1e-12,
    require_matching_scale = require_matching_scale,
    max_loading_error_abs = max_loading_error_abs
  )
}

make_unit_if_replay_from_V <- function(
  V,
  n_units = 20L,
  batch_size = 5L,
  scale_id = "n_minus_1_centered",
  nuisance_fit_id = "nuis_fixture_v1"
) {
  m <- nrow(V)
  R <- chol(V)
  # Construct an exact centered n_units x m matrix Phi whose (n_units - 1)^{-1} Phi'Phi == V
  U_raw <- matrix(0, nrow = n_units, ncol = m)
  for (j in seq_len(m)) {
    U_raw[j, j] <- 1
    U_raw[n_units, j] <- -1
  }
  U_centered <- scale(U_raw, center = TRUE, scale = FALSE)
  Q_u <- qr.Q(qr(U_centered))[, seq_len(m), drop = FALSE]
  Phi <- sqrt(n_units - 1) * (Q_u %*% R)
  unit_ids <- sprintf("u_%03d", seq_len(n_units))

  cb <- function(consume_batch) {
    idx_starts <- seq(1L, n_units, by = batch_size)
    for (st in idx_starts) {
      en <- min(st + batch_size - 1L, n_units)
      consume_batch(list(
        unit_id = unit_ids[st:en],
        phi = Phi[st:en, , drop = FALSE],
        weight = rep(1, en - st + 1L),
        nuisance_fit_id = nuisance_fit_id
      ))
    }
  }

  gmm_if_replay_spec(
    callback = cb,
    input_row_level = "unit",
    aggregate_to_clusters = FALSE,
    n_input_rows = n_units,
    n_output_rows = n_units,
    n_positive_contributing_rows = n_units,
    centered = TRUE,
    weighted = TRUE,
    weight_type = "unit_uniform",
    normalization_factor = 1 / (n_units - 1),
    scale_id = scale_id,
    nuisance_fit_id = nuisance_fit_id,
    expected_unit_ids = unit_ids
  )
}

make_default_provenance <- function(
  n_units = 20L,
  scale_id = "n_minus_1_centered",
  nuisance_fit_id = "nuis_fixture_v1",
  reduction = list(mode = "none")
) {
  gmm_covariance_provenance(
    row_level = "unit",
    n_units = n_units,
    n_clusters = NA_real_,
    n_positive_contributing_rows = n_units,
    centered = TRUE,
    weighted = TRUE,
    weight_type = "unit_uniform",
    normalization_factor = 1 / (n_units - 1),
    scale_id = scale_id,
    nuisance_fit_id = nuisance_fit_id,
    reduction = reduction
  )
}

build_fixture_inputs <- function(fixture_name) {
  if (fixture_name == "hand_n3_att4_lone45_a") {
    # Cohort g=2, ref=1 has N3 post-nuisance at t=3 (unqualified), ATT(2,2) qualified,
    # and ATT(2,4) isolated via target_post_times; cohort g=3, ref=1 has lone 4->5 link.
    links <- data.frame(
      g = c(2, 2, 2, 3),
      s = c(1, 1, 2, 4),
      t = c(2, 3, 3, 5),
      ref = c(1, 1, 1, 1),
      functional_key = c("k_nuis", "k_nuis", "k_nuis", "k_att"),
      s_qualified = c(TRUE, TRUE, TRUE, TRUE),
      t_qualified = c(TRUE, FALSE, FALSE, TRUE),
      stringsAsFactors = FALSE
    )
    design <- construct_chained_design(
      links,
      target_post_times = data.frame(g = 2, time = 4)
    )
    mhat <- c(2, 5, 4, 7)
    V <- matrix(c(
      3, 1, 1, 1,
      1, 4, 1, 0,
      1, 1, 3, 1,
      1, 0, 1, 2
    ), nrow = 4L, byrow = TRUE)
    targets <- list(
      gmm_target_spec("N3_nuisance", 1L, 1.0),
      gmm_target_spec("ATT_2_2", 2L, 1.0),
      gmm_target_spec("ATT22_minus_N3", c(1L, 2L), c(-1.0, 1.0)),
      gmm_target_spec("lone_45_contrast", c(4L, 5L), c(-1.0, 1.0)),
      gmm_target_spec("isolated_ATT4", 3L, 1.0),
      gmm_target_spec("unanchored_ATT_3_4", 4L, 1.0),
      gmm_target_spec("unanchored_ATT_3_5", 5L, 1.0),
      gmm_target_spec("cross_unanchored_2_to_5", c(2L, 5L), c(-1.0, 1.0))
    )
    return(list(design = design, mhat = mhat, V = V, targets = targets))
  }

  if (fixture_name == "hand_n3_att4_lone45_single_cohort") {
    # Single cohort g=2, ref=1 with N3 post-nuisance at t=3, lone 4->5 post-nuisance contrast,
    # estimable ATT(2,2), and isolated ATT(2,4)
    links <- data.frame(
      g = c(2, 2, 2, 2),
      s = c(1, 1, 2, 4),
      t = c(2, 3, 3, 5),
      ref = c(1, 1, 1, 1),
      functional_key = c("k_n3", "k_n3", "k_n3", "k_45"),
      s_qualified = c(TRUE, TRUE, TRUE, FALSE),
      t_qualified = c(TRUE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    )
    design <- construct_chained_design(
      links,
      target_post_times = data.frame(g = 2, time = 4)
    )
    mhat <- c(2, 5, 4, 7)
    V <- matrix(c(
      3, 1, 1, 1,
      1, 4, 1, 0,
      1, 1, 3, 1,
      1, 0, 1, 2
    ), nrow = 4L, byrow = TRUE)
    targets <- list(
      gmm_target_spec("N3_nuisance", 1L, 1.0),
      gmm_target_spec("ATT_2_2", 4L, 1.0),
      gmm_target_spec("ATT22_minus_N3", c(1L, 4L), c(-1.0, 1.0)),
      gmm_target_spec("lone_45_contrast", c(2L, 3L), c(-1.0, 1.0)),
      gmm_target_spec("unanchored_4", 2L, 1.0),
      gmm_target_spec("unanchored_5", 3L, 1.0),
      gmm_target_spec("isolated_ATT4", 5L, 1.0)
    )
    return(list(design = design, mhat = mhat, V = V, targets = targets))
  }

  if (fixture_name == "transfer_w21_nonzero") {
    links <- data.frame(
      g = c(2, 2, 2, 2),
      s = c(0, 0, 1, 2),
      t = c(1, 2, 2, 3),
      ref = c(1, 1, 1, 1),
      functional_key = rep("k1", 4L),
      stringsAsFactors = FALSE
    )
    design <- construct_chained_design(links)
    mhat <- c(-1, 3, 2, 4)
    V <- matrix(c(
      3,  1, 0,  1,
      1,  4, 1, -1,
      0,  1, 2,  1,
      1, -1, 1,  3
    ), nrow = 4L, byrow = TRUE)
    targets <- list(
      gmm_target_spec("p1", 1L, 1.0),
      gmm_target_spec("p2", 2L, 1.0),
      gmm_target_spec("p3", 3L, 1.0),
      gmm_target_spec("p3_minus_p1", c(1L, 3L), c(-1.0, 1.0))
    )
    return(list(design = design, mhat = mhat, V = V, targets = targets))
  }

  if (fixture_name %in% c("spd_multi_cohort_opt", "spd_multi_cohort_id")) {
    # Links in exact canonical order:
    # l1: g=2, 0->1
    # l2: g=2, 1->2
    # l3: g=2, 0->2
    # l4: g=2, 2->3
    # l5: g=2, 1->3
    # l6: g=4, 4->5
    links <- data.frame(
      g = c(2, 2, 2, 2, 2, 4),
      s = c(0, 1, 0, 2, 1, 4),
      t = c(1, 2, 2, 3, 3, 5),
      ref = c(1, 1, 1, 1, 1, 1),
      functional_key = rep("k1", 6L),
      stringsAsFactors = FALSE
    )
    design <- construct_chained_design(
      links,
      target_post_times = data.frame(g = 2, time = 4)
    )
    mhat <- c(-2, 3, 6, 2, 5, 4)
    V <- matrix(c(
      4, 1, 1, 0, 0, 1,
      1, 3, 1, 1, 0, 0,
      1, 1, 5, 0, 1, 1,
      0, 1, 0, 3, 1, 0,
      0, 0, 1, 1, 4, 1,
      1, 0, 1, 0, 1, 3
    ), nrow = 6L, byrow = TRUE)
    targets <- list(
      gmm_target_spec("p1_pre", 1L, 1.0),
      gmm_target_spec("p2_att22", 2L, 1.0),
      gmm_target_spec("p3_att23", 3L, 1.0),
      gmm_target_spec("p6_minus_p5", c(5L, 6L), c(-1.0, 1.0)),
      gmm_target_spec("agg_target", c(2L, 3L, 5L, 6L), c(0.5, 0.5, 1.0, -1.0))
    )
    return(list(design = design, mhat = mhat, V = V, targets = targets))
  }

  stop("Unknown fixture: ", fixture_name)
}

test_that("1. Hand N3 nuisance, ATT(4) isolated, and lone 4->5 contrast preserve estimable targets and set unsupported levels to NA", {
  ref_df <- utils::read.csv(find_gmm_ref_file("fresh_exact_kkt_reference.csv"), stringsAsFactors = FALSE)

  for (fx_name in c("hand_n3_att4_lone45_a", "hand_n3_att4_lone45_single_cohort")) {
    fx <- build_fixture_inputs(fx_name)
    replay <- make_unit_if_replay_from_V(fx$V)
    prov <- make_default_provenance()
    res <- solve_chained_gmm_core(
      design = fx$design,
      mhat = fx$mhat,
      V = fx$V,
      provenance = prov,
      targets = fx$targets,
      weighting = "optimal",
      if_replay = replay,
      accuracy = default_test_accuracy("optimal"),
      bounds = default_test_bounds(),
      gauge_rule = "lowest",
      request_param_cov = TRUE
    )

    # No rows of W_G may be dropped
    expect_equal(nrow(res$internal_gauge$Lambda), fx$design$incidence$n_links)
    expect_equal(res$internal_gauge$gauge$n_links, 4L)
    expect_equal(res$internal_gauge$gauge$rank_k, 3L)

    # First 4 targets are estimable; remaining level/cross-component targets are non-estimable (NA)
    expect_identical(res$targets$estimable[1:4], rep(TRUE, 4L))
    expect_true(all(!res$targets$estimable[5:nrow(res$targets)]))
    expect_true(all(is.na(res$targets$estimate[5:nrow(res$targets)])))
    expect_true(all(is.na(res$targets$se_sandwich[5:nrow(res$targets)])))
    expect_true(all(is.na(res$targets$se_model[5:nrow(res$targets)])))

    # Parameter table: only connected_to_ref (level_estimable) parameters have level estimates;
    # gauge roots and unanchored/isolated vertices must be NA!
    est_p <- which(res$parameters$level_estimable)
    unest_p <- which(!res$parameters$level_estimable)
    expect_length(est_p, 2L)
    expect_true(all(!is.na(res$parameters$level_estimate[est_p])))
    expect_true(all(is.na(res$parameters$level_estimate[unest_p])))
    expect_true(any(res$parameters$is_gauge_root))
    expect_true(all(is.na(res$parameters$level_estimate[res$parameters$is_gauge_root])))

    # Compare estimable targets against exact Fraction reference
    sub_tau <- ref_df[ref_df$fixture == fx_name & ref_df$key == "tau", , drop = FALSE]
    expect_equal(res$targets$estimate[1:4], as.double(sub_tau$double), tolerance = 1e-14)
  }

  # Pure isolated ATT(4) without any 4->5 link: N3(3) and ATT(2) identified; ATT(4) isolated -> NA
  links_iso <- data.frame(
    g = c(2, 2, 2),
    s = c(1, 1, 2),
    t = c(2, 3, 3),
    ref = c(1, 1, 1),
    functional_key = c("k_n3", "k_n3", "k_n3"),
    s_qualified = c(TRUE, TRUE, TRUE),
    t_qualified = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  des_iso <- construct_chained_design(links_iso, target_post_times = data.frame(g = 2, time = 4))
  V3 <- matrix(c(3.0, 0.5, -0.2, 0.5, 2.5, 0.3, -0.2, 0.3, 2.0), nrow = 3L, byrow = TRUE)
  m3 <- c(-0.9, 1.2, 1.8)
  t_iso <- list(
    gmm_target_spec("N3_nuisance", 1L, 1.0),
    gmm_target_spec("ATT_2_2", 2L, 1.0),
    gmm_target_spec("ATT_2_4_isolated", 3L, 1.0)
  )
  res_iso <- chained_gmm_core(
    design = des_iso,
    mhat = m3,
    V = V3,
    provenance = make_default_provenance(),
    targets = t_iso,
    weighting = "optimal",
    if_replay = make_unit_if_replay_from_V(V3),
    accuracy = default_test_accuracy("optimal"),
    bounds = default_test_bounds(),
    gauge_rule = "lowest"
  )
  expect_identical(res_iso$targets$estimable, c(TRUE, TRUE, FALSE))
  expect_identical(res_iso$targets$status, c("estimable", "estimable", "isolated_unsupported"))
  expect_true(is.na(res_iso$targets$estimate[3L]))
  expect_true(is.na(res_iso$parameters$level_estimate[3L]))
})

test_that("2. Correlated block residual transfer works for both W_21 = 0 and W_21 != 0 in joint solve", {
  ref_df <- utils::read.csv(find_gmm_ref_file("fresh_exact_kkt_reference.csv"), stringsAsFactors = FALSE)

  # Case A: W_21 = 0 (`hand_n3_att4_lone45_a`: block 1 overidentified on links 1..3, block 2 on link 4 has W_21 = 0)
  fx_a <- build_fixture_inputs("hand_n3_att4_lone45_a")
  V_joint <- fx_a$V
  V_blockdiag <- V_joint
  V_blockdiag[1:3, 4] <- 0
  V_blockdiag[4, 1:3] <- 0

  res_joint_a <- solve_chained_gmm_core(
    design = fx_a$design,
    mhat = fx_a$mhat,
    V = V_joint,
    provenance = make_default_provenance(),
    targets = fx_a$targets[1:4],
    weighting = "optimal",
    if_replay = make_unit_if_replay_from_V(V_joint),
    accuracy = default_test_accuracy("optimal"),
    bounds = default_test_bounds(),
    gauge_rule = "lowest"
  )
  res_bdiag_a <- solve_chained_gmm_core(
    design = fx_a$design,
    mhat = fx_a$mhat,
    V = V_blockdiag,
    provenance = make_default_provenance(),
    targets = fx_a$targets[1:4],
    weighting = "optimal",
    if_replay = make_unit_if_replay_from_V(V_blockdiag),
    accuracy = default_test_accuracy("optimal"),
    bounds = default_test_bounds(),
    gauge_rule = "lowest"
  )

  sub_tau_a <- ref_df[ref_df$fixture == "hand_n3_att4_lone45_a" & ref_df$key == "tau", , drop = FALSE]
  sub_tau_blk <- ref_df[ref_df$fixture == "hand_n3_att4_lone45_a" & ref_df$key == "tau_blockwise", , drop = FALSE]
  expect_equal(res_joint_a$targets$estimate, as.double(sub_tau_a$double), tolerance = 1e-14)
  expect_equal(res_bdiag_a$targets$estimate, as.double(sub_tau_blk$double), tolerance = 1e-14)
  # When V_12 != 0 and r_1 != 0, joint estimate of lone_45_contrast (27/4 = 6.75) differs from blockwise (7.0)
  # AND its model variance shrinks strictly below V[4,4] = 2.0!
  expect_equal(res_joint_a$targets$estimate[4L], 6.75, tolerance = 1e-14)
  expect_equal(res_bdiag_a$targets$estimate[4L], 7.0, tolerance = 1e-14)
  expect_lt(res_joint_a$raw_Z1tZ1[4L, 4L], res_bdiag_a$raw_Z1tZ1[4L, 4L] - 0.1)

  # Case B: W_21 != 0 (`transfer_w21_nonzero`)
  fx_b <- build_fixture_inputs("transfer_w21_nonzero")
  res_joint_b <- solve_chained_gmm_core(
    design = fx_b$design,
    mhat = fx_b$mhat,
    V = fx_b$V,
    provenance = make_default_provenance(),
    targets = fx_b$targets,
    weighting = "optimal",
    if_replay = make_unit_if_replay_from_V(fx_b$V),
    accuracy = default_test_accuracy("optimal"),
    bounds = default_test_bounds(),
    gauge_rule = "lowest"
  )
  sub_tau_b <- ref_df[ref_df$fixture == "transfer_w21_nonzero" & ref_df$key == "tau", , drop = FALSE]
  expect_equal(res_joint_b$targets$estimate, as.double(sub_tau_b$double), tolerance = 1e-14)

  # Compare against no-covariance-transfer formula p3_no_cov_transfer = m4 + p2_blk1
  p3_no_cov <- as.double(ref_df$double[ref_df$fixture == "transfer_w21_nonzero" & ref_df$key == "p3_no_cov_transfer"])
  v21_y1 <- as.double(ref_df$double[ref_df$fixture == "transfer_w21_nonzero" & ref_df$key == "v21_y1"])
  expect_equal(res_joint_b$targets$estimate[3L], p3_no_cov - v21_y1, tolerance = 1e-14)
  expect_gt(abs(res_joint_b$targets$estimate[3L] - p3_no_cov), 0.1)
})

test_that("3. Gauge invariance holds across lowest, highest, custom, and centered minimum-norm gauges", {
  fx <- build_fixture_inputs("spd_multi_cohort_opt")
  replay <- make_unit_if_replay_from_V(fx$V)
  prov <- make_default_provenance()
  bounds <- default_test_bounds(allow_materialize_if = TRUE, max_materialized_if_entries = 5000L)
  acc <- default_test_accuracy("optimal")

  res_low <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets, weighting = "optimal", if_replay = replay,
    accuracy = acc, bounds = bounds, gauge_rule = "lowest",
    materialize_target_if = TRUE
  )
  res_high <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets, weighting = "optimal", if_replay = replay,
    accuracy = acc, bounds = bounds, gauge_rule = "highest",
    materialize_target_if = TRUE
  )
  # Unanchored components in spd_multi_cohort_opt are {4} (isolated) and {5, 6}:
  # pick custom roots c(4L, 6L)
  res_cust <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets, weighting = "optimal", if_replay = replay,
    accuracy = acc, bounds = bounds, gauge_rule = "custom",
    custom_gauge_roots = c(4L, 6L),
    materialize_target_if = TRUE
  )

  # Estimable target estimates, SEs, loadings Lambda, projected IFs Psi, KKT multiplier y, and J_stat are gauge-invariant
  expect_equal(res_low$targets$estimate, res_high$targets$estimate, tolerance = 1e-14)
  expect_equal(res_low$targets$estimate, res_cust$targets$estimate, tolerance = 1e-14)
  expect_equal(res_low$targets$se_sandwich, res_high$targets$se_sandwich, tolerance = 1e-14)
  expect_equal(res_low$targets$se_criterion_conditional, res_high$targets$se_criterion_conditional, tolerance = 1e-14)
  expect_equal(res_low$internal_gauge$Lambda, res_high$internal_gauge$Lambda, tolerance = 1e-14)
  expect_equal(res_low$internal_gauge$Lambda, res_cust$internal_gauge$Lambda, tolerance = 1e-14)
  expect_equal(res_low$if_projection$materialized_target_if, res_high$if_projection$materialized_target_if, tolerance = 1e-14)
  expect_equal(res_low$if_projection$materialized_target_if, res_cust$if_projection$materialized_target_if, tolerance = 1e-14)
  expect_equal(res_low$internal_gauge$y_multiplier, res_high$internal_gauge$y_multiplier, tolerance = 1e-14)
  expect_equal(res_low$raw_quadratic_residual, res_high$raw_quadratic_residual, tolerance = 1e-14)

  # Internal gauge coordinates differ on unanchored components, but map to the identical centered minimum-norm vector
  expect_gt(max(abs(res_low$internal_gauge$theta_gauge_coords - res_high$internal_gauge$theta_gauge_coords)), 0.1)
  cent_low <- transform_gmm_gauge_coordinates(res_low$internal_gauge$gauge, res_low$internal_gauge$theta_gauge_coords, "centered")
  cent_high <- transform_gmm_gauge_coordinates(res_high$internal_gauge$gauge, res_high$internal_gauge$theta_gauge_coords, "centered")
  cent_cust <- transform_gmm_gauge_coordinates(res_cust$internal_gauge$gauge, res_cust$internal_gauge$theta_gauge_coords, "centered")
  expect_equal(cent_low, cent_high, tolerance = 1e-14)
  expect_equal(cent_low, cent_cust, tolerance = 1e-14)

  # Sparse BFS path witnesses satisfy W' w_j = c_j for every estimable target
  W_full <- matrix(0, nrow = fx$design$incidence$n_links, ncol = fx$design$incidence$n_parameters)
  for (l in seq_len(fx$design$incidence$n_links)) {
    if (fx$design$incidence$start[l] > 0L) W_full[l, fx$design$incidence$start[l]] <- -1
    if (fx$design$incidence$end[l] > 0L) W_full[l, fx$design$incidence$end[l]] <- 1
  }
  for (j in seq_along(fx$targets)) {
    wit <- res_low$witnesses[[j]]
    w_vec <- numeric(nrow(W_full))
    w_vec[wit$link_index] <- wit$weight
    c_expected <- numeric(ncol(W_full))
    c_expected[fx$targets[[j]]$param_index] <- fx$targets[[j]]$coeff
    expect_equal(drop(crossprod(W_full, w_vec)), c_expected, tolerance = 1e-15)
  }
})

test_that("4. All 5 fresh exact Fraction KKT reference fixtures match solve_chained_gmm_core to machine precision", {
  ref_df <- utils::read.csv(find_gmm_ref_file("fresh_exact_kkt_reference.csv"), stringsAsFactors = FALSE)
  fixtures <- unique(ref_df$fixture)
  expect_length(fixtures, 5L)

  for (fx_name in fixtures) {
    fx <- build_fixture_inputs(fx_name)
    wt <- if (fx_name == "spd_multi_cohort_id") "identity" else "optimal"
    res <- solve_chained_gmm_core(
      design = fx$design,
      mhat = fx$mhat,
      V = fx$V,
      provenance = make_default_provenance(),
      targets = fx$targets,
      weighting = wt,
      if_replay = make_unit_if_replay_from_V(fx$V),
      accuracy = default_test_accuracy(wt),
      bounds = default_test_bounds(),
      gauge_rule = "lowest"
    )
    sub <- ref_df[ref_df$fixture == fx_name, , drop = FALSE]

    th_ref <- as.double(sub$double[sub$key == "theta"])
    expect_equal(unname(res$internal_gauge$theta_gauge_coords), th_ref, tolerance = 1e-13)

    tau_ref <- as.double(sub$double[sub$key == "tau"])
    est_idx <- which(res$targets$estimable)
    expect_equal(res$targets$estimate[est_idx], tau_ref, tolerance = 1e-13)

    y_ref <- as.double(sub$double[sub$key == "y"])
    expect_equal(unname(res$internal_gauge$y_multiplier), y_ref, tolerance = 1e-13)

    lam_sub <- sub[sub$key == "Lambda", , drop = FALSE]
    for (r in seq_len(nrow(lam_sub))) {
      expect_equal(
        res$internal_gauge$Lambda[lam_sub$row[r], lam_sub$col[r]],
        as.double(lam_sub$double[r]),
        tolerance = 1e-13
      )
    }

    mc_sub <- sub[sub$key == "Mc", , drop = FALSE]
    for (r in seq_len(nrow(mc_sub))) {
      expect_equal(
        res$internal_gauge$Mc[mc_sub$row[r], mc_sub$col[r]],
        as.double(mc_sub$double[r]),
        tolerance = 1e-13
      )
    }

    mcov_sub <- sub[sub$key == "model_cov", , drop = FALSE]
    for (r in seq_len(nrow(mcov_sub))) {
      expect_equal(
        res$raw_Z1tZ1[mcov_sub$row[r], mcov_sub$col[r]],
        as.double(mcov_sub$double[r]),
        tolerance = 1e-13
      )
    }

    j_ref <- as.double(sub$double[sub$key == "J_stat"])
    expect_equal(res$raw_quadratic_residual, j_ref, tolerance = 1e-13)
  }
})

test_that("5. Identity vs optimal weighting differ, identity uses sandwich SEs, and never reduces or auto-falls back", {
  fx <- build_fixture_inputs("spd_multi_cohort_opt")
  replay <- make_unit_if_replay_from_V(fx$V)
  prov <- make_default_provenance()

  res_opt <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets, weighting = "optimal", if_replay = replay,
    accuracy = default_test_accuracy("optimal"), bounds = default_test_bounds(),
    gauge_rule = "lowest"
  )
  res_id <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets, weighting = "identity", if_replay = replay,
    accuracy = default_test_accuracy("identity"), bounds = default_test_bounds(),
    gauge_rule = "lowest"
  )

  # Point estimates and loadings differ between optimal and identity weighting
  expect_gt(max(abs(res_opt$targets$estimate - res_id$targets$estimate)), 0.05)
  expect_gt(max(abs(res_opt$internal_gauge$Lambda - res_id$internal_gauge$Lambda)), 0.05)

  # Under identity weighting: model_cov is NULL, se_model is NA, J_stat_optimal is NA,
  # while sandwich_cov equals Lambda_id' V Lambda_id!
  expect_null(res_id$model_cov)
  expect_true(all(is.na(res_id$targets$se_model)))
  expect_true(is.na(res_id$J_stat_optimal))
  expect_false(res_id$scale_matches_inference)
  expected_id_sand <- t(res_id$internal_gauge$Lambda) %*% fx$V %*% res_id$internal_gauge$Lambda
  expect_equal(res_id$sandwich_cov, expected_id_sand, tolerance = 1e-13)
  expect_equal(res_id$targets$se_sandwich, unname(sqrt(diag(expected_id_sand))), tolerance = 1e-13)
  # Optimal SEs are strictly smaller than identity sandwich SEs on overidentified targets
  expect_lt(res_opt$targets$se_sandwich[1L], res_id$targets$se_sandwich[1L])

  # No automatic fallback to identity when optimal weighting fails PD/Cholesky
  V_bad <- fx$V
  V_bad[1, 1] <- -1.0
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = V_bad, provenance = prov,
      targets = fx$targets, weighting = "optimal", if_replay = replay,
      accuracy = default_test_accuracy("optimal"), bounds = default_test_bounds(),
      gauge_rule = "lowest"
    ),
    class = "cdid_gmm_pd_error"
  )

  # Reduction is refused under identity weighting
  prov_red <- make_default_provenance(
    reduction = list(
      mode = "pre_reduced",
      context_id = "ctx1",
      cohort = "2",
      operator_id = "op_dr",
      vertex_map_id = "vmap1",
      risk_set_pair_vector_id = "pair_s1",
      shared_nuisance_object_id = "nuis_shared_1",
      retained_link_ids = fx$design$incidence$link_ids,
      original_n_links = 7L,
      certified_by = "structural_key_equality"
    )
  )
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov_red,
      targets = fx$targets, weighting = "identity", if_replay = replay,
      accuracy = default_test_accuracy("identity"), bounds = default_test_bounds(),
      gauge_rule = "lowest"
    ),
    class = "cdid_gmm_reduction_error"
  )
})

test_that("6. Streamed IF projection and weighted cluster centering across batch boundaries match unbatched oracle and refuse invalid streams", {
  fx <- build_fixture_inputs("hand_n3_att4_lone45_a")
  # 11 units across 4 clusters (sizes 3, 4, 2, 2) with cluster weights
  cluster_ids <- c("cl_A", "cl_A", "cl_A", "cl_B", "cl_B", "cl_B", "cl_B", "cl_C", "cl_C", "cl_D", "cl_D")
  unit_ids <- sprintf("unit_%02d", 1:11)
  cluster_w_map <- c(cl_A = 1.5, cl_B = 0.8, cl_C = 2.2, cl_D = 1.1)
  unit_w <- unname(cluster_w_map[cluster_ids])

  # 5 clusters so M - 1 = 4 >= m_R = 4!
  cluster_ids <- c("cl_A", "cl_A", "cl_A", "cl_B", "cl_B", "cl_B", "cl_C", "cl_C", "cl_D", "cl_D", "cl_E")
  cluster_w_map <- c(cl_A = 1.5, cl_B = 0.8, cl_C = 2.2, cl_D = 1.1, cl_E = 1.3)
  unit_w <- unname(cluster_w_map[cluster_ids])

  Phi_units <- matrix(c(
     1.2, -0.4,  0.5,  0.3,
    -0.3,  0.8, -0.2,  0.1,
     0.5,  0.2,  0.7, -0.6,
    -0.9, -0.5,  0.4,  0.8,
     0.4,  1.1, -0.8, -0.2,
     0.2, -0.7,  0.3,  0.5,
    -0.6,  0.3,  0.9, -0.4,
     1.0, -0.9, -0.5,  0.7,
    -0.4,  0.6,  0.2, -0.9,
    -0.8,  0.4, -0.6,  0.2,
     0.7, -0.2,  0.1, -0.5
  ), nrow = 11L, byrow = TRUE)

  # Exact cluster-aggregated weighted centered covariance in moment space
  u_cl <- c("cl_A", "cl_B", "cl_C", "cl_D", "cl_E")
  Phi_cl <- rowsum(Phi_units, cluster_ids, reorder = FALSE)[u_cl, , drop = FALSE]
  w_cl <- unname(cluster_w_map[u_cl])
  mu_cl <- colSums(Phi_cl * w_cl) / sum(w_cl)
  Phi_cl_cent <- sweep(Phi_cl, 2L, mu_cl, `-`)
  alpha_inf <- 1 / (length(u_cl) - 1)
  V_cl <- alpha_inf * (t(Phi_cl_cent) %*% (Phi_cl_cent * w_cl))
  V_cl <- 0.5 * (V_cl + t(V_cl))

  prov_cl <- gmm_covariance_provenance(
    row_level = "cluster",
    n_units = 11L,
    n_clusters = 5L,
    n_positive_contributing_rows = 5L,
    centered = TRUE,
    weighted = TRUE,
    weight_type = "cluster_multiplier",
    normalization_factor = alpha_inf,
    scale_id = "cluster_centered_w",
    nuisance_fit_id = "nuis_cl_1"
  )

  # Test across batch sizes B_b = 2 (clusters split across multiple batch boundaries!), 3, and 11
  for (b_sz in c(2L, 3L, 11L)) {
    cb_cl <- function(consume_batch) {
      starts <- seq(1L, 11L, by = b_sz)
      for (st in starts) {
        en <- min(st + b_sz - 1L, 11L)
        consume_batch(list(
          unit_id = unit_ids[st:en],
          cluster_id = cluster_ids[st:en],
          phi = Phi_units[st:en, , drop = FALSE],
          weight = unit_w[st:en],
          nuisance_fit_id = "nuis_cl_1"
        ))
      }
    }
    rep_cl <- gmm_if_replay_spec(
      callback = cb_cl,
      input_row_level = "unit",
      aggregate_to_clusters = TRUE,
      n_input_rows = 11L,
      n_output_rows = 5L,
      n_positive_contributing_rows = 5L,
      centered = TRUE,
      weighted = TRUE,
      weight_type = "cluster_multiplier",
      normalization_factor = alpha_inf,
      scale_id = "cluster_centered_w",
      nuisance_fit_id = "nuis_cl_1",
      expected_unit_ids = unit_ids,
      expected_cluster_ids = u_cl
    )
    res_cl <- solve_chained_gmm_core(
      design = fx$design,
      mhat = fx$mhat,
      V = V_cl,
      provenance = prov_cl,
      targets = fx$targets[1:4],
      weighting = "optimal",
      if_replay = rep_cl,
      accuracy = default_test_accuracy("optimal"),
      bounds = default_test_bounds(max_batch_rows = b_sz),
      gauge_rule = "lowest"
    )
    Lam <- res_cl$internal_gauge$Lambda
    expected_sand <- t(Lam) %*% V_cl %*% Lam
    expect_equal(res_cl$sandwich_cov, expected_sand, tolerance = 1e-14)
    # Actual covariance equality is proved by this fixture, not by the metadata.
    expect_equal(res_cl$raw_Z1tZ1, expected_sand, tolerance = 1e-13)
    expect_null(res_cl$model_cov)
    expect_lt(res_cl$variance_gap$abs, 1e-13)
  }

  # Refusal when cluster reappears after closure (unsorted clusters)
  bad_cl_ids <- c("cl_A", "cl_B", "cl_A", "cl_B", "cl_B", "cl_B", "cl_C", "cl_C", "cl_D", "cl_D", "cl_E")
  rep_unsorted <- gmm_if_replay_spec(
    callback = function(consume_batch) {
      consume_batch(list(unit_id = unit_ids[1:2], cluster_id = bad_cl_ids[1:2], phi = Phi_units[1:2, ], weight = c(1.5, 0.8)))
      consume_batch(list(unit_id = unit_ids[3:11], cluster_id = bad_cl_ids[3:11], phi = Phi_units[3:11, ], weight = unit_w[3:11]))
    },
    input_row_level = "unit",
    aggregate_to_clusters = TRUE,
    n_input_rows = 11L,
    n_output_rows = 5L,
    n_positive_contributing_rows = 5L,
    centered = TRUE,
    weighted = TRUE,
    weight_type = "cluster_multiplier",
    normalization_factor = alpha_inf,
    scale_id = "cluster_centered_w",
    nuisance_fit_id = "nuis_cl_1"
  )
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = V_cl, provenance = prov_cl,
      targets = fx$targets[1:4], weighting = "optimal", if_replay = rep_unsorted,
      accuracy = default_test_accuracy("optimal"), bounds = default_test_bounds(),
      gauge_rule = "lowest"
    ),
    class = "cdid_gmm_replay_error"
  )

  # Refusal when input_row_level == "cluster" and aggregate_to_clusters == TRUE (double summation)
  expect_error(
    gmm_if_replay_spec(
      callback = function(consume_batch) NULL,
      input_row_level = "cluster",
      aggregate_to_clusters = TRUE,
      n_input_rows = 5L,
      n_output_rows = 5L,
      n_positive_contributing_rows = 5L,
      centered = TRUE,
      weighted = TRUE,
      weight_type = "cluster_multiplier",
      normalization_factor = alpha_inf,
      scale_id = "cluster_centered_w",
      nuisance_fit_id = "nuis_cl_1"
    ),
    class = "cdid_gmm_replay_error"
  )
})

test_that("7. Stored indefinite (stored_6c_V) and truly singular (singular_V) matrices pass base::chol but are REFUSED by numerical guard", {
  for (mat_name in c("stored_6c_V", "singular_V")) {
    V_mat <- read_stored_3x3_matrix(mat_name)

    # 1. Confirm base::chol succeeds with positive diagonal pivots!
    R_floating <- base::chol(V_mat)
    expect_true(is.matrix(R_floating))
    expect_true(all(diag(R_floating) > 0))
    expect_equal(min(diag(R_floating)), 2^-30)

    # 2. Test on 3-link design: both must be REFUSED by solve_chained_gmm_core
    links_3 <- data.frame(
      g = c(2, 2, 2),
      s = c(0, 0, 1),
      t = c(1, 2, 2),
      ref = c(1, 1, 1),
      functional_key = rep("k1", 3L),
      stringsAsFactors = FALSE
    )
    des_3 <- construct_chained_design(links_3)
    targets_3 <- list(gmm_target_spec("p1", 1L, 1.0), gmm_target_spec("p2", 2L, 1.0))
    replay_3 <- gmm_if_replay_spec(
      callback = function(consume_batch) {
        consume_batch(list(unit_id = sprintf("u%d", 1:10), phi = matrix(0.1, nrow = 10L, ncol = 3L), weight = rep(1, 10L)))
      },
      input_row_level = "unit",
      aggregate_to_clusters = FALSE,
      n_input_rows = 10L,
      n_output_rows = 10L,
      n_positive_contributing_rows = 10L,
      centered = TRUE,
      weighted = TRUE,
      weight_type = "unit_uniform",
      normalization_factor = 1 / 9,
      scale_id = "s1",
      nuisance_fit_id = "n1"
    )
    prov_3 <- make_default_provenance(n_units = 10L, scale_id = "s1", nuisance_fit_id = "n1")

    err <- expect_error(
      solve_chained_gmm_core(
        design = des_3,
        mhat = c(-0.2, 0.4, 0.6),
        V = V_mat,
        provenance = prov_3,
        targets = targets_3,
        weighting = "optimal",
        if_replay = replay_3,
        accuracy = default_test_accuracy("optimal"),
        bounds = default_test_bounds(),
        gauge_rule = "lowest"
      ),
      class = "cdid_gmm_pd_error"
    )
    expect_true("cholesky_factor_residual_exceeds_certified_min_eigenvalue" %in% err$diagnostics$reasons)
    expect_lte(err$diagnostics$lambda_min_cert, 0)
  }
})

test_that("8. Missing/invalid controls, bounds, overflow, scale mismatch, untrusted metadata, and uncertified reductions refuse", {
  fx <- build_fixture_inputs("hand_n3_att4_lone45_a")
  replay <- make_unit_if_replay_from_V(fx$V)
  prov <- make_default_provenance()
  acc <- default_test_accuracy("optimal")
  bnds <- default_test_bounds()

  # Missing required control field
  acc_missing <- acc
  acc_missing$max_kkt_residual_rel <- NULL
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
      targets = fx$targets, weighting = "optimal", if_replay = replay,
      accuracy = acc_missing, bounds = bnds, gauge_rule = "lowest"
    ),
    class = "cdid_gmm_control_error"
  )

  # Uncertified method name (chol_only / rcond_only / ordinary_indicator) refused
  acc_bad_method <- acc
  acc_bad_method$method <- "chol_and_rcond_only"
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
      targets = fx$targets, weighting = "optimal", if_replay = replay,
      accuracy = acc_bad_method, bounds = bnds, gauge_rule = "lowest"
    ),
    class = "cdid_gmm_accuracy_error"
  )

  # Untrusted metadata PD flags on V or provenance refused
  V_flagged <- fx$V
  attr(V_flagged, "is_pd") <- TRUE
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = V_flagged, provenance = prov,
      targets = fx$targets, weighting = "optimal", if_replay = replay,
      accuracy = acc, bounds = bnds, gauge_rule = "lowest"
    ),
    class = "cdid_gmm_covariance_error"
  )
  prov_flagged <- prov
  prov_flagged$certified_pd <- TRUE
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov_flagged,
      targets = fx$targets, weighting = "optimal", if_replay = replay,
      accuracy = acc, bounds = bnds, gauge_rule = "lowest"
    ),
    class = "cdid_gmm_covariance_error"
  )

  # Sample rank bound m_R > M - 1 refused before Cholesky
  prov_rank_low <- make_default_provenance(n_units = 4L) # centered = TRUE => M - 1 = 3 < m_R = 4
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov_rank_low,
      targets = fx$targets, weighting = "optimal", if_replay = replay,
      accuracy = acc, bounds = bnds, gauge_rule = "lowest"
    ),
    class = "cdid_gmm_rank_bound_error"
  )

  # Integer product overflow and preflight budget refusal
  expect_error(gmm_safe_product(1e9, 1e9, 1e15, "overflow_check"), class = "cdid_gmm_bound_error")
  bnds_tiny <- default_test_bounds(max_matrix_entries = 10L)
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
      targets = fx$targets, weighting = "optimal", if_replay = replay,
      accuracy = acc, bounds = bnds_tiny, gauge_rule = "lowest"
    ),
    class = "cdid_gmm_bound_error"
  )

  # Scale mismatch between V provenance and IF replay
  prov_diff_scale <- prov
  prov_diff_scale$normalization_factor <- 1.0 # differs from 1/19
  expect_error(
    solve_chained_gmm_core(
      design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov_diff_scale,
      targets = fx$targets[1:4], weighting = "optimal", if_replay = replay,
      accuracy = default_test_accuracy("optimal", target_error_unit = "abs", require_matching_scale = TRUE),
      bounds = bnds, gauge_rule = "lowest"
    ),
    class = "cdid_gmm_scale_error"
  )
  res_unmatched <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov_diff_scale,
    targets = fx$targets[1:4], weighting = "optimal", if_replay = replay,
    accuracy = default_test_accuracy("optimal", target_error_unit = "abs", require_matching_scale = FALSE),
    bounds = bnds, gauge_rule = "lowest"
  )
  expect_false(res_unmatched$scale_matches_inference)
  expect_null(res_unmatched$model_cov)
  expect_true(all(is.na(res_unmatched$targets$se_model)))
  expect_true(is.na(res_unmatched$J_stat_optimal))
  expect_true(all(is.finite(res_unmatched$targets$se_sandwich)))

  # Unimplemented IPW auto-detection and numeric-agreement reduction modes refuse
  for (bad_mode in c("auto_detect_ipw", "numeric_agreement", "class_label_only", "sample_rank")) {
    p_bad <- make_default_provenance(reduction = list(mode = bad_mode))
    expect_error(
      solve_chained_gmm_core(
        design = fx$design, mhat = fx$mhat, V = fx$V, provenance = p_bad,
        targets = fx$targets[1:4], weighting = "optimal", if_replay = replay,
        accuracy = acc, bounds = bnds, gauge_rule = "lowest"
      ),
      class = "cdid_gmm_reduction_error"
    )
  }

  # Spanning forest helper restarts BFS across disconnected date components within one certified class
  cat_disc <- data.frame(
    s = c(1, 1, 2, 10, 10, 11),
    t = c(2, 3, 3, 11, 12, 12),
    context_id = rep("ctx1", 6L),
    g = rep(20, 6L),
    operator_id = rep("dr_att", 6L),
    vertex_map_id = rep("vmap1", 6L),
    pair_s_id = rep("ps_shared", 6L),
    nuisance_object_id = rep("nuis_shared", 6L),
    positive_rates = rep(TRUE, 6L),
    regular_logit = rep(TRUE, 6L),
    stringsAsFactors = FALSE
  )
  # Conditional graph algebra requires an explicit bound; metadata is no class proof.
  forest <- build_certified_class_forest(cat_disc, weighting = "optimal", max_matrix_entries = 50000L)
  expect_identical(forest$retained_links, c(1L, 2L, 4L, 5L))
  expect_equal(unname(forest$E[3L, ]), c(-1L, 1L, 0L, 0L))
  expect_equal(unname(forest$E[6L, ]), c(0L, 0L, -1L, 1L))
})

test_that("9. Exact IEEE 754 binary64 dyadic coefficient sum handles extreme exponents without false zero cancellation", {
  # Pairwise binary64 addition ((2^100 + 0.1) - 2^100) equals 0, while exact sum is +0.1
  x_cat <- c(2^100, 0.1, -2^100)
  expect_equal((x_cat[1] + x_cat[2]) + x_cat[3], 0)
  ds_cat <- exact_binary64_dyadic_sum(x_cat)
  expect_false(ds_cat$is_zero)
  expect_identical(ds_cat$sign, 1L)

  # Extreme exponent span: 2^1020 and subnormal 2^-1074 (where even 80-bit extended precision sum() is 0!)
  min_subnormal <- hex16_to_double("0000000000000001") # 2^-1074
  x_ext_pos <- c(2^1020, min_subnormal, -2^1020)
  expect_equal(sum(x_ext_pos), 0)
  ds_ext_pos <- exact_binary64_dyadic_sum(x_ext_pos)
  expect_false(ds_ext_pos$is_zero)
  expect_identical(ds_ext_pos$sign, 1L)

  x_ext_neg <- c(2^1020, -min_subnormal, -2^1020)
  expect_equal(sum(x_ext_neg), 0)
  ds_ext_neg <- exact_binary64_dyadic_sum(x_ext_neg)
  expect_false(ds_ext_neg$is_zero)
  expect_identical(ds_ext_neg$sign, -1L)

  # Exact zero across extreme exponents and permutations
  x_zero <- c(2^1020, min_subnormal, 0.1, -2^1020, -0.1, -min_subnormal)
  expect_true(exact_binary64_dyadic_sum(x_zero)$is_zero)
  expect_true(exact_binary64_dyadic_sum(rev(x_zero))$is_zero)

  # Verify target estimability integration on unanchored component {4, 5} of hand_n3_att4_lone45_a
  fx <- build_fixture_inputs("hand_n3_att4_lone45_a")
  t_extreme <- list(
    # Sum on unanchored {4, 5} is 2^100 + 0.1 - 2^100 = 0.1 != 0 -> MUST be non-estimable (NA)
    gmm_target_spec("false_zero_unanchored", c(4L, 4L, 5L), c(2^100, 0.1, -2^100)),
    # Exact zero sum on {4, 5}: 2^100 - 2^100 + (-1 + 1) = 0 -> estimable
    gmm_target_spec("true_zero_unanchored", c(4L, 4L, 4L, 5L), c(2^100, -2^100, -1.0, 1.0))
  )
  res_ext <- solve_chained_gmm_core(
    design = fx$design,
    mhat = fx$mhat,
    V = fx$V,
    provenance = make_default_provenance(),
    targets = t_extreme,
    weighting = "optimal",
    if_replay = make_unit_if_replay_from_V(fx$V),
    accuracy = default_test_accuracy("optimal"),
    bounds = default_test_bounds(),
    gauge_rule = "lowest"
  )
  expect_identical(res_ext$targets$estimable, c(FALSE, TRUE))
  expect_identical(res_ext$targets$status, c("disconnected_component", "estimable"))
  expect_true(is.na(res_ext$targets$estimate[1L]))
  expect_equal(res_ext$targets$estimate[2L], 6.75, tolerance = 1e-14)
})

test_that("10. Deliberate mutations (wrong whitening transpose, wrong KKT sign, missing adjoint residual term) are caught and refused", {
  fx <- build_fixture_inputs("hand_n3_att4_lone45_a")
  replay <- make_unit_if_replay_from_V(fx$V)
  prov <- make_default_provenance()
  acc <- default_test_accuracy("optimal")
  bnds <- default_test_bounds()

  res_ok <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets[1:4], weighting = "optimal", if_replay = replay,
    accuracy = acc, bounds = bnds, gauge_rule = "lowest"
  )
  ig <- res_ok$internal_gauge
  W_G <- materialize_gmm_gauge_matrix(ig$gauge, bnds$max_matrix_entries)
  A_ok <- backsolve(ig$R_chol, W_G, transpose = TRUE)
  t_info <- validate_gmm_targets(fx$design, ig$gauge, fx$targets[1:4], bnds)

  # Mutation A: Wrong whitening transpose R^{-1} W_G instead of R^{-T} W_G
  A_wrong <- backsolve(ig$R_chol, W_G, transpose = FALSE)
  expect_error(
    certify_gmm_numerical_guard(
      V_R = fx$V, R_chol = ig$R_chol, W_G = W_G, A_whitened = A_wrong,
      U_qr = ig$U_qr, piv = ig$pivot, mhat_R = fx$mhat,
      theta = ig$theta_gauge_coords, y = ig$y_multiplier,
      C_gauge = t_info$C_gauge, Lambda = ig$Lambda, Mc = ig$Mc, Z1 = ig$Z1,
      weighting = "optimal", accuracy = acc
    ),
    class = "cdid_gmm_qr_error"
  )

  # Mutation B: Wrong KKT multiplier sign (y_wrong = -y)
  expect_error(
    certify_gmm_numerical_guard(
      V_R = fx$V, R_chol = ig$R_chol, W_G = W_G, A_whitened = A_ok,
      U_qr = ig$U_qr, piv = ig$pivot, mhat_R = fx$mhat,
      theta = ig$theta_gauge_coords, y = -ig$y_multiplier,
      C_gauge = t_info$C_gauge, Lambda = ig$Lambda, Mc = ig$Mc, Z1 = ig$Z1,
      weighting = "optimal", accuracy = acc
    ),
    class = "cdid_gmm_accuracy_error"
  )

  # Mutation C: Perturbed approximate adjoint and primal solution where ordinary indicator
  # |Lambda_1' r1| + |Mc_1' r2| == 0 by orthogonality, while actual target error c_1' delta_theta != 0
  # is carried entirely by the missing term -s_1' K^{-1} r!
  # Perturb Lambda[, 1] slightly so g1 = W_G' Lambda_bad[, 1] = c_1 + (0, 1e-2, 0),
  # and choose delta_theta = (1e-4, -1e-2, 0) so g1' delta_theta = 1e-4 - 1e-4 = 0
  # (making Lambda_bad[, 1]' r1 = -g1' delta_theta = 0), while actual error c_1' delta_theta = 1e-4!
  Lambda_bad <- ig$Lambda
  # Link 1 has row [0, 1, 0] in W_G, so adding 1e-2 to Lambda_bad[1, 1] shifts W_G' Lambda_bad[, 1] by (0, 1e-2, 0)
  Lambda_bad[1L, 1L] <- Lambda_bad[1L, 1L] + 1e-2
  delta_th <- c(1e-4, -1e-2, 0)
  theta_pert <- ig$theta_gauge_coords + delta_th

  acc_relax_res <- acc
  acc_relax_res$max_kkt_residual_rel <- 0.05
  acc_relax_res$max_loading_residual_rel <- 0.05
  acc_relax_res$target_error_unit <- "abs"
  acc_relax_res$max_target_error_abs <- 1e-6

  err_mut_c <- expect_error(
    certify_gmm_numerical_guard(
      V_R = fx$V, R_chol = ig$R_chol, W_G = W_G, A_whitened = A_ok,
      U_qr = ig$U_qr, piv = ig$pivot, mhat_R = fx$mhat,
      theta = theta_pert, y = ig$y_multiplier,
      C_gauge = t_info$C_gauge, Lambda = Lambda_bad, Mc = ig$Mc, Z1 = ig$Z1,
      weighting = "optimal", accuracy = acc_relax_res
    ),
    class = "cdid_gmm_accuracy_error"
  )
  # Confirm that the ordinary computed indicator on target 1 was ~0 (< 1e-14, well below budget 1e-6)
  # while actual error is 1e-4 and target_error_cert (via adjoint_missing_term_cert) caught it!
  expect_lt(err_mut_c$diagnostics$target_error_est[1L], 1e-14)
  expect_gt(err_mut_c$diagnostics$adjoint_missing_term_cert[1L], 1e-4)
  expect_gt(err_mut_c$diagnostics$target_error_cert[1L], 1e-4)
})

# D12 addition: verify conditional point acceptance, outward G_lambda export, explicit loading-error
# budgeting (NULL absence state vs explicit budget), and separation from pending complete inference.
test_that("11. D12 guard exports outward G_lambda, separates point/loading/complete inference, and accepts zero-dual exact-fit candidates", {
  fx <- build_fixture_inputs("hand_n3_att4_lone45_a")
  replay <- make_unit_if_replay_from_V(fx$V)
  prov <- make_default_provenance()
  bnds <- default_test_bounds()

  # Default NULL max_loading_error_abs means not supplied (absence state):
  # point is certified when G_c <= budget, while loading certification remains incomplete.
  acc_no_load <- default_test_accuracy("optimal", max_loading_error_abs = NULL)
  res_no_load <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets, weighting = "optimal", if_replay = replay,
    accuracy = acc_no_load, bounds = bnds, gauge_rule = "lowest"
  )
  expect_true(res_no_load$point_certified)
  expect_false(res_no_load$loading_certified)
  expect_false(res_no_load$loading_budget_supplied)
  expect_identical(res_no_load$loading_validation_status, "incomplete_missing_caller_budget")
  expect_false(res_no_load$complete_inference_certified)
  expect_identical(res_no_load$targets$complete_inference_certified, rep(FALSE, nrow(res_no_load$targets)))
  expect_identical(
    res_no_load$full_inference_status,
    "pending_unvalidated_raw_root_covariance_projection_multiplier"
  )
  est_idx <- which(res_no_load$targets$estimable)
  unest_idx <- which(!res_no_load$targets$estimable)
  expect_true(all(is.finite(res_no_load$targets$G_lambda[est_idx])))
  expect_identical(res_no_load$targets$G_lambda, res_no_load$targets$loading_error_cert)
  expect_true(all(is.na(res_no_load$targets$G_lambda[unest_idx])))
  expect_true(all(res_no_load$targets$point_certified[est_idx]))
  expect_true(all(!res_no_load$targets$loading_certified[est_idx]))
  expect_true(all(!res_no_load$targets$complete_inference_certified[est_idx]))

  # Supplying an explicit caller loading budget certifies stored-array loadings (still pending full inference)
  acc_with_load <- default_test_accuracy("optimal", max_loading_error_abs = 1e-10)
  res_with_load <- solve_chained_gmm_core(
    design = fx$design, mhat = fx$mhat, V = fx$V, provenance = prov,
    targets = fx$targets, weighting = "optimal", if_replay = replay,
    accuracy = acc_with_load, bounds = bnds, gauge_rule = "lowest"
  )
  expect_true(res_with_load$point_certified)
  expect_true(res_with_load$loading_certified)
  expect_true(res_with_load$loading_budget_supplied)
  expect_identical(res_with_load$loading_validation_status, "certified_conditional_stored_array")
  expect_false(res_with_load$complete_inference_certified)
  expect_identical(res_with_load$targets$complete_inference_certified, rep(FALSE, nrow(res_with_load$targets)))

  # Corrupt loading Lambda while keeping primal (theta, y) and point estimate accurate:
  # When y == 0 (exact-fit / zero-dual regime), corrupting Lambda leaves r = 0 and G_c ~ 0,
  # so point-only succeeds when max_loading_error_abs is NULL, but explicit loading budget refuses!
  ig <- res_with_load$internal_gauge
  W_G <- materialize_gmm_gauge_matrix(ig$gauge, bnds$max_matrix_entries)
  A_ok <- backsolve(ig$R_chol, W_G, transpose = TRUE)
  t_info <- validate_gmm_targets(fx$design, ig$gauge, fx$targets[1:4], bnds)
  mhat_exact_fit <- drop(W_G %*% ig$theta_gauge_coords)
  y_zero <- rep(0.0, nrow(W_G))
  Lambda_corrupt <- ig$Lambda
  Lambda_corrupt[1L, 1L] <- Lambda_corrupt[1L, 1L] + 1e-3

  g_pt_only <- certify_gmm_numerical_guard(
    V_R = fx$V, R_chol = ig$R_chol, W_G = W_G, A_whitened = A_ok,
    U_qr = ig$U_qr, piv = ig$pivot, mhat_R = mhat_exact_fit,
    theta = ig$theta_gauge_coords, y = y_zero,
    C_gauge = t_info$C_gauge, Lambda = Lambda_corrupt, Mc = ig$Mc, Z1 = ig$Z1,
    weighting = "optimal", accuracy = acc_no_load
  )
  expect_true(g_pt_only$point_certified)
  expect_false(g_pt_only$loading_certified)
  expect_true(g_pt_only$targets$loading_rel_exceeds_threshold[1L])
  expect_gt(g_pt_only$targets$G_lambda[1L], 1e-3)

  err_load <- expect_error(
    certify_gmm_numerical_guard(
      V_R = fx$V, R_chol = ig$R_chol, W_G = W_G, A_whitened = A_ok,
      U_qr = ig$U_qr, piv = ig$pivot, mhat_R = mhat_exact_fit,
      theta = ig$theta_gauge_coords, y = y_zero,
      C_gauge = t_info$C_gauge, Lambda = Lambda_corrupt, Mc = ig$Mc, Z1 = ig$Z1,
      weighting = "optimal", accuracy = acc_with_load
    ),
    class = "cdid_gmm_accuracy_error"
  )
  expect_true("certified_loading_error_abs_exceeded" %in% err_load$diagnostics$reasons)
})

