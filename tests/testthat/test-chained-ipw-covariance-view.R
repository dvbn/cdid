require_cdid_namespace <- function() {
  if (!isNamespaceLoaded("cdid") && !requireNamespace("cdid", quietly = TRUE)) {
    stop("Actual installed 'cdid' namespace is required; no fallback is permitted in packaged tests.", call. = FALSE)
  }
  ns <- asNamespace("cdid")
  if (!isNamespace(ns) || !environmentIsLocked(ns)) {
    stop("'cdid' namespace must be a locked R namespace.", call. = FALSE)
  }
  ns
}

# Require installed namespace immediately at file load (never skip or fallback)
require_cdid_namespace()

test_that("private compact-view IPW covariance transfer entrypoints exist in locked cdid namespace and are unexported", {
  ns <- require_cdid_namespace()
  expect_true(isNamespace(ns))
  expect_true(environmentIsLocked(ns))

  priv_fns <- c(
    "ipw_cov_transfer_view_require_deps",
    "ipw_cov_transfer_validate_frame_owner_header",
    "ipw_cov_transfer_validate_links_view_headers",
    "ipw_cov_transfer_validate_links_view_content",
    "preflight_ipw_covariance_transfer_view",
    "certify_ipw_covariance_transfer_view"
  )
  exports <- getNamespaceExports("cdid")
  for (fn in priv_fns) {
    expect_true(exists(fn, envir = ns, mode = "function", inherits = FALSE), info = fn)
    expect_true(bindingIsLocked(fn, ns), info = fn)
    expect_identical(environment(get(fn, envir = ns, inherits = FALSE)), ns, info = fn)
    expect_false(fn %in% exports, info = fn)
  }
})

test_that("certify_ipw_covariance_transfer_view matches dense wrapper bitwise on 14-unit 3-link dataset for unit GMM and cluster inference across chunk/block sizes", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))

  for (id_mode in c("character", "numeric")) {
    b14 <- fx$build_14unit_3link_view_and_dense(
      ns, r_ball = 1e-6, chunk_rows = 5L, operator_block_rows = 4L, id_mode = id_mode
    )

    res_v_gmm <- ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view,
      candidate_rows = b14$candidate_rows,
      frame_owner = b14$frame_owner,
      C = b14$C,
      op = b14$op_gmm,
      V_candidate = b14$V_cand_gmm,
      control = b14$cov_ctrl,
      reader_control = b14$o_ctrl
    )
    res_d_gmm <- ns$certify_ipw_covariance_transfer(
      links = b14$links_dense,
      candidate_rows = b14$candidate_rows,
      frame = b14$frame_dense,
      C = b14$C,
      op = b14$op_gmm,
      V_candidate = b14$V_cand_gmm,
      control = b14$cov_ctrl
    )

    expect_identical(res_v_gmm$storage_mode, "compact_view")
    expect_true(res_v_gmm$certified)
    expect_true(res_v_gmm$enclosure_certified)
    expect_true(res_v_gmm$pd_transferred)
    expect_identical(res_v_gmm$status, "CERTIFIED")
    expect_false(res_v_gmm$inconclusive)
    expect_false(res_v_gmm$separation_certified)
    expect_false(res_v_gmm$sampling_certified)
    expect_false(res_v_gmm$identification_certified)
    expect_false(res_v_gmm$raw_panel_eligibility_certified)
    expect_false(res_v_gmm$public_gate_adopted)
    expect_false(res_v_gmm$singularity_proved)
    expect_false(res_v_gmm$nonidentification_proved)
    expect_false(res_v_gmm$rank_deficiency_proved)
    ns$validate_ipw_covariance_transfer(res_v_gmm, require_certified = TRUE)

    # Bitwise numerical equivalence between compact view and dense covariance transfer
    expect_identical(res_v_gmm$enclosure$lo, res_d_gmm$enclosure$lo)
    expect_identical(res_v_gmm$enclosure$hi, res_d_gmm$enclosure$hi)
    expect_identical(res_v_gmm$epsilon_acc, res_d_gmm$epsilon_acc)
    expect_identical(res_v_gmm$nu_eval, res_d_gmm$nu_eval)
    expect_identical(res_v_gmm$nu_root, res_d_gmm$nu_root)
    expect_identical(res_v_gmm$nu_total, res_d_gmm$nu_total)
    expect_identical(res_v_gmm$epsilon_V, res_d_gmm$epsilon_V)
    expect_identical(res_v_gmm$eigen_bounds, res_d_gmm$eigen_bounds)
    expect_identical(res_v_gmm$work$actual_4bd_interval_ops, res_d_gmm$work$actual_4bd_interval_ops)
    expect_identical(res_v_gmm$work$actual_component_interval_ops, res_d_gmm$work$actual_component_interval_ops)
    expect_identical(res_v_gmm$work$total_actual_interval_ops, res_d_gmm$work$total_actual_interval_ops)

    for (j in 1:2) {
      for (k in 1:2) {
        expect_lte(res_v_gmm$enclosure$lo[j, k], b14$ref_exact_V_gmm[j, k])
        expect_gte(res_v_gmm$enclosure$hi[j, k], b14$ref_exact_V_gmm[j, k])
      }
    }

    # Verify scratch reuse across links l1 and l2 sharing prob12_view, and release before l3 (prob3_view)
    expect_identical(res_v_gmm$n_unique_nuisance_problems, 2L)
    expect_identical(res_v_gmm$work$scratch_allocations, 3L)
    expect_identical(res_v_gmm$work$scratch_reuse_count, 3L)
    expect_gt(res_v_gmm$work$reused_shared_4bd_interval_ops, 0)
    expect_gt(res_v_gmm$work$actual_reader_work_ops, 0)
    expect_gt(res_v_gmm$work$cache_copy_work_ops_charged, 0)
    expect_true(res_v_gmm$lifetime$callback_aliases_released_before_links)
    expect_true(res_v_gmm$lifetime$link_results_released_each_step)
    expect_true(res_v_gmm$lifetime$scratch_released_on_problem_or_chunk_end)
    expect_true(res_v_gmm$lifetime$expected_a_id_chunk_bounded)
    expect_lte(res_v_gmm$lifetime$max_expected_a_id_length, 5L)

    # Deduplicated existing 4bd retained entries
    expected_4bd <- as.double(b14$prob12_view$preflight$total_retained_entries) +
      as.double(b14$prob3_view$preflight$total_retained_entries) +
      as.double(b14$links_view$l1$preflight$total_retained_entries) +
      as.double(b14$links_view$l2$preflight$total_retained_entries) +
      as.double(b14$links_view$l3$preflight$total_retained_entries)
    expect_identical(res_v_gmm$preflight$existing_4bd_retained_entries, expected_4bd)

    # Cluster-level inference covariance transfer (retaining zero-IF cluster "zero")
    res_v_cl <- ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view,
      candidate_rows = b14$candidate_rows,
      frame_owner = b14$frame_owner,
      C = b14$C,
      op = b14$op_cl,
      V_candidate = b14$V_cand_cl,
      control = b14$cov_ctrl,
      reader_control = b14$o_ctrl
    )
    res_d_cl <- ns$certify_ipw_covariance_transfer(
      links = b14$links_dense,
      candidate_rows = b14$candidate_rows,
      frame = b14$frame_dense,
      C = b14$C,
      op = b14$op_cl,
      V_candidate = b14$V_cand_cl,
      control = b14$cov_ctrl
    )
    expect_true(res_v_cl$certified)
    expect_identical(res_v_cl$enclosure$lo, res_d_cl$enclosure$lo)
    expect_identical(res_v_cl$enclosure$hi, res_d_cl$enclosure$hi)
    expect_identical(res_v_cl$epsilon_V, res_d_cl$epsilon_V)
    expect_identical(res_v_cl$cluster_summary$n_declared_clusters, 5L)
    expect_identical(res_v_cl$cluster_summary$n_zero_if_clusters, 1L)
    for (j in 1:2) {
      for (k in 1:2) {
        expect_lte(res_v_cl$enclosure$lo[j, k], b14$ref_exact_V_cl[j, k])
        expect_gte(res_v_cl$enclosure$hi[j, k], b14$ref_exact_V_cl[j, k])
      }
    }

    # Chunk-size and operator-block-size invariance
    for (chk in c(1L, 3L, 14L)) {
      ctrl_chk <- ns$ipw_covariance_transfer_control(
        max_chunk_rows = chk,
        max_operator_block_rows = 2L,
        max_cluster_rows = 14L,
        max_matrix_entries = 15000L,
        max_work_ops = 1000000L,
        max_output_entries = 10000L,
        link_control = b14$link_ctrl
      )
      r_chk <- ns$certify_ipw_covariance_transfer_view(
        links = b14$links_view,
        candidate_rows = b14$candidate_rows,
        frame_owner = b14$frame_owner,
        C = b14$C,
        op = b14$op_cl,
        V_candidate = b14$V_cand_cl,
        control = ctrl_chk,
        reader_control = b14$o_ctrl
      )
      expect_true(r_chk$certified, info = paste(id_mode, chk))
      expect_lte(r_chk$lifetime$max_expected_a_id_length, chk)
      for (j in 1:2) {
        for (k in 1:2) {
          expect_lte(r_chk$enclosure$lo[j, k], b14$ref_exact_V_cl[j, k])
          expect_gte(r_chk$enclosure$hi[j, k], b14$ref_exact_V_cl[j, k])
        }
      }
    }
  }
})

test_that("compact view preserves complete nonobserver rate IF terms (singular observer-only Q vs PD full V) and 10-unit exact/rootball fixtures", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))

  sq <- fx$build_singular_observer_q_view_and_dense(ns, r_ball = 1e-6)
  expect_lt(abs(det(sq$Q_obs)), 1e-17)
  expect_gt(det(sq$V_full), 0.045)

  res_sq_v <- ns$certify_ipw_covariance_transfer_view(
    links = sq$links_view,
    candidate_rows = sq$candidate_rows_full,
    frame_owner = sq$frame_owner,
    C = sq$C,
    op = sq$op,
    V_candidate = sq$V_full,
    control = sq$cov_ctrl,
    reader_control = sq$o_ctrl
  )
  res_sq_d <- ns$certify_ipw_covariance_transfer(
    links = sq$links_dense,
    candidate_rows = sq$candidate_rows_full,
    frame = sq$frame_dense,
    C = sq$C,
    op = sq$op,
    V_candidate = sq$V_full,
    control = sq$cov_ctrl
  )
  expect_true(res_sq_v$certified)
  expect_true(res_sq_v$pd_transferred)
  expect_identical(res_sq_v$enclosure$lo, res_sq_d$enclosure$lo)
  expect_identical(res_sq_v$enclosure$hi, res_sq_d$enclosure$hi)
  expect_identical(res_sq_v$eigen_bounds, res_sq_d$eigen_bounds)
  # Distinct nuisance problems (la and lb have different control cohorts K) -> 0 scratch allocations/reuses
  expect_identical(res_sq_v$n_unique_nuisance_problems, 2L)
  expect_identical(res_sq_v$work$scratch_allocations, 0L)
  expect_identical(res_sq_v$work$scratch_reuse_count, 0L)

  # Singular observer-only Q_obs supplied as V_candidate is INCONCLUSIVE, never certified PD
  res_q_v <- ns$certify_ipw_covariance_transfer_view(
    links = sq$links_view,
    candidate_rows = sq$candidate_rows_full,
    frame_owner = sq$frame_owner,
    C = sq$C,
    op = sq$op,
    V_candidate = sq$Q_obs,
    control = sq$cov_ctrl,
    reader_control = sq$o_ctrl
  )
  expect_false(res_q_v$certified)
  expect_false(res_q_v$pd_transferred)
  expect_identical(res_q_v$status, "INCONCLUSIVE")
  expect_false(res_q_v$singularity_proved)
  expect_false(res_q_v$nonidentification_proved)

  # 10-unit exact rational and nonzero-beta root-ball fixtures
  f10 <- fx$build_10unit_enclosure_fixtures_view_and_dense(ns)
  res_10_v <- ns$certify_ipw_covariance_transfer_view(
    links = f10$fixture1$links_view,
    candidate_rows = f10$fixture1$candidate_rows,
    frame_owner = f10$fixture1$frame_owner,
    C = diag(2L),
    op = f10$op_10,
    V_candidate = f10$fixture1$V_cand,
    control = f10$cov_ctrl,
    reader_control = f10$o_ctrl
  )
  res_10_d <- ns$certify_ipw_covariance_transfer(
    links = f10$fixture1$links_dense,
    candidate_rows = f10$fixture1$candidate_rows,
    frame = f10$fixture1$frame_dense,
    C = diag(2L),
    op = f10$op_10,
    V_candidate = f10$fixture1$V_cand,
    control = f10$cov_ctrl
  )
  expect_true(res_10_v$certified)
  expect_identical(res_10_v$enclosure$lo, res_10_d$enclosure$lo)
  expect_identical(res_10_v$enclosure$hi, res_10_d$enclosure$hi)
  expect_identical(res_10_v$epsilon_V, res_10_d$epsilon_V)

  res_f2_v <- ns$certify_ipw_covariance_transfer_view(
    links = f10$fixture2$links_view,
    candidate_rows = f10$fixture2$candidate_rows,
    frame_owner = f10$fixture2$frame_owner,
    C = matrix(1.0, 1L, 1L),
    op = f10$op_10,
    V_candidate = f10$fixture2$V_cand,
    control = f10$cov_ctrl,
    reader_control = f10$o_ctrl
  )
  res_f2_d <- ns$certify_ipw_covariance_transfer(
    links = f10$fixture2$links_dense,
    candidate_rows = f10$fixture2$candidate_rows,
    frame = f10$fixture2$frame_dense,
    C = matrix(1.0, 1L, 1L),
    op = f10$op_10,
    V_candidate = f10$fixture2$V_cand,
    control = f10$cov_ctrl
  )
  expect_true(res_f2_v$certified)
  expect_identical(res_f2_v$enclosure$lo, res_f2_d$enclosure$lo)
  expect_identical(res_f2_v$enclosure$hi, res_f2_d$enclosure$hi)
  expect_lte(abs(f10$fixture2$V_star - f10$fixture2$V_cand[1L, 1L]), res_f2_v$epsilon_V)
})

test_that("frozen coherent24 6-link panel matches dense covariance transfer bitwise on 6x6 Omega and 3x3 signed C_target covariance", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))
  c24 <- fx$build_coherent24_sixlink_view_and_dense(ns, chunk_rows = 5L, operator_block_rows = 4L)

  op_24 <- ns$ipw_covariance_operator_spec(
    purpose = "GMM",
    row_level = "unit",
    omega = rep(1.0, c24$n),
    centered = TRUE,
    kappa = 1.0 / c24$n
  )

  res_omega_v <- ns$certify_ipw_covariance_transfer_view(
    links = c24$links_view,
    candidate_rows = c24$callback,
    frame_owner = c24$frame_owner,
    C = c24$C_id,
    op = op_24,
    V_candidate = c24$Omega_exact,
    control = c24$cov_ctrl,
    reader_control = c24$o_ctrl
  )
  res_omega_d <- ns$certify_ipw_covariance_transfer(
    links = c24$links_dense,
    candidate_rows = c24$callback,
    frame = c24$frame_dense,
    C = c24$C_id,
    op = op_24,
    V_candidate = c24$Omega_exact,
    control = c24$cov_ctrl
  )
  expect_true(res_omega_v$certified)
  expect_true(res_omega_v$pd_transferred)
  expect_identical(res_omega_v$enclosure$lo, res_omega_d$enclosure$lo)
  expect_identical(res_omega_v$enclosure$hi, res_omega_d$enclosure$hi)
  expect_identical(res_omega_v$epsilon_V, res_omega_d$epsilon_V)
  expect_identical(res_omega_v$eigen_bounds, res_omega_d$eigen_bounds)
  for (j in 1:6) {
    for (k in 1:6) {
      expect_lte(res_omega_v$enclosure$lo[j, k], c24$Omega_exact[j, k])
      expect_gte(res_omega_v$enclosure$hi[j, k], c24$Omega_exact[j, k])
    }
  }

  # 3x3 finite-sample target covariance Sigma_target = (1 / n) * C_target %*% Omega %*% t(C_target) under signed C_target
  op_24_inf <- ns$ipw_covariance_operator_spec(
    purpose = "inference",
    row_level = "unit",
    omega = rep(1.0, c24$n),
    centered = TRUE,
    kappa = 1.0 / (c24$n * c24$n)
  )
  res_tgt_v <- ns$certify_ipw_covariance_transfer_view(
    links = c24$links_view,
    candidate_rows = c24$callback,
    frame_owner = c24$frame_owner,
    C = c24$C_target,
    op = op_24_inf,
    V_candidate = c24$target_cov_exact,
    control = c24$cov_ctrl,
    reader_control = c24$o_ctrl
  )
  res_tgt_d <- ns$certify_ipw_covariance_transfer(
    links = c24$links_dense,
    candidate_rows = c24$callback,
    frame = c24$frame_dense,
    C = c24$C_target,
    op = op_24_inf,
    V_candidate = c24$target_cov_exact,
    control = c24$cov_ctrl
  )
  expect_true(res_tgt_v$certified)
  expect_true(res_tgt_v$pd_transferred)
  expect_identical(res_tgt_v$enclosure$lo, res_tgt_d$enclosure$lo)
  expect_identical(res_tgt_v$enclosure$hi, res_tgt_d$enclosure$hi)
  expect_identical(res_tgt_v$epsilon_V, res_tgt_d$epsilon_V)
  expect_identical(res_tgt_v$eigen_bounds, res_tgt_d$eigen_bounds)
  for (j in 1:3) {
    for (k in 1:3) {
      expect_lte(res_tgt_v$enclosure$lo[j, k], c24$target_cov_exact[j, k])
      expect_gte(res_tgt_v$enclosure$hi[j, k], c24$target_cov_exact[j, k])
    }
  }

  # Signed C_target vs |C_target| witness: replacing C_target with abs(C_target) changes point covariance V
  V_abs_C <- (abs(c24$C_target) %*% c24$Omega_exact %*% t(abs(c24$C_target))) / c24$n
  V_abs_C <- 0.5 * (V_abs_C + t(V_abs_C))
  expect_gt(max(abs(V_abs_C - c24$target_cov_exact)), 0.01)
})

test_that("weighted96 shared-mask 3-link triangle is an honest INCONCLUSIVE singularity diagnostic while 2-link subset is CERTIFIED PD", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))
  w96 <- fx$build_weighted96_shared_mask_view(ns, chunk_rows = 24L, operator_block_rows = 16L)

  # On the shared mask, phi_12 + phi_23 - phi_13 == 0 identically across all 96 units
  expect_lt(max(abs(w96$phi_mat[, 1L] + w96$phi_mat[, 2L] - w96$phi_mat[, 3L])), 1e-13)

  op_96 <- ns$ipw_covariance_operator_spec(
    purpose = "GMM",
    row_level = "unit",
    omega = rep(1.0, w96$n),
    centered = TRUE,
    kappa = 1.0 / w96$n
  )

  # 1. 3-link triangle: enclosure_certified == TRUE, certified == FALSE, pd_transferred == FALSE, status == "INCONCLUSIVE"
  res_3 <- ns$certify_ipw_covariance_transfer_view(
    links = w96$links_3,
    candidate_rows = function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(
        phi = w96$phi_mat[idx, , drop = FALSE],
        id = w96$id[idx],
        link_ids = c("l12", "l23", "l13"),
        a_index = idx
      )
    },
    frame_owner = w96$frame_owner,
    C = diag(1.0, 3L),
    op = op_96,
    V_candidate = w96$V_3,
    control = w96$cov_ctrl,
    reader_control = w96$o_ctrl
  )
  expect_true(res_3$enclosure_certified)
  expect_false(res_3$certified)
  expect_false(res_3$pd_transferred)
  expect_identical(res_3$status, "INCONCLUSIVE")
  expect_true(res_3$inconclusive)
  expect_true("candidate_lambda_min_cert_nonpositive" %in% res_3$inconclusive_reasons)
  expect_false(res_3$singularity_proved)
  expect_false(res_3$nonidentification_proved)
  expect_false(res_3$rank_deficiency_proved)
  expect_identical(res_3$n_unique_nuisance_problems, 1L)
  expect_identical(res_3$work$scratch_allocations, 4L)
  expect_identical(res_3$work$scratch_reuse_count, 8L)

  # 2. 2-link subset (l12, l23): strictly PD and CERTIFIED
  res_2 <- ns$certify_ipw_covariance_transfer_view(
    links = w96$links_2,
    candidate_rows = function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(
        phi = w96$phi_mat[idx, 1:2, drop = FALSE],
        id = w96$id[idx],
        link_ids = c("l12", "l23"),
        a_index = idx
      )
    },
    frame_owner = w96$frame_owner,
    C = diag(1.0, 2L),
    op = op_96,
    V_candidate = w96$V_2,
    control = w96$cov_ctrl,
    reader_control = w96$o_ctrl
  )
  expect_true(res_2$certified)
  expect_true(res_2$enclosure_certified)
  expect_true(res_2$pd_transferred)
  expect_identical(res_2$status, "CERTIFIED")
  expect_identical(res_2$work$scratch_allocations, 4L)
  expect_identical(res_2$work$scratch_reuse_count, 4L)
  expect_gt(res_2$eigen_bounds$lambda_min_star_lower, 0.25)
})

test_that("owner reference identity, forged dense field, preflight cap, callback contract, and purpose-separation refusals occur before outcome reads", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))
  b14 <- fx$build_14unit_3link_view_and_dense(ns, r_ball = 1e-6, chunk_rows = 5L, operator_block_rows = 4L)

  make_traced_reader_ctrl <- function(base_ctrl = b14$o_ctrl) {
    te <- new.env(parent = emptyenv())
    te$calls <- 0L
    rc <- base_ctrl
    attr(rc, "trace_env") <- te
    list(rc = rc, te = te)
  }

  # 1. Distinct frame_owner with identical unit count, IDs, weights, cohorts, and Z_a is refused with 0 outcome reads
  panel_clone <- ns$construct_ipw_panel_owner(b14$pan_df, "id", "time", "y", "g", b14$o_ctrl)
  frame_clone <- ns$construct_ipw_frame_owner(panel_clone, id = b14$id, Z = b14$Z, weights = b14$w, control = b14$o_ctrl)
  tr1 <- make_traced_reader_ctrl()
  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view,
      candidate_rows = b14$candidate_rows,
      frame_owner = frame_clone,
      C = b14$C,
      op = b14$op_gmm,
      V_candidate = b14$V_cand_gmm,
      control = b14$cov_ctrl,
      reader_control = tr1$rc
    ),
    class = "cdid_ipw_owner_mismatch_error"
  )
  expect_identical(tr1$te$calls, 0L)

  # 2. Dense frame spec or dense link enclosure passed to compact-view adapter is refused with 0 outcome reads
  tr2 <- make_traced_reader_ctrl()
  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view,
      candidate_rows = b14$candidate_rows,
      frame_owner = b14$frame_dense,
      C = b14$C,
      op = b14$op_gmm,
      V_candidate = b14$V_cand_gmm,
      control = b14$cov_ctrl,
      reader_control = tr2$rc
    ),
    class = "cdid_ipw_covariance_transfer_frame_error"
  )
  expect_identical(tr2$te$calls, 0L)

  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = b14$links_dense,
      candidate_rows = b14$candidate_rows,
      frame_owner = b14$frame_owner,
      C = b14$C,
      op = b14$op_gmm,
      V_candidate = b14$V_cand_gmm,
      control = b14$cov_ctrl,
      reader_control = tr2$rc
    ),
    class = "cdid_ipw_enclosure_link_error"
  )
  expect_identical(tr2$te$calls, 0L)

  # 3. Forged dense fields inside a sealed compact link enclosure or problem are refused with 0 outcome reads
  forged_ref <- new.env(parent = emptyenv())
  lk_forged <- get("link_enclosure", envir = b14$links_view$l1$ref, inherits = FALSE)
  lk_forged$endpoints$y_left_obs <- rep(0.0, b14$n)
  lk_forged$ref <- forged_ref
  assign("link_enclosure", lk_forged, envir = forged_ref)
  lockBinding("link_enclosure", forged_ref)
  lockEnvironment(forged_ref)
  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = list(l1 = lk_forged, l2 = b14$links_view$l2, l3 = b14$links_view$l3),
      candidate_rows = b14$candidate_rows,
      frame_owner = b14$frame_owner,
      C = b14$C,
      op = b14$op_gmm,
      V_candidate = b14$V_cand_gmm,
      control = b14$cov_ctrl,
      reader_control = tr2$rc
    ),
    class = "cdid_ipw_enclosure_integrity_error"
  )
  expect_identical(tr2$te$calls, 0L)

  # 4. Tight preflight bounds (matrix, work, output, cluster, reader) refuse BEFORE any callback or outcome read
  cb_calls <- 0L
  counting_cb <- function(a_start, a_end) {
    cb_calls <<- cb_calls + 1L
    b14$callback(a_start, a_end)
  }

  tight_me_ctrl <- ns$ipw_covariance_transfer_control(
    max_chunk_rows = 5L, max_operator_block_rows = 4L, max_cluster_rows = 14L,
    max_matrix_entries = 600L, max_work_ops = 1000000L, max_output_entries = 10000L,
    link_control = b14$link_ctrl
  )
  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view, candidate_rows = counting_cb, frame_owner = b14$frame_owner,
      C = b14$C, op = b14$op_gmm, V_candidate = b14$V_cand_gmm,
      control = tight_me_ctrl, reader_control = tr2$rc
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )
  expect_identical(cb_calls, 0L)
  expect_identical(tr2$te$calls, 0L)

  tight_wo_ctrl <- ns$ipw_covariance_transfer_control(
    max_chunk_rows = 5L, max_operator_block_rows = 4L, max_cluster_rows = 14L,
    max_matrix_entries = 15000L, max_work_ops = 40000L, max_output_entries = 10000L,
    link_control = b14$link_ctrl
  )
  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view, candidate_rows = counting_cb, frame_owner = b14$frame_owner,
      C = b14$C, op = b14$op_gmm, V_candidate = b14$V_cand_gmm,
      control = tight_wo_ctrl, reader_control = tr2$rc
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )
  expect_identical(cb_calls, 0L)
  expect_identical(tr2$te$calls, 0L)

  tight_cl_ctrl <- ns$ipw_covariance_transfer_control(
    max_chunk_rows = 5L, max_operator_block_rows = 4L, max_cluster_rows = 3L,
    max_matrix_entries = 15000L, max_work_ops = 1000000L, max_output_entries = 10000L,
    link_control = b14$link_ctrl
  )
  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view, candidate_rows = counting_cb, frame_owner = b14$frame_owner,
      C = b14$C, op = b14$op_cl, V_candidate = b14$V_cand_cl,
      control = tight_cl_ctrl, reader_control = tr2$rc
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )
  expect_identical(cb_calls, 0L)
  expect_identical(tr2$te$calls, 0L)

  # 5. Callback refusals (warning, oversized phi, bad ID, bad link_ids, bad a_index, nonfinite phi) refuse BEFORE outcome reads
  bad_callbacks <- list(
    warn = list(
      cb = function(a_start, a_end) { warning("cb warn"); b14$callback(a_start, a_end) },
      cls = "cdid_ipw_covariance_transfer_callback_error"
    ),
    oversize = list(
      cb = function(a_start, a_end) {
        out <- b14$callback(a_start, a_end)
        out$phi <- matrix(0.0, nrow = 14L, ncol = 3L)
        out
      },
      cls = "cdid_ipw_covariance_transfer_bound_error"
    ),
    bad_id = list(
      cb = function(a_start, a_end) {
        out <- b14$callback(a_start, a_end)
        out$id <- rev(out$id)
        out
      },
      cls = "cdid_ipw_covariance_transfer_id_error"
    ),
    bad_lids = list(
      cb = function(a_start, a_end) {
        out <- b14$callback(a_start, a_end)
        out$link_ids <- c("l2", "l1", "l3")
        out
      },
      cls = "cdid_ipw_covariance_transfer_callback_error"
    ),
    bad_a_idx = list(
      cb = function(a_start, a_end) {
        out <- b14$callback(a_start, a_end)
        out$a_index <- out$a_index + 1L
        out
      },
      cls = "cdid_ipw_covariance_transfer_callback_error"
    ),
    nonfinite_phi = list(
      cb = function(a_start, a_end) {
        out <- b14$callback(a_start, a_end)
        out$phi[1L, 1L] <- NA_real_
        out
      },
      cls = "cdid_ipw_covariance_transfer_candidate_error"
    ),
    forged_flag = list(
      cb = function(a_start, a_end) {
        out <- b14$callback(a_start, a_end)
        out$certified <- TRUE
        out
      },
      cls = "cdid_ipw_covariance_transfer_metadata_error"
    )
  )

  for (nm in names(bad_callbacks)) {
    tr_b <- make_traced_reader_ctrl()
    expect_error(
      ns$certify_ipw_covariance_transfer_view(
        links = b14$links_view,
        candidate_rows = bad_callbacks[[nm]]$cb,
        frame_owner = b14$frame_owner,
        C = b14$C,
        op = b14$op_gmm,
        V_candidate = b14$V_cand_gmm,
        control = b14$cov_ctrl,
        reader_control = tr_b$rc
      ),
      class = bad_callbacks[[nm]]$cls,
      info = nm
    )
    expect_identical(tr_b$te$calls, 0L, info = nm)
  }

  # 6. Purpose separation & D11: clustered operator or non-unit omega for purpose = "GMM" is refused
  expect_error(
    ns$ipw_covariance_operator_spec(
      purpose = "GMM",
      row_level = "cluster",
      omega = b14$omega_cl,
      centered = TRUE,
      kappa = b14$kappa,
      cluster_id = b14$labels_cl,
      cluster_order = b14$order_cl
    ),
    class = "cdid_ipw_covariance_transfer_d11_error"
  )
  op_bad_d11 <- ns$ipw_covariance_operator_spec(
    purpose = "GMM",
    row_level = "unit",
    omega = rep(0.5, b14$n),
    centered = TRUE,
    kappa = b14$kappa
  )
  expect_error(
    ns$certify_ipw_covariance_transfer_view(
      links = b14$links_view,
      candidate_rows = b14$candidate_rows,
      frame_owner = b14$frame_owner,
      C = b14$C,
      op = op_bad_d11,
      V_candidate = b14$V_cand_gmm,
      control = b14$cov_ctrl,
      reader_control = tr2$rc
    ),
    class = "cdid_ipw_covariance_transfer_d11_error"
  )
  expect_identical(tr2$te$calls, 0L)
})


test_that("optional indices are exact before outcome reads and integer phi uses the same numerical contract", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))
  b <- fx$build_14unit_3link_view_and_dense(ns, id_mode = "numeric")
  te <- new.env(parent = emptyenv()); te$calls <- 0L
  rc <- b$o_ctrl; attr(rc, "trace_env") <- te
  for (change in fx$invalid_index_cases) {
    cb_calls <- 0L
    cb <- function(s, t) {
      cb_calls <<- cb_calls + 1L
      x <- b$callback(s, t); x$a_index <- change(seq.int(s, t)); x
    }
    expect_error(ns$certify_ipw_covariance_transfer_view(b$links_view, cb, b$frame_owner,
      b$C, b$op_gmm, b$V_cand_gmm, b$cov_ctrl, rc),
      class = "cdid_ipw_covariance_transfer_callback_error")
    expect_identical(cb_calls, 1L)
    expect_identical(te$calls, 0L)
  }
  cb <- function(s, t) {
    x <- b$callback(s, t)
    x$a_index <- as.double(seq.int(s, t))
    x$phi <- matrix(0L, t - s + 1L, 3L); x
  }
  rv <- ns$certify_ipw_covariance_transfer_view(b$links_view, cb, b$frame_owner,
    b$C, b$op_gmm, b$V_cand_gmm, b$cov_ctrl, b$o_ctrl)
  rd <- ns$certify_ipw_covariance_transfer(b$links_dense, cb, b$frame_dense,
    b$C, b$op_gmm, b$V_cand_gmm, b$cov_ctrl)
  expect_identical(rv$enclosure, rd$enclosure)
  expect_identical(rv$epsilon_V, rd$epsilon_V)
})

test_that("resource preflight reserves distinct references and header bytes with no callback or outcome read", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))
  b <- fx$build_14unit_3link_view_and_dense(ns)
  links <- fx$resource_links(ns, b, paste0("distinct", seq_len(12L)), distinct = TRUE)
  fh <- ns$ipw_cov_transfer_validate_frame_owner_header(b$frame_owner)
  lh <- ns$ipw_cov_transfer_validate_links_view_headers(links, fh)
  expect_identical(lh$n_unique_problems, 12L)
  expect_identical(lh$n_unique_nuisance_views, 12L)
  expect_equal(lh$reference_comparisons, 12 * 11)
  expect_true(all(lh$first_problem))
  te <- new.env(parent = emptyenv()); te$calls <- 0L
  rc <- b$o_ctrl; attr(rc, "trace_env") <- te
  calls <- 0L
  cb <- function(s, t) { calls <<- calls + 1L; stop("preflight must refuse") }
  long <- fx$resource_links(ns, b, paste0(strrep("x", 10000L), 1:3))
  expect_error(ns$certify_ipw_covariance_transfer_view(long, cb, b$frame_owner,
    b$C, b$op_gmm, b$V_cand_gmm, b$cov_ctrl, rc),
    class = "cdid_ipw_covariance_transfer_bound_error")
  expect_identical(calls, 0L); expect_identical(te$calls, 0L)
  # This also checks the entrypoint's stage cap, not only an isolated helper.
  tight <- b$cov_ctrl; tight$max_matrix_entries <- 1L
  expect_error(ns$certify_ipw_covariance_transfer_view(b$links_view, cb, b$frame_owner,
    b$C, b$op_gmm, b$V_cand_gmm, tight, rc),
    class = "cdid_ipw_covariance_transfer_bound_error")
  expect_identical(calls, 0L); expect_identical(te$calls, 0L)
})

test_that("canonical candidate spec and bare callback keep the same header reservation and tight refusal", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_view_fixtures.R"))
  b <- fx$build_14unit_3link_view_and_dense(ns)
  a <- ns$ipw_cov_transfer_view_header_budget(b$links_view, b$frame_owner,
    b$candidate_rows, b$C, b$op_gmm, b$V_cand_gmm, b$cov_ctrl)
  f <- ns$ipw_cov_transfer_view_header_budget(b$links_view, b$frame_owner,
    b$callback, b$C, b$op_gmm, b$V_cand_gmm, b$cov_ctrl)
  expect_identical(a, f)
  z <- ns$certify_ipw_covariance_transfer_view(b$links_view, b$candidate_rows,
    b$frame_owner, b$C, b$op_gmm, b$V_cand_gmm, b$cov_ctrl, b$o_ctrl)
  tight <- b$cov_ctrl; tight$max_work_ops <- floor(z$preflight$combined_total_work_ops_bound - 500)
  calls <- 0L; te <- new.env(parent = emptyenv()); te$calls <- 0L
  rc <- b$o_ctrl; attr(rc, "trace_env") <- te
  cb <- function(s, t) { calls <<- calls + 1L; b$callback(s, t) }
  expect_error(ns$certify_ipw_covariance_transfer_view(b$links_view, cb,
    b$frame_owner, b$C, b$op_gmm, b$V_cand_gmm, tight, rc),
    class = "cdid_ipw_covariance_transfer_bound_error")
  expect_identical(calls, 0L); expect_identical(te$calls, 0L)
})
