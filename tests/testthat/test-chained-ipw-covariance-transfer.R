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

test_that("private IPW covariance transfer entrypoints exist in locked cdid namespace and are unexported", {
  ns <- require_cdid_namespace()
  expect_true(isNamespace(ns))
  expect_true(environmentIsLocked(ns))

  priv_fns <- c(
    "ipw_covariance_transfer_control",
    "validate_ipw_covariance_transfer_control",
    "ipw_covariance_frame_spec",
    "ipw_covariance_operator_spec",
    "ipw_candidate_rows_spec",
    "preflight_ipw_covariance_transfer",
    "certify_ipw_covariance_transfer",
    "validate_ipw_covariance_transfer",
    "ipw_cov_transfer_stream_kernel"
  )
  exports <- getNamespaceExports("cdid")
  for (fn in priv_fns) {
    expect_true(exists(fn, envir = ns, mode = "function", inherits = FALSE), info = fn)
    expect_true(bindingIsLocked(fn, ns), info = fn)
    expect_identical(environment(get(fn, envir = ns, inherits = FALSE)), ns, info = fn)
    expect_false(fn %in% exports, info = fn)
  }
})

test_that("streaming Gram and covariance-transfer kernel encloses all 4 frozen reference fixtures across chunk and operator block sizes", {
  ns <- require_cdid_namespace()
  fx_path <- testthat::test_path("fixtures", "ipw_covariance_transfer_fixtures.R")
  expect_true(file.exists(fx_path))
  fx <- dget(fx_path)

  for (f in fx$operator_fixtures) {
    n <- nrow(f$phi)
    L <- ncol(f$phi)
    d <- nrow(f$C)
    is_cl <- !(f$name %in% c("complete_unit_GMM", "singular_observer_Q_full_PD"))
    unit_to_cl <- if (is_cl) match(f$labels, f$cluster_order) else NULL
    op_norm <- list(
      purpose = f$purpose,
      row_level = if (is_cl) "cluster" else "unit",
      omega = as.double(f$omega),
      centered = isTRUE(f$centered),
      kappa = as.double(f$kappa),
      n_clusters = length(f$cluster_order),
      n_positive_weight_rows = sum(f$omega > 0),
      n_zero_unit_declared_clusters = 0L,
      unit_to_cluster = unit_to_cl,
      cluster_id_key = f$labels,
      cluster_order_key = f$cluster_order
    )

    fetch_fn <- function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(
        certified = TRUE,
        Phi_block = f$phi[idx, , drop = FALSE],
        E_eval_block = f$E_eval[idx, , drop = FALSE],
        E_root_block = f$E_root[idx, , drop = FALSE],
        chunk_4bd_ops = 0,
        max_direct_star_dist = max(f$E_eval[idx, , drop = FALSE] + f$E_root[idx, , drop = FALSE])
      )
    }

    for (b_in in unique(c(1L, 2L, n))) {
      for (b_op in unique(c(1L, 2L, length(f$cluster_order)))) {
        res <- ns$ipw_cov_transfer_stream_kernel(
          n = n, L = L, d = d,
          b_in = b_in, b_op = b_op,
          fetch_chunk = fetch_fn,
          C_dbl = f$C,
          op_norm = op_norm,
          V_cand_dbl = f$V_candidate,
          max_matrix_entries = 10000L,
          preflight = NULL
        )
        expect_true(res$enclosure_certified, info = paste(f$name, b_in, b_op))
        for (j in seq_len(d)) {
          for (k in seq_len(d)) {
            expect_lte(res$enclosure$lo[j, k], f$ref_exact_V[j, k])
            expect_gte(res$enclosure$hi[j, k], f$ref_exact_V[j, k])
          }
        }
        expect_gte(res$epsilon_acc + 1e-18, f$ref_epsilon_acc)
        expect_gte(res$nu_eval + 1e-18, f$ref_nu_eval)
        expect_gte(res$nu_root + 1e-18, f$ref_nu_root)
        expect_gte(res$epsilon_V + 1e-18, f$ref_epsilon_V)

        # Verify separated quadratic error budget identity
        eb <- res$error_budget
        expect_gte(eb$quad_total_term, eb$quad_eval_term + eb$quad_cross_term + eb$quad_root_term - 1e-18)
        expect_gte(eb$epsilon_V, eb$epsilon_acc + eb$linear_total_term + eb$quad_total_term - 1e-18)

        if (f$name == "cluster_rounding_cancellation") {
          # In cluster_rounding_cancellation, 2^53 + 1 - 2^53 cancels in float to 0, so V_candidate is far from exact V
          # and epsilon_acc > 1.75 exceeds lambda_min_cert(V_candidate), producing honest INCONCLUSIVE PD transfer
          # while still certifying the covariance interval enclosure!
          expect_false(res$certified)
          expect_true(res$enclosure_certified)
          expect_false(res$pd_transferred)
          expect_identical(res$status, "INCONCLUSIVE")
          expect_false(res$singularity_proved)
          expect_false(res$nonidentification_proved)
        } else {
          expect_true(res$certified, info = paste(f$name, b_in, b_op))
          expect_true(res$pd_transferred, info = paste(f$name, b_in, b_op))
          expect_identical(res$status, "CERTIFIED")
          expect_gt(res$eigen_bounds$lambda_min_star_lower, 0.0)
        }
      }
    }
  }
})

test_that("end-to-end certify_ipw_covariance_transfer certifies 14-unit 3-link sealed 4bd dataset for both unit GMM and cluster inference", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_transfer_fixtures.R"))
  b14 <- fx$build_14unit_3link_4bd(ns, r_ball = 1e-6, max_chunk_rows = 5L, max_operator_block_rows = 4L)
  f_gmm_ref <- fx$operator_fixtures[[1L]]
  f_cl_ref  <- fx$operator_fixtures[[2L]]

  # 1. Unit-level GMM covariance transfer
  res_gmm <- ns$certify_ipw_covariance_transfer(
    links = b14$links,
    candidate_rows = b14$candidate_rows,
    frame = b14$frame,
    C = b14$C,
    op = b14$op_gmm,
    V_candidate = b14$V_cand_gmm,
    control = b14$control
  )
  expect_true(res_gmm$certified)
  expect_true(res_gmm$enclosure_certified)
  expect_true(res_gmm$pd_transferred)
  expect_identical(res_gmm$status, "CERTIFIED")
  expect_false(res_gmm$inconclusive)
  expect_false(res_gmm$separation_certified)
  expect_false(res_gmm$sampling_certified)
  expect_false(res_gmm$identification_certified)
  expect_false(res_gmm$raw_panel_eligibility_certified)
  expect_false(res_gmm$public_gate_adopted)
  expect_false(res_gmm$singularity_proved)
  expect_false(res_gmm$nonidentification_proved)
  expect_false(res_gmm$rank_deficiency_proved)
  ns$validate_ipw_covariance_transfer(res_gmm, require_certified = TRUE)

  expect_identical(res_gmm$n_a_units, 14.0)
  expect_identical(res_gmm$n_links, 3.0)
  expect_identical(res_gmm$n_moments, 2.0)
  expect_identical(res_gmm$n_unique_nuisance_problems, 2L)
  expect_identical(res_gmm$frame$embedding, "identity")
  expect_false(res_gmm$frame$larger_frame_transport_applied)
  expect_false(res_gmm$frame$a_rows_rescaled)
  expect_false(res_gmm$frame$zero_rows_invented)
  expect_false(res_gmm$operator$case_weights_reapplied)

  for (j in 1:2) {
    for (k in 1:2) {
      expect_lte(res_gmm$enclosure$lo[j, k], f_gmm_ref$ref_exact_V[j, k])
      expect_gte(res_gmm$enclosure$hi[j, k], f_gmm_ref$ref_exact_V[j, k])
    }
  }
  expect_gt(res_gmm$nu_eval, 0.0)
  expect_gt(res_gmm$nu_root, 0.0)
  expect_gt(res_gmm$error_budget$quad_cross_term, 0.0)
  expect_gt(res_gmm$eigen_bounds$lambda_min_star_lower, 0.09)
  expect_lte(res_gmm$work$actual_4bd_interval_ops, res_gmm$preflight$total_4bd_stream_work_ops_bound)
  expect_lte(res_gmm$work$actual_component_interval_ops, res_gmm$preflight$component_work_ops_bound)
  expect_lte(res_gmm$work$total_actual_interval_ops, res_gmm$preflight$total_work_ops_bound)

  # Verify existing 4bd retained storage deduplicates shared nuisance problem prob12
  expected_4bd_entries <- as.double(b14$prob12$preflight$total_retained_entries) +
    as.double(b14$prob3$preflight$total_retained_entries) +
    as.double(b14$links$l1$preflight$total_retained_entries) +
    as.double(b14$links$l2$preflight$total_retained_entries) +
    as.double(b14$links$l3$preflight$total_retained_entries)
  expect_identical(res_gmm$preflight$existing_4bd_retained_entries, expected_4bd_entries)

  # 2. Cluster-level inference covariance transfer (retaining zero-IF cluster "zero" with omega = 7/5)
  res_cl <- ns$certify_ipw_covariance_transfer(
    links = b14$links,
    candidate_rows = b14$candidate_rows,
    frame = b14$frame,
    C = b14$C,
    op = b14$op_cl,
    V_candidate = b14$V_cand_cl,
    control = b14$control
  )
  expect_true(res_cl$certified)
  expect_true(res_cl$enclosure_certified)
  expect_true(res_cl$pd_transferred)
  expect_identical(res_cl$status, "CERTIFIED")
  ns$validate_ipw_covariance_transfer(res_cl, require_certified = TRUE)

  for (j in 1:2) {
    for (k in 1:2) {
      expect_lte(res_cl$enclosure$lo[j, k], f_cl_ref$ref_exact_V[j, k])
      expect_gte(res_cl$enclosure$hi[j, k], f_cl_ref$ref_exact_V[j, k])
    }
  }
  expect_identical(res_cl$cluster_summary$n_declared_clusters, 5L)
  expect_identical(res_cl$cluster_summary$n_zero_if_clusters, 1L)
  expect_true(res_cl$work$weight_rounding_occurred)
  expect_gt(res_cl$eigen_bounds$lambda_min_star_lower, 0.0008)

  # 3. Verify chunk-size invariance of enclosure validity across max_chunk_rows in {1, 3, 14}
  for (chk in c(1L, 3L, 14L)) {
    ctrl_chk <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = chk,
      max_operator_block_rows = 2L,
      max_cluster_rows = 14L,
      max_matrix_entries = 10000L,
      max_work_ops = 1000000L,
      max_output_entries = 10000L,
      link_control = b14$link_control
    )
    r_chk <- ns$certify_ipw_covariance_transfer(
      links = b14$links,
      candidate_rows = b14$candidate_rows,
      frame = b14$frame,
      C = b14$C,
      op = b14$op_cl,
      V_candidate = b14$V_cand_cl,
      control = ctrl_chk
    )
    expect_true(r_chk$certified, info = paste("chunk size", chk))
    for (j in 1:2) {
      for (k in 1:2) {
        expect_lte(r_chk$enclosure$lo[j, k], f_cl_ref$ref_exact_V[j, k])
        expect_gte(r_chk$enclosure$hi[j, k], f_cl_ref$ref_exact_V[j, k])
      }
    }
  }
})

test_that("real 4bd sealed 2-link dataset proves singular observer-only Q vs strictly PD full complete-IF covariance V", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_transfer_fixtures.R"))
  sq <- fx$build_singular_observer_q_4bd(ns, r_ball = 1e-6)

  # Observer-only matrix Q_obs is singular (det(Q_obs) == 0), whereas full V_full is strictly PD
  expect_lt(abs(det(sq$Q_obs)), 1e-17)
  expect_gt(det(sq$V_full), 0.045)

  # Full complete-IF covariance transfer succeeds and certifies strict PD
  res_full <- ns$certify_ipw_covariance_transfer(
    links = sq$links,
    candidate_rows = sq$candidate_rows_full,
    frame = sq$frame,
    C = sq$C,
    op = sq$op,
    V_candidate = sq$V_full,
    control = sq$control
  )
  expect_true(res_full$certified)
  expect_true(res_full$pd_transferred)
  expect_identical(res_full$status, "CERTIFIED")
  expect_gt(res_full$eigen_bounds$lambda_min_star_lower, 0.033)

  # If a caller supplies the singular observer-only Q_obs as V_candidate, PD transfer is refused (INCONCLUSIVE)
  res_q_cand <- ns$certify_ipw_covariance_transfer(
    links = sq$links,
    candidate_rows = sq$candidate_rows_full,
    frame = sq$frame,
    C = sq$C,
    op = sq$op,
    V_candidate = sq$Q_obs,
    control = sq$control
  )
  expect_false(res_q_cand$certified)
  expect_false(res_q_cand$pd_transferred)
  expect_identical(res_q_cand$status, "INCONCLUSIVE")
  expect_false(res_q_cand$singularity_proved)
  expect_false(res_q_cand$nonidentification_proved)
})

test_that("10-unit 2-link exact rational and nonzero-beta root-ball fixtures certify covariance transfer", {
  ns <- require_cdid_namespace()
  ref_4bd <- dget(testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R"))
  f1 <- ref_4bd$fixture1_exact_rational
  f2 <- ref_4bd$fixture2_nonzero_beta_rootball

  cc <- ns$logit_certificate_control(18L, 40, 16L, 4L, 2500L, 500000L)
  pc <- ns$ipw_nuisance_problem_control(4L, 5000L, 2500L, 500000L)
  lc <- ns$ipw_link_enclosure_control(10L, 5000L, 2500L, 500000L, 5000L)
  ctrl <- ns$ipw_covariance_transfer_control(
    max_chunk_rows = 4L,
    max_operator_block_rows = 4L,
    max_cluster_rows = 10L,
    max_matrix_entries = 5000L,
    max_work_ops = 500000L,
    max_output_entries = 5000L,
    link_control = lc
  )

  prob1 <- ns$construct_certified_ipw_nuisance_problem(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G,
    pair_observed = f1$pair_observed, cohorts = f1$cohorts,
    beta0 = f1$beta_star_float, F = f1$F_float, r = 1e-5,
    certificate_control = cc, problem_control = pc
  )
  lk1 <- ns$enclose_ipw_link(prob1, f1$id, f1$yL_1_float, f1$yR_1_float, f1$link1$lambda_beta_float, lc, link_id = "link1")
  lk2 <- ns$enclose_ipw_link(prob1, f1$id, f1$yL_2_float, f1$yR_2_float, f1$link2$lambda_beta_float, lc, link_id = "link2")

  phi_10 <- cbind(f1$link1$phi_float, f1$link2$phi_float)
  V_10 <- crossprod(scale(phi_10, center = TRUE, scale = FALSE)) / 10.0
  V_10 <- 0.5 * (V_10 + t(V_10))

  cand_10 <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(id = f1$id[idx], link_ids = c("link1", "link2"), phi = phi_10[idx, , drop = FALSE])
    },
    link_ids = c("link1", "link2")
  )
  frame_10 <- ns$ipw_covariance_frame_spec(id = f1$id, normalizer_n = 10, embedding = "identity", w = f1$w_float, G = f1$G)
  op_10 <- ns$ipw_covariance_operator_spec(
    purpose = "GMM", row_level = "unit", omega = rep(1.0, 10L), centered = TRUE, kappa = 0.1
  )

  res_10 <- ns$certify_ipw_covariance_transfer(
    links = list(link1 = lk1, link2 = lk2),
    candidate_rows = cand_10,
    frame = frame_10,
    C = diag(2L),
    op = op_10,
    V_candidate = V_10,
    control = ctrl
  )
  expect_true(res_10$certified)
  expect_true(res_10$pd_transferred)
  expect_identical(res_10$frame$id_type, "numeric")
  expect_gt(res_10$eigen_bounds$lambda_min_star_lower, 0.35)

  # Fixture 2 (nonzero beta0 with active target cohort rate loading lambda_{r,1} != 0)
  prob2 <- ns$construct_certified_ipw_nuisance_problem(
    id = f2$id, Z = f2$Z_float, w = f2$w_float, G = f2$G,
    pair_observed = f2$pair_observed, cohorts = f2$cohorts,
    beta0 = f2$beta0_perturbed_float, F = f2$F_float, r = f2$r_ball,
    certificate_control = cc, problem_control = pc
  )
  lk2_single <- ns$enclose_ipw_link(
    prob2, f2$id, f2$yL_float, f2$yR_float,
    f2$at_beta0_perturbed$lambda_beta, lc, link_id = "f2_link"
  )
  phi_f2 <- matrix(f2$at_beta0_perturbed$phi, ncol = 1L)
  V_f2 <- matrix(sum((phi_f2 - mean(phi_f2))^2) / 10.0, 1L, 1L)
  cand_f2 <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(id = f2$id[idx], link_ids = "f2_link", phi = phi_f2[idx, , drop = FALSE])
    },
    link_ids = "f2_link"
  )
  frame_f2 <- ns$ipw_covariance_frame_spec(id = f2$id, normalizer_n = 10, embedding = "identity", w = f2$w_float, G = f2$G)
  res_f2 <- ns$certify_ipw_covariance_transfer(
    links = list(f2_link = lk2_single),
    candidate_rows = cand_f2,
    frame = frame_f2,
    C = matrix(1.0, 1L, 1L),
    op = op_10,
    V_candidate = V_f2,
    control = ctrl
  )
  expect_true(res_f2$certified)
  expect_true(res_f2$pd_transferred)
  # True covariance at beta_star must lie inside [V_candidate - epsilon_V, V_candidate + epsilon_V]
  phi_f2_star <- f2$at_beta_star$phi
  V_f2_star <- sum((phi_f2_star - mean(phi_f2_star))^2) / 10.0
  expect_lte(abs(V_f2_star - V_f2[1L, 1L]), res_f2$epsilon_V)
})
