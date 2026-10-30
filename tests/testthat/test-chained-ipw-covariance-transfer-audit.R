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

make_source_mutant <- function(fn, target, replacement, ns) {
  src <- paste(deparse(fn, width.cutoff = 500L), collapse = "\n")
  if (!grepl(target, src, fixed = TRUE)) {
    stop(paste0("Source mutant target not found: ", target), call. = FALSE)
  }
  mut_src <- sub(target, replacement, src, fixed = TRUE)
  if (identical(mut_src, src)) {
    stop(paste0("Source mutant produced identical code for: ", target), call. = FALSE)
  }
  eval(parse(text = mut_src), envir = ns)
}

enclosure_contains_matrix <- function(enc, mat) {
  all(enc$lo <= mat & mat <= enc$hi)
}

test_that("all 15 source mutants of ipw_cov_transfer_stream_kernel and certify_ipw_covariance_transfer are actively rejected", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_transfer_fixtures.R"))
  f_gmm <- fx$operator_fixtures[[1L]]
  f_cl  <- fx$operator_fixtures[[2L]]
  f_rc  <- fx$operator_fixtures[[3L]]
  b14   <- fx$build_14unit_3link_4bd(ns, r_ball = 1e-6, max_chunk_rows = 5L, max_operator_block_rows = 2L)

  orig_kern <- ns$ipw_cov_transfer_stream_kernel
  orig_cert <- ns$certify_ipw_covariance_transfer

  op_norm_cl <- list(
    purpose = "inference",
    row_level = "cluster",
    omega = as.double(f_cl$omega),
    centered = TRUE,
    kappa = as.double(f_cl$kappa),
    n_clusters = length(f_cl$cluster_order),
    n_positive_weight_rows = sum(f_cl$omega > 0),
    n_zero_unit_declared_clusters = 0L,
    unit_to_cluster = match(f_cl$labels, f_cl$cluster_order),
    cluster_id_key = f_cl$labels,
    cluster_order_key = f_cl$cluster_order
  )
  fetch_cl <- function(a_start, a_end) {
    idx <- seq.int(a_start, a_end)
    list(
      certified = TRUE,
      Phi_block = f_cl$phi[idx, , drop = FALSE],
      E_eval_block = f_cl$E_eval[idx, , drop = FALSE],
      E_root_block = f_cl$E_root[idx, , drop = FALSE],
      chunk_4bd_ops = 0,
      max_direct_star_dist = 0
    )
  }

  op_norm_rc <- list(
    purpose = "inference",
    row_level = "cluster",
    omega = as.double(f_rc$omega),
    centered = FALSE,
    kappa = as.double(f_rc$kappa),
    n_clusters = length(f_rc$cluster_order),
    n_positive_weight_rows = sum(f_rc$omega > 0),
    n_zero_unit_declared_clusters = 0L,
    unit_to_cluster = match(f_rc$labels, f_rc$cluster_order),
    cluster_id_key = f_rc$labels,
    cluster_order_key = f_rc$cluster_order
  )
  fetch_rc <- function(a_start, a_end) {
    idx <- seq.int(a_start, a_end)
    list(
      certified = TRUE,
      Phi_block = f_rc$phi[idx, , drop = FALSE],
      E_eval_block = f_rc$E_eval[idx, , drop = FALSE],
      E_root_block = f_rc$E_root[idx, , drop = FALSE],
      chunk_4bd_ops = 0,
      max_direct_star_dist = 0
    )
  }

  # 1. Mutant: map_as_exact_float (omits transform residual res_U on Uhat = Phi %*% t(C))
  mut_map <- make_source_mutant(
    orig_kern,
    "u_iv <- gmm_iv_add(gmm_iv(Uhat[u, j]), c(res_U$lo[u, j], res_U$hi[u, j]))",
    "u_iv <- gmm_iv(Uhat[u, j])",
    ns
  )
  # Construct 1-row uncentered unit operator where Phi %*% t(C) has catastrophic float cancellation:
  # Phi = [2^53, 1, -2^53], C = [1, 1, 1] -> exact U = 1, float Uhat = (2^53 + 1) - 2^53 = 0
  op_1row <- list(
    purpose = "GMM", row_level = "unit", omega = 1.0, centered = FALSE, kappa = 1.0,
    n_clusters = 1L, n_positive_weight_rows = 1L, n_zero_unit_declared_clusters = 0L,
    unit_to_cluster = NULL, cluster_id_key = NULL, cluster_order_key = NULL
  )
  fetch_1row <- function(s, e) {
    list(
      certified = TRUE,
      Phi_block = matrix(c(2^53, 1.0, -2^53), nrow = 1L, ncol = 3L),
      E_eval_block = matrix(0.0, 1L, 3L),
      E_root_block = matrix(0.0, 1L, 3L),
      chunk_4bd_ops = 0, max_direct_star_dist = 0
    )
  }
  C_1row <- matrix(c(1.0, 1.0, 1.0), nrow = 1L, ncol = 3L)
  exact_V_1row <- matrix(1.0, 1L, 1L)
  good_1row <- orig_kern(1L, 3L, 1L, 1L, 1L, fetch_1row, C_1row, op_1row, matrix(1.0, 1L, 1L), 1000L)
  bad_1row  <- mut_map(1L, 3L, 1L, 1L, 1L, fetch_1row, C_1row, op_1row, matrix(1.0, 1L, 1L), 1000L)
  expect_true(enclosure_contains_matrix(good_1row$enclosure, exact_V_1row))
  expect_false(enclosure_contains_matrix(bad_1row$enclosure, exact_V_1row))

  # 2. Mutant: wrong_C_transpose (uses C_dbl instead of t(C_dbl))
  mut_tC <- make_source_mutant(orig_kern, "tC <- t(C_dbl)", "tC <- C_dbl", ns)
  fetch_tC <- function(s, e) {
    list(
      certified = TRUE,
      Phi_block = matrix(c(1.0, 0.0), nrow = 1L, ncol = 2L),
      E_eval_block = matrix(0.0, 1L, 2L),
      E_root_block = matrix(0.0, 1L, 2L),
      chunk_4bd_ops = 0, max_direct_star_dist = 0
    )
  }
  C_nonsym <- matrix(c(1.0, 2.0, 0.0, 1.0), nrow = 2L, ncol = 2L) # C[,1] = (1, 2), C[1,] = (1, 0)
  exact_V_tC <- tcrossprod(c(1.0, 2.0))
  good_tC <- orig_kern(1L, 2L, 2L, 1L, 1L, fetch_tC, C_nonsym, op_1row, exact_V_tC, 1000L)
  bad_tC  <- mut_tC(1L, 2L, 2L, 1L, 1L, fetch_tC, C_nonsym, op_1row, exact_V_tC, 1000L)
  expect_true(enclosure_contains_matrix(good_tC$enclosure, exact_V_tC))
  expect_false(enclosure_contains_matrix(bad_tC$enclosure, exact_V_tC))

  # 3. Mutant: cluster_as_exact_float (ignores outward cluster interval sum, using floating cluster_z_hat)
  mut_cl_float <- make_source_mutant(
    orig_kern,
    "z_sum_iv <- gmm_iv_add(c(cluster_z_lo[cu, j], cluster_z_hi[cu, j]), c(U_lo[u, j], U_hi[u, j]))",
    "z_sum_iv <- c(zh_new, zh_new)",
    ns
  )
  bad_cl_float <- mut_cl_float(4L, 2L, 2L, 2L, 2L, fetch_rc, f_rc$C, op_norm_rc, f_rc$V_candidate, 10000L)
  expect_false(enclosure_contains_matrix(bad_cl_float$enclosure, f_rc$ref_exact_V))

  # 4. Mutant: sqrt_weight_as_exact (sets E_mat[r, j] <- 0, ignoring sqrt(omega) and cluster radius)
  mut_sqrt_w <- make_source_mutant(
    orig_kern,
    "E_mat[r, j] <- ipw_cov_transfer_iv_point_dist(y_rj_iv, yh_rj)",
    "E_mat[r, j] <- 0",
    ns
  )
  bad_sqrt_w <- mut_sqrt_w(4L, 2L, 2L, 2L, 2L, fetch_rc, f_rc$C, op_norm_rc, f_rc$V_candidate, 10000L)
  expect_false(enclosure_contains_matrix(bad_sqrt_w$enclosure, f_rc$ref_exact_V))

  # 5. Mutant: drop_gram_widening (omits D_jk widening of point Gram interval)
  mut_no_widen <- make_source_mutant(
    orig_kern,
    "widened_jk <- gmm_iv_add(c(gram_pt$lo[j, k], gram_pt$hi[j, k]), c(-D_jk, D_jk))",
    "widened_jk <- c(gram_pt$lo[j, k], gram_pt$hi[j, k])",
    ns
  )
  bad_no_widen <- mut_no_widen(4L, 2L, 2L, 2L, 2L, fetch_rc, f_rc$C, op_norm_rc, f_rc$V_candidate, 10000L)
  expect_false(enclosure_contains_matrix(bad_no_widen$enclosure, f_rc$ref_exact_V))

  # 6. Mutant: wrong_gram_transpose (computes Yhat %*% t(Yhat) instead of t(Yhat) %*% Yhat on a 2x2 block)
  mut_gram_tr <- make_source_mutant(
    orig_kern,
    "gram_pt <- gmm_iv_mat_residual(t(Yhat), Yhat, zero_d_d, max_matrix_entries = max_matrix_entries)",
    "gram_pt <- gmm_iv_mat_residual(Yhat, t(Yhat), zero_d_d, max_matrix_entries = max_matrix_entries)",
    ns
  )
  bad_gram_tr <- mut_gram_tr(4L, 2L, 2L, 4L, 2L, fetch_rc, f_rc$C, op_norm_rc, f_rc$V_candidate, 10000L)
  expect_false(enclosure_contains_matrix(bad_gram_tr$enclosure, f_rc$ref_exact_V))

  # 7. Mutant: omit_global_center (skips global centering subtraction t_iv * t_iv / W_iv)
  mut_no_cent <- make_source_mutant(orig_kern, "if (centered) {", "if (FALSE) {", ns)
  bad_no_cent <- mut_no_cent(14L, 3L, 2L, 5L, 2L, fetch_cl, f_cl$C, op_norm_cl, f_cl$V_candidate, 10000L)
  expect_false(enclosure_contains_matrix(bad_no_cent$enclosure, f_cl$ref_exact_V))

  # 8. Mutant: center_per_batch (centers each operator block of size b_op = 2 instead of globally)
  mut_batch_cent <- make_source_mutant(
    orig_kern,
    "n_operator_blocks_seen <<- n_operator_blocks_seen + 1L",
    "for (jj in seq_len(d)) for (kk in seq_len(d)) S_iv[jj, kk, ] <<- gmm_iv_sub(S_iv[jj, kk, ], gmm_iv_div(gmm_iv_mul(t_iv[jj, ], t_iv[kk, ]), W_iv)); t_iv[,] <<- 0; W_iv <<- c(0, 0); n_operator_blocks_seen <<- n_operator_blocks_seen + 1L",
    ns
  )
  mut_batch_cent <- make_source_mutant(mut_batch_cent, "if (centered) {", "if (FALSE) {", ns)
  mut_batch_cent <- make_source_mutant(mut_batch_cent, "if (!is.finite(W_iv[1L]) || W_iv[1L] <= 0) {", "if (FALSE) {", ns)
  bad_batch_cent <- mut_batch_cent(14L, 3L, 2L, 5L, 2L, fetch_cl, f_cl$C, op_norm_cl, f_cl$V_candidate, 10000L)
  expect_false(enclosure_contains_matrix(bad_batch_cent$enclosure, f_cl$ref_exact_V))

  # 9. Mutant: drop_zero_cluster_mass (omits zero-IF cluster "zero" weight omega = 7/5 from W_iv)
  mut_drop_zero_cl <- make_source_mutant(
    orig_kern,
    "W_iv <<- gmm_iv_add(W_iv, w_iv)",
    "if (any(z_hat_blk[r, ] != 0)) W_iv <<- gmm_iv_add(W_iv, w_iv)",
    ns
  )
  bad_drop_zero_cl <- mut_drop_zero_cl(14L, 3L, 2L, 5L, 2L, fetch_cl, f_cl$C, op_norm_cl, f_cl$V_candidate, 10000L)
  expect_false(enclosure_contains_matrix(bad_drop_zero_cl$enclosure, f_cl$ref_exact_V))

  # 10. Mutant: double_caseweights (multiplies complete IF rows by case weights w a second time)
  mut_dbl_w <- make_source_mutant(
    orig_cert,
    "E_root_block[, l] <- sep_rows$root_error",
    "Phi_block[, l] <- Phi_block[, l] * link_hdr$problems[[l]]$raw$w[a_seq]; E_root_block[, l] <- sep_rows$root_error",
    ns
  )
  bad_dbl_w <- mut_dbl_w(
    links = b14$links, candidate_rows = b14$candidate_rows, frame = b14$frame,
    C = b14$C, op = b14$op_cl, V_candidate = b14$V_cand_cl, control = b14$control
  )
  expect_false(enclosure_contains_matrix(bad_dbl_w$enclosure, f_cl$ref_exact_V))

  # 11. Mutant: drop_nonobserver (zeroes out non-observer rows !in_observed)
  mut_drop_nonobs <- make_source_mutant(
    orig_cert,
    "E_root_block[, l] <- sep_rows$root_error",
    "Phi_block[, l] <- ifelse(link_hdr$problems[[l]]$masks$in_observed[a_seq], Phi_block[, l], 0); E_root_block[, l] <- sep_rows$root_error",
    ns
  )
  bad_drop_nonobs <- mut_drop_nonobs(
    links = b14$links, candidate_rows = b14$candidate_rows, frame = b14$frame,
    C = b14$C, op = b14$op_cl, V_candidate = b14$V_cand_cl, control = b14$control
  )
  expect_false(enclosure_contains_matrix(bad_drop_nonobs$enclosure, f_cl$ref_exact_V))

  # 12. Mutant: drop_eval_root_cross (omits 2 * nu_eval * nu_root from quadratic error term)
  mut_no_cross <- make_source_mutant(
    orig_kern,
    "quad_cross_term <- gmm_iv_mul(gmm_iv(2), gmm_iv_mul(gmm_iv(nu_eval), gmm_iv(nu_root)))[2L]",
    "quad_cross_term <- 0",
    ns
  )
  mut_no_cross <- make_source_mutant(
    mut_no_cross,
    "gmm_iv_mul(gmm_iv(nu_total), gmm_iv(nu_total))[2L],",
    "0,",
    ns
  )
  # Scalar 2-unit witness: phi_cand = (0, 0), V_cand = 0, E_eval = (1, 1), E_root = (1, 1) -> phi_star = (-2, 2) -> V_star = 4
  op_2u <- list(
    purpose = "GMM", row_level = "unit", omega = c(1.0, 1.0), centered = TRUE, kappa = 0.5,
    n_clusters = 2L, n_positive_weight_rows = 2L, n_zero_unit_declared_clusters = 0L,
    unit_to_cluster = NULL, cluster_id_key = NULL, cluster_order_key = NULL
  )
  fetch_cross_witness <- function(s, e) {
    idx <- seq.int(s, e)
    list(
      certified = TRUE,
      Phi_block = matrix(0.0, nrow = length(idx), ncol = 1L),
      E_eval_block = matrix(1.0, nrow = length(idx), ncol = 1L),
      E_root_block = matrix(1.0, nrow = length(idx), ncol = 1L),
      chunk_4bd_ops = 0, max_direct_star_dist = 2.0
    )
  }
  good_cross <- orig_kern(2L, 1L, 1L, 2L, 2L, fetch_cross_witness, matrix(1.0, 1L, 1L), op_2u, matrix(0.0, 1L, 1L), 1000L)
  bad_cross  <- mut_no_cross(2L, 1L, 1L, 2L, 2L, fetch_cross_witness, matrix(1.0, 1L, 1L), op_2u, matrix(0.0, 1L, 1L), 1000L)
  expect_gte(good_cross$epsilon_V, 4.0) # Covers true V_star = 4
  expect_lt(bad_cross$epsilon_V, 4.0)   # Underestimates true V_star = 4 (gives 2 < 4)

  # 13. Mutant: drop_quad (omits (nu_eval + nu_root)^2 from perturbation_total)
  mut_no_quad <- make_source_mutant(
    orig_kern,
    "perturbation_total <- gmm_iv_add(gmm_iv(linear_total_term), gmm_iv(quad_total_term))[2L]",
    "perturbation_total <- linear_total_term",
    ns
  )
  bad_quad <- mut_no_quad(2L, 1L, 1L, 2L, 2L, fetch_cross_witness, matrix(1.0, 1L, 1L), op_2u, matrix(0.0, 1L, 1L), 1000L)
  expect_lt(bad_quad$epsilon_V, 4.0) # Gives 0 < 4

  # 14. Mutant: use_candidate_norm_only (uses ||V_candidate||_F without + epsilon_acc in v_upper)
  mut_cand_norm_only <- make_source_mutant(
    orig_kern,
    "v_upper <- gmm_iv_add(gmm_iv(v_cand_norm_upper), gmm_iv(epsilon_acc))[2L]",
    "v_upper <- v_cand_norm_upper",
    ns
  )
  # Witness: Phi_cand = (-2, 2) so V_exact(Phi_cand) = 4, but V_candidate = 0 (so epsilon_acc = 4, v_cand_norm = 0).
  # With E_root = (2, 2), Phi_star = (-4, 4) so V(Phi_star) = 16.
  # True budget: eps_acc(4) + 2*sqrt(4)*2(8) + 2^2(4) = 16. Mutant budget: 4 + 0 + 4 = 8 < 16!
  fetch_vupper_witness <- function(s, e) {
    idx <- seq.int(s, e)
    phi_w <- matrix(c(-2.0, 2.0)[idx], ncol = 1L)
    list(
      certified = TRUE,
      Phi_block = phi_w,
      E_eval_block = matrix(0.0, nrow = length(idx), ncol = 1L),
      E_root_block = matrix(2.0, nrow = length(idx), ncol = 1L),
      chunk_4bd_ops = 0, max_direct_star_dist = 2.0
    )
  }
  good_vup <- orig_kern(2L, 1L, 1L, 2L, 2L, fetch_vupper_witness, matrix(1.0, 1L, 1L), op_2u, matrix(0.0, 1L, 1L), 1000L)
  bad_vup  <- mut_cand_norm_only(2L, 1L, 1L, 2L, 2L, fetch_vupper_witness, matrix(1.0, 1L, 1L), op_2u, matrix(0.0, 1L, 1L), 1000L)
  expect_gte(good_vup$epsilon_V, 16.0)
  expect_lt(bad_vup$epsilon_V, 16.0)

  # 15. Mutant: relabel_direct_star_as_root (uses sep_rows$direct_star_error instead of sep_rows$root_error)
  mut_relabel <- make_source_mutant(
    orig_cert,
    "E_root_block[, l] <- sep_rows$root_error",
    "E_root_block[, l] <- sep_rows$direct_star_error",
    ns
  )
  # When phi_candidate is perturbed by +0.5 on in-risk row 1 (G=2), direct_star_error includes eval_error (~0.5)
  # and differs from the pure root perturbation radius sep_rows$root_error (~1e-5)!
  phi_shifted <- b14$phi_cand
  phi_shifted[1L, ] <- phi_shifted[1L, ] + 0.5
  cand_shifted <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(
        id = b14$id[idx],
        link_ids = c("l1", "l2", "l3"),
        phi = phi_shifted[idx, , drop = FALSE]
      )
    },
    link_ids = c("l1", "l2", "l3")
  )
  good_sep <- orig_cert(b14$links, cand_shifted, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control)
  bad_sep  <- mut_relabel(b14$links, cand_shifted, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control)
  expect_lt(good_sep$nu_root, 1e-3)
  expect_gt(bad_sep$nu_root, 0.05) # Mutant conflated evaluation shift 0.5 into nu_root!
})

test_that("malformed, unsealed, order-mismatched, larger-frame, D11, callback, and preflight controls refuse before execution", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_transfer_fixtures.R"))
  b14 <- fx$build_14unit_3link_4bd(ns, r_ball = 1e-6, max_chunk_rows = 5L, max_operator_block_rows = 4L)

  # 1. Missing required arguments in certify_ipw_covariance_transfer and ipw_covariance_transfer_control
  expect_error(
    ns$certify_ipw_covariance_transfer(b14$links, b14$candidate_rows, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm),
    class = "cdid_ipw_covariance_transfer_input_error"
  )
  expect_error(
    ns$ipw_covariance_transfer_control(max_chunk_rows = 5L),
    class = "cdid_ipw_covariance_transfer_control_error"
  )

  # 2. Tampered / unsealed 4bd link object refused
  bad_l1 <- b14$links$l1
  bad_l1$link_id <- "tampered_l1"
  expect_error(
    ns$certify_ipw_covariance_transfer(
      list(tampered_l1 = bad_l1, l2 = b14$links$l2, l3 = b14$links$l3),
      b14$candidate_rows, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_enclosure_integrity_error"
  )

  # 3. Decorative shortcut attributes on frame, op, C, V_candidate, or control are refused
  frame_dec <- b14$frame
  attr(frame_dec, "same_frame") <- TRUE
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, frame_dec, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_metadata_error"
  )
  V_dec <- b14$V_cand_gmm
  attr(V_dec, "certified_pd") <- TRUE
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, b14$frame, b14$C, b14$op_gmm, V_dec, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_metadata_error"
  )

  # 4. Larger-frame embedding, normalizer mismatch, or shortcut fields in frame are refused
  expect_error(
    ns$ipw_covariance_frame_spec(id = b14$id, normalizer_n = 14, embedding = "pad_zeros"),
    class = "cdid_ipw_covariance_transfer_frame_error"
  )
  frame_bad_norm <- list(id = c(b14$id, "unit_14"), normalizer_n = 15, embedding = "identity")
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, frame_bad_norm, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_frame_error"
  )
  frame_shortcut <- list(id = b14$id, normalizer_n = 14, embedding = "identity", allow_rescale = TRUE)
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, frame_shortcut, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_metadata_error"
  )

  # 5. Reordered or mismatched unit IDs, case weights w, or cohorts G between frame and source links are refused
  frame_reordered <- ns$ipw_covariance_frame_spec(
    id = rev(b14$id), normalizer_n = 14, embedding = "identity", w = b14$w, G = b14$G
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, frame_reordered, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_id_error"
  )
  frame_bad_w <- ns$ipw_covariance_frame_spec(
    id = b14$id, normalizer_n = 14, embedding = "identity", w = rev(b14$w), G = b14$G
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, frame_bad_w, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_frame_error"
  )
  frame_bad_G <- ns$ipw_covariance_frame_spec(
    id = b14$id, normalizer_n = 14, embedding = "identity", w = b14$w, G = rev(b14$G)
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, frame_bad_G, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_frame_error"
  )

  # 6. Links from nuisance problems with different w, G, or id order on the same count n=14 are refused
  prob_diff_w <- ns$construct_certified_ipw_nuisance_problem(
    id = b14$id, Z = b14$Z, w = rep(1.0, 14L), G = b14$G,
    pair_observed = b14$S, cohorts = c(2L, 5L),
    beta0 = c(0.0, 0.0), F = diag(2L), r = 0.05,
    certificate_control = b14$certificate_control, problem_control = b14$problem_control
  )
  l3_diff_w <- ns$enclose_ipw_link(prob_diff_w, b14$id, b14$y[, 1L], b14$y[, 3L], c(0.0, 0.0), b14$link_control, link_id = "l3")
  expect_error(
    ns$certify_ipw_covariance_transfer(
      list(l1 = b14$links$l1, l2 = b14$links$l2, l3 = l3_diff_w),
      b14$candidate_rows, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_frame_error"
  )

  # 7. D11 separation controls: GMM purpose forbids cluster row_level, cluster_id, or non-unit omega
  expect_error(
    ns$ipw_covariance_operator_spec(
      purpose = "GMM", row_level = "cluster",
      omega = b14$omega_cl, centered = TRUE, kappa = b14$kappa,
      cluster_id = b14$labels_cl, cluster_order = b14$order_cl
    ),
    class = "cdid_ipw_covariance_transfer_d11_error"
  )
  op_gmm_bad_omega <- ns$ipw_covariance_operator_spec(
    purpose = "GMM", row_level = "unit",
    omega = rep(2.0, 14L), centered = TRUE, kappa = b14$kappa
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, b14$frame, b14$C, op_gmm_bad_omega, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_d11_error"
  )
  expect_error(
    ns$ipw_covariance_operator_spec(
      purpose = "GMM", row_level = "unit",
      omega = rep(1.0, 14L), centered = TRUE, kappa = b14$kappa,
      cluster_id = b14$id, cluster_order = b14$id
    ),
    class = "cdid_ipw_covariance_transfer_d11_error"
  )

  # 8. Undeclared cluster IDs or duplicate cluster_order entries are refused
  op_undecl <- ns$ipw_covariance_operator_spec(
    purpose = "inference", row_level = "cluster",
    omega = b14$omega_cl[1:4], centered = TRUE, kappa = b14$kappa,
    cluster_id = b14$labels_cl, cluster_order = b14$order_cl[1:4] # omits "b"
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, b14$candidate_rows, b14$frame, b14$C, op_undecl, b14$V_cand_cl, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_operator_error"
  )

  # 9. Preflight resource bounds precede content scans and callback invocation (callback_calls == 0)
  cb_counter <- 0L
  cand_counted <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      cb_counter <<- cb_counter + 1L
      b14$candidate_rows$callback(a_start, a_end)
    },
    link_ids = c("l1", "l2", "l3")
  )
  ctrl_tiny_work <- ns$ipw_covariance_transfer_control(
    max_chunk_rows = 5L, max_operator_block_rows = 4L, max_cluster_rows = 14L,
    max_matrix_entries = 10000L, max_work_ops = 50L, max_output_entries = 10000L,
    link_control = b14$link_control
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, cand_counted, b14$frame, b14$C, b14$op_cl, b14$V_cand_cl, ctrl_tiny_work
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )
  expect_identical(cb_counter, 0L)

  ctrl_tiny_cl <- ns$ipw_covariance_transfer_control(
    max_chunk_rows = 5L, max_operator_block_rows = 2L, max_cluster_rows = 3L, # < 5 clusters
    max_matrix_entries = 10000L, max_work_ops = 1000000L, max_output_entries = 10000L,
    link_control = b14$link_control
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, cand_counted, b14$frame, b14$C, b14$op_cl, b14$V_cand_cl, ctrl_tiny_cl
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )
  expect_identical(cb_counter, 0L)

  # Character temporary byte budget enforced
  ctrl_tiny_mem <- ns$ipw_covariance_transfer_control(
    max_chunk_rows = 5L, max_operator_block_rows = 2L, max_cluster_rows = 5L,
    max_matrix_entries = 250L, max_work_ops = 1000000L, max_output_entries = 10000L,
    link_control = b14$link_control
  )
  long_id_frame <- ns$ipw_covariance_frame_spec(
    id = paste0(b14$id, strrep("x", 400L)), normalizer_n = 14, embedding = "identity"
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, cand_counted, long_id_frame, b14$C, b14$op_cl, b14$V_cand_cl, ctrl_tiny_mem
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )
  expect_identical(cb_counter, 0L)

  # 10. Callback controls: oversized phi refused before copy, warning refused, wrong ID order refused
  cand_oversized <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      list(
        id = b14$id,
        link_ids = c("l1", "l2", "l3"),
        phi = matrix(0.0, nrow = 14L, ncol = 3L) # 14 > max_chunk_rows = 5
      )
    },
    link_ids = c("l1", "l2", "l3")
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, cand_oversized, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )

  cand_warn <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      warning("suspicious candidate fallback")
      b14$candidate_rows$callback(a_start, a_end)
    },
    link_ids = c("l1", "l2", "l3")
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, cand_warn, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_callback_error"
  )

  cand_bad_order <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      out <- b14$candidate_rows$callback(a_start, a_end)
      out$id <- rev(out$id)
      out
    },
    link_ids = c("l1", "l2", "l3")
  )
  expect_error(
    ns$certify_ipw_covariance_transfer(
      b14$links, cand_bad_order, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
    ),
    class = "cdid_ipw_covariance_transfer_id_error"
  )

  # 11. Honest INCONCLUSIVE propagation:
  # (a) Uncertified sealed 4bd link (with NULL endpoints from sep_try refusal)
  xmax <- .Machine$double.xmax
  l3_inc <- ns$enclose_ipw_link(
    b14$prob3, b14$id, b14$y[, 1L], b14$y[, 3L],
    c(0.0, 0.0), b14$link_control, link_id = "l3",
    candidates = list(mu = -xmax)
  )
  expect_false(l3_inc$certified)
  res_inc_link <- ns$certify_ipw_covariance_transfer(
    list(l1 = b14$links$l1, l2 = b14$links$l2, l3 = l3_inc),
    b14$candidate_rows, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
  )
  expect_false(res_inc_link$certified)
  expect_identical(res_inc_link$status, "INCONCLUSIVE")
  expect_identical(res_inc_link$inconclusive_type, "source_link_inconclusive")
  expect_false(res_inc_link$singularity_proved)
  expect_false(res_inc_link$nonidentification_proved)
  ns$validate_ipw_covariance_transfer(res_inc_link, require_certified = FALSE)
  expect_error(
    ns$validate_ipw_covariance_transfer(res_inc_link, require_certified = TRUE),
    class = "cdid_ipw_covariance_transfer_inconclusive"
  )

  # (b) Overflowing candidate rows in callback chunk -> honest INCONCLUSIVE
  cand_overflow <- ns$ipw_candidate_rows_spec(
    callback = function(a_start, a_end) {
      out <- b14$candidate_rows$callback(a_start, a_end)
      out$phi[1L, 1L] <- xmax
      out
    },
    link_ids = c("l1", "l2", "l3")
  )
  res_inc_chunk <- ns$certify_ipw_covariance_transfer(
    b14$links, cand_overflow, b14$frame, b14$C, b14$op_gmm, b14$V_cand_gmm, b14$control
  )
  expect_false(res_inc_chunk$certified)
  expect_identical(res_inc_chunk$status, "INCONCLUSIVE")
  expect_false(res_inc_chunk$singularity_proved)
  expect_false(res_inc_chunk$nonidentification_proved)
})
