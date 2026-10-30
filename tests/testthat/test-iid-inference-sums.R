# SPDX-License-Identifier: GPL-3.0-only
#
# Deterministic regression tests for the unclustered plain-double row-sum fast
# path in `compute_chained_ipw_inference` and `aggregate_chained_ipw_targets`
# (`R/chained_ipw_pipeline.R`).

resolve_iid_sums_test_env <- function() {
  if (isNamespaceLoaded("cdid")) {
    return(asNamespace("cdid"))
  }
  env <- parent.frame()
  while (!identical(env, emptyenv())) {
    if (exists("compute_chained_ipw_inference", envir = env, inherits = FALSE)) {
      return(env)
    }
    env <- parent.env(env)
  }
  .GlobalEnv
}

make_iid_sums_test_budget <- function() {
  list(
    max_units = 64L,
    max_clusters = 64L,
    max_targets = 16L,
    max_draws = 64L,
    max_if_batch_rows = 16L,
    max_if_batch_cells = 256L,
    max_sign_batch_rows = 16L,
    max_sign_batch_cells = 256L,
    max_workspace_cells = 16384L,
    max_workspace_bytes = 131072L
  )
}

extract_aggregate_center_sums_and_decomposition <- function(agg_fn, n_val, cluster_full_val = NULL) {
  fn_env <- new.env(parent = environment(agg_fn))
  fn_env$n <- as.integer(n_val)
  fn_env$cluster_full <- cluster_full_val
  stmts <- as.list(body(agg_fn))
  hit <- which(vapply(stmts, function(s) {
    is.call(s) && identical(s[[1L]], as.symbol("<-")) &&
      identical(s[[2L]], as.symbol("center_sums"))
  }, logical(1L)))
  if (length(hit) != 1L) {
    stop("Unique `center_sums` definition not found in `aggregate_chained_ipw_targets`.")
  }
  eval(stmts[[hit]], envir = fn_env)
  eval_decomp <- function(X, W, P, D, av, am, M, R, K, labels, ids_all, cohorts, shares_known = TRUE) {
    fn_env$X <- X
    fn_env$W <- W
    fn_env$P <- P
    fn_env$D <- D
    fn_env$av <- av
    fn_env$am <- am
    fn_env$M <- M
    fn_env$R <- R
    fn_env$K <- K
    fn_env$labels <- labels
    fn_env$ids_all <- ids_all
    fn_env$cohorts <- cohorts
    fn_env$shares_known <- shares_known
    eval(stmts[[hit + 1L]], envir = fn_env)
    eval(stmts[[hit + 2L]], envir = fn_env)
    fn_env$decomposition
  }
  list(center_sums = fn_env$center_sums, eval_decomposition = eval_decomp)
}

test_that("compute_chained_ipw_inference unclustered fast path preserves 0.0 + Psi, signed zero, metadata, dispatch, flags, and refusals", {
  env <- resolve_iid_sums_test_env()
  infer_fn <- get("compute_chained_ipw_inference", envir = env)
  supp_fn <- get("ipw_pipeline_supported_inference", envir = env)
  mb <- make_iid_sums_test_budget()

  rng_kind_before <- RNGkind()
  seed_existed_before <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  seed_before <- if (seed_existed_before) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  assert_rng_unchanged <- function() {
    expect_identical(RNGkind(), rng_kind_before)
    expect_identical(exists(".Random.seed", envir = globalenv(), inherits = FALSE), seed_existed_before)
    if (seed_existed_before) {
      expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE), seed_before)
    }
  }

  neg_zero <- -0.0
  expect_identical(1.0 / neg_zero, -Inf)
  expect_identical(1.0 / 0.0, Inf)

  # 1. Tiny matrix (n = 1, J = 1) with -0.0, dimnames, and arbitrary attributes.
  psi_1x1 <- matrix(neg_zero, nrow = 1L, ncol = 1L, dimnames = list("u1", "t1"))
  attr(psi_1x1, "irrelevant_attr") <- list(tag = "discard_me")
  expect_identical(1.0 / psi_1x1[1L, 1L], -Inf)
  res_1x1 <- infer_fn(
    estimates = c(t1 = 0.0),
    if_full = psi_1x1,
    cluster_full = NULL,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expect_identical(res_1x1$raw_cluster_sums, matrix(0.0, 1L, 1L))
  expect_identical(1.0 / res_1x1$raw_cluster_sums, matrix(Inf, 1L, 1L))
  expect_identical(1.0 / res_1x1$cluster_if, matrix(Inf, 1L, 1L))
  expect_identical(attributes(res_1x1$raw_cluster_sums), list(dim = c(1L, 1L)))
  expect_identical(attributes(res_1x1$cluster_if), list(dim = c(1L, 1L)))
  expect_false(res_1x1$clustered)
  expect_identical(res_1x1$n_total, 1L)
  expect_identical(res_1x1$n_clusters, 1L)
  expect_identical(res_1x1$cband_crit_val, 0.0)

  # 2. Tiny matrix (n = 1, J = 3) with mixed -0.0 / +0.0 and attributes.
  psi_1x3 <- matrix(c(neg_zero, 0.0, neg_zero), nrow = 1L, ncol = 3L,
                    dimnames = list("u1", c("t1", "t2", "t3")))
  attr(psi_1x3, "custom_meta") <- "drop_me"
  expect_identical(as.double(1.0 / psi_1x3), c(-Inf, Inf, -Inf))
  res_1x3 <- infer_fn(
    estimates = c(t1 = 0.0, t2 = 0.0, t3 = 0.0),
    if_full = psi_1x3,
    cluster_full = NULL,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expect_identical(res_1x3$raw_cluster_sums, matrix(0.0, 1L, 3L))
  expect_identical(1.0 / res_1x3$raw_cluster_sums, matrix(Inf, 1L, 3L))
  expect_identical(1.0 / res_1x3$cluster_if, matrix(Inf, 1L, 3L))
  expect_identical(attributes(res_1x3$raw_cluster_sums), list(dim = c(1L, 3L)))
  expect_identical(attributes(res_1x3$cluster_if), list(dim = c(1L, 3L)))

  # 3. Tiny matrix (n = 4, J = 1) with -0.0 / +0.0 vs literal formulas and singleton clusters.
  psi_4x1 <- matrix(c(neg_zero, 2.0, 0.0, -2.0), nrow = 4L, ncol = 1L,
                    dimnames = list(paste0("u", 1:4), "t1"))
  attr(psi_4x1, "extra_attr") <- 99L
  expect_identical(unname(1.0 / psi_4x1[c(1L, 3L), 1L]), c(-Inf, Inf))
  res_4x1_iid <- infer_fn(
    estimates = c(t1 = 1.5),
    if_full = psi_4x1,
    cluster_full = NULL,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  res_4x1_cl1 <- infer_fn(
    estimates = c(t1 = 1.5),
    if_full = psi_4x1,
    cluster_full = 1:4,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expected_raw_4x1 <- matrix(c(0.0, 2.0, 0.0, -2.0), nrow = 4L, ncol = 1L)
  expected_V_4x1 <- matrix(8.0 / 16.0, nrow = 1L, ncol = 1L, dimnames = list("t1", "t1"))
  expect_identical(res_4x1_iid$raw_cluster_sums, expected_raw_4x1)
  expect_identical(1.0 / res_4x1_iid$raw_cluster_sums[c(1L, 3L), 1L], c(Inf, Inf))
  expect_identical(1.0 / res_4x1_iid$cluster_if[c(1L, 3L), 1L], c(Inf, Inf))
  expect_identical(attributes(res_4x1_iid$raw_cluster_sums), list(dim = c(4L, 1L)))
  expect_identical(attributes(res_4x1_iid$cluster_if), list(dim = c(4L, 1L)))
  expect_identical(res_4x1_iid$vcov, expected_V_4x1)
  expect_identical(res_4x1_iid$Omega, expected_V_4x1 * 4.0)
  expect_identical(res_4x1_iid$se, c(t1 = sqrt(0.5)))
  expect_identical(res_4x1_iid$raw_cluster_sums, res_4x1_cl1$raw_cluster_sums)
  expect_identical(1.0 / res_4x1_iid$raw_cluster_sums, 1.0 / res_4x1_cl1$raw_cluster_sums)
  expect_identical(1.0 / res_4x1_iid$cluster_if, 1.0 / res_4x1_cl1$cluster_if)
  expect_identical(res_4x1_iid$vcov, res_4x1_cl1$vcov)
  expect_identical(res_4x1_iid$se, res_4x1_cl1$se)
  expect_false(res_4x1_iid$clustered)
  expect_true(res_4x1_cl1$clustered)

  # 4. Tiny matrix (n = 4, J = 3) with -0.0 / +0.0, dimnames, and attributes.
  psi_4x3 <- rbind(
    c(neg_zero, 1.0, -3.0),
    c(2.0, neg_zero, 1.0),
    c(-2.0, 3.0, 0.0),
    c(0.0, -4.0, 2.0)
  )
  dimnames(psi_4x3) <- list(paste0("unit_", 1:4), c("t1", "t2", "t3"))
  attr(psi_4x3, "arbitrary_attr") <- c("a", "b")
  expect_identical(1.0 / psi_4x3[1L, 1L], -Inf)
  expect_identical(1.0 / psi_4x3[4L, 1L], Inf)
  expect_identical(1.0 / psi_4x3[2L, 2L], -Inf)
  expect_identical(1.0 / psi_4x3[3L, 3L], Inf)

  est_3 <- c(t1 = 1.25, t2 = -2.5, t3 = 0.75)
  res_4x3_iid <- infer_fn(
    estimates = est_3,
    if_full = psi_4x3,
    cluster_full = NULL,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  res_4x3_cl1 <- infer_fn(
    estimates = est_3,
    if_full = psi_4x3,
    cluster_full = 1:4,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expected_raw_4x3 <- rbind(
    c(0.0, 1.0, -3.0),
    c(2.0, 0.0, 1.0),
    c(-2.0, 3.0, 0.0),
    c(0.0, -4.0, 2.0)
  )
  expected_V_4x3 <- matrix(
    c(8.0, -6.0, 2.0,
      -6.0, 26.0, -11.0,
      2.0, -11.0, 14.0) / 16.0,
    nrow = 3L, ncol = 3L,
    dimnames = list(names(est_3), names(est_3))
  )
  expected_se_4x3 <- stats::setNames(sqrt(c(8.0, 26.0, 14.0) / 16.0), names(est_3))
  zcrit <- stats::qnorm(0.975)

  expect_identical(res_4x3_iid$raw_cluster_sums, expected_raw_4x3)
  expect_identical(res_4x3_iid$cluster_if, expected_raw_4x3)
  expect_identical(1.0 / res_4x3_iid$raw_cluster_sums[c(1L, 4L), 1L], c(Inf, Inf))
  expect_identical(1.0 / res_4x3_iid$raw_cluster_sums[2L, 2L], Inf)
  expect_identical(1.0 / res_4x3_iid$raw_cluster_sums[3L, 3L], Inf)
  expect_identical(1.0 / res_4x3_iid$cluster_if[c(1L, 4L), 1L], c(Inf, Inf))
  expect_identical(1.0 / res_4x3_iid$cluster_if[2L, 2L], Inf)
  expect_identical(1.0 / res_4x3_iid$cluster_if[3L, 3L], Inf)
  expect_identical(attributes(res_4x3_iid$raw_cluster_sums), list(dim = c(4L, 3L)))
  expect_identical(attributes(res_4x3_iid$cluster_if), list(dim = c(4L, 3L)))
  expect_identical(res_4x3_iid$vcov, expected_V_4x3)
  expect_identical(res_4x3_iid$Omega, expected_V_4x3 * 4.0)
  expect_identical(res_4x3_iid$se, expected_se_4x3)
  expect_identical(res_4x3_iid$analytical_se, unname(expected_se_4x3))
  expect_equal(res_4x3_iid$ci_lower, est_3 - zcrit * expected_se_4x3, tolerance = 1e-15)
  expect_equal(res_4x3_iid$ci_upper, est_3 + zcrit * expected_se_4x3, tolerance = 1e-15)
  for (fld in c("estimates", "se", "analytical_se", "z_stat", "p_value",
                "ci_lower", "ci_upper", "cband_crit_val", "cband_lower",
                "cband_upper", "vcov", "Omega", "cluster_if", "raw_cluster_sums")) {
    expect_identical(res_4x3_iid[[fld]], res_4x3_cl1[[fld]], info = fld)
  }

  # Supported-subvector wrapper preserves NA scatter and unclustered fast path.
  res_supp <- supp_fn(
    estimates = c(est_3, t4 = NA_real_),
    if_full = cbind(psi_4x3, t4 = NA_real_),
    available = c(TRUE, TRUE, TRUE, FALSE),
    cluster_full = NULL,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expect_identical(res_supp$compact$raw_cluster_sums, expected_raw_4x3)
  expect_identical(1.0 / res_supp$raw_cluster_sums[1L, 1L], Inf)
  expect_true(all(is.na(res_supp$raw_cluster_sums[, 4L])))

  # 5. Non-double (integer) and custom-classed matrix inputs use the else loop branch.
  psi_int <- matrix(c(0L, 2L, -2L, 0L, 1L, 0L, 3L, -4L, -3L, 1L, 0L, 2L),
                    nrow = 4L, ncol = 3L, dimnames = dimnames(psi_4x3))
  expect_identical(typeof(psi_int), "integer")
  res_int <- infer_fn(
    estimates = est_3,
    if_full = psi_int,
    cluster_full = NULL,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expect_identical(res_int$raw_cluster_sums, expected_raw_4x3)
  expect_identical(res_int$vcov, expected_V_4x3)

  row_index_calls <- 0L
  as_matrix_custom <- function(x, ...) x
  subset_custom <- function(x, i, j, ..., drop = TRUE) {
    if (!missing(i) && missing(j)) {
      row_index_calls <<- row_index_calls + 1L
    }
    NextMethod("[")
  }
  registerS3method("as.matrix", "cdid_test_custom_mat", as_matrix_custom, envir = environment(infer_fn))
  registerS3method("[", "cdid_test_custom_mat", subset_custom, envir = environment(infer_fn))
  psi_cls <- structure(psi_4x3, class = "cdid_test_custom_mat")
  expect_true(is.object(psi_cls))
  res_cls <- infer_fn(
    estimates = est_3,
    if_full = psi_cls,
    cluster_full = NULL,
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expect_identical(row_index_calls, 4L)
  expect_identical(res_cls$raw_cluster_sums, expected_raw_4x3)
  expect_identical(1.0 / res_cls$raw_cluster_sums[1L, 1L], Inf)
  expect_identical(res_cls$vcov, expected_V_4x3)

  # 6. Fixed supplied sign matrices across inference flag laws (no RNG draw).
  signs_8x4 <- rbind(
    c(-1.0, -1.0,  1.0,  1.0),
    c(-1.0,  1.0, -1.0,  1.0),
    c(-1.0,  1.0,  1.0, -1.0),
    c( 1.0, -1.0, -1.0,  1.0),
    c( 1.0, -1.0,  1.0, -1.0),
    c( 1.0,  1.0, -1.0, -1.0),
    c(-1.0, -1.0, -1.0,  1.0),
    c( 1.0,  1.0,  1.0, -1.0)
  )
  raw_draws <- (signs_8x4 %*% expected_raw_4x3) / sqrt(4.0)
  iqr_norm <- stats::qnorm(0.75) - stats::qnorm(0.25)
  iqr_scale <- apply(raw_draws, 2L, function(col) {
    s <- sort(col)
    (s[ceiling(0.75 * 8L)] - s[ceiling(0.25 * 8L)]) / iqr_norm
  })
  boot_se <- stats::setNames(iqr_scale * sqrt(4.0) / 4.0, names(est_3))
  sup_t <- apply(abs(sweep(raw_draws, 2L, iqr_scale, "/")), 1L, max)
  crit_expected <- sort(sup_t)[ceiling(0.95 * 8L)]

  res_tf <- infer_fn(est_3, psi_4x3, cluster_full = NULL, cband = FALSE, bstrap = TRUE,
                     biters = 8L, fixed_signs = signs_8x4, multiplier_budget = mb,
                     inference_accuracy = NULL, endpoint = "cells")
  res_tt <- infer_fn(est_3, psi_4x3, cluster_full = NULL, cband = TRUE, bstrap = TRUE,
                     biters = 8L, fixed_signs = signs_8x4, multiplier_budget = mb,
                     inference_accuracy = NULL, endpoint = "cells")
  res_ft_cells <- infer_fn(est_3, psi_4x3, cluster_full = NULL, cband = TRUE, bstrap = FALSE,
                           biters = 8L, fixed_signs = signs_8x4, multiplier_budget = mb,
                           inference_accuracy = NULL, endpoint = "cells")
  res_ft_bins <- infer_fn(est_3, psi_4x3, cluster_full = NULL, cband = TRUE, bstrap = FALSE,
                          biters = 8L, fixed_signs = signs_8x4, multiplier_budget = mb,
                          inference_accuracy = NULL, endpoint = "aggte_bins")
  res_tt_cl1 <- infer_fn(est_3, psi_4x3, cluster_full = 1:4, cband = TRUE, bstrap = TRUE,
                         biters = 8L, fixed_signs = signs_8x4, multiplier_budget = mb,
                         inference_accuracy = NULL, endpoint = "cells")

  expect_equal(res_tf$se, boot_se, tolerance = 1e-15)
  expect_identical(res_tf$cband_crit_val, zcrit)
  expect_true(res_tf$bootstrap)
  expect_false(res_tf$simultaneous)
  expect_identical(res_tf$n_draws, 8L)

  expect_equal(res_tt$se, boot_se, tolerance = 1e-15)
  expect_equal(res_tt$cband_crit_val, crit_expected, tolerance = 1e-15)
  expect_true(res_tt$bootstrap)
  expect_true(res_tt$simultaneous)
  expect_identical(res_tt$se, res_tt_cl1$se)
  expect_identical(res_tt$cband_crit_val, res_tt_cl1$cband_crit_val)

  expect_identical(res_ft_cells$se, expected_se_4x3)
  expect_identical(res_ft_cells$cband_crit_val, zcrit)
  expect_null(res_ft_cells$multiplier)
  expect_identical(res_ft_cells$n_draws, NA_integer_)

  expect_identical(res_ft_bins$se, expected_se_4x3)
  expect_equal(res_ft_bins$cband_crit_val, crit_expected, tolerance = 1e-15)
  expect_false(res_ft_bins$bootstrap)
  expect_true(res_ft_bins$simultaneous)
  expect_identical(res_ft_bins$n_draws, 8L)

  # 7. Multi-unit clustered (G < n) path and refusal behavior remain unchanged.
  cl_2 <- c(20L, 10L, 20L, 10L)
  res_cl2 <- infer_fn(
    estimates = est_3,
    if_full = psi_4x3,
    cluster_full = cl_2,
    eligible_full = c(TRUE, FALSE, TRUE, FALSE),
    multiplier_budget = mb,
    inference_accuracy = NULL
  )
  expected_raw_cl2 <- rbind(
    (0.0 + psi_4x3[1L, ]) + psi_4x3[3L, ],
    (0.0 + psi_4x3[2L, ]) + psi_4x3[4L, ]
  )
  dimnames(expected_raw_cl2) <- NULL
  expected_cent_cl2 <- sweep(expected_raw_cl2, 2L, colMeans(expected_raw_cl2), "-")
  expected_V_cl2 <- crossprod(expected_cent_cl2) / 16.0
  dimnames(expected_V_cl2) <- list(names(est_3), names(est_3))
  expect_identical(res_cl2$raw_cluster_sums, expected_raw_cl2)
  expect_identical(res_cl2$cluster_if, expected_cent_cl2)
  expect_identical(res_cl2$vcov, expected_V_cl2)
  expect_true(res_cl2$clustered)
  expect_identical(res_cl2$n_clusters, 2L)
  expect_identical(res_cl2$n_zero_eligible_clusters, 1L)

  # Budget and input refusals
  mb_units <- mb; mb_units$max_units <- 3L
  expect_error(infer_fn(est_3, psi_4x3, multiplier_budget = mb_units, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_bound_error")
  mb_cells <- mb; mb_cells$max_workspace_cells <- 131L # infer_entries = 4*4*3 + 4*4*3 + 4*3*3 = 132
  expect_error(infer_fn(est_3, psi_4x3, multiplier_budget = mb_cells, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_bound_error")
  mb_bytes <- mb; mb_bytes$max_workspace_bytes <- 1055L # 8 * 132 = 1056
  expect_error(infer_fn(est_3, psi_4x3, multiplier_budget = mb_bytes, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_bound_error")
  mb_draws <- mb; mb_draws$max_draws <- 7L
  expect_error(infer_fn(est_3, psi_4x3, cband = TRUE, biters = 8L, fixed_signs = signs_8x4,
                        multiplier_budget = mb_draws, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_bound_error")
  expect_error(infer_fn(est_3[1:2], psi_4x3, multiplier_budget = mb, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_input_error")
  expect_error(infer_fn(est_3, psi_4x3, cluster_full = c(1L, 1L, 1L, 1L),
                        multiplier_budget = mb, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_cluster_error")
  expect_error(infer_fn(est_3, psi_4x3, cluster_full = c(1L, 2L, NA_integer_, 2L),
                        multiplier_budget = mb, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_cluster_error")
  expect_error(infer_fn(est_3, psi_4x3, cband = TRUE, biters = 8L, btype = "mammen",
                        multiplier_budget = mb, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_bootstrap_error")
  expect_error(infer_fn(c(t1 = 1.0), matrix(0.0, 4L, 1L), multiplier_budget = mb, inference_accuracy = NULL),
               class = "cdid_ipw_pipeline_scale_error")

  assert_rng_unchanged()
})

test_that("aggregate_chained_ipw_targets local center_sums fast path preserves 0.0 + z, signed zero, metadata, dispatch, and Omega decomposition", {
  env <- resolve_iid_sums_test_env()
  agg_fn <- get("aggregate_chained_ipw_targets", envir = env)

  rng_kind_before <- RNGkind()
  seed_existed_before <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  seed_before <- if (seed_existed_before) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL

  neg_zero <- -0.0
  # 1. Tiny matrices (n = 1, J = 1), (n = 1, J = 3), (n = 4, J = 1), (n = 4, J = 3)
  ext_1 <- extract_aggregate_center_sums_and_decomposition(agg_fn, 1L, NULL)
  cs_1 <- ext_1$center_sums
  z_1x1 <- matrix(neg_zero, 1L, 1L, dimnames = list("u1", "c1"))
  attr(z_1x1, "extra") <- "drop"
  out_1x1 <- cs_1(z_1x1)
  expect_identical(out_1x1, matrix(0.0, 1L, 1L))
  expect_identical(1.0 / out_1x1, matrix(Inf, 1L, 1L))
  expect_identical(attributes(out_1x1), list(dim = c(1L, 1L)))

  z_1x3 <- matrix(c(neg_zero, 0.0, neg_zero), 1L, 3L, dimnames = list("u1", c("a", "b", "c")))
  attr(z_1x3, "extra") <- 1L
  out_1x3 <- cs_1(z_1x3)
  expect_identical(1.0 / out_1x3, matrix(Inf, 1L, 3L))
  expect_identical(attributes(out_1x3), list(dim = c(1L, 3L)))

  ext_4_iid <- extract_aggregate_center_sums_and_decomposition(agg_fn, 4L, NULL)
  ext_4_cl1 <- extract_aggregate_center_sums_and_decomposition(agg_fn, 4L, 1:4)
  ext_4_cl2 <- extract_aggregate_center_sums_and_decomposition(agg_fn, 4L, c(1L, 1L, 2L, 2L))
  cs_4_iid <- ext_4_iid$center_sums
  cs_4_cl1 <- ext_4_cl1$center_sums
  cs_4_cl2 <- ext_4_cl2$center_sums

  z_4x1 <- matrix(c(neg_zero, 2.0, 0.0, -2.0), 4L, 1L, dimnames = list(paste0("u", 1:4), "c1"))
  attr(z_4x1, "meta") <- list(x = 1)
  out_4x1_iid <- cs_4_iid(z_4x1)
  out_4x1_cl1 <- cs_4_cl1(z_4x1)
  expect_identical(out_4x1_iid, out_4x1_cl1)
  # In IEEE-754, (-0.0) - (+0.0) is -0.0, whereas ((+0.0) + (-0.0)) - (+0.0) is +0.0.
  expect_identical(1.0 / out_4x1_iid[c(1L, 3L), 1L], c(Inf, Inf))
  expect_identical(1.0 / out_4x1_cl1[c(1L, 3L), 1L], c(Inf, Inf))
  expect_identical(attributes(out_4x1_iid), list(dim = c(4L, 1L)))

  X_4x3 <- rbind(
    c(neg_zero,  1.0, -2.0),
    c( 2.0, neg_zero,  0.0),
    c(-2.0, -1.0,  2.0),
    c( 0.0,  0.0, neg_zero)
  )
  dimnames(X_4x3) <- list(paste0("u", 1:4), c("a", "b", "c"))
  attr(X_4x3, "custom") <- "strip"
  xc_iid <- cs_4_iid(X_4x3)
  xc_cl1 <- cs_4_cl1(X_4x3)
  expect_identical(xc_iid, xc_cl1)
  expect_identical(1.0 / xc_iid[c(1L, 4L), 1L], c(Inf, Inf))
  expect_identical(1.0 / xc_iid[c(2L, 4L), 2L], c(Inf, Inf))
  expect_identical(1.0 / xc_iid[c(2L, 4L), 3L], c(Inf, Inf))
  expect_identical(attributes(xc_iid), list(dim = c(4L, 3L)))

  # Integer and custom-classed inputs take the else branch in center_sums.
  X_int <- matrix(as.integer(X_4x3), 4L, 3L, dimnames = dimnames(X_4x3))
  expect_identical(cs_4_iid(X_int), xc_iid)
  cs_row_calls <- 0L
  subset_cs <- function(x, i, j, ..., drop = TRUE) {
    if (!missing(i) && missing(j)) cs_row_calls <<- cs_row_calls + 1L
    NextMethod("[")
  }
  registerS3method("[", "cdid_test_cs_mat", subset_cs, envir = environment(agg_fn))
  X_cls <- structure(X_4x3, class = "cdid_test_cs_mat")
  expect_identical(cs_4_iid(X_cls), xc_iid)
  expect_identical(cs_row_calls, 4L)

  # 2. Complete share/cross IF decomposition AST block (unclustered vs singleton-clustered vs multi-unit clustered).
  th <- c(a = 2.0, b = 5.0, c = 8.0)
  gs <- c(3, 3, 4, Inf)
  ws <- c(1.0, 3.0, 2.0, 2.0) # sum = 8, mean = 2; p_3 = 0.5, p_4 = 0.25
  p <- c("3" = 0.5, "4" = 0.25)
  P <- cbind(
    "3" = ws * (as.double(gs == 3) - p[1L]) / mean(ws),
    "4" = ws * (as.double(gs == 4) - p[2L]) / mean(ws)
  )
  den <- 1.25
  W <- matrix(c(0.4, 0.4, 0.2), nrow = 1L, dimnames = list("ATT_simple", names(th)))
  est_simple <- sum(W * th) # 4.4
  D <- matrix(c((th[1L] - est_simple) + (th[2L] - est_simple), th[3L] - est_simple) / den,
              nrow = 1L, dimnames = list("ATT_simple", c("3", "4")))
  fixed_if <- X_4x3 %*% t(W)
  share_if <- P %*% t(D)
  total_if <- fixed_if + share_if

  dec_iid <- ext_4_iid$eval_decomposition(X_4x3, W, P, D, 1:3, 1L, 1L, 3L, 2L,
                                          "ATT_simple", names(th), c(3, 4), TRUE)
  dec_cl1 <- ext_4_cl1$eval_decomposition(X_4x3, W, P, D, 1:3, 1L, 1L, 3L, 2L,
                                          "ATT_simple", names(th), c(3, 4), TRUE)
  dec_cl2 <- ext_4_cl2$eval_decomposition(X_4x3, W, P, D, 1:3, 1L, 1L, 3L, 2L,
                                          "ATT_simple", names(th), c(3, 4), TRUE)
  expect_identical(dec_iid, dec_cl1)

  # Independent total-IF variance identity:
  tc_iid <- cs_4_iid(total_if)
  fc_iid <- cs_4_iid(fixed_if)
  sc_iid <- cs_4_iid(share_if)
  expect_equal(unname(dec_iid$Omega_total), crossprod(tc_iid) / 4.0, tolerance = 1e-15)
  expect_equal(unname(dec_iid$Omega_fixed), crossprod(fc_iid) / 4.0, tolerance = 1e-15)
  expect_equal(unname(dec_iid$Omega_share), crossprod(sc_iid) / 4.0, tolerance = 1e-15)
  expect_equal(unname(dec_iid$Omega_cross), (crossprod(fc_iid, sc_iid) + crossprod(sc_iid, fc_iid)) / 4.0, tolerance = 1e-15)
  expect_gt(abs(dec_iid$Omega_cross[1L, 1L]), 0.1)

  # Clustered (G = 2) decomposition identity:
  tc_cl2 <- cs_4_cl2(total_if)
  expect_equal(unname(dec_cl2$Omega_total), crossprod(tc_cl2) / 4.0, tolerance = 1e-15)

  expect_identical(RNGkind(), rng_kind_before)
  expect_identical(exists(".Random.seed", envir = globalenv(), inherits = FALSE), seed_existed_before)
  if (seed_existed_before) {
    expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE), seed_before)
  }
})

test_that("aggregate_chained_ipw_targets unclustered vs clustered full entrypoint share and cross IF decomposition", {
  env <- resolve_iid_sums_test_env()
  agg_fn <- get("aggregate_chained_ipw_targets", envir = env)
  mb <- make_iid_sums_test_budget()

  neg_zero <- -0.0
  th <- c(a = 2.0, b = 5.0, c = 8.0)
  X <- rbind(
    c(neg_zero,  1.0, -2.0),
    c( 2.0, neg_zero,  0.0),
    c(-2.0, -1.0,  2.0),
    c( 0.0,  0.0, neg_zero)
  )
  dimnames(X) <- list(paste0("u", 1:4), names(th))
  attr(X, "irrelevant") <- "drop_me"
  gs <- c(3, 3, 4, Inf)
  ws <- c(1.0, 3.0, 2.0, 2.0)
  meta <- data.frame(
    g = c(3, 3, 4),
    time = c(3, 4, 4),
    role = "post",
    target_id = names(th),
    support_status = "available",
    stringsAsFactors = FALSE
  )

  res_iid <- agg_fn(
    theta_hat = th, target_if_full = X, target_meta = meta,
    full_frame_g = gs, full_frame_w = ws, type = "simple",
    cluster_full = NULL, multiplier_budget = mb, inference_accuracy = NULL
  )
  res_cl1 <- agg_fn(
    theta_hat = th, target_if_full = X, target_meta = meta,
    full_frame_g = gs, full_frame_w = ws, type = "simple",
    cluster_full = 1:4, multiplier_budget = mb, inference_accuracy = NULL
  )
  expect_identical(res_iid$estimates, res_cl1$estimates)
  expect_identical(res_iid$if_full, res_cl1$if_full)
  expect_identical(res_iid$if_fixed_full, res_cl1$if_fixed_full)
  expect_identical(res_iid$if_share_full, res_cl1$if_share_full)
  expect_identical(res_iid$decomposition, res_cl1$decomposition)
  expect_identical(res_iid$se, res_cl1$se)
  expect_identical(res_iid$se_fixed_shares_only, res_cl1$se_fixed_shares_only)
  expect_equal(unname(res_iid$decomposition$Omega_total), 4.0 * unname(res_iid$inference$vcov), tolerance = 1e-15)
  expect_gt(abs(res_iid$decomposition$Omega_cross[1L, 1L]), 0.1)
})
