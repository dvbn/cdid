# SPDX-License-Identifier: GPL-3.0-only
#
# Focused end-to-end tests for CP02 Chained IPW-GMM Pipeline Integration
# (`R/chained_ipw_pipeline.R`) and the callback `srcref` dispatch repair
# (`R/chained_ipw_covariance_view.R:1212`, `R/chained_ipw_covariance_transfer.R:2184`).

parse_test_rational <- function(x) {
  vapply(as.character(x), function(s) {
    parts <- strsplit(s, "/", fixed = TRUE)[[1L]]
    if (length(parts) == 1L) as.double(parts[1L]) else as.double(parts[1L]) / as.double(parts[2L])
  }, numeric(1L), USE.NAMES = FALSE)
}

parse_test_rational_matrix <- function(m) {
  out <- apply(m, c(1L, 2L), parse_test_rational)
  storage.mode(out) <- "double"
  out
}

make_test_pipeline_control <- function(
  require_matching_scale = TRUE,
  max_target_error_abs = 1e-4,
  max_loading_error_abs = 1e-4,
  max_work_ops = 500000000L
) {
  lnk_ctrl <- ipw_link_enclosure_control(
    max_chunk_rows = 64L,
    max_retained_entries = 2000000L,
    max_matrix_entries = 2000000L,
    max_work_ops = max_work_ops,
    max_output_entries = 500000L
  )
  ipw_pipeline_control(
    owner_control = ipw_owner_control(
      max_chunk_rows = 64L,
      max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops,
      max_output_entries = 500000L
    ),
    catalog_bounds = list(
      max_candidates = 500L,
      max_rows = 500L,
      max_pool_entries = 5000L,
      max_scratch = 100000L,
      max_work = max_work_ops
    ),
    class_bounds = list(
      max_links = 500L,
      max_classes = 50L,
      max_class_dates = 50L,
      max_E_entries = 5000L,
      max_mask_entries = 100000L
    ),
    fitted_control = ipw_fitted_view_control(
      max_chunk_rows = 64L,
      max_observer_units = 2000L,
      max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops,
      max_output_entries = 500000L
    ),
    point_control = ipw_point_view_control(
      max_chunk_rows = 64L,
      max_observer_units = 2000L,
      max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops,
      max_output_entries = 500000L
    ),
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    certificate_control = logit_certificate_control(
      taylor_degree = 12L,
      max_abs_eta = 50.0,
      max_squarings = 16L,
      max_chunk_rows = 64L,
      max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops
    ),
    problem_control = ipw_nuisance_problem_control(
      max_chunk_rows = 64L,
      max_retained_entries = 2000000L,
      max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops
    ),
    link_control = lnk_ctrl,
    covariance_control = ipw_covariance_transfer_control(
      max_chunk_rows = 64L,
      max_operator_block_rows = 64L,
      max_cluster_rows = 2000L,
      max_matrix_entries = 2000000L,
      max_work_ops = max_work_ops,
      max_output_entries = 1000000L,
      link_control = lnk_ctrl
    ),
    gmm_accuracy = gmm_accuracy_control(
      method = "certified_eigen_chol_kkt",
      target_error_unit = "abs",
      max_target_error_abs = max_target_error_abs,
      max_target_error_se = 1e-2,
      max_orthogonality_defect = 1e-8,
      min_certified_lambda_rel = 1e-10,
      max_factor_residual_rel = 1e-10,
      min_rcond_R = 1e-10,
      min_rcond_U = 1e-10,
      min_sigma_rel_A = 1e-10,
      max_kkt_residual_rel = 1e-10,
      max_loading_residual_rel = 1e-10,
      require_matching_scale = require_matching_scale,
      max_loading_error_abs = max_loading_error_abs
    ),
    gmm_bounds = gmm_bounds_control(
      max_matrix_entries = 2000000L,
      max_batch_rows = 64L,
      max_targets = 100L,
      max_open_clusters = 2000L,
      allow_materialize_if = TRUE,
      max_materialized_if_entries = 1000000L
    ),
    multiplier_budget = list(max_units=2000,max_clusters=2000,max_targets=100,max_draws=2000,
      max_if_batch_rows=64,max_if_batch_cells=6400,max_sign_batch_rows=32,
      max_sign_batch_cells=64000,max_workspace_cells=1000000,max_workspace_bytes=8000000),
    inference_accuracy = list(max_estimate_error_abs=1e-3,max_se_error_abs=1e-2,
      max_vcov_error_abs=1e-2,max_multiplier_error_abs=1e-1,max_band_error_abs=1e-1)
  )
}

resolve_test_cdid_env <- function() {
  if (isNamespaceLoaded("cdid")) {
    return(asNamespace("cdid"))
  }
  env <- parent.frame()
  while (!identical(env, emptyenv())) {
    if (exists("run_chained_ipw_pipeline", envir = env, inherits = FALSE)) {
      return(env)
    }
    env <- parent.env(env)
  }
  .GlobalEnv
}

locate_test_fixture <- function(filename) {
  candidates <- c(
    tryCatch(testthat::test_path("fixtures", filename), error = function(e) character(0L)),
    file.path("fixtures", filename),
    file.path("tests", "testthat", "fixtures", filename),
    file.path("..", "..", "tests", "testthat", "fixtures", filename)
  )
  for (p in candidates) {
    if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  }
  stop("Fixture file not found: ", filename)
}

locate_optional_workspace_file <- function(rel_path) {
  candidates <- c(
    rel_path,
    file.path("..", "..", rel_path),
    file.path("..", "..", "..", rel_path)
  )
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  }
  NULL
}

test_that("callback srcref metadata is accepted on both function and spec-list forms while forbidden shortcut attributes are refused", {
  fixture_path <- locate_test_fixture("ipw_covariance_view_fixtures.R")
  fx_builders <- source(fixture_path, local = TRUE, keep.source = TRUE)$value
  fx <- fx_builders$build_14unit_3link_view_and_dense(ns = resolve_test_cdid_env())

  # Construct callbacks explicitly parsed with keep.source = TRUE (carrying `srcref`) and keep.source = FALSE
  cb_text <- "function(a_start, a_end) { idx <- seq.int(a_start, a_end); list(id = fx$id[idx], link_ids = c('l1', 'l2', 'l3'), phi = fx$phi_cand[idx, , drop = FALSE], a_index = idx) }"
  cb_with_srcref <- eval(parse(text = cb_text, keep.source = TRUE))
  cb_without_srcref <- eval(parse(text = cb_text, keep.source = FALSE))

  expect_false(is.null(attr(cb_with_srcref, "srcref", exact = TRUE)))
  expect_true(is.null(attr(cb_without_srcref, "srcref", exact = TRUE)))

  spec_with_srcref <- ipw_candidate_rows_spec(
    callback = cb_with_srcref,
    link_ids = c("l1", "l2", "l3")
  )
  spec_without_srcref <- ipw_candidate_rows_spec(
    callback = cb_without_srcref,
    link_ids = c("l1", "l2", "l3")
  )

  # 1. Compact view entrypoint (`certify_ipw_covariance_transfer_view`)
  res_view_fn_src <- certify_ipw_covariance_transfer_view(
    links = fx$links_view,
    candidate_rows = cb_with_srcref,
    frame_owner = fx$frame_owner,
    C = fx$C,
    op = fx$op_gmm,
    V_candidate = fx$V_cand_gmm,
    control = fx$cov_ctrl,
    reader_control = fx$o_ctrl
  )
  res_view_fn_nosrc <- certify_ipw_covariance_transfer_view(
    links = fx$links_view,
    candidate_rows = cb_without_srcref,
    frame_owner = fx$frame_owner,
    C = fx$C,
    op = fx$op_gmm,
    V_candidate = fx$V_cand_gmm,
    control = fx$cov_ctrl,
    reader_control = fx$o_ctrl
  )
  res_view_spec_src <- certify_ipw_covariance_transfer_view(
    links = fx$links_view,
    candidate_rows = spec_with_srcref,
    frame_owner = fx$frame_owner,
    C = fx$C,
    op = fx$op_gmm,
    V_candidate = fx$V_cand_gmm,
    control = fx$cov_ctrl,
    reader_control = fx$o_ctrl
  )
  res_view_spec_nosrc <- certify_ipw_covariance_transfer_view(
    links = fx$links_view,
    candidate_rows = spec_without_srcref,
    frame_owner = fx$frame_owner,
    C = fx$C,
    op = fx$op_gmm,
    V_candidate = fx$V_cand_gmm,
    control = fx$cov_ctrl,
    reader_control = fx$o_ctrl
  )

  expect_true(res_view_fn_src$certified)
  expect_equal(res_view_fn_src$epsilon_V, res_view_fn_nosrc$epsilon_V, tolerance = 0)
  expect_equal(res_view_fn_src$epsilon_V, res_view_spec_src$epsilon_V, tolerance = 0)
  expect_equal(res_view_fn_src$epsilon_V, res_view_spec_nosrc$epsilon_V, tolerance = 0)

  # 2. Dense entrypoint (`certify_ipw_covariance_transfer`)
  res_dense_fn_src <- certify_ipw_covariance_transfer(
    links = fx$links_dense,
    candidate_rows = cb_with_srcref,
    frame = fx$frame_dense,
    C = fx$C,
    op = fx$op_gmm,
    V_candidate = fx$V_cand_gmm,
    control = fx$cov_ctrl
  )
  res_dense_fn_nosrc <- certify_ipw_covariance_transfer(
    links = fx$links_dense,
    candidate_rows = cb_without_srcref,
    frame = fx$frame_dense,
    C = fx$C,
    op = fx$op_gmm,
    V_candidate = fx$V_cand_gmm,
    control = fx$cov_ctrl
  )
  res_dense_spec_src <- certify_ipw_covariance_transfer(
    links = fx$links_dense,
    candidate_rows = spec_with_srcref,
    frame = fx$frame_dense,
    C = fx$C,
    op = fx$op_gmm,
    V_candidate = fx$V_cand_gmm,
    control = fx$cov_ctrl
  )
  expect_true(res_dense_fn_src$certified)
  expect_equal(res_dense_fn_src$epsilon_V, res_dense_fn_nosrc$epsilon_V, tolerance = 0)
  expect_equal(res_dense_fn_src$epsilon_V, res_dense_spec_src$epsilon_V, tolerance = 0)
  expect_equal(res_view_fn_src$epsilon_V, res_dense_fn_src$epsilon_V, tolerance = 0)

  # 3. Forbidden shortcut attributes are still refused on both function and spec-list forms
  cb_bad <- cb_with_srcref
  attr(cb_bad, "certified") <- TRUE
  expect_error(
    certify_ipw_covariance_transfer_view(
      links = fx$links_view,
      candidate_rows = cb_bad,
      frame_owner = fx$frame_owner,
      C = fx$C,
      op = fx$op_gmm,
      V_candidate = fx$V_cand_gmm,
      control = fx$cov_ctrl,
      reader_control = fx$o_ctrl
    ),
    class = "cdid_ipw_covariance_transfer_metadata_error"
  )
  expect_error(
    certify_ipw_covariance_transfer(
      links = fx$links_dense,
      candidate_rows = cb_bad,
      frame = fx$frame_dense,
      C = fx$C,
      op = fx$op_gmm,
      V_candidate = fx$V_cand_gmm,
      control = fx$cov_ctrl
    ),
    class = "cdid_ipw_covariance_transfer_metadata_error"
  )

  spec_bad <- spec_with_srcref
  attr(spec_bad, "extra_attr") <- "forbidden"
  expect_error(
    certify_ipw_covariance_transfer_view(
      links = fx$links_view,
      candidate_rows = spec_bad,
      frame_owner = fx$frame_owner,
      C = fx$C,
      op = fx$op_gmm,
      V_candidate = fx$V_cand_gmm,
      control = fx$cov_ctrl,
      reader_control = fx$o_ctrl
    ),
    class = "cdid_ipw_covariance_transfer_metadata_error"
  )

  # 4. Also verify the four bare-callback predicates recorded in `context/installed_gate_failure.md`
  #    (lines 402, 592, 727, 789 of `test-chained-ipw-covariance-view.R`):
  w96 <- fx_builders$build_weighted96_shared_mask_view(resolve_test_cdid_env(), chunk_rows = 24L, operator_block_rows = 16L)
  op_96 <- ipw_covariance_operator_spec(
    purpose = "GMM", row_level = "unit", omega = rep(1.0, w96$n), centered = TRUE, kappa = 1.0 / w96$n
  )
  res_w96_3 <- certify_ipw_covariance_transfer_view(
    links = w96$links_3,
    candidate_rows = function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(phi = w96$phi_mat[idx, , drop = FALSE], id = w96$id[idx], link_ids = c("l12", "l23", "l13"), a_index = idx)
    },
    frame_owner = w96$frame_owner, C = diag(1.0, 3L), op = op_96,
    V_candidate = w96$V_3, control = w96$cov_ctrl, reader_control = w96$o_ctrl
  )
  expect_true(res_w96_3$enclosure_certified)
  expect_false(res_w96_3$certified)
  expect_identical(res_w96_3$status, "INCONCLUSIVE")

  res_w96_2 <- certify_ipw_covariance_transfer_view(
    links = w96$links_2,
    candidate_rows = function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(phi = w96$phi_mat[idx, 1:2, drop = FALSE], id = w96$id[idx], link_ids = c("l12", "l23"), a_index = idx)
    },
    frame_owner = w96$frame_owner, C = diag(1.0, 2L), op = op_96,
    V_candidate = w96$V_2, control = w96$cov_ctrl, reader_control = w96$o_ctrl
  )
  expect_true(res_w96_2$certified)
  expect_true(res_w96_2$pd_transferred)

  tight_cl_ctrl <- ipw_covariance_transfer_control(
    max_chunk_rows = 5L, max_operator_block_rows = 4L, max_cluster_rows = 3L,
    max_matrix_entries = 15000L, max_work_ops = 1000000L, max_output_entries = 10000L,
    link_control = fx$link_ctrl
  )
  expect_error(
    certify_ipw_covariance_transfer_view(
      links = fx$links_view,
      candidate_rows = function(s, t) fx$callback(s, t),
      frame_owner = fx$frame_owner, C = fx$C, op = fx$op_cl, V_candidate = fx$V_cand_cl,
      control = tight_cl_ctrl, reader_control = fx$o_ctrl
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )

  expect_error(
    certify_ipw_covariance_transfer_view(
      links = fx$links_view,
      candidate_rows = function(s, t) {
        x <- fx$callback(s, t)
        x$a_index <- x$a_index + 1L
        x
      },
      frame_owner = fx$frame_owner, C = fx$C, op = fx$op_gmm, V_candidate = fx$V_cand_gmm,
      control = fx$cov_ctrl, reader_control = fx$o_ctrl
    ),
    class = "cdid_ipw_covariance_transfer_callback_error"
  )

  tight_work <- fx$cov_ctrl
  tight_work$max_work_ops <- floor(res_view_spec_src$preflight$combined_total_work_ops_bound - 500)
  expect_error(
    certify_ipw_covariance_transfer_view(
      links = fx$links_view,
      candidate_rows = function(s, t) fx$callback(s, t),
      frame_owner = fx$frame_owner, C = fx$C, op = fx$op_gmm, V_candidate = fx$V_cand_gmm,
      control = tight_work, reader_control = fx$o_ctrl
    ),
    class = "cdid_ipw_covariance_transfer_bound_error"
  )
})

test_that("rotating80 exact binary-representable fixture matches reference equations and certifies complete Theorem 12 transfer", {
  fx80 <- dget(locate_test_fixture("rotating80_fixture.R"))
  df80 <- fx80$panel

  # When the workspace root `reference/rotating80/` directory is present on disk,
  # cross-check that `rotating80_fixture.R` matches the raw reference files identically.
  raw_p80 <- locate_optional_workspace_file("reference/rotating80/panel.csv")
  if (!is.null(raw_p80)) {
    expect_equal(df80, utils::read.csv(raw_p80, stringsAsFactors = FALSE))
  }

  ctrl_opt <- make_test_pipeline_control(require_matching_scale = TRUE)
  ctrl_id  <- make_test_pipeline_control(require_matching_scale = FALSE)

  res_exact <- run_chained_ipw_pipeline(
    data = df80, idname = "id", tname = "t", yname = "y_exactfit", gname = "g",
    xformla = ~ x, weightsname = "weight", control = ctrl_opt, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", weighting = "optimal",
    run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6), cband = TRUE, biters = 200L, seed = 42L
  )

  res_pert_opt <- run_chained_ipw_pipeline(
    data = df80, idname = "id", tname = "t", yname = "y_pairperturbed", gname = "g",
    xformla = ~ x, weightsname = "weight", control = ctrl_opt, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", weighting = "optimal",
    run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6), cband = TRUE, biters = 200L, seed = 42L
  )

  res_pert_id <- run_chained_ipw_pipeline(
    data = df80, idname = "id", tname = "t", yname = "y_pairperturbed", gname = "g",
    xformla = ~ x, weightsname = "weight", control = ctrl_id, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", weighting = "identity",
    run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6), cband = TRUE, biters = 200L, seed = 42L
  )

  Omega_ref <- fx80$Omega_ref
  Sigma_ref <- fx80$Sigma_ref
  C_opt_ref <- fx80$C_opt_ref
  m_exact_ref <- fx80$m_exact_ref
  theta_exact_ref <- fx80$theta_exact_ref
  m_pert_ref <- fx80$m_pert_ref
  theta_pert_opt_ref <- fx80$theta_pert_opt_ref
  theta_pert_id_ref <- fx80$theta_pert_id_ref

  # Five disjoint observer masks forbid structural redundancy reduction
  expect_false(res_exact$reduction_summary$reduction_active)
  expect_equal(res_exact$reduction_summary$n_links_all, 5L)
  expect_equal(res_exact$reduction_summary$n_links_used, 5L)

  # Point moments, covariance, optimal & identity target estimates, and target IF covariance
  expect_equal(unname(res_exact$covariance_stream$mhat), m_exact_ref, tolerance = 1e-12)
  expect_equal(unname(res_exact$covariance_stream$V), Omega_ref, tolerance = 1e-12)
  expect_equal(unname(res_exact$reporting$estimates[res_exact$reporting$target_meta$role != "reference"]), theta_exact_ref, tolerance = 1e-12)
  expect_equal(unname(res_exact$reporting$Omega_a[res_exact$reporting$target_meta$role != "reference", res_exact$reporting$target_meta$role != "reference"]), Sigma_ref, tolerance = 1e-12)

  expect_equal(unname(res_pert_opt$covariance_stream$mhat), m_pert_ref, tolerance = 1e-12)
  expect_equal(unname(res_pert_opt$reporting$estimates[res_pert_opt$reporting$target_meta$role != "reference"]), theta_pert_opt_ref, tolerance = 1e-12)
  expect_equal(unname(res_pert_id$reporting$estimates[res_pert_id$reporting$target_meta$role != "reference"]), theta_pert_id_ref, tolerance = 1e-12)
  expect_equal(unname(res_pert_opt$reporting$Omega_a[res_pert_opt$reporting$target_meta$role != "reference", res_pert_opt$reporting$target_meta$role != "reference"]), Sigma_ref, tolerance = 1e-12)

  # Optimal and identity post-treatment ATTs differ on pair_perturbed
  expect_gt(min(abs(res_pert_opt$reporting$estimates[res_pert_opt$reporting$target_meta$role == "post"] - res_pert_id$reporting$estimates[res_pert_id$reporting$target_meta$role == "post"])), 0.05)

  # Verify unit-level IF rows against `reference/rotating80/if_rows.csv`
  # Note: in `cdid_chained_design`, parameter 1 is `level(1) - level(2) = -m_1`,
  # while post-treatment parameters 2..4 are `level(t) - level(2)` (`t = 3, 4, 5`).
  phi_ref_mat <- fx80$phi_ref_mat
  target_if_ref <- phi_ref_mat %*% t(C_opt_ref)
  target_if_ref[, 1L] <- -target_if_ref[, 1L]
  expect_equal(unname(res_exact$reporting$if_full[,res_exact$reporting$target_meta$role != "reference",drop=FALSE]), unname(target_if_ref), tolerance = 1e-12)

  # Complete Theorem 12 (`thm:gmm`) numerical transfer certificates
  expect_true(res_exact$numerical_transfer$certified)
  expect_true(res_pert_opt$numerical_transfer$certified)
  expect_true(res_pert_id$numerical_transfer$certified)

  gt_opt <- res_pert_opt$numerical_transfer$gmm_transfer
  gt_id  <- res_pert_id$numerical_transfer$gmm_transfer
  expect_true(all(is.finite(gt_opt$total_target_error_cert)))
  expect_true(all(gt_opt$total_target_error_cert < 1e-4))
  expect_true(all(is.finite(gt_opt$total_se_error_cert)))
  expect_true(all(gt_opt$total_se_error_cert < 1e-4))
  # Genuine formula discriminator: optimal weighting has positive changing-covariance loading error `delta_C_stat > 0`,
  # whereas identity weighting has `delta_C_stat == 0`.
  expect_true(all(gt_opt$delta_C_stat > 0))
  expect_true(all(gt_id$delta_C_stat == 0))
  expect_true(all(gt_opt$total_loading_error_cert > gt_opt$delta_C_num))
})

test_that("panel24 deterministic fitted reference matches rational moments and post-ATT targets", {
  df24 <- utils::read.csv(locate_test_fixture("fitted_panel24.csv"), stringsAsFactors = FALSE)
  ctrl_opt <- make_test_pipeline_control(require_matching_scale = TRUE)

  res24 <- run_chained_ipw_pipeline(
    data = df24, idname = "id", tname = "time", yname = "y", gname = "g",
    xformla = ~ 1, control = ctrl_opt, target_cohorts = 3,
    control_group = "nevertreated", base_period = "universal", weighting = "optimal",
    run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6)
  )

  expect_equal(unname(res24$covariance_stream$mhat), c(-1.5, -0.5, 1.0, 1.0, 0.0, 0.0), tolerance = 1e-12)
  expect_equal(unname(res24$reporting$estimates[res24$reporting$target_meta$role != "reference"]), c(1.5, 1.0, 1.0, 1.0), tolerance = 1e-12)
  expect_true(res24$numerical_transfer$certified)
})

test_that("shared-mask structural redundancy reduction, D13 not-yet-treated reuse, anticipation, base_period, clustering, and aggregations work end-to-end", {
  set.seed(20261016L)
  n_u <- 48L
  u_ids <- seq_len(n_u)
  u_g <- rep(c(3, 4, Inf), each = 16L)
  u_x <- rep(c(-1.0, 1.0, -0.5, 0.5), times = 12L)
  u_w <- rep(c(1.0, 1.5, 0.8, 1.2), times = 12L)
  u_cl <- rep(seq_len(12L), each = 4L)
  times_grid <- 1:4

  panel_list <- vector("list", n_u * length(times_grid))
  idx_r <- 0L
  for (i in seq_len(n_u)) {
    for (tt in times_grid) {
      idx_r <- idx_r + 1L
      eff <- if (u_g[i] == 3 && tt >= 3) 2.0 + 0.5 * (tt - 3) else if (u_g[i] == 4 && tt >= 4) 1.5 else 0.0
      y_val <- 0.25 * i + 0.5 * tt + 0.3 * u_x[i] * tt + eff + 0.2 * rnorm(1L)
      panel_list[[idx_r]] <- data.frame(
        id = u_ids[i], t = tt, g = u_g[i], x = u_x[i], weight = u_w[i], cluster = u_cl[i], y = y_val,
        stringsAsFactors = FALSE
      )
    }
  }
  df_bal <- do.call(rbind, panel_list)

  # Full frame has 52 units: 48 panel units + 4 frame-only singleton units in cluster 13
  f_id_full <- c(u_ids, 49:52)
  f_Z_full <- stats::model.matrix(~ x, data = data.frame(x = c(u_x, -1.0, 1.0, -0.5, 0.5)))
  f_w_full <- c(u_w, 1.0, 1.0, 1.0, 1.0)
  f_g_full <- c(u_g, 3, 4, Inf, Inf)
  cl_full <- c(u_cl, rep(13L, 4L))

  ctrl_opt <- make_test_pipeline_control(require_matching_scale = TRUE)
  ctrl_id  <- make_test_pipeline_control(require_matching_scale = FALSE)

  # 1. Optimal weighting with structural redundancy reduction vs unreduced cycle singularity refusal
  res_red <- run_chained_ipw_pipeline(
    data = df_bal, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, weightsname = "weight", cluster_vec = cl_full,
    frame_id = f_id_full, frame_Z = f_Z_full, frame_weights = f_w_full, frame_g = f_g_full,
    allow_frame_only = TRUE, control = ctrl_opt,
    control_group = "notyettreated", base_period = "varying", weighting = "optimal",
    reduce_structural_redundancy = TRUE, run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6),
    cband = TRUE, biters = 200L, seed = 99L
  )

  expect_true(res_red$reduction_summary$reduction_active)
  expect_lt(res_red$reduction_summary$n_links_used, res_red$reduction_summary$n_links_all)
  # CP08 P3: built-in optimal IPW certifies its classes outcome-free and never fits
  # dropped rows, so the fitted-versus-expanded self-test has no value on this
  # route; the accepted post-fit path still evaluates it on the same retained rows.
  expect_true(is.na(res_red$reduction_summary$max_moment_expansion_diff))
  expect_identical(res_red$reduction_summary$n_rows_fitted, res_red$reduction_summary$n_links_used)
  res_post <- run_chained_ipw_pipeline(
    data = df_bal, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, weightsname = "weight", cluster_vec = cl_full,
    frame_id = f_id_full, frame_Z = f_Z_full, frame_weights = f_w_full, frame_g = f_g_full,
    allow_frame_only = TRUE, control = ctrl_opt,
    control_group = "notyettreated", base_period = "varying", weighting = "optimal",
    reduce_structural_redundancy = TRUE, run_numerical_transfer = FALSE,
    structural_reduction_stage = "post_point_fit"
  )
  expect_lt(res_post$reduction_summary$max_moment_expansion_diff, 1e-12)
  expect_identical(res_post$reduction_summary$kept_indices, res_red$reduction_summary$kept_indices)
  expect_identical(res_post$gmm_fit$targets$estimate, res_red$gmm_fit$targets$estimate)
  expect_true(res_red$numerical_transfer$certified)

  # Unreduced optimal GMM on cyclic shared-mask panel has rank-deficient covariance and is refused
  expect_error(
    run_chained_ipw_pipeline(
      data = df_bal, idname = "id", tname = "t", yname = "y", gname = "g",
      xformla = ~ x, weightsname = "weight", cluster_vec = cl_full,
      frame_id = f_id_full, frame_Z = f_Z_full, frame_weights = f_w_full, frame_g = f_g_full,
      allow_frame_only = TRUE, control = ctrl_opt,
      control_group = "notyettreated", base_period = "varying", weighting = "optimal",
      reduce_structural_redundancy = FALSE, run_numerical_transfer = FALSE
    ),
    class = "cdid_gmm_pd_error"
  )

  # 2. Identity weighting keeps all rows even when `reduce_structural_redundancy = TRUE`
  res_id <- run_chained_ipw_pipeline(
    data = df_bal, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, weightsname = "weight", cluster_vec = cl_full,
    frame_id = f_id_full, frame_Z = f_Z_full, frame_weights = f_w_full, frame_g = f_g_full,
    allow_frame_only = TRUE, control = ctrl_id,
    control_group = "notyettreated", base_period = "varying", weighting = "identity",
    reduce_structural_redundancy = TRUE, run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6)
  )
  expect_true(res_id$numerical_transfer$gmm_transfer$complete_inference_certified)
  expect_true(res_id$numerical_transfer$covariance_transfer$enclosure_certified)
  expect_false(res_id$numerical_transfer$covariance_transfer$pd_transferred)
  expect_false(res_id$reduction_summary$reduction_active)
  expect_equal(res_id$reduction_summary$n_links_used, res_id$reduction_summary$n_links_all)

  # 3. `base_period = "universal"` vs `"varying"` reporting map
  res_univ <- run_chained_ipw_pipeline(
    data = df_bal, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, weightsname = "weight", cluster_vec = cl_full,
    frame_id = f_id_full, frame_Z = f_Z_full, frame_weights = f_w_full, frame_g = f_g_full,
    allow_frame_only = TRUE, control = ctrl_opt,
    control_group = "notyettreated", base_period = "universal", weighting = "optimal",
    reduce_structural_redundancy = TRUE, run_numerical_transfer = FALSE
  )
  post_mask <- res_red$reporting$target_meta$role == "post"
  # Post-treatment estimates and IFs are identical under `varying` and `universal`
  expect_equal(res_red$reporting$estimates[post_mask], res_univ$reporting$estimates[res_univ$reporting$target_meta$role == "post"], tolerance = 1e-12)
  # For cohort g = 4, pre-treatment period t = 1 has base_t = 3 under universal and base_t = 2 under varying
  g4_t2 <- which(res_red$reporting$target_meta$g == 4 & res_red$reporting$target_meta$time == 2)
  u4_t2 <- which(res_univ$reporting$target_meta$g == 4 & res_univ$reporting$target_meta$time == 2)
  expect_true(length(g4_t2)>0L)
  expect_true(all(res_univ$reporting$target_meta$base_t[u4_t2] == 3))
  expect_true(all(res_red$reporting$target_meta$base_t[g4_t2] == 1))
  expect_false(any(res_red$reporting$target_meta$time == 1))

  # 4. Anticipation (`anticipation = 1L`) shifts the reference boundary and D13 control pool
  res_ant1 <- run_chained_ipw_pipeline(
    data = df_bal, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, weightsname = "weight",
    control = ctrl_opt, anticipation = 1L,
    control_group = "nevertreated", base_period = "universal", weighting = "optimal",
    reduce_structural_redundancy = TRUE, run_numerical_transfer = FALSE
  )
  expect_true(all(res_ant1$reporting$target_meta$ref <= res_ant1$reporting$target_meta$g - 2))

  # 5. Scale invariance under positive constant rescaling of sampling weights (w -> 3 * w)
  df_bal_3w <- df_bal
  df_bal_3w$weight <- 3.0 * df_bal$weight
  res_3w <- run_chained_ipw_pipeline(
    data = df_bal_3w, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, weightsname = "weight", cluster_vec = cl_full,
    frame_id = f_id_full, frame_Z = f_Z_full, frame_weights = 3.0 * f_w_full, frame_g = f_g_full,
    allow_frame_only = TRUE, control = ctrl_opt,
    control_group = "notyettreated", base_period = "varying", weighting = "optimal",
    reduce_structural_redundancy = TRUE, run_numerical_transfer = FALSE, cband = TRUE, biters = 200L, seed = 99L
  )
  expect_equal(res_red$reporting$estimates, res_3w$reporting$estimates, tolerance = 1e-12)
  expect_equal(res_red$reporting$se, res_3w$reporting$se, tolerance = 1e-12)
  expect_equal(res_red$aggregations$simple$estimates, res_3w$aggregations$simple$estimates, tolerance = 1e-12)
  expect_equal(res_red$aggregations$simple$se, res_3w$aggregations$simple$se, tolerance = 1e-12)

  # 6. Genuine formula-fault discriminators on full-frame transport, clustering, and cohort-share aggregation
  # (a) Full-frame scale factor is n_total / n_A = 52 / 48, and unscaled embedding would change V(theta_hat) by (48/52)^2
  expect_equal(res_red$full_frame_transport$scale_factor, 52.0 / 48.0)
  Phi_A <- res_red$gmm_fit$if_projection$materialized_target_if
  Phi_full <- res_red$full_frame_transport$Phi_full
  expect_equal(
    crossprod(Phi_full) / (52.0^2),
    crossprod(Phi_A) / (48.0^2),
    tolerance = 1e-12
  )
  # Frame-only rows (49:52) have zero target IF but nonzero cohort-share IF
  expect_true(all(Phi_full[49:52, ] == 0.0))
  expect_true(any(res_red$aggregations$simple$if_share_full[49:52, ] != 0.0))
  # Zero-eligible cluster 13 is retained in cluster count
  expect_equal(res_red$reporting$inference$n_clusters, 13L)
  expect_equal(res_red$reporting$inference$n_zero_eligible_clusters, 1L)
  # Nonzero cross term between target IFs and cohort-share IFs in aggregation decomposition
  expect_gt(max(abs(res_red$aggregations$simple$decomposition$Omega_cross)), 1e-6)
  expect_false(isTRUE(all.equal(
    res_red$aggregations$simple$se,
    res_red$aggregations$simple$se_fixed_shares_only,
    tolerance = 1e-6
  )))
})

test_that("gapped/nonadjacent panels without any balanced subset or adjacent pairs are estimated via chained comparisons", {
  # Construct a panel over t in {1, 3, 5} where no unit observes all 3 dates:
  # half the units observe only (1, 3) and half observe only (3, 5).
  set.seed(20261017L)
  n_u <- 32L
  rows_list <- list()
  for (i in seq_len(n_u)) {
    g_i <- if (i <= 16L) 5 else Inf
    x_i <- if (i %% 2L == 0L) 1.0 else -1.0
    obs_dates <- if (i %% 4L %in% c(1L, 2L)) c(1, 3) else c(3, 5)
    for (tt in obs_dates) {
      eff <- if (g_i == 5 && tt >= 5) 2.5 else 0.0
      y_val <- 0.1 * i + 0.4 * tt + 0.2 * x_i * tt + eff + 0.15 * rnorm(1L)
      rows_list[[length(rows_list) + 1L]] <- data.frame(
        id = i, t = tt, g = g_i, x = x_i, y = y_val, stringsAsFactors = FALSE
      )
    }
  }
  df_gap <- do.call(rbind, rows_list)
  ctrl_opt <- make_test_pipeline_control(require_matching_scale = TRUE)

  res_gap <- run_chained_ipw_pipeline(
    data = df_gap, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, control = ctrl_opt, target_cohorts = 5,
    control_group = "nevertreated", base_period = "universal", weighting = "optimal",
    run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6)
  )

  expect_true(res_gap$numerical_transfer$certified)
  expect_equal(nrow(res_gap$reporting$target_meta), 3L)
  expect_true(all(is.finite(res_gap$reporting$estimates)))
  expect_true(all(res_gap$reporting$se[res_gap$reporting$target_meta$role != "reference"] > 0))
})

test_that("non-estimable requested targets, singular covariance, and accuracy/resource budget violations are explicitly refused", {
  df80 <- dget(locate_test_fixture("rotating80_fixture.R"))$panel
  ctrl_opt <- make_test_pipeline_control(require_matching_scale = TRUE)

  # 1. Non-estimable requested post target (t = 5 is in calendar via g = Inf, but g = 3 has no observations at t = 5)
  #    Contract change (approved D5, October 16 2026): the transitional private entry point refused this
  #    whole fit with class `cdid_gmm_target_estimability_error` ("Requested post-treatment target(s) are
  #    not structurally estimable from the qualified comparison graph: (3, 5) [isolated_unsupported]").
  #    That refusal was support-only: the two supported cells (3,3) and (3,4) and their joint GMM were
  #    discarded because one requested cell lacked catalog support. D5 keeps the declared ledger, returns
  #    the unsupported cell as NA with a catalog-relative reason (never a population non-identification
  #    claim), and preserves every supported level, IF, covariance and certificate. Genuine refusals
  #    (singular covariance, accuracy, resources) below are unchanged.
  df80_no_g3_t5 <- df80[!(df80$g == 3 & df80$t == 5), ]
  res_partial <- run_chained_ipw_pipeline(
    data = df80_no_g3_t5, idname = "id", tname = "t", yname = "y_exactfit", gname = "g",
    xformla = ~ x, control = ctrl_opt, target_cohorts = 3,
    requested_post = data.frame(g = 3, t = 5),
    target_post_times = data.frame(g = 3, t = c(3, 4, 5)),
    control_group = "nevertreated", base_period = "universal", weighting = "optimal",
    run_numerical_transfer = FALSE
  )
  expect_identical(res_partial$support_status, "partial_support")
  led <- res_partial$support$ledger
  expect_identical(led$cell_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)", "ATT(g=3,t=5)"))
  expect_identical(led$support_status, c("estimable", "estimable", "isolated_unsupported"))
  expect_false(any(grepl("nonidentif", led$support_reason)))
  post <- which(res_partial$reporting$target_meta$role == "post")
  expect_equal(unname(res_partial$reporting$estimates[post]), c(1.0, 1.0, NA_real_), tolerance = 1e-12)
  expect_true(all(is.na(res_partial$reporting$if_full[, post[3L]])))
  expect_true(all(is.finite(res_partial$reporting$if_full[, post[1:2]])))
  expect_identical(res_partial$gmm_fit$targets$status[res_partial$gmm_fit$targets$target_id == "ATT(g=3,t=5)"], "isolated_unsupported")
  expect_true(is.na(res_partial$aggregations$simple$estimates[["ATT_simple"]]))
  expect_match(res_partial$aggregations$simple$cell_status$reason, "ATT\\(g=3,t=5\\)")
  expect_false(res_partial$support$whole_request_complete)

  # 2. Singular moment covariance under optimal weighting (constant outcomes within each cohort/pair -> V = 0)
  df_zero_var <- df80
  df_zero_var$y_const <- ifelse(df_zero_var$g == 3 & df_zero_var$t >= 3, 1.0, 0.0)
  expect_error(
    run_chained_ipw_pipeline(
      data = df_zero_var, idname = "id", tname = "t", yname = "y_const", gname = "g",
      xformla = ~ 1, control = ctrl_opt, target_cohorts = 3,
      control_group = "nevertreated", base_period = "universal", weighting = "optimal",
      run_numerical_transfer = FALSE
    ),
    class = "cdid_gmm_pd_error"
  )

  # 3. Deliberate accuracy-budget refusal (`max_target_error_abs = 1e-18` refuses Theorem 12 transfer bound)
  ctrl_tight_acc <- make_test_pipeline_control(
    require_matching_scale = TRUE,
    max_target_error_abs = 1e-18
  )
  expect_error(
    run_chained_ipw_pipeline(
      data = df80, idname = "id", tname = "t", yname = "y_pairperturbed", gname = "g",
      xformla = ~ x, weightsname = "weight", control = ctrl_tight_acc, target_cohorts = 3,
      control_group = "nevertreated", base_period = "universal", weighting = "optimal",
      run_numerical_transfer = TRUE, root_witness = function(fv) list(F=unname(fv$factors$H_chol),r=1e-6)
    ),
    class = "cdid_gmm_accuracy_error"
  )

  # 4. Deliberate resource-bound refusal (`max_work_ops = 100L`)
  ctrl_tight_work <- make_test_pipeline_control(
    require_matching_scale = TRUE,
    max_work_ops = 100L
  )
  expect_error(
    run_chained_ipw_pipeline(
      data = df80, idname = "id", tname = "t", yname = "y_exactfit", gname = "g",
      xformla = ~ x, weightsname = "weight", control = ctrl_tight_work, target_cohorts = 3,
      control_group = "nevertreated", base_period = "universal", weighting = "optimal",
      run_numerical_transfer = FALSE
    ),
    class = "cdid_ipw_owner_bound_error"
  )
})

test_that("declared fixed signs and semidefinite inference preserve positive marginals and zero-eligible clusters", {
  ctrl <- make_test_pipeline_control()
  X <- rbind(c(-2,4),c(-1,2),c(1,-2),c(2,-4),c(0,0))
  th <- c(a=1,b=2); clusters <- c(1L,1L,2L,2L,3L)
  signs <- unname(as.matrix(expand.grid(rep(list(c(-1.0,1.0)),3L))))
  storage.mode(signs) <- "double"
  enc <- list(theta=list(lo=matrix(th,1L),hi=matrix(th,1L)),rows=list(lo=X,hi=X))
  before_kind <- RNGkind(); set.seed(123L); before_seed <- .Random.seed
  z <- compute_chained_ipw_inference(th,X,clusters,cband=TRUE,biters=8L,
    multiplier_budget=ctrl$multiplier_budget,inference_accuracy=ctrl$inference_accuracy,
    fixed_signs=signs,enclosure=enc,eligible_full=c(TRUE,TRUE,TRUE,TRUE,FALSE))
  S <- rbind(c(-3,6),c(3,-6),c(0,0))
  raw <- signs %*% S / sqrt(3.0)
  scale <- apply(raw,2L,function(x)diff(stats::quantile(x,c(.25,.75),type=1,names=FALSE)))/
    diff(stats::qnorm(c(.25,.75)))
  bse <- scale*sqrt(3.0)/5.0
  critical <- stats::quantile(apply(abs(sweep(raw,2L,scale,"/")),1L,max),.95,type=1,names=FALSE)
  expect_equal(unname(z$vcov),crossprod(S)/25.0,tolerance=1e-15)
  expect_equal(unname(z$se),bse,tolerance=1e-15)
  expect_equal(z$cband_crit_val,critical,tolerance=1e-15)
  expect_equal(unname(z$cband_lower),th-critical*bse,tolerance=1e-15,ignore_attr=TRUE)
  expect_true(all(z$analytical_se>0))
  expect_equal(qr(z$vcov)$rank,1L)
  expect_equal(z$n_zero_eligible_clusters,1L)
  expect_true(z$numerical_certificate$certified)
  expect_identical(.Random.seed,before_seed)
  expect_identical(RNGkind(),before_kind)
  expect_error(compute_chained_ipw_inference(th,X,clusters,cband=TRUE,biters=8L,btype="normal",
    multiplier_budget=ctrl$multiplier_budget,inference_accuracy=ctrl$inference_accuracy),
    class="cdid_ipw_pipeline_bootstrap_error")
  bad <- ctrl$multiplier_budget;bad$max_draws<-7
  expect_error(compute_chained_ipw_inference(th,X,clusters,cband=TRUE,biters=8L,fixed_signs=signs,
    multiplier_budget=bad,inference_accuracy=ctrl$inference_accuracy),class="cdid_ipw_pipeline_bound_error")
})

test_that("released bstrap = FALSE, cband = TRUE bin law with declared fixed signs: analytic SE certificate, multiplier scales apart, native interval product", {
  ctrl <- make_test_pipeline_control()
  X <- rbind(c(-2,4),c(-1,2),c(1,-2),c(2,-4),c(0,0))
  th <- c(a=1,b=2); clusters <- c(1L,1L,2L,2L,3L)
  signs <- unname(as.matrix(expand.grid(rep(list(c(-1.0,1.0)),3L))))
  storage.mode(signs) <- "double"
  enc <- list(theta=list(lo=matrix(th,1L),hi=matrix(th,1L)),rows=list(lo=X,hi=X))
  args <- list(th,X,clusters,cband=TRUE,biters=8L,multiplier_budget=ctrl$multiplier_budget,
    inference_accuracy=ctrl$inference_accuracy,fixed_signs=signs,enclosure=enc,eligible_full=c(TRUE,TRUE,TRUE,TRUE,FALSE))
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  tt <- do.call(compute_chained_ipw_inference, args)                                   # private law (bstrap = NULL)
  ft <- do.call(compute_chained_ipw_inference, c(args, list(bstrap=FALSE, endpoint="aggte_bins")))
  ff <- do.call(compute_chained_ipw_inference, c(args[names(args) != "cband"], list(cband=FALSE, bstrap=FALSE)))
  expect_false(exists(".Random.seed", envir = globalenv()))
  S <- rbind(c(-3,6),c(3,-6),c(0,0))
  raw <- signs %*% S / sqrt(3.0)
  scale <- apply(raw,2L,function(x)diff(stats::quantile(x,c(.25,.75),type=1,names=FALSE)))/diff(stats::qnorm(c(.25,.75)))
  critical <- stats::quantile(apply(abs(sweep(raw,2L,scale,"/")),1L,max),.95,type=1,names=FALSE)
  # analytic reported SEs (getSE), sup-t constant from the draws' own IQR scale (mboot), band = est +/- crit * analytic SE
  expect_equal(unname(ft$se),sqrt(colSums(S^2))/5.0,tolerance=1e-15)
  expect_identical(ft$se, ff$se); expect_identical(unname(ft$se), unname(ft$analytical_se))
  expect_equal(ft$cband_crit_val,critical,tolerance=1e-15); expect_identical(ft$cband_crit_val, tt$cband_crit_val)
  expect_equal(unname(ft$cband_lower),th-critical*sqrt(colSums(S^2))/5.0,tolerance=1e-15,ignore_attr=TRUE)
  expect_false(ft$bootstrap); expect_true(ft$simultaneous); expect_true(tt$bootstrap)
  expect_identical(ft$flag_law, "released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band")
  # certificate: each reported quantity against its own law, both scales retained, budgets unchanged
  cf <- ft$numerical_certificate
  expect_true(cf$certified); expect_identical(names(cf$errors), names(ff$numerical_certificate$errors))
  expect_identical(cf$se, ff$numerical_certificate$se)
  expect_identical(cf$critical_value, tt$numerical_certificate$critical_value)
  expect_true(all(cf$multiplier_scale[1L,] <= scale & scale <= cf$multiplier_scale[2L,]))
  pure <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "chained_gmm_core_pure_r_reference.R"), envir = pure)
  factor <- pure$gmm_iv_div(pure$gmm_iv_sqrt(pure$gmm_iv(3.0)), pure$gmm_iv(5.0))
  expected_se_interval <- vapply(seq_len(ncol(cf$multiplier_scale)), function(j)
    pure$gmm_iv_mul(cf$multiplier_scale[, j], factor), numeric(2L))
  expect_identical(expected_se_interval, tt$numerical_certificate$se)
  expect_identical(cf$errors[["max_multiplier_error_abs"]], tt$numerical_certificate$errors[["max_multiplier_error_abs"]])
  expect_true(all(cf$band$lo <= ft$cband_lower & ft$cband_lower <= cf$band$hi))
  expect_null(tt$numerical_certificate$multiplier_scale); expect_null(ff$numerical_certificate$multiplier_scale)
  # the certificate's sign-by-cluster products run on the native kernel, bit-identical to the frozen R recurrence
  ref_env <- new.env(parent = baseenv())
  sys.source(testthat::test_path("reference", "cp02_interval_product_R_reference.R"), envir = ref_env)
  C <- ft$numerical_certificate$cluster
  nat <- ipw_pipeline_interval_product(signs, signs, C$lo, C$hi)
  ref <- ref_env$ipw_pipeline_interval_product(signs, signs, C$lo, C$hi)
  bits <- function(x) writeBin(as.double(x), raw(), size = 8L, endian = "little")
  expect_identical(bits(nat$lo), bits(ref$lo)); expect_identical(bits(nat$hi), bits(ref$hi))
  expect_true(all(nat$lo <= raw * sqrt(3.0) & raw * sqrt(3.0) <= nat$hi))
})

test_that("weighted aggregation uses full-frame share IF and cross terms under the same cluster operator", {
  ctrl <- make_test_pipeline_control()
  th <- c(a=2,b=4,c=7)
  X <- rbind(c(1,-1,2),c(-1,2,-1),c(2,0,-2),c(-2,-1,1),c(0,0,0),c(0,0,0))
  gs <- c(3,3,4,4,Inf,3); ws <- c(1,2,3,1,2,4)
  meta <- data.frame(g=c(3,3,4),time=c(3,4,4),role="post",target_id=names(th))
  cl <- c(1L,1L,2L,2L,3L,3L)
  z <- aggregate_chained_ipw_targets(th,X,meta,gs,ws,type="simple",cluster_full=cl,
    multiplier_budget=ctrl$multiplier_budget,inference_accuracy=ctrl$inference_accuracy)
  mass <- c(sum(ws[gs==3]),sum(ws[gs==4]));p<-mass/sum(ws)
  w<-p[c(1,1,2)]/sum(p[c(1,1,2)]);T<-sum(w*th)
  D<-c(sum(th[1:2]-T),th[3]-T)/sum(p[c(1,1,2)])
  P<-cbind(ws*(as.double(gs==3)-p[1]),ws*(as.double(gs==4)-p[2]))/mean(ws)
  expected<-as.double(X%*%w+P%*%D)
  sums<-rbind(colSums(cbind(X,P)[1:2,,drop=FALSE]),colSums(cbind(X,P)[3:4,,drop=FALSE]),colSums(cbind(X,P)[5:6,,drop=FALSE]))
  centered<-sweep(sums,2L,colMeans(sums),"-")
  total<-as.double(centered%*%c(w,D))
  expect_equal(z$cohort_shares,stats::setNames(p,c("3","4")),tolerance=1e-15)
  expect_equal(z$if_full[,1L],expected,tolerance=1e-14)
  expect_equal(unname(z$inference$vcov),matrix(sum(total^2)/36.0),tolerance=1e-14)
  expect_equal(unname(z$decomposition$Omega_total),6*z$inference$vcov,tolerance=1e-14,ignore_attr=TRUE)
  expect_gt(max(abs(z$decomposition$Omega_cross)),.01)
})

test_that("direct did-style auxiliary pretrends can use early-clean cohorts without changing post GMM", {
  ctrl <- make_test_pipeline_control(require_matching_scale=FALSE)
  gs <- rep(c(3,4,Inf),each=16L); v <- rep(c(-1.0,1.0),24L)
  dat <- expand.grid(id=1:48,t=1:4)
  dat$g <- gs[dat$id]
  dat$y <- as.double(dat$g==3 & dat$t==2)*10.0 + .125*v[dat$id]*dat$t
  args <- list(data=dat,idname="id",tname="t",yname="y",gname="g",control=ctrl,
    target_cohorts=c(3,4),control_group="notyettreated",weighting="identity",
    root_witness=function(fv)list(F=unname(fv$factors$H_chol),r=1e-6),run_numerical_transfer=TRUE)
  varying <- do.call(run_chained_ipw_pipeline,c(args,list(base_period="varying")))
  universal <- do.call(run_chained_ipw_pipeline,c(args,list(base_period="universal")))
  rv <- which(varying$reporting$target_meta$g==4 & varying$reporting$target_meta$time==2)
  ru <- which(universal$reporting$target_meta$g==4 & universal$reporting$target_meta$time==2)
  # Accepted released did2.5.1 reporting oracle: g4 report2/from1 uses controls
  # g3 and never; universal report2/from3 uses only never. See supplied receipt.
  expect_equal(unname(varying$reporting$estimates[rv]),-5.0,tolerance=1e-12)
  expect_equal(unname(universal$reporting$estimates[ru]),0.0,tolerance=1e-12)
  direct_if <- ifelse(gs==4,.375*v,ifelse(gs==3,-7.5-.1875*v,7.5-.1875*v))
  expect_equal(unname(varying$reporting$if_full[,rv]),direct_if,tolerance=1e-12)
  descriptor <- varying$auxiliary_reporting$selections$varying
  spec <- descriptor[[which(vapply(descriptor,function(x)x$g==4 && x$t==2,logical(1L)))]]
  expect_equal(spec$K,c(3,Inf))
  expect_true(any(!varying$metadata$catalog$rows$original))
  expect_equal(varying$covariance_stream,universal$covariance_stream,tolerance=0)
  expect_equal(varying$gmm_fit$targets,universal$gmm_fit$targets,tolerance=0)
  expect_equal(varying$gmm_fit$internal_gauge$Lambda,universal$gmm_fit$internal_gauge$Lambda,tolerance=0)
  expect_equal(varying$gmm_fit$if_projection$materialized_target_if,
    universal$gmm_fit$if_projection$materialized_target_if,tolerance=0)
  expect_equal(varying$auxiliary_reporting$post_GMM_rows_added,0L)
  expect_true(varying$numerical_transfer$gmm_transfer$complete_inference_certified)
  expect_true(varying$numerical_transfer$covariance_transfer$enclosure_certified)
  expect_false(varying$numerical_transfer$covariance_transfer$pd_transferred)
  expect_true(all(varying$reporting$numerical_certificate$se[1L,]>0))
  expect_true(all(direct_if>=varying$reporting$numerical_certificate$rows$lo[,rv] &
    direct_if<=varying$reporting$numerical_certificate$rows$hi[,rv]))
})

test_that("explicit frames preserve known singleton cohorts and refuse partial frame input", {
  ctrl <- make_test_pipeline_control()
  d <- dget(locate_test_fixture("rotating80_fixture.R"))$panel
  one <- d[1L,,drop=FALSE];one$id<-20000L;one$x<-0.0;one$weight<-1.0
  d <- rbind(d,one)
  unit <- d[!duplicated(d$id),,drop=FALSE]
  Z <- stats::model.matrix(~x,unit)
  z <- run_chained_ipw_pipeline(d,"id","t","y_exactfit","g",ctrl,xformla=~x,
    frame_id=unit$id,frame_Z=Z,frame_weights=unit$weight,frame_g=NULL,
    target_cohorts=3,control_group="nevertreated",base_period="universal",run_numerical_transfer=FALSE)
  expect_equal(z$full_frame_transport$n_total,81L)
  expect_true(20000L %in% z$frame_owner$singletons$id)
  expect_equal(unname(z$aggregations$simple$cohort_shares),sum(unit$weight[unit$g==3])/sum(unit$weight),tolerance=1e-15)
  expect_error(run_chained_ipw_pipeline(d,"id","t","y_exactfit","g",ctrl,frame_id=unit$id,
    run_numerical_transfer=FALSE),class="cdid_ipw_pipeline_input_error")
})
