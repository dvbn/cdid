# Private D52 Chained Pure-Imputation backend tests.
# Uses locked installed cdid namespace functions and real sealed owners
# (panel, frame, reader, design_owner, link views) constructed from the frozen
# independent reference CSVs (tests/testthat/fixtures/imputation_finite_v1).

imp_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, envir = ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}

imp_fx_dir <- testthat::test_path("fixtures", "imputation_finite_v1")
imp_loader <- new.env(parent = asNamespace("cdid"))
sys.source(file.path(imp_fx_dir, "load_fixture.R"), envir = imp_loader)
sys.source(file.path(imp_fx_dir, "diagnostics.R"), envir = imp_loader)

imp_specs <- imp_loader$imp_fixture_specs
imp_perm_order <- c(seq.int(36L, 2L, -2L), seq.int(35L, 1L, -2L))
imp_mx <- function(x) max(c(0, abs(x)))

# Inherited CP06 QR refit convention, fixed before this replay.
# See tests/testthat/test-chained-or-wls.R:30-36 for the refit convention; this is not a public accuracy guarantee. Independent-reference
# thresholds, point comparisons, and covariance comparisons remain unchanged.
imp_qr_refit_difference <- function(actual, expected) {
  imp_mx(actual - expected) / max(1, imp_mx(expected))
}

imp_load <- function(fixture_id, order = NULL, max_chunk_rows = 13L) {
  imp_loader$load_imputation_finite_reference(
    imp_fx_dir,
    fixture_id = fixture_id,
    order = order,
    max_chunk_rows = max_chunk_rows
  )
}

imp_fit <- function(fx,
                    est_method = "dr",
                    imputation_pooling = fx$imputation_pooling,
                    views = fx$views,
                    control = fx$control,
                    glm_control = fx$glm_control,
                    rank_tol = fx$rank_tol,
                    outcome_recipe = fx$outcome_recipe,
                    treatment_columns = fx$treatment_columns,
                    design_owner = fx$design_owner,
                    X_frame = NULL,
                    pair_outcomes = NULL,
                    selection_xformla = NULL,
                    selection_model = NULL,
                    selection_pooling = NULL) {
  imp_ns("fit_imputation_backend")(
    reader = fx$reader,
    link_views = views,
    control = control,
    glm_control = glm_control,
    rank_tol = rank_tol,
    outcome_recipe = outcome_recipe,
    treatment_columns = treatment_columns,
    est_method = est_method,
    imputation_pooling = imputation_pooling,
    design_owner = design_owner,
    X_frame = X_frame,
    pair_outcomes = pair_outcomes,
    selection_xformla = selection_xformla,
    selection_model = selection_model,
    selection_pooling = selection_pooling
  )
}

imp_collect <- function(fit, chunk_rows = 11L, link_batch_size = 2L,
                        scale = "frame", details = TRUE) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  N <- st$N
  L <- length(fit$links)
  phi <- matrix(NA_real_, N, L)
  phi_frame <- matrix(NA_real_, N, L)
  batches <- list()

  info <- imp_ns("replay_imputation_if")(
    fit,
    chunk_rows = chunk_rows,
    link_batch_size = link_batch_size,
    callback = function(b) {
      stopifnot(identical(b$id, st$id[b$rows]))
      phi[b$rows, b$links] <<- b$phi
      phi_frame[b$rows, b$links] <<- if (scale == "frame") b$phi else b$phi / st$lambda
      batches[[length(batches) + 1L]] <<- b
    },
    scale = scale,
    details = details
  )
  list(phi = phi, phi_frame = phi_frame, info = info, batches = batches)
}

imp_cluster_cov <- function(phi_frame, cluster_id, n_A) {
  N <- nrow(phi_frame)
  L <- ncol(phi_frame)
  clusters <- sort(unique(cluster_id))
  G_cl <- matrix(0, length(clusters), L)
  for (i in seq_along(clusters)) {
    rows <- which(cluster_id == clusters[i])
    G_cl[i, ] <- colSums(phi_frame[rows, , drop = FALSE])
  }
  V_unit <- crossprod(phi_frame) / (N^2)
  V_cl <- crossprod(G_cl) / (N^2)
  Omega_cl <- N * V_cl
  # Core psi = (n_A / N) * phi_frame with 1 / n_A^2 cluster covariance normalization
  G_core <- (n_A / N) * G_cl
  V_core_cl <- crossprod(G_core) / (n_A^2)
  list(
    unit_coefficient = V_unit,
    cluster_coefficient = V_cl,
    cluster_sqrtN = Omega_cl,
    core_cluster_coefficient = V_core_cl
  )
}

testthat::test_that("packaged imputation reference hashes and masked panel atoms hold", {
  imp_loader$imp_fixture_verify_hashes(imp_fx_dir)
  fx <- imp_load("F1_sep_int_x1x2")
  testthat::expect_identical(fx$panel$n_a_units, 24L)
  testthat::expect_identical(fx$reader$state$N, 36L)
  testthat::expect_identical(fx$reader$state$n_A, 24L)
  testthat::expect_true(all(c(10, 22, 34) %in% fx$panel$a_id))
  testthat::expect_true(all(is.na(fx$Y_obs[fx$R == 0])))
  imp_ns("validate_imputation_design_owner")(fx$design_owner, fx$reader)
})

testthat::test_that("all 9 fixtures and 108 method-links match points, nuisances, IF rows, 8-channel components, predictions, projections, and covariances", {
  total_method_links <- 0L
  total_nuisance <- 0L
  total_uif <- 0L
  total_comp <- 0L
  total_cov <- 0L
  total_pred <- 0L
  total_proj <- 0L
  total_proj_if <- 0L

  # The fixed projection is read from the independently pinned CSV below.

  for (fid in imp_specs$fixture_id) {
    fx <- imp_load(fid, max_chunk_rows = 11L)
    exp <- fx$expectations
    N <- fx$reader$state$N
    n_A <- fx$reader$state$n_A
    link_ids <- fx$catalog$links$link_id

    # 1. Fit real sealed-owner backend for all 3 methods on exact Delta_e AND on panel Y_obs
    fits_delta <- list()
    fits_y <- list()
    for (m in c("ipw", "reg", "dr")) {
      fits_delta[[m]] <- imp_fit(fx, est_method = m, pair_outcomes = fx$pair_outcomes_delta)
      imp_ns("validate_imputation_backend")(fits_delta[[m]])
      expected_groups <- switch(fx$imputation_pooling, separate = 9L, within_cohort = 3L, across_cohorts = 1L)
      testthat::expect_identical(length(fits_delta[[m]]$groups), expected_groups)
      expected_p <- if (fx$spec$or_slopes == "") 0L else 2L
      testthat::expect_true(all(vapply(fits_delta[[m]]$groups, `[[`, integer(1L), "p") == expected_p))
      testthat::expect_identical(fits_delta[[m]]$accounting$n_imputation_fits, expected_groups)
      testthat::expect_identical(fits_delta[[m]]$accounting$observation_logit_models, 0L)
      testthat::expect_identical(fits_delta[[m]]$accounting$n_or_fits, if (m == "ipw") 0L else 1L)
      testthat::expect_identical(fits_delta[[m]]$accounting$n_treatment_fits, if (m == "reg") 0L else 4L)
      fits_y[[m]] <- imp_fit(fx, est_method = m, pair_outcomes = NULL)
      imp_ns("validate_imputation_backend")(fits_y[[m]])

      # Also fit pure-R table helper on Y_obs and verify exact machine-precision parity with sealed owner
      tb_y <- imp_loader$fit_imputation_table(
        X = fx$X_full,
        Y = fx$Y_obs,
        obs = fx$R,
        g = fx$raw$cohort,
        w = fx$raw$w,
        calendar = fx$catalog$calendar$time,
        links = fx$catalog$link_list,
        outcome_recipe = fx$outcome_recipe,
        treatment_columns = fx$treatment_columns,
        est_method = m,
        imputation_pooling = fx$imputation_pooling,
        glm_control = fx$glm_control,
        rank_tol = fx$rank_tol
      )
      pts_owner_y <- vapply(seq_along(fx$views), function(j) imp_ns("chained_full_link_estimate")(fits_y[[m]], j), numeric(1))
      testthat::expect_lt(imp_mx(pts_owner_y - tb_y$estimates), 2e-14)
      col_y <- imp_collect(fits_y[[m]], chunk_rows = 9L, link_batch_size = 2L, scale = "frame", details = FALSE)
      testthat::expect_lt(imp_mx(col_y$phi - tb_y$phi), 2e-14)

      # Real owner on 8-decimal Y_obs vs 10-decimal Delta_e reference
      testthat::expect_lt(imp_mx(pts_owner_y - exp$point_method(m, "estimate")), 5e-8)
      testthat::expect_lt(imp_mx(col_y$phi - exp$uif_method_mat(m, "analytic")), 1.5e-7)

      # Exact match on Delta_e for all 108 method-links
      fit_m <- fits_delta[[m]]
      pts_m <- vapply(seq_along(fx$views), function(j) imp_ns("chained_full_link_estimate")(fit_m, j), numeric(1))
      testthat::expect_lt(imp_mx(pts_m - exp$point_method(m, "estimate")), 1e-14)
      testthat::expect_lt(imp_mx(vapply(fit_m$links, `[[`, numeric(1), "T_gl") - exp$point_method(m, "T")), 1e-14)
      testthat::expect_lt(imp_mx(vapply(fit_m$links, `[[`, numeric(1), "d_plus") - exp$point_method(m, "d_plus")), 1e-14)
      if (m %in% c("ipw", "dr")) {
        testthat::expect_lt(imp_mx(vapply(fit_m$links, `[[`, numeric(1), "C_ol") - exp$point_method(m, "C")), 1e-14)
        testthat::expect_lt(imp_mx(vapply(fit_m$links, `[[`, numeric(1), "d_minus") - exp$point_method(m, "d_minus")), 1e-14)
      }
      if (m %in% c("reg", "dr")) {
        testthat::expect_lt(imp_mx(vapply(fit_m$links, `[[`, numeric(1), "B_gl") - exp$point_method(m, "Bg")), 1e-14)
      }
      if (m == "dr") {
        testthat::expect_lt(imp_mx(vapply(fit_m$links, `[[`, numeric(1), "B_ol") - exp$point_method(m, "Bo")), 1e-14)
      }
      total_method_links <- total_method_links + length(pts_m)

      # Replay IF on frame and core scales + private producer dispatch
      col_frame <- imp_collect(fit_m, chunk_rows = 9L, link_batch_size = 2L, scale = "frame", details = FALSE)
      col_core <- imp_collect(fit_m, chunk_rows = 11L, link_batch_size = 3L, scale = "core", details = FALSE)
      testthat::expect_identical(col_frame$info$normalization_count, 36L)
      testthat::expect_identical(col_core$info$normalization_count, 24L)
      testthat::expect_lt(imp_mx(col_core$phi - (n_A / N) * col_frame$phi), 1e-15)
      testthat::expect_lt(imp_mx(col_frame$phi - exp$uif_method_mat(m, "analytic")), 3e-14)
      testthat::expect_lt(imp_mx(col_frame$phi - exp$uif_method_mat(m, "fd_central_1e6")), 1e-8)

      # Private producer dispatch (chained_full_if_budget, chained_full_if_chunk, cp04_plan_rows)
      planned_B <- imp_ns("cp04_plan_rows")(fit_m, list(1:2, 3:4))
      testthat::expect_true(is.numeric(planned_B) && length(planned_B) == 1L &&
        is.finite(planned_B) && planned_B == 11 && planned_B == floor(planned_B))
      b_info <- imp_ns("chained_full_if_budget")(fit_m, 7L, 1:4)
      testthat::expect_true(b_info$output > 0 && b_info$scratch > 0 && b_info$work > 0)
      ch_chunk <- imp_ns("chained_full_if_chunk")(fit_m, 1:7, 1:4, "frame")
      testthat::expect_lt(imp_mx(ch_chunk$phi - col_frame$phi[1:7, , drop = FALSE]), 1e-15)

      # Streamed covariance vs reference covariance.csv
      cov_stream <- imp_ns("stream_imputation_covariance")(fit_m)
      V_unit_stream <- cov_stream$covariance / n_A
      cl_covs <- imp_cluster_cov(col_frame$phi, fx$raw$cluster_id, n_A)
      testthat::expect_lt(imp_mx(V_unit_stream - cl_covs$unit_coefficient), 1e-15)
      for (ctype in names(cl_covs)) {
        ref_c <- exp$cov_matrix(ctype, m)
        scale_div <- if (ctype == "cluster_sqrtN") N else 1.0
        testthat::expect_lt(imp_mx(cl_covs[[ctype]] - ref_c) / scale_div, 1e-14)
        # Complete joint covariance is counted once per fixture below.
      }

      # Fixed-loading projections (projection.csv & projection_if.csv)
      for (p_idx in 0:1) {
        p_row <- exp$projection_df[exp$projection_df$method == m & exp$projection_df$coordinate == p_idx, , drop = FALSE]
        testthat::expect_equal(nrow(p_row), 1L)
        w_vec <- as.double(strsplit(p_row$fixed_loading_row, ";", fixed = TRUE)[[1L]])
        testthat::expect_length(w_vec, 4L)
        est_proj <- sum(w_vec * pts_m)
        phi_proj <- as.vector(col_frame$phi %*% w_vec)
        var_unit_proj <- as.double(t(w_vec) %*% cl_covs$unit_coefficient %*% w_vec)
        var_cl_proj <- as.double(t(w_vec) %*% cl_covs$cluster_coefficient %*% w_vec)
        var_sqrtN_proj <- as.double(t(w_vec) %*% cl_covs$cluster_sqrtN %*% w_vec)
        var_core_proj <- as.double(t(w_vec) %*% cl_covs$core_cluster_coefficient %*% w_vec)

        testthat::expect_lt(abs(est_proj - as.double(p_row$value)), 1e-14)
        for (ctype in names(cl_covs)) {
          got <- switch(ctype, unit_coefficient = var_unit_proj,
            cluster_coefficient = var_cl_proj, cluster_sqrtN = var_sqrtN_proj,
            core_cluster_coefficient = var_core_proj)
          ref_var <- as.double(t(w_vec) %*% exp$cov_matrix(ctype, m) %*% w_vec)
          testthat::expect_lt(abs(got - ref_var) / if (ctype == "cluster_sqrtN") N else 1, 1e-14)
        }
        total_proj <- total_proj + 1L

        pif_sub <- exp$projection_if_df[exp$projection_if_df$method == m & exp$projection_if_df$coordinate == p_idx, , drop = FALSE]
        idx_u <- match(fx$raw$unit_id, pif_sub$unit_id)
        testthat::expect_lt(imp_mx(phi_proj - as.double(pif_sub$psi[idx_u])), 3e-14)
        total_proj_if <- total_proj_if + length(phi_proj)
      }
    }

    # 2. Full nuisance, unit_if, 8-channel components, predictions, and sparse adjoint verification on DR fit
    ref_sum <- imp_loader$imputation_reference_summary(fits_delta$dr, cluster_id = fx$raw$cluster_id,
      fixed_test_functions = list(one = rep(1, N), x1 = fx$raw$x1, x2 = fx$raw$x2))
    detail_batches <- imp_collect(fits_delta$dr, chunk_rows = 11L, link_batch_size = 4L, details = TRUE)$batches
    for (ctype in names(ref_sum$covariances)) {
      ref_joint <- exp$cov_matrix(ctype)
      got_joint <- ref_sum$covariances[[ctype]][rownames(ref_joint), colnames(ref_joint), drop = FALSE]
      testthat::expect_lt(imp_mx(got_joint - ref_joint) / if (ctype == "cluster_sqrtN") N else 1, 1e-14)
      total_cov <- total_cov + length(ref_joint)
    }

    # Nuisance scalars
    nuis_names <- names(exp$nuisance)
    got_nuis <- vapply(nuis_names, function(nm) ref_sum$nuisance[[nm]], numeric(1))
    testthat::expect_lt(imp_mx(got_nuis - exp$nuisance), 1e-14)
    total_nuisance <- total_nuisance + length(nuis_names)

    # All unit_if coordinates (including repeat-observer IDs 10, 22, 34 and off-A units)
    uif_names <- unique(exp$unit_if_df$coordinate)
    for (nm in uif_names) {
      got_u <- ref_sum$unit_if[, nm]
      ref_u <- exp$uif_vec(nm)
      testthat::expect_lt(imp_mx(got_u - ref_u), 3e-14)
      total_uif <- total_uif + length(got_u)
    }

    # 8-channel components and sparse adjoint identity
    for (j in 1:4) {
      lid <- link_ids[j]
      ld <- lapply(detail_batches[[1L]]$details$links[[j]], identity)
      for (field in names(ld)) {
        if (field == "link_id" || is.null(ld[[field]])) next
        vals <- lapply(detail_batches, function(ch) ch$details$links[[j]][[field]])
        ld[[field]] <- if (is.matrix(vals[[1L]])) do.call(rbind, vals) else unlist(vals, use.names = FALSE)
      }
      testthat::expect_lt(imp_mx(ld$j_T_adjoint - ld$j_T_forward), 1e-14)
      testthat::expect_lt(imp_mx(ld$j_C_adjoint - ld$j_C_forward), 1e-14)
      for (m in c("ipw", "reg", "dr")) {
        comp_got <- ld[[paste0("components_", m)]]
        comp_ref <- exp$comp_method_link(m, lid)
        testthat::expect_lt(imp_mx(comp_got - comp_ref), 3e-14)
        testthat::expect_lt(imp_mx(rowSums(comp_got) - ld[[m]]), 1e-15)
        total_comp <- total_comp + nrow(comp_got)
      }
    }

    # Cell predictions on all N = 36 F rows
    for (cid in seq_along(fits_delta$dr$cells)) {
      cf <- fits_delta$dr$cells[[cid]]
      gf <- fits_delta$dr$groups[[cf$group_id]]
      Z_imp <- if (gf$p > 0L) fx$X_full[, gf$kept_cols_F, drop = FALSE] else matrix(0.0, N, 0L)
      jhat <- if (gf$p > 0L) cf$alpha + as.vector((cf$L * Z_imp) %*% gf$beta) else rep(cf$alpha, N)
      psub <- exp$predictions_df[exp$predictions_df$cell == cf$name, , drop = FALSE]
      idx_p <- match(fx$raw$unit_id, psub$unit_id)
      testthat::expect_lt(imp_mx(jhat - as.double(psub$j_hat[idx_p])), 1e-14)
      total_pred <- total_pred + length(jhat)
    }
  }

  testthat::expect_identical(total_method_links, 108L)
  testthat::expect_identical(total_nuisance, 1315L)
  testthat::expect_identical(total_uif, 47340L)
  testthat::expect_identical(total_comp, 3888L)
  testthat::expect_identical(total_cov, 11664L)
  testthat::expect_identical(total_pred, 2916L)
  testthat::expect_identical(total_proj, 54L)
  testthat::expect_identical(total_proj_if, 1944L)
})

testthat::test_that("intercept-only fixtures F3, F6, F9 satisfy exact IPW == REG == DR parity for points and IFs", {
  for (fid in c("F3_sep_int_only", "F6_within_int_only", "F9_across_int_only")) {
    fx <- imp_load(fid)
    fit_ipw <- imp_fit(fx, est_method = "ipw")
    fit_reg <- imp_fit(fx, est_method = "reg")
    fit_dr <- imp_fit(fx, est_method = "dr")

    pts_ipw <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_ipw, j), numeric(1))
    pts_reg <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_reg, j), numeric(1))
    pts_dr <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_dr, j), numeric(1))
    testthat::expect_lt(imp_mx(pts_ipw - pts_reg), 1e-14)
    testthat::expect_lt(imp_mx(pts_dr - pts_ipw), 1e-14)

    phi_ipw <- imp_collect(fit_ipw, chunk_rows = 12L, details = FALSE)$phi
    phi_reg <- imp_collect(fit_reg, chunk_rows = 12L, details = FALSE)$phi
    phi_dr <- imp_collect(fit_dr, chunk_rows = 12L, details = FALSE)$phi
    testthat::expect_lt(imp_mx(phi_ipw - phi_reg), 1e-14)
    testthat::expect_lt(imp_mx(phi_dr - phi_ipw), 1e-14)
  }
})

testthat::test_that("chunk size, supplied row order permutation, and positive weight scaling leave estimates, IFs, and covariances invariant", {
  fx_std <- imp_load("F4_within_int_x1x2", max_chunk_rows = 36L)
  fit_std <- imp_fit(fx_std, est_method = "dr")
  col_std <- imp_collect(fit_std, chunk_rows = 36L, link_batch_size = 4L, details = FALSE)
  cov_std <- imp_ns("stream_imputation_covariance")(fit_std)$covariance

  # 1. Chunk-size invariance (1, 5, 13, 36)
  for (cr in c(1L, 5L, 13L, 36L)) {
    fx_c <- imp_load("F4_within_int_x1x2", max_chunk_rows = cr)
    fit_c <- imp_fit(fx_c, est_method = "dr")
    pts_c <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_c, j), numeric(1))
    pts_s <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_std, j), numeric(1))
    testthat::expect_lt(imp_mx(pts_c - pts_s), 1e-14)
    col_c <- imp_collect(fit_c, chunk_rows = cr, link_batch_size = 1L, details = FALSE)
    testthat::expect_lt(imp_qr_refit_difference(col_c$phi, col_std$phi), 1e-10)
    cov_c <- imp_ns("stream_imputation_covariance")(fit_c)$covariance
    testthat::expect_lt(imp_mx(cov_c - cov_std), 1e-14)
  }

  # 2. Supplied row-order permutation invariance
  fx_perm <- imp_load("F4_within_int_x1x2", order = imp_perm_order, max_chunk_rows = 7L)
  fit_perm <- imp_fit(fx_perm, est_method = "dr")
  pts_perm <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_perm, j), numeric(1))
  pts_std <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_std, j), numeric(1))
  testthat::expect_lt(imp_mx(pts_perm - pts_std), 1e-14)
  col_perm <- imp_collect(fit_perm, chunk_rows = 7L, link_batch_size = 2L, details = FALSE)
  reidx <- match(fx_std$raw$unit_id, fx_perm$raw$unit_id)
  testthat::expect_lt(imp_qr_refit_difference(col_perm$phi[reidx, , drop = FALSE], col_std$phi), 1e-10)
  cov_perm <- imp_ns("stream_imputation_covariance")(fit_perm)$covariance
  testthat::expect_lt(imp_mx(cov_perm - cov_std), 1e-14)

  # 3. Positive weight-scaling invariance (w -> 3.5 * w)
  tb_1 <- imp_loader$fit_imputation_table(
    X = fx_std$X_full, Y = fx_std$Y_obs, obs = fx_std$R,
    g = fx_std$raw$cohort, w = fx_std$raw$w,
    calendar = fx_std$catalog$calendar$time, links = fx_std$catalog$link_list,
    outcome_recipe = fx_std$outcome_recipe, treatment_columns = fx_std$treatment_columns,
    est_method = "dr", imputation_pooling = "within_cohort"
  )
  tb_scale <- imp_loader$fit_imputation_table(
    X = fx_std$X_full, Y = fx_std$Y_obs, obs = fx_std$R,
    g = fx_std$raw$cohort, w = 3.5 * fx_std$raw$w,
    calendar = fx_std$catalog$calendar$time, links = fx_std$catalog$link_list,
    outcome_recipe = fx_std$outcome_recipe, treatment_columns = fx_std$treatment_columns,
    est_method = "dr", imputation_pooling = "within_cohort"
  )
  testthat::expect_lt(imp_mx(tb_1$estimates - tb_scale$estimates), 1e-14)
  testthat::expect_lt(imp_mx(tb_1$phi - tb_scale$phi), 1e-14)
})

testthat::test_that("repeat-observer IDs 10, 22, 34 with no retained pair S_e = 0 and off-A units contribute nonzero IF entries without entering Stage 1 WLS", {
  fx <- imp_load("F1_sep_int_x1x2")
  fit <- imp_fit(fx, est_method = "dr", pair_outcomes = fx$pair_outcomes_delta)
  col <- imp_collect(fit, chunk_rows = 12L, details = FALSE)
  batches <- imp_collect(fit, chunk_rows = 12L, link_batch_size = 4L, details = TRUE)$batches

  no_pair_rows <- which(fx$raw$unit_id %in% c(10L, 22L, 34L))
  off_a_rows <- which(fx$A == 0L)
  testthat::expect_identical(length(no_pair_rows), 3L)
  testthat::expect_identical(length(off_a_rows), 12L)

  # Stage 1 score is strictly zero on units 10, 22, 34 and all off-A units
  scores <- do.call(rbind, lapply(batches, function(ch) ch$details$imputation_cell_scores))
  testthat::expect_equal(imp_mx(scores[c(no_pair_rows, off_a_rows), , drop = FALSE]), 0)
  # Group slope scores use these same observed-pair residuals, hence also vanish.


  # Complete F-row IF is nonzero on units 10, 22, 34 and off-A units
  testthat::expect_gt(imp_mx(col$phi[no_pair_rows, , drop = FALSE]), 0.05)
  testthat::expect_gt(imp_mx(col$phi[off_a_rows, , drop = FALSE]), 0.05)
})

testthat::test_that("method dependency isolation holds: ipw skips OR, reg skips treatment logit, and none fits observation logits or forces selection_xformla", {
  fx <- imp_load("F1_sep_int_x1x2")

  # 1. Unused selection_xformla is never forced; uses_selection_models is FALSE
  fit_ipw <- imp_fit(fx, est_method = "ipw", selection_xformla = stop("selection_xformla must not be forced"),
    selection_model = stop("selection_model must not be forced"),
    selection_pooling = stop("selection_pooling must not be forced"))
  fit_reg <- imp_fit(fx, est_method = "reg",
    treatment_columns = stop("unused treatment recipe must not be forced"),
    glm_control = stop("unused treatment fit control must not be forced"), selection_xformla = stop("selection_xformla must not be forced"),
    selection_model = stop("selection_model must not be forced"),
    selection_pooling = stop("selection_pooling must not be forced"))
  fit_dr <- imp_fit(fx, est_method = "dr", selection_xformla = stop("selection_xformla must not be forced"),
    selection_model = stop("selection_model must not be forced"),
    selection_pooling = stop("selection_pooling must not be forced"))

  testthat::expect_false(fit_ipw$uses_selection_models)
  testthat::expect_identical(fit_ipw$accounting$n_selection_fits, 0L)
  testthat::expect_null(fit_ipw$shared_or)
  testthat::expect_identical(length(fit_ipw$or_blocks), 0L)
  testthat::expect_identical(length(fit_ipw$treatments), 4L)

  testthat::expect_false(fit_reg$uses_selection_models)
  testthat::expect_identical(fit_reg$accounting$n_selection_fits, 0L)
  testthat::expect_false(is.null(fit_reg$shared_or))
  testthat::expect_identical(length(fit_reg$or_blocks), 3L)
  testthat::expect_identical(length(fit_reg$treatments), 0L)

  testthat::expect_false(fit_dr$uses_selection_models)
  testthat::expect_identical(fit_dr$accounting$n_selection_fits, 0L)
  testthat::expect_false(is.null(fit_dr$shared_or))
  testthat::expect_identical(length(fit_dr$or_blocks), 3L)
  testthat::expect_identical(length(fit_dr$treatments), 4L)

  # 2. When treatment covariates cause complete separation in treatment propensity logit,
  # est_method = "reg" still succeeds while "ipw" and "dr" refuse with cdid_cp04_separation_error
  X_sep <- cbind(fx$X_full, x_sep = as.double(fx$raw$cohort == 2005))
  testthat::expect_error(
    imp_loader$fit_imputation_table(
      X = X_sep, Y = fx$Y_obs, obs = fx$R, g = fx$raw$cohort, w = fx$raw$w,
      calendar = fx$catalog$calendar$time, links = fx$catalog$link_list,
      outcome_recipe = fx$outcome_recipe,
      treatment_columns = c("(Intercept)", "x1", "x2", "x_sep"),
      est_method = "ipw", imputation_pooling = "separate"
    ),
    class = "cdid_cp04_separation_error"
  )
  fit_reg_sep <- imp_loader$fit_imputation_table(
    X = X_sep, Y = fx$Y_obs, obs = fx$R, g = fx$raw$cohort, w = fx$raw$w,
    calendar = fx$catalog$calendar$time, links = fx$catalog$link_list,
    outcome_recipe = fx$outcome_recipe,
    treatment_columns = c("(Intercept)", "x1", "x2", "x_sep"),
    est_method = "reg", imputation_pooling = "separate"
  )
  pts_reg_base <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_reg, j), numeric(1))
  testthat::expect_lt(imp_mx(fit_reg_sep$estimates - pts_reg_base), 2e-14)
})

testthat::test_that("design resolution refuses mismatched probability Z on reader unless sealed design_owner or X_frame is supplied", {
  fx <- imp_load("F1_sep_int_x1x2")

  # Reader in fx has 1-column Q_prob while frame$Z_a has 3 columns -> omitting design_owner & X_frame must error
  testthat::expect_error(
    imp_fit(fx, est_method = "dr", design_owner = NULL, X_frame = NULL),
    class = "cdid_imputation_design_error"
  )

  # Supplying X_frame directly succeeds and matches design_owner fit identically
  fit_do <- imp_fit(fx, est_method = "dr", design_owner = fx$design_owner, X_frame = NULL)
  fit_xf <- imp_fit(fx, est_method = "dr", design_owner = NULL, X_frame = fx$X_full)
  pts_do <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_do, j), numeric(1))
  pts_xf <- vapply(1:4, function(j) imp_ns("chained_full_link_estimate")(fit_xf, j), numeric(1))
  testthat::expect_lt(imp_mx(pts_do - pts_xf), 1e-15)

  # Perturbing X_frame on an A row so it disagrees with frame$Z_a is refused
  X_bad <- fx$X_full
  a_row <- which(fx$A == 1L)[1L]
  X_bad[a_row, "x1"] <- X_bad[a_row, "x1"] + 1.0
  testthat::expect_error(
    imp_ns("construct_imputation_design_owner")(frame_owner = fx$frame, reader = fx$reader, id = fx$reader$state$id, X = X_bad, control = fx$control),
    class = "cdid_imputation_identity_error"
  )
})

testthat::test_that("group-local retained basis p_{imp,r} preserves predictions on valid cohort aliases and refuses off-S_e alias breakdowns on F", {
  fx <- imp_load("F1_sep_int_x1x2")

  # 1. Exact duplicate column x1_dup = x1 and linear combination x_comb = 2 * x1 - 0.5 * x2 on all F rows
  X_aug <- cbind(
    fx$X_full,
    x1_dup = fx$X_full[, "x1"],
    x_comb = 2.0 * fx$X_full[, "x1"] - 0.5 * fx$X_full[, "x2"]
  )
  rec_aug <- imp_ns("or_pooled_outcome_recipe")(
    c("(Intercept)", "x1", "x2", "x1_dup", "x_comb"),
    TRUE,
    c("x1", "x2", "x1_dup", "x_comb")
  )
  tb_base <- imp_loader$fit_imputation_table(
    X = fx$X_full, Y = fx$Y_obs, obs = fx$R, g = fx$raw$cohort, w = fx$raw$w,
    calendar = fx$catalog$calendar$time, links = fx$catalog$link_list,
    outcome_recipe = fx$outcome_recipe, treatment_columns = fx$treatment_columns,
    est_method = "dr", imputation_pooling = "separate"
  )
  tb_aug <- imp_loader$fit_imputation_table(
    X = X_aug, Y = fx$Y_obs, obs = fx$R, g = fx$raw$cohort, w = fx$raw$w,
    calendar = fx$catalog$calendar$time, links = fx$catalog$link_list,
    outcome_recipe = rec_aug, treatment_columns = fx$treatment_columns,
    est_method = "dr", imputation_pooling = "separate"
  )
  testthat::expect_lt(imp_mx(tb_aug$estimates - tb_base$estimates), 1e-13)
  testthat::expect_lt(imp_mx(tb_aug$phi - tb_base$phi), 1e-13)
  testthat::expect_true(all(vapply(tb_aug$core$groups, `[[`, integer(1), "p") == 2L))

  # 2. Group-local rank heterogeneity: in cohort 0, set x2 = x1 on ALL cohort-0 rows of F,
  # while cohorts 2005 and 2010 keep linearly independent (x1, x2).
  X_het <- fx$X_full
  c0_rows <- which(fx$raw$cohort == 0)
  X_het[c0_rows, "x2"] <- X_het[c0_rows, "x1"]
  tb_het <- imp_loader$fit_imputation_table(
    X = X_het, Y = fx$Y_obs, obs = fx$R, g = fx$raw$cohort, w = fx$raw$w,
    calendar = fx$catalog$calendar$time, links = fx$catalog$link_list,
    outcome_recipe = fx$outcome_recipe, treatment_columns = fx$treatment_columns,
    est_method = "dr", imputation_pooling = "separate"
  )
  p_by_group <- vapply(tb_het$core$groups, `[[`, integer(1), "p")
  # Cohort 0 cells have p_{imp,r} == 1 while cohorts 2005 and 2010 have p_{imp,r} == 2
  testthat::expect_identical(unname(p_by_group[1:3]), c(1L, 1L, 1L))
  testthat::expect_identical(unname(p_by_group[4:9]), rep(2L, 6L))

  # 3. Alias that holds on observed pairs S_{ie} = 1 of cohort 0, e12, but breaks on an
  # unobserved F row (S_{i,e12} = 0) of cohort 0 -> MUST refuse with cdid_or_wls_rank_error
  # and carry is_identifiability_claim == FALSE
  X_break <- X_het
  c0_unobs <- which(fx$raw$cohort == 0 & fx$raw$S_e12 == 0L)[1L]
  X_break[c0_unobs, "x2"] <- X_break[c0_unobs, "x1"] + 0.5
  err <- tryCatch(
    imp_loader$fit_imputation_table(
      X = X_break, Y = fx$Y_obs, obs = fx$R, g = fx$raw$cohort, w = fx$raw$w,
      calendar = fx$catalog$calendar$time, links = fx$catalog$link_list,
      outcome_recipe = fx$outcome_recipe, treatment_columns = fx$treatment_columns,
      est_method = "dr", imputation_pooling = "separate"
    ),
    cdid_or_wls_rank_error = identity
  )
  testthat::expect_s3_class(err, "cdid_or_wls_rank_error")
  testthat::expect_false(isTRUE(err$is_identifiability_claim))
  testthat::expect_match(conditionMessage(err), "does not establish population ATT nonidentification")
})

testthat::test_that("resource bounds are enforced before allocations with typed cdid_ipw_owner_bound_error", {
  fx <- imp_load("F1_sep_int_x1x2", max_chunk_rows = 11L)
  fit <- imp_fit(fx, est_method = "dr")

  # Replay rows exceeding max_chunk_rows
  testthat::expect_error(
    imp_ns("evaluate_imputation_if_chunk")(fit, 1:12, 1:4),
    class = "cdid_ipw_owner_bound_error"
  )

  # Tight max_retained_entries refuses fit preflight
  ctrl_tight_ret <- imp_ns("ipw_owner_control")(
    max_chunk_rows = 11L,
    max_retained_entries = 100L,
    max_matrix_entries = 2e7,
    max_work_ops = 1e10,
    max_output_entries = 1e7
  )
  testthat::expect_error(
    imp_fit(fx, control = ctrl_tight_ret),
    class = "cdid_ipw_owner_bound_error"
  )

  # Tight max_work_ops refuses fit preflight
  ctrl_tight_work <- imp_ns("ipw_owner_control")(
    max_chunk_rows = 11L,
    max_retained_entries = 1e7,
    max_matrix_entries = 2e7,
    max_work_ops = 500L,
    max_output_entries = 1e7
  )
  testthat::expect_error(
    imp_fit(fx, control = ctrl_tight_work),
    class = "cdid_ipw_owner_bound_error"
  )
})
