local({
  ref_env <- new.env(parent = baseenv())
  sys.source(
    testthat::test_path("reference", "static_observation_adapter_reference.R"),
    envir = ref_env
  )

  # Deterministic fixture with N = 60 original units in non-sorted character ID
  # order, including zero-outcome units, single-outcome units, and multi-outcome
  # units across overlapping fixed role groups (G and C) and varying designs.
  make_adapter_fixture <- function() {
    set.seed(20261104L)
    N <- 60L
    # Arbitrary non-sequential original unit IDs
    unit_id <- sprintf("id_%03d", sample.int(300L, N, replace = FALSE))

    # Fixed role masks:
    # - g1_G (switcher role G) and g1_C (stayer role C) overlap on several units
    # - g2_C has a third overlapping population
    # - last 6 units are outside all masks; first 6 role units have 0 outcomes;
    #   next 6 role units have only 1 outcome (date t1 only -> pair O = 0).
    A_g1_G <- logical(N); A_g1_G[c(1:28, 35:40)] <- TRUE
    A_g1_C <- logical(N); A_g1_C[c(15:48)] <- TRUE
    A_g2_C <- logical(N); A_g2_C[c(10:32, 42:54)] <- TRUE

    # Date observation indicators on t1, t2, t3
    O_t1 <- rbinom(N, 1L, 0.78) == 1L
    O_t2 <- rbinom(N, 1L, 0.74) == 1L
    O_t3 <- rbinom(N, 1L, 0.72) == 1L

    # Enforce explicit 0-outcome units (1:5) and single-outcome units (6:10)
    O_t1[1:5] <- FALSE; O_t2[1:5] <- FALSE; O_t3[1:5] <- FALSE
    O_t1[6:10] <- TRUE; O_t2[6:10] <- FALSE; O_t3[6:10] <- FALSE

    # Ensure balanced interior support across cells
    O_t1[11:16] <- TRUE; O_t2[11:16] <- TRUE; O_t3[11:16] <- TRUE
    O_t1[17:20] <- TRUE; O_t2[17:20] <- FALSE; O_t3[17:20] <- TRUE
    O_t1[21:24] <- FALSE; O_t2[21:24] <- TRUE; O_t3[21:24] <- FALSE

    O_p12 <- O_t1 & O_t2
    O_p13 <- O_t1 & O_t3

    make_X <- function(seed_offset) {
      set.seed(20261104L + seed_offset)
      X <- cbind(
        x1 = rnorm(N, 0, 0.55),
        x2 = runif(N, -0.7, 0.7)
      )
      rownames(X) <- unit_id
      X
    }

    X_sel_base <- make_X(1L)
    X_sel_alt  <- make_X(2L)
    X_imp_base <- make_X(3L)
    X_imp_alt  <- make_X(4L)
    Z_stayer   <- cbind(z1 = X_sel_base[, 1L] + 0.2, z2 = X_sel_base[, 2L] - 0.1)
    rownames(Z_stayer) <- unit_id

    dy_12_G <- 0.45 + 0.6 * X_imp_base[, 1L] - 0.35 * X_imp_base[, 2L] + rnorm(N, 0, 0.25)
    dy_13_G <- 0.80 + 1.1 * X_imp_alt[, 1L]  - 0.50 * X_imp_alt[, 2L]  + rnorm(N, 0, 0.30)
    dy_12_C <- -0.15 + 0.4 * X_imp_base[, 1L] + 0.25 * X_imp_base[, 2L] + rnorm(N, 0, 0.25)
    dy_13_C <- -0.25 + 0.7 * X_imp_alt[, 1L]  + 0.45 * X_imp_alt[, 2L]  + rnorm(N, 0, 0.30)
    dy_12_g2 <- 0.10 - 0.3 * X_imp_base[, 1L] + 0.55 * X_imp_base[, 2L] + rnorm(N, 0, 0.25)

    names(dy_12_G) <- names(dy_13_G) <- names(dy_12_C) <- names(dy_13_C) <- names(dy_12_g2) <- unit_id

    dates_g1_G <- list(
      t1 = list(observed = O_t1, selection_design = X_sel_base),
      t2 = list(observed = O_t2, selection_design = X_sel_alt),
      t3 = list(observed = O_t3, selection_design = X_sel_base)
    )
    dates_g1_C <- list(
      t1 = list(observed = O_t1, selection_design = X_sel_alt),
      t2 = list(observed = O_t2, selection_design = X_sel_base),
      t3 = list(observed = O_t3, selection_design = X_sel_alt)
    )
    dates_g2_C <- list(
      t1 = list(observed = O_t1, selection_design = X_sel_base),
      t2 = list(observed = O_t2, selection_design = X_sel_base)
    )

    groups <- list(
      g1_G = list(
        role = "G",
        target = "2014",
        prior_state = 0L,
        mask = A_g1_G,
        Z = Z_stayer,
        dates = dates_g1_G,
        pairs = list(
          p12 = list(
            start = "t1", end = "t2", gap = 1L,
            observed = O_p12, dy = dy_12_G,
            selection_design = X_sel_base,
            imputation_design = X_imp_base
          ),
          p13 = list(
            start = "t1", end = "t3", gap = 2L,
            observed = O_p13, dy = dy_13_G,
            selection_design = X_sel_alt,
            imputation_design = X_imp_alt
          )
        )
      ),
      g1_C = list(
        role = "C",
        target = "2014",
        prior_state = 0L,
        mask = A_g1_C,
        Z = Z_stayer,
        dates = dates_g1_C,
        pairs = list(
          p12 = list(
            start = "t1", end = "t2", gap = 1L,
            observed = O_p12, dy = dy_12_C,
            selection_design = X_sel_base,
            imputation_design = X_imp_base
          ),
          p13 = list(
            start = "t1", end = "t3", gap = 2L,
            observed = O_p13, dy = dy_13_C,
            selection_design = X_sel_alt,
            imputation_design = X_imp_alt
          )
        )
      ),
      g2_C = list(
        role = "C",
        target = "2016",
        prior_state = 1L,
        mask = A_g2_C,
        Z = Z_stayer,
        dates = dates_g2_C,
        pairs = list(
          p12 = list(
            start = "t1", end = "t2", gap = 1L,
            observed = O_p12, dy = dy_12_g2,
            selection_design = X_sel_base,
            imputation_design = X_imp_base
          )
        )
      )
    )

    # Upstream fitted treatment tilt h on union of g1_G and g1_C for testing
    # non-trivial h + h_jacobian + h_coef_if provider.
    A_trt <- A_g1_G | A_g1_C
    y_trt <- as.double(A_g1_G)
    X_trt <- cbind(1, X_sel_base)
    colnames(X_trt) <- c("alpha", "x1", "x2")
    fit_trt_fn <- function(w_u = rep.int(1, N)) {
      rows <- which(A_trt)
      f <- suppressWarnings(stats::glm.fit(
        X_trt[rows, , drop = FALSE], y_trt[rows], weights = w_u[rows],
        family = stats::binomial(),
        control = stats::glm.control(epsilon = 1e-14, maxit = 200)
      ))
      unname(f$coefficients)
    }
    gamma_hat <- fit_trt_fn()
    eta_trt <- drop(X_trt %*% gamma_hat)
    p_trt <- stats::plogis(eta_trt)
    h_vec <- stats::setNames(exp(eta_trt), unit_id)
    h_jac_mat <- matrix(0, nrow = N, ncol = 3L,
                        dimnames = list(unit_id, colnames(X_trt)))
    h_jac_mat[A_trt, ] <- h_vec[A_trt] * X_trt[A_trt, , drop = FALSE]
    H_trt <- crossprod(
      X_trt[A_trt, , drop = FALSE],
      (p_trt[A_trt] * (1 - p_trt[A_trt])) * X_trt[A_trt, , drop = FALSE]
    )
    S_trt <- matrix(0, nrow = N, ncol = 3L)
    S_trt[A_trt, ] <- (y_trt[A_trt] - p_trt[A_trt]) * X_trt[A_trt, , drop = FALSE]
    L_trt <- N * t(solve(H_trt, t(S_trt)))
    dimnames(L_trt) <- list(unit_id, colnames(X_trt))

    h_fitter <- function(w_u) {
      g_w <- fit_trt_fn(w_u)
      exp(drop(X_trt %*% g_w))
    }

    list(
      N = N,
      unit_id = unit_id,
      groups = groups,
      Z = Z_stayer,
      h_one = stats::setNames(rep.int(1, N), unit_id),
      h_vec = h_vec,
      h_jacobian = list(trt_block = h_jac_mat),
      h_coef_if_list = list(trt_block = L_trt),
      h_coef_if_fn = function(ids) {
        pos <- match(ids, unit_id)
        list(trt_block = L_trt[pos, , drop = FALSE])
      },
      h_fitter = h_fitter,
      glm_control = list(maxit = 100L, epsilon = 1e-12, trace = FALSE),
      rank_tol = 1e-10
    )
  }

  testthat::test_that("all 4 routes, pair & date_product models, and 3 poolings match stacked glm/lm and literal ratio/pseudo-change reference", {
    fx <- make_adapter_fixture()

    for (route in c("baseline", "weighting", "imputation", "augmented")) {
      sel_models <- if (route %in% c("weighting", "augmented")) {
        c("pair_logit", "independent_period_logit")
      } else {
        list(NULL)
      }
      poolings <- c("separate", "within_group", "pooled")
      for (s_mod in sel_models) {
        for (pool in poolings) {
          ad <- static_observation_adapter_fit(
            unit_id = fx$unit_id,
            groups = fx$groups,
            route = route,
            selection_model = s_mod,
            selection_pooling = if (!is.null(s_mod)) pool else NULL,
            imputation_pooling = if (route %in% c("imputation", "augmented")) pool else NULL,
            imputation_scale = if (route %in% c("imputation", "augmented")) "per_period" else NULL,
            glm_control = if (!is.null(s_mod)) fx$glm_control else NULL,
            rank_tol = if (route %in% c("imputation", "augmented")) fx$rank_tol else NULL,
            max_chunk_rows = 17L
          )

          # Check that only needed nuisance fits are constructed
          if (route %in% c("baseline", "imputation")) {
            testthat::expect_null(ad$selection_fit)
          } else {
            testthat::expect_s3_class(ad$selection_fit, "cdid_static_nuisance_fit")
          }
          if (route %in% c("baseline", "weighting")) {
            testthat::expect_null(ad$imputation_fit)
          } else {
            testthat::expect_s3_class(ad$imputation_fit, "cdid_static_nuisance_fit")
          }

          # Compare fitted nuisance coefficients against independent stacked glm.fit / lm.wfit
          ref_sp <- ref_env$adapter_ref_build_specs(
            N = fx$N, groups = fx$groups, route = route,
            selection_model = s_mod, selection_pooling = pool,
            imputation_pooling = pool, imputation_scale = "per_period"
          )
          if (!is.null(ad$selection_fit)) {
            ref_sel <- ref_env$adapter_ref_fit_all(
              ref_sp$sel_cells, ref_sp$sel_partition, "logit", rep.int(1, fx$N)
            )
            testthat::expect_equal(
              unname(ad$selection_fit$coefficients),
              ref_sel$coefficients,
              tolerance = 1e-10
            )
          }
          if (!is.null(ad$imputation_fit)) {
            ref_imp <- ref_env$adapter_ref_fit_all(
              ref_sp$imp_cells, ref_sp$imp_partition, "wls", rep.int(1, fx$N)
            )
            testthat::expect_equal(
              unname(ad$imputation_fit$coefficients),
              ref_imp$coefficients,
              tolerance = 1e-12
            )
          }

          # Compare mean (with estimated tilt h) and cov (with Z) against independent literal reference
          for (edge in list(c("g1_G", "p12"), c("g1_G", "p13"), c("g1_C", "p13"))) {
            g_k <- edge[1L]
            p_k <- edge[2L]
            m_obj <- ad$mean(
              group = g_k, pair = p_k,
              h = fx$h_vec, h_jacobian = fx$h_jacobian, h_coef_if = fx$h_coef_if_fn
            )
            c_obj <- ad$cov(group = g_k, pair = p_k, Z = fx$Z)

            ref_m <- ref_env$adapter_ref_edge_functional(
              N = fx$N, groups = fx$groups, group_key = g_k, pair_key = p_k,
              route = route, selection_model = s_mod, selection_pooling = pool,
              imputation_pooling = pool, imputation_scale = "per_period",
              kind = "mean", h_vec = fx$h_vec
            )
            ref_c <- ref_env$adapter_ref_edge_functional(
              N = fx$N, groups = fx$groups, group_key = g_k, pair_key = p_k,
              route = route, selection_model = s_mod, selection_pooling = pool,
              imputation_pooling = pool, imputation_scale = "per_period",
              kind = "cov", Z = fx$Z
            )

            testthat::expect_equal(unname(m_obj$estimate), unname(ref_m), tolerance = 1e-11)
            testthat::expect_equal(unname(c_obj$estimate), unname(ref_c), tolerance = 1e-11)

            # Mean of full-N IF is ~0 for centered direct IF, and finite everywhere
            if_m <- m_obj$eval_if()
            if_c <- c_obj$eval_if()
            testthat::expect_identical(names(if_m), fx$unit_id)
            testthat::expect_identical(rownames(if_c), fx$unit_id)
            testthat::expect_true(all(is.finite(if_m)))
            testthat::expect_true(all(is.finite(if_c)))
          }
        }
      }
    }
  })

  testthat::test_that("eval_if matches full-primitive static_obs_mean/cov and is invariant to row chunking and arbitrary ID order", {
    fx <- make_adapter_fixture()

    for (s_mod in c("pair_logit", "independent_period_logit")) {
      for (pool in c("separate", "within_group", "pooled")) {
        ad <- static_observation_adapter_fit(
          unit_id = fx$unit_id,
          groups = fx$groups,
          route = "augmented",
          selection_model = s_mod,
          selection_pooling = pool,
          imputation_pooling = pool,
          imputation_scale = "per_period",
          glm_control = fx$glm_control,
          rank_tol = fx$rank_tol,
          max_chunk_rows = 13L
        )

        m_obj <- ad$mean(
          "g1_G", "p13",
          h = fx$h_vec, h_jacobian = fx$h_jacobian, h_coef_if = fx$h_coef_if_list
        )
        c_obj <- ad$cov("g1_C", "p13", Z = fx$Z)

        # Call accepted full primitives directly on full N for comparison
        A_G <- fx$groups$g1_G$mask
        inA_G <- which(A_G)
        pr_G <- fx$groups$g1_G$pairs$p13
        q_pred_G <- if (s_mod == "pair_logit") {
          static_nuisance_q_predict(ad$selection_fit, fx$unit_id[inA_G], pair = "g1_G::p13")
        } else {
          static_nuisance_q_predict(ad$selection_fit, fx$unit_id[inA_G], dates = c("g1_G::t1", "g1_G::t3"))
        }
        j_pred_G <- static_nuisance_j_predict(ad$imputation_fit, fx$unit_id[inA_G], cell = "g1_G::p13")

        q_full <- rep(NA_real_, fx$N); names(q_full) <- fx$unit_id; q_full[inA_G] <- q_pred_G$value
        j_full <- rep(NA_real_, fx$N); names(j_full) <- fx$unit_id; j_full[inA_G] <- j_pred_G$value
        q_jac <- static_obs_adapter_block_jacobians(
          ad$selection_fit, "selection", inA_G, q_pred_G, fx$N, fx$unit_id
        )
        j_jac <- static_obs_adapter_block_jacobians(
          ad$imputation_fit, "imputation", inA_G, j_pred_G, fx$N, fx$unit_id
        )
        nuis_if_all <- c(fx$h_coef_if_list, ad$nuisance_if(fx$unit_id))
        prim_m <- static_obs_mean(
          id = fx$unit_id, A = A_G, h = fx$h_vec, route = "augmented",
          O = pr_G$observed, dy = pr_G$dy, q = q_full, j = j_full,
          h_jacobian = fx$h_jacobian, q_jacobian = q_jac, j_jacobian = j_jac,
          coef_if = nuis_if_all[names(m_obj$loadings)]
        )

        parts_m <- m_obj$eval_if(include_components = TRUE)
        testthat::expect_equal(parts_m$direct_if, prim_m$direct_if, tolerance = 1e-14)
        testthat::expect_equal(parts_m$influence, prim_m$influence, tolerance = 1e-14)

        # Row-chunk replay and reversed/permuted ID replay must equal full-N IF slice exactly
        perm_ids <- rev(fx$unit_id)
        testthat::expect_equal(m_obj$eval_if(perm_ids), prim_m$influence[perm_ids], tolerance = 1e-14)
        chunk_if <- c(
          m_obj$eval_if(fx$unit_id[1:11]),
          m_obj$eval_if(fx$unit_id[12:37]),
          m_obj$eval_if(fx$unit_id[38:60])
        )
        testthat::expect_equal(chunk_if, prim_m$influence, tolerance = 1e-14)

        # Same check for covariance
        c_full <- c_obj$eval_if()
        c_chunks <- rbind(
          c_obj$eval_if(fx$unit_id[1:19]),
          c_obj$eval_if(fx$unit_id[20:41]),
          c_obj$eval_if(fx$unit_id[42:60])
        )
        testthat::expect_equal(c_chunks, c_full, tolerance = 1e-14)
        testthat::expect_equal(c_obj$eval_if(perm_ids), c_full[perm_ids, , drop = FALSE], tolerance = 1e-14)

        # Closure memory check: moment closure must not retain full-N q_tmp/j_tmp/res objects
        env_m <- environment(m_obj$eval_if)
        testthat::expect_false(exists("res", envir = env_m, inherits = FALSE))
        testthat::expect_false(exists("q_tmp", envir = env_m, inherits = FALSE))
        testthat::expect_false(exists("j_tmp", envir = env_m, inherits = FALSE))
      }
    }
  })

  testthat::test_that("imputation_scale per_period vs per_change: separate equivalence and pooled difference across unequal gaps", {
    fx <- make_adapter_fixture()

    # 1. Under separate pooling, per_period and per_change are exact slope reparameterizations:
    #    predictions, estimates, and complete IFs match; gap-2 slopes scale by 1/2.
    ad_sep_pp <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "augmented",
      selection_model = "pair_logit", selection_pooling = "separate",
      imputation_pooling = "separate", imputation_scale = "per_period",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 25L
    )
    ad_sep_pc <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "augmented",
      selection_model = "pair_logit", selection_pooling = "separate",
      imputation_pooling = "separate", imputation_scale = "per_change",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 25L
    )

    m_sep_pp <- ad_sep_pp$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
    m_sep_pc <- ad_sep_pc$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
    c_sep_pp <- ad_sep_pp$cov("g1_C", "p13")
    c_sep_pc <- ad_sep_pc$cov("g1_C", "p13")

    testthat::expect_equal(m_sep_pp$estimate, m_sep_pc$estimate, tolerance = 1e-12)
    testthat::expect_equal(m_sep_pp$eval_if(), m_sep_pc$eval_if(), tolerance = 1e-12)
    testthat::expect_equal(c_sep_pp$estimate, c_sep_pc$estimate, tolerance = 1e-12)
    testthat::expect_equal(c_sep_pp$eval_if(), c_sep_pc$eval_if(), tolerance = 1e-12)

    # Gap = 2 cell slopes under per_change are 2x the per_period slopes
    b_pp <- ad_sep_pp$imputation_fit$block_fits[["g1_G::p13"]]$beta
    b_pc <- ad_sep_pc$imputation_fit$block_fits[["g1_G::p13"]]$beta
    testthat::expect_equal(2 * b_pp, b_pc, tolerance = 1e-12)

    # 2. Under within_group or pooled pooling across gap=1 (p12) and gap=2 (p13),
    #    per_period and per_change impose genuinely different restrictions.
    for (pool in c("within_group", "pooled")) {
      ad_pool_pp <- static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = fx$groups, route = "augmented",
        selection_model = "pair_logit", selection_pooling = pool,
        imputation_pooling = pool, imputation_scale = "per_period",
        glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 25L
      )
      ad_pool_pc <- static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = fx$groups, route = "augmented",
        selection_model = "pair_logit", selection_pooling = pool,
        imputation_pooling = pool, imputation_scale = "per_change",
        glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 25L
      )
      m_pool_pp <- ad_pool_pp$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
      m_pool_pc <- ad_pool_pc$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
      testthat::expect_gt(abs(m_pool_pp$estimate - m_pool_pc$estimate), 1e-3)
      testthat::expect_gt(max(abs(m_pool_pp$eval_if() - m_pool_pc$eval_if())), 1e-2)
    }
  })

  testthat::test_that("complete original-unit IF matches unit-weight-perturbation refit finite differences and catches mutations", {
    fx <- make_adapter_fixture()
    probe_units <- c(1L, 7L, 14L, 22L, 38L, 57L) # includes 0-outcome, 1-outcome, overlapping, and out-of-role units
    probe_ids <- fx$unit_id[probe_units]

    # Case 1: augmented + independent_period_logit + within_group + estimated h
    ad_date <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "augmented",
      selection_model = "independent_period_logit", selection_pooling = "within_group",
      imputation_pooling = "within_group", imputation_scale = "per_period",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 20L
    )
    m_date <- ad_date$mean(
      "g1_G", "p13",
      h = fx$h_vec, h_jacobian = fx$h_jacobian, h_coef_if = fx$h_coef_if_fn
    )
    analytic_if_m <- m_date$eval_if(probe_ids)

    fd_if_m <- ref_env$adapter_ref_fd_if(
      units = probe_units,
      eval_fn = function(w_u) {
        ref_env$adapter_ref_edge_functional(
          N = fx$N, groups = fx$groups, group_key = "g1_G", pair_key = "p13",
          route = "augmented", selection_model = "independent_period_logit",
          selection_pooling = "within_group", imputation_pooling = "within_group",
          imputation_scale = "per_period", kind = "mean",
          h_fitter = fx$h_fitter, w_unit = w_u
        )
      },
      N = fx$N, eps = 1e-5
    )
    testthat::expect_equal(unname(analytic_if_m), unname(fd_if_m), tolerance = 1e-5)

    # Case 2: weighting + pooled + cov
    ad_wt <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "weighting",
      selection_model = "pair_logit", selection_pooling = "pooled",
      glm_control = fx$glm_control, max_chunk_rows = 20L
    )
    c_wt <- ad_wt$cov("g1_C", "p12", Z = fx$Z)
    analytic_if_c <- c_wt$eval_if(probe_ids)

    fd_if_c <- ref_env$adapter_ref_fd_if(
      units = probe_units,
      eval_fn = function(w_u) {
        ref_env$adapter_ref_edge_functional(
          N = fx$N, groups = fx$groups, group_key = "g1_C", pair_key = "p12",
          route = "weighting", selection_model = "pair_logit",
          selection_pooling = "pooled", kind = "cov", Z = fx$Z, w_unit = w_u
        )
      },
      N = fx$N, eps = 1e-5
    )
    testthat::expect_equal(unname(analytic_if_c), unname(fd_if_c), tolerance = 1e-5)

    # Mutation 1a: losing a shared date parameter contribution in independent_period_logit
    # (zeroing out the start-date alpha or dropping one date's slope contribution)
    # deviates from FD by > 2000x the 1e-5 FD tolerance on augmented and > 0.15 on weighting.
    mut_loadings <- m_date$loadings
    mut_loadings[["selection|g1_G"]]["g1_G|alpha|g1_G::t1"] <- 0
    mut_contrib <- static_obs_contract(
      mut_loadings,
      c(fx$h_coef_if_fn(probe_ids), ad_date$nuisance_if(probe_ids))[names(mut_loadings)],
      probe_ids
    )
    mut_if_date <- m_date$eval_if(probe_ids, include_components = TRUE)$direct_if + drop(mut_contrib)
    testthat::expect_gt(max(abs(unname(mut_if_date) - unname(fd_if_m))), 0.015)

    # Also check weighting + independent_period_logit + separate: p12 and p13 share date cell g1_G::t1
    ad_wt_date <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "weighting",
      selection_model = "independent_period_logit", selection_pooling = "separate",
      glm_control = fx$glm_control, max_chunk_rows = 20L
    )
    m_wt_12 <- ad_wt_date$mean("g1_G", "p12", h = fx$h_one, h_jacobian = list())
    m_wt_13 <- ad_wt_date$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
    testthat::expect_true("selection|g1_G::t1" %in% names(m_wt_12$loadings))
    testthat::expect_true("selection|g1_G::t1" %in% names(m_wt_13$loadings))
    fd_wt_13 <- ref_env$adapter_ref_fd_if(
      units = probe_units,
      eval_fn = function(w_u) {
        ref_env$adapter_ref_edge_functional(
          N = fx$N, groups = fx$groups, group_key = "g1_G", pair_key = "p13",
          route = "weighting", selection_model = "independent_period_logit",
          selection_pooling = "separate", kind = "mean",
          h_vec = fx$h_one, w_unit = w_u
        )
      },
      N = fx$N, eps = 1e-5
    )
    testthat::expect_equal(unname(m_wt_13$eval_if(probe_ids)), unname(fd_wt_13), tolerance = 1e-5)
    # Dropping the shared date cell loading selection|g1_G::t1 fails FD by > 0.10 (>10,000x FD tol)
    mut_wt_load <- m_wt_13$loadings["selection|g1_G::t3"]
    mut_wt_if <- m_wt_13$eval_if(probe_ids, include_components = TRUE)$direct_if +
      drop(static_obs_contract(mut_wt_load, ad_wt_date$nuisance_if(probe_ids)[names(mut_wt_load)], probe_ids))
    testthat::expect_gt(max(abs(unname(mut_wt_if) - unname(fd_wt_13))), 0.10)

    # Mutation 2: wrong N normalization (normalizing by n_role = sum(A_g1_G) instead of full N)
    n_role <- sum(fx$groups$g1_G$mask)
    mut_if_norm <- analytic_if_m * (n_role / fx$N)
    testthat::expect_gt(max(abs(unname(mut_if_norm) - unname(fd_if_m))), 0.25)

    # Mutation 3: subset-normalized direct IF (recomputing S and M on requested subset instead of frozen full-N S, M)
    parts_sub <- m_date$eval_if(probe_ids, include_components = TRUE)
    inA_G <- which(fx$groups$g1_G$mask)
    j_all <- static_nuisance_j_predict(ad_date$imputation_fit, fx$unit_id[inA_G], "g1_G::p13")$value
    q_all <- static_nuisance_q_predict(ad_date$selection_fit, fx$unit_id[inA_G], dates = c("g1_G::t1", "g1_G::t3"))$value
    # Subset of probe_ids inside A_g1_G
    sub_in_A <- intersect(probe_ids, fx$unit_id[inA_G])
    pos_in_A <- match(sub_in_A, fx$unit_id[inA_G])
    pos_in_N <- match(sub_in_A, fx$unit_id)
    O_sub <- fx$groups$g1_G$pairs$p13$observed[pos_in_N]
    phi_sub <- ifelse(O_sub, 1 / ipw_clip_floor(q_all[pos_in_A]), 0)
    resid_sub <- ifelse(O_sub, fx$groups$g1_G$pairs$p13$dy[pos_in_N] - j_all[pos_in_A], 0)
    Y_sub <- j_all[pos_in_A] + phi_sub * resid_sub
    w_sub <- fx$h_vec[pos_in_N]
    S_wrong <- sum(w_sub)
    M_wrong <- sum(w_sub * Y_sub) / S_wrong
    direct_wrong <- length(probe_ids) * w_sub * (Y_sub - M_wrong) / S_wrong
    testthat::expect_gt(max(abs(direct_wrong - parts_sub$direct_if[sub_in_A])), 0.2)
  })

  testthat::test_that("nonfinite entries outside relevant population are ignored, identical-mask groups stay distinct, and pair_callback works", {
    fx <- make_adapter_fixture()
    grp_dirty <- fx$groups

    # Contaminate outside-role rows on masks/designs/h/Z and unobserved rows on dy
    for (g_k in names(grp_dirty)) {
      A <- grp_dirty[[g_k]]$mask
      out_A <- which(!A)
      grp_dirty[[g_k]]$Z[out_A, ] <- NaN
      for (d_k in names(grp_dirty[[g_k]]$dates)) {
        grp_dirty[[g_k]]$dates[[d_k]]$observed[out_A] <- NA
        grp_dirty[[g_k]]$dates[[d_k]]$selection_design[out_A, ] <- Inf
      }
      for (p_k in names(grp_dirty[[g_k]]$pairs)) {
        O <- grp_dirty[[g_k]]$pairs[[p_k]]$observed
        unobs <- which(!(A & O))
        grp_dirty[[g_k]]$pairs[[p_k]]$observed[out_A] <- NA
        grp_dirty[[g_k]]$pairs[[p_k]]$dy[unobs] <- -Inf
        grp_dirty[[g_k]]$pairs[[p_k]]$selection_design[out_A, ] <- NA_real_
        grp_dirty[[g_k]]$pairs[[p_k]]$imputation_design[out_A, ] <- NaN
      }
    }
    h_dirty <- fx$h_vec
    h_dirty[!fx$groups$g1_G$mask] <- NA_real_

    ad_clean <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "augmented",
      selection_model = "independent_period_logit", selection_pooling = "within_group",
      imputation_pooling = "within_group", imputation_scale = "per_period",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 16L
    )
    ad_dirty <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = grp_dirty, route = "augmented",
      selection_model = "independent_period_logit", selection_pooling = "within_group",
      imputation_pooling = "within_group", imputation_scale = "per_period",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 16L
    )

    m_clean <- ad_clean$mean("g1_G", "p13", h = fx$h_vec, h_jacobian = fx$h_jacobian, h_coef_if = fx$h_coef_if_fn)
    m_dirty <- ad_dirty$mean("g1_G", "p13", h = h_dirty, h_jacobian = fx$h_jacobian, h_coef_if = fx$h_coef_if_fn)
    c_clean <- ad_clean$cov("g1_C", "p12")
    c_dirty <- ad_dirty$cov("g1_C", "p12")

    testthat::expect_equal(m_dirty$estimate, m_clean$estimate, tolerance = 1e-15)
    testthat::expect_equal(m_dirty$eval_if(), m_clean$eval_if(), tolerance = 1e-15)
    testthat::expect_equal(c_dirty$estimate, c_clean$estimate, tolerance = 1e-15)
    testthat::expect_equal(c_dirty$eval_if(), c_clean$eval_if(), tolerance = 1e-15)

    # Identical-mask groups must NOT be deduplicated: two groups with identical masks and data
    # still get distinct cell intercepts and keys in pooled fit
    grp_dup <- list(
      gA = fx$groups$g1_G,
      gB = fx$groups$g1_G
    )
    ad_dup <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = grp_dup, route = "augmented",
      selection_model = "pair_logit", selection_pooling = "pooled",
      imputation_pooling = "pooled", imputation_scale = "per_period",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 30L
    )
    testthat::expect_length(ad_dup$selection_fit$cells, 4L)
    testthat::expect_length(ad_dup$imputation_fit$cells, 4L)

    # Pair callback interface produces identical results to inline pair records
    grp_cb <- lapply(fx$groups, function(g) {
      g$pairs <- lapply(g$pairs, function(p) {
        p$observed <- NULL
        p$dy <- NULL
        p
      })
      g
    })
    cb_fn <- function(g_key, p_key) {
      orig <- fx$groups[[g_key]]$pairs[[p_key]]
      list(observed = orig$observed, dy = orig$dy)
    }
    ad_cb <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = grp_cb, route = "augmented",
      selection_model = "pair_logit", selection_pooling = "within_group",
      imputation_pooling = "within_group", imputation_scale = "per_period",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 20L,
      pair_callback = cb_fn
    )
    ad_inline <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "augmented",
      selection_model = "pair_logit", selection_pooling = "within_group",
      imputation_pooling = "within_group", imputation_scale = "per_period",
      glm_control = fx$glm_control, rank_tol = fx$rank_tol, max_chunk_rows = 20L
    )
    testthat::expect_equal(ad_cb$selection_fit$coefficients, ad_inline$selection_fit$coefficients, tolerance = 1e-15)
    testthat::expect_equal(ad_cb$imputation_fit$coefficients, ad_inline$imputation_fit$coefficients, tolerance = 1e-15)
  })

  testthat::test_that("clipping applies once to final product q and fit/channel/input refusals are enforced", {
    fx <- make_adapter_fixture()

    # Check that raw date_product probabilities below 1e-6 are floored once in static_obs_mean/cov
    # even when individual date probabilities q_s, q_u > 1e-6
    ad <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "weighting",
      selection_model = "independent_period_logit", selection_pooling = "separate",
      glm_control = fx$glm_control, max_chunk_rows = 30L
    )
    m_wt <- ad$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
    testthat::expect_equal(m_wt$kink$eps, ipw_clip_eps)
    testthat::expect_equal(m_wt$kink$n_q_floored, 0L)
    testthat::expect_true(grepl("conditional date independence", m_wt$diagnostics$assumption, fixed = TRUE))

    # Construct a fixture where individual date probabilities satisfy q_s > 1e-6 and q_u > 1e-6,
    # while their product satisfies q_s * q_u < 1e-6 so the final-product floor activates once
    grp_floor <- fx$groups
    inA_G <- which(grp_floor$g1_G$mask)
    obs_13 <- which(grp_floor$g1_G$mask & grp_floor$g1_G$pairs$p13$observed)
    unobs_t1 <- setdiff(inA_G, which(grp_floor$g1_G$dates$t1$observed))
    unobs_t3 <- setdiff(inA_G, which(grp_floor$g1_G$dates$t3$observed))
    u_ext <- obs_13[1L]
    # Construct a 2-covariate design where x2 keeps p ~ 0.7 on y=0 high-x1 rows
    # (driving beta_1 < -1), while u_ext (y=1) has x1 > 0, x2 = 0 so individual
    # q_s, q_u in (1e-6, 1e-3) and product q_s * q_u < 1e-6.
    make_floor_X <- function(obs_vec) {
      X <- matrix(0, nrow = fx$N, ncol = 1L, dimnames = list(fx$unit_id, "x1"))
      obs_in <- inA_G[obs_vec[inA_G]]
      unobs_in <- inA_G[!obs_vec[inA_G]]
      X[obs_in, 1L] <- 0.0
      X[unobs_in, 1L] <- 0.35
      X[u_ext, 1L] <- 0.95
      X
    }
    grp_floor$g1_G$dates$t1$selection_design <- make_floor_X(grp_floor$g1_G$dates$t1$observed)
    grp_floor$g1_G$dates$t3$selection_design <- make_floor_X(grp_floor$g1_G$dates$t3$observed)
    ad_fl <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = grp_floor, route = "weighting",
      selection_model = "independent_period_logit", selection_pooling = "separate",
      glm_control = fx$glm_control, max_chunk_rows = 30L
    )
    raw_pred <- static_nuisance_q_predict(
      ad_fl$selection_fit, fx$unit_id[u_ext], dates = c("g1_G::t1", "g1_G::t3")
    )
    testthat::expect_gt(raw_pred$q_dates[1L, "start"], ipw_clip_eps)
    testthat::expect_gt(raw_pred$q_dates[1L, "end"], ipw_clip_eps)
    testthat::expect_lt(unname(raw_pred$value), ipw_clip_eps)
    m_fl <- ad_fl$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
    testthat::expect_equal(m_fl$kink$n_q_floored, 1L)
    ref_fl <- ref_env$adapter_ref_edge_functional(
      N = fx$N, groups = grp_floor, group_key = "g1_G", pair_key = "p13",
      route = "weighting", selection_model = "independent_period_logit",
      selection_pooling = "separate", kind = "mean", h_vec = fx$h_one
    )
    testthat::expect_equal(unname(m_fl$estimate), unname(ref_fl), tolerance = 1e-10)
    testthat::expect_equal(
      m_fl$eval_if(fx$unit_id[1:15]),
      m_fl$eval_if()[1:15],
      tolerance = 1e-14
    )

    # Refusal: all-success selection cell (never silently replaced by known q = 1)
    grp_all1 <- fx$groups
    grp_all1$g2_C$pairs$p12$observed[] <- TRUE
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = grp_all1, route = "weighting",
        selection_model = "pair_logit", selection_pooling = "separate",
        glm_control = fx$glm_control, max_chunk_rows = 30L
      ),
      class = "cdid_error"
    )

    # Refusal: all-failure selection cell
    grp_all0 <- fx$groups
    grp_all0$g2_C$pairs$p12$observed[] <- FALSE
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = grp_all0, route = "weighting",
        selection_model = "pair_logit", selection_pooling = "separate",
        glm_control = fx$glm_control, max_chunk_rows = 30L
      ),
      class = "cdid_error"
    )

    # Refusal: missing or non-positive gap
    grp_nogap <- fx$groups
    grp_nogap$g1_G$pairs$p12$gap <- NULL
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = grp_nogap, route = "baseline",
        max_chunk_rows = 30L
      ),
      class = "cdid_static_obs_adapter_error"
    )

    # Refusal: reserved namespace in h_jacobian or missing h_coef_if
    testthat::expect_error(
      ad$mean("g1_G", "p13", h = fx$h_vec,
              h_jacobian = list("selection|g1_G" = fx$h_jacobian$trt_block),
              h_coef_if = list("selection|g1_G" = fx$h_coef_if_list$trt_block)),
      class = "cdid_static_obs_channel_error"
    )
    testthat::expect_error(
      ad$mean("g1_G", "p13", h = fx$h_vec, h_jacobian = fx$h_jacobian, h_coef_if = NULL),
      class = "cdid_static_obs_channel_error"
    )

    # Refusal: wrong-type or unknown IDs in eval_if
    testthat::expect_error(m_wt$eval_if(1:5), class = "cdid_static_nuisance_input_error")
    testthat::expect_error(m_wt$eval_if(c("unknown_id")), class = "cdid_static_nuisance_input_error")

    # Refusal: omitted, NULL, or non-positive max_chunk_rows; valid explicit positive bound accepted
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = fx$groups, route = "weighting",
        selection_model = "independent_period_logit", selection_pooling = "separate",
        glm_control = fx$glm_control
      ),
      class = "cdid_static_nuisance_input_error"
    )
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = fx$groups, route = "baseline"
      ),
      class = "cdid_static_obs_adapter_error"
    )
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = fx$groups, route = "weighting",
        selection_model = "independent_period_logit", selection_pooling = "separate",
        glm_control = fx$glm_control, max_chunk_rows = NULL
      ),
      class = "cdid_static_nuisance_input_error"
    )
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = fx$groups, route = "baseline",
        max_chunk_rows = NULL
      ),
      class = "cdid_static_obs_adapter_error"
    )
    testthat::expect_error(
      static_observation_adapter_fit(
        unit_id = fx$unit_id, groups = fx$groups, route = "weighting",
        selection_model = "independent_period_logit", selection_pooling = "separate",
        glm_control = fx$glm_control, max_chunk_rows = 0L
      ),
      class = "cdid_static_nuisance_input_error"
    )
    ad_chunk11 <- static_observation_adapter_fit(
      unit_id = fx$unit_id, groups = fx$groups, route = "weighting",
      selection_model = "independent_period_logit", selection_pooling = "separate",
      glm_control = fx$glm_control, max_chunk_rows = 11L
    )
    m_wt_11 <- ad_chunk11$mean("g1_G", "p13", h = fx$h_one, h_jacobian = list())
    testthat::expect_identical(ad_chunk11$max_chunk_rows, 11L)
    testthat::expect_equal(m_wt_11$estimate, m_wt$estimate, tolerance = 1e-15)
    testthat::expect_equal(m_wt_11$eval_if(), m_wt$eval_if(), tolerance = 1e-14)
  })

  testthat::test_that("operation-local with_if_batch shares q/j nuisance IFs once per ordered batch and cleans up on exit or error", {
    fx <- make_adapter_fixture()
    ad <- static_observation_adapter_fit(
      unit_id = fx$unit_id,
      groups = fx$groups,
      route = "augmented",
      selection_model = "independent_period_logit",
      selection_pooling = "within_group",
      imputation_pooling = "within_group",
      imputation_scale = "per_period",
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol,
      max_chunk_rows = 20L
    )

    h_calls <- 0L
    counted_h_cb <- function(ids) {
      h_calls <<- h_calls + 1L
      fx$h_coef_if_fn(ids)
    }

    m_g12 <- ad$mean("g1_G", "p12", h = fx$h_vec, h_jacobian = fx$h_jacobian, h_coef_if = counted_h_cb)
    m_g13 <- ad$mean("g1_G", "p13", h = fx$h_vec, h_jacobian = fx$h_jacobian, h_coef_if = counted_h_cb)
    m_c12 <- ad$mean("g1_C", "p12", h = fx$h_one, h_jacobian = list())
    c_c12 <- ad$cov("g1_C", "p12", Z = fx$Z)
    c_c13 <- ad$cov("g1_C", "p13", Z = fx$Z)
    c_g2  <- ad$cov("g2_C", "p12", Z = fx$Z)

    ns <- environment(static_observation_adapter_fit)
    orig_nuisance_if <- get("static_nuisance_if", envir = ns)
    was_locked <- base::bindingIsLocked("static_nuisance_if", ns)
    if (was_locked) base::unlockBinding("static_nuisance_if", ns)
    on.exit({
      if (base::bindingIsLocked("static_nuisance_if", ns)) base::unlockBinding("static_nuisance_if", ns)
      assign("static_nuisance_if", orig_nuisance_if, envir = ns)
      if (was_locked) base::lockBinding("static_nuisance_if", ns)
    }, add = TRUE)

    sel_calls <- 0L
    imp_calls <- 0L
    assign("static_nuisance_if", function(fit, ids) {
      if (identical(fit$kind, "logit")) {
        sel_calls <<- sel_calls + 1L
      } else if (identical(fit$kind, "wls")) {
        imp_calls <<- imp_calls + 1L
      }
      orig_nuisance_if(fit, ids)
    }, envir = ns)
    if (was_locked) base::lockBinding("static_nuisance_if", ns)

    # Non-sorted batch including 0-outcome (1, 2), 1-outcome (7), overlapping (18, 35, 14, 22), and out-of-role (58) units
    batch_pos_idx <- c(18L, 1L, 7L, 35L, 2L, 58L, 14L, 22L)
    batch_ids <- fx$unit_id[batch_pos_idx]

    eval_six <- function(ids) {
      list(
        m_g12 = m_g12$eval_if(ids, include_components = TRUE),
        m_g13 = m_g13$eval_if(ids, include_components = TRUE),
        m_c12 = m_c12$eval_if(ids),
        c_c12 = c_c12$eval_if(ids),
        c_c13 = c_c13$eval_if(ids),
        c_g2  = c_g2$eval_if(ids)
      )
    }

    # 1. Standalone uncached evaluation computes q and j IFs on every moment call (6 + 6)
    sel_calls <- 0L; imp_calls <- 0L; h_calls <- 0L
    uncached <- eval_six(batch_ids)
    testthat::expect_identical(sel_calls, 6L)
    testthat::expect_identical(imp_calls, 6L)
    testthat::expect_identical(h_calls, 2L)

    # 2. Admitted with_if_batch evaluates static_nuisance_if at most once for q and once for j
    #    across all 6 moments while still evaluating custom h_coef_if on every moment call
    sel_calls <- 0L; imp_calls <- 0L; h_calls <- 0L
    cached <- ad$with_if_batch(batch_ids, function() eval_six(batch_ids), enabled = TRUE)
    testthat::expect_identical(sel_calls, 1L)
    testthat::expect_identical(imp_calls, 1L)
    testthat::expect_identical(h_calls, 2L)
    testthat::expect_equal(cached, uncached, tolerance = 1e-14)
    testthat::expect_identical(cached$m_g12$influence, uncached$m_g12$influence)
    testthat::expect_identical(cached$c_c13, uncached$c_c13)

    # Holder is cleared immediately upon normal return; no per-edge matrices retained
    env_ad <- environment(ad$with_if_batch)
    testthat::expect_false(env_ad$batch_active)
    testthat::expect_null(env_ad$batch_pos)
    testthat::expect_null(env_ad$batch_sel_if)
    testthat::expect_null(env_ad$batch_imp_if)
    for (m_env in list(environment(m_g12$eval_if), environment(c_c12$eval_if))) {
      testthat::expect_false(exists("sel_if", envir = m_env, inherits = FALSE))
      testthat::expect_false(exists("imp_if", envir = m_env, inherits = FALSE))
      testthat::expect_false(exists("full_if", envir = m_env, inherits = FALSE))
    }

    # 3. Next batch / replay pass calculates anew once per fit (no cross-pass retention)
    sel_calls <- 0L; imp_calls <- 0L; h_calls <- 0L
    cached_replay <- ad$with_if_batch(batch_ids, function() eval_six(batch_ids), enabled = TRUE)
    testthat::expect_identical(sel_calls, 1L)
    testthat::expect_identical(imp_calls, 1L)
    testthat::expect_identical(h_calls, 2L)
    testthat::expect_identical(cached_replay, cached)

    # 4. Declined reuse (enabled = FALSE) and oversized request (> max_chunk_rows) leave holder inactive
    sel_calls <- 0L; imp_calls <- 0L; h_calls <- 0L
    declined <- ad$with_if_batch(batch_ids, function() {
      testthat::expect_false(env_ad$batch_active)
      eval_six(batch_ids)
    }, enabled = FALSE)
    testthat::expect_identical(sel_calls, 6L)
    testthat::expect_identical(imp_calls, 6L)
    testthat::expect_identical(h_calls, 2L)
    testthat::expect_equal(declined, uncached, tolerance = 1e-14)

    oversize_ids <- fx$unit_id[1:25]
    sel_calls <- 0L; imp_calls <- 0L
    over_res <- ad$with_if_batch(oversize_ids, function() {
      testthat::expect_false(env_ad$batch_active)
      m_c12$eval_if(oversize_ids)
    }, enabled = TRUE)
    testthat::expect_identical(sel_calls, 2L)
    testthat::expect_identical(imp_calls, 2L)
    testthat::expect_equal(over_res, m_c12$eval_if(oversize_ids), tolerance = 1e-14)

    # 5. Exact ordered-row miss (permuted order, subset), nested scope, and invalid row errors
    sel_calls <- 0L; imp_calls <- 0L
    ad$with_if_batch(batch_ids, function() {
      r1 <- m_g12$eval_if(batch_ids)
      testthat::expect_identical(sel_calls, 1L)
      testthat::expect_identical(imp_calls, 1L)

      # Permuted order of the same unit set does NOT hit and does NOT replace outer holder
      r_perm <- m_g13$eval_if(rev(batch_ids))
      testthat::expect_identical(sel_calls, 2L)
      testthat::expect_identical(imp_calls, 2L)
      testthat::expect_identical(env_ad$batch_pos, as.integer(batch_pos_idx))
      testthat::expect_equal(r_perm, uncached$m_g13$influence[rev(batch_ids)], tolerance = 1e-14)

      # Subset does NOT hit and nested with_if_batch declines without replacing outer holder
      r_sub <- ad$with_if_batch(batch_ids[1:4], function() {
        c_c12$eval_if(batch_ids[1:4])
      }, enabled = TRUE)
      testthat::expect_identical(sel_calls, 3L)
      testthat::expect_identical(imp_calls, 3L)
      testthat::expect_identical(env_ad$batch_pos, as.integer(batch_pos_idx))
      testthat::expect_equal(r_sub, uncached$c_c12[batch_ids[1:4], , drop = FALSE], tolerance = 1e-14)

      # Outer batch still hits after non-matching and nested calls
      r_hit <- c_c13$eval_if(batch_ids)
      testthat::expect_identical(sel_calls, 3L)
      testthat::expect_identical(imp_calls, 3L)
      testthat::expect_equal(r_hit, uncached$c_c13, tolerance = 1e-14)

      # Invalid row keys still raise cdid_static_nuisance_input_error
      testthat::expect_error(m_c12$eval_if(seq_along(batch_ids)), class = "cdid_static_nuisance_input_error")
      testthat::expect_error(m_c12$eval_if(c(batch_ids[1L], batch_ids[1L])), class = "cdid_static_nuisance_input_error")
      testthat::expect_error(m_c12$eval_if("unknown_unit_id"), class = "cdid_static_nuisance_input_error")
      r1
    }, enabled = TRUE)

    # 6. Forced callback error inside with_if_batch cleans up populated holder via on.exit
    m_err <- ad$mean(
      "g1_G", "p13",
      h = fx$h_vec,
      h_jacobian = fx$h_jacobian,
      h_coef_if = function(ids) stop("forced custom h_coef_if callback error")
    )
    sel_calls <- 0L; imp_calls <- 0L
    testthat::expect_error(
      ad$with_if_batch(batch_ids, function() m_err$eval_if(batch_ids), enabled = TRUE),
      "forced custom h_coef_if callback error"
    )
    # Confirm q and j matrices were actually computed before h_coef_if threw, and then cleaned up
    testthat::expect_identical(sel_calls, 1L)
    testthat::expect_identical(imp_calls, 1L)
    testthat::expect_false(env_ad$batch_active)
    testthat::expect_null(env_ad$batch_pos)
    testthat::expect_null(env_ad$batch_sel_if)
    testthat::expect_null(env_ad$batch_imp_if)
  })
})
