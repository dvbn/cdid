# ==============================================================================
# CP02 Custom-IPW Callback Bridge Tests (`test-cp02-custom-ipw.R`)
# ==============================================================================
# Discriminating checks for the CP02 custom-IPW callback seam on the D23
# candidate against the independent 174-case reference
# (`reference/cp02_custom_ipw_reference.R`), pooled-NYT with non-saturated
# covariates (nonzero direct + indirect propensity-rate and treated-rate
# sensitivities), full-frame row embedding, shared/reduction ownership,
# malformed/uncentered callback outputs, checked preflight before callback
# execution, and public request/pipeline/reporting/no-refit `aggte` seams.
# ==============================================================================

if (!exists("cp02_custom_reference_runs", inherits = TRUE)) {
  cp02_custom_reference_runs <- new.env(parent = emptyenv())
  cp02_custom_reference_runs$envs <- list()
}

cp02_load_reference_helpers <- function() {
  env <- new.env(parent = globalenv())
  env$counts <- list(callback_calls = 0L, newton_calls = 0L, newton_iterations = 0L)
  ref_path <- testthat::test_path("reference", "cp02_custom_ipw_reference.R")
  if (!file.exists(ref_path)) stop("Supply the pinned reference/cp02_custom_ipw_reference.R test file; no replacement reference is generated.")
  exprs <- parse(ref_path)
  for (nm in c("newton", "jacobian", "custom_ipw", "compose", "full_stack")) {
    take <- vapply(exprs, function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name(nm)), logical(1L))
    stopifnot(sum(take) == 1L)
    eval(exprs[[which(take)]], envir = env)
  }
  # Diagnostic registry for the affected runner; includes failed test attempts.
  cp02_custom_reference_runs$envs <- c(cp02_custom_reference_runs$envs, list(env))
  env
}

cp02_custom_clone <- function(overrides = list()) {
  ns <- environment(chained_public_att_gt)
  e <- new.env(parent = ns)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, ns)
    if (is.function(f) && identical(environment(f), ns)) {
      environment(f) <- e
      assign(nm, f, e)
    }
  }
  list2env(overrides, e)
  e
}

test_that("CP02 custom-IPW matches the 174-check independent reference across all four panel24 regimes at 1e-10", {
  ref_env <- cp02_load_reference_helpers()
  panel <- utils::read.csv(testthat::test_path("fixtures", "panel24.csv"))
  ids <- sort(unique(panel$id))
  n <- length(ids)
  calendar <- sort(unique(panel$time))
  first <- match(ids, panel$id)
  G0 <- panel$g[first]
  sign <- panel$sign[first]
  Y <- matrix(NA_real_, n, length(calendar))
  Y[cbind(match(panel$id, ids), match(panel$time, calendar))] <- panel$y
  observed <- !is.na(Y)

  regimes <- list(
    NT_unit = list(
      G = G0,
      w = rep(1, n),
      X = cbind(intercept = rep(1, n)),
      cg = "nevertreated"
    ),
    NT_weighted_covariate = list(
      G = G0,
      w = 1 + (ids %% 5) / 4,
      X = cbind(intercept = 1, sign = sign),
      cg = "nevertreated"
    ),
    NYT_unit = list(
      G = replace(G0, ids %in% c(4, 8, 16), 4),
      w = rep(1, n),
      X = cbind(intercept = rep(1, n)),
      cg = "notyettreated"
    ),
    NYT_weighted = list(
      G = replace(G0, ids %in% c(4, 8, 16), 4),
      w = 1 + (ids %% 5) / 4,
      X = cbind(intercept = rep(1, n)),
      cg = "notyettreated"
    )
  )

  omega <- function(phi) crossprod(sweep(phi, 2, colMeans(phi))) / nrow(phi)

  for (regime in names(regimes)) {
    fx <- regimes[[regime]]
    G <- fx$G
    w <- fx$w
    X <- fx$X
    dat <- panel
    dat$g <- G[match(dat$id, ids)]
    pair <- construct_observed_pairs(dat, "id", "time", "y", "g")
    ctx <- construct_ipw_context(pair, ids, X, w)
    expect_identical(ctx$a_id, ids)
    expect_identical(ctx$n_a_units, n)

    occ <- do.call(rbind, lapply(utils::combn(calendar, 2, simplify = FALSE), function(e) {
      data.frame(
        cohort = sort(unique(G)),
        s = e[1L],
        t = e[2L],
        count = vapply(
          sort(unique(G)),
          function(g) sum(G == g & observed[, e[1L]] & observed[, e[2L]]),
          numeric(1L)
        )
      )
    }))
    cat <- build_comparison_catalog(
      calendar, sort(unique(G)), 3, 0, occ,
      fx$cg, "fixed-baseline-reference", NULL,
      list(
        max_candidates = 1000, max_rows = 1000, max_pool_entries = 10000,
        max_scratch = 1e7, max_work = 1e8
      )
    )
    rows <- cat$rows
    pool <- function(j) {
      cat$spec$cohorts[cat$pools$entries[seq.int(
        cat$pools$start[j], length.out = cat$pools$size[j]
      )]]
    }
    L <- nrow(rows)
    expect_gt(L, 0L)

    custom_fits <- vector("list", L)
    builtin_fits <- vector("list", L)
    ref_composed <- matrix(0, n, L)
    ref_stack <- matrix(0, n, L)
    builtin_phi <- matrix(0, n, L)
    ref_mu <- numeric(L)
    stack_mu <- numeric(L)
    builtin_mu <- numeric(L)
    raw_point_errors <- numeric()
    z_perm <- NULL

    for (j in seq_len(L)) {
      K <- pool(rows$pool_id[j])
      cohorts <- c(3, K)
      gi <- match(G, cohorts)
      s <- rows$s[j]
      t <- rows$t[j]
      S <- observed[, s] & observed[, t]
      obs <- which(S & !is.na(gi))
      dy <- numeric(n)
      dy[obs] <- Y[obs, t] - Y[obs, s]
      I <- outer(G, cohorts, "==") * 1
      rate <- colSums(I * w * S) / colSums(I * w)

      # Independent reference with reversed observer order vs production canonical A-unit order
      obs_rev <- rev(obs)
      ans_rev <- ref_env$custom_ipw(
        Y[obs_rev, t], Y[obs_rev, s], as.numeric(G[obs_rev] == 3),
        X[obs_rev, , drop = FALSE], w[obs_rev] / rate[gi[obs_rev]], TRUE
      )
      z <- ref_env$compose(ans_rev, obs_rev, gi, w, S, rate, n)
      st <- ref_env$full_stack(G, S, dy, X, w, cohorts)

      # Production custom-IPW link fit (canonical increasing A-unit order)
      cb_before <- ref_env$counts$callback_calls
      cfit <- fit_custom_ipw_link(
        ctx, 3, K, s, t,
        custom_est_method = ref_env$custom_ipw,
        link_id = paste0("link_", j),
        functional_key = rows$functional_key[j],
        max_entries = 1e5
      )
      expect_identical(ref_env$counts$callback_calls, cb_before + 1L)
      expect_true(cfit$accounting$observer_released)
      expect_false(cfit$preflight$peak_rss_promised)
      expect_false(cfit$shared_across_dates)
      expect_false(cfit$nuisance_mask_verified)

      # Built-in IPW link fit for comparison
      bfit <- fit_ipw_link(
        ctx, 3, K, s, t,
        stats::glm.control(epsilon = 1e-12, maxit = 50, trace = FALSE),
        link_id = paste0("link_", j),
        functional_key = rows$functional_key[j]
      )
      replay_ipw_if(ctx, bfit, 7, function(batch) builtin_phi[batch$a_index, j] <<- batch$phi[, 1L])

      custom_fits[[j]] <- cfit
      builtin_fits[[j]] <- bfit
      ref_composed[, j] <- z$phi
      ref_stack[, j] <- st$phi
      ref_mu[j] <- z$estimate
      stack_mu[j] <- st$estimate
      builtin_mu[j] <- bfit$estimate

      expect_equal(unname(cfit$estimate), unname(z$estimate), tolerance = 1e-10, info = paste(regime, "point vs compose", j))
      expect_equal(unname(cfit$estimate), unname(st$estimate), tolerance = 1e-10, info = paste(regime, "point vs full_stack", j))
      expect_equal(unname(cfit$estimate), unname(bfit$estimate), tolerance = 1e-10, info = paste(regime, "point vs builtin", j))
      expect_equal(unname(cfit$b_rate), unname(z$sensitivity), tolerance = 1e-10, info = paste(regime, "b_rate vs compose", j))
      expect_equal(unname(cfit$b_rate), unname(st$sensitivity), tolerance = 1e-10, info = paste(regime, "b_rate vs full_stack", j))

      if (j == 1L && fx$cg == "notyettreated") {
        ans_perm <- ans_rev
        in_g <- which(G[obs_rev] == 3)
        ans_perm$att.inf.func[in_g] <- rev(ans_perm$att.inf.func[in_g])
        z_perm <- ref_env$compose(ans_perm, obs_rev, gi, w, S, rate, n)
        expect_equal(unname(z_perm$sensitivity), unname(st$sensitivity), tolerance = 1e-10)
        expect_gt(max(abs(z_perm$phi - st$phi)), 1e-4)
      }
      if (fx$cg == "notyettreated" && length(K) > 1L) {
        expect_gt(max(abs(z$base - st$phi)), 1e-4)
        raw <- ref_env$custom_ipw(
          Y[obs_rev, t], Y[obs_rev, s], as.numeric(G[obs_rev] == 3),
          X[obs_rev, , drop = FALSE], w[obs_rev], TRUE
        )
        raw_point_errors <- c(raw_point_errors, abs(raw$ATT - st$estimate))
      }
    }

    if (length(raw_point_errors)) {
      expect_gt(max(raw_point_errors), 1e-4)
    }

    # Replay & materialization via production custom-IPW helpers
    mat_custom <- materialize_custom_ipw_if(
      ctx, custom_fits,
      max_entries = 1e5,
      chunk_size = 7,
      include_components = TRUE
    )
    cov_custom <- stream_custom_ipw_covariance(ctx, custom_fits, chunk_size = 7)
    cov_builtin <- stream_ipw_covariance(ctx, builtin_fits, chunk_size = 7)

    expect_equal(unname(mat_custom$estimates), unname(ref_mu), tolerance = 1e-10, info = paste(regime, "estimates"))
    expect_equal(unname(mat_custom$phi_matrix), unname(ref_composed), tolerance = 1e-10, info = paste(regime, "joint IF vs compose"))
    expect_equal(unname(mat_custom$phi_matrix), unname(ref_stack), tolerance = 1e-10, info = paste(regime, "joint IF vs full_stack"))
    expect_equal(unname(mat_custom$phi_matrix), unname(builtin_phi), tolerance = 1e-10, info = paste(regime, "joint IF vs builtin"))
    expect_equal(
      unname(mat_custom$components$base + mat_custom$components$rate),
      unname(mat_custom$phi_matrix),
      tolerance = 1e-14,
      info = paste(regime, "base + rate decomposition")
    )
    expect_equal(unname(cov_custom$covariance), unname(omega(ref_stack)), tolerance = 1e-10, info = paste(regime, "cov vs full_stack"))
    expect_equal(unname(cov_custom$covariance), unname(cov_builtin$covariance), tolerance = 1e-10, info = paste(regime, "cov vs builtin"))

    if (fx$cg == "notyettreated") {
      wrong_within <- mat_custom$phi_matrix
      wrong_within[, 1L] <- z_perm$phi
      expect_gt(max(abs(omega(wrong_within) - cov_custom$covariance)), 1e-4)
    }

    # Full-frame row embedding (N = n + 2, reversed A-unit order + 2 non-A rows)
    N <- n + 2L
    frame_ids <- c(-101, rev(ids), -102)
    map <- match(ids, frame_ids)
    embedded <- matrix(0, N, L)
    embedded[map, ] <- (N / n) * mat_custom$phi_matrix
    transported <- transport_if_to_full_frame(mat_custom$phi_matrix, N, map)$Phi_full
    expect_equal(
      unname(c(colMeans(mat_custom$phi_matrix), colMeans(transported))),
      rep(0, 2L * L),
      tolerance = 1e-10,
      info = paste(regime, "sample-root IF centering")
    )
    expect_equal(unname(transported), unname(embedded), tolerance = 1e-10, info = paste(regime, "full-frame transport"))
    expect_equal(
      unname(omega(transported) / N),
      unname(cov_custom$covariance / n),
      tolerance = 1e-10,
      info = paste(regime, "full-frame variance normalization")
    )
    expect_gt(max(abs(embedded * (n / N) - transported)), 1e-4)

    if (isNamespace(environment(build_custom_ipw_catalog_links))) {
    # Actual owner/view-backed bridge: all D10/D13 rows, with the same independent
    # full equation-stack reference as above (no pipeline or constructor doubles).
    own <- ipw_owner_control(7, 2e6, 2e6, 5e8, 5e5)
    pt <- ipw_point_view_control(7, 2000, 2e6, 2e6, 5e8, 5e5)
    po <- construct_ipw_panel_owner(dat, "id", "time", "y", "g", own)
    fo <- construct_ipw_frame_owner(po, ids, X, w, own)
    meta <- build_owner_comparison_catalog(po, 3, 0, fx$cg, "fixed-baseline-reference",
      control = own, frame_owner = fo, build_design = TRUE,
      declared_qualified = data.frame(g = 3, t = 3:5), target_post_times = 3:5,
      catalog_bounds = list(max_candidates = 1e5, max_rows = 1e5,
        max_pool_entries = 1e6, max_scratch = 1e7, max_work = 1e7))
    expect_identical(meta$catalog$rows, cat$rows)
    cb_before <- ref_env$counts$callback_calls
    links <- build_custom_ipw_catalog_links(meta, fo, own, pt, ref_env$custom_ipw)
    expect_identical(ref_env$counts$callback_calls, cb_before + L)
    view_phi <- matrix(0, n, L)
    replay_custom_ipw_if_view(links$point_links, pt, callback_schema = TRUE,
      callback = function(batch) view_phi[batch$a_index, ] <<- batch$phi)
    view_cov <- stream_custom_ipw_covariance_view(links$point_links, pt)
    expect_equal(unname(view_phi), unname(ref_stack), tolerance = 1e-10)
    expect_equal(unname(view_cov$covariance), unname(omega(ref_stack)), tolerance = 1e-10)
    expect_equal(unname(view_cov$estimates), unname(stack_mu), tolerance = 1e-10)
    cert <- certify_ipw_structural_classes_view(links$point_links, meta$design,
      match(seq_len(L), meta$design$links$catalog_row), "optimal",
      list(max_links = 500, max_classes = 50, max_class_dates = 50,
           max_E_entries = 5000, max_mask_entries = 100000), own)
    expect_equal(unname(apply_ipw_reduction_map(cert, view_cov$estimates[cert$retained_index])),
                 unname(stack_mu), tolerance = 1e-10)
    expect_equal(unname(t(apply_ipw_reduction_map(cert, t(view_phi[, cert$retained_index, drop = FALSE])))),
                 unname(ref_stack), tolerance = 1e-10)

    }

    # Native GMM is an installed check. Source helpers above run in every regime.
    if (!isNamespace(environment(solve_chained_gmm_core))) next

    # Fixed-system GMM with streamed custom-IPW IF replay
    dq <- data.frame(g = 3, t = 3:5)
    design <- build_comparison_catalog_design(cat, dq, 3:5)
    W <- materialize_chained_design(design, 1e5)$W
    p <- ncol(W)
    V <- cov_custom$covariance
    KKT <- rbind(cbind(V, W), cbind(t(W), matrix(0, p, p)))
    sol <- tryCatch(solve(KKT, rbind(diag(L), matrix(0, p, L))), error = identity)

    prov <- gmm_covariance_provenance(
      "unit", n, n, n, TRUE, FALSE, "none", 1 / n,
      "custom-ipw-test", "custom-ipw-stack"
    )
    replay <- gmm_if_replay_spec(
      function(consume_batch) {
        replay_custom_ipw_if(ctx, custom_fits, chunk_size = 7, callback = function(b) {
          consume_batch(list(unit_id = b$id, phi = b$phi))
        })
      },
      "unit", FALSE, n, n, n, TRUE, FALSE, "none", 1 / n,
      "custom-ipw-test", "custom-ipw-stack",
      expected_unit_ids = ids
    )
    carried <- chained_public_carried_settings()$gmm_guard
    accuracy <- do.call(gmm_accuracy_control, c(
      list(
        method = "conventional_chol_qr",
        target_error_unit = "abs",
        max_target_error_abs = chained_public_staging_envelope(),
        max_target_error_se = NA_real_
      ),
      carried,
      list(
        require_matching_scale = TRUE,
        max_loading_error_abs = chained_public_staging_envelope()
      )
    ))
    bounds <- gmm_bounds_control(1e7, n, 100, n, TRUE, 1e5)
    C <- diag(p)
    colnames(C) <- paste0("parameter_", seq_len(p))

    fitg <- tryCatch(
      solve_chained_gmm_core(
        design, unname(mat_custom$estimates), V, prov, C, "optimal", replay,
        accuracy, bounds, "lowest",
        request_param_cov = TRUE,
        materialize_target_if = TRUE
      ),
      error = identity
    )

    if (identical(regime, "NT_weighted_covariate")) {
      # D5: genuine finite-fixture covariance/Cholesky refusal is preserved (no pinv/ridge/moment drop)
      expect_s3_class(fitg, "cdid_gmm_pd_error")
    } else {
      expect_false(inherits(sol, "error"))
      expect_false(inherits(fitg, "error"))
      B <- sol[L + seq_len(p), , drop = FALSE]
      theta <- drop(B %*% stack_mu)
      psi <- ref_stack %*% t(B)
      Sigma <- omega(psi) / n
      expect_equal(unname(fitg$targets$estimate), unname(theta), tolerance = 1e-10, info = paste(regime, "GMM point"))
      expect_equal(
        unname(fitg$if_projection$materialized_target_if),
        unname(psi),
        tolerance = 1e-10,
        info = paste(regime, "GMM joint IF")
      )
      expect_equal(
        unname(fitg$targets$se_sandwich / sqrt(n)),
        unname(sqrt(diag(Sigma))),
        tolerance = 1e-10,
        info = paste(regime, "GMM sampling SE")
      )
      expect_equal(
        unname(fitg$if_projection$sandwich_cov),
        unname(Sigma * n),
        tolerance = 1e-10,
        info = paste(regime, "GMM core asymptotic covariance")
      )
    }
  }
})

test_that("CP02 custom-IPW pooled-NYT covariate checks verify nonzero indirect/treated-rate sensitivities, 1-col matrix IF, and shared/reduction ownership", {
  ref_env <- cp02_load_reference_helpers()
  panel <- utils::read.csv(testthat::test_path("fixtures", "panel24.csv"))
  ids <- sort(unique(panel$id))
  n <- length(ids)
  calendar <- sort(unique(panel$time))
  first <- match(ids, panel$id)
  G0 <- panel$g[first]
  sign <- panel$sign[first]
  Y <- matrix(NA_real_, n, length(calendar))
  Y[cbind(match(panel$id, ids), match(panel$time, calendar))] <- panel$y
  observed <- !is.na(Y)

  G_e <- replace(G0, ids %in% c(4, 8, 16), 4)
  w_e <- 1 + (ids %% 5) / 4
  X_e <- cbind(intercept = 1, x = sign + (ids %% 3 - 1) / 4)
  dat_e <- panel
  dat_e$g <- G_e[match(dat_e$id, ids)]
  pair_e <- construct_observed_pairs(dat_e, "id", "time", "y", "g")
  ctx_e <- construct_ipw_context(pair_e, ids, X_e, w_e)

  K_e <- c(4, Inf)
  cohorts_e <- c(3, K_e)
  gi_e <- match(G_e, cohorts_e)
  S_e <- observed[, 1L] & observed[, 2L]
  obs_e <- which(S_e & !is.na(gi_e))
  dy_e <- numeric(n)
  dy_e[obs_e] <- Y[obs_e, 2L] - Y[obs_e, 1L]

  ref_e <- ref_env$full_stack(G_e, S_e, dy_e, X_e, w_e, cohorts_e)
  ctrl_e <- stats::glm.control(epsilon = 1e-12, maxit = 50, trace = FALSE)
  fit_builtin_e <- fit_ipw_link(ctx_e, 3, K_e, 1, 2, ctrl_e)
  phi_builtin_e <- materialize_ipw_if(ctx_e, fit_builtin_e, max_entries = 1000, chunk_size = 7)$phi

  sess_e <- construct_ipw_shared_session(ctx_e, max_nuisance = 4, max_links = 4, max_unit_entries = 1000)
  nu_e <- construct_ipw_nuisance(sess_e, 3, K_e, 1, 2, ctrl_e)
  decl_e <- declare_ipw_shared_link(sess_e, nu_e, 3, K_e, 1, 2, ctrl_e)
  fit_sh_e <- fit_ipw_link_shared(sess_e, decl_e)
  phi_shared_e <- materialize_ipw_if(ctx_e, fit_sh_e, max_entries = 1000, chunk_size = 7)$phi

  # Custom-IPW fit with vector att.inf.func and with 1-column matrix att.inf.func (plus forwarded `...`)
  seen_dots <- NULL
  cb_vec <- function(y1, y0, D, covariates, i.weights, inffunc, bonus_arg = NULL) {
    seen_dots <<- bonus_arg
    ref_env$custom_ipw(y1, y0, D, covariates, i.weights, inffunc)
  }
  cb_mat <- function(y1, y0, D, covariates, i.weights, inffunc) {
    ans <- ref_env$custom_ipw(y1, y0, D, covariates, i.weights, inffunc)
    list(ATT = ans$ATT, att.inf.func = matrix(ans$att.inf.func, ncol = 1L))
  }

  fit_c_vec <- fit_custom_ipw_link(
    ctx_e, 3, K_e, 1, 2,
    custom_est_method = cb_vec,
    custom_dots = list(bonus_arg = "forwarded_ok"),
    link_id = "nyt_e"
  )
  expect_identical(seen_dots, "forwarded_ok")
  fit_c_mat <- fit_custom_ipw_link(
    ctx_e, 3, K_e, 1, 2,
    custom_est_method = cb_mat,
    link_id = "nyt_e_mat"
  )
  phi_c_vec <- materialize_custom_ipw_if(ctx_e, fit_c_vec, max_entries = 1000, chunk_size = 7)$phi
  phi_c_mat <- materialize_custom_ipw_if(ctx_e, fit_c_mat, max_entries = 1000, chunk_size = 7)$phi

  expect_equal(unname(fit_c_vec$estimate), unname(ref_e$estimate), tolerance = 1e-10)
  expect_equal(unname(fit_c_vec$estimate), unname(fit_sh_e$estimate), tolerance = 1e-10)
  expect_equal(unname(fit_c_vec$b_rate), unname(ref_e$sensitivity), tolerance = 1e-10)
  expect_equal(unname(phi_c_vec), unname(ref_e$phi), tolerance = 1e-10)
  expect_equal(unname(phi_c_vec), unname(phi_builtin_e), tolerance = 1e-10)
  expect_equal(unname(phi_c_vec), unname(phi_shared_e), tolerance = 1e-10)
  expect_identical(fit_c_mat$estimate, fit_c_vec$estimate)
  expect_identical(fit_c_mat$b_rate, fit_c_vec$b_rate)
  expect_identical(phi_c_mat, phi_c_vec)

  # Direct + indirect rate sensitivity decomposition (nonzero indirect and nonzero treated-rate b_{ell, g})
  k_e <- length(cohorts_e)
  p_e <- ncol(X_e)
  G_r_sens <- -ref_e$J[k_e + p_e + 1L, seq_len(k_e)] / ref_e$J[k_e + p_e + 1L, k_e + p_e + 1L] +
    ref_e$J[k_e + p_e + 2L, seq_len(k_e)] / ref_e$J[k_e + p_e + 2L, k_e + p_e + 2L]
  indirect_sens <- -drop(
    (ref_e$J[k_e + p_e + 2L, k_e + seq_len(p_e)] / ref_e$J[k_e + p_e + 2L, k_e + p_e + 2L]) %*%
      solve(ref_e$J[k_e + seq_len(p_e), k_e + seq_len(p_e)], ref_e$J[k_e + seq_len(p_e), seq_len(k_e)])
  )
  expect_equal(unname(G_r_sens + indirect_sens), unname(fit_c_vec$b_rate), tolerance = 1e-10)
  expect_gt(max(abs(indirect_sens)), 1e-4)
  expect_gt(abs(fit_c_vec$b_rate[1L]), 1e-4)

  if (isNamespace(environment(fit_custom_ipw_link_view))) {
  own <- ipw_owner_control(7, 2e6, 2e6, 5e8, 5e5)
  pt <- ipw_point_view_control(7, 2000, 2e6, 2e6, 5e8, 5e5)
  po <- construct_ipw_panel_owner(dat_e, "id", "time", "y", "g", own)
  fo <- construct_ipw_frame_owner(po, ids, X_e, w_e, own)
  lv <- construct_ipw_link_view(fo, 1, 2, 3, K_e, control = own, link_id = "nonsaturated")
  cv <- fit_custom_ipw_link_view(fo, lv, pt, ref_env$custom_ipw)
  phi_view <- numeric(n)
  replay_custom_ipw_if_view(cv, pt, callback = function(batch) phi_view[batch$a_index] <<- batch$phi[, 1L])
  expect_equal(unname(cv$b_rate), unname(ref_e$sensitivity), tolerance = 1e-10)
  expect_equal(unname(phi_view), unname(ref_e$phi), tolerance = 1e-10)
  expect_true(all(cv$obs_a_idx %in% seq_len(n)))
  rate_only <- setdiff(which(G_e %in% cohorts_e), cv$obs_a_idx)
  expect_gt(length(rate_only), 0L)
  expect_gt(max(abs(phi_view[rate_only])), 1e-4)

  }

  # Unweighted propensity with weighted outcome means fails on non-saturated covariates,
  # but coincides on intercept-only covariates (confirming the refined TeX characterization).
  cb_unweighted_prop <- function(y1, y0, D, covariates, i.weights, inffunc) {
    X <- covariates
    d <- D
    w <- i.weights
    dy <- y1 - y0
    b <- ref_env$newton(X, d, rep(1, length(d)))
    p <- ncol(X)
    odds <- exp(drop(X %*% b))
    a <- w * d
    c_w <- w * (1 - d) * odds
    T_m <- sum(a * dy) / sum(a)
    M_m <- sum(c_w * dy) / sum(c_w)
    eq <- function(th) {
      eta <- drop(X %*% th[seq_len(p)])
      cbind(
        X * (d - 1 / (1 + exp(-eta))),
        a * (dy - th[p + 1L]),
        w * (1 - d) * exp(eta) * (dy - th[p + 2L])
      )
    }
    theta <- c(b, T_m, M_m)
    xi <- t(solve(ref_env$jacobian(theta, eq), -t(eq(theta))))
    list(ATT = T_m - M_m, att.inf.func = xi[, p + 1L] - xi[, p + 2L])
  }
  fit_unw_cov <- fit_custom_ipw_link(ctx_e, 3, K_e, 1, 2, custom_est_method = cb_unweighted_prop)
  expect_gt(abs(fit_unw_cov$estimate - ref_e$estimate), 1e-4)
  expect_gt(max(abs(fit_unw_cov$b_rate - ref_e$sensitivity)), 1e-4)

  ctx_int <- construct_ipw_context(pair_e, ids, cbind(intercept = rep(1, n)), w_e)
  fit_w_int <- fit_custom_ipw_link(ctx_int, 3, K_e, 1, 2, custom_est_method = ref_env$custom_ipw)
  fit_u_int <- fit_custom_ipw_link(ctx_int, 3, K_e, 1, 2, custom_est_method = cb_unweighted_prop)
  expect_equal(fit_u_int$estimate, fit_w_int$estimate, tolerance = 1e-12)

  # Shared/reduction ownership and seal enforcement
  expect_error(
    declare_ipw_shared_link(sess_e, nu_e, 3, K_e, 2, 3, ctrl_e),
    class = "cdid_ipw_shared_mask_error"
  )
  expect_error(
    replay_ipw_if(ctx_e, fit_c_vec, 7, function(b) NULL),
    class = "cdid_ipw_input_error"
  )
  tampered <- fit_c_vec
  tampered$estimate <- tampered$estimate + 1
  expect_error(
    validate_custom_ipw_link(tampered),
    class = "cdid_custom_ipw_seal_error"
  )
})

test_that("CP02 custom-IPW refuses malformed callback outputs, never centers away an uncentered IF, and enforces preflight before calls", {
  ref_env <- cp02_load_reference_helpers()
  panel <- utils::read.csv(testthat::test_path("fixtures", "panel24.csv"))
  ids <- sort(unique(panel$id))
  n <- length(ids)
  pair <- construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- construct_ipw_context(pair, ids, cbind(intercept = rep(1, n)), rep(1, n))

  bad_callbacks <- list(
    not_a_list = list(fn = function(...) 1.0, cls = "cdid_custom_ipw_output_error"),
    data_frame = list(fn = function(...) data.frame(ATT = 0, att.inf.func = 0), cls = "cdid_custom_ipw_output_error"),
    missing_fields = list(fn = function(...) list(ATT = 0), cls = "cdid_custom_ipw_output_error"),
    na_att = list(fn = function(y1, ...) list(ATT = NA_real_, att.inf.func = rep(0, length(y1))), cls = "cdid_custom_ipw_output_error"),
    inf_att = list(fn = function(y1, ...) list(ATT = Inf, att.inf.func = rep(0, length(y1))), cls = "cdid_custom_ipw_output_error"),
    vec_att = list(fn = function(y1, ...) list(ATT = c(0, 1), att.inf.func = rep(0, length(y1))), cls = "cdid_custom_ipw_output_error"),
    wrong_len_if = list(fn = function(y1, ...) list(ATT = 0, att.inf.func = rep(0, length(y1) - 1L)), cls = "cdid_custom_ipw_dimension_error"),
    two_col_if = list(fn = function(y1, ...) list(ATT = 0, att.inf.func = matrix(0, length(y1), 2L)), cls = "cdid_custom_ipw_dimension_error"),
    nonfinite_if = list(fn = function(y1, ...) list(ATT = 0, att.inf.func = c(Inf, rep(0, length(y1) - 1L))), cls = "cdid_custom_ipw_output_error")
  )
  for (nm in names(bad_callbacks)) {
    spec <- bad_callbacks[[nm]]
    expect_error(
      fit_custom_ipw_link(ctx, 3, Inf, 1, 2, custom_est_method = spec$fn, link_id = nm),
      class = spec$cls,
      info = nm
    )
  }

  # An intentionally uncentered IF is NEVER silently centered away:
  # adding constant delta = 0.25 to psi shifts b_{ell, k} by -0.25 * (m_k / m) / r_k and preserves exact psi_obs.
  fit_base <- fit_custom_ipw_link(ctx, 3, Inf, 1, 2, custom_est_method = ref_env$custom_ipw)
  delta_shift <- 0.25
  cb_uncentered <- function(y1, y0, D, covariates, i.weights, inffunc) {
    ans <- ref_env$custom_ipw(y1, y0, D, covariates, i.weights, inffunc)
    ans$att.inf.func <- ans$att.inf.func + delta_shift
    ans
  }
  fit_unc <- fit_custom_ipw_link(ctx, 3, Inf, 1, 2, custom_est_method = cb_uncentered)
  expect_equal(fit_unc$psi_obs, fit_base$psi_obs + delta_shift, tolerance = 1e-15)
  expect_equal(fit_unc$diagnostics$callback_if_mean,
               fit_base$diagnostics$callback_if_mean + delta_shift, tolerance = 1e-15)
  # This deliberately invalid IF probes preservation/diagnostics, not estimator validity.
  phi_unc <- materialize_custom_ipw_if(ctx, fit_unc, max_entries = 1000, chunk_size = 7)$phi
  expect_equal(mean(phi_unc), delta_shift, tolerance = 1e-14)
  obs_share <- fit_base$counts$by_cohort_pair / fit_base$counts$n_observed_risk_units
  expect_equal(
    unname(fit_unc$b_rate),
    unname(fit_base$b_rate - delta_shift * obs_share / fit_base$masses$rates),
    tolerance = 1e-14
  )

  # Preflight bound check runs BEFORE invoking the user callback
  called <- 0L
  cb_spy <- function(...) {
    called <<- called + 1L
    list(ATT = 0, att.inf.func = rep(0, 10))
  }
  expect_error(
    fit_custom_ipw_link(ctx, 3, Inf, 1, 2, custom_est_method = cb_spy, max_entries = 5),
    class = "cdid_ipw_bound_error"
  )
  expect_identical(called, 0L)

  expect_error(
    materialize_custom_ipw_if(ctx, fit_base, max_entries = 5, chunk_size = 7),
    class = "cdid_ipw_bound_error"
  )
})

test_that("CP02 source request and no-refit aggte seams forward custom est_method and ... using pipeline/constructor test doubles", {
  z <- readRDS(testthat::test_path("fixtures", "cp02_result_adapter_N200.rds"))
  r <- z$reporting
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L

  custom_fn <- function(y1, y0, D, covariates, i.weights, inffunc = TRUE, tune_lambda = 0) {
    list(ATT = 0, att.inf.func = rep(0, length(y1)))
  }

  spy <- cp02_custom_clone(list(
    run_chained_ipw_pipeline = function(...) {
      a <- list(...)
      calls$n <- calls$n + 1L
      calls$args <- a
      expect_identical(a$custom_est_method, custom_fn)
      expect_identical(a$custom_dots, list(tune_lambda = 0.75, label = "bridge"))
      expect_false(a$run_numerical_transfer)
      structure(
        list(
          reporting = r,
          panel_owner = list(
            a_id = z$frame$id,
            n_a_units = 200L,
            all_times = z$calendar,
            cohorts = c(3, 4, 5, Inf)
          ),
          numerical_transfer = NULL,
          support = z$support,
          reduction_summary = list(
            weighting = "optimal",
            reduce_structural_redundancy = TRUE,
            reduction_active = FALSE,
            n_links_all = 10L,
            n_links_used = 10L,
            max_moment_expansion_diff = 0
          )
        ),
        class = "cdid_ipw_pipeline_result"
      )
    },
    chained_did_AGGTEobj = function(aggregation, ...) list(aggregation = aggregation, chained = list()),
    chained_did_MP = function(reporting, full_id, call, params, units, ...) {
      structure(
        list(
          att = reporting$estimates,
          inffunc = reporting$if_full,
          DIDparams = params,
          chained = list(
            frame = units,
            available = reporting$available,
            target_id = reporting$target_meta$target_id,
            target_meta = reporting$target_meta
          )
        ),
        class = c("cdid_MP", "MP")
      )
    }
  ))

  d <- expand.grid(id = z$frame$id, period = z$calendar)
  d$G <- z$frame$G[d$id]
  d$x <- z$frame$x[d$id]
  d$y <- d$id + d$period

  dot_evals <- 0L
  mp <- spy$att_gt(
    "y", "period", "id", "G", ~x, d,
    est_method = custom_fn,
    bstrap = FALSE, cband = FALSE,
    tune_lambda = { dot_evals <- dot_evals + 1L; 0.75 }, label = "bridge"
  )
  expect_identical(calls$n, 1L)
  expect_identical(dot_evals, 1L)
  expect_identical(mp$DIDparams$est_method, custom_fn)
  expect_identical(mp$DIDparams$extra_args, list(tune_lambda = 0.75, label = "bridge"))
  expect_identical(mp$chained$request$est_method, "custom")
  expect_match(mp$chained$public$entry, "custom IPW callback", fixed = TRUE)
  expect_true(mp$chained$reduction_summary$reduce_structural_redundancy)
  expect_false(mp$chained$reduction_summary$reduction_active)

  for (agg_type in c("simple", "group", "dynamic", "calendar")) {
    agg <- spy$chained_aggte_request(mp, type = agg_type)
    expect_identical(calls$n, 1L) # no refit on any aggregation
    expect_equal(
      unname(agg$aggregation$estimates),
      unname(z$aggregations[[agg_type]]$estimates),
      tolerance = 1e-14,
      info = agg_type
    )
    expect_equal(
      unname(agg$aggregation$if_full),
      unname(z$aggregations[[agg_type]]$if_full),
      tolerance = 1e-14,
      info = agg_type
    )
  }
})

test_that("CP02 installed public rotating80 pipeline matches built-in IPW with custom est_method and ... forwarding", {
  skip_if_not_installed("cdid")
  skip_if_not_installed("did")
  ref_env <- cp02_load_reference_helpers()
  z <- source(testthat::test_path("fixtures", "rotating80_fixture.R"))$value
  d <- z$panel
  d$g[is.infinite(d$g)] <- 0
  d$y <- d$y_exactfit

  seen_extra <- integer()
  custom_cb <- function(y1, y0, D, covariates, i.weights, inffunc = TRUE, tag_val = 0L) {
    seen_extra <<- c(seen_extra, as.integer(tag_val))
    ref_env$custom_ipw(y1, y0, D, covariates, i.weights, inffunc)
  }

  seed_before <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  # Preserve the original unit-weight/varying case; add universal reporting with
  # nonuniform weights, permuted full frame, non-A rows and inference-only clusters.
  frame <- d[!duplicated(d$id), c("id", "g", "x")]
  frame$sample_weight <- 1 + seq_len(nrow(frame)) %% 5 / 4
  frame$cluster <- (seq_len(nrow(frame)) - 1L) %/% 4L
  extra <- frame[1:2, ]; extra$id <- c(-101L, -102L)
  frame <- rbind(extra, frame[nrow(frame):1L, ])
  for (bp in c("varying", "universal")) {
    extra_args <- if (bp == "varying") list() else list(frame = frame,
      weightsname = "sample_weight", clustervars = "cluster")
    common <- c(list(yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~x,
      data = d, base_period = bp, bstrap = FALSE, cband = bp == "universal"), extra_args)
    mp_ipw <- do.call(cdid::att_gt, c(common, list(est_method = "ipw")))
    mp_cus <- do.call(cdid::att_gt, c(common, list(est_method = custom_cb, tag_val = 42L)))
    expect_gt(length(seen_extra), 0L)
    expect_true(all(seen_extra == 42L))
    expect_identical(mp_cus$DIDparams$est_method, custom_cb)
    expect_identical(mp_cus$DIDparams$extra_args, list(tag_val = 42L))
    for (field in c("att", "se", "inffunc", "V_analytical", "W", "Wpval")) {
      expect_equal(unname(mp_cus[[field]]), unname(mp_ipw[[field]]), tolerance = 1e-10, info = paste(bp, field))
    }
    for (field in c("status", "declared_cells", "n_declared_pretrend", "n_available_pretrend", "cells")) {
      expect_identical(mp_cus$chained$pretest[[field]], mp_ipw$chained$pretest[[field]])
    }
    expect_identical(mp_cus$chained$reduction_summary, mp_ipw$chained$reduction_summary)
    if (bp == "universal") {
      expect_identical(rownames(mp_cus$inffunc), as.character(frame$id))
      expect_true(all(mp_cus$inffunc[1:2, ] == 0))
    }
    cb_count_before <- length(seen_extra)
    for (agg_type in c("simple", "group", "dynamic", "calendar")) {
      a_ipw <- cdid::aggte(mp_ipw, type = agg_type, bstrap = FALSE, cband = FALSE)
      a_cus <- cdid::aggte(mp_cus, type = agg_type, bstrap = FALSE, cband = FALSE)
      for (field in c("overall.att", "overall.se", "egt", "att.egt", "se.egt", "inf.function")) {
        expect_equal(a_cus[[field]], a_ipw[[field]], tolerance = 1e-10, info = paste(bp, agg_type, field))
      }
    }
    expect_identical(length(seen_extra), cb_count_before)
  }
  seed_after <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  expect_identical(seed_after, seed_before)
})

test_that("original-IPW callbacks preserve balanced-panel point, design and complete-IF forest identities", {
  skip_if_not_installed("cdid")
  ref <- cp02_load_reference_helpers()
  p <- utils::read.csv(testthat::test_path("fixtures", "panel24.csv"))
  p <- p[p$id %in% 21:24 & p$time %in% c(1, 2, 4), ]
  expect_false(anyNA(p$y))
  own <- ipw_owner_control(3, 2e6, 2e6, 5e8, 5e5)
  pt <- ipw_point_view_control(3, 2000, 2e6, 2e6, 5e8, 5e5)
  po <- construct_ipw_panel_owner(p, "id", "time", "y", "g", own)
  fo <- construct_ipw_frame_owner(po, 21:24, matrix(1, 4, 1), rep(1, 4), own)
  meta <- build_owner_comparison_catalog(po, 3, 0, "nevertreated", "fixed-baseline-reference",
    control = own, frame_owner = fo, build_design = TRUE,
    declared_qualified = data.frame(g = 3, t = 4), target_post_times = 4,
    catalog_bounds = list(max_candidates = 1e5, max_rows = 1e5,
      max_pool_entries = 1e6, max_scratch = 1e7, max_work = 1e7))
  links <- build_custom_ipw_catalog_links(meta, fo, own, pt, ref$custom_ipw)
  cert <- certify_ipw_structural_classes_view(links$point_links, meta$design,
    match(seq_along(links$point_links), meta$design$links$catalog_row), "optimal",
    list(max_links = 500, max_classes = 50, max_class_dates = 50,
         max_E_entries = 5000, max_mask_entries = 100000), own)
  expect_gt(length(cert$dropped_index), 0L)
  phi <- matrix(0, 4, length(links$point_links))
  replay_custom_ipw_if_view(links$point_links, pt,
    callback = function(b) phi[b$a_index, ] <<- b$phi)
  mu <- vapply(links$point_links, `[[`, numeric(1), "estimate")
  W <- materialize_chained_design(meta$design, 1e5)$W
  keep <- cert$retained_index
  expect_equal(unname(apply_ipw_reduction_map(cert, mu[keep])), unname(mu), tolerance = 1e-10)
  expect_equal(unname(t(apply_ipw_reduction_map(cert, t(phi[, keep, drop = FALSE])))),
    unname(phi), tolerance = 1e-10)
  expect_equal(unname(apply_ipw_reduction_map(cert, W[keep, , drop = FALSE])), unname(W), tolerance = 1e-10)
  # Checked whole-catalog simultaneous retention refuses before the first callback.
  before <- ref$counts$callback_calls
  expect_error(build_custom_ipw_catalog_links(meta, fo, own, pt, ref$custom_ipw,
    retained_pipeline_entries = pt$max_matrix_entries), class = "cdid_ipw_owner_bound_error")
  expect_identical(ref$counts$callback_calls, before)
})
