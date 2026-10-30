# Private CP07 Chained Doubly Robust (DR) backend tests (DR_SPEC.md, TASK.txt).
# Uses ONLY locked installed cdid namespace functions and real sealed owners
# (panel, frame, reader, link views) constructed from the frozen independent
# reference CSVs (tests/testthat/fixtures/dr_finite_v1).

dr_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, envir = ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}

dr_fx_dir <- testthat::test_path("fixtures", "dr_finite_v1")
dr_loader <- new.env(parent = globalenv())
sys.source(file.path(dr_fx_dir, "load_fixture.R"), envir = dr_loader)

dr_specs <- dr_loader$dr_fixture_csv(dr_fx_dir, "fixture_specs.csv")
dr_quals <- dr_loader$dr_fixture_csv(dr_fx_dir, "qualifications.csv")
dr_pop <- dr_loader$dr_fixture_csv(dr_fx_dir, "population_examples.csv")
dr_perm_order <- c(seq.int(36L, 2L, -2L), seq.int(35L, 1L, -2L))
dr_mx <- function(x) max(c(0, abs(x)))

dr_load <- function(fixture_id, order = NULL, max_chunk_rows = 13L) {
  dr_loader$load_dr_finite_reference(
    dr_fx_dir,
    fixture_id = fixture_id,
    order = order,
    max_chunk_rows = max_chunk_rows
  )
}

dr_fit <- function(fx,
                   views = fx$views,
                   control = fx$control,
                   glm_control = fx$glm_control,
                   rank_tol = fx$rank_tol,
                   outcome_recipe = fx$outcome_recipe,
                   treatment_columns = fx$treatment_columns,
                   selection = fx$selection,
                   selection_model = fx$selection_model,
                   selection_pooling = fx$selection_pooling,
                   sel_fit = NULL) {
  dr_ns("fit_dr_backend")(
    reader = fx$reader,
    link_views = views,
    control = control,
    glm_control = glm_control,
    rank_tol = rank_tol,
    outcome_recipe = outcome_recipe,
    treatment_columns = treatment_columns,
    selection = selection,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    sel_fit = sel_fit
  )
}

dr_collect <- function(fit, chunk_rows = 11L, link_batch_size = 2L,
                       scale = "frame", details = TRUE) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  N <- st$N
  L <- length(fit$links)
  p <- fit$shared$p
  G <- length(fit$regressions)
  J <- length(st$cohorts)
  na_mat <- function(k) matrix(NA_real_, N, k)

  out <- list(
    phi = na_mat(L),
    psi_dir = na_mat(L),
    psi_or = na_mat(L),
    psi_ps = na_mat(L),
    psi_sel = na_mat(L),
    psi_ipw = na_mat(L),
    psi_or_estimator = na_mat(L),
    beta = na_mat(p),
    beta_score = na_mat(p),
    block_scores = lapply(seq_len(G), function(b) na_mat(p)),
    alpha = na_mat(G),
    tau = na_mat(J),
    eta = na_mat(J),
    lambda = rep(NA_real_, N)
  )

  out$info <- dr_ns("replay_dr_if")(
    fit,
    chunk_rows = chunk_rows,
    link_batch_size = link_batch_size,
    callback = function(b) {
      stopifnot(identical(b$id, st$id[b$rows]))
      out$phi[b$rows, b$links] <<- b$phi
      if (!details) return(invisible(NULL))
      d <- b$details
      out$psi_dir[b$rows, b$links] <<- d$psi_dir
      out$psi_or[b$rows, b$links] <<- d$psi_or
      out$psi_ps[b$rows, b$links] <<- d$psi_ps
      out$psi_sel[b$rows, b$links] <<- d$psi_sel
      out$psi_ipw[b$rows, b$links] <<- d$psi_ipw
      out$psi_or_estimator[b$rows, b$links] <<- d$psi_or_estimator
      if (p > 0L) {
        out$beta[b$rows, ] <<- d$beta
        out$beta_score[b$rows, ] <<- d$beta_score
        for (k in seq_along(d$beta_regression_ids)) {
          gid <- d$beta_regression_ids[k]
          out$block_scores[[gid]][b$rows, ] <<- d$beta_block_scores[[k]]
        }
      }
      out$alpha[b$rows, d$regression_ids] <<- d$alpha
      out$tau[b$rows, ] <<- d$tau
      out$eta[b$rows, ] <<- d$eta
      out$lambda[b$rows] <<- d$lambda
    },
    scale = scale,
    details = details
  )
  out
}

dr_fit_table <- function(fx, pair_outcomes = NULL) {
  if (fx$selection == "none") {
    dr_ns("fit_dr_baseline_table")(
      X_or = fx$X_or,
      V = fx$V_trt,
      Y = fx$Y_obs,
      obs = fx$R,
      g = fx$raw$cohort,
      w = fx$raw$w,
      calendar = fx$catalog$calendar$time,
      links = fx$catalog$link_list,
      outcome_recipe = fx$outcome_recipe,
      pair_outcomes = pair_outcomes,
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    )
  } else {
    dr_ns("fit_dr_selection_table")(
      Z = fx$Q,
      X_or = fx$X_or,
      V = fx$V_trt,
      Y = fx$Y_obs,
      obs = fx$R,
      g = fx$raw$cohort,
      w = fx$raw$w,
      U = fx$U,
      rho = fx$rho,
      universe = fx$spec$U,
      calendar = fx$catalog$calendar$time,
      links = fx$catalog$link_list,
      outcome_recipe = fx$outcome_recipe,
      pair_outcomes = pair_outcomes,
      selection_model = fx$selection_model,
      selection_pooling = fx$selection_pooling,
      has_intercept = fx$has_intercept,
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol
    )
  }
}

# Map backend regression blocks to frozen reference block IDs (b1, b2, b3).
dr_block_ids <- function(fx, fit) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  bl <- fx$catalog$blocks
  bl_ctl <- vapply(
    strsplit(bl$K, "|", fixed = TRUE),
    function(x) paste(sort(as.double(x)), collapse = ","),
    character(1)
  )
  vapply(fit$regression_specs, function(sp) {
    pr <- fit$support_specs[[sp$support_ids[1L]]]$pair
    ctl <- paste(sort(as.double(st$cohorts[sp$cohorts])), collapse = ",")
    idx <- which(bl$s == pr$s & bl$t == pr$t & bl_ctl == ctl)
    stopifnot(length(idx) == 1L)
    bl$block_id[idx]
  }, character(1))
}

# Complete-unit-first cluster covariance helper: sums complete unit IF rows
# within each cluster before forming outer products.
dr_cluster_cov <- function(phi, cluster_id) {
  N <- nrow(phi)
  L <- ncol(phi)
  clusters <- sort(unique(cluster_id))
  G_cl <- matrix(0, length(clusters), L)
  for (i in seq_along(clusters)) {
    rows <- which(cluster_id == clusters[i])
    G_cl[i, ] <- colSums(phi[rows, , drop = FALSE])
  }
  Omega_cl <- crossprod(G_cl) / N
  V_cl <- Omega_cl / N
  list(V_cluster = V_cl, Omega_cluster = Omega_cl)
}

# Independent accepted functional evaluated at the owner's actual masked-Y
# differences. Original stored-Delta reference targets remain separate below.
dr_expect_actual_y <- function(fx, fit, out) {
  all_values <- dr_loader$dr_fixture_csv(dr_fx_dir, "actual_y_values.csv")
  ref <- all_values[all_values$fixture_id == fx$fixture_id, , drop = FALSE]
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  ids <- fx$raw$unit_id
  links <- vapply(fit$links, `[[`, character(1), "link_id")
  beta_cols <- if (fit$shared$p > 0L) paste0("beta_", fit$shared$basis$kept_columns) else character(0)
  blocks <- dr_block_ids(fx, fit)
  cluster <- dr_cluster_cov(out$phi, fx$raw$cluster_id)
  matrices <- list(V_unit = crossprod(out$phi) / st$N^2,
                   V_cluster = cluster$V_cluster, Omega_cluster = cluster$Omega_cluster)
  components <- c(psi_dr = "phi", psi_dir = "psi_dir", psi_or = "psi_or",
                  psi_ps = "psi_ps", psi_sel = "psi_sel", psi_ipw = "psi_ipw",
                  beta_if = "beta", tau_if = "tau", eta_if = "eta")
  actual <- vapply(seq_len(nrow(ref)), function(i) {
    row <- ref[i, , drop = FALSE]
    if (row$group %in% c("mu_dr", "mu_or", "mu_ipw")) {
      return(fit$links[[match(row$key1, links)]][[row$group]])
    }
    if (row$group == "beta_hat") return(unname(fit$shared$beta[match(row$key1, beta_cols)]))
    if (row$group == "alpha_hat") return(fit$regressions[[match(row$key1, blocks)]]$alpha)
    if (row$group %in% names(components)) {
      cols <- if (row$group == "beta_if") beta_cols else if (row$group == "tau_if") paste0("tau_", st$cohorts) else
        if (row$group == "eta_if") paste0("eta_", st$cohorts) else links
      return(out[[components[[row$group]]]][match(as.integer(row$key1), ids), match(row$key2, cols)])
    }
    matrices[[row$group]][match(row$key1, links), match(row$key2, links)]
  }, numeric(1))
  testthat::expect_true(nrow(ref) > 0L && all(is.finite(actual)))
  for (group in unique(ref$group)) {
    error <- max(abs(actual[ref$group == group] - ref$value[ref$group == group]))
    limit <- if (group %in% c("V_unit", "V_cluster", "Omega_cluster")) 1e-15 else
      if (group %in% c("tau_if", "eta_if")) 1e-14 else 2e-14
    if (group == "Omega_cluster") error <- error / st$N
    testthat::expect_lt(error, limit)
  }
  nrow(ref)
}

testthat::test_that("CSV covariance expectations retain scientific link coordinates when reordered", {
  atoms <- dr_loader$dr_fixture_csv(dr_fx_dir, "atoms.csv")
  links <- dr_loader$dr_fixture_csv(dr_fx_dir, "catalog_links.csv")
  reversed <- rev(links$link_id)
  for (fid in dr_specs$fixture_id) {
    original <- dr_loader$dr_load_expectations(dr_fx_dir, fid, atoms$unit_id, links$link_id)
    reordered <- dr_loader$dr_load_expectations(dr_fx_dir, fid, atoms$unit_id, reversed)
    for (name in c("V_unit", "V_cluster", "Omega_cluster")) {
      testthat::expect_identical(reordered[[name]], original[[name]][reversed, reversed, drop = FALSE])
    }
  }
})

testthat::test_that("F1 baseline DR point estimates, 4-term IF rows, nuisances, and covariances match the frozen reference", {
  fx <- dr_load("F1_baseline_none_intercept_slope")
  testthat::expect_identical(fx$panel$n_a_units, 24L)
  testthat::expect_identical(fx$reader$state$N, 36L)
  testthat::expect_identical(fx$reader$state$n_A, 24L)
  testthat::expect_true(all(c(10, 22, 34) %in% fx$panel$a_id))

  fit <- dr_fit(fx)
  dr_ns("validate_dr_baseline_backend")(fit)
  dr_ns("validate_dr_backend")(fit)

  out_frame <- dr_collect(fit, chunk_rows = 11L, link_batch_size = 2L, scale = "frame", details = TRUE)
  out_core <- dr_collect(fit, chunk_rows = 13L, link_batch_size = 4L, scale = "core", details = FALSE)
  testthat::expect_identical(out_frame$info$normalization_count, 36L)
  testthat::expect_identical(out_core$info$normalization_count, 24L)
  testthat::expect_lt(dr_mx(out_core$phi - (24 / 36) * out_frame$phi), 1e-15)

  # Exact machine-precision match between real owner pipeline and table helper on Y_obs
  tb_y <- dr_fit_table(fx)
  testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_dr") - tb_y$estimates), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$phi - tb_y$phi), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$psi_dir - tb_y$psi_dir), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$psi_or - tb_y$psi_or), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$psi_ps - tb_y$psi_ps), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$psi_sel - tb_y$psi_sel), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$psi_ipw - tb_y$psi_ipw), 1e-14)

  # Exact machine-precision match against frozen reference on observed Delta_e
  tb_d <- dr_fit_table(fx, pair_outcomes = fx$pair_outcomes_delta)
  exp <- fx$expectations
  testthat::expect_lt(dr_mx(tb_d$estimates - exp$mu_dr), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$estimates_or - exp$mu_or), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$estimates_ipw - exp$mu_ipw), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$shared$beta - exp$beta_hat), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$phi - exp$psi_dr), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$phi_core - exp$psi_dr_core), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$psi_dir - exp$psi_dir), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$psi_or - exp$psi_or), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$psi_ps - exp$psi_ps), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$psi_sel - exp$psi_sel), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$psi_ipw - exp$psi_ipw), 1e-14)
  testthat::expect_lt(dr_mx(tb_d$beta_if - exp$beta_if), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$tau - exp$tau_if), 1e-14)
  testthat::expect_lt(dr_mx(out_frame$eta - exp$eta_if), 1e-14)

  # Direct comparison of real owner fit (using Y_t - Y_s from 8-decimal Y_obs) against reference
  bid <- dr_block_ids(fx, fit)
  testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_dr") - exp$mu_dr), 5e-9)
  testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_or") - exp$mu_or), 5e-9)
  testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_ipw") - exp$mu_ipw), 5e-9)
  testthat::expect_lt(dr_mx(fit$shared$beta - exp$beta_hat), 5e-9)
  testthat::expect_lt(dr_mx(vapply(fit$regressions, `[[`, numeric(1), "alpha") - exp$alpha_hat[bid]), 5e-9)
  testthat::expect_lt(dr_mx(out_frame$phi - exp$psi_dr), 1.5e-7)
})

testthat::test_that("all 8 frozen fixtures F1-F8 match points, 1152 unit-link IF rows, 2088 nuisance rows, and covariances", {
  total_uif_rows <- 0L
  total_nif_rows <- 0L
  total_actual_y_values <- 0L

  for (fid in dr_specs$fixture_id) {
    fx <- dr_load(fid, order = dr_perm_order, max_chunk_rows = 11L)
    fit <- dr_fit(fx)
    dr_ns("validate_dr_backend")(fit)

    out <- dr_collect(fit, chunk_rows = 9L, link_batch_size = 2L, scale = "frame", details = TRUE)
    out_core <- dr_collect(fit, chunk_rows = 11L, link_batch_size = 3L, scale = "core", details = FALSE)
    exp <- fx$expectations
    N <- fx$reader$state$N
    n_A <- fx$reader$state$n_A
    L <- length(fit$links)
    p <- fit$shared$p
    bid <- dr_block_ids(fx, fit)
    total_actual_y_values <- total_actual_y_values + dr_expect_actual_y(fx, fit, out)

    total_uif_rows <- total_uif_rows + nrow(out$phi) * ncol(out$phi)
    total_nif_rows <- total_nif_rows + nrow(out$tau) * ncol(out$tau) +
      nrow(out$eta) * ncol(out$eta) + nrow(out$beta) * ncol(out$beta)

    # 1. Normalization & 4-term decomposition identities
    testthat::expect_identical(out$info$normalization_count, 36L)
    testthat::expect_identical(out_core$info$normalization_count, 24L)
    testthat::expect_lt(dr_mx(out_core$phi - (n_A / N) * out$phi), 1e-15)
    testthat::expect_lt(
      dr_mx(out$phi - (out$psi_dir + out$psi_or + out$psi_ps + out$psi_sel)),
      1e-14
    )
    testthat::expect_lt(dr_mx(colMeans(out$phi)), 1e-13)
    for (ell in seq_len(L)) {
      lk <- fit$links[[ell]]
      testthat::expect_lt(abs(lk$mu_dr - lk$mu_dr_decomp), 1e-14)
      testthat::expect_lt(abs(lk$d_plus - exp$link_summaries$d_plus[ell]), 1e-14)
      testthat::expect_lt(abs(lk$d_minus - exp$link_summaries$d_minus[ell]), 1e-14)
    }

    # 2. Global beta score is the sum over ALL 3 deduplicated regression blocks
    if (p > 0L) {
      sum_blocks <- Reduce(`+`, out$block_scores)
      testthat::expect_lt(dr_mx(out$beta_score - sum_blocks), 1e-14)
      testthat::expect_identical(out$info$batches, as.integer(ceiling(N / 9L) * ceiling(L / 2L)))
    }

    # 3. Off-A support mask discipline across routes and universes
    off_A <- which(fx$raw$A == 0L)
    b_minus_a <- which(fx$raw$B == 1L & fx$raw$A == 0L)
    f_minus_b <- which(fx$raw$F == 1L & fx$raw$B == 0L)
    if (fx$selection == "none" || fx$spec$U == "A") {
      testthat::expect_equal(dr_mx(out$phi[off_A, , drop = FALSE]), 0)
    } else if (fx$spec$U == "B") {
      testthat::expect_gt(dr_mx(out$phi[b_minus_a, , drop = FALSE]), 1e-3)
      testthat::expect_equal(dr_mx(out$phi[f_minus_b, , drop = FALSE]), 0)
    } else {
      testthat::expect_gt(dr_mx(out$phi[b_minus_a, , drop = FALSE]), 1e-3)
      testthat::expect_gt(dr_mx(out$phi[f_minus_b, , drop = FALSE]), 1e-3)
    }

    # 4. Exact match between real owner pipeline and table helper on Y_obs
    tb_y <- dr_fit_table(fx)
    testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_dr") - tb_y$estimates), 2e-14)
    testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_or") - tb_y$estimates_or), 2e-14)
    testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_ipw") - tb_y$estimates_ipw), 2e-14)
    testthat::expect_lt(dr_mx(out$phi - tb_y$phi), 2e-14)
    testthat::expect_lt(dr_mx(out$psi_dir - tb_y$psi_dir), 2e-14)
    testthat::expect_lt(dr_mx(out$psi_or - tb_y$psi_or), 2e-14)
    testthat::expect_lt(dr_mx(out$psi_ps - tb_y$psi_ps), 2e-14)
    testthat::expect_lt(dr_mx(out$psi_sel - tb_y$psi_sel), 2e-14)
    testthat::expect_lt(dr_mx(out$psi_ipw - tb_y$psi_ipw), 2e-14)

    # 5. Exact machine-precision match against frozen reference on observed Delta_e
    tb_d <- dr_fit_table(fx, pair_outcomes = fx$pair_outcomes_delta)
    testthat::expect_lt(dr_mx(tb_d$estimates - exp$mu_dr), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$estimates_or - exp$mu_or), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$estimates_ipw - exp$mu_ipw), 2e-14)
    if (p > 0L) {
      testthat::expect_lt(dr_mx(tb_d$shared$beta - exp$beta_hat), 2e-14)
      testthat::expect_lt(dr_mx(tb_d$beta_if - exp$beta_if), 2e-14)
    }
    if (fx$outcome_recipe$intercept) {
      testthat::expect_lt(
        dr_mx(vapply(tb_d$regressions, `[[`, numeric(1), "alpha") - exp$alpha_hat[bid]),
        2e-14
      )
    }
    testthat::expect_lt(dr_mx(tb_d$phi - exp$psi_dr), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$phi_core - exp$psi_dr_core), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$psi_dir - exp$psi_dir), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$psi_or - exp$psi_or), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$psi_ps - exp$psi_ps), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$psi_sel - exp$psi_sel), 2e-14)
    testthat::expect_lt(dr_mx(tb_d$psi_ipw - exp$psi_ipw), 2e-14)
    testthat::expect_lt(dr_mx(out$tau - exp$tau_if), 1e-14)
    testthat::expect_lt(dr_mx(out$eta - exp$eta_if), 1e-14)

    # 6. Direct comparison of real owner fit (Y_obs) against reference
    testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_dr") - exp$mu_dr), 5e-9)
    testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_or") - exp$mu_or), 5e-9)
    testthat::expect_lt(dr_mx(vapply(fit$links, `[[`, numeric(1), "mu_ipw") - exp$mu_ipw), 5e-9)
    testthat::expect_lt(dr_mx(out$phi - exp$psi_dr), 1.5e-7)

    # 7. Streamed unit covariance and complete-unit-first cluster covariance
    cov_stream <- dr_ns("stream_dr_covariance")(fit)
    V_unit_owner <- cov_stream$covariance / n_A
    V_unit_phi <- crossprod(out$phi) / (N^2)
    testthat::expect_lt(dr_mx(V_unit_owner - V_unit_phi), 1e-15)
    testthat::expect_lt(dr_mx(crossprod(tb_d$phi) / (N^2) - exp$V_unit), 1e-15)
    testthat::expect_lt(dr_mx(V_unit_owner - exp$V_unit), 5e-9)

    cl_d <- dr_cluster_cov(tb_d$phi, fx$raw$cluster_id)
    testthat::expect_lt(dr_mx(cl_d$V_cluster - exp$V_cluster), 1e-15)
    # Compare the same covariance on its existing V scale; retain Omega = N V.
    testthat::expect_lt(dr_mx((cl_d$Omega_cluster - exp$Omega_cluster) / N), 1e-15)
    testthat::expect_lt(dr_mx(cl_d$Omega_cluster / N - cl_d$V_cluster), 1e-15)
    testthat::expect_lt(dr_mx(exp$Omega_cluster / N - exp$V_cluster), 1e-15)
    testthat::expect_gt(dr_mx(cl_d$Omega_cluster - N * exp$V_unit), 0.05)
  }

  testthat::expect_identical(total_uif_rows, 1152L)
  testthat::expect_identical(total_nif_rows, 2088L)
  testthat::expect_identical(total_actual_y_values, 9514L)
})

testthat::test_that("intercept-only DR equals IPW identically across Routes 1, 2, and 3 including IF rows and covariances", {
  # Route 1: F2_baseline_none_intercept_only
  fx2 <- dr_load("F2_baseline_none_intercept_only")
  fit2 <- dr_fit(fx2)
  out2 <- dr_collect(fit2, chunk_rows = 12L, link_batch_size = 2L, scale = "frame", details = TRUE)
  cov2 <- dr_ns("stream_dr_covariance")(fit2)
  testthat::expect_identical(fit2$shared$p, 0L)
  testthat::expect_lt(
    dr_mx(vapply(fit2$links, `[[`, numeric(1), "mu_dr") - vapply(fit2$links, `[[`, numeric(1), "mu_ipw")),
    1e-15
  )
  testthat::expect_equal(dr_mx(out2$psi_or), 0)
  testthat::expect_lt(dr_mx(out2$phi - out2$psi_ipw), 3e-15)

  # Route 2 (pair_logit) and Route 3 (independent_period_logit) intercept-only DR == IPW parity
  for (fid in c("F3_pair_logit_F_separate", "F7_date_product_F_within_cohort")) {
    fx <- dr_load(fid)
    icpt_recipe <- dr_ns("or_pooled_outcome_recipe")("(Intercept)", TRUE, character(0))
    fit_dr <- dr_fit(fx, outcome_recipe = icpt_recipe)
    out_dr <- dr_collect(fit_dr, chunk_rows = 13L, link_batch_size = 2L, scale = "frame", details = TRUE)

    # Also fit standalone CP04 selection IPW consumer on the exact same precomputed sel_fit
    fit_ipw <- dr_ns("fit_cp04_selection_models")(
      fx$reader,
      fx$views,
      control = fx$control,
      glm_control = fx$glm_control,
      selection_model = fx$selection_model,
      selection_pooling = fx$selection_pooling,
      sel_fit = fit_dr$selection_backend
    )
    ipw_pts <- vapply(seq_along(fx$views), function(j) dr_ns("chained_full_link_estimate")(fit_ipw, j), numeric(1))
    N <- fx$reader$state$N
    ipw_phi <- matrix(NA_real_, N, length(fx$views))
    for (s in seq.int(1L, N, by = 11L)) {
      rr <- seq.int(s, min(N, s + 10L))
      ipw_phi[rr, ] <- dr_ns("evaluate_cp04_selection_if_chunk")(fit_ipw, rr, seq_along(fx$views), scale = "frame")$phi
    }

    testthat::expect_lt(dr_mx(vapply(fit_dr$links, `[[`, numeric(1), "mu_dr") - ipw_pts), 1e-15)
    testthat::expect_equal(dr_mx(out_dr$psi_or), 0)
    # Existing deterministic IF convention: selection-public.R:686-703.
    testthat::expect_lt(dr_mx(out_dr$phi - out_dr$psi_ipw), 1e-13)
    testthat::expect_lt(dr_mx(out_dr$phi - ipw_phi), 1e-13)

    cov_dr <- dr_ns("stream_dr_covariance")(fit_dr)
    cov_ipw <- dr_ns("cp04_stream_covariance")(fit_ipw, seq_along(fx$views))
    testthat::expect_lt(dr_mx(cov_dr$covariance - cov_ipw$covariance), 1e-14)
  }
})

testthat::test_that("generic 4-term nonzero witnesses, saturated binary-x1 identity, and float64/binary128 qualifications hold", {
  # 1. Verify retained float64 eps=1e-6 failure vs binary128 quad reference in qualifications.csv
  sm <- dr_quals[dr_quals$metric_group == "summary_metrics", ]
  q_val <- function(nm) as.double(sm$value[sm$metric_name == nm])
  testthat::expect_gt(q_val("max_dr_if_fd_err_1e6_float64"), 1e-8)
  testthat::expect_lt(abs(q_val("max_dr_if_fd_err_1e6_float64") - 1.1957433265763484e-8), 1e-20)
  testthat::expect_lt(q_val("max_dr_if_fd_err_1e6_quad"), 1e-8)
  testthat::expect_lt(abs(q_val("max_dr_if_fd_err_1e6_quad") - 1.2576606422953773e-12), 1e-24)
  testthat::expect_lt(q_val("max_dr_if_fd_err_1e5_float64"), 2e-9)

  # 2. Verify generic 4-term nonzero witnesses on F5, F7, F8
  for (fid in c("F5_pair_logit_A_across_no_intercept_or",
                "F7_date_product_F_within_cohort",
                "F8_date_product_F_across_no_intercept_sel")) {
    fx <- dr_load(fid)
    tb_d <- dr_fit_table(fx, pair_outcomes = fx$pair_outcomes_delta)
    w_sub <- dr_quals[dr_quals$metric_group == "four_term_nonzero_witnesses" & dr_quals$fixture_id == fid, ]
    w_val <- function(nm) as.double(w_sub$value[w_sub$metric_name == nm])

    testthat::expect_lt(abs(dr_mx(tb_d$psi_dir) - w_val("max_over_links_dir")), 1e-13)
    testthat::expect_lt(abs(dr_mx(tb_d$psi_or) - w_val("max_over_links_or")), 1e-13)
    testthat::expect_lt(abs(dr_mx(tb_d$psi_ps) - w_val("max_over_links_ps")), 1e-13)
    testthat::expect_lt(abs(dr_mx(tb_d$psi_sel) - w_val("max_over_links_sel")), 1e-13)

    # Omitting any of the 4 terms fails against the binary128 FD derivative
    fd_quad <- fx$expectations$fd_1e6_quad
    testthat::expect_lt(dr_mx(tb_d$phi - fd_quad), 2e-12)
    testthat::expect_gt(dr_mx((tb_d$phi - tb_d$psi_dir) - fd_quad), 3.0)
    testthat::expect_gt(dr_mx((tb_d$phi - tb_d$psi_or) - fd_quad), 1e-3)
    testthat::expect_gt(dr_mx((tb_d$phi - tb_d$psi_ps) - fd_quad), 1.0)
    testthat::expect_gt(dr_mx((tb_d$phi - tb_d$psi_sel) - fd_quad), 0.3)
  }

  # 3. Verify retained failed parent prediction 1 on saturated binary-x1 designs (F1, F3, F4, F6)
  pred_sub <- dr_quals[dr_quals$metric_group == "retained_failed_predictions", ]
  testthat::expect_identical(
    pred_sub$value[pred_sub$metric_name == "PARENT_PLAN_PRED_1_SATURATED_BINARY_X1_OR_LOADING:original_predicate_status"],
    "FAIL"
  )
  testthat::expect_identical(
    pred_sub$value[pred_sub$metric_name == "PARENT_PLAN_PRED_1_SATURATED_BINARY_X1_OR_LOADING:corrected_identity_status"],
    "PASS"
  )
  for (fid in c("F1_baseline_none_intercept_slope",
                "F3_pair_logit_F_separate",
                "F4_pair_logit_B_within_cohort",
                "F6_date_product_F_separate")) {
    fx <- dr_load(fid)
    fit <- dr_fit(fx)
    out <- dr_collect(fit, chunk_rows = 12L, link_batch_size = 2L, scale = "frame", details = TRUE)
    d_diffs <- vapply(fit$links, function(lk) dr_mx(lk$d_bar_diff), numeric(1))
    testthat::expect_lt(max(d_diffs), 1e-14)
    testthat::expect_lt(dr_mx(out$psi_or), 1e-15)
    testthat::expect_gt(dr_mx(out$beta), 0.1)
  }

  # 4. Verify population examples CSV (P1-P6) including P3B retained failed prediction 2
  p3b <- dr_pop[dr_pop$case_id == "P3B_PS_ONLY_EQUAL_LAWS_INTERCEPT_CONSISTENCY", ]
  p3a <- dr_pop[dr_pop$case_id == "P3A_PS_ONLY_DISTINCT_LAWS", ]
  testthat::expect_lt(abs(as.double(p3b$or_bias)), 1e-14)
  testthat::expect_gt(abs(as.double(p3a$or_bias)), 0.19)
  testthat::expect_lt(abs(as.double(p3a$dr_bias)), 1e-14)
})

testthat::test_that("no-intercept, exact alias, and zero-slope OR branches work and preserve invariants", {
  # 1. Exact alias in OR slopes on both baseline (F1) and conditional_mtar (F7)
  for (fid in c("F1_baseline_none_intercept_slope", "F7_date_product_F_within_cohort")) {
    fx <- dr_load(fid)
    fit_base <- dr_fit(fx)
    out_base <- dr_collect(fit_base, chunk_rows = 12L, link_batch_size = 2L, scale = "frame", details = TRUE)
    cov_base <- dr_ns("stream_dr_covariance")(fit_base)

    # Augment frame$Z_a with exact aliased columns x1_dup = x1 and x_comb = 2 * x1
    Z_aug <- cbind(
      `(Intercept)` = 1,
      x1 = fx$raw$x1,
      x2 = fx$raw$x2,
      x1_dup = fx$raw$x1,
      x_comb = 2 * fx$raw$x1
    )
    Z_aug[fx$raw$A == 0L, ] <- NA_real_
    frame_aug <- dr_ns("construct_ipw_frame_owner")(
      fx$panel,
      as.double(fx$raw$unit_id),
      Z_aug,
      weights = as.double(fx$raw$w),
      control = fx$control,
      frame_g = as.double(fx$raw$cohort),
      allow_frame_only = TRUE,
      design_on_A_only = TRUE
    )
    just <- if (fx$spec$U %in% c("B", "A")) "known transport under the frozen independent fixed atom law" else NULL
    reader_aug <- dr_ns("construct_cp04_frame_reader")(
      frame_aug,
      as.double(fx$raw$unit_id),
      fx$Q,
      as.double(fx$raw$cohort),
      as.double(fx$raw$w),
      fx$control,
      U = fx$U,
      rho = fx$rho,
      universe = fx$spec$U,
      justification = just,
      has_intercept = fx$has_intercept
    )
    views_aug <- lapply(seq_len(nrow(fx$catalog$links)), function(i) {
      lk <- fx$catalog$links[i, , drop = FALSE]
      dr_ns("construct_ipw_link_view")(
        frame_aug,
        as.double(lk$s),
        as.double(lk$t),
        g = as.double(lk$g),
        K = as.double(dr_loader$dr_split_pipe(lk$K)),
        control = fx$control,
        link_id = lk$link_id,
        functional_key = lk$functional_key
      )
    })
    or_slopes_base <- fx$outcome_recipe$slope_columns
    or_slopes_aug <- c(or_slopes_base, "x1_dup", "x_comb")
    rec_aug <- dr_ns("or_pooled_outcome_recipe")(
      c("(Intercept)", or_slopes_aug),
      TRUE,
      or_slopes_aug
    )
    fit_aug <- dr_ns("fit_dr_backend")(
      reader_aug,
      views_aug,
      control = fx$control,
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol,
      outcome_recipe = rec_aug,
      treatment_columns = fx$treatment_columns,
      selection = fx$selection,
      selection_model = fx$selection_model,
      selection_pooling = fx$selection_pooling
    )
    out_aug <- dr_collect(fit_aug, chunk_rows = 12L, link_batch_size = 2L, scale = "frame", details = TRUE)
    cov_aug <- dr_ns("stream_dr_covariance")(fit_aug)

    testthat::expect_identical(fit_aug$shared$p, fit_base$shared$p)
    testthat::expect_identical(fit_aug$accounting$n_aliased_slopes, 2L)
    testthat::expect_identical(fit_aug$shared$basis$aliased_columns, c("x1_dup", "x_comb"))
    testthat::expect_lt(
      dr_mx(vapply(fit_aug$links, `[[`, numeric(1), "mu_dr") - vapply(fit_base$links, `[[`, numeric(1), "mu_dr")),
      1e-14
    )
    testthat::expect_lt(dr_mx(out_aug$phi - out_base$phi), 1e-14)
    testthat::expect_lt(dr_mx(out_aug$psi_dir - out_base$psi_dir), 1e-14)
    testthat::expect_lt(dr_mx(out_aug$psi_or - out_base$psi_or), 1e-14)
    testthat::expect_lt(dr_mx(out_aug$psi_ps - out_base$psi_ps), 1e-14)
    testthat::expect_lt(dr_mx(out_aug$psi_sel - out_base$psi_sel), 1e-14)
    testthat::expect_lt(dr_mx(cov_aug$covariance - cov_base$covariance), 1e-14)

    # Zero-slope column (x_zero = 0) aliased with intercept reduces to intercept-only DR == IPW
    Z_zero <- cbind(`(Intercept)` = 1, x1 = fx$raw$x1, x2 = fx$raw$x2, x_zero = 0)
    Z_zero[fx$raw$A == 0L, ] <- NA_real_
    frame_zero <- dr_ns("construct_ipw_frame_owner")(
      fx$panel,
      as.double(fx$raw$unit_id),
      Z_zero,
      weights = as.double(fx$raw$w),
      control = fx$control,
      frame_g = as.double(fx$raw$cohort),
      allow_frame_only = TRUE,
      design_on_A_only = TRUE
    )
    reader_zero <- dr_ns("construct_cp04_frame_reader")(
      frame_zero,
      as.double(fx$raw$unit_id),
      fx$Q,
      as.double(fx$raw$cohort),
      as.double(fx$raw$w),
      fx$control,
      U = fx$U,
      rho = fx$rho,
      universe = fx$spec$U,
      justification = just,
      has_intercept = fx$has_intercept
    )
    views_zero <- lapply(seq_len(nrow(fx$catalog$links)), function(i) {
      lk <- fx$catalog$links[i, , drop = FALSE]
      dr_ns("construct_ipw_link_view")(
        frame_zero,
        as.double(lk$s),
        as.double(lk$t),
        g = as.double(lk$g),
        K = as.double(dr_loader$dr_split_pipe(lk$K)),
        control = fx$control,
        link_id = lk$link_id,
        functional_key = lk$functional_key
      )
    })
    rec_zero <- dr_ns("or_pooled_outcome_recipe")(c("(Intercept)", "x_zero"), TRUE, "x_zero")
    fit_zero <- dr_ns("fit_dr_backend")(
      reader_zero,
      views_zero,
      control = fx$control,
      glm_control = fx$glm_control,
      rank_tol = fx$rank_tol,
      outcome_recipe = rec_zero,
      treatment_columns = fx$treatment_columns,
      selection = fx$selection,
      selection_model = fx$selection_model,
      selection_pooling = fx$selection_pooling
    )
    out_zero <- dr_collect(fit_zero, chunk_rows = 12L, link_batch_size = 2L, scale = "frame", details = TRUE)
    testthat::expect_identical(fit_zero$shared$p, 0L)
    testthat::expect_identical(fit_zero$accounting$n_aliased_slopes, 1L)
    testthat::expect_lt(
      dr_mx(vapply(fit_zero$links, `[[`, numeric(1), "mu_dr") - vapply(fit_zero$links, `[[`, numeric(1), "mu_ipw")),
      1e-15
    )
    testthat::expect_lt(dr_mx(out_zero$phi - out_zero$psi_ipw), 1e-13)
  }
})

testthat::test_that("case-weight scaling, unit ID permutation, and chunk/batch sizes leave points, IFs, and covariances invariant", {
  for (fid in c("F1_baseline_none_intercept_slope",
                "F4_pair_logit_B_within_cohort",
                "F5_pair_logit_A_across_no_intercept_or",
                "F7_date_product_F_within_cohort",
                "F8_date_product_F_across_no_intercept_sel")) {
    fx_id <- dr_load(fid, order = seq_len(36L), max_chunk_rows = 36L)
    fx_perm <- dr_load(fid, order = dr_perm_order, max_chunk_rows = 7L)

    fit_id <- dr_fit(fx_id)
    fit_perm <- dr_fit(fx_perm)

    out_id <- dr_collect(fit_id, chunk_rows = 36L, link_batch_size = 4L, scale = "frame", details = TRUE)
    out_perm_1 <- dr_collect(fit_perm, chunk_rows = 1L, link_batch_size = 1L, scale = "frame", details = TRUE)
    out_perm_5 <- dr_collect(fit_perm, chunk_rows = 5L, link_batch_size = 2L, scale = "frame", details = TRUE)

    # Chunk-size and link-batch invariance
    testthat::expect_lt(dr_mx(out_perm_1$phi - out_perm_5$phi), 1e-15)
    testthat::expect_lt(dr_mx(out_perm_1$psi_dir - out_perm_5$psi_dir), 1e-15)
    testthat::expect_lt(dr_mx(out_perm_1$psi_or - out_perm_5$psi_or), 1e-15)
    testthat::expect_lt(dr_mx(out_perm_1$psi_ps - out_perm_5$psi_ps), 1e-15)
    testthat::expect_lt(dr_mx(out_perm_1$psi_sel - out_perm_5$psi_sel), 1e-15)

    # Row permutation invariance (aligned by unit_id)
    reord <- match(fx_id$raw$unit_id, fx_perm$raw$unit_id)
    testthat::expect_lt(
      dr_mx(vapply(fit_id$links, `[[`, numeric(1), "mu_dr") - vapply(fit_perm$links, `[[`, numeric(1), "mu_dr")),
      1e-14
    )
    testthat::expect_lt(dr_mx(out_id$phi - out_perm_5$phi[reord, , drop = FALSE]), 1e-14)
    testthat::expect_lt(dr_mx(out_id$psi_dir - out_perm_5$psi_dir[reord, , drop = FALSE]), 1e-13)
    testthat::expect_lt(dr_mx(out_id$psi_or - out_perm_5$psi_or[reord, , drop = FALSE]), 1e-14)
    testthat::expect_lt(dr_mx(out_id$psi_ps - out_perm_5$psi_ps[reord, , drop = FALSE]), 1e-12)
    testthat::expect_lt(dr_mx(out_id$psi_sel - out_perm_5$psi_sel[reord, , drop = FALSE]), 1e-14)
    testthat::expect_lt(dr_mx(out_id$tau - out_perm_5$tau[reord, , drop = FALSE]), 1e-14)
    testthat::expect_lt(dr_mx(out_id$eta - out_perm_5$eta[reord, , drop = FALSE]), 1e-14)

    cov_id <- dr_ns("stream_dr_covariance")(fit_id)
    cov_perm <- dr_ns("stream_dr_covariance")(fit_perm)
    testthat::expect_lt(dr_mx(cov_id$covariance - cov_perm$covariance), 1e-14)

    # Positive case-weight rescaling invariance (w -> 3.75 * w)
    w_scaled <- 3.75 * as.double(fx_id$raw$w)
    Z_frame <- cbind(`(Intercept)` = 1, x1 = fx_id$raw$x1, x2 = fx_id$raw$x2)[, colnames(fx_id$frame$Z_a), drop = FALSE]
    Z_frame[fx_id$raw$A == 0L, ] <- NA_real_
    frame_sc <- dr_ns("construct_ipw_frame_owner")(
      fx_id$panel,
      as.double(fx_id$raw$unit_id),
      Z_frame,
      weights = w_scaled,
      control = fx_id$control,
      frame_g = as.double(fx_id$raw$cohort),
      allow_frame_only = TRUE,
      design_on_A_only = TRUE
    )
    just <- if (fx_id$spec$U %in% c("B", "A")) "known transport under the frozen independent fixed atom law" else NULL
    reader_sc <- dr_ns("construct_cp04_frame_reader")(
      frame_sc,
      as.double(fx_id$raw$unit_id),
      fx_id$Q,
      as.double(fx_id$raw$cohort),
      w_scaled,
      fx_id$control,
      U = fx_id$U,
      rho = fx_id$rho,
      universe = fx_id$spec$U,
      justification = just,
      has_intercept = fx_id$has_intercept
    )
    views_sc <- lapply(seq_len(nrow(fx_id$catalog$links)), function(i) {
      lk <- fx_id$catalog$links[i, , drop = FALSE]
      dr_ns("construct_ipw_link_view")(
        frame_sc,
        as.double(lk$s),
        as.double(lk$t),
        g = as.double(lk$g),
        K = as.double(dr_loader$dr_split_pipe(lk$K)),
        control = fx_id$control,
        link_id = lk$link_id,
        functional_key = lk$functional_key
      )
    })
    fit_sc <- dr_ns("fit_dr_backend")(
      reader_sc,
      views_sc,
      control = fx_id$control,
      glm_control = fx_id$glm_control,
      rank_tol = fx_id$rank_tol,
      outcome_recipe = fx_id$outcome_recipe,
      treatment_columns = fx_id$treatment_columns,
      selection = fx_id$selection,
      selection_model = fx_id$selection_model,
      selection_pooling = fx_id$selection_pooling
    )
    out_sc <- dr_collect(fit_sc, chunk_rows = 12L, link_batch_size = 2L, scale = "frame", details = TRUE)
    testthat::expect_lt(
      dr_mx(vapply(fit_sc$links, `[[`, numeric(1), "mu_dr") - vapply(fit_id$links, `[[`, numeric(1), "mu_dr")),
      1e-14
    )
    testthat::expect_lt(dr_mx(out_sc$phi - out_id$phi), 1e-14)
  }
})

testthat::test_that("single-fit nuisance deduplication, precomputed selection reuse, and resource refusal guards hold", {
  # 1. Selection fitted ONCE, shared OR fitted ONCE (3 deduplicated blocks for 4 links),
  #    distinct treatment roots fitted ONCE (4 distinct (pair, g, K) for 5 links when l1 is duplicated).
  fx7 <- dr_load("F7_date_product_F_within_cohort")
  dup_view <- dr_ns("construct_ipw_link_view")(
    fx7$frame,
    2000,
    2002,
    g = 2005,
    K = 0,
    control = fx7$control,
    link_id = "l1_dup_consumer",
    functional_key = "cmpcat:l1_dup_consumer|g=2005|K=0|pair=2000->2002"
  )
  views5 <- c(fx7$views, list(l1_dup_consumer = dup_view))

  sel_fit <- dr_ns("fit_selection_models")(
    reader = fx7$reader,
    link_views = views5,
    selection_model = fx7$selection_model,
    selection_pooling = fx7$selection_pooling,
    control = fx7$control,
    glm_control = fx7$glm_control
  )
  # Active date cells on F7/F8 are restricted to the 9 retained endpoint/cohort cells (dates 2000, 2002, 2005)
  testthat::expect_identical(length(sel_fit$cells), 9L)
  testthat::expect_identical(length(sel_fit$blocks), 3L)

  fit_fresh <- dr_fit(fx7, views = views5)
  fit_reused <- dr_fit(fx7, views = views5, sel_fit = sel_fit)
  testthat::expect_identical(fit_reused$selection_backend, sel_fit)
  testthat::expect_identical(fit_fresh$accounting$n_selection_cells, 9L)
  testthat::expect_identical(fit_fresh$accounting$n_selection_blocks, 3L)
  testthat::expect_identical(fit_fresh$accounting$n_regression_fits, 1L)
  testthat::expect_identical(fit_fresh$accounting$n_regression_blocks, 3L)
  testthat::expect_identical(fit_fresh$accounting$n_treatment_fits, 4L)
  testthat::expect_identical(length(fit_fresh$links), 5L)
  testthat::expect_identical(fit_fresh$links[[1L]]$model, fit_fresh$links[[5L]]$model)
  testthat::expect_identical(fit_fresh$accounting$full_F_masks_per_model, 0L)
  testthat::expect_false(fit_fresh$accounting$full_F_by_link_IF_retained)
  testthat::expect_false(fit_fresh$accounting$inverse_matrices_retained)
  testthat::expect_false(fit_fresh$accounting$dense_dummy_regression)

  out_fresh <- dr_collect(fit_fresh, chunk_rows = 11L, link_batch_size = 3L, scale = "frame", details = TRUE)
  out_reused <- dr_collect(fit_reused, chunk_rows = 11L, link_batch_size = 3L, scale = "frame", details = TRUE)
  testthat::expect_identical(out_fresh$phi, out_reused$phi)
  testthat::expect_identical(out_fresh$phi[, 1L], out_fresh$phi[, 5L])

  # 2. Precomputed selection adapter preservation on OR and IPW consumers
  or_fresh <- dr_ns("fit_or_selection_models")(
    fx7$reader, fx7$views, fx7$control, fx7$glm_control, fx7$rank_tol, fx7$outcome_recipe,
    selection_model = fx7$selection_model, selection_pooling = fx7$selection_pooling
  )
  or_reused <- dr_ns("fit_or_selection_models")(
    fx7$reader, fx7$views, fx7$control, fx7$glm_control, fx7$rank_tol, fx7$outcome_recipe,
    selection_model = fx7$selection_model, selection_pooling = fx7$selection_pooling,
    sel_fit = fit_fresh$selection_backend
  )
  testthat::expect_identical(or_reused$selection_backend, fit_fresh$selection_backend)
  testthat::expect_identical(
    vapply(or_fresh$links, `[[`, numeric(1), "mu"),
    vapply(or_reused$links, `[[`, numeric(1), "mu")
  )

  # 3. Pre-allocation resource refusal guards
  fx1 <- dr_load("F1_baseline_none_intercept_slope")
  ctrl_tiny_ret <- dr_ns("ipw_owner_control")(
    max_chunk_rows = 13L, max_retained_entries = 500L,
    max_matrix_entries = 2e7, max_work_ops = 1e10, max_output_entries = 1e7
  )
  testthat::expect_error(
    dr_fit(fx1, control = ctrl_tiny_ret),
    class = "cdid_ipw_owner_bound_error"
  )
  ctrl_tiny_work <- dr_ns("ipw_owner_control")(
    max_chunk_rows = 13L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1000L, max_output_entries = 1e7
  )
  testthat::expect_error(
    dr_fit(fx1, control = ctrl_tiny_work),
    class = "cdid_ipw_owner_bound_error"
  )
  fit1 <- dr_fit(fx1)
  testthat::expect_error(
    dr_ns("evaluate_dr_if_chunk")(fit1, seq_len(20L), 1:4),
    class = "cdid_dr_baseline_input_error"
  )
  testthat::expect_error(
    dr_ns("replay_dr_if")(fit1, 20L, 2L, function(b) {}),
    class = "cdid_dr_baseline_input_error"
  )
  testthat::expect_error(
    dr_ns("stream_dr_covariance")(fit1, retained = 2e7),
    class = "cdid_ipw_owner_bound_error"
  )

  # 4. Explicit solver/rank controls required on private backend; public est_method = "dr"
  #    is now activated by D48 CP07 and passes the scope gate while unsupported est_method strings refuse
  testthat::expect_error(
    dr_ns("fit_dr_baseline")(fx1$reader, fx1$views, control = fx1$control, rank_tol = 1e-10, outcome_recipe = fx1$outcome_recipe),
    class = "cdid_ipw_control_error"
  )
  testthat::expect_error(
    dr_ns("fit_dr_baseline")(fx1$reader, fx1$views, control = fx1$control, glm_control = fx1$glm_control, outcome_recipe = fx1$outcome_recipe),
    class = "cdid_dr_baseline_input_error"
  )
  testthat::expect_error(
    dr_ns("att_gt")(
      yname = "y", tname = "time", idname = "id", gname = "g",
      data = stop("data must not be evaluated"),
      est_method = "dr", mode = "chained"
    ),
    "data must not be evaluated"
  )
  e_bad_em <- testthat::expect_error(
    dr_ns("att_gt")(
      yname = "y", tname = "time", idname = "id", gname = "g",
      data = stop("data must not be evaluated"),
      est_method = "stdipw", mode = "chained"
    ),
    class = "cdid_chained_request_est_method"
  )
  testthat::expect_no_match(conditionMessage(e_bad_em), "data must not be evaluated")

  # 5. Full-calendar A (including units 10, 22, 34 with no retained pair) vs retained-pair union
  obs_catalog_only <- fx1$R
  obs_catalog_only[, 4L] <- 0
  tb_d_full_A <- dr_fit_table(fx1, pair_outcomes = fx1$pair_outcomes_delta)
  tb_d_cat_A <- dr_ns("fit_dr_baseline_table")(
    X_or = fx1$X_or, V = fx1$V_trt, Y = fx1$Y_obs, obs = obs_catalog_only,
    g = fx1$raw$cohort, w = fx1$raw$w, calendar = fx1$catalog$calendar$time,
    links = fx1$catalog$link_list, outcome_recipe = fx1$outcome_recipe,
    pair_outcomes = fx1$pair_outcomes_delta, glm_control = fx1$glm_control, rank_tol = fx1$rank_tol
  )
  testthat::expect_gt(dr_mx(tb_d_full_A$estimates - tb_d_cat_A$estimates), 1e-4)
  testthat::expect_gt(dr_mx(tb_d_full_A$phi - tb_d_cat_A$phi), 0.04)
})
