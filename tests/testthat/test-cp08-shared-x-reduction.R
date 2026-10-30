# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 stage 6, shared-X extension: the post-fit complete-support structural
# reduction map (cp04_complete_reduction) for structurally common fixed-X
# outcome recipes (p = 1, 2) on baseline OR and DR, including the nonconstant
# shared DR treatment propensity. Every expected value comes from the unchanged
# generic routes of the same build: the unreduced pipeline (reduce = FALSE) on
# the same complete-support fit, the general (non-complete) route forced by
# test-local instrumentation, and the full-catalog IF replay of every original
# link. The reduction never supplies its own truth; no formula of the reducer
# is duplicated here.
#
# STATUS: UNEXECUTED at authoring time (no R, build, install or fit was run in
# the authoring job). The root sequences the independent numerical gates.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures (test-local copies; testthat files do not share definitions)
# ---------------------------------------------------------------------------

# Balanced nonadditive panel with a binary unit covariate `x` and a smooth
# second covariate `x2`, one treated cohort 3 and a never-treated pool.
sx_nonadd_panel <- function(N = 41L, T = 5L, n_g = 18L, seed = 41L, id_perm = FALSE, rows_perm = FALSE) {
  set.seed(seed)
  g <- rep(0, N)
  g[sample.int(N, n_g)] <- 3
  ids <- if (id_perm) sample(1000L + seq_len(5L * N), N) else seq_len(N)
  xi <- stats::rbinom(N, 1, 0.5)
  ue <- stats::rnorm(N)
  d <- data.frame(id = rep(ids, each = T), period = rep(seq_len(T), N))
  u <- match(d$id, ids)
  d$G <- g[u]
  d$x <- xi[u]
  eff <- ifelse(d$G > 0 & d$period >= d$G, 1 + 0.25 * (d$period - d$G), 0)
  d$y <- ue[u] + 0.1 * d$period + 0.35 * sin(0.7 * u * d$period) + 0.2 * d$x * (d$period^2) + eff + stats::rnorm(nrow(d))
  d$cl <- (u - 1L) %/% 4L + 1L
  d$x2 <- cos(d$id * 0.31)
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

# Deterministic small panel: explicit treated units and explicit x = 1 units,
# so both cohorts carry both x values (identified OR slope, nonconstant PS).
sx_small_panel <- function(N = 12L, T = 5L, g_units = c(1L, 3L, 5L, 7L, 9L, 11L),
                           x_units = c(1L, 2L, 3L, 5L, 8L), seed = 55L) {
  set.seed(seed)
  g <- rep(0, N); g[g_units] <- 3
  x <- as.numeric(seq_len(N) %in% x_units)
  ue <- stats::rnorm(N)
  d <- data.frame(id = rep(seq_len(N), each = T), period = rep(seq_len(T), N))
  u <- d$id
  d$G <- g[u]
  d$x <- x[u]
  eff <- ifelse(d$G > 0 & d$period >= d$G, 1 + 0.25 * (d$period - d$G), 0)
  d$y <- ue[u] + 0.1 * d$period + 0.35 * sin(0.7 * u * d$period) + 0.2 * d$x * (d$period^2) + eff + stats::rnorm(nrow(d))
  d$cl <- (u - 1L) %/% 2L + 1L
  d
}

sx_run_pipeline <- function(d, est = "reg", chunk = 7L, reduce = TRUE, weighting = "optimal",
                            pt_assumption = "post", pt_start = NULL, xformla = NULL,
                            weightsname = NULL, control_group = "nevertreated",
                            target_cohorts = NULL, frame = NULL, bstrap = FALSE,
                            cband = FALSE, biters = 99L, seed = 101L,
                            aggregations = c("simple", "group", "dynamic", "calendar"),
                            condition_trigger = NULL, ridge_relative = NULL, ns = api) {
  ctrl_list <- if (is.null(chunk)) list() else list(chunk_rows = chunk)
  if (!is.null(condition_trigger)) ctrl_list$condition_trigger <- condition_trigger
  if (!is.null(ridge_relative)) ctrl_list$ridge_relative <- ridge_relative
  if (!length(ctrl_list)) ctrl_list <- NULL
  resolved <- api$chained_public_resolve_control(ctrl_list, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  prep <- api$chained_request_prepare_inputs(
    data = d, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xformla, weightsname = weightsname, clustervars = "cl", anticipation = 0,
    frame = frame, outcome_regression = identical(est, "reg"),
    doubly_robust = identical(est, "dr"), or_resource = resolved$profile
  )
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = biters, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  tc <- if (is.null(target_cohorts)) prep$target_cohorts else target_cohorts
  ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = built$control, xformla = prep$xformla, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = tc, anticipation = 0L, control_group = control_group, base_period = "varying",
    weighting = weighting, reduce_structural_redundancy = reduce, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = aggregations, cband = cband, alp = 0.05, biters = biters, seed = seed,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = bstrap, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = pt_assumption, pt_start = pt_start, cp04 = prep$cp04
  )
}

# The unchanged general (non-complete-support) route of the same build.
sx_general_ns <- function() {
  e <- new.env(parent = api)
  for (nm in c("run_chained_ipw_pipeline", "build_cp04_public_registry", "or_public_fit", "dr_public_fit",
               "fit_or_baseline", "fit_dr_baseline")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  e$or_baseline_complete_facts <- function(...) NULL
  e
}

sx_dense_E <- function(cert) {
  E <- matrix(0, cert$n_links, length(cert$retained_index))
  E[cbind(cert$retained_index, seq_along(cert$retained_index))] <- 1
  if (length(cert$E$i)) {
    rpos <- match(cert$E$j, cert$retained_index)
    for (k in seq_along(cert$E$i)) {
      E[cert$E$i[k], rpos[k]] <- E[cert$E$i[k], rpos[k]] + cert$E$x[k]
    }
  }
  E
}

sx_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, sx_strip); attributes(x) <- a }
  x
}
# Route-only differences: the complete-support record and resource bounds.
sx_backend_view <- function(b) {
  x <- unclass(b)[setdiff(names(b), c("ref", "reader_ref", "complete_support", "accounting"))]
  list(fields = sx_strip(x), accounting = b$accounting[setdiff(names(b$accounting),
    c("complete_support", "work_ops_bound", "fit_scratch_bound", "retained_bound", "owned_entries"))])
}
sx_binding <- function(est) if (est == "dr") "dr_baseline_backend" else "or_baseline_backend"
sx_reseal <- function(b, est) api$ipw_owner_seal(b, sx_binding(est))

# Full-catalog IF replay of every original link on the fitted backend (the
# actual generic per-link IF evaluation; nothing is reconstructed here).
sx_full_if <- function(backend, M, N, chunk) {
  H <- matrix(0, N, M)
  pos <- 0L
  api$cp04_replay_links(backend, seq_len(M), function(batch) {
    idx <- seq.int(pos + 1L, pos + nrow(batch$phi))
    H[idx, ] <<- batch$phi
    pos <<- pos + nrow(batch$phi)
  }, chunk_rows = chunk)
  expect_identical(pos, N)
  H
}

sx_reduction_args <- function(res) {
  M <- length(res$point_links_all)
  list(backend = res$nuisance_registry$cp04, points = res$point_links_all, design = res$design_all,
    dli = match(seq_len(M), res$design_all$links$catalog_row), fo = res$frame_owner,
    cb = res$controls$pipeline_control$class_bounds, oc = res$controls$pipeline_control$owner_control, M = M)
}
sx_reduce <- function(a, backend = a$backend, points = a$points, oc = a$oc, cb = a$cb, fn = api$cp04_complete_reduction, ...) {
  fn(backend, points, a$design, a$dli, a$fo, cb, oc, ...)
}

# ---------------------------------------------------------------------------
# 1. Shared-X p = 1 and p = 2, OR and DR: the 10 original moments and the full
#    original-unit IFs are reconstructed from the 4 forest rows; the global
#    beta and nonconstant PS IF terms are retained; fits, auxiliary IDs and
#    diagnostics are untouched; public results equal both generic routes.
# ---------------------------------------------------------------------------
test_that("shared-X p = 1 and p = 2 OR/DR reduction reconstructs all original moments and full IFs from the forest rows", {
  dp <- sx_nonadd_panel(41L, 5L, 18L, seed = 41L, id_perm = TRUE, rows_perm = TRUE)
  # Keep the original binary-x fixture, and use its existing continuous x2
  # column for a genuinely non-saturated p = 1 DR beta-IF omission witness.
  for (xf in list(~x, ~x2, ~x + x2)) {
    p_exp <- length(all.vars(xf))
    for (est in c("reg", "dr")) {
      info <- paste(est, deparse(xf))
      res_r <- sx_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, reduce = TRUE)
      res_u <- sx_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, reduce = FALSE)
      res_g <- sx_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, reduce = TRUE, ns = sx_general_ns())

      cert <- res_r$structural_certificate
      expect_false(is.null(cert), info = info)
      expect_false(inherits(cert, "cdid_ipw_structural_certificate"), info = info)
      expect_true(api$ipw_owner_seal_ok(cert, "cp04_complete_reduction"), info = info)
      expect_identical(cert$n_links, 10L, info = info)
      expect_identical(length(cert$retained_index), 4L, info = info)
      expect_identical(length(cert$dropped_index), 6L, info = info)
      expect_identical(cert$link_kind, rep(if (est == "dr") "complete_support_DR" else "complete_support_OR", 10L), info = info)
      expect_identical(res_r$gmm_fit$provenance$reduction$mode, "sparse_class_certificate", info = info)
      expect_true(res_r$reduction_summary$reduction_active, info = info)
      expect_identical(res_r$reduction_summary$n_links_used, 4L, info = info)
      # The general route stays generic: no complete-support facts, no map.
      expect_null(res_g$nuisance_registry$cp04$complete_support, info = info)
      expect_null(res_g$structural_certificate, info = info)
      expect_null(res_u$structural_certificate, info = info)

      # The fitted backend carries the shared-X facts: p kept slopes, a global
      # beta of that length and, for DR, a nonconstant shared propensity.
      b_r <- res_r$nuisance_registry$cp04
      b_u <- res_u$nuisance_registry$cp04
      expect_identical(b_r$shared$p, p_exp, info = info)
      expect_identical(length(b_r$shared$beta), p_exp, info = info)
      expect_identical(b_r$outcome_recipe$p, p_exp, info = info)
      if (est == "dr") {
        expect_identical(b_r$treatment_recipe$p, p_exp + 1L, info = info)
        coef1 <- b_r$treatments[[1L]]$coefficients
        expect_identical(length(coef1), p_exp + 1L, info = info)
        expect_gt(max(abs(coef1[-1L])), 1e-3)
        expect_identical(length(unique(lapply(b_r$treatments, `[[`, "coefficients"))), 1L, info = info)
      }

      # Every original fit, point link, auxiliary index/report and public diagnostic is preserved bitwise.
      expect_identical(sx_backend_view(b_r), sx_backend_view(b_u), info = info)
      expect_identical(res_r$point_links_all, res_u$point_links_all, info = info)
      expect_identical(res_r$nuisance_registry$cp04_points, res_u$nuisance_registry$cp04_points, info = info)
      expect_identical(res_r$nuisance_registry$cp04_aux_index, res_u$nuisance_registry$cp04_aux_index, info = info)
      expect_identical(sx_strip(res_r$nuisance_registry$cp04_aux), sx_strip(res_u$nuisance_registry$cp04_aux), info = info)
      expect_identical(sx_strip(res_r$auxiliary_reporting), sx_strip(res_u$auxiliary_reporting), info = info)
      expect_identical(res_r$reporting$target_meta, res_u$reporting$target_meta, info = info)
      expect_identical(res_r$reporting$target_meta, res_g$reporting$target_meta, info = info)
      diag_r <- api$cp04_public_diagnostics(res_r, res_r$controls$pipeline_control)
      diag_u <- api$cp04_public_diagnostics(res_u, res_u$controls$pipeline_control)
      expect_identical(sx_strip(diag_r), sx_strip(diag_u), info = info)

      # Actual fitted moments and the full original-unit IFs of ALL 10 links
      # (generic per-link replay on the fitted backend; core scale = frame scale, n_A == N).
      reader <- res_r$nuisance_registry$cp04_reader
      N <- reader$state$N
      m_all <- vapply(res_r$point_links_all, `[[`, numeric(1L), "estimate")
      H_all <- sx_full_if(b_r, 10L, N, 7L)
      E <- sx_dense_E(cert)
      m_R <- m_all[cert$retained_index]
      H_R <- H_all[, cert$retained_index, drop = FALSE]
      V_all <- crossprod(sweep(H_all, 2L, colMeans(H_all), "-")) / N
      V_R <- crossprod(sweep(H_R, 2L, colMeans(H_R), "-")) / N
      expect_lt(max(abs(m_all - drop(E %*% m_R))), 1e-12)
      expect_lt(max(abs(H_all - H_R %*% t(E))), 1e-12)
      expect_lt(max(abs(V_all - E %*% V_R %*% t(E))), 1e-12)
      expect_lt(res_r$reduction_summary$max_moment_expansion_diff, 1e-12)
      expect_equal(api$cp04_expand_reduction_map(cert, m_R), drop(E %*% m_R), tolerance = 1e-12)
      expect_equal(t(api$cp04_expand_reduction_map(cert, t(H_R))), H_R %*% t(E), tolerance = 1e-12)
      # The reduced covariance stream is exactly the retained columns of the full-catalog stream.
      expect_equal(unname(res_r$covariance_stream$V), unname(V_R), tolerance = 1e-12)
      expect_equal(unname(res_u$covariance_stream$V), unname(V_all), tolerance = 1e-12)

      # Global-beta / PS omission mutants: the reconstructed IF equals the full
      # IF (above) and is far from any IF that drops the global slope term or
      # the propensity score term. Both terms are nontrivial on this design.
      if (est == "dr") {
        psi_or <- psi_ps <- psi_sel <- psi_dir <- matrix(0, N, 10L)
        catalog_cp04 <- vapply(res_r$point_links_all, `[[`, integer(1L), "cp04_index")
        api$replay_dr_baseline_if(b_r, 7L, 10L, function(b) {
          dst <- match(b$links, catalog_cp04)
          take <- which(!is.na(dst))
          if (!length(take)) return(invisible(NULL))
          psi_or[b$rows, dst[take]] <<- b$details$psi_or[, take, drop = FALSE]
          psi_ps[b$rows, dst[take]] <<- b$details$psi_ps[, take, drop = FALSE]
          psi_sel[b$rows, dst[take]] <<- b$details$psi_sel[, take, drop = FALSE]
          psi_dir[b$rows, dst[take]] <<- b$details$psi_dir[, take, drop = FALSE]
        }, scale = "frame", details = TRUE)
        expect_identical(max(abs(psi_sel)), 0)
        expect_lt(max(abs(H_all - (psi_dir + psi_or + psi_ps + psi_sel))), 1e-12)
        # Intercept + binary x is a saturated PS: its weighted covariate
        # contrast, hence the DR beta loading, vanishes in exact arithmetic.
        # Do not claim an omission witness from solver roundoff. The original
        # binary case retains every numeric reconstruction/output check;
        # continuous p = 1 and p = 2 provide the nonzero beta witnesses.
        saturated_binary_ps <- identical(all.vars(xf), "x")
        H_no_beta <- H_all - psi_or
        H_no_ps <- H_all - psi_ps
        if (!saturated_binary_ps) {
          expect_gt(max(abs(psi_or)), 1e-6)
          expect_gt(max(abs(H_R %*% t(E) - H_no_beta)), 1e-6)
        }
        expect_gt(max(abs(psi_ps)), 1e-6)
        expect_gt(max(abs(H_R %*% t(E) - H_no_ps)), 1e-6)
        # The PS is nonconstant across units: the fitted odds vary.
        fo <- res_r$frame_owner
        Xt <- fo$Z_a[reader$state$a_index, b_r$treatment_recipe$index, drop = FALSE]
        odds <- exp(drop(Xt %*% b_r$treatments[[1L]]$coefficients))
        expect_gt(max(odds) - min(odds), 1e-3)
      } else {
        # OR: the detailed replay exposes xi_beta (the global slope IF). Setting
        # xi_beta to zero changes each link's IF by + xi_beta (c_link - mean_d_block)
        # (both fitted descriptors of that link's block/treated contrast).
        xi_beta <- matrix(0, N, p_exp)
        api$replay_or_baseline_if(b_r, 7L, 10L, function(b) {
          xi_beta[b$rows, ] <<- b$details$beta
        }, scale = "frame", details = TRUE)
        expect_gt(max(abs(xi_beta)), 1e-6)
        H_no_beta <- H_all
        for (j in seq_len(10L)) {
          l <- b_r$links[[res_r$point_links_all[[j]]$cp04_index]]
          rg <- b_r$regressions[[l$regression]]
          H_no_beta[, j] <- H_all[, j] + drop(xi_beta %*% (l$c - rg$mean_d))
        }
        expect_gt(max(abs(H_no_beta - H_all)), 1e-6)
        expect_gt(max(abs(H_R %*% t(E) - H_no_beta)), 1e-6)
      }

      # Public results: reduced == unreduced complete route == general route at the existing tolerances.
      for (other in list(res_u, res_g)) {
        expect_equal(res_r$reporting$estimates, other$reporting$estimates, tolerance = 1e-12, info = info)
        expect_equal(res_r$reporting$se, other$reporting$se, tolerance = 1e-12, info = info)
        expect_equal(res_r$reporting$analytical_se, other$reporting$analytical_se, tolerance = 1e-12, info = info)
        expect_equal(res_r$reporting$if_full, other$reporting$if_full, tolerance = 1e-12, info = info)
        expect_equal(res_r$reporting$vcov, other$reporting$vcov, tolerance = 1e-12, info = info)
        for (ag in c("simple", "group", "dynamic", "calendar")) {
          expect_equal(res_r$aggregations[[ag]]$overall$estimate, other$aggregations[[ag]]$overall$estimate, tolerance = 1e-12, info = info)
          expect_equal(res_r$aggregations[[ag]]$overall$se, other$aggregations[[ag]]$overall$se, tolerance = 1e-12, info = info)
          expect_equal(res_r$aggregations[[ag]]$estimates, other$aggregations[[ag]]$estimates, tolerance = 1e-12, info = info)
          expect_equal(res_r$aggregations[[ag]]$se, other$aggregations[[ag]]$se, tolerance = 1e-12, info = info)
        }
      }
      # Full E'E metric semantics with the original M on the reduced route.
      ch <- res_r$covariance_handling
      expect_identical(ch$dimensions$M_full, 10L, info = info)
      expect_identical(ch$dimensions$L_retained, 4L, info = info)
      expect_identical(ch$gram_provenance, "sparse_structural_reduction_E_t_E", info = info)
      expect_equal(ch$trace_B, res_u$covariance_handling$trace_B, tolerance = 1e-12, info = info)
      expect_equal(ch$scale, res_u$covariance_handling$scale, tolerance = 1e-12, info = info)
      # Established D26 behavior: removing exact cycle redundancies leaves
      # a well-conditioned retained covariance; the singular full covariance
      # triggers its ordinary handler. Full-M scale/trace and scientific
      # outputs above still agree. The forced-ridge block below checks equal
      # lambda when both representations invoke the Tikhonov handler.
      expect_identical(ch$handler_used, "none", info = info)
      expect_identical(ch$trigger_reason, "well_conditioned_pd", info = info)
      expect_identical(res_u$covariance_handling$handler_used, "tikhonov", info = info)
    }
  }
})

# ---------------------------------------------------------------------------
# 2. 190 -> 19: the shared-X T = 20 panel replays and streams only the 19
#    forest columns while every one of the 190 fitted moments is reconstructed.
# ---------------------------------------------------------------------------
test_that("shared-X 190-to-19 moment and covariance column reduction with full-catalog reconstruction", {
  d20 <- sx_nonadd_panel(45L, 20L, 20L, seed = 77L)
  for (est in c("reg", "dr")) {
    e_cnt <- new.env(parent = api)
    for (nm in c("run_chained_ipw_pipeline", "cp04_replay_links")) {
      f <- get(nm, api); environment(f) <- e_cnt; assign(nm, f, e_cnt)
    }
    cov_fn <- get("cp04_stream_covariance", api); environment(cov_fn) <- e_cnt
    cov_cols <- integer()
    rep_cols <- integer()
    e_cnt$cp04_stream_covariance <- function(backend, links, retained = 0) {
      cov_cols <<- c(cov_cols, length(links))
      cov_fn(backend, links, retained)
    }
    e_cnt$cp04_replay_links <- function(backend, links, callback, scale = "core", chunk_rows = NULL) {
      rep_cols <<- c(rep_cols, length(links))
      api$cp04_replay_links(backend, links, callback, scale = scale, chunk_rows = chunk_rows)
    }
    res20 <- sx_run_pipeline(d20, est = est, chunk = 50L, reduce = TRUE, xformla = ~x, aggregations = "simple", ns = e_cnt)
    expect_identical(res20$reduction_summary$n_links_all, 190L)
    expect_identical(res20$reduction_summary$n_links_used, 19L)
    expect_true(res20$reduction_summary$reduction_active)
    expect_identical(res20$nuisance_registry$cp04$shared$p, 1L)
    if (est == "dr") expect_identical(res20$nuisance_registry$cp04$treatment_recipe$p, 2L)
    expect_identical(cov_cols, 19L)
    # 2 covariance passes (19 cols each) + 2 GMM IF projection passes (19 cols each).
    expect_identical(rep_cols, rep(19L, 4L))
    # All 190 fitted moments are reconstructed from the 19 (the pipeline's own
    # fitted-versus-expanded self-test on this route) and the map is sealed.
    cert <- res20$structural_certificate
    expect_identical(cert$n_links, 190L)
    expect_identical(length(cert$retained_index), 19L)
    expect_lt(res20$reduction_summary$max_moment_expansion_diff, 1e-12)
    expect_identical(res20$covariance_handling$dimensions$M_full, 190L)
    expect_identical(res20$covariance_handling$dimensions$L_retained, 19L)
    # The full-catalog IF of all 190 links (generic replay) telescopes exactly.
    b <- res20$nuisance_registry$cp04
    N <- res20$nuisance_registry$cp04_reader$state$N
    H_all <- sx_full_if(b, 190L, N, min(50L, b$control$max_chunk_rows))
    E <- sx_dense_E(cert)
    expect_lt(max(abs(H_all - H_all[, cert$retained_index, drop = FALSE] %*% t(E))), 1e-12)
    # Unreduced generic reference on the same fit: public results agree.
    res20_u <- sx_run_pipeline(d20, est = est, chunk = 50L, reduce = FALSE, xformla = ~x, aggregations = "simple")
    expect_null(res20_u$structural_certificate)
    expect_identical(res20_u$reduction_summary$n_links_used, 190L)
    expect_equal(res20$reporting$estimates, res20_u$reporting$estimates, tolerance = 1e-12)
    expect_equal(res20$reporting$se, res20_u$reporting$se, tolerance = 1e-12)
    expect_equal(res20$reporting$if_full, res20_u$reporting$if_full, tolerance = 1e-12)
    expect_equal(res20$aggregations$simple$overall$estimate, res20_u$aggregations$simple$overall$estimate, tolerance = 1e-12)
    expect_equal(res20$aggregations$simple$overall$se, res20_u$aggregations$simple$overall$se, tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# 3. Ridge-active full-metric behavior with X: identical Tikhonov lambda and
#    G = E'E keep reduced and unreduced criteria equal.
# ---------------------------------------------------------------------------
test_that("shared-X ridge-active witness preserves the full-catalog E'E metric and refusal semantics", {
  d_ridge <- sx_small_panel()
  for (est in c("reg", "dr")) {
    res_r <- sx_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = TRUE, xformla = ~x, condition_trigger = 0.99)
    res_u <- sx_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = FALSE, xformla = ~x, condition_trigger = 0.99)
    expect_false(is.null(res_r$structural_certificate))
    expect_identical(res_r$nuisance_registry$cp04$shared$p, 1L)
    ch_r <- res_r$covariance_handling
    ch_u <- res_u$covariance_handling
    expect_identical(ch_r$handler_used, "tikhonov")
    expect_identical(ch_u$handler_used, "tikhonov")
    expect_identical(ch_r$dimensions$M_full, 10L)
    expect_identical(ch_r$dimensions$L_retained, 4L)
    expect_false(ch_r$gram_is_identity)
    expect_false(ch_r$gram_is_diagonal)
    expect_equal(ch_r$trace_B, ch_u$trace_B, tolerance = 1e-12)
    expect_equal(ch_r$scale, ch_u$scale, tolerance = 1e-12)
    expect_equal(ch_r$lambda, ch_u$lambda, tolerance = 1e-12)
    expect_equal(res_r$reporting$estimates, res_u$reporting$estimates, tolerance = 1e-12)
    expect_equal(res_r$reporting$se, res_u$reporting$se, tolerance = 1e-12)
    expect_equal(res_r$reporting$if_full, res_u$reporting$if_full, tolerance = 1e-12)
    expect_equal(res_r$reporting$vcov, res_u$reporting$vcov, tolerance = 1e-12)
    for (ag in c("simple", "group", "dynamic", "calendar")) {
      expect_equal(res_r$aggregations[[ag]]$overall$estimate, res_u$aggregations[[ag]]$overall$estimate, tolerance = 1e-12)
      expect_equal(res_r$aggregations[[ag]]$overall$se, res_u$aggregations[[ag]]$overall$se, tolerance = 1e-12)
    }
    res_r_big <- sx_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = TRUE, xformla = ~x, condition_trigger = 0.99, ridge_relative = 0.5)
    res_u_big <- sx_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = FALSE, xformla = ~x, condition_trigger = 0.99, ridge_relative = 0.5)
    expect_equal(res_r_big$reporting$estimates, res_u_big$reporting$estimates, tolerance = 1e-12)
    expect_equal(res_r_big$reporting$se, res_u_big$reporting$se, tolerance = 1e-12)
  }
  # Constant non-zero treatment effect with zero within-cohort variance refuses identically with X.
  d_zero1 <- sx_small_panel()
  d_zero1$y <- ifelse(d_zero1$G == 3 & d_zero1$period >= 3, 2 * (d_zero1$period - 2), 0)
  for (est in c("reg", "dr")) {
    err_r <- tryCatch(sx_run_pipeline(d_zero1, est = est, chunk = 7L, reduce = TRUE, xformla = ~x), error = identity)
    err_u <- tryCatch(sx_run_pipeline(d_zero1, est = est, chunk = 7L, reduce = FALSE, xformla = ~x), error = identity)
    expect_s3_class(err_r, "cdid_ipw_pipeline_scale_error")
    expect_identical(class(err_r), class(err_u))
  }
})

# ---------------------------------------------------------------------------
# 4. Exact ineligibility: designs that are not structurally common fixed-X
#    complete-roster fits stay generic BEFORE any work; sealed-fact mutants of
#    the new bindings return NULL; no numerical comparison is reached; the
#    full map work is composed before construction.
# ---------------------------------------------------------------------------
test_that("shared-X exact ineligibility negatives, binding mutants and work composition", {
  d0 <- sx_nonadd_panel(41L, 5L, 18L, seed = 41L)
  dm <- d0[!(d0$id == 5L & d0$period == 3L), ]
  dw <- d0; dw$w <- 1 + (dw$id %% 3) / 4
  d2 <- d0; d2$G[d2$id %in% c(4L, 9L, 14L, 19L)] <- 4
  fr <- d0[!duplicated(d0$id), c("id", "G", "x", "cl")]
  fr <- rbind(fr[order(fr$id), ], data.frame(id = 9999L, G = 0, x = 0, cl = 99L))

  neg_cases <- list(
    missing_date_or = function() sx_run_pipeline(dm, est = "reg", xformla = ~x),
    missing_date_dr = function() sx_run_pipeline(dm, est = "dr", xformla = ~x),
    varying_weights_or = function() sx_run_pipeline(dw, est = "reg", xformla = ~x, weightsname = "w"),
    varying_weights_dr = function() sx_run_pipeline(dw, est = "dr", xformla = ~x, weightsname = "w"),
    frame_only_unit = function() sx_run_pipeline(d0, est = "reg", xformla = ~x, frame = fr),
    multiple_targets = function() sx_run_pipeline(d2, est = "dr", xformla = ~x, control_group = "notyettreated"),
    identity_weighting = function() sx_run_pipeline(d0, est = "dr", xformla = ~x, weighting = "identity"),
    pt_all = function() sx_run_pipeline(d0, est = "dr", xformla = ~x, pt_assumption = "all", pt_start = 1)
  )
  for (nm in names(neg_cases)) {
    res <- neg_cases[[nm]]()
    expect_null(res$structural_certificate, info = nm)
    expect_false(res$reduction_summary$reduction_active, info = nm)
    expect_identical(res$reduction_summary$n_links_used, res$reduction_summary$n_links_all, info = nm)
  }

  for (est in c("reg", "dr")) {
    res <- sx_run_pipeline(d0, est = est, chunk = 7L, xformla = ~x + x2, reduce = TRUE)
    a <- sx_reduction_args(res)
    cert <- res$structural_certificate
    expect_false(is.null(cert), info = est)
    expect_identical(a$backend$shared$p, 2L, info = est)
    expect_identical(sx_reduce(a)$E, cert$E, info = est)

    # Existing bound/refusal semantics are preserved on the shared-X fit.
    expect_null(sx_reduce(a, backend = NULL))
    bad_pl <- a$points; bad_pl[[2L]]$link_id <- "wrong_id"
    expect_null(sx_reduce(a, points = bad_pl))
    cb_links <- a$cb; cb_links$max_links <- 5L
    expect_error(sx_reduce(a, cb = cb_links), class = "cdid_ipw_bound_error")
    cb_e <- a$cb; cb_e$max_E_entries <- 2L
    expect_error(sx_reduce(a, cb = cb_e), class = "cdid_ipw_bound_error")
    oc_work <- a$oc; oc_work$max_work_ops <- a$backend$accounting$work_ops_bound + 10
    expect_error(sx_reduce(a, oc = oc_work), class = "cdid_ipw_owner_bound_error")

    # The full map work is composed before any design/forest/path construction.
    e <- new.env(parent = api)
    fn <- api$cp04_complete_reduction; environment(fn) <- e
    calls <- c(binding_scan = 0L, design = 0L, forest = 0L, path = 0L)
    e$lapply <- function(...) { calls["binding_scan"] <<- calls["binding_scan"] + 1L; base::lapply(...) }
    e$validate_ipw_design_for_classes <- function(...) { calls["design"] <<- calls["design"] + 1L; api$validate_ipw_design_for_classes(...) }
    e$incidence_forest <- function(...) { calls["forest"] <<- calls["forest"] + 1L; api$incidence_forest(...) }
    e$ipw_forest_path <- function(...) { calls["path"] <<- calls["path"] + 1L; api$ipw_forest_path(...) }
    expect_error(sx_reduce(a, oc = oc_work, fn = fn), class = "cdid_ipw_owner_bound_error")
    expect_identical(unname(calls), rep(0L, 4L))
    expect_identical(sx_reduce(a, fn = fn)$E, cert$E)
    expect_gt(calls[["binding_scan"]], 0L)
    expect_gt(calls[["forest"]], 0L)
    expect_gt(calls[["path"]], 0L)

    # No numerical-near-equality predicate licenses the shared-X map: the
    # copied helper's abs trap is never reached on valid p = 2 input.
    e$abs <- function(...) stop("numerical eligibility comparison reached")
    expect_identical(sx_reduce(a, fn = fn)$E, cert$E)
    tiny_pair <- a$points
    tiny_pair[[1L]]$pair$s <- tiny_pair[[1L]]$pair$s + .Machine$double.eps
    expect_null(sx_reduce(a, points = tiny_pair, fn = fn))

    # Sealed-fact mutants of the new structural bindings (each re-sealed, so
    # only the binding itself decides): every one returns NULL, exactly.
    b0 <- a$backend
    idx1 <- a$points[[1L]]$cp04_index
    l1 <- b0$links[[idx1]]
    mut <- list()
    m1 <- b0; m1$supports[[l1$support_id]]$sumZ <- m1$supports[[l1$support_id]]$sumZ + 1; mut$target_sumZ <- m1
    csid <- b0$regression_specs[[l1$regression]]$support_ids[1L]
    m2 <- b0; m2$supports[[csid]]$sumZ <- m2$supports[[csid]]$sumZ + 1; mut$control_sumZ <- m2
    m3 <- b0; m3$regressions[[l1$regression]]$L <- m3$regressions[[l1$regression]]$L + 1; mut$block_L <- m3
    m4 <- b0; m4$supports[[l1$support_id]]$L <- m4$supports[[l1$support_id]]$L + 1; mut$support_L <- m4
    m5 <- b0; m5$shared$fit$rank <- 1L; mut$rank <- m5
    m6 <- b0; m6$shared$basis$frozen_during_differentiation <- FALSE; mut$frozen <- m6
    m7 <- b0; m7$shared$keep <- m7$shared$keep[1L]; mut$keep_dim <- m7
    m8 <- b0; m8$shared$beta <- m8$shared$beta[1L]; mut$beta_dim <- m8
    m9 <- b0; m9$complete_support$p <- 0L; mut$facts_p <- m9
    m10 <- b0; m10$regressions[[l1$regression]]$mean_d <- m10$regressions[[l1$regression]]$mean_d[1L]; mut$mean_d_dim <- m10
    m11 <- b0; m11$links[[idx1]]$c <- m11$links[[idx1]]$c[1L]; mut$link_c_dim <- m11
    m12 <- b0; m12$shared$control_support_ids <- m12$shared$control_support_ids[1L]; mut$control_ids <- m12
    m13 <- b0; m13$supports[[l1$support_id]]$rate <- 1 - .Machine$double.eps; mut$rate <- m13
    if (est == "dr") {
      tid <- l1$model
      n1 <- b0; n1$treatments[[tid]]$coefficients[2L] <- n1$treatments[[tid]]$coefficients[2L] + 0.5; mut$ps_coef <- n1
      n2 <- b0; n2$treatments[[tid]]$partials$beta <- n2$treatments[[tid]]$partials$beta[, 1L, drop = FALSE]; mut$ps_partials_dim <- n2
      n3 <- b0; n3$treatments[[tid]]$cross <- n3$treatments[[tid]]$cross[1L, , drop = FALSE]; mut$ps_cross_dim <- n3
      n4 <- b0; n4$treatments[[tid]]$d_bar_diff <- n4$treatments[[tid]]$d_bar_diff[1L]; mut$ps_dbar_dim <- n4
      n5 <- b0; n5$treatment_recipe$p <- 1L; mut$trt_recipe_p <- n5
      n6 <- b0; n6$treatment_recipe$index <- n6$treatment_recipe$index[c(2L, 1L, 3L)]; mut$trt_recipe_index <- n6
      n7 <- b0; n7$treatments[[tid]]$H_chol <- n7$treatments[[tid]]$H_chol * 2; mut$ps_H <- n7
    }
    for (nm in names(mut)) {
      expect_null(sx_reduce(a, backend = sx_reseal(mut[[nm]], est)), info = paste(est, nm))
    }
    # A deliberate mutant that keeps every binding intact still produces the map:
    # the mutant list above is not vacuous.
    expect_identical(sx_reduce(a, backend = sx_reseal(b0, est))$E, cert$E, info = est)
  }
})

# ---------------------------------------------------------------------------
# 5. The intercept-only route is unchanged: same map, same E, same NULL/refusal
#    outcomes as before the extension on the p = 0 fixture.
# ---------------------------------------------------------------------------
test_that("intercept-only p = 0 reduction is unchanged by the shared-X extension", {
  d <- sx_nonadd_panel(41L, 5L, 18L, seed = 31L)
  for (est in c("reg", "dr")) {
    res_r <- sx_run_pipeline(d, est = est, chunk = 7L, reduce = TRUE)
    res_u <- sx_run_pipeline(d, est = est, chunk = 7L, reduce = FALSE)
    cert <- res_r$structural_certificate
    expect_false(is.null(cert))
    expect_identical(res_r$nuisance_registry$cp04$shared$p, 0L)
    expect_identical(cert$n_links, 10L)
    expect_identical(length(cert$retained_index), 4L)
    expect_equal(res_r$reporting$estimates, res_u$reporting$estimates, tolerance = 1e-12)
    expect_equal(res_r$reporting$se, res_u$reporting$se, tolerance = 1e-12)
    a <- sx_reduction_args(res_r)
    # The unchanged p = 0 work composition: the published preflight term is the
    # only term (no shared-X descriptor charge is added for the intercept route).
    M <- a$M
    D <- min(length(a$backend$complete_support$dates), 2 * M)
    W <- 256 * M * (D + ceiling(log2(M + 1))) + 64 * M * min(M, max(0, D - 1))
    B <- a$backend$accounting$work_ops_bound
    cap <- a$oc; cap$max_work_ops <- B + W - 1
    expect_error(sx_reduce(a, oc = cap), class = "cdid_ipw_owner_bound_error")
    cap$max_work_ops <- B + W
    expect_identical(sx_reduce(a, oc = cap)$E, cert$E)
  }
})
