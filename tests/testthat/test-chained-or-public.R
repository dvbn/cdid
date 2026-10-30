# Root runs this file against the real installed namespace and native DLL.
# These deterministic algebra fixtures are not a coverage simulation.
orp_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  f <- get(name, ns, inherits = FALSE)
  stopifnot(is.function(f), identical(environment(f), ns)); f
}
orp_path <- function(directory, name) testthat::test_path("fixtures", directory, name)
orp_raw <- read.csv(orp_path("cp06_or_baseline", "raw.csv"))
orp_raw$cohort <- as.numeric(replace(orp_raw$cohort, orp_raw$cohort == "NT", "Inf"))
sys.source(testthat::test_path("reference", "cp06_public_or_reference.R"), envir = environment())
sys.source(testthat::test_path("reference", "cp06_pooled_public_or_reference.R"), envir = environment())

orp_fixture <- function(route = "none", dates = 1:4) {
  d <- orp_raw; N <- nrow(d)
  y_cols <- if (all(dates %in% 1:4)) paste0("Y", dates) else paste0("Y", seq_along(dates))
  Y <- as.matrix(d[y_cols])
  ids <- paste0("original-unit-", d$id)
  g_mapped <- if (length(dates) == 4L && !identical(as.numeric(dates), as.numeric(1:4))) {
    ifelse(is.infinite(d$cohort), 0,
      ifelse(d$cohort == 3, dates[3L], ifelse(d$cohort == 4, dates[4L], max(dates) + 60)))
  } else {
    ifelse(is.infinite(d$cohort), 0, d$cohort)
  }
  data <- data.frame(id = rep(ids, each = length(dates)), time = rep(dates, N),
    y = as.vector(t(Y)), g = rep(g_mapped, each = length(dates)))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  # Retain one all-NA panel identity and the rest as external frame-only units.
  data <- data[!data$id %in% ids[zero[-1L]], ]
  frame <- data.frame(id = ids, g = g_mapped, x = d$x,
    duplicate = d$x, combination = 2 + 3 * d$x, const_col = rep(5, N),
    w = d$w, cl = paste0("cluster-", (d$id - 1L) %/% 3L), stringsAsFactors = FALSE)
  frame[["(Intercept)"]] <- d$x + 1.5
  if (route %in% c("B", "A")) {
    tab <- read.csv(orp_path("cp06_or_corollary1", "transport.csv"))
    tab <- tab[tab$route == route, ]; tab$cohort <- as.numeric(replace(tab$cohort, tab$cohort == "NT", "Inf"))
    frame$rho <- tab$rho[match(paste(d$cohort, d$x), paste(tab$cohort, tab$x))]
  }
  args <- list(yname = "y", tname = "time", idname = "id", gname = "g", data = data, frame = frame,
    xformla = ~x, weightsname = "w", est_method = "reg", control_group = "nevertreated",
    bstrap = FALSE, cband = FALSE, gmm_weighting = "identity")
  if (route != "none") {
    args$selection <- "corollary1"; args$selection_universe <- route; args$selection_xformla <- ~x
    if (route != "F") {
      args$selection_rho <- "rho"; args$selection_justification <- "frozen fixture population-law transport"
    }
  }
  list(args = args, frame = frame, data = data)
}

orp_fit <- function(fx, ...) do.call(cdid::att_gt, utils::modifyList(fx$args, list(...)))

orp_backend <- function(fx, ...) {
  args <- utils::modifyList(fx$args, list(...))
  control_group <- if (is.null(args$control_group)) "nevertreated" else args$control_group
  anticipation <- if (is.null(args$anticipation)) 0L else as.integer(args$anticipation)
  base_period <- if (is.null(args$base_period)) "varying" else args$base_period
  pt_assumption <- if (is.null(args$pt_assumption)) "post" else args$pt_assumption
  pt_start <- args$pt_start
  selection <- if (is.null(args$selection)) "none" else args$selection
  selection_universe <- if (is.null(args$selection_universe)) "F" else args$selection_universe
  gmm_weighting <- if (is.null(args$gmm_weighting)) "identity" else args$gmm_weighting
  bstrap <- if (is.null(args$bstrap)) FALSE else args$bstrap
  cband <- if (is.null(args$cband)) FALSE else args$cband
  biters <- if (is.null(args$biters)) 1000L else as.integer(args$biters)
  alp <- if (is.null(args$alp)) 0.05 else args$alp

  resolved <- orp_ns("chained_public_resolve_control")(args[["control"]],
    gmm_weighting = gmm_weighting, gmm_covariance = "tikhonov")
  prep <- orp_ns("chained_request_prepare_inputs")(
    data = args$data, yname = args$yname, tname = args$tname, idname = args$idname, gname = args$gname,
    xformla = args$xformla, weightsname = args$weightsname, clustervars = args$clustervars,
    anticipation = anticipation, frame = args$frame, target_cohorts = args$target_cohorts,
    selection = selection, selection_universe = selection_universe,
    selection_xformla = args$selection_xformla, selection_rho = args$selection_rho,
    selection_justification = args$selection_justification,
    outcome_regression = TRUE, or_resource = resolved$profile
  )
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- orp_ns("chained_public_dimensions")(prep, biters = biters, clustered = !is.null(prep$cluster_full))
  built <- orp_ns("chained_public_pipeline_control")(dims, resolved,
    gmm_weighting = gmm_weighting, gmm_covariance = "tikhonov")
  res <- orp_ns("run_chained_ipw_pipeline")(
    data = prep$data_internal, idname = args$idname, tname = args$tname, yname = args$yname, gname = args$gname,
    control = built$control, xformla = prep$xformla, weightsname = NULL,
    cluster_vec = prep$cluster_full, target_cohorts = prep$target_cohorts,
    anticipation = anticipation, control_group = control_group,
    base_period = base_period, weighting = gmm_weighting,
    run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = character(0L), cband = cband, alp = alp, biters = biters,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w,
    frame_g = prep$frame_g_arg, allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE,
    bstrap = bstrap, caller_rng = TRUE, gmm_covariance = "tikhonov",
    pt_assumption = pt_assumption, pt_start = pt_start, cp04 = prep$cp04
  )
  res$nuisance_registry$cp04
}

orp_state <- new.env(parent = emptyenv())

testthat::test_that("public baseline and F/B/A match independent pooled dense-functional reference and full-F nuisance/share rows", {
  for (route in c("none", "F", "B", "A")) {
    fx <- orp_fixture(route); fit <- orp_fit(fx)
    raw_r <- orp_raw
    if (route %in% c("B", "A")) raw_r$rho <- fx$frame$rho
    ref <- cp06_pooled_public_or_reference(raw_r, route = route, intercept = TRUE,
      slope_mat = cbind(x = raw_r$x), Q_mat = cbind(`(Intercept)` = 1, x = raw_r$x),
      dates = 1:4, control_group = "nevertreated")

    testthat::expect_equal(fit$att, ref$att, tolerance = 1e-10)
    testthat::expect_equal(unname(fit$inffunc), unname(ref$inffunc), tolerance = 1e-7)
    testthat::expect_equal(unname(fit$V_analytical), unname(ref$V_analytical), tolerance = 1e-7)
    testthat::expect_equal(fit$se, sqrt(diag(ref$V_analytical) / fit$n), tolerance = 1e-8)
    testthat::expect_identical(rownames(fit$inffunc), fx$frame$id)
    testthat::expect_equal(fit$n, nrow(fx$frame))
    testthat::expect_gt(sum(fit$chained$frame$n_obs == 0), 0)
    testthat::expect_gt(sum(fit$chained$frame$n_obs == 1), 0)

    diag <- fit$chained$nuisance_diagnostics
    testthat::expect_identical(diag$method, "reg")
    testthat::expect_identical(diag$outcome_model, "shared one-period outcome-change control-mean regression")
    testthat::expect_identical(diag$fit_structure, "single_joint_fit_with_block_local_intercepts")
    testthat::expect_identical(diag$duration_convention, "stored-calendar transition count L_e = i_t - i_s")
    testthat::expect_equal(diag$coefficients$shared_slopes, ref$beta, tolerance = 1e-10)
    testthat::expect_equal(unname(diag$coefficients$local_intercepts), ref$alpha, tolerance = 1e-10)
    testthat::expect_true(diag$coefficients$has_local_intercepts)
    testthat::expect_equal(diag$fitted_catalog$mu, ref$mu, tolerance = 1e-10)
    testthat::expect_equal(diag$accounting$n_treatment_fits, 0)
    testthat::expect_equal(diag$accounting$n_regression_fits, 1)
    testthat::expect_equal(diag$accounting$n_regression_blocks, length(ref$alpha))
    testthat::expect_equal(diag$accounting$n_local_intercepts, length(ref$alpha))
    testthat::expect_equal(diag$accounting$n_shared_slopes, length(ref$beta))
    testthat::expect_equal(diag$accounting$n_selection_fits == 0, route == "none")
    testthat::expect_null(diag$records$treatment)
    testthat::expect_true(!is.null(diag$shared$fit$R))
    testthat::expect_equal(dim(diag$shared$fit$R), rep(length(ref$beta), 2L))
    testthat::expect_true(all(vapply(diag$records$regression, function(r) is.null(r$qr), logical(1))))
    testthat::expect_true(all(vapply(diag$records$regression, `[[`, character(1), "role") ==
      "local_intercept_block_of_single_joint_fit"))
    testthat::expect_no_match(fit$chained$selection$assumptions, "Synthetic treatment")
    testthat::expect_match(fit$chained$selection$outcome_model, "shared one-period outcome-change")
    testthat::expect_equal(fit$chained$reduction_summary$n_links_all, fit$chained$reduction_summary$n_links_used)
    testthat::expect_identical(fit$DIDparams$est_method, "reg")
    testthat::expect_identical(fit$DIDparams$xformla, fx$args$xformla)
    testthat::expect_match(fit$chained$public$entry, 'est_method = "reg"', fixed = TRUE)

    # Verify global shared beta, local block alpha, link mu, and full-F share/nuisance rows
    b <- orp_backend(fx)
    N <- nrow(raw_r); L <- length(b$links)
    eval_fn <- if (route == "none") orp_ns("evaluate_or_baseline_if_chunk") else orp_ns("evaluate_or_corollary1_if_chunk")
    det <- eval_fn(b, seq_len(N), seq_len(L), scale = "frame", details = TRUE)
    ord_u <- match(fx$frame$id, det$id)
    testthat::expect_setequal(det$id, fx$frame$id)
    testthat::expect_identical(det$id[ord_u], fx$frame$id)
    testthat::expect_identical(b$shared$basis$kept_columns, colnames(ref$phi_beta))
    testthat::expect_equal(unname(det$phi[ord_u, , drop = FALSE]), unname(ref$phi_mu), tolerance = 1e-7)
    testthat::expect_equal(unname(det$details$beta[ord_u, , drop = FALSE]), unname(ref$phi_beta), tolerance = 1e-7)
    testthat::expect_equal(unname(det$details$alpha[ord_u, , drop = FALSE]), unname(ref$phi_alpha), tolerance = 1e-7)
    testthat::expect_equal(unname(det$details$tau[ord_u, , drop = FALSE]), unname(ref$phi_tau), tolerance = 1e-7)
    testthat::expect_equal(unname(det$details$eta[ord_u, , drop = FALSE]), unname(ref$phi_eta), tolerance = 1e-7)
    if (route %in% c("B", "A")) {
      na_cols <- which(det$details$upstream_ids <= b$upstream$n_a)
      coh_idx <- b$upstream$cohort[det$details$upstream_ids[na_cols]]
      testthat::expect_equal(unname(det$details$upstream[ord_u, na_cols, drop = FALSE]),
        unname(ref$phi_a_U[, coh_idx, drop = FALSE]), tolerance = 1e-7)
    }
    # Streaming a single link batch must still evaluate all global beta blocks
    det_single <- eval_fn(b, seq_len(N), 1L, scale = "frame", details = FALSE)
    testthat::expect_equal(det_single$diagnostics$n_global_beta_blocks, length(b$regression_specs))
    testthat::expect_equal(unname(det_single$phi[ord_u, 1L]), unname(ref$phi_mu[, 1L]), tolerance = 1e-7)

    off_A <- which(orp_raw$A == 0)
    if (route %in% c("F", "B")) {
      testthat::expect_gt(max(abs(fit$inffunc[off_A, ])), 1e-4)
      # Nuisance effects can cancel in beta on this fixture. The complete ATT
      # rows must retain the independently computed nonzero off-A effects.
      testthat::expect_equal(unname(fit$inffunc[off_A, , drop = FALSE]),
        unname(ref$inffunc[off_A, , drop = FALSE]), tolerance = 1e-7)
    } else {
      testthat::expect_equal(max(abs(fit$inffunc[off_A, ])), 0, tolerance = 0)
      testthat::expect_equal(max(abs(det$details$beta[ord_u[off_A], ])), 0, tolerance = 0)
    }

    # Verify all four aggregations and their complete F-row IFs (including cohort-share tau derivatives)
    for (type in c("simple", "group", "calendar", "dynamic")) {
      agg <- cdid::aggte(fit, type = type, bstrap = FALSE, cband = FALSE)
      testthat::expect_equal(agg$overall.att, unname(ref$agg_overall[[type]]), tolerance = 1e-10)
      actual_if <- switch(type, simple = agg$inf.function$simple.att,
        group = agg$inf.function$selective.inf.func, calendar = agg$inf.function$calendar.inf.func,
        dynamic = agg$inf.function$dynamic.inf.func)
      testthat::expect_equal(unname(as.vector(actual_if)), unname(ref$phi_agg_overall[, type]), tolerance = 1e-7)
    }

    assign(route, fit, orp_state)
  }
})

testthat::test_that("Corollary1 optimal weighting retains every moment and complete auxiliary rows", {
  for (route in c("F", "B", "A")) {
    base <- get(route, orp_state); fit <- orp_fit(orp_fixture(route), gmm_weighting = "optimal")
    pre <- which(base$chained$role == "pre")
    testthat::expect_equal(fit$att[pre], base$att[pre], tolerance = 1e-11)
    testthat::expect_equal(fit$inffunc[, pre], base$inffunc[, pre], tolerance = 1e-10)
    testthat::expect_equal(fit$chained$reduction_summary$n_links_all, fit$chained$reduction_summary$n_links_used)
    testthat::expect_identical(fit$chained$covariance_handling$weighting_canonical, "optimal")
    testthat::expect_false(fit$chained$numerical_transfer$certified)
  }
})

testthat::test_that("global exact aliases and pooled control prediction-preserving removal preserve point estimates and complete IF rows", {
  for (route in c("none", "F", "B", "A")) {
    fx <- orp_fixture(route); base <- get(route, orp_state)
    aliased <- orp_fit(fx, xformla = ~x + duplicate + combination)
    testthat::expect_equal(aliased$att, base$att, tolerance = 1e-11)
    testthat::expect_equal(aliased$inffunc, base$inffunc, tolerance = 1e-10)
    testthat::expect_equal(aliased$V_analytical, base$V_analytical, tolerance = 1e-11)
    d <- aliased$chained$outcome_design
    testthat::expect_identical(d$original_columns, c("(Intercept)", "x", "duplicate", "combination"))
    testthat::expect_identical(d$aliased_columns, c("duplicate", "combination"))
    testthat::expect_identical(d$retained_columns, c("(Intercept)", "x"))
    testthat::expect_true(d$intercept)
    testthat::expect_identical(d$slope_columns, "x")
    testthat::expect_true(d$frozen_during_differentiation)
    testthat::expect_identical(aliased$chained$selection$formula, base$chained$selection$formula)

    # Constant slope collapsing into intercept retains local block intercepts + duration-scaled slope space
    const_alias <- orp_fit(fx, xformla = ~x + const_col)
    testthat::expect_equal(const_alias$att, base$att, tolerance = 1e-11)
    testthat::expect_equal(const_alias$inffunc, base$inffunc, tolerance = 1e-10)
    testthat::expect_equal(const_alias$V_analytical, base$V_analytical, tolerance = 1e-11)
    dc <- const_alias$chained$outcome_design
    testthat::expect_true(dc$intercept)
    testthat::expect_identical(dc$intercept_column, "(Intercept)")
    testthat::expect_identical(dc$aliased_columns, "const_col")
    testthat::expect_identical(dc$retained_columns, c("(Intercept)", "x"))
    testthat::expect_identical(dc$slope_columns, "x")
  }

  # Pooled-control exact alias (non-aliased on common A due to cohort 4 on dates 1:3, aliased on pooled controls,
  # and reconstructing treated cohort 3) is removed with preserved predictions and IFs.
  # Cohort 4 is beyond this three-date panel and is not a requested target.
  fx_sub <- orp_fixture("none", dates = 1:3)
  fx_sub$args$frame$x_c3_nt <- ifelse(fx_sub$frame$g == 4, fx_sub$frame$x^2 + 1, 2 * fx_sub$frame$x)
  fit_sub_base <- orp_fit(fx_sub, xformla = ~x)
  fit_sub_alias <- orp_fit(fx_sub, xformla = ~x + x_c3_nt)
  testthat::expect_setequal(unique(fit_sub_alias$group), 3)
  testthat::expect_identical(fit_sub_alias$chained$outcome_design$retained_columns, c("(Intercept)", "x", "x_c3_nt"))
  testthat::expect_identical(fit_sub_alias$chained$nuisance_diagnostics$basis$pooled_control$aliased_columns, "x_c3_nt")
  testthat::expect_identical(fit_sub_alias$chained$nuisance_diagnostics$basis$pooled_control$kept_columns, "x")
  testthat::expect_equal(fit_sub_alias$att, fit_sub_base$att, tolerance = 1e-11)
  testthat::expect_equal(fit_sub_alias$inffunc, fit_sub_base$inffunc, tolerance = 1e-10)
  testthat::expect_equal(fit_sub_alias$V_analytical, fit_sub_base$V_analytical, tolerance = 1e-11)
})

testthat::test_that("two-date public OR and all aggregations match independent normalized means and share IFs", {
  fx <- orp_fixture(dates = c(2, 4)); fx$args$xformla <- ~1
  ref <- cp06_mean_reference(orp_raw)
  for (weighting in c("identity", "optimal")) {
    fit <- orp_fit(fx, gmm_weighting = weighting)
    tm <- fit$chained$target_meta; j <- which(tm$role == "post")
    testthat::expect_equal(tm$g[j], c(3, 4))
    testthat::expect_equal(fit$att[j], ref$theta, tolerance = 1e-10)
    testthat::expect_equal(unname(fit$inffunc[, j]), ref$phi, tolerance = 1e-10)
    testthat::expect_equal(unname(fit$V_analytical[j, j] / fit$n), ref$covariance, tolerance = 1e-11)
    testthat::expect_equal(fit$se[j], sqrt(diag(ref$covariance)), tolerance = 1e-11)
    for (type in c("simple", "group", "calendar", "dynamic")) {
      a <- cdid::aggte(fit, type = type, bstrap = FALSE, cband = FALSE)
      ov <- if (type == "dynamic") ref$dynamic else ref$overall
      oi <- if (type == "dynamic") ref$dynamic_if else ref$overall_if
      testthat::expect_equal(a$overall.att, ov, tolerance = 1e-10)
      actual_if <- switch(type, simple = a$inf.function$simple.att,
        group = a$inf.function$selective.inf.func, calendar = a$inf.function$calendar.inf.func,
        dynamic = a$inf.function$dynamic.inf.func)
      testthat::expect_equal(unname(actual_if), oi, tolerance = 1e-10)
      testthat::expect_equal(a$overall.se, sqrt(sum((oi - mean(oi))^2)) / nrow(orp_raw), tolerance = 1e-11)
      testthat::expect_match(a$chained$estimator, "outcome regression/control-mean", fixed = TRUE)
    }
  }
})

testthat::test_that("NT/NYT PT-post/all preserve post/base reporting, match independent pooled reference, and bind the full pretest policy", {
  for (route in c("none", "F", "B", "A")) for (cg in c("nevertreated", "notyettreated")) {
    fx <- orp_fixture(route)
    raw_r <- orp_raw
    if (route %in% c("B", "A")) raw_r$rho <- fx$frame$rho
    for (pt in c("post", "all")) {
      extra <- list(control_group = cg, pt_assumption = pt)
      if (pt == "all") extra$pt_start <- 1
      v <- do.call(orp_fit, c(list(fx), extra))
      u <- do.call(orp_fit, c(list(fx, base_period = "universal"), extra))
      ref_v <- cp06_pooled_public_or_reference(raw_r, route = route, intercept = TRUE,
        slope_mat = cbind(x = raw_r$x), Q_mat = cbind(`(Intercept)` = 1, x = raw_r$x),
        dates = 1:4, control_group = cg, pt_assumption = pt,
        pt_start = if (pt == "all") 1 else NULL, base_period = "varying")
      ref_u <- cp06_pooled_public_or_reference(raw_r, route = route, intercept = TRUE,
        slope_mat = cbind(x = raw_r$x), Q_mat = cbind(`(Intercept)` = 1, x = raw_r$x),
        dates = 1:4, control_group = cg, pt_assumption = pt,
        pt_start = if (pt == "all") 1 else NULL, base_period = "universal")
      testthat::expect_equal(v$att, ref_v$att, tolerance = 1e-10)
      testthat::expect_equal(unname(v$inffunc), unname(ref_v$inffunc), tolerance = 1e-7)
      testthat::expect_equal(u$att, ref_u$att, tolerance = 1e-10)
      testthat::expect_equal(unname(u$inffunc), unname(ref_u$inffunc), tolerance = 1e-7)
      iv <- which(v$chained$role == "post"); iu <- which(u$chained$role == "post")
      testthat::expect_equal(v$att[iv], u$att[iu], tolerance = 1e-11)
      testthat::expect_equal(unname(v$inffunc[, iv]), unname(u$inffunc[, iu]), tolerance = 1e-10)
      for (fit in list(v, u)) {
        pre <- which(fit$chained$role == "pre"); record <- fit$chained$pretest
        testthat::expect_setequal(record$cells, fit$chained$target_id[pre])
        if (length(pre) && all(fit$chained$available[pre])) {
          V <- fit$V_analytical[pre, pre, drop = FALSE]
          F <- tryCatch(chol(V), error = function(e) NULL)
          if (!is.null(F) && rcond(F, triangular = TRUE) >= 1e-10) {
            expected <- fit$n * sum(fit$att[pre] * solve(V, fit$att[pre]))
            testthat::expect_identical(record$status, "computed")
            testthat::expect_equal(fit$W, expected, tolerance = 1e-9)
            testthat::expect_equal(record$df, length(pre))
          } else testthat::expect_identical(record$status, "unavailable")
        } else if (length(pre)) testthat::expect_identical(record$status, "unavailable")
        else testthat::expect_identical(record$status, "not_applicable")
      }
    }
  }
})

testthat::test_that("formula routing handles ~1, ~0+x, user-named `(Intercept)` covariates, and interval/calendar labels", {
  # ~1 intercept-only
  fx1 <- orp_fixture("none")
  fit1 <- orp_fit(fx1, xformla = ~1)
  ref1 <- cp06_pooled_public_or_reference(orp_raw, route = "none", intercept = TRUE,
    slope_mat = NULL, dates = 1:4, control_group = "nevertreated")
  testthat::expect_true(fit1$chained$outcome_design$intercept)
  testthat::expect_length(fit1$chained$outcome_design$slope_columns, 0L)
  testthat::expect_equal(fit1$att, ref1$att, tolerance = 1e-10)
  testthat::expect_equal(unname(fit1$inffunc), unname(ref1$inffunc), tolerance = 1e-7)

  # ~0 + x no-intercept raw uncentered shared-slope fitting across baseline and F/B/A
  for (route in c("none", "F", "B", "A")) {
    fx0 <- orp_fixture(route)
    raw_r <- orp_raw
    if (route %in% c("B", "A")) raw_r$rho <- fx0$frame$rho
    fit0 <- orp_fit(fx0, xformla = ~0 + x)
    ref0 <- cp06_pooled_public_or_reference(raw_r, route = route, intercept = FALSE,
      slope_mat = cbind(x = raw_r$x), Q_mat = cbind(`(Intercept)` = 1, x = raw_r$x),
      dates = 1:4, control_group = "nevertreated")
    diag0 <- fit0$chained$nuisance_diagnostics
    testthat::expect_false(fit0$chained$outcome_design$intercept)
    testthat::expect_identical(fit0$chained$outcome_design$slope_columns, "x")
    testthat::expect_false(diag0$coefficients$has_local_intercepts)
    testthat::expect_equal(diag0$accounting$n_local_intercepts, 0)
    testthat::expect_true(all(diag0$coefficients$local_intercepts == 0))
    testthat::expect_equal(diag0$coefficients$shared_slopes, ref0$beta, tolerance = 1e-10)
    testthat::expect_equal(fit0$att, ref0$att, tolerance = 1e-10)
    testthat::expect_equal(unname(fit0$inffunc), unname(ref0$inffunc), tolerance = 1e-7)
  }

  # User-named `(Intercept)` covariate when the formula has no intercept (~0 + `(Intercept)` and ~0 + `(Intercept)` + x)
  fx_icpt <- orp_fixture("none")
  fit_u1 <- orp_fit(fx_icpt, xformla = ~0 + `(Intercept)`)
  ref_u1 <- cp06_pooled_public_or_reference(orp_raw, route = "none", intercept = FALSE,
    slope_mat = cbind(`(Intercept)` = fx_icpt$frame[["(Intercept)"]]),
    dates = 1:4, control_group = "nevertreated")
  testthat::expect_false(fit_u1$chained$outcome_design$intercept)
  testthat::expect_identical(fit_u1$chained$outcome_design$slope_columns, "(Intercept)")
  testthat::expect_equal(fit_u1$chained$nuisance_diagnostics$coefficients$shared_slopes, ref_u1$beta, tolerance = 1e-10)
  testthat::expect_equal(fit_u1$att, ref_u1$att, tolerance = 1e-10)
  testthat::expect_equal(unname(fit_u1$inffunc), unname(ref_u1$inffunc), tolerance = 1e-7)

  fit_u2 <- orp_fit(fx_icpt, xformla = ~0 + `(Intercept)` + x)
  ref_u2 <- cp06_pooled_public_or_reference(orp_raw, route = "none", intercept = FALSE,
    slope_mat = cbind(`(Intercept)` = fx_icpt$frame[["(Intercept)"]], x = fx_icpt$frame$x),
    dates = 1:4, control_group = "nevertreated")
  testthat::expect_false(fit_u2$chained$outcome_design$intercept)
  testthat::expect_identical(fit_u2$chained$outcome_design$slope_columns, c("(Intercept)", "x"))
  testthat::expect_equal(fit_u2$chained$nuisance_diagnostics$coefficients$shared_slopes, ref_u2$beta, tolerance = 1e-10)
  testthat::expect_equal(fit_u2$att, ref_u2$att, tolerance = 1e-10)
  testthat::expect_equal(unname(fit_u2$inffunc), unname(ref_u2$inffunc), tolerance = 1e-7)

  # Interval/calendar labels: pair design uses stored-calendar transition count L_e = i_t - i_s, not raw t - s
  cal_dates <- c(10, 25, 40, 100)
  fx_cal <- orp_fixture("F", dates = cal_dates)
  testthat::expect_false(anyNA(fx_cal$frame$g))
  testthat::expect_equal(fx_cal$frame$g[orp_raw$cohort == 6],
    rep(max(cal_dates) + 60, sum(orp_raw$cohort == 6)))
  testthat::expect_gt(min(fx_cal$frame$g[orp_raw$cohort == 6]), max(cal_dates))
  fit_cal <- orp_fit(fx_cal, control_group = "notyettreated")
  fit_std <- orp_fit(orp_fixture("F", dates = 1:4), control_group = "notyettreated")
  diag_cal <- fit_cal$chained$nuisance_diagnostics
  testthat::expect_true(all(diag_cal$fitted_catalog$duration_steps %in% c(1, 2, 3)))
  testthat::expect_false(any((diag_cal$fitted_catalog$t - diag_cal$fitted_catalog$s) %in% c(1, 2, 3)))
  testthat::expect_equal(diag_cal$coefficients$shared_slopes,
    fit_std$chained$nuisance_diagnostics$coefficients$shared_slopes, tolerance = 1e-12)
  testthat::expect_equal(unname(diag_cal$coefficients$local_intercepts),
    unname(fit_std$chained$nuisance_diagnostics$coefficients$local_intercepts), tolerance = 1e-12)
  testthat::expect_equal(fit_cal$att, fit_std$att, tolerance = 1e-12)
  testthat::expect_equal(unname(fit_cal$inffunc), unname(fit_std$inffunc), tolerance = 1e-11)
})

testthat::test_that("local-deficient/global-supported pooled fit succeeds while full control deficiency and near-collinearity refuse", {
  # Local-deficient / global-supported case: x_local_def is constant (0) on never-treated (g == 0),
  # so blocks with control pool {Inf} are locally deficient, while blocks with control pool {4, Inf}
  # under NYT provide global support for the shared slope.
  fx_loc <- orp_fixture("none")
  fx_loc$args$frame$x_local_def <- ifelse(fx_loc$frame$g == 0, 0, fx_loc$frame$x)
  fit_loc <- orp_fit(fx_loc, xformla = ~x_local_def, control_group = "notyettreated")
  raw_loc <- orp_raw
  raw_loc$x_local_def <- fx_loc$args$frame$x_local_def
  ref_loc <- cp06_pooled_public_or_reference(raw_loc, route = "none", intercept = TRUE,
    slope_mat = cbind(x_local_def = raw_loc$x_local_def),
    dates = 1:4, control_group = "notyettreated")
  testthat::expect_equal(fit_loc$chained$nuisance_diagnostics$coefficients$shared_slopes, ref_loc$beta, tolerance = 1e-10)
  testthat::expect_equal(unname(fit_loc$chained$nuisance_diagnostics$coefficients$local_intercepts), ref_loc$alpha, tolerance = 1e-10)
  testthat::expect_equal(fit_loc$att, ref_loc$att, tolerance = 1e-10)
  testthat::expect_equal(unname(fit_loc$inffunc), unname(ref_loc$inffunc), tolerance = 1e-7)

  # Under nevertreated, the same covariate is globally rank-deficient on controls and varies on treated -> refused
  e <- testthat::expect_error(orp_fit(fx_loc, xformla = ~x_local_def, control_group = "nevertreated"),
    class = "cdid_or_wls_rank_error")
  testthat::expect_match(conditionMessage(e), "computation failure")
  testthat::expect_identical(e$diagnostics$request_context$input$outcome_design$retained_columns, c("(Intercept)", "x_local_def"))

  # Pooled-control alias whose null direction fails to reconstruct active treated rows -> refused
  fx_bad <- orp_fixture("none")
  fx_bad$args$frame$x_bad_tr <- ifelse(fx_bad$frame$g == 0, 2 * fx_bad$frame$x, fx_bad$frame$x^2 + 1)
  testthat::expect_error(orp_fit(fx_bad, xformla = ~x + x_bad_tr), class = "cdid_or_wls_rank_error")

  # Near-collinear design -> refused
  near <- orp_fixture(); near$args$frame$near <- near$frame$x + 1e-9 * (seq_len(nrow(near$frame)) %% 5 - 2)
  e <- testthat::expect_error(orp_fit(near, xformla = ~x + near), class = "cdid_or_wls_rank_error")
  testthat::expect_length(e$diagnostics$request_context$input$outcome_design$aliased_columns, 0L)

  # Empty outcome formulas remain refused on the public path
  testthat::expect_error(orp_fit(orp_fixture(), xformla = ~0), class = "cdid_chained_request_error")
  testthat::expect_error(orp_fit(orp_fixture(), xformla = ~-1), class = "cdid_chained_request_error")

  # Lazy estimator refusals (D48 CP07 activated est_method = "dr"; unsupported strings still refuse lazily)
  for (args in list(list(est_method = "stdipw"), list(est_method = "reg", selection = "custom"),
      list(est_method = function(...) stop("callback forced"), selection = "corollary1",
        selection_est_method = "imputation"))) {
    e <- testthat::expect_error(do.call(cdid::att_gt, c(list(yname = "y", tname = "time", idname = "id", gname = "g",
      data = quote(stop("data forced"))), args)), class = "cdid_chained_request_error")
    testthat::expect_no_match(conditionMessage(e), "data forced|callback forced")
  }
  for (method in list("dr", function(...) stop("callback forced"))) {
    testthat::expect_error(
      cdid::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
        data = stop("supported data reached"), est_method = method, selection = "corollary1"),
      "supported data reached"
    )
  }
})

testthat::test_that("original identity, clustering, multiplier and anticipation use existing consumers", {
  fx <- orp_fixture(); fit <- get("none", orp_state)
  clustered <- orp_fit(fx, clustervars = c("id", "cl"))
  testthat::expect_equal(clustered$att, fit$att, tolerance = 1e-12)
  testthat::expect_equal(clustered$inffunc, fit$inffunc, tolerance = 1e-12)
  available <- which(fit$chained$available)
  sums <- rowsum(fit$inffunc[, available, drop = FALSE], fx$frame$cl, reorder = FALSE)
  expected <- crossprod(scale(sums, scale = FALSE)) / fit$n^2
  testthat::expect_equal(unname(clustered$V_analytical[available, available] / clustered$n), unname(expected), tolerance = 1e-11)
  ord <- c(seq.int(2, nrow(fx$frame), 2), seq.int(1, nrow(fx$frame), 2))
  shuffled <- fx; shuffled$args$frame <- fx$frame[ord, ]; shuffled$args$data <- fx$data[nrow(fx$data):1, ]
  s <- orp_fit(shuffled)
  testthat::expect_equal(s$att, fit$att, tolerance = 1e-10)
  testthat::expect_identical(rownames(s$inffunc), fx$frame$id[ord])
  testthat::expect_equal(unname(s$inffunc), unname(fit$inffunc[ord, ]), tolerance = 1e-10)
  # Unknown external cohorts are preserved for baseline; linked panel cohorts
  # still define all regressions, and no new population-share assumption is made.
  unknown <- fx; unknown$args$frame <- fx$frame[, setdiff(names(fx$frame), "g")]
  un <- orp_fit(unknown)
  testthat::expect_equal(un$att, fit$att, tolerance = 1e-10)
  testthat::expect_equal(un$inffunc, fit$inffunc, tolerance = 1e-10)
  testthat::expect_true(anyNA(un$chained$frame$cohort))
  # Multiplier draws are inference from a fixed finite fixture, not a simulation.
  set.seed(1006); boot <- orp_fit(fx, bstrap = TRUE, cband = TRUE, biters = 199L)
  testthat::expect_equal(boot$att, fit$att, tolerance = 1e-12)
  testthat::expect_equal(boot$inffunc, fit$inffunc, tolerance = 1e-12)
  testthat::expect_true(boot$chained$inference$bstrap)
  testthat::expect_equal(boot$chained$inference$biters, 199L)
  for (route in c("F", "B", "A")) {
    set.seed(1006); z <- orp_fit(orp_fixture(route), bstrap = TRUE, cband = TRUE, biters = 199L)
    testthat::expect_equal(z$att, get(route, orp_state)$att, tolerance = 1e-12)
    testthat::expect_equal(z$inffunc, get(route, orp_state)$inffunc, tolerance = 1e-12)
    testthat::expect_true(z$chained$inference$bstrap)
  }
  anticipated <- orp_fit(fx, anticipation = 1L)
  testthat::expect_equal(anticipated$DIDparams$anticipation, 1)
  tm <- anticipated$chained$target_meta
  testthat::expect_equal(tm$base_t[tm$role == "post"], tm$g[tm$role == "post"] - 2)
  # Selection designs are unchanged: an outcome alias must not cure a singular logit.
  testthat::expect_error(orp_fit(orp_fixture("F"), selection_xformla = ~x + I(2 * x)), class = "cdid_error")
})

testthat::test_that("public print/summary/tidy/glance/plot labels identify control-mean OR", {
  fit <- get("none", orp_state)
  testthat::expect_match(paste(capture.output(print(fit)), collapse = "\n"), "outcome regression/control-mean", fixed = TRUE)
  testthat::expect_match(paste(capture.output(summary(fit)), collapse = "\n"), "outcome regression/control-mean", fixed = TRUE)
  testthat::expect_true(all(cdid::tidy(fit)$estimator == fit$chained$estimator))
  testthat::expect_identical(cdid::glance(fit)$estimator, fit$chained$estimator)
  testthat::expect_identical(cdid::ggdid(fit)$labels$caption, fit$chained$estimator)
  a <- cdid::aggte(fit, type = "group", bstrap = FALSE, cband = FALSE)
  testthat::expect_true(all(cdid::tidy(a)$estimator == a$chained$estimator))
  testthat::expect_identical(cdid::glance(a)$estimator, a$chained$estimator)
  testthat::expect_identical(cdid::ggdid(a)$labels$caption, a$chained$estimator)
})

testthat::test_that("public OR resource planning preserves explicit chunks and finite defaults", {
  fx <- orp_fixture("F")
  automatic <- orp_fit(fx)
  fixed <- orp_fit(fx, control = list(chunk_rows = 7L))
  testthat::expect_equal(automatic$att, fixed$att, tolerance = 1e-11)
  testthat::expect_equal(automatic$inffunc, fixed$inffunc, tolerance = 1e-10)
  testthat::expect_equal(automatic$V_analytical, fixed$V_analytical, tolerance = 1e-11)
  testthat::expect_equal(automatic$chained$nuisance_diagnostics$accounting$automatic_chunk_rows, TRUE)
  testthat::expect_equal(fixed$chained$nuisance_diagnostics$accounting$automatic_chunk_rows, FALSE)
  testthat::expect_equal(fixed$chained$resource_accounting$point_producer$chunk_rows, 7, tolerance = 0)
  profile <- orp_ns("chained_public_profile_defaults")()
  testthat::expect_equal(profile$numeric_workspace_bytes, 1073741824, tolerance = 0)
  testthat::expect_identical(profile$counted_work_limit, Inf)
  testthat::expect_equal(profile$workers, 1, tolerance = 0)
  # Existing baseline route has the same complete IF under explicit partitions.
  baseline <- orp_fit(orp_fixture(), control = list(chunk_rows = 7L))
  testthat::expect_equal(baseline$att, get("none", orp_state)$att, tolerance = 1e-11)
  testthat::expect_equal(baseline$inffunc, get("none", orp_state)$inffunc, tolerance = 1e-10)
  testthat::expect_equal(baseline$V_analytical, get("none", orp_state)$V_analytical, tolerance = 1e-11)
})

testthat::test_that("automatic OR row planning reserves the final diagnostic copy", {
  fit <- orp_fit(orp_fixture("F"))
  resource <- fit$chained$resource_accounting
  owned <- fit$chained$nuisance_diagnostics$accounting$owned_entries
  bound <- 4 * owned + 4096
  testthat::expect_equal(resource$final_diagnostic_copy_entries, 5 * bound, tolerance = 0)
  testthat::expect_lte(max(resource$simultaneous_entries, resource$initial_envelope_entries) +
    resource$final_diagnostic_copy_entries, fit$chained$public$numeric_workspace_entries)
  baseline <- orp_fit(orp_fixture())
  testthat::expect_equal(baseline$chained$resource_accounting$final_diagnostic_copy_entries, 0, tolerance = 0)
})
