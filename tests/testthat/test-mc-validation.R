# Strict installed-namespace tests for private R/mc_validation.R under D7/D14/D15/D16.
# Requires an actual installed locked cdid namespace at file load with scoped
# formal/body/environment bindings and missing-namespace refusal active.
# No source(), sys.source(), or dyn.load() fallback for R/mc_validation.R is permitted.

verify_strict_cdid_namespace <- function(pkg = "cdid", env_override = NULL) {
  ns <- if (is.null(env_override)) {
    if (!isNamespaceLoaded(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Strict packaged tests require an actual installed locked `%s` namespace at file load; missing-namespace refusal active (no source/dyn.load fallback).",
                   pkg), call. = FALSE)
    }
    asNamespace(pkg)
  } else {
    env_override
  }

  if (!isNamespace(ns) || !identical(unname(getNamespaceName(ns)), "cdid")) {
    stop("Target environment is not the `cdid` package namespace.", call. = FALSE)
  }
  if (!identical(unname(getNamespaceVersion(ns)), "1.1.0.0")) {
    stop("`cdid` namespace version must be 1.1.0.0.", call. = FALSE)
  }
  if (!environmentIsLocked(ns) || !environmentIsLocked(parent.env(ns))) {
    stop("`cdid` namespace and its imports environment must be locked.",
         call. = FALSE)
  }

  pkg_path <- getNamespaceInfo(ns, "path")
  if (!is.character(pkg_path) || length(pkg_path) != 1L || !nzchar(pkg_path) ||
      !file.exists(file.path(pkg_path, "DESCRIPTION")) ||
      !file.exists(file.path(pkg_path, "NAMESPACE"))) {
    stop("`cdid` namespace path must contain valid package metadata (DESCRIPTION and NAMESPACE).",
         call. = FALSE)
  }

  expected_exports <- c("aggte", "att_gt", "conditional_did_pretest",
                        "ggdid", "glance", "tidy")
  actual_exports <- getNamespaceExports(ns)
  if (!setequal(actual_exports, expected_exports)) {
    stop("`cdid` namespace exports do not match the immutable NAMESPACE specification.",
         call. = FALSE)
  }

  private_fns <- list(
    mc_allowed_metrics = character(),
    mc_allowed_looks = character(),
    mc_is_whole_number = c("x"),
    mc_is_scalar_flag = c("x"),
    mc_is_scalar_string = c("x"),
    mc_is_scalar_finite_number = c("x"),
    mc_clopper_pearson_ci = c("k", "d", "alpha_tail"),
    mc_wilson_ci = c("k", "d", "conf_level"),
    mc_classify_interval_decision = c("metric", "lower", "upper", "nominal",
                                      "supported", "unsupported_reasons"),
    mc_continuous_bias_ci = c("t_errors", "n_attempted", "premises_supported",
                              "unsupported_reasons", "z"),
    mc_continuous_se_ratio_ci = c("t_errors", "se_values", "n_attempted",
                                  "premises_supported", "unsupported_reasons", "z"),
    mc_resolve_premise_flag = c("premises", "specific_key", "fallback_key"),
    mc_check_metric_premises = c("metric", "specification"),
    validate_mc_specification = c("specification", "look"),
    validate_mc_records = c("records", "specification", "look"),
    mc_evaluate_coverage_metric = c("records", "specification", "metric_name"),
    mc_evaluate_null_size_metric = c("records", "specification"),
    mc_nonfinite_error_metric = c("metric", "m", "n_attempted", "premises", "n_zero_se"),
    mc_evaluate_single_metric = c("metric", "records", "specification"),
    mc_compute_cell_status = c("metric_statuses"),
    mc_summarize_failures_and_events = c("records", "specification"),
    mc_specifications_compatible = c("spec_2000", "spec_8000"),
    mc_records_prefix_identical = c("prefix_df", "full_df"),
    evaluate_private_mc_cell = c("records", "specification", "look",
                                 "prior_assessment"),
    mc_calibration_outer_seed = c("cell_index", "outer_index"),
    mc_calibration_continuation_decision = c("assessment_2000"),
    prepare_mc_calibration_plumbing = character()
  )

  if (any(names(private_fns) %in% actual_exports)) {
    stop("Private `R/mc_validation.R` symbols must never be exported by `cdid`.",
         call. = FALSE)
  }

  for (fn_name in names(private_fns)) {
    if (!exists(fn_name, envir = ns, inherits = FALSE)) {
      stop(sprintf("Required private function `%s` is missing from `cdid` namespace.",
                   fn_name), call. = FALSE)
    }
    if (!bindingIsLocked(fn_name, ns)) {
      stop(sprintf("Binding `%s` in `cdid` namespace must be locked.", fn_name),
           call. = FALSE)
    }
    fn <- get(fn_name, envir = ns, inherits = FALSE)
    if (typeof(fn) != "closure" || !is.function(fn)) {
      stop(sprintf("`%s` in `cdid` namespace must be a closure.", fn_name),
           call. = FALSE)
    }
    if (!identical(environment(fn), ns)) {
      stop(sprintf("Closure `%s` environment must be identical to the `cdid` namespace.",
                   fn_name), call. = FALSE)
    }
    actual_formals <- names(formals(fn))
    expected_formals <- private_fns[[fn_name]]
    if (length(expected_formals) == 0L) {
      if (!is.null(actual_formals)) {
        stop(sprintf("`%s` must take zero arguments.", fn_name), call. = FALSE)
      }
    } else if (!identical(actual_formals, expected_formals)) {
      stop(sprintf("Formals of `%s` do not match expected signature.", fn_name),
           call. = FALSE)
    }
  }

  eval_fn <- get("evaluate_private_mc_cell", envir = ns, inherits = FALSE)
  if (!is.null(formals(eval_fn)$prior_assessment)) {
    stop("`evaluate_private_mc_cell` default for `prior_assessment` must be NULL.",
         call. = FALSE)
  }

  # Verify scoped body AST against R/mc_validation.R when source tree is readable
  src_candidates <- c("R/mc_validation.R", file.path("..", "..", "R", "mc_validation.R"))
  src_hit <- src_candidates[file.exists(src_candidates)]
  if (length(src_hit) > 0L) {
    parsed_env <- new.env(parent = baseenv())
    sys.source(src_hit[[1L]], envir = parsed_env, keep.source = FALSE)
    for (fn_name in names(private_fns)) {
      ns_fn <- get(fn_name, envir = ns, inherits = FALSE)
      src_fn <- get(fn_name, envir = parsed_env, inherits = FALSE)
      if (!identical(formals(ns_fn), formals(src_fn)) ||
          !identical(deparse(body(ns_fn)), deparse(body(src_fn)))) {
        stop(sprintf("Installed namespace function `%s` formals/body diverge from `R/mc_validation.R`.",
                     fn_name), call. = FALSE)
      }
    }
  }

  ns
}

cdid_ns <- verify_strict_cdid_namespace("cdid")

ref_fixture_candidates <- c(
  "fixtures/mc_validation_reference.R",
  file.path("tests", "testthat", "fixtures", "mc_validation_reference.R")
)
ref_fixture_hit <- ref_fixture_candidates[file.exists(ref_fixture_candidates)]
if (length(ref_fixture_hit) == 0L) {
  stop("Fixture file `mc_validation_reference.R` not found.", call. = FALSE)
}
sys.source(ref_fixture_hit[[1L]], envir = environment())

# Independent fixture arithmetic, including the conditional joint subset and
# fixed critical 4.0. This uses the unchanged external reference calculator,
# never the package's interval or classification functions.
cp02_mc_expect_ratio_reference <- function(metric, records, specification) {
  joint <- records$point_computable & records$se_computable
  ref <- ref_sample_moments(records$estimate[joint] - specification$truth, records$se[joint])
  expected_status <- if (ref$ci_r[["lower"]] >= 0.95 && ref$ci_r[["upper"]] <= 1.05) "PASS" else
    if (ref$ci_r[["upper"]] < 0.95 || ref$ci_r[["lower"]] > 1.05) "FAIL" else "INCONCLUSIVE"
  expect_identical(metric$z, 4.0)
  expect_equal(metric$estimate, ref$rho_hat, tolerance = 1e-15)
  expect_equal(metric$emp_sd, ref$sd, tolerance = 1e-15)
  expect_equal(metric$if_var, ref$vr, tolerance = 1e-15)
  expect_equal(metric$mcse_log, ref$mcse_r, tolerance = 1e-15)
  expect_equal(metric$interval, ref$ci_r, tolerance = 1e-15)
  expect_identical(metric$status, expected_status)
  ref
}

evaluate_private_mc_cell <- get("evaluate_private_mc_cell", envir = cdid_ns, inherits = FALSE)
mc_clopper_pearson_ci <- get("mc_clopper_pearson_ci", envir = cdid_ns, inherits = FALSE)
mc_wilson_ci <- get("mc_wilson_ci", envir = cdid_ns, inherits = FALSE)
mc_classify_interval_decision <- get("mc_classify_interval_decision", envir = cdid_ns, inherits = FALSE)
mc_continuous_bias_ci <- get("mc_continuous_bias_ci", envir = cdid_ns, inherits = FALSE)
mc_continuous_se_ratio_ci <- get("mc_continuous_se_ratio_ci", envir = cdid_ns, inherits = FALSE)
prepare_mc_calibration_plumbing <- get("prepare_mc_calibration_plumbing", envir = cdid_ns, inherits = FALSE)
mc_calibration_outer_seed <- get("mc_calibration_outer_seed", envir = cdid_ns, inherits = FALSE)
mc_calibration_continuation_decision <- get("mc_calibration_continuation_decision", envir = cdid_ns, inherits = FALSE)

test_that("strict locked cdid namespace binding and missing-namespace refusal work", {
  expect_true(isNamespace(cdid_ns))
  expect_identical(unname(getNamespaceName(cdid_ns)), "cdid")
  expect_true(environmentIsLocked(cdid_ns))
  expect_false("evaluate_private_mc_cell" %in% getNamespaceExports(cdid_ns))
  expect_identical(environment(evaluate_private_mc_cell), cdid_ns)

  # Refusal on missing package
  expect_error(
    verify_strict_cdid_namespace("cdid_missing_pkg_999999"),
    "missing-namespace refusal active"
  )
  # Refusal on wrong installed namespace (stats)
  expect_error(
    verify_strict_cdid_namespace(env_override = asNamespace("stats")),
    "not the `cdid` package namespace"
  )
  # Refusal on unlocked environment
  fake_env <- new.env(parent = emptyenv())
  expect_error(
    verify_strict_cdid_namespace(env_override = fake_env),
    "not the `cdid` package namespace"
  )
})

test_that("evaluator and plumbing touch zero RNG state and perform no draws or fits", {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  seed_before <- if (had_seed) get(".Random.seed", envir = genv, inherits = FALSE) else NULL

  spec <- make_mc_spec(n_seeds = 2000L)
  rec <- make_mc_records(spec, look = 2000L)
  res <- evaluate_private_mc_cell(rec, spec, look = 2000L)
  plumb <- prepare_mc_calibration_plumbing()

  expect_identical(res$look, 2000L)
  expect_identical(plumb$outer_draws_executed, 0L)
  expect_identical(plumb$fits_executed, 0L)
  expect_identical(exists(".Random.seed", envir = genv, inherits = FALSE), had_seed)
  if (had_seed) {
    expect_identical(get(".Random.seed", envir = genv, inherits = FALSE), seed_before)
  }
})

test_that("frozen 91 exact checks, 17 algebraic controls, and root replay 1217 match", {
  frozen_path <- find_project_file("reference/approved/exact_output_frozen.json")
  replay_path <- find_project_file("reference/root_replay_1217/exact_output.json")
  receipt_path <- find_project_file("reference/root_replay_1217/receipt.json")
  expect_false(is.na(frozen_path))
  expect_false(is.na(replay_path))
  expect_false(is.na(receipt_path))

  frozen_bytes <- readBin(frozen_path, what = "raw", n = file.info(frozen_path)$size)
  replay_bytes <- readBin(replay_path, what = "raw", n = file.info(replay_path)$size)
  expect_identical(frozen_bytes, replay_bytes)

  frozen_txt <- paste(readLines(frozen_path, warn = FALSE), collapse = "\n")
  expect_match(frozen_txt, '"exact_checks": 91', fixed = TRUE)
  expect_match(frozen_txt, '"rejecting_controls": 17', fixed = TRUE)
  expect_match(frozen_txt, '"fits": 0', fixed = TRUE)
  expect_match(frozen_txt, '"outer_calibration_draws": 0', fixed = TRUE)
  expect_match(frozen_txt, '"status": "PASS"', fixed = TRUE)

  # Exact three-row sample moments witness from exact_reference.py:
  # rows = [(-4/5, 1), (1/5, 2), (6/5, 3)]
  t3 <- c(-4 / 5, 1 / 5, 6 / 5)
  s3 <- c(1, 2, 3)
  ref3 <- ref_sample_moments(t3, s3)
  expect_equal(ref3$mu, 1 / 5, tolerance = 1e-15)
  expect_equal(ref3$v, 1, tolerance = 1e-15)
  expect_equal(ref3$nu, 2, tolerance = 1e-15)
  expect_equal(ref3$pb_bar, 1 / 30, tolerance = 1e-15)
  expect_equal(ref3$mcse2_b, 301 / 900, tolerance = 1e-15)
  expect_equal(ref3$pr_bar, 1 / 6, tolerance = 1e-15)
  expect_equal(ref3$mcse2_r, 1 / 9, tolerance = 1e-15)

  b3 <- mc_continuous_bias_ci(t3, n_attempted = 2000L)
  r3 <- mc_continuous_se_ratio_ci(t3, s3, n_attempted = 2000L)
  r3_legacy <- mc_continuous_se_ratio_ci(t3, s3, n_attempted = 2000L,
                                         z = stats::qnorm(399 / 400))
  expect_equal(b3$raw_bias, 1 / 5, tolerance = 1e-15)
  expect_equal(b3$emp_var, 1, tolerance = 1e-15)
  expect_equal(b3$emp_sd, 1, tolerance = 1e-15)
  expect_equal(b3$estimate, 1 / 5, tolerance = 1e-15)
  expect_equal(b3$psi_bar, 1 / 30, tolerance = 1e-15)
  expect_equal(b3$if_var, 301 / 300, tolerance = 1e-15)
  expect_equal(b3$mcse^2, 301 / 900, tolerance = 1e-15)
  expect_equal(b3$z, stats::qnorm(399 / 400), tolerance = 1e-15)
  expect_equal(b3$interval, ref3$ci_b, tolerance = 1e-15)
  expect_identical(b3$denominator, "3/2000")

  expect_equal(r3$mean_se, 2, tolerance = 1e-15)
  expect_equal(r3$emp_var, 1, tolerance = 1e-15)
  expect_equal(r3$emp_sd, 1, tolerance = 1e-15)
  expect_equal(r3$estimate, 2, tolerance = 1e-15)
  expect_equal(r3$psi_bar, 1 / 6, tolerance = 1e-15)
  expect_equal(r3$if_var, 1 / 3, tolerance = 1e-15)
  expect_equal(r3$mcse_log^2, 1 / 9, tolerance = 1e-15)
  expect_identical(r3$z, 4.0)
  expect_equal(r3$interval, ref3$ci_r, tolerance = 1e-15)
  expect_equal(r3_legacy$z, stats::qnorm(399 / 400), tolerance = 1e-15)
  expect_equal(r3_legacy$interval, ref3$ci_r_legacy_qnorm, tolerance = 1e-15)
  expect_identical(r3$denominator, "3/2000")

  # Also verify through full evaluate_private_mc_cell at look = 2000 with 3 computable rows
  # and a 4th point-only row T = 9 (se_computable = FALSE):
  # Point-only T = c(-1, 0, 1, 9) has variance 251/12, while joint 3 rows have T = c(-4/5, 1/5, 6/5)
  # or T = c(-1, 0, 1) with variance 1.
  spec_sub <- make_mc_spec(n_seeds = 2000L, truth = 0,
                           required_metrics = c("bias", "se_ratio"))
  rec_sub <- make_mc_records(spec_sub, look = 2000L)
  rec_sub$point_computable <- rep(FALSE, 2000L)
  rec_sub$se_computable <- rep(FALSE, 2000L)
  rec_sub$interval_computable <- rep(FALSE, 2000L)
  rec_sub$band_computable <- rep(FALSE, 2000L)
  rec_sub$test_computable <- rep(FALSE, 2000L)
  rec_sub$failure_type <- rep("estimation_failed", 2000L)
  rec_sub$estimate <- rep(999, 2000L)
  rec_sub$se <- rep(999, 2000L)

  # First 3 rows are both point and SE computable with exact 3-row witness
  rec_sub$point_computable[1:3] <- TRUE
  rec_sub$se_computable[1:3] <- TRUE
  rec_sub$estimate[1:3] <- t3
  rec_sub$se[1:3] <- s3
  rec_sub$failure_type[1:3] <- "inference_failed"

  res_3row <- evaluate_private_mc_cell(rec_sub, spec_sub, look = 2000L)
  expect_equal(res_3row$metrics$bias$psi_bar, 1 / 30, tolerance = 1e-15)
  expect_equal(res_3row$metrics$bias$mcse^2, 301 / 900, tolerance = 1e-15)
  expect_equal(res_3row$metrics$se_ratio$psi_bar, 1 / 6, tolerance = 1e-15)
  expect_equal(res_3row$metrics$se_ratio$mcse_log^2, 1 / 9, tolerance = 1e-15)

  # Now set 1:3 to T = (-1, 0, 1), S = (1, 2, 3) and row 4 to point-only T = 9 (se_computable = FALSE)
  rec_sub$estimate[1:3] <- c(-1, 0, 1)
  rec_sub$se[1:3] <- c(1, 2, 3)
  rec_sub$point_computable[4] <- TRUE
  rec_sub$se_computable[4] <- FALSE
  rec_sub$estimate[4] <- 9
  rec_sub$se[4] <- NA_real_
  rec_sub$failure_type[4] <- "inference_failed"

  res_pt_only <- evaluate_private_mc_cell(rec_sub, spec_sub, look = 2000L)
  expect_identical(res_pt_only$metrics$bias$m_computable, 4L)
  expect_identical(res_pt_only$metrics$bias$denominator, "4/2000")
  expect_equal(res_pt_only$metrics$bias$emp_var, 251 / 12, tolerance = 1e-14)
  expect_identical(res_pt_only$metrics$se_ratio$m_computable, 3L)
  expect_identical(res_pt_only$metrics$se_ratio$denominator, "3/2000")
  expect_equal(res_pt_only$metrics$se_ratio$emp_var, 1, tolerance = 1e-15)
  expect_equal(res_pt_only$metrics$se_ratio$mean_se, 2, tolerance = 1e-15)
  expect_equal(res_pt_only$metrics$se_ratio$estimate, 2, tolerance = 1e-15)
})

test_that("skewed correlated atom law, population Jacobian/polynomial derivatives, and omit-derivative discriminators hold", {
  # Z = (-1, 0, 2) with probabilities (1/3, 1/2, 1/6), S = 1 + (Z^2 - 1)/4
  zs <- c(-1, 0, 2)
  qs <- c(1 / 3, 1 / 2, 1 / 6)
  expect_equal(
    vapply(1:4, function(j) sum(qs * zs^j), numeric(1L)),
    c(0, 1, 1, 3),
    tolerance = 1e-15
  )

  # Population complete frequency Jacobian and polynomial quotient derivative checks
  expected_derivs <- list(
    `0_1` = c(db = 3, dr = -3 / 4),
    `0_0` = c(db = 0, dr = 0),
    `0.2_1` = c(db = 11 / 4, dr = -3 / 4),
    `0.2_0` = c(db = 0, dr = 0),
    `-0.2_1` = c(db = 13 / 4, dr = -3 / 4),
    `-0.2_0` = c(db = 0, dr = 0)
  )
  for (bias_val in c(0, 1 / 5, -1 / 5)) {
    pi_val <- 1 / 2
    mu_val <- bias_val
    var_val <- 1
    nu_val <- 1
    p <- pi_val
    a <- pi_val * mu_val
    b <- pi_val * (var_val + mu_val^2)
    c_tot <- pi_val * nu_val

    for (d_flag in c(1, 0)) {
      t_atom <- if (d_flag == 1) bias_val + 3 / 2 else 99
      s_atom <- if (d_flag == 1) 5 / 4 else 99
      psi_b <- (t_atom - mu_val) - mu_val * ((t_atom - mu_val)^2 - 1) / 2
      psi_r <- (s_atom - nu_val) / nu_val - ((t_atom - mu_val)^2 - 1) / 2
      db <- d_flag * psi_b / pi_val
      dr <- d_flag * psi_r / pi_val
      key <- paste0(bias_val, "_", d_flag)
      expect_equal(c(db = db, dr = dr), expected_derivs[[key]], tolerance = 1e-15)

      # Second independent polynomial quotient derivation
      dp <- d_flag - p
      da <- d_flag * t_atom - a
      db2 <- d_flag * t_atom^2 - b
      dc <- d_flag * s_atom - c_tot
      n0 <- c_tot^2
      n1 <- 2 * c_tot * dc
      n2 <- dc^2
      d0 <- b * p - a^2
      d1 <- db2 * p + b * dp - 2 * a * da
      d2 <- db2 * dp - da^2
      big_a <- n1 * d0 - n0 * d1
      expect_equal(dr, big_a / (2 * n0 * d0), tolerance = 1e-15)
      if (bias_val != 0) {
        expect_equal(db, (2 * a * da * d0 - a^2 * d1) / (2 * bias_val * d0^2), tolerance = 1e-15)
      }
      for (h in c(1 / 10, 1 / 100, 1 / 1000)) {
        den <- d0 + d1 * h + d2 * h^2
        gh <- (n0 + n1 * h + n2 * h^2) / den
        rem <- h * ((n2 * d0 - n0 * d2) * d0 - big_a * d1 - big_a * d2 * h) / (d0^2 * den)
        expect_true(abs((gh - n0 / d0) - h * (big_a / d0^2) - h * rem) < 1e-14)
        expect_true(abs((gh - n0 / d0) / h - big_a / d0^2 - rem) < 1e-11)
      }
    }

    pb_pop <- zs - bias_val * (zs^2 - 1) / 2
    pr_pop <- -(zs^2 - 1) / 4
    expect_equal(sum(qs * pb_pop), 0, tolerance = 1e-15)
    expect_equal(sum(qs * pr_pop), 0, tolerance = 1e-15)
    vb_pop <- sum(qs * pb_pop^2)
    vr_pop <- sum(qs * pr_pop^2)
    expect_equal(vb_pop, 1 + bias_val^2 / 2 - bias_val, tolerance = 1e-15)
    expect_equal(vr_pop, 1 / 8, tolerance = 1e-15)
    expect_equal(2 * vb_pop / 2000, vb_pop / 1000, tolerance = 1e-15)
  }

  # Repeat the exact 6-atom cycle (2 copies of -1, 3 copies of 0, 1 copy of 2)
  # across m = 1200 computable rows out of N = 2000 attempts.
  z6 <- c(-1, -1, 0, 0, 0, 2)
  s6 <- 1 + (z6^2 - 1) / 4
  z_rep <- rep(z6, 200L)
  s_rep <- rep(s6, 200L)

  for (bias_val in c(0, 1 / 5, -1 / 5)) {
    t_rep <- z_rep + bias_val
    b_out <- mc_continuous_bias_ci(t_rep, n_attempted = 2000L)
    r_out <- mc_continuous_se_ratio_ci(t_rep, s_rep, n_attempted = 2000L)
    r_out_legacy <- mc_continuous_se_ratio_ci(t_rep, s_rep, n_attempted = 2000L,
                                              z = stats::qnorm(399 / 400))

    # Sample variance uses m - 1 = 1199, while population variance on the 6-atom cycle is 1
    m <- length(t_rep)
    s2_emp <- m / (m - 1)
    # Verify exact agreement with independent reference formula (which includes SD derivative in beta
    # and both mean-SE and SD derivatives with covariance/skewness in log-rho)
    ref_out <- ref_sample_moments(t_rep, s_rep)
    expect_equal(b_out$estimate, ref_out$beta_hat, tolerance = 1e-15)
    expect_equal(b_out$if_var, ref_out$vb, tolerance = 1e-15)
    expect_equal(b_out$mcse, ref_out$mcse_b, tolerance = 1e-15)
    expect_equal(b_out$z, ref_out$z_b, tolerance = 1e-15)
    expect_equal(b_out$interval, ref_out$ci_b, tolerance = 1e-15)

    expect_equal(r_out$estimate, ref_out$rho_hat, tolerance = 1e-15)
    expect_equal(r_out$if_var, ref_out$vr, tolerance = 1e-15)
    expect_equal(r_out$mcse_log, ref_out$mcse_r, tolerance = 1e-15)
    expect_identical(r_out$z, 4.0)
    expect_equal(r_out$interval, ref_out$ci_r, tolerance = 1e-15)
    expect_equal(r_out_legacy$z, ref_out$z_b, tolerance = 1e-15)
    expect_equal(r_out_legacy$interval, ref_out$ci_r_legacy_qnorm, tolerance = 1e-15)

    # Discriminate against wrong formulas:
    # 1) Omitting SD derivative in psi_beta (using (T - Tbar)/s only, whose variance is 1):
    if (bias_val != 0) {
      wrong_psi_b <- (t_rep - mean(t_rep)) / sqrt(s2_emp)
      wrong_vb <- stats::var(wrong_psi_b)
      expect_true(abs(b_out$if_var - wrong_vb) > 0.05)
    }
    # 2) Omitting mean-SE derivative in psi_logrho (using -((T - Tbar)^2 - s2)/(2*s2) only):
    wrong_psi_r_no_se <- -((t_rep - mean(t_rep))^2 - s2_emp) / (2 * s2_emp)
    expect_true(abs(r_out$if_var - stats::var(wrong_psi_r_no_se)) > 0.10)
    # 3) Normal-kurtosis / independent shortcut Var(psi_logrho) = 1/2 instead of ~ 1/8:
    expect_true(abs(r_out$if_var - 0.5) > 0.30)
    # 4) Using attempted N = 2000 instead of conditional m = 1200 in MCSE:
    expect_true(abs(r_out$mcse_log - sqrt(r_out$if_var / 2000)) > 1e-3)
  }
})

test_that("all 5 degenerate/undefined sample moment controls return INCONCLUSIVE", {
  # 1) m = 0
  b0 <- mc_continuous_bias_ci(numeric(0), n_attempted = 2000L)
  r0 <- mc_continuous_se_ratio_ci(numeric(0), numeric(0), n_attempted = 2000L)
  expect_identical(b0$status, "INCONCLUSIVE")
  expect_false(b0$supported)
  expect_true("m_lt_2" %in% b0$unsupported_reasons)
  expect_identical(r0$status, "INCONCLUSIVE")
  expect_false(r0$supported)
  expect_true("m_lt_2" %in% r0$unsupported_reasons)

  # 2) m = 1
  b1 <- mc_continuous_bias_ci(1, n_attempted = 2000L)
  r1 <- mc_continuous_se_ratio_ci(1, 1, n_attempted = 2000L)
  expect_identical(b1$status, "INCONCLUSIVE")
  expect_true("m_lt_2" %in% b1$unsupported_reasons)
  expect_identical(r1$status, "INCONCLUSIVE")
  expect_true("m_lt_2" %in% r1$unsupported_reasons)

  # 3) Zero empirical SD: [(1, 1), (1, 2)]
  b_zsd <- mc_continuous_bias_ci(c(1, 1), n_attempted = 2000L)
  r_zsd <- mc_continuous_se_ratio_ci(c(1, 1), c(1, 2), n_attempted = 2000L)
  expect_identical(b_zsd$status, "INCONCLUSIVE")
  expect_true("zero_or_nonpositive_empirical_sd" %in% b_zsd$unsupported_reasons)
  expect_identical(r_zsd$status, "INCONCLUSIVE")
  expect_true("zero_or_nonpositive_empirical_sd" %in% r_zsd$unsupported_reasons)

  # 4) Zero mean SE: [(-1, 0), (0, 0), (1, 0)]
  # Individual zero SEs are valid inputs (not rejected as malformed),
  # bias is computable, and se_ratio is INCONCLUSIVE due to nonpositive mean SE.
  b_zse <- mc_continuous_bias_ci(c(-1, 0, 1), n_attempted = 2000L)
  r_zse <- mc_continuous_se_ratio_ci(c(-1, 0, 1), c(0, 0, 0), n_attempted = 2000L)
  expect_true(b_zse$supported)
  expect_identical(r_zse$status, "INCONCLUSIVE")
  expect_false(r_zse$supported)
  expect_identical(r_zse$n_zero_se, 3L)
  expect_true("nonpositive_or_nonfinite_mean_se" %in% r_zse$unsupported_reasons)

  # 5) Degenerate plug-in IF variance: [(-1, 7), (0, 4), (1, 7)]
  # Here psi_logrho = (1/6, 1/6, 1/6) is constant, so centered IF variance is 0!
  r_deg <- mc_continuous_se_ratio_ci(c(-1, 0, 1), c(7, 4, 7), n_attempted = 2000L)
  expect_identical(r_deg$status, "INCONCLUSIVE")
  expect_false(r_deg$supported)
  expect_equal(r_deg$psi_bar, 1 / 6, tolerance = 1e-15)
  expect_equal(r_deg$if_var, 0, tolerance = 1e-15)
  expect_true("zero_or_nonpositive_influence_variance" %in% r_deg$unsupported_reasons)

  # Also degenerate bias IF variance when T = (-1, 1):
  b_deg <- mc_continuous_bias_ci(c(-1, 1), n_attempted = 2000L)
  # Wait: for m=2, T=(-1, 1), Tbar=0, s^2=2, s=sqrt(2), psi_beta = (-1/sqrt(2), 1/sqrt(2)), v_hat = 1 > 0;
  # whereas for S=(1, 1), T=(-1, 1), Sbar=1, s^2=2, psi_logrho = (1/4, 1/4) -> constant!
  r_deg2 <- mc_continuous_se_ratio_ci(c(-1, 1), c(1, 1), n_attempted = 2000L)
  expect_identical(r_deg2$status, "INCONCLUSIVE")
  expect_equal(r_deg2$if_var, 0, tolerance = 1e-15)
  expect_true("zero_or_nonpositive_influence_variance" %in% r_deg2$unsupported_reasons)
})

test_that("finite individual zero SE is retained in evaluate_private_mc_cell when mean SE > 0", {
  spec <- make_mc_spec(n_seeds = 2000L, required_metrics = c("bias", "se_ratio"))
  rec <- make_mc_records(spec, look = 2000L)
  # Set 50 rows to have finite zero reported SE while se_computable == TRUE
  rec$se[1:50] <- 0.0
  # Rescale remaining SEs so overall mean SE equals empirical SD
  emp_sd_all <- stats::sd(rec$estimate - spec$truth)
  rec$se[51:2000] <- rec$se[51:2000] * (emp_sd_all * 2000 / sum(rec$se[51:2000]))

  res <- evaluate_private_mc_cell(rec, spec, look = 2000L)
  expect_identical(res$metrics$se_ratio$m_computable, 2000L)
  expect_identical(res$metrics$se_ratio$n_zero_se, 50L)
  expect_true(res$metrics$se_ratio$supported)
  expect_equal(res$metrics$se_ratio$estimate, 1.0, tolerance = 1e-12)
  expect_identical(res$metrics$se_ratio$status, "INCONCLUSIVE")
  ref_zero <- cp02_mc_expect_ratio_reference(res$metrics$se_ratio, rec, spec)
  # Retain the old vectors: their qnorm interval passed, but critical 4.0
  # crosses the same unchanged D7 margins. Zero SEs remain computable.
  expect_true(ref_zero$ci_r_legacy_qnorm[["lower"]] >= 0.95 &&
    ref_zero$ci_r_legacy_qnorm[["upper"]] <= 1.05)
  expect_true(ref_zero$ci_r[["lower"]] < 0.95 || ref_zero$ci_r[["upper"]] > 1.05)

  # Now if ALL 2000 computable SEs are 0, record validation still succeeds (0 SE is valid input),
  # and se_ratio returns INCONCLUSIVE due to nonpositive mean SE
  rec_all_zero_se <- rec
  rec_all_zero_se$se <- rep(0.0, 2000L)
  res_zero <- evaluate_private_mc_cell(rec_all_zero_se, spec, look = 2000L)
  expect_identical(res_zero$metrics$se_ratio$status, "INCONCLUSIVE")
  expect_identical(res_zero$metrics$se_ratio$n_zero_se, 2000L)
  expect_true("nonpositive_or_nonfinite_mean_se" %in% res_zero$metrics$se_ratio$unsupported_reasons)
})

test_that("Clopper-Pearson 99.5% boundaries, binomial tail recurrences, and pinned stats match", {
  # d = 0 is undefined / NA
  expect_identical(mc_clopper_pearson_ci(0L, 0L), c(lower = NA_real_, upper = NA_real_))

  # d = 1 boundaries
  cp_1_0 <- mc_clopper_pearson_ci(0L, 1L, alpha_tail = 1 / 400)
  cp_1_1 <- mc_clopper_pearson_ci(1L, 1L, alpha_tail = 1 / 400)
  expect_equal(cp_1_0, c(lower = 0, upper = 399 / 400), tolerance = 1e-15)
  expect_equal(cp_1_1, c(lower = 1 / 400, upper = 1), tolerance = 1e-15)

  # Exact binomial tail equations and correspondence with stats::binom.test(conf.level = 199/200)
  for (d in c(2L, 5L, 10L, 20L, 80L, 2000L, 8000L)) {
    cp_0 <- mc_clopper_pearson_ci(0L, d)
    expect_identical(cp_0[["lower"]], 0)
    expect_equal((1 - cp_0[["upper"]])^d, 1 / 400, tolerance = 1e-12)

    cp_d <- mc_clopper_pearson_ci(d, d)
    expect_identical(cp_d[["upper"]], 1)
    expect_equal(cp_d[["lower"]]^d, 1 / 400, tolerance = 1e-12)
  }

  for (n in c(2L, 5L, 10L)) {
    for (k in 0L:(n - 1L)) {
      expect_equal(
        stats::pbinom(k, n, 3 / 50),
        ref_binom_cdf_recurrence(n, k, 3 / 50),
        tolerance = 1e-15
      )
      cp_nk <- mc_clopper_pearson_ci(k, n)
      bt_nk <- as.numeric(stats::binom.test(k, n, conf.level = 199 / 200)$conf.int)
      expect_equal(unname(cp_nk), bt_nk, tolerance = 1e-14)
      if (k > 0L) {
        expect_equal(stats::pbinom(k - 1L, n, cp_nk[["lower"]], lower.tail = FALSE),
                     1 / 400, tolerance = 1e-12)
      }
      expect_equal(stats::pbinom(k, n, cp_nk[["upper"]]),
                   1 / 400, tolerance = 1e-12)
    }
  }

  # Pinned stats reference check (99% CP interval at 1900/2000, 3800/4000, 7600/8000)
  expect_equal(
    unname(mc_clopper_pearson_ci(1900L, 2000L, alpha_tail = 0.005)),
    c(0.9360299778, 0.96147225),
    tolerance = 1e-2
  )
})

test_that("D15 null size conditional vs joint false-PASS witnesses (2000/20/20 and 8000/80/80), zero computability, and null premise absence", {
  # N = 2000, m_test = 20, K = 20: joint rejection is 20/2000 = 1% (joint CP upper < 0.06),
  # but conditional rejection is 20/20 = 100% (conditional CP lower > 0.06) -> MUST FAIL!
  spec_2000 <- make_mc_spec(n_seeds = 8000L, nominal_size = 0.05,
                            null_premise = TRUE, required_metrics = "null_size")
  rec_2000 <- make_mc_records(spec_2000, look = 2000L)
  rec_2000$test_computable <- c(rep(TRUE, 20L), rep(FALSE, 1980L))
  rec_2000$rejected <- c(rep(TRUE, 20L), rep(FALSE, 1980L))
  rec_2000$failure_type <- c(rep("none", 20L), rep("test_failed", 1980L))

  res_2000 <- evaluate_private_mc_cell(rec_2000, spec_2000, look = 2000L)
  ns_2000 <- res_2000$metrics$null_size
  expect_identical(ns_2000$m_computable, 20L)
  expect_identical(ns_2000$n_failed, 1980L)
  expect_identical(ns_2000$k_rejected_computable, 20L)
  expect_equal(ns_2000$joint_rejection_rate, 0.01, tolerance = 1e-15)
  expect_true(ns_2000$joint_cp_995_interval[["upper"]] < 0.06)
  expect_false(ns_2000$joint_can_confer_pass)
  expect_equal(ns_2000$estimate, 1.0, tolerance = 1e-15)
  expect_true(ns_2000$interval[["lower"]] > 0.06)
  expect_identical(ns_2000$status, "FAIL")
  expect_identical(res_2000$cell_status, "FAIL")
  expect_identical(ns_2000$denominator, "20/2000")

  # N = 8000, m_test = 80, K = 80 witness:
  # To evaluate at 8000 via evaluate_private_mc_cell, start from an INCONCLUSIVE 2000 prefix
  # (e.g., m_test = 20, K = 1 where 99.5% CP interval [0.000125, 0.317] straddles 0.06)
  # and then test both the direct 8000/80/80 witness and the 2000->8000 extension!
  rec_8000_direct <- make_mc_records(spec_2000, look = 8000L)
  rec_8000_direct$test_computable <- c(rep(TRUE, 80L), rep(FALSE, 7920L))
  rec_8000_direct$rejected <- c(rep(TRUE, 80L), rep(FALSE, 7920L))
  rec_8000_direct$failure_type <- c(rep("none", 80L), rep("test_failed", 7920L))
  ns_8000_fn <- get("mc_evaluate_null_size_metric", envir = cdid_ns, inherits = FALSE)
  ns_8000 <- ns_8000_fn(rec_8000_direct, spec_2000)
  expect_identical(ns_8000$m_computable, 80L)
  expect_identical(ns_8000$n_failed, 7920L)
  expect_equal(ns_8000$joint_rejection_rate, 0.01, tolerance = 1e-15)
  expect_true(ns_8000$joint_cp_995_interval[["upper"]] < 0.06)
  expect_true(ns_8000$interval[["lower"]] > 0.06)
  expect_identical(ns_8000$status, "FAIL")
  expect_identical(ns_8000$denominator, "80/8000")

  # Zero test computability (m_test = 0) at N = 2000:
  # Joint rejection is 0/2000 = 0 (joint CP upper < 0.06), but conditional size is UNDEFINED -> INCONCLUSIVE!
  rec_zero_test <- rec_2000
  rec_zero_test$test_computable <- rep(FALSE, 2000L)
  rec_zero_test$rejected <- rep(FALSE, 2000L)
  rec_zero_test$failure_type <- rep("test_failed", 2000L)
  res_zero_test <- evaluate_private_mc_cell(rec_zero_test, spec_2000, look = 2000L)
  expect_identical(res_zero_test$metrics$null_size$m_computable, 0L)
  expect_identical(res_zero_test$metrics$null_size$joint_rejection_count, 0L)
  expect_true(res_zero_test$metrics$null_size$joint_cp_995_interval[["upper"]] < 0.06)
  expect_identical(res_zero_test$metrics$null_size$status, "INCONCLUSIVE")
  expect_false(res_zero_test$metrics$null_size$supported)
  expect_true("zero_test_computability" %in% res_zero_test$metrics$null_size$unsupported_reasons)
  expect_identical(res_zero_test$metrics$null_size$denominator, "0/2000")

  # Null premise absence (null_premise = FALSE or NULL with null_holds = FALSE):
  # Even when truth == 0 and observed rejections are 0/2000, size is NOT inferred from data -> INCONCLUSIVE!
  spec_no_null <- make_mc_spec(
    n_seeds = 2000L,
    truth = 0.0,
    nominal_size = 0.05,
    null_premise = FALSE,
    premises = ref_default_premises(list(null_holds = FALSE)),
    required_metrics = "null_size"
  )
  rec_no_null <- make_mc_records(spec_no_null, look = 2000L, n_rejected = 0L)
  res_no_null <- evaluate_private_mc_cell(rec_no_null, spec_no_null, look = 2000L)
  expect_identical(res_no_null$metrics$null_size$status, "INCONCLUSIVE")
  expect_false(res_no_null$metrics$null_size$supported)
  expect_true("null_premise_not_declared_true" %in% res_no_null$metrics$null_size$unsupported_reasons)
})

test_that("coverage failure retention, separate 5 events, and false event values work", {
  # 49/50 conditional coverage vs 49/100 all-attempt coverage witness:
  cp_cond_50 <- mc_clopper_pearson_ci(49L, 50L)
  cp_all_100 <- mc_clopper_pearson_ci(49L, 100L)
  expect_true(cp_all_100[["upper"]] < 0.94)
  expect_true(cp_cond_50[["upper"]] > 0.94)

  # Full N = 2000 cell with 5 distinct event counts and false-event values that must NOT leak:
  spec <- make_mc_spec(n_seeds = 2000L)
  rec <- make_mc_records(spec, look = 2000L, n_covered = 2000L,
                         n_band_covered = 2000L, n_rejected = 20L)

  # Make computable counts distinct:
  # point: 1950, se: 1920 (with joint 1900), interval: 1000, band: 1960, test: 1850
  rec$point_computable[1951:2000] <- FALSE
  rec$se_computable[c(1901:1950, 1971:2000)] <- FALSE
  rec$interval_computable[1001:2000] <- FALSE
  rec$band_computable[1961:2000] <- FALSE
  rec$test_computable[1851:2000] <- FALSE

  # Rows 1901:1950 are point_computable = TRUE but se_computable = FALSE:
  # give them large finite point estimates so point-only SD != joint point/SE SD
  rec$estimate[1901:1950] <- 25.0

  any_fail <- !(rec$point_computable & rec$se_computable &
                  rec$interval_computable & rec$band_computable &
                  rec$test_computable)
  rec$failure_type[any_fail] <- "stage_partial_failure"

  # Put extreme / contradictory values on FALSE-event rows to verify they are ignored and never zero-filled:
  rec$estimate[!rec$point_computable] <- 1e9
  rec$se[!rec$se_computable] <- -999
  rec$covered[!rec$interval_computable] <- TRUE
  rec$band_covered[!rec$band_computable] <- TRUE
  rec$rejected[!rec$test_computable] <- TRUE

  res <- evaluate_private_mc_cell(rec, spec, look = 2000L)
  expect_identical(res$metrics$bias$m_computable, 1950L)
  expect_identical(res$metrics$bias$denominator, "1950/2000")
  expect_identical(res$metrics$se_ratio$m_computable, 1900L)
  expect_identical(res$metrics$se_ratio$denominator, "1900/2000")
  expect_equal(
    res$metrics$se_ratio$emp_sd,
    stats::sd(rec$estimate[rec$point_computable & rec$se_computable] - spec$truth),
    tolerance = 1e-15
  )
  expect_true(abs(res$metrics$se_ratio$emp_sd - res$metrics$bias$emp_sd) > 2.0)
  expect_identical(res$metrics$coverage$m_computable, 1000L)
  expect_identical(res$metrics$coverage$k_valid_and_covered, 1000L)
  expect_equal(res$metrics$coverage$conditional_coverage_rate, 1.0, tolerance = 1e-15)
  expect_equal(res$metrics$coverage$estimate, 0.5, tolerance = 1e-15)
  expect_identical(res$metrics$coverage$status, "FAIL")
  expect_identical(res$metrics$band_coverage$m_computable, 1960L)
  expect_identical(res$metrics$band_coverage$k_valid_and_covered, 1960L)
  expect_equal(res$metrics$band_coverage$estimate, 1960 / 2000, tolerance = 1e-15)
  expect_identical(res$metrics$band_coverage$status, "PASS")
  expect_identical(res$metrics$null_size$m_computable, 1850L)
  expect_identical(res$metrics$null_size$k_rejected_computable, 20L)
  expect_identical(res$metrics$null_size$denominator, "1850/2000")
  expect_identical(res$metrics$null_size$status, "PASS")
  expect_identical(res$cell_status, "FAIL")
})

test_that("t3 vs t5 fourth-moment premise and other external premises are never inferred from finite data", {
  # Keep the original finite 2000-row data. Bias lies inside [-0.10, 0.10],
  # while critical 4.0 widens the SE-ratio interval across [0.95, 1.05].
  # A fourth-moment premise is still external (FALSE as in t3).
  spec_t3 <- make_mc_spec(
    n_seeds = 2000L,
    required_metrics = c("bias", "se_ratio"),
    premises = ref_default_premises(list(
      conditional_fourth_moment_T_point = FALSE,
      conditional_fourth_moment_T_joint = FALSE
    ))
  )
  rec_t3 <- make_mc_records(spec_t3, look = 2000L, bias_shift = 0.0, se_scale = 1.0)
  res_t3 <- evaluate_private_mc_cell(rec_t3, spec_t3, look = 2000L)

  expect_true(res_t3$metrics$bias$arithmetic_defined)
  expect_true(all(is.finite(res_t3$metrics$bias$descriptive_interval)))
  expect_true(res_t3$metrics$bias$descriptive_interval[["lower"]] >= -0.10 &&
                res_t3$metrics$bias$descriptive_interval[["upper"]] <= 0.10)
  expect_false(res_t3$metrics$bias$supported)
  expect_true(all(is.na(res_t3$metrics$bias$interval)))
  expect_identical(res_t3$metrics$bias$status, "INCONCLUSIVE")
  expect_true("conditional_fourth_moment_T_point" %in% res_t3$metrics$bias$unsupported_reasons)

  expect_true(res_t3$metrics$se_ratio$arithmetic_defined)
  expect_false(res_t3$metrics$se_ratio$supported)
  expect_identical(res_t3$metrics$se_ratio$status, "INCONCLUSIVE")
  expect_true("conditional_fourth_moment_T_joint" %in% res_t3$metrics$se_ratio$unsupported_reasons)
  expect_identical(res_t3$cell_status, "INCONCLUSIVE")

  # Same exact data under t5: bias passes, while the wider SE-ratio interval
  # is arithmetically supported but INCONCLUSIVE; the premise is not inferred.
  spec_t5 <- make_mc_spec(
    n_seeds = 2000L,
    required_metrics = c("bias", "se_ratio"),
    premises = ref_default_premises()
  )
  res_t5 <- evaluate_private_mc_cell(rec_t3, spec_t5, look = 2000L)
  expect_identical(res_t5$metrics$bias$status, "PASS")
  expect_identical(res_t5$metrics$se_ratio$status, "INCONCLUSIVE")
  expect_identical(res_t5$cell_status, "INCONCLUSIVE")
  ref_t5 <- cp02_mc_expect_ratio_reference(res_t5$metrics$se_ratio, rec_t3, spec_t5)
  expect_true(res_t5$metrics$se_ratio$supported)
  expect_true(ref_t5$ci_r_legacy_qnorm[["lower"]] >= 0.95 &&
    ref_t5$ci_r_legacy_qnorm[["upper"]] <= 1.05)
  expect_equal(res_t3$metrics$se_ratio$descriptive_interval, ref_t5$ci_r, tolerance = 1e-15)

  # Missing calibration support or missing second moment of S also yields INCONCLUSIVE:
  spec_no_cal <- make_mc_spec(
    n_seeds = 2000L,
    required_metrics = c("bias", "se_ratio"),
    premises = ref_default_premises(list(
      calibration_supported_bias = FALSE,
      conditional_second_moment_S_joint = FALSE
    ))
  )
  res_no_cal <- evaluate_private_mc_cell(rec_t3, spec_no_cal, look = 2000L)
  expect_identical(res_no_cal$metrics$bias$status, "INCONCLUSIVE")
  expect_true("calibration_supported_bias" %in% res_no_cal$metrics$bias$unsupported_reasons)
  expect_identical(res_no_cal$metrics$se_ratio$status, "INCONCLUSIVE")
  expect_true("conditional_second_moment_S_joint" %in% res_no_cal$metrics$se_ratio$unsupported_reasons)
})

test_that("CI threshold equalities obey strict lower-bound coverage rule and inclusive bias/ratio/size rules", {
  # Coverage and band_coverage: floor = 0.95 - 0.01 = 0.94
  # Strict lower bound > 0.94 required for PASS; equality at lower or upper is INCONCLUSIVE!
  for (m_cov in c("coverage", "band_coverage")) {
    expect_identical(mc_classify_interval_decision(m_cov, 0.94, 0.97, nominal = 0.95), "INCONCLUSIVE")
    expect_identical(mc_classify_interval_decision(m_cov, 0.94 + 1e-12, 0.97, nominal = 0.95), "PASS")
    expect_identical(mc_classify_interval_decision(m_cov, 0.91, 0.94, nominal = 0.95), "INCONCLUSIVE")
    expect_identical(mc_classify_interval_decision(m_cov, 0.91, 0.94 - 1e-12, nominal = 0.95), "FAIL")
    expect_identical(mc_classify_interval_decision(m_cov, 0.93, 0.95, nominal = 0.95), "INCONCLUSIVE")
  }

  # Null size: ceiling = 0.05 + 0.01 = 0.06
  expect_identical(mc_classify_interval_decision("null_size", 0.02, 0.06, nominal = 0.05), "PASS")
  expect_identical(mc_classify_interval_decision("null_size", 0.02, 0.06 + 1e-12, nominal = 0.05), "INCONCLUSIVE")
  expect_identical(mc_classify_interval_decision("null_size", 0.06, 0.08, nominal = 0.05), "INCONCLUSIVE")
  expect_identical(mc_classify_interval_decision("null_size", 0.06 + 1e-12, 0.08, nominal = 0.05), "FAIL")

  # Signed normalized bias: inclusive [-0.10, 0.10]
  expect_identical(mc_classify_interval_decision("bias", -0.10, 0.10), "PASS")
  expect_identical(mc_classify_interval_decision("bias", -0.10, 0.00), "PASS")
  expect_identical(mc_classify_interval_decision("bias", 0.00, 0.10), "PASS")
  expect_identical(mc_classify_interval_decision("bias", 0.10, 0.15), "INCONCLUSIVE")
  expect_identical(mc_classify_interval_decision("bias", -0.15, -0.10), "INCONCLUSIVE")
  expect_identical(mc_classify_interval_decision("bias", 0.10 + 1e-12, 0.15), "FAIL")
  expect_identical(mc_classify_interval_decision("bias", -0.15, -0.10 - 1e-12), "FAIL")

  # SE ratio: inclusive [0.95, 1.05]
  expect_identical(mc_classify_interval_decision("se_ratio", 0.95, 1.05), "PASS")
  expect_identical(mc_classify_interval_decision("se_ratio", 0.95, 1.00), "PASS")
  expect_identical(mc_classify_interval_decision("se_ratio", 1.00, 1.05), "PASS")
  expect_identical(mc_classify_interval_decision("se_ratio", 1.05, 1.10), "INCONCLUSIVE")
  expect_identical(mc_classify_interval_decision("se_ratio", 0.90, 0.95), "INCONCLUSIVE")
  expect_identical(mc_classify_interval_decision("se_ratio", 1.05 + 1e-12, 1.10), "FAIL")
  expect_identical(mc_classify_interval_decision("se_ratio", 0.90, 0.95 - 1e-12), "FAIL")

  # Full evaluate_private_mc_cell exact IEEE-754 boundary equality witnesses at look = 2000:
  l_1800 <- stats::qbeta(1 / 400, 1800L, 201L)
  u_1800 <- stats::qbeta(399 / 400, 1801L, 200L)
  expect_identical((l_1800 + 0.01) - 0.01, l_1800)
  expect_identical((u_1800 + 0.01) - 0.01, u_1800)

  spec_cov_eq_l <- make_mc_spec(n_seeds = 2000L, nominal_coverage = l_1800 + 0.01,
                                required_metrics = c("coverage", "band_coverage"))
  rec_cov_eq_l <- make_mc_records(spec_cov_eq_l, look = 2000L,
                                  n_covered = 1800L, n_band_covered = 1800L)
  res_cov_eq_l <- evaluate_private_mc_cell(rec_cov_eq_l, spec_cov_eq_l, look = 2000L)
  expect_identical(res_cov_eq_l$metrics$coverage$interval[["lower"]],
                   spec_cov_eq_l$nominal_coverage - 0.01)
  expect_identical(res_cov_eq_l$metrics$coverage$status, "INCONCLUSIVE")
  expect_identical(res_cov_eq_l$metrics$band_coverage$status, "INCONCLUSIVE")
  expect_identical(res_cov_eq_l$cell_status, "INCONCLUSIVE")

  spec_cov_eq_u <- make_mc_spec(n_seeds = 2000L, nominal_coverage = u_1800 + 0.01,
                                required_metrics = c("coverage", "band_coverage"))
  res_cov_eq_u <- evaluate_private_mc_cell(rec_cov_eq_l, spec_cov_eq_u, look = 2000L)
  expect_identical(res_cov_eq_u$metrics$coverage$interval[["upper"]],
                   spec_cov_eq_u$nominal_coverage - 0.01)
  expect_identical(res_cov_eq_u$metrics$coverage$status, "INCONCLUSIVE")
  expect_identical(res_cov_eq_u$metrics$band_coverage$status, "INCONCLUSIVE")

  u_50 <- stats::qbeta(399 / 400, 51L, 1950L)
  l_50 <- stats::qbeta(1 / 400, 50L, 1951L)
  expect_identical((u_50 - 0.01) + 0.01, u_50)
  expect_identical((l_50 - 0.01) + 0.01, l_50)
  spec_sz_eq_u <- make_mc_spec(n_seeds = 2000L, nominal_size = u_50 - 0.01,
                               required_metrics = "null_size")
  rec_sz_50 <- make_mc_records(spec_sz_eq_u, look = 2000L, n_rejected = 50L)
  res_sz_eq_u <- evaluate_private_mc_cell(rec_sz_50, spec_sz_eq_u, look = 2000L)
  expect_identical(res_sz_eq_u$metrics$null_size$interval[["upper"]],
                   spec_sz_eq_u$nominal_size + 0.01)
  expect_identical(res_sz_eq_u$metrics$null_size$status, "PASS")

  spec_sz_eq_l <- make_mc_spec(n_seeds = 2000L, nominal_size = l_50 - 0.01,
                               required_metrics = "null_size")
  res_sz_eq_l <- evaluate_private_mc_cell(rec_sz_50, spec_sz_eq_l, look = 2000L)
  expect_identical(res_sz_eq_l$metrics$null_size$interval[["lower"]],
                   spec_sz_eq_l$nominal_size + 0.01)
  expect_identical(res_sz_eq_l$metrics$null_size$status, "INCONCLUSIVE")

  # Full-cell bias and se_ratio PASS, INCONCLUSIVE, and FAIL regimes:
  spec_br <- make_mc_spec(n_seeds = 2000L, required_metrics = c("bias", "se_ratio"))
  rec_br_fail <- make_mc_records(spec_br, look = 2000L, bias_shift = 0.25, se_scale = 1.25)
  res_br_fail <- evaluate_private_mc_cell(rec_br_fail, spec_br, look = 2000L)
  expect_identical(res_br_fail$metrics$bias$status, "FAIL")
  expect_identical(res_br_fail$metrics$se_ratio$status, "FAIL")
  expect_identical(res_br_fail$cell_status, "FAIL")
})

test_that("all-fail N=2000 and record/specification/attempt/endpoint corruptions behave as specified", {
  spec <- make_mc_spec(n_seeds = 2000L)
  rec_allfail <- make_mc_records(spec, look = 2000L)
  rec_allfail$point_computable <- rep(FALSE, 2000L)
  rec_allfail$se_computable <- rep(FALSE, 2000L)
  rec_allfail$interval_computable <- rep(FALSE, 2000L)
  rec_allfail$band_computable <- rep(FALSE, 2000L)
  rec_allfail$test_computable <- rep(FALSE, 2000L)
  rec_allfail$estimate <- rep(NA_real_, 2000L)
  rec_allfail$se <- rep(NA_real_, 2000L)
  rec_allfail$covered <- rep(NA, 2000L)
  rec_allfail$band_covered <- rep(NA, 2000L)
  rec_allfail$rejected <- rep(NA, 2000L)
  rec_allfail$failure_type <- rep("generation_failed", 2000L)
  rec_allfail$failure_message <- rep("all attempts failed", 2000L)

  res_af <- evaluate_private_mc_cell(rec_allfail, spec, look = 2000L)
  expect_identical(res_af$metrics$coverage$status, "FAIL")
  expect_identical(res_af$metrics$band_coverage$status, "FAIL")
  expect_identical(res_af$metrics$null_size$status, "INCONCLUSIVE")
  expect_identical(res_af$metrics$bias$status, "INCONCLUSIVE")
  expect_identical(res_af$metrics$se_ratio$status, "INCONCLUSIVE")
  expect_identical(res_af$cell_status, "FAIL")
  expect_false(res_af$can_extend_to_8000)

  # Intermediate or development lengths are refused
  rec_valid <- make_mc_records(spec, look = 2000L)
  for (bad_look in c(200L, 1000L, 1999L, 2001L, 4000L, 5000L, 10000L)) {
    expect_error(
      evaluate_private_mc_cell(rec_valid, spec, look = bad_look),
      "look.*must be 2000L or 8000L"
    )
  }

  # Duplicated, missing, or reordered attempts
  rec_dup_att <- rec_valid
  rec_dup_att$attempt[2] <- 1L
  expect_error(evaluate_private_mc_cell(rec_dup_att, spec, look = 2000L), "attempt")

  rec_reord_att <- rec_valid
  rec_reord_att$attempt[1:2] <- c(2L, 1L)
  expect_error(evaluate_private_mc_cell(rec_reord_att, spec, look = 2000L), "attempt")

  rec_miss_att <- rec_valid
  rec_miss_att$attempt[2000] <- 2001L
  expect_error(evaluate_private_mc_cell(rec_miss_att, spec, look = 2000L), "attempt")

  # Duplicated or reordered seeds
  rec_reord_seed <- rec_valid
  rec_reord_seed$seed[1:2] <- rec_valid$seed[2:1]
  expect_error(evaluate_private_mc_cell(rec_reord_seed, spec, look = 2000L), "seed")

  # Mismatched truth, target, cell_id
  rec_bad_truth <- rec_valid
  rec_bad_truth$truth[10] <- spec$truth + 1
  expect_error(evaluate_private_mc_cell(rec_bad_truth, spec, look = 2000L), "truth")

  rec_bad_target <- rec_valid
  rec_bad_target$target[10] <- "wrong_target"
  expect_error(evaluate_private_mc_cell(rec_bad_target, spec, look = 2000L), "target")

  rec_bad_cell <- rec_valid
  rec_bad_cell$cell_id[10] <- "wrong_cell"
  expect_error(evaluate_private_mc_cell(rec_bad_cell, spec, look = 2000L), "cell_id")

  # Computed flag TRUE with invalid value is an input refusal, never silently turned into failure
  rec_bad_est <- rec_valid
  rec_bad_est$estimate[1] <- NA_real_
  expect_error(evaluate_private_mc_cell(rec_bad_est, spec, look = 2000L), "point_computable")

  rec_inf_est <- rec_valid
  rec_inf_est$estimate[1] <- Inf
  expect_error(evaluate_private_mc_cell(rec_inf_est, spec, look = 2000L), "point_computable")

  rec_neg_se <- rec_valid
  rec_neg_se$se[1] <- -1e-6
  expect_error(evaluate_private_mc_cell(rec_neg_se, spec, look = 2000L), "se_computable")

  rec_na_se <- rec_valid
  rec_na_se$se[1] <- NA_real_
  expect_error(evaluate_private_mc_cell(rec_na_se, spec, look = 2000L), "se_computable")

  rec_na_cov <- rec_valid
  rec_na_cov$covered[1] <- NA
  expect_error(evaluate_private_mc_cell(rec_na_cov, spec, look = 2000L), "interval_computable")

  # Endpoint binding check when lower and upper are present
  rec_endpt <- rec_valid
  rec_endpt$lower <- rec_endpt$estimate - 1.96 * rec_endpt$se
  rec_endpt$upper <- rec_endpt$estimate + 1.96 * rec_endpt$se
  rec_endpt$covered <- rec_endpt$lower <= spec$truth & spec$truth <= rec_endpt$upper
  expect_silent(evaluate_private_mc_cell(rec_endpt, spec, look = 2000L))
  rec_endpt$covered[1] <- !rec_endpt$covered[1]
  expect_error(evaluate_private_mc_cell(rec_endpt, spec, look = 2000L), "lower <= truth")

  # Merged event columns refused
  spec_merged_ev <- spec
  spec_merged_ev$events$se <- "point_computable"
  expect_error(evaluate_private_mc_cell(rec_valid, spec_merged_ev, look = 2000L), "separate fields")
})

test_that("D14 two-look extension (2000 -> 8000) preserves conclusive decisions, reevaluates inconclusive only, and enforces prefix/status guards", {
  spec <- make_mc_spec(
    n_seeds = 8000L,
    required_metrics = c("coverage", "bias", "se_ratio")
  )
  # At look 2000:
  # - coverage has 1930/2000 = 0.965 -> L = 0.9518 > 0.94 -> PASS
  # - se_ratio has se_scale = 1.0 but critical 4.0 crosses the D7 margins -> INCONCLUSIVE
  # - bias has bias_shift = 0.085 -> CI straddles 0.10 -> INCONCLUSIVE
  # Therefore cell_status at 2000 is INCONCLUSIVE; both continuous metrics extend.
  rec_8000 <- make_mc_records(spec, look = 8000L, bias_shift = 0.02, se_scale = 1.0,
                              n_covered = 7720L)
  rec_2000 <- make_mc_records(spec, look = 2000L, bias_shift = 0.085, se_scale = 1.0,
                              n_covered = 1930L)
  rec_8000[1:2000, ] <- rec_2000

  # Make rows 2001:8000 have terrible coverage (0 covered!) to prove that look-2000 conclusive
  # coverage PASS is PRESERVED and NOT overturned at look 8000!
  rec_8000$covered[2001:8000] <- FALSE

  res_2000 <- evaluate_private_mc_cell(rec_2000, spec, look = 2000L)
  expect_identical(res_2000$metrics$coverage$status, "PASS")
  expect_identical(res_2000$metrics$se_ratio$status, "INCONCLUSIVE")
  ref_2000 <- cp02_mc_expect_ratio_reference(res_2000$metrics$se_ratio, rec_2000, spec)
  expect_true(ref_2000$ci_r_legacy_qnorm[["lower"]] >= 0.95 &&
    ref_2000$ci_r_legacy_qnorm[["upper"]] <= 1.05)
  expect_identical(res_2000$metrics$bias$status, "INCONCLUSIVE")
  expect_identical(res_2000$cell_status, "INCONCLUSIVE")
  expect_true(res_2000$can_extend_to_8000)

  res_8000 <- evaluate_private_mc_cell(rec_8000, spec, look = 8000L,
                                       prior_assessment = res_2000)
  expect_identical(res_8000$look, 8000L)
  expect_identical(res_8000$metrics$coverage$status, "PASS")
  expect_identical(res_8000$metrics$coverage$decision_look, 2000L)
  expect_true(res_8000$metrics$coverage$preserved_from_look_2000)
  expect_false(res_8000$metrics$coverage$reevaluated_at_8000)
  expect_identical(res_8000$metrics$coverage$look_8000_descriptive_only$status, "FAIL")

  expect_identical(res_8000$metrics$se_ratio$status, "PASS")
  cp02_mc_expect_ratio_reference(res_8000$metrics$se_ratio, rec_8000, spec)
  expect_identical(res_8000$metrics$se_ratio$decision_look, 8000L)
  expect_false(res_8000$metrics$se_ratio$preserved_from_look_2000)
  expect_true(res_8000$metrics$se_ratio$reevaluated_at_8000)

  expect_identical(res_8000$metrics$bias$status, "PASS")
  expect_identical(res_8000$metrics$bias$decision_look, 8000L)
  expect_false(res_8000$metrics$bias$preserved_from_look_2000)
  expect_true(res_8000$metrics$bias$reevaluated_at_8000)

  expect_identical(res_8000$cell_status, "PASS")
  expect_false(res_8000$can_extend_to_8000)
  expect_false(res_8000$unresolved_at_cap)
  expect_identical(res_8000$history$preserved_metrics, "coverage")
  expect_identical(res_8000$history$reevaluated_metrics, c("bias", "se_ratio"))

  # Unresolved at cap (8000) remains INCONCLUSIVE
  rec_8000_inc <- make_mc_records(spec, look = 8000L, bias_shift = 0.095, se_scale = 1.0,
                                  n_covered = 7720L)
  rec_8000_inc[1:2000, ] <- rec_2000
  res_8000_inc <- evaluate_private_mc_cell(rec_8000_inc, spec, look = 8000L,
                                           prior_assessment = res_2000)
  expect_identical(res_8000_inc$metrics$bias$status, "INCONCLUSIVE")
  expect_identical(res_8000_inc$cell_status, "INCONCLUSIVE")
  expect_true(res_8000_inc$unresolved_at_cap)
  expect_false(res_8000_inc$can_extend_to_8000)

  # Guard 1: look = 8000 without prior_assessment
  expect_error(
    evaluate_private_mc_cell(rec_8000, spec, look = 8000L, prior_assessment = NULL),
    "prior_assessment"
  )

  # Guard 2: look = 2000 with non-NULL prior_assessment
  expect_error(
    evaluate_private_mc_cell(rec_2000, spec, look = 2000L, prior_assessment = res_2000),
    "prior_assessment.*must be NULL"
  )

  # The original zero-bias fixture also has an INCONCLUSIVE SE-ratio interval
  # under critical 4.0. Preserve it and verify its legitimate extension.
  # The conclusive-cell refusal is checked on a new deterministic case below.
  rec_2000_pass <- make_mc_records(spec, look = 2000L, bias_shift = 0.0, se_scale = 1.0,
                                   n_covered = 1930L)
  res_2000_pass <- evaluate_private_mc_cell(rec_2000_pass, spec, look = 2000L)
  expect_identical(res_2000_pass$cell_status, "INCONCLUSIVE")
  expect_true(res_2000_pass$can_extend_to_8000)
  cp02_mc_expect_ratio_reference(res_2000_pass$metrics$se_ratio, rec_2000_pass, spec)
  rec_8000_from_pass <- rec_8000
  rec_8000_from_pass[1:2000, ] <- rec_2000_pass
  res_8000_from_pass <- evaluate_private_mc_cell(rec_8000_from_pass, spec, look = 8000L,
    prior_assessment = res_2000_pass)
  cp02_mc_expect_ratio_reference(res_8000_from_pass$metrics$se_ratio, rec_8000_from_pass, spec)
  expect_identical(res_8000_from_pass$cell_status, "PASS")
  expect_identical(res_8000_from_pass$metrics$se_ratio$decision_look, 8000L)
  expect_true(res_8000_from_pass$metrics$se_ratio$reevaluated_at_8000)

  # Guard 4: FAIL + INCONCLUSIVE cell at 2000 has cell_status == "FAIL" and CANNOT extend!
  rec_2000_fail_inc <- rec_2000
  rec_2000_fail_inc$covered[1:500] <- FALSE # coverage = 1430/2000 = 0.715 -> FAIL, while bias is INCONCLUSIVE
  res_2000_fail_inc <- evaluate_private_mc_cell(rec_2000_fail_inc, spec, look = 2000L)
  expect_identical(res_2000_fail_inc$metrics$coverage$status, "FAIL")
  expect_identical(res_2000_fail_inc$metrics$bias$status, "INCONCLUSIVE")
  expect_identical(res_2000_fail_inc$cell_status, "FAIL")
  expect_false(res_2000_fail_inc$can_extend_to_8000)
  rec_8000_from_fail_inc <- rec_8000
  rec_8000_from_fail_inc[1:2000, ] <- rec_2000_fail_inc
  expect_error(
    evaluate_private_mc_cell(rec_8000_from_fail_inc, spec, look = 8000L,
                             prior_assessment = res_2000_fail_inc),
    "cell-INCONCLUSIVE"
  )

  # Guard 5: Forged decorative status on prior_assessment is caught by deterministic reproduction
  forged_prior <- res_2000_fail_inc
  forged_prior$cell_status <- "INCONCLUSIVE"
  expect_error(
    evaluate_private_mc_cell(rec_8000_from_fail_inc, spec, look = 8000L,
                             prior_assessment = forged_prior),
    "deterministic reproduction"
  )

  # Guard 6: Prefix modification in first 2000 rows of records at 8000
  rec_8000_tampered <- rec_8000
  rec_8000_tampered$estimate[2000] <- rec_8000_tampered$estimate[2000] + 1e-6
  expect_error(
    evaluate_private_mc_cell(rec_8000_tampered, spec, look = 8000L,
                             prior_assessment = res_2000),
    "records_prefix"
  )

  # Guard 7: Changed original seeds in specification at 8000
  spec_bad_seeds <- spec
  spec_bad_seeds$seeds[1:2] <- spec_bad_seeds$seeds[2:1]
  rec_8000_bad_seeds <- rec_8000
  rec_8000_bad_seeds$seed[1:2] <- spec_bad_seeds$seeds[1:2]
  expect_error(
    evaluate_private_mc_cell(rec_8000_bad_seeds, spec_bad_seeds, look = 8000L,
                             prior_assessment = res_2000),
    "specification"
  )

  # Guard 8: Changed specification truth/target/metrics at 8000
  spec_changed <- spec
  spec_changed$required_metrics <- c("coverage", "bias")
  expect_error(
    evaluate_private_mc_cell(rec_8000, spec_changed, look = 8000L,
                             prior_assessment = res_2000),
    "specification"
  )
})

test_that("D14 light-tailed fixed records preserve conclusive ratio decisions and refuse all-conclusive extension", {
  spec <- make_mc_spec(n_seeds = 8000L, required_metrics = c("coverage", "bias", "se_ratio"))
  # No random draws: exactly representable two-point errors and reported-SE
  # atoms. At 2000, SD=sqrt(2000/1999), mean SE=1, and the centered
  # ratio IF is +/-1/64, so its variance is (1/64)^2*2000/1999.
  # Binary atoms avoid treating an ideal decimal shift/normalization as exact.
  light_records <- function(look, shift) {
    rec <- make_mc_records(spec, look = look, bias_shift = 0, se_scale = 1)
    rec$estimate <- spec$truth + shift + rep(c(-1, 1), length.out = look)
    rec$se <- rep(c(63/64, 63/64, 65/64, 65/64), length.out = look)
    rec
  }
  first <- light_records(2000L, 3/32)
  full <- light_records(8000L, 1/64)
  full$se[2001:8000] <- 2 * full$se[2001:8000]
  full$covered[2001:8000] <- FALSE
  full[1:2000, ] <- first
  ref_first <- ref_sample_moments(first$estimate - spec$truth, first$se)
  ref_full <- ref_sample_moments(full$estimate - spec$truth, full$se)
  expect_equal(ref_first$sd, sqrt(2000 / 1999), tolerance = 1e-15)
  expect_equal(ref_first$vr, (1/64)^2 * 2000 / 1999, tolerance = 1e-15)
  expect_equal(ref_first$mcse_r, 1 / (64 * sqrt(1999)), tolerance = 1e-15)
  expect_equal(ref_first$ci_r, sqrt(1999 / 2000) * c(lower = exp(-4 / (64 * sqrt(1999))),
    upper = exp(4 / (64 * sqrt(1999)))), tolerance = 1e-15)
  expect_true(ref_first$ci_r[["lower"]] >= .95 && ref_first$ci_r[["upper"]] <= 1.05)
  expect_true(ref_first$ci_b[["lower"]] < .10 && ref_first$ci_b[["upper"]] > .10)
  expect_true(ref_full$ci_b[["lower"]] >= -.10 && ref_full$ci_b[["upper"]] <= .10)
  expect_true(ref_full$ci_r[["lower"]] > 1.05)

  res_first <- evaluate_private_mc_cell(first, spec, look = 2000L)
  cp02_mc_expect_ratio_reference(res_first$metrics$se_ratio, first, spec)
  expect_equal(res_first$metrics$bias$interval, ref_first$ci_b, tolerance = 1e-15)
  expect_identical(res_first$metrics$coverage$status, "PASS")
  expect_identical(res_first$metrics$se_ratio$status, "PASS")
  expect_identical(res_first$metrics$bias$status, "INCONCLUSIVE")
  expect_identical(res_first$cell_status, "INCONCLUSIVE")
  expect_true(res_first$can_extend_to_8000)
  res_full <- evaluate_private_mc_cell(full, spec, look = 8000L, prior_assessment = res_first)
  # Both conclusive first-look metrics stay frozen, even when their full-8000
  # descriptive intervals fail. Only the previously inconclusive bias extends.
  for (metric in c("coverage", "se_ratio")) {
    expect_identical(res_full$metrics[[metric]]$status, "PASS")
    expect_identical(res_full$metrics[[metric]]$interval, res_first$metrics[[metric]]$interval)
    expect_identical(res_full$metrics[[metric]]$decision_look, 2000L)
    expect_true(res_full$metrics[[metric]]$preserved_from_look_2000)
    expect_false(res_full$metrics[[metric]]$reevaluated_at_8000)
    expect_identical(res_full$metrics[[metric]]$look_8000_descriptive_only$status, "FAIL")
  }
  cp02_mc_expect_ratio_reference(res_full$metrics$se_ratio$look_8000_descriptive_only, full, spec)
  expect_equal(res_full$metrics$bias$interval, ref_full$ci_b, tolerance = 1e-15)
  expect_identical(res_full$metrics$bias$status, "PASS")
  expect_identical(res_full$metrics$bias$decision_look, 8000L)
  expect_false(res_full$metrics$bias$preserved_from_look_2000)
  expect_true(res_full$metrics$bias$reevaluated_at_8000)
  expect_identical(res_full$cell_status, "PASS")
  expect_false(res_full$can_extend_to_8000)
  expect_false(res_full$unresolved_at_cap)
  expect_identical(res_full$history$preserved_metrics, c("coverage", "se_ratio"))
  expect_identical(res_full$history$reevaluated_metrics, "bias")

  conclusive <- light_records(2000L, 0)
  ref_conclusive <- ref_sample_moments(conclusive$estimate - spec$truth, conclusive$se)
  expect_true(ref_conclusive$ci_b[["lower"]] >= -.10 && ref_conclusive$ci_b[["upper"]] <= .10)
  res_conclusive <- evaluate_private_mc_cell(conclusive, spec, look = 2000L)
  cp02_mc_expect_ratio_reference(res_conclusive$metrics$se_ratio, conclusive, spec)
  expect_identical(res_conclusive$cell_status, "PASS")
  expect_false(res_conclusive$can_extend_to_8000)
  full[1:2000, ] <- conclusive
  expect_error(evaluate_private_mc_cell(full, spec, look = 8000L,
    prior_assessment = res_conclusive), "cell-INCONCLUSIVE")
})

test_that("four-law 50-outer calibration plumbing and cell-rule continuation overlay match specifications", {
  plumb <- prepare_mc_calibration_plumbing()
  expect_identical(plumb$status, "PROPOSED_PLUMBING_ONLY_UNLAUNCHED")
  expect_false(plumb$execution_now)
  expect_identical(plumb$outer_draws_executed, 0L)
  expect_identical(plumb$inner_attempts_executed, 0L)
  expect_identical(plumb$fits_executed, 0L)
  expect_identical(plumb$historical_plan_sha256,
                   "d7f055c8d201146feb2c55900f16d7f0605e063f09869080a1687d4a8a44b683")
  expect_false(plumb$historical_plan_modified)
  expect_identical(plumb$outer_draws_per_cell, 50L)
  expect_identical(plumb$cell_count, 4L)
  expect_identical(plumb$outer_draws_total, 200L)
  expect_identical(plumb$inner_attempt_assessments, c(2000L, 8000L))
  expect_identical(plumb$batch_rows_maximum, 1000L)
  expect_identical(plumb$maximum_inner_attempts_upper_bound, 1600000L)
  expect_identical(plumb$maximum_primitive_scalar_draws_Z_plus_U, 3200000L)
  expect_identical(plumb$maximum_outer_look_rows, 400L)
  expect_equal(plumb$expected_marginal_misses_per_law_at_995, 0.25, tolerance = 1e-15)
  expect_false(plumb$can_certify_995_calibration)
  expect_false(plumb$can_certify_fitted_estimator)
  expect_false(plumb$bootstrap_precision_limits_adopted)
  expect_false(plumb$failure_ceiling_adopted)

  # Analytic influence variances for the first three laws (a=1/8, b=1/10):
  expect_equal(plumb$cells[[1]]$Var_psi_beta, 201 / 200, tolerance = 1e-15)
  expect_equal(plumb$cells[[1]]$Var_psi_logrho, 9 / 32, tolerance = 1e-15)
  expect_equal(plumb$cells[[2]]$Var_psi_beta, 51 / 50, tolerance = 1e-15)
  expect_equal(plumb$cells[[2]]$Var_psi_logrho, 9 / 8, tolerance = 1e-15)
  expect_equal(plumb$cells[[3]]$Var_psi_beta, 41 / 50, tolerance = 1e-15)
  expect_equal(plumb$cells[[3]]$Var_psi_logrho, 9 / 8, tolerance = 1e-15)
  expect_equal(plumb$cells[[4]]$computability, 1 / 2, tolerance = 1e-15)
  expect_equal(plumb$cells[[4]]$coverage_all_attempt, 19 / 40, tolerance = 1e-15)

  # Seed formula check
  expect_identical(mc_calibration_outer_seed(1L, 1L), 14332001L)
  expect_identical(mc_calibration_outer_seed(4L, 50L), 14362050L)
  expect_identical(nrow(plumb$outer_seed_grid), 200L)
  expect_false(any(plumb$outer_seed_grid$executed))

  # Dependent-computability cell 4 overlay check: low all-attempt coverage (19/40 = 47.5%)
  # causes coverage FAIL at look 2000, so cell-level continuation rule stops at 2000 even if bias is INCONCLUSIVE.
  spec_c4 <- make_mc_spec(n_seeds = 2000L, required_metrics = c("coverage", "bias"))
  rec_c4 <- make_mc_records(spec_c4, look = 2000L, bias_shift = 0.085, n_covered = 950L)
  res_c4 <- evaluate_private_mc_cell(rec_c4, spec_c4, look = 2000L)
  cont_c4 <- mc_calibration_continuation_decision(res_c4)
  expect_identical(res_c4$metrics$coverage$status, "FAIL")
  expect_identical(res_c4$metrics$bias$status, "INCONCLUSIVE")
  expect_identical(cont_c4$cell_status_2000, "FAIL")
  expect_false(cont_c4$continue_to_8000)
  expect_identical(cont_c4$stop_reason, "cell_fail_plus_inconclusive_stop_at_2000_no_extension")
})

test_that("independent audit: scalar records cannot recycle across matrix columns", {
  spec <- make_mc_spec(n_seeds = 8000L, required_metrics = "bias")
  rec <- make_mc_records(spec, look = 2000L)
  rec$point_computable[1001:2000] <- FALSE
  rec$failure_type[1001:2000] <- "point_not_computed"
  rec$estimate <- I(cbind(rep(1, 2000L), rep(-1, 2000L)))
  expect_error(evaluate_private_mc_cell(rec, spec, 2000L), "scalar-per-attempt")
  base_rec <- make_mc_records(spec, look = 2000L)
  cols <- c("attempt", "seed", "estimate", "se", "covered", "band_covered", "rejected",
    "point_computable", "se_computable", "interval_computable", "band_computable", "test_computable",
    "failure_type", "failure_message", "n_warnings", "warnings")
  for (nm in cols) {
    bad <- base_rec; bad[[nm]] <- I(cbind(bad[[nm]], bad[[nm]]))
    expect_error(evaluate_private_mc_cell(bad, spec, 2000L), "scalar-per-attempt")
  }
  bad <- base_rec; bad$extra <- bad$estimate; names(bad)[ncol(bad)] <- "estimate"
  expect_error(evaluate_private_mc_cell(bad, spec, 2000L), "unique column names")
})

test_that("independent audit: finite stored inputs with unrepresentable errors are INCONCLUSIVE", {
  spec <- make_mc_spec(n_seeds = 8000L, truth = -1e308, required_metrics = c("bias", "se_ratio"))
  rec <- make_mc_records(spec, look = 2000L)
  rec$estimate <- rep(c(-1e308, 1e308), 1000L)
  rec$se[1L] <- 0
  res <- evaluate_private_mc_cell(rec, spec, 2000L)
  expect_identical(res$cell_status, "INCONCLUSIVE")
  for (nm in c("bias", "se_ratio")) {
    expect_identical(res$metrics[[nm]]$status, "INCONCLUSIVE")
    expect_false(res$metrics[[nm]]$arithmetic_defined)
    expect_true("nonfinite_error_arithmetic" %in% res$metrics[[nm]]$unsupported_reasons)
    expect_identical(res$metrics[[nm]]$m_computable, 2000L)
    expect_identical(res$metrics[[nm]]$denominator, "2000/2000")
  }
  expect_identical(res$metrics$se_ratio$n_zero_se, 1L)
  expect_identical(unname(res$failure_and_event_summary$computable_counts["point"]), 2000L)
  imax <- .Machine$integer.max
  spec$truth <- -imax; rec$truth[] <- -imax
  rec$estimate <- rep(c(-imax, imax), 1000L)
  res_int <- evaluate_private_mc_cell(rec, spec, 2000L)
  expect_identical(res_int$metrics$bias$status, "FAIL")
  expect_true(res_int$metrics$bias$arithmetic_defined)
  expect_equal(res_int$metrics$bias$raw_bias, as.double(imax), tolerance = 0)
  # Invalid stored computed inputs still refuse, rather than becoming arithmetic diagnostics.
  rec$estimate[1L] <- Inf
  expect_error(evaluate_private_mc_cell(rec, spec, 2000L), "point_computable")
})

test_that("independent audit: entire specification binding and early size guards are retained", {
  spec <- make_mc_spec(n_seeds = 8000L, required_metrics = c("bias", "null_size"))
  spec$premises$calibration_supported_bias <- FALSE
  spec$binding_note <- "fixed opaque source binding"
  rec <- make_mc_records(spec, look = 8000L, n_rejected = 0L)
  prior <- evaluate_private_mc_cell(rec[1:2000, ], spec, 2000L)
  expect_identical(prior$cell_status, "INCONCLUSIVE")
  changed <- spec; changed$binding_note <- "changed opaque binding"
  expect_error(evaluate_private_mc_cell(rec, changed, 8000L, prior), "specification.*identical")
  expect_identical(evaluate_private_mc_cell(rec, spec, 8000L, prior)$cell_status, "INCONCLUSIVE")
  duplicate <- spec; duplicate$premises <- c(duplicate$premises, list(iid_attempts = FALSE))
  expect_error(evaluate_private_mc_cell(rec[1:2000, ], duplicate, 2000L), "named list")
  # Spy clone changes only the scalar scanner to expose guard order; no namespace mutation.
  ce <- new.env(parent = cdid_ns)
  original_scan <- get("mc_is_whole_number", cdid_ns)
  ce$mc_is_whole_number <- function(x) {
    if (length(x) == 2L || length(x) > 8000L) stop("whole_scan_tripwire")
    original_scan(x)
  }
  validator <- get("validate_mc_specification", cdid_ns)
  environment(validator) <- ce
  expect_error(validator(spec, c(2000L, 8000L)), "look.*2000L or 8000L")
  bad_seeds <- spec; bad_seeds$seeds <- seq_len(8001L)
  expect_error(validator(bad_seeds, 2000L), "seeds.*sequence")
})

test_that("D27 conservative asymptotic SE-ratio critical 4.0 endpoints, strict nesting, both-look invariants, truthful metadata, and formula-mutation failure controls hold", {
  z_old <- stats::qnorm(399 / 400)
  z_new <- 4.0
  lim_cont <- 2 * stats::pnorm(4.0) - 1
  expected_conf_type <- "conservative_asymptotic_mc_interval_fixed_critical_4.0_required_calibration_target_99.5pct_limiting_2Phi4_minus_1_not_exact_99.5pct_no_exact_finite_M_guarantee"

  # 1) Default formal argument in mc_continuous_se_ratio_ci is 4.0; bias remains qnorm(399/400)
  expect_identical(eval(formals(mc_continuous_se_ratio_ci)$z), 4.0)
  expect_equal(eval(formals(mc_continuous_bias_ci)$z), z_old, tolerance = 1e-15)

  # 2) Exact endpoints, strict nesting, and log-half-width ratio 4 / qnorm(399/400)
  t_wit <- c(-4 / 5, 1 / 5, 6 / 5)
  s_wit <- c(1, 2, 3)
  r_new <- mc_continuous_se_ratio_ci(t_wit, s_wit, n_attempted = 2000L)
  r_old <- mc_continuous_se_ratio_ci(t_wit, s_wit, n_attempted = 2000L, z = z_old)
  b_wit <- mc_continuous_bias_ci(t_wit, n_attempted = 2000L)

  expect_identical(r_new$z, 4.0)
  expect_equal(r_new$estimate, 2.0, tolerance = 1e-15)
  expect_equal(r_new$log_estimate, log(2.0), tolerance = 1e-15)
  expect_equal(r_new$mcse_log, 1 / 3, tolerance = 1e-15)
  expect_equal(r_new$interval, c(lower = 2 * exp(-4 / 3), upper = 2 * exp(4 / 3)), tolerance = 1e-15)
  expect_equal(r_old$interval, c(lower = 2 * exp(-z_old / 3), upper = 2 * exp(z_old / 3)), tolerance = 1e-15)

  # Strict nesting and half-width ratio on log scale
  expect_true(r_new$interval[["lower"]] < r_old$interval[["lower"]])
  expect_true(r_new$interval[["upper"]] > r_old$interval[["upper"]])
  log_hw_new <- 0.5 * (log(r_new$interval[["upper"]]) - log(r_new$interval[["lower"]]))
  log_hw_old <- 0.5 * (log(r_old$interval[["upper"]]) - log(r_old$interval[["lower"]]))
  expect_equal(log_hw_new / log_hw_old, 4.0 / z_old, tolerance = 1e-15)

  # Unchanged ratio, log_estimate, mean_se, emp_sd, emp_var, psi, psi_bar, if_var, mcse_log, denominator
  for (fld in c("estimate", "log_estimate", "mean_se", "n_zero_se", "t_mean_joint",
                "emp_sd", "emp_var", "psi", "psi_bar", "if_var", "mcse_log",
                "m_computable", "n_attempted", "n_failed", "computability_rate",
                "denominator", "conditioning_event")) {
    expect_equal(r_new[[fld]], r_old[[fld]], tolerance = 1e-15)
  }
  expect_equal(b_wit$z, z_old, tolerance = 1e-15)

  # 3) Truthful method metadata on supported, unsupported, degenerate, nonfinite-error, and plumbing outputs
  check_se_ratio_meta <- function(obj) {
    expect_identical(obj$z, 4.0)
    expect_identical(obj$confidence_type, expected_conf_type)
    expect_equal(obj$required_calibration_target, 199 / 200, tolerance = 1e-15)
    expect_equal(obj$limiting_normal_containment, lim_cont, tolerance = 1e-15)
    expect_true(abs(obj$limiting_normal_containment - 0.995) > 0.0049)
    expect_false(obj$exact_finite_m_guarantee)
    expect_identical(obj$method_metadata$interval_class, "conservative_asymptotic_mc_interval")
    expect_identical(obj$method_metadata$critical_value, 4.0)
    expect_identical(obj$method_metadata$default_critical_value, 4.0)
    expect_equal(obj$method_metadata$required_calibration_target, 199 / 200, tolerance = 1e-15)
    expect_equal(obj$method_metadata$limiting_normal_containment, lim_cont, tolerance = 1e-15)
    expect_false(obj$method_metadata$limiting_containment_equals_target)
    expect_false(obj$method_metadata$exact_finite_m_guarantee)
    expect_false(obj$method_metadata$data_dependent_critical_selection)
    expect_identical(obj$method_metadata$looks_applied, c(2000L, 8000L))
  }
  check_se_ratio_meta(r_new)
  check_se_ratio_meta(mc_continuous_se_ratio_ci(t_wit, s_wit, n_attempted = 2000L, premises_supported = FALSE))
  check_se_ratio_meta(mc_continuous_se_ratio_ci(1.0, 1.0, n_attempted = 2000L))

  nf_fn <- get("mc_nonfinite_error_metric", envir = cdid_ns, inherits = FALSE)
  nf_se <- nf_fn("se_ratio", 2000L, 2000L, list(supported = TRUE, unsupported_reasons = character()))
  nf_bi <- nf_fn("bias", 2000L, 2000L, list(supported = TRUE, unsupported_reasons = character()))
  check_se_ratio_meta(nf_se)
  expect_equal(nf_bi$z, z_old, tolerance = 1e-15)

  plumb <- prepare_mc_calibration_plumbing()
  expect_true(grepl("4\\.0\\*MCSE_logrho", plumb$per_look_continuous_intervals))
  expect_identical(plumb$se_ratio_interval_metadata$critical_value, 4.0)
  expect_false(plumb$se_ratio_interval_metadata$exact_finite_m_guarantee)
  expect_false(plumb$se_ratio_interval_metadata$limiting_containment_equals_target)

  # 4) Both-look use (look = 2000 and look = 8000) where z = 4.0 governs both looks:
  # Choose se_scale so that at look = 2000, legacy qnorm(399/400) would PASS (upper <= 1.05),
  # whereas D27 critical 4.0 straddles 1.05 -> INCONCLUSIVE at look = 2000, and then at look = 8000
  # D27 critical 4.0 shrinks inside [0.95, 1.05] -> PASS!
  spec_2look <- make_mc_spec(n_seeds = 8000L, required_metrics = c("bias", "se_ratio"))
  base_2000 <- make_mc_records(spec_2look, look = 2000L, bias_shift = 0.0, se_scale = 1.0)
  base_r2000 <- mc_continuous_se_ratio_ci(base_2000$estimate - spec_2look$truth, base_2000$se, n_attempted = 2000L)
  mcse_2000 <- base_r2000$mcse_log
  # Pick target rho_hat at look 2000 so log(rho_hat) + z_old * mcse_2000 < log(1.05) < log(rho_hat) + 4.0 * mcse_2000
  target_log_rho_2000 <- log(1.05) - 0.5 * (z_old + 4.0) * mcse_2000
  scale_2000 <- exp(target_log_rho_2000) / base_r2000$estimate

  rec_2000_d27 <- make_mc_records(spec_2look, look = 2000L, bias_shift = 0.0, se_scale = scale_2000)
  legacy_check_2000 <- mc_continuous_se_ratio_ci(
    rec_2000_d27$estimate - spec_2look$truth,
    rec_2000_d27$se,
    n_attempted = 2000L,
    z = z_old
  )
  expect_identical(legacy_check_2000$status, "PASS")

  eval_2000_d27 <- evaluate_private_mc_cell(rec_2000_d27, spec_2look, look = 2000L)
  expect_identical(eval_2000_d27$metrics$bias$status, "PASS")
  expect_identical(eval_2000_d27$metrics$se_ratio$z, 4.0)
  expect_identical(eval_2000_d27$metrics$se_ratio$status, "INCONCLUSIVE")
  expect_identical(eval_2000_d27$cell_status, "INCONCLUSIVE")
  expect_true(eval_2000_d27$can_extend_to_8000)

  rec_8000_d27 <- make_mc_records(spec_2look, look = 8000L, bias_shift = 0.0, se_scale = 1.0)
  rec_8000_d27[1:2000, ] <- rec_2000_d27
  eval_8000_d27 <- evaluate_private_mc_cell(rec_8000_d27, spec_2look, look = 8000L,
                                            prior_assessment = eval_2000_d27)
  expect_identical(eval_8000_d27$metrics$se_ratio$z, 4.0)
  expect_identical(eval_8000_d27$metrics$se_ratio$decision_look, 8000L)
  expect_true(eval_8000_d27$metrics$se_ratio$reevaluated_at_8000)
  expect_identical(eval_8000_d27$metrics$se_ratio$status, "PASS")
  expect_identical(eval_8000_d27$cell_status, "PASS")

  # 5) Formula-mutation failure controls: any mutation to z, m-1 divisor, psi_logrho derivatives,
  # MCSE denominator m, or bias critical value MUST fail our reference checks.
  t_mut <- rep(c(-1, -1, 0, 0, 0, 2), 200L) + 0.1
  s_mut <- rep(1 + (c(-1, -1, 0, 0, 0, 2)^2 - 1) / 4, 200L)
  tru_r <- mc_continuous_se_ratio_ci(t_mut, s_mut, n_attempted = 2000L)
  tru_b <- mc_continuous_bias_ci(t_mut, n_attempted = 2000L)

  # Mutation A: reverting se_ratio critical to qnorm(399/400)
  mut_a_ci <- exp(tru_r$log_estimate + c(-1, 1) * z_old * tru_r$mcse_log)
  expect_true(max(abs(unname(tru_r$interval) - mut_a_ci)) > 1e-2)

  # Mutation B: using m instead of m - 1 in empirical SD or IF variance
  m_mut <- length(t_mut)
  mut_b_rho <- mean(s_mut) / sqrt(sum((t_mut - mean(t_mut))^2) / m_mut)
  expect_true(abs(tru_r$estimate - mut_b_rho) > 1e-4)

  # Mutation C: omitting mean-SE derivative in psi_logrho
  s2_mut <- stats::var(t_mut)
  mut_c_psi <- -((t_mut - mean(t_mut))^2 - s2_mut) / (2 * s2_mut)
  mut_c_mcse <- sqrt(stats::var(mut_c_psi) / m_mut)
  expect_true(abs(tru_r$mcse_log - mut_c_mcse) > 1e-3)

  # Mutation D: using n_attempted (2000) instead of m_computable (1200) in MCSE_log
  mut_d_mcse <- sqrt(tru_r$if_var / 2000)
  expect_true(abs(tru_r$mcse_log - mut_d_mcse) > 1e-3)

  # Mutation E: accidentally changing bias critical value to 4.0
  mut_e_ci <- tru_b$estimate + c(-1, 1) * 4.0 * tru_b$mcse
  expect_true(max(abs(unname(tru_b$interval) - mut_e_ci)) > 1e-2)
})
