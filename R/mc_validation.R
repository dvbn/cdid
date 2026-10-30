# Private Monte Carlo validation evaluator under approved D7/D14/D15/D16 and D19/D27 SE-ratio amendment.
# Unexported internal implementation for cdid; no public API, RNG draws or fits.

mc_allowed_metrics <- function() {
  c("coverage", "band_coverage", "null_size", "bias", "se_ratio")
}

mc_allowed_looks <- function() {
  c(2000L, 8000L)
}

mc_is_whole_number <- function(x) {
  is.numeric(x) && !anyNA(x) && all(is.finite(x)) &&
    all(x == round(x)) && all(abs(x) <= .Machine$integer.max)
}

mc_is_scalar_flag <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

mc_is_scalar_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

mc_is_scalar_finite_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

mc_clopper_pearson_ci <- function(k, d, alpha_tail = 1 / 400) {
  if (length(k) != 1L || !mc_is_whole_number(k) || k < 0) {
    stop("`k` must be a single non-negative whole number.", call. = FALSE)
  }
  if (length(d) != 1L || !mc_is_whole_number(d) || d < 0) {
    stop("`d` must be a single non-negative whole number.", call. = FALSE)
  }
  if (!mc_is_scalar_finite_number(alpha_tail) ||
      alpha_tail <= 0 || alpha_tail >= 0.5) {
    stop("`alpha_tail` must be a single finite number in (0, 0.5).",
         call. = FALSE)
  }
  k <- as.integer(k)
  d <- as.integer(d)
  if (d == 0L) {
    if (k != 0L) {
      stop("`k` cannot exceed `d` when `d` is zero.", call. = FALSE)
    }
    return(c(lower = NA_real_, upper = NA_real_))
  }
  if (k > d) {
    stop("`k` must not exceed `d`.", call. = FALSE)
  }
  lower <- if (k == 0L) 0 else stats::qbeta(alpha_tail, k, d - k + 1L)
  upper <- if (k == d) 1 else stats::qbeta(1 - alpha_tail, k + 1L, d - k)
  c(lower = lower, upper = upper)
}

mc_wilson_ci <- function(k, d, conf_level = 0.99) {
  if (length(k) != 1L || !mc_is_whole_number(k) || k < 0) {
    stop("`k` must be a single non-negative whole number.", call. = FALSE)
  }
  if (length(d) != 1L || !mc_is_whole_number(d) || d < 0) {
    stop("`d` must be a single non-negative whole number.", call. = FALSE)
  }
  if (!mc_is_scalar_finite_number(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single finite number in (0, 1).",
         call. = FALSE)
  }
  k <- as.integer(k)
  d <- as.integer(d)
  if (d == 0L) {
    if (k != 0L) {
      stop("`k` cannot exceed `d` when `d` is zero.", call. = FALSE)
    }
    return(c(lower = NA_real_, upper = NA_real_))
  }
  if (k > d) {
    stop("`k` must not exceed `d`.", call. = FALSE)
  }
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  p <- k / d
  z2_d <- (z * z) / d
  denom <- 1 + z2_d
  center <- (p + 0.5 * z2_d) / denom
  half <- (z * sqrt((p * (1 - p)) / d + (z2_d / (4 * d)))) / denom
  c(lower = max(0, center - half), upper = min(1, center + half))
}

mc_classify_interval_decision <- function(metric, lower, upper, nominal = NULL,
                                          supported = TRUE,
                                          unsupported_reasons = character()) {
  if (!mc_is_scalar_string(metric) || !(metric %in% mc_allowed_metrics())) {
    stop("`metric` must be a single supported metric name.", call. = FALSE)
  }
  if (!mc_is_scalar_flag(supported)) {
    stop("`supported` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(unsupported_reasons) || anyNA(unsupported_reasons)) {
    stop("`unsupported_reasons` must be a character vector without NA.",
         call. = FALSE)
  }
  if (!supported || length(unsupported_reasons) > 0L ||
      !mc_is_scalar_finite_number(lower) ||
      !mc_is_scalar_finite_number(upper) || lower > upper) {
    return("INCONCLUSIVE")
  }
  if (metric %in% c("coverage", "band_coverage")) {
    if (!mc_is_scalar_finite_number(nominal) || nominal <= 0 || nominal >= 1) {
      stop("`nominal` coverage level must be in (0, 1).", call. = FALSE)
    }
    floor_val <- nominal - 0.01
    # D7 strict lower-bound rule: PASS if L > nominal - 0.01, FAIL if U < nominal - 0.01,
    # equality at either boundary is INCONCLUSIVE.
    if (lower > floor_val) {
      return("PASS")
    }
    if (upper < floor_val) {
      return("FAIL")
    }
    return("INCONCLUSIVE")
  }
  if (metric == "null_size") {
    if (!mc_is_scalar_finite_number(nominal) || nominal <= 0 || nominal >= 1) {
      stop("`nominal` size level must be in (0, 1).", call. = FALSE)
    }
    ceiling_val <- nominal + 0.01
    # D7/D15 null size rule: accuracy region [0, nominal + 0.01].
    # PASS if U <= nominal + 0.01; FAIL if L > nominal + 0.01; overlap INCONCLUSIVE.
    if (upper <= ceiling_val) {
      return("PASS")
    }
    if (lower > ceiling_val) {
      return("FAIL")
    }
    return("INCONCLUSIVE")
  }
  if (metric == "bias") {
    # D7/D16 signed normalized bias inclusive region [-0.10, 0.10].
    if (lower >= -0.10 && upper <= 0.10) {
      return("PASS")
    }
    if (lower > 0.10 || upper < -0.10) {
      return("FAIL")
    }
    return("INCONCLUSIVE")
  }
  if (metric == "se_ratio") {
    # D7/D16 SE ratio inclusive region [0.95, 1.05].
    if (lower >= 0.95 && upper <= 1.05) {
      return("PASS")
    }
    if (lower > 1.05 || upper < 0.95) {
      return("FAIL")
    }
    return("INCONCLUSIVE")
  }
  stop("Unhandled metric in `mc_classify_interval_decision`.", call. = FALSE)
}

mc_continuous_bias_ci <- function(t_errors, n_attempted = length(t_errors),
                                  premises_supported = TRUE,
                                  unsupported_reasons = character(),
                                  z = stats::qnorm(399 / 400)) {
  if (!is.numeric(t_errors) || any(!is.finite(t_errors))) {
    stop("`t_errors` must be a numeric vector of finite values on D_point.",
         call. = FALSE)
  }
  if (length(n_attempted) != 1L || !mc_is_whole_number(n_attempted) ||
      n_attempted < length(t_errors)) {
    stop("`n_attempted` must be a single whole number >= length(t_errors).",
         call. = FALSE)
  }
  if (!mc_is_scalar_flag(premises_supported)) {
    stop("`premises_supported` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(unsupported_reasons) || anyNA(unsupported_reasons)) {
    stop("`unsupported_reasons` must be a character vector without NA.",
         call. = FALSE)
  }
  if (!mc_is_scalar_finite_number(z) || z <= 0) {
    stop("`z` must be a single positive finite number.", call. = FALSE)
  }

  m <- length(t_errors)
  n_att <- as.integer(n_attempted)
  reasons <- unique(unsupported_reasons)
  if (!premises_supported && length(reasons) == 0L) {
    reasons <- c(reasons, "required_external_premises_unsupported")
  }

  empty_out <- function(reasons_vec, tbar_val = NA_real_, s_val = NA_real_,
                        s2_val = NA_real_, beta_val = NA_real_,
                        psi_vec = numeric(), psibar_val = NA_real_,
                        vhat_val = NA_real_, mcse_val = NA_real_,
                        desc_ci = c(lower = NA_real_, upper = NA_real_),
                        arith_ok = FALSE) {
    list(
      metric = "bias",
      status = "INCONCLUSIVE",
      supported = FALSE,
      premises_supported = premises_supported,
      arithmetic_defined = arith_ok,
      unsupported_reasons = unique(reasons_vec),
      estimate = beta_val,
      raw_bias = tbar_val,
      emp_sd = s_val,
      emp_var = s2_val,
      psi = psi_vec,
      psi_bar = psibar_val,
      if_var = vhat_val,
      mcse = mcse_val,
      z = z,
      interval = c(lower = NA_real_, upper = NA_real_),
      descriptive_interval = desc_ci,
      threshold = c(lower = -0.10, upper = 0.10),
      region_inclusive = c(lower = TRUE, upper = TRUE),
      m_computable = m,
      n_attempted = n_att,
      n_failed = n_att - m,
      computability_rate = if (n_att > 0L) m / n_att else NA_real_,
      denominator = paste0(m, "/", n_att),
      conditioning_event = "point_computable",
      confidence_type = "asymptotic_99.5pct_per_look_delta_method_not_exact_finite_99pct"
    )
  }

  if (m < 2L) {
    return(empty_out(c(reasons, "m_lt_2")))
  }

  tbar <- sum(t_errors) / m
  centered_t <- t_errors - tbar
  s2 <- sum(centered_t * centered_t) / (m - 1L)
  if (!is.finite(s2) || s2 <= 0) {
    return(empty_out(c(reasons, "zero_or_nonpositive_empirical_sd"),
                     tbar_val = tbar, s2_val = s2))
  }
  s <- sqrt(s2)
  beta_hat <- tbar / s
  if (!is.finite(s) || s <= 0 || !is.finite(beta_hat)) {
    return(empty_out(c(reasons, "nonfinite_standardized_bias"),
                     tbar_val = tbar, s_val = s, s2_val = s2,
                     beta_val = beta_hat))
  }

  psi_beta <- centered_t / s - beta_hat * ((centered_t * centered_t - s2) / (2 * s2))
  if (any(!is.finite(psi_beta))) {
    return(empty_out(c(reasons, "nonfinite_plugin_influence"),
                     tbar_val = tbar, s_val = s, s2_val = s2,
                     beta_val = beta_hat, psi_vec = psi_beta))
  }
  psibar <- sum(psi_beta) / m
  centered_psi <- psi_beta - psibar
  v_hat <- sum(centered_psi * centered_psi) / (m - 1L)
  if (!is.finite(v_hat) || v_hat <= 0) {
    return(empty_out(c(reasons, "zero_or_nonpositive_influence_variance"),
                     tbar_val = tbar, s_val = s, s2_val = s2,
                     beta_val = beta_hat, psi_vec = psi_beta,
                     psibar_val = psibar, vhat_val = v_hat))
  }

  mcse <- sqrt(v_hat / m)
  ci_raw <- c(lower = beta_hat - z * mcse, upper = beta_hat + z * mcse)
  if (!is.finite(mcse) || mcse <= 0 || any(!is.finite(ci_raw))) {
    return(empty_out(c(reasons, "nonfinite_ci_endpoints"),
                     tbar_val = tbar, s_val = s, s2_val = s2,
                     beta_val = beta_hat, psi_vec = psi_beta,
                     psibar_val = psibar, vhat_val = v_hat, mcse_val = mcse))
  }

  if (!premises_supported || length(reasons) > 0L) {
    return(empty_out(reasons,
                     tbar_val = tbar, s_val = s, s2_val = s2,
                     beta_val = beta_hat, psi_vec = psi_beta,
                     psibar_val = psibar, vhat_val = v_hat, mcse_val = mcse,
                     desc_ci = ci_raw, arith_ok = TRUE))
  }

  status <- mc_classify_interval_decision("bias", ci_raw[["lower"]],
                                          ci_raw[["upper"]],
                                          supported = TRUE,
                                          unsupported_reasons = character())
  list(
    metric = "bias",
    status = status,
    supported = TRUE,
    premises_supported = TRUE,
    arithmetic_defined = TRUE,
    unsupported_reasons = character(),
    estimate = beta_hat,
    raw_bias = tbar,
    emp_sd = s,
    emp_var = s2,
    psi = psi_beta,
    psi_bar = psibar,
    if_var = v_hat,
    mcse = mcse,
    z = z,
    interval = ci_raw,
    descriptive_interval = ci_raw,
    threshold = c(lower = -0.10, upper = 0.10),
    region_inclusive = c(lower = TRUE, upper = TRUE),
    m_computable = m,
    n_attempted = n_att,
    n_failed = n_att - m,
    computability_rate = if (n_att > 0L) m / n_att else NA_real_,
    denominator = paste0(m, "/", n_att),
    conditioning_event = "point_computable",
    confidence_type = "asymptotic_99.5pct_per_look_delta_method_not_exact_finite_99pct"
  )
}

mc_continuous_se_ratio_ci <- function(t_errors, se_values,
                                      n_attempted = length(t_errors),
                                      premises_supported = TRUE,
                                      unsupported_reasons = character(),
                                      z = 4.0) {
  if (!is.numeric(t_errors) || any(!is.finite(t_errors))) {
    stop("`t_errors` must be a numeric vector of finite values on D_point & D_se.",
         call. = FALSE)
  }
  if (!is.numeric(se_values) || length(se_values) != length(t_errors) ||
      any(!is.finite(se_values)) || any(se_values < 0)) {
    stop("`se_values` must be a finite non-negative numeric vector matching `t_errors`.",
         call. = FALSE)
  }
  if (length(n_attempted) != 1L || !mc_is_whole_number(n_attempted) ||
      n_attempted < length(t_errors)) {
    stop("`n_attempted` must be a single whole number >= length(t_errors).",
         call. = FALSE)
  }
  if (!mc_is_scalar_flag(premises_supported)) {
    stop("`premises_supported` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(unsupported_reasons) || anyNA(unsupported_reasons)) {
    stop("`unsupported_reasons` must be a character vector without NA.",
         call. = FALSE)
  }
  if (!mc_is_scalar_finite_number(z) || z <= 0) {
    stop("`z` must be a single positive finite number.", call. = FALSE)
  }

  m <- length(t_errors)
  n_att <- as.integer(n_attempted)
  n_zero_se <- sum(se_values == 0)
  reasons <- unique(unsupported_reasons)
  if (!premises_supported && length(reasons) == 0L) {
    reasons <- c(reasons, "required_external_premises_unsupported")
  }
  req_cal_target <- 199 / 200
  lim_containment <- 2 * stats::pnorm(z) - 1
  conf_type <- "conservative_asymptotic_mc_interval_fixed_critical_4.0_required_calibration_target_99.5pct_limiting_2Phi4_minus_1_not_exact_99.5pct_no_exact_finite_M_guarantee"
  method_meta <- list(
    interval_class = "conservative_asymptotic_mc_interval",
    critical_value = z,
    default_critical_value = 4.0,
    required_calibration_target = req_cal_target,
    limiting_normal_containment = lim_containment,
    limiting_containment_equals_target = FALSE,
    exact_finite_m_guarantee = FALSE,
    data_dependent_critical_selection = FALSE,
    looks_applied = c(2000L, 8000L)
  )

  empty_out <- function(reasons_vec, tbar_val = NA_real_, sbar_val = NA_real_,
                        s_val = NA_real_, s2_val = NA_real_,
                        rho_val = NA_real_, log_rho_val = NA_real_,
                        psi_vec = numeric(), psibar_val = NA_real_,
                        vhat_val = NA_real_, mcse_val = NA_real_,
                        desc_ci = c(lower = NA_real_, upper = NA_real_),
                        arith_ok = FALSE) {
    list(
      metric = "se_ratio",
      status = "INCONCLUSIVE",
      supported = FALSE,
      premises_supported = premises_supported,
      arithmetic_defined = arith_ok,
      unsupported_reasons = unique(reasons_vec),
      estimate = rho_val,
      log_estimate = log_rho_val,
      mean_se = sbar_val,
      n_zero_se = n_zero_se,
      t_mean_joint = tbar_val,
      emp_sd = s_val,
      emp_var = s2_val,
      psi = psi_vec,
      psi_bar = psibar_val,
      if_var = vhat_val,
      mcse_log = mcse_val,
      z = z,
      interval = c(lower = NA_real_, upper = NA_real_),
      descriptive_interval = desc_ci,
      threshold = c(lower = 0.95, upper = 1.05),
      region_inclusive = c(lower = TRUE, upper = TRUE),
      m_computable = m,
      n_attempted = n_att,
      n_failed = n_att - m,
      computability_rate = if (n_att > 0L) m / n_att else NA_real_,
      denominator = paste0(m, "/", n_att),
      conditioning_event = "point_computable AND se_computable",
      confidence_type = conf_type,
      required_calibration_target = req_cal_target,
      limiting_normal_containment = lim_containment,
      exact_finite_m_guarantee = FALSE,
      method_metadata = method_meta
    )
  }

  if (m < 2L) {
    return(empty_out(c(reasons, "m_lt_2")))
  }

  tbar <- sum(t_errors) / m
  sbar <- sum(se_values) / m
  centered_t <- t_errors - tbar
  s2 <- sum(centered_t * centered_t) / (m - 1L)

  bad_reasons <- character()
  if (!is.finite(s2) || s2 <= 0) {
    bad_reasons <- c(bad_reasons, "zero_or_nonpositive_empirical_sd")
  }
  if (!is.finite(sbar) || sbar <= 0) {
    bad_reasons <- c(bad_reasons, "nonpositive_or_nonfinite_mean_se")
  }
  if (length(bad_reasons) > 0L) {
    return(empty_out(c(reasons, bad_reasons),
                     tbar_val = tbar, sbar_val = sbar, s2_val = s2))
  }

  s <- sqrt(s2)
  rho_hat <- sbar / s
  log_rho_hat <- log(rho_hat)
  if (!is.finite(s) || s <= 0 || !is.finite(rho_hat) || rho_hat <= 0 ||
      !is.finite(log_rho_hat)) {
    return(empty_out(c(reasons, "nonfinite_se_ratio"),
                     tbar_val = tbar, sbar_val = sbar, s_val = s, s2_val = s2,
                     rho_val = rho_hat, log_rho_val = log_rho_hat))
  }

  psi_logrho <- (2 * s2 * (se_values - sbar) - sbar * (centered_t * centered_t - s2)) / (2 * sbar * s2)
  if (any(!is.finite(psi_logrho))) {
    return(empty_out(c(reasons, "nonfinite_plugin_influence"),
                     tbar_val = tbar, sbar_val = sbar, s_val = s, s2_val = s2,
                     rho_val = rho_hat, log_rho_val = log_rho_hat,
                     psi_vec = psi_logrho))
  }
  psibar <- sum(psi_logrho) / m
  centered_psi <- psi_logrho - psibar
  v_hat <- sum(centered_psi * centered_psi) / (m - 1L)
  if (!is.finite(v_hat) || v_hat <= 0) {
    return(empty_out(c(reasons, "zero_or_nonpositive_influence_variance"),
                     tbar_val = tbar, sbar_val = sbar, s_val = s, s2_val = s2,
                     rho_val = rho_hat, log_rho_val = log_rho_hat,
                     psi_vec = psi_logrho, psibar_val = psibar,
                     vhat_val = v_hat))
  }

  mcse_log <- sqrt(v_hat / m)
  ci_raw <- c(
    lower = exp(log_rho_hat - z * mcse_log),
    upper = exp(log_rho_hat + z * mcse_log)
  )
  if (!is.finite(mcse_log) || mcse_log <= 0 || any(!is.finite(ci_raw))) {
    return(empty_out(c(reasons, "nonfinite_ci_endpoints"),
                     tbar_val = tbar, sbar_val = sbar, s_val = s, s2_val = s2,
                     rho_val = rho_hat, log_rho_val = log_rho_hat,
                     psi_vec = psi_logrho, psibar_val = psibar,
                     vhat_val = v_hat, mcse_val = mcse_log))
  }

  if (!premises_supported || length(reasons) > 0L) {
    return(empty_out(reasons,
                     tbar_val = tbar, sbar_val = sbar, s_val = s, s2_val = s2,
                     rho_val = rho_hat, log_rho_val = log_rho_hat,
                     psi_vec = psi_logrho, psibar_val = psibar,
                     vhat_val = v_hat, mcse_val = mcse_log,
                     desc_ci = ci_raw, arith_ok = TRUE))
  }

  status <- mc_classify_interval_decision("se_ratio", ci_raw[["lower"]],
                                          ci_raw[["upper"]],
                                          supported = TRUE,
                                          unsupported_reasons = character())
  list(
    metric = "se_ratio",
    status = status,
    supported = TRUE,
    premises_supported = TRUE,
    arithmetic_defined = TRUE,
    unsupported_reasons = character(),
    estimate = rho_hat,
    log_estimate = log_rho_hat,
    mean_se = sbar,
    n_zero_se = n_zero_se,
    t_mean_joint = tbar,
    emp_sd = s,
    emp_var = s2,
    psi = psi_logrho,
    psi_bar = psibar,
    if_var = v_hat,
    mcse_log = mcse_log,
    z = z,
    interval = ci_raw,
    descriptive_interval = ci_raw,
    threshold = c(lower = 0.95, upper = 1.05),
    region_inclusive = c(lower = TRUE, upper = TRUE),
    m_computable = m,
    n_attempted = n_att,
    n_failed = n_att - m,
    computability_rate = if (n_att > 0L) m / n_att else NA_real_,
    denominator = paste0(m, "/", n_att),
    conditioning_event = "point_computable AND se_computable",
    confidence_type = conf_type,
    required_calibration_target = req_cal_target,
    limiting_normal_containment = lim_containment,
    exact_finite_m_guarantee = FALSE,
    method_metadata = method_meta
  )
}

mc_resolve_premise_flag <- function(premises, specific_key, fallback_key = NULL) {
  val <- premises[[specific_key]]
  if (is.null(val) && !is.null(fallback_key)) {
    val <- premises[[fallback_key]]
  }
  if (is.null(val)) {
    return(FALSE)
  }
  if (!mc_is_scalar_flag(val)) {
    stop(sprintf("Premise `%s` must be TRUE or FALSE when supplied.",
                 specific_key), call. = FALSE)
  }
  val
}

mc_check_metric_premises <- function(metric, specification) {
  premises <- specification$premises
  missing_or_false <- character()

  check_flag <- function(label, specific_key, fallback_key = NULL) {
    ok <- mc_resolve_premise_flag(premises, specific_key, fallback_key)
    if (!ok) {
      missing_or_false <<- c(missing_or_false, label)
    }
  }

  check_flag("iid_attempts", "iid_attempts")
  check_flag("fixed_predeclared_events", "fixed_predeclared_events")
  check_flag("no_hidden_peeks", "no_hidden_peeks")

  if (metric == "coverage") {
    # All-attempt valid-and-covered binomial interval requires iid attempts,
    # fixed predeclared events and no hidden peeks.
  } else if (metric == "band_coverage") {
    # Same binomial requirements for simultaneous band valid-and-covered.
  } else if (metric == "null_size") {
    null_ok <- FALSE
    if (!is.null(specification$null_premise)) {
      if (!mc_is_scalar_flag(specification$null_premise)) {
        stop("`specification$null_premise` must be TRUE or FALSE.",
             call. = FALSE)
      }
      null_ok <- specification$null_premise
    } else {
      null_ok <- mc_resolve_premise_flag(premises, "null_holds", "null_premise")
    }
    if (!null_ok) {
      missing_or_false <- c(missing_or_false, "null_premise_not_declared_true")
    }
    check_flag("positive_test_event_probability",
               "positive_test_event_probability",
               "positive_event_probability")
  } else if (metric == "bias") {
    check_flag("positive_point_event_probability",
               "positive_point_event_probability",
               "positive_event_probability")
    check_flag("positive_conditional_sd_point",
               "positive_conditional_sd_point",
               "positive_conditional_sd")
    check_flag("conditional_fourth_moment_T_point",
               "conditional_fourth_moment_T_point",
               "conditional_fourth_moment_T")
    check_flag("nondegenerate_if_variance_bias",
               "nondegenerate_if_variance_bias",
               "nondegenerate_if_variance")
    check_flag("calibration_supported_bias",
               "calibration_supported_bias",
               "calibration_supported")
  } else if (metric == "se_ratio") {
    check_flag("positive_joint_event_probability",
               "positive_joint_event_probability",
               "positive_event_probability")
    check_flag("positive_conditional_sd_joint",
               "positive_conditional_sd_joint",
               "positive_conditional_sd")
    check_flag("conditional_fourth_moment_T_joint",
               "conditional_fourth_moment_T_joint",
               "conditional_fourth_moment_T")
    check_flag("conditional_second_moment_S_joint",
               "conditional_second_moment_S_joint",
               "conditional_second_moment_S")
    check_flag("positive_conditional_mean_se_joint",
               "positive_conditional_mean_se_joint",
               "positive_conditional_mean_se")
    check_flag("nondegenerate_if_variance_se_ratio",
               "nondegenerate_if_variance_se_ratio",
               "nondegenerate_if_variance")
    check_flag("calibration_supported_se_ratio",
               "calibration_supported_se_ratio",
               "calibration_supported")
  }

  list(
    supported = length(missing_or_false) == 0L,
    unsupported_reasons = unique(missing_or_false)
  )
}

validate_mc_specification <- function(specification, look) {
  if (length(look) != 1L || !mc_is_whole_number(look) ||
      !(as.integer(look) %in% mc_allowed_looks())) {
    stop("`look` must be 2000L or 8000L; intermediate or development lengths are not eligible for acceptance assessment.",
         call. = FALSE)
  }
  look <- as.integer(look)

  if (!is.list(specification) || is.data.frame(specification) ||
      is.null(names(specification)) || anyNA(names(specification)) ||
      any(!nzchar(names(specification))) || anyDuplicated(names(specification))) {
    stop("`specification` must be a named list.", call. = FALSE)
  }
  if (!mc_is_scalar_string(specification$cell_id)) {
    stop("`specification$cell_id` must be a single non-empty character string.",
         call. = FALSE)
  }
  if (!mc_is_scalar_string(specification$target)) {
    stop("`specification$target` must be a single non-empty character string.",
         call. = FALSE)
  }
  if (!mc_is_scalar_finite_number(specification$truth)) {
    stop("`specification$truth` must be a single finite numeric value.",
         call. = FALSE)
  }

  seeds <- specification$seeds
  if (!(length(seeds) %in% c(look, 8000L)) ||
      length(seeds) < look || !mc_is_whole_number(seeds) ||
      anyDuplicated(as.integer(seeds))) {
    stop("`specification$seeds` must be a duplicate-free whole-number sequence of length `look` (or 8000).",
         call. = FALSE)
  }

  req_metrics <- specification$required_metrics
  if (!is.character(req_metrics) || length(req_metrics) < 1L ||
      anyNA(req_metrics) || anyDuplicated(req_metrics) ||
      !all(req_metrics %in% mc_allowed_metrics())) {
    stop("`specification$required_metrics` must be a non-empty duplicate-free subset of allowed metrics.",
         call. = FALSE)
  }

  if (any(c("coverage", "band_coverage") %in% req_metrics) ||
      !is.null(specification$nominal_coverage)) {
    if (!mc_is_scalar_finite_number(specification$nominal_coverage) ||
        specification$nominal_coverage <= 0 ||
        specification$nominal_coverage >= 1) {
      stop("`specification$nominal_coverage` must be a single finite number in (0, 1).",
           call. = FALSE)
    }
  }

  if (("null_size" %in% req_metrics) || !is.null(specification$nominal_size)) {
    if (!mc_is_scalar_finite_number(specification$nominal_size) ||
        specification$nominal_size <= 0 ||
        specification$nominal_size >= 1) {
      stop("`specification$nominal_size` must be a single finite number in (0, 1).",
           call. = FALSE)
    }
  }

  if (!is.null(specification$null_premise) &&
      !mc_is_scalar_flag(specification$null_premise)) {
    stop("`specification$null_premise` must be TRUE, FALSE, or NULL.",
         call. = FALSE)
  }

  required_event_keys <- c("point", "se", "interval", "band", "test")
  events <- specification$events
  if (!is.list(events) || anyNA(names(events)) ||
      any(!nzchar(names(events))) || anyDuplicated(names(events)) ||
      !all(required_event_keys %in% names(events))) {
    stop("`specification$events` must be a named list declaring separate `point`, `se`, `interval`, `band`, and `test` event column names.",
         call. = FALSE)
  }
  event_cols <- vapply(required_event_keys, function(k) {
    val <- events[[k]]
    if (!mc_is_scalar_string(val)) {
      stop(sprintf("`specification$events$%s` must be a single non-empty column name.", k),
           call. = FALSE)
    }
    val
  }, character(1L))
  if (anyDuplicated(event_cols)) {
    stop("Point, SE, interval, band and test computability events must be separate fields.",
         call. = FALSE)
  }

  if (!is.list(specification$premises) || is.data.frame(specification$premises) ||
      (length(specification$premises) > 0L &&
       (is.null(names(specification$premises)) || anyNA(names(specification$premises)) ||
        any(!nzchar(names(specification$premises))) || anyDuplicated(names(specification$premises))))) {
    stop("`specification$premises` must be a named list of declared mathematical and calibration premises.",
         call. = FALSE)
  }
  for (p_name in names(specification$premises)) {
    p_val <- specification$premises[[p_name]]
    if (!is.null(p_val) && !mc_is_scalar_flag(p_val)) {
      stop(sprintf("Premise `%s` in `specification$premises` must be TRUE or FALSE.",
                   p_name), call. = FALSE)
    }
  }

  if (!is.null(specification$provenance) &&
      (!is.list(specification$provenance) && !mc_is_scalar_string(specification$provenance))) {
    stop("`specification$provenance` must be a list or non-empty character string when supplied.",
         call. = FALSE)
  }

  invisible(TRUE)
}

validate_mc_records <- function(records, specification, look) {
  validate_mc_specification(specification, look)
  look <- as.integer(look)

  if (!is.data.frame(records)) {
    stop("`records` must be a data.frame.", call. = FALSE)
  }
  if (nrow(records) != look) {
    stop(sprintf("`nrow(records)` (%d) must equal `look` (%d).",
                 nrow(records), look), call. = FALSE)
  }

  event_cols <- vapply(c("point", "se", "interval", "band", "test"),
                       function(k) specification$events[[k]],
                       character(1L))
  required_cols <- c("attempt", "seed", "estimate", "se", "covered",
                     "band_covered", "rejected", unname(event_cols),
                     "failure_type", "failure_message", "n_warnings",
                     "warnings")
  missing_cols <- setdiff(required_cols, names(records))
  if (length(missing_cols) > 0L) {
    stop(sprintf("`records` is missing required columns: %s.",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Each required field supplies one scalar per attempted replication.
  # Matrix/array columns otherwise recycle logical subsets and change denominators.
  if (anyDuplicated(names(records))) {
    stop("`records` must have unique column names.", call. = FALSE)
  }
  scalar_cols <- unique(c(required_cols, intersect(names(records),
    c("rep", "cell_id", "target", "truth", "nominal_coverage", "nominal_size",
      "lower", "upper", "band_lower", "band_upper"))))
  for (col_nm in scalar_cols) {
    value <- records[[col_nm]]
    if (!is.atomic(value) || !is.null(dim(value)) || length(value) != look) {
      stop(sprintf("Record column `%s` must be a scalar-per-attempt vector of length look, without dimensions.",
                   col_nm), call. = FALSE)
    }
  }

  # Attempt order and completeness
  if (!mc_is_whole_number(records$attempt) ||
      !identical(as.integer(records$attempt), seq_len(look))) {
    stop("`records$attempt` must be strictly ordered 1..look with no missing, duplicated, or reordered attempts.",
         call. = FALSE)
  }
  if ("rep" %in% names(records)) {
    if (!mc_is_whole_number(records$rep) ||
        !identical(as.integer(records$rep), seq_len(look))) {
      stop("`records$rep` must match 1..look when present.", call. = FALSE)
    }
  }

  # Seed binding
  expected_seeds <- as.integer(specification$seeds[seq_len(look)])
  if (!mc_is_whole_number(records$seed) ||
      anyDuplicated(as.integer(records$seed)) ||
      !identical(as.integer(records$seed), expected_seeds)) {
    stop("`records$seed` must match `specification$seeds[1:look]` in exact order without duplicates.",
         call. = FALSE)
  }

  # Optional cell/target/truth/nominal binding columns in records
  if ("cell_id" %in% names(records)) {
    if (!is.character(records$cell_id) || anyNA(records$cell_id) ||
        !all(records$cell_id == specification$cell_id)) {
      stop("`records$cell_id` does not match `specification$cell_id`.",
           call. = FALSE)
    }
  }
  if ("target" %in% names(records)) {
    if (!is.character(records$target) || anyNA(records$target) ||
        !all(records$target == specification$target)) {
      stop("`records$target` does not match `specification$target`.",
           call. = FALSE)
    }
  }
  if ("truth" %in% names(records)) {
    if (!is.numeric(records$truth) || any(!is.finite(records$truth)) ||
        !all(records$truth == specification$truth)) {
      stop("`records$truth` does not match `specification$truth`.",
           call. = FALSE)
    }
  }
  if ("nominal_coverage" %in% names(records) &&
      !is.null(specification$nominal_coverage)) {
    if (!is.numeric(records$nominal_coverage) ||
        any(!is.finite(records$nominal_coverage)) ||
        !all(records$nominal_coverage == specification$nominal_coverage)) {
      stop("`records$nominal_coverage` does not match `specification$nominal_coverage`.",
           call. = FALSE)
    }
  }
  if ("nominal_size" %in% names(records) &&
      !is.null(specification$nominal_size)) {
    if (!is.numeric(records$nominal_size) ||
        any(!is.finite(records$nominal_size)) ||
        !all(records$nominal_size == specification$nominal_size)) {
      stop("`records$nominal_size` does not match `specification$nominal_size`.",
           call. = FALSE)
    }
  }

  # Validate event columns are logical without NA
  for (ev_key in names(event_cols)) {
    col_nm <- event_cols[[ev_key]]
    ev_vec <- records[[col_nm]]
    if (!is.logical(ev_vec) || anyNA(ev_vec)) {
      stop(sprintf("Event column `%s` (%s) must be logical without NA.",
                   col_nm, ev_key), call. = FALSE)
    }
  }

  d_point <- records[[event_cols[["point"]]]]
  d_se <- records[[event_cols[["se"]]]]
  d_interval <- records[[event_cols[["interval"]]]]
  d_band <- records[[event_cols[["band"]]]]
  d_test <- records[[event_cols[["test"]]]]

  # Validate value types and computed-flag consistency (never silently turn invalid values into failures)
  if (!is.numeric(records$estimate)) {
    stop("`records$estimate` must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(records$estimate[d_point]))) {
    stop("Rows with `point_computable == TRUE` must have finite `estimate`.",
         call. = FALSE)
  }

  if (!is.numeric(records$se)) {
    stop("`records$se` must be numeric.", call. = FALSE)
  }
  # Retain finite zero reported SE if declared computed; negative or non-finite SE is invalid input.
  if (any(!is.finite(records$se[d_se]) | records$se[d_se] < 0)) {
    stop("Rows with `se_computable == TRUE` must have finite non-negative `se` (finite zero SE is retained; negative/non-finite SE is invalid).",
         call. = FALSE)
  }

  if (!is.logical(records$covered)) {
    stop("`records$covered` must be logical.", call. = FALSE)
  }
  if (anyNA(records$covered[d_interval])) {
    stop("Rows with `interval_computable == TRUE` must have non-NA logical `covered`.",
         call. = FALSE)
  }
  if (any(c("lower", "upper") %in% names(records))) {
    if (!all(c("lower", "upper") %in% names(records)) ||
        !is.numeric(records$lower) || !is.numeric(records$upper)) {
      stop("When interval endpoints are supplied, both `lower` and `upper` must be numeric columns.",
           call. = FALSE)
    }
    if (any(!is.finite(records$lower[d_interval]) |
            !is.finite(records$upper[d_interval]) |
            records$lower[d_interval] > records$upper[d_interval])) {
      stop("Rows with `interval_computable == TRUE` must have finite `lower <= upper` endpoints.",
           call. = FALSE)
    }
    expected_cov <- records$lower[d_interval] <= specification$truth &
      specification$truth <= records$upper[d_interval]
    if (!identical(records$covered[d_interval], expected_cov)) {
      stop("`records$covered` does not match `lower <= truth & truth <= upper` on `interval_computable == TRUE` rows.",
           call. = FALSE)
    }
  }

  if (!is.logical(records$band_covered)) {
    stop("`records$band_covered` must be logical.", call. = FALSE)
  }
  if (anyNA(records$band_covered[d_band])) {
    stop("Rows with `band_computable == TRUE` must have non-NA logical `band_covered`.",
         call. = FALSE)
  }
  if (any(c("band_lower", "band_upper") %in% names(records))) {
    if (!all(c("band_lower", "band_upper") %in% names(records)) ||
        !is.numeric(records$band_lower) || !is.numeric(records$band_upper)) {
      stop("When band endpoints are supplied, both `band_lower` and `band_upper` must be numeric columns.",
           call. = FALSE)
    }
    if (any(!is.finite(records$band_lower[d_band]) |
            !is.finite(records$band_upper[d_band]) |
            records$band_lower[d_band] > records$band_upper[d_band])) {
      stop("Rows with `band_computable == TRUE` must have finite `band_lower <= band_upper` endpoints.",
           call. = FALSE)
    }
    expected_band_cov <- records$band_lower[d_band] <= specification$truth &
      specification$truth <= records$band_upper[d_band]
    if (!identical(records$band_covered[d_band], expected_band_cov)) {
      stop("`records$band_covered` does not match `band_lower <= truth & truth <= band_upper` on `band_computable == TRUE` rows.",
           call. = FALSE)
    }
  }

  if (!is.logical(records$rejected)) {
    stop("`records$rejected` must be logical.", call. = FALSE)
  }
  if (anyNA(records$rejected[d_test])) {
    stop("Rows with `test_computable == TRUE` must have non-NA logical `rejected`.",
         call. = FALSE)
  }

  # Failure types, messages and warnings
  if (!is.character(records$failure_type)) {
    stop("`records$failure_type` must be a character vector.", call. = FALSE)
  }
  if (!is.character(records$failure_message)) {
    stop("`records$failure_message` must be a character vector.", call. = FALSE)
  }
  any_failed_event <- !(d_point & d_se & d_interval & d_band & d_test)
  no_failure_labels <- c("none", "complete", "")
  if (any(is.na(records$failure_type[any_failed_event]) |
          records$failure_type[any_failed_event] %in% no_failure_labels)) {
    stop("Rows with any false computability event must record a non-empty failure category in `records$failure_type`.",
         call. = FALSE)
  }
  all_events_ok <- !any_failed_event
  bad_ok_label <- !is.na(records$failure_type[all_events_ok]) &
    !(records$failure_type[all_events_ok] %in% c("none", "complete"))
  if (any(bad_ok_label)) {
    stop("Rows where all computability events are TRUE must not carry a failure category in `records$failure_type`.",
         call. = FALSE)
  }

  if (!mc_is_whole_number(records$n_warnings) || any(records$n_warnings < 0)) {
    stop("`records$n_warnings` must be non-negative whole numbers without NA.",
         call. = FALSE)
  }
  if (!is.character(records$warnings)) {
    stop("`records$warnings` must be a character vector.", call. = FALSE)
  }
  has_warn <- records$n_warnings > 0L
  if (any(is.na(records$warnings[has_warn]) | !nzchar(records$warnings[has_warn]))) {
    stop("Rows with `n_warnings > 0` must retain non-empty warning text in `records$warnings`.",
         call. = FALSE)
  }

  invisible(TRUE)
}

mc_evaluate_coverage_metric <- function(records, specification,
                                        metric_name = c("coverage", "band_coverage")) {
  metric_name <- match.arg(metric_name)
  n_att <- nrow(records)
  nominal <- specification$nominal_coverage
  floor_val <- nominal - 0.01

  ev_col <- if (metric_name == "coverage") {
    specification$events$interval
  } else {
    specification$events$band
  }
  cov_col <- if (metric_name == "coverage") "covered" else "band_covered"
  endpt_cols <- if (metric_name == "coverage") {
    c("lower", "upper")
  } else {
    c("band_lower", "band_upper")
  }

  d_ev <- records[[ev_col]]
  cov_vec <- records[[cov_col]]
  # All-attempt valid-and-covered indicator: false event rows contribute 0.
  valid_and_covered <- d_ev & (cov_vec %in% TRUE)
  k_all <- sum(valid_and_covered)
  m_comp <- sum(d_ev)
  n_fail <- n_att - m_comp

  prem <- mc_check_metric_premises(metric_name, specification)
  cp_all <- mc_clopper_pearson_ci(k_all, n_att, alpha_tail = 1 / 400)
  wilson_all <- mc_wilson_ci(k_all, n_att, conf_level = 0.99)
  cp_cond <- if (m_comp > 0L) {
    mc_clopper_pearson_ci(k_all, m_comp, alpha_tail = 1 / 400)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }

  status <- mc_classify_interval_decision(
    metric = metric_name,
    lower = cp_all[["lower"]],
    upper = cp_all[["upper"]],
    nominal = nominal,
    supported = prem$supported,
    unsupported_reasons = prem$unsupported_reasons
  )

  list(
    metric = metric_name,
    status = status,
    supported = prem$supported,
    premises_supported = prem$supported,
    arithmetic_defined = n_att > 0L,
    unsupported_reasons = prem$unsupported_reasons,
    estimate = k_all / n_att,
    k_valid_and_covered = k_all,
    m_computable = m_comp,
    n_attempted = n_att,
    n_failed = n_fail,
    computability_rate = m_comp / n_att,
    conditional_coverage_rate = if (m_comp > 0L) k_all / m_comp else NA_real_,
    conditional_cp_995_interval = cp_cond,
    interval = if (prem$supported) cp_all else c(lower = NA_real_, upper = NA_real_),
    descriptive_interval = cp_all,
    descriptive_wilson_99_interval = wilson_all,
    nominal = nominal,
    threshold = c(lower_strict = floor_val, upper = 1.0),
    strict_lower_bound_required = TRUE,
    denominator = paste0(n_att, "/", n_att),
    computable_denominator = paste0(m_comp, "/", n_att),
    conditioning_event = paste0("all_attempts_valid_and_covered (", ev_col, " AND ", cov_col, ")"),
    endpoints_verified = all(endpt_cols %in% names(records)),
    confidence_type = "exact_central_99.5pct_clopper_pearson_per_look_at_least_99pct_across_two_looks_per_metric"
  )
}

mc_evaluate_null_size_metric <- function(records, specification) {
  n_att <- nrow(records)
  nominal <- specification$nominal_size
  ceiling_val <- nominal + 0.01
  ev_col <- specification$events$test
  d_test <- records[[ev_col]]
  rej_vec <- records$rejected

  m_test <- sum(d_test)
  n_fail <- n_att - m_test
  valid_and_rejected <- d_test & (rej_vec %in% TRUE)
  k_rej <- sum(valid_and_rejected)

  # All-N joint rejection is always computed and reported separately,
  # and can NEVER confer D7 null-size PASS.
  joint_rate <- k_rej / n_att
  joint_cp_995 <- mc_clopper_pearson_ci(k_rej, n_att, alpha_tail = 1 / 400)
  joint_wilson_99 <- mc_wilson_ci(k_rej, n_att, conf_level = 0.99)

  prem <- mc_check_metric_premises("null_size", specification)
  reasons <- prem$unsupported_reasons
  if (m_test == 0L) {
    reasons <- unique(c(reasons, "zero_test_computability"))
  }
  supported <- prem$supported && (m_test > 0L)

  cond_cp_995 <- if (m_test > 0L) {
    mc_clopper_pearson_ci(k_rej, m_test, alpha_tail = 1 / 400)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  cond_wilson_99 <- if (m_test > 0L) {
    mc_wilson_ci(k_rej, m_test, conf_level = 0.99)
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }

  status <- mc_classify_interval_decision(
    metric = "null_size",
    lower = cond_cp_995[["lower"]],
    upper = cond_cp_995[["upper"]],
    nominal = nominal,
    supported = supported,
    unsupported_reasons = reasons
  )

  failure_breakdown <- table(records$failure_type[!d_test], useNA = "ifany")
  failure_counts <- stats::setNames(as.integer(failure_breakdown),
                                    names(failure_breakdown))

  list(
    metric = "null_size",
    status = status,
    supported = supported,
    premises_supported = prem$supported,
    arithmetic_defined = m_test > 0L,
    unsupported_reasons = reasons,
    estimate = if (m_test > 0L) k_rej / m_test else NA_real_,
    k_rejected_computable = k_rej,
    m_computable = m_test,
    n_attempted = n_att,
    n_failed = n_fail,
    computability_rate = m_test / n_att,
    failure_categories_noncomputable = failure_counts,
    joint_rejection_count = k_rej,
    joint_rejection_rate = joint_rate,
    joint_cp_995_interval = joint_cp_995,
    joint_wilson_99_interval = joint_wilson_99,
    joint_can_confer_pass = FALSE,
    interval = if (supported) cond_cp_995 else c(lower = NA_real_, upper = NA_real_),
    descriptive_interval = cond_cp_995,
    descriptive_wilson_99_interval = cond_wilson_99,
    nominal = nominal,
    threshold = c(lower = 0.0, upper_inclusive = ceiling_val),
    denominator = paste0(m_test, "/", n_att),
    conditioning_event = paste0("conditional_on_", ev_col, "_under_declared_null"),
    confidence_type = "exact_conditional_99.5pct_clopper_pearson_per_look_at_least_99pct_across_two_looks_per_metric"
  )
}

mc_nonfinite_error_metric <- function(metric, m, n_attempted, premises,
                                      n_zero_se = 0L) {
  # Finite stored estimates/truth can have a nonrepresentable binary64 difference.
  # This is derived-arithmetic INCONCLUSIVE, not an invalid stored input or fallback.
  out <- list(
    metric = metric, status = "INCONCLUSIVE", supported = FALSE,
    premises_supported = premises$supported, arithmetic_defined = FALSE,
    unsupported_reasons = unique(c(premises$unsupported_reasons,
                                  "nonfinite_error_arithmetic")),
    estimate = NA_real_, emp_sd = NA_real_, emp_var = NA_real_,
    psi = numeric(), psi_bar = NA_real_, if_var = NA_real_,
    z = stats::qnorm(399 / 400),
    interval = c(lower = NA_real_, upper = NA_real_),
    descriptive_interval = c(lower = NA_real_, upper = NA_real_),
    region_inclusive = c(lower = TRUE, upper = TRUE),
    m_computable = m, n_attempted = n_attempted, n_failed = n_attempted - m,
    computability_rate = if (n_attempted > 0L) m / n_attempted else NA_real_,
    denominator = paste0(m, "/", n_attempted),
    confidence_type = "asymptotic_99.5pct_per_look_delta_method_not_exact_finite_99pct"
  )
  if (metric == "bias") {
    out$raw_bias <- NA_real_; out$mcse <- NA_real_
    out$threshold <- c(lower = -0.10, upper = 0.10)
    out$conditioning_event <- "point_computable"
  } else {
    out$z <- 4.0
    out$log_estimate <- NA_real_; out$mean_se <- NA_real_
    out$n_zero_se <- n_zero_se; out$t_mean_joint <- NA_real_; out$mcse_log <- NA_real_
    out$threshold <- c(lower = 0.95, upper = 1.05)
    out$conditioning_event <- "point_computable AND se_computable"
    out$confidence_type <- "conservative_asymptotic_mc_interval_fixed_critical_4.0_required_calibration_target_99.5pct_limiting_2Phi4_minus_1_not_exact_99.5pct_no_exact_finite_M_guarantee"
    out$required_calibration_target <- 199 / 200
    out$limiting_normal_containment <- 2 * stats::pnorm(4.0) - 1
    out$exact_finite_m_guarantee <- FALSE
    out$method_metadata <- list(
      interval_class = "conservative_asymptotic_mc_interval",
      critical_value = 4.0,
      default_critical_value = 4.0,
      required_calibration_target = 199 / 200,
      limiting_normal_containment = 2 * stats::pnorm(4.0) - 1,
      limiting_containment_equals_target = FALSE,
      exact_finite_m_guarantee = FALSE,
      data_dependent_critical_selection = FALSE,
      looks_applied = c(2000L, 8000L)
    )
  }
  out
}

mc_evaluate_single_metric <- function(metric, records, specification) {
  if (metric %in% c("coverage", "band_coverage")) {
    return(mc_evaluate_coverage_metric(records, specification, metric))
  }
  if (metric == "null_size") {
    return(mc_evaluate_null_size_metric(records, specification))
  }
  n_att <- nrow(records)
  d_point <- records[[specification$events$point]]
  if (metric == "bias") {
    prem <- mc_check_metric_premises("bias", specification)
    t_err <- as.double(records$estimate[d_point]) - as.double(specification$truth)
    if (any(!is.finite(t_err))) {
      return(mc_nonfinite_error_metric("bias", sum(d_point), n_att, prem))
    }
    return(mc_continuous_bias_ci(
      t_errors = t_err,
      n_attempted = n_att,
      premises_supported = prem$supported,
      unsupported_reasons = prem$unsupported_reasons
    ))
  }
  if (metric == "se_ratio") {
    prem <- mc_check_metric_premises("se_ratio", specification)
    d_se <- records[[specification$events$se]]
    # Identical joint point/SE subset for both mean SE and empirical SD
    d_joint <- d_point & d_se
    t_err_joint <- as.double(records$estimate[d_joint]) - as.double(specification$truth)
    se_joint <- records$se[d_joint]
    if (any(!is.finite(t_err_joint))) {
      return(mc_nonfinite_error_metric("se_ratio", sum(d_joint), n_att, prem,
                                      n_zero_se = sum(se_joint == 0)))
    }
    return(mc_continuous_se_ratio_ci(
      t_errors = t_err_joint,
      se_values = se_joint,
      n_attempted = n_att,
      premises_supported = prem$supported,
      unsupported_reasons = prem$unsupported_reasons
    ))
  }
  stop(sprintf("Unsupported metric `%s`.", metric), call. = FALSE)
}

mc_compute_cell_status <- function(metric_statuses) {
  if (!is.character(metric_statuses) || length(metric_statuses) < 1L ||
      anyNA(metric_statuses) ||
      !all(metric_statuses %in% c("PASS", "FAIL", "INCONCLUSIVE"))) {
    stop("`metric_statuses` must be a non-empty vector of PASS/FAIL/INCONCLUSIVE.",
         call. = FALSE)
  }
  if (all(metric_statuses == "PASS")) {
    return("PASS")
  }
  if (any(metric_statuses == "FAIL")) {
    return("FAIL")
  }
  "INCONCLUSIVE"
}

mc_summarize_failures_and_events <- function(records, specification) {
  n_att <- nrow(records)
  ev <- specification$events
  d_point <- records[[ev$point]]
  d_se <- records[[ev$se]]
  d_joint <- d_point & d_se
  d_interval <- records[[ev$interval]]
  d_band <- records[[ev$band]]
  d_test <- records[[ev$test]]

  any_fail <- !(d_point & d_se & d_interval & d_band & d_test)
  ft_tab <- table(records$failure_type[any_fail], useNA = "ifany")

  list(
    n_attempted = n_att,
    computable_counts = c(
      point = sum(d_point),
      se = sum(d_se),
      joint_point_se = sum(d_joint),
      interval = sum(d_interval),
      band = sum(d_band),
      test = sum(d_test)
    ),
    failure_counts_by_event = c(
      point = n_att - sum(d_point),
      se = n_att - sum(d_se),
      joint_point_se = n_att - sum(d_joint),
      interval = n_att - sum(d_interval),
      band = n_att - sum(d_band),
      test = n_att - sum(d_test)
    ),
    denominators = c(
      point = paste0(sum(d_point), "/", n_att),
      se = paste0(sum(d_se), "/", n_att),
      joint_point_se = paste0(sum(d_joint), "/", n_att),
      interval = paste0(sum(d_interval), "/", n_att),
      band = paste0(sum(d_band), "/", n_att),
      test = paste0(sum(d_test), "/", n_att)
    ),
    n_rows_with_any_event_failure = sum(any_fail),
    failure_type_counts = stats::setNames(as.integer(ft_tab), names(ft_tab)),
    n_rows_with_warnings = sum(records$n_warnings > 0L),
    total_warning_count = sum(as.integer(records$n_warnings))
  )
}

mc_specifications_compatible <- function(spec_2000, spec_8000) {
  # Preserve the entire supplied specification, including opaque bindings.
  # Only the explicitly permitted seed-sequence extension differs.
  old_fixed <- spec_2000[setdiff(names(spec_2000), "seeds")]
  new_fixed <- spec_8000[setdiff(names(spec_8000), "seeds")]
  if (!identical(old_fixed, new_fixed)) {
    return(FALSE)
  }
  seeds_2000 <- as.integer(spec_2000$seeds)
  seeds_8000 <- as.integer(spec_8000$seeds)
  if (length(seeds_8000) != 8000L) {
    return(FALSE)
  }
  if (!identical(seeds_2000[seq_len(2000L)], seeds_8000[seq_len(2000L)])) {
    return(FALSE)
  }
  if (length(seeds_2000) == 8000L && !identical(seeds_2000, seeds_8000)) {
    return(FALSE)
  }
  TRUE
}

mc_records_prefix_identical <- function(prefix_df, full_df) {
  if (!is.data.frame(prefix_df) || !is.data.frame(full_df) ||
      nrow(prefix_df) != 2000L || nrow(full_df) < 2000L ||
      !identical(names(prefix_df), names(full_df))) {
    return(FALSE)
  }
  sub_df <- full_df[seq_len(2000L), , drop = FALSE]
  rownames(prefix_df) <- NULL
  rownames(sub_df) <- NULL
  identical(prefix_df, sub_df)
}

evaluate_private_mc_cell <- function(records, specification, look,
                                     prior_assessment = NULL) {
  validate_mc_records(records, specification, look)
  look <- as.integer(look)
  req_metrics <- specification$required_metrics

  if (look == 2000L) {
    if (!is.null(prior_assessment)) {
      stop("`prior_assessment` must be NULL at `look = 2000`.", call. = FALSE)
    }
    metrics_out <- stats::setNames(
      lapply(req_metrics, function(m) {
        res <- mc_evaluate_single_metric(m, records, specification)
        res$decision_look <- 2000L
        res$preserved_from_look_2000 <- FALSE
        res$reevaluated_at_8000 <- FALSE
        res
      }),
      req_metrics
    )
    metric_statuses <- vapply(metrics_out, function(x) x$status, character(1L))
    cell_status <- mc_compute_cell_status(metric_statuses)
    prefix_copy <- records
    rownames(prefix_copy) <- NULL

    return(structure(
      list(
        cell_id = specification$cell_id,
        target = specification$target,
        truth = specification$truth,
        look = 2000L,
        n_attempted = 2000L,
        required_metrics = req_metrics,
        metric_statuses = metric_statuses,
        cell_status = cell_status,
        can_extend_to_8000 = identical(cell_status, "INCONCLUSIVE"),
        unresolved_at_cap = FALSE,
        metrics = metrics_out,
        failure_and_event_summary = mc_summarize_failures_and_events(records, specification),
        specification = specification,
        records_prefix = prefix_copy,
        history = list(
          look_2000 = NULL,
          prefix_verified = TRUE,
          external_provenance_premise = "Stateless evaluator cannot detect unreported external peeks/redraws; prespecification and no-hidden-history are external premises."
        ),
        certification_scope = "Conditional on declared computability events and mathematical/calibration premises; per-metric confidence only, never simultaneous across metrics/targets/cells and not a fitted-estimator certification."
      ),
      class = "mc_cell_assessment"
    ))
  }

  # look == 8000L
  if (length(specification$seeds) != 8000L) {
    stop("At `look = 8000`, `specification$seeds` must have length 8000.",
         call. = FALSE)
  }
  if (is.null(prior_assessment) || !inherits(prior_assessment, "mc_cell_assessment") ||
      !is.list(prior_assessment) || !identical(prior_assessment$look, 2000L)) {
    stop("At `look = 8000`, `prior_assessment` must be a valid `look = 2000` `mc_cell_assessment` object.",
         call. = FALSE)
  }
  if (!mc_specifications_compatible(prior_assessment$specification, specification)) {
    stop("At `look = 8000`, `specification` must be identical to `prior_assessment$specification` with the exact same first 2000 seed prefix.",
         call. = FALSE)
  }
  if (!mc_records_prefix_identical(prior_assessment$records_prefix, records)) {
    stop("At `look = 8000`, the first 2000 rows of `records` must be byte-for-byte identical to `prior_assessment$records_prefix`.",
         call. = FALSE)
  }

  # Deterministically reproduce look 2000 on the first 2000 rows rather than trusting a decorative status.
  reproduced_2000 <- evaluate_private_mc_cell(
    records = records[seq_len(2000L), , drop = FALSE],
    specification = prior_assessment$specification,
    look = 2000L,
    prior_assessment = NULL
  )
  if (!identical(prior_assessment$cell_status, reproduced_2000$cell_status) ||
      !identical(prior_assessment$metric_statuses, reproduced_2000$metric_statuses)) {
    stop("Supplied `prior_assessment` status does not match the deterministic reproduction of the first 2000 rows.",
         call. = FALSE)
  }
  for (m_nm in req_metrics) {
    p_m <- prior_assessment$metrics[[m_nm]]
    r_m <- reproduced_2000$metrics[[m_nm]]
    if (is.null(p_m) || !identical(p_m$status, r_m$status) ||
        !identical(p_m$interval, r_m$interval) ||
        !identical(p_m$denominator, r_m$denominator)) {
      stop(sprintf("Supplied `prior_assessment` metric `%s` does not match deterministic look-2000 reproduction.",
                   m_nm), call. = FALSE)
    }
  }

  if (!identical(reproduced_2000$cell_status, "INCONCLUSIVE")) {
    stop(sprintf("Extension to 8000 is permitted only from a cell-INCONCLUSIVE assessment at 2000; prior cell status was `%s`.",
                 reproduced_2000$cell_status), call. = FALSE)
  }

  # Preserve conclusive first-look metric decisions; reevaluate only previously inconclusive metrics.
  metrics_out <- stats::setNames(
    lapply(req_metrics, function(m) {
      first_look_m <- reproduced_2000$metrics[[m]]
      full_8000_m <- mc_evaluate_single_metric(m, records, specification)
      if (first_look_m$status %in% c("PASS", "FAIL")) {
        preserved <- first_look_m
        preserved$decision_look <- 2000L
        preserved$preserved_from_look_2000 <- TRUE
        preserved$reevaluated_at_8000 <- FALSE
        preserved$look_8000_descriptive_only <- full_8000_m
        return(preserved)
      }
      full_8000_m$decision_look <- 8000L
      full_8000_m$preserved_from_look_2000 <- FALSE
      full_8000_m$reevaluated_at_8000 <- TRUE
      full_8000_m$look_2000_prior <- first_look_m
      full_8000_m
    }),
    req_metrics
  )

  metric_statuses <- vapply(metrics_out, function(x) x$status, character(1L))
  cell_status <- mc_compute_cell_status(metric_statuses)

  structure(
    list(
      cell_id = specification$cell_id,
      target = specification$target,
      truth = specification$truth,
      look = 8000L,
      n_attempted = 8000L,
      required_metrics = req_metrics,
      metric_statuses = metric_statuses,
      cell_status = cell_status,
      can_extend_to_8000 = FALSE,
      unresolved_at_cap = identical(cell_status, "INCONCLUSIVE"),
      metrics = metrics_out,
      failure_and_event_summary = mc_summarize_failures_and_events(records, specification),
      specification = specification,
      records_prefix = reproduced_2000$records_prefix,
      history = list(
        look_2000 = reproduced_2000,
        prefix_verified = TRUE,
        preserved_metrics = names(Filter(function(x) isTRUE(x$preserved_from_look_2000), metrics_out)),
        reevaluated_metrics = names(Filter(function(x) isTRUE(x$reevaluated_at_8000), metrics_out)),
        external_provenance_premise = "Stateless evaluator cannot detect unreported external peeks/redraws; prespecification and no-hidden-history are external premises."
      ),
      certification_scope = "Conditional on declared computability events and mathematical/calibration premises; per-metric confidence only, never simultaneous across metrics/targets/cells and not a fitted-estimator certification."
    ),
    class = "mc_cell_assessment"
  )
}

mc_calibration_outer_seed <- function(cell_index, outer_index) {
  if (length(cell_index) != 1L || !mc_is_whole_number(cell_index) ||
      !(as.integer(cell_index) %in% 1L:4L)) {
    stop("`cell_index` must be an integer in 1..4.", call. = FALSE)
  }
  if (length(outer_index) != 1L || !mc_is_whole_number(outer_index) ||
      !(as.integer(outer_index) %in% 1L:50L)) {
    stop("`outer_index` must be an integer in 1..50.", call. = FALSE)
  }
  as.integer(14322000L + 10000L * as.integer(cell_index) + as.integer(outer_index))
}

mc_calibration_continuation_decision <- function(assessment_2000) {
  if (!inherits(assessment_2000, "mc_cell_assessment") ||
      !identical(assessment_2000$look, 2000L)) {
    stop("`assessment_2000` must be a `look = 2000` `mc_cell_assessment` object.",
         call. = FALSE)
  }
  st <- assessment_2000$cell_status
  m_st <- assessment_2000$metric_statuses
  continue_flag <- identical(st, "INCONCLUSIVE")
  stop_reason <- if (continue_flag) {
    "cell_inconclusive_extend_to_8000"
  } else if (identical(st, "PASS")) {
    "cell_pass_stop_at_2000"
  } else if (any(m_st == "FAIL") && any(m_st == "INCONCLUSIVE")) {
    "cell_fail_plus_inconclusive_stop_at_2000_no_extension"
  } else {
    "cell_fail_stop_at_2000"
  }
  list(
    cell_status_2000 = st,
    continue_to_8000 = continue_flag,
    frozen_conclusive_metrics = names(m_st[m_st %in% c("PASS", "FAIL")]),
    metrics_to_reevaluate = if (continue_flag) names(m_st[m_st == "INCONCLUSIVE"]) else character(),
    stop_reason = stop_reason
  )
}

prepare_mc_calibration_plumbing <- function() {
  cells <- list(
    list(
      index = 1L,
      law = "Z~N(0,1)",
      computability = 1,
      beta = 1 / 10,
      rho = 1,
      EZ3 = 0,
      EZ4 = 3,
      Var_psi_beta = 201 / 200,
      Var_psi_logrho = 9 / 32,
      coverage_conditional = 19 / 20,
      coverage_all_attempt = 19 / 20,
      size_conditional = 1 / 20,
      joint_rejection = 1 / 20,
      role = "Light-tailed Gaussian calibration law"
    ),
    list(
      index = 2L,
      law = "Z=sqrt(3/5)*t_5",
      computability = 1,
      beta = 1 / 10,
      rho = 1,
      EZ3 = 0,
      EZ4 = 9,
      Var_psi_beta = 51 / 50,
      Var_psi_logrho = 9 / 8,
      coverage_conditional = 19 / 20,
      coverage_all_attempt = 19 / 20,
      size_conditional = 1 / 20,
      joint_rejection = 1 / 20,
      role = "Heavy-tailed finite-fourth-moment standardized t5 calibration law (t3 has infinite fourth moment and is a negative control)"
    ),
    list(
      index = 3L,
      law = "Z=Exp(rate=1)-1",
      computability = 1,
      beta = 1 / 10,
      rho = 1,
      EZ3 = 2,
      EZ4 = 9,
      Var_psi_beta = 41 / 50,
      Var_psi_logrho = 9 / 8,
      coverage_conditional = 19 / 20,
      coverage_all_attempt = 19 / 20,
      size_conditional = 1 / 20,
      joint_rejection = 1 / 20,
      role = "Skewed centered exponential calibration law"
    ),
    list(
      index = 4L,
      law = "Z~N(0,1); D=1{Z>0}; T=Z-sqrt(2/pi), S=sqrt(1-2/pi); point/SE/test/interval computability all use D",
      computability = 1 / 2,
      beta = 0,
      rho = 1,
      mean_Z_given_D = sqrt(2 / pi),
      variance_Z_given_D = 1 - 2 / pi,
      central_third = sqrt(2 / pi) * (4 / pi - 1),
      central_fourth = 3 - 4 / pi - 12 / (pi * pi),
      coverage_conditional = 19 / 20,
      coverage_all_attempt = 19 / 40,
      size_conditional = 1 / 20,
      joint_rejection = 1 / 40,
      role = "Predeclared dependent-computability diagnostic; low all-attempt coverage remains a failure of that coverage criterion, never relabeled after results."
    )
  )

  grid_rows <- expand.grid(
    outer_index = seq_len(50L),
    cell_index = seq_len(4L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid_rows <- grid_rows[, c("cell_index", "outer_index")]
  grid_rows$outer_seed <- 14322000L + 10000L * grid_rows$cell_index + grid_rows$outer_index
  grid_rows$initial_look <- 2000L
  grid_rows$cap_look <- 8000L
  grid_rows$batch_rows <- 1000L
  grid_rows$executed <- FALSE

  list(
    status = "PROPOSED_PLUMBING_ONLY_UNLAUNCHED",
    historical_plan_sha256 = "d7f055c8d201146feb2c55900f16d7f0605e063f09869080a1687d4a8a44b683",
    historical_plan_modified = FALSE,
    controlling_authority = "Approved D14/D15/D16 at commit740718e and root clarification: explicit cell status; coverage PASS strictly above its floor",
    execution_now = FALSE,
    outer_draws_executed = 0L,
    inner_attempts_executed = 0L,
    fits_executed = 0L,
    outer_draws_per_cell = 50L,
    cell_count = 4L,
    outer_draws_total = 200L,
    outer_seed_formula = "14322000 + 10000*cell_index + outer_index; cell_index=1..4, outer_index=1..50",
    inner_attempt_assessments = c(2000L, 8000L),
    batch_rows_maximum = 1000L,
    chunks_initial_2000 = 2L,
    chunks_extension_8000 = 6L,
    maximum_inner_attempts_upper_bound = 1600000L,
    maximum_primitive_scalar_draws_Z_plus_U = 3200000L,
    maximum_outer_look_rows = 400L,
    per_look_binomial_confidence = "199/200 central Clopper-Pearson; each tail=1/400",
    per_look_continuous_intervals = "beta_hat +/- qnorm(399/400)*MCSE_beta; exp(log(rho_hat) +/- 4.0*MCSE_logrho) [conservative asymptotic MC interval, critical 4.0, required calibration target 99.5%, limiting normal containment 2*Phi(4)-1 != 99.5%, no exact finite-M guarantee]",
    se_ratio_interval_metadata = list(
      interval_class = "conservative_asymptotic_mc_interval",
      critical_value = 4.0,
      required_calibration_target = 199 / 200,
      limiting_normal_containment = 2 * stats::pnorm(4.0) - 1,
      limiting_containment_equals_target = FALSE,
      exact_finite_m_guarantee = FALSE,
      data_dependent_critical_selection = FALSE,
      looks_applied = c(2000L, 8000L)
    ),
    continuation_rule = "Only previous cell status INCONCLUSIVE permits 8000; all metrics PASS means cell PASS, any metric FAIL means cell FAIL. Preserve first-look conclusive metric decisions and reevaluate previously inconclusive metrics only. FAIL+INCONCLUSIVE must not extend.",
    precision_limit = "At 99.5% nominal, only 0.25 expected marginal misses per law with 50 outer draws. Counts/interval-hit uncertainty are plumbing diagnostics, not validation of nominal tail calibration or estimator accuracy.",
    expected_marginal_misses_per_law_at_995 = 0.25,
    can_certify_995_calibration = FALSE,
    can_certify_fitted_estimator = FALSE,
    bootstrap_precision_limits_adopted = FALSE,
    failure_ceiling_adopted = FALSE,
    cells = cells,
    outer_seed_grid = grid_rows
  )
}
