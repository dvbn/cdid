# Independent reference calculators and deterministic fixture builders for
# test-mc-validation.R. This file never sources or defines R/mc_validation.R.

find_project_file <- function(rel_path) {
  candidates <- c(
    rel_path,
    file.path("..", "..", rel_path),
    file.path("..", rel_path),
    file.path("fixtures", "mc_validation", rel_path),
    file.path("tests", "testthat", "fixtures", "mc_validation", rel_path)
  )
  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  NA_character_
}

ref_default_events <- function() {
  list(
    point = "point_computable",
    se = "se_computable",
    interval = "interval_computable",
    band = "band_computable",
    test = "test_computable"
  )
}

ref_default_premises <- function(overrides = list()) {
  base_p <- list(
    iid_attempts = TRUE,
    fixed_predeclared_events = TRUE,
    no_hidden_peeks = TRUE,
    null_holds = TRUE,
    positive_point_event_probability = TRUE,
    positive_joint_event_probability = TRUE,
    positive_test_event_probability = TRUE,
    positive_conditional_sd_point = TRUE,
    positive_conditional_sd_joint = TRUE,
    conditional_fourth_moment_T_point = TRUE,
    conditional_fourth_moment_T_joint = TRUE,
    conditional_second_moment_S_joint = TRUE,
    positive_conditional_mean_se_joint = TRUE,
    nondegenerate_if_variance_bias = TRUE,
    nondegenerate_if_variance_se_ratio = TRUE,
    calibration_supported_bias = TRUE,
    calibration_supported_se_ratio = TRUE
  )
  for (nm in names(overrides)) {
    base_p[[nm]] <- overrides[[nm]]
  }
  base_p
}

make_mc_spec <- function(n_seeds = 8000L,
                         first_seed = 700001L,
                         cell_id = "cell_baseline_01",
                         target = "att_g2_t2",
                         truth = 0.0,
                         nominal_coverage = 0.95,
                         nominal_size = 0.05,
                         null_premise = TRUE,
                         required_metrics = c("coverage", "band_coverage",
                                              "null_size", "bias", "se_ratio"),
                         events = ref_default_events(),
                         premises = ref_default_premises(),
                         provenance = list(source = "deterministic_unit_fixture")) {
  list(
    cell_id = cell_id,
    target = target,
    truth = truth,
    seeds = as.integer(first_seed + seq_len(n_seeds) - 1L),
    nominal_coverage = nominal_coverage,
    nominal_size = nominal_size,
    null_premise = null_premise,
    required_metrics = required_metrics,
    events = events,
    premises = premises,
    provenance = provenance
  )
}

# Deterministic cyclic error/SE pattern generator (zero RNG usage)
make_mc_records <- function(spec, look = 2000L,
                            bias_shift = 0.0,
                            se_scale = 1.0,
                            n_covered = as.integer(round(0.965 * look)),
                            n_band_covered = as.integer(round(0.965 * look)),
                            n_rejected = as.integer(round(0.035 * look))) {
  look <- as.integer(look)
  seeds <- as.integer(spec$seeds[seq_len(look)])
  # Deterministic 3-point standardized cycle (-1, 0, 1) with mean 0 and variance 1:
  # Repeat block of 4: (-1, 0, 1, 0) * sqrt(2) has mean 0 and population variance 1;
  # vary SE slightly so plug-in IF is non-degenerate.
  idx <- seq_len(look)
  cycle4 <- c(-1, 0, 1, 0) * sqrt(2)
  z_det <- cycle4[((idx - 1L) %% 4L) + 1L]
  # Add tiny deterministic modulation across 8-period so sample SD and IF variance are strictly positive
  mod8 <- c(-0.2, 0.1, 0.2, -0.1, 0.15, -0.05, -0.15, 0.05)[((idx - 1L) %% 8L) + 1L]
  raw_t <- z_det + mod8
  raw_t <- (raw_t - mean(raw_t)) / stats::sd(raw_t)
  t_err <- bias_shift + raw_t

  se_mod <- c(0.92, 1.08, 0.95, 1.05, 0.90, 1.10, 0.97, 1.03)[((idx - 1L) %% 8L) + 1L]
  se_vec <- se_scale * (se_mod / mean(se_mod))

  est_vec <- spec$truth + t_err
  cov_vec <- idx <= as.integer(n_covered)
  band_cov_vec <- idx <= as.integer(n_band_covered)
  rej_vec <- idx <= as.integer(n_rejected)

  data.frame(
    attempt = idx,
    seed = seeds,
    cell_id = rep(spec$cell_id, look),
    target = rep(spec$target, look),
    truth = rep(spec$truth, look),
    estimate = est_vec,
    se = se_vec,
    covered = cov_vec,
    band_covered = band_cov_vec,
    rejected = rej_vec,
    point_computable = rep(TRUE, look),
    se_computable = rep(TRUE, look),
    interval_computable = rep(TRUE, look),
    band_computable = rep(TRUE, look),
    test_computable = rep(TRUE, look),
    failure_type = rep("none", look),
    failure_message = rep(NA_character_, look),
    n_warnings = rep(0L, look),
    warnings = rep(NA_character_, look),
    stringsAsFactors = FALSE
  )
}

ref_binom_cdf_recurrence <- function(n, k, p) {
  mass <- (1 - p)^n
  out <- mass
  if (k >= 1L) {
    for (j in seq_len(k) - 1L) {
      mass <- mass * ((n - j) / (j + 1)) * (p / (1 - p))
      out <- out + mass
    }
  }
  out
}

ref_sample_moments <- function(t_vec, s_vec) {
  m <- length(t_vec)
  stopifnot(m >= 2L, length(s_vec) == m)
  mu <- sum(t_vec) / m
  nu <- sum(s_vec) / m
  v <- sum((t_vec - mu)^2) / (m - 1L)
  stopifnot(v > 0, nu > 0)
  sd_val <- sqrt(v)
  beta_hat <- mu / sd_val
  rho_hat <- nu / sd_val
  pb <- (t_vec - mu) / sd_val - beta_hat * ((t_vec - mu)^2 - v) / (2 * v)
  pr <- (s_vec - nu) / nu - ((t_vec - mu)^2 - v) / (2 * v)
  pb_bar <- sum(pb) / m
  pr_bar <- sum(pr) / m
  vb <- sum((pb - pb_bar)^2) / (m - 1L)
  vr <- sum((pr - pr_bar)^2) / (m - 1L)
  stopifnot(vb > 0, vr > 0)
  z_b <- stats::qnorm(399 / 400)
  z_r <- 4.0
  mcse_b <- sqrt(vb / m)
  mcse_r <- sqrt(vr / m)
  list(
    mu = mu,
    v = v,
    sd = sd_val,
    nu = nu,
    beta_hat = beta_hat,
    rho_hat = rho_hat,
    pb = pb,
    pr = pr,
    pb_bar = pb_bar,
    pr_bar = pr_bar,
    vb = vb,
    vr = vr,
    mcse2_b = vb / m,
    mcse2_r = vr / m,
    mcse_b = mcse_b,
    mcse_r = mcse_r,
    z_b = z_b,
    z_r = z_r,
    ci_b = c(lower = beta_hat - z_b * mcse_b, upper = beta_hat + z_b * mcse_b),
    ci_r = c(lower = exp(log(rho_hat) - z_r * mcse_r),
             upper = exp(log(rho_hat) + z_r * mcse_r)),
    ci_r_legacy_qnorm = c(lower = exp(log(rho_hat) - z_b * mcse_r),
                          upper = exp(log(rho_hat) + z_b * mcse_r))
  )
}
