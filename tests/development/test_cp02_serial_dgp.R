# Deterministic verification suite for candidate/inst/development/cp02_serial_dgp.R
#
# Strict validation-preparation constraints:
#   - NO RNG initialization or draws (.Random.seed must never exist).
#   - NO stochastic DGP constructor calls (sim_draw, cp02_baseline_variant_draw,
#     cp02_serial_dgp_draw with live draws).
#   - NO estimator fits, package inference, bootstrap, calibration, or MC runs.
#   - Uses only hand-built fake accepted rotating draw objects, supplied
#     deterministic innovations, tiny enumeration, and independent finite
#     cohort/X population arithmetic.

assert_no_rng_touched <- function(context_label) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stop("RNG invariant violated during ", context_label,
         ": .Random.seed was initialized in globalenv().", call. = FALSE)
  }
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) {
    stop("Assertion failed: ", msg, call. = FALSE)
  }
}

assert_error <- function(expr, pattern = NULL, label = "expected error") {
  err <- tryCatch({
    force(expr)
    NULL
  }, error = function(e) e)
  if (is.null(err)) {
    stop("Expected an error for [", label, "], but expression succeeded.",
         call. = FALSE)
  }
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(err), ignore.case = TRUE)) {
    stop("Error message for [", label, "] ('", conditionMessage(err),
         "') did not match expected pattern '", pattern, "'.", call. = FALSE)
  }
  invisible(conditionMessage(err))
}

assert_no_rng_touched("startup")

source("inst/development/seed.R")
source("inst/development/dgp.R")
source("inst/development/cp02_baseline_variants_dgp.R")
source("inst/development/cp02_serial_dgp.R")

assert_no_rng_touched("sourcing reference and candidate files")

# ---------------------------------------------------------------------------
# 1. Catalog, specification metadata, original sim_spec preservation, and
#    strict refusal of undeclared rho / law / design overrides
# ---------------------------------------------------------------------------
test_catalog_spec_and_refusals <- function() {
  cat_df <- cp02_serial_dgp_catalog()
  assert_true(is.data.frame(cat_df) && nrow(cat_df) == 2L,
              "catalog has exactly 2 rows")
  assert_true(identical(cat_df$law, c("AR05", "AR1")),
              "catalog laws are AR05 and AR1")
  assert_true(identical(cat_df$design_table_law_index, c(20L, 21L)),
              "catalog design table law indices are 20 and 21")
  assert_true(identical(cat_df$rho, c(0.5, 1.0)),
              "catalog rho values are 0.5 and 1.0")
  assert_true(identical(cat_df$innovation_variance, c(1.0, 1.0)),
              "catalog innovation_variance is 1.0")
  assert_true(identical(cat_df$stationary, c(TRUE, FALSE)),
              "catalog stationarity flags are TRUE for AR05 and FALSE for AR1")
  assert_true(identical(cat_df$pattern_family, c("rotating", "rotating")),
              "catalog pattern_family is rotating")
  assert_true(identical(cat_df$x_trend, c(0.0, 0.0)),
              "catalog x_trend is 0.0")

  # Original sim_spec() default x_trend = 0.5 must remain unaltered.
  orig_spec <- sim_spec()
  assert_true(identical(orig_spec$x_trend, 0.5),
              "original sim_spec() default x_trend == 0.5 is unaltered")

  s05 <- cp02_serial_dgp_spec("AR05")
  s1 <- cp02_serial_dgp_spec("AR1")
  assert_true(identical(s05, cp02_serial_dgp_spec(law = "AR05", rho = 0.5)),
              "AR05 spec matches explicit (law = 'AR05', rho = 0.5)")
  assert_true(identical(s05, cp02_serial_dgp_spec(rho = 0.5)),
              "AR05 spec matches rho = 0.5 lookup")
  assert_true(identical(s1, cp02_serial_dgp_spec(law = "AR1", rho = 1.0)),
              "AR1 spec matches explicit (law = 'AR1', rho = 1.0)")
  assert_true(identical(s1, cp02_serial_dgp_spec(rho = 1.0)),
              "AR1 spec matches rho = 1.0 lookup")

  assert_true(identical(s05$target_spec, sim_spec(x_trend = 0)),
              "AR05 target_spec equals sim_spec(x_trend = 0)")
  assert_true(identical(s1$target_spec, sim_spec(x_trend = 0)),
              "AR1 target_spec equals sim_spec(x_trend = 0)")
  assert_true(identical(s05$base_variant_spec$variant, "rotating") &&
                identical(s1$base_variant_spec$variant, "rotating"),
              "Both laws bind to rotating base variant")
  assert_true(identical(s05$campaign_fit_spec$control_group, "NT") &&
                identical(s05$campaign_fit_spec$xformla, "~1") &&
                identical(s05$campaign_fit_spec$weighting, "T"),
              "Campaign metadata records NT / ~1 / T")

  # Clear refusal of undeclared laws, rho values, and design overrides.
  assert_error(cp02_serial_dgp_spec(), "Must specify", "empty spec call")
  for (bad_law in list("AR0", "AR025", "AR09", "AR2", "AR", "A", "B", "", NA_character_, c("AR05", "AR1"), 5)) {
    assert_error(cp02_serial_dgp_spec(law = bad_law), "Undeclared",
                 paste("bad law:", paste(bad_law, collapse = ",")))
  }
  for (bad_rho in list(0, -0.5, 0.25, 0.49, 0.51, 0.75, 0.99, 1.01, -1, 2, NA_real_, NaN, Inf, -Inf, c(0.5, 1), "0.5")) {
    assert_error(cp02_serial_dgp_spec(rho = bad_rho), "Undeclared",
                 paste("bad rho:", paste(bad_rho, collapse = ",")))
  }
  assert_error(cp02_serial_dgp_spec(law = "AR05", rho = 1.0), "Inconsistent",
               "mismatched AR05 and rho = 1")
  assert_error(cp02_serial_dgp_spec(law = "AR1", rho = 0.5), "Inconsistent",
               "mismatched AR1 and rho = 0.5")
  assert_error(cp02_serial_dgp_spec("AR05", x_trend = 0.5), "Undeclared design override",
               "extra x_trend override")
  assert_error(cp02_serial_dgp_spec("AR1", pattern_family = "balanced"), "Undeclared design override",
               "extra pattern_family override")
  assert_error(cp02_serial_dgp_spec("AR05", sd_eps = 2), "Undeclared design override",
               "extra sd_eps override")

  # Refusal of tampered spec objects.
  bad_s1 <- s1
  bad_s1$stationary <- TRUE
  assert_error(cp02_validate_serial_spec(bad_s1), "Invalid start/stationarity",
               "AR1 marked stationary")

  bad_s05_x <- s05
  bad_s05_x$target_spec$x_trend <- 0.5
  assert_error(cp02_validate_serial_spec(bad_s05_x), "Undeclared `target_spec`",
               "tampered x_trend = 0.5")

  bad_s05_pat <- s05
  bad_s05_pat$base_variant_spec <- cp02_baseline_variant_spec("balanced", base_spec = sim_spec(x_trend = 0))
  assert_error(cp02_validate_serial_spec(bad_s05_pat), "Undeclared base variant",
               "tampered balanced variant")

  bad_s05_var <- s05
  bad_s05_var$innovation_variance <- 2.0
  assert_error(cp02_validate_serial_spec(bad_s05_var), "Undeclared innovation",
               "tampered innovation variance")

  # Fixed metadata must not redirect seeds or mislabel the transformed law.
  bad_s05_index <- s05
  bad_s05_index$design_table_law_index <- 21L
  assert_error(cp02_serial_dgp_seed(bad_s05_index, "development"),
               "Undeclared serial DGP metadata", "AR05 in AR1 seed namespace")
  for (field in c("status", "start_distribution", "recurrence_equation",
                  "eps_scale_semantics", "legacy_provenance", "patterns",
                  "campaign_fit_spec", "draw_order", "stream_note")) {
    bad_metadata <- s05
    bad_metadata[[field]] <- NULL
    assert_error(cp02_validate_serial_spec(bad_metadata),
                 "Undeclared serial DGP metadata",
                 paste("missing declared metadata", field))
  }
  bad_s05_campaign <- s05
  bad_s05_campaign$campaign_fit_spec$xformla <- "~x"
  assert_error(cp02_validate_serial_spec(bad_s05_campaign),
               "Undeclared serial DGP metadata", "tampered no-covariate fit")
  bad_s05_order <- s05
  bad_s05_order$draw_order <- "second generator"
  assert_error(cp02_validate_serial_spec(bad_s05_order),
               "Undeclared serial DGP metadata", "tampered draw order")

  # Refusal in cp02_serial_dgp_draw before any constructor/RNG call.
  assert_error(cp02_serial_dgp_draw("AR_BAD", n = 10, seed = 1), "Undeclared",
               "draw with undeclared law")
  assert_error(cp02_serial_dgp_draw(rho = 0.25, n = 10, seed = 1), "Undeclared",
               "draw with undeclared rho")
  assert_error(cp02_serial_dgp_draw(0.25, n = 10, seed = 1), "Undeclared",
               "draw with positional undeclared rho")
  assert_error(cp02_serial_dgp_draw("AR05", n = 10, seed = 1, x_trend = 0.5), "Undeclared design override",
               "draw with extra design override")
  assert_error(cp02_serial_dgp_draw("AR05", n = 0, seed = 1), "whole number >= 1",
               "draw with n = 0")
  assert_error(cp02_serial_dgp_draw("AR05", n = 10, seed = NA_real_), "seed",
               "draw with NA seed")
  assert_error(cp02_serial_dgp_draw("AR1", n = 10, seed = 1, n_clusters = 0), "n_clusters",
               "draw with n_clusters = 0")

  # Refusal across all other entry points for undeclared rho or design overrides.
  assert_error(cp02_serial_loading_matrix(rho = 0.3), "Undeclared", "loading_matrix undeclared rho")
  assert_error(cp02_serial_covariance_matrix(rho = 0.3), "Undeclared", "covariance_matrix undeclared rho")
  assert_error(cp02_serial_noise_moments(rho = 0.3), "Undeclared", "noise_moments undeclared rho")
  assert_error(cp02_serial_dgp_truth(rho = 0.3), "Undeclared", "dgp_truth undeclared rho")
  assert_error(cp02_serial_transform_innovations(rep(0, 5), rho = 0.3), "Undeclared", "transform_innovations undeclared rho")

  # Deterministic seed formula check (Section 5 of CP02_FINAL_MC_DESIGN.md:
  # AR05 = law 20, AR1 = law 21).
  assert_true(identical(cp02_serial_dgp_seed("AR05", "development", rung = 1L, rep = 1L),
                        620100001L),
              "AR05 development rung 1 rep 1 seed == 620100001")
  assert_true(identical(cp02_serial_dgp_seed("AR1", "development", rung = 2L, rep = 200L),
                        621200200L),
              "AR1 development rung 2 rep 200 seed == 621200200")
  assert_true(identical(cp02_serial_dgp_seed("AR05", "production", rung = 3L, rep = 2000L),
                        720302000L),
              "AR05 production rung 3 rep 2000 seed == 720302000")
  assert_true(identical(cp02_serial_dgp_seed("AR1", "inference_production", rung = 3L, rep = 8000L),
                        871308000L),
              "AR1 inference_production rung 3 rep 8000 seed == 871308000")
  assert_error(cp02_serial_dgp_seed("AR05", "development", rung = 4L, rep = 1L),
               "rung", "seed with invalid rung 4")
  assert_error(cp02_serial_dgp_seed("AR05", "development", rung = 1L, rep = 201L),
               "rep", "development seed with rep > 200")

  assert_no_rng_touched("test_catalog_spec_and_refusals")
}

# ---------------------------------------------------------------------------
# 2. Recurrence vs independently formed lower-triangular loading matrix and
#    analytic finite-horizon covariance (rho = 0.5 stationary, rho = 1 RW)
# ---------------------------------------------------------------------------
test_recurrence_loading_and_covariance <- function() {
  periods <- 1:5
  T_len <- 5L

  # Independently construct 5x5 lower-triangular loading matrix for AR05.
  rho05 <- 0.5
  s1_05 <- 1 / sqrt(1 - rho05^2)
  L05_ind <- matrix(0, nrow = T_len, ncol = T_len,
                    dimnames = list(as.character(periods), as.character(periods)))
  for (r in 1:T_len) {
    L05_ind[r, 1L] <- (rho05^(r - 1L)) * s1_05
    if (r >= 2L) {
      for (c in 2:r) {
        L05_ind[r, c] <- rho05^(r - c)
      }
    }
  }
  # Independently construct 5x5 stationary covariance rho^|s-t| / (1 - rho^2).
  Sigma05_ind <- outer(periods, periods, function(s, t) {
    (rho05^abs(s - t)) / (1 - rho05^2)
  })
  dimnames(Sigma05_ind) <- list(as.character(periods), as.character(periods))

  # Independently construct 5x5 lower-triangular loading matrix for AR1 (unit lower triangular).
  L1_ind <- matrix(0, nrow = T_len, ncol = T_len,
                   dimnames = list(as.character(periods), as.character(periods)))
  for (r in 1:T_len) {
    L1_ind[r, 1:r] <- 1.0
  }
  # Independently construct 5x5 random-walk covariance min(s, t).
  Sigma1_ind <- outer(periods, periods, pmin) * 1.0
  dimnames(Sigma1_ind) <- list(as.character(periods), as.character(periods))

  # Verify L %*% t(L) == Sigma independently.
  assert_true(max(abs(L05_ind %*% t(L05_ind) - Sigma05_ind)) < 1e-15,
              "AR05 L_ind %*% t(L_ind) equals rho^|s-t|/(1-rho^2)")
  assert_true(max(abs(L1_ind %*% t(L1_ind) - Sigma1_ind)) == 0,
              "AR1 L_ind %*% t(L_ind) equals min(s, t) exactly")

  # Verify candidate functions match independent matrices.
  assert_true(max(abs(cp02_serial_loading_matrix("AR05") - L05_ind)) < 1e-15,
              "cp02_serial_loading_matrix('AR05') matches L05_ind")
  assert_true(max(abs(cp02_serial_loading_matrix("AR1") - L1_ind)) == 0,
              "cp02_serial_loading_matrix('AR1') matches L1_ind")
  assert_true(max(abs(cp02_serial_covariance_matrix("AR05") - Sigma05_ind)) < 1e-15,
              "cp02_serial_covariance_matrix('AR05') matches Sigma05_ind")
  assert_true(max(abs(cp02_serial_covariance_matrix("AR1") - Sigma1_ind)) == 0,
              "cp02_serial_covariance_matrix('AR1') matches Sigma1_ind")

  # Verify noise moments (period variances and pair-difference variances).
  nm05 <- cp02_serial_noise_moments("AR05")
  nm1 <- cp02_serial_noise_moments("AR1")
  assert_true(max(abs(nm05$period_moments$var_idiosyncratic - rep(4 / 3, 5L))) < 1e-15,
              "AR05 marginal period variances all equal 4/3")
  assert_true(max(abs(nm05$period_moments$fourth_central_idiosyncratic - rep(16 / 3, 5L))) < 1e-14,
              "AR05 fourth central moments all equal 3 * (4/3)^2 = 16/3")
  assert_true(max(abs(nm1$period_moments$var_idiosyncratic - (1:5))) == 0,
              "AR1 marginal period variances equal 1:5")
  assert_true(max(abs(nm1$period_moments$fourth_central_idiosyncratic - 3 * ((1:5)^2))) == 0,
              "AR1 fourth central moments equal 3 * t^2")

  pairs <- utils::combn(5L, 2L)
  expected_diff_var_05 <- 2 * (1 - 0.5^(pairs[2L, ] - pairs[1L, ])) / (1 - 0.5^2)
  expected_diff_var_1 <- as.numeric(pairs[2L, ] - pairs[1L, ])
  assert_true(max(abs(nm05$pair_difference_moments$var_diff_idiosyncratic - expected_diff_var_05)) < 1e-15,
              "AR05 pair-difference variances match 2*(1 - rho^(t-s))/(1 - rho^2)")
  assert_true(max(abs(nm1$pair_difference_moments$var_diff_idiosyncratic - expected_diff_var_1)) == 0,
              "AR1 pair-difference variances match (t - s)")

  # Verify step-by-step recurrence against matrix multiplication Z_mat %*% t(L_ind)
  # on (a) the 5x5 identity basis I_5 and (b) a 32-unit deterministic innovation matrix.
  I5 <- diag(5L)
  assert_true(max(abs(cp02_serial_transform_innovations(I5, "AR05") - t(L05_ind))) < 1e-15,
              "AR05 recurrence on I_5 equals t(L05_ind)")
  assert_true(max(abs(cp02_serial_transform_innovations(I5, "AR1") - t(L1_ind))) == 0,
              "AR1 recurrence on I_5 equals t(L1_ind)")

  # Deterministic 32 x 5 innovation matrix with varied exact rational entries.
  z_vec <- ((seq_len(160L) * 17L + 5L) %% 29L - 14L) / 7.0
  Z_mat <- matrix(z_vec, nrow = 32L, ncol = 5L, byrow = TRUE)

  E05_mat <- cp02_serial_transform_innovations(Z_mat, "AR05")
  E05_vec <- cp02_serial_transform_innovations(z_vec, "AR05")
  E05_expected_mat <- Z_mat %*% t(L05_ind)
  assert_true(max(abs(E05_mat - E05_expected_mat)) < 1e-14,
              "AR05 matrix recurrence matches Z_mat %*% t(L05_ind)")
  assert_true(max(abs(E05_vec - as.vector(t(E05_expected_mat)))) < 1e-14,
              "AR05 vector recurrence matches row-by-unit unpacked Z_mat %*% t(L05_ind)")

  E1_mat <- cp02_serial_transform_innovations(Z_mat, "AR1")
  E1_vec <- cp02_serial_transform_innovations(z_vec, "AR1")
  E1_expected_mat <- Z_mat %*% t(L1_ind)
  assert_true(max(abs(E1_mat - E1_expected_mat)) < 1e-14,
              "AR1 matrix recurrence matches Z_mat %*% t(L1_ind)")
  assert_true(max(abs(E1_vec - as.vector(t(E1_expected_mat)))) < 1e-14,
              "AR1 vector recurrence matches row-by-unit unpacked Z_mat %*% t(L1_ind)")

  v05 <- cp02_serial_verify_recurrence_and_covariance("AR05", eps_raw = Z_mat, eps_transformed = E05_mat)
  v1 <- cp02_serial_verify_recurrence_and_covariance("AR1", eps_raw = Z_mat, eps_transformed = E1_mat)
  assert_true(isTRUE(v05$verified) && isTRUE(v1$verified),
              "verification helper passes true AR05 and AR1 recurrence and covariance")

  assert_no_rng_touched("test_recurrence_loading_and_covariance")
}

# ---------------------------------------------------------------------------
# 3. Deliberate wrong-start (AR05) and rho1-stationary (AR1) mutant rejection
# ---------------------------------------------------------------------------
test_mutant_rejection_wrong_start_and_rho1_stationary <- function() {
  z_vec <- ((seq_len(80L) * 13L + 3L) %% 23L - 11L) / 5.0
  Z_mat <- matrix(z_vec, nrow = 16L, ncol = 5L, byrow = TRUE)

  # Mutant 1: Wrong-start AR05 mutant (starts eps_1 ~ N(0, 1) without 1/sqrt(1 - rho^2)).
  L05_wrong_start <- cp02_serial_loading_matrix("AR05")
  L05_wrong_start[, 1L] <- 0.5^(0:4)
  Sigma05_wrong_start <- L05_wrong_start %*% t(L05_wrong_start)
  E05_wrong_start <- Z_mat %*% t(L05_wrong_start)

  # Verify that the wrong-start diagonal is (1, 1.25, 1.3125, 1.328125, 1.33203125) != 4/3.
  assert_true(abs(Sigma05_wrong_start[1L, 1L] - 1.0) < 1e-15 &&
                abs(Sigma05_wrong_start[1L, 1L] - 4 / 3) > 0.33,
              "AR05 wrong-start mutant has non-stationary period-1 variance 1 != 4/3")
  assert_error(
    cp02_serial_verify_recurrence_and_covariance(
      "AR05", cov_matrix = Sigma05_wrong_start
    ),
    "Supplied `cov_matrix` rejected",
    "AR05 wrong-start mutant covariance rejection"
  )
  assert_error(
    cp02_serial_verify_recurrence_and_covariance(
      "AR05", eps_raw = Z_mat, eps_transformed = E05_wrong_start
    ),
    "Transformed innovations rejected",
    "AR05 wrong-start mutant innovation recurrence rejection"
  )
  bad_spec_05 <- cp02_serial_dgp_spec("AR05")
  bad_spec_05$loading_matrix <- L05_wrong_start
  assert_error(
    cp02_validate_serial_spec(bad_spec_05),
    "loading_matrix",
    "AR05 wrong-start loading_matrix spec rejection"
  )

  # Mutant 2a: rho = 1 naive stationary formula 1 / sqrt(1 - rho^2) -> Inf / NaN.
  rho1 <- 1.0
  s1_stationary_at_1 <- 1.0 / sqrt(1.0 - rho1^2) # Inf
  E1_stationary_inf <- matrix(0.0, nrow = nrow(Z_mat), ncol = 5L)
  E1_stationary_inf[, 1L] <- Z_mat[, 1L] * s1_stationary_at_1
  for (tt in 2:5) {
    E1_stationary_inf[, tt] <- rho1 * E1_stationary_inf[, tt - 1L] + Z_mat[, tt]
  }
  Sigma1_stationary_inf <- matrix(s1_stationary_at_1, nrow = 5L, ncol = 5L)
  assert_error(
    cp02_serial_verify_recurrence_and_covariance(
      "AR1", cov_matrix = Sigma1_stationary_inf
    ),
    "non-finite",
    "AR1 stationary 1/sqrt(1-rho^2) non-finite covariance rejection"
  )
  assert_error(
    cp02_serial_verify_recurrence_and_covariance(
      "AR1", eps_raw = Z_mat, eps_transformed = E1_stationary_inf
    ),
    "non-finite",
    "AR1 stationary 1/sqrt(1-rho^2) non-finite innovation rejection"
  )

  # Mutant 2b: rho = 1 finite wrong-stationary-start mutant (starts with AR05
  # stationary variance 4/3 at t = 1 and accumulates unit-root steps, giving Var = t + 1/3 != t).
  L1_wrong_stationary <- cp02_serial_loading_matrix("AR1")
  L1_wrong_stationary[, 1L] <- 1.0 / sqrt(1.0 - 0.5^2)
  Sigma1_wrong_stationary <- L1_wrong_stationary %*% t(L1_wrong_stationary)
  E1_wrong_stationary <- Z_mat %*% t(L1_wrong_stationary)
  assert_error(
    cp02_serial_verify_recurrence_and_covariance(
      "AR1", cov_matrix = Sigma1_wrong_stationary
    ),
    "Supplied `cov_matrix` rejected",
    "AR1 finite wrong-stationary-start covariance rejection"
  )
  assert_error(
    cp02_serial_verify_recurrence_and_covariance(
      "AR1", eps_raw = Z_mat, eps_transformed = E1_wrong_stationary
    ),
    "Transformed innovations rejected",
    "AR1 finite wrong-stationary-start recurrence rejection"
  )
  bad_spec_1 <- cp02_serial_dgp_spec("AR1")
  bad_spec_1$loading_matrix <- L1_wrong_stationary
  assert_error(
    cp02_validate_serial_spec(bad_spec_1),
    "loading_matrix",
    "AR1 wrong-start loading_matrix spec rejection"
  )

  assert_no_rng_touched("test_mutant_rejection_wrong_start_and_rho1_stationary")
}

# ---------------------------------------------------------------------------
# 4. Deterministic 32-unit enumeration fake accepted rotating draw object:
#    transformation, non-mutation, literal mask retention, treatment-effect
#    identity, and oracle/observed outcome & endpoint consistency
# ---------------------------------------------------------------------------
build_deterministic_rotating_base_draw <- function() {
  target_spec <- sim_spec(x_trend = 0)
  vspec <- cp02_baseline_variant_spec(variant = "rotating", base_spec = target_spec)
  pat <- vspec$patterns
  periods <- target_spec$periods
  T_len <- length(periods)
  codes <- c(0, target_spec$cohorts)

  # Complete 32-cell enumeration: 2 (x in 0,1) x 4 (G in 0,3,4,5) x 4 (rotating patterns 1..4)
  grid <- expand.grid(
    pidx = 1:4,
    G = codes,
    x = c(0L, 1L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n <- nrow(grid) # 32
  m <- n * T_len   # 160

  x <- as.integer(grid$x)
  G <- as.numeric(grid$G)
  pidx <- as.integer(grid$pidx)
  z <- ((seq_len(n) * 7L - 3L) %% 19L - 9L) / 4.0
  unit_effect <- ((seq_len(n) * 11L + 1L) %% 17L - 8L) / 5.0
  n_cl <- target_spec$n_clusters
  cluster <- as.integer(((seq_len(n) - 1L) %% n_cl) + 1L)
  cluster_obs_flag <- rep(TRUE, n_cl)
  cluster_shocks <- matrix(0.0, nrow = n_cl, ncol = T_len)
  sample_weight <- rep(1.0, n)
  norm_weight <- rep(1.0, n)

  tilt <- cp02_baseline_tilt_table(vspec)
  cell_idx <- cbind(x + 1L, match(G, codes))
  inclusion_prob_unit <- unname(tilt$inclusion_prob_matrix[cell_idx])
  p_target_unit <- unname(tilt$p_target_matrix[cell_idx])
  p_sample_unit <- unname(tilt$p_sample_matrix[cell_idx])

  id <- rep(seq_len(n), each = T_len)
  period <- rep(periods, times = n)
  t_pos <- match(period, periods)
  Gi <- G[id]
  xi <- x[id]
  gk <- as.character(Gi)

  eps_raw <- ((seq_len(m) * 23L + 7L) %% 31L - 15L) / 6.0
  eps_scale <- rep(1.0, m)
  eps <- eps_scale * eps_raw
  c_shock <- rep(0.0, m)
  shift <- rep(0.0, m)

  y0 <- unname(target_spec$cohort_level[gk]) +
    unname(target_spec$time_effect[as.character(period)]) +
    target_spec$x_trend * xi * (period - periods[1L]) +
    target_spec$z_coef * z[id] +
    unit_effect[id] +
    c_shock +
    eps +
    shift

  treated <- Gi > 0 & period >= Gi
  effect <- numeric(m)
  effect[treated] <- unname(target_spec$theta[gk[treated]]) +
    target_spec$delta * (period[treated] - Gi[treated]) +
    target_spec$kappa * xi[treated]
  y_g <- y0 + effect

  seen <- matrix(FALSE, nrow = nrow(pat), ncol = T_len)
  for (j in seq_len(nrow(pat))) {
    dates_j <- c(pat$date1[j], pat$date2[j])
    seen[j, match(dates_j, periods)] <- TRUE
  }
  observed <- seen[cbind(pidx[id], t_pos)]
  n_obs_unit <- rowSums(seen)[pidx]
  eligible_ge2_unit <- n_obs_unit >= 2L
  unit_cl_obs <- cluster_obs_flag[cluster]

  units <- data.frame(
    id = seq_len(n),
    G = G,
    x = x,
    z = z,
    cluster = cluster,
    cluster_obs_flag = unit_cl_obs,
    sample_weight = sample_weight,
    norm_weight = norm_weight,
    pattern = pat$pattern[pidx],
    pattern_type = pat$type[pidx],
    n_obs = as.integer(n_obs_unit),
    eligible_ge2 = eligible_ge2_unit,
    stringsAsFactors = FALSE
  )

  panel <- data.frame(
    id = id,
    period = period,
    G = Gi,
    x = xi,
    z = z[id],
    cluster = cluster[id],
    cluster_obs_flag = unit_cl_obs[id],
    sample_weight = sample_weight[id],
    norm_weight = norm_weight[id],
    pattern = units$pattern[id],
    pattern_type = units$pattern_type[id],
    n_obs = units$n_obs[id],
    eligible_ge2 = units$eligible_ge2[id],
    observed = observed,
    y = ifelse(observed, y_g, NA_real_),
    stringsAsFactors = FALSE
  )

  y_mat <- matrix(panel$y, nrow = n, ncol = T_len, byrow = TRUE)
  d1 <- pat$date1[pidx]
  d2 <- pat$date2[pidx]
  d1_pos <- match(d1, periods)
  d2_pos <- match(d2, periods)
  y1_val <- y_mat[cbind(seq_len(n), d1_pos)]
  y2_val <- y_mat[cbind(seq_len(n), d2_pos)]
  endpoints <- data.frame(
    id = seq_len(n),
    G = G,
    x = x,
    cluster = cluster,
    sample_weight = sample_weight,
    pattern = units$pattern,
    pattern_type = units$pattern_type,
    n_obs = units$n_obs,
    eligible_ge2 = units$eligible_ge2,
    date1 = d1,
    date2 = d2,
    y_date1 = y1_val,
    y_date2 = y2_val,
    dy = y2_val - y1_val,
    stringsAsFactors = FALSE
  )

  cl_ids <- seq_len(n_cl)
  cl_n_units <- as.integer(tabulate(cluster, nbins = n_cl))
  cl_n_elig <- as.integer(tabulate(cluster[eligible_ge2_unit], nbins = n_cl))
  cl_n_sing <- as.integer(tabulate(cluster[n_obs_unit == 1L], nbins = n_cl))
  cl_n_none <- as.integer(tabulate(cluster[n_obs_unit == 0L], nbins = n_cl))
  clusters_df <- data.frame(
    cluster = cl_ids,
    cluster_obs_flag = cluster_obs_flag,
    n_units = cl_n_units,
    n_eligible_ge2 = cl_n_elig,
    n_singleton = cl_n_sing,
    n_none = cl_n_none,
    zero_eligible = cl_n_elig == 0L,
    stringsAsFactors = FALSE
  )

  counts <- as.data.frame(
    table(
      G = factor(G, levels = codes),
      pattern = factor(units$pattern, levels = pat$pattern)
    ),
    responseName = "units",
    stringsAsFactors = FALSE
  )
  counts$G <- as.numeric(as.character(counts$G))

  colnames(cluster_shocks) <- paste0("shock_t", periods)
  oracle_clusters <- cbind(
    data.frame(cluster = cl_ids, cluster_obs_flag = cluster_obs_flag),
    as.data.frame(cluster_shocks)
  )

  settings <- vspec
  settings$n <- n
  settings$n_clusters_realized <- n_cl
  settings$seed <- 620100001L
  settings$rng <- list(
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection",
    r_version = R.version.string,
    draw_order = c(
      "runif(n): binary X",
      "runif(n): cohort G | X",
      "rnorm(n): covariate Z",
      "rnorm(n): unit effect U",
      "sample.int(n_clusters, n, replace = TRUE): cluster labels",
      "rnorm(n * T): Gaussian period innovations",
      "runif(n): base observation pattern"
    )
  )

  list(
    data = list(
      panel = panel,
      units = units,
      endpoints = endpoints,
      clusters = clusters_df,
      counts = counts
    ),
    oracle = list(
      panel = data.frame(
        id = id,
        period = period,
        y0 = y0,
        y_g = y_g,
        effect = effect,
        eps = eps,
        eps_raw = eps_raw,
        eps_scale = eps_scale,
        cluster_shock = c_shock,
        shift = shift
      ),
      units = data.frame(
        id = seq_len(n),
        unit_effect = unit_effect,
        base_pattern = pat$pattern[pidx],
        inclusion_prob = inclusion_prob_unit,
        p_target_cell = p_target_unit,
        p_sample_cell = p_sample_unit
      ),
      clusters = oracle_clusters,
      truth = cp02_baseline_variant_truth(vspec)
    ),
    settings = settings
  )
}

test_tiny_enumeration_fake_draw_transform <- function() {
  base_draw <- build_deterministic_rotating_base_draw()
  base_snapshot <- unserialize(serialize(base_draw, NULL))

  for (law in c("AR05", "AR1")) {
    out <- cp02_serial_transform_draw(base_draw, law)

    # (a) Original input not mutated.
    assert_true(identical(base_draw, base_snapshot),
                paste(law, "did not mutate original base_draw"))

    # (b) Independent rotating mask and unit/cluster/count metadata retained literally.
    assert_true(identical(out$data$units, base_draw$data$units),
                paste(law, "retains data$units literally"))
    assert_true(identical(out$data$clusters, base_draw$data$clusters),
                paste(law, "retains data$clusters literally"))
    assert_true(identical(out$data$counts, base_draw$data$counts),
                paste(law, "retains data$counts literally"))
    assert_true(identical(out$oracle$units, base_draw$oracle$units),
                paste(law, "retains oracle$units literally"))
    assert_true(identical(out$oracle$clusters, base_draw$oracle$clusters),
                paste(law, "retains oracle$clusters literally"))

    panel_meta_cols <- c("id", "period", "G", "x", "z", "cluster",
                         "cluster_obs_flag", "sample_weight", "norm_weight",
                         "pattern", "pattern_type", "n_obs", "eligible_ge2",
                         "observed")
    assert_true(identical(out$data$panel[, panel_meta_cols],
                          base_draw$data$panel[, panel_meta_cols]),
                paste(law, "retains all data$panel non-outcome columns and mask literally"))

    end_meta_cols <- c("id", "G", "x", "cluster", "sample_weight", "pattern",
                       "pattern_type", "n_obs", "eligible_ge2", "date1", "date2")
    assert_true(identical(out$data$endpoints[, end_meta_cols],
                          base_draw$data$endpoints[, end_meta_cols]),
                paste(law, "retains all data$endpoints metadata and dates literally"))

    # (c) Unchanged treatment-effect identity.
    assert_true(identical(out$oracle$panel$effect, base_draw$oracle$panel$effect),
                paste(law, "retains oracle$panel$effect literally"))
    assert_true(max(abs((out$oracle$panel$y_g - out$oracle$panel$y0) -
                          base_draw$oracle$panel$effect)) < 1e-15,
                paste(law, "preserves y_g - y0 == effect across all 160 rows"))
    pre_or_nt <- out$data$panel$G == 0 | out$data$panel$period < out$data$panel$G
    assert_true(all(out$oracle$panel$y_g[pre_or_nt] == out$oracle$panel$y0[pre_or_nt]),
                paste(law, "has y_g == y0 on all pre-treatment and never-treated cells"))

    # (d) Oracle innovations, serial errors, y0, and y_g updated consistently.
    assert_true(identical(out$oracle$panel$eps_raw, base_draw$oracle$panel$eps_raw),
                paste(law, "retains oracle$panel$eps_raw literally"))
    L_ind <- cp02_serial_loading_matrix(law)
    Z_mat <- matrix(base_draw$oracle$panel$eps_raw, nrow = 32L, ncol = 5L, byrow = TRUE)
    expected_eps <- as.vector(t(Z_mat %*% t(L_ind)))
    assert_true(max(abs(out$oracle$panel$eps - expected_eps)) < 1e-14,
                paste(law, "oracle$panel$eps matches lower-triangular loading transform"))

    delta_eps <- expected_eps - base_draw$oracle$panel$eps
    assert_true(max(abs((out$oracle$panel$y0 - base_draw$oracle$panel$y0) - delta_eps)) < 1e-14,
                paste(law, "oracle$panel$y0 shifts by exact serial error difference"))
    assert_true(max(abs((out$oracle$panel$y_g - base_draw$oracle$panel$y_g) - delta_eps)) < 1e-14,
                paste(law, "oracle$panel$y_g shifts by exact serial error difference"))

    # (e) Observed panel outcome y and endpoints y_date1, y_date2, dy consistency.
    obs <- out$data$panel$observed
    assert_true(all(is.na(out$data$panel$y[!obs])),
                paste(law, "data$panel$y is NA whenever !observed"))
    assert_true(max(abs(out$data$panel$y[obs] - out$oracle$panel$y_g[obs])) == 0,
                paste(law, "data$panel$y equals oracle$panel$y_g whenever observed"))

    for (i in seq_len(32L)) {
      rows_i <- which(out$data$panel$id == i)
      obs_rows_i <- rows_i[out$data$panel$observed[rows_i]]
      assert_true(length(obs_rows_i) == 2L,
                  paste(law, "unit", i, "has exactly 2 observed rows"))
      d1 <- out$data$endpoints$date1[i]
      d2 <- out$data$endpoints$date2[i]
      assert_true(identical(as.integer(out$data$panel$period[obs_rows_i]), c(d1, d2)),
                  paste(law, "unit", i, "observed periods match date1, date2"))
      assert_true(identical(out$data$endpoints$y_date1[i], out$data$panel$y[obs_rows_i[1L]]) &&
                    identical(out$data$endpoints$y_date2[i], out$data$panel$y[obs_rows_i[2L]]),
                  paste(law, "unit", i, "endpoint outcomes match panel$y at date1, date2"))
      assert_true(abs(out$data$endpoints$dy[i] -
                        (out$data$endpoints$y_date2[i] - out$data$endpoints$y_date1[i])) < 1e-15,
                  paste(law, "unit", i, "endpoint dy equals y_date2 - y_date1"))
    }

    # (f) Truth and settings metadata on transformed draw.
    assert_true(identical(out$settings$law, law) &&
                  identical(out$settings$rho, if (law == "AR05") 0.5 else 1.0) &&
                  identical(out$settings$innovation_variance, 1.0) &&
                  identical(out$settings$pattern_family, "rotating"),
                paste(law, "transformed settings metadata is truthful"))
  }

  # Refusal when base_draw has undeclared design or corrupted structure.
  bad_draw_x <- base_draw
  bad_draw_x$settings$target_spec$x_trend <- 0.5
  assert_error(cp02_serial_transform_draw(bad_draw_x, "AR05"),
               "base_draw", "refuse base_draw with x_trend = 0.5")

  bad_draw_pat <- base_draw
  bad_draw_pat$settings$variant <- "balanced"
  assert_error(cp02_serial_transform_draw(bad_draw_pat, "AR1"),
               "base_draw", "refuse base_draw with balanced variant")

  bad_draw_mask <- base_draw
  bad_draw_mask$data$panel$observed[1L] <- !bad_draw_mask$data$panel$observed[1L]
  assert_error(cp02_serial_transform_draw(bad_draw_mask, "AR05"),
               "observed", "refuse base_draw with corrupted rotating mask")

  bad_draw_y0 <- base_draw
  bad_draw_y0$oracle$panel$y0[1L] <- bad_draw_y0$oracle$panel$y0[1L] + 1.0
  assert_error(cp02_serial_transform_draw(bad_draw_y0, "AR1"),
               "y0", "refuse base_draw with inconsistent oracle y0")

  assert_no_rng_touched("test_tiny_enumeration_fake_draw_transform")
}

# ---------------------------------------------------------------------------
# 5. Independent finite cohort/X population arithmetic verification of
#    analytic ATT, contrasts, and all four aggregations (never calling sim_truth)
# ---------------------------------------------------------------------------
test_independent_finite_population_truth_arithmetic <- function() {
  # Finite (X, G) population law for Law A / AR05 / AR1:
  # X in {0, 1} with P(X = 1) = 1/2, P(X = 0) = 1/2.
  # G in {0, 3, 4, 5}:
  #   P(G | X = 0) = (2/5, 1/5, 1/5, 1/5)
  #   P(G | X = 1) = (1/10, 3/10, 3/10, 3/10)
  p_x0_g <- c("0" = 0.5 * 0.40, "3" = 0.5 * 0.20, "4" = 0.5 * 0.20, "5" = 0.5 * 0.20)
  p_x1_g <- c("0" = 0.5 * 0.10, "3" = 0.5 * 0.30, "4" = 0.5 * 0.30, "5" = 0.5 * 0.30)
  pi_g <- p_x0_g + p_x1_g
  mean_x_g <- p_x1_g / pi_g

  assert_true(max(abs(pi_g - rep(0.25, 4L))) < 1e-15,
              "independent pi_g equals (1/4, 1/4, 1/4, 1/4)")
  assert_true(max(abs(mean_x_g - c("0" = 0.20, "3" = 0.60, "4" = 0.60, "5" = 0.60))) < 1e-15,
              "independent E[X | G = g] equals (1/5, 3/5, 3/5, 3/5)")

  # Finite-stratum potential outcome difference E[Y_it(g) - Y_it(0) | G = g]:
  # For x in {0, 1}, tau(g, t, x) = 1{t >= g} * (theta_g + delta * (t - g) + kappa * x)
  theta_map <- c("3" = 1.0, "4" = 2.0, "5" = 1.5)
  delta_val <- 0.5
  kappa_val <- 1.0

  att_finite_stratum <- function(g, tt) {
    if (tt < g) return(0.0)
    g_chr <- as.character(g)
    w0 <- p_x0_g[g_chr] / pi_g[g_chr]
    w1 <- p_x1_g[g_chr] / pi_g[g_chr]
    tau_x0 <- unname(theta_map[g_chr]) + delta_val * (tt - g) + kappa_val * 0.0
    tau_x1 <- unname(theta_map[g_chr]) + delta_val * (tt - g) + kappa_val * 1.0
    unname(w0 * tau_x0 + w1 * tau_x1)
  }

  # Exact rational truths for the 6 post-treatment levels in declared order:
  # (3,3)=8/5, (3,4)=21/10, (3,5)=13/5, (4,4)=13/5, (4,5)=31/10, (5,5)=21/10.
  expected_post_levels <- data.frame(
    group = c(3, 3, 3, 4, 4, 5),
    t = c(3, 4, 5, 4, 5, 5),
    e = c(0, 1, 2, 0, 1, 0),
    att = c(8 / 5, 21 / 10, 13 / 5, 13 / 5, 31 / 10, 21 / 10)
  )
  for (r in seq_len(nrow(expected_post_levels))) {
    g <- expected_post_levels$group[r]
    tt <- expected_post_levels$t[r]
    assert_true(abs(att_finite_stratum(g, tt) - expected_post_levels$att[r]) < 1e-15,
                paste("finite-stratum ATT matches rational for", g, tt))
  }

  # Exact rational truths for the 3 shared-date contrasts:
  # ATT(3,4)-ATT(4,4) = 21/10 - 13/5 = -1/2
  # ATT(3,5)-ATT(4,5) = 13/5 - 31/10 = -1/2
  # ATT(4,5)-ATT(5,5) = 31/10 - 21/10 = 1
  expected_contrasts <- c(-1 / 2, -1 / 2, 1.0)

  # Exact rational truths for the 4 population aggregations using pi_g:
  w_post <- pi_g[as.character(expected_post_levels$group)]
  expected_simple <- sum(w_post * expected_post_levels$att) / sum(w_post)
  assert_true(abs(expected_simple - 47 / 20) < 1e-15,
              "independent simple aggregate equals 47/20 = 2.35")

  expected_group_by <- c(
    mean(expected_post_levels$att[expected_post_levels$group == 3]),
    mean(expected_post_levels$att[expected_post_levels$group == 4]),
    mean(expected_post_levels$att[expected_post_levels$group == 5])
  )
  assert_true(max(abs(expected_group_by - c(21 / 10, 57 / 20, 21 / 10))) < 1e-15,
              "independent group by-cohort ATTs equal (21/10, 57/20, 21/10)")
  wg <- pi_g[c("3", "4", "5")]
  expected_group_overall <- sum(wg * expected_group_by) / sum(wg)
  assert_true(abs(expected_group_overall - 47 / 20) < 1e-15,
              "independent group overall equals 47/20 = 2.35")

  expected_dyn_by <- c(
    "-4" = 0.0, "-3" = 0.0, "-2" = 0.0, "-1" = 0.0,
    "0" = sum(pi_g[c("3", "4", "5")] * c(8 / 5, 13 / 5, 21 / 10)) / sum(pi_g[c("3", "4", "5")]),
    "1" = sum(pi_g[c("3", "4")] * c(21 / 10, 31 / 10)) / sum(pi_g[c("3", "4")]),
    "2" = unname((pi_g["3"] * (13 / 5)) / pi_g["3"])
  )
  assert_true(max(abs(expected_dyn_by - c(0, 0, 0, 0, 21 / 10, 13 / 5, 13 / 5))) < 1e-15,
              "independent dynamic event-time curve matches exact rationals")
  expected_dyn_overall <- mean(expected_dyn_by[c("0", "1", "2")])
  assert_true(abs(expected_dyn_overall - 73 / 30) < 1e-15,
              "independent dynamic post overall equals 73/30")

  expected_cal_by <- c(
    "3" = unname((pi_g["3"] * (8 / 5)) / pi_g["3"]),
    "4" = sum(pi_g[c("3", "4")] * c(21 / 10, 13 / 5)) / sum(pi_g[c("3", "4")]),
    "5" = sum(pi_g[c("3", "4", "5")] * c(13 / 5, 31 / 10, 21 / 10)) / sum(pi_g[c("3", "4", "5")])
  )
  assert_true(max(abs(expected_cal_by - c(8 / 5, 47 / 20, 13 / 5))) < 1e-15,
              "independent calendar curve equals (8/5, 47/20, 13/5)")
  expected_cal_overall <- mean(expected_cal_by)
  assert_true(abs(expected_cal_overall - 131 / 60) < 1e-15,
              "independent calendar overall equals 131/60")

  for (law in c("AR05", "AR1")) {
    tr <- cp02_serial_dgp_truth(law)
    assert_true(max(abs(tr$cohort$pi - unname(pi_g))) < 1e-15,
                paste(law, "truth cohort pi matches independent pi_g"))
    assert_true(max(abs(tr$cohort$mean_x - unname(mean_x_g))) < 1e-15,
                paste(law, "truth cohort mean_x matches independent E[X | G]"))

    for (r in seq_len(nrow(tr$att_gt))) {
      g <- tr$att_gt$group[r]
      tt <- tr$att_gt$t[r]
      exp_att <- att_finite_stratum(g, tt)
      assert_true(abs(tr$att_gt$att[r] - exp_att) < 1e-15 &&
                    abs(tr$att_gt$causal_att[r] - exp_att) < 1e-15 &&
                    abs(tr$att_gt$reporting_did_truth[r] - exp_att) < 1e-15,
                  paste(law, "att_gt row", r, "matches independent finite-stratum ATT"))
    }
    assert_true(all(tr$att_gt$recoverable_under_pattern[tr$att_gt$post]),
                paste(law, "all 6 post levels are recoverable under rotating pattern"))

    assert_true(identical(as.numeric(tr$post_levels$group), expected_post_levels$group) &&
                  identical(as.numeric(tr$post_levels$t), expected_post_levels$t) &&
                  max(abs(tr$post_levels$att - expected_post_levels$att)) < 1e-15,
                paste(law, "6 post levels match declared order and exact rationals"))

    assert_true(max(abs(tr$shared_date_contrasts$truth - expected_contrasts)) < 1e-15,
                paste(law, "3 shared-date contrasts match (-1/2, -1/2, 1)"))

    assert_true(abs(tr$simple - expected_simple) < 1e-15,
                paste(law, "simple aggregate matches 47/20"))
    assert_true(max(abs(tr$group$by$att - expected_group_by)) < 1e-15 &&
                  abs(tr$group$overall - expected_group_overall) < 1e-15,
                paste(law, "group aggregate matches independent finite arithmetic"))
    assert_true(max(abs(tr$dynamic$by$att - unname(expected_dyn_by))) < 1e-15 &&
                  abs(tr$dynamic$overall - expected_dyn_overall) < 1e-15,
                paste(law, "dynamic aggregate matches independent finite arithmetic"))
    assert_true(max(abs(tr$calendar$by$att - unname(expected_cal_by))) < 1e-15 &&
                  abs(tr$calendar$overall - expected_cal_overall) < 1e-15,
                paste(law, "calendar aggregate matches independent finite arithmetic"))

    # All 24 non-reference pre-reporting operators have causal_att == 0 and
    # reporting_did_truth == 0 under x_trend = 0 and no period-1 shift, and
    # pair_observed_in_pattern is TRUE iff |base_period - t| == 1 under rotating.
    assert_true(nrow(tr$pre_reporting_operators) == 24L &&
                  all(tr$pre_reporting_operators$causal_att == 0) &&
                  all(tr$pre_reporting_operators$reporting_did_truth == 0) &&
                  identical(tr$pre_reporting_operators$pair_observed_in_pattern,
                            abs(tr$pre_reporting_operators$base_period -
                                  tr$pre_reporting_operators$t) == 1L),
                paste(law, "all 24 pre-reporting operator truths are 0 with rotating pair availability"))
  }

  assert_no_rng_touched("test_independent_finite_population_truth_arithmetic")
}

# ---------------------------------------------------------------------------
# 6. Structural verification of future caller RNG restoration via the
#    accepted single with_seed constructor (without calling any RNG now)
# ---------------------------------------------------------------------------
test_structural_single_with_seed_contract_no_rng_executed <- function() {
  # Verify with_seed body structure: saves RNGkind(), checks exists('.Random.seed'),
  # registers on.exit() restoration, and sets L'Ecuyer-CMRG / Inversion / Rejection.
  ws_syms <- all.names(body(with_seed), unique = FALSE)
  assert_true("on.exit" %in% ws_syms &&
                "RNGkind" %in% ws_syms &&
                "set.seed" %in% ws_syms &&
                "exists" %in% ws_syms &&
                "assign" %in% ws_syms &&
                "rm" %in% ws_syms,
              "with_seed body contains RNGkind/set.seed/on.exit restoration primitives")

  # Verify cp02_baseline_variant_draw calls with_seed exactly once and never
  # calls set.seed or RNGkind directly.
  bvd_syms <- all.names(body(cp02_baseline_variant_draw), unique = FALSE)
  assert_true(sum(bvd_syms == "with_seed") == 1L &&
                !("set.seed" %in% bvd_syms) &&
                !("RNGkind" %in% bvd_syms),
              "cp02_baseline_variant_draw delegates RNG isolation to a single with_seed() call")

  # Verify NO function in cp02_serial_dgp.R calls set.seed, RNGkind, with_seed,
  # rnorm, runif, rbinom, rexp, sample, sample.int, or sim_draw.
  serial_fns <- list(
    cp02_serial_dgp_catalog = cp02_serial_dgp_catalog,
    cp02_resolve_serial_law = cp02_resolve_serial_law,
    cp02_serial_loading_matrix = cp02_serial_loading_matrix,
    cp02_serial_covariance_matrix = cp02_serial_covariance_matrix,
    cp02_serial_dgp_seed = cp02_serial_dgp_seed,
    cp02_serial_dgp_spec = cp02_serial_dgp_spec,
    cp02_validate_serial_spec = cp02_validate_serial_spec,
    cp02_as_serial_spec = cp02_as_serial_spec,
    cp02_serial_noise_moments = cp02_serial_noise_moments,
    cp02_serial_dgp_truth = cp02_serial_dgp_truth,
    cp02_serial_transform_innovations = cp02_serial_transform_innovations,
    cp02_serial_verify_recurrence_and_covariance = cp02_serial_verify_recurrence_and_covariance,
    cp02_validate_base_rotating_draw = cp02_validate_base_rotating_draw,
    cp02_serial_transform_draw = cp02_serial_transform_draw,
    cp02_serial_dgp_draw = cp02_serial_dgp_draw
  )
  forbidden_rng_syms <- c("set.seed", "RNGkind", "with_seed", "rnorm", "runif",
                          "rbinom", "rexp", "rt", "sample", "sample.int", "sim_draw")
  for (fn_name in names(serial_fns)) {
    syms <- all.names(body(serial_fns[[fn_name]]), unique = TRUE)
    hit <- intersect(syms, forbidden_rng_syms)
    assert_true(length(hit) == 0L,
                paste(fn_name, "contains zero direct RNG/with_seed calls"))
  }

  # Verify cp02_serial_dgp_draw calls cp02_baseline_variant_draw exactly once
  # and cp02_serial_transform_draw exactly once.
  draw_syms <- all.names(body(cp02_serial_dgp_draw), unique = FALSE)
  assert_true(sum(draw_syms == "cp02_baseline_variant_draw") == 1L &&
                sum(draw_syms == "cp02_serial_transform_draw") == 1L,
              "cp02_serial_dgp_draw makes a single cp02_baseline_variant_draw call then cp02_serial_transform_draw")

  assert_no_rng_touched("test_structural_single_with_seed_contract_no_rng_executed")
}

# Run all deterministic test suites.
test_catalog_spec_and_refusals()
test_recurrence_loading_and_covariance()
test_mutant_rejection_wrong_start_and_rho1_stationary()
test_tiny_enumeration_fake_draw_transform()
test_independent_finite_population_truth_arithmetic()
test_structural_single_with_seed_contract_no_rng_executed()

assert_no_rng_touched("completion of all deterministic tests")

cat("ALL DETERMINISTIC CP02 SERIAL DGP TESTS PASSED\n")
cat("RNG_TOUCHED=FALSE (.Random.seed exists: ",
    exists(".Random.seed", envir = globalenv(), inherits = FALSE), ")\n", sep = "")
cat("LIVE_RNG_RESTORATION_EXECUTION=PENDING_ACTUAL_HOST_CHECK\n")
