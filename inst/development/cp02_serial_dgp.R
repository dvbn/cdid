# Development simulation support for cdid (not package API).
# Requires seed.R (with_seed), dgp.R (sim_spec, sim_truth), and
# cp02_baseline_variants_dgp.R (cp02_baseline_variant_spec,
# cp02_baseline_variant_truth, cp02_baseline_variant_draw,
# cp02_check_variant_n).
# Does not clone or modify seed.R, dgp.R, or cp02_baseline_variants_dgp.R.
#
# Scoped CP02 serial-error DGP bindings AR05 and AR1 from
# reference/CP02_FINAL_MC_DESIGN.md:
#   - AR05 : rotating two-consecutive-date windows on periods 1:5, cohorts
#            3,4,5 plus never-treated (0), accepted sim_spec default parameters
#            with x_trend = 0 (law A outcome specification), Gaussian
#            innovation variance 1, rho = 0.5 with stationary start
#            eps_{i,1} ~ N(0, 1 / (1 - rho^2)) and recurrence
#            eps_{i,t} = rho * eps_{i,t-1} + nu_{i,t}, nu_{i,t} ~ iid N(0, 1).
#   - AR1  : same rotating x_trend = 0 population and masks, Gaussian
#            innovation variance 1, rho = 1 with unit-variance start
#            eps_{i,1} ~ N(0, 1) and cumulative innovations
#            eps_{i,t} = eps_{i,t-1} + nu_{i,t} (finite five-date random walk,
#            never stationary AR(1) at rho = 1).
#
# Reuses cp02_baseline_variant_draw() for the rotating x_trend = 0 base draw
# inside a single with_seed() call and applies a deterministic row-by-unit
# lower-triangular AR transform to oracle$panel$eps_raw, updating only outcome
# errors (oracle$panel$eps, oracle$panel$y0, oracle$panel$y_g) and observed
# outcomes (data$panel$y, data$endpoints$y_date1, y_date2, dy) while preserving
# masks, identifiers, unit/cluster metadata, treatment effects, and accepted
# population truth.

cp02_serial_dgp_catalog <- function() {
  data.frame(
    law = c("AR05", "AR1"),
    design_table_law_index = c(20L, 21L),
    rho = c(0.5, 1.0),
    innovation_family = c("gaussian", "gaussian"),
    innovation_variance = c(1.0, 1.0),
    sd_innovation = c(1.0, 1.0),
    start_convention = c(
      "stationary_normal_0_1_over_1_minus_rho2",
      "random_walk_normal_0_1_start"
    ),
    period1_error_variance = c(1.0 / (1.0 - 0.5^2), 1.0),
    stationary = c(TRUE, FALSE),
    process_type = c("stationary_ar1", "finite_horizon_random_walk"),
    pattern_family = c("rotating", "rotating"),
    x_trend = c(0.0, 0.0),
    campaign_control_group = c("NT", "NT"),
    campaign_xformla = c("~1", "~1"),
    campaign_weighting = c("T", "T"),
    stringsAsFactors = FALSE
  )
}

cp02_resolve_serial_law <- function(law = NULL, rho = NULL) {
  if (is.null(law) && is.null(rho)) {
    stop("Must specify declared serial DGP `law` ('AR05' or 'AR1') ",
         "or `rho` (0.5 or 1).", call. = FALSE)
  }
  if (!is.null(law) && is.numeric(law) && is.null(rho)) {
    rho <- law
    law <- NULL
  }
  if (!is.null(rho)) {
    if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || !is.finite(rho)) {
      stop("Undeclared serial correlation parameter `rho`: must be a single ",
           "finite number equal to 0.5 (AR05) or 1 (AR1).", call. = FALSE)
    }
    if (rho == 0.5) {
      rho_law <- "AR05"
    } else if (rho == 1.0) {
      rho_law <- "AR1"
    } else {
      stop("Undeclared serial correlation parameter `rho` = ", rho,
           ": only rho = 0.5 (AR05) and rho = 1 (AR1) are declared.",
           call. = FALSE)
    }
  } else {
    rho_law <- NULL
  }

  if (!is.null(law)) {
    if (!is.character(law) || length(law) != 1L || is.na(law) ||
        !(law %in% c("AR05", "AR1"))) {
      stop("Undeclared serial DGP `law`: only 'AR05' (rho = 0.5) and ",
           "'AR1' (rho = 1) are declared in CP02_FINAL_MC_DESIGN.md.",
           call. = FALSE)
    }
    if (!is.null(rho_law) && !identical(law, rho_law)) {
      stop("Inconsistent `law` ('", law, "') and `rho` (", rho,
           "): AR05 requires rho = 0.5 and AR1 requires rho = 1.",
           call. = FALSE)
    }
    return(law)
  }
  rho_law
}

cp02_serial_loading_matrix <- function(law_or_spec = NULL, rho = NULL, ...) {
  if (length(list(...)) > 0L) {
    stop("Undeclared design override in `cp02_serial_loading_matrix`.",
         call. = FALSE)
  }
  law <- if (is.list(law_or_spec)) {
    sspec <- cp02_validate_serial_spec(law_or_spec)
    if (!is.null(rho)) cp02_resolve_serial_law(law = sspec$law, rho = rho)
    sspec$law
  } else {
    cp02_resolve_serial_law(law = law_or_spec, rho = rho)
  }
  T_len <- 5L
  L <- matrix(0.0, nrow = T_len, ncol = T_len,
              dimnames = list(as.character(seq_len(T_len)),
                              as.character(seq_len(T_len))))
  if (law == "AR05") {
    rho_val <- 0.5
    s1 <- 1.0 / sqrt(1.0 - rho_val^2)
    for (r in seq_len(T_len)) {
      L[r, 1L] <- (rho_val^(r - 1L)) * s1
      if (r >= 2L) {
        for (c in 2L:r) {
          L[r, c] <- rho_val^(r - c)
        }
      }
    }
  } else {
    # AR1: finite five-date random walk, never stationary AR(1) at rho = 1.
    for (r in seq_len(T_len)) {
      L[r, seq_len(r)] <- 1.0
    }
  }
  L
}

cp02_serial_covariance_matrix <- function(law_or_spec = NULL, rho = NULL, ...) {
  if (length(list(...)) > 0L) {
    stop("Undeclared design override in `cp02_serial_covariance_matrix`.",
         call. = FALSE)
  }
  law <- if (is.list(law_or_spec)) {
    sspec <- cp02_validate_serial_spec(law_or_spec)
    if (!is.null(rho)) cp02_resolve_serial_law(law = sspec$law, rho = rho)
    sspec$law
  } else {
    cp02_resolve_serial_law(law = law_or_spec, rho = rho)
  }
  periods <- 1:5
  grid_s <- matrix(periods, nrow = 5L, ncol = 5L)
  grid_t <- matrix(periods, nrow = 5L, ncol = 5L, byrow = TRUE)
  if (law == "AR05") {
    rho_val <- 0.5
    Sigma <- (rho_val^abs(grid_s - grid_t)) / (1.0 - rho_val^2)
  } else {
    Sigma <- pmin(grid_s, grid_t) * 1.0
  }
  dimnames(Sigma) <- list(as.character(periods), as.character(periods))
  Sigma
}

cp02_serial_dgp_seed <- function(
    law_or_spec = NULL,
    stage = c("development", "production",
              "inference_development", "inference_production"),
    rung = 1L,
    rep = 1L,
    rho = NULL) {
  sspec <- cp02_as_serial_spec(law_or_spec = law_or_spec, rho = rho)
  if (!is.character(stage) || length(stage) != 1L || is.na(stage) ||
      !(stage %in% c("development", "production",
                     "inference_development", "inference_production"))) {
    stop("`stage` must be one of 'development', 'production', ",
         "'inference_development', or 'inference_production'.",
         call. = FALSE)
  }
  if (!is.numeric(rung) || length(rung) != 1L || is.na(rung) ||
      !(rung %in% c(1, 2, 3))) {
    stop("`rung` must be 1 (N=2400), 2 (N=9600), or 3 (N=38400).",
         call. = FALSE)
  }
  max_rep <- if (stage %in% c("development", "inference_development")) 200L else 8000L
  if (!is.numeric(rep) || length(rep) != 1L || is.na(rep) || !is.finite(rep) ||
      rep != round(rep) || rep < 1L || rep > max_rep) {
    stop("`rep` must be a whole number in 1..", max_rep, " for stage '",
         stage, "'.", call. = FALSE)
  }
  base_map <- c(
    development = 600000000L,
    production = 700000000L,
    inference_development = 800000000L,
    inference_production = 850000000L
  )
  as.integer(
    unname(base_map[stage]) +
      1000000L * sspec$design_table_law_index +
      100000L * as.integer(rung) +
      as.integer(rep)
  )
}

cp02_serial_dgp_spec <- function(law = NULL, rho = NULL, ...) {
  extra_args <- list(...)
  if (length(extra_args) > 0L) {
    stop("Undeclared design override in `cp02_serial_dgp_spec`: AR05 and AR1 ",
         "are fixed to the rotating x_trend = 0 Gaussian unit-iid design.",
         call. = FALSE)
  }
  law <- cp02_resolve_serial_law(law = law, rho = rho)
  cat_df <- cp02_serial_dgp_catalog()
  cat_row <- cat_df[cat_df$law == law, , drop = FALSE]

  # Law A outcome specification: accepted sim_spec() defaults with x_trend = 0.
  target_spec <- sim_spec(x_trend = 0)
  base_variant_spec <- cp02_baseline_variant_spec(
    variant = "rotating",
    base_spec = target_spec
  )

  rho_val <- cat_row$rho
  L <- cp02_serial_loading_matrix(law)
  Sigma <- cp02_serial_covariance_matrix(law)

  list(
    family = "cp02_serial_error_dgp",
    law = law,
    design_table_law_index = as.integer(cat_row$design_table_law_index),
    status = paste(
      "Development validation DGP binding for CP02_FINAL_MC_DESIGN.md;",
      "actual stochastic constructor/binding checks and campaign execution",
      "are separately frozen and admitted by root."
    ),
    rho = rho_val,
    innovation_family = "gaussian",
    innovation_variance = 1.0,
    sd_innovation = 1.0,
    start_convention = cat_row$start_convention,
    start_distribution = if (law == "AR05") {
      "Normal(0, 1 / (1 - rho^2)) = Normal(0, 4/3)"
    } else {
      "Normal(0, 1)"
    },
    recurrence_equation = if (law == "AR05") {
      "eps_1 = z_1 / sqrt(1 - rho^2); eps_t = 0.5 * eps_{t-1} + z_t for t = 2..5"
    } else {
      "eps_1 = z_1; eps_t = eps_{t-1} + z_t = sum_{s <= t} z_s for t = 2..5"
    },
    stationary = isTRUE(cat_row$stationary),
    process_type = cat_row$process_type,
    eps_scale_semantics = paste(
      "oracle$panel$eps_scale records the per-period Gaussian innovation",
      "standard deviation (1.0); marginal period error standard deviations",
      "sqrt(diag(covariance_matrix)) are recorded in",
      "oracle$truth$noise$period_moments$sd_idiosyncratic."
    ),
    legacy_provenance = list(
      read_only_source = "simulation_audit/sims/common/dgp_functions.R:110-123",
      full_source_sha256 = "fd4baa42a8a1bfb63e851d353201890934f077cce3ecb5f4697ee927e32c1df9",
      role = "equation/initialization convention only; no attrition repair or code port"
    ),
    target_spec = target_spec,
    base_variant_spec = base_variant_spec,
    pattern_family = "rotating",
    patterns = base_variant_spec$patterns,
    x_trend = 0.0,
    campaign_fit_spec = list(
      control_group = "NT",
      xformla = "~1",
      weighting = "T",
      selection = "none",
      pt_assumption = "post",
      alpha = 0.05,
      baselines = c("varying", "universal"),
      rungs_n = c(rung1 = 2400L, rung2 = 9600L, rung3 = 38400L),
      development_attempts_per_rung = 200L,
      production_initial_attempts = 2000L,
      production_extension_attempts = 8000L
    ),
    loading_matrix = L,
    covariance_matrix = Sigma,
    draw_order = c(
      "runif(n): binary X",
      "runif(n): cohort G | X",
      "rnorm(n): covariate Z",
      "rnorm(n): unit effect U",
      "sample.int(n_clusters, n, replace = TRUE): cluster labels",
      "rnorm(n * T): Gaussian period innovations eps_raw (row-by-unit t = 1..5)",
      "runif(n): base observation pattern (rotating 1-2, 2-3, 3-4, 4-5)",
      paste0(
        "deterministic row-by-unit AR transform (", law, ", rho = ", rho_val,
        ", ", cat_row$start_convention,
        ") of oracle$panel$eps_raw -> oracle$panel$eps (no extra RNG draws)"
      )
    ),
    stream_note = paste(
      "Calls cp02_baseline_variant_draw(rotating, x_trend = 0) once inside a",
      "single with_seed(seed, ...) block, then applies a deterministic",
      "row-by-unit AR transform to oracle$panel$eps_raw. No second seeded",
      "generator, resampling, or extra RNG consumption."
    )
  )
}

cp02_validate_serial_spec <- function(sspec) {
  if (!is.list(sspec) || !identical(sspec$family, "cp02_serial_error_dgp")) {
    stop("Invalid serial DGP specification: expected an object from ",
         "`cp02_serial_dgp_spec()`.", call. = FALSE)
  }
  law <- cp02_resolve_serial_law(law = sspec$law, rho = sspec$rho)
  ref <- cp02_serial_dgp_spec(law = law)

  if (!identical(sspec$stationary, ref$stationary) ||
      !identical(sspec$start_convention, ref$start_convention) ||
      !identical(sspec$process_type, ref$process_type)) {
    stop("Invalid start/stationarity metadata in serial DGP specification ",
         "for law '", law, "'.", call. = FALSE)
  }
  if (!is.numeric(sspec$innovation_variance) ||
      length(sspec$innovation_variance) != 1L ||
      !identical(as.numeric(sspec$innovation_variance), 1.0) ||
      !is.numeric(sspec$sd_innovation) ||
      length(sspec$sd_innovation) != 1L ||
      !identical(as.numeric(sspec$sd_innovation), 1.0) ||
      !identical(sspec$innovation_family, "gaussian")) {
    stop("Undeclared innovation specification: AR05 and AR1 require ",
         "Gaussian innovations with variance 1.", call. = FALSE)
  }
  if (!identical(sspec$target_spec, ref$target_spec)) {
    stop("Undeclared `target_spec` modification: AR05 and AR1 require ",
         "accepted `sim_spec(x_trend = 0)` defaults.", call. = FALSE)
  }
  if (!identical(sspec$base_variant_spec, ref$base_variant_spec) ||
      !identical(sspec$pattern_family, "rotating") ||
      !identical(as.numeric(sspec$x_trend), 0.0)) {
    stop("Undeclared base variant design: AR05 and AR1 require the ",
         "rotating independent two-consecutive-date design with x_trend = 0.",
         call. = FALSE)
  }
  # These fixed fields are propagated into seeds and saved draw metadata.
  # Refuse changes rather than accepting a different law namespace or label.
  metadata_fields <- c(
    "design_table_law_index", "status", "start_distribution",
    "recurrence_equation", "eps_scale_semantics", "legacy_provenance",
    "patterns", "campaign_fit_spec", "draw_order", "stream_note"
  )
  if (!identical(sspec[metadata_fields], ref[metadata_fields])) {
    stop("Undeclared serial DGP metadata modification for law '", law,
         "': law index, start, patterns, campaign, and draw-order metadata ",
         "must match the declared specification.", call. = FALSE)
  }
  if (!is.matrix(sspec$loading_matrix) ||
      !identical(dim(sspec$loading_matrix), c(5L, 5L)) ||
      any(!is.finite(sspec$loading_matrix)) ||
      max(abs(sspec$loading_matrix - ref$loading_matrix)) > 1e-12) {
    stop("Serial DGP `loading_matrix` does not match the analytic ",
         "lower-triangular loading matrix for law '", law, "'.",
         call. = FALSE)
  }
  if (!is.matrix(sspec$covariance_matrix) ||
      !identical(dim(sspec$covariance_matrix), c(5L, 5L)) ||
      any(!is.finite(sspec$covariance_matrix)) ||
      max(abs(sspec$covariance_matrix - ref$covariance_matrix)) > 1e-12) {
    stop("Serial DGP `covariance_matrix` does not match the analytic ",
         "finite-horizon covariance matrix for law '", law, "'.",
         call. = FALSE)
  }
  sspec
}

cp02_as_serial_spec <- function(law_or_spec = NULL, rho = NULL, ...) {
  if (length(list(...)) > 0L) {
    stop("Undeclared design override in serial DGP helper: AR05 and AR1 ",
         "are fixed to the rotating x_trend = 0 Gaussian unit-iid design.",
         call. = FALSE)
  }
  if (is.list(law_or_spec)) {
    sspec <- cp02_validate_serial_spec(law_or_spec)
    if (!is.null(rho)) {
      cp02_resolve_serial_law(law = sspec$law, rho = rho)
    }
    return(sspec)
  }
  cp02_serial_dgp_spec(law = law_or_spec, rho = rho)
}

cp02_serial_noise_moments <- function(law_or_spec = NULL, rho = NULL, ...) {
  sspec <- cp02_as_serial_spec(law_or_spec = law_or_spec, rho = rho, ...)
  periods <- sspec$target_spec$periods
  L <- cp02_serial_loading_matrix(sspec)
  Sigma <- cp02_serial_covariance_matrix(sspec)
  var_t <- unname(diag(Sigma))
  sd_t <- sqrt(var_t)

  period_df <- data.frame(
    period = periods,
    sd_idiosyncratic = sd_t,
    var_idiosyncratic = var_t,
    var_cluster_shock = rep(0.0, length(periods)),
    var_total_noise = var_t,
    cond_mean = rep(0.0, length(periods)),
    third_central_idiosyncratic = rep(0.0, length(periods)),
    fourth_central_idiosyncratic = 3.0 * (var_t^2),
    within_cluster_cov_same_period = rep(0.0, length(periods)),
    stringsAsFactors = FALSE
  )

  pairs <- utils::combn(length(periods), 2L)
  s_idx <- pairs[1L, ]
  t_idx <- pairs[2L, ]
  cov_st <- Sigma[cbind(s_idx, t_idx)]
  var_diff_idio <- var_t[t_idx] + var_t[s_idx] - 2.0 * cov_st

  pair_df <- data.frame(
    date1 = periods[s_idx],
    date2 = periods[t_idx],
    cov_idiosyncratic = unname(cov_st),
    cond_mean_diff = rep(0.0, ncol(pairs)),
    var_diff_idiosyncratic = unname(var_diff_idio),
    var_diff_cluster = rep(0.0, ncol(pairs)),
    var_diff_total = unname(var_diff_idio),
    third_central_diff = rep(0.0, ncol(pairs)),
    skewness_diff_idiosyncratic = rep(0.0, ncol(pairs)),
    fourth_central_diff_total = unname(3.0 * (var_diff_idio^2)),
    within_cluster_cov_diff = rep(0.0, ncol(pairs)),
    stringsAsFactors = FALSE
  )

  list(
    law = sspec$law,
    rho = sspec$rho,
    innovation_family = sspec$innovation_family,
    innovation_variance = sspec$innovation_variance,
    start_convention = sspec$start_convention,
    stationary = sspec$stationary,
    finite_moments_all_orders = TRUE,
    loading_matrix = L,
    covariance_matrix = Sigma,
    max_loading_cov_diff = max(abs(L %*% t(L) - Sigma)),
    period_moments = period_df,
    pair_difference_moments = pair_df
  )
}

cp02_serial_dgp_truth <- function(law_or_spec = NULL, rho = NULL, ...) {
  sspec <- cp02_as_serial_spec(law_or_spec = law_or_spec, rho = rho, ...)
  # Population cohort shares, conditional X means, post-treatment ATT(g,t),
  # contrasts, and all four population aggregates equal the accepted A
  # population (sim_spec(x_trend = 0)), never sample or survivor shares.
  base_truth <- cp02_baseline_variant_truth(sspec$base_variant_spec)
  post_cells <- base_truth$att_gt[base_truth$att_gt$post, , drop = FALSE]
  post_order <- order(post_cells$group, post_cells$t)
  post_cells <- post_cells[post_order, , drop = FALSE]
  rownames(post_cells) <- NULL

  att_lookup <- function(g, tt) {
    post_cells$att[post_cells$group == g & post_cells$t == tt]
  }
  contrasts_df <- data.frame(
    contrast = c("ATT(3,4)-ATT(4,4)", "ATT(3,5)-ATT(4,5)", "ATT(4,5)-ATT(5,5)"),
    group1 = c(3, 3, 4),
    group2 = c(4, 4, 5),
    t = c(4, 5, 5),
    truth = c(
      att_lookup(3, 4) - att_lookup(4, 4),
      att_lookup(3, 5) - att_lookup(4, 5),
      att_lookup(4, 5) - att_lookup(5, 5)
    ),
    stringsAsFactors = FALSE
  )

  base_truth$law <- sspec$law
  base_truth$rho <- sspec$rho
  base_truth$post_levels <- post_cells[, c("group", "t", "e", "att")]
  base_truth$shared_date_contrasts <- contrasts_df
  base_truth$noise <- cp02_serial_noise_moments(sspec)
  base_truth$serial_covariance <- sspec$covariance_matrix
  base_truth$serial_loading <- sspec$loading_matrix
  base_truth
}

cp02_serial_transform_innovations <- function(
    eps_raw,
    law_or_spec = NULL,
    rho = NULL,
    ...) {
  sspec <- cp02_as_serial_spec(law_or_spec = law_or_spec, rho = rho, ...)
  T_len <- length(sspec$target_spec$periods)
  if (!is.numeric(eps_raw) || length(eps_raw) < T_len ||
      any(is.na(eps_raw)) || any(!is.finite(eps_raw))) {
    stop("`eps_raw` must be a finite numeric vector or n x 5 matrix.",
         call. = FALSE)
  }

  was_matrix <- is.matrix(eps_raw)
  if (was_matrix) {
    if (ncol(eps_raw) != T_len || nrow(eps_raw) < 1L) {
      stop("Matrix `eps_raw` must have 5 columns (periods 1:5) and >= 1 row.",
           call. = FALSE)
    }
    Z_mat <- eps_raw
    n <- nrow(Z_mat)
  } else {
    if (length(eps_raw) %% T_len != 0L) {
      stop("Vector `eps_raw` length must be a multiple of 5 (periods 1:5).",
           call. = FALSE)
    }
    n <- length(eps_raw) %/% T_len
    Z_mat <- matrix(eps_raw, nrow = n, ncol = T_len, byrow = TRUE)
  }

  # Deterministic row-by-unit recurrence matching legacy_AR_excerpt.txt:
  # innovation variance = 1;
  # AR05 (rho = 0.5): eps[, 1] = z[, 1] / sqrt(1 - rho^2),
  #                   eps[, t] = rho * eps[, t - 1] + z[, t] for t = 2..5.
  # AR1  (rho = 1.0): eps[, 1] = z[, 1],
  #                   eps[, t] = eps[, t - 1] + z[, t] for t = 2..5.
  E_mat <- matrix(0.0, nrow = n, ncol = T_len)
  if (sspec$law == "AR05") {
    rho_val <- 0.5
    E_mat[, 1L] <- Z_mat[, 1L] / sqrt(1.0 - rho_val^2)
    for (tt in 2L:T_len) {
      E_mat[, tt] <- rho_val * E_mat[, tt - 1L] + Z_mat[, tt]
    }
  } else {
    E_mat[, 1L] <- Z_mat[, 1L]
    for (tt in 2L:T_len) {
      E_mat[, tt] <- E_mat[, tt - 1L] + Z_mat[, tt]
    }
  }

  if (was_matrix) {
    colnames(E_mat) <- as.character(sspec$target_spec$periods)
    return(E_mat)
  }
  as.vector(t(E_mat))
}

cp02_serial_verify_recurrence_and_covariance <- function(
    law_or_spec = NULL,
    eps_raw = NULL,
    eps_transformed = NULL,
    cov_matrix = NULL,
    tol = 1e-12,
    rho = NULL,
    ...) {
  sspec <- cp02_as_serial_spec(law_or_spec = law_or_spec, rho = rho, ...)
  periods <- 1:5
  T_len <- length(periods)

  # Independently form lower-triangular loading matrix L and analytic Sigma.
  L_ind <- matrix(0.0, nrow = T_len, ncol = T_len)
  if (sspec$law == "AR05") {
    rho_val <- 0.5
    s1 <- 1.0 / sqrt(1.0 - rho_val^2)
    for (r in seq_len(T_len)) {
      L_ind[r, 1L] <- (rho_val^(r - 1L)) * s1
      if (r >= 2L) {
        for (c in 2L:r) {
          L_ind[r, c] <- rho_val^(r - c)
        }
      }
    }
    Sigma_ind <- outer(periods, periods, function(s, t) {
      (rho_val^abs(s - t)) / (1.0 - rho_val^2)
    })
  } else {
    for (r in seq_len(T_len)) {
      L_ind[r, seq_len(r)] <- 1.0
    }
    Sigma_ind <- outer(periods, periods, pmin) * 1.0
  }

  cov_from_L <- L_ind %*% t(L_ind)
  max_L_cov_err <- max(abs(cov_from_L - Sigma_ind))
  if (!is.finite(max_L_cov_err) || max_L_cov_err > tol) {
    stop("Analytic loading covariance L %*% t(L) failed to match target ",
         "covariance for law '", sspec$law, "'.", call. = FALSE)
  }

  if (!is.null(cov_matrix)) {
    if (!is.matrix(cov_matrix) || !identical(dim(cov_matrix), c(T_len, T_len)) ||
        any(is.na(cov_matrix)) || any(!is.finite(cov_matrix))) {
      stop("Supplied `cov_matrix` is non-finite or has wrong dimensions.",
           call. = FALSE)
    }
    max_supplied_cov_err <- max(abs(cov_matrix - Sigma_ind))
    if (max_supplied_cov_err > tol) {
      stop("Supplied `cov_matrix` rejected for law '", sspec$law,
           "': max abs diff = ", signif(max_supplied_cov_err, 6), ".",
           call. = FALSE)
    }
  } else {
    max_supplied_cov_err <- max(abs(sspec$covariance_matrix - Sigma_ind))
    if (max_supplied_cov_err > tol) {
      stop("Spec `covariance_matrix` rejected for law '", sspec$law, "'.",
           call. = FALSE)
    }
  }

  max_recurrence_err <- 0.0
  if (!is.null(eps_raw)) {
    Z_mat <- if (is.matrix(eps_raw)) {
      eps_raw
    } else {
      matrix(eps_raw, ncol = T_len, byrow = TRUE)
    }
    expected_mat <- Z_mat %*% t(L_ind)
    actual_mat <- if (is.null(eps_transformed)) {
      cp02_serial_transform_innovations(Z_mat, sspec)
    } else if (is.matrix(eps_transformed)) {
      eps_transformed
    } else {
      matrix(eps_transformed, ncol = T_len, byrow = TRUE)
    }
    if (!identical(dim(actual_mat), dim(expected_mat)) ||
        any(is.na(actual_mat)) || any(!is.finite(actual_mat))) {
      stop("Transformed innovations are non-finite or have wrong dimensions.",
           call. = FALSE)
    }
    max_recurrence_err <- max(abs(actual_mat - expected_mat))
    if (max_recurrence_err > tol) {
      stop("Transformed innovations rejected against lower-triangular loading ",
           "matrix for law '", sspec$law, "': max abs diff = ",
           signif(max_recurrence_err, 6), ".", call. = FALSE)
    }
  }

  list(
    law = sspec$law,
    rho = sspec$rho,
    loading_matrix = L_ind,
    covariance_matrix = Sigma_ind,
    max_loading_cov_diff = max_L_cov_err,
    max_supplied_cov_diff = max_supplied_cov_err,
    max_recurrence_diff = max_recurrence_err,
    verified = TRUE
  )
}

cp02_validate_base_rotating_draw <- function(base_draw) {
  if (!is.list(base_draw) ||
      !all(c("data", "oracle", "settings") %in% names(base_draw))) {
    stop("`base_draw` must be a list with `data`, `oracle`, and `settings`.",
         call. = FALSE)
  }
  st <- base_draw$settings
  ref_target <- sim_spec(x_trend = 0)
  if (!is.list(st) ||
      !identical(st$variant, "rotating") ||
      !identical(st$pattern_family, "rotating") ||
      !identical(st$innovation_family, "gaussian") ||
      !identical(st$sampling_scheme, "iid_target") ||
      !identical(st$cluster_mode, "unit_iid") ||
      isTRUE(st$has_period1_shift) ||
      isTRUE(st$post_null) ||
      !identical(st$target_spec, ref_target) ||
      !identical(as.numeric(st$time_scale), rep(1.0, 5L))) {
    stop("`base_draw` settings do not match the required rotating ",
         "x_trend = 0 Gaussian unit-iid iid-target base specification.",
         call. = FALSE)
  }

  n <- cp02_check_variant_n(st$n, cp02_baseline_variant_spec(
    variant = "rotating", base_spec = ref_target
  ))
  T_len <- 5L
  m <- n * T_len

  d_panel <- base_draw$data$panel
  d_units <- base_draw$data$units
  d_end <- base_draw$data$endpoints
  o_panel <- base_draw$oracle$panel
  o_units <- base_draw$oracle$units

  if (!is.data.frame(d_panel) || !is.data.frame(d_units) ||
      !is.data.frame(d_end) || !is.data.frame(o_panel) ||
      !is.data.frame(o_units) ||
      nrow(d_units) != n || nrow(d_end) != n || nrow(o_units) != n ||
      nrow(d_panel) != m || nrow(o_panel) != m) {
    stop("`base_draw` table dimensions are inconsistent with `n` = ", n, ".",
         call. = FALSE)
  }

  expected_id_unit <- seq_len(n)
  expected_id_panel <- rep(expected_id_unit, each = T_len)
  expected_period_panel <- rep(1:5, times = n)
  if (!identical(as.integer(d_units$id), expected_id_unit) ||
      !identical(as.integer(d_end$id), expected_id_unit) ||
      !identical(as.integer(o_units$id), expected_id_unit) ||
      !identical(as.integer(d_panel$id), expected_id_panel) ||
      !identical(as.integer(o_panel$id), expected_id_panel) ||
      !identical(as.numeric(d_panel$period), as.numeric(expected_period_panel)) ||
      !identical(as.numeric(o_panel$period), as.numeric(expected_period_panel))) {
    stop("`base_draw` unit/period ordering does not match canonical ",
         "row-by-unit (id = 1..n, period = 1..5) layout.", call. = FALSE)
  }

  valid_patterns <- c("1-2", "2-3", "3-4", "4-5")
  if (!all(d_units$pattern %in% valid_patterns) ||
      !all(d_units$pattern_type == "pair") ||
      !all(d_units$n_obs == 2L) ||
      !all(d_units$eligible_ge2) ||
      !identical(d_end$pattern, d_units$pattern) ||
      !identical(as.integer(d_end$date2 - d_end$date1), rep(1L, n)) ||
      !all(d_end$date1 %in% 1:4)) {
    stop("`base_draw` observation masks do not match rotating two-date ",
         "consecutive windows on periods 1:5.", call. = FALSE)
  }

  expected_obs <- (d_panel$period == d_end$date1[d_panel$id]) |
    (d_panel$period == d_end$date2[d_panel$id])
  if (!identical(as.logical(d_panel$observed), expected_obs)) {
    stop("`base_draw$data$panel$observed` is inconsistent with ",
         "`endpoints$date1` and `endpoints$date2`.", call. = FALSE)
  }

  if (!is.numeric(o_panel$eps_raw) || length(o_panel$eps_raw) != m ||
      any(is.na(o_panel$eps_raw)) || any(!is.finite(o_panel$eps_raw)) ||
      !identical(as.numeric(o_panel$eps_scale), rep(1.0, m)) ||
      max(abs(o_panel$eps - o_panel$eps_raw)) > 1e-12 ||
      !identical(as.numeric(o_panel$cluster_shock), rep(0.0, m)) ||
      !identical(as.numeric(o_panel$shift), rep(0.0, m))) {
    stop("`base_draw$oracle$panel` innovations/scales/shocks/shifts are ",
         "inconsistent with the unit-iid Gaussian no-shift rotating design.",
         call. = FALSE)
  }

  gk <- as.character(d_panel$G)
  expected_y0_non_error <- unname(ref_target$cohort_level[gk]) +
    unname(ref_target$time_effect[as.character(d_panel$period)]) +
    ref_target$x_trend * d_panel$x * (d_panel$period - 1.0) +
    ref_target$z_coef * d_panel$z +
    o_units$unit_effect[d_panel$id] +
    o_panel$cluster_shock +
    o_panel$shift
  if (max(abs((o_panel$y0 - o_panel$eps) - expected_y0_non_error)) > 1e-12) {
    stop("`base_draw$oracle$panel$y0` does not match the accepted Law A ",
         "(x_trend = 0) untreated potential outcome equation.", call. = FALSE)
  }

  treated <- d_panel$G > 0 & d_panel$period >= d_panel$G
  expected_effect <- numeric(m)
  expected_effect[treated] <- unname(ref_target$theta[gk[treated]]) +
    ref_target$delta * (d_panel$period[treated] - d_panel$G[treated]) +
    ref_target$kappa * d_panel$x[treated]
  if (max(abs(o_panel$effect - expected_effect)) > 1e-12 ||
      max(abs(o_panel$y_g - (o_panel$y0 + o_panel$effect))) > 1e-12) {
    stop("`base_draw$oracle$panel` treatment effects or `y_g` are ",
         "inconsistent with the accepted outcome equation.", call. = FALSE)
  }

  invisible(base_draw)
}

cp02_serial_transform_draw <- function(
    base_draw,
    law_or_spec = NULL,
    rho = NULL,
    ...) {
  sspec <- cp02_as_serial_spec(law_or_spec = law_or_spec, rho = rho, ...)
  cp02_validate_base_rotating_draw(base_draw)

  n <- base_draw$settings$n
  periods <- sspec$target_spec$periods
  T_len <- length(periods)

  eps_raw <- base_draw$oracle$panel$eps_raw
  eps_serial <- cp02_serial_transform_innovations(eps_raw, sspec)

  # Replace only outcome errors; preserve non-error potential outcome
  # components and treatment effects identically.
  y0_non_error <- base_draw$oracle$panel$y0 - base_draw$oracle$panel$eps
  y0_new <- y0_non_error + eps_serial
  effect <- base_draw$oracle$panel$effect
  y_g_new <- y0_new + effect

  out <- base_draw
  out$oracle$panel$eps <- eps_serial
  out$oracle$panel$y0 <- y0_new
  out$oracle$panel$y_g <- y_g_new
  out$oracle$truth <- cp02_serial_dgp_truth(sspec)

  observed <- out$data$panel$observed
  out$data$panel$y <- ifelse(observed, y_g_new, NA_real_)

  y_mat <- matrix(out$data$panel$y, nrow = n, ncol = T_len, byrow = TRUE)
  d1_pos <- match(out$data$endpoints$date1, periods)
  d2_pos <- match(out$data$endpoints$date2, periods)
  y1_val <- ifelse(is.na(d1_pos), NA_real_, y_mat[cbind(seq_len(n), d1_pos)])
  y2_val <- ifelse(is.na(d2_pos), NA_real_, y_mat[cbind(seq_len(n), d2_pos)])
  out$data$endpoints$y_date1 <- y1_val
  out$data$endpoints$y_date2 <- y2_val
  out$data$endpoints$dy <- ifelse(
    !is.na(y1_val) & !is.na(y2_val),
    y2_val - y1_val,
    NA_real_
  )

  settings <- out$settings
  settings$family <- sspec$family
  settings$law <- sspec$law
  settings$design_table_law_index <- sspec$design_table_law_index
  settings$status <- sspec$status
  settings$rho <- sspec$rho
  settings$innovation_family <- sspec$innovation_family
  settings$innovation_variance <- sspec$innovation_variance
  settings$sd_innovation <- sspec$sd_innovation
  settings$start_convention <- sspec$start_convention
  settings$start_distribution <- sspec$start_distribution
  settings$recurrence_equation <- sspec$recurrence_equation
  settings$stationary <- sspec$stationary
  settings$process_type <- sspec$process_type
  settings$eps_scale_semantics <- sspec$eps_scale_semantics
  settings$legacy_provenance <- sspec$legacy_provenance
  settings$campaign_fit_spec <- sspec$campaign_fit_spec
  settings$loading_matrix <- sspec$loading_matrix
  settings$covariance_matrix <- sspec$covariance_matrix
  settings$stream_note <- sspec$stream_note
  if (is.list(settings$rng)) {
    settings$rng$draw_order <- sspec$draw_order
  }
  out$settings <- settings
  out
}

cp02_serial_dgp_draw <- function(
    law_or_spec = NULL,
    n,
    seed,
    n_clusters = NULL,
    rho = NULL,
    ...) {
  sspec <- cp02_as_serial_spec(law_or_spec = law_or_spec, rho = rho, ...)
  n <- cp02_check_variant_n(n, sspec$base_variant_spec)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("`seed` must be a single whole number in integer range.",
         call. = FALSE)
  }
  if (!is.null(n_clusters)) {
    if (!is.numeric(n_clusters) || length(n_clusters) != 1L ||
        !is.finite(n_clusters) || n_clusters != round(n_clusters) ||
        n_clusters < 1) {
      stop("`n_clusters` must be a positive whole number.", call. = FALSE)
    }
  }
  # Single seeded base constructor call (executes inside one with_seed(seed, ...)
  # block and restores the caller's RNG state on exit), followed by the pure
  # deterministic AR outcome-error transformation.
  base_draw <- cp02_baseline_variant_draw(
    vspec = sspec$base_variant_spec,
    n = n,
    seed = seed,
    n_clusters = n_clusters
  )
  cp02_serial_transform_draw(base_draw = base_draw, law_or_spec = sspec)
}
