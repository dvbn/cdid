# Development simulation support for cdid (not package API).
# Requires seed.R (with_seed) and dgp.R (sim_spec, sim_truth).
# Does not clone or modify dgp.R, seed.R, or cp02_continuous_dgp.R.
#
# Scoped CP02 baseline-validation DGP variants (proposal/development only,
# not production cell approval, sample-size/role approval, or inference
# certification). Mathematical and engineering contract:
# docs/theory/cp02_baseline_variants_dgp.tex.
#
# Reuses sim_spec() and sim_truth() for the binary-(X, G) target population and
# with_seed() for deterministic RNG isolation. Covers the missing baseline
# validation variants in plan/CP02_FINAL_MC_DESIGN.md:
#   1. balanced                        : 5-period balanced binary-X Gaussian
#   2. rotating                        : 4 consecutive two-period windows
#   3. general_pairs                   : 10 two-date pairs + singleton{3} + none
#   4. fixed_horizon_bridge            : {1,2}/{1,4} (0.4 each) + {3}/none (0.1
#                                        each) with explicit cohort-specific
#                                        period-1 shifts on BOTH potential paths
#   5. skewed_innovations              : centered-exponential period errors with
#                                        time-varying scale
#   6. nonconstant_independent_weights : bounded weights independent of (G,X,Y)
#   7. unequal_prob_iid                : tilted sampling P_s(x,g) on finite
#                                        binary(X,G) cells with w = 1/pi(x,g)
#   8. cluster_dependence              : fixed 10 units/cluster with independent
#                                        common cluster-by-period shocks
#   9. cluster_zero_eligible           : fixed 10 units/cluster + independent
#                                        cluster-level observation flag creating
#                                        singleton-only / zero-eligible clusters
#  10. true_post_null                  : post-effect parameters zero in sim_spec

cp02_baseline_variants_catalog <- function() {
  data.frame(
    variant = c(
      "balanced",
      "rotating",
      "general_pairs",
      "fixed_horizon_bridge",
      "skewed_innovations",
      "nonconstant_independent_weights",
      "unequal_prob_iid",
      "cluster_dependence",
      "cluster_zero_eligible",
      "true_post_null"
    ),
    pattern_family = c(
      "balanced",
      "rotating",
      "general_pairs",
      "fixed_horizon_bridge",
      "general_pairs",
      "general_pairs",
      "general_pairs",
      "general_pairs",
      "general_pairs",
      "general_pairs"
    ),
    innovation_family = c(
      rep("gaussian", 4L),
      "centered_exponential",
      rep("gaussian", 5L)
    ),
    sampling_scheme = c(
      rep("iid_target", 5L),
      "independent_weights",
      "unequal_prob_iid",
      rep("iid_target", 3L)
    ),
    cluster_mode = c(
      rep("unit_iid", 7L),
      "cluster_shock",
      "cluster_zero_eligible",
      "unit_iid"
    ),
    period1_shift = c(
      FALSE, FALSE, FALSE, TRUE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    ),
    post_null = c(rep(FALSE, 9L), TRUE),
    stringsAsFactors = FALSE
  )
}

cp02_baseline_patterns <- function(
    pattern_family = c("general_pairs", "balanced", "rotating",
                       "fixed_horizon_bridge")) {
  pattern_family <- match.arg(pattern_family)
  if (pattern_family == "balanced") {
    return(data.frame(
      pattern = "1-2-3-4-5",
      type = "full",
      date1 = 1L,
      date2 = 5L,
      n_obs = 5L,
      eligible_ge2 = TRUE,
      slots = 1L,
      prob = 1.0,
      stringsAsFactors = FALSE
    ))
  }
  if (pattern_family == "rotating") {
    return(data.frame(
      pattern = c("1-2", "2-3", "3-4", "4-5"),
      type = rep("pair", 4L),
      date1 = 1:4,
      date2 = 2:5,
      n_obs = rep(2L, 4L),
      eligible_ge2 = rep(TRUE, 4L),
      slots = rep(1L, 4L),
      prob = rep(0.25, 4L),
      stringsAsFactors = FALSE
    ))
  }
  if (pattern_family == "general_pairs") {
    pairs <- utils::combn(5L, 2L)
    return(data.frame(
      pattern = c(paste(pairs[1L, ], pairs[2L, ], sep = "-"), "3", "none"),
      type = c(rep("pair", 10L), "singleton", "none"),
      date1 = c(pairs[1L, ], 3L, NA_integer_),
      date2 = c(pairs[2L, ], NA_integer_, NA_integer_),
      n_obs = c(rep(2L, 10L), 1L, 0L),
      eligible_ge2 = c(rep(TRUE, 10L), FALSE, FALSE),
      slots = c(rep(4L, 10L), 5L, 5L),
      prob = c(rep(0.08, 10L), 0.10, 0.10),
      stringsAsFactors = FALSE
    ))
  }
  # fixed_horizon_bridge: {1,2} w.p. 0.4, {1,4} w.p. 0.4, {3} w.p. 0.1, none 0.1
  data.frame(
    pattern = c("1-2", "1-4", "3", "none"),
    type = c("pair", "pair", "singleton", "none"),
    date1 = c(1L, 1L, 3L, NA_integer_),
    date2 = c(2L, 4L, NA_integer_, NA_integer_),
    n_obs = c(2L, 2L, 1L, 0L),
    eligible_ge2 = c(TRUE, TRUE, FALSE, FALSE),
    slots = c(4L, 4L, 1L, 1L),
    prob = c(0.40, 0.40, 0.10, 0.10),
    stringsAsFactors = FALSE
  )
}

# Exact slot-based inverse-CDF map from U(0, 1) to pattern row index.
cp02_baseline_pattern_index <- function(
    v,
    pattern_family = c("general_pairs", "balanced", "rotating",
                       "fixed_horizon_bridge")) {
  pattern_family <- match.arg(pattern_family)
  stopifnot(is.numeric(v), all(is.finite(v)), all(v >= 0), all(v <= 1))
  pat <- cp02_baseline_patterns(pattern_family)
  total_slots <- sum(pat$slots)
  slot <- pmin(floor(total_slots * v), total_slots - 1L)
  cum_slots <- cumsum(pat$slots)
  1L + rowSums(outer(slot, cum_slots[-length(cum_slots)], ">="))
}

# Exact inverse-CDF map for standardized centered rate-1 exponential E - 1.
cp02_centered_exp_quantile <- function(u) {
  stopifnot(is.numeric(u), all(is.finite(u)), all(u >= 0), all(u < 1))
  (-log1p(-u)) - 1
}

cp02_default_inclusion_prob <- function(cohorts = c(3, 4, 5)) {
  codes <- as.character(c(0, cohorts))
  stopifnot(length(codes) == 4L)
  # Dyadic inverses w = 1/pi in {1.25, 2.00, 2.50}; canonical E_target[pi] = 0.50.
  pi_mat <- rbind(
    x0 = c(0.40, 0.50, 0.80, 0.40),
    x1 = c(0.80, 0.40, 0.50, 0.50)
  )
  colnames(pi_mat) <- codes
  pi_mat
}

cp02_default_period1_shift <- function(cohorts = c(3, 4, 5), active = FALSE) {
  codes <- c(0, cohorts)
  stopifnot(length(codes) == 4L)
  if (!active) {
    return(data.frame(
      G = codes,
      intercept = rep(0, 4L),
      slope_x = rep(0, 4L)
    ))
  }
  data.frame(
    G = codes,
    intercept = c(0, 2, -1, 1),
    slope_x = c(0, 0.7, 0.2, -0.4)
  )
}

cp02_baseline_variant_spec <- function(
    variant = c(
      "general_pairs",
      "balanced",
      "rotating",
      "fixed_horizon_bridge",
      "skewed_innovations",
      "nonconstant_independent_weights",
      "unequal_prob_iid",
      "cluster_dependence",
      "cluster_zero_eligible",
      "true_post_null"
    ),
    pattern_family = NULL,
    innovation_family = NULL,
    time_scale = NULL,
    sampling_scheme = NULL,
    inclusion_prob = NULL,
    indep_weight_bounds = c(0.5, 2.0),
    cluster_mode = NULL,
    cluster_size = 10L,
    sd_cluster = NULL,
    cluster_obs_prob = NULL,
    period1_shift = NULL,
    post_null = NULL,
    base_spec = NULL) {
  variant <- match.arg(variant)
  cat_row <- cp02_baseline_variants_catalog()
  cat_row <- cat_row[cat_row$variant == variant, , drop = FALSE]

  if (is.null(pattern_family)) pattern_family <- cat_row$pattern_family
  pattern_family <- match.arg(
    pattern_family,
    c("general_pairs", "balanced", "rotating", "fixed_horizon_bridge")
  )
  if (is.null(innovation_family)) innovation_family <- cat_row$innovation_family
  innovation_family <- match.arg(
    innovation_family,
    c("gaussian", "centered_exponential")
  )
  if (is.null(sampling_scheme)) sampling_scheme <- cat_row$sampling_scheme
  sampling_scheme <- match.arg(
    sampling_scheme,
    c("iid_target", "unequal_prob_iid", "independent_weights")
  )
  if (is.null(cluster_mode)) cluster_mode <- cat_row$cluster_mode
  cluster_mode <- match.arg(
    cluster_mode,
    c("unit_iid", "cluster_shock", "cluster_zero_eligible")
  )
  if (is.null(post_null)) post_null <- isTRUE(cat_row$post_null)
  stopifnot(is.logical(post_null), length(post_null) == 1L, !is.na(post_null))

  if (is.null(base_spec)) {
    target_spec <- if (post_null) {
      sim_spec(theta = c(0, 0, 0), delta = 0, kappa = 0)
    } else {
      sim_spec()
    }
  } else {
    target_spec <- if (post_null) {
      sim_spec(
        periods = base_spec$periods,
        cohorts = base_spec$cohorts,
        p_x = base_spec$p_x,
        cohort_prob = base_spec$cohort_prob,
        theta = rep(0, length(base_spec$cohorts)),
        delta = 0,
        kappa = 0,
        cohort_level = unname(base_spec$cohort_level),
        time_effect = unname(base_spec$time_effect),
        x_trend = base_spec$x_trend,
        z_coef = base_spec$z_coef,
        sd_unit = base_spec$sd_unit,
        sd_eps = base_spec$sd_eps,
        n_clusters = base_spec$n_clusters
      )
    } else {
      base_spec
    }
  }

  periods <- target_spec$periods
  cohorts <- target_spec$cohorts
  codes <- c(0, cohorts)
  stopifnot(identical(as.numeric(periods), as.numeric(1:5)),
            identical(as.numeric(cohorts), c(3, 4, 5)))

  if (is.null(time_scale)) {
    time_scale <- if (innovation_family == "centered_exponential") {
      (2 + periods) / 4
    } else {
      rep(1, length(periods))
    }
  }
  time_scale <- as.numeric(time_scale)
  stopifnot(length(time_scale) == length(periods),
            all(is.finite(time_scale)), all(time_scale > 0))
  names(time_scale) <- as.character(periods)

  if (is.null(period1_shift)) {
    shift_df <- cp02_default_period1_shift(cohorts, active = cat_row$period1_shift)
  } else if (is.logical(period1_shift) && length(period1_shift) == 1L) {
    shift_df <- cp02_default_period1_shift(cohorts, active = period1_shift)
  } else {
    stopifnot(
      is.data.frame(period1_shift),
      identical(as.numeric(period1_shift$G), as.numeric(codes)),
      is.numeric(period1_shift$intercept),
      is.numeric(period1_shift$slope_x),
      all(is.finite(period1_shift$intercept)),
      all(is.finite(period1_shift$slope_x))
    )
    shift_df <- period1_shift
  }
  has_period1_shift <- any(shift_df$intercept != 0 | shift_df$slope_x != 0)

  if (is.null(inclusion_prob)) {
    if (sampling_scheme == "unequal_prob_iid") {
      inclusion_prob <- cp02_default_inclusion_prob(cohorts)
    } else {
      inclusion_prob <- matrix(1, nrow = 2L, ncol = length(codes),
                               dimnames = list(c("x0", "x1"), codes))
    }
  } else {
    stopifnot(is.matrix(inclusion_prob),
              identical(dim(inclusion_prob), c(2L, length(codes))),
              all(is.finite(inclusion_prob)),
              all(inclusion_prob > 0),
              all(inclusion_prob <= 1))
    dimnames(inclusion_prob) <- list(c("x0", "x1"), codes)
  }

  indep_weight_bounds <- as.numeric(indep_weight_bounds)
  stopifnot(length(indep_weight_bounds) == 2L,
            all(is.finite(indep_weight_bounds)),
            indep_weight_bounds[1L] > 0,
            indep_weight_bounds[2L] > indep_weight_bounds[1L])

  cluster_size <- as.integer(cluster_size)
  stopifnot(length(cluster_size) == 1L, cluster_size == 10L)

  if (is.null(sd_cluster)) {
    sd_cluster <- if (cluster_mode == "unit_iid") 0 else 0.75
  }
  stopifnot(is.numeric(sd_cluster), length(sd_cluster) == 1L,
            is.finite(sd_cluster), sd_cluster >= 0)
  if (cluster_mode != "unit_iid" && sd_cluster <= 0) {
    stop("Genuine cluster dependence requires `sd_cluster > 0`.", call. = FALSE)
  }

  if (is.null(cluster_obs_prob)) {
    cluster_obs_prob <- if (cluster_mode == "cluster_zero_eligible") 0.70 else 1.0
  }
  stopifnot(is.numeric(cluster_obs_prob), length(cluster_obs_prob) == 1L,
            is.finite(cluster_obs_prob), cluster_obs_prob > 0,
            cluster_obs_prob <= 1)
  if (cluster_mode == "cluster_zero_eligible" && cluster_obs_prob >= 1) {
    stop("`cluster_zero_eligible` requires `cluster_obs_prob < 1`.",
         call. = FALSE)
  }

  # Finite-logit validation requires positive cell odds. Refuse a sampling
  # declaration that would silently change the target without inverse weights.
  if (any(target_spec$cohort_prob <= 0)) {
    stop("Binary-X finite-logit validation requires positive (X,G) cells.",
         call. = FALSE)
  }
  if (sampling_scheme != "unequal_prob_iid" && any(inclusion_prob != 1)) {
    stop("Nontrivial `inclusion_prob` requires `sampling_scheme = unequal_prob_iid`.",
         call. = FALSE)
  }
  patterns <- cp02_baseline_patterns(pattern_family)
  if (cluster_mode == "cluster_zero_eligible" &&
      !any(patterns$type == "singleton")) {
    # This row is a flag-induced mask, never an additional base-pattern draw.
    patterns <- rbind(patterns, data.frame(
      pattern = "3", type = "singleton", date1 = 3L, date2 = NA_integer_,
      n_obs = 1L, eligible_ge2 = FALSE, slots = 0L, prob = 0,
      stringsAsFactors = FALSE
    ))
  }
  patterns$prob_unconditional <- patterns$prob
  if (cluster_mode == "cluster_zero_eligible") {
    collapsing <- patterns$type %in% c("pair", "full")
    singleton <- which(patterns$type == "singleton")
    patterns$prob_unconditional[collapsing] <-
      cluster_obs_prob * patterns$prob[collapsing]
    patterns$prob_unconditional[singleton] <- patterns$prob[singleton] +
      (1 - cluster_obs_prob) * sum(patterns$prob[collapsing])
  }

  list(
    family = "cp02_baseline_validation_variants",
    variant = variant,
    status = paste(
      "Development validation DGP proposal only; not production cell",
      "approval, sample-size/role approval, public default, or",
      "estimator inference certification."
    ),
    target_spec = target_spec,
    pattern_family = pattern_family,
    patterns = patterns,
    pattern_probability_note = paste(
      "patterns$prob is the independent base-pattern law;",
      "patterns$prob_unconditional includes the independent cluster flag.",
      "Inactive-cluster full/pair masks become singleton {3}; none stays none."
    ),
    innovation_family = innovation_family,
    time_scale = time_scale,
    period1_shift = shift_df,
    has_period1_shift = has_period1_shift,
    period1_shift_both_paths = TRUE,
    sampling_scheme = sampling_scheme,
    inclusion_prob = inclusion_prob,
    indep_weight_bounds = stats::setNames(indep_weight_bounds,
                                          c("lower", "upper")),
    cluster_mode = cluster_mode,
    cluster_size = cluster_size,
    sd_cluster = sd_cluster,
    cluster_obs_prob = cluster_obs_prob,
    post_null = post_null,
    default_fitting_eligibility = "n_obs >= 2L (singletons and none preserved in frame/panel)",
    assumptions_and_scope = list(
      truth_invariant = paste(
        "Target population ATT(g,t) and population aggregation truths are",
        "defined solely by sim_truth(target_spec). Realized, selected, or",
        "tilted sample counts never define truth."
      ),
      saturated_binary_x_logit = paste(
        "Binary X in {0,1} with intercept is saturated, so ~ x and ~ x + z",
        "are exact logits for all NT and NYT control pools (including early",
        "NYT pools). By contrast, cp02_continuous_dgp.R remains unchanged",
        "and guarantees correct logits only for its declared future-control",
        "operators, not early NYT pools containing G=3."
      ),
      bridge_semantics = paste(
        "In fixed_horizon_bridge, cohort-specific period-1 shifts s_G(X) are",
        "added to BOTH potential outcome paths Y_it(0) and Y_it(g). Causal",
        "pre-treatment effects remain 0, and conditional PT holds from",
        "baseline t >= 2 onward. Pre-period DiD reporting estimands involving",
        "period 1 are nonzero and operator/baseline-specific."
      ),
      weighted_sampling_scope = paste(
        "Unequal-probability iid sampling tilts finite (X,G) cells by",
        "bounded positive inclusion probabilities pi(x,g) with inverse",
        "weights w(x,g) = 1/pi(x,g), satisfying E_s[w f]/E_s[w] = E_target[f]",
        "exactly on finite cells. This does not claim arbitrary complex",
        "survey-design or endogenous selection-correction validity."
      ),
      cluster_scope = paste(
        "Cluster variants draw independent clusters of fixed 10 units with",
        "shared cluster-by-period shocks eta_{c,t} ~ N(0, sd_cluster^2) and",
        "an independent cluster observation flag B_c ~ Bernoulli(cluster_obs_prob).",
        "Asymptotic validity requires growing cluster count C; small-C draws",
        "are mechanical verification fixtures only."
      ),
      deferred_variants = paste(
        "Anticipation and irregular-calendar variants are deferred to their",
        "own focused calendar/anticipation check; final feature completeness",
        "is not claimed here."
      )
    )
  )
}

cp02_baseline_tilt_table <- function(vspec) {
  spec <- vspec$target_spec
  codes <- c(0, spec$cohorts)
  p_target_mat <- rbind(
    x0 = (1 - spec$p_x) * spec$cohort_prob[1L, ],
    x1 = spec$p_x * spec$cohort_prob[2L, ]
  )
  colnames(p_target_mat) <- as.character(codes)
  pi_mat <- vspec$inclusion_prob
  mean_pi_target <- sum(p_target_mat * pi_mat)
  p_sample_mat <- (pi_mat * p_target_mat) / mean_pi_target
  w_mat <- 1 / pi_mat
  mean_w_sample <- sum(p_sample_mat * w_mat)
  norm_w_mat <- w_mat / mean_w_sample

  cells <- expand.grid(
    x = c(0L, 1L),
    G = codes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  idx <- cbind(cells$x + 1L, match(cells$G, codes))
  cells$p_target <- p_target_mat[idx]
  cells$pi_include <- pi_mat[idx]
  cells$p_sample <- p_sample_mat[idx]
  cells$inv_weight <- w_mat[idx]
  cells$norm_weight <- norm_w_mat[idx]
  # Exact term-by-term contribution to E_s[w f] / E_s[w]
  cells$recovered_p_target <- (cells$p_sample * cells$inv_weight) / mean_w_sample
  cells$abs_recovery_diff <- abs(cells$recovered_p_target - cells$p_target)

  pi_g_target <- colSums(p_target_mat)
  mean_x_target <- p_target_mat[2L, ] / pi_g_target
  pi_g_sample_unweighted <- colSums(p_sample_mat)
  mean_x_sample_unweighted <- p_sample_mat[2L, ] / pi_g_sample_unweighted

  pi_g_recovered <- vapply(seq_along(codes), function(j) {
    sum(p_sample_mat[, j] * w_mat[, j]) / mean_w_sample
  }, numeric(1))
  mean_x_recovered <- vapply(seq_along(codes), function(j) {
    (p_sample_mat[2L, j] * w_mat[2L, j]) / sum(p_sample_mat[, j] * w_mat[, j])
  }, numeric(1))

  cohort_comparison <- data.frame(
    G = codes,
    pi_target = unname(pi_g_target),
    pi_sample_unweighted = unname(pi_g_sample_unweighted),
    pi_sample_weighted = unname(pi_g_recovered),
    mean_x_target = unname(mean_x_target),
    mean_x_sample_unweighted = unname(mean_x_sample_unweighted),
    mean_x_sample_weighted = unname(mean_x_recovered),
    stringsAsFactors = FALSE
  )

  list(
    cells = cells,
    p_target_matrix = p_target_mat,
    inclusion_prob_matrix = pi_mat,
    p_sample_matrix = p_sample_mat,
    inv_weight_matrix = w_mat,
    mean_pi_target = mean_pi_target,
    mean_w_sample = mean_w_sample,
    inverse_weight_table_note = paste(
      "inv_weight, norm_weight and mean_w_sample describe inverse inclusion",
      "weights 1/pi, not the independently drawn random W scheme."
    ),
    mean_actual_sample_weight = if (vspec$sampling_scheme == "independent_weights") {
      mean(vspec$indep_weight_bounds)
    } else {
      mean_w_sample
    },
    identity_identity_error = max(cells$abs_recovery_diff),
    cohort_comparison = cohort_comparison
  )
}

cp02_baseline_propensity_table <- function(vspec) {
  spec <- vspec$target_spec
  codes <- c(0, spec$cohorts)
  p_target_mat <- rbind(
    x0 = (1 - spec$p_x) * spec$cohort_prob[1L, ],
    x1 = spec$p_x * spec$cohort_prob[2L, ]
  )
  colnames(p_target_mat) <- as.character(codes)

  tilt <- cp02_baseline_tilt_table(vspec)
  laws <- list(
    target = p_target_mat,
    sample_unweighted = tilt$p_sample_matrix,
    sample_weighted = p_target_mat
  )
  rows <- list()
  for (law in names(laws)) {
    p_law <- laws[[law]]
    for (g in spec$cohorts) {
      for (tt in spec$periods) {
        for (ctl in c("NT", "NYT")) {
          pool <- if (ctl == "NT") {
            0
          } else {
            c(0, setdiff(spec$cohorts[spec$cohorts > tt], g))
          }
          pool_chr <- as.character(pool)
          p0_g <- p_law["x0", as.character(g)]
          p1_g <- p_law["x1", as.character(g)]
          p0_c <- sum(p_law["x0", pool_chr])
          p1_c <- sum(p_law["x1", pool_chr])
          alpha <- log(p0_g / p0_c)
          beta_x <- log(p1_g / p1_c) - alpha
          rows[[length(rows) + 1L]] <- data.frame(
            population_law = law,
            pool_rule = "cohorts untreated at t, excluding target cohort",
            group = g,
            t = tt,
            post = tt >= g,
            control = ctl,
            pool = paste(pool, collapse = ","),
            pool_size = length(pool),
            contains_early_g3 = 3 %in% pool,
            intercept = unname(alpha),
            slope_x = unname(beta_x),
            slope_z = 0,
            prob_x0 = unname(p0_g / (p0_g + p0_c)),
            prob_x1 = unname(p1_g / (p1_g + p1_c)),
            correct_x_only = TRUE,
            correct_x_plus_z = TRUE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

cp02_baseline_pre_reporting_table <- function(vspec) {
  spec <- vspec$target_spec
  codes <- c(0, spec$cohorts)
  p_target_mat <- rbind(
    x0 = (1 - spec$p_x) * spec$cohort_prob[1L, ],
    x1 = spec$p_x * spec$cohort_prob[2L, ]
  )
  colnames(p_target_mat) <- as.character(codes)
  pi_g <- colSums(p_target_mat)
  px_given_g <- sweep(p_target_mat, 2L, pi_g, "/")

  shift_df <- vspec$period1_shift
  s_mat <- rbind(
    x0 = shift_df$intercept,
    x1 = shift_df$intercept + shift_df$slope_x
  )
  colnames(s_mat) <- as.character(shift_df$G)

  p_sample_mat <- cp02_baseline_tilt_table(vspec)$p_sample_matrix
  sample_did <- function(g, pool, sign) {
    sum(vapply(c("x0", "x1"), function(x) {
      px_g <- p_sample_mat[x, as.character(g)] /
        sum(p_sample_mat[, as.character(g)])
      pool_shift <- sum(p_sample_mat[x, as.character(pool)] *
                          s_mat[x, as.character(pool)]) /
        sum(p_sample_mat[x, as.character(pool)])
      px_g * (s_mat[x, as.character(g)] - pool_shift)
    }, numeric(1))) * sign
  }

  pat <- vspec$patterns
  pair_rows <- pat[pat$type %in% c("full", "pair"), , drop = FALSE]
  has_pair <- function(b, tt) {
    if (vspec$pattern_family == "balanced") return(TRUE)
    lo <- min(b, tt)
    hi <- max(b, tt)
    any(pair_rows$date1 == lo & pair_rows$date2 == hi)
  }

  out <- list()
  for (g in spec$cohorts) {
    pre_ts <- spec$periods[spec$periods < g]
    for (base_mode in c("varying", "universal")) {
      for (tt in pre_ts) {
        b <- if (base_mode == "varying") tt - 1 else g - 1
        if (!(b %in% spec$periods) || b == tt) next
        for (ctl in c("NT", "NYT")) {
          max_t <- max(b, tt)
          pool <- if (ctl == "NT") {
            0
          } else {
            c(0, setdiff(spec$cohorts[spec$cohorts > max_t], g))
          }
          pool_chr <- as.character(pool)
          w_pool_x0 <- p_target_mat["x0", pool_chr] /
            sum(p_target_mat["x0", pool_chr])
          w_pool_x1 <- p_target_mat["x1", pool_chr] /
            sum(p_target_mat["x1", pool_chr])
          s_pool_x0 <- sum(w_pool_x0 * s_mat["x0", pool_chr])
          s_pool_x1 <- sum(w_pool_x1 * s_mat["x1", pool_chr])
          diff_s_x0 <- s_mat["x0", as.character(g)] - s_pool_x0
          diff_s_x1 <- s_mat["x1", as.character(g)] - s_pool_x1
          mean_diff_s <- px_given_g["x0", as.character(g)] * diff_s_x0 +
            px_given_g["x1", as.character(g)] * diff_s_x1
          period1_sign <- as.numeric(tt == spec$periods[1L]) -
            as.numeric(b == spec$periods[1L])
          did_estimand <- period1_sign * mean_diff_s
          out[[length(out) + 1L]] <- data.frame(
            group = g,
            t = tt,
            base_period_mode = base_mode,
            base_period = b,
            control = ctl,
            pool = paste(pool, collapse = ","),
            causal_att = 0,
            involves_period1 = (b == spec$periods[1L] || tt == spec$periods[1L]),
            truth_law = "target",
            reporting_did_truth = unname(did_estimand),
            reporting_did_truth_sample_unweighted = sample_did(g, pool, period1_sign),
            reporting_did_truth_sample_weighted = unname(did_estimand),
            prob_x0_target = unname(p_target_mat["x0", as.character(g)] /
              sum(p_target_mat["x0", as.character(c(g, pool))])),
            prob_x1_target = unname(p_target_mat["x1", as.character(g)] /
              sum(p_target_mat["x1", as.character(c(g, pool))])),
            prob_x0_sample_unweighted = unname(p_sample_mat["x0", as.character(g)] /
              sum(p_sample_mat["x0", as.character(c(g, pool))])),
            prob_x1_sample_unweighted = unname(p_sample_mat["x1", as.character(g)] /
              sum(p_sample_mat["x1", as.character(c(g, pool))])),
            pair_observed_in_pattern = has_pair(b, tt),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, out)
}

cp02_baseline_bridge_d13_table <- function(vspec) {
  spec <- vspec$target_spec
  base_truth <- sim_truth(spec)
  codes <- c(0, spec$cohorts)
  p_target_mat <- rbind(
    x0 = (1 - spec$p_x) * spec$cohort_prob[1L, ],
    x1 = spec$p_x * spec$cohort_prob[2L, ]
  )
  colnames(p_target_mat) <- as.character(codes)
  pi_g <- colSums(p_target_mat)
  px_given_g <- sweep(p_target_mat, 2L, pi_g, "/")

  shift_df <- vspec$period1_shift
  s_mat <- rbind(
    x0 = shift_df$intercept,
    x1 = shift_df$intercept + shift_df$slope_x
  )
  colnames(s_mat) <- as.character(shift_df$G)

  att_34 <- base_truth$att_gt$att[
    base_truth$att_gt$group == 3 & base_truth$att_gt$t == 4
  ]

  rows <- list()
  for (ctl in c("NT", "NYT")) {
    # Fixed-horizon control pool for (g=3, t=4) must be untreated at t=4:
    pool <- if (ctl == "NT") 0 else c(0, 5)
    pool_chr <- as.character(pool)
    w_pool_x0 <- p_target_mat["x0", pool_chr] / sum(p_target_mat["x0", pool_chr])
    w_pool_x1 <- p_target_mat["x1", pool_chr] / sum(p_target_mat["x1", pool_chr])
    s_pool_x0 <- sum(w_pool_x0 * s_mat["x0", pool_chr])
    s_pool_x1 <- sum(w_pool_x1 * s_mat["x1", pool_chr])
    shift_bias <- px_given_g["x0", "3"] * (s_mat["x0", "3"] - s_pool_x0) +
      px_given_g["x1", "3"] * (s_mat["x1", "3"] - s_pool_x1)
    link_14_did <- att_34 - shift_bias
    link_12_did <- -shift_bias
    bridged_24_did <- link_14_did - link_12_did
    rows[[length(rows) + 1L]] <- data.frame(
      group = 3,
      t = 4,
      baseline = 2,
      anchor = 1,
      control = ctl,
      pool = paste(pool, collapse = ","),
      link_14_did = unname(link_14_did),
      link_12_did = unname(link_12_did),
      period1_shift_bias_cancelled = unname(shift_bias),
      bridged_att = unname(bridged_24_did),
      target_att = unname(att_34),
      abs_cancellation_diff = abs(unname(bridged_24_did - att_34)),
      stringsAsFactors = FALSE
    )
  }

  post_cells <- base_truth$att_gt[base_truth$att_gt$post, c("group", "t", "att")]
  post_cells$baseline <- post_cells$group - 1
  post_cells$recoverable_in_bridge <- post_cells$group == 3 & post_cells$t == 4
  post_cells$reason <- ifelse(
    post_cells$recoverable_in_bridge,
    "Recoverable via D13 fixed-horizon cancellation (1->4) - (1->2) at common anchor 1",
    ifelse(
      post_cells$t %in% c(3, 5),
      "Target period has no observed pair links in {1-2, 1-4}; t=3 is singleton-only and t=5 is unobserved",
      "Canonical baseline r_g = g - 1 = 3 has no observed pair links in {1-2, 1-4}"
    )
  )
  rownames(post_cells) <- NULL

  list(
    att34_cancellation = do.call(rbind, rows),
    post_cell_recoverability = post_cells
  )
}

cp02_baseline_noise_moments <- function(vspec) {
  spec <- vspec$target_spec
  periods <- spec$periods
  sigma_t <- spec$sd_eps * unname(vspec$time_scale)
  var_cl <- if (vspec$cluster_mode == "unit_iid") 0 else vspec$sd_cluster^2

  period_df <- data.frame(
    period = periods,
    sd_idiosyncratic = sigma_t,
    var_idiosyncratic = sigma_t^2,
    var_cluster_shock = rep(var_cl, length(periods)),
    var_total_noise = sigma_t^2 + var_cl,
    cond_mean = rep(0, length(periods)),
    third_central_idiosyncratic = if (vspec$innovation_family == "centered_exponential") {
      2 * sigma_t^3
    } else {
      rep(0, length(periods))
    },
    fourth_central_idiosyncratic = if (vspec$innovation_family == "centered_exponential") {
      9 * sigma_t^4
    } else {
      3 * sigma_t^4
    },
    within_cluster_cov_same_period = rep(var_cl, length(periods))
  )

  pairs <- utils::combn(length(periods), 2L)
  s_idx <- pairs[1L, ]
  t_idx <- pairs[2L, ]
  sig_s <- sigma_t[s_idx]
  sig_t <- sigma_t[t_idx]
  var_diff_idio <- sig_t^2 + sig_s^2
  var_diff_total <- var_diff_idio + 2 * var_cl
  third_diff_idio <- if (vspec$innovation_family == "centered_exponential") {
    2 * (sig_t^3 - sig_s^3)
  } else {
    rep(0, ncol(pairs))
  }
  fourth_diff_idio <- if (vspec$innovation_family == "centered_exponential") {
    9 * (sig_t^4 + sig_s^4) + 6 * sig_t^2 * sig_s^2
  } else {
    3 * var_diff_idio^2
  }
  fourth_diff_total <- fourth_diff_idio +
    6 * var_diff_idio * (2 * var_cl) + 3 * (2 * var_cl)^2

  pair_df <- data.frame(
    date1 = periods[s_idx],
    date2 = periods[t_idx],
    cond_mean_diff = rep(0, ncol(pairs)),
    var_diff_idiosyncratic = var_diff_idio,
    var_diff_cluster = rep(2 * var_cl, ncol(pairs)),
    var_diff_total = var_diff_total,
    third_central_diff = third_diff_idio,
    skewness_diff_idiosyncratic = ifelse(
      var_diff_idio > 0,
      third_diff_idio / (var_diff_idio^(3 / 2)),
      0
    ),
    fourth_central_diff_total = fourth_diff_total,
    within_cluster_cov_diff = rep(2 * var_cl, ncol(pairs))
  )

  list(
    innovation_family = vspec$innovation_family,
    finite_moments_all_orders = TRUE,
    period_moments = period_df,
    pair_difference_moments = pair_df
  )
}

cp02_baseline_variant_truth <- function(vspec) {
  # Reuse sim_truth(target_spec) directly for population cohort shares/means,
  # post-treatment ATT(g,t), and all four population aggregates.
  base_truth <- sim_truth(vspec$target_spec)
  cells <- base_truth$att_gt
  cells$causal_att <- cells$att
  if (vspec$has_period1_shift) {
    cells$reporting_did_truth <- ifelse(cells$post, cells$att, NA_real_)
  } else {
    cells$reporting_did_truth <- cells$att
  }
  if (vspec$pattern_family == "fixed_horizon_bridge") {
    cells$recoverable_under_pattern <- ifelse(
      cells$post,
      cells$group == 3 & cells$t == 4,
      NA
    )
  } else {
    cells$recoverable_under_pattern <- ifelse(cells$post, TRUE, NA)
  }

  pre_ops <- cp02_baseline_pre_reporting_table(vspec)
  bridge_d13 <- if (vspec$pattern_family == "fixed_horizon_bridge") {
    cp02_baseline_bridge_d13_table(vspec)
  } else {
    NULL
  }

  list(
    cohort = base_truth$cohort,
    att_gt = cells,
    simple = base_truth$simple,
    group = base_truth$group,
    dynamic = base_truth$dynamic,
    calendar = base_truth$calendar,
    pre_reporting_operators = pre_ops,
    bridge_d13 = bridge_d13,
    propensity = cp02_baseline_propensity_table(vspec),
    tilt = cp02_baseline_tilt_table(vspec),
    noise = cp02_baseline_noise_moments(vspec),
    pre_note = if (vspec$has_period1_shift) {
      paste(
        "Causal pre-treatment effects are 0 because period-1 shifts s_G(X)",
        "enter both Y_it(0) and Y_it(g). However, pre-period DiD reporting",
        "estimands involving period 1 are nonzero and depend on the control",
        "pool C and baseline b; see truth$pre_reporting_operators."
      )
    } else {
      paste(
        "No period-1 shift is active; conditional parallel trends hold over",
        "all periods 1..5 and both causal pre-effects and pre-period DiD",
        "reporting estimands are 0."
      )
    }
  )
}

cp02_check_variant_n <- function(n, vspec) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n != round(n) || n < 1) {
    stop("`n` must be a single whole number >= 1.", call. = FALSE)
  }
  if (length(vspec$target_spec$periods) * n > .Machine$integer.max) {
    stop("`n` too large: T * n rows must not exceed .Machine$integer.max.",
         call. = FALSE)
  }
  n_int <- as.integer(n)
  if (vspec$cluster_mode != "unit_iid" && (n_int %% vspec$cluster_size != 0L)) {
    stop("Cluster variants require `n` to be a multiple of fixed ",
         "cluster_size = ", vspec$cluster_size, ".", call. = FALSE)
  }
  n_int
}

cp02_baseline_variant_draw <- function(vspec, n, seed, n_clusters = NULL) {
  n <- cp02_check_variant_n(n, vspec)
  if (!is.null(n_clusters)) {
    if (!is.numeric(n_clusters) || length(n_clusters) != 1L ||
        !is.finite(n_clusters) || n_clusters != round(n_clusters) ||
        n_clusters < 1) {
      stop("`n_clusters` must be a positive whole number.", call. = FALSE)
    }
    n_clusters <- as.integer(n_clusters)
    if (vspec$cluster_mode != "unit_iid" &&
        n != n_clusters * vspec$cluster_size) {
      stop("Inconsistent `n` (", n, ") and `n_clusters` (", n_clusters,
           ") for fixed cluster_size = ", vspec$cluster_size, ".",
           call. = FALSE)
    }
  }

  spec <- vspec$target_spec
  periods <- spec$periods
  T_len <- length(periods)
  codes <- c(0, spec$cohorts)
  m <- n * T_len

  tilt <- cp02_baseline_tilt_table(vspec)
  p_sample_mat <- tilt$p_sample_matrix
  p_x1_sample <- sum(p_sample_mat["x1", ])
  cond_g_given_x <- sweep(p_sample_mat, 1L, rowSums(p_sample_mat), "/")
  cum_g <- t(apply(cond_g_given_x, 1L, cumsum))

  pat <- vspec$patterns
  singleton_row <- which(pat$type == "singleton")[1L]

  with_seed(seed, {
    u_x <- stats::runif(n)
    x <- as.integer(u_x < p_x1_sample)
    u_g <- stats::runif(n)
    g_idx <- 1L + rowSums(u_g > cum_g[x + 1L, -ncol(cum_g), drop = FALSE])
    G <- codes[g_idx]

    z <- stats::rnorm(n)
    unit_effect <- stats::rnorm(n, 0, spec$sd_unit)

    if (vspec$cluster_mode == "unit_iid") {
      n_cl <- if (is.null(n_clusters)) spec$n_clusters else n_clusters
      cluster <- sample.int(n_cl, n, replace = TRUE)
      cluster_shocks <- matrix(0, nrow = n_cl, ncol = T_len)
      cluster_obs_flag <- rep(TRUE, n_cl)
    } else {
      n_cl <- n %/% vspec$cluster_size
      cluster <- rep(seq_len(n_cl), each = vspec$cluster_size)
      cluster_shocks <- matrix(
        stats::rnorm(n_cl * T_len, 0, vspec$sd_cluster),
        nrow = n_cl,
        ncol = T_len
      )
      if (vspec$cluster_mode == "cluster_zero_eligible") {
        u_cl_obs <- stats::runif(n_cl)
        cluster_obs_flag <- u_cl_obs < vspec$cluster_obs_prob
      } else {
        cluster_obs_flag <- rep(TRUE, n_cl)
      }
    }

    if (vspec$innovation_family == "gaussian") {
      eps_raw <- stats::rnorm(m)
    } else {
      u_eps <- stats::runif(m)
      eps_raw <- cp02_centered_exp_quantile(u_eps)
    }

    u_pat <- stats::runif(n)
    pidx_base <- cp02_baseline_pattern_index(u_pat, vspec$pattern_family)

    if (vspec$sampling_scheme == "independent_weights") {
      u_w <- stats::runif(n)
      sample_weight <- vspec$indep_weight_bounds["lower"] +
        diff(vspec$indep_weight_bounds) * u_w
    } else if (vspec$sampling_scheme == "unequal_prob_iid") {
      cell_idx <- cbind(x + 1L, match(G, codes))
      sample_weight <- unname(tilt$inv_weight_matrix[cell_idx])
    } else {
      sample_weight <- rep(1, n)
    }
  })

  pidx <- pidx_base
  unit_cl_obs <- cluster_obs_flag[cluster]
  if (vspec$cluster_mode == "cluster_zero_eligible") {
    # Clusters with cluster_obs_flag == FALSE have follow-up masked so pair
    # and full patterns collapse to singleton {3}, yielding 0 eligible units.
    mask_to_singleton <- !unit_cl_obs & (pat$type[pidx] %in% c("pair", "full"))
    pidx[mask_to_singleton] <- singleton_row
  }

  norm_weight <- sample_weight / mean(sample_weight)
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

  shift_idx <- match(Gi, vspec$period1_shift$G)
  shift <- ifelse(
    period == periods[1L],
    vspec$period1_shift$intercept[shift_idx] +
      vspec$period1_shift$slope_x[shift_idx] * xi,
    0
  )

  eps_scale <- spec$sd_eps * unname(vspec$time_scale[t_pos])
  eps <- eps_scale * eps_raw
  c_shock <- cluster_shocks[cbind(cluster[id], t_pos)]

  y0 <- unname(spec$cohort_level[gk]) +
    unname(spec$time_effect[as.character(period)]) +
    spec$x_trend * xi * (period - periods[1L]) +
    spec$z_coef * z[id] +
    unit_effect[id] +
    c_shock +
    eps +
    shift

  treated <- Gi > 0 & period >= Gi
  effect <- numeric(m)
  effect[treated] <- unname(spec$theta[gk[treated]]) +
    spec$delta * (period[treated] - Gi[treated]) +
    spec$kappa * xi[treated]
  # Period-1 shift enters y0 and therefore enters BOTH y0 and y_g = y0 + effect.
  y_g <- y0 + effect

  seen <- matrix(FALSE, nrow = nrow(pat), ncol = T_len)
  for (j in seq_len(nrow(pat))) {
    if (pat$type[j] == "full") {
      seen[j, ] <- TRUE
    } else {
      dates_j <- stats::na.omit(c(pat$date1[j], pat$date2[j]))
      if (length(dates_j) > 0L) {
        seen[j, match(dates_j, periods)] <- TRUE
      }
    }
  }
  observed <- seen[cbind(pidx[id], t_pos)]
  n_obs_unit <- rowSums(seen)[pidx]
  eligible_ge2_unit <- n_obs_unit >= 2L

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
  y1_val <- ifelse(is.na(d1_pos), NA_real_, y_mat[cbind(seq_len(n), d1_pos)])
  y2_val <- ifelse(is.na(d2_pos), NA_real_, y_mat[cbind(seq_len(n), d2_pos)])
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
    dy = ifelse(!is.na(y1_val) & !is.na(y2_val), y2_val - y1_val, NA_real_),
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
  settings$seed <- seed
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
      if (vspec$cluster_mode == "unit_iid") {
        "sample.int(n_clusters, n, replace = TRUE): cluster labels"
      } else {
        c(
          "rnorm(n_clusters * T): cluster-by-period shocks",
          if (vspec$cluster_mode == "cluster_zero_eligible") {
            "runif(n_clusters): cluster observation flag"
          }
        )
      },
      if (vspec$innovation_family == "gaussian") {
        "rnorm(n * T): Gaussian period innovations"
      } else {
        "runif(n * T): centered-exponential period innovations"
      },
      "runif(n): base observation pattern",
      if (vspec$sampling_scheme == "independent_weights") {
        "runif(n): independent bounded weights"
      }
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
        base_pattern = pat$pattern[pidx_base],
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
