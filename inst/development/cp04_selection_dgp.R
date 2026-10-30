# CP04 Corollary 1 validation-fixture DGP family (development only, not package API).
# Pure base R; no external file or package dependencies. Never modifies inputs,
# never fits cdid, and never repairs or redraws finite-sample support failures.

cp04_with_seed <- function(seed, expr) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("`seed` must be a finite integer scalar.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  force(expr)
}

cp04_expit <- function(eta) {
  1 / (1 + exp(-eta))
}

cp04_dgp_catalog <- function() {
  data.frame(
    variant = c(
      "baseline_general_cont",
      "baseline_general_binary",
      "single_pair_witness",
      "consecutive_rotation_cont",
      "gapped_reentry_cont",
      "weighted_general_cont",
      "clustered_general_cont",
      "misspecified_selection_logit"
    ),
    history_family = c(
      "general_pairs",
      "general_pairs",
      "single_pair",
      "consecutive_rotation",
      "gapped_reentry",
      "general_pairs",
      "general_pairs",
      "general_pairs"
    ),
    x_type = c(
      "continuous_bounded",
      "binary",
      "binary",
      "continuous_bounded",
      "continuous_bounded",
      "continuous_bounded",
      "continuous_bounded",
      "continuous_bounded"
    ),
    weight_mode = c(
      "unit",
      "unit",
      "unit",
      "unit",
      "unit",
      "weighted_x",
      "unit",
      "unit"
    ),
    cluster_mode = c(
      "unit_iid",
      "unit_iid",
      "unit_iid",
      "unit_iid",
      "unit_iid",
      "unit_iid",
      "cluster_dependent",
      "unit_iid"
    ),
    selection_spec = c(
      rep("correct_linear", 7L),
      "misspecified_quadratic"
    ),
    valid_nominal_regime = c(rep(TRUE, 7L), FALSE),
    description = c(
      "5-period general exact-two pairs + singleton + zero-outcome; continuous X ~ U(-1,1), G indep X, exact cohort-specific pair logits.",
      "5-period general exact-two pairs + singleton + zero-outcome; binary X in {0,1}, saturated selection and treatment logits.",
      "2-period single-pair selection-bias witness matching population_reference.R (p(0)=1/4, p(1)=3/4, a=1/2, ATT=0.50, naive=0.75).",
      "5-period consecutive rotation exact-two pairs (1-2, 2-3, 3-4, 4-5) + singleton + zero-outcome; non-adjacent pairs unsupported directly.",
      "5-period gapped/reentry exact-two pairs (1-3, 1-4, 1-5, 2-4, 2-5, 3-5) + singleton + zero-outcome; adjacent pairs unsupported directly.",
      "5-period general pairs with unequal fixed positive sampling weights w(X)=1+0.5*X; analytic target E_P[X]=1/6.",
      "5-period general pairs with 5 units/cluster and common cluster-by-period outcome shocks preserving MAR and PT.",
      "NEGATIVE CONTROL (never valid-nominal): true quadratic selection logit in X with quadratic X trend, fitted with linear ~x formula."
    ),
    stringsAsFactors = FALSE
  )
}

cp04_build_pair_table <- function(periods, history_family) {
  pairs_mat <- utils::combn(periods, 2L)
  all_pairs <- data.frame(
    pair = paste(pairs_mat[1L, ], pairs_mat[2L, ], sep = "-"),
    s = as.numeric(pairs_mat[1L, ]),
    t = as.numeric(pairs_mat[2L, ]),
    gap = as.numeric(pairs_mat[2L, ] - pairs_mat[1L, ]),
    stringsAsFactors = FALSE
  )
  if (history_family == "single_pair") {
    all_pairs$supported <- all_pairs$s == periods[1L] & all_pairs$t == periods[2L]
  } else if (history_family == "consecutive_rotation") {
    all_pairs$supported <- all_pairs$gap == 1
  } else if (history_family == "gapped_reentry") {
    all_pairs$supported <- all_pairs$gap >= 2
  } else if (history_family == "general_pairs") {
    all_pairs$supported <- rep(TRUE, nrow(all_pairs))
  } else {
    stop("Unknown `history_family`.", call. = FALSE)
  }
  all_pairs
}

cp04_dgp_spec <- function(
    variant = "baseline_general_cont",
    history_family = NULL,
    x_type = NULL,
    weight_mode = NULL,
    cluster_mode = NULL,
    selection_spec = NULL,
    has_never_treated = TRUE,
    cluster_size = 5L,
    sd_cluster = 0.5) {
  cat_df <- cp04_dgp_catalog()
  if (!is.character(variant) || length(variant) != 1L || is.na(variant) ||
      !variant %in% cat_df$variant) {
    stop("`variant` must be one of cp04_dgp_catalog()$variant.", call. = FALSE)
  }
  row <- cat_df[cat_df$variant == variant, , drop = FALSE]
  if (is.null(history_family)) history_family <- row$history_family
  if (is.null(x_type)) x_type <- row$x_type
  if (is.null(weight_mode)) weight_mode <- row$weight_mode
  if (is.null(cluster_mode)) cluster_mode <- row$cluster_mode
  if (is.null(selection_spec)) selection_spec <- row$selection_spec

  history_family <- match.arg(
    history_family,
    c("general_pairs", "single_pair", "consecutive_rotation", "gapped_reentry")
  )
  x_type <- match.arg(x_type, c("continuous_bounded", "binary"))
  weight_mode <- match.arg(weight_mode, c("unit", "weighted_x"))
  cluster_mode <- match.arg(cluster_mode, c("unit_iid", "cluster_dependent"))
  selection_spec <- match.arg(selection_spec, c("correct_linear", "misspecified_quadratic"))
  stopifnot(
    is.logical(has_never_treated), length(has_never_treated) == 1L, !is.na(has_never_treated),
    is.numeric(cluster_size), length(cluster_size) == 1L, cluster_size >= 2L,
    cluster_size == round(cluster_size),
    is.numeric(sd_cluster), length(sd_cluster) == 1L, is.finite(sd_cluster), sd_cluster >= 0
  )
  if (history_family == "single_pair" && !has_never_treated) {
    stop("`single_pair` requires `has_never_treated = TRUE`.", call. = FALSE)
  }

  if (history_family == "single_pair") {
    periods <- c(1, 2)
    cohorts <- c(2)
    singleton_date <- 1L
  } else {
    periods <- as.numeric(1:5)
    cohorts <- c(3, 4, 5)
    singleton_date <- 3L
  }
  cohort_codes <- if (has_never_treated) c(0, cohorts) else cohorts
  n_cohorts <- length(cohort_codes)
  tau_H <- stats::setNames(rep(1 / n_cohorts, n_cohorts), as.character(cohort_codes))
  # Because G is independent of X and w(X) is a common function of X across cohorts,
  # the target weighted cohort shares tau_P = P(G=k) equal tau_H identically.
  tau_P <- tau_H

  weight_slope <- if (weight_mode == "weighted_x") 0.5 else 0.0
  p_x_binary <- 0.5

  if (x_type == "continuous_bounded") {
    x_support <- c(-1.0, 1.0)
    mean_w_H <- 1.0
    mean_x_H <- 0.0
    mean_x_P <- weight_slope / 3.0
    mean_q2_P <- 0.0 # int_{-1}^1 0.5*(1 + lambda*x)*(3*x^2 - 1) dx = 0 identically
  } else {
    x_support <- c(0.0, 1.0)
    w0 <- 1.0
    w1 <- 1.0 + weight_slope
    mean_w_H <- (1.0 - p_x_binary) * w0 + p_x_binary * w1
    mean_x_H <- p_x_binary
    mean_x_P <- (p_x_binary * w1) / mean_w_H
    mean_q2_P <- (1.0 - mean_x_P) * (-1.0) + mean_x_P * 2.0
  }

  all_pairs <- cp04_build_pair_table(periods, history_family)
  active_pairs <- all_pairs[all_pairs$supported, , drop = FALSE]
  M <- nrow(active_pairs)

  histories <- rbind(
    data.frame(
      history = active_pairs$pair,
      history_type = rep("pair", M),
      date1 = as.integer(active_pairs$s),
      date2 = as.integer(active_pairs$t),
      n_obs = rep(2L, M),
      eligible_ge1 = rep(TRUE, M),
      eligible_ge2 = rep(TRUE, M),
      stringsAsFactors = FALSE
    ),
    data.frame(
      history = c(as.character(singleton_date), "none"),
      history_type = c("singleton", "none"),
      date1 = c(as.integer(singleton_date), NA_integer_),
      date2 = c(NA_integer_, NA_integer_),
      n_obs = c(1L, 0L),
      eligible_ge1 = c(TRUE, FALSE),
      eligible_ge2 = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    )
  )
  rownames(histories) <- NULL
  singleton_split <- if (history_family == "single_pair") 2 / 3 else 0.5

  # Fixed selection logit parameters p_{e,k}(X) = expit(alpha_{e,k} + beta_{e,k}*X + gamma_{e,k}*(3*X^2 - 1)).
  if (history_family == "single_pair") {
    c_shift <- stats::setNames(c(0.0, 0.0), as.character(cohort_codes))
    d_slope <- stats::setNames(c(2 * log(3), 2 * log(3)), as.character(cohort_codes))
    g_quad <- stats::setNames(c(0.0, 0.0), as.character(cohort_codes))
    alpha_base <- -log(3)
    pair_u <- 0.0
    pair_v <- 1.0
  } else {
    c_map <- c("0" = -0.15, "3" = 0.15, "4" = 0.05, "5" = -0.05)
    d_map <- c("0" = -0.35, "3" = 0.40, "4" = 0.25, "5" = -0.20)
    g_map <- if (selection_spec == "misspecified_quadratic") {
      c("0" = -0.25, "3" = 0.25, "4" = 0.15, "5" = -0.15)
    } else {
      c("0" = 0.0, "3" = 0.0, "4" = 0.0, "5" = 0.0)
    }
    c_shift <- c_map[as.character(cohort_codes)]
    d_slope <- d_map[as.character(cohort_codes)]
    g_quad <- g_map[as.character(cohort_codes)]
    alpha_base <- if (selection_spec == "misspecified_quadratic") {
      -3.40
    } else if (M == 4L) {
      -1.90
    } else if (M == 6L) {
      -2.35
    } else {
      -2.90
    }
    s_m <- if (M > 1L) {
      (seq_len(M) - (M + 1) / 2) / ((M - 1) / 2)
    } else {
      0.0
    }
    pair_u <- 0.10 * s_m
    pair_v <- 1.0 + 0.15 * s_m
  }

  param_rows <- list()
  for (k_chr in as.character(cohort_codes)) {
    k_num <- as.numeric(k_chr)
    act_idx <- 0L
    for (j in seq_len(nrow(all_pairs))) {
      sup <- all_pairs$supported[j]
      if (sup) {
        act_idx <- act_idx + 1L
        alp <- unname(alpha_base + c_shift[k_chr] + pair_u[act_idx])
        bet <- unname(d_slope[k_chr] * pair_v[act_idx])
        gam <- unname(g_quad[k_chr])
        if (x_type == "binary") {
          eta_vals <- c(alp - gam, alp + bet + 2 * gam)
          eta_sup <- max(eta_vals)
          eta_inf <- min(eta_vals)
        } else {
          eta_sup <- alp + abs(bet) + max(-gam, 2 * gam)
          eta_inf <- alp - abs(bet) + min(-gam, 2 * gam)
        }
        p_up <- cp04_expit(eta_sup)
        p_lo <- cp04_expit(eta_inf)
      } else {
        alp <- NA_real_
        bet <- NA_real_
        gam <- NA_real_
        eta_sup <- -Inf
        eta_inf <- -Inf
        p_up <- 0.0
        p_lo <- 0.0
      }
      param_rows[[length(param_rows) + 1L]] <- data.frame(
        pair = all_pairs$pair[j],
        s = all_pairs$s[j],
        t = all_pairs$t[j],
        G = k_num,
        supported = sup,
        alpha = alp,
        beta = bet,
        gamma_quad = gam,
        eta_sup = eta_sup,
        p_upper_bound = p_up,
        p_min_bound = p_lo,
        stringsAsFactors = FALSE
      )
    }
  }
  pair_logit_params <- do.call(rbind, param_rows)
  rownames(pair_logit_params) <- NULL

  cohort_bounds <- do.call(rbind, lapply(cohort_codes, function(k) {
    sub_df <- pair_logit_params[pair_logit_params$G == k & pair_logit_params$supported, , drop = FALSE]
    sum_up <- sum(sub_df$p_upper_bound)
    sum_lo <- sum(sub_df$p_min_bound)
    if (!is.finite(sum_up) || sum_up >= 1.0) {
      stop("Analytic pair probability upper bound sum_e p_e,k(X) must be < 1.", call. = FALSE)
    }
    data.frame(
      G = k,
      n_supported_pairs = nrow(sub_df),
      sum_p_upper_bound = sum_up,
      sum_p_lower_bound = sum_lo,
      min_residual_prob = 1.0 - sum_up,
      stringsAsFactors = FALSE
    )
  }))

  if (history_family == "single_pair") {
    cohort_level <- stats::setNames(c(0.0, 0.0), as.character(cohort_codes))
    time_effect <- stats::setNames(c(0.0, 0.0), as.character(periods))
    theta <- stats::setNames(c(0.5), as.character(cohorts))
    delta <- 0.0
    kappa <- 1.0
    x_trend <- 0.0
    x_quad_trend <- 0.0
    z_coef <- 0.0
    sd_unit <- 0.5
    sd_eps <- 0.5
  } else {
    full_levels <- c("0" = 0.0, "3" = 0.5, "4" = -0.5, "5" = 1.0)
    cohort_level <- full_levels[as.character(cohort_codes)]
    time_effect <- stats::setNames(0.25 * periods, as.character(periods))
    theta <- stats::setNames(c(1.0, 2.0, 1.5), as.character(cohorts))
    delta <- 0.5
    kappa <- 1.0
    x_trend <- 0.5
    x_quad_trend <- if (selection_spec == "misspecified_quadratic") 0.4 else 0.0
    z_coef <- 0.5
    sd_unit <- 0.5
    sd_eps <- 0.5
  }

  valid_nominal <- (selection_spec == "correct_linear")
  anticipation <- 0
  refs <- stats::setNames(
    vapply(cohorts, function(g) max(periods[periods + anticipation < g]), numeric(1L)),
    as.character(cohorts)
  )

  list(
    variant = variant,
    history_family = history_family,
    x_type = x_type,
    weight_mode = weight_mode,
    cluster_mode = cluster_mode,
    selection_spec = selection_spec,
    valid_nominal_regime = valid_nominal,
    periods = periods,
    cohorts = cohorts,
    cohort_codes = cohort_codes,
    never_treated_code = if (has_never_treated) 0 else NA_real_,
    has_never_treated = has_never_treated,
    anticipation = anticipation,
    reference = refs,
    tau_H = tau_H,
    tau_P = tau_P,
    x_support = x_support,
    p_x_binary = p_x_binary,
    weight_slope = weight_slope,
    mean_w_H = mean_w_H,
    mean_x_H = mean_x_H,
    mean_x_P = mean_x_P,
    mean_q2_P = mean_q2_P,
    all_pairs = all_pairs,
    histories = histories,
    singleton_date = singleton_date,
    singleton_split = singleton_split,
    pair_logit_params = pair_logit_params,
    cohort_bounds = cohort_bounds,
    cohort_level = cohort_level,
    time_effect = time_effect,
    theta = theta,
    delta = delta,
    kappa = kappa,
    x_trend = x_trend,
    x_quad_trend = x_quad_trend,
    z_coef = z_coef,
    sd_unit = sd_unit,
    sd_eps = sd_eps,
    sd_cluster = if (cluster_mode == "cluster_dependent") sd_cluster else 0.0,
    cluster_size = as.integer(cluster_size),
    formulas = list(
      default_selection_xformla = "~ x",
      default_treatment_xformla = "~ x",
      true_selection_xformla = if (selection_spec == "correct_linear") {
        "~ x"
      } else {
        "~ x + I(3 * x^2 - 1)"
      },
      true_treatment_xformla = "~ 1",
      redundant_zero_treatment_covariates = "~ x"
    ),
    regime_metadata = list(
      variant = variant,
      valid_nominal_regime = valid_nominal,
      regime_role = if (valid_nominal) {
        "valid_nominal_corollary1_fixture"
      } else {
        "negative_control_misspecified_selection_logit"
      },
      supported_controls = if (has_never_treated) c("NT", "NYT") else c("NYT"),
      supported_pt_modes = c("PT-post", "PT-all"),
      default_selection = "none",
      corollary1_fitting_universe = "F",
      corollary1_models = "fitted_actual_pair_cohort_logits",
      probability_role = "fixture_oracle_only",
      unsupported_regimes = c(
        "uncorrected_selection_none_is_biased",
        "marginal_period_product_selection_probabilities_invalid",
        "observer_only_selection_logit_dropping_A0_invalid",
        "postdraw_minimum_two_repair_or_redraw_prohibited",
        if (!valid_nominal) "linear_selection_logit_misspecified_identifies_pseudo_functional_only"
      )
    )
  )
}

cp04_eval_weight <- function(spec, x) {
  stopifnot(is.numeric(x), all(is.finite(x)))
  w <- 1.0 + spec$weight_slope * x
  if (any(w <= 0)) {
    stop("Sampling weight w(x) must be strictly positive on support.", call. = FALSE)
  }
  w
}

cp04_eval_pair_prob <- function(spec, G, x, pair = NULL) {
  n <- max(length(G), length(x))
  stopifnot(
    n >= 1L,
    length(G) %in% c(1L, n),
    length(x) %in% c(1L, n),
    is.numeric(G), all(is.finite(G)), all(G %in% spec$cohort_codes),
    is.numeric(x), all(is.finite(x))
  )
  G <- rep_len(as.numeric(G), n)
  x <- rep_len(as.numeric(x), n)
  pairs_vec <- spec$all_pairs$pair
  if (!is.null(pair)) {
    stopifnot(is.character(pair), length(pair) == 1L, pair %in% pairs_vec)
    pairs_vec <- pair
  }
  out <- matrix(0.0, nrow = n, ncol = length(pairs_vec))
  colnames(out) <- pairs_vec
  ptab <- spec$pair_logit_params
  q2 <- 3.0 * x^2 - 1.0
  for (j in seq_along(pairs_vec)) {
    pr <- pairs_vec[j]
    sub_p <- ptab[ptab$pair == pr, , drop = FALSE]
    idx <- match(G, sub_p$G)
    sup <- sub_p$supported[idx]
    if (any(sup)) {
      eta <- sub_p$alpha[idx[sup]] +
        sub_p$beta[idx[sup]] * x[sup] +
        sub_p$gamma_quad[idx[sup]] * q2[sup]
      out[sup, j] <- cp04_expit(eta)
    }
  }
  if (!is.null(pair)) {
    return(as.numeric(out[, 1L]))
  }
  out
}

cp04_eval_history_probs <- function(spec, G, x) {
  p_all <- cp04_eval_pair_prob(spec, G, x)
  active_names <- spec$all_pairs$pair[spec$all_pairs$supported]
  p_active <- p_all[, active_names, drop = FALSE]
  sum_active <- rowSums(p_active)
  residual <- 1.0 - sum_active
  if (any(residual <= 0)) {
    stop("Analytic bound violated: sum of active pair probabilities >= 1.", call. = FALSE)
  }
  p_sing <- spec$singleton_split * residual
  p_none <- (1.0 - spec$singleton_split) * residual
  out <- cbind(p_active, p_sing, p_none)
  colnames(out) <- spec$histories$history
  out
}

cp04_eval_potential <- function(
    spec,
    G,
    period,
    x,
    z = 0.0,
    unit_effect = 0.0,
    cluster_shock = 0.0,
    eps = 0.0) {
  args <- list(G, period, x, z, unit_effect, cluster_shock, eps)
  n <- max(lengths(args))
  stopifnot(
    n >= 1L,
    all(lengths(args) %in% c(1L, n)),
    all(vapply(args, function(a) is.numeric(a) && all(is.finite(a)), logical(1L)))
  )
  G <- rep_len(as.numeric(G), n)
  period <- rep_len(as.numeric(period), n)
  x <- rep_len(as.numeric(x), n)
  z <- rep_len(as.numeric(z), n)
  unit_effect <- rep_len(as.numeric(unit_effect), n)
  cluster_shock <- rep_len(as.numeric(cluster_shock), n)
  eps <- rep_len(as.numeric(eps), n)
  stopifnot(all(G %in% spec$cohort_codes), all(period %in% spec$periods))

  gk <- as.character(G)
  pk <- as.character(period)
  t1 <- spec$periods[1L]
  q2 <- 3.0 * x^2 - 1.0
  trend <- (spec$x_trend * x + spec$x_quad_trend * q2) * (period - t1)

  y0 <- unname(spec$cohort_level[gk]) +
    unname(spec$time_effect[pk]) +
    trend +
    spec$z_coef * z +
    unit_effect +
    cluster_shock +
    eps

  effect <- numeric(n)
  post <- G > 0 & period >= G
  if (any(post)) {
    effect[post] <- unname(spec$theta[gk[post]]) +
      spec$delta * (period[post] - G[post]) +
      spec$kappa * (x[post] - spec$mean_x_P)
  }
  data.frame(
    y0 = y0,
    y_g = y0 + effect,
    effect = effect,
    trend = trend,
    eps = eps,
    cluster_shock = cluster_shock,
    stringsAsFactors = FALSE
  )
}

# One-dimensional population means: established base-R integration or exact
# binary enumeration. No population fitting or pseudo-functional optimization.
cp04_mean <- function(spec, f, weighted = TRUE) {
  integrand <- function(x) f(x) * if (weighted) cp04_eval_weight(spec, x) / spec$mean_w_H else 1
  if (spec$x_type == "binary") {
    sum(c(1 - spec$p_x_binary, spec$p_x_binary) * integrand(c(0, 1)))
  } else {
    stats::integrate(function(x) integrand(x) / 2, -1, 1,
                     rel.tol = 1e-12, abs.tol = 1e-13)$value
  }
}

# All qualified controls are untreated at both endpoints. Shared untreated
# trends and common cohort X laws give mu(g,s,t,K)=ATT(g,t)-ATT(g,s), including
# zero pre-comparisons. Retain actual D10/D13 pools and horizon provenance.
cp04_comparisons <- function(spec, eta, control) {
  rows <- list()
  effect <- function(g, t) ifelse(t >= g, spec$theta[as.character(g)] + spec$delta * (t - g), 0)
  for (g in spec$cohorts) for (j in seq_len(nrow(spec$all_pairs))) {
    pair <- spec$all_pairs[j, ]
    cutoffs <- if (control == "NT") pair$t else
      unique(c(pair$t, spec$periods[spec$periods >= max(g, pair$t)]))
    pools <- list()
    for (h in cutoffs) {
      K <- if (spec$has_never_treated) 0 else numeric()
      if (control == "NYT") K <- c(K, spec$cohorts[spec$cohorts > g & spec$cohorts > h + spec$anticipation])
      if (!length(K)) next
      key <- paste(K, collapse = ",")
      pools[[key]] <- c(pools[[key]], h)
    }
    for (key in names(pools)) {
      K <- as.numeric(strsplit(key, ",", fixed = TRUE)[[1L]])
      row <- data.frame(group = g, s = pair$s, t = pair$t, pair = pair$pair,
                        pair_supported = pair$supported, reference = unname(spec$reference[as.character(g)]),
                        pool = key, original = pair$t %in% pools[[key]],
                        beta_s = unname(effect(g, pair$s)), beta_t = unname(effect(g, pair$t)),
                        mu = unname(effect(g, pair$t) - effect(g, pair$s)),
                        true_treat_intercept = log(eta[as.character(g)] / sum(eta[as.character(K)])),
                        true_treat_slope_x = 0, row.names = NULL)
      row$horizons <- list(if (control == "NYT") pools[[key]][pools[[key]] >= g] else numeric())
      rows[[length(rows) + 1L]] <- row
    }
  }
  all <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(control = control, rows = all,
       supported_rows = if (nrow(all)) all[all$pair_supported, , drop = FALSE] else all)
}

# Availability under the PT-post comparison graph: pre-dates retain their pool
# identity; post-dates share ATT coordinates. A reference endpoint anchors zero.
# Reachability replaces the old generic operator-coordinate/incidence/QR stack.
cp04_available <- function(spec, rows, g) {
  dates <- spec$periods[spec$periods >= g]
  if (is.null(rows) || !nrow(rows)) return(rep(FALSE, length(dates)))
  rows <- rows[rows$group == g, , drop = FALSE]
  key <- function(t, K) ifelse(t == spec$reference[as.character(g)], "reference",
                              ifelse(t >= g, paste0("ATT:", t), paste0(K, ":pre:", t)))
  left <- key(rows$s, rows$pool)
  right <- key(rows$t, rows$pool)
  reached <- "reference"
  repeat {
    next_set <- unique(c(reached, right[left %in% reached], left[right %in% reached]))
    if (length(next_set) == length(reached)) break
    reached <- next_set
  }
  paste0("ATT:", dates) %in% reached
}

cp04_dgp_truth <- function(spec) {
  ptab <- spec$pair_logit_params
  ptab$mean_p_H <- ptab$mean_p_P <- ptab$cond_pair_given_A_P <- 0
  for (r in which(ptab$supported)) {
    f <- function(x) cp04_eval_pair_prob(spec, ptab$G[r], x, ptab$pair[r])
    ptab$mean_p_H[r] <- cp04_mean(spec, f, weighted = FALSE)
    ptab$mean_p_P[r] <- cp04_mean(spec, f)
  }
  rates <- do.call(rbind, lapply(spec$cohort_codes, function(k) {
    idx <- ptab$G == k
    a_H <- sum(ptab$mean_p_H[idx])
    a_P <- sum(ptab$mean_p_P[idx])
    data.frame(G = k, tau_H = unname(spec$tau_H[as.character(k)]),
               tau_P = unname(spec$tau_P[as.character(k)]), mean_x_P = spec$mean_x_P,
               a_H = a_H, a_P = a_P, singleton_rate_P = spec$singleton_split * (1 - a_P),
               zero_obs_rate_P = (1 - spec$singleton_split) * (1 - a_P))
  }))
  rates$sum_p_upper_bound <- spec$cohort_bounds$sum_p_upper_bound
  rates$min_residual_prob <- spec$cohort_bounds$min_residual_prob
  rates$eta_Q <- rates$tau_P * rates$a_P / sum(rates$tau_P * rates$a_P)
  a <- stats::setNames(rates$a_P, rates$G)
  eta <- stats::setNames(rates$eta_Q, rates$G)
  ptab$cond_pair_given_A_P <- ptab$mean_p_P / a[as.character(ptab$G)]
  cohort <- data.frame(group = rates$G, pi = rates$tau_P, mean_x = rates$mean_x_P,
                       repeat_rate_a = rates$a_P, eta_Q = rates$eta_Q,
                       contract_group = ifelse(rates$G == 0, Inf, rates$G))
  cells <- expand.grid(t = spec$periods, group = spec$cohorts)[, 2:1]
  cells$post <- cells$t >= cells$group
  cells$e <- cells$t - cells$group
  cells$att <- ifelse(cells$post, spec$theta[as.character(cells$group)] + spec$delta * cells$e, 0)
  comparison <- list(NT = if (spec$has_never_treated) cp04_comparisons(spec, eta, "NT") else NULL,
                     NYT = cp04_comparisons(spec, eta, "NYT"))
  for (ctl in names(comparison)) {
    cells[[paste0("available_", ctl)]] <- NA
    rows <- comparison[[ctl]]$supported_rows
    for (g in spec$cohorts) {
      idx <- cells$group == g & cells$post
      cells[[paste0("available_", ctl)]][idx] <- cp04_available(spec, rows, g)
    }
  }

  # Explicit population cell weights, independent of histories or repeat shares.
  tau <- unname(spec$tau_P[as.character(cells$group)])
  post <- cells$post
  weights <- cells[, c("group", "t", "e", "post")]
  weights$simple <- ifelse(post, tau / sum(tau[post]), 0)
  weights$group <- ifelse(post, tau / sum(spec$tau_P[as.character(spec$cohorts)]) /
                           (max(spec$periods) - cells$group + 1), 0)
  weights$dynamic <- weights$calendar <- numeric(nrow(cells))
  dynamic_by <- data.frame(e = sort(unique(cells$e)), att = NA_real_)
  dynamic_weights <- calendar_weights <- numeric(nrow(cells))
  for (e in dynamic_by$e) {
    idx <- cells$e == e
    dynamic_weights[idx] <- tau[idx] / sum(tau[idx])
    dynamic_by$att[dynamic_by$e == e] <- sum(dynamic_weights[idx] * cells$att[idx])
    if (e >= 0) weights$dynamic[idx] <- dynamic_weights[idx] / sum(dynamic_by$e >= 0)
  }
  calendar_by <- data.frame(t = spec$periods[spec$periods >= min(spec$cohorts)], att = NA_real_)
  for (t in calendar_by$t) {
    idx <- cells$t == t & post
    calendar_weights[idx] <- tau[idx] / sum(tau[idx])
    calendar_by$att[calendar_by$t == t] <- sum(calendar_weights[idx] * cells$att[idx])
    weights$calendar[idx] <- calendar_weights[idx] / nrow(calendar_by)
  }
  group_by <- data.frame(group = spec$cohorts,
                         weight_tau_P = unname(spec$tau_P[as.character(spec$cohorts)]) /
                           sum(spec$tau_P[as.character(spec$cohorts)]),
                         att = vapply(spec$cohorts, function(g) mean(cells$att[cells$group == g & post]), numeric(1)))
  list(cohort = cohort, att_gt = cells, aggregate_weights = weights,
       simple = sum(weights$simple * cells$att),
       group = list(by = group_by, overall = sum(weights$group * cells$att)),
       dynamic = list(by = dynamic_by, overall = sum(weights$dynamic * cells$att),
                      weights = data.frame(cells[, c("group", "t", "e")], weight = dynamic_weights)),
       calendar = list(by = calendar_by, overall = sum(weights$calendar * cells$att),
                       weights = data.frame(cells[, c("group", "t")], weight = calendar_weights)),
       selection = list(pair_models = ptab, cohort_repeat_rates = rates,
                        eval_pair_prob = function(s, t, G, x) cp04_eval_pair_prob(spec, G, x, paste(s, t, sep = "-")),
                        eval_repeat_prob = function(G, x) rowSums(cp04_eval_pair_prob(spec, G, x))),
       comparison = comparison, formulas = spec$formulas, regime_metadata = spec$regime_metadata)
}

cp04_validate_draw_inputs <- function(spec, n, seed) {
  stopifnot(is.list(spec), !is.null(spec$variant), !is.null(spec$periods))
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != round(n) || (n * length(spec$periods)) > .Machine$integer.max) {
    stop("`n` must be a positive integer with n * length(periods) in integer range.", call. = FALSE)
  }
  if (spec$cluster_mode == "cluster_dependent" && (as.integer(n) %% spec$cluster_size != 0L)) {
    stop("For `cluster_dependent`, `n` must be a multiple of `cluster_size`.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("`seed` must be a finite integer scalar.", call. = FALSE)
  }
  invisible(TRUE)
}

cp04_dgp_draw <- function(spec = cp04_dgp_spec(), n = 120L, seed = 20260401L) {
  cp04_validate_draw_inputs(spec, n, seed)
  n <- as.integer(n)
  periods <- spec$periods
  T_len <- length(periods)
  codes <- spec$cohort_codes
  cum_tau <- cumsum(unname(spec$tau_H))
  n_cl <- if (spec$cluster_mode == "cluster_dependent") n %/% spec$cluster_size else n

  # Single forward RNG pass; no redraws or post-draw support repairs.
  cp04_with_seed(seed, {
    u_x <- stats::runif(n)
    x <- if (spec$x_type == "binary") {
      as.numeric(u_x < spec$p_x_binary)
    } else {
      -1.0 + 2.0 * u_x
    }
    u_g <- stats::runif(n)
    g_idx <- 1L + rowSums(outer(u_g, cum_tau[-length(cum_tau)], ">"))
    G <- codes[g_idx]

    z <- stats::rnorm(n)
    unit_effect <- stats::rnorm(n, mean = 0.0, sd = spec$sd_unit)
    cluster_shocks <- if (spec$cluster_mode == "cluster_dependent") {
      matrix(stats::rnorm(n_cl * T_len, mean = 0.0, sd = spec$sd_cluster), nrow = n_cl, ncol = T_len)
    } else {
      matrix(0.0, nrow = n_cl, ncol = T_len)
    }
    eps_mat <- matrix(stats::rnorm(n * T_len, mean = 0.0, sd = spec$sd_eps), nrow = n, ncol = T_len, byrow = TRUE)
    u_sel <- stats::runif(n)
  })

  cluster <- if (spec$cluster_mode == "cluster_dependent") {
    rep(seq_len(n_cl), each = spec$cluster_size)
  } else {
    seq_len(n)
  }

  sample_weight <- cp04_eval_weight(spec, x)
  norm_weight <- sample_weight / mean(sample_weight)

  # Mutually exclusive history draw from exact conditional probabilities
  h_probs <- cp04_eval_history_probs(spec, G, x)
  cum_h <- t(apply(h_probs, 1L, cumsum))
  h_idx <- 1L + rowSums(u_sel > cum_h[, -ncol(cum_h), drop = FALSE])
  hist_df <- spec$histories[h_idx, , drop = FALSE]

  obs_mask <- matrix(FALSE, nrow = n, ncol = T_len)
  colnames(obs_mask) <- as.character(periods)
  for (i in seq_len(n)) {
    d_i <- stats::na.omit(c(hist_df$date1[i], hist_df$date2[i]))
    if (length(d_i) > 0L) {
      obs_mask[i, match(d_i, periods)] <- TRUE
    }
  }

  # Actual pair indicators S_{i,e} = S_{i,s} * S_{i,t}
  n_pairs <- nrow(spec$all_pairs)
  pair_indicators <- matrix(0L, nrow = n, ncol = n_pairs)
  colnames(pair_indicators) <- spec$all_pairs$pair
  for (j in seq_len(n_pairs)) {
    s_col <- match(spec$all_pairs$s[j], periods)
    t_col <- match(spec$all_pairs$t[j], periods)
    pair_indicators[, j] <- as.integer(obs_mask[, s_col] & obs_mask[, t_col])
  }

  frame_df <- data.frame(
    id = seq_len(n),
    cluster = as.integer(cluster),
    G = G,
    x = x,
    z = z,
    sample_weight = sample_weight,
    norm_weight = norm_weight,
    history = hist_df$history,
    history_type = hist_df$history_type,
    date1 = hist_df$date1,
    date2 = hist_df$date2,
    n_obs = as.integer(rowSums(obs_mask)),
    eligible_ge1 = rowSums(obs_mask) >= 1L,
    eligible_ge2 = rowSums(obs_mask) >= 2L,
    stringsAsFactors = FALSE
  )

  # Long panel construction (unit-major, period-minor)
  id_long <- rep(seq_len(n), each = T_len)
  period_long <- rep(periods, times = n)
  t_pos_long <- rep(seq_len(T_len), times = n)
  c_shock_long <- cluster_shocks[cbind(cluster[id_long], t_pos_long)]
  eps_long <- as.vector(t(eps_mat))
  obs_long <- as.vector(t(obs_mask))

  pot_long <- cp04_eval_potential(
    spec = spec,
    G = G[id_long],
    period = period_long,
    x = x[id_long],
    z = z[id_long],
    unit_effect = unit_effect[id_long],
    cluster_shock = c_shock_long,
    eps = eps_long
  )

  panel_df <- data.frame(
    id = id_long,
    period = period_long,
    cluster = frame_df$cluster[id_long],
    G = frame_df$G[id_long],
    x = frame_df$x[id_long],
    z = frame_df$z[id_long],
    sample_weight = frame_df$sample_weight[id_long],
    norm_weight = frame_df$norm_weight[id_long],
    history = frame_df$history[id_long],
    history_type = frame_df$history_type[id_long],
    n_obs = frame_df$n_obs[id_long],
    eligible_ge1 = frame_df$eligible_ge1[id_long],
    eligible_ge2 = frame_df$eligible_ge2[id_long],
    observed = obs_long,
    y = ifelse(obs_long, pot_long$y_g, NA_real_),
    stringsAsFactors = FALSE
  )
  observed_df <- panel_df[panel_df$observed, , drop = FALSE]
  rownames(observed_df) <- NULL

  # Cluster summary table across all frame clusters
  cl_ids <- seq_len(n_cl)
  cl_n_units <- as.integer(tabulate(frame_df$cluster, nbins = n_cl))
  cl_n_elig2 <- as.integer(tabulate(frame_df$cluster[frame_df$eligible_ge2], nbins = n_cl))
  cl_n_sing <- as.integer(tabulate(frame_df$cluster[frame_df$n_obs == 1L], nbins = n_cl))
  cl_n_none <- as.integer(tabulate(frame_df$cluster[frame_df$n_obs == 0L], nbins = n_cl))
  clusters_df <- data.frame(
    cluster = cl_ids,
    n_units = cl_n_units,
    n_eligible_ge2 = cl_n_elig2,
    n_singleton = cl_n_sing,
    n_none = cl_n_none,
    zero_eligible = cl_n_elig2 == 0L,
    stringsAsFactors = FALSE
  )

  counts_df <- as.data.frame(
    table(
      G = factor(frame_df$G, levels = codes),
      history = factor(frame_df$history, levels = spec$histories$history)
    ),
    responseName = "units",
    stringsAsFactors = FALSE
  )
  counts_df$G <- as.numeric(counts_df$G)

  # Oracle pair probabilities and Corollary 1 multipliers
  truth_obj <- cp04_dgp_truth(spec)
  pair_prob_mat <- cp04_eval_pair_prob(spec, G, x)
  a_prob_unit <- rowSums(pair_prob_mat)
  a_P_map <- stats::setNames(
    truth_obj$selection$cohort_repeat_rates$a_P,
    as.character(truth_obj$selection$cohort_repeat_rates$G)
  )
  a_rate_cohort <- unname(a_P_map[as.character(G)])
  oracle_b_obs <- rep(NA_real_, n)
  is_pair <- frame_df$eligible_ge2
  if (any(is_pair)) {
    pair_col <- match(frame_df$history[is_pair], colnames(pair_prob_mat))
    p_obs <- pair_prob_mat[cbind(which(is_pair), pair_col)]
    oracle_b_obs[is_pair] <- a_rate_cohort[is_pair] / p_obs
  }

  # Support diagnostics (never repaired or redrawn)
  cp_rows <- list()
  for (k in codes) {
    in_k <- frame_df$G == k
    for (j in seq_len(n_pairs)) {
      pr <- spec$all_pairs$pair[j]
      sup <- spec$all_pairs$supported[j]
      cnt <- sum(pair_indicators[in_k, j])
      cp_rows[[length(cp_rows) + 1L]] <- data.frame(
        G = k,
        pair = pr,
        supported = sup,
        realized_count = as.integer(cnt),
        empirical_support_failure = sup && (cnt == 0L),
        stringsAsFactors = FALSE
      )
    }
  }
  cp_counts_df <- do.call(rbind, cp_rows)
  empty_sup <- cp_counts_df[cp_counts_df$empirical_support_failure, , drop = FALSE]
  empty_labels <- if (nrow(empty_sup)) {
    paste0("G=", empty_sup$G, ":pair=", empty_sup$pair)
  } else {
    character(0L)
  }

  settings <- spec
  settings$n <- n
  settings$n_clusters_realized <- n_cl
  settings$seed <- as.integer(seed)
  settings$rng <- list(
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection",
    r_version = R.version.string,
    draw_order = c(
      "runif(n): covariate X",
      "runif(n): cohort G (independent of X)",
      "rnorm(n): covariate Z",
      "rnorm(n): unit effect U",
      if (spec$cluster_mode == "cluster_dependent") "rnorm(n_clusters * T): cluster-by-period shocks",
      "rnorm(n * T): period innovations eps",
      "runif(n): mutually exclusive history selection U_sel"
    ),
    postdraw_repair = FALSE
  )

  list(
    data = list(
      frame = frame_df,
      units = frame_df,
      panel = panel_df,
      observed = observed_df,
      pair_indicators = pair_indicators,
      observation_mask = obs_mask,
      clusters = clusters_df,
      counts = counts_df
    ),
    oracle = list(
      panel = cbind(data.frame(id = id_long, period = period_long), pot_long),
      units = data.frame(
        id = seq_len(n),
        unit_effect = unit_effect,
        a_prob_unit = a_prob_unit,
        a_rate_cohort = a_rate_cohort,
        oracle_b_observed = oracle_b_obs,
        stringsAsFactors = FALSE
      ),
      pair_probabilities = pair_prob_mat,
      truth = truth_obj
    ),
    support_diagnostics = list(
      postdraw_repair = FALSE,
      rng_passes = 1L,
      n_frame = n,
      n_eligible_ge1 = sum(frame_df$eligible_ge1),
      n_eligible_ge2 = sum(frame_df$eligible_ge2),
      n_singleton = sum(frame_df$n_obs == 1L),
      n_zero_outcome = sum(frame_df$n_obs == 0L),
      cohort_pair_counts = cp_counts_df,
      n_empty_supported_cohort_pair_cells = nrow(empty_sup),
      empty_supported_cells = empty_labels
    ),
    settings = settings
  )
}
