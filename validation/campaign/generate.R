# Portable reliability campaign: data generation and population truth helpers.
# Inert when sourced or parsed; executes no RNG or estimation on load.
expit_num <- function(x) 1.0 / (1.0 + exp(-x))
# Explicit finite dispatch table for all 64 recipes in contract/RECIPES.csv.
# Unknown recipe IDs are refused and never fall through to a base law.
make_recipe_specs <- function() {
  spec <- function(stream_key, stream_index, family, core_sweep = FALSE, N = 2000L, T = 5L, calendar = seq_len(T), p = 5L,
                   xformla = "~X1+X2+X3+X4", cohort_law = "default_0345", untreated_law = "default", effect_law = "default",
                   obs_law = "comp_arbitrary", error_law = "normal", het_law = "hom", cluster_C = 0L, derived = "none",
                   truth_law = "default_0345") {
    list(
      stream_key = stream_key, stream_index = as.integer(stream_index), family = family, core_sweep = isTRUE(core_sweep),
      N = as.integer(N), T = as.integer(T), calendar = as.double(calendar), p = as.integer(p), xformla = xformla, cohort_law = cohort_law,
      untreated_law = untreated_law, effect_law = effect_law, obs_law = obs_law, error_law = error_law, het_law = het_law,
      cluster_C = as.integer(cluster_C), derived = derived, truth_law = truth_law
    )
  }; x20 <- "~X1+X2+X3+X4+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20"; xv <- "~X1+V+X3+X4"; list(
    COMP_adjacent = spec("COMP_adjacent", 1L, "COMPOSITION", TRUE, obs_law = "comp_adjacent"),
    COMP_arbitrary = spec("COMP_arbitrary", 2L, "COMPOSITION", TRUE, obs_law = "comp_arbitrary"),
    CC_BASE = spec("CC_BASE", 3L, "SELECTED_CC", TRUE, obs_law = "cc_uniform"),
    CC_COMPOSITION = spec("CC_COMPOSITION", 4L, "SELECTED_CC", TRUE, obs_law = "cc_comp"),
    BAL_BASE = spec("BAL_BASE", 5L, "SHARED_BALANCED", TRUE, obs_law = "complete"),
    C_MCAR_DATE = spec("C_MCAR_DATE", 6L, "SHARED_MCAR", TRUE, obs_law = "mcar_date_65"),
    C_MCAR_PAIR = spec("C_MCAR_PAIR", 7L, "SHARED_MCAR", TRUE, obs_law = "mcar_pair"),
    COMP_LADDER_01 = spec("COMP_LADDER_01", 8L, "COMPOSITION_LADDER", FALSE, N = 500L, T = 5L, p = 5L, obs_law = "comp_arbitrary"),
    COMP_LADDER_02 = spec("COMP_LADDER_02", 9L, "COMPOSITION_LADDER", FALSE, N = 2000L, T = 5L, p = 20L, xformla = x20, obs_law = "comp_arbitrary", het_law = "het_x1"),
    COMP_LADDER_03 = spec("COMP_LADDER_03", 10L, "COMPOSITION_LADDER", FALSE, N = 8000L, T = 5L, p = 5L, obs_law = "comp_arbitrary", error_law = "exp"),
    COMP_LADDER_04 = spec("COMP_LADDER_04", 11L, "COMPOSITION_LADDER", FALSE, N = 32000L, T = 5L, p = 5L, obs_law = "comp_arbitrary"),
    COMP_LADDER_05 = spec("COMP_LADDER_05", 12L, "COMPOSITION_LADDER", FALSE, N = 8000L, T = 8L, p = 20L, xformla = x20, obs_law = "comp_arbitrary", error_law = "t5", het_law = "het_x1"),
    COMP_LADDER_06 = spec("COMP_LADDER_06", 13L, "COMPOSITION_LADDER", FALSE, N = 2000L, T = 8L, p = 5L, obs_law = "comp_arbitrary", error_law = "exp", het_law = "het_x1"),
    COMP_LADDER_07 = spec("COMP_LADDER_07", 14L, "COMPOSITION_LADDER", FALSE, N = 2000L, T = 12L, p = 5L, obs_law = "comp_arbitrary"),
    COMP_LADDER_08 = spec("COMP_LADDER_08", 15L, "COMPOSITION_LADDER", FALSE, N = 8000L, T = 12L, p = 20L, xformla = x20, obs_law = "comp_arbitrary", error_law = "t5"),
    CC_LADDER_01 = spec("CC_LADDER_01", 16L, "SELECTED_CC", FALSE, N = 500L, T = 5L, p = 5L, obs_law = "cc_uniform", het_law = "het_x1"),
    CC_LADDER_02 = spec("CC_LADDER_02", 17L, "SELECTED_CC", FALSE, N = 2000L, T = 5L, p = 20L, xformla = x20, obs_law = "cc_uniform", error_law = "exp", het_law = "het_x1"),
    CC_LADDER_03 = spec("CC_LADDER_03", 18L, "SELECTED_CC", FALSE, N = 8000L, T = 5L, p = 20L, xformla = x20, obs_law = "cc_uniform", error_law = "t5", het_law = "het_x1"),
    CC_LADDER_04 = spec("CC_LADDER_04", 19L, "SELECTED_CC", FALSE, N = 32000L, T = 5L, p = 5L, obs_law = "cc_uniform", het_law = "het_x1"),
    COMP_LATE = spec("COMP_LATE", 20L, "COMPOSITION", FALSE, N = 8000L, T = 12L, p = 5L, cohort_law = "late_0_3_7_11", obs_law = "comp_arbitrary", truth_law = "late_0_3_7_11"),
    DATE_BASE = spec("DATE_BASE", 21L, "DATES", FALSE, obs_law = "date_base_35_65", derived = "rho_date_base"),
    DATE_HET = spec("DATE_HET", 22L, "DATES", FALSE, obs_law = "date_het"),
    DATE_COHORT_SLOPES = spec("DATE_COHORT_SLOPES", 23L, "DATES", FALSE, obs_law = "date_cohort_slopes"),
    REPRESENTATIVE_ROSTERS = spec("REPRESENTATIVE_ROSTERS", 24L, "TRANSPORT", FALSE, obs_law = "representative_rosters"),
    DATE_CONSTANT_EFFECT = spec("DATE_CONSTANT_EFFECT", 25L, "DATES", FALSE, effect_law = "constant_effect_no_x1", obs_law = "date_base_35_65", truth_law = "constant_effect_no_x1"),
    DATE_INTERACTION_Q = spec("DATE_INTERACTION_Q", 26L, "DATES", FALSE, obs_law = "date_interaction_q"),
    DATE_DISTINCT_0 = spec("DATE_DISTINCT_0", 27L, "DATES", FALSE, obs_law = "date_distinct"),
    DATE_DISTINCT_1 = spec("DATE_DISTINCT_1", 28L, "DATES", FALSE, untreated_law = "distinct_x5_shift", obs_law = "date_distinct"),
    DATE_WRONG_Q_J = spec("DATE_WRONG_Q_J", 29L, "DATES", FALSE, untreated_law = "wrong_qj_x1x5_shift", effect_law = "wrong_qj_x1x5_effect", obs_law = "date_interaction_q"),
    INNOVATION_SELECTION = spec("INNOVATION_SELECTION", 30L, "MEAN_MTAR", FALSE, obs_law = "innovation_selection"),
    CONT_PS_only = spec("CONT_PS_only", 31L, "TREATMENT_ROBUSTNESS", FALSE, xformla = xv, cohort_law = "cont_linear_v", untreated_law = "cont_quad_v", obs_law = "complete", truth_law = "cont_single_g3"),
    CONT_OR_only = spec("CONT_OR_only", 32L, "TREATMENT_ROBUSTNESS", FALSE, xformla = xv, cohort_law = "cont_quad_v", untreated_law = "cont_linear_v", obs_law = "complete", truth_law = "cont_single_g3"),
    CONT_both_wrong = spec("CONT_both_wrong", 33L, "TREATMENT_ROBUSTNESS", FALSE, xformla = xv, cohort_law = "cont_quad_v", untreated_law = "cont_quad_v", obs_law = "complete", truth_law = "cont_single_g3"),
    CONT_both_correct = spec("CONT_both_correct", 34L, "TREATMENT_ROBUSTNESS", FALSE, xformla = xv, cohort_law = "cont_linear_v", untreated_law = "cont_quad_v", obs_law = "complete", truth_law = "cont_single_g3"),
    EARLY_DRIFT = spec("EARLY_DRIFT", 35L, "PT", FALSE, cohort_law = "equiprob_0_4", untreated_law = "early_drift_g4", obs_law = "complete", truth_law = "early_drift_g4"),
    PT_POST_INVALID = spec("PT_POST_INVALID", 36L, "PT", FALSE, untreated_law = "pt_post_invalid", obs_law = "complete"),
    ANTICIPATION_ONE = spec("ANTICIPATION_ONE", 37L, "ANTICIPATION", FALSE, effect_law = "anticipation_one", obs_law = "complete", truth_law = "anticipation_one"),
    GAP_LABELS = spec("GAP_LABELS", 38L, "CALENDAR", FALSE, calendar = c(2000, 2002, 2007, 2010, 2020), cohort_law = "gap_labels_0345", untreated_law = "default", effect_law = "gap_index_effect", obs_law = "comp_arbitrary", truth_law = "gap_labels"),
    DATED_WEIGHTS = spec("DATED_WEIGHTS", 39L, "WEIGHTS_REFERENCE", FALSE, obs_law = "comp_arbitrary", derived = "dated_x2_and_w", truth_law = "dated_weights"),
    W_STATIC = spec("W_STATIC", 40L, "WEIGHTS_REFERENCE", FALSE, obs_law = "complete", derived = "static_w", truth_law = "w_static"),
    CLUSTER_N2000_C100 = spec("CLUSTER_N2000_C100", 41L, "CLUSTER", FALSE, N = 2000L, obs_law = "comp_arbitrary", het_law = "het_x1", cluster_C = 100L, derived = "cluster"),
    CLUSTER_N8000_C400 = spec("CLUSTER_N8000_C400", 42L, "CLUSTER", FALSE, N = 8000L, obs_law = "comp_arbitrary", het_law = "het_x1", cluster_C = 400L, derived = "cluster"),
    CLUSTER_N2000_C20 = spec("CLUSTER_N2000_C20", 43L, "CLUSTER", FALSE, N = 2000L, obs_law = "comp_arbitrary", het_law = "het_x1", cluster_C = 20L, derived = "cluster"),
    CLUSTER_N2000_C5 = spec("CLUSTER_N2000_C5", 44L, "CLUSTER", FALSE, N = 2000L, obs_law = "comp_arbitrary", het_law = "het_x1", cluster_C = 5L, derived = "cluster"),
    MOMENT_t3 = spec("MOMENT_t3", 45L, "ERROR_MOMENTS", FALSE, obs_law = "comp_arbitrary", error_law = "t3"),
    MOMENT_t1.5 = spec("MOMENT_t1.5", 46L, "ERROR_MOMENTS", FALSE, obs_law = "comp_arbitrary", error_law = "t1.5"),
    TEST_null = spec("TEST_null", 47L, "PRETEST", FALSE, cohort_law = "equiprob_0_3", obs_law = "complete", truth_law = "equiprob_0_3"),
    TEST_conditional_alternative = spec("TEST_conditional_alternative", 48L, "PRETEST", FALSE, cohort_law = "equiprob_0_3", untreated_law = "cond_alt_g3", obs_law = "complete", truth_law = "equiprob_0_3"),
    GMM_SPARSE = spec("GMM_SPARSE", 49L, "GMM_COVARIANCE", FALSE, N = 500L, T = 12L, obs_law = "comp_arbitrary"),
    GMM_ZERO_SCALE = spec("GMM_ZERO_SCALE", 50L, "GMM_COVARIANCE", FALSE, untreated_law = "zero_noise", effect_law = "theta_only", obs_law = "complete", error_law = "zero", truth_law = "theta_only"),
    BAL_ALIAS = spec("BAL_BASE", 5L, "FORMULA", FALSE, obs_law = "complete", derived = "xdup"),
    F_NO_INTERCEPT = spec("F_NO_INTERCEPT", 51L, "FORMULA", FALSE, untreated_law = "no_intercept_drift", obs_law = "complete"),
    DATE_NO_INTERCEPT = spec("DATE_NO_INTERCEPT", 52L, "DATES", FALSE, untreated_law = "no_intercept_drift", effect_law = "half_x1_only", obs_law = "date_base_35_65", truth_law = "half_x1_only"),
    MCAR_HALF = spec("MCAR_HALF", 53L, "DATES", FALSE, obs_law = "mcar_date_50"),
    E_TWO_DATES_A = spec("E_TWO_DATES_A", 54L, "REFUSAL", FALSE, N = 2000L, T = 2L, cohort_law = "equiprob_0_2", effect_law = "two_dates_g2", obs_law = "mcar_date_65", truth_law = "two_dates_g2"),
    E_OR_NEAR_COLLINEAR = spec("BAL_BASE", 5L, "REFUSAL", FALSE, obs_law = "complete", derived = "near_collinear_v"),
    DATED_MISSING = spec("BAL_BASE", 5L, "REFUSAL", FALSE, obs_law = "complete", derived = "dated_x2_missing_t1"),
    E_CHANGING_WEIGHT_NULL = spec("BAL_BASE", 5L, "REFUSAL", FALSE, obs_law = "complete", derived = "dated_w_only"),
    E_OVERLAP = spec("E_OVERLAP", 55L, "REFUSAL", FALSE, cohort_law = "overlap_x1_determines_g3", obs_law = "complete", truth_law = "overlap_g3"),
    E_SELECTION_SEPARATED = spec("E_SELECTION_SEPARATED", 56L, "REFUSAL", FALSE, obs_law = "selection_separated_x1"),
    E_PREDICTION_SUPPORT = spec("E_PREDICTION_SUPPORT", 57L, "REFUSAL", FALSE, cohort_law = "support_x3_determines_g3", obs_law = "complete", truth_law = "equiprob_0_3"),
    SUPPORT_DISCONNECTED = spec("SUPPORT_DISCONNECTED", 58L, "SUPPORT", FALSE, obs_law = "disconnected_12_45"),
    AR1_MCAR_DATE = spec("C_MCAR_DATE", 6L, "SERIAL_ERROR", TRUE, obs_law = "mcar_date_65", error_law = "ar1"),
    ATTRITION_MTAR = spec("ATTRITION_MTAR", 59L, "MONOTONE_ATTRITION", TRUE, obs_law = "attrition_mtar")
  )
}; RECIPE_SPECS <- make_recipe_specs(); get_recipe_spec <- function(recipe_id) {
  if (!is.character(recipe_id) || length(recipe_id) != 1L || is.na(recipe_id) ||
      !nzchar(recipe_id) || !(recipe_id %in% names(RECIPE_SPECS))) {
    stop("Unknown or unsupported recipe_id: ", as.character(recipe_id), " (refusing fallback to any base law).", call. = FALSE)
  }; RECIPE_SPECS[[recipe_id]]
}
# RNG stream index n = ((s - 1)*100 + (r - 1))*5 + p, p in 1..5:
# 1 = population, 2 = observation, 3 = FIT bootstrap, 4 = aggregation, 5 = pretest.
rng_stream_index <- function(stream_index, rep, purpose) {
  s <- as.integer(stream_index); r <- as.integer(rep); p <- as.integer(purpose)
  stopifnot(s >= 1L, s <= 59L, r >= 1L, r <= 100L, p >= 1L, p <= 5L); ((s - 1L) * 100L + (r - 1L)) * 5L + p
}; expand_rng_streams <- function(n_streams = 59L, n_reps = 100L, n_purposes = 5L, seed = 61012026L) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(as.integer(seed))
  total <- as.integer(n_streams) * as.integer(n_reps) * as.integer(n_purposes); states <- vector("list", total)
  cur <- get(".Random.seed", envir = .GlobalEnv); for (i in seq_len(total)) {
    cur <- parallel::nextRNGStream(cur); states[[i]] <- as.integer(cur)
  }; states
}; restore_rng_state <- function(state_vec) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); assign(".Random.seed", as.integer(state_vec), envir = .GlobalEnv); invisible(NULL)
}; invert_discrete_matrix <- function(u, cum_probs) {
  M <- ncol(cum_probs); if (M == 1L) return(rep(1L, length(u))); 1L + rowSums(u > cum_probs[, -M, drop = FALSE])
}; invert_discrete_vector <- function(u, probs) {
  cp <- cumsum(probs); cp[length(cp)] <- 1.0; findInterval(u, cp, rightmost.closed = TRUE, left.open = TRUE) + 1L
}; lex_edges <- function(T, kind = c("arbitrary", "adjacent")) {
  kind <- match.arg(kind); if (identical(kind, "adjacent")) {
    cbind(s = seq_len(T - 1L), t = seq_len(T - 1L) + 1L)
  } else {
    idx <- which(upper.tri(matrix(0, T, T)), arr.ind = TRUE); idx <- idx[order(idx[, 1L], idx[, 2L]), , drop = FALSE]
    cbind(s = as.integer(idx[, 1L]), t = as.integer(idx[, 2L]))
  }
}; generate_recipe_draw <- function(recipe_id, pop_state, obs_state) {
  sp <- get_recipe_spec(recipe_id)
  N <- sp$N; T <- sp$T; p <- sp$p; cal <- sp$calendar
  # 1. Population stream: covariates (1..p), U_G, alpha, innovations, cluster shocks
  restore_rng_state(pop_state); covs <- vector("list", p); cov_names <- paste0("X", seq_len(p))
  is_cont_v <- sp$cohort_law %in% c("cont_linear_v", "cont_quad_v"); if (is_cont_v) cov_names[2L] <- "V"; for (j in seq_len(p)) {
    if (j == 1L) {
      covs[[j]] <- ifelse(stats::runif(N) <= 0.5, 0.0, 1.0)
    } else if (j == 2L) {
      covs[[j]] <- if (is_cont_v) stats::runif(N, -1.0, 1.0) else stats::rnorm(N)
    } else if (j == 3L || j == 5L) {
      covs[[j]] <- ifelse(stats::runif(N) <= 0.5, -0.5, 0.5)
    } else if (j == 4L) {
      covs[[j]] <- stats::rexp(N, rate = 1.0) - 1.0
    } else {
      rem <- (j - 6L) %% 3L
      covs[[j]] <- if (rem == 0L) stats::rnorm(N) else if (rem == 1L) ifelse(stats::runif(N) <= 0.5, -0.5, 0.5) else stats::rexp(N, rate = 1.0) - 1.0
    }
  }; names(covs) <- cov_names; X1 <- covs[["X1"]]; X3 <- covs[["X3"]]; X5 <- covs[["X5"]]; u_G <- stats::runif(N)
  if (sp$cohort_law %in% c("default_0345", "late_0_3_7_11", "gap_labels_0345")) {
    codes_pos <- if (identical(sp$cohort_law, "late_0_3_7_11")) c(0L, 3L, 7L, 11L) else c(0L, 3L, 4L, 5L)
    q0 <- c(0.4, 0.7, 0.9, 1.0); q1 <- c(0.2, 0.4, 0.7, 1.0)
    cp <- ifelse(matrix(X1 == 1.0, N, 4L), matrix(q1, N, 4L, byrow = TRUE), matrix(q0, N, 4L, byrow = TRUE))
    g_rank <- invert_discrete_matrix(u_G, cp); G_pos <- codes_pos[g_rank]
    G_label <- if (identical(sp$cohort_law, "gap_labels_0345")) ifelse(G_pos == 0L, 0.0, cal[pmax(G_pos, 1L)]) else as.double(G_pos)
  } else if (identical(sp$cohort_law, "cont_linear_v") || identical(sp$cohort_law, "cont_quad_v")) {
    V <- covs[["V"]]; p3 <- if (identical(sp$cohort_law, "cont_linear_v")) expit_num(0.2 + 0.6 * V) else expit_num(0.2 + 0.6 * (V^2))
    G_pos <- ifelse(u_G <= 1.0 - p3, 0L, 3L); G_label <- as.double(G_pos); g_rank <- ifelse(G_pos == 0L, 1L, 2L)
  } else if (sp$cohort_law %in% c("equiprob_0_4", "equiprob_0_3", "equiprob_0_2")) {
    tr_code <- if (identical(sp$cohort_law, "equiprob_0_4")) 4L else if (identical(sp$cohort_law, "equiprob_0_3")) 3L else 2L
    G_pos <- ifelse(u_G <= 0.5, 0L, tr_code); G_label <- as.double(G_pos); g_rank <- ifelse(G_pos == 0L, 1L, 2L)
  } else if (identical(sp$cohort_law, "overlap_x1_determines_g3")) {
    G_pos <- ifelse(X1 == 1.0, 3L, 0L); G_label <- as.double(G_pos); g_rank <- ifelse(G_pos == 0L, 1L, 2L)
  } else if (identical(sp$cohort_law, "support_x3_determines_g3")) {
    G_pos <- ifelse(X3 > 0.0, 3L, 0L); G_label <- as.double(G_pos); g_rank <- ifelse(G_pos == 0L, 1L, 2L)
  } else {
    stop("Unhandled cohort_law: ", sp$cohort_law, call. = FALSE)
  }; alpha <- ifelse(stats::runif(N) <= 0.5, -1.0, 1.0); U_mat <- if (sp$error_law %in% c("normal", "ar1")) {
    matrix(stats::rnorm(N * T), nrow = N, ncol = T)
  } else if (identical(sp$error_law, "exp")) {
    matrix(stats::rexp(N * T, rate = 1.0) - 1.0, nrow = N, ncol = T)
  } else if (identical(sp$error_law, "t5")) {
    matrix(sqrt(3.0 / 5.0) * stats::rt(N * T, df = 5.0), nrow = N, ncol = T)
  } else if (identical(sp$error_law, "t3")) {
    matrix(sqrt(1.0 / 3.0) * stats::rt(N * T, df = 3.0), nrow = N, ncol = T)
  } else if (identical(sp$error_law, "t1.5")) {
    matrix(stats::rt(N * T, df = 1.5), nrow = N, ncol = T)
  } else if (identical(sp$error_law, "zero")) {
    matrix(0.0, nrow = N, ncol = T)
  } else {
    stop("Unhandled error_law: ", sp$error_law, call. = FALSE)
  }; h_vec <- if (identical(sp$het_law, "het_x1")) sqrt(1.0 + 0.5 * X1) else rep(1.0, N); id_vec <- seq_len(N)
  cluster_vec <- if (sp$cluster_C > 0L) ((id_vec - 1L) %% sp$cluster_C) + 1L else NULL; if (sp$cluster_C > 0L) {
    v_ct <- matrix(stats::rnorm(sp$cluster_C * T), nrow = sp$cluster_C, ncol = T)
    eps_mat <- 0.5 * h_vec * (sqrt(0.2) * v_ct[cluster_vec, , drop = FALSE] + sqrt(0.8) * U_mat)
  } else if (identical(sp$error_law, "ar1")) {
    eps_mat <- matrix(0.0, nrow = N, ncol = T); eps_mat[, 1L] <- 0.5 * h_vec * U_mat[, 1L]
    for (k in 2L:T) eps_mat[, k] <- 0.7 * eps_mat[, k - 1L] + 0.5 * h_vec * sqrt(0.51) * U_mat[, k]
  } else {
    eps_mat <- 0.5 * h_vec * U_mat
  }; bz <- 0.3 * X1; b_other <- 0.2 / sqrt(as.double(p)); for (j in setdiff(seq_len(p), c(1L, 5L))) bz <- bz + b_other * covs[[j]]
  Y_mat <- matrix(0.0, nrow = N, ncol = T); for (k in seq_len(T)) {
    drift <- if (identical(sp$untreated_law, "no_intercept_drift")) 0.0 else 0.1 * k; y0_k <- alpha + drift + k * bz + eps_mat[, k]
    if (identical(sp$untreated_law, "distinct_x5_shift")) {
      y0_k <- y0_k + 0.4 * k * X5
    } else if (identical(sp$untreated_law, "wrong_qj_x1x5_shift")) {
      y0_k <- y0_k + 0.6 * k * X1 * X5
    } else if (identical(sp$untreated_law, "cont_quad_v")) {
      y0_k <- y0_k + 0.6 * k * (covs[["V"]]^2)
    } else if (identical(sp$untreated_law, "cont_linear_v")) {
      y0_k <- y0_k + 0.6 * k * covs[["V"]]
    } else if (identical(sp$untreated_law, "early_drift_g4")) {
      y0_k <- y0_k + 0.4 * as.double(G_pos == 4L) * max(3L - k, 0L)
    } else if (identical(sp$untreated_law, "pt_post_invalid")) {
      y0_k <- y0_k + 0.3 * as.double(G_pos > 0L) * pmax(k - G_pos + 1L, 0L)
    } else if (identical(sp$untreated_law, "cond_alt_g3")) {
      y0_k <- y0_k + 0.4 * as.double(G_pos == 3L) * (2.0 * X1 - 1.0) * k
    } else if (!sp$untreated_law %in% c("default", "no_intercept_drift", "zero_noise")) {
      stop("Unhandled untreated_law: ", sp$untreated_law, call. = FALSE)
    }; treated_k <- (G_pos > 0L) & (k >= G_pos); theta_g <- 1.0 + 0.1 * (G_pos - 3.0)
    eff_k <- if (sp$effect_law %in% c("default", "gap_index_effect")) {
      ifelse(treated_k, theta_g + 0.25 * (k - G_pos) + 0.5 * X1, 0.0)
    } else if (identical(sp$effect_law, "constant_effect_no_x1")) {
      ifelse(treated_k, theta_g + 0.25 * (k - G_pos), 0.0)
    } else if (identical(sp$effect_law, "wrong_qj_x1x5_effect")) {
      ifelse(treated_k, theta_g + 0.25 * (k - G_pos) + 0.5 * X1 + 0.8 * X1 * X5, 0.0)
    } else if (identical(sp$effect_law, "anticipation_one")) {
      ifelse(treated_k, theta_g + 0.25 * (k - G_pos) + 0.5 * X1, ifelse(G_pos > 0L & k == (G_pos - 1L), 0.3, 0.0))
    } else if (identical(sp$effect_law, "theta_only")) {
      ifelse(treated_k, theta_g, 0.0)
    } else if (identical(sp$effect_law, "half_x1_only")) {
      ifelse(treated_k, 0.5 * X1, 0.0)
    } else if (identical(sp$effect_law, "two_dates_g2")) {
      ifelse(G_pos == 2L & k >= 2L, 1.0 + 0.5 * X1, 0.0)
    } else {
      stop("Unhandled effect_law: ", sp$effect_law, call. = FALSE)
    }; Y_mat[, k] <- y0_k + eff_k
  }
  # 2. Observation stream: draw masks in listed order
  restore_rng_state(obs_state); obs_mat <- matrix(FALSE, nrow = N, ncol = T); lam_4 <- c(0.0, 0.6, 0.4, 0.2)
  if (sp$obs_law %in% c("comp_adjacent", "comp_arbitrary")) {
    edges <- lex_edges(T, if (identical(sp$obs_law, "comp_adjacent")) "adjacent" else "arbitrary")
    K <- nrow(edges); h_e <- (T + 1.0 - edges[, "s"] - edges[, "t"]) / (T - 2.0); probs <- (1.0 + outer(lam_4[g_rank] * alpha, h_e)) / K
    cp <- t(apply(probs, 1L, cumsum)); cp[, K] <- 1.0; chosen <- invert_discrete_matrix(stats::runif(N), cp)
    obs_mat[cbind(id_vec, edges[chosen, "s"])] <- TRUE; obs_mat[cbind(id_vec, edges[chosen, "t"])] <- TRUE
  } else if (identical(sp$obs_law, "cc_uniform") || identical(sp$obs_law, "cc_comp")) {
    edges <- lex_edges(T, "arbitrary"); K <- nrow(edges); u_c <- stats::runif(N); u_e <- stats::runif(N)
    is_cc <- u_c <= ifelse(X1 == 1.0, 0.6, 0.2); if (identical(sp$obs_law, "cc_uniform")) {
      chosen <- invert_discrete_vector(u_e, rep(1.0 / K, K))
    } else {
      h_e <- (T + 1.0 - edges[, "s"] - edges[, "t"]) / (T - 2.0); probs <- (1.0 + outer(lam_4[g_rank] * alpha, h_e)) / K
      cp <- t(apply(probs, 1L, cumsum)); cp[, K] <- 1.0; chosen <- invert_discrete_matrix(u_e, cp)
    }; obs_mat[is_cc, ] <- TRUE; not_cc <- which(!is_cc)
    obs_mat[cbind(not_cc, edges[chosen[not_cc], "s"])] <- TRUE; obs_mat[cbind(not_cc, edges[chosen[not_cc], "t"])] <- TRUE
  } else if (identical(sp$obs_law, "complete")) {
    obs_mat[,] <- TRUE
  } else if (identical(sp$obs_law, "mcar_date_65") || identical(sp$obs_law, "mcar_date_50")) {
    rate <- if (identical(sp$obs_law, "mcar_date_65")) 0.65 else 0.50; obs_mat <- matrix(stats::runif(N * T), nrow = N, ncol = T) <= rate
  } else if (identical(sp$obs_law, "mcar_pair")) {
    edges <- lex_edges(T, "arbitrary"); K <- nrow(edges); chosen <- invert_discrete_vector(stats::runif(N), rep(1.0 / K, K))
    obs_mat[cbind(id_vec, edges[chosen, "s"])] <- TRUE; obs_mat[cbind(id_vec, edges[chosen, "t"])] <- TRUE
  } else if (identical(sp$obs_law, "date_base_35_65")) {
    U_obs <- matrix(stats::runif(N * T), nrow = N, ncol = T); obs_mat <- U_obs <= matrix(ifelse(X1 == 1.0, 0.65, 0.35), nrow = N, ncol = T)
  } else if (identical(sp$obs_law, "date_het")) {
    U_obs <- matrix(stats::runif(N * T), nrow = N, ncol = T)
    obs_mat <- U_obs <= expit_num(outer(0.6 * X1, c(-1.0, -0.6, -0.2, 0.2, 0.6), "+"))
  } else if (identical(sp$obs_law, "date_cohort_slopes")) {
    U_obs <- matrix(stats::runif(N * T), nrow = N, ncol = T); r_i <- expit_num(-0.5 + c(0.2, 0.6, 1.0, 0.8)[g_rank] * X1)
    obs_mat <- U_obs <= matrix(r_i, nrow = N, ncol = T)
  } else if (identical(sp$obs_law, "representative_rosters")) {
    edges <- lex_edges(T, "arbitrary"); K <- nrow(edges); u_type <- stats::runif(N); u_sub <- stats::runif(N)
    htype <- invert_discrete_vector(u_type, c(0.2, 0.4, 0.2, 0.2)); obs_mat[htype == 1L, ] <- TRUE; pair_idx <- which(htype == 2L)
    if (length(pair_idx)) {
      chosen_e <- invert_discrete_vector(u_sub[pair_idx], rep(1.0 / K, K))
      obs_mat[cbind(pair_idx, edges[chosen_e, "s"])] <- TRUE; obs_mat[cbind(pair_idx, edges[chosen_e, "t"])] <- TRUE
    }; sing_idx <- which(htype == 3L); if (length(sing_idx)) {
      chosen_d <- invert_discrete_vector(u_sub[sing_idx], rep(1.0 / T, T)); obs_mat[cbind(sing_idx, chosen_d)] <- TRUE
    }
  } else if (identical(sp$obs_law, "date_interaction_q") || identical(sp$obs_law, "date_distinct")) {
    U_obs <- matrix(stats::runif(N * T), nrow = N, ncol = T)
    inter <- if (identical(sp$obs_law, "date_interaction_q")) 1.2 * X1 * X5 else 0.0
    obs_mat <- U_obs <= matrix(expit_num(-0.5 + 0.6 * X1 + 0.8 * X5 + inter), nrow = N, ncol = T)
  } else if (identical(sp$obs_law, "innovation_selection")) {
    edges <- lex_edges(T, "arbitrary"); K <- nrow(edges); probs <- matrix(1.0 / K, nrow = N, ncol = K)
    e23 <- which(edges[, "s"] == 2L & edges[, "t"] == 3L); e45 <- which(edges[, "s"] == 4L & edges[, "t"] == 5L)
    g3_idx <- which(G_pos == 3L); sgn <- sign(eps_mat[g3_idx, 3L] - eps_mat[g3_idx, 2L])
    probs[g3_idx, e23] <- (1.0 + 0.6 * sgn) / K; probs[g3_idx, e45] <- (1.0 - 0.6 * sgn) / K
    cp <- t(apply(probs, 1L, cumsum)); cp[, K] <- 1.0; chosen <- invert_discrete_matrix(stats::runif(N), cp)
    obs_mat[cbind(id_vec, edges[chosen, "s"])] <- TRUE; obs_mat[cbind(id_vec, edges[chosen, "t"])] <- TRUE
  } else if (identical(sp$obs_law, "selection_separated_x1")) {
    obs_mat <- matrix(X1 == 1.0, nrow = N, ncol = T)
  } else if (identical(sp$obs_law, "disconnected_12_45")) {
    first_pair <- stats::runif(N) <= 0.5; obs_mat[first_pair, c(1L, 2L)] <- TRUE; obs_mat[!first_pair, c(4L, 5L)] <- TRUE
  } else if (identical(sp$obs_law, "attrition_mtar")) {
    u_att <- stats::runif(N); obs_mat[, 1L] <- TRUE; a_att <- c(1.5, 1.0, 0.5, 0.0); gam_i <- lam_4[g_rank]
    for (k in 2L:T) obs_mat[, k] <- u_att <= expit_num(a_att[k - 1L] + 0.6 * X1 + gam_i * alpha)
  } else {
    stop("Unhandled obs_law: ", sp$obs_law, call. = FALSE)
  }
  # 3. Build PANEL (N * T rows, unit-major) and F (N rows)
  id_long <- rep(id_vec, each = T); k_long <- rep(seq_len(T), times = N); t_long <- cal[k_long]; G_long <- G_label[id_long]
  y_long <- as.double(t(ifelse(obs_mat, Y_mat, NA_real_)))
  panel_list <- list(id = as.integer(id_long), t = as.double(t_long), G = as.double(G_long))
  frame_list <- list(id = as.integer(id_vec), G = as.double(G_label)); for (nm in cov_names) {
    if (nm == "X2" && sp$derived %in% c("dated_x2_and_w", "dated_x2_missing_t1")) {
      x2_dated <- covs[["X2"]][id_long] + 0.1 * t_long; if (identical(sp$derived, "dated_x2_missing_t1")) x2_dated[k_long == 1L] <- NA_real_
      panel_list[["X2"]] <- as.double(x2_dated)
    } else {
      panel_list[[nm]] <- as.double(covs[[nm]][id_long]); frame_list[[nm]] <- as.double(covs[[nm]])
    }
  }; if (sp$derived %in% c("dated_x2_and_w", "dated_w_only")) {
    panel_list[["w"]] <- as.double(1.0 + 0.1 * t_long * X1[id_long])
  } else if (identical(sp$derived, "static_w")) {
    w_stat <- as.double(1.0 + 0.5 * X1); panel_list[["w"]] <- w_stat[id_long]; frame_list[["w"]] <- w_stat
  } else if (identical(sp$derived, "cluster")) {
    panel_list[["cluster"]] <- as.integer(cluster_vec[id_long]); frame_list[["cluster"]] <- as.integer(cluster_vec)
  } else if (identical(sp$derived, "xdup")) {
    xdup <- as.double(2.0 * covs[["X2"]]); panel_list[["Xdup"]] <- xdup[id_long]; frame_list[["Xdup"]] <- xdup
  } else if (identical(sp$derived, "near_collinear_v")) {
    v_col <- as.double(covs[["X2"]] + 1e-9 * covs[["X4"]]); panel_list[["V"]] <- v_col[id_long]; frame_list[["V"]] <- v_col
  } else if (identical(sp$derived, "rho_date_base")) {
    r0 <- 0.35; r1 <- 0.65; uB0 <- 1.0 - (1.0 - r0)^T; uB1 <- 1.0 - (1.0 - r1)^T
    uA0 <- uB0 - T * r0 * (1.0 - r0)^(T - 1L); uA1 <- uB1 - T * r1 * (1.0 - r1)^(T - 1L); mu_pop <- c(1.0 / 3.0, 0.4, 0.6, 0.75)[g_rank]
    uA_x <- ifelse(X1 == 1.0, uA1, uA0); uB_x <- ifelse(X1 == 1.0, uB1, uB0)
    frame_list[["rho_A"]] <- as.double(uA_x / ((1.0 - mu_pop) * uA0 + mu_pop * uA1))
    frame_list[["rho_B"]] <- as.double(uB_x / ((1.0 - mu_pop) * uB0 + mu_pop * uB1)); frame_list[["rho_one"]] <- rep(1.0, N)
  }; panel_list[["y"]] <- y_long; PANEL <- as.data.frame(panel_list, stringsAsFactors = FALSE)
  F_df <- as.data.frame(frame_list, stringsAsFactors = FALSE); list(
    PANEL = PANEL, F = F_df, meta = list(
      recipe_id = recipe_id, stream_key = sp$stream_key, stream_index = sp$stream_index, n_units = N, n_dates = T, n_rows = nrow(PANEL),
      n_observed_rows = sum(!is.na(PANEL$y)), n_complete_units = sum(rowSums(obs_mat) == T),
      n_pair_eligible_units = sum(rowSums(obs_mat) >= 2L)
    )
  )
}; recipe_population_masses <- function(recipe_id, weightsname = NULL, fix_weights = NULL) {
  sp <- get_recipe_spec(recipe_id)
  weighted <- !is.null(weightsname) && nzchar(weightsname); if (sp$cohort_law %in% c("default_0345", "late_0_3_7_11", "gap_labels_0345")) {
    codes_pos <- if (identical(sp$cohort_law, "late_0_3_7_11")) c(0L, 3L, 7L, 11L) else c(0L, 3L, 4L, 5L)
    labels <- if (identical(sp$cohort_law, "gap_labels_0345")) c(0.0, sp$calendar[codes_pos[-1L]]) else as.double(codes_pos)
    pi_base <- c(0.30, 0.25, 0.25, 0.20); mu_base <- c(1.0 / 3.0, 0.4, 0.6, 0.75)
    pi_raw <- if (weighted && identical(recipe_id, "W_STATIC")) {
      pi_base * (1.0 + 0.5 * mu_base)
    } else if (weighted && identical(recipe_id, "DATED_WEIGHTS")) {
      pi_base * (1.0 + 0.1 * mu_base)
    } else {
      pi_base
    }; list(codes_pos = codes_pos, labels = labels, pi_raw = pi_raw, mu_base = mu_base)
  } else if (sp$cohort_law %in% c("cont_linear_v", "cont_quad_v", "equiprob_0_3", "support_x3_determines_g3")) {
    p3 <- if (sp$cohort_law %in% c("cont_linear_v", "cont_quad_v")) {
      stats::integrate(function(v) expit_num(0.2 + 0.6 * if (sp$cohort_law == "cont_quad_v") v^2 else v),
                       -1, 1, rel.tol = 1e-12)$value / 2
    } else 0.5
    list(codes_pos = c(0L, 3L), labels = c(0.0, 3.0), pi_raw = c(1 - p3, p3), mu_base = c(0.5, 0.5))
  } else if (identical(sp$cohort_law, "overlap_x1_determines_g3")) {
    list(codes_pos = c(0L, 3L), labels = c(0.0, 3.0), pi_raw = c(0.5, 0.5), mu_base = c(0.0, 1.0))
  } else if (identical(sp$cohort_law, "equiprob_0_4")) {
    list(codes_pos = c(0L, 4L), labels = c(0.0, 4.0), pi_raw = c(0.5, 0.5), mu_base = c(0.5, 0.5))
  } else if (identical(sp$cohort_law, "equiprob_0_2")) {
    list(codes_pos = c(0L, 2L), labels = c(0.0, 2.0), pi_raw = c(0.5, 0.5), mu_base = c(0.5, 0.5))
  } else {
    stop("Unhandled cohort_law in recipe_population_masses: ", sp$cohort_law, call. = FALSE)
  }
}
compute_cell_att_truth <- function(recipe_id, g_pos, t_pos, mu_g, anticipation = 0L, weightsname = NULL, fix_weights = NULL, cal = NULL) {
  sp <- get_recipe_spec(recipe_id)
  if (is.null(cal)) cal <- sp$calendar; weighted <- !is.null(weightsname) && nzchar(weightsname); theta_g <- 1.0 + 0.1 * (g_pos - 3.0)
  if (t_pos < g_pos) {
    if (identical(recipe_id, "ANTICIPATION_ONE") && as.integer(anticipation) == 1L && t_pos == (g_pos - 1L)) return(0.3); return(0.0)
  }; mu_eff <- if (weighted && identical(recipe_id, "W_STATIC")) {
    1.5 * mu_g / (1.0 + 0.5 * mu_g)
  } else if (weighted && identical(recipe_id, "DATED_WEIGHTS")) {
    anchor <- if (identical(fix_weights, "base_period")) as.double(g_pos - 1L - as.integer(anticipation)) else 1.0
    mu_g * (1.0 + 0.1 * anchor) / (1.0 + 0.1 * anchor * mu_g)
  } else {
    mu_g
  }; if (identical(sp$truth_law, "constant_effect_no_x1")) theta_g + 0.25 * (t_pos - g_pos)
  else if (identical(sp$truth_law, "theta_only")) theta_g
  else if (identical(sp$truth_law, "half_x1_only")) 0.5 * mu_eff
  else if (identical(sp$truth_law, "two_dates_g2")) 1.0 + 0.5 * mu_eff
  else theta_g + 0.25 * (t_pos - g_pos) + 0.5 * mu_eff
}; compute_aggte_truth <- function(gt_df, masses_df, calendar, type = c("simple", "group", "dynamic", "calendar"),
                                balance_e = NULL, min_e = -Inf, max_e = Inf) {
  type <- match.arg(type); post <- gt_df[gt_df$q_t >= gt_df$q_g, , drop = FALSE]
  w_map <- stats::setNames(masses_df$pi_raw, as.character(masses_df$label)); if (identical(type, "simple")) {
    sub <- post[post$q_t <= (post$q_g + max_e), , drop = FALSE]; w <- w_map[as.character(sub$group)]
    list(overall = sum(w * sub$primary_truth) / sum(w), egt = numeric(0), att_egt = numeric(0))
  } else if (identical(type, "group")) {
    g_list <- sort(unique(post$group))
    att_g <- vapply(g_list, function(g) mean(post$primary_truth[post$group == g & post$q_t <= (post$q_g + max_e)]), numeric(1L))
    wg <- w_map[as.character(g_list)]; list(overall = sum(wg * att_g) / sum(wg), egt = g_list, att_egt = att_g)
  } else if (identical(type, "dynamic")) {
    d1 <- min(calendar); dT <- max(calendar); df <- gt_df
    df$val <- ifelse(!is.na(df$primary_truth), df$primary_truth, df$auxiliary_truth)
    df <- df[!is.na(df$val), , drop = FALSE]; df$e <- df$t - df$group; if (!is.null(balance_e) && is.finite(balance_e)) {
      B <- as.double(balance_e); df <- df[(dT - df$group) >= B & df$e >= (B - dT + d1) & df$e <= B, , drop = FALSE]
    }; df <- df[df$e >= min_e & df$e <= max_e, , drop = FALSE]; e_list <- sort(unique(df$e)); att_e <- vapply(e_list, function(ev) {
      sub <- df[df$e == ev, , drop = FALSE]; w <- w_map[as.character(sub$group)]; sum(w * sub$val) / sum(w)
    }, numeric(1L)); nonneg <- e_list >= 0
    list(overall = if (any(nonneg)) mean(att_e[nonneg]) else NA_real_, egt = e_list, att_egt = att_e)
  } else {
    t_list <- sort(unique(post$t)); att_t <- vapply(t_list, function(tv) {
      sub <- post[post$t == tv, , drop = FALSE]; w <- w_map[as.character(sub$group)]; sum(w * sub$primary_truth) / sum(w)
    }, numeric(1L)); list(overall = mean(att_t), egt = t_list, att_egt = att_t)
  }
}
