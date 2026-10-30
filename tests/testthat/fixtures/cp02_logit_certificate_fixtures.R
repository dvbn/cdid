# Deterministic reference fixtures for private CP02 raw weighted-logit root certificate tests.
# All fixtures are constructed from exact dyadic or small-integer rationals without
# calling any optimizer, glm(), or random number generator.

make_logit_cert_default_control <- function(
  taylor_degree = 18L,
  max_abs_eta = 700.0,
  max_squarings = 11L,
  max_chunk_rows = 128L,
  max_matrix_entries = 100000L,
  max_work_ops = 10000000L
) {
  logit_certificate_control(
    taylor_degree = taylor_degree,
    max_abs_eta = max_abs_eta,
    max_squarings = max_squarings,
    max_chunk_rows = max_chunk_rows,
    max_matrix_entries = max_matrix_entries,
    max_work_ops = max_work_ops
  )
}

# 1. Exact rational beta0 = 0 problem with non-dyadic cohort rates r_j = Pi_j / M_j
#    and outside-risk + pair non-observer units.
#    n = 8 eligible A-units:
#      Cohort 2 (target, D = 1):
#        unit 1: w = 1,   z = (1, 0), pair_observed = TRUE
#        unit 2: w = 2,   z = (1, 1), pair_observed = TRUE
#        unit 3: w = 4,   z = (1, 9), pair_observed = FALSE  (non-observer in R)
#        -> M_1 = 7, Pi_1 = 3, r_1 = 3/7, c_1 = 7/3, c_2 = 14/3
#      Cohort 5 (control 1, D = 0):
#        unit 4: w = 2,   z = (1, 0), pair_observed = TRUE
#        unit 5: w = 1,   z = (1, 1), pair_observed = TRUE
#        unit 6: w = 3,   z = (1, -5), pair_observed = FALSE (non-observer in R)
#        -> M_2 = 6, Pi_2 = 3, r_2 = 1/2, c_4 = 4, c_5 = 2
#      Cohort Inf (control 2, D = 0):
#        unit 7: w = 3,   z = (1, 1), pair_observed = TRUE
#        -> M_3 = 3, Pi_3 = 3, r_3 = 1, c_7 = 3
#      Outside-risk unit (G = 9 not in cohorts = c(2, 5, Inf)):
#        unit 8: w = 5,   z = (1, 100), pair_observed = TRUE (outside-risk observer)
#
#    At beta0 = (0, 0):
#      eta_i = 0, o_i = 1, p_i = 1/2, 1 - p_i = 1/2, v_i = 1/4.
#      Sum of c_i over O:
#        target (D=1): c_1 + c_2 = 7/3 + 14/3 = 7
#        control (D=0): c_4 + c_5 + c_7 = 4 + 2 + 3 = 9
#      Exact score g(0) = (1/8) * sum_{i in O} c_i (1/2 - D_i) z_i:
#        g_1(0) = (1/8) * (7 * (-1/2) + 9 * (1/2)) = (1/8) * 1 = 1/8 = 0.125
#        g_2(0) = (1/8) * (c_2 * (-1/2) * 1 + (c_5 + c_7) * (1/2) * 1)
#               = (1/8) * ((14/3)*(-1/2) + 5*(1/2)) = (1/8) * (1/6) = 1/48
#      Exact Hessian H(0) = (1/8) * (1/4) * sum_{i in O} c_i z_i z_i':
#        sum_O c_i = 16       -> H_{11}(0) = 16 / 32 = 1/2 = 0.5
#        sum_O c_i x_i = 14/3 + 2 + 3 = 29/3 -> H_{12}(0) = H_{21}(0) = H_{22}(0) = 29 / 96
make_exact_beta0_zero_fixture <- function() {
  Z <- rbind(
    c(1.0, 0.0),
    c(1.0, 1.0),
    c(1.0, 9.0),
    c(1.0, 0.0),
    c(1.0, 1.0),
    c(1.0, -5.0),
    c(1.0, 1.0),
    c(1.0, 100.0)
  )
  colnames(Z) <- c("intercept", "x")
  w <- c(1.0, 2.0, 4.0, 2.0, 1.0, 3.0, 3.0, 5.0)
  G <- c(2.0, 2.0, 2.0, 5.0, 5.0, 5.0, Inf, 9.0)
  pair_observed <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
  cohorts <- c(2.0, 5.0, Inf)
  beta0 <- c(0.0, 0.0)
  # Upper-triangular F with exact M = F'F:
  # F = [[1, 1/2], [0, 1/2]] -> M = F'F = [[1, 1/2], [1/2, 1/2]]
  F_mat <- matrix(c(1.0, 0.0, 0.5, 0.5), nrow = 2L, ncol = 2L)

  list(
    Z = Z,
    w = w,
    G = G,
    pair_observed = pair_observed,
    cohorts = cohorts,
    beta0 = beta0,
    F = F_mat,
    r = 0.5,
    exact = list(
      n = 8.0,
      M_num = c(7, 6, 3),
      M_den = c(1, 1, 1),
      Pi_num = c(3, 3, 3),
      Pi_den = c(1, 1, 1),
      r_num = c(3, 1, 1),
      r_den = c(7, 2, 1),
      q_num = c(7, 3, 3),
      q_den = c(8, 4, 8),
      g_num = c(1, 1),
      g_den = c(8, 48),
      H_num = matrix(c(1, 29, 29, 29), nrow = 2L, ncol = 2L),
      H_den = matrix(c(2, 96, 96, 96), nrow = 2L, ncol = 2L)
    )
  )
}

# 2. The 16-unit nuisance-A fixture from docs/theory/cp02_ipw_gmm_numerical_propagation.tex
#    (Section 6 / ipw_gmm_propagation_reference.py) with stored root
#    beta0 = c(-0.725796364095429, 0.12361921796210888).
make_propagation_16unit_fixture <- function() {
  G <- c(2, 2, 2, 2, 2, 5, 5, 5, 5, 5, Inf, Inf, Inf, Inf, 9, 9)
  w <- c(1.0, 1.5, 1.0, 0.75, 1.25, 1.0, 2.0, 1.0, 0.5, 1.0, 1.0, 1.25, 1.0, 0.75, 1.0, 2.0)
  x <- c(0, 1, 2, 1, 3, 0, 1, 2, 1, 3, 0, 2, 1, 3, 1, 0)
  pair_observed <- c(
    TRUE, TRUE, TRUE, FALSE, TRUE,
    TRUE, TRUE, FALSE, TRUE, TRUE,
    TRUE, TRUE, FALSE, TRUE,
    TRUE, FALSE
  )
  Z <- cbind(intercept = rep(1.0, 16L), x = as.double(x))
  cohorts <- c(2.0, 5.0, Inf)
  # Stored beta0 from ipw_gmm_propagation_reference.py (rounded high-precision root for nuisance A)
  beta0 <- c(-0.725796364095429, 0.12361921796210888)
  # Upper-triangular factor F from ipw_gmm_propagation_reference.py (30-bit dyadic Cholesky of H(beta0))
  # H(beta0) at beta0 is [[0.20073486, 0.27885915], [0.27885915, 0.67432350]]
  H_mid <- matrix(c(0.2007348617, 0.2788591507, 0.2788591507, 0.6743235022), 2L, 2L)
  F_mat <- round(chol(H_mid) * 2^30) / 2^30
  F_mat[2L, 1L] <- 0.0

  list(
    Z = Z,
    w = w,
    G = G,
    pair_observed = pair_observed,
    cohorts = cohorts,
    beta0 = beta0,
    F = F_mat,
    r = 0.01,
    beta_star_hp = c(-0.725796364095429, 0.12361921796210888)
  )
}

# 3. Overlap fixture OV from docs/theory/cp02_finite_logit_validation.tex (Section 3, Example)
#    embedded as an 8-unit raw A-unit dataset (6 observers + 2 outside-risk units so n = 8,
#    with r_j = 1 so c_i = w_i matches OV's exact c_i).
#    At stored odds u0 = (11/8, 9/8), beta0 = (log(11/8), log(9/8)).
make_finite_logit_ov_fixture <- function() {
  Z <- rbind(
    c(1.0, 0.0),
    c(1.0, 0.0),
    c(1.0, 1.0),
    c(1.0, 1.0),
    c(1.0, 2.0),
    c(1.0, 2.0),
    c(1.0, 0.0),
    c(1.0, 1.0)
  )
  colnames(Z) <- c("intercept", "x")
  w <- c(1.0, 2.0, 1.5, 1.0, 1.0, 2.5, 1.0, 1.0)
  G <- c(5.0, 2.0, 5.0, 2.0, 5.0, 2.0, 9.0, 9.0)
  pair_observed <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  cohorts <- c(2.0, 5.0)
  beta0 <- c(log(11 / 8), log(9 / 8))
  F_mat <- chol(matrix(c(972 / 4096, 957 / 4096, 957 / 4096, 1602 / 4096), 2L, 2L))
  F_mat <- round(F_mat * 4096) / 4096
  F_mat[2L, 1L] <- 0.0

  list(
    Z = Z,
    w = w,
    G = G,
    pair_observed = pair_observed,
    cohorts = cohorts,
    beta0 = beta0,
    F = F_mat,
    r = (34 / 256) / 2.66015625,
    beta_star_hp = c(0.28969593300018094, 0.14197023953543122)
  )
}

# 4. Unsaturated 3-cell fixture with non-zero cell score residuals (1, -2, 1) at exact root
#    beta* = (0, log(2)) (odds = (1, 2, 4) at x = 0, 1, 2).
make_unsaturated_3cell_fixture <- function() {
  Z <- rbind(
    c(1.0, 0.0),
    c(1.0, 1.0),
    c(1.0, 2.0),
    c(1.0, 0.0),
    c(1.0, 0.0),
    c(1.0, 1.0),
    c(1.0, 2.0),
    c(1.0, 3.0)
  )
  colnames(Z) <- c("intercept", "x")
  w <- c(8.0, 9.0, 12.5, 29.5, 14.0, 12.0, 5.0, 4.0)
  G <- c(2.0, 2.0, 2.0, 2.0, 5.0, 5.0, 5.0, 99.0)
  pair_observed <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
  cohorts <- c(2.0, 5.0)
  beta0 <- c(0.0, log(2.0))
  H_exact <- matrix(c(569 / 240, 61 / 30, 61 / 30, 97 / 30), 2L, 2L)
  F_mat <- chol(H_exact)
  F_mat[2L, 1L] <- 0.0

  list(
    Z = Z,
    w = w,
    G = G,
    pair_observed = pair_observed,
    cohorts = cohorts,
    beta0 = beta0,
    F = F_mat,
    r = 0.01,
    H_exact = H_exact
  )
}

# 5. Quasi-separated fixture QS from docs/theory/cp02_finite_logit_validation.tex:
#    No X = 1 control (rows: (1,0) D=0, (1,0) D=1, (1,1) D=1, (1,1) D=1).
make_quasi_separated_fixture <- function(slope = log(2^20)) {
  Z <- rbind(
    c(1.0, 0.0),
    c(1.0, 0.0),
    c(1.0, 1.0),
    c(1.0, 1.0)
  )
  colnames(Z) <- c("intercept", "x")
  w <- c(1.0, 1.0, 1.0, 1.0)
  G <- c(5.0, 2.0, 2.0, 2.0)
  pair_observed <- c(TRUE, TRUE, TRUE, TRUE)
  cohorts <- c(2.0, 5.0)
  beta0 <- c(0.0, as.double(slope))
  p_x1 <- 1.0 / (1.0 + exp(-slope))
  v_x1 <- p_x1 * (1.0 - p_x1)
  H_approx <- matrix(c(
    (2 * 0.25 + 2 * v_x1) / 4.0,
    (2 * v_x1) / 4.0,
    (2 * v_x1) / 4.0,
    (2 * v_x1) / 4.0
  ), 2L, 2L)
  F_mat <- chol(H_approx)
  F_mat[2L, 1L] <- 0.0

  list(
    Z = Z,
    w = w,
    G = G,
    pair_observed = pair_observed,
    cohorts = cohorts,
    beta0 = beta0,
    F = F_mat,
    r = 0.1
  )
}

# 6. Raw weights vs cached rounded w/r discrepancy fixtures:
#    (a) Small integer normal weights (Pi_1 = 5, M_1 = 29):
#        unit 1 (observer, G=2):     w_1 = 5
#        unit 2 (non-observer, G=2): w_2 = 24
#        -> Exact M_1 = 29, Pi_1 = 5, exact c_1 = 5 * 29 / 5 = 29.
#        -> Cached rounded r_1_fl = fl(5 / 29), and fl(5 / r_1_fl) = 29 - 2^-48 < 29!
#        unit 3 (observer, G=5):     w_3 = 29 - 2^-48
#        -> Surrogate cached-w/r score at beta0 = 0 is g_cached(0) = 0.0 (falsely claiming
#           beta0 = 0 is the exact root), whereas exact raw score is g_raw(0) = -2^-49 / 3 < 0
#           and true root is beta*_raw = log(29 / (29 - 2^-48)) > 0.
#    (b) Subnormal cached rate precision loss (Pi_1 = 3 * 2^-540, M_1 = 3 * 2^-540 + 7 * 2^500):
#        Cached r_1_fl = fl(Pi_1 / M_1) underflows into subnormal range (~ (3/7)*2^-1040,
#        losing 18 bits = 262,144 ULPs of precision), so cached c_1_fl = fl(w_1 / r_1_fl)
#        and g_cached(0) lie strictly outside the raw-weight outward enclosures c_iv and g_iv!
make_raw_vs_cached_wr_fixture <- function() {
  w3 <- 29.0 - 2^(-48)
  Z <- matrix(c(1.0, 1.0, 1.0), nrow = 3L, ncol = 1L, dimnames = list(NULL, "intercept"))
  w <- c(5.0, 24.0, w3)
  G <- c(2.0, 2.0, 5.0)
  pair_observed <- c(TRUE, FALSE, TRUE)
  cohorts <- c(2.0, 5.0)
  beta0 <- c(0.0)
  F_mat <- matrix(sqrt((58.0 - 2^(-48)) / 12.0), nrow = 1L, ncol = 1L)

  r1_cached <- 5.0 / 29.0
  c1_cached <- 5.0 / r1_cached
  g0_cached <- (c1_cached * (-0.5) + w3 * 0.5) / 3.0
  g0_exact_num <- -2^(-49)
  g0_exact <- g0_exact_num / 3.0
  beta_star_raw <- log(29.0 / w3)

  # Subnormal rate case (3 units at O(1) mass scale):
  w_sub <- c(3.0 * 2^(-1040), 7.0, 7.0)
  r1_sub_cached <- w_sub[1L] / (w_sub[1L] + w_sub[2L])
  c1_sub_cached <- w_sub[1L] / r1_sub_cached
  g0_sub_cached <- (c1_sub_cached * (-0.5) + w_sub[3L] * 0.5) / 3.0
  F_sub <- matrix(sqrt(14.0 / 12.0), nrow = 1L, ncol = 1L)

  list(
    Z = Z,
    w = w,
    G = G,
    pair_observed = pair_observed,
    cohorts = cohorts,
    beta0 = beta0,
    F = F_mat,
    r1_cached = r1_cached,
    c1_cached = c1_cached,
    c1_exact = 29.0,
    w3 = w3,
    g0_cached = g0_cached,
    g0_exact = g0_exact,
    beta_star_raw = beta_star_raw,
    subnormal_rate = list(
      Z = Z,
      w = w_sub,
      G = G,
      pair_observed = pair_observed,
      cohorts = cohorts,
      beta0 = beta0,
      F = F_sub,
      r1_cached = r1_sub_cached,
      c1_cached = c1_sub_cached,
      g0_cached = g0_sub_cached
    )
  )
}
