# Tests for R/static_observation_moments.R (private static observation kernel).
# Expected values come from hand arithmetic (written in comments) or from the
# independent reference in reference/static_observation_reference.R, which
# never calls production code. Refit finite-difference fixtures validate the
# complete original-unit IF (direct term + every keyed nuisance loading).

# Production functions are NOT sourced here. They are taken from the
# environment this file runs in: under the package test runner that is a clone
# of the package namespace (or an isolated test environment holding the frozen
# helpers and the candidate module).
# `prod` is the closure environment of the kernel (the namespace, or the
# runner's isolated environment). Only the independent reference is sourced,
# from the package-safe test path, into an environment that cannot see `prod`.
if (!exists("static_obs_mean", mode = "function")) {
  stop("static_obs_mean is not bound: run via the package tests or run_static_observation_checks.R.")
}
prod <- environment(static_obs_mean)
ref <- new.env(parent = globalenv())
sys.source(testthat::test_path("reference", "static_observation_reference.R"), envir = ref)

smean <- prod$static_obs_mean
scov <- prod$static_obs_cov
evidence <- function(...) cat("EVIDENCE:", ..., "\n")
rel_err <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
FD_TOL <- 1e-6

# ---------------------------------------------------------------------------
# Hand fixture. Units 101..106; unit 105 is outside the role (all NA).
#   A  = 1 1 1 1 0 1
#   O  = 1 1 0 1 NA 0
#   dy = 2 4 NA 6 NA NA
#   h  = 1 2 1 1 NA 3
#   q  = .5 .5 NA .25 NA NA   -> phi = 2 2 . 4
#   j  = 1 3 5 7 NA 2
#   Z  = 0 1 2 3 NA 4
hid <- c(101, 102, 103, 104, 105, 106)
hA <- c(1, 1, 1, 1, 0, 1)
hO <- c(1, 1, 0, 1, NA, 0)
hdy <- c(2, 4, NA, 6, NA, NA)
hh <- c(1, 2, 1, 1, NA, 3)
hq <- c(0.5, 0.5, NA, 0.25, NA, NA)
hj <- c(1, 3, 5, 7, NA, 2)
hZ <- c(0, 1, 2, 3, NA, 4)

test_that("reference file is independent of production code", {
  prod_names <- ls(prod)
  for (nm in ls(ref)) {
    obj <- get(nm, ref)
    if (is.function(obj)) {
      expect_length(intersect(all.names(body(obj)), prod_names), 0L)
    }
  }
  expect_false(any(ls(ref) %in% prod_names))
})

test_that("baseline mean, direct IF and covariance match hand arithmetic", {
  r <- smean(hid, hA, hh, "baseline", O = hO, dy = hdy, h_jacobian = list())
  # rows 101,102,104: w = 1,2,1; S = 4; M = (2 + 8 + 6) / 4 = 4
  expect_equal(r$estimate, 4)
  expect_equal(r$denominator_sum, 4)
  # IF = n w (Y - M) / S = 6 * (1*-2, 2*0, 1*2) / 4
  expect_equal(unname(r$direct_if), c(-3, 0, 0, 3, 0, 0))
  expect_identical(names(r$direct_if), as.character(hid))
  expect_equal(r$influence, r$direct_if)          # no channels -> full = direct
  rc <- scov(hid, hA, hZ, "baseline", O = hO, dy = hdy)
  # Z = 0,1,3; Zbar = 4/3; C = [(-4/3)2 + (-1/3)4 + (5/3)6] / 3 = 2
  expect_equal(unname(rc$estimate), 2)
  expect_equal(unname(rc$weighted_design_mean), 4 / 3)
  # direct IF = n w ((Z - Zbar)(Y - M1) - C) / S, S = 3, M1 = 4
  expect_equal(unname(rc$direct_if[, 1]),
               6 * c((-4 / 3) * -2 - 2, (-1 / 3) * 0 - 2, 0, (5 / 3) * 2 - 2, 0, 0) / 3)
})

test_that("weighting mean/covariance and q-loading match hand arithmetic", {
  ones <- matrix(1, 6, 1, dimnames = list(NULL, "a"))
  r <- smean(hid, hA, hh, "weighting", O = hO, dy = hdy, q = hq,
             h_jacobian = list(), q_jacobian = list(eta = ones))
  # w = h phi = 2,4,4; S = 10; M = (4 + 16 + 24) / 10 = 4.4
  expect_equal(r$estimate, 4.4)
  # J_q = sum h dphi (dy - M) / S, dphi = -1/q^2 = -4,-4,-16:
  #   1*-4*-2.4 + 2*-4*-0.4 + 1*-16*1.6 = 9.6 + 3.2 - 25.6 = -12.8; / 10
  expect_equal(unname(r$loadings$eta), -1.28)
  rc <- scov(hid, hA, hZ, "weighting", O = hO, dy = hdy, q = hq, q_jacobian = list())
  # w = phi = 2,2,4; S = 8; Zbar = 14/8 = 1.75
  # C = [2(-1.75)2 + 2(-0.75)4 + 4(1.25)6] / 8 = 17/8
  expect_equal(unname(rc$estimate), 17 / 8)
  expect_equal(unname(rc$weighted_design_mean), 1.75)
})

test_that("imputation and augmented match hand arithmetic, incl. eq:aug-relation", {
  ones <- matrix(1, 6, 1, dimnames = list(NULL, "a"))
  zcol <- matrix(hZ, 6, 1, dimnames = list(NULL, "b"))
  ri <- smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(),
              j_jacobian = list(th = cbind(a = 1, b = hZ)))
  # A rows 101,102,103,104,106: h = 1,2,1,1,3; S = 8; M = (1+6+5+7+6)/8
  expect_equal(ri$estimate, 25 / 8)
  # J_j = sum h b / S: intercept 8/8 = 1; slope (0+2+2+3+12)/8 = 19/8
  expect_equal(unname(ri$loadings$th), c(1, 19 / 8))
  ra <- smean(hid, hA, hh, "augmented", O = hO, dy = hdy, q = hq, j = hj,
              h_jacobian = list(), q_jacobian = list(), j_jacobian = list(th = ones))
  # Y = 3, 5, 5, 3, 2; M = (3 + 10 + 5 + 3 + 6) / 8 = 27/8
  expect_equal(ra$estimate, 27 / 8)
  expect_equal(ra$components$imputation_part, 25 / 8)
  # correction = (1*2*1 + 2*2*1 + 1*4*(-1)) / 8 = 0.25 (divided by sum A h)
  expect_equal(ra$components$correction_part, 0.25)
  # J_j = sum h (1 - O phi) / S = (1(-1) + 2(-1) + 1 + 1(-3) + 3) / 8 = -2/8
  expect_equal(unname(ra$loadings$th), -0.25)
  rci <- scov(hid, hA, hZ, "imputation", j = hj, j_jacobian = list())
  # Zbar = 2; C = [(-2)1 + (-1)3 + 0 + 7 + 2*2] / 5 = 6/5
  expect_equal(unname(rci$estimate), 6 / 5)
  rca <- scov(hid, hA, hZ, "augmented", O = hO, dy = hdy, q = hq, j = hj,
              q_jacobian = list(), j_jacobian = list())
  # C = [(-2)3 + (-1)5 + 0 + 3 + 2*2] / 5 = -4/5
  expect_equal(unname(rca$estimate), -0.8)
})

test_that("augmented and weighting use genuinely different denominators", {
  w <- smean(hid, hA, hh, "weighting", O = hO, dy = hdy, q = hq,
             h_jacobian = list(), q_jacobian = list())
  a0 <- smean(hid, hA, hh, "augmented", O = hO, dy = hdy, q = hq, j = rep(0, 6),
              h_jacobian = list(), q_jacobian = list(), j_jacobian = list())
  # j = 0: Horvitz-Thompson sum h O phi dy / sum A h = (4 + 16 + 24) / 8 = 5.5
  expect_equal(a0$estimate, 5.5)
  expect_equal(w$estimate, 4.4)            # Hajek: / sum A O h phi = 10
  # normalized-correction alternative is NOT used:
  ra <- smean(hid, hA, hh, "augmented", O = hO, dy = hdy, q = hq, j = hj,
              h_jacobian = list(), q_jacobian = list(), j_jacobian = list())
  normalized_alt <- 25 / 8 + (2 + 4 - 4) / 10
  expect_false(isTRUE(all.equal(ra$estimate, normalized_alt)))
  expect_equal(ra$estimate, 25 / 8 + 0.25)
})

test_that("floor is applied once to the joint q, not to each marginal", {
  qs <- c(5e-4, 0.5, 0.5); qu <- c(1e-3, 1, 1)       # product 5e-7 < 1e-6 < marginals
  qprod <- qs * qu
  expect_true(all(qs > 1e-6) && all(qu > 1e-6) && qprod[1] < 1e-6)
  dy3 <- c(1, 3, 5)
  r <- smean(1:3, c(1, 1, 1), c(1, 1, 1), "weighting", O = c(1, 1, 1), dy = dy3, q = qprod,
             h_jacobian = list(), q_jacobian = list())
  # product floor: phi = 1e6, 2, 2
  expect_equal(r$estimate, (1e6 * 1 + 2 * 3 + 2 * 5) / (1e6 + 4))
  marginal_floor <- (2e6 * 1 + 2 * 3 + 2 * 5) / (2e6 + 4)   # 1/(5e-7) = 2e6
  expect_false(isTRUE(all.equal(r$estimate, marginal_floor, tolerance = 1e-12)))
  expect_equal(r$kink$n_q_floored, 1L)
  # Jacobian row of the floored unit is not read (derivative zero there).
  J <- matrix(c(NA, 1, 2), 3, 1)
  r2 <- smean(1:3, c(1, 1, 1), c(1, 1, 1), "weighting", O = c(1, 1, 1), dy = dy3, q = qprod,
              h_jacobian = list(), q_jacobian = list(e = J))
  M <- r2$estimate
  S <- 1e6 + 4
  # rows 2,3: h dphi (dy - M) J / S with dphi = -1/0.25 = -4
  expect_equal(unname(r2$loadings$e), (-4 * (3 - M) * 1 + -4 * (5 - M) * 2) / S)
})

test_that("exact floor boundary is reported (zero derivative), not refused", {
  q <- c(1e-6, 0.5, 0.5)
  r <- smean(1:3, c(1, 1, 1), c(1, 1, 1), "weighting", O = c(1, 1, 1), dy = c(1, 3, 5), q = q,
             h_jacobian = list(), q_jacobian = list(e = matrix(c(NA, 1, 1), 3, 1)))
  expect_equal(r$kink$n_q_at_floor, 1L)
  expect_equal(r$kink$n_q_floored, 0L)
  expect_true(nzchar(r$kink$qualification))
  expect_true(is.finite(r$loadings$e))
})

test_that("unread NAs outside the role / observed intersection do not poison results", {
  filled <- function(x, v) { x[is.na(x)] <- v; x }
  for (route in c("baseline", "weighting", "imputation", "augmented")) {
    args <- list(id = hid, A = hA, h = hh, route = route, h_jacobian = list(hk = cbind(c(1, 2, 3, 4, NA, 6))))
    if (route != "imputation") { args$O <- hO; args$dy <- hdy }
    if (route %in% c("weighting", "augmented")) {
      args$q <- hq
      args$q_jacobian <- list(qk = cbind(c(1, -1, NA, 2, NA, NA)))
    }
    if (route %in% c("imputation", "augmented")) {
      args$j <- hj
      args$j_jacobian <- list(jk = cbind(c(1, 1, 1, 1, NA, 1)))
    }
    if (route %in% c("baseline", "weighting")) {
      args$h <- replace(hh, 3, NA)            # unit 103 unobserved: h not read
      args$h_jacobian$hk[3, 1] <- NA
    }
    r_na <- do.call(smean, args)
    args2 <- args
    args2$h <- filled(args$h, 99)
    args2$h_jacobian$hk <- filled(args$h_jacobian$hk, 99)
    if (!is.null(args2$O)) args2$O <- filled(args$O, 1)
    if (!is.null(args2$dy)) args2$dy <- filled(args$dy, 1e9)
    if (!is.null(args2$q)) { args2$q <- filled(args$q, 0.123); args2$q_jacobian$qk <- filled(args$q_jacobian$qk, 77) }
    if (!is.null(args2$j)) { args2$j <- filled(args$j, -55); args2$j_jacobian$jk <- filled(args$j_jacobian$jk, 5) }
    args2$A <- hA
    # unit 105 has A = 0, so every filler on it is unread; for baseline/weighting
    # unit 103/106 are unobserved and their dy/q/h are unread.
    r_fill <- do.call(smean, args2)
    expect_true(all(is.finite(r_na$direct_if)), info = route)
    expect_equal(r_na$estimate, r_fill$estimate, info = route)
    expect_equal(r_na$direct_if, r_fill$direct_if, info = route)
    expect_equal(r_na$loadings, r_fill$loadings, info = route)
    expect_length(r_na$direct_if, 6L)          # no original unit dropped
  }
  rc <- scov(hid, hA, hZ, "augmented", O = hO, dy = hdy, q = hq, j = hj,
             q_jacobian = list(), j_jacobian = list())
  expect_true(all(is.finite(rc$direct_if)))
  expect_equal(unname(rc$direct_if[5, 1]), 0)
})

test_that("imputation is available for a role with no observed pair", {
  r <- smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(), j_jacobian = list())
  expect_equal(r$n_observed_rows, 0L)
  expect_equal(r$estimate, 25 / 8)
})

# ---------------------------------------------------------------------------
# Augmented route with A O = 0 throughout the role: Y = j, weight
# A h, denominator sum(A h) (mean) / sum(A) (covariance), zero correction.
# The moment, direct IF, h/j loadings and full IF equal imputation; the q
# loading is exactly zero; dy, q and q-Jacobian rows are never read.
test_that("augmented with no observed pair: review fixture reduces to imputation (hand arithmetic)", {
  id4 <- c(10L, 20L, 30L, 40L); A4 <- c(1L, 1L, 0L, 1L); O4 <- c(0L, 0L, NA, 0L)
  h4 <- c(1, 2, NA, 3); j4 <- c(2, 5, NA, -1); Z4 <- c(0, 1, NA, 4)
  na4 <- rep(NA_real_, 4)
  ra <- smean(id4, A4, h4, "augmented", O = O4, dy = na4, q = na4, j = j4,
              h_jacobian = list(), q_jacobian = list(), j_jacobian = list())
  # S = 1 + 2 + 3 = 6; M = (2 + 10 - 3) / 6 = 1.5
  expect_equal(ra$estimate, 1.5)
  expect_equal(ra$denominator_sum, 6)
  expect_equal(ra$n_observed_rows, 0L)
  expect_equal(ra$components$correction_part, 0)
  # direct IF = n h (j - M) / S = 4 (0.5, 7, ., -7.5) / 6
  expect_equal(unname(ra$direct_if), c(1 / 3, 14 / 3, 0, -5))
  expect_equal(ra$influence, ra$direct_if)
  rca <- scov(id4, A4, Z4, "augmented", O = O4, dy = na4, q = na4, j = j4,
              q_jacobian = list(), j_jacobian = list())
  # Zbar = 5/3, M1 = 2; C = [(-5/3)2 + (-2/3)5 + (7/3)(-1)] / 3 = -3
  expect_equal(unname(rca$estimate), -3)
  expect_equal(unname(rca$weighted_design_mean), 5 / 3)
  # direct IF = n {(Z - Zbar)(j - M1) - C} / S = 4 (3, 1, ., -4) / 3
  expect_equal(unname(rca$direct_if[, 1]), c(4, 4 / 3, 0, -16 / 3))
  ri <- smean(id4, A4, h4, "imputation", j = j4, h_jacobian = list(), j_jacobian = list())
  rci <- scov(id4, A4, Z4, "imputation", j = j4, j_jacobian = list())
  expect_equal(ra$direct_if, ri$direct_if)
  expect_equal(rca$direct_if, rci$direct_if)
})

test_that("augmented with no observed pair: mean/covariance, full IF and h/j contraction equal imputation", {
  O0 <- c(0, 0, 0, 0, NA, 0)
  na6 <- rep(NA_real_, 6)
  idx <- c(1, 2, 3, 4, 6)                               # role rows (unit 105 outside)
  Jh <- cbind(g1 = c(1, -1, 2, 0.5, NA, 1), g2 = c(0, 1, 1, -2, NA, 3))
  Jj <- cbind(t1 = c(1, 1, 1, 1, NA, 1), t2 = c(0, 1, 2, 3, NA, 4))
  Jq <- matrix(NA_real_, 6, 1, dimnames = list(NULL, "e1"))   # no row is read
  rn <- as.character(hid)
  Lg <- matrix(c(1, -2, 0.5, 3, -1, 2, 0, 1, -1, 2, 4, -3), 6, 2, dimnames = list(rn, c("g1", "g2")))
  Lt <- matrix(c(2, 0, -1, 1, 3, -2, 1, 1, 0, -1, 2, 5), 6, 2, dimnames = list(rn, c("t1", "t2")))
  Le <- matrix(c(9, -7, 5, 3, -11, 13), 6, 1, dimnames = list(rn, "e1"))

  ra <- smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = hj,
              h_jacobian = list(gamma = Jh), q_jacobian = list(eta = Jq), j_jacobian = list(theta = Jj),
              coef_if = list(gamma = Lg, eta = Le, theta = Lt))
  ri <- smean(hid, hA, hh, "imputation", j = hj,
              h_jacobian = list(gamma = Jh), j_jacobian = list(theta = Jj),
              coef_if = list(gamma = Lg, theta = Lt))
  # hand values: h = 1,2,1,1,3 on the role, S = 8, M = 25/8 (as for imputation)
  expect_equal(ra$estimate, 25 / 8)
  expect_equal(ra$denominator_sum, 8)
  expect_equal(ra$n_observed_rows, 0L)
  expect_equal(ra$n_weight_rows, 5L)
  expect_equal(ra$components$imputation_part, 25 / 8)
  expect_equal(ra$components$correction_part, 0)
  expect_equal(ra$kink$n_q_read, 0L)
  # zero q loading of the right dimension
  expect_identical(unname(ra$loadings$eta), 0)
  # h loading: sum (j - M) Jh / S ; j loading: sum h (1 - O phi) Jj / S = sum h Jj / S
  Jgam <- colSums((hj[idx] - 25 / 8) * Jh[idx, ]) / 8
  Jth <- colSums(hh[idx] * Jj[idx, ]) / 8
  expect_equal(ra$loadings$gamma, Jgam)
  expect_equal(ra$loadings$theta, Jth)
  dir_m <- numeric(6); dir_m[idx] <- 6 * hh[idx] * (hj[idx] - 25 / 8) / 8
  expect_equal(unname(ra$direct_if), dir_m)
  expect_equal(unname(ra$influence), unname(dir_m + drop(Lg %*% Jgam) + drop(Lt %*% Jth)))
  expect_identical(names(ra$influence), rn)
  # equality with imputation (value, direct IF, h/j loadings, full original-unit IF)
  expect_equal(ra$estimate, ri$estimate)
  expect_equal(ra$direct_if, ri$direct_if)
  expect_equal(ra$loadings$gamma, ri$loadings$gamma)
  expect_equal(ra$loadings$theta, ri$loadings$theta)
  expect_equal(ra$influence, ri$influence)
  # independent reference functional (augmented with O = 0 uses j)
  expect_equal(ra$estimate, ref$ref_moments("augmented", hA, O0, na6, hh, na6, hj)$M)

  rca <- scov(hid, hA, hZ, "augmented", O = O0, dy = na6, q = na6, j = hj,
              q_jacobian = list(eta = Jq), j_jacobian = list(theta = Jj),
              coef_if = list(eta = Le, theta = Lt))
  rci <- scov(hid, hA, hZ, "imputation", j = hj, j_jacobian = list(theta = Jj),
              coef_if = list(theta = Lt))
  # Zbar = 2 over the five role rows, C = 6/5 (imputation hand fixture above)
  expect_equal(unname(rca$estimate), 6 / 5)
  expect_equal(unname(rca$weighted_design_mean), 2)
  expect_equal(rca$denominator_sum, 5)
  expect_identical(dim(rca$loadings$eta), c(1L, 1L))
  expect_identical(max(abs(rca$loadings$eta)), 0)
  JthC <- colSums((hZ[idx] - 2) * Jj[idx, ]) / 5
  expect_equal(drop(rca$loadings$theta), JthC)
  expect_equal(rca$estimate, rci$estimate)
  expect_equal(rca$weighted_design_mean, rci$weighted_design_mean)
  expect_equal(rca$direct_if, rci$direct_if)
  expect_equal(rca$loadings$theta, rci$loadings$theta)
  expect_equal(rca$influence, rci$influence)
  expect_equal(unname(rca$estimate), unname(ref$ref_moments("augmented", hA, O0, na6, NULL, na6, hj, hZ)$C))
})

test_that("augmented with no observed pair: loadings are directional derivatives of the reference", {
  O0 <- c(0, 0, 0, 0, NA, 0)
  na6 <- rep(NA_real_, 6)
  Jh <- cbind(c(1, -1, 2, 0.5, NA, 1), c(0, 1, 1, -2, NA, 3))
  Jj <- cbind(c(1, 1, 1, 1, NA, 1), c(0, 1, 2, 3, NA, 4))
  hpar <- function(th) hh + drop(Jh %*% th)              # linear in theta
  jpar <- function(th) hj + drop(Jj %*% th)
  delta <- 1e-6
  r <- smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = hj,
             h_jacobian = list(theta = Jh), q_jacobian = list(), j_jacobian = list(theta = Jj))
  rc <- scov(hid, hA, hZ, "augmented", O = O0, dy = na6, q = na6, j = hj,
             q_jacobian = list(), j_jacobian = list(theta = Jj))
  fdm <- fdc <- numeric(2)
  for (l in 1:2) {
    e <- c(0, 0); e[l] <- delta
    fdm[l] <- (ref$ref_moments("augmented", hA, O0, na6, hpar(e), na6, jpar(e))$M -
               ref$ref_moments("augmented", hA, O0, na6, hpar(-e), na6, jpar(-e))$M) / (2 * delta)
    fdc[l] <- (ref$ref_moments("augmented", hA, O0, na6, NULL, na6, jpar(e), hZ)$C -
               ref$ref_moments("augmented", hA, O0, na6, NULL, na6, jpar(-e), hZ)$C) / (2 * delta)
  }
  expect_lt(rel_err(unname(r$loadings$theta), fdm), 1e-6)
  expect_lt(rel_err(unname(drop(rc$loadings$theta)), fdc), 1e-6)
  evidence(sprintf("no-observed augmented directional mean rel.err=%.2e cov rel.err=%.2e",
                   rel_err(unname(r$loadings$theta), fdm), rel_err(unname(drop(rc$loadings$theta)), fdc)))
})

test_that("augmented with no observed pair reads dy/q/q-Jacobian on no row; required data and contracts kept", {
  O0 <- c(0, 0, 0, 0, NA, 0)
  na6 <- rep(NA_real_, 6)
  Jq <- matrix(NA_real_, 6, 1)
  base <- smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = hj,
                h_jacobian = list(), q_jacobian = list(eta = Jq), j_jacobian = list())
  # values that would be refused if read (q outside [0, 1], nonfinite dy/Jacobian)
  bad <- smean(hid, hA, hh, "augmented", O = replace(O0, 5, 7), dy = rep(Inf, 6), q = rep(7, 6),
               j = replace(hj, 5, NaN), h_jacobian = list(),
               q_jacobian = list(eta = matrix(Inf, 6, 1)), j_jacobian = list())
  expect_equal(bad$estimate, base$estimate)
  expect_equal(bad$direct_if, base$direct_if)
  expect_equal(bad$loadings, base$loadings)
  bc <- scov(hid, hA, replace(hZ, 5, NaN), "augmented", O = O0, dy = rep(-Inf, 6), q = rep(-3, 6),
             j = hj, q_jacobian = list(eta = matrix(NaN, 6, 1)), j_jacobian = list())
  gc <- scov(hid, hA, hZ, "augmented", O = O0, dy = na6, q = na6, j = hj,
             q_jacobian = list(eta = Jq), j_jacobian = list())
  expect_equal(bc$estimate, gc$estimate)
  expect_equal(bc$direct_if, gc$direct_if)
  expect_equal(bc$loadings, gc$loadings)
  # rows that ARE read stay required: O and j on every role row, h on the weight rows
  ie <- "cdid_static_obs_input_error"
  expect_error(smean(hid, hA, hh, "augmented", O = replace(O0, 3, NA), dy = na6, q = na6, j = hj,
                     h_jacobian = list(), q_jacobian = list(), j_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = replace(hj, 6, NA),
                     h_jacobian = list(), q_jacobian = list(), j_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, replace(hh, 3, NA), "augmented", O = O0, dy = na6, q = na6, j = hj,
                     h_jacobian = list(), q_jacobian = list(), j_jacobian = list()), class = ie)
  expect_error(scov(hid, hA, replace(hZ, 2, NA), "augmented", O = O0, dy = na6, q = na6, j = hj,
                    q_jacobian = list(), j_jacobian = list()), class = ie)
  # full-role denominator and empty role are still refused
  se <- "cdid_static_obs_support_error"
  expect_error(smean(hid, hA, c(0, 0, 0, 0, NA, 0), "augmented", O = O0, dy = na6, q = na6, j = hj,
                     h_jacobian = list(), q_jacobian = list(), j_jacobian = list()), class = se)
  expect_error(smean(hid, rep(0, 6), hh, "augmented", O = na6, dy = na6, q = na6, j = hj,
                     h_jacobian = list(), q_jacobian = list(), j_jacobian = list()), class = se)
  # nuisance contracts unchanged: q/j/h Jacobian lists required (list() = known)
  ce <- "cdid_static_obs_channel_error"
  expect_error(smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = hj,
                     h_jacobian = list(), j_jacobian = list()), class = ce)
  expect_error(smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = hj,
                     h_jacobian = list(), q_jacobian = list()), class = ce)
  expect_error(scov(hid, hA, hZ, "augmented", O = O0, dy = na6, q = na6, j = hj,
                    j_jacobian = list()), class = ce)
  # the zero q loading still needs its coefficient IF key (exact key set)
  rn <- as.character(hid)
  Le <- matrix(1, 6, 1, dimnames = list(rn, NULL))
  expect_error(smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = hj,
                     h_jacobian = list(), q_jacobian = list(eta = Jq), j_jacobian = list(),
                     coef_if = list()), class = ce)
  full <- smean(hid, hA, hh, "augmented", O = O0, dy = na6, q = na6, j = hj,
                h_jacobian = list(), q_jacobian = list(eta = Jq), j_jacobian = list(),
                coef_if = list(eta = Le))
  expect_equal(full$influence, full$direct_if)
  # dy stays a mandatory argument of the route (a vector, even if no row is read)
  expect_error(smean(hid, hA, hh, "augmented", O = O0, q = na6, j = hj,
                     h_jacobian = list(), q_jacobian = list(), j_jacobian = list()), class = ie)
})

test_that("insufficient support, zero denominators and nonfinite required data are refused", {
  cls <- "cdid_static_obs_support_error"
  expect_error(smean(hid, hA, hh, "baseline", O = c(0, 0, 0, 0, NA, 0), dy = hdy, h_jacobian = list()), class = cls)
  # observed support is required by baseline/weighting only; the
  # augmented route with no observed pair reduces to imputation (tested below)
  expect_error(smean(hid, hA, hh, "weighting", O = c(0, 0, 0, 0, NA, 0), dy = hdy, q = hq,
                     h_jacobian = list(), q_jacobian = list()), class = cls)
  expect_error(scov(hid, hA, hZ, "baseline", O = c(0, 0, 0, 0, NA, 0), dy = hdy), class = cls)
  expect_error(scov(hid, hA, hZ, "weighting", O = c(0, 0, 0, 0, NA, 0), dy = hdy, q = hq,
                    q_jacobian = list()), class = cls)
  expect_error(smean(hid, hA, c(0, 0, 5, 0, NA, 5), "weighting", O = hO, dy = hdy, q = hq,
                     h_jacobian = list(), q_jacobian = list()), class = cls)
  expect_error(smean(hid, hA, c(0, 0, 0, 0, NA, 0), "imputation", j = hj,
                     h_jacobian = list(), j_jacobian = list()), class = cls)
  ie <- "cdid_static_obs_input_error"
  expect_error(smean(hid, hA, hh, "baseline", O = hO, dy = replace(hdy, 2, NA), h_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, hh, "baseline", O = replace(hO, 3, NA), dy = hdy, h_jacobian = list()), class = ie)
  expect_error(smean(hid, replace(hA, 5, NA), hh, "baseline", O = hO, dy = hdy, h_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, replace(hh, 1, -1), "baseline", O = hO, dy = hdy, h_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, replace(hh, 2, Inf), "baseline", O = hO, dy = hdy, h_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, hh, "weighting", O = hO, dy = hdy, q = replace(hq, 1, 1.2),
                     h_jacobian = list(), q_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, hh, "weighting", O = hO, dy = hdy, q = replace(hq, 4, NA),
                     h_jacobian = list(), q_jacobian = list()), class = ie)
  expect_error(smean(hid, hA, hh, "imputation", j = replace(hj, 6, NaN),
                     h_jacobian = list(), j_jacobian = list()), class = ie)
  expect_error(scov(hid, hA, replace(hZ, 3, NA), "imputation", j = hj, j_jacobian = list()), class = ie)
  expect_error(smean(c(1, 1, 2, 3, 4, 5), hA, hh, "baseline", O = hO, dy = hdy, h_jacobian = list()), class = ie)
})

test_that("routes are explicit: invalid route and unused inputs are refused", {
  re <- "cdid_static_obs_route_error"
  expect_error(smean(hid, hA, hh, "ipw", O = hO, dy = hdy, h_jacobian = list()), class = re)
  expect_error(smean(hid, hA, hh, "baseline", O = hO, dy = hdy, q = hq, h_jacobian = list()), class = re)
  expect_error(smean(hid, hA, hh, "imputation", O = hO, j = hj, h_jacobian = list(), j_jacobian = list()), class = re)
  expect_error(smean(hid, hA, hh, "weighting", O = hO, dy = hdy, q = hq, j = hj,
                     h_jacobian = list(), q_jacobian = list()), class = re)
  ce <- "cdid_static_obs_channel_error"
  expect_error(smean(hid, hA, hh, "baseline", O = hO, dy = hdy), class = ce)          # h_jacobian NULL
  expect_error(smean(hid, hA, hh, "weighting", O = hO, dy = hdy, q = hq, h_jacobian = list()), class = ce)
})

test_that("channel keys: repeats refused, cross-nuisance sharing summed, IF keys/order checked", {
  ce <- "cdid_static_obs_channel_error"
  oe <- "cdid_static_obs_order_error"
  J1 <- cbind(a = c(1, 2, 3, 4, NA, 6))
  dupl <- list(J1, J1); names(dupl) <- c("k", "k")
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = dupl, j_jacobian = list()), class = ce)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(J1), j_jacobian = list()), class = ce)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1[1:5, , drop = FALSE]),
                     j_jacobian = list()), class = ce)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = replace(J1, 2, NA)),
                     j_jacobian = list()), class = ce)
  # same key in two nuisances = same parameter: loadings are summed
  J2 <- cbind(a = c(0.5, -1, 2, 0, NA, 1))
  sep <- smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(kh = J1), j_jacobian = list(kj = J2))
  shr <- smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list(k = J2))
  expect_equal(shr$loadings$k, sep$loadings$kh + sep$loadings$kj)
  # equal arrays under different keys are distinct parameters: never deduplicated
  twin <- smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k1 = J1), j_jacobian = list(k2 = J1))
  expect_setequal(names(twin$loadings), c("k1", "k2"))
  L1 <- matrix(1, 6, 1, dimnames = list(as.character(hid), "a"))
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k1 = J1), j_jacobian = list(k2 = J1),
                     coef_if = list(k1 = L1)), class = ce)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1),
                     j_jacobian = list(k = cbind(J2, J2))), class = ce)
  # coefficient IFs: keys must match exactly; rownames = id in order
  L <- matrix(c(1, -2, 3, 0.5, 4, -1), 6, 1, dimnames = list(as.character(hid), "a"))
  ok <- smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list(k = J2),
              coef_if = list(k = L))
  expect_equal(unname(ok$influence), unname(ok$direct_if + drop(L %*% shr$loadings$k)))
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(kh = J1), j_jacobian = list(kj = J2),
                     coef_if = list(kh = L)), class = ce)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list(),
                     coef_if = list(k = L, extra = L)), class = ce)
  Lrev <- L[6:1, , drop = FALSE]
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list(),
                     coef_if = list(k = Lrev)), class = oe)
  Lnon <- unname(L)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list(),
                     coef_if = list(k = Lnon)), class = oe)
  Lw <- cbind(L, b = 1)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list(),
                     coef_if = list(k = Lw)), class = ce)
  dupL <- list(L, L); names(dupL) <- c("k", "k")
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list(),
                     coef_if = dupL), class = ce)
  # labelled vectors / Jacobian rows in another order are refused
  expect_error(smean(hid, hA, stats::setNames(hh, rev(hid)), "imputation", j = hj,
                     h_jacobian = list(), j_jacobian = list()), class = oe)
  Jr <- J1; rownames(Jr) <- rev(hid)
  expect_error(smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = Jr), j_jacobian = list()), class = oe)
  # loadings without coefficient IFs: no full IF claimed
  noif <- smean(hid, hA, hh, "imputation", j = hj, h_jacobian = list(k = J1), j_jacobian = list())
  expect_null(noif$influence)
})

# ---------------------------------------------------------------------------
# Directional derivatives of the loadings (conditional algebra), with clipping
# active on some rows and inactive on others, away from the kink. One shared
# parameter theta enters q, j and h; the reference is re-evaluated at theta +/-.
test_that("loadings equal directional derivatives of the independent reference", {
  set.seed(5)
  n <- 14
  id <- sample(1e5, n)
  A <- c(0, rep(1, n - 1))
  O <- c(NA, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1)
  x <- round(stats::rnorm(n), 2)
  B <- cbind(1, x)
  off <- c(0, -22, -22, rep(0, n - 3))          # rows 2,3: q ~ 1e-9 (floored)
  dyv <- ifelse(O %in% 1, round(stats::rnorm(n, 1), 2), NA)
  Z <- cbind(x, x^2)
  th0 <- c(0.2, -0.4)
  pr <- function(th) {
    list(q = 1 / (1 + exp(-(off + drop(B %*% th)))), j = drop(B %*% th) * 1.5 + 0.3,
         h = exp(0.3 * drop(B %*% th)))
  }
  p0 <- pr(th0)
  expect_true(all(p0$q[2:3] < 1e-7) && all(p0$q[-(1:3)] > 1e-3))
  Jq <- p0$q * (1 - p0$q) * B
  Jj <- 1.5 * B
  Jh <- 0.3 * p0$h * B
  delta <- 1e-6
  for (route in c("baseline", "weighting", "imputation", "augmented")) {
    for (kind in c("mean", "cov")) {
      args <- list(id = id, A = A, route = route)
      if (kind == "mean") { args$h <- p0$h; args$h_jacobian <- list(theta = Jh) } else args$Z <- Z
      if (route != "imputation") { args$O <- O; args$dy <- dyv }
      if (route %in% c("weighting", "augmented")) { args$q <- p0$q; args$q_jacobian <- list(theta = Jq) }
      if (route %in% c("imputation", "augmented")) { args$j <- p0$j; args$j_jacobian <- list(theta = Jj) }
      r <- do.call(if (kind == "mean") smean else scov, args)
      fd <- sapply(1:2, function(l) {
        e <- c(0, 0); e[l] <- delta
        pp <- pr(th0 + e); pm <- pr(th0 - e)
        hp <- if (kind == "mean") pp$h else NULL
        hm <- if (kind == "mean") pm$h else NULL
        rp <- ref$ref_moments(route, A, O, dyv, hp, pp$q, pp$j, Z)
        rm <- ref$ref_moments(route, A, O, dyv, hm, pm$q, pm$j, Z)
        if (kind == "mean") (rp$M - rm$M) / (2 * delta) else (rp$C - rm$C) / (2 * delta)
      })
      an <- r$loadings$theta
      if (kind == "mean") fd <- drop(fd)
      if (route == "baseline" && kind == "cov") {
        # no nuisance enters the baseline covariance: no channel, zero derivative
        expect_length(r$loadings, 0L)
        expect_identical(max(abs(fd)), 0)
        an <- 0 * fd
      } else {
        expect_identical(dim(as.matrix(an)), dim(as.matrix(fd)))
      }
      expect_lt(rel_err(unname(an), unname(fd)), 1e-6)
      v0 <- ref$ref_moments(route, A, O, dyv, if (kind == "mean") p0$h else NULL, p0$q, p0$j, Z)
      if (kind == "mean") expect_equal(r$estimate, v0$M) else expect_equal(unname(r$estimate), unname(v0$C))
      if (kind == "cov") expect_equal(unname(r$weighted_design_mean), unname(v0$Zbar))
      evidence(sprintf("directional %s/%s rel.err=%.2e", route, kind, rel_err(unname(an), unname(fd))))
    }
  }
})

# ---------------------------------------------------------------------------
# Refit finite-difference fixtures: every nuisance is refitted under unit
# empirical weights nu; n * dF/dnu_i is compared with the analytic full IF.
call_kernel <- function(kernel, kind, fx, inp, cell, route, prod_q = NULL) {
  args <- list(id = fx$id, A = fx$A, route = route)
  cif <- list()
  if (kind == "mean") {
    args$h <- inp$pred$h
    args$h_jacobian <- list(gamma = inp$J_h)
    cif$gamma <- inp$L_gam
  } else {
    args$Z <- fx$Z
  }
  if (route != "imputation") {
    O <- fx$Ok[[cell]]; O[fx$A == 0L] <- NA          # unread outside the role
    dy <- fx$dyk[[cell]]; dy[fx$A == 0L] <- NA
    args$O <- O; args$dy <- dy
  }
  if (route %in% c("weighting", "augmented")) {
    if (is.null(prod_q)) {
      q <- inp$pred$q[[cell]]
      q[fx$Ok[[cell]] == 0L | fx$A == 0L] <- NA       # q read only on observed role rows
      args$q <- q
      args$q_jacobian <- list(eta_pool = inp$J_q[[cell]])
      cif$eta_pool <- inp$L_eta
    } else if (prod_q == "separate") {
      args$q <- inp$pred$qprod
      args$q_jacobian <- list(eta_date_s = inp$J_prod_s, eta_date_u = inp$J_prod_u)
      cif$eta_date_s <- inp$L_s
      cif$eta_date_u <- inp$L_u
    } else {
      args$q <- inp$pred$qprod_sh
      args$q_jacobian <- list(eta_date = inp$J_prod_sh)
      cif$eta_date <- inp$L_sh
    }
  }
  if (route %in% c("imputation", "augmented")) {
    args$j <- inp$pred$j[[cell]]
    args$j_jacobian <- list(theta_pool = inp$J_j[[cell]])
    cif$theta_pool <- inp$L_th
  }
  args$coef_if <- cif
  do.call(kernel, args)
}

fx0 <- ref$ref_make_fixture(seed = 11L)
inp0 <- ref$ref_inputs(fx0)
fd0 <- ref$ref_fd_if(fx0)
fxF <- ref$ref_make_fixture(seed = 11L, pair_offset = -16, date_offset = -9)
inpF <- ref$ref_inputs(fxF)
fdF <- ref$ref_fd_if(fxF)

test_that("refit FD: full IF of mean (with fitted tilt) and covariance, all routes, both pairs", {
  for (cell in 1:2) for (route in c("baseline", "weighting", "imputation", "augmented")) {
    rm <- call_kernel(smean, "mean", fx0, inp0, cell, route)
    key <- paste0("M_", route, "_k", cell)
    expect_equal(rm$estimate, unname(fd0$value[key]), tolerance = 1e-12)
    e <- rel_err(rm$influence, fd0$fd[, key])
    expect_lt(e, FD_TOL)
    expect_identical(names(rm$influence), as.character(fx0$id))
    rc <- call_kernel(scov, "cov", fx0, inp0, cell, route)
    kc <- paste0(c("C1_", "C2_"), route, "_k", cell)
    expect_equal(unname(rc$estimate), unname(fd0$value[kc]), tolerance = 1e-12)
    ec <- rel_err(unname(rc$influence), unname(fd0$fd[, kc]))
    expect_lt(ec, FD_TOL)
    evidence(sprintf("refit FD k%d %-10s mean rel.err=%.2e  cov rel.err=%.2e", cell, route, e, ec))
  }
})

test_that("zero/single-outcome units and switchers enter through full means and nuisances", {
  zero_single <- as.character(fx0$id[17:20])   # stayers with 0 or 1 recorded outcome
  switchers <- as.character(fx0$id[1:10])
  ri <- call_kernel(smean, "mean", fx0, inp0, 1L, "imputation")
  ra <- call_kernel(smean, "mean", fx0, inp0, 1L, "augmented")
  rw <- call_kernel(smean, "mean", fx0, inp0, 1L, "weighting")
  # direct term nonzero for no-pair role units in imputation/augmented
  expect_true(all(ri$direct_if[zero_single] != 0))
  expect_true(all(ra$direct_if[zero_single] != 0))
  # weighting: direct term zero for them, but the IF is nonzero via q and h
  expect_true(all(rw$direct_if[zero_single] == 0))
  expect_true(all(rw$influence[zero_single] != 0))
  expect_lt(max(abs(rw$influence[zero_single] - fd0$fd[zero_single, "M_weighting_k1"])), 1e-6)
  # switchers (outside the role) enter only via the treatment tilt fit
  expect_true(all(rw$direct_if[switchers] == 0))
  expect_true(all(abs(rw$influence[switchers]) > 0))
  expect_lt(max(abs(rw$influence[switchers] - fd0$fd[switchers, "M_weighting_k1"])), 1e-6)
})

test_that("shared pooled channel: one contraction across pairs equals sum of per-pair IFs", {
  r1 <- call_kernel(smean, "mean", fx0, inp0, 1L, "augmented")
  r2 <- call_kernel(smean, "mean", fx0, inp0, 2L, "augmented")
  keys <- names(r1$loadings)
  expect_setequal(keys, names(r2$loadings))
  diff_load <- lapply(stats::setNames(keys, keys), function(k) r1$loadings[[k]] - r2$loadings[[k]])
  once <- r1$direct_if - r2$direct_if +
    drop(prod$static_obs_contract(diff_load, list(gamma = inp0$L_gam, eta_pool = inp0$L_eta,
                                                  theta_pool = inp0$L_th), fx0$id))
  expect_equal(once, r1$influence - r2$influence)
  expect_lt(rel_err(once, fd0$fd[, "M_augmented_k1"] - fd0$fd[, "M_augmented_k2"]), FD_TOL)
})

test_that("refit FD with the floor active (pair logit and date product) away from kinks", {
  for (cell in 1:2) for (route in c("weighting", "augmented")) {
    rm <- call_kernel(smean, "mean", fxF, inpF, cell, route)
    expect_gt(rm$kink$n_q_floored, 0L)
    expect_equal(rm$kink$n_q_at_floor, 0L)
    key <- paste0("M_", route, "_k", cell)
    e <- rel_err(rm$influence, fdF$fd[, key])
    expect_lt(e, FD_TOL)
    rc <- call_kernel(scov, "cov", fxF, inpF, cell, route)
    kc <- paste0(c("C1_", "C2_"), route, "_k", cell)
    ec <- rel_err(unname(rc$influence), unname(fdF$fd[, kc]))
    expect_lt(ec, FD_TOL)
    evidence(sprintf("floor-active FD k%d %-9s floored=%d mean rel.err=%.2e cov rel.err=%.2e",
                     cell, route, rm$kink$n_q_floored, e, ec))
  }
  obs <- fxF$A == 1L & fxF$Ok[[1]] == 1L
  qs <- inpF$pred$qs[obs]; qu <- inpF$pred$qu[obs]
  expect_true(any(qs * qu < 1e-6 & qs > 1e-6 & qu > 1e-6))   # product floor only
  for (route in c("weighting", "augmented")) {
    rs <- call_kernel(smean, "mean", fxF, inpF, 1L, route, prod_q = "separate")
    e1 <- rel_err(rs$influence, fdF$fd[, paste0("Mprod_", route)])
    expect_lt(e1, FD_TOL)
    expect_equal(rs$estimate, unname(fdF$value[paste0("Mprod_", route)]), tolerance = 1e-12)
    rcs <- call_kernel(scov, "cov", fxF, inpF, 1L, route, prod_q = "separate")
    e2 <- rel_err(unname(rcs$influence[, 1]), fdF$fd[, paste0("C1prod_", route)])
    expect_lt(e2, FD_TOL)
    rsh <- call_kernel(smean, "mean", fxF, inpF, 1L, route, prod_q = "shared")
    e3 <- rel_err(rsh$influence, fdF$fd[, paste0("Mprodsh_", route)])
    expect_lt(e3, FD_TOL)
    rcsh <- call_kernel(scov, "cov", fxF, inpF, 1L, route, prod_q = "shared")
    e4 <- rel_err(unname(rcsh$influence[, 2]), fdF$fd[, paste0("C2prodsh_", route)])
    expect_lt(e4, FD_TOL)
    evidence(sprintf("date-product FD %-9s separate mean=%.2e cov=%.2e shared mean=%.2e cov=%.2e",
                     route, e1, e2, e3, e4))
  }
})

# ---------------------------------------------------------------------------
# Deliberate formula perturbations: a mutated copy of the production code must
# fail the same meaningful checks that the production code passes.
mutate <- function(edits) {
  env <- new.env(parent = prod)
  for (fname in names(edits)) {
    txt <- paste(deparse(get(fname, prod), width.cutoff = 500L), collapse = "\n")
    from <- edits[[fname]][1]
    hits <- gregexpr(from, txt, fixed = TRUE)[[1]]
    if (sum(hits > 0) != 1L) stop("mutation target not found exactly once: ", from)
    f <- eval(parse(text = sub(from, edits[[fname]][2], txt, fixed = TRUE)))
    environment(f) <- env
    assign(fname, f, envir = env)
  }
  for (w in c("static_obs_mean", "static_obs_cov", "static_obs_core", "static_obs_floor_inverse")) {
    if (!exists(w, envir = env, inherits = FALSE)) {
      f <- get(w, prod); environment(f) <- env; assign(w, f, envir = env)
    }
  }
  env
}

test_that("perturbation 1: derivative not zeroed below the floor fails the refit FD check", {
  mt <- mutate(list(static_obs_floor_inverse = c("act <- ipw_clip_floor_active(q)",
                                                 "act <- as.double(q >= 0)")))
  good <- call_kernel(smean, "mean", fxF, inpF, 1L, "weighting")
  bad <- call_kernel(mt$static_obs_mean, "mean", fxF, inpF, 1L, "weighting")
  eg <- rel_err(good$influence, fdF$fd[, "M_weighting_k1"])
  eb <- rel_err(bad$influence, fdF$fd[, "M_weighting_k1"])
  evidence(sprintf("perturbation 1 (floor derivative): production %.2e, mutant %.2e", eg, eb))
  expect_lt(eg, FD_TOL)
  expect_gt(eb, 1e3 * FD_TOL)
  expect_equal(bad$estimate, good$estimate)       # the value alone cannot detect it
})

test_that("perturbation 2: augmented j-loading without (1 - O phi) fails the refit FD check", {
  mt <- mutate(list(static_obs_core = c("else 1 - phi_w", "else rep(1, nw)")))
  good <- call_kernel(smean, "mean", fx0, inp0, 1L, "augmented")
  bad <- call_kernel(mt$static_obs_mean, "mean", fx0, inp0, 1L, "augmented")
  eg <- rel_err(good$influence, fd0$fd[, "M_augmented_k1"])
  eb <- rel_err(bad$influence, fd0$fd[, "M_augmented_k1"])
  evidence(sprintf("perturbation 2 (aug j-loading): production %.2e, mutant %.2e", eg, eb))
  expect_lt(eg, FD_TOL)
  expect_gt(eb, 1e3 * FD_TOL)
})

test_that("perturbation 3: augmented correction renormalized by inverse-weight mass fails", {
  mt <- mutate(list(static_obs_core = c("Y <- j_w + phi_w * resid_w",
                                        "Y <- j_w + phi_w * resid_w * sum(h_w)/sum(h_w * phi_w)")))
  good <- call_kernel(smean, "mean", fx0, inp0, 1L, "augmented")
  bad <- call_kernel(mt$static_obs_mean, "mean", fx0, inp0, 1L, "augmented")
  ref_value <- unname(fd0$value["M_augmented_k1"])
  evidence(sprintf("perturbation 3 (aug renormalized): production |diff|=%.2e, mutant |diff|=%.2e",
                   abs(good$estimate - ref_value), abs(bad$estimate - ref_value)))
  expect_equal(good$estimate, ref_value, tolerance = 1e-12)
  expect_gt(abs(bad$estimate - ref_value), 1e-3)
})

test_that("complete-data baseline with unit tilt reduces to the role mean of dy", {
  fx <- fx0
  dy_full <- (fx$Z[, 1] + 1:fx$n / 7)                # arbitrary complete changes
  r <- smean(fx$id, fx$A, rep(1, fx$n), "baseline", O = rep(1, fx$n), dy = dy_full, h_jacobian = list())
  expect_equal(r$estimate, mean(dy_full[fx$A == 1L]))
  expect_equal(r$denominator_mean, mean(fx$A))
})
