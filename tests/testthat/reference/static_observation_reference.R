# Independent mathematical reference and refit fixtures for the static
# observation kernel tests. NOTHING here calls the production kernel or its
# helpers (no static_obs_*, ipw_clip_*, cdid_*). The floor 1e-6 is written as
# a literal here on purpose, so a change of the production constant is caught.
#
# Contents
#   ref_moments()      unit-loop evaluation of M_hat(h) and C_hat with optional
#                      unit empirical weights nu (accumulators, uncentered
#                      covariance algebra E[ZY] - E[Z]E[Y]; different from the
#                      production centered form)
#   ref_logit_newton(), ref_wls()   weighted refits (own Newton / normal eqs)
#   ref_make_fixture() small panel with zero/single-outcome units, stayer role,
#                      pooled pair logit (shared slope, local intercepts),
#                      pooled imputation LS, treatment logit, date logits
#   ref_pipeline()     refit every nuisance under unit weights nu, return the
#                      functionals of interest
#   ref_fd_if()        n-times central-difference derivative in nu_i
#   ref_inputs()       predictions, Jacobians and coefficient IFs at nu = 1
#                      from closed-form M-estimator formulas

ref_eps <- 1e-6

ref_inv_floor <- function(q) if (q > ref_eps) 1 / q else 1 / ref_eps

ref_moments <- function(route, A, O, dy, h, q, j, Z = NULL, nu = NULL) {
  n <- length(A)
  if (is.null(nu)) nu <- rep(1, n)
  if (!is.null(Z)) Z <- as.matrix(Z)
  d <- if (is.null(Z)) 0L else ncol(Z)
  sw <- 0; swy <- 0                 # tilted (mean) accumulators
  s1 <- 0; s1y <- 0                 # untilted (covariance) accumulators
  s1z <- numeric(d); s1zy <- numeric(d)
  for (i in seq_len(n)) {
    if (A[i] != 1) next
    if (route %in% c("baseline", "weighting") && O[i] != 1) next
    om <- if (route == "weighting") ref_inv_floor(q[i]) else 1
    y <- if (route %in% c("baseline", "weighting")) {
      dy[i]
    } else if (route == "imputation") {
      j[i]
    } else if (O[i] == 1) {
      j[i] + ref_inv_floor(q[i]) * (dy[i] - j[i])
    } else {
      j[i]
    }
    if (!is.null(h)) {
      sw <- sw + nu[i] * om * h[i]
      swy <- swy + nu[i] * om * h[i] * y
    }
    s1 <- s1 + nu[i] * om
    s1y <- s1y + nu[i] * om * y
    if (d) {
      s1z <- s1z + nu[i] * om * Z[i, ]
      s1zy <- s1zy + nu[i] * om * Z[i, ] * y
    }
  }
  list(
    M = if (!is.null(h)) swy / sw else NA_real_,
    M1 = s1y / s1,
    Zbar = if (d) s1z / s1 else NULL,
    C = if (d) s1zy / s1 - (s1z / s1) * (s1y / s1) else NULL
  )
}

ref_logit_newton <- function(y, X, w, offset = 0) {
  b <- rep(0, ncol(X))
  for (it in 1:200) {
    eta <- offset + drop(X %*% b)
    p <- 1 / (1 + exp(-eta))
    g <- crossprod(X, w * (y - p))
    H <- crossprod(X, X * (w * p * (1 - p)))
    step <- solve(H, g)
    b <- b + drop(step)
    if (max(abs(step)) < 1e-14) break
  }
  if (max(abs(step)) > 1e-10) stop("reference logit did not converge")
  b
}

ref_wls <- function(y, X, w) {
  drop(solve(crossprod(X, X * w), crossprod(X, w * y)))
}

# Panel fixture. Role = stayers (A = 1); switchers have A = 0 and enter only
# through the treatment logit. Dates 1..3; cells k1 = (2,3), k2 = (1,3).
ref_make_fixture <- function(seed = 11L, n = 40L, pair_offset = 0, date_offset = 0) {
  set.seed(seed)
  id <- 1000L + sample.int(500L, n)              # non-sequential original ids
  x1 <- round(stats::rnorm(n), 3)
  x2 <- round(stats::runif(n, -1, 1), 3)
  flag <- integer(n)
  A <- as.integer(seq_len(n) > 10L)              # units 1..10 switchers
  G <- 1L - A
  flag[c(11L, 12L, 13L, 14L, 15L, 16L)] <- 1L    # remote stayers (offset units)
  Y <- matrix(NA_real_, n, 3L)
  base <- 1 + 0.5 * x1 - 0.3 * x2
  for (s in 1:3) Y[, s] <- base + 0.4 * s + 0.2 * s * x2 + round(stats::rnorm(n, sd = 0.5), 3)
  Y <- round(Y, 3)
  Ob <- matrix(stats::runif(n * 3L) < stats::plogis(0.8 + 0.6 * x1), n, 3L)
  # forced patterns: zero-outcome and single-outcome stayers, observed remote units
  Ob[17L, ] <- FALSE; Ob[18L, ] <- FALSE              # zero outcomes
  Ob[19L, ] <- c(FALSE, FALSE, TRUE)                  # single outcome
  Ob[20L, ] <- c(TRUE, FALSE, FALSE)                  # single outcome
  Ob[11L, ] <- TRUE; Ob[12L, ] <- TRUE                # remote but fully observed
  Ob[13L, ] <- FALSE; Ob[14L, ] <- c(TRUE, FALSE, FALSE)
  Ob[15L, ] <- c(FALSE, TRUE, FALSE); Ob[16L, ] <- c(FALSE, FALSE, TRUE)
  Ob[21L, ] <- c(TRUE, TRUE, FALSE); Ob[22L, ] <- c(FALSE, TRUE, TRUE)
  Ob[23L, ] <- c(TRUE, FALSE, TRUE); Ob[24L, ] <- c(TRUE, TRUE, TRUE)
  Ob[25L, ] <- c(FALSE, TRUE, TRUE)
  Yobs <- Y
  Yobs[!Ob] <- NA_real_
  cells <- list(k1 = c(2L, 3L), k2 = c(1L, 3L))
  Ok <- lapply(cells, function(k) as.integer(Ob[, k[1]] & Ob[, k[2]]))
  dyk <- lapply(cells, function(k) Yobs[, k[2]] - Yobs[, k[1]])  # NA unless both observed
  list(
    n = n, id = id, x1 = x1, x2 = x2, flag = flag, A = A, G = G, Ob = Ob,
    Yobs = Yobs, cells = cells, Ok = Ok, dyk = dyk,
    Z = cbind(z1 = x1, z2 = x2 + 0.5 * x1^2),
    Xt = cbind(1, x1, x2),                               # treatment logit design
    pair_offset = pair_offset, date_offset = date_offset
  )
}

# Unit-level design of the pooled pair logit for cell index c (1 or 2):
# local intercepts (2) + shared slope on x1.
ref_pair_design <- function(fx, c) {
  cbind(a1 = rep(as.numeric(c == 1L), fx$n), a2 = rep(as.numeric(c == 2L), fx$n), b = fx$x1)
}
ref_imp_design <- function(fx, c) {
  cbind(a1 = rep(as.numeric(c == 1L), fx$n), a2 = rep(as.numeric(c == 2L), fx$n),
        d1 = fx$x1, d2 = fx$x2)
}
ref_date_design <- function(fx) cbind(a = rep(1, fx$n), b = fx$x1)

ref_fits <- function(fx, nu) {
  n <- fx$n
  role <- which(fx$A == 1L)
  gam <- ref_logit_newton(fx$G, fx$Xt, nu)
  # pooled pair logit over cells k1, k2 of the role (stacked)
  Bs <- rbind(ref_pair_design(fx, 1L)[role, ], ref_pair_design(fx, 2L)[role, ])
  ys <- c(fx$Ok[[1]][role], fx$Ok[[2]][role])
  ws <- c(nu[role], nu[role])
  offs <- rep(fx$pair_offset * fx$flag[role], 2L)
  eta <- ref_logit_newton(ys, Bs, ws, offs)
  # pooled imputation LS over role x observed rows of k1, k2
  r1 <- role[fx$Ok[[1]][role] == 1L]
  r2 <- role[fx$Ok[[2]][role] == 1L]
  Bj <- rbind(ref_imp_design(fx, 1L)[r1, ], ref_imp_design(fx, 2L)[r2, ])
  yj <- c(fx$dyk[[1]][r1], fx$dyk[[2]][r2])
  th <- ref_wls(yj, Bj, c(nu[r1], nu[r2]))
  # separate date logits for dates 2 and 3 on the role
  Dd <- ref_date_design(fx)[role, ]
  doff <- fx$date_offset * fx$flag[role]
  eta_s <- ref_logit_newton(as.numeric(fx$Ob[role, 2L]), Dd, nu[role], doff)
  eta_u <- ref_logit_newton(as.numeric(fx$Ob[role, 3L]), Dd, nu[role], doff)
  # shared date block: local intercepts per date, shared slope on x1
  Dsh <- rbind(cbind(1, 0, fx$x1[role]), cbind(0, 1, fx$x1[role]))
  eta_sh <- ref_logit_newton(c(as.numeric(fx$Ob[role, 2L]), as.numeric(fx$Ob[role, 3L])),
                             Dsh, c(nu[role], nu[role]), c(doff, doff))
  list(gam = gam, eta = eta, th = th, eta_s = eta_s, eta_u = eta_u, eta_sh = eta_sh)
}

ref_predictions <- function(fx, f) {
  lg <- function(u) 1 / (1 + exp(-u))
  h <- exp(drop(fx$Xt %*% f$gam))
  q <- lapply(1:2, function(c) lg(fx$pair_offset * fx$flag + drop(ref_pair_design(fx, c) %*% f$eta)))
  j <- lapply(1:2, function(c) drop(ref_imp_design(fx, c) %*% f$th))
  D <- ref_date_design(fx)
  doff <- fx$date_offset * fx$flag
  qs <- lg(doff + drop(D %*% f$eta_s))
  qu <- lg(doff + drop(D %*% f$eta_u))
  qs_sh <- lg(doff + f$eta_sh[1] + f$eta_sh[3] * fx$x1)
  qu_sh <- lg(doff + f$eta_sh[2] + f$eta_sh[3] * fx$x1)
  list(h = h, q = q, j = j, qs = qs, qu = qu, qprod = qs * qu,
       qs_sh = qs_sh, qu_sh = qu_sh, qprod_sh = qs_sh * qu_sh)
}

# All functionals of interest under unit weights nu (everything refitted).
ref_pipeline <- function(fx, nu) {
  f <- ref_fits(fx, nu)
  p <- ref_predictions(fx, f)
  out <- c()
  for (c in 1:2) {
    O <- fx$Ok[[c]]; dy <- fx$dyk[[c]]
    for (route in c("baseline", "weighting", "imputation", "augmented")) {
      r <- ref_moments(route, fx$A, O, dy, p$h, p$q[[c]], p$j[[c]], fx$Z, nu)
      out[paste0("M_", route, "_k", c)] <- r$M
      out[paste0("C1_", route, "_k", c)] <- r$C[1]
      out[paste0("C2_", route, "_k", c)] <- r$C[2]
    }
  }
  # date-product routes on k1 = (2,3)
  O <- fx$Ok[[1]]; dy <- fx$dyk[[1]]
  for (route in c("weighting", "augmented")) {
    r <- ref_moments(route, fx$A, O, dy, p$h, p$qprod, p$j[[1]], fx$Z, nu)
    out[paste0("Mprod_", route)] <- r$M
    out[paste0("C1prod_", route)] <- r$C[1]
    r <- ref_moments(route, fx$A, O, dy, p$h, p$qprod_sh, p$j[[1]], fx$Z, nu)
    out[paste0("Mprodsh_", route)] <- r$M
    out[paste0("C2prodsh_", route)] <- r$C[2]
  }
  out
}

ref_fd_if <- function(fx, step = 1e-4) {
  n <- fx$n
  base <- ref_pipeline(fx, rep(1, n))
  out <- matrix(NA_real_, n, length(base), dimnames = list(as.character(fx$id), names(base)))
  for (i in seq_len(n)) {
    up <- dn <- rep(1, n)
    up[i] <- 1 + step
    dn[i] <- 1 - step
    out[i, ] <- n * (ref_pipeline(fx, up) - ref_pipeline(fx, dn)) / (2 * step)
  }
  list(value = base, fd = out)
}

# Predictions, prediction-row Jacobians and empirical original-unit coefficient
# IFs (n x p, rownames = id) at nu = 1, from closed-form M-estimator formulas:
#   L_i = n * Hsum^{-1} * score_i (logit: Hsum = sum p(1-p) x x'; LS: sum x x').
ref_inputs <- function(fx) {
  n <- fx$n
  ids <- as.character(fx$id)
  nu <- rep(1, n)
  f <- ref_fits(fx, nu)
  p <- ref_predictions(fx, f)
  role <- which(fx$A == 1L)
  lab <- function(m, cn) {
    m <- as.matrix(m)
    dimnames(m) <- list(ids, cn)
    m
  }
  # treatment logit
  pt <- 1 / (1 + exp(-drop(fx$Xt %*% f$gam)))
  Ht <- crossprod(fx$Xt, fx$Xt * (pt * (1 - pt)))
  L_gam <- lab(n * (fx$Xt * (fx$G - pt)) %*% solve(Ht), c("g0", "g1", "g2"))
  J_h <- lab(p$h * fx$Xt, c("g0", "g1", "g2"))
  # pooled pair logit
  Hq <- matrix(0, 3, 3)
  Sq <- matrix(0, n, 3)
  for (c in 1:2) {
    B <- ref_pair_design(fx, c)
    qq <- p$q[[c]]
    Hq <- Hq + crossprod(B[role, ], B[role, ] * (qq[role] * (1 - qq[role])))
    Sq[role, ] <- Sq[role, ] + B[role, ] * (fx$Ok[[c]][role] - qq[role])
  }
  L_eta <- lab(n * Sq %*% solve(Hq), c("a1", "a2", "b"))
  J_q <- lapply(1:2, function(c) lab(p$q[[c]] * (1 - p$q[[c]]) * ref_pair_design(fx, c), c("a1", "a2", "b")))
  # pooled imputation LS
  Gj <- matrix(0, 4, 4)
  Sj <- matrix(0, n, 4)
  for (c in 1:2) {
    B <- ref_imp_design(fx, c)
    rr <- role[fx$Ok[[c]][role] == 1L]
    Gj <- Gj + crossprod(B[rr, ])
    Sj[rr, ] <- Sj[rr, ] + B[rr, ] * (fx$dyk[[c]][rr] - p$j[[c]][rr])
  }
  L_th <- lab(n * Sj %*% solve(Gj), c("a1", "a2", "d1", "d2"))
  J_j <- lapply(1:2, function(c) lab(ref_imp_design(fx, c), c("a1", "a2", "d1", "d2")))
  # date logits (separate blocks)
  D <- ref_date_design(fx)
  date_if <- function(qd, Od) {
    H <- crossprod(D[role, ], D[role, ] * (qd[role] * (1 - qd[role])))
    S <- matrix(0, n, 2)
    S[role, ] <- D[role, ] * (Od[role] - qd[role])
    lab(n * S %*% solve(H), c("a", "b"))
  }
  L_s <- date_if(p$qs, as.numeric(fx$Ob[, 2L]))
  L_u <- date_if(p$qu, as.numeric(fx$Ob[, 3L]))
  # chain rule for the product q = qs * qu, written by the caller
  J_prod_s <- lab(p$qu * p$qs * (1 - p$qs) * D, c("a", "b"))
  J_prod_u <- lab(p$qs * p$qu * (1 - p$qu) * D, c("a", "b"))
  # shared date block: the two date contributions are summed into one key
  Dsh_s <- cbind(1, 0, fx$x1)
  Dsh_u <- cbind(0, 1, fx$x1)
  Hsh <- crossprod(Dsh_s[role, ], Dsh_s[role, ] * (p$qs_sh[role] * (1 - p$qs_sh[role]))) +
    crossprod(Dsh_u[role, ], Dsh_u[role, ] * (p$qu_sh[role] * (1 - p$qu_sh[role])))
  Ssh <- matrix(0, n, 3)
  Ssh[role, ] <- Dsh_s[role, ] * (fx$Ob[role, 2L] - p$qs_sh[role]) +
    Dsh_u[role, ] * (fx$Ob[role, 3L] - p$qu_sh[role])
  L_sh <- lab(n * Ssh %*% solve(Hsh), c("as", "au", "b"))
  J_prod_sh <- lab(p$qu_sh * p$qs_sh * (1 - p$qs_sh) * Dsh_s +
                     p$qs_sh * p$qu_sh * (1 - p$qu_sh) * Dsh_u, c("as", "au", "b"))
  list(
    pred = p, fits = f,
    J_h = J_h, L_gam = L_gam,
    J_q = J_q, L_eta = L_eta,
    J_j = J_j, L_th = L_th,
    J_prod_s = J_prod_s, J_prod_u = J_prod_u, L_s = L_s, L_u = L_u,
    J_prod_sh = J_prod_sh, L_sh = L_sh
  )
}
