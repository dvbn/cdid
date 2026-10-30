# Independent shared-nuisance reference for the CP02 weighted-pair IPW link.
# Plain dense arrays over A-units; NO package function, glm.fit or kernel cache.
# Implements the accepted weighted-pair IF (cp02_weighted_pair_if.tex) with
# mechanical case weights w, a local Newton logit solve, and the per-date class
# potentials of Prop. 5 (cp02_gmm_computation.tex). Newton settings below are
# test-fixture choices only and define no production solver policy.

ref_shared_nuisance <- function(S, G, Z, w, g, K, tol = 1e-15, maxit = 60L) {
  n <- length(G)
  cohorts <- c(g, K)
  I <- vapply(cohorts, function(k) as.numeric(G == k), numeric(n))
  if (!is.matrix(I)) I <- matrix(I, nrow = n)
  R <- rowSums(I)
  D <- I[, 1L]
  Cc <- R - D
  q <- colSums(I * w) / n
  r <- colSums(I * w * S) / colSums(I * w)
  rG <- ifelse(R > 0, drop(I %*% r), 1)
  pw <- w * R * S / rG
  beta <- numeric(ncol(Z))
  for (it in seq_len(maxit)) {
    p <- 1 / (1 + exp(-drop(Z %*% beta)))
    grad <- colSums(Z * (pw * (D - p)))
    info <- crossprod(Z, Z * (pw * p * (1 - p)))
    step <- solve(info, grad)
    beta <- beta + step
    if (max(abs(step)) <= tol * (1 + max(abs(beta)))) break
  }
  odds <- exp(drop(Z %*% beta))
  p <- odds / (1 + odds)
  H <- crossprod(Z, Z * (pw * p * (1 - p))) / n
  L <- vapply(seq_along(cohorts), function(j)
    -colSums(Z * (w * I[, j] * S * (D - p))) / n / r[j]^2, numeric(ncol(Z)))
  if (!is.matrix(L)) L <- matrix(L, nrow = ncol(Z))
  xi_r <- I * w * (S - rep(r, each = n)) / rep(q, each = n)
  u <- Z * (w * R * S * (D - p) / rG)
  xi_beta <- t(solve(H, t(u + xi_r %*% t(L))))
  list(n = n, S = S, w = w, I = I, D = D, Cc = Cc, rG = rG, r = r, q = q,
       beta = beta, odds = odds, H = H, L = L, xi_r = xi_r, xi_beta = xi_beta, Z = Z)
}

# Accepted complete link IF for outcome-difference vector `delta` (0 if S = 0).
ref_link_if <- function(nu, delta) {
  S <- nu$S
  a_p <- nu$w * nu$D * S / nu$r[1L]
  a_m <- nu$w * nu$Cc * S * nu$odds / nu$rG
  dp <- mean(a_p)
  dm <- mean(a_m)
  Tm <- sum(a_p * delta) / sum(a_p)
  Mm <- sum(a_m * delta) / sum(a_m)
  v <- colSums(nu$Z * (a_m * (delta - Mm))) / nu$n / dm
  m <- ncol(nu$I)
  gamma <- c(0, vapply(seq_len(m)[-1L], function(j)
    mean(nu$w * nu$I[, j] * S * nu$odds * (delta - Mm)) / (nu$r[j]^2 * dm), numeric(1L)))
  phi <- a_p / dp * (delta - Tm) - a_m / dm * (delta - Mm) +
    drop(nu$xi_r %*% gamma) - drop(nu$xi_beta %*% v)
  list(estimate = Tm - Mm, treated_mean = Tm, control_mean = Mm, phi = phi)
}

# Prop. 5 potential of one date: identical algebra with Y_u in place of Delta.
ref_date_potential <- function(nu, y_u) {
  y <- ifelse(nu$S == 1, y_u, 0)
  out <- ref_link_if(nu, y)
  list(lambda = out$estimate, psi = out$phi)
}
