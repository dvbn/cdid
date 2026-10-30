# Independent company CP03-B reference. Base R only; no production helpers.
# Rows below are literal date-pair differences, in the stated fixture order.
cp03_company_dense_fit <- function(W, m, V, Q = solve(V)) {
  A <- solve(crossprod(W, Q %*% W), crossprod(W, Q))
  list(point = drop(A %*% m), Lambda = t(A), covariance = A %*% V %*% t(A))
}

cp03_company_dense_reference <- function() {
  W <- rbind(c(-1, 1, 0), c(0, -1, 0), c(0, 0, 1))
  H <- rbind(c(1, 0), c(0, 0), c(0, 1))
  WH <- rbind(c(-1, 0), c(0, 0), c(0, 1))
  V <- rbind(c(4, -1, 0), c(-1, 3, -1), c(0, -1, 2))
  m <- c(3/2, 1/2, 5/2)
  zero <- cp03_company_dense_fit(WH, m, V)
  dropped <- cp03_company_dense_fit(WH[c(1, 3), ], m[c(1, 3)], V[c(1, 3), c(1, 3)])
  identity <- cp03_company_dense_fit(WH, m, V, diag(3))

  # Full rows 12,13,23,34; retained rows 12,13,34. Row 23 = -12 + 13.
  E <- rbind(c(1, 0, 0), c(0, 1, 0), c(-1, 1, 0), c(0, 0, 1))
  B <- rbind(c(1, 0, 0), c(1, 1, 0), c(0, 1, 0), c(0, 0, 1))
  Vbasis <- rbind(c(2, 0, 1), c(0, 2, 1), c(1, 1, 1))
  Vfull <- B %*% Vbasis %*% t(B)
  mfull <- drop(B %*% c(1, 1, 2))
  Wfull <- rbind(c(-1, 0), c(-1, 0), c(0, 0), c(0, 1))
  keep <- c(1L, 2L, 4L)
  Vr <- Vfull[keep, keep]
  Wr <- Wfull[keep, ]
  mr <- mfull[keep]
  lambda <- (1/4) * sum(diag(Vfull)) / 4
  Q <- solve(Vfull + diag(lambda, 4))
  Qr <- crossprod(E, Q %*% E)
  ridge <- cp03_company_dense_fit(Wfull, mfull, Vfull, Q)
  ridge_reduced <- cp03_company_dense_fit(Wr, mr, Vr, Qr)
  wrong_lambda <- (1/4) * sum(diag(Vr)) / 3
  wrong_ridge <- cp03_company_dense_fit(Wr, mr, Vr, solve(Vr + diag(wrong_lambda, 3)))

  # Four fixed units, links 12,13,23,34, balanced intercept-only shared IPW.
  # Treated endpoint differences minus their treated means, and negative
  # control differences minus their control means; multiply each by 2.
  shared_phi <- rbind(c(0, 0, 0, -1), c(0, 0, 0, 1),
                      c(2, 2, 0, 1), c(-2, -2, 0, -1))
  # The fixture's custom callback returns psi=(2D-1)*dy minus its sample mean.
  custom_phi <- rbind(c(1, 5/2, 3/2, 2), c(1, 5/2, 3/2, 3),
                      c(0, -3/2, -3/2, -2), c(-2, -7/2, -3/2, -3))
  # Columns: anticipation2015, ATT2017, ATT2020; rows 2010-12,12-15,12-17,17-20.
  irregular <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, -1, 1))
  # Columns: private1, ATT4, ATT5; rows 12,23,24,35; anchors 2 and 3.
  multi_anchor <- rbind(c(-1, 0, 0), c(0, 0, 0), c(0, 1, 0), c(0, 0, 1))
  # Rows 13,24: the second component is grounded at 2 only after restriction.
  newly_post <- rbind(c(-1, 0, 0), c(0, -1, 1))
  newly_all <- rbind(c(-1, 0), c(0, 1))
  list(W = W, H = H, WH = WH, V = V, m = m, zero = zero, dropped = dropped,
       identity = identity, E = E, G = crossprod(E), Wfull = Wfull,
       Vfull = Vfull, mfull = mfull, lambda = lambda, ridge = ridge,
       ridge_reduced = ridge_reduced, wrong_ridge = wrong_ridge,
       shared_phi = shared_phi, custom_phi = custom_phi,
       irregular = irregular, multi_anchor = multi_anchor,
       newly_post = newly_post, newly_all = newly_all)
}

cp03_company_reference_check <- function() {
  r <- cp03_company_dense_reference()
  near <- function(x, y) stopifnot(max(abs(x - y)) < 1e-13)
  stopifnot(identical(r$W %*% r$H, r$WH))
  near(r$zero$point, c(-5/3, 8/3))
  near(r$zero$covariance, rbind(c(11/3, 1/3), c(1/3, 5/3)))
  near(r$zero$Lambda, rbind(c(-1, 0), c(-1/3, 1/3), c(0, 1)))
  near(r$dropped$covariance[2, 2], 2)
  near(r$identity$point, c(-3/2, 5/2))
  near(r$lambda, 9/16)
  near(r$ridge$point, r$ridge_reduced$point)
  near(r$ridge$covariance, r$ridge_reduced$covariance)
  near(crossprod(r$E, r$ridge$Lambda), r$ridge_reduced$Lambda)
  stopifnot(max(abs(r$wrong_ridge$point - r$ridge$point)) > 0.02)
  near(t(r$E %*% t(r$shared_phi[, c(1, 2, 4)])), r$shared_phi)
  near(t(r$E %*% t(r$custom_phi[, c(1, 2, 4)])), r$custom_phi)
  invisible(r)
}
