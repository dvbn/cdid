# Independent finite U-functional refits. No production code/IF/Jacobians.
cp04_finite_u <- function(ref, U, rho, cw = rep(1, ref$N), Z = ref$Z, V = ref$V) {
  N <- ref$N; ww <- ref$w * cw; nc <- sum(cw)
  I <- ref$I; A <- ref$A; m <- colSums(ww * I); mu <- colSums(ww * U * I)
  ma <- colSums(ww * A * I); a <- ma / mu
  al <- list(); pr <- matrix(NA_real_, N, 2L)
  for (e in 1:2) for (k in 1:3) {
    rr <- which(U & ref$dat$k == k)
    alpha <- ref$fit_logit(Z[rr, , drop = FALSE], ref$S[rr,e], ww[rr])
    al[[length(al) + 1L]] <- alpha
    pr[rr,e] <- plogis(drop(Z[rr, , drop = FALSE] %*% alpha))
  }
  be <- list(); T <- M <- dp <- dm <- numeric(3L)
  for (l in 1:3) {
    link <- ref$links[[l]]; e <- link$e
    rr <- which(ref$S[,e] == 1 & ref$dat$k %in% c(1L, link$K))
    D <- as.double(ref$dat$k[rr] == 1L)
    fw <- ww[rr] * a[ref$dat$k[rr]] / (rho[rr] * pr[rr,e])
    beta <- ref$fit_logit(V[rr, , drop = FALSE], D, fw); be[[l]] <- beta
    od <- exp(drop(V[rr, , drop = FALSE] %*% beta))
    pair <- ref$pairs[[e]]; dy <- ref$Y[rr,pair[2]] - ref$Y[rr,pair[1]]
    wp <- fw * D; wm <- fw * (1-D) * od
    dp[l] <- sum(wp)/nc; dm[l] <- sum(wm)/nc
    T[l] <- sum(wp * dy)/sum(wp); M[l] <- sum(wm * dy)/sum(wm)
  }
  list(vector = c(a, unlist(al), unlist(be), T, M, dp, dm), phi_point = T - M,
    a = a, alpha = al, beta = be, T = T, M = M, dp = dp, dm = dm)
}
