# Independent finite weighted-mean reference: two dates, an intercept-only
# control mean, NT controls. Each cohort's common pair rate cancels exactly.
# No production helper, fitted value, random draw or did equality is used.
cp06_mean_reference <- function(raw, dates = c(2, 4), cohorts = c(3, 4)) {
  N <- nrow(raw); S <- !is.na(raw[[paste0("Y", dates[1])]]) & !is.na(raw[[paste0("Y", dates[2])]])
  dy <- raw[[paste0("Y", dates[2])]] - raw[[paste0("Y", dates[1])]]
  mean_if <- function(g) {
    hit <- S & raw$cohort == g; mass <- sum(raw$w[hit])
    mu <- sum(raw$w[hit]*dy[hit])/mass
    phi <- numeric(N); phi[hit] <- N*raw$w[hit]*(dy[hit]-mu)/mass
    list(mu = mu, phi = phi)
  }
  control <- mean_if(Inf); targets <- lapply(cohorts, mean_if)
  theta <- vapply(targets, function(z) z$mu-control$mu, numeric(1))
  phi <- vapply(targets, function(z) z$phi-control$phi, numeric(N))
  mass <- vapply(cohorts, function(g) sum(raw$w[raw$cohort == g]), numeric(1))
  shares <- mass/sum(mass); overall <- sum(shares*theta)
  overall_if <- as.vector(phi%*%shares)
  for (j in seq_along(cohorts)) {
    hit <- raw$cohort == cohorts[j]
    overall_if[hit] <- overall_if[hit]+N*raw$w[hit]*(theta[j]-overall)/sum(mass)
  }
  list(theta = theta, phi = phi, overall = overall, overall_if = overall_if,
    dynamic = mean(theta), dynamic_if = rowMeans(phi),
    covariance = crossprod(scale(phi, scale = FALSE))/N^2)
}
