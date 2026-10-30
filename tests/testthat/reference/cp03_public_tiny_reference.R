# Independent raw-data oracle: no cdid functions or fitted production values.
# Potential outcomes: Y0_it = id/4 + t/2 + epsilon_it; Y1_it = Y0_it + 3/2
# from g=4. Errors have mean zero on symmetric support {-2,-1,-.4,0,.4,1,2};
# this is one fixed realization, not a simulation/coverage experiment.
cp03_public_tiny_data <- function() {
  eps <- rbind(c(0,.4,-1,2), c(1,-1,2,0), c(-1,2,0,1), c(2,0,1,-1),
    c(1,0,2,-1), c(0,2,-1,1), c(2,-1,1,0), c(-1,1,0,2))
  d <- expand.grid(id = 1:8, t = 1:4)
  d$g <- ifelse(d$id <= 4, 4, 0)
  d$y0 <- d$id/4 + d$t/2 + eps[cbind(d$id, d$t)]
  d$tau <- ifelse(d$g == 4 & d$t >= 4, 1.5, 0)
  d$y <- d$y0 + d$tau
  d$weight <- 1
  d
}

cp03_public_tiny_reference <- function(weighting = c("identity", "two_step", "ridge"), alpha = 1e-4) {
  weighting <- match.arg(weighting)
  d <- cp03_public_tiny_data(); Y <- matrix(d$y, 8, 4)
  pairs <- rbind(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
  D <- Y[, pairs[,2], drop=FALSE] - Y[, pairs[,1], drop=FALSE]
  m <- colMeans(D[1:4,,drop=FALSE]) - colMeans(D[5:8,,drop=FALSE])
  phi <- rbind(2*sweep(D[1:4,,drop=FALSE], 2, colMeans(D[1:4,,drop=FALSE]), "-"),
    -2*sweep(D[5:8,,drop=FALSE], 2, colMeans(D[5:8,,drop=FALSE]), "-"))
  V <- crossprod(phi)/8
  # Columns: the free pre-start differential beta_1 and ATT_4. c=2, b=3.
  W <- rbind(c(-1,0), c(-1,0), c(-1,1), c(0,0), c(0,1), c(0,1))
  E <- rbind(c(1,0,0), c(1,1,0), c(1,1,1), c(0,1,0), c(0,1,1), c(0,0,1))
  basis <- c(1,4,6)
  if (weighting == "two_step") {
    A <- W[basis,,drop=FALSE]; Q <- solve(V[basis,basis]); mu <- m[basis]
    Lr <- Q %*% A %*% solve(crossprod(A,Q%*%A))
    Lambda <- matrix(0,6,2); Lambda[basis,] <- Lr
  } else {
    Q <- if (weighting == "identity") diag(6) else solve(V + alpha*sum(diag(V))/6*diag(6))
    Lambda <- Q %*% W %*% solve(crossprod(W,Q%*%W))
  }
  IF <- phi %*% Lambda
  list(data=d, pairs=pairs, W=W, E=E, m=m, phi=phi, V=V, Lambda=Lambda,
    point=drop(crossprod(Lambda,m)), IF=IF, covariance=crossprod(IF)/64,
    lambda=if (weighting == "ridge") alpha*sum(diag(V))/6 else 0)
}
