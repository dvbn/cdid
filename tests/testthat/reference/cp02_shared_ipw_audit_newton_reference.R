# Independent Newton logit and complete weighted-pair IF reference for the
# packaged shared-IPW audit tests. Uses only base R and stats; it is sourced into an
# environment whose parent is baseenv() and never calls cdid code. Expected
# values are computed here, never by the candidate.
independent_fit <- function(fx,s,t,w=fx$w) {
  n <- fx$n; K <- c(4,7,Inf)
  I <- outer(fx$G,K,'==')*1
  S <- as.numeric(fx$obs[,s]&fx$obs[,t]&rowSums(I)>0)
  D <- I[,1]; C <- rowSums(I)-D
  q <- colSums(I*w)/n; r <- colSums(I*w*S)/colSums(I*w)
  rG <- drop(I%*%r); rG[rG==0] <- 1
  pw <- w*S/rG; Z <- fx$Z; beta <- numeric(ncol(Z))
  for(k in 1:50) {
    p <- stats::plogis(drop(Z%*%beta))
    step <- solve(crossprod(Z,Z*(pw*p*(1-p))),colSums(Z*(pw*(D-p))))
    beta <- beta+step
    if(max(abs(step))<1e-13) break
  }
  p <- stats::plogis(drop(Z%*%beta)); odds <- exp(drop(Z%*%beta))
  H <- crossprod(Z,Z*(pw*p*(1-p)))/n
  L <- vapply(seq_along(K),function(j) -colSums(Z*(w*I[,j]*S*(D-p)))/n/r[j]^2,numeric(ncol(Z)))
  delta <- fx$Y[,t]-fx$Y[,s]
  a <- w*D*S/r[1]; b <- w*C*S*odds/rG
  T <- sum(a*delta)/sum(a); M <- sum(b*delta)/sum(b)
  Gbeta <- -colSums(Z*(b*(delta-M)))/sum(b)
  Gr <- c(0,vapply(2:3,function(j)sum(w*I[,j]*S*odds*(delta-M))/r[j]^2/sum(b),numeric(1)))
  lb <- solve(H,Gbeta); lr <- Gr+drop(crossprod(L,lb))
  xir <- I*w*(S-rep(r,each=n))/rep(q,each=n)
  components <- list(treated=a/mean(a)*(delta-T),control=-b/mean(b)*(delta-M),beta_score=pw*(D-p)*drop(Z%*%lb),rate=drop(xir%*%lr))
  list(estimate=T-M,beta=beta,components=components,phi=Reduce('+',components),r=r,H=H,L=L,lb=lb,lr=lr)
}
