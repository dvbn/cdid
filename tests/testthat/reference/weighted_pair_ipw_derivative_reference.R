# Independent normalized weighted-pair IPW M-estimation derivative reference.
# A finite-sample derivative check, not population identification or MC evidence.
# Test-only glm controls do not select a production tolerance or solver policy.
set.seed(141026)
n <- 240L
x <- rnorm(n)
Z <- cbind(intercept = 1, x = x)
G <- sample(c(3,4,5,Inf), n, replace = TRUE, prob = c(.32,.26,.15,.27))
K <- c(4,Inf)
risk <- G %in% c(3,K)
D <- as.numeric(G == 3)
C <- as.numeric(G %in% K)
patterns <- list(c(1,2,3),c(1,2),c(1,3),c(2,3),c(1,4),c(3,4))
patt <- sample(seq_along(patterns),n,replace=TRUE)
observed <- matrix(FALSE,n,4)
for (i in seq_len(n)) observed[i,patterns[[patt[i]]]] <- TRUE
stopifnot(all(rowSums(observed)>=2))
alpha <- rnorm(n)
pre <- ifelse(G==3,2,ifelse(G==4,-1,ifelse(G==5,.5,1)))
Y <- cbind(alpha,alpha+pre+.3*x+rnorm(n),
           alpha+pre+1+.8*x+rnorm(n),alpha+pre+2+.5*x+rnorm(n))
Y[G==3,3] <- Y[G==3,3]+1+.4*x[G==3]
cohorts <- c(3,K)
gi <- match(G,cohorts)
fit_link <- function(ell, case=rep(1,n), with_if=FALSE) {
  S <- as.numeric(observed[,ell[1]] & observed[,ell[2]])
  delta <- ifelse(S==1,Y[,ell[2]]-Y[,ell[1]],0)
  group_mass <- vapply(cohorts,function(k) sum(case[G==k]),numeric(1))
  rates <- vapply(cohorts,function(k) sum(case[G==k]*S[G==k]),numeric(1))/group_mass
  rr <- rep(1,n);rr[risk] <- rates[gi[risk]]
  use <- risk & S==1
  fit <- withCallingHandlers(
    glm.fit(Z[use,,drop=FALSE],D[use],weights=case[use]/rr[use],
            family=binomial(),control=glm.control(epsilon=1e-12,maxit=100)),
    warning=function(e) {
      if (grepl("non-integer #successes",conditionMessage(e),fixed=TRUE))
        invokeRestart("muffleWarning") else stop(e)
    })
  stopifnot(fit$converged,all(is.finite(fit$coefficients)))
  odds <- exp(drop(Z %*% fit$coefficients));prob <- odds/(1+odds)
  a <- D*S/rr;b <- C*S*odds/rr
  Tmean <- sum(case*a*delta)/sum(case*a)
  Cmean <- sum(case*b*delta)/sum(case*b)
  result <- list(value=Tmean-Cmean)
  if (!with_if) return(result)
  stopifnot(all(case==1))
  q <- group_mass/n
  rate_if <- vapply(seq_along(cohorts),function(j)
    as.numeric(G==cohorts[j])*(S-rates[j])/q[j],numeric(n))
  score <- Z*as.numeric(risk*S/rr*(D-prob))
  H <- crossprod(Z,Z*as.numeric(risk*S/rr*prob*(1-prob)))/n
  A_beta_rate <- vapply(seq_along(cohorts),function(j)
    -colMeans(Z*as.numeric((G==cohorts[j])*S/rates[j]^2*(D-prob))),numeric(ncol(Z)))
  beta_if <- t(solve(H,t(score+rate_if %*% t(A_beta_rate))))
  dbeta <- -colSums(Z*as.numeric(b*(delta-Cmean)))/sum(b)
  drate <- c(0,vapply(seq_along(K),function(j) {
    k <- K[j];rj <- rates[j+1L]
    sum(as.numeric(G==k)*S*odds*(delta-Cmean))/(rj^2*sum(b))
  },numeric(1)))
  direct <- a/mean(a)*(delta-Tmean)-b/mean(b)*(delta-Cmean)
  phi <- direct+drop(rate_if %*% drate)+drop(beta_if %*% dbeta)
  beta_if_omitted_rate <- t(solve(H,t(score)))
  wrong <- direct+drop(rate_if %*% drate)+drop(beta_if_omitted_rate %*% dbeta)
  stopifnot(max(abs(phi[!risk]))==0,abs(mean(phi))<1e-10)
  list(value=result$value,phi=phi,wrong=wrong,S=S,
       rate_partial=drate,beta_partial=dbeta)
}
pairs <- list(c(1,2),c(1,3))
fits <- lapply(pairs,fit_link,with_if=TRUE)
phi <- do.call(cbind,lapply(fits,`[[`,"phi"))
finite <- function(h) vapply(seq_len(n),function(i) {
  plus <- minus <- rep(1,n);plus[i]<-1+h;minus[i]<-1-h
  n*vapply(pairs,function(ell)
    (fit_link(ell,plus)$value-fit_link(ell,minus)$value)/(2*h),numeric(1))
},numeric(2))
d1 <- t(finite(1e-4));d2 <- t(finite(5e-5))
err <- max(abs(phi-d2));resolution <- max(abs(d1-d2))
wrong_err <- max(abs(do.call(cbind,lapply(fits,`[[`,"wrong"))-d2))
stopifnot(err<1e-6,resolution<1e-6,wrong_err>1e-4,
          max(abs((phi[,2]-phi[,1])-(d2[,2]-d2[,1])))<1e-6,
          any(abs(fits[[1]]$phi[fits[[1]]$S==0 & risk])>1e-4))
covariance <- crossprod(sweep(phi,2,colMeans(phi)))/n
stopifnot(abs(covariance[1,2])>1e-4)
cat("Weighted-pair IPW derivative; 240 units, all observed at least twice.\n")
cat("Pairs12/13; propensity fitted on inverse-rate-weighted observed pairs.\n")
cat("Maximum complete-IF versus refitted case-weight derivative error:",err,"\n")
cat("Maximum change across derivative step sizes:",resolution,"\n")
cat("Error omitting rate effects inside propensity-score Jacobian:",wrong_err,"\n")
cat("Joint link covariance off-diagonal:",covariance[1,2],"\n")
cat("Pair non-observer rate/nuisance contributions detected; outside-risk rows zero.\n")
cat("PASS: full signed normalization/rate/propensity IF and bridge difference derivative.\n")
cat("This checks differentiation only, not statistical identification, inference or MC coverage.\n")
