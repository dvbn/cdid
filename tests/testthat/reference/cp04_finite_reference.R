# CP04 deterministic finite-functional reference. Base R only; no RNG or package code.
# Test-local independent refitter; sourced into a separate environment by the installed tests.

cat("CP04 direct frame-logit independent finite reference\n")
cat("R:", R.version.string, "\nNo RNG, sampled DGP, package source, or inference simulation.\n")
tol <- 3e-7                         # fixed before execution; all coordinates
root_tol <- 2e-11
patterns <- rbind(c(1,1,1,1), c(1,1,0,0), c(1,0,1,0),
                  c(0,1,0,1), c(0,0,1,1), c(1,0,0,1),
                  c(1,0,0,0), c(0,1,0,0), c(0,0,0,0))
dat <- do.call(rbind, lapply(1:3, function(k) do.call(rbind,
  lapply(1:3, function(j) {
    counts <- c(2+(k+j)%%2, 1+j%%2, 1+(k==1), 1+(k==2),
                1+(j==2), 1+(k==3), 1+k, 1+(j==3), 1+2*(k-1))
    do.call(rbind, lapply(1:9, function(h)
      data.frame(k=k, x=j-2, h=h, r=seq_len(counts[h]))))
  }))))
N <- nrow(dat); ids <- seq_len(N)
w <- 0.6 + (1+(7*ids)%%11)/9
obs <- patterns[dat$h,,drop=FALSE]; A <- as.numeric(rowSums(obs)>=2)
Y <- sapply(1:4, function(t) 0.2*dat$k + t*(dat$x+0.31*dat$k) +
              t^2*(0.07*dat$x^2+0.09*dat$k) +
              (t%%2)*(dat$r-1)*0.13 + t*dat$h*0.025)
Y[obs==0] <- NA_real_               # unobserved differences cannot be used
Z <- cbind(intercept=1, x=dat$x); V <- Z
pairs <- list(c(1L,2L), c(2L,4L))
S <- sapply(pairs, function(e) obs[,e[1]]*obs[,e[2]])
links <- list(list(e=1L,K=c(2L,3L),name="NYT12"),
              list(e=2L,K=c(2L,3L),name="NYT24"),
              list(e=1L,K=3L,name="NT12"))
I <- sapply(1:3, function(k) as.numeric(dat$k==k))
qsolve <- function(H, B) qr.solve(H, B, tol=1e-12)
fit_logit <- function(X, D, ww) {
  # Actual finite equation solve; deliberately independent of analytic IF code.
  f <- withCallingHandlers(stats::glm.fit(x=X, y=D, weights=ww,
    family=stats::binomial(),
    control=stats::glm.control(epsilon=1e-13,maxit=60)), warning=function(e) {
      # Analytic case weights are intentionally fractional; no other warning is ignored.
      if (identical(conditionMessage(e),"non-integer #successes in a binomial glm!"))
        invokeRestart("muffleWarning")
    })
  stopifnot(f$converged, f$rank==ncol(X), all(is.finite(f$coefficients)))
  pr <- plogis(drop(X %*% f$coefficients))
  stopifnot(all(pr>0 & pr<1), max(abs(crossprod(X,ww*(D-pr))))<root_tol)
  unname(f$coefficients)
}
finite <- function(cw=rep(1,N)) {
  # No analytic derivatives, Jacobians, IFs, or supplied private fitter are called.
  ww <- w*cw; nc <- sum(cw)
  m <- colSums(ww*I); ma <- colSums(ww*A*I); a <- ma/m
  al <- array(0,c(2,3,2)); pr <- matrix(0,N,2)
  for (e in 1:2) for (k in 1:3) {
    ii <- which(dat$k==k)
    al[,k,e] <- fit_logit(Z[ii,,drop=FALSE],S[ii,e],ww[ii])
    pr[ii,e] <- plogis(drop(Z[ii,,drop=FALSE] %*% al[,k,e]))
  }
  out <- list(a=a,al=al,pr=pr, beta=list(), T=numeric(3), M=numeric(3),
              dp=numeric(3),dm=numeric(3))
  for (l in 1:3) {
    e <- links[[l]]$e; C <- as.numeric(dat$k %in% links[[l]]$K)
    D <- I[,1]; R <- C+D; ii <- which(R*S[,e]==1)
    b <- a[dat$k]/pr[,e]
    beta <- fit_logit(V[ii,,drop=FALSE],D[ii],ww[ii]*b[ii])
    odds <- exp(drop(V %*% beta)); dy <- rep(0,N)
    jj <- which(S[,e]==1); pp <- pairs[[e]]
    dy[jj] <- Y[jj,pp[2]]-Y[jj,pp[1]]
    wp <- ww*D*S[,e]*b; wm <- ww*C*S[,e]*b*odds
    out$beta[[l]] <- beta
    out$dp[l] <- sum(wp)/nc; out$dm[l] <- sum(wm)/nc
    out$T[l] <- sum(wp*dy)/sum(wp); out$M[l] <- sum(wm*dy)/sum(wm)
  }
  out$r <- colSums(ww*S[,1]*I)/ma
  out$pi <- colSums(ww*S[,1]*I)/m
  out$eta <- ma/sum(ma); out$tau <- m/sum(m); out$lambda <- sum(cw*A)/nc
  out$vector <- c(a=a,alpha=as.vector(al),beta=unlist(out$beta),
                  T=out$T,M=out$M,dp=out$dp,dm=out$dm,
                  eta=out$eta,tau=out$tau,lambda=out$lambda,r=out$r,pi=out$pi,
                  dp_A=out$dp/out$lambda,dm_A=out$dm/out$lambda)
  out
}
f <- finite()
# Analytic side: only baseline roots enter; all derivatives are explicit below.
mk <- colMeans(w*I); mak <- colMeans(w*A*I)
xa <- w*I*sweep(matrix(A,N,3),2,f$a,"-")
xa <- sweep(xa,2,mk,"/")
xal <- matrix(0,N,12)
for (e in 1:2) for (k in 1:3) {
  ix <- (e-1)*6+(k-1)*2+1:2; pp <- f$pr[,e]
  H <- crossprod(Z,w*I[,k]*pp*(1-pp)*Z)/N
  uu <- w*I[,k]*(S[,e]-pp)*Z
  xal[,ix] <- t(qsolve(H,t(uu)))
}
xb <- matrix(0,N,6); xT <- xM <- xdp <- xdm <- matrix(0,N,3)
direct <- matrix(0,N,3); sel_forcing <- matrix(0,N,3)
bad <- matrix(0,N,3)
for (l in 1:3) {
  e <- links[[l]]$e; D <- I[,1]; C <- as.numeric(dat$k %in% links[[l]]$K)
  R <- C+D; pp <- f$pr[,e]; b <- f$a[dat$k]/pp
  h <- plogis(drop(V %*% f$beta[[l]])); odds <- exp(drop(V %*% f$beta[[l]]))
  dy <- rep(0,N); jj <- which(S[,e]==1); ee <- pairs[[e]]
  dy[jj] <- Y[jj,ee[2]]-Y[jj,ee[1]]
  wp <- w*D*S[,e]*b; wm <- w*C*S[,e]*b*odds
  uu <- w*R*S[,e]*b*(D-h)*V
  H <- crossprod(V,w*R*S[,e]*b*h*(1-h)*V)/N
  Ca <- sapply(1:3,function(k) colMeans(I[,k]*uu)/f$a[k])
  Cal <- matrix(0,2,12)
  da <- sweep(I,2,f$a,"/"); dal <- matrix(0,N,12)
  for (k in 1:3) {
    ix <- (e-1)*6+(k-1)*2+1:2
    dal[,ix] <- -I[,k]*(1-pp)*Z
    Cal[,ix] <- crossprod(uu,-I[,k]*(1-pp)*Z)/N
  }
  xbeta <- t(qsolve(H,t(uu+xa%*%t(Ca)+xal%*%t(Cal))))
  xbeta_bad <- t(qsolve(H,t(uu+xa%*%t(Ca))))
  xb[,(l-1)*2+1:2] <- xbeta
  # Separate derivative of each actual empirical numerator/denominator.
  plus_a <- colMeans(wp*(dy-f$T[l])*da)/f$dp[l]
  minus_a <- colMeans(wm*(dy-f$M[l])*da)/f$dm[l]
  plus_al <- colMeans(wp*(dy-f$T[l])*dal)/f$dp[l]
  minus_al <- colMeans(wm*(dy-f$M[l])*dal)/f$dm[l]
  vb <- colMeans(wm*(dy-f$M[l])*V)/f$dm[l]
  xT[,l] <- wp*(dy-f$T[l])/f$dp[l] + xa%*%plus_a+xal%*%plus_al
  xM[,l] <- wm*(dy-f$M[l])/f$dm[l] + xa%*%minus_a+xal%*%minus_al+xbeta%*%vb
  xdp[,l] <- wp-f$dp[l]+xa%*%colMeans(wp*da)+xal%*%colMeans(wp*dal)
  xdm[,l] <- wm-f$dm[l]+xa%*%colMeans(wm*da)+xal%*%colMeans(wm*dal)+xbeta%*%colMeans(wm*V)
  direct[,l] <- wp*(dy-f$T[l])/f$dp[l]-wm*(dy-f$M[l])/f$dm[l]
  sel_forcing[,l] <- rowSums(abs(xal%*%t(Cal)))
  bad[,l] <- xT[,l]-(wm*(dy-f$M[l])/f$dm[l]+xa%*%minus_a+xal%*%minus_al+xbeta_bad%*%vb)
}
xeta <- w*A*sweep(I,2,f$eta,"-")/mean(w*A)
xtau <- w*sweep(I,2,f$tau,"-")/mean(w)
xr <- sweep(w*I*(S[,1]-A*rep(f$r,each=N)),2,mak,"/")
xpi <- sweep(w*I*sweep(matrix(S[,1],N,3),2,f$pi,"-"),2,mk,"/")
analytic <- cbind(xa,xal,xb,xT,xM,xdp,xdm,xeta,xtau,A-f$lambda,xr,xpi,
                 xdp/f$lambda-outer(A-f$lambda,f$dp/f$lambda^2),
                 xdm/f$lambda-outer(A-f$lambda,f$dm/f$lambda^2))
phi <- xT-xM
stopifnot(ncol(analytic)==length(f$vector), nrow(analytic)==N)
err <- matrix(0,2,3,dimnames=list(c("1e-4","5e-5"),c("stack","link","normalizers")))
last <- NULL
for (j in 1:2) {
  step <- c(1e-4,5e-5)[j]
  fd <- matrix(0,N,length(f$vector))
  for (i in seq_len(N)) {
    cp <- cm <- rep(1,N); cp[i] <- 1+step; cm[i] <- 1-step
    fd[i,] <- N*(finite(cp)$vector-finite(cm)$vector)/(2*step)
  }
  # a:3 alpha:12 beta:6 T:3 M:3 dp:3 dm:3 eta:3 tau:3 lambda:1 r:3 pi:3 dp_A:3 dm_A:3
  fdphi <- fd[,22:24,drop=FALSE]-fd[,25:27,drop=FALSE]
  err[j,] <- c(max(abs(fd-analytic)), max(abs(fdphi-phi)),
               max(abs(fd[,c(28:33,47:52)]-analytic[,c(28:33,47:52)])))
  stopifnot(all(err[j,]<tol))
  if (!is.null(last)) stopifnot(max(abs(fd-last))<tol)
  last <- fd
}
cat("Frame N =",N,"; repeated n_A =",sum(A),"; stack dimension =",ncol(analytic),"\n")
cat("Fixed tolerances: derivative",tol,"; logit score",root_tol,"\n")
cat("Weighted cohort repeat rates:",format(f$a),"\n")
cat("Link estimates (NYT12, NYT24, NT12):",format(f$T-f$M),"\n")
cat("Max absolute errors versus independent complete case refits:\n"); print(err)
center_error <- max(abs(colMeans(analytic)))
stopifnot(center_error<root_tol)
cat("Uncentered analytic IF mean, max absolute:",center_error,"\n")
stopifnot(max(sel_forcing)>1e-3, max(abs(xdp))>1e-3, max(abs(xdm))>1e-3)
stopifnot(max(abs(phi[A==0,]))>1e-3, max(abs(phi[S[,1]==0,1]))>1e-3)
cat("Max absolute A=0 link IF:",max(abs(phi[A==0,])),"\n")
cat("Max selection -> treatment forcing:",max(sel_forcing),"\n")
stopifnot(max(abs(phi[dat$k==2,3]))<root_tol) # NT has no cohort-5 role
cat("NT unused cohort link rows are zero; relevant NT nonobserver rows retained.\n")
lambda <- sum(A)/N; psi <- lambda*phi
vF <- crossprod(scale(phi,scale=FALSE))/N^2
vA <- crossprod(scale(psi,scale=FALSE))/sum(A)^2
stopifnot(max(abs(vF-vA))<1e-13)
covFD <- crossprod(scale(last[,22:24]-last[,25:27],scale=FALSE))/N^2
stopifnot(max(abs(vF-covFD))<tol)
stopifnot(abs(vF[1,3])>1e-5)
cat("Shared-pair/shared-selection NYT12-NT12 covariance:",vF[1,3],"\n")
cat("Frame versus n_A variance-scaling error:",max(abs(vF-vA)),"\n")
# All frame clusters, including zero-eligible clusters; no inference simulation.
cl <- ifelse(A==0, 1000L+ids, 1L+(ids%%7L))
cs <- rowsum(phi,cl); ca <- rowsum(psi,cl)
stopifnot(max(abs(crossprod(scale(cs,scale=FALSE))/N^2-
                     crossprod(scale(ca,scale=FALSE))/sum(A)^2))<1e-13)
cat("Cluster scaling checked with",nrow(cs),"clusters, including",sum(A==0),"zero-eligible clusters.\n")
# Existing marginal pi/r representation cancels to a, with its derivative.
stopifnot(max(abs(f$pi/f$r-f$a))<1e-13)
cancel_if <- sweep(xpi,2,f$r,"/")-sweep(xr,2,f$pi/f$r^2,"*")
stopifnot(max(abs(cancel_if-xa))<1e-12)
cat("Marginal-rate cancellation and derivative checked.\n")
# Representative-B known-probability cancellation. This does not assert finite
# equality of two separately fitted, nonsaturated logits on F and B.
agB <- c(.61,.72,.83); bg <- c(.8,.6,.5)
pb <- c(.22,.43,.29,.48,.31,.57); kg <- rep(1:3,each=2)
ratioB <- agB[kg]/pb
ratioF <- function(bb) (agB*bb)[kg]/(pb*bb[kg])
stopifnot(max(abs(ratioF(bg)-ratioB))<1e-14)
for (k in 1:3) {
  bp <- bm <- bg; bp[k] <- bp[k]+1e-5; bm[k] <- bm[k]-1e-5
  stopifnot(max(abs((ratioF(bp)-ratioF(bm))/2e-5))<1e-10)
}
xx <- rep(c(-1,1),3); vv <- cbind(1,xx); dd <- as.numeric(kg==1)
known <- function(rr) {
  be <- fit_logit(vv,dd,rr*rep(c(1,1.2),3))
  ww <- rr*rep(c(1,1.2),3); od <- exp(drop(vv%*%be)); yy <- kg+xx*kg/3
  c(be,sum(ww*dd),sum(ww*(1-dd)*od),
    sum(ww*dd*yy)/sum(ww*dd)-sum(ww*(1-dd)*od*yy)/sum(ww*(1-dd)*od))
}
stopifnot(max(abs(known(ratioF(bg))-known(ratioB)))<1e-12)
cat("Representative-B b<1 cancellation: weights, score root, both normalizers and link checked.\n")
# Negative controls must fail at the unchanged positive-check tolerance.
mutant_gap <- max(abs(bad-phi))
zero_nonobserver_gap <- max(abs(phi[A==0,]))
wrong_sign_gap <- max(abs((xT+xM)-phi))
wrong_scale_gap <- max(abs(phi-psi))
stopifnot(min(mutant_gap,zero_nonobserver_gap,wrong_sign_gap,wrong_scale_gap)>100*tol)
cat("Rejected omission of selection -> treatment Jacobian; max IF error:",mutant_gap,"\n")
cat("Rejected zero A=0 IF, positive control sign and unconverted scale; gaps:",
    zero_nonobserver_gap,wrong_sign_gap,wrong_scale_gap,"\n")
cat("PASS: finite derivatives, both normalizers, NYT mixture, NT, shares/rates, shared covariance and scaling.\n")
cat("Limits: deterministic fixed data only; no identification, raw-data adapter, statistical coverage, package or public-API validation.\n")
