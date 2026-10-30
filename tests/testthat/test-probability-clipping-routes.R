# Bounded D93 active-engine supplement. Independent finite refits use explicit
# cell masses and ordinary WLS normal equations; production tables are subjects,
# never oracles. FD 1e-6 is the existing new-test convention. Identity-weight
# public checks are wiring checks, not derivatives of two-step GMM or coverage.
d93r_fn <- function(n) getFromNamespace(n, "cdid")
d93r_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
d93r_eq <- function(a, b, tol = 1e-6) {
  testthat::expect_equal(unname(a), unname(b), tolerance = tol)
}

# Same 32 atoms as the independent IPW reference. z is an outcome predictor
# distinct from the treatment/observation predictor x, avoiding saturated-OR
# cancellation of the changed propensity and selection loadings.
d93r_atoms <- function(bound = TRUE) {
  d <- expand.grid(r1 = 0:1, r2 = 0:1, r3 = 0:1, x = 0:1, k = 1:2)
  p <- c(.3, 4e-4)[d$x + 1L]; p3 <- c(.4, .7)[d$x + 1L]
  mass <- matrix(c(.7, 1.2, 1.5, if (bound) 2e-7 else 2e-5), 2L)
  w <- mass[cbind(d$k, d$x + 1L)] * ifelse(d$r1, p, 1-p) *
    ifelse(d$r2, p, 1-p) * ifelse(d$r3, p3, 1-p3)
  R <- as.matrix(d[c("r1", "r2", "r3")]); storage.mode(R) <- "double"
  dy <- 1 + d$k + 2*d$x + .7*d$r3 + .6*d$r1*d$r3 + .3*d$k*d$x
  Y <- cbind(0, dy, 2*dy + d$r1); Y[R == 0] <- NA_real_
  list(N = nrow(d), x = d$x, z = d$r3 + .2*d$r1, w = nrow(d)*w/sum(w),
    R = R, Y = Y, g = ifelse(d$k == 1, 3, Inf), K = Inf)
}

# Reuse the author's empirical-rate pattern, with moderate control weights:
# r_4 is below eps while every fitted propensity remains below the upper cap.
d93r_rates <- function() {
  pattern <- c("123","12","13","123","12","23","2","12","12",
    "13","23","123","12","13","12","123","23","123")
  R <- t(vapply(pattern, function(p) as.double(vapply(1:3,
    function(s) grepl(as.character(s), p), logical(1))), numeric(3)))
  g <- c(rep(3,7),rep(4,4),rep(Inf,6),5)
  x <- c(0,0,0,1,1,1,0,0,1,0,1,0,0,0,1,1,1,0)
  w <- c(1,1.5,.7,1.2,.8,1.1,1,1e-7,1e-8,1,1,1,2,1,.7,1.3,1.5,1)
  Y <- outer(1:18,1:3,function(i,t) cos(1.7*i+t)+.25*t*(g[i]==3))
  Y[R == 0] <- NA_real_
  list(N=18L,x=x,z=sin(1:18),w=w,R=R,Y=Y,g=g,K=c(4,Inf))
}

d93r_X <- function(f) cbind("(Intercept)" = 1, x = as.double(f$x), z = f$z)
d93r_links <- function(f) list(list(s=1,t=2,g=3,K=f$K,link_id="D93route"))
d93r_wls <- function(X, y, w) as.vector(solve(crossprod(X, w*X), crossprod(X, w*y)))

# Complete independent refit of ONE actual pair, including nuisance roots,
# both normalizers, WLS and the final reconstructed response.
d93r_ref <- function(f, route, cw=rep(1,f$N), model="pair_logit") {
  N <- f$N; w <- f$w*cw; X <- d93r_X(f); O <- X[,c(1,3),drop=FALSE]
  V <- X[,1:2,drop=FALSE]; A <- rowSums(f$R)>=2
  S <- f$R[,1]*f$R[,2]; cohorts <- c(3,f$K); k <- match(f$g,cohorts)
  risk <- !is.na(k); D <- as.double(f$g==3)
  y <- numeric(N); y[S==1] <- f$Y[S==1,2]-f$Y[S==1,1]
  a <- r <- numeric(length(cohorts)); q <- rep(1,N)
  for (j in seq_along(cohorts)) {
    on <- risk & k %in% j
    a[j] <- sum(w[on]*A[on])/sum(w[on])
    r[j] <- sum(w[on]*S[on])/sum(w[on]*A[on])
    for (v in 0:1) {
      cell <- on & f$x==v
      if (!any(cell)) next
      q[cell] <- if (model=="pair_logit") sum(w[cell]*S[cell])/sum(w[cell]) else
        prod(colSums(w[cell]*f$R[cell,1:2,drop=FALSE])/sum(w[cell]))
    }
  }
  h <- numeric(N)
  if (route=="dr_baseline") h[risk] <- w[risk]*S[risk]/pmax(r,1e-6)[k[risk]] else
    h[risk] <- w[risk]*S[risk]*a[k[risk]]/pmax(q[risk],1e-6)
  response <- y
  if (route %in% c("imputation","augmented")) {
    for (j in seq_along(cohorts)) {
      on <- which(k %in% j); obs <- on[S[on]==1]
      b <- d93r_wls(O[obs,,drop=FALSE],y[obs],w[obs])
      response[on] <- as.vector(O[on,,drop=FALSE]%*%b)
    }
    if (route=="augmented") response <- response+S*(y-response)/pmax(q,1e-6)
    h[] <- 0; h[risk] <- w[risk]*a[k[risk]]
  }
  control <- which(risk & D==0)
  or <- d93r_wls(O[control,,drop=FALSE], response[control], h[control])
  m <- as.vector(O%*%or)
  if (route=="or") {
    dp <- sum(h*D)/N
    T <- sum(h*D*(response-m))/sum(h*D)
    return(list(value=c(mu=T,T=T,dp=dp,or0=or[1],orz=or[2]),q=q,h=h))
  }
  beta <- numeric(2); p <- odds <- numeric(N)
  for (v in 0:1) {
    cell <- risk & f$x==v
    mt <- sum(h[cell]*D[cell]); mc <- sum(h[cell]*(1-D[cell]))
    beta[v+1] <- log(mt/mc)
    p[cell] <- mt/(mt+mc)
    odds[cell] <- if (p[cell][1] < 1-1e-6) mt/mc else (1-1e-6)/(1-(1-1e-6))
  }
  beta <- c(beta[1],beta[2]-beta[1])
  wp <- h*D; wm <- h*(1-D)*odds
  res <- if (route %in% c("dr_baseline","dr_selection","augmented")) response-m else response
  dp <- sum(wp)/N; dm <- sum(wm)/N
  T <- sum(wp*res)/sum(wp); M <- sum(wm*res)/sum(wm)
  # Imputation and augmented details expose unresidualized T,C,Bg,Bo.
  TT <- sum(wp*response)/sum(wp); CC <- sum(wm*response)/sum(wm)
  Bg <- sum(wp*m)/sum(wp); Bo <- sum(wm*m)/sum(wm)
  list(value=c(mu=T-M,T=T,M=M,dp=dp,dm=dm,b0=beta[1],bx=beta[2],
    or0=or[1],orz=or[2],TT=TT,CC=CC,Bg=Bg,Bo=Bo),
    q=q,p=p,h=h,wp=wp,wm=wm,y=y,res=res,rate=r,k=k,S=S,D=D,V=V,
    response=response,m=m,act=as.double(p < 1-1e-6))
}

d93r_fd <- function(f, route, model="pair_logit") {
  h <- 1e-5; base <- d93r_ref(f,route,model=model)$value
  ans <- t(vapply(seq_len(f$N),function(i) {
    cp <- cm <- rep(1,f$N); cp[i] <- 1+h; cm[i] <- 1-h
    f$N*(d93r_ref(f,route,cp,model)$value-d93r_ref(f,route,cm,model)$value)/(2*h)
  },base))
  colnames(ans) <- names(base); ans
}

d93r_owner <- function(f, links=d93r_links(f), outcome_only=FALSE) {
  ctrl <- d93r_fn("ipw_owner_control")(max_chunk_rows=7L,max_retained_entries=1e7,
    max_matrix_entries=2e7,max_work_ops=1e10,max_output_entries=1e7)
  long <- data.frame(id=rep(seq_len(f$N),each=3L),time=rep(1:3,f$N),
    y=as.vector(t(f$Y)),g=rep(f$g,each=3L))
  panel <- d93r_fn("construct_ipw_panel_owner")(long,"id","time","y","g",ctrl)
  X <- d93r_X(f)
  # OR requires exactly the common outcome columns in frame Z_a. Selection
  # still uses its distinct (Intercept, x) design in the reader below.
  frame_X <- if(outcome_only) X[,c(1,3),drop=FALSE] else X
  frame <- d93r_fn("construct_ipw_frame_owner")(panel,seq_len(f$N),frame_X,
    weights=f$w,control=ctrl,frame_g=f$g,allow_frame_only=TRUE)
  reader <- d93r_fn("construct_cp04_frame_reader")(frame,seq_len(f$N),X[,1:2,drop=FALSE],
    f$g,f$w,ctrl,has_intercept=TRUE)
  views <- lapply(links,function(l) d93r_fn("construct_ipw_link_view")(frame,l$s,l$t,
    g=l$g,K=l$K,control=ctrl,link_id=l$link_id,functional_key=paste0("cmpcat:",l$link_id)))
  recipe <- d93r_fn("or_pooled_outcome_recipe")(c("(Intercept)","z"),TRUE,"z")
  list(ctrl=ctrl,frame=frame,reader=reader,views=views,recipe=recipe,X=X)
}

d93r_gram <- function(phi, want, f) {
  a <- sweep(as.matrix(phi),2,colMeans(as.matrix(phi)))
  b <- sweep(as.matrix(want),2,colMeans(as.matrix(want)))
  d93r_eq(crossprod(a)/f$N^2,crossprod(b)/f$N^2)
  cl <- rep(1:4,length.out=f$N)
  d93r_eq(crossprod(rowsum(a,cl))/f$N^2,crossprod(rowsum(b,cl))/f$N^2)
}

testthat::test_that("weighted OR forcing/pullback has an independent clipped WLS refit", {
  f <- d93r_atoms(); o <- d93r_ref(f,"or"); fd <- d93r_fd(f,"or")
  own <- d93r_owner(f,outcome_only=TRUE)
  testthat::expect_true(any(o$q < 1e-6) && any(o$q > 1e-6))
  b <- d93r_fn("fit_or_corollary1")(own$reader,own$views,own$ctrl,d93r_gc,1e-10,own$recipe)
  phi <- numeric(f$N)
  d93r_fn("replay_or_corollary1_if")(b,7L,1L,function(z) phi[z$rows] <<- z$phi[,1])
  d93r_eq(b$links[[1]]$mu,o$value["mu"],1e-8)
  d93r_eq(phi,fd[,"mu"]); d93r_gram(phi,fd[,"mu"],f)
  d93r_eq(b$shared$beta,o$value["orz"],1e-8)
  # Table helper is a second implementation under test, never an expectation.
  tab <- d93r_fn("fit_or_selection_table")(own$X[,1:2],own$X[,c(1,3)],f$Y,f$R,f$g,f$w,
    calendar=1:3,links=d93r_links(f),outcome_recipe=own$recipe,has_intercept=TRUE)
  d93r_eq(tab$phi[,1],fd[,"mu"])
  d93r_eq(tab$beta_if[,1],fd[,"orz"])
  # A frozen observation-weight regression omits the selection pullback.
  # Its case derivative differs: independently recompute WLS with fixed h/cw.
  h0 <- o$h; h <- 1e-5; O <- own$X[,c(1,3)]; y <- numeric(f$N)
  S <- f$R[,1]*f$R[,2]; y[S==1] <- f$Y[S==1,2]-f$Y[S==1,1]
  frozen <- vapply(seq_len(f$N),function(i) {
    calc <- function(cw) {
      hw <- h0*cw; on <- which(f$g==Inf)
      bt <- d93r_wls(O[on,,drop=FALSE],y[on],hw[on])
      sum(hw*(f$g==3)*(y-as.vector(O%*%bt)))/sum(hw*(f$g==3))
    }
    cp <- cm <- rep(1,f$N); cp[i]<-1+h; cm[i]<-1-h
    f$N*(calc(cp)-calc(cm))/(2*h)
  },numeric(1))
  testthat::expect_gt(max(abs(frozen-fd[,"mu"])),1e-6)
})

# Independent DR partials at the raw saturated treatment root; the rate
# derivative uses a strict mask, but score/H and raw rate IF remain unmasked.
d93r_partials <- function(f,o) {
  N <- f$N; dp <- o$value["dp"]; dm <- o$value["dm"]
  score <- (o$h*(o$D-o$p))*o$V
  H <- crossprod(o$V,(o$h*o$p*(1-o$p))*o$V)/N
  rr <- cbind(o$wp*(o$res-o$value["T"])/dp,o$wm*(o$res-o$value["M"])/dm,o$wp,o$wm)
  C <- matrix(0,2,length(o$rate)); U <- matrix(0,4,length(o$rate))
  for (k in seq_along(o$rate)) {
    on <- which(o$k %in% k); scale <- if(o$rate[k]>1e-6) -1/(N*o$rate[k]) else 0
    C[,k] <- colSums(score[on,,drop=FALSE])*scale
    U[,k] <- colSums(rr[on,,drop=FALSE])*scale
  }
  V <- matrix(0,4,2); V[2,] <- colSums((o$wm*o$act*(o$res-o$value["M"]))*o$V)/(N*dm)
  V[4,] <- colSums((o$wm*o$act)*o$V)/N
  list(H=H,C=C,U=U,V=V,score=score,rr=rr)
}

testthat::test_that("DR rate/propensity partials and complete IF agree with independent refits", {
  for (route in c("dr_baseline","dr_selection")) {
    f <- if(route=="dr_baseline") d93r_rates() else d93r_atoms()
    o <- d93r_ref(f,route); fd <- d93r_fd(f,route); own <- d93r_owner(f)
    args <- list(reader=own$reader,link_views=own$views,control=own$ctrl,glm_control=d93r_gc,
      rank_tol=1e-10,outcome_recipe=own$recipe,treatment_columns=c("(Intercept)","x"))
    b <- do.call(d93r_fn(if(route=="dr_baseline") "fit_dr_baseline" else "fit_dr_corollary1"),args)
    ev <- if(route=="dr_baseline") "evaluate_dr_baseline_if_chunk" else "evaluate_dr_corollary1_if_chunk"
    phi <- numeric(f$N); norm <- matrix(0,f$N,2); beta <- matrix(0,f$N,2)
    for (start in seq.int(1L,f$N,by=7L)) {
      rows <- seq.int(start,min(f$N,start+6L))
      z <- d93r_fn(ev)(b,rows,1L,details=TRUE); tr <- z$details$models[[1]]
      phi[rows] <- z$phi[,1]; norm[rows,] <- tr$coordinates[,3:4]; beta[rows,] <- tr$beta
      d93r_eq(d93r_fn(ev)(b,rows,1L,scale="core")$phi[,1],sum(rowSums(f$R)>=2)/f$N*z$phi[,1])
    }
    fit <- b$treatments[[1]]
    d93r_eq(fit$means[c("T","M","dp","dm")],o$value[c("T","M","dp","dm")],1e-8)
    d93r_eq(phi,fd[,"mu"]); d93r_eq(beta,fd[,c("b0","bx")])
    d93r_eq(norm,sweep(fd[,c("dp","dm")],2,o$value[c("dp","dm")],"-"))
    H <- crossprod(o$V,(o$h*o$p*(1-o$p))*o$V)/f$N
    d93r_eq(crossprod(fit$H_chol),H,1e-10)
    testthat::expect_lt(max(abs(colSums((o$h*(o$D-o$p))*o$V))),1e-8)
    d93r_gram(phi,fd[,"mu"],f)
    if(route=="dr_baseline") {
      pp <- d93r_partials(f,o)
      d93r_eq(fit$cross,pp$C,1e-10); d93r_eq(fit$partials$upstream,pp$U,1e-10)
      d93r_eq(fit$partials$beta,pp$V,1e-10)
      testthat::expect_lt(o$rate[2],1e-6); testthat::expect_true(all(o$p[o$h>0] < 1-1e-6))
      testthat::expect_identical(unname(fit$cross[,2]),c(0,0))
    } else {
      testthat::expect_true(any(o$p > 1-1e-6) && any(o$p[o$h>0] < 1-1e-6))
      # The unscaled score mass is only 1.1e-7 in this rare cell. Measure the
      # omitted raw-score forcing in coefficient units using independent H.
      clipped_score <- (o$h*(o$D-pmin(o$p,1-1e-6)))*o$V
      testthat::expect_gt(max(abs(colMeans(clipped_score)%*%solve(H))),1e-6)
    }
  }
})

testthat::test_that("imputation/augmented wdot, final inverse response and complete rows match finite refits", {
  for (route in c("imputation","augmented")) for(bound in c(TRUE,FALSE)) {
    f <- d93r_atoms(bound); model <- if(route=="augmented") "independent_period_logit" else "pair_logit"
    o <- d93r_ref(f,route,model=model); fd <- d93r_fd(f,route,model); own <- d93r_owner(f)
    args <- list(reader=own$reader,link_views=own$views,control=own$ctrl,glm_control=d93r_gc,
      rank_tol=1e-10,outcome_recipe=own$recipe,treatment_columns=c("(Intercept)","x"),
      X_frame=own$X,est_method=if(route=="imputation") "ipw" else "dr")
    if(route=="augmented") args <- c(args,list(selection_model=model,selection_pooling="separate"))
    b <- do.call(d93r_fn(paste0("fit_",route,"_backend")),args)
    lf <- b$links[[1]]; fit <- b$treatments[[1]]
    d93r_eq(c(lf$mu,lf$d_plus,lf$d_minus),o$value[c("mu","dp","dm")],1e-8)
    wantdot <- colSums((o$wm*o$act)*o$V)/f$N
    d93r_eq(lf$denminus_gamma,wantdot,1e-10)
    if(bound) testthat::expect_gt(max(abs(colSums(o$wm*o$V)/f$N-wantdot)),1e-6)
    testthat::expect_identical(any(o$p > 1-1e-6),bound)
    H <- crossprod(o$V,(o$h*o$p*(1-o$p))*o$V)/f$N
    d93r_eq(crossprod(fit$H_chol),H,1e-10)
    phi <- numeric(f$N); norm <- matrix(0,f$N,2); coords <- matrix(0,f$N,4)
    for(start in seq.int(1L,f$N,by=7L)) {
      rows <- seq.int(start,min(f$N,start+6L))
      z <- d93r_fn(paste0("evaluate_",route,"_if_chunk"))(b,rows,1L,details=TRUE)
      ld <- z$details$links[[1]]; phi[rows] <- z$phi[,1]
      norm[rows,] <- cbind(ld$denplus,ld$denminus)
      coords[rows,] <- cbind(ld$T,ld$C,ld$Bg,ld$Bo)
    }
    d93r_eq(phi,fd[,"mu"]); d93r_eq(norm,fd[,c("dp","dm")])
    if(route=="augmented") d93r_eq(coords,fd[,c("TT","CC","Bg","Bo")]) else
      d93r_eq(coords[,1:2],fd[,c("TT","CC")])
    d93r_gram(cbind(phi,norm),cbind(fd[,"mu"],fd[,c("dp","dm")]),f)
    if(route=="augmented") {
      testthat::expect_true(any(o$q < 1e-6) && any(o$q > 1e-6))
      # Confirm the actual cell reconstruction uses the FINAL product floor.
      sel <- b$selection_backend
      for(cid in seq_along(b$cells)) {
        cf <- b$cells[[cid]]; gf <- b$groups[[cf$group_id]]
        rows <- which(f$g==c(3,Inf)[cf$cohort_index])
        pr <- function(rr,pair,read) list(S=f$R[rr,1]*f$R[rr,2],dy=o$y[rr])
        z <- d93r_fn("augmented_cell_rows")(sel,b$up_map,own$X,pr,cf,gf,rows)
        d93r_eq(z$Z,o$response[rows],1e-8)
      }
    }
  }
})

# Equal-weight, complete-support fixture with a KNOWN raw logit root (0,1).
# Choose the design by its two score identities, independently of a fitter.
# The single high-x control lies on a fixed side of the cap; the other cells
# have ordinary probabilities. This keeps native eligibility unchanged.
d93r_complete <- function(bound) {
  po <- 1 - (if(bound) 5e-7 else 2e-6)
  xo <- stats::qlogis(po)
  equation <- function(pc) {
    q <- 100*pc+po; pt <- 1-q/12
    -100*pc*stats::qlogis(pc)+q*stats::qlogis(pt)-po*xo
  }
  pc <- stats::uniroot(equation,c(.02,.05),tol=.Machine$double.eps)$root
  pt <- 1-(100*pc+po)/12
  x <- c(rep(stats::qlogis(pc),100),rep(stats::qlogis(pt),12),xo)
  z <- c(rep(c(0,1),50),rep(c(0,1),6),.5)
  g <- c(rep(Inf,100),rep(3,12),Inf)
  dy <- 1+.4*z+.2*sin(seq_along(x))+.3*(g==3)
  # Both adjacent changes have the same slope and residuals; their joint OR
  # score is exactly the one-pair WLS score for every case-weight perturbation.
  list(N=length(x),x=x,z=z,w=rep(1,length(x)),R=matrix(1,length(x),3),
    Y=cbind(0,dy,2*dy+.25),g=g,K=Inf)
}

d93r_complete_ref <- function(f) {
  N <- f$N; X <- d93r_X(f); V <- X[,1:2]; O <- X[,c(1,3)]
  D <- as.double(f$g==3); y <- f$Y[,2]-f$Y[,1]
  p <- stats::plogis(f$x); pc <- pmin(p,1-1e-6); odds <- pc/(1-pc)
  score <- (D-p)*V; H <- crossprod(V,(p*(1-p))*V)/N
  xb <- score%*%solve(H)
  on <- which(D==0); bo <- d93r_wls(O[on,,drop=FALSE],y[on],rep(1,length(on)))
  m <- as.vector(O%*%bo); es <- y-m
  os <- (1-D)*es*O; Ho <- crossprod(O,(1-D)*O)/N; xbo <- os%*%solve(Ho)
  wp <- D; wm <- (1-D)*odds; dp <- mean(wp); dm <- mean(wm)
  T <- sum(wp*es)/sum(wp); M <- sum(wm*es)/sum(wm)
  Vb <- colSums((wm*(p<1-1e-6)*(es-M))*V)/(N*dm)
  Lo <- -colSums((wp/dp-wm/dm)*O)/N
  phi <- wp*(es-T)/dp-wm*(es-M)/dm-as.vector(xb%*%Vb)+as.vector(xbo%*%Lo)
  norm <- cbind(wp-dp,wm-dm+as.vector(xb%*%(colSums((wm*(p<1-1e-6))*V)/N)))
  list(phi=phi,norm=norm,H=H,score=score,xb=xb,mu=T-M,p=p)
}

testthat::test_that("owned native DR parity and capped-control summary decline are separate cases", {
  # At this extreme but finite root, deviance can stagnate at roughly 1e-13
  # under R 4.5.0. Keep the independent score/IF checks strict while allowing
  # glm.fit to declare convergence before rounding noise exhausts maxit.
  fixture_gc <- utils::modifyList(d93r_gc, list(epsilon = 1e-12))
  for(bound in c(FALSE,TRUE)) {
    f <- d93r_complete(bound); ref <- d93r_complete_ref(f)
    links <- c(d93r_links(f),list(list(s=2,t=3,g=3,K=Inf,link_id="D93second")))
    own <- d93r_owner(f,links)
    testthat::expect_lt(max(abs(colSums(ref$score))),1e-10)
    testthat::expect_identical(any(ref$p[f$g==Inf] > 1-1e-6),bound)
    b <- d93r_fn("fit_dr_baseline")(own$reader,own$views,own$ctrl,fixture_gc,1e-10,
      own$recipe,treatment_columns=c("(Intercept)","x"))
    testthat::expect_false(is.null(b$complete_support))
    count <- b$accounting$complete_support$n_treatment_summary_batched_models
    if(bound) d93r_eq(count,0,0) else testthat::expect_gt(count,0)
    testthat::expect_true(all(is.finite(c(b$accounting$retained_bound,
      b$accounting$fit_scratch_bound,b$accounting$work_ops_bound))))
    testthat::expect_lte(b$accounting$retained_bound,own$ctrl$max_retained_entries)
    testthat::expect_lte(b$accounting$fit_scratch_bound,own$ctrl$max_matrix_entries)
    testthat::expect_lte(b$accounting$work_ops_bound,own$ctrl$max_work_ops)
    if(bound) {
      d93r_eq(b$accounting$complete_support$n_treatment_logit_fits,1,0)
      d93r_eq(b$accounting$complete_support$n_treatment_logit_reuses,1,0)
      # A finite cap just below the completed reservation must reach and refuse
      # the structural-decline reservation, before the remaining model runs.
      tight <- own$ctrl; tight$max_work_ops <- b$accounting$work_ops_bound-1
      err <- tryCatch(d93r_fn("fit_dr_baseline")(own$reader,own$views,tight,fixture_gc,
        1e-10,own$recipe,treatment_columns=c("(Intercept)","x")),error=identity)
      testthat::expect_s3_class(err,"cdid_ipw_owner_bound_error")
      testthat::expect_match(conditionMessage(err),
        "structural-decline original-loop work|original-loop fit plus original-factor validation")
    }
    phi <- numeric(f$N)
    for(start in seq.int(1L,f$N,by=7L)) {
      rows <- seq.int(start,min(f$N,start+6L))
      native <- d93r_fn("evaluate_dr_baseline_if_chunk")(b,rows,1:2,details=FALSE)
      R <- d93r_fn("evaluate_dr_baseline_if_chunk")(b,rows,1:2,details=TRUE)
      d93r_eq(native$phi,R$phi,1e-8)
      d93r_eq(native$phi[,1],ref$phi[rows],1e-8)
      d93r_eq(R$details$models[[1]]$coordinates[,3:4,drop=FALSE],ref$norm[rows,,drop=FALSE],1e-8)
      d93r_eq(R$details$models[[1]]$beta,ref$xb[rows,,drop=FALSE],1e-8)
      phi[rows] <- native$phi[,1]
    }
    d93r_eq(b$links[[1]]$mu,ref$mu,1e-8); d93r_gram(phi,ref$phi,f)
    d93r_eq(crossprod(b$treatments[[1]]$H_chol),ref$H,1e-10)
  }
})

testthat::test_that("native treatment-summary rate floor is actually active with uncapped propensities", {
  # Owned complete-support eligibility implies rates==1. To exercise the C
  # rate mask itself, supply the actual observer rows/weights of the independent
  # incomplete-rate fixture at the kernel boundary. This is kernel parity,
  # not a claim that an incomplete roster qualifies for the owned native route.
  f <- d93r_rates(); o <- d93r_ref(f,"dr_baseline"); pp <- d93r_partials(f,o)
  testthat::expect_lt(o$rate[2],1e-6)
  testthat::expect_true(all(o$p[o$h>0] < 1-1e-6))
  rows <- as.integer(which(o$h>0)); X <- d93r_X(f); N <- f$N
  co <- lapply(seq_along(o$rate),function(k) as.integer(which(o$k[rows] %in% k)))
  ss <- vapply(co,function(rr) colSums(pp$score[rows[rr],,drop=FALSE]),numeric(2))
  or <- o$value[c("or0","orz")]
  res <- .Call(d93r_fn("or_baseline_native_kernel")("C_gmm_cp08_dr_treatment_batch"),
    as.double(c(rep(0,N),o$y,rep(0,N))),as.integer(c(1,N+1,2*N+1)),
    rep(seq_len(N),3L),seq_len(3L*N),rows,seq_len(N),o$V[rows,,drop=FALSE],
    X[rows,3,drop=FALSE],as.double(or[2]),o$wp[rows],o$wm[rows],co,ss,
    colSums(o$wm*o$V),as.double(c(sum(o$wp),sum(o$wm),o$value[c("dp","dm")],N)),
    matrix(c(1L,2L),2L),matrix(c(1,or[1],o$rate),2+length(o$rate)),
    as.integer(c(7,isTRUE(unname(capabilities("long.double"))))),2e7)
  testthat::expect_identical(res$status[1],0L)
  d93r_eq(res$means[1:4,1],o$value[c("T","M","dp","dm")],1e-8)
  d93r_eq(matrix(res$cross[,1],2),pp$C,1e-10)
  d93r_eq(matrix(res$U_dr[,1],4),pp$U,1e-10)
  d93r_eq(matrix(res$V_dr[,1],4),pp$V,1e-10)
  testthat::expect_identical(unname(matrix(res$U_dr[,1],4)[,2]),rep(0,4))
  # Omitting this mask is observable on the same supplied native inputs.
  unmasked <- -colSums(pp$rr[o$k %in% 2,,drop=FALSE])/(N*o$rate[2])
  testthat::expect_gt(max(abs(unmasked)),1e-6)
})

testthat::test_that("finite clipping never rescues genuine raw-root and support refusals", {
  f <- d93r_atoms(); own <- d93r_owner(f)
  # The raw MLE design is rank deficient even though bounded weights exist.
  X <- own$X[,c(1,1),drop=FALSE]; colnames(X) <- c("one","duplicate")
  testthat::expect_error(d93r_fn("fit_dr_baseline_table")(own$X[,c(1,3)],X,f$Y,f$R,f$g,f$w,
    calendar=1:3,links=d93r_links(f),outcome_recipe=own$recipe),class="cdid_cp04_rank_error")
  zero <- f$R; zero[f$g==Inf,1] <- 0
  testthat::expect_error(d93r_fn("fit_dr_baseline_table")(own$X[,c(1,3)],own$X[,1:2],
    f$Y,zero,f$g,f$w,calendar=1:3,links=d93r_links(f),outcome_recipe=own$recipe),
    class="cdid_ipw_fit_error")
  testthat::expect_error(d93r_fn("cp04_logit")(cbind(1,c(-2,-1,1,2)),c(0,0,1,1),
    rep(1,4),d93r_gc,4,list()),class="cdid_cp04_separation_error")
  testthat::expect_error(d93r_fn("cp04_logit")(cbind(1,c(-2,-1,0,0,1,2)),c(0,1,0,1,1,1),
    rep(1,6),list(epsilon=1e-13,maxit=1L,trace=FALSE),6,list()),class="cdid_cp04_convergence_error")
  bad <- own$X[,1:2]; bad[1,2] <- Inf
  testthat::expect_error(d93r_fn("cp04_logit")(bad,as.double(f$g==3),f$w,
    d93r_gc,f$N,list()),class="cdid_ipw_fit_error")
})

testthat::test_that("public att_gt and conditional pretest activate clipped paths", {
  f <- d93r_atoms(); N <- f$N
  panel <- data.frame(id=rep(seq_len(N),each=3L),time=rep(1:3,N),y=as.vector(t(f$Y)),
    g=rep(ifelse(is.infinite(f$g),0,f$g),each=3L))
  frame <- data.frame(id=seq_len(N),g=ifelse(is.infinite(f$g),0,f$g),x=f$x,w=f$w)
  a <- list(yname="y",tname="time",idname="id",gname="g",data=panel,frame=frame,
    xformla=~x,weightsname="w",control_group="nevertreated",est_method="dr",mode="chained",
    selection="conditional_mtar",selection_model="pair_logit",selection_pooling="separate",
    selection_xformla=~x,bstrap=FALSE,cband=FALSE)
  fit <- do.call(d93r_fn("att_gt"),c(a,list(gmm_weighting="identity")))
  testthat::expect_s3_class(fit,"MP")
  tr <- fit$chained$nuisance_diagnostics$records$treatment
  testthat::expect_gt(sum(vapply(tr,function(t) t$clipping[["p_capped_ctrl"]],numeric(1))),0)
  testthat::expect_true(all(is.finite(fit$att)))
  # The public conditional API always draws multipliers. Root runs this small
  # activation check; no draw is made during parse/review, no size claim follows.
  pre <- do.call(d93r_fn("conditional_did_pretest"),c(a,list(biters=19L)))
  testthat::expect_s3_class(pre,"MP.TEST")
  testthat::expect_identical(pre$chained$status,"available")
  testthat::expect_true(is.finite(pre$CvM))
})
