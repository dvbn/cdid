# Fixture, packaged-reference loaders and scoped spy for the CP02 shared-IPW
# audit tests (test-chained-ipw-shared-audit.R, -independent.R). Ported from
# an independent audit's loader: every production call goes to the loaded
# cdid namespace; no production R file is sourced. Defines functions only.

audit_fixture <- function(weighted=TRUE) {
  n <- 96L; G <- rep(c(4,7,Inf,9), each=24L)
  x <- rep(seq(-1.4,1.4,length.out=24L),4L)
  Z <- cbind(one=1,x=x)
  w <- if(weighted) 0.7 + ((seq_len(n)*13L)%%29L)/20 else rep(1,n)
  obs <- matrix(FALSE,n,8L)
  for(i in seq_len(n)) {
    k <- ((i-1L)%%24L)+1L
    if(G[i] != 9 && k<=18L) obs[i,1:6] <- TRUE else {
      obs[i,7:8] <- TRUE
      obs[i,1L+(k%%6L)] <- TRUE
    }
  }
  Y <- outer(seq_len(n)/17,1:8,function(a,b) sin(a*b)+a+b^2/13) + outer(x,1:8,function(a,b) a*b/3)
  pan <- data.frame(id=rep(seq_len(n),each=8L),time=rep(as.double(1:8),n),y=as.vector(t(ifelse(obs,Y,NA_real_))),g=rep(G,each=8L))
  idx <- cdid:::construct_observed_pairs(pan,'id','time','y','g')
  ctx <- cdid:::construct_ipw_context(idx,id=seq_len(n),Z=Z,weights=w)
  list(n=n,G=G,Z=Z,w=w,obs=obs,Y=Y,panel=pan,index=idx,context=ctx)
}
audit_control <- stats::glm.control(epsilon=1e-12,maxit=100,trace=FALSE)
audit_bounds <- list(max_links=100,max_classes=20,max_class_dates=20,max_E_entries=500,max_mask_entries=100000)
audit_pairs <- list(c(1,2),c(2,3),c(1,3))
audit_setup <- function(fx=audit_fixture()) {
  ses <- cdid:::construct_ipw_shared_session(fx$context,4,40,1e5)
  nu <- cdid:::construct_ipw_nuisance(ses,4,c(7,Inf),1,2,audit_control)
  decl <- lapply(audit_pairs,function(p) cdid:::declare_ipw_shared_link(ses,nu,4,c(7,Inf),p[1],p[2],audit_control,functional_key='A'))
  des <- cdid:::build_chained_design(data.frame(g=4,s=c(1,2,1),t=c(2,3,3),ref=3,functional_key='A'))
  list(fx=fx,session=ses,nuisance=nu,declarations=decl,design=des)
}
audit_cert <- function(z,des=z$design,bounds=audit_bounds) cdid:::certify_ipw_structural_classes(z$session,z$declarations,des,1:3,'optimal',bounds)

# Packaged external references, bound by SHA-256 before use. A missing, empty
# or altered reference is an error, never a skip.
audit_reference_sha <- c(
  chained_ipw_link_original_6550dc.R = "6550dc592451cb7e9b7f45573d083fa14cacc1f3a09421a3330d31997183b0b5",
  cp02_shared_ipw_audit_newton_reference.R = "2e5647c6b0530e83f7701aae00684fddbf5f72ac37e6000c71c6a5bfffd94b72"
)
audit_load_reference <- function(file, env) {
  path <- testthat::test_path("reference", file)
  if (!file.exists(path)) stop("missing packaged reference: ", file, call. = FALSE)
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required to bind packaged references.", call. = FALSE)
  if (!identical(digest::digest(file = path, algo = "sha256"), audit_reference_sha[[file]])) {
    stop("packaged reference bytes differ: ", file, call. = FALSE)
  }
  # Only the immutable original kernel receives the one matched-start edit.
  # The independent Newton/IF reference is sourced without adaptation.
  if (file == "chained_ipw_link_original_6550dc.R") {
    ipw_reference_source(path, envir = env)
  } else {
    sys.source(path, envir = env)
  }
  env
}
audit_original_kernel <- function() {
  env <- new.env(parent = baseenv())
  env$cdid_abort <- cdid:::cdid_abort
  audit_load_reference("chained_ipw_link_original_6550dc.R", env)
}
audit_newton_reference <- function() {
  audit_load_reference("cp02_shared_ipw_audit_newton_reference.R", new.env(parent = baseenv()))
}

# Scoped spy: isolated clones of the actual cdid namespace closures, enclosed
# by a child environment holding counting wrappers that delegate to the real
# namespace functions. The namespace is never edited. Binding is proved on
# every construction: each clone has the namespace body and formals, and every
# hooked name is a namespace closure with the signature the wrapper forwards.
audit_spy_hooks <- list(ipw_glm_fit_internal = c("Z_use", "D_use", "fit_weights", "ctrl"),
                        ipw_pair_outcomes = c("context", "hit"),
                        query_pair_mask_internal = NULL,
                        ipw_forest_path = c("a", "b", "forest", "depth", "count_only"))
audit_spy <- function() {
  ns <- asNamespace("cdid")
  e <- new.env(parent = ns)
  cloned <- character(0L)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, envir = ns, inherits = FALSE)
    if (is.function(f) && identical(environment(f), ns)) {
      environment(f) <- e
      assign(nm, f, envir = e)
      cloned <- c(cloned, nm)
    }
  }
  for (nm in cloned) {
    f <- get(nm, envir = ns, inherits = FALSE)
    g <- get(nm, envir = e, inherits = FALSE)
    if (!identical(body(g), body(f)) || !identical(formals(g), formals(f)) || !identical(environment(g), e)) {
      stop("spy clone is not bound to cdid namespace function ", nm, call. = FALSE)
    }
  }
  if (!all(names(audit_spy_hooks) %in% cloned)) stop("spy hook target missing from cdid namespace", call. = FALSE)
  real <- mget(names(audit_spy_hooks), envir = ns, inherits = FALSE)
  for (nm in names(audit_spy_hooks)) {
    if (!is.null(audit_spy_hooks[[nm]]) && !identical(names(formals(real[[nm]])), audit_spy_hooks[[nm]])) {
      stop("spy hook signature differs from cdid namespace function ", nm, call. = FALSE)
    }
  }
  count <- new.env(); count$glm<-0L; count$out<-list(); count$mask<-0L; count$paths<-0L
  e$ipw_glm_fit_internal <- function(...) {count$glm<-count$glm+1L;real$ipw_glm_fit_internal(...)}
  e$ipw_pair_outcomes <- function(context,hit) {count$out[[length(count$out)+1L]]<-if(length(hit$obs_s)) c(context$observations$time[hit$obs_s[1]],context$observations$time[hit$obs_t[1]]) else numeric();real$ipw_pair_outcomes(context,hit)}
  e$query_pair_mask_internal <- function(...) {count$mask<-count$mask+1L;real$query_pair_mask_internal(...)}
  e$ipw_forest_path <- function(a,b,forest,depth,count_only) {if(!count_only)count$paths<-count$paths+1L;real$ipw_forest_path(a,b,forest,depth,count_only)}
  list(e=e,count=count,cloned=cloned,real=real)
}
