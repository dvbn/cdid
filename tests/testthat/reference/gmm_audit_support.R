load_core <- function(path='R/chained_gmm_core.R',env=globalenv()) {
  for (p in c('R/forward.R','R/chained_pairs.R','R/chained_design.R',path)) sys.source(p,env)
}
audit_fixture <- function() {
  links <- data.frame(g=c(2,2,2,3),s=c(1,1,2,4),t=c(2,3,3,5),ref=rep(1,4),
    functional_key=c('fresh_A','fresh_A','fresh_A','fresh_B'),s_qualified=rep(TRUE,4),t_qualified=c(TRUE,FALSE,FALSE,TRUE))
  design <- construct_chained_design(links,target_post_times=data.frame(g=2,time=4))
  B <- matrix(c(2,0,0,0,1,3,0,0,-1,1,2,0,1,-1,1,2),4,4,byrow=TRUE)
  V <- tcrossprod(B)+diag(4)
  targets <- list(gmm_target_spec('nuis',1,1),gmm_target_spec('level',2,1),
    gmm_target_spec('difference',c(1,2),c(-1,1)),gmm_target_spec('unanch_contrast',c(4,5),c(-1,1)))
  Phi <- t(vapply(0:7,function(i)c(i-3,(i*i)%%7-2,(2*i)%%5-1,(-1)^i/2),numeric(4)))
  list(design=design,V=V,mhat=c(7/4,19/8,-5/2,13/4),targets=targets,Phi=Phi)
}
audit_bounds <- function() gmm_bounds_control(50000,8,16,1,TRUE,1000)
# D12 control plumbing: allow callers to pass optional max_loading_error_abs (default NULL = not supplied)
# without selecting any hidden default or modifying existing control values.
audit_accuracy <- function(weighting='optimal', max_loading_error_abs=NULL) gmm_accuracy_control(
  if(weighting=='optimal') 'certified_eigen_chol_kkt' else 'certified_qr_kkt','abs',1e-10,NA,
  1e-11,if(weighting=='optimal') 1e-8 else NA,1e-12,if(weighting=='optimal') 1e-6 else NA,
  1e-6,1e-6,1e-12,1e-12,FALSE,max_loading_error_abs=max_loading_error_abs)
audit_replay <- function(Phi,cluster=FALSE,batch_size=3,mutate=NULL) {
  uid <- paste0('id',0:7); wi <- rep(1:4,each=2); cid <- rep(paste0('c',0:3),each=2)
  pass <- 0L
  cb <- function(consume) {
    pass <<- pass+1L; P <- Phi
    if (!is.null(mutate)) P <- mutate(P,pass)
    for (st in seq(1,8,by=batch_size)) {
      ix <- st:min(8,st+batch_size-1)
      consume(list(unit_id=uid[ix],cluster_id=cid[ix],phi=P[ix,,drop=FALSE],weight=wi[ix],nuisance_fit_id='frozen_assumed'))
    }
  }
  gmm_if_replay_spec(cb,'unit',cluster,8,if(cluster)4 else 8,if(cluster)4 else 8,
    TRUE,TRUE,'case_weights',1/8,'explicit_scale','frozen_assumed',uid,if(cluster)paste0('c',0:3) else NULL)
}
audit_provenance <- function() gmm_covariance_provenance('unit',8,NA,8,TRUE,TRUE,'case_weights',1/8,'criterion_scale','frozen_assumed')
audit_solve <- function(fx,weighting='optimal',gauge='lowest',targets=fx$targets, replay=audit_replay(fx$Phi),bounds=audit_bounds(),max_loading_error_abs=NULL,accuracy=audit_accuracy(weighting,max_loading_error_abs=max_loading_error_abs)) {
  solve_chained_gmm_core(fx$design,fx$mhat,fx$V,audit_provenance(),targets,weighting,replay,
    accuracy,bounds,gauge,request_param_cov=TRUE,materialize_target_if=TRUE)
}
