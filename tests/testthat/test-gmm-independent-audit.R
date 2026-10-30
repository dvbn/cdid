# Independent finite references and adversarial checks; no statistical certification.
sys.source(file.path('reference','gmm_audit_support.R'), envir=environment())
codex_ref <- read.csv(file.path('reference','codex_independent_kkt.csv'),stringsAsFactors=FALSE)
ref_mat <- function(fixture, field) {
  d <- codex_ref[codex_ref$fixture==fixture & codex_ref$field==field,]
  stopifnot(nrow(d)>0)
  z <- matrix(0,max(d$row),max(d$col)); z[cbind(d$row,d$col)] <- as.numeric(d$numerator)/as.numeric(d$denominator); z
}

test_that('Independent fresh exact joint KKT, loadings, offdiagonals and projected IF covariance', {
  fx <- audit_fixture()
  for (wt in c('optimal','identity')) {
    res <- audit_solve(fx,wt)
    ig <- res$internal_gauge
    W <- materialize_gmm_gauge_matrix(ig$gauge,50000)
    expect_equal(unname(W),ref_mat(wt,'W'),tolerance=0)
    expect_equal(unname(fx$V),ref_mat(wt,'V'),tolerance=0)
    expect_equal(unname(ig$theta_gauge_coords),drop(ref_mat(wt,'theta')),tolerance=1e-13)
    expect_equal(unname(ig$y_multiplier),drop(ref_mat(wt,'y')),tolerance=1e-13)
    expect_equal(unname(ig$Lambda),ref_mat(wt,'Lambda'),tolerance=1e-13)
    expect_equal(unname(ig$Mc),ref_mat(wt,'Mc'),tolerance=1e-13)
    expect_equal(res$targets$estimate,drop(ref_mat(wt,'tau')),tolerance=1e-13)
    expect_equal(unname(res$raw_Z1tZ1),ref_mat(wt,'criterion_cov'),tolerance=1e-13)
    expect_equal(unname(res$sandwich_cov),ref_mat(wt,'sandwich'),tolerance=1e-13)
    expect_equal(unname(res$if_projection$materialized_target_if),ref_mat(wt,'projected_IF'),tolerance=1e-13)
    expect_true(any(abs(res$sandwich_cov[row(res$sandwich_cov)!=col(res$sandwich_cov)])>0.1))
    expect_false(res$model_covariance_provenance_certified)
    expect_null(res$model_cov)
    expect_true(is.na(res$J_stat_optimal))
    other <- audit_solve(fx,wt,'highest')
    expect_equal(other$targets$estimate,res$targets$estimate,tolerance=1e-13)
    expect_equal(unname(other$internal_gauge$Lambda),unname(ig$Lambda),tolerance=1e-13)
    expect_true(all(is.na(res$parameters$level_estimate[!res$parameters$level_estimable])))
    for (bs in c(1,3,8)) {
      rr <- audit_solve(fx,wt,replay=audit_replay(fx$Phi,TRUE,bs))
      expect_equal(unname(rr$sandwich_cov),ref_mat(wt,'cluster_sandwich'),tolerance=1e-13)
    }
    # Correlation removal changes the joint answer, including the disconnected contrast.
    fxb <- fx; fxb$V[4,1:3] <- fxb$V[1:3,4] <- 0
    if (wt=='optimal') expect_gt(abs(audit_solve(fxb,wt)$targets$estimate[4]-res$targets$estimate[4]),0.01)
  }
})

test_that('Independent exact dyadic target predicates and representability refusal', {
  fx <- audit_fixture(); g <- build_gmm_gauge(fx$design,'lowest'); b <- audit_bounds()
  mn <- 2^-1074
  expect_false(exact_binary64_dyadic_sum(c(2^1020,mn,-2^1020))$is_zero)
  expect_identical(exact_binary64_dyadic_sum(c(2^1020,-mn,-2^1020))$sign,-1L)
  for (v in list(c(2^1020,mn,-2^1020,-mn),c(0.1,-0.1),c(-0,0))) expect_true(exact_binary64_dyadic_sum(v)$is_zero)
  # Exact sum of duplicates is 1, despite ordinary long-double summation losing it.
  tt <- list(gmm_target_spec('exact_cancel',c(1,1,1),c(2^100,1,-2^100)))
  ti <- validate_gmm_targets(fx$design,g,tt,b)
  expect_equal(ti$c_vectors[[1]][1],1,tolerance=0)
  expect_equal(audit_solve(fx,targets=tt)$targets$estimate,audit_solve(fx)$targets$estimate[1],tolerance=1e-13)
  tt <- list(gmm_target_spec('inexact_duplicate',c(1,1),c(2^100,0.1)))
  expect_error(validate_gmm_targets(fx$design,g,tt,b),class='cdid_gmm_target_error')
  for (a in c(mn,2^1020)) {
    t <- list(gmm_target_spec('dyadic_contrast',c(4,5),c(-a,a)))
    expect_true(validate_gmm_targets(fx$design,g,t,b)$estimable)
  }
  for (t in list(gmm_target_spec('unanch_level',4,1),gmm_target_spec('isolated',3,1),
    gmm_target_spec('hidden_subnormal',c(4,4,5),c(2^1020,mn,-2^1020))))
    expect_false(validate_gmm_targets(fx$design,g,list(t),b)$estimable)
  expect_error(gmm_target_spec('empty',numeric(0),numeric(0)),class='cdid_gmm_target_error')
  expect_error(gmm_target_spec('overflow_index',2^32,1),class='cdid_gmm_target_error')
  M <- matrix(1,5,2,dimnames=list(NULL,c('same','same')))
  expect_error(validate_gmm_targets(fx$design,g,M,b),class='cdid_gmm_target_error')
})

test_that('Independent replay checksum collision, second-pass validation, counts and overflow refuse', {
  fx <- audit_fixture(); b <- audit_bounds(); p <- audit_provenance()
  # A collision in the old coordinate-weighted checksum is now detected byte-for-byte.
  pass <- 0L
  cb <- function(consume) {
    pass <<- pass+1L; P <- matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE)
    if(pass==2L) P[1,] <- P[1,]+c(2,-1)
    consume(list(unit_id=c('a','b','c'),phi=P))
  }
  r <- gmm_if_replay_spec(cb,'unit',FALSE,3,3,3,TRUE,FALSE,'uniform',1,'s','frozen_assumed')
  expect_error(project_gmm_if_stream(r,diag(2),2,b,p),class='cdid_gmm_replay_error')
  for (mut in list(function(P,pass) {if(pass==2) P[1,1]<-Inf;P},function(P,pass) {if(pass==2) P<-P[1:7,,drop=FALSE];P})) {
    rp <- audit_replay(fx$Phi,mutate=mut)
    expect_error(audit_solve(fx,replay=rp),class='cdid_gmm_replay_error')
  }
  rp <- audit_replay(fx$Phi); rp$expected_unit_ids[1] <- 'wrong'
  expect_error(audit_solve(fx,replay=rp),class='cdid_gmm_replay_error')
  rp <- audit_replay(fx$Phi); rp$n_input_rows <- rp$n_output_rows <- 9
  expect_error(audit_solve(fx,replay=rp),class='cdid_gmm_replay_error')
  rp <- audit_replay(fx$Phi); rp$callback <- function(consume) consume(list(unit_id=rep('dup',8),phi=fx$Phi,weight=rep(1:4,each=2)))
  expect_error(audit_solve(fx,replay=rp),class='cdid_gmm_replay_error')
  rp <- audit_replay(fx$Phi); rp$callback <- function(consume) consume(list(unit_id=paste0('id',0:7),phi=matrix(1e308,8,4),weight=rep(1:4,each=2)))
  expect_error(audit_solve(fx,replay=rp),class='cdid_gmm_replay_error')
  rp <- audit_replay(fx$Phi); rp$callback <- function(consume) NULL
  expect_error(audit_solve(fx,replay=rp),class='cdid_gmm_replay_error')
})

test_that('Independent preflight, graph metadata and conditional class labels cannot certify', {
  fx <- audit_fixture(); b <- audit_bounds(); b$max_matrix_entries <- 10
  called <- FALSE; rp <- audit_replay(fx$Phi); rp$callback <- function(consume) {called <<- TRUE;stop('called')}
  expect_error(audit_solve(fx,replay=rp,bounds=b),class='cdid_gmm_bound_error'); expect_false(called)
  b <- audit_bounds(); rp$n_input_rows <- rp$n_output_rows <- 100000
  expect_error(audit_solve(fx,replay=rp,bounds=b),class='cdid_gmm_bound_error'); expect_false(called)
  fx2 <- fx; fx2$design$parameters$level_estimable[5] <- TRUE
  expect_error(audit_solve(fx2),class='cdid_gmm_design_error')
  fx2 <- fx; dimnames(fx2$V) <- list(paste0('wrong',1:4),NULL)
  expect_error(audit_solve(fx2),class='cdid_gmm_covariance_error')
  expect_error(validate_gmm_reduction_spec(list(mode='pre_reduced',certified_by='made_up'), 'optimal'),class='cdid_gmm_reduction_error')
  ca <- data.frame(s=c(1,1,2,10,10,11),t=c(2,3,3,11,12,12),context_id=rep('a',6),g=20,
    operator_id='b',vertex_map_id='v',pair_s_id='p',nuisance_object_id='n',positive_rates=TRUE,regular_logit=TRUE)
  expect_error(build_certified_class_forest(ca,'optimal'),class='cdid_gmm_bound_error')
  f <- build_certified_class_forest(ca,'optimal',50000)
  expect_false(f$structural_proof_verified); expect_identical(f$retained_links,c(1L,2L,4L,5L))
  expect_equal(unname(f$E[3,]),c(-1L,1L,0L,0L),tolerance=0)
  ca$context_id[1:3] <- 'a|20'; ca$operator_id[1:3] <- 'b'
  ca$context_id[4:6] <- 'a'; ca$operator_id[4:6] <- '20|b'
  expect_equal(length(unique(build_certified_class_forest(ca,'optimal',50000)$class_key)),2)
})

test_that('Independent outward operations enclose subnormal/huge norms and refuse bad candidates', {
  for (x in c(2^-1074,2^-600,1,2^500,2^1000)) {
    n <- gmm_iv_norm(x)
    expect_lte(n[1],x); expect_gte(n[2],x)
  }
  expect_identical(gmm_next(0),2^-1074)
  expect_identical(gmm_next(1,FALSE),1-2^-53)
  expect_identical(gmm_next(1),1+2^-52)
  expect_error(gmm_iv_mul(gmm_iv(1e308),gmm_iv(1e308)),class='cdid_gmm_accuracy_error')
  fx <- audit_fixture(); res <- audit_solve(fx); ig <- res$internal_gauge; W <- materialize_gmm_gauge_matrix(ig$gauge,50000)
  ti <- validate_gmm_targets(fx$design,ig$gauge,fx$targets,audit_bounds())
  a <- list(V_R=fx$V,R_chol=ig$R_chol,W_G=W,A_whitened=backsolve(ig$R_chol,W,transpose=TRUE),
    U_qr=ig$U_qr,piv=ig$pivot,mhat_R=fx$mhat,theta=ig$theta_gauge_coords,y=ig$y_multiplier,
    C_gauge=ti$C_gauge,Lambda=ig$Lambda,Mc=ig$Mc,Z1=ig$Z1,weighting='optimal',accuracy=audit_accuracy())
  for (kind in c('pivot','transpose','sign','crosscov')) {
    bad <- a
    if (kind=='pivot') bad$theta[bad$piv] <- rev(bad$theta[bad$piv])
    if (kind=='transpose') bad$A_whitened <- backsolve(ig$R_chol,W)
    if (kind=='sign') bad$y <- -bad$y
    if (kind=='crosscov') bad$V_R[1,4] <- bad$V_R[4,1] <- 0
    expect_error(do.call(certify_gmm_numerical_guard,bad),class='cdid_gmm_refusal')
  }
})

test_that('Internal recovery distinguishes criterion algebra from unverified model provenance', {
  fx <- audit_fixture(); fx$Phi[,] <- 0
  p <- audit_provenance(); p$scale_id <- 'explicit_scale'
  a <- audit_accuracy(); a$require_matching_scale <- TRUE
  res <- solve_chained_gmm_core(fx$design,fx$mhat,fx$V,p,fx$targets,'optimal',
    audit_replay(fx$Phi),a,audit_bounds(),'lowest')
  expect_true(res$declared_scale_matches_inference)
  expect_false(res$scale_matches_inference)
  expect_false(res$model_covariance_provenance_certified)
  expect_null(res$model_cov)
  expect_true(all(is.na(res$targets$se_model)))
  expect_true(is.na(res$J_stat_optimal))
  expect_equal(unname(res$sandwich_cov),matrix(0,4,4),tolerance=0)
  expect_gt(max(abs(res$raw_Z1tZ1)),0.1)
  expect_true(all(res$targets$se_criterion_conditional>0))
  expect_gt(res$raw_quadratic_residual,0)
  fx <- audit_fixture(); fx$mhat <- fx$mhat*2^520
  fx$targets <- list(gmm_target_spec('zero_target',1,0))
  expect_error(audit_solve(fx,'identity'),class='cdid_gmm_accuracy_error')
})

test_that('Native C backend registration, symbol forcing, arithmetic probes and no-fallback refusal', {
  nb <- gmm_native_backend_info()
  expect_identical(nb$backend, 'native_c_ieee754_binary64')
  expect_true(isTRUE(nb$c_registration))
  expect_true(isTRUE(nb$force_symbols))
  expect_identical(nb$fallback_permitted, FALSE)
  expect_identical(nb$flt_radix, 2L)
  expect_identical(nb$dbl_mant_dig, 53L)
  expect_identical(nb$rounding_mode, 'FE_TONEAREST')
  expect_true(isTRUE(nb$gradual_underflow_verified))
  expect_true(isTRUE(nb$volatile_fp_stores))
  expect_true(isTRUE(nb$dll_loaded))
  # October 16, 2026: the measured post-fit interval-product kernel
  # (`C_gmm_iv_mat_product`, 5 arguments) is the fifth registered hotspot.
  expect_identical(
    nb$compiled_hotspots,
    c('C_gmm_iv_dot', 'C_gmm_iv_norm', 'C_gmm_iv_mat_residual', 'C_gmm_iv_eigen_residual',
      'C_gmm_iv_mat_product')
  )

  dll <- base::getLoadedDLLs()[['cdid']]
  expect_false(is.null(dll))
  reg <- base::getDLLRegisteredRoutines(dll, addNames = TRUE)[['.Call']]
  # Exact current registration table in src/cdid.c (21 routines). The generic
  # GMM resolver still authenticates its original nine plus each route's arity.
  expected_registry <- c(
    'C_gmm_native_status' = 0L,
    'C_gmm_test_arithmetic_probe' = 5L,
    'C_gmm_sexp_aliases' = 2L,
    'C_gmm_next' = 2L,
    'C_gmm_iv_dot' = 3L,
    'C_gmm_iv_norm' = 3L,
    'C_gmm_iv_mat_residual' = 4L,
    'C_gmm_iv_eigen_residual' = 4L,
    'C_gmm_iv_mat_product' = 5L,
    'C_gmm_ipw_if_chunk' = 18L,
    'C_gmm_cp08_support_chunk' = 12L,
    'C_gmm_cp08_pair_rows' = 9L,
    'C_gmm_cp08_fit_pass1' = 15L,
    'C_gmm_cp08_fit_pass2' = 15L,
    'C_gmm_cp08_or_if_chunk' = 17L,
    'C_gmm_cp08_dr_if_prepare' = 14L,
    'C_gmm_cp08_dr_if_assemble' = 17L,
    'C_gmm_if_unit_batch' = 10L,
    'C_gmm_cp08_sparse_pair_rows' = 6L,
    'C_gmm_cp08_dr_treatment_batch' = 19L,
    'C_gmm_cp08_payload' = 3L
  )
  expect_identical(sort(names(reg)), sort(names(expected_registry)))
  expect_identical(vapply(reg[names(expected_registry)], function(x) as.integer(x$numParameters), 0L), expected_registry)
  expect_identical(as.integer(reg[['C_gmm_iv_mat_product']]$numParameters), 5L)
  expect_identical(as.integer(reg[['C_gmm_ipw_if_chunk']]$numParameters), 18L)

  syms <- gmm_require_native_dll()
  expect_true(all(c('C_gmm_iv_mat_product') %in% names(syms)))
  expect_true(.Call(syms$C_gmm_test_arithmetic_probe, -1L, 2^-1074, 1.0, 1.0 + 2^-51, 0.0))
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, 999L, 2^-1074, 1.0, 1.0 + 2^-51, 0.0),
               class = 'cdid_gmm_native_error')
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, 0.0, 1.0, 1.0 + 2^-51, 0.0),
               class = 'cdid_gmm_native_error')
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, 2^-1074, 1.0 + 2^-52, 1.0 + 2^-51, 0.0),
               class = 'cdid_gmm_native_error')
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, 2^-1074, 1.0, 1.0 + 2^-52, 0.0),
               class = 'cdid_gmm_native_error')
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, 2^-1074, 1.0, 1.0 + 2^-51, 2^-104),
               class = 'cdid_gmm_native_error')

  fx <- audit_fixture()
  # A process option cannot disable or silently switch this explicit backend.
  old_opt <- options(cdid.gmm.disable_native = TRUE)
  on.exit(options(old_opt), add = TRUE)
  expect_true(all(is.finite(gmm_iv_dot(c(1, 2), c(3, 4)))))
  expect_true(all(is.finite(gmm_iv_norm(c(1, 2)))))
  expect_true(all(is.finite(gmm_iv_mat_residual(diag(2), diag(2), diag(2))$norm)))
  expect_true(all(is.finite(unlist(gmm_iv_mat_product(diag(2), diag(2), diag(2), diag(2))))))
  expect_true(certify_symmetric_eigen_bounds(diag(2))$valid)
  expect_true(audit_solve(fx, 'optimal')$numerical_guard$certified)
  options(old_opt)

  # Check an actually missing backend in a separate process. Never invalidate an
  # installed package namespace by manually unloading its registered DLL.
  probe <- tempfile(fileext='.R')
  on.exit(unlink(probe),add=TRUE)
  writeLines(c(paste('cdid_gmm_abort <-',paste(deparse(cdid_gmm_abort),collapse='\n')),
    paste('gmm_require_native_dll <-',paste(deparse(gmm_require_native_dll),collapse='\n')),
    "e<-tryCatch(gmm_require_native_dll(),error=identity)",
    "stopifnot(inherits(e,'cdid_gmm_native_error'),'cdid' %in% names(getLoadedDLLs())==FALSE)",
    "cat('missing DLL refused without autoload or fallback\\n')"),probe)
  output <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(probe)),stdout=TRUE,stderr=TRUE)
  expect_null(attr(output,'status'))
  expect_match(paste(output,collapse=' '),'missing DLL refused')
  expect_identical(gmm_native_backend_info()$backend, 'native_c_ieee754_binary64')
})

test_that('Native C interval endpoints match frozen pure-R reference bit-for-bit across exhaustive and adversarial inputs', {
  ref_path <- file.path('reference', 'chained_gmm_core_pure_r_reference.R')
  expect_true(file.exists(ref_path))
  expect_identical(as.integer(file.info(ref_path)$size), 106617L)
  expect_identical(unname(tools::md5sum(ref_path)), '21307628a8060a6c223dcabfc0c74600')

  pure_r <- new.env(parent = globalenv())
  for (fn in c('construct_chained_design', 'root_path_certificate')) {
    assign(fn, get(fn, envir = parent.env(environment())), envir = pure_r)
  }
  sys.source(ref_path, envir = pure_r)

  raw_bits <- function(x) writeBin(as.double(x), raw(), size = 8L, endian = 'little')
  syms <- gmm_require_native_dll()

  vals <- c(
    0.0, -0.0, 2^-1074, -2^-1074, 2^-1073, -2^-1073,
    2^-1023, -2^-1023, 2^-1022, -2^-1022,
    2^-600, -2^-600, 0.1, -0.1, 1.0 - 2^-53, -(1.0 - 2^-53),
    1.0, -1.0, 1.0 + 2^-52, -(1.0 + 2^-52),
    3.141592653589793, -3.141592653589793,
    2^500, -2^500, 2^511, -2^511,
    .Machine$double.xmax, -.Machine$double.xmax, Inf, -Inf
  )
  for (v in vals) {
    for (up in c(TRUE, FALSE)) {
      c_nx <- .Call(syms$C_gmm_next, v, up)
      r_nx <- pure_r$gmm_next(v, up)
      expect_identical(raw_bits(c_nx), raw_bits(r_nx))
    }
  }

  fin_vals <- vals[is.finite(vals) & abs(vals) <= 2^511]
  for (i in seq_along(fin_vals)) {
    for (j in seq_along(fin_vals)) {
      vx <- fin_vals[i]; vy <- fin_vals[j]
      expect_identical(raw_bits(gmm_iv_dot(vx, vy)), raw_bits(pure_r$gmm_iv_dot(vx, vy)))
      lo_v <- min(vx, vy); hi_v <- max(vx, vy)
      expect_identical(raw_bits(gmm_iv_norm(lo_v, hi_v)), raw_bits(pure_r$gmm_iv_norm(lo_v, hi_v)))
    }
    expect_identical(raw_bits(gmm_iv_norm(fin_vals[i])), raw_bits(pure_r$gmm_iv_norm(fin_vals[i])))
  }

  adv_vecs <- list(
    list(x = c(0, 0, 0, 0), y = c(1, -1, 2^-1074, 2^500)),
    list(x = c(2^500, 2^-1074, -2^500, 1), y = c(1, 1, 1, -2^-1074)),
    list(x = c(2^-1074, -2^-1074, 2^-1022, -2^-1022, 1 - 2^-53), y = c(1 + 2^-52, 2^500, -2^500, 0.1, -0.1)),
    list(x = sin(1:11) * c(rep(2^-600, 5), 1, rep(2^400, 5)), y = cos(1:11) * c(rep(2^400, 5), 1, rep(2^-600, 5)))
  )
  for (av in adv_vecs) {
    expect_identical(raw_bits(gmm_iv_dot(av$x, av$y)), raw_bits(pure_r$gmm_iv_dot(av$x, av$y)))
    expect_identical(raw_bits(gmm_iv_norm(av$x)), raw_bits(pure_r$gmm_iv_norm(av$x)))
    lo_a <- pmin(av$x, av$y); hi_a <- pmax(av$x, av$y)
    expect_identical(raw_bits(gmm_iv_norm(lo_a, hi_a)), raw_bits(pure_r$gmm_iv_norm(lo_a, hi_a)))
  }

  X_adv <- matrix(c(2^-1074, -1.0, 2^200, 0.1, -2^-1022, 1 + 2^-52,
                    -0.1, 2^-500, 1 - 2^-53, 0.0, -2^200, 3.5), nrow = 3, ncol = 4)
  Y_adv <- matrix(c(1.0, -2^200, 0.25, 2^-1074, -3.0, 0.5, 2^-1022, -1.0), nrow = 4, ncol = 2)
  D_adv <- X_adv %*% Y_adv
  res_c <- gmm_iv_mat_residual(X_adv, Y_adv, D_adv)
  res_r <- pure_r$gmm_iv_mat_residual(X_adv, Y_adv, D_adv)
  expect_identical(raw_bits(res_c$lo), raw_bits(res_r$lo))
  expect_identical(raw_bits(res_c$hi), raw_bits(res_r$hi))
  expect_identical(raw_bits(res_c$norm), raw_bits(res_r$norm))

  S_adv <- crossprod(matrix(c(2, -1, 0.5, 2^-30, 0, 3, -0.25, 1, -1, 0.5, 4, -0.5, 0.2, -0.3, 0.7, 1.5), 4, 4))
  for (j in 1:4) for (i in 1:j) S_adv[j, i] <- S_adv[i, j]
  eb_c <- certify_symmetric_eigen_bounds(S_adv, 'S_adv')
  eb_r <- pure_r$certify_symmetric_eigen_bounds(S_adv, 'S_adv')
  for (fld in c('lambda_min_est', 'lambda_max_est', 'lambda_min_cert', 'lambda_max_cert',
                'delta_Q_est', 'delta_Q_cert', 'eps_eig_est', 'eps_eig_cert', 'pert_cert')) {
    expect_identical(raw_bits(eb_c[[fld]]), raw_bits(eb_r[[fld]]))
  }

  fx <- audit_fixture()
  for (wt in c('optimal', 'identity')) {
    for (sc in c(2^-6, 1.0, 2^6)) {
      fxs <- fx
      fxs$mhat <- fx$mhat * sc
      fxs$V <- fx$V * (sc^2)
      fxs$Phi <- fx$Phi * sc
      acc <- audit_accuracy()
      if (wt == 'identity') {
        acc$target_error_unit <- 'abs'
        acc$max_target_error_se <- NA_real_
        acc$min_certified_lambda_rel <- NA_real_
        acc$min_rcond_R <- NA_real_
        acc$require_matching_scale <- FALSE
      }
      sol_c <- solve_chained_gmm_core(
        fxs$design, fxs$mhat, fxs$V, audit_provenance(), fxs$targets, wt,
        audit_replay(fxs$Phi), acc, audit_bounds(), 'lowest',
        request_param_cov = TRUE, materialize_target_if = TRUE
      )
      sol_r <- pure_r$solve_chained_gmm_core(
        fxs$design, fxs$mhat, fxs$V, audit_provenance(), fxs$targets, wt,
        audit_replay(fxs$Phi), acc, audit_bounds(), 'lowest',
        request_param_cov = TRUE, materialize_target_if = TRUE
      )
      expect_identical(sol_c$numerical_guard$backend, 'native_c_ieee754_binary64')
      expect_identical(sol_c$provenance$native_backend$backend, 'native_c_ieee754_binary64')
      for (col in c('estimate', 'se_sandwich', 'se_criterion_conditional', 'target_error_est', 'target_error_cert')) {
        expect_identical(raw_bits(sol_c$targets[[col]]), raw_bits(sol_r$targets[[col]]))
      }
      expect_identical(raw_bits(sol_c$sandwich_cov), raw_bits(sol_r$sandwich_cov))
      expect_identical(raw_bits(sol_c$raw_Z1tZ1), raw_bits(sol_r$raw_Z1tZ1))
      expect_identical(raw_bits(sol_c$raw_quadratic_residual), raw_bits(sol_r$raw_quadratic_residual))
      for (fld in c('theta_gauge_coords', 'y_multiplier', 'Lambda', 'Mc', 'Z1', 'U_qr', 'param_cov_gauge')) {
        expect_identical(raw_bits(sol_c$internal_gauge[[fld]]), raw_bits(sol_r$internal_gauge[[fld]]))
      }
      for (sec in c('cholesky', 'whitening_and_qr', 'primal_kkt', 'targets')) {
        nms <- names(sol_r$numerical_guard[[sec]])
        for (nm in nms) {
          val_r <- sol_r$numerical_guard[[sec]][[nm]]
          if (is.numeric(val_r)) {
            expect_identical(raw_bits(sol_c$numerical_guard[[sec]][[nm]]), raw_bits(val_r))
          }
        }
      }
    }
  }
})

test_that('Native R and C entrypoints validate dimensions, types, finiteness, aliasing, bounds and overflow before allocation', {
  syms <- gmm_require_native_dll()
  v <- c(1.0, 2.0, 3.0)
  w <- c(4.0, 5.0, 6.0)
  M2 <- matrix(c(1.0, 2.0, 3.0, 4.0), 2, 2)
  M2b <- matrix(c(5.0, 6.0, 7.0, 8.0), 2, 2)
  M2c <- matrix(c(9.0, 10.0, 11.0, 12.0), 2, 2)

  # Read-only aliases are valid and do not allocate defensive copies.
  pure <- new.env(parent=globalenv())
  sys.source(file.path('reference','chained_gmm_core_pure_r_reference.R'),pure)
  expect_identical(gmm_iv_dot(v,v),pure$gmm_iv_dot(v,v))
  expect_identical(.Call(syms$C_gmm_iv_dot,v,v,NULL),pure$gmm_iv_dot(v,v))
  expect_identical(gmm_iv_norm(v,v),pure$gmm_iv_norm(v,v))
  expect_identical(.Call(syms$C_gmm_iv_norm,v,v,NULL),pure$gmm_iv_norm(v,v))
  expect_identical(gmm_iv_mat_residual(M2,M2,M2c),pure$gmm_iv_mat_residual(M2,M2,M2c))
  expect_identical(gmm_iv_mat_residual(M2,M2b,M2),pure$gmm_iv_mat_residual(M2,M2b,M2))
  expect_identical(.Call(syms$C_gmm_iv_mat_residual,M2,M2,M2c,NULL),pure$gmm_iv_mat_residual(M2,M2,M2c))
  E <- .Call(syms$C_gmm_iv_eigen_residual,M2,M2,c(1.0,2.0),NULL)
  lo <- hi <- matrix(0,2,2)
  for(i in 1:2) for(j in 1:2) {
    z<-pure$gmm_iv_sub(pure$gmm_iv_dot(M2[i,],M2[,j]),pure$gmm_iv_mul(pure$gmm_iv(M2[i,j]),pure$gmm_iv(c(1,2)[j])))
    lo[i,j]<-z[1];hi[i,j]<-z[2]
  }
  expect_identical(E,list(lo=lo,hi=hi,norm=pure$gmm_iv_norm(lo,hi)))
  expect_error(gmm_iv_residual(M2,c(1.0,2.0),c(1.0,2.0)[1:2]),regexp=NA)
  u_alias <- c(1.0,2.0)
  expect_identical(gmm_iv_residual(M2,u_alias,u_alias),pure$gmm_iv_residual(M2,u_alias,u_alias))

  # Dimension, empty, and conformability refusals
  expect_error(gmm_iv_dot(numeric(0), numeric(0)), class = 'cdid_gmm_bound_error')
  expect_error(gmm_iv_dot(c(1.0, 2.0), c(1.0, 2.0, 3.0)), class = 'cdid_gmm_bound_error')
  expect_error(.Call(syms$C_gmm_iv_dot, numeric(0), numeric(0), NULL), class = 'cdid_gmm_bound_error')
  expect_error(.Call(syms$C_gmm_iv_dot, c(1.0, 2.0), c(1.0, 2.0, 3.0), NULL), class = 'cdid_gmm_bound_error')
  expect_error(gmm_iv_norm(numeric(0)), class = 'cdid_gmm_bound_error')
  expect_error(gmm_iv_norm(c(1.0, 2.0), c(1.0, 2.0, 3.0)), class = 'cdid_gmm_bound_error')
  expect_error(gmm_iv_norm(M2, c(1.0, 2.0, 3.0, 4.0)), class = 'cdid_gmm_bound_error')
  expect_error(gmm_iv_mat_residual(M2, matrix(1.0, 3, 2), M2c), class = 'cdid_gmm_bound_error')
  expect_error(.Call(syms$C_gmm_iv_eigen_residual, M2, M2b, c(1.0, 2.0, 3.0), NULL), class = 'cdid_gmm_bound_error')

  # Caller max_matrix_entries budget refusals before indexing/allocation
  expect_error(gmm_iv_dot(v, w, max_matrix_entries = 4), class = 'cdid_gmm_bound_error')
  expect_error(.Call(syms$C_gmm_iv_dot, v, w, 4), class = 'cdid_gmm_bound_error')
  expect_error(gmm_iv_norm(v, w, max_matrix_entries = 5), class = 'cdid_gmm_bound_error')
  expect_error(gmm_iv_mat_residual(M2, M2b, M2c, max_matrix_entries = 7), class = 'cdid_gmm_bound_error')
  expect_error(.Call(syms$C_gmm_iv_eigen_residual, M2, M2b, c(1.0, 2.0), 7), class = 'cdid_gmm_bound_error')
  expect_error(certify_symmetric_eigen_bounds(diag(3), max_matrix_entries = 20), class = 'cdid_gmm_bound_error')

  # Type, object, finiteness, and inverted-interval refusals
  expect_error(.Call(syms$C_gmm_iv_dot, 1:3, w, NULL), class = 'cdid_gmm_accuracy_error')
  expect_error(gmm_iv_dot(structure(v, class = 'custom_obj'), w), class = 'cdid_gmm_accuracy_error')
  expect_error(gmm_iv_dot(c(1.0, NA_real_, 3.0), w), class = 'cdid_gmm_accuracy_error')
  expect_error(gmm_iv_dot(c(1.0, Inf, 3.0), w), class = 'cdid_gmm_accuracy_error')
  expect_error(gmm_iv_norm(c(2.0, 1.0), c(1.0, 2.0)), class = 'cdid_gmm_accuracy_error')
  expect_error(.Call(syms$C_gmm_iv_norm, c(2.0, 1.0), c(1.0, 2.0), NULL), class = 'cdid_gmm_accuracy_error')
  expect_error(.Call(syms$C_gmm_next, NA_real_, TRUE), class = 'cdid_gmm_accuracy_error')
  expect_error(.Call(syms$C_gmm_next, 1.0, NA), class = 'cdid_gmm_accuracy_error')

  # Validated arithmetic overflow refusals in compiled C hotspots
  huge <- c(.Machine$double.xmax, .Machine$double.xmax)
  huge2 <- c(.Machine$double.xmax, .Machine$double.xmax)
  expect_error(gmm_iv_dot(huge, huge2), class = 'cdid_gmm_accuracy_error')
  expect_error(gmm_iv_norm(huge), class = 'cdid_gmm_accuracy_error')
  Mh1 <- matrix(.Machine$double.xmax, 2, 2)
  Mh2 <- matrix(.Machine$double.xmax, 2, 2)
  expect_error(gmm_iv_mat_residual(Mh1, Mh2, M2c), class = 'cdid_gmm_accuracy_error')
  expect_error(.Call(syms$C_gmm_iv_eigen_residual, Mh1, Mh2, c(1.0, 1.0), NULL), class = 'cdid_gmm_accuracy_error')
})

# D12 addition: verify explicit loading-budget control validation (NULL absence state vs invalid NA/nonpositive),
# outward G_lambda inclusion of both adjoint blocks s1 and s2, and non-certification of complete inference.
test_that('Independent D12 loading budget control, adjoint s1/s2 discrimination in G_lambda, and pending inference status', {
  fx <- audit_fixture()
  # NULL is an explicit absence state: point certified, loading incomplete, complete inference FALSE
  s_null <- audit_solve(fx, 'optimal', max_loading_error_abs = NULL)
  expect_true(s_null$point_certified)
  expect_false(s_null$loading_certified)
  expect_false(s_null$loading_budget_supplied)
  expect_identical(s_null$loading_validation_status, 'incomplete_missing_caller_budget')
  expect_false(s_null$complete_inference_certified)

  # Explicit positive loading budget certifies stored-array loadings while complete inference remains FALSE
  s_load <- audit_solve(fx, 'optimal', max_loading_error_abs = 1e-10)
  expect_true(s_load$point_certified)
  expect_true(s_load$loading_certified)
  expect_true(s_load$loading_budget_supplied)
  expect_identical(s_load$loading_validation_status, 'certified_conditional_stored_array')
  expect_false(s_load$complete_inference_certified)

  # Invalid max_loading_error_abs values (NA, negative, Inf, vector) refuse with cdid_gmm_control_error
  for (bad_lb in list(NA_real_, -1e-10, Inf, c(1e-10, 1e-10), '1e-10')) {
    acc_bad <- audit_accuracy('optimal')
    acc_bad$max_loading_error_abs <- bad_lb
    expect_error(validate_gmm_accuracy_control(acc_bad, 'optimal'), class = 'cdid_gmm_control_error')
  }

  # Too-tight loading budget refuses even when point budget passes
  expect_error(audit_solve(fx, 'optimal', max_loading_error_abs = 1e-30), class = 'cdid_gmm_accuracy_error')
})
