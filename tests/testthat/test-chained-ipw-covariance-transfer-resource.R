require_cdid_namespace <- function() {
  if (!isNamespaceLoaded("cdid") && !requireNamespace("cdid", quietly = TRUE)) {
    stop("Actual installed 'cdid' namespace is required; no fallback is permitted in packaged tests.", call. = FALSE)
  }
  ns <- asNamespace("cdid")
  if (!isNamespace(ns) || !environmentIsLocked(ns)) stop("'cdid' namespace must be a locked R namespace.", call. = FALSE)
  ns
}
require_cdid_namespace()

test_that("covariance output cap includes retained S t and every declared cluster array", {
  ns <- require_cdid_namespace()
  fx <- dget(testthat::test_path("fixtures", "ipw_covariance_transfer_fixtures.R"))
  b <- fx$build_14unit_3link_4bd(ns, r_ball = 1e-6, max_chunk_rows = 4L, max_operator_block_rows = 5L)
  ctrl <- b$control
  ctrl$max_cluster_rows <- 100L
  op <- ns$ipw_covariance_operator_spec("inference", "cluster", rep(1,100), TRUE, b$kappa,
    as.double(seq_along(b$id)), as.double(seq_len(100)))
  hdr <- ns$ipw_cov_transfer_validate_links_headers(b$links)
  pre <- ns$preflight_ipw_covariance_transfer(14,3,2,"cluster",100,hdr$p_vec,hdr$m_vec,hdr$certificate_controls,ctrl)
  # Independently counted:6*d²+d*L+4*d+3*Ncluster*d+64*L+256 =1086.
  expect_identical(pre$output_entries,1086)
  expect_identical(pre$cluster_buffer_entries,1000)
  expect_identical(pre$b_by_d_entries,14*4*2+12*5*2)
  expect_identical(pre$d_by_L_entries,4*2*3)
  ctrl$max_output_entries <- pre$output_entries
  ans <- ns$certify_ipw_covariance_transfer(b$links,b$candidate_rows,b$frame,b$C,op,diag(2),ctrl)
  arrays <- list();count <- 0
  visit <- function(x) {
    if (is.numeric(x)||is.logical(x)) {
      if (!any(vapply(arrays,function(y)ns$gmm_sexp_aliases(x,y),logical(1)))) {
        arrays[[length(arrays)+1L]] <<- x
        count <<- count+length(x)
      }
    } else if (is.list(x)) for(y in x)visit(y)
  }
  visit(ans)
  expect_lte(count,ctrl$max_output_entries)
  expect_identical(dim(ans$S_unnormalized_iv),c(2L,2L,2L))
  expect_identical(dim(ans$t_iv),c(2L,2L))
  expect_identical(ans$cluster_summary$n_zero_unit_declared_clusters,86L)
  expect_true(all(ans$cluster_summary$cluster_point[15:100,]==0))
  expect_true(all(ans$cluster_summary$cluster_interval_lo[15:100,]==0))
  expect_true(all(ans$cluster_summary$cluster_interval_hi[15:100,]==0))
  expect_identical(ans$preflight$total_with_existing_4bd_entries,
    ans$preflight$existing_4bd_retained_entries+ans$preflight$concurrent_entries)
  calls <- 0L
  cb <- ns$ipw_candidate_rows_spec(function(s,e){calls <<- calls+1L;b$candidate_rows$callback(s,e)},b$candidate_rows$link_ids)
  ctrl$max_output_entries <- pre$output_entries-1
  expect_error(ns$certify_ipw_covariance_transfer(b$links,cb,b$frame,b$C,op,diag(2),ctrl),class="cdid_ipw_covariance_transfer_bound_error")
  expect_identical(calls,0L)
  ctrl <- b$control;ctrl$max_matrix_entries <- 100L
  expect_error(ns$certify_ipw_covariance_transfer(rep(list(NULL),1000),cb,b$frame,b$C,b$op_gmm,diag(2),ctrl),
    regexp="link header vectors",class="cdid_ipw_covariance_transfer_bound_error")
  expect_identical(calls,0L)
})

# The packaged entrypoint below requires the real installed namespace. The pure
# check function also permits an explicitly disclosed source-module audit.
covariance_lifetime_checks <- function(ns, fx) {
  checks <- list(); details <- list()
  ck <- function(name, value) checks[[name]] <<- isTRUE(value)
  cc <- ns$logit_certificate_control(18L, 40, 16L, 14L, 5000L, 500000L)
  lc <- ns$ipw_link_enclosure_control(1000L, 1e7, 1e7, 1e12, 1e7)
  cases <- list(
    cluster = list(N=100, M=10000, d=2, b=1, q=1, L=1, base=121447, old_live=140000),
    unit = list(N=2000, M=2000, d=1, b=1000, q=1, L=100, base=389400, old_live=400000)
  )
  for (name in names(cases)) {
    x <- cases[[name]]
    ctrl <- ns$ipw_covariance_transfer_control(x$b,x$q,x$M,1e7,1e12,1e7,lc)
    pre <- ns$preflight_ipw_covariance_transfer(x$N,x$L,x$d,name,x$M,
      rep(3,x$L),rep(3,x$L),rep(list(cc),x$L),ctrl)
    # Independent payload sum:3bL+b,14bd+12qd,4dL,24d²+4d+256,
    # unit7N+66L / cluster8N+2M+66L,5Md if clustered,H=53b+104.
    expected <- 3*x$b*x$L+x$b+14*x$b*x$d+12*x$q*x$d+4*x$d*x$L+
      24*x$d^2+4*x$d+256+53*x$b+104+66*x$L+
      if(name=='unit')7*x$N else 8*x$N+2*x$M+5*x$M*x$d
    ck(paste0(name,'_exact_preflight'),identical(pre$base_component_concurrent_entries,expected)&&expected==x$base)
    ck(paste0(name,'_original_liveness_counterexample'),x$old_live>expected)
    ck(paste0(name,'_single_chunk_charge'),pre$b_by_L_entries==3*x$b*x$L+x$b)
    tight <- ctrl;tight$max_matrix_entries <- expected-1
    err <- tryCatch(ns$preflight_ipw_covariance_transfer(x$N,x$L,x$d,name,x$M,
      rep(3,x$L),rep(3,x$L),rep(list(cc),x$L),tight),error=identity)
    ck(paste0(name,'_unchanged_cap_refusal'),inherits(err,'cdid_ipw_covariance_transfer_bound_error'))
    details[[name]] <- list(preflight=pre,original_live_lower_bound=x$old_live)
  }
  # Check absence of old chunk aliases at the precise next-fetch boundary.
  dead <- c('ch_data','Phi_block','E_eval_block','E_root_block','Uhat','res_U',
    'U_lo','U_hi','zero_b_d','res_E_eval','res_E_root','ae_eval_blk','ae_root_blk')
  seen <- list()
  fetch <- function(s,e) {
    seen[[length(seen)+1L]] <<- intersect(dead,ls(parent.frame(),all.names=TRUE))
    z <- matrix(c(-1,1),nrow=e-s+1L,ncol=1L)
    list(certified=TRUE,Phi_block=z,E_eval_block=z*0,E_root_block=z*0,
      chunk_4bd_ops=0,max_direct_star_dist=0)
  }
  op <- list(row_level='unit',n_clusters=4L,omega=rep(1,4),centered=TRUE,kappa=.25)
  ans <- ns$ipw_cov_transfer_stream_kernel(4L,1L,1L,2L,1L,fetch,matrix(1),op,matrix(1),10000L)
  ck('kernel_two_chunks_completed',ans$enclosure_certified&&length(seen)==2L)
  ck('kernel_previous_chunk_released',all(lengths(seen)==0L))
  details$kernel_boundary_bindings <- seen
  # Spy only on the consumed-result lifetime; call the accepted 4bd helper with
  # identical arguments. No fitter, alternative arithmetic, or namespace mask.
  b <- fx$build_14unit_3link_4bd(ns,r_ball=1e-6,max_chunk_rows=4L,max_operator_block_rows=2L)
  retained <- list(); probe <- new.env(parent=environment(ns$certify_ipw_covariance_transfer))
  probe$enclose_ipw_link_if_chunk <- function(...) {
    retained[[length(retained)+1L]] <<- intersect(c('ch_l','sep_rows'),ls(parent.frame(),all.names=TRUE))
    ns$enclose_ipw_link_if_chunk(...)
  }
  cert <- ns$certify_ipw_covariance_transfer;environment(cert) <- probe
  res <- cert(b$links,b$candidate_rows,b$frame,b$C,b$op_gmm,b$V_cand_gmm,b$control)
  ck('link_spy_complete_no_fit',res$enclosure_certified&&length(retained)==12L)
  ck('previous_link_result_released',all(lengths(retained)==0L))
  details$link_boundary_bindings <- retained
  # Evaluate the actual scalar diagnostic block against independent finite
  # endpoint fixtures. Signed zeros/subnormal/positive-overflow rows discriminate.
  exprs <- as.list(body(ns$ipw_cov_transfer_stream_kernel))[-1L]
  idx <- which(vapply(exprs,function(e)is.call(e)&&identical(e[[1L]],as.name('<-'))&&
    identical(e[[2L]],as.name('n_zero_if_clusters')),logical(1L)))
  ck('scalar_scan_present',length(idx)==1L)
  out_list <- exprs[[length(exprs)]]
  cs <- out_list[['cluster_summary']][[3L]][[2L]]
  diag_expr <- cs[['n_zero_if_clusters']]
  ck('returned_count_is_scalar_binding',identical(diag_expr,as.name('n_zero_if_clusters')))
  if(length(idx)==1L) {
    s <- new.env(parent=environment(ns$ipw_cov_transfer_stream_kernel))
    s$is_cluster <- TRUE;s$N_c <- 5L;s$d <- 2L;s$b_op <- 2L
    tiny <- .Machine$double.xmin*.Machine$double.eps
    s$cluster_z_lo <- matrix(c(0,-tiny,-1e308,0,0,-0.0,0,1e308,0,tiny),5L,2L)
    s$cluster_z_hi <- s$cluster_z_lo
    # Rows1/4 are all zero. Row3's old abs-sum overflows, but is nonzero.
    eval(exprs[[idx]],s);eval(exprs[[idx+1L]],s)
    ck('finite_zero_count_exact',identical(s$n_zero_if_clusters,2L))
    ck('finite_old_count_identical',identical(s$n_zero_if_clusters,
      sum(rowSums(abs(s$cluster_z_lo)+abs(s$cluster_z_hi))==0)))
  } else {
    ck('finite_zero_count_exact',FALSE);ck('finite_old_count_identical',FALSE)
  }
  list(checks=checks,details=details)
}

require_covariance_lifetime_namespace <- function() {
  if(!requireNamespace('cdid',quietly=TRUE))stop("Actual installed cdid namespace required; no source fallback.",call.=FALSE)
  ns <- asNamespace('cdid')
  stopifnot(isNamespace(ns),environmentIsLocked(ns))
  ns
}

test_that('covariance chunk lifetimes and scalar zero diagnostic honor the declared payload envelope', {
  ns <- require_covariance_lifetime_namespace()
  fx <- dget(testthat::test_path('fixtures','ipw_covariance_transfer_fixtures.R'))
  got <- covariance_lifetime_checks(ns,fx)
  for(name in names(got$checks))expect_true(got$checks[[name]],info=name)
})
