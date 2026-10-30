# D26 fixed-matrix algebra only: independent dense full-catalog references,
# literal closed forms, actual automatic core outputs and resource boundaries.
# No nuisance/data fits, RNG, DLL substitution, bootstrap or fixture retuning.
d26_accuracy <- function(weighting = "optimal", alpha = 1e-4, trigger = 1e-10, strategy = TRUE) {
  gmm_accuracy_control(method = "conventional_chol_qr", target_error_unit = "abs",
    max_target_error_abs = .Machine$double.xmax, max_target_error_se = NA_real_,
    max_orthogonality_defect = 1e-8, min_certified_lambda_rel = min(trigger, 1e-10),
    max_factor_residual_rel = 1e-10, min_rcond_R = 1e-10, min_rcond_U = 1e-10,
    min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10,
    max_loading_residual_rel = 1e-10, require_matching_scale = FALSE,
    covariance_strategy = if (strategy) gmm_covariance_strategy(condition_trigger = trigger,
      ridge_relative = alpha, weighting_requested = weighting,
      weighting_canonical = resolve_gmm_weighting_choice(weighting)$canonical) else NULL)
}
d26_bounds <- function(entries = 1e7) gmm_bounds_control(entries, 4, 20, 4, TRUE, 1000)
d26_provenance <- function(n = 4, V = NULL, internal = FALSE, reduction = list(mode = "none")) {
  p <- gmm_covariance_provenance("unit", n, n, n, TRUE, FALSE, "unweighted", 1/n,
    "D26:fixed", "D26:fixed", reduction)
  if (internal) p$complete_unit_if_gram <- ipw_seal(list(covariance = V, n_units = p$n_units,
    nuisance_fit_id = p$nuisance_fit_id), "complete_unit_if_gram")
  p
}
# Independent reference: dense full-space normal equations are fine for these
# tiny fixed matrices. Production uses Cholesky/QR; no production helper here.
d26_reference <- function(H, m, V, Q) {
  A <- solve(crossprod(H, Q %*% H), crossprod(H, Q))
  list(point = drop(A %*% m), A = A, sandwich = A %*% V %*% t(A),
    bread = solve(crossprod(H, Q %*% H)))
}
d26_handler <- function(H, m, V, G = diag(nrow(V)), M = nrow(V), alpha = 1e-4,
                        provenance = d26_provenance(), trigger = 1e-10) {
  k <- ncol(H); L <- nrow(H)
  gram <- assemble_gmm_gram_matrix(L, "optimal", list(G = G, n_links_all = M), 1e7)
  a <- d26_accuracy(alpha = alpha, trigger = trigger)
  b <- diagnose_gmm_transformed_covariance(V, gram, provenance, a$covariance_strategy, 1e7)
  h <- solve_gmm_covariance_handler(V, m, H,
    list(n_links = L, rank_k = k, link_ids = paste0("m", seq_len(L)), keep_param_ids = paste0("p", seq_len(k))),
    list(n_targets = k, estimable = rep(TRUE, k), estimable_indices = seq_len(k),
      C_gauge = diag(k)), gram, b, a$covariance_strategy, a, 1e7)
  list(point = unname(h$theta_keep), A = t(unname(h$Lambda)),
    sandwich = crossprod(h$Lambda, V %*% h$Lambda), bread = h$raw_Z1tZ1,
    lambda = h$lambda, diagnosis = b, solve = h)
}
d26_core <- function(V = diag(c(1, 0)), m = c(3, 2), weighting = "optimal",
                     alpha = 1e-4, trigger = 1e-10, Phi = cbind(c(1,1,-1,-1), 0),
                     clustered = FALSE, strategy = TRUE, n_positive = nrow(Phi),
                     accuracy = NULL, targets = NULL, reduction = list(mode = "none")) {
  des <- build_chained_design(data.frame(g = c(2,2), s = c(1,1), t = c(2,2), ref = c(1,1),
    functional_key = c("a", "b")))
  n <- nrow(Phi)
  p <- d26_provenance(n, reduction = reduction); p$n_positive_contributing_rows <- n_positive
  replay <- gmm_if_replay_spec(function(consume) {
    for (i in seq_len(n)) consume(list(unit_id = i, phi = Phi[i,,drop=FALSE],
      cluster_id = if (clustered) ceiling(i/2) else NULL, nuisance_fit_id = "D26:fixed"))
  }, "unit", clustered, n, if (clustered) n/2 else n, if (clustered) n/2 else n,
  TRUE, FALSE, "unweighted", 1/n, "D26:fixed", "D26:fixed", seq_len(n))
  if (is.null(targets)) targets <- list(gmm_target_spec("theta", 1L, 1))
  if (is.null(accuracy)) accuracy <- d26_accuracy(weighting, alpha, trigger, strategy)
  solve_chained_gmm_core(des, m, V, p, targets, weighting, replay, accuracy,
    d26_bounds(), "lowest", request_param_cov = TRUE, materialize_target_if = TRUE)
}

test_that("D26 retains supplied closed forms and actual sandwich rather than bread", {
  H <- matrix(c(1,1), 2); V <- diag(c(1,0)); m <- c(3,2)
  r <- d26_handler(H,m,V,alpha=2) # s=1/2, lambda=1
  expect_equal(r$lambda, 1, tolerance=1e-14)
  expect_equal(r$point, 7/3, tolerance=1e-14)
  expect_equal(as.double(r$A), c(1/3,2/3), tolerance=1e-14)
  expect_equal(as.double(r$sandwich), 1/9, tolerance=1e-14)
  expect_equal(as.double(r$bread), 2/3, tolerance=1e-14)
  expect_equal(drop(r$A %*% H), 1, tolerance=1e-14)
  z <- d26_handler(matrix(c(0,1),2),m,V,alpha=2)
  expect_equal(z$point, 2, tolerance=1e-14)
  expect_equal(as.double(z$A), c(0,1), tolerance=1e-14)
  expect_equal(as.double(z$sandwich), 0, tolerance=1e-14)
  E <- rbind(c(1,0),c(1,0),c(0,1)); Vr <- diag(c(1,2))
  f <- d26_reference(E %*% H, E %*% m, E %*% Vr %*% t(E),
    solve(E %*% Vr %*% t(E) + diag(3)))
  reduced <- d26_handler(H,m,Vr,crossprod(E),3,alpha=3/4)
  naive <- d26_handler(H,m,Vr,alpha=2/3)
  expect_equal(f$point, 8/3, tolerance=1e-14)
  expect_equal(reduced$point, 8/3, tolerance=1e-14)
  expect_equal(naive$point, 13/5, tolerance=1e-14)
  expect_equal(unname(f$A %*% E), reduced$A, tolerance=1e-14)
  expect_equal(unname(reduced$sandwich), matrix(2/3), tolerance=1e-14)
  # A plausible formula mutation (forgetting reconstruction G) is rejected by
  # the literal independent 8/3 target, not a fixture from the tested helper.
  expect_gt(abs(naive$point - 8/3), .05)
})

test_that("D26 sparse reconstruction handles signed nondiagonal G and permuted kept rows", {
  E <- rbind(c(1,0),c(0,1),c(1,1))
  sparse <- list(n_links=3,retained_index=c(1L,2L),dropped_index=3L,
    E=list(i=c(3L,3L),j=c(1L,2L),x=c(1,1)),kept_indices=c(2L,1L))
  g <- assemble_gmm_gram_matrix(2,"optimal",sparse,1e7)
  expect_equal(g$G, matrix(c(2,1,1,2),2), tolerance=1e-14)
  expect_identical(g$kept_indices,c(2L,1L)); expect_false(g$is_diagonal_gram)
  H <- matrix(c(1,1),2); m <- c(3,2); V <- diag(c(1,0))
  r <- d26_handler(H,m,V,crossprod(E),3,alpha=1.5)
  Q <- crossprod(E,solve(E %*% V %*% t(E)+diag(3),E))
  expect_equal(Q,matrix(c(2/3,1/3,1/3,5/3),2),tolerance=1e-14)
  expect_equal(r$diagnosis$scale_s,2/3,tolerance=1e-14)
  expect_equal(r$point,7/3,tolerance=1e-14)
  # Denominator mutation tr(G) changes lambda and fails the unchanged formula.
  expect_gt(abs(sum(diag(V %*% crossprod(E)))/sum(diag(crossprod(E))) - 2/3), .1)
  # Unequal signed rows ensure column permutation actually changes G.
  sparse$E$x <- c(-2,1)
  Ep <- rbind(c(0,1),c(1,0),c(1,-2))
  gp <- assemble_gmm_gram_matrix(2,"optimal",sparse,1e7)
  expect_equal(gp$G,crossprod(Ep),tolerance=1e-14)
  expect_equal(crossprod(gp$C),crossprod(Ep),tolerance=1e-14)
  sparse$E$j[1] <- 1.5
  expect_error(assemble_gmm_gram_matrix(2,"optimal",sparse,1e7),class="cdid_gmm_reduction_error")
})

test_that("D26 full/reduced and equivalent-basis ridge loadings and sandwiches agree", {
  E <- rbind(c(1,0),c(0,1),c(-1,1),c(2,1)); G <- crossprod(E)
  H <- rbind(c(1,0),c(1,1)); m <- c(3,2); V <- diag(c(1,0)); alpha <- .4
  lambda <- alpha*sum(diag(E %*% V %*% t(E)))/nrow(E)
  f <- d26_reference(E%*%H,E%*%m,E%*%V%*%t(E),solve(E%*%V%*%t(E)+diag(lambda,nrow(E))))
  r <- d26_handler(H,m,V,G,nrow(E),alpha)
  expect_equal(r$point,f$point,tolerance=1e-14)
  expect_equal(r$A,unname(f$A%*%E),tolerance=1e-14)
  expect_equal(unname(r$sandwich),unname(f$sandwich),tolerance=1e-14)
  expect_equal(r$A%*%H,diag(2),tolerance=1e-14)
  S <- rbind(c(1,1),c(0,1)); E2 <- E%*%solve(S)
  q <- d26_handler(S%*%H,drop(S%*%m),S%*%V%*%t(S),crossprod(E2),nrow(E),alpha)
  expect_equal(q$point,r$point,tolerance=1e-14)
  expect_equal(q$A%*%S,r$A,tolerance=1e-14)
  expect_equal(unname(q$sandwich),unname(r$sandwich),tolerance=1e-14)
  expect_equal(q$lambda,r$lambda,tolerance=1e-14)
})

test_that("D26 automatic singular/near-singular/zero outputs and ordinary nonactivation are real", {
  f <- d26_core(alpha=2)
  expect_true(f$covariance_handling$handler_invoked)
  expect_equal(f$targets$estimate,7/3,tolerance=1e-14)
  expect_equal(as.double(f$internal_gauge$Lambda),c(1/3,2/3),tolerance=1e-14)
  expect_equal(as.double(f$sandwich_cov),1/9,tolerance=1e-14)
  expect_equal(as.double(f$raw_Z1tZ1),2/3,tolerance=1e-14)
  expect_equal(unname(f$if_projection$materialized_target_if[,1]),c(1,1,-1,-1)/3,tolerance=1e-14)
  expect_true(all(f$targets$estimable)); expect_false(f$point_certified)
  def <- d26_core(); lam <- 1e-4/2
  expect_equal(def$targets$estimate,2+lam/(1+2*lam),tolerance=1e-14)
  expect_equal(def$covariance_handling$lambda,lam,tolerance=1e-14)
  near <- d26_core(V=diag(c(1,1e-12)))
  expect_true(near$covariance_handling$handler_invoked)
  ref <- d26_reference(matrix(c(1,1),2),c(3,2),diag(c(1,1e-12)),
    solve(diag(c(1,1e-12))+diag(1e-4*(1+1e-12)/2,2)))
  expect_equal(near$targets$estimate,ref$point,tolerance=1e-14)
  expect_equal(t(unname(near$internal_gauge$Lambda)),ref$A,tolerance=1e-14)
  pd <- d26_core(V=diag(c(1,2))); old <- d26_core(V=diag(c(1,2)),strategy=FALSE)
  expect_false(pd$covariance_handling$handler_invoked)
  expect_identical(pd$targets,old$targets)
  expect_identical(pd$internal_gauge$Lambda,old$internal_gauge$Lambda)
  expect_identical(pd$raw_quadratic_residual,old$raw_quadratic_residual)
  # Equality belongs to the ordinary branch, if Cholesky succeeds.
  for (rho in c(1e-10/2,1e-10,2e-10)) {
    x <- d26_core(V=diag(c(1,rho)))
    expect_identical(x$covariance_handling$handler_invoked,rho<1e-10)
  }
  # Known rank bound activates even if a roundoff perturbation looks PD.
  rank <- d26_core(V=diag(2),n_positive=2)
  expect_true(rank$covariance_handling$handler_invoked)
  expect_identical(rank$covariance_handling$trigger_reason,"known_unit_covariance_rank_bound")
  z <- d26_core(V=matrix(0,2,2),Phi=matrix(0,4,2))
  expect_equal(z$targets$estimate,2.5,tolerance=1e-14)
  expect_equal(as.double(z$sandwich_cov),0,tolerance=1e-14)
  expect_true(z$covariance_handling$zero_covariance)
  expect_identical(z$covariance_handling$zero_covariance_branch,"scale_free_gram_limit")
  expect_false(z$complete_inference_certified)
  budget <- list(max_units=4,max_clusters=4,max_targets=1,max_draws=1,max_if_batch_rows=4,
    max_if_batch_cells=4,max_sign_batch_rows=1,max_sign_batch_cells=4,max_workspace_cells=1e5,max_workspace_bytes=8e5)
  expect_error(ipw_pipeline_supported_inference(z$targets$estimate,matrix(0,4,1),TRUE,
    multiplier_budget=budget,inference_accuracy=NULL,bstrap=FALSE),class="cdid_ipw_pipeline_scale_error")
})

test_that("D26 relative ridge and explicit identity preserve common positive scale invariance", {
  H <- matrix(c(1,1),2); m <- c(3,2); V <- diag(c(1,0)); G <- matrix(c(2,1,1,2),2)
  r <- d26_handler(H,m,V,G,3,alpha=.1)
  for (a in c(1e-100,.125,8,1e100)) {
    q <- d26_handler(H,m,a*V,G,3,alpha=.1)
    expect_equal(q$point,r$point,tolerance=1e-14)
    expect_equal(q$A,r$A,tolerance=1e-14)
    expect_equal(q$diagnosis$rho,r$diagnosis$rho,tolerance=1e-14)
    expect_equal(q$lambda/a,r$lambda,tolerance=1e-14)
    expect_equal(q$A%*%H,matrix(1),tolerance=1e-14)
    expect_equal(q$sandwich/a,r$sandwich,tolerance=1e-14)
    id <- d26_core(V=a*V,weighting="identity")
    expect_false(id$covariance_handling$handler_invoked)
    expect_equal(id$targets$estimate,2.5,tolerance=1e-14)
    expect_equal(as.double(id$internal_gauge$Lambda),c(.5,.5),tolerance=1e-14)
  }
  z <- d26_handler(H,m,matrix(0,2,2),G,3)
  full_id <- d26_reference(H,m,matrix(0,2,2),G)
  expect_equal(z$point,full_id$point,tolerance=1e-14)
  expect_equal(z$A,full_id$A,tolerance=1e-14)
})

test_that("D26 clustering changes only the actual IF sandwich", {
  u <- d26_core(alpha=2); cl <- d26_core(alpha=2,clustered=TRUE)
  expect_identical(cl$targets$estimate,u$targets$estimate)
  expect_identical(cl$internal_gauge$Lambda,u$internal_gauge$Lambda)
  expect_identical(cl$covariance_handling$lambda,u$covariance_handling$lambda)
  expect_equal(as.double(cl$sandwich_cov),2/9,tolerance=1e-14)
  expect_equal(unname(cl$if_projection$materialized_target_if[,1]),c(2,-2)/3,tolerance=1e-14)
  expect_false(isTRUE(all.equal(cl$sandwich_cov,cl$raw_Z1tZ1)))
})

test_that("D26 numerical negatives require actual Gram provenance and the exact D27 rule", {
  # Recorded arithmetic perturbation of the Gram from Phi=(1,1): epsilon
  # changes the off diagonal and leaves a small negative eigenvalue.
  V <- matrix(c(1,1+2^-50,1+2^-50,1),2)
  p <- d26_provenance(V=V,internal=TRUE)
  g <- assemble_gmm_gram_matrix(2,"optimal",NULL,1e7); s <- gmm_covariance_strategy()
  run <- function(v=V,prov=p) diagnose_gmm_transformed_covariance(v,g,prov,s,1e7)
  b <- run()
  expect_true(b$diagnosed_numerical_negative)
  expect_identical(b$negative_eigenvalue_diagnosis,"consistent_with_roundoff")
  expect_equal(b$b_bound,100*2*.Machine$double.eps*sqrt(sum(V^2)),tolerance=1e-14)
  expect_lte(b$r_eig,b$b_bound); expect_gte(b$d_min,-(b$b_bound+b$r_eig))
  r <- d26_handler(matrix(c(1,1),2),c(3,2),V,provenance=p)
  expect_equal(r$lambda,1e-4,tolerance=1e-14)
  expect_identical(r$diagnosis$B,V) # original B, never a clipped reconstruction
  expect_error(run(prov=d26_provenance()),"provenance",class="cdid_gmm_pd_error")
  fake <- d26_provenance(); fake$gram_provenance_certified <- TRUE
  expect_error(run(prov=fake),"provenance",class="cdid_gmm_pd_error")
  expect_error(run(v=diag(2)),"provenance",class="cdid_gmm_pd_error")
  material <- matrix(c(1,2,2,1),2)
  expect_error(run(material,d26_provenance(V=material,internal=TRUE)),"indefinite",class="cdid_gmm_pd_error")
  # Fault only the local eigensolver seam; production has no caller attributes
  # for falsifying residuals/factor failures, and no fake namespace/DLL.
  faulty <- diagnose_gmm_transformed_covariance
  local <- new.env(parent=environment(faulty))
  local$eigen <- function(...) { z <- base::eigen(...); z$vectors[1,1] <- 0; z }
  environment(faulty) <- local
  expect_error(faulty(V,g,p,s,1e7),"residual",class="cdid_gmm_pd_error")
  expect_error(d26_handler(matrix(c(1,1),2),c(3,2),V,alpha=1e-20,provenance=p),
    "factorization",class="cdid_gmm_pd_error")
  strict <- d26_accuracy(); strict$method <- "certified_eigen_chol_kkt"
  expect_error(solve_gmm_covariance_handler(V,c(3,2),matrix(c(1,1),2),
    list(n_links=2,rank_k=1),list(estimable_indices=1L,C_gauge=matrix(1)),g,b,s,strict,1e7),
    "certification",class="cdid_gmm_pd_error")
})

test_that("D26 validates inputs, deferral, finite controls and adopted defaults strictly", {
  expect_error(resolve_gmm_covariance_choice("moore_penrose"),"deferred",class="cdid_gmm_control_error")
  for (x in list(NULL,NA_character_,c("tikhonov","x"),"tik",structure("tikhonov",class="x")))
    expect_error(resolve_gmm_covariance_choice(x),class="cdid_gmm_control_error")
  for (x in list(list(mp_cutoff=.1),list(moore_penrose_cutoff=.1),list(unknown=1),
    list(ridge_relative=0),list(ridge_relative=Inf),list(condition_trigger=1),
    list(condition_trigger=NA_real_),list(condition_trigger=c(.1,.2))))
    expect_error(resolve_gmm_covariance_tuning(x),class="cdid_gmm_control_error")
  for (x in list(NULL,NA_character_,c("identity","optimal"),factor("optimal")))
    expect_error(resolve_gmm_weighting_choice(x),class="cdid_gmm_control_error")
  alias <- d26_core(weighting="two_step")
  expect_equal(alias$targets$estimate,d26_core()$targets$estimate,tolerance=1e-14)
  expect_identical(alias$covariance_handling$weighting_requested,"two_step")
  for (V in list(matrix(Inf,2,2),matrix(NaN,2,2),matrix(c(1,0,1,1),2),diag(c(1,-1))))
    expect_error(d26_core(V),class="cdid_gmm_error")
  expect_error(d26_core(Phi=matrix(Inf,4,2)),class="cdid_gmm_replay_error")
  expect_error(d26_handler(matrix(c(1,1),2),c(1e200,-1e200),diag(c(1,0)),alpha=2),
    "quadratic residual overflowed",class="cdid_gmm_accuracy_error")
  expect_error(d26_core(reduction=list(mode="sample_rank")),class="cdid_gmm_reduction_error")
  expect_error(d26_core(weighting="identity",reduction=list(mode="none",reduction_active=TRUE)),class="cdid_gmm_reduction_error")
  a <- chained_public_resolve_control(NULL)
  e <- chained_public_resolve_control(list(numeric_workspace_bytes=1073741824,counted_work_limit=100000000000,workers=1))
  expect_identical(a$profile$counted_work_limit,Inf)
  expect_identical(e$profile$counted_work_limit,100000000000)
  expect_identical(a$profile[c("numeric_workspace_bytes","workers")],e$profile[c("numeric_workspace_bytes","workers")])
  expect_identical(a$covariance_strategy,e$covariance_strategy)
  expect_identical(a$profile_status,"ADOPTED_D78_OPTIONAL_WORK_CAP_D27_MEMORY_DEFAULTS")
  for (v in list(list(workers=0),list(counted_work_limit=-Inf),list(numeric_workspace_bytes=7),
    list(mp_cutoff=.1),list(ridge_relative=0),list(chunk_rows=matrix(1)),
    setNames(list(1),NA_character_))) expect_error(chained_public_resolve_control(v),class="cdid_chained_public_control_error")
  dims <- chained_public_dimensions_from_counts(4,4,8,2,2,1,1,2,4)
  p <- chained_public_pipeline_control(dims,a)
  pe <- chained_public_pipeline_control(dims,e)
  expect_identical(p$control,chained_public_pipeline_control(dims,chained_public_resolve_control(list(counted_work_limit=Inf)))$control)
  expect_identical(pe$control$owner_control$max_work_ops,100000000000)
  pi <- chained_public_pipeline_control(dims,chained_public_resolve_control(NULL,"identity"))
  expect_false(pi$control$gmm_accuracy$require_matching_scale)
  expect_equal(d26_core(weighting="identity",accuracy=pi$control$gmm_accuracy)$targets$estimate,2.5,tolerance=1e-14)
})

test_that("D26 allocation and counted-work boundaries include simultaneous buffers", {
  bounds <- d26_bounds(); s <- gmm_covariance_strategy()
  pre <- function(b=bounds,strategy=s) preflight_gmm_allocations(2,1,1,"optimal",b,
    n_out_rows=4,n_input_rows=4,n_parameters=1,request_param_cov=TRUE,
    materialize_target_if=TRUE,certify_error=FALSE,covariance_strategy=strategy,n_dropped_triplets=3)
  p <- pre(); old <- pre(strategy=NULL)
  # Independent enumeration of the specified live-buffer reservations, including
  # copied/transient matrices. Literal dimensions 2,1,1 and three triplets.
  expected <- c(2*4,4*4,12*3,4*4,8*4,40*2,6*4,6*4,6*2,8*2,6,4,16*4,4)
  names(expected) <- c("D26_reconstruction_G_C","D26_reconstruction_scratch","D26_sparse_triplets_indices",
    "D26_transformed_covariance_products","D26_eigensystem_matrix_workspace","D26_eigenvalues_linear_workspace",
    "D26_eigen_residual_norm_buffers","D26_ridge_factor_residual_buffers","D26_design_solve_buffers",
    "D26_loading_solve_buffers","D26_coefficient_solve_buffers","D26_design_gram_workspace",
    "D26_vector_workspace","D26_gram_provenance_snapshot")
  expect_identical(p$extra_payload_entries[names(expected)],expected)
  expect_equal(p$total_payload_entries-old$total_payload_entries,sum(expected),tolerance=0)
  expect_equal(p$covariance_strategy_work_ops,160*8+64*4*3+64*2*3+16*3*6,tolerance=0)
  b <- bounds; b$max_matrix_entries <- p$total_payload_entries
  expect_identical(pre(b)$total_payload_entries,p$total_payload_entries)
  b$max_matrix_entries <- b$max_matrix_entries-1
  expect_error(pre(b),class="cdid_gmm_bound_error")
  work <- p$covariance_strategy_work_ops
  expect_equal(ipw_owner_safe_sum(c(17,work),17+work,"D26 boundary"),17+work,tolerance=0)
  expect_error(ipw_owner_safe_sum(c(17,work),16+work,"D26 boundary"),class="cdid_ipw_owner_bound_error")
  expect_error(assemble_gmm_gram_matrix(2,"optimal",NULL,23),class="cdid_gmm_bound_error")
  expect_error(diagnose_gmm_transformed_covariance(diag(2),assemble_gmm_gram_matrix(2,"optimal"),
    d26_provenance(),s,151),class="cdid_gmm_bound_error")
  expect_error(preflight_gmm_allocations(2^30,1,1,"optimal",bounds,covariance_strategy=s),class="cdid_gmm_bound_error")
})

test_that("D26 automatic reconstructed ridge and zero PSD preserve the full catalog", {
  red <- list(mode="sparse_class_certificate",n_links=3,retained_index=c(1L,2L),
    dropped_index=3L,kept_indices=c(2L,1L),
    E=list(i=c(3L,3L),j=c(1L,2L),x=c(-2,1)))
  E <- rbind(c(0,1),c(1,0),c(1,-2)); H <- matrix(c(1,1),2); m <- c(3,2)
  V <- diag(c(1,0)); Vf <- E%*%V%*%t(E); lambda <- 1e-4*sum(diag(Vf))/3
  f <- d26_core(reduction=red)
  r <- d26_reference(E%*%H,E%*%m,Vf,solve(Vf+diag(lambda,3)))
  expect_equal(f$targets$estimate,r$point,tolerance=1e-14)
  expect_equal(t(unname(f$internal_gauge$Lambda)),unname(r$A%*%E),tolerance=1e-14)
  expect_equal(unname(f$sandwich_cov),unname(r$sandwich),tolerance=1e-14)
  expect_equal(f$covariance_handling$lambda,lambda,tolerance=1e-14)
  expect_identical(f$covariance_handling$dimensions$M_full,3L)
  expect_identical(f$covariance_handling$kept_indices,c(2L,1L))
  z <- d26_core(V=matrix(0,2,2),Phi=matrix(0,4,2),reduction=red)
  id <- d26_reference(E%*%H,E%*%m,matrix(0,3,3),diag(3))
  expect_equal(z$targets$estimate,id$point,tolerance=1e-14)
  expect_equal(t(unname(z$internal_gauge$Lambda)),unname(id$A%*%E),tolerance=1e-14)
  expect_equal(as.double(z$sandwich_cov),0,tolerance=1e-14)
  # Structural reduction alone does not activate the handler on remaining PD.
  pd <- d26_core(V=diag(2),reduction=red)
  expect_false(pd$covariance_handling$handler_invoked)
})

test_that("D26 ridge keeps D5 anchored levels, unanchored contrasts and unavailable levels", {
  des <- build_chained_design(data.frame(g=c(2,3),s=c(1,4),t=c(2,5),ref=c(1,1),
    functional_key=c("a","b")),target_post_times=data.frame(g=2,time=6))
  pt <- des$parameters
  idx <- function(g,t) pt$param_index[pt$g==g & pt$time==t]
  targets <- list(gmm_target_spec("level",idx(2,2),1),
    gmm_target_spec("contrast",c(idx(3,4),idx(3,5)),c(-1,1)),
    gmm_target_spec("unanchored4",idx(3,4),1),gmm_target_spec("unanchored5",idx(3,5),1),
    gmm_target_spec("isolated",idx(2,6),1))
  x <- c(1,1,-1,-1); Phi <- cbind(x,2*x); V <- unname(crossprod(Phi)/4)
  p <- d26_provenance(V=V,internal=TRUE)
  replay <- gmm_if_replay_spec(function(consume) consume(list(unit_id=1:4,phi=Phi,
    nuisance_fit_id="D26:fixed")),"unit",FALSE,4,4,4,TRUE,FALSE,"unweighted",.25,
    "D26:fixed","D26:fixed",1:4)
  run <- function(w) solve_chained_gmm_core(des,c(3,7),V,p,targets,w,replay,d26_accuracy(w),
    d26_bounds(),"lowest",materialize_target_if=TRUE)
  r <- run("optimal"); id <- run("identity")
  expect_identical(r$targets$estimable,c(TRUE,TRUE,FALSE,FALSE,FALSE))
  expect_identical(r$targets$estimable,id$targets$estimable)
  expect_identical(r$targets$status,id$targets$status)
  expect_identical(r$targets$support_reason,id$targets$support_reason)
  expect_equal(r$targets$estimate[1:2],c(3,7),tolerance=1e-12)
  expect_true(all(is.na(r$targets$estimate[3:5])))
  expect_true(all(is.na(r$parameters$level_estimate[!r$parameters$level_estimable])))
  expect_identical(ncol(r$if_projection$materialized_target_if),2L)
  expect_equal(unname(r$if_projection$materialized_target_if),unname(Phi),tolerance=1e-12)
  expect_length(r$witnesses,5)
})

test_that("D26 complete IFs retain aggregation share terms and the full pretrend-null rule", {
  # Fixed target IF stack, as from the preceding exactly identified system.
  x <- c(1,1,-1,-1); X <- cbind(x,2*x)
  meta <- data.frame(target_id=c("a","b"),g=c(2,3),time=c(2,3),
    role="post",available=TRUE,support_status="estimable")
  dims <- chained_public_dimensions_from_counts(4,4,8,3,2,2,1,2,4)
  budget <- chained_public_pipeline_control(dims,chained_public_resolve_control(NULL))$control$multiplier_budget
  for (type in c("simple","group","dynamic","calendar")) {
    a <- aggregate_chained_ipw_targets(c(3,7),X,meta,c(2,2,3,3),rep(1,4),type=type,
      multiplier_budget=budget,inference_accuracy=NULL,bstrap=FALSE,cband=FALSE)
    # Simple, group and event-time outer weights depend on cohort shares.
    # Calendar averages the two single-cohort date bins with fixed outer weights.
    expected <- if (type=="calendar") 1.5*x else -.5*x
    expect_equal(a$overall$estimate,5,tolerance=1e-14)
    expect_equal(as.double(a$overall$if_full),expected,tolerance=1e-14)
    expect_equal(a$overall$se,sqrt(sum(expected^2)/16),tolerance=1e-14)
    expect_null(a$numerical_certificate)
  }
  pmeta <- data.frame(target_id=c("p1","ref","p2"),role=c("pre","reference","pre"),
    available=TRUE,support_status="supported")
  report <- list(target_meta=pmeta,estimates=c(.3,0,-.2),vcov=diag(c(.04,0,.09)))
  full <- chained_pretest_record(report)
  expect_identical(full$status,"computed"); expect_identical(full$df,2L)
  expect_equal(full$W,.3^2/.04+.2^2/.09,tolerance=1e-14)
  report$target_meta$available[3] <- FALSE
  missing <- chained_pretest_record(report)
  expect_identical(missing$status,"unavailable")
  expect_identical(missing$cells,full$cells); expect_null(missing$W)
  report$target_meta <- pmeta; report$vcov[c(1,3),c(1,3)] <- 1
  singular <- chained_pretest_record(report)
  expect_identical(singular$status,"unavailable")
  expect_identical(singular$cells,full$cells); expect_null(singular$W)
})

# Fixed 8x8 witness of differing LAPACK eigenvalue paths at the shared trigger.
# On platforms routing below the boundary, verify the declared ridge instead.
test_that("D28 ordinary routing is not contradicted by a second eigenvalue path", {
  L <- 8L; s <- 1.3
  Q <- qr.Q(qr(outer(seq_len(L),seq_len(L),function(i,j) sin(s*i*j+i)+cos(j/(s+1)))))
  d <- c(1, exp(seq(log(.9),log(1e-3),length.out=L-2L)),1e-10*(1-1e-7))
  V <- Q %*% (d*t(Q)); V <- (V+t(V))/2
  des <- build_chained_design(data.frame(g=rep(2,L),s=1,t=2,ref=1,
                                        functional_key=paste0("boundary",seq_len(L))))
  dimnames(V) <- list(des$incidence$link_ids,des$incidence$link_ids)
  n <- L+5L; Phi <- matrix(0,n,L,dimnames=list(NULL,des$incidence$link_ids))
  p <- gmm_covariance_provenance("unit",n,n,n,TRUE,FALSE,"unweighted",1/n,
                                "boundary","boundary",list(mode="none"))
  replay <- gmm_if_replay_spec(function(consume) consume(list(unit_id=seq_len(n),
    phi=Phi,nuisance_fit_id="boundary")),"unit",FALSE,n,n,n,TRUE,FALSE,"unweighted",1/n,
    "boundary","boundary",seq_len(n))
  acc <- gmm_accuracy_control(method="conventional_chol_qr",target_error_unit="abs",
    max_target_error_abs=.Machine$double.xmax,max_target_error_se=NA_real_,
    max_orthogonality_defect=1e-8,min_certified_lambda_rel=1e-10,max_factor_residual_rel=1e-10,
    min_rcond_R=1e-10,min_rcond_U=1e-10,min_sigma_rel_A=1e-10,max_kkt_residual_rel=1e-10,
    max_loading_residual_rel=1e-10,require_matching_scale=FALSE,
    covariance_strategy=gmm_covariance_strategy())
  gram <- assemble_gmm_gram_matrix(L,"optimal",list(mode="none"),1e7)
  diagnosis <- diagnose_gmm_transformed_covariance(V,gram,p,acc$covariance_strategy,1e7)
  m <- seq_len(L)/L
  fit <- solve_chained_gmm_core(des,m,V,p,list(gmm_target_spec("theta",1L,1)),"optimal",replay,
    acc,gmm_bounds_control(1e7,n,20,n,TRUE,1e5),"lowest")
  expect_identical(fit$covariance_handling$handler_invoked,diagnosis$triggered)
  lambda <- if(diagnosis$triggered) 1e-4*sum(diag(V))/L else 0
  q1 <- solve(V+diag(lambda,L),rep(1,L)); ref <- sum(q1*m)/sum(q1)
  expect_equal(unname(fit$targets$estimate),unname(ref),tolerance=1e-10)
  expect_true(all(fit$targets$estimable))
})
