# D23 conventional public contract. Deterministic unit fixtures only. No DGP,
# bootstrap, benchmark or saved-input whole fit. Saved N200 fields test behavior
# parity; the unchanged Newton and rotating80 references test formulas separately.
d23_file <- function(kind, name) {
  testthat::test_path(kind, name)
}
d23_controls <- function(n = 200, p = 2) {
  dims <- chained_public_dimensions_from_counts(n, n, 5*n, 5, 4, 3, p, 32, n)
  chained_public_pipeline_control(dims, chained_public_resolve_control(NULL))
}
# Scoped clones preserve production bodies; substitutes test only composition
# seams without masquerading as an installed/native execution.
d23_clone <- function(overrides = list()) {
  ns <- environment(chained_public_att_gt); e <- new.env(parent = ns)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, ns)
    if (is.function(f) && identical(environment(f), ns)) {
      environment(f) <- e; assign(nm, f, e)
    }
  }
  list2env(overrides, e); e
}

test_that("D23 public controls defer accuracy and retain adopted D27 resource overrides", {
  out <- d23_controls(); policy <- out$record$numerical_policy
  expect_identical(policy$decision, "D23")
  expect_false(policy$complete_per_fit_certification_requested)
  expect_false(policy$complete_per_fit_certification_performed)
  expect_identical(policy$relative_accuracy_status, "deferred_by_D23")
  expect_null(out$control$inference_accuracy)
  expect_identical(out$control$gmm_accuracy$method, "conventional_chol_qr")
  expect_identical(out$record$profile_status, "ADOPTED_D78_OPTIONAL_WORK_CAP_D27_MEMORY_DEFAULTS")
  for (v in list(1e-6, 1e-2, NULL, NA_real_))
    expect_error(chained_public_resolve_control(list(relative_accuracy = v)), class = "cdid_chained_public_accuracy_deferred")
  r <- chained_public_resolve_control(list(chunk_rows = 7, numeric_workspace_bytes = 2^30, counted_work_limit = 2e9))
  expect_identical(r$chunk_rows, 7); expect_identical(r$profile$counted_work_limit, 2e9)
  expect_identical(r$profile$numeric_workspace_bytes, 2^30)
})

test_that("D23 saved-output parity preserves N200 analytical inference and every aggregation map without enclosures", {
  z <- readRDS(d23_file("fixtures", "cp02_result_adapter_N200.rds")); r <- z$reporting; fr <- z$frame
  budget <- d23_controls()$control$multiplier_budget
  inf <- ipw_pipeline_supported_inference(r$estimates, r$if_full, r$available,
    multiplier_budget = budget, inference_accuracy = NULL, bstrap = FALSE, caller_rng = TRUE)
  for (f in c("estimates", "se", "analytical_se", "vcov", "ci_lower", "ci_upper", "cband_lower", "cband_upper"))
    expect_equal(unname(inf[[f]]), unname(r$inference[[f]]), tolerance = 1e-14, info = f)
  expect_null(inf$numerical_certificate)
  expect_match(inf$band_family, "supported_ids")
  G <- replace(fr$G, fr$G == 0, Inf)
  for (type in c("simple", "group", "dynamic", "calendar")) {
    a <- aggregate_chained_ipw_targets(r$estimates, r$if_full, r$target_meta, G, fr$sample_weight,
      type = type, available = r$available, multiplier_budget = budget, inference_accuracy = NULL,
      bstrap = FALSE, caller_rng = TRUE, calendar = z$calendar, cohort_grid = sort(unique(G)))
    old <- z$aggregations[[type]]
    for (f in c("estimates", "se", "if_full", "weights", "share_gradient"))
      expect_equal(unname(a[[f]]), unname(old[[f]]), tolerance = 1e-14, info = paste(type, f))
    expect_equal(a$overall$estimate, old$overall$estimate, tolerance = 1e-14)
    expect_equal(a$overall$se, old$overall$se, tolerance = 1e-14)
    expect_null(a$numerical_certificate); expect_null(a$overall$numerical_certificate)
    expect_null(a$weight_enclosure); expect_null(a$if_shares_enclosure)
  }
  bad_if <- r$if_full; bad_if[1,1] <- Inf
  expect_error(ipw_pipeline_supported_inference(r$estimates,bad_if,r$available,
    multiplier_budget=budget,inference_accuracy=NULL,bstrap=FALSE),class="cdid_ipw_pipeline_integrity_error")
  # Independent cluster-sum covariance formula (distinct from saved parity).
  cl <- fr$cluster; S <- rowsum(r$if_full, cl, reorder = FALSE)
  S <- sweep(S, 2, colMeans(S), "-"); V <- crossprod(S)/nrow(fr)^2
  ci <- ipw_pipeline_supported_inference(r$estimates, r$if_full, r$available, cl,
    multiplier_budget = budget, inference_accuracy = NULL, bstrap = FALSE)
  expect_equal(unname(ci$vcov), unname(V), tolerance = 1e-14)
  expect_identical(ci$estimates, inf$estimates) # D11: clustering changes inference only.
  # FALSE/TRUE cells retain analytic inference and consume no RNG draws.
  seed_before <- get0(".Random.seed", globalenv(), inherits = FALSE)
  ft <- ipw_pipeline_supported_inference(r$estimates, r$if_full, r$available, cband = TRUE,
    multiplier_budget = budget, inference_accuracy = NULL, bstrap = FALSE, caller_rng = TRUE)
  expect_identical(ft$se, inf$se); expect_null(ft$multiplier)
  expect_identical(get0(".Random.seed", globalenv(), inherits = FALSE), seed_before)
})

test_that("D23 D21 computes the full declared null and preserves unavailable restrictions", {
  meta <- data.frame(target_id = c("p1", "ref", "p2", "post"), role = c("pre", "reference", "pre", "post"),
    available = TRUE, support_status = "supported")
  V <- diag(c(.04, 0, .09, .1)); V[1,3] <- V[3,1] <- .01
  r <- list(base_period = "universal", target_meta = meta, estimates = c(.3, 0, -.2, 2), vcov = V)
  a <- chained_pretest_record(r)
  expect_identical(a$status, "computed"); expect_identical(a$df, 2L)
  expect_identical(a$cells, c("p1", "p2")); expect_identical(a$structural_reference_cells, "ref")
  expect_equal(a$W, .0109/.0035, tolerance = 1e-13)
  expect_equal(a$Wpval, exp(-a$W/2), tolerance = 1e-13)
  expect_false(a$certified); expect_null(a$W_interval); expect_null(a$arithmetic$error_bound)
  before <- r; r$target_meta$available[3] <- FALSE
  b <- chained_pretest_record(r)
  expect_identical(b$status, "unavailable"); expect_identical(b$df, 2L)
  expect_identical(b$cells, a$cells); expect_null(b$W)
  expect_identical(r$estimates, before$estimates)
  r <- before; r$vcov[c(1,3),c(1,3)] <- 1
  expect_identical(chained_pretest_record(r)$status, "unavailable")
  r <- before; r$vcov[1,1] <- Inf
  expect_identical(chained_pretest_record(r)$status, "unavailable")
  r <- before; r$vcov <- NULL
  expect_identical(chained_pretest_record(r)$status, "unavailable")
  r <- before; r$vcov <- diag(2)
  expect_identical(chained_pretest_record(r)$status, "unavailable")
  r <- before; r$target_meta$role[c(1,3)] <- "post"
  expect_identical(chained_pretest_record(r)$status, "not_applicable")
})

test_that("D23 zero-dual fixture remains accurate despite its large stationarity ratio", {
  tab <- read.delim(d23_file("reference", "zero_dual_stored_binary64.tsv"), stringsAsFactors = FALSE)
  pull <- function(nm) {
    t <- tab[tab$name == nm, ]; v <- as.numeric(t$hex)
    if (nrow(t) && nzchar(t$dim[1])) matrix(v, 3, 3) else v
  }
  args <- lapply(c("V_R", "R_chol", "W_G", "A_whitened", "U_qr", "piv", "mhat_R", "theta", "y", "C_gauge", "Lambda", "Mc", "Z1"), pull)
  names(args) <- c("V_R", "R_chol", "W_G", "A_whitened", "U_qr", "piv", "mhat_R", "theta", "y", "C_gauge", "Lambda", "Mc", "Z1")
  acc <- d23_controls()$control$gmm_accuracy
  # Deliberately tiny unused certificate budgets: no conservative-bound refusal.
  acc$max_target_error_abs <- acc$max_loading_error_abs <- 1e-300
  run <- function(a = args, w = "optimal") do.call(certify_gmm_numerical_guard, c(a, list(weighting = w, accuracy = acc)))
  g <- run()
  expect_equal(g$computed_target_values, c(-2,2,5), tolerance = 1e-14)
  expect_gt(g$primal_kkt$stat_rel_est, .5)
  expect_false(g$certified); expect_false(g$point_certified)
  expect_true(all(is.na(g$targets$target_error_cert)))
  expect_true(all(is.finite(g$targets$se_criterion_est)))
  expect_null(g$propagation_bounds)
  id_args <- args; id_args$A_whitened <- id_args$W_G
  expect_true(all(is.na(run(id_args, "identity")$targets$se_criterion_est)))
  bad <- args; bad$theta[1] <- Inf
  expect_error(run(bad), class = "cdid_gmm_accuracy_error")
  bad <- args; bad$R_chol[1,1] <- 0
  expect_error(run(bad), class = "cdid_gmm_pd_error")
  bad <- args; bad$V_R[1,1] <- -1
  expect_error(run(bad), class = "cdid_gmm_pd_error")
  bad <- args; bad$U_qr[3,3] <- 0
  expect_error(run(bad), class = "cdid_gmm_qr_error")
  # No interval guard buffer reservation for work that did not run.
  bounds <- d23_controls()$control$gmm_bounds
  strict <- preflight_gmm_allocations(3,3,3,"optimal",bounds)
  plain <- preflight_gmm_allocations(3,3,3,"optimal",bounds,certify_error = FALSE)
  expect_identical(unname(plain$extra_payload_entries["validated_guard_buffers"]), 0)
  expect_lt(plain$total_payload_entries, strict$total_payload_entries)
})

test_that("D23 exact P1 obstruction survives a converged glm flag, with finite counterexamples", {
  # Exact supplied P1 rows (checkdata.csv SHA256 a245bcd0...ba21c4):
  # 9 treated x=0, 9 treated x=1, 15 controls x=0; positive weights 53/6,11.
  D <- c(rep(1,18), rep(0,15)); Z <- cbind(1, c(rep(0,9), rep(1,9), rep(0,15)))
  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 50L)
  old <- ipw_glm_fit_internal(Z, D, c(rep(53/6,18),rep(11,15)), ctrl)$fit
  expect_true(old$converged); expect_false(old$boundary)
  expect_true(all(is.finite(old$coefficients)))
  run <- function(z = Z, gi = ifelse(D == 1,1L,2L), control = ctrl, de = NULL) ipw_nuisance_observer_fit_worker(
    list(n = 639, m_cohorts = 2L, group_mass = c(159,165)), seq_len(nrow(z)), gi, rep(1,nrow(z)), z, control, diag_env = de)
  de <- new.env(parent = emptyenv())
  e <- expect_error(run(de = de), class = "cdid_ipw_separation_error")
  expect_match(conditionMessage(e), "no finite root")
  expect_identical(e$diagnostics$separation_witness$coordinate, 2L)
  expect_identical(e$diagnostics$separation_witness$strict_rows, 10:18)
  expect_identical(length(e$diagnostics$separation_witness$zero_rows), 24L)
  expect_identical(de$error_diagnostics$stage, "coordinate_separation_check")
  expect_identical(de$error_diagnostics$separation_witness, e$diagnostics$separation_witness)
  # One control at x=1 gives both outcomes at both x values: finite root,
  # including after a negative scale (test +/- signed coordinate directions).
  finite <- Z; finite[33,2] <- 1
  expect_null(ipw_logit_coordinate_separation(finite,D))
  expect_true(run(finite)$fit_info$converged)
  # Zero-weight rows neither mask positive-weight separation nor trigger false positives.
  expect_identical(ipw_logit_coordinate_separation(finite, D, c(rep(1, 32), 0))$coordinate, 2L)
  expect_null(ipw_logit_coordinate_separation(Z, D, c(rep(1, 9), rep(0, 9), rep(1, 15))))
  finite[,2] <- -finite[,2]
  expect_null(ipw_logit_coordinate_separation(finite,D))
  expect_true(run(finite)$fit_info$converged)
  negative <- Z; negative[,2] <- -negative[,2]
  expect_identical(ipw_logit_coordinate_separation(negative,D)$sign, -1)
  allzero <- matrix(0,33,1)
  expect_null(ipw_logit_coordinate_separation(allzero,D)) # no strict product
  expect_error(run(finite, control = stats::glm.control(epsilon = 1e-15, maxit = 1)), class = "cdid_ipw_convergence_error")
  # Hessian curvature still refuses unusable fits independently of transfer.
  expect_error(ipw_nuisance_hessian_stage(list(Z_use = cbind(1,1), fit_weights = 1, prob_use = .5), list(n=1)),
    class = "cdid_ipw_rank_error")
})

test_that("D23 finite dense point and complete IF match the immutable kernel plus ONE matched-start adaptation and unchanged independent reference", {
  # Existing deterministic shared-audit fixture; no random draws.
  n <- 96L; G <- rep(c(4,7,Inf,9),each=24L); x <- rep(seq(-1.4,1.4,length.out=24),4)
  Z <- cbind(one=1,x=x); w <- .7 + ((seq_len(n)*13L)%%29L)/20
  obs <- matrix(FALSE,n,8)
  for(i in seq_len(n)) {
    k <- ((i-1L)%%24L)+1L
    if(G[i] != 9 && k<=18L) obs[i,1:6] <- TRUE else {obs[i,7:8] <- TRUE; obs[i,1L+(k%%6L)] <- TRUE}
  }
  Y <- outer(seq_len(n)/17,1:8,function(a,b) sin(a*b)+a+b^2/13) + outer(x,1:8,function(a,b) a*b/3)
  panel <- data.frame(id=rep(seq_len(n),each=8),time=rep(as.double(1:8),n),y=as.vector(t(ifelse(obs,Y,NA_real_))),g=rep(G,each=8))
  ctx <- construct_ipw_context(construct_observed_pairs(panel,"id","time","y","g"),seq_len(n),Z,weights=w)
  ctrl <- stats::glm.control(epsilon=1e-12,maxit=100)
  f <- fit_ipw_link(ctx,4,c(7,Inf),1,2,ctrl)
  X <- materialize_ipw_if(ctx,f,max_entries=1000,chunk_size=13,include_components=TRUE)$phi_matrix
  old <- new.env(parent=environment(fit_ipw_link))
  ipw_reference_source(d23_file("reference","chained_ipw_link_original_6550dc.R"),old)
  # The frozen kernel accepts three fields; native-only fields must be defaults.
  legacy3 <- ctrl[c("epsilon", "maxit", "trace")]
  expect_identical(do.call(stats::glm.control, legacy3), ctrl)
  fo <- old$fit_ipw_link(ctx,4,c(7,Inf),1,2,legacy3)
  expect_identical(f$estimate,fo$estimate)
  expect_equal(X,old$materialize_ipw_if(ctx,fo,max_entries=1000,chunk_size=13)$phi_matrix,tolerance=1e-13)
  ref <- new.env(parent=baseenv())
  sys.source(d23_file("reference","cp02_shared_ipw_audit_newton_reference.R"),ref)
  expected <- ref$independent_fit(list(n=n,G=G,Z=Z,w=w,obs=obs,Y=Y),1,2)
  expect_equal(f$estimate,expected$estimate,tolerance=1e-13)
  expect_equal(as.double(X),expected$phi,tolerance=1e-13)
  expect_equal(as.double(stream_ipw_covariance(ctx,f,13)$covariance),var(expected$phi)*(n-1)/n,tolerance=1e-13)
})

test_that("D23 public execution and result metadata request no transfer and aggte retains that policy", {
  z <- readRDS(d23_file("fixtures","cp02_result_adapter_N200.rds")); r <- z$reporting
  calls <- new.env(); calls$n <- 0L
  spy <- d23_clone(list(
    run_chained_ipw_pipeline = function(...) {
      a <- list(...); calls$n <- calls$n+1L; calls$args <- a
      expect_false(a$run_numerical_transfer); expect_null(a$root_witness)
      expect_null(a$control$inference_accuracy)
      structure(list(reporting=r, panel_owner=list(a_id=z$frame$id, n_a_units=200L,all_times=z$calendar,
        cohorts=c(3,4,5,Inf)), numerical_transfer=NULL, support=z$support, reduction_summary=z$reduction_summary),class="cdid_ipw_pipeline_result")
    },
    chained_did_AGGTEobj = function(aggregation, ...) list(aggregation=aggregation,chained=list()),
    chained_did_MP = function(reporting, full_id, call, params, units, ...) {
      structure(list(att=reporting$estimates,inffunc=reporting$if_full,DIDparams=params,
        chained=list(frame=units,available=reporting$available,target_id=reporting$target_meta$target_id,target_meta=reporting$target_meta)),class=c("cdid_MP","MP"))
    }))
  # Request-composition spy only: saved reporting values pass through unchanged.
  d <- expand.grid(id=z$frame$id,period=z$calendar); d$G <- z$frame$G[d$id]; d$x <- z$frame$x[d$id]; d$y <- d$id+d$period
  mp <- spy$att_gt("y","period","id","G",~x,d,est_method="ipw",bstrap=FALSE,cband=FALSE)
  expect_identical(calls$n,1L)
  expect_identical(mp$att,r$estimates); expect_identical(mp$inffunc,r$if_full)
  expect_false(mp$chained$numerical_transfer$requested); expect_false(mp$chained$numerical_transfer$performed)
  expect_false(mp$chained$numerical_transfer$certified)
  expect_null(mp$chained$numerical_policy$reported_cells_d20)
  expect_null(mp$chained$aggregation_state$enclosure)
  expect_null(mp$chained$aggregation_state$inference_accuracy)
  expect_identical(mp$chained$pretest$status,"computed"); expect_false(mp$chained$pretest$certified)
  expect_false(mp$chained$request$run_numerical_transfer)
  agg <- spy$chained_aggte_request(mp,type="group")
  expect_identical(calls$n,1L) # actual aggregation consumer, no refit
  expect_null(agg$aggregation$numerical_certificate)
  expect_identical(agg$chained$numerical_policy,mp$chained$numerical_policy)
  expect_identical(agg$chained$numerical_transfer,mp$chained$numerical_transfer)
  expect_equal(agg$aggregation$if_full,z$aggregations$group$if_full,tolerance=1e-14)
  # Keep this source-only test focused on policy metadata; the upstream
  # estimate table remains covered by the installed presentation tests.
  spy$chained_result_upstream_summary <- function(...) character(0L)
  spy$chained_result_banner <- function(...) invisible(NULL)
  shown <- mp
  shown$chained$inference <- list(bstrap=FALSE,cband=FALSE,clustered=FALSE,n_total=200L)
  shown$chained$pretest <- mp$chained$pretest
  shown$n <- shown$chained$inference$n_total # the source-only constructor spy omits MP metadata
  txt <- paste(capture.output(spy$summary.cdid_MP(shown)),collapse="\n")
  expect_no_match(txt,"Numerical policy:",fixed=TRUE)
  expect_identical(shown$chained$public$numerical_policy$decision,"D23")
  expect_no_match(txt,"Numerical acceptance (D20)",fixed=TRUE)
  expect_no_match(txt,"certified enclosure of W",fixed=TRUE)

})

test_that("D23 ordinary delegation stays lazy and the final dr default stays intact", {
  seen <- NULL
  spy <- d23_clone(list(forward_call=function(fun,formal_names,env,override=list(),dots=TRUE) {
    seen <<- list(fun=fun,names=formal_names); "delegated"
  }))
  expect_identical(spy$att_gt("y","t","id","g",data=stop("forced data"),mode="ordinary",control=stop("forced control")),"delegated")
  expect_identical(seen$fun,quote(did::att_gt)); expect_identical(seen$names,upstream_att_gt_formals)
  expect_identical(formals(att_gt)$est_method,"dr")
  expect_error(att_gt("y","t","id","g",data=stop("forced data")),"forced data")
  expect_error(att_gt("y","t","id","g",data=stop("forced data"),est_method="stdipw"),class="cdid_chained_unavailable")
})

test_that("D23 installed public rotating80 pipeline preserves independent point and full IF formulas", {
  skip_if_not_installed("cdid"); skip_if_not_installed("did")
  z <- source(d23_file("fixtures","rotating80_fixture.R"))$value; d <- z$panel
  d$g[is.infinite(d$g)] <- 0; d$y <- d$y_exactfit
  mp <- att_gt("y","t","id","g",~x,d,est_method="ipw",bstrap=FALSE,cband=FALSE)
  post <- which(mp$chained$target_meta$role=="post")
  expect_equal(unname(mp$att[post]),c(1,1,1),tolerance=1e-12)
  # Independent reference targets order (pre, post3, post4, post5).
  ref_if <- z$phi_ref_mat %*% t(z$C_opt_ref)
  expect_equal(unname(mp$inffunc),unname(ref_if),tolerance=1e-12)
  expect_null(mp$chained$aggregation_state$enclosure)
  expect_false(mp$chained$numerical_transfer$requested)
  expect_identical(mp$chained$resource_accounting$d5_arrays$ledger_enclosure_entries,0)
  expect_identical(mp$chained$resource_accounting$gmm_consumer$extra_payload_entries[["validated_guard_buffers"]],0)
})

test_that("D23 conventional core handoff preserves rotating80 independent point/IF/covariance references", {
  z <- source(d23_file("fixtures","rotating80_fixture.R"))$value
  des <- build_chained_design(data.frame(g=3,s=c(1,2,2,3,4),t=c(2,3,4,4,5),ref=2,functional_key="A"))
  ids <- des$incidence$link_ids; n <- 80L
  Phi <- z$phi_ref_mat; colnames(Phi) <- ids
  targets <- lapply(seq_len(nrow(des$parameters)),function(j) gmm_target_spec(paste0("t",j),j,if(j==1) -1 else 1))
  prov <- gmm_covariance_provenance(row_level="unit",n_units=n,n_clusters=n,n_positive_contributing_rows=n,
    centered=TRUE,weighted=FALSE,weight_type="unweighted",normalization_factor=1/n,scale_id="d23:fixed",nuisance_fit_id="d23:fixed",reduction=list(mode="none"))
  replay <- gmm_if_replay_spec(callback=function(consume) consume(list(unit_id=seq_len(n),phi=Phi,weight=rep(1,n),nuisance_fit_id="d23:fixed")),
    input_row_level="unit",aggregate_to_clusters=FALSE,n_input_rows=n,n_output_rows=n,n_positive_contributing_rows=n,
    centered=TRUE,weighted=FALSE,weight_type="unweighted",normalization_factor=1/n,scale_id="d23:fixed",nuisance_fit_id="d23:fixed",expected_unit_ids=seq_len(n))
  ctrl <- d23_controls(n)$control
  run <- function(V=z$Omega_ref,m=z$m_pert_ref) solve_chained_gmm_core(des,m,V,prov,targets,"optimal",replay,
    ctrl$gmm_accuracy,ctrl$gmm_bounds,"lowest",materialize_target_if=TRUE)
  f <- run()
  expect_equal(f$targets$estimate,z$theta_pert_opt_ref,tolerance=1e-12)
  expect_equal(unname(f$if_projection$materialized_target_if),unname(Phi%*%t(z$C_opt_ref)),tolerance=1e-12)
  # Core sandwich_cov is the IF covariance Omega; reported vcov is Omega/n.
  expect_equal(unname(f$sandwich_cov),unname(z$Sigma_ref),tolerance=1e-12)
  expect_false(f$point_certified); expect_false(f$loading_certified)
  expect_true(all(is.na(f$targets$target_error_cert)))
  expect_identical(unname(f$allocations$extra_payload_entries["validated_guard_buffers"]),0)
  # Valid singular automatic outputs are checked independently in D26's fixed
  # matrix tests; a materially indefinite covariance still refuses here.
  expect_error(run(V=-diag(5)),class="cdid_gmm_pd_error")
  expect_error(run(m=c(Inf,1,1,0,0)),class="cdid_gmm_design_error")
})
