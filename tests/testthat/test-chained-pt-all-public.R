# CP03 bounded software integration, not MC qualification. No fixture RNG.
# The only bootstrap equality check replays 1000 prescribed Rademacher rows.
pt_ref <- new.env(parent = baseenv())
sys.source(testthat::test_path("reference", "cp03_public_tiny_reference.R"), pt_ref)
pt_counts <- new.env(parent = emptyenv())
pt_counts$public_requests <- pt_counts$public_completed <- pt_counts$pipeline_fits <- pt_counts$aggregations <- 0L
pt_counts$fixed_multiplier_rows <- pt_counts$nuisance_fits_in_completed_requests <- 0L
pt_fit <- function(d, ...) {
  pt_counts$public_requests <- pt_counts$public_requests + 1L
  out <- do.call(cdid::att_gt, c(list(yname="y", tname="t", idname="id", gname="g", data=d,
    est_method="ipw", bstrap=FALSE, cband=FALSE), list(...)))
  pt_counts$public_completed <- pt_counts$public_completed+1L
  pt_counts$nuisance_fits_in_completed_requests <- pt_counts$nuisance_fits_in_completed_requests+length(out$chained$nuisance_diagnostics$records)
  out
}
pt_control <- function(d, weighting="optimal", advanced=NULL, biters=1000) {
  prep <- chained_request_prepare_inputs(d, NULL, "y", "t", "id", "g", ~1)
  chained_public_pipeline_control(chained_public_dimensions(prep, biters),
    chained_public_resolve_control(advanced, gmm_weighting=weighting), gmm_weighting=weighting)$control
}
pt_pipeline <- function(d, weighting="optimal", advanced=NULL, ...) {
  pt_counts$pipeline_fits <- pt_counts$pipeline_fits + 1L
  internal <- d; internal$g[internal$g == 0] <- Inf
  out <- do.call(run_chained_ipw_pipeline, c(list(data=internal, idname="id", tname="t", yname="y", gname="g",
    control=pt_control(d, weighting, advanced), weighting=weighting, control_group="nevertreated",
    run_numerical_transfer=FALSE, bstrap=FALSE, cband=FALSE, aggregations=character()), list(...)))
  pt_counts$nuisance_fits_in_completed_requests <- pt_counts$nuisance_fits_in_completed_requests+length(out$auxiliary_reporting$fitted_views)
  out
}
pt_agg <- function(mp, type) {
  pt_counts$aggregations <- pt_counts$aggregations + 1L
  args <- list(MP=mp, type=type, bstrap=FALSE, cband=FALSE)
  # Hold post targets fixed for all four types; calendar already reports post.
  if (type != "calendar") args$min_e <- 0
  do.call(cdid::aggte, args)
}
pt_post_equal <- function(a,b) {
  ia <- which(a$chained$role == "post"); ib <- which(b$chained$role == "post")
  expect_identical(a$chained$target_id[ia], b$chained$target_id[ib])
  expect_identical(a$att[ia], b$att[ib])
  expect_identical(a$inffunc[,ia,drop=FALSE], b$inffunc[,ib,drop=FALSE])
  pt_simple_leaf_num_ok <- function(actual, expected, tol = 1e-12, entries = 12L) {
    identical(typeof(actual), "double") &&
      identical(typeof(expected), "double") &&
      identical(class(actual), class(expected)) &&
      identical(is.object(actual), is.object(expected)) &&
      identical(attributes(actual), attributes(expected)) &&
      identical(length(actual), entries) &&
      identical(length(expected), entries) &&
      all(is.finite(actual)) &&
      all(is.finite(expected)) &&
      identical(actual == 0, expected == 0) &&
      identical(1 / actual[actual == 0], 1 / expected[expected == 0]) &&
      identical(sign(actual), sign(expected)) &&
      all(abs(actual - expected) <= tol * pmax(1, abs(actual), abs(expected)))
  }
  pt_simple_inf_fn_ok <- function(ax, ay, tol = 1e-12) {
    if (!is.list(ax) || !is.list(ay) ||
        !identical(class(ax), class(ay)) ||
        !identical(is.object(ax), is.object(ay)) ||
        !identical(attributes(ax), attributes(ay)) ||
        !identical(length(ax), length(ay)) ||
        !identical(which(names(ax) == "simple.att"), which(names(ay) == "simple.att")) ||
        length(which(names(ax) == "simple.att")) != 1L) {
      return(FALSE)
    }
    ax_rest <- ax; ax_rest[["simple.att"]] <- 0
    ay_rest <- ay; ay_rest[["simple.att"]] <- 0
    identical(ax_rest, ay_rest) &&
      pt_simple_leaf_num_ok(ax[["simple.att"]], ay[["simple.att"]], tol = tol)
  }
  nyt_post_ctx <- identical(ia, c(3L, 4L, 5L, 10L)) &&
    identical(ib, c(4L, 5L, 6L, 12L)) &&
    identical(nrow(a$inffunc), 12L) &&
    identical(nrow(b$inffunc), 12L) &&
    exists("cg", parent.frame(), inherits = FALSE) &&
    identical(get("cg", parent.frame(), inherits = FALSE), "notyettreated") &&
    exists("pt", parent.frame(), inherits = FALSE) &&
    get("pt", parent.frame(), inherits = FALSE) %in% c("post", "all") &&
    (!is.character(a$chained$request$control_group) ||
      identical(a$chained$request$control_group, "notyettreated")) &&
    (!is.character(b$chained$request$control_group) ||
      identical(b$chained$request$control_group, "notyettreated")) &&
    (!is.character(a$chained$request$pt_assumption) ||
      a$chained$request$pt_assumption %in% c("post", "all")) &&
    (!is.character(b$chained$request$pt_assumption) ||
      identical(a$chained$request$pt_assumption, b$chained$request$pt_assumption)) &&
    (!is.character(a$chained$request$base_period) ||
      identical(a$chained$request$base_period, "varying")) &&
    (!is.character(b$chained$request$base_period) ||
      identical(b$chained$request$base_period, "universal"))
  # Measured NYT post/all covariance only; every other context stays exact.
  if (nyt_post_ctx) {
    expect_true(pt_simple_leaf_num_ok(a$V_analytical[ia,ia,drop=FALSE],
      b$V_analytical[ib,ib,drop=FALSE], entries = 16L))
  } else {
    expect_identical(a$V_analytical[ia,ia,drop=FALSE], b$V_analytical[ib,ib,drop=FALSE])
  }
  for (type in c("simple","group","dynamic","calendar")) {
    x <- pt_agg(a,type); y <- pt_agg(b,type)
    nyt_simple <- nyt_post_ctx && identical(type, "simple")
    if (nyt_simple) {
      expect_false(pt_simple_leaf_num_ok(
        y$chained$if_fixed_full + 2e-12 * pmax(1, abs(y$chained$if_fixed_full)),
        y$chained$if_fixed_full
      ))
      expect_false(pt_simple_leaf_num_ok(as.vector(y$chained$if_fixed_full), y$chained$if_fixed_full))
      expect_false(pt_simple_leaf_num_ok(replace(y$chained$if_fixed_full, 1L, 0), y$chained$if_fixed_full))
      expect_false(pt_simple_leaf_num_ok(-y$chained$if_fixed_full, y$chained$if_fixed_full))
      expect_false(pt_simple_inf_fn_ok(c(y$inf.function, list(extra = 1)), y$inf.function))
    }
    for (nm in c("overall.att","overall.se","att.egt","se.egt","inf.function")) {
      if (nyt_simple && identical(nm, "inf.function")) {
        expect_true(pt_simple_inf_fn_ok(x[[nm]], y[[nm]]), info=type)
      } else {
        expect_identical(x[[nm]], y[[nm]], info=type)
      }
    }
    for (nm in c("if_fixed_full","if_share_full","se_fixed_shares_only")) {
      if (nyt_simple && identical(nm, "if_fixed_full")) {
        expect_true(pt_simple_leaf_num_ok(x$chained[[nm]], y$chained[[nm]]), info=type)
      } else {
        expect_identical(x$chained[[nm]],y$chained[[nm]],info=type)
      }
    }
    expect_identical(x$chained$pt,a$chained$pt)
  }
}

# New deterministic NYT fixture, with explicit potential outcomes and truth.
pt_nyt_data <- function() {
  d <- expand.grid(id=1:12,t=1:6); d$g <- rep(c(4,6,0),each=4)[d$id]
  d$y0 <- d$id/4 + d$t/2 + sin(d$id*d$t)/5
  d$tau <- ifelse(d$g != 0 & d$t >= d$g, ifelse(d$g==4,1.5,2),0)
  d$y <- d$y0+d$tau; d
}


test_that("D33 validates explicit scalar start before data, then exact calendar membership", {
  expect_identical(names(formals(cdid::att_gt))[match("pt_assumption",names(formals(cdid::att_gt)))+1L],"pt_start")
  expect_null(formals(cdid::att_gt)$pt_start)
  for (bad in list(NULL, numeric(), c(1,2), NA_real_, NaN, Inf, -Inf, "2", TRUE, list(2), matrix(2), as.Date("2000-01-01"))) {
    expect_error(cdid::att_gt(yname="y",tname="t",idname="id",gname="g",data=stop("forced data"),
      est_method="ipw",pt_assumption="all",pt_start=bad),class="cdid_chained_request_pt_start")
  }
  expect_error(cdid::att_gt(yname="y",tname="t",idname="id",gname="g",data=stop("forced data"),
    est_method="ipw",pt_start=2),class="cdid_chained_request_pt_start")
  expect_error(cdid::att_gt(yname="y",tname="t",idname="id",gname="g",data=stop("forced data"),
    est_method="ipw",pt_assumption="all",pt_start=2),"forced data")
  d <- pt_ref$cp03_public_tiny_data()
  for (bad in c(0,2+2^-40,5)) expect_error(pt_fit(d,pt_assumption="all",pt_start=bad),class="cdid_chained_request_pt_start")
  # Development path shares both validations, before any fit/control use.
  expect_error(chained_att_gt_request(yname="y",tname="t",idname="id",gname="g",data=stop("forced data"),
    est_method="ipw",pt_assumption="all"),class="cdid_chained_request_pt_start")
})

test_that("NT raw-data points, complete IF and sandwich agree with independent tiny matrices", {
  d <- pt_ref$cp03_public_tiny_data()
  for (mode in c("identity","two_step","ridge")) {
    ref <- pt_ref$cp03_public_tiny_reference(mode)
    fit <- pt_pipeline(d,weighting=if(mode=="ridge") "optimal" else mode,
      advanced=if(mode=="ridge") list(condition_trigger=.99) else NULL,
      pt_assumption="all",pt_start=2)
    # Pair ordering is specified independently, not extracted as oracle values.
    idx <- match(paste(fit$design_all$links$s,fit$design_all$links$t),paste(ref$pairs[,1],ref$pairs[,2]))
    expect_false(anyNA(idx)); expect_length(idx,6L)
    expect_equal(unname(materialize_chained_design(fit$design_all,1000)$W),ref$W[idx,,drop=FALSE],tolerance=0)
    if (mode=="identity") {
      # Identity keeps all rows fitted; direct extraction remains valid here.
      expect_equal(unname(vapply(fit$point_links_all, `[[`, numeric(1L), "estimate")),ref$m[idx],tolerance=1e-12)
    } else {
      # Optimal/ridge use the outcome-free reduction. E-expanded moments still
      # face the independent raw-data matrix truth for EVERY original pair.
      expect_true(fit$reduction_summary$reduction_active)
      expect_equal(unname(as.double(fit$reduction_summary$expanded_m)),ref$m[idx],tolerance=1e-12)
      expect_equal(unname(vapply(fit$point_links_used, `[[`, numeric(1L), "estimate")),
        ref$m[idx][fit$reduction_summary$kept_indices],tolerance=1e-12)
      expect_identical(fit$reduction_summary$reduction_stage,"outcome_free_before_point_fits")
      expect_equal(fit$reduction_summary$n_rows_fitted,length(fit$reduction_summary$kept_indices))
      expect_identical(fit$reduction_summary$max_moment_expansion_diff,NA_real_)
      expect_identical(fit$reduction_summary$expansion_self_test,
        "not_evaluated: dropped rows are never fitted (Prop. 5/Cor. 6)")
      expect_gt(length(fit$reduction_summary$dropped_indices),0L)
      for (m in fit$reduction_summary$dropped_indices) {
        record <- fit$point_links_all[[m]]
        expect_s3_class(record,"cdid_ipw_reduced_link_record")
        expect_identical(names(record),c("link_id","catalog_row","link_view_ref","nuisance_view_ref",
          "fitted_nuisance_ref","expansion","point_fit_performed","outcomes_read","status","scope"))
        expect_null(record[["estimate",exact=TRUE]])
        expect_false(record$point_fit_performed); expect_false(record$outcomes_read)
        expect_identical(record$status,"structurally_reduced_outcome_free")
        expect_match(record$scope,"never fitted, replayed or read",fixed=TRUE)
      }
      # Fit every original row via the accepted private reference on the same
      # tiny fixture, preserving the independently fitted expansion self-test.
      reference_fit <- pt_pipeline(d,weighting=if(mode=="ridge") "optimal" else mode,
        advanced=if(mode=="ridge") list(condition_trigger=.99) else NULL,
        pt_assumption="all",pt_start=2,structural_reduction_stage="post_point_fit")
      expect_identical(reference_fit$reduction_summary$reduction_stage,"after_point_fits")
      expect_equal(reference_fit$reduction_summary$n_rows_fitted,6L)
      expect_identical(reference_fit$reduction_summary$expansion_self_test,"evaluated_on_fitted_rows")
      all_m <- unname(vapply(reference_fit$point_links_all, `[[`, numeric(1L), "estimate"))
      expect_equal(all_m,ref$m[idx],tolerance=1e-12)
      expect_equal(unname(as.double(reference_fit$reduction_summary$expanded_m)),all_m,tolerance=1e-12)
      expect_lt(reference_fit$reduction_summary$max_moment_expansion_diff,1e-12)
      for (nm in c("support","covariance_stream","full_frame_transport","reporting","basis","aggregations"))
        expect_equal(fit[[nm]],reference_fit[[nm]],tolerance=1e-10,info=nm)
      for (nm in c("targets","parameters","internal_gauge","sandwich_cov","model_cov","if_projection","numerical_guard"))
        expect_equal(fit$gmm_fit[[nm]],reference_fit$gmm_fit[[nm]],tolerance=1e-10,info=nm)
      expect_identical(fit$covariance_handling$handler_invoked,reference_fit$covariance_handling$handler_invoked)
      expect_equal(fit$covariance_handling$lambda,reference_fit$covariance_handling$lambda,tolerance=1e-10)
    }
    p <- which(fit$gmm_fit$parameters$role=="post")
    target <- match(fit$gmm_fit$parameters$param_id[p],fit$gmm_fit$targets$target_id)
    expect_length(target,1L)
    expect_equal(unname(fit$gmm_fit$targets$estimate[target]),ref$point[2],tolerance=1e-10)
    expect_equal(unname(fit$full_frame_transport$Phi_full[,target]),ref$IF[,2],tolerance=1e-10)
    expect_equal(unname(fit$gmm_fit$sandwich_cov[target,target]/8),ref$covariance[2,2],tolerance=1e-10)
    expect_gt(max(abs(ref$phi[,4])),0)
    expect_true(any(ref$V[4,-4] != 0)) # informative zero-loading comparison
    expect_equal(colSums(ref$IF),c(0,0),tolerance=1e-12)
    if (mode=="identity") {
      expect_false(fit$reduction_summary$reduction_active)
      expect_equal(unname(fit$gmm_fit$internal_gauge$Lambda[,target]),ref$Lambda[idx,2],tolerance=1e-12)
    }
    if (mode=="ridge") {
      expect_true(fit$covariance_handling$handler_invoked)
      expect_equal(fit$covariance_handling$lambda,ref$lambda,tolerance=1e-10)
      full <- pt_pipeline(d,advanced=list(condition_trigger=.99),pt_assumption="all",pt_start=2,reduce_structural_redundancy=FALSE)
      expect_equal(full$gmm_fit$targets$estimate,fit$gmm_fit$targets$estimate,tolerance=1e-10)
      expect_equal(full$full_frame_transport$Phi_full,fit$full_frame_transport$Phi_full,tolerance=1e-10)
      expect_equal(unname(full$gmm_fit$internal_gauge$Lambda[,target]),ref$Lambda[idx,2],tolerance=1e-10)
    }
  }
})

test_that("NT and fixed-horizon NYT keep post inference and four fixed-target aggregations under either base period", {
  for (cg in c("nevertreated","notyettreated")) {
    d <- if (cg=="nevertreated") pt_ref$cp03_public_tiny_data() else pt_nyt_data()
    for (pt in c("post","all")) {
      a <- pt_fit(d,control_group=cg,pt_assumption=pt,pt_start=if(pt=="all") 2 else NULL,base_period="varying")
      b <- pt_fit(d,control_group=cg,pt_assumption=pt,pt_start=if(pt=="all") 2 else NULL,base_period="universal")
      pt_post_equal(a,b)
      if (pt=="all") {
        expect_identical(a$chained$request$pt_start,2)
        expect_identical(a$chained$pt$pt_start,2)
        expect_gt(a$chained$pt$n_restricted_pre_parameters,0L)
        expect_match(a$chained$pretest$interpretation,"auxiliary NYT pool")
        # Neither base period substitutes mechanical restricted zeros for direct pre comparisons.
        for (mp in list(a,b)) {
          pre <- which(mp$chained$role=="pre" & mp$chained$available)
          expect_gt(length(pre),0L); expect_gt(max(abs(mp$att[pre])),1e-8)
          expect_gt(max(abs(mp$inffunc[,pre,drop=FALSE])),0)
        }
      }
    }
  }
  # Same source data, different controls: inspect actual NYT catalog via its public path dimensions.
  d <- pt_nyt_data(); z <- d; z$g[z$g==0] <- Inf
  pt_counts$pipeline_fits <- pt_counts$pipeline_fits+1L
  nyt <- run_chained_ipw_pipeline(z,"id","t","y","g",pt_control(d),
    control_group="notyettreated",pt_assumption="all",pt_start=2,
    run_numerical_transfer=FALSE,bstrap=FALSE,cband=FALSE,aggregations=character())
  pt_counts$nuisance_fits_in_completed_requests <- pt_counts$nuisance_fits_in_completed_requests+length(nyt$auxiliary_reporting$fitted_views)
  rows <- nyt$metadata$catalog$rows
  expect_true(any(rows$original)); expect_true(any(!rows$original & rows$n_horizons>0))
  early <- nyt$design_all$links$s==1 & nyt$design_all$links$t==2 & nyt$design_all$links$g==4
  expect_identical(sort(nyt$design_all$links$control_pool[early]),c("{6,Inf}","{Inf}"))
  expect_identical(nyt$design_all$links$catalog_row,rows$row_id)
  expect_gt(nyt$design_all$pt_all_spec$n_zero_loading_links,0L)
  tm <- nyt$reporting$target_meta
  expect_true(any(tm$pt_window_tag=="in_pt_window" & !tm$pt_operator_covered,na.rm=TRUE))
})

test_that("the explicit development request carries the same PT-all fit and request settings", {
  d <- pt_ref$cp03_public_tiny_data()
  public <- pt_fit(d,pt_assumption="all",pt_start=2)
  pt_counts$development_requests <- 1L
  dev <- chained_att_gt_request(yname="y",tname="t",idname="id",gname="g",data=d,
    est_method="ipw",pt_assumption="all",pt_start=2,control=pt_control(d),
    run_numerical_transfer=FALSE,bstrap=FALSE,cband=FALSE)
  pt_counts$development_completed <- 1L
  pt_counts$nuisance_fits_in_completed_requests <- pt_counts$nuisance_fits_in_completed_requests+length(dev$chained$nuisance_diagnostics$records)
  expect_identical(dev$att,public$att)
  expect_identical(dev$inffunc,public$inffunc)
  expect_identical(dev$V_analytical,public$V_analytical)
  expect_identical(dev$chained$request$pt_start,2)
  expect_identical(dev$chained$pt,public$chained$pt)
})

test_that("D34 global identity map has exact post point, IF and covariance equivalence", {
  d <- pt_ref$cp03_public_tiny_data()
  for (weight in c("identity","two_step")) {
    post <- pt_fit(d,gmm_weighting=weight)
    all <- pt_fit(d,gmm_weighting=weight,pt_assumption="all",pt_start=4)
    pt_post_equal(all,post)
    expect_true(all$chained$pt$restriction_is_identity)
    expect_true(all$chained$pt$no_extra_restrictions)
    expect_equal(all$chained$pt$by_cohort$n_extra_calendar_dates,0)
    expect_match(paste(capture.output(print(all)),collapse="\n"),"No additional restrictions",fixed=TRUE)
    agg <- pt_agg(all,"group")
    expect_match(paste(capture.output(summary(agg)),collapse="\n"),"No additional restrictions",fixed=TRUE)
  }
})

test_that("stored windows, singleton records and full-frame transport remain distinct from represented restrictions", {
  d <- pt_ref$cp03_public_tiny_data()
  # Dates 1 and 2 exist only in singleton records. The sole retained pair is 3->4.
  keep <- d[d$t %in% c(3,4),]
  single <- d[d$id %in% c(1,5) & d$t %in% c(1,2),]
  single <- single[c(1,4),]; single$id <- c(9L,10L)
  sparse <- rbind(keep,single)
  fit <- pt_fit(sparse,pt_assumption="all",pt_start=2,gmm_weighting="identity")
  expect_identical(fit$chained$calendar,as.double(1:4))
  expect_equal(fit$n,10)
  expect_equal(sum(fit$chained$frame$n_obs==1),2)
  expect_equal(fit$chained$pt$by_cohort$n_extra_calendar_dates,1)
  expect_equal(fit$chained$pt$by_cohort$n_restricted_pre_parameters,0)
  expect_false(fit$chained$pt$no_extra_calendar_window)
  expect_true(fit$chained$pt$no_extra_restrictions)
  expect_match(paste(capture.output(summary(fit)),collapse="\n"),"No additional restrictions",fixed=TRUE)

  # Separate complete balanced moments, one singleton and one frame-only unit.
  one <- d[1,,drop=FALSE]; one$id <- 9L
  panel <- rbind(d,one)
  frame <- data.frame(id=1:10,g=c(rep(4,4),rep(0,4),4,0),cluster=c(rep(1:4,each=2),5,5))
  a <- pt_fit(panel,frame=frame,pt_assumption="all",pt_start=2)
  b <- pt_fit(panel,frame=frame,pt_assumption="all",pt_start=2,clustervars="cluster")
  expect_identical(a$att,b$att); expect_identical(a$inffunc,b$inffunc)
  expect_equal(a$n,10); expect_true(a$chained$frame$frame_only[10])
  post <- which(a$chained$role=="post")
  expect_equal(unname(a$inffunc[9:10,post]),c(0,0))
  ref <- pt_ref$cp03_public_tiny_reference("two_step")
  expect_equal(unname(a$inffunc[1:8,post]),ref$IF[,2]*10/8,tolerance=1e-10)
  expect_equal(b$chained$inference$n_zero_eligible_clusters,1)
})

test_that("a restricted bridge is newly identified without an unrestricted post fit, and D21 stays unavailable", {
  a <- pt_ref$cp03_public_tiny_data(); b <- a; b$id <- b$id+8L
  d <- rbind(a[a$t %in% c(2,4),],b[b$t %in% c(1,3),])
  extra <- d[1,,drop=FALSE]; extra$id <- 17L; extra$g <- 0; extra$t <- 5
  extra$y0 <- 17/4+5/2; extra$tau <- 0; extra$y <- extra$y0
  d <- rbind(d,extra)
  all <- pt_fit(d,pt_assumption="all",pt_start=2,gmm_weighting="identity")
  post <- which(all$chained$role=="post" & all$t==4)
  unsupported <- which(all$chained$role=="post" & all$t==5)
  expect_false(all$chained$available[unsupported]); expect_true(is.na(all$att[unsupported]))
  expect_true(all(is.na(all$inffunc[,unsupported])))
  expect_true(all(all$chained$available[post]))
  expect_true(all(is.finite(all$att[post])))
  expect_identical(all$chained$pretest$status,"unavailable")
  expect_null(all$W); expect_null(all$Wpval)
  expect_match(all$chained$pretest$reason,"no available-only subset",fixed=TRUE)
  # Build the unrestricted design for support diagnosis only: no post solve.
  full <- pt_pipeline(d,pt_assumption="all",pt_start=2,weighting="identity")
  cat <- full$metadata$catalog
  des <- build_comparison_catalog_design(cat,data.frame(g=4,t=4),data.frame(g=4,t=4))
  expect_false(des$post_targets$level_estimable)
  expect_identical(full$support$ledger$support_status,c("estimable","isolated_unsupported"))
  expect_true(full$design_all$post_targets$level_estimable[full$design_all$post_targets$time==4])
  expect_false(full$design_all$post_targets$level_estimable[full$design_all$post_targets$time==5])
})

test_that("irregular calendar with actual anticipation retains mixed cohorts and free anticipation effects", {
  cal <- c(2000,2002,2005,2008,2009,2014)
  d <- expand.grid(id=1:24,t=cal); d$g <- rep(c(2002,2009,0),each=8)[d$id]
  # Explicit potential outcomes: deterministic Y0 and heterogeneous timing;
  # tau=3/4 during the one-year anticipation interval, 3/2 after adoption.
  d$y0 <- d$id/4+d$t/500+sin(d$id*d$t)/10
  d$tau <- ifelse(d$g>0 & d$t>=d$g,1.5,ifelse(d$g>0 & d$t<d$g & d$t+1>=d$g,.75,0))
  d$y <- d$y0+d$tau
  fit <- pt_fit(d,anticipation=1,pt_assumption="all",pt_start=2002,gmm_weighting="identity",base_period="universal")
  expect_identical(fit$chained$target_cohorts,c(2002,2009))
  by <- fit$chained$pt$by_cohort
  expect_identical(by$ref,c(2000,2005))
  expect_identical(by$no_extra_restrictions,c(TRUE,FALSE))
  expect_equal(by$n_extra_calendar_dates,c(0,1))
  tags <- fit$chained$pretest$window_tags
  ant <- match("pre(g=2009,t=2008)",tags$target_id)
  expect_identical(tags$pt_window_tag[ant],"anticipation")
  expect_true(is.finite(fit$att[match(tags$target_id[ant],fit$chained$target_id)]))
  expect_match(fit$chained$pt$assumptions,"including dates before the start",fixed=TRUE)
})

test_that("identity IF equals a raw weighted finite difference and fixed draws replay the existing bootstrap law", {
  d <- pt_ref$cp03_public_tiny_data()
  center <- pt_fit(d,gmm_weighting="identity",pt_assumption="all",pt_start=2,weightsname="weight")
  h <- 1e-5; plus <- minus <- d
  plus$weight[plus$id==1] <- 1+h; minus$weight[minus$id==1] <- 1-h
  a <- pt_fit(plus,gmm_weighting="identity",pt_assumption="all",pt_start=2,weightsname="weight")
  b <- pt_fit(minus,gmm_weighting="identity",pt_assumption="all",pt_start=2,weightsname="weight")
  post <- which(center$chained$role=="post")
  expect_equal((a$att[post]-b$att[post])/(2*h),center$inffunc[1,post]/8,tolerance=1e-8)
  patterns <- unname(as.matrix(expand.grid(rep(list(c(-1,1)),8L))))
  signs <- patterns[rep(seq_len(256L),length.out=1000L),,drop=FALSE]
  z <- d; z$g[z$g==0] <- Inf; ctrl <- pt_control(d,"identity")
  run <- function(pt) {
    pt_counts$pipeline_fits <- pt_counts$pipeline_fits+1L
    pt_counts$fixed_multiplier_rows <- pt_counts$fixed_multiplier_rows+1000L
    out <- run_chained_ipw_pipeline(z,"id","t","y","g",ctrl,weighting="identity",control_group="nevertreated",
      pt_assumption=pt,pt_start=if(pt=="all") 4 else NULL,run_numerical_transfer=FALSE,
      bstrap=TRUE,cband=TRUE,biters=1000L,fixed_signs=signs,aggregations=character())
    pt_counts$nuisance_fits_in_completed_requests <- pt_counts$nuisance_fits_in_completed_requests+length(out$auxiliary_reporting$fitted_views)
    out
  }
  set.seed(3103); before <- .Random.seed
  x <- run("post"); y <- run("all")
  expect_identical(.Random.seed,before)
  expect_identical(x$reporting$estimates,y$reporting$estimates)
  expect_identical(x$reporting$if_full,y$reporting$if_full)
  expect_identical(x$reporting$se,y$reporting$se)
  expect_identical(x$reporting$cband_crit_val,y$reporting$cband_crit_val)
  expect_identical(y$controls$btype,"rademacher"); expect_identical(y$controls$biters,1000L)
  path <- Sys.getenv("CP03_PUBLIC_ARTIFACT_DIR")
  if (nzchar(path)) saveRDS(as.list(pt_counts),file.path(path,"public_fixture_counts.rds"))
})
