# Independent shared-IPW audit regressions: actual-incidence W != E W_R, intrinsic roles and
# parameter maps, sealed certificates, cached-fit and storage-counter swaps,
# and forged nuisance cohort/N/risk-count metadata must all refuse. Every call
# goes to the loaded cdid namespace; fixtures live in helper-cp02-shared-audit.R.

if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}

test_that('malformed actual incidence cannot obtain a structural range certificate', {
  z <- audit_setup()
  expect_s3_class(audit_cert(z),'cdid_ipw_structural_certificate')
  cases <- list(
    incidence_end = list(
      edit = function(d) {d$incidence$end[2] <- d$incidence$start[2]; d},
      class = 'cdid_ipw_design_error',
      message = 'Actual incidence differs from the intrinsic link and parameter maps.'),
    n_parameters = list(
      edit = function(d) {d$incidence$n_parameters <- d$incidence$n_parameters+1L; d},
      class = 'cdid_ipw_design_error',
      message = 'Actual incidence differs from the intrinsic link and parameter maps.'),
    start_parameter = list(
      edit = function(d) {d$links$start_param_index[2] <- 0L; d},
      class = 'cdid_ipw_design_error',
      message = 'Actual incidence differs from the intrinsic link and parameter maps.'),
    vertex_parameter = list(
      edit = function(d) {d$vertices$param_index[1] <- 99L; d},
      class = 'cdid_ipw_class_vertex_error',
      message = "Endpoint vertex is not the intrinsic qualified vertex of the link's cohort, date and operator."),
    parameter_date = list(
      edit = function(d) {d$parameters$time[1] <- d$parameters$time[1]+1; d},
      class = 'cdid_ipw_design_error',
      message = 'Parameter rows differ from the intrinsic qualified vertices.'),
    incidence_ids = list(
      edit = function(d) {d$incidence$param_ids <- rev(d$incidence$param_ids); d},
      class = 'cdid_ipw_design_error',
      message = 'Actual incidence differs from the intrinsic link and parameter maps.'))
  for (case in cases) {
    expect_error(audit_cert(z,case$edit(z$design)), class=case$class,
      regexp=case$message, fixed=TRUE)
  }
})

test_that('qualified and unqualified endpoint roles are intrinsic and complete', {
  z <- audit_setup()
  for(edit in list(
    function(d) {v <- d$links$start_vertex_index[1]; d$vertices$role[v] <- 'reference'; d$vertices$functional_key[v] <- 'foreign'; d},
    function(d) {v <- d$links$start_vertex_index[1]; d$vertices$role[v] <- NA_character_; d},
    function(d) {d$links$ref <- 2; d})) {
    expect_error(audit_cert(z,edit(z$design)),class='cdid_ipw_certificate_error')
  }
  pp <- list(c(3,4),c(4,5),c(3,5))
  ds <- lapply(pp,function(p) cdid:::declare_ipw_shared_link(z$session,z$nuisance,4,c(7,Inf),p[1],p[2],audit_control,functional_key='A'))
  L <- data.frame(g=4,s=c(3,4,3),t=c(4,5,5),ref=3,functional_key='A',s_qualified=c(TRUE,FALSE,TRUE),t_qualified=FALSE)
  d <- cdid:::build_chained_design(L)
  good <- cdid:::certify_ipw_structural_classes(z$session,ds,d,1:3,'optimal',audit_bounds)
  expect_identical(good$dropped_index,2L)
  v <- d$links$end_vertex_index[1]
  d$vertices$role[v] <- 'post'
  expect_error(cdid:::certify_ipw_structural_classes(z$session,ds,d,1:3,'optimal',audit_bounds),class='cdid_ipw_class_vertex_error')
})

test_that('issued certificates and cache counters cannot silently change meaning', {
  z <- audit_setup()
  cert <- audit_cert(z)
  bad <- cert; bad$E$x <- -bad$E$x
  expect_error(cdid:::apply_ipw_reduction_map(bad,c(2,5)),class='cdid_ipw_certificate_error')
  fake <- structure(list(n_links=3L,retained_index=c(1L,3L),E=list(i=2L,j=1L,x=1L)),class='cdid_ipw_structural_certificate')
  expect_error(cdid:::apply_ipw_reduction_map(fake,c(2,5)),class='cdid_ipw_certificate_error')
  expect_invisible(cdid:::verify_ipw_structural_certificate(cert,z$session,z$declarations,z$design,1:3))
  f1 <- cdid:::fit_ipw_link_shared(z$session,z$declarations[[1]])
  f2 <- cdid:::fit_ipw_link_shared(z$session,z$declarations[[2]])
  z$session$ref$fit_refs[[1]] <- f2$link_ref
  expect_error(cdid:::fit_ipw_link_shared(z$session,z$declarations[[1]]),class='cdid_ipw_session_error')
  z$session$ref$unit_entries <- 0
  expect_error(cdid:::validate_ipw_shared_session(z$session),class='cdid_ipw_session_error')
})

test_that('intrinsic nuisance metadata derives from g K and the frozen context', {
  for(edit in list(
    function(nu,z) {nu$g<-8;nu},
    function(nu,z) {nu$stage$n<-nu$stage$n/2;nu$hessian<-cdid:::ipw_nuisance_hessian_stage(cdid:::ipw_stage_from_nuisance(nu,z$fx$context),nu$cohort_labels);nu},
    function(nu,z) {nu$mask$n_risk_a_units<-1L;nu$mask$n_pair_nonobservers<-0L;nu},
    function(nu,z) {nu$cohort_labels<-rev(nu$cohort_labels);nu})) {
    z<-audit_setup();nu<-z$nuisance;nu$ref<-NULL;nu<-edit(nu,z);nu$serial<-2L
    nu<-cdid:::ipw_seal(nu,'nuisance');z$session$ref$nuisance_refs[[2]]<-nu$ref
    z$session$ref$unit_entries<-z$session$ref$unit_entries+nu$unit_entries
    d<-cdid:::declare_ipw_shared_link(z$session,nu,nu$g,nu$K,1,2,audit_control,functional_key='A')
    des<-cdid:::build_chained_design(data.frame(g=nu$g,s=1,t=2,ref=3,functional_key='A'))
    expect_error(cdid:::certify_ipw_structural_classes(z$session,list(d),des,1L,'optimal',audit_bounds),class='cdid_ipw_nuisance_integrity_error')
  }
})
