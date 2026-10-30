if (!('cdid' %in% loadedNamespaces())) stop('Load the installed candidate namespace.')
raw_ns<-asNamespace('cdid')
ref_path<-testthat::test_path('reference','raw_logit_independent_reference.R')
if(!identical(digest::digest(file=ref_path,algo='sha256'),'7fda102445b1be51d13c963efeffd4ed9d3ff54d60fe5206266b7df60dec0b7e'))stop('Frozen raw reference binding differs.')
sys.source(ref_path,envir=environment())
test_that('Frozen pre-access Fraction fixtures bind actual installed raw-logit arithmetic', {
 control<-get('logit_certificate_control',raw_ns)(18,700,11,3,1000000,10000000)
 run<-get('certify_chained_logit_root',raw_ns)
 include<-function(iv,bracket) {
  expect_true(is.numeric(iv)&&length(iv)==2L&&all(is.finite(iv)))
  expect_lte(iv[1L],bracket[1L]);expect_gte(iv[2L],bracket[2L])
 }
 for(f in raw_logit_independent_cases) {
  args<-f$input;args$control<-control
  if(!f$expected$valid_domain) {
   expect_error(do.call(run,args),class='cdid_logit_certificate_error')
   next
  }
  z<-do.call(run,args);e<-f$expected
  expect_identical(z$dimensions$n_a_units,e$n)
  expect_false(z$separation_certified);expect_false(z$sampling_certified)
  for(name in c('M_cohort','Pi_cohort','r_cohort','g','whitened_g'))
   for(j in seq_along(e[[name]]))include(z$enclosures[[name]][j,],e[[name]][[j]])
  for(name in c('H','M','whitened_H_minus_I'))
   for(j in seq_along(e[[name]]))for(k in seq_along(e[[name]][[j]]))
    include(z$enclosures[[name]][j,k,],e[[name]][[j]][[k]])
  for(name in c('t_iv','delta_iv','a_iv'))include(z$enclosures[[name]],e[[name]])
  if(f$name %in% c('balanced_n16_nondiagonal','balanced_n32_nonobservers','weighted_n14_intercept'))expect_true(z$certified)
  if(f$name %in% c('asymmetric_weighted_nonroot_small_radius','rounded_rate_cancellation',
   'rank_deficient_observers','strict_x_boundary','strict_q_boundary','nonzero_beta_exact_predictors','x_above_one'))
   expect_identical(z$status,'INCONCLUSIVE')
  expect_lte(z$work$actual_interval_ops,z$work$preflight_max_interval_ops)
 }
})
