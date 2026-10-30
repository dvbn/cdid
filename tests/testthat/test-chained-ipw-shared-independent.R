# Independent shared-IPW checks. Expected values come only from the packaged
# immutable 6550dc kernel plus ONE matched-start adaptation and the unchanged
# independent Newton/IF reference, both hash-bound by their loaders.
# Spy tests use isolated clones of
# the actual cdid namespace closures (audit_spy); all other calls use cdid:::.

if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}

test_that('fresh weighted and unit-weight fixtures match the immutable kernel plus ONE matched-start adaptation bitwise and every independent IF component', {
  original <- audit_original_kernel()
  independent_fit <- audit_newton_reference()$independent_fit
  # The frozen kernel accepts three fields; native-only fields must be defaults.
  legacy3 <- audit_control[c("epsilon", "maxit", "trace")]
  expect_identical(do.call(stats::glm.control, legacy3), audit_control)
  for(weighted in c(FALSE,TRUE)) {
    fx <- audit_fixture(weighted); z <- audit_setup(fx)
    oldctx <- original$construct_ipw_context(fx$index,id=seq_len(fx$n),Z=fx$Z,weights=fx$w)
    for(j in 1:3) {
      p <- audit_pairs[[j]]
      f <- cdid:::fit_ipw_link_shared(z$session,z$declarations[[j]])
      old <- original$fit_ipw_link(oldctx,4,c(7,Inf),p[1],p[2],legacy3,functional_key='A')
      field <- setdiff(names(old),c('context_ref','diagnostics'))
      expect_identical(unclass(f)[field],unclass(old)[field])
      expect_identical(names(f$diagnostics),c(names(old$diagnostics),'clipping'))
      expect_identical(f$diagnostics[names(old$diagnostics)],old$diagnostics)
      expect_true(is.list(f$diagnostics$clipping) && identical(f$diagnostics$clipping$eps,1e-6))
      expect_identical(as.integer(f$diagnostics$clipping$p_n),f$counts$n_observed_risk_units)
      expect_identical(as.integer(f$diagnostics$clipping$r_n),f$cohorts$n_cohorts)
      a <- cdid:::materialize_ipw_if(fx$context,f,fx$n,11,TRUE)
      o <- original$materialize_ipw_if(oldctx,old,fx$n,11,TRUE)
      expect_identical(a$components,o$components)
      expect_identical(a$phi,o$phi)
      ref <- independent_fit(fx,p[1],p[2])
      expect_lt(abs(f$estimate-ref$estimate),1e-11)
      expect_lt(max(abs(f$coefficients-ref$beta)),1e-10)
      for(nm in names(ref$components)) {
        expect_length(a$components[[nm]],fx$n)
        expect_lt(max(abs(drop(a$components[[nm]])-ref$components[[nm]])),1e-9)
      }
      expect_lt(max(abs(a$phi-ref$phi)),1e-9)
      h <- 1e-5
      fd <- vapply(seq_len(fx$n),function(i) {
        wp<-wm<-fx$w; wp[i]<-wp[i]*(1+h);wm[i]<-wm[i]*(1-h)
        fx$n*(independent_fit(fx,p[1],p[2],wp)$estimate-independent_fit(fx,p[1],p[2],wm)$estimate)/(2*h)
      },numeric(1))
      expect_lt(max(abs(fd-a$phi)),2e-6)
      non <- fx$G!=9 & !fx$obs[,1]
      expect_true(any(abs(a$phi[non])>1e-5))
      expect_identical(max(abs(a$phi[fx$G==9])),0)
      # Independently break each complete IF component; every mutation matters.
      for(nm in names(ref$components)) expect_gt(max(abs(a$phi-(ref$phi-ref$components[[nm]]))),1e-7)
    }
    ord <- rev(seq_len(fx$n))
    context <- cdid:::construct_ipw_context(fx$index,id=ord,Z=fx$Z[ord,,drop=FALSE],weights=fx$w[ord])
    f <- cdid:::fit_ipw_link(context,4,c(7,Inf),1,2,audit_control)
    expect_identical(f$estimate,cdid:::fit_ipw_link(fx$context,4,c(7,Inf),1,2,audit_control)$estimate)
  }
})

test_that('one nuisance fit, retained-only outcomes/replay and disconnected sparse forest are observed directly', {
  fx<-audit_fixture(); sp<-audit_spy(); e<-sp$e
  ses<-e$construct_ipw_shared_session(fx$context,2,20,1e5)
  nu<-e$construct_ipw_nuisance(ses,4,c(7,Inf),1,2,audit_control)
  expect_identical(sp$count$glm,1L);expect_length(sp$count$out,0)
  pp<-list(c(1,2),c(2,3),c(1,3),c(4,5),c(5,6),c(4,6))
  dl<-lapply(pp,function(p)e$declare_ipw_shared_link(ses,nu,4,c(7,Inf),p[1],p[2],audit_control,functional_key='A'))
  des<-e$build_chained_design(data.frame(g=4,s=c(1,2,1,4,5,4),t=c(2,3,3,5,6,6),ref=3,functional_key='A'))
  cert<-e$certify_ipw_structural_classes(ses,dl,des,1:6,'optimal',audit_bounds)
  expect_identical(cert$classes[[1]]$roots,c(1,4)); expect_identical(cert$dropped_index,c(2L,5L))
  expect_identical(sp$count$glm,1L);expect_length(sp$count$out,0)
  fit<-lapply(dl[cert$retained_index],function(d)e$fit_ipw_link_shared(ses,d))
  e$replay_ipw_if(fx$context,fit,13,function(batch) invisible(NULL))
  expect_identical(sp$count$glm,1L)
  expect_true(all(vapply(sp$count$out,function(p)any(vapply(pp[cert$retained_index],function(q)identical(as.double(q),p),logical(1))),logical(1))))
  bad<-audit_bounds;bad$max_E_entries<-0
  sp$count$paths<-0L
  expect_error(e$certify_ipw_structural_classes(ses,dl,des,1:6,'optimal',bad),class='cdid_ipw_bound_error')
  expect_identical(sp$count$paths,0L)
})

test_that('bounds and malformed/overflow contexts refuse before mask allocation or GLM', {
  fx<-audit_fixture();sp<-audit_spy();e<-sp$e
  tiny<-e$construct_ipw_shared_session(fx$context,1,1,1)
  expect_error(e$construct_ipw_nuisance(tiny,4,c(7,Inf),1,2,audit_control),class='cdid_ipw_bound_error')
  expect_identical(sp$count$mask,0L);expect_identical(sp$count$glm,0L)
  mut<-fx$context;mut$w_a[1]<-mut$w_a[1]*2
  expect_error(e$fit_ipw_link(mut,4,c(7,Inf),1,2,audit_control),class='cdid_ipw_context_mismatch_error')
  expect_identical(sp$count$glm,0L)
  pan<-fx$panel
  pan$y<-ifelse(pan$time==1,-1e308,1e308)
  expect_error(e$construct_observed_pairs(pan,'id','time','y','g'),class='cdid_pairs_outcome_error')
  # Deliberately mutated compact input reaches the old query guard even though
  # the normal pair constructor already refuses this overflowing outcome pair.
  idx<-fx$index; idx$observations$y<-ifelse(idx$observations$time==1,-1e308,1e308)
  ctx<-e$construct_ipw_context(idx,id=seq_len(fx$n),Z=fx$Z,weights=fx$w)
  expect_error(e$fit_ipw_link(ctx,4,c(7,Inf),1,2,audit_control),class='cdid_ipw_outcome_error')
  expect_identical(sp$count$glm,0L)
  nu<-e$construct_ipw_nuisance(e$construct_ipw_shared_session(ctx,1,1,1e5),4,c(7,Inf),1,2,audit_control)
  expect_s3_class(nu,'cdid_ipw_nuisance')
  expect_identical(sp$count$glm,1L)
})

test_that('audit spy clones and packaged references are bound to the loaded cdid namespace, not a snapshot', {
  ns <- asNamespace("cdid")
  sp <- audit_spy()
  expect_identical(parent.env(sp$e), ns)
  # Compare every clone while keeping the assertion count independent of
  # unrelated package modules added to the namespace.
  unhooked <- setdiff(sp$cloned, names(audit_spy_hooks))
  expect_true(all(vapply(unhooked, function(nm)
    identical(environment(get(nm, envir = ns, inherits = FALSE)), ns), logical(1L))))
  expect_identical(lapply(unhooked, function(nm) body(get(nm, envir = sp$e, inherits = FALSE))),
                   lapply(unhooked, function(nm) body(get(nm, envir = ns, inherits = FALSE))))
  expect_identical(lapply(unhooked, function(nm) formals(get(nm, envir = sp$e, inherits = FALSE))),
                   lapply(unhooked, function(nm) formals(get(nm, envir = ns, inherits = FALSE))))
  for (nm in names(audit_spy_hooks)) {
    expect_identical(sp$real[[nm]], get(nm, envir = ns, inherits = FALSE))
    expect_identical(sp$real[[nm]], utils::getFromNamespace(nm, "cdid"))
  }
  for (nm in c("construct_ipw_shared_session", "construct_ipw_nuisance", "declare_ipw_shared_link",
               "certify_ipw_structural_classes", "fit_ipw_link_shared", "fit_ipw_link", "replay_ipw_if",
               "materialize_ipw_if", "construct_ipw_context", "construct_observed_pairs", "build_chained_design")) {
    expect_true(nm %in% sp$cloned)
    expect_identical(environment(utils::getFromNamespace(nm, "cdid")), ns)
  }
  ref <- audit_newton_reference()
  expect_identical(ls(ref, all.names = TRUE), "independent_fit")
  expect_identical(parent.env(environment(ref$independent_fit)), baseenv())
  original <- audit_original_kernel()
  expect_identical(parent.env(original), baseenv())
  expect_false(exists("fit_ipw_link_shared", envir = original, inherits = FALSE))
})
