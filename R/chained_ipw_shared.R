# Private CP02 shared-nuisance weighted-pair IPW layer and structural classes.
# One outcome-free nuisance object (rates, logit, H factor) is fitted ONCE for an
# exact risk-set pair-observation mask and explicitly reused by every link whose
# own mask is identical. Classes of Prop. 5 / Cor. 6 (gmm_computation.tex) are
# certified from actual object identity, recomputed masks and the intrinsic
# design vertex map; labels, counts, rates or numeric agreement never certify.
# All state lives in a caller-explicit bounded session; there is no hidden cache.
# No public estimator, catalog search, GMM solver or numerical acceptance level.

ipw_seal_ok <- function(obj, binding) {
  ref <- obj$ref
  is.environment(ref) && environmentIsLocked(ref) &&
    exists(binding, envir = ref, inherits = FALSE) &&
    bindingIsLocked(binding, ref) && identical(obj, get(binding, envir = ref, inherits = FALSE))
}

ipw_seal <- function(obj, binding) {
  ref <- new.env(parent = emptyenv())
  obj$ref <- ref
  assign(binding, obj, envir = ref)
  lockEnvironment(ref, bindings = TRUE)
  obj
}

ipw_serial_ok <- function(k, n) {
  is.integer(k) && length(k) == 1L && !is.na(k) && k >= 1L && k <= n
}

construct_ipw_shared_session <- function(context, max_nuisance, max_links, max_unit_entries) {
  validate_ipw_context(context)
  if (missing(max_nuisance) || missing(max_links) || missing(max_unit_entries)) {
    cdid_abort(
      "Supply explicit `max_nuisance`, `max_links` and `max_unit_entries` session bounds.",
      c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
    )
  }
  bounds <- list(
    max_nuisance = ipw_bound(max_nuisance, "max_nuisance"),
    max_links = ipw_bound(max_links, "max_links"),
    max_unit_entries = ipw_bound(max_unit_entries, "max_unit_entries")
  )
  ref <- new.env(parent = emptyenv())
  ref$context <- context
  ref$bounds <- bounds
  ref$nuisance_refs <- list()
  ref$declaration_refs <- list()
  ref$fit_refs <- list()
  ref$unit_entries <- 0
  lockEnvironment(ref)
  lockBinding("context", ref)
  lockBinding("bounds", ref)
  structure(list(ref = ref), class = "cdid_ipw_shared_session")
}

validate_ipw_shared_session <- function(session) {
  fields <- c("context", "bounds", "nuisance_refs", "declaration_refs", "fit_refs", "unit_entries")
  ref <- if (is.list(session)) session$ref else NULL
  ok <- inherits(session, "cdid_ipw_shared_session") && is.list(session) &&
    identical(names(session), "ref") && is.environment(ref) && environmentIsLocked(ref) &&
    setequal(ls(ref, all.names = TRUE), fields) &&
    bindingIsLocked("context", ref) && bindingIsLocked("bounds", ref)
  if (ok) {
    b <- ref$bounds
    ok <- is.list(b) && identical(names(b), c("max_nuisance", "max_links", "max_unit_entries")) &&
      is.list(ref$nuisance_refs) && is.list(ref$declaration_refs) && is.list(ref$fit_refs) &&
      length(ref$nuisance_refs) <= b$max_nuisance &&
      length(ref$declaration_refs) <= b$max_links &&
      length(ref$fit_refs) == length(ref$declaration_refs) &&
      ipw_scalar_numeric(ref$unit_entries) && ref$unit_entries >= 0 &&
      ref$unit_entries <= b$max_unit_entries
  }
  if (!ok) {
    cdid_abort("Expected a well-formed `cdid_ipw_shared_session`.",
               c("cdid_ipw_session_error", "cdid_ipw_input_error"))
  }
  validate_ipw_context(ref$context)
  actual_entries <- 0
  for (k in seq_along(ref$nuisance_refs)) {
    nr <- ref$nuisance_refs[[k]]
    if (!is.environment(nr) || !exists("nuisance", envir = nr, inherits = FALSE)) {
      cdid_abort("Session nuisance registry is malformed.", c("cdid_ipw_session_error", "cdid_ipw_input_error"))
    }
    nu <- get("nuisance", envir = nr, inherits = FALSE)
    validate_ipw_nuisance(nu, ref)
    if (!identical(nu$serial, k) || !identical(nu$unit_entries, 8 * as.double(length(nu$mask$observed_a_index)))) {
      cdid_abort("Session nuisance storage accounting is malformed.", c("cdid_ipw_session_error", "cdid_ipw_input_error"))
    }
    actual_entries <- actual_entries + nu$unit_entries
  }
  if (!identical(ref$unit_entries, actual_entries)) {
    cdid_abort("Session unit storage counter differs from its issued nuisance objects.",
               c("cdid_ipw_session_error", "cdid_ipw_input_error"))
  }
  ref
}

validate_ipw_nuisance <- function(nuisance, sref) {
  if (!inherits(nuisance, "cdid_ipw_nuisance") || !is.list(nuisance) || !ipw_seal_ok(nuisance, "nuisance")) {
    cdid_abort(
      "Shared nuisance object is malformed or differs from its frozen snapshot; it cannot be reused.",
      c("cdid_ipw_nuisance_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  k <- nuisance$serial
  if (!identical(nuisance$session_ref, sref) || !ipw_serial_ok(k, length(sref$nuisance_refs)) ||
      !identical(sref$nuisance_refs[[k]], nuisance$ref)) {
    cdid_abort(
      "Shared nuisance object was not issued by this session; it cannot be reused.",
      c("cdid_ipw_nuisance_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  if (!identical(nuisance$context_ref, sref$context$ref)) {
    cdid_abort(
      "Shared nuisance object belongs to a different `cdid_ipw_context`.",
      c("cdid_ipw_context_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  invisible(nuisance)
}

# Outcome-free nuisance construction: rates, weighted logit, mean denominators,
# Hessian factor and every score-rate column, for the EXACT risk-set
# pair-observation mask of `(s, t)`. Outcomes are never read here.
construct_ipw_nuisance <- function(session, g, K, s, t, control) {
  sref <- validate_ipw_shared_session(session)
  context <- sref$context
  ctrl <- validate_ipw_glm_control(control)
  coh <- validate_ipw_cohorts(g, K)
  desc <- make_date_pair_descriptor(s, t, context)
  if (length(sref$nuisance_refs) >= sref$bounds$max_nuisance) {
    cdid_abort("Session nuisance cache bound reached; no object is evicted or refitted.",
               c("cdid_ipw_bound_error", "cdid_ipw_cache_bound_error"))
  }
  sup <- ipw_cohort_support(context, coh)
  # Preflight: eight retained risk-unit vectors, bounded before any mask query.
  upper <- 8 * sum(as.double(sup$group_a_count))
  if (upper > sref$bounds$max_unit_entries - sref$unit_entries) {
    cdid_abort("Nuisance unit storage would exceed the session `max_unit_entries` bound.",
               c("cdid_ipw_bound_error", "cdid_ipw_cache_bound_error"))
  }

  hit <- query_pair_mask_internal(
    context, desc$s, desc$t,
    ds = desc$s_calendar_index, dt = desc$t_calendar_index
  )
  risk <- ipw_pair_risk_set(context, coh$cohorts, hit)
  st <- ipw_nuisance_fit_stage(context, sup, risk, ctrl)
  cohort_labels <- ipw_cohort_labels(coh)
  hs <- ipw_nuisance_hessian_stage(st, cohort_labels)
  if (any(!is.finite(hs$A_beta_r))) {
    cdid_abort(
      "Nonfinite propensity or outcome partial derivatives.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }

  stage <- st[c("n", "m_cohorts", "p", "risk_gi", "risk_w", "pair_count", "pair_mass",
                "rates", "fit_weights", "beta", "odds_use", "prob_use", "is_t",
                "a_obs", "b_obs", "sum_a", "mean_a", "sum_b", "mean_b", "fit_info")]
  n_obs <- length(risk$risk_a_idx)
  serial <- length(sref$nuisance_refs) + 1L
  nu <- structure(
    list(
      serial = serial,
      g = coh$g,
      K = coh$K,
      cohorts = coh$cohorts,
      cohort_labels = cohort_labels,
      control = ctrl,
      defining_pair = list(s = desc$s, t = desc$t,
                           s_calendar_index = desc$s_calendar_index,
                           t_calendar_index = desc$t_calendar_index),
      support = sup,
      mask = list(
        observed_a_index = risk$risk_a_idx,
        n_risk_a_units = sum(sup$group_a_count),
        n_observed_risk_units = n_obs,
        n_pair_nonobservers = sum(sup$group_a_count) - n_obs
      ),
      stage = stage,
      hessian = hs,
      unit_entries = 8 * as.double(n_obs),
      context_ref = context$ref,
      session_ref = sref
    ),
    class = "cdid_ipw_nuisance"
  )
  nu <- ipw_seal(nu, "nuisance")
  sref$nuisance_refs[[serial]] <- nu$ref
  sref$unit_entries <- sref$unit_entries + nu$unit_entries
  nu
}

# Recompute a link's risk-set mask from calendar postings and require EXACT
# identity with the nuisance mask on every A-unit of {g, K}.
ipw_shared_mask_check <- function(context, nuisance, desc) {
  hit <- query_pair_mask_internal(
    context, desc$s, desc$t,
    ds = desc$s_calendar_index, dt = desc$t_calendar_index
  )
  gi <- match(context$a_g[hit$a_index], nuisance$cohorts)
  in_risk <- which(!is.na(gi))
  if (!identical(hit$a_index[in_risk], nuisance$mask$observed_a_index)) {
    cdid_abort(
      "Link risk-set pair-observation mask differs from the shared nuisance mask; reuse is refused (no refit, no sample change).",
      c("cdid_ipw_shared_mask_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  list(hit = hit, in_risk = in_risk, n_obs_pair_all = length(hit$a_index))
}

validate_ipw_shared_declaration <- function(declaration, sref) {
  if (!inherits(declaration, "cdid_ipw_shared_link_declaration") || !is.list(declaration) ||
      !ipw_seal_ok(declaration, "declaration")) {
    cdid_abort(
      "Shared link declaration is malformed or differs from its frozen snapshot.",
      c("cdid_ipw_declaration_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  k <- declaration$serial
  if (!identical(declaration$session_ref, sref) ||
      !ipw_serial_ok(k, length(sref$declaration_refs)) ||
      !identical(sref$declaration_refs[[k]], declaration$ref)) {
    cdid_abort(
      "Shared link declaration was not issued by this session.",
      c("cdid_ipw_declaration_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  nu_ref <- declaration$nuisance_ref
  if (!is.environment(nu_ref) || !exists("nuisance", envir = nu_ref, inherits = FALSE)) {
    cdid_abort("Shared link declaration lacks its nuisance object.",
               c("cdid_ipw_declaration_error", "cdid_ipw_shared_incompatible_error"))
  }
  nu <- get("nuisance", envir = nu_ref, inherits = FALSE)
  validate_ipw_nuisance(nu, sref)
  if (!identical(nu$ref, nu_ref)) {
    cdid_abort("Shared link declaration references a foreign nuisance environment.",
               c("cdid_ipw_declaration_error", "cdid_ipw_shared_incompatible_error"))
  }
  nu
}

# Explicit reuse request: refuses unless context, target g, EXACT ordered K,
# validated glm.control and the full risk-set mask are identical. Outcome-free.
declare_ipw_shared_link <- function(session, nuisance, g, K, s, t, control,
                                    link_id = NULL, functional_key = NULL) {
  sref <- validate_ipw_shared_session(session)
  context <- sref$context
  validate_ipw_nuisance(nuisance, sref)
  ctrl <- validate_ipw_glm_control(control)
  coh <- validate_ipw_cohorts(g, K)
  if (!identical(coh$g, nuisance$g) || !identical(coh$K, nuisance$K) ||
      !identical(ctrl, nuisance$control)) {
    cdid_abort(
      "Target `g`, exact ordered `K` and validated `glm.control` must be identical to the shared nuisance object.",
      c("cdid_ipw_shared_incompatible_error", "cdid_ipw_input_error")
    )
  }
  desc <- make_date_pair_descriptor(s, t, context)
  validate_ipw_link_labels(link_id, functional_key)
  if (length(sref$declaration_refs) >= sref$bounds$max_links) {
    cdid_abort("Session link bound reached.", c("cdid_ipw_bound_error", "cdid_ipw_cache_bound_error"))
  }
  chk <- ipw_shared_mask_check(context, nuisance, desc)
  serial <- length(sref$declaration_refs) + 1L
  decl <- structure(
    list(
      serial = serial,
      pair = list(s = desc$s, t = desc$t,
                  s_calendar_index = desc$s_calendar_index,
                  t_calendar_index = desc$t_calendar_index,
                  link_id = link_id, functional_key = functional_key),
      g = nuisance$g,
      K = nuisance$K,
      n_observed_pair_units = chk$n_obs_pair_all,
      nuisance_serial = nuisance$serial,
      nuisance_ref = nuisance$ref,
      context_ref = context$ref,
      session_ref = sref
    ),
    class = "cdid_ipw_shared_link_declaration"
  )
  decl <- ipw_seal(decl, "declaration")
  sref$declaration_refs[[serial]] <- decl$ref
  sref$fit_refs[serial] <- list(NULL)
  decl
}

# Rebuild the kernel stage list from the stored nuisance (no refit).
ipw_stage_from_nuisance <- function(nuisance, context) {
  st <- nuisance$stage
  st$risk_a_idx <- nuisance$mask$observed_a_index
  st$Z_use <- context$Z_a[st$risk_a_idx, , drop = FALSE]
  st$D_use <- as.double(st$risk_gi == 1L)
  st$is_c <- !st$is_t
  st
}

# Outcome stage for a declared link: reads only this link's pair outcomes and
# reuses the shared nuisance object and its single H factor.
fit_ipw_link_shared <- function(session, declaration) {
  sref <- validate_ipw_shared_session(session)
  context <- sref$context
  nu <- validate_ipw_shared_declaration(declaration, sref)
  k <- declaration$serial
  existing <- sref$fit_refs[[k]]
  if (!is.null(existing)) {
    if (!is.environment(existing) || !exists("fit", envir = existing, inherits = FALSE)) {
      cdid_abort("Session fit registry is malformed.", c("cdid_ipw_session_error", "cdid_ipw_input_error"))
    }
    fit <- get("fit", envir = existing, inherits = FALSE)
    validate_ipw_shared_link_fit(fit, context)
    if (!identical(fit$declaration_ref, declaration$ref) || !identical(fit$nuisance_ref, nu$ref)) {
      cdid_abort("Session fit registry does not match the requested declaration.",
                 c("cdid_ipw_session_error", "cdid_ipw_input_error"))
    }
    return(fit)
  }
  desc <- declaration$pair
  chk <- ipw_shared_mask_check(context, nu, desc)
  pq <- ipw_pair_outcomes(context, chk$hit)
  risk_dy <- pq$dy[chk$in_risk]
  st <- ipw_stage_from_nuisance(nu, context)
  means <- ipw_link_means(st, risk_dy)
  coh <- list(g = nu$g, K = nu$K, cohorts = nu$cohorts)
  part <- ipw_link_partials(st, nu$hessian, risk_dy, means, nu$support$q, nu$cohort_labels)
  fit <- ipw_assemble_link_fit(context, desc, coh, nu$support, chk$n_obs_pair_all, st,
                               nu$hessian, means, part, nu$control,
                               desc$link_id, desc$functional_key)
  fit$shared <- list(nuisance_serial = nu$serial, declaration_serial = k)
  fit$nuisance_ref <- nu$ref
  fit$declaration_ref <- declaration$ref
  class(fit) <- c("cdid_ipw_shared_link_fit", "cdid_ipw_link_fit")
  ref <- new.env(parent = emptyenv())
  fit$link_ref <- ref
  ref$fit <- fit
  lockEnvironment(ref, bindings = TRUE)
  sref$fit_refs[[k]] <- ref
  fit
}

validate_ipw_shared_link_fit <- function(fit, context = NULL) {
  ref <- if (is.list(fit)) fit$link_ref else NULL
  ok <- inherits(fit, "cdid_ipw_shared_link_fit") && is.environment(ref) &&
    environmentIsLocked(ref) && exists("fit", envir = ref, inherits = FALSE) &&
    bindingIsLocked("fit", ref) && identical(fit, get("fit", envir = ref, inherits = FALSE)) &&
    is.environment(fit$declaration_ref) && is.environment(fit$nuisance_ref) &&
    exists("declaration", envir = fit$declaration_ref, inherits = FALSE) &&
    exists("nuisance", envir = fit$nuisance_ref, inherits = FALSE)
  if (ok) {
    decl <- get("declaration", envir = fit$declaration_ref, inherits = FALSE)
    nu <- get("nuisance", envir = fit$nuisance_ref, inherits = FALSE)
    ok <- ipw_seal_ok(decl, "declaration") && ipw_seal_ok(nu, "nuisance") &&
      identical(decl$nuisance_ref, fit$nuisance_ref) &&
      identical(fit$context_ref, nu$context_ref) &&
      identical(fit$pair, decl$pair)
  }
  if (!ok) {
    cdid_abort(
      "Shared link fit is malformed or differs from its frozen snapshot.",
      c("cdid_ipw_shared_fit_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  if (!is.null(context) && !identical(fit$context_ref, context$ref)) {
    cdid_abort("Shared link fit belongs to a different `cdid_ipw_context`.",
               c("cdid_ipw_context_mismatch_error", "cdid_ipw_context_error"))
  }
  invisible(fit)
}

# Recompute every outcome-free nuisance quantity from the frozen context, the
# stored mask and the stored coefficients; require bitwise identity and the
# finite-fit prerequisites. This does not certify glm convergence/asymptotics.
verify_ipw_nuisance_internals <- function(nuisance, context) {
  fail <- function(msg) {
    cdid_abort(paste0("Nuisance integrity check failed: ", msg),
               c("cdid_ipw_nuisance_integrity_error", "cdid_ipw_certificate_error"))
  }
  mask <- nuisance$mask$observed_a_index
  n_a <- context$n_a_units
  if (!is.integer(mask) || length(mask) < 1L || anyNA(mask) || any(mask < 1L | mask > n_a) ||
      (length(mask) > 1L && any(diff(mask) <= 0L))) fail("mask is not a strictly increasing A-unit index.")
  coh <- validate_ipw_cohorts(nuisance$g, nuisance$K)
  if (!identical(nuisance$cohorts, coh$cohorts) ||
      !identical(nuisance$cohort_labels, ipw_cohort_labels(coh)) ||
      !identical(nuisance$control, validate_ipw_glm_control(nuisance$control))) {
    fail("cohort, ordered control or covariate-map metadata is inconsistent.")
  }
  st <- ipw_stage_from_nuisance(nuisance, context)
  if (!identical(st$n, as.double(context$n_a_units)) ||
      !identical(st$m_cohorts, length(coh$cohorts)) || !identical(st$p, context$n_covariates)) {
    fail("nuisance dimensions differ from the frozen context and declared cohorts.")
  }
  chk <- tryCatch({
    sup <- ipw_cohort_support(context, coh)
    gi <- match(context$a_g[mask], nuisance$cohorts)
    w <- context$w_a[mask]
    pr <- ipw_pair_rates(gi, w, sup)
    fw <- w / ipw_clip_floor(pr$rates)[gi]
    lv <- ipw_logit_values(st$Z_use, st$beta)
    ow <- ipw_obs_weights(fw, lv$odds, gi, sup$n)
    st2 <- st
    st2$fit_weights <- fw
    st2$prob_use <- lv$prob
    hs <- ipw_nuisance_hessian_stage(st2, nuisance$cohort_labels)
    list(sup = sup, gi = gi, w = w, pr = pr, fw = fw, lv = lv, ow = ow, hs = hs)
  }, cdid_error = function(e) fail(conditionMessage(e)))
  same <- identical(chk$sup, nuisance$support) && identical(chk$gi, st$risk_gi) &&
    identical(chk$w, st$risk_w) && identical(chk$pr$rates, st$rates) &&
    identical(chk$pr$pair_count, st$pair_count) && identical(chk$pr$pair_mass, st$pair_mass) &&
    identical(chk$fw, st$fit_weights) && identical(chk$lv$odds, st$odds_use) &&
    identical(chk$lv$prob, st$prob_use) && identical(chk$ow$is_t, st$is_t) &&
    identical(chk$ow$a_obs, st$a_obs) && identical(chk$ow$b_obs, st$b_obs) &&
    identical(chk$ow$sum_a, st$sum_a) && identical(chk$ow$sum_b, st$sum_b) &&
    identical(chk$ow$mean_a, st$mean_a) && identical(chk$ow$mean_b, st$mean_b) &&
    identical(chk$hs, nuisance$hessian) &&
    identical(nuisance$mask$n_observed_risk_units, length(mask)) &&
    identical(nuisance$mask$n_risk_a_units, sum(chk$sup$group_a_count)) &&
    identical(nuisance$mask$n_pair_nonobservers, sum(chk$sup$group_a_count) - length(mask))
  if (!same) fail("stored nuisance quantities are not those implied by its context, mask and coefficients.")
  pos <- all(is.finite(st$beta)) && all(st$rates > 0 & st$rates <= 1) &&
    all(st$pair_count >= 1L) && is.finite(st$sum_a) && st$sum_a > 0 &&
    is.finite(st$sum_b) && st$sum_b > 0 && all(is.finite(nuisance$hessian$H_chol)) &&
    all(diag(nuisance$hessian$H_chol) > 0) && all(st$prob_use > 0 & st$prob_use < 1)
  if (!pos) fail("positive rates/denominators or finite algebra prerequisites fail.")
  invisible(TRUE)
}

validate_ipw_certificate_bounds <- function(bounds) {
  req <- c("max_links", "max_classes", "max_class_dates", "max_E_entries", "max_mask_entries")
  if (!is.list(bounds) || is.object(bounds) || length(bounds) != length(req) ||
      is.null(names(bounds)) || anyDuplicated(names(bounds)) || !setequal(names(bounds), req)) {
    cdid_abort(paste0("`bounds` must be a list with exactly: ", paste(req, collapse = ", "), "."),
               c("cdid_ipw_bound_error", "cdid_pairs_bound_error"))
  }
  out <- lapply(req, function(nm) ipw_bound(bounds[[nm]], nm, allow_zero = nm == "max_E_entries"))
  names(out) <- req
  out
}

validate_ipw_design_for_classes <- function(design, n_links, design_link_index) {
  bad <- function(msg) cdid_abort(msg, c("cdid_ipw_design_error", "cdid_ipw_certificate_error"))
  if (!inherits(design, "cdid_chained_design") || !is.list(design) ||
      !is.data.frame(design$links) || !is.data.frame(design$vertices) || !is.list(design$incidence) ||
      anyDuplicated(names(design)) || anyDuplicated(names(design$links)) ||
      anyDuplicated(names(design$vertices)) || anyDuplicated(names(design$incidence))) {
    bad("Expected a `cdid_chained_design` object.")
  }
  L <- design$links
  V <- design$vertices
  need_l <- c("g", "s", "t", "functional_key", "start_vertex_index", "end_vertex_index")
  need_v <- c("g", "time", "role", "functional_key")
  if (!all(need_l %in% names(L)) || !all(need_v %in% names(V)) || nrow(L) != n_links ||
      !identical(design$incidence$n_links, n_links) ||
      !is.numeric(L$g) || !is.numeric(L$s) || !is.numeric(L$t) || !is.character(L$functional_key) ||
      !is.integer(L$start_vertex_index) || !is.integer(L$end_vertex_index) ||
      anyNA(L$start_vertex_index) || anyNA(L$end_vertex_index) ||
      any(L$start_vertex_index < 1L | L$start_vertex_index > nrow(V)) ||
      any(L$end_vertex_index < 1L | L$end_vertex_index > nrow(V)) ||
      !is.numeric(V$g) || !is.numeric(V$time) || !is.character(V$role)) {
    bad("Malformed design links/vertices, or design rows do not match the supplied links one-to-one.")
  }
  has_cal <- !is.null(design$calendar_spec) || !is.null(attr(L, "calendar_spec", exact = TRUE))
  cspec <- NULL
  if (has_cal) {
    if (!is.list(design$calendar_spec) || !inherits(design$calendar_spec, "cdid_calendar_spec")) {
      bad("Missing or malformed `design$calendar_spec` on calendar-adapted design.")
    }
    cspec <- tryCatch(
      validate_calendar_spec(
        design$calendar_spec$calendar,
        design$calendar_spec$cohorts,
        design$calendar_spec$anticipation,
        pt_assumption = if (is.null(design$calendar_spec$pt_assumption)) "post" else design$calendar_spec$pt_assumption,
        pt_start = design$calendar_spec$pt_start
      ),
      cdid_error = function(e) bad(conditionMessage(e))
    )
    if (!identical(design$calendar_spec, cspec) ||
        (!is.null(attr(L, "calendar_spec", exact = TRUE)) && !identical(attr(L, "calendar_spec", exact = TRUE), cspec))) {
      bad("`design$calendar_spec` differs from its canonical validated specification.")
    }
  }
  tryCatch(validate_design_links(L, calendar_spec = cspec), cdid_error = function(e) {
    if (inherits(e, "cdid_design_duplicate_link_error")) {
      cdid_abort(conditionMessage(e), c("cdid_ipw_class_duplicate_error", "cdid_ipw_certificate_error"))
    }
    bad(conditionMessage(e))
  })
  valid_roles <- if (has_cal) {
    c("reference", "pre", "anticipation_nuisance", "anticipation", "post_nuisance", "post")
  } else {
    c("reference", "pre", "post_nuisance", "post")
  }
  if (!is.character(V$functional_key) || anyNA(V$g) || anyNA(V$time) || anyNA(V$role) ||
      any(!is.finite(V$g) | !is.finite(V$time)) ||
      any(!V$role %in% valid_roles)) {
    bad("Malformed intrinsic design vertex coordinates or roles.")
  }
  dli <- design_link_index
  if (!is.atomic(dli) || !is.null(dim(dli)) || is.object(dli) || !is.numeric(dli) ||
      length(dli) != n_links || anyNA(dli) || any(!is.finite(dli) | dli != floor(dli)) ||
      !identical(sort(as.integer(dli)), seq_len(n_links))) {
    bad("`design_link_index` must be a permutation mapping every supplied link to one design row.")
  }
  as.integer(dli)
}

# The actual W is read from incidence, not from the displayed vertex indices.
# Check the compact parameter and incidence maps against those intrinsic indices
# before issuing a range certificate; no dense W or reconstructed graph is needed.
validate_ipw_design_incidence <- function(design) {
  bad <- function(msg) cdid_abort(msg, c("cdid_ipw_class_vertex_error", "cdid_ipw_design_error",
                                        "cdid_ipw_certificate_error"))
  L <- design$links
  V <- design$vertices
  P <- design$parameters
  inc <- design$incidence
  need_v <- c("vertex_index", "param_index", "vertex_id")
  need_p <- c("g", "time", "functional_key", "role", "vertex_index", "param_index", "param_id")
  need_l <- c("start_param_index", "end_param_index", "link_id")
  if (!all(need_v %in% names(V)) || !all(need_l %in% names(L)) ||
      !is.data.frame(P) || anyDuplicated(names(P)) || !all(need_p %in% names(P)) ||
      !is.integer(V$vertex_index) || !identical(V$vertex_index, seq_len(nrow(V))) ||
      !is.integer(V$param_index) || !is.character(V$vertex_id) || anyNA(V$vertex_id) ||
      any(!nzchar(V$vertex_id)) || anyDuplicated(V$vertex_id)) {
    bad("Malformed intrinsic vertex-to-parameter map.")
  }
  v_tuple <- paste(match(V$g, sort(unique(V$g))), match(V$time, sort(unique(V$time))),
                   V$role, ifelse(is.na(V$functional_key), "", paste0("fk:", V$functional_key)), sep = ":")
  if (anyDuplicated(v_tuple)) {
    bad("Duplicate intrinsic vertex coordinates/roles in `design$vertices`.")
  }
  shared_v <- V$role %in% c("reference", "anticipation", "post")
  if (any(!is.na(V$functional_key[shared_v])) ||
      any(is.na(V$functional_key[!shared_v]) | !nzchar(V$functional_key[!shared_v]))) {
    bad("Vertex `functional_key` disagrees with its shared/private role.")
  }
  sv <- L$start_vertex_index; tv <- L$end_vertex_index
  if (any(V$g[sv] != L$g | V$g[tv] != L$g | V$time[sv] != L$s | V$time[tv] != L$t |
          (!shared_v[sv] & V$functional_key[sv] != L$functional_key) |
          (!shared_v[tv] & V$functional_key[tv] != L$functional_key))) {
    bad("Link endpoints differ from their intrinsic cohort, date or operator vertices.")
  }
  if (!is.null(design$calendar_spec)) {
    cspec <- design$calendar_spec
    if (!inherits(cspec, "cdid_calendar_spec") || !is.list(cspec)) {
      bad("Design `calendar_spec` is not a valid `cdid_calendar_spec` object.")
    }
    cspec_val <- tryCatch(
      validate_calendar_spec(
        cspec$calendar, cspec$cohorts, cspec$anticipation,
        pt_assumption = if (is.null(cspec$pt_assumption)) "post" else cspec$pt_assumption,
        pt_start = cspec$pt_start
      ),
      cdid_error = function(e) bad(conditionMessage(e))
    )
    if (!identical(cspec, cspec_val)) {
      bad("Design `calendar_spec` differs from its canonical validated representation.")
    }
    target_cohorts <- sort(unique(L$g))
    used_vertices <- unique(c(L$start_vertex_index, L$end_vertex_index))
    unlinked_vertices <- setdiff(seq_len(nrow(V)), used_vertices)
    if (any(!V$g %in% target_cohorts) || any(!V$g %in% cspec$cohorts) ||
        any(!V$time %in% cspec$calendar)) {
      bad("Design vertex cohort or date lies outside `design$calendar_spec`.")
    }
    if (any(!V$role[unlinked_vertices] %in% c("reference", "post"))) {
      bad("An unlinked vertex must be a cohort reference or requested post target.")
    }
    reference_vertices <- which(V$role == "reference")
    if (length(reference_vertices) != length(target_cohorts) ||
        anyDuplicated(V$g[reference_vertices]) ||
        !setequal(V$g[reference_vertices], target_cohorts)) {
      bad("Every declared target cohort requires its unique intrinsic reference vertex.")
    }
    post_vertices <- which(V$role == "post")
    targets <- design$post_targets
    if (!is.data.frame(targets) || anyDuplicated(names(targets)) ||
        !all(c("g", "time", "param_id", "param_index") %in% names(targets)) ||
        nrow(targets) != length(post_vertices) ||
        !identical(targets$g, V$g[post_vertices]) ||
        !identical(targets$time, V$time[post_vertices]) ||
        !identical(targets$param_id, V$vertex_id[post_vertices]) ||
        !identical(targets$param_index, V$param_index[post_vertices])) {
      bad("Post target rows differ from their intrinsic causal vertices.")
    }
    for (v in seq_len(nrow(V))) {
      r_v <- tryCatch(calendar_reference_date(V$g[v], cspec$calendar, cspec$anticipation),
                      cdid_error = function(e) bad(conditionMessage(e)))
      cmp_v <- tryCatch(exact_cmp_sum(V$time[v], cspec$anticipation, V$g[v]),
                        cdid_error = function(e) bad(conditionMessage(e)))
      rv <- V$role[v]
      ok_v <- if (V$time[v] == r_v) {
        identical(rv, "reference")
      } else if (cmp_v < 0L) {
        identical(rv, "pre")
      } else if (V$time[v] < V$g[v]) {
        rv %in% c("anticipation", "anticipation_nuisance")
      } else {
        rv %in% c("post", "post_nuisance")
      }
      if (!ok_v) bad("Design vertex role disagrees with its calendar date and anticipation window.")
      # Isolated requested post targets are desired coordinates, not qualified
      # operator endpoints. Their lack of eligible controls must remain visible.
      if (rv %in% c("anticipation", "post") && v %in% used_vertices) {
        cmp_max <- tryCatch(exact_cmp_sum(V$time[v], cspec$anticipation, cspec$max_cohort),
                            cdid_error = function(e) bad(conditionMessage(e)))
        if (cmp_max >= 0L) bad("Qualified causal vertex violates the latest-cohort clean-control boundary.")
      }
    }
  }
  is_pt_all <- !is.null(design$calendar_spec) && identical(design$calendar_spec$pt_assumption, "all")
  is_anchor_v <- if (is_pt_all) {
    V$role == "reference" | (V$role == "pre" & V$time >= design$calendar_spec$pt_start)
  } else {
    V$role == "reference"
  }
  if ("anchored_value" %in% names(V)) {
    expected_av <- ifelse(is_anchor_v, 0, NA_real_)
    if (!identical(as.double(V$anchored_value), as.double(expected_av))) {
      bad("Vertex `anchored_value` disagrees with its reference or PT-all anchoring status.")
    }
  }
  if (!is.null(design$calendar_spec)) {
    pt_start_val <- if (is_pt_all) design$calendar_spec$pt_start else Inf
    for (side in 1:2) {
      rc_col <- if (side == 1L) "s_coordinate_role" else "t_coordinate_role"
      if (rc_col %in% names(L)) {
        v_idx <- if (side == 1L) L$start_vertex_index else L$end_vertex_index
        rv <- V$role[v_idx]
        uv <- V$time[v_idx]
        rc_val <- L[[rc_col]]
        ok_rc <- (rv == "reference" & rc_val == "reference") |
          (rv == "pre" & rc_val %in% c("pre", "pre_nuisance")) |
          (rv == "pre" & is_pt_all & uv >= pt_start_val & rc_val == "pre_pt_all_zero") |
          (rv == "anticipation" & rc_val %in% c("anticipation", "anticipation_att")) |
          (rv == "anticipation_nuisance" & rc_val == "anticipation_nuisance") |
          (rv == "post" & rc_val %in% c("post", "post_att")) |
          (rv == "post_nuisance" & rc_val == "post_nuisance")
        if (!all(ok_rc)) bad("Design link endpoint coordinate role disagrees with its intrinsic vertex role or PT-all window.")
      }
    }
  }
  pv <- which(!is_anchor_v)
  vp <- rep(NA_integer_, nrow(V))
  vp[pv] <- seq_along(pv)
  if (!identical(V$param_index, vp) || nrow(P) != length(pv) ||
      !identical(P$param_id, V$vertex_id[pv]) ||
      !all(vapply(c("g", "time", "functional_key", "role", "vertex_index", "param_index"),
                  function(nm) identical(P[[nm]], V[[nm]][pv]), logical(1L)))) {
    bad("Parameter rows differ from the intrinsic qualified vertices.")
  }
  sp <- vp[L$start_vertex_index]
  ep <- vp[L$end_vertex_index]
  if (!identical(L$start_param_index, sp) || !identical(L$end_param_index, ep) ||
      !identical(inc$n_parameters, length(pv)) ||
      !identical(inc$start, ifelse(is.na(sp), 0L, sp)) ||
      !identical(inc$end, ifelse(is.na(ep), 0L, ep)) ||
      !identical(inc$link_ids, L$link_id) || !identical(inc$param_ids, P$param_id)) {
    bad("Actual incidence differs from the intrinsic link and parameter maps.")
  }
  invisible(TRUE)
}

# Signed forest path s -> t: mu_(s,t) = sum sign * mu_e over forest edges e.
ipw_forest_path <- function(a, b, forest, depth, count_only) {
  k <- 0
  edges <- signs <- integer(0L)
  if (!count_only) {
    len <- depth[a] + depth[b]
    edges <- signs <- integer(len)
  }
  x <- a
  y <- b
  while (x != y) {
    k <- k + 1
    if (depth[x] >= depth[y]) {
      if (!count_only) {
        edges[k] <- forest$parent_link[x]
        signs[k] <- -forest$parent_sign[x]
      }
      x <- forest$parent_vertex[x]
    } else {
      if (!count_only) {
        edges[k] <- forest$parent_link[y]
        signs[k] <- forest$parent_sign[y]
      }
      y <- forest$parent_vertex[y]
    }
    if (x == 0L || y == 0L) {
      cdid_abort("Forest path left its component.", c("cdid_ipw_certificate_error"))
    }
  }
  if (count_only) return(k)
  list(edges = edges[seq_len(k)], signs = signs[seq_len(k)])
}

ipw_forest_depth <- function(forest) {
  nv <- length(forest$parent_vertex)
  depth <- rep(-1L, nv)
  depth[forest$parent_vertex == 0L] <- 0L
  for (v in seq_len(nv)) {
    path <- integer(0L)
    cur <- v
    while (depth[cur] < 0L) {
      path <- c(cur, path)
      cur <- forest$parent_vertex[cur]
    }
    d <- depth[cur]
    for (u in path) {
      d <- d + 1L
      depth[u] <- d
    }
  }
  depth
}

# Structural class certificate (Prop. 5 / Cor. 6). `links` holds shared fits,
# shared declarations (unfitted, outcome-free) and/or ordinary fits (unclassed,
# always retained). Dropped rows are never evaluated here.
certify_ipw_structural_classes <- function(session, links, design, design_link_index,
                                           weighting, bounds) {
  sref <- validate_ipw_shared_session(session)
  context <- sref$context
  if (!is.character(weighting) || length(weighting) != 1L || is.na(weighting) ||
      !weighting %in% c("optimal", "identity")) {
    cdid_abort("`weighting` must be exactly \"optimal\" or \"identity\".", c("cdid_ipw_input_error"))
  }
  if (!identical(weighting, "optimal")) {
    cdid_abort(
      "Structural reduction is licensed only for optimal weighting; identity weighting keeps every row.",
      c("cdid_ipw_reduction_weighting_error", "cdid_ipw_certificate_error")
    )
  }
  b <- validate_ipw_certificate_bounds(bounds)
  if (!is.list(links) || is.object(links) || length(links) < 1L) {
    cdid_abort("`links` must be a nonempty plain list.", c("cdid_ipw_input_error", "cdid_ipw_certificate_error"))
  }
  n_links <- length(links)
  if (n_links > b$max_links) {
    cdid_abort("Number of links exceeds `max_links`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  # 1. Actual object identity: sealed, session-issued, same frozen context.
  kind <- integer(n_links)
  decls <- vector("list", n_links)
  for (ell in seq_len(n_links)) {
    obj <- links[[ell]]
    if (inherits(obj, "cdid_ipw_shared_link_fit")) {
      validate_ipw_shared_link_fit(obj, context)
      decl <- get("declaration", envir = obj$declaration_ref, inherits = FALSE)
      validate_ipw_shared_declaration(decl, sref)
      if (!identical(sref$fit_refs[[decl$serial]], obj$link_ref)) {
        cdid_abort("Shared link fit was not issued by this session.",
                   c("cdid_ipw_shared_fit_error", "cdid_ipw_certificate_error"))
      }
      kind[ell] <- 2L
      decls[[ell]] <- decl
    } else if (inherits(obj, "cdid_ipw_shared_link_declaration")) {
      validate_ipw_shared_declaration(obj, sref)
      kind[ell] <- 1L
      decls[[ell]] <- obj
    } else if (inherits(obj, "cdid_ipw_link_fit") && is.list(obj)) {
      if (!identical(obj$context_ref, context$ref)) {
        cdid_abort("Ordinary link fit belongs to a different `cdid_ipw_context`.",
                   c("cdid_ipw_context_mismatch_error", "cdid_ipw_certificate_error"))
      }
      kind[ell] <- 0L
    } else {
      cdid_abort("Every element of `links` must be a shared fit, shared declaration or ordinary fit.",
                 c("cdid_ipw_input_error", "cdid_ipw_certificate_error"))
    }
  }
  shared <- which(kind > 0L)
  mask_work <- as.double(n_links) * as.double(context$n_a_units)
  if (mask_work > b$max_mask_entries) {
    cdid_abort("Mask recomputation exceeds `max_mask_entries`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  link_g <- link_s <- link_t <- numeric(n_links)
  link_fk <- vector("list", n_links)
  link_K <- vector("list", n_links)
  for (ell in seq_len(n_links)) {
    if (kind[ell] > 0L) {
      link_g[ell] <- decls[[ell]]$g
      link_s[ell] <- decls[[ell]]$pair$s
      link_t[ell] <- decls[[ell]]$pair$t
      link_fk[ell] <- list(decls[[ell]]$pair$functional_key)
      link_K[[ell]] <- as.double(decls[[ell]]$K)
    } else {
      link_g[ell] <- links[[ell]]$cohorts$g
      link_s[ell] <- links[[ell]]$pair$s
      link_t[ell] <- links[[ell]]$pair$t
      link_fk[ell] <- list(links[[ell]]$pair$functional_key)
      link_K[[ell]] <- as.double(links[[ell]]$cohorts$K)
    }
  }

  # 2. Intrinsic per-date vertex map from the design (never from labels).
  dli <- validate_ipw_design_for_classes(design, n_links, design_link_index)
  DL <- design$links
  DV <- design$vertices
  cspec <- design$calendar_spec
  has_cal <- !is.null(cspec)
  is_pt_all <- FALSE
  vs <- DL$start_vertex_index[dli]
  vt <- DL$end_vertex_index[dli]
  fk <- DL$functional_key[dli]
  map_bad <- function(msg) cdid_abort(msg, c("cdid_ipw_class_vertex_error", "cdid_ipw_certificate_error"))
  if (has_cal) {
    if (!inherits(cspec, "cdid_calendar_spec") || !is.list(cspec)) {
      map_bad("Design `calendar_spec` is not a valid `cdid_calendar_spec` object.")
    }
    cspec_val <- tryCatch(
      validate_calendar_spec(
        cspec$calendar, cspec$cohorts, cspec$anticipation,
        pt_assumption = if (is.null(cspec$pt_assumption)) "post" else cspec$pt_assumption,
        pt_start = cspec$pt_start
      ),
      cdid_error = function(e) map_bad(conditionMessage(e))
    )
    if (!identical(cspec, cspec_val)) {
      map_bad("Design `calendar_spec` differs from its canonical validated representation.")
    }
    is_pt_all <- identical(cspec$pt_assumption, "all")
  }
  # One typed parse per distinct design pool string (CP08 scaling action 1).
  dl_pools <- if (has_cal && "control_pool" %in% names(DL)) design_pool_cache(DL$control_pool) else NULL
  for (ell in seq_len(n_links)) {
    j <- dli[ell]
    if (!identical(as.double(DL$g[j]), link_g[ell]) || !identical(as.double(DL$s[j]), link_s[ell]) ||
        !identical(as.double(DL$t[j]), link_t[ell])) {
      map_bad("Design row cohort/dates differ from the link they are mapped to.")
    }
    if (has_cal) {
      delta <- cspec$anticipation
      K_ell <- link_K[[ell]]
      ref_clean <- tryCatch(
        all(vapply(K_ell, function(k) exact_cmp_sum(DL$ref[j], delta, k) < 0L, logical(1L))),
        cdid_error = function(e) map_bad(conditionMessage(e))
      )
      if (any(!K_ell %in% cspec$cohorts) || any(K_ell == link_g[ell]) || !ref_clean) {
        map_bad("Link control pool K is outside the design cohort universe or not clean at the cohort reference period.")
      }
      if ("control_pool" %in% names(DL)) {
        cp_vals <- dl_pools$values(j)
        if (is.null(cp_vals) || !setequal(as.double(cp_vals), K_ell)) {
          map_bad("Design row `control_pool` differs from the link's actual control pool K.")
        }
      }
    }
    for (side in 1:2) {
      v <- if (side == 1L) vs[ell] else vt[ell]
      u <- if (side == 1L) link_s[ell] else link_t[ell]
      qual_col <- if (side == 1L) "s_qualified" else "t_qualified"
      role_col <- if (side == 1L) "s_coordinate_role" else "t_coordinate_role"
      if (!has_cal) {
        qualified <- if (qual_col %in% names(DL)) DL[[qual_col]][j] else (u == DL$ref[j] || u >= link_g[ell])
        role <- if (u == DL$ref[j]) "reference" else if (u < link_g[ell]) "pre" else if (qualified) "post" else "post_nuisance"
        key <- if (role %in% c("pre", "post_nuisance")) fk[ell] else NA_character_
        expected_anchored <- (u == DL$ref[j])
      } else {
        cmp_ug <- tryCatch(exact_cmp_sum(u, delta, link_g[ell]),
                           cdid_error = function(e) map_bad(conditionMessage(e)))
        qualified <- if (qual_col %in% names(DL)) DL[[qual_col]][j] else (u == DL$ref[j] || cmp_ug >= 0L)
        if (cmp_ug >= 0L && qualified) {
          u_clean <- tryCatch(
            exact_cmp_sum(u, delta, cspec$max_cohort) < 0L &&
              all(vapply(K_ell, function(k) exact_cmp_sum(u, delta, k) < 0L, logical(1L))),
            cdid_error = function(e) map_bad(conditionMessage(e))
          )
          if (!u_clean) {
            map_bad("Qualified causal endpoint is not clean under the link's actual control pool K and calendar anticipation.")
          }
        }
        role <- if (u == DL$ref[j]) {
          "reference"
        } else if (cmp_ug < 0L) {
          "pre"
        } else if (u < link_g[ell]) {
          if (qualified) "anticipation" else "anticipation_nuisance"
        } else {
          if (qualified) "post" else "post_nuisance"
        }
        key <- if (role %in% c("pre", "anticipation_nuisance", "post_nuisance")) fk[ell] else NA_character_
        expected_anchored <- (u == DL$ref[j]) || (is_pt_all && cmp_ug < 0L && u >= cspec$pt_start)
        if (role_col %in% names(DL)) {
          rc_val <- DL[[role_col]][j]
          ok_rc <- if (u == DL$ref[j]) {
            identical(rc_val, "reference")
          } else if (cmp_ug < 0L) {
            rc_val %in% c("pre", "pre_nuisance") || (is_pt_all && u >= cspec$pt_start && identical(rc_val, "pre_pt_all_zero"))
          } else if (u < link_g[ell]) {
            if (qualified) rc_val %in% c("anticipation", "anticipation_att") else identical(rc_val, "anticipation_nuisance")
          } else {
            if (qualified) rc_val %in% c("post", "post_att") else identical(rc_val, "post_nuisance")
          }
          if (!isTRUE(ok_rc)) {
            map_bad("Design row endpoint coordinate role disagrees with its calendar region or PT-all window.")
          }
        }
      }
      if (!identical(as.double(DV$time[v]), u) || !identical(as.double(DV$g[v]), link_g[ell]) ||
          !identical(DV$role[v], role) || !identical(DV$functional_key[v], key) ||
          !identical(is.na(DV$param_index[v]), expected_anchored)) {
        map_bad("Endpoint vertex is not the intrinsic qualified vertex of the link's cohort, date and operator.")
      }
    }
    if (vs[ell] == vt[ell]) map_bad("Loop link in design.")
    if (!is.null(link_fk[[ell]]) && !identical(link_fk[[ell]], fk[ell])) {
      cdid_abort("Link functional key differs from its design operator key.",
                 c("cdid_ipw_class_operator_error", "cdid_ipw_certificate_error"))
    }
  }

  validate_ipw_design_incidence(design)

  # 3. Recompute every shared mask and nuisance internal from the frozen context.
  nu_refs <- list()
  nu_gid <- integer(n_links)
  for (ell in shared) {
    r <- decls[[ell]]$nuisance_ref
    hit <- 0L
    for (k in seq_along(nu_refs)) if (identical(nu_refs[[k]], r)) { hit <- k; break }
    if (hit == 0L) {
      nu_refs[[length(nu_refs) + 1L]] <- r
      hit <- length(nu_refs)
      verify_ipw_nuisance_internals(get("nuisance", envir = r, inherits = FALSE), context)
    }
    nu_gid[ell] <- hit
  }
  for (ell in shared) {
    nu <- get("nuisance", envir = nu_refs[[nu_gid[ell]]], inherits = FALSE)
    chk <- ipw_shared_mask_check(context, nu, decls[[ell]]$pair)
    if (!identical(chk$n_obs_pair_all, decls[[ell]]$n_observed_pair_units)) {
      cdid_abort("Declared pair count differs from recomputed postings.",
                 c("cdid_ipw_declaration_error", "cdid_ipw_certificate_error"))
    }
  }

  # 4. Classes: (verified nuisance object, design operator key). Bounded.
  fk_levels <- sort(unique(fk), method = "radix")
  fk_id <- match(fk, fk_levels)
  ckey_n <- nu_gid[shared]
  ckey_f <- fk_id[shared]
  ord <- order(ckey_n, ckey_f, method = "radix")
  starts <- if (length(shared)) which(c(TRUE, diff(ckey_n[ord]) != 0L | diff(ckey_f[ord]) != 0L)) else integer(0L)
  n_classes <- length(starts)
  if (n_classes > b$max_classes) {
    cdid_abort("Number of structural classes exceeds `max_classes`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }
  ends <- c(starts[-1L] - 1L, length(shared))
  plans <- vector("list", n_classes)
  total_E <- 0
  for (ci in seq_len(n_classes)) {
    members <- shared[ord[seq.int(starts[ci], ends[ci])]]
    members <- members[order(link_s[members], link_t[members], method = "radix")]
    dates <- sort(unique(c(link_s[members], link_t[members])))
    if (length(dates) > b$max_class_dates) {
      cdid_abort("Class date count exceeds `max_class_dates`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
    }
    a_idx <- match(link_s[members], dates)
    b_idx <- match(link_t[members], dates)
    if (anyDuplicated(paste(a_idx, b_idx))) {
      cdid_abort("Duplicate (s, t) link within one structural class.",
                 c("cdid_ipw_class_duplicate_error", "cdid_ipw_certificate_error"))
    }
    vmap <- rep(NA_integer_, length(dates))
    for (side in 1:2) {
      di <- if (side == 1L) a_idx else b_idx
      vv <- if (side == 1L) vs[members] else vt[members]
      for (h in seq_along(di)) {
        if (is.na(vmap[di[h]])) vmap[di[h]] <- vv[h] else if (vmap[di[h]] != vv[h]) {
          map_bad("Class links map one date to different qualified vertices; no common vertex map.")
        }
      }
    }
    if (anyDuplicated(vmap)) map_bad("Class vertex map is not injective over dates.")
    forest <- incidence_forest(length(dates), a_idx, b_idx, integer(0L))
    tree <- sort(forest$parent_link[forest$parent_link > 0L])
    dropped <- setdiff(seq_along(members), tree)
    depth <- ipw_forest_depth(forest)
    n_e <- 0
    for (h in dropped) n_e <- n_e + ipw_forest_path(a_idx[h], b_idx[h], forest, depth, TRUE)
    total_E <- total_E + n_e
    plans[[ci]] <- list(members = members, dates = dates, a_idx = a_idx, b_idx = b_idx,
                        vmap = vmap, forest = forest, depth = depth, tree = tree, dropped = dropped)
  }
  if (total_E > b$max_E_entries) {
    cdid_abort("Reduction map entries exceed `max_E_entries`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
  }

  # 5. Retained rows built from the same exact nuisance inputs (context, g,
  # ordered K, control, risk-set mask) -- one object, duplicate objects or
  # ordinary refits -- satisfy exact IF relations along date cycles. Such rows
  # must stay acyclic across classes; otherwise the singularity is unexplained.
  src_keys <- list()
  src_id <- integer(n_links)
  for (ell in seq_len(n_links)) {
    if (kind[ell] > 0L) {
      nu <- get("nuisance", envir = nu_refs[[nu_gid[ell]]], inherits = FALSE)
      key <- list(nu$g, nu$K, nu$control, nu$mask$observed_a_index)
    } else {
      f <- links[[ell]]
      coh <- validate_ipw_cohorts(f$cohorts$g, f$cohorts$K)
      desc <- make_date_pair_descriptor(f$pair$s, f$pair$t, context)
      hit <- query_pair_mask_internal(context, desc$s, desc$t,
                                      ds = desc$s_calendar_index, dt = desc$t_calendar_index)
      gi <- match(context$a_g[hit$a_index], coh$cohorts)
      key <- list(coh$g, coh$K, validate_ipw_glm_control(f$diagnostics$control), hit$a_index[!is.na(gi)])
    }
    hit_k <- 0L
    for (k in seq_along(src_keys)) if (identical(src_keys[[k]], key)) { hit_k <- k; break }
    if (hit_k == 0L) {
      src_keys[[length(src_keys) + 1L]] <- key
      hit_k <- length(src_keys)
    }
    src_id[ell] <- hit_k
  }
  dropped_all <- unlist(lapply(plans, function(pl) pl$members[pl$dropped]), use.names = FALSE)
  retained_all <- setdiff(seq_len(n_links), dropped_all)
  for (k in seq_along(src_keys)) {
    rk <- retained_all[src_id[retained_all] == k]
    if (length(rk) < 2L) next
    dk <- sort(unique(c(link_s[rk], link_t[rk])))
    parent <- seq_along(dk)
    find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
    for (ell in rk) {
      x <- find(match(link_s[ell], dk))
      y <- find(match(link_t[ell], dk))
      if (x == y) {
        cdid_abort(
          "Retained links with identical nuisance inputs (same object across operators, duplicate objects or ordinary refits) form a date cycle; this singularity is not certified by any class.",
          c("cdid_ipw_reduction_singular_error", "cdid_ipw_certificate_error")
        )
      }
      parent[x] <- y
    }
  }

  # 6. Sparse E for dropped rows only (allocation bounded above).
  E_i <- E_j <- E_x <- integer(total_E)
  pos <- 0L
  classes <- vector("list", n_classes)
  for (ci in seq_len(n_classes)) {
    pl <- plans[[ci]]
    for (h in pl$dropped) {
      path <- ipw_forest_path(pl$a_idx[h], pl$b_idx[h], pl$forest, pl$depth, FALSE)
      o <- order(pl$members[path$edges])
      idx <- pos + seq_along(o)
      E_i[idx] <- pl$members[h]
      E_j[idx] <- pl$members[path$edges[o]]
      E_x[idx] <- path$signs[o]
      pos <- pos + length(o)
    }
    first <- pl$members[1L]
    nu <- get("nuisance", envir = nu_refs[[nu_gid[first]]], inherits = FALSE)
    classes[[ci]] <- list(
      nuisance_serial = nu$serial,
      functional_key = fk[first],
      g = nu$g,
      K = nu$K,
      dates = pl$dates,
      vertex_index = pl$vmap,
      link_index = pl$members,
      forest_link_index = pl$members[pl$tree],
      dropped_link_index = pl$members[pl$dropped],
      roots = pl$dates[pl$forest$roots],
      n_components = pl$forest$n_components,
      cycle_rank = length(pl$members) - (length(pl$dates) - pl$forest$n_components),
      n_risk_a_units = nu$mask$n_risk_a_units,
      n_observers = nu$mask$n_observed_risk_units,
      n_pair_nonobservers = nu$mask$n_pair_nonobservers,
      n_outside_risk_a_units = context$n_a_units - nu$mask$n_risk_a_units
    )
  }
  ord_e <- order(E_i, E_j, method = "radix")
  dropped_index <- sort(dropped_all)
  certificate <- structure(
    list(
      weighting = "optimal",
      n_links = n_links,
      link_kind = c("ordinary_unclassed", "shared_declaration", "shared_fit")[kind + 1L],
      retained_index = setdiff(seq_len(n_links), dropped_index),
      dropped_index = dropped_index,
      E = list(i = E_i[ord_e], j = E_j[ord_e], x = E_x[ord_e]),
      classes = classes,
      design_link_index = dli,
      nuisance_refs = nu_refs,
      object_refs = lapply(seq_len(n_links), function(ell) {
        if (kind[ell] == 2L) links[[ell]]$link_ref else if (kind[ell] == 1L) links[[ell]]$ref else NULL
      }),
      context_ref = context$ref,
      bounds = b,
      retained_covariance = "not_certified: PD/accuracy of the retained covariance must be validated by the caller (Cor. 6); otherwise refuse",
      scope = "exact real-algebra structural identities of Prop. 5 / Cor. 6 for optimal GMM; no numerical acceptance, sampling or catalog claim"
    ),
    class = "cdid_ipw_structural_certificate"
  )
  ipw_seal(certificate, "certificate")
}

# Re-derive the certificate from the actual objects; any difference refuses.
verify_ipw_structural_certificate <- function(certificate, session, links, design, design_link_index) {
  if (!inherits(certificate, "cdid_ipw_structural_certificate") || !is.list(certificate) ||
      !is.list(certificate$bounds)) {
    cdid_abort("Expected `cdid_ipw_structural_certificate`.", c("cdid_ipw_certificate_error", "cdid_ipw_input_error"))
  }
  again <- certify_ipw_structural_classes(session, links, design, design_link_index,
                                          "optimal", certificate$bounds)
  if (!ipw_seal_ok(certificate, "certificate") ||
      !identical(certificate[names(certificate) != "ref"], again[names(again) != "ref"])) {
    cdid_abort("Certificate does not match the actual objects, masks and vertex map.",
               c("cdid_ipw_certificate_mismatch_error", "cdid_ipw_certificate_error"))
  }
  invisible(TRUE)
}

# Expand retained rows to all rows: x = E x_R (vector or row-matrix input).
apply_ipw_reduction_map <- function(certificate, x_retained) {
  if (!inherits(certificate, "cdid_ipw_structural_certificate") || !is.list(certificate) ||
      !ipw_seal_ok(certificate, "certificate")) {
    cdid_abort("Expected an unchanged issued `cdid_ipw_structural_certificate`.",
               c("cdid_ipw_certificate_error", "cdid_ipw_input_error"))
  }
  vec <- is.null(dim(x_retained))
  X <- if (vec) matrix(x_retained, ncol = 1L) else x_retained
  n_r <- length(certificate$retained_index)
  if (!is.numeric(X) || !is.matrix(X) || nrow(X) != n_r || any(!is.finite(X))) {
    cdid_abort("`x_retained` must be finite with one row per retained link.", c("cdid_ipw_input_error"))
  }
  out <- matrix(0, certificate$n_links, ncol(X))
  out[certificate$retained_index, ] <- X
  rpos <- match(certificate$E$j, certificate$retained_index)
  for (k in seq_along(certificate$E$i)) {
    i <- certificate$E$i[k]
    out[i, ] <- out[i, ] + certificate$E$x[k] * X[rpos[k], ]
  }
  if (vec) drop(out) else out
}
