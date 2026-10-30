# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 IPW native batching on the shared fixed-X design: the trusted batch route
# (eligibility, compact (4 + p + m) x L descriptors, the 18-argument kernel's
# p-term loading product, the generalized charges) for p >= 2 propensity
# columns (p counts the intercept; `~x1` is p = 2), with the p = 1 route
# unchanged. Expected values come from the unchanged general route of the same
# installed build (eligibility declined in a re-bound namespace copy), from the
# raw owner arrays and from the four reported IF components; the native kernel
# never supplies its own truth. Floating outputs compare at the established
# 1e-12 tolerance; identities, orders, shapes and statuses compare exactly.
#
# UNEXECUTED at authoring time (no R in the authoring job). Fixtures are
# private to this file; the helpers mirror those of
# test-cp08-native-batch.R (nb_expect_match, nb_reseal) and
# test-cp08-ipw-construction.R (cc_counted, cc_public_ctrl, cc_strip_env).
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Private fixtures (deterministic, no RNG)
# ---------------------------------------------------------------------------

# Balanced long panel: one target cohort g = 4 and the never-treated pool, two
# nonconstant unit-level covariates, constant and unit-varying weights.
# `extra` moves a block of unit positions to an unreferenced cohort 6 (off the
# nuisance's cohorts); `drop` removes "id t" rows.
xb_panel <- function(N = 41L, T = 6L, extra = NULL, drop = NULL) {
  u <- seq_len(N)
  x1 <- (u %% 5L) / 4 - 0.5
  x2 <- ((3L * u) %% 7L) / 6 - 0.5
  g <- ifelse(((u %% 3L) == 0L) | (x1 > 0 & (u %% 2L) == 0L), 4, Inf)
  if (!is.null(extra)) g[extra] <- 6
  d <- expand.grid(t = seq_len(T), id = u, KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$x1 <- x1[d$id]
  d$x2 <- x2[d$id]
  d$y <- 0.3 * d$id + 0.2 * d$t + ifelse(d$t >= d$g, 1, 0) + 0.5 * d$x1 * d$t - 0.25 * d$x2 +
    ((7 * d$id + 3 * d$t) %% 5) / 10
  d$w_const <- 1
  d$w_vary <- 1 + (d$id %% 4L) / 8
  if (!is.null(drop)) d <- d[!(paste(d$id, d$t) %in% drop), ]
  d
}
xb_formula <- function(p) switch(as.character(p), "1" = ~1, "2" = ~x1, "3" = ~x1 + x2,
  stop("p must be 1, 2 or 3"))
# Public pipeline control for `n_covariates = p` (same construction as cc_public_ctrl).
xb_ctrl <- function(d, p, chunk = 7L, weighting = "optimal") {
  gs <- unique(d$g)
  n_frame <- length(unique(d$id))
  dims <- api$chained_public_dimensions_from_counts(n_frame = n_frame, n_panel_units = n_frame,
    n_rows = nrow(d), n_dates = length(unique(d$t)), n_cohorts = length(gs), n_targets = sum(is.finite(gs)),
    n_covariates = p, biters = 99, n_clusters = n_frame)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(list(chunk_rows = chunk),
    gmm_weighting = weighting), gmm_weighting = weighting)$control
}
xb_run <- function(ns, d, ctrl, p, ...) ns$run_chained_ipw_pipeline(
  data = d, idname = "id", tname = "t", yname = "y", gname = "g", control = ctrl, xformla = xb_formula(p),
  control_group = "nevertreated", base_period = "varying", weighting = "optimal", run_numerical_transfer = FALSE,
  root_witness = NULL, structural_reduction_stage = "outcome_free", ...)

# Re-bound copy of the namespace (as cc_counted): native resolution stays in
# the genuine namespace; `edits` replace internal bindings; `wrap` names are
# counted on every call.
xb_ns <- function(edits = list(), wrap = character()) {
  e <- new.env(parent = api)
  counts <- new.env(parent = emptyenv())
  genuine <- "gmm_require_native_dll"
  for (nm in ls(api, all.names = TRUE)) {
    f <- get(nm, envir = api, inherits = FALSE)
    if (is.function(f) && !is.primitive(f) && identical(environment(f), api) && !nm %in% genuine) {
      environment(f) <- e
      assign(nm, f, envir = e)
    }
  }
  for (nm in genuine) assign(nm, get(nm, envir = api, inherits = FALSE), envir = e)
  for (nm in names(edits)) assign(nm, edits[[nm]], envir = e)
  make_wrap <- function(nm, inner) {
    force(nm); force(inner)
    function(...) {
      assign(nm, get(nm, envir = counts, inherits = FALSE) + 1L, envir = counts)
      inner(...)
    }
  }
  for (nm in wrap) {
    if (!exists(nm, envir = e, inherits = FALSE)) stop("No internal function `", nm, "` to count.", call. = FALSE)
    assign(nm, 0L, envir = counts)
    assign(nm, make_wrap(nm, get(nm, envir = e, inherits = FALSE)), envir = e)
  }
  e$.counts <- counts
  e
}
xb_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0L)
# The general route of the same build: the structural eligibility declines.
xb_general <- function() xb_ns(list(ipw_pipeline_batch_eligibility = function(...) NULL))

xb_strip_env <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.list(x)) { a <- attributes(x); y <- lapply(x, xb_strip_env); attributes(y) <- a; return(y) }
  x
}
# Same strict schema/leaf predicates as nb_expect_match, one reported assertion.
xb_expect_match <- function(a, b, tol = 1e-12, info = NULL) {
  mismatch <- function(x, y, path) {
    if (!identical(class(x), class(y)) || !identical(typeof(x), typeof(y)) ||
        !identical(names(x), names(y)) || !identical(dim(x), dim(y)) ||
        !identical(dimnames(x), dimnames(y))) return(paste(path, "schema"))
    if (is.list(x)) {
      if (length(x) != length(y)) return(paste(path, "length"))
      return(unlist(lapply(seq_along(x), function(i) mismatch(x[[i]], y[[i]],
        paste0(path, "$", if (is.null(names(x))) i else names(x)[i]))), use.names = FALSE))
    }
    same <- if (typeof(x) == "double") length(waldo::compare(x, y, tolerance = tol)) == 0L else identical(x, y)
    if (same) character() else paste(path, "value")
  }
  bad <- mismatch(a, b, "root")
  expect_true(length(bad) == 0L, info = paste(c(info, head(bad, 12L)), collapse = "\n"))
}
xb_reseal <- function(x, binding, edit) { x$ref <- NULL; x <- edit(x); api$ipw_owner_seal(x, binding) }
xb_point_fields <- c("link_id", "estimate", "means", "coefficients", "loadings", "partials", "masses", "counts")
xb_record <- function(pl) xb_strip_env(pl[intersect(xb_point_fields, names(pl))])
xb_fields <- c("support", "basis", "aggregations", "gmm_fit", "structural_certificate", "design_all", "design_used")
# Collect every replayed batch of the plan or of the general view.
xb_collect <- function(fn, ..., include_components = TRUE) {
  out <- list()
  fn(..., callback = function(b) { force(b); out[[length(out) + 1L]] <<- b }, include_components = include_components)
  out
}
xb_compare_batches <- function(bb, gb, tol = 1e-12, info = NULL) {
  expect_identical(length(bb), length(gb), info = info)
  for (i in seq_along(bb)) {
    b <- bb[[i]]; g <- gb[[i]]
    expect_identical(b$id, g$id, info = paste(info, "ids", i))
    expect_identical(b$a_index, g$a_index, info = paste(info, "a_index", i))
    expect_identical(b$link_ids, g$link_ids, info = paste(info, "link_ids", i))
    expect_identical(colnames(b$phi), b$link_ids, info = paste(info, "phi columns", i))
    xb_expect_match(b$phi, g$phi, tol = tol, info = paste(info, "phi", i))
    for (nm in c("rate", "beta_score", "treated", "control")) {
      xb_expect_match(b$components[[nm]], g$components[[nm]], tol = tol, info = paste(info, nm, i))
    }
    # The kernel sums exactly ((treated + control) + beta_score) + rate, the
    # general route's order; R's elementwise sums reproduce it bitwise.
    expect_identical(b$phi, ((b$components$treated + b$components$control) + b$components$beta_score) +
      b$components$rate, info = paste(info, "four-term order", i))
  }
}
xb_compare_cov <- function(cp, cg, tol = 1e-12, info = NULL) {
  for (nm in c("estimates", "mean_phi", "covariance")) xb_expect_match(cp[[nm]], cg[[nm]], tol = tol, info = paste(info, nm))
  expect_identical(cp$n_a_units, cg$n_a_units)
  expect_equal(sqrt(diag(cp$covariance)), sqrt(diag(cg$covariance)), tolerance = tol)
}

# ---------------------------------------------------------------------------
# Steps 1-5: nonconstant p = 2 / p = 3 designs through the whole pipeline
# ---------------------------------------------------------------------------

test_that("shared fixed-X p = 2 and p = 3: the batch route is selected and equals the general route", {
  d <- xb_panel()
  for (p in c(2L, 3L)) for (chunk in c(7L, 16L)) {
    ctrl <- xb_ctrl(d, p, chunk = chunk)
    res_b <- xb_run(api, d, ctrl, p)
    res_g <- xb_run(xb_general(), d, ctrl, p)
    info <- paste("p =", p, "chunk =", chunk)
    pb <- res_b$resource_accounting$point_batch
    expect_false(is.null(pb), info = info)
    expect_null(res_g$resource_accounting$point_batch)
    expect_identical(pb$if_plan$p_covariates, p)
    expect_identical(pb$if_plan$kernel, "C_gmm_ipw_if_chunk")
    expect_identical(pb$if_plan$normalizations_performed, 1L)
    expect_false(isTRUE(pb$if_plan$full_n_by_L_retained))
    # The batch fitter reused the one fitted propensity: no refit, no solver call.
    expect_identical(pb$point_fit$propensity_solver_calls, 0L)
    expect_false(isTRUE(pb$point_fit$preflight$propensity_refit_performed))
    expect_identical(pb$point_fit$preflight$p_covariates, as.double(p))
    # Every structural, support, basis, GMM and design field agrees with the
    # general route (doubles at 1e-12, everything else exactly).
    for (f in xb_fields) xb_expect_match(xb_strip_env(res_b[[f]]), xb_strip_env(res_g[[f]]), info = paste(info, f))
    expect_identical(res_b$structural_certificate$retained_index, res_g$structural_certificate$retained_index)
    expect_identical(names(res_b$point_links_used), names(res_g$point_links_used))
    # Point records of the retained links: means, the one fitted propensity's
    # coefficients, the p loadings and rate coefficients, partials, masses, counts.
    for (j in seq_along(res_b$point_links_used)) {
      xb_expect_match(xb_record(res_b$point_links_used[[j]]), xb_record(res_g$point_links_used[[j]]), info = paste(info, "point link", j))
      expect_length(res_b$point_links_used[[j]]$loadings$loading_beta, p)
      expect_identical(res_b$point_links_used[[j]]$coefficients, res_b$nuisance_registry$fitted_views[[1L]]$coefficients)
    }
    # Full-F covariance (`V`), link estimates (`mhat`), mean IF and standard errors.
    cs_b <- xb_strip_env(res_b$covariance_stream); cs_g <- xb_strip_env(res_g$covariance_stream)
    expect_identical(names(cs_b), names(cs_g))
    for (nm in c("mhat", "V", "mean_phi", "n_a_units")) {
      xb_expect_match(cs_b[[nm]], cs_g[[nm]], info = paste(info, "covariance_stream", nm))
    }
    expect_identical(dimnames(cs_b$V), dimnames(cs_g$V))
    expect_equal(sqrt(diag(cs_b$V)), sqrt(diag(cs_g$V)), tolerance = 1e-12)
    # Check actual full-F reporting, not only the A-frame moment covariance.
    expect_true(all(c("reporting", "full_frame_transport") %in% names(res_b)))
    expect_identical(res_b$reporting$available, res_g$reporting$available)
    available <- which(res_b$reporting$available)
    expect_gt(length(available), 0L)
    expect_true(is.matrix(res_b$reporting$if_full) && nrow(res_b$reporting$if_full) == length(unique(d$id)))
    expect_true(all(is.finite(res_b$reporting$if_full[, available, drop = FALSE])))
    for (nm in c("estimates", "if_full", "vcov", "se", "analytical_se")) {
      expect_false(is.null(res_b$reporting[[nm]]), info = paste(info, nm))
      xb_expect_match(res_b$reporting[[nm]], res_g$reporting[[nm]], info = paste(info, "full-F", nm))
    }
  }
})

# ---------------------------------------------------------------------------
# Steps 2-4: plan descriptors, all four IF terms, partial and zero-observer chunks
# ---------------------------------------------------------------------------

test_that("p = 2 plan: four IF components, full IF, IDs and order equal the general replay on partial and zero-observer chunks", {
  # Unit positions 15..21 form the third chunk of 7 and sit in cohort 6, outside
  # the nuisance's cohorts {4, Inf}: a chunk with zero observers. 41 = 5 * 7 + 6
  # leaves a partial last chunk.
  d <- xb_panel(extra = 15:21)
  p <- 2L
  ctrl <- xb_ctrl(d, p, chunk = 7L)
  res_b <- xb_run(api, d, ctrl, p, target_cohorts = 4)
  res_g <- xb_run(xb_general(), d, ctrl, p, target_cohorts = 4)
  expect_false(is.null(res_b$resource_accounting$point_batch))
  expect_null(res_g$resource_accounting$point_batch)
  expect_length(res_b$nuisance_registry$fitted_views, 1L)
  fv <- res_b$nuisance_registry$fitted_views[[1L]]
  links <- res_b$point_links_used
  pc <- ctrl$point_control
  po <- res_b$panel_owner
  L <- length(links)
  m <- as.integer(fv$cohorts$m_cohorts)
  plan <- api$ipw_point_if_batch_plan(links, pc, fv)
  expect_identical(dim(plan$link_par), c(4L + p + m, L))
  expect_identical(dim(plan$link_dates), c(2L, L))
  expect_identical(plan$batch$p_covariates, p)
  expect_identical(plan$batch$m_cohorts, m)
  expect_identical(plan$link_ids, names(links))
  expect_identical(plan$coefficients, fv$coefficients)
  expect_length(plan$coefficients, p)
  # Descriptor layout, row by row, from the retained point records.
  for (ell in seq_len(L)) {
    lk <- links[[ell]]
    expect_identical(plan$link_par[, ell], unname(c(lk$means$treated_mean, lk$means$control_mean, lk$means$mean_a,
      lk$means$mean_b, as.double(lk$loadings$loading_beta), as.double(lk$loadings$rate_coef_s1))))
  }
  # Zero-observer chunk 3 and partial chunk 6, from the owner's cohort index.
  a_chunks <- lapply(seq.int(1L, 41L, by = 7L), function(s) seq.int(s, min(41L, s + 6L)))
  n_obs_chunk <- vapply(a_chunks, function(a) sum(!is.na(api$ipw_owner_view_cohort_index(plan$cohorts,
    po$a_cohort_index[a]))), integer(1))
  expect_identical(n_obs_chunk[3L], 0L)
  expect_true(all(n_obs_chunk[-3L] > 0L))
  expect_identical(length(a_chunks[[6L]]), 6L)

  bb <- xb_collect(api$replay_ipw_point_if_plan, plan = plan, control = pc)
  gb <- xb_collect(api$replay_ipw_point_if_view, point_links = links, control = pc)
  expect_length(bb, 6L)
  xb_compare_batches(bb, gb, info = "p = 2 components")
  expect_true(all(bb[[3L]]$phi == 0))
  for (nm in c("rate", "beta_score", "treated", "control")) expect_true(all(bb[[3L]]$components[[nm]] == 0))
  expect_identical(nrow(bb[[6L]]$phi), 6L)
  # Off-risk rows (cohort 6 positions elsewhere would be zero too): every row of
  # a positive chunk belongs to an observer here, every zero row to a nonobserver.
  for (i in c(1L, 2L, 4L, 5L, 6L)) {
    gi <- api$ipw_owner_view_cohort_index(plan$cohorts, po$a_cohort_index[a_chunks[[i]]])
    expect_true(all(bb[[i]]$phi[is.na(gi), ] == 0))
    expect_identical(bb[[i]]$id, po$a_id[a_chunks[[i]]])
  }
  # IF-only replay (the production path) equals the component replay's phi.
  b0 <- xb_collect(api$replay_ipw_point_if_plan, plan = plan, control = pc, include_components = FALSE)
  for (i in seq_along(b0)) { expect_identical(b0[[i]]$phi, bb[[i]]$phi); expect_null(b0[[i]]$components) }
  # Two-pass covariance through the plan equals the general two-pass stream.
  xb_compare_cov(api$stream_ipw_point_covariance_plan(plan, pc), api$stream_ipw_point_covariance_view(links, pc),
    info = "p = 2 covariance")
})

test_that("p = 3 plan equals the general replay and covariance", {
  d <- xb_panel()
  p <- 3L
  ctrl <- xb_ctrl(d, p, chunk = 16L)
  res_b <- xb_run(api, d, ctrl, p)
  expect_false(is.null(res_b$resource_accounting$point_batch))
  fv <- res_b$nuisance_registry$fitted_views[[1L]]
  links <- res_b$point_links_used
  pc <- ctrl$point_control
  plan <- api$ipw_point_if_batch_plan(links, pc, fv)
  expect_identical(nrow(plan$link_par), 4L + p + as.integer(fv$cohorts$m_cohorts))
  bb <- xb_collect(api$replay_ipw_point_if_plan, plan = plan, control = pc)
  gb <- xb_collect(api$replay_ipw_point_if_view, point_links = links, control = pc)
  expect_length(bb, 3L)
  xb_compare_batches(bb, gb, info = "p = 3 components")
  expect_identical(nrow(bb[[3L]]$phi), 9L)
  xb_compare_cov(api$stream_ipw_point_covariance_plan(plan, pc), api$stream_ipw_point_covariance_view(links, pc),
    info = "p = 3 covariance")
})

# ---------------------------------------------------------------------------
# Step 6: reference failure when a term or the loading order is wrong
# ---------------------------------------------------------------------------

test_that("omitting the beta-score or rate term, or permuting the p loadings, breaks agreement with the general route", {
  d <- xb_panel()
  for (p in c(2L, 3L)) {
    ctrl <- xb_ctrl(d, p, chunk = 7L)
    res_b <- xb_run(api, d, ctrl, p)
    fv <- res_b$nuisance_registry$fitted_views[[1L]]
    links <- res_b$point_links_used
    pc <- ctrl$point_control
    plan <- api$ipw_point_if_batch_plan(links, pc, fv)
    good <- xb_collect(api$replay_ipw_point_if_plan, plan = plan, control = pc)
    gb <- xb_collect(api$replay_ipw_point_if_view, point_links = links, control = pc)
    xb_compare_batches(good, gb, info = paste("p =", p, "intact"))
    beta_rows <- 4L + seq_len(p)
    rate_rows <- 4L + p + seq_len(nrow(plan$link_par) - 4L - p)
    mutants <- list(
      no_beta_score = function(x) { x$link_par[beta_rows, ] <- 0; x },
      no_rate = function(x) { x$link_par[rate_rows, ] <- 0; x },
      permuted_loadings = function(x) { x$link_par[beta_rows, ] <- x$link_par[rev(beta_rows), ]; x })
    for (nm in names(mutants)) {
      mut <- xb_reseal(plan, "point_if_batch_plan", mutants[[nm]])
      mb <- xb_collect(api$replay_ipw_point_if_plan, plan = mut, control = pc)
      phi_m <- do.call(rbind, lapply(mb, `[[`, "phi")); phi_g <- do.call(rbind, lapply(gb, `[[`, "phi"))
      if (nm == "no_rate") {
        # This complete two-cohort fit can have a zero rate loading. Check its
        # arithmetic here; a nonzero descriptor witness below detects omission.
        three <- do.call(rbind, lapply(good, function(b) (b$components$treated + b$components$control) + b$components$beta_score))
        xb_expect_match(phi_m, three, info = paste("p =", p, "original zero-rate arithmetic"))
      } else expect_gt(max(abs(phi_m - phi_g)), 1e-8, label = paste("p =", p, nm))
      # Only the altered term moves; the two outcome terms are untouched.
      for (i in seq_along(mb)) {
        expect_identical(mb[[i]]$components$treated, good[[i]]$components$treated)
        expect_identical(mb[[i]]$components$control, good[[i]]$components$control)
        if (nm == "no_beta_score") expect_true(all(mb[[i]]$components$beta_score == 0))
        if (nm == "no_rate") expect_true(all(mb[[i]]$components$rate == 0))
      }
    }
    # The general reference keeps all four terms: dropping either one from its
    # own components no longer reproduces its IF.
    for (i in seq_along(gb)) {
      three <- (gb[[i]]$components$treated + gb[[i]]$components$control) + gb[[i]]$components$rate
      expect_gt(max(abs(three - gb[[i]]$phi)), 1e-8)
    }
    # Explicit synthetic descriptor witness, not a newly fitted estimator:
    # assign a known nonzero observer rate coefficient to the SAME linked
    # point records, and independently replay those descriptors generically.
    rate_links <- lapply(links, function(lk) xb_reseal(lk, "point_link_view", function(x) {
      x$loadings$rate_coef_s1[] <- 0.25; x }))
    names(rate_links) <- names(links)
    rate_plan <- api$ipw_point_if_batch_plan(rate_links, pc, fv)
    rate_good <- xb_collect(api$replay_ipw_point_if_plan, plan = rate_plan, control = pc)
    rate_ref <- xb_collect(api$replay_ipw_point_if_view, point_links = rate_links, control = pc)
    xb_compare_batches(rate_good, rate_ref, info = paste("p =", p, "synthetic nonzero rate descriptor"))
    rate_mut <- xb_reseal(rate_plan, "point_if_batch_plan", function(x) { x$link_par[rate_rows, ] <- 0; x })
    rate_bad <- xb_collect(api$replay_ipw_point_if_plan, plan = rate_mut, control = pc)
    expect_gt(max(abs(do.call(rbind, lapply(rate_bad, `[[`, "phi")) -
      do.call(rbind, lapply(rate_ref, `[[`, "phi")))), 1e-8)
    for (i in seq_along(rate_ref)) {
      expect_true(all(rate_ref[[i]]$components$rate == 0.25))
      expect_true(all(rate_bad[[i]]$components$rate == 0))
    }
  }
})

# ---------------------------------------------------------------------------
# Step 3: kernel shape validation before any endpoint read; k = 0 chunks
# ---------------------------------------------------------------------------

test_that("malformed design and descriptor shapes refuse with typed bound errors before endpoint reads; k = 0 is valid", {
  d <- xb_panel()
  p <- 2L
  ctrl <- xb_ctrl(d, p, chunk = 7L)
  res_b <- xb_run(api, d, ctrl, p)
  fv <- res_b$nuisance_registry$fitted_views[[1L]]
  links <- res_b$point_links_used
  pc <- ctrl$point_control
  plan <- api$ipw_point_if_batch_plan(links, pc, fv)
  po <- res_b$panel_owner; fo <- res_b$frame_owner; post <- po$postings
  m <- as.integer(fv$cohorts$m_cohorts)
  # The chunk's shared observer arrays, formed as ipw_point_if_plan_chunk forms them.
  a_seq <- 1:7
  gi <- api$ipw_owner_view_cohort_index(plan$cohorts, po$a_cohort_index[a_seq])
  pos_k <- which(!is.na(gi)); gj_k <- gi[pos_k]; k <- length(pos_k)
  expect_gt(k, 0L)
  w_k <- fo$w_a[a_seq][pos_k]; r_k <- plan$rates[gj_k]; D_k <- as.double(gj_k == 1L)
  Z_k <- fo$Z_a[a_seq[pos_k], , drop = FALSE]
  expect_identical(dim(Z_k), c(k, p))
  odds_k <- exp(drop(Z_k %*% plan$coefficients)); prob_k <- odds_k / (1 + odds_k)
  sb_k <- (w_k / r_k) * (D_k - prob_k)
  kern <- api$ipw_point_if_native_kernel()
  call <- function(...) {
    a <- list(po$observations$y, post$date_start, post$post_a_index, post$post_obs_index, 1L, 7L, pos_k, gj_k,
      w_k, r_k, odds_k, sb_k, Z_k, plan$link_dates, plan$link_par, plan$dimnames, TRUE, plan$max_output_entries)
    edits <- list(...); for (j in names(edits)) a[[as.integer(j)]] <- edits[[j]]
    do.call(.Call, c(list(kern), a))
  }
  good <- call()
  expect_identical(good$status, c(0L, 0L, 0L))
  expect_identical(names(good), c("phi", "rate", "beta_score", "treated", "control", "status"))
  gb <- xb_collect(api$replay_ipw_point_if_view, point_links = links, control = pc)
  xb_expect_match(good$phi, gb[[1L]]$phi, info = "direct kernel chunk 1 vs general")
  xb_expect_match(good$beta_score, gb[[1L]]$components$beta_score, info = "direct kernel beta_score vs general")
  # Independent p-term check from the raw arrays: beta_score = score_base * (Z %*% loading_beta).
  for (ell in seq_along(links)) {
    ref <- sb_k * drop(Z_k %*% links[[ell]]$loadings$loading_beta)
    expect_equal(good$beta_score[pos_k, ell], unname(ref), tolerance = 1e-12)
  }
  # A nonfinite endpoint at an observer row is a status, not an abort ...
  y_nan <- po$observations$y
  y_nan[post$post_obs_index[post$date_start[plan$link_dates[2L, 1L]] - 1L + a_seq[pos_k[1L]]]] <- NaN
  expect_identical(call(`1` = y_nan)$status[1L], 1L)
  # ... but every malformed shape aborts first, so that endpoint is never read.
  bad <- list(
    rows_short = list(`13` = Z_k[-1L, , drop = FALSE]),
    dimless_kp = list(`13` = as.vector(Z_k)),
    zero_columns = list(`13` = matrix(numeric(0), nrow = k, ncol = 0L)),
    rows_long = list(`13` = rbind(Z_k, 0)),
    wrong_dim_type = list(`13` = structure(as.vector(Z_k), dim = c(k, p, 1L))),
    par_too_short = list(`15` = plan$link_par[seq_len(4L + p), , drop = FALSE]),
    par_beta_row_dropped = list(`15` = plan$link_par[-5L, , drop = FALSE]),
    par_columns = list(`15` = plan$link_par[, -1L, drop = FALSE]),
    par_vector = list(`15` = as.vector(plan$link_par)),
    extra_column_design = list(`13` = cbind(Z_k, 1)),
    dates_rows = list(`14` = plan$link_dates[1L, , drop = FALSE]),
    cap_below_output = list(`18` = 7 * length(links) - 1))
  for (nm in names(bad)) {
    expect_error(do.call(call, c(list(`1` = y_nan), bad[[nm]])), class = "cdid_gmm_bound_error", info = nm)
  }
  # k = 0: an empty observer set with a 0 x p design (or a dimless empty
  # vector, p = 1 layout with the matching descriptors) is a valid all-zero chunk.
  empty <- call(`7` = integer(0), `8` = integer(0), `9` = numeric(0), `10` = numeric(0), `11` = numeric(0),
    `12` = numeric(0), `13` = matrix(numeric(0), nrow = 0L, ncol = p))
  expect_identical(empty$status, c(0L, 0L, 0L))
  expect_identical(dim(empty$phi), c(7L, length(links)))
  expect_true(all(empty$phi == 0) && all(empty$rate == 0) && all(empty$beta_score == 0))
  expect_identical(colnames(empty$phi), plan$link_ids)
  # Zero rows allow an extreme column count without allocating a large array.
  # Descriptor validation must use subtraction-safe arithmetic before any read.
  expect_error(call(`7` = integer(0), `8` = integer(0), `9` = numeric(0), `10` = numeric(0), `11` = numeric(0),
    `12` = numeric(0), `13` = matrix(numeric(0), nrow = 0L, ncol = .Machine$integer.max)),
    class = "cdid_gmm_bound_error")
  expect_error(call(`7` = integer(0), `8` = integer(0), `9` = numeric(0), `10` = numeric(0), `11` = numeric(0),
    `12` = numeric(0), `13` = matrix(numeric(0), nrow = 0L, ncol = p), `15` = plan$link_par[seq_len(4L + p), , drop = FALSE]),
    class = "cdid_gmm_bound_error")
  # Plan-level malformed loading sizes: a point link whose loading has p - 1
  # entries is refused by the plan before any descriptor is stored.
  short <- xb_reseal(links[[1L]], "point_link_view", function(x) { x$loadings$loading_beta <- x$loadings$loading_beta[-1L]; x })
  links_short <- links; links_short[[1L]] <- short
  expect_error(api$ipw_point_if_batch_plan(links_short, pc, fv), class = "cdid_error")
  wrong_names <- xb_reseal(links[[1L]], "point_link_view", function(x) {
    names(x$loadings$loading_beta) <- rev(names(x$loadings$loading_beta)); x })
  links_wn <- links; links_wn[[1L]] <- wrong_names
  expect_error(api$ipw_point_if_batch_plan(links_wn, pc, fv), class = "cdid_error")
  # Positional, unnamed loadings are legitimate under the existing compact
  # point-link validator; the new plan must not strengthen that contract.
  links_un <- lapply(links, function(lk) xb_reseal(lk, "point_link_view", function(x) {
    names(x$loadings$loading_beta) <- NULL; x }))
  names(links_un) <- names(links)
  plan_un <- api$ipw_point_if_batch_plan(links_un, pc, fv)
  xb_compare_batches(xb_collect(api$replay_ipw_point_if_plan, plan = plan_un, control = pc),
    xb_collect(api$replay_ipw_point_if_view, point_links = links_un, control = pc), info = "unnamed positional loadings")
})

# ---------------------------------------------------------------------------
# Step 5: charges at the adjacent caps refuse before any kernel call or callback
# ---------------------------------------------------------------------------

test_that("adjacent work and matrix caps refuse the p = 2 plan before any chunk, kernel call or callback", {
  d <- xb_panel()
  p <- 2L
  ctrl <- xb_ctrl(d, p, chunk = 7L)
  res_b <- xb_run(api, d, ctrl, p)
  fv <- res_b$nuisance_registry$fitted_views[[1L]]
  links <- res_b$point_links_used
  pc <- ctrl$point_control
  e <- xb_ns(wrap = "ipw_point_if_plan_chunk")
  plan <- e$ipw_point_if_batch_plan(links, pc, fv)
  expect_identical(xb_n(e, "ipw_point_if_plan_chunk"), 0L)
  m <- plan$batch$m_cohorts
  L <- plan$batch$n_links
  B <- as.double(plan$batch$effective_chunk_rows)
  # The descriptor charge is the generalized (2 + 4 + p + m + 16) * L plus IDs/header;
  # the shared design charge is B * p; the traversal bound grows with p.
  expect_gte(plan$batch$plan_retained_entries, (2 + 4 + p + m + 16) * L)
  expect_gte(plan$batch$chunk_shared_entries, 24 * B + B * p)
  expect_gte(plan$batch$traversal_work_ops_bound, 41 * (40 + 4 * (p - 1)) * L)
  # One below the traversal work bound: refused as a typed bound error with no chunk.
  W <- plan$batch$traversal_work_ops_bound
  pc_w <- pc; pc_w$max_work_ops <- W - 1
  pc_w <- api$validate_ipw_point_view_control(pc_w)
  err <- tryCatch(e$ipw_point_if_batch_plan(links, pc_w, fv), error = identity)
  expect_s3_class(err, "cdid_ipw_owner_bound_error")
  if (identical(err$diagnostics$label, "batch IF plan traversal work")) expect_identical(err$diagnostics$requested, W)
  expect_identical(xb_n(e, "ipw_point_if_plan_chunk"), 0L)
  # Component output: the four extra B x L matrices are charged with the
  # retained plan before any chunk; one below refuses, exactly at admits.
  pre <- plan$preflight
  extra <- 4 * pre$b_by_L_entries
  need <- pre$peak_concurrent_entries + 2 * extra + plan$batch$additional_live_entries
  calls <- 0L
  for (cap in c(need - 1, need)) {
    tight <- xb_reseal(plan, "point_if_batch_plan", function(x) {
      x$control$max_matrix_entries <- cap; x$control_in <- x$control; x })
    calls <- 0L
    r <- tryCatch(e$replay_ipw_point_if_plan(tight, tight$control, function(b) calls <<- calls + 1L, include_components = TRUE),
      error = identity)
    if (cap < need) {
      expect_s3_class(r, "cdid_ipw_owner_bound_error")
      expect_identical(calls, 0L)
    } else {
      expect_false(inherits(r, "error"))
      expect_identical(calls, as.integer(pre$n_chunks))
    }
  }
  # The IF-only production path at the same tightened cap is unaffected by the
  # component charge (its preflight is the original one).
  tight <- xb_reseal(plan, "point_if_batch_plan", function(x) {
    x$control$max_matrix_entries <- need - 1; x$control_in <- x$control; x })
  calls <- 0L
  expect_false(inherits(tryCatch(e$replay_ipw_point_if_plan(tight, tight$control, function(b) calls <<- calls + 1L),
    error = identity), "error"))
  expect_identical(calls, as.integer(pre$n_chunks))
})

# ---------------------------------------------------------------------------
# Steps 1 and 6: exact structural negatives stay generic; p = 1 unchanged
# ---------------------------------------------------------------------------

test_that("p = 1 intercept-only batching is unchanged and equals the general route", {
  d <- xb_panel()
  ctrl <- xb_ctrl(d, 1L, chunk = 7L)
  res_b <- xb_run(api, d, ctrl, 1L)
  res_g <- xb_run(xb_general(), d, ctrl, 1L)
  pb <- res_b$resource_accounting$point_batch
  expect_false(is.null(pb))
  expect_identical(pb$if_plan$p_covariates, 1L)
  for (f in xb_fields) xb_expect_match(xb_strip_env(res_b[[f]]), xb_strip_env(res_g[[f]]), info = paste("p = 1", f))
  fv <- res_b$nuisance_registry$fitted_views[[1L]]
  links <- res_b$point_links_used; pc <- ctrl$point_control
  plan <- api$ipw_point_if_batch_plan(links, pc, fv)
  m <- as.integer(fv$cohorts$m_cohorts)
  # The original 2 x L and (5 + m) x L layout, the original 40 * L per-unit term.
  expect_identical(dim(plan$link_par), c(5L + m, length(links)))
  expect_identical(plan$link_par[5L, ], unname(vapply(links, function(lk) unname(lk$loadings$loading_beta[1L]), numeric(1))))
  bb <- xb_collect(api$replay_ipw_point_if_plan, plan = plan, control = pc)
  gb <- xb_collect(api$replay_ipw_point_if_view, point_links = links, control = pc)
  xb_compare_batches(bb, gb, info = "p = 1")
  xb_compare_cov(api$stream_ipw_point_covariance_plan(plan, pc), api$stream_ipw_point_covariance_view(links, pc),
    info = "p = 1 covariance")
})

test_that("differing nuisances, nonconstant weights, incomplete dates and malformed designs keep the general route", {
  d <- xb_panel()
  p <- 2L
  # Unit-varying (time-invariant) weights: structurally ineligible; general route.
  ctrl <- xb_ctrl(d, p, chunk = 7L)
  res_w <- xb_run(api, d, ctrl, p, weightsname = "w_vary")
  expect_null(res_w$resource_accounting$point_batch)
  expect_null(api$ipw_pipeline_batch_eligibility(res_w$nuisance_registry, res_w$structural_certificate$retained_index,
    res_w$frame_owner, ctrl$point_control))
  # Constant weights supplied by name: eligible (the weight value is irrelevant).
  res_c <- xb_run(api, d, ctrl, p, weightsname = "w_const")
  expect_false(is.null(res_c$resource_accounting$point_batch))
  # Two target cohorts: the retained links span two fitted nuisances.
  d2 <- xb_panel(extra = c(5L, 11L, 23L, 29L, 35L))
  ctrl2 <- xb_ctrl(d2, p, chunk = 7L)
  res_2 <- tryCatch(xb_run(api, d2, ctrl2, p), error = identity)
  if (!inherits(res_2, "error")) {
    expect_null(res_2$resource_accounting$point_batch)
    expect_gt(length(res_2$nuisance_registry$fitted_views), 1L)
  } else {
    expect_s3_class(res_2, "cdid_error")
  }
  # An incomplete date (one unit's date-2 row missing): the posting fact fails;
  # legacy clones agree across all fields while production raw-IPW agrees on
  # all scientific outputs at 1e-12.
  d3 <- xb_panel(drop = "3 2")
  ctrl3 <- xb_ctrl(d3, p, chunk = 7L)
  legacy_ns <- xb_ns()
  general_ns <- xb_general()
  res_3b <- tryCatch(xb_run(api, d3, ctrl3, p), error = identity)
  res_3l <- tryCatch(xb_run(legacy_ns, d3, ctrl3, p), error = identity)
  res_3g <- tryCatch(xb_run(general_ns, d3, ctrl3, p), error = identity)
  if (inherits(res_3g, "error")) {
    expect_identical(class(res_3b), class(res_3g))
    expect_identical(conditionMessage(res_3b), conditionMessage(res_3g))
    expect_identical(class(res_3l), class(res_3g))
    expect_identical(conditionMessage(res_3l), conditionMessage(res_3g))
  } else {
    expect_null(res_3l$resource_accounting$point_batch)
    for (f in xb_fields) xb_expect_match(xb_strip_env(res_3l[[f]]), xb_strip_env(res_3g[[f]]), info = paste("incomplete legacy", f))
    expect_identical(res_3l$reporting, res_3g$reporting)
    for (nm in c("mhat", "V", "mean_phi"))
      expect_identical(res_3l$covariance_stream[[nm]], res_3g$covariance_stream[[nm]], info = paste("incomplete legacy", nm))
    expect_null(res_3b$resource_accounting$point_batch)
    expect_s3_class(res_3b$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
    for (f in c("support", "basis", "aggregations", "design_all", "design_used", "reporting")) {
      xb_expect_match(xb_strip_env(res_3b[[f]]), xb_strip_env(res_3g[[f]]), info = paste("incomplete", f))
    }
    xb_expect_match(res_3b$covariance_stream$mhat, res_3g$covariance_stream$mhat, info = "incomplete covariance_stream mhat")
    xb_expect_match(res_3b$covariance_stream$V, res_3g$covariance_stream$V, info = "incomplete covariance_stream V")
    expect_equal(unname(res_3b$covariance_stream$mean_phi), unname(res_3g$covariance_stream$mean_phi), tolerance = 1e-12)
    xb_expect_match(res_3b$gmm_fit$targets, res_3g$gmm_fit$targets, info = "incomplete gmm_fit targets")
    for (nm in c("sandwich_cov", "materialized_target_if", "n_batches", "max_batch_rows_seen")) {
      xb_expect_match(res_3b$gmm_fit$if_projection[[nm]], res_3g$gmm_fit$if_projection[[nm]], info = paste("incomplete if_projection", nm))
    }
    for (nm in c("n_links", "retained_index", "dropped_index", "E", "design_link_index", "link_kind", "classes")) {
      xb_expect_match(xb_strip_env(res_3b$structural_certificate[[nm]]), xb_strip_env(res_3g$structural_certificate[[nm]]),
        info = paste("incomplete structural_certificate", nm))
    }
    expect_null(res_3g$nuisance_registry$raw_ipw)
    # Each outcome-free declaration owns this route's registry link views.
    # Raw records and legacy sealed views have distinct private schemas.
    raw_schema <- c("link_id", "functional_key", "nuisance_id", "nuisance_ref", "frame_ref", "panel_ref",
      "pair", "defining_pair", "cohorts", "nuisance_mask_verified", "support", "control", "preflight", "ref")
    legacy_schema <- c("link_id", "functional_key", "frame_ref", "panel_ref", "nuisance_ref",
      "nuisance_mask_verified", "pair", "cohorts", "support", "loading_beta_candidate", "candidates",
      "control", "preflight", "accounting", "work", "ref")
    for (route in c("raw", "legacy")) {
      rr <- if (route == "raw") res_3b else res_3g
      reg <- rr$nuisance_registry; cert <- rr$structural_certificate
      expect_identical(rr$reduction_summary$reduction_stage, "outcome_free_before_point_fits")
      expect_identical(cert$link_kind, rep("shared_view_declaration", length(reg$link_views)))
      expect_identical(cert$context_ref, rr$frame_owner$ref)
      refs <- lapply(reg$link_views, `[[`, "ref")
      if (route == "legacy") names(refs) <- NULL # legacy core appends in catalog order
      expect_identical(cert$object_refs, refs)
      idx <- if (route == "raw") seq_along(reg$nuisance_views) else unique(reg$link_to_nuisance_idx)
      nrefs <- lapply(reg$nuisance_views[idx], `[[`, "ref")
      if (route == "legacy") names(nrefs) <- NULL
      expect_identical(cert$nuisance_refs, nrefs)
      for (j in seq_along(reg$link_views)) {
        lv <- reg$link_views[[j]]
        expect_identical(names(lv), if (route == "raw") raw_schema else legacy_schema)
        expect_identical(get("link_view", cert$object_refs[[j]], inherits = FALSE), lv)
        expect_identical(lv$nuisance_ref, reg$nuisance_views[[reg$link_to_nuisance_idx[j]]]$ref)
        expect_identical(lv$frame_ref, rr$frame_owner$ref)
        expect_identical(lv$panel_ref, rr$frame_owner$panel_ref)
        expect_true(lv$nuisance_mask_verified)
      }
    }
    # Every remaining certificate field keeps the existing cross-route contract.
    rest <- setdiff(names(res_3b$structural_certificate), c("nuisance_refs", "object_refs", "context_ref", "ref"))
    expect_identical(rest, setdiff(names(res_3g$structural_certificate), c("nuisance_refs", "object_refs", "context_ref", "ref")))
    xb_expect_match(xb_strip_env(res_3b$structural_certificate[rest]),
      xb_strip_env(res_3g$structural_certificate[rest]), info = "incomplete remaining certificate")
    expect_identical(res_3b$covariance_handling$handler_used, res_3g$covariance_handling$handler_used)
    penalty_identity <- c("weighting_requested", "weighting_canonical", "handler_selected", "handler_invoked",
                          "condition_trigger", "ridge_relative", "target_availability")
    expect_identical(res_3b$covariance_handling[penalty_identity], res_3g$covariance_handling[penalty_identity])
    expect_equal(res_3b$covariance_handling$lambda, res_3g$covariance_handling$lambda, tolerance = 1e-12)
  }
  # Structural design facts on an eligible registry: a frame copy whose design
  # shape or column count disagrees with its owner is declined, never admitted
  # by its values.
  res_b <- xb_run(api, d, ctrl, p)
  reg <- res_b$nuisance_registry; kept <- res_b$structural_certificate$retained_index; fo <- res_b$frame_owner
  expect_false(is.null(api$ipw_pipeline_batch_eligibility(reg, kept, fo, ctrl$point_control)))
  fo_cols <- fo; fo_cols$Z_a <- fo$Z_a[, 1L, drop = FALSE]
  expect_null(api$ipw_pipeline_batch_eligibility(reg, kept, fo_cols, ctrl$point_control))
  fo_count <- fo; fo_count$n_covariates <- 2.5
  expect_null(api$ipw_pipeline_batch_eligibility(reg, kept, fo_count, ctrl$point_control))
  fo_int <- fo; storage.mode(fo_int$Z_a) <- "integer"
  expect_null(api$ipw_pipeline_batch_eligibility(reg, kept, fo_int, ctrl$point_control))
  fo_vec <- fo; fo_vec$Z_a <- as.vector(fo$Z_a)
  expect_null(api$ipw_pipeline_batch_eligibility(reg, kept, fo_vec, ctrl$point_control))
  fo_w <- fo; fo_w$w_a[2L] <- fo_w$w_a[2L] + 0.5
  expect_null(api$ipw_pipeline_batch_eligibility(reg, kept, fo_w, ctrl$point_control))
  # Mixed functional keys or an unverified mask on one retained link decline.
  reg_key <- reg; reg_key$link_views[[kept[1L]]]$functional_key <- "other"
  expect_null(api$ipw_pipeline_batch_eligibility(reg_key, kept, fo, ctrl$point_control))
  reg_mask <- reg; reg_mask$link_views[[kept[1L]]]$nuisance_mask_verified <- FALSE
  expect_null(api$ipw_pipeline_batch_eligibility(reg_mask, kept, fo, ctrl$point_control))
})

# Same operation, not a guessed BLAS reduction order: finite cancellation and
# overflowing products retain the generic R primitive's result/refusal class.
test_that("general-p loading products preserve exact R primitive arithmetic and overflow classification", {
  d <- xb_panel(); ctrl <- xb_ctrl(d, 3L); res <- xb_run(api, d, ctrl, 3L)
  plan <- api$ipw_point_if_batch_plan(res$point_links_used, ctrl$point_control, res$nuisance_registry$fitted_views[[1L]])
  po <- res$panel_owner; post <- po$postings; kern <- api$ipw_point_if_native_kernel()
  cases <- list(decimal = list(z = c(.1,.2,.3), l = c(.3,.2,-.1)),
    huge_cancel = list(z = c(1e308,-1e308,1), l = c(1,1,1)),
    intermediate_overflow = list(z = c(1e308,1e308,-1e308), l = c(1,1,1)),
    product_overflow_cancel = list(z = c(1e200,1e200,1), l = c(1e200,-1e200,1)))
  for (nm in names(cases)) {
    z <- matrix(cases[[nm]]$z, nrow = 1L); l <- cases[[nm]]$l
    par <- matrix(c(0,0,1,1,l,0,0), ncol = 1L); before <- unserialize(serialize(list(z,par),NULL))
    ref <- unname(drop(z %*% l))
    r <- .Call(kern, rep(0.0,length(po$observations$y)), post$date_start, post$post_a_index, post$post_obs_index,
      1L,1L,1L,1L,1.0,1.0,1.0,1.0,z,plan$link_dates[,1L,drop=FALSE],par,list(NULL,"probe"),TRUE,1e6)
    expect_identical(colnames(r$beta_score), "probe", info = nm)
    expect_identical(unname(r$beta_score[1L,1L]), ref, info = nm)
    if (is.finite(ref)) expect_identical(r$status,c(0L,0L,0L), info = nm) else {
      expect_identical(r$status,c(3L,1L,1L), info = nm)
      expect_error(api$ipw_point_if_plan_refuse(r$status), class = "cdid_ipw_numeric_error")
    }
    expect_identical(list(z,par), before, info = paste(nm,"borrowed inputs"))
  }
})
