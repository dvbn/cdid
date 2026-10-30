# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 balanced-panel performance, stage 2: the trusted batch route (shared
# observer state for retained point fits, one sealed plan, registered kernel
# `C_gmm_ipw_if_chunk` for the two covariance and two GMM projection passes).
# Every expected value comes from the unchanged general route of the same
# build (scalar point fits, `query_ipw_view_chunk` replays, the post-fit
# pipeline stage) or from formulas written here over raw panel rows; the batch
# route never supplies its own truth.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

bb_owner_ctrl <- function(chunk) api$ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = 1e7,
  max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
bb_point_ctrl <- function(chunk) api$ipw_point_view_control(max_chunk_rows = chunk, max_observer_units = 1e4,
  max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
bb_fit_ctrl <- function(chunk) api$ipw_fitted_view_control(max_chunk_rows = chunk, max_observer_units = 1e4,
  max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
bb_glm <- stats::glm.control(epsilon = 1e-12, maxit = 50)

# Balanced long panel: g = 3 treated, Inf never treated, 41 units (a partial
# last chunk for chunk sizes 7 and 16), five dates.
bb_panel <- function(N = 41L, T = 5L, seed = 3L) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:2] <- c(3, Inf)
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$y <- 0.2 * d$id + 0.4 * d$t + ifelse(d$t >= d$g, 1.25, 0) + stats::rnorm(nrow(d))
  d
}

bb_fixture <- function(d = bb_panel(), chunk = 7L, weight = 1, pairs = list(c(1, 2), c(2, 3), c(1, 4), c(3, 5), c(4, 5))) {
  oc <- bb_owner_ctrl(chunk)
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", oc)
  u <- d[!duplicated(d$id), ]; u <- u[order(u$id), ]
  fo <- api$construct_ipw_frame_owner(po, u$id, matrix(1, nrow(u), 1L, dimnames = list(NULL, "(Intercept)")),
    weights = rep(weight, nrow(u)), control = oc)
  nv <- api$construct_ipw_nuisance_view(fo, 3, Inf, 1, 2, oc)
  fv <- api$fit_ipw_nuisance_view(nv, bb_glm, bb_fit_ctrl(chunk))
  lvs <- lapply(seq_along(pairs), function(k) api$construct_ipw_link_view(fo, pairs[[k]][1], pairs[[k]][2],
    nuisance_view = nv, control = oc, link_id = paste0("L", k), functional_key = "A"))
  pc <- bb_point_ctrl(chunk)
  scalar <- lapply(lvs, function(lv) api$fit_ipw_point_link_view(fv, lv, pc))
  list(d = d, po = po, fo = fo, nv = nv, fv = fv, lvs = lvs, pc = pc, oc = oc, chunk = chunk, scalar = scalar)
}

bb_snapshot_fields <- c("link_id", "functional_key", "nuisance_id", "estimate", "treated_mean", "control_mean",
  "pair", "defining_pair", "shared_across_dates", "nuisance_mask_verified", "cohorts", "dimensions", "counts",
  "masses", "means", "coefficients", "beta", "factors", "partials", "loadings", "diagnostics",
  "point_summary_disclosure", "glm_control", "control")
# Prove private certificate refs against this route's sealed registry before
# comparing across routes (whose reference environments are necessarily distinct).
# Positional identity is required: labels never establish identity or order.
bb_certificate_ref_map <- function(certificate, pipeline) {
  bad <- function() stop("Certificate nuisance reference mapping differs from the sealed registry/model.", call. = FALSE)
  reg <- pipeline$nuisance_registry
  if (!inherits(certificate, "cdid_ipw_structural_certificate") ||
      !isTRUE(api$ipw_seal_ok(certificate, "certificate")) ||
      !identical(certificate, pipeline$structural_certificate) ||
      !identical(certificate$context_ref, pipeline$frame_owner$ref)) bad()
  order <- unique(reg$link_to_nuisance_idx)
  if (!is.integer(order) || !length(order) || anyNA(order) ||
      any(order < 1L | order > length(reg$nuisance_views)) ||
      length(reg$fitted_views) != length(reg$nuisance_views) ||
      !identical(names(reg$nuisance_views), reg$nuisance_keys) ||
      !identical(names(reg$fitted_views), reg$nuisance_keys)) bad()
  expected <- unname(lapply(reg$nuisance_views[order], `[[`, "ref"))
  # Exact list equality preserves order, identity and every attribute except
  # the one optional private names attribute proved by registry keys below.
  refs <- certificate$nuisance_refs
  if (!identical(unname(refs), expected) ||
      (!is.null(names(refs)) && !identical(names(refs), unname(reg$nuisance_keys[order])))) bad()
  for (k in seq_along(order)) {
    j <- order[k]; nv <- reg$nuisance_views[[j]]; fv <- reg$fitted_views[[j]]
    api$validate_ipw_nuisance_view(nv, frame_owner = pipeline$frame_owner)
    api$validate_ipw_fitted_nuisance_view(fv, nuisance_view = nv, frame_owner = pipeline$frame_owner)
    if (!identical(get("nuisance_view", envir = refs[[k]], inherits = FALSE), nv) ||
        !identical(get("fitted_nuisance_view", envir = fv$ref, inherits = FALSE), fv) ||
        !identical(fv$nuisance_view_ref, refs[[k]])) bad()
  }
  if (length(pipeline$point_links_used) != length(certificate$retained_index)) bad()
  for (h in seq_along(pipeline$point_links_used)) {
    row <- certificate$retained_index[h]; j <- reg$link_to_nuisance_idx[row]
    pl <- pipeline$point_links_used[[h]]; lv <- reg$link_views[[row]]
    if (!identical(pl$nuisance_view_ref, reg$nuisance_views[[j]]$ref) ||
        !identical(pl$fitted_nuisance_ref, reg$fitted_views[[j]]$ref) ||
        !identical(lv$nuisance_ref, reg$nuisance_views[[j]]$ref) ||
        !identical(pl$link_view_ref, lv$ref)) bad()
  }
  unname(order)
}

bb_strip <- function(x, certificate_pipeline = NULL) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (inherits(x, "cdid_ipw_structural_certificate")) {
    # Route-specific link kind and object references (as in the accepted
    # early-reduction tests); every other certificate field is compared.
    if (!is.null(certificate_pipeline)) {
      bb_certificate_ref_map(x, certificate_pipeline)
      # Descriptor assembly preserves nuisance registry names through lapply;
      # post-fit assembly appends refs to an unnamed list. Normalize only this
      # proved attribute on a comparison copy; the original seal stays intact.
      names(x$nuisance_refs) <- NULL
    }
    x <- unclass(x)
    return(lapply(x[setdiff(names(x), c("link_kind", "object_refs", "ref"))],
      bb_strip, certificate_pipeline = certificate_pipeline))
  }
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, bb_strip,
    certificate_pipeline = certificate_pipeline); attributes(x) <- a }
  x
}
bb_snap <- function(pl) bb_strip(unclass(pl)[bb_snapshot_fields])

# Independent four-term IF of one link over A-units `a`, from raw panel rows.
bb_reference_phi <- function(fx, pl, a) {
  ids <- fx$po$a_id[a]
  gcoh <- fx$po$a_g[a]
  gi <- ifelse(gcoh == 3, 1L, ifelse(is.infinite(gcoh), 2L, NA_integer_))
  yat <- function(id, t) vapply(id, function(i) fx$d$y[fx$d$id == i & fx$d$t == t], numeric(1L))
  phi <- numeric(length(a)); comp <- list(rate = phi, beta_score = phi, treated = phi, control = phi)
  k <- which(!is.na(gi))
  if (length(k)) {
    g <- gi[k]; w <- fx$fo$w_a[a[k]]
    r <- pl$masses$rates[g]
    eta <- drop(matrix(1, length(k), 1L) %*% pl$coefficients)
    odds <- exp(eta); prob <- odds / (1 + odds)
    sb <- (w / r) * (as.double(g == 1L) - prob)
    zl <- drop(matrix(1, length(k), 1L) %*% pl$loadings$loading_beta)
    dy <- yat(ids[k], pl$pair$t) - yat(ids[k], pl$pair$s)
    tr <- ifelse(g == 1L, (w / r) / pl$means$mean_a * (dy - pl$means$treated_mean), 0)
    co <- ifelse(g == 1L, 0, -(w * odds / r) / pl$means$mean_b * (dy - pl$means$control_mean))
    rate <- w * pl$loadings$rate_coef_s1[g]
    comp$rate[k] <- rate; comp$beta_score[k] <- sb * zl; comp$treated[k] <- tr; comp$control[k] <- co
    phi[k] <- tr + co + sb * zl + rate
  }
  list(phi = phi, components = comp)
}

bb_kernel_call <- function(plan, a_s, a_e, ..., kernel = api$ipw_point_if_native_kernel()) {
  # The plan chunk's own argument list, with optional replacements.
  a_seq <- seq.int(as.integer(a_s), as.integer(a_e))
  fo <- plan$frame_owner; po <- plan$panel_owner
  gi <- api$ipw_owner_view_cohort_index(plan$cohorts, po$a_cohort_index[a_seq])
  pos <- which(!is.na(gi)); g <- gi[pos]; w <- fo$w_a[a_seq][pos]; r <- plan$rates[g]
  Z <- fo$Z_a[a_seq[pos], , drop = FALSE]; odds <- exp(drop(Z %*% plan$coefficients))
  sb <- (w / r) * (as.double(g == 1L) - odds / (1 + odds))
  args <- list(y = po$observations$y, date_start = po$postings$date_start, post_a = po$postings$post_a_index,
    post_obs = po$postings$post_obs_index, a_start = as.integer(a_s), b_len = length(a_seq), pos = pos, g = g,
    w = w, r = r, odds = odds, sb = sb, z = Z, dates = plan$link_dates, par = plan$link_par,
    dimnames = plan$dimnames, components = FALSE, max = plan$max_output_entries)
  repl <- list(...)
  for (nm in names(repl)) args[nm] <- list(repl[[nm]])
  do.call(.Call, c(list(kernel), unname(args)))
}

fx <- bb_fixture()

test_that("batch point fits equal the scalar fitter field by field and record the shared state truthfully", {
  for (case in list(list(fx = fx), list(fx = bb_fixture(chunk = 16L, weight = 1.5)))) {
    f <- case$fx
    tr <- new.env(parent = emptyenv())
    batch <- api$fit_ipw_point_link_view_batch(f$fv, f$lvs, f$pc, trace_env = tr)
    expect_identical(length(batch), length(f$scalar))
    for (i in seq_along(batch)) {
      expect_identical(bb_snap(batch[[i]]), bb_snap(f$scalar[[i]]))
      expect_silent(api$validate_ipw_point_link_view(batch[[i]], fitted_view = f$fv, link_view = f$lvs[[i]]))
      # Revalidation recomputes through the scalar fitter and compares every snapshot field.
      expect_identical(bb_snap(api$revalidate_ipw_point_link_view(batch[[i]], f$fv, f$lvs[[i]], f$pc)), bb_snap(f$scalar[[i]]))
      expect_true(batch[[i]]$accounting$shared_observer_state)
      expect_false(batch[[i]]$accounting$observer_buffer_released)
      expect_identical(batch[[i]]$work$propensity_solver_calls, 0L)
      expect_identical(batch[[i]]$work$chunks_queried, 0L)
    }
    n_obs <- f$fv$counts$n_observed_risk_units
    expect_true(tr$shared_state_released)
    expect_false(tr$shared_state_active)
    expect_identical(tr$links_fitted, length(batch))
    expect_identical(tr$endpoints_read, 2 * n_obs * length(batch))
    expect_identical(tr$propensity_solver_calls, 0L)
    expect_lte(tr$preflight$peak_concurrent_entries, f$pc$max_matrix_entries)
  }
})

test_that("plan chunks equal the general replay and an independent raw-row formula, components included", {
  pl <- fx$scalar
  plan <- api$ipw_point_if_batch_plan(pl, fx$pc, fx$fv)
  n <- fx$fo$n_a_units
  for (rng in list(c(1L, 7L), c(8L, 14L), c(36L, n), c(n, n), c(5L, 9L))) {
    gen <- api$evaluate_ipw_point_if_view_chunk(pl, rng[1], rng[2], fx$pc, include_components = TRUE, use_nuisance_cache = FALSE)
    got <- api$ipw_point_if_plan_chunk(plan, rng[1], rng[2], api$ipw_point_if_native_kernel(), include_components = TRUE)
    expect_identical(got$phi, gen$phi)
    expect_identical(got$components, gen$components)
    expect_identical(got$id, gen$id)
    expect_identical(got$a_index, gen$a_index)
    for (j in seq_along(pl)) {
      ref <- bb_reference_phi(fx, pl[[j]], seq.int(rng[1], rng[2]))
      expect_identical(unname(got$phi[, j]), ref$phi)
      for (cn in names(ref$components)) expect_identical(unname(got$components[[cn]][, j]), ref$components[[cn]])
    }
  }
  # Whole replays: every callback batch is identical, for full and partial chunks.
  gen <- list(); bat <- list()
  api$replay_ipw_point_if_view(pl, fx$pc, callback_schema = TRUE, callback = function(b) gen[[length(gen) + 1L]] <<- b)
  info <- api$replay_ipw_point_if_plan(plan, fx$pc, callback = function(b) bat[[length(bat) + 1L]] <<- b)
  expect_identical(bat, gen)
  expect_identical(info$n_batches, as.integer(ceiling(n / fx$chunk)))
  expect_identical(info$normalizations_performed, 0L)
  cv_gen <- api$stream_ipw_point_covariance_view(pl, fx$pc)
  cv_bat <- api$stream_ipw_point_covariance_plan(plan, fx$pc)
  for (nm in c("covariance", "mean_phi", "estimates", "n_a_units")) expect_identical(cv_bat[[nm]], cv_gen[[nm]])
  rec <- plan$batch
  expect_identical(rec$normalizations_performed, 1L)
  expect_lte(rec$producer_workspace_entries, rec$general_producer_workspace_entries)
  expect_identical(rec$additional_live_entries, rec$plan_retained_entries)
  expect_false(rec$full_n_by_L_retained)
})

bb_public_ctrl <- function(d, n_frame, chunk = NULL, tuning = NULL, weighting = "optimal", n_covariates = 1) {
  tuning <- c(tuning, if (!is.null(chunk)) list(chunk_rows = chunk))
  gs <- unique(d$g)
  dims <- api$chained_public_dimensions_from_counts(n_frame = n_frame, n_panel_units = length(unique(d$id)),
    n_rows = nrow(d), n_dates = length(unique(d$t)), n_cohorts = length(gs), n_targets = sum(is.finite(gs)),
    n_covariates = n_covariates, biters = 99, n_clusters = n_frame)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(tuning, gmm_weighting = weighting),
    gmm_weighting = weighting)$control
}
bb_run <- function(d, ctrl, stage = "outcome_free", ...) api$run_chained_ipw_pipeline(data = d, idname = "id", tname = "t",
  yname = "y", gname = "g", control = ctrl, xformla = ~1, control_group = "nevertreated", base_period = "varying",
  weighting = "optimal", run_numerical_transfer = FALSE, root_witness = NULL, structural_reduction_stage = stage, ...)
bb_complete_fields <- c("support_status", "support", "reporting", "basis", "aggregations", "covariance_stream",
  "gmm_fit", "covariance_handling", "auxiliary_reporting", "full_frame_transport", "design_all", "design_used",
  "structural_certificate", "controls", "metadata", "panel_owner", "frame_owner", "selection")

test_that("eligible pipelines equal the accepted post-fit route in every output, partial chunks, IDs, clusters and ridge included", {
  d <- bb_panel(N = 41L, T = 6L, seed = 5L)
  set.seed(9L)
  perm <- d; perm$id <- c(sample(500:900, 41L))[perm$id]; perm <- perm[sample(nrow(perm)), ]
  cl <- rep(seq_len(11L), length.out = 41L)
  cases <- list(
    plain = list(d = d, ctrl = bb_public_ctrl(d, 41L, chunk = 7), extra = list()),
    permuted_clustered = list(d = perm, ctrl = bb_public_ctrl(perm, 41L, chunk = 16),
      extra = list(cluster_vec = cl, aggregations = c("simple", "group", "dynamic", "calendar"), cband = TRUE,
        biters = 99L, seed = 7L)),
    ridge = list(d = d, ctrl = bb_public_ctrl(d, 41L, tuning = list(condition_trigger = 1 - 1e-12, ridge_relative = 0.01)),
      extra = list()))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    e <- do.call(bb_run, c(list(cs$d, cs$ctrl, "outcome_free"), cs$extra))
    l <- do.call(bb_run, c(list(cs$d, cs$ctrl, "post_point_fit"), cs$extra))
    expect_false(is.null(e$resource_accounting$point_batch), info = nm)
    expect_null(l$resource_accounting$point_batch, info = nm)
    # This fixture has exactly the same sole nuisance key on both routes.
    # Account for the actual original names before any comparison normalization.
    expect_identical(e$nuisance_registry$nuisance_keys, "op1:mask1", info = nm)
    expect_identical(l$nuisance_registry$nuisance_keys, "op1:mask1", info = nm)
    expect_identical(names(e$structural_certificate$nuisance_refs), "op1:mask1", info = nm)
    expect_null(names(l$structural_certificate$nuisance_refs), info = nm)
    expect_identical(bb_certificate_ref_map(e$structural_certificate, e), 1L, info = nm)
    expect_identical(bb_certificate_ref_map(l$structural_certificate, l), 1L, info = nm)
    # A same-label ref from the other valid route is stale here; even a newly
    # sealed certificate cannot use it. A wrong private label also refuses.
    if (nm == "plain") {
      wrong <- e$structural_certificate
      wrong$nuisance_refs[[1L]] <- l$structural_certificate$nuisance_refs[[1L]]
      wrong <- api$ipw_seal(wrong, "certificate")
      other <- e; other$structural_certificate <- wrong
      expect_error(bb_strip(wrong, other), "Certificate nuisance reference mapping")
      wrong <- e$structural_certificate
      names(wrong$nuisance_refs) <- "op1:mask2"
      wrong <- api$ipw_seal(wrong, "certificate")
      other <- e; other$structural_certificate <- wrong
      expect_error(bb_strip(wrong, other), "Certificate nuisance reference mapping")
    }
    for (f in bb_complete_fields) expect_identical(bb_strip(e[[f]], e), bb_strip(l[[f]], l), info = paste(nm, f))
    expect_true(api$ipw_seal_ok(e$structural_certificate, "certificate"), info = nm)
    expect_true(api$ipw_seal_ok(l$structural_certificate, "certificate"), info = nm)
    expect_identical(lapply(e$point_links_used, bb_snap), lapply(l$point_links_used, bb_snap), info = nm)
    pb <- e$resource_accounting$point_batch
    expect_true(pb$point_fit$shared_state_released)
    expect_identical(pb$point_fit$propensity_solver_calls, 0L)
    expect_identical(pb$if_plan$kernel, "C_gmm_ipw_if_chunk")
    expect_identical(pb$if_plan$effective_chunk_rows, e$resource_accounting$point_producer$effective_chunk_rows)
  }
  expect_identical(e$covariance_handling$handler_used, "tikhonov")
})

test_that("eligibility is decided from actual facts; every other request keeps the general route unchanged", {
  d <- bb_panel(N = 30L, T = 5L, seed = 7L)
  ctrl <- bb_public_ctrl(d, 30L, chunk = 8)
  expect_false(is.null(bb_run(d, ctrl)$resource_accounting$point_batch))
  variants <- list(
    missing_row = d[!(d$id == 4 & d$t == 3), ],
    explicit_na = within(d, y[id == 6 & t == 5] <- NA),
    second_target = within(d, g[id %in% c(5, 9, 13)] <- 4))
  for (nm in names(variants)) {
    v <- variants[[nm]]
    vc <- bb_public_ctrl(v, 30L, chunk = 8)
    e <- bb_run(v, vc)
    l <- bb_run(v, vc, "post_point_fit")
    expect_null(e$resource_accounting$point_batch, info = nm)
    for (f in bb_complete_fields) {
      if (f == "auxiliary_reporting" && !is.null(e$nuisance_registry$raw_ipw)) {
        ae <- e[[f]]; al <- l[[f]]
        expect_identical(ae$descriptors, al$descriptors, info = nm)
        expect_identical(ae$descriptor_status, al$descriptor_status, info = nm)
        supported <- which(ae$descriptor_status$status == "supported")
        expect_gt(length(supported), 0L, label = nm)
        expect_identical(length(ae$point_links), length(ae$descriptors), info = nm)
        expect_identical(length(al$point_links), length(al$descriptors), info = nm)
        # Raw assembly emits fitted catalog rows first, then supported auxiliary
        # descriptors in declaration order (not a fixed second list slot).
        catalog_order <- seq_len(nrow(e$metadata$catalog$rows))
        early_order <- sort(e$structural_certificate$retained_index)
        expect_identical(e$metadata$catalog$rows, l$metadata$catalog$rows, info = nm)
        expect_identical(e$reduction_summary$n_rows_fitted, length(early_order), info = nm)
        expect_identical(l$reduction_summary$n_rows_fitted, length(catalog_order), info = nm)
        for (k in seq_along(supported)) {
          j <- supported[k]
          desc <- ae$descriptors[[j]]
          for (pl in list(ae$point_links[[j]], al$point_links[[j]])) {
            expect_identical(pl$pair$s, desc$s, info = nm)
            expect_identical(pl$pair$t, desc$t, info = nm)
            expect_identical(pl$cohorts$g, as.double(desc$g), info = nm)
            expect_identical(pl$cohorts$K, as.double(desc$K), info = nm)
          }
          expect_identical(ae$point_links[[j]]$cp04_index,
            as.integer(length(early_order) + k), info = nm)
          expect_identical(al$point_links[[j]]$cp04_index,
            as.integer(length(catalog_order) + k), info = nm)
          al$point_links[[j]]$cp04_index <- ae$point_links[[j]]$cp04_index
        }
        owned_diff <- l$nuisance_registry$raw_ipw$accounting$owned_entries -
          e$nuisance_registry$raw_ipw$accounting$owned_entries
        expect_identical(al$preflight$producer_workspace_entries - ae$preflight$producer_workspace_entries,
          owned_diff, info = nm)
        expect_identical(al$preflight$peak_concurrent_entries - ae$preflight$peak_concurrent_entries,
          owned_diff, info = nm)
        al$preflight$producer_workspace_entries <- ae$preflight$producer_workspace_entries
        al$preflight$peak_concurrent_entries <- ae$preflight$peak_concurrent_entries
        expect_identical(bb_strip(ae), bb_strip(al), info = paste(nm, f))
      } else {
        expect_identical(bb_strip(e[[f]]), bb_strip(l[[f]]), info = paste(nm, f))
      }
    }
  }
  base <- list(data = d, idname = "id", tname = "t", yname = "y", gname = "g", control = ctrl, control_group = "nevertreated",
    run_numerical_transfer = FALSE, root_witness = NULL)
  u <- d[!duplicated(d$id), ]; u <- u[order(u$id), ]
  arg_variants <- list(shared_fixed_design = list(frame_id = u$id, frame_Z = cbind(x = 1 + (u$id %% 2))),
    varying_weights = list(frame_id = u$id, frame_Z = matrix(1, nrow(u), 1L), frame_weights = 1 + (u$id %% 3) / 4),
    pt_all = list(pt_assumption = "all", pt_start = 1),
    identity = list(weighting = "identity", control = bb_public_ctrl(d, 30L, chunk = 8, weighting = "identity")))
  for (nm in names(arg_variants)) {
    args <- arg_variants[[nm]]
    e <- do.call(api$run_chained_ipw_pipeline, c(base[setdiff(names(base), names(args))], args))
    if (nm == "shared_fixed_design") expect_false(is.null(e$resource_accounting$point_batch)) else
      expect_null(e$resource_accounting$point_batch)
  }
  # The predicate itself: positive on the balanced registry, NULL on each
  # broken fact, with no refusal raised.
  reg_of <- function(dd, cohorts = Inf, Z = NULL, w = NULL) {
    po <- api$construct_ipw_panel_owner(dd, "id", "t", "y", "g", bb_owner_ctrl(8L))
    uu <- dd[!duplicated(dd$id), ]; uu <- uu[order(uu$id), ]
    fo <- api$construct_ipw_frame_owner(po, uu$id, if (is.null(Z)) matrix(1, nrow(uu), 1L) else Z, weights = w, control = bb_owner_ctrl(8L))
    meta <- api$build_owner_comparison_catalog(po, targets = 3, anticipation = 0, control_group = "nevertreated",
      conventions = "cp08-batch", control = bb_owner_ctrl(8L), catalog_bounds = list(max_candidates = 1e4, max_rows = 1e4,
      max_pool_entries = 1e5, max_scratch = 1e7, max_work = 1e9), frame_owner = fo)
    reg <- api$build_pipeline_nuisance_and_link_views(meta, fo, bb_owner_ctrl(8L), bb_fit_ctrl(8L), bb_point_ctrl(8L),
      bb_glm, fit_point_links = FALSE)
    list(reg = reg, fo = fo, M = length(reg$link_views))
  }
  ok <- reg_of(d)
  facts <- api$ipw_pipeline_batch_eligibility(ok$reg, seq_len(ok$M), ok$fo, bb_point_ctrl(8L))
  expect_false(is.null(facts))
  expect_identical(facts$fitted_view, ok$reg$fitted_views[[1L]])
  fixed_design <- reg_of(d, Z = matrix(2, 30L, 1L))
  expect_false(is.null(api$ipw_pipeline_batch_eligibility(fixed_design$reg, seq_len(fixed_design$M),
    fixed_design$fo, bb_point_ctrl(8L))))
  for (bad in list(reg_of(variants$missing_row), reg_of(d, w = rep(c(1, 2), 15L)))) {
    expect_null(api$ipw_pipeline_batch_eligibility(bad$reg, seq_len(bad$M), bad$fo, bb_point_ctrl(8L)))
  }
  forged <- ok$reg; forged$link_to_nuisance_idx[2L] <- 2L
  expect_null(api$ipw_pipeline_batch_eligibility(forged, seq_len(ok$M), ok$fo, bb_point_ctrl(8L)))
})

test_that("the plan re-verifies seal, controls and every live binding at boundaries and after callbacks", {
  pl <- fx$scalar
  plan <- api$ipw_point_if_batch_plan(pl, fx$pc, fx$fv)
  rebind <- function(env, nm, value) {
    old <- get(nm, envir = env)
    unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
    function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
  }
  calls <- 0L
  count <- function(b) calls <<- calls + 1L
  # Forged or altered plans and changed controls refuse before any chunk.
  bad_plan <- plan; bad_plan$link_par[1L, 1L] <- 0
  expect_error(api$replay_ipw_point_if_plan(bad_plan, fx$pc, count), class = "cdid_ipw_point_view_batch_error")
  forged <- unclass(plan); forged$ref <- NULL
  expect_error(api$replay_ipw_point_if_plan(forged, fx$pc, count), class = "cdid_ipw_point_view_batch_error")
  expect_error(api$replay_ipw_point_if_plan(plan, bb_point_ctrl(8L), count), class = "cdid_ipw_control_mismatch_error")
  expect_identical(calls, 0L)
  # A swapped owner, nuisance, link or point-link binding refuses before any chunk.
  fo2 <- fx$fo; fo2$sum_w_a <- fo2$sum_w_a + 1
  pl2 <- pl[[2L]]; pl2$estimate <- pl2$estimate + 1
  for (sw in list(list(env = fx$fo$ref, nm = "frame_owner", value = fo2),
                  list(env = fx$fv$ref, nm = "fitted_nuisance_view", value = unclass(fx$fv)),
                  list(env = fx$lvs[[3L]]$ref, nm = "link_view", value = fx$lvs[[1L]]),
                  list(env = pl[[2L]]$ref, nm = "point_link_view", value = pl2))) {
    restore <- rebind(sw$env, sw$nm, sw$value)
    expect_error(api$replay_ipw_point_if_plan(plan, fx$pc, count), class = "cdid_ipw_point_view_batch_error", info = sw$nm)
    restore()
  }
  expect_identical(calls, 0L)
  # A callback that rebinds an owner is caught before the next chunk is read.
  restore <- NULL
  kernel_calls <- 0L
  e <- new.env(parent = api)
  for (nm in c("replay_ipw_point_if_plan", "ipw_point_if_plan_chunk")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$ipw_point_if_native_kernel <- function() { kernel_calls <<- kernel_calls + 1L; api$ipw_point_if_native_kernel() }
  expect_error(e$replay_ipw_point_if_plan(plan, fx$pc, function(b) {
    if (is.null(restore)) restore <<- rebind(fx$nv$ref, "nuisance_view", unclass(fx$nv))
  }), class = "cdid_ipw_point_view_batch_error")
  restore()
  expect_identical(kernel_calls, 1L)
  invisible(api$replay_ipw_point_if_plan(plan, fx$pc, count))
  expect_identical(calls, as.integer(ceiling(fx$fo$n_a_units / fx$chunk)))
})

test_that("the kernel refuses malformed shapes and positions, maps data failures to the general route's conditions, never writes inputs", {
  plan <- api$ipw_point_if_batch_plan(fx$scalar, fx$pc, fx$fv)
  base <- bb_kernel_call(plan, 1L, 7L)
  expect_identical(base$status, c(0L, 0L, 0L))
  expect_identical(base$phi, api$ipw_point_if_plan_chunk(plan, 1L, 7L, api$ipw_point_if_native_kernel())$phi)
  pos <- which(!is.na(api$ipw_owner_view_cohort_index(plan$cohorts, fx$po$a_cohort_index[1:7])))
  bound <- "cdid_gmm_bound_error"
  expect_error(bb_kernel_call(plan, 1L, 7L, pos = as.double(pos)), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, pos = rev(pos)), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, pos = replace(pos, 1L, 99L)), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, b_len = 3L), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, g = replace(rep(1L, length(pos)), 1L, 3L)), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, dates = replace(plan$link_dates, 1L, 99L)), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, par = plan$link_par[1:5, , drop = FALSE]), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, z = cbind(1, rep(1, length(pos)))), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, w = 1), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, components = NA), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, dimnames = list(NULL, "x")), class = bound)
  expect_error(bb_kernel_call(plan, 1L, 7L, max = 7 * ncol(plan$link_par) - 1), class = bound)
  # Posting slots that do not hold their A-unit (shifted or truncated
  # postings) stop before any endpoint is read.
  shifted <- fx$po$postings$date_start + 1L
  expect_identical(bb_kernel_call(plan, 1L, 7L, date_start = shifted)$status[1L], 4L)
  expect_identical(bb_kernel_call(plan, 1L, 7L, post_obs = fx$po$postings$post_obs_index[1:3],
    post_a = fx$po$postings$post_a_index[1:3])$status[1L], 4L)
  expect_error(api$ipw_point_if_plan_refuse(c(4L, 1L, 1L)), class = "cdid_ipw_point_view_batch_error")
  # Endpoint, raw-difference and IF failures carry the general route's conditions.
  y <- fx$po$observations$y
  slot <- function(t, a) fx$po$postings$post_obs_index[fx$po$postings$date_start[t] - 1L + a]
  a1 <- pos[1L]
  y_inf <- y; y_inf[slot(2L, a1)] <- Inf
  st1 <- bb_kernel_call(plan, 1L, 7L, y = y_inf)$status
  expect_identical(st1, c(1L, 1L, a1))
  expect_error(api$ipw_point_if_plan_refuse(st1), "Observed risk-set outcome endpoints contain NA, NaN, or Inf.",
    fixed = TRUE, class = "cdid_ipw_owner_endpoint_error")
  y_big <- y; y_big[slot(2L, a1)] <- 1e308; y_big[slot(1L, a1)] <- -1e308
  st2 <- bb_kernel_call(plan, 1L, 7L, y = y_big)$status
  expect_identical(st2, c(2L, 1L, a1))
  expect_error(api$ipw_point_if_plan_refuse(st2), "Pair difference overflow.", fixed = TRUE, class = "cdid_ipw_outcome_error")
  odds_inf <- rep(Inf, length(pos))
  st3 <- bb_kernel_call(plan, 1L, 7L, odds = odds_inf)$status
  expect_identical(st3[1:2], c(3L, 1L))
  expect_error(api$ipw_point_if_plan_refuse(st3), "Nonfinite IF contribution encountered during replay.", fixed = TRUE,
    class = "cdid_ipw_numeric_error")
  # Borrowed inputs are read only.
  expect_identical(fx$po$observations$y, y)
  expect_identical(plan$link_par, api$ipw_point_if_batch_plan(fx$scalar, fx$pc, fx$fv)$link_par)
})

test_that("a raw-difference overflow refuses identically on the batch and scalar point fits", {
  d <- bb_panel(N = 30L, T = 5L, seed = 7L)
  d$y[d$id == 2 & d$t == 3] <- 1e308; d$y[d$id == 2 & d$t == 4] <- -1e308
  oc <- bb_owner_ctrl(8L)
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", oc)
  u <- d[!duplicated(d$id), ]; u <- u[order(u$id), ]
  fo <- api$construct_ipw_frame_owner(po, u$id, matrix(1, nrow(u), 1L), control = oc)
  nv <- api$construct_ipw_nuisance_view(fo, 3, Inf, 1, 2, oc)
  fv <- api$fit_ipw_nuisance_view(nv, bb_glm, bb_fit_ctrl(8L))
  lvs <- lapply(list(c(1, 2), c(3, 4)), function(p) api$construct_ipw_link_view(fo, p[1], p[2], nuisance_view = nv, control = oc))
  es <- tryCatch(api$fit_ipw_point_link_view(fv, lvs[[2L]], bb_point_ctrl(8L)), error = identity)
  eb <- tryCatch(api$fit_ipw_point_link_view_batch(fv, lvs, bb_point_ctrl(8L)), error = identity)
  expect_s3_class(es, "cdid_ipw_outcome_error")
  expect_identical(class(eb), class(es))
  expect_identical(conditionMessage(eb), conditionMessage(es))
})

test_that("preflights refuse before any observer state, endpoint read or kernel call", {
  e <- new.env(parent = api)
  for (nm in c("fit_ipw_point_link_view_batch", "ipw_point_if_batch_plan", "replay_ipw_point_if_plan")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  touched <- 0L
  e$ipw_owner_view_cohort_index <- function(...) { touched <<- touched + 1L; api$ipw_owner_view_cohort_index(...) }
  e$ipw_point_if_native_kernel <- function() { touched <<- touched + 1L; api$ipw_point_if_native_kernel() }
  tr <- new.env(parent = emptyenv())
  low <- fx$pc; low$max_matrix_entries <- 2000
  expect_error(e$fit_ipw_point_link_view_batch(fx$fv, fx$lvs, low, trace_env = tr), class = "cdid_ipw_owner_bound_error")
  expect_null(tr$preflight_completed)
  few <- fx$pc; few$max_observer_units <- 3
  expect_error(e$fit_ipw_point_link_view_batch(fx$fv, fx$lvs, few), class = "cdid_ipw_point_view_bound_error")
  expect_error(e$ipw_point_if_batch_plan(fx$scalar, low, fx$fv), class = "cdid_ipw_owner_bound_error")
  expect_identical(touched, 0L)
  # A link outside the nuisance's mask class refuses during validation, before any read.
  di <- bb_panel()[-3L, ] # unit 1 (g = 3) misses date 3: equal masks, incomplete dates
  oc <- bb_owner_ctrl(7L)
  po <- api$construct_ipw_panel_owner(di, "id", "t", "y", "g", oc)
  u <- di[!duplicated(di$id), ]; u <- u[order(u$id), ]
  fo <- api$construct_ipw_frame_owner(po, u$id, matrix(1, nrow(u), 1L), control = oc)
  nv <- api$construct_ipw_nuisance_view(fo, 3, Inf, 1, 3, oc)
  fv <- api$fit_ipw_nuisance_view(nv, bb_glm, bb_fit_ctrl(7L))
  lvs <- lapply(list(c(1, 3), c(3, 5)), function(p) api$construct_ipw_link_view(fo, p[1], p[2], nuisance_view = nv, control = oc))
  expect_error(e$fit_ipw_point_link_view_batch(fv, lvs, bb_point_ctrl(7L)), class = "cdid_ipw_point_view_batch_error")
  expect_identical(touched, 0L)
  expect_identical(length(lapply(lvs, function(lv) api$fit_ipw_point_link_view(fv, lv, bb_point_ctrl(7L)))), 2L)
})

bb_shlib <- function(code, name, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  src <- file.path(dir, paste0(name, ".c")); writeLines(code, src)
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", shQuote(basename(src))),
    stdout = TRUE, stderr = TRUE))
  so <- file.path(dir, paste0(name, .Platform$dynlib.ext))
  if (!is.null(attr(out, "status")) || !file.exists(so)) stop("R CMD SHLIB failed:\n", paste(out, collapse = "\n"), call. = FALSE)
  so
}
# A build location that admits loading shared objects (`tempdir()` may be
# mounted noexec), probed once as the native-overhead tests do.
bb_exec_base <- local({
  found <- NULL
  function() {
    if (!is.null(found)) return(found)
    for (base in unique(c(tempdir(), getwd()))) {
      d <- tempfile("bb_probe_", tmpdir = base)
      so <- tryCatch(bb_shlib(c("#include <R.h>", "#include <Rinternals.h>",
        "SEXP C_bb_probe(void) { return Rf_ScalarLogical(1); }"), "bb_probe", d), error = function(e) NULL)
      ok <- !is.null(so) && isTRUE(tryCatch({ dyn.load(so); dyn.unload(so); TRUE }, error = function(e) FALSE))
      unlink(d, recursive = TRUE, force = TRUE)
      if (ok) { found <<- base; return(base) }
    }
    stop("No build directory admits loading a compiled shared object.", call. = FALSE)
  }
})
bb_compile <- function(code, name) bb_shlib(code, name, tempfile(paste0("bb_", name, "_"), tmpdir = bb_exec_base()))

test_that("the kernel's own live arithmetic guard refuses changed rounding and subnormal modes", {
  so <- bb_compile(c("#include <R.h>", "#include <Rinternals.h>", "#include <fenv.h>", "#include <stdint.h>",
    "SEXP C_bb_set_rmode(SEXP m) { int v = INTEGER(m)[0]; int rc = -1;",
    "  if (v == 0) rc = fesetround(FE_TONEAREST);",
    "#ifdef FE_UPWARD", "  else if (v == 1) rc = fesetround(FE_UPWARD);", "#endif",
    "  return Rf_ScalarInteger(rc); }",
    "SEXP C_bb_set_fz(SEXP e) { int en = LOGICAL(e)[0];",
    "#if defined(__aarch64__)",
    "  uint64_t fpcr = 0; __asm__ __volatile__(\"mrs %0, fpcr\" : \"=r\"(fpcr));",
    "  if (en) fpcr |= (1ULL << 24); else fpcr &= ~(1ULL << 24);",
    "  __asm__ __volatile__(\"msr fpcr, %0\" : : \"r\"(fpcr)); return Rf_ScalarLogical(1);",
    "#elif defined(__x86_64__) || defined(__SSE2__)",
    "  unsigned int mx = 0; __asm__ __volatile__(\"stmxcsr %0\" : \"=m\"(mx));",
    "  if (en) mx |= 0x8040U; else mx &= ~0x8040U;",
    "  __asm__ __volatile__(\"ldmxcsr %0\" : : \"m\"(mx)); return Rf_ScalarLogical(1);",
    "#else", "  return Rf_ScalarLogical(0);", "#endif", "}"), "bb_fp_helper")
  helper <- dyn.load(so)
  set_rmode <- getNativeSymbolInfo("C_bb_set_rmode", PACKAGE = helper)$address
  set_fz <- getNativeSymbolInfo("C_bb_set_fz", PACKAGE = helper)$address
  on.exit({ try(.Call(set_rmode, 0L), silent = TRUE); try(.Call(set_fz, FALSE), silent = TRUE); try(dyn.unload(so), silent = TRUE)
    unlink(dirname(so), recursive = TRUE, force = TRUE) }, add = TRUE)
  plan <- api$ipw_point_if_batch_plan(fx$scalar, fx$pc, fx$fv)
  kernel <- api$ipw_point_if_native_kernel()
  ok <- bb_kernel_call(plan, 1L, 7L, kernel = kernel)
  rm_guard_re <- paste0(
    "(FE_TONEAREST is required|",
    "Gradual underflow check failed \\(FTZ/DAZ or subnormal flushing detected\\)\\.|",
    "Round-to-nearest ties-to-even \\(RNTE\\) addition/subtraction probe failed\\.)"
  )
  if (identical(.Call(set_rmode, 1L), 0L)) {
    expect_error(bb_kernel_call(plan, 1L, 7L, kernel = kernel), rm_guard_re, class = "cdid_gmm_native_error")
    expect_error(api$replay_ipw_point_if_plan(plan, fx$pc, function(b) NULL), class = "cdid_gmm_native_error")
    .Call(set_rmode, 0L)
  }
  if (isTRUE(.Call(set_fz, TRUE))) {
    expect_error(bb_kernel_call(plan, 1L, 7L, kernel = kernel), "Gradual underflow check failed", class = "cdid_gmm_native_error")
    .Call(set_fz, FALSE)
  }
  expect_identical(bb_kernel_call(plan, 1L, 7L, kernel = kernel), ok)
})

test_that("a backend without the kernel, or with another arity, refuses; unload is refused safely; no R fallback", {
  ref <- testthat::test_path("reference", "cp02_native_original", "src", "cdid.c")
  expect_true(file.exists(ref))
  code <- readLines(ref)
  old_so <- bb_compile(code, "cdid")
  stub <- c(paste0("SEXP C_gmm_ipw_if_chunk(", paste0("SEXP a", 1:17, collapse = ", "), ") { return R_NilValue; }"))
  at <- grep("^static const R_CallMethodDef CallEntries", code)
  arity <- c(code[seq_len(at - 1L)], stub, code[at:length(code)])
  arity <- sub("  {NULL, NULL, 0}", "  {\"C_gmm_ipw_if_chunk\", (DL_FUNC) &C_gmm_ipw_if_chunk, 17},\n  {NULL, NULL, 0}", arity, fixed = TRUE)
  arity_so <- bb_compile(arity, "cdid")
  mock <- function(so) {
    ns <- new.env(hash = TRUE, parent = baseenv())
    ns$.__NAMESPACE__. <- new.env(hash = TRUE, parent = baseenv())
    ns$.__NAMESPACE__.$spec <- c(name = "cdid_bb_mock", version = "0.0.0")
    ns$.__NAMESPACE__.$DLLs <- list(cdid = dyn.load(so))
    for (nm in ls(api, all.names = TRUE)) {
      if (startsWith(nm, ".__")) next
      val <- get(nm, envir = api, inherits = FALSE)
      if (is.function(val) && identical(environment(val), api)) environment(val) <- ns
      else if (nm %in% c("gmm_native_dll_cache", "ipw_point_if_native_cache")) val <- new.env(parent = emptyenv())
      assign(nm, val, envir = ns)
    }
    lockEnvironment(ns, bindings = TRUE)
    ns
  }
  m_old <- mock(old_so)
  on.exit({ try(dyn.unload(old_so), silent = TRUE)
    unlink(c(dirname(old_so), dirname(arity_so)), recursive = TRUE, force = TRUE) }, add = TRUE)
  expect_length(m_old$gmm_require_native_dll(), 9L)
  expect_error(m_old$ipw_point_if_native_kernel(), "C_gmm_ipw_if_chunk", class = "cdid_gmm_native_error")
  expect_error(m_old$replay_ipw_point_if_plan(api$ipw_point_if_batch_plan(fx$scalar, fx$pc, fx$fv), fx$pc, function(b) NULL),
    class = "cdid_gmm_native_error")
  m_arity <- mock(arity_so)
  expect_error(m_arity$ipw_point_if_native_kernel(), "C_gmm_ipw_if_chunk", class = "cdid_gmm_native_error")
  dyn.unload(arity_so)
  expect_error(m_arity$ipw_point_if_native_kernel(), class = "cdid_gmm_native_error")
  expect_null(m_arity$ipw_point_if_native_cache$address)
  # The current backend resolves the kernel and refreshes it with the resolver's table.
  k1 <- api$ipw_point_if_native_kernel()
  expect_identical(api$ipw_point_if_native_kernel(), k1)
  expect_true(api$gmm_sexp_aliases(api$ipw_point_if_native_cache$syms, api$gmm_require_native_dll()))
})

test_that("a running plan replay responds to SIGINT", {
  if (identical(.Platform$OS.type, "windows")) testthat::skip("POSIX SIGINT child probe (forked child) is unavailable on Windows.")
  plan <- api$ipw_point_if_batch_plan(fx$scalar, fx$pc, fx$fv)
  ready <- tempfile("bb_ready_"); result <- tempfile("bb_result_")
  job <- parallel::mcparallel({
    res <- tryCatch({
      writeLines(as.character(Sys.getpid()), ready)
      repeat api$replay_ipw_point_if_plan(plan, fx$pc, function(b) NULL)
      "NO_INTERRUPT"
    }, interrupt = function(e) "CAUGHT_INTERRUPT")
    writeLines(res, result)
    TRUE
  }, mc.set.seed = FALSE)
  on.exit(try(tools::pskill(job$pid, tools::SIGKILL), silent = TRUE), add = TRUE)
  deadline <- Sys.time() + 20
  while (!file.exists(ready) && Sys.time() < deadline) Sys.sleep(0.02)
  expect_true(file.exists(ready))
  Sys.sleep(0.2)
  expect_true(tools::pskill(job$pid, tools::SIGINT))
  deadline <- Sys.time() + 20
  while (!file.exists(result) && Sys.time() < deadline) Sys.sleep(0.02)
  expect_identical(readLines(result, warn = FALSE)[1L], "CAUGHT_INTERRUPT")
  expect_identical(parallel::mccollect(job, wait = TRUE, timeout = 20)[[1L]], TRUE)
})


test_that("optional component matrices are preflighted with the live plan before any callback or chunk allocation", {
  # Repeated pairs with distinct link IDs keep this witness small (41 units),
  # while making the four extra B x L matrices exceed the non-component
  # reservation. No statistical value is used as a resource reference.
  f <- bb_fixture(chunk = 41L, pairs = rep(list(c(1, 2)), 128L))
  plan <- api$ipw_point_if_batch_plan(f$scalar, f$pc, f$fv)
  ctrl <- f$pc
  ctrl$max_matrix_entries <- plan$preflight$peak_concurrent_entries + plan$batch$additional_live_entries
  plan <- api$ipw_point_if_batch_plan(f$scalar, ctrl, f$fv)
  calls <- 0L
  expect_error(api$replay_ipw_point_if_plan(plan, ctrl, function(batch) calls <<- calls + 1L,
    include_components = TRUE), class = "cdid_ipw_owner_bound_error")
  expect_identical(calls, 0L)
  expect_error(api$ipw_point_if_plan_chunk(plan, 1L, 41L, api$ipw_point_if_native_kernel(),
    include_components = TRUE), class = "cdid_ipw_owner_bound_error")
  # Under an ample bound the existing numerical component results still run,
  # and replay's receipt describes all five matrices and their live allowance.
  full <- api$ipw_point_if_batch_plan(fx$scalar, fx$pc, fx$fv)
  info <- api$replay_ipw_point_if_plan(full, fx$pc, function(batch) NULL, include_components = TRUE)
  expect_true(info$preflight$include_components)
  expect_identical(info$preflight$n_output_matrices_per_chunk, 5L)
  expect_identical(info$preflight$chunk_output_matrix_entries, 5 * full$preflight$b_by_L_entries)
  expect_identical(info$preflight$peak_concurrent_entries,
    full$preflight$peak_concurrent_entries + 8 * full$preflight$b_by_L_entries)
})
