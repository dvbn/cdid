# CP08 IPW construction: unique nuisance/fitted validation, single seal of
# proved complete-support links, and constructor-body reuse through the
# invocation-local bound context.
#
# UNEXECUTED at authoring time (no R in the authoring job). Every expected
# value comes from the unchanged standalone/general path of the same installed
# build (no-context constructors, public validators), never from the private
# construction path itself. Fixtures are private to this file.

if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
if (!requireNamespace("testthat", quietly = TRUE)) stop("testthat required.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Private fixtures (deterministic, no RNG)
# ---------------------------------------------------------------------------

cc_ctrl <- function(chunk = 7L, retained = 1e7, matrix = 1e7, work = 1e10) api$ipw_owner_control(
  max_chunk_rows = chunk, max_retained_entries = retained, max_matrix_entries = matrix, max_work_ops = work,
  max_output_entries = 1e7)

# Long panel; g = Inf codes never treated. Nonadditive outcome.
cc_panel <- function(N = 41L, T = 20L, drop = NULL, na = NULL) {
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- ifelse(d$id %% 3 == 0, 3, Inf)
  d$y <- 0.3 * d$id + 0.2 * d$t + ifelse(d$t >= d$g, 1, 0) + ((7 * d$id + 3 * d$t) %% 5) / 10
  if (!is.null(drop)) d <- d[!(paste(d$id, d$t) %in% drop), ]
  if (!is.null(na)) d$y[paste(d$id, d$t) %in% na] <- NA
  d
}

cc_owners <- function(d, chunk = 7L, weights = NULL) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", cc_ctrl(chunk))
  u <- d[!duplicated(d$id), ]; u <- u[order(u$id), ]
  w <- if (is.null(weights)) NULL else weights(u$id)
  fo <- api$construct_ipw_frame_owner(po, u$id, matrix(1, nrow(u), 1L), weights = w, control = cc_ctrl(chunk))
  list(d = d, po = po, fo = fo)
}

# Units without RNG: permuted ids, cohorts cycling through `cohorts`, one
# covariate and nonconstant weights (dyadic eighths, or non-dyadic tenths plus
# one third).
cc_units <- function(N = 41L, char_ids = FALSE, cohorts = c(4, 6, Inf), dyadic = TRUE) {
  u <- seq_len(N)
  id <- (u * 37L) %% 1009L
  data.frame(id = if (char_ids) sprintf("id-%04d", id) else id, g = cohorts[u %% length(cohorts) + 1L],
    w = if (dyadic) 1 + (u %% 7L) / 8 else 1 + (u %% 7L) / 10 + 1 / 3, x = (u %% 5L) / 4, stringsAsFactors = FALSE)
}
# Long rows in a permuted order; `drop` removes and `na` blanks "unit-position date" rows.
cc_rs_panel <- function(units, T = 6L, drop = NULL, na = NULL, shift = 0) {
  k <- rep(seq_len(nrow(units)), each = T)
  d <- data.frame(id = units$id[k], t = rep(seq_len(T), nrow(units)), g = units$g[k], stringsAsFactors = FALSE)
  d$y <- 0.3 * k + 0.2 * d$t + ifelse(d$t >= d$g, 1, 0) + ((7 * k + 3 * d$t) %% 5) / 10 + shift
  key <- paste(k, d$t)
  d$y[key %in% na] <- NA
  d <- d[!(key %in% drop), ]
  d[order((seq_len(nrow(d)) * 13L) %% nrow(d)), ]
}
cc_rs_owners <- function(units, d, ctrl = cc_ctrl()) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  fo <- api$construct_ipw_frame_owner(po, units$id, cbind("(Intercept)" = 1, x = units$x), weights = units$w,
    control = ctrl)
  list(units = units, d = d, po = po, fo = fo)
}

cc_pairs <- function(T) {
  out <- list()
  for (s in seq_len(T - 1L)) for (t in seq.int(s + 1L, T)) out[[length(out) + 1L]] <- c(s = s, t = t)
  out
}

cc_public_ctrl <- function(d, n_frame, chunk = NULL, tuning = NULL, weighting = "optimal", n_covariates = 1) {
  tuning <- c(tuning, if (!is.null(chunk)) list(chunk_rows = chunk))
  gs <- unique(d$g)
  dims <- api$chained_public_dimensions_from_counts(n_frame = n_frame, n_panel_units = length(unique(d$id)),
    n_rows = nrow(d), n_dates = length(unique(d$t)), n_cohorts = length(gs), n_targets = sum(is.finite(gs)),
    n_covariates = n_covariates, biters = 99, n_clusters = n_frame)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(tuning, gmm_weighting = weighting),
    gmm_weighting = weighting)$control
}
cc_run <- function(ns, d, ctrl, stage = "outcome_free", weighting = "optimal", ...) ns$run_chained_ipw_pipeline(
  data = d, idname = "id", tname = "t", yname = "y", gname = "g", control = ctrl, xformla = ~1,
  control_group = "nevertreated", base_period = "varying", weighting = weighting, run_numerical_transfer = FALSE,
  root_witness = NULL, structural_reduction_stage = stage, ...)
cc_complete_fields <- c("support_status", "support", "reporting", "basis", "aggregations", "covariance_stream",
  "gmm_fit", "covariance_handling", "auxiliary_reporting", "full_frame_transport", "design_all", "design_used",
  "structural_certificate", "controls", "metadata", "panel_owner", "frame_owner", "selection")

# ---------------------------------------------------------------------------
# Instrumentation: a test-local rebinding of the namespace in which named
# internal functions count their actual calls (recursive ones included).
# ---------------------------------------------------------------------------

cc_counted <- function(wrap = character()) {
  e <- new.env(parent = api)
  counts <- new.env(parent = emptyenv())
  # Keep native resolution in the genuine authenticated
  # namespace, as in the existing CP08 selective-counter helpers. Rebinding
  # the resolver without its own namespace cache repeats DLL discovery on
  # every seal/interrupt query. The genuine resolver still checks live native
  # arithmetic on every call; no guard is replaced by a stand-in.
  genuine <- "gmm_require_native_dll"
  for (nm in ls(api, all.names = TRUE)) {
    f <- get(nm, envir = api, inherits = FALSE)
    if (is.function(f) && !is.primitive(f) && identical(environment(f), api) && !nm %in% genuine) {
      environment(f) <- e
      assign(nm, f, envir = e)
    }
  }
  for (nm in genuine) assign(nm, get(nm, envir = api, inherits = FALSE), envir = e)
  stopifnot(all(vapply(genuine, function(nm)
    identical(get(nm, envir = e, inherits = FALSE), get(nm, envir = api, inherits = FALSE)), logical(1L))))
  make_wrap <- function(nm, inner) {
    force(nm); force(inner)
    function(...) {
      assign(nm, get(nm, envir = counts, inherits = FALSE) + 1, envir = counts)
      inner(...)
    }
  }
  for (nm in wrap) {
    if (!exists(nm, envir = e, inherits = FALSE)) stop("No internal function `", nm, "` to count.", call. = FALSE)
    assign(nm, 0, envir = counts)
    assign(nm, make_wrap(nm, get(nm, envir = e, inherits = FALSE)), envir = e)
  }
  # Seals are counted per binding name ("seal:link_view", "seal:nuisance_view", ...).
  # Link views bound to a nuisance view are additionally counted under
  # "seal:link_view:nuisance" (nuisance-free CP04/auxiliary views are not).
  seal_inner <- get("ipw_owner_seal", envir = e, inherits = FALSE)
  assign("ipw_owner_seal", function(obj, binding) {
    bump <- function(key) assign(key, get0(key, envir = counts, inherits = FALSE, ifnotfound = 0) + 1, envir = counts)
    bump(paste0("seal:", binding))
    if (identical(binding, "link_view") && is.list(obj) && !is.null(obj$nuisance_ref)) bump("seal:link_view:nuisance")
    seal_inner(obj, binding)
  }, envir = e)
  # Auxiliary reporting constructs its own nuisance-bound links after the main
  # registry. Keep the total counter above and separately measure the actual
  # main-registry constructor boundary, without guessing auxiliary counts.
  build_inner <- get("build_pipeline_nuisance_and_link_views", envir = e, inherits = FALSE)
  assign("build_pipeline_nuisance_and_link_views", function(...) {
    before <- get0("seal:link_view:nuisance", envir = counts, inherits = FALSE, ifnotfound = 0)
    out <- build_inner(...)
    after <- get0("seal:link_view:nuisance", envir = counts, inherits = FALSE, ifnotfound = 0)
    assign("seal:link_view:nuisance:registry", after - before, envir = counts)
    out
  }, envir = e)
  e$.counts <- counts
  e
}
cc_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0)
cc_snapshot <- function(e) mget(ls(e$.counts, all.names = TRUE), envir = e$.counts)
cc_delta <- function(e, before, nm) cc_n(e, nm) - (if (is.null(before[[nm]])) 0 else before[[nm]])

# Explicitly decline the new reuse paths in a re-bound copy of the same build.
# Only a fresh constructor object declines the fact: its sealed verification
# still uses the accepted posting fact, reproducing the prior two-seal body.
cc_generic <- function() {
  e <- cc_counted()
  for (nm in c("ipw_owner_context_frame_checked", "ipw_owner_context_nuisance_checked",
               "ipw_link_view_nuisance_checked_has", "ipw_owner_immutable_view_snapshot")) {
    assign(nm, function(...) FALSE, envir = e)
  }
  original_fact <- get("ipw_view_mask_compatibility_fact", envir = e, inherits = FALSE)
  assign("ipw_view_mask_compatibility_fact", function(nuisance_view, link_view, ...) {
    out <- original_fact(nuisance_view, link_view, ...)
    if (is.null(link_view$ref)) { out$decided <- FALSE; out$result <- NULL }
    out
  }, envir = e)
  e
}

# Drop only the object's own fresh seal environment; every other reference
# (frame, panel, nuisance) must stay identical between the compared objects.
cc_strip_ref <- function(x) { x$ref <- NULL; x }
# Replace every environment anywhere in a result by NULL (separate runs own
# fresh environments); everything else, including attributes, is kept.
cc_strip_env <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.list(x)) {
    a <- attributes(x)
    y <- lapply(x, cc_strip_env)
    attributes(y) <- a
    return(y)
  }
  x
}
cc_same_error <- function(a, b) {
  expect_true(inherits(a, "error")); expect_true(inherits(b, "error"))
  expect_identical(class(a), class(b))
  expect_identical(conditionMessage(a), conditionMessage(b))
}

# ---------------------------------------------------------------------------
# Step 2: single seal of proved complete-support links
# ---------------------------------------------------------------------------

test_that("complete support: every fresh link gets one final seal, no stream, identical content to the standalone constructor", {
  T <- 20L
  pairs <- cc_pairs(T)
  expect_length(pairs, 190L)
  for (chunk in c(7L, 16L)) {
    o <- cc_owners(cc_panel(41L, T), chunk)
    ctrl <- cc_ctrl(chunk)
    generic <- cc_generic()
    e <- cc_counted(c("verify_ipw_view_mask_compatibility", "ipw_view_mask_compatibility_fact",
      "ipw_owner_query_pair_mask", "validate_ipw_nuisance_view", "validate_ipw_frame_owner",
      "ipw_owner_chunked_pair_support", "ipw_owner_context_link_support"))
    ctx <- e$ipw_owner_construction_context(o$fo)
    e$ipw_owner_context_check_frame(ctx, o$fo)
    expect_identical(ctx$checked_frame_ref, o$fo$ref)
    nv <- e$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl,
      nuisance_id = "op3:mask1", construction_context = ctx)
    nv_ref <- api$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl,
      nuisance_id = "op3:mask1")
    expect_true(identical(cc_strip_ref(nv), cc_strip_ref(nv_ref), num.eq = FALSE))
    before <- cc_snapshot(e)
    links <- vector("list", length(pairs))
    for (j in seq_along(pairs)) {
      s <- pairs[[j]][["s"]]; t <- pairs[[j]][["t"]]
      id <- paste0("L", s, "_", t)
      v <- e$construct_ipw_link_view(o$fo, s, t, nuisance_view = nv, control = ctrl, link_id = id,
        functional_key = "att", construction_context = ctx)
      # Explicit legacy body: no fresh-object fact shortcut, original seal/re-seal.
      w <- generic$construct_ipw_link_view(o$fo, s, t, nuisance_view = nv, control = ctrl, link_id = id,
        functional_key = "att")
      expect_true(identical(cc_strip_ref(v), cc_strip_ref(w), num.eq = FALSE))
      expect_true(isTRUE(v$nuisance_mask_verified))
      expect_identical(v$nuisance_ref, nv$ref)
      expect_identical(v$frame_ref, o$fo$ref)
      expect_identical(names(v), names(w))
      expect_true(api$ipw_owner_seal_ok(v, "link_view"))
      expect_identical(get("link_view", envir = v$ref, inherits = FALSE)$support, v$support)
      expect_invisible(api$validate_ipw_link_view(v, frame_owner = o$fo, nuisance_view = nv))
      expect_identical(api$verify_ipw_view_mask_compatibility(nv, v, ctrl)$compatible, TRUE)
      links[[j]] <- v
    }
    M <- length(pairs)
    expect_identical(cc_delta(e, before, "seal:link_view"), as.numeric(M))
    expect_identical(cc_n(generic, "seal:link_view"), 2 * as.numeric(M))
    expect_identical(cc_delta(e, before, "ipw_view_mask_compatibility_fact"), as.numeric(M))
    expect_identical(cc_delta(e, before, "verify_ipw_view_mask_compatibility"), 0)
    expect_identical(cc_delta(e, before, "ipw_owner_query_pair_mask"), 0)
    expect_identical(cc_delta(e, before, "ipw_owner_chunked_pair_support"), 0)
    expect_identical(cc_delta(e, before, "ipw_owner_context_link_support"), as.numeric(M))
    expect_identical(ctx$support_reuses, as.numeric(M))
    # Step 3: the checked frame and the checked nuisance template are validated
    # once per invocation, not once per link.
    expect_identical(cc_delta(e, before, "validate_ipw_nuisance_view"), 0)
    expect_identical(cc_delta(e, before, "validate_ipw_frame_owner"), 0)
    refs <- lapply(links, `[[`, "ref")
    expect_length(unique(refs), M)
  }
})

test_that("undecided posting fact keeps the original path: seal, streamed verification, verified re-seal", {
  units <- cc_units(41L)
  # Unit position 1 is cohort 6 (outside g = 4 / K = Inf). Dropping its date-3
  # row makes date 3 incomplete while the cohort-filtered masks stay equal.
  for (chunk in c(7L, 16L)) {
    o <- cc_rs_owners(units, cc_rs_panel(units, T = 6L, drop = "1 3"), cc_ctrl(chunk))
    ctrl <- cc_ctrl(chunk)
    e <- cc_counted(c("verify_ipw_view_mask_compatibility", "ipw_view_mask_compatibility_fact",
      "ipw_owner_query_pair_mask"))
    ctx <- e$ipw_owner_construction_context(o$fo)
    e$ipw_owner_context_check_frame(ctx, o$fo)
    nv <- e$construct_ipw_nuisance_view(o$fo, g = 4, K = Inf, s = 1, t = 2, control = ctrl,
      nuisance_id = "op4:mask1", construction_context = ctx)
    before <- cc_snapshot(e)
    v <- e$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl, link_id = "L3_4",
      functional_key = "att", construction_context = ctx)
    w <- api$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl, link_id = "L3_4",
      functional_key = "att")
    expect_true(identical(cc_strip_ref(v), cc_strip_ref(w), num.eq = FALSE))
    expect_true(isTRUE(v$nuisance_mask_verified))
    # Two seals (the original sequence), one streamed verification, streamed reads.
    expect_identical(cc_delta(e, before, "seal:link_view"), 2)
    expect_identical(cc_delta(e, before, "verify_ipw_view_mask_compatibility"), 1)
    expect_gt(cc_delta(e, before, "ipw_owner_query_pair_mask"), 0)
    # Identical pair: decided by the fact even on an incomplete date.
    before <- cc_snapshot(e)
    nv34 <- e$construct_ipw_nuisance_view(o$fo, g = 4, K = Inf, s = 3, t = 4, control = ctrl,
      nuisance_id = "op4:mask34", construction_context = ctx)
    v2 <- e$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv34, control = ctrl, link_id = "L3_4b",
      functional_key = "att", construction_context = ctx)
    w2 <- api$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv34, control = ctrl, link_id = "L3_4b",
      functional_key = "att")
    expect_true(identical(cc_strip_ref(v2), cc_strip_ref(w2), num.eq = FALSE))
    expect_identical(cc_delta(e, before, "seal:link_view"), 1)
    expect_identical(cc_delta(e, before, "verify_ipw_view_mask_compatibility"), 0)
    # A link whose masks actually differ refuses on the original streamed path
    # with the original class, whether or not a context is bound.
    units2 <- cc_units(41L)
    o2 <- cc_rs_owners(units2, cc_rs_panel(units2, T = 6L, drop = "3 3"), cc_ctrl(chunk)) # unit 3 is cohort 4
    ctx2 <- e$ipw_owner_construction_context(o2$fo)
    e$ipw_owner_context_check_frame(ctx2, o2$fo)
    nv2 <- e$construct_ipw_nuisance_view(o2$fo, g = 4, K = Inf, s = 1, t = 2, control = ctrl,
      nuisance_id = "op4:mask1", construction_context = ctx2)
    a <- tryCatch(e$construct_ipw_link_view(o2$fo, 3, 4, nuisance_view = nv2, control = ctrl,
      construction_context = ctx2), error = identity)
    b <- tryCatch(api$construct_ipw_link_view(o2$fo, 3, 4, nuisance_view = nv2, control = ctrl), error = identity)
    cc_same_error(a, b)
    expect_s3_class(a, "cdid_ipw_owner_mask_mismatch_error")
  }
})

test_that("mask-compatibility bounds and typed refusals precede the posting fact and are unchanged", {
  o <- cc_owners(cc_panel(41L, 6L), 7L)
  ctrl <- cc_ctrl(7L)
  nv <- api$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl, nuisance_id = "n1")
  v <- api$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl)
  step <- min(o$fo$n_a_units, 7)
  m <- nv$cohorts$m_cohorts
  n_a <- o$fo$n_a_units
  need_matrix <- 64 * step + 32 * m + 256
  need_work <- 128 * n_a * (1 + ceiling(log2(m + 1))) + 128 * ceiling(n_a / step) * (1 + ceiling(log2(n_a + 1)))
  e <- cc_counted(c("ipw_owner_pair_masks_equal_fact"))
  ctx <- e$ipw_owner_construction_context(o$fo)
  e$ipw_owner_context_check_frame(ctx, o$fo)
  for (ctrl_bad in list(
    api$ipw_owner_control(max_chunk_rows = 7L, max_retained_entries = 1e7, max_matrix_entries = need_matrix - 1,
      max_work_ops = 1e10, max_output_entries = 1e7),
    api$ipw_owner_control(max_chunk_rows = 7L, max_retained_entries = 1e7, max_matrix_entries = 1e7,
      max_work_ops = need_work - 1, max_output_entries = 1e7))) {
    before <- cc_snapshot(e)
    a <- tryCatch(e$verify_ipw_view_mask_compatibility(nv, v, ctrl_bad, construction_context = ctx), error = identity)
    b <- tryCatch(api$verify_ipw_view_mask_compatibility(nv, v, ctrl_bad), error = identity)
    cc_same_error(a, b)
    expect_s3_class(a, "cdid_ipw_owner_bound_error")
    expect_identical(cc_delta(e, before, "ipw_owner_pair_masks_equal_fact"), 0)
  }
  # At the bound the fact decides without any stream.
  ctrl_edge <- api$ipw_owner_control(max_chunk_rows = 7L, max_retained_entries = 1e7, max_matrix_entries = need_matrix,
    max_work_ops = need_work, max_output_entries = 1e7)
  r <- api$verify_ipw_view_mask_compatibility(nv, v, ctrl_edge)
  expect_identical(r$compatible, TRUE)
  expect_identical(r$chunks_compared, as.integer(ceiling(n_a / step)))
  expect_identical(r$endpoints_read, 0L)
  # Cohort and owner refusals keep their classes and order.
  nv6 <- api$construct_ipw_nuisance_view(o$fo, g = 3, K = c(Inf), s = 1, t = 3, control = ctrl, nuisance_id = "n2")
  v_other <- api$construct_ipw_link_view(o$fo, 1, 3, g = 3, K = Inf, control = ctrl)
  expect_identical(api$verify_ipw_view_mask_compatibility(nv6, v_other, ctrl)$compatible, TRUE)
  o_b <- cc_owners(cc_panel(41L, 6L), 7L)
  v_foreign <- api$construct_ipw_link_view(o_b$fo, 3, 4, g = 3, K = Inf, control = ctrl)
  a <- tryCatch(api$verify_ipw_view_mask_compatibility(nv, v_foreign, ctrl), error = identity)
  expect_s3_class(a, "cdid_ipw_owner_mismatch_error")
  # Malformed/unsealed link objects are refused before any guard.
  bad <- cc_strip_ref(v)
  a <- tryCatch(api$verify_ipw_view_mask_compatibility(nv, bad, ctrl), error = identity)
  expect_s3_class(a, "cdid_ipw_owner_integrity_error")
})

# ---------------------------------------------------------------------------
# Step 3: constructor-body reuse through the bound context; forgery/staleness
# ---------------------------------------------------------------------------

test_that("stale, forged, relabelled and foreign objects never reuse a context check", {
  o <- cc_owners(cc_panel(41L, 8L), 7L)
  ctrl <- cc_ctrl(7L)
  e <- cc_counted(c("validate_ipw_nuisance_view", "validate_ipw_frame_owner"))
  ctx <- e$ipw_owner_construction_context(o$fo)
  e$ipw_owner_context_check_frame(ctx, o$fo)
  nv <- e$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl,
    nuisance_id = "op3:mask1", construction_context = ctx)

  # Stale copy of the nuisance view (same `ref`, different object): refused
  # with the original integrity error, identically to the no-context path.
  stale <- nv; stale$support$q[1L] <- stale$support$q[1L]
  a <- tryCatch(e$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = stale, control = ctrl,
    construction_context = ctx), error = identity)
  b <- tryCatch(api$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = stale, control = ctrl), error = identity)
  cc_same_error(a, b); expect_s3_class(a, "cdid_ipw_owner_integrity_error")
  # Forged support values with the same label and reference.
  forged <- nv; forged$support$rates[] <- 1
  a <- tryCatch(e$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = forged, control = ctrl,
    construction_context = ctx), error = identity)
  expect_s3_class(a, "cdid_ipw_owner_integrity_error")
  # Same label, different sealed object (standalone construction): not the
  # template; validated in full and accepted on its own merits.
  twin <- api$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl, nuisance_id = "op3:mask1")
  expect_false(e$ipw_owner_context_nuisance_checked(ctx, twin, o$fo))
  expect_true(e$ipw_owner_context_nuisance_checked(ctx, nv, o$fo))
  before <- cc_snapshot(e)
  v_twin <- e$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = twin, control = ctrl, construction_context = ctx)
  expect_gte(cc_delta(e, before, "validate_ipw_nuisance_view"), 1)
  expect_identical(v_twin$nuisance_ref, twin$ref)
  # Stale frame owner copy: refused as before, context or not.
  fo_stale <- o$fo; fo_stale$n_a_units <- fo_stale$n_a_units
  a <- tryCatch(e$construct_ipw_link_view(fo_stale, 3, 4, nuisance_view = nv, control = ctrl,
    construction_context = ctx), error = identity)
  b <- tryCatch(api$construct_ipw_link_view(fo_stale, 3, 4, nuisance_view = nv, control = ctrl), error = identity)
  cc_same_error(a, b); expect_s3_class(a, "cdid_ipw_owner_integrity_error")
  a <- tryCatch(e$construct_ipw_nuisance_view(fo_stale, g = 3, K = Inf, s = 1, t = 2, control = ctrl,
    construction_context = ctx), error = identity)
  b <- tryCatch(api$construct_ipw_nuisance_view(fo_stale, g = 3, K = Inf, s = 1, t = 2, control = ctrl),
    error = identity)
  cc_same_error(a, b)
  # Foreign frame owner with this context: full validation, original mismatch.
  o2 <- cc_owners(cc_panel(41L, 8L), 7L)
  a <- tryCatch(e$construct_ipw_link_view(o2$fo, 3, 4, nuisance_view = nv, control = ctrl,
    construction_context = ctx), error = identity)
  b <- tryCatch(api$construct_ipw_link_view(o2$fo, 3, 4, nuisance_view = nv, control = ctrl), error = identity)
  cc_same_error(a, b); expect_s3_class(a, "cdid_ipw_owner_mismatch_error")
  expect_false(e$ipw_owner_context_frame_checked(ctx, o2$fo))
  # Context without a completed frame check (shape of the CP04 variant): the
  # constructors validate the frame and nuisance per call, outputs unchanged.
  ctx_plain <- e$ipw_owner_construction_context(o$fo)
  expect_false(e$ipw_owner_context_frame_checked(ctx_plain, o$fo))
  expect_false(e$ipw_owner_context_nuisance_checked(ctx_plain, nv, o$fo))
  before <- cc_snapshot(e)
  v <- e$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl, construction_context = ctx_plain)
  expect_gte(cc_delta(e, before, "validate_ipw_frame_owner"), 1)
  expect_gte(cc_delta(e, before, "validate_ipw_nuisance_view"), 1)
  w <- api$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl)
  expect_true(identical(cc_strip_ref(v), cc_strip_ref(w), num.eq = FALSE))
  # Repeated nuisance ID in one invocation is ambiguous: no template reuse.
  dup <- e$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 3, control = ctrl,
    nuisance_id = "op3:mask1", construction_context = ctx)
  expect_false(e$ipw_owner_context_nuisance_checked(ctx, nv, o$fo))
  expect_false(e$ipw_owner_context_nuisance_checked(ctx, dup, o$fo))
  before <- cc_snapshot(e)
  v <- e$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl, construction_context = ctx)
  expect_gte(cc_delta(e, before, "validate_ipw_nuisance_view"), 1)
  # The context never leaks a skip into the public validators.
  expect_error(api$validate_ipw_link_view(cc_strip_ref(v)), class = "cdid_ipw_owner_integrity_error")
  expect_error(api$validate_ipw_nuisance_view(stale), class = "cdid_ipw_owner_integrity_error")
  expect_error(api$validate_ipw_frame_owner(fo_stale), class = "cdid_ipw_owner_integrity_error")
})

test_that("cohort/control/operator/date changes and malformed endpoints stay on their original paths", {
  o <- cc_owners(cc_panel(41L, 8L), 7L)
  ctrl <- cc_ctrl(7L)
  e <- cc_counted(character())
  ctx <- e$ipw_owner_construction_context(o$fo)
  e$ipw_owner_context_check_frame(ctx, o$fo)
  nv <- e$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl,
    nuisance_id = "op3:mask1", construction_context = ctx)
  cases <- list(
    list(s = 3, t = 4, g = 3, K = c(Inf, 3)),       # K contains g
    list(s = 3, t = 4, g = Inf, K = 3),             # different operator
    list(s = 4, t = 3, g = 3, K = Inf),             # reversed endpoints
    list(s = 3, t = 99, g = 3, K = Inf),            # missing date
    list(s = 3, t = 4, g = 3, K = c(Inf, 7)))       # unsupported cohort
  for (x in cases) {
    a <- tryCatch(e$construct_ipw_link_view(o$fo, x$s, x$t, g = x$g, K = x$K, nuisance_view = nv, control = ctrl,
      construction_context = ctx), error = identity)
    b <- tryCatch(api$construct_ipw_link_view(o$fo, x$s, x$t, g = x$g, K = x$K, nuisance_view = nv, control = ctrl),
      error = identity)
    if (inherits(b, "error")) cc_same_error(a, b) else {
      expect_true(identical(cc_strip_ref(a), cc_strip_ref(b), num.eq = FALSE))
    }
  }
  # Nuisance-free links (g, K supplied) with and without the context.
  for (x in list(list(s = 3, t = 4, g = 3, K = Inf), list(s = 2, t = 7, g = 3, K = Inf))) {
    a <- tryCatch(e$construct_ipw_link_view(o$fo, x$s, x$t, g = x$g, K = x$K, control = ctrl,
      construction_context = ctx), error = identity)
    b <- tryCatch(api$construct_ipw_link_view(o$fo, x$s, x$t, g = x$g, K = x$K, control = ctrl), error = identity)
    if (inherits(b, "error")) cc_same_error(a, b) else {
      expect_true(identical(cc_strip_ref(a), cc_strip_ref(b), num.eq = FALSE))
      expect_null(a$nuisance_ref)
      expect_true(isTRUE(a$nuisance_mask_verified))
    }
  }
  # A false verification flag or a dropped reference is refused by the validator.
  v <- api$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl)
  unverified <- cc_strip_ref(v); unverified$nuisance_mask_verified <- FALSE
  unverified <- api$ipw_owner_seal(unverified, "link_view")
  expect_error(api$validate_ipw_link_view(unverified, frame_owner = o$fo, nuisance_view = nv),
    class = "cdid_ipw_owner_mask_unverified_error")
  nv_other <- api$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl, nuisance_id = "x")
  expect_error(api$validate_ipw_link_view(v, frame_owner = o$fo, nuisance_view = nv_other),
    class = "cdid_ipw_owner_mismatch_error")
})

test_that("mutable context markers and self-populated checked sets never replace the original owner predicates", {
  o <- cc_owners(cc_panel(41L, 6L), 7L)
  ctrl <- cc_ctrl(7L)
  e <- cc_counted()
  nv <- api$construct_ipw_nuisance_view(o$fo, g = 3, K = Inf, s = 1, t = 2, control = ctrl)

  # Own seals alone do not establish the validity of a referenced panel.
  # This context was never checked, but a caller has populated its marker.
  bad_frame <- cc_strip_ref(o$fo)
  bad_frame$panel_ref <- new.env(parent = emptyenv())
  bad_frame <- api$ipw_owner_seal(bad_frame, "frame_owner")
  spoof <- e$ipw_owner_construction_context(bad_frame)
  spoof$checked_frame_ref <- bad_frame$ref
  expect_false(e$ipw_owner_context_frame_checked(spoof, bad_frame))
  a <- tryCatch(e$construct_ipw_nuisance_view(bad_frame, g = 3, K = Inf, s = 1, t = 2,
    control = ctrl, construction_context = spoof), error = identity)
  b <- tryCatch(api$construct_ipw_nuisance_view(bad_frame, g = 3, K = Inf, s = 1, t = 2,
    control = ctrl), error = identity)
  cc_same_error(a, b)
  expect_s3_class(a, "cdid_ipw_owner_integrity_error")

  # The set contains the exact aliased sealed object, but its panel binding
  # is invalid. Membership cannot replace that original nuisance predicate.
  bad_nv <- cc_strip_ref(nv)
  bad_nv$panel_ref <- new.env(parent = emptyenv())
  bad_nv <- api$ipw_owner_seal(bad_nv, "nuisance_view")
  fake_set <- new.env(parent = emptyenv())
  fake_set$frame_ref <- o$fo$ref
  fake_set$refs <- list(bad_nv$ref)
  fake_set$objs <- list(bad_nv)
  expect_false(e$ipw_link_view_nuisance_checked_has(fake_set, bad_nv$ref, bad_nv, o$fo$ref))
  lv <- api$construct_ipw_link_view(o$fo, 3, 4, nuisance_view = nv, control = ctrl)
  bad_link <- cc_strip_ref(lv); bad_link$nuisance_ref <- bad_nv$ref
  bad_link <- api$ipw_owner_seal(bad_link, "link_view")
  a <- tryCatch(e$ipw_link_view_validate_body(bad_link, frame_owner = o$fo, nuisance_view = bad_nv,
    nuisance_checked = fake_set), error = identity)
  b <- tryCatch(api$validate_ipw_link_view(bad_link, frame_owner = o$fo, nuisance_view = bad_nv), error = identity)
  cc_same_error(a, b)
  expect_s3_class(a, "cdid_ipw_owner_mismatch_error")

  # A locked foreign container can hold a valid view without being its own
  # immutable snapshot. Reference-group memoization must decline that case.
  holder <- new.env(parent = emptyenv())
  assign("nuisance_view", nv, envir = holder)
  lockEnvironment(holder, bindings = TRUE)
  expect_false(e$ipw_owner_immutable_view_snapshot(holder, nv, "nuisance_view"))
  expect_true(e$ipw_owner_immutable_view_snapshot(nv$ref, nv, "nuisance_view"))
})

test_that("non-dyadic weights: context links equal the original chunked scan bitwise, including the typed rate refusal", {
  classes <- list(list(g = 4, K = Inf), list(g = 4, K = c(6, Inf)), list(g = 4, K = c(Inf, 6)),
    list(g = 6, K = c(Inf, 4)))
  for (char_ids in c(FALSE, TRUE)) for (chunk in c(7L, 16L)) for (rule in c(TRUE, FALSE)) {
    units <- cc_units(41L, char_ids, dyadic = FALSE)
    o <- cc_rs_owners(units, cc_rs_panel(units, T = 6L), cc_ctrl(chunk))
    ctrl <- cc_ctrl(chunk)
    e <- cc_counted(c("ipw_owner_chunked_pair_support"))
    ctx <- e$ipw_owner_construction_context(o$fo)
    e$ipw_owner_context_check_frame(ctx, o$fo)
    for (cls in classes) {
      nv <- tryCatch(e$construct_ipw_nuisance_view(o$fo, g = cls$g, K = cls$K, s = 1, t = 2, control = ctrl,
        nuisance_id = paste0("n", cls$g, paste(cls$K, collapse = "_")), require_support = rule,
        construction_context = ctx), error = identity)
      nv_ref <- tryCatch(api$construct_ipw_nuisance_view(o$fo, g = cls$g, K = cls$K, s = 1, t = 2, control = ctrl,
        nuisance_id = paste0("n", cls$g, paste(cls$K, collapse = "_")), require_support = rule), error = identity)
      if (inherits(nv_ref, "error")) { cc_same_error(nv, nv_ref); next }
      expect_true(identical(cc_strip_ref(nv), cc_strip_ref(nv_ref), num.eq = FALSE))
      for (pr in cc_pairs(6L)) {
        v <- tryCatch(e$construct_ipw_link_view(o$fo, pr[["s"]], pr[["t"]], nuisance_view = nv, control = ctrl,
          require_support = rule, construction_context = ctx), error = identity)
        w <- tryCatch(api$construct_ipw_link_view(o$fo, pr[["s"]], pr[["t"]], nuisance_view = nv, control = ctrl,
          require_support = rule), error = identity)
        if (inherits(w, "error")) {
          expect_s3_class(w, "cdid_ipw_rate_error")
          cc_same_error(v, w)
        } else {
          expect_true(identical(cc_strip_ref(v), cc_strip_ref(w), num.eq = FALSE))
          expect_identical(v$support$by_cohort_pair_mass, w$support$by_cohort_pair_mass)
          expect_identical(v$support$rates, w$support$rates)
        }
      }
    }
  }
})

# ---------------------------------------------------------------------------
# Step 1 and the whole pipeline: unique validation, 190 bindings, full output
# ---------------------------------------------------------------------------

test_that("pipeline: generic versus instrumented runs agree on every field; one seal per link; unique validation", {
  d <- cc_panel(41L, 20L)
  n_frame <- length(unique(d$id))
  for (chunk in c(7L, 16L)) {
    ctrl <- cc_public_ctrl(d, n_frame, chunk = chunk)
    e <- cc_counted(c("validate_ipw_fitted_nuisance_view", "validate_ipw_nuisance_view",
      "ipw_link_view_validate_body", "validate_ipw_link_view", "verify_ipw_view_mask_compatibility",
      "certify_ipw_structural_classes_from_views", "ipw_view_structural_classes_core",
      "preflight_ipw_early_structural_reduction"))
    assign("ipw_pipeline_compact_basis_facts", function(...) NULL, envir = e)
    # Snapshot the counters around the structural certification.
    trace <- new.env(parent = emptyenv())
    certify_inner <- get("certify_ipw_structural_classes_from_views", envir = e, inherits = FALSE)
    assign("certify_ipw_structural_classes_from_views", function(registry, ...) {
      trace$before <- cc_snapshot(e)
      trace$K <- length(registry$fitted_views)
      trace$M <- length(registry$link_views)
      out <- certify_inner(registry, ...)
      trace$after <- cc_snapshot(e)
      out
    }, envir = e)
    ref <- cc_run(local({ z <- cc_generic(); z$ipw_pipeline_compact_basis_facts <- function(...) NULL; z }), d, ctrl, stage = "outcome_free")
    res <- cc_run(e, d, ctrl, stage = "outcome_free")
    expect_true(all(cc_complete_fields %in% names(ref)))
    expect_identical(names(res), names(ref))
    for (f in cc_complete_fields) {
      expect_true(identical(cc_strip_env(res[[f]]), cc_strip_env(ref[[f]]), num.eq = FALSE), info = f)
    }
    cert <- res$structural_certificate
    M <- cert$n_links
    expect_identical(M, 190L)
    expect_length(cert$object_refs, 190L)
    expect_true(all(vapply(cert$object_refs, is.environment, logical(1L))))
    expect_length(unique(cert$object_refs), 190L)
    for (r in cert$object_refs) {
      lv <- get("link_view", envir = r, inherits = FALSE)
      expect_true(api$ipw_owner_seal_ok(lv, "link_view"))
      expect_true(isTRUE(lv$nuisance_mask_verified))
      expect_identical(lv$frame_ref, cert$context_ref)
    }
    expect_identical(length(cert$retained_index) + length(cert$dropped_index), 190L)
    expect_identical(sort(c(cert$retained_index, cert$dropped_index)), seq_len(190L))
    expect_identical(cert$E, ref$structural_certificate$E)
    expect_identical(cert$classes, ref$structural_certificate$classes)
    expect_identical(cert$retained_index, ref$structural_certificate$retained_index)
    # One final seal per eligible fresh link on a complete panel.
    expect_identical(cc_n(e, "seal:link_view:nuisance:registry"), as.numeric(M))
    # Unique-pair validation inside certification: the recursive per-link
    # nuisance validation is gone, the per-link body still ran for every link
    # in both the wrapper and the core (and once more in the preflight).
    K <- trace$K
    expect_identical(trace$M, 190L)
    d_body <- trace$after[["ipw_link_view_validate_body"]] - trace$before[["ipw_link_view_validate_body"]]
    d_nv <- trace$after[["validate_ipw_nuisance_view"]] - trace$before[["validate_ipw_nuisance_view"]]
    d_fv <- trace$after[["validate_ipw_fitted_nuisance_view"]] - trace$before[["validate_ipw_fitted_nuisance_view"]]
    d_stream <- trace$after[["verify_ipw_view_mask_compatibility"]] - trace$before[["verify_ipw_view_mask_compatibility"]]
    expect_identical(d_body, 3 * 190)
    expect_lte(d_nv, 6 * K)
    expect_lt(d_nv, 190)
    expect_identical(d_fv, 3 * K)
    expect_identical(d_stream, 0)
    expect_gte(cc_n(e, "ipw_view_structural_classes_core"), 1)
  }
})

test_that("pipeline: permuted ids/rows, partial chunks, missing dates and NA outcomes agree with the generic run", {
  units <- cc_units(41L, cohorts = c(4, 6, Inf))
  for (variant in list(list(drop = NULL, na = NULL), list(drop = c("1 3", "9 11"), na = NULL),
    list(drop = NULL, na = c("4 2", "22 9")))) {
    d <- cc_rs_panel(units, T = 12L, drop = variant$drop, na = variant$na)
    n_frame <- length(unique(d$id))
    for (chunk in c(7L, 16L)) {
      ctrl <- cc_public_ctrl(d, n_frame, chunk = chunk)
      e <- cc_counted(c("verify_ipw_view_mask_compatibility"))
      assign("ipw_pipeline_compact_basis_facts", function(...) NULL, envir = e)
      ref <- tryCatch(cc_run(local({ z <- cc_generic(); z$ipw_pipeline_compact_basis_facts <- function(...) NULL; z }), d, ctrl, stage = "outcome_free"), error = identity)
      res <- tryCatch(cc_run(e, d, ctrl, stage = "outcome_free"), error = identity)
      if (inherits(ref, "error")) { cc_same_error(res, ref); next }
      for (f in cc_complete_fields) {
        expect_true(identical(cc_strip_env(res[[f]]), cc_strip_env(ref[[f]]), num.eq = FALSE), info = f)
      }
      M <- res$structural_certificate$n_links
      expect_true(all(vapply(res$structural_certificate$object_refs, function(r)
        api$ipw_owner_seal_ok(get("link_view", envir = r, inherits = FALSE), "link_view"), logical(1L))))
      # Complete panels: one seal per link and no streamed verification; the
      # incomplete variants keep the streamed path (at least one seal per link).
      if (is.null(variant$drop) && is.null(variant$na)) {
        expect_identical(cc_n(e, "seal:link_view:nuisance:registry"), as.numeric(M))
      } else {
        expect_gte(cc_n(e, "seal:link_view:nuisance:registry"), as.numeric(M))
      }
    }
  }
})

test_that("pipeline: identity weighting and post-point-fit stage keep their existing behaviour", {
  d <- cc_panel(41L, 8L)
  n_frame <- length(unique(d$id))
  for (weighting in c("optimal", "identity")) for (stage in c("outcome_free", "post_point_fit")) {
    ctrl <- cc_public_ctrl(d, n_frame, chunk = 7L, weighting = weighting)
    e <- cc_counted(character())
    ref <- tryCatch(cc_run(cc_generic(), d, ctrl, stage = stage, weighting = weighting), error = identity)
    res <- tryCatch(cc_run(e, d, ctrl, stage = stage, weighting = weighting), error = identity)
    if (inherits(ref, "error")) { cc_same_error(res, ref); next }
    for (f in intersect(cc_complete_fields, names(ref))) {
      expect_true(identical(cc_strip_env(res[[f]]), cc_strip_env(ref[[f]]), num.eq = FALSE), info = f)
    }
  }
})
