# CP08 compact canonical basis: catalog/mapping reuse on the design pass and the
# early structural certificate from canonical catalog descriptors (before any
# link view) on the complete-roster built-in baseline IPW route.
#
# UNEXECUTED at authoring time (no R in the authoring job). Every expected
# value comes from the unchanged accepted routes of the same installed build:
# the general view-registry outcome-free route (forced in a disclosed test-local
# re-bound namespace copy that declines the compact predicate) and the accepted
# post-point-fit reference (`structural_reduction_stage = "post_point_fit"`,
# genuine namespace). Counter instrumentation lives only in test-local clones
# (`cb_counted`, `cb_general`); the genuine public reference gates run on the
# installed namespace. Fixtures are private to this file.

if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
if (!requireNamespace("testthat", quietly = TRUE)) stop("testthat required.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Private fixtures (deterministic, no RNG)
# ---------------------------------------------------------------------------

cb_ctrl <- function(chunk = 7L, retained = 1e7, matrix = 1e7, work = 1e10, output = 1e7) api$ipw_owner_control(
  max_chunk_rows = chunk, max_retained_entries = retained, max_matrix_entries = matrix, max_work_ops = work,
  max_output_entries = output)

# Complete balanced long panel: cohort 3 (id %% 3 == 0) against never treated
# (g = Inf), one unit-level covariate x, nonadditive outcome. `dates` replaces
# the integer calendar by literal stored values; `drop` removes rows; `id_perm`
# and `rows_perm` permute identifiers and row order without RNG.
cb_panel <- function(N = 41L, T = 20L, dates = NULL, drop = NULL, id_perm = FALSE, rows_perm = FALSE) {
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- ifelse(d$id %% 3 == 0, 3, Inf)
  d$y <- 0.3 * d$id + 0.2 * d$t + ifelse(d$t >= d$g, 1, 0) + ((7 * d$id + 3 * d$t) %% 5) / 10
  d$x <- (d$id %% 5) / 4
  d$w <- 1 + (d$id %% 4) / 4
  if (!is.null(drop)) d <- d[!(paste(d$id, d$t) %in% drop), ]
  if (!is.null(dates)) {
    stopifnot(length(dates) == T)
    d$g <- ifelse(is.finite(d$g), dates[3L], Inf)
    d$t <- dates[d$t]
  }
  if (id_perm) d$id <- (d$id * 37L) %% 1009L
  if (rows_perm) d <- d[order((seq_len(nrow(d)) * 13L) %% nrow(d)), ]
  rownames(d) <- NULL
  d
}

cb_public_ctrl <- function(d, n_frame, chunk = NULL, weighting = "optimal", n_covariates = 1) {
  tuning <- if (!is.null(chunk)) list(chunk_rows = chunk) else NULL
  gs <- unique(d$g)
  dims <- api$chained_public_dimensions_from_counts(n_frame = n_frame, n_panel_units = length(unique(d$id)),
    n_rows = nrow(d), n_dates = length(unique(d$t)), n_cohorts = length(gs), n_targets = sum(is.finite(gs)),
    n_covariates = n_covariates, biters = 99, n_clusters = n_frame)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(tuning, gmm_weighting = weighting),
    gmm_weighting = weighting)$control
}

cb_run <- function(ns, d, ctrl, stage = "outcome_free", weighting = "optimal", xformla = ~1,
                   control_group = "nevertreated", ...) ns$run_chained_ipw_pipeline(
  data = d, idname = "id", tname = "t", yname = "y", gname = "g", control = ctrl, xformla = xformla,
  control_group = control_group, base_period = "varying", weighting = weighting, run_numerical_transfer = FALSE,
  root_witness = NULL, structural_reduction_stage = stage, ...)

cb_fields <- c("support_status", "support", "reporting", "basis", "aggregations", "covariance_stream",
  "gmm_fit", "covariance_handling", "auxiliary_reporting", "full_frame_transport", "design_all", "design_used",
  "controls", "metadata", "panel_owner", "frame_owner", "selection")
cb_cert_fields <- c("weighting", "n_links", "retained_index", "dropped_index", "E", "classes", "design_link_index", "bounds")

# ---------------------------------------------------------------------------
# Instrumentation (disclosed test-local clones of the namespace)
# ---------------------------------------------------------------------------

cb_counted <- function(wrap = character()) {
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
  seal_inner <- get("ipw_owner_seal", envir = e, inherits = FALSE)
  assign("ipw_owner_seal", function(obj, binding) {
    key <- paste0("seal:", binding)
    assign(key, get0(key, envir = counts, inherits = FALSE, ifnotfound = 0) + 1, envir = counts)
    seal_inner(obj, binding)
  }, envir = e)
  e$.counts <- counts
  e
}
cb_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0)
cb_snapshot <- function(e) mget(ls(e$.counts, all.names = TRUE), envir = e$.counts)

# The accepted general view-registry outcome-free route of the same build:
# only the compact predicate is declined; everything else is the installed code.
cb_general <- function(wrap = character()) {
  e <- cb_counted(wrap)
  assign("ipw_pipeline_compact_basis_facts", function(...) NULL, envir = e)
  e
}

cb_strip_env <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.list(x)) {
    a <- attributes(x)
    y <- lapply(x, cb_strip_env)
    attributes(y) <- a
    return(y)
  }
  x
}
cb_same <- function(a, b, info) {
  # Only the GMM provenance's duplicate certificate representation differs:
  # descriptors versus views, and optional list names on reference groups.
  # Actual kinds and class/ref dimensions are checked before normalization;
  # every numerical/provenance/class/forest/metric field remains compared.
  if (is.list(a) && is.list(b) && is.list(a$provenance$reduction$class_certificate) &&
      is.list(b$provenance$reduction$class_certificate) &&
      !identical(a$provenance$reduction$class_certificate$link_kind, b$provenance$reduction$class_certificate$link_kind)) {
    ca <- a$provenance$reduction$class_certificate; cb <- b$provenance$reduction$class_certificate
    kinds <- c("canonical_descriptor_declaration", "shared_view_declaration")
    expect_true(length(ca$link_kind) > 0L && all(ca$link_kind %in% kinds), info = info)
    expect_true(length(cb$link_kind) > 0L && all(cb$link_kind %in% kinds), info = info)
    expect_identical(length(ca$link_kind), length(cb$link_kind), info = info)
    expect_identical(length(ca$nuisance_refs), length(cb$nuisance_refs), info = info)
    expect_identical(ca$link_kind, a$provenance$weight_covariance$reduction$class_certificate$link_kind, info = info)
    expect_identical(cb$link_kind, b$provenance$weight_covariance$reduction$class_certificate$link_kind, info = info)
    for (which in c("reduction", "weight_covariance")) {
      if (which == "reduction") {
        a$provenance$reduction$class_certificate$link_kind <- b$provenance$reduction$class_certificate$link_kind <- NULL
        a$provenance$reduction$class_certificate$nuisance_refs <- unname(a$provenance$reduction$class_certificate$nuisance_refs)
        b$provenance$reduction$class_certificate$nuisance_refs <- unname(b$provenance$reduction$class_certificate$nuisance_refs)
      } else {
        a$provenance$weight_covariance$reduction$class_certificate$link_kind <- b$provenance$weight_covariance$reduction$class_certificate$link_kind <- NULL
        a$provenance$weight_covariance$reduction$class_certificate$nuisance_refs <- unname(a$provenance$weight_covariance$reduction$class_certificate$nuisance_refs)
        b$provenance$weight_covariance$reduction$class_certificate$nuisance_refs <- unname(b$provenance$weight_covariance$reduction$class_certificate$nuisance_refs)
      }
    }
  }
  expect_true(identical(cb_strip_env(a), cb_strip_env(b), num.eq = FALSE), info = info)
}
cb_err <- function(expr) tryCatch({ expr; NULL }, error = identity)
cb_expect_class <- function(e, cls, info = NULL) {
  expect_true(inherits(e, "error"), info = info)
  expect_true(inherits(e, cls), info = paste(info, paste(class(e), collapse = "/")))
}

# ---------------------------------------------------------------------------
# 1. Parity: compact canonical basis versus the general early route (bitwise)
#    and the accepted post-point-fit reference (existing tolerances)
# ---------------------------------------------------------------------------

test_that("compact canonical basis equals the general early route and the post-fit reference on every output", {
  d <- cb_panel(41L, 20L)
  n_frame <- length(unique(d$id))
  for (xf in list(~1, ~x)) for (chunk in c(7L, 16L)) {
    info <- paste(deparse(xf), chunk)
    ctrl <- cb_public_ctrl(d, n_frame, chunk = chunk, n_covariates = length(all.vars(xf)) + 1)
    res <- cb_run(api, d, ctrl, xformla = xf)
    ref_g <- cb_run(cb_general(), d, ctrl, xformla = xf)
    ref_p <- cb_run(api, d, ctrl, stage = "post_point_fit", xformla = xf)

    expect_identical(res$reduction_summary$reduction_stage, "outcome_free_compact_canonical_basis", info = info)
    expect_identical(ref_g$reduction_summary$reduction_stage, "outcome_free_before_point_fits", info = info)
    expect_identical(ref_p$reduction_summary$reduction_stage, "after_point_fits", info = info)
    expect_identical(names(res), names(ref_g), info = info)
    for (f in cb_fields) cb_same(res[[f]], ref_g[[f]], paste(info, f))

    # The same certificate: kept indices, E, classes, design map, bounds.
    cert <- res$structural_certificate
    expect_true(api$ipw_seal_ok(cert, "certificate"), info = info)
    for (f in cb_cert_fields) cb_same(cert[[f]], ref_g$structural_certificate[[f]], paste(info, "cert", f))
    for (f in cb_cert_fields) cb_same(cert[[f]], ref_p$structural_certificate[[f]], paste(info, "cert/post", f))
    expect_identical(cert$n_links, 190L, info = info)
    expect_identical(length(cert$retained_index), 19L, info = info)
    expect_identical(cert$link_kind, rep("canonical_descriptor_declaration", 190L), info = info)
    expect_identical(ref_g$structural_certificate$link_kind, rep("shared_view_declaration", 190L), info = info)
    expect_identical(cert$context_ref, res$frame_owner$ref, info = info)
    expect_length(cert$nuisance_refs, 1L)
    # Every original row's certificate object is its sealed canonical descriptor
    # bound to the actual owners, metadata, nuisance view and fitted model.
    expect_length(cert$object_refs, 190L)
    expect_length(unique(cert$object_refs), 190L)
    fv <- res$nuisance_registry$fitted_views[[1L]]
    for (ell in seq_len(190L)) {
      desc <- get("link_descriptor", envir = cert$object_refs[[ell]], inherits = FALSE)
      expect_true(api$ipw_owner_seal_ok(desc, "link_descriptor"))
      expect_identical(desc$catalog_row, ell)
      expect_identical(desc$frame_ref, res$frame_owner$ref)
      expect_identical(desc$metadata_ref, res$metadata$ref)
      expect_identical(desc$nuisance_view_ref, cert$nuisance_refs[[1L]])
      expect_identical(desc$fitted_nuisance_ref, fv$ref)
      expect_identical(desc$link_id, paste0("link:r", ell))
    }
    expect_identical(res$reduction_summary$kept_indices, ref_g$reduction_summary$kept_indices, info = info)
    expect_identical(res$reduction_summary$expanded_m, ref_g$reduction_summary$expanded_m, info = info)
    expect_identical(res$reduction_summary$n_rows_fitted, 19L, info = info)
    expect_identical(res$reduction_summary$n_link_views_constructed, 19L, info = info)
    expect_equal(ref_g$reduction_summary$n_link_views_constructed, 190, info = info)

    # Retained point links are genuine fits bitwise equal to the general route;
    # dropped rows are outcome-free compact records with the certified expansion.
    expect_length(res$point_links_all, 190L)
    expect_identical(names(res$point_links_all), names(ref_g$point_links_all), info = info)
    for (m in cert$retained_index) {
      expect_true(inherits(res$point_links_all[[m]], "cdid_ipw_point_link_view"), info = info)
      cb_same(res$point_links_all[[m]], ref_g$point_links_all[[m]], paste(info, "point", m))
    }
    for (m in cert$dropped_index) {
      rec <- res$point_links_all[[m]]
      ref_rec <- ref_g$point_links_all[[m]]
      expect_true(inherits(rec, "cdid_ipw_reduced_link_record"), info = info)
      expect_identical(rec$status, "structurally_reduced_outcome_free_compact", info = info)
      expect_false(rec$point_fit_performed); expect_false(rec$link_view_constructed); expect_false(rec$outcomes_read)
      expect_null(rec$link_view_ref)
      expect_identical(rec$descriptor_ref, cert$object_refs[[m]])
      expect_identical(rec$expansion, ref_rec$expansion, info = info)
      expect_identical(rec$link_id, ref_rec$link_id, info = info)
      expect_identical(rec$catalog_row, ref_rec$catalog_row, info = info)
    }

    # Points, expanded moments, full original-unit IF, V, SE and the active ridge
    # against the accepted post-point-fit reference at the existing tolerance.
    expect_equal(res$covariance_stream$mhat, ref_p$covariance_stream$mhat, tolerance = 1e-12, info = info)
    expect_equal(res$covariance_stream$V, ref_p$covariance_stream$V, tolerance = 1e-12, info = info)
    expect_equal(res$reduction_summary$expanded_m, ref_p$reduction_summary$expanded_m, tolerance = 1e-12, info = info)
    expect_equal(res$gmm_fit$targets$estimate, ref_p$gmm_fit$targets$estimate, tolerance = 1e-12, info = info)
    expect_equal(res$reporting$estimates, ref_p$reporting$estimates, tolerance = 1e-12, info = info)
    expect_equal(res$reporting$se, ref_p$reporting$se, tolerance = 1e-12, info = info)
    expect_equal(res$reporting$vcov, ref_p$reporting$vcov, tolerance = 1e-12, info = info)
    expect_equal(res$reporting$if_full, ref_p$reporting$if_full, tolerance = 1e-12, info = info)
    expect_equal(res$full_frame_transport$Phi_full, ref_p$full_frame_transport$Phi_full, tolerance = 1e-12, info = info)
    expect_identical(nrow(res$full_frame_transport$Phi_full), n_frame, info = info)
    expect_identical(cb_strip_env(res$covariance_handling), cb_strip_env(ref_p$covariance_handling), info = info)
    expect_identical(res$gmm_fit$provenance$reduction$kept_indices, ref_p$gmm_fit$provenance$reduction$kept_indices, info = info)
    expect_identical(res$gmm_fit$provenance$reduction$n_links_all, 190L, info = info)
    expect_identical(res$gmm_fit$provenance$nuisance_fit_id, ref_p$gmm_fit$provenance$nuisance_fit_id, info = info)
    expect_true(is.na(res$reduction_summary$max_moment_expansion_diff), info = info)
    expect_true(is.finite(ref_p$reduction_summary$max_moment_expansion_diff), info = info)
    expect_equal(as.double(ref_p$reduction_summary$expanded_m),
      unname(vapply(ref_p$point_links_all, `[[`, numeric(1L), "estimate")), tolerance = 1e-12, info = info)
  }
})

# ---------------------------------------------------------------------------
# 2. Actual constructor / canonical-fact / certification counts and ordering
# ---------------------------------------------------------------------------

test_that("compact route: 190 preserved, heavy views only 19 plus auxiliary, one catalog/mapping construction, caps before allocation", {
  d <- cb_panel(41L, 20L)
  n_frame <- length(unique(d$id))
  wrap <- c("construct_ipw_link_view", "construct_ipw_nuisance_view", "fit_ipw_nuisance_view",
    "build_comparison_catalog", "pipeline_metadata_resolve_catalog_mappings", "derive_owner_sparse_occupancy",
    "build_pipeline_nuisance_and_link_views", "certify_ipw_structural_classes_from_views",
    "preflight_ipw_early_structural_reduction", "certify_ipw_structural_classes_from_descriptors",
    "ipw_pipeline_compact_basis_preflight", "ipw_pipeline_compact_basis_facts", "verify_ipw_view_mask_compatibility",
    "ipw_link_view_validate_body", "fingerprint_ipw_view_mask", "ipw_owner_date_complete",
    "build_owner_comparison_catalog")
  for (xf in list(~1, ~x)) for (chunk in c(7L, 16L)) {
    info <- paste(deparse(xf), chunk)
    ctrl <- cb_public_ctrl(d, n_frame, chunk = chunk, n_covariates = length(all.vars(xf)) + 1)
    e <- cb_counted(wrap)
    g <- cb_general(wrap)
    trace <- new.env(parent = emptyenv())
    # The compact preflight must run before any view, fit or seal of this registry.
    pre_inner <- get("ipw_pipeline_compact_basis_preflight", envir = e, inherits = FALSE)
    assign("ipw_pipeline_compact_basis_preflight", function(...) {
      trace$at_preflight <- cb_snapshot(e)
      pre_inner(...)
    }, envir = e)
    # Capture the first catalog object to prove the design pass reused it.
    build_inner <- get("build_owner_comparison_catalog", envir = e, inherits = FALSE)
    assign("build_owner_comparison_catalog", function(...) {
      out <- build_inner(...)
      if (is.null(trace$first_meta)) trace$first_meta <- out
      out
    }, envir = e)
    registry_inner <- get("build_pipeline_compact_basis_registry", envir = e, inherits = FALSE)
    assign("build_pipeline_compact_basis_registry", function(...) {
      trace$before_registry <- cb_snapshot(e)
      out <- registry_inner(...)
      trace$after_registry <- cb_snapshot(e)
      out
    }, envir = e)
    res <- cb_run(e, d, ctrl, xformla = xf)
    ref <- cb_run(g, d, ctrl, xformla = xf)
    for (f in cb_fields) cb_same(res[[f]], ref[[f]], paste(info, f))

    # Ordering: preflight before any allocation of the route.
    at <- trace$at_preflight
    expect_false(is.null(at), info = info)
    for (nm in c("construct_ipw_link_view", "construct_ipw_nuisance_view", "fit_ipw_nuisance_view",
      "certify_ipw_structural_classes_from_descriptors")) expect_identical(at[[nm]], 0, info = paste(info, nm))
    expect_true(is.null(at[["seal:link_view"]]) || at[["seal:link_view"]] == 0, info = info)
    expect_true(is.null(at[["seal:nuisance_view"]]) || at[["seal:nuisance_view"]] == 0, info = info)
    expect_true(is.null(at[["seal:fitted_nuisance_view"]]) || at[["seal:fitted_nuisance_view"]] == 0, info = info)
    expect_identical(at[["ipw_pipeline_compact_basis_facts"]], 1, info = info)

    # One catalog, one occupancy and one CSR mapping construction on both routes;
    # the design metadata reused the first sealed catalog.
    for (ns in list(e, g)) {
      expect_identical(cb_n(ns, "build_owner_comparison_catalog"), 2, info = info)
      expect_identical(cb_n(ns, "build_comparison_catalog"), 1, info = info)
      expect_identical(cb_n(ns, "pipeline_metadata_resolve_catalog_mappings"), 1, info = info)
      expect_identical(cb_n(ns, "derive_owner_sparse_occupancy"), 1, info = info)
    }
    expect_true(isTRUE(res$metadata$diagnostics$catalog_reuse$reused), info = info)
    expect_identical(res$metadata$catalog_source_ref, trace$first_meta$ref, info = info)
    expect_true(api$ipw_owner_seal_ok(trace$first_meta, "pipeline_metadata"), info = info)
    expect_identical(res$metadata$catalog, trace$first_meta$catalog, info = info)
    expect_identical(res$metadata$row_mapping, trace$first_meta$row_mapping, info = info)
    expect_identical(res$metadata$operator_mapping, trace$first_meta$operator_mapping, info = info)
    expect_identical(res$metadata$occupancy_ref, trace$first_meta$occupancy_ref, info = info)
    expect_false(is.null(res$metadata$design), info = info)
    expect_null(trace$first_meta$design, info = info)
    expect_identical(res$reduction_summary$catalog_constructions, 1L, info = info)

    # Compact route: no view registry, no view certificate, no streamed mask,
    # one descriptor certificate; the general route keeps its own counts.
    expect_identical(cb_n(e, "build_pipeline_nuisance_and_link_views"), 0, info = info)
    expect_identical(cb_n(e, "certify_ipw_structural_classes_from_views"), 0, info = info)
    expect_identical(cb_n(e, "preflight_ipw_early_structural_reduction"), 0, info = info)
    expect_identical(cb_n(e, "certify_ipw_structural_classes_from_descriptors"), 1, info = info)
    expect_identical(cb_n(e, "ipw_pipeline_compact_basis_preflight"), 1, info = info)
    expect_identical(trace$after_registry[["verify_ipw_view_mask_compatibility"]] -
      trace$before_registry[["verify_ipw_view_mask_compatibility"]], 0, info = info)
    expect_identical(cb_n(e, "fingerprint_ipw_view_mask"), 0, info = info)
    expect_true(cb_n(e, "ipw_owner_date_complete") >= 20, info = info)
    expect_identical(cb_n(g, "build_pipeline_nuisance_and_link_views"), 1, info = info)
    expect_identical(cb_n(g, "certify_ipw_structural_classes_from_views"), 1, info = info)
    expect_identical(cb_n(g, "certify_ipw_structural_classes_from_descriptors"), 0, info = info)

    # Heavy views: 19 retained basis rows in the main registry (190 on the
    # general route); the auxiliary work is unchanged on both routes, so the
    # whole-run link constructions differ by exactly the 171 dropped rows and
    # the nuisance constructions/fits agree.
    cw <- res$nuisance_registry$construction_work
    expect_identical(cw$route, "compact_canonical_basis", info = info)
    expect_identical(cw$n_rows, 190L, info = info)
    expect_identical(cw$link_view_constructions, 19L, info = info)
    expect_identical(cw$heavy_link_views_constructed, 19L, info = info)
    expect_identical(cw$canonical_descriptor_rows, 190L, info = info)
    expect_identical(cw$dropped_rows_without_link_view, 171L, info = info)
    expect_identical(cw$n_distinct_nuisance_views, 1L, info = info)
    expect_identical(ref$nuisance_registry$construction_work$n_rows, 190L, info = info)
    expect_identical(cb_n(e, "construct_ipw_link_view"), cb_n(g, "construct_ipw_link_view") - 171, info = info)
    expect_identical(cb_n(e, "construct_ipw_nuisance_view"), cb_n(g, "construct_ipw_nuisance_view"), info = info)
    expect_identical(cb_n(e, "fit_ipw_nuisance_view"), cb_n(g, "fit_ipw_nuisance_view"), info = info)
    expect_identical(cb_n(e, "seal:link_view"), cb_n(g, "seal:link_view") - 171, info = info)
    expect_identical(cb_n(e, "seal:link_descriptor"), 190, info = info)
    expect_identical(cb_n(e, "seal:nuisance_view"), cb_n(g, "seal:nuisance_view"), info = info)
    Q_aux <- res$resource_accounting$auxiliary_plan$n_descriptors_union
    expect_identical(cb_n(e, "construct_ipw_link_view"), 19 + Q_aux, info = info)
    expect_identical(cb_n(e, "seal:link_view"), 19 + Q_aux, info = info)
    expect_identical(cb_n(e, "fit_ipw_nuisance_view"), 1, info = info)
    expect_identical(length(Filter(function(v) inherits(v, "cdid_ipw_link_view"), res$nuisance_registry$link_views)), 19L, info = info)
    expect_identical(length(res$nuisance_registry$link_views), 190L, info = info)
    expect_identical(names(res$nuisance_registry$link_views), names(ref$nuisance_registry$link_views), info = info)
    expect_identical(res$nuisance_registry$nuisance_keys, ref$nuisance_registry$nuisance_keys, info = info)
    cb_same(res$nuisance_registry$nuisance_views[[1L]], ref$nuisance_registry$nuisance_views[[1L]], paste(info, "nuisance view"))
    cb_same(res$nuisance_registry$fitted_views[[1L]], ref$nuisance_registry$fitted_views[[1L]], paste(info, "fitted view"))

    # Mask facts: every catalog date complete, re-read in the certificate itself.
    mf <- res$nuisance_registry$compact_basis$mask_facts
    expect_true(api$ipw_seal_ok(mf, "mask_facts"), info = info)
    expect_identical(mf$date_indices_complete, 1:20, info = info)
    expect_identical(mf$n_date_facts_reread, 20L, info = info)
    expect_identical(mf$n_rows, 190, info = info)

    # Resource accounting: the reservation is the smaller compact one (no
    # original link views), checked against the unchanged caps.
    es <- res$resource_accounting$early_structural
    expect_identical(es$route, "compact_canonical_basis", info = info)
    expect_identical(es$original_link_view_entries, 0, info = info)
    expect_identical(es$retained_basis_bound, 19, info = info)
    expect_identical(es$mask_streams_charged, 0, info = info)
    ref_es <- ref$resource_accounting$early_structural
    expect_lt(es$additional_entries, ref_es$additional_entries)
    expect_lt(res$resource_accounting$initial_envelope_entries, ref$resource_accounting$initial_envelope_entries)
    expect_identical(res$controls$pipeline_control, ref$controls$pipeline_control, info = info)
  }
})

# ---------------------------------------------------------------------------
# 3. Catalog reuse: exact arguments and sealed owners; stale/forged/rebound refuse
# ---------------------------------------------------------------------------

test_that("catalog reuse binds exact original arguments and sealed owners; stale, forged, rebound and foreign sources refuse", {
  d <- cb_panel(41L, 8L)
  n_frame <- length(unique(d$id))
  pub <- cb_public_ctrl(d, n_frame, chunk = 7L)
  oc <- pub$owner_control; cb <- pub$catalog_bounds
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", oc)
  u <- d[!duplicated(d$id), ]; u <- u[order(u$id), ]
  fo <- api$construct_ipw_frame_owner(po, u$id, matrix(1, nrow(u), 1L), weights = NULL, control = oc)
  fo2 <- api$construct_ipw_frame_owner(po, u$id, cbind(1, u$x), weights = NULL, control = oc)
  tg <- 3
  build <- function(ns = api, ...) ns$build_owner_comparison_catalog(panel_owner = po, targets = tg, anticipation = 0L,
    control_group = "nevertreated", conventions = "D10_D13_admissible_frame", requested_post = NULL,
    control = oc, catalog_bounds = cb, frame_owner = fo, ...)
  meta1 <- build(build_design = FALSE)
  expect_true(api$ipw_owner_seal_ok(meta1, "pipeline_metadata"))
  expect_false(meta1$diagnostics$catalog_reuse$reused)
  expect_identical(meta1$construction$targets, tg)
  # Design inputs derived as the pipeline derives them from the first catalog.
  c_rows <- meta1$catalog$rows
  s_c <- c_rows$s_role %in% c("anticipation_conditional", "post_conditional")
  t_c <- c_rows$t_role %in% c("anticipation_conditional", "post_conditional")
  dq <- unique(data.frame(g = as.double(c(c_rows$g[s_c], c_rows$g[t_c])), t = as.double(c(c_rows$s[s_c], c_rows$t[t_c])),
    stringsAsFactors = FALSE))
  dq <- dq[order(dq$g, dq$t), , drop = FALSE]; rownames(dq) <- NULL
  tpt <- data.frame(g = 3, t = as.double(3:8), stringsAsFactors = FALSE)
  build <- function(ns = api, ...) ns$build_owner_comparison_catalog(panel_owner = po, targets = tg, anticipation = 0L,
    control_group = "nevertreated", conventions = "D10_D13_admissible_frame", requested_post = NULL,
    control = oc, catalog_bounds = cb, frame_owner = fo, declared_qualified = dq, ...)
  e <- cb_counted(c("build_comparison_catalog", "pipeline_metadata_resolve_catalog_mappings", "derive_owner_sparse_occupancy"))
  meta2 <- build(ns = e, occupancy = meta1$occupancy_object, target_post_times = tpt, build_design = TRUE, metadata_source = meta1)
  meta2_ref <- build(occupancy = meta1$occupancy_object, target_post_times = tpt, build_design = TRUE)
  expect_identical(cb_n(e, "build_comparison_catalog"), 0)
  expect_identical(cb_n(e, "pipeline_metadata_resolve_catalog_mappings"), 0)
  expect_identical(cb_n(e, "derive_owner_sparse_occupancy"), 0)
  expect_true(api$ipw_owner_seal_ok(meta2, "pipeline_metadata"))
  expect_true(meta2$diagnostics$catalog_reuse$reused)
  expect_identical(meta2$catalog_source_ref, meta1$ref)
  expect_identical(meta2$catalog, meta1$catalog)
  expect_identical(meta2$row_mapping, meta1$row_mapping)
  expect_identical(meta2$occupancy_ref, meta1$occupancy_ref)
  for (f in c("catalog", "row_mapping", "operator_mapping", "design_links", "design", "occupancy", "cell_indices",
    "owner_summary", "control", "construction")) cb_same(meta2[[f]], meta2_ref[[f]], f)
  expect_identical(meta2$diagnostics$catalog_reuse$borrowed_catalog_mapping_entries,
    meta2_ref$diagnostics$catalog_reuse$owned_catalog_mapping_entries)
  expect_identical(meta2$diagnostics$catalog_reuse$owned_design_entries, meta2_ref$diagnostics$catalog_reuse$owned_design_entries)
  expect_gt(meta2$diagnostics$catalog_reuse$owned_design_entries, 0)
  # A source without the occupancy argument: the source's sealed occupancy is validated and reused.
  meta3 <- build(target_post_times = tpt, build_design = TRUE, metadata_source = meta1)
  cb_same(meta3$design, meta2$design, "design without occupancy arg")
  expect_identical(meta3$occupancy_ref, meta1$occupancy_ref)

  # Different exact arguments refuse (type, convention and anticipation), typed.
  cb_expect_class(cb_err(api$build_owner_comparison_catalog(panel_owner = po, targets = 3L, anticipation = 0L,
    control_group = "nevertreated", conventions = "D10_D13_admissible_frame", control = oc, catalog_bounds = cb,
    frame_owner = fo, declared_qualified = dq, build_design = TRUE, target_post_times = tpt, metadata_source = meta1)),
    "cdid_pipeline_metadata_reuse_error", "targets type")
  cb_expect_class(cb_err(api$build_owner_comparison_catalog(panel_owner = po, targets = tg, anticipation = 0L,
    control_group = "nevertreated", conventions = "other", control = oc, catalog_bounds = cb,
    frame_owner = fo, declared_qualified = dq, build_design = TRUE, target_post_times = tpt, metadata_source = meta1)),
    "cdid_pipeline_metadata_reuse_error", "conventions")
  cb_expect_class(cb_err(api$build_owner_comparison_catalog(panel_owner = po, targets = tg, anticipation = 1L,
    control_group = "nevertreated", conventions = "D10_D13_admissible_frame", control = oc, catalog_bounds = cb,
    frame_owner = fo, declared_qualified = dq, build_design = TRUE, target_post_times = tpt, metadata_source = meta1)),
    "cdid_pipeline_metadata_reuse_error", "anticipation")
  # Stale and forged sources refuse through the original seal validator.
  stale <- meta1; stale$catalog$rows$s[1L] <- stale$catalog$rows$s[1L] + 0.5
  cb_expect_class(cb_err(build(target_post_times = tpt, build_design = TRUE, metadata_source = stale)),
    "cdid_pipeline_metadata_integrity_error", "stale")
  forged <- meta1; forged$ref <- new.env(parent = emptyenv())
  cb_expect_class(cb_err(build(target_post_times = tpt, build_design = TRUE, metadata_source = forged)),
    "cdid_pipeline_metadata_integrity_error", "forged")
  nocons <- meta1; nocons$construction <- NULL
  cb_expect_class(cb_err(build(target_post_times = tpt, build_design = TRUE, metadata_source = nocons)),
    "cdid_pipeline_metadata_integrity_error", "construction removed")
  cb_expect_class(cb_err(build(target_post_times = tpt, build_design = TRUE, metadata_source = list(a = 1))),
    "cdid_pipeline_metadata_integrity_error", "malformed")
  # Foreign frame and rebound occupancy refuse.
  cb_expect_class(cb_err(api$build_owner_comparison_catalog(panel_owner = po, targets = tg, anticipation = 0L,
    control_group = "nevertreated", conventions = "D10_D13_admissible_frame", control = oc, catalog_bounds = cb,
    frame_owner = fo2, declared_qualified = dq, build_design = TRUE, target_post_times = tpt, metadata_source = meta1)),
    "cdid_pipeline_metadata_owner_mismatch_error", "foreign frame")
  occ2 <- api$derive_owner_sparse_occupancy(panel_owner = po, control = oc, frame_owner = fo)
  cb_expect_class(cb_err(build(occupancy = occ2, target_post_times = tpt, build_design = TRUE, metadata_source = meta1)),
    "cdid_pipeline_metadata_owner_mismatch_error", "rebound occupancy")
  # Malformed arguments keep their original typed refusals with or without a source.
  a <- cb_err(build(target_post_times = tpt, build_design = TRUE, metadata_source = meta1, pt_assumption = "bogus"))
  b <- cb_err(build(occupancy = meta1$occupancy_object, target_post_times = tpt, build_design = TRUE, pt_assumption = "bogus"))
  expect_true(inherits(a, "error")); expect_identical(class(a), class(b)); expect_identical(conditionMessage(a), conditionMessage(b))
  a <- cb_err(api$build_owner_comparison_catalog(panel_owner = po, targets = "three", anticipation = 0L,
    control_group = "nevertreated", conventions = "D10_D13_admissible_frame", control = oc, catalog_bounds = cb,
    frame_owner = fo, declared_qualified = dq, build_design = TRUE, target_post_times = tpt, metadata_source = meta1))
  b <- cb_err(api$build_owner_comparison_catalog(panel_owner = po, targets = "three", anticipation = 0L,
    control_group = "nevertreated", conventions = "D10_D13_admissible_frame", control = oc, catalog_bounds = cb,
    frame_owner = fo, declared_qualified = dq, build_design = TRUE, target_post_times = tpt, occupancy = meta1$occupancy_object))
  expect_true(inherits(a, "error")); expect_identical(class(a), class(b)); expect_identical(conditionMessage(a), conditionMessage(b))
})

# ---------------------------------------------------------------------------
# 4. Compact certificate: genuine facts only; mutations refuse; other contexts
#    keep their routes; permutations and floating dates agree
# ---------------------------------------------------------------------------

test_that("compact certificate refuses stale/forged owner, metadata, nuisance and mask mutations; incomplete and other contexts keep the general route", {
  d <- cb_panel(41L, 12L)
  n_frame <- length(unique(d$id))
  ctrl <- cb_public_ctrl(d, n_frame, chunk = 7L)
  res <- cb_run(api, d, ctrl)
  expect_identical(res$reduction_summary$reduction_stage, "outcome_free_compact_canonical_basis")
  M <- res$structural_certificate$n_links
  expect_identical(M, 66L)
  reg <- res$nuisance_registry
  cbnd <- res$controls$pipeline_control$class_bounds
  oc <- res$controls$pipeline_control$owner_control
  dli <- match(seq_len(M), res$design_all$links$catalog_row)
  certify <- function(md = res$metadata, fo = res$frame_owner, nvs = reg$nuisance_views, fvs = reg$fitted_views,
                      idx = rep(1L, M), design = res$design_all, bounds = cbnd, control = oc)
    api$certify_ipw_structural_classes_from_descriptors(md, fo, nvs, fvs, idx, design, dli, bounds, control)
  again <- certify()
  for (f in cb_cert_fields) cb_same(again$certificate[[f]], res$structural_certificate[[f]], f)
  expect_identical(again$n_date_facts, 12L)
  expect_length(again$descriptors, M)

  # Facts: the route predicate on the genuine objects.
  facts <- api$ipw_pipeline_compact_basis_facts(res$metadata, res$frame_owner, oc, "post", FALSE, n_frame)
  expect_false(is.null(facts))
  expect_identical(facts$g, 3); expect_identical(facts$K, Inf); expect_identical(facts$n_rows, M)
  expect_identical(facts$date_indices, 1:12); expect_identical(facts$first_row, 1L)
  expect_null(api$ipw_pipeline_compact_basis_facts(res$metadata, res$frame_owner, oc, "all", FALSE, n_frame))
  expect_null(api$ipw_pipeline_compact_basis_facts(res$metadata, res$frame_owner, oc, "post", TRUE, n_frame))
  expect_null(api$ipw_pipeline_compact_basis_facts(res$metadata, res$frame_owner, oc, "post", FALSE, n_frame + 1L))

  # Caps are checked before any allocation: a tight retained or work cap refuses
  # the preflight itself with the owner bound error, constructing nothing.
  e <- cb_counted(c("construct_ipw_link_view", "construct_ipw_nuisance_view", "fit_ipw_nuisance_view"))
  tight_ret <- api$ipw_owner_control(max_chunk_rows = 7, max_retained_entries = 2000, max_matrix_entries = 1e7,
    max_work_ops = 1e10, max_output_entries = 1e7)
  cb_expect_class(cb_err(e$ipw_pipeline_compact_basis_preflight(facts, res$design_all, cbnd, tight_ret, 0, 0)),
    "cdid_ipw_owner_bound_error", "tight retained")
  tight_work <- api$ipw_owner_control(max_chunk_rows = 7, max_retained_entries = 1e7, max_matrix_entries = 1e7,
    max_work_ops = 5000, max_output_entries = 1e7)
  cb_expect_class(cb_err(e$ipw_pipeline_compact_basis_preflight(facts, res$design_all, cbnd, tight_work, 0, 0)),
    "cdid_ipw_owner_bound_error", "tight work")
  ok <- api$ipw_pipeline_compact_basis_preflight(facts, res$design_all, cbnd, oc, 0, 0)
  expect_identical(ok$retained_basis_bound, 11)
  cb_expect_class(cb_err(api$ipw_pipeline_compact_basis_preflight(facts, res$design_all, cbnd, oc,
    oc$max_retained_entries - 1, 0)), "cdid_ipw_owner_bound_error", "adjacent live envelope")
  for (nm in c("construct_ipw_link_view", "construct_ipw_nuisance_view", "fit_ipw_nuisance_view")) expect_identical(cb_n(e, nm), 0)

  # Stale/forged nuisance, fitted model, metadata and owner mutations refuse.
  nv_bad <- reg$nuisance_views[[1L]]; nv_bad$support$n_risk_a_units <- 0L
  cb_expect_class(cb_err(certify(nvs = list(nv_bad))), "cdid_ipw_owner_integrity_error", "stale nuisance")
  fv_bad <- reg$fitted_views[[1L]]; fv_bad$coefficients[1L] <- fv_bad$coefficients[1L] + 1
  cb_expect_class(cb_err(certify(fvs = list(fv_bad))), "cdid_ipw_fitted_view_integrity_error", "forged fit")
  md_bad <- res$metadata; md_bad$catalog$rows$s[2L] <- md_bad$catalog$rows$s[2L] + 0.5
  cb_expect_class(cb_err(certify(md = md_bad)), "cdid_pipeline_metadata_integrity_error", "stale metadata")
  md_rows <- res$metadata; md_rows$row_mapping$functional_key[1L] <- "forged"
  cb_expect_class(cb_err(certify(md = md_rows)), "cdid_pipeline_metadata_integrity_error", "forged mapping")
  res_x <- cb_run(api, d, cb_public_ctrl(d, n_frame, chunk = 7L, n_covariates = 2), xformla = ~x)
  foreign <- cb_err(certify(fo = res_x$frame_owner))
  expect_true(inherits(foreign, "error"))
  expect_true(inherits(foreign, "cdid_pipeline_metadata_owner_mismatch_error") || inherits(foreign, "cdid_ipw_owner_mismatch_error"),
    info = paste(class(foreign), collapse = "/"))
  cb_expect_class(cb_err(certify(nvs = res_x$nuisance_registry$nuisance_views, fvs = res_x$nuisance_registry$fitted_views)),
    "cdid_ipw_owner_mismatch_error", "foreign nuisance")
  cb_expect_class(cb_err(certify(idx = rep(2L, M))), "cdid_ipw_input_error", "row map out of range")
  bad_design <- res$design_all; bad_design$links$s[3L] <- bad_design$links$s[3L] + 0.25
  cb_expect_class(cb_err(certify(design = bad_design)), "cdid_ipw_certificate_error", "mutated design")

  # Incomplete roster: two legacy-route clones agree bitwise on the general
  # early route (`ipw_pipeline_compact_basis_facts` returns NULL on incomplete
  # dates), the compact certificate refuses because the mask facts are not
  # established (never a silent fallback inside), and the production namespace
  # takes the compact standalone raw-IPW backend and agrees at 1e-12 on all
  # stable scientific outputs (while `z_stat == estimates / se` holds on each
  # route, including the two deterministic unit-change links with SE ~ 1e-15).
  cb_match_1e12 <- function(a, b, info = NULL) {
    a <- cb_strip_env(a); b <- cb_strip_env(b)
    if (is.numeric(a) && is.numeric(b) && !is.object(a) && !is.object(b)) {
      expect_identical(dim(a), dim(b), info = info)
      expect_identical(dimnames(a), dimnames(b), info = info)
      expect_identical(names(a), names(b), info = info)
      expect_identical(is.na(a), is.na(b), info = info)
      ok <- !is.na(a)
      if (any(ok)) {
        expect_identical(is.infinite(a[ok]), is.infinite(b[ok]), info = info)
        fin <- ok & is.finite(a) & is.finite(b)
        if (any(fin)) {
          sc <- pmax(1, abs(a[fin]), abs(b[fin]))
          expect_true(max(abs(a[fin] - b[fin]) / sc) <= 1e-12, info = info)
        }
      }
      return(invisible(NULL))
    }
    if (is.data.frame(a) && is.data.frame(b)) {
      expect_identical(names(a), names(b), info = info)
      expect_identical(nrow(a), nrow(b), info = info)
      for (nm in names(a)) cb_match_1e12(a[[nm]], b[[nm]], paste(info, nm))
      return(invisible(NULL))
    }
    if (is.list(a) && is.list(b)) {
      expect_identical(class(a), class(b), info = info)
      expect_identical(names(a), names(b), info = info)
      expect_identical(length(a), length(b), info = info)
      for (i in seq_along(a)) cb_match_1e12(a[[i]], b[[i]], paste(info, names(a)[i]))
      return(invisible(NULL))
    }
    expect_identical(a, b, info = info)
  }
  d_inc <- cb_panel(41L, 12L, drop = c("4 5"))
  n_inc <- length(unique(d_inc$id))
  ctrl_inc <- cb_public_ctrl(d_inc, n_inc, chunk = 7L)
  legacy_ns <- cb_counted()
  general_ns <- cb_general()
  legacy_inc <- cb_run(legacy_ns, d_inc, ctrl_inc)
  ref_inc <- cb_run(general_ns, d_inc, ctrl_inc)
  expect_identical(legacy_inc$reduction_summary$reduction_stage, "outcome_free_before_point_fits")
  expect_equal(legacy_inc$reduction_summary$n_link_views_constructed, legacy_inc$structural_certificate$n_links)
  for (f in cb_fields) cb_same(legacy_inc[[f]], ref_inc[[f]], paste("incomplete legacy same-route", f))

  res_inc <- cb_run(api, d_inc, ctrl_inc)
  expect_s3_class(res_inc$nuisance_registry$raw_ipw, "cdid_ipw_raw_backend")
  expect_null(ref_inc$nuisance_registry$raw_ipw)
  expect_identical(res_inc$reduction_summary$reduction_stage, "outcome_free_before_point_fits")
  expect_equal(res_inc$reduction_summary$n_link_views_constructed, res_inc$structural_certificate$n_links)
  expect_null(res_inc$resource_accounting$point_general)
  expect_false(is.null(ref_inc$resource_accounting$point_general))
  expect_identical(res_inc$resource_accounting$point_producer$retained_moment_if$status, "retained")
  cb_same(res_inc$metadata$catalog$rows, ref_inc$metadata$catalog$rows, "incomplete catalog rows")
  cb_same(res_inc$support, ref_inc$support, "incomplete support")
  cb_same(res_inc$basis, ref_inc$basis, "incomplete basis")
  cb_same(res_inc$design_all, ref_inc$design_all, "incomplete design_all")
  cb_same(res_inc$design_used, ref_inc$design_used, "incomplete design_used")
  cb_same(res_inc$reporting$target_meta, ref_inc$reporting$target_meta, "incomplete reporting target_meta")
  cb_same(res_inc$reporting$available, ref_inc$reporting$available, "incomplete reporting available")
  cb_same(res_inc$reporting$support, ref_inc$reporting$support, "incomplete reporting support")
  cb_same(res_inc$structural_certificate$retained_index, ref_inc$structural_certificate$retained_index, "incomplete retained_index")
  cb_same(res_inc$structural_certificate$dropped_index, ref_inc$structural_certificate$dropped_index, "incomplete dropped_index")
  cb_same(res_inc$structural_certificate$classes, ref_inc$structural_certificate$classes, "incomplete classes")
  cb_same(res_inc$structural_certificate$E, ref_inc$structural_certificate$E, "incomplete E")
  cb_same(res_inc$structural_certificate$weighting, ref_inc$structural_certificate$weighting, "incomplete weighting")
  target_identity <- c("target_id", "estimable", "status", "support_reason")
  cb_same(res_inc$gmm_fit$targets[target_identity], ref_inc$gmm_fit$targets[target_identity], "incomplete target identification")
  penalty_identity <- c("weighting_requested", "weighting_canonical", "handler_selected", "handler_invoked",
                        "condition_trigger", "ridge_relative", "target_availability")
  cb_same(res_inc$covariance_handling[penalty_identity], ref_inc$covariance_handling[penalty_identity], "incomplete weighting and penalty policy")
  cb_same(res_inc$covariance_handling$handler_used, ref_inc$covariance_handling$handler_used, "incomplete handler_used")
  cb_expect_incomplete_lambda_relative <- function(x, y) {
    schema_ok <- typeof(x) == "double" && typeof(y) == "double" &&
      length(x) == 1L && length(y) == 1L &&
      identical(attributes(x), attributes(y), num.eq = FALSE)
    expect_true(schema_ok, info = "incomplete lambda")
    if (!schema_ok) return(invisible(NULL))
    finite_ok <- all(is.finite(x)) && all(is.finite(y))
    expect_true(finite_ok, info = "incomplete lambda")
    if (!finite_ok) return(invisible(NULL))
    expect_identical(x == 0, y == 0, info = "incomplete lambda zero status")
    expect_identical(sign(x), sign(y), info = "incomplete lambda sign")
    expect_identical(1 / x[x == 0], 1 / y[y == 0], info = "incomplete lambda zero sign")
    if (x != 0 && y != 0)
      expect_true(abs(x - y) / max(abs(x), abs(y)) <= 1e-12, info = "incomplete lambda relative error")
  }
  cb_expect_incomplete_lambda_relative(res_inc$covariance_handling$lambda, ref_inc$covariance_handling$lambda)
  expect_identical(res_inc$reporting$z_stat, res_inc$reporting$estimates / res_inc$reporting$se)
  expect_identical(ref_inc$reporting$z_stat, ref_inc$reporting$estimates / ref_inc$reporting$se)
  # FIXED_INPUT_PARITY.json also records the same ill-conditioned z values in
  # inference/compact and in the dynamic/calendar aggregation copies. Validate
  # each copy's route-local algebra; compare every other aggregation field.
  cb_stable_inference <- function(x, info) {
    expect_identical(x$z_stat, x$estimates / x$se, info = info)
    expect_identical(x$compact$z_stat, x$compact$estimates / x$compact$se, info = paste(info, "compact"))
    x$z_stat <- NULL
    x$compact$z_stat <- NULL
    x
  }
  invisible(cb_stable_inference(res_inc$reporting$inference, "incomplete raw reporting inference"))
  invisible(cb_stable_inference(ref_inc$reporting$inference, "incomplete legacy reporting inference"))
  cb_same(res_inc$covariance_stream$n_a_units, ref_inc$covariance_stream$n_a_units, "incomplete n_a_units")
  cb_match_1e12(res_inc$covariance_stream$mhat, ref_inc$covariance_stream$mhat, "incomplete mhat")
  cb_match_1e12(res_inc$covariance_stream$V, ref_inc$covariance_stream$V, "incomplete V")
  cb_match_1e12(as.double(res_inc$covariance_stream$mean_phi), as.double(ref_inc$covariance_stream$mean_phi), "incomplete mean_phi")
  cb_match_1e12(res_inc$gmm_fit$targets, ref_inc$gmm_fit$targets, "incomplete gmm_fit targets")
  cb_match_1e12(res_inc$gmm_fit$if_projection$sandwich_cov, ref_inc$gmm_fit$if_projection$sandwich_cov, "incomplete sandwich_cov")
  cb_match_1e12(res_inc$gmm_fit$if_projection$materialized_target_if, ref_inc$gmm_fit$if_projection$materialized_target_if, "incomplete materialized_target_if")
  cb_match_1e12(res_inc$full_frame_transport, ref_inc$full_frame_transport, "incomplete full_frame_transport")
  for (nm in c("estimates", "se", "analytical_se", "p_value", "ci_lower", "ci_upper",
               "cband_crit_val", "cband_lower", "cband_upper", "vcov", "Omega_a", "if_full")) {
    cb_match_1e12(res_inc$reporting[[nm]], ref_inc$reporting[[nm]], paste("incomplete reporting", nm))
  }
  agg_inc <- res_inc$aggregations; agg_ref <- ref_inc$aggregations
  for (nm in c("dynamic", "calendar")) {
    agg_inc[[nm]]$inference <- cb_stable_inference(agg_inc[[nm]]$inference, paste("incomplete raw aggregation", nm))
    agg_ref[[nm]]$inference <- cb_stable_inference(agg_ref[[nm]]$inference, paste("incomplete legacy aggregation", nm))
  }
  cb_match_1e12(agg_inc, agg_ref, "incomplete aggregations")
  cb_match_1e12(res_inc$auxiliary_reporting$estimates, ref_inc$auxiliary_reporting$estimates, "incomplete auxiliary estimates")
  cb_same(res_inc$auxiliary_reporting$selections, ref_inc$auxiliary_reporting$selections, "incomplete auxiliary selections")
  cb_same(res_inc$auxiliary_reporting$descriptors, ref_inc$auxiliary_reporting$descriptors, "incomplete auxiliary descriptors")
  cb_same(res_inc$auxiliary_reporting$descriptor_status, ref_inc$auxiliary_reporting$descriptor_status, "incomplete auxiliary descriptor_status")
  expect_null(api$ipw_pipeline_compact_basis_facts(res_inc$metadata, res_inc$frame_owner, oc, "post", FALSE, n_inc))
  M_inc <- res_inc$structural_certificate$n_links
  cb_expect_class(cb_err(api$certify_ipw_structural_classes_from_descriptors(res_inc$metadata, res_inc$frame_owner,
    res_inc$nuisance_registry$nuisance_views[1L], res_inc$nuisance_registry$fitted_views[1L], rep(1L, M_inc),
    res_inc$design_all, match(seq_len(M_inc), res_inc$design_all$links$catalog_row), cbnd, oc)),
    "cdid_ipw_shared_mask_error", "incomplete mask facts")

  d_two <- d; d_two$w <- 2
  two <- cb_run(api, d_two, ctrl, weightsname = "w")
  two_ref <- cb_run(cb_general(), d_two, ctrl, weightsname = "w")
  expect_identical(two$reduction_summary$reduction_stage, "outcome_free_before_point_fits")
  for (f in cb_fields) cb_same(two[[f]], two_ref[[f]], paste("constant nonunit", f))

  # Nonconstant weights, identity weighting, unreduced requests and the
  # post-fit stage keep their existing routes unchanged by the compact code.
  for (case in list(list(args = list(weightsname = "w"), stage = "outcome_free_before_point_fits"),
                    list(args = list(weighting = "identity"), stage = "none"),
                    list(args = list(reduce_structural_redundancy = FALSE), stage = "after_point_fits"),
                    list(args = list(stage = "post_point_fit"), stage = "after_point_fits"))) {
    ctrl_c <- cb_public_ctrl(d, n_frame, chunk = 7L, weighting = if (is.null(case$args$weighting)) "optimal" else case$args$weighting)
    r1 <- tryCatch(do.call(cb_run, c(list(ns = api, d = d, ctrl = ctrl_c), case$args)), error = identity)
    r2 <- tryCatch(do.call(cb_run, c(list(ns = cb_general(), d = d, ctrl = ctrl_c), case$args)), error = identity)
    if (inherits(r1, "error")) {
      expect_true(inherits(r2, "error")); expect_identical(class(r1), class(r2)); expect_identical(conditionMessage(r1), conditionMessage(r2))
      next
    }
    expect_identical(r1$reduction_summary$reduction_stage, case$stage, info = case$stage)
    for (f in intersect(cb_fields, names(r1))) cb_same(r1[[f]], r2[[f]], paste(case$stage, f))
  }
})

test_that("compact canonical basis: row/ID permutations and adjacent floating dates agree with the general early route", {
  dates <- 0.1 * (1:12) + 1e-9 * (1:12)
  for (variant in list(list(id_perm = TRUE, rows_perm = TRUE, dates = NULL),
                       list(id_perm = FALSE, rows_perm = TRUE, dates = dates),
                       list(id_perm = TRUE, rows_perm = FALSE, dates = dates))) {
    d <- cb_panel(41L, 12L, dates = variant$dates, id_perm = variant$id_perm, rows_perm = variant$rows_perm)
    n_frame <- length(unique(d$id))
    for (chunk in c(7L, 16L)) {
      info <- paste(variant$id_perm, variant$rows_perm, !is.null(variant$dates), chunk)
      ctrl <- cb_public_ctrl(d, n_frame, chunk = chunk)
      res <- tryCatch(cb_run(api, d, ctrl), error = identity)
      ref <- tryCatch(cb_run(cb_general(), d, ctrl), error = identity)
      if (inherits(ref, "error")) {
        expect_true(inherits(res, "error"), info = info)
        expect_identical(class(res), class(ref), info = info); expect_identical(conditionMessage(res), conditionMessage(ref), info = info)
        next
      }
      expect_identical(res$reduction_summary$reduction_stage, "outcome_free_compact_canonical_basis", info = info)
      for (f in cb_fields) cb_same(res[[f]], ref[[f]], paste(info, f))
      for (f in cb_cert_fields) cb_same(res$structural_certificate[[f]], ref$structural_certificate[[f]], paste(info, "cert", f))
      # Literal calendar values are carried exactly (no formatted key) on both routes.
      cls <- res$structural_certificate$classes[[1L]]
      expect_identical(cls$dates, sort(unique(c(res$metadata$catalog$rows$s, res$metadata$catalog$rows$t))), info = info)
      expect_identical(cls$dates, ref$structural_certificate$classes[[1L]]$dates, info = info)
      expect_identical(res$reduction_summary$n_link_views_constructed, length(res$structural_certificate$retained_index), info = info)
    }
  }
})
