# Owner storage counted once: the persistent prefix absorbs the fitted backend's
# owned entries (and, on CP04 routes, the reader's), and already holds the frame
# and panel owner storage; compositions that add the backend's owned + borrowed
# entries remove that storage once, as recorded in a sealed record.
osc_ns <- function(name) get(name, asNamespace("cdid"), inherits = FALSE)

osc_data <- function(N = 60L, TT = 8L) {
  set.seed(11)
  g <- sample(c(0, 4, 6), N, TRUE)
  d <- expand.grid(t = seq_len(TT), id = seq_len(N)); d$g <- g[d$id]
  d$y <- stats::ave(stats::rnorm(nrow(d)), d$id, FUN = cumsum) + ifelse(d$g > 0 & d$t >= d$g, 1, 0)
  d$y[stats::runif(nrow(d)) < 0.2] <- NA
  d
}

# Runs one fit and captures the pipeline's registry and its prefix at the first
# stream plan (after every composition below has been formed).
osc_fit <- function(est_method, weighting = "within", extra = list()) {
  seen <- new.env(parent = emptyenv())
  plan <- osc_ns("cp04_stream_plan")
  testthat::local_mocked_bindings(
    cp04_stream_plan = function(backend, links, retained = 0) {
      pf <- parent.frame()
      if (is.null(seen$registry) && exists("nuis_res", envir = pf, inherits = FALSE)) {
        seen$registry <- get("nuis_res", pf)
        for (nm in c("combined_initial", "owners_live", "or_preprocess_entries", "cp04_retained_with_diagnostics",
                     "retained_with_diagnostics", "cp04_planner_retained", "frame_owner", "panel_owner",
                     "final_diagnostic_copy_entries", "gmm_preflight", "d5_arrays"))
          if (exists(nm, envir = pf, inherits = FALSE)) assign(nm, get(nm, pf), envir = seen)
      }
      plan(backend, links, retained)
    }, .package = "cdid")
  args <- c(list(yname = "y", tname = "t", idname = "id", gname = "g", data = osc_data(), xformla = ~1,
    est_method = est_method, gmm_weighting = weighting, control_group = "nevertreated", base_period = "universal",
    allow_unbalanced_panel = TRUE, bstrap = FALSE, cband = FALSE), extra)
  fit <- suppressWarnings(do.call(cdid::att_gt, args))
  list(fit = fit, seen = seen)
}

testthat::test_that("the CP04 record holds backend, reader and owner storage once each", {
  x <- osc_fit("dr")
  reg <- x$seen$registry; b <- reg$cp04; reader <- reg$cp04_reader
  testthat::expect_true(inherits(b, "cdid_dr_baseline"))
  rec <- osc_ns("cp04_owner_overlap_checked")(reg)
  testthat::expect_false(is.null(rec))
  owners <- osc_ns("ipw_point_view_owner_entries")(x$seen$frame_owner, x$seen$panel_owner)
  testthat::expect_identical(rec$backend_owned_entries, b$accounting$owned_entries)
  testthat::expect_identical(rec$reader_owned_entries, reader$accounting$owned_entries)
  testthat::expect_identical(rec$owner_entries, reader$accounting$borrowed_owner_entries)
  testthat::expect_equal(rec$owner_entries, owners)
  testthat::expect_identical(rec$overlap_entries,
    b$accounting$owned_entries + reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries)
  # The owner storage is part of the prefix's owner term; backend and reader were added to it.
  testthat::expect_lte(rec$owner_entries, x$seen$owners_live - x$seen$or_preprocess_entries)
  testthat::expect_lte(rec$prefix_entries, x$seen$combined_initial)
  # Each composition removes exactly the external prefix and the recorded storage.
  r <- x$seen$retained_with_diagnostics
  testthat::expect_identical(x$seen$cp04_retained_with_diagnostics,
    osc_ns("cp04_registry_retained_entries")(reg, r) - rec$overlap_entries)
  d5 <- x$seen$d5_arrays
  planner <- x$seen$final_diagnostic_copy_entries + x$seen$combined_initial + x$seen$gmm_preflight$total_payload_entries +
    d5$ledger_if_entries + d5$basis_if_entries + d5$ledger_enclosure_entries + d5$basis_enclosure_entries +
    d5$retained_diagnostic_entries
  testthat::expect_identical(x$seen$cp04_planner_retained,
    osc_ns("cp04_registry_retained_entries")(reg, planner) - rec$overlap_entries)
  pc <- x$fit$chained$resource_accounting$memory_prefix_composition
  testthat::expect_identical(pc$simultaneous_entries_before_prefix_dedup -
    x$fit$chained$resource_accounting$simultaneous_entries, pc$embedded_external_entries + rec$overlap_entries)
})

testthat::test_that("the compact raw IPW record holds backend and owner storage once each", {
  x <- osc_fit("ipw", "optimal")
  reg <- x$seen$registry; b <- reg$raw_ipw
  testthat::expect_true(inherits(b, "cdid_ipw_raw_backend"))
  rec <- osc_ns("cp04_owner_overlap_checked")(reg)
  owners <- osc_ns("ipw_point_view_owner_entries")(x$seen$frame_owner, x$seen$panel_owner)
  testthat::expect_identical(rec$backend_borrowed_entries, owners)
  testthat::expect_identical(rec$overlap_entries, b$accounting$owned_entries + owners)
  testthat::expect_identical(rec$reader_owned_entries, 0)
  testthat::expect_identical(x$seen$cp04_retained_with_diagnostics, x$seen$retained_with_diagnostics - rec$overlap_entries)
  ra <- x$fit$chained$resource_accounting
  testthat::expect_identical(ra$memory_prefix_composition$simultaneous_entries_before_prefix_dedup - ra$simultaneous_entries,
    ra$memory_prefix_composition$embedded_external_entries + rec$overlap_entries)
})

testthat::test_that("storage outside a verified identity stays charged", {
  x <- osc_fit("dr")
  reg <- x$seen$registry; f <- x$seen$frame_owner; p <- x$seen$panel_owner
  rec_fn <- function(...) get("record", envir = osc_ns("cp04_owner_overlap_record")(...), inherits = FALSE)
  owners_in_prefix <- x$seen$owners_live - x$seen$or_preprocess_entries
  base <- rec_fn(reg, x$seen$combined_initial, f, p, owners_in_prefix)
  # Borrowed entries that are not exactly reader + owners + design + external:
  # only the backend's own storage is removed.
  bad <- reg; bad$cp04$accounting$borrowed_entries <- bad$cp04$accounting$borrowed_entries + 1
  r1 <- rec_fn(bad, x$seen$combined_initial, f, p, owners_in_prefix)
  testthat::expect_identical(r1$overlap_entries, reg$cp04$accounting$owned_entries)
  testthat::expect_identical(r1$reader_owned_entries, 0)
  # Owners the prefix does not cover stay charged.
  r2 <- rec_fn(reg, x$seen$combined_initial, f, p, base$owner_entries - 1)
  testthat::expect_identical(r2$owner_entries, 0)
  testthat::expect_identical(r2$overlap_entries, base$overlap_entries - base$owner_entries)
  # Owners of another frame stay charged.
  f2 <- f; f2$ref <- new.env()
  testthat::expect_identical(rec_fn(reg, x$seen$combined_initial, f2, p, owners_in_prefix)$owner_entries, 0)
})

testthat::test_that("the record is checked against its backend and its prefix", {
  x <- osc_fit("dr")
  reg <- x$seen$registry; composed <- osc_ns("cp04_composed_retained_entries")
  rec <- osc_ns("cp04_owner_overlap_checked")(reg)
  testthat::expect_error(composed(reg, rec$prefix_entries - 1), class = "cdid_cp04_integrity_error")
  testthat::expect_identical(composed(reg, rec$prefix_entries + 10),
    osc_ns("cp04_registry_retained_entries")(reg, rec$prefix_entries + 10) - rec$overlap_entries)
  moved <- reg; moved$cp04$accounting$owned_entries <- moved$cp04$accounting$owned_entries + 1
  testthat::expect_error(composed(moved, rec$prefix_entries), class = "cdid_cp04_integrity_error")
  copy <- rec; copy$overlap_entries <- 0
  forged <- reg; forged$prefix_owner_storage <- new.env(parent = emptyenv())
  assign("record", copy, envir = forged$prefix_owner_storage)
  testthat::expect_error(composed(forged, rec$prefix_entries), class = "cdid_cp04_integrity_error")
  # The diagnostic copy walks the returned result without the record, so the
  # record adds no counted entries of its own.
  plain2 <- reg; plain2$prefix_owner_storage <- NULL
  pay <- osc_ns("cp04_payload"); strip <- osc_ns("cp04_without_owner_record")
  testthat::expect_identical(unname(pay(strip(list(a = 1, nuisance_registry = reg)))["entries"]),
    unname(pay(list(a = 1, nuisance_registry = plain2))["entries"]))
  # Without a record nothing beyond the external prefix is removed.
  plain <- reg; plain$prefix_owner_storage <- NULL
  testthat::expect_identical(composed(plain, rec$prefix_entries), osc_ns("cp04_registry_retained_entries")(plain, rec$prefix_entries))
})

testthat::test_that("the post-pipeline diagnostic copy removes the owner storage, and only it, once", {
  seen <- new.env(parent = emptyenv())
  retained_fn <- osc_ns("cp04_public_diagnostic_retained")
  testthat::local_mocked_bindings(
    cp04_public_diagnostic_retained = function(res, public_object = NULL) {
      out <- retained_fn(res, public_object)
      if (is.null(seen$res)) { seen$res <- res; seen$public <- public_object; seen$out <- out }
      out
    }, .package = "cdid")
  x <- osc_fit("dr")
  res <- seen$res
  testthat::expect_false(is.null(res))
  phase <- res$resource_accounting$memory_phases
  testthat::expect_false(is.null(phase))
  rec <- osc_ns("cp04_owner_overlap_checked")(res$nuisance_registry)
  testthat::expect_gt(rec$owner_entries, 0)
  # The owner part differs from the whole recorded storage on this fit.
  testthat::expect_lt(rec$owner_entries, rec$overlap_entries)
  walk <- unname(osc_ns("cp04_payload")(list(osc_ns("cp04_without_owner_record")(res), seen$public))["entries"])
  # Before the change: persistent prefix plus the walk; after: less the owner storage exactly.
  testthat::expect_identical(seen$out, phase$persistent_entries - rec$owner_entries + walk)
  testthat::expect_identical(phase$persistent_entries + walk - seen$out, rec$owner_entries)
})

testthat::test_that("a record whose backend owned figure differs is refused", {
  x <- osc_fit("ipw", "optimal")
  reg <- x$seen$registry
  rec <- osc_ns("cp04_owner_overlap_checked")(reg)
  moved <- reg; moved$raw_ipw$accounting$owned_entries <- moved$raw_ipw$accounting$owned_entries + 1
  testthat::expect_error(osc_ns("cp04_owner_overlap_checked")(moved), class = "cdid_cp04_integrity_error")
  testthat::expect_error(osc_ns("cp04_composed_retained_entries")(moved, rec$prefix_entries), class = "cdid_cp04_integrity_error")
  testthat::expect_error(osc_ns("cp04_owner_overlap_entries")(moved), class = "cdid_cp04_integrity_error")
})
