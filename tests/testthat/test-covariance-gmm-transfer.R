v61_ns <- function(nm) getFromNamespace(nm, "cdid")
v61_e <- new.env(parent = asNamespace("cdid"))
for (x in parse(file = testthat::test_path("test-augmented-adjoint-replay.R")))
  if (is.call(x) && identical(x[[1L]], as.name("<-"))) eval(x, v61_e)

v61_filled <- function() {
  a <- v61_e$v57_aug$rap_case()
  fit <- v61_e$v57_aug$rap_fit(a, imp = "separate", model = "pair_logit", mode = "pooled")
  links <- rev(seq_along(fit$links)); N <- nrow(a$Y)
  raw <- matrix(0, N, length(links))
  v61_ns("cp04_replay_links")(fit, links, function(b) raw[b$rows, ] <<- b$phi,
    scale = "core", chunk_rows = 7L)
  reference <- v61_ns("cp04_stream_covariance")(fit, links)
  caller <- unname(v61_ns("cp04_payload")(list(raw, reference))["entries"])
  covariance <- v61_ns("cp04_stream_plan")(fit, links, caller)
  point <- v61_ns("cp04_preflight_links")(fit, links, caller, chunk_rows = 7L)
  d5 <- list(ledger_if_entries = 0, basis_if_entries = 0,
    ledger_enclosure_entries = 0, basis_enclosure_entries = 0)
  inv <- v61_ns("cp04_retained_if_phase_inventory")("covariance", fit, covariance, point,
    caller, caller, d5, length(links), length(fit$links))
  store <- v61_ns("cp04_retained_if_phase_reserve")(fit, links, "covariance", covariance$chunk_rows,
    covariance, retained = inv$retained, live = inv$live, work = 0, control = fit$control)
  testthat::expect_true(is.environment(store))
  actual <- v61_ns("cp04_stream_covariance")(v61_ns("cp04_retained_if_bind")(fit, store), links, inv$retained)
  testthat::expect_equal(actual$covariance, reference$covariance, tolerance = 1e-12)
  testthat::expect_identical(store$producer_passes, 1L)
  testthat::expect_identical(store$slice_passes, 1L)
  list(a = a, fit = fit, links = links, raw = raw, reference = reference, covariance = covariance,
    point = point, d5 = d5, caller = caller, store = store)
}

testthat::test_that("one covariance matrix serves original GMM chunks with no second producer", {
  z <- v61_filled(); store <- z$store; fit <- z$fit; links <- z$links
  on.exit(v61_ns("cp04_retained_if_release")(store), add = TRUE)
  inv <- v61_ns("cp04_retained_if_phase_inventory")("gmm", fit, z$covariance, z$point,
    z$caller, z$caller, z$d5, length(links), length(fit$links))
  B <- z$point$chunk_rows
  testthat::expect_false(isTRUE(B == store$B))
  original <- store$phi; raw_ref <- store$raw_ref
  peak <- max(B, store$B)
  boundary <- v61_ns("cp04_retained_if_boundary_cost")(fit, links, peak)
  entries <- 1024 + max(store$entries,
    v61_ns("cp04_retained_if_entries")(store$N, store$L, peak) + boundary$entries)
  total <- inv$retained + fit$accounting$owned_entries + fit$accounting$borrowed_entries + entries
  names_before <- ls(store, all.names = TRUE)
  tight <- v61_ns("ipw_owner_control")(fit$control$max_chunk_rows, total - 1, total - 1,
    fit$control$max_work_ops, fit$control$max_output_entries)
  testthat::expect_null(v61_ns("cp04_retained_if_covariance_to_gmm")(store, fit, links, B,
    retained = inv$retained, work = 0, control = tight))
  testthat::expect_identical(ls(store, all.names = TRUE), names_before)
  testthat::expect_true(v61_ns("gmm_sexp_aliases")(store$phi, original))
  exact <- v61_ns("ipw_owner_control")(fit$control$max_chunk_rows, total, total,
    fit$control$max_work_ops, fit$control$max_output_entries)
  transferred <- v61_ns("cp04_retained_if_covariance_to_gmm")(store, fit, links, B,
    retained = inv$retained, work = 0, control = exact)
  testthat::expect_identical(transferred, store)
  testthat::expect_true(v61_ns("gmm_sexp_aliases")(store$phi, original))
  testthat::expect_identical(store$raw_ref, raw_ref)
  testthat::expect_equal(store$gmm_transfer$admission_total, total, tolerance = 0)
  testthat::expect_lte(store$gmm_transfer$admission_total, exact$max_matrix_entries)
  testthat::expect_error(v61_ns("cp04_retained_if_check")(store, fit, rev(links), B, "core"),
    class = "cdid_cp04_integrity_error")
  testthat::expect_error(v61_ns("cp04_retained_if_check")(store, fit, links, B - 1, "core"),
    class = "cdid_cp04_integrity_error")
  producer <- v61_ns("cp04_retained_if_bind")(fit, store)
  for (pass in 1:2) {
    actual <- matrix(0, nrow(z$raw), ncol(z$raw)); visited <- integer(0L)
    v61_ns("cp04_replay_links")(producer, links, function(b) {
      testthat::expect_identical(b$rows, seq.int(min(b$rows), max(b$rows)))
      testthat::expect_lte(length(b$rows), B)
      visited <<- c(visited, b$rows)
      actual[b$rows, ] <<- b$phi
    }, scale = "core", chunk_rows = B)
    testthat::expect_identical(visited, seq_len(nrow(z$raw)))
    testthat::expect_equal(actual, z$raw, tolerance = 1e-12)
    testthat::expect_identical(actual == 0, z$raw == 0)
    testthat::expect_equal(crossprod(actual), crossprod(z$raw), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(actual))), sqrt(diag(crossprod(z$raw))), tolerance = 1e-12)
  }
  testthat::expect_identical(store$producer_passes, 1L)
  testthat::expect_identical(store$slice_passes, 3L)
  testthat::expect_true(v61_ns("gmm_sexp_aliases")(store$phi, original))
  record <- v61_ns("cp04_retained_if_release")(store)
  testthat::expect_identical(record$mode, "covariance_to_gmm")
  testthat::expect_identical(record$producer_passes, 1L)
  testthat::expect_identical(record$slice_passes, 3L)
  testthat::expect_equal(record$covariance_chunk_rows, z$covariance$chunk_rows, tolerance = 0)
  testthat::expect_equal(record$gmm_chunk_rows, B, tolerance = 0)
  testthat::expect_identical(store$state, "released")
  testthat::expect_false(any(c("phi", "raw_ref", "backend_snapshot", "gmm_transfer") %in% ls(store)))
})

testthat::test_that("transferred covariance storage releases on an actual GMM slice callback error", {
  z <- v61_filled(); store <- z$store
  on.exit(v61_ns("cp04_retained_if_release")(store), add = TRUE)
  inv <- v61_ns("cp04_retained_if_phase_inventory")("gmm", z$fit, z$covariance, z$point,
    z$caller, z$caller, z$d5, length(z$links), length(z$fit$links))
  transferred <- v61_ns("cp04_retained_if_covariance_to_gmm")(store, z$fit, z$links, z$point$chunk_rows,
    retained = inv$retained, work = 0, control = z$fit$control)
  testthat::expect_identical(transferred, store)
  testthat::expect_error(v61_ns("cp04_replay_links")(v61_ns("cp04_retained_if_bind")(z$fit, store),
    z$links, function(b) stop("actual GMM callback failure"), chunk_rows = z$point$chunk_rows),
    "actual GMM callback failure", fixed = TRUE)
  testthat::expect_identical(store$state, "released")
  testthat::expect_false(any(c("phi", "raw_ref", "backend_snapshot", "gmm_transfer") %in% ls(store)))
})


testthat::test_that("historical replay peaks do not become final diagnostic live reservations", {
  f <- v61_ns("run_chained_ipw_pipeline")
  matches <- list()
  walk <- function(x) {
    if (!is.call(x)) return(invisible(NULL))
    if (identical(x[[1L]], as.name("<-")) &&
        identical(x[[2L]], as.name("replay_peak_entries"))) matches[[length(matches) + 1L]] <<- x
    for (i in seq_along(x)[-1L]) if (!identical(x[[i]], quote(expr = ))) walk(x[[i]])
    invisible(NULL)
  }
  walk(body(f))
  initial <- Filter(function(x) identical(x[[3L]], as.name("combined_live")), matches)
  transition <- Filter(function(x) "transferred" %in% all.names(x[[3L]]), matches)
  testthat::expect_length(initial, 1L); testthat::expect_length(transition, 1L)
  e <- new.env(parent = baseenv()); e$combined_live <- 100
  eval(initial[[1L]], e)
  e$transferred <- list(gmm_transfer = list(admission_total = 150))
  eval(transition[[1L]], e)
  testthat::expect_equal(e$combined_live, 100, tolerance = 0)
  testthat::expect_equal(e$replay_peak_entries, 150, tolerance = 0)
  # Keep the actual diagnostic current-prefix gate; a released replay matrix
  # is represented only by the historical scalar, not retained payload.
  res <- list(resource_accounting = list(simultaneous_entries = e$combined_live,
    initial_envelope_entries = 90, historical_replay_peak_entries = e$replay_peak_entries))
  retained <- v61_ns("cp04_public_diagnostic_retained")(res)
  testthat::expect_equal(retained, 100, tolerance = 0)
  control <- v61_ns("ipw_owner_control")(1L, 140, 140, 1000, 100)
  testthat::expect_true(v61_ns("cp04_budget")(control, c(retained, 5), 20, 0, 5,
    "released replay diagnostic current-phase witness"))
  testthat::expect_error(v61_ns("cp04_budget")(control, c(e$replay_peak_entries, 5), 20, 0, 5,
    "incorrect historical-as-live witness"), class = "cdid_ipw_owner_bound_error")
})
