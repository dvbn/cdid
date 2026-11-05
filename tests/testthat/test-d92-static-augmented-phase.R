# Existing D92 stores, partitions and scientific references; only the ordinary
# static augmented gate and its already-live diagnostic payload are new.
aps <- function(name) getFromNamespace(name, "cdid")
aps_helpers <- function(file, wanted) {
  e <- new.env(parent = environment()); e$api <- asNamespace("cdid")
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1), envir = e, inherits = FALSE)))
  e
}
aps_aug <- aps_helpers("test-augmented-outcome-pooling.R",
  c("rap_ns", "rap_gc", "rap_case", "rap_fit", "rap_point", "rap_reference"))
aps_phase <- aps_helpers("test-d92-phase-retained-if.R", c("d92_distinct_partition_ns"))
aps_public <- aps_helpers("test-d92-shared-if-reuse.R",
  c("d92_route_panel", "d92_run_public", "d92_strip", "d92_find_store", "d92_expect_parity"))

testthat::test_that("static augmented phase admission adds exactly the ordinary route", {
  eligible <- aps("cp04_retained_if_phase_static_augmented")
  # Stop at the unchanged domain boundary to inspect eligibility independently
  # of fitting. No fake backend is ever admitted to a store or IF producer.
  e <- new.env(parent = asNamespace("cdid"))
  reserve <- aps("cp04_retained_if_phase_reserve"); environment(reserve) <- e
  e$cp04_retained_if_validate_backend <- function(backend) invisible(TRUE)
  e$cp04_backend_domain <- function(backend) stop("original phase domain reached")
  old_eligible <- function(b) inherits(b, "cdid_ipw_raw_backend") ||
    aps("cp04_retained_if_phase_baseline")(b) ||
    (inherits(b, "cdid_augmented_backend") && !is.null(b$reference_ref) && !is.null(b$layout))
  for (cls in c("cdid_augmented_backend", "cdid_imputation_backend", "cdid_ipw_raw_backend",
      "cdid_or_baseline", "cdid_dr_baseline", "cdid_cp04_corollary1", "cdid_or_corollary1",
      "cdid_dr_corollary1", "cdid_reference_ipw_backend", "cdid_cp04_custom_backend")) {
    for (ref in c(FALSE, TRUE)) for (layout in c(FALSE, TRUE)) for (pair in c(FALSE, TRUE)) {
      b <- structure(list(reference_ref = if (ref) new.env(parent = emptyenv()) else NULL,
        layout = if (layout) list(frame_id = 1L) else NULL,
        pair_outcomes = if (pair) list(private = TRUE) else NULL, outcome_pooling = "separate"), class = cls)
      added <- cls == "cdid_augmented_backend" && !ref && !layout && !pair
      testthat::expect_identical(eligible(b), added)
      call <- function(transfer = FALSE) reserve(b, 1L, "covariance", 1L,
        list(blocks = list(1L)), 0, 0, 0, NULL, numerical_transfer = transfer)
      if (old_eligible(b) || added) {
        testthat::expect_error(call(), "original phase domain reached", fixed = TRUE)
      } else testthat::expect_null(call())
      testthat::expect_null(call(TRUE))
    }
  }
  testthat::expect_false(eligible(NULL))
  testthat::expect_null(reserve(NULL, 1L, "covariance", 1L, NULL, 0, 0, 0, NULL))
  testthat::expect_error(reserve(NULL, 1L, "bad", 1L, NULL, 0, 0, 0, NULL),
    class = "cdid_cp04_input_error")
})

testthat::test_that("static augmented phase stores preserve all F rows and exact capacity boundaries", {
  for (model in c("pair_logit", "independent_period_logit")) {
    a <- aps_aug$rap_case()
    fit <- aps_aug$rap_fit(a, imp = "across_cohorts", sel = "within_cohort", model = model, mode = "pooled")
    testthat::expect_true(aps("cp04_retained_if_phase_static_augmented")(fit))
    testthat::expect_equal(aps_aug$rap_point(fit), aps_aug$rap_reference(a, fit), tolerance = 2e-8)
    links <- c(3L, 1L, 2L); N <- nrow(a$Y); B <- fit$control$max_chunk_rows
    raw <- matrix(0, N, length(links))
    aps("cp04_replay_links")(fit, links, function(z) raw[z$rows, ] <<- z$phi, chunk_rows = B)
    reference <- aps("cp04_stream_covariance")(fit, links)
    caller <- unname(aps("cp04_payload")(list(raw, reference))["entries"])
    plan <- aps("cp04_stream_plan")(fit, links, caller)
    point <- aps("cp04_preflight_links")(fit, links, caller, chunk_rows = B)
    d5 <- list(ledger_if_entries = 17, basis_if_entries = 19,
      ledger_enclosure_entries = 0, basis_enclosure_entries = 0)
    inv <- aps("cp04_retained_if_phase_inventory")("covariance", fit, plan, point,
      caller, caller, d5, length(links), length(fit$links))
    testthat::expect_length(plan$blocks, 1L)
    testthat::expect_equal(plan$chunk_rows, B)
    testthat::expect_equal(inv$d5_persistent_entries, 36)
    boundary <- aps("cp04_retained_if_boundary_cost")(fit, links, B)
    entries <- aps("cp04_retained_if_entries")(N, length(links), B) + boundary$entries
    work <- aps("cp04_retained_if_work")(N, length(links)) + (32 + 8 * (ceiling(N / B) + 1)) * boundary$work
    cap <- inv$live + entries
    exact <- aps("ipw_owner_control")(B, cap, cap, fit$control$max_work_ops, fit$control$max_output_entries)
    reserve <- function(phase, live = inv$live, control = exact, active = NULL) {
      aps("cp04_retained_if_phase_reserve")(fit, links, phase, B, plan,
        retained = inv$retained, live = live, work = 0, control = control, active_store = active)
    }
    testthat::expect_null(reserve("covariance", live = inv$live + 1))
    testthat::expect_null(reserve("gmm", live = inv$live + 1))
    multi <- plan; multi$blocks <- list(1L, 2:3)
    testthat::expect_null(aps("cp04_retained_if_phase_reserve")(fit, links, "covariance", B,
      multi, 0, 0, 0, exact))
    cov_store <- reserve("covariance")
    testthat::expect_true(is.environment(cov_store))
    testthat::expect_identical(cov_store$entries, entries)
    testthat::expect_identical(cov_store$work, work)
    testthat::expect_equal(cov_store$admission_total, cap)
    testthat::expect_equal(cov_store$B, plan$chunk_rows)
    producer <- aps("cp04_retained_if_bind")(fit, cov_store)
    cached <- aps("cp04_stream_covariance")(producer, links, inv$retained)
    for (nm in c("estimates", "mean_phi", "covariance"))
      testthat::expect_equal(cached[[nm]], reference[[nm]], tolerance = 1e-12)
    testthat::expect_equal(unname(cov_store$phi), raw, tolerance = 1e-12)
    testthat::expect_equal(dim(cov_store$phi), c(N, length(links)))
    testthat::expect_identical(cov_store$producer_passes, 1L)
    testthat::expect_identical(cov_store$slice_passes, 1L)
    # The established raw-frequency reference includes zero/single-outcome F
    # units and refits the nuisances, independently of either replay route.
    reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
    for (i in c(which(rowSums(a$obs) == 0L)[1L], which(rowSums(a$obs) == 1L)[1L])) {
      direction <- rep(-1, N); direction[i] <- N - 1
      numerical <- (aps_aug$rap_reference(a, fit, 1 + 1e-6 * direction) -
        aps_aug$rap_reference(a, fit, 1 - 1e-6 * direction)) / 2e-6
      testthat::expect_equal(unname(cov_store$phi[i, ]) / reader$state$lambda,
        tail(numerical, length(fit$links))[links], tolerance = 1e-5)
    }
    testthat::expect_error(reserve("gmm", active = cov_store), class = "cdid_cp04_integrity_error")
    aps("cp04_retained_if_release")(cov_store); producer <- NULL
    testthat::expect_false(exists("phi", cov_store, inherits = FALSE))
    gmm_store <- reserve("gmm", active = cov_store)
    testthat::expect_true(is.environment(gmm_store))
    testthat::expect_identical(gmm_store$entries, entries)
    producer <- aps("cp04_retained_if_bind")(fit, gmm_store)
    for (pass in 1:2) aps("cp04_replay_links")(producer, links, function(z)
      testthat::expect_equal(unname(z$phi), raw[z$rows, , drop = FALSE], tolerance = 1e-12), chunk_rows = B)
    testthat::expect_null(aps("cp04_retained_if_auxiliary_admit")(gmm_store, fit, links,
      retained = 0, live = 0, work = 0, control = exact))
    testthat::expect_identical(gmm_store$producer_passes, 1L)
    testthat::expect_identical(gmm_store$slice_passes, 1L)
    aps("cp04_retained_if_release")(gmm_store); producer <- NULL
    bad <- fit; bad$links[[1L]]$mu <- bad$links[[1L]]$mu + 1
    testthat::expect_error(aps("cp04_retained_if_phase_reserve")(bad, links, "covariance", B,
      plan, 0, 0, 0, exact), class = "cdid_augmented_integrity_error")
    failed <- reserve("covariance")
    testthat::expect_error(aps("cp04_replay_links")(aps("cp04_retained_if_bind")(fit, failed),
      links, function(z) stop("phase callback failed"), chunk_rows = B), "phase callback failed", fixed = TRUE)
    testthat::expect_identical(failed$state, "released")
    testthat::expect_false(exists("phi", failed, inherits = FALSE))
  }
})

aps_pipeline_ns <- function(gmm = TRUE, stream = FALSE, fail_at = NULL, allow_transfer = TRUE) {
  api <- asNamespace("cdid")
  e <- aps_phase$d92_distinct_partition_ns(B_cov = 7L, B_gmm = 5L)
  trace <- new.env(parent = emptyenv()); trace$stores <- list(); trace$inventory <- list(); trace$rows <- list()
  trace$transition_attempted <- FALSE; trace$transferred <- FALSE
  trace$max_live_stores <- 0L; trace$streamed_passes <- 0L; trace$partitions <- integer()
  pipeline <- e$run_chained_ipw_pipeline
  e$run_chained_ipw_pipeline <- function(...) {
    out <- pipeline(...); trace$result <- out; out
  }
  replay <- aps("cp04_retained_if_replay"); environment(replay) <- e
  e$cp04_retained_if_replay <- function(store, backend, links, callback, scale, B, pre) {
    trace$partitions <- c(trace$partitions, B)
    replay(store, backend, links, callback, scale, B, pre)
  }
  replay_links <- aps("cp04_replay_links"); environment(replay_links) <- e
  e$cp04_replay_links <- function(backend, links, callback, scale = "core", chunk_rows = NULL) {
    if (!inherits(backend, "cdid_cp04_retained_if_bound"))
      trace$streamed_passes <- trace$streamed_passes + 1L
    replay_links(backend, links, callback, scale, chunk_rows)
  }
  # Force only the whole-lifetime decline. Every phase's original charged
  # inventory, planner probe and one-store-at-a-time guard remain live.
  e$cp04_retained_if_reserve <- function(...) NULL
  inventory <- aps("cp04_retained_if_phase_inventory")
  e$cp04_retained_if_phase_inventory <- function(phase, backend, cov_plan, point_preflight,
      retained_with_diagnostics, cp04_planner_retained, d5_arrays, L_used, M_all) {
    f <- parent.frame()
    actual <- unname(aps("cp04_payload")(get("cp04_diagnostic_selection", f))["entries"])
    original <- get("cp04_retained_with_diagnostics", f)
    if (phase == "covariance") {
      testthat::expect_equal(retained_with_diagnostics - original, actual)
      testthat::expect_gt(actual, 0)
      testthat::expect_gt(nrow(get("cp04_diagnostic_selection", f)$active_cells), 0)
    } else {
      testthat::expect_identical(retained_with_diagnostics, original)
      testthat::expect_identical(cp04_planner_retained, get("cp04_planner_retained", f))
      testthat::expect_gte(get("final_diagnostic_copy_entries", f), 5 * actual)
      if (!trace$transition_attempted) {
        # The transition inventories the filled covariance buffer, before
        # deciding whether the same buffer can serve the GMM partition.
        testthat::expect_identical(trace$stores$covariance$state, "filled")
        testthat::expect_true(exists("phi", trace$stores$covariance, inherits = FALSE))
        testthat::expect_identical(trace$stores$covariance$producer_passes, 1L)
        testthat::expect_identical(trace$stores$covariance$slice_passes, 1L)
      } else {
        testthat::expect_false(trace$transferred)
        testthat::expect_identical(trace$stores$covariance$state, "released")
        testthat::expect_false(exists("phi", trace$stores$covariance, inherits = FALSE))
      }
    }
    out <- inventory(phase, backend, cov_plan, point_preflight, retained_with_diagnostics,
      cp04_planner_retained, d5_arrays, L_used, M_all)
    before <- inventory(phase, backend, cov_plan, point_preflight, original,
      cp04_planner_retained, d5_arrays, L_used, M_all)
    for (nm in c("retained", "live")) testthat::expect_equal(out[[nm]] - before[[nm]],
      if (phase == "covariance") actual else 0)
    trace$inventory[[phase]] <- list(out = out, before = before, diagnostic_entries = actual,
      covariance_partition = cov_plan$chunk_rows, gmm_partition = point_preflight$chunk_rows)
    out
  }
  reserve <- e$cp04_retained_if_phase_reserve
  e$cp04_retained_if_phase_reserve <- function(backend, links, phase, ...) {
    if (stream || (phase == "gmm" && !gmm)) return(NULL)
    if (phase == "gmm") {
      testthat::expect_identical(trace$stores$covariance$state, "released")
      testthat::expect_false(exists("phi", trace$stores$covariance, inherits = FALSE))
    }
    out <- reserve(backend, links, phase = phase, ...)
    testthat::expect_true(is.environment(out))
    trace$stores[[phase]] <- out
    live <- sum(vapply(trace$stores, function(s) !identical(s$state, "released"), logical(1L)))
    trace$max_live_stores <- max(trace$max_live_stores, live)
    testthat::expect_identical(live, 1L)
    out
  }
  transfer <- aps("cp04_retained_if_covariance_to_gmm")
  e$cp04_retained_if_covariance_to_gmm <- function(store, backend, links, chunk_rows, retained, work, control) {
    trace$transition_attempted <- TRUE
    testthat::expect_identical(store, trace$stores$covariance)
    testthat::expect_identical(store$state, "filled")
    testthat::expect_identical(store$producer_passes, 1L)
    testthat::expect_identical(store$slice_passes, 1L)
    if (!gmm || !allow_transfer) return(NULL)
    original <- store$phi; raw_ref <- store$raw_ref
    out <- transfer(store, backend, links, chunk_rows, retained, work, control)
    testthat::expect_identical(out, store)
    testthat::expect_true(api$gmm_sexp_aliases(out$phi, original))
    testthat::expect_identical(out$raw_ref, raw_ref)
    testthat::expect_true(api$gmm_sexp_aliases(out$phi, out$raw_ref$phi))
    testthat::expect_identical(out$producer_passes, 1L)
    testthat::expect_identical(out$slice_passes, 1L)
    testthat::expect_lte(out$gmm_transfer$admission_total,
      min(control$max_retained_entries, control$max_matrix_entries))
    trace$transferred <- TRUE
    out
  }
  e$chained_full_if_chunk <- function(backend, rows, links, scale, plan = NULL) {
    trace$rows[[length(trace$rows) + 1L]] <- rows
    api$chained_full_if_chunk(backend, rows, links, scale, plan)
  }
  e$gmm_covariance_provenance <- function(...) {
    if (!stream) {
      store <- trace$stores$covariance
      if (trace$transferred) {
        testthat::expect_true(trace$transferred)
        testthat::expect_identical(store$state, "filled")
        testthat::expect_true(api$gmm_sexp_aliases(store$phi, store$raw_ref$phi))
        testthat::expect_identical(store$producer_passes, 1L)
      } else testthat::expect_identical(store$state, "released")
    }
    api$gmm_covariance_provenance(...)
  }
  e$build_ipw_auxiliary_reporting <- function(...) {
    for (store in trace$stores) {
      testthat::expect_identical(store$state, "released")
      for (nm in c("phi", "raw_ref", "backend_snapshot", "reader_snapshot", "gmm_transfer"))
        testthat::expect_false(exists(nm, store, inherits = FALSE))
    }
    api$build_ipw_auxiliary_reporting(...)
  }
  if (!is.null(fail_at)) assign(fail_at, function(...) stop("pipeline phase failed"), e)
  e$trace <- trace
  e
}

testthat::test_that("static augmented pipeline funds its actual diagnostic and releases each phase", {
  d <- aps_public$d92_route_panel()
  formula <- ~x
  run <- function(e) aps_public$d92_run_public(d, chunk = 11L, ns = e,
    xformla = formula, weightsname = "w", est_method = "dr", selection = "conditional_mtar",
    selection_est_method = "dr", selection_model = "pair_logit", selection_pooling = "within_cohort",
    imputation_pooling = "across_cohorts")
  retained <- aps_pipeline_ns(); streamed <- aps_pipeline_ns(stream = TRUE)
  # Streamed inventory still runs; it has no covariance store to inspect.
  streamed$cp04_retained_if_phase_inventory <- aps("cp04_retained_if_phase_inventory")
  out <- run(retained); ref <- run(streamed)
  a <- retained$trace$result; b <- streamed$trace$result
  aps_public$d92_expect_parity(a, b)
  testthat::expect_equal(out$att, ref$att, tolerance = 1e-12)
  testthat::expect_equal(out$inffunc, ref$inffunc, tolerance = 1e-12)
  testthat::expect_identical(a$selection, b$selection)
  testthat::expect_identical(a$support, b$support)
  testthat::expect_identical(a$reporting$target_meta, b$reporting$target_meta)
  record <- a$resource_accounting$point_producer$retained_moment_if
  testthat::expect_identical(record$mode, "covariance_to_gmm")
  testthat::expect_equal(c(covariance = record$covariance_chunk_rows, gmm = record$gmm_chunk_rows), c(covariance = 7, gmm = 5))
  testthat::expect_identical(record$producer_passes, 1L)
  testthat::expect_identical(record$slice_passes, 3L)
  testthat::expect_identical(retained$trace$streamed_passes, 0L)
  testthat::expect_identical(retained$trace$max_live_stores, 1L)
  testthat::expect_identical(names(retained$trace$stores), "covariance")
  testthat::expect_equal(retained$trace$partitions, c(7, 7, 5, 5))
  testthat::expect_equal(record$n_links, length(a$point_links_used))
  testthat::expect_identical(a$reduction_summary, b$reduction_summary)
  N <- a$nuisance_registry$cp04_reader$state$N
  testthat::expect_length(retained$trace$rows, ceiling(N / 7))
  testthat::expect_identical(unlist(retained$trace$rows, use.names = FALSE), seq_len(N))
  for (store in retained$trace$stores) {
    testthat::expect_identical(store$state, "released")
    for (nm in c("phi", "raw_ref", "backend_snapshot", "reader_snapshot", "gmm_transfer"))
      testthat::expect_false(exists(nm, store, inherits = FALSE))
  }
  testthat::expect_false(aps_public$d92_find_store(a))
  testthat::expect_false(aps_public$d92_find_store(out))
  # Declining only transfer retains the original two independent stores.
  # The covariance allocation is released before the GMM store is admitted.
  separate <- aps_pipeline_ns(allow_transfer = FALSE)
  separate_out <- run(separate)
  separate_result <- separate$trace$result
  aps_public$d92_expect_parity(separate_result, b)
  testthat::expect_equal(separate_out$att, ref$att, tolerance = 1e-12)
  testthat::expect_equal(separate_out$inffunc, ref$inffunc, tolerance = 1e-12)
  separate_record <- separate_result$resource_accounting$point_producer$retained_moment_if
  testthat::expect_identical(separate_record$mode, "phase_local")
  testthat::expect_equal(separate_record$chunk_rows, c(covariance = 7, gmm = 5))
  testthat::expect_identical(separate_record$producer_passes, 2L)
  testthat::expect_identical(separate_record$slice_passes, 2L)
  testthat::expect_identical(separate_record$streamed_passes, 0L)
  testthat::expect_identical(separate_record$simultaneous_stores, 1L)
  testthat::expect_identical(separate$trace$max_live_stores, 1L)
  testthat::expect_false(separate$trace$transferred)
  testthat::expect_identical(names(separate$trace$stores), c("covariance", "gmm"))
  testthat::expect_false(identical(separate$trace$stores$covariance, separate$trace$stores$gmm))
  testthat::expect_equal(separate$trace$partitions, c(7, 7, 5, 5))
  testthat::expect_length(separate$trace$rows, ceiling(N / 7) + ceiling(N / 5))
  testthat::expect_identical(unlist(separate$trace$rows, use.names = FALSE), c(seq_len(N), seq_len(N)))
  for (store in separate$trace$stores) {
    testthat::expect_identical(store$state, "released")
    for (nm in c("phi", "raw_ref", "backend_snapshot", "reader_snapshot", "gmm_transfer"))
      testthat::expect_false(exists(nm, store, inherits = FALSE))
  }
  testthat::expect_false(aps_public$d92_find_store(separate_result))
  testthat::expect_false(aps_public$d92_find_store(separate_out))
  # A legitimate GMM decline keeps the covariance-only gain and streams the
  # two original GMM passes at their original partition.
  cov_only <- aps_pipeline_ns(gmm = FALSE)
  run(cov_only)
  aps_public$d92_expect_parity(cov_only$trace$result, b)
  cov_record <- cov_only$trace$result$resource_accounting$point_producer$retained_moment_if
  testthat::expect_identical(cov_record$covariance$status, "retained")
  testthat::expect_identical(cov_record$gmm$status, "streamed")
  testthat::expect_identical(cov_record$streamed_passes, 2L)
  for (failure in c("cp04_stream_covariance", "solve_chained_gmm_core")) {
    failed <- aps_pipeline_ns(fail_at = failure)
    testthat::expect_error(run(failed), "pipeline phase failed", fixed = TRUE)
    testthat::expect_true(length(failed$trace$stores) > 0L)
    for (store in failed$trace$stores) {
      testthat::expect_identical(store$state, "released")
      testthat::expect_false(exists("phi", store, inherits = FALSE))
    }
  }
})
