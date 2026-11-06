# Whole-store lifetime accounting only; existing numerical tolerances stay fixed.
alt_ns <- function(name) getFromNamespace(name, "cdid")
alt_backend <- function(cls, reference = NULL, layout = NULL, complete = NULL, general = FALSE)
  structure(list(reference_ref = reference, layout = layout, complete_support = complete,
    uses_selection_models = general), class = cls)

testthat::test_that("the shared auxiliary gate preserves the existing backend predicate", {
  eligible <- alt_ns("cp04_retained_if_auxiliary_eligible")
  old <- function(b) is.list(b) && inherits(b, c("cdid_or_baseline", "cdid_dr_baseline", "cdid_ipw_raw_backend")) &&
    is.null(b$reference_ref) && !is.null(b$complete_support)
  classes <- c("cdid_or_baseline", "cdid_dr_baseline", "cdid_ipw_raw_backend",
    "cdid_cp04_corollary1", "cdid_or_corollary1", "cdid_dr_corollary1",
    "cdid_cp04_custom_backend", "cdid_reference_ipw_backend", "cdid_imputation_backend", "cdid_augmented_backend")
  for (cls in classes) for (complete in c(FALSE, TRUE)) for (reference in c(FALSE, TRUE)) {
    b <- alt_backend(cls, reference = if (reference) new.env() else NULL,
      complete = if (complete) list(dates = 1:4) else NULL)
    testthat::expect_identical(eligible(b), old(b))
  }
  testthat::expect_false(eligible(NULL))
  # Actual admission shares the same predicate, before any numeric planner.
  b <- alt_backend("cdid_dr_baseline", complete = list(dates = 1:4)); b$ref <- new.env()
  store <- new.env(parent = emptyenv()); store$kind <- "cp04_retained_if"
  store$state <- "reserved"; store$backend_ref <- b$ref
  testthat::local_mocked_bindings(cp04_retained_if_auxiliary_eligible = function(backend) FALSE, .package = "cdid")
  testthat::expect_null(alt_ns("cp04_retained_if_auxiliary_admit")(
    store, b, 1L, retained = 0, live = 0, work = 0, control = NULL))
})

testthat::test_that("only unextended non-transfer static conditional stores omit future output", {
  charge <- alt_ns("cp04_retained_if_auxiliary_output_entries")
  N <- 10000; Q <- 37L; R <- 20L; extra <- 8 * N * (Q + R)
  for (cls in c("cdid_cp04_corollary1", "cdid_or_corollary1", "cdid_dr_corollary1")) {
    for (general in c(FALSE, TRUE)) {
      b <- alt_backend(cls, general = general)
      testthat::expect_identical(charge(b, FALSE, N, Q, R), 0)
      testthat::expect_identical(charge(b, TRUE, N, Q, R), extra)
      b$reference_ref <- new.env()
      testthat::expect_identical(charge(b, FALSE, N, Q, R), extra)
      b$reference_ref <- NULL; b$layout <- list(frame_id = 1:3)
      testthat::expect_identical(charge(b, FALSE, N, Q, R), extra)
    }
  }
  for (cls in c("cdid_or_baseline", "cdid_dr_baseline", "cdid_ipw_raw_backend",
      "cdid_cp04_custom_backend", "cdid_reference_ipw_backend", "cdid_imputation_backend", "cdid_augmented_backend"))
    for (complete in c(FALSE, TRUE)) testthat::expect_identical(
      charge(alt_backend(cls, complete = if (complete) list(dates = 1:4) else NULL), FALSE, N, Q, R), extra)
  testthat::expect_identical(charge(NULL, FALSE, N, Q, R), extra)
  # Future eligibility expansion restores the original output reservation.
  testthat::local_mocked_bindings(cp04_retained_if_auxiliary_eligible = function(backend) TRUE, .package = "cdid")
  testthat::expect_identical(charge(alt_backend("cdid_dr_corollary1", general = TRUE), FALSE, N, Q, R), extra)
})

# Reuse the established public missing-panel fixture, including zero- and
# single-outcome units. The forced-stream reference is the unchanged installed
# producer, not the store being tested.
alt_fixture <- function() {
  d <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_baseline", "raw.csv"))
  g <- as.numeric(replace(d$cohort, d$cohort == "NT", "Inf")); g[is.infinite(g)] <- 0
  Y <- as.matrix(d[paste0("Y", 1:4)]); ids <- paste0("unit-", d$id)
  data <- data.frame(id = rep(ids, each = 4L), time = rep(1:4, nrow(d)),
    y = as.vector(t(Y)), g = rep(g, each = 4L))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  data <- data[!data$id %in% ids[zero[-1L]], , drop = FALSE]
  list(data = data, frame = data.frame(id = ids, g = g, x = d$x, w = d$w),
    zero = zero, single = which(rowSums(!is.na(Y)) == 1L))
}
alt_run <- function(fx, pooling, stream = FALSE) {
  reserve <- alt_ns("cp04_retained_if_reserve"); release <- alt_ns("cp04_retained_if_release")
  report <- alt_ns("build_ipw_auxiliary_reporting"); trace <- new.env(parent = emptyenv())
  trace$events <- character(); trace$args <- NULL; trace$entries <- NULL
  testthat::local_mocked_bindings(
    cp04_retained_if_reserve = function(backend, links, chunk_rows, cov_plan, retained, gmm_retained, live, work, control) {
      caller <- parent.frame()
      trace$args <- list(backend = backend, links = links, chunk_rows = chunk_rows, cov_plan = cov_plan,
        retained = retained, gmm_retained = gmm_retained, live = live, work = work, control = control)
      trace$extra <- 8 * get("n_total", caller) *
        (length(get("aux_sel_bp", caller)) + nrow(get("ledger", caller)))
      trace$common <- get("combined_live", caller) + get("final_diagnostic_copy_entries", caller) +
        max(0, get("cp04_work", caller)$peak + backend$accounting$owned_entries +
          backend$accounting$borrowed_entries - get("point_preflight", caller)$producer_workspace_entries)
      ans <- if (stream) NULL else reserve(backend, links, chunk_rows, cov_plan, retained, gmm_retained, live, work, control)
      if (!is.null(ans)) trace$entries <- ans$entries
      ans
    },
    cp04_retained_if_release = function(store) {
      if (is.environment(store) && !identical(store$state, "released")) {
        trace$events <- c(trace$events, "release")
        testthat::expect_false(exists("target_if_A", parent.frame(), inherits = FALSE))
        testthat::expect_null(store$auxiliary)
        trace$raw <- store$phi
      }
      release(store)
    },
    build_ipw_auxiliary_reporting = function(...) {
      trace$events <- c(trace$events, "report")
      report(...)
    }, .package = "cdid")
  fit <- cdid::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    data = fx$data, frame = fx$frame, xformla = ~x, weightsname = "w", est_method = "dr",
    outcome_pooling = "separate", selection = "conditional_mtar", selection_xformla = ~x,
    selection_model = "pair_logit", selection_pooling = pooling, selection_universe = "F",
    control_group = "nevertreated", gmm_weighting = "identity", control = list(chunk_rows = 7L),
    bstrap = FALSE, cband = FALSE)
  list(fit = fit, trace = trace)
}

testthat::test_that("the whole store admits at its true lifetime cap and releases before future output", {
  fx <- alt_fixture(); reserve <- alt_ns("cp04_retained_if_reserve")
  for (pooling in c("separate", "within_cohort")) {
    got <- alt_run(fx, pooling); ref <- alt_run(fx, pooling, stream = TRUE)
    tr <- got$trace; args <- tr$args
    testthat::expect_gt(tr$extra, 0)
    testthat::expect_equal(args$live, tr$common)
    testthat::expect_equal(tr$common + tr$extra - args$live, tr$extra)
    testthat::expect_gt(tr$entries, 0)
    testthat::expect_identical(tr$events, c("release", "report"))
    testthat::expect_identical(nrow(tr$raw), nrow(fx$frame))
    testthat::expect_true(all(is.finite(tr$raw[c(fx$zero, fx$single), , drop = FALSE])))
    # Exact adjacent main-store capacity; original planner probes still run.
    cap <- args$live + tr$entries
    args$control$max_matrix_entries <- cap; args$control$max_retained_entries <- cap
    admitted <- do.call(reserve, args)
    testthat::expect_true(is.environment(admitted))
    testthat::expect_identical(admitted$B, args$chunk_rows)
    testthat::expect_equal(admitted$admission_total, cap)
    alt_ns("cp04_retained_if_release")(admitted)
    args$control$max_matrix_entries <- cap - 1
    testthat::expect_null(do.call(reserve, args))
    args$control$max_matrix_entries <- cap; args$live <- args$live + tr$extra
    testthat::expect_null(do.call(reserve, args))
    # Full raw core IF against the unchanged ordinary producer, then all
    # public original-F target rows, covariance, points and status fields.
    ordinary <- matrix(NA_real_, nrow(tr$raw), ncol(tr$raw))
    alt_ns("cp04_replay_links")(tr$args$backend, tr$args$links, function(b)
      ordinary[b$rows, ] <<- b$phi, scale = "core", chunk_rows = tr$args$chunk_rows)
    testthat::expect_equal(unname(tr$raw), ordinary, tolerance = 1e-12)
    testthat::expect_equal(crossprod(unname(tr$raw)), crossprod(ordinary), tolerance = 1e-12)
    for (field in c("att", "se", "V_analytical", "inffunc"))
      testthat::expect_equal(got$fit[[field]], ref$fit[[field]], tolerance = 1e-12)
    testthat::expect_identical(rownames(got$fit$inffunc), fx$frame$id)
    testthat::expect_identical(got$fit$chained$available, ref$fit$chained$available)
    testthat::expect_identical(cdid::tidy(got$fit)[c("available", "support.status", "support.reason")],
      cdid::tidy(ref$fit)[c("available", "support.status", "support.reason")])
  }
})
