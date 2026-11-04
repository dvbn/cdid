# Unchanged detailed replay and independent raw-row nuisance refits are the
# references. No expected value is produced by the projected helper itself.
dsp <- function(name) getFromNamespace(name, "cdid")
dsp_case <- function(intercept = TRUE, slopes = TRUE, alias = FALSE, chunk_rows = 7L,
                     copies = 1L, long_ids = FALSE) {
  j <- rep(1:12, 4); group <- rep(1:4, each = 12); N <- length(j) + 2L
  id <- seq_len(N); g <- c(c(3, 4, 5, Inf)[group], 3, Inf)
  x <- c((j - 6.5) / 5, -.7, .9); z <- c(sin(j) / 2 + j / 30, .2, -.3)
  w <- c((1 + (j %% 4) / 10) * (1 + group * x[1:48] / 30), 1.2, 1.7)
  Y <- vapply(1:5, function(t) c(j / 7 + group * t / 5 + (.2 * t + .03 * t^2) * x[1:48] +
    sin(j + group * t) / 9 + as.double(t >= g[1:48]) * (t + x[1:48] * group / 6), NA_real_, NA_real_), numeric(N))
  Y[c(which(j == 1L), which(j == 4L)), 2L] <- NA_real_
  Y[which(j == 3L), 4L] <- NA_real_
  Y[N, 1L] <- 4.1 # one single-outcome unit; row N-1 has no outcomes
  X <- matrix(numeric(), N, 0L)
  if (intercept) X <- cbind(`(Intercept)` = rep(1, N))
  if (slopes) X <- cbind(X, x = x)
  if (alias) X <- cbind(X, alias = 3 + 2 * x)
  V <- cbind(`(Intercept)` = 1, z = z)
  Xa <- cbind(`(Intercept)` = 1, x = x, z = z)
  if (alias) Xa <- cbind(Xa, alias = X[, "alias"])
  control <- dsp("ipw_owner_control")(chunk_rows, 1e8, 2e8, 1e13, 1e8)
  data <- data.frame(id = rep(id, each = 5L), t = rep(1:5, N), g = rep(g, each = 5L), y = as.vector(t(Y)))
  # Exercise both absent rows and stored NA outcomes without changing A.
  data <- data[!(data$id == 2L & data$t == 3L), ]
  Y[2L, 3L] <- NA_real_
  panel <- dsp("construct_ipw_panel_owner")(data, "id", "t", "y", "g", control)
  frame <- dsp("construct_ipw_frame_owner")(panel, id, Xa, w, control)
  reader <- dsp("construct_cp04_frame_reader")(frame, id,
    matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, control)
  defs <- list(list(s = 1, t = 2, g = 3, K = c(4, 5, Inf)),
    list(s = 1, t = 3, g = 3, K = c(4, 5, Inf)), list(s = 2, t = 3, g = 3, K = c(4, 5, Inf)),
    list(s = 1, t = 3, g = 4, K = c(5, Inf)), list(s = 1, t = 3, g = 5, K = Inf),
    list(s = 2, t = 5, g = 3, K = c(4, Inf)), list(s = 1, t = 3, g = 3, K = Inf),
    list(s = 1, t = 3, g = 4, K = Inf), list(s = 1, t = 3, g = 3, K = c(4, 5, Inf)))
  defs <- rep(defs, copies)
  views <- lapply(seq_along(defs), function(i) {
    d <- defs[[i]]
    key <- paste0("projected-", i, if (long_ids) strrep("long-id", 80L) else "")
    dsp("construct_ipw_link_view")(frame, d$s, d$t, d$g, d$K, control = control,
      link_id = key, functional_key = paste0("cmpcat:", key))
  })
  recipe <- dsp("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  fit <- dsp("fit_dr_baseline")(reader, views, control,
    list(epsilon = 1e-13, maxit = 60L, trace = FALSE), 1e-10, recipe,
    treatment_columns = colnames(V), outcome_pooling = "separate")
  list(fit = fit, reader = reader, Y = Y, X = X, V = V, w = w, g = g, defs = defs,
    intercept = intercept, N = N)
}
dsp_eval <- function(a, rows, links, details = FALSE, scale = "frame")
  dsp("evaluate_dr_separate_if_chunk")(a$fit, rows, links, scale, details)

# Directly refit each requested raw control WLS and pair-risk logit under
# perturbed original-unit weights, recomputing the empirical A-to-pair rates.
# Only base stats fits and the fixture's own data/definitions enter this truth.
dsp_raw <- function(a, frequency = rep(1, a$N)) {
  w <- a$w * frequency; A <- rowSums(!is.na(a$Y)) >= 2L
  vapply(a$defs, function(d) {
    obs <- !is.na(a$Y[, d$s]) & !is.na(a$Y[, d$t]); cohorts <- c(d$g, d$K)
    rates <- vapply(cohorts, function(k) sum(w[a$g == k & obs]) / sum(w[a$g == k & A]), numeric(1))
    rows <- which(obs & a$g %in% cohorts)
    q <- w[rows] / rates[match(a$g[rows], cohorts)]
    dy <- a$Y[rows, d$t] - a$Y[rows, d$s]
    D <- as.double(a$g[rows] == d$g)
    design <- a$X[rows, , drop = FALSE]
    if (a$intercept) design[, -1L] <- design[, -1L, drop = FALSE] * (d$t - d$s) else
      design <- design * (d$t - d$s)
    or_fit <- stats::lm.wfit(design[D == 0, , drop = FALSE], dy[D == 0], q[D == 0], tol = 1e-10)
    stopifnot(all(is.finite(or_fit$coefficients)))
    residual <- dy - as.vector(design %*% or_fit$coefficients)
    trt <- suppressWarnings(stats::glm.fit(a$V[rows, , drop = FALSE], D, weights = q,
      start = rep(0, ncol(a$V)), family = stats::binomial(),
      control = stats::glm.control(epsilon = 1e-13, maxit = 60L), intercept = FALSE))
    stopifnot(trt$converged, all(is.finite(trt$coefficients)))
    prob <- stats::plogis(as.vector(a$V[rows, , drop = FALSE] %*% trt$coefficients))
    stopifnot(all(prob > 1e-6 & prob < 1 - 1e-6))
    weighted.mean(residual, q * D) - weighted.mean(residual, q * (1 - D) * prob / (1 - prob))
  }, numeric(1))
}

testthat::test_that("ordinary incomplete separate DR preserves full coefficient effects and frozen operators", {
  for (options in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE))) {
    a <- do.call(dsp_case, options); fit <- a$fit; L <- length(fit$links)
    testthat::expect_null(fit$separate_complete)
    frozen <- list(groups = fit$outcome_groups, regressions = fit$regressions, treatments = fit$treatments)
    full <- matrix(0, a$N, L); detailed <- full
    for (first in seq.int(1L, a$N, by = 7L)) {
      rows <- seq.int(first, min(a$N, first + 6L))
      ordinary <- dsp_eval(a, rows, seq_len(L))
      old <- dsp_eval(a, rows, seq_len(L), TRUE)
      full[rows, ] <- ordinary$phi; detailed[rows, ] <- old$phi
      testthat::expect_equal(ordinary$phi, old$phi, tolerance = 1e-10)
      testthat::expect_identical(ordinary$id, old$id)
      testthat::expect_identical(ordinary$normalization_count, old$normalization_count)
      # Detailed replay intentionally has larger output/scratch/work bounds.
      # Compare the semantic counts exactly and check each output's own bound.
      bound_fields <- c("output_entries_bound", "concurrent_scratch_bound", "work_ops_bound")
      semantic_fields <- setdiff(names(old$diagnostics), bound_fields)
      testthat::expect_identical(names(ordinary$diagnostics), names(old$diagnostics))
      testthat::expect_identical(ordinary$diagnostics[semantic_fields], old$diagnostics[semantic_fields])
      testthat::expect_lte(unname(dsp("cp04_payload")(ordinary)["entries"]), ordinary$diagnostics$output_entries_bound)
      testthat::expect_lte(unname(dsp("cp04_payload")(old)["entries"]), old$diagnostics$output_entries_bound)
    }
    for (rows in list(7:1, c(50L, 49L, 2L, 1L), 23:19)) for (links in list(c(9L, 3L, 1L), c(8L, 5L, 7L))) {
      ordinary <- dsp_eval(a, rows, links)
      old <- dsp_eval(a, rows, links, TRUE)
      testthat::expect_equal(ordinary$phi, old$phi, tolerance = 1e-10)
      testthat::expect_equal(unname(ordinary$phi), full[rows, links, drop = FALSE], tolerance = 1e-10)
      core <- dsp_eval(a, rows, links, scale = "core")
      testthat::expect_equal(core$phi, a$reader$state$lambda * ordinary$phi, tolerance = 1e-10)
    }
    testthat::expect_equal(crossprod(full) / a$N, crossprod(detailed) / a$N, tolerance = 1e-10)
    testthat::expect_identical(list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments), frozen)
  }
})

testthat::test_that("projected incomplete DR matches independent raw case-weight score derivatives", {
  a <- dsp_case(); fit <- a$fit; L <- length(fit$links)
  testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), dsp_raw(a), tolerance = 1e-9)
  actual <- matrix(0, a$N, L)
  for (first in seq.int(1L, a$N, by = 7L)) {
    rows <- seq.int(first, min(a$N, first + 6L))
    actual[rows, ] <- dsp_eval(a, rows, seq_len(L))$phi
  }
  numerical <- t(vapply(seq_len(a$N), function(i) {
    direction <- rep(-1, a$N); direction[i] <- a$N - 1
    (dsp_raw(a, 1 + 1e-6 * direction) - dsp_raw(a, 1 - 1e-6 * direction)) / 2e-6
  }, numeric(L)))
  testthat::expect_equal(actual, numerical, tolerance = 2e-6)
  testthat::expect_equal(crossprod(actual) / a$N, crossprod(numerical) / a$N, tolerance = 2e-6)
  # Explicitly retain zero- and single-outcome rows in the reference universe.
  testthat::expect_identical(nrow(actual), a$N)
  testthat::expect_equal(actual[49:50, ], numerical[49:50, ], tolerance = 2e-6)
})

testthat::test_that("ordinary baseline projected replay avoids the detailed OR worker", {
  a <- dsp_case(); original <- dsp("or_baseline_separate_if_core")
  n <- 0L
  testthat::local_mocked_bindings(or_baseline_separate_if_core = function(...) {
    n <<- n + 1L; original(...)
  }, .package = "cdid")
  dsp_eval(a, 7:1, c(9L, 3L, 1L))
  testthat::expect_identical(n, 0L)
  dsp_eval(a, 7:1, c(9L, 3L, 1L), TRUE)
  testthat::expect_gt(n, 0L)
})

testthat::test_that("frozen separate DR replay operators match V8 and detailed IFs without chunk solves", {
  for (args in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE))) {
    a <- do.call(dsp_case, c(args, list(chunk_rows = 50L)))
    fit <- a$fit; links <- rev(seq_along(fit$links)); rows <- rev(seq_len(a$N))
    frozen <- list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, reader_ref = fit$reader_ref)
    plan <- dsp("dr_separate_baseline_replay_plan")(fit, links, 50L)
    testthat::expect_true(is.list(plan))
    testthat::expect_lte(unname(dsp("cp04_payload")(plan)["entries"]), plan$entries)
    testthat::expect_lte(plan$entries, plan$retained_bound)
    testthat::expect_lte(plan$reservation, 16 * plan$min_B * length(links))
    run <- function(plan) {
      pooled <- dsp("or_pooled_wls_gram_solve"); plain <- dsp("or_wls_gram_solve"); chol <- dsp("cp04_solve")
      solves <- 0L
      testthat::local_mocked_bindings(
        or_pooled_wls_gram_solve = function(...) { solves <<- solves + 1L; pooled(...) },
        or_wls_gram_solve = function(...) { solves <<- solves + 1L; plain(...) },
        cp04_solve = function(...) { solves <<- solves + 1L; chol(...) }, .package = "cdid")
      out <- dsp("chained_full_if_chunk")(fit, rows, links, "core", plan)
      list(out = out, solves = solves)
    }
    planned <- run(plan); original <- run(NULL)
    detailed <- dsp_eval(a, rows, links, TRUE, "core")
    testthat::expect_identical(planned$solves, 0L)
    testthat::expect_gt(original$solves, 0L)
    testthat::expect_equal(planned$out, original$out, tolerance = 1e-10)
    testthat::expect_equal(planned$out$phi, detailed$phi, tolerance = 1e-10)
    testthat::expect_equal(crossprod(planned$out$phi), crossprod(original$out$phi), tolerance = 1e-10)
    testthat::expect_identical(list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, reader_ref = fit$reader_ref), frozen)
    broken <- plan; broken$eta[1L, 1L] <- broken$eta[1L, 1L] + 1
    testthat::expect_error(dsp("chained_full_if_chunk")(fit, rows, links, "core", broken),
      class = "cdid_dr_baseline_integrity_error")
    testthat::expect_error(dsp("chained_full_if_chunk")(fit, rows, rev(links), "core", plan),
      class = "cdid_dr_baseline_integrity_error")
    testthat::expect_error(dsp("chained_full_if_chunk")(fit, rows[-1L], links, "core", plan),
      class = "cdid_dr_baseline_integrity_error")
    # The one-row tail cannot fund the retained descriptor, so the whole
    # replay declines without selecting a different row partition.
    testthat::expect_null(dsp("chained_full_if_plan")(fit, links, 49L))
    testthat::expect_null(dsp("chained_full_if_plan")(fit, 1L, 1L))
  }
  foreign <- dsp_case(chunk_rows = 50L)
  testthat::expect_error(dsp("chained_full_if_chunk")(foreign$fit, rows, links, "core", plan),
    class = "cdid_dr_baseline_integrity_error")
})

testthat::test_that("one compact plan serves duplicate consumers and the admitted short tail", {
  a <- dsp_case(chunk_rows = 50L, copies = 40L, long_ids = TRUE)
  fit <- a$fit; links <- rev(seq_along(fit$links))
  plan <- dsp("dr_separate_baseline_replay_plan")(fit, links, 16L)
  testthat::expect_true(is.list(plan))
  testthat::expect_equal(plan$min_B, 2)
  testthat::expect_lte(plan$reservation, 16 * 2L * length(links))
  testthat::expect_lte(length(plan$tids), 9L)
  testthat::expect_lte(length(plan$gids), 9L)
  # Labels stay in the existing backend/output; no label or fit snapshots
  # are replicated in this retained descriptor.
  testthat::expect_false(any(vapply(plan[setdiff(names(plan), "kind")], is.character, logical(1L))))
  testthat::expect_false(any(c("backend_snapshot", "reader_snapshot", "fit", "labels") %in% names(plan)))
  direct <- matrix(0, a$N, length(links)); replayed <- direct
  for (first in seq.int(1L, a$N, by = 16L)) {
    rows <- seq.int(first, min(a$N, first + 15L))
    out <- dsp("chained_full_if_chunk")(fit, rows, links, "core", plan)
    old <- dsp_eval(a, rows, links, FALSE, "core")
    testthat::expect_equal(out, old, tolerance = 1e-10)
    direct[rows, ] <- old$phi
  }
  subset <- seq.int(1L, length(links), by = 2L)
  partial_plan <- dsp("dr_separate_baseline_replay_plan")(fit, links[subset], 50L)
  testthat::expect_true(is.list(partial_plan))
  partial <- dsp("chained_full_if_chunk")(fit, rev(seq_len(a$N)), links[subset], "core", partial_plan)
  testthat::expect_equal(unname(partial$phi), direct[rev(seq_len(a$N)), subset, drop = FALSE], tolerance = 1e-10)
  create <- dsp("dr_separate_baseline_replay_plan"); created <- 0L
  testthat::local_mocked_bindings(dr_separate_baseline_replay_plan = function(...) {
    created <<- created + 1L; create(...)
  }, .package = "cdid")
  dsp("cp04_replay_links")(fit, links, function(out) replayed[out$rows, ] <<- out$phi,
    scale = "core", chunk_rows = 16L)
  testthat::expect_identical(created, 1L)
  testthat::expect_equal(replayed, direct, tolerance = 1e-10)
  pre <- dsp("dr_baseline_if_budget")(fit, 16L, links, FALSE)
  owned <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  retained <- fit$control$max_matrix_entries - owned - pre$output - pre$scratch
  admitted <- dsp("cp04_preflight_links")(fit, links, retained, chunk_rows = 16L)
  testthat::expect_equal(admitted$chunk_rows, 16)
  testthat::expect_error(dsp("cp04_preflight_links")(fit, links, retained + 1,
    chunk_rows = 16L), class = "cdid_ipw_owner_bound_error")
})

dsp_phase_or <- function(a) {
  frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
  views <- lapply(seq_along(a$defs), function(i) {
    d <- a$defs[[i]]
    dsp("construct_ipw_link_view")(frame, d$s, d$t, d$g, d$K, control = a$fit$control,
      link_id = paste0("phase-or-", i), functional_key = paste0("cmpcat:phase-or-", i))
  })
  recipe <- dsp("or_pooled_outcome_recipe")(colnames(frame$Z_a), TRUE,
    setdiff(colnames(frame$Z_a), "(Intercept)"))
  dsp("fit_or_baseline")(a$reader, views, a$fit$control, 1e-10,
    recipe, outcome_pooling = "separate")
}

testthat::test_that("baseline phase eligibility excludes pooled, selection, reference and selected complete routes", {
  # Eligibility-only metadata, never validated or evaluated as fitted owners.
  meta <- structure(list(outcome_pooling = "separate"), class = "cdid_dr_baseline")
  eligible <- dsp("cp04_retained_if_phase_baseline")
  testthat::expect_true(eligible(meta))
  for (field in c("uses_selection_models", "reference_ref", "separate_complete", "outcome_pooling")) {
    other <- meta
    other[[field]] <- switch(field, uses_selection_models = TRUE,
      reference_ref = new.env(parent = emptyenv()), separate_complete = list(selected = TRUE),
      outcome_pooling = "pooled")
    testthat::expect_false(eligible(other))
  }
  class(meta) <- "cdid_imputation_backend"
  testthat::expect_false(eligible(meta))
})

testthat::test_that("static incomplete OR and DR phase stores preserve original-unit IFs and release before reuse", {
  a <- dsp_case()
  fits <- list(dr = a$fit, or = dsp_phase_or(a))
  testthat::expect_lt(a$reader$state$lambda, 1)
  links <- c(6L, 3L, 1L, 5L, 2L, 4L)
  d5 <- list(ledger_if_entries = 0, basis_if_entries = 0,
    ledger_enclosure_entries = 0, basis_enclosure_entries = 0)
  for (kind in names(fits)) {
    fit <- fits[[kind]]; ctrl <- fit$control
    frozen <- list(ref = fit$ref, reader_ref = fit$reader_ref, regressions = fit$regressions,
      outcome_groups = fit$outcome_groups, treatments = fit$treatments)
    testthat::expect_true(dsp("cp04_retained_if_phase_baseline")(fit))
    raw <- matrix(0, a$N, length(links))
    dsp("cp04_replay_links")(fit, links, function(b) raw[b$rows, ] <<- b$phi,
      scale = "core", chunk_rows = 7L)
    reference <- dsp("cp04_stream_covariance")(fit, links)
    # The independent reference outputs remain live in this test caller.
    caller <- unname(dsp("cp04_payload")(list(raw, reference))["entries"])
    cov_plan <- dsp("cp04_stream_plan")(fit, links, caller)
    point <- dsp("cp04_preflight_links")(fit, links, caller, chunk_rows = 7L)
    testthat::expect_length(cov_plan$blocks, 1L)
    testthat::expect_equal(cov_plan$chunk_rows, 7)
    cov_inv <- dsp("cp04_retained_if_phase_inventory")("covariance", fit, cov_plan, point,
      caller, caller, d5, length(links), length(fit$links))
    entries <- dsp("cp04_retained_if_entries")(a$N, length(links), 7L) +
      dsp("cp04_retained_if_boundary_cost")(fit, links, 7L)$entries
    small <- dsp("ipw_owner_control")(7L, cov_inv$live + entries - 1,
      cov_inv$live + entries - 1, ctrl$max_work_ops, ctrl$max_output_entries)
    # The original phase buffers fit this finite cap, while its optional
    # cache does not. No sealed owner or its control is mutated.
    testthat::expect_true(dsp("cp04_budget")(small,
      c(cov_inv$retained, fit$accounting$owned_entries, fit$accounting$borrowed_entries),
      c(cov_plan$peak, cov_inv$gram_entries), 0, 0, "uncached phase witness"))
    reserve <- function(phase, inv, control = ctrl, active = NULL, transfer = FALSE) {
      dsp("cp04_retained_if_phase_reserve")(fit, links, phase, 7L, cov_plan,
        retained = inv$retained, live = inv$live, work = 0, control = control,
        numerical_transfer = transfer, active_store = active)
    }
    testthat::expect_null(reserve("covariance", cov_inv, small))
    testthat::expect_null(reserve("covariance", cov_inv, transfer = TRUE))
    too_small <- dsp("ipw_owner_control")(7L, cov_inv$live - 1, cov_inv$live - 1,
      ctrl$max_work_ops, ctrl$max_output_entries)
    testthat::expect_error(dsp("cp04_budget")(too_small,
      c(cov_inv$retained, fit$accounting$owned_entries, fit$accounting$borrowed_entries),
      c(cov_plan$peak, cov_inv$gram_entries), 0, 0, "original overbudget phase"),
      class = "cdid_ipw_owner_bound_error")
    cov_store <- reserve("covariance", cov_inv)
    testthat::expect_true(is.environment(cov_store))
    producer <- dsp("cp04_retained_if_bind")(fit, cov_store)
    cached <- dsp("cp04_stream_covariance")(producer, links, cov_inv$retained)
    testthat::expect_equal(cached$estimates, reference$estimates, tolerance = 1e-10)
    testthat::expect_equal(cached$mean_phi, reference$mean_phi, tolerance = 1e-10)
    testthat::expect_equal(cached$covariance, reference$covariance, tolerance = 1e-10)
    testthat::expect_equal(unname(cov_store$phi), raw, tolerance = 1e-10)
    testthat::expect_identical(cov_store$producer_passes, 1L)
    testthat::expect_identical(cov_store$slice_passes, 1L)
    testthat::expect_equal(unname(cov_store$phi[49:50, , drop = FALSE]),
      raw[49:50, , drop = FALSE], tolerance = 1e-10)
    testthat::expect_error(dsp("cp04_retained_if_check")(cov_store, fit, rev(links), 7L, "core"),
      class = "cdid_cp04_integrity_error")
    testthat::expect_error(dsp("cp04_retained_if_check")(cov_store, fit, links, 6L, "core"),
      class = "cdid_cp04_integrity_error")
    foreign <- fits[[if (kind == "dr") "or" else "dr"]]
    testthat::expect_error(dsp("cp04_retained_if_check")(cov_store, foreign, links, 7L, "core"),
      class = "cdid_cp04_integrity_error")
    gmm_inv <- dsp("cp04_retained_if_phase_inventory")("gmm", fit, cov_plan, point,
      caller, caller, d5, length(links), length(fit$links))
    testthat::expect_error(reserve("gmm", gmm_inv, active = cov_store),
      class = "cdid_cp04_integrity_error")
    cov_record <- dsp("cp04_retained_if_release")(cov_store)
    producer <- NULL
    testthat::expect_identical(cov_store$state, "released")
    testthat::expect_false(any(c("phi", "raw_ref", "backend_snapshot", "reader_snapshot") %in% ls(cov_store)))
    gmm_store <- reserve("gmm", gmm_inv, active = cov_store)
    testthat::expect_true(is.environment(gmm_store))
    producer <- dsp("cp04_retained_if_bind")(fit, gmm_store)
    for (pass in 1:2) dsp("cp04_replay_links")(producer, links, function(b) {
      testthat::expect_equal(unname(b$phi), raw[b$rows, , drop = FALSE], tolerance = 1e-10)
    }, scale = "core", chunk_rows = 7L)
    testthat::expect_identical(gmm_store$producer_passes, 1L)
    testthat::expect_identical(gmm_store$slice_passes, 1L)
    gmm_record <- dsp("cp04_retained_if_release")(gmm_store)
    producer <- NULL
    testthat::expect_false(any(c("phi", "raw_ref", "backend_snapshot", "reader_snapshot") %in% ls(gmm_store)))
    record <- dsp("cp04_retained_if_phase_record")(cov_record, gmm_record)
    testthat::expect_identical(record$mode, "phase_local")
    testthat::expect_identical(record$simultaneous_stores, 1L)
    testthat::expect_identical(record$producer_passes, 2L)
    testthat::expect_identical(record$slice_passes, 2L)
    testthat::expect_identical(list(ref = fit$ref, reader_ref = fit$reader_ref, regressions = fit$regressions,
      outcome_groups = fit$outcome_groups, treatments = fit$treatments), frozen)
  }
})


# V28 exact original V11 planned transport, independent of the new C worker.
dsp28_native <- function(a, plan, rows, data = plan$data, op = plan$v11, y = NULL, cap = NULL) {
  fit <- a$fit; reader <- a$reader; st <- reader$state
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(y)) y <- panel$observations$y
  if (is.null(cap)) cap <- dsp("dr_separate_sparse_live")(length(rows), length(plan$sids),
    length(plan$dates), length(plan$tids), length(plan$links))
  .Call(dsp("or_baseline_native_kernel")("C_gmm_dr_separate_sparse_if"), y,
    panel$postings$post_a_index, panel$postings$post_obs_index, as.integer(rows), st$cohort_index,
    st$a_index, st$w, plan$bounds, frame$Z_a, data, op, as.double(cap))
}

testthat::test_that("native projected separate DR matches frozen V11, local bases and raw derivatives", {
  for (args in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE))) {
    a <- do.call(dsp_case, c(args, list(chunk_rows = 50L))); fit <- a$fit
    links <- rev(seq_along(fit$links)); rows <- rev(seq_len(a$N))
    old <- dsp("dr_separate_baseline_replay_plan")(fit, links, 50L)
    plan <- dsp("chained_full_if_plan")(fit, links, 50L)
    testthat::expect_identical(plan$kind, "dr_separate_sparse_plan")
    testthat::expect_true(dsp("ipw_owner_seal_ok")(plan, "dr_separate_sparse_plan"))
    # Actual borrowed identity, not merely equal coefficients or fitted weights.
    for (field in c("or_adj", "or_center", "or_up", "eta", "prop_up"))
      testthat::expect_true(dsp("gmm_sexp_aliases")(plan$v11[[field]],
        get("dr_separate_baseline_replay_plan", plan$v11$ref, inherits = FALSE)[[field]]))
    frozen <- list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, state = a$reader$state)
    direct <- dsp28_native(a, plan, rows)
    testthat::expect_identical(direct$status[1L], 0L)
    out <- dsp("chained_full_if_chunk")(fit, rows, links, "frame", plan)
    ref <- dsp("chained_full_if_chunk")(fit, rows, links, "frame", old)
    testthat::expect_equal(out, ref, tolerance = 1e-10)
    testthat::expect_equal(unname(out$phi), direct$phi, tolerance = 1e-10)
    testthat::expect_equal(out$phi, dsp_eval(a, rows, links, TRUE)$phi, tolerance = 1e-10)
    testthat::expect_equal(dsp("chained_full_if_chunk")(fit, rows, links, "core", plan)$phi,
      a$reader$state$lambda * out$phi, tolerance = 1e-10)
    testthat::expect_equal(crossprod(out$phi) / a$N, crossprod(ref$phi) / a$N, tolerance = 1e-10)
    testthat::expect_equal(out$phi[match(49:50, rows), ], ref$phi[match(49:50, rows), ], tolerance = 1e-10)
    testthat::expect_identical(list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, state = a$reader$state), frozen)
    pre <- dsp("dr_baseline_if_budget")(fit, 50L, links, FALSE)
    testthat::expect_lte(plan$entries, plan$retained_bound)
    testthat::expect_gte(plan$reservation, plan$v11$reservation)
    testthat::expect_lte(plan$reservation + dsp("dr_separate_sparse_inventory")(50L,
      length(plan$sids), length(plan$dates), length(plan$tids), length(links), plan$metadata), pre$scratch)
    if (!length(args)) {
      numerical <- t(vapply(seq_len(a$N), function(i) {
        direction <- rep(-1, a$N); direction[i] <- a$N - 1
        (dsp_raw(a, 1 + 1e-6 * direction) - dsp_raw(a, 1 - 1e-6 * direction)) / 2e-6
      }, numeric(length(links))))
      testthat::expect_equal(unname(out$phi[match(seq_len(a$N), rows), ]),
        numerical[, links, drop = FALSE], tolerance = 2e-6)
      # Deliberately remove one actual rate contraction from the old operators.
      no_rate <- plan$v11; no_rate$or_up[] <- 0; no_rate$prop_up[] <- 0
      testthat::expect_gt(max(abs(dsp28_native(a, plan, rows, op = no_rate)$phi - direct$phi)), 1e-6)
      # The centering channel is an ALL-F scalar, including both off-A units.
      shifted <- plan$v11; shift <- seq_along(shifted$or_center) / 100
      shifted$or_center <- shifted$or_center + shift
      expected <- direct$phi
      for (j in seq_along(links)) {
        t <- match(fit$links[[links[j]]]$model, plan$tids)
        expected[, j] <- expected[, j] + shift[t]
      }
      testthat::expect_equal(dsp28_native(a, plan, rows, op = shifted)$phi, expected, tolerance = 1e-10)
      testthat::expect_gt(max(abs(expected[match(49:50, rows), ] - direct$phi[match(49:50, rows), ])), 1e-6)
    }
  }
})

testthat::test_that("projected native DR preserves finite-phi negative-logit and unused-coordinate contracts", {
  a <- dsp_case(chunk_rows = 50L); fit <- a$fit; links <- seq_along(fit$links); rows <- seq_len(a$N)
  # A deliberate frozen-operator numeric witness, NOT another fitted model:
  # raw eta overflows negatively on observed rows, while p=odds=0 and phi is
  # finite under the CURRENT V11 projected arithmetic.
  changed <- fit; changed$ref <- NULL
  for (j in seq_along(changed$treatments)) changed$treatments[[j]]$coefficients[] <- -.Machine$double.xmax
  changed <- dsp("ipw_owner_seal")(changed, "dr_baseline_backend")
  b <- a; b$fit <- changed
  old <- dsp("dr_separate_baseline_replay_plan")(changed, links, 50L)
  plan <- dsp("chained_full_if_plan")(changed, links, 50L)
  testthat::expect_identical(plan$kind, "dr_separate_sparse_plan")
  frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
  eta <- as.vector(frame$Z_a[, changed$treatment_recipe$index, drop = FALSE] %*% changed$treatments[[1L]]$coefficients)
  testthat::expect_true(any(is.infinite(eta) & eta < 0))
  ref <- dsp("chained_full_if_chunk")(changed, rows, links, "frame", old)
  testthat::expect_true(all(is.finite(ref$phi)))
  testthat::expect_equal(dsp("chained_full_if_chunk")(changed, rows, links, "frame", plan), ref, tolerance = 1e-10)
  # Fields not consumed by ordinary V11 may not acquire new finite guards.
  changed <- fit; changed$ref <- NULL
  for (j in seq_along(changed$treatments)) {
    changed$treatments[[j]]$partials$beta[c(1L, 3L, 4L), ] <- Inf
    changed$treatments[[j]]$partials$upstream[3:4, ] <- Inf
  }
  changed <- dsp("ipw_owner_seal")(changed, "dr_baseline_backend")
  old <- dsp("dr_separate_baseline_replay_plan")(changed, links, 50L)
  plan <- dsp("chained_full_if_plan")(changed, links, 50L)
  testthat::expect_identical(plan$kind, "dr_separate_sparse_plan")
  ref <- dsp("chained_full_if_chunk")(changed, rows, links, "frame", old)
  testthat::expect_true(all(is.finite(ref$phi)))
  testthat::expect_equal(dsp("chained_full_if_chunk")(changed, rows, links, "frame", plan), ref, tolerance = 1e-10)
  # p=0 returns its zero outcome projection before touching those operators.
  b <- dsp_case(chunk_rows = 50L, slopes = FALSE); links <- seq_along(b$fit$links)
  old <- dsp("dr_separate_baseline_replay_plan")(b$fit, links, 50L); old$ref <- NULL
  old$or_center[] <- Inf; old$or_up[] <- Inf
  old <- dsp("ipw_owner_seal")(old, "dr_separate_baseline_replay_plan")
  plan <- dsp("dr_separate_sparse_plan")(b$fit, links, 50L, old)
  testthat::expect_identical(plan$kind, "dr_separate_sparse_plan")
  ref <- dsp("chained_full_if_chunk")(b$fit, rows, links, "frame", old)
  testthat::expect_true(all(is.finite(ref$phi)))
  testthat::expect_equal(dsp("chained_full_if_chunk")(b$fit, rows, links, "frame", plan), ref, tolerance = 1e-10)
})

testthat::test_that("projected native DR covers subsets, short tails and exact source/capacity refusals", {
  a <- dsp_case(chunk_rows = 50L, copies = 40L, long_ids = TRUE); fit <- a$fit
  links <- rev(seq_along(fit$links)); plan <- dsp("chained_full_if_plan")(fit, links, 16L)
  testthat::expect_identical(plan$kind, "dr_separate_sparse_plan")
  testthat::expect_equal(plan$min_B, 2)
  testthat::expect_lte(length(plan$tids), 9L); testthat::expect_lte(length(plan$gids), 9L)
  testthat::expect_false(any(c("backend_snapshot", "reader_snapshot", "fit", "labels") %in% names(plan)))
  full <- matrix(0, a$N, length(links))
  for (first in seq.int(1L, a$N, by = 16L)) {
    rows <- seq.int(first, min(a$N, first + 15L))
    out <- dsp("chained_full_if_chunk")(fit, rows, links, "core", plan)
    old <- dsp("chained_full_if_chunk")(fit, rows, links, "core", plan$v11)
    testthat::expect_equal(out, old, tolerance = 1e-10); full[rows, ] <- out$phi
  }
  sub <- seq.int(1L, length(links), by = 3L)
  partial <- dsp("chained_full_if_plan")(fit, links[sub], 50L)
  testthat::expect_identical(partial$kind, "dr_separate_sparse_plan")
  rows <- rev(seq_len(a$N))
  testthat::expect_equal(unname(dsp("chained_full_if_chunk")(fit, rows, links[sub], "core", partial)$phi),
    full[rows, sub, drop = FALSE], tolerance = 1e-10)
  broken <- plan; broken$v11$or_adj[1L] <- broken$v11$or_adj[1L] + .1
  testthat::expect_error(dsp("chained_full_if_chunk")(fit, 1:16, links, "core", broken),
    class = "cdid_dr_baseline_integrity_error")
  broken <- plan; broken$kind <- "bad-kind"
  testthat::expect_error(dsp("chained_full_if_chunk")(fit, 1:16, links, "core", broken),
    class = "cdid_dr_baseline_integrity_error")
  testthat::expect_error(dsp("chained_full_if_chunk")(fit, 1:15, links, "core", plan),
    class = "cdid_dr_baseline_integrity_error")
  testthat::expect_error(dsp("chained_full_if_chunk")(fit, 1:16, rev(links), "core", plan),
    class = "cdid_dr_baseline_integrity_error")
  broken <- fit; broken$treatments[[1L]]$coefficients[1L] <- broken$treatments[[1L]]$coefficients[1L] + .1
  testthat::expect_error(dsp("chained_full_if_chunk")(broken, 1:16, links, "core", plan),
    class = "cdid_dr_baseline_integrity_error")
  invalid <- plan$data; invalid$keep[1L] <- 0L
  testthat::expect_error(dsp28_native(a, plan, 1:16, invalid), class = "cdid_gmm_bound_error")
  testthat::expect_error(dsp28_native(a, plan, 1:16, cap = 16 * length(links) - 1), class = "cdid_gmm_bound_error")
  # Source endpoint and raw-difference refusals retain precedence even p=0.
  p0 <- dsp_case(chunk_rows = 50L, slopes = FALSE); ll <- seq_along(p0$fit$links)
  p <- dsp("chained_full_if_plan")(p0$fit, ll, 50L)
  testthat::expect_identical(p$kind, "dr_separate_sparse_plan")
  panel <- get("panel_owner", p0$reader$panel_ref, inherits = FALSE)
  y <- panel$observations$y; y[1L] <- Inf
  testthat::expect_identical(dsp28_native(p0, p, 1:50, y = y)$status[1L], 2L)
  s <- p0$fit$support_specs[[p$sids[1L]]]; st <- p0$reader$state
  cohort_rows <- which(st$cohort_index == s$cohort)
  observed <- which(!is.na(p0$Y[cohort_rows, s$pair$s_calendar_index]) & !is.na(p0$Y[cohort_rows, s$pair$t_calendar_index]))[1L]
  row <- cohort_rows[observed]; ai <- st$a_index[row]
  ix <- vapply(c(s$pair$s_calendar_index, s$pair$t_calendar_index), function(d) {
    r <- seq.int(panel$postings$date_start[d], panel$postings$date_end[d])
    panel$postings$post_obs_index[r[match(ai, panel$postings$post_a_index[r])]]
  }, integer(1))
  y <- panel$observations$y; y[ix] <- c(-.Machine$double.xmax, .Machine$double.xmax)
  testthat::expect_identical(dsp28_native(p0, p, 1:50, y = y)$status[1L], 3L)
  # Final-phi failure is still R's original cp04_numeric error and message.
  bad <- p; bad$ref <- NULL; bad$v11 <- p$v11; bad$v11$ref <- NULL
  bad$v11$eta[] <- NaN
  bad$v11 <- dsp("ipw_owner_seal")(bad$v11, "dr_separate_baseline_replay_plan")
  bad <- dsp("ipw_owner_seal")(bad, "dr_separate_sparse_plan")
  err <- tryCatch(dsp("chained_full_if_chunk")(p0$fit, 1:50, ll, "frame", bad), error = identity)
  testthat::expect_s3_class(err, "cdid_cp04_numeric_error")
  testthat::expect_identical(conditionMessage(err), "Nonfinite complete separate DR original-unit IF.")
  # Optional added lifetime cannot create a caller refusal. A source-bound
  # inventory pressure probe returns the exact existing V11 plan instead.
  v11 <- dsp("dr_separate_baseline_replay_plan")(fit, links, 16L)
  original <- dsp("dr_separate_sparse_inventory")
  pre <- dsp("dr_baseline_if_budget")(fit, 16L, links)
  testthat::local_mocked_bindings(dr_separate_sparse_inventory = function(...) pre$scratch + 1,
    .package = "cdid")
  declined <- dsp("dr_separate_sparse_plan")(fit, links, 16L, v11)
  testthat::expect_identical(declined, v11)
  expected_declined <- full[1:16, , drop = FALSE]
  dimnames(expected_declined) <- list(NULL,
    vapply(fit$links[links], `[[`, character(1), "link_id"))
  testthat::expect_equal(dsp("chained_full_if_chunk")(fit, 1:16, links, "core", declined)$phi,
    expected_declined, tolerance = 1e-10)
})
