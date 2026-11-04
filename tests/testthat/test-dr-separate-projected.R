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
    plan <- dsp("chained_full_if_plan")(fit, links, 50L)
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
  plan <- dsp("chained_full_if_plan")(fit, links, 16L)
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
  partial_plan <- dsp("chained_full_if_plan")(fit, links[subset], 50L)
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
