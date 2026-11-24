# Memory-phase accounting of the direct same-pair `within` working covariance
# (`R/chained_cp04_memory_phases.R`, `R/chained_ipw_pipeline.R`). The phase
# route must admit fits that the whole-path envelope overcounts, return the
# same estimates as the whole-path route, and charge each `within` phase at
# least the vector heap that the phase actually allocates.

wmp_ns <- function(name) getFromNamespace(name, "cdid")

wmp_panel <- function(n = 60L, n_t = 16L, seed = 7L) {
  set.seed(seed)
  g <- sample(c(0, 4, 7, 10), n, TRUE)
  d <- expand.grid(t = seq_len(n_t), id = seq_len(n))
  d$g <- g[d$id]
  d$y <- stats::ave(stats::rnorm(nrow(d)), d$id, FUN = cumsum) + ifelse(d$g > 0 & d$t >= d$g, 1, 0)
  d
}

# `whole = TRUE` restores the whole-path envelope for `within` only.
wmp_fit <- function(d, bytes = NULL, whole = FALSE, method = "dr", law = "nevertreated") {
  if (whole) {
    eligible <- wmp_ns("cp04_pipeline_phase_eligible")
    testthat::local_mocked_bindings(cp04_pipeline_phase_eligible = function(cp04, transfer, pretest, imputation,
        augmented, custom, weighting, strategy) {
      if (identical(weighting, "within")) FALSE else
        eligible(cp04, transfer, pretest, imputation, augmented, custom, weighting, strategy)
    }, .package = "cdid")
  }
  ctl <- list(chunk_rows = 16L)
  if (!is.null(bytes)) ctl$numeric_workspace_bytes <- bytes
  cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~1, data = d,
    bstrap = FALSE, cband = FALSE, base_period = "universal", est_method = method,
    gmm_weighting = "within", control_group = law, control = ctl)
}

testthat::test_that("phase eligibility admits every CP04 within route and leaves optimal unchanged", {
  eligible <- wmp_ns("cp04_pipeline_phase_eligible")
  args <- list(cp04 = list(dr_method = "baseline", outcome_pooling = "separate"), transfer = FALSE,
    pretest = FALSE, imputation = FALSE, augmented = FALSE, custom = NULL, weighting = "within",
    strategy = list(method = "tikhonov"))
  testthat::expect_true(do.call(eligible, args))
  for (field in list(list(outcome_pooling = "pooled"), list(selection_est_method = "ipw"),
                     list(canonical_assumption = "conditional_mtar"), list(reference_source = list()))) {
    a <- args; a$cp04 <- utils::modifyList(a$cp04, field)
    testthat::expect_true(do.call(eligible, a))
    a$weighting <- "optimal"
    testthat::expect_identical(do.call(eligible, a), is.null(a$cp04$reference_source))
  }
  for (flag in c("imputation", "augmented")) {
    a <- args; a[[flag]] <- TRUE
    testthat::expect_true(do.call(eligible, a))
    a$weighting <- "optimal"
    testthat::expect_false(do.call(eligible, a))
  }
  for (flag in c("transfer", "pretest")) {
    a <- args; a[[flag]] <- TRUE
    testthat::expect_false(do.call(eligible, a))
  }
  a <- args; a["strategy"] <- list(NULL)
  testthat::expect_false(do.call(eligible, a))
  a <- args; a$weighting <- "identity"
  testthat::expect_false(do.call(eligible, a))
  a <- args; a["cp04"] <- list(NULL)
  testthat::expect_false(do.call(eligible, a))
})

testthat::test_that("phase accounting admits a within fit the whole-path envelope refuses, with identical results", {
  d <- wmp_panel()
  # Between the smallest admitted workspace of the phase route (about 224 MB
  # for this design) and that of the whole-path envelope (about 250 MB).
  bytes <- 236e6
  testthat::expect_error(wmp_fit(d, bytes, whole = TRUE), class = "cdid_ipw_owner_bound_error")
  admitted <- wmp_fit(d, bytes)
  reference <- wmp_fit(d, whole = TRUE)
  phase <- admitted$chained$resource_accounting$memory_phases
  testthat::expect_false(is.null(phase))
  testthat::expect_null(reference$chained$resource_accounting$memory_phases)
  within <- admitted$chained$resource_accounting$within_phases
  testthat::expect_gt(within$reduction_entries, 0)
  testthat::expect_gt(within$covariance_entries, 0)
  testthat::expect_lte(max(unlist(within)), bytes / 8)
  testthat::expect_lt(admitted$chained$resource_accounting$simultaneous_entries,
    reference$chained$resource_accounting$simultaneous_entries)
  testthat::expect_identical(admitted$att, reference$att)
  testthat::expect_identical(admitted$se, reference$se)
  testthat::expect_identical(admitted$inffunc, reference$inffunc)
  testthat::expect_identical(admitted$V_analytical, reference$V_analytical)
  testthat::expect_identical(admitted$chained$available, reference$chained$available)
  testthat::expect_identical(admitted$chained$covariance_handling, reference$chained$covariance_handling)
})

# The reachable vector heap is sampled by a full collection before every
# top-level statement of the two construction functions and of the forest
# assembly they call; the peak above the heap in use at entry must not exceed
# the corresponding phase charge. Temporaries inside a single statement are
# covered by the charge derivations in the source, not by this sampling.
testthat::test_that("each within phase charge bounds the reachable heap of its construction", {
  d <- wmp_panel(n = 80L, n_t = 12L, seed = 11L)
  seen <- new.env(parent = emptyenv())
  ns <- asNamespace("cdid")
  forest_entries <- wmp_ns("cp04_pipeline_within_forest_entries")
  covariance_entries <- wmp_ns("cp04_pipeline_within_covariance_entries")
  sample_heap <- function(phase) {
    used <- gc()["Vcells", "used"]
    if (is.null(seen[[phase]])) seen[[phase]] <- c(start = used, peak = used) else
      seen[[phase]]["peak"] <- max(seen[[phase]]["peak"], used)
  }
  # Sample from the first statement that allocates a link-sized array onward.
  steps <- function(name, from) {
    b <- as.list(body(get(name, ns)))
    first <- which(vapply(b, function(x) any(grepl(from, deparse(x), fixed = TRUE)), logical(1)))[1L]
    seq.int(first, length(b))
  }
  starts <- list(within_mask_spanning_forest = "masks_m <-", build_mask_forest_and_gram = "kept_chunks <-",
    compute_pipeline_within_covariance = "Omega_scaled <-")
  tracers <- list(within_mask_spanning_forest = "forest", build_mask_forest_and_gram = "forest",
    compute_pipeline_within_covariance = "covariance")
  for (fn in names(tracers)) {
    probe <- bquote(.(sample_heap)(.(tracers[[fn]])))
    suppressMessages(trace(fn, where = ns, at = steps(fn, starts[[fn]]), print = FALSE, tracer = probe, exit = probe))
  }
  on.exit(for (fn in names(tracers)) suppressMessages(untrace(fn, where = ns)), add = TRUE)
  testthat::local_mocked_bindings(
    cp04_pipeline_within_forest_entries = function(M, n_a, n_cal) {
      seen$forest_charge <- forest_entries(M, n_a, n_cal); seen$forest_charge
    },
    cp04_pipeline_within_covariance_entries = function(L, n_a, n_cal, p) {
      seen$covariance_charge <- covariance_entries(L, n_a, n_cal, p); seen$covariance_charge
    },
    .package = "cdid")
  for (route in list(c("dr", "nevertreated"), c("reg", "notyettreated"))) {
    rm(list = ls(seen), envir = seen)
    fit <- wmp_fit(d, method = route[1], law = route[2])
    testthat::expect_false(is.null(fit$chained$resource_accounting$memory_phases))
    testthat::expect_true(all(c("forest", "covariance", "forest_charge", "covariance_charge") %in% ls(seen)))
    forest_peak <- unname(seen$forest["peak"] - seen$forest["start"])
    covariance_peak <- unname(seen$covariance["peak"] - seen$covariance["start"])
    testthat::expect_gt(covariance_peak, 0)
    testthat::expect_lte(forest_peak, seen$forest_charge)
    testthat::expect_lte(covariance_peak, seen$covariance_charge)
  }
})

# Every term of every within phase charge, against an independent statement
# of the documented formulas.
testthat::test_that("within phase charges keep every counted term", {
  forest <- wmp_ns("cp04_pipeline_within_forest_entries")
  covariance <- wmp_ns("cp04_pipeline_within_covariance_entries")
  bundle <- wmp_ns("cp04_pipeline_within_bundle_entries")
  stream <- wmp_ns("cp04_pipeline_within_stream_entries")
  set.seed(31L)
  for (i in seq_len(25L)) {
    M <- sample.int(3000L, 1L); L <- sample.int(M, 1L); M_u <- M + sample.int(50L, 1L) - 1L
    n_a <- sample.int(5000L, 1L); n_cal <- sample.int(40L, 1L); p <- sample.int(12L, 1L)
    n_g <- sample.int(8L, 1L); b <- sample.int(n_a, 1L)
    testthat::expect_identical(forest(M, n_a, n_cal),
      (n_cal + M + 8) * n_a + 6 * M^2 + 64 * M + 1024)
    testthat::expect_identical(covariance(L, n_a, n_cal, p),
      (n_cal + 8) * n_a + 2 * n_a * p + 16 * n_a * (p + 1) + 4 * p^2 + L * n_a + 10 * L^2 + 16 * L + 1024)
    testthat::expect_identical(bundle(M, M_u, n_a, n_cal, n_g, p),
      2 * n_a * n_cal + 4 * M_u * n_a + 16 * n_a * n_cal * (p + 1) + 4 * p^2 + n_g * (n_cal^2 + n_a) +
        7 * M * n_a + 4 * n_a + 6 * M^2 + 64 * M + 64 * M_u + 1024)
    testthat::expect_identical(stream(L, b, n_a, M), 6 * L^2 + 4 * b * L + 4 * n_a + M + 16 * L + 1024)
  }
})

# The direct route's covariance phase composes the persistent prefix, the
# backend, the empirical covariance stream, any retained store and the
# construction charge; the general route's phases compose the prefix, the
# backend, the loading charge, the log and the stream charge, and the log is
# released before the solve.
testthat::test_that("within phase checks compose every live array", {
  seen <- new.env(parent = emptyenv()); seen$checks <- list()
  check <- wmp_ns("ipw_pipeline_within_phase_check")
  stream_cov <- wmp_ns("cp04_stream_covariance")
  make_log <- wmp_ns("level_loading_log")
  solve <- wmp_ns("solve_chained_gmm_core")
  block_rows <- wmp_ns("stream_block_rows")
  reserve <- wmp_ns("cp04_retained_if_reserve")
  testthat::local_mocked_bindings(
    ipw_pipeline_within_phase_check = function(nuis_res, retained, phase_entries, label, registry_prefix = TRUE) {
      total <- check(nuis_res, retained, phase_entries, label, registry_prefix)
      b <- nuis_res$cp04
      # A retained value that extends the prefix loses the recorded external
      # prefix and, once, the owner storage the prefix already holds.
      base <- if (registry_prefix) wmp_ns("cp04_composed_retained_entries")(nuis_res, retained) else retained
      overlap <- if (registry_prefix) wmp_ns("cp04_owner_overlap_entries")(nuis_res) else 0
      seen$checks[[length(seen$checks) + 1L]] <- list(label = label, phase = phase_entries, total = total,
        expected = base + b$accounting$owned_entries + b$accounting$borrowed_entries + sum(phase_entries),
        overlap = overlap, owned = b$accounting$owned_entries, registry_prefix = registry_prefix,
        max = b$control$max_matrix_entries)
      total
    },
    cp04_stream_covariance = function(...) { out <- stream_cov(...); seen$stream <- out; out },
    cp04_retained_if_reserve = function(...) { out <- reserve(...); seen$store <- out; out },
    level_loading_log = function(n_cal, n_a, L, limit = Inf) { seen$limit <- limit; make_log(n_cal, n_a, L, limit) },
    stream_block_rows = function(n_rows, per_row, max_entries, live_entries, cap = n_rows) {
      seen$blocks[[length(seen$blocks) + 1L]] <- c(per_row = per_row, live = live_entries)
      block_rows(n_rows, per_row, max_entries, live_entries, cap)
    },
    solve_chained_gmm_core = function(...) {
      bundle <- get0("same_pair_within_bundle", envir = parent.frame(), inherits = FALSE)
      seen$log_at_solve <- if (is.null(bundle)) "none" else if (is.null(bundle$U_log)) "released" else "live"
      solve(...)
    },
    .package = "cdid")
  d <- wmp_panel(n = 60L, n_t = 10L)
  # Direct route.
  seen$checks <- list(); seen$store <- NULL; fit <- wmp_fit(d)
  cov_check <- Filter(function(x) grepl("within covariance phase", x$label), seen$checks)
  testthat::expect_length(cov_check, 1L)
  cc <- cov_check[[1L]]
  testthat::expect_identical(cc$phase[1L], unname(wmp_ns("cp04_payload")(seen$stream)["entries"]))
  testthat::expect_gt(cc$phase[1L], length(seen$stream$covariance) - 1)
  # The retained moment IF store admitted on this route is charged in full.
  testthat::expect_false(is.null(seen$store))
  testthat::expect_gt(seen$store$entries, 0)
  testthat::expect_identical(cc$phase[2L], seen$store$entries)
  for (x in seen$checks) { testthat::expect_identical(x$total, x$expected); testthat::expect_lte(x$total, x$max)
    if (x$registry_prefix) testthat::expect_gte(x$overlap, x$owned) }
  testthat::expect_identical(seen$log_at_solve, "none")
  # General route (pooled outcome models).
  seen$checks <- list(); seen$limit <- NULL; seen$blocks <- list()
  ctl <- list(chunk_rows = 16L)
  fit <- cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g", xformla = ~1, data = d,
    bstrap = FALSE, cband = FALSE, base_period = "universal", est_method = "dr", outcome_pooling = "pooled",
    gmm_weighting = "within", control_group = "nevertreated", control = ctl)
  labels <- vapply(seen$checks, `[[`, character(1), "label")
  testthat::expect_true(all(c("CP04 within loading phase simultaneous entries",
    "CP04 within covariance phase simultaneous entries") %in% labels))
  load <- seen$checks[[which(labels == "CP04 within loading phase simultaneous entries")]]
  covp <- seen$checks[[which(labels == "CP04 within covariance phase simultaneous entries")]]
  for (x in seen$checks) { testthat::expect_identical(x$total, x$expected); testthat::expect_lte(x$total, x$max) }
  testthat::expect_identical(seen$limit, load$max - load$total)
  rec <- fit$chained$resource_accounting$within_phases
  testthat::expect_identical(covp$phase[1L], rec$loading_log_entries)
  testthat::expect_gt(rec$loading_log_entries, 0)
  testthat::expect_identical(seen$log_at_solve, "released")
  # The working-covariance blocks count the log and 5 L^2 as live.
  L <- as.double(fit$chained$covariance_handling$dimensions$L)
  testthat::expect_gt(L, 0)
  gram <- Filter(function(x) x[["per_row"]] == 4 * L, seen$blocks)
  testthat::expect_length(gram, 1L)
  testthat::expect_identical(unname(gram[[1L]][["live"]]), rec$loading_log_entries + 5 * L^2)
})

# One small design per general-route family: the whole-path envelope refuses
# at the stated workspace, the phase route admits it with identical results.
wmp_panel_random <- function(n = 100L, n_t = 12L, seed = 7L) {
  set.seed(seed)
  g <- sample(c(0, 4, 7, 10), n, TRUE)
  d <- expand.grid(t = seq_len(n_t), id = seq_len(n)); d$g <- g[d$id]
  d$xv <- round(stats::rnorm(n)[d$id] + 0.1 * d$t + stats::rnorm(nrow(d), sd = 0.3), 3)
  d$y <- stats::ave(stats::rnorm(nrow(d)), d$id, FUN = cumsum) + ifelse(d$g > 0 & d$t >= d$g, 1, 0)
  d$y[stats::runif(nrow(d)) < 0.25] <- NA
  d
}

wmp_panel_fixed <- function(n = 60L, n_t = 16L, seed = 7L) {
  set.seed(seed)
  g <- sample(c(0, 4, 7, 10), n, TRUE)
  d <- expand.grid(t = seq_len(n_t), id = seq_len(n)); d$g <- g[d$id]
  d$xv <- round(stats::rnorm(n)[d$id] + 0.1 * d$t + stats::rnorm(nrow(d), sd = 0.3), 3)
  d$y <- stats::ave(stats::rnorm(nrow(d)), d$id, FUN = cumsum) + ifelse(d$g > 0 & d$t >= d$g, 1, 0)
  d$y[((d$id * 7 + d$t * 3) %% 11) == 0] <- NA
  d
}

wmp_family_fit <- function(d, family, bytes = NULL, whole = FALSE) {
  if (whole) {
    eligible <- wmp_ns("cp04_pipeline_phase_eligible")
    testthat::local_mocked_bindings(cp04_pipeline_phase_eligible = function(cp04, transfer, pretest, imputation,
        augmented, custom, weighting, strategy) {
      if (identical(weighting, "within")) FALSE else
        eligible(cp04, transfer, pretest, imputation, augmented, custom, weighting, strategy)
    }, .package = "cdid")
  }
  ctl <- list(chunk_rows = 16L)
  if (!is.null(bytes)) ctl$numeric_workspace_bytes <- bytes
  args <- list(yname = "y", tname = "t", idname = "id", gname = "g", data = d, xformla = ~1,
    allow_unbalanced_panel = TRUE, bstrap = FALSE, cband = FALSE, base_period = "universal",
    gmm_weighting = "within", control_group = "nevertreated", control = ctl)
  args[names(family$args)] <- family$args
  if (!is.null(args$selection)) args$selection_xformla <- ~1
  suppressWarnings(do.call(cdid::att_gt, args))
}

testthat::test_that("phase accounting admits every general within family the whole-path envelope refuses", {
  # Workspaces between the smallest admitted workspace of the phase route
  # and that of the whole-path envelope on each design (bytes).
  families <- list(
    pooled_dr = list(panel = "fixed", bytes = 223e6, args = list(est_method = "dr", outcome_pooling = "pooled")),
    imputation_ipw = list(panel = "fixed", bytes = 169e6,
      args = list(est_method = "ipw", selection = "conditional_mtar", selection_est_method = "imputation")),
    reference_dr = list(panel = "fixed", bytes = 473e6, args = list(est_method = "dr", xformla = ~xv)),
    selection_reg = list(panel = "random", bytes = 76e6,
      args = list(est_method = "reg", selection = "conditional_mtar", selection_est_method = "ipw")),
    selection_ipw = list(panel = "random", bytes = 90e6,
      args = list(est_method = "ipw", selection = "conditional_mtar", selection_est_method = "ipw")),
    augmented_reg = list(panel = "random", bytes = 128e6,
      args = list(est_method = "reg", selection = "conditional_mtar", selection_est_method = "dr")))
  for (nm in names(families)) {
    f <- families[[nm]]
    d <- if (identical(f$panel, "fixed")) wmp_panel_fixed() else wmp_panel_random()
    testthat::expect_error(wmp_family_fit(d, f, f$bytes, whole = TRUE), class = "cdid_ipw_owner_bound_error",
      info = nm)
    admitted <- wmp_family_fit(d, f, f$bytes)
    reference <- wmp_family_fit(d, f, whole = TRUE)
    ra <- admitted$chained$resource_accounting
    testthat::expect_false(is.null(ra$memory_phases), info = nm)
    testthat::expect_gt(ra$within_phases$loading_entries, 0)
    testthat::expect_gt(ra$within_phases$covariance_entries, 0)
    testthat::expect_lte(max(unlist(ra$within_phases)), f$bytes / 8)
    for (field in c("att", "se", "inffunc", "V_analytical"))
      testthat::expect_identical(admitted[[field]], reference[[field]], info = paste(nm, field))
    testthat::expect_identical(admitted$chained$available, reference$chained$available, info = nm)
    testthat::expect_identical(admitted$chained$covariance_handling, reference$chained$covariance_handling, info = nm)
  }
})

# The standalone raw IPW backend (incomplete panel, static covariate).
wmp_raw_panel <- function(N = 100L, n_t = 12L) {
  set.seed(7)
  g <- sample(c(0, 4, 7, 10), N, TRUE)
  d <- expand.grid(t = seq_len(n_t), id = seq_len(N)); d$g <- g[d$id]
  d$xs <- round(stats::rnorm(N)[d$id], 3)
  d$y <- stats::ave(stats::rnorm(nrow(d)), d$id, FUN = cumsum) + ifelse(d$g > 0 & d$t >= d$g, 1, 0) + 0.3 * d$xs
  d$y[stats::runif(nrow(d)) < 0.25] <- NA
  d
}

wmp_raw_fit <- function(d, w, bytes = NULL, whole = FALSE) {
  if (whole) testthat::local_mocked_bindings(ipw_raw_phase_eligible = function(...) FALSE, .package = "cdid")
  ctl <- list(chunk_rows = 16L)
  if (!is.null(bytes)) ctl$numeric_workspace_bytes <- bytes
  suppressWarnings(cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d, xformla = ~xs,
    allow_unbalanced_panel = TRUE, bstrap = FALSE, cband = FALSE, base_period = "universal", gmm_weighting = w,
    control_group = "nevertreated", est_method = "ipw", control = ctl))
}

testthat::test_that("raw IPW phase eligibility and deferred terms", {
  eligible <- wmp_ns("ipw_raw_phase_eligible")
  strategy <- list(method = "tikhonov")
  testthat::expect_true(eligible(TRUE, FALSE, FALSE, NULL, "within", strategy))
  testthat::expect_true(eligible(TRUE, FALSE, FALSE, NULL, "optimal", strategy))
  testthat::expect_false(eligible(TRUE, FALSE, FALSE, NULL, "identity", strategy))
  testthat::expect_false(eligible(FALSE, FALSE, FALSE, NULL, "within", strategy))
  testthat::expect_false(eligible(TRUE, TRUE, FALSE, NULL, "within", strategy))
  testthat::expect_false(eligible(TRUE, FALSE, TRUE, NULL, "within", strategy))
  testthat::expect_false(eligible(TRUE, FALSE, FALSE, function(...) NULL, "within", strategy))
  testthat::expect_false(eligible(TRUE, FALSE, FALSE, NULL, "within", NULL))
  defer <- wmp_ns("ipw_raw_phase_defer")
  phase <- list(kind = "cp04_pipeline_memory_phases", whole_path_estimate = 1e6,
    deferred_covariance_entries = 100, deferred_replay_entries = 7)
  out <- defer(phase, 1600, 320)
  testthat::expect_identical(out$whole_path_estimate, 1e6 + 1600 + 320)
  testthat::expect_identical(out$deferred_covariance_entries, 100 + 1600)
  testthat::expect_identical(out$deferred_replay_entries, 7 + 320)
})

testthat::test_that("raw IPW optimal and within are admitted by phase where the whole-path envelope refuses", {
  d <- wmp_raw_panel()
  # Between the smallest admitted workspace of the phase record (74.37 MB on
  # this design, both weightings) and that of the whole-path envelope (76.20 MB).
  bytes <- 75.3e6
  for (w in c("optimal", "within")) {
    testthat::expect_error(wmp_raw_fit(d, w, bytes, whole = TRUE), class = "cdid_ipw_owner_bound_error", info = w)
    admitted <- wmp_raw_fit(d, w, bytes)
    reference <- wmp_raw_fit(d, w, whole = TRUE)
    ra <- admitted$chained$resource_accounting
    testthat::expect_false(is.null(ra$memory_phases), info = w)
    testthat::expect_null(reference$chained$resource_accounting$memory_phases)
    # The retained-basis covariance is deferred, not part of the persistent prefix.
    L <- as.double(admitted$chained$covariance_handling$dimensions$L)
    testthat::expect_gte(ra$memory_phases$deferred_covariance_entries, 16 * L^2)
    testthat::expect_lt(ra$simultaneous_entries, reference$chained$resource_accounting$simultaneous_entries)
    # The simultaneous count falls by exactly the deferred terms (including the
    # retained-basis covariance and replay), less the phase metadata record.
    mp <- ra$memory_phases
    testthat::expect_identical(reference$chained$resource_accounting$simultaneous_entries - ra$simultaneous_entries,
      mp$whole_path_estimate - mp$persistent_entries)
    for (field in c("att", "se", "inffunc", "V_analytical"))
      testthat::expect_identical(admitted[[field]], reference[[field]], info = paste(w, field))
    testthat::expect_identical(admitted$chained$available, reference$chained$available)
    testthat::expect_identical(admitted$chained$covariance_handling, reference$chained$covariance_handling)
  }
  # identity keeps the whole-path envelope, as on every CP04 backend.
  testthat::expect_null(wmp_raw_fit(d, "identity")$chained$resource_accounting$memory_phases)
})

# The raw IPW `within` route checks its working-covariance phase with the raw
# backend, the empirical covariance stream and the construction charge.
testthat::test_that("raw IPW within covariance phase composes every live array", {
  seen <- new.env(parent = emptyenv()); seen$checks <- list()
  check <- wmp_ns("ipw_pipeline_within_phase_check")
  stream_cov <- wmp_ns("cp04_stream_covariance")
  covariance_entries <- wmp_ns("cp04_pipeline_within_covariance_entries")
  testthat::local_mocked_bindings(
    ipw_pipeline_within_phase_check = function(nuis_res, retained, phase_entries, label, registry_prefix = TRUE) {
      total <- check(nuis_res, retained, phase_entries, label, registry_prefix)
      raw <- nuis_res$raw_ipw
      # The prefix holds the raw backend's owned entries and the owner storage it
      # borrows; both are removed once from a retained value that extends it.
      base <- if (registry_prefix) wmp_ns("cp04_composed_retained_entries")(nuis_res, retained) else retained
      seen$checks[[length(seen$checks) + 1L]] <- list(label = label, phase = phase_entries, total = total,
        overlap = if (registry_prefix) wmp_ns("cp04_owner_overlap_entries")(nuis_res) else 0,
        cp04 = nuis_res$cp04, owned = raw$accounting$owned_entries,
        expected = base + raw$accounting$owned_entries + raw$accounting$borrowed_entries + sum(phase_entries),
        max = raw$control$max_matrix_entries)
      total
    },
    cp04_pipeline_within_covariance_entries = function(L, n_a, n_cal, p) {
      seen$formula <- covariance_entries(L, n_a, n_cal, p); seen$formula
    },
    cp04_stream_covariance = function(...) { out <- stream_cov(...); seen$stream <- out; out },
    .package = "cdid")
  fit <- wmp_raw_fit(wmp_raw_panel(), "within")
  ra <- fit$chained$resource_accounting
  testthat::expect_false(is.null(ra$memory_phases))
  cov <- Filter(function(x) identical(x$label, "CP04 within covariance phase simultaneous entries"), seen$checks)
  testthat::expect_length(cov, 1L)
  cc <- cov[[1L]]
  testthat::expect_null(cc$cp04)
  testthat::expect_gt(cc$owned, 0)
  testthat::expect_identical(cc$total, cc$expected)
  testthat::expect_lte(cc$total, cc$max)
  testthat::expect_identical(cc$phase[1L], unname(wmp_ns("cp04_payload")(seen$stream)["entries"]))
  testthat::expect_identical(cc$phase[length(cc$phase)], seen$formula)
  testthat::expect_identical(ra$within_phases$covariance_entries, cc$total)
})
