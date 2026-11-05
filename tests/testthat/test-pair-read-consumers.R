# Production constructors and independent raw-data oracles are loaded as
# definitions only; no stored expected values or old test statements change.
v51_ns <- function(name) getFromNamespace(name, "cdid")
v51_helpers <- function(file, wanted) {
  e <- new.env(parent = asNamespace("cdid"))
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE)))
  e
}
v51_dr <- v51_helpers("test-dr-outcome-pooling.R", c("drp_ns", "drp_gc", "drp_case", "drp_fit", "drp_point", "drp_reference"))
v51_or <- v51_helpers("test-or-selection-outcome-pooling.R", c("osop_ns", "osop_gc", "osop_case", "osop_fit"))
v51_rec <- v51_helpers("test-reconstruction-score-cache.R", c("rsc_ns", "rsc_case", "rsc_fit"))
v51_rcp <- v51_helpers("test-reconstructed-outcome-pooling.R",
  c("rcp_ns", "rcp_gc", "rcp_case", "rcp_fit", "rcp_point", "rcp_reference"))
v51_trace <- function(run, enabled = TRUE) {
  factory <- v51_ns("cp04_pair_read_operation"); date <- v51_ns("cp04_date_rows")
  operations <- list(); read_only <- outcomes <- dates <- 0L
  testthat::local_mocked_bindings(cp04_pair_read_operation = function(...) {
    op <- if (enabled) factory(...) else NULL
    if (!is.null(op)) {
      pair <- op$pair
      op$pair <- function(r_sub, descriptor, read_outcomes = FALSE) {
        if (read_outcomes) outcomes <<- outcomes + 1L else read_only <<- read_only + 1L
        pair(r_sub, descriptor, read_outcomes)
      }
      operations[[length(operations) + 1L]] <<- op
    }
    op
  }, cp04_date_rows = function(...) { dates <<- dates + 1L; date(...) }, .package = "cdid")
  value <- run()
  list(value = value, operations = operations, counts = lapply(operations, function(x) x$counts()),
    read_only = read_only, outcomes = outcomes, dates = dates)
}
v51_full <- function(fit, fn, enabled = TRUE, details = FALSE, links = 3:1) {
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE); N <- reader$state$N
  v51_trace(function() {
    phi <- matrix(0, N, length(links)); budgets <- list()
    for (first in seq.int(1L, N, by = fit$control$max_chunk_rows)) {
      rows <- rev(seq.int(first, min(N, first + fit$control$max_chunk_rows - 1L)))
      out <- v51_ns(fn)(fit, rows, links, details = details)
      phi[rows, ] <- out$phi; budgets[[length(budgets) + 1L]] <- out$diagnostics
    }
    list(phi = phi, budgets = budgets)
  }, enabled)
}
v51_equal <- function(fast, old, positive = TRUE) {
  testthat::expect_identical(fast$value, old$value)
  testthat::expect_identical(crossprod(fast$value$phi), crossprod(old$value$phi))
  testthat::expect_identical(sqrt(diag(crossprod(fast$value$phi))), sqrt(diag(crossprod(old$value$phi))))
  if (positive) {
    testthat::expect_gt(length(fast$operations), 0L)
    testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["queries"]), numeric(1L))), 0)
    testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["cache_hits"]), numeric(1L))), 0)
    testthat::expect_gt(fast$outcomes, 0L)
    for (j in seq_along(fast$operations)) {
      op <- fast$operations[[j]]
      testthat::expect_lte(fast$counts[[j]]["searches"], op$chunk_rows * op$n_dates)
      testthat::expect_lte(op$reservation, op$funding_entries)
      testthat::expect_error(op$pair(1L, list()), class = "cdid_cp04_integrity_error")
    }
  }
}

testthat::test_that("general IPW and OR share outcome and selection pair positions on every original row", {
  for (model in c("pair_logit", "independent_period_logit")) {
    a <- v51_dr$drp_case()
    fit <- v51_ns("fit_cp04_selection_models")(a$reader, a$views, a$control, v51_dr$drp_gc,
      selection_model = model, selection_pooling = "within_cohort")
    fast <- v51_full(fit, "evaluate_cp04_if_chunk"); old <- v51_full(fit, "evaluate_cp04_if_chunk", FALSE)
    v51_equal(fast, old)
    if (model == "pair_logit") testthat::expect_gt(fast$read_only, 0L) else {
      testthat::expect_identical(fast$read_only, 0L)
      testthat::expect_gt(fast$dates, 0L)
      testthat::expect_identical(fast$dates, old$dates)
    }
    testthat::expect_gt(max(abs(fast$value$phi[!a$A, , drop = FALSE])), 1e-4)
  }
  for (mode in c("pooled", "separate")) {
    a <- v51_or$osop_case()
    fit <- v51_or$osop_fit(a, pooling = "within_cohort", mode = mode, general = TRUE)
    fast <- v51_full(fit, "evaluate_or_corollary1_if_chunk")
    v51_equal(fast, v51_full(fit, "evaluate_or_corollary1_if_chunk", FALSE))
    testthat::expect_gt(fast$read_only, 0L)
    detailed <- v51_full(fit, "evaluate_or_corollary1_if_chunk", details = TRUE)
    testthat::expect_length(detailed$operations, 0L)
    testthat::expect_equal(fast$value$phi, detailed$value$phi, tolerance = 1e-10)
  }
  a <- v51_or$osop_case(); fit <- v51_or$osop_fit(a, mode = "pooled")
  testthat::expect_false(isTRUE(fit$uses_selection_models))
  v51_equal(v51_full(fit, "evaluate_or_corollary1_if_chunk"),
    v51_full(fit, "evaluate_or_corollary1_if_chunk", FALSE))
})

testthat::test_that("pooled and separate conditional DR keep full-F nuisance derivatives and date singleton responses", {
  for (mode in c("pooled", "separate")) for (model in c("pair_logit", "independent_period_logit")) {
    a <- v51_dr$drp_case()
    fit <- v51_dr$drp_fit(a, "conditional_mtar", model, "within_cohort", mode)
    fast <- v51_full(fit, "evaluate_dr_if_chunk")
    old <- v51_full(fit, "evaluate_dr_if_chunk", FALSE)
    v51_equal(fast, old)
    rich <- v51_full(fit, "evaluate_dr_if_chunk", details = TRUE)
    testthat::expect_length(rich$operations, 0L)
    testthat::expect_equal(fast$value$phi, rich$value$phi, tolerance = 1e-10)
    if (model == "pair_logit") {
      testthat::expect_gt(fast$read_only, 0L)
      offset <- length(v51_dr$drp_point(fit)) - 2 * length(fit$links) - 4 * length(fit$treatments)
      N <- nrow(a$Y)
      for (i in c(1L, N)) {
        delta <- rep(-1, N); delta[i] <- N - 1
        derivative <- (v51_dr$drp_reference(a, fit, 1 + 1e-6 * delta) -
          v51_dr$drp_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6
        testthat::expect_equal(fast$value$phi[i, ], derivative[offset + 3:1], tolerance = 1e-5)
      }
    } else {
      testthat::expect_identical(fast$read_only, 0L)
      testthat::expect_gt(fast$dates, 0L)
      testthat::expect_identical(fast$dates, old$dates)
    }
  }
})

testthat::test_that("all pure methods and augmented pooling modes use the same owned ordinary pair reader", {
  specs <- list(c("imputation", "ipw", "pooled"), c("imputation", "reg", "pooled"),
    c("imputation", "reg", "separate"), c("imputation", "dr", "pooled"),
    c("imputation", "dr", "separate"), c("augmented", "dr", "pooled"), c("augmented", "dr", "separate"))
  for (spec in specs) {
    a <- v51_rec$rsc_case(); fit <- v51_rec$rsc_fit(a, spec[1L], spec[2L], spec[3L])
    fn <- if (spec[1L] == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk"
    fast <- v51_full(fit, fn); old <- v51_full(fit, fn, FALSE)
    small_reg <- identical(spec, c("imputation", "reg", "separate"))
    v51_equal(fast, old, positive = !small_reg)
    if (small_reg) {
      # Actual root inventory: this original48-row fixture has no treatment
      # rows. Its140B slice cannot fund4120+73B at either B32 or tail16.
      # Keep its complete science checks and require truthful safe decline;
      # the existing larger fixture below supplies the positive same-route gate.
      testthat::expect_length(fast$operations, 0L)
      testthat::expect_length(v51_ns("imputation_if_needed_sets")(fit, 3:1, "reg")$tids, 0L)
      for (B in c(32L, 16L)) testthat::expect_lt(140 * B, 4096 + 64 * B + 8 * 3 + 3 * B * 3)
    }
    rich <- v51_full(fit, fn, details = TRUE)
    testthat::expect_length(rich$operations, 0L)
    testthat::expect_equal(fast$value$phi, rich$value$phi, tolerance = 1e-10)
    st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
    testthat::expect_gt(max(abs(fast$value$phi[st$a_index == 0L, , drop = FALSE])), 1e-4)
  }
})

testthat::test_that("existing larger finite fixture positively funds the same pure reg separate route", {
  a <- v51_rcp$rcp_case()
  testthat::expect_gte(nrow(a$Y), 62L)
  # Explicit private test batch only: retain the fixture's data, weights,
  # formulas, other controls and every numeric capacity/tolerance unchanged.
  a$control$max_chunk_rows <- 64L
  fit <- v51_rcp$rcp_fit(a, "reg", "within_cohort", "separate")
  testthat::expect_identical(fit$est_method, "reg")
  testthat::expect_identical(fit$outcome_pooling, "separate")
  fast <- v51_full(fit, "evaluate_imputation_if_chunk")
  old <- v51_full(fit, "evaluate_imputation_if_chunk", FALSE)
  # Retain every original query/hit/outcome/final-source assertion here.
  v51_equal(fast, old)
  testthat::expect_true(any(vapply(fast$operations, function(op) op$chunk_rows >= 62L, logical(1L))))
  rich <- v51_full(fit, "evaluate_imputation_if_chunk", details = TRUE)
  testthat::expect_length(rich$operations, 0L)
  testthat::expect_equal(fast$value$phi, rich$value$phi, tolerance = 1e-10)
  testthat::expect_gt(max(abs(fast$value$phi[!a$A, , drop = FALSE])), 1e-4)
  offset <- length(v51_rcp$rcp_point(fit)) - length(fit$links); N <- nrow(a$Y)
  for (i in c(1L, N)) {
    delta <- rep(-1, N); delta[i] <- N - 1
    derivative <- (v51_rcp$rcp_reference(a, fit, 1 + 1e-6 * delta) -
      v51_rcp$rcp_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6
    testthat::expect_equal(fast$value$phi[i, ], derivative[offset + 3:1], tolerance = 1e-5)
  }
})

testthat::test_that("capacity-only decline leaves original rows/budgets and excluded readers unchanged", {
  a <- v51_rec$rsc_case(); fit <- v51_rec$rsc_fit(a, "imputation", "dr", "separate")
  # The old one-row partition admits its old buffers but cannot fund the full
  # position-operation fixed header. It must stay a truthful old-reader call.
  declined <- v51_trace(function() v51_ns("evaluate_imputation_if_chunk")(fit, 1L, 1L))
  testthat::expect_length(declined$operations, 0L)
  testthat::expect_identical(declined$value,
    v51_trace(function() v51_ns("evaluate_imputation_if_chunk")(fit, 1L, 1L), FALSE)$value)
  # Explicit exclusion is before any new source permission or pair query.
  pre <- v51_ns("imputation_if_budget")(fit, 32L, 1:3)
  reader <- a$reader
  testthat::expect_null(v51_ns("cp04_pair_consumer_operation")(fit, reader, 1:32, 1:3, pre, "reconstruction", TRUE))
  retained <- fit; retained$pair_outcomes <- list()
  testthat::expect_null(v51_ns("cp04_pair_consumer_operation")(retained, reader, 1:32, 1:3, pre, "reconstruction", FALSE))
  layout <- fit; layout$layout <- list()
  testthat::expect_null(v51_ns("cp04_pair_consumer_operation")(layout, reader, 1:32, 1:3, pre, "reconstruction", FALSE))
  # No operation changes the exact existing numeric admission boundary.
  owned <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  cap <- fit; cap["ref"] <- list(NULL); cap$control$max_matrix_entries <- owned + pre$scratch
  cap <- v51_ns("ipw_owner_seal")(cap, "imputation_backend")
  v51_equal(v51_full(cap, "evaluate_imputation_if_chunk"),
    v51_full(cap, "evaluate_imputation_if_chunk", FALSE))
  low <- cap; low["ref"] <- list(NULL); low$control$max_matrix_entries <- low$control$max_matrix_entries - 1
  low <- v51_ns("ipw_owner_seal")(low, "imputation_backend")
  testthat::expect_error(v51_ns("evaluate_imputation_if_chunk")(low, 1:32, 1:3), class = "cdid_ipw_owner_bound_error")
})

testthat::test_that("actual public conditional DR uses shared positions and keeps inference/status exactly", {
  e <- new.env(parent = asNamespace("cdid"))
  e$sp_raw <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_baseline", "raw.csv"))
  e$sp_raw$cohort <- as.numeric(replace(e$sp_raw$cohort, e$sp_raw$cohort == "NT", "Inf"))
  e$sp_trans <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_corollary1", "transport.csv"))
  e$sp_trans$cohort <- as.numeric(replace(e$sp_trans$cohort, e$sp_trans$cohort == "NT", "Inf"))
  for (x in parse(file = testthat::test_path("test-selection-public.R"))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name("sp_public_fixture"))) eval(x, e)
  args <- e$sp_public_fixture("dr", selection_pooling = "within_cohort")$args
  args$outcome_pooling <- "pooled"; args$control <- list(chunk_rows = 23L)
  fast <- v51_trace(function() do.call(cdid::att_gt, args))
  old <- v51_trace(function() do.call(cdid::att_gt, args), FALSE)
  testthat::expect_gt(length(fast$operations), 0L)
  testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["queries"]), numeric(1L))), 0)
  testthat::expect_gt(sum(vapply(fast$counts, function(x) unname(x["cache_hits"]), numeric(1L))), 0)
  testthat::expect_gt(fast$read_only, 0L)
  for (nm in c("att", "inffunc", "V_analytical", "se")) testthat::expect_identical(fast$value[[nm]], old$value[[nm]])
  for (nm in c("available", "target_meta", "support_status", "support_reason", "covariance_handling"))
    testthat::expect_identical(fast$value$chained[[nm]], old$value$chained[[nm]])
})

testthat::test_that("owned consumer successful boundaries still validate original sources and abort after refusal", {
  a <- v51_rec$rsc_case(); fit <- v51_rec$rsc_fit(a, "imputation", "dr", "pooled")
  factory <- v51_ns("cp04_pair_read_operation"); ops <- list(); aborted <- 0L
  original <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
  restore <- function() {
    ref <- a$reader$frame_ref; unlockBinding("frame_owner", ref)
    assign("frame_owner", original, ref); lockBinding("frame_owner", ref)
  }
  on.exit(restore(), add = TRUE)
  testthat::local_mocked_bindings(cp04_pair_read_operation = function(...) {
    op <- factory(...)
    if (!is.null(op)) {
      finish <- op$finish; abort <- op$abort
      op$finish <- function() {
        ref <- a$reader$frame_ref; other <- original; other$unit_weights_all[1L] <- other$unit_weights_all[1L] * 2
        unlockBinding("frame_owner", ref); assign("frame_owner", other, ref); lockBinding("frame_owner", ref)
        finish()
      }
      op$abort <- function() { aborted <<- aborted + 1L; abort() }
      ops[[length(ops) + 1L]] <<- op
    }
    op
  }, .package = "cdid")
  testthat::expect_error(v51_ns("evaluate_imputation_if_chunk")(fit, 1:32, 1:3), class = "cdid_cp04_integrity_error")
  testthat::expect_length(ops, 1L)
  testthat::expect_gt(aborted, 0L)
  restore()
})
