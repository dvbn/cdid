# V43 reuses the accepted independent finite-law fixture and raw nuisance
# refits. No other test expressions or fixture outputs are executed/copied.
v43_ns <- function(name) getFromNamespace(name, "cdid")
v43_fixture <- local({
  e <- new.env(parent = environment())
  wanted <- c("drp_ns", "drp_gc", "drp_case", "drp_reference", "drp_fit", "drp_point")
  for (x in parse(file = testthat::test_path("test-dr-outcome-pooling.R"))) {
    if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
        as.character(x[[2L]]) %in% wanted) eval(x, e)
  }
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE)))
  e
})
v43_case <- function(model = "pair_logit", pooling = "separate", universe = "F", zero = FALSE) {
  a <- v43_fixture$drp_case(universe)
  frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
  if (zero) {
    a$reader <- v43_ns("construct_cp04_frame_reader")(frame, seq_len(nrow(a$Y)),
      matrix(numeric(), nrow(a$Y), 0L), a$g, a$w, a$control,
      U = a$U, rho = a$rho, universe = universe,
      justification = "existing finite-law zero-column selection witness", has_intercept = FALSE)
  }
  originals <- a$views
  # Genuine comparison views, distinct consumer identities, three actual
  # models/blocks. Repeated consumers fund the optional index preparation in
  # the old row-independent map slab without adding any nuisance fit.
  a$views <- lapply(seq_len(2100L), function(i) {
    v <- originals[[1L + (i - 1L) %% 3L]]
    id <- paste0("v43-", i, if (i <= 3L) strrep("\u03bb\u6f22", 40L) else "")
    v43_ns("construct_ipw_link_view")(frame, v$pair$s, v$pair$t, v$cohorts$g, v$cohorts$K,
      control = a$control, link_id = id, functional_key = paste0("cmpcat:", id))
  })
  fit <- v43_fixture$drp_fit(a, "conditional_mtar", model, pooling)
  list(a = a, fit = fit)
}
v43_reseal <- function(fit, control = fit$control, treatments = fit$treatments) {
  x <- fit; x["ref"] <- list(NULL); x$control <- control; x$treatments <- treatments
  v43_ns("ipw_owner_seal")(x, "dr_corollary1_backend")
}
v43_rebind <- function(env, name, value) {
  old <- get(name, env, inherits = FALSE)
  unlockBinding(name, env); assign(name, value, env); lockBinding(name, env)
  function() { unlockBinding(name, env); assign(name, old, env); lockBinding(name, env) }
}

testthat::test_that("issued conditional metadata reuses exact full/tail budgets and original full-F arithmetic", {
  cases <- list(list(model = "pair_logit", pooling = "separate", universe = "F"),
    list(model = "independent_period_logit", pooling = "across_cohorts", universe = "F"))
  check <- function(spec) {
    data <- do.call(v43_case, spec); a <- data$a; fit <- data$fit
    N <- nrow(a$Y); links <- rev(seq_along(fit$links)); B <- 11L
    attr(links, "source_notes") <- strrep("\u03bb\u6f22", 40L)
    original_inventory <- v43_ns("dr_conditional_separate_ordinary_inventory")
    original_budget <- v43_ns("cp04_budget")
    original_operation <- v43_ns("selection_prediction_operation")
    inventories <- integer(); guard_labels <- character(); operations <- 0L
    testthat::local_mocked_bindings(
      dr_conditional_separate_ordinary_inventory = function(backend, B, ...) {
        inventories <<- c(inventories, B); original_inventory(backend, B, ...)
      }, cp04_budget = function(ctrl, retained, scratch, work, output, label) {
        guard_labels <<- c(guard_labels, label)
        original_budget(ctrl, retained, scratch, work, output, label)
      }, selection_prediction_operation = function(...) {
        operations <<- operations + 1L; original_operation(...)
      }, .package = "cdid")
    plan <- v43_ns("chained_full_if_plan")(fit, links, B)
    testthat::expect_type(plan, "list")
    testthat::expect_identical(plan$kind, "dr_conditional_replay_preparation")
    testthat::expect_identical(sort(inventories), c(8L, 11L))
    testthat::expect_identical(plan$full, v43_ns("dr_separate_if_envelope")(fit, B, links, FALSE))
    testthat::expect_identical(plan$tail, v43_ns("dr_separate_if_envelope")(fit, 8L, links, FALSE))
    testthat::expect_lte(plan$entries, plan$retained_bound)
    testthat::expect_lte(unname(v43_ns("cp04_payload")(plan)["entries"]), plan$entries)
    testthat::expect_lte(plan$reservation, floor(plan$maps / 2))
    testthat::expect_equal(plan$reservation, plan$retained_bound + plan$construction_peak)
    testthat::expect_length(fit$treatments, 3L)
    testthat::expect_length(fit$outcome_groups$groups, 3L)
    testthat::expect_true(any(rowSums(a$obs) == 0L))
    testthat::expect_true(any(rowSums(a$obs) == 1L))
    inventories <- integer(); guard_labels <- character(); operations <- 0L
    actual <- matrix(0, N, length(links))
    batches <- 0L
    for (first in seq.int(1L, N, by = B)) {
      rows <- rev(seq.int(first, min(N, first + B - 1L)))
      out <- v43_ns("chained_full_if_chunk")(fit, rows, links, "frame", plan)
      actual[rows, ] <- out$phi; batches <- batches + 1L
    }
    testthat::expect_length(inventories, 0L)
    testthat::expect_identical(operations, batches)
    for (label in c("separate DR dependency metadata", "selection IF chunk", "separate DR IF chunk"))
      testthat::expect_identical(sum(guard_labels == label), batches)
    # Both partitions compare whole result/labels/diagnostics, not only Phi.
    for (rows in list(11:1, 96:89)) testthat::expect_identical(
      v43_ns("chained_full_if_chunk")(fit, rows, links, "frame", plan),
      v43_ns("evaluate_dr_corollary1_if_chunk")(fit, rows, links, "frame"))
    representatives <- 1:3
    scalar <- detailed <- matrix(0, N, 3L)
    for (first in seq.int(1L, N, by = B)) {
      rows <- seq.int(first, min(N, first + B - 1L))
      scalar[rows, ] <- v43_ns("evaluate_dr_corollary1_if_chunk")(fit, rows, representatives)$phi
      detailed[rows, ] <- v43_ns("evaluate_dr_corollary1_if_chunk")(fit, rows, representatives, details = TRUE)$phi
    }
    target <- match(representatives, links)
    testthat::expect_equal(actual[, target], scalar, tolerance = 1e-10)
    testthat::expect_equal(actual[, target], detailed, tolerance = 1e-10)
    testthat::expect_equal(crossprod(actual[, target]), crossprod(scalar), tolerance = 1e-10)
    offset <- length(v43_fixture$drp_point(fit)) - 2 * length(fit$links) - 4 * length(fit$treatments)
    for (i in c(1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      derivative <- (v43_fixture$drp_reference(a, fit, 1 + 1e-6 * delta) -
        v43_fixture$drp_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6
      testthat::expect_equal(actual[i, target], derivative[offset + representatives], tolerance = 1e-5)
    }
    # One-row tail has exactly the old one-row record: metadata funding does
    # not depend on B or suppress original tail guards/zero-observer rows.
    one_tail <- v43_ns("chained_full_if_plan")(fit, links, 5L)
    testthat::expect_identical(one_tail$tail_B, 1L)
    testthat::expect_identical(one_tail$tail, v43_ns("dr_separate_if_envelope")(fit, 1L, links, FALSE))
    testthat::expect_identical(v43_ns("chained_full_if_chunk")(fit, N, links, "core", one_tail),
      v43_ns("evaluate_dr_corollary1_if_chunk")(fit, N, links, "core"))
    subset <- links[-seq_len(60L)]
    subset_plan <- v43_ns("chained_full_if_plan")(fit, subset, B)
    testthat::expect_identical(subset_plan$kind, "dr_conditional_replay_preparation")
    testthat::expect_identical(v43_ns("chained_full_if_chunk")(fit, 11:1, subset, "frame", subset_plan),
      v43_ns("evaluate_dr_corollary1_if_chunk")(fit, 11:1, subset))
    # Fresh source operation even after a caller callback changes a sealed
    # selection binding; the prepared integer indexes grant no permission.
    bad <- fit$selection_backend; bad$blocks[[1L]]$alpha[1L] <- Inf
    restore <- v43_rebind(fit$selection_backend$ref, "selection_models", bad)
    error <- tryCatch(v43_ns("chained_full_if_chunk")(fit, 1:11, links, "frame", plan), error = identity)
    restore()
    testthat::expect_s3_class(error, "cdid_cp04_integrity_error")
    during_chunk <- function() {
      restore <- NULL; original_pair <- v43_ns("cp04_pair_rows")
      on.exit(if (!is.null(restore)) restore(), add = TRUE)
      testthat::local_mocked_bindings(cp04_pair_rows = function(...) {
        value <- original_pair(...)
        if (is.null(restore)) restore <<- v43_rebind(fit$selection_backend$ref, "selection_models", bad)
        value
      }, .package = "cdid")
      tryCatch(v43_ns("chained_full_if_chunk")(fit, 1:11, links, "frame", plan), error = identity)
    }
    testthat::expect_s3_class(during_chunk(), "cdid_cp04_integrity_error")
    testthat::expect_error(v43_ns("chained_full_if_chunk")(fit, 1:10, links, "frame", plan),
      class = "cdid_dr_corollary1_integrity_error")
    testthat::expect_error(v43_ns("chained_full_if_chunk")(fit, 1:11, rev(links), "frame", plan),
      class = "cdid_dr_corollary1_integrity_error")
    changed <- plan; changed$support_pos[1L] <- 0L
    testthat::expect_error(v43_ns("chained_full_if_chunk")(fit, 1:11, links, "frame", changed),
      class = "cdid_dr_corollary1_integrity_error")
    other <- v43_reseal(fit)
    testthat::expect_error(v43_ns("chained_full_if_chunk")(other, 1:11, links, "frame", plan),
      class = "cdid_dr_corollary1_integrity_error")
    testthat::expect_error(v43_ns("dr_conditional_replay_upstream_index")(plan,
      fit$selection_backend, fit$up_map, rev(plan$full$uids), FALSE), class = "cdid_cp04_integrity_error")
    noisy <- links; attr(noisy, "source_notes") <- strrep("\u03bb\u6f22", 100000L)
    testthat::expect_null(v43_ns("chained_full_if_plan")(fit, noisy, B))
  }
  for (spec in cases) check(spec)
})

testthat::test_that("unfunded metadata keeps old evaluation, caps, dated and noneligible paths", {
  a <- v43_fixture$drp_case("B"); fit <- v43_fixture$drp_fit(a, "conditional_mtar", pooling = "within_cohort")
  links <- 1:3
  seals <- 0L; original <- v43_ns("ipw_owner_seal")
  testthat::local_mocked_bindings(ipw_owner_seal = function(...) {
    seals <<- seals + 1L; original(...)
  }, .package = "cdid")
  testthat::expect_null(v43_ns("chained_full_if_plan")(fit, links, 11L))
  testthat::expect_identical(seals, 0L)
  testthat::expect_identical(v43_ns("chained_full_if_chunk")(fit, 1:11, links, "frame"),
    v43_ns("evaluate_dr_corollary1_if_chunk")(fit, 1:11, links))
  pre <- v43_ns("dr_separate_if_envelope")(fit, 11L, links, FALSE)
  owned <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  ctrl <- fit$control; ctrl$max_matrix_entries <- owned + pre$scratch
  admitted <- v43_reseal(fit, ctrl)
  testthat::expect_identical(v43_ns("chained_full_if_chunk")(admitted, 1:11, links, "frame"),
    v43_ns("evaluate_dr_corollary1_if_chunk")(fit, 1:11, links))
  ctrl$max_matrix_entries <- ctrl$max_matrix_entries - 1
  refused <- v43_reseal(fit, ctrl)
  testthat::expect_error(v43_ns("chained_full_if_chunk")(refused, 1:11, links, "frame"),
    class = "cdid_ipw_owner_bound_error")
  tr <- fit$treatments; tr[[1L]]$partials$beta[2L, 1L] <- Inf
  bad <- v43_reseal(fit, treatments = tr)
  actual <- tryCatch(v43_ns("chained_full_if_chunk")(bad, 1:11, links, "frame"), error = identity)
  reference <- tryCatch(v43_ns("evaluate_dr_corollary1_if_chunk")(bad, 1:11, links), error = identity)
  testthat::expect_s3_class(actual, "cdid_cp04_numeric_error")
  testthat::expect_identical(class(actual), class(reference))
  testthat::expect_identical(conditionMessage(actual), conditionMessage(reference))
  for (mode in c("pooled", "separate")) {
    dated <- v43_fixture$drp_case(dated = TRUE, extra_reference = TRUE)
    b <- v43_fixture$drp_fit(dated, "conditional_mtar", "independent_period_logit", "across_cohorts", mode)
    plan <- v43_ns("chained_full_if_plan")(b, links, 11L)
    testthat::expect_null(plan)
    testthat::expect_identical(v43_ns("chained_full_if_chunk")(b, 1:11, links, "frame", plan),
      v43_ns("evaluate_dr_corollary1_if_chunk")(b, 1:11, links))
    testthat::expect_identical(v43_ns("evaluate_dr_corollary1_if_chunk")(b, 1:11, links, details = TRUE, plan = plan),
      v43_ns("evaluate_dr_corollary1_if_chunk")(b, 1:11, links, details = TRUE))
  }
})


testthat::test_that("theta-free selection still funds both retained leaf-budget cell maps", {
  data <- v43_case(zero = TRUE); fit <- data$fit
  links <- seq_along(fit$links)
  plan <- v43_ns("chained_full_if_plan")(fit, links, 11L)
  testthat::expect_identical(plan$kind, "dr_conditional_replay_preparation")
  testthat::expect_identical(plan$active_count, 0)
  testthat::expect_identical(plan$cell_count, 0)
  testthat::expect_length(plan$upstream_index$active_cells, 0L)
  testthat::expect_gt(length(plan$selection_full$active_cells), 0L)
  testthat::expect_identical(plan$selection_full$active_cells, plan$selection_tail$active_cells)
  testthat::expect_lte(unname(v43_ns("cp04_payload")(plan)["entries"]), plan$retained_bound)
  testthat::expect_lte(plan$reservation, floor(plan$maps / 2))
  for (rows in list(11:1, 96:89)) testthat::expect_identical(
    v43_ns("chained_full_if_chunk")(fit, rows, links, "frame", plan),
    v43_ns("evaluate_dr_corollary1_if_chunk")(fit, rows, links))
})
