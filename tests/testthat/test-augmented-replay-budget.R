v63_ns <- function(nm) getFromNamespace(nm, "cdid")
v63_e <- new.env(parent = asNamespace("cdid"))
for (x in parse(file = testthat::test_path("test-augmented-adjoint-replay.R")))
  if (is.call(x) && identical(x[[1L]], as.name("<-"))) eval(x, v63_e)

testthat::test_that("issued augmented full chunks reuse exact preflight and still check actual budgets", {
  a <- v63_e$v57_wide_case()
  fit <- v63_e$v57_aug$rap_fit(a, imp = "separate", model = "pair_logit", mode = "pooled")
  links <- rev(seq_along(fit$links)); B <- fit$control$max_chunk_rows; N <- nrow(a$Y)
  plan <- v63_ns("chained_full_if_plan")(fit, links, B)
  on.exit(v63_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
  testthat::expect_type(plan$replay_preflight, "list")
  fresh <- v63_ns("augmented_if_budget")(fit, B, links, FALSE)
  testthat::expect_identical(plan$replay_preflight, fresh)
  testthat::expect_identical(v63_ns("augmented_adjoint_replay_preflight")(plan, fit, links, B), fresh)
  testthat::expect_lte(unname(v63_ns("cp04_payload")(fresh)["entries"]), plan$replay_preflight_retained_bound)
  testthat::expect_lte(plan$retained_bound + plan$construction_peak, plan$funded_entries)
  testthat::expect_null(v63_ns("augmented_adjoint_replay_preflight")(plan, fit, rev(links), B))
  testthat::expect_null(v63_ns("augmented_adjoint_replay_preflight")(plan, fit, links, N %% B))
  testthat::expect_null(v63_ns("augmented_adjoint_replay_preflight")(NULL, fit, links, B))
  budget <- v63_ns("augmented_if_budget"); arithmetic <- v63_ns("cp04_budget")
  inventories <- checked <- 0L
  testthat::local_mocked_bindings(augmented_if_budget = function(...) {
    inventories <<- inventories + 1L; budget(...)
  }, cp04_budget = function(ctrl, retained, scratch, work, output, label) {
    if (identical(label, "augmented observation-correction IF chunk")) checked <<- checked + 1L
    arithmetic(ctrl, retained, scratch, work, output, label)
  }, .package = "cdid")
  old <- v63_e$v57_collect(fit, links, NULL)
  old_inventories <- inventories; old_checked <- checked
  inventories <- checked <- 0L
  actual <- v63_e$v57_collect(fit, links, plan)
  testthat::expect_identical(old_inventories, as.integer(ceiling(N / B)))
  testthat::expect_identical(inventories, as.integer(N %% B != 0L))
  testthat::expect_identical(old_checked, as.integer(ceiling(N / B)))
  testthat::expect_identical(checked, as.integer(ceiling(N / B)))
  testthat::expect_equal(actual, old, tolerance = 1e-12)
  testthat::expect_identical(actual == 0, old == 0)
  testthat::expect_equal(crossprod(actual), crossprod(old), tolerance = 1e-12)
  testthat::expect_equal(sqrt(diag(crossprod(actual))), sqrt(diag(crossprod(old))), tolerance = 1e-12)
  testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
  # Details still constructs and checks its own different full preflight.
  inventories <- checked <- 0L
  detail <- v63_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, details = TRUE, plan = plan)
  testthat::expect_identical(inventories, 1L); testthat::expect_identical(checked, 1L)
  testthat::expect_equal(unname(detail$phi), unname(actual[1:11, , drop = FALSE]), tolerance = 1e-12)
})

testthat::test_that("cached preflight never replaces real plan, model or control integrity", {
  a <- v63_e$v57_wide_case(); fit <- v63_e$v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  links <- rev(seq_along(fit$links)); B <- fit$control$max_chunk_rows
  plan <- v63_ns("chained_full_if_plan")(fit, links, B)
  on.exit(v63_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
  changed <- plan; changed$replay_preflight$output <- 0
  testthat::expect_null(v63_ns("augmented_adjoint_replay_preflight")(changed, fit, links, B))
  testthat::expect_error(v63_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, plan = changed),
    class = "cdid_augmented_integrity_error")
  changed <- fit; changed$control$max_output_entries <- 1
  testthat::expect_error(v63_ns("evaluate_augmented_if_chunk")(changed, 1:11, links, plan = plan),
    class = "cdid_augmented_integrity_error")
  v63_ns("augmented_adjoint_replay_release")(plan)
  testthat::expect_null(v63_ns("augmented_adjoint_replay_preflight")(plan, fit, links, B))
  # Released optional metadata falls back to the original numerical route.
  actual <- v63_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, plan = plan)
  reference <- v63_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, plan = NULL)
  testthat::expect_equal(actual$phi, reference$phi, tolerance = 1e-12)
})


testthat::test_that("decorated link resources always use the original full inventory", {
  a <- v63_e$v57_wide_case(); fit <- v63_e$v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  links <- rev(seq_along(fit$links)); B <- fit$control$max_chunk_rows
  plan <- v63_ns("chained_full_if_plan")(fit, links, B)
  on.exit(v63_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
  named <- stats::setNames(links, paste0("actual decorated link ", seq_along(links)))
  testthat::expect_null(v63_ns("augmented_adjoint_replay_preflight")(plan, fit, named, B))
  budget <- v63_ns("augmented_if_budget"); seen <- list()
  testthat::local_mocked_bindings(augmented_if_budget = function(backend, B, links, details = FALSE) {
    seen[[length(seen) + 1L]] <<- links
    budget(backend, B, links, details)
  }, .package = "cdid")
  actual <- v63_ns("evaluate_augmented_if_chunk")(fit, 1:11, named, plan = plan)
  testthat::expect_length(seen, 1L); testthat::expect_identical(seen[[1L]], named)
  reference <- v63_ns("evaluate_augmented_if_chunk")(fit, 1:11, named, plan = NULL)
  testthat::expect_equal(actual$phi, reference$phi, tolerance = 1e-12)
})
