v59_ns <- function(nm) getFromNamespace(nm, "cdid")
v59_e <- new.env(parent = asNamespace("cdid"))
for (x in parse(file = testthat::test_path("test-augmented-adjoint-replay.R"))) if (is.call(x) && identical(x[[1L]], as.name("<-"))) eval(x, v59_e)

v59_collect <- function(fit, enabled = TRUE) {
  links <- rev(seq_along(fit$links)); B <- fit$control$max_chunk_rows
  plan <- v59_ns("chained_full_if_plan")(fit, links, B)
  on.exit(v59_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
  testthat::expect_type(plan, "list")
  access <- v59_ns("augmented_adjoint_replay_access")
  hits <- misses <- tapped_chunks <- duplicate_models <- step_predictions_saved <- 0L
  testthat::local_mocked_bindings(augmented_adjoint_replay_access = function(...) {
    op <- access(...); if (is.null(op)) return(NULL)
    # Disabled numerical reference retains the original independent model rows.
    if (!enabled) op$treatment_representatives <- op$treatment_ids
    get_old <- op$get; tapped <- FALSE
    op$get <- function(...) {
      out <- get_old(...)
      if (!tapped) for (core in rev(sys.frames())) {
        if (!exists("coarse_operation", core, inherits = FALSE) || !exists("eval_saved_adjoint", core, inherits = FALSE)) next
        native <- get("coarse_operation", core, inherits = FALSE)
        if (is.null(native)) next
        current <- native$contract_current; prepared <- FALSE
        native$contract_current <- function(reader, adj, ids) {
          val <- current(reader, adj, ids)
          if (!prepared && exists("odds_slot", core, inherits = FALSE)) {
            slot <- get("odds_slot", core, inherits = FALSE)
            nt <- length(get("sets", core)$tids); b <- get("B", core); pt <- fit$treatment_recipe$p
            testthat::expect_equal(dim(get("odds_cache", core)), c(b, nt))
            testthat::expect_lte(length(get("odds_cache", core)) + length(get("odds_seen", core)), 2 * b * nt * pt)
            trt <- get("trt_if", core)
            tids <- get("sets", core)$tids
            reps <- get("adjoint_replay", core)$treatment_representatives
            duplicate <- which(tids != reps)
            duplicate_models <<- duplicate_models + length(duplicate)
            for (pos in duplicate) {
              testthat::expect_identical(trt[[tids[pos]]], trt[[reps[pos]]])
              tf <- fit$treatments[[tids[pos]]]
              if (any(get("ci", core) %in% tf$all_cohort_indices))
                step_predictions_saved <<- step_predictions_saved + 1L
            }
            testthat::expect_true(all(vapply(trt[get("sets", core)$tids], function(z)
              identical(names(z), c("direct", "a", "gamma")), logical(1L))))
            assign("odds_slot", function(tid, cohort, ix, coefficients) {
              k <- slot(tid, cohort, ix, coefficients)
              if (k) {
                if (all(get("odds_seen", core)[ix, k])) hits <<- hits + 1L else misses <<- misses + 1L
                if (length(ix) > 1L) testthat::expect_identical(slot(tid, cohort, rev(ix), coefficients), 0L)
                testthat::expect_identical(slot(tid, cohort, ix, coefficients + 1), 0L)
              }
              if (enabled) k else 0L
            }, core)
            prepared <<- TRUE; tapped_chunks <<- tapped_chunks + 1L
          }
          val
        }
        assign("coarse_operation", native, core); tapped <<- TRUE; break
      }
      out
    }
    op
  }, .package = "cdid")
  phi <- v59_e$v57_collect(fit, links, plan)
  list(phi = phi, hits = hits, misses = misses, chunks = tapped_chunks,
    duplicate_models = duplicate_models, step_predictions_saved = step_predictions_saved)
}

testthat::test_that("ordinary treatment odds reuse the exact original link rows and clipped predictor", {
  for (model in c("pair_logit", "independent_period_logit")) {
    a <- v59_e$v57_wide_case(); fit <- v59_e$v57_aug$rap_fit(a, imp = "separate", model = model, mode = "pooled")
    old <- v59_collect(fit, FALSE); new <- v59_collect(fit, TRUE)
    testthat::expect_gt(new$chunks, 0L); testthat::expect_gt(new$hits, 0L); testthat::expect_gt(new$misses, 0L)
    testthat::expect_equal(new$phi, old$phi, tolerance = 1e-12)
    testthat::expect_identical(new$phi == 0, old$phi == 0)
    testthat::expect_equal(crossprod(new$phi), crossprod(old$phi), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(new$phi))), sqrt(diag(crossprod(old$phi))), tolerance = 1e-12)
    testthat::expect_gt(max(abs(new$phi[!a$A, , drop = FALSE])), 1e-4)
    N <- nrow(a$Y); offset <- length(v59_e$v57_aug$rap_point(fit)) - length(fit$links); links <- rev(seq_along(fit$links))
    for (i in c(1L, N)) {
      d <- rep(-1, N); d[i] <- N - 1
      ref <- (v59_e$v57_aug$rap_reference(a, fit, 1 + 1e-6 * d) - v59_e$v57_aug$rap_reference(a, fit, 1 - 1e-6 * d)) / 2e-6
      testthat::expect_equal(new$phi[i, ], ref[offset + links], tolerance = 1e-5)
    }
    detail <- v59_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, details = TRUE)
    testthat::expect_true(all(vapply(detail$details$trt_if[!vapply(detail$details$trt_if,is.null,logical(1L))], function(z)
      identical(names(z), c("direct", "a", "gamma", "score_direct", "score_a")), logical(1L))))
    testthat::expect_identical(colnames(detail$phi), vapply(fit$links[links], `[[`, character(1L), "link_id"))
    testthat::expect_equal(unname(new$phi[1:11, ]), unname(detail$phi), tolerance = 1e-12)
  }
})

testthat::test_that("transparent development instrumentation preserves scientific output", {
  a <- v59_e$v57_wide_case(); fit <- v59_e$v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  links <- rev(seq_along(fit$links)); baseline <- v59_e$v57_collect(fit, links)
  original <- v59_ns("ipw_logit_values"); calls <- 0L
  testthat::local_mocked_bindings(ipw_logit_values = function(...) { calls <<- calls + 1L; original(...) }, .package = "cdid")
  plan <- v59_ns("chained_full_if_plan")(fit, links, fit$control$max_chunk_rows)
  actual <- v59_e$v57_collect(fit, links, plan)
  v59_ns("augmented_adjoint_replay_release")(plan)
  testthat::expect_gt(calls, 0L)
  testthat::expect_equal(actual, baseline, tolerance = 1e-12)
})


testthat::test_that("repeated exact keys avoid real link prediction calls", {
  a <- v59_e$v57_wide_case(); fit <- v59_e$v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  original <- v59_ns("ipw_logit_values"); solve <- v59_ns("cp04_solve"); calls <- solves <- 0L
  # This is transparent DEVELOPMENT instrumentation, not a production
  # authenticity requirement. Both paths call the same installed arithmetic.
  testthat::local_mocked_bindings(ipw_logit_values = function(...) {
    calls <<- calls + 1L; original(...)
  }, cp04_solve = function(...) { solves <<- solves + 1L; solve(...) }, .package = "cdid")
  old <- v59_collect(fit, FALSE); old_calls <- calls; old_solves <- solves
  calls <- solves <- 0L; new <- v59_collect(fit, TRUE); new_calls <- calls; new_solves <- solves
  testthat::expect_gt(new$hits, 0L)
  testthat::expect_lt(new_calls, old_calls)
  testthat::expect_gt(new$duplicate_models, 0L)
  testthat::expect_equal(old_calls - new_calls, new$hits + new$step_predictions_saved, tolerance = 0)
  testthat::expect_equal(old_solves - new_solves, 2L * new$duplicate_models, tolerance = 0)
  testthat::expect_equal(new$phi, old$phi, tolerance = 1e-12)
})


testthat::test_that("producer treatment families retain every original model and exact fitted operand", {
  a <- v59_e$v57_wide_case(); fit <- v59_e$v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  links <- rev(seq_along(fit$links)); before <- fit$treatments
  plan <- v59_ns("chained_full_if_plan")(fit, links, fit$control$max_chunk_rows)
  on.exit(v59_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
  testthat::expect_type(plan, "list")
  testthat::expect_lt(length(unique(plan$treatment_representatives)), length(plan$treatment_ids))
  fields <- c("target_cohort_index", "all_cohort_indices", "coefficients", "H_chol", "a_forcing", "column_names")
  for (pos in seq_along(plan$treatment_ids)) {
    tid <- plan$treatment_ids[pos]; representative <- plan$treatment_representatives[pos]
    testthat::expect_true(representative %in% plan$treatment_ids[seq_len(pos)])
    for (nm in fields) testthat::expect_identical(fit$treatments[[tid]][[nm]], fit$treatments[[representative]][[nm]])
  }
  payload <- v59_ns("cp04_payload")(list(ids = plan$treatment_ids, representatives = plan$treatment_representatives))
  testthat::expect_lte(unname(payload["entries"]), plan$treatment_retained_bound)
  testthat::expect_lte(plan$treatment_retained_bound + plan$treatment_construction_peak, plan$treatment_funded_entries)
  testthat::expect_lte(plan$treatment_funded_entries, 1024 + 256 * (length(plan$treatment_ids) + length(fit$treatments)))
  testthat::expect_identical(fit$treatments, before)
  new <- v59_e$v57_collect(fit, links, plan)
  old <- v59_e$v57_collect(fit, links, NULL)
  testthat::expect_equal(new, old, tolerance = 1e-12)
  testthat::expect_identical(new == 0, old == 0)
})
