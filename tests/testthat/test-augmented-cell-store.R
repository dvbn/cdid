v65_ns <- function(nm) getFromNamespace(nm, "cdid")
v65_e <- new.env(parent = asNamespace("cdid"))
for (x in parse(file = testthat::test_path("test-augmented-adjoint-replay.R")))
  if (is.call(x) && identical(x[[1L]], as.name("<-"))) eval(x, v65_e)

v65_reservation <- function(fit, links, caller = 0, control = fit$control) {
  covariance <- v65_ns("cp04_stream_plan")(fit, links, caller)
  point <- v65_ns("cp04_preflight_links")(fit, links, caller)
  d5 <- list(ledger_if_entries = 0, basis_if_entries = 0, ledger_enclosure_entries = 0, basis_enclosure_entries = 0)
  inv <- v65_ns("cp04_retained_if_phase_inventory")("covariance", fit, covariance, point,
    caller, caller, d5, length(links), length(fit$links))
  store <- v65_ns("cp04_retained_if_phase_reserve")(fit, links, "covariance", covariance$chunk_rows,
    covariance, retained = inv$retained, live = inv$live, work = 0, control = control)
  list(store = store, covariance = covariance, inv = inv)
}

v65_collect <- function(fit, links, B, plan = NULL) {
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  out <- matrix(0, N, length(links))
  for (first in seq.int(1L, N, by = B)) {
    rows <- seq.int(first, min(N, first + B - 1L))
    out[rows, ] <- v65_ns("evaluate_augmented_if_chunk")(fit, rows, links, plan = plan)$phi
  }
  out
}

v65_score_reference <- function(fit, rows, cid) {
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE); st <- reader$state
  X <- v65_ns("augmented_backend_X_F")(fit); cf <- fit$cells[[cid]]; gf <- fit$groups[[cf$group_id]]
  ix <- which(st$cohort_index[rows] == cf$cohort_index)
  if (!length(ix)) return(NULL)
  hit <- v65_ns("cp04_pair_rows")(reader, rows[ix], cf$pair, TRUE); obs <- which(hit$S == 1)
  if (!length(obs)) return(NULL)
  rr <- ix[obs]
  d <- if (gf$p) cf$L * X[rows[rr], gf$kept_cols_F, drop = FALSE] else matrix(0, length(rr), 0L)
  j <- if (gf$p) cf$alpha + as.vector(d %*% gf$beta) else rep(cf$alpha, length(rr))
  list(rr = rr, d_rr = d, s_rr = st$w[rows[rr]] * (hit$dy[obs] - j))
}

testthat::test_that("consumer-funded ragged cells replace repeated construction while full-F science stays", {
  run_model <- function(model) {
    a <- v65_e$v57_wide_case(); fit <- v65_e$v57_aug$rap_fit(a, imp = "separate", model = model, mode = "pooled")
    links <- rev(seq_along(fit$links)); reservation <- v65_reservation(fit, links)
    store <- reservation$store; testthat::expect_true(is.environment(store))
    B <- reservation$covariance$chunk_rows
    plan <- v65_ns("chained_full_if_plan")(fit, links, B, phase_store = store)
    on.exit(v65_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
    on.exit(v65_ns("cp04_retained_if_release")(store), add = TRUE)
    cell <- plan$store$cell_replay
    testthat::expect_type(cell, "list")
    info <- cell$accounting
    testthat::expect_lte(info$actual_entries, info$retained_entries)
    testthat::expect_lte(info$admission_total, info$cap)
    testthat::expect_lte(info$build_rows, fit$control$max_chunk_rows)
    testthat::expect_equal(info$admission_total,
      store$admission_total + 2 * info$retained_entries + info$construction_peak, tolerance = 0)
    N <- nrow(a$Y); pre <- v65_ns("augmented_if_budget")(fit, B, links, FALSE)
    testthat::expect_lt(info$build_batches, ceiling(N / B) * length(pre$cids))
    wanted <- sort(unique(c(pre$cids, pre$zids)))
    expected_cohort <- sum(vapply(fit$cells[wanted], function(cf) {
      st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
      st$cohort_end[cf$cohort_index] - st$cohort_start[cf$cohort_index] + 1L
    }, numeric(1L)))
    testthat::expect_equal(info$cohort_rows, expected_cohort, tolerance = 0)
    for (first in unique(c(1L, 1L + B, 1L + B * floor((N - 1L) / B)))) {
      rows <- seq.int(first, min(N, first + B - 1L))
      testthat::expect_true(cell$matches(rows))
      for (cid in pre$cids) testthat::expect_equal(cell$score(cid, rows),
        v65_score_reference(fit, rows, cid), tolerance = 1e-12)
    }
    testthat::expect_false(cell$matches(rev(seq_len(B))))
    reference <- v65_collect(fit, links, B)
    worker <- v65_ns("augmented_cell_rows"); calls <- 0L
    testthat::local_mocked_bindings(augmented_cell_rows = function(...) {
      calls <<- calls + 1L; worker(...)
    }, .package = "cdid")
    actual <- v65_collect(fit, links, B, plan)
    testthat::expect_identical(calls, 0L)
    testthat::expect_equal(actual, reference, tolerance = 1e-12)
    testthat::expect_identical(actual == 0, reference == 0)
    testthat::expect_equal(crossprod(actual), crossprod(reference), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(actual))), sqrt(diag(crossprod(reference))), tolerance = 1e-12)
    testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
    # Unsupported private row shape selects old callbacks before arithmetic.
    calls <- 0L; rows <- rev(seq_len(B))
    arbitrary <- v65_ns("evaluate_augmented_if_chunk")(fit, rows, links, plan = plan)
    old <- v65_ns("evaluate_augmented_if_chunk")(fit, rows, links)
    testthat::expect_gt(calls, 0L)
    testthat::expect_equal(arbitrary$phi, old$phi, tolerance = 1e-12)
    offset <- length(v65_e$v57_aug$rap_point(fit)) - length(fit$links)
    for (i in c(1L, N)) {
      d <- rep(-1, N); d[i] <- N - 1
      ref <- (v65_e$v57_aug$rap_reference(a, fit, 1 + 1e-6 * d) -
        v65_e$v57_aug$rap_reference(a, fit, 1 - 1e-6 * d)) / 2e-6
      testthat::expect_equal(actual[i, ], ref[offset + links], tolerance = 1e-5)
    }
    v65_ns("augmented_adjoint_replay_release")(plan)
    testthat::expect_error(cell$score(pre$cids[1L], seq_len(B)), class = "cdid_augmented_integrity_error")
    testthat::expect_null(plan$store$cell_replay)
  }
  for (model in c("pair_logit", "independent_period_logit")) run_model(model)
})

testthat::test_that("insufficient genuine consumer headroom declines cell building before outcomes", {
  a <- v65_e$v57_wide_case(); fit <- v65_e$v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  links <- rev(seq_along(fit$links)); reservation <- v65_reservation(fit, links); store <- reservation$store
  plan <- v65_ns("chained_full_if_plan")(fit, links, reservation$covariance$chunk_rows, phase_store = store)
  cell <- plan$store$cell_replay; testthat::expect_type(cell, "list")
  cap <- cell$accounting$admission_total - 1
  v65_ns("augmented_adjoint_replay_release")(plan); v65_ns("cp04_retained_if_release")(store)
  control <- v65_ns("ipw_owner_control")(fit$control$max_chunk_rows, cap, cap,
    fit$control$max_work_ops, fit$control$max_output_entries)
  small <- v65_reservation(fit, links, control = control); testthat::expect_true(is.environment(small$store))
  reads <- 0L; read <- v65_ns("cp04_pair_rows")
  testthat::local_mocked_bindings(cp04_pair_rows = function(...) { reads <<- reads + 1L; read(...) }, .package = "cdid")
  declined <- v65_ns("chained_full_if_plan")(fit, links, small$covariance$chunk_rows, phase_store = small$store)
  on.exit(v65_ns("augmented_adjoint_replay_release")(declined), add = TRUE)
  on.exit(v65_ns("cp04_retained_if_release")(small$store), add = TRUE)
  testthat::expect_type(declined, "list"); testthat::expect_null(declined$store$cell_replay)
  testthat::expect_identical(reads, 0L)
  testthat::expect_false(exists("producer_cell_accounting", small$store, inherits = FALSE))
  actual <- v65_collect(fit, links, small$covariance$chunk_rows, declined)
  reference <- v65_collect(fit, links, small$covariance$chunk_rows)
  testthat::expect_equal(actual, reference, tolerance = 1e-12)
})


testthat::test_that("ragged primitive retains p0 matrices and the original probability floor", {
  # Existing deterministic finite-law generator, only its supported recipe
  # role is changed. No RNG or hand-edited outcome/weight fixture.
  original_case <- v65_e$v57_aug$rap_case
  v65_e$v57_aug$rap_case <- function() original_case(slopes = FALSE)
  on.exit(v65_e$v57_aug$rap_case <- original_case, add = TRUE)
  a <- v65_e$v57_wide_case()
  fit <- v65_e$v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  # Same frozen-fit low-intercept branch witness as the existing replay test.
  selection <- fit$selection_backend
  for (j in seq_along(selection$blocks)) {
    bl <- selection$blocks[[j]]
    if (identical(bl$branch, "cp04_separate_pair")) {
      bl$coefficients[] <- 0; bl$coefficients[1L] <- -30
    } else { bl$alpha[] <- -30; bl$beta[] <- 0 }
    selection$blocks[[j]] <- bl
  }
  selection["ref"] <- list(NULL)
  selection <- v65_ns("ipw_owner_seal")(selection, "selection_models")
  helper <- new.env(parent = asNamespace("cdid")); helper$v45_ns <- v65_ns
  for (x in parse(file = testthat::test_path("test-augmented-replay-operation.R")))
    if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
        identical(x[[2L]], as.name("v45_reseal"))) eval(x, helper)
  fit <- helper$v45_reseal(fit, selection)
  links <- rev(seq_along(fit$links)); reservation <- v65_reservation(fit, links)
  on.exit(v65_ns("cp04_retained_if_release")(reservation$store), add = TRUE)
  B <- reservation$covariance$chunk_rows
  pre <- v65_ns("augmented_if_budget")(fit, B, links, FALSE)
  cell <- v65_ns("augmented_cell_replay_store")(fit, pre, reservation$store)
  testthat::expect_type(cell, "list"); on.exit(cell$release(), add = TRUE)
  testthat::expect_true(all(vapply(fit$groups, `[[`, integer(1L), "p") == 0L))
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  X <- v65_ns("augmented_backend_X_F")(fit)
  pair_reader <- function(rows, pair, read_outcomes = TRUE) v65_ns("cp04_pair_rows")(reader, rows, pair, read_outcomes)
  positive <- empty <- FALSE
  for (first in seq.int(1L, nrow(a$Y), by = B)) {
    rows <- seq.int(first, min(nrow(a$Y), first + B - 1L))
    for (cid in pre$cids) {
      score <- cell$score(cid, rows); ref <- v65_score_reference(fit, rows, cid)
      testthat::expect_equal(score, ref, tolerance = 1e-12)
      if (is.null(score)) empty <- TRUE else {
        positive <- TRUE; testthat::expect_equal(ncol(score$d_rr), 0L, tolerance = 0)
      }
    }
  }
  testthat::expect_true(positive); testthat::expect_true(empty)
  cid <- pre$zids[1L]; cf <- fit$cells[[cid]]; gf <- fit$groups[[cf$group_id]]
  for (first in seq.int(1L, nrow(a$Y), by = B)) {
    rows <- seq.int(first, min(nrow(a$Y), first + B - 1L))
    ix <- which(reader$state$cohort_index[rows] == cf$cohort_index)
    if (!length(ix)) next
    reference <- v65_ns("augmented_cell_rows")(selection, fit$up_map, X, pair_reader, cf, gf, rows[ix], FALSE)
    testthat::expect_equal(cell$Z(cid, rows), list(ix = ix, Z = reference$Z), tolerance = 1e-12)
    if (length(reference$q)) {
      testthat::expect_true(any(reference$q < 1e-6))
      testthat::expect_identical(v65_ns("ipw_clip_floor")(min(reference$q)), 1e-6)
    }
  }
})
