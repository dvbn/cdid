amr_ns <- function(name) getFromNamespace(name, "cdid")
amr_helpers <- function(file, wanted) {
  e <- new.env(parent = asNamespace("cdid"))
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE))); e
}
amr_aug <- amr_helpers("test-augmented-outcome-pooling.R",
  c("rap_ns", "rap_gc", "rap_case", "rap_fit", "rap_point", "rap_reference"))
amr_v57 <- amr_helpers("test-augmented-adjoint-replay.R",
  c("v57_ns", "v57_collect"))

amr_wide_case <- function(intercept = TRUE, slopes = TRUE) {
  a <- amr_aug$rap_case(intercept = intercept, slopes = slopes); original <- nrow(a$Y)
  take <- rep(seq_len(original), 7L)
  for (nm in c("X", "V", "X_full", "Q", "Y", "obs")) a[[nm]] <- a[[nm]][take, , drop = FALSE]
  for (nm in c("w", "U", "A", "rho")) a[[nm]] <- a[[nm]][take]
  labels <- c(seq(11, by = 12, length.out = 20L), Inf)
  a$g <- rep(labels, each = original / 3L)
  N <- length(take); id <- seq_len(N); T <- length(a$calendar)
  long <- data.frame(id = rep(id, each = T), time = rep(a$calendar, N),
    y = as.vector(t(a$Y)), g = rep(a$g, each = T))
  panel <- amr_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", a$control)
  Xa <- a$X_full; Xa[!a$A, ] <- NA_real_
  frame <- amr_ns("construct_ipw_frame_owner")(panel, id, Xa, a$w, a$control,
    frame_g = a$g, design_on_A_only = TRUE)
  a$reader <- amr_ns("construct_cp04_frame_reader")(frame, id, a$Q, a$g, a$w, a$control,
    U = a$U, rho = a$rho, universe = "F", justification = "Existing finite-law rows copied to a wider cohort catalogue")
  pairs <- utils::combn(a$calendar, 2L); a$views <- list()
  for (g in labels[is.finite(labels)]) for (j in seq_len(ncol(pairs))) {
    k <- length(a$views) + 1L
    a$views[[k]] <- amr_ns("construct_ipw_link_view")(frame, pairs[1L, j], pairs[2L, j], g, Inf,
      control = a$control, link_id = paste0("amr-existing-grid-", k),
      functional_key = paste0("cmpcat:amr-existing-grid-", k))
  }
  a
}

amr_reservation <- function(fit, links, caller = 0, control = fit$control) {
  covariance <- amr_ns("cp04_stream_plan")(fit, links, caller)
  point <- amr_ns("cp04_preflight_links")(fit, links, caller)
  d5 <- list(ledger_if_entries = 0, basis_if_entries = 0, ledger_enclosure_entries = 0, basis_enclosure_entries = 0)
  inv <- amr_ns("cp04_retained_if_phase_inventory")("covariance", fit, covariance, point,
    caller, caller, d5, length(links), length(fit$links))
  store <- amr_ns("cp04_retained_if_phase_reserve")(fit, links, "covariance", covariance$chunk_rows,
    covariance, retained = inv$retained, live = inv$live, work = 0, control = control)
  list(store = store, covariance = covariance, inv = inv)
}

amr_canonical_collect <- function(fit, links, B, plan = NULL) {
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  out <- matrix(0, N, length(links))
  for (first in seq.int(1L, N, by = B)) {
    rows <- seq.int(first, min(N, first + B - 1L))
    out[rows, ] <- amr_ns("evaluate_augmented_if_chunk")(fit, rows, links, plan = plan)$phi
  }
  out
}

amr_counted_collect <- function(fit, links, B, plan) {
  access <- amr_ns("augmented_adjoint_replay_access")
  sparse_ids <- amr_ns("reconstruction_sparse_adjoint_ids")
  cell_worker <- amr_ns("augmented_cell_rows")
  logit_worker <- amr_ns("ipw_logit_values")
  puts <- hits <- T_contracts <- C_contracts <- cell_calls <- logit_calls <- 0L
  tap_coarse <- function() {
    for (e in rev(sys.frames())) {
      if (!exists("eval_saved_adjoint", e, inherits = FALSE) ||
          !exists("coarse_operation", e, inherits = FALSE)) next
      actual <- get("coarse_operation", e, inherits = FALSE)
      if (is.null(actual) || isTRUE(attr(actual, "amr_tapped"))) next
      contract <- actual$contract; core <- e
      actual$contract <- function(view, state, adjoint, ids) {
        out <- contract(view, state, adjoint, ids)
        if (!is.null(out) && exists("lf", core, inherits = FALSE)) {
          lf <- get("lf", core, inherits = FALSE)
          if (identical(adjoint, lf$adjoint_T)) T_contracts <<- T_contracts + 1L
          if (identical(adjoint, lf$adjoint_C)) C_contracts <<- C_contracts + 1L
        }
        out
      }
      current_contract <- actual$contract_current
      if (!is.null(current_contract)) actual$contract_current <- function(reader, adjoint, ids) {
        out <- current_contract(reader, adjoint, ids)
        if (!is.null(out) && exists("lf", core, inherits = FALSE)) {
          lf <- get("lf", core, inherits = FALSE)
          if (identical(adjoint, lf$adjoint_T)) T_contracts <<- T_contracts + 1L
          if (identical(adjoint, lf$adjoint_C)) C_contracts <<- C_contracts + 1L
        }
        out
      }
      attr(actual, "amr_tapped") <- TRUE
      assign("coarse_operation", actual, e)
      break
    }
  }
  testthat::local_mocked_bindings(
    augmented_adjoint_replay_access = function(...) {
      op <- access(...); if (is.null(op)) return(NULL)
      old_get <- op$get; old_put <- op$put
      op$get <- function(...) {
        x <- old_get(...); if (!is.null(x)) hits <<- hits + 1L
        tap_coarse()
        x
      }
      op$put <- function(...) { puts <<- puts + 1L; old_put(...) }
      op
    },
    reconstruction_sparse_adjoint_ids = function(...) {
      tap_coarse()
      sparse_ids(...)
    },
    augmented_cell_rows = function(...) {
      cell_calls <<- cell_calls + 1L
      cell_worker(...)
    },
    ipw_logit_values = function(...) {
      logit_calls <<- logit_calls + 1L
      logit_worker(...)
    },
    .package = "cdid"
  )
  phi <- amr_canonical_collect(fit, links, B, plan)
  list(phi = phi, puts = puts, hits = hits, T_contracts = T_contracts, C_contracts = C_contracts,
    cell_calls = cell_calls, logit_calls = logit_calls)
}

testthat::test_that("augmented IPW activates cell store, treatment/odds reuse, and native coarse contraction without OR", {
  for (model in c("pair_logit", "independent_period_logit")) {
    for (imp in c("separate", "within_cohort", "across_cohorts")) {
      a <- amr_wide_case()
      fit <- amr_aug$rap_fit(a, imp = imp, sel = imp, model = model, mode = "separate", est_method = "ipw")
      testthat::expect_identical(fit$est_method, "ipw")
      testthat::expect_identical(fit$outcome_pooling, "inactive")
      testthat::expect_null(fit$shared_or)
      testthat::expect_identical(length(fit$or_blocks), 0L)
      links <- rev(seq_along(fit$links))
      res <- amr_reservation(fit, links); store <- res$store; B <- res$covariance$chunk_rows
      plan <- amr_ns("chained_full_if_plan")(fit, links, B, phase_store = store)
      on.exit(amr_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
      on.exit(amr_ns("cp04_retained_if_release")(store), add = TRUE)
      testthat::expect_type(plan, "list")
      testthat::expect_identical(plan$kind, "augmented_adjoint_replay_plan")
      testthat::expect_null(plan$shared)
      testthat::expect_identical(plan$bids, integer(0L))
      testthat::expect_identical(plan$q, 0L)
      testthat::expect_identical(plan$alpha, integer(0L))
      testthat::expect_identical(length(plan$store$entries), 0L)
      testthat::expect_true(plan$retained_bound + plan$construction_peak <= plan$funded_entries)
      testthat::expect_true(plan$treatment_retained_bound + plan$treatment_construction_peak <= plan$treatment_funded_entries)
      testthat::expect_lt(length(unique(plan$treatment_representatives)), length(plan$treatment_ids))
      testthat::expect_type(plan$replay_preflight, "list")
      testthat::expect_type(plan$store$cell_replay, "list")
      testthat::expect_lte(plan$store$cell_replay$accounting$actual_entries, plan$store$cell_replay$accounting$retained_entries)
      testthat::expect_lte(plan$store$cell_replay$accounting$admission_total, plan$store$cell_replay$accounting$cap)

      logit_ref_calls <- 0L; logit_worker <- amr_ns("ipw_logit_values")
      reference <- testthat::with_mocked_bindings(
        amr_canonical_collect(fit, links, B, NULL),
        ipw_logit_values = function(...) { logit_ref_calls <<- logit_ref_calls + 1L; logit_worker(...) },
        .package = "cdid"
      )
      counted <- amr_counted_collect(fit, links, B, plan)
      actual <- counted$phi
      testthat::expect_identical(counted$puts, 0L)
      testthat::expect_identical(counted$hits, 0L)
      testthat::expect_identical(counted$cell_calls, 0L)
      testthat::expect_gt(counted$T_contracts, 0L)
      testthat::expect_gt(counted$C_contracts, 0L)
      testthat::expect_lt(counted$logit_calls, logit_ref_calls)
      testthat::expect_equal(actual, reference, tolerance = 1e-12)
      testthat::expect_identical(actual == 0, reference == 0)
      testthat::expect_equal(crossprod(actual), crossprod(reference), tolerance = 1e-12)
      testthat::expect_equal(sqrt(diag(crossprod(actual))), sqrt(diag(crossprod(reference))), tolerance = 1e-12)
      testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)

      if (identical(imp, "separate")) {
        N <- nrow(a$Y); offset <- length(amr_aug$rap_point(fit)) - length(fit$links)
        for (i in c(1L, N)) {
          d <- rep(-1, N); d[i] <- N - 1
          deriv <- (amr_aug$rap_reference(a, fit, 1 + 1e-6 * d) -
            amr_aug$rap_reference(a, fit, 1 - 1e-6 * d)) / 2e-6
          testthat::expect_equal(actual[i, ], deriv[offset + links], tolerance = 1e-5)
        }
      }
      amr_ns("augmented_adjoint_replay_release")(plan)
      amr_ns("cp04_retained_if_release")(store)
    }
  }
})

testthat::test_that("augmented OR and DR reuse cell store and preflight across both outcome poolings without cross-block beta keys", {
  specs <- list(
    list(method = "reg", mode = "separate", imp = "separate", sel = "separate", model = "pair_logit"),
    list(method = "reg", mode = "separate", imp = "within_cohort", sel = "across_cohorts", model = "independent_period_logit"),
    list(method = "reg", mode = "pooled", imp = "within_cohort", sel = "within_cohort", model = "pair_logit"),
    list(method = "reg", mode = "pooled", imp = "across_cohorts", sel = "separate", model = "independent_period_logit"),
    list(method = "dr", mode = "separate", imp = "separate", sel = "separate", model = "pair_logit"),
    list(method = "dr", mode = "separate", imp = "across_cohorts", sel = "within_cohort", model = "independent_period_logit"),
    list(method = "dr", mode = "pooled", imp = "within_cohort", sel = "across_cohorts", model = "pair_logit"),
    list(method = "dr", mode = "pooled", imp = "across_cohorts", sel = "across_cohorts", model = "independent_period_logit")
  )
  for (spec in specs) {
    a <- amr_wide_case()
    fit <- amr_aug$rap_fit(a, imp = spec$imp, sel = spec$sel, model = spec$model,
      mode = spec$mode, est_method = spec$method)
    links <- rev(seq_along(fit$links))
    res <- amr_reservation(fit, links); store <- res$store; B <- res$covariance$chunk_rows
    plan <- amr_ns("chained_full_if_plan")(fit, links, B, phase_store = store)
    on.exit(amr_ns("augmented_adjoint_replay_release")(plan), add = TRUE)
    on.exit(amr_ns("cp04_retained_if_release")(store), add = TRUE)
    testthat::expect_type(plan, "list")
    testthat::expect_true(plan$retained_bound + plan$construction_peak <= plan$funded_entries)
    testthat::expect_type(plan$replay_preflight, "list")
    testthat::expect_type(plan$store$cell_replay, "list")
    if (identical(spec$method, "reg")) {
      testthat::expect_null(fit$treatments)
      testthat::expect_identical(plan$treatment_ids, integer(0L))
      testthat::expect_identical(plan$treatment_representatives, integer(0L))
    } else {
      testthat::expect_lt(length(unique(plan$treatment_representatives)), length(plan$treatment_ids))
    }
    if (identical(spec$mode, "separate")) {
      testthat::expect_identical(plan$q, 0L)
      testthat::expect_identical(plan$alpha, integer(0L))
      testthat::expect_identical(length(plan$store$entries), 0L)
      testthat::expect_true(isTRUE(plan$replay_preflight$stage1_forward$enabled))
    } else {
      testthat::expect_identical(plan$q, fit$shared_or$p)
      testthat::expect_gt(length(plan$store$entries), 0L)
    }

    reference <- amr_canonical_collect(fit, links, B, NULL)
    counted <- amr_counted_collect(fit, links, B, plan)
    actual <- counted$phi
    testthat::expect_identical(counted$cell_calls, 0L)
    if (identical(spec$mode, "separate")) {
      testthat::expect_identical(counted$puts, 0L)
      testthat::expect_identical(counted$hits, 0L)
      testthat::expect_null(plan$store$seed)
    } else {
      testthat::expect_identical(counted$puts, as.integer(length(plan$store$entries)))
      testthat::expect_gt(counted$hits, counted$puts)
      testthat::expect_gt(counted$T_contracts, 0L)
      if (identical(spec$method, "dr")) testthat::expect_gt(counted$C_contracts, 0L) else
        testthat::expect_identical(counted$C_contracts, 0L)
    }
    testthat::expect_equal(actual, reference, tolerance = 1e-12)
    testthat::expect_identical(actual == 0, reference == 0)
    testthat::expect_equal(crossprod(actual), crossprod(reference), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(actual))), sqrt(diag(crossprod(reference))), tolerance = 1e-12)
    testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)

    if (identical(spec$imp, "separate")) {
      N <- nrow(a$Y); offset <- length(amr_aug$rap_point(fit)) - length(fit$links)
      for (i in c(1L, N)) {
        d <- rep(-1, N); d[i] <- N - 1
        deriv <- (amr_aug$rap_reference(a, fit, 1 + 1e-6 * d) -
          amr_aug$rap_reference(a, fit, 1 - 1e-6 * d)) / 2e-6
        testthat::expect_equal(actual[i, ], deriv[offset + links], tolerance = 1e-5)
      }
    }
    amr_ns("augmented_adjoint_replay_release")(plan)
    amr_ns("cp04_retained_if_release")(store)
  }
})

testthat::test_that("intercept-only, no-intercept, dated recipes, permutations, memory declines, and covariance/SE agree", {
  for (basis in list(c(FALSE, TRUE), c(TRUE, FALSE))) {
    a <- amr_wide_case(intercept = basis[1L], slopes = basis[2L])
    for (cfg in list(list(m = "ipw", mode = "separate"), list(m = "reg", mode = "separate"),
        list(m = "reg", mode = "pooled"), list(m = "dr", mode = "separate"), list(m = "dr", mode = "pooled"))) {
      fit <- amr_aug$rap_fit(a, imp = "separate", mode = cfg$mode, est_method = cfg$m)
      links <- seq_along(fit$links)
      res <- amr_reservation(fit, links); store <- res$store; B <- res$covariance$chunk_rows
      plan <- amr_ns("chained_full_if_plan")(fit, links, B, phase_store = store)
      testthat::expect_type(plan, "list")
      testthat::expect_type(plan$store$cell_replay, "list")
      if (identical(cfg$mode, "pooled") && !identical(cfg$m, "ipw")) {
        if (!basis[1L]) testthat::expect_identical(plan$alpha, integer(0L))
        if (!basis[2L]) testthat::expect_identical(plan$q, 0L)
      } else {
        testthat::expect_identical(plan$q, 0L)
        testthat::expect_identical(plan$alpha, integer(0L))
      }
      actual <- amr_canonical_collect(fit, links, B, plan)
      reference <- amr_canonical_collect(fit, links, B, NULL)
      testthat::expect_equal(actual, reference, tolerance = 1e-12)
      amr_ns("augmented_adjoint_replay_release")(plan)
      amr_ns("cp04_retained_if_release")(store)
    }
  }

  # Dated recipes decline common replay plan across all methods and poolings.
  dated <- amr_aug$rap_case(dated = TRUE, route_kind = "augmented")
  for (cfg in list(list(m = "ipw", mode = "separate"), list(m = "reg", mode = "separate"),
      list(m = "dr", mode = "separate"), list(m = "dr", mode = "pooled"))) {
    dfit <- amr_aug$rap_fit(dated, imp = "separate", mode = cfg$mode, est_method = cfg$m)
    testthat::expect_null(amr_ns("chained_full_if_plan")(dfit, 1:3, dfit$control$max_chunk_rows))
    testthat::expect_equal(amr_aug$rap_point(dfit), amr_aug$rap_reference(dated, dfit), tolerance = 2e-8)
  }

  # Permutations, non-canonical row fallback, constrained memory decline, and stream_augmented_covariance agreement.
  a <- amr_wide_case()
  for (cfg in list(list(m = "ipw", mode = "separate"), list(m = "reg", mode = "separate"),
      list(m = "dr", mode = "separate"), list(m = "dr", mode = "pooled"))) {
    fit <- amr_aug$rap_fit(a, imp = "separate", mode = cfg$mode, est_method = cfg$m)
    links <- rev(seq_along(fit$links))
    res <- amr_reservation(fit, links); store <- res$store; B <- res$covariance$chunk_rows
    plan <- amr_ns("chained_full_if_plan")(fit, links, B, phase_store = store)
    cell <- plan$store$cell_replay
    testthat::expect_type(cell, "list")
    cap <- cell$accounting$admission_total - 1
    # Non-canonical reversed rows decline cell_replay and match un-planned evaluation.
    calls <- 0L; worker <- amr_ns("augmented_cell_rows")
    rev_rows <- rev(seq_len(B))
    rev_actual <- testthat::with_mocked_bindings(
      amr_ns("evaluate_augmented_if_chunk")(fit, rev_rows, links, plan = plan),
      augmented_cell_rows = function(...) { calls <<- calls + 1L; worker(...) },
      .package = "cdid"
    )
    rev_ref <- amr_ns("evaluate_augmented_if_chunk")(fit, rev_rows, links)
    testthat::expect_gt(calls, 0L)
    testthat::expect_identical(rev_actual$id, rev_ref$id)
    testthat::expect_equal(rev_actual$phi, rev_ref$phi, tolerance = 1e-12)
    amr_ns("augmented_adjoint_replay_release")(plan)
    amr_ns("cp04_retained_if_release")(store)

    # Constrained consumer headroom declines cell_replay before pair reads while keeping plan active.
    ctrl <- amr_ns("ipw_owner_control")(fit$control$max_chunk_rows, cap, cap,
      fit$control$max_work_ops, fit$control$max_output_entries)
    small <- amr_reservation(fit, links, control = ctrl)
    reads <- 0L; read_fn <- amr_ns("cp04_pair_rows")
    declined <- testthat::with_mocked_bindings(
      amr_ns("chained_full_if_plan")(fit, links, small$covariance$chunk_rows, phase_store = small$store),
      cp04_pair_rows = function(...) { reads <<- reads + 1L; read_fn(...) },
      .package = "cdid"
    )
    testthat::expect_type(declined, "list")
    testthat::expect_null(declined$store$cell_replay)
    testthat::expect_identical(reads, 0L)
    testthat::expect_equal(amr_canonical_collect(fit, links, small$covariance$chunk_rows, declined),
      amr_canonical_collect(fit, links, small$covariance$chunk_rows, NULL), tolerance = 1e-12)
    amr_ns("augmented_adjoint_replay_release")(declined)
    amr_ns("cp04_retained_if_release")(small$store)

    # Exact matrix-entry budget boundary check.
    pre <- amr_ns("augmented_if_budget")(fit, B, links, FALSE)
    need <- sum(c(fit$accounting$owned_entries, fit$accounting$borrowed_entries, pre$scratch))
    below <- fit; below$control$max_matrix_entries <- need - 1
    below <- amr_ns("ipw_owner_seal")(below, "augmented_backend")
    old_err <- tryCatch(amr_ns("augmented_if_budget")(below, B, links, FALSE), error = identity)
    new_err <- tryCatch(amr_ns("chained_full_if_plan")(below, links, B), error = identity)
    testthat::expect_s3_class(old_err, "error")
    testthat::expect_identical(class(new_err), class(old_err))
    testthat::expect_identical(conditionMessage(new_err), conditionMessage(old_err))
    exact <- fit; exact$control$max_matrix_entries <- need
    exact <- amr_ns("ipw_owner_seal")(exact, "augmented_backend")
    exact_plan <- amr_ns("chained_full_if_plan")(exact, links, B)
    testthat::expect_type(exact_plan, "list")
    amr_ns("augmented_adjoint_replay_release")(exact_plan)

    # Full covariance stream agreement (covariance, SE, status/retained record).
    cov_reservation <- amr_reservation(fit, links)
    bound_fit <- amr_ns("cp04_retained_if_bind")(fit, cov_reservation$store)
    cov_res <- amr_ns("cp04_stream_covariance")(bound_fit, links)
    rel <- amr_ns("cp04_retained_if_release")(cov_reservation$store)
    cov_direct <- amr_ns("stream_augmented_covariance")(fit, links)
    ref_core <- amr_ns("evaluate_augmented_if_chunk")(fit, seq_len(B), links, scale = "core")
    full_core <- matrix(0, nrow(a$Y), length(links))
    for (first in seq.int(1L, nrow(a$Y), by = B)) {
      rr <- seq.int(first, min(nrow(a$Y), first + B - 1L))
      full_core[rr, ] <- amr_ns("evaluate_augmented_if_chunk")(fit, rr, links, scale = "core")$phi
    }
    n_A <- ref_core$normalization_count
    centered_core <- sweep(full_core, 2L, colSums(full_core) / nrow(a$Y), FUN = "-")
    expected_cov <- crossprod(centered_core) / n_A
    testthat::expect_equal(unname(cov_res$estimates),
      vapply(fit$links[links], `[[`, numeric(1L), "mu"), tolerance = 1e-12)
    testthat::expect_equal(cov_res$covariance, cov_direct$covariance, tolerance = 1e-12)
    testthat::expect_equal(cov_res$mean_phi, cov_direct$mean_phi, tolerance = 1e-12)
    testthat::expect_equal(unname(cov_res$covariance), unname(expected_cov), tolerance = 1e-12)
    testthat::expect_equal(unname(sqrt(diag(cov_res$covariance))),
      unname(sqrt(diag(expected_cov))), tolerance = 1e-12)
    testthat::expect_identical(rel$status, "retained")
    testthat::expect_identical(rel$producer_passes, 1L)
    testthat::expect_identical(rel$slice_passes, 1L)
    testthat::expect_type(rel$producer_cell_store, "list")
  }
})
