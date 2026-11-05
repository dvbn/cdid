v57_ns <- function(name) getFromNamespace(name, "cdid")
v57_helpers <- function(file, wanted) {
  e <- new.env(parent = asNamespace("cdid"))
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE))); e
}
v57_aug <- v57_helpers("test-augmented-outcome-pooling.R",
  c("rap_ns", "rap_gc", "rap_case", "rap_fit", "rap_point", "rap_reference"))
v57_old <- v57_helpers("test-reconstruction-coarse-pullback.R", c("v53_ns", "v53_scalar", "v53_run", "v53_native"))

# Enlarge ONLY the existing deterministic finite grid's cohort catalogue.
# Every outcome/design/weight/observation row is copied from rap_case(), with
# new original-unit IDs and cohort labels. No RNG, redraw or fitted output.
v57_wide_case <- function() {
  a <- v57_aug$rap_case(); original <- nrow(a$Y)
  take <- rep(seq_len(original), 7L)
  for (nm in c("X", "V", "X_full", "Q", "Y", "obs")) a[[nm]] <- a[[nm]][take, , drop = FALSE]
  for (nm in c("w", "U", "A", "rho")) a[[nm]] <- a[[nm]][take]
  labels <- c(seq(11, by = 12, length.out = 20L), Inf)
  a$g <- rep(labels, each = original / 3L)
  N <- length(take); id <- seq_len(N); T <- length(a$calendar)
  long <- data.frame(id = rep(id, each = T), time = rep(a$calendar, N),
    y = as.vector(t(a$Y)), g = rep(a$g, each = T))
  panel <- v57_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", a$control)
  Xa <- a$X_full; Xa[!a$A, ] <- NA_real_
  frame <- v57_ns("construct_ipw_frame_owner")(panel, id, Xa, a$w, a$control,
    frame_g = a$g, design_on_A_only = TRUE)
  a$reader <- v57_ns("construct_cp04_frame_reader")(frame, id, a$Q, a$g, a$w, a$control,
    U = a$U, rho = a$rho, universe = "F", justification = "Existing finite-law rows copied to a wider cohort catalogue")
  pairs <- utils::combn(a$calendar, 2L); a$views <- list()
  for (g in labels[is.finite(labels)]) for (j in seq_len(ncol(pairs))) {
    k <- length(a$views) + 1L
    a$views[[k]] <- v57_ns("construct_ipw_link_view")(frame, pairs[1L, j], pairs[2L, j], g, Inf,
      control = a$control, link_id = paste0("v57-existing-grid-", k),
      functional_key = paste0("cmpcat:v57-existing-grid-", k))
  }
  a
}

v57_collect <- function(fit, links, plan = NULL) {
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  B <- fit$control$max_chunk_rows
  out <- matrix(0, N, length(links))
  for (first in seq.int(1L, N, by = B)) {
    rows <- rev(seq.int(first, min(N, first + B - 1L)))
    out[rows, ] <- v57_ns("evaluate_augmented_if_chunk")(fit, rows, links, plan = plan)$phi
  }
  out
}

# Observe methods without changing any namespace numerical/source helper.
v57_counted_collect <- function(fit, links, plan) {
  access <- v57_ns("augmented_adjoint_replay_access"); puts <- hits <- T_contracts <- C_contracts <- 0L
  testthat::local_mocked_bindings(augmented_adjoint_replay_access = function(...) {
    op <- access(...); if (is.null(op)) return(NULL)
    old_get <- op$get; old_put <- op$put; tapped <- FALSE
    op$get <- function(...) {
      x <- old_get(...); if (!is.null(x)) hits <<- hits + 1L
      # Observe the already-issued core-local method. No namespace/source or
      # numerical helper is replaced, and the original method decides use.
      if (!tapped) for (e in rev(sys.frames())) {
        if (!exists("eval_saved_adjoint", e, inherits = FALSE) ||
            !exists("coarse_operation", e, inherits = FALSE)) next
        actual <- get("coarse_operation", e, inherits = FALSE)
        if (is.null(actual)) next
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
        assign("coarse_operation", actual, e); tapped <<- TRUE
        break
      }
      x
    }
    op$put <- function(...) { puts <<- puts + 1L; old_put(...) }
    op
  }, .package = "cdid")
  phi <- v57_collect(fit, links, plan)
  list(phi = phi, puts = puts, hits = hits, T_contracts = T_contracts, C_contracts = C_contracts)
}

# Rebuild only a fixed loading with the literal OLD R group/scalar-solve loop;
# all-null score records remove row data without changing any group operation.
v57_original_beta <- function(fit, plan, qq) {
  cells <- fit$cells; groups <- fit$groups; nc <- length(cells)
  lam <- numeric(nc); nu <- lapply(cells, function(cf) numeric(groups[[cf$group_id]]$p))
  for (bid in plan$bids) {
    bf <- fit$shared_or$blocks[[bid]]
    for (k in seq_along(bf$control_cell_ids)) {
      cid <- bf$control_cell_ids[k]
      if (fit$outcome_recipe$intercept) lam[cid] <- lam[cid] + bf$j_force_beta_alpha[qq, k]
      nu[[cid]] <- nu[[cid]] + bf$j_force_beta_beta[[k]][qq, ]
    }
  }
  ids <- unique(unlist(lapply(fit$shared_or$blocks[plan$bids], `[[`, "control_cell_ids")))
  active <- unique(vapply(cells[ids], `[[`, integer(1L), "group_id"))
  pos <- integer(nc); pos[] <- seq_len(nc)
  state <- list(cells = cells, groups = groups, N = plan$N, B = 1L,
    intercept = fit$outcome_recipe$intercept)
  v57_old$v53_scalar(list(cache = vector("list", nc), position = pos), state, lam, nu,
    active, as.integer(ids), vector("list", length(groups)))
}

testthat::test_that("one augmented producer reuses only exact original successful OR adjoints", {
  for (model in c("pair_logit", "independent_period_logit")) {
    a <- v57_wide_case(); fit <- v57_aug$rap_fit(a, imp = "separate", sel = "separate", model = model, mode = "pooled")
    links <- rev(seq_along(fit$links)); B <- fit$control$max_chunk_rows
    plan <- v57_ns("chained_full_if_plan")(fit, links, B)
    testthat::expect_type(plan, "list")
    testthat::expect_identical(plan$kind, "augmented_adjoint_replay_plan")
    testthat::expect_true(plan$retained_bound + plan$construction_peak <= plan$funded_entries)
    testthat::expect_true(plan$min_B == nrow(a$Y) %% B)
    testthat::expect_true(all(vapply(plan$store$entries, is.null, logical(1L))))
    testthat::expect_null(plan$store$seed)
    expected <- v57_collect(fit, links)
    counted <- v57_counted_collect(fit, links, plan)
    actual <- counted$phi; puts <- counted$puts; hits <- counted$hits
    testthat::expect_equal(actual, expected, tolerance = 1e-12)
    testthat::expect_equal(crossprod(actual), crossprod(expected), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(actual))), sqrt(diag(crossprod(expected))), tolerance = 1e-12)
    testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
    testthat::expect_identical(puts, as.integer(length(plan$store$entries)))
    testthat::expect_gt(hits, puts)
    testthat::expect_gt(counted$T_contracts, 0L)
    testthat::expect_gt(counted$C_contracts, 0L)
    testthat::expect_false(any(vapply(plan$store$entries, is.null, logical(1L))))
    testthat::expect_type(plan$store$seed, "list")
    testthat::expect_true(v57_ns("ipw_owner_seal_ok")(plan$store$seed, "augmented_adjoint_seed"))
    original <- v57_original_beta(fit, plan, 1L)
    testthat::expect_equal(plan$store$entries[[1L]]$adjoint, original$adjoint, tolerance = 1e-12)
    testthat::expect_equal(plan$store$seed$templates, original$templates, tolerance = 1e-12)
    # Measure the actual new owned payload; fitted group/cell source values are
    # existing borrowed owners, so replace ONLY those alias values with NULL
    # while keeping every actual list slot/name and current numeric record.
    seed_owned <- length(plan$store$seed$templates) + 256
    for (record in plan$store$seed$templates) if (!is.null(record)) {
      owned <- record; owned["group"] <- list(NULL)
      owned["cells"] <- list(rep(list(NULL), length(record$cells)))
      seed_owned <- seed_owned + unname(v57_ns("cp04_payload")(owned)["entries"])
    }
    entries_owned <- sum(vapply(plan$store$entries, function(x)
      unname(v57_ns("cp04_payload")(x)["entries"]), numeric(1L)))
    testthat::expect_lte(seed_owned + entries_owned + 4096 + 32 * length(plan$store$entries), plan$retained_bound)
    for (j in seq_along(plan$alpha)) {
      adj <- plan$store$entries[[plan$q + j]]$adjoint
      testthat::expect_lte(unname(v57_ns("cp04_payload")(adj)["entries"]), plan$entry_bounds[plan$q + j])
      if (identical(adj$kind, "compact_reconstruction_adjoint")) {
        cids <- fit$shared_or$blocks[[plan$alpha[j]]]$control_cell_ids
        testthat::expect_identical(adj$cell_ids, sort(unique(as.integer(cids))))
        testthat::expect_identical(adj$group_ids,
          sort(unique(vapply(fit$cells[cids], `[[`, integer(1L), "group_id"))))
      }
    }
    # Existing independent frequency refits move EVERY original observation,
    # including zero/singletons, and refit every nuisance at the original law.
    N <- nrow(a$Y); offset <- length(v57_aug$rap_point(fit)) - length(fit$links)
    for (i in c(1L, N)) {
      d <- rep(-1, N); d[i] <- N - 1
      deriv <- (v57_aug$rap_reference(a, fit, 1 + 1e-6 * d) -
        v57_aug$rap_reference(a, fit, 1 - 1e-6 * d)) / 2e-6
      testthat::expect_equal(actual[i, ], deriv[offset + links], tolerance = 1e-5)
    }
    # Fresh rows are evaluated on every chunk; reuse is not a row-score cache.
    first <- v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, plan = plan)$phi
    next_rows <- v57_ns("evaluate_augmented_if_chunk")(fit, 12:22, links, plan = plan)$phi
    testthat::expect_equal(first, v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, links)$phi, tolerance = 1e-12)
    testthat::expect_equal(next_rows, v57_ns("evaluate_augmented_if_chunk")(fit, 12:22, links)$phi, tolerance = 1e-12)
    testthat::expect_gt(max(abs(first - next_rows)), 1e-8)
    v57_ns("augmented_adjoint_replay_release")(plan)
    testthat::expect_false(plan$store$active)
    testthat::expect_null(plan$store$entries); testthat::expect_null(plan$store$seed)
    testthat::expect_equal(v57_collect(fit, links, plan), expected, tolerance = 1e-12)
  }
})

testthat::test_that("adjoint replay keeps partition/source refusals and original first-walk failures", {
  a <- v57_wide_case(); fit <- v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  B <- fit$control$max_chunk_rows; links <- rev(seq_along(fit$links))
  plan <- v57_ns("augmented_adjoint_replay_plan")(fit, links, B)
  testthat::expect_type(plan, "list")
  first <- v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, plan = plan)
  testthat::expect_true(any(!vapply(plan$store$entries, is.null, logical(1L))))
  testthat::expect_error(v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, rev(links), plan = plan), class = "cdid_cp04_integrity_error")
  testthat::expect_error(v57_ns("evaluate_augmented_if_chunk")(fit, 1:10, links, plan = plan), class = "cdid_cp04_integrity_error")
  foreign <- v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  testthat::expect_error(v57_ns("evaluate_augmented_if_chunk")(foreign, 1:11, links, plan = plan), class = "cdid_cp04_integrity_error")
  # Changed helper selects the unchanged route before any saved scalar solve
  # is trusted. No numerical exception is caught as an alternate estimator.
  helper <- v57_ns("reconstruction_zero_template")
  testthat::local_mocked_bindings(reconstruction_zero_template = function(...) helper(...), .package = "cdid")
  pre <- v57_ns("augmented_if_budget")(fit, 11L, links, FALSE)
  testthat::expect_null(v57_ns("augmented_adjoint_replay_access")(plan, fit, links, 11L, pre))
  testthat::expect_equal(v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, plan = plan)$phi,
    v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, links)$phi, tolerance = 1e-12)
})

testthat::test_that("ordinary adjoint memo eligibility leaves details and other models unchanged", {
  a <- v57_aug$rap_case()
  for (spec in list(c("within_cohort", "pooled"), c("separate", "separate"))) {
    fit <- v57_aug$rap_fit(a, imp = spec[1L], mode = spec[2L])
    testthat::expect_null(v57_ns("chained_full_if_plan")(fit, 1:3, fit$control$max_chunk_rows))
  }
  a <- v57_wide_case(); fit <- v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  links <- rev(seq_along(fit$links))
  plan <- v57_ns("chained_full_if_plan")(fit, links, fit$control$max_chunk_rows)
  testthat::expect_type(plan, "list")
  detailed <- v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, details = TRUE, plan = plan)
  testthat::expect_equal(detailed, v57_ns("evaluate_augmented_if_chunk")(fit, 1:11, links, details = TRUE), tolerance = 1e-12)
  testthat::expect_true(all(vapply(plan$store$entries, is.null, logical(1L))))
  dated <- v57_aug$rap_case(dated = TRUE, route_kind = "augmented")
  dated_fit <- v57_aug$rap_fit(dated, imp = "separate", mode = "pooled")
  testthat::expect_null(v57_ns("chained_full_if_plan")(dated_fit, 1:3, dated_fit$control$max_chunk_rows))
})


testthat::test_that("guarded coarse contraction accepts saved dense and structural compact records", {
  a <- v57_aug$rap_case(); fit <- v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  run <- v57_old$v53_run(fit, "augmented")
  testthat::expect_gt(run$used, 0L)
  testthat::expect_type(run$mutant, "list")
  x <- run$mutant$sample; api <- v57_old$v53_native(x)
  ptr <- api$create(); on.exit(api$close(ptr), add = TRUE)
  testthat::expect_type(ptr, "externalptr")
  ref <- do.call(v57_old$v53_scalar, x)
  contract <- function(adj, view = x$view, state = x$state, ids = x$ids)
    api$kernels[[2L]](ptr, view, state, adj, NULL, NULL, ids, NULL)
  testthat::expect_equal(contract(ref$adjoint), ref$phi, tolerance = 1e-12)
  structural <- sort(unique(unlist(lapply(x$state$groups[x$active], `[[`, "cell_ids"))))
  compact <- v57_ns("reconstruction_pack_plain_adjoint")(ref$adjoint$tilde_lambda,
    ref$adjoint$tilde_nu, structural, x$state$cells)
  testthat::expect_equal(contract(compact), ref$phi, tolerance = 1e-12)
  # This nonzero original-oracle direction must distinguish an incorrect memo.
  wrong <- ref$adjoint; wrong$tilde_lambda[run$mutant$cid] <- wrong$tilde_lambda[run$mutant$cid] + 1
  testthat::expect_gt(max(abs(contract(wrong) - ref$phi)), 1e-8)
  bad <- x$state; bad$B <- bad$B + 1L
  testthat::expect_null(contract(ref$adjoint, state = bad))
  decorated <- ref$adjoint; attr(decorated$tilde_lambda, "private") <- TRUE
  testthat::expect_null(contract(decorated))
  # Provided invalid centered inputs and changed source behavior keep the old
  # score preflight; the original scalar reference owns the resulting refusal.
  api$close(ptr); testthat::expect_null(contract(ref$adjoint))
})


testthat::test_that("memo admission and release preserve the original explicit payload limits", {
  a <- v57_wide_case(); fit <- v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  B <- fit$control$max_chunk_rows; links <- rev(seq_along(fit$links))
  plan <- v57_ns("chained_full_if_plan")(fit, links, B)
  testthat::expect_type(plan, "list")
  # Check the actual selected full and one-row inventory, not a named nominal
  # slab that might be absent when the old envelope remains selected.
  for (b in c(1L, B, plan$min_B)) {
    pre <- v57_ns("augmented_if_budget")(fit, b, links, FALSE)
    sp <- v57_ns("selection_if_budget")(fit$selection_backend, b, pre$q_ids, details = FALSE)
    inv <- v57_ns("reconstruction_ordinary_inventory")(fit, b, links, pre, sp)
    testthat::expect_identical(pre$scratch, inv$scratch)
  }
  pre <- v57_ns("augmented_if_budget")(fit, B, links, FALSE)
  need <- sum(c(fit$accounting$owned_entries, fit$accounting$borrowed_entries, pre$scratch))
  below <- fit; below$control$max_matrix_entries <- need - 1
  below <- v57_ns("ipw_owner_seal")(below, "augmented_backend")
  old <- tryCatch(v57_ns("augmented_if_budget")(below, B, links, FALSE), error = identity)
  new <- tryCatch(v57_ns("chained_full_if_plan")(below, links, B), error = identity)
  testthat::expect_s3_class(old, "error")
  testthat::expect_identical(class(new), class(old))
  testthat::expect_identical(conditionMessage(new), conditionMessage(old))
  exact <- fit; exact$control$max_matrix_entries <- need
  exact <- v57_ns("ipw_owner_seal")(exact, "augmented_backend")
  testthat::expect_type(v57_ns("chained_full_if_plan")(exact, links, B), "list")
  # Callback failure releases numeric memo entries before escaping the issued
  # producer. Original callback/source/error order is otherwise unchanged.
  factory <- v57_ns("chained_full_if_plan"); issued <- NULL
  testthat::local_mocked_bindings(chained_full_if_plan = function(...) {
    issued <<- factory(...); issued
  }, .package = "cdid")
  testthat::expect_error(v57_ns("cp04_replay_links")(fit, links, function(x) stop("V57 callback sentinel"),
    chunk_rows = B), "V57 callback sentinel", fixed = TRUE)
  testthat::expect_type(issued, "list")
  testthat::expect_false(issued$store$active)
  testthat::expect_null(issued$store$entries); testthat::expect_null(issued$store$seed)
})


testthat::test_that("the original tiny catalogue safely retains its genuinely selected legacy route", {
  a <- v57_aug$rap_case(); fit <- v57_aug$rap_fit(a, imp = "separate", mode = "pooled")
  links <- 3:1; B <- fit$control$max_chunk_rows
  pre <- v57_ns("augmented_if_budget")(fit, B, links, FALSE)
  sp <- v57_ns("selection_if_budget")(fit$selection_backend, B, pre$q_ids, details = FALSE)
  inv <- v57_ns("reconstruction_ordinary_inventory")(fit, B, links, pre, sp)
  testthat::expect_false(identical(pre$scratch, inv$scratch))
  plan <- v57_ns("chained_full_if_plan")(fit, links, B)
  testthat::expect_null(plan)
  actual <- v57_collect(fit, links, plan)
  expected <- v57_collect(fit, links, NULL)
  testthat::expect_identical(actual, expected)
  testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
  N <- nrow(a$Y); offset <- length(v57_aug$rap_point(fit)) - length(fit$links)
  for (i in c(1L, N)) {
    d <- rep(-1, N); d[i] <- N - 1
    derivative <- (v57_aug$rap_reference(a, fit, 1 + 1e-6 * d) -
      v57_aug$rap_reference(a, fit, 1 - 1e-6 * d)) / 2e-6
    testthat::expect_equal(actual[i, ], derivative[offset + links], tolerance = 1e-5)
  }
})
