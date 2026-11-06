# Exact row-invocation score reuse; independent nuisance oracles remain in the
# accepted reconstruction/augmented outcome-pooling test files.
rsc_ns <- function(name) getFromNamespace(name, "cdid")
rsc_case <- function() {
  g0 <- expand.grid(R3 = 0:1, R2 = 0:1, R1 = 0:1, x = c(-1, 1), k = 1:3)
  N <- nrow(g0); id <- seq_len(N); g <- c(3, 4, Inf)[g0$k]
  observed <- as.matrix(g0[c("R1", "R2", "R3")]) != 0; A <- rowSums(observed) >= 2L
  w <- 1 + g0$R1 / 10 + g0$x / 20
  X <- cbind(`(Intercept)` = 1, x = g0$x)
  Y <- vapply(1:3, function(d) d + g0$k * d / 4 +
    (d + .13 * d^2 + .07 * g0$k * d) * g0$x + sin(id + d) / 8, numeric(N))
  Y[!observed] <- NA_real_
  control <- rsc_ns("ipw_owner_control")(32L, 1e8, 2e8, 1e13, 1e8)
  data <- data.frame(id = rep(id, each = 3L), time = rep(1:3, N),
    y = as.vector(t(Y)), g = rep(g, each = 3L))
  panel <- rsc_ns("construct_ipw_panel_owner")(data, "id", "time", "y", "g", control)
  Xa <- X; Xa[!A, ] <- NA_real_
  frame <- rsc_ns("construct_ipw_frame_owner")(panel, id, Xa, w, control,
    frame_g = g, design_on_A_only = TRUE)
  reader <- rsc_ns("construct_cp04_frame_reader")(frame, id, X, g, w, control,
    U = rep(TRUE, N), rho = rep(1, N), universe = "F")
  defs <- list(list(s = 1, t = 3, g = 3, K = c(4, Inf)),
    list(s = 2, t = 3, g = 3, K = Inf), list(s = 1, t = 3, g = 4, K = Inf))
  views <- lapply(seq_along(defs), function(j) {
    d <- defs[[j]]
    rsc_ns("construct_ipw_link_view")(frame, d$s, d$t, d$g, d$K, control = control,
      link_id = paste0("score-cache-", j), functional_key = paste0("cmpcat:score-cache-", j))
  })
  list(reader = reader, views = views, X = X, control = control,
    recipe = rsc_ns("or_pooled_outcome_recipe")(colnames(X), TRUE, "x"))
}
rsc_fit <- function(a, kind, method = "dr", mode = "separate") {
  args <- list(reader = a$reader, link_views = a$views, control = a$control,
    glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE), rank_tol = 1e-10,
    outcome_recipe = a$recipe, treatment_columns = c("(Intercept)", "x"), X_frame = a$X,
    imputation_pooling = "within_cohort", est_method = method, outcome_pooling = mode)
  if (kind == "augmented") args <- c(args, list(selection_model = "independent_period_logit", selection_pooling = "within_cohort"))
  do.call(rsc_ns(if (kind == "augmented") "fit_augmented_backend" else "fit_imputation_backend"), args)
}
rsc_run <- function(fit, kind, rows, links, cached, details = TRUE) {
  factory <- rsc_ns("reconstruction_cached_cell_reader")
  calls <- integer(length(fit$cells))
  cache_cids <- integer(0L)
  counted_factory <- function(read, plan, n_cells) {
    force(read)
    cache_cids <<- plan$cids
    counted <- function(cid) {
      calls[cid] <<- calls[cid] + 1L
      read(cid)
    }
    if (cached) factory(counted, plan, n_cells) else counted
  }
  testthat::local_mocked_bindings(reconstruction_cached_cell_reader = counted_factory, .package = "cdid")
  fn <- rsc_ns(if (kind == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk")
  list(output = fn(fit, rows, links, details = details), calls = calls, cache_cids = cache_cids)
}
rsc_forward_run <- function(fit, kind, rows, links, forward, details = FALSE) {
  planner <- rsc_ns("reconstruction_stage1_forward_plan")
  pooled_solve <- rsc_ns("or_pooled_wls_gram_solve")
  plain_solve <- rsc_ns("or_wls_gram_solve")
  gram_calls <- 0L; plans <- list()
  testthat::local_mocked_bindings(
    reconstruction_stage1_forward_plan = function(...) {
      plan <- if (forward) planner(...) else list(enabled = FALSE, gids = integer(0L), entries = 0)
      plans[[length(plans) + 1L]] <<- plan
      plan
    },
    or_pooled_wls_gram_solve = function(...) {
      gram_calls <<- gram_calls + 1L
      pooled_solve(...)
    },
    or_wls_gram_solve = function(...) {
      gram_calls <<- gram_calls + 1L
      plain_solve(...)
    }, .package = "cdid")
  fn <- rsc_ns(if (kind == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk")
  out <- fn(fit, rows, links, details = details)
  list(output = out, gram_calls = gram_calls, plans = plans)
}

testthat::test_that("score reader caches exact records and NULL only within its own row invocation", {
  columns <- c("(Intercept)", paste0("long-", strrep("x", 400L)))
  meta <- list(shared_or = list(outcome_pooling = "separate"), control = list(max_chunk_rows = 4L),
    cells = list(list(group_id = 1L), list(group_id = 1L), list(group_id = 1L)),
    groups = list(list(p = 2L, kept_columns = columns)))
  X <- matrix(1, 4L, 2L, dimnames = list(NULL, columns))
  plan <- rsc_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, X, 100L)
  forward <- rsc_ns("reconstruction_stage1_forward_plan")(meta, 4L, 1L, X, 100L, FALSE)
  testthat::expect_true(forward$enabled)
  joint <- rsc_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, X, 100L,
    reserved_entries = forward$entries)
  testthat::expect_lte(forward$entries + joint$entries, 32 * 4L * 100L)
  testthat::expect_lt(length(joint$cids), length(plan$cids))
  testthat::expect_identical(rsc_ns("reconstruction_stage1_forward_plan")(
    meta, 12L, 1L, X, 100L, FALSE), forward)
  testthat::expect_false(rsc_ns("reconstruction_stage1_forward_plan")(meta, 4L, 1L, X, 1L, FALSE)$enabled)
  testthat::expect_false(rsc_ns("reconstruction_stage1_forward_plan")(meta, 4L, 1L, X, 100L, TRUE)$enabled)
  calls <- integer(3L)
  raw <- function(cid) {
    calls[cid] <<- calls[cid] + 1L
    if (cid == 2L) return(NULL)
    list(rr = 1:2, d_rr = X[1:2, , drop = FALSE], s_rr = c(cid, -cid))
  }
  reader <- rsc_ns("reconstruction_cached_cell_reader")(raw, plan, 3L)
  records <- lapply(c(3L, 2L, 1L, 2L, 3L, 1L), reader)
  testthat::expect_identical(calls, rep(1L, 3L))
  testthat::expect_null(records[[2L]]); testthat::expect_null(records[[4L]])
  testthat::expect_identical(records[[1L]], records[[5L]])
  testthat::expect_identical(records[[3L]], records[[6L]])
  testthat::expect_lte(unname(rsc_ns("cp04_payload")(records[1:3])["entries"]), plan$entries)
  fresh <- rsc_ns("reconstruction_cached_cell_reader")(raw, plan, 3L)
  fresh(2L); testthat::expect_identical(calls, c(1L, 2L, 1L))
  named <- X; rownames(named) <- paste0("named-", 1:4)
  no <- rsc_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, named, 100L)
  testthat::expect_false(no$enabled)
  testthat::expect_false(rsc_ns("reconstruction_stage1_forward_plan")(meta, 4L, 1L, named, 100L, FALSE)$enabled)
  testthat::expect_identical(rsc_ns("reconstruction_cached_cell_reader")(raw, no, 3L), raw)
  meta$shared_or$outcome_pooling <- "pooled"
  testthat::expect_true(rsc_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, X, 100L)$enabled)
  testthat::expect_false(rsc_ns("reconstruction_stage1_forward_plan")(meta, 4L, 1L, X, 100L, FALSE)$enabled)
  meta$shared_or <- NULL
  testthat::expect_true(rsc_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, X, 100L)$enabled)
  testthat::expect_false(rsc_ns("reconstruction_stage1_forward_plan")(meta, 4L, 1L, X, 100L, FALSE)$enabled)
})

testthat::test_that("pure and augmented complete IF channels reuse each funded score at most once", {
  a <- rsc_case()
  for (case in list(c("imputation", "reg"), c("imputation", "dr"), c("augmented", "dr"))) {
    kind <- case[1L]; fit <- rsc_fit(a, kind, case[2L])
    frozen <- list(cells = fit$cells, groups = fit$groups, shared_or = fit$shared_or,
      treatments = fit$treatments, reader_ref = fit$reader_ref)
    for (rows in list(1:32, c(33L, 1L, 18L))) {
      cached <- rsc_run(fit, kind, rows, c(3L, 1L), TRUE)
      old <- rsc_run(fit, kind, rows, c(3L, 1L), FALSE)
      testthat::expect_equal(cached$output, old$output, tolerance = 1e-12)
      if (length(rows) == 32L) {
        testthat::expect_gt(length(cached$cache_cids), 0L)
        testthat::expect_lte(max(cached$calls[cached$cache_cids]), 1L)
        testthat::expect_gt(sum(old$calls), sum(cached$calls))
      } else {
        # This batch cannot fund even the cache header inside its existing
        # scratch slice. It keeps the exact original reader and row partition.
        testthat::expect_length(cached$cache_cids, 0L)
        testthat::expect_identical(cached$calls, old$calls)
      }
    }
    full <- rsc_run(fit, kind, 1:32, 1:3, TRUE, FALSE)$output
    partial <- rsc_run(fit, kind, 32:1, c(3L, 1L), TRUE, FALSE)$output
    testthat::expect_equal(partial$phi, full$phi[32:1, c(3L, 1L), drop = FALSE], tolerance = 1e-10)
    budget <- rsc_ns(if (kind == "augmented") "augmented_if_budget" else "imputation_if_budget")
    pre <- budget(fit, 32L, 1:3, FALSE)
    testthat::expect_true(pre$score_cache$enabled || pre$stage1_forward$enabled)
    testthat::expect_true(all(pre$score_cache$cids %in% pre$cids))
    optional_entries <- pre$score_cache$entries + pre$stage1_forward$entries
    testthat::expect_lte(optional_entries, 32 * 32L * 3L)
    testthat::expect_lte(unname(rsc_ns("cp04_payload")(full)["entries"]), pre$output)
    # Preserve the original finding witness: adding an unfunded cache would
    # make the original admitted buffers fail. The implemented bounded cache
    # instead leaves every numeric envelope and public row plan unchanged.
    owned <- c(fit$accounting$owned_entries, fit$accounting$borrowed_entries)
    work <- c(fit$accounting$work_ops_bound, pre$work)
    limit <- ceiling(sum(owned) + pre$scratch)
    bounded <- rsc_ns("ipw_owner_control")(32L, fit$control$max_retained_entries,
      limit, fit$control$max_work_ops, fit$control$max_output_entries)
    testthat::expect_true(rsc_ns("cp04_budget")(bounded, owned,
      pre$scratch, work, pre$output, "funded resource witness"))
    testthat::expect_error(rsc_ns("cp04_budget")(bounded, owned, pre$scratch + optional_entries,
      work, pre$output, "unfunded resource witness"), class = "cdid_ipw_owner_bound_error")
    retained <- fit$control$max_matrix_entries - sum(owned) - pre$output - pre$scratch
    plan <- rsc_ns("cp04_preflight_links")(fit, 1:3, retained, chunk_rows = 32L)
    original <- function() {
      testthat::local_mocked_bindings(reconstruction_score_cache_plan = function(...) {
        list(enabled = FALSE, cids = integer(0L), entries = 0)
      }, reconstruction_stage1_forward_plan = function(...) {
        list(enabled = FALSE, gids = integer(0L), entries = 0)
      }, .package = "cdid")
      list(pre = budget(fit, 32L, 1:3, FALSE),
        plan = rsc_ns("cp04_preflight_links")(fit, 1:3, retained, chunk_rows = 32L))
    }
    uncached <- original()
    for (field in c("output", "scratch", "work")) testthat::expect_identical(pre[[field]], uncached$pre[[field]])
    testthat::expect_identical(plan, uncached$plan)
    testthat::expect_identical(plan$chunk_rows, 32)
    testthat::expect_error(rsc_ns("cp04_preflight_links")(fit, 1:3, retained + 1,
      chunk_rows = 32L), class = "cdid_ipw_owner_bound_error")

    # Ordinary forward mode reuses the already accepted detailed score/Gram
    # formulas; the declined reference retains the exact V6 adjoint route.
    reused <- rsc_forward_run(fit, kind, 1:32, 1:3, TRUE)
    adjoint <- rsc_forward_run(fit, kind, 1:32, 1:3, FALSE)
    detailed <- rsc_forward_run(fit, kind, 1:32, 1:3, TRUE, TRUE)
    testthat::expect_length(reused$plans, 1L)
    testthat::expect_true(all(vapply(reused$plans, `[[`, logical(1L), "enabled")))
    testthat::expect_equal(reused$output, adjoint$output, tolerance = 1e-10)
    testthat::expect_equal(reused$output$phi, detailed$output$phi, tolerance = 1e-10)
    testthat::expect_lt(reused$gram_calls, adjoint$gram_calls)
    testthat::expect_false(any(vapply(detailed$plans, `[[`, logical(1L), "enabled")))
    testthat::expect_identical(pre$stage1_forward$gids, pre$gids)
    declined <- rsc_forward_run(fit, kind, c(33L, 1L, 18L), c(3L, 1L), TRUE)
    legacy <- rsc_forward_run(fit, kind, c(33L, 1L, 18L), c(3L, 1L), FALSE)
    testthat::expect_false(any(vapply(declined$plans, `[[`, logical(1L), "enabled")))
    testthat::expect_identical(declined$output, legacy$output)
    testthat::expect_identical(declined$gram_calls, legacy$gram_calls)
    testthat::expect_identical(list(cells = fit$cells, groups = fit$groups, shared_or = fit$shared_or,
      treatments = fit$treatments, reader_ref = fit$reader_ref), frozen)
  }
})


# V36 changes only score-read reuse for generic ordinary replay. Gram calls
# and all old scalar adjoint/finite checks remain the accepted reference path.
rsc36_run <- function(fit, kind, rows, links, cached, details = FALSE) {
  factory <- rsc_ns("reconstruction_cached_cell_reader")
  pooled <- rsc_ns("or_pooled_wls_gram_solve"); plain <- rsc_ns("or_wls_gram_solve")
  calls <- nulls <- integer(length(fit$cells)); plans <- list(); grams <- 0L
  testthat::local_mocked_bindings(
    reconstruction_cached_cell_reader = function(read, plan, n_cells) {
      force(read)
      plans[[length(plans) + 1L]] <<- plan
      counted <- function(cid) {
        calls[cid] <<- calls[cid] + 1L
        record <- read(cid)
        if (is.null(record)) nulls[cid] <<- nulls[cid] + 1L
        record
      }
      if (cached) factory(counted, plan, n_cells) else counted
    },
    or_pooled_wls_gram_solve = function(...) {
      grams <<- grams + 1L
      pooled(...)
    },
    or_wls_gram_solve = function(...) {
      grams <<- grams + 1L
      plain(...)
    }, .package = "cdid")
  fn <- rsc_ns(if (kind == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk")
  output <- fn(fit, rows, links, details = details)
  list(output = output, calls = calls, nulls = nulls, plans = plans, grams = grams)
}

testthat::test_that("generic ordinary IPW and pooled scores reuse records without forward solves or new bounds", {
  a <- rsc_case()
  for (case in list(c("imputation", "ipw"), c("imputation", "reg"),
      c("imputation", "dr"), c("augmented", "dr"))) {
    kind <- case[1L]; fit <- rsc_fit(a, kind, case[2L], "pooled")
    if (case[2L] == "ipw") testthat::expect_null(fit$shared_or) else
      testthat::expect_false(rsc_ns("reconstruction_separate")(fit$shared_or))
    frozen <- list(cells = fit$cells, groups = fit$groups, shared_or = fit$shared_or,
      treatments = fit$treatments, reader_ref = fit$reader_ref)
    budget <- rsc_ns(if (kind == "augmented") "augmented_if_budget" else "imputation_if_budget")
    evaluate <- rsc_ns(if (kind == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk")
    pre <- budget(fit, 32L, 1:3, FALSE)
    testthat::expect_true(pre$score_cache$enabled)
    testthat::expect_null(pre$stage1_forward)
    testthat::expect_true(all(pre$score_cache$cids %in% pre$cids))
    cache_slice <- if (kind == "augmented") 3040 else 2016
    testthat::expect_lte(pre$score_cache$entries, cache_slice * 32L * 3L)
    detailed_pre <- budget(fit, 32L, 1:3, TRUE)
    testthat::expect_null(detailed_pre$score_cache)
    testthat::expect_null(detailed_pre$stage1_forward)
    actual <- old <- matrix(0, nrow(a$X), 3L)
    saved_grams <- 0L; null_reused <- FALSE
    for (rows in list(1:32, 33:48)) {
      cached <- rsc36_run(fit, kind, rows, 1:3, TRUE)
      raw <- rsc36_run(fit, kind, rows, 1:3, FALSE)
      detailed <- rsc36_run(fit, kind, rows, 1:3, TRUE, TRUE)
      testthat::expect_equal(cached$output, raw$output, tolerance = 1e-12)
      testthat::expect_equal(cached$output$phi, detailed$output$phi, tolerance = 1e-10)
      testthat::expect_identical(cached$grams, raw$grams)
      saved_grams <- saved_grams + cached$grams
      testthat::expect_length(cached$plans, 1L)
      ids <- cached$plans[[1L]]$cids
      if (length(ids)) {
        testthat::expect_true(all(cached$calls[ids] <= 1L))
        null_reused <- null_reused || any(raw$nulls[ids] > cached$nulls[ids])
      }
      if (length(rows) == 32L) testthat::expect_gt(sum(raw$calls), sum(cached$calls))
      actual[rows, ] <- cached$output$phi; old[rows, ] <- raw$output$phi
    }
    if (case[2L] == "ipw") testthat::expect_identical(saved_grams, 0L) else
      testthat::expect_gt(saved_grams, 0L)
    testthat::expect_true(null_reused)
    testthat::expect_equal(crossprod(actual) / nrow(actual)^2,
      crossprod(old) / nrow(old)^2, tolerance = 1e-10)
    reordered <- rsc36_run(fit, kind, 32:1, c(3L, 1L), TRUE)
    testthat::expect_equal(unname(reordered$output$phi),
      unname(actual[32:1, c(3L, 1L), drop = FALSE]), tolerance = 1e-10)
    testthat::expect_identical(dimnames(reordered$output$phi), list(NULL,
      vapply(fit$links[c(3L, 1L)], `[[`, character(1L), "link_id")))
    # V39 can fund the one-row tail from the larger generic slab. Numeric
    # output, Gram work and partition remain equal to the old raw reader.
    tiny <- rsc36_run(fit, kind, 48L, 1:3, TRUE)
    tiny_raw <- rsc36_run(fit, kind, 48L, 1:3, FALSE)
    testthat::expect_lte(tiny$plans[[1L]]$entries, cache_slice * 3L)
    testthat::expect_equal(tiny$output, tiny_raw$output, tolerance = 1e-12)
    testthat::expect_identical(tiny$grams, tiny_raw$grams)
    testthat::expect_true(all(tiny$calls <= tiny_raw$calls))
    old_plan <- function() {
      testthat::local_mocked_bindings(reconstruction_score_cache_plan = function(...) {
        list(enabled = FALSE, cids = integer(0L), entries = 0)
      }, .package = "cdid")
      budget(fit, 32L, 1:3, FALSE)
    }
    old_pre <- old_plan()
    for (field in c("cids", "gids", "bids", "tids", "output", "scratch", "work"))
      testthat::expect_identical(pre[[field]], old_pre[[field]])
    owned <- c(fit$accounting$owned_entries, fit$accounting$borrowed_entries)
    retained <- fit$control$max_matrix_entries - sum(owned) - pre$output - pre$scratch
    admitted <- rsc_ns("cp04_preflight_links")(fit, 1:3, retained, chunk_rows = 32L)
    testthat::expect_identical(admitted$chunk_rows, 32)
    testthat::expect_error(rsc_ns("cp04_preflight_links")(fit, 1:3, retained + 1,
      chunk_rows = 32L), class = "cdid_ipw_owner_bound_error")
    testthat::expect_lte(unname(rsc_ns("cp04_payload")(reordered$output)["entries"]), pre$output)
    broken <- fit; broken$cells[[1L]]$alpha <- broken$cells[[1L]]$alpha + 1
    testthat::expect_error(evaluate(broken, 1:32, 1:3),
      class = if (kind == "augmented") "cdid_augmented_integrity_error" else "cdid_imputation_integrity_error")
    testthat::expect_identical(list(cells = fit$cells, groups = fit$groups, shared_or = fit$shared_or,
      treatments = fit$treatments, reader_ref = fit$reader_ref), frozen)
  }
})

testthat::test_that("cache bookkeeping includes the full cell-position map with sparse selection and long names", {
  nc <- 10000L
  columns <- c("(Intercept)", paste0("long-", strrep("column", 80L)))
  meta <- list(shared_or = NULL, control = list(max_chunk_rows = 4L),
    cells = rep(list(list(group_id = 1L)), nc), groups = list(list(p = 2L, kept_columns = columns)))
  X <- matrix(1, 4L, 2L, dimnames = list(NULL, columns))
  plan <- rsc_ns("reconstruction_score_cache_plan")(meta, 4L, c(2L, 9999L), X, 100L)
  testthat::expect_true(plan$enabled)
  raw <- function(cid) if (cid == 2L) NULL else list(rr = 1:4, d_rr = X, s_rr = rep(0, 4))
  reader <- rsc_ns("reconstruction_cached_cell_reader")(raw, plan, nc)
  for (cid in plan$cids) reader(cid)
  env <- environment(reader)
  testthat::expect_length(env$position, nc)
  testthat::expect_identical(which(env$position > 0L), plan$cids)
  testthat::expect_length(env$cache, length(plan$cids))
  testthat::expect_true(all(env$seen))
  # Character names/records are in plan$entries; the global map/slots and
  # planning vectors use unchanged generic metadata, even for two selected IDs.
  metadata <- 16384 + 1024 * (length(plan$cids) + length(meta$groups)) + 64 * nc * 2
  payload <- rsc_ns("cp04_payload")(list(position = env$position, seen = env$seen,
    cache = env$cache, plan = plan))
  testthat::expect_lte(unname(payload["entries"]) + 1024, metadata + plan$entries)
  named <- X; rownames(named) <- paste0("row-", 1:4)
  no <- rsc_ns("reconstruction_score_cache_plan")(meta, 4L, c(2L, 9999L), named, 100L)
  testthat::expect_false(no$enabled)
  testthat::expect_identical(rsc_ns("reconstruction_cached_cell_reader")(raw, no, nc), raw)
  unfunded <- rsc_ns("reconstruction_score_cache_plan")(meta, 1L, c(2L, 9999L), X, 1L)
  testthat::expect_false(unfunded$enabled)
  testthat::expect_identical(rsc_ns("reconstruction_cached_cell_reader")(raw, unfunded, nc), raw)
})
