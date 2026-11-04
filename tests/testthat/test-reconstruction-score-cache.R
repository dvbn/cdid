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
rsc_fit <- function(a, kind, method = "dr") {
  args <- list(reader = a$reader, link_views = a$views, control = a$control,
    glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE), rank_tol = 1e-10,
    outcome_recipe = a$recipe, treatment_columns = c("(Intercept)", "x"), X_frame = a$X,
    imputation_pooling = "within_cohort", est_method = method, outcome_pooling = "separate")
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

testthat::test_that("score reader caches exact records and NULL only within its own row invocation", {
  columns <- c("(Intercept)", paste0("long-", strrep("x", 400L)))
  meta <- list(shared_or = list(outcome_pooling = "separate"), control = list(max_chunk_rows = 4L),
    cells = list(list(group_id = 1L), list(group_id = 1L), list(group_id = 1L)),
    groups = list(list(p = 2L, kept_columns = columns)))
  X <- matrix(1, 4L, 2L, dimnames = list(NULL, columns))
  plan <- rsc_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, X, 100L)
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
  testthat::expect_identical(rsc_ns("reconstruction_cached_cell_reader")(raw, no, 3L), raw)
  meta$shared_or$outcome_pooling <- "pooled"
  testthat::expect_false(rsc_ns("reconstruction_score_cache_plan")(meta, 4L, 1:3, X, 100L)$enabled)
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
    testthat::expect_true(pre$score_cache$enabled)
    testthat::expect_true(all(pre$score_cache$cids %in% pre$cids))
    testthat::expect_lte(pre$score_cache$entries, 32 * 32L * 3L)
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
    testthat::expect_error(rsc_ns("cp04_budget")(bounded, owned, pre$scratch + pre$score_cache$entries,
      work, pre$output, "unfunded resource witness"), class = "cdid_ipw_owner_bound_error")
    retained <- fit$control$max_matrix_entries - sum(owned) - pre$output - pre$scratch
    plan <- rsc_ns("cp04_preflight_links")(fit, 1:3, retained, chunk_rows = 32L)
    original <- function() {
      testthat::local_mocked_bindings(reconstruction_score_cache_plan = function(...) {
        list(enabled = FALSE, cids = integer(0L), entries = 0)
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
    testthat::expect_identical(list(cells = fit$cells, groups = fit$groups, shared_or = fit$shared_or,
      treatments = fit$treatments, reader_ref = fit$reader_ref), frozen)
  }
})
