v39_ns <- function(name) getFromNamespace(name, "cdid")
v39_helpers <- function(file, names) {
  e <- new.env(parent = asNamespace("cdid"))
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% names) eval(x, e)
  stopifnot(all(vapply(names, exists, logical(1), envir = e, inherits = FALSE)))
  e
}
v39_pure <- v39_helpers("test-reconstructed-outcome-pooling.R",
  c("rcp_ns", "rcp_gc", "rcp_case", "rcp_fit", "rcp_reference", "rcp_point"))
v39_aug <- v39_helpers("test-augmented-outcome-pooling.R",
  c("rap_ns", "rap_gc", "rap_case", "rap_fit", "rap_reference", "rap_point"))
v39_bound_old <- new.env(parent = asNamespace("cdid"))
sys.source(testthat::test_path("fixtures", "reconstruction-v39-bound-reference.R"), v39_bound_old)

# Scalar reference retains the exact current equations but never reuses a
# zero template or prunes a record. It also uses the original arbitrary-input
# payload packer. No solver/RHS/formula/tolerance is changed by this adapter.
v39_env <- function(scalar) {
  e <- new.env(parent = asNamespace("cdid"))
  for (name in c("fit_imputation_backend", "imputation_fit_core", "fit_augmented_backend",
      "evaluate_imputation_if_chunk", "evaluate_imputation_layout_if_chunk", "imputation_evaluate_if_chunk_core",
      "evaluate_augmented_if_chunk", "evaluate_augmented_layout_if_chunk", "augmented_evaluate_if_chunk_core")) {
    f <- v39_ns(name); environment(f) <- e; assign(name, f, e)
  }
  if (scalar) {
    e$reconstruction_zero_template <- function(...) NULL
    e$reconstruction_complete_cache_witness <- function(...) NULL
    e$reconstruction_pack_plain_adjoint <- v39_ns("reconstruction_pack_adjoint")
    e$identical <- function(x, y) {
      if (is.list(x) && all(c("rr", "d_rr", "s_rr") %in% names(x))) FALSE else base::identical(x, y)
    }
  }
  e
}
v39_fit <- function(a, kind, method = "dr", imp = "separate", mode = "separate",
                    sel_model = "pair_logit", scalar = FALSE, corrupt_after = 0L) {
  e <- v39_env(scalar); solves <- packs <- bound_calls <- 0L; records <- list()
  for (name in c("or_pooled_wls_gram_solve", "or_wls_gram_solve")) local({
    original <- v39_ns(name)
    assign(name, function(...) {
      if (exists("rhs_mat", parent.frame(), inherits = FALSE) && exists("tilde_lambda", parent.frame(), inherits = FALSE))
        solves <<- solves + 1L
      original(...)
    }, e)
  })
  pack <- if (scalar) v39_ns("reconstruction_pack_adjoint") else v39_ns("reconstruction_pack_plain_adjoint")
  e$reconstruction_pack_plain_adjoint <- function(lambda, nu, cids, cells) {
    packs <<- packs + 1L
    dense <- list(tilde_lambda = lambda, tilde_nu = nu)
    testthat::expect_equal(unname(v39_ns("cp04_payload")(dense)["entries"]),
      length(lambda) + 9 * length(nu) + sum(lengths(nu)) + 55)
    out <- pack(lambda, nu, cids, cells)
    records[[length(records) + 1L]] <<- list(dense = dense, cids = cids)
    if (corrupt_after == packs) {
      f <- parent.frame()
      while (!exists("cells_fitted", f, inherits = FALSE) && !identical(f, emptyenv())) f <- parent.env(f)
      if (exists("cells_fitted", f, inherits = FALSE)) {
        x <- get("cells_fitted", f, inherits = FALSE)
        # First arm active -> first inactive later; second-arm mutation also
        # invalidates a zero template from an earlier successful inactive arm.
        cid <- if (corrupt_after == 1L) cids[1L] else tail(cids, 1L)
        x[[cid]]$mean_d[] <- Inf
        assign("cells_fitted", x, envir = f)
      }
    }
    out
  }
  for (name in c("imputation_retained_core_bounds", "augmented_retained_core_bounds")) local({
    original <- v39_ns(name); old_name <- sub("bounds$", "bound", name)
    assign(name, function(...) {
      args <- list(...); value <- original(...); bound_calls <<- bound_calls + 1L
      for (dense in c(FALSE, TRUE)) testthat::expect_identical(value[[if (dense) "dense" else "compact"]],
        do.call(v39_bound_old[[old_name]], c(args, list(dense_adjoints = dense))))
      value
    }, e)
  })
  args <- list(reader = a$reader, link_views = a$views, control = a$control,
    glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE), rank_tol = 1e-10,
    outcome_recipe = a$recipe, treatment_columns = c("(Intercept)", "treatment_v"),
    reference = a$reference, imputation_pooling = imp, outcome_pooling = mode)
  if (kind == "augmented") args <- c(args, list(selection_model = sel_model,
    selection_pooling = "across_cohorts", X_frame = if (is.null(a$reference)) a$X_full else NULL)) else
    args <- c(args, list(est_method = method, X_frame = if (is.null(a$reference)) a$X_full else NULL))
  fit <- do.call(e[[if (kind == "augmented") "fit_augmented_backend" else "fit_imputation_backend"]], args)
  list(fit = fit, solves = solves, packs = packs, records = records, bound_calls = bound_calls)
}

v39_core <- function(fit, kind, rows, links, scalar = FALSE, cache = "full", damage = NULL) {
  e <- v39_env(scalar); factory <- v39_ns("reconstruction_cached_cell_reader"); products <- 0L
  e[["%*%"]] <- function(x, y) { products <<- products + 1L; base::`%*%`(x, y) }
  e$reconstruction_cached_cell_reader <- function(read, plan, n_cells) {
    force(read)
    if (cache == "none") plan <- list(enabled = FALSE, cids = integer(0L), entries = 0)
    if (cache == "partial" && isTRUE(plan$enabled)) plan$cids <- fit$links[[1L]]$target_cell_id
    wrapped <- factory(read, plan, n_cells); calls <- integer(n_cells)
    function(cid) {
      value <- wrapped(cid); calls[cid] <<- calls[cid] + 1L
      if (!is.null(damage) && cid == fit$links[[1L]]$target_cell_id && !is.null(value)) {
        if (damage == "late" && calls[cid] >= 2L) value$s_rr[] <- Inf
        if (damage == "late_mean" && calls[cid] >= 2L) {
          f <- parent.frame()
          while (!exists("core_fit", f, inherits = FALSE) && !identical(f, emptyenv())) f <- parent.env(f)
          if (exists("core_fit", f, inherits = FALSE)) {
            x <- get("core_fit", f, inherits = FALSE); x$cells[[cid]]$mean_d[] <- Inf
            assign("core_fit", x, envir = f)
          }
        }
        if (damage == "score") value$s_rr[] <- Inf
        if (damage == "centered") value$d_rr[] <- -.Machine$double.xmax
      }
      value
    }
  }
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  st <- reader$state; pre <- v39_ns(if (kind == "augmented") "augmented_if_budget" else "imputation_if_budget")(
    fit, length(rows), links, FALSE)
  # Keep the evaluator's cache declaration consistent with this test reader.
  # Changing only the factory's local plan does not change its caller's pre.
  if (cache == "none") pre$score_cache <- list(enabled = FALSE, cids = integer(0L), entries = 0)
  if (cache == "partial" && isTRUE(pre$score_cache$enabled))
    pre$score_cache$cids <- fit$links[[1L]]$target_cell_id
  if (identical(damage, "centered")) fit$cells[[fit$links[[1L]]$target_cell_id]]$mean_d[] <- .Machine$double.xmax
  pair_reader <- function(r, pair, read_outcomes = TRUE) v39_ns("cp04_pair_rows")(reader, r, pair, read_outcomes)
  out <- if (kind == "augmented") e$augmented_evaluate_if_chunk_core(st = st,
    X_F = v39_ns("augmented_backend_X_F")(fit), pair_reader = pair_reader, backend = fit,
    rows = rows, links = links, scale = "frame", details = FALSE, sets = pre) else
    e$imputation_evaluate_if_chunk_core(st = st, X_F = v39_ns("imputation_backend_X_F")(fit, reader),
      pair_reader = pair_reader, core_fit = fit, recipe = fit$outcome_recipe, or_frame_index = fit$or_frame_index,
      trt_recipe = fit$treatment_recipe, est_method = fit$est_method, rows = rows, links = links,
      scale = "frame", details = FALSE, sets = pre)
  attr(out, "v39_product_calls") <- products
  out
}

testthat::test_that("active scalar reconstruction preserves old bounds, fitted adjoints and full-F derivatives", {
  specs <- list(list(kind = "imputation", method = "ipw", imp = "separate", mode = "separate", dated = FALSE),
    list(kind = "imputation", method = "reg", imp = "within_cohort", mode = "separate", dated = FALSE),
    list(kind = "imputation", method = "dr", imp = "across_cohorts", mode = "pooled", dated = TRUE),
    list(kind = "augmented", method = "dr", imp = "separate", mode = "pooled", dated = FALSE),
    list(kind = "augmented", method = "dr", imp = "within_cohort", mode = "separate", dated = FALSE))
  for (spec in specs) {
    kind <- spec$kind
    a <- if (kind == "augmented") v39_aug$rap_case(route_kind = "augmented") else
      v39_pure$rcp_case(dated = spec$dated, route_kind = "imputation")
    args <- list(a = a, kind = kind, method = spec$method, imp = spec$imp, mode = spec$mode,
      sel_model = if (spec$mode == "pooled") "independent_period_logit" else "pair_logit")
    current <- do.call(v39_fit, args); old <- do.call(v39_fit, c(args, list(scalar = TRUE)))
    f <- current$fit; g <- old$fit; N <- nrow(a$Y); links <- rev(seq_along(f$links))
    testthat::expect_identical(current$records, old$records)
    testthat::expect_equal(current$fit$links, old$fit$links, tolerance = 0)
    for (field in c("retained_bound", "fit_scratch_bound", "work_ops_bound"))
      testthat::expect_identical(f$accounting[[field]], g$accounting[[field]])
    testthat::expect_lte(current$solves, old$solves)
    if (spec$imp == "separate") testthat::expect_lt(current$solves, old$solves)
    testthat::expect_equal(current$bound_calls, 1L)
    evaluate <- v39_ns(if (kind == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk")
    raw_env <- v39_env(TRUE); reference <- raw_env[[if (kind == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk"]]
    actual <- plain <- matrix(0, N, length(links))
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(first, min(N, first + 10L)))
      actual[rows, ] <- evaluate(f, rows, links)$phi
      plain[rows, ] <- reference(g, rows, links)$phi
      if (!spec$dated) testthat::expect_equal(actual[rows, , drop = FALSE],
        unname(evaluate(f, rows, links, details = TRUE)$phi), tolerance = 1e-10)
    }
    testthat::expect_equal(actual, plain, tolerance = 1e-10)
    testthat::expect_equal(crossprod(actual), crossprod(plain), tolerance = 1e-10)
    point <- if (kind == "augmented") v39_aug$rap_point else v39_pure$rcp_point
    raw <- if (kind == "augmented") v39_aug$rap_reference else v39_pure$rcp_reference
    testthat::expect_equal(point(f), raw(a, f), tolerance = 2e-8)
    for (i in c(1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      deriv <- (raw(a, f, 1 + 1e-6 * delta) - raw(a, f, 1 - 1e-6 * delta)) / 2e-6
      offset <- length(point(f)) - length(links)
      testthat::expect_equal(actual[i, ], deriv[offset + links], tolerance = 1e-5)
    }
    if (!spec$dated) for (cache in c("full", "partial", "none")) {
      current_core <- v39_core(f, kind, 11:1, links, cache = cache)
      old_core <- v39_core(f, kind, 11:1, links, scalar = TRUE, cache = cache)
      testthat::expect_equal(current_core$phi_frame, old_core$phi_frame, tolerance = 1e-10)
      if (spec$method == "ipw") {
        if (cache == "full") testthat::expect_lt(attr(current_core, "v39_product_calls"), attr(old_core, "v39_product_calls"))
        if (cache == "none") testthat::expect_identical(attr(current_core, "v39_product_calls"), attr(old_core, "v39_product_calls"))
      }
    }
  }
})

testthat::test_that("inactive scalar reuse waits for the first actual zero and invalidates changed source inputs", {
  a <- v39_pure$rcp_case(route_kind = "imputation")
  for (at in c(1L, 2L)) {
    actual <- tryCatch(v39_fit(a, "imputation", corrupt_after = at), error = identity)
    old <- tryCatch(v39_fit(a, "imputation", scalar = TRUE, corrupt_after = at), error = identity)
    testthat::expect_s3_class(actual, "error")
    testthat::expect_identical(class(actual), class(old))
    testthat::expect_identical(conditionMessage(actual), conditionMessage(old))
  }
  result <- v39_fit(a, "imputation", method = "ipw")$fit
  for (cache in c("full", "partial", "none")) for (damage in c("late", "late_mean", "score", "centered")) {
    actual <- tryCatch(v39_core(result, "imputation", 1:11, 1:3, cache = cache, damage = damage), error = identity)
    old <- tryCatch(v39_core(result, "imputation", 1:11, 1:3, scalar = TRUE, cache = cache, damage = damage), error = identity)
    testthat::expect_s3_class(actual, "error")
    testthat::expect_identical(class(actual), class(old))
    testthat::expect_identical(conditionMessage(actual), conditionMessage(old))
  }
})

testthat::test_that("plain constructor size and generic cache funding stay inside exact original slabs", {
  cells <- lapply(1:9, function(j) list(group_id = as.integer(j)))
  for (p in c(0L, 1L, 4L)) {
    lambda <- numeric(9L); nu <- lapply(cells, function(c) numeric(p))
    testthat::expect_identical(v39_ns("reconstruction_pack_plain_adjoint")(lambda, nu, c(1L, 3L), cells),
      v39_ns("reconstruction_pack_adjoint")(lambda, nu, c(1L, 3L), cells))
    testthat::expect_equal(unname(v39_ns("cp04_payload")(list(tilde_lambda = lambda, tilde_nu = nu))["entries"]),
      9 + 9 * 9 + sum(lengths(nu)) + 55)
    scratch <- v39_ns("reconstruction_adjoint_construction_scratch")(9L, 9L, p, 9L)
    testthat::expect_lte(scratch, 128 * 9L * (p + 2)^2 + 8192)
  }
  # Decorated private callers retain the old arbitrary-input entry point.
  lambda <- setNames(numeric(9L), paste0("id", 1:9)); nu <- lapply(cells, function(c) setNames(0, "slope"))
  decorated <- v39_ns("reconstruction_pack_adjoint")(lambda, nu, 1:9, cells)
  testthat::expect_lte(unname(v39_ns("cp04_payload")(decorated)["entries"]),
    unname(v39_ns("cp04_payload")(list(tilde_lambda = lambda, tilde_nu = nu))["entries"]))
})
