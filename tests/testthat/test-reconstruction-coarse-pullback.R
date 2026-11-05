v53_ns <- function(name) getFromNamespace(name, "cdid")
v53_helpers <- function(file, wanted) {
  e <- new.env(parent = asNamespace("cdid"))
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE))); e
}
v53_pure <- v53_helpers("test-reconstructed-outcome-pooling.R",
  c("rcp_ns", "rcp_gc", "rcp_case", "rcp_fit", "rcp_point", "rcp_reference"))
v53_aug <- v53_helpers("test-augmented-outcome-pooling.R",
  c("rap_ns", "rap_gc", "rap_case", "rap_fit", "rap_point", "rap_reference"))

# Exact old group loop, including original first-hit and V50 suffix calls.
v53_scalar <- function(view, state, lambda, nu, active, ids, templates, centered = NULL) {
  cells <- state$cells; groups <- state$groups; N <- state$N
  tl <- numeric(length(cells)); tn <- vector("list", length(groups)); next_gid <- 1L
  solves <- hits <- suffix <- 0L
  for (gid in seq_along(groups)) {
    if (gid < next_gid) next
    gf <- groups[[gid]]
    cached <- v53_ns("reconstruction_zero_template")(templates, gid, gid %in% active, gf, cells, N)
    if (!is.null(cached)) {
      hits <- hits + 1L; tl[gf$cell_ids] <- cached$lambda; tn[[gid]] <- cached$nu
      run <- if (gid < length(groups)) v53_ns("reconstruction_zero_template_run")(
        templates, groups, cells, active, N, gid + 1L, v53_ns("reconstruction_zero_template")) else NULL
      if (!is.null(run)) {
        tl[run$cell_ids] <- run$lambda; tn[run$group_ids] <- run$nu
        suffix <- suffix + length(run$group_ids); next_gid <- run$next_gid
      }
      next
    }
    rhs <- numeric(gf$p)
    for (cid in gf$cell_ids) {
      cf <- cells[[cid]]; lam <- if (state$intercept) lambda[cid] else 0
      if (state$intercept) tl[cid] <- lam / cf$mass
      if (gf$p) rhs <- rhs + nu[[cid]] - lam * cf$mean_d
    }
    tn[[gid]] <- if (!gf$p) numeric() else {
      solves <- solves + 1L
      as.double(if (state$intercept) v53_ns("or_pooled_wls_gram_solve")(gf$fit, rhs, normalizer = N) else
        v53_ns("or_wls_gram_solve")(gf$fit, rhs, normalizer = N))
    }
    if (!gid %in% active) templates[[gid]] <- list(lambda = tl[gf$cell_ids], nu = tn[[gid]], group = gf, cells = cells[gf$cell_ids], N = N)
  }
  out <- numeric(state$B)
  for (cid in ids) {
    rec <- view$cache[[view$position[cid]]]; if (is.null(rec)) next
    cf <- cells[[cid]]; gf <- groups[[cf$group_id]]
    lam <- if (state$intercept) tl[cid] else 0.0; val <- rep(lam, length(rec$rr))
    if (gf$p && length(tn[[cf$group_id]])) {
      d <- if (state$intercept) { if (!is.null(centered)) centered[[view$position[cid]]] else sweep(rec$d_rr, 2L, cf$mean_d, FUN = "-") } else rec$d_rr
      val <- val + as.vector(d %*% tn[[cf$group_id]])
    }
    out[rec$rr] <- out[rec$rr] + rec$s_rr * val
  }
  list(phi = out, templates = templates, adjoint = list(tilde_lambda = tl, tilde_nu = tn), counts = c(solves, hits, suffix))
}

v53_run <- function(fit, kind, enabled = TRUE) {
  factory <- v53_ns("reconstruction_coarse_operation")
  issued <- calls <- used <- 0L; suffix <- 0L; samples <- operations <- list(); mutant <- NULL
  testthat::local_mocked_bindings(reconstruction_coarse_operation = function(view, state, centered = NULL) {
    op <- if (enabled) factory(view, state, centered) else NULL
    if (!is.null(op)) {
      issued <<- issued + 1L; pullback <- op$pullback
      op$pullback <- function(view, state, lambda, nu, groups, ids, templates) {
        out <- pullback(view, state, lambda, nu, groups, ids, templates); calls <<- calls + 1L
        if (!is.null(out)) {
          used <<- used + 1L; suffix <<- suffix + out$counts[3L]
          ref <- v53_scalar(view, state, lambda, nu, groups, ids, templates, centered)
          testthat::expect_identical(out$counts, ref$counts)
          testthat::expect_identical(names(out), names(ref))
          testthat::expect_equal(out$adjoint, ref$adjoint, tolerance = 1e-12)
          testthat::expect_equal(out$phi, ref$phi, tolerance = 1e-12)
          testthat::expect_equal(out$templates, ref$templates, tolerance = 1e-12)
          input <- list(view = view, state = state, lambda = lambda, nu = nu,
            active = groups, ids = ids, templates = templates, centered = centered)
          if (!length(samples)) samples[[1L]] <<- input
          # Select a nonzero test direction using ONLY the independent old-R
          # oracle on existing observed cells. The first treated-only chunk
          # legitimately has zero response for every control-cell loading.
          if (is.null(mutant)) for (cid in ids) {
            rec <- view$cache[[view$position[cid]]]
            if (is.null(rec) || !state$cells[[cid]]$group_id %in% groups || !any(rec$s_rr != 0)) next
            wrong <- input; wrong$lambda[cid] <- wrong$lambda[cid] + 1
            effect <- max(abs(do.call(v53_scalar, wrong)$phi - ref$phi))
            if (is.finite(effect) && effect > 1e-8) {
              mutant <<- list(sample = input, cid = cid, effect = effect)
              break
            }
          }
        }
        out
      }
      operations[[length(operations) + 1L]] <<- op
    }
    op
  }, .package = "cdid")
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  out <- matrix(0, N, length(fit$links)); fn <- v53_ns(if (kind == "augmented") "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk")
  for (first in seq.int(1L, N, by = fit$control$max_chunk_rows)) {
    rows <- rev(seq.int(first, min(N, first + fit$control$max_chunk_rows - 1L)))
    out[rows, ] <- fn(fit, rows, 3:1)$phi
  }
  for (op in operations) testthat::expect_null(op$pullback(NULL, NULL, NULL, NULL, NULL, NULL, NULL))
  list(phi = out, issued = issued, calls = calls, used = used, suffix = suffix, samples = samples, mutant = mutant)
}

testthat::test_that("coarse scalar assembly preserves original solves/cache decisions and complete-F inference", {
  specs <- list(list(kind = "imputation", method = "reg", imp = "separate"),
    list(kind = "imputation", method = "dr", imp = "within_cohort"),
    list(kind = "imputation", method = "dr", imp = "across_cohorts"),
    list(kind = "augmented", imp = "within_cohort", model = "pair_logit"),
    list(kind = "augmented", imp = "across_cohorts", model = "independent_period_logit"))
  suffix_total <- 0L
  for (spec in specs) {
    a <- if (spec$kind == "augmented") v53_aug$rap_case() else v53_pure$rcp_case()
    fit <- if (spec$kind == "augmented") v53_aug$rap_fit(a, imp = spec$imp, sel = "within_cohort", model = spec$model, mode = "pooled") else
      v53_pure$rcp_fit(a, spec$method, spec$imp, "pooled")
    fast <- v53_run(fit, spec$kind); old <- v53_run(fit, spec$kind, FALSE)
    testthat::expect_gt(fast$used, 0L); testthat::expect_gt(fast$issued, 0L)
    testthat::expect_gt(fast$used, fast$issued)
    testthat::expect_lte(fast$issued, ceiling(nrow(a$Y) / fit$control$max_chunk_rows))
    testthat::expect_equal(fast$phi, old$phi, tolerance = 1e-12)
    testthat::expect_equal(crossprod(fast$phi), crossprod(old$phi), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(fast$phi))), sqrt(diag(crossprod(old$phi))), tolerance = 1e-12)
    testthat::expect_gt(max(abs(fast$phi[!a$A, , drop = FALSE])), 1e-4)
    suffix_total <- suffix_total + fast$suffix
    point <- if (spec$kind == "augmented") v53_aug$rap_point(fit) else v53_pure$rcp_point(fit)
    reference <- if (spec$kind == "augmented") v53_aug$rap_reference else v53_pure$rcp_reference
    offset <- length(point) - length(fit$links); N <- nrow(a$Y)
    for (i in c(1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      derivative <- (reference(a, fit, 1 + 1e-6 * delta) - reference(a, fit, 1 - 1e-6 * delta)) / 2e-6
      testthat::expect_equal(fast$phi[i, ], derivative[offset + 3:1], tolerance = 1e-5)
    }
  }
  testthat::expect_gt(suffix_total, 0L)
})

v53_native <- function(sample) {
  kernels <- v53_ns("reconstruction_coarse_kernels")()
  defaults <- get("defaults", environment(v53_ns("reconstruction_coarse_operation")), inherits = FALSE)
  current <- lapply(names(defaults), function(nm) get(nm, asNamespace("cdid"), inherits = TRUE))
  names(current) <- names(defaults)
  testthat::expect_true(all(vapply(names(defaults), function(nm) identical(current[[nm]], defaults[[nm]]), logical(1L))))
  functions <- c(current, list(.owner = asNamespace("cdid")))
  list(kernels = kernels, functions = functions,
    create = function(cap = 1e8, view = sample$view, state = sample$state)
      .Call(kernels[[1L]], view, state, sample$centered, functions, cap),
    call = function(ptr, x = sample) .Call(kernels[[2L]], ptr, x$view, x$state, x$lambda,
      x$nu, x$active, x$ids, x$templates),
    close = function(ptr) .Call(kernels[[2L]], ptr, NULL, NULL, NULL, NULL, NULL, NULL, NULL))
}

testthat::test_that("native coarse boundary keeps first-inactive arithmetic, active zeros, caps and refusals", {
  a <- v53_pure$rcp_case(); fit <- v53_pure$rcp_fit(a, "reg", "separate", "pooled")
  run <- v53_run(fit, "imputation"); testthat::expect_gt(run$used, 0L)
  testthat::expect_type(run$mutant, "list")
  testthat::expect_gt(run$mutant$effect, 1e-8)
  x <- run$mutant$sample; api <- v53_native(x); ptr <- api$create(); on.exit(api$close(ptr), add = TRUE)
  testthat::expect_type(ptr, "externalptr")
  # Exact helper semantics allow serialized closure copies; descriptor binding
  # must still use the actual current closure objects, not the stale copies.
  copied <- unserialize(serialize(api$functions$or_wls_gram_solve, NULL))
  testthat::expect_identical(copied, api$functions$or_wls_gram_solve)
  testthat::expect_false(v53_ns("gmm_sexp_aliases")(copied, api$functions$or_wls_gram_solve))
  stale_functions <- api$functions; stale_functions$or_wls_gram_solve <- copied
  testthat::expect_null(.Call(api$kernels[[1L]], x$view, x$state, x$centered, stale_functions, 1e8))
  fresh <- x; fresh$templates <- vector("list", length(x$state$groups))
  got <- api$call(ptr, fresh); ref <- do.call(v53_scalar, fresh)
  testthat::expect_equal(got, ref, tolerance = 1e-12)
  testthat::expect_identical(got$counts[1L], as.integer(sum(vapply(x$state$groups, `[[`, integer(1L), "p") > 0L)))
  zero <- fresh; zero$lambda[] <- 0; zero$nu <- lapply(zero$nu, function(v) v * 0)
  zero_result <- api$call(ptr, zero)
  testthat::expect_equal(zero_result, do.call(v53_scalar, zero), tolerance = 1e-12)
  testthat::expect_identical(zero_result$counts[1L], got$counts[1L])
  # A coefficient perturbation is detected by the original scalar oracle.
  wrong <- fresh; cid <- run$mutant$cid; wrong$lambda[cid] <- wrong$lambda[cid] + 1
  testthat::expect_gt(max(abs(api$call(ptr, wrong)$phi - ref$phi)), 1e-8)
  testthat::expect_equal(api$call(ptr, wrong)$phi, do.call(v53_scalar, wrong)$phi, tolerance = 1e-12)
  overflow <- fresh; overflow$lambda[cid] <- Inf
  actual <- tryCatch(api$call(ptr, overflow), error = identity)
  expected <- tryCatch(do.call(v53_scalar, overflow), error = identity)
  testthat::expect_s3_class(actual, "error"); testthat::expect_identical(class(actual), class(expected))
  testthat::expect_identical(conditionMessage(actual), conditionMessage(expected))
  changed <- x; changed$state$B <- changed$state$B + 1L; testthat::expect_null(api$call(ptr, changed))
  unseen <- x$view; unseen$seen[1L] <- FALSE; testthat::expect_null(api$create(view = unseen))
  decorated <- fresh; attr(decorated$nu[[cid]], "private") <- TRUE; testthat::expect_null(api$call(ptr, decorated))
  foreign <- setdiff(seq_along(x$state$cells), x$state$cids)
  if (length(foreign)) { bad <- x; bad$ids <- as.integer(foreign[1L]); testthat::expect_null(api$call(ptr, bad)) }
  peak <- max(vapply(x$state$groups, `[[`, integer(1L), "p")); sump <- sum(vapply(x$state$groups, `[[`, integer(1L), "p"))
  bound <- 2048 + 4096 + 4 * length(x$state$cells) + 128 * length(x$state$groups) + 4 * sump +
    32 * peak^2 + x$state$B * (4 * peak + 8)
  testthat::expect_null(api$create(bound - 1)); exact <- api$create(bound)
  testthat::expect_type(exact, "externalptr"); testthat::expect_equal(api$call(exact), do.call(v53_scalar, x), tolerance = 1e-12)
  api$close(exact); testthat::expect_null(api$call(exact)); testthat::expect_silent(api$close(exact))
  badmass <- x$state; badmass$cells[[cid]]$mass <- NA_integer_; testthat::expect_null(api$create(state = badmass))
})

testthat::test_that("p0 no-intercept and already-fitted IPW paths retain their scalar contracts", {
  for (spec in list(c(TRUE, FALSE), c(FALSE, TRUE))) {
    a <- v53_pure$rcp_case(intercept = spec[1L], slopes = spec[2L])
    fit <- v53_pure$rcp_fit(a, "reg", "separate", "pooled")
    fast <- v53_run(fit, "imputation"); old <- v53_run(fit, "imputation", FALSE)
    if (spec[1L]) testthat::expect_gt(fast$used, 0L) else {
      # Pooled, one-slope, no-intercept OR has only ONE scalar beta-j
      # loading, before the first full traversal creates its witness. There
      # is no later eligible loading: retain the original scalar contract.
      testthat::expect_false(fit$outcome_recipe$intercept)
      testthat::expect_identical(fit$shared_or$p, 1L)
      testthat::expect_identical(fast$issued, 0L)
      testthat::expect_identical(fast$used, 0L)
    }
    testthat::expect_equal(fast$phi, old$phi, tolerance = 1e-12)
  }
  a <- v53_pure$rcp_case(); fit <- v53_pure$rcp_fit(a, "ipw", "within_cohort", "pooled")
  fast <- v53_run(fit, "imputation"); old <- v53_run(fit, "imputation", FALSE)
  testthat::expect_identical(fast$issued, 0L)
  testthat::expect_identical(fast$phi, old$phi)
  # Instrumenting a private helper does not disable the prepared operation.
  fit <- v53_pure$rcp_fit(a, "reg", "separate", "pooled")
  helper <- v53_ns("reconstruction_zero_template_run"); calls <- 0L
  testthat::local_mocked_bindings(reconstruction_zero_template_run = function(...) { calls <<- calls + 1L; helper(...) }, .package = "cdid")
  custom <- v53_run(fit, "imputation"); reference <- v53_run(fit, "imputation", FALSE)
  testthat::expect_gt(custom$used, 0L); testthat::expect_gt(calls, 0L)
  testthat::expect_equal(custom$phi, reference$phi, tolerance = 1e-12)
})
