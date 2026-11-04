# Use the installed package namespace; no route fixtures or extra dependencies.
outcome_blocks_private <- function(name) getFromNamespace(name, "cdid")

outcome_blocks_recipe <- function(intercept = TRUE, slopes = "x") {
  outcome_blocks_private("or_pooled_outcome_recipe")(
    c(if (intercept) "(Intercept)", slopes), intercept, slopes)
}

outcome_blocks_rows <- function(recipe, z, y, w = rep(1, length(y)), L = 1, scale = 1) {
  z <- matrix(z, nrow = length(y), ncol = recipe$p)
  Z <- matrix(0, length(y), length(recipe$columns), dimnames = list(NULL, recipe$columns))
  if (recipe$intercept) Z[, recipe$intercept_index] <- 1
  if (recipe$p) Z[, recipe$slope_index] <- z
  acc <- outcome_blocks_private("or_pooled_support_start")(recipe, L)
  acc <- outcome_blocks_private("or_pooled_support_add")(acc, recipe, Z, y, w)
  acc$scale <- scale
  list(Z = Z, y = y, w = w, scale = scale, support = acc)
}

outcome_blocks_case <- function(recipe, rows, blocks, links) {
  visits <- new.env(parent = emptyenv()); visits$ids <- integer()
  visit_rows <- function(sid, visit) {
    visits$ids <- c(visits$ids, as.integer(sid))
    Z <- rows[[sid]]$Z
    # Exercise bounded chunk replay, including a one-row final chunk.
    for (start in seq.int(1L, nrow(Z), by = 2L))
      visit(Z[seq.int(start, min(start + 1L, nrow(Z))), , drop = FALSE])
    invisible(NULL)
  }
  list(recipe = recipe, rows = rows, supports = lapply(rows, `[[`, "support"),
    blocks = blocks, links = links, visit_rows = visit_rows, visits = visits)
}

outcome_blocks_fit <- function(case, mode, N = 97L) {
  outcome_blocks_private("or_outcome_fit_blocks")(case$recipe, case$supports,
    case$blocks, case$links, 1e-7, N, case$visit_rows, mode)
}

# Independent dense reference used only on tiny fixtures. One indicator per
# selected block plus common L*z slopes constructs the pooled joint objective;
# selecting a singleton constructs its independent control WLS reference.
outcome_blocks_reference <- function(case, bids) {
  recipe <- case$recipe
  pieces <- lapply(seq_along(bids), function(j) {
    b <- case$blocks[[bids[j]]]; rows <- case$rows[b$support_ids]
    Z <- do.call(rbind, lapply(rows, `[[`, "Z"))
    d <- b$L * Z[, recipe$slope_index, drop = FALSE]
    if (recipe$intercept) {
      indicators <- matrix(0, nrow(Z), length(bids)); indicators[, j] <- 1
      d <- cbind(indicators, d)
    }
    list(X = d, y = unlist(lapply(rows, `[[`, "y"), use.names = FALSE),
      w = unlist(lapply(rows, function(r) r$w * r$scale), use.names = FALSE))
  })
  X <- do.call(rbind, lapply(pieces, `[[`, "X"))
  y <- unlist(lapply(pieces, `[[`, "y"), use.names = FALSE)
  w <- unlist(lapply(pieces, `[[`, "w"), use.names = FALSE)
  fit <- stats::lm.wfit(X, y, w, tol = 1e-7)
  list(coefficients = unname(fit$coefficients), rss = sum(w * fit$residuals^2),
    X = X, y = y, w = w)
}

outcome_blocks_heterogeneous <- function(intercept = TRUE) {
  recipe <- outcome_blocks_recipe(intercept)
  rows <- list(
    outcome_blocks_rows(recipe, c(-2, 0, 1), c(-5, 1, 5), c(1, 2, 1), L = 2, scale = 3),
    outcome_blocks_rows(recipe, c(-1, 1, 3, 4), c(9, 3, -2, -6), c(2, 1, 3, 1), scale = 0.4),
    outcome_blocks_rows(recipe, c(2, 3), c(8, 10), c(2, 1), L = 2, scale = 0.2),
    outcome_blocks_rows(recipe, c(0.5, 1.5), c(0, 0), L = 2, scale = 123),
    outcome_blocks_rows(recipe, c(0, 2), c(0, 0), scale = 456))
  outcome_blocks_case(recipe, rows,
    list(list(name = "long", support_ids = c(3L, 1L), L = 2),
      list(name = "short", support_ids = 2L, L = 1)),
    list(list(regression = 2L, support_id = 5L), list(regression = 1L, support_id = 4L)))
}

outcome_blocks_alias_case <- function(duplicate = FALSE, failing = FALSE) {
  recipe <- outcome_blocks_recipe(TRUE, c("x", "alias")); x <- c(-2, -1, 0, 1, 2)
  rows <- list(
    outcome_blocks_rows(recipe, cbind(x, 4 - x), 2 + 3 * x, L = 1),
    outcome_blocks_rows(recipe, cbind(c(-1, 2), c(1, 7)), c(0, 0), L = 2),
    # Entirely unused support: including this scale in alias checks would hide
    # the failing required prediction below, besides widening the row domain.
    outcome_blocks_rows(recipe, cbind(1e14, -1e14), 0),
    outcome_blocks_rows(recipe, cbind(x, 3 + 2 * x), 1 + 4 * x, L = 2),
    outcome_blocks_rows(recipe, cbind(c(-1, 2), c(5, 2)), c(0, 0)),
    outcome_blocks_rows(recipe, cbind(1, if (failing) 3.1 else 3), 0))
  links <- list(list(regression = 2L, support_id = 5L),
    list(regression = 1L, support_id = 2L), list(regression = 2L, support_id = 6L))
  if (duplicate) links <- c(links, links[c(2L, 1L, 2L)])
  outcome_blocks_case(recipe, rows,
    list(list(name = "first", support_ids = 4L, L = 2),
      list(name = "second", support_ids = 1L, L = 1)), links)
}

# separate fits versus independent per-block WLS, scales and full N.
testthat::test_that("separate blocks retain heterogeneous control slopes and normalized operators", {
  case <- outcome_blocks_heterogeneous(); fit <- outcome_blocks_fit(case, "separate")
  block_fit <- outcome_blocks_private("or_outcome_block_fit")
  predict <- outcome_blocks_private("or_outcome_block_predict")
  gram <- outcome_blocks_private("or_outcome_block_gram_solve")
  testthat::expect_identical(fit$block_map$group_id, 1:2)
  for (b in 1:2) {
    ref <- outcome_blocks_reference(case, b); f <- block_fit(fit, b)
    testthat::expect_equal(c(f$alpha, f$beta), ref$coefficients, tolerance = 1e-10)
    testthat::expect_equal(f$weight_sum, sum(ref$w), tolerance = 1e-10)
    testthat::expect_equal(f$mass, sum(ref$w) / 97, tolerance = 1e-10)
    testthat::expect_equal(f$mean_d, sum(ref$w * ref$X[, 2L]) / sum(ref$w), tolerance = 1e-10)
    testthat::expect_equal(f$mean_y, sum(ref$w * ref$y) / sum(ref$w), tolerance = 1e-10)
    dc <- ref$X[, 2L] - sum(ref$w * ref$X[, 2L]) / sum(ref$w)
    testthat::expect_equal(gram(fit, b, 2), 97 * 2 / sum(ref$w * dc^2), tolerance = 1e-10)
    sid <- case$links[[match(b, vapply(case$links, `[[`, integer(1), "regression"))]]$support_id
    Z <- case$rows[[sid]]$Z
    expected <- as.vector(cbind(1, case$blocks[[b]]$L * Z[, 2L]) %*% ref$coefficients)
    testthat::expect_equal(predict(fit, b, Z), expected, tolerance = 1e-10)
  }
  testthat::expect_gt(abs(block_fit(fit, 1)$beta - block_fit(fit, 2)$beta), 1)
  refs <- lapply(1:2, function(b) outcome_blocks_reference(case, b))
  testthat::expect_equal(fit$objective, sum(vapply(refs, `[[`, numeric(1), "rss")) / 97,
    tolerance = 1e-10)
  twice <- outcome_blocks_fit(case, "separate", N = 194L)
  testthat::expect_equal(gram(twice, 1, 2), 2 * gram(fit, 1, 2), tolerance = 1e-10)
  testthat::expect_equal(twice$objective, fit$objective / 2, tolerance = 1e-10)
  testthat::expect_equal(block_fit(twice, 1)$mass, block_fit(fit, 1)$mass / 2, tolerance = 1e-10)
})

# independently assembled pooled joint fit and exact old-helper result.
testthat::test_that("pooled mode delegates the unchanged joint fit", {
  case <- outcome_blocks_heterogeneous(); fit <- outcome_blocks_fit(case, "pooled")
  ref <- outcome_blocks_reference(case, 1:2)
  old <- outcome_blocks_private("or_pooled_fit_blocks")(case$recipe, case$supports,
    case$blocks, case$links, 1e-7, 97L, case$visit_rows)
  testthat::expect_identical(fit$groups[[1L]], old)
  testthat::expect_identical(fit$block_map$group_id, c(1L, 1L))
  testthat::expect_identical(fit$block_map$local_block_id, 1:2)
  alpha <- vapply(fit$groups[[1L]]$regressions, `[[`, numeric(1), "alpha")
  testthat::expect_equal(c(alpha, fit$groups[[1L]]$beta), ref$coefficients, tolerance = 1e-10)
  testthat::expect_equal(fit$rss, ref$rss, tolerance = 1e-10)
  testthat::expect_equal(fit$objective, ref$rss / 97, tolerance = 1e-10)
  f <- outcome_blocks_private("or_outcome_block_fit")(fit, "short")
  testthat::expect_identical(f$group_id, 1L)
  testthat::expect_identical(f$local_block_id, 2L)
  testthat::expect_equal(f$alpha, ref$coefficients[2L], tolerance = 1e-10)
})

# multiple retained coordinates, matrix RHS and a non-leading intercept.
testthat::test_that("block accessors preserve recipe order and multi-column operators", {
  recipe <- outcome_blocks_private("or_pooled_outcome_recipe")(
    c("x", "(Intercept)", "v"), TRUE, c("v", "x"))
  x <- c(-2, -1, 0, 1, 2); v <- c(1, 0, 2, -1, 1)
  rows <- list(
    outcome_blocks_rows(recipe, cbind(x, v), 2 + 3 * x - v, c(1, 2, 1, 3, 1), L = 2),
    outcome_blocks_rows(recipe, cbind(x, v), 4 - x + 2 * v, c(2, 1, 2, 1, 3), scale = 0.3))
  case <- outcome_blocks_case(recipe, rows,
    list(list(name = "left", support_ids = 1L, L = 2),
      list(name = "right", support_ids = 2L, L = 1)), list())
  rhs <- matrix(c(1, 2, -1, 3), 2L, 2L)
  for (mode in c("separate", "pooled")) {
    fit <- outcome_blocks_fit(case, mode)
    for (b in 1:2) {
      bids <- if (mode == "pooled") 1:2 else b
      ref <- outcome_blocks_reference(case, bids)
      f <- outcome_blocks_private("or_outcome_block_fit")(fit, b)
      testthat::expect_identical(f$keep, 1:2)
      testthat::expect_identical(f$basis$kept_columns, c("x", "v"))
      testthat::expect_equal(f$beta, tail(ref$coefficients, 2L), tolerance = 1e-10)
      alpha <- ref$coefficients[match(b, bids)]
      d <- case$blocks[[b]]$L * rows[[b]]$Z[, c(1L, 3L), drop = FALSE]
      testthat::expect_equal(outcome_blocks_private("or_outcome_block_predict")(fit, b, rows[[b]]$Z),
        as.vector(alpha + d %*% tail(ref$coefficients, 2L)), tolerance = 1e-10)
      centered <- ref$X[, length(bids) + seq_len(2L), drop = FALSE]
      for (j in seq_along(bids)) {
        ix <- which(ref$X[, j] == 1)
        center <- colSums(ref$w[ix] * centered[ix, , drop = FALSE]) / sum(ref$w[ix])
        centered[ix, ] <- sweep(centered[ix, , drop = FALSE], 2L, center)
      }
      expected <- unname(97 * solve(crossprod(centered, ref$w * centered), rhs))
      testthat::expect_equal(outcome_blocks_private("or_outcome_block_gram_solve")(fit, b, rhs),
        expected, tolerance = 1e-10)
    }
  }
})

# no-intercept separate and pooled numeric references and Gram scaling.
testthat::test_that("no-intercept recipes keep uncentered fits and add no alpha", {
  case <- outcome_blocks_heterogeneous(FALSE)
  block_fit <- outcome_blocks_private("or_outcome_block_fit")
  for (mode in c("separate", "pooled")) {
    fit <- outcome_blocks_fit(case, mode)
    for (b in 1:2) {
      ref <- outcome_blocks_reference(case, if (mode == "pooled") 1:2 else b)
      f <- block_fit(fit, b)
      testthat::expect_identical(f$alpha, 0)
      testthat::expect_equal(f$beta, ref$coefficients, tolerance = 1e-10)
      testthat::expect_identical(f$p, 1L)
      testthat::expect_identical(f$coefficient_dimension, 1L)
      testthat::expect_equal(outcome_blocks_private("or_outcome_block_gram_solve")(fit, b, 1),
        97 / sum(ref$w * ref$X[, 1L]^2), tolerance = 1e-10)
      Z <- case$rows[[case$blocks[[b]]$support_ids[1L]]]$Z
      testthat::expect_equal(outcome_blocks_private("or_outcome_block_predict")(fit, b, Z),
        as.vector(case$blocks[[b]]$L * Z %*% ref$coefficients), tolerance = 1e-10)
    }
  }
})

# intercept-only keeps zero-width slopes and the existing mean solve.
testthat::test_that("intercept-only fits preserve means and zero-width operators", {
  recipe <- outcome_blocks_recipe(TRUE, character())
  rows <- list(outcome_blocks_rows(recipe, numeric(), c(2, 4, 8), c(1, 2, 1), scale = 3),
    outcome_blocks_rows(recipe, numeric(), c(5, 7), c(2, 1), L = 3, scale = 0.5))
  case <- outcome_blocks_case(recipe, rows,
    list(list(name = "a", support_ids = 1L, L = 1),
      list(name = "b", support_ids = 2L, L = 3)), list())
  for (mode in c("separate", "pooled")) {
    fit <- outcome_blocks_fit(case, mode)
    for (b in 1:2) {
      ref <- outcome_blocks_reference(case, b)
      f <- outcome_blocks_private("or_outcome_block_fit")(fit, b)
      testthat::expect_identical(f$p, 0L)
      testthat::expect_length(f$beta, 0L)
      testthat::expect_equal(f$alpha, ref$coefficients[1L], tolerance = 1e-10)
      testthat::expect_equal(outcome_blocks_private("or_outcome_block_predict")(fit, b, rows[[b]]$Z),
        rep(ref$coefficients[1L], nrow(rows[[b]]$Z)), tolerance = 1e-10)
      testthat::expect_equal(outcome_blocks_private("or_outcome_block_gram_solve")(
        fit, b, matrix(numeric(), 0L, 3L)), matrix(numeric(), 0L, 3L))
      testthat::expect_length(outcome_blocks_private("or_outcome_block_predict")(
        fit, b, rows[[b]]$Z[FALSE, , drop = FALSE]), 0L)
    }
  }
})

# modeled-calendar duration, rather than numeric date-label gaps.
testthat::test_that("calendar durations scale the fitted design and predictions", {
  recipe <- outcome_blocks_recipe(); x <- c(-2, -1, 1, 3)
  pair <- list(s_calendar_index = 2L, t_calendar_index = 5L,
    s = 1990, t = 2070, time_gap = 80)
  L <- outcome_blocks_private("or_pooled_transitions")(pair)
  testthat::expect_identical(L, 3)
  rows <- list(outcome_blocks_rows(recipe, x, 4 + 6 * x, L = L),
    outcome_blocks_rows(recipe, x, 2 - 4 * x, L = 2))
  case <- outcome_blocks_case(recipe, rows,
    list(list(name = "irregular", support_ids = 1L, L = L),
      list(name = "two", support_ids = 2L, L = 2)), list())
  fit <- outcome_blocks_fit(case, "separate")
  testthat::expect_equal(outcome_blocks_private("or_outcome_block_fit")(fit, 1)$beta, 2,
    tolerance = 1e-10)
  testthat::expect_equal(outcome_blocks_private("or_outcome_block_fit")(fit, 2)$beta, -2,
    tolerance = 1e-10)
  for (b in 1:2) testthat::expect_equal(
    outcome_blocks_private("or_outcome_block_predict")(fit, b, rows[[b]]$Z), rows[[b]]$y,
    tolerance = 1e-10)
})

# local aliases, original support/link IDs and required row domains.
testthat::test_that("separate alias checks use only local controls and consuming treated rows", {
  case <- outcome_blocks_alias_case(); fit <- outcome_blocks_fit(case, "separate")
  testthat::expect_identical(fit$group_map[[1L]]$support_ids, c(2L, 4L))
  testthat::expect_identical(fit$group_map[[2L]]$support_ids, c(1L, 5L, 6L))
  testthat::expect_identical(fit$group_map[[1L]]$link_ids, 2L)
  testthat::expect_identical(fit$group_map[[2L]]$link_ids, c(1L, 3L))
  testthat::expect_identical(case$visits$ids, c(2L, 4L, 1L, 5L, 6L))
  for (b in 1:2) {
    f <- outcome_blocks_private("or_outcome_block_fit")(fit, b)
    testthat::expect_identical(f$keep, 1L)
    testthat::expect_identical(f$basis$aliased_columns, "alias")
    testthat::expect_identical(rownames(f$basis$block_constants), case$blocks[[b]]$name)
    sid <- case$blocks[[b]]$support_ids
    r <- case$rows[[sid]]; X <- cbind(1, case$blocks[[b]]$L * r$Z[, 2L])
    ref <- stats::lm.wfit(X, r$y, r$w * r$scale, tol = 1e-7)
    testthat::expect_equal(c(f$alpha, f$beta), unname(ref$coefficients), tolerance = 1e-10)
    for (l in case$links) if (l$regression == b) {
      Z <- case$rows[[l$support_id]]$Z
      testthat::expect_equal(outcome_blocks_private("or_outcome_block_predict")(fit, b, Z),
        as.vector(cbind(1, case$blocks[[b]]$L * Z[, 2L]) %*% ref$coefficients), tolerance = 1e-10)
    }
  }
  testthat::expect_equal(fit$groups[[1L]]$basis$block_constants[1L, 1L], 3, tolerance = 1e-10)
  testthat::expect_equal(fit$groups[[2L]]$basis$block_constants[1L, 1L], 4, tolerance = 1e-10)
})

# a later link of original block 2 must not be lost in local remapping.
testthat::test_that("every consuming prediction domain must preserve a local exact alias", {
  case <- outcome_blocks_alias_case(failing = TRUE)
  testthat::expect_error(outcome_blocks_fit(case, "separate"),
    "required treated predictions", class = "cdid_or_wls_rank_error")
  testthat::expect_identical(case$visits$ids, c(2L, 4L, 1L, 5L, 6L))
})

# link duplication changes the mapping only, not fits or replay work.
testthat::test_that("duplicated consuming links deterministically reuse the caller block", {
  one <- outcome_blocks_alias_case(); many <- outcome_blocks_alias_case(duplicate = TRUE)
  fit_one <- outcome_blocks_fit(one, "separate"); fit_many <- outcome_blocks_fit(many, "separate")
  testthat::expect_identical(fit_many$groups, fit_one$groups)
  testthat::expect_identical(fit_many$block_map, fit_one$block_map)
  testthat::expect_identical(many$visits$ids, one$visits$ids)
  testthat::expect_identical(fit_many$group_map[[1L]]$link_ids, c(2L, 4L, 6L))
  testthat::expect_identical(fit_many$group_map[[2L]]$link_ids, c(1L, 3L, 5L))
  testthat::expect_identical(fit_many$retained, fit_one$retained)
})

# equal fitted numbers are insufficient to merge distinct source fits.
testthat::test_that("distinct caller blocks remain separate when their fitted numbers coincide", {
  recipe <- outcome_blocks_recipe(); r <- outcome_blocks_rows(recipe, c(-1, 0, 2), c(1, 3, 7))
  case <- outcome_blocks_case(recipe, list(r, r),
    list(list(name = "context-a", support_ids = 1L, L = 1),
      list(name = "context-b", support_ids = 2L, L = 1)), list())
  fit <- outcome_blocks_fit(case, "separate")
  testthat::expect_length(fit$groups, 2L)
  testthat::expect_identical(fit$block_map$group_id, 1:2)
  testthat::expect_equal(fit$groups[[1L]]$beta, fit$groups[[2L]]$beta, tolerance = 1e-10)
})

# no-intercept exact aliases keep a zero constant and refuse violations.
testthat::test_that("no-intercept local aliases require prediction support without a constant", {
  recipe <- outcome_blocks_recipe(FALSE, c("x", "alias")); x <- c(-2, -1, 1, 2)
  control <- outcome_blocks_rows(recipe, cbind(x, 2 * x), 3 * x, L = 2)
  treated <- outcome_blocks_rows(recipe, cbind(c(1, 3), c(2, 6)), c(0, 0), L = 2)
  blocks <- list(list(name = "uncentered", support_ids = 2L, L = 2))
  links <- list(list(regression = 1L, support_id = 1L))
  fit <- outcome_blocks_fit(outcome_blocks_case(recipe, list(treated, control), blocks, links), "separate")
  f <- outcome_blocks_private("or_outcome_block_fit")(fit, 1)
  ref <- stats::lm.wfit(matrix(2 * x, ncol = 1L), 3 * x, rep(1, 4), tol = 1e-7)
  testthat::expect_equal(f$beta, unname(ref$coefficients), tolerance = 1e-10)
  testthat::expect_equal(f$basis$block_constants[1L, 1L], 0)
  testthat::expect_equal(outcome_blocks_private("or_outcome_block_predict")(fit, 1, treated$Z),
    c(3, 9), tolerance = 1e-10)
  bad <- outcome_blocks_rows(recipe, cbind(c(1, 3), c(3, 7)), c(0, 0), L = 2)
  testthat::expect_error(outcome_blocks_fit(
    outcome_blocks_case(recipe, list(bad, control), blocks, links), "separate"),
    "required treated predictions", class = "cdid_or_wls_rank_error")
})

# block-local dimensions may differ, including a zero-slope mean fit.
testthat::test_that("local retained dimensions can differ across separate blocks", {
  recipe <- outcome_blocks_recipe()
  rows <- list(outcome_blocks_rows(recipe, c(2, 2, 2), c(1, 3, 5)),
    outcome_blocks_rows(recipe, c(-1, 0, 2), c(0, 2, 5)))
  case <- outcome_blocks_case(recipe, rows,
    list(list(name = "constant", support_ids = 1L, L = 1),
      list(name = "varying", support_ids = 2L, L = 1)), list())
  fit <- outcome_blocks_fit(case, "separate")
  p <- outcome_blocks_private("or_outcome_block_p")
  testthat::expect_identical(c(p(fit, 1), p(fit, 2)), c(0L, 1L))
  testthat::expect_equal(outcome_blocks_private("or_outcome_block_predict")(fit, 1, rows[[1L]]$Z),
    c(3, 3, 3), tolerance = 1e-10)
  testthat::expect_identical(fit$retained[[1L]]$operator_dimension, c(0L, 0L))
  testthat::expect_identical(fit$retained[[2L]]$operator_dimension, c(1L, 1L))
})

# retain existing rank/roundoff policy, without a new residual-df rule.
testthat::test_that("exactly determined fits remain allowed and near aliases remain subject to rank checks", {
  recipe <- outcome_blocks_recipe()
  rows <- list(outcome_blocks_rows(recipe, c(0, 1), c(2, 5)))
  blocks <- list(list(name = "zero-df", support_ids = 1L, L = 1))
  fit <- outcome_blocks_fit(outcome_blocks_case(recipe, rows, blocks, list()), "separate")
  testthat::expect_equal(c(outcome_blocks_private("or_outcome_block_fit")(fit, 1)$alpha,
    fit$groups[[1L]]$beta), c(2, 3), tolerance = 1e-10)
  recipe <- outcome_blocks_recipe(TRUE, c("x", "near")); x <- c(-2, -1, 0, 1, 2)
  rows <- list(outcome_blocks_rows(recipe, cbind(x, x + 1e-9 * x^2), c(1, 2, 1, 3, 4)))
  testthat::expect_error(outcome_blocks_fit(outcome_blocks_case(recipe, rows, blocks, list()), "separate"),
    "Final pooled WLS slope rank", class = "cdid_or_wls_rank_error")
  recipe <- outcome_blocks_recipe(FALSE)
  rows <- list(outcome_blocks_rows(recipe, c(0, 0), c(1, 2)))
  testthat::expect_error(outcome_blocks_fit(outcome_blocks_case(recipe, rows, blocks, list()), "separate"),
    "no retained columns", class = "cdid_or_wls_rank_error")
})

# metadata preflight describes dimensions without fitting or scanning.
testthat::test_that("preflight retains original IDs and rejects invalid local-map inputs", {
  case <- outcome_blocks_alias_case()
  layout <- outcome_blocks_private("or_outcome_group_layout")
  plan <- layout(case$recipe, case$supports, case$blocks, case$links, "separate")
  testthat::expect_identical(case$visits$ids, integer())
  testthat::expect_identical(plan$group_map[[2L]]$block_ids, 2L)
  testthat::expect_equal(plan$group_map[[2L]]$slope_dimension_bound, 2)
  testthat::expect_equal(plan$group_map[[2L]]$coefficient_dimension_bound, 3)
  testthat::expect_error(layout(case$recipe, case$supports, case$blocks, case$links, "sep"),
    "explicit private outcome mode")
  bad <- case$links; bad[[1L]]$regression <- 3L
  testthat::expect_error(layout(case$recipe, case$supports, case$blocks, bad, "separate"),
    "original block ID")
})
