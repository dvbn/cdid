v64_ns <- function(nm) getFromNamespace(nm, "cdid")
v64_e <- new.env(parent = asNamespace("cdid"))
for (x in parse(file = testthat::test_path("test-augmented-adjoint-replay.R")))
  if (is.call(x) && identical(x[[1L]], as.name("<-"))) eval(x, v64_e)

# Restore the literal original Stage3 loop by deleting ONLY the four new
# reuse statements. Every original fit/check/forcing expression is retained.
v64_original_core <- function() {
  f <- v64_ns("imputation_fit_core"); removed <- 0L
  rewrite <- function(x) {
    if (!is.call(x)) return(x)
    if (identical(x[[1L]], as.name("{"))) {
      out <- list(x[[1L]])
      for (item in as.list(x)[-1L]) {
        if (is.call(item) && identical(item[[1L]], as.name("<-")) &&
            identical(item[[2L]], as.name("numeric_fields"))) { removed <<- removed + 1L; next }
        out[length(out) + 1L] <- list(rewrite(item))
      }
      return(as.call(out))
    }
    if (identical(x[[1L]], as.name("for")) && identical(x[[2L]], as.name("tid")) &&
        identical(x[[3L]], quote(seq_len(n_trt)))) {
      items <- as.list(x[[4L]])
      testthat::expect_identical(items[[3L]], quote(representative <- 0L))
      testthat::expect_identical(items[[4L]][[1L]], as.name("if"))
      testthat::expect_identical(items[[4L]][[2L]], quote(tid > 1L))
      testthat::expect_identical(items[[5L]][[2L]], as.name("representative"))
      x[[4L]] <- as.call(items[-c(3L, 4L, 5L)]); removed <<- removed + 3L
      return(x)
    }
    for (i in seq_along(x)[-1L]) if (!identical(x[[i]], quote(expr = ))) x[i] <- list(rewrite(x[[i]]))
    x
  }
  body(f) <- rewrite(body(f)); testthat::expect_identical(removed, 4L); f
}

v64_fit_count <- function(a, model, original = NULL) {
  logit <- v64_ns("cp04_logit"); identities <- list()
  testthat::local_mocked_bindings(cp04_logit = function(X, y, w, glm_control, N, identity) {
    if (is.list(identity) && !is.null(identity$target_cohort_index) && !is.null(identity$all_cohort_indices))
      identities[[length(identities) + 1L]] <<- identity
    logit(X, y, w, glm_control, N, identity)
  }, .package = "cdid")
  if (!is.null(original)) testthat::local_mocked_bindings(imputation_fit_core = original, .package = "cdid")
  fit <- v64_e$v57_aug$rap_fit(a, imp = "separate", model = model, mode = "pooled")
  list(fit = fit, identities = identities)
}

testthat::test_that("identical treatment data fit once while all original models and full-F IF remain", {
  original <- v64_original_core()
  for (model in c("pair_logit", "independent_period_logit")) {
    a <- v64_e$v57_wide_case()
    old <- v64_fit_count(a, model, original); actual <- v64_fit_count(a, model)
    keys <- lapply(actual$fit$treatments, function(tf)
      list(target = tf$target_cohort_index, risk = tf$all_cohort_indices))
    testthat::expect_identical(length(old$identities), length(old$fit$treatments))
    testthat::expect_identical(length(actual$identities), length(unique(keys)))
    testthat::expect_lt(length(actual$identities), length(old$identities))
    testthat::expect_identical(actual$identities[[1L]], old$identities[[1L]])
    testthat::expect_identical(length(actual$fit$treatments), length(old$fit$treatments))
    testthat::expect_equal(actual$fit$treatments, old$fit$treatments, tolerance = 1e-12)
    testthat::expect_equal(v64_e$v57_aug$rap_point(actual$fit), v64_e$v57_aug$rap_point(old$fit), tolerance = 1e-12)
    links <- rev(seq_along(actual$fit$links))
    raw <- v64_e$v57_collect(old$fit, links); phi <- v64_e$v57_collect(actual$fit, links)
    testthat::expect_equal(phi, raw, tolerance = 1e-12)
    testthat::expect_identical(phi == 0, raw == 0)
    testthat::expect_equal(crossprod(phi), crossprod(raw), tolerance = 1e-12)
    testthat::expect_equal(sqrt(diag(crossprod(phi))), sqrt(diag(crossprod(raw))), tolerance = 1e-12)
    testthat::expect_gt(max(abs(phi[!a$A, , drop = FALSE])), 1e-4)
    N <- nrow(a$Y); offset <- length(v64_e$v57_aug$rap_point(actual$fit)) - length(links)
    for (i in c(1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      ref <- (v64_e$v57_aug$rap_reference(a, actual$fit, 1 + 1e-6 * delta) -
        v64_e$v57_aug$rap_reference(a, actual$fit, 1 - 1e-6 * delta)) / 2e-6
      testthat::expect_equal(phi[i, ], ref[offset + links], tolerance = 1e-5)
    }
  }
})
