v47_ns <- function(name) getFromNamespace(name, "cdid")
v47_fixture <- local({
  e <- new.env(parent = environment())
  wanted <- c("rap_ns", "rap_gc", "rap_case", "rap_reference", "rap_point", "rap_fit")
  for (x in parse(file = testthat::test_path("test-augmented-outcome-pooling.R"))) {
    if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
        as.character(x[[2L]]) %in% wanted) eval(x, e)
  }
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE)))
  e
})
v47_rebind <- function(env, name, value) {
  old <- get(name, env, inherits = FALSE)
  unlockBinding(name, env); assign(name, value, env); lockBinding(name, env)
  function() { unlockBinding(name, env); assign(name, old, env); lockBinding(name, env) }
}
v47_fit <- function(a, spec = list(), default = FALSE, fail = NULL, clip = FALSE, pair_outcomes = NULL) {
  factory <- v47_ns("selection_prediction_operation"); cell <- v47_ns("augmented_cell_rows")
  predict <- v47_ns("predict_pair_selection"); derivative <- v47_ns("selection_pair_cohort_log_deriv")
  pure <- v47_ns("imputation_fit_core"); selector <- v47_ns("fit_selection_models")
  ops <- list(); owned <- active <- logical(); pairs <- derivs <- finishes <- aborts <- integer()
  expected <- list(); probabilities <- numeric(); clipped_derivatives <- 0L
  predicts <- derivatives <- pure_calls <- or_calls <- 0L; restore <- NULL
  on.exit(if (!is.null(restore)) restore(), add = TRUE)
  testthat::local_mocked_bindings(
    predict_pair_selection = function(...) { predicts <<- predicts + 1L; predict(...) },
    selection_pair_cohort_log_deriv = function(...) { derivatives <<- derivatives + 1L; derivative(...) },
    selection_prediction_operation = function(backend, expected_sources = NULL, source_reader = NULL, source_reference = NULL) {
      op <- factory(backend, expected_sources, source_reader, source_reference); j <- length(ops) + 1L
      owned[j] <<- !is.null(source_reader); active[j] <<- TRUE
      pairs[j] <<- 0L; derivs[j] <<- 0L; finishes[j] <<- 0L; aborts[j] <<- 0L
      expected[j] <<- list(expected_sources)
      pair <- op$pair; deriv <- op$deriv; finish <- op$finish; abort <- op$abort
      op$pair <- function(...) {
        out <- pair(...); pairs[j] <<- pairs[j] + 1L
        probabilities <<- c(probabilities, out$prob); out
      }
      op$deriv <- function(...) {
        out <- deriv(...); derivs[j] <<- derivs[j] + 1L
        probabilities <<- c(probabilities, out$pred$prob)
        low <- which(out$pred$prob <= 1e-6)
        if (length(low) && all(out$T_logit[low, , drop = FALSE] == 0))
          clipped_derivatives <<- clipped_derivatives + 1L
        out
      }
      op$finish <- function() { finishes[j] <<- finishes[j] + 1L; value <- finish(); active[j] <<- FALSE; value }
      op$abort <- function() { aborts[j] <<- aborts[j] + 1L; active[j] <<- FALSE; abort() }
      ops[[j]] <<- op; op
    },
    imputation_fit_core = function(...) {
      pure_calls <<- pure_calls + 1L
      testthat::expect_false(any(active[owned]))
      if (identical(fail, "glm_warning")) warning("v47 original imputation warning", call. = FALSE)
      pure(...)
    }, .package = "cdid")
  outcome <- v47_ns("or_outcome_fit_blocks"); pooled <- v47_ns("or_pooled_fit_blocks")
  visit_fitter <- function(original, ...) {
    or_calls <<- or_calls + 1L
    testthat::expect_false(any(active[owned]))
    if (identical(fail, "or_error") && any(owned)) stop("v47 original OR error", call. = FALSE)
    out <- original(...)
    if (identical(fail, "or_source") && any(owned) && is.null(restore)) {
      frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
      other <- frame; other$unit_weights_all[1L] <- other$unit_weights_all[1L] * 2
      restore <<- v47_rebind(a$reader$frame_ref, "frame_owner", other)
    }
    out
  }
  testthat::local_mocked_bindings(
    or_outcome_fit_blocks = function(...) visit_fitter(outcome, ...),
    or_pooled_fit_blocks = function(...) visit_fitter(pooled, ...), .package = "cdid")
  if (default) testthat::local_mocked_bindings(augmented_cell_rows =
    function(..., pair_predictor = NULL, log_derivative = NULL)
      cell(..., pair_predictor = v47_ns("predict_pair_selection"),
        log_derivative = v47_ns("selection_pair_cohort_log_deriv")), .package = "cdid")
  if (clip) testthat::local_mocked_bindings(fit_selection_models = function(...) {
    b <- selector(...)
    for (j in seq_along(b$blocks)) {
      bl <- b$blocks[[j]]
      if (identical(bl$branch, "cp04_separate_pair")) {
        bl$coefficients[] <- 0; bl$coefficients[1L] <- -30
      } else { bl$alpha[] <- -30; bl$beta[] <- 0 }
      b$blocks[[j]] <- bl
    }
    b["ref"] <- list(NULL); v47_ns("ipw_owner_seal")(b, "selection_models")
  }, .package = "cdid")
  run <- function() {
    if (is.null(pair_outcomes)) return(do.call(v47_fixture$rap_fit, c(list(a = a), spec)))
    choice <- utils::modifyList(list(imp = "separate", sel = "separate", model = "pair_logit", mode = "separate"), spec)
    v47_ns("fit_augmented_backend")(a$reader, a$views, a$control, v47_fixture$rap_gc, 1e-10,
      a$recipe, treatment_columns = c("(Intercept)", "treatment_v"), X_frame = a$X_full,
      pair_outcomes = pair_outcomes, imputation_pooling = choice$imp,
      selection_model = choice$model, selection_pooling = choice$sel, outcome_pooling = choice$mode)
  }
  result <- tryCatch(withCallingHandlers(run(), warning = function(w) {
    if (identical(fail, "glm_warning")) stop(paste0("v47 warning handler: ", conditionMessage(w)), call. = FALSE)
  }), error = identity)
  list(result = result, ops = ops, owned = owned, active = active, pairs = pairs, derivs = derivs,
    finishes = finishes, aborts = aborts, expected = expected, predicts = predicts, derivatives = derivatives,
    pure_calls = pure_calls, or_calls = or_calls, probabilities = probabilities,
    clipped_derivatives = clipped_derivatives)
}

testthat::test_that("augmented fitting shares only owned scans with exact nuisance and full-F IF results", {
  specs <- list(list(imp = "separate", sel = "separate", model = "pair_logit", mode = "separate"),
    list(imp = "within_cohort", sel = "across_cohorts", model = "pair_logit", mode = "pooled"),
    list(imp = "across_cohorts", sel = "within_cohort", model = "independent_period_logit", mode = "separate"),
    list(imp = "separate", sel = "across_cohorts", model = "independent_period_logit", mode = "pooled"))
  for (spec in specs) {
    a <- v47_fixture$rap_case(); fast <- v47_fit(a, spec); old <- v47_fit(a, spec, default = TRUE)
    fit <- fast$result; reference <- old$result
    testthat::expect_s3_class(fit, "cdid_augmented_backend")
    testthat::expect_s3_class(reference, "cdid_augmented_backend")
    testthat::expect_identical(fast$owned, c(TRUE, TRUE))
    testthat::expect_false(any(fast$active))
    testthat::expect_identical(fast$finishes, c(1L, 1L))
    testthat::expect_true(all(fast$aborts >= 1L))
    testthat::expect_null(fast$expected[[1L]])
    testthat::expect_identical(fast$expected[[2L]], fast$ops[[1L]]$sources)
    testthat::expect_true(all(fast$derivs > 0L))
    testthat::expect_gt(fast$pairs[1L], 0L)
    testthat::expect_identical(fast$predicts, 0L)
    testthat::expect_identical(fast$derivatives, 0L)
    testthat::expect_gt(old$predicts + old$derivatives, 0L)
    testthat::expect_identical(fast$pure_calls, 1L)
    testthat::expect_gt(fast$or_calls, 0L)
    for (op in fast$ops) testthat::expect_error(op$pair(1L, 1L), class = "cdid_cp04_integrity_error")
    for (field in c("cells", "groups", "treatments", "links", "shared_or", "accounting"))
      testthat::expect_identical(fit[[field]], reference[[field]])
    testthat::expect_identical(fit$selection_backend$coefficients, reference$selection_backend$coefficients)
    testthat::expect_equal(v47_fixture$rap_point(fit), v47_fixture$rap_reference(a, fit), tolerance = 2e-8)
    N <- nrow(a$Y); actual <- scalar <- matrix(0, N, 3L)
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(first, min(N, first + 10L)))
      x <- v47_ns("evaluate_augmented_if_chunk")(fit, rows, 1:3)
      y <- v47_ns("evaluate_augmented_if_chunk")(reference, rows, 1:3)
      testthat::expect_identical(x, y)
      actual[rows, ] <- x$phi; scalar[rows, ] <- y$phi
    }
    testthat::expect_identical(crossprod(actual), crossprod(scalar))
    testthat::expect_identical(sqrt(diag(crossprod(actual))) / N, sqrt(diag(crossprod(scalar))) / N)
    for (i in c(1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      derivative <- tail((v47_fixture$rap_reference(a, fit, 1 + 1e-6 * delta) -
        v47_fixture$rap_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6, 3L)
      testthat::expect_equal(actual[i, ], derivative, tolerance = 1e-5)
    }
  }
})

testthat::test_that("fitting preserves clipped derivative masks, warning errors and callback source boundaries", {
  a <- v47_fixture$rap_case()
  for (model in c("pair_logit", "independent_period_logit")) {
    fast <- v47_fit(a, list(model = model), clip = TRUE)
    old <- v47_fit(a, list(model = model), default = TRUE, clip = TRUE)
    testthat::expect_s3_class(fast$result, "cdid_augmented_backend")
    testthat::expect_true(any(fast$probabilities < 1e-6))
    testthat::expect_gt(fast$clipped_derivatives, 0L)
    testthat::expect_identical(v47_fixture$rap_point(fast$result), v47_fixture$rap_point(old$result))
    testthat::expect_identical(fast$result$links, old$result$links)
  }
  for (kind in c("glm_warning", "or_error")) {
    fast <- v47_fit(a, fail = kind); old <- v47_fit(a, default = TRUE, fail = kind)
    testthat::expect_s3_class(fast$result, "condition")
    testthat::expect_identical(class(fast$result), class(old$result))
    testthat::expect_identical(conditionMessage(fast$result), conditionMessage(old$result))
    testthat::expect_false(any(fast$active))
    if (kind == "glm_warning") testthat::expect_length(fast$ops, 0L) else
      testthat::expect_identical(fast$finishes, 1L)
  }
  changed <- v47_fit(a, fail = "or_source")
  testthat::expect_s3_class(changed$result, "cdid_cp04_integrity_error")
  testthat::expect_false(any(changed$active))
  testthat::expect_identical(changed$finishes, 1L)
})


testthat::test_that("reference and retained private pair fits retain standalone operation validation", {
  a <- v47_fixture$rap_case(dated = TRUE, route_kind = "augmented")
  fit <- v47_fit(a, list(imp = "within_cohort", sel = "across_cohorts"))
  old <- v47_fit(a, list(imp = "within_cohort", sel = "across_cohorts"), default = TRUE)
  testthat::expect_s3_class(fit$result, "cdid_augmented_backend")
  testthat::expect_false(any(fit$owned))
  testthat::expect_gt(fit$predicts + fit$derivatives, 0L)
  testthat::expect_identical(v47_fixture$rap_point(fit$result), v47_fixture$rap_point(old$result))
  a <- v47_fixture$rap_case(); pair_outcomes <- list()
  for (v in a$views) pair_outcomes[[paste0(v$pair$s, ":", v$pair$t)]] <-
    a$Y[, match(v$pair$t, a$calendar)] - a$Y[, match(v$pair$s, a$calendar)]
  fit <- v47_fit(a, pair_outcomes = pair_outcomes)
  old <- v47_fit(a, default = TRUE, pair_outcomes = pair_outcomes)
  testthat::expect_s3_class(fit$result, "cdid_augmented_backend")
  testthat::expect_false(any(fit$owned))
  testthat::expect_gt(fit$predicts + fit$derivatives, 0L)
  testthat::expect_identical(v47_fixture$rap_point(fit$result), v47_fixture$rap_point(old$result))
})
