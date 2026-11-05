# Existing deterministic raw nuisance/refit reference, loaded without running
# another test file's expressions or changing fixture inputs/outputs.
v45_ns <- function(name) getFromNamespace(name, "cdid")
v45_fixture <- local({
  e <- new.env(parent = environment())
  wanted <- c("rap_ns", "rap_gc", "rap_case", "rap_reference", "rap_point", "rap_fit")
  for (x in parse(file = testthat::test_path("test-augmented-outcome-pooling.R"))) {
    if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) &&
        as.character(x[[2L]]) %in% wanted) eval(x, e)
  }
  stopifnot(all(vapply(wanted, exists, logical(1L), envir = e, inherits = FALSE)))
  e
})
v45_run <- function(fit, rows, links = 1:3, default = FALSE, details = FALSE, capture_error = FALSE) {
  validate <- v45_ns("validate_selection_models")
  predict <- v45_ns("predict_pair_selection")
  factory <- v45_ns("selection_prediction_operation")
  core <- v45_ns("augmented_evaluate_if_chunk_core")
  validations <- predictions <- 0L; ops <- list(); pairs <- finishes <- aborts <- integer(); probabilities <- numeric(); owned <- logical()
  testthat::local_mocked_bindings(
    validate_selection_models = function(...) { validations <<- validations + 1L; validate(...) },
    predict_pair_selection = function(...) { predictions <<- predictions + 1L; predict(...) },
    selection_prediction_operation = function(...) {
      op <- factory(...); j <- length(ops) + 1L
      pairs[j] <<- 0L; finishes[j] <<- 0L; aborts[j] <<- 0L
      owned[j] <<- "source_reader" %in% names(list(...))
      pair <- op$pair; finish <- op$finish; abort <- op$abort
      op$pair <- function(...) {
        value <- pair(...); pairs[j] <<- pairs[j] + 1L
        probabilities <<- c(probabilities, value$prob); value
      }
      op$finish <- function() { finishes[j] <<- finishes[j] + 1L; finish() }
      op$abort <- function() { aborts[j] <<- aborts[j] + 1L; abort() }
      ops[[j]] <<- op; op
    }, .package = "cdid")
  if (default) testthat::local_mocked_bindings(augmented_evaluate_if_chunk_core =
    function(..., cell_worker = NULL) core(..., cell_worker = v45_ns("augmented_cell_rows")), .package = "cdid")
  run <- function() v45_ns("evaluate_augmented_if_chunk")(fit, rows, links, details = details)
  out <- if (capture_error) tryCatch(run(), error = identity) else run()
  list(out = out, validations = validations, predictions = predictions, ops = ops,
    pairs = pairs, finishes = finishes, aborts = aborts, probabilities = probabilities, owned = owned)
}
v45_rebind <- function(env, name, value) {
  old <- get(name, env, inherits = FALSE)
  unlockBinding(name, env); assign(name, value, env); lockBinding(name, env)
  function() { unlockBinding(name, env); assign(name, old, env); lockBinding(name, env) }
}
v45_reseal <- function(fit, selection) {
  b <- fit; b["ref"] <- list(NULL); b$selection_backend <- selection
  v45_ns("ipw_owner_seal")(b, "augmented_backend")
}

testthat::test_that("owned augmented replay reuses prediction validation with identical complete-F science", {
  specs <- list(list(imp = "separate", sel = "separate", model = "pair_logit", mode = "separate"),
    list(imp = "within_cohort", sel = "across_cohorts", model = "pair_logit", mode = "pooled"),
    list(imp = "across_cohorts", sel = "within_cohort", model = "independent_period_logit", mode = "separate"),
    list(imp = "separate", sel = "across_cohorts", model = "independent_period_logit", mode = "pooled"))
  for (spec in specs) {
    a <- v45_fixture$rap_case(); fit <- do.call(v45_fixture$rap_fit, c(list(a = a), spec))
    N <- nrow(a$Y); actual <- reference <- matrix(0, N, 3L)
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(first, min(N, first + 10L)))
      fast <- v45_run(fit, rows); old <- v45_run(fit, rows, default = TRUE)
      testthat::expect_identical(fast$out, old$out)
      testthat::expect_identical(fast$predictions, 0L)
      testthat::expect_identical(length(fast$ops), 2L) # wrapper + unchanged upstream
      testthat::expect_identical(fast$owned, c(TRUE, FALSE))
      testthat::expect_identical(old$validations - fast$validations, old$predictions)
      testthat::expect_equal(sum(fast$pairs), old$predictions)
      testthat::expect_identical(fast$finishes, c(1L, 1L))
      testthat::expect_true(all(fast$aborts >= 1L))
      for (op in fast$ops) testthat::expect_error(op$cell(1L, 1L), class = "cdid_cp04_integrity_error")
      actual[rows, ] <- fast$out$phi; reference[rows, ] <- old$out$phi
    }
    testthat::expect_identical(actual, reference)
    testthat::expect_identical(crossprod(actual), crossprod(reference))
    testthat::expect_identical(sqrt(diag(crossprod(actual))) / N, sqrt(diag(crossprod(reference))) / N)
    testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
    testthat::expect_equal(v45_fixture$rap_point(fit), v45_fixture$rap_reference(a, fit), tolerance = 2e-8)
    for (i in c(1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      derivative <- tail((v45_fixture$rap_reference(a, fit, 1 + 1e-6 * delta) -
        v45_fixture$rap_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6, 3L)
      testthat::expect_equal(actual[i, ], derivative, tolerance = 1e-5)
    }
    rich <- v45_run(fit, 11:1, c(3L, 1L), details = TRUE)
    testthat::expect_gt(rich$predictions, 0L)
    testthat::expect_false(any(rich$owned))
    testthat::expect_identical(colnames(rich$out$phi), vapply(fit$links[c(3L, 1L)], `[[`, character(1L), "link_id"))
    testthat::expect_equal(unname(rich$out$phi), actual[11:1, c(3L, 1L), drop = FALSE], tolerance = 1e-10)
  }
})

testthat::test_that("standalone cell/default private core, clipping and output/source boundaries remain checked", {
  a <- v45_fixture$rap_case(); fit <- v45_fixture$rap_fit(a)
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  X <- v45_ns("augmented_backend_X_F")(fit); cf <- fit$cells[[1L]]; gf <- fit$groups[[cf$group_id]]
  rows <- head(which(reader$state$cohort_index == cf$cohort_index), 11L)
  calls <- 0L
  pair_reader <- function(rows, pair, read_outcomes = TRUE) {
    calls <<- calls + 1L; v45_ns("cp04_pair_rows")(reader, rows, pair, read_outcomes)
  }
  cell <- v45_ns("augmented_cell_rows")
  deriv <- cell(fit$selection_backend, fit$up_map, X, pair_reader, cf, gf, rows, TRUE)
  testthat::expect_identical(deriv, cell(fit$selection_backend, fit$up_map, X, pair_reader, cf, gf, rows, TRUE,
    pair_predictor = function(...) stop("derivative branch must not use predictor")))
  testthat::expect_gt(calls, 0L)
  pre <- v45_ns("augmented_if_budget")(fit, length(rows), 1:3, FALSE)
  core <- v45_ns("augmented_evaluate_if_chunk_core")(reader$state, X, pair_reader, fit, rows, 1:3, sets = pre)
  testthat::expect_equal(core$phi, v45_run(fit, rows)$out$phi, tolerance = 1e-10)
  # Frozen-fit primitive branch witness: the same existing nuisance records
  # are resealed at a finite low intercept; no DGP/redraw/refit/control changes.
  for (model in c("pair_logit", "independent_period_logit")) {
    b <- if (model == "pair_logit") fit else v45_fixture$rap_fit(a, model = model)
    selection <- b$selection_backend
    for (j in seq_along(selection$blocks)) {
      bl <- selection$blocks[[j]]
      if (identical(bl$branch, "cp04_separate_pair")) {
        bl$coefficients[] <- 0; bl$coefficients[1L] <- -30
      } else { bl$alpha[] <- -30; bl$beta[] <- 0 }
      selection$blocks[[j]] <- bl
    }
    selection["ref"] <- list(NULL); selection <- v45_ns("ipw_owner_seal")(selection, "selection_models")
    b <- v45_reseal(b, selection)
    fast <- v45_run(b, rows); old <- v45_run(b, rows, default = TRUE)
    testthat::expect_true(any(fast$probabilities < 1e-6))
    testthat::expect_identical(v45_ns("ipw_clip_floor")(min(fast$probabilities)), 1e-6)
    testthat::expect_identical(fast$out, old$out)
  }
  during <- function() {
    restore <- NULL; original <- v45_ns("augmented_cell_rows")
    bad <- fit$selection_backend; bad$coefficients[1L] <- Inf
    on.exit(if (!is.null(restore)) restore(), add = TRUE)
    testthat::local_mocked_bindings(augmented_cell_rows = function(...) {
      value <- original(...)
      if (is.null(restore)) restore <<- v45_rebind(fit$selection_backend$ref, "selection_models", bad)
      value
    }, .package = "cdid")
    v45_run(fit, rows, capture_error = TRUE)
  }
  changed <- during()
  testthat::expect_s3_class(changed$out, "cdid_cp04_integrity_error")
  testthat::expect_identical(changed$finishes[1L], 1L)
  testthat::expect_gte(changed$aborts[1L], 1L)
  fail_output <- function() {
    payload <- v45_ns("cp04_payload")
    testthat::local_mocked_bindings(cp04_payload = function(x, ...) {
      value <- payload(x, ...)
      if (is.list(x) && identical(names(x), c("rows", "id", "links", "phi", "scale", "normalization_count", "diagnostics")))
        value["entries"] <- pre$output + 1
      value
    }, .package = "cdid")
    v45_run(fit, rows, capture_error = TRUE)
  }
  failed <- fail_output()
  testthat::expect_s3_class(failed$out, "cdid_augmented_bound_error")
  testthat::expect_identical(failed$finishes[1L], 0L)
  testthat::expect_gte(failed$aborts[1L], 1L)
})


testthat::test_that("excluded layout and retained-pair wrappers retain the default cell validation path", {
  a <- v45_fixture$rap_case()
  pair_outcomes <- list()
  for (v in a$views) {
    s <- match(v$pair$s, a$calendar); t <- match(v$pair$t, a$calendar)
    pair_outcomes[[paste0(v$pair$s, ":", v$pair$t)]] <- a$Y[, t] - a$Y[, s]
  }
  retained <- v45_ns("fit_augmented_backend")(a$reader, a$views, a$control, v45_fixture$rap_gc,
    1e-10, a$recipe, treatment_columns = c("(Intercept)", "treatment_v"),
    X_frame = a$X_full, pair_outcomes = pair_outcomes, outcome_pooling = "separate")
  value <- v45_run(retained, 11:1)
  testthat::expect_gt(value$predictions, 0L)
  testthat::expect_false(any(value$owned))
  testthat::expect_identical(value$out, v45_run(retained, 11:1, default = TRUE)$out)
  dated <- v45_fixture$rap_case(dated = TRUE, route_kind = "augmented")
  layout <- v45_fixture$rap_fit(dated, imp = "within_cohort", sel = "across_cohorts")
  value <- v45_run(layout, 11:1)
  testthat::expect_gt(value$predictions, 0L)
  testthat::expect_false(any(value$owned))
  testthat::expect_identical(value$out, v45_run(layout, 11:1, default = TRUE)$out)
})

testthat::test_that("public augmented ATT full-F IF covariance SE and handler status match default cell replay", {
  e <- new.env(parent = as.environment("package:stats"))
  sys.source(testthat::test_path("reference", "imputation_public_reference.R"), e, keep.source = FALSE)
  e$imputation_public_helpers(helper_dir = testthat::test_path("reference"))
  raw <- e$imputation_public_raw(utils::read.csv(testthat::test_path("fixtures", "aug_public_v1", "raw.csv")))
  args <- e$imputation_public_fixture(raw)$args
  args$selection_est_method <- "dr"
  for (spec in list(list(selection_model = "pair_logit", outcome_pooling = "separate"),
      list(selection_model = "independent_period_logit", outcome_pooling = "pooled"))) {
    request <- utils::modifyList(args, spec)
    fast <- do.call(cdid::att_gt, request)
    default_fit <- function() {
      core <- v45_ns("augmented_evaluate_if_chunk_core")
      testthat::local_mocked_bindings(augmented_evaluate_if_chunk_core =
        function(..., cell_worker = NULL) core(..., cell_worker = v45_ns("augmented_cell_rows")), .package = "cdid")
      do.call(cdid::att_gt, request)
    }
    old <- default_fit()
    for (nm in c("att", "inffunc", "V_analytical", "se")) testthat::expect_identical(fast[[nm]], old[[nm]])
    for (nm in c("available", "target_meta", "support_status", "support_reason", "covariance_handling"))
      testthat::expect_identical(fast$chained[[nm]], old$chained[[nm]])
    testthat::expect_true(is.list(fast$chained$covariance_handling))
  }
})
