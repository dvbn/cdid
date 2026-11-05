# Focused context lifetime/check-count witnesses. Existing independent selection
# frequency, dense-GLM, full-F IF and dated-reference fixtures remain the
# numerical anchors; these fixtures and tolerances are not regenerated here.
spc_ns <- function(n) getFromNamespace(n, "cdid")
spc_case <- function(universe = "F", zero = FALSE) {
  N <- 28L; id <- seq_len(N); g <- rep(c(3, Inf), each = 14L)
  x <- rep(-2:4, 4L); w <- 1 + (id %% 3L) / 10
  Y <- vapply(1:3, function(d) d + d * x / 3 + sin(id + d) / 10 +
    as.double(g == 3) * d, numeric(N))
  observed <- matrix(TRUE, N, 3L)
  for (d in 1:3) observed[id %% 4L == d - 1L, d] <- FALSE
  observed[c(13L, 27L), ] <- FALSE; observed[c(13L, 27L), 1L] <- TRUE
  observed[c(14L, 28L), ] <- FALSE
  ctrl <- spc_ns("ipw_owner_control")(5L, 1e8, 1e8, 1e14, 1e8)
  data <- data.frame(id = rep(id, each = 3L), time = rep(1:3, N),
    y = as.vector(t(Y)), g = rep(g, each = 3L))
  data <- data[as.vector(t(observed)), ]
  panel <- spc_ns("construct_ipw_panel_owner")(data, "id", "time", "y", "g", ctrl)
  X <- cbind(`(Intercept)` = 1, x = x)
  frame <- spc_ns("construct_ipw_frame_owner")(panel, id, X, w, ctrl,
    frame_g = g, allow_frame_only = TRUE)
  U <- switch(universe, F = rep(TRUE, N), B = rowSums(observed) > 0L, A = rowSums(observed) >= 2L)
  Z <- matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept"))
  if (zero) Z <- matrix(numeric(), N, 0L)
  reader <- spc_ns("construct_cp04_frame_reader")(frame, id, Z, g, w, ctrl, U = U,
    universe = universe, justification = if (universe == "F") NULL else "fixed finite-law context witness",
    has_intercept = !zero)
  views <- lapply(1:2, function(j) spc_ns("construct_ipw_link_view")(frame,
    j, 3, 3, Inf, control = ctrl, link_id = paste0("prediction-context-", j),
    functional_key = paste0("cmpcat:prediction-context-", j)))
  list(reader = reader, views = views, frame = frame, panel = panel, control = ctrl,
    observed = observed, w = w, U = U, g = g,
    recipe = spc_ns("or_pooled_outcome_recipe")(colnames(X), TRUE, "x"))
}
spc_fit <- function(a, model = "pair_logit", pooling = "separate")
  spc_ns("fit_selection_models")(a$reader, a$views, selection_model = model,
    selection_pooling = pooling, control = a$control,
    glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE))
spc_rebind <- function(env, binding, value) {
  old <- get(binding, env, inherits = FALSE)
  unlockBinding(binding, env); assign(binding, value, env); lockBinding(binding, env)
  function() { unlockBinding(binding, env); assign(binding, old, env); lockBinding(binding, env) }
}

testthat::test_that("bounded predictions match independent empirical frequencies across model/pooling/universe branches", {
  for (universe in c("F", "B", "A")) for (model in c("pair_logit", "independent_period_logit")) {
    if (model == "independent_period_logit" && universe != "F") next
    a <- spc_case(universe)
    for (pooling in c("separate", "within_cohort", "across_cohorts")) {
      fit <- spc_fit(a, model, pooling); op <- spc_ns("selection_prediction_operation")(fit)
      up_map <- spc_ns("build_selection_upstream_map")(fit)
      for (q in seq_along(fit$pair_cohorts)) {
        pc <- fit$pair_cohorts[[q]]; k <- a$reader$state$cohorts[pc$cohort]
        U <- a$U & a$g == k
        rows <- rev(which(U))
        s <- pc$pair$s_calendar_index; t <- pc$pair$t_calendar_index
        mean_response <- function(v) sum(a$w[U] * v[U]) / sum(a$w[U])
        ps <- mean_response(a$observed[, s]); pt <- mean_response(a$observed[, t])
        expected <- if (model == "pair_logit") mean_response(a$observed[, s] * a$observed[, t]) else ps * pt
        pred <- op$pair(rows = rows, pair_cohort_id = q)
        testthat::expect_equal(pred$prob, rep(expected, length(rows)), tolerance = 1e-12)
        testthat::expect_equal(pred$weight, a$w[rows] * a$reader$state$a[pc$cohort] / expected, tolerance = 1e-12)
        testthat::expect_identical(pred, spc_ns("predict_pair_selection")(fit, rows, q))
        if (model == "independent_period_logit") {
          testthat::expect_equal(pred$prob_s, rep(ps, length(rows)), tolerance = 1e-12)
          testthat::expect_equal(pred$prob_t, rep(pt, length(rows)), tolerance = 1e-12)
          testthat::expect_identical(pred, op$period(rows, q))
        }
        testthat::expect_identical(op$deriv(up_map, q, rows, pred = pred),
          spc_ns("selection_pair_cohort_log_deriv")(fit, up_map, q, rows, pred = pred))
      }
      # All F rows include singleton, zero-observer and frame-only units.
      rows <- 1:5
      testthat::expect_identical(op$upstream(up_map, rows, seq_len(up_map$D), TRUE),
        spc_ns("evaluate_selection_upstream_chunk")(fit, up_map, rows, seq_len(up_map$D), TRUE))
      op$finish()
      testthat::expect_error(op$cell(1L, 1L), class = "cdid_cp04_integrity_error")
    }
  }
})

testthat::test_that("zero-column and standalone argument/refusal behavior is retained", {
  a <- spc_case(zero = TRUE)
  for (model in c("pair_logit", "independent_period_logit")) {
    fit <- spc_fit(a, model, "across_cohorts"); op <- spc_ns("selection_prediction_operation")(fit)
    testthat::expect_identical(op$cell(c(5, 1, 5), 1L), rep(0.5, 3L))
    testthat::expect_identical(op$pair(5:1, 1L)$prob, rep(if (model == "pair_logit") 0.5 else 0.25, 5L))
    for (rows in list(integer(), c(1, NA), 0L, 29L, matrix(1L), 1.5)) {
      testthat::expect_error(op$cell(rows, 1L), class = "cdid_cp04_input_error")
      testthat::expect_error(spc_ns("predict_selection_cell")(fit, rows, 1L), class = "cdid_cp04_input_error")
    }
    testthat::expect_error(op$pair(1L, 0L), class = "cdid_cp04_input_error")
    if (model == "pair_logit") testthat::expect_error(op$period(1L, 1L), class = "cdid_cp04_unsupported_error")
    op$abort()
    testthat::expect_error(op$pair(1L, 1L), class = "cdid_cp04_integrity_error")
  }
})

testthat::test_that("entry, exit and restart reject mutation and foreign exact sources", {
  a <- spc_case(); fit <- spc_fit(a, "independent_period_logit", "within_cohort")
  bad <- fit; bad$blocks[[1L]]$alpha[1L] <- bad$blocks[[1L]]$alpha[1L] + 0.1
  testthat::expect_error(spc_ns("selection_prediction_operation")(bad), class = "cdid_cp04_integrity_error")
  other <- spc_case()
  testthat::expect_error(spc_ns("selection_prediction_operation")(fit, source_reader = other$reader),
    class = "cdid_cp04_integrity_error")
  for (sw in list(list(env = fit$ref, binding = "selection_models", value = bad),
      list(env = a$reader$ref, binding = "cp04_reader", value = other$reader),
      list(env = a$frame$ref, binding = "frame_owner", value = other$frame),
      list(env = a$panel$ref, binding = "panel_owner", value = other$panel))) {
    op <- spc_ns("selection_prediction_operation")(fit)
    op$pair(1:5, 1L)
    restore <- spc_rebind(sw$env, sw$binding, sw$value)
    error <- tryCatch(op$finish(), error = identity)
    restore()
    testthat::expect_s3_class(error, "cdid_cp04_integrity_error")
    testthat::expect_error(op$pair(1L, 1L), class = "cdid_cp04_integrity_error")
  }
  op <- spc_ns("selection_prediction_operation")(fit); saved <- op$sources; op$finish(); op <- NULL
  fresh <- spc_ns("selection_prediction_operation")(fit, expected_sources = saved)
  fresh$finish()
  restore <- spc_rebind(a$frame$ref, "frame_owner", other$frame)
  error <- tryCatch(spc_ns("selection_prediction_operation")(fit, expected_sources = saved), error = identity)
  restore()
  testthat::expect_s3_class(error, "cdid_error")
})

testthat::test_that("one full source validation serves a bounded operation and standalone upstream chunk", {
  a <- spc_case(); fit <- spc_fit(a, "independent_period_logit", "across_cohorts")
  up_map <- spc_ns("build_selection_upstream_map")(fit)
  original <- spc_ns("validate_selection_models"); calls <- 0L
  testthat::local_mocked_bindings(validate_selection_models = function(backend) {
    calls <<- calls + 1L; original(backend)
  }, .package = "cdid")
  op <- spc_ns("selection_prediction_operation")(fit)
  for (j in 1:4) {
    op$pair(5:1, 1L); op$cell(1:5, 1L)
    op$upstream(up_map, 1:5, seq_len(up_map$D), TRUE)
  }
  op$finish()
  testthat::expect_identical(calls, 1L)
  calls <- 0L
  spc_ns("evaluate_selection_upstream_chunk")(fit, up_map, 1:5, seq_len(up_map$D), TRUE)
  testthat::expect_identical(calls, 1L)
})

testthat::test_that("successful exit keeps live DLL verification and closes on failure", {
  a <- spc_case(); fit <- spc_fit(a)
  op <- spc_ns("selection_prediction_operation")(fit)
  op$pair(1:5, 1L)
  testthat::local_mocked_bindings(gmm_require_native_dll = function()
    stop("live selection-context DLL verification failure"), .package = "cdid")
  testthat::expect_error(op$finish(), "live selection-context DLL verification failure", fixed = TRUE)
  testthat::expect_error(op$cell(1L, 1L), class = "cdid_cp04_integrity_error")
})

testthat::test_that("stored observer predictions give the same treatment derivatives as recomputation", {
  a <- spc_case()
  gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
  for (model in c("pair_logit", "independent_period_logit")) {
    sel <- spc_fit(a, model, "within_cohort")
    fit <- function(mode) spc_ns("fit_dr_corollary1")(a$reader, a$views, a$control, gc, 1e-10,
      a$recipe, treatment_columns = "(Intercept)", selection_model = model,
      selection_pooling = "within_cohort", sel_fit = sel, outcome_pooling = mode)
    for (mode in c("pooled", "separate")) {
      reused <- fit(mode)
      recompute <- function() {
        original <- spc_ns("selection_prediction_operation")
        testthat::local_mocked_bindings(selection_prediction_operation = function(...) {
          op <- original(...); derivative <- op$deriv
          # Deliberately recompute the same raw pair/endpoints in each
          # derivative chunk, rather than accepting collector pp/w buffers.
          op$deriv <- function(up_map, q, rows, pred = NULL, Z = NULL, a_vec = NULL)
            derivative(up_map, q, rows, pred = NULL, Z = Z, a_vec = a_vec)
          op
        }, .package = "cdid")
        fit(mode)
      }
      repeated <- recompute()
      testthat::expect_identical(reused$treatments, repeated$treatments)
      testthat::expect_identical(reused$links, repeated$links)
      fn <- spc_ns("evaluate_dr_corollary1_if_chunk")
      for (details in c(FALSE, TRUE)) for (rows in list(1:5, 25:28))
        testthat::expect_identical(fn(reused, rows, 2:1, details = details),
          fn(repeated, rows, 2:1, details = details))
    }
  }
})

testthat::test_that("treatment GLM executes after closing prediction permission and restart catches handler mutation", {
  a <- spc_case(); sel <- spc_fit(a, "independent_period_logit", "within_cohort")
  original_op <- spc_ns("selection_prediction_operation"); original_glm <- spc_ns("cp04_logit")
  last <- NULL; observed_closed <- 0L
  testthat::local_mocked_bindings(selection_prediction_operation = function(...) {
    last <<- original_op(...); last
  }, cp04_logit = function(...) {
    if (!is.null(last)) {
      testthat::expect_error(last$pair(1L, 1L), class = "cdid_cp04_integrity_error")
      observed_closed <<- observed_closed + 1L
    }
    original_glm(...)
  }, .package = "cdid")
  gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
  spc_ns("fit_dr_corollary1")(a$reader, a$views, a$control, gc, 1e-10,
    a$recipe, treatment_columns = "(Intercept)", selection_model = "independent_period_logit",
    selection_pooling = "within_cohort", sel_fit = sel, outcome_pooling = "separate")
  testthat::expect_gt(observed_closed, 0L)
  # Equivalent warning/user-handler boundary witness: a closed context grants
  # no permission to restart on a rebound owner, even when that new owner is
  # independently valid. The actual GLM numerical result above is unchanged.
  saved <- last$sources; last$abort()
  other <- spc_case(); restore <- spc_rebind(a$frame$ref, "frame_owner", other$frame)
  error <- tryCatch(original_op(sel, expected_sources = saved), error = identity)
  restore()
  testthat::expect_s3_class(error, "cdid_error")
})

testthat::test_that("source promises resolve before entry and active snapshots never grant permission", {
  a <- spc_case(); other <- spc_case(); fit <- spc_fit(a)
  source_reader <- function() {
    restore <- spc_rebind(a$frame$ref, "frame_owner", other$frame)
    list(value = a$reader, restore = restore)
  }
  changed <- NULL
  error <- tryCatch(spc_ns("selection_prediction_operation")(fit,
    source_reader = { changed <- source_reader(); changed$value }), error = identity)
  changed$restore()
  testthat::expect_s3_class(error, "cdid_error")
  op <- spc_ns("selection_prediction_operation")(fit); src <- op$sources; op$finish()
  active <- new.env(parent = emptyenv()); called <- FALSE
  for (nm in ls(src)) {
    if (nm == "reader") makeActiveBinding(nm, function(value) { called <<- TRUE; a$reader }, active) else
      assign(nm, get(nm, src, inherits = FALSE), active)
    lockBinding(nm, active)
  }
  lockEnvironment(active, bindings = TRUE)
  testthat::expect_error(spc_ns("selection_prediction_operation")(fit, expected_sources = active),
    class = "cdid_cp04_integrity_error")
  testthat::expect_false(called)
})
