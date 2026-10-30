# Fixtures and comparisons for ordinary-mode parity with did.

parity_spec <- function() {
  sim_spec(periods = 1:4, cohorts = c(3, 4), p_x = 0.5,
           cohort_prob = rbind(c(0.4, 0.3, 0.3), c(0.2, 0.4, 0.4)),
           theta = c(1, 2), cohort_level = c(0, 0.5, -0.5), n_clusters = 12L)
}

# Balanced panel with a time-invariant sampling-weight column.
parity_panel <- function(n = 240, seed = 101) {
  d <- sim_draw(parity_spec(), n = n, panel = TRUE, seed = seed)
  d$y0 <- NULL
  d$w <- 0.5 + (d$id %% 4) / 2
  d
}

# Unbalanced panel: a fixed set of rows removed.
parity_unbalanced <- function() {
  d <- parity_panel()
  d[-c(2, 7, 15, 40, 41, 99, 230, 517), ]
}

parity_rc <- function(n = 150, seed = 202) {
  d <- sim_draw(parity_spec(), n = n, panel = FALSE, seed = seed)
  d$y0 <- NULL
  d$w <- 0.5 + (d$id %% 3) / 2
  d
}

# Custom estimators following did's documented est_method contracts. `shift`
# arrives only through att_gt()'s `...`, so it makes forwarding observable.
custom_panel <- function(y1, y0, D, covariates, i.weights = NULL,
                         inffunc = TRUE, shift = 0, ...) {
  w <- if (is.null(i.weights)) rep(1, length(D)) else i.weights
  w <- w / mean(w)
  dy <- y1 - y0
  p1 <- mean(w * D)
  p0 <- mean(w * (1 - D))
  m1 <- sum(w * D * dy) / sum(w * D)
  m0 <- sum(w * (1 - D) * dy) / sum(w * (1 - D))
  inf <- w * D * (dy - m1) / p1 - w * (1 - D) * (dy - m0) / p0
  list(ATT = m1 - m0 + shift, att.inf.func = inf)
}

custom_rc <- function(y, post, D, covariates, i.weights = NULL,
                      inffunc = TRUE, shift = 0, ...) {
  w <- if (is.null(i.weights)) rep(1, length(D)) else i.weights
  w <- w / mean(w)
  cell <- function(d, p) {
    k <- w * (D == d) * (post == p)
    m <- sum(k * y) / sum(k)
    list(m = m, inf = k * (y - m) / mean(k))
  }
  c11 <- cell(1, 1); c10 <- cell(1, 0); c01 <- cell(0, 1); c00 <- cell(0, 0)
  list(ATT = (c11$m - c10$m) - (c01$m - c00$m) + shift,
       att.inf.func = (c11$inf - c10$inf) - (c01$inf - c00$inf))
}

# Evaluate, collecting warnings and messages in order (errors propagate).
with_conditions <- function(expr) {
  warnings <- character()
  messages <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, warnings = warnings, messages = messages)
}

# did creates a ~1 formula inside its own preprocessing frame when xformla
# is absent or literally NULL. Its environment retains the whole call frame,
# including delegated call text, but ~1 has no variable lookup. User-supplied
# formulas keep their original environments, including explicitly supplied ~1.
did_default_formula <- function(x) {
  dp <- x$DIDparams
  f <- dp$xformla
  if (!inherits(f, "formula") || length(f) != 2L ||
      !identical(f[[2L]], 1) || is.null(environment(f))) return(FALSE)
  call_args <- as.list(dp$call)
  absent_or_null <- !("xformla" %in% names(call_args)) ||
    is.null(call_args[["xformla"]])
  absent_or_null && identical(parent.env(environment(f)), asNamespace("did"))
}

# Normalize only recorded calls, data.table bookkeeping pointers and the
# proven generated-intercept environment. All statistical fields remain exact.
normalize_did <- function(x, default_formula = did_default_formula(x)) {
  strip <- function(o) {
    if (inherits(o, "data.table")) {
      attr(o, ".internal.selfref") <- NULL
      return(o)
    }
    if (is.list(o) && !is.data.frame(o) && length(o)) {
      o[] <- lapply(o, strip)
    }
    o
  }
  if (default_formula) {
    f <- x$DIDparams$xformla
    stopifnot(inherits(f, "formula"), length(f) == 2L, identical(f[[2L]], 1),
              identical(parent.env(environment(f)), asNamespace("did")))
    environment(f) <- baseenv()
    x$DIDparams$xformla <- f
  }
  if (is.list(x$DIDparams)) x$DIDparams$call <- NULL
  if (!is.null(x$call)) x$call <- NULL
  strip(x)
}

mp_fields <- c("group", "t", "att", "se", "c", "V_analytical", "inffunc",
               "n", "W", "Wpval", "alp")
aggte_fields <- c("overall.att", "overall.se", "type", "egt", "att.egt",
                  "se.egt", "crit.val.egt", "inf.function", "min_e", "max_e",
                  "balance_e")

expect_same_result <- function(ours, ref, fields) {
  expect_identical(class(ours$value), class(ref$value))
  for (f in fields) {
    expect_identical(ours$value[[f]], ref$value[[f]], label = paste0("ours$", f))
  }
  expect_identical(ours$warnings, ref$warnings)
  expect_identical(ours$messages, ref$messages)
  default_formula <- did_default_formula(ref$value)
  expect_equal(normalize_did(ours$value, default_formula),
               normalize_did(ref$value, default_formula), tolerance = 0)
}

# Evaluate the same att_gt(...) expression with cdid (adding `mode`, unless
# NULL) and with did::att_gt, from the caller's environment, after resetting
# the seed each time. Returns both results invisibly.
expect_att_gt_parity <- function(expr, mode = "ordinary", seed = 20261013) {
  call <- substitute(expr)
  env <- parent.frame()
  ours_call <- call
  ours_call[[1L]] <- quote(cdid::att_gt)
  if (!is.null(mode)) ours_call$mode <- mode
  ref_call <- call
  ref_call[[1L]] <- quote(did::att_gt)
  set.seed(seed)
  ours <- with_conditions(eval(ours_call, env))
  ours_rng <- .Random.seed
  set.seed(seed)
  ref <- with_conditions(eval(ref_call, env))
  expect_identical(ours_rng, .Random.seed)
  expect_same_result(ours, ref, mp_fields)
  invisible(list(ours = ours$value, ref = ref$value))
}

expect_aggte_parity <- function(pair, ..., seed = 20261014) {
  set.seed(seed)
  ours <- with_conditions(cdid::aggte(pair$ours, ...))
  set.seed(seed)
  ref <- with_conditions(did::aggte(pair$ref, ...))
  expect_same_result(ours, ref, aggte_fields)
  invisible(list(ours = ours$value, ref = ref$value))
}

# Same error message from cdid (ordinary mode) and did.
expect_same_error <- function(expr, mode = "ordinary") {
  call <- substitute(expr)
  env <- parent.frame()
  ours_call <- call
  ours_call[[1L]] <- quote(cdid::att_gt)
  if (!is.null(mode)) ours_call$mode <- mode
  ref_call <- call
  ref_call[[1L]] <- quote(did::att_gt)
  e_ours <- tryCatch(suppressWarnings(suppressMessages(eval(ours_call, env))),
                     error = function(e) e)
  e_ref <- tryCatch(suppressWarnings(suppressMessages(eval(ref_call, env))),
                    error = function(e) e)
  expect_s3_class(e_ref, "error")
  expect_s3_class(e_ours, "error")
  expect_identical(conditionMessage(e_ours), conditionMessage(e_ref))
}
