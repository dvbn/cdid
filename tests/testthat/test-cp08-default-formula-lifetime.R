# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 default-formula lifetime: when the caller supplies no `xformla`, the
# request layer creates the intercept-only formula ~1. A literal ~1 written
# inside chained_request_prepare_inputs() carried that preparation frame (data,
# data_internal, cp04, ...) as its environment and travelled into DIDparams/MP
# through `prep$xformla` (saved no-X r1 results serialized to ~85-90 MB versus
# ~2-3 MB with an explicit ~x). The package default is now built by the private
# constructor chained_request_default_formula() in the base environment.
#
# These tests check (i) the constructor's exact ~1 AST, terms/intercept and
# design metadata, (ii) that no-X preparations and genuine no-X public results
# reach no preparation frame, (iii) that every caller-supplied formula keeps its
# exact identity and environment, and (iv) that static, transform, factor,
# context-dependent (refused), dated-covariate, dated-weight and supplied-frame
# preparation routes keep their existing outputs and errors. Fixtures follow the
# pc_panel/pc_prepare patterns of test-cp08-prepared-context.R.
#
# UNEXECUTED at authoring time (no R in the authoring job). Expected values come
# from base R (model.matrix/terms of a literal ~1), from the explicit-formula
# route of the same installed build, and from the source of the request layer.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

# Benchmark-shaped complete balanced panel in public coding (G = 0 never
# treated): one treated cohort 3, unit covariates x (integer-valued double) and
# z, a factor f, a unit weight w and a unit column of ones.
dfl_panel <- function(N = 40L, T = 6L, seed = 17L, ids = seq_len(N), rows_perm = FALSE, second = FALSE,
                      singleton = FALSE, integer_ids = FALSE) {
  set.seed(seed)
  g <- rep(0, N); g[seq(1L, N, by = 2L)] <- 3
  if (second) g[seq(2L, N, by = 4L)] <- 4
  d <- data.frame(id = rep(ids, each = T), period = rep(seq_len(T), N))
  u <- rep(seq_len(N), each = T)
  d$G <- g[u]
  x <- stats::rbinom(N, 3, 0.5); z <- round(stats::rnorm(N), 3)
  d$x <- as.double(x[u]); d$z <- z[u]
  d$f <- factor(c("a", "b")[1 + (u %% 2L)])
  d$w <- 1 + (u %% 3L) / 4
  d$one <- 1
  d$y <- 0.3 * u + 0.1 * d$period + 0.2 * d$x * d$period + ifelse(d$G > 0 & d$period >= d$G, 1.5, 0) +
    0.35 * sin(0.7 * u * d$period) + stats::rnorm(nrow(d))
  if (singleton) d$y[d$id == ids[N] & d$period > 1L] <- NA_real_
  if (integer_ids) d$id <- as.integer(d$id)
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

dfl_resolved <- function(chunk = 7L) api$chained_public_resolve_control(list(chunk_rows = chunk))

dfl_prepare <- function(d, est = "ipw", xformla = NULL, frame = NULL, weightsname = NULL, fix_weights = NULL,
                        resolved = dfl_resolved()) {
  api$chained_request_prepare_inputs(data = d, frame = frame, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xformla, weightsname = weightsname, anticipation = 0,
    outcome_regression = identical(est, "reg"), doubly_robust = identical(est, "dr"), or_resource = resolved$profile,
    fix_weights = fix_weights)
}

# Bounded reachability walk over list elements, attributes, call elements,
# closure environments and formula environments. Environments are collected,
# never descended (the shape of the saved-evidence walk). Lists are values in R,
# so the walk terminates; the node bound guards against pathological inputs.
dfl_environments <- function(x, max_nodes = 500000L) {
  envs <- list(); n <- 0L
  add_env <- function(e) {
    for (k in seq_along(envs)) if (identical(envs[[k]], e)) return(invisible(NULL))
    envs[[length(envs) + 1L]] <<- e
    invisible(NULL)
  }
  walk <- function(v) {
    n <<- n + 1L
    if (n > max_nodes) stop("reachability walk exceeded its node bound", call. = FALSE)
    if (is.environment(v)) return(add_env(v))
    if (is.function(v)) {
      if (!is.primitive(v)) add_env(environment(v))
      return(invisible(NULL))
    }
    a <- attributes(v)
    if (!is.null(a)) for (k in seq_along(a)) walk(a[[k]])
    if (is.list(v)) {
      for (k in seq_along(v)) walk(v[[k]])
    } else if (is.call(v)) {
      for (k in seq_along(v)) {
        if (identical(v[[k]], quote(expr = ))) next
        walk(v[[k]])
      }
    }
    invisible(NULL)
  }
  walk(x)
  envs
}

# A chained_request_prepare_inputs() evaluation frame is recognisable by its
# local closure `covariate_gate` together with its local `data_internal`
# binding (R02's release sets the invocation maps to NULL but keeps the names).
dfl_is_preparation_frame <- function(e) {
  !identical(e, baseenv()) && !identical(e, emptyenv()) && !identical(e, globalenv()) && !isNamespace(e) &&
    exists("covariate_gate", envir = e, inherits = FALSE) && exists("data_internal", envir = e, inherits = FALSE)
}

dfl_expect_no_preparation_frame <- function(x, label) {
  envs <- dfl_environments(x)
  hits <- vapply(envs, dfl_is_preparation_frame, logical(1L))
  expect_false(any(hits), info = paste0(label, ": ", sum(hits), " chained_request_prepare_inputs() frame(s) reachable"))
  invisible(envs)
}

dfl_strip_env <- function(t) { attr(t, ".Environment") <- NULL; t }

# ---------------------------------------------------------------------------
# 1. The private constructor: exact ~1, base environment, unchanged metadata
# ---------------------------------------------------------------------------

test_that("the private default-formula constructor builds the exact ~1 intercept formula in the base environment", {
  ctor <- api$chained_request_default_formula
  expect_true(is.function(ctor))
  f <- ctor()
  expect_s3_class(f, "formula")
  expect_identical(environment(f), baseenv())
  expect_identical(length(f), 2L)
  expect_identical(all.vars(f), character(0))
  expect_identical(deparse(f), "~1")
  # Exact AST and attributes of a literal ~1 once that literal is moved to the base environment.
  lit <- ~1
  environment(lit) <- baseenv()
  expect_identical(f, lit)
  expect_true(is.double(f[[2L]])); expect_identical(f[[2L]], 1)
  expect_identical(as.character(f[[1L]]), "~")
  # terms/intercept/response metadata equal the literal's apart from the environment.
  tf <- stats::terms(f); tl <- stats::terms(~1)
  expect_identical(attr(tf, "intercept"), 1L)
  expect_identical(attr(tf, "response"), 0L)
  expect_identical(attr(tf, "term.labels"), character(0))
  expect_identical(attr(tf, ".Environment"), baseenv())
  expect_identical(dfl_strip_env(tf), dfl_strip_env(tl))
  # The intercept-only design on a frame is identical to the literal's design.
  d <- dfl_panel(N = 9L, T = 3L, seed = 3L); first <- d[!duplicated(d$id), ]
  expect_identical(stats::model.matrix(f, data = first), stats::model.matrix(~1, data = first))
  expect_identical(stats::model.matrix(f, data = d), stats::model.matrix(~1, data = d))
  # One fresh, equal formula per call; the constructor keeps no state.
  g <- ctor()
  expect_identical(f, g)
  expect_identical(environment(g), baseenv())
  expect_identical(names(formals(ctor)), NULL)
})

# ---------------------------------------------------------------------------
# 2. No-X preparation on the static IPW/OR/DR routes
# ---------------------------------------------------------------------------

test_that("no-X preparations on IPW/OR/DR static routes carry the base-environment default, reach no preparation frame and equal an explicit ~1 field for field", {
  d <- dfl_panel(seed = 17L)
  fenv <- new.env(parent = baseenv())
  one <- eval(quote(~1), fenv)                       # a caller's own ~1 with a caller-owned environment
  for (est in c("ipw", "reg", "dr")) {
    a <- dfl_prepare(d, est = est)                   # xformla = NULL: the package default
    b <- dfl_prepare(d, est = est, xformla = one)    # the caller's explicit ~1
    expect_s3_class(a$xformla, "formula")
    expect_identical(environment(a$xformla), baseenv())
    expect_identical(deparse(a$xformla), "~1")
    expect_identical(all.vars(a$xformla), character(0))
    dfl_expect_no_preparation_frame(a, paste("prepared input", est, "no X"))
    # The caller's ~1 is never replaced by the package default.
    expect_identical(b$xformla, one)
    expect_identical(environment(b$xformla), fenv)
    expect_false(identical(environment(b$xformla), baseenv()))
    # Apart from the formula environment and the sealed record (sealed per call),
    # the two preparations are identical: design, cp04, units, data, basis, routes.
    aa <- a[setdiff(names(a), "prepared")]; bb <- b[setdiff(names(b), "prepared")]
    environment(bb$xformla) <- baseenv()
    expect_identical(aa, bb)
    expect_identical(a$frame_Z, b$frame_Z)
    expect_identical(a$cp04, b$cp04)
    expect_identical(a$or_basis, b$or_basis)
    if (identical(est, "ipw"))
      expect_identical(a$frame_Z, matrix(1, length(a$frame_ids), 1L, dimnames = list(NULL, "(Intercept)")))
    # The prepared-intake eligibility of the intercept-only design is unchanged on both inputs.
    expect_s3_class(a$prepared, "cdid_request_prepared_intake")
    expect_s3_class(b$prepared, "cdid_request_prepared_intake")
  }
})

# ---------------------------------------------------------------------------
# 3. Caller-supplied formulas keep identity and environment on every route
# ---------------------------------------------------------------------------

test_that("caller-supplied formulas keep their exact identity and environment on static, transform, factor, dated and refused routes", {
  d <- dfl_panel(seed = 13L); first <- d[!duplicated(d$id), ]
  fenv <- new.env(parent = globalenv())
  # Static direct columns and a factor run on every route; row-wise transforms
  # decline the direct route and run the NT route on IPW (existing patterns).
  shapes <- list(list(ex = quote(~x), est = c("ipw", "reg", "dr")),
                 list(ex = quote(~x + z), est = c("ipw", "reg", "dr")),
                 list(ex = quote(~f), est = c("ipw", "reg", "dr")),
                 list(ex = quote(~I(x^2)), est = "ipw"),
                 list(ex = quote(~log(1 + x)), est = "ipw"))
  for (sh in shapes) {
    xf <- eval(sh$ex, fenv)
    expect_identical(environment(xf), fenv)
    for (est in sh$est) {
      p <- dfl_prepare(d, est = est, xformla = xf)
      expect_identical(p$xformla, xf)
      expect_identical(environment(p$xformla), fenv)
      expect_true(api$gmm_sexp_aliases(p$xformla, xf))
      dfl_expect_no_preparation_frame(p, paste("prepared input", est, deparse(sh$ex)))
    }
    # Independent reference of the IPW frame design: the first-observed frame's model matrix.
    p <- dfl_prepare(d, est = "ipw", xformla = xf)
    m <- stats::model.matrix(xf, data = first); rownames(m) <- NULL; attr(m, "assign") <- NULL
    expect_identical(p$frame_Z, m)
  }
  # Context-dependent calls keep the existing typed refusal on every route; the
  # caller's object is untouched.
  for (ex in list(quote(~scale(z)), quote(~poly(z, 2)))) {
    xf <- eval(ex, fenv)
    for (est in c("ipw", "reg", "dr")) {
      e <- tryCatch(dfl_prepare(d, est = est, xformla = xf), error = identity)
      expect_s3_class(e, "cdid_chained_request_covariate_error")
    }
    expect_identical(environment(xf), fenv)
    expect_identical(xf, eval(ex, fenv))
  }
  # A genuinely changing covariate still becomes the dated reference route with the caller's formula object.
  dv <- d; dv$x[dv$id == d$id[1L] & dv$period == 2L] <- 9
  xf <- eval(quote(~x), fenv)
  for (est in c("ipw", "reg")) {
    p <- dfl_prepare(dv, est = est, xformla = xf)
    expect_identical(p$weight_reduction, "reference_recipes")
    expect_identical(p$xformla, xf)
    expect_identical(environment(p$xformla), fenv)
    expect_null(p$prepared)
    dfl_expect_no_preparation_frame(p, paste("dated covariate", est))
  }
  # Caller data are unchanged by preparation.
  expect_identical(d, dfl_panel(seed = 13L))
})

# ---------------------------------------------------------------------------
# 4. Default formula on supplied-frame, unit-weight and dated-weight routes
# ---------------------------------------------------------------------------

test_that("the default formula on supplied-frame, unit-weight and dated-weight routes stays in the base environment and reaches no preparation frame", {
  d <- dfl_panel(seed = 11L)
  fr <- data.frame(id = rev(unique(d$id)), G = rev(d$G[!duplicated(d$id)]), x = rev(d$x[!duplicated(d$id)]))
  for (est in c("ipw", "dr")) {
    p <- dfl_prepare(d, est = est, frame = fr)
    expect_identical(environment(p$xformla), baseenv())
    expect_identical(deparse(p$xformla), "~1")
    expect_identical(p$frame_ids, fr$id)
    dfl_expect_no_preparation_frame(p, paste("supplied frame", est))
    if (identical(est, "ipw"))
      expect_identical(p$frame_Z, matrix(1, nrow(fr), 1L, dimnames = list(NULL, "(Intercept)")))
  }
  # Unit weights (all one) keep the static law and the record with the default formula.
  p1 <- dfl_prepare(d, est = "reg", weightsname = "one")
  expect_identical(environment(p1$xformla), baseenv())
  expect_identical(p1$weight_reduction, "static")
  expect_s3_class(p1$prepared, "cdid_request_prepared_intake")
  dfl_expect_no_preparation_frame(p1, "unit weights reg")
  # Genuinely changing data weights with an explicit anchor (existing pattern):
  # not the static law, no record, default formula still in the base environment.
  dw <- d; dw$w[dw$period == 3L] <- dw$w[dw$period == 3L] + 1
  p2 <- api$chained_request_prepare_inputs(data = dw, yname = "y", tname = "period", idname = "id", gname = "G",
    weightsname = "w", fix_weights = "first_period")
  expect_identical(environment(p2$xformla), baseenv())
  expect_identical(deparse(p2$xformla), "~1")
  expect_false(identical(p2$weight_reduction, "static"))
  expect_null(p2$prepared)
  dfl_expect_no_preparation_frame(p2, "dated weights first_period")
  # The same changing weights without an anchor keep their existing typed refusal.
  expect_error(api$chained_request_prepare_inputs(data = dw, yname = "y", tname = "period", idname = "id", gname = "G",
    weightsname = "w"), class = "cdid_chained_request_weight_error")
})

# ---------------------------------------------------------------------------
# 5. Genuine public no-X results (IPW/OR/DR)
# ---------------------------------------------------------------------------

test_that("genuine public no-X IPW/OR/DR results reach no preparation frame, carry the base-environment default in DIDparams and equal the explicit ~1 fit", {
  skip_if_not_installed("did")
  d <- dfl_panel(N = 36L, T = 5L, seed = 31L)
  fenv <- new.env(parent = baseenv())
  one <- eval(quote(~1), fenv)
  fit <- function(est, xformla) cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d,
    xformla = xformla, est_method = est, control_group = "nevertreated", bstrap = FALSE, cband = FALSE,
    base_period = "varying", mode = "chained", pt_assumption = "post")
  for (est in c("ipw", "reg", "dr")) {
    mp0 <- fit(est, NULL)
    mp1 <- fit(est, one)
    expect_s3_class(mp0, "cdid_MP")
    expect_s3_class(mp1, "cdid_MP")
    # DIDparams carries the package default in the base environment, never a preparation frame.
    expect_s3_class(mp0$DIDparams$xformla, "formula")
    expect_identical(environment(mp0$DIDparams$xformla), baseenv())
    expect_identical(deparse(mp0$DIDparams$xformla), "~1")
    dfl_expect_no_preparation_frame(mp0, paste("public MP", est, "no X"))
    # The caller's ~1 is retained as the caller's own object.
    expect_identical(mp1$DIDparams$xformla, one)
    expect_identical(environment(mp1$DIDparams$xformla), fenv)
    # Identical science: IDs, point estimates, SEs, full influence functions, covariance.
    expect_identical(mp0$group, mp1$group); expect_identical(mp0$t, mp1$t)
    expect_identical(mp0$att, mp1$att)
    expect_identical(mp0$se, mp1$se)
    expect_identical(mp0$inffunc, mp1$inffunc)
    expect_identical(mp0$V_analytical, mp1$V_analytical)
    expect_true(any(is.finite(mp0$att) & is.finite(mp0$se)))
    expect_identical(mp0$chained$recoding, mp1$chained$recoding)
    expect_identical(mp0$chained$frame, mp1$chained$frame)
    # Serialized state: with no captured preparation frame the no-X result is no
    # larger than the explicit-~1 result (whose formula carries one small empty
    # environment); the saved pre-fix evidence was ~40x larger.
    expect_lte(length(serialize(mp0, NULL)), length(serialize(mp1, NULL)) + 1024L)
  }
  # Caller data are unchanged by the public requests.
  expect_identical(d, dfl_panel(N = 36L, T = 5L, seed = 31L))
})
