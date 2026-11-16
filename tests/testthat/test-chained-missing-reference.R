# Issue #16: att_gt(missing_reference = "error" | "exclude"). Installed-namespace
# tests against the independent base-R reference reference/issue16_reference.R
# (masked link points, identity-GMM cells, D61 masses with eligibility and
# central-difference original-unit influence functions; no cdid code).

i16_installed <- function() isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE)

i16_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  get(name, ns, inherits = FALSE)
}

i16_ref <- new.env(parent = baseenv())
sys.source(testthat::test_path("reference", "issue16_reference.R"), envir = i16_ref, keep.source = FALSE)

# Deterministic panel (no RNG): dates 1..4, cohorts 3 (b = 2), 4 (b = 3) and
# never treated; x changes every date; weights change for even ids; NA outcomes
# keep their covariate and weight rows.
i16_panel <- function(n = 60, dates = 1:4, cohorts = c(3, 4, 0)) {
  id <- seq_len(n); g <- rep(cohorts, length.out = n)
  d <- expand.grid(time = dates, id = id)[, c("id", "time")]
  ti <- match(d$time, dates)
  d$G <- g[d$id]
  d$x <- round(sin(1.7 * d$id) + 0.25 * ti + 0.3 * cos(2.3 * d$id + 1.1 * ti), 2)
  bw <- round(0.6 + (d$id %% 7) / 7, 3)
  d$w <- round(bw * (1 + 0.15 * ti * ((d$id %% 2) == 0)), 3)
  te <- ifelse(d$G > 0 & d$time >= d$G, 1 + 0.5 * (d$time - d$G), 0)
  d$y <- round(cos(d$id) + 0.4 * d$time + 0.8 * d$x + te + 0.5 * sin(5.1 * d$id + 2.7 * ti), 3)
  d$y[((d$id * 37 + ti * 11) %% 13) < 2] <- NA
  d
}

i16_fit <- function(d, missing_reference = NULL, est_method = "ipw", fix_weights = "base_period", ...) {
  args <- list(yname = "y", tname = "time", idname = "id", gname = "G", data = d, xformla = ~x, weightsname = "w",
    est_method = est_method, fix_weights = fix_weights, bstrap = FALSE, cband = FALSE, gmm_weighting = "identity", ...)
  if (!is.null(missing_reference)) args$missing_reference <- missing_reference
  if (identical(missing_reference, "exclude")) suppressWarnings(do.call(i16_ns("att_gt"), args)) else
    do.call(i16_ns("att_gt"), args)
}

i16_post <- function(fit) which(fit$t >= fit$group)

i16_close <- function(a, b, tol) {
  testthat::expect_true(all(is.finite(a)) && all(is.finite(b)))
  testthat::expect_lte(max(abs(a - b)), tol * max(1, max(abs(b))))
}

# Link points, post cells and their complete original-unit IF against the reference.
i16_compare <- function(fit, d, fix_weights = "base_period", eps = 1e-4) {
  S <- i16_ref$i16_setup(d, fix_weights = fix_weights)
  R <- i16_ref$i16_if(S, eps = eps)
  links <- fit$chained$nuisance_diagnostics$fitted_catalog
  cat_rows <- which(links$role == "catalog")
  want <- paste(vapply(S$links, `[[`, numeric(1), "g"), S$cal[vapply(S$links, `[[`, integer(1), "is")],
    S$cal[vapply(S$links, `[[`, integer(1), "it")])
  m <- match(paste(links$cohort[cat_rows], links$s[cat_rows], links$t[cat_rows]), want)
  testthat::expect_false(anyNA(m))
  testthat::expect_setequal(m, seq_along(want))
  testthat::expect_equal(links$mu[cat_rows], R$base$mu[m], tolerance = 1e-10)
  post <- i16_post(fit)
  mc <- match(paste(fit$group[post], fit$t[post]), paste(R$base$cells$g, R$base$cells$t))
  av <- fit$chained$available[post]
  testthat::expect_false(anyNA(mc[av]))
  testthat::expect_equal(unname(fit$att[post][av]), R$base$cells$att[mc[av]], tolerance = 1e-10)
  rows <- match(fit$chained$frame$id, S$ids)
  i16_close(unname(fit$inffunc[, post[av], drop = FALSE]), R$phi_cells[rows, mc[av], drop = FALSE], 3e-7)
  invisible(list(S = S, R = R, rows = rows))
}

i16_agg_compare <- function(fit, R, rows, na.rm = FALSE) {
  agg <- i16_ns("aggte")
  types <- c("simple", "group", "dynamic", "calendar")
  for (j in seq_along(types)) {
    a <- agg(fit, type = types[j], na.rm = na.rm)
    testthat::expect_equal(a$overall.att, unname(R$base$agg[types[j]]), tolerance = 1e-10)
    inf <- switch(types[j], simple = a$inf.function$simple.att, group = a$inf.function$selective.inf.func,
      dynamic = a$inf.function$dynamic.inf.func, calendar = a$inf.function$calendar.inf.func)
    i16_close(as.numeric(inf), R$phi_agg[rows, j], 3e-7)
  }
}

testthat::test_that("API: default warning, strict error and validation before data", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  att <- i16_ns("att_gt")
  testthat::expect_identical(formals(att)$missing_reference, "exclude")
  testthat::expect_identical(utils::tail(names(formals(att)), 1L), "missing_reference")
  testthat::expect_true("missing_reference" %in% i16_ns("cdid_chained_only_formals"))
  base <- list(yname = "y", tname = "time", idname = "id", gname = "G", data = quote(stop("data must not be evaluated")),
    est_method = "ipw")
  for (bad in list("drop", NA_character_, c("error", "exclude"), TRUE)) {
    e <- testthat::expect_error(do.call(att, c(base, list(missing_reference = bad))),
      class = "cdid_chained_request_missing_reference_error")
    testthat::expect_no_match(conditionMessage(e), "data must not be evaluated")
  }
  for (ok in c("error", "exclude"))
    testthat::expect_error(do.call(att, c(base, list(missing_reference = ok))), "data must not be evaluated")
  # The time-invariant preset and switching models refuse the opt-in before data.
  testthat::expect_error(do.call(att, c(base, list(missing_reference = "exclude", control = list(preset = "cohort_time_fe")))),
    class = "cdid_chained_request_missing_reference_error")
  sb <- list(yname = "y", tname = "time", idname = "id", dname = "D", data = quote(stop("data must not be evaluated")))
  testthat::expect_error(do.call(att, c(sb, list(missing_reference = "exclude"))), class = "cdid_static_routing_error")
  testthat::expect_error(do.call(att, c(sb, list(missing_reference = "error"))), "data must not be evaluated")
  # Ordinary mode neither evaluates nor forwards the chained-only argument.
  od <- list(yname = "y", tname = "time", idname = "id", gname = "G", data = quote(stop("data must not be evaluated")),
    mode = "ordinary", missing_reference = quote(stop("missing_reference must not be evaluated")))
  e <- tryCatch(do.call(att, od), error = identity)
  testthat::expect_no_match(conditionMessage(e), "missing_reference must not be evaluated")
  # The default reports the narrowed sample. Strict mode still identifies the missing value.
  d <- i16_panel()
  u <- d$id[d$G == 0 & d$time == 2 & !is.na(d$y)][1L]
  dx <- d; dx$x[dx$id == u & dx$time == 2] <- NA
  e <- testthat::expect_error(i16_fit(dx, missing_reference = "error"), class = "cdid_reference_missing_error")
  testthat::expect_match(conditionMessage(e), paste0("id ", u, ", date 2, variable `x`"), fixed = TRUE)
  testthat::expect_warning(f <- i16_fit(dx), "Excluded units with missing reference covariates", fixed = TRUE)
  testthat::expect_identical(f$chained$reference_eligibility$status, "restricted")
  testthat::expect_identical(nrow(f$inffunc), length(unique(dx$id)))
})

testthat::test_that("unrestricted parity: complete reference values keep the existing path exactly", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  agg <- i16_ns("aggte")
  d <- i16_panel()
  for (em in list("ipw", "dr")) {
    a <- i16_fit(d, est_method = em, missing_reference = "error")
    b <- i16_fit(d, est_method = em, missing_reference = "exclude")
    testthat::expect_identical(b$att, a$att)
    testthat::expect_identical(b$inffunc, a$inffunc)
    testthat::expect_identical(b$V_analytical, a$V_analytical)
    testthat::expect_identical(b$chained$reference_inputs, a$chained$reference_inputs)
    testthat::expect_null(a$chained$reference_eligibility)
    testthat::expect_null(b$chained$reference_eligibility)
    for (tp in c("simple", "dynamic")) testthat::expect_identical(agg(b, type = tp)$overall.att, agg(a, type = tp)$overall.att)
  }
  i16_compare(i16_fit(d, missing_reference = "exclude"), d)
})

testthat::test_that("a never-treated unit missing x at year 2 but known at year 3 enters only target 4's context", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  d <- i16_panel()
  u <- d$id[d$G == 0 & d$time == 2 & !is.na(d$y) & d$id %in% d$id[d$time == 3 & !is.na(d$y)]][1L]
  d$x[d$id == u & d$time == 2] <- NA
  testthat::expect_false(is.na(d$x[d$id == u & d$time == 3]))
  fit <- i16_fit(d, missing_reference = "exclude")
  el <- fit$chained$reference_eligibility
  testthat::expect_identical(el$status, "restricted")
  testthat::expect_equal(el$targets$n_never_treated_excluded[el$targets$g == 3], 1)
  testthat::expect_equal(el$targets$n_never_treated_excluded[el$targets$g == 4], 0)
  testthat::expect_equal(as.numeric(el$restricted_targets), 3)
  cmp <- i16_compare(fit, d)
  # Identity GMM decouples targets: the unit's IF row vanishes for every
  # target-3 cell and is active for target 4; its row and id are retained.
  post <- i16_post(fit); iu <- match(u, fit$chained$frame$id)
  testthat::expect_identical(nrow(fit$inffunc), length(unique(d$id)))
  testthat::expect_lt(max(abs(fit$inffunc[iu, post[fit$group[post] == 3]])), 1e-12)
  testthat::expect_true(any(fit$inffunc[iu, post[fit$group[post] == 4]] != 0))
  # Target 4 keeps its date-3 recipe. Its date-2 pretrend comparison cannot
  # use the unknown covariate because that unit is eligible in the date-3 context.
  rt <- fit$chained$reference_inputs$recipes$recipes
  testthat::expect_true("eligibility_reference_date" %in% names(rt))
  testthat::expect_true(any(rt$design_date == 3 & is.na(rt$eligibility_reference_date)))
  a <- el$auxiliary_unavailable
  testthat::expect_true(any(a$g == 4L & a$s == 2L))
  testthat::expect_true(any(rt$design_date == 2 & rt$eligibility_reference_date %in% 2))
  testthat::expect_output(print(fit), "Restricted comparisons for cohort", fixed = TRUE)
  # Negative control: the unrestricted reference (unit kept in E_3) differs.
  d_full <- i16_panel()
  R0 <- i16_ref$i16_eval(i16_ref$i16_setup(d_full))
  testthat::expect_gt(max(abs(R0$cells$att[R0$cells$g == 3] - cmp$R$base$cells$att[cmp$R$base$cells$g == 3])), 1e-6)
})

testthat::test_that("known reference covariate-only rows keep a unit; an absent reference row excludes it", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  d <- i16_panel()
  u <- d$id[d$G == 0 & d$time == 2 & is.na(d$y)][1L]
  testthat::expect_false(is.na(u))
  # Y = NA with x supplied at b_3 still supplies the required covariate.
  a <- i16_fit(d, missing_reference = "error"); b <- i16_fit(d, missing_reference = "exclude")
  testthat::expect_identical(b$att, a$att)
  testthat::expect_null(b$chained$reference_eligibility)
  # The same unit without its date-2 row: unknown at b_3, excluded from E_3 only.
  dd <- d[!(d$id == u & d$time == 2), ]
  fit <- i16_fit(dd, missing_reference = "exclude")
  el <- fit$chained$reference_eligibility
  testthat::expect_identical(el$status, "restricted")
  testthat::expect_identical(el$contexts$n_excluded[el$contexts$reference_date == 2], 1L)
  testthat::expect_false(any(el$targets$g == 4L & el$targets$restricted))
  i16_compare(fit, dd)
})

testthat::test_that("restricted cohort shares use m_g = H(V I_g E_g) with their IF channel", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  d <- i16_panel()
  tr <- d$id[d$G == 3 & d$time == 2 & !is.na(d$y)][1:2]
  nt <- d$id[d$G == 0 & d$time == 2 & !is.na(d$y)][1L]
  d$x[d$id %in% c(tr, nt) & d$time == 2] <- NA
  fit <- i16_fit(d, missing_reference = "exclude")
  el <- fit$chained$reference_eligibility
  testthat::expect_equal(el$targets$n_treated_eligible[el$targets$g == 3], el$targets$n_treated[el$targets$g == 3] - 2)
  cmp <- i16_compare(fit, d)
  i16_agg_compare(fit, cmp$R, cmp$rows)
  testthat::expect_equal(fit$chained$aggregation_state$share_eligibility$cohort, 3)
  # Negative control: the same cells with unrestricted masses H(V I_g) differ.
  S <- cmp$S; cells <- cmp$R$base$cells
  m_full <- vapply(S$targets, function(g) sum(S$V[S$G == g]) / S$N, numeric(1))
  mg <- m_full[match(cells$g, S$targets)]
  testthat::expect_gt(abs(sum(mg * cells$att) / sum(mg) - cmp$R$base$agg[["simple"]]), 1e-8)
})

testthat::test_that("a target whose eligible support vanishes is unavailable with its reason", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  agg <- i16_ns("aggte")
  d <- i16_panel()
  d$x[d$G == 4 & d$time == 3] <- NA
  fit <- i16_fit(d, missing_reference = "exclude")
  post <- i16_post(fit); g4 <- post[fit$group[post] == 4]
  testthat::expect_true(length(g4) > 0L)
  testthat::expect_false(any(fit$chained$available[g4]))
  testthat::expect_true(all(is.na(fit$att[g4])))
  testthat::expect_true(all(grepl("eligibility", fit$chained$support_reason[g4], fixed = TRUE)))
  testthat::expect_equal(fit$chained$reference_eligibility$targets$n_treated_eligible[
    fit$chained$reference_eligibility$targets$g == 4], 0)
  testthat::expect_true(all(fit$chained$available[post[fit$group[post] == 3]]))
  cmp <- i16_compare(fit, d)
  testthat::expect_true(is.na(agg(fit, type = "simple")$overall.att))
  a <- agg(fit, type = "simple", na.rm = TRUE)
  testthat::expect_equal(a$overall.att, unname(cmp$R$base$agg["simple"]), tolerance = 1e-10)
})

testthat::test_that("auxiliary-only missing values mark the affected pretrend reports unavailable; E_g is unchanged", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  d <- i16_panel()
  u <- d$id[d$G == 0 & d$time == 1 & !is.na(d$y) & d$id %in% d$id[d$time == 2 & !is.na(d$y)]][1L]
  dx <- d; dx$x[dx$id == u & dx$time == 1] <- NA # date 1 is no target's reference date
  testthat::expect_error(i16_fit(dx, missing_reference = "error"), class = "cdid_reference_missing_error")
  full <- i16_fit(d)
  fit <- i16_fit(dx, missing_reference = "exclude")
  testthat::expect_identical(fit$chained$reference_eligibility$status, "auxiliary_unavailable")
  post <- i16_post(fit); pf <- match(paste(fit$group[post], fit$t[post]), paste(full$group, full$t))
  testthat::expect_identical(unname(fit$att[post]), unname(full$att[pf]))
  testthat::expect_identical(unname(fit$inffunc[, post]), unname(full$inffunc[, pf]))
  au <- fit$chained$reference_eligibility$auxiliary_unavailable
  testthat::expect_gt(nrow(au), 0L)
  testthat::expect_true(all(au$s == 1))
  ds <- fit$chained$nuisance_diagnostics$auxiliary_descriptor_status
  testthat::expect_true(all(ds$status[au$descriptor] == "unsupported_pair"))
  testthat::expect_true(all(grepl("missing_reference", ds$reason[au$descriptor], fixed = TRUE)))
  # Full joint test rule: no joint pretrend test with an unavailable report.
  testthat::expect_false(identical(fit$chained$pretest$status, "computed"))
  # Accepted selection case: the same auxiliary-only gap under conditional MTAR.
  sel <- i16_fit(dx, missing_reference = "exclude", selection = "conditional_mtar")
  testthat::expect_gt(nrow(sel$chained$reference_eligibility$auxiliary_unavailable), 0L)
  testthat::expect_true(all(is.finite(sel$att[i16_post(sel)])))
  testthat::expect_error(i16_fit(dx, missing_reference = "error", selection = "conditional_mtar"), class = "cdid_reference_missing_error")
})

testthat::test_that("restricted custom callbacks and first-date anchors read the same contexts", {
  testthat::skip_if_not(i16_installed(), "installed cdid required")
  d <- i16_panel()
  u <- d$id[d$G == 0 & d$time == 2 & !is.na(d$y)][1L]
  d$x[d$id == u & d$time == 2] <- NA
  i16_compare(i16_fit(d, missing_reference = "exclude", fix_weights = "first_period"), d, fix_weights = "first_period")
  custom <- function(y1, y0, D, covariates, i.weights = NULL, inffunc = TRUE, ...) {
    w <- i.weights / mean(i.weights); dy <- y1 - y0
    m1 <- sum(w * D * dy) / sum(w * D); m0 <- sum(w * (1 - D) * dy) / sum(w * (1 - D))
    list(ATT = m1 - m0, att.inf.func = w * D * (dy - m1) / mean(w * D) - w * (1 - D) * (dy - m0) / mean(w * (1 - D)))
  }
  fit <- i16_fit(d, missing_reference = "exclude", est_method = custom)
  post <- i16_post(fit); iu <- match(u, fit$chained$frame$id)
  testthat::expect_lt(max(abs(fit$inffunc[iu, post[fit$group[post] == 3]])), 1e-12)
  testthat::expect_identical(fit$chained$reference_eligibility$status, "restricted")
})
