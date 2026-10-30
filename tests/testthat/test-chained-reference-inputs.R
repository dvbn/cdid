# CP08 D60/D61 reference-date covariates, fixed-anchor changing weights and
# first-calendar aggregation masses. Installed-namespace tests: the public
# att_gt()/aggte() results are compared with the independent base-R references
# reference/d60_or_reference.R (shared OR) and reference/d60_routes_reference.R
# (baseline IPW/custom/DR, conditional-selection IPW/reg/DR, pure imputation,
# augmented observation DR; no cdid code), whose influence functions are
# central differences of original-unit case weights at steps 1e-4 and 5e-5.
# Static calls must stay identical under every accepted fix_weights value.

rid_installed <- function() isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE)

rid_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  get(name, ns, inherits = FALSE)
}

rid_ref <- new.env(parent = baseenv())
sys.source(testthat::test_path("reference", "d60_or_reference.R"), envir = rid_ref, keep.source = FALSE)
sys.source(testthat::test_path("reference", "d60_routes_reference.R"), envir = rid_ref, keep.source = FALSE)

# Deterministic panel (no RNG): irregular calendar, changing x, weights that
# change for every even id (all cohorts), and Y = NA rows that keep x and w.
rid_panel <- function(n = 48, dates = c(1, 2, 4, 5, 7), cohorts = c(4, 5, 0),
                      vary_x = TRUE, vary_w = TRUE, drop_y = TRUE) {
  id <- seq_len(n); g <- rep(cohorts, length.out = n)
  d <- expand.grid(time = dates, id = id)[, c("id", "time")]
  ti <- match(d$time, dates)
  d$G <- g[d$id]
  bx <- round(sin(1.7 * id), 2)
  d$x <- if (vary_x) round(bx[d$id] + 0.25 * ti + 0.3 * cos(2.3 * d$id + 1.1 * ti), 2) else bx[d$id]
  bw <- round(0.6 + (id %% 7) / 7, 3)
  d$w <- if (vary_w) round(bw[d$id] * (1 + 0.15 * ti * ((d$id %% 2) == 0)), 3) else bw[d$id]
  te <- ifelse(d$G > 0 & d$time >= d$G, 1 + 0.5 * (d$time - d$G), 0)
  d$y <- round(cos(d$id) + 0.4 * d$time + 0.8 * d$x + te + 0.5 * sin(5.1 * d$id + 2.7 * ti), 3)
  if (drop_y) d$y[(d$id * 7 + ti * 3) %% 9 == 0] <- NA
  d
}

rid_fit <- function(d, xformla = ~x, est_method = "reg", control_group = "nevertreated", ...) {
  rid_ns("att_gt")(yname = "y", tname = "time", idname = "id", gname = "G", data = d, xformla = xformla,
    weightsname = "w", est_method = est_method, control_group = control_group, bstrap = FALSE, cband = FALSE,
    gmm_weighting = "identity", ...)
}

rid_post <- function(fit) which(fit$t >= fit$group)

# |a - b| <= tol * max(1, max|b|) entrywise maximum.
rid_close <- function(a, b, tol) {
  testthat::expect_true(all(is.finite(a)) && all(is.finite(b)))
  testthat::expect_lte(max(abs(a - b)), tol * max(1, max(abs(b))))
}

rid_compare <- function(fit, S, steps = c(1e-4, 5e-5)) {
  post <- rid_post(fit)
  for (h in steps) {
    R <- rid_ref$d60_ref_if(S, eps = h)
    m <- match(paste(fit$group[post], fit$t[post]), paste(R$base$cells$g, R$base$cells$t))
    testthat::expect_false(anyNA(m))
    testthat::expect_equal(unname(fit$att[post]), R$base$cells$att[m], tolerance = 1e-10)
    rid_close(unname(fit$inffunc[, post, drop = FALSE]), R$phi_cells[, m, drop = FALSE], 3e-7)
  }
  invisible(R)
}

# Route fixture (no RNG): every cohort misses outcomes at every date, so the
# selection and treatment logits have finite roots; rho is a fixed known unit
# ratio (constant within unit, read from the first-observed frame row).
rid_route_panel <- function(n = 36, dates = c(1, 2, 4, 5), cohorts = c(4, 5, 0)) {
  d <- rid_panel(n = n, dates = dates, cohorts = cohorts, drop_y = FALSE)
  ti <- match(d$time, dates)
  d$y[((d$id * 37 + ti * 11) %% 13) < 3] <- NA
  d$rho <- round(0.6 + 0.3 * (d$id %% 4) / 4, 3)
  d
}

# did-style callback with an exact influence function (regression adjustment
# on the supplied covariates, DRDID reg_did_panel form): the custom route must
# hand it z(C^r) at each context's recipe date and i.weights = W_g / r_G.
rid_custom_reg <- function(y1, y0, D, covariates, i.weights = NULL, inffunc = TRUE, ...) {
  w <- i.weights / mean(i.weights); dy <- y1 - y0; X <- as.matrix(covariates)
  c0 <- (1 - D) * w
  M <- crossprod(X, X * c0) / length(D)
  beta <- solve(crossprod(X, X * c0), crossprod(X, c0 * dy))
  fit <- drop(X %*% beta)
  p1 <- mean(w * D)
  att <- sum(w * D * (dy - fit)) / sum(w * D)
  mx1 <- colMeans(w * D * X) / p1
  list(ATT = att, att.inf.func = w * D * (dy - fit - att) / p1 - drop((c0 * (dy - fit) * X) %*% solve(M, mx1)))
}

# Public route arguments for one reference route specification.
rid_route_args <- function(rt) {
  sel <- list(selection = "conditional_mtar", selection_model = rt$sel_model, selection_pooling = rt$sel_pooling)
  switch(rt$kind,
    baseline = list(est_method = if (identical(rt$est, "custom")) rid_custom_reg else rt$est),
    corollary1 = c(list(est_method = rt$est), sel, if (!is.null(rt$universe) && rt$universe != "F")
      list(selection_universe = rt$universe, selection_rho = "rho",
        selection_justification = "fixture: known unit-level normalized ratio rho")),
    imputation = list(est_method = rt$est, selection = "conditional_mtar", selection_est_method = "imputation",
      imputation_pooling = rt$imp_pooling),
    augmented = c(list(est_method = "dr", selection_est_method = "dr", imputation_pooling = rt$imp_pooling), sel))
}

# Post cells, complete original-unit IF and the D61 simple/dynamic aggregates
# of one public fit against the independent route reference.
rid_route_compare <- function(fit, S, steps = c(1e-4, 5e-5)) {
  agg <- rid_ns("aggte"); post <- rid_post(fit)
  for (h in steps) {
    R <- rid_ref$d60r_if(S, eps = h)
    m <- match(paste(fit$group[post], fit$t[post]), paste(R$base$cells$g, R$base$cells$t))
    testthat::expect_false(anyNA(m))
    testthat::expect_equal(unname(fit$att[post]), R$base$cells$att[m], tolerance = 1e-10)
    rid_close(unname(fit$inffunc[, post, drop = FALSE]), R$phi_cells[, m, drop = FALSE], 3e-7)
    for (tp in c("simple", "dynamic")) {
      a <- agg(fit, type = tp)
      testthat::expect_equal(a$overall.att, unname(R$base$agg[tp]), tolerance = 1e-10)
      rid_close(as.numeric(if (tp == "simple") a$inf.function$simple.att else a$inf.function$dynamic.inf.func),
        R$phi_agg[, match(tp, c("simple", "group", "dynamic", "calendar"))], 3e-7)
    }
  }
  invisible(R)
}

testthat::test_that("fix_weights contract: varying refused before data; first/base accepted; ordinary mode unchanged", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  att <- rid_ns("att_gt")
  base <- list(yname = "y", tname = "time", idname = "id", gname = "G", data = quote(stop("data must not be evaluated")),
    est_method = "reg")
  e <- testthat::expect_error(do.call(att, c(base, list(fix_weights = "varying"))),
    class = "cdid_chained_request_unsupported_weights")
  testthat::expect_no_match(conditionMessage(e), "data must not be evaluated")
  testthat::expect_match(conditionMessage(e), "deferred", fixed = TRUE)
  for (bad in list("anchor", NA_character_, c("first_period", "base_period"), 1))
    testthat::expect_error(do.call(att, c(base, list(fix_weights = bad))), class = "cdid_chained_request_unsupported_weights")
  for (ok in list(NULL, "first_period", "base_period"))
    testthat::expect_error(do.call(att, c(base, list(fix_weights = ok))), "data must not be evaluated")
  testthat::expect_identical(rid_ns("upstream_att_gt_formals")[12L], "fix_weights")
})

testthat::test_that("static calls are identical under NULL, base_period and first_period on every route", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  agg <- rid_ns("aggte")
  d <- rid_panel(vary_x = FALSE, vary_w = FALSE)
  custom <- function(y1, y0, D, covariates, i.weights = NULL, inffunc = TRUE, ...) {
    w <- i.weights / mean(i.weights); dy <- y1 - y0
    m1 <- sum(w * D * dy) / sum(w * D); m0 <- sum(w * (1 - D) * dy) / sum(w * (1 - D))
    list(ATT = m1 - m0, att.inf.func = w * D * (dy - m1) / mean(w * D) - w * (1 - D) * (dy - m0) / mean(w * (1 - D)))
  }
  routes <- list(ipw = list(est_method = "ipw"), reg = list(est_method = "reg"), dr = list(est_method = "dr"),
    custom = list(est_method = custom))
  for (nm in names(routes)) {
    fits <- lapply(list(NULL, "base_period", "first_period"), function(f) do.call(rid_fit, c(list(d), routes[[nm]], list(fix_weights = f))))
    for (k in 2:3) {
      testthat::expect_identical(fits[[k]]$att, fits[[1L]]$att)
      testthat::expect_identical(fits[[k]]$inffunc, fits[[1L]]$inffunc)
      testthat::expect_identical(fits[[k]]$V_analytical, fits[[1L]]$V_analytical)
      testthat::expect_identical(fits[[k]]$chained$frame$weight, fits[[1L]]$chained$frame$weight)
      for (tp in c("simple", "dynamic")) {
        a1 <- agg(fits[[1L]], type = tp); ak <- agg(fits[[k]], type = tp)
        testthat::expect_identical(ak$overall.att, a1$overall.att)
        testthat::expect_identical(ak$overall.se, a1$overall.se)
      }
    }
    testthat::expect_null(fits[[1L]]$chained$reference_inputs)
    testthat::expect_identical(fits[[2L]]$chained$reference_inputs$reduction, "static")
  }
})

testthat::test_that("unit-level x with first_period is the existing static system on w at the first calendar date", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_panel(vary_x = FALSE, vary_w = TRUE)
  d1 <- d; d1$w <- d$w[d$time == 1][match(d$id, d$id[d$time == 1])]
  for (em in c("ipw", "reg", "dr")) {
    a <- rid_fit(d, est_method = em, fix_weights = "first_period")
    b <- rid_fit(d1, est_method = em)
    testthat::expect_identical(a$att, b$att)
    testthat::expect_identical(a$inffunc, b$inffunc)
    testthat::expect_identical(a$chained$reference_inputs$reduction, "first_period_static")
  }
  # NULL is did's per-cell weight and stays refused for changing weights.
  e <- testthat::expect_error(rid_fit(d, est_method = "reg"), class = "cdid_chained_request_unsupported_weights")
  testthat::expect_match(conditionMessage(e), "first_period", fixed = TRUE)
})

testthat::test_that("dated OR matches the independent reference: both anchors, NT/NYT, anticipation, irregular calendar, PT-all", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_panel()
  cases <- list(
    list(fix = "base_period", cg = "nevertreated", ant = 0, f = ~x),
    list(fix = "first_period", cg = "notyettreated", ant = 0, f = ~x),
    list(fix = "base_period", cg = "notyettreated", ant = 1, f = ~x),
    list(fix = "first_period", cg = "nevertreated", ant = 0, f = ~x + a),
    list(fix = "base_period", cg = "notyettreated", ant = 0, f = ~scale(x)),
    list(fix = "base_period", cg = "nevertreated", ant = 0, f = ~x, pt = "all"))
  d$a <- round(sin(d$id), 2)
  for (cs in cases) {
    dd <- if (cs$ant > 0) rid_panel(n = 60, dates = c(1, 2, 4, 5, 7, 8), cohorts = c(5, 7, 0)) else d
    if (cs$ant > 0) dd$a <- round(sin(dd$id), 2)
    extra <- if (identical(cs$pt, "all")) list(pt_assumption = "all", pt_start = min(dd$time)) else list()
    fit <- do.call(rid_fit, c(list(dd, xformla = cs$f, control_group = cs$cg, anticipation = cs$ant, fix_weights = cs$fix), extra))
    S <- rid_ref$d60_ref_setup(dd, formula = cs$f, fix_weights = cs$fix, control_group = cs$cg, anticipation = cs$ant,
      pt_assumption = if (is.null(cs$pt)) "post" else cs$pt)
    R <- rid_compare(fit, S)
    rec <- fit$chained$reference_inputs
    testthat::expect_identical(rec$reduction, "reference_recipes")
    testthat::expect_identical(rec$share_law, "first_calendar_masses")
    # Post links of target g extract at b_g; weights at t_1 or b_g for the whole context.
    rt <- rec$recipes$recipes; vr <- rec$recipes$view_recipe
    links <- fit$chained$nuisance_diagnostics$fitted_catalog
    cal <- sort(unique(dd$time))
    for (j in which(links$role == "catalog")) {
      b <- max(cal[cal + cs$ant < links$cohort[j]])
      testthat::expect_equal(rt$design_date[vr[j]], b)
      testthat::expect_equal(rt$weight_anchor[vr[j]], if (cs$fix == "first_period") cal[1L] else b)
    }
    for (j in which(links$role == "auxiliary")) testthat::expect_equal(rt$design_date[vr[j]], links$s[j])
  }
})

testthat::test_that("static designs keep the static design and alias map; one canonical recipe never duplicates blocks", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  ds <- rid_panel(vary_x = FALSE)
  # Unit-level x with base-period weights: one recipe per target weight law.
  fit <- rid_fit(ds, fix_weights = "base_period", control_group = "notyettreated")
  rid_compare(fit, rid_ref$d60_ref_setup(ds, fix_weights = "base_period", control_group = "notyettreated"), steps = 1e-4)
  testthat::expect_true(all(is.na(fit$chained$reference_inputs$recipes$recipes$design_date)))
  # first_period with frame-only units (no first-date weight, never read by these
  # cells): the single recipe equals the static fit, block for block.
  extra <- data.frame(id = c(101L, 102L), G = c(0, 0))
  fr <- rbind(data.frame(id = unique(ds$id), G = ds$G[match(unique(ds$id), ds$id)], x = ds$x[match(unique(ds$id), ds$id)]),
    data.frame(extra, x = NA_real_))
  a <- rid_fit(ds, frame = fr, fix_weights = "first_period", control_group = "notyettreated")
  testthat::expect_identical(a$chained$reference_inputs$reduction, "reference_recipes")
  testthat::expect_identical(nrow(a$chained$reference_inputs$recipes$recipes), 1L)
  dc <- ds; dc$w <- ds$w[ds$time == 1][match(ds$id, ds$id[ds$time == 1])]
  fr_s <- fr; fr_s$w <- c(dc$w[match(unique(ds$id), dc$id)], 1, 1)
  b <- rid_fit(dc, frame = fr_s, control_group = "notyettreated")
  testthat::expect_identical(a$att, b$att)
  testthat::expect_identical(a$inffunc, b$inffunc)
  testthat::expect_identical(a$chained$nuisance_diagnostics$accounting$n_regression_blocks,
    b$chained$nuisance_diagnostics$accounting$n_regression_blocks)
  rid_compare(a, rid_ref$d60_ref_setup(ds, fix_weights = "first_period", control_group = "notyettreated",
    frame_extra = extra), steps = 1e-4)
  # A character covariate measured every period: levels and contrasts freeze on the union.
  dd <- rid_panel(); dd$f_t <- ifelse(dd$x > 1, "hi", "lo")
  fit <- rid_fit(dd, xformla = ~x + f_t, fix_weights = "base_period", control_group = "notyettreated")
  rid_compare(fit, rid_ref$d60_ref_setup(dd, formula = ~x + f_t, fix_weights = "base_period",
    control_group = "notyettreated"), steps = 1e-4)
})

testthat::test_that("D61 aggregation uses first-calendar masses with share channels; never-treated V is not demanded", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  agg <- rid_ns("aggte")
  types <- c("simple", "group", "dynamic", "calendar")
  inf_of <- function(a, tp) as.numeric(switch(tp, simple = a$inf.function$simple.att,
    group = a$inf.function$selective.inf.func, dynamic = a$inf.function$dynamic.inf.func,
    calendar = a$inf.function$calendar.inf.func))
  d <- rid_panel()
  d_nt <- d[!(d$G == 0 & d$id %in% c(3, 6, 9) & d$time == 1), ]
  for (dd in list(d, d_nt)) {
    fit <- rid_fit(dd, fix_weights = "base_period")
    S <- rid_ref$d60_ref_setup(dd, fix_weights = "base_period")
    for (h in c(1e-4, 5e-5)) {
      R <- rid_ref$d60_ref_if(S, eps = h)
      for (j in seq_along(types)) {
        a <- agg(fit, type = types[j])
        testthat::expect_equal(a$overall.att, unname(R$base$agg[types[j]]), tolerance = 1e-10)
        rid_close(inf_of(a, types[j]), R$phi_agg[, j], 3e-7)
      }
    }
  }
  # A treated unit without a first-calendar weight: every aggregate reading its
  # cohort's mass is refused with the id; mass-free bins stay available.
  d_tr <- d[!(d$id == 1 & d$time == 1), ]
  fit <- rid_fit(d_tr, fix_weights = "base_period")
  a <- agg(fit, type = "simple")
  testthat::expect_true(is.na(a$overall.att))
  testthat::expect_match(a$chained$cell$reason[1L], "missing_first_period_weight", fixed = TRUE)
  testthat::expect_match(a$chained$cell$reason[1L], "ids 1", fixed = TRUE)
  g <- agg(fit, type = "group", na.rm = TRUE)
  testthat::expect_true(all(is.finite(g$att.egt)))
  testthat::expect_true(is.na(g$overall.att))
})

testthat::test_that("post estimates are invariant to base_period; equal anchors under constant weights", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_panel()
  a <- rid_fit(d, fix_weights = "base_period", base_period = "varying")
  b <- rid_fit(d, fix_weights = "base_period", base_period = "universal")
  pa <- rid_post(a); pb <- match(paste(a$group[pa], a$t[pa]), paste(b$group, b$t))
  testthat::expect_identical(unname(a$att[pa]), unname(b$att[pb]))
  testthat::expect_identical(unname(a$inffunc[, pa]), unname(b$inffunc[, pb]))
  dc <- rid_panel(vary_w = FALSE)
  e1 <- rid_fit(dc, fix_weights = "first_period"); e2 <- rid_fit(dc, fix_weights = "base_period")
  testthat::expect_identical(e1$att, e2$att)
  testthat::expect_identical(e1$inffunc, e2$inffunc)
})

testthat::test_that("negative controls: wrong anchors, summed X, merged or duplicated blocks, unit loss and omitted channels are detected", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_panel()
  agg <- rid_ns("aggte")
  fit <- rid_fit(d, fix_weights = "base_period", control_group = "notyettreated")
  post <- rid_post(fit)
  gap <- function(variant) {
    S <- rid_ref$d60_ref_setup(d, fix_weights = "base_period", control_group = "notyettreated", variant = variant)
    R <- rid_ref$d60_ref_if(S, eps = 1e-4)
    m <- match(paste(fit$group[post], fit$t[post]), paste(R$base$cells$g, R$base$cells$t))
    list(att = max(abs(unname(fit$att[post]) - R$base$cells$att[m])),
         inf = if (nrow(R$phi_cells) == nrow(fit$inffunc)) max(abs(unname(fit$inffunc[, post]) - R$phi_cells[, m])) else Inf,
         R = R)
  }
  testthat::expect_lt(gap("none")$att, 1e-10)
  for (v in c("own_anchor", "summed_x", "merge_recipes", "duplicate_aux")) testthat::expect_gt(gap(v)$att, 1e-6)
  testthat::expect_gt(gap("frozen_beta")$inf, 1e-4)
  dl <- d[!(d$id %in% c(3, 4) & d$time == 1), ]
  fl <- rid_fit(dl, fix_weights = "base_period", control_group = "notyettreated")
  Sl <- rid_ref$d60_ref_setup(dl, fix_weights = "base_period", control_group = "notyettreated", variant = "drop_no_t1")
  Rl <- rid_ref$d60_ref_eval(Sl)
  pl <- rid_post(fl)
  testthat::expect_gt(max(abs(unname(fl$att[pl]) - Rl$cells$att[match(paste(fl$group[pl], fl$t[pl]), paste(Rl$cells$g, Rl$cells$t))])), 1e-6)
  testthat::expect_identical(nrow(fl$inffunc), length(unique(dl$id)))
  Sf <- rid_ref$d60_ref_setup(d, fix_weights = "base_period", control_group = "notyettreated", variant = "frozen_shares")
  Rf <- rid_ref$d60_ref_if(Sf, eps = 1e-4)
  a <- agg(fit, type = "simple")
  testthat::expect_gt(max(abs(as.numeric(a$inf.function$simple.att) - Rf$phi_agg[, 1L])), 1e-4)
})

testthat::test_that("missing required reference values are refused with ids, never filled; every route now runs", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_panel()
  # A never-treated control observing (2, 4) needs x at b_4 = 2 for target 4.
  nt <- d$id[d$G == 0 & d$time == 2 & !is.na(d$y)]
  u <- nt[nt %in% d$id[d$time == 4 & !is.na(d$y)]][1L]
  dx <- d; dx$x[dx$id == u & dx$time == 2] <- NA
  e <- testthat::expect_error(rid_fit(dx, fix_weights = "base_period"), class = "cdid_reference_missing_error")
  testthat::expect_match(conditionMessage(e), paste0("id ", u, ", date 2, variable `x`"), fixed = TRUE)
  testthat::expect_match(conditionMessage(e), "support (", fixed = TRUE)
  dw <- d; dw$w[dw$id == u & dw$time == 2] <- NA
  e <- testthat::expect_error(rid_fit(dw, fix_weights = "base_period"), class = "cdid_reference_missing_error")
  testthat::expect_match(conditionMessage(e), "sampling weights at the context anchor", fixed = TRUE)
  # Every other route consumes the sealed recipe-indexed reader (no pending refusal).
  for (em in c("ipw", "dr")) {
    a <- rid_fit(d, est_method = em, fix_weights = "first_period")
    testthat::expect_identical(a$chained$reference_inputs$reduction, "reference_recipes")
    testthat::expect_true(all(is.finite(a$att[rid_post(a)])))
    b <- rid_fit(rid_panel(vary_x = FALSE), est_method = em, fix_weights = "base_period")
    testthat::expect_identical(b$chained$reference_inputs$reduction, "reference_recipes")
  }
  # Pure imputation predicts on EVERY roster row of each cohort: x at b_g is
  # demanded even for a unit whose outcomes never form a pair there.
  dr <- rid_route_panel()
  u <- dr$id[dr$G == 0 & dr$time == 2 & is.na(dr$y)][1L]
  dx <- dr; dx$x[dx$id == u & dx$time == 2] <- NA
  e <- testthat::expect_error(do.call(rid_fit, c(list(dx), rid_route_args(list(kind = "imputation", est = "dr",
    imp_pooling = "separate")), list(fix_weights = "base_period"))), class = "cdid_reference_missing_error")
  testthat::expect_match(conditionMessage(e), paste0("id ", u, ", date 2, variable `x`, prediction roster of support ("),
    fixed = TRUE)
  # A separate dated selection formula is demanded on the whole selection roster
  # at the recipe date (its own frozen design; xformla stays complete).
  dz <- dr; dz$z <- round(dz$x - 0.1 * dz$time, 2); dz$z[dz$id == u & dz$time == 2] <- NA
  e <- testthat::expect_error(do.call(rid_fit, c(list(dz), rid_route_args(list(kind = "corollary1", est = "ipw",
    sel_model = "pair_logit", sel_pooling = "separate")), list(fix_weights = "base_period", selection_xformla = ~z))),
    class = "cdid_reference_missing_error")
  testthat::expect_match(conditionMessage(e), paste0("id ", u, ", date 2, variable `z`, selection fitting roster of support ("),
    fixed = TRUE)
  # Frame supply cannot override changing data values.
  fr <- data.frame(id = unique(d$id), G = d$G[match(unique(d$id), d$id)], x = 0)
  testthat::expect_error(rid_fit(d, frame = fr, fix_weights = "base_period"), class = "cdid_chained_request_covariate_error")
})


testthat::test_that("dated constructors use canonical rows when the static constructor fails first", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_panel(n = 24)
  # The final date is never a covariate extraction date in this fixed union.
  # Its absent x must not be passed to poly() or treated as required input.
  d$x[d$time == 7] <- NA_real_
  for (anchor in c("base_period", "first_period")) {
    a <- rid_fit(d, xformla = ~poly(x, 2), fix_weights = anchor)
    S <- rid_ref$d60_ref_setup(d, formula = ~poly(x, 2), fix_weights = anchor)
    rid_compare(a, S, steps = 1e-4)
  }
  # All first-observed frame rows have one level; required dated rows have both
  # levels, with within-block variation. The static contrasts error is not a
  # reason to reject this valid dated design.
  d <- rid_panel(n = 24)
  d$f <- ifelse(d$time == 1 | d$id %% 2 == 0, "low", "high")
  a <- rid_fit(d, xformla = ~x + f, fix_weights = "base_period")
  S <- rid_ref$d60_ref_setup(d, formula = ~x + f, fix_weights = "base_period")
  rid_compare(a, S, steps = 1e-4)
})

testthat::test_that("every dated route matches the independent reference: anchors, NT/NYT, pooling, universes, anticipation, PT-all", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_route_panel()
  cases <- list(
    list(fix = "base_period", cg = "nevertreated", rt = list(kind = "baseline", est = "ipw")),
    list(fix = "first_period", cg = "notyettreated", rt = list(kind = "baseline", est = "custom", custom = rid_custom_reg)),
    list(fix = "base_period", cg = "notyettreated", rt = list(kind = "baseline", est = "dr")),
    list(fix = "first_period", cg = "nevertreated",
      rt = list(kind = "corollary1", est = "ipw", sel_model = "pair_logit", sel_pooling = "separate")),
    list(fix = "base_period", cg = "nevertreated",
      rt = list(kind = "corollary1", est = "reg", sel_model = "pair_logit", sel_pooling = "within_cohort", universe = "B")),
    list(fix = "base_period", cg = "notyettreated",
      rt = list(kind = "corollary1", est = "dr", sel_model = "independent_period_logit", sel_pooling = "across_cohorts")),
    list(fix = "first_period", cg = "notyettreated",
      rt = list(kind = "corollary1", est = "ipw", sel_model = "pair_logit", sel_pooling = "across_cohorts", universe = "A")),
    list(fix = "base_period", cg = "nevertreated", rt = list(kind = "imputation", est = "ipw", imp_pooling = "within_cohort")),
    list(fix = "first_period", cg = "notyettreated", rt = list(kind = "imputation", est = "reg", imp_pooling = "across_cohorts")),
    list(fix = "base_period", cg = "nevertreated", rt = list(kind = "imputation", est = "dr", imp_pooling = "separate")),
    list(fix = "base_period", cg = "nevertreated", rt = list(kind = "augmented", est = "dr", sel_model = "pair_logit",
      sel_pooling = "across_cohorts", imp_pooling = "separate")),
    list(fix = "first_period", cg = "notyettreated", rt = list(kind = "augmented", est = "dr",
      sel_model = "independent_period_logit", sel_pooling = "within_cohort", imp_pooling = "within_cohort")),
    list(fix = "base_period", cg = "notyettreated", ant = 1,
      rt = list(kind = "imputation", est = "dr", imp_pooling = "within_cohort")),
    list(fix = "first_period", cg = "nevertreated", pt = "all",
      rt = list(kind = "corollary1", est = "dr", sel_model = "pair_logit", sel_pooling = "within_cohort")))
  for (cs in cases) {
    ant <- if (is.null(cs$ant)) 0 else cs$ant
    dd <- if (ant > 0) rid_route_panel(dates = c(1, 2, 4, 5, 7), cohorts = c(5, 7, 0)) else d
    rt <- cs$rt
    if (!is.null(rt$universe)) rt$rho <- dd$rho[match(unique(dd$id), dd$id)]
    extra <- if (identical(cs$pt, "all")) list(pt_assumption = "all", pt_start = min(dd$time)) else list()
    fit <- do.call(rid_fit, c(list(dd), rid_route_args(rt),
      list(fix_weights = cs$fix, control_group = cs$cg, anticipation = ant), extra))
    S <- rid_ref$d60r_setup(dd, rt, fix_weights = cs$fix, control_group = cs$cg, anticipation = ant,
      pt_assumption = if (is.null(cs$pt)) "post" else cs$pt)
    rid_route_compare(fit, S)
    rec <- fit$chained$reference_inputs
    testthat::expect_identical(rec$reduction, "reference_recipes")
    testthat::expect_identical(rec$share_law, "first_calendar_masses")
    testthat::expect_identical(nrow(fit$inffunc), length(unique(dd$id)))
  }
})

testthat::test_that("route negative controls: first-date X, own-cohort anchors and cross-recipe pooling are detected", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_route_panel()
  gap <- function(fit, S) {
    post <- rid_post(fit); R <- rid_ref$d60r_eval(S)
    max(abs(unname(fit$att[post]) - R$cells$att[match(paste(fit$group[post], fit$t[post]), paste(R$cells$g, R$cells$t))]))
  }
  for (rt in list(list(kind = "imputation", est = "dr", imp_pooling = "across_cohorts"),
                  list(kind = "corollary1", est = "ipw", sel_model = "pair_logit", sel_pooling = "across_cohorts"))) {
    fit <- do.call(rid_fit, c(list(d), rid_route_args(rt), list(fix_weights = "base_period", control_group = "notyettreated")))
    testthat::expect_lt(gap(fit, rid_ref$d60r_setup(d, rt, control_group = "notyettreated")), 1e-10)
    for (v in c("static_x", "own_anchor", "pool_across_recipes"))
      testthat::expect_gt(gap(fit, rid_ref$d60r_setup(d, rt, control_group = "notyettreated", variant = v)), 1e-6)
  }
})

testthat::test_that("selection, imputation and augmented routes keep static reductions and base_period invariance", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_route_panel()
  first <- function(v) d[[v]][d$time == 1][match(d$id, d$id[d$time == 1])]
  ds <- d; ds$x <- first("x"); ds$w <- first("w")
  dx <- d; dx$x <- first("x")
  routes <- list(
    list(kind = "corollary1", est = "dr", sel_model = "independent_period_logit", sel_pooling = "within_cohort"),
    list(kind = "imputation", est = "dr", imp_pooling = "separate"),
    list(kind = "augmented", est = "dr", sel_model = "pair_logit", sel_pooling = "separate", imp_pooling = "separate"))
  for (rt in routes) {
    args <- rid_route_args(rt)
    fits <- lapply(list(NULL, "base_period", "first_period"), function(f) do.call(rid_fit, c(list(ds), args, list(fix_weights = f))))
    for (k in 2:3) {
      testthat::expect_identical(fits[[k]]$att, fits[[1L]]$att)
      testthat::expect_identical(fits[[k]]$inffunc, fits[[1L]]$inffunc)
    }
    # Unit-level x with first_period is the static system on w at t_1.
    a <- do.call(rid_fit, c(list(dx), args, list(fix_weights = "first_period")))
    b <- do.call(rid_fit, c(list(ds), args))
    testthat::expect_identical(a$chained$reference_inputs$reduction, "first_period_static")
    testthat::expect_identical(a$att, b$att)
    testthat::expect_identical(a$inffunc, b$inffunc)
    # One union registry: post cells do not depend on the reporting base period.
    v <- do.call(rid_fit, c(list(d), args, list(fix_weights = "base_period", base_period = "varying")))
    u <- do.call(rid_fit, c(list(d), args, list(fix_weights = "base_period", base_period = "universal")))
    pv <- rid_post(v); pu <- match(paste(v$group[pv], v$t[pv]), paste(u$group, u$t))
    testthat::expect_identical(unname(v$att[pv]), unname(u$att[pu]))
    testthat::expect_identical(unname(v$inffunc[, pv]), unname(u$inffunc[, pu]))
  }
})

testthat::test_that("dated augmented guards charge live sources and every recipe copy", {
  testthat::skip_if_not(rid_installed(), "installed cdid required")
  d <- rid_route_panel()
  resolved <- rid_ns("chained_public_resolve_control")(NULL, gmm_weighting = "identity")
  prep <- rid_ns("chained_request_prepare_inputs")(data = d, yname = "y", tname = "time", idname = "id", gname = "G",
    xformla = ~x, weightsname = "w", selection = "conditional_mtar", selection_model = "pair_logit",
    selection_pooling = "separate", selection_est_method = "dr", imputation_pooling = "separate",
    fix_weights = "base_period", or_resource = resolved$profile)
  prep$cp04$automatic_chunk_rows <- TRUE
  dims <- rid_ns("chained_public_dimensions")(prep, biters = 1000, clustered = FALSE)
  produced <- rid_ns("chained_public_pipeline_control")(dims, resolved, gmm_weighting = "identity")
  res <- rid_ns("run_chained_ipw_pipeline")(data = prep$data_internal, idname = "id", tname = "time", yname = "y", gname = "G",
    control = produced$control, xformla = prep$xformla, target_cohorts = prep$target_cohorts,
    control_group = "notyettreated", base_period = "varying", weighting = "identity", run_numerical_transfer = FALSE,
    aggregations = character(), cband = FALSE, bstrap = FALSE, caller_rng = TRUE, cp04 = prep$cp04,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, gmm_covariance = "tikhonov")
  b <- res$nuisance_registry$cp04
  ref <- rid_ns("chained_backend_reference")(b)
  source <- unname(rid_ns("cp04_payload")(prep$cp04$reference_source)["entries"])
  testthat::expect_equal(ref$accounting$source_entries, source)
  testthat::expect_gt(ref$accounting$owned_entries, source)
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  m <- max(b$layout$by_orig$count)
  testthat::expect_gt(m, 1)
  testthat::expect_gt(b$layout$N_layout, st$N)
  e0 <- rid_ns("augmented_if_envelope")(b, 0L)
  e1 <- rid_ns("augmented_if_envelope")(b, 1L)
  G <- length(st$cohorts); p_or <- b$outcome_recipe$p; p_trt <- b$treatment_recipe$p
  n_alpha <- length(unique(vapply(b$links, `[[`, integer(1), "or_block_id")))
  width <- 3 * G + 24 * n_alpha + 40 * p_or + length(e1$uids) +
    length(e1$tids) * (6 * p_trt + 8) + 48 * length(b$links) + 8
  # This fails when B counts only original rows while phi_l retains all copies.
  testthat::expect_equal(e1$output - e0$output, m * (width + 16 + frame$id_byte_metadata$max))
  rows <- seq_len(min(st$N, b$control$max_chunk_rows))
  ev <- rid_ns("evaluate_augmented_if_chunk")(b, rows, seq_along(b$links), scale = "frame")
  pre <- rid_ns("augmented_if_envelope")(b, length(rows))
  testthat::expect_equal(ev$diagnostics$output_entries_bound, pre$output)
  testthat::expect_lte(unname(rid_ns("cp04_payload")(ev)["entries"]), pre$output)
  # The layout count enters the fitting reservation, never empirical normalizers.
  n <- b$layout$N_layout; nc <- length(b$cells); nb <- b$accounting$n_or_blocks; nt <- length(b$treatments)
  floor <- 256 * n * ((nc + nb + 1) * (p_or + 2)^2 * (1 + ceiling(log2(n + 1))) +
    nt * (p_trt + 2)^2 * b$glm_control$maxit + length(b$links) * (p_or + p_trt + G + 4))
  testthat::expect_gte(b$accounting$work_ops_bound, floor)
  testthat::expect_equal(ev$normalization_count, length(unique(d$id)))
})
