# Public integration tests for D48 CP07: chained doubly robust (DR) estimation
# via public att_gt(..., est_method = "dr", mode = "chained") and aggte().

drp_has_installed_cdid <- function() {
  if (exists("or_public_helpers", inherits = TRUE)) return(FALSE)
  isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE)
}

drp_ns <- function(name) {
  if (exists("or_public_helpers", inherits = TRUE)) {
    env <- get("or_public_helpers", inherits = TRUE)
    stopifnot(is.environment(env), exists(name, envir = env, inherits = FALSE))
    fn <- get(name, envir = env, inherits = FALSE)
    stopifnot(is.function(fn))
    return(fn)
  }
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}

sys.source(testthat::test_path("reference", "cp06_pooled_public_or_reference.R"), envir = environment())
sys.source(testthat::test_path("reference", "selection_public_reference.R"), envir = environment())
sys.source(testthat::test_path("reference", "cp07_public_dr_reference.R"), envir = environment())

drp_raw <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_baseline", "raw.csv"))
drp_raw$cohort <- as.numeric(replace(drp_raw$cohort, drp_raw$cohort == "NT", "Inf"))
drp_trans <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_corollary1", "transport.csv"))
drp_trans$cohort <- as.numeric(replace(drp_trans$cohort, drp_trans$cohort == "NT", "Inf"))

drp_public_fixture <- function(route = "none",
                               selection = if (identical(route, "none")) "none" else "conditional_mtar",
                               selection_model = "pair_logit",
                               selection_pooling = "separate",
                               xformla = ~x,
                               selection_xformla = ~x,
                               est_method = "dr") {
  d <- drp_raw
  N <- nrow(d)
  dates <- 1:4
  Y <- as.matrix(d[paste0("Y", dates)])
  ids <- paste0("unit-", d$id)
  g_mapped <- ifelse(is.infinite(d$cohort), 0, d$cohort)
  data <- data.frame(
    id = rep(ids, each = length(dates)),
    time = rep(dates, N),
    y = as.vector(t(Y)),
    g = rep(g_mapped, each = length(dates)),
    stringsAsFactors = FALSE
  )
  # Retain one all-NA unit in panel while making the remaining all-NA units frame-only
  zero <- which(rowSums(!is.na(Y)) == 0L)
  data <- data[!data$id %in% ids[zero[-1L]], , drop = FALSE]
  frame <- data.frame(
    id = ids,
    g = g_mapped,
    x = d$x,
    w = d$w,
    cl = paste0("cl-", (d$id - 1L) %/% 3L),
    stringsAsFactors = FALSE
  )
  frame[["(Intercept)"]] <- d$x + 1.5
  if (route %in% c("B", "A")) {
    tb <- drp_trans[drp_trans$route == route, ]
    frame$rho <- tb$rho[match(paste(d$cohort, d$x), paste(tb$cohort, tb$x))]
  }
  args <- list(
    yname = "y", tname = "time", idname = "id", gname = "g",
    data = data, frame = frame, xformla = xformla, weightsname = "w",
    est_method = est_method, control_group = "nevertreated",
    bstrap = FALSE, cband = FALSE, gmm_weighting = "identity"
  )
  if (!identical(route, "none")) {
    args$selection <- selection
    args$selection_model <- selection_model
    args$selection_pooling <- selection_pooling
    args$selection_xformla <- selection_xformla
    args$selection_universe <- route
    if (route %in% c("B", "A")) {
      args$selection_rho <- "rho"
      args$selection_justification <- "frozen population-law normalized transport ratio"
    }
  }
  list(args = args, frame = frame, data = data, raw = d)
}

testthat::test_that("pure-table DR consumers match cp07_public_dr_reference across baseline MTAR and all 12 conditional MTAR modes for NT and NYT", {
  fit_b_tab <- drp_ns("fit_dr_baseline_table")
  fit_s_tab <- drp_ns("fit_dr_selection_table")
  mk_recipe <- drp_ns("or_pooled_outcome_recipe")
  Y <- as.matrix(drp_raw[paste0("Y", 1:4)])
  obs <- as.matrix(drp_raw[paste0("R", 1:4)])
  Z <- cbind(`(Intercept)` = 1, x = drp_raw$x)
  attr(Z, "intercept") <- TRUE; attr(Z, "assign") <- c(0L, 1L)
  V <- Z
  recipe_or <- mk_recipe(colnames(V), TRUE, "x")

  for (cg in c("nevertreated", "notyettreated")) {
    cat_info <- cp06_pooled_public_catalog(drp_raw, dates = 1:4, control_group = cg)
    links_df <- cat_info$links_df
    links <- lapply(seq_len(nrow(links_df)), function(i) {
      list(
        link_id = paste0(links_df$role[i], "_", i),
        s = links_df$s[i], t = links_df$t[i], g = links_df$g[i],
        K = as.numeric(strsplit(links_df$controls[i], "|", fixed = TRUE)[[1L]])
      )
    })

    # Baseline MTAR
    tab_b <- fit_b_tab(
      X_or = V, V = V, Y = Y, obs = obs, g = drp_raw$cohort, w = drp_raw$w,
      calendar = 1:4, links = links, outcome_recipe = recipe_or
    )
    ref_b <- cp07_public_dr_reference(drp_raw, route = "none", control_group = cg)
    testthat::expect_lt(max(abs(tab_b$estimates - ref_b$mu)), 1e-11)
    testthat::expect_lt(max(abs(tab_b$phi - ref_b$phi_mu)), 3e-7)

    # All 12 conditional MTAR modes
    for (sm in c("pair_logit", "independent_period_logit")) {
      for (sp in c("separate", "within_cohort", "across_cohorts")) {
        routes <- if (identical(sm, "independent_period_logit")) "F" else c("F", "B", "A")
        for (route in routes) {
          U <- if (route == "F") rep(TRUE, nrow(drp_raw)) else if (route == "B") (drp_raw$B == 1) else (drp_raw$A == 1)
          tb <- drp_trans[drp_trans$route == route, ]
          raw_r <- drp_raw
          raw_r$rho <- tb$rho[match(paste(drp_raw$cohort, drp_raw$x), paste(tb$cohort, tb$x))]

          tab_s <- fit_s_tab(
            Z = Z, X_or = V, V = V, Y = Y, obs = obs, g = raw_r$cohort, w = raw_r$w,
            U = U, rho = raw_r$rho, universe = route, calendar = 1:4, links = links,
            outcome_recipe = recipe_or, selection_model = sm, selection_pooling = sp
          )
          ref_s <- cp07_public_dr_reference(
            raw_r, route = route, selection_model = sm, selection_pooling = sp, control_group = cg
          )
          testthat::expect_lt(max(abs(tab_s$estimates - ref_s$mu)), 1e-11)
          testthat::expect_lt(max(abs(tab_s$phi - ref_s$phi_mu)), 3e-7)
        }
      }
    }
  }
})

testthat::test_that("public att_gt default est_method = 'dr' and explicit 'dr' match cp07_public_dr_reference across baseline MTAR and all 12 conditional MTAR modes", {
  testthat::skip_if_not(drp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- drp_ns("att_gt")
  agg_fn <- drp_ns("aggte")

  # 1. Baseline MTAR: default est_method vs explicit est_method = "dr"
  fx_none <- drp_public_fixture(route = "none")
  args_default <- fx_none$args
  args_default$est_method <- NULL
  fit_def <- do.call(att_fn, args_default)
  fit_dr <- do.call(att_fn, fx_none$args)

  testthat::expect_identical(fit_def$att, fit_dr$att)
  testthat::expect_identical(fit_def$se, fit_dr$se)
  testthat::expect_identical(fit_def$V_analytical, fit_dr$V_analytical)
  testthat::expect_identical(fit_def$inffunc, fit_dr$inffunc)
  testthat::expect_identical(fit_dr$DIDparams$est_method, "dr")
  testthat::expect_identical(fit_dr$chained$estimator, "chained doubly robust GMM")
  testthat::expect_identical(fit_dr$chained$selection$canonical_assumption, "none")
  testthat::expect_identical(fit_dr$chained$selection$outcome_model, "shared one-period outcome-change control-mean regression")
  testthat::expect_identical(fit_dr$chained$selection$treatment_model, "risk-set treatment propensity logit")
  diag <- fit_dr$chained$nuisance_diagnostics
  testthat::expect_false(diag$accounting$full_F_by_link_IF_retained)
  testthat::expect_equal(diag$accounting$n_regression_fits, 1L)
  testthat::expect_equal(diag$accounting$n_regression_blocks, length(unique(diag$fitted_catalog$regression_block)))
  testthat::expect_equal(diag$accounting$n_treatment_fits, length(unique(diag$fitted_catalog$treatment_model)))
  testthat::expect_length(diag$records$treatment, diag$accounting$n_treatment_fits)
  bound <- 4 * diag$accounting$owned_entries + 4096
  resource <- fit_dr$chained$resource_accounting
  testthat::expect_equal(resource$final_diagnostic_copy_entries, 5 * bound, tolerance = 0)
  testthat::expect_lte(max(resource$simultaneous_entries, resource$initial_envelope_entries) +
    resource$final_diagnostic_copy_entries, fit_dr$chained$public$numeric_workspace_entries)

  ref_none <- cp07_public_dr_reference(drp_raw, route = "none", control_group = "nevertreated")
  testthat::expect_equal(fit_dr$att, ref_none$att, tolerance = 1e-10)
  testthat::expect_equal(unname(fit_dr$inffunc), unname(ref_none$inffunc), tolerance = 1e-7)
  testthat::expect_equal(unname(fit_dr$V_analytical), unname(ref_none$V_analytical), tolerance = 1e-7)
  testthat::expect_equal(fit_dr$se, sqrt(diag(ref_none$V_analytical) / fit_dr$n), tolerance = 1e-8)
  testthat::expect_identical(rownames(fit_dr$inffunc), fx_none$frame$id)
  testthat::expect_equal(fit_dr$n, nrow(fx_none$frame))

  for (tp in c("simple", "group", "dynamic", "calendar")) {
    ag <- agg_fn(fit_dr, type = tp)
    testthat::expect_equal(ag$overall.att, unname(ref_none$agg_overall[tp]), tolerance = 1e-10)
    if_ag <- as.numeric(switch(tp,
      simple = ag$inf.function$simple.att,
      group = ag$inf.function$selective.inf.func,
      calendar = ag$inf.function$calendar.inf.func,
      dynamic = ag$inf.function$dynamic.inf.func
    ))
    testthat::expect_equal(if_ag, as.numeric(ref_none$phi_agg_overall[, tp]), tolerance = 1e-7)
    se_ref <- sqrt(sum((if_ag - mean(if_ag))^2) / (fit_dr$n^2))
    testthat::expect_equal(ag$overall.se, se_ref, tolerance = 1e-8)
  }

  # 2. All 12 conditional MTAR modes (pair_logit F/B/A x 3 poolings; independent_period_logit F x 3 poolings)
  for (sm in c("pair_logit", "independent_period_logit")) {
    for (sp in c("separate", "within_cohort", "across_cohorts")) {
      routes <- if (identical(sm, "independent_period_logit")) "F" else c("F", "B", "A")
      for (route in routes) {
        fx <- drp_public_fixture(
          route = route, selection = "conditional_mtar",
          selection_model = sm, selection_pooling = sp
        )
        fit <- do.call(att_fn, fx$args)
        raw_r <- drp_raw
        if (route %in% c("B", "A")) raw_r$rho <- fx$frame$rho
        ref <- cp07_public_dr_reference(
          raw_r, route = route, selection_model = sm, selection_pooling = sp,
          control_group = "nevertreated"
        )

        testthat::expect_equal(fit$att, ref$att, tolerance = 1e-10)
        testthat::expect_equal(unname(fit$inffunc), unname(ref$inffunc), tolerance = 1e-7)
        testthat::expect_equal(unname(fit$V_analytical), unname(ref$V_analytical), tolerance = 1e-7)
        testthat::expect_equal(fit$se, sqrt(diag(ref$V_analytical) / fit$n), tolerance = 1e-8)

        sel_rec <- fit$chained$selection
        testthat::expect_identical(sel_rec$canonical_assumption, "conditional_mtar")
        testthat::expect_identical(sel_rec$selection_model, sm)
        testthat::expect_identical(sel_rec$selection_pooling, sp)
        testthat::expect_identical(sel_rec$fitting_universe, route)
        testthat::expect_identical(sel_rec$outcome_model, "shared one-period outcome-change control-mean regression")
        testthat::expect_identical(sel_rec$treatment_model, "risk-set treatment propensity logit")
        testthat::expect_equal(sel_rec$n_active_cells, ref$n_active_cells)
        testthat::expect_equal(sel_rec$n_selection_blocks, ref$n_selection_blocks)
        testthat::expect_equal(unname(sel_rec$empirical_repeat_rates), unname(ref$a_U), tolerance = 1e-12)
        diag <- fit$chained$nuisance_diagnostics
        testthat::expect_false(diag$accounting$full_F_by_link_IF_retained)
        testthat::expect_equal(diag$accounting$n_regression_fits, 1L)
        testthat::expect_equal(diag$accounting$n_regression_blocks, length(unique(diag$fitted_catalog$regression_block)))
        testthat::expect_equal(diag$accounting$n_treatment_fits, length(unique(diag$fitted_catalog$treatment_model)))
        testthat::expect_length(diag$records$treatment, diag$accounting$n_treatment_fits)
        bound <- 4 * diag$accounting$owned_entries + 4096
        resource <- fit$chained$resource_accounting
        testthat::expect_equal(resource$final_diagnostic_copy_entries, 5 * bound, tolerance = 0)
        testthat::expect_lte(max(resource$simultaneous_entries, resource$initial_envelope_entries) +
          resource$final_diagnostic_copy_entries, fit$chained$public$numeric_workspace_entries)

        for (tp in c("simple", "group", "dynamic", "calendar")) {
          ag <- agg_fn(fit, type = tp)
          testthat::expect_equal(ag$overall.att, unname(ref$agg_overall[tp]), tolerance = 1e-10)
          if_ag <- as.numeric(switch(tp,
            simple = ag$inf.function$simple.att,
            group = ag$inf.function$selective.inf.func,
            calendar = ag$inf.function$calendar.inf.func,
            dynamic = ag$inf.function$dynamic.inf.func
          ))
          testthat::expect_equal(if_ag, as.numeric(ref$phi_agg_overall[, tp]), tolerance = 1e-7)
        }
      }
    }
  }
})

testthat::test_that("canonical selection = 'conditional_mtar' has exact parity with legacy alias 'corollary1' under DR across F/B/A", {
  testthat::skip_if_not(drp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- drp_ns("att_gt")
  agg_fn <- drp_ns("aggte")

  for (route in c("F", "B", "A")) {
    fx_c <- drp_public_fixture(
      route = route, selection = "conditional_mtar",
      selection_model = "pair_logit", selection_pooling = "separate"
    )
    fx_a <- drp_public_fixture(
      route = route, selection = "corollary1",
      selection_model = "pair_logit", selection_pooling = "separate"
    )
    fit_c <- do.call(att_fn, fx_c$args)
    fit_a <- do.call(att_fn, fx_a$args)

    testthat::expect_identical(fit_c$att, fit_a$att)
    testthat::expect_identical(fit_c$se, fit_a$se)
    testthat::expect_identical(fit_c$V_analytical, fit_a$V_analytical)
    testthat::expect_identical(fit_c$inffunc, fit_a$inffunc)
    testthat::expect_identical(fit_c$chained$selection$canonical_assumption, "conditional_mtar")
    testthat::expect_identical(fit_a$chained$selection$canonical_assumption, "conditional_mtar")
    testthat::expect_identical(fit_c$chained$selection$requested_selection, "conditional_mtar")
    testthat::expect_identical(fit_a$chained$selection$requested_selection, "corollary1")

    for (tp in c("simple", "group", "dynamic", "calendar")) {
      ag_c <- agg_fn(fit_c, type = tp)
      ag_a <- agg_fn(fit_a, type = tp)
      testthat::expect_identical(ag_c$overall.att, ag_a$overall.att)
      testthat::expect_identical(ag_c$overall.se, ag_a$overall.se)
      testthat::expect_identical(ag_c$att.egt, ag_a$att.egt)
      testthat::expect_identical(ag_c$se.egt, ag_a$se.egt)
    }
  }
})

testthat::test_that("intercept-only DR and IPW have exact point, complete IF, covariance, and aggte identity in the same correction mode", {
  testthat::skip_if_not(drp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- drp_ns("att_gt")
  agg_fn <- drp_ns("aggte")

  modes <- list(
    list(route = "none", sm = "pair_logit", sp = "separate"),
    list(route = "F", sm = "pair_logit", sp = "separate"),
    list(route = "B", sm = "pair_logit", sp = "within_cohort"),
    list(route = "A", sm = "pair_logit", sp = "across_cohorts"),
    list(route = "F", sm = "independent_period_logit", sp = "within_cohort")
  )
  for (cfg in modes) {
    fx_dr <- drp_public_fixture(
      route = cfg$route, selection_model = cfg$sm, selection_pooling = cfg$sp,
      xformla = ~1, selection_xformla = ~x, est_method = "dr"
    )
    fx_ipw <- drp_public_fixture(
      route = cfg$route, selection_model = cfg$sm, selection_pooling = cfg$sp,
      xformla = ~1, selection_xformla = ~x, est_method = "ipw"
    )
    fit_dr <- do.call(att_fn, fx_dr$args)
    fit_ipw <- do.call(att_fn, fx_ipw$args)

    testthat::expect_equal(fit_dr$att, fit_ipw$att, tolerance = 1e-12)
    testthat::expect_equal(unname(fit_dr$inffunc), unname(fit_ipw$inffunc), tolerance = 1e-11)
    testthat::expect_equal(unname(fit_dr$V_analytical), unname(fit_ipw$V_analytical), tolerance = 1e-11)
    testthat::expect_equal(fit_dr$se, fit_ipw$se, tolerance = 1e-11)

    for (tp in c("simple", "group", "dynamic", "calendar")) {
      ag_dr <- agg_fn(fit_dr, type = tp)
      ag_ipw <- agg_fn(fit_ipw, type = tp)
      testthat::expect_equal(ag_dr$overall.att, ag_ipw$overall.att, tolerance = 1e-12)
      testthat::expect_equal(ag_dr$overall.se, ag_ipw$overall.se, tolerance = 1e-11)
      testthat::expect_equal(ag_dr$att.egt, ag_ipw$att.egt, tolerance = 1e-12)
      testthat::expect_equal(ag_dr$se.egt, ag_ipw$se.egt, tolerance = 1e-11)
    }
  }
})

testthat::test_that("public DR supports NYT, PT-all, universal base_period, anticipation, formula variants, basis invariance, permutations, chunk invariance, clustering, and bootstrap", {
  testthat::skip_if_not(drp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- drp_ns("att_gt")
  agg_fn <- drp_ns("aggte")

  # 1. NYT, PT-post/PT-all, varying/universal base_period, and anticipation
  for (cfg in list(
    list(route = "none", sm = "pair_logit", sp = "separate"),
    list(route = "F", sm = "independent_period_logit", sp = "within_cohort"),
    list(route = "B", sm = "pair_logit", sp = "across_cohorts")
  )) {
    fx <- drp_public_fixture(
      route = cfg$route, selection_model = cfg$sm, selection_pooling = cfg$sp
    )
    raw_r <- drp_raw
    if (cfg$route %in% c("B", "A")) raw_r$rho <- fx$frame$rho

    for (bp in c("varying", "universal")) {
      for (pt in c("post", "all")) {
        ps <- if (pt == "all") 2 else NULL
        fit <- do.call(att_fn, utils::modifyList(fx$args, list(
          control_group = "notyettreated", base_period = bp,
          pt_assumption = pt, pt_start = ps
        )))
        ref <- cp07_public_dr_reference(
          raw_r, route = cfg$route, selection_model = cfg$sm, selection_pooling = cfg$sp,
          control_group = "notyettreated", base_period = bp,
          pt_assumption = pt, pt_start = ps
        )
        testthat::expect_equal(fit$att, ref$att, tolerance = 1e-10)
        avail <- !is.na(ref$att)
        testthat::expect_equal(
          unname(fit$inffunc[, avail, drop = FALSE]),
          unname(ref$inffunc[, avail, drop = FALSE]),
          tolerance = 1e-7
        )
        testthat::expect_equal(
          unname(fit$V_analytical[avail, avail, drop = FALSE]),
          unname(ref$V_analytical[avail, avail, drop = FALSE]),
          tolerance = 1e-7
        )
      }
    }

    # Anticipation = 1
    fit_ant <- do.call(att_fn, utils::modifyList(fx$args, list(
      control_group = "notyettreated", anticipation = 1L
    )))
    ref_ant <- cp07_public_dr_reference(
      raw_r, route = cfg$route, selection_model = cfg$sm, selection_pooling = cfg$sp,
      control_group = "notyettreated", anticipation = 1L
    )
    testthat::expect_equal(fit_ant$att, ref_ant$att, tolerance = 1e-10)
    avail_ant <- !is.na(ref_ant$att)
    testthat::expect_equal(
      unname(fit_ant$inffunc[, avail_ant, drop = FALSE]),
      unname(ref_ant$inffunc[, avail_ant, drop = FALSE]),
      tolerance = 1e-7
    )
  }

  # 2. Formula variants (~0 + x, ~0) and separate xformla vs selection_xformla
  fx_f <- drp_public_fixture(
    route = "F", selection_model = "independent_period_logit", selection_pooling = "across_cohorts"
  )
  # An empty common outcome/treatment design is an existing public refusal;
  # zero-slope OR is the supported intercept-only ~1 case checked above.
  testthat::expect_error(do.call(att_fn, utils::modifyList(fx_f$args, list(
    xformla = ~0, selection_xformla = ~0 + x
  ))), class = "cdid_chained_request_covariate_error")
  for (v in list(
    list(xf = ~0 + x, sf = ~x, or_icpt = FALSE, X_sl = cbind(x = drp_raw$x), V_tr = cbind(x = drp_raw$x), sel_icpt = TRUE, sel_zero = FALSE, Q_sl = cbind(x = drp_raw$x)),
    list(xf = ~x, sf = ~0, or_icpt = TRUE, X_sl = cbind(x = drp_raw$x), V_tr = cbind(`(Intercept)` = 1, x = drp_raw$x), sel_icpt = FALSE, sel_zero = TRUE, Q_sl = NULL)
  )) {
    fit_v <- do.call(att_fn, utils::modifyList(fx_f$args, list(
      xformla = v$xf, selection_xformla = v$sf
    )))
    ref_v <- cp07_public_dr_reference(
      drp_raw, route = "F", selection_model = "independent_period_logit",
      selection_pooling = "across_cohorts",
      sel_intercept = v$sel_icpt, sel_zero_col = v$sel_zero, Q_slope = v$Q_sl,
      or_intercept = v$or_icpt, X_slope = v$X_sl, V_treat = v$V_tr
    )
    testthat::expect_equal(fit_v$att, ref_v$att, tolerance = 1e-10)
    testthat::expect_equal(unname(fit_v$inffunc), unname(ref_v$inffunc), tolerance = 1e-7)
  }

  # 3. Equivalent shifted/scaled basis with intercept vs non-equivalent no-intercept shift
  fx_b <- drp_public_fixture(route = "none")
  fx_b$frame$x_shift <- 2.5 * fx_b$frame$x - 1.25
  fit_x <- do.call(att_fn, fx_b$args)
  fit_shift_icpt <- do.call(att_fn, utils::modifyList(fx_b$args, list(
    frame = fx_b$frame, xformla = ~x_shift
  )))
  testthat::expect_equal(fit_shift_icpt$att, fit_x$att, tolerance = 1e-10)
  testthat::expect_equal(unname(fit_shift_icpt$inffunc), unname(fit_x$inffunc), tolerance = 1e-8)

  fit_noicpt_x <- do.call(att_fn, utils::modifyList(fx_b$args, list(
    frame = fx_b$frame, xformla = ~0 + x
  )))
  fit_noicpt_shift <- do.call(att_fn, utils::modifyList(fx_b$args, list(
    frame = fx_b$frame, xformla = ~0 + x_shift
  )))
  testthat::expect_gt(max(abs(fit_noicpt_x$att - fit_noicpt_shift$att)), 1e-3)

  # 4. Row/ID permutations and chunk_rows invariance
  set.seed(4807L)
  perm_frame <- fx_f$frame[sample.int(nrow(fx_f$frame)), , drop = FALSE]
  perm_data <- fx_f$data[sample.int(nrow(fx_f$data)), , drop = FALSE]
  fit_orig <- do.call(att_fn, fx_f$args)
  fit_perm <- do.call(att_fn, utils::modifyList(fx_f$args, list(
    data = perm_data, frame = perm_frame, control = list(chunk_rows = 7L)
  )))
  fit_chunk23 <- do.call(att_fn, utils::modifyList(fx_f$args, list(control = list(chunk_rows = 23L))))

  testthat::expect_equal(fit_perm$att, fit_orig$att, tolerance = 1e-12)
  testthat::expect_equal(fit_perm$se, fit_orig$se, tolerance = 1e-12)
  testthat::expect_equal(unname(fit_perm$V_analytical), unname(fit_orig$V_analytical), tolerance = 1e-12)
  testthat::expect_equal(
    unname(fit_perm$inffunc[fx_f$frame$id, , drop = FALSE]),
    unname(fit_orig$inffunc),
    tolerance = 1e-12
  )
  testthat::expect_equal(fit_chunk23$att, fit_orig$att, tolerance = 1e-13)
  testthat::expect_equal(unname(fit_chunk23$inffunc), unname(fit_orig$inffunc), tolerance = 1e-13)

  # 5. Clustering (inference-only) and multiplier bootstrap
  fit_cl <- do.call(att_fn, utils::modifyList(fx_f$args, list(clustervars = "cl")))
  testthat::expect_identical(fit_cl$att, fit_orig$att)
  testthat::expect_identical(fit_cl$inffunc, fit_orig$inffunc)
  S_cl <- rowsum(fit_orig$inffunc, fx_f$frame$cl, reorder = FALSE)
  S_cl <- sweep(S_cl, 2L, colMeans(S_cl), "-")
  V_cl <- crossprod(S_cl) / nrow(fx_f$frame)
  testthat::expect_equal(unname(fit_cl$V_analytical), unname(V_cl), tolerance = 1e-12)
  testthat::expect_equal(fit_cl$se, unname(sqrt(diag(V_cl) / nrow(fx_f$frame))), tolerance = 1e-12)

  set.seed(20261021L)
  fit_bs <- do.call(att_fn, utils::modifyList(fx_f$args, list(
    clustervars = "cl", bstrap = TRUE, cband = TRUE, biters = 64L
  )))
  testthat::expect_identical(fit_bs$att, fit_orig$att)
  testthat::expect_true(all(is.finite(fit_bs$se) & fit_bs$se > 0))
  testthat::expect_true(is.finite(fit_bs$c) && fit_bs$c > 0)
  ag_bs <- agg_fn(fit_bs, type = "dynamic")
  testthat::expect_true(is.finite(ag_bs$overall.att) && is.finite(ag_bs$overall.se) && ag_bs$overall.se > 0)

  # 6. Result reporting adapters: tidy / glance / ggdid
  td <- generics::tidy(fit_orig)
  gl <- generics::glance(fit_orig)
  testthat::expect_true("estimator" %in% names(td))
  testthat::expect_true(all(td$estimator == "chained doubly robust GMM"))
  testthat::expect_identical(gl$estimator, "chained doubly robust GMM")
  td_ag <- generics::tidy(ag_bs)
  gl_ag <- generics::glance(ag_bs)
  testthat::expect_true("estimator" %in% names(td_ag))
  testthat::expect_true(all(td_ag$estimator == "chained doubly robust GMM"))
  testthat::expect_identical(gl_ag$estimator, "chained doubly robust GMM")
})

testthat::test_that("public DR enforces pre-allocation resource bounds and rank-deficiency refusals truthfully", {
  testthat::skip_if_not(drp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- drp_ns("att_gt")

  fx <- drp_public_fixture(route = "F", selection_model = "pair_logit", selection_pooling = "separate")

  # Tight retained entries bound refuses before allocation
  ctrl_ret <- list(chunk_rows = 13L, numeric_workspace_bytes = 8 * 200L)
  testthat::expect_error(
    do.call(att_fn, utils::modifyList(fx$args, list(control = ctrl_ret))),
    class = "cdid_ipw_owner_bound_error"
  )

  # Tight work_ops bound refuses before allocation
  ctrl_work <- list(chunk_rows = 13L, counted_work_limit = 500L)
  testthat::expect_error(
    do.call(att_fn, utils::modifyList(fx$args, list(control = ctrl_work))),
    class = "cdid_ipw_owner_bound_error"
  )

  # Collinear covariate design refuses on rank-deficient treatment propensity logit
  fx_col <- fx
  fx_col$frame$x_dup <- 2 * fx_col$frame$x + 1
  testthat::expect_error(
    do.call(att_fn, utils::modifyList(fx_col$args, list(
      frame = fx_col$frame, xformla = ~x + x_dup
    ))),
    class = "cdid_error"
  )
})
