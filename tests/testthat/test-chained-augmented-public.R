# Public integration tests for the D52 augmented observation correction
# (`selection_est_method = "dr"` with `est_method = "dr"`) through
# cdid::att_gt()/aggte().
#
# Section B compares the installed public path with the frozen independent
# public oracle (fixtures/aug_public_v1: byte copies of the accepted
# public_reference transport, eight configurations on the accepted 144-unit
# fixture, identity GMM, analytic inference). Its criteria are the fixed public
# criteria of tools/verify_installed_imputation_public.R: cell, aggregate and
# primitive link points max-abs < 1e-11; complete original-F IF rows max-abs
# < 3e-7; V_analytical all.equal 1e-7; se and aggregate se all.equal 1e-8.
# Cells map by (group, t) and IF rows by original F identity. The other
# sections are gate, wiring, metadata, invariance and refusal checks; no
# expected value anywhere in this file is produced by the backend under test.

augp_has_installed_cdid <- function() isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE)

augp_val <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  get(name, envir = ns, inherits = FALSE)
}

augp_dir <- testthat::test_path("fixtures", "aug_public_v1")
# The accepted independent public helpers (unchanged) in their own environment.
augp_ref <- new.env(parent = as.environment("package:stats"))
sys.source(testthat::test_path("reference", "imputation_public_reference.R"), envir = augp_ref, keep.source = FALSE)
augp_ref$imputation_public_helpers(helper_dir = testthat::test_path("reference"))
augp_raw <- augp_ref$imputation_public_raw(utils::read.csv(file.path(augp_dir, "raw.csv")))
# Deterministic selection-only covariate of the oracle (validate_reference.R).
augp_raw$V <- ((5 * augp_raw$id) %% 11 - 5) / 5
augp_mx <- function(x) max(c(0, abs(x)))
augp_key <- function(g, t) paste(as.numeric(g), as.numeric(t), sep = "|")

# The eight oracle configurations (public_reference/validate_reference.R).
# "date_intercept_only_NT" is intercept-only in its SELECTION formula only: the
# oracle's `args$X_slope <- NULL` removes an absent list element, so its outcome
# and treatment designs keep x (saved metadata: treatment design dimension 2).
augp_configs <- list(
  pair_sep_jsep_NT_NULL = list(model = "pair_logit", qpool = "separate", jpool = "separate",
    cg = "nevertreated", q = "NULL"),
  pair_within_jacross_NYT_xV = list(model = "pair_logit", qpool = "within_cohort", jpool = "across_cohorts",
    cg = "notyettreated", q = "xV"),
  pair_across_jwithin_NYT_V = list(model = "pair_logit", qpool = "across_cohorts", jpool = "within_cohort",
    cg = "notyettreated", q = "V"),
  date_sep_jwithin_NYT_xV = list(model = "independent_period_logit", qpool = "separate", jpool = "within_cohort",
    cg = "notyettreated", q = "xV"),
  date_within_jsep_NT_V = list(model = "independent_period_logit", qpool = "within_cohort", jpool = "separate",
    cg = "nevertreated", q = "V"),
  date_across_jacross_NYT_NULL = list(model = "independent_period_logit", qpool = "across_cohorts",
    jpool = "across_cohorts", cg = "notyettreated", q = "NULL"),
  pair_irregular_all_universal_cluster_xV = list(model = "pair_logit", qpool = "across_cohorts",
    jpool = "within_cohort", cg = "notyettreated", q = "xV", irregular = TRUE),
  date_intercept_only_NT = list(model = "independent_period_logit", qpool = "within_cohort",
    jpool = "across_cohorts", cg = "nevertreated", q = "intercept_only")
)

# Public request of one configuration plus the oracle blueprint (link identities).
augp_request <- function(cfg) {
  dates <- if (isTRUE(cfg$irregular)) c(1, 2, 3, 5) else 1:4
  fx <- augp_ref$imputation_public_fixture(augp_raw, dates = dates)
  fx$frame$V <- augp_raw$V
  args <- fx$args
  args$frame <- fx$frame
  args$selection <- "conditional_mtar"
  args$selection_est_method <- "dr"
  args$est_method <- "dr"
  args$imputation_pooling <- cfg$jpool
  args$selection_model <- cfg$model
  args$selection_pooling <- cfg$qpool
  args$control_group <- cfg$cg
  # NULL inherits xformla = ~x (the oracle's explicit Q_slope = bp$Z).
  sel <- switch(cfg$q, `NULL` = NULL, xV = ~x + V, V = ~V, intercept_only = ~1)
  if (!is.null(sel)) args$selection_xformla <- sel
  bp_args <- list(raw = augp_raw, imputation_pooling = cfg$jpool, control_group = cfg$cg)
  if (isTRUE(cfg$irregular)) {
    args$pt_assumption <- "all"
    args$pt_start <- 2
    args$base_period <- "universal"
    args$clustervars <- "cl"
    bp_args <- c(bp_args, list(dates = dates, pt_assumption = "all", pt_start = 2, base_period = "universal"))
  }
  list(args = args, fx = fx, dates = dates, bp = do.call(augp_ref$imputation_public_blueprint, bp_args))
}

augp_state <- new.env(parent = emptyenv())
augp_oracle_fit <- function(name) {
  if (is.null(augp_state[[name]])) {
    rq <- augp_request(augp_configs[[name]])
    augp_state[[name]] <- list(rq = rq, fit = do.call(cdid::att_gt, rq$args),
      ref = readRDS(file.path(augp_dir, "expected", paste0(name, ".rds"))))
  }
  augp_state[[name]]
}

# Overrides are scalar/formula/list arguments only (never data frames).
augp_base_args <- function(...) utils::modifyList(augp_request(augp_configs$pair_sep_jsep_NT_NULL)$args, list(...))
augp_fit <- function(...) do.call(cdid::att_gt, augp_base_args(...))

augp_refuse <- function(cls, ..., entry = cdid::att_gt) {
  base <- list(yname = "y", tname = "time", idname = "id", gname = "g",
    data = quote(stop("data must not be evaluated")))
  e <- testthat::expect_error(do.call(entry, utils::modifyList(base, list(...))), class = cls[[1L]])
  for (k in cls[-1L]) testthat::expect_s3_class(e, k)
  testthat::expect_no_match(conditionMessage(e), "must not be evaluated|was evaluated")
  invisible(e)
}

# ------------------------------------------------------------------------------
# 0. Frozen oracle inputs (no cdid needed)
# ------------------------------------------------------------------------------

testthat::test_that("the frozen public oracle files are intact and the request mapping rebuilds its fixture", {
  man <- utils::read.csv(file.path(augp_dir, "manifest.csv"), colClasses = "character")
  testthat::expect_identical(nrow(man), 10L)
  paths <- file.path(augp_dir, man$file)
  testthat::expect_true(all(file.exists(paths)))
  testthat::expect_identical(as.character(file.size(paths)), man$bytes)
  testthat::expect_identical(unname(tools::md5sum(paths)), man$md5)
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE))
    testthat::expect_identical(unname(tools::sha256sum(paths)), man$sha256)
  # The saved oracle fixture equals the one rebuilt here from raw.csv and the V rule.
  pf <- readRDS(file.path(augp_dir, "public_fixture.rds"))
  rq <- augp_request(augp_configs$pair_sep_jsep_NT_NULL)
  testthat::expect_equal(rq$fx$frame, pf$frame)
  testthat::expect_equal(rq$fx$data, pf$data)
  testthat::expect_identical(nrow(rq$fx$frame), 144L)
  testthat::expect_identical(qr(cbind(1, augp_raw$x, augp_raw$V))$rank, 3L)
  for (name in names(augp_configs)) {
    cfg <- augp_configs[[name]]
    ref <- readRDS(file.path(augp_dir, "expected", paste0(name, ".rds")))
    rq <- augp_request(cfg)
    meta <- ref$metadata
    testthat::expect_identical(c(meta$selection_est_method, meta$est_method, meta$selection),
      c("dr", "dr", "conditional_mtar"))
    testthat::expect_identical(c(meta$selection_model, meta$selection_pooling, meta$imputation_pooling),
      c(cfg$model, cfg$qpool, cfg$jpool))
    testthat::expect_identical(meta$selection_recipe, cfg$q)
    testthat::expect_identical(meta$configuration$cg, cfg$cg)
    testthat::expect_identical(meta$treatment_design_dimension, 2L)
    # Blueprint identities used below: primitive links, cells and IF row identities.
    testthat::expect_identical(length(ref$results$dr$mu), nrow(rq$bp$links))
    testthat::expect_identical(rownames(ref$results$dr$inffunc), rq$fx$frame$id)
    testthat::expect_identical(colnames(ref$results$dr$phi_agg_overall), names(ref$results$dr$agg_overall))
    testthat::expect_false(anyDuplicated(augp_key(ref$results$dr$cells$g, ref$results$dr$cells$time)) > 0L)
  }
})

# ------------------------------------------------------------------------------
# A. Request gates (before any data promise is forced)
# ------------------------------------------------------------------------------

testthat::test_that("augmented gates refuse unsupported combinations before the data promise", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  uns <- "cdid_chained_request_unsupported_selection"
  # Treatment IPW/reg with the augmented correction: refused, nothing substituted,
  # and neither argument group is forced.
  for (m in c("ipw", "reg")) {
    e <- augp_refuse(c(uns, "cdid_chained_unavailable"), est_method = m, selection = "conditional_mtar",
      selection_est_method = "dr", imputation_pooling = quote(stop("imputation_pooling was evaluated")),
      selection_model = quote(stop("selection_model was evaluated")))
    testthat::expect_match(conditionMessage(e), "not implemented", fixed = TRUE)
    testthat::expect_match(conditionMessage(e), 'est_method = "dr"', fixed = TRUE)
  }
  augp_refuse(uns, est_method = function(...) NULL, selection = "conditional_mtar", selection_est_method = "dr")
  e <- augp_refuse(uns, selection_est_method = "dr")
  testthat::expect_match(conditionMessage(e), "conditional_mtar", fixed = TRUE)
  for (u in c("B", "A")) augp_refuse(c(uns, "cdid_cp04_unsupported_error"), selection = "conditional_mtar",
    selection_est_method = "dr", selection_universe = u, selection_justification = "declared representativeness")
  augp_refuse("cdid_chained_request_input_error", selection = "conditional_mtar", selection_est_method = "dr",
    selection_rho = "rho")
  augp_refuse(c(uns, "cdid_cp04_unsupported_error"), selection = "corollary1", selection_est_method = "dr",
    selection_model = "independent_period_logit")
  augp_refuse(c(uns, "cdid_cp04_unsupported_error"), selection = "conditional_mtar", selection_est_method = "dr",
    selection_model = "independent_period_logit", selection_justification = "not allowed for date products")
  for (bad in list("pooled", "", NA_character_, c("separate", "within_cohort"), 1L)) {
    augp_refuse(uns, selection = "conditional_mtar", selection_est_method = "dr", imputation_pooling = bad)
    augp_refuse(uns, selection = "conditional_mtar", selection_est_method = "dr", selection_pooling = bad)
    augp_refuse(uns, selection = "conditional_mtar", selection_est_method = "dr", selection_model = bad)
  }
  # The treatment-method refusal keeps its precedence.
  augp_refuse(c("cdid_chained_request_est_method", "cdid_chained_unavailable"), est_method = "stdipw",
    selection_est_method = "dr")
  dev <- augp_val("chained_att_gt_request")
  augp_refuse(c(uns, "cdid_chained_unavailable"), entry = dev, est_method = "ipw", selection = "conditional_mtar",
    selection_est_method = "dr")
})

testthat::test_that("augmented gates pass every axis combination and force both argument groups", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  base <- list(yname = "y", tname = "time", idname = "id", gname = "g",
    data = quote(stop("data must not be evaluated")))
  reach <- function(...) testthat::expect_error(do.call(cdid::att_gt, c(base, list(...))), "data must not be evaluated")
  for (sm in c("pair_logit", "independent_period_logit")) for (sp in c("separate", "within_cohort", "across_cohorts"))
    for (ip in c("separate", "within_cohort", "across_cohorts")) {
      reach(est_method = "dr", selection = "conditional_mtar", selection_est_method = "dr",
        selection_model = sm, selection_pooling = sp, imputation_pooling = ip)
    }
  reach(selection = "corollary1", selection_est_method = "dr")
  reach(selection = "conditional_mtar", selection_est_method = "dr", selection_justification = "supplied representative roster")
  # Both groups are active under "dr": an unforceable promise in either is evaluated.
  testthat::expect_error(do.call(cdid::att_gt, c(base, list(selection = "conditional_mtar", selection_est_method = "dr",
    imputation_pooling = quote(stop("imputation_pooling was evaluated"))))), "imputation_pooling was evaluated")
  testthat::expect_error(do.call(cdid::att_gt, c(base, list(selection = "conditional_mtar", selection_est_method = "dr",
    selection_pooling = quote(stop("selection_pooling was evaluated"))))), "selection_pooling was evaluated")
})

testthat::test_that("ordinary mode neither evaluates nor forwards the augmented arguments", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  d <- parity_panel()
  res <- cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d, bstrap = FALSE,
    mode = "ordinary", selection_est_method = "dr", selection = stop("selection must not be evaluated"),
    imputation_pooling = stop("imputation_pooling must not be evaluated"))
  ref <- did::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d, bstrap = FALSE)
  testthat::expect_false(inherits(res, "cdid_MP"))
  testthat::expect_identical(res$att, ref$att)
  testthat::expect_identical(res$se, ref$se)
  testthat::expect_identical(res$inffunc, ref$inffunc)
})

# ------------------------------------------------------------------------------
# B. Frozen independent public oracle (eight configurations)
# ------------------------------------------------------------------------------

testthat::test_that("installed public results match the frozen augmented public oracle", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  n_compared <- 0L
  for (name in names(augp_configs)) {
    obj <- augp_oracle_fit(name)
    fit <- obj$fit
    ref <- obj$ref$results$dr
    meta <- obj$ref$metadata
    cfg <- augp_configs[[name]]
    frame_ids <- obj$rq$fx$frame$id
    N <- length(frame_ids)
    testthat::expect_s3_class(fit, "cdid_MP")
    testthat::expect_equal(fit$n, N)
    # Cells by (group, t); IF rows by original F identity.
    pk <- augp_key(fit$group, fit$t)
    rk <- augp_key(ref$cells$g, ref$cells$time)
    j <- match(rk, pk)
    testthat::expect_identical(length(pk), length(rk))
    testthat::expect_false(anyNA(j))
    IF <- as.matrix(fit$inffunc)
    testthat::expect_identical(nrow(IF), N)
    testthat::expect_setequal(rownames(IF), frame_ids)
    ri <- match(frame_ids, rownames(IF))
    av <- !is.na(ref$att)
    testthat::expect_identical(av, !is.na(fit$att[j]))
    a <- which(av)
    testthat::expect_lt(augp_mx(fit$att[j[a]] - ref$att[a]), 1e-11)
    testthat::expect_lt(augp_mx(IF[ri, j[a], drop = FALSE] - ref$inffunc[, a, drop = FALSE]), 3e-7)
    V <- as.matrix(fit$V_analytical)
    testthat::expect_equal(unname(V[j[a], j[a], drop = FALSE]), unname(ref$V_analytical[a, a, drop = FALSE]),
      tolerance = 1e-7)
    nonref <- a[ref$cells$role[a] != "reference"]
    testthat::expect_equal(unname(fit$se[j[nonref]]), unname(ref$se[nonref]), tolerance = 1e-8)
    structural <- a[ref$cells$role[a] == "reference"]
    if (length(structural)) testthat::expect_true(all(is.na(fit$se[j[structural]]) | fit$se[j[structural]] == 0))
    n_compared <- n_compared + length(a)
    # Primitive link points by (role, cohort, s, t, control-cohort set), auxiliary rows included.
    cat_pub <- fit$chained$nuisance_diagnostics$fitted_catalog
    links <- obj$rq$bp$links
    ckey <- function(x) vapply(as.character(x), function(s) {
      z <- as.numeric(strsplit(s, "[|,]")[[1L]]); z[z == 0] <- Inf; paste(sort(z), collapse = ",")
    }, "")
    lk <- function(role, g, s, t, controls) paste(role, g, s, t, ckey(controls), sep = "|")
    jj <- match(lk(links$role, links$g, links$s, links$t, links$controls),
      lk(cat_pub$role, cat_pub$cohort, cat_pub$s, cat_pub$t, cat_pub$controls))
    testthat::expect_identical(nrow(cat_pub), nrow(links))
    testthat::expect_false(anyNA(jj))
    testthat::expect_lt(augp_mx(cat_pub$mu[jj] - ref$mu), 1e-11)
    # The fitted nuisance structure is the oracle's: cells, pooling groups,
    # observation-probability roots and their dimensions, treatment designs.
    nd <- fit$chained$nuisance_diagnostics
    testthat::expect_identical(nd$method, "augmented")
    testthat::expect_identical(nrow(nd$cells), nrow(obj$rq$bp$cells))
    testthat::expect_identical(length(nd$groups), length(obj$rq$bp$groups))
    testthat::expect_equal(nd$support$observation_probability_models, length(meta$q_design_dimensions))
    testthat::expect_equal(nd$records$selection$n_parameters, sum(meta$q_design_dimensions))
    testthat::expect_equal(nd$records$selection$n_slope_columns,
      switch(cfg$q, `NULL` = 1, xV = 2, V = 1, intercept_only = 0))
    testthat::expect_identical(fit$chained$selection$selection_formula_inherited, cfg$q == "NULL")
    testthat::expect_identical(nd$support$design_columns, c("(Intercept)", "x"))
    testthat::expect_true(all(lengths(nd$coefficients$treatment) == meta$treatment_design_dimension))
    # Aggregates with their complete share paths.
    for (tp in c("simple", "group", "dynamic", "calendar")) {
      ag <- cdid::aggte(fit, type = tp)
      testthat::expect_lt(abs(ag$overall.att - unname(ref$agg_overall[tp])), 1e-11)
      ifv <- switch(tp, simple = ag$inf.function$simple.att, group = ag$inf.function$selective.inf.func,
        dynamic = ag$inf.function$dynamic.inf.func, calendar = ag$inf.function$calendar.inf.func)
      ifv <- as.matrix(ifv)
      rows <- if (!is.null(rownames(ifv))) rownames(ifv) else rownames(IF)
      testthat::expect_identical(dim(ifv), c(N, 1L))
      testthat::expect_lt(augp_mx(ifv[match(frame_ids, rows), 1L] - ref$phi_agg_overall[, tp]), 3e-7)
      testthat::expect_equal(ag$overall.se, unname(ref$se_agg[tp]), tolerance = 1e-8)
      n_compared <- n_compared + 1L
    }
  }
  testthat::expect_gte(n_compared, 8L * 4L + 8L * 5L)
})

# ------------------------------------------------------------------------------
# C. Metadata, diagnostics and reporting
# ------------------------------------------------------------------------------

testthat::test_that("both method axes, the selection formula and every fit count are reported", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  obj <- augp_oracle_fit("pair_within_jacross_NYT_xV")
  fit <- obj$fit
  N <- nrow(obj$rq$fx$frame)
  label <- augp_val("chained_result_estimator")
  testthat::expect_identical(fit$chained$estimator, label("dr", "dr"))
  testthat::expect_match(fit$chained$estimator, "(augmented observation correction)", fixed = TRUE)
  testthat::expect_identical(fit$DIDparams$est_method, "dr")
  req <- fit$chained$request
  testthat::expect_identical(c(req$selection_est_method, req$est_method, req$imputation_pooling),
    c("dr", "dr", "across_cohorts"))
  testthat::expect_identical(c(req$selection_model, req$selection_pooling), c("pair_logit", "within_cohort"))
  testthat::expect_identical(req$inactive_arguments, character(0))
  sel <- fit$chained$selection
  testthat::expect_identical(sel$method, "augmented")
  testthat::expect_identical(c(sel$selection_est_method, sel$est_method, sel$canonical_assumption, sel$fitting_universe),
    c("dr", "dr", "conditional_mtar", "F"))
  testthat::expect_identical(c(sel$selection_model, sel$selection_pooling, sel$imputation_pooling),
    c("pair_logit", "within_cohort", "across_cohorts"))
  testthat::expect_false(sel$selection_formula_inherited)
  testthat::expect_identical(deparse(sel$formula), "~x + V")
  testthat::expect_identical(deparse(sel$outcome_formula), "~x")
  testthat::expect_length(sel$inactive_arguments, 0L)
  testthat::expect_gt(sel$models$n_selection_fits, 0)
  testthat::expect_equal(sel$models$observation_logit_models, sel$models$n_selection_blocks)
  testthat::expect_match(sel$assumptions, "common conditioning set", fixed = TRUE)
  testthat::expect_match(sel$assumptions, "sufficient, not necessary", fixed = TRUE)
  nd <- fit$chained$nuisance_diagnostics
  testthat::expect_identical(nd$method, "augmented")
  testthat::expect_identical(c(nd$selection_est_method, nd$est_method, nd$imputation_pooling),
    c("dr", "dr", "across_cohorts"))
  testthat::expect_identical(nd$support$design_source, "X_frame")
  testthat::expect_identical(nd$support$design_columns, c("(Intercept)", "x"))
  testthat::expect_identical(nd$support$selection_design_columns, c("(Intercept)", "x", "V"))
  testthat::expect_true(nd$support$observation_probabilities_fitted)
  testthat::expect_equal(nd$support$design_rows, N)
  testthat::expect_equal(nd$support$n_full_frame_units, N)
  testthat::expect_equal(nd$accounting$n_or_fits, 1)
  testthat::expect_equal(nd$accounting$observation_logit_models, nd$support$observation_probability_models)
  testthat::expect_true(nd$accounting$automatic_chunk_rows)
  testthat::expect_false(nd$accounting$full_F_by_link_IF_retained)
  testthat::expect_identical(nd$counterfactual_or$response, "reconstructed augmented response Z")
  testthat::expect_true(all(c("selection_pair_cohort", "roster_mass", "probability_weighted_mass",
    "mass_discrepancy", "imputation_total", "augmented_total") %in% names(nd$cells)))
  testthat::expect_true(all(c("catalog", "auxiliary") %in% nd$fitted_catalog$role))
  testthat::expect_true(all(nd$fitted_catalog$target_cell %in% nd$cells$cell_id))
  od <- fit$chained$outcome_design
  testthat::expect_false(od$common_A_alias_reduction)
  testthat::expect_identical(od$columns, c("(Intercept)", "x"))
  testthat::expect_identical(od$selection_columns, c("(Intercept)", "x", "V"))
  testthat::expect_equal(od$n_rows, N)
  ra <- fit$chained$resource_accounting
  testthat::expect_equal(ra$final_diagnostic_copy_entries, 5 * (4 * nd$accounting$owned_entries + 4096), tolerance = 0)
  testthat::expect_lte(max(ra$simultaneous_entries, ra$initial_envelope_entries) + ra$final_diagnostic_copy_entries,
    fit$chained$public$numeric_workspace_entries)
  testthat::expect_true(all(cdid::tidy(fit)$estimator == fit$chained$estimator))
  testthat::expect_identical(cdid::glance(fit)$estimator, fit$chained$estimator)
  txt <- paste(utils::capture.output(summary(fit)), collapse = "\n")
  testthat::expect_match(txt, "augmented correction under selection", fixed = TRUE)
  testthat::expect_match(txt, "selection_model = \"pair_logit\"", fixed = TRUE)
  testthat::expect_match(txt, "from selection_xformla", fixed = TRUE)
  testthat::expect_match(fit$chained$public$entry, "selection_est_method = \"dr\"", fixed = TRUE)
})

# ------------------------------------------------------------------------------
# D. Formulas, NULL inheritance, legacy alias and the active axes
# ------------------------------------------------------------------------------

testthat::test_that("NULL selection_xformla inherits xformla and the legacy alias is the direct-pair route", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  ref <- augp_oracle_fit("pair_sep_jsep_NT_NULL")$fit
  explicit <- augp_fit(selection_xformla = ~x)
  for (nm in c("att", "se", "V_analytical", "inffunc")) testthat::expect_identical(explicit[[nm]], ref[[nm]])
  testthat::expect_true(ref$chained$selection$selection_formula_inherited)
  testthat::expect_false(explicit$chained$selection$selection_formula_inherited)
  alias <- augp_fit(selection = "corollary1")
  for (nm in c("att", "se", "V_analytical", "inffunc")) testthat::expect_identical(alias[[nm]], ref[[nm]])
  testthat::expect_true(alias$chained$selection$legacy_alias)
  testthat::expect_identical(alias$chained$selection$requested_selection, "corollary1")
  testthat::expect_identical(alias$chained$selection$selection_model, "pair_logit")
})

testthat::test_that("distinct formulas and model axes retain saturated identities and material V contrasts", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  ref <- augp_oracle_fit("pair_sep_jsep_NT_NULL")$fit
  difference <- function(fit, baseline = ref) max(abs(fit$att - baseline$att), na.rm = TRUE)
  differs <- function(fit, baseline = ref) difference(fit, baseline) > 1e-8
  # A selection-only covariate enters the observation logits only: never the
  # imputation/counterfactual OR design or the treatment logit.
  fV <- augp_fit(selection_xformla = ~x + V)
  testthat::expect_true(differs(fV))
  ndV <- fV$chained$nuisance_diagnostics
  testthat::expect_identical(ndV$support$design_columns, c("(Intercept)", "x"))
  testthat::expect_identical(ndV$records$selection$formula_columns, c("(Intercept)", "x", "V"))
  testthat::expect_false("V" %in% fV$chained$outcome_design$columns)
  testthat::expect_false(any(vapply(ndV$treatment, function(tr) "V" %in% names(tr$coefficients), logical(1))))
  testthat::expect_true(differs(augp_fit(selection_xformla = ~V)))
  # With binary x and selection ~x, the finite score spans give the same
  # points in these fixed cases. The accepted independent point helper confirms
  # those equalities; a model change need not move points. Selection-only V
  # supplies material contrasts on these same original observations.
  testthat::expect_lt(difference(augp_fit(selection_model = "independent_period_logit")), 1e-11)
  testthat::expect_lt(difference(augp_fit(selection_pooling = "across_cohorts")), 1e-11)
  testthat::expect_lt(difference(augp_fit(imputation_pooling = "within_cohort")), 1e-11)
  testthat::expect_true(differs(augp_fit(selection_model = "independent_period_logit", selection_xformla = ~x + V), fV))
  testthat::expect_true(differs(augp_fit(selection_pooling = "across_cohorts", selection_xformla = ~x + V), fV))
  testthat::expect_true(differs(augp_fit(imputation_pooling = "within_cohort", selection_xformla = ~x + V), fV))
  # The V-based augmented correction differs from both one-arm corrections on
  # the same data; pure imputation still never forces inactive probability arguments.
  imp <- augp_fit(selection_est_method = "imputation",
    selection_model = quote(stop("inactive selection_model was evaluated")),
    selection_pooling = quote(stop("inactive selection_pooling was evaluated")),
    selection_xformla = quote(stop("inactive selection_xformla was evaluated")))
  testthat::expect_identical(imp$chained$selection$method, "imputation")
  testthat::expect_equal(imp$chained$nuisance_diagnostics$support$observation_probability_models, 0)
  testthat::expect_lt(difference(imp), 1e-11)
  testthat::expect_true(differs(imp, fV))
  w <- augp_fit(selection_est_method = "ipw", imputation_pooling = quote(stop("imputation_pooling was evaluated")))
  testthat::expect_false(identical(w$chained$selection$method, "augmented"))
  testthat::expect_lt(difference(w), 1e-11)
  wV <- augp_fit(selection_est_method = "ipw", selection_xformla = ~x + V,
    imputation_pooling = quote(stop("imputation_pooling was evaluated")))
  testthat::expect_true(differs(wV, fV))
  # Intercept-only formulas for both arms run and are reported as such.
  io <- augp_fit(xformla = ~1, selection_xformla = ~1)
  testthat::expect_identical(io$chained$outcome_design$columns, "(Intercept)")
  testthat::expect_equal(io$chained$nuisance_diagnostics$records$selection$n_slope_columns, 0)
  testthat::expect_equal(io$chained$nuisance_diagnostics$counterfactual_or$rank, 0L)
})

# ------------------------------------------------------------------------------
# E. Full F, frame-only identities, weights, order and chunking
# ------------------------------------------------------------------------------

testthat::test_that("both designs are required on every F row, and frame-only rows carry IF mass", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  ref <- augp_oracle_fit("pair_sep_jsep_NT_NULL")
  fit <- ref$fit
  args <- ref$rq$args
  frame_only <- setdiff(args$frame$id, args$data$id)
  testthat::expect_gt(length(frame_only), 0L)
  testthat::expect_true(all(frame_only %in% rownames(fit$inffunc)))
  av <- fit$chained$available
  testthat::expect_true(any(fit$inffunc[frame_only, av, drop = FALSE] != 0))
  bad_v <- args
  bad_v$frame$V[bad_v$frame$id == frame_only[1L]] <- NA
  bad_v$selection_xformla <- ~x + V
  testthat::expect_error(do.call(cdid::att_gt, bad_v), class = "cdid_chained_request_covariate_error")
  bad_x <- args
  bad_x$frame$x[bad_x$frame$id == frame_only[1L]] <- NA
  testthat::expect_error(do.call(cdid::att_gt, bad_x), class = "cdid_chained_request_covariate_error")
  unknown <- args
  unknown$frame$g <- NULL
  e <- testthat::expect_error(do.call(cdid::att_gt, unknown), class = "cdid_chained_request_input_error")
  testthat::expect_match(conditionMessage(e), "known cohorts", fixed = TRUE)
})

testthat::test_that("row order, global weight scale and chunk size leave augmented results unchanged", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  obj <- augp_oracle_fit("date_sep_jwithin_NYT_xV")
  ref <- obj$fit
  args <- obj$rq$args
  ord <- c(seq.int(2L, nrow(args$frame), 2L), seq.int(1L, nrow(args$frame), 2L))
  sh <- args
  sh$frame <- args$frame[ord, , drop = FALSE]
  sh$data <- args$data[rev(seq_len(nrow(args$data))), , drop = FALSE]
  s <- do.call(cdid::att_gt, sh)
  testthat::expect_identical(rownames(s$inffunc), args$frame$id[ord])
  testthat::expect_equal(s$att, ref$att, tolerance = 1e-10)
  testthat::expect_equal(unname(s$inffunc[args$frame$id, , drop = FALSE]), unname(ref$inffunc), tolerance = 1e-10)
  # Degree-zero normalization: a common weight scale cancels.
  sc <- args
  sc$frame$w <- 3 * args$frame$w
  scaled <- do.call(cdid::att_gt, sc)
  testthat::expect_equal(scaled$att, ref$att, tolerance = 1e-10)
  testthat::expect_equal(unname(scaled$inffunc), unname(ref$inffunc), tolerance = 1e-10)
  # An explicit chunk size replaces automatic row planning without changing values.
  ch_args <- args
  ch_args$control <- list(chunk_rows = 7L)
  ch <- do.call(cdid::att_gt, ch_args)
  testthat::expect_false(ch$chained$nuisance_diagnostics$accounting$automatic_chunk_rows)
  testthat::expect_equal(ch$att, ref$att, tolerance = 1e-10)
  testthat::expect_equal(unname(ch$inffunc), unname(ref$inffunc), tolerance = 1e-10)
})

# ------------------------------------------------------------------------------
# F. Inference consumers: clusters, bands, aggregation, pretest, GMM weighting
# ------------------------------------------------------------------------------

testthat::test_that("clustering, multiplier bands, aggregation and the pretest use the existing consumers", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  obj <- augp_oracle_fit("pair_across_jwithin_NYT_V")
  fit <- obj$fit
  args <- obj$rq$args
  N <- nrow(args$frame)
  cl_args <- args
  cl_args$clustervars <- "cl"
  fit_cl <- do.call(cdid::att_gt, cl_args)
  # Clustering changes inference only.
  testthat::expect_identical(fit_cl$att, fit$att)
  testthat::expect_identical(fit_cl$inffunc, fit$inffunc)
  av <- fit$chained$available
  sums <- rowsum(fit$inffunc[, av, drop = FALSE], args$frame$cl[match(rownames(fit$inffunc), args$frame$id)],
    reorder = FALSE)
  V_cl <- crossprod(sweep(sums, 2L, colMeans(sums), "-")) / N
  testthat::expect_equal(unname(fit_cl$V_analytical[av, av]), unname(V_cl), tolerance = 1e-12)
  testthat::expect_true(fit_cl$chained$inference$clustered)
  # Multiplier draws use the caller RNG stream.
  bs_args <- cl_args
  bs_args$bstrap <- TRUE
  bs_args$cband <- TRUE
  bs_args$biters <- 64L
  set.seed(52711L)
  b1 <- do.call(cdid::att_gt, bs_args)
  set.seed(52711L)
  b2 <- do.call(cdid::att_gt, bs_args)
  testthat::expect_identical(b1$att, fit$att)
  testthat::expect_identical(b1$se, b2$se)
  testthat::expect_identical(b1$c, b2$c)
  testthat::expect_true(is.finite(b1$c) && b1$c > 0)
  for (tp in c("simple", "group", "dynamic", "calendar")) {
    ag <- cdid::aggte(b1, type = tp)
    testthat::expect_s3_class(ag, "cdid_AGGTEobj")
    testthat::expect_true(is.finite(ag$overall.att) && is.finite(ag$overall.se) && ag$overall.se > 0)
    testthat::expect_identical(ag$chained$estimator, fit$chained$estimator)
  }
  # The D21 full declared pretrend test is computed on the public path.
  testthat::expect_identical(fit$chained$pretest$status, "computed")
  testthat::expect_true(is.finite(fit$W) && is.finite(fit$Wpval))
})

testthat::test_that("default feasible two-step weighting and its explicit alias run the augmented consumers", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  # Consumer smoke/self-consistency only: the independent oracle is identity-GMM.
  def <- augp_fit(gmm_weighting = NULL)
  alias <- augp_fit(gmm_weighting = "two_step", gmm_covariance = "tikhonov")
  for (nm in c("att", "se", "V_analytical", "inffunc")) testthat::expect_identical(alias[[nm]], def[[nm]])
  testthat::expect_identical(def$chained$request$gmm_weighting, "optimal")
  testthat::expect_identical(def$chained$covariance_handling$handler_selected, "tikhonov")
  testthat::expect_true(any(def$chained$available))
  testthat::expect_true(all(is.finite(def$inffunc[, def$chained$available, drop = FALSE])))
  testthat::expect_identical(def$chained$nuisance_diagnostics$method, "augmented")
})

# ------------------------------------------------------------------------------
# G. Legitimate failures: no fallback to another correction
# ------------------------------------------------------------------------------

testthat::test_that("an all-success observation cell refuses under the augmented correction (no fallback)", {
  testthat::skip_if_not(augp_has_installed_cdid(), "Installed cdid namespace required.")
  args <- augp_base_args()
  # Balanced panel on the panel roster: every pair is observed by every unit, so
  # each observation logit has no finite root.
  bal <- args
  bal$frame <- args$frame[args$frame$id %in% args$data$id, , drop = FALSE]
  k <- match(args$data$id, bal$frame$id)
  fill <- 0.2 * args$data$time + 0.1 * bal$frame$x[k] + 0.05 * sin(1.7 * k + args$data$time)
  bal$data$y <- ifelse(is.na(args$data$y), fill, args$data$y)
  e <- testthat::expect_error(do.call(cdid::att_gt, bal), class = "cdid_ipw_fit_error")
  testthat::expect_s3_class(e, "cdid_error")
  # Pure imputation fits no observation logit and runs on the same balanced data.
  imp_args <- bal
  imp_args$selection_est_method <- "imputation"
  imp <- do.call(cdid::att_gt, imp_args)
  testthat::expect_identical(imp$chained$selection$method, "imputation")
  testthat::expect_equal(imp$chained$nuisance_diagnostics$accounting$observation_logit_models, 0)
})
