# Deterministic installed-namespace and independent reference tests for D48 Step 2
# Section B: private and public selection-model connection for chained IPW and
# shared one-period outcome-change OR (CONNECTION_BRIEF.md and DATE_SELECTION_SPEC.md).

sp_has_installed_cdid <- function() {
  if (exists("selection_public_helpers", inherits = TRUE)) return(FALSE)
  isNamespaceLoaded("cdid") || requireNamespace("cdid", quietly = TRUE)
}

sp_ns <- function(name) {
  if (exists("selection_public_helpers", inherits = TRUE)) {
    env <- get("selection_public_helpers", inherits = TRUE)
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

sp_val <- function(name) {
  if (exists("selection_public_helpers", inherits = TRUE)) {
    env <- get("selection_public_helpers", inherits = TRUE)
    stopifnot(is.environment(env), exists(name, envir = env, inherits = FALSE))
    return(get(name, envir = env, inherits = FALSE))
  }
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  get(name, ns, inherits = FALSE)
}

sys.source(testthat::test_path("reference", "cp06_pooled_public_or_reference.R"), envir = environment())
sys.source(testthat::test_path("reference", "selection_public_reference.R"), envir = environment())

sp_fx <- sp_ref_build_fixture()
sp_raw <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_baseline", "raw.csv"))
sp_raw$cohort <- as.numeric(replace(sp_raw$cohort, sp_raw$cohort == "NT", "Inf"))
sp_trans <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_corollary1", "transport.csv"))
sp_trans$cohort <- as.numeric(replace(sp_trans$cohort, sp_trans$cohort == "NT", "Inf"))

# Align nuisance coordinates by their scientific cell/block identities. The
# independent dense reference uses traversal order; the accepted backend sorts
# by cohort and calendar endpoints. Neither order changes the numerical test.
sp_selection_reference_order <- function(selection_fit, reference_point) {
  actual_keys <- vapply(selection_fit$cells, function(cell) {
    if (cell$type == "date") paste(cell$cohort, cell$calendar_index, sep = ":") else
      paste(cell$cohort, cell$pair$s_calendar_index, cell$pair$t_calendar_index, sep = ":")
  }, character(1))
  reference_keys <- vapply(reference_point$sel_cells, function(cell) {
    if (cell$kind == "date") paste(cell$cohort, cell$t_idx, sep = ":") else
      paste(cell$cohort, cell$s_idx, cell$t_idx, sep = ":")
  }, character(1))
  stopifnot(setequal(actual_keys, reference_keys))
  actual_block_keys <- vapply(selection_fit$blocks, function(block)
    paste(sort(actual_keys[block$cell_ids]), collapse = "|"), character(1))
  reference_block_keys <- vapply(reference_point$sel_block_groups, function(group)
    paste(sort(reference_keys[group]), collapse = "|"), character(1))
  blocks <- match(reference_block_keys, actual_block_keys)
  stopifnot(!anyNA(blocks), !anyDuplicated(blocks))
  widths <- vapply(selection_fit$blocks, function(block) length(block$beta), integer(1))
  ends <- cumsum(widths)
  beta <- unlist(lapply(blocks, function(m) {
    if (!widths[m]) integer(0) else seq.int(ends[m] - widths[m] + 1L, ends[m])
  }), use.names = FALSE)
  list(alpha = match(reference_keys[unlist(reference_point$sel_block_groups)], actual_keys), beta = beta)
}

sp_public_fixture <- function(est_method = "ipw", route = "F",
                              selection = "conditional_mtar",
                              selection_model = "pair_logit",
                              selection_pooling = "separate",
                              selection_xformla = ~x) {
  d <- sp_raw
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
    tb <- sp_trans[sp_trans$route == route, ]
    frame$rho <- tb$rho[match(paste(d$cohort, d$x), paste(tb$cohort, tb$x))]
  }
  args <- list(
    yname = "y", tname = "time", idname = "id", gname = "g",
    data = data, frame = frame, xformla = ~x, weightsname = "w",
    est_method = est_method, control_group = "nevertreated",
    bstrap = FALSE, cband = FALSE, gmm_weighting = "identity",
    selection = selection,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    selection_xformla = selection_xformla,
    selection_universe = route
  )
  if (route %in% c("B", "A")) {
    args$selection_rho <- "rho"
    args$selection_justification <- "frozen population-law normalized transport ratio"
  }
  list(args = args, frame = frame, data = data)
}

# ------------------------------------------------------------------------------
# 1. Public surface, normalize_selection_declaration, and scope-gate refusals
# ------------------------------------------------------------------------------

testthat::test_that("att_gt signature and chained-only formals include selection_model and selection_pooling with unchanged defaults", {
  att_fn <- sp_ns("att_gt")
  up_formals <- sp_val("upstream_att_gt_formals")
  ch_formals <- sp_val("cdid_chained_only_formals")
  f <- formals(att_fn)
  n_up <- length(up_formals)
  testthat::expect_identical(names(f)[seq_len(n_up)], up_formals)
  testthat::expect_identical(names(f)[n_up + 1L], "...")
  testthat::expect_identical(names(f)[-seq_len(n_up + 1L)], ch_formals)
  testthat::expect_true(all(c("selection", "selection_model", "selection_pooling") %in% ch_formals))
  testthat::expect_identical(eval(f$selection), "none")
  testthat::expect_identical(eval(f$selection_model), "pair_logit")
  testthat::expect_identical(eval(f$selection_pooling), "separate")
  testthat::expect_identical(names(f)[length(f)], "control")
  testthat::expect_null(f$control)
})

testthat::test_that("scope gate validates selection, selection_model, selection_pooling, and universe/transport compatibility before forcing data", {
  att_fn <- sp_ns("att_gt")
  base <- list(
    yname = "y", tname = "time", idname = "id", gname = "g",
    data = quote(stop("data must not be evaluated")),
    est_method = "ipw"
  )
  refuse <- function(cls, ...) {
    e <- testthat::expect_error(do.call(att_fn, utils::modifyList(base, list(...))), class = cls)
    testthat::expect_no_match(conditionMessage(e), "data must not be evaluated")
    invisible(e)
  }

  # Invalid selection / selection_model / selection_pooling values
  for (bad_sel in list("ipw", "smar", "corollary2", "", NA_character_, c("none", "corollary1"))) {
    refuse("cdid_chained_request_unsupported_selection", selection = bad_sel)
  }
  for (bad_mod in list("probit", "pair", "", NA_character_, 1L)) {
    refuse("cdid_chained_request_unsupported_selection", selection = "conditional_mtar", selection_model = bad_mod)
  }
  for (bad_pool in list("pooled", "global", "", NA_character_, 1L)) {
    refuse("cdid_chained_request_unsupported_selection", selection = "conditional_mtar", selection_pooling = bad_pool)
  }

  # selection = "none" refuses any non-default selection_model, selection_pooling, formula, universe, or transport
  refuse("cdid_chained_request_unsupported_selection", selection = "none", selection_model = "independent_period_logit")
  refuse("cdid_chained_request_unsupported_selection", selection = "none", selection_pooling = "within_cohort")
  refuse("cdid_chained_request_unsupported_selection", selection = "none", selection_pooling = "across_cohorts")
  refuse("cdid_chained_request_unsupported_selection", selection = "none", selection_xformla = ~x)
  refuse("cdid_chained_request_unsupported_selection", selection = "none", selection_universe = "B")
  refuse("cdid_chained_request_unsupported_selection", selection = "none", selection_rho = "rho")
  refuse("cdid_chained_request_unsupported_selection", selection = "none", selection_justification = "why")

  # Legacy alias selection = "corollary1" is direct-pair only: refuses independent_period_logit
  refuse("cdid_chained_request_unsupported_selection", selection = "corollary1", selection_model = "independent_period_logit")

  # independent_period_logit is supported only on U = F with rho = NULL and justification = NULL
  for (u in c("B", "A")) {
    refuse("cdid_chained_request_unsupported_selection",
           selection = "conditional_mtar", selection_model = "independent_period_logit",
           selection_universe = u, selection_justification = "breaks date factorization")
  }
  refuse("cdid_chained_request_unsupported_selection",
         selection = "conditional_mtar", selection_model = "independent_period_logit",
         selection_universe = "F", selection_rho = "rho")
  refuse("cdid_chained_request_unsupported_selection",
         selection = "conditional_mtar", selection_model = "independent_period_logit",
         selection_universe = "F", selection_justification = "unneeded on F")

  # pair_logit on B/A still requires justification, and on F refuses rho
  refuse("cdid_chained_request_input_error",
         selection = "conditional_mtar", selection_model = "pair_logit",
         selection_universe = "B")
  refuse("cdid_chained_request_input_error",
         selection = "conditional_mtar", selection_model = "pair_logit",
         selection_universe = "F", selection_rho = "rho")

  # D64 allows custom probability correction; imputation/augmented callbacks
  # still have no complete forcing contract and must refuse before data.
  for (axis in c("imputation", "dr")) {
    refuse("cdid_chained_request_unsupported_selection",
           est_method = function(...) stop("callback forced"), selection = "conditional_mtar",
           selection_est_method = axis)
  }
  for (method in list("reg", "dr", function(...) stop("callback forced"))) {
    testthat::expect_error(
      do.call(att_fn, utils::modifyList(base, list(est_method = method, selection = "conditional_mtar"))),
      "data must not be evaluated"
    )
  }

  # Valid combinations pass the scope gate and reach the data promise
  for (sm in c("pair_logit", "independent_period_logit")) {
    for (sp in c("separate", "within_cohort", "across_cohorts")) {
      testthat::expect_error(
        do.call(att_fn, c(base, list(selection = "conditional_mtar", selection_model = sm, selection_pooling = sp))),
        "data must not be evaluated"
      )
    }
  }
  for (sp in c("separate", "within_cohort", "across_cohorts")) {
    testthat::expect_error(
      do.call(att_fn, c(base, list(selection = "corollary1", selection_model = "pair_logit", selection_pooling = sp))),
      "data must not be evaluated"
    )
  }
})

testthat::test_that("cp04_public_design and chained_request_prepare_inputs preserve intercept/assign metadata and support ~0 selection formulas", {
  prep_fn <- sp_ns("chained_request_prepare_inputs")
  fx <- sp_public_fixture("ipw", route = "F", selection = "conditional_mtar",
                          selection_model = "independent_period_logit", selection_pooling = "across_cohorts")

  # 1. Standard intercept + covariate formula ~ x
  p1 <- prep_fn(
    data = fx$data, yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~x, weightsname = "w", frame = fx$frame,
    selection = "conditional_mtar", selection_model = "independent_period_logit",
    selection_pooling = "across_cohorts", selection_xformla = ~x
  )
  testthat::expect_identical(p1$cp04$selection_requested, "conditional_mtar")
  testthat::expect_identical(p1$cp04$canonical_assumption, "conditional_mtar")
  testthat::expect_false(p1$cp04$legacy_alias)
  testthat::expect_identical(p1$cp04$selection_model, "independent_period_logit")
  testthat::expect_identical(p1$cp04$selection_pooling, "across_cohorts")
  testthat::expect_true(attr(p1$cp04$Z, "intercept"))
  testthat::expect_identical(attr(p1$cp04$Z, "assign"), c(0L, 1L))

  # 2. Intercept-only ~ 1
  p_icpt <- prep_fn(
    data = fx$data, yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~x, weightsname = "w", frame = fx$frame,
    selection = "conditional_mtar", selection_model = "pair_logit",
    selection_pooling = "within_cohort", selection_xformla = ~1
  )
  testthat::expect_equal(ncol(p_icpt$cp04$Z), 1L)
  testthat::expect_true(attr(p_icpt$cp04$Z, "intercept"))
  testthat::expect_identical(attr(p_icpt$cp04$Z, "assign"), 0L)

  # 3. No-intercept ~ 0 + x and backticked `(Intercept)` slope column
  p_noicpt <- prep_fn(
    data = fx$data, yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~x, weightsname = "w", frame = fx$frame,
    selection = "conditional_mtar", selection_model = "pair_logit",
    selection_pooling = "across_cohorts", selection_xformla = ~ 0 + `(Intercept)` + x
  )
  testthat::expect_equal(ncol(p_noicpt$cp04$Z), 2L)
  testthat::expect_false(attr(p_noicpt$cp04$Z, "intercept"))
  testthat::expect_identical(attr(p_noicpt$cp04$Z, "assign"), c(1L, 2L))
  testthat::expect_identical(colnames(p_noicpt$cp04$Z), c("(Intercept)", "x"))

  # 4. Zero-column ~ 0 selection formula is accepted for selection_xformla, while ~ 0 xformla still refuses
  p_zero <- prep_fn(
    data = fx$data, yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~x, weightsname = "w", frame = fx$frame,
    selection = "conditional_mtar", selection_model = "independent_period_logit",
    selection_pooling = "separate", selection_xformla = ~0
  )
  testthat::expect_equal(ncol(p_zero$cp04$Z), 0L)
  testthat::expect_false(attr(p_zero$cp04$Z, "intercept"))
  testthat::expect_identical(attr(p_zero$cp04$Z, "assign"), integer(0))
  testthat::expect_error(
    prep_fn(
      data = fx$data, yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~0, weightsname = "w", frame = fx$frame,
      selection = "conditional_mtar", selection_xformla = ~x
    ),
    class = "cdid_chained_request_covariate_error"
  )
})

# ------------------------------------------------------------------------------
# 2. Pure-table consumer verification across all 3 date modes on F and all 9
#    pair modes on F/B/A against independent two-step central FD reference
# ------------------------------------------------------------------------------

testthat::test_that("IPW and shared-slope OR consumers match independent two-step FD reference across all 3 date modes on F and all 9 pair modes on F/B/A", {
  fit_ipw_tab <- sp_ns("fit_cp04_selection_table")
  fit_or_tab <- sp_ns("fit_or_selection_table")
  mk_recipe <- sp_ns("or_pooled_outcome_recipe")
  recipe_or <- mk_recipe(colnames(sp_fx$V), TRUE, c("x1", "x2"))
  L <- length(sp_fx$links)

  for (sm in c("pair_logit", "independent_period_logit")) {
    for (sp in c("separate", "within_cohort", "across_cohorts")) {
      univs <- if (identical(sm, "independent_period_logit")) "F" else c("F", "B", "A")
      for (univ in univs) {
        U <- sp_fx$U[[univ]]
        rho <- sp_fx$rho[[univ]]

        # Independent two-step (1e-4, 5e-5) central FD for IPW
        ref_ipw <- sp_ref_central_fd(
          sp_fx$N,
          function(cw) sp_ref_ipw_point(sp_fx, cw, sm, sp, univ)
        )
        tab_ipw <- fit_ipw_tab(
          sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
          U = U, rho = rho, universe = univ, calendar = sp_ref_calendar,
          links = sp_fx$links, selection_model = sm, selection_pooling = sp
        )
        testthat::expect_lt(ref_ipw$step_diff, sp_ref_tol)
        testthat::expect_lt(max(abs(tab_ipw$selection$a - ref_ipw$point$a)), 1e-12)
        sel_order <- sp_selection_reference_order(tab_ipw$selection, ref_ipw$point)
        testthat::expect_lt(max(abs(tab_ipw$selection$alpha[sel_order$alpha] - ref_ipw$point$sel_alpha)), 2e-11)
        testthat::expect_lt(max(abs(tab_ipw$selection$beta[sel_order$beta] - ref_ipw$point$sel_beta)), 2e-11)
        testthat::expect_lt(max(abs(tab_ipw$estimates - ref_ipw$point$mu)), 1e-11)
        testthat::expect_lt(max(abs(tab_ipw$phi - ref_ipw$fd_phi)), sp_ref_tol)
        testthat::expect_lt(max(abs(colMeans(tab_ipw$phi))), sp_ref_root_tol)

        # Verify every intermediate IPW coordinate IF (trt_beta, T, M, dp, dm)
        trt_beta_if <- do.call(cbind, lapply(tab_ipw$treatments, `[[`, "beta_if"))
        T_if <- do.call(cbind, lapply(tab_ipw$treatments, function(f) f$coordinates_if[, "T"]))
        M_if <- do.call(cbind, lapply(tab_ipw$treatments, function(f) f$coordinates_if[, "M"]))
        dp_if <- do.call(cbind, lapply(tab_ipw$treatments, function(f) f$coordinates_if[, "dp"]))
        dm_if <- do.call(cbind, lapply(tab_ipw$treatments, function(f) f$coordinates_if[, "dm"]))
        n_sel <- length(ref_ipw$point$a) + length(ref_ipw$point$sel_alpha) + length(ref_ipw$point$sel_beta)
        n_tb <- ncol(trt_beta_if)
        testthat::expect_lt(max(abs(trt_beta_if - ref_ipw$fd_stack[, n_sel + seq_len(n_tb), drop = FALSE])), sp_ref_tol)
        testthat::expect_lt(max(abs(T_if - ref_ipw$fd_stack[, n_sel + n_tb + seq_len(L), drop = FALSE])), sp_ref_tol)
        testthat::expect_lt(max(abs(M_if - ref_ipw$fd_stack[, n_sel + n_tb + L + seq_len(L), drop = FALSE])), sp_ref_tol)
        testthat::expect_lt(max(abs(dp_if - ref_ipw$fd_stack[, n_sel + n_tb + 2L * L + seq_len(L), drop = FALSE])), sp_ref_tol)
        testthat::expect_lt(max(abs(dm_if - ref_ipw$fd_stack[, n_sel + n_tb + 3L * L + seq_len(L), drop = FALSE])), sp_ref_tol)

        # Independent two-step (1e-4, 5e-5) central FD for shared-slope OR
        ref_or <- sp_ref_central_fd(
          sp_fx$N,
          function(cw) sp_ref_or_point(sp_fx, cw, sm, sp, univ)
        )
        tab_or <- fit_or_tab(
          sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
          U = U, rho = rho, universe = univ, calendar = sp_ref_calendar,
          links = sp_fx$links, outcome_recipe = recipe_or,
          selection_model = sm, selection_pooling = sp
        )
        testthat::expect_lt(ref_or$step_diff, sp_ref_tol)
        testthat::expect_lt(max(abs(tab_or$shared$beta - ref_or$point$or_beta)), 1e-11)
        or_order <- vapply(tab_or$links[ref_or$point$or_block_first_link], `[[`, integer(1), "regression")
        or_al <- vapply(tab_or$regressions, `[[`, numeric(1), "alpha")[or_order]
        testthat::expect_lt(max(abs(or_al - ref_or$point$or_alpha)), 1e-11)
        testthat::expect_lt(max(abs(tab_or$estimates - ref_or$point$mu)), 1e-11)
        testthat::expect_lt(max(abs(tab_or$phi - ref_or$fd_phi)), sp_ref_tol)
        testthat::expect_lt(max(abs(colMeans(tab_or$phi))), sp_ref_root_tol)
        G_or <- length(or_al)
        p_or <- length(ref_or$point$or_beta)
        testthat::expect_lt(max(abs(tab_or$alpha_if[, or_order, drop = FALSE] - ref_or$fd_stack[, n_sel + seq_len(G_or), drop = FALSE])), sp_ref_tol)
        testthat::expect_lt(max(abs(tab_or$beta_if - ref_or$fd_stack[, n_sel + G_or + seq_len(p_or), drop = FALSE])), sp_ref_tol)
      }
    }
  }
})

testthat::test_that("intercept-only (~1), no-intercept (~0 + z1 + z2), and zero-column (~0) selection designs match independent FD reference", {
  fit_ipw_tab <- sp_ns("fit_cp04_selection_table")
  fit_or_tab <- sp_ns("fit_or_selection_table")
  mk_recipe <- sp_ns("or_pooled_outcome_recipe")
  recipe_or <- mk_recipe(colnames(sp_fx$V), TRUE, c("x1", "x2"))

  Z_icpt <- sp_fx$Z[, "(Intercept)", drop = FALSE]
  attr(Z_icpt, "intercept") <- TRUE
  attr(Z_icpt, "assign") <- 0L

  Z_noicpt <- sp_fx$Z[, c("z1", "z2"), drop = FALSE]
  attr(Z_noicpt, "intercept") <- FALSE
  attr(Z_noicpt, "assign") <- c(1L, 2L)

  Z_zero <- sp_fx$Z[, integer(0), drop = FALSE]
  attr(Z_zero, "intercept") <- FALSE
  attr(Z_zero, "assign") <- integer(0)

  variants <- list(
    list(name = "icpt_only", Z = Z_icpt, has_icpt = TRUE, zero_col = FALSE, slopes = character(0)),
    list(name = "no_icpt", Z = Z_noicpt, has_icpt = FALSE, zero_col = FALSE, slopes = c("z1", "z2")),
    list(name = "zero_col", Z = Z_zero, has_icpt = FALSE, zero_col = TRUE, slopes = character(0))
  )

  for (v in variants) {
    for (sm in c("pair_logit", "independent_period_logit")) {
      for (sp in c("separate", "within_cohort", "across_cohorts")) {
        ref_ipw <- sp_ref_central_fd(
          sp_fx$N,
          function(cw) sp_ref_ipw_point(sp_fx, cw, sm, sp, "F", has_intercept = v$has_icpt, zero_col = v$zero_col, slope_cols = v$slopes)
        )
        tab_ipw <- fit_ipw_tab(
          v$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
          U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
          links = sp_fx$links, selection_model = sm, selection_pooling = sp,
          has_intercept = v$has_icpt
        )
        testthat::expect_lt(max(abs(tab_ipw$estimates - ref_ipw$point$mu)), 1e-11)
        testthat::expect_lt(max(abs(tab_ipw$phi - ref_ipw$fd_phi)), sp_ref_tol)

        ref_or <- sp_ref_central_fd(
          sp_fx$N,
          function(cw) sp_ref_or_point(sp_fx, cw, sm, sp, "F", has_intercept = v$has_icpt, zero_col = v$zero_col, slope_cols = v$slopes)
        )
        tab_or <- fit_or_tab(
          v$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
          U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
          links = sp_fx$links, outcome_recipe = recipe_or,
          selection_model = sm, selection_pooling = sp,
          has_intercept = v$has_icpt
        )
        testthat::expect_lt(max(abs(tab_or$estimates - ref_or$point$mu)), 1e-11)
        testthat::expect_lt(max(abs(tab_or$phi - ref_or$fd_phi)), sp_ref_tol)
      }
    }
  }
})

# ------------------------------------------------------------------------------
# 3. Omission witnesses: endpoint, within-cohort other-cell, cross-cohort, off-A
# ------------------------------------------------------------------------------

testthat::test_that("omitting date endpoints, within-cohort shared cells, cross-cohort shared slopes, or off-A units fails by > 100 * tol", {
  fit_ipw_tab <- sp_ns("fit_cp04_selection_table")
  fit_or_tab <- sp_ns("fit_or_selection_table")
  mk_recipe <- sp_ns("or_pooled_outcome_recipe")
  recipe_or <- mk_recipe(colnames(sp_fx$V), TRUE, c("x1", "x2"))

  # A. Date endpoint omission witness on independent_period_logit (F)
  ipw_date <- fit_ipw_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = sp_fx$links, selection_model = "independent_period_logit", selection_pooling = "separate"
  )
  or_date <- fit_or_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = sp_fx$links, outcome_recipe = recipe_or,
    selection_model = "independent_period_logit", selection_pooling = "separate"
  )
  # Link 1 (NYT12) uses pair (1, 2): s_calendar_index = 1, t_calendar_index = 2.
  # Identify upstream columns belonging to date s=1 vs date t=2.
  s_cells <- which(vapply(ipw_date$selection$cells, function(c) c$calendar_index == 1L, logical(1)))
  t_cells <- which(vapply(ipw_date$selection$cells, function(c) c$calendar_index == 2L, logical(1)))
  s_uids <- unlist(ipw_date$up_map$theta_coords[s_cells])
  t_uids <- unlist(ipw_date$up_map$theta_coords[t_cells])
  f1 <- ipw_date$treatments[[1L]]
  up_f1 <- ipw_date$upstream_if[, f1$upstream_ids, drop = FALSE]
  s_pos <- which(f1$upstream_ids %in% s_uids)
  t_pos <- which(f1$upstream_ids %in% t_uids)
  omit_ipw_uids <- function(pos) {
    up_mod <- up_f1
    up_mod[, pos] <- 0
    xb_mod <- t(sp_ns("cp04_solve")(f1$H_chol, t(up_mod %*% t(f1$cross))))
    coord_delta <- (up_f1 - up_mod) %*% t(f1$partials$upstream) +
      (t(sp_ns("cp04_solve")(f1$H_chol, t(up_f1 %*% t(f1$cross)))) - xb_mod) %*% t(f1$partials$beta)
    max(abs(coord_delta[, 1L] - coord_delta[, 2L]))
  }
  testthat::expect_gt(omit_ipw_uids(s_pos), 100 * sp_ref_tol)
  testthat::expect_gt(omit_ipw_uids(t_pos), 100 * sp_ref_tol)

  # Same endpoint omission check for OR on Link 1
  s_uids_or <- unlist(or_date$up_map$theta_coords[s_cells])
  t_uids_or <- unlist(or_date$up_map$theta_coords[t_cells])
  testthat::expect_gt(max(abs(or_date$upstream_if[, s_uids_or, drop = FALSE] %*%
                                t(or_date$shared$forcing_beta[, s_uids_or, drop = FALSE]))), 100 * sp_ref_tol)
  testthat::expect_gt(max(abs(or_date$upstream_if[, t_uids_or, drop = FALSE] %*%
                                t(or_date$shared$forcing_beta[, t_uids_or, drop = FALSE]))), 100 * sp_ref_tol)

  # B. Within-cohort other-cell omission witness:
  # Under within_cohort pooling, pair (2, 4) in cohort g=3 shares slope beta_1 with pair (1, 2) in cohort g=3.
  # Zeroing out the cell score of pair (2, 4) changes the IF of link NT12 (which only evaluates pair (1, 2)).
  ipw_wc <- fit_ipw_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = sp_fx$links, selection_model = "pair_logit", selection_pooling = "within_cohort"
  )
  bl1 <- ipw_wc$selection$blocks[[1L]]
  cs1 <- ipw_wc$selection$cell_scores[, bl1$cell_ids, drop = FALSE]
  # Cell 1 in block 1 is pair (1,2) in g=3; Cell 2 in block 1 is pair (2,4) in g=3
  cs1_omit_other <- cs1
  cs1_omit_other[, -1L] <- 0
  Z_sl <- sp_fx$Z[, c("z1", "z2"), drop = FALSE]
  bif_full <- sp_ns("evaluate_selection_block_if")(bl1, Z_sl, cs1, sp_fx$N)
  bif_omit <- sp_ns("evaluate_selection_block_if")(bl1, Z_sl, cs1_omit_other, sp_fx$N)
  tc1 <- ipw_wc$up_map$theta_coords[[1L]]
  up_wc_mod <- ipw_wc$upstream_if
  up_wc_mod[, tc1] <- bif_omit$xi_theta
  f3 <- ipw_wc$treatments[[3L]] # NT12 (pair 1,2, cohorts g=3 and Inf)
  u3 <- f3$upstream_ids
  up_diff <- ipw_wc$upstream_if[, u3, drop = FALSE] - up_wc_mod[, u3, drop = FALSE]
  xb_diff <- t(sp_ns("cp04_solve")(f3$H_chol, t(up_diff %*% t(f3$cross))))
  phi3_diff <- up_diff %*% (f3$partials$upstream[1L, ] - f3$partials$upstream[2L, ]) +
    xb_diff %*% (f3$partials$beta[1L, ] - f3$partials$beta[2L, ])
  testthat::expect_gt(max(abs(phi3_diff)), 100 * sp_ref_tol)

  # Same within-cohort other-cell omission witness for OR on two NT links (NT12 and NT13)
  links_nt_two <- sp_fx$links[c(1L, 2L)] # NYT12 (1,2) and NYT24 (2,4)
  or_wc <- fit_or_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = links_nt_two, outcome_recipe = recipe_or,
    selection_model = "pair_logit", selection_pooling = "within_cohort"
  )
  # Omitting pair (2,4) selection scores in cohort k=2 (g=5) changes link 1 (pair 1,2) OR IF
  bl2_or <- or_wc$selection$blocks[[2L]]
  cs2_or <- or_wc$selection$cell_scores[, bl2_or$cell_ids, drop = FALSE]
  cs2_or_omit <- cs2_or; cs2_or_omit[, -1L] <- 0
  bif2_omit <- sp_ns("evaluate_selection_block_if")(bl2_or, Z_sl, cs2_or_omit, sp_fx$N)
  tc2_or <- or_wc$up_map$theta_coords[[2L]]
  delta_up_or <- matrix(0, sp_fx$N, or_wc$up_map$D)
  delta_up_or[, tc2_or] <- or_wc$upstream_if[, tc2_or, drop = FALSE] - bif2_omit$xi_theta
  delta_rhs_beta <- delta_up_or %*% t(or_wc$shared$forcing_beta)
  testthat::expect_gt(max(abs(delta_rhs_beta)), 100 * sp_ref_tol)

  # C. Cross-cohort omission witness under across_cohorts pooling:
  # Link NT12 uses only cohorts g=3 (k=1) and g=Inf (k=3). Under separate and within_cohort,
  # cohort k=2 (g=5) has identically zero IF on NT12 in IPW (and in OR when fitted on NT-only links).
  # Under across_cohorts, cohort k=2's selection scores couple into the shared selection slope beta,
  # making phi[k == 2, "NT12"] strictly nonzero (> 100 * tol).
  ipw_ac <- fit_ipw_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = sp_fx$links, selection_model = "pair_logit", selection_pooling = "across_cohorts"
  )
  ipw_ac_date <- fit_ipw_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = sp_fx$links, selection_model = "independent_period_logit", selection_pooling = "across_cohorts"
  )
  testthat::expect_equal(max(abs(ipw_wc$phi[sp_fx$dat$k == 2L, "NT12"])), 0)
  testthat::expect_gt(max(abs(ipw_ac$phi[sp_fx$dat$k == 2L, "NT12"])), 100 * sp_ref_tol)
  testthat::expect_gt(max(abs(ipw_ac_date$phi[sp_fx$dat$k == 2L, "NT12"])), 100 * sp_ref_tol)

  links_nt_only <- sp_fx$links[c(3L, 4L)] # NT12 (g=3, K=Inf) and NT13 (g=5, K=Inf)
  or_nt_wc <- fit_or_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = links_nt_only, outcome_recipe = recipe_or,
    selection_model = "pair_logit", selection_pooling = "within_cohort"
  )
  or_nt_ac <- fit_or_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$F, rho = sp_fx$rho$F, universe = "F", calendar = sp_ref_calendar,
    links = links_nt_only, outcome_recipe = recipe_or,
    selection_model = "pair_logit", selection_pooling = "across_cohorts"
  )
  testthat::expect_equal(max(abs(or_nt_wc$phi[sp_fx$dat$k == 2L, "NT12"])), 0)
  testthat::expect_gt(max(abs(or_nt_ac$phi[sp_fx$dat$k == 2L, "NT12"])), 100 * sp_ref_tol)

  # D. Off-A omission witness: singletons (B \ A) and zero-observers (F \ B) carry nonzero IF on F and B
  sing_rows <- which(sp_fx$B == 1 & sp_fx$A == 0)
  zero_rows <- which(sp_fx$B == 0)
  testthat::expect_gt(max(abs(ipw_date$phi[sing_rows, ])), 100 * sp_ref_tol)
  testthat::expect_gt(max(abs(ipw_date$phi[zero_rows, ])), 100 * sp_ref_tol)
  testthat::expect_gt(max(abs(or_date$phi[sing_rows, ])), 100 * sp_ref_tol)
  testthat::expect_gt(max(abs(or_date$phi[zero_rows, ])), 100 * sp_ref_tol)

  ipw_b <- fit_ipw_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$B, rho = sp_fx$rho$B, universe = "B", calendar = sp_ref_calendar,
    links = sp_fx$links, selection_model = "pair_logit", selection_pooling = "within_cohort"
  )
  or_b <- fit_or_tab(
    sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
    U = sp_fx$U$B, rho = sp_fx$rho$B, universe = "B", calendar = sp_ref_calendar,
    links = sp_fx$links, outcome_recipe = recipe_or,
    selection_model = "pair_logit", selection_pooling = "within_cohort"
  )
  testthat::expect_gt(max(abs(ipw_b$phi[sing_rows, ])), 100 * sp_ref_tol)
  testthat::expect_equal(max(abs(ipw_b$phi[zero_rows, ])), 0)
  testthat::expect_gt(max(abs(or_b$phi[sing_rows, ])), 100 * sp_ref_tol)
  testthat::expect_equal(max(abs(or_b$phi[zero_rows, ])), 0)
})

# ------------------------------------------------------------------------------
# 4. Installed-namespace private owners, batching, covariance streaming, and
#    full-F cluster totals (including zero-eligible A=0 clusters)
# ------------------------------------------------------------------------------

testthat::test_that("installed owners match FD reference across all 12 modes, preserve legacy separate/pair_logit consumer, and are batch/cluster invariant", {
  testthat::skip_if_not(sp_has_installed_cdid(), "Real installed cdid owners required.")
  ctrl <- sp_ns("ipw_owner_control")(
    max_chunk_rows = 23L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7
  )
  gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)

  private_long <- sp_fx$long
  private_long$g[private_long$g == 0] <- Inf
  panel <- sp_ns("construct_ipw_panel_owner")(
    private_long, "id", "time", "y", "g", ctrl
  )
  frame <- sp_ns("construct_ipw_frame_owner")(
    panel, sp_fx$dat$id, sp_fx$V, weights = sp_fx$dat$w,
    control = ctrl, frame_g = sp_fx$dat$g, allow_frame_only = TRUE
  )
  testthat::expect_gt(frame$n_frame_only_units, 0L)
  testthat::expect_gt(panel$n_singleton_units, 0L)
  testthat::expect_gt(panel$n_unobserved_units, 0L)

  views <- lapply(sp_fx$links, function(lk) {
    sp_ns("construct_ipw_link_view")(
      frame, lk$s, lk$t, g = lk$g, K = lk$K, control = ctrl,
      link_id = lk$link_id, functional_key = paste0("cmpcat:sp:", lk$link_id)
    )
  })
  recipe_or <- sp_ns("or_pooled_outcome_recipe")(colnames(frame$Z_a), TRUE, c("x1", "x2"))

  # Cluster assignment with explicit zero-eligible (A=0) clusters
  cl_zero_a <- ifelse(sp_fx$A == 0, 1000L + sp_fx$dat$id, sp_fx$dat$cluster)

  for (sm in c("pair_logit", "independent_period_logit")) {
    for (sp in c("separate", "within_cohort", "across_cohorts")) {
      univs <- if (identical(sm, "independent_period_logit")) "F" else c("F", "B", "A")
      for (univ in univs) {
        U <- sp_fx$U[[univ]]
        rho <- sp_fx$rho[[univ]]
        selection_Z <- sp_fx$Z
        selection_Z[!U, ] <- NA_real_
        reader <- sp_ns("construct_cp04_frame_reader")(
          frame, sp_fx$dat$id, selection_Z, sp_fx$dat$g, sp_fx$dat$w, ctrl,
          U = U, rho = rho, universe = univ,
          justification = if (univ == "F") NULL else "known normalized population-law transport ratio"
        )

        # Legacy consumer dispatch check: positive-column pair_logit + separate uses legacy consumer
        is_legacy <- identical(sm, "pair_logit") && identical(sp, "separate")
        testthat::expect_identical(
          sp_ns("selection_use_legacy_consumer")(reader, sm, sp),
          is_legacy
        )

        # Fit IPW and OR via private entry points
        fit_ipw <- sp_ns("fit_cp04_corollary1")(
          reader, views, ctrl, gc, selection_model = sm, selection_pooling = sp
        )
        fit_or <- sp_ns("fit_or_corollary1")(
          reader, views, ctrl, gc, rank_tol = 1e-10, outcome_recipe = recipe_or,
          selection_model = sm, selection_pooling = sp
        )
        testthat::expect_identical(isTRUE(fit_ipw$uses_selection_models), !is_legacy)
        testthat::expect_identical(isTRUE(fit_or$uses_selection_models), !is_legacy)

        # Compare against independent table/FD reference
        tab_ipw <- sp_ns("fit_cp04_selection_table")(
          sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
          U = U, rho = rho, universe = univ, calendar = sp_ref_calendar,
          links = sp_fx$links, selection_model = sm, selection_pooling = sp
        )
        tab_or <- sp_ns("fit_or_selection_table")(
          sp_fx$Z, sp_fx$V, sp_fx$Y, sp_fx$obs, sp_fx$dat$g, sp_fx$dat$w,
          U = U, rho = rho, universe = univ, calendar = sp_ref_calendar,
          links = sp_fx$links, outcome_recipe = recipe_or,
          selection_model = sm, selection_pooling = sp
        )

        # Replay IFs across different chunk sizes and link batch sizes
        collect_ipw <- function(ch, lb, scale = "frame", details = FALSE) {
          phi <- matrix(NA_real_, sp_fx$N, length(views))
          sp_ns("replay_cp04_if")(fit_ipw, ch, lb, function(b) {
            phi[b$rows, b$links] <<- b$phi
          }, scale = scale, details = details)
          phi
        }
        collect_or <- function(ch, lb, scale = "frame", details = FALSE) {
          phi <- matrix(NA_real_, sp_fx$N, length(views))
          sp_ns("replay_or_corollary1_if")(fit_or, ch, lb, function(b) {
            phi[b$rows, b$links] <<- b$phi
          }, scale = scale, details = details)
          phi
        }

        phi_ipw_23_4 <- collect_ipw(23L, 4L, "frame", details = TRUE)
        phi_ipw_7_1 <- collect_ipw(7L, 1L, "frame", details = FALSE)
        phi_ipw_core <- collect_ipw(13L, 2L, "core", details = FALSE)

        phi_or_23_4 <- collect_or(23L, 4L, "frame", details = TRUE)
        phi_or_7_1 <- collect_or(7L, 1L, "frame", details = FALSE)
        phi_or_core <- collect_or(13L, 2L, "core", details = FALSE)

        testthat::expect_lt(max(abs(phi_ipw_23_4 - tab_ipw$phi)), 1e-12)
        testthat::expect_lt(max(abs(phi_ipw_7_1 - phi_ipw_23_4)), 1e-13)
        testthat::expect_lt(max(abs(phi_or_23_4 - tab_or$phi)), 1e-12)
        testthat::expect_lt(max(abs(phi_or_7_1 - phi_or_23_4)), 1e-13)

        lambda <- reader$state$lambda
        testthat::expect_lt(max(abs(phi_ipw_core - lambda * phi_ipw_23_4)), 1e-13)
        testthat::expect_lt(max(abs(phi_or_core - lambda * phi_or_23_4)), 1e-13)

        # Streamed covariance matches centered crossprod on core scale
        cov_ipw <- sp_ns("cp04_stream_covariance")(fit_ipw, seq_along(views))
        cov_or <- sp_ns("cp04_stream_covariance")(fit_or, seq_along(views))
        testthat::expect_lt(max(abs(cov_ipw$estimates - tab_ipw$estimates)), 1e-12)
        testthat::expect_lt(max(abs(cov_or$estimates - tab_or$estimates)), 1e-12)
        v_ipw_expected <- crossprod(scale(phi_ipw_core, scale = FALSE)) / reader$state$n_A
        v_or_expected <- crossprod(scale(phi_or_core, scale = FALSE)) / reader$state$n_A
        testthat::expect_lt(max(abs(cov_ipw$covariance - v_ipw_expected)), 1e-12)
        testthat::expect_lt(max(abs(cov_or$covariance - v_or_expected)), 1e-12)

        # Full-F cluster totals including zero-eligible (A=0) clusters match between frame and core scales
        cs_ipw_f <- rowsum(phi_ipw_23_4, cl_zero_a)
        cs_ipw_a <- rowsum(phi_ipw_core, cl_zero_a)
        testthat::expect_lt(
          max(abs(crossprod(scale(cs_ipw_f, scale = FALSE)) / sp_fx$N^2 -
                    crossprod(scale(cs_ipw_a, scale = FALSE)) / reader$state$n_A^2)),
          1e-13
        )
        cs_or_f <- rowsum(phi_or_23_4, cl_zero_a)
        cs_or_a <- rowsum(phi_or_core, cl_zero_a)
        testthat::expect_lt(
          max(abs(crossprod(scale(cs_or_f, scale = FALSE)) / sp_fx$N^2 -
                    crossprod(scale(cs_or_a, scale = FALSE)) / reader$state$n_A^2)),
          1e-13
        )
      }
    }
  }
})

# ------------------------------------------------------------------------------
# 5. Public att_gt / aggte integration: canonical vs legacy alias parity, all 12
#    modes for IPW and OR, NT/NYT, PT-post/PT-all, base_period, clustering, bootstrap
# ------------------------------------------------------------------------------

testthat::test_that("pure-table consumers on cp06_or_baseline/raw.csv + transport.csv match sp_ref_public_reference across all 12 modes for NT and NYT", {
  fit_ipw_tab <- sp_ns("fit_cp04_selection_table")
  fit_or_tab <- sp_ns("fit_or_selection_table")
  mk_recipe <- sp_ns("or_pooled_outcome_recipe")
  Y <- as.matrix(sp_raw[paste0("Y", 1:4)])
  obs <- as.matrix(sp_raw[paste0("R", 1:4)])
  Z <- cbind(`(Intercept)` = 1, x = sp_raw$x)
  attr(Z, "intercept") <- TRUE; attr(Z, "assign") <- c(0L, 1L)
  V <- Z
  recipe_or <- mk_recipe(colnames(V), TRUE, "x")

  for (cg in c("nevertreated", "notyettreated")) {
    cat_info <- cp06_pooled_public_catalog(sp_raw, dates = 1:4, control_group = cg)
    links_df <- cat_info$links_df
    links <- lapply(seq_len(nrow(links_df)), function(i) {
      list(
        link_id = paste0(links_df$role[i], "_", i),
        s = links_df$s[i], t = links_df$t[i], g = links_df$g[i],
        K = as.numeric(strsplit(links_df$controls[i], "|", fixed = TRUE)[[1L]])
      )
    })
    for (sm in c("pair_logit", "independent_period_logit")) {
      for (sp in c("separate", "within_cohort", "across_cohorts")) {
        routes <- if (identical(sm, "independent_period_logit")) "F" else c("F", "B", "A")
        for (route in routes) {
          U <- if (route == "F") rep(TRUE, nrow(sp_raw)) else if (route == "B") (sp_raw$B == 1) else (sp_raw$A == 1)
          tb <- sp_trans[sp_trans$route == route, ]
          raw_r <- sp_raw
          raw_r$rho <- tb$rho[match(paste(sp_raw$cohort, sp_raw$x), paste(tb$cohort, tb$x))]

          tab_ipw <- fit_ipw_tab(Z, V, Y, obs, raw_r$cohort, raw_r$w, U, raw_r$rho, route, 1:4, links, sm, sp)
          ref_ipw <- sp_ref_public_reference(raw_r, est_method = "ipw", selection_model = sm, selection_pooling = sp, route = route, control_group = cg)
          testthat::expect_lt(max(abs(tab_ipw$estimates - ref_ipw$mu)), 1e-11)
          testthat::expect_lt(max(abs(tab_ipw$phi - ref_ipw$phi_mu)), sp_ref_tol)

          tab_or <- fit_or_tab(Z, V, Y, obs, raw_r$cohort, raw_r$w, U, raw_r$rho, route, 1:4, links, recipe_or, sm, sp)
          ref_or <- sp_ref_public_reference(raw_r, est_method = "reg", selection_model = sm, selection_pooling = sp, route = route, control_group = cg)
          testthat::expect_lt(max(abs(tab_or$estimates - ref_or$mu)), 1e-11)
          testthat::expect_lt(max(abs(tab_or$phi - ref_or$phi_mu)), sp_ref_tol)
        }
      }
    }
  }
})

testthat::test_that("canonical selection = 'conditional_mtar' with pair_logit/separate has exact parity with legacy alias selection = 'corollary1'", {
  testthat::skip_if_not(sp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- sp_ns("att_gt")
  agg_fn <- sp_ns("aggte")

  for (em in c("ipw", "reg")) {
    for (route in c("F", "B", "A")) {
      fx_canon <- sp_public_fixture(
        est_method = em, route = route,
        selection = "conditional_mtar", selection_model = "pair_logit", selection_pooling = "separate"
      )
      fx_alias <- sp_public_fixture(
        est_method = em, route = route,
        selection = "corollary1", selection_model = "pair_logit", selection_pooling = "separate"
      )
      fit_canon <- do.call(att_fn, fx_canon$args)
      fit_alias <- do.call(att_fn, fx_alias$args)

      testthat::expect_identical(fit_canon$att, fit_alias$att)
      testthat::expect_identical(fit_canon$se, fit_alias$se)
      testthat::expect_identical(fit_canon$V_analytical, fit_alias$V_analytical)
      testthat::expect_identical(fit_canon$inffunc, fit_alias$inffunc)
      testthat::expect_identical(fit_canon$chained$selection$canonical_assumption, "conditional_mtar")
      testthat::expect_identical(fit_alias$chained$selection$canonical_assumption, "conditional_mtar")
      testthat::expect_identical(fit_canon$chained$selection$requested_selection, "conditional_mtar")
      testthat::expect_identical(fit_alias$chained$selection$requested_selection, "corollary1")
      testthat::expect_identical(fit_canon$chained$selection$selection_model, "pair_logit")
      testthat::expect_identical(fit_canon$chained$selection$selection_pooling, "separate")

      for (tp in c("simple", "group", "dynamic", "calendar")) {
        ag_c <- agg_fn(fit_canon, type = tp)
        ag_a <- agg_fn(fit_alias, type = tp)
        testthat::expect_identical(ag_c$overall.att, ag_a$overall.att)
        testthat::expect_identical(ag_c$overall.se, ag_a$overall.se)
        testthat::expect_identical(ag_c$att.egt, ag_a$att.egt)
        testthat::expect_identical(ag_c$se.egt, ag_a$se.egt)
      }
    }
  }
})

testthat::test_that("public att_gt and aggte match independent public reference across all 3 date modes on F and all 9 pair modes on F/B/A for both IPW and OR", {
  testthat::skip_if_not(sp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- sp_ns("att_gt")
  agg_fn <- sp_ns("aggte")

  for (em in c("ipw", "reg")) {
    for (sm in c("pair_logit", "independent_period_logit")) {
      for (sp in c("separate", "within_cohort", "across_cohorts")) {
        routes <- if (identical(sm, "independent_period_logit")) "F" else c("F", "B", "A")
        for (route in routes) {
          fx <- sp_public_fixture(
            est_method = em, route = route,
            selection = "conditional_mtar", selection_model = sm, selection_pooling = sp
          )
          fit <- do.call(att_fn, fx$args)
          raw_r <- sp_raw
          if (route %in% c("B", "A")) raw_r$rho <- fx$frame$rho
          ref <- sp_ref_public_reference(
            raw_r, est_method = em, selection_model = sm, selection_pooling = sp,
            route = route, control_group = "nevertreated"
          )

          testthat::expect_equal(fit$att, ref$att, tolerance = 1e-10)
          testthat::expect_equal(unname(fit$inffunc), unname(ref$inffunc), tolerance = 1e-7)
          testthat::expect_equal(unname(fit$V_analytical), unname(ref$V_analytical), tolerance = 1e-7)
          testthat::expect_equal(fit$se, sqrt(diag(ref$V_analytical) / fit$n), tolerance = 1e-8)
          testthat::expect_identical(rownames(fit$inffunc), fx$frame$id)
          testthat::expect_equal(fit$n, nrow(fx$frame))

          sel_rec <- fit$chained$selection
          testthat::expect_identical(sel_rec$canonical_assumption, "conditional_mtar")
          testthat::expect_identical(sel_rec$selection_model, sm)
          testthat::expect_identical(sel_rec$selection_pooling, sp)
          testthat::expect_identical(sel_rec$fitting_universe, route)
          testthat::expect_equal(sel_rec$n_active_cells, ref$n_active_cells)
          testthat::expect_equal(sel_rec$n_selection_blocks, ref$n_selection_blocks)
          testthat::expect_equal(unname(sel_rec$empirical_repeat_rates), unname(ref$a_U), tolerance = 1e-12)

          for (tp in c("simple", "group", "dynamic", "calendar")) {
            ag <- agg_fn(fit, type = tp)
            testthat::expect_equal(ag$overall.att, unname(ref$agg_overall[tp]), tolerance = 1e-10)
            if_ag <- as.numeric(switch(tp, simple = ag$inf.function$simple.att,
              group = ag$inf.function$selective.inf.func, calendar = ag$inf.function$calendar.inf.func,
              dynamic = ag$inf.function$dynamic.inf.func))
            testthat::expect_equal(if_ag, as.numeric(ref$phi_agg_overall[, tp]), tolerance = 1e-7)
            se_ref <- sqrt(sum((if_ag - mean(if_ag))^2) / (fit$n^2))
            testthat::expect_equal(ag$overall.se, se_ref, tolerance = 1e-8)
          }
        }
      }
    }
  }
})

testthat::test_that("public NYT, PT-all, universal base_period, formula variants (~1, ~0+x, ~0), clustering, and multiplier bootstrap work across selection modes", {
  testthat::skip_if_not(sp_has_installed_cdid(), "Real installed cdid owners required.")
  att_fn <- sp_ns("att_gt")
  agg_fn <- sp_ns("aggte")

  # 1. Not-yet-treated, PT-all (pt_start = 2), and universal base_period on pooled date and pair modes
  for (em in c("ipw", "reg")) {
    for (cfg in list(
      list(sm = "independent_period_logit", sp = "within_cohort", route = "F"),
      list(sm = "pair_logit", sp = "across_cohorts", route = "B")
    )) {
      fx <- sp_public_fixture(
        est_method = em, route = cfg$route,
        selection = "conditional_mtar", selection_model = cfg$sm, selection_pooling = cfg$sp
      )
      raw_r <- sp_raw
      if (cfg$route %in% c("B", "A")) raw_r$rho <- fx$frame$rho

      for (bp in c("varying", "universal")) {
        for (pt in c("post", "all")) {
          ps <- if (pt == "all") 2 else NULL
          fit <- do.call(att_fn, utils::modifyList(fx$args, list(
            control_group = "notyettreated", base_period = bp,
            pt_assumption = pt, pt_start = ps
          )))
          ref <- sp_ref_public_reference(
            raw_r, est_method = em, selection_model = cfg$sm, selection_pooling = cfg$sp,
            route = cfg$route, control_group = "notyettreated",
            base_period = bp, pt_assumption = pt, pt_start = ps
          )
          testthat::expect_equal(fit$att, ref$att, tolerance = 1e-10)
          avail <- !is.na(ref$att)
          testthat::expect_equal(unname(fit$inffunc[, avail, drop = FALSE]),
                                 unname(ref$inffunc[, avail, drop = FALSE]), tolerance = 1e-7)
          testthat::expect_equal(unname(fit$V_analytical[avail, avail, drop = FALSE]),
                                 unname(ref$V_analytical[avail, avail, drop = FALSE]), tolerance = 1e-7)
        }
      }
    }
  }

  # 2. Formula variants (~1, ~0 + x, ~0) on public att_gt
  for (em in c("ipw", "reg")) {
    fx <- sp_public_fixture(
      est_method = em, route = "F",
      selection = "conditional_mtar", selection_model = "independent_period_logit",
      selection_pooling = "across_cohorts"
    )
    for (v in list(
      list(fml = ~1, icpt = TRUE, zero = FALSE, Q_sl = NULL),
      list(fml = ~0 + x, icpt = FALSE, zero = FALSE, Q_sl = cbind(x = sp_raw$x)),
      list(fml = ~0, icpt = FALSE, zero = TRUE, Q_sl = NULL)
    )) {
      fit_v <- do.call(att_fn, utils::modifyList(fx$args, list(selection_xformla = v$fml)))
      ref_v <- sp_ref_public_reference(
        sp_raw, est_method = em, selection_model = "independent_period_logit",
        selection_pooling = "across_cohorts", route = "F",
        sel_intercept = v$icpt, sel_zero_col = v$zero, Q_slope = v$Q_sl
      )
      testthat::expect_equal(fit_v$att, ref_v$att, tolerance = 1e-10)
      testthat::expect_equal(unname(fit_v$inffunc), unname(ref_v$inffunc), tolerance = 1e-7)
    }
  }

  # 3. Clustering (inference-only) and multiplier bootstrap on full F
  for (em in c("ipw", "reg")) {
    fx <- sp_public_fixture(
      est_method = em, route = "F",
      selection = "conditional_mtar", selection_model = "independent_period_logit",
      selection_pooling = "within_cohort"
    )
    fit_uncl <- do.call(att_fn, fx$args)
    fit_cl <- do.call(att_fn, utils::modifyList(fx$args, list(clustervars = "cl")))
    # Point estimates and unit IF rows are unchanged by clustering
    testthat::expect_identical(fit_cl$att, fit_uncl$att)
    testthat::expect_identical(fit_cl$inffunc, fit_uncl$inffunc)
    S_cl <- rowsum(fit_uncl$inffunc, fx$frame$cl, reorder = FALSE)
    S_cl <- sweep(S_cl, 2L, colMeans(S_cl), "-")
    V_cl <- crossprod(S_cl) / nrow(fx$frame)
    testthat::expect_equal(unname(fit_cl$V_analytical), unname(V_cl), tolerance = 1e-12)
    testthat::expect_equal(fit_cl$se, unname(sqrt(diag(V_cl) / nrow(fx$frame))), tolerance = 1e-12)

    # Multiplier bootstrap and simultaneous confidence bands
    set.seed(20261021L)
    fit_bs <- do.call(att_fn, utils::modifyList(fx$args, list(
      clustervars = "cl", bstrap = TRUE, cband = TRUE, biters = 64L
    )))
    testthat::expect_identical(fit_bs$att, fit_uncl$att)
    testthat::expect_true(all(is.finite(fit_bs$se) & fit_bs$se > 0))
    testthat::expect_true(is.finite(fit_bs$c) && fit_bs$c > 0)
    ag_bs <- agg_fn(fit_bs, type = "dynamic")
    testthat::expect_true(is.finite(ag_bs$overall.att) && is.finite(ag_bs$overall.se) && ag_bs$overall.se > 0)
  }
})

testthat::test_that("composed OR admission charges pooled selection cells outside the requested link supports", {
  ctrl <- sp_ns("ipw_owner_control")(
    max_chunk_rows = 23L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7
  )
  long <- sp_fx$long
  long$g[long$g == 0] <- Inf
  panel <- sp_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  frame <- sp_ns("construct_ipw_frame_owner")(panel, sp_fx$dat$id, sp_fx$V,
    weights = sp_fx$dat$w, control = ctrl, frame_g = sp_fx$dat$g, allow_frame_only = TRUE)
  reader <- sp_ns("construct_cp04_frame_reader")(frame, sp_fx$dat$id, sp_fx$Z,
    sp_fx$dat$g, sp_fx$dat$w, ctrl)
  views <- lapply(seq_along(c(3, 5)), function(j) {
    sp_ns("construct_ipw_link_view")(frame, 1, 2, g = c(3, 5)[j], K = Inf,
      control = ctrl, link_id = paste0("selection-budget-", j),
      functional_key = paste0("cmpcat:selection-budget:", j))
  })
  recipe <- sp_ns("or_pooled_outcome_recipe")(colnames(frame$Z_a), TRUE, c("x1", "x2"))
  receipts <- lapply(1:2, function(n) {
    backend <- sp_ns("fit_or_corollary1")(reader, views[seq_len(n)], ctrl,
      list(epsilon = 1e-13, maxit = 60L, trace = FALSE), rank_tol = 1e-10,
      outcome_recipe = recipe, selection_model = "pair_logit", selection_pooling = "across_cohorts")
    envelopes <- lapply(1:3, function(B) sp_ns("or_corollary1_if_envelope")(backend, B, 1L))
    testthat::expect_equal(envelopes[[3]]$scratch - envelopes[[2]]$scratch,
      envelopes[[2]]$scratch - envelopes[[1]]$scratch)
    testthat::expect_equal(envelopes[[3]]$work - envelopes[[2]]$work,
      envelopes[[2]]$work - envelopes[[1]]$work)
    admission <- sp_ns("cp04_preflight_links")(backend, 1L, chunk_rows = 23L)
    score_cells <- integer(0)
    sp_ns("cp04_replay_links")(backend, 1L, function(batch) {
      score_cells <<- c(score_cells, batch$diagnostics$n_selection_score_evaluations)
    }, chunk_rows = 23L)
    testthat::expect_true(all(score_cells == n + 1L))
    testthat::expect_false(backend$accounting$full_F_by_link_IF_retained)
    list(pre = envelopes[[1]], admission = admission)
  })
  # Same selected link and two outcome supports, with one extra pooled score cell.
  testthat::expect_identical(length(receipts[[1]]$pre$sids), length(receipts[[2]]$pre$sids))
  testthat::expect_gt(receipts[[2]]$pre$scratch, receipts[[1]]$pre$scratch)
  testthat::expect_gt(receipts[[2]]$pre$work, receipts[[1]]$pre$work)
  testthat::expect_gt(receipts[[2]]$admission$work_ops_bound, receipts[[1]]$admission$work_ops_bound)
})

testthat::test_that("selection IPW owner helper requires explicit GLM control and preserves an explicit fit", {
  fit_fn <- sp_ns("fit_cp04_selection_models")
  testthat::expect_identical(formals(fit_fn)["glm_control"], alist(glm_control = ))
  ctrl <- sp_ns("ipw_owner_control")(
    max_chunk_rows = 23L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7
  )
  gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
  long <- sp_fx$long
  long$g[long$g == 0] <- Inf
  panel <- sp_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  frame <- sp_ns("construct_ipw_frame_owner")(panel, sp_fx$dat$id, sp_fx$V,
    weights = sp_fx$dat$w, control = ctrl, frame_g = sp_fx$dat$g, allow_frame_only = TRUE)
  reader <- sp_ns("construct_cp04_frame_reader")(frame, sp_fx$dat$id, sp_fx$Z,
    sp_fx$dat$g, sp_fx$dat$w, ctrl)
  views <- lapply(sp_fx$links, function(lk) sp_ns("construct_ipw_link_view")(
    frame, lk$s, lk$t, g = lk$g, K = lk$K, control = ctrl,
    link_id = lk$link_id, functional_key = paste0("cmpcat:explicit-control:", lk$link_id)))
  error <- testthat::expect_error(fit_fn(reader, views, control = ctrl), class = "cdid_ipw_control_error")
  testthat::expect_identical(conditionMessage(error),
    "Supply an explicit `glm.control` list; no default solver policy is adopted.")
  testthat::expect_error(fit_fn(reader, views, control = ctrl, glm_control = NULL),
    class = "cdid_ipw_control_error")
  fit <- fit_fn(reader, views, control = ctrl, glm_control = gc,
    selection_model = "independent_period_logit", selection_pooling = "within_cohort")
  testthat::expect_identical(fit$glm_control, sp_ns("validate_ipw_glm_control")(gc))
  reference <- sp_ref_ipw_point(sp_fx, selection_model = "independent_period_logit",
    selection_pooling = "within_cohort", universe = "F")
  point <- vapply(seq_along(views), function(j) sp_ns("chained_full_link_estimate")(fit, j), numeric(1))
  testthat::expect_lt(max(abs(point - reference$mu)), 1e-11)
  batch <- sp_ns("evaluate_cp04_if_chunk")(fit, 1:23, seq_along(views), details = TRUE)
  testthat::expect_true(all(is.finite(batch$phi)))
  testthat::expect_identical(batch$normalization_count, reader$state$N)
})
