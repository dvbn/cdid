# CP08 / D78 (2026-10-24): the public estimated-work cap is OPT-IN, and the
# upstream work prefix handed to the baseline OR/DR fits is charged exactly once.
#
# Expected public policy is stated here independently of the production code:
#   * control = NULL, list() and list(counted_work_limit = Inf) carry NO user
#     work budget; the public profile shows Inf; memory/worker defaults are the
#     unchanged D27 values (1 GiB, one worker);
#   * an explicit positive finite exact whole number is a user cap that reaches
#     every private stage control unchanged and refuses with that exact bound;
#   * zero, negative, -Inf, NA/NaN, non-integral, out-of-exact-range, non-numeric,
#     dimensioned, classed or non-scalar values refuse;
#   * the finite private constructors keep refusing Inf; only the public seam
#     hands them the exact-count arithmetic ceiling 2^53 - 1, recorded apart from
#     the public Inf (never relabelled as the user's cap);
#   * the pipeline prefix P handed to the CP04 registry is counted ONCE in the
#     composed total: for a cumulative baseline OR/DR backend the total is
#     P + registry preparation + own work (its cumulative bound), for a backend
#     with another convention the unchanged legacy composition P + bound holds.
# No DGP draw, RNG or bootstrap. Fixture-backed blocks skip when the installed
# fixture directories are absent; the ordinary-mode block skips without did.

owc_small_panel <- function() {
  d <- expand.grid(id = 1:6, period = 1:4)
  d$G <- c(0, 0, 3, 3, 4, 4)[d$id]; d$x <- c(0.5, -1, 1, 0, 2, -2)[d$id]; d$y <- d$id + d$period
  d
}
owc_dims <- function() {
  prep <- chained_request_prepare_inputs(owc_small_panel(), NULL, "y", "period", "id", "G", ~x)
  list(prep = prep, dims = chained_public_dimensions(prep, biters = 50))
}
owc_ceiling <- 2^53 - 1
owc_stage_controls <- c("owner_control", "fitted_control", "point_control", "problem_control", "link_control", "certificate_control")

test_that("D78: NULL, an empty list and an explicit Inf carry no user work budget; memory/worker defaults are unchanged", {
  r0 <- chained_public_resolve_control(NULL)
  re <- chained_public_resolve_control(list())
  ri <- chained_public_resolve_control(list(counted_work_limit = Inf))
  for (r in list(r0, re, ri)) {
    expect_identical(r$profile$counted_work_limit, Inf)
    expect_identical(r$profile$numeric_workspace_bytes, 1073741824)      # D27 memory default unchanged
    expect_identical(r$profile$workers, 1)                               # D27 worker default unchanged
    expect_identical(r$profile_status, "ADOPTED_D78_OPTIONAL_WORK_CAP_D27_MEMORY_DEFAULTS")
    expect_null(r$chunk_rows)
  }
  # default versus explicit opt-out: same profile, different provenance
  expect_identical(r0$profile, ri$profile)
  expect_identical(r0$overrides, character(0L)); expect_identical(re$overrides, character(0L))
  expect_identical(ri$overrides, "counted_work_limit")
  expect_match(r0$source, "optional work cap disabled", fixed = TRUE); expect_match(re$source, "optional work cap disabled", fixed = TRUE)
  # explicit finite exact whole numbers are user caps and are kept exactly
  expect_identical(chained_public_resolve_control(list(counted_work_limit = 5e9))$profile$counted_work_limit, 5e9)
  expect_identical(chained_public_resolve_control(list(counted_work_limit = 7L))$profile$counted_work_limit, 7)
  expect_identical(chained_public_resolve_control(list(counted_work_limit = 1))$profile$counted_work_limit, 1)
  expect_identical(chained_public_resolve_control(list(counted_work_limit = 2^53 - 1))$profile$counted_work_limit, 2^53 - 1)
  expect_identical(chained_public_resolve_control(list(counted_work_limit = 1e11))$overrides, "counted_work_limit")
  # a finite cap does not touch the other defaults
  rf <- chained_public_resolve_control(list(counted_work_limit = 1e11))
  expect_identical(rf$profile$numeric_workspace_bytes, 1073741824); expect_identical(rf$profile$workers, 1)
  expect_identical(rf$covariance_strategy, r0$covariance_strategy); expect_identical(rf$glm, r0$glm)
})

test_that("D78: invalid work budgets refuse by shape, class, type and value", {
  refuse <- function(x) expect_error(chained_public_resolve_control(list(counted_work_limit = x)),
    class = "cdid_chained_public_control_error")
  refuse(0); refuse(-1); refuse(-1e9); refuse(-Inf)
  refuse(NA_real_); refuse(NaN); refuse(NA); refuse(NA_integer_)
  refuse(1.5); refuse(0.5); refuse(1e9 + 0.25)
  refuse(2^53); refuse(2^53 + 2); refuse(1e300)
  refuse("Inf"); refuse("1e11"); refuse(TRUE); refuse(NULL); refuse(list(1e9)); refuse(as.complex(1e9))
  refuse(c(1e9, 1e9)); refuse(numeric(0)); refuse(matrix(1e9)); refuse(array(1e9, c(1, 1, 1)))
  refuse(structure(1e9, class = "foo")); refuse(as.Date("2026-10-24")); refuse(factor("1e9"))
  # the independent guards are untouched
  expect_error(chained_public_resolve_control(list(bogus = 1)), class = "cdid_chained_public_control_error")
  expect_error(chained_public_resolve_control(list(numeric_workspace_bytes = Inf)), class = "cdid_chained_public_control_error")
  expect_error(chained_public_resolve_control(list(workers = Inf)), class = "cdid_chained_public_control_error")
  expect_error(chained_public_resolve_control(list(chunk_rows = Inf)), class = "cdid_chained_public_control_error")
  expect_error(chained_public_resolve_control(list(relative_accuracy = 1e-5)), class = "cdid_chained_public_accuracy_deferred")
})

test_that("D78 seam: a disabled cap hands every private stage control the 2^53-1 arithmetic ceiling and records it apart from the public Inf", {
  d <- owc_dims()
  for (ctl in list(NULL, list(), list(counted_work_limit = Inf))) {
    out <- chained_public_pipeline_control(d$dims, chained_public_resolve_control(ctl))
    rec <- out$record; ctrl <- out$control
    expect_identical(rec$counted_work_limit, Inf)                       # transparent public value
    expect_identical(rec$profile$counted_work_limit, Inf)
    expect_identical(rec$work_cap$status, "disabled")
    expect_identical(rec$work_cap$user_limit, Inf)
    expect_identical(rec$work_cap$effective_private_max_work_ops, owc_ceiling)
    expect_identical(rec$work_cap$internal_arithmetic_ceiling, owc_ceiling)
    expect_match(rec$work_cap$meaning, "not a resource budget")
    expect_identical(rec$admission$counted_work_limit, Inf)
    expect_identical(rec$admission$work_cap_status, "disabled")
    expect_identical(rec$admission$private_max_work_ops, owc_ceiling)
    for (nm in owc_stage_controls) expect_identical(ctrl[[nm]]$max_work_ops, owc_ceiling)
    cc <- ctrl$covariance_control
    expect_true(is.null(cc$max_work_ops) || identical(cc$max_work_ops, owc_ceiling))
    expect_identical(ctrl$catalog_bounds$max_work, owc_ceiling)
    expect_identical(ctrl$class_bounds, list(max_links = 36, max_classes = 36, max_class_dates = 4, max_E_entries = 144, max_mask_entries = 216))
    expect_match(rec$derived_bounds$bound[2L], "disabled")              # the ceiling row is labelled as such
    expect_identical(rec$derived_bounds$value[2L], owc_ceiling)
    expect_true(is.data.frame(rec$derived_bounds) && nrow(rec$derived_bounds) == 16L)
  }
  # default and explicit opt-out produce byte-identical private controls, with distinct provenance
  dflt <- chained_public_pipeline_control(d$dims, chained_public_resolve_control(NULL))
  optout <- chained_public_pipeline_control(d$dims, chained_public_resolve_control(list(counted_work_limit = Inf)))
  expect_identical(dflt$control, optout$control)
  expect_match(dflt$record$work_cap$source, "default"); expect_match(optout$record$work_cap$source, "opt-out")
})

test_that("D78 seam: an explicit finite cap reaches every private stage control unchanged and leaves every other control as under the disabled cap", {
  d <- owc_dims()
  cap <- 123456789
  fin <- chained_public_pipeline_control(d$dims, chained_public_resolve_control(list(counted_work_limit = cap)))
  dis <- chained_public_pipeline_control(d$dims, chained_public_resolve_control(NULL))
  expect_identical(fin$record$counted_work_limit, cap)
  expect_identical(fin$record$profile$counted_work_limit, cap)
  expect_identical(fin$record$work_cap$status, "enabled")
  expect_identical(fin$record$work_cap$user_limit, cap)
  expect_identical(fin$record$work_cap$effective_private_max_work_ops, cap)
  expect_identical(fin$record$work_cap$internal_arithmetic_ceiling, owc_ceiling)
  expect_identical(fin$record$admission$counted_work_limit, cap); expect_identical(fin$record$admission$work_cap_status, "enabled")
  for (nm in owc_stage_controls) expect_identical(fin$control[[nm]]$max_work_ops, cap)
  expect_identical(fin$control$catalog_bounds$max_work, cap)
  expect_match(fin$record$derived_bounds$bound[2L], "explicit"); expect_identical(fin$record$derived_bounds$value[2L], cap)
  # the small ceiling-independent part of the policy is untouched by the cap choice
  expect_identical(fin$record$chunk_rows, dis$record$chunk_rows)
  expect_identical(fin$record$numeric_workspace_entries, dis$record$numeric_workspace_entries)
  expect_identical(fin$control$owner_control[c("max_chunk_rows", "max_retained_entries", "max_matrix_entries", "max_output_entries")],
                   dis$control$owner_control[c("max_chunk_rows", "max_retained_entries", "max_matrix_entries", "max_output_entries")])
  expect_identical(fin$control$multiplier_budget, dis$control$multiplier_budget)
  expect_identical(fin$control$gmm_accuracy, dis$control$gmm_accuracy)
  expect_identical(fin$control$gmm_bounds, dis$control$gmm_bounds)
  expect_identical(fin$control$class_bounds, dis$control$class_bounds)
  expect_identical(fin$control$glm_control, dis$control$glm_control)
  # the exact ceiling is itself an admissible explicit cap; the record then says "enabled"
  edge <- chained_public_pipeline_control(d$dims, chained_public_resolve_control(list(counted_work_limit = 2^53 - 1)))
  expect_identical(edge$record$work_cap$status, "enabled"); expect_identical(edge$control$owner_control$max_work_ops, owc_ceiling)
})

test_that("the finite private constructors keep refusing Inf directly; the exact-sum helper treats Inf as a no-op bound", {
  expect_error(ipw_owner_control(max_chunk_rows = 4, max_retained_entries = 1e6, max_matrix_entries = 1e6,
    max_work_ops = Inf, max_output_entries = 1e6), class = "cdid_ipw_owner_control_error")
  expect_error(ipw_owner_int_bound(Inf, "max_work_ops"), class = "cdid_ipw_owner_bound_error")
  expect_error(ipw_owner_int_bound(2^53, "max_work_ops"), class = "cdid_ipw_owner_bound_error")
  expect_error(chained_public_count(Inf, "counted_work_limit"), class = "cdid_chained_public_dimension_error")
  ok <- ipw_owner_control(max_chunk_rows = 4, max_retained_entries = 1e6, max_matrix_entries = 1e6,
    max_work_ops = 2^53 - 1, max_output_entries = 1e6)
  expect_identical(ok$max_work_ops, 2^53 - 1)
  expect_error(ipw_owner_control(max_chunk_rows = 4, max_retained_entries = 1e6, max_matrix_entries = 1e6,
    max_work_ops = 2^53, max_output_entries = 1e6), class = "cdid_ipw_owner_bound_error")
  expect_identical(chained_public_work_arithmetic_ceiling(), 2^53 - 1)
  # the read-only OR-basis / reference-source consumers compare their charges against
  # the public profile value directly: Inf is a no-op bound, a finite bound still refuses
  expect_identical(ipw_owner_safe_sum(c(3, 4), Inf, "x"), 7)
  expect_error(ipw_owner_safe_sum(c(3, 4), 6, "x"), class = "cdid_ipw_owner_bound_error")
  expect_identical(ipw_owner_safe_sum(c(3, 4), 7, "x"), 7)
})

test_that("an explicit finite cap refuses the accepted pipeline at the owner preflight with the exact bound, before any nuisance fit", {
  d <- owc_dims(); prep <- d$prep
  capped <- chained_public_pipeline_control(d$dims, chained_public_resolve_control(list(counted_work_limit = 1024)))
  expect_identical(capped$record$counted_work_limit, 1024); expect_identical(capped$control$owner_control$max_work_ops, 1024)
  e <- tryCatch(run_chained_ipw_pipeline(prep$data_internal, "id", "period", "y", "G", capped$control, xformla = ~x,
    target_cohorts = prep$target_cohorts, control_group = "nevertreated", run_numerical_transfer = FALSE,
    root_witness = chained_public_automatic_witness, aggregations = character(0L), frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
    frame_weights = prep$frame_w, frame_g = prep$frame_g_arg, bstrap = FALSE, caller_rng = TRUE), error = function(e) e)
  expect_s3_class(e, "cdid_ipw_owner_bound_error")
  expect_match(conditionMessage(e), "exceeds explicit bound")
  expect_identical(e$diagnostics$bound, 1024)                            # the user's finite cap, no relaxation or ceiling substitution
})

test_that("registry work ownership: the upstream prefix is charged once for cumulative baseline backends; other conventions keep the legacy composition", {
  compose <- cp04_registry_work_accounting
  # hand-specified: pipeline prefix P = 1000; registry preparation 50 (prepared = 1050); backend own work 300
  cumulative <- function(cls, bound = 1350, own = 300, ext = 1050) structure(list(accounting = list(
    work_ops_bound = bound, own_work_ops = own, external_work_ops = ext, work_ops_convention = "cumulative")), class = cls)
  a <- compose(1000, 1050, cumulative("cdid_or_baseline"), 1e6)
  expect_identical(a$composed_total, 1350)                 # P + 50 + 300, not P + 1350
  expect_identical(a$prefix, 1000); expect_identical(a$prepared, 1050); expect_identical(a$registry_preparation, 50)
  expect_identical(a$backend_own_work_ops, 300); expect_identical(a$backend_convention, "cumulative")
  expect_identical(compose(1000, 1050, cumulative("cdid_dr_baseline"), 1e6)$composed_total, 1350)
  expect_identical(compose(0, 0, cumulative("cdid_or_baseline", bound = 300, own = 300, ext = 0), 1e6)$composed_total, 300)
  # Hand-independent components: P=1000, D=20, C=30, own=300.
  # The most-specific selection subclass wins over a cumulative OR parent.
  P <- 1000; D <- 20; C <- 30; own <- 300
  for (cls in list("cdid_cp04_corollary1",
      c("cdid_cp04_selection_backend", "cdid_cp04_corollary1"),
      c("cdid_or_selection_backend", "cdid_or_corollary1"), "cdid_dr_corollary1")) {
    inc <- structure(list(accounting = list(work_ops_bound = own)), class = cls)
    got <- compose(P, P + D + C, inc, P + D + C + own)
    expect_identical(got$backend_convention, "incremental")
    expect_identical(got$registry_preparation, D + C)
    expect_identical(got$composed_total, 1350)
    expect_error(compose(P, P + D + C, inc, 1349), class = "cdid_ipw_owner_bound_error")
  }
  old_or <- structure(list(accounting = list(work_ops_bound = 1350, external_work_ops = 1050)),
    class = "cdid_or_corollary1")
  expect_identical(compose(1000,1050,old_or,1350)$composed_total,1350)
  expect_identical(compose(1000,1050,old_or,1350)$backend_own_work_ops,300)
  expect_identical(compose(1000,1050,old_or,1350)$backend_convention,"cumulative")
  expect_error(compose(1000,1050,old_or,1349),class="cdid_ipw_owner_bound_error")
  # Truly unresolved wrappers keep the old composition; no convention guessed.
  incremental <- structure(list(accounting = list(work_ops_bound = 300)), class = "cdid_unresolved_wrapper")
  b <- compose(1000,1050,incremental,1e6)
  expect_identical(b$composed_total,1300);expect_identical(b$backend_convention,"legacy_unlabelled")
  # no backend: the prepared prefix
  expect_identical(compose(1000, 1050, NULL, 1e6)$composed_total, 1050)
  # a finite cap is never bypassed
  expect_error(compose(1000, 1050, cumulative("cdid_or_baseline"), 1349), class = "cdid_ipw_owner_bound_error")
  expect_identical(compose(1000, 1050, cumulative("cdid_or_baseline"), 1350)$composed_total, 1350)
  expect_error(compose(1000, 1050, incremental, 1299), class = "cdid_ipw_owner_bound_error")
  expect_error(compose(1000, 1050, NULL, 1049), class = "cdid_ipw_owner_bound_error")
  # inconsistent declarations refuse instead of guessing; no negative increments
  expect_error(compose(1000, 1050, cumulative("cdid_or_baseline", ext = 1000), 1e6), class = "cdid_cp04_integrity_error")   # saw another prefix
  expect_error(compose(1000, 1050, cumulative("cdid_or_baseline", bound = 1400), 1e6), class = "cdid_cp04_integrity_error") # own + ext != bound
  expect_error(compose(1000, 1050, cumulative("cdid_cp04_corollary1"), 1e6), class = "cdid_cp04_integrity_error")        # only baseline OR/DR are cumulative
  expect_error(compose(1000, 1050, cumulative("cdid_or_baseline", own = -1, bound = 1049), 1e6), class = "cdid_ipw_owner_bound_error")
  expect_error(compose(1000, 999, cumulative("cdid_or_baseline"), 1e6), class = "cdid_cp04_integrity_error")               # prepared below prefix
  expect_error(compose(Inf, 1050, incremental, 1e6), class = "cdid_ipw_owner_bound_error")
  # a legacy backend without a recorded bound keeps the former `c(prefix, NULL)` composition, labelled
  none <- compose(1000, 1050, structure(list(accounting = list()), class = "cdid_or_baseline"), 1e6)
  expect_identical(none$composed_total, 1000); expect_identical(none$backend_convention, "legacy_no_recorded_work_ops_bound")
})

owc_or_fixture <- function() {
  dir <- testthat::test_path("fixtures", "cp06_pooled_or")
  if (!dir.exists(dir) || !file.exists(file.path(dir, "load_fixture.R"))) return(NULL)
  loader <- new.env(parent = globalenv())
  sys.source(file.path(dir, "load_fixture.R"), envir = loader)
  list(dir = dir, loader = loader)
}

test_that("baseline OR fit: own work is prefix-independent and the cumulative bound equals the prefix plus own work (counted once)", {
  fxl <- owc_or_fixture()
  skip_if(is.null(fxl), "pooled OR fixture directory is not installed")
  order <- c(seq.int(108L, 2L, -2L), seq.int(107L, 1L, -2L))
  fx <- fxl$loader$load_pooled_or_reference(fxl$dir, "NYT_baseline", order)
  ctrl <- ipw_owner_control(max_chunk_rows = 23L, max_retained_entries = 1e7, max_matrix_entries = 2e7,
    max_work_ops = 1e10, max_output_entries = 1e7)
  recipe <- or_pooled_outcome_recipe(colnames(fx$frame$Z_a), fx$outcome_recipe$intercept, fx$outcome_recipe$shared_columns)
  W <- 123456
  f0 <- fit_or_baseline(fx$reader, fx$views, ctrl, rank_tol = 1e-10, outcome_recipe = recipe)
  fW <- fit_or_baseline(fx$reader, fx$views, ctrl, rank_tol = 1e-10, outcome_recipe = recipe, external_work_ops = W)
  for (f in list(f0, fW)) {
    expect_identical(f$accounting$work_ops_convention, "cumulative")
    expect_true(is.numeric(f$accounting$own_work_ops) && f$accounting$own_work_ops > 0)
  }
  expect_identical(f0$accounting$external_work_ops, 0); expect_identical(fW$accounting$external_work_ops, W)
  expect_identical(fW$accounting$own_work_ops, f0$accounting$own_work_ops)          # own work does not depend on the prefix
  expect_identical(f0$accounting$work_ops_bound, f0$accounting$own_work_ops)        # zero prefix: cumulative == own
  expect_identical(fW$accounting$work_ops_bound, f0$accounting$own_work_ops + W)    # the prefix appears exactly once
  # estimates are untouched by the prefix
  expect_identical(fW$shared$beta, f0$shared$beta)
  est <- function(f) vapply(seq_along(f$links), function(j) chained_full_link_estimate(f, j), numeric(1))
  expect_identical(est(fW), est(f0))
  # the registry composition on this real backend charges W once; the former
  # pipeline composition (W + cumulative bound) would have charged it twice
  wa <- cp04_registry_work_accounting(W, W, fW, ctrl$max_work_ops)
  expect_identical(wa$composed_total, fW$accounting$work_ops_bound)
  expect_identical(wa$composed_total, W + f0$accounting$own_work_ops)
  expect_identical((W + fW$accounting$work_ops_bound) - wa$composed_total, W)
  # a cap between the single and the double charge admits the fit under the corrected composition only
  expect_identical(cp04_registry_work_accounting(W, W, fW, wa$composed_total)$composed_total, wa$composed_total)
  expect_error(ipw_owner_safe_sum(c(W, fW$accounting$work_ops_bound), wa$composed_total, "former double composition"),
    class = "cdid_ipw_owner_bound_error")
})

owc_dr_fixture <- function() {
  dir <- testthat::test_path("fixtures", "dr_finite_v1")
  if (!dir.exists(dir) || !file.exists(file.path(dir, "load_fixture.R"))) return(NULL)
  loader <- new.env(parent = globalenv())
  sys.source(file.path(dir, "load_fixture.R"), envir = loader)
  list(dir = dir, loader = loader)
}

test_that("baseline DR fit: own work is prefix-independent and the cumulative bound equals the prefix plus own work (counted once)", {
  fxl <- owc_dr_fixture()
  skip_if(is.null(fxl), "DR finite fixture directory is not installed")
  fx <- fxl$loader$load_dr_finite_reference(fxl$dir, fixture_id = "F1_baseline_none_intercept_slope", order = NULL, max_chunk_rows = 13L)
  ctrl <- ipw_owner_control(max_chunk_rows = fx$control$max_chunk_rows, max_retained_entries = fx$control$max_retained_entries,
    max_matrix_entries = fx$control$max_matrix_entries, max_work_ops = 1e12, max_output_entries = fx$control$max_output_entries)
  W <- 123456
  fit <- function(w) fit_dr_baseline(fx$reader, fx$views, control = ctrl, glm_control = fx$glm_control, rank_tol = fx$rank_tol,
    outcome_recipe = fx$outcome_recipe, treatment_columns = fx$treatment_columns, external_work_ops = w)
  f0 <- fit(0); fW <- fit(W)
  for (f in list(f0, fW)) {
    expect_identical(f$accounting$work_ops_convention, "cumulative")
    expect_true(is.numeric(f$accounting$own_work_ops) && f$accounting$own_work_ops > 0)
  }
  expect_identical(f0$accounting$external_work_ops, 0); expect_identical(fW$accounting$external_work_ops, W)
  expect_identical(fW$accounting$own_work_ops, f0$accounting$own_work_ops)
  expect_identical(f0$accounting$work_ops_bound, f0$accounting$own_work_ops)
  expect_identical(fW$accounting$work_ops_bound, f0$accounting$own_work_ops + W)
  expect_identical(fW$shared$beta, f0$shared$beta)
  est <- function(f) vapply(seq_along(f$links), function(j) chained_full_link_estimate(f, j), numeric(1))
  expect_identical(est(fW), est(f0))
  wa <- cp04_registry_work_accounting(W, W, fW, ctrl$max_work_ops)
  expect_identical(wa$composed_total, fW$accounting$work_ops_bound)
  expect_identical(wa$composed_total, W + f0$accounting$own_work_ops)
  # a different-convention wrapper of the same family is untouched: the Corollary 1
  # DR route records own work only and composes as prefix + bound
  legacy <- structure(list(accounting = list(work_ops_bound = f0$accounting$own_work_ops)), class = "cdid_dr_corollary1")
  expect_identical(cp04_registry_work_accounting(W, W, legacy, ctrl$max_work_ops)$composed_total, W + f0$accounting$own_work_ops)
})

test_that("ordinary did mode delegates before evaluating the chained control", {
  skip_if_not_installed("did")
  env <- new.env()
  utils::data("mpdta", package = "did", envir = env)
  res <- att_gt(yname = "lemp", tname = "year", idname = "countyreal", gname = "first.treat", data = env$mpdta,
    mode = "ordinary", bstrap = FALSE, cband = FALSE,
    control = stop("the chained `control` promise must not be forced in ordinary mode"))
  expect_s3_class(res, "MP")
  expect_false(inherits(res, "cdid_MP"))
})
