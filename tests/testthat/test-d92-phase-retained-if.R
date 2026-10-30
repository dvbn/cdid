# SPDX-License-Identifier: GPL-3.0-only
#
# D92 phase-local CP04 retained moment IF: when the accepted whole-lifetime
# store (`cp04_retained_if_reserve`) declines for a one-block, non-transfer,
# dated augmented backend (`cdid_augmented_backend` with `reference_ref` and
# `layout`), independently consider a covariance-only store at the original
# covariance partition (`cov_plan$chunk_rows`) and a GMM-only store at the
# original GMM partition (`point_preflight$chunk_rows`). Each phase fills on
# its first evaluation and slices on its second evaluation at its own original
# partition. The covariance store and every closure/reference are released
# before allocating GMM consumers; the GMM store is released before target IF
# transport, auxiliary reporting and final diagnostic copying. Both stores are
# never live simultaneously.
#
# STATUS: Test bodies UNEXECUTED in this adaptation job; only static Base R
# parse/AST checks ran. The root sequences the independent numerical gates.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures (test-local copies; testthat files do not share definitions)
# ---------------------------------------------------------------------------

d92_ctrl <- function(chunk = 7L, work = 1e10, matrix = 2e7, retained = matrix, output = matrix) {
  api$ipw_owner_control(
    max_chunk_rows = chunk,
    max_retained_entries = retained,
    max_matrix_entries = matrix,
    max_work_ops = work,
    max_output_entries = output
  )
}
d92_glm <- list(epsilon = 1e-12, maxit = 50L, trace = FALSE)

# Unbalanced panel with 2 treated cohorts (3, 5) and never-treated (0 / Inf),
# time-varying covariate `x` and time-varying positive sampling weight `w`
# supplied at every unit's first calendar period (period == 1).
d92_dated_panel <- function(N = 48L, T = 6L, seed = 43L, keep_reference = FALSE) {
  set.seed(seed)
  g <- rep(c(3, 5, 0), length.out = N)
  xi <- stats::rnorm(N)
  wi <- 0.8 + 0.4 * stats::runif(N)
  d <- expand.grid(period = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "period")]
  u <- d$id
  d$G <- g[u]
  d$x <- xi[u] + 0.15 * cos(0.5 * u + d$period)
  d$w <- wi[u] * (1 + 0.05 * (d$period - 1L))
  eff <- ifelse(d$G > 0 & d$period >= d$G, 0.9 + 0.2 * (d$period - d$G), 0)
  d$y <- 0.25 * u + 0.3 * d$period + 0.4 * d$x + eff + stats::rnorm(nrow(d))
  d$cl <- (u - 1L) %/% 4L + 1L
  # The original default drops both dates and witnesses missing reference X.
  # Keeping the already generated period-4 X/w (with y still missing) gets to
  # the independent-date logit's correct finite-root refusal on this panel.
  if (keep_reference) d$y[d$period == 4L & (d$id %% 7L == 0L)] <- NA_real_
  drop_rows <- which(d$period %in% (if (keep_reference) 6L else c(4L, 6L)) & (d$id %% 7L == 0L))
  if (length(drop_rows)) d <- d[-drop_rows, , drop = FALSE]
  rownames(d) <- NULL
  d
}


# Positive cache fixture: exact test-local copy of d92_route_panel from
# test-d92-shared-if-reuse.R (producer SHA256
# 76e01846f15f4c469d4472dc65971b011e18b8f889d58e6702df80caa1dc60bd).
# Reuse its exposed arguments verified in DATED_FIXTURE_FINDING.md; only the
# time-column name changes. This supplies finite-root preconditions for cache
# lifecycle checks, not an independent statistical or numerical target.
d92_route_panel <- function(n = 36L, dates = c(1, 2, 4, 5), cohorts = c(4, 5, 0),
                            vary_x = FALSE, vary_w = FALSE, drop_y = TRUE) {
  id <- seq_len(n)
  g <- rep(cohorts, length.out = n)
  d <- expand.grid(time = dates, id = id, KEEP.OUT.ATTRS = FALSE)[, c("id", "time")]
  ti <- match(d$time, dates)
  d$G <- g[d$id]
  bx <- round(sin(1.7 * id), 2)
  d$x <- if (vary_x) round(bx[d$id] + 0.25 * ti + 0.3 * cos(2.3 * d$id + 1.1 * ti), 2) else bx[d$id]
  bw <- round(0.6 + (id %% 7L) / 7, 3)
  d$w <- if (vary_w) round(bw[d$id] * (1 + 0.15 * ti * ((d$id %% 2L) == 0L)), 3) else bw[d$id]
  te <- ifelse(d$G > 0 & d$time >= d$G, 1 + 0.5 * (d$time - d$G), 0)
  d$y <- round(cos(d$id) + 0.4 * d$time + 0.8 * d$x + te + 0.5 * sin(5.1 * d$id + 2.7 * ti), 3)
  if (drop_y) d$y[((d$id * 37L + ti * 11L) %% 13L) < 3L] <- NA_real_
  d$rho <- round(0.6 + 0.3 * (d$id %% 4L) / 4, 3)
  d$cl <- (d$id - 1L) %/% 3L + 1L
  d
}

d92_positive_dated_panel <- function() {
  d <- d92_route_panel(n = 48L, dates = as.double(1:6), cohorts = c(3, 5, 0),
                       vary_x = TRUE, vary_w = TRUE, drop_y = TRUE)
  names(d)[names(d) == "time"] <- "period"
  d
}

d92_pairs <- function(T = 6L) {
  pr <- list()
  for (s in seq_len(T - 1L)) for (t in seq.int(s + 1L, T)) pr[[length(pr) + 1L]] <- c(s, t)
  pr
}

# Direct sealed owners, reference reader and dated augmented backend.
# Finite private check budget admits the legacy planner before the one-block
# phase plan is considered; the old 1e10 refusal is tested explicitly below.
d92_dated_fixture <- function(d = d92_positive_dated_panel(), chunk = 7L, ctrl = d92_ctrl(chunk, work = 1e11),
                              xf = ~x, sf = ~x, fix_weights = "first_period",
                              imputation_pooling = "across_cohorts",
                              selection_model = "independent_period_logit",
                              selection_pooling = "within_cohort",
                              automatic_chunk_rows = TRUE) {
  u <- d[!duplicated(d$id), c("id", "G")]
  u <- u[order(u$id), , drop = FALSE]
  N <- nrow(u)
  fg <- ifelse(u$G == 0, Inf, as.double(u$G))
  d_int <- d
  d_int$G <- ifelse(d_int$G == 0, Inf, as.double(d_int$G))
  po <- api$construct_ipw_panel_owner(d_int, "id", "period", "y", "G", ctrl)
  cal <- as.double(po$all_times)
  pos_unit <- match(d_int$id, u$id)
  first <- match(u$id, d_int$id)
  in_panel <- seq_along(u$id)
  wres <- api$chained_reference_weights(
    d_int, NULL, "w", d_int$id, d_int$period, cal, u$id,
    in_panel, integer(0), first, pos_unit, frame_supplied = FALSE
  )
  vars <- unique(c(all.vars(xf), all.vars(sf)))
  prov <- api$chained_reference_classify(
    vars, d_int, NULL, FALSE, d_int$id, d_int$period, "period", u$id, pos_unit, "xformla"
  )
  api$chained_reference_check_classes(prov, "xformla")
  formulas <- list(
    x = list(formula = xf, dated = TRUE, static = NULL),
    s = list(formula = sf, dated = TRUE, static = NULL)
  )
  src <- api$chained_reference_source(
    formulas, prov, d_int, NULL, FALSE, d_int$id, d_int$period, "period",
    u$id, cal, 0, wres, api$chained_reference_fix_weights(fix_weights),
    api$chained_reference_route("augmented", "F"), resource = NULL
  )
  Z_car <- matrix(1, N, 1L, dimnames = list(NULL, "carrier_intercept"))
  V_car <- matrix(1, N, 1L, dimnames = list(NULL, "reader_intercept"))
  wt <- as.double(src$V)
  fo <- api$construct_ipw_frame_owner(
    po, u$id, Z_car, weights = wt, control = ctrl, frame_g = fg,
    allow_frame_only = FALSE, design_on_A_only = TRUE
  )
  reader <- api$construct_cp04_frame_reader(
    fo, u$id, V_car, fg, wt, ctrl, U = rep(TRUE, N), rho = rep(1, N),
    universe = "F", allow_unknown_off_A = FALSE, reference_carrier = TRUE
  )
  prs <- d92_pairs(length(cal))
  views <- list()
  for (g_target in c(3, 5)) {
    K_ctrl <- if (g_target == 3) c(5, Inf) else Inf
    for (k in seq_along(prs)) {
      p <- prs[[k]]
      idx <- length(views) + 1L
      views[[idx]] <- api$construct_ipw_link_view(
        fo, cal[p[1L]], cal[p[2L]], g = g_target, K = K_ctrl, control = ctrl,
        link_id = paste0("L_g", g_target, "_", p[1L], "_", p[2L]),
        functional_key = paste0("cmpcat:d92:g", g_target, ":", p[1L], "-", p[2L])
      )
    }
  }
  ref <- api$build_chained_reference_reader(src, reader, views, length(views), ctrl, retained = 0, work = 0)
  rr <- api$chained_reference_outcome_recipe(ref)
  fit <- api$fit_augmented_backend(
    reader, views, ctrl, d92_glm, 1e-10,
    outcome_recipe = rr$recipe, treatment_columns = rr$columns,
    est_method = "dr", imputation_pooling = imputation_pooling,
    selection_model = selection_model, selection_pooling = selection_pooling,
    retained_external_entries = 0, external_work_ops = 0,
    automatic_chunk_rows = automatic_chunk_rows, reference = ref
  )
  list(d = d, d_int = d_int, po = po, fo = fo, reader = reader, views = views,
       src = src, ref = ref, rr = rr, fit = fit, ctrl = ctrl, chunk = chunk)
}

# Evaluate complete original-unit IF at a specific row partition B.
# Never compare across different row partitions on a recipe-cohort layout.
d92_phi_at_partition <- function(fit, links = seq_along(fit$links), B = fit$control$max_chunk_rows, scale = "core") {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  # The producer carries unit IDs separately and uses NULL row labels;
  # column labels follow the independently fitted link declaration order.
  out <- matrix(0, st$N, length(links),
    dimnames = list(NULL, vapply(fit$links[links], `[[`, character(1L), "link_id")))
  for (s in seq.int(1L, st$N, by = B)) {
    rr <- seq.int(s, min(st$N, s + B - 1L))
    out[rr, ] <- api$evaluate_augmented_if_chunk(fit, rr, links, scale = scale, details = FALSE)$phi
  }
  out
}

d92_prepare_dated_augmented <- function(d, xf = ~x, sf = NULL, weightsname = "w",
                                        fix_weights = "first_period",
                                        imputation_pooling = "across_cohorts",
                                        selection_model = "independent_period_logit",
                                        selection_pooling = "within_cohort",
                                        clustered = FALSE, profile = NULL) {
  prep <- api$chained_request_prepare_inputs(
    data = d, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = NULL, weightsname = NULL, clustervars = if (clustered) "cl" else NULL,
    anticipation = 0, frame = NULL, outcome_regression = FALSE,
    doubly_robust = TRUE, or_resource = profile
  )
  sf_eff <- if (is.null(sf)) xf else sf
  pos_unit <- match(d$id, prep$frame_ids)
  first <- match(prep$frame_ids, d$id)
  in_panel <- seq_along(prep$frame_ids)
  wres <- api$chained_reference_weights(
    d, NULL, weightsname, d$id, d$period, prep$calendar, prep$frame_ids,
    in_panel, integer(0), first, pos_unit, frame_supplied = FALSE
  )
  vars <- unique(c(all.vars(xf), all.vars(sf_eff)))
  prov <- api$chained_reference_classify(
    vars, d, NULL, FALSE, d$id, d$period, "period", prep$frame_ids, pos_unit, "xformla"
  )
  api$chained_reference_check_classes(prov, "xformla")
  formulas <- list(
    x = list(formula = xf, dated = TRUE, static = NULL),
    s = list(formula = sf_eff, dated = TRUE, static = NULL)
  )
  src <- api$chained_reference_source(
    formulas, prov, d, NULL, FALSE, d$id, d$period, "period",
    prep$frame_ids, prep$calendar, 0, wres,
    api$chained_reference_fix_weights(fix_weights),
    api$chained_reference_route("augmented", "F"),
    resource = profile
  )
  sel_norm <- api$chained_request_augmented_declaration(
    "conditional_mtar", selection_model, selection_pooling, "F", NULL, NULL, FALSE
  )
  cp04 <- api$chained_request_reference_selection_record(
    "augmented", xf, sf_eff, sel_norm, NULL, FALSE, prep$n_frame,
    imputation_pooling, selection_formula_inherited = is.null(sf)
  )
  cp04$reference_source <- src
  prep$reference_source <- src
  prep$cp04 <- cp04
  prep$frame_Z <- matrix(1, prep$n_frame, 1L, dimnames = list(NULL, "carrier_intercept"))
  if (!is.null(src$V)) prep$frame_w <- as.double(src$V)
  if (is.null(prep$frame_g_arg)) prep$frame_g_arg <- prep$frame_g
  prep
}

d92_run_dated_augmented <- function(d = d92_positive_dated_panel(), chunk = NULL, reduce = TRUE,
                                    weighting = "optimal", control_group = "notyettreated",
                                    base_period = "universal",
                                    aggregations = c("simple", "group", "dynamic", "calendar"),
                                    owner_caps = NULL, clustered = FALSE, ns = api) {
  ctrl_list <- if (is.null(chunk)) NULL else list(chunk_rows = chunk)
  resolved <- api$chained_public_resolve_control(ctrl_list, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  prep <- d92_prepare_dated_augmented(d, clustered = clustered, profile = resolved$profile)
  prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = 99L, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  if (!is.null(owner_caps)) for (nm in names(owner_caps)) built$control$owner_control[[nm]] <- owner_caps[[nm]]
  res <- ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = built$control, xformla = ~1, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = prep$target_cohorts, anticipation = 0L, control_group = control_group,
    base_period = base_period, weighting = weighting, reduce_structural_redundancy = reduce,
    run_numerical_transfer = FALSE, root_witness = NULL, aggregations = aggregations,
    cband = FALSE, alp = 0.05, biters = 99L, seed = 101L,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w,
    frame_g = prep$frame_g_arg, allow_frame_only = prep$n_frame_only > 0L,
    na.rm = FALSE, bstrap = FALSE, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = "post", pt_start = NULL, cp04 = prep$cp04
  )
  res$nuisance_diagnostics <- ns$cp04_public_diagnostics(res, built$control)
  res
}

# Test-local installed functions issue distinct original partitions. This
# instruments the planners, not the IF producer or its arithmetic; the separate
# capacity tests below exercise the real planner and adjacent cap refusals.
d92_distinct_partition_ns <- function(B_cov = 9L, B_gmm = 7L, base_ns = api) {
  e <- new.env(parent = base_ns)
  for (nm in c("att_gt", "chained_public_att_gt", "cp04_public_diagnostics",
               "run_chained_ipw_pipeline", "cp04_stream_plan", "cp04_stream_plan_legacy",
               "cp04_stream_plan_full_block", "cp04_preflight_links",
               "cp04_retained_if_reserve", "cp04_retained_if_phase_reserve",
               "cp04_stream_covariance", "cp04_replay_links")) {
    f <- get(nm, base_ns); environment(f) <- e; assign(nm, f, e)
  }
  e$cp04_stream_plan_legacy <- function(backend, links, retained = 0) {
    res <- base_ns$cp04_stream_plan_legacy(backend, links, retained)
    res$chunk_rows <- min(res$chunk_rows, as.integer(B_cov))
    res
  }
  e$cp04_preflight_links <- function(backend, links, retained = 0, chunk_rows = NULL) {
    cr <- if (is.null(chunk_rows)) as.integer(B_gmm) else chunk_rows
    base_ns$cp04_preflight_links(backend, links, retained = retained, chunk_rows = cr)
  }
  e
}

d92_stream_ns <- function(base_ns = api) {
  e <- new.env(parent = base_ns)
  for (nm in names(as.list(base_ns, all.names = TRUE))) {
    if (exists(nm, base_ns, inherits = FALSE)) assign(nm, get(nm, base_ns, inherits = FALSE), e)
  }
  f <- get("run_chained_ipw_pipeline", base_ns); environment(f) <- e; assign("run_chained_ipw_pipeline", f, e)
  e$cp04_retained_if_reserve <- function(...) NULL
  e$cp04_retained_if_phase_reserve <- function(...) NULL
  e
}

d92_phase_filter_ns <- function(allow_phases = c("covariance", "gmm"), base_ns = api) {
  e <- new.env(parent = base_ns)
  for (nm in names(as.list(base_ns, all.names = TRUE))) {
    if (exists(nm, base_ns, inherits = FALSE)) assign(nm, get(nm, base_ns, inherits = FALSE), e)
  }
  f <- get("run_chained_ipw_pipeline", base_ns); environment(f) <- e; assign("run_chained_ipw_pipeline", f, e)
  e$cp04_retained_if_reserve <- function(...) NULL
  e$cp04_retained_if_phase_reserve <- function(backend, links, phase, ...) {
    if (!phase %in% allow_phases) return(NULL)
    base_ns$cp04_retained_if_phase_reserve(backend, links, phase = phase, ...)
  }
  e
}

d92_find_store <- function(x, seen = new.env(parent = emptyenv())) {
  if (inherits(x, "cdid_cp04_retained_if_bound")) return(TRUE)
  if (is.environment(x)) {
    if (identical(x$kind, "cp04_retained_if")) return(TRUE)
    key <- format(x)
    if (exists(key, seen, inherits = FALSE)) return(FALSE)
    assign(key, TRUE, seen)
    for (nm in ls(x, all.names = TRUE)) if (d92_find_store(get(nm, x, inherits = FALSE), seen)) return(TRUE)
    return(FALSE)
  }
  if (is.list(x)) for (el in x) if (d92_find_store(el, seen)) return(TRUE)
  FALSE
}

d92_equal <- function(a, b) {
  expect_false(d92_find_store(a))
  expect_false(d92_find_store(b))
  expect_true(is.matrix(a$reporting$if_full))
  expect_true(is.matrix(a$covariance_stream$V))
  for (nm in c("mhat", "V", "mean_phi"))
    expect_equal(a$covariance_stream[[nm]], b$covariance_stream[[nm]], tolerance = 1e-12)
  expect_equal(a$gmm_fit$targets, b$gmm_fit$targets, tolerance = 1e-12)
  for (nm in c("sandwich_cov", "materialized_target_if"))
    expect_equal(a$gmm_fit$if_projection[[nm]], b$gmm_fit$if_projection[[nm]], tolerance = 1e-12)
  expect_equal(a$covariance_handling, b$covariance_handling, tolerance = 1e-12)
  expect_equal(a$full_frame_transport$Phi_full, b$full_frame_transport$Phi_full, tolerance = 1e-12)
  for (nm in c("estimates", "se", "analytical_se", "if_full", "vcov"))
    expect_equal(a$reporting[[nm]], b$reporting[[nm]], tolerance = 1e-12)
  expect_identical(a$reporting$target_meta, b$reporting$target_meta)
  expect_equal(a$auxiliary_reporting, b$auxiliary_reporting, tolerance = 1e-12)
  expect_equal(a$nuisance_diagnostics, b$nuisance_diagnostics, tolerance = 1e-12)
  expect_equal(a$selection, b$selection, tolerance = 1e-12)
  pa <- api$chained_pretest_record(a$reporting)
  pb <- api$chained_pretest_record(b$reporting)
  expect_identical(pa$status, pb$status)
  expect_equal(pa, pb, tolerance = 1e-12)
  expect_identical(sort(names(a$aggregations)), sort(names(b$aggregations)))
  for (ag in names(a$aggregations)) for (nm in c("estimates", "se", "if_full", "overall"))
    expect_equal(a$aggregations[[ag]][[nm]], b$aggregations[[ag]][[nm]], tolerance = 1e-12)
}

# ---------------------------------------------------------------------------
# 1. Complete original-unit IF, covariance, SE and public outputs versus forced
#    streaming at 1e-12; both original callback partitions and order (B_cov != B_gmm).
# ---------------------------------------------------------------------------

test_that("the original dated panel refuses missing reference X and separated date logits", {
  missing_ref <- tryCatch(d92_dated_fixture(d = d92_dated_panel()), error = identity)
  expect_s3_class(missing_ref, "cdid_reference_missing_error")
  expect_match(conditionMessage(missing_ref), "date 4, variable `x`", fixed = TRUE)

  separated <- tryCatch(d92_dated_fixture(d = d92_dated_panel(keep_reference = TRUE)), error = identity)
  expect_s3_class(separated, "cdid_cp04_separation_error")
  expect_match(conditionMessage(separated), "Logit has no finite root", fixed = TRUE)
})

test_that("dated augmented phase-local stores match forced streaming at 1e-12 across distinct covariance and GMM partitions", {
  phase_ns <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L)
  fx <- d92_dated_fixture(chunk = 11L)
  b <- fx$fit
  expect_s3_class(b, "cdid_augmented_backend")
  expect_false(is.null(b$reference_ref))
  expect_false(is.null(b$layout))
  links <- seq_along(b$links)
  L <- length(links)
  N <- get("cp04_reader", b$reader_ref, inherits = FALSE)$state$N

  # Covariance partition B_cov = 9, GMM partition B_gmm = 7.
  B_cov <- 9L
  B_gmm <- 7L
  cov_plan <- phase_ns$cp04_stream_plan(b, links, retained = 0)
  expect_identical(length(cov_plan$blocks), 1L)

  # Whole-lifetime store declines when covariance and GMM row partitions differ.
  expect_null(api$cp04_retained_if_reserve(b, links, B_gmm, cov_plan, 0, 0, 1000, 1000, fx$ctrl))

  # Phase 1: covariance store fills at B_cov on mean pass and slices at B_cov on Gram pass.
  store_cov <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  expect_true(is.environment(store_cov))
  expect_identical(store_cov$phase, "covariance")
  expect_identical(store_cov$max_slice_passes, 1L)
  expect_identical(store_cov$B, B_cov)
  expect_identical(store_cov$state, "reserved")

  cov_batches <- list()
  prod_cov <- api$cp04_retained_if_bind(b, store_cov)
  api$cp04_replay_links(prod_cov, links, function(z) {
    cov_batches[[length(cov_batches) + 1L]] <<- list(rows = z$rows, id = z$id, links = z$links, phi = z$phi)
  }, scale = "core", chunk_rows = B_cov)
  expect_identical(store_cov$state, "filled")
  expect_identical(store_cov$producer_passes, 1L)
  expect_identical(store_cov$slice_passes, 0L)
  expect_identical(length(cov_batches), as.integer(ceiling(N / B_cov)))
  expect_identical(store_cov$id, sort(unique(fx$d$id)))
  expect_equal(store_cov$phi, d92_phi_at_partition(b, links, B = B_cov, scale = "core"), tolerance = 1e-12)
  cov_original <- list()
  api$cp04_replay_links(b, links, function(z) {
    cov_original[[length(cov_original) + 1L]] <<- list(rows = z$rows, id = z$id, links = z$links, phi = z$phi)
  }, scale = "core", chunk_rows = B_cov)
  expect_equal(cov_batches, cov_original, tolerance = 1e-12)
  expect_identical(unlist(lapply(cov_batches, `[[`, "rows"), use.names = FALSE), seq_len(N))

  cov_slices <- list()
  api$cp04_replay_links(prod_cov, links, function(z) {
    cov_slices[[length(cov_slices) + 1L]] <<- list(rows = z$rows, id = z$id, links = z$links, phi = z$phi)
  }, scale = "core", chunk_rows = B_cov)
  expect_identical(store_cov$producer_passes, 1L)
  expect_identical(store_cov$slice_passes, 1L)
  expect_identical(length(cov_slices), length(cov_batches))
  for (k in seq_along(cov_batches)) {
    expect_identical(cov_slices[[k]]$rows, cov_batches[[k]]$rows)
    expect_identical(cov_slices[[k]]$id, cov_batches[[k]]$id)
    expect_identical(cov_slices[[k]]$links, cov_batches[[k]]$links)
    expect_equal(cov_slices[[k]]$phi, cov_batches[[k]]$phi, tolerance = 1e-12)
  }

  rec_cov <- api$cp04_retained_if_release(store_cov)
  rm(prod_cov)
  expect_identical(store_cov$state, "released")
  expect_false(exists("phi", store_cov, inherits = FALSE))
  expect_identical(rec_cov$phase, "covariance")
  expect_identical(rec_cov$chunk_rows, B_cov)
  expect_identical(rec_cov$producer_passes, 1L)
  expect_identical(rec_cov$slice_passes, 1L)

  # Phase 2: GMM store fills at B_gmm on pass 1 and slices at B_gmm on pass 2.
  store_gmm <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B_gmm,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000,
    control = fx$ctrl, active_store = store_cov
  )
  expect_true(is.environment(store_gmm))
  expect_identical(store_gmm$phase, "gmm")
  expect_identical(store_gmm$max_slice_passes, 1L)
  expect_identical(store_gmm$B, B_gmm)

  gmm_batches <- list()
  prod_gmm <- api$cp04_retained_if_bind(b, store_gmm)
  api$cp04_replay_links(prod_gmm, links, function(z) {
    gmm_batches[[length(gmm_batches) + 1L]] <<- list(rows = z$rows, id = z$id, links = z$links, phi = z$phi)
  }, scale = "core", chunk_rows = B_gmm)
  expect_identical(store_gmm$state, "filled")
  expect_identical(store_gmm$producer_passes, 1L)
  expect_identical(store_gmm$slice_passes, 0L)
  expect_identical(length(gmm_batches), as.integer(ceiling(N / B_gmm)))
  expect_identical(store_gmm$id, sort(unique(fx$d$id)))
  expect_equal(store_gmm$phi, d92_phi_at_partition(b, links, B = B_gmm, scale = "core"), tolerance = 1e-12)
  gmm_original <- list()
  api$cp04_replay_links(b, links, function(z) {
    gmm_original[[length(gmm_original) + 1L]] <<- list(rows = z$rows, id = z$id, links = z$links, phi = z$phi)
  }, scale = "core", chunk_rows = B_gmm)
  expect_equal(gmm_batches, gmm_original, tolerance = 1e-12)
  expect_identical(unlist(lapply(gmm_batches, `[[`, "rows"), use.names = FALSE), seq_len(N))

  gmm_slices <- list()
  api$cp04_replay_links(prod_gmm, links, function(z) {
    gmm_slices[[length(gmm_slices) + 1L]] <<- list(rows = z$rows, id = z$id, links = z$links, phi = z$phi)
  }, scale = "core", chunk_rows = B_gmm)
  expect_identical(store_gmm$producer_passes, 1L)
  expect_identical(store_gmm$slice_passes, 1L)
  for (k in seq_along(gmm_batches)) {
    expect_identical(gmm_slices[[k]]$rows, gmm_batches[[k]]$rows)
    expect_identical(gmm_slices[[k]]$id, gmm_batches[[k]]$id)
    expect_identical(gmm_slices[[k]]$links, gmm_batches[[k]]$links)
    expect_equal(gmm_slices[[k]]$phi, gmm_batches[[k]]$phi, tolerance = 1e-12)
  }

  rec_gmm <- api$cp04_retained_if_release(store_gmm)
  rm(prod_gmm)
  expect_identical(store_gmm$state, "released")
  expect_false(exists("phi", store_gmm, inherits = FALSE))

  comp <- api$cp04_retained_if_phase_record(rec_cov, rec_gmm)
  expect_identical(comp$status, "retained")
  expect_identical(comp$mode, "phase_local")
  expect_identical(comp$producer_passes, 2L)
  expect_identical(comp$slice_passes, 2L)
  expect_identical(comp$streamed_passes, 0L)
  expect_identical(comp$chunk_rows, c(covariance = as.double(B_cov), gmm = as.double(B_gmm)))
  expect_identical(comp$simultaneous_stores, 1L)

  # Full pipeline parity with distinct covariance (B_cov = 9) and GMM (B_gmm = 7) partitions.
  ns_split <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L, base_ns = api)
  ns_stream <- d92_stream_ns(base_ns = ns_split)
  res_ret <- d92_run_dated_augmented(fx$d, chunk = NULL, ns = ns_split)
  res_str <- d92_run_dated_augmented(fx$d, chunk = NULL, ns = ns_stream)
  d92_equal(res_ret, res_str)

  rm_rec <- res_ret$resource_accounting$point_producer$retained_moment_if
  expect_identical(rm_rec$status, "retained")
  expect_identical(rm_rec$mode, "phase_local")
  expect_identical(rm_rec$producer_passes, 2L)
  expect_identical(rm_rec$slice_passes, 2L)
  expect_identical(rm_rec$streamed_passes, 0L)
  expect_identical(rm_rec$chunk_rows, c(covariance = 9, gmm = 7))
  expect_identical(rm_rec$simultaneous_stores, 1L)
  expect_identical(rm_rec$covariance$status, "retained")
  expect_equal(rm_rec$covariance$chunk_rows, 9, tolerance = 0)
  expect_identical(rm_rec$gmm$status, "retained")
  expect_equal(rm_rec$gmm$chunk_rows, 7, tolerance = 0)
  expect_identical(res_str$resource_accounting$point_producer$retained_moment_if$status, "streamed")
})

# ---------------------------------------------------------------------------
# 2. Independent covariance/GMM admission and exact one-entry, one-work-op and
#    planner-probe boundaries.
# ---------------------------------------------------------------------------

test_that("covariance and GMM phase-local stores enforce exact one-entry, one-work-op and planner boundaries independently", {
  old_cap <- d92_dated_fixture(chunk = 9L, ctrl = d92_ctrl(9L, work = 1e10))
  expect_error(api$cp04_stream_plan(old_cap$fit, seq_along(old_cap$fit$links), retained = 0),
    "CP04 covariance cumulative work", class = "cdid_ipw_owner_bound_error")
  phase_ns <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L)
  fx <- d92_dated_fixture(chunk = 9L, automatic_chunk_rows = TRUE)
  b <- fx$fit
  links <- seq_along(b$links)
  L <- length(links)
  st <- get("cp04_reader", b$reader_ref, inherits = FALSE)$state
  N <- st$N
  B_cov <- 9L
  B_gmm <- 7L
  cov_plan <- api$cp04_stream_plan(b, links, retained = 0)

  # Non-dated or non-transfer gates decline cleanly.
  expect_null(phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000,
    control = fx$ctrl, numerical_transfer = TRUE
  ))
  multi_cov <- cov_plan
  multi_cov$blocks <- list(1L, seq.int(2L, L))
  expect_null(phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = multi_cov, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  ))

  # Exact covariance phase store costs at B_cov.
  b_cov <- api$cp04_retained_if_boundary_cost(b, links, B_cov)
  entries_cov <- api$cp04_retained_if_entries(N, L, B_cov) + b_cov$entries
  checks_cov <- 32 + 8 * (ceiling(N / B_cov) + 1)
  work_cov <- api$cp04_retained_if_work(N, L) + checks_cov * b_cov$work
  cap <- min(fx$ctrl$max_retained_entries, fx$ctrl$max_matrix_entries)

  # Exact one-entry boundary for covariance store.
  s_cov_ok <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = 0, live = cap - entries_cov, work = 1000, control = fx$ctrl
  )
  expect_true(is.environment(s_cov_ok))
  expect_identical(s_cov_ok$entries, entries_cov)
  expect_identical(s_cov_ok$work, work_cov)
  api$cp04_retained_if_release(s_cov_ok)

  expect_null(phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = 0, live = cap - entries_cov + 1, work = 1000, control = fx$ctrl
  ))

  # Exact one-work-op boundary for covariance store.
  s_cov_wok <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = 0, live = 1000, work = fx$ctrl$max_work_ops - work_cov, control = fx$ctrl
  )
  expect_true(is.environment(s_cov_wok))
  api$cp04_retained_if_release(s_cov_wok)

  expect_null(phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = 0, live = 1000, work = fx$ctrl$max_work_ops - work_cov + 1, control = fx$ctrl
  ))

  # Exact planner-probe boundary for covariance store: when retained_probe
  # exceeds the 1-block covariance envelope by 1 entry, admission declines.
  own <- b$accounting$owned_entries + b$accounting$borrowed_entries
  env_cov <- api$augmented_if_envelope(b, B_cov, links, FALSE)
  consumer_cov <- if (L > fx$ctrl$max_chunk_rows) 6 * B_cov * L + 2 * L * L + 16 * L + 256 + b_cov$entries else 0
  room_cov <- fx$ctrl$max_matrix_entries - own - (env_cov$output + env_cov$scratch) - 4 * L * L - consumer_cov
  s_cov_plan_ok <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = room_cov - entries_cov, live = 1000, work = 1000, control = fx$ctrl
  )
  expect_true(is.environment(s_cov_plan_ok))
  api$cp04_retained_if_release(s_cov_plan_ok)

  expect_null(phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = room_cov - entries_cov + 1, live = 1000, work = 1000, control = fx$ctrl
  ))

  # Exact GMM phase store costs at B_gmm (and planner boundary at B_cov).
  b_gmm <- api$cp04_retained_if_boundary_cost(b, links, B_gmm)
  entries_gmm <- api$cp04_retained_if_entries(N, L, B_gmm) + b_gmm$entries
  checks_gmm <- 32 + 8 * (ceiling(N / B_gmm) + 1)
  work_gmm <- api$cp04_retained_if_work(N, L) + checks_gmm * b_gmm$work
  room_gmm <- fx$ctrl$max_matrix_entries - own - (env_cov$output + env_cov$scratch)

  same_partition_ns <- d92_distinct_partition_ns(B_cov = B_cov, B_gmm = B_cov)
  s_gmm_plan_ok <- same_partition_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = room_gmm - entries_cov, live = 1000, work = 1000, control = fx$ctrl
  )
  expect_true(is.environment(s_gmm_plan_ok))
  api$cp04_retained_if_release(s_gmm_plan_ok)

  expect_null(same_partition_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B_cov,
    cov_plan = cov_plan, retained = room_gmm - entries_cov + 1, live = 1000, work = 1000, control = fx$ctrl
  ))

  # Exact one-entry and one-work-op boundaries for GMM store at B_gmm.
  s_gmm_ok <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B_gmm,
    cov_plan = cov_plan, retained = 0, live = cap - entries_gmm,
    work = fx$ctrl$max_work_ops - work_gmm, control = fx$ctrl
  )
  expect_true(is.environment(s_gmm_ok))
  expect_identical(s_gmm_ok$entries, entries_gmm)
  expect_identical(s_gmm_ok$work, work_gmm)
  api$cp04_retained_if_release(s_gmm_ok)

  expect_null(phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B_gmm,
    cov_plan = cov_plan, retained = 0, live = cap - entries_gmm + 1,
    work = fx$ctrl$max_work_ops - work_gmm, control = fx$ctrl
  ))
  expect_null(phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B_gmm,
    cov_plan = cov_plan, retained = 0, live = cap - entries_gmm,
    work = fx$ctrl$max_work_ops - work_gmm + 1, control = fx$ctrl
  ))

  # Optional probes may swallow only finite requested > bound capacity errors.
  for (phase in c("covariance", "gmm")) {
    probe_ns <- new.env(parent = phase_ns)
    reserve <- api$cp04_retained_if_phase_reserve
    environment(reserve) <- probe_ns
    reserve_probe <- function() reserve(b, links, phase = phase,
      chunk_rows = if (phase == "covariance") B_cov else B_gmm,
      cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl)
    probe_name <- if (phase == "covariance") "cp04_stream_plan" else "cp04_preflight_links"
    assign(probe_name, function(...) api$ipw_owner_safe_sum(2, 1, "finite probe capacity"), probe_ns)
    expect_null(reserve_probe())
    assign(probe_name, function(...) api$ipw_owner_int_bound(NA_real_, "invalid probe count"), probe_ns)
    expect_error(reserve_probe(), class = "cdid_ipw_owner_control_error")
    assign(probe_name, function(...) api$ipw_owner_safe_sum(c(2^53 - 1, 1)), probe_ns)
    expect_error(reserve_probe(), class = "cdid_ipw_owner_bound_error")
    assign(probe_name, function(...) api$ipw_owner_check_matrix_dims(2^31, 1), probe_ns)
    expect_error(reserve_probe(), class = "cdid_ipw_owner_matrix_dim_error")
    assign(probe_name, function(...) api$cp04_abort("probe integrity", "integrity"), probe_ns)
    expect_error(reserve_probe(), class = "cdid_cp04_integrity_error")
    assign(probe_name, function(...) stop("probe numerical failure"), probe_ns)
    expect_error(reserve_probe(), "probe numerical failure")
  }

  # Independent pipeline admission: covariance-only admitted vs GMM-only admitted vs forced stream.
  ns_split <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L, base_ns = api)
  ns_cov_only <- d92_phase_filter_ns("covariance", base_ns = ns_split)
  ns_gmm_only <- d92_phase_filter_ns("gmm", base_ns = ns_split)
  ns_stream <- d92_stream_ns(base_ns = ns_split)

  res_cov_only <- d92_run_dated_augmented(fx$d, chunk = NULL, ns = ns_cov_only)
  res_gmm_only <- d92_run_dated_augmented(fx$d, chunk = NULL, ns = ns_gmm_only)
  res_stream <- d92_run_dated_augmented(fx$d, chunk = NULL, ns = ns_stream)

  d92_equal(res_cov_only, res_stream)
  d92_equal(res_gmm_only, res_stream)

  rec_co <- res_cov_only$resource_accounting$point_producer$retained_moment_if
  expect_identical(rec_co$mode, "phase_local")
  expect_identical(rec_co$producer_passes, 1L)
  expect_identical(rec_co$slice_passes, 1L)
  expect_identical(rec_co$streamed_passes, 2L)
  expect_identical(rec_co$covariance$status, "retained")
  expect_identical(rec_co$gmm$status, "streamed")

  rec_go <- res_gmm_only$resource_accounting$point_producer$retained_moment_if
  expect_identical(rec_go$mode, "phase_local")
  expect_identical(rec_go$producer_passes, 1L)
  expect_identical(rec_go$slice_passes, 1L)
  expect_identical(rec_go$streamed_passes, 2L)
  expect_identical(rec_go$covariance$status, "streamed")
  expect_identical(rec_go$gmm$status, "retained")
})

# ---------------------------------------------------------------------------
# 3. No simultaneous stores and strict phase lifetime / release discipline.
# ---------------------------------------------------------------------------

test_that("covariance store is released before GMM consumers and GMM store is released before transport, auxiliary and diagnostics", {
  phase_ns <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L)
  fx <- d92_dated_fixture(chunk = 9L)
  b <- fx$fit
  links <- seq_along(b$links)
  cov_plan <- api$cp04_stream_plan(b, links, retained = 0)

  store_cov <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "covariance", chunk_rows = 9L,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  # Refuses reserving a second store while active_store is not released.
  expect_error(
    phase_ns$cp04_retained_if_phase_reserve(
      b, links, phase = "gmm", chunk_rows = 7L,
      cov_plan = cov_plan, retained = 0, live = 1000, work = 1000,
      control = fx$ctrl, active_store = store_cov
    ),
    class = "cdid_cp04_integrity_error"
  )
  # Refuses reserving on an already-bound carrier.
  bound_cov <- api$cp04_retained_if_bind(b, store_cov)
  expect_error(
    phase_ns$cp04_retained_if_phase_reserve(
      bound_cov, links, phase = "gmm", chunk_rows = 7L,
      cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
    ),
    class = "cdid_cp04_integrity_error"
  )
  rm(bound_cov)
  api$cp04_retained_if_release(store_cov)

  # Wrong-phase use in cp04_stream_covariance is refused and releases the store.
  same_partition_ns <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 9L)
  store_gmm <- same_partition_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = 9L,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  expect_error(
    api$cp04_stream_covariance(api$cp04_retained_if_bind(b, store_gmm), links, retained = 0),
    class = "cdid_cp04_integrity_error"
  )
  expect_identical(store_gmm$state, "released")

  # Phase-local store refuses a second slice pass and refuses auxiliary extension.
  store_gmm2 <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = 7L,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  prod2 <- api$cp04_retained_if_bind(b, store_gmm2)
  api$cp04_replay_links(prod2, links, function(z) NULL, scale = "core", chunk_rows = 7L)
  api$cp04_replay_links(prod2, links, function(z) NULL, scale = "core", chunk_rows = 7L)
  expect_null(api$cp04_retained_if_auxiliary_admit(
    store_gmm2, b, links, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  ))
  expect_error(
    api$cp04_replay_links(prod2, links, function(z) NULL, scale = "core", chunk_rows = 7L),
    class = "cdid_cp04_integrity_error"
  )
  expect_identical(store_gmm2$state, "released")

  # Pipeline lifetime trace: verify covariance store is released before
  # solve_chained_gmm_core and GMM store is released before transport,
  # build_ipw_auxiliary_reporting and cp04_public_diagnostics.
  ns_trace <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L, base_ns = api)
  created_stores <- list()
  trace_reserve <- ns_trace$cp04_retained_if_phase_reserve
  ns_trace$cp04_retained_if_phase_reserve <- function(backend, links, phase, ...) {
    if (identical(phase, "gmm") && !is.null(created_stores$covariance)) {
      expect_identical(created_stores$covariance$state, "released")
      expect_false(exists("phi", created_stores$covariance, inherits = FALSE))
    }
    st_out <- trace_reserve(backend, links, phase = phase, ...)
    if (!is.null(st_out)) created_stores[[phase]] <<- st_out
    st_out
  }
  gmm_entry_states <- NULL
  ns_trace$solve_chained_gmm_core <- function(...) {
    gmm_entry_states <<- list(
      cov_state = created_stores$covariance$state,
      cov_has_phi = exists("phi", created_stores$covariance, inherits = FALSE),
      gmm_state = created_stores$gmm$state
    )
    api$solve_chained_gmm_core(...)
  }
  aux_entry_states <- NULL
  ns_trace$build_ipw_auxiliary_reporting <- function(...) {
    if (!is.null(created_stores$gmm)) {
      aux_entry_states <<- list(
        cov_state = created_stores$covariance$state,
        gmm_state = created_stores$gmm$state,
        gmm_has_phi = exists("phi", created_stores$gmm, inherits = FALSE)
      )
    }
    api$build_ipw_auxiliary_reporting(...)
  }
  diag_entry_states <- NULL
  ns_trace$cp04_public_diagnostics <- function(...) {
    diag_entry_states <<- list(
      cov_state = created_stores$covariance$state,
      gmm_state = created_stores$gmm$state,
      gmm_has_phi = exists("phi", created_stores$gmm, inherits = FALSE)
    )
    api$cp04_public_diagnostics(...)
  }

  res <- d92_run_dated_augmented(fx$d, chunk = NULL, ns = ns_trace)
  expect_identical(gmm_entry_states$cov_state, "released")
  expect_false(gmm_entry_states$cov_has_phi)
  expect_identical(gmm_entry_states$gmm_state, "reserved")
  expect_identical(aux_entry_states$cov_state, "released")
  expect_identical(aux_entry_states$gmm_state, "released")
  expect_false(aux_entry_states$gmm_has_phi)
  expect_identical(diag_entry_states$cov_state, "released")
  expect_identical(diag_entry_states$gmm_state, "released")
  expect_false(diag_entry_states$gmm_has_phi)
  for (out in created_stores) for (nm in c("phi", "raw_ref", "backend_snapshot", "reference_snapshot"))
    expect_false(exists(nm, out, inherits = FALSE))
  expect_false(d92_find_store(res))

  # Pipeline on.exit cleanup also covers failures before a phase's first replay.
  for (failed_phase in c("covariance", "gmm")) {
    ns_error <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L)
    error_stores <- list()
    error_reserve <- ns_error$cp04_retained_if_phase_reserve
    ns_error$cp04_retained_if_phase_reserve <- function(backend, links, phase, ...) {
      out <- error_reserve(backend, links, phase = phase, ...)
      if (!is.null(out)) error_stores[[phase]] <<- out
      out
    }
    fail_at <- if (failed_phase == "covariance") "cp04_stream_covariance" else "solve_chained_gmm_core"
    assign(fail_at, function(...) stop("phase consumer failure"), ns_error)
    expect_error(d92_run_dated_augmented(fx$d, ns = ns_error), "phase consumer failure")
    expect_true(is.environment(error_stores[[failed_phase]]))
    for (out in error_stores) {
      expect_identical(out$state, "released")
      for (nm in c("phi", "raw_ref", "backend_snapshot", "reference_snapshot"))
        expect_false(exists(nm, out, inherits = FALSE))
    }
  }
})

# ---------------------------------------------------------------------------
# 4. Owner/raw mutation checks before and after every callback including the
#    last callback, and cleanup on errors.
# ---------------------------------------------------------------------------

test_that("phase-local stores verify raw/owner identity before and after every callback including the last and release on errors", {
  phase_ns <- d92_distinct_partition_ns(B_cov = 9L, B_gmm = 7L)
  fx <- d92_dated_fixture(chunk = 9L)
  b <- fx$fit
  links <- seq_along(b$links)
  N <- get("cp04_reader", b$reader_ref, inherits = FALSE)$state$N
  B <- 7L
  n_chunks <- as.integer(ceiling(N / B))
  cov_plan <- api$cp04_stream_plan(b, links, retained = 0)

  # 1. Mutation during the LAST callback of the fill pass is caught after that callback.
  store_last_fill <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  prod_last_fill <- api$cp04_retained_if_bind(b, store_last_fill)
  seen_fill <- 0L
  expect_error(
    api$cp04_replay_links(prod_last_fill, links, function(z) {
      seen_fill <<- seen_fill + 1L
      if (seen_fill == n_chunks) store_last_fill$phi <- matrix(0, 1L, 1L)
    }, scale = "core", chunk_rows = B),
    class = "cdid_cp04_integrity_error"
  )
  expect_identical(seen_fill, n_chunks)
  expect_identical(store_last_fill$state, "released")
  expect_false(exists("phi", store_last_fill, inherits = FALSE))

  # 2. Copy-on-write raw replacement during the LAST callback of the slice pass is caught after that callback.
  store_last_slice <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  prod_last_slice <- api$cp04_retained_if_bind(b, store_last_slice)
  api$cp04_replay_links(prod_last_slice, links, function(z) NULL, scale = "core", chunk_rows = B)
  expect_identical(store_last_slice$state, "filled")
  seen_slice <- 0L
  expect_error(
    api$cp04_replay_links(prod_last_slice, links, function(z) {
      seen_slice <<- seen_slice + 1L
      if (seen_slice == n_chunks) store_last_slice$phi[1L] <- store_last_slice$phi[1L] + 0
    }, scale = "core", chunk_rows = B),
    class = "cdid_cp04_integrity_error"
  )
  expect_identical(seen_slice, n_chunks)
  expect_identical(store_last_slice$state, "released")
  expect_false(exists("phi", store_last_slice, inherits = FALSE))

  # Replacement between passes is rejected before any slice callback runs.
  store_before_slice <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  prod_before_slice <- api$cp04_retained_if_bind(b, store_before_slice)
  api$cp04_replay_links(prod_before_slice, links, function(z) NULL, scale = "core", chunk_rows = B)
  store_before_slice$phi <- store_before_slice$phi + 0
  seen_before_slice <- 0L
  expect_error(
    api$cp04_replay_links(prod_before_slice, links, function(z) {
      seen_before_slice <<- seen_before_slice + 1L
    }, scale = "core", chunk_rows = B),
    class = "cdid_cp04_integrity_error"
  )
  expect_identical(seen_before_slice, 0L)
  expect_identical(store_before_slice$state, "released")

  # 3. Owner / reference_reader mutation on the last callback is caught and releases the store.
  fx_mut <- d92_dated_fixture(chunk = 9L)
  b_mut <- fx_mut$fit
  cov_mut <- api$cp04_stream_plan(b_mut, links, retained = 0)
  store_owner <- phase_ns$cp04_retained_if_phase_reserve(
    b_mut, links, phase = "covariance", chunk_rows = 9L,
    cov_plan = cov_mut, retained = 0, live = 1000, work = 1000, control = fx_mut$ctrl
  )
  prod_owner <- api$cp04_retained_if_bind(b_mut, store_owner)
  n_cov_chunks <- as.integer(ceiling(N / 9L))
  seen_cov <- 0L
  expect_error(
    api$cp04_replay_links(prod_owner, links, function(z) {
      seen_cov <<- seen_cov + 1L
      if (seen_cov == n_cov_chunks) {
        unlockBinding("reference_reader", b_mut$reference_ref)
        ref_copy <- get("reference_reader", b_mut$reference_ref, inherits = FALSE)
        ref_copy$w <- ref_copy$w + 0
        assign("reference_reader", ref_copy, b_mut$reference_ref)
      }
    }, scale = "core", chunk_rows = 9L),
    class = "cdid_reference_integrity_error"
  )
  expect_identical(seen_cov, n_cov_chunks)
  expect_identical(store_owner$state, "released")
  expect_false(exists("phi", store_owner, inherits = FALSE))

  # 4. Callback error during fill or slice releases the phase-local store via on.exit.
  store_err <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  prod_err <- api$cp04_retained_if_bind(b, store_err)
  expect_error(
    api$cp04_replay_links(prod_err, links, function(z) stop("consumer failure"), scale = "core", chunk_rows = B),
    "consumer failure"
  )
  expect_identical(store_err$state, "released")
  expect_false(exists("phi", store_err, inherits = FALSE))
  expect_false(exists("backend_snapshot", store_err, inherits = FALSE))

  store_slice_err <- phase_ns$cp04_retained_if_phase_reserve(
    b, links, phase = "gmm", chunk_rows = B,
    cov_plan = cov_plan, retained = 0, live = 1000, work = 1000, control = fx$ctrl
  )
  prod_slice_err <- api$cp04_retained_if_bind(b, store_slice_err)
  api$cp04_replay_links(prod_slice_err, links, function(z) NULL, scale = "core", chunk_rows = B)
  expect_error(
    api$cp04_replay_links(prod_slice_err, links, function(z) stop("slice consumer failure"), scale = "core", chunk_rows = B),
    "slice consumer failure"
  )
  expect_identical(store_slice_err$state, "released")
  for (nm in c("phi", "raw_ref", "backend_snapshot", "reference_snapshot"))
    expect_false(exists(nm, store_slice_err, inherits = FALSE))
})
