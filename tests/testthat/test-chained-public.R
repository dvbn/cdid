# Public CP02 candidate surface: att_gt()/aggte() forwarding and dispatch, the
# scope gate, the D20 scale-aware numerical acceptance and the public control
# argument. Deterministic helper tests on fixed matrices, the frozen N200
# reporting fields and typed refusals. Ordinary delegation uses the accepted
# parity fixture draws/fits only; the other blocks have no DGP draw or fit. Parts that need
# did 2.5.1 (ordinary delegation, parity) skip honestly when did is absent and
# run on the installed host.

pub_locate_fixture <- function(filename) {
  candidates <- c(tryCatch(testthat::test_path("fixtures", filename), error = function(e) character(0L)),
    file.path("fixtures", filename), file.path("tests", "testthat", "fixtures", filename), file.path("..", "..", "tests", "testthat", "fixtures", filename))
  for (p in candidates) if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  stop("Fixture file not found: ", filename)
}
pub_budget <- function() list(max_units = 2000, max_clusters = 2000, max_targets = 100, max_draws = 2000, max_if_batch_rows = 64,
  max_if_batch_cells = 6400, max_sign_batch_rows = 32, max_sign_batch_cells = 64000, max_workspace_cells = 1000000, max_workspace_bytes = 8000000)
pub_abs <- function() list(max_estimate_error_abs = 1e-3, max_se_error_abs = 1e-2, max_vcov_error_abs = 1e-2, max_multiplier_error_abs = 1e-1, max_band_error_abs = 1e-1)
pub_rel <- function(level = 1e-6) list(max_relative_error_lower_se = level)
# Ten units in five clusters (the request-path fixture): distinct nonzero cluster sums.
pub_fixed_if <- function() { X <- rbind(c(-2, 4), c(1, 1), c(-1, 2), c(1, -3), c(1, -2), c(2, 1), c(2, -4), c(-1, 2), c(0.5, 1), c(-3.5, -3)); colnames(X) <- c("a", "b"); X }
pub_exact_enclosure <- function(th, X) list(theta = list(lo = matrix(th, 1L), hi = matrix(th, 1L)), rows = list(lo = X, hi = X))
# Explicit fixture metadata for the universal-reference operator: roles are
# declared independently of numerical zeros and target names.
pub_reference_meta <- function(ids, reference) data.frame(target_id = ids,
  role = ifelse(seq_along(ids) %in% reference, "reference", "post"),
  support_status = ifelse(seq_along(ids) %in% reference, "deterministic_reference_zero", "estimable"),
  time = ifelse(seq_along(ids) %in% reference, 1, 2), ref = 1, base_t = 1,
  stringsAsFactors = FALSE)
pub_capture_warnings <- function(expr) { w <- character(); val <- withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }); list(value = val, warnings = w) }

# ------------------------------------------------------------------------------
# A. Public surface and scope gate (no did needed)
# ------------------------------------------------------------------------------

test_that("att_gt appends control = NULL after the chained-only extensions; did's formals and order are intact", {
  f <- formals(att_gt)
  n_up <- length(upstream_att_gt_formals)
  expect_identical(names(f)[seq_len(n_up)], upstream_att_gt_formals)
  expect_identical(names(f)[n_up + 1L], "...")
  expect_identical(names(f)[-seq_len(n_up + 1L)], cdid_chained_only_formals)
  expect_identical(names(f)[length(f)], "control"); expect_null(f$control)
  expect_identical(f$est_method, "dr"); expect_identical(f$bstrap, TRUE); expect_identical(f$cband, TRUE); expect_identical(f$biters, 1000)
  expect_identical(names(formals(aggte)), upstream_aggte_formals)          # aggte signature unchanged
})

test_that("the chained scope gate refuses before the data promise is forced, with the familiar classes", {
  base <- list(yname = "y", tname = "period", idname = "id", gname = "G", data = quote(stop("data must not be evaluated")))
  refuse <- function(cls, ...) {
    e <- expect_error(do.call(att_gt, c(base, list(...))), class = cls)
    expect_no_match(conditionMessage(e), "data must not be evaluated")
    invisible(e)
  }
  # D48 CP07 activated the "dr" default; unsupported methods are never substituted and refuse with the CP01 class retained
  e <- refuse("cdid_chained_request_est_method", est_method = "stdipw")
  expect_s3_class(e, "cdid_chained_unavailable")
  expect_match(conditionMessage(e), "est_method = \"dr\"", fixed = TRUE)
  expect_match(conditionMessage(e), "est_method = \"ipw\"", fixed = TRUE); expect_match(conditionMessage(e), "mode = \"ordinary\"", fixed = TRUE)
  expect_match(conditionMessage(e), "without silent replacement", fixed = TRUE)
  for (m in list("stdipw", NA_character_, c("ipw", "dr"))) refuse("cdid_chained_unavailable", est_method = m)
  refuse("cdid_chained_request_unsupported_selection", est_method = "ipw", selection = "ipw")
  refuse("cdid_chained_request_unsupported_selection", est_method = function(...) NULL, selection = "ipw",
         unknown = quote(stop("extra must not be evaluated")))
  refuse("cdid_chained_request_unsupported_selection", est_method = "ipw", selection_xformla = ~x)
  refuse("cdid_chained_request_pt_start", est_method = "ipw", pt_assumption = "all")
  # D60: "varying" stays refused; "first_period"/"base_period" pass the gate and reach the data.
  refuse("cdid_chained_request_unsupported_weights", est_method = "ipw", fix_weights = "varying")
  for (fw in c("first_period", "base_period"))
    expect_error(do.call(att_gt, c(base, list(est_method = "ipw", fix_weights = fw))), "data must not be evaluated")
  # CP08 inherited options: point estimates only, declared bootstrap workers and
  # progress pass the gate; undeclared workers and invalid flags still refuse.
  for (ok in list(list(compute_inffunc = FALSE), list(pl = TRUE), list(cores = 2), list(print_details = TRUE),
                  list(pl = TRUE, cores = 2, control = list(workers = 2))))
    expect_error(do.call(att_gt, c(base, list(est_method = "ipw"), ok)), "data must not be evaluated")
  refuse("cdid_chained_request_parallel_error", est_method = "ipw", pl = TRUE, cores = 2)
  refuse("cdid_chained_request_input_error", est_method = "ipw", compute_inffunc = NA)
  refuse("cdid_chained_request_input_error", est_method = "ipw", cores = 0)
  refuse("cdid_chained_request_input_error", est_method = "ipw", unknown = quote(stop("extra must not be evaluated")))
  refuse("cdid_chained_request_input_error", est_method = "ipw", shift = 1)                      # unknown dots are never ignored for est_method = "ipw"
  refuse("cdid_chained_request_input_error", est_method = "ipw", control_group = "everyone")
  refuse("cdid_chained_panel", est_method = "ipw", panel = FALSE, mode = "chained")
  refuse("cdid_chained_balanced", est_method = "ipw", allow_unbalanced_panel = FALSE)
  refuse("cdid_chained_balanced", allow_unbalanced_panel = FALSE)                                # before the est_method gate, as in CP01
  refuse("cdid_mode_error", est_method = "ipw", mode = "Chained")
  # the public control is validated before the data is touched
  refuse("cdid_chained_public_control_error", est_method = "ipw", control = list(bogus = 1))
  refuse("cdid_chained_public_control_error", est_method = "ipw", control = list(workers = 0))
  refuse("cdid_chained_public_control_error", est_method = "ipw", control = list(workers = 1.5))
  refuse("cdid_chained_public_control_error", est_method = "ipw", control = structure(list(), class = "cdid_ipw_pipeline_control"))
  # a passing gate reaches the data (which this stub refuses): nothing else stops the request
  expect_error(do.call(att_gt, base), "data must not be evaluated")
  expect_error(do.call(att_gt, c(base, list(est_method = "dr"))), "data must not be evaluated")
  expect_error(do.call(att_gt, c(base, list(est_method = "ipw"))), "data must not be evaluated")
  expect_error(do.call(att_gt, c(base, list(est_method = "reg"))), "data must not be evaluated")
  expect_error(do.call(att_gt, c(base, list(est_method = function(...) NULL, shift = 1))), "data must not be evaluated")
  # D23: reject the deferred override explicitly before forcing data.
  refuse("cdid_chained_public_accuracy_deferred", est_method = "ipw", control = list(relative_accuracy = 1e-5))
  expect_error(do.call(att_gt, c(base, list(est_method = "ipw", bstrap = FALSE, cband = TRUE))), "data must not be evaluated")
})

test_that("aggte dispatches cdid_MP objects to the chained no-refit path and leaves others to did", {
  fake <- structure(list(), class = c("cdid_MP", "MP"))
  e <- expect_error(aggte(fake, type = "group"), class = "cdid_chained_request_input_error")
  expect_match(conditionMessage(e), "aggregation state")
  expect_error(aggte(fake, type = "nonsense"), class = "cdid_chained_request_input_error")
  if (!requireNamespace("did", quietly = TRUE)) {
    # without did the ordinary branch fails only at did's own namespace, proving no chained validation ran
    e2 <- tryCatch(aggte(structure(list(), class = "MP")), error = function(e) e)
    expect_false(inherits(e2, "cdid_error"))
  }
})

test_that("ordinary mode neither evaluates nor forwards control or any chained-only argument", {
  skip_if_not_installed("did")
  d <- parity_panel()
  ours <- with_conditions(cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d, bstrap = FALSE,
    mode = "ordinary", control = stop("forced control"), frame = stop("forced frame"),
    pt_start = stop("forced pt_start")))
  ref <- with_conditions(did::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d, bstrap = FALSE))
  expect_same_result(ours, ref, mp_fields)
  ours2 <- with_conditions(cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d, bstrap = FALSE,
    mode = "ordinary", control = list(bogus = 1)))
  expect_same_result(ours2, ref, mp_fields)
  # ordinary aggregation stays did's own, lazily forwarded
  set.seed(4); a <- cdid::aggte(ours$value, type = "simple"); set.seed(4); b <- did::aggte(ref$value, type = "simple")
  expect_identical(a$overall.att, b$overall.att); expect_identical(a$overall.se, b$overall.se)
})

# ------------------------------------------------------------------------------
# B. D20 scale-aware numerical acceptance (copied native DLL; no did)
# ------------------------------------------------------------------------------

test_that("inference-accuracy budgets accept the legacy absolute form and the D20 level; everything else refuses", {
  expect_identical(validate_ipw_pipeline_inference_accuracy(pub_abs()), pub_abs())
  expect_identical(validate_ipw_pipeline_inference_accuracy(pub_rel()), pub_rel())
  both <- c(pub_abs(), pub_rel()); expect_identical(validate_ipw_pipeline_inference_accuracy(both), both)
  expect_true(ipw_pipeline_accuracy_is_relative(pub_rel())); expect_false(ipw_pipeline_accuracy_is_relative(pub_abs()))
  refuse <- function(x) expect_error(validate_ipw_pipeline_inference_accuracy(x), class = "cdid_ipw_pipeline_control_error")
  refuse(pub_abs()[-1L]); refuse(list(max_relative_error_lower_se = 0)); refuse(list(max_relative_error_lower_se = 1))
  refuse(list(max_relative_error_lower_se = -1e-6)); refuse(list(max_relative_error_lower_se = Inf)); refuse(c(pub_rel(), list(extra = 1)))
  refuse(list(max_relative_error_lower_se = 1e-6, max_se_error_abs = -1)); refuse(structure(pub_rel(), class = "x"))
  # D23: the public control requests no certification; explicit private D20
  # budgets and their mathematical tests above remain research functionality.
  d <- expand.grid(id = 1:6, period = 1:4); d$G <- c(0, 0, 3, 3, 4, 4)[d$id]; d$x <- c(0.5, -1, 1, 0, 2, -2)[d$id]; d$y <- d$id + d$period
  prep <- chained_request_prepare_inputs(d, NULL, "y", "period", "id", "G", ~x)
  pc <- chained_public_pipeline_control(chained_public_dimensions(prep, 50), chained_public_resolve_control(NULL))
  expect_null(pc$control$inference_accuracy)
})

test_that("under the D20 level the reported numbers are unchanged, the certificate gains normalized ratios, and refusals are typed", {
  X <- pub_fixed_if(); th <- c(a = 0.5, b = -0.2); cl <- rep(1:5, each = 2); enc <- pub_exact_enclosure(th, X)
  run <- function(acc, ...) compute_chained_ipw_inference(th, X, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = acc, enclosure = enc, ...)
  leg <- run(pub_abs(), cband = FALSE, bstrap = FALSE); rel <- run(pub_rel(), cband = FALSE, bstrap = FALSE)
  # estimates, SEs, covariance, intervals and the absolute error record are identical: D20 only accepts or refuses
  for (f in c("estimates", "se", "analytical_se", "vcov", "ci_lower", "ci_upper", "cband_crit_val", "cband_lower", "cband_upper", "band_family"))
    expect_identical(rel[[f]], leg[[f]], info = f)
  expect_identical(rel$numerical_certificate$errors, leg$numerical_certificate$errors)
  expect_null(leg$numerical_certificate$d20); d20 <- rel$numerical_certificate$d20
  expect_identical(d20$decision, "D20"); expect_identical(d20$level, 1e-6); expect_true(d20$certified)
  expect_identical(names(d20$worst_ratios), c("estimate", "se", "ci_lower", "ci_upper", "band_lower", "band_upper", "covariance"))
  expect_true(all(d20$worst_ratios <= 1e-6)); expect_identical(names(d20$lower_se_scale), c("a", "b"))
  expect_identical(unname(d20$lower_se_scale), unname(rel$numerical_certificate$se[1L, ]))
  expect_true(all(d20$lower_se_scale > 0)); expect_identical(d20$structural_zero_targets, character(0L))
  expect_identical(dim(d20$ratios_by_target), c(6L, 2L)); expect_identical(dim(d20$covariance_ratios), c(2L, 2L))
  expect_match(rel$numerical_certificate$scope, "D20")
  # the ratios are the recorded errors over the lower scales (outward division), hand-checked
  expect_true(all(abs(d20$ratios_by_target["estimate", ] - d20$errors_by_target["estimate", ] / d20$lower_se_scale) <= 1e-15 * (1 + d20$ratios_by_target["estimate", ])))
  # bootstrap law: scale is the IQR-SE interval's lower endpoint; raw draws use n/sqrt(G) x analytic lower scale
  set.seed(5); tt <- run(pub_rel(), cband = TRUE, bstrap = TRUE, caller_rng = TRUE)
  set.seed(5); tt_leg <- run(pub_abs(), cband = TRUE, bstrap = TRUE, caller_rng = TRUE)
  expect_identical(tt$se, tt_leg$se); expect_identical(tt$cband_crit_val, tt_leg$cband_crit_val)
  cb <- tt$numerical_certificate
  expect_identical(unname(cb$d20$lower_se_scale), unname(cb$se[1L, ])); expect_false(isTRUE(all.equal(cb$se, cb$analytic_se)))
  expect_equal(unname(cb$d20$raw_draw_scale_lower), unname(10 / sqrt(5) * cb$analytic_se[1L, ]), tolerance = 1e-15)
  expect_identical(length(cb$multiplier_error_by_target), 2L); expect_true(max(cb$d20$multiplier_ratios) <= 1e-6)
  # widened estimate enclosure: refusal with the D20 class, diagnostics name the worst family, nothing relaxed
  wide <- enc; wide$theta$lo <- wide$theta$lo - 5e-4; wide$theta$hi <- wide$theta$hi + 5e-4
  e <- expect_error(compute_chained_ipw_inference(th, X, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = pub_rel(),
    enclosure = wide, cband = FALSE, bstrap = FALSE), class = "cdid_ipw_pipeline_d20_error")
  expect_s3_class(e, "cdid_ipw_pipeline_accuracy_error"); expect_match(conditionMessage(e), "exceeds the level 1e-06")
  expect_true(e$diagnostics$worst_ratios[["estimate"]] > 1e-6); expect_identical(e$diagnostics$level, 1e-6)
  # the legacy absolute budgets accept that same enclosure (5e-4 < 1e-3 abs): D20 is the stricter public rule, not a relabelled toy budget
  expect_true(compute_chained_ipw_inference(th, X, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = pub_abs(),
    enclosure = wide, cband = FALSE, bstrap = FALSE)$numerical_certificate$certified)
  # an explicit larger level (recorded by the caller) accepts; the default never moves
  expect_true(compute_chained_ipw_inference(th, X, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = pub_rel(1e-2),
    enclosure = wide, cband = FALSE, bstrap = FALSE)$numerical_certificate$d20$certified)
  # a nonzero constant IF column has no positive scale: the existing scale refusal stays (rank alone never refuses)
  Xc <- cbind(X, c = 1); thc <- c(th, c = 0.7)
  expect_error(compute_chained_ipw_inference(thc, Xc, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = pub_rel(),
    enclosure = pub_exact_enclosure(thc, Xc), cband = FALSE, bstrap = FALSE), class = "cdid_ipw_pipeline_scale_error")
})

test_that("exact structural reference zeros are identities under D20: no scale is inverted, their boxes must contain zero", {
  X <- pub_fixed_if(); cl <- rep(1:5, each = 2)
  Z <- cbind(X, ref = 0); thz <- c(a = 0.5, b = -0.2, ref = 0)
  dz <- compute_chained_ipw_inference(thz, Z, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = pub_rel(),
    enclosure = pub_exact_enclosure(thz, Z), cband = FALSE, bstrap = FALSE, target_meta = pub_reference_meta(names(thz), 3L))
  d20 <- dz$numerical_certificate$d20
  expect_identical(d20$structural_zero_targets, "ref"); expect_identical(unname(d20$lower_se_scale["ref"]), 0)
  expect_true(all(d20$ratios_by_target[, "ref"] == 0)); expect_true(all(d20$covariance_ratios["ref", ] == 0))
  expect_true(all(d20$worst_ratios <= 1e-6)); expect_identical(dz$numerical_certificate$deterministic, 3L)
  # A zero computed column whose exact enclosure is not a zero identity must refuse.
  uncertain <- pub_exact_enclosure(thz, Z)
  uncertain$theta$lo[3L] <- -1e-12; uncertain$theta$hi[3L] <- 1e-12
  expect_error(compute_chained_ipw_inference(thz, Z, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = pub_rel(),
    enclosure = uncertain, cband = FALSE, bstrap = FALSE, target_meta = pub_reference_meta(names(thz), 3L)), class = "cdid_ipw_pipeline_d20_error")
  uncertain_if <- pub_exact_enclosure(thz, Z)
  uncertain_if$rows$lo[1L, 3L] <- -1e-12; uncertain_if$rows$hi[1L, 3L] <- 1e-12
  expect_error(compute_chained_ipw_inference(thz, Z, cl, 0.05, biters = 32L, multiplier_budget = pub_budget(), inference_accuracy = pub_rel(),
    enclosure = uncertain_if, cband = FALSE, bstrap = FALSE, target_meta = pub_reference_meta(names(thz), 3L)), class = "cdid_ipw_pipeline_d20_error")
  # a reference-only selection is deterministic end to end (no quantile, no draw)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  det <- compute_chained_ipw_inference(c(ref = 0), matrix(0, 10L, 1L, dimnames = list(NULL, "ref")), cl, 0.05, biters = 32L, cband = TRUE, bstrap = TRUE,
    caller_rng = TRUE, multiplier_budget = pub_budget(), inference_accuracy = pub_rel(), enclosure = pub_exact_enclosure(0, matrix(0, 10L, 1L)), target_meta = pub_reference_meta("ref", 1L))
  expect_false(exists(".Random.seed", envir = globalenv())); expect_true(det$numerical_certificate$d20$certified)
  expect_identical(det$numerical_certificate$d20$structural_zero_targets, "ref"); expect_identical(det$band_family, "deterministic_exact_zero_over_certified_ids")
})

test_that("frozen N200 reporting fields: analytic families meet 1e-6 by two orders of magnitude; the multiplier band on this fixture does not", {
  z <- readRDS(pub_locate_fixture("cp02_result_adapter_N200.rds")); rep_ <- z$reporting; fr <- z$frame; G <- fr$G; G[G == 0] <- Inf
  enc <- list(theta = rep_$numerical_certificate$theta, rows = rep_$numerical_certificate$rows)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  # reported cells (analytic, no draw): same SEs as the frozen record, D20 ratios ~1e-8 (the decision document's N200 magnitude)
  ci <- ipw_pipeline_supported_inference(rep_$estimates, rep_$if_full, rep_$available, NULL, 0.05, FALSE, 1000L, "rademacher", NULL, pub_budget(), pub_rel(), enc,
    NULL, fr$eligible_ge2, bstrap = FALSE, caller_rng = TRUE)
  expect_equal(unname(ci$se), unname(rep_$se), tolerance = 1e-14); expect_true(ci$numerical_certificate$certified)
  expect_true(max(ci$numerical_certificate$d20$worst_ratios) < 1e-7); expect_true(all(ci$numerical_certificate$d20$lower_se_scale > 0))
  # the four analytic aggregation families and their overalls, without refit, no RNG use
  for (typ in c("simple", "group", "dynamic", "calendar")) {
    a <- aggregate_chained_ipw_targets(rep_$estimates, rep_$if_full, rep_$target_meta, G, fr$sample_weight, type = typ, multiplier_budget = pub_budget(),
      inference_accuracy = pub_rel(), available = rep_$available, enclosure = enc, bstrap = FALSE, caller_rng = TRUE, calendar = z$calendar, cohort_grid = sort(unique(G)))
    o <- z$aggregations[[typ]]
    expect_identical(a$status, o$status, info = typ); expect_identical(unname(a$estimates), unname(o$estimates), info = typ)
    expect_equal(unname(a$se), unname(o$se), tolerance = 1e-14, info = typ); expect_equal(a$overall$se, o$overall$se, tolerance = 1e-14, info = typ)
    expect_true(max(a$numerical_certificate$d20$worst_ratios) < 1e-7, info = typ)
    expect_true(max(a$overall$numerical_certificate$d20$worst_ratios) < 1e-7, info = typ)
    expect_true(max(a$fixed_share_numerical_certificate$d20$worst_ratios) < 1e-7, info = typ)
  }
  expect_false(exists(".Random.seed", envir = globalenv()))
  # Known consequence recorded for the author (verification/review.md): the sup-t
  # critical value's certified enclosure on this fixture has relative width about
  # 1e-6 (built from the raw-draw and IQR-scale boxes), so band endpoints carry
  # about 1.7e-6 of an SE under bstrap = TRUE, cband = TRUE and are refused at the
  # default level; the IQR standard errors themselves pass. An explicit recorded
  # override accepts the same computation; nothing is relaxed automatically.
  set.seed(3)
  e <- expect_error(aggregate_chained_ipw_targets(rep_$estimates, rep_$if_full, rep_$target_meta, G, fr$sample_weight, type = "group", biters = 40L,
    multiplier_budget = pub_budget(), inference_accuracy = pub_rel(), available = rep_$available, enclosure = enc, bstrap = TRUE, cband = TRUE, caller_rng = TRUE,
    calendar = z$calendar, cohort_grid = sort(unique(G))), class = "cdid_ipw_pipeline_d20_error")
  expect_true(names(which.max(e$diagnostics$worst_ratios)) %in% c("band_lower", "band_upper"))
  expect_true(e$diagnostics$worst_ratios[["se"]] <= 1e-6); expect_true(max(e$diagnostics$multiplier_ratios) <= 1e-6)
  expect_true(e$diagnostics$worst_ratios[["band_upper"]] < 1e-5)
  set.seed(3)
  ok <- aggregate_chained_ipw_targets(rep_$estimates, rep_$if_full, rep_$target_meta, G, fr$sample_weight, type = "group", biters = 40L,
    multiplier_budget = pub_budget(), inference_accuracy = pub_rel(1e-5), available = rep_$available, enclosure = enc, bstrap = TRUE, cband = TRUE, caller_rng = TRUE,
    calendar = z$calendar, cohort_grid = sort(unique(G)))
  expect_true(ok$numerical_certificate$d20$certified); expect_identical(ok$numerical_certificate$d20$level, 1e-5)
  expect_true(ok$inference$bootstrap); expect_true(ok$inference$simultaneous)
  # the bootstrap IQR SEs alone (cband = FALSE) pass the default level on this fixture
  set.seed(3)
  tf <- aggregate_chained_ipw_targets(rep_$estimates, rep_$if_full, rep_$target_meta, G, fr$sample_weight, type = "group", biters = 40L,
    multiplier_budget = pub_budget(), inference_accuracy = pub_rel(), available = rep_$available, enclosure = enc, bstrap = TRUE, cband = FALSE, caller_rng = TRUE,
    calendar = z$calendar, cohort_grid = sort(unique(G)))
  expect_true(tf$numerical_certificate$d20$certified); expect_true(tf$inference$bootstrap)
})

test_that("the pipeline resolves the D20 list through its control validator and the released flag laws are untouched", {
  d <- expand.grid(id = 1:6, period = 1:4); d$G <- c(0, 0, 3, 3, 4, 4)[d$id]; d$x <- c(0.5, -1, 1, 0, 2, -2)[d$id]; d$y <- d$id + d$period
  prep <- chained_request_prepare_inputs(d, NULL, "y", "period", "id", "G", ~x)
  pc <- chained_public_pipeline_control(chained_public_dimensions(prep, 50), chained_public_resolve_control(NULL))
  v <- validate_ipw_pipeline_control(pc$control, weighting = "optimal")
  expect_null(v$inference_accuracy)
  d20_ctrl <- pc$control; d20_ctrl$inference_accuracy <- pub_rel(); d20_ctrl$gmm_accuracy$method <- "certified_eigen_chol_kkt"
  expect_identical(validate_ipw_pipeline_control(d20_ctrl, weighting = "optimal")$inference_accuracy, pub_rel())
  expect_identical(ipw_pipeline_flag_law(), c(private = "private_bootstrap_iff_cband", explicit = "explicit_bstrap_and_cband",
    cells_false_true = "released_bstrap_false_cband_true_cells_analytic_pointwise_no_draws",
    aggte_bins_false_true = "released_bstrap_false_cband_true_aggte_bins_analytic_se_multiplier_sup_t_band", empty = "empty_selection_no_draws"))
  ft <- ipw_pipeline_resolve_bootstrap_flags(FALSE, TRUE, endpoint = "cells"); expect_false(ft$draws); expect_false(ft$simultaneous)
  fb <- ipw_pipeline_resolve_bootstrap_flags(FALSE, TRUE, endpoint = "aggte_bins"); expect_true(fb$draws); expect_true(fb$simultaneous); expect_false(fb$bootstrap)
})

test_that("D20's certified lower SE normalization discriminates an upper-endpoint mutant", {
  # Independent arithmetic: 0.1005/1 > 0.1, whereas 0.1005/1.01 < 0.1.
  # The other errors (SE <= .005, covariance zero) meet either scale. No RNG/fit.
  cert <- list(theta = list(lo = matrix(0.3995, 1L), hi = matrix(0.6005, 1L)),
    rows = list(lo = matrix(c(-1, 1), 2L), hi = matrix(c(-1, 1), 2L)),
    se = matrix(c(1, 1.01), 2L), vcov = list(lo = matrix(1.005^2), hi = matrix(1.005^2)),
    band = list(lo = matrix(0.3995, 1L), hi = matrix(0.6005, 1L)),
    band_upper = list(lo = matrix(0.3995, 1L), hi = matrix(0.6005, 1L)), deterministic = integer())
  args <- list(cert = cert, estimates = c(a = 0.5), se = c(a = 1.005), vcov = matrix(1.005^2), level = 0.1,
    err_est = 0.1005, err_se = 0.005, err_v = matrix(0), band_lower = 0.5, band_upper = 0.5)
  expect_gt(0.1005 / 1, 0.1); expect_lt(0.1005 / 1.01, 0.1)
  expect_error(do.call(ipw_pipeline_enforce_d20, args), class = "cdid_ipw_pipeline_d20_error")
  mutant <- ipw_pipeline_enforce_d20
  original <- paste(deparse(body(mutant)), collapse = "\n")
  changed <- sub("cert$se[1L, ]", "cert$se[2L, ]", original, fixed = TRUE)
  expect_false(identical(changed, original))
  body(mutant) <- parse(text = changed)[[1L]]
  expect_true(do.call(mutant, args)$d20$certified)
})


test_that("D20 accepts requested outputs without inventing a raw-intermediate cutoff", {
  cert <- list(theta = list(lo = matrix(0.5), hi = matrix(0.5)),
    rows = list(lo = matrix(c(-1,1),2L),hi = matrix(c(-1,1),2L)),
    se = matrix(c(1,1),2L),vcov = list(lo = matrix(1),hi = matrix(1)),
    band = list(lo = matrix(0.5),hi = matrix(0.5)),band_upper = list(lo = matrix(0.5),hi = matrix(0.5)),
    deterministic = integer(),multiplier_error_by_target = 1e-5,raw_draw_scale_lower = 1)
  z <- ipw_pipeline_enforce_d20(cert,c(a=0.5),1,matrix(1),1e-6,0,0,matrix(0),0.5,0.5)
  expect_true(z$d20$certified)
  expect_gt(z$d20$multiplier_ratios[1L],1e-6)
  expect_false("multiplier" %in% z$d20$acceptance_families)
  expect_false(z$d20$intermediate_diagnostics$multiplier_ratio_is_acceptance_cutoff)
})


test_that("D20 needs authoritative reference role/operator plus exact boxes, while private zeros stay unchanged", {
  X <- matrix(0,4L,1L,dimnames=list(NULL,"nonreference"))
  E <- pub_exact_enclosure(0,X)
  run <- function(meta=NULL,acc=pub_rel(),enc=E) compute_chained_ipw_inference(c(nonreference=0),X,NULL,.05,FALSE,32L,
    multiplier_budget=pub_budget(),inference_accuracy=acc,enclosure=enc,bstrap=FALSE,caller_rng=TRUE,target_meta=meta)
  if(exists(".Random.seed",globalenv(),inherits=FALSE)) rm(".Random.seed",envir=globalenv())
  expect_error(run(),class="cdid_ipw_pipeline_d20_error")
  post <- pub_reference_meta("nonreference",integer())
  expect_error(run(post),class="cdid_ipw_pipeline_d20_error")
  known <- pub_reference_meta("nonreference",1L) # name is deliberately not evidence of role
  good <- run(known)
  expect_true(good$numerical_certificate$d20$certified)
  expect_identical(good$numerical_certificate$structural_reference,1L)
  expect_identical(good$numerical_certificate$d20$structural_zero_targets,"nonreference")
  wrong_operator <- known;wrong_operator$time <- 2
  expect_error(run(wrong_operator),class="cdid_ipw_pipeline_d20_error")
  missing_operator <- known;missing_operator$support_status <- "estimable"
  expect_error(run(missing_operator),class="cdid_ipw_pipeline_d20_error")
  wide <- E;wide$theta$lo[1] <- -1e-12;wide$theta$hi[1] <- 1e-12
  expect_error(run(known,enc=wide),class="cdid_ipw_pipeline_d20_error")
  legacy <- run(acc=pub_abs())
  expect_identical(legacy$se,c(nonreference=0));expect_null(legacy$numerical_certificate$d20)
  expect_null(legacy$numerical_certificate$structural_reference)
  expect_identical(legacy$numerical_certificate$deterministic,1L)
  expect_false(exists(".Random.seed",globalenv(),inherits=FALSE))
  # The declared support selector carries reference metadata by row, not name.
  EE <- list(theta=list(lo=matrix(c(NA,0),1L),hi=matrix(c(NA,0),1L)),
    rows=list(lo=cbind(NA_real_,X),hi=cbind(NA_real_,X)))
  tm <- pub_reference_meta(c("missing","nonreference"),2L)
  zz <- ipw_pipeline_supported_inference(c(missing=NA_real_,nonreference=0),cbind(NA_real_,X),c(FALSE,TRUE),
    multiplier_budget=pub_budget(),inference_accuracy=pub_rel(),enclosure=EE,bstrap=FALSE,target_meta=tm)
  expect_identical(zz$numerical_certificate$structural_reference,2L)
  expect_identical(zz$numerical_certificate$d20$structural_zero_targets,"nonreference")
  expect_true(is.na(zz$se[1L]));expect_identical(unname(zz$se[2L]),0)
})

test_that("dynamic reference composition requires nonempty known source identities and exact input boxes", {
  X <- cbind(r3=0,p3=pub_fixed_if()[,1L],r4=0,p4=pub_fixed_if()[,2L])
  th <- c(r3=0,p3=.5,r4=0,p4=1)
  tm <- data.frame(target_id=names(th),g=c(3,3,4,4),time=c(2,3,3,4),ref=c(2,2,3,3),base_t=c(2,2,3,3),
    role=c("reference","post","reference","post"),support_status=c("deterministic_reference_zero","estimable","deterministic_reference_zero","estimable"))
  E <- pub_exact_enclosure(th,X)
  run <- function(meta=tm,enc=E,acc=pub_rel()) aggregate_chained_ipw_targets(th,X,meta,rep(c(3,4),each=5L),rep(1,10L),
    type="dynamic",multiplier_budget=pub_budget(),inference_accuracy=acc,enclosure=enc,bstrap=FALSE,cband=FALSE)
  if(exists(".Random.seed",globalenv(),inherits=FALSE)) rm(".Random.seed",envir=globalenv())
  a <- run()
  # Independent identity: every source at horizon -1 is Y_ref-Y_ref; arbitrary
  # cohort weights and their derivatives therefore multiply exact zero.
  expect_identical(unname(a$estimates[1L]),0);expect_identical(unname(a$se[1L]),0)
  expect_identical(a$cell_keys,c(-1,0))
  expect_identical(a$cell_used_target_ids[[1L]],c("r3","r4"))
  expect_identical(a$numerical_certificate$structural_reference,1L)
  expect_identical(a$numerical_certificate$d20$structural_zero_targets,"event_-1")
  expect_true(all(a$numerical_certificate$theta$lo[1L]==0));expect_true(all(a$numerical_certificate$rows$hi[,1L]==0))
  expect_true(a$numerical_certificate$d20$certified);expect_true(a$fixed_share_numerical_certificate$d20$certified)
  not_ref <- tm;not_ref$role[c(1L,3L)] <- "pre"
  expect_error(run(not_ref),class="cdid_ipw_pipeline_d20_error")
  one_not_ref <- tm;one_not_ref$role[3L] <- "pre"
  expect_error(run(one_not_ref),class="cdid_ipw_pipeline_d20_error")
  wider <- E;wider$rows$lo[,1L] <- -1e-12;wider$rows$hi[,1L] <- 1e-12
  expect_error(run(enc=wider),class="cdid_ipw_pipeline_d20_error")
  legacy <- run(not_ref,acc=pub_abs())
  expect_identical(legacy$estimates,a$estimates);expect_identical(legacy$se,a$se)
  expect_null(legacy$numerical_certificate$structural_reference)
  expect_false(exists(".Random.seed",globalenv(),inherits=FALSE))
  # No empty source family can prove an aggregate identity.
  fake <- pub_reference_meta("x",integer());fake$role <- "reference_aggregate"
  fake$support_status <- "deterministic_reference_zero";fake$reference_operator <- "weighted_reference_identity"
  fake$reference_sources <- list(tm[FALSE,,drop=FALSE])
  expect_error(compute_chained_ipw_inference(c(x=0),matrix(0,4L,1L),multiplier_budget=pub_budget(),
    inference_accuracy=pub_rel(),enclosure=pub_exact_enclosure(0,matrix(0,4L,1L)),bstrap=FALSE,target_meta=fake),
    class="cdid_ipw_pipeline_d20_error")
})
