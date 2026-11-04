# Public composition checks reuse the accepted finite fixture builder. The
# independent numerical references and complete IF tests remain private tests.
opub_ns <- function(name) getFromNamespace(name, "cdid")
opub_fixture <- function() {
  env <- new.env(parent = as.environment("package:stats"))
  sys.source(testthat::test_path("reference", "imputation_public_reference.R"), env)
  raw <- utils::read.csv(testthat::test_path("fixtures", "aug_public_v1", "raw.csv"))
  args <- env$imputation_public_fixture(raw)$args
  args$outcome_pooling <- NULL
  args
}

testthat::test_that("outcome choice is appended and inactive promises remain lazy", {
  testthat::expect_identical(tail(names(formals(cdid::conditional_did_pretest)), 2L), c("control", "outcome_pooling"))
  for (method in c("reg", "dr")) testthat::expect_error(cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g", data = stop("forced data"),
    est_method = method, outcome_pooling = "invalid"), class = "cdid_chained_request_outcome_pooling")
  for (method in list("ipw", function(...) stop("callback must not run")))
    testthat::expect_error(cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g",
      data = stop("forced data"), est_method = method, outcome_pooling = stop("forced outcome")), "forced data")
  testthat::expect_error(cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g",
    data = stop("ordinary data"), mode = "ordinary", outcome_pooling = stop("forced outcome")), "ordinary data")
  cpt <- cdid::conditional_did_pretest(yname = "y", tname = "t", idname = "id", gname = "g",
    data = stop("custom data"), est_method = function(...) NULL, outcome_pooling = stop("forced outcome"))
  testthat::expect_identical(cpt$chained$status, "not_applicable")
})

testthat::test_that("public default explicitly routes separate coefficients across accepted families", {
  routes <- list(list(est_method = "reg", selection = "none", selection_est_method = "ipw"),
    list(est_method = "dr", selection = "none", selection_est_method = "ipw", control_group = "notyettreated"),
    list(est_method = "reg", selection_est_method = "ipw", selection_pooling = "within_cohort"),
    list(est_method = "dr", selection_est_method = "ipw", selection_model = "independent_period_logit", selection_pooling = "within_cohort"),
    list(est_method = "dr", selection_est_method = "imputation", imputation_pooling = "within_cohort"),
    list(est_method = "dr", selection_est_method = "dr", imputation_pooling = "across_cohorts", selection_pooling = "within_cohort"))
  for (route in routes) {
    args <- utils::modifyList(opub_fixture(), route)
    omitted <- do.call(cdid::att_gt, args)
    separate <- do.call(cdid::att_gt, c(args, list(outcome_pooling = "separate")))
    pooled <- do.call(cdid::att_gt, c(args, list(outcome_pooling = "pooled")))
    testthat::expect_equal(omitted$att, separate$att, tolerance = 1e-10)
    testthat::expect_equal(omitted$inffunc, separate$inffunc, tolerance = 1e-10)
    testthat::expect_identical(omitted$DIDparams$outcome_pooling, "separate")
    testthat::expect_identical(omitted$chained$request$outcome_pooling, "separate")
    testthat::expect_identical(pooled$chained$request$outcome_pooling, "pooled")
    diag <- omitted$chained$nuisance_diagnostics
    testthat::expect_identical(diag$outcome_pooling, "separate")
    testthat::expect_null(diag$coefficients$shared_slopes)
    testthat::expect_gt(length(diag$coefficients$block_slopes), 0L)
  }
})

testthat::test_that("D11 clustering leaves outcome points and unit IF rows unchanged", {
  args <- utils::modifyList(opub_fixture(), list(selection = "none", selection_est_method = "ipw",
    est_method = "dr", gmm_weighting = "optimal"))
  unit <- do.call(cdid::att_gt, args)
  cluster <- do.call(cdid::att_gt, c(args, list(clustervars = "cl")))
  testthat::expect_equal(unit$att, cluster$att, tolerance = 1e-10)
  testthat::expect_equal(unit$inffunc, cluster$inffunc, tolerance = 1e-10)
})
