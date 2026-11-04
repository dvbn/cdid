# Public composition and support edges; no new data-generating experiment.
oped_ns <- function(name) getFromNamespace(name, "cdid")
oped_fixture <- function(dates = 1:4) {
  e <- new.env(parent = as.environment("package:stats"))
  sys.source(testthat::test_path("reference", "imputation_public_reference.R"), e)
  raw <- utils::read.csv(testthat::test_path("fixtures", "aug_public_v1", "raw.csv"))
  fx <- e$imputation_public_fixture(raw, dates)
  fx$args$outcome_pooling <- NULL
  fx$args$selection <- "none"; fx$args$selection_est_method <- "ipw"
  fx
}
oped_fit <- function(fx, overrides = list()) do.call(cdid::att_gt, utils::modifyList(fx$args, overrides))

testthat::test_that("public separate routing covers both control laws and PT systems for OR and DR", {
  fx <- oped_fixture()
  cases <- expand.grid(control_group = c("nevertreated", "notyettreated"),
    pt_assumption = c("post", "all"), est_method = c("reg", "dr"), stringsAsFactors = FALSE)
  # Both reporting base families appear without multiplying unrelated axes.
  cases$base_period <- ifelse(cases$pt_assumption == "all", "universal", "varying")
  for (i in seq_len(nrow(cases))) {
    cs <- as.list(cases[i, ]); if (cs$pt_assumption == "all") cs$pt_start <- 1
    a <- oped_fit(fx, cs)
    b <- oped_fit(fx, c(cs, list(outcome_pooling = "separate", control = list(chunk_rows = 7L))))
    info <- paste(unlist(cs), collapse = "/")
    testthat::expect_equal(a$att, b$att, tolerance = 1e-9, info = info)
    testthat::expect_equal(a$inffunc, b$inffunc, tolerance = 1e-9, info = info)
    testthat::expect_equal(a$V_analytical, b$V_analytical, tolerance = 1e-9, info = info)
    testthat::expect_identical(a$chained$available, b$chained$available)
    testthat::expect_identical(a$chained$request$outcome_pooling, "separate")
    testthat::expect_identical(a$chained$request$control_group, cs$control_group)
    testthat::expect_identical(a$chained$request$pt_assumption, cs$pt_assumption)
    testthat::expect_identical(a$DIDparams$base_period, cs$base_period)
    testthat::expect_identical(rownames(a$inffunc), fx$frame$id)
    good <- which(a$chained$available)
    X <- a$inffunc[, good, drop = FALSE]
    X <- sweep(X, 2L, colMeans(X), "-")
    testthat::expect_equal(unname(a$V_analytical[good, good, drop = FALSE]), unname(crossprod(X) / nrow(X)),
      tolerance = 1e-8, info = info)
    nd <- a$chained$nuisance_diagnostics
    testthat::expect_equal(nd$accounting$n_regression_fits, nd$accounting$n_regression_blocks)
    testthat::expect_equal(length(nd$outcome_blocks), nd$accounting$n_regression_blocks)
    if (cs$base_period == "universal") {
      reference <- which(a$chained$role == "reference")
      testthat::expect_gt(length(reference), 0L)
      testthat::expect_true(all(a$att[reference] == 0))
    }
  }
})

testthat::test_that("separate dated recipes retain weight anchors and chunk-invariant original-unit rows", {
  d <- utils::read.csv(testthat::test_path("fixtures", "cpt", "fBd.csv"))
  d$w <- d$w * (1 + d$time / 10)
  for (cs in list(list(est_method = "reg", fix_weights = "first_period"),
    list(est_method = "dr", fix_weights = "base_period"))) {
    args <- c(list(yname = "y", tname = "time", idname = "id", gname = "G", data = d,
      xformla = ~x + z, weightsname = "w", control_group = "notyettreated",
      bstrap = FALSE, cband = FALSE, gmm_weighting = "identity"), cs)
    a <- do.call(cdid::att_gt, args)
    b <- do.call(cdid::att_gt, c(args, list(outcome_pooling = "separate", control = list(chunk_rows = 7L))))
    testthat::expect_equal(a$att, b$att, tolerance = 1e-9)
    testthat::expect_equal(a$inffunc, b$inffunc, tolerance = 1e-9)
    testthat::expect_identical(nrow(a$inffunc), length(unique(d$id)))
    ref <- a$chained$reference_inputs
    testthat::expect_identical(ref$reduction, "reference_recipes")
    testthat::expect_identical(a$DIDparams$outcome_pooling, "separate")
    rec <- ref$recipes$recipes
    testthat::expect_true(all(is.finite(rec$weight_anchor)))
    if (cs$fix_weights == "first_period") testthat::expect_true(all(rec$weight_anchor == min(d$time)))
    testthat::expect_equal(length(a$chained$nuisance_diagnostics$basis$blocks),
      a$chained$nuisance_diagnostics$accounting$n_regression_blocks)
  }
})

testthat::test_that("separate local aliases preserve supported predictions and refuse unsupported predictions", {
  fx <- oped_fixture(1:3); fx$args$est_method <- "reg"
  fx$args$frame$local_alias <- ifelse(fx$frame$g == 4, fx$frame$x^2 + 1, 2 * fx$frame$x)
  base <- oped_fit(fx)
  alias <- oped_fit(fx, list(xformla = ~x + local_alias))
  testthat::expect_equal(alias$att, base$att, tolerance = 1e-9)
  testthat::expect_equal(alias$inffunc, base$inffunc, tolerance = 1e-9)
  testthat::expect_true(any(vapply(alias$chained$nuisance_diagnostics$basis$blocks,
    function(b) "local_alias" %in% b$aliased_columns, logical(1L))))
  bad <- oped_fixture(); bad$args$est_method <- "reg"
  bad$args$frame$bad <- ifelse(bad$frame$g == 0, 2 * bad$frame$x, bad$frame$x^2 + 1)
  testthat::expect_error(oped_fit(bad, list(xformla = ~x + bad)), class = "cdid_or_wls_rank_error")
  # Empty outcome designs retain their existing refusal, including separate mode.
  testthat::expect_error(oped_fit(bad, list(xformla = ~0)), class = "cdid_chained_request_error")
})

testthat::test_that("separate public intercept-only DR equals IPW and no-intercept OR retains its formula", {
  fx <- oped_fixture()
  dr <- oped_fit(fx, list(est_method = "dr", xformla = ~1))
  ipw <- oped_fit(fx, list(est_method = "ipw", xformla = ~1, outcome_pooling = quote(stop("inactive outcome choice"))))
  testthat::expect_equal(dr$att, ipw$att, tolerance = 1e-9)
  testthat::expect_equal(dr$inffunc, ipw$inffunc, tolerance = 1e-9)
  testthat::expect_equal(dr$V_analytical, ipw$V_analytical, tolerance = 1e-9)
  no_intercept <- oped_fit(fx, list(est_method = "reg", xformla = ~0 + x))
  nd <- no_intercept$chained$nuisance_diagnostics
  testthat::expect_false(nd$coefficients$has_local_intercepts)
  testthat::expect_true(all(nd$coefficients$local_intercepts == 0))
  testthat::expect_identical(no_intercept$DIDparams$outcome_pooling, "separate")
})

testthat::test_that("empty post catalogs preserve separate auxiliary fits and public request metadata", {
  fx <- oped_fixture(); fx$args$est_method <- "reg"
  post <- fx$data$g > 0 & fx$data$time >= fx$data$g
  fx$args$data$y[post] <- NA_real_
  mp <- oped_fit(fx)
  testthat::expect_identical(mp$DIDparams$outcome_pooling, "separate")
  testthat::expect_identical(mp$chained$request$outcome_pooling, "separate")
  testthat::expect_true(all(is.na(mp$att[mp$t >= mp$group])))
  diag <- mp$chained$nuisance_diagnostics
  testthat::expect_gt(length(diag$outcome_blocks), 0L)
  testthat::expect_true(all(diag$fitted_catalog$role == "auxiliary"))
  args <- fx$args
  for (nm in c("imputation_pooling", "gmm_weighting")) args[[nm]] <- NULL
  args$biters <- 9L
  test <- do.call(cdid::conditional_did_pretest, args)
  testthat::expect_identical(test$chained$status, "available")
  testthat::expect_identical(test$chained$request$outcome_pooling, "separate")
  testthat::expect_true(is.finite(test$CvM))
})
