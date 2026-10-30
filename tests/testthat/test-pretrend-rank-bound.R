test_that("known centered cluster rank bounds precede numerical Cholesky", {
  r <- readRDS(test_path("fixtures", "pretrend_cluster5_reporting.rds"))
  result <- chained_pretest_record(r)
  expect_identical(result$status, "unavailable")
  expect_identical(result$df, 6L)
  expect_equal(result$covariance_rank_upper, 4)
  expect_null(result$W)
  expect_null(result$Wpval)
  expect_match(result$reason, "known centered covariance rank bound")
  expect_length(result$cells, 6L)
  # The refusal is based on the count, including a roundoff-PD matrix.
  r$inference$vcov <- diag(length(r$estimates))
  expect_identical(chained_pretest_record(r)$status, "unavailable")
})

test_that("unit rank bounds and ordinary sufficient-row tests remain distinct", {
  r <- list(target_meta=data.frame(role=c('pre','pre'),target_id=c('a','b'),available=TRUE),
    estimates=c(.1,.2),base_period='varying',
    inference=list(vcov=diag(2),clustered=FALSE,n_total=2))
  expect_identical(chained_pretest_record(r)$status, "unavailable")
  r$inference$n_total <- 3
  result <- chained_pretest_record(r)
  expect_identical(result$status, "computed")
  expect_equal(result$W, .05)
  r$inference$n_total <- NULL
  expect_identical(chained_pretest_record(r)$status, "computed")
})
