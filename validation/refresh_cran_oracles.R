#!/usr/bin/env Rscript

# Run from the source root with the intended cdid build installed.
library(cdid)
library(testthat)

references <- list()
capture_env <- new.env(parent = asNamespace("cdid"))
capture_env$cran_regression_oracle <- function(case, compute) {
  stopifnot(!case %in% names(references))
  reference <- force(compute)
  references[[case]] <<- reference
  reference
}

testthat::test_dir("tests/testthat", filter = "^cran-regression$",
                   env = capture_env, package = "cdid", load_package = "installed",
                   stop_on_failure = TRUE, stop_on_warning = TRUE)
stopifnot(length(references) == 11L)
saveRDS(references, "tests/testthat/fixtures/cran-regression-oracles.rds",
        version = 2, compress = "xz")
print(sessionInfo())
