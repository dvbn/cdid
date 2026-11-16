library(testthat)
library(cdid)

# Ordinary R CMD check always executes the fixed bounded regression filter in
# tests/testthat/test-cran-regression.R and
# tests/testthat/test-cran-regression-static.R without environment-variable
# branching. The anchored union selects exactly those two files.
# Use saved independent outputs for the deterministic finite refits. The
# unfiltered suite recomputes them and runs the same comparison assertions.
references <- readRDS(file.path("testthat", "fixtures", "cran-regression-oracles.rds"))
check_env <- new.env(parent = asNamespace("cdid"))
check_env$cran_regression_oracle <- function(case, compute) {
  if (is.null(references[[case]])) stop("Missing CRAN regression oracle: ", case)
  references[[case]]
}
test_check("cdid", filter = "^(cran-regression|cran-regression-static)$", env = check_env)
