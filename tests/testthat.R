library(testthat)
library(cdid)

# Ordinary R CMD check always executes the fixed bounded regression filter in
# tests/testthat/test-cran-regression.R and
# tests/testthat/test-cran-regression-static.R without environment-variable
# branching. All 185 comprehensive cohort and static development and Monte Carlo
# test files remain packaged in tests/testthat/ and are run unfiltered
# (filter = NULL) in CI and via validation/run_full_tests.R.
test_check("cdid", filter = "cran-regression")
