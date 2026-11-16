# `cdid` Test Suites

`cdid` ships two levels of automated verification in `tests/testthat/`:

1. **Bounded `R CMD check` regression suite (`tests/testthat.R` -> `tests/testthat/test-cran-regression.R` and `tests/testthat/test-cran-regression-static.R`)**
   Standard `R CMD check` runs `testthat::test_check("cdid", filter = "^(cran-regression|cran-regression-static)$")` unconditionally without environment-variable branching; the anchored filter selects exactly these two files. This suite exercises all supported public estimators (`att_gt` with `gname` and `dname`, `aggte`, `conditional_did_pretest`), treatment adjustments (`"ipw"`, `"reg"`, `"dr"`, and custom IPW callbacks), missing-outcome adjustments (`selection = "none"`, `"conditional_mtar"` / `"corollary1"` with `"ipw"`, `"imputation"`, and `"dr"`), outcome regression poolings (`"separate"`, `"pooled"`, and static cell partitions), selection logit models (`"pair_logit"` and `"independent_period_logit"`), cohort and static pooling modes (`"separate"`, `"within_cohort"`, `"across_cohorts"`, `"within_group"`, `"pooled"`), static imputation scales (`"per_period"` and `"per_change"`), static transition selectors (`switchers` and `switch_times`), control groups (`"nevertreated"` and `"notyettreated"`), parallel-trends assumptions (`"post"` and `"all"`), GMM weighting and Tikhonov regularization, time-varying covariates and weight anchors (`fix_weights`), analytical and multiplier-bootstrap inference, ordinary-mode delegation to `did`, summary/reference citations, and typed refusal gates against independent finite-sample references and frozen oracle fixtures.

The bounded suite uses representative independent references, not an independent numerical reference for every option combination. Additional combinations and full nuisance-derivative checks remain in the mandatory complete suite.

For finite-difference checks on fixed input, the bounded driver reads saved output from the independent oracle. Every public fit and result comparison still runs. The same tests recompute the reference in the unfiltered suite. To regenerate the saved outputs after an intentional fixture or oracle change, install the source and run `Rscript validation/refresh_cran_oracles.R` from the checkout. The generator runs the live reference tests and writes only their independent oracle outputs after all assertions pass.

2. **Unfiltered complete installed test suite (`filter = NULL`)**
   The complete suite runs all 191 packaged test files, including the bounded tests and the full cohort, static switching, cohort-by-date, cross-pair, and missing-reference checks. From a GitHub checkout, use the convenience script to execute the entire unfiltered suite against an installed `cdid` build and the packaged test directory after `R CMD check`:

   ```sh
   Rscript validation/run_full_tests.R \
     --lib-loc=check/cdid.Rcheck \
     --test-dir=check/cdid.Rcheck/tests/testthat \
     --output-dir=full-tests
   ```

   The direct R command also works with the tests supplied in the CRAN source archive:

   ```r
   testthat::test_dir(
     "tests/testthat",
     filter = NULL,
     package = "cdid",
     load_package = "installed",
     stop_on_failure = TRUE,
     stop_on_warning = TRUE
   )
   ```

   Every four-platform CI job in `.github/workflows/R-CMD-check.yaml` runs both stages against the identical built source archive: first `R CMD check --as-cran` (with zero errors, warnings, or notes), and then the unfiltered complete installed `test_dir(..., filter = NULL)` gate.
