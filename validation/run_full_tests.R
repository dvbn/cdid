#!/usr/bin/env Rscript
# Execute every testthat file against an installed cdid build. The native
# DLL lifecycle tests run separately from fits that hold native finalizers.
#
# Standard usage after R CMD check:
#   Rscript validation/run_full_tests.R \
#     --lib-loc=check/cdid.Rcheck \
#     --test-dir=check/cdid.Rcheck/tests/testthat \
#     --output-dir=full-tests
#
# Standard usage against an already installed cdid in .libPaths():
#   Rscript validation/run_full_tests.R

command <- commandArgs()
args <- commandArgs(trailingOnly = TRUE)
native_only <- "--native-only" %in% args
native_file <- "test-cp02-native-overhead.R"
native_filter <- "^cp02-native-overhead$"

parse_flag <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[length(hit)])
}

default_lib <- if (dir.exists(file.path("check", "cdid.Rcheck", "cdid"))) {
  file.path("check", "cdid.Rcheck")
} else {
  NULL
}
default_test_dir <- if (dir.exists(file.path("check", "cdid.Rcheck", "tests", "testthat"))) {
  file.path("check", "cdid.Rcheck", "tests", "testthat")
} else {
  file.path("tests", "testthat")
}

lib_loc <- parse_flag("--lib-loc", default_lib)
test_dir <- parse_flag("--test-dir", default_test_dir)
output_dir <- parse_flag("--output-dir", "full-tests")

if (!is.null(lib_loc)) {
  if (!nzchar(lib_loc)) stop("--lib-loc must name a library containing cdid.")
  lib_loc <- normalizePath(lib_loc, winslash = "/", mustWork = TRUE)
  stopifnot(dir.exists(file.path(lib_loc, "cdid")))
  .libPaths(unique(c(lib_loc, .libPaths())))
}
test_dir <- normalizePath(test_dir, winslash = "/", mustWork = TRUE)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

library(testthat)
if (is.null(lib_loc)) library(cdid) else library(cdid, lib.loc = lib_loc)
package_path <- normalizePath(getNamespaceInfo(asNamespace("cdid"), "path"), winslash = "/", mustWork = TRUE)
if (!is.null(lib_loc)) {
  stopifnot(identical(package_path,
                      normalizePath(file.path(lib_loc, "cdid"), winslash = "/", mustWork = TRUE)))
}

# Match testthat's file discovery and verify that the complete archive contains
# the primary estimator and nuisance comparison suites.
test_files <- sort(list.files(test_dir, pattern = "^test.*\\.[rR]$", full.names = FALSE))
required_test_files <- c("test-cran-regression.R", "test-cran-regression-static.R",
                         "test-nuisance-models.R", "test-cross-pair-chaining.R",
                         "test-cohort-time-fe.R", "test-cohort-time-nuisance-scores.R",
                         "test-static-nuisance-capture.R",
                         "test-chained-missing-reference.R")
expected_test_count <- 191L
if (length(test_files) != expected_test_count || !all(required_test_files %in% test_files)) {
  stop(sprintf("Expected all %d packaged testthat files in '%s', found %d.",
               expected_test_count, test_dir, length(test_files)))
}
stopifnot(native_file %in% test_files)
if (native_only) test_files <- native_file

log_file <- file.path(output_dir, "full-testthat.log")
summary_file <- file.path(output_dir, "full-test-summary.txt")
results_rds <- file.path(output_dir, "full-test-results.rds")

list_rep <- testthat::ListReporter$new()
summary_rep <- testthat::SummaryReporter$new(show_praise = FALSE, file = log_file)
multi_rep <- testthat::MultiReporter$new(reporters = list(summary_rep, list_rep))

started_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
t0 <- proc.time()
run_error <- NULL
native_error <- NULL
native_results <- NULL

if (!native_only) {
  self <- sub("^--file=", "", grep("^--file=", command, value = TRUE))
  stopifnot(length(self) == 1L, file.exists(self))
  native_dir <- tempfile("native-overhead-", tmpdir = output_dir)
  native_args <- c("--vanilla", shQuote(self), "--native-only",
    shQuote(paste0("--lib-loc=", dirname(package_path))),
    shQuote(paste0("--test-dir=", test_dir)),
    shQuote(paste0("--output-dir=", native_dir)))
  native_log <- file.path(output_dir, "native-overhead-console.log")
  native_status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    native_args, stdout = native_log, stderr = native_log))
  native_rds <- file.path(native_dir, "full-test-results.rds")
  result <- if (file.exists(native_rds)) tryCatch(readRDS(native_rds), error = identity) else NULL
  valid_result <- is.data.frame(result) && identical(unique(result$file), native_file)
  if (valid_result) native_results <- result
  if (!identical(native_status, 0L) || !valid_result) {
    native_error <- sprintf("Native DLL test process failed (status %s). See %s.",
      native_status, native_log)
  }
}

tryCatch(
  {
    testthat::test_dir(
      path = test_dir,
      filter = native_filter,
      invert = !native_only,
      package = "cdid",
      load_package = "installed",
      reporter = multi_rep,
      stop_on_failure = TRUE,
      stop_on_warning = TRUE
    )
  },
  error = function(e) {
    run_error <<- e
  }
)

elapsed_sec <- as.numeric((proc.time() - t0)[["elapsed"]])
finished_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

res_list <- list_rep$get_results()
res_df <- as.data.frame(res_list)
if (!is.null(native_results)) res_df <- rbind(res_df, native_results)
saveRDS(res_df, file = results_rds)

# ListReporter has one row per test block; file is the basename passed to
# start_file(), not the context label. Compare identities, not just counts.
executed_files <- sort(unique(res_df$file), na.last = TRUE)
missing_files <- setdiff(test_files, executed_files)
unexpected_files <- setdiff(executed_files, test_files)
n_files_executed <- length(executed_files)
n_passed <- if (nrow(res_df)) sum(res_df$passed) else 0L
n_failed <- if (nrow(res_df)) sum(res_df$failed) else 0L
n_warnings <- if (nrow(res_df)) sum(res_df$warning) else 0L
n_skipped <- if (nrow(res_df)) sum(as.integer(res_df$skipped)) else 0L
n_errors <- if (nrow(res_df)) sum(as.integer(res_df$error)) else 0L

status_ok <- is.null(run_error) && is.null(native_error) &&
  n_failed == 0L &&
  n_warnings == 0L &&
  n_errors == 0L &&
  !anyNA(executed_files) &&
  setequal(executed_files, test_files)

summary_lines <- c(
  "cdid Unfiltered Full Installed Test Suite Summary",
  sprintf("status: %s", if (status_ok) "PASS" else "FAIL"),
  sprintf("package_version: %s", as.character(getNamespaceVersion("cdid"))),
  sprintf("package_path: %s", package_path),
  sprintf("requested_library: %s", if (is.null(lib_loc)) "default .libPaths()" else lib_loc),
  sprintf("test_dir: %s", test_dir),
  if (native_only) "scope: isolated native DLL tests" else "scope: all files across two processes",
  sprintf("started_utc: %s", started_utc),
  sprintf("finished_utc: %s", finished_utc),
  sprintf("elapsed_seconds: %.2f", elapsed_sec),
  sprintf("elapsed_minutes: %.2f", elapsed_sec / 60),
  sprintf("packaged_test_files: %d", length(test_files)),
  sprintf("executed_test_files: %d", n_files_executed),
  sprintf("missing_test_files: %s", paste(missing_files, collapse = ", ")),
  sprintf("unexpected_test_files: %s", paste(unexpected_files, collapse = ", ")),
  sprintf("passed_assertions: %d", n_passed),
  sprintf("failed_assertions: %d", n_failed),
  sprintf("warnings: %d", n_warnings),
  sprintf("errors: %d", n_errors),
  sprintf("skipped_tests: %d", n_skipped),
  sprintf("log_file: %s", log_file)
)
if (!is.null(native_error)) summary_lines <- c(summary_lines, native_error)
if (!is.null(run_error)) {
  summary_lines <- c(summary_lines, sprintf("error_message: %s", conditionMessage(run_error)))
}

writeLines(summary_lines, con = summary_file)
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")

if (!status_ok) {
  if (!is.null(native_error)) stop(native_error, call. = FALSE)
  if (!is.null(run_error)) {
    stop(run_error)
  }
  stop(sprintf(
    "Unfiltered test suite failed: failed=%d, warnings=%d, errors=%d, executed_files=%d/%d.",
    n_failed, n_warnings, n_errors, n_files_executed, length(test_files)
  ))
}
