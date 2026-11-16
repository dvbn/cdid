source_commit <- "956537a5e05ff6d727473eea5edd6650be2fa711"
archive_sha <- "41da46bf563a25ae28deeef9595a23037c094dbbef125319fc410fff019cb644"
archive <- file.path("source-archive", "cdid_1.3.0.0.tar.gz")
stopifnot(identical(readLines(file.path("source-archive", "SOURCE_COMMIT")), source_commit),
          identical(digest::digest(file = archive, algo = "sha256"), archive_sha))
dir.create("windows-diagnostic")
report <- normalizePath("windows-diagnostic", winslash = "/", mustWork = TRUE)
Sys.setenv(CDID_WINDOWS_PROBE_REPORT = report)
lib <- tempfile("cdid-windows-lib-")
dir.create(lib)
check <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", paste0("--library=", shQuote(lib)), shQuote(archive)),
  stdout = file.path(report, "install.log"), stderr = file.path(report, "install.log")
)
stopifnot(identical(check, 0L))
.libPaths(c(lib, .libPaths()))
library(cdid)
stopifnot(packageVersion("cdid") == "1.3.0.0",
          identical(normalizePath(dirname(find.package("cdid"))), normalizePath(lib)))
writeLines(c(paste("Source commit", source_commit), paste("Archive SHA256", archive_sha),
             capture.output(sessionInfo())), file.path(report, "environment.txt"))

# Capture the original linker's untruncated output from the published helper.
original <- system2("git", c("show", paste0(source_commit,
  ":tests/testthat/test-cp02-native-overhead.R")), stdout = TRUE)
helper <- new.env()
for (expr in parse(text = original)) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name("noh_compile_dll"))) eval(expr, helper)
}
before <- tryCatch(helper$noh_compile_dll(readLines("src/cdid.c"), "cdid",
  tempdir(), source_dir = normalizePath("src")), error = identity)
writeLines(if (inherits(before, "error")) conditionMessage(before) else "Original test DLL linked.",
           file.path(report, "original-link.txt"))

insert_probe <- function(path, marker, lines) {
  source <- readLines(path)
  at <- which(trimws(source) == marker)
  stopifnot(length(at) == 1L)
  writeLines(append(source, lines, after = at - 1L), path)
}

# Measure the operands and original-unit coordinates before the unchanged test.
insert_probe("tests/testthat/test-augmented-replay-operation.R",
  "testthat::expect_identical(actual == 0, expected == 0)", c(
    '        mismatches <- which(xor(actual == 0, expected == 0))',
    '        if (length(mismatches)) {',
    '          rows <- (mismatches - 1L) %% NROW(actual) + 1L',
    '          cols <- (mismatches - 1L) %/% NROW(actual) + 1L',
    '          keys <- rownames(actual)',
    '          vals <- data.frame(index = mismatches, row = rows, column = cols,',
    '            row_key = if (is.null(keys)) rows else keys[rows],',
    '            actual_hex = sprintf("%a", as.numeric(actual[mismatches])),',
    '            expected_hex = sprintf("%a", as.numeric(expected[mismatches])))',
    '          if (identical(nm, "inffunc")) {',
    '            vals$cohort <- fast$group[cols]',
    '            vals$date <- fast$t[cols]',
    '          }',
    '          write.csv(vals, file.path(Sys.getenv("CDID_WINDOWS_PROBE_REPORT"),',
    '            paste0("zero-mask-", nm, ".csv")), row.names = FALSE)',
    '        }'
  )
)
pretest <- "tests/testthat/test-chained-conditional-pretest.R"
insert_probe(pretest, "set.seed(11); r1 <- call_it(); u1 <- stats::runif(1L)", c(
  '  serial_samples <- list()',
  '  worker_samples <- list()',
  '  native_batch <- getFromNamespace("cdid_fast_cpt_batch_cvmb", "cdid")',
  '  worker_apply <- parallel::parLapply',
  '  testthat::local_mocked_bindings(cdid_fast_cpt_batch_cvmb = function(E, cells, N) {',
  '    value <- native_batch(E, cells, N)',
  '    serial_samples[[length(serial_samples) + 1L]] <<- list(E = E, value = value)',
  '    value',
  '  }, .package = "cdid")',
  '  testthat::local_mocked_bindings(parLapply = function(cl, X, fun, ...) {',
  '    value <- worker_apply(cl, X, fun, ...)',
  '    worker_samples[[length(worker_samples) + 1L]] <<- list(signs = X, value = value, args = list(...))',
  '    value',
  '  }, .package = "parallel")'
))
insert_probe(pretest, 'expect_identical(r1$CvMb, r2$CvMb)', c(
  '  serial_batches <- ceiling(length(r1$CvMb) / r1$chained$inference$rng$batch_rows)',
  '  serial_signs <- lapply(serial_samples[seq_len(serial_batches)], `[[`, "E")',
  '  serial_seed <- .Random.seed'
))
insert_probe(pretest, 'expect_identical(r6$CvMb, r1$CvMb)', c(
  '  probe_dir <- Sys.getenv("CDID_WINDOWS_PROBE_REPORT")',
  '  worker_signs <- do.call(c, lapply(worker_samples, `[[`, "signs"))',
  '  pure_batch <- getFromNamespace("cpt_batch_cvmb", "cdid")',
  '  local_worker_draws <- unlist(lapply(worker_samples, function(s)',
  '    lapply(s$signs, function(E) do.call(pure_batch, c(list(E = E), s$args)))))',
  '  stopifnot(length(r1$CvMb) == length(r6$CvMb))',
  '  write.csv(data.frame(draw = seq_along(r1$CvMb), serial_hex = sprintf("%a", r1$CvMb),',
  '    worker_hex = sprintf("%a", r6$CvMb), delta = r6$CvMb - r1$CvMb),',
  '    file.path(probe_dir, "pretest-draws.csv"), row.names = FALSE)',
  '  writeLines(c(sprintf("sign_batches_identical: %s", identical(serial_signs, worker_signs)),',
  '    sprintf("rng_state_identical: %s", identical(serial_seed, .Random.seed)),',
  '    sprintf("max_absolute_draw_delta: %.17g", max(abs(r6$CvMb - r1$CvMb))),',
  '    sprintf("R_caller_matches_PSOCK: %s", identical(local_worker_draws, r6$CvMb)),',
  '    sprintf("CvM: %s %s", format(r1$CvM, digits = 17), format(r6$CvM, digits = 17)),',
  '    sprintf("pvalue: %s %s", format(r1$CvMpval, digits = 17), format(r6$CvMpval, digits = 17)),',
  '    sprintf("critical: %s %s", format(r1$CvMcval, digits = 17), format(r6$CvMcval, digits = 17)),',
  '    sprintf("exceedances_identical: %s", identical(r1$CvMb >= r1$CvM, r6$CvMb >= r6$CvM))),',
  '    file.path(probe_dir, "pretest-diagnostics.txt"))'
))

reporter <- testthat::ListReporter$new()
tests <- testthat::test_dir("tests/testthat",
  filter = "^(augmented-replay-operation|chained-conditional-pretest|cp02-native-overhead)$",
  package = "cdid", load_package = "installed",
  reporter = testthat::MultiReporter$new(reporters = list(
    reporter, testthat::SummaryReporter$new(file = file.path(report, "tests.log")))),
  stop_on_failure = FALSE, stop_on_warning = FALSE
)
rows <- as.data.frame(tests)
cols <- c("file", "test", "failed", "error", "warning", "skipped", "passed")
write.csv(rows[cols], file.path(report, "tests.csv"), row.names = FALSE)
stopifnot(length(unique(rows$file)) == 3L)
print(aggregate(cbind(failed, error, warning, skipped, passed) ~ file, rows, sum))
if (any(rows$failed > 0L | rows$error | rows$warning > 0L)) quit(status = 1L)
