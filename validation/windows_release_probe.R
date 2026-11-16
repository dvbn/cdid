source_commit <- "956537a5e05ff6d727473eea5edd6650be2fa711"
archive_sha <- "41da46bf563a25ae28deeef9595a23037c094dbbef125319fc410fff019cb644"
archive <- file.path("source-archive", "cdid_1.3.0.0.tar.gz")
stopifnot(identical(readLines(file.path("source-archive", "SOURCE_COMMIT")), source_commit),
          identical(digest::digest(file = archive, algo = "sha256"), archive_sha))
dir.create("windows-diagnostic")
report <- normalizePath("windows-diagnostic", winslash = "/", mustWork = TRUE)
runner <- normalizePath("validation/run_full_tests.R", winslash = "/", mustWork = TRUE)
files <- c("R/chained_augmented.R", "R/chained_conditional_engine.R",
  "man/conditional_did_pretest.Rd", "src/cdid_fast_impaug.c",
  "tests/testthat/test-chained-conditional-pretest.R",
  "tests/testthat/test-cp02-native-overhead.R")
# Keep the archive's frozen reference bytes. Read the six patched files
# from Git blobs rather than from a Windows checkout with CRLF conversion.
staging <- tempfile("cdid-source-")
dir.create(staging)
untar(archive, exdir = staging)
source <- file.path(staging, "cdid")
for (file in files) {
  stopifnot(file.exists(file.path(source, file)))
  status <- system2("git", c("show", shQuote(paste0("HEAD:", file))),
    stdout = file.path(source, file))
  stopifnot(identical(status, 0L))
}
writeLines(c(paste("Original archive SHA256", archive_sha),
  paste("Diagnostic revision", Sys.getenv("GITHUB_SHA")),
  paste(unname(tools::md5sum(file.path(source, files))), files)),
  file.path(report, "source-files.txt"))
lib <- tempfile("cdid-library-")
dir.create(lib)
log <- file.path(report, "install.log")
status <- system2(file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", paste0("--library=", shQuote(lib)), shQuote(source)),
  stdout = log, stderr = log)
stopifnot(identical(status, 0L))
native <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(runner), "--native-only",
    shQuote(paste0("--lib-loc=", lib)),
    shQuote(paste0("--test-dir=", file.path(source, "tests/testthat"))),
    shQuote(paste0("--output-dir=", file.path(report, "native-overhead")))),
  stdout = file.path(report, "native-console.log"), stderr = file.path(report, "native-console.log")))
writeLines(as.character(native), file.path(report, "native-status.txt"))
.libPaths(c(lib, .libPaths()))
library(cdid)
stopifnot(identical(normalizePath(dirname(find.package("cdid"))), normalizePath(lib)))
writeLines(capture.output(sessionInfo()), file.path(report, "environment.txt"))
reporter <- testthat::ListReporter$new()
tests <- testthat::test_dir(file.path(source, "tests/testthat"),
  filter = "^(augmented-replay-operation|augmented-outcome-pooling|chained-conditional-pretest|cran-regression(-static)?)$",
  package = "cdid", load_package = "installed",
  reporter = testthat::MultiReporter$new(reporters = list(
    reporter, testthat::SummaryReporter$new(file = file.path(report, "tests.log")))),
  stop_on_failure = FALSE, stop_on_warning = FALSE
)
rows <- as.data.frame(tests)
cols <- c("file", "test", "failed", "error", "warning", "skipped", "passed")
write.csv(rows[cols], file.path(report, "tests.csv"), row.names = FALSE)
stopifnot(length(unique(rows$file)) == 5L)
print(aggregate(cbind(failed, error, warning, skipped, passed) ~ file, rows, sum))
if (native != 0L || any(rows$failed > 0L | rows$error | rows$warning > 0L)) quit(status = 1L)
