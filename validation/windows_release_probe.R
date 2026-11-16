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
link_log <- file.path(report, "original-build.log")
original <- gsub("stdout = TRUE, stderr = TRUE", "stdout = link_log, stderr = link_log", original, fixed = TRUE)
original <- gsub('!is.null(attr(out, "status"))', '!identical(out, 0L)', original, fixed = TRUE)
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

# Run the DLL unload/reload tests in a process without live fit objects.
native_args <- c("--vanilla", "validation/run_full_tests.R", "--native-only",
  shQuote(paste0("--lib-loc=", normalizePath(lib, winslash = "/"))),
  "--test-dir=tests/testthat",
  shQuote(paste0("--output-dir=", file.path(report, "native-overhead"))))
native_log <- file.path(report, "native-console.log")
native_status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), native_args,
  stdout = native_log, stderr = native_log))
writeLines(as.character(native_status), file.path(report, "native-status.txt"))

insert_probe("tests/testthat/test-augmented-replay-operation.R",
  "fast <- do.call(cdid::att_gt, request)", c(
    '    contexts <- list()',
    '    native <- v45_ns("cdid_ia_fast_aug_if_chunk_core")',
    '    testthat::local_mocked_bindings(cdid_ia_fast_aug_if_chunk_core =',
    '      function(st, X_F, backend, rows, links, scale, pre, pair_op = NULL) {',
    '        value <- native(st, X_F, backend, rows, links, scale, pre, pair_op)',
    '        if (!is.null(value) && spec$selection_model == "pair_logit" && !length(contexts))',
    '          contexts[[1L]] <<- list(st = st, X_F = X_F, backend = backend)',
    '        value',
    '      }, .package = "cdid")'
  )
)
insert_probe("tests/testthat/test-augmented-replay-operation.R", "old <- default_fit()", c(
  '    if (length(contexts)) {',
  '      st <- contexts[[1L]]$st',
  '      backend <- contexts[[1L]]$backend',
  '      X_F <- contexts[[1L]]$X_F',
  '      trace <- list()',
  '      for (g in c(3, 4)) {',
  '        links <- which(vapply(backend$links, function(link) {',
  '          st$cohorts[link$target_cohort_index] == g && link$pair$s == 1 && link$pair$t == 2',
  '        }, logical(1)))',
  '        if (!length(links)) next',
  '        for (ell in links) {',
  '          link <- backend$links[[ell]]',
  '          bid <- link$or_block_id',
  '          block <- backend$shared_or$blocks[[bid]]',
  '          active <- st$cohort_index %in% c(link$target_cohort_index, link$control_cohort_indices)',
  '          rr <- c(which(!active)[1L], which(active)[1L])',
  '          pre <- v45_ns("augmented_if_needed_sets")(backend, ell, FALSE)',
  '          a <- native(st, X_F, backend, rr, ell, "frame", pre)',
  '          b <- v45_ns("evaluate_augmented_if_chunk")(backend, rr, ell, details = TRUE)',
  '          if (is.null(a)) next',
  '          fit <- v45_ns("reconstruction_or_block")(backend$shared_or, bid)$fit',
  '          trace[[length(trace) + 1L]] <- list(cohort = g, id = st$id[rr],',
  '            native = sprintf("%a", drop(a$phi)), ordinary = sprintf("%a", drop(b$phi)),',
  '            components = b$details$links[[1L]]$components, mean_alpha = block$score_mean_alpha,',
  '            mean_beta = block$score_mean_beta, mass = block$mass, WLS = fit, mean_d = block$mean_d,',
  '            d_plus = link$d_bar_plus_or, d_minus = link$d_bar_minus_or)',
  '        }',
  '      }',
  '      writeLines(capture.output(dput(trace)), file.path(Sys.getenv("CDID_WINDOWS_PROBE_REPORT"),',
  '        "augmented-primitive.txt"))',
  '    }'
))

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
reporter <- testthat::ListReporter$new()
tests <- testthat::test_dir("tests/testthat",
  filter = "^(augmented-replay-operation|chained-conditional-pretest)$",
  package = "cdid", load_package = "installed",
  reporter = testthat::MultiReporter$new(reporters = list(
    reporter, testthat::SummaryReporter$new(file = file.path(report, "tests.log")))),
  stop_on_failure = FALSE, stop_on_warning = FALSE
)
rows <- as.data.frame(tests)
cols <- c("file", "test", "failed", "error", "warning", "skipped", "passed")
write.csv(rows[cols], file.path(report, "tests.csv"), row.names = FALSE)
stopifnot(length(unique(rows$file)) == 2L)
print(aggregate(cbind(failed, error, warning, skipped, passed) ~ file, rows, sum))
if (native_status != 0L || any(rows$failed > 0L | rows$error | rows$warning > 0L)) quit(status = 1L)
