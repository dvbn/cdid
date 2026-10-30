# Build the existing cdid documentation as a private local pkgdown preview.
# Rscript --vanilla build_site.R --source DIR --library DIR --logo FILE.png --output NEW_DIR
# Requires existing tooling; installs and deploys nothing. Source files are read only.

fail <- function(...) stop(paste0(...), call. = FALSE)

parse_args <- function(args) {
  keys <- c("source", "library", "logo", "output")
  opt <- list()
  i <- 1L
  while (i <= length(args)) {
    if (grepl("^--[a-z]+=", args[[i]])) {
      key <- sub("^--([a-z]+)=.*$", "\\1", args[[i]])
      value <- sub("^--[a-z]+=", "", args[[i]])
      i <- i + 1L
    } else if (grepl("^--[a-z]+$", args[[i]]) && i < length(args) &&
               !startsWith(args[[i + 1L]], "--")) {
      key <- substring(args[[i]], 3L)
      value <- args[[i + 1L]]
      i <- i + 2L
    } else {
      fail("Unexpected or incomplete argument: ", args[[i]])
    }
    if (!key %in% keys) fail("Unknown option --", key)
    if (!is.null(opt[[key]])) fail("--", key, " is given twice")
    if (!nzchar(value)) fail("--", key, " is empty")
    opt[[key]] <- value
  }
  missing <- setdiff(keys, names(opt))
  if (length(missing)) fail("Missing option(s): ", paste0("--", missing, collapse = ", "))
  opt
}

this_dir <- function() {
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(file) != 1L) fail("Run this file with Rscript.")
  dirname(normalizePath(gsub("~+~", " ", file, fixed = TRUE), mustWork = TRUE))
}

existing <- function(path, what, directory = FALSE) {
  if (!file.exists(path)) fail(what, " not found: ", path)
  if (directory && !dir.exists(path)) fail(what, " must be a directory: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

inside <- function(path, dir) startsWith(paste0(path, "/"), paste0(dir, "/"))
is_link <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}

python_call <- function(python, helper, args) {
  status <- system2(python, shQuote(c("-I", helper, args)))
  if (!identical(status, 0L)) fail("Website preparation or validation failed; see the error above.")
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  site_dir <- this_dir()
  src <- existing(opt$source, "--source", TRUE)
  lib <- existing(opt$library, "--library", TRUE)
  for (p in c("did", "cdid")) existing(file.path(lib, p, "DESCRIPTION"), paste(p, "in --library"))
  logo <- existing(opt$logo, "--logo")
  if (dir.exists(logo)) fail("--logo must be a PNG file")
  for (p in c("DESCRIPTION", "NAMESPACE", "README.md", "NEWS.md", "docs/README.md"))
    existing(file.path(src, p), paste("Source file", p))

  # Resolve the existing parent before checking containment; refuse dangling links too.
  out <- file.path(existing(dirname(opt$output), "Parent of --output", TRUE), basename(opt$output))
  if (file.exists(out) || is_link(out)) fail("--output already exists; choose a new directory: ", out)
  if (basename(out) %in% c(".", "..", "")) fail("--output must name a new directory")
  if (inside(out, src) || inside(out, lib)) fail("--output must be outside --source and --library")
  stage_root <- file.path(dirname(out), ".site-build")
  if (inside(stage_root, src) || inside(stage_root, lib) || inside(out, stage_root))
    fail("Choose --output with its parent outside --source, --library and .site-build")
  if (is_link(stage_root) || (file.exists(stage_root) && !dir.exists(stage_root)))
    fail(".site-build must be an ordinary directory")

  python <- Sys.which("python3")
  if (!nzchar(python)) fail("python3 is not on PATH")
  .libPaths(c(lib, .libPaths()))
  for (p in c("pkgdown", "rmarkdown", "knitr", "bslib", "ragg", "jsonlite"))
    if (!requireNamespace(p, quietly = TRUE)) fail("Package not installed: ", p)
  if (packageVersion("pkgdown") != "2.2.0") fail("This builder requires pkgdown 2.2.0")
  if (!rmarkdown::pandoc_available("2.0")) fail("pandoc is not available")
  message("pkgdown ", packageVersion("pkgdown"), ", rmarkdown ", packageVersion("rmarkdown"),
          ", knitr ", packageVersion("knitr"), ", bslib ", packageVersion("bslib"),
          ", ragg ", packageVersion("ragg"), ", pandoc ", rmarkdown::pandoc_version())

  if (!dir.exists(stage_root) && !dir.create(stage_root)) fail("Cannot create .site-build")
  stage <- file.path(stage_root, paste0("stage-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid()))
  built <- FALSE
  on.exit(if (!built) {
    if (dir.exists(stage)) message("Stage kept for inspection: ", stage)
    if (file.exists(out)) message("Partial site left in ", out, "; choose a new output for the next build.")
  })
  python_call(python, file.path(site_dir, "stage_site.py"),
              c("--source", src, "--stage", stage, "--config", file.path(site_dir, "_pkgdown.yml"), "--logo", logo))
  pkg_dir <- file.path(stage, "pkg")
  manifest_file <- file.path(stage, "site-manifest.json")
  manifest <- jsonlite::read_json(manifest_file, simplifyVector = TRUE)

  # A pkgdown object carries the explicit destination through every build step.
  # This also checks nested article names and the explicitly indexed reexports topic.
  pkg <- pkgdown::as_pkgdown(pkg_dir, override = list(destination = out))
  if (!identical(as.character(pkg$dst_path), out)) fail("pkgdown destination differs from --output")
  article_names <- vapply(manifest$pages$stage[startsWith(manifest$pages$stage, "pkg/vignettes/")],
                          function(x) sub("\\.Rmd$", "", sub("^pkg/vignettes/", "", x)), "")
  if (!setequal(pkg$vignettes$name, article_names)) fail("pkgdown did not discover every staged article")
  if (!all(manifest$reference_topics %in% pkg$topics$name)) fail("pkgdown did not discover every help topic")

  # Suppress pkgdown's automatic remote favicon service. No favicon generation is needed.
  old_ci <- Sys.getenv("CI", unset = NA_character_)
  Sys.setenv(CI = "true")
  on.exit(if (is.na(old_ci)) Sys.unsetenv("CI") else Sys.setenv(CI = old_ci), add = TRUE)
  pkgdown::build_site(pkg, examples = FALSE, run_dont_run = FALSE, seed = NULL,
                      preview = FALSE, devel = FALSE, new_process = FALSE, install = FALSE)

  # Copy recorded assets explicitly; do not rely on rmarkdown resource discovery.
  for (i in seq_len(nrow(manifest$figures))) {
    dst <- file.path(out, manifest$figures$site[i])
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(stage, manifest$figures$stage[i]), dst, overwrite = TRUE))
      fail("Could not copy figure: ", manifest$figures$site[i])
  }
  for (name in manifest$logos$site)
    if (!file.copy(file.path(pkg_dir, name), file.path(out, name), overwrite = TRUE)) fail("Could not copy logo: ", name)
  python_call(python, file.path(site_dir, "check_site.py"),
              c("--output", out, "--manifest", manifest_file,
                "--host-path", src, "--host-path", lib, "--host-path", stage,
                "--host-path", normalizePath("~", winslash = "/", mustWork = FALSE)))
  built <- TRUE
  unlink(stage, recursive = TRUE)
  message("Static website checks passed. Browser review remains: ", file.path(out, "index.html"))
}

main()
