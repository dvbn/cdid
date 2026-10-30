#!/usr/bin/env Rscript
# Dependency setup only; run with Rscript --vanilla. No estimation is run.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(project = NULL, library = NULL, repository = "https://cloud.r-project.org")
  mode <- "check-only"
  flags <- character()
  while (length(args)) {
    key <- args[1L]; args <- args[-1L]
    if (key %in% c("--install", "--check-only", "--plan")) {
      flags <- c(flags, key); mode <- sub("^--", "", key)
    } else {
      name <- sub("^--", "", key)
      if (!startsWith(key, "--") || !name %in% names(opts) || !length(args))
        stop("Use --project PATH [--library PATH] [--repository HTTPS_URL] ",
             "and one of --install / --check-only / --plan.")
      opts[[name]] <- args[1L]; args <- args[-1L]
    }
  }
  if (length(flags) > 1L || is.null(opts$project)) stop("Supply --project and at most one mode.")
  canon <- function(p) {
    p <- path.expand(p)
    if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", p)) p <- file.path(getwd(), p)
    if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
    parent <- dirname(p)
    if (identical(parent, p)) stop("Cannot resolve path: ", p)
    file.path(canon(parent), basename(p))
  }
  project <- canon(opts$project)
  if (!dir.exists(project) ||
      !all(file.exists(file.path(project, c("DESCRIPTION", "validation/setup.R")))))
    stop("--project must be the cdid checkout root containing validation/setup.R.")
  identity <- read.dcf(file.path(project, "DESCRIPTION"), fields = c("Package", "Version"))[1L, ]
  if (!identical(unname(identity[["Package"]]), "cdid") ||
      !identical(unname(identity[["Version"]]), "1.1.0.0"))
    stop("--project must select a cdid 1.1.0.0 source checkout.")
  receipt_dir <- file.path(project, "validation", "setup-receipts")
  if (!dir.exists(receipt_dir) && !dir.create(receipt_dir))
    stop("Cannot create project-local receipt directory.")
  stem <- paste0(format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC"), "-", Sys.getpid())
  attempt <- file.path(receipt_dir, gsub("[^0-9A-Za-z.-]", "", stem))
  suffix <- 0L
  while (file.exists(attempt)) {
    suffix <- suffix + 1L; attempt <- file.path(receipt_dir, paste0(stem, "-", suffix))
  }
  if (!dir.create(attempt)) stop("Cannot create a fresh receipt directory.")
  receipt <- list(status = "failed", mode = mode, started_utc = stem, warnings = character(),
                  minimum_R = "4.1.0", tested_R_baseline = "4.5.2", project = project)
  on.exit({
    receipt$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    receipt$session <- capture.output(sessionInfo())
    dput(receipt, file = file.path(attempt, "receipt.dput"))
    message("Receipt: ", attempt)
  }, add = TRUE)
  tryCatch(withCallingHandlers({
    if (getRversion() < "4.1.0") stop("R >= 4.1.0 is required; 4.5.2 is the tested baseline.")
    lib <- canon(if (is.null(opts$library)) file.path(project, "validation", ".lib") else opts$library)
    inside <- function(a, b) {
      if (.Platform$OS.type == "windows") { a <- tolower(a); b <- tolower(b) }
      identical(a, b) || startsWith(paste0(a, "/"), paste0(b, "/"))
    }
    user_libs <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1L]]
    protected <- unique(vapply(c(.Library, .Library.site, user_libs[nzchar(user_libs)],
                                 file.path(project, ".runtime")), canon, ""))
    if (any(vapply(protected, function(p) inside(lib, p) || inside(p, lib), logical(1))) ||
        inside(project, lib)) stop("Choose a dedicated library outside base/user/site libraries, project .runtime, and source root.")
    if (file.exists(lib) && !dir.exists(lib)) stop("Library path is not a directory.")
    source_dirs <- file.path(project, c("R", "src", "inst", "man", "tests", "data", "exec", "po", "tools", ".git"))
    if (any(vapply(source_dirs, function(p) inside(lib, canon(p)), logical(1))))
      stop("Library must be separate from package source and Git metadata directories.")
    repo <- sub("/+$", "", opts$repository)
    if (!grepl("^https://[^/@?#[:space:]]+(/[^?#[:space:]]*)?$", repo))
      stop("Repository must be an HTTPS CRAN repository URL without credentials/query/fragment.")
    receipt$library <- lib; receipt$repository <- repo
    .libPaths(c(if (dir.exists(lib)) lib, .libPaths()))
    fields <- c("Package", "Version", "Depends", "Imports", "LinkingTo")
    description <- function(path) as.list(read.dcf(path, fields = fields)[1L, ])
    cdid <- description(file.path(project, "DESCRIPTION"))
    pins <- c(DRDID = "1.3.0", did = "2.5.1", cdid = "1.1.0.0")
    if (!identical(cdid$Package, "cdid") || !identical(cdid$Version, pins[["cdid"]]))
      stop("Expected selected cdid source version 1.1.0.0.")
    if (!all(file.exists(file.path(project, c("NAMESPACE", "R", "src")))))
      stop("A complete cdid source checkout, including compiled src, is required.")
    native_files <- list.files(file.path(project, "src"), recursive = TRUE, all.files = TRUE)
    if (any(grepl("\\.(o|obj|so|dll|dylib|a)$", native_files, ignore.case = TRUE)))
      stop("Source checkout contains native build outputs; select a clean source checkout.")
    vendor <- file.path(project, "validation", "vendor")
    archives <- read.csv(file.path(vendor, "SOURCES.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("package", "version", "license", "file", "source_url", "archive_url", "bytes", "md5",
                  "sha256", "depends", "imports", "linking_to")
    if (!identical(names(archives), required) || nrow(archives) != 2L ||
        !identical(archives$package, names(pins)[1:2]) ||
        !identical(archives$version, unname(pins[1:2])) || any(archives$license != "GPL-3") ||
        any(archives$file != paste0(archives$package, "_", archives$version, ".tar.gz")) ||
        any(archives$source_url != paste0("https://cran.r-project.org/src/contrib/", archives$file)) ||
        any(archives$archive_url != paste0("https://cran.r-project.org/src/contrib/Archive/",
                                         archives$package, "/", archives$file)) ||
        any(!grepl("^[0-9a-f]{32}$", archives$md5)) || any(!grepl("^[0-9a-f]{64}$", archives$sha256)))
      stop("Invalid pinned source manifest.")
    archive_paths <- file.path(vendor, archives$file)
    if (any(!file.exists(archive_paths))) stop("Missing bundled source archives; restore the original files listed in validation/vendor/SOURCES.csv.")
    actual_md5 <- unname(tools::md5sum(archive_paths))
    receipt$archives <- transform(archives, actual_md5 = actual_md5)
    if (anyNA(actual_md5) || any(actual_md5 != archives$md5) ||
        any(file.info(archive_paths)$size != archives$bytes)) stop("Bundled archive MD5/size verification failed.")
    sources <- setNames(lapply(seq_len(nrow(archives)), function(i) {
      setNames(as.list(archives[i, c("package", "version", "depends", "imports", "linking_to")]), fields)
    }), archives$package)
    sources$cdid <- cdid
    # Only package source trees enter the fingerprint, never libraries, receipts or .git.
    manifest <- function(root, entries, output) {
      files <- unique(unlist(lapply(file.path(root, entries), function(p) {
        if (dir.exists(p)) list.files(p, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
        else if (file.exists(p)) p else character()
      })))
      files <- sort(files[!file.info(files)$isdir])
      data <- data.frame(file = substring(files, nchar(root) + 2L), md5 = unname(tools::md5sum(files)))
      if (anyNA(data$md5)) stop("Cannot hash all source/installed files.")
      write.csv(data, output, row.names = FALSE)
      unname(tools::md5sum(output))
    }
    entries <- c("DESCRIPTION", "NAMESPACE", "R", "src", "man", "data", "inst", "tests", "exec", "po", "tools",
                 "configure", "configure.win", "cleanup", "cleanup.win", "LICENSE", "LICENSE.note",
                 "LICENCE", "LICENCE.note", "NEWS", "NEWS.md", "README", "README.md", "INDEX", ".Rbuildignore")
    git <- Sys.which("git")
    if (nzchar(git)) {
      git_run <- function(args) {
        result <- system2(git, c("-C", shQuote(project), args), stdout = TRUE, stderr = TRUE)
        if (!is.null(attr(result, "status"))) stop("Git inspection failed; use a clean Git clone of the cdid repository.")
        result
      }
      if (!identical(canon(git_run(c("rev-parse", "--show-toplevel"))), project))
        stop("--project must select the root of the cdid Git clone.")
      if (length(git_run(c("status", "--porcelain", "--untracked-files=no"))) ||
          length(git_run(c("ls-files", "--others", "--", shQuote(entries)))))
        stop("Checkout has changed tracked files or untracked package sources; use a clean checkout.")
      receipt$git_revision <- git_run(c("rev-parse", "HEAD"))
    } else if (file.exists(file.path(project, ".git"))) stop("Git is required to inspect this checkout.")
    else receipt$git_revision <- "unavailable; source file hashes recorded"
    receipt$source_md5 <- manifest(project, entries, file.path(attempt, "source-files.csv"))
    # Parse ordinary dependency fields, preserving all declared version constraints.
    parse_deps <- function(d) {
      values <- unlist(d[c("Depends", "Imports", "LinkingTo")], use.names = FALSE)
      terms <- trimws(unlist(strsplit(paste(values[!is.na(values)], collapse = ","), ",", fixed = TRUE)))
      terms <- terms[nzchar(terms)]
      result <- data.frame(package = character(), op = character(), version = character())
      for (term in terms) {
        match <- regmatches(term, regexec("^([A-Za-z][A-Za-z0-9.]*)[[:space:]]*(\\(([><=!]+)[[:space:]]*([0-9][0-9A-Za-z.-]*)\\))?$", term))[[1L]]
        if (!length(match) || (!match[4L] %in% c("", ">=", "<=", "==", ">", "<", "!=")))
          stop("Unsupported dependency declaration: ", term)
        result <- rbind(result, data.frame(package = match[2L], op = match[4L], version = match[5L]))
      }
      unique(result)
    }
    satisfies <- function(version, req) {
      all(vapply(seq_len(nrow(req)), function(i) {
        op <- req$op[i]; if (!nzchar(op)) return(TRUE)
        cmp <- utils::compareVersion(version, req$version[i])
        switch(op, ">=" = cmp >= 0L, "<=" = cmp <= 0L, "==" = cmp == 0L,
               ">" = cmp > 0L, "<" = cmp < 0L, "!=" = cmp != 0L)
      }, logical(1)))
    }
    suggested <- as.list(read.dcf(file.path(project, "DESCRIPTION"), fields = "Suggests")[1L, ])
    validation <- parse_deps(list(Imports = suggested[[1L]]))
    validation <- validation[validation$package %in% c("testthat", "digest", "ggplot2"), ]
    if (!setequal(validation$package, c("testthat", "digest", "ggplot2")))
      stop("Expected existing validation dependencies in cdid Suggests.")
    roots <- unique(do.call(rbind, c(lapply(sources, parse_deps), list(validation))))
    installed <- installed.packages(lib.loc = .libPaths(), fields = fields[-(1:2)])
    installed <- installed[!duplicated(installed[, "Package"]), , drop = FALSE]
    receipt$initial_versions <- installed[, c("Package", "Version", "LibPath"), drop = FALSE]
    available <- NULL; chosen <- list(); unresolved <- character()
    # Recompute reachable requirements when a selected dependency version changes.
    for (pass in seq_len(100L)) {
      req <- roots; todo <- unique(req$package); seen <- character()
      while (length(todo)) {
        p <- todo[1L]; todo <- todo[-1L]
        if (p %in% seen) next
        seen <- c(seen, p); d <- if (p %in% names(sources)) sources[[p]] else chosen[[p]]
        if (!is.null(d)) {
          extra <- parse_deps(d); req <- unique(rbind(req, extra))
          todo <- unique(c(todo, setdiff(extra$package, seen)))
        }
      }
      next_chosen <- list(); unresolved <- character()
      for (p in unique(req$package)) {
        bounds <- req[req$package == p, ]
        if (p == "R") {
          if (!satisfies(as.character(getRversion()), bounds)) stop("R cannot satisfy declared dependency bounds.")
          next
        }
        if (p %in% names(pins)) {
          if (!satisfies(pins[[p]], bounds)) stop("Declared dependencies conflict with exact pin: ", p)
          next
        }
        d <- if (p %in% rownames(installed)) as.list(installed[p, fields]) else NULL
        if (!is.null(d) && satisfies(d$Version, bounds)) {
          d$origin <- "installed"
        } else if (mode != "install") {
          unresolved <- c(unresolved, p); d <- NULL
        } else {
          if (is.null(available)) available <- available.packages(repos = repo, type = "source")
          if (!p %in% rownames(available)) stop("No source candidate for ", p, " on this R/repository.")
          d <- as.list(available[p, fields]); d$origin <- "repository"
          if (!satisfies(d$Version, bounds)) stop("Repository source cannot satisfy bounds for ", p)
        }
        if (!is.null(d)) next_chosen[[p]] <- d
      }
      stable <- identical(next_chosen, chosen)
      chosen <- next_chosen
      if (stable) break
      if (pass == 100L) stop("Dependency requirements did not converge; review repository compatibility.")
    }
    write.csv(req, file.path(attempt, "requirements.csv"), row.names = FALSE)
    ordinary <- names(chosen)[vapply(chosen, function(d) d$origin == "repository", logical(1))]
    pin_ok <- vapply(names(pins), function(p) {
      path <- file.path(lib, p, "DESCRIPTION")
      file.exists(path) && identical(description(path)$Version, pins[[p]]) &&
        identical(description(path)$Package, p) && identical(dirname(canon(dirname(path))), lib)
    }, logical(1))
    receipt$required_pins <- pins
    receipt$planned <- c(if (length(c(unresolved, ordinary))) paste("CRAN source:", unique(c(unresolved, ordinary))),
                         if (any(!pin_ok[1:2])) paste("Bundled source:", names(pins)[1:2][!pin_ok[1:2]]),
                         "cdid: install current checkout source when --install is selected")
    message(paste(receipt$planned, collapse = "\n"))
    if (length(unresolved)) message("Unavailable locally: ", paste(unresolved, collapse = ", "),
                                   ". Transitive requirements of missing sources require repository metadata during --install.")
    rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
    need_tools <- mode == "install" || length(unresolved) || any(!pin_ok)
    if (need_tools) {
      for (key in c("CC", "CXX")) {
        config <- system2(rbin, c("CMD", "config", key), stdout = TRUE, stderr = TRUE)
        if (!is.null(attr(config, "status")) || !length(config)) stop("Cannot inspect R compiler: ", key)
        compiler <- scan(text = config[1L], what = "", quiet = TRUE)[1L]
        if (is.na(compiler) || !nzchar(Sys.which(compiler))) stop("Required compiler unavailable: ", key)
      }
      if (!nzchar(Sys.which("make"))) stop("Source installation requires make (Rtools on Windows).")
    }
    if (mode != "install" && (length(unresolved) || any(!pin_ok)))
      stop("Preflight incomplete: dependencies/pins need installation. No network or installation was attempted.")
    if (mode == "install") {
      if (!dir.exists(lib) && !dir.create(lib, recursive = TRUE)) stop("Cannot create explicit library.")
      .libPaths(c(lib, .libPaths()))
      receipt$phase <- "ordinary dependencies"
      check_destination <- function(p) {
        dest <- file.path(lib, p)
        link <- Sys.readlink(dest)
        if ((!is.na(link) && nzchar(link)) || (file.exists(dest) && !identical(dirname(canon(dest)), lib)))
          stop("Refusing to replace a package linked outside the explicit library: ", p)
      }
      for (p in unique(c(ordinary, names(pins)))) check_destination(p)
      done <- character()
      install_one <- function(p, active = character()) {
        if (p %in% done || !p %in% ordinary) return(invisible(NULL))
        if (p %in% active) stop("Dependency install cycle involving ", p)
        for (dep in parse_deps(chosen[[p]])$package) install_one(dep, c(active, p))
        install.packages(p, lib = lib, repos = repo, type = "source", dependencies = FALSE, available = available)
        d <- description(file.path(lib, p, "DESCRIPTION"))
        if (!identical(d$Version, chosen[[p]]$Version)) stop("Source install failed/version changed: ", p)
        done <<- c(done, p)
      }
      for (p in ordinary) install_one(p)
    }
    # digest is an existing validation dependency, not a bootstrap requirement for MD5.
    if (!requireNamespace("digest", quietly = TRUE)) stop("digest unavailable for SHA-256 verification.")
    actual_sha256 <- vapply(archive_paths, function(p) digest::digest(file = p, algo = "sha256"), "")
    receipt$archives$actual_sha256 <- unname(actual_sha256)
    if (any(actual_sha256 != archives$sha256)) stop("Bundled archive SHA-256 verification failed.")
    if (mode == "install") {
      source_install <- function(path, package) {
        receipt$phase <<- paste("source install", package)
        log <- file.path(attempt, paste0(package, "-install.log"))
        status <- system2(rbin, c("CMD", "INSTALL", shQuote(paste0("--library=", lib)), shQuote(path)),
                          stdout = log, stderr = log)
        if (status != 0L) stop("Source installation failed: ", package, "; see ", log)
      }
      for (i in seq_len(2L)) if (!pin_ok[i]) {
        # Recheck immediately before consuming a bundled source.
        if (unname(tools::md5sum(archive_paths[i])) != archives$md5[i]) stop("Archive MD5 changed before install.")
        source_install(archive_paths[i], archives$package[i])
      }
      # Compile a fresh source copy so the selected checkout stays clean for repeat checks.
      stage <- file.path(attempt, "cdid-source")
      if (!dir.create(stage)) stop("Cannot create fresh source staging directory.")
      for (entry in entries[file.exists(file.path(project, entries))])
        if (!file.copy(file.path(project, entry), stage, recursive = TRUE, copy.date = TRUE))
          stop("Could not stage cdid source: ", entry)
      staged_md5 <- manifest(stage, entries, file.path(attempt, "staged-source-files.csv"))
      if (!identical(staged_md5, receipt$source_md5)) stop("Checkout source changed during setup/staging.")
      source_install(stage, "cdid")
    }
    current <- installed.packages(lib.loc = .libPaths())
    current <- current[!duplicated(current[, "Package"]), , drop = FALSE]
    for (p in setdiff(unique(req$package), "R")) {
      if (!p %in% rownames(current) || !satisfies(current[p, "Version"], req[req$package == p, ]))
        stop("Installed dependency does not satisfy declared bounds: ", p)
    }
    # Check pinned namespace versions, installation identity, and cdid's minimum-version guard.
    versions <- do.call(rbind, lapply(names(pins), function(p) {
      ns <- loadNamespace(p, lib.loc = c(lib, .libPaths()))
      actual <- as.character(getNamespaceVersion(ns))
      path <- canon(getNamespaceInfo(ns, "path"))
      if (!identical(actual, pins[[p]]) || !identical(path, canon(file.path(lib, p))) ||
          !identical(dirname(path), lib))
        stop("Namespace version/library identity mismatch: ", p)
      data.frame(package = p, version = actual, library = lib)
    }))
    guard <- get("check_upstream_baseline", envir = asNamespace("cdid"), inherits = FALSE)
    guard()
    receipt$versions <- versions
    receipt$installed_versions <- installed.packages(lib.loc = .libPaths())[, c("Package", "Version", "LibPath"), drop = FALSE]
    receipt$source_installed_this_attempt <- mode == "install"
    receipt$status <- if (mode == "install") "setup_verified" else "preflight_verified_no_install"
    receipt$phase <- "complete"
    message(if (mode == "install") "Setup verified." else
      "Installed versions, namespaces and guard verified. Preflight only: checkout source was not installed.")
  }, warning = function(w) {
    receipt$warnings <<- c(receipt$warnings, conditionMessage(w))
  }), error = function(e) {
    receipt$error <<- conditionMessage(e)
    stop(conditionMessage(e), call. = FALSE)
  })
}
tryCatch(main(), error = function(e) { message("SETUP FAILED: ", conditionMessage(e)); quit(status = 1L) })
