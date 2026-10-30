# The oldest upstream versions cdid accepts, and the versions it has been
# tested against. Loading refuses a version below the minimum. Later versions
# are accepted, but they are untested: a larger version number is not evidence
# of compatibility. A version is added to `tested_upstream` only after the
# ordinary-mode parity tests pass against it.
minimum_upstream <- c(did = "2.5.1", DRDID = "1.3.0")
tested_upstream <- c(did = "2.5.1", DRDID = "1.3.0")

namespace_path <- function(pkg) {
  normalizePath(getNamespaceInfo(pkg, "path"), winslash = "/",
                mustWork = FALSE)
}

# TRUE where `version` is at least `minimum`, compared as package versions
# (so "2.10.0" is later than "2.5.1"); NA where either is missing.
upstream_at_least <- function(version, minimum) {
  minimum <- rep_len(as.character(minimum), length(version))
  ok <- rep(NA, length(version))
  known <- !is.na(version) & !is.na(minimum)
  if (any(known)) {
    ok[known] <- package_version(version[known]) >= package_version(minimum[known])
  }
  ok
}

# One row per package: the loaded version, the tested version, whether loading
# is permitted (`ok`, the version is at least `minimum`), and whether the
# loaded version is the tested one (`matches_tested`; NA for cdid itself).
upstream_version_table <- function(package, version, path) {
  minimum <- unname(minimum_upstream[package])
  tested <- unname(tested_upstream[package])
  data.frame(
    package = package,
    version = version,
    tested = tested,
    ok = is.na(minimum) | upstream_at_least(version, minimum) %in% TRUE,
    path = path,
    minimum = minimum,
    matches_tested = version == tested,
    stringsAsFactors = FALSE
  )
}

# Private diagnostic used by the namespace baseline check.
cdid_dependencies <- function() {
  # Explicitly use the declared dependency; versions below still describe the
  # loaded namespaces, including any namespace already loaded by did.
  loadNamespace("DRDID")
  pkgs <- c("cdid", names(minimum_upstream))
  version <- vapply(pkgs, function(p) {
    as.character(getNamespaceVersion(p))
  }, character(1))
  upstream_version_table(
    pkgs, unname(version),
    vapply(pkgs, namespace_path, character(1), USE.NAMES = FALSE)
  )
}

check_upstream_baseline <- function(deps = cdid_dependencies()) {
  deps <- deps[!is.na(deps$minimum), , drop = FALSE]
  bad <- deps[!deps$ok, , drop = FALSE]
  if (nrow(bad) == 0L) {
    return(invisible(TRUE))
  }
  found <- paste0(bad$package, " ", bad$version, " (loaded from ", bad$path,
                  ")", collapse = "; ")
  wanted <- paste0(bad$package, " ", bad$minimum, " or later",
                   collapse = " and ")
  tested <- paste0(names(tested_upstream), " ", tested_upstream,
                   collapse = " and ")
  cdid_abort(
    paste0("cdid ", cdid_version(), " requires ", wanted,
           ", but this session loaded ", found, ". Update ",
           paste(bad$package, collapse = " and "), ", for example with ",
           "install.packages(), then start a fresh R session before loading ",
           "cdid. cdid has been tested with ", tested, "."),
    "cdid_untested_upstream"
  )
}

.onLoad <- function(libname, pkgname) {
  check_upstream_baseline()
}
