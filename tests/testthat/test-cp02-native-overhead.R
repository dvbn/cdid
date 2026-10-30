# Focused tests for CP02 native DLL lookup caching, owner/symbol/lifecycle/FP
# guards, redundant owner-validation avoidance, and C-level interrupt polling
# (`cdid/src/cdid.c`, `cdid/R/chained_gmm_core.R`, `cdid/R/chained_ipw_owner.R`).

noh_raw_bits <- function(x) {
  writeBin(as.double(x), raw(), size = 8L, endian = "little")
}

noh_locate_pkg_root <- function() {
  cands <- c(
    tryCatch(testthat::test_path("..", ".."), error = function(e) character(0L)),
    tryCatch(testthat::test_path("..", "..", "00_pkg_src", "cdid"), error = function(e) character(0L)),
    file.path("..", "..", "00_pkg_src", "cdid"),
    file.path("..", ".."),
    "cdid",
    "."
  )
  for (p in cands) {
    if (nzchar(p) &&
        file.exists(file.path(p, "R", "chained_gmm_core.R")) &&
        file.exists(file.path(p, "R", "chained_ipw_owner.R"))) {
      return(normalizePath(p, mustWork = TRUE))
    }
  }
  NULL
}

noh_locate_c_source <- function() {
  root <- noh_locate_pkg_root()
  if (!is.null(root) && file.exists(file.path(root, "src", "cdid.c"))) {
    return(normalizePath(file.path(root, "src", "cdid.c"), mustWork = TRUE))
  }
  cands <- c(
    tryCatch(testthat::test_path("..", "..", "src", "cdid.c"), error = function(e) character(0L)),
    file.path("..", "..", "src", "cdid.c"),
    file.path("src", "cdid.c")
  )
  for (p in cands) {
    if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  }
  NULL
}

noh_locate_reference_original <- function() {
  root <- testthat::test_path("reference", "cp02_native_original")
  # Portable base-R byte guards; canonical SHA256 pins are recorded in packaging verification.
  expected <- c(
    "R/chained_gmm_core.R" = "553365a02f1d390d6e4d1fa0f9129736",
    "R/chained_ipw_owner.R" = "8402158457e7f6b88a1df0bb5600d553",
    "src/cdid.c" = "92ceb4294e897e084465d1d7f2631873"
  )
  paths <- file.path(root, names(expected))
  if (!all(file.exists(paths))) {
    stop("Missing packaged CP02 original reference file(s): ",
      paste(names(expected)[!file.exists(paths)], collapse = ", "), call. = FALSE)
  }
  if (!identical(unname(tools::md5sum(paths)), unname(expected))) {
    stop("Wrong packaged CP02 original reference bytes.", call. = FALSE)
  }
  normalizePath(root, mustWork = TRUE)
}

noh_require_sources <- function(pkg_root, c_src) {
  if (is.null(pkg_root) && is.null(c_src)) {
    testthat::skip("CP02 source-only test requires current R/C sources, unavailable in installed binary tests.")
  }
  if (is.null(pkg_root) || is.null(c_src)) {
    stop("Incomplete current R/C source tree for CP02 source-only test.", call. = FALSE)
  }
  invisible(NULL)
}

# Create a temporary directory on an executable filesystem (falling back to a
# workspace-local directory if `tempdir()` is mounted `noexec`).
noh_make_exec_dir <- function(prefix = "noh_build_") {
  base_candidates <- unique(c(
    tempdir(),
    if (!is.null(noh_locate_pkg_root())) dirname(noh_locate_pkg_root()) else character(0L),
    getwd()
  ))
  for (b in base_candidates) {
    if (!dir.exists(b)) next
    d <- tempfile(pattern = prefix, tmpdir = b)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(d)) next
    probe_c <- file.path(d, "probe.c")
    writeLines(c(
      "#include <R.h>",
      "#include <Rinternals.h>",
      "SEXP C_noh_probe(void) { return Rf_ScalarLogical(1); }"
    ), probe_c)
    old <- setwd(d)
    out <- suppressWarnings(system2(
      file.path(R.home("bin"), "R"),
      c("CMD", "SHLIB", "probe.c"),
      stdout = TRUE, stderr = TRUE
    ))
    setwd(old)
    so <- file.path(d, paste0("probe", .Platform$dynlib.ext))
    loaded <- if (file.exists(so)) {
      tryCatch({
        dll <- dyn.load(so)
        dyn.unload(so)
        inherits(dll, "DLLInfo")
      }, error = function(e) FALSE)
    } else {
      FALSE
    }
    unlink(c(probe_c, file.path(d, "probe.o"), so), force = TRUE)
    if (isTRUE(loaded)) {
      return(normalizePath(d, mustWork = TRUE))
    }
    unlink(d, recursive = TRUE, force = TRUE)
  }
  stop("Unable to create an executable temporary directory for R CMD SHLIB.", call. = FALSE)
}

noh_compile_dll <- function(c_lines, dll_name, parent_dir) {
  sub_dir <- tempfile(pattern = paste0("dll_", dll_name, "_"), tmpdir = parent_dir)
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  src <- file.path(sub_dir, paste0(dll_name, ".c"))
  writeLines(c_lines, src)
  old <- setwd(sub_dir)
  on.exit(setwd(old), add = TRUE)
  out <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", shQuote(basename(src))),
    stdout = TRUE, stderr = TRUE
  )
  so <- file.path(sub_dir, paste0(dll_name, .Platform$dynlib.ext))
  if (!is.null(attr(out, "status")) || !file.exists(so)) {
    stop("R CMD SHLIB failed for ", dll_name, ":\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  normalizePath(so, mustWork = TRUE)
}

noh_populate_env_from_pkg <- function(env, pkg_root = NULL, keep.source = FALSE) {
  if (!is.null(pkg_root)) {
    r_files <- sort(list.files(file.path(pkg_root, "R"), pattern = "\\.[rR]$", full.names = TRUE))
    for (p in r_files) {
      sys.source(p, envir = env, keep.source = keep.source)
    }
  } else if (isNamespaceLoaded("cdid")) {
    ns <- asNamespace("cdid")
    for (nm in ls(ns, all.names = TRUE)) {
      val <- get(nm, envir = ns, inherits = FALSE)
      if (is.function(val) && identical(environment(val), ns)) {
        environment(val) <- env
      } else if (identical(nm, "gmm_native_dll_cache") && is.environment(val)) {
        val <- new.env(parent = emptyenv())
      }
      assign(nm, val, envir = env)
    }
  } else {
    stop("Neither source tree nor loaded `cdid` namespace is available.", call. = FALSE)
  }
  env
}

noh_load_source_owner_env <- function(pkg_root = NULL, keep.source = FALSE) {
  env <- new.env(parent = baseenv())
  noh_populate_env_from_pkg(env, pkg_root = pkg_root, keep.source = keep.source)
}

# Construct a synthetic namespace-like environment (`isNamespace(ns) == TRUE`)
# with its own `.__NAMESPACE__.$DLLs` registry for deterministic testing of
# `isNamespace(owner)` branches without modifying installed packages.
noh_make_mock_namespace <- function(name, dll_info = NULL, pkg_root = NULL) {
  ns <- new.env(hash = TRUE, parent = baseenv())
  ns$.__NAMESPACE__. <- new.env(hash = TRUE, parent = baseenv())
  ns$.__NAMESPACE__.$spec <- c(name = name, version = "0.0.0")
  dll_list <- list()
  if (!is.null(dll_info)) {
    dll_list[["cdid"]] <- dll_info
  }
  ns$.__NAMESPACE__.$DLLs <- dll_list
  noh_populate_env_from_pkg(ns, pkg_root = pkg_root, keep.source = FALSE)
  ns
}

noh_capture_installed_cdid_dll <- function() {
  if (!isNamespaceLoaded("cdid")) return(NULL)
  ns <- asNamespace("cdid")
  dll <- tryCatch(getNamespaceInfo(ns, "DLLs")[["cdid"]], error = function(e) NULL)
  if (is.null(dll) || !inherits(dll, "DLLInfo")) return(NULL)
  dll[["path"]]
}

noh_unload_all_cdid_dlls <- function() {
  loaded <- base::getLoadedDLLs()
  for (dll in rev(loaded)) {
    if (identical(dll[["name"]], "cdid")) {
      try(dyn.unload(dll[["path"]]), silent = TRUE)
    }
  }
}

noh_restore_installed_cdid_dll <- function(orig_path) {
  noh_unload_all_cdid_dlls()
  if (!is.null(orig_path) && nzchar(orig_path) && file.exists(orig_path)) {
    dll <- dyn.load(orig_path)
    if (isNamespaceLoaded("cdid")) {
      ns <- asNamespace("cdid")
      namespace_info <- get(".__NAMESPACE__.", envir = ns, inherits = FALSE)
      namespace_info$DLLs[["cdid"]] <- dll
    }
  }
}

test_that("1. Native symbol table cache follows current owner across same-basename DLLs and verifies warm hits in C", {
  pkg_root <- noh_locate_pkg_root()
  c_src <- noh_locate_c_source()
  noh_require_sources(pkg_root, c_src)
  expect_false(is.null(pkg_root))
  expect_false(is.null(c_src))

  orig_dll_path <- noh_capture_installed_cdid_dll()
  work_dir <- noh_make_exec_dir("noh_test1_")
  on.exit({
    noh_restore_installed_cdid_dll(orig_dll_path)
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  noh_unload_all_cdid_dlls()

  c_lines <- readLines(c_src)
  so_a <- noh_compile_dll(c_lines, "cdid", work_dir)
  so_b <- noh_compile_dll(c_lines, "cdid", work_dir)
  expect_false(identical(so_a, so_b))

  dll_a <- dyn.load(so_a)
  dll_b <- dyn.load(so_b)
  expect_identical(dll_a[["name"]], "cdid")
  expect_identical(dll_b[["name"]], "cdid")
  expect_false(identical(dll_a[["info"]], dll_b[["info"]]))

  req_names <- c(
    "C_gmm_native_status",
    "C_gmm_test_arithmetic_probe",
    "C_gmm_sexp_aliases",
    "C_gmm_next",
    "C_gmm_iv_dot",
    "C_gmm_iv_norm",
    "C_gmm_iv_mat_residual",
    "C_gmm_iv_eigen_residual",
    "C_gmm_iv_mat_product"
  )

  # Direct symbol addresses from dll_a and dll_b are distinct for all 9 routines
  addrs_a <- lapply(req_names, function(nm) getNativeSymbolInfo(nm, PACKAGE = dll_a)$address)
  addrs_b <- lapply(req_names, function(nm) getNativeSymbolInfo(nm, PACKAGE = dll_b)$address)
  names(addrs_a) <- names(addrs_b) <- req_names
  for (nm in req_names) {
    expect_false(identical(addrs_a[[nm]], addrs_b[[nm]]), info = nm)
  }

  # Namespace-owner selection: two locked mock namespaces bound to dll_a and dll_b
  ns_a <- noh_make_mock_namespace("cdid_ns_a", dll_info = dll_a, pkg_root = pkg_root)
  ns_b <- noh_make_mock_namespace("cdid_ns_b", dll_info = dll_b, pkg_root = pkg_root)
  lockEnvironment(ns_a, bindings = TRUE)
  lockEnvironment(ns_b, bindings = TRUE)
  expect_true(isNamespace(ns_a) && environmentIsLocked(ns_a))
  expect_true(isNamespace(ns_b) && environmentIsLocked(ns_b))

  syms_a1 <- ns_a$gmm_require_native_dll()
  syms_b1 <- ns_b$gmm_require_native_dll()
  syms_a2 <- ns_a$gmm_require_native_dll()
  syms_b2 <- ns_b$gmm_require_native_dll()

  expect_identical(names(syms_a1), req_names)
  expect_identical(names(syms_b1), req_names)
  for (nm in req_names) {
    expect_identical(syms_a1[[nm]], addrs_a[[nm]], info = nm)
    expect_identical(syms_a2[[nm]], addrs_a[[nm]], info = nm)
    expect_identical(syms_b1[[nm]], addrs_b[[nm]], info = nm)
    expect_identical(syms_b2[[nm]], addrs_b[[nm]], info = nm)
    expect_false(identical(syms_a2[[nm]], syms_b2[[nm]]), info = nm)
  }
  expect_identical(ns_a$gmm_native_backend_info()$dll_path, dll_a[["path"]])
  expect_identical(ns_b$gmm_native_backend_info()$dll_path, dll_b[["path"]])

  # Bounded cache structure: only `owner`, `info`, `syms`
  expect_identical(sort(ls(ns_a$gmm_native_dll_cache, all.names = TRUE)), c("info", "owner", "syms"))
  expect_identical(ns_a$gmm_native_dll_cache$owner, ns_a)
  expect_identical(ns_a$gmm_native_dll_cache$info, dll_a[["info"]])

  # Source-owner selection (`isNamespace(owner) == FALSE`): follows `base::getLoadedDLLs()[["cdid"]]`
  src_env <- noh_load_source_owner_env(pkg_root)
  expect_false(isNamespace(src_env))
  syms_src_a <- src_env$gmm_require_native_dll()
  for (nm in req_names) {
    expect_identical(syms_src_a[[nm]], addrs_a[[nm]], info = nm)
  }
  # Unloading dll_a makes dll_b the active `base::getLoadedDLLs()[["cdid"]]`
  dyn.unload(so_a)
  syms_src_b <- src_env$gmm_require_native_dll()
  for (nm in req_names) {
    expect_identical(syms_src_b[[nm]], addrs_b[[nm]], info = nm)
  }
})

test_that("2. Unload/reload, missing-DLL, deparsed standalone probe, and invalid-DLL mutants refuse with typed cdid_gmm_native_error", {
  pkg_root <- noh_locate_pkg_root()
  c_src <- noh_locate_c_source()
  noh_require_sources(pkg_root, c_src)
  orig_dll_path <- noh_capture_installed_cdid_dll()
  work_dir <- noh_make_exec_dir("noh_test2_")
  on.exit({
    noh_restore_installed_cdid_dll(orig_dll_path)
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  noh_unload_all_cdid_dlls()
  src_env <- noh_load_source_owner_env(pkg_root)

  # 2a. No DLL loaded: every consumer refuses with typed cdid_gmm_native_error
  v1 <- c(1.0, 2.0); v2 <- c(3.0, 4.0)
  m1 <- matrix(c(2.0, 0.5, 0.5, 2.0), 2L, 2L)
  consumers <- list(
    resolver = function(e) e$gmm_require_native_dll(),
    info = function(e) e$gmm_native_backend_info(),
    alias = function(e) e$gmm_sexp_aliases(v1, v1),
    dot = function(e) e$gmm_iv_dot(v1, v2),
    norm = function(e) e$gmm_iv_norm(v1, v1),
    mat_product = function(e) e$gmm_iv_mat_product(m1, m1, m1, m1),
    mat_residual = function(e) e$gmm_iv_mat_residual(m1, m1, m1),
    residual = function(e) e$gmm_iv_residual(m1, v1, v2),
    eigen_bounds = function(e) e$certify_symmetric_eigen_bounds(m1),
    interrupt_check = function(e) e$ipw_owner_check_interrupt()
  )
  for (nm in names(consumers)) {
    err <- expect_error(consumers[[nm]](src_env), class = "cdid_gmm_native_error", label = nm)
    expect_s3_class(err, "cdid_gmm_accuracy_error")
    expect_match(conditionMessage(err), "Production native C backend \\(`cdid` DLL\\) is not loaded", info = nm)
  }

  # 2b. Deparsed standalone probe: evaluate deparse(gmm_require_native_dll) in an
  # isolated environment with ONLY `cdid_gmm_abort` (no cache binding, no autoload)
  standalone_env <- new.env(parent = baseenv())
  standalone_env$cdid_gmm_abort <- src_env$cdid_gmm_abort
  standalone_fn <- eval(parse(text = deparse(src_env$gmm_require_native_dll)), envir = standalone_env)
  environment(standalone_fn) <- standalone_env
  err_standalone <- expect_error(standalone_fn(), class = "cdid_gmm_native_error")
  expect_identical(
    conditionMessage(err_standalone),
    "Production native C backend (`cdid` DLL) is not loaded; automatic pure-R fallback is forbidden."
  )

  # 2c. Load two valid DLLs, warm namespace ns_a on dll_a, then unload dll_a while dll_b stays loaded
  c_lines <- readLines(c_src)
  so_a <- noh_compile_dll(c_lines, "cdid", work_dir)
  so_b <- noh_compile_dll(c_lines, "cdid", work_dir)
  dll_a <- dyn.load(so_a)
  dll_b <- dyn.load(so_b)

  # Standalone function also succeeds when a valid DLL is loaded, without cache binding
  syms_standalone <- standalone_fn()
  expect_length(syms_standalone, 9L)

  ns_a <- noh_make_mock_namespace("cdid_ns_unload", dll_info = dll_a, pkg_root = pkg_root)
  lockEnvironment(ns_a, bindings = TRUE)
  expect_length(ns_a$gmm_require_native_dll(), 9L)
  expect_true(ns_a$gmm_sexp_aliases(v1, v1))

  # Unload dll_a: ns_a still holds the stale DLLInfo object in its locked namespace,
  # and dll_b is still loaded globally. Every consumer in ns_a must refuse with
  # controlled cdid_gmm_native_error (neither crashing nor hijacking dll_b nor
  # leaking a raw simpleError).
  dyn.unload(so_a)
  for (nm in names(consumers)) {
    err_u <- expect_error(consumers[[nm]](ns_a), class = "cdid_gmm_native_error", label = nm)
    expect_s3_class(err_u, "cdid_gmm_accuracy_error")
    expect_match(conditionMessage(err_u), "Loaded `cdid` DLL is missing required registered `.Call` routines", info = nm)
  }
  expect_null(ns_a$gmm_native_dll_cache$syms)
  expect_null(ns_a$gmm_native_dll_cache$info)

  dyn.unload(so_b)
  for (nm in names(consumers)) {
    expect_error(consumers[[nm]](src_env), class = "cdid_gmm_native_error", label = nm)
  }

  # 2d. Mutant A: missing `C_gmm_iv_mat_product` registration
  mut_missing <- c_lines[!grepl("\\{\"C_gmm_iv_mat_product\"", c_lines)]
  so_missing <- noh_compile_dll(mut_missing, "cdid", work_dir)
  dll_missing <- dyn.load(so_missing)
  err_m1 <- expect_error(src_env$gmm_require_native_dll(), class = "cdid_gmm_native_error")
  expect_match(conditionMessage(err_m1), "missing required registered `.Call` routines")
  dyn.unload(so_missing)

  # 2e. Mutant B: wrong registered arity on `C_gmm_iv_dot` (2 instead of 3)
  mut_arity <- sub(
    "{\"C_gmm_iv_dot\", (DL_FUNC) &C_gmm_iv_dot, 3},",
    "{\"C_gmm_iv_dot\", (DL_FUNC) &C_gmm_iv_dot, 2},",
    c_lines, fixed = TRUE
  )
  so_arity <- noh_compile_dll(mut_arity, "cdid", work_dir)
  dll_arity <- dyn.load(so_arity)
  err_m2 <- expect_error(src_env$gmm_require_native_dll(), class = "cdid_gmm_native_error")
  expect_match(conditionMessage(err_m2), "Native symbol registration/forcing verification failed for `C_gmm_iv_dot`")
  dyn.unload(so_arity)

  # 2f. Mutant C: bad status payload (`force_symbols = FALSE`)
  mut_status <- sub(
    "SET_STRING_ELT(names, 2, Rf_mkChar(\"force_symbols\"));\n  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(1));",
    "SET_STRING_ELT(names, 2, Rf_mkChar(\"force_symbols\"));\n  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(0));",
    paste(c_lines, collapse = "\n"), fixed = TRUE
  )
  so_status <- noh_compile_dll(strsplit(mut_status, "\n", fixed = TRUE)[[1L]], "cdid", work_dir)
  dll_status <- dyn.load(so_status)
  err_m3 <- expect_error(src_env$gmm_require_native_dll(), class = "cdid_gmm_native_error")
  expect_match(conditionMessage(err_m3), "Native C backend arithmetic/registration status check failed")
  dyn.unload(so_status)

  # 2g. Reload valid DLL and verify clean recovery
  dll_ok <- dyn.load(so_a)
  expect_length(src_env$gmm_require_native_dll(), 9L)
  expect_true(src_env$gmm_sexp_aliases(v1, v1))
})

test_that("3. Live hardware FP rounding and subnormal-flush changes are caught on warm cache hits across all consumers", {
  pkg_root <- noh_locate_pkg_root()
  c_src <- noh_locate_c_source()
  ref_root <- noh_locate_reference_original()
  noh_require_sources(pkg_root, c_src)
  orig_dll_path <- noh_capture_installed_cdid_dll()
  work_dir <- noh_make_exec_dir("noh_test3_")
  on.exit({
    noh_restore_installed_cdid_dll(orig_dll_path)
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  noh_unload_all_cdid_dlls()
  so_prod <- noh_compile_dll(readLines(c_src), "cdid", work_dir)
  dll_prod <- dyn.load(so_prod)
  src_env <- noh_load_source_owner_env(pkg_root)

  # Compile a small foreign C helper to toggle hardware FP rounding mode and FZ/FTZ state
  fp_helper_c <- c(
    "#include <R.h>",
    "#include <Rinternals.h>",
    "#include <fenv.h>",
    "#include <stdint.h>",
    "SEXP C_noh_set_rmode(SEXP mode_sexp) {",
    "  int m = INTEGER(mode_sexp)[0];",
    "  int rc = -1;",
    "  if (m == 0) rc = fesetround(FE_TONEAREST);",
    "#ifdef FE_UPWARD",
    "  else if (m == 1) rc = fesetround(FE_UPWARD);",
    "#endif",
    "#ifdef FE_DOWNWARD",
    "  else if (m == 2) rc = fesetround(FE_DOWNWARD);",
    "#endif",
    "#ifdef FE_TOWARDZERO",
    "  else if (m == 3) rc = fesetround(FE_TOWARDZERO);",
    "#endif",
    "  return Rf_ScalarInteger(rc);",
    "}",
    "SEXP C_noh_set_fz(SEXP enable_sexp) {",
    "  int en = LOGICAL(enable_sexp)[0];",
    "#if defined(__aarch64__)",
    "  uint64_t fpcr = 0;",
    "  __asm__ __volatile__(\"mrs %0, fpcr\" : \"=r\"(fpcr));",
    "  if (en) fpcr |= (1ULL << 24); else fpcr &= ~(1ULL << 24);",
    "  __asm__ __volatile__(\"msr fpcr, %0\" : : \"r\"(fpcr));",
    "  return Rf_ScalarLogical(1);",
    "#elif defined(__x86_64__) || defined(__SSE2__)",
    "  unsigned int mxcsr = 0;",
    "  __asm__ __volatile__(\"stmxcsr %0\" : \"=m\"(mxcsr));",
    "  if (en) mxcsr |= 0x8040U; else mxcsr &= ~0x8040U;",
    "  __asm__ __volatile__(\"ldmxcsr %0\" : : \"m\"(mxcsr));",
    "  return Rf_ScalarLogical(1);",
    "#else",
    "  return Rf_ScalarLogical(0);",
    "#endif",
    "}"
  )
  so_fp <- noh_compile_dll(fp_helper_c, "noh_fp_helper", work_dir)
  dll_fp <- dyn.load(so_fp)
  sym_rmode <- getNativeSymbolInfo("C_noh_set_rmode", PACKAGE = dll_fp)$address
  sym_fz <- getNativeSymbolInfo("C_noh_set_fz", PACKAGE = dll_fp)$address
  on.exit({
    try(.Call(sym_rmode, 0L), silent = TRUE)
    try(.Call(sym_fz, FALSE), silent = TRUE)
    try(dyn.unload(so_fp), silent = TRUE)
  }, add = TRUE)

  # Warm the symbol cache under valid FE_TONEAREST state
  syms <- src_env$gmm_require_native_dll()
  expect_identical(src_env$gmm_native_dll_cache$info, dll_prod[["info"]])

  # An older valid backend has the same nine routines/status but an identity-only
  # alias routine. It must still enforce live FP checks when selected by its
  # namespace, even while the current candidate DLL is also loaded.
  so_old <- noh_compile_dll(readLines(file.path(ref_root, "src", "cdid.c")), "cdid", work_dir)
  dll_old <- dyn.load(so_old)
  old_ns <- noh_make_mock_namespace("cdid_old_backend", dll_info = dll_old, pkg_root = pkg_root)
  lockEnvironment(old_ns, bindings = TRUE)
  old_syms <- old_ns$gmm_require_native_dll()
  expect_identical(old_syms$C_gmm_next,
    getNativeSymbolInfo("C_gmm_next", PACKAGE = dll_old)$address)
  expect_true(old_ns$gmm_sexp_aliases(old_ns, old_ns))

  v1 <- c(1.0, 2.0); v2 <- c(3.0, 4.0)
  m1 <- matrix(c(2.0, 0.5, 0.5, 2.0), 2L, 2L)
  consumers <- list(
    resolver = function() src_env$gmm_require_native_dll(),
    info = function() src_env$gmm_native_backend_info(),
    alias = function() src_env$gmm_sexp_aliases(v1, v1),
    dot = function() src_env$gmm_iv_dot(v1, v2),
    norm = function() src_env$gmm_iv_norm(v1, v1),
    mat_product = function() src_env$gmm_iv_mat_product(m1, m1, m1, m1),
    mat_residual = function() src_env$gmm_iv_mat_residual(m1, m1, m1),
    residual = function() src_env$gmm_iv_residual(m1, v1, v2),
    eigen_bounds = function() src_env$certify_symmetric_eigen_bounds(m1),
    interrupt_check = function() src_env$ipw_owner_check_interrupt()
  )
  consumers$original_backend_resolver <- function() old_ns$gmm_require_native_dll()
  consumers$original_backend_alias <- function() old_ns$gmm_sexp_aliases(v1, v1)
  consumers$original_backend_interrupt <- function() old_ns$ipw_owner_check_interrupt()

  # Toggle FE_UPWARD (1L), FE_DOWNWARD (2L), and FE_TOWARDZERO (3L) while cache is warm
  rm_guard_re <- paste0(
    "(FE_TONEAREST is required|",
    "Gradual underflow check failed \\(FTZ/DAZ or subnormal flushing detected\\)\\.|",
    "Round-to-nearest ties-to-even \\(RNTE\\) addition/subtraction probe failed\\.)"
  )
  for (rm_code in c(1L, 2L, 3L)) {
    rc <- .Call(sym_rmode, rm_code)
    if (identical(rc, 0L)) {
      for (nm in names(consumers)) {
        err_fp <- expect_error(consumers[[nm]](), class = "cdid_gmm_native_error", label = paste0("rm_", rm_code, "_", nm))
        expect_match(conditionMessage(err_fp), rm_guard_re, info = nm)
      }
      .Call(sym_rmode, 0L)
      expect_true(src_env$gmm_sexp_aliases(v1, v1))
    }
  }

  # Toggle hardware subnormal flushing (FZ / FTZ+DAZ) while cache is warm
  fz_supported <- isTRUE(.Call(sym_fz, TRUE))
  if (fz_supported) {
    for (nm in names(consumers)) {
      err_fz <- expect_error(consumers[[nm]](), class = "cdid_gmm_native_error", label = paste0("fz_", nm))
      expect_match(conditionMessage(err_fz), "Gradual underflow check failed", info = nm)
    }
    .Call(sym_fz, FALSE)
    expect_true(src_env$gmm_sexp_aliases(v1, v1))
  }

  # Also verify `C_gmm_test_arithmetic_probe` synthetic fault detection
  min_sub <- 2^-1074
  odd_tie <- 1.0 + 2.0 * (2^-52)
  expect_true(.Call(syms$C_gmm_test_arithmetic_probe, -1L, min_sub, 1.0, odd_tie, 0.0))
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, 0.0, 1.0, odd_tie, 0.0), class = "cdid_gmm_native_error")
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, min_sub, 1.0 + 2^-52, odd_tie, 0.0), class = "cdid_gmm_native_error")
  expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, min_sub, 1.0, odd_tie, 2^-104), class = "cdid_gmm_native_error")
})

test_that("4. Owner seal identity, mutated/mismatched owner refusals, ALTREP compactness, and outcome filtering are preserved", {
  pkg_root <- noh_locate_pkg_root()
  c_src <- noh_locate_c_source()
  noh_require_sources(pkg_root, c_src)
  orig_dll_path <- noh_capture_installed_cdid_dll()
  work_dir <- noh_make_exec_dir("noh_test4_")
  on.exit({
    noh_restore_installed_cdid_dll(orig_dll_path)
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  noh_unload_all_cdid_dlls()
  so_prod <- noh_compile_dll(readLines(c_src), "cdid", work_dir)
  dll_prod <- dyn.load(so_prod)
  ns <- noh_make_mock_namespace("cdid", dll_info = dll_prod, pkg_root = pkg_root)
  lockEnvironment(ns, bindings = TRUE)

  # ALTREP compactness preserved by `gmm_sexp_aliases`
  alt_seq <- seq_len(2^31 + 5)
  before_insp <- capture.output(.Internal(inspect(alt_seq, 0L, 1L)))[1L]
  expect_true(ns$gmm_sexp_aliases(alt_seq, alt_seq))
  after_insp <- capture.output(.Internal(inspect(alt_seq, 0L, 1L)))[1L]
  expect_match(before_insp, "compact")
  expect_match(after_insp, "compact")

  ctrl <- ns$ipw_owner_control(
    max_chunk_rows = 3L,
    max_retained_entries = 100000L,
    max_matrix_entries = 100000L,
    max_work_ops = 500000L,
    max_output_entries = 100000L
  )
  df <- rbind(
    data.frame(id = c(1, 1, 1), time = c(1, 2, 3), y = c(1.0, 2.0, 3.5), g = c(2, 2, 2)),
    data.frame(id = c(2, 2, 2), time = c(1, 2, 3), y = c(1.5, NA_real_, 4.0), g = c(2, 2, 2)),
    data.frame(id = c(3, 3, 3), time = c(1, 2, 3), y = c(999.0, -999.0, 2.5), g = c(3, 3, 3)),
    data.frame(id = c(4, 4, 4), time = c(1, 2, 3), y = c(2.0, 3.0, 4.5), g = c(Inf, Inf, Inf))
  )
  p1 <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl)
  p2 <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl)
  Z <- cbind(1.0, c(0.1, 0.2, 0.3, 0.4))
  f1 <- ns$construct_ipw_frame_owner(p1, id = 1:4, Z = Z, control = ctrl)
  f2 <- ns$construct_ipw_frame_owner(p1, id = 1:4, Z = Z, control = ctrl)
  nu1 <- ns$construct_ipw_nuisance_view(f1, g = 2, K = c(Inf), s = 1, t = 2, control = ctrl)
  nu2 <- ns$construct_ipw_nuisance_view(f1, g = 2, K = c(Inf), s = 1, t = 2, control = ctrl)
  lk1 <- ns$construct_ipw_link_view(f1, s = 1, t = 2, nuisance_view = nu1, control = ctrl)
  lk_unverified <- ns$construct_ipw_link_view(f1, s = 1, t = 2, nuisance_view = nu1, control = ctrl, verify_nuisance_mask = FALSE)

  # Valid owner/view validations succeed
  expect_silent(ns$validate_ipw_frame_owner(f1, panel_owner = p1))
  expect_silent(ns$validate_ipw_nuisance_view(nu1, frame_owner = f1))
  expect_silent(ns$validate_ipw_link_view(lk1, frame_owner = f1, nuisance_view = nu1))

  # Mutated copies and mismatched valid snapshots are refused
  p1_mut <- p1; p1_mut$unit_id_key <- as.double(p1$unit_id_key) + 0
  f1_mut <- f1; f1_mut$w_a <- as.double(f1$w_a) + 0
  nu1_mut <- nu1; nu1_mut$defining_pair <- as.list(nu1$defining_pair)
  expect_error(ns$validate_ipw_frame_owner(f1, panel_owner = p1_mut), class = "cdid_ipw_owner_integrity_error")
  expect_error(ns$validate_ipw_frame_owner(f1, panel_owner = p2), class = "cdid_ipw_owner_mismatch_error")
  expect_error(ns$validate_ipw_nuisance_view(nu1, frame_owner = f1_mut), class = "cdid_ipw_owner_integrity_error")
  expect_error(ns$validate_ipw_nuisance_view(nu1, frame_owner = f2), class = "cdid_ipw_owner_mismatch_error")
  expect_error(ns$validate_ipw_link_view(lk1, frame_owner = f1_mut), class = "cdid_ipw_owner_integrity_error")
  expect_error(ns$validate_ipw_link_view(lk1, frame_owner = f2), class = "cdid_ipw_owner_mismatch_error")
  expect_error(ns$validate_ipw_link_view(lk1, frame_owner = f1, nuisance_view = nu1_mut), class = "cdid_ipw_owner_integrity_error")
  expect_error(ns$validate_ipw_link_view(lk1, frame_owner = f1, nuisance_view = nu2), class = "cdid_ipw_owner_mismatch_error")
  expect_error(ns$validate_ipw_link_view(lk_unverified, frame_owner = f1, nuisance_view = nu1), class = "cdid_ipw_owner_mask_unverified_error")

  # Bound refusal when queried chunk length (4) exceeds `max_chunk_rows` (3)
  expect_error(
    ns$query_ipw_view_chunk(lk1, 1L, 4L, ctrl, read_outcomes = FALSE),
    class = "cdid_ipw_owner_bound_error"
  )

  # Outcome read filtering: unit 3 (g = 3, outside risk set {2, Inf}) has y = (999, -999) at (1, 2)
  # and must NEVER be read when querying or streaming `lk1`
  ctrl4 <- ns$ipw_owner_control(
    max_chunk_rows = 4L,
    max_retained_entries = 100000L,
    max_matrix_entries = 100000L,
    max_work_ops = 500000L,
    max_output_entries = 100000L
  )
  tr <- new.env(parent = emptyenv())
  ch_no_out <- ns$query_ipw_view_chunk(lk1, 1L, 4L, ctrl4, read_outcomes = FALSE, trace_env = tr)
  expect_identical(tr$calls, 0L)
  expect_null(ch_no_out$observed_risk$y_left_obs)

  ch_with_out <- ns$query_ipw_view_chunk(lk1, 1L, 4L, ctrl4, read_outcomes = TRUE, enclose_differences = TRUE, trace_env = tr)
  expect_identical(ch_with_out$observed_risk$a_index, c(1L, 4L))
  expect_identical(ch_with_out$observed_risk$y_left_obs, c(1.0, 2.0))
  expect_identical(ch_with_out$observed_risk$y_right_obs, c(2.0, 3.0))
  expect_false(any(c(7L, 8L) %in% tr$pos_s_index))

  # Multi-batch stream with `max_chunk_rows = 3L` (batches [1,3] and [4,4])
  streamed_obs <- integer(0L)
  st_summary <- ns$stream_ipw_view_chunks(
    lk1, control = ctrl, read_outcomes = TRUE,
    callback = function(ch) streamed_obs <<- c(streamed_obs, ch$observed_risk$a_index)
  )
  expect_identical(st_summary$n_batches, 2L)
  expect_identical(streamed_obs, c(1L, 4L))
})

test_that("5. Actual OS SIGINT in throwaway child processes is caught by C-backed ipw_owner_check_interrupt and stream_ipw_view_chunks", {
  pkg_root <- noh_locate_pkg_root()
  c_src <- noh_locate_c_source()
  noh_require_sources(pkg_root, c_src)
  if (identical(.Platform$OS.type, "windows")) {
    testthat::skip("POSIX SIGINT child probe requires forked parallel::mcparallel/mccollect and POSIX SIGINT delivery, which are unavailable on Windows.")
  }
  orig_dll_path <- noh_capture_installed_cdid_dll()
  work_dir <- noh_make_exec_dir("noh_test5_")
  on.exit({
    noh_restore_installed_cdid_dll(orig_dll_path)
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  so_prod <- noh_compile_dll(readLines(c_src), "cdid", work_dir)
  rscript_bin <- file.path(R.home("bin"), "Rscript")

  run_child_sigint_probe <- function(mode_label, loop_lines) {
    ready_file <- file.path(work_dir, paste0("child_ready_", mode_label, ".txt"))
    result_file <- file.path(work_dir, paste0("child_result_", mode_label, ".txt"))
    script_file <- file.path(work_dir, paste0("child_worker_", mode_label, ".R"))
    stdout_file <- file.path(work_dir, paste0("child_stdout_", mode_label, ".txt"))
    stderr_file <- file.path(work_dir, paste0("child_stderr_", mode_label, ".txt"))
    job <- NULL
    child_exit <- NULL
    child_pending <- function() {
      if (is.null(job) || !is.null(child_exit)) return(FALSE)
      done <- parallel::mccollect(job, wait = FALSE)
      if (!is.null(done)) {
        child_exit <<- done
        return(FALSE)
      }
      TRUE
    }
    child_output <- function() {
      paste(vapply(c(stdout_file, stderr_file), function(p) {
        paste0(basename(p), ":\n", if (file.exists(p)) paste(readLines(p, warn = FALSE), collapse = "\n") else "<absent>")
      }, character(1L)), collapse = "\n")
    }
    on.exit({
      if (child_pending()) {
        tools::pskill(job$pid, tools::SIGTERM)
        deadline <- Sys.time() + 2
        while (child_pending() && Sys.time() < deadline) Sys.sleep(0.02)
        if (child_pending()) {
          tools::pskill(job$pid, tools::SIGKILL)
          deadline <- Sys.time() + 2
          while (child_pending() && Sys.time() < deadline) Sys.sleep(0.02)
        }
      }
    }, add = TRUE)

    writeLines(c(
      sprintf("pkg_root <- %s", deparse(pkg_root)),
      sprintf("so_prod <- %s", deparse(so_prod)),
      sprintf("ready_file <- %s", deparse(ready_file)),
      sprintf("result_file <- %s", deparse(result_file)),
      "dyn.load(so_prod)",
      "for (p in sort(list.files(file.path(pkg_root, 'R'), pattern = '\\\\.[rR]$', full.names = TRUE))) {",
      "  sys.source(p, envir = globalenv(), keep.source = FALSE)",
      "}",
      loop_lines
    ), script_file)

    # Use R's base parallel child handle rather than an asynchronous shell.
    # mc.set.seed = FALSE preserves the parent's RNG state. Each child unloads
    # inherited cdid DLLs before loading the one standalone backend under test.
    job <- parallel::mcparallel(local({
      out <- file(stdout_file, open = "wt")
      err <- file(stderr_file, open = "wt")
      sink(out)
      sink(err, type = "message")
      on.exit({
        sink(type = "message")
        sink()
        close(err)
        close(out)
      }, add = TRUE)
      noh_unload_all_cdid_dlls()
      sys.source(script_file, envir = globalenv(), keep.source = FALSE)
      TRUE
    }), mc.set.seed = FALSE)
    child_pid <- job$pid
    deadline <- Sys.time() + 15
    while (!file.exists(ready_file) && child_pending() && Sys.time() < deadline) Sys.sleep(0.02)
    expect_true(file.exists(ready_file), label = mode_label)
    if (!file.exists(ready_file)) stop("Child did not become ready.\n", child_output(), call. = FALSE)
    expect_true(is.finite(child_pid) && child_pid > 1L && child_pid != Sys.getpid(), label = mode_label)
    expect_identical(as.integer(trimws(readLines(ready_file, warn = FALSE)[1L])), child_pid)
    if (!child_pending()) stop("Ready child already exited.\n", child_output(), call. = FALSE)

    # Send exactly one SIGINT only to this identity-checked throwaway child.
    expect_true(tools::pskill(child_pid, tools::SIGINT), label = mode_label)
    deadline_res <- Sys.time() + 15
    while (!file.exists(result_file) && child_pending() && Sys.time() < deadline_res) Sys.sleep(0.02)
    expect_true(file.exists(result_file), label = mode_label)
    if (!file.exists(result_file)) stop("Child produced no interrupt result.\n", child_output(), call. = FALSE)
    while (child_pending() && Sys.time() < deadline_res) Sys.sleep(0.02)
    expect_false(child_pending(), label = mode_label)
    expect_identical(child_exit[[1L]], TRUE)
    readLines(result_file, warn = FALSE)[1L]
  }

  # Probe A: Direct `ipw_owner_check_interrupt()` loop
  out_direct <- run_child_sigint_probe("direct", c(
    "ipw_owner_check_interrupt()",
    "iters <- 0L",
    "res <- tryCatch({",
    "  writeLines(as.character(Sys.getpid()), ready_file)",
    "  repeat { iters <- iters + 1L; ipw_owner_check_interrupt() }",
    "  'NO_INTERRUPT'",
    "}, interrupt = function(e) {",
    "  paste0('CAUGHT_INTERRUPT:', paste(class(e), collapse = ','), ':iters=', iters)",
    "})",
    "writeLines(res, result_file)"
  ))
  expect_match(out_direct, "^CAUGHT_INTERRUPT:interrupt,condition:iters=[0-9]+$")

  # Probe B: Active `stream_ipw_view_chunks(...)` traversal
  out_stream <- run_child_sigint_probe("stream", c(
    "ctrl <- ipw_owner_control(max_chunk_rows = 2L, max_retained_entries = 100000L, max_matrix_entries = 100000L, max_work_ops = 500000L, max_output_entries = 100000L)",
    "df <- rbind(",
    "  data.frame(id = c(1,1), time = c(1,2), y = c(1.0,2.0), g = c(2,2)),",
    "  data.frame(id = c(2,2), time = c(1,2), y = c(1.5,2.5), g = c(Inf,Inf)),",
    "  data.frame(id = c(3,3), time = c(1,2), y = c(2.0,3.5), g = c(2,2)),",
    "  data.frame(id = c(4,4), time = c(1,2), y = c(0.5,1.5), g = c(Inf,Inf))",
    ")",
    "p <- construct_ipw_panel_owner(df, 'id', 'time', 'y', 'g', control = ctrl)",
    "f <- construct_ipw_frame_owner(p, id = 1:4, Z = cbind(1.0, 1:4), control = ctrl)",
    "lk <- construct_ipw_link_view(f, s = 1, t = 2, g = 2, K = Inf, control = ctrl)",
    "iters <- 0L",
    "res <- tryCatch({",
    "  writeLines(as.character(Sys.getpid()), ready_file)",
    "  repeat {",
    "    iters <- iters + 1L",
    "    stream_ipw_view_chunks(lk, control = ctrl, read_outcomes = TRUE, callback = function(ch) NULL)",
    "  }",
    "  'NO_INTERRUPT'",
    "}, interrupt = function(e) {",
    "  paste0('CAUGHT_STREAM_INTERRUPT:', paste(class(e), collapse = ','), ':iters=', iters)",
    "})",
    "writeLines(res, result_file)"
  ))
  expect_match(out_stream, "^CAUGHT_STREAM_INTERRUPT:interrupt,condition:iters=[0-9]+$")
})

test_that("6. Structural delta and bit-identical equivalence against bundled original", {
  pkg_root <- noh_locate_pkg_root()
  c_src <- noh_locate_c_source()
  ref_root <- noh_locate_reference_original()
  noh_require_sources(pkg_root, c_src)

  orig_dll_path <- noh_capture_installed_cdid_dll()
  work_dir <- noh_make_exec_dir("noh_test6_")
  on.exit({
    noh_restore_installed_cdid_dll(orig_dll_path)
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  # 6a. Structural checks across both keep.source modes
  for (ks in c(FALSE, TRUE)) {
    e_orig_core <- new.env(parent = baseenv())
    e_cand_core <- new.env(parent = baseenv())
    sys.source(file.path(ref_root, "R", "chained_gmm_core.R"), envir = e_orig_core, keep.source = ks)
    sys.source(file.path(pkg_root, "R", "chained_gmm_core.R"), envir = e_cand_core, keep.source = ks)

    orig_core_names <- ls(e_orig_core, all.names = TRUE)
    cand_core_names <- ls(e_cand_core, all.names = TRUE)
    # The cache-stage contract remains historical. The separately accepted
    # CP08 GMM unit accumulator adds exactly these three helpers (1421031):
    # gmm_unit_confirmation.json and sharedX_unit_integration.json, with the
    # latter's verified corrected_source_pins.json binding this exact core.
    cp08_unit_core_added_bindings <- c("gmm_if_unit_batch", "gmm_if_unit_batch_eligible",
      "gmm_require_unit_batch_symbol")
    expect_identical(setdiff(cand_core_names, orig_core_names),
      sort(c("gmm_native_dll_cache", cp08_unit_core_added_bindings)))
    expect_identical(setdiff(orig_core_names, cand_core_names), character(0L))
    # D28 adds these reviewed controls and changes these six core bodies
    # (the solve alias shares one body). The original native-cache delta and
    # the accepted CP03 gauge-body delta remain separate. Pin reviewed bytes;
    # test-chained-gmm-weighting.R supplies independent fixed-matrix references.
    d28_formal_additions <- list(
      chained_gmm_core = list(gmm_covariance = NULL, covariance_tuning = NULL, reduction_spec = NULL),
      gmm_accuracy_control = list(covariance_strategy = NULL),
      preflight_gmm_allocations = list(covariance_strategy = NULL, n_dropped_triplets = 0L),
      solve_chained_gmm_core = list(gmm_covariance = NULL, covariance_tuning = NULL, reduction_spec = NULL))
    d28_changed_bodies <- c("chained_gmm_core", "gmm_accuracy_control", "preflight_gmm_allocations",
      "solve_chained_gmm_core", "validate_gmm_accuracy_control", "validate_gmm_reduction_spec")
    # CP03: multiple anchors, zero-loading rows and intrinsic validation;
    # only this additional body changed, with no additional formal changes.
    cp03_changed_bodies <- "build_gmm_gauge"
    # CP08 changes duplicate-ID bookkeeping and, in accepted 1421031, routes
    # eligible unweighted unit replays through the two-pass native accumulator.
    # The exact changed-body set is unchanged; the current pin binds both deltas.
    cp08_changed_bodies <- "project_gmm_if_stream"
    actual_changed_formals <- orig_core_names[vapply(orig_core_names, function(nm)
      !identical(formals(get(nm, e_cand_core)), formals(get(nm, e_orig_core))), logical(1))]
    actual_changed_bodies <- orig_core_names[vapply(orig_core_names, function(nm)
      !identical(body(utils::removeSource(get(nm, e_cand_core))),
        body(utils::removeSource(get(nm, e_orig_core)))), logical(1))]
    expect_identical(actual_changed_formals, names(d28_formal_additions))
    expect_identical(actual_changed_bodies,
      sort(c("gmm_require_native_dll", d28_changed_bodies, cp03_changed_bodies, cp08_changed_bodies)))
    # Historical D28-only core SHA256 (before the accepted CP03 gauge delta):
    # 5f1fd263e702fecfe8ccc04e9005d83a3e67ce265bc5be418df7fe38e0e23b7e
    # Historical CP08 cache-stage whole-core pin (not the current source):
    # cb0b987738d7629d7cb5eedc8db0020fd35f2c5ce444d1efdb92a2f9b1f8bb34
    # Earlier CP03 acceptance remains in cp03_backend_acceptance.json.
    # Historical sharedX_unit_integration.json pin:
    # 93adc0e289ba55f84b4cd8c1c79d8127f0d5f95d58f423cc6f61146cbaa02a38.
    # Current accepted pin (post_validation_fixes_20261029/CONFIRMED_PINS.json):
    expect_identical(digest::digest(file = file.path(pkg_root, "R", "chained_gmm_core.R"), algo = "sha256"),
      "c1b5f022ac8d7771c5ef6eb8c1e0b1148182ab244668408cb124090ffe74ea3f")
    for (nm in orig_core_names) {
      vo <- get(nm, envir = e_orig_core, inherits = FALSE)
      vc <- get(nm, envir = e_cand_core, inherits = FALSE)
      additions <- d28_formal_additions[[nm]]
      if (is.null(additions)) {
        expect_identical(formals(vc), formals(vo), info = paste0("core_formals_", ks, "_", nm))
      } else {
        expect_identical(names(formals(vc)), c(names(formals(vo)), names(additions)),
          info = paste0("core_added_formals_", ks, "_", nm))
        expect_identical(as.pairlist(formals(vc)[names(formals(vo))]), formals(vo),
          info = paste0("core_legacy_formals_", ks, "_", nm))
        expect_identical(formals(vc)[names(additions)], additions,
          info = paste0("core_added_defaults_", ks, "_", nm))
      }
      if (!(nm %in% c("gmm_require_native_dll", d28_changed_bodies, cp03_changed_bodies, cp08_changed_bodies))) {
        expect_identical(body(vc), body(vo), info = paste0("core_body_", ks, "_", nm))
      }
    }

    e_orig_own <- new.env(parent = baseenv())
    e_cand_own <- new.env(parent = baseenv())
    sys.source(file.path(ref_root, "R", "chained_ipw_owner.R"), envir = e_orig_own, keep.source = ks)
    sys.source(file.path(pkg_root, "R", "chained_ipw_owner.R"), envir = e_cand_own, keep.source = ks)

    orig_own_names <- ls(e_orig_own, all.names = TRUE)
    cand_own_names <- ls(e_cand_own, all.names = TRUE)
    # Accepted CP08 P1 adds exactly four mask helpers and one constant. Equal
    # fingerprints still require an exact mask comparison; no generic
    # allowance for added bindings or changed legacy formals/bodies is made.
    p1_owner_added_bindings <- c("fingerprint_ipw_view_mask", "ipw_owner_exact_mod_sum",
      "ipw_owner_mask_fingerprint_modulus", "ipw_owner_mask_mismatch_message", "ipw_owner_pair_mask_fingerprint")
    # Exact named CP08 delta, independently reviewed with installed references.
    # Whole-file pins below bind the bodies of every added/changed helper.
    cp08_owner_added_bindings <- c(
      "cdid_dependency_check_done", "cdid_dependency_check_record", "cdid_dependency_checks",
      "ipw_owner_construction_context", "ipw_owner_construction_context_for",
      "ipw_owner_context_fingerprint_key", "ipw_owner_context_link_support", "ipw_owner_context_record_nuisance",
      "ipw_owner_date_complete", "ipw_owner_pair_complete", "ipw_owner_pair_masks_equal_fact",
      "ipw_owner_posting_date_complete", "ipw_owner_query_complete_pair_mask", "ipw_owner_sorted_lookup_vectorizable",
      "ipw_owner_validated_control", "ipw_owner_validated_preflight", "ipw_owner_validated_view_match",
      "ipw_owner_validated_view_record")
    # CP08 stage 5: nuisance-free link support records of one CP04 registry.
    cp08_stage5_owner_added_bindings <- c("ipw_owner_context_link_admits", "ipw_owner_context_link_key",
      "ipw_owner_context_link_lookup", "ipw_owner_context_record_link", "ipw_owner_link_support_context",
      "ipw_owner_pair_counts_full")
    # Later accepted CP08 source, bound by common_context_integration.json's
    # exact owner pin and confirmation (last owner change 734983c). These are
    # explicit deltas beyond stage 5, not a wildcard for new private helpers.
    # Prepared intake shares header validation and sealed representation with
    # the standalone panel/frame constructors; public formals are unchanged.
    cp08_prepared_owner_added_bindings <- c("construct_ipw_frame_owner_prepared",
      "construct_ipw_panel_owner_prepared", "ipw_frame_owner_assemble", "ipw_frame_owner_headers",
      "ipw_panel_owner_assemble", "ipw_panel_owner_headers")
    # Checked sets retain exact immutable nuisance snapshots within one call;
    # the public link validator still invokes the full body with a NULL set.
    cp08_checked_owner_added_bindings <- c("ipw_link_view_nuisance_checked_add",
      "ipw_link_view_nuisance_checked_has", "ipw_link_view_nuisance_checked_set",
      "ipw_link_view_validate_body", "ipw_owner_immutable_view_snapshot")
    # Constructor context reuse rechecks seals/owners; the mask-fact helper
    # retains owner, work-bound and ordered-cohort guards before posting facts.
    cp08_context_owner_added_bindings <- c("ipw_owner_context_check_frame",
      "ipw_owner_context_frame_checked", "ipw_owner_context_nuisance_checked",
      "ipw_view_mask_compatibility_fact")
    expect_identical(cand_own_names, sort(c(orig_own_names, p1_owner_added_bindings, cp08_owner_added_bindings,
      cp08_stage5_owner_added_bindings, cp08_prepared_owner_added_bindings,
      cp08_checked_owner_added_bindings, cp08_context_owner_added_bindings)))
    expect_identical(e_cand_own$ipw_owner_mask_fingerprint_modulus, 67108859)
    modified_owner_fns <- c(
      "ipw_owner_check_interrupt",
      "validate_ipw_frame_owner",
      "validate_ipw_nuisance_view",
      "validate_ipw_link_view",
      "validate_common_frame_owner"
    )
    # CP04: these three constructor aliases add only the opt-in A-only design
    # validation control. Pin the accepted owner bytes; the installed CP04
    # finite-reference tests cover the TRUE route, and 6b retains the baseline.
    cp04_owner_formal_additions <- list(
      chained_ipw_frame_owner = list(design_on_A_only = FALSE),
      construct_common_frame_owner = list(design_on_A_only = FALSE),
      construct_ipw_frame_owner = list(design_on_A_only = FALSE))
    cp04_owner_changed_bodies <- names(cp04_owner_formal_additions)
    cp08_owner_formal_additions <- list(
      construct_ipw_link_view = list(construction_context = NULL),
      construct_ipw_nuisance_view = list(construction_context = NULL),
      declare_ipw_link_view = list(construction_context = NULL),
      declare_ipw_nuisance_view = list(construction_context = NULL),
      query_ipw_link_view_chunk = list(validated = NULL),
      query_ipw_nuisance_view_chunk = list(validated = NULL),
      query_ipw_view_chunk = list(validated = NULL),
      read_ipw_view_chunk = list(validated = NULL),
      verify_ipw_view_mask_compatibility = list(construction_context = NULL))
    owner_formal_additions <- c(cp04_owner_formal_additions, cp08_owner_formal_additions)
    # Panel constructor/header/assembler extraction changes these three alias
    # bodies only. Frame aliases, view constructors, mask verifier and link
    # validator already occur in the exact historical changed-body groups;
    # their later prepared/context/checked-set deltas are bound by the new pin.
    cp08_prepared_owner_changed_bodies <- c("chained_ipw_panel_owner",
      "construct_common_panel_owner", "construct_ipw_panel_owner")
    cp08_owner_changed_bodies <- c(setdiff(names(cp08_owner_formal_additions), "verify_ipw_view_mask_compatibility"),
      "ipw_owner_match_sorted_double", "ipw_owner_require_deps", cp08_prepared_owner_changed_bodies)
    actual_changed_formals <- orig_own_names[vapply(orig_own_names, function(nm)
      !identical(formals(get(nm, e_cand_own)), formals(get(nm, e_orig_own))), logical(1))]
    expect_identical(actual_changed_formals, sort(names(owner_formal_additions)))
    # The earlier native-cache validator/interrupt and CP04 constructor deltas
    # remain exact. P1 changes only the verify-mask body's literal message to
    # its helper call among existing bindings; no existing formal changes.
    # Complete-support rate rule changes only ipw_owner_chunked_pair_support to
    # set rates[pair_count == coh_desc$group_a_count] <- 1 when require_support = TRUE.
    p1_owner_changed_bodies <- "verify_ipw_view_mask_compatibility"
    cp08_full_support_owner_changed_bodies <- "ipw_owner_chunked_pair_support"
    actual_changed_bodies <- orig_own_names[vapply(orig_own_names, function(nm)
      !identical(body(utils::removeSource(get(nm, e_cand_own))),
        body(utils::removeSource(get(nm, e_orig_own)))), logical(1))]
    expect_identical(actual_changed_bodies,
      sort(c(modified_owner_fns, cp04_owner_changed_bodies, p1_owner_changed_bodies,
        cp08_owner_changed_bodies, cp08_full_support_owner_changed_bodies)))
    # Historical owner pin includes the CP08 reviewed optimization delta and
    # stage-5 link support records:
    # 5b0419b32155244b7aa37a3197a007fac58f7b7e71cbf0448905b550c787ea01
    # (accepted stage-4 owner SHA:
    # 39929ee64e80584894e16c094e4647fea997dba2f4b4d070f3c337233d2b8ca8).
    # Historical P1 bytes: accepted_P1.json, status P1_ACCEPTED_EXACT_SOURCE_INTEGRATED,
    # files$R/chained_ipw_owner.R (P23 leaves this owner unchanged). Historical
    # accepted CP04 owner SHA: 518a9b09057389dbdb03f0be5f01f97b07630595ce66cd1612ffd8aaa05071c8.
    # Historical common_context_integration.json owner SHA:
    # 8b27a3b1ceac5643731803f657943d45ae9e1112d8141247102ca9741e57fa63.
    # Current accepted owner pin including the complete-support exact rate = 1 delta:
    expect_identical(digest::digest(file = file.path(pkg_root, "R", "chained_ipw_owner.R"), algo = "sha256"),
      "b68487d7d16cc3546e6a4beb4b2f0cdf9e53cee6526691ba46405e4db1260e78")
    for (nm in orig_own_names) {
      vo <- get(nm, envir = e_orig_own, inherits = FALSE)
      vc <- get(nm, envir = e_cand_own, inherits = FALSE)
      if (is.function(vo)) {
        additions <- owner_formal_additions[[nm]]
        if (is.null(additions)) {
          expect_identical(formals(vc), formals(vo), info = paste0("own_formals_", ks, "_", nm))
        } else {
          expect_identical(names(formals(vc)), c(names(formals(vo)), names(additions)),
            info = paste0("own_added_formals_", ks, "_", nm))
          expect_identical(as.pairlist(formals(vc)[names(formals(vo))]), formals(vo),
            info = paste0("own_legacy_formals_", ks, "_", nm))
          expect_identical(formals(vc)[names(additions)], additions,
            info = paste0("own_added_defaults_", ks, "_", nm))
        }
        if (nm %in% p1_owner_changed_bodies) {
          # The message is unchanged. CP08 also adds the exact posting-fact
          # branch, whose complete body is bound by the reviewed owner hash.
          message <- "Link risk-set pair-observation mask differs from the nuisance view mask on chunk; reuse is refused."
          expect_identical(e_cand_own$ipw_owner_mask_mismatch_message(), message)
        } else if (nm %in% cp08_full_support_owner_changed_bodies) {
          b_cand <- as.list(body(utils::removeSource(vc)))
          b_orig <- as.list(body(utils::removeSource(vo)))
          rate_expr <- parse(
            text = "rates <- ifelse(coh_desc$group_a_mass > 0, pair_mass / coh_desc$group_a_mass, NA_real_)",
            keep.source = FALSE
          )[[1L]]
          delta_expr <- parse(
            text = "if (isTRUE(require_support)) {\n  rates[pair_count == coh_desc$group_a_count] <- 1\n  ops <- ops + 4 * m\n}",
            keep.source = FALSE
          )[[1L]]
          idx_delta <- which(vapply(b_cand, identical, logical(1L), rate_expr)) + 1L
          expect_identical(length(idx_delta), 1L, info = paste0("own_body_rate_idx_", ks, "_", nm))
          expect_identical(b_cand[[idx_delta]], delta_expr, info = paste0("own_body_delta_", ks, "_", nm))
          expect_identical(b_cand[-idx_delta], b_orig, info = paste0("own_body_rest_", ks, "_", nm))
        } else if (!(nm %in% c(modified_owner_fns, cp04_owner_changed_bodies, cp08_owner_changed_bodies))) {
          expect_identical(body(vc), body(vo), info = paste0("own_body_", ks, "_", nm))
        }
      } else {
        expect_identical(vc, vo, info = paste0("own_val_", ks, "_", nm))
      }
    }
  }

  # 6b. Bit-for-bit numerical equivalence between original and candidate C + R implementations
  noh_unload_all_cdid_dlls()
  so_orig <- noh_compile_dll(readLines(file.path(ref_root, "src", "cdid.c")), "cdid", work_dir)
  so_cand <- noh_compile_dll(readLines(c_src), "cdid", work_dir)
  dll_orig <- dyn.load(so_orig)
  dll_cand <- dyn.load(so_cand)

  ns_orig <- noh_make_mock_namespace("cdid_orig", dll_info = dll_orig, pkg_root = pkg_root)
  sys.source(file.path(ref_root, "R", "chained_gmm_core.R"), envir = ns_orig, keep.source = FALSE)
  sys.source(file.path(ref_root, "R", "chained_ipw_owner.R"), envir = ns_orig, keep.source = FALSE)
  lockEnvironment(ns_orig, bindings = TRUE)

  ns_cand <- noh_make_mock_namespace("cdid_cand", dll_info = dll_cand, pkg_root = pkg_root)
  lockEnvironment(ns_cand, bindings = TRUE)

  # Compare interval primitives and matrix products bit-for-bit
  xlo <- matrix(c(1, -2.5, 2^-40, 0, -0, 3.125), 2L, 3L)
  xhi <- xlo + matrix(c(2^-50, 2^-45, 2^-40, 0, 0, 2^-48), 2L, 3L)
  ylo <- matrix(c(0.5, -1.25, 2.0, 0, 1/3, -1/7), 3L, 2L)
  yhi <- ylo + matrix(c(2^-48, 2^-48, 2^-48, 0, 2^-50, 2^-50), 3L, 2L)

  prod_o <- ns_orig$gmm_iv_mat_product(xlo, xhi, ylo, yhi)
  prod_c <- ns_cand$gmm_iv_mat_product(xlo, xhi, ylo, yhi)
  expect_identical(noh_raw_bits(prod_c$lo), noh_raw_bits(prod_o$lo))
  expect_identical(noh_raw_bits(prod_c$hi), noh_raw_bits(prod_o$hi))

  # Compare owner/view chunk queries and enclosures bit-for-bit
  ctrl_o <- ns_orig$ipw_owner_control(4L, 100000L, 100000L, 500000L, 100000L)
  ctrl_c <- ns_cand$ipw_owner_control(4L, 100000L, 100000L, 500000L, 100000L)
  df <- rbind(
    data.frame(id = c(1, 1, 1), time = c(1, 2, 3), y = c(1.1, 2.3, 3.7), g = c(2, 2, 2)),
    data.frame(id = c(2, 2, 2), time = c(1, 2, 3), y = c(0.4, 1.9, 2.8), g = c(2, 2, 2)),
    data.frame(id = c(3, 3, 3), time = c(1, 2, 3), y = c(1.5, 2.1, 3.2), g = c(Inf, Inf, Inf)),
    data.frame(id = c(4, 4, 4), time = c(1, 2, 3), y = c(0.8, 1.4, 2.6), g = c(Inf, Inf, Inf))
  )
  Z <- cbind(1.0, c(-0.5, 0.2, -0.1, 0.7))
  po <- ns_orig$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl_o)
  pc <- ns_cand$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl_c)
  fo <- ns_orig$construct_ipw_frame_owner(po, id = 1:4, Z = Z, control = ctrl_o)
  fc <- ns_cand$construct_ipw_frame_owner(pc, id = 1:4, Z = Z, control = ctrl_c)
  nuo <- ns_orig$construct_ipw_nuisance_view(fo, g = 2, K = Inf, s = 1, t = 2, control = ctrl_o)
  nuc <- ns_cand$construct_ipw_nuisance_view(fc, g = 2, K = Inf, s = 1, t = 2, control = ctrl_c)
  lko <- ns_orig$construct_ipw_link_view(fo, s = 1, t = 3, nuisance_view = nuo, control = ctrl_o)
  lkc <- ns_cand$construct_ipw_link_view(fc, s = 1, t = 3, nuisance_view = nuc, control = ctrl_c)

  cho <- ns_orig$query_ipw_view_chunk(lko, 1L, 4L, ctrl_o, read_outcomes = TRUE, include_design = TRUE, enclose_differences = TRUE)
  chc <- ns_cand$query_ipw_view_chunk(lkc, 1L, 4L, ctrl_c, read_outcomes = TRUE, include_design = TRUE, enclose_differences = TRUE)
  expect_identical(noh_raw_bits(chc$observed_risk$y_left_obs), noh_raw_bits(cho$observed_risk$y_left_obs))
  expect_identical(noh_raw_bits(chc$observed_risk$y_right_obs), noh_raw_bits(cho$observed_risk$y_right_obs))
  expect_identical(noh_raw_bits(chc$observed_risk$Delta_iv_obs), noh_raw_bits(cho$observed_risk$Delta_iv_obs))
  expect_identical(chc$counts, cho$counts)
  expect_identical(chc$read_diagnostics, cho$read_diagnostics)
  expect_identical(chc$work, cho$work)
})
