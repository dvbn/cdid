# Independent bit-for-bit checks of the native directed interval matrix product
# (`C_gmm_iv_mat_product`, R wrapper `gmm_iv_mat_product`, pipeline entry
# `ipw_pipeline_interval_product`) against the frozen historical pure-R
# recurrence in reference/cp02_interval_product_R_reference.R. Everything here
# is deterministic: small explicit interval matrices, no RNG, no fit, no did
# namespace. In a plain helper environment without a loaded `cdid` DLL the
# copied C source is compiled once with R's own compiler settings and loaded as
# the actual registered DLL; an installed package run uses its namespace DLL.

ipn_raw_bits <- function(x) writeBin(as.double(x), raw(), size = 8L, endian = "little")

ipn_reference_env <- function() {
  path <- testthat::test_path("reference", "cp02_interval_product_R_reference.R")
  if (!file.exists(path)) stop("missing frozen R reference: ", path, call. = FALSE)
  # Frozen bytes: any edit of the independent reference is an error, never a skip.
  expect_identical(as.integer(file.info(path)$size), 3035L)
  expect_identical(unname(tools::md5sum(path)), "aebaf353099aa228f04783356fae6486")
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  env
}

ipn_locate_source <- function() {
  cands <- c(tryCatch(testthat::test_path("..", "..", "src", "cdid.c"), error = function(e) character(0L)),
             file.path("..", "..", "src", "cdid.c"), file.path("src", "cdid.c"))
  for (p in cands) if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = TRUE))
  NULL
}

# Compile a C source text into a loadable DLL named `name` using R CMD SHLIB
# (the existing compiler settings) and load it; returns the DLLInfo.
ipn_build_dll <- function(c_text, name) {
  dir <- tempfile(paste0("ipn_", name, "_")); dir.create(dir)
  src <- file.path(dir, paste0(name, ".c"))
  writeLines(c_text, src)
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  out <- system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", shQuote(basename(src))), stdout = TRUE, stderr = TRUE)
  so <- file.path(dir, paste0(name, .Platform$dynlib.ext))
  if (!is.null(attr(out, "status")) || !file.exists(so)) stop("R CMD SHLIB failed for ", name, ":\n", paste(out, collapse = "\n"), call. = FALSE)
  dyn.load(so)
}

ipn_production_syms <- function() {
  if (!("cdid" %in% names(base::getLoadedDLLs()))) {
    src <- ipn_locate_source()
    if (is.null(src)) stop("no loaded `cdid` DLL and no src/cdid.c to compile; the native kernel cannot be bound.", call. = FALSE)
    ipn_build_dll(readLines(src), "cdid")
  }
  gmm_require_native_dll()
}

# Deterministic interval-matrix cases (lo <= hi everywhere; all finite).
ipn_cases <- function() {
  xmax <- .Machine$double.xmax
  signs <- matrix(c(1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1), 6L, 5L)
  cl <- matrix(c(0.1, -2.5, 3, 2^-20, 7, -0.3, 1 / 3, -1 / 7, 11, 0, 2.75, -4, 0.5, 1e-3, -9), 5L, 3L)
  cancel_x <- matrix(c(2^53, 1, -2^53, -2^53, 1, 2^53), 2L, 3L, byrow = TRUE)
  cancel_y <- matrix(c(1, 1, 1, 0.5, -0.5, 0.25), 3L, 2L)
  bnd <- c(2^-1074, -2^-1074, 2^-1022, -(1 - 2^-53), 1 + 2^-52, 2^500, -2^500, 2^511, 0.1)
  bnd_x <- matrix(bnd, 3L, 3L)
  bnd_y <- matrix(c(1, 2^-1074, -1, 2^-600, 1 - 2^-53, 0, 2^-1022, -2^-1074, 3), 3L, 3L)
  zero_x <- rbind(c(0, 0, 0), c(1, 0, 2), c(0, 0, -0), c(-3, 4, 0))
  zero_y <- cbind(c(0, 0, 0), c(1, -1, 2), c(0, -0, 0), c(0.5, 0, 0))
  wide_x <- matrix(sin(1:12) * c(2^-40, 1, 2^40, 3, -2, 0.125, 2^-30, 7, -11, 0.3, 2^20, -2^-20), 4L, 3L)
  wide_y <- matrix(cos(1:6) * c(1, 2^-3, -5, 2^33, 0.7, -2^-33), 3L, 2L)
  rad <- 2^-40
  list(
    signs_by_cluster = list(xlo = signs, xhi = signs, ylo = cl, yhi = cl),
    cancellation = list(xlo = cancel_x, xhi = cancel_x, ylo = cancel_y, yhi = cancel_y),
    subnormal_boundary = list(xlo = bnd_x, xhi = bnd_x, ylo = bnd_y, yhi = bnd_y),
    zero_shortcuts_and_signed_zero = list(xlo = zero_x, xhi = zero_x, ylo = zero_y, yhi = zero_y),
    nondegenerate_boxes = list(xlo = wide_x - rad * abs(wide_x), xhi = wide_x + rad * abs(wide_x),
                               ylo = wide_y - rad, yhi = wide_y + rad),
    mixed_degenerate_rows = list(xlo = rbind(signs[1:2, ], 0), xhi = rbind(signs[1:2, ], 0),
                                 ylo = cbind(cl[, 1], 0, cl[, 3]), yhi = cbind(cl[, 1] + 2^-45, 0, cl[, 3] + 2^-45)),
    single_entry = list(xlo = matrix(-2^-1074), xhi = matrix(2^-1074), ylo = matrix(xmax / 4), yhi = matrix(xmax / 4))
  )
}

test_that("native interval product matches the frozen R recurrence bit for bit on deterministic cases", {
  ref <- ipn_reference_env()
  syms <- ipn_production_syms()
  cases <- ipn_cases()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    r <- ref$ipw_pipeline_interval_product(cs$xlo, cs$xhi, cs$ylo, cs$yhi)
    n_wrap <- gmm_iv_mat_product(cs$xlo, cs$xhi, cs$ylo, cs$yhi)
    n_call <- .Call(syms$C_gmm_iv_mat_product, cs$xlo, cs$xhi, cs$ylo, cs$yhi, NULL)
    n_pipe <- ipw_pipeline_interval_product(cs$xlo, cs$xhi, cs$ylo, cs$yhi)
    expect_identical(names(n_wrap), c("lo", "hi"), info = nm)
    expect_identical(dim(n_wrap$lo), c(nrow(cs$xlo), ncol(cs$ylo)), info = nm)
    expect_identical(dim(n_wrap$hi), dim(r$hi), info = nm)
    expect_identical(ipn_raw_bits(n_wrap$lo), ipn_raw_bits(r$lo), info = nm)
    expect_identical(ipn_raw_bits(n_wrap$hi), ipn_raw_bits(r$hi), info = nm)
    expect_identical(ipn_raw_bits(n_call$lo), ipn_raw_bits(r$lo), info = nm)
    expect_identical(ipn_raw_bits(n_call$hi), ipn_raw_bits(r$hi), info = nm)
    expect_identical(ipn_raw_bits(n_pipe$lo), ipn_raw_bits(r$lo), info = nm)
    expect_identical(ipn_raw_bits(n_pipe$hi), ipn_raw_bits(r$hi), info = nm)
    expect_true(all(n_wrap$lo <= n_wrap$hi), info = nm)
    # Enclosure of the exact product is implied by outward rounding; checked on
    # the exactly representable sign-by-cluster case at double precision.
    if (identical(nm, "signs_by_cluster")) {
      exact <- cs$xlo %*% cs$ylo
      expect_true(all(n_wrap$lo <= exact + 1e-12 & n_wrap$hi >= exact - 1e-12))
      expect_true(any(n_wrap$lo < n_wrap$hi))
    }
    if (identical(nm, "zero_shortcuts_and_signed_zero")) {
      # Shortcut rows/columns stay +0.0 in both endpoints (R reference law).
      expect_identical(ipn_raw_bits(n_wrap$lo[1L, ]), ipn_raw_bits(rep(0.0, 4L)))
      expect_identical(ipn_raw_bits(n_wrap$hi[3L, ]), ipn_raw_bits(rep(0.0, 4L)))
      expect_identical(ipn_raw_bits(n_wrap$lo[, 1L]), ipn_raw_bits(rep(0.0, 4L)))
      expect_identical(ipn_raw_bits(n_wrap$hi[, 3L]), ipn_raw_bits(rep(0.0, 4L)))
    }
  }
  # Production default `yhi = ylo` (degenerate second operand) through the pipeline entry.
  cs <- cases$signs_by_cluster
  r1 <- ref$ipw_pipeline_interval_product(cs$xlo, cs$xhi, cs$ylo)
  n1 <- ipw_pipeline_interval_product(cs$xlo, cs$xhi, cs$ylo)
  expect_identical(ipn_raw_bits(n1$lo), ipn_raw_bits(r1$lo)); expect_identical(ipn_raw_bits(n1$hi), ipn_raw_bits(r1$hi))
  # Integer +/-1 sign matrices are coerced exactly (same bits as the double matrix).
  si <- cs$xlo; storage.mode(si) <- "integer"
  n2 <- gmm_iv_mat_product(si, si, cs$ylo, cs$yhi)
  expect_identical(ipn_raw_bits(n2$lo), ipn_raw_bits(r1$lo)); expect_identical(ipn_raw_bits(n2$hi), ipn_raw_bits(r1$hi))
  # Transposed-operand pattern used by the covariance certificate: t(C) %*% C.
  C <- cases$nondegenerate_boxes
  rt <- ref$ipw_pipeline_interval_product(t(C$ylo), t(C$yhi), C$ylo, C$yhi)
  nt <- ipw_pipeline_interval_product(t(C$ylo), t(C$yhi), C$ylo, C$yhi)
  expect_identical(ipn_raw_bits(nt$lo), ipn_raw_bits(rt$lo)); expect_identical(ipn_raw_bits(nt$hi), ipn_raw_bits(rt$hi))
})

test_that("native interval product keeps the typed failures of the recurrence and the native entry contract", {
  ref <- ipn_reference_env()
  syms <- ipn_production_syms()
  xmax <- .Machine$double.xmax
  H <- matrix(xmax, 2L, 2L)
  # Overflow inside the recurrence: same class and message as the frozen reference.
  e_ref <- tryCatch(ref$ipw_pipeline_interval_product(H, H, H, H), error = identity)
  e_nat <- tryCatch(gmm_iv_mat_product(H, H, H, H), error = identity)
  e_call <- tryCatch(.Call(syms$C_gmm_iv_mat_product, H, H, H, H, NULL), error = identity)
  expect_s3_class(e_ref, "cdid_gmm_accuracy_error"); expect_s3_class(e_nat, "cdid_gmm_accuracy_error"); expect_s3_class(e_call, "cdid_gmm_accuracy_error")
  expect_identical(conditionMessage(e_nat), conditionMessage(e_ref))
  expect_identical(conditionMessage(e_call), conditionMessage(e_ref))
  # A single overflowing product (2^600 * 2^600) with otherwise harmless entries.
  O <- matrix(c(2^600, 1, 1, 1), 2L, 2L)
  e_ref2 <- tryCatch(ref$ipw_pipeline_interval_product(O, O, O, O), error = identity)
  e_nat2 <- tryCatch(ipw_pipeline_interval_product(O, O, O, O), error = identity)
  expect_s3_class(e_ref2, "cdid_gmm_accuracy_error"); expect_s3_class(e_nat2, "cdid_gmm_accuracy_error")
  expect_identical(conditionMessage(e_nat2), conditionMessage(e_ref2))
  # Nonfinite operands: refused as accuracy errors by both (the reference through its interval check).
  I <- matrix(c(Inf, 1, 1, 1), 2L, 2L); U <- matrix(1, 2L, 2L)
  expect_error(ref$ipw_pipeline_interval_product(I, I, U, U), class = "cdid_gmm_accuracy_error")
  expect_error(gmm_iv_mat_product(I, I, U, U), class = "cdid_gmm_accuracy_error")
  expect_error(.Call(syms$C_gmm_iv_mat_product, U, U, I, I, NULL), class = "cdid_gmm_accuracy_error")
  N <- matrix(c(NA_real_, 1, 1, 1), 2L, 2L)
  expect_error(gmm_iv_mat_product(N, N, U, U), class = "cdid_gmm_accuracy_error")
  expect_error(.Call(syms$C_gmm_iv_mat_product, N, N, U, U, NULL), class = "cdid_gmm_accuracy_error")
  # Invalid operands under the established native contract (before any allocation).
  M23 <- matrix(1, 2L, 3L); M32 <- matrix(1, 3L, 2L); M22 <- matrix(1, 2L, 2L)
  expect_error(gmm_iv_mat_product(M22, M22, M32, M32), class = "cdid_gmm_bound_error")       # ncol(x) != nrow(y)
  expect_error(.Call(syms$C_gmm_iv_mat_product, M22, M22, M32, M32, NULL), class = "cdid_gmm_bound_error")
  expect_error(gmm_iv_mat_product(M23, M22, M32, M32), class = "cdid_gmm_bound_error")       # dim(xlo) != dim(xhi)
  expect_error(.Call(syms$C_gmm_iv_mat_product, M23, M23, M32, M23, NULL), class = "cdid_gmm_bound_error")
  expect_error(gmm_iv_mat_product(c(1, 2), c(1, 2), M22, M22), class = "cdid_gmm_accuracy_error")  # not a matrix
  expect_error(.Call(syms$C_gmm_iv_mat_product, c(1, 2), c(1, 2), M22, M22, NULL), class = "cdid_gmm_accuracy_error")
  Mi <- matrix(1L, 2L, 2L)
  expect_error(.Call(syms$C_gmm_iv_mat_product, Mi, Mi, M22, M22, NULL), class = "cdid_gmm_accuracy_error")  # integer storage at the C boundary
  expect_error(gmm_iv_mat_product(structure(M22, class = "custom_obj"), M22, M22, M22), class = "cdid_gmm_accuracy_error")
  expect_error(gmm_iv_mat_product(matrix(1, 0L, 2L), matrix(1, 0L, 2L), M22, M22), class = "cdid_gmm_bound_error")
  inv <- M22; inv[1L, 1L] <- 2
  expect_error(gmm_iv_mat_product(inv, M22, M22, M22), class = "cdid_gmm_accuracy_error")      # xlo > xhi
  expect_error(.Call(syms$C_gmm_iv_mat_product, M22, M22, inv, M22, NULL), class = "cdid_gmm_accuracy_error")  # ylo > yhi
  # Allocation controls: the interval output (lo + hi) and every operand are
  # charged against an explicit `max_matrix_entries` before allocation.
  expect_error(gmm_iv_mat_product(M23, M23, M32, M32, max_matrix_entries = 7), class = "cdid_gmm_bound_error")   # 2 * (2 x 2) = 8 > 7
  expect_error(.Call(syms$C_gmm_iv_mat_product, M23, M23, M32, M32, 7), class = "cdid_gmm_bound_error")
  expect_error(gmm_iv_mat_product(M23, M23, M32, M32, max_matrix_entries = 5), class = "cdid_gmm_bound_error")   # operand 6 > 5
  ok <- gmm_iv_mat_product(M23, M23, M32, M32, max_matrix_entries = 8)
  expect_identical(ipn_raw_bits(ok$lo), ipn_raw_bits(ref$ipw_pipeline_interval_product(M23, M23, M32, M32)$lo))
  expect_error(gmm_iv_mat_product(M23, M23, M32, M32, max_matrix_entries = 0), class = "cdid_gmm_bound_error")
  expect_error(.Call(syms$C_gmm_iv_mat_product, M23, M23, M32, M32, -1), class = "cdid_gmm_bound_error")
  # Registration of the new entry (narrow): five arguments, listed hotspot.
  expect_identical(as.integer(getNativeSymbolInfo("C_gmm_iv_mat_product", PACKAGE = base::getLoadedDLLs()[["cdid"]])$numParameters), 5L)
  expect_true("C_gmm_iv_mat_product" %in% .Call(syms$C_gmm_native_status)$compiled_hotspots)
})

test_that("the independent checks discriminate a kernel error: deliberate accumulate-without-rounding mutation", {
  ref <- ipn_reference_env()
  cases <- ipn_cases()
  # (1) R-level mutant of the recurrence: the per-term product is unchanged but
  # the accumulation drops outward rounding. Same inputs, same bit comparison.
  mutant_r <- function(xlo, xhi, ylo, yhi) {
    outlo <- outhi <- matrix(0.0, nrow(xlo), ncol(ylo))
    for (i in seq_len(nrow(xlo))) for (j in seq_len(ncol(ylo))) {
      if (all(ylo[, j] == 0.0 & yhi[, j] == 0.0) || all(xlo[i, ] == 0.0 & xhi[i, ] == 0.0)) next
      z <- c(0, 0)
      for (k in seq_len(ncol(xlo))) {
        p <- ref$gmm_iv_mul(c(xlo[i, k], xhi[i, k]), c(ylo[k, j], yhi[k, j]))
        z <- c(z[1L] + p[1L], z[2L] + p[2L])
      }
      outlo[i, j] <- z[1L]; outhi[i, j] <- z[2L]
    }
    list(lo = outlo, hi = outhi)
  }
  discriminated <- 0L
  for (nm in c("signs_by_cluster", "cancellation", "nondegenerate_boxes")) {
    cs <- cases[[nm]]
    r <- ref$ipw_pipeline_interval_product(cs$xlo, cs$xhi, cs$ylo, cs$yhi)
    m <- mutant_r(cs$xlo, cs$xhi, cs$ylo, cs$yhi)
    same <- identical(ipn_raw_bits(m$lo), ipn_raw_bits(r$lo)) && identical(ipn_raw_bits(m$hi), ipn_raw_bits(r$hi))
    if (!same) discriminated <- discriminated + 1L
  }
  expect_identical(discriminated, 3L)
  # The mutant still agrees where no accumulation occurs (shortcut-only output),
  # so the discrimination above is attributable to the accumulation step.
  z <- cases$zero_shortcuts_and_signed_zero
  z1 <- list(xlo = z$xlo[c(1L, 3L), , drop = FALSE], xhi = z$xhi[c(1L, 3L), , drop = FALSE], ylo = z$ylo, yhi = z$yhi)
  expect_identical(ipn_raw_bits(mutant_r(z1$xlo, z1$xhi, z1$ylo, z1$yhi)$hi),
                   ipn_raw_bits(ref$ipw_pipeline_interval_product(z1$xlo, z1$xhi, z1$ylo, z1$yhi)$hi))
})

test_that("a compiled C mutant of the kernel's accumulation step is discriminated by the same bit comparison", {
  src <- ipn_locate_source()
  if (is.null(src)) skip("src/cdid.c is not reachable from this test directory (installed-package run); the C-level mutant compile needs the source tree. The R-level mutation check above ran.")
  ref <- ipn_reference_env()
  txt <- readLines(src)
  marker <- "z = gmm_iv_add_c(z, prod, &st); /* IPN_ACCUMULATE */"
  init <- "void R_init_cdid(DllInfo *dll)"
  expect_identical(sum(grepl(marker, txt, fixed = TRUE)), 1L)
  expect_identical(sum(grepl(init, txt, fixed = TRUE)), 2L)   # prototype + definition
  mut <- sub(marker, "z.lo = gmm_fp_add(z.lo, prod.lo); z.hi = gmm_fp_add(z.hi, prod.hi); /* IPN_MUTATED */", txt, fixed = TRUE)
  mut <- gsub(init, "void R_init_cdid_ipn_mutant(DllInfo *dll)", mut, fixed = TRUE)
  expect_identical(sum(grepl("IPN_MUTATED", mut, fixed = TRUE)), 1L)
  dll_m <- ipn_build_dll(mut, "cdid_ipn_mutant")
  on.exit(try(dyn.unload(dll_m[["path"]]), silent = TRUE), add = TRUE)
  sym_m <- getNativeSymbolInfo("C_gmm_iv_mat_product", PACKAGE = dll_m)$address
  cases <- ipn_cases()
  discriminated <- 0L
  for (nm in c("signs_by_cluster", "cancellation", "nondegenerate_boxes")) {
    cs <- cases[[nm]]
    r <- ref$ipw_pipeline_interval_product(cs$xlo, cs$xhi, cs$ylo, cs$yhi)
    m <- .Call(sym_m, cs$xlo, cs$xhi, cs$ylo, cs$yhi, NULL)
    same <- identical(ipn_raw_bits(m$lo), ipn_raw_bits(r$lo)) && identical(ipn_raw_bits(m$hi), ipn_raw_bits(r$hi))
    if (!same) discriminated <- discriminated + 1L
  }
  expect_identical(discriminated, 3L)
  # Shortcut-only inputs still agree: the compiled mutant differs only in accumulation.
  z <- cases$zero_shortcuts_and_signed_zero
  z1 <- list(xlo = z$xlo[c(1L, 3L), , drop = FALSE], xhi = z$xhi[c(1L, 3L), , drop = FALSE], ylo = z$ylo, yhi = z$yhi)
  m0 <- .Call(sym_m, z1$xlo, z1$xhi, z1$ylo, z1$yhi, NULL)
  expect_identical(ipn_raw_bits(m0$lo), ipn_raw_bits(ref$ipw_pipeline_interval_product(z1$xlo, z1$xhi, z1$ylo, z1$yhi)$lo))
  # The production kernel (unmutated) still passes the same comparison afterwards.
  syms <- ipn_production_syms()
  cs <- cases$signs_by_cluster
  r <- ref$ipw_pipeline_interval_product(cs$xlo, cs$xhi, cs$ylo, cs$yhi)
  p <- .Call(syms$C_gmm_iv_mat_product, cs$xlo, cs$xhi, cs$ylo, cs$yhi, NULL)
  expect_identical(ipn_raw_bits(p$lo), ipn_raw_bits(r$lo)); expect_identical(ipn_raw_bits(p$hi), ipn_raw_bits(r$hi))
})
