# Frozen independent TEST reference (October 16, 2026): the historical pure-R
# directed interval matrix product `ipw_pipeline_interval_product` exactly as it
# stood in R/chained_ipw_pipeline.R before the native kernel, together with the
# scalar binary64 primitives it used (R/chained_gmm_core.R: gmm_next, gmm_iv,
# gmm_iv_check, gmm_iv_add, gmm_iv_neg, gmm_iv_sub, gmm_iv_mul, cdid_gmm_abort).
# It exists only so tests can compare the native kernel `C_gmm_iv_mat_product`
# bit for bit against the historical recurrence. It is never a production
# fallback and is not loaded by the package. Do not edit: the test pins its
# bytes.

cdid_gmm_abort <- function(message, class, diagnostics = list()) {
  stop(structure(
    class = unique(c(class, "cdid_gmm_refusal", "cdid_gmm_error", "cdid_error", "error", "condition")),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

# Scalar outward rounding. Assumptions: IEEE binary64, round-to-nearest ties-to-even,
# correctly rounded scalar +,-,*,/,sqrt, gradual underflow, no FTZ/DAZ/reassociation.
gmm_next <- function(x, up = TRUE) {
  if (length(x) != 1L || is.na(x)) cdid_gmm_abort("Invalid interval endpoint.", "cdid_gmm_accuracy_error")
  if (x == 0) return(if (up) 2^-1074 else -2^-1074)
  if (is.infinite(x)) {
    if ((x > 0 && up) || (x < 0 && !up)) return(x)
    return(if (x > 0) .Machine$double.xmax else -.Machine$double.xmax)
  }
  b <- as.integer(writeBin(as.double(x), raw(), size = 8L, endian = "little"))
  inc <- (x > 0) == up
  for (i in seq_len(8L)) {
    if (inc) {
      if (b[i] < 255L) { b[i] <- b[i] + 1L; break }
      b[i] <- 0L
    } else {
      if (b[i] > 0L) { b[i] <- b[i] - 1L; break }
      b[i] <- 255L
    }
  }
  readBin(as.raw(b), "double", n = 1L, size = 8L, endian = "little")
}
gmm_iv <- function(x) c(x, x)
gmm_iv_check <- function(z) {
  if (length(z) != 2L || any(!is.finite(z)) || z[1L] > z[2L])
    cdid_gmm_abort("Validated scalar arithmetic overflowed or became inconclusive.", "cdid_gmm_accuracy_error")
  z
}
gmm_iv_add <- function(a, b) {
  if (identical(unname(a), c(0, 0))) return(b)
  if (identical(unname(b), c(0, 0))) return(a)
  gmm_iv_check(c(gmm_next(a[1L] + b[1L], FALSE), gmm_next(a[2L] + b[2L])))
}
gmm_iv_neg <- function(a) -rev(a)
gmm_iv_sub <- function(a, b) gmm_iv_add(a, gmm_iv_neg(b))
gmm_iv_mul <- function(a, b) {
  if (all(a == 0) || all(b == 0)) return(c(0, 0))
  z <- c(a[1L]*b[1L], a[1L]*b[2L], a[2L]*b[1L], a[2L]*b[2L])
  gmm_iv_check(c(gmm_next(min(z), FALSE), gmm_next(max(z))))
}

ipw_pipeline_interval_product <- function(xlo, xhi, ylo, yhi = ylo) {
  outlo <- outhi <- matrix(0.0, nrow(xlo), ncol(ylo))
  for (i in seq_len(nrow(xlo))) for (j in seq_len(ncol(ylo))) {
    if (all(ylo[,j]==0.0 & yhi[,j]==0.0) || all(xlo[i,]==0.0 & xhi[i,]==0.0)) next
    z <- gmm_iv(0.0)
    for (k in seq_len(ncol(xlo))) {
      z <- gmm_iv_add(z, gmm_iv_mul(c(xlo[i,k], xhi[i,k]), c(ylo[k,j], yhi[k,j])))
    }
    outlo[i,j] <- z[1L]; outhi[i,j] <- z[2L]
  }
  list(lo = outlo, hi = outhi)
}
