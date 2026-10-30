#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("fp-contract=off")
#endif
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#pragma STDC FP_CONTRACT OFF
#endif

#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <fenv.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)) || \
    (defined(__ASSOC_MATH__) && (__ASSOC_MATH__ > 0)) || \
    (defined(__RECIPROCAL_MATH__) && (__RECIPROCAL_MATH__ > 0))
#error "Fast-math / associative-math / finite-math-only compilation is forbidden for validated interval arithmetic."
#endif

#if !defined(FLT_RADIX) || (FLT_RADIX != 2) || \
    !defined(DBL_MANT_DIG) || (DBL_MANT_DIG != 53) || \
    !defined(DBL_MAX_EXP) || (DBL_MAX_EXP != 1024) || \
    !defined(DBL_MIN_EXP) || (DBL_MIN_EXP != -1021)
#error "Unsupported floating-point model: IEEE-754 binary64 (FLT_RADIX=2, DBL_MANT_DIG=53) is required."
#endif

typedef char gmm_static_assert_double_64bit[(sizeof(double) == 8 && sizeof(uint64_t) == 8) ? 1 : -1];

#define GMM_MAX_EXACT_INT 9007199254740991ULL /* 2^53 - 1 */

typedef struct {
  double lo;
  double hi;
} gmm_iv_t;

typedef enum {
  GMM_IV_OK = 0,
  GMM_IV_ERR_OVERFLOW = 1,
  GMM_IV_ERR_DIV_ZERO = 2,
  GMM_IV_ERR_NEG_SQRT = 3,
  GMM_IV_ERR_ENDPOINT = 4
} gmm_iv_status_t;

static void gmm_c_abort(int nprotect, const char *message, const char *cls) {
  if (nprotect > 0) {
    UNPROTECT(nprotect);
  }
  SEXP cond = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("message"));
  SET_STRING_ELT(names, 1, Rf_mkChar("call"));
  SET_STRING_ELT(names, 2, Rf_mkChar("diagnostics"));
  Rf_setAttrib(cond, R_NamesSymbol, names);
  SET_VECTOR_ELT(cond, 0, Rf_mkString(message));
  SET_VECTOR_ELT(cond, 1, R_NilValue);
  SET_VECTOR_ELT(cond, 2, Rf_allocVector(VECSXP, 0));

  SEXP classes;
  if (strcmp(cls, "cdid_gmm_native_error") == 0) {
    classes = PROTECT(Rf_allocVector(STRSXP, 7));
    SET_STRING_ELT(classes, 0, Rf_mkChar("cdid_gmm_native_error"));
    SET_STRING_ELT(classes, 1, Rf_mkChar("cdid_gmm_accuracy_error"));
    SET_STRING_ELT(classes, 2, Rf_mkChar("cdid_gmm_refusal"));
    SET_STRING_ELT(classes, 3, Rf_mkChar("cdid_gmm_error"));
    SET_STRING_ELT(classes, 4, Rf_mkChar("cdid_error"));
    SET_STRING_ELT(classes, 5, Rf_mkChar("error"));
    SET_STRING_ELT(classes, 6, Rf_mkChar("condition"));
  } else {
    classes = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_STRING_ELT(classes, 0, Rf_mkChar(cls));
    SET_STRING_ELT(classes, 1, Rf_mkChar("cdid_gmm_refusal"));
    SET_STRING_ELT(classes, 2, Rf_mkChar("cdid_gmm_error"));
    SET_STRING_ELT(classes, 3, Rf_mkChar("cdid_error"));
    SET_STRING_ELT(classes, 4, Rf_mkChar("error"));
    SET_STRING_ELT(classes, 5, Rf_mkChar("condition"));
  }
  Rf_setAttrib(cond, R_ClassSymbol, classes);

  SEXP stop_call = PROTECT(Rf_lang2(Rf_install("stop"), cond));
  Rf_eval(stop_call, R_BaseEnv);
  UNPROTECT(4);
  Rf_error("%s", message);
}

static void gmm_c_abort_status(int nprotect, gmm_iv_status_t st) {
  switch (st) {
    case GMM_IV_ERR_DIV_ZERO:
      gmm_c_abort(nprotect, "Interval division crosses zero.", "cdid_gmm_accuracy_error");
      break;
    case GMM_IV_ERR_NEG_SQRT:
      gmm_c_abort(nprotect, "Negative interval square root.", "cdid_gmm_accuracy_error");
      break;
    case GMM_IV_ERR_ENDPOINT:
      gmm_c_abort(nprotect, "Invalid interval endpoint.", "cdid_gmm_accuracy_error");
      break;
    case GMM_IV_ERR_OVERFLOW:
    default:
      gmm_c_abort(nprotect, "Validated scalar arithmetic overflowed or became inconclusive.", "cdid_gmm_accuracy_error");
      break;
  }
}

/* Explicit volatile stores prevent FMA contraction, reassociation, and excess precision. */
static inline double gmm_fp_add(double x, double y) {
  volatile double vx = x;
  volatile double vy = y;
  volatile double vz = vx + vy;
  return vz;
}

static inline double gmm_fp_sub(double x, double y) {
  volatile double vx = x;
  volatile double vy = y;
  volatile double vz = vx - vy;
  return vz;
}

static inline double gmm_fp_mul(double x, double y) {
  volatile double vx = x;
  volatile double vy = y;
  volatile double vz = vx * vy;
  return vz;
}

static inline double gmm_fp_div(double x, double y) {
  volatile double vx = x;
  volatile double vy = y;
  volatile double vz = vx / vy;
  return vz;
}

static inline double gmm_fp_sqrt(double x) {
  volatile double vx = x;
  volatile double vz = sqrt(vx);
  return vz;
}

static inline uint64_t gmm_double_to_bits(double x) {
  volatile double vx = x;
  double tmp = vx;
  uint64_t u = 0;
  memcpy(&u, &tmp, sizeof(uint64_t));
  return u;
}

static inline double gmm_next_c(double x, int up) {
  volatile double vx = x;
  if (isnan(vx)) {
    return R_NaN;
  }
  if (vx == 0.0) {
    volatile double z = up ? nextafter(0.0, INFINITY) : nextafter(0.0, -INFINITY);
    return z;
  }
  if (isinf(vx)) {
    if ((vx > 0.0 && up) || (vx < 0.0 && !up)) {
      return vx;
    }
    volatile double z = (vx > 0.0) ? DBL_MAX : -DBL_MAX;
    return z;
  }
  volatile double z = nextafter(vx, up ? INFINITY : -INFINITY);
  return z;
}

static int gmm_verify_arithmetic_values(
    int rmode,
    double min_sub_pos,
    double even_tie_pos,
    double odd_tie_pos,
    double fma_residual,
    const char **err_msg) {
#if defined(FE_TONEAREST)
  if (rmode != FE_TONEAREST) {
    *err_msg = "Unsupported floating-point rounding mode: FE_TONEAREST is required.";
    return 0;
  }
#endif
  if (gmm_double_to_bits(1.0) != 0x3ff0000000000000ULL ||
      gmm_double_to_bits(-2.0) != 0xc000000000000000ULL) {
    *err_msg = "Unsupported IEEE-754 binary64 bit representation.";
    return 0;
  }
  if (!(min_sub_pos > 0.0) ||
      gmm_double_to_bits(min_sub_pos) != 0x0000000000000001ULL ||
      gmm_fp_add(min_sub_pos, 0.0) != min_sub_pos ||
      gmm_fp_mul(min_sub_pos, 1.0) != min_sub_pos ||
      gmm_fp_div(min_sub_pos, 2.0) != 0.0 ||
      gmm_next_c(min_sub_pos, 0) != 0.0) {
    *err_msg = "Gradual underflow check failed (FTZ/DAZ or subnormal flushing detected).";
    return 0;
  }
  double min_sub_neg = gmm_next_c(0.0, 0);
  if (!(min_sub_neg < 0.0) ||
      gmm_double_to_bits(min_sub_neg) != 0x8000000000000001ULL ||
      gmm_fp_add(min_sub_pos, min_sub_neg) != 0.0) {
    *err_msg = "Negative subnormal stepping or cancellation check failed.";
    return 0;
  }
  double half_min_norm = gmm_fp_mul(DBL_MIN, 0.5);
  double max_sub = gmm_fp_sub(DBL_MIN, min_sub_pos);
  if (!(half_min_norm > 0.0 && half_min_norm < DBL_MIN) ||
      !(max_sub > 0.0 && max_sub < DBL_MIN) ||
      gmm_fp_add(max_sub, min_sub_pos) != DBL_MIN) {
    *err_msg = "Normal-to-subnormal gradual underflow transition check failed.";
    return 0;
  }
  double u = 0x1.0p-53;
  double ulp1 = 0x1.0p-52;
  double tiny = 0x1.0p-105;
  if (even_tie_pos != 1.0 ||
      odd_tie_pos != gmm_fp_add(1.0, 2.0 * ulp1) ||
      gmm_fp_sub(-1.0, u) != -1.0 ||
      gmm_fp_sub(-gmm_fp_add(1.0, ulp1), u) != -gmm_fp_add(1.0, 2.0 * ulp1) ||
      gmm_fp_add(1.0, gmm_fp_add(u, tiny)) != gmm_fp_add(1.0, ulp1) ||
      gmm_fp_sub(-1.0, gmm_fp_add(u, tiny)) != -gmm_fp_add(1.0, ulp1)) {
    *err_msg = "Round-to-nearest ties-to-even (RNTE) addition/subtraction probe failed.";
    return 0;
  }
  if (fma_residual != 0.0) {
    *err_msg = "Unwanted FMA contraction or non-RNTE multiplication detected.";
    return 0;
  }
  if (gmm_fp_div(1.0, 3.0) != 0x1.5555555555555p-2 ||
      gmm_fp_sqrt(2.0) != 0x1.6a09e667f3bcdp+0) {
    *err_msg = "Round-to-nearest ties-to-even division/sqrt probe failed.";
    return 0;
  }
  return 1;
}

static void gmm_require_valid_arithmetic(int nprotect) {
  int rmode = 0;
#if defined(FE_TONEAREST)
  rmode = fegetround();
#endif
  double u = 0x1.0p-53;
  double ulp1 = 0x1.0p-52;
  double min_sub_pos = gmm_next_c(0.0, 1);
  double even_tie_pos = gmm_fp_add(1.0, u);
  double odd_tie_pos = gmm_fp_add(gmm_fp_add(1.0, ulp1), u);
  double a_fma = gmm_fp_add(1.0, ulp1);
  double prod = gmm_fp_mul(a_fma, a_fma);
  double fma_res = gmm_fp_sub(prod, gmm_fp_add(1.0, 2.0 * ulp1));
  const char *err_msg = NULL;
  if (!gmm_verify_arithmetic_values(rmode, min_sub_pos, even_tie_pos, odd_tie_pos, fma_res, &err_msg)) {
    gmm_c_abort(nprotect, err_msg, "cdid_gmm_native_error");
  }
}

static inline int gmm_buffers_overlap(const double *a, R_xlen_t na, const double *b, R_xlen_t nb) {
  if (a == NULL || b == NULL || na <= 0 || nb <= 0) return 0;
  if (a == b) return 1;
  uintptr_t ua = (uintptr_t) a;
  uintptr_t ub = (uintptr_t) b;
  uintptr_t bytes_a = (uintptr_t) na * sizeof(double);
  uintptr_t bytes_b = (uintptr_t) nb * sizeof(double);
  uintptr_t ea = (ua <= UINTPTR_MAX - bytes_a) ? (ua + bytes_a) : UINTPTR_MAX;
  uintptr_t eb = (ub <= UINTPTR_MAX - bytes_b) ? (ub + bytes_b) : UINTPTR_MAX;
  return (ua < eb) && (ub < ea);
}

static uint64_t gmm_checked_mul_u64(int nprotect, uint64_t a, uint64_t b, double max_entries, const char *label) {
  if (a > 0ULL && b > GMM_MAX_EXACT_INT / a) {
    gmm_c_abort(nprotect, "Entry count overflows exact IEEE binary64 integer range.", "cdid_gmm_bound_error");
  }
  uint64_t prod = a * b;
  if (prod > GMM_MAX_EXACT_INT) {
    gmm_c_abort(nprotect, "Allocation product exceeds exact integer range.", "cdid_gmm_bound_error");
  }
  if (R_finite(max_entries) && (double) prod > max_entries) {
    char buf[256];
    snprintf(buf, sizeof(buf),
             "Allocation for %s (%.0f entries) exceeds explicit `max_matrix_entries` (%.0f).",
             label, (double) prod, max_entries);
    gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
  }
  return prod;
}

static double gmm_parse_max_entries(int nprotect, SEXP max_entries_sexp) {
  if (max_entries_sexp == R_NilValue) {
    return R_PosInf;
  }
  if (Rf_isObject(max_entries_sexp) ||
      !(TYPEOF(max_entries_sexp) == REALSXP || TYPEOF(max_entries_sexp) == INTSXP) ||
      Rf_getAttrib(max_entries_sexp, R_DimSymbol) != R_NilValue ||
      XLENGTH(max_entries_sexp) != 1) {
    gmm_c_abort(nprotect, "`max_matrix_entries` must be an explicit finite integer bound in [1, 2^53 - 1].",
                "cdid_gmm_bound_error");
  }
  double val = (TYPEOF(max_entries_sexp) == REALSXP)
                   ? REAL(max_entries_sexp)[0]
                   : ((INTEGER(max_entries_sexp)[0] == NA_INTEGER) ? R_NaReal : (double) INTEGER(max_entries_sexp)[0]);
  if (!R_finite(val) || val != floor(val) || val < 1.0 || val > (double) GMM_MAX_EXACT_INT) {
    gmm_c_abort(nprotect, "`max_matrix_entries` must be an explicit finite integer bound in [1, 2^53 - 1].",
                "cdid_gmm_bound_error");
  }
  return val;
}

static void gmm_validate_real_vector_header(
    int nprotect, SEXP x, const char *name, double max_entries, R_xlen_t *out_len, const double **out_ptr) {
  if (x == R_NilValue || Rf_isObject(x) || TYPEOF(x) != REALSXP ||
      Rf_getAttrib(x, R_DimSymbol) != R_NilValue) {
    char buf[160];
    snprintf(buf, sizeof(buf), "`%s` must be a non-object atomic double vector.", name);
    gmm_c_abort(nprotect, buf, "cdid_gmm_accuracy_error");
  }
  R_xlen_t n = XLENGTH(x);
  if (n < 1) {
    gmm_c_abort(nprotect, "Invalid empty vector length in native interval guard.",
                "cdid_gmm_bound_error");
  }
  (void) gmm_checked_mul_u64(nprotect, (uint64_t) n, 1ULL, max_entries, name);
  *out_len = n;
  *out_ptr = NULL;
}

static void gmm_validate_real_matrix_header(
    int nprotect, SEXP x, const char *name, double max_entries,
    int *out_nr, int *out_nc, R_xlen_t *out_len, const double **out_ptr) {
  if (x == R_NilValue || Rf_isObject(x) || TYPEOF(x) != REALSXP) {
    char buf[160];
    snprintf(buf, sizeof(buf), "`%s` must be a non-object double matrix.", name);
    gmm_c_abort(nprotect, buf, "cdid_gmm_accuracy_error");
  }
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (dim == R_NilValue || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) {
    char buf[160];
    snprintf(buf, sizeof(buf), "`%s` must have a 2D integer `dim` attribute.", name);
    gmm_c_abort(nprotect, buf, "cdid_gmm_accuracy_error");
  }
  int nr = INTEGER(dim)[0];
  int nc = INTEGER(dim)[1];
  if (nr == NA_INTEGER || nc == NA_INTEGER || nr < 1 || nc < 1) {
    char buf[160];
    snprintf(buf, sizeof(buf), "`%s` must have positive finite matrix dimensions.", name);
    gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
  }
  uint64_t prod = gmm_checked_mul_u64(nprotect, (uint64_t) nr, (uint64_t) nc, max_entries, name);
  if ((R_xlen_t) prod != XLENGTH(x)) {
    gmm_c_abort(nprotect, "Matrix dimension product disagrees with storage length.",
                "cdid_gmm_bound_error");
  }
  *out_nr = nr;
  *out_nc = nc;
  *out_len = (R_xlen_t) prod;
  *out_ptr = NULL;
}

static void gmm_validate_real_finite(int nprotect, const double *px, R_xlen_t n, const char *name) {
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    if (!R_finite(px[i])) {
      char buf[160];
      snprintf(buf, sizeof(buf), "Nonfinite element in `%s`.", name);
      gmm_c_abort(nprotect, buf, "cdid_gmm_accuracy_error");
    }
  }
}

static inline gmm_iv_t gmm_iv_check_c(double lo, double hi, gmm_iv_status_t *st) {
  gmm_iv_t out = {lo, hi};
  if (*st != GMM_IV_OK) return out;
  if (!R_finite(lo) || !R_finite(hi) || lo > hi) {
    *st = GMM_IV_ERR_OVERFLOW;
  }
  return out;
}

static inline gmm_iv_t gmm_iv_add_c(gmm_iv_t a, gmm_iv_t b, gmm_iv_status_t *st) {
  if (*st != GMM_IV_OK) return a;
  if (a.lo == 0.0 && a.hi == 0.0) return b;
  if (b.lo == 0.0 && b.hi == 0.0) return a;
  double s_lo = gmm_fp_add(a.lo, b.lo);
  double s_hi = gmm_fp_add(a.hi, b.hi);
  if (isnan(s_lo) || isnan(s_hi)) {
    *st = GMM_IV_ERR_ENDPOINT;
    return a;
  }
  double lo = gmm_next_c(s_lo, 0);
  double hi = gmm_next_c(s_hi, 1);
  return gmm_iv_check_c(lo, hi, st);
}

static inline gmm_iv_t gmm_iv_neg_c(gmm_iv_t a) {
  volatile double n_lo = -a.hi;
  volatile double n_hi = -a.lo;
  gmm_iv_t out = {n_lo, n_hi};
  return out;
}

static inline gmm_iv_t gmm_iv_sub_c(gmm_iv_t a, gmm_iv_t b, gmm_iv_status_t *st) {
  return gmm_iv_add_c(a, gmm_iv_neg_c(b), st);
}

static inline gmm_iv_t gmm_iv_mul_c(gmm_iv_t a, gmm_iv_t b, gmm_iv_status_t *st) {
  gmm_iv_t zero_iv = {0.0, 0.0};
  if (*st != GMM_IV_OK) return zero_iv;
  if ((a.lo == 0.0 && a.hi == 0.0) || (b.lo == 0.0 && b.hi == 0.0)) {
    return zero_iv;
  }
  double z0 = gmm_fp_mul(a.lo, b.lo);
  double z1 = gmm_fp_mul(a.lo, b.hi);
  double z2 = gmm_fp_mul(a.hi, b.lo);
  double z3 = gmm_fp_mul(a.hi, b.hi);
  if (isnan(z0) || isnan(z1) || isnan(z2) || isnan(z3)) {
    *st = GMM_IV_ERR_ENDPOINT;
    return zero_iv;
  }
  double zmin = z0;
  if (z1 < zmin) zmin = z1;
  if (z2 < zmin) zmin = z2;
  if (z3 < zmin) zmin = z3;
  double zmax = z0;
  if (z1 > zmax) zmax = z1;
  if (z2 > zmax) zmax = z2;
  if (z3 > zmax) zmax = z3;
  double lo = gmm_next_c(zmin, 0);
  double hi = gmm_next_c(zmax, 1);
  return gmm_iv_check_c(lo, hi, st);
}

static inline gmm_iv_t gmm_iv_div_c(gmm_iv_t a, gmm_iv_t b, gmm_iv_status_t *st) {
  gmm_iv_t zero_iv = {0.0, 0.0};
  if (*st != GMM_IV_OK) return zero_iv;
  if (b.lo <= 0.0 && b.hi >= 0.0) {
    *st = GMM_IV_ERR_DIV_ZERO;
    return zero_iv;
  }
  double z0 = gmm_fp_div(a.lo, b.lo);
  double z1 = gmm_fp_div(a.lo, b.hi);
  double z2 = gmm_fp_div(a.hi, b.lo);
  double z3 = gmm_fp_div(a.hi, b.hi);
  if (isnan(z0) || isnan(z1) || isnan(z2) || isnan(z3)) {
    *st = GMM_IV_ERR_ENDPOINT;
    return zero_iv;
  }
  double zmin = z0;
  if (z1 < zmin) zmin = z1;
  if (z2 < zmin) zmin = z2;
  if (z3 < zmin) zmin = z3;
  double zmax = z0;
  if (z1 > zmax) zmax = z1;
  if (z2 > zmax) zmax = z2;
  if (z3 > zmax) zmax = z3;
  double lo = gmm_next_c(zmin, 0);
  double hi = gmm_next_c(zmax, 1);
  return gmm_iv_check_c(lo, hi, st);
}

static inline gmm_iv_t gmm_iv_sqrt_c(gmm_iv_t a, gmm_iv_status_t *st) {
  gmm_iv_t zero_iv = {0.0, 0.0};
  if (*st != GMM_IV_OK) return zero_iv;
  if (a.lo < 0.0) {
    *st = GMM_IV_ERR_NEG_SQRT;
    return zero_iv;
  }
  double lo = 0.0;
  if (a.lo != 0.0) {
    double s_lo = gmm_fp_sqrt(a.lo);
    double n_lo = gmm_next_c(s_lo, 0);
    lo = (n_lo > 0.0) ? n_lo : 0.0;
  }
  double hi = 0.0;
  if (a.hi != 0.0) {
    double s_hi = gmm_fp_sqrt(a.hi);
    hi = gmm_next_c(s_hi, 1);
  }
  return gmm_iv_check_c(lo, hi, st);
}

static inline gmm_iv_t gmm_iv_abs_c(gmm_iv_t a) {
  double a0 = fabs(a.lo);
  double a1 = fabs(a.hi);
  double lo = (a.lo <= 0.0 && a.hi >= 0.0) ? 0.0 : ((a0 < a1) ? a0 : a1);
  double hi = (a0 > a1) ? a0 : a1;
  gmm_iv_t out = {lo, hi};
  return out;
}

static inline gmm_iv_t gmm_iv_dot_strided(
    const double *x, R_xlen_t inc_x,
    const double *y, R_xlen_t inc_y,
    R_xlen_t n, gmm_iv_status_t *st) {
  gmm_iv_t z = {0.0, 0.0};
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    double xv = x[i * inc_x];
    double yv = y[i * inc_y];
    gmm_iv_t xi = {xv, xv};
    gmm_iv_t yi = {yv, yv};
    gmm_iv_t prod = gmm_iv_mul_c(xi, yi, st);
    if (*st != GMM_IV_OK) return z;
    z = gmm_iv_add_c(z, prod, st);
    if (*st != GMM_IV_OK) return z;
  }
  return z;
}

static inline gmm_iv_t gmm_iv_norm_raw(
    const double *lo, const double *hi, R_xlen_t n, gmm_iv_status_t *st) {
  gmm_iv_t zero_iv = {0.0, 0.0};
  if (*st != GMM_IV_OK) return zero_iv;
  double scale = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    double alo = fabs(lo[i]);
    double ahi = fabs(hi[i]);
    double amax = (alo > ahi) ? alo : ahi;
    if (amax > scale) scale = amax;
  }
  if (scale == 0.0) {
    return zero_iv;
  }
  gmm_iv_t scale_iv = {scale, scale};
  gmm_iv_t s = {0.0, 0.0};
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    gmm_iv_t elem = {lo[i], hi[i]};
    gmm_iv_t abs_elem = gmm_iv_abs_c(elem);
    gmm_iv_t t = gmm_iv_div_c(abs_elem, scale_iv, st);
    if (*st != GMM_IV_OK) return zero_iv;
    if (t.lo <= 0.0) t.lo = 0.0;
    gmm_iv_t sq = gmm_iv_mul_c(t, t, st);
    if (*st != GMM_IV_OK) return zero_iv;
    if (sq.lo <= 0.0) sq.lo = 0.0;
    s = gmm_iv_add_c(s, sq, st);
    if (*st != GMM_IV_OK) return zero_iv;
  }
  if (s.lo <= 0.0) s.lo = 0.0;
  gmm_iv_t root_s = gmm_iv_sqrt_c(s, st);
  if (*st != GMM_IV_OK) return zero_iv;
  return gmm_iv_mul_c(scale_iv, root_s, st);
}

/* Registered .Call entrypoints */
SEXP C_gmm_native_status(void);
SEXP C_gmm_test_arithmetic_probe(
    SEXP rmode_sexp, SEXP min_sub_sexp, SEXP even_tie_sexp, SEXP odd_tie_sexp, SEXP fma_res_sexp);
SEXP C_gmm_sexp_aliases(SEXP a_sexp, SEXP b_sexp);
SEXP C_gmm_next(SEXP x_sexp, SEXP up_sexp);
SEXP C_gmm_iv_dot(SEXP x_sexp, SEXP y_sexp, SEXP max_entries_sexp);
SEXP C_gmm_iv_norm(SEXP lo_sexp, SEXP hi_sexp, SEXP max_entries_sexp);
SEXP C_gmm_iv_mat_residual(SEXP X_sexp, SEXP Y_sexp, SEXP D_sexp, SEXP max_entries_sexp);
SEXP C_gmm_iv_eigen_residual(SEXP S_sexp, SEXP Q_sexp, SEXP lam_sexp, SEXP max_entries_sexp);
SEXP C_gmm_iv_mat_product(SEXP xlo_sexp, SEXP xhi_sexp, SEXP ylo_sexp, SEXP yhi_sexp, SEXP max_entries_sexp);
SEXP C_gmm_ipw_if_chunk(
    SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
    SEXP a_start_sexp, SEXP b_len_sexp, SEXP obs_pos_sexp, SEXP obs_g_sexp,
    SEXP obs_w_sexp, SEXP obs_r_sexp, SEXP obs_odds_sexp, SEXP obs_score_sexp, SEXP obs_z_sexp,
    SEXP link_dates_sexp, SEXP link_par_sexp, SEXP dimnames_sexp, SEXP components_sexp,
    SEXP max_entries_sexp);
SEXP C_gmm_cp08_support_chunk(
    SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
    SEXP order_sexp, SEXP a_sexp, SEXP w_sexp, SEXP z_sexp, SEXP ctl_sexp,
    SEXP dates_sexp, SEXP long_double_sexp, SEXP max_entries_sexp);
SEXP C_gmm_cp08_pair_rows(
    SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
    SEXP rows_sexp, SEXP a_sexp, SEXP dates_sexp, SEXP chunk_sexp, SEXP max_entries_sexp);
SEXP C_gmm_cp08_sparse_pair_rows(
    SEXP a_sexp, SEXP post_a_sexp, SEXP post_obs_sexp, SEXP bounds_sexp,
    SEXP read_sexp, SEXP max_entries_sexp);
SEXP C_gmm_cp08_or_if_chunk(
    SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
    SEXP rows_sexp, SEXP ci_sexp, SEXP a_sexp, SEXP w_sexp, SEXP dates_sexp,
    SEXP sup_sexp, SEXP sup_par_sexp, SEXP blk_sexp, SEXP blk_par_sexp,
    SEXP lnk_sexp, SEXP lnk_par_sexp, SEXP dimnames_sexp, SEXP max_entries_sexp);
SEXP C_gmm_cp08_dr_if_prepare(
    SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
    SEXP rows_sexp, SEXP ci_sexp, SEXP a_sexp, SEXP w_sexp, SEXP dates_sexp,
    SEXP sup_sexp, SEXP sup_par_sexp, SEXP mdl_sexp, SEXP mdl_cross_sexp,
    SEXP max_entries_sexp);
SEXP C_gmm_cp08_dr_if_assemble(
    SEXP rows_sexp, SEXP ci_sexp, SEXP a_sexp, SEXP w_sexp, SEXP gather_sexp,
    SEXP sup_sexp, SEXP sup_par_sexp, SEXP mdl_sexp, SEXP mdl_par_sexp,
    SEXP xi_sexp, SEXP rr_sexp, SEXP wp_sexp, SEXP wm_sexp, SEXP lnk_sexp,
    SEXP n_ok_sexp, SEXP dimnames_sexp, SEXP max_entries_sexp);
SEXP C_gmm_cp08_dr_treatment_batch(
    SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
    SEXP rows_sexp, SEXP a_sexp, SEXP X_trt_sexp, SEXP Z_or_sexp, SEXP beta_sexp,
    SEXP wp_sexp, SEXP wm_sexp, SEXP cohort_rows_sexp, SEXP score_sums_sexp,
    SEXP wm_design_sexp, SEXP norm_sexp, SEXP pairs_sexp, SEXP model_par_sexp,
    SEXP ctl_sexp, SEXP max_entries_sexp);
SEXP C_gmm_if_unit_batch(
    SEXP pass_sexp, SEXP psi_sexp, SEXP w_sexp, SEXP flags_sexp, SEXP counts_sexp,
    SEXP weighted_sum_sexp, SEXP raw_sum_sexp, SEXP psi_mean_sexp, SEXP S_sexp,
    SEXP max_entries_sexp);
SEXP C_gmm_cp08_payload(SEXP x_sexp, SEXP scan_attrs_sexp, SEXP call_env_sexp);
void R_init_cdid(DllInfo *dll);

SEXP C_gmm_native_status(void) {
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 10));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
  SEXP hotspots = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(hotspots, 0, Rf_mkChar("C_gmm_iv_dot"));
  SET_STRING_ELT(hotspots, 1, Rf_mkChar("C_gmm_iv_norm"));
  SET_STRING_ELT(hotspots, 2, Rf_mkChar("C_gmm_iv_mat_residual"));
  SET_STRING_ELT(hotspots, 3, Rf_mkChar("C_gmm_iv_eigen_residual"));
  SET_STRING_ELT(hotspots, 4, Rf_mkChar("C_gmm_iv_mat_product"));

  SET_STRING_ELT(names, 0, Rf_mkChar("backend"));
  SET_VECTOR_ELT(out, 0, Rf_mkString("native_c_ieee754_binary64"));
  SET_STRING_ELT(names, 1, Rf_mkChar("c_registration"));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(1));
  SET_STRING_ELT(names, 2, Rf_mkChar("force_symbols"));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(1));
  SET_STRING_ELT(names, 3, Rf_mkChar("fallback_permitted"));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(0));
  SET_STRING_ELT(names, 4, Rf_mkChar("flt_radix"));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(FLT_RADIX));
  SET_STRING_ELT(names, 5, Rf_mkChar("dbl_mant_dig"));
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(DBL_MANT_DIG));
  SET_STRING_ELT(names, 6, Rf_mkChar("rounding_mode"));
  SET_VECTOR_ELT(out, 6, Rf_mkString("FE_TONEAREST"));
  SET_STRING_ELT(names, 7, Rf_mkChar("gradual_underflow_verified"));
  SET_VECTOR_ELT(out, 7, Rf_ScalarLogical(1));
  SET_STRING_ELT(names, 8, Rf_mkChar("volatile_fp_stores"));
  SET_VECTOR_ELT(out, 8, Rf_ScalarLogical(1));
  SET_STRING_ELT(names, 9, Rf_mkChar("compiled_hotspots"));
  SET_VECTOR_ELT(out, 9, hotspots);

  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(3);
  return out;
}

SEXP C_gmm_test_arithmetic_probe(
    SEXP rmode_sexp, SEXP min_sub_sexp, SEXP even_tie_sexp, SEXP odd_tie_sexp, SEXP fma_res_sexp) {
  if (TYPEOF(rmode_sexp) != INTSXP || XLENGTH(rmode_sexp) != 1 ||
      TYPEOF(min_sub_sexp) != REALSXP || XLENGTH(min_sub_sexp) != 1 ||
      TYPEOF(even_tie_sexp) != REALSXP || XLENGTH(even_tie_sexp) != 1 ||
      TYPEOF(odd_tie_sexp) != REALSXP || XLENGTH(odd_tie_sexp) != 1 ||
      TYPEOF(fma_res_sexp) != REALSXP || XLENGTH(fma_res_sexp) != 1) {
    gmm_c_abort(0, "Invalid probe arguments.", "cdid_gmm_control_error");
  }
  int rmode = INTEGER(rmode_sexp)[0];
#if defined(FE_TONEAREST)
  if (rmode < 0) rmode = FE_TONEAREST;
#endif
  const char *err_msg = NULL;
  if (!gmm_verify_arithmetic_values(
          rmode, REAL(min_sub_sexp)[0], REAL(even_tie_sexp)[0],
          REAL(odd_tie_sexp)[0], REAL(fma_res_sexp)[0], &err_msg)) {
    gmm_c_abort(0, err_msg, "cdid_gmm_native_error");
  }
  return Rf_ScalarLogical(1);
}

/* Object identity diagnostic only. All production operands are read-only.
   Do not force ALTREP data merely to diagnose identity. */
SEXP C_gmm_sexp_aliases(SEXP a_sexp, SEXP b_sexp) {
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  return Rf_ScalarLogical(a_sexp != R_NilValue && a_sexp == b_sexp);
}

SEXP C_gmm_next(SEXP x_sexp, SEXP up_sexp) {
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  if (Rf_isObject(x_sexp) || TYPEOF(x_sexp) != REALSXP ||
      Rf_getAttrib(x_sexp, R_DimSymbol) != R_NilValue || XLENGTH(x_sexp) != 1) {
    gmm_c_abort(0, "Invalid interval endpoint.", "cdid_gmm_accuracy_error");
  }
  if (Rf_isObject(up_sexp) || TYPEOF(up_sexp) != LGLSXP ||
      Rf_getAttrib(up_sexp, R_DimSymbol) != R_NilValue || XLENGTH(up_sexp) != 1 ||
      LOGICAL(up_sexp)[0] == NA_LOGICAL) {
    gmm_c_abort(0, "Invalid interval direction.", "cdid_gmm_accuracy_error");
  }
  double x = REAL(x_sexp)[0];
  if (ISNA(x) || ISNAN(x)) {
    gmm_c_abort(0, "Invalid interval endpoint.", "cdid_gmm_accuracy_error");
  }
  double out = gmm_next_c(x, LOGICAL(up_sexp)[0] != 0);
  return Rf_ScalarReal(out);
}

SEXP C_gmm_iv_dot(SEXP x_sexp, SEXP y_sexp, SEXP max_entries_sexp) {
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  R_xlen_t nx = 0, ny = 0;
  const double *px = NULL, *py = NULL;
  gmm_validate_real_vector_header(0, x_sexp, "x", max_entries, &nx, &px);
  gmm_validate_real_vector_header(0, y_sexp, "y", max_entries, &ny, &py);
  if (nx != ny) {
    gmm_c_abort(0, "`x` and `y` in `gmm_iv_dot` must have identical lengths.", "cdid_gmm_bound_error");
  }
  (void) gmm_checked_mul_u64(0, 2ULL, (uint64_t) nx, max_entries, "gmm_iv_dot input pair (x + y)");
  px = REAL(x_sexp); py = REAL(y_sexp);
  gmm_validate_real_finite(0, px, nx, "x");
  gmm_validate_real_finite(0, py, ny, "y");

  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  double *pout = REAL(out);
  if (gmm_buffers_overlap(pout, 2, px, nx) || gmm_buffers_overlap(pout, 2, py, ny)) {
    gmm_c_abort(1, "Allocated output aliases input in `gmm_iv_dot`.", "cdid_gmm_bound_error");
  }

  gmm_iv_status_t st = GMM_IV_OK;
  gmm_iv_t z = gmm_iv_dot_strided(px, 1, py, 1, nx, &st);
  if (st != GMM_IV_OK) {
    gmm_c_abort_status(1, st);
  }
  pout[0] = z.lo;
  pout[1] = z.hi;
  UNPROTECT(1);
  return out;
}

SEXP C_gmm_iv_norm(SEXP lo_sexp, SEXP hi_sexp, SEXP max_entries_sexp) {
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  if (lo_sexp == R_NilValue || Rf_isObject(lo_sexp) || TYPEOF(lo_sexp) != REALSXP) {
    gmm_c_abort(0, "`lo` in `gmm_iv_norm` must be a non-object double vector or matrix.",
                "cdid_gmm_accuracy_error");
  }
  SEXP dim_lo = Rf_getAttrib(lo_sexp, R_DimSymbol);
  R_xlen_t n = XLENGTH(lo_sexp);
  if (n < 1) {
    gmm_c_abort(0, "Invalid empty input length in `gmm_iv_norm`.", "cdid_gmm_bound_error");
  }
  if (dim_lo != R_NilValue) {
    if (TYPEOF(dim_lo) != INTSXP || XLENGTH(dim_lo) != 2 ||
        INTEGER(dim_lo)[0] == NA_INTEGER || INTEGER(dim_lo)[1] == NA_INTEGER ||
        INTEGER(dim_lo)[0] < 1 || INTEGER(dim_lo)[1] < 1) {
      gmm_c_abort(0, "Invalid matrix dimensions on `lo` in `gmm_iv_norm`.", "cdid_gmm_bound_error");
    }
    uint64_t prod = gmm_checked_mul_u64(0, (uint64_t) INTEGER(dim_lo)[0], (uint64_t) INTEGER(dim_lo)[1],
                                        max_entries, "gmm_iv_norm(lo)");
    if ((R_xlen_t) prod != n) {
      gmm_c_abort(0, "Dimension product on `lo` disagrees with length.", "cdid_gmm_bound_error");
    }
  } else {
    (void) gmm_checked_mul_u64(0, (uint64_t) n, 1ULL, max_entries, "gmm_iv_norm(lo)");
  }

  const double *plo = NULL;
  const double *phi = NULL;
  if (hi_sexp != R_NilValue) {
    if (Rf_isObject(hi_sexp) || TYPEOF(hi_sexp) != REALSXP) {
      gmm_c_abort(0, "`hi` in `gmm_iv_norm` must be a non-object double vector or matrix.",
                  "cdid_gmm_accuracy_error");
    }
    if (XLENGTH(hi_sexp) != n) {
      gmm_c_abort(0, "`lo` and `hi` in `gmm_iv_norm` must have identical lengths.", "cdid_gmm_bound_error");
    }
    SEXP dim_hi = Rf_getAttrib(hi_sexp, R_DimSymbol);
    if ((dim_lo == R_NilValue) != (dim_hi == R_NilValue)) {
      gmm_c_abort(0, "`lo` and `hi` in `gmm_iv_norm` must have matching dimension structure.",
                  "cdid_gmm_bound_error");
    }
    if (dim_lo != R_NilValue) {
      if (TYPEOF(dim_hi) != INTSXP || XLENGTH(dim_hi) != 2 ||
          INTEGER(dim_hi)[0] != INTEGER(dim_lo)[0] ||
          INTEGER(dim_hi)[1] != INTEGER(dim_lo)[1]) {
        gmm_c_abort(0, "`lo` and `hi` matrix dimensions in `gmm_iv_norm` do not match.",
                    "cdid_gmm_bound_error");
      }
    }
    (void) gmm_checked_mul_u64(0, 2ULL, (uint64_t) n, max_entries, "gmm_iv_norm(lo + hi)");

  }

  plo = REAL(lo_sexp); phi = hi_sexp == R_NilValue ? plo : REAL(hi_sexp);
  gmm_validate_real_finite(0, plo, n, "lo");
  if (hi_sexp != R_NilValue) {
    for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
      if (!R_finite(phi[i])) {
        gmm_c_abort(0, "Nonfinite element in `hi` for `gmm_iv_norm`.", "cdid_gmm_accuracy_error");
      }
      if (plo[i] > phi[i]) {
        gmm_c_abort(0, "Inverted interval endpoints (`lo > hi`) in `gmm_iv_norm`.", "cdid_gmm_accuracy_error");
      }
    }
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  double *pout = REAL(out);
  if (gmm_buffers_overlap(pout, 2, plo, n) || (phi != plo && gmm_buffers_overlap(pout, 2, phi, n))) {
    gmm_c_abort(1, "Allocated output aliases input in `gmm_iv_norm`.", "cdid_gmm_bound_error");
  }

  gmm_iv_status_t st = GMM_IV_OK;
  gmm_iv_t z = gmm_iv_norm_raw(plo, phi, n, &st);
  if (st != GMM_IV_OK) {
    gmm_c_abort_status(1, st);
  }
  pout[0] = z.lo;
  pout[1] = z.hi;
  UNPROTECT(1);
  return out;
}

SEXP C_gmm_iv_mat_residual(SEXP X_sexp, SEXP Y_sexp, SEXP D_sexp, SEXP max_entries_sexp) {
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  int nr_x = 0, nc_x = 0, nr_y = 0, nc_y = 0, nr_d = 0, nc_d = 0;
  R_xlen_t len_x = 0, len_y = 0, len_d = 0;
  const double *px = NULL, *py = NULL, *pd = NULL;

  gmm_validate_real_matrix_header(0, X_sexp, "X", max_entries, &nr_x, &nc_x, &len_x, &px);
  gmm_validate_real_matrix_header(0, Y_sexp, "Y", max_entries, &nr_y, &nc_y, &len_y, &py);
  gmm_validate_real_matrix_header(0, D_sexp, "D", max_entries, &nr_d, &nc_d, &len_d, &pd);

  if (nc_x != nr_y || nr_d != nr_x || nc_d != nc_y) {
    gmm_c_abort(0, "Conformability error in `gmm_iv_mat_residual`: require ncol(X)==nrow(Y), dim(D)==c(nrow(X),ncol(Y)).",
                "cdid_gmm_bound_error");
  }
  (void) gmm_checked_mul_u64(0, 2ULL, (uint64_t) len_d, max_entries,
                                         "gmm_iv_mat_residual interval output (lo + hi)");

  px = REAL(X_sexp); py = REAL(Y_sexp); pd = REAL(D_sexp);
  gmm_validate_real_finite(0, px, len_x, "X");
  gmm_validate_real_finite(0, py, len_y, "Y");
  gmm_validate_real_finite(0, pd, len_d, "D");

  SEXP lo_sexp = PROTECT(Rf_allocMatrix(REALSXP, nr_x, nc_y));
  SEXP hi_sexp = PROTECT(Rf_allocMatrix(REALSXP, nr_x, nc_y));
  SEXP norm_sexp = PROTECT(Rf_allocVector(REALSXP, 2));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  const int nprot = 5;

  double *plo = REAL(lo_sexp);
  double *phi = REAL(hi_sexp);
  double *pnorm = REAL(norm_sexp);

  if (gmm_buffers_overlap(plo, len_d, phi, len_d) ||
      gmm_buffers_overlap(plo, len_d, px, len_x) ||
      gmm_buffers_overlap(plo, len_d, py, len_y) ||
      gmm_buffers_overlap(plo, len_d, pd, len_d) ||
      gmm_buffers_overlap(phi, len_d, px, len_x) ||
      gmm_buffers_overlap(phi, len_d, py, len_y) ||
      gmm_buffers_overlap(phi, len_d, pd, len_d)) {
    gmm_c_abort(nprot, "Allocated output buffers alias input buffers in `gmm_iv_mat_residual`.",
                "cdid_gmm_bound_error");
  }

  gmm_iv_status_t st = GMM_IV_OK;
  R_xlen_t nr = (R_xlen_t) nr_x;
  R_xlen_t nk = (R_xlen_t) nc_x;
  R_xlen_t nc = (R_xlen_t) nc_y;

  for (R_xlen_t i = 0; i < nr; ++i) {
    if ((i % 128) == 0) R_CheckUserInterrupt();
    for (R_xlen_t j = 0; j < nc; ++j) {
      if (j > 0 && (j % 4096) == 0) R_CheckUserInterrupt();
      R_xlen_t idx = i + j * nr;
      gmm_iv_t dot_ij = gmm_iv_dot_strided(px + i, nr, py + j * nk, 1, nk, &st);
      if (st != GMM_IV_OK) gmm_c_abort_status(nprot, st);
      gmm_iv_t d_ij = {pd[idx], pd[idx]};
      gmm_iv_t z = gmm_iv_sub_c(dot_ij, d_ij, &st);
      if (st != GMM_IV_OK) gmm_c_abort_status(nprot, st);
      plo[idx] = z.lo;
      phi[idx] = z.hi;
    }
  }

  gmm_iv_t norm_iv = gmm_iv_norm_raw(plo, phi, len_d, &st);
  if (st != GMM_IV_OK) {
    gmm_c_abort_status(nprot, st);
  }
  pnorm[0] = norm_iv.lo;
  pnorm[1] = norm_iv.hi;

  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("lo"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("hi"));
  SET_STRING_ELT(names_sexp, 2, Rf_mkChar("norm"));
  SET_VECTOR_ELT(out_sexp, 0, lo_sexp);
  SET_VECTOR_ELT(out_sexp, 1, hi_sexp);
  SET_VECTOR_ELT(out_sexp, 2, norm_sexp);
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);

  UNPROTECT(nprot);
  return out_sexp;
}

SEXP C_gmm_iv_eigen_residual(SEXP S_sexp, SEXP Q_sexp, SEXP lam_sexp, SEXP max_entries_sexp) {
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  int nr_s = 0, nc_s = 0, nr_q = 0, nc_q = 0;
  R_xlen_t len_s = 0, len_q = 0, len_lam = 0;
  const double *ps = NULL, *pq = NULL, *plam = NULL;

  gmm_validate_real_matrix_header(0, S_sexp, "S", max_entries, &nr_s, &nc_s, &len_s, &ps);
  gmm_validate_real_matrix_header(0, Q_sexp, "Q", max_entries, &nr_q, &nc_q, &len_q, &pq);
  gmm_validate_real_vector_header(0, lam_sexp, "lam", max_entries, &len_lam, &plam);

  if (nr_s != nc_s || nr_q != nc_q || nr_s != nr_q || len_lam != (R_xlen_t) nr_s) {
    gmm_c_abort(0, "Conformability error in `C_gmm_iv_eigen_residual`: require m x m `S`, `Q` and length-m `lam`.",
                "cdid_gmm_bound_error");
  }
  (void) gmm_checked_mul_u64(0, 2ULL, (uint64_t) len_s, max_entries,
                                         "C_gmm_iv_eigen_residual interval output (lo + hi)");

  ps = REAL(S_sexp); pq = REAL(Q_sexp); plam = REAL(lam_sexp);
  gmm_validate_real_finite(0, ps, len_s, "S");
  gmm_validate_real_finite(0, pq, len_q, "Q");
  gmm_validate_real_finite(0, plam, len_lam, "lam");

  SEXP lo_sexp = PROTECT(Rf_allocMatrix(REALSXP, nr_s, nc_s));
  SEXP hi_sexp = PROTECT(Rf_allocMatrix(REALSXP, nr_s, nc_s));
  SEXP norm_sexp = PROTECT(Rf_allocVector(REALSXP, 2));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  const int nprot = 5;

  double *plo = REAL(lo_sexp);
  double *phi = REAL(hi_sexp);
  double *pnorm = REAL(norm_sexp);

  if (gmm_buffers_overlap(plo, len_s, phi, len_s) ||
      gmm_buffers_overlap(plo, len_s, ps, len_s) ||
      gmm_buffers_overlap(plo, len_s, pq, len_q) ||
      gmm_buffers_overlap(plo, len_s, plam, len_lam) ||
      gmm_buffers_overlap(phi, len_s, ps, len_s) ||
      gmm_buffers_overlap(phi, len_s, pq, len_q) ||
      gmm_buffers_overlap(phi, len_s, plam, len_lam)) {
    gmm_c_abort(nprot, "Allocated output buffers alias input buffers in `C_gmm_iv_eigen_residual`.",
                "cdid_gmm_bound_error");
  }

  gmm_iv_status_t st = GMM_IV_OK;
  R_xlen_t m = (R_xlen_t) nr_s;

  for (R_xlen_t i = 0; i < m; ++i) {
    if ((i % 128) == 0) R_CheckUserInterrupt();
    for (R_xlen_t j = 0; j < m; ++j) {
      R_xlen_t idx = i + j * m;
      gmm_iv_t dot_ij = gmm_iv_dot_strided(ps + i, m, pq + j * m, 1, m, &st);
      if (st != GMM_IV_OK) gmm_c_abort_status(nprot, st);
      gmm_iv_t q_ij = {pq[idx], pq[idx]};
      gmm_iv_t l_j = {plam[j], plam[j]};
      gmm_iv_t ql_ij = gmm_iv_mul_c(q_ij, l_j, &st);
      if (st != GMM_IV_OK) gmm_c_abort_status(nprot, st);
      gmm_iv_t z = gmm_iv_sub_c(dot_ij, ql_ij, &st);
      if (st != GMM_IV_OK) gmm_c_abort_status(nprot, st);
      plo[idx] = z.lo;
      phi[idx] = z.hi;
    }
  }

  gmm_iv_t norm_iv = gmm_iv_norm_raw(plo, phi, len_s, &st);
  if (st != GMM_IV_OK) {
    gmm_c_abort_status(nprot, st);
  }
  pnorm[0] = norm_iv.lo;
  pnorm[1] = norm_iv.hi;

  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("lo"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("hi"));
  SET_STRING_ELT(names_sexp, 2, Rf_mkChar("norm"));
  SET_VECTOR_ELT(out_sexp, 0, lo_sexp);
  SET_VECTOR_ELT(out_sexp, 1, hi_sexp);
  SET_VECTOR_ELT(out_sexp, 2, norm_sexp);
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);

  UNPROTECT(nprot);
  return out_sexp;
}

/* Directed interval matrix product (CP02 pipeline `ipw_pipeline_interval_product`).
   Given interval matrices X = [xlo, xhi] (nr x nk) and Y = [ylo, yhi] (nk x nc), the
   entry (i, j) is the same-order recurrence of the frozen R reference:
     z_0 = [0, 0];  z_k = gmm_iv_add_c(z_{k-1}, gmm_iv_mul_c(X[i,k], Y[k,j]))  for k = 1..nk,
   with the R reference's two zero shortcuts (an all-zero row i of X or an all-zero
   column j of Y leaves the output entry at [0, 0]). Only the accepted primitives
   gmm_iv_mul_c / gmm_iv_add_c are used: no reassociation, no BLAS, no fused
   operations, no clipping. Failures carry the primitives' typed status. */
SEXP C_gmm_iv_mat_product(SEXP xlo_sexp, SEXP xhi_sexp, SEXP ylo_sexp, SEXP yhi_sexp, SEXP max_entries_sexp) {
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  int nr_xl = 0, nc_xl = 0, nr_xh = 0, nc_xh = 0, nr_yl = 0, nc_yl = 0, nr_yh = 0, nc_yh = 0;
  R_xlen_t len_xl = 0, len_xh = 0, len_yl = 0, len_yh = 0;
  const double *pxl = NULL, *pxh = NULL, *pyl = NULL, *pyh = NULL;

  gmm_validate_real_matrix_header(0, xlo_sexp, "xlo", max_entries, &nr_xl, &nc_xl, &len_xl, &pxl);
  gmm_validate_real_matrix_header(0, xhi_sexp, "xhi", max_entries, &nr_xh, &nc_xh, &len_xh, &pxh);
  gmm_validate_real_matrix_header(0, ylo_sexp, "ylo", max_entries, &nr_yl, &nc_yl, &len_yl, &pyl);
  gmm_validate_real_matrix_header(0, yhi_sexp, "yhi", max_entries, &nr_yh, &nc_yh, &len_yh, &pyh);

  if (nr_xl != nr_xh || nc_xl != nc_xh) {
    gmm_c_abort(0, "`xlo` and `xhi` in `gmm_iv_mat_product` must have identical dimensions.", "cdid_gmm_bound_error");
  }
  if (nr_yl != nr_yh || nc_yl != nc_yh) {
    gmm_c_abort(0, "`ylo` and `yhi` in `gmm_iv_mat_product` must have identical dimensions.", "cdid_gmm_bound_error");
  }
  if (nc_xl != nr_yl) {
    gmm_c_abort(0, "Conformability error in `gmm_iv_mat_product`: require ncol(x) == nrow(y).", "cdid_gmm_bound_error");
  }
  uint64_t out_len_u = gmm_checked_mul_u64(0, (uint64_t) nr_xl, (uint64_t) nc_yl, max_entries,
                                           "gmm_iv_mat_product output (nrow(x) x ncol(y))");
  (void) gmm_checked_mul_u64(0, 2ULL, out_len_u, max_entries, "gmm_iv_mat_product interval output (lo + hi)");

  pxl = REAL(xlo_sexp); pxh = REAL(xhi_sexp); pyl = REAL(ylo_sexp); pyh = REAL(yhi_sexp);
  gmm_validate_real_finite(0, pxl, len_xl, "xlo");
  gmm_validate_real_finite(0, pxh, len_xh, "xhi");
  gmm_validate_real_finite(0, pyl, len_yl, "ylo");
  gmm_validate_real_finite(0, pyh, len_yh, "yhi");
  for (R_xlen_t i = 0; i < len_xl; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    if (pxl[i] > pxh[i]) {
      gmm_c_abort(0, "Inverted interval endpoints (`lo > hi`) in `gmm_iv_mat_product` operand `x`.", "cdid_gmm_accuracy_error");
    }
  }
  for (R_xlen_t i = 0; i < len_yl; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    if (pyl[i] > pyh[i]) {
      gmm_c_abort(0, "Inverted interval endpoints (`lo > hi`) in `gmm_iv_mat_product` operand `y`.", "cdid_gmm_accuracy_error");
    }
  }

  R_xlen_t nr = (R_xlen_t) nr_xl;
  R_xlen_t nk = (R_xlen_t) nc_xl;
  R_xlen_t nc = (R_xlen_t) nc_yl;
  R_xlen_t out_len = (R_xlen_t) out_len_u;

  SEXP lo_sexp = PROTECT(Rf_allocMatrix(REALSXP, nr_xl, nc_yl));
  SEXP hi_sexp = PROTECT(Rf_allocMatrix(REALSXP, nr_xl, nc_yl));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
  const int nprot = 4;

  double *plo = REAL(lo_sexp);
  double *phi = REAL(hi_sexp);

  if (gmm_buffers_overlap(plo, out_len, phi, out_len) ||
      gmm_buffers_overlap(plo, out_len, pxl, len_xl) || gmm_buffers_overlap(plo, out_len, pxh, len_xh) ||
      gmm_buffers_overlap(plo, out_len, pyl, len_yl) || gmm_buffers_overlap(plo, out_len, pyh, len_yh) ||
      gmm_buffers_overlap(phi, out_len, pxl, len_xl) || gmm_buffers_overlap(phi, out_len, pxh, len_xh) ||
      gmm_buffers_overlap(phi, out_len, pyl, len_yl) || gmm_buffers_overlap(phi, out_len, pyh, len_yh)) {
    gmm_c_abort(nprot, "Allocated output buffers alias input buffers in `gmm_iv_mat_product`.", "cdid_gmm_bound_error");
  }

  /* The R reference initialises both outputs to +0.0 and skips shortcut entries. */
  for (R_xlen_t idx = 0; idx < out_len; ++idx) {
    plo[idx] = 0.0;
    phi[idx] = 0.0;
  }

  /* Zero shortcuts of the R reference, evaluated once per row of X / column of Y:
     all(xlo[i, ] == 0 & xhi[i, ] == 0) and all(ylo[, j] == 0 & yhi[, j] == 0).
     Transient flags only (nr + nc ints); no full-unit/moment Cartesian allocation. */
  int *x_zero_row = (int *) R_alloc((size_t) nr, sizeof(int));
  int *y_zero_col = (int *) R_alloc((size_t) nc, sizeof(int));
  for (R_xlen_t i = 0; i < nr; ++i) {
    int all_zero = 1;
    for (R_xlen_t k = 0; k < nk; ++k) {
      R_xlen_t idx = i + k * nr;
      if (!(pxl[idx] == 0.0 && pxh[idx] == 0.0)) { all_zero = 0; break; }
    }
    x_zero_row[i] = all_zero;
  }
  for (R_xlen_t j = 0; j < nc; ++j) {
    int all_zero = 1;
    for (R_xlen_t k = 0; k < nk; ++k) {
      R_xlen_t idx = k + j * nk;
      if (!(pyl[idx] == 0.0 && pyh[idx] == 0.0)) { all_zero = 0; break; }
    }
    y_zero_col[j] = all_zero;
  }

  gmm_iv_status_t st = GMM_IV_OK;
  for (R_xlen_t i = 0; i < nr; ++i) {
    if ((i % 128) == 0) R_CheckUserInterrupt();
    for (R_xlen_t j = 0; j < nc; ++j) {
      if (j > 0 && (j % 4096) == 0) R_CheckUserInterrupt();
      if (y_zero_col[j] || x_zero_row[i]) continue;
      gmm_iv_t z = {0.0, 0.0};
      for (R_xlen_t k = 0; k < nk; ++k) {
        R_xlen_t ix = i + k * nr;
        R_xlen_t iy = k + j * nk;
        gmm_iv_t xi = {pxl[ix], pxh[ix]};
        gmm_iv_t yi = {pyl[iy], pyh[iy]};
        gmm_iv_t prod = gmm_iv_mul_c(xi, yi, &st);
        if (st != GMM_IV_OK) gmm_c_abort_status(nprot, st);
        z = gmm_iv_add_c(z, prod, &st); /* IPN_ACCUMULATE */
        if (st != GMM_IV_OK) gmm_c_abort_status(nprot, st);
      }
      R_xlen_t idx = i + j * nr;
      plo[idx] = z.lo;
      phi[idx] = z.hi;
    }
  }

  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("lo"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("hi"));
  SET_VECTOR_ELT(out_sexp, 0, lo_sexp);
  SET_VECTOR_ELT(out_sexp, 1, hi_sexp);
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);

  UNPROTECT(nprot);
  return out_sexp;
}

/* Bounded all-links complete point-IF chunk of the trusted IPW batch plan
   (`ipw_point_if_plan_chunk`, R/chained_ipw_point_view.R). One call writes the
   B x L IF of every retained link of one fitted nuisance for the A-units
   a_start .. a_start + B - 1, in place of L chunk queries and R column
   assemblies. R supplies the chunk's shared observer bases exactly as the
   generic replay's one-nuisance chunk cache forms them (w, r = rates[g], and
   odds/score base from R's own exp and matrix product) and the sealed links'
   compact descriptors (dates, treated/control means, normalizers, the p beta
   loadings and rate_coef_s1, (4 + p + m) x L; p counts every propensity design
   column including the intercept). Every risk row must be a pair observer at
   both dates of every link: posting slot date_start[d] - 1 + a is verified to
   hold A-unit a, for all links, before any outcome endpoint is read; off-risk
   rows and nonobservers are never read. Per link, in the generic replay's
   order: endpoint finiteness, the raw y_t - y_s and its finiteness, then the
   four signed terms ((treated + control) + beta_score) + rate with explicit
   volatile operations, then IF finiteness. The beta loading factor is R's
   unchanged single-column recurrence 0 + lb * z for p = 1. For p > 1,
   evaluate the same base R matrix-product primitive as generic replay on
   borrowed Z and one private loading vector. This preserves the active R
   matprod implementation, its accumulation and nonfinite behavior; it is
   unconditional arithmetic, with no caught-failure alternate path. Data-dependent failures
   return status (code, link, row) so R raises the generic route's exact
   condition; malformed arguments abort here. Borrowed inputs are never written. */
#define GMM_IPW_IF_OK 0
#define GMM_IPW_IF_ENDPOINT 1
#define GMM_IPW_IF_DIFFERENCE 2
#define GMM_IPW_IF_NONFINITE 3
#define GMM_IPW_IF_POSTING 4

static R_xlen_t gmm_ipw_vector_length(int nprotect, SEXP x, int type, const char *name) {
  if (x == R_NilValue || Rf_isObject(x) || TYPEOF(x) != type) {
    char buf[200];
    snprintf(buf, sizeof(buf), "`%s` in `C_gmm_ipw_if_chunk` must be a non-object %s vector.",
             name, type == INTSXP ? "integer" : "double");
    gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
  }
  return XLENGTH(x);
}

static int gmm_ipw_int_scalar(int nprotect, SEXP x, const char *name, int min_val) {
  if (gmm_ipw_vector_length(nprotect, x, INTSXP, name) != 1 ||
      INTEGER(x)[0] == NA_INTEGER || INTEGER(x)[0] < min_val) {
    char buf[200];
    snprintf(buf, sizeof(buf), "`%s` in `C_gmm_ipw_if_chunk` must be one integer >= %d.", name, min_val);
    gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
  }
  return INTEGER(x)[0];
}

static void gmm_ipw_matrix_dims(int nprotect, SEXP x, int type, const char *name, int *nr, int *nc) {
  R_xlen_t len = gmm_ipw_vector_length(nprotect, x, type, name);
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (dim == R_NilValue || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
      INTEGER(dim)[0] == NA_INTEGER || INTEGER(dim)[1] == NA_INTEGER ||
      INTEGER(dim)[0] < 1 || INTEGER(dim)[1] < 1 ||
      (R_xlen_t) INTEGER(dim)[0] * (R_xlen_t) INTEGER(dim)[1] != len) {
    char buf[200];
    snprintf(buf, sizeof(buf), "`%s` in `C_gmm_ipw_if_chunk` must be a matrix with positive dimensions.", name);
    gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
  }
  *nr = INTEGER(dim)[0];
  *nc = INTEGER(dim)[1];
}

static inline double gmm_fp_neg(double x) {
  volatile double vx = x;
  volatile double vz = -vx;
  return vz;
}

SEXP C_gmm_ipw_if_chunk(
    SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
    SEXP a_start_sexp, SEXP b_len_sexp, SEXP obs_pos_sexp, SEXP obs_g_sexp,
    SEXP obs_w_sexp, SEXP obs_r_sexp, SEXP obs_odds_sexp, SEXP obs_score_sexp, SEXP obs_z_sexp,
    SEXP link_dates_sexp, SEXP link_par_sexp, SEXP dimnames_sexp, SEXP components_sexp,
    SEXP max_entries_sexp) {
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);

  R_xlen_t ny = gmm_ipw_vector_length(0, y_sexp, REALSXP, "y");
  R_xlen_t n_dates = gmm_ipw_vector_length(0, date_start_sexp, INTSXP, "date_start");
  R_xlen_t n_post = gmm_ipw_vector_length(0, post_a_sexp, INTSXP, "post_a_index");
  if (gmm_ipw_vector_length(0, post_obs_sexp, INTSXP, "post_obs_index") != n_post) {
    gmm_c_abort(0, "`post_a_index` and `post_obs_index` in `C_gmm_ipw_if_chunk` differ in length.", "cdid_gmm_bound_error");
  }
  int a_start = gmm_ipw_int_scalar(0, a_start_sexp, "a_start", 1);
  int b_len = gmm_ipw_int_scalar(0, b_len_sexp, "b_len", 1);
  if ((int64_t) a_start + (int64_t) b_len - 1 > (int64_t) INT_MAX) {
    gmm_c_abort(0, "Chunk range in `C_gmm_ipw_if_chunk` exceeds the integer A-unit range.", "cdid_gmm_bound_error");
  }
  R_xlen_t k = gmm_ipw_vector_length(0, obs_pos_sexp, INTSXP, "obs_pos");
  if (k > (R_xlen_t) b_len ||
      gmm_ipw_vector_length(0, obs_g_sexp, INTSXP, "obs_g") != k ||
      gmm_ipw_vector_length(0, obs_w_sexp, REALSXP, "obs_w") != k ||
      gmm_ipw_vector_length(0, obs_r_sexp, REALSXP, "obs_r") != k ||
      gmm_ipw_vector_length(0, obs_odds_sexp, REALSXP, "obs_odds") != k ||
      gmm_ipw_vector_length(0, obs_score_sexp, REALSXP, "obs_score_base") != k) {
    gmm_c_abort(0, "Observer arrays in `C_gmm_ipw_if_chunk` must share one length no larger than `b_len`.",
                "cdid_gmm_bound_error");
  }
  /* Shared observer design: a k x p column-major double matrix (p >= 1, k may be
     0); a dimensionless vector is the one-column (p = 1) design. Its length must
     be the checked product k * p before any descriptor or row is read. */
  R_xlen_t z_len = gmm_ipw_vector_length(0, obs_z_sexp, REALSXP, "obs_z");
  SEXP z_dim = Rf_getAttrib(obs_z_sexp, R_DimSymbol);
  int p = 1;
  if (z_dim != R_NilValue) {
    if (TYPEOF(z_dim) != INTSXP || XLENGTH(z_dim) != 2 ||
        INTEGER(z_dim)[0] == NA_INTEGER || INTEGER(z_dim)[1] == NA_INTEGER ||
        INTEGER(z_dim)[0] != (int) k || INTEGER(z_dim)[1] < 1) {
      gmm_c_abort(0, "`obs_z` in `C_gmm_ipw_if_chunk` must be a k x p design matrix (p >= 1 columns).",
                  "cdid_gmm_bound_error");
    }
    p = INTEGER(z_dim)[1];
  }
  uint64_t z_len_u = gmm_checked_mul_u64(0, (uint64_t) k, (uint64_t) p, max_entries,
                                         "C_gmm_ipw_if_chunk shared observer design (k x p)");
  if (z_len < 0 || (uint64_t) z_len != z_len_u) {
    gmm_c_abort(0, "`obs_z` in `C_gmm_ipw_if_chunk` must hold exactly k x p entries.", "cdid_gmm_bound_error");
  }
  int nr_dates = 0, n_links = 0, n_par = 0, n_links_par = 0;
  gmm_ipw_matrix_dims(0, link_dates_sexp, INTSXP, "link_dates", &nr_dates, &n_links);
  gmm_ipw_matrix_dims(0, link_par_sexp, REALSXP, "link_par", &n_par, &n_links_par);
  if (nr_dates != 2 || n_par < 5 || p > n_par - 5 || n_links_par != n_links) {
    gmm_c_abort(0, "Link descriptors in `C_gmm_ipw_if_chunk` must be 2 x L dates and (4 + p + m) x L parameters.",
                "cdid_gmm_bound_error");
  }
  int m = n_par - 4 - p;
  if (Rf_isObject(components_sexp) || TYPEOF(components_sexp) != LGLSXP || XLENGTH(components_sexp) != 1 ||
      LOGICAL(components_sexp)[0] == NA_LOGICAL) {
    gmm_c_abort(0, "`include_components` in `C_gmm_ipw_if_chunk` must be TRUE or FALSE.", "cdid_gmm_bound_error");
  }
  int with_components = LOGICAL(components_sexp)[0] != 0;
  if (dimnames_sexp != R_NilValue) {
    if (TYPEOF(dimnames_sexp) != VECSXP || XLENGTH(dimnames_sexp) != 2 ||
        VECTOR_ELT(dimnames_sexp, 0) != R_NilValue ||
        TYPEOF(VECTOR_ELT(dimnames_sexp, 1)) != STRSXP ||
        XLENGTH(VECTOR_ELT(dimnames_sexp, 1)) != (R_xlen_t) n_links) {
      gmm_c_abort(0, "`dimnames` in `C_gmm_ipw_if_chunk` must be NULL or list(NULL, <L link IDs>).",
                  "cdid_gmm_bound_error");
    }
  }
  uint64_t out_len_u = gmm_checked_mul_u64(0, (uint64_t) b_len, (uint64_t) n_links, max_entries,
                                           "C_gmm_ipw_if_chunk output (B x L)");
  (void) gmm_checked_mul_u64(0, with_components ? 5ULL : 1ULL, out_len_u, max_entries,
                             "C_gmm_ipw_if_chunk outputs (IF plus optional components)");

  const double *py = REAL(y_sexp);
  const int *pds = INTEGER(date_start_sexp);
  const int *ppa = INTEGER(post_a_sexp);
  const int *ppo = INTEGER(post_obs_sexp);
  const int *opos = INTEGER(obs_pos_sexp);
  const int *og = INTEGER(obs_g_sexp);
  const double *ow = REAL(obs_w_sexp);
  const double *orr = REAL(obs_r_sexp);
  const double *oodds = REAL(obs_odds_sexp);
  const double *osb = REAL(obs_score_sexp);
  const double *oz = REAL(obs_z_sexp);
  const int *pld = INTEGER(link_dates_sexp);
  const double *plp = REAL(link_par_sexp);

  /* One private p-vector and one k-vector result, reused/released per link.
     R composes B + p entries with the producer before any allocation. */
  if (p > 1 && k > 0) {
    (void) gmm_checked_mul_u64(0, 1ULL, (uint64_t) k + (uint64_t) p, max_entries,
                               "C_gmm_ipw_if_chunk exact matrix-product scratch (k + p)");
  }

  for (R_xlen_t i = 0; i < k; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    if (opos[i] == NA_INTEGER || opos[i] < 1 || opos[i] > b_len || (i > 0 && opos[i] <= opos[i - 1]) ||
        og[i] == NA_INTEGER || og[i] < 1 || og[i] > m) {
      gmm_c_abort(0, "Observer positions in `C_gmm_ipw_if_chunk` must strictly increase within the chunk, with cohort indices in 1..m.",
                  "cdid_gmm_bound_error");
    }
  }
  for (int ell = 0; ell < n_links; ++ell) {
    for (int side = 0; side < 2; ++side) {
      int d = pld[2 * (R_xlen_t) ell + side];
      if (d == NA_INTEGER || d < 1 || (R_xlen_t) d > n_dates) {
        gmm_c_abort(0, "Link calendar indices in `C_gmm_ipw_if_chunk` are outside the posting dates.",
                    "cdid_gmm_bound_error");
      }
    }
  }

  R_xlen_t out_len = (R_xlen_t) out_len_u;
  int n_out = with_components ? 5 : 1;
  SEXP mats[5];
  int nprot = 0;
  for (int j = 0; j < n_out; ++j) {
    mats[j] = PROTECT(Rf_allocMatrix(REALSXP, b_len, n_links));
    ++nprot;
    if (dimnames_sexp != R_NilValue) Rf_setAttrib(mats[j], R_DimNamesSymbol, dimnames_sexp);
  }
  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 3));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, n_out + 1));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, n_out + 1));
  nprot += 3;
  SEXP beta_vec = R_NilValue, dot_call = R_NilValue;
  if (p > 1 && k > 0) {
    beta_vec = PROTECT(Rf_allocVector(REALSXP, p));
    dot_call = PROTECT(Rf_lang3(Rf_install("%*%"), obs_z_sexp, beta_vec));
    nprot += 2;
  }
  double *out[5] = {NULL, NULL, NULL, NULL, NULL};
  for (int j = 0; j < n_out; ++j) {
    out[j] = REAL(mats[j]);
    if (gmm_buffers_overlap(out[j], out_len, py, ny) || gmm_buffers_overlap(out[j], out_len, ow, k) ||
        gmm_buffers_overlap(out[j], out_len, plp, (R_xlen_t) n_par * n_links)) {
      gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_ipw_if_chunk`.", "cdid_gmm_bound_error");
    }
    for (R_xlen_t idx = 0; idx < out_len; ++idx) out[j][idx] = 0.0;
  }
  int *st = INTEGER(status_sexp);
  st[0] = GMM_IPW_IF_OK;
  st[1] = 0;
  st[2] = 0;

  /* 1. Complete-date posting fact for every link and observer before any read. */
  for (int ell = 0; ell < n_links && st[0] == GMM_IPW_IF_OK; ++ell) {
    R_CheckUserInterrupt();
    for (int side = 0; side < 2 && st[0] == GMM_IPW_IF_OK; ++side) {
      int64_t base = (int64_t) pds[pld[2 * (R_xlen_t) ell + side] - 1] - 1;
      for (R_xlen_t i = 0; i < k; ++i) {
        if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
        int64_t a = (int64_t) a_start - 1 + opos[i];
        int64_t slot = base + a;
        if (pds[pld[2 * (R_xlen_t) ell + side] - 1] == NA_INTEGER || slot < 1 || slot > (int64_t) n_post ||
            (int64_t) ppa[slot - 1] != a || ppo[slot - 1] == NA_INTEGER || ppo[slot - 1] < 1 ||
            (R_xlen_t) ppo[slot - 1] > ny) {
          st[0] = GMM_IPW_IF_POSTING;
          st[1] = ell + 1;
          st[2] = opos[i];
          break;
        }
      }
    }
  }

  /* 2. Links in order: endpoints, raw differences, four-term IF, IF finiteness. */
  for (int ell = 0; ell < n_links && st[0] == GMM_IPW_IF_OK; ++ell) {
    R_CheckUserInterrupt();
    const double *par = plp + (R_xlen_t) n_par * ell;
    int64_t base_s = (int64_t) pds[pld[2 * (R_xlen_t) ell] - 1] - 1;
    int64_t base_t = (int64_t) pds[pld[2 * (R_xlen_t) ell + 1] - 1] - 1;
    double *col = out[0] + (R_xlen_t) b_len * ell;
    for (R_xlen_t i = 0; i < k; ++i) {
      if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
      int64_t a = (int64_t) a_start - 1 + opos[i];
      double y_s = py[ppo[base_s + a - 1] - 1];
      double y_t = py[ppo[base_t + a - 1] - 1];
      if (!R_FINITE(y_s) || !R_FINITE(y_t)) {
        st[0] = GMM_IPW_IF_ENDPOINT;
        st[1] = ell + 1;
        st[2] = opos[i];
        break;
      }
    }
    if (st[0] != GMM_IPW_IF_OK) break;
    for (R_xlen_t i = 0; i < k; ++i) {
      if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
      int64_t a = (int64_t) a_start - 1 + opos[i];
      double dy = gmm_fp_sub(py[ppo[base_t + a - 1] - 1], py[ppo[base_s + a - 1] - 1]);
      col[opos[i] - 1] = dy;
      if (!R_FINITE(dy)) {
        st[0] = GMM_IPW_IF_DIFFERENCE;
        st[1] = ell + 1;
        st[2] = opos[i];
        break;
      }
    }
    if (st[0] != GMM_IPW_IF_OK) break;
    double treated_mean = par[0], control_mean = par[1], mean_a = par[2], mean_b = par[3];
    const double *loading_beta = par + 4;
    const double *rate_coef_s1 = par + 4 + p;
    const double *z_load_all = NULL;
    if (p > 1 && k > 0) {
      /* Match generic drop(Z_k %*% loading_beta), including R's selected
         primitive/BLAS reduction and overflow. Inputs remain borrowed/read-only.
         No result is kept across links or after this native call. */
      for (int j = 0; j < p; ++j) REAL(beta_vec)[j] = loading_beta[j];
      SEXP dot_result = PROTECT(Rf_eval(dot_call, R_BaseEnv));
      z_load_all = REAL(dot_result);
    }
    for (R_xlen_t i = 0; i < k; ++i) {
      if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
      R_xlen_t r = (R_xlen_t) opos[i] - 1;
      int g = og[i];
      double w = ow[i];
      double dy = col[r];
      double rate = gmm_fp_mul(w, rate_coef_s1[g - 1]);
      /* Keep every p = 1 operation unchanged. General p uses exactly the
         matrix product evaluated above, rather than predicting BLAS order. */
      double z_load = p == 1 ? gmm_fp_add(0.0, gmm_fp_mul(loading_beta[0], oz[i])) : z_load_all[i];
      double score = gmm_fp_mul(osb[i], z_load);
      double treated = 0.0;
      double control = 0.0;
      if (g == 1) {
        treated = gmm_fp_mul(gmm_fp_div(gmm_fp_div(w, orr[i]), mean_a), gmm_fp_sub(dy, treated_mean));
      } else {
        double weight = gmm_fp_div(gmm_fp_mul(w, oodds[i]), orr[i]);
        control = gmm_fp_mul(gmm_fp_div(gmm_fp_neg(weight), mean_b), gmm_fp_sub(dy, control_mean));
      }
      double phi = gmm_fp_add(gmm_fp_add(gmm_fp_add(treated, control), score), rate);
      col[r] = phi;
      if (with_components) {
        R_xlen_t idx = (R_xlen_t) b_len * ell + r;
        out[1][idx] = rate;
        out[2][idx] = score;
        out[3][idx] = treated;
        out[4][idx] = control;
      }
      if (!R_FINITE(phi)) {
        st[0] = GMM_IPW_IF_NONFINITE;
        st[1] = ell + 1;
        st[2] = opos[i];
        break;
      }
    }
    if (p > 1 && k > 0) UNPROTECT(1); /* this link's matrix-product result */
  }

  static const char *out_names[6] = {"phi", "rate", "beta_score", "treated", "control", "status"};
  for (int j = 0; j < n_out; ++j) {
    SET_VECTOR_ELT(out_sexp, j, mats[j]);
    SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  }
  SET_VECTOR_ELT(out_sexp, n_out, status_sexp);
  SET_STRING_ELT(names_sexp, n_out, Rf_mkChar(out_names[5]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* CP08 stage 4: native complete-support batches of the eligible baseline OR/DR
   route (R/chained_or_baseline.R, R/chained_dr_backend.R). R selects the route
   from the fitted complete-support facts, preflights every allocation against
   the unchanged caller caps and passes the validated sealed owners' arrays and
   compact descriptors. Each floating-point operation is one explicit volatile
   operation in the order of the R expression it replaces. R's matrix products
   are emulated in the reference BLAS order (zero start, left-to-right terms);
   a product that would accumulate more than one term with two nonzero factors
   is refused (CP08_PRODUCT), the only case in which an optimized or fused BLAS
   could round it differently. R's sum() is replicated with its long double
   accumulator when R uses one. Posting slots are checked before any endpoint
   is read. Data-dependent failures return (code, stage, index, row) in the
   original route's order so R raises that route's exact condition; malformed
   arguments abort here. Borrowed inputs are never written. */
#define CP08_OK 0
#define CP08_SLOT 1
#define CP08_ENDPOINT 2
#define CP08_DIFFERENCE 3
#define CP08_INTERCEPT 4
#define CP08_OR_LINK_IF 5
#define CP08_DR_FRAME_IF 7
#define CP08_DR_LINK_IF 8
#define CP08_PRODUCT 9
#define CP08_GATHER 10
#define CP08_SLOPE_RHS 11
#define CP08_SLOPE_GRAM 12
#define CP08_DR_TRT_MEANS 14
#define CP08_DR_TRT_PARTIALS 15
/* Probability-clipping epsilon (did 2.5.1 / DRDID 1.3.0), the value of the R
   constant ipw_clip_eps: a support rate r <= eps enters weights as max(r, eps)
   and has no rate derivative. */
#define CP08_CLIP_EPS 1e-6

static R_xlen_t cp08_len(int nprotect, SEXP x, int type, const char *fn, const char *name) {
  if (x == R_NilValue || Rf_isObject(x) || TYPEOF(x) != type) {
    char buf[256];
    snprintf(buf, sizeof(buf), "`%s` in `%s` must be a non-object %s vector.", name, fn,
             type == INTSXP ? "integer" : (type == REALSXP ? "double" : "logical"));
    gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
  }
  return XLENGTH(x);
}

static void cp08_dims(int nprotect, SEXP x, int type, const char *fn, const char *name,
                      int nr_req, int nc_req, int *nr, int *nc) {
  R_xlen_t len = cp08_len(nprotect, x, type, fn, name);
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  int ok = dim != R_NilValue && TYPEOF(dim) == INTSXP && XLENGTH(dim) == 2 &&
    INTEGER(dim)[0] != NA_INTEGER && INTEGER(dim)[1] != NA_INTEGER &&
    INTEGER(dim)[0] >= 1 && INTEGER(dim)[1] >= 1 &&
    (R_xlen_t) INTEGER(dim)[0] * (R_xlen_t) INTEGER(dim)[1] == len &&
    (nr_req < 1 || INTEGER(dim)[0] == nr_req) && (nc_req < 1 || INTEGER(dim)[1] == nc_req);
  if (!ok) {
    char buf[256];
    snprintf(buf, sizeof(buf), "`%s` in `%s` must be a matrix of the required positive dimensions.", name, fn);
    gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
  }
  *nr = INTEGER(dim)[0];
  *nc = INTEGER(dim)[1];
}

static void cp08_malformed(int nprotect, const char *fn, const char *what) {
  char buf[320];
  snprintf(buf, sizeof(buf), "Malformed arguments to `%s`: %s.", fn, what);
  gmm_c_abort(nprotect, buf, "cdid_gmm_bound_error");
}

static void cp08_range(int nprotect, const int *x, R_xlen_t n, int lo, int hi, const char *fn, const char *what) {
  for (R_xlen_t i = 0; i < n; ++i) {
    if (x[i] == NA_INTEGER || x[i] < lo || x[i] > hi) cp08_malformed(nprotect, fn, what);
  }
}

static void cp08_dimnames(int nprotect, SEXP dimnames_sexp, int n_cols, const char *fn) {
  if (dimnames_sexp != R_NilValue &&
      (TYPEOF(dimnames_sexp) != VECSXP || XLENGTH(dimnames_sexp) != 2 ||
       (VECTOR_ELT(dimnames_sexp, 0) != R_NilValue && TYPEOF(VECTOR_ELT(dimnames_sexp, 0)) != VECSXP) ||
       TYPEOF(VECTOR_ELT(dimnames_sexp, 1)) != STRSXP ||
       XLENGTH(VECTOR_ELT(dimnames_sexp, 1)) != (R_xlen_t) n_cols)) {
    cp08_malformed(nprotect, fn, "dimnames must be NULL or list(NULL_or_slope, <one link ID per column>)");
  }
}

static void cp08_set_dimnames(SEXP phi_sexp, SEXP dimnames_sexp) {
  if (dimnames_sexp == R_NilValue) return;
  if (VECTOR_ELT(dimnames_sexp, 0) == R_NilValue) {
    Rf_setAttrib(phi_sexp, R_DimNamesSymbol, dimnames_sexp);
    return;
  }
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, R_NilValue);
  SET_VECTOR_ELT(dn, 1, VECTOR_ELT(dimnames_sexp, 1));
  Rf_setAttrib(phi_sexp, R_DimNamesSymbol, dn);
  UNPROTECT(1);
}

/* R's sum() of doubles: accumulate in R's LDOUBLE, then saturate and round. */
static double cp08_ld_value(long double s) {
  if (s > DBL_MAX) return R_PosInf;
  if (s < -DBL_MAX) return R_NegInf;
  return (double) s;
}

/* Bounded row chunk of the sealed reader and panel owner plus its gathered
   endpoints: gather column c holds calendar date dates[c], each (row, date)
   posting fact is decided once and each endpoint is read at most once. */
typedef struct {
  const double *y;
  R_xlen_t ny;
  const int *date_start;
  R_xlen_t n_dates;
  const int *post_a;
  const int *post_obs;
  R_xlen_t n_post;
  const int *rows;
  int B;
  R_xlen_t n_frame;
  const int *ci;
  const int *a;
  const double *w;
  const int *dates;
  int nd;
  double *Y;            /* B x nd endpoints, NaN until read */
  unsigned char *state; /* 0 unknown, 1 slot holds a, 2 slot fails, 3 finite read, 4 nonfinite read */
  double slot_checks;
  double reads;
} cp08_chunk_t;

static void cp08_chunk_args(int nprotect, const char *fn, cp08_chunk_t *g, SEXP rows_sexp, SEXP ci_sexp,
                            SEXP a_sexp, SEXP w_sexp) {
  R_xlen_t B = cp08_len(nprotect, rows_sexp, INTSXP, fn, "rows");
  R_xlen_t n_frame = cp08_len(nprotect, ci_sexp, INTSXP, fn, "cohort_index");
  if (cp08_len(nprotect, a_sexp, INTSXP, fn, "a_index") != n_frame ||
      cp08_len(nprotect, w_sexp, REALSXP, fn, "w") != n_frame) {
    cp08_malformed(nprotect, fn, "reader cohort, A-index and weight vectors differ in length");
  }
  if (B < 1 || B > INT_MAX || n_frame < 1 || n_frame > INT_MAX) {
    cp08_malformed(nprotect, fn, "empty or oversized row chunk or frame");
  }
  g->rows = INTEGER(rows_sexp);
  g->B = (int) B;
  g->n_frame = n_frame;
  cp08_range(nprotect, g->rows, B, 1, (int) n_frame, fn, "chunk rows outside the frame");
  g->ci = INTEGER(ci_sexp);
  g->a = INTEGER(a_sexp);
  g->w = REAL(w_sexp);
  g->y = NULL;
  g->ny = 0;
  g->date_start = NULL;
  g->n_dates = 0;
  g->post_a = NULL;
  g->post_obs = NULL;
  g->n_post = 0;
  g->dates = NULL;
  g->nd = 0;
  g->Y = NULL;
  g->state = NULL;
  g->slot_checks = 0.0;
  g->reads = 0.0;
}

static void cp08_panel_args(int nprotect, const char *fn, cp08_chunk_t *g, SEXP y_sexp, SEXP date_start_sexp,
                            SEXP post_a_sexp, SEXP post_obs_sexp, SEXP dates_sexp) {
  g->ny = cp08_len(nprotect, y_sexp, REALSXP, fn, "y");
  g->n_dates = cp08_len(nprotect, date_start_sexp, INTSXP, fn, "date_start");
  g->n_post = cp08_len(nprotect, post_a_sexp, INTSXP, fn, "post_a_index");
  if (cp08_len(nprotect, post_obs_sexp, INTSXP, fn, "post_obs_index") != g->n_post) {
    cp08_malformed(nprotect, fn, "posting A-unit and observation vectors differ in length");
  }
  R_xlen_t nd = cp08_len(nprotect, dates_sexp, INTSXP, fn, "dates");
  if (nd < 1 || nd > INT_MAX || g->n_dates < 1 || g->n_dates > INT_MAX) {
    cp08_malformed(nprotect, fn, "empty or oversized gather dates or calendar");
  }
  g->dates = INTEGER(dates_sexp);
  g->nd = (int) nd;
  cp08_range(nprotect, g->dates, nd, 1, (int) g->n_dates, fn, "gather dates outside the calendar");
  g->y = REAL(y_sexp);
  g->date_start = INTEGER(date_start_sexp);
  g->post_a = INTEGER(post_a_sexp);
  g->post_obs = INTEGER(post_obs_sexp);
}

/* Posting fact of chunk row i at gather column c: slot date_start - 1 + a
   holds A-unit a >= 1 (the check of or_baseline_complete_pair_rows). */
static int cp08_slot(cp08_chunk_t *g, int i, int c) {
  unsigned char *st = g->state + (R_xlen_t) c * g->B + i;
  if (*st == 0) {
    int a = g->a[g->rows[i] - 1];
    int ds = g->date_start[g->dates[c] - 1];
    int ok = 0;
    if (a != NA_INTEGER && a >= 1 && ds != NA_INTEGER) {
      int64_t slot = (int64_t) ds - 1 + (int64_t) a;
      ok = slot >= 1 && slot <= (int64_t) g->n_post && g->post_a[slot - 1] == a;
    }
    *st = (unsigned char) (ok ? 1 : 2);
    g->slot_checks += 1.0;
  }
  return *st != 2;
}

/* Endpoint of a checked slot (ipw_owner_read_y_positions): finite or not. */
static int cp08_read(cp08_chunk_t *g, int i, int c) {
  R_xlen_t idx = (R_xlen_t) c * g->B + i;
  unsigned char *st = g->state + idx;
  if (*st == 1) {
    int a = g->a[g->rows[i] - 1];
    int64_t slot = (int64_t) g->date_start[g->dates[c] - 1] - 1 + (int64_t) a;
    int o = g->post_obs[slot - 1];
    double v = R_NaN;
    if (o != NA_INTEGER && o >= 1 && (R_xlen_t) o <= g->ny) v = g->y[o - 1];
    g->reads += 1.0;
    if (R_FINITE(v)) {
      g->Y[idx] = v;
      *st = 3;
    } else {
      *st = 4;
    }
  }
  return *st == 3;
}

/* Raw right-minus-left outcome difference y_t - y_s of chunk row i. */
static inline double cp08_dy(const cp08_chunk_t *g, int i, int cs, int ct) {
  return gmm_fp_sub(g->Y[(R_xlen_t) ct * g->B + i], g->Y[(R_xlen_t) cs * g->B + i]);
}

static inline int cp08_member(const cp08_chunk_t *g, int i, const int *ks, int nk) {
  int c = g->ci[g->rows[i] - 1];
  for (int j = 0; j < nk; ++j) {
    if (c == ks[j]) return 1;
  }
  return 0;
}

/* One pair-row query of the original route on the chunk rows whose cohort is
   in ks: every posting slot at both dates, then (read) every endpoint, then
   every raw difference. Returns the failing row (code set) or -1. */
static int cp08_pair_rows(cp08_chunk_t *g, const int *ks, int nk, int cs, int ct, int read, int *code) {
  for (int i = 0; i < g->B; ++i) {
    if (!cp08_member(g, i, ks, nk)) continue;
    if (!cp08_slot(g, i, cs) || !cp08_slot(g, i, ct)) {
      *code = CP08_SLOT;
      return i;
    }
  }
  if (!read) return -1;
  for (int i = 0; i < g->B; ++i) {
    if (!cp08_member(g, i, ks, nk)) continue;
    if (!cp08_read(g, i, cs) || !cp08_read(g, i, ct)) {
      *code = CP08_ENDPOINT;
      return i;
    }
  }
  for (int i = 0; i < g->B; ++i) {
    if (!cp08_member(g, i, ks, nk)) continue;
    if (!R_FINITE(cp08_dy(g, i, cs, ct))) {
      *code = CP08_DIFFERENCE;
      return i;
    }
  }
  return -1;
}

/* Rate IF of one support on chunk row i: w * (S - A * rate) / massA on its
   cohort rows (S = 1 once its slots hold), zero elsewhere. */
static inline double cp08_up(const cp08_chunk_t *g, int i, int k, double rate, double mass_a) {
  int r = g->rows[i] - 1;
  if (g->ci[r] != k) return 0.0;
  double A = g->a[r] > 0 ? 1.0 : 0.0;
  return gmm_fp_div(gmm_fp_mul(g->w[r], gmm_fp_sub(1.0, gmm_fp_mul(A, rate))), mass_a);
}

/* Reference-order accumulation of one product term: z + c * u. Counts terms
   whose two factors are both nonzero. */
static inline double cp08_term(double z, double c, double u, int *nonzero) {
  if (c != 0.0 && u != 0.0) ++(*nonzero);
  return gmm_fp_add(z, gmm_fp_mul(c, u));
}

static void cp08_status(int *st, int code, int stage, int index, int row) {
  st[0] = code;
  st[1] = stage;
  st[2] = index;
  st[3] = row;
}

/* Fit-time support preparation: one cohort-row chunk (positions start..end of
   cohort_order) of one (pair, cohort) support, as or_baseline_visit with the
   complete-date reader and or_pooled_support_add at p = 0: slots of both pair
   dates, endpoints and raw differences, the declared intercept column, then
   R's sums of w[A], w[obs] and w[obs] * dy and the augmented rows
   [1 * sqrt(w), dy * sqrt(w)] that or_wls_summary compresses. */
SEXP C_gmm_cp08_support_chunk(SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
                              SEXP order_sexp, SEXP a_sexp, SEXP w_sexp, SEXP z_sexp, SEXP ctl_sexp,
                              SEXP dates_sexp, SEXP long_double_sexp, SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_support_chunk";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  R_xlen_t ny = cp08_len(0, y_sexp, REALSXP, fn, "y");
  R_xlen_t n_dates = cp08_len(0, date_start_sexp, INTSXP, fn, "date_start");
  R_xlen_t n_post = cp08_len(0, post_a_sexp, INTSXP, fn, "post_a_index");
  if (cp08_len(0, post_obs_sexp, INTSXP, fn, "post_obs_index") != n_post) {
    cp08_malformed(0, fn, "posting A-unit and observation vectors differ in length");
  }
  R_xlen_t n_order = cp08_len(0, order_sexp, INTSXP, fn, "cohort_order");
  R_xlen_t n_frame = cp08_len(0, a_sexp, INTSXP, fn, "a_index");
  if (cp08_len(0, w_sexp, REALSXP, fn, "w") != n_frame || n_frame < 1 || n_frame > INT_MAX) {
    cp08_malformed(0, fn, "weights and A-indices differ in length or are empty");
  }
  int z_nr = 0, z_nc = 0;
  cp08_dims(0, z_sexp, REALSXP, fn, "Z_a", 0, 0, &z_nr, &z_nc);
  if (cp08_len(0, ctl_sexp, INTSXP, fn, "ctl") != 3 || cp08_len(0, dates_sexp, INTSXP, fn, "dates") != 2 ||
      cp08_len(0, long_double_sexp, LGLSXP, fn, "long_double") != 1 ||
      LOGICAL(long_double_sexp)[0] == NA_LOGICAL) {
    cp08_malformed(0, fn, "ctl must be (start, end, intercept column), dates two calendar indices, long_double one flag");
  }
  const int *ctl = INTEGER(ctl_sexp);
  const int *dd = INTEGER(dates_sexp);
  int start = ctl[0], end = ctl[1], zcol = ctl[2];
  if (start == NA_INTEGER || end == NA_INTEGER || start < 1 || end < start || (R_xlen_t) end > n_order ||
      zcol == NA_INTEGER || zcol < 1 || zcol > z_nc) {
    cp08_malformed(0, fn, "chunk positions or intercept column outside their arrays");
  }
  for (int s = 0; s < 2; ++s) {
    if (dd[s] == NA_INTEGER || dd[s] < 1 || (R_xlen_t) dd[s] > n_dates) {
      cp08_malformed(0, fn, "pair dates outside the calendar");
    }
  }
  int n = end - start + 1;
  (void) gmm_checked_mul_u64(0, (uint64_t) n, 2ULL, max_entries, "C_gmm_cp08_support_chunk augmented rows (n x 2)");
  const int *ord = INTEGER(order_sexp) + (start - 1);
  cp08_range(0, ord, n, 1, (int) n_frame, fn, "cohort rows outside the frame");

  const double *py = REAL(y_sexp);
  const int *pds = INTEGER(date_start_sexp);
  const int *ppa = INTEGER(post_a_sexp);
  const int *ppo = INTEGER(post_obs_sexp);
  const int *pa = INTEGER(a_sexp);
  const double *pw = REAL(w_sexp);
  const double *pz = REAL(z_sexp) + (R_xlen_t) (zcol - 1) * z_nr;
  int ld = LOGICAL(long_double_sexp)[0] != 0;

  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP sums_sexp = PROTECT(Rf_allocVector(REALSXP, 3));
  SEXP count_sexp = PROTECT(Rf_allocVector(INTSXP, 1));
  SEXP aug_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, 2));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 2));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 5));
  const int nprot = 7;
  int *st = INTEGER(status_sexp);
  double *sums = REAL(sums_sexp);
  double *aug = REAL(aug_sexp);
  double *counts = REAL(counts_sexp);
  if (gmm_buffers_overlap(aug, 2 * (R_xlen_t) n, py, ny) || gmm_buffers_overlap(aug, 2 * (R_xlen_t) n, pw, n_frame) ||
      gmm_buffers_overlap(aug, 2 * (R_xlen_t) n, pz, z_nr)) {
    gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_cp08_support_chunk`.", "cdid_gmm_bound_error");
  }
  for (R_xlen_t idx = 0; idx < 2 * (R_xlen_t) n; ++idx) aug[idx] = 0.0;
  cp08_status(st, CP08_OK, 0, 0, 0);
  sums[0] = sums[1] = sums[2] = 0.0;
  INTEGER(count_sexp)[0] = 0;
  counts[0] = counts[1] = 0.0;

  /* 1. Both posting slots of every chunk row before any endpoint is read. */
  for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    int a = pa[ord[i] - 1];
    for (int s = 0; s < 2; ++s) {
      int ok = 0;
      if (a != NA_INTEGER && a >= 1 && pds[dd[s] - 1] != NA_INTEGER) {
        int64_t slot = (int64_t) pds[dd[s] - 1] - 1 + (int64_t) a;
        ok = slot >= 1 && slot <= (int64_t) n_post && ppa[slot - 1] == a;
      }
      counts[0] += 1.0;
      if (!ok) {
        cp08_status(st, CP08_SLOT, 1, 1, i + 1);
        break;
      }
    }
  }
  /* 2. Endpoints, then 3. raw right-minus-left differences (held in column 2). */
  for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    int a = pa[ord[i] - 1];
    double v[2];
    for (int s = 0; s < 2; ++s) {
      int o = ppo[(int64_t) pds[dd[s] - 1] - 1 + (int64_t) a - 1];
      v[s] = (o != NA_INTEGER && o >= 1 && (R_xlen_t) o <= ny) ? py[o - 1] : R_NaN;
      counts[1] += 1.0;
    }
    if (!R_FINITE(v[0]) || !R_FINITE(v[1])) {
      cp08_status(st, CP08_ENDPOINT, 1, 1, i + 1);
      break;
    }
    aug[n + i] = gmm_fp_sub(v[1], v[0]);
  }
  for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
    if (!R_FINITE(aug[n + i])) cp08_status(st, CP08_DIFFERENCE, 1, 1, i + 1);
  }
  /* 4. The declared intercept column equals one on every observed pair row. */
  for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
    int a = pa[ord[i] - 1];
    if (a > z_nr) {
      cp08_malformed(nprot, fn, "A-unit outside the design rows");
    }
    if (!(pz[a - 1] == 1.0)) cp08_status(st, CP08_INTERCEPT, 1, 1, i + 1);
  }
  /* 5. R's sums and the augmented rows. */
  if (st[0] == CP08_OK) {
    long double lA = 0.0L, lW = 0.0L, lWD = 0.0L;
    double dA = 0.0, dW = 0.0, dWD = 0.0;
    int n_a = 0;
    for (int i = 0; i < n; ++i) {
      int r = ord[i] - 1;
      double wi = pw[r];
      double dy = aug[n + i];
      if (pa[r] > 0) {
        if (ld) lA += (long double) wi; else dA = gmm_fp_add(dA, wi);
        ++n_a;
      }
      if (ld) lW += (long double) wi; else dW = gmm_fp_add(dW, wi);
      double prod = gmm_fp_mul(wi, dy);
      if (ld) lWD += (long double) prod; else dWD = gmm_fp_add(dWD, prod);
      double s = gmm_fp_sqrt(wi);
      aug[i] = gmm_fp_mul(1.0, s);
      aug[n + i] = gmm_fp_mul(dy, s);
    }
    sums[0] = ld ? cp08_ld_value(lA) : dA;
    sums[1] = ld ? cp08_ld_value(lW) : dW;
    sums[2] = ld ? cp08_ld_value(lWD) : dWD;
    INTEGER(count_sexp)[0] = n_a;
  }

  static const char *out_names[5] = {"status", "sums", "n_A", "augmented", "counts"};
  SET_VECTOR_ELT(out_sexp, 0, status_sexp);
  SET_VECTOR_ELT(out_sexp, 1, sums_sexp);
  SET_VECTOR_ELT(out_sexp, 2, count_sexp);
  SET_VECTOR_ELT(out_sexp, 3, aug_sexp);
  SET_VECTOR_ELT(out_sexp, 4, counts_sexp);
  for (int j = 0; j < 5; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* Fit-time treatment-model reads of the eligible DR route: the raw differences
   y_t - y_s of one model's frame rows (increasing), checked as the original
   scan's or_baseline_complete_pair_rows calls on its frame-row chunks
   ((row - 1) / chunk): every slot of a chunk, then its endpoints, then its
   differences, chunk by chunk. */
SEXP C_gmm_cp08_pair_rows(SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
                          SEXP rows_sexp, SEXP a_sexp, SEXP dates_sexp, SEXP chunk_sexp, SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_pair_rows";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  R_xlen_t ny = cp08_len(0, y_sexp, REALSXP, fn, "y");
  R_xlen_t n_dates = cp08_len(0, date_start_sexp, INTSXP, fn, "date_start");
  R_xlen_t n_post = cp08_len(0, post_a_sexp, INTSXP, fn, "post_a_index");
  if (cp08_len(0, post_obs_sexp, INTSXP, fn, "post_obs_index") != n_post) {
    cp08_malformed(0, fn, "posting A-unit and observation vectors differ in length");
  }
  R_xlen_t n = cp08_len(0, rows_sexp, INTSXP, fn, "rows");
  R_xlen_t n_frame = cp08_len(0, a_sexp, INTSXP, fn, "a_index");
  if (n_frame < 1 || n_frame > INT_MAX || n > n_frame || cp08_len(0, dates_sexp, INTSXP, fn, "dates") != 2 ||
      cp08_len(0, chunk_sexp, INTSXP, fn, "chunk") != 1) {
    cp08_malformed(0, fn, "rows exceed the frame, or dates/chunk have the wrong length");
  }
  const int *rows = INTEGER(rows_sexp);
  const int *dd = INTEGER(dates_sexp);
  int chunk = INTEGER(chunk_sexp)[0];
  cp08_range(0, rows, n, 1, (int) n_frame, fn, "rows outside the frame");
  for (R_xlen_t i = 1; i < n; ++i) {
    if (rows[i] <= rows[i - 1]) cp08_malformed(0, fn, "rows must strictly increase");
  }
  cp08_range(0, dd, 2, 1, (int) n_dates, fn, "pair dates outside the calendar");
  if (chunk == NA_INTEGER || chunk < 1) cp08_malformed(0, fn, "chunk must be a positive count");
  (void) gmm_checked_mul_u64(0, (uint64_t) (n > 0 ? n : 1), 1ULL, max_entries, "C_gmm_cp08_pair_rows differences");
  const double *py = REAL(y_sexp);
  const int *pds = INTEGER(date_start_sexp);
  const int *ppa = INTEGER(post_a_sexp);
  const int *ppo = INTEGER(post_obs_sexp);
  const int *pa = INTEGER(a_sexp);

  SEXP dy_sexp = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 2));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  const int nprot = 5;
  double *dy = REAL(dy_sexp);
  double *counts = REAL(counts_sexp);
  int *st = INTEGER(status_sexp);
  if (gmm_buffers_overlap(dy, n, py, ny)) {
    gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_cp08_pair_rows`.", "cdid_gmm_bound_error");
  }
  for (R_xlen_t i = 0; i < n; ++i) dy[i] = 0.0;
  counts[0] = counts[1] = 0.0;
  cp08_status(st, CP08_OK, 0, 0, 0);
  R_xlen_t lo = 0;
  while (lo < n && st[0] == CP08_OK) {
    R_CheckUserInterrupt();
    int c = (rows[lo] - 1) / chunk;
    R_xlen_t hi = lo;
    while (hi < n && (rows[hi] - 1) / chunk == c) ++hi;
    for (R_xlen_t i = lo; i < hi && st[0] == CP08_OK; ++i) {
      int a = pa[rows[i] - 1];
      for (int s = 0; s < 2; ++s) {
        int ok = 0;
        if (a != NA_INTEGER && a >= 1 && pds[dd[s] - 1] != NA_INTEGER) {
          int64_t slot = (int64_t) pds[dd[s] - 1] - 1 + (int64_t) a;
          ok = slot >= 1 && slot <= (int64_t) n_post && ppa[slot - 1] == a;
        }
        counts[0] += 1.0;
        if (!ok) {
          cp08_status(st, CP08_SLOT, 2, c + 1, (int) i + 1);
          break;
        }
      }
    }
    for (R_xlen_t i = lo; i < hi && st[0] == CP08_OK; ++i) {
      int a = pa[rows[i] - 1];
      double v[2];
      for (int s = 0; s < 2; ++s) {
        int o = ppo[(int64_t) pds[dd[s] - 1] - 1 + (int64_t) a - 1];
        v[s] = (o != NA_INTEGER && o >= 1 && (R_xlen_t) o <= ny) ? py[o - 1] : R_NaN;
        counts[1] += 1.0;
      }
      if (!R_FINITE(v[0]) || !R_FINITE(v[1])) {
        cp08_status(st, CP08_ENDPOINT, 2, c + 1, (int) i + 1);
        break;
      }
      dy[i] = gmm_fp_sub(v[1], v[0]);
    }
    for (R_xlen_t i = lo; i < hi && st[0] == CP08_OK; ++i) {
      if (!R_FINITE(dy[i])) cp08_status(st, CP08_DIFFERENCE, 2, c + 1, (int) i + 1);
    }
    lo = hi;
  }
  static const char *out_names[3] = {"dy", "status", "counts"};
  SET_VECTOR_ELT(out_sexp, 0, dy_sexp);
  SET_VECTOR_ELT(out_sexp, 1, status_sexp);
  SET_VECTOR_ELT(out_sexp, 2, counts_sexp);
  for (int j = 0; j < 3; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* The scalar search of ipw_owner_match_sorted_double over one posting
   subrange [lo, hi] (1-based positions): mid = lo + floor((hi - lo) / 2), stop
   at the first equal value, descend by the comparison of that value with the
   target. Duplicate or unsorted postings therefore yield the position of the
   same path. An NA posting value is where the R comparison would fail; it is
   refused as malformed instead of being ordered. Returns 0 when absent. */
static int cp08_sparse_search(int nprotect, const char *fn, SEXP pa, int64_t lo, int64_t hi, int target) {
  while (lo <= hi) {
    int64_t mid = lo + (hi - lo) / 2;
    int value = INTEGER_ELT(pa, (R_xlen_t) mid - 1);
    if (value == NA_INTEGER) cp08_malformed(nprotect, fn, "a visited posting A-index is NA");
    if (value == target) return (int) mid;
    if (value < target) lo = mid + 1; else hi = mid - 1;
  }
  return 0;
}

/* General-route pair-row join of cp04_pair_rows on the owner's sparse calendar
   postings (missing panels; no complete-date slot arithmetic). For every query
   A-index a > 0 (0 marks a frame row off A) the scalar search above runs over
   the s-date posting subrange bounds[0..1] and, on a hit, over the t-date
   subrange bounds[2..3]. Observed local rows (1-based query positions) are
   returned in query order; when read is TRUE their two posting observation
   positions follow, which R hands to ipw_owner_read_y_positions exactly as
   before. The postings are borrowed and only visited positions are read; the
   only scratch is three ints per query row, reserved against max_entries. */
SEXP C_gmm_cp08_sparse_pair_rows(SEXP a_sexp, SEXP post_a_sexp, SEXP post_obs_sexp, SEXP bounds_sexp,
                                 SEXP read_sexp, SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_sparse_pair_rows";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  R_xlen_t n = cp08_len(0, a_sexp, INTSXP, fn, "a_index");
  R_xlen_t n_post = cp08_len(0, post_a_sexp, INTSXP, fn, "post_a_index");
  if (cp08_len(0, post_obs_sexp, INTSXP, fn, "post_obs_index") != n_post) {
    cp08_malformed(0, fn, "posting A-unit and observation vectors differ in length");
  }
  if (n > INT_MAX || cp08_len(0, bounds_sexp, INTSXP, fn, "bounds") != 4 ||
      cp08_len(0, read_sexp, LGLSXP, fn, "read") != 1 || LOGICAL(read_sexp)[0] == NA_LOGICAL) {
    cp08_malformed(0, fn, "oversized query, or bounds/read have the wrong length or are NA");
  }
  if (n_post < 1 || n_post > INT_MAX) cp08_malformed(0, fn, "empty or oversized postings");
  const int *bd = INTEGER(bounds_sexp);
  cp08_range(0, bd, 4, 1, (int) n_post, fn, "posting bounds outside the postings");
  if (bd[0] > bd[1] || bd[2] > bd[3]) cp08_malformed(0, fn, "a posting subrange ends before it starts");
  uint64_t scratch_entries = gmm_checked_mul_u64(0, (uint64_t) (n > 0 ? n : 1), 3ULL, max_entries,
                             "C_gmm_cp08_sparse_pair_rows row scratch (n x 3)");
  int read = LOGICAL(read_sexp)[0] != 0;
  /* Element access borrows compact ALTREP postings without expanding the
     whole sparse table. The query is checked only after its allocation bound. */
  for (R_xlen_t i = 0; i < n; ++i) {
    int a = INTEGER_ELT(a_sexp, i);
    if (a == NA_INTEGER || a < 0) cp08_malformed(0, fn, "query A-indices must be nonnegative");
  }

  size_t cap = (size_t) (n > 0 ? n : 1);
  int *obs = (int *) R_alloc(cap, sizeof(int));
  int *ms = (int *) R_alloc(cap, sizeof(int));
  int *mt = (int *) R_alloc(cap, sizeof(int));
  R_xlen_t k = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    int a = INTEGER_ELT(a_sexp, i);
    if (a <= 0) continue;
    int s = cp08_sparse_search(0, fn, post_a_sexp, bd[0], bd[1], a);
    if (s == 0) continue;
    int t = cp08_sparse_search(0, fn, post_a_sexp, bd[2], bd[3], a);
    if (t == 0) continue;
    obs[k] = (int) (i + 1);
    ms[k] = s;
    mt[k] = t;
    ++k;
  }

  /* Scratch remains live while exact-length results are allocated. */
  uint64_t output_entries = gmm_checked_mul_u64(0, (uint64_t) k,
    read ? 3ULL : 1ULL, max_entries, "C_gmm_cp08_sparse_pair_rows output");
  if (R_finite(max_entries) && scratch_entries > (uint64_t) max_entries - output_entries)
    gmm_c_abort(0, "C_gmm_cp08_sparse_pair_rows live scratch plus output exceeds the explicit entry bound.",
      "cdid_gmm_bound_error");

  SEXP obs_sexp = PROTECT(Rf_allocVector(INTSXP, k));
  SEXP ps_sexp = PROTECT(read ? Rf_allocVector(INTSXP, k) : R_NilValue);
  SEXP pt_sexp = PROTECT(read ? Rf_allocVector(INTSXP, k) : R_NilValue);
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  const int nprot = 5;
  int *po = INTEGER(obs_sexp);
  for (R_xlen_t j = 0; j < k; ++j) po[j] = obs[j];
  if (read) {
    int *qs = INTEGER(ps_sexp);
    int *qt = INTEGER(pt_sexp);
    for (R_xlen_t j = 0; j < k; ++j) {
      qs[j] = INTEGER_ELT(post_obs_sexp, (R_xlen_t) ms[j] - 1);
      qt[j] = INTEGER_ELT(post_obs_sexp, (R_xlen_t) mt[j] - 1);
    }
  }
  static const char *out_names[3] = {"obs", "pos_s", "pos_t"};
  SET_VECTOR_ELT(out_sexp, 0, obs_sexp);
  SET_VECTOR_ELT(out_sexp, 1, ps_sexp);
  SET_VECTOR_ELT(out_sexp, 2, pt_sexp);
  for (int j = 0; j < 3; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* Checked addition is local to the new summary's simultaneous-entry preflight. */
static uint64_t cp08_dr_summary_add(int nprotect, uint64_t a, uint64_t b,
                                  double max_entries, const char *label) {
  if (a > UINT64_MAX - b)
    gmm_c_abort(nprotect, "DR summary allocation sum overflowed.", "cdid_gmm_bound_error");
  return gmm_checked_mul_u64(nprotect, 1ULL, a + b, max_entries, label);
}

/* Complete-roster DR treatment-summary batch for models sharing the first
   model's prepared outcome-independent invariants (share$prep and share$inputs).
   For each model t = 1..M_batch in model order:
   1. Checks both posting slots of every observer chunk (rows[i] - 1) / chunk,
      reads both endpoints and subtracts raw y_t - y_s BEFORE any weighted sum.
   2. Computes T_ipw, M_ipw, mu_ipw, residualized r_obs = dy - (alpha + L Z_or beta),
      T_dr, M_dr, mu_dr, d_bar_plus, d_bar_minus, d_bar_diff and mu_dr_decomp,
      then checks finiteness (CP08_DR_TRT_MEANS).
   3. Computes rate cross-derivatives C, upstream partials U_dr/U_ipw and beta
      partials V_dr/V_ipw directly without N x 4 or N x p_trt per-model matrices,
      then checks finiteness (CP08_DR_TRT_PARTIALS). A floored support rate
      (r <= CP08_CLIP_EPS) has zero rate scale. The R caller declines this batch
      whenever a control odds is capped, so V uses wm without an odds mask. */
SEXP C_gmm_cp08_dr_treatment_batch(SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
                                   SEXP rows_sexp, SEXP a_sexp, SEXP X_trt_sexp, SEXP Z_or_sexp,
                                   SEXP beta_sexp, SEXP wp_sexp, SEXP wm_sexp, SEXP cohort_rows_sexp,
                                   SEXP score_sums_sexp, SEXP wm_design_sexp, SEXP norm_sexp,
                                   SEXP pairs_sexp, SEXP model_par_sexp, SEXP ctl_sexp,
                                   SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_dr_treatment_batch";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  R_xlen_t ny = cp08_len(0, y_sexp, REALSXP, fn, "y");
  R_xlen_t n_dates = cp08_len(0, date_start_sexp, INTSXP, fn, "date_start");
  R_xlen_t n_post = cp08_len(0, post_a_sexp, INTSXP, fn, "post_a_index");
  if (cp08_len(0, post_obs_sexp, INTSXP, fn, "post_obs_index") != n_post) {
    cp08_malformed(0, fn, "posting A-unit and observation vectors differ in length");
  }
  R_xlen_t n_x = cp08_len(0, rows_sexp, INTSXP, fn, "rows");
  R_xlen_t n_frame = cp08_len(0, a_sexp, INTSXP, fn, "a_index");
  if (n_frame < 1 || n_frame > INT_MAX || n_x < 1 || n_x > n_frame || n_dates < 1 || n_dates > INT_MAX) {
    cp08_malformed(0, fn, "rows exceed the frame or calendar is empty");
  }
  int n = (int) n_x;
  const int *rows = INTEGER(rows_sexp);
  cp08_range(0, rows, n_x, 1, (int) n_frame, fn, "rows outside the frame");
  for (R_xlen_t i = 1; i < n_x; ++i) {
    if (rows[i] <= rows[i - 1]) cp08_malformed(0, fn, "rows must strictly increase");
  }
  int x_nr = 0, p_trt = 0;
  cp08_dims(0, X_trt_sexp, REALSXP, fn, "X_trt", n, 0, &x_nr, &p_trt);
  R_xlen_t z_len = cp08_len(0, Z_or_sexp, REALSXP, fn, "Z_or");
  SEXP z_dim = Rf_getAttrib(Z_or_sexp, R_DimSymbol);
  if (z_dim == R_NilValue || TYPEOF(z_dim) != INTSXP || XLENGTH(z_dim) != 2 ||
      INTEGER(z_dim)[0] != n || INTEGER(z_dim)[1] == NA_INTEGER || INTEGER(z_dim)[1] < 0 ||
      (R_xlen_t) n * (R_xlen_t) INTEGER(z_dim)[1] != z_len) {
    cp08_malformed(0, fn, "Z_or must be an n x p_or double matrix with p_or >= 0");
  }
  int p_or = INTEGER(z_dim)[1];
  if (cp08_len(0, beta_sexp, REALSXP, fn, "beta") != (R_xlen_t) p_or) {
    cp08_malformed(0, fn, "beta length must equal ncol(Z_or)");
  }
  const double *beta = REAL(beta_sexp);
  for (int j = 0; j < p_or; ++j) {
    if (!R_FINITE(beta[j])) cp08_malformed(0, fn, "beta must be finite");
  }
  if (cp08_len(0, wp_sexp, REALSXP, fn, "wp") != n_x ||
      cp08_len(0, wm_sexp, REALSXP, fn, "wm") != n_x) {
    cp08_malformed(0, fn, "wp and wm must have length n");
  }
  if (cohort_rows_sexp == R_NilValue || Rf_isObject(cohort_rows_sexp) ||
      TYPEOF(cohort_rows_sexp) != VECSXP || XLENGTH(cohort_rows_sexp) < 1 ||
      XLENGTH(cohort_rows_sexp) > INT_MAX) {
    cp08_malformed(0, fn, "cohort_rows must be a nonempty list of integer row vectors");
  }
  int m = (int) XLENGTH(cohort_rows_sexp);
  for (int k = 0; k < m; ++k) {
    SEXP rr_k = VECTOR_ELT(cohort_rows_sexp, k);
    R_xlen_t n_rr = cp08_len(0, rr_k, INTSXP, fn, "cohort_rows[[k]]");
    if (n_rr > n_x) cp08_malformed(0, fn, "cohort_rows[[k]] exceeds n");
    const int *rk = INTEGER(rr_k);
    cp08_range(0, rk, n_rr, 1, n, fn, "cohort_rows[[k]] outside 1..n");
    for (R_xlen_t i = 1; i < n_rr; ++i) {
      if (rk[i] <= rk[i - 1]) cp08_malformed(0, fn, "cohort_rows[[k]] must strictly increase");
    }
  }
  int ss_nr = 0, ss_nc = 0;
  cp08_dims(0, score_sums_sexp, REALSXP, fn, "score_sums", p_trt, m, &ss_nr, &ss_nc);
  if (cp08_len(0, wm_design_sexp, REALSXP, fn, "wm_design") != (R_xlen_t) p_trt) {
    cp08_malformed(0, fn, "wm_design length must equal p_trt");
  }
  if (cp08_len(0, norm_sexp, REALSXP, fn, "norm") != 5) {
    cp08_malformed(0, fn, "norm must be (sum_wp, sum_wm, dp, dm, N)");
  }
  const double *norm = REAL(norm_sexp);
  double sum_wp = norm[0], sum_wm = norm[1], dp = norm[2], dm = norm[3], N_dbl = norm[4];
  if (!R_FINITE(sum_wp) || sum_wp <= 0.0 || !R_FINITE(sum_wm) || sum_wm <= 0.0 ||
      !R_FINITE(dp) || dp <= 0.0 || !R_FINITE(dm) || dm <= 0.0 ||
      !R_FINITE(N_dbl) || N_dbl < 1.0) {
    cp08_malformed(0, fn, "normalizers and N must be positive finite");
  }
  int pr_nr = 0, M_batch = 0, mp_nr = 0, mp_nc = 0;
  cp08_dims(0, pairs_sexp, INTSXP, fn, "pairs", 2, 0, &pr_nr, &M_batch);
  cp08_dims(0, model_par_sexp, REALSXP, fn, "model_par", 2 + m, M_batch, &mp_nr, &mp_nc);
  const int *pairs = INTEGER(pairs_sexp);
  cp08_range(0, pairs, 2 * (R_xlen_t) M_batch, 1, (int) n_dates, fn, "pair dates outside the calendar");
  const double *model_par = REAL(model_par_sexp);
  for (int t = 0; t < M_batch; ++t) {
    const double *mp_t = model_par + (R_xlen_t) (2 + m) * t;
    if (!R_FINITE(mp_t[0]) || !R_FINITE(mp_t[1])) {
      cp08_malformed(0, fn, "model L and alpha must be finite");
    }
    for (int k = 0; k < m; ++k) {
      if (!R_FINITE(mp_t[2 + k]) || mp_t[2 + k] == 0.0) {
        cp08_malformed(0, fn, "model support rates must be nonzero finite");
      }
    }
  }
  if (cp08_len(0, ctl_sexp, INTSXP, fn, "ctl") != 2) {
    cp08_malformed(0, fn, "ctl must be (chunk_rows, long_double)");
  }
  const int *ctl = INTEGER(ctl_sexp);
  int chunk = ctl[0], ld = ctl[1];
  if (chunk == NA_INTEGER || chunk < 1 || ld == NA_INTEGER || (ld != 0 && ld != 1)) {
    cp08_malformed(0, fn, "chunk_rows must be positive and long_double 0 or 1");
  }

  uint64_t scratch_entries = gmm_checked_mul_u64(0, 2ULL, (uint64_t) n, max_entries,
    "C_gmm_cp08_dr_treatment_batch raw-pair scratch (2 x n)");
  uint64_t per_model_out = cp08_dr_summary_add(0, 12ULL,
    gmm_checked_mul_u64(0, 3ULL, (uint64_t) p_or, max_entries, "C_gmm_cp08_dr_treatment_batch d_bar width"),
    max_entries, "C_gmm_cp08_dr_treatment_batch means + d_bar width");
  per_model_out = cp08_dr_summary_add(0, per_model_out,
    gmm_checked_mul_u64(0, (uint64_t) p_trt, (uint64_t) m, max_entries, "C_gmm_cp08_dr_treatment_batch cross width"),
    max_entries, "C_gmm_cp08_dr_treatment_batch cross summary width");
  per_model_out = cp08_dr_summary_add(0, per_model_out,
    gmm_checked_mul_u64(0, 8ULL, (uint64_t) (m + p_trt), max_entries, "C_gmm_cp08_dr_treatment_batch partials width"),
    max_entries, "C_gmm_cp08_dr_treatment_batch per-model summary width");
  uint64_t out_entries = gmm_checked_mul_u64(0, (uint64_t) M_batch, per_model_out, max_entries,
    "C_gmm_cp08_dr_treatment_batch summary outputs");
  (void) cp08_dr_summary_add(0, scratch_entries, out_entries, max_entries,
    "C_gmm_cp08_dr_treatment_batch total native entries");

  if ((uint64_t) p_trt * (uint64_t) m > INT_MAX ||
      4ULL * (uint64_t) m > INT_MAX || 4ULL * (uint64_t) p_trt > INT_MAX)
    cp08_malformed(0, fn, "summary row dimension exceeds integer matrix dimensions");

  const double *py = REAL(y_sexp);
  const int *pds = INTEGER(date_start_sexp);
  const int *ppa = INTEGER(post_a_sexp);
  const int *ppo = INTEGER(post_obs_sexp);
  const int *pa = INTEGER(a_sexp);
  const double *px_trt = REAL(X_trt_sexp);
  const double *pz_or = REAL(Z_or_sexp);
  const double *pwp = REAL(wp_sexp);
  const double *pwm = REAL(wm_sexp);
  const double *pscore_sums = REAL(score_sums_sexp);
  const double *pwm_design = REAL(wm_design_sexp);

  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP means_sexp = PROTECT(Rf_allocMatrix(REALSXP, 12, M_batch));
  SEXP d_bar_plus_sexp = PROTECT(Rf_allocMatrix(REALSXP, p_or, M_batch));
  SEXP d_bar_minus_sexp = PROTECT(Rf_allocMatrix(REALSXP, p_or, M_batch));
  SEXP d_bar_diff_sexp = PROTECT(Rf_allocMatrix(REALSXP, p_or, M_batch));
  SEXP cross_sexp = PROTECT(Rf_allocMatrix(REALSXP, p_trt * m, M_batch));
  SEXP U_dr_sexp = PROTECT(Rf_allocMatrix(REALSXP, 4 * m, M_batch));
  SEXP U_ipw_sexp = PROTECT(Rf_allocMatrix(REALSXP, 4 * m, M_batch));
  SEXP V_dr_sexp = PROTECT(Rf_allocMatrix(REALSXP, 4 * p_trt, M_batch));
  SEXP V_ipw_sexp = PROTECT(Rf_allocMatrix(REALSXP, 4 * p_trt, M_batch));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 3));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 11));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 11));
  const int nprot = 13;

  int *st = INTEGER(status_sexp);
  double *means = REAL(means_sexp);
  double *d_bar_plus = REAL(d_bar_plus_sexp);
  double *d_bar_minus = REAL(d_bar_minus_sexp);
  double *d_bar_diff = REAL(d_bar_diff_sexp);
  double *cross_out = REAL(cross_sexp);
  double *U_dr_out = REAL(U_dr_sexp);
  double *U_ipw_out = REAL(U_ipw_sexp);
  double *V_dr_out = REAL(V_dr_sexp);
  double *V_ipw_out = REAL(V_ipw_sexp);
  double *counts = REAL(counts_sexp);
  double *dy_buf = (double *) R_alloc((size_t) n, sizeof(double));
  double *r_buf = (double *) R_alloc((size_t) n, sizeof(double));

  if (gmm_buffers_overlap(dy_buf, n_x, py, ny) ||
      gmm_buffers_overlap(r_buf, n_x, py, ny) ||
      gmm_buffers_overlap(means, 12 * (R_xlen_t) M_batch, py, ny) ||
      gmm_buffers_overlap(cross_out, (R_xlen_t) (p_trt * m) * M_batch, pscore_sums, (R_xlen_t) p_trt * m) ||
      gmm_buffers_overlap(V_dr_out, (R_xlen_t) (4 * p_trt) * M_batch, px_trt, (R_xlen_t) n * p_trt)) {
    gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_cp08_dr_treatment_batch`.", "cdid_gmm_bound_error");
  }

  cp08_status(st, CP08_OK, 0, 0, 0);
  for (R_xlen_t idx = 0; idx < 12 * (R_xlen_t) M_batch; ++idx) means[idx] = 0.0;
  for (R_xlen_t idx = 0; idx < (R_xlen_t) p_or * M_batch; ++idx) {
    d_bar_plus[idx] = 0.0;
    d_bar_minus[idx] = 0.0;
    d_bar_diff[idx] = 0.0;
  }
  for (R_xlen_t idx = 0; idx < (R_xlen_t) (p_trt * m) * M_batch; ++idx) cross_out[idx] = 0.0;
  for (R_xlen_t idx = 0; idx < (R_xlen_t) (4 * m) * M_batch; ++idx) {
    U_dr_out[idx] = 0.0;
    U_ipw_out[idx] = 0.0;
  }
  for (R_xlen_t idx = 0; idx < (R_xlen_t) (4 * p_trt) * M_batch; ++idx) {
    V_dr_out[idx] = 0.0;
    V_ipw_out[idx] = 0.0;
  }
  counts[0] = counts[1] = counts[2] = 0.0;

  double denom_v = gmm_fp_mul(N_dbl, dm);

  for (int t = 0; t < M_batch && st[0] == CP08_OK; ++t) {
    R_CheckUserInterrupt();
    int dd[2] = {pairs[2 * (R_xlen_t) t], pairs[2 * (R_xlen_t) t + 1]};
    const double *mp_t = model_par + (R_xlen_t) (2 + m) * t;
    double L_reg = mp_t[0];
    double alpha_part = mp_t[1];
    const double *rates = mp_t + 2;

    /* 1. Raw pair differences y_t - y_s in original observer chunk order. */
    R_xlen_t lo = 0;
    while (lo < n_x && st[0] == CP08_OK) {
      if (lo > 0 && (lo % 4096) == 0) R_CheckUserInterrupt();
      int c = (rows[lo] - 1) / chunk;
      R_xlen_t hi = lo;
      while (hi < n_x && (rows[hi] - 1) / chunk == c) ++hi;
      for (R_xlen_t i = lo; i < hi && st[0] == CP08_OK; ++i) {
        int a = pa[rows[i] - 1];
        for (int s = 0; s < 2; ++s) {
          int ok = 0;
          if (a != NA_INTEGER && a >= 1 && pds[dd[s] - 1] != NA_INTEGER) {
            int64_t slot = (int64_t) pds[dd[s] - 1] - 1 + (int64_t) a;
            ok = slot >= 1 && slot <= (int64_t) n_post && ppa[slot - 1] == a;
          }
          counts[0] += 1.0;
          if (!ok) {
            cp08_status(st, CP08_SLOT, t + 1, c + 1, (int) i + 1);
            break;
          }
        }
      }
      for (R_xlen_t i = lo; i < hi && st[0] == CP08_OK; ++i) {
        int a = pa[rows[i] - 1];
        double v[2];
        for (int s = 0; s < 2; ++s) {
          int o = ppo[(int64_t) pds[dd[s] - 1] - 1 + (int64_t) a - 1];
          v[s] = (o != NA_INTEGER && o >= 1 && (R_xlen_t) o <= ny) ? py[o - 1] : R_NaN;
          counts[1] += 1.0;
        }
        if (!R_FINITE(v[0]) || !R_FINITE(v[1])) {
          cp08_status(st, CP08_ENDPOINT, t + 1, c + 1, (int) i + 1);
          break;
        }
        dy_buf[i] = gmm_fp_sub(v[1], v[0]);
      }
      for (R_xlen_t i = lo; i < hi && st[0] == CP08_OK; ++i) {
        if (!R_FINITE(dy_buf[i])) {
          cp08_status(st, CP08_DIFFERENCE, t + 1, c + 1, (int) i + 1);
          break;
        }
      }
      lo = hi;
    }
    if (st[0] != CP08_OK) break;

    /* 2. IPW and residualized DR means and covariate contrasts. */
    long double l_wp_dy = 0.0L, l_wm_dy = 0.0L, l_wp_r = 0.0L, l_wm_r = 0.0L;
    double d_wp_dy = 0.0, d_wm_dy = 0.0, d_wp_r = 0.0, d_wm_r = 0.0;
    for (R_xlen_t i = 0; i < n_x; ++i) {
      double dy_i = dy_buf[i];
      double db = 0.0;
      for (int j = 0; j < p_or; ++j) {
        double d_ij = gmm_fp_mul(L_reg, pz_or[i + (R_xlen_t) j * n_x]);
        db = gmm_fp_add(db, gmm_fp_mul(d_ij, beta[j]));
      }
      double m_hat_i = gmm_fp_add(alpha_part, db);
      double r_i = gmm_fp_sub(dy_i, m_hat_i);
      r_buf[i] = r_i;
      double wp_i = pwp[i], wm_i = pwm[i];
      double p_wp_dy = gmm_fp_mul(wp_i, dy_i);
      double p_wm_dy = gmm_fp_mul(wm_i, dy_i);
      double p_wp_r = gmm_fp_mul(wp_i, r_i);
      double p_wm_r = gmm_fp_mul(wm_i, r_i);
      if (ld) {
        l_wp_dy += (long double) p_wp_dy;
        l_wm_dy += (long double) p_wm_dy;
        l_wp_r += (long double) p_wp_r;
        l_wm_r += (long double) p_wm_r;
      } else {
        d_wp_dy = gmm_fp_add(d_wp_dy, p_wp_dy);
        d_wm_dy = gmm_fp_add(d_wm_dy, p_wm_dy);
        d_wp_r = gmm_fp_add(d_wp_r, p_wp_r);
        d_wm_r = gmm_fp_add(d_wm_r, p_wm_r);
      }
    }
    double T_ipw = gmm_fp_div(ld ? cp08_ld_value(l_wp_dy) : d_wp_dy, sum_wp);
    double M_ipw = gmm_fp_div(ld ? cp08_ld_value(l_wm_dy) : d_wm_dy, sum_wm);
    double mu_ipw = gmm_fp_sub(T_ipw, M_ipw);
    double T_dr = gmm_fp_div(ld ? cp08_ld_value(l_wp_r) : d_wp_r, sum_wp);
    double M_dr = gmm_fp_div(ld ? cp08_ld_value(l_wm_r) : d_wm_r, sum_wm);
    double mu_dr = gmm_fp_sub(T_dr, M_dr);

    double *dbp = d_bar_plus + (R_xlen_t) p_or * t;
    double *dbm = d_bar_minus + (R_xlen_t) p_or * t;
    double *dbd = d_bar_diff + (R_xlen_t) p_or * t;
    long double l_decomp = 0.0L;
    double d_decomp = 0.0;
    for (int j = 0; j < p_or; ++j) {
      const double *z_col = pz_or + (R_xlen_t) j * n_x;
      long double l_p = 0.0L, l_m = 0.0L;
      double d_p = 0.0, d_m = 0.0;
      for (R_xlen_t i = 0; i < n_x; ++i) {
        double d_ij = gmm_fp_mul(L_reg, z_col[i]);
        double pp = gmm_fp_mul(pwp[i], d_ij);
        double pm = gmm_fp_mul(pwm[i], d_ij);
        if (ld) {
          l_p += (long double) pp;
          l_m += (long double) pm;
        } else {
          d_p = gmm_fp_add(d_p, pp);
          d_m = gmm_fp_add(d_m, pm);
        }
      }
      dbp[j] = gmm_fp_div(ld ? cp08_ld_value(l_p) : d_p, sum_wp);
      dbm[j] = gmm_fp_div(ld ? cp08_ld_value(l_m) : d_m, sum_wm);
      dbd[j] = gmm_fp_sub(dbp[j], dbm[j]);
      double prod_j = gmm_fp_mul(dbd[j], beta[j]);
      if (ld) l_decomp += (long double) prod_j; else d_decomp = gmm_fp_add(d_decomp, prod_j);
    }
    double decomp_sub = (p_or > 0) ? (ld ? cp08_ld_value(l_decomp) : d_decomp) : 0.0;
    double mu_dr_decomp = gmm_fp_sub(mu_ipw, decomp_sub);

    int means_ok = R_FINITE(T_ipw) && R_FINITE(M_ipw) && R_FINITE(mu_ipw) &&
                   R_FINITE(T_dr) && R_FINITE(M_dr) && R_FINITE(mu_dr) &&
                   R_FINITE(mu_dr_decomp);
    for (int j = 0; j < p_or && means_ok; ++j) {
      if (!R_FINITE(dbp[j]) || !R_FINITE(dbm[j])) means_ok = 0;
    }
    if (!means_ok) {
      cp08_status(st, CP08_DR_TRT_MEANS, t + 1, 1, 1);
      break;
    }

    double *mt = means + 12 * (R_xlen_t) t;
    mt[0] = T_dr;
    mt[1] = M_dr;
    mt[2] = dp;
    mt[3] = dm;
    mt[4] = T_dr;
    mt[5] = M_dr;
    mt[6] = mu_dr;
    mt[7] = T_ipw;
    mt[8] = M_ipw;
    mt[9] = mu_ipw;
    mt[10] = T_dr;
    mt[11] = mu_dr_decomp;

    /* 3. Analytic cross derivatives C, U_dr, U_ipw, V_dr, V_ipw. */
    double *C_t = cross_out + (R_xlen_t) (p_trt * m) * t;
    double *U_dr_t = U_dr_out + (R_xlen_t) (4 * m) * t;
    double *U_ipw_t = U_ipw_out + (R_xlen_t) (4 * m) * t;
    for (int k = 0; k < m; ++k) {
      SEXP rr_k = VECTOR_ELT(cohort_rows_sexp, k);
      R_xlen_t n_rr = XLENGTH(rr_k);
      if (n_rr == 0) continue;
      const int *rr = INTEGER(rr_k);
      double r_k = rates[k];
      double scale_k = (r_k > CP08_CLIP_EPS) ? gmm_fp_div(-1.0, gmm_fp_mul(N_dbl, r_k)) : 0.0;
      const double *ss_k = pscore_sums + (R_xlen_t) p_trt * k;
      double *C_tk = C_t + (R_xlen_t) p_trt * k;
      for (int j = 0; j < p_trt; ++j) {
        C_tk[j] = gmm_fp_mul(ss_k[j], scale_k);
      }
      long double ldr0 = 0.0L, ldr1 = 0.0L, ldr2 = 0.0L, ldr3 = 0.0L;
      long double lip0 = 0.0L, lip1 = 0.0L, lip2 = 0.0L, lip3 = 0.0L;
      double ddr0 = 0.0, ddr1 = 0.0, ddr2 = 0.0, ddr3 = 0.0;
      double dip0 = 0.0, dip1 = 0.0, dip2 = 0.0, dip3 = 0.0;
      for (R_xlen_t idx = 0; idx < n_rr; ++idx) {
        R_xlen_t i = (R_xlen_t) (rr[idx] - 1);
        double wp_i = pwp[i], wm_i = pwm[i];
        double rdr0 = gmm_fp_div(gmm_fp_mul(wp_i, gmm_fp_sub(r_buf[i], T_dr)), dp);
        double rdr1 = gmm_fp_div(gmm_fp_mul(wm_i, gmm_fp_sub(r_buf[i], M_dr)), dm);
        double rip0 = gmm_fp_div(gmm_fp_mul(wp_i, gmm_fp_sub(dy_buf[i], T_ipw)), dp);
        double rip1 = gmm_fp_div(gmm_fp_mul(wm_i, gmm_fp_sub(dy_buf[i], M_ipw)), dm);
        if (ld) {
          ldr0 += (long double) rdr0;
          ldr1 += (long double) rdr1;
          ldr2 += (long double) wp_i;
          ldr3 += (long double) wm_i;
          lip0 += (long double) rip0;
          lip1 += (long double) rip1;
          lip2 += (long double) wp_i;
          lip3 += (long double) wm_i;
        } else {
          ddr0 = gmm_fp_add(ddr0, rdr0);
          ddr1 = gmm_fp_add(ddr1, rdr1);
          ddr2 = gmm_fp_add(ddr2, wp_i);
          ddr3 = gmm_fp_add(ddr3, wm_i);
          dip0 = gmm_fp_add(dip0, rip0);
          dip1 = gmm_fp_add(dip1, rip1);
          dip2 = gmm_fp_add(dip2, wp_i);
          dip3 = gmm_fp_add(dip3, wm_i);
        }
      }
      double *u_dr_k = U_dr_t + 4 * (R_xlen_t) k;
      double *u_ip_k = U_ipw_t + 4 * (R_xlen_t) k;
      u_dr_k[0] = gmm_fp_mul(ld ? cp08_ld_value(ldr0) : ddr0, scale_k);
      u_dr_k[1] = gmm_fp_mul(ld ? cp08_ld_value(ldr1) : ddr1, scale_k);
      u_dr_k[2] = gmm_fp_mul(ld ? cp08_ld_value(ldr2) : ddr2, scale_k);
      u_dr_k[3] = gmm_fp_mul(ld ? cp08_ld_value(ldr3) : ddr3, scale_k);
      u_ip_k[0] = gmm_fp_mul(ld ? cp08_ld_value(lip0) : dip0, scale_k);
      u_ip_k[1] = gmm_fp_mul(ld ? cp08_ld_value(lip1) : dip1, scale_k);
      u_ip_k[2] = gmm_fp_mul(ld ? cp08_ld_value(lip2) : dip2, scale_k);
      u_ip_k[3] = gmm_fp_mul(ld ? cp08_ld_value(lip3) : dip3, scale_k);
    }

    double *V_dr_t = V_dr_out + (R_xlen_t) (4 * p_trt) * t;
    double *V_ipw_t = V_ipw_out + (R_xlen_t) (4 * p_trt) * t;
    for (int j = 0; j < p_trt; ++j) {
      const double *x_col = px_trt + (R_xlen_t) j * n_x;
      long double l_vdr = 0.0L, l_vip = 0.0L;
      double d_vdr = 0.0, d_vip = 0.0;
      for (R_xlen_t i = 0; i < n_x; ++i) {
        double wm_i = pwm[i];
        double x_ij = x_col[i];
        double w_r = gmm_fp_mul(wm_i, gmm_fp_sub(r_buf[i], M_dr));
        double w_dy = gmm_fp_mul(wm_i, gmm_fp_sub(dy_buf[i], M_ipw));
        double p_dr = gmm_fp_mul(w_r, x_ij);
        double p_ip = gmm_fp_mul(w_dy, x_ij);
        if (ld) {
          l_vdr += (long double) p_dr;
          l_vip += (long double) p_ip;
        } else {
          d_vdr = gmm_fp_add(d_vdr, p_dr);
          d_vip = gmm_fp_add(d_vip, p_ip);
        }
      }
      double v4 = gmm_fp_div(pwm_design[j], N_dbl);
      double *vdr_col = V_dr_t + 4 * (R_xlen_t) j;
      double *vip_col = V_ipw_t + 4 * (R_xlen_t) j;
      vdr_col[0] = 0.0;
      vdr_col[1] = gmm_fp_div(ld ? cp08_ld_value(l_vdr) : d_vdr, denom_v);
      vdr_col[2] = 0.0;
      vdr_col[3] = v4;
      vip_col[0] = 0.0;
      vip_col[1] = gmm_fp_div(ld ? cp08_ld_value(l_vip) : d_vip, denom_v);
      vip_col[2] = 0.0;
      vip_col[3] = v4;
    }

    int part_ok = 1;
    for (int idx = 0; idx < p_trt * m && part_ok; ++idx) {
      if (!R_FINITE(C_t[idx])) part_ok = 0;
    }
    for (int idx = 0; idx < 4 * m && part_ok; ++idx) {
      if (!R_FINITE(U_dr_t[idx]) || !R_FINITE(U_ipw_t[idx])) part_ok = 0;
    }
    for (int idx = 0; idx < 4 * p_trt && part_ok; ++idx) {
      if (!R_FINITE(V_dr_t[idx]) || !R_FINITE(V_ipw_t[idx])) part_ok = 0;
    }
    if (!part_ok) {
      cp08_status(st, CP08_DR_TRT_PARTIALS, t + 1, 1, 1);
      break;
    }
    counts[2] += 1.0;
  }

  static const char *out_names[11] = {
    "status", "means", "d_bar_plus", "d_bar_minus", "d_bar_diff",
    "cross", "U_dr", "U_ipw", "V_dr", "V_ipw", "counts"
  };
  SET_VECTOR_ELT(out_sexp, 0, status_sexp);
  SET_VECTOR_ELT(out_sexp, 1, means_sexp);
  SET_VECTOR_ELT(out_sexp, 2, d_bar_plus_sexp);
  SET_VECTOR_ELT(out_sexp, 3, d_bar_minus_sexp);
  SET_VECTOR_ELT(out_sexp, 4, d_bar_diff_sexp);
  SET_VECTOR_ELT(out_sexp, 5, cross_sexp);
  SET_VECTOR_ELT(out_sexp, 6, U_dr_sexp);
  SET_VECTOR_ELT(out_sexp, 7, U_ipw_sexp);
  SET_VECTOR_ELT(out_sexp, 8, V_dr_sexp);
  SET_VECTOR_ELT(out_sexp, 9, V_ipw_sexp);
  SET_VECTOR_ELT(out_sexp, 10, counts_sexp);
  for (int j = 0; j < 11; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* Pass 1 of the common-design OR fit: visits every support s = 1..ns in
   specs order across bounded cohort_order chunks of size chunk_rows, checking
   posting slots, reading endpoints, subtracting raw differences y_t - y_s
   BEFORE summation, checking the declared intercept column, and accumulating
   per-support (repeat_mass, weight_mass, sumY, n_repeat, n_pair), per-support
   slope moments (sumZ, colmax) when p > 0, and the control-cohort collapsed
   response ystar_i = sum_b (v_b L_b / D) dy_ib when p > 0. */
SEXP C_gmm_cp08_fit_pass1(SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
                          SEXP order_sexp, SEXP a_sexp, SEXP w_sexp, SEXP z_sexp, SEXP sup_sexp,
                          SEXP cols_sexp, SEXP ystar_w_sexp, SEXP ctl_sexp, SEXP long_double_sexp,
                          SEXP max_entries_sexp, SEXP validate_sexp) {
  static const char *fn = "C_gmm_cp08_fit_pass1";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  R_xlen_t ny = cp08_len(0, y_sexp, REALSXP, fn, "y");
  R_xlen_t n_dates = cp08_len(0, date_start_sexp, INTSXP, fn, "date_start");
  R_xlen_t n_post = cp08_len(0, post_a_sexp, INTSXP, fn, "post_a_index");
  if (cp08_len(0, post_obs_sexp, INTSXP, fn, "post_obs_index") != n_post) {
    cp08_malformed(0, fn, "posting A-unit and observation vectors differ in length");
  }
  R_xlen_t n_order = cp08_len(0, order_sexp, INTSXP, fn, "cohort_order");
  R_xlen_t n_frame = cp08_len(0, a_sexp, INTSXP, fn, "a_index");
  if (cp08_len(0, w_sexp, REALSXP, fn, "w") != n_frame || n_frame < 1 || n_frame > INT_MAX) {
    cp08_malformed(0, fn, "weights and A-indices differ in length or are empty");
  }
  int z_nr = 0, z_nc = 0, sup_nr = 0, ns = 0;
  cp08_dims(0, z_sexp, REALSXP, fn, "Z_a", 0, 0, &z_nr, &z_nc);
  cp08_dims(0, sup_sexp, INTSXP, fn, "sup", 4, 0, &sup_nr, &ns);
  R_xlen_t p_x = cp08_len(0, cols_sexp, INTSXP, fn, "cols");
  if (p_x > z_nc || cp08_len(0, ystar_w_sexp, REALSXP, fn, "ystar_w") != (R_xlen_t) ns ||
      cp08_len(0, ctl_sexp, INTSXP, fn, "ctl") != 4 ||
      cp08_len(0, long_double_sexp, LGLSXP, fn, "long_double") != 1 ||
      LOGICAL(long_double_sexp)[0] == NA_LOGICAL) {
    cp08_malformed(0, fn, "cols, ystar_w, ctl or long_double have invalid lengths");
  }
  if (validate_sexp != R_NilValue && !Rf_isFunction(validate_sexp))
    cp08_malformed(0, fn, "original-factor validator must be NULL or a function");
  int p = (int) p_x;
  const int *cols = INTEGER(cols_sexp);
  if (p > 0) cp08_range(0, cols, p, 1, z_nc, fn, "slope columns outside Z_a");
  const int *ctl = INTEGER(ctl_sexp);
  int icpt_col = ctl[0], chunk = ctl[1], ctrl_first = ctl[2], ctrl_last = ctl[3];
  if (icpt_col == NA_INTEGER || icpt_col < 1 || icpt_col > z_nc ||
      chunk == NA_INTEGER || chunk < 1 ||
      ctrl_first == NA_INTEGER || ctrl_last == NA_INTEGER ||
      ctrl_first < 1 || ctrl_last < ctrl_first || (R_xlen_t) ctrl_last > n_order) {
    cp08_malformed(0, fn, "ctl (icpt_col, chunk_rows, ctrl_first, ctrl_last) out of range");
  }
  const int *sup = INTEGER(sup_sexp);
  const double *ystar_w = REAL(ystar_w_sexp);
  const int *p_ord = INTEGER(order_sexp);
  for (int s = 0; s < ns; ++s) {
    int first = sup[4 * (R_xlen_t) s], last = sup[4 * (R_xlen_t) s + 1];
    int ds = sup[4 * (R_xlen_t) s + 2], dt = sup[4 * (R_xlen_t) s + 3];
    if (first == NA_INTEGER || last == NA_INTEGER || first < 1 || last < first || (R_xlen_t) last > n_order ||
        ds == NA_INTEGER || dt == NA_INTEGER || ds < 1 || (R_xlen_t) ds > n_dates ||
        dt < 1 || (R_xlen_t) dt > n_dates || !R_FINITE(ystar_w[s])) {
      cp08_malformed(0, fn, "support range, calendar dates or ystar weight out of range");
    }
    if (ystar_w[s] != 0.0 && (first != ctrl_first || last != ctrl_last)) {
      cp08_malformed(0, fn, "control support range disagrees with ctrl_first..ctrl_last");
    }
    cp08_range(0, p_ord + (first - 1), last - first + 1, 1, (int) n_frame, fn, "support rows outside the frame");
  }
  int n_c = ctrl_last - ctrl_first + 1;
  int max_b = chunk;
  (void) gmm_checked_mul_u64(0, (uint64_t) max_b, 1ULL, max_entries, "C_gmm_cp08_fit_pass1 chunk differences");
  (void) gmm_checked_mul_u64(0, 3ULL, (uint64_t) ns, max_entries, "C_gmm_cp08_fit_pass1 support sums");
  if (p > 0) {
    (void) gmm_checked_mul_u64(0, (uint64_t) p, (uint64_t) ns, max_entries, "C_gmm_cp08_fit_pass1 support slope moments");
    (void) gmm_checked_mul_u64(0, (uint64_t) n_c, 1ULL, max_entries, "C_gmm_cp08_fit_pass1 ystar");
  }

  const double *py = REAL(y_sexp);
  const int *pds = INTEGER(date_start_sexp);
  const int *ppa = INTEGER(post_a_sexp);
  const int *ppo = INTEGER(post_obs_sexp);
  const int *pa = INTEGER(a_sexp);
  const double *pw = REAL(w_sexp);
  const double *pz = REAL(z_sexp);
  const double *pz_icpt = pz + (R_xlen_t) (icpt_col - 1) * z_nr;
  int ld = LOGICAL(long_double_sexp)[0] != 0;

  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP sums_sexp = PROTECT(Rf_allocMatrix(REALSXP, 3, ns));
  SEXP n_counts_sexp = PROTECT(Rf_allocMatrix(INTSXP, 2, ns));
  SEXP sumZ_sexp = PROTECT(Rf_allocMatrix(REALSXP, p, ns));
  SEXP colmax_sexp = PROTECT(Rf_allocMatrix(REALSXP, p, ns));
  SEXP ystar_sexp = PROTECT(Rf_allocVector(REALSXP, p > 0 ? n_c : 0));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 3));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 7));
  const int nprot = 9;

  int *st = INTEGER(status_sexp);
  double *sums = REAL(sums_sexp);
  int *n_counts = INTEGER(n_counts_sexp);
  double *sumZ = REAL(sumZ_sexp);
  double *colmax = REAL(colmax_sexp);
  double *ystar = REAL(ystar_sexp);
  double *counts = REAL(counts_sexp);
  double *dy_buf = (double *) R_alloc((size_t) max_b, sizeof(double));

  cp08_status(st, CP08_OK, 0, 0, 0);
  for (R_xlen_t idx = 0; idx < 3 * (R_xlen_t) ns; ++idx) sums[idx] = 0.0;
  for (R_xlen_t idx = 0; idx < 2 * (R_xlen_t) ns; ++idx) n_counts[idx] = 0;
  for (R_xlen_t idx = 0; idx < (R_xlen_t) p * ns; ++idx) {
    sumZ[idx] = 0.0;
    colmax[idx] = 0.0;
  }
  for (int i = 0; i < (p > 0 ? n_c : 0); ++i) ystar[i] = 0.0;
  counts[0] = counts[1] = counts[2] = 0.0;

  for (int s_idx = 0; s_idx < ns && st[0] == CP08_OK; ++s_idx) {
    int first = sup[4 * (R_xlen_t) s_idx], last = sup[4 * (R_xlen_t) s_idx + 1];
    int dd[2] = {sup[4 * (R_xlen_t) s_idx + 2], sup[4 * (R_xlen_t) s_idx + 3]};
    double yw = ystar_w[s_idx];
    double rep_mass = 0.0, wt_mass = 0.0, sum_y = 0.0;
    int n_rep = 0, n_pr = 0, original_factor_checked = 0;
    double *sz = sumZ + (R_xlen_t) p * s_idx;
    double *cm = colmax + (R_xlen_t) p * s_idx;

    for (int start = first, c_idx = 1; start <= last && st[0] == CP08_OK; start += chunk, ++c_idx) {
      R_CheckUserInterrupt();
      int end = start + chunk - 1;
      if (end > last) end = last;
      int n = end - start + 1;
      const int *ord = p_ord + (start - 1);
      counts[2] += 1.0;

      /* 1. Both posting slots of every chunk row before any endpoint is read. */
      for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
        int a = pa[ord[i] - 1];
        for (int side = 0; side < 2; ++side) {
          int ok = 0;
          if (a != NA_INTEGER && a >= 1 && pds[dd[side] - 1] != NA_INTEGER) {
            int64_t slot = (int64_t) pds[dd[side] - 1] - 1 + (int64_t) a;
            ok = slot >= 1 && slot <= (int64_t) n_post && ppa[slot - 1] == a;
          }
          counts[0] += 1.0;
          if (!ok) {
            cp08_status(st, CP08_SLOT, s_idx + 1, c_idx, i + 1);
            break;
          }
        }
      }
      /* 2. Endpoints, then 3. raw right-minus-left differences BEFORE summation. */
      for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
        int a = pa[ord[i] - 1];
        double v[2];
        for (int side = 0; side < 2; ++side) {
          int o = ppo[(int64_t) pds[dd[side] - 1] - 1 + (int64_t) a - 1];
          v[side] = (o != NA_INTEGER && o >= 1 && (R_xlen_t) o <= ny) ? py[o - 1] : R_NaN;
          counts[1] += 1.0;
        }
        if (!R_FINITE(v[0]) || !R_FINITE(v[1])) {
          cp08_status(st, CP08_ENDPOINT, s_idx + 1, c_idx, i + 1);
          break;
        }
        dy_buf[i] = gmm_fp_sub(v[1], v[0]);
      }
      for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
        if (!R_FINITE(dy_buf[i])) {
          cp08_status(st, CP08_DIFFERENCE, s_idx + 1, c_idx, i + 1);
          break;
        }
      }
      /* 4. The declared intercept column equals one on every observed pair row. */
      for (int i = 0; i < n && st[0] == CP08_OK; ++i) {
        int a = pa[ord[i] - 1];
        if (a < 1 || a > z_nr) cp08_malformed(nprot, fn, "A-unit outside the design rows");
        if (!(pz_icpt[a - 1] == 1.0)) {
          cp08_status(st, CP08_INTERCEPT, s_idx + 1, c_idx, i + 1);
          break;
        }
      }
      if (st[0] != CP08_OK) break;

      /* 5. Chunk moments and ystar accumulation. */
      long double lA = 0.0L, lW = 0.0L, lWD = 0.0L;
      double dA = 0.0, dW = 0.0, dWD = 0.0, response_norm = 0.0;
      int n_a = 0;
      for (int i = 0; i < n; ++i) {
        int r = ord[i] - 1;
        double wi = pw[r];
        double dy = dy_buf[i];
        /* Match a nonfinite raw weighted QR column without squaring it:
           stable QR can represent a finite norm whose square overflows. */
        double weighted_dy = gmm_fp_mul(sqrt(wi), dy);
        response_norm = hypot(response_norm, weighted_dy);
        int raw_square_overflow = !R_FINITE(gmm_fp_mul(weighted_dy, weighted_dy));
        for (int j = 0; j < p && !raw_square_overflow; ++j) {
          double d = gmm_fp_mul((double) (dd[1] - dd[0]),
            pz[(R_xlen_t) (cols[j] - 1) * z_nr + pa[r] - 1]);
          double wd = gmm_fp_mul(sqrt(wi), d);
          raw_square_overflow = !R_FINITE(gmm_fp_mul(wd, wd));
        }
        if (validate_sexp != R_NilValue && !original_factor_checked && raw_square_overflow) {
          SEXP sid_sexp = PROTECT(Rf_ScalarInteger(s_idx + 1));
          SEXP check_call = PROTECT(Rf_lang2(validate_sexp, sid_sexp));
          Rf_eval(check_call, R_GlobalEnv);
          UNPROTECT(2);
          original_factor_checked = 1;
        }
        if (!R_FINITE(weighted_dy) || !R_FINITE(response_norm)) {
          cp08_status(st, 13, s_idx + 1, c_idx, i + 1);
          break;
        }
        if (pa[r] > 0) {
          if (ld) lA += (long double) wi; else dA = gmm_fp_add(dA, wi);
          ++n_a;
        }
        if (ld) lW += (long double) wi; else dW = gmm_fp_add(dW, wi);
        double prod = gmm_fp_mul(wi, dy);
        if (ld) lWD += (long double) prod; else dWD = gmm_fp_add(dWD, prod);
        if (p > 0 && yw != 0.0) {
          int c_pos = (start - ctrl_first) + i;
          ystar[c_pos] = gmm_fp_add(ystar[c_pos], gmm_fp_mul(yw, dy));
        }
      }
      if (st[0] != CP08_OK) break;
      rep_mass = gmm_fp_add(rep_mass, ld ? cp08_ld_value(lA) : dA);
      wt_mass = gmm_fp_add(wt_mass, ld ? cp08_ld_value(lW) : dW);
      sum_y = gmm_fp_add(sum_y, ld ? cp08_ld_value(lWD) : dWD);
      if (!R_FINITE(rep_mass) || !R_FINITE(wt_mass) || !R_FINITE(sum_y)) {
        cp08_status(st, 13, s_idx + 1, c_idx, 1);
        break;
      }
      n_rep += n_a;
      n_pr += n;

      for (int j = 0; j < p; ++j) {
        const double *pz_col = pz + (R_xlen_t) (cols[j] - 1) * z_nr;
        long double lZ = 0.0L;
        double dZ = 0.0;
        double cmax = cm[j];
        for (int i = 0; i < n; ++i) {
          int r = ord[i] - 1;
          int a = pa[r];
          double zij = pz_col[a - 1];
          double wz = gmm_fp_mul(pw[r], zij);
          if (ld) lZ += (long double) wz; else dZ = gmm_fp_add(dZ, wz);
          double az = fabs(zij);
          if (az > cmax) cmax = az;
        }
        sz[j] = gmm_fp_add(sz[j], ld ? cp08_ld_value(lZ) : dZ);
        cm[j] = cmax;
      }
    }
    sums[3 * (R_xlen_t) s_idx] = rep_mass;
    sums[3 * (R_xlen_t) s_idx + 1] = wt_mass;
    sums[3 * (R_xlen_t) s_idx + 2] = sum_y;
    n_counts[2 * (R_xlen_t) s_idx] = n_rep;
    n_counts[2 * (R_xlen_t) s_idx + 1] = n_pr;
  }

  static const char *out_names[7] = {"status", "sums", "n_counts", "sumZ", "colmax", "ystar", "counts"};
  SET_VECTOR_ELT(out_sexp, 0, status_sexp);
  SET_VECTOR_ELT(out_sexp, 1, sums_sexp);
  SET_VECTOR_ELT(out_sexp, 2, n_counts_sexp);
  SET_VECTOR_ELT(out_sexp, 3, sumZ_sexp);
  SET_VECTOR_ELT(out_sexp, 4, colmax_sexp);
  SET_VECTOR_ELT(out_sexp, 5, ystar_sexp);
  SET_VECTOR_ELT(out_sexp, 6, counts_sexp);
  for (int j = 0; j < 7; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* Pass 2 of the common-design OR fit: after alpha_b and beta are solved prior
   to Pass 2, visits each control block b = 1..ng across bounded chunks of the
   never-treated cohort, computing from raw differences dy_ib = y_it - y_is:
   r_ib = (dy_ib - alpha_b) - L_b (z_i' beta), raw-weighted forcing moments
   g1_b = sum w_i r_ib, g_slope_b = sum w_i (L_b z_i) r_ib, e2_b = sum w_i r_ib^2,
   and (when p_keep > 0) K_b = sum (v_b w_i) ((dy_ib - L_b ystar_i) - e_bar_b)^2
   accumulated strictly as non-negative squares. */
SEXP C_gmm_cp08_fit_pass2(SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
                          SEXP order_sexp, SEXP a_sexp, SEXP w_sexp, SEXP z_sexp, SEXP blk_sexp,
                          SEXP blk_par_sexp, SEXP cols_sexp, SEXP beta_sexp, SEXP ystar_sexp,
                          SEXP ctl_sexp, SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_fit_pass2";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  R_xlen_t ny = cp08_len(0, y_sexp, REALSXP, fn, "y");
  R_xlen_t n_dates = cp08_len(0, date_start_sexp, INTSXP, fn, "date_start");
  R_xlen_t n_post = cp08_len(0, post_a_sexp, INTSXP, fn, "post_a_index");
  if (cp08_len(0, post_obs_sexp, INTSXP, fn, "post_obs_index") != n_post) {
    cp08_malformed(0, fn, "posting A-unit and observation vectors differ in length");
  }
  R_xlen_t n_order = cp08_len(0, order_sexp, INTSXP, fn, "cohort_order");
  R_xlen_t n_frame = cp08_len(0, a_sexp, INTSXP, fn, "a_index");
  if (cp08_len(0, w_sexp, REALSXP, fn, "w") != n_frame || n_frame < 1 || n_frame > INT_MAX) {
    cp08_malformed(0, fn, "weights and A-indices differ in length or are empty");
  }
  int z_nr = 0, z_nc = 0, blk_nr = 0, ng = 0, tmp = 0;
  cp08_dims(0, z_sexp, REALSXP, fn, "Z_a", 0, 0, &z_nr, &z_nc);
  cp08_dims(0, blk_sexp, INTSXP, fn, "blk", 2, 0, &blk_nr, &ng);
  cp08_dims(0, blk_par_sexp, REALSXP, fn, "blk_par", 4, ng, &blk_nr, &tmp);
  R_xlen_t p_x = cp08_len(0, cols_sexp, INTSXP, fn, "cols");
  if (p_x > z_nc || cp08_len(0, beta_sexp, REALSXP, fn, "beta") != p_x ||
      cp08_len(0, ctl_sexp, INTSXP, fn, "ctl") != 3) {
    cp08_malformed(0, fn, "cols, beta or ctl have invalid lengths");
  }
  int p = (int) p_x;
  const int *cols = INTEGER(cols_sexp);
  if (p > 0) cp08_range(0, cols, p, 1, z_nc, fn, "slope columns outside Z_a");
  const int *ctl = INTEGER(ctl_sexp);
  int ctrl_first = ctl[0], ctrl_last = ctl[1], chunk = ctl[2];
  if (ctrl_first == NA_INTEGER || ctrl_last == NA_INTEGER || chunk == NA_INTEGER ||
      ctrl_first < 1 || ctrl_last < ctrl_first || (R_xlen_t) ctrl_last > n_order || chunk < 1) {
    cp08_malformed(0, fn, "ctl (ctrl_first, ctrl_last, chunk_rows) out of range");
  }
  int n_c = ctrl_last - ctrl_first + 1;
  R_xlen_t ny_star = cp08_len(0, ystar_sexp, REALSXP, fn, "ystar");
  if ((p > 0 && ny_star != (R_xlen_t) n_c) || (p == 0 && ny_star != 0 && ny_star != (R_xlen_t) n_c)) {
    cp08_malformed(0, fn, "ystar length does not match control cohort");
  }
  const int *blk = INTEGER(blk_sexp);
  cp08_range(0, blk, 2 * (R_xlen_t) ng, 1, (int) n_dates, fn, "block dates outside calendar");
  const int *p_ord = INTEGER(order_sexp);
  cp08_range(0, p_ord + (ctrl_first - 1), n_c, 1, (int) n_frame, fn, "control rows outside the frame");
  (void) gmm_checked_mul_u64(0, 3ULL, (uint64_t) ng, max_entries, "C_gmm_cp08_fit_pass2 block moments");
  if (p > 0) {
    (void) gmm_checked_mul_u64(0, (uint64_t) p, (uint64_t) ng, max_entries, "C_gmm_cp08_fit_pass2 block slope moments");
  }

  const double *py = REAL(y_sexp);
  const int *pds = INTEGER(date_start_sexp);
  const int *ppa = INTEGER(post_a_sexp);
  const int *ppo = INTEGER(post_obs_sexp);
  const int *pa = INTEGER(a_sexp);
  const double *pw = REAL(w_sexp);
  const double *pz = REAL(z_sexp);
  const double *bp = REAL(blk_par_sexp);
  const double *beta = REAL(beta_sexp);
  const double *ystar = REAL(ystar_sexp);

  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP moments_sexp = PROTECT(Rf_allocMatrix(REALSXP, 3, ng));
  SEXP g_slope_sexp = PROTECT(Rf_allocMatrix(REALSXP, p, ng));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 3));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 4));
  const int nprot = 6;

  int *st = INTEGER(status_sexp);
  double *moments = REAL(moments_sexp);
  double *g_slope = REAL(g_slope_sexp);
  double *counts = REAL(counts_sexp);

  cp08_status(st, CP08_OK, 0, 0, 0);
  for (R_xlen_t idx = 0; idx < 3 * (R_xlen_t) ng; ++idx) moments[idx] = 0.0;
  for (R_xlen_t idx = 0; idx < (R_xlen_t) p * ng; ++idx) g_slope[idx] = 0.0;
  counts[0] = counts[1] = counts[2] = 0.0;

  for (int b = 0; b < ng && st[0] == CP08_OK; ++b) {
    int ds = blk[2 * (R_xlen_t) b], dt = blk[2 * (R_xlen_t) b + 1];
    double L_b = bp[4 * (R_xlen_t) b], v_b = bp[4 * (R_xlen_t) b + 1];
    double alpha_b = bp[4 * (R_xlen_t) b + 2], e_bar_b = bp[4 * (R_xlen_t) b + 3];
    double g1 = 0.0, e2 = 0.0, K_b = 0.0;
    double *gs = g_slope + (R_xlen_t) p * b;

    for (int start = ctrl_first, c_idx = 1; start <= ctrl_last && st[0] == CP08_OK; start += chunk, ++c_idx) {
      R_CheckUserInterrupt();
      int end = start + chunk - 1;
      if (end > ctrl_last) end = ctrl_last;
      int n = end - start + 1;
      const int *ord = p_ord + (start - 1);
      counts[2] += 1.0;

      for (int i = 0; i < n; ++i) {
        int r = ord[i] - 1;
        int a = pa[r];
        if (a < 1 || a > z_nr) cp08_malformed(nprot, fn, "A-unit outside the design rows");
        int64_t ss = (int64_t) pds[ds - 1] - 1 + (int64_t) a;
        int64_t tt = (int64_t) pds[dt - 1] - 1 + (int64_t) a;
        counts[0] += 2.0;
        if (pds[ds - 1] == NA_INTEGER || pds[dt - 1] == NA_INTEGER ||
            ss < 1 || ss > (int64_t) n_post || tt < 1 || tt > (int64_t) n_post ||
            ppa[ss - 1] != a || ppa[tt - 1] != a) {
          cp08_status(st, CP08_SLOT, b + 1, c_idx, i + 1);
          break;
        }
        int os = ppo[ss - 1];
        int ot = ppo[tt - 1];
        double vs = (os != NA_INTEGER && os >= 1 && (R_xlen_t) os <= ny) ? py[os - 1] : R_NaN;
        double vt = (ot != NA_INTEGER && ot >= 1 && (R_xlen_t) ot <= ny) ? py[ot - 1] : R_NaN;
        counts[1] += 2.0;
        if (!R_FINITE(vs) || !R_FINITE(vt)) {
          cp08_status(st, CP08_ENDPOINT, b + 1, c_idx, i + 1);
          break;
        }
        double dy = gmm_fp_sub(vt, vs);
        if (!R_FINITE(dy)) {
          cp08_status(st, CP08_DIFFERENCE, b + 1, c_idx, i + 1);
          break;
        }
        double wi = pw[r];
        double db = 0.0;
        for (int j = 0; j < p; ++j) {
          double dj = gmm_fp_mul(L_b, pz[(R_xlen_t) (cols[j] - 1) * z_nr + (a - 1)]);
          db = gmm_fp_add(db, gmm_fp_mul(dj, beta[j]));
        }
        double res = gmm_fp_sub(gmm_fp_sub(dy, alpha_b), db);
        double wr = gmm_fp_mul(wi, res);
        g1 = gmm_fp_add(g1, wr);
        e2 = gmm_fp_add(e2, gmm_fp_mul(wr, res));
        for (int j = 0; j < p; ++j) {
          double dj = gmm_fp_mul(L_b, pz[(R_xlen_t) (cols[j] - 1) * z_nr + (a - 1)]);
          gs[j] = gmm_fp_add(gs[j], gmm_fp_mul(wr, dj));
        }
        if (p > 0) {
          int c_pos = (start - ctrl_first) + i;
          double eb = gmm_fp_sub(dy, gmm_fp_mul(L_b, ystar[c_pos]));
          double dev = gmm_fp_sub(eb, e_bar_b);
          double q_ib = gmm_fp_mul(v_b, wi);
          K_b = gmm_fp_add(K_b, gmm_fp_mul(q_ib, gmm_fp_mul(dev, dev)));
        }
      }
    }
    moments[3 * (R_xlen_t) b] = g1;
    moments[3 * (R_xlen_t) b + 1] = e2;
    moments[3 * (R_xlen_t) b + 2] = K_b;
  }

  static const char *out_names[4] = {"status", "moments", "g_slope", "counts"};
  SET_VECTOR_ELT(out_sexp, 0, status_sexp);
  SET_VECTOR_ELT(out_sexp, 1, moments_sexp);
  SET_VECTOR_ELT(out_sexp, 2, g_slope_sexp);
  SET_VECTOR_ELT(out_sexp, 3, counts_sexp);
  for (int j = 0; j < 4; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

typedef struct {
  int p;
  int n_a;
  const double *z_a;
  const int *cols;
  const double *beta;
  const double *R_mat;
  const int *pivot;
  double normalizer;
  int ng_beta;
  const int *beta_blk;
  const double *beta_par;
  int n_ctrl;
  const int *ctrl_sup;
  const double *forcing_beta;
  SEXP tail_sexp;
  double *xi_beta;
} cp08_slope_t;

static inline double cp08_z_elem(const cp08_slope_t *sl, int a, int k) {
  int col = sl->cols != NULL ? (sl->cols[k] - 1) : k;
  return sl->z_a[(R_xlen_t) col * sl->n_a + (a - 1)];
}

static int cp08_or_slope_if(int nprotect, const char *fn, cp08_chunk_t *g, int ns, const int *sup,
                            const double *sp, SEXP slope_sexp, double max_entries, int *st,
                            cp08_slope_t *out) {
  if (slope_sexp == R_NilValue || TYPEOF(slope_sexp) != VECSXP || XLENGTH(slope_sexp) != 10) {
    cp08_malformed(nprotect, fn, "slope descriptor must be a 10-element list");
  }
  int n_a = 0, p = 0, z_nc = 0, r_nr = 0, r_nc = 0, bp_nr = 0, ng_beta = 0, fb_nr = 0, n_ctrl = 0;
  SEXP z_elt = VECTOR_ELT(slope_sexp, 0);
  const double *z_a = NULL;
  const int *cols = NULL;
  if (z_elt != R_NilValue && TYPEOF(z_elt) == VECSXP && XLENGTH(z_elt) == 2) {
    cp08_dims(nprotect, VECTOR_ELT(z_elt, 0), REALSXP, fn, "slope$Z_a", 0, 0, &n_a, &z_nc);
    p = (int) cp08_len(nprotect, VECTOR_ELT(z_elt, 1), INTSXP, fn, "slope$cols");
    if (p < 1 || p > z_nc) cp08_malformed(nprotect, fn, "slope$cols length out of range");
    cols = INTEGER(VECTOR_ELT(z_elt, 1));
    cp08_range(nprotect, cols, p, 1, z_nc, fn, "slope$cols outside Z_a columns");
    z_a = REAL(VECTOR_ELT(z_elt, 0));
  } else {
    cp08_dims(nprotect, z_elt, REALSXP, fn, "slope$z_a", 0, 0, &n_a, &p);
    z_a = REAL(z_elt);
  }
  if (cp08_len(nprotect, VECTOR_ELT(slope_sexp, 1), REALSXP, fn, "slope$beta") != (R_xlen_t) p) {
    cp08_malformed(nprotect, fn, "slope$beta length mismatch");
  }
  cp08_dims(nprotect, VECTOR_ELT(slope_sexp, 2), REALSXP, fn, "slope$R", p, p, &r_nr, &r_nc);
  if (cp08_len(nprotect, VECTOR_ELT(slope_sexp, 3), INTSXP, fn, "slope$pivot") != (R_xlen_t) p) {
    cp08_malformed(nprotect, fn, "slope$pivot length mismatch");
  }
  const int *pivot = INTEGER(VECTOR_ELT(slope_sexp, 3));
  cp08_range(nprotect, pivot, p, 1, p, fn, "slope$pivot outside 1..p");
  if (cp08_len(nprotect, VECTOR_ELT(slope_sexp, 4), REALSXP, fn, "slope$normalizer") != 1) {
    cp08_malformed(nprotect, fn, "slope$normalizer must be a scalar");
  }
  double normalizer = REAL(VECTOR_ELT(slope_sexp, 4))[0];
  if (!R_FINITE(normalizer) || normalizer <= 0.0) {
    cp08_malformed(nprotect, fn, "slope$normalizer must be positive finite");
  }
  ng_beta = (int) cp08_len(nprotect, VECTOR_ELT(slope_sexp, 5), INTSXP, fn, "slope$beta_blk");
  if (ng_beta < 1) cp08_malformed(nprotect, fn, "slope$beta_blk must be nonempty");
  const int *beta_blk = INTEGER(VECTOR_ELT(slope_sexp, 5));
  cp08_range(nprotect, beta_blk, ng_beta, 1, ns, fn, "slope$beta_blk outside supports");
  cp08_dims(nprotect, VECTOR_ELT(slope_sexp, 6), REALSXP, fn, "slope$beta_par", 2 + p, ng_beta, &bp_nr, &ng_beta);
  n_ctrl = (int) cp08_len(nprotect, VECTOR_ELT(slope_sexp, 7), INTSXP, fn, "slope$ctrl_sup");
  if (n_ctrl < 1) cp08_malformed(nprotect, fn, "slope$ctrl_sup must be nonempty");
  const int *ctrl_sup = INTEGER(VECTOR_ELT(slope_sexp, 7));
  cp08_range(nprotect, ctrl_sup, n_ctrl, 1, ns, fn, "slope$ctrl_sup outside supports");
  cp08_dims(nprotect, VECTOR_ELT(slope_sexp, 8), REALSXP, fn, "slope$forcing_beta", p, n_ctrl, &fb_nr, &n_ctrl);

  int B = g->B;
  uint64_t bp_len = gmm_checked_mul_u64(nprotect, (uint64_t) B, (uint64_t) p, max_entries,
                                        "cp08_or_slope_if slope IF (B x p)");
  double *h_beta = (double *) R_alloc((size_t) bp_len, sizeof(double));
  double *xi_beta = (double *) R_alloc((size_t) bp_len, sizeof(double));
  double *v_sol = (double *) R_alloc((size_t) p, sizeof(double));
  double *y_sol = (double *) R_alloc((size_t) p, sizeof(double));
  double *x_sol = (double *) R_alloc((size_t) p, sizeof(double));
  for (R_xlen_t idx = 0; idx < (R_xlen_t) bp_len; ++idx) {
    h_beta[idx] = 0.0;
    xi_beta[idx] = 0.0;
  }

  const double *beta = REAL(VECTOR_ELT(slope_sexp, 1));
  const double *R_mat = REAL(VECTOR_ELT(slope_sexp, 2));
  const double *beta_par = REAL(VECTOR_ELT(slope_sexp, 6));
  const double *forcing_beta = REAL(VECTOR_ELT(slope_sexp, 8));
  out->p = p;
  out->n_a = n_a;
  out->z_a = z_a;
  out->cols = cols;

  for (int kb = 0; kb < ng_beta; ++kb) {
    R_CheckUserInterrupt();
    int c = beta_blk[kb] - 1;
    const int *s = sup + 3 * (R_xlen_t) c;
    double v = fabs(sp[3 * (R_xlen_t) c + 2]);
    const double *bpar = beta_par + (R_xlen_t) (2 + p) * kb;
    double L_b = bpar[0], alpha_b = bpar[1];
    const double *mean_d = bpar + 2;
    for (int i = 0; i < B; ++i) {
      int r = g->rows[i] - 1;
      if (g->ci[r] != s[0]) continue;
      int a = g->a[r];
      if (a < 1 || a > n_a) cp08_malformed(nprotect, fn, "A-unit outside Z_a rows");
      double dy = cp08_dy(g, i, s[1] - 1, s[2] - 1);
      if (!R_FINITE(dy)) {
        cp08_status(st, CP08_GATHER, 2, kb + 1, i + 1);
        return 0;
      }
      double q = gmm_fp_mul(g->w[r], v);
      double db = 0.0;
      for (int j = 0; j < p; ++j) {
        double dj = gmm_fp_mul(L_b, cp08_z_elem(out, a, j));
        db = gmm_fp_add(db, gmm_fp_mul(dj, beta[j]));
      }
      double res = gmm_fp_sub(gmm_fp_sub(dy, alpha_b), db);
      double qr_val = gmm_fp_mul(q, res);
      for (int j = 0; j < p; ++j) {
        double dj = gmm_fp_mul(L_b, cp08_z_elem(out, a, j));
        double term = gmm_fp_mul(qr_val, gmm_fp_sub(dj, mean_d[j]));
        h_beta[(R_xlen_t) j * B + i] = gmm_fp_add(h_beta[(R_xlen_t) j * B + i], term);
      }
    }
  }

  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < B; ++i) {
      double uf = 0.0;
      for (int k = 0; k < n_ctrl; ++k) {
        int c = ctrl_sup[k] - 1;
        const int *s = sup + 3 * (R_xlen_t) c;
        const double *sv = sp + 3 * (R_xlen_t) c;
        double u_ik = cp08_up(g, i, s[0], sv[0], sv[1]);
        uf = gmm_fp_add(uf, gmm_fp_mul(u_ik, forcing_beta[(R_xlen_t) k * p + j]));
      }
      h_beta[(R_xlen_t) j * B + i] = gmm_fp_add(h_beta[(R_xlen_t) j * B + i], uf);
      if (!R_FINITE(h_beta[(R_xlen_t) j * B + i])) {
        cp08_status(st, CP08_SLOPE_RHS, 2, j + 1, i + 1);
        return 0;
      }
    }
  }

  for (int i = 0; i < B; ++i) {
    for (int r = 0; r < p; ++r) {
      v_sol[r] = h_beta[(R_xlen_t) (pivot[r] - 1) * B + i];
    }
    for (int r = 0; r < p; ++r) {
      double s = v_sol[r];
      for (int c = 0; c < r; ++c) {
        s = gmm_fp_sub(s, gmm_fp_mul(R_mat[(R_xlen_t) r * p + c], y_sol[c]));
      }
      y_sol[r] = gmm_fp_div(s, R_mat[(R_xlen_t) r * p + r]);
    }
    for (int r = p - 1; r >= 0; --r) {
      double s = y_sol[r];
      for (int c = r + 1; c < p; ++c) {
        s = gmm_fp_sub(s, gmm_fp_mul(R_mat[(R_xlen_t) c * p + r], x_sol[c]));
      }
      x_sol[r] = gmm_fp_div(s, R_mat[(R_xlen_t) r * p + r]);
    }
    for (int r = 0; r < p; ++r) {
      double val = gmm_fp_mul(x_sol[r], normalizer);
      if (!R_FINITE(val)) {
        cp08_status(st, CP08_SLOPE_GRAM, 2, r + 1, i + 1);
        return 0;
      }
      xi_beta[(R_xlen_t) (pivot[r] - 1) * B + i] = val;
    }
  }

  out->p = p;
  out->n_a = n_a;
  out->z_a = z_a;
  out->beta = beta;
  out->R_mat = R_mat;
  out->pivot = pivot;
  out->normalizer = normalizer;
  out->ng_beta = ng_beta;
  out->beta_blk = beta_blk;
  out->beta_par = beta_par;
  out->n_ctrl = n_ctrl;
  out->ctrl_sup = ctrl_sup;
  out->forcing_beta = forcing_beta;
  out->tail_sexp = VECTOR_ELT(slope_sexp, 9);
  out->xi_beta = xi_beta;
  return 1;
}

/* Eligible non-detailed OR IF chunk (evaluate_or_baseline_if_chunk with the
   complete-date reader, or_pooled_if_core at p >= 0 with intercepts). sup is
   3 x ns (cohort, gather columns of s and t), sup_par 3 x ns (rate, massA, v)
   in pre$sids order; blk the control-support position of every requested
   block, blk_par 3 x ng (alpha, mass, forcing_alpha) in pre$gids order; lnk
   2 x nl (target-support and block positions), lnk_par 3 x nl (mu, d,
   rate_loading) in the requested link order. When p > 0, dimnames[[1]] holds
   the slope descriptor list. */
SEXP C_gmm_cp08_or_if_chunk(SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
                            SEXP rows_sexp, SEXP ci_sexp, SEXP a_sexp, SEXP w_sexp, SEXP dates_sexp,
                            SEXP sup_sexp, SEXP sup_par_sexp, SEXP blk_sexp, SEXP blk_par_sexp,
                            SEXP lnk_sexp, SEXP lnk_par_sexp, SEXP dimnames_sexp, SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_or_if_chunk";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  cp08_chunk_t g;
  cp08_chunk_args(0, fn, &g, rows_sexp, ci_sexp, a_sexp, w_sexp);
  cp08_panel_args(0, fn, &g, y_sexp, date_start_sexp, post_a_sexp, post_obs_sexp, dates_sexp);
  int nr = 0, ns = 0, ng = 0, nl = 0, tmp = 0;
  cp08_dims(0, sup_sexp, INTSXP, fn, "sup", 3, 0, &nr, &ns);
  cp08_dims(0, sup_par_sexp, REALSXP, fn, "sup_par", 3, ns, &nr, &tmp);
  ng = (int) cp08_len(0, blk_sexp, INTSXP, fn, "blk");
  if (ng < 1) cp08_malformed(0, fn, "no requested block");
  cp08_dims(0, blk_par_sexp, REALSXP, fn, "blk_par", 3, ng, &nr, &tmp);
  cp08_dims(0, lnk_sexp, INTSXP, fn, "lnk", 2, 0, &nr, &nl);
  cp08_dims(0, lnk_par_sexp, REALSXP, fn, "lnk_par", 3, nl, &nr, &tmp);
  cp08_dimnames(0, dimnames_sexp, nl, fn);
  SEXP slope_sexp = (dimnames_sexp != R_NilValue) ? VECTOR_ELT(dimnames_sexp, 0) : R_NilValue;
  const int *sup = INTEGER(sup_sexp);
  const double *sp = REAL(sup_par_sexp);
  const int *blk = INTEGER(blk_sexp);
  const double *bp = REAL(blk_par_sexp);
  const int *lnk = INTEGER(lnk_sexp);
  const double *lp = REAL(lnk_par_sexp);
  for (int j = 0; j < ns; ++j) {
    cp08_range(0, sup + 3 * (R_xlen_t) j + 1, 2, 1, g.nd, fn, "support gather columns outside the dates");
    if (sup[3 * (R_xlen_t) j] == NA_INTEGER) cp08_malformed(0, fn, "support cohort is NA");
  }
  cp08_range(0, blk, ng, 1, ns, fn, "block control supports outside the requested supports");
  for (int j = 0; j < nl; ++j) {
    cp08_range(0, lnk + 2 * (R_xlen_t) j, 1, 1, ns, fn, "link supports outside the requested supports");
    cp08_range(0, lnk + 2 * (R_xlen_t) j + 1, 1, 1, ng, fn, "link blocks outside the requested blocks");
  }
  int B = g.B;
  uint64_t phi_len = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) nl, max_entries, "C_gmm_cp08_or_if_chunk IF (B x L)");
  uint64_t alpha_len = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) ng, max_entries, "C_gmm_cp08_or_if_chunk block IF columns");
  uint64_t gather_len = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) g.nd, max_entries, "C_gmm_cp08_or_if_chunk endpoint gather");

  SEXP phi_sexp = PROTECT(Rf_allocMatrix(REALSXP, B, nl));
  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 4));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  const int nprot = 5;
  double *phi = REAL(phi_sexp);
  if (gmm_buffers_overlap(phi, (R_xlen_t) phi_len, g.y, g.ny) || gmm_buffers_overlap(phi, (R_xlen_t) phi_len, g.w, g.n_frame)) {
    gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_cp08_or_if_chunk`.", "cdid_gmm_bound_error");
  }
  cp08_set_dimnames(phi_sexp, dimnames_sexp);
  for (R_xlen_t idx = 0; idx < (R_xlen_t) phi_len; ++idx) phi[idx] = 0.0;
  double *alpha = (double *) R_alloc((size_t) alpha_len, sizeof(double));
  g.Y = (double *) R_alloc((size_t) gather_len, sizeof(double));
  g.state = (unsigned char *) R_alloc((size_t) gather_len, sizeof(unsigned char));
  for (R_xlen_t idx = 0; idx < (R_xlen_t) gather_len; ++idx) {
    g.Y[idx] = R_NaN;
    g.state[idx] = 0;
  }
  int *st = INTEGER(status_sexp);
  cp08_status(st, CP08_OK, 0, 0, 0);

  /* 1. Rate stage: every requested support in order. */
  for (int j = 0; j < ns && st[0] == CP08_OK; ++j) {
    R_CheckUserInterrupt();
    const int *s = sup + 3 * (R_xlen_t) j;
    int code = CP08_OK;
    int row = cp08_pair_rows(&g, s, 1, s[1] - 1, s[2] - 1, 1, &code);
    if (row >= 0) cp08_status(st, code, 1, j + 1, row + 1);
  }
  cp08_slope_t sl;
  sl.p = 0;
  const double *blk_lm = NULL;
  const double *lnk_c = NULL;
  if (st[0] == CP08_OK && slope_sexp != R_NilValue) {
    if (cp08_or_slope_if(nprot, fn, &g, ns, sup, sp, slope_sexp, max_entries, st, &sl)) {
      if (sl.tail_sexp == R_NilValue || TYPEOF(sl.tail_sexp) != VECSXP || XLENGTH(sl.tail_sexp) != 2) {
        cp08_malformed(nprot, fn, "OR slope tail must be list(blk_L_meand, lnk_c)");
      }
      int r1 = 0, c1 = 0, r2 = 0, c2 = 0;
      cp08_dims(nprot, VECTOR_ELT(sl.tail_sexp, 0), REALSXP, fn, "blk_L_meand", 1 + sl.p, ng, &r1, &c1);
      cp08_dims(nprot, VECTOR_ELT(sl.tail_sexp, 1), REALSXP, fn, "lnk_c", sl.p, nl, &r2, &c2);
      blk_lm = REAL(VECTOR_ELT(sl.tail_sexp, 0));
      lnk_c = REAL(VECTOR_ELT(sl.tail_sexp, 1));
    }
  }
  /* 2. Block intercept IF columns. */
  for (int b = 0; b < ng && st[0] == CP08_OK; ++b) {
    R_CheckUserInterrupt();
    int c = blk[b] - 1;
    const int *s = sup + 3 * (R_xlen_t) c;
    double rate = sp[3 * (R_xlen_t) c], mass_a = sp[3 * (R_xlen_t) c + 1], v = sp[3 * (R_xlen_t) c + 2];
    double alpha_b = bp[3 * (R_xlen_t) b], mass = bp[3 * (R_xlen_t) b + 1], fa = bp[3 * (R_xlen_t) b + 2];
    double L_b = (sl.p > 0) ? blk_lm[(R_xlen_t) (1 + sl.p) * b] : 0.0;
    const double *mean_d = (sl.p > 0) ? (blk_lm + (R_xlen_t) (1 + sl.p) * b + 1) : NULL;
    double *col = alpha + (R_xlen_t) b * B;
    for (int i = 0; i < B; ++i) {
      double h = 0.0;
      int r = g.rows[i] - 1;
      if (g.ci[r] == s[0]) {
        double q = gmm_fp_mul(g.w[r], v);
        double db = 0.0;
        if (sl.p > 0) {
          int a = g.a[r];
          for (int k = 0; k < sl.p; ++k) {
            double dk = gmm_fp_mul(L_b, cp08_z_elem(&sl, a, k));
            db = gmm_fp_add(db, gmm_fp_mul(dk, sl.beta[k]));
          }
        }
        double res = gmm_fp_sub(gmm_fp_sub(cp08_dy(&g, i, s[1] - 1, s[2] - 1), alpha_b), db);
        h = gmm_fp_add(0.0, gmm_fp_mul(q, res));
      }
      int nz = 0;
      double forced = cp08_term(0.0, fa, cp08_up(&g, i, s[0], rate, mass_a), &nz);
      double xi_md = 0.0;
      for (int k = 0; k < sl.p; ++k) {
        xi_md = gmm_fp_add(xi_md, gmm_fp_mul(sl.xi_beta[(R_xlen_t) k * B + i], mean_d[k]));
      }
      col[i] = gmm_fp_sub(gmm_fp_div(gmm_fp_add(h, forced), mass), xi_md);
    }
  }
  /* 3. Link columns in the requested order, then IF finiteness. */
  for (int j = 0; j < nl && st[0] == CP08_OK; ++j) {
    R_CheckUserInterrupt();
    int c = lnk[2 * (R_xlen_t) j] - 1, b = lnk[2 * (R_xlen_t) j + 1] - 1;
    const int *s = sup + 3 * (R_xlen_t) c;
    double rate = sp[3 * (R_xlen_t) c], mass_a = sp[3 * (R_xlen_t) c + 1], v = sp[3 * (R_xlen_t) c + 2];
    double alpha_b = bp[3 * (R_xlen_t) b];
    double L_b = (sl.p > 0) ? blk_lm[(R_xlen_t) (1 + sl.p) * b] : 0.0;
    const double *lc = (sl.p > 0) ? (lnk_c + (R_xlen_t) sl.p * j) : NULL;
    double mu = lp[3 * (R_xlen_t) j], d = lp[3 * (R_xlen_t) j + 1], rl = lp[3 * (R_xlen_t) j + 2];
    const double *acol = alpha + (R_xlen_t) b * B;
    double *out = phi + (R_xlen_t) j * B;
    for (int i = 0; i < B; ++i) {
      double direct = 0.0;
      int r = g.rows[i] - 1;
      if (g.ci[r] == s[0]) {
        double q = gmm_fp_mul(g.w[r], v);
        double db = 0.0;
        if (sl.p > 0) {
          int a = g.a[r];
          for (int k = 0; k < sl.p; ++k) {
            double dk = gmm_fp_mul(L_b, cp08_z_elem(&sl, a, k));
            db = gmm_fp_add(db, gmm_fp_mul(dk, sl.beta[k]));
          }
        }
        double res = gmm_fp_sub(gmm_fp_sub(cp08_dy(&g, i, s[1] - 1, s[2] - 1), alpha_b), db);
        direct = gmm_fp_div(gmm_fp_mul(q, gmm_fp_sub(res, mu)), d);
      }
      int nz = 0;
      double load = cp08_term(0.0, rl, cp08_up(&g, i, s[0], rate, mass_a), &nz);
      double xi_c = 0.0;
      for (int k = 0; k < sl.p; ++k) {
        xi_c = gmm_fp_add(xi_c, gmm_fp_mul(sl.xi_beta[(R_xlen_t) k * B + i], lc[k]));
      }
      out[i] = gmm_fp_add(gmm_fp_sub(gmm_fp_sub(direct, acol[i]), xi_c), load);
    }
  }
  for (int j = 0; j < nl && st[0] == CP08_OK; ++j) {
    const double *out = phi + (R_xlen_t) j * B;
    for (int i = 0; i < B; ++i) {
      if (!R_FINITE(out[i])) {
        cp08_status(st, CP08_OR_LINK_IF, 3, j + 1, i + 1);
        break;
      }
    }
  }
  double *counts = REAL(counts_sexp);
  counts[0] = g.slot_checks;
  counts[1] = g.reads;
  counts[2] = (double) g.nd;
  counts[3] = (double) phi_len;
  static const char *out_names[3] = {"phi", "status", "counts"};
  SET_VECTOR_ELT(out_sexp, 0, phi_sexp);
  SET_VECTOR_ELT(out_sexp, 1, status_sexp);
  SET_VECTOR_ELT(out_sexp, 2, counts_sexp);
  for (int j = 0; j < 3; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* Eligible non-detailed DR IF chunk, part 1 (evaluate_dr_baseline_if_chunk
   with the complete-date reader up to each model's solve): the rate stage's
   slots (and endpoints/differences when sp[2] < 0 signals p_or > 0), then
   every treatment model in order: slots, endpoints and raw differences of its
   cohorts' rows, and the upstream product upstream %*% t(cross) (p_trt >= 1)
   in reference order. mdl is (2 + 2 m) x nt: gather columns of s and t, the
   model's m support positions, then its m cohorts; mdl_cross is (p_trt * m) x nt.
   Returns the gathered endpoints, the products of the models whose reads passed
   and whether they are bitwise identical. */
SEXP C_gmm_cp08_dr_if_prepare(SEXP y_sexp, SEXP date_start_sexp, SEXP post_a_sexp, SEXP post_obs_sexp,
                              SEXP rows_sexp, SEXP ci_sexp, SEXP a_sexp, SEXP w_sexp, SEXP dates_sexp,
                              SEXP sup_sexp, SEXP sup_par_sexp, SEXP mdl_sexp, SEXP mdl_cross_sexp,
                              SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_dr_if_prepare";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  cp08_chunk_t g;
  cp08_chunk_args(0, fn, &g, rows_sexp, ci_sexp, a_sexp, w_sexp);
  cp08_panel_args(0, fn, &g, y_sexp, date_start_sexp, post_a_sexp, post_obs_sexp, dates_sexp);
  int nr = 0, ns = 0, nt = 0, tmp = 0, mr = 0, mcr = 0;
  cp08_dims(0, sup_sexp, INTSXP, fn, "sup", 3, 0, &nr, &ns);
  cp08_dims(0, sup_par_sexp, REALSXP, fn, "sup_par", 3, ns, &nr, &tmp);
  cp08_dims(0, mdl_sexp, INTSXP, fn, "mdl", 0, 0, &mr, &nt);
  if (mr < 4 || (mr % 2) != 0) cp08_malformed(0, fn, "model descriptors must be (2 + 2 m) x nt");
  int m = (mr - 2) / 2;
  cp08_dims(0, mdl_cross_sexp, REALSXP, fn, "mdl_cross", 0, nt, &mcr, &tmp);
  if (mcr < m || (mcr % m) != 0) cp08_malformed(0, fn, "mdl_cross must be (p_trt * m) x nt");
  int p_trt = mcr / m;
  const int *sup = INTEGER(sup_sexp);
  const double *sp = REAL(sup_par_sexp);
  const int *mdl = INTEGER(mdl_sexp);
  const double *mc = REAL(mdl_cross_sexp);
  for (int j = 0; j < ns; ++j) {
    cp08_range(0, sup + 3 * (R_xlen_t) j + 1, 2, 1, g.nd, fn, "support gather columns outside the dates");
    if (sup[3 * (R_xlen_t) j] == NA_INTEGER) cp08_malformed(0, fn, "support cohort is NA");
  }
  for (int t = 0; t < nt; ++t) {
    const int *md = mdl + (R_xlen_t) mr * t;
    cp08_range(0, md, 2, 1, g.nd, fn, "model gather columns outside the dates");
    cp08_range(0, md + 2, m, 1, ns, fn, "model supports outside the requested supports");
    for (int k = 0; k < m; ++k) {
      if (md[2 + m + k] == NA_INTEGER) cp08_malformed(0, fn, "model cohort is NA");
    }
  }
  int B = g.B;
  uint64_t gather_len = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) g.nd, max_entries, "C_gmm_cp08_dr_if_prepare endpoint gather");
  uint64_t bp_trt = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) p_trt, max_entries, "C_gmm_cp08_dr_if_prepare model upstream product");
  uint64_t prod_len = gmm_checked_mul_u64(0, bp_trt, (uint64_t) nt, max_entries, "C_gmm_cp08_dr_if_prepare upstream products");

  SEXP gather_sexp = PROTECT(Rf_allocMatrix(REALSXP, B, g.nd));
  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP same_sexp = PROTECT(Rf_allocVector(LGLSXP, 1));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 4));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 5));
  int nprot = 6;
  g.Y = REAL(gather_sexp);
  if (gmm_buffers_overlap(g.Y, (R_xlen_t) gather_len, g.y, g.ny)) {
    gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_cp08_dr_if_prepare`.", "cdid_gmm_bound_error");
  }
  g.state = (unsigned char *) R_alloc((size_t) gather_len, sizeof(unsigned char));
  for (R_xlen_t idx = 0; idx < (R_xlen_t) gather_len; ++idx) {
    g.Y[idx] = R_NaN;
    g.state[idx] = 0;
  }
  double *prods = (double *) R_alloc((size_t) prod_len, sizeof(double));
  int *st = INTEGER(status_sexp);
  cp08_status(st, CP08_OK, 0, 0, 0);

  /* 1. Rate stage: slots (and reads when p_or > 0) of every requested support, in order. */
  int core_read = (ns > 0 && sp[2] < 0.0) ? 1 : 0;
  for (int j = 0; j < ns && st[0] == CP08_OK; ++j) {
    R_CheckUserInterrupt();
    const int *s = sup + 3 * (R_xlen_t) j;
    int code = CP08_OK;
    int row = cp08_pair_rows(&g, s, 1, s[1] - 1, s[2] - 1, core_read, &code);
    if (row >= 0) cp08_status(st, code, 1, j + 1, row + 1);
  }
  /* 2. Treatment models in order. */
  int n_ok = 0;
  for (int t = 0; t < nt && st[0] == CP08_OK; ++t) {
    R_CheckUserInterrupt();
    const int *md = mdl + (R_xlen_t) mr * t;
    int code = CP08_OK;
    int row = cp08_pair_rows(&g, md + 2 + m, m, md[0] - 1, md[1] - 1, 1, &code);
    if (row >= 0) {
      cp08_status(st, code, 2, t + 1, row + 1);
      break;
    }
    double *pc = prods + (R_xlen_t) t * bp_trt;
    for (int jp = 0; jp < p_trt && st[0] == CP08_OK; ++jp) {
      double *pc_col = pc + (R_xlen_t) jp * B;
      const double *mc_col = mc + (R_xlen_t) mcr * t + (R_xlen_t) m * jp;
      for (int i = 0; i < B; ++i) {
        double z = 0.0;
        int nz = 0;
        for (int k = 0; k < m; ++k) {
          const int *s = sup + 3 * (R_xlen_t) (md[2 + k] - 1);
          const double *p = sp + 3 * (R_xlen_t) (md[2 + k] - 1);
          z = cp08_term(z, mc_col[k], cp08_up(&g, i, s[0], p[0], p[1]), &nz);
        }
        if (nz > 1) {
          cp08_status(st, CP08_PRODUCT, 2, t + 1, i + 1);
          break;
        }
        pc_col[i] = z;
      }
    }
    if (st[0] != CP08_OK) break;
    ++n_ok;
  }
  /* 3. Bitwise identity of every passed model's product with the first. */
  int same = 1;
  for (int t = 1; t < n_ok && same; ++t) {
    if (memcmp(prods, prods + (R_xlen_t) t * bp_trt, sizeof(double) * (size_t) bp_trt) != 0) same = 0;
  }
  int n_cols = (n_ok < 1 ? 1 : (same ? 1 : n_ok)) * p_trt;
  SEXP prod_sexp = PROTECT(Rf_allocMatrix(REALSXP, B, n_cols));
  ++nprot;
  double *po = REAL(prod_sexp);
  for (R_xlen_t idx = 0; idx < (R_xlen_t) B * n_cols; ++idx) po[idx] = n_ok < 1 ? R_NaN : prods[idx];
  LOGICAL(same_sexp)[0] = same;
  double *counts = REAL(counts_sexp);
  counts[0] = g.slot_checks;
  counts[1] = g.reads;
  counts[2] = (double) g.nd;
  counts[3] = (double) n_ok;
  static const char *out_names[5] = {"status", "gather", "upstream", "identical", "counts"};
  SET_VECTOR_ELT(out_sexp, 0, status_sexp);
  SET_VECTOR_ELT(out_sexp, 1, gather_sexp);
  SET_VECTOR_ELT(out_sexp, 2, prod_sexp);
  SET_VECTOR_ELT(out_sexp, 3, same_sexp);
  SET_VECTOR_ELT(out_sexp, 4, counts_sexp);
  for (int j = 0; j < 5; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* Eligible non-detailed DR IF chunk, part 2, after R's unchanged cp04_solve
   of beta_score + upstream %*% t(cross): for the first n_ok models in order,
   the four terms psi_dir + psi_or + psi_ps + psi_sel in their addition order
   with direct = [wp (r - T) / dp, wm (r - M) / dm, -dp + wp, -dm + wm] on the
   replay rows rr (r = dy - (alpha + L_b z_i' beta)), psi_or = -xi_beta %*% d_bar_diff,
   psi_ps = -(V[2, ] %*% xi), psi_sel = upstream %*% (U[1, ] - U[2, ]), and the
   frame coordinates (direct + upstream %*% t(U)) + xi %*% t(V) whose finiteness
   is checked with xi and psi; then, when every model passed, the link columns
   in the requested order and their finiteness. mdl as in part 1; mdl_par is
   (5 + 4 * p_trt + 4 * m) x nt: alpha, T_dr, M_dr, dp, dm, V (4 x p_trt by col),
   U (4 x m by col). When p_or > 0, dimnames[[1]] holds the slope descriptor. */
SEXP C_gmm_cp08_dr_if_assemble(SEXP rows_sexp, SEXP ci_sexp, SEXP a_sexp, SEXP w_sexp, SEXP gather_sexp,
                               SEXP sup_sexp, SEXP sup_par_sexp, SEXP mdl_sexp, SEXP mdl_par_sexp,
                               SEXP xi_sexp, SEXP rr_sexp, SEXP wp_sexp, SEXP wm_sexp, SEXP lnk_sexp,
                               SEXP n_ok_sexp, SEXP dimnames_sexp, SEXP max_entries_sexp) {
  static const char *fn = "C_gmm_cp08_dr_if_assemble";
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);
  cp08_chunk_t g;
  cp08_chunk_args(0, fn, &g, rows_sexp, ci_sexp, a_sexp, w_sexp);
  int B = g.B;
  int nr = 0, nd = 0, ns = 0, nt = 0, tmp = 0, mr = 0, xc = 0, kpr = 0;
  cp08_dims(0, gather_sexp, REALSXP, fn, "gather", B, 0, &nr, &nd);
  g.Y = REAL(gather_sexp);
  g.nd = nd;
  cp08_dims(0, sup_sexp, INTSXP, fn, "sup", 3, 0, &nr, &ns);
  cp08_dims(0, sup_par_sexp, REALSXP, fn, "sup_par", 3, ns, &nr, &tmp);
  cp08_dims(0, mdl_sexp, INTSXP, fn, "mdl", 0, 0, &mr, &nt);
  if (mr < 4 || (mr % 2) != 0) cp08_malformed(0, fn, "model descriptors must be (2 + 2 m) x nt");
  int m = (mr - 2) / 2;
  cp08_dims(0, mdl_par_sexp, REALSXP, fn, "mdl_par", 0, nt, &kpr, &tmp);
  if (kpr < 9 + 4 * m || ((kpr - 5 - 4 * m) % 4) != 0) {
    cp08_malformed(0, fn, "mdl_par must be (5 + 4 * p_trt + 4 * m) x nt");
  }
  int p_trt = (kpr - 5 - 4 * m) / 4;
  int kp = kpr;
  cp08_dims(0, xi_sexp, REALSXP, fn, "xi", B, 0, &nr, &xc);
  R_xlen_t n_rr = cp08_len(0, rr_sexp, INTSXP, fn, "rr");
  if (cp08_len(0, wp_sexp, REALSXP, fn, "wp") != n_rr || cp08_len(0, wm_sexp, REALSXP, fn, "wm") != n_rr ||
      n_rr > (R_xlen_t) B) {
    cp08_malformed(0, fn, "replay rows and weights differ in length or exceed the chunk");
  }
  R_xlen_t nl_x = cp08_len(0, lnk_sexp, INTSXP, fn, "lnk");
  if (cp08_len(0, n_ok_sexp, INTSXP, fn, "n_ok") != 1 || nl_x < 1 || nl_x > INT_MAX) {
    cp08_malformed(0, fn, "n_ok must be one count and links nonempty");
  }
  int nl = (int) nl_x;
  int n_ok = INTEGER(n_ok_sexp)[0];
  if (n_ok == NA_INTEGER || n_ok < 0 || n_ok > nt || (xc != p_trt && xc != p_trt * n_ok)) {
    cp08_malformed(0, fn, "model count or xi columns outside the models");
  }
  cp08_dimnames(0, dimnames_sexp, nl, fn);
  SEXP slope_sexp = (dimnames_sexp != R_NilValue) ? VECTOR_ELT(dimnames_sexp, 0) : R_NilValue;
  const int *sup = INTEGER(sup_sexp);
  const double *sp = REAL(sup_par_sexp);
  const int *mdl = INTEGER(mdl_sexp);
  const double *mp = REAL(mdl_par_sexp);
  const double *xi = REAL(xi_sexp);
  const int *rr = INTEGER(rr_sexp);
  const double *wp = REAL(wp_sexp);
  const double *wm = REAL(wm_sexp);
  const int *lnk = INTEGER(lnk_sexp);
  for (int j = 0; j < ns; ++j) {
    if (sup[3 * (R_xlen_t) j] == NA_INTEGER) cp08_malformed(0, fn, "support cohort is NA");
  }
  for (int t = 0; t < nt; ++t) {
    const int *md = mdl + (R_xlen_t) mr * t;
    cp08_range(0, md, 2, 1, nd, fn, "model gather columns outside the gather");
    cp08_range(0, md + 2, m, 1, ns, fn, "model supports outside the requested supports");
  }
  cp08_range(0, rr, n_rr, 1, B, fn, "replay rows outside the chunk");
  for (R_xlen_t p = 1; p < n_rr; ++p) {
    if (rr[p] <= rr[p - 1]) cp08_malformed(0, fn, "replay rows must strictly increase");
  }
  cp08_range(0, lnk, nl, 1, nt, fn, "link models outside the requested models");
  uint64_t psi_len = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) nt, max_entries, "C_gmm_cp08_dr_if_assemble model IF columns");
  uint64_t phi_len = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) nl, max_entries, "C_gmm_cp08_dr_if_assemble IF (B x L)");

  SEXP phi_sexp = PROTECT(Rf_allocMatrix(REALSXP, B, nl));
  SEXP status_sexp = PROTECT(Rf_allocVector(INTSXP, 4));
  SEXP counts_sexp = PROTECT(Rf_allocVector(REALSXP, 2));
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  const int nprot = 5;
  double *phi = REAL(phi_sexp);
  if (gmm_buffers_overlap(phi, (R_xlen_t) phi_len, xi, (R_xlen_t) B * xc) ||
      gmm_buffers_overlap(phi, (R_xlen_t) phi_len, g.Y, (R_xlen_t) B * nd)) {
    gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_cp08_dr_if_assemble`.", "cdid_gmm_bound_error");
  }
  cp08_set_dimnames(phi_sexp, dimnames_sexp);
  for (R_xlen_t idx = 0; idx < (R_xlen_t) phi_len; ++idx) phi[idx] = 0.0;
  double *psi = (double *) R_alloc((size_t) psi_len, sizeof(double));
  double *upv = (double *) R_alloc((size_t) m, sizeof(double));
  int *st = INTEGER(status_sexp);
  cp08_status(st, CP08_OK, 0, 0, 0);

  cp08_slope_t sl;
  sl.p = 0;
  const double *mdl_or_par = NULL;
  if (slope_sexp != R_NilValue) {
    if (cp08_or_slope_if(nprot, fn, &g, ns, sup, sp, slope_sexp, max_entries, st, &sl)) {
      int r1 = 0, c1 = 0;
      cp08_dims(nprot, sl.tail_sexp, REALSXP, fn, "mdl_or_par", 1 + sl.p, nt, &r1, &c1);
      mdl_or_par = REAL(sl.tail_sexp);
    }
  }

  int done = 0;
  for (int t = 0; t < n_ok && st[0] == CP08_OK; ++t) {
    R_CheckUserInterrupt();
    const int *md = mdl + (R_xlen_t) mr * t;
    const double *par = mp + (R_xlen_t) kp * t;
    double alpha = par[0], T = par[1], M = par[2], dp = par[3], dm = par[4];
    const double *V = par + 5;
    const double *U = par + 5 + 4 * p_trt;
    const double *xc_t = xi + (R_xlen_t) (xc == p_trt ? 0 : t * p_trt) * B;
    double L_t = (sl.p > 0) ? mdl_or_par[(R_xlen_t) (1 + sl.p) * t] : 0.0;
    const double *d_bar_diff = (sl.p > 0) ? (mdl_or_par + (R_xlen_t) (1 + sl.p) * t + 1) : NULL;
    double *out = psi + (R_xlen_t) t * B;
    R_xlen_t p = 0;
    int bad_row = -1;
    for (int i = 0; i < B && st[0] == CP08_OK; ++i) {
      double d1 = 0.0, d2 = 0.0, d3 = gmm_fp_neg(dp), d4 = gmm_fp_neg(dm);
      if (p < n_rr && rr[p] - 1 == i) {
        double dy = cp08_dy(&g, i, md[0] - 1, md[1] - 1);
        if (!R_FINITE(dy)) {
          cp08_status(st, CP08_GATHER, 2, t + 1, i + 1);
          break;
        }
        double db = 0.0;
        if (sl.p > 0) {
          int r_idx = g.rows[i] - 1;
          int a = g.a[r_idx];
          for (int k = 0; k < sl.p; ++k) {
            double dk = gmm_fp_mul(L_t, cp08_z_elem(&sl, a, k));
            db = gmm_fp_add(db, gmm_fp_mul(dk, sl.beta[k]));
          }
        }
        double m_hat = gmm_fp_add(alpha, db);
        double r = gmm_fp_sub(dy, m_hat);
        d1 = gmm_fp_div(gmm_fp_mul(wp[p], gmm_fp_sub(r, T)), dp);
        d2 = gmm_fp_div(gmm_fp_mul(wm[p], gmm_fp_sub(r, M)), dm);
        d3 = gmm_fp_add(d3, wp[p]);
        d4 = gmm_fp_add(d4, wm[p]);
        ++p;
      }
      int x_finite = 1;
      for (int jp = 0; jp < p_trt; ++jp) {
        if (!R_FINITE(xc_t[(R_xlen_t) jp * B + i])) x_finite = 0;
      }
      for (int k = 0; k < m; ++k) {
        const int *s = sup + 3 * (R_xlen_t) (md[2 + k] - 1);
        const double *sv = sp + 3 * (R_xlen_t) (md[2 + k] - 1);
        upv[k] = cp08_up(&g, i, s[0], sv[0], sv[1]);
      }
      int nz = 0;
      double psi_dir = gmm_fp_sub(d1, d2);
      double psi_or = 0.0;
      if (sl.p > 0) {
        double dot_or = 0.0;
        for (int k = 0; k < sl.p; ++k) {
          dot_or = gmm_fp_add(dot_or, gmm_fp_mul(sl.xi_beta[(R_xlen_t) k * B + i], d_bar_diff[k]));
        }
        psi_or = gmm_fp_neg(dot_or);
      }
      double dot_ps = 0.0;
      for (int jp = 0; jp < p_trt; ++jp) {
        dot_ps = gmm_fp_add(dot_ps, gmm_fp_mul(xc_t[(R_xlen_t) jp * B + i], V[4 * (R_xlen_t) jp + 1]));
      }
      double psi_ps = gmm_fp_neg(dot_ps);
      double psi_sel = 0.0;
      for (int k = 0; k < m; ++k) {
        psi_sel = cp08_term(psi_sel, gmm_fp_sub(U[4 * (R_xlen_t) k], U[4 * (R_xlen_t) k + 1]), upv[k], &nz);
      }
      if (nz > 1) {
        cp08_status(st, CP08_PRODUCT, 2, t + 1, i + 1);
        break;
      }
      double psi_dr = gmm_fp_add(gmm_fp_add(gmm_fp_add(psi_dir, psi_or), psi_ps), psi_sel);
      double dcol[4] = {d1, d2, d3, d4};
      int finite = x_finite && R_FINITE(psi_dr);
      for (int cc = 0; cc < 4; ++cc) {
        double z = 0.0;
        nz = 0;
        for (int k = 0; k < m; ++k) z = cp08_term(z, U[4 * (R_xlen_t) k + cc], upv[k], &nz);
        if (nz > 1) {
          cp08_status(st, CP08_PRODUCT, 2, t + 1, i + 1);
          break;
        }
        double vx = 0.0;
        for (int jp = 0; jp < p_trt; ++jp) {
          vx = gmm_fp_add(vx, gmm_fp_mul(xc_t[(R_xlen_t) jp * B + i], V[4 * (R_xlen_t) jp + cc]));
        }
        double coord = gmm_fp_add(gmm_fp_add(dcol[cc], z), vx);
        if (!R_FINITE(coord)) finite = 0;
      }
      if (st[0] != CP08_OK) break;
      out[i] = psi_dr;
      if (!finite && bad_row < 0) bad_row = i;
    }
    if (st[0] != CP08_OK) break;
    if (bad_row >= 0) {
      cp08_status(st, CP08_DR_FRAME_IF, 2, t + 1, bad_row + 1);
      break;
    }
    ++done;
  }
  if (st[0] == CP08_OK && done == nt) {
    for (int j = 0; j < nl; ++j) {
      const double *src = psi + (R_xlen_t) (lnk[j] - 1) * B;
      double *dst = phi + (R_xlen_t) j * B;
      for (int i = 0; i < B; ++i) dst[i] = src[i];
    }
    for (int j = 0; j < nl && st[0] == CP08_OK; ++j) {
      const double *col = phi + (R_xlen_t) j * B;
      for (int i = 0; i < B; ++i) {
        if (!R_FINITE(col[i])) {
          cp08_status(st, CP08_DR_LINK_IF, 3, j + 1, i + 1);
          break;
        }
      }
    }
  }
  double *counts = REAL(counts_sexp);
  counts[0] = (double) done;
  counts[1] = (st[0] == CP08_OK && done == nt) ? (double) nl : 0.0;
  static const char *out_names[3] = {"phi", "status", "counts"};
  SET_VECTOR_ELT(out_sexp, 0, phi_sexp);
  SET_VECTOR_ELT(out_sexp, 1, status_sexp);
  SET_VECTOR_ELT(out_sexp, 2, counts_sexp);
  for (int j = 0; j < 3; ++j) SET_STRING_ELT(names_sexp, j, Rf_mkChar(out_names[j]));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* Bounded unit-row IF projection batch accumulator for the unweighted
   nonclustered unit-row branch of `project_gmm_if_stream` (`R/chained_gmm_core.R`).
   Replaces the per-row R `flush` slice/vector/tcrossprod/S-copy loops across both
   complete producer passes while preserving `phi %*% Lambda` in R, exact
   projected-row/ID/weight/batch replay checks, row-by-row addition and refusal
   order (`nout > no`, row finiteness, duplicate/metadata flag, then scalar
   `sum_w + w`, `weighted_sum + w * vec`, `raw_sum + vec` in pass 1 and
   `z = vec - psi_mean`, `S + w * (z_j * z_k)` in pass 2). Inputs are never
   modified in place. All simultaneously live batch/state/returned buffers are
   reserved against explicit `max_matrix_entries` before allocation. */
SEXP C_gmm_if_unit_batch(
    SEXP pass_sexp, SEXP psi_sexp, SEXP w_sexp, SEXP flags_sexp, SEXP counts_sexp,
    SEXP weighted_sum_sexp, SEXP raw_sum_sexp, SEXP psi_mean_sexp, SEXP S_sexp,
    SEXP max_entries_sexp) {
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  if (max_entries_sexp == R_NilValue) {
    gmm_c_abort(0, "`max_matrix_entries` must be an explicit finite integer bound in [1, 2^53 - 1].",
                "cdid_gmm_bound_error");
  }
  double max_entries = gmm_parse_max_entries(0, max_entries_sexp);

  if (pass_sexp == R_NilValue || Rf_isObject(pass_sexp) || TYPEOF(pass_sexp) != INTSXP ||
      Rf_getAttrib(pass_sexp, R_DimSymbol) != R_NilValue || XLENGTH(pass_sexp) != 1 ||
      INTEGER(pass_sexp)[0] == NA_INTEGER || (INTEGER(pass_sexp)[0] != 1 && INTEGER(pass_sexp)[0] != 2)) {
    gmm_c_abort(0, "`pass` in `C_gmm_if_unit_batch` must be integer 1L or 2L.", "cdid_gmm_bound_error");
  }
  int pass = INTEGER(pass_sexp)[0];

  int B = 0, J = 0;
  R_xlen_t psi_len = 0;
  const double *psi_unused = NULL;
  gmm_validate_real_matrix_header(0, psi_sexp, "psi", max_entries, &B, &J, &psi_len, &psi_unused);

  R_xlen_t w_len = 0;
  const double *w_unused = NULL;
  gmm_validate_real_vector_header(0, w_sexp, "w", max_entries, &w_len, &w_unused);
  if (w_len != (R_xlen_t) B) {
    gmm_c_abort(0, "`w` length in `C_gmm_if_unit_batch` must equal `nrow(psi)`.", "cdid_gmm_bound_error");
  }

  if (flags_sexp == R_NilValue || Rf_isObject(flags_sexp) || TYPEOF(flags_sexp) != LGLSXP ||
      Rf_getAttrib(flags_sexp, R_DimSymbol) != R_NilValue) {
    gmm_c_abort(0, "`flags` in `C_gmm_if_unit_batch` must be a non-object logical vector.",
                "cdid_gmm_accuracy_error");
  }
  if (XLENGTH(flags_sexp) != (R_xlen_t) B) {
    gmm_c_abort(0, "`flags` length in `C_gmm_if_unit_batch` must equal `nrow(psi)`.", "cdid_gmm_bound_error");
  }
  (void) gmm_checked_mul_u64(0, (uint64_t) B, 1ULL, max_entries, "flags");

  R_xlen_t counts_len = 0;
  const double *counts_unused = NULL;
  gmm_validate_real_vector_header(0, counts_sexp, "counts", max_entries, &counts_len, &counts_unused);
  if (counts_len != 4) {
    gmm_c_abort(0, "`counts` in `C_gmm_if_unit_batch` must have length 4 (`c(nout, no, npos, sum_w)`).",
                "cdid_gmm_bound_error");
  }

  uint64_t bj = gmm_checked_mul_u64(0, (uint64_t) B, (uint64_t) J, max_entries, "unit batch psi (B x J)");
  uint64_t b2 = gmm_checked_mul_u64(0, 2ULL, (uint64_t) B, max_entries, "unit batch weights/flags (2 x B)");

  if (pass == 1) {
    if (psi_mean_sexp != R_NilValue || S_sexp != R_NilValue) {
      gmm_c_abort(0, "Pass 1 of `C_gmm_if_unit_batch` requires `psi_mean` and `S` to be NULL.",
                  "cdid_gmm_bound_error");
    }
    R_xlen_t ws_len = 0, rs_len = 0;
    const double *ws_unused = NULL, *rs_unused = NULL;
    gmm_validate_real_vector_header(0, weighted_sum_sexp, "weighted_sum", max_entries, &ws_len, &ws_unused);
    gmm_validate_real_vector_header(0, raw_sum_sexp, "raw_sum", max_entries, &rs_len, &rs_unused);
    if (ws_len != (R_xlen_t) J || rs_len != (R_xlen_t) J) {
      gmm_c_abort(0, "`weighted_sum` and `raw_sum` lengths in `C_gmm_if_unit_batch` must equal `ncol(psi)`.",
                  "cdid_gmm_bound_error");
    }
    uint64_t j4 = gmm_checked_mul_u64(0, 4ULL, (uint64_t) J, max_entries, "unit batch pass-1 sums (4 x J)");
    (void) gmm_checked_mul_u64(0, 1ULL, bj + b2 + j4 + 8ULL, max_entries,
                               "unit batch pass-1 simultaneous buffers");
  } else {
    if (weighted_sum_sexp != R_NilValue || raw_sum_sexp != R_NilValue) {
      gmm_c_abort(0, "Pass 2 of `C_gmm_if_unit_batch` requires `weighted_sum` and `raw_sum` to be NULL.",
                  "cdid_gmm_bound_error");
    }
    R_xlen_t pm_len = 0, S_len_hdr = 0;
    int S_nr = 0, S_nc = 0;
    const double *pm_unused = NULL, *S_unused = NULL;
    gmm_validate_real_vector_header(0, psi_mean_sexp, "psi_mean", max_entries, &pm_len, &pm_unused);
    if (pm_len != (R_xlen_t) J) {
      gmm_c_abort(0, "`psi_mean` length in `C_gmm_if_unit_batch` must equal `ncol(psi)`.",
                  "cdid_gmm_bound_error");
    }
    gmm_validate_real_matrix_header(0, S_sexp, "S", max_entries, &S_nr, &S_nc, &S_len_hdr, &S_unused);
    if (S_nr != J || S_nc != J) {
      gmm_c_abort(0, "`S` dimensions in `C_gmm_if_unit_batch` must be `ncol(psi) x ncol(psi)`.",
                  "cdid_gmm_bound_error");
    }
    uint64_t j2 = gmm_checked_mul_u64(0, 2ULL, (uint64_t) J, max_entries, "unit batch pass-2 centering (2 x J)");
    uint64_t jj = gmm_checked_mul_u64(0, (uint64_t) J, (uint64_t) J, max_entries, "unit batch covariance (J x J)");
    uint64_t jj2 = gmm_checked_mul_u64(0, 2ULL, jj, max_entries, "unit batch pass-2 covariance buffers (2 x J x J)");
    (void) gmm_checked_mul_u64(0, 1ULL, bj + b2 + j2 + jj2 + 8ULL, max_entries,
                               "unit batch pass-2 simultaneous buffers");
  }

  const double *psi = REAL(psi_sexp);
  const double *w = REAL(w_sexp);
  for (int i = 0; i < B; ++i) {
    if (!R_FINITE(w[i]) || w[i] < 0.0) {
      gmm_c_abort(0, "Invalid replay weights.", "cdid_gmm_replay_error");
    }
  }
  const int *flags = LOGICAL(flags_sexp);
  for (int i = 0; i < B; ++i) {
    if (flags[i] == NA_LOGICAL) {
      gmm_c_abort(0, "`flags` in `C_gmm_if_unit_batch` cannot contain NA.", "cdid_gmm_bound_error");
    }
  }
  const double *counts_in = REAL(counts_sexp);
  double d_nout = counts_in[0];
  double d_no = counts_in[1];
  double d_npos = counts_in[2];
  double sum_w = counts_in[3];
  if (!R_FINITE(d_nout) || d_nout != floor(d_nout) || d_nout < 0.0 ||
      !R_FINITE(d_no) || d_no != floor(d_no) || d_no < 1.0 || d_no > (double) INT_MAX ||
      d_nout > d_no ||
      !R_FINITE(d_npos) || d_npos != floor(d_npos) || d_npos < 0.0 ||
      d_npos > (pass == 1 ? d_nout : d_no) ||
      !R_FINITE(sum_w) || sum_w < 0.0) {
    gmm_c_abort(0, "Invalid `counts` state in `C_gmm_if_unit_batch`.", "cdid_gmm_bound_error");
  }
  int nout = (int) d_nout;
  int no = (int) d_no;
  int npos = (int) d_npos;

  if (pass == 1) {
    const double *ws_in = REAL(weighted_sum_sexp);
    const double *rs_in = REAL(raw_sum_sexp);
    for (int j = 0; j < J; ++j) {
      if (!R_FINITE(ws_in[j]) || !R_FINITE(rs_in[j])) {
        gmm_c_abort(0, "Replay mean/weight accumulation overflowed.", "cdid_gmm_replay_error");
      }
    }

    SEXP counts_out_sexp = PROTECT(Rf_allocVector(REALSXP, 4));
    SEXP ws_out_sexp = PROTECT(Rf_allocVector(REALSXP, J));
    SEXP rs_out_sexp = PROTECT(Rf_allocVector(REALSXP, J));
    SEXP psi_dn = Rf_getAttrib(psi_sexp, R_DimNamesSymbol);
    SEXP col_nms = (psi_dn != R_NilValue && TYPEOF(psi_dn) == VECSXP && XLENGTH(psi_dn) == 2)
                       ? VECTOR_ELT(psi_dn, 1)
                       : R_NilValue;
    SEXP ws_nms = Rf_getAttrib(weighted_sum_sexp, R_NamesSymbol);
    SEXP rs_nms = Rf_getAttrib(raw_sum_sexp, R_NamesSymbol);
    if (ws_nms != R_NilValue) {
      Rf_setAttrib(ws_out_sexp, R_NamesSymbol, ws_nms);
    } else if (col_nms != R_NilValue && TYPEOF(col_nms) == STRSXP && XLENGTH(col_nms) == (R_xlen_t) J) {
      Rf_setAttrib(ws_out_sexp, R_NamesSymbol, col_nms);
    }
    if (rs_nms != R_NilValue) {
      Rf_setAttrib(rs_out_sexp, R_NamesSymbol, rs_nms);
    } else if (col_nms != R_NilValue && TYPEOF(col_nms) == STRSXP && XLENGTH(col_nms) == (R_xlen_t) J) {
      Rf_setAttrib(rs_out_sexp, R_NamesSymbol, col_nms);
    }
    SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
    const int nprot = 5;
    double *counts_out = REAL(counts_out_sexp);
    double *ws_out = REAL(ws_out_sexp);
    double *rs_out = REAL(rs_out_sexp);
    if (gmm_buffers_overlap(counts_out, 4, psi, psi_len) ||
        gmm_buffers_overlap(counts_out, 4, w, B) ||
        gmm_buffers_overlap(counts_out, 4, counts_in, 4) ||
        gmm_buffers_overlap(counts_out, 4, ws_in, J) ||
        gmm_buffers_overlap(counts_out, 4, rs_in, J) ||
        gmm_buffers_overlap(ws_out, J, psi, psi_len) ||
        gmm_buffers_overlap(ws_out, J, w, B) ||
        gmm_buffers_overlap(ws_out, J, counts_in, 4) ||
        gmm_buffers_overlap(ws_out, J, ws_in, J) ||
        gmm_buffers_overlap(ws_out, J, rs_in, J) ||
        gmm_buffers_overlap(rs_out, J, psi, psi_len) ||
        gmm_buffers_overlap(rs_out, J, w, B) ||
        gmm_buffers_overlap(rs_out, J, counts_in, 4) ||
        gmm_buffers_overlap(rs_out, J, ws_in, J) ||
        gmm_buffers_overlap(rs_out, J, rs_in, J) ||
        gmm_buffers_overlap(ws_out, J, rs_out, J)) {
      gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_if_unit_batch`.", "cdid_gmm_bound_error");
    }
    for (int j = 0; j < J; ++j) {
      ws_out[j] = ws_in[j];
      rs_out[j] = rs_in[j];
    }
    for (int i = 0; i < B; ++i) {
      if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
      if (nout >= no) {
        gmm_c_abort(nprot, "Replay exceeds declared output count.", "cdid_gmm_replay_error");
      }
      ++nout;
      for (int j = 0; j < J; ++j) {
        if (!R_FINITE(psi[i + (R_xlen_t) j * B])) {
          gmm_c_abort(nprot, "Cluster sum overflowed.", "cdid_gmm_replay_error");
        }
      }
      if (flags[i] != 0) {
        gmm_c_abort(nprot, "Closed cluster reappeared.", "cdid_gmm_replay_error");
      }
      double wi = w[i];
      if (wi > 0.0) ++npos;
      sum_w = gmm_fp_add(sum_w, wi);
      int ok = R_FINITE(sum_w);
      for (int j = 0; j < J; ++j) {
        double v = psi[i + (R_xlen_t) j * B];
        double ws = gmm_fp_add(ws_out[j], gmm_fp_mul(wi, v));
        double rs = gmm_fp_add(rs_out[j], v);
        ws_out[j] = ws;
        rs_out[j] = rs;
        if (!R_FINITE(ws) || !R_FINITE(rs)) ok = 0;
      }
      if (!ok) {
        gmm_c_abort(nprot, "Replay mean/weight accumulation overflowed.", "cdid_gmm_replay_error");
      }
    }
    counts_out[0] = (double) nout;
    counts_out[1] = (double) no;
    counts_out[2] = (double) npos;
    counts_out[3] = sum_w;
    SET_STRING_ELT(names_sexp, 0, Rf_mkChar("counts"));
    SET_STRING_ELT(names_sexp, 1, Rf_mkChar("weighted_sum"));
    SET_STRING_ELT(names_sexp, 2, Rf_mkChar("raw_sum"));
    SET_VECTOR_ELT(out_sexp, 0, counts_out_sexp);
    SET_VECTOR_ELT(out_sexp, 1, ws_out_sexp);
    SET_VECTOR_ELT(out_sexp, 2, rs_out_sexp);
    Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
    UNPROTECT(nprot);
    return out_sexp;
  }

  R_xlen_t S_len = (R_xlen_t) J * (R_xlen_t) J;
  const double *pm_in = REAL(psi_mean_sexp);
  const double *S_in = REAL(S_sexp);
  for (int j = 0; j < J; ++j) {
    if (!R_FINITE(pm_in[j])) {
      gmm_c_abort(0, "Replay mean overflowed.", "cdid_gmm_replay_error");
    }
  }
  for (R_xlen_t idx = 0; idx < S_len; ++idx) {
    if (!R_FINITE(S_in[idx])) {
      gmm_c_abort(0, "Replay covariance overflowed.", "cdid_gmm_replay_error");
    }
  }

  SEXP counts_out_sexp = PROTECT(Rf_allocVector(REALSXP, 4));
  SEXP S_out_sexp = PROTECT(Rf_allocMatrix(REALSXP, J, J));
  SEXP S_dn = Rf_getAttrib(S_sexp, R_DimNamesSymbol);
  if (S_dn != R_NilValue) {
    Rf_setAttrib(S_out_sexp, R_DimNamesSymbol, S_dn);
  }
  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
  const int nprot = 4;
  double *counts_out = REAL(counts_out_sexp);
  double *S_out = REAL(S_out_sexp);
  double *z = (double *) R_alloc((size_t) J, sizeof(double));
  if (gmm_buffers_overlap(counts_out, 4, psi, psi_len) ||
      gmm_buffers_overlap(counts_out, 4, w, B) ||
      gmm_buffers_overlap(counts_out, 4, counts_in, 4) ||
      gmm_buffers_overlap(counts_out, 4, pm_in, J) ||
      gmm_buffers_overlap(counts_out, 4, S_in, S_len) ||
      gmm_buffers_overlap(S_out, S_len, psi, psi_len) ||
      gmm_buffers_overlap(S_out, S_len, w, B) ||
      gmm_buffers_overlap(S_out, S_len, counts_in, 4) ||
      gmm_buffers_overlap(S_out, S_len, pm_in, J) ||
      gmm_buffers_overlap(S_out, S_len, S_in, S_len)) {
    gmm_c_abort(nprot, "Allocated output aliases input in `C_gmm_if_unit_batch`.", "cdid_gmm_bound_error");
  }
  for (R_xlen_t idx = 0; idx < S_len; ++idx) {
    S_out[idx] = S_in[idx];
  }
  for (int i = 0; i < B; ++i) {
    if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
    if (nout >= no) {
      gmm_c_abort(nprot, "Replay exceeds declared output count.", "cdid_gmm_replay_error");
    }
    ++nout;
    for (int j = 0; j < J; ++j) {
      if (!R_FINITE(psi[i + (R_xlen_t) j * B])) {
        gmm_c_abort(nprot, "Cluster sum overflowed.", "cdid_gmm_replay_error");
      }
    }
    if (flags[i] != 0) {
      gmm_c_abort(nprot, "Replay output metadata changed.", "cdid_gmm_replay_error");
    }
    double wi = w[i];
    for (int j = 0; j < J; ++j) {
      z[j] = gmm_fp_sub(psi[i + (R_xlen_t) j * B], pm_in[j]);
    }
    int ok = 1;
    for (int k = 0; k < J; ++k) {
      double zk = z[k];
      double *S_col = S_out + (R_xlen_t) k * J;
      for (int j = 0; j < J; ++j) {
        double zjk = gmm_fp_mul(z[j], zk);
        double wzjk = gmm_fp_mul(wi, zjk);
        double s_new = gmm_fp_add(S_col[j], wzjk);
        S_col[j] = s_new;
        if (!R_FINITE(s_new)) ok = 0;
      }
    }
    if (!ok) {
      gmm_c_abort(nprot, "Replay covariance overflowed.", "cdid_gmm_replay_error");
    }
  }
  counts_out[0] = (double) nout;
  counts_out[1] = (double) no;
  counts_out[2] = (double) npos;
  counts_out[3] = sum_w;
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("counts"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("S"));
  SET_VECTOR_ELT(out_sexp, 0, counts_out_sexp);
  SET_VECTOR_ELT(out_sexp, 1, S_out_sexp);
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(nprot);
  return out_sexp;
}

/* ---------------------------------------------------------------------------
   CP08 payload accounting: the exact charge of the recursive R walk
   cp04_payload(x, scan_attrs) in R/chained_cp04_corollary1.R, without building
   the equivalent attributes()/list recursion graph in R.
     * environments and NULL charge nothing (and their attributes are not read);
     * every other node charges length(x) + 8, where length(x) is R's length()
       (base xlength, or the dispatched method when x carries the object bit);
     * lists and pairlists (is.list) recurse into every element with attribute
       scanning on; every alias occurrence is charged again, never deduplicated;
     * character vectors add ceil(bytes / 2) entries and the bytes of each
       non-NA string (nchar(type = "bytes"); NA_character_ counts zero);
     * when scan_attrs is on and attributes(x) is non-NULL, that attribute list
       is charged ONCE as a list node (its own names are not scanned) and each
       attribute value is walked with attribute scanning on. The values are the
       ones attributes(x) exposes: compact data.frame row names expand to
       seq_len(n), tagged pairlists expose their tag names, one-dimensional
       arrays expose dimnames[[1]] as names.
   Counts accumulate in exact 64-bit integers and are refused beyond 2^53 - 1;
   the interpreter stack and user interrupts are checked while walking.
   Borrowed inputs are never written. */
typedef struct {
  uint64_t entries;
  uint64_t char_bytes;
  uint64_t nodes;
  SEXP call_env;
  SEXP scan_attrs;
} cp08_payload_acc_t;

static void cp08_payload_add(cp08_payload_acc_t *acc, uint64_t entries, uint64_t bytes) {
  if (entries > GMM_MAX_EXACT_INT - acc->entries || bytes > GMM_MAX_EXACT_INT - acc->char_bytes) {
    gmm_c_abort(0, "Payload count overflows exact IEEE binary64 integer range (2^53 - 1).", "cdid_gmm_bound_error");
  }
  acc->entries += entries;
  acc->char_bytes += bytes;
}

static void cp08_payload_tick(cp08_payload_acc_t *acc) {
  acc->nodes += 1ULL;
  if ((acc->nodes % 4096ULL) == 0ULL) R_CheckUserInterrupt();
}

/* R-level length(x): base xlength(), or the dispatched method for an object
   (evaluated on the quoted node so no language object is ever executed). */
static uint64_t cp08_payload_length(SEXP x, SEXP call_env) {
  if (!Rf_isObject(x)) return (uint64_t) Rf_xlength(x);
  SEXP quoted = PROTECT(Rf_lang2(Rf_install("quote"), x));
  SEXP call = PROTECT(Rf_lang2(Rf_install("length"), quoted));
  SEXP val = PROTECT(Rf_eval(call, call_env));
  double len = -1.0;
  if (TYPEOF(val) == INTSXP && XLENGTH(val) == 1) {
    len = (INTEGER(val)[0] == NA_INTEGER) ? -1.0 : (double) INTEGER(val)[0];
  } else if (TYPEOF(val) == REALSXP && XLENGTH(val) == 1) {
    len = REAL(val)[0];
  }
  UNPROTECT(3);
  if (!R_FINITE(len) || len < 0.0 || len != floor(len) || len > (double) GMM_MAX_EXACT_INT) {
    gmm_c_abort(0, "A `length()` method returned an invalid count during payload accounting.", "cdid_gmm_bound_error");
  }
  return (uint64_t) len;
}

static void cp08_payload_walk(SEXP x, int scan_attrs, cp08_payload_acc_t *acc);

/* The names vector attributes() derives from the tags of a pairlist/call:
   one string per node (blank for untagged nodes), no attributes of its own. */
static void cp08_payload_tag_names(SEXP x, cp08_payload_acc_t *acc) {
  uint64_t len = 0, entries = 0, bytes = 0;
  for (SEXP p = x; p != R_NilValue; p = CDR(p)) {
    len += 1ULL;
    if (TYPEOF(TAG(p)) == SYMSXP) {
      uint64_t nb = (uint64_t) LENGTH(PRINTNAME(TAG(p)));
      bytes += nb;
      entries += (nb + 1ULL) / 2ULL;
    }
  }
  cp08_payload_tick(acc);
  cp08_payload_add(acc, len + 8ULL, 0ULL);
  cp08_payload_add(acc, entries, bytes);
}

static int cp08_payload_any_tag(SEXP x) {
  for (SEXP p = x; p != R_NilValue; p = CDR(p)) {
    if (TYPEOF(TAG(p)) == SYMSXP) return 1;
  }
  return 0;
}

typedef SEXP (*cp08_attrib_visitor_fn)(SEXP tag, SEXP value, void *data);

static int cp08_payload_has_raw_attrib(SEXP x) {
#if R_VERSION >= R_Version(4, 5, 0)
  return ANY_ATTRIB(x) != 0;
#else
  return ATTRIB(x) != R_NilValue;
#endif
}

static SEXP cp08_payload_each_attrib(SEXP x, cp08_attrib_visitor_fn cb, void *data) {
  PROTECT(x);
#if R_VERSION >= R_Version(4, 6, 0)
  SEXP out = R_mapAttrib(x, cb, data);
#else
  SEXP out = NULL;
  for (SEXP a = ATTRIB(x); a != R_NilValue; a = CDR(a)) {
    PROTECT(a);
    SEXP tag = PROTECT(TAG(a));
    SEXP value = PROTECT(CAR(a));
    out = cb(tag, value, data);
    UNPROTECT(3);
    if (out != NULL) break;
  }
#endif
  UNPROTECT(1);
  return out;
}

static SEXP cp08_payload_find_attrib_cb(SEXP tag, SEXP value, void *data) {
  SEXP name = (SEXP) data;
  return (tag == name) ? value : NULL;
}

static SEXP cp08_payload_attrib_node(SEXP x, SEXP name) {
  SEXP found = cp08_payload_each_attrib(x, cp08_payload_find_attrib_cb, (void *) name);
  return (found != NULL) ? found : R_NilValue;
}

/* One attribute value exactly as do_attributes() obtains it through getAttrib(). */
static void cp08_payload_attr_value(SEXP x, SEXP tag, SEXP value, cp08_payload_acc_t *acc) {
  if (TYPEOF(tag) != SYMSXP) {
    cp08_payload_walk(value, 1, acc);
    return;
  }
  if (tag == R_RowNamesSymbol && TYPEOF(value) == INTSXP && !Rf_inherits(value, "factor") &&
      XLENGTH(value) == 2 && INTEGER(value)[0] == NA_INTEGER) {
    /* compact c(NA, n) row names: getAttrib() returns seq_len(|n|) with no attributes */
    int n = INTEGER(value)[1];
    uint64_t len = 0ULL;
    if (n != NA_INTEGER && n != 0) len = (uint64_t) (n < 0 ? -(int64_t) n : (int64_t) n);
    cp08_payload_tick(acc);
    cp08_payload_add(acc, len + 8ULL, 0ULL);
    return;
  }
  if (tag == R_NamesSymbol) {
    int type = TYPEOF(x);
    int vector_like = Rf_isVector(x) || type == LISTSXP || type == LANGSXP;
    if (vector_like) {
      SEXP dim = cp08_payload_attrib_node(x, R_DimSymbol);
      if (TYPEOF(dim) == INTSXP && XLENGTH(dim) == 1) {
        SEXP dn = cp08_payload_attrib_node(x, R_DimNamesSymbol);
        if (dn != R_NilValue && TYPEOF(dn) == VECSXP && XLENGTH(dn) >= 1) {
          cp08_payload_walk(VECTOR_ELT(dn, 0), 1, acc);
          return;
        }
      }
    }
    if (type == LISTSXP || type == LANGSXP || type == DOTSXP) {
      if (cp08_payload_any_tag(x)) cp08_payload_tag_names(x, acc);
      return; /* no tagged node: getAttrib() yields NULL, which charges nothing */
    }
  }
  cp08_payload_walk(value, 1, acc);
}

static SEXP cp08_payload_count_attrib_cb(SEXP tag, SEXP value, void *data) {
  (void) tag;
  (void) value;
  uint64_t *nvalues = (uint64_t *) data;
  *nvalues += 1ULL;
  return NULL;
}

typedef struct {
  SEXP x;
  cp08_payload_acc_t *acc;
} cp08_payload_attr_visit_t;

static SEXP cp08_payload_visit_attrib_cb(SEXP tag, SEXP value, void *data) {
  cp08_payload_attr_visit_t *visit = (cp08_payload_attr_visit_t *) data;
  PROTECT(tag);
  PROTECT(value);
  cp08_payload_attr_value(visit->x, tag, value, visit->acc);
  UNPROTECT(2);
  return NULL;
}

static void cp08_payload_attributes(SEXP x, cp08_payload_acc_t *acc) {
  uint64_t nvalues = 0ULL;
  int names_from_tags = (TYPEOF(x) == LISTSXP) && cp08_payload_any_tag(x);
  if (names_from_tags) nvalues += 1ULL;
  cp08_payload_each_attrib(x, cp08_payload_count_attrib_cb, &nvalues);
  if (nvalues == 0ULL) return; /* attributes(x) is NULL */
  cp08_payload_tick(acc);
  cp08_payload_add(acc, nvalues + 8ULL, 0ULL); /* the attribute list node, scan_attrs = FALSE */
  if (names_from_tags) cp08_payload_tag_names(x, acc);
  cp08_payload_attr_visit_t visit = {x, acc};
  cp08_payload_each_attrib(x, cp08_payload_visit_attrib_cb, &visit);
}

static void cp08_payload_walk(SEXP x, int scan_attrs, cp08_payload_acc_t *acc) {
  R_CheckStack();
  cp08_payload_tick(acc);
  if (x == R_NilValue || TYPEOF(x) == ENVSXP) return;
  uint64_t len = cp08_payload_length(x, acc->call_env);
  cp08_payload_add(acc, len + 8ULL, 0ULL);
  switch (TYPEOF(x)) {
    case VECSXP: {
      R_xlen_t n = XLENGTH(x);
      for (R_xlen_t i = 0; i < n; ++i) cp08_payload_walk(VECTOR_ELT(x, i), 1, acc);
      break;
    }
    case LISTSXP: {
      for (SEXP p = x; p != R_NilValue; p = CDR(p)) cp08_payload_walk(CAR(p), 1, acc);
      break;
    }
    case STRSXP: {
      R_xlen_t n = XLENGTH(x);
      uint64_t entries = 0ULL, bytes = 0ULL;
      for (R_xlen_t i = 0; i < n; ++i) {
        if (i > 0 && (i % 4096) == 0) R_CheckUserInterrupt();
        SEXP s = STRING_ELT(x, i);
        if (s != NA_STRING) {
          uint64_t nb = (uint64_t) LENGTH(s);
          bytes += nb;
          entries += (nb + 1ULL) / 2ULL;
        }
      }
      cp08_payload_add(acc, entries, bytes);
      break;
    }
    default:
      break;
  }
  if (scan_attrs < 0 && TYPEOF(acc->scan_attrs) == LGLSXP &&
      XLENGTH(acc->scan_attrs) == 1 && LOGICAL(acc->scan_attrs)[0] != NA_LOGICAL) {
    if (LOGICAL(acc->scan_attrs)[0]) cp08_payload_attributes(x, acc);
  } else if (scan_attrs < 0) {
    /* Preserve the old R condition, including scalar coercion and NA with no
       attributes, AFTER this node and its children have been charged. */
    int has_attrs = cp08_payload_has_raw_attrib(x) ||
      (TYPEOF(x) == LISTSXP && cp08_payload_any_tag(x));
    SEXP quoted = PROTECT(Rf_lang2(Rf_install("quote"), acc->scan_attrs));
    SEXP yes = PROTECT(Rf_ScalarLogical(1));
    SEXP no = PROTECT(Rf_ScalarLogical(0));
    SEXP both = PROTECT(Rf_lang3(Rf_install("&&"), quoted, has_attrs ? yes : no));
    SEXP branch = PROTECT(Rf_lang4(Rf_install("if"), both, yes, no));
    SEXP value = PROTECT(Rf_eval(branch, acc->call_env));
    int do_attrs = LOGICAL(value)[0];
    UNPROTECT(6);
    if (do_attrs) cp08_payload_attributes(x, acc);
  } else if (scan_attrs) cp08_payload_attributes(x, acc);
}

SEXP C_gmm_cp08_payload(SEXP x_sexp, SEXP scan_attrs_sexp, SEXP call_env_sexp) {
  R_CheckUserInterrupt();
  gmm_require_valid_arithmetic(0);
  if (TYPEOF(call_env_sexp) != ENVSXP)
    gmm_c_abort(0, "Payload accounting requires its R call environment.", "cdid_gmm_bound_error");
  cp08_payload_acc_t acc = {0ULL, 0ULL, 0ULL, call_env_sexp, scan_attrs_sexp};
  cp08_payload_walk(x_sexp, -1, &acc);
  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, 2));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
  REAL(out_sexp)[0] = (double) acc.entries;
  REAL(out_sexp)[1] = (double) acc.char_bytes;
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("entries"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("char_bytes"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  UNPROTECT(2);
  return out_sexp;
}

static const R_CallMethodDef CallEntries[] = {
  {"C_gmm_native_status", (DL_FUNC) &C_gmm_native_status, 0},
  {"C_gmm_test_arithmetic_probe", (DL_FUNC) &C_gmm_test_arithmetic_probe, 5},
  {"C_gmm_sexp_aliases", (DL_FUNC) &C_gmm_sexp_aliases, 2},
  {"C_gmm_next", (DL_FUNC) &C_gmm_next, 2},
  {"C_gmm_iv_dot", (DL_FUNC) &C_gmm_iv_dot, 3},
  {"C_gmm_iv_norm", (DL_FUNC) &C_gmm_iv_norm, 3},
  {"C_gmm_iv_mat_residual", (DL_FUNC) &C_gmm_iv_mat_residual, 4},
  {"C_gmm_iv_eigen_residual", (DL_FUNC) &C_gmm_iv_eigen_residual, 4},
  {"C_gmm_iv_mat_product", (DL_FUNC) &C_gmm_iv_mat_product, 5},
  {"C_gmm_ipw_if_chunk", (DL_FUNC) &C_gmm_ipw_if_chunk, 18},
  {"C_gmm_cp08_support_chunk", (DL_FUNC) &C_gmm_cp08_support_chunk, 12},
  {"C_gmm_cp08_pair_rows", (DL_FUNC) &C_gmm_cp08_pair_rows, 9},
  {"C_gmm_cp08_sparse_pair_rows", (DL_FUNC) &C_gmm_cp08_sparse_pair_rows, 6},
  {"C_gmm_cp08_fit_pass1", (DL_FUNC) &C_gmm_cp08_fit_pass1, 15},
  {"C_gmm_cp08_fit_pass2", (DL_FUNC) &C_gmm_cp08_fit_pass2, 15},
  {"C_gmm_cp08_or_if_chunk", (DL_FUNC) &C_gmm_cp08_or_if_chunk, 17},
  {"C_gmm_cp08_dr_if_prepare", (DL_FUNC) &C_gmm_cp08_dr_if_prepare, 14},
  {"C_gmm_cp08_dr_if_assemble", (DL_FUNC) &C_gmm_cp08_dr_if_assemble, 17},
  {"C_gmm_cp08_dr_treatment_batch", (DL_FUNC) &C_gmm_cp08_dr_treatment_batch, 19},
  {"C_gmm_if_unit_batch", (DL_FUNC) &C_gmm_if_unit_batch, 10},
  {"C_gmm_cp08_payload", (DL_FUNC) &C_gmm_cp08_payload, 3},
  {NULL, NULL, 0}
};

void R_init_cdid(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
