#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("fp-contract=off")
#endif
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#pragma STDC FP_CONTRACT OFF
#endif

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
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
void R_init_cdid(DllInfo *dll);

SEXP C_gmm_native_status(void) {
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
  return Rf_ScalarLogical(a_sexp != R_NilValue && a_sexp == b_sexp);
}

SEXP C_gmm_next(SEXP x_sexp, SEXP up_sexp) {
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
  {NULL, NULL, 0}
};

void R_init_cdid(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
