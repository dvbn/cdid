/* Existing complete witnessed score records only. This worker does NOT read
   outcomes, fit/solve factors, assemble RHS/loadings or change finite guards. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)) || \
    (defined(__ASSOC_MATH__) && (__ASSOC_MATH__ > 0)) || \
    (defined(__RECIPROCAL_MATH__) && (__RECIPROCAL_MATH__ > 0))
#error "Unsafe floating-point compilation is forbidden for score pullbacks."
#endif
extern void cdid_pair_operation_arithmetic(void);

static int sp_plain(SEXP x, int type) {
  return TYPEOF(x) == type && !OBJECT(x) && !ALTREP(x) && ATTRIB(x) == R_NilValue;
}
static int sp_list(SEXP x) {
  if (TYPEOF(x) != VECSXP || OBJECT(x) || ALTREP(x)) return 0;
  SEXP a = ATTRIB(x), n = Rf_getAttrib(x, R_NamesSymbol);
  return a != R_NilValue && TAG(a) == R_NamesSymbol && CDR(a) == R_NilValue &&
    sp_plain(n, STRSXP) && XLENGTH(n) == XLENGTH(x);
}
static SEXP sp_field(SEXP x, const char *name) {
  SEXP n = Rf_getAttrib(x, R_NamesSymbol);
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
    if (STRING_ELT(n, i) != NA_STRING && !strcmp(CHAR(STRING_ELT(n, i)), name)) return VECTOR_ELT(x, i);
  return R_NilValue;
}
static int sp_scalar(SEXP x, int minimum, int *out) {
  if (!((sp_plain(x, INTSXP) || sp_plain(x, REALSXP)) && XLENGTH(x) == 1)) return 0;
  double d = TYPEOF(x) == INTSXP ? INTEGER(x)[0] : REAL(x)[0];
  if (!R_FINITE(d) || d < minimum || d > INT_MAX || d != floor(d)) return 0;
  *out = (int) d; return 1;
}
static int sp_matrix(SEXP x, int n, int p) {
  if (TYPEOF(x) != REALSXP || OBJECT(x) || ALTREP(x)) return 0;
  for (SEXP a = ATTRIB(x); a != R_NilValue; a = CDR(a))
    if (TAG(a) != R_DimSymbol && TAG(a) != R_DimNamesSymbol) return 0;
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  return sp_plain(d, INTSXP) && XLENGTH(d) == 2 && INTEGER(d)[0] == n && INTEGER(d)[1] == p &&
    XLENGTH(x) == (R_xlen_t) n * p;
}
/* Constructors explicitly name mean_d with the fitted kept-column names.
   These labels do not change sweep values; other decorations stay scalar. */
static int sp_mean(SEXP x, int p, SEXP columns) {
  if (TYPEOF(x) != REALSXP || OBJECT(x) || ALTREP(x) || XLENGTH(x) != p) return 0;
  if (ATTRIB(x) == R_NilValue) return 1;
  SEXP a = ATTRIB(x), names = Rf_getAttrib(x, R_NamesSymbol);
  return TAG(a) == R_NamesSymbol && CDR(a) == R_NilValue &&
    sp_plain(names, STRSXP) && sp_plain(columns, STRSXP) && XLENGTH(names) == p && XLENGTH(columns) == p &&
    R_compute_identical(names, columns, IDENT_USE_CLOENV);
}
static int sp_ids(SEXP x, R_xlen_t upper, int sorted) {
  if (!sp_plain(x, INTSXP)) return 0;
  int old = 0;
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
    int v = INTEGER(x)[i];
    if (v < 1 || v > upper || (sorted && v <= old)) return 0;
    old = v;
  }
  return 1;
}
static int sp_find(SEXP ids, int id) {
  R_xlen_t lo = 0, hi = XLENGTH(ids);
  while (lo < hi) {
    R_xlen_t mid = lo + (hi - lo) / 2; int key = INTEGER(ids)[mid];
    if (key == id) return (int) mid;
    if (key < id) lo = mid + 1; else hi = mid;
  }
  return -1;
}
static double sp_add(double a, double b) { volatile double x = a + b; return x; }
static double sp_sub(double a, double b) { volatile double x = a - b; return x; }
static double sp_mul(double a, double b) { volatile double x = a * b; return x; }

SEXP C_gmm_reconstruction_score_pullback(SEXP view, SEXP state, SEXP adj,
                                        SEXP ids, SEXP centered, SEXP cap) {
  if (!sp_list(view) || !sp_list(state) || !sp_list(adj) ||
      !((sp_plain(cap, REALSXP) || sp_plain(cap, INTSXP)) && XLENGTH(cap) == 1)) return R_NilValue;
  double limit = TYPEOF(cap) == REALSXP ? REAL(cap)[0] : INTEGER(cap)[0];
  if (!R_FINITE(limit) || limit <= 0) return R_NilValue;
  SEXP cells = sp_field(state, "cells"), groups = sp_field(state, "groups");
  SEXP position = sp_field(view, "position"), records = sp_field(view, "cache");
  SEXP intercept = sp_field(state, "intercept");
  int B = 0;
  if (!sp_plain(cells, VECSXP) || !sp_plain(groups, VECSXP) ||
      XLENGTH(cells) > INT_MAX || XLENGTH(groups) > INT_MAX ||
      !sp_plain(position, INTSXP) || XLENGTH(position) != XLENGTH(cells) ||
      !sp_plain(records, VECSXP) || !sp_scalar(sp_field(state, "B"), 1, &B) ||
      !sp_plain(intercept, LGLSXP) || XLENGTH(intercept) != 1 || LOGICAL(intercept)[0] == NA_LOGICAL ||
      !sp_ids(ids, XLENGTH(cells), 0) || XLENGTH(ids) > XLENGTH(cells) ||
      (centered != R_NilValue && (!sp_plain(centered, VECSXP) || XLENGTH(centered) != XLENGTH(records))))
    return R_NilValue;
  int has_intercept = LOGICAL(intercept)[0];
  SEXP kind = sp_field(adj, "kind"), lambda = sp_field(adj, "tilde_lambda"), nus = sp_field(adj, "tilde_nu");
  int compact = sp_plain(kind, STRSXP) && XLENGTH(kind) == 1 && STRING_ELT(kind, 0) != NA_STRING &&
    !strcmp(CHAR(STRING_ELT(kind, 0)), "compact_reconstruction_adjoint");
  SEXP cids = compact ? sp_field(adj, "cell_ids") : R_NilValue;
  SEXP gids = compact ? sp_field(adj, "group_ids") : R_NilValue;
  if (!sp_plain(lambda, REALSXP) || !sp_plain(nus, VECSXP)) return R_NilValue;
  if (compact) {
    if (!sp_ids(cids, XLENGTH(cells), 1) || !sp_ids(gids, XLENGTH(groups), 1) ||
        XLENGTH(lambda) != XLENGTH(cids) || XLENGTH(nus) != XLENGTH(gids)) return R_NilValue;
  } else if (XLENGTH(lambda) != XLENGTH(cells) || XLENGTH(nus) != XLENGTH(groups)) return R_NilValue;

  int peak_p = 0;
  /* Admission reads only plain metadata/operands; no numeric product, callback,
     new finiteness refusal or owned workspace precedes complete admission. */
  for (R_xlen_t j = 0; j < XLENGTH(ids); ++j) {
    int cid = INTEGER(ids)[j], slot = INTEGER(position)[cid - 1];
    if (slot < 1 || slot > XLENGTH(records)) return R_NilValue;
    SEXP rec = VECTOR_ELT(records, slot - 1), cf = VECTOR_ELT(cells, cid - 1);
    if (!sp_list(cf)) return R_NilValue;
    int gid = 0, p = 0;
    if (!sp_scalar(sp_field(cf, "group_id"), 1, &gid) || gid > XLENGTH(groups)) return R_NilValue;
    SEXP gf = VECTOR_ELT(groups, gid - 1);
    if (!sp_list(gf) || !sp_scalar(sp_field(gf, "p"), 0, &p)) return R_NilValue;
    if (p > peak_p) peak_p = p;
    if (rec == R_NilValue) continue;
    if (!sp_list(rec)) return R_NilValue;
    SEXP rr = sp_field(rec, "rr"), score = sp_field(rec, "s_rr"), d = sp_field(rec, "d_rr");
    if (!sp_ids(rr, B, 1) || XLENGTH(rr) > B || !sp_plain(score, REALSXP) ||
        XLENGTH(score) != XLENGTH(rr) || !sp_matrix(d, (int) XLENGTH(rr), p)) return R_NilValue;
    if (has_intercept && p && !sp_mean(sp_field(cf, "mean_d"), p, sp_field(gf, "kept_columns")))
      return R_NilValue;
    int gi = compact ? sp_find(gids, gid) : gid - 1;
    if (gi >= 0 && (!sp_plain(VECTOR_ELT(nus, gi), REALSXP) || XLENGTH(VECTOR_ELT(nus, gi)) != p))
      return R_NilValue;
    if (centered != R_NilValue) {
      SEXP cd = VECTOR_ELT(centered, slot - 1);
      /* Old R chooses this supplied entry, even when NULL. A missing p>0
         entry must return to that original error, never recompute a sweep. */
      if ((has_intercept && p && cd == R_NilValue) ||
          (cd != R_NilValue && !sp_matrix(cd, (int) XLENGTH(rr), p))) return R_NilValue;
    }
  }
  /* One output B, current/replacement centered matrices <=2Bp, base matprod
     result/temporary and R call/index/header overlap. No R_alloc in the loop.
     This entire peak fits the REPLACED scalar-adjoint compartment, separate
     from V49 position state, score/centered records and V50 forcing records. */
  double required = 1024.0 + (double) B * (4.0 * peak_p + 8.0) + 8.0 * peak_p;
  double funded = 4096.0 + 8.0 * B * (peak_p + 4.0) + 8.0 * peak_p * peak_p;
  if (!R_FINITE(required) || !R_FINITE(funded) || required > funded || required > limit ||
      (double) B * peak_p > R_XLEN_T_MAX) return R_NilValue;
  cdid_pair_operation_arithmetic();
  SEXP out = PROTECT(Rf_allocVector(REALSXP, B));
  for (int i = 0; i < B; ++i) REAL(out)[i] = 0.0;
  SEXP dot_call = PROTECT(Rf_lang3(Rf_install("%*%"), R_NilValue, R_NilValue));
  for (R_xlen_t j = 0; j < XLENGTH(ids); ++j) {
    if (j && !(j % 256)) R_CheckUserInterrupt();
    int cid = INTEGER(ids)[j], slot = INTEGER(position)[cid - 1];
    SEXP rec = VECTOR_ELT(records, slot - 1);
    if (rec == R_NilValue) continue;
    SEXP cf = VECTOR_ELT(cells, cid - 1); int gid = 0, p = 0;
    sp_scalar(sp_field(cf, "group_id"), 1, &gid);
    sp_scalar(sp_field(VECTOR_ELT(groups, gid - 1), "p"), 0, &p);
    SEXP rr = sp_field(rec, "rr"), score = sp_field(rec, "s_rr"), d = sp_field(rec, "d_rr");
    int n = (int) XLENGTH(rr), ci = compact ? sp_find(cids, cid) : cid - 1;
    double lam = has_intercept && ci >= 0 ? REAL(lambda)[ci] : 0.0;
    int gi = compact ? sp_find(gids, gid) : gid - 1, np = 0;
    SEXP nu = gi >= 0 ? VECTOR_ELT(nus, gi) : R_NilValue;
    if (nu == R_NilValue) {
      nu = PROTECT(Rf_allocVector(REALSXP, p)); ++np;
      for (int k = 0; k < p; ++k) REAL(nu)[k] = 0.0;
    }
    SEXP product = R_NilValue;
    if (p) {
      SEXP centered_d = has_intercept && centered != R_NilValue ? VECTOR_ELT(centered, slot - 1) : R_NilValue;
      if (centered_d == R_NilValue) {
        if (has_intercept) {
          centered_d = PROTECT(Rf_allocMatrix(REALSXP, n, p)); ++np;
          SEXP mean = sp_field(cf, "mean_d");
          for (int k = 0; k < p; ++k) for (int i = 0; i < n; ++i)
            REAL(centered_d)[i + (R_xlen_t) n * k] = sp_sub(REAL(d)[i + (R_xlen_t) n * k], REAL(mean)[k]);
        } else centered_d = d;
      }
      SETCADR(dot_call, centered_d); SETCADDR(dot_call, nu);
      product = PROTECT(Rf_eval(dot_call, R_BaseEnv)); ++np;
      /* Same base one-vector product and as.vector numeric values as R. */
      if (TYPEOF(product) != REALSXP || XLENGTH(product) != n)
        Rf_error("Unexpected base matrix-product result in reconstruction score pullback.");
    }
    for (int i = 0; i < n; ++i) {
      double val = p ? sp_add(lam, REAL(product)[i]) : lam;
      double contribution = sp_mul(REAL(score)[i], val);
      int row = INTEGER(rr)[i] - 1;
      REAL(out)[row] = sp_add(REAL(out)[row], contribution);
    }
    SETCADR(dot_call, R_NilValue); SETCADDR(dot_call, R_NilValue);
    UNPROTECT(np);
  }
  UNPROTECT(2); return out;
}
