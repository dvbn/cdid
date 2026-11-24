#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("fp-contract=off")
#endif
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#pragma STDC FP_CONTRACT OFF
#endif

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <R.h>
#include <Rinternals.h>
#include "cdid_r_api.h"
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
/* Older R BLAS and LAPACK headers use C int for dimensions. */
#ifndef BLAS_INT
#define BLAS_INT int
#endif
#include <R_ext/Utils.h>
#ifndef FCONE
#define FCONE
#endif
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#define CDID_SEL_CLIP_EPS (1e-6)
#define CDID_SEL_CLIP_CAP (1.0 - 1e-6)
#define CDID_SEL_CLIP_CAP_ODDS ((1.0 - 1e-6) / (1.0 - (1.0 - 1e-6)))

static inline SEXP cdid_sel_elt(SEXP list, const char *nm) {
  if (TYPEOF(list) != VECSXP) return R_NilValue;
  SEXP nms = Rf_getAttrib(list, R_NamesSymbol);
  if (TYPEOF(nms) != STRSXP) return R_NilValue;
  R_xlen_t n = XLENGTH(list);
  if (XLENGTH(nms) < n) n = XLENGTH(nms);
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP s = STRING_ELT(nms, i);
    if (s != NA_STRING && strcmp(CHAR(s), nm) == 0) return VECTOR_ELT(list, i);
  }
  return R_NilValue;
}

static inline double cdid_sel_num_at(SEXP v, R_xlen_t idx) {
  if (TYPEOF(v) == REALSXP) return REAL(v)[idx];
  if (TYPEOF(v) == INTSXP) {
    int iv = INTEGER(v)[idx];
    return (iv == NA_INTEGER) ? R_NaReal : (double) iv;
  }
  return R_NaReal;
}

static inline int cdid_sel_int_at(SEXP v, R_xlen_t idx) {
  if (TYPEOF(v) == INTSXP) return INTEGER(v)[idx];
  if (TYPEOF(v) == REALSXP) {
    double dv = REAL(v)[idx];
    if (!R_FINITE(dv) || dv != floor(dv) || dv < (double) INT_MIN || dv > (double) INT_MAX) return NA_INTEGER;
    return (int) dv;
  }
  return NA_INTEGER;
}

static inline void cdid_sel_clip_logit_scalar(double eta, double *prob_out, double *odds_out, double *act_out) {
  double odds, prob, active;
  if (R_FINITE(eta) && eta > 40.0) {
    odds = 1.0;
    prob = 1.0 / (1.0 + exp(-eta));
  } else {
    odds = exp(eta);
    prob = odds / (1.0 + odds);
  }
  if (R_FINITE(prob) && prob < CDID_SEL_CLIP_CAP) {
    active = 1.0;
  } else if (ISNAN(prob)) {
    active = R_NaReal;
  } else {
    active = 0.0;
    odds = CDID_SEL_CLIP_CAP_ODDS;
  }
  *prob_out = prob;
  *odds_out = odds;
  *act_out = active;
}

/* 0-based exact binary search on strictly increasing double table[0..n-1] with linear-next hint */
static inline int cdid_sel_bsearch_dbl(double target, const double *tab, int n, int *hint) {
  if (n <= 0 || ISNAN(target)) return -1;
  int h = *hint;
  if (h >= 0 && h + 1 < n && tab[h + 1] == target) {
    *hint = h + 1;
    return h + 1;
  }
  if (h >= 0 && h < n && tab[h] == target) {
    return h;
  }
  int lo = 0, hi = n - 1;
  while (lo <= hi) {
    int mid = lo + ((hi - lo) >> 1);
    double val = tab[mid];
    if (val == target) {
      *hint = mid;
      return mid;
    }
    if (val < target) lo = mid + 1;
    else hi = mid - 1;
  }
  return -1;
}

/* 0-based exact binary search on sorted int sub-array tab[lo0..hi0] (inclusive, 0-based) */
static inline int cdid_sel_bsearch_int_range(int target, const int *tab, int lo0, int hi0, int *hint) {
  if (lo0 > hi0) return -1;
  int h = *hint;
  if (h >= lo0 && h + 1 <= hi0 && tab[h + 1] == target) {
    *hint = h + 1;
    return h + 1;
  }
  int lo = lo0, hi = hi0;
  while (lo <= hi) {
    int mid = lo + ((hi - lo) >> 1);
    int val = tab[mid];
    if (val == target) {
      *hint = mid;
      return mid;
    }
    if (val < target) lo = mid + 1;
    else hi = mid - 1;
  }
  return -1;
}

static inline int cdid_sel_all_finite(const double *x, R_xlen_t n) {
  if (n > 0 && x == NULL) return 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!R_FINITE(x[i])) return 0;
  }
  return 1;
}

static inline void cdid_dd_add(double *hi, double *lo, double x) {
  double s = *hi + x;
  double v = s - *hi;
  double e = (*hi - (s - v)) + (x - v);
  *hi = s;
  *lo += e;
}

static inline void cdid_dd_prod_add(double *hi, double *lo, double a, double b) {
  double p = a * b;
  double err = fma(a, b, -p);
  cdid_dd_add(hi, lo, p);
  *lo += err;
}

/* Compensated Householder QR WLS solver for row-permutation-invariant block selection fits */
static int cdid_sel_qr_wls_solve_dd(
    double *W_X,
    double *w_z,
    double *R_out,
    double *beta_out,
    int n,
    int p,
    double rel_tol) {
  if (n < p || p < 1) return 0;
  for (int step = 0; step < p; ++step) {
    double *col_s = W_X + (size_t) step * (size_t) n;
    double alpha = col_s[step];
    double xnorm = 0.0;
    if (step + 1 < n) {
      double scale = 0.0;
      for (int i = step + 1; i < n; ++i) {
        double ax = fabs(col_s[i]);
        if (ax > scale) scale = ax;
      }
      if (scale > 0.0) {
        double inv_sc = 1.0 / scale;
        double s_hi = 0.0, s_lo = 0.0;
        for (int i = step + 1; i < n; ++i) {
          double u = col_s[i] * inv_sc;
          cdid_dd_prod_add(&s_hi, &s_lo, u, u);
        }
        xnorm = scale * sqrt(s_hi + s_lo);
      }
    }
    double norm_full = hypot(alpha, xnorm);
    if (!R_FINITE(norm_full) || !(norm_full > 0.0)) return 0;
    if (xnorm > 0.0) {
      double beta = (alpha >= 0.0) ? -norm_full : norm_full;
      double tau = (beta - alpha) / beta;
      double inv_denom = 1.0 / (alpha - beta);
      for (int i = step + 1; i < n; ++i) col_s[i] *= inv_denom;
      col_s[step] = beta;

      for (int j = step + 1; j < p; ++j) {
        double *col_j = W_X + (size_t) j * (size_t) n;
        double d_hi = col_j[step], d_lo = 0.0;
        for (int i = step + 1; i < n; ++i) cdid_dd_prod_add(&d_hi, &d_lo, col_s[i], col_j[i]);
        double tdot = tau * (d_hi + d_lo);
        col_j[step] -= tdot;
        for (int i = step + 1; i < n; ++i) col_j[i] = fma(-tdot, col_s[i], col_j[i]);
      }

      if (w_z != NULL) {
        double dz_hi = w_z[step], dz_lo = 0.0;
        for (int i = step + 1; i < n; ++i) cdid_dd_prod_add(&dz_hi, &dz_lo, col_s[i], w_z[i]);
        double tdot_z = tau * (dz_hi + dz_lo);
        w_z[step] -= tdot_z;
        for (int i = step + 1; i < n; ++i) w_z[i] = fma(-tdot_z, col_s[i], w_z[i]);
      }
    }
  }

  memset(R_out, 0, (size_t) p * (size_t) p * sizeof(double));
  double max_d = 0.0, min_d = DBL_MAX;
  for (int j = 0; j < p; ++j) {
    const double *col_j = W_X + (size_t) j * (size_t) n;
    for (int i = 0; i <= j; ++i) {
      if (!R_FINITE(col_j[i])) return 0;
      R_out[i + j * p] = col_j[i];
    }
    double ad = fabs(col_j[j]);
    if (ad > max_d) max_d = ad;
    if (ad < min_d) min_d = ad;
  }
  if (!(max_d > 0.0) || !(min_d > rel_tol * max_d)) return 0;

  if (w_z != NULL && beta_out != NULL) {
    for (int i = p - 1; i >= 0; --i) {
      double s_hi = w_z[i], s_lo = 0.0;
      for (int k = i + 1; k < p; ++k) cdid_dd_prod_add(&s_hi, &s_lo, -R_out[i + k * p], beta_out[k]);
      double b_i = (s_hi + s_lo) / R_out[i + i * p];
      if (!R_FINITE(b_i)) return 0;
      beta_out[i] = b_i;
    }
  }
  return 1;
}

/* Unpivoted QR WLS solver via R's Cdqrls matching stats::glm.fit to 0 ULP */
static int cdid_sel_qr_wls_solve(
    double *W_X,
    double *w_z,
    double *R_out,
    double *beta_out,
    int n,
    int p,
    double dqr_tol) {
  if (n < p || p < 1 || p > 256 || w_z == NULL || beta_out == NULL) return 0;
  int ny = 1;
  int rank = 0;
  int jpvt[256];
  double qraux[256];
  double work[512];
  double b_piv[256];
  const void *vmax = vmaxget();
  double *rsd = (double *) R_alloc((size_t) n, sizeof(double));
  double *qty = (double *) R_alloc((size_t) n, sizeof(double));
  for (int j = 0; j < p; ++j) {
    jpvt[j] = j + 1;
    qraux[j] = 0.0;
    b_piv[j] = 0.0;
    work[2 * j] = 0.0;
    work[2 * j + 1] = 0.0;
  }
  F77_CALL(dqrls)(W_X, &n, &p, w_z, &ny, &dqr_tol, b_piv, rsd, qty, &rank, jpvt, qraux, work);
  vmaxset(vmax);
  if (rank < p) return 0;
  for (int j = 0; j < p; ++j) {
    if (jpvt[j] != j + 1 || !R_FINITE(b_piv[j])) return 0;
    beta_out[j] = b_piv[j];
  }

  memset(R_out, 0, (size_t) p * (size_t) p * sizeof(double));
  double max_d = 0.0;
  for (int j = 0; j < p; ++j) {
    const double *col_j = W_X + (size_t) j * (size_t) n;
    for (int i = 0; i <= j; ++i) {
      if (!R_FINITE(col_j[i])) return 0;
      R_out[i + j * p] = col_j[i];
    }
    double ad = fabs(col_j[j]);
    if (ad > max_d) max_d = ad;
  }
  if (!(max_d > 0.0)) return 0;
  return 1;
}

static int cdid_sel_chol_upper(const double *H, double *R_out, int p, double rel_tol) {
  if (p < 1) return 0;
  double max_diag = 0.0;
  for (int j = 0; j < p; ++j) {
    double hd = H[j + j * p];
    if (!R_FINITE(hd)) return 0;
    if (hd > max_diag) max_diag = hd;
  }
  if (!(max_diag > 0.0)) return 0;
  for (int idx = 0; idx < p * p; ++idx) {
    if (!R_FINITE(H[idx])) return 0;
    R_out[idx] = H[idx];
  }
  int p_b = p, info = 0;
  F77_CALL(dpotrf)("U", &p_b, R_out, &p_b, &info FCONE);
  if (info != 0) return 0;
  double min_allowed = rel_tol * sqrt(max_diag);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i <= j; ++i) {
      if (!R_FINITE(R_out[i + j * p])) return 0;
    }
    for (int i = j + 1; i < p; ++i) {
      R_out[i + j * p] = 0.0;
    }
    double rjj = R_out[j + j * p];
    if (!(rjj > 0.0) || !(rjj > min_allowed)) return 0;
  }
  return 1;
}

/* Solve (R^T R) x = rhs * scale_norm for B rows (rhs is B x p col-major, out is B x p col-major).
   If pivot1 != NULL (1-based length p), solves in pivot order: P^T R^T R P x = rhs * scale_norm. */
static int cdid_sel_rtr_solve_rows(
    const double *R_up,
    const int *pivot1,
    const double *rhs,
    double *out,
    int B,
    int p,
    double scale_norm,
    double *tmp_p) {
  (void) tmp_p;
  if (p == 0 || B == 0) return 1;
  if (!R_FINITE(scale_norm)) return 0;
  for (int j = 0; j < p; ++j) {
    int pj = pivot1 ? (pivot1[j] - 1) : j;
    if (pj < 0 || pj >= p) return 0;
    double d = R_up[j + j * p];
    if (!R_FINITE(d) || d == 0.0) return 0;
    for (int i = 0; i <= j; ++i) {
      if (!R_FINITE(R_up[i + j * p])) return 0;
    }
  }
  const void *vm = vmaxget();
  double *L_low = (double *) R_alloc((size_t) p * (size_t) p, sizeof(double));
  double *work = (double *) R_alloc((size_t) p * (size_t) B, sizeof(double));
  for (int c = 0; c < p; ++c) {
    for (int r = 0; r < p; ++r) {
      L_low[r + (size_t) c * (size_t) p] = (r >= c) ? R_up[c + (size_t) r * (size_t) p] : 0.0;
    }
  }
  for (int j = 0; j < p; ++j) {
    int pj = pivot1 ? (pivot1[j] - 1) : j;
    const double *rhs_col = rhs + (size_t) pj * (size_t) B;
    for (int r = 0; r < B; ++r) {
      double v = rhs_col[r];
      if (!R_FINITE(v)) { vmaxset(vm); return 0; }
      work[j + (size_t) r * (size_t) p] = v;
    }
  }
  BLAS_INT pb = (BLAS_INT) p;
  BLAS_INT Bb = (BLAS_INT) B;
  double one = 1.0;
  F77_CALL(dtrsm)("L", "L", "N", "N", &pb, &Bb, &one, L_low, &pb, work, &pb FCONE FCONE FCONE FCONE);
  F77_CALL(dtrsm)("L", "U", "N", "N", &pb, &Bb, &one, R_up, &pb, work, &pb FCONE FCONE FCONE FCONE);
  for (int j = 0; j < p; ++j) {
    int pj = pivot1 ? (pivot1[j] - 1) : j;
    double *out_col = out + (size_t) pj * (size_t) B;
    for (int r = 0; r < B; ++r) {
      double val = work[j + (size_t) r * (size_t) p] * scale_norm;
      if (!R_FINITE(val)) { vmaxset(vm); return 0; }
      out_col[r] = val;
    }
  }
  vmaxset(vm);
  return 1;
}

/* ---------------------------------------------------------------------------
   1. C_cdid_sel_ref_lookup(keys_sexp, table_sexp)
   --------------------------------------------------------------------------- */
SEXP C_cdid_sel_ref_lookup(SEXP keys_sexp, SEXP table_sexp) {
  if (Rf_isObject(keys_sexp) || Rf_isObject(table_sexp)) return R_NilValue;
  int tk = TYPEOF(keys_sexp);
  int tt = TYPEOF(table_sexp);
  if ((tk != REALSXP && tk != INTSXP) || (tt != REALSXP && tt != INTSXP)) return R_NilValue;
  R_xlen_t m = XLENGTH(keys_sexp);
  R_xlen_t n_xl = XLENGTH(table_sexp);
  if (n_xl > (R_xlen_t) INT_MAX) return R_NilValue;
  int n = (int) n_xl;

  const void *vmax = vmaxget();
  const double *tab = NULL;
  if (tt == REALSXP) {
    tab = REAL(table_sexp);
    for (int i = 0; i < n; ++i) {
      if (ISNAN(tab[i]) || (i > 0 && !(tab[i] > tab[i - 1]))) {
        vmaxset(vmax);
        return R_NilValue;
      }
    }
  } else {
    const int *ti = INTEGER(table_sexp);
    double *td = (n > 0) ? (double *) R_alloc((size_t) n, sizeof(double)) : NULL;
    for (int i = 0; i < n; ++i) {
      if (ti[i] == NA_INTEGER || (i > 0 && !(ti[i] > ti[i - 1]))) {
        vmaxset(vmax);
        return R_NilValue;
      }
      td[i] = (double) ti[i];
    }
    tab = td;
  }

  if (tk == REALSXP) {
    const double *kv = REAL(keys_sexp);
    for (R_xlen_t i = 0; i < m; ++i) {
      if (ISNAN(kv[i])) { vmaxset(vmax); return R_NilValue; }
    }
    SEXP out_s = PROTECT(Rf_allocVector(INTSXP, m));
    int *out = INTEGER(out_s);
    int hint = -1;
    for (R_xlen_t i = 0; i < m; ++i) {
      int pos = cdid_sel_bsearch_dbl(kv[i], tab, n, &hint);
      out[i] = (pos >= 0) ? (pos + 1) : NA_INTEGER;
    }
    UNPROTECT(1);
    vmaxset(vmax);
    return out_s;
  } else {
    const int *kv = INTEGER(keys_sexp);
    for (R_xlen_t i = 0; i < m; ++i) {
      if (kv[i] == NA_INTEGER) { vmaxset(vmax); return R_NilValue; }
    }
    SEXP out_s = PROTECT(Rf_allocVector(INTSXP, m));
    int *out = INTEGER(out_s);
    int hint = -1;
    for (R_xlen_t i = 0; i < m; ++i) {
      int pos = cdid_sel_bsearch_dbl((double) kv[i], tab, n, &hint);
      out[i] = (pos >= 0) ? (pos + 1) : NA_INTEGER;
    }
    UNPROTECT(1);
    vmaxset(vmax);
    return out_s;
  }
}

/* ---------------------------------------------------------------------------
   2. C_cdid_sel_ref_slice(mode_sexp, rows_sexp, base_sexp, table_sexp, data_sexp, cols_sexp)
      mode 1: data_sexp is REALSXP vector -> returns REALSXP vector[m]
      mode 2: data_sexp is REALSXP matrix -> returns m x ncol matrix (with colnames)
      mode 3: data_sexp is REALSXP matrix, cols_sexp is 1-based INTSXP -> returns m x k matrix
      base_sexp: c(base_offset, N_max) as REALSXP length 2
   --------------------------------------------------------------------------- */
SEXP C_cdid_sel_ref_slice(
    SEXP mode_sexp,
    SEXP rows_sexp,
    SEXP base_sexp,
    SEXP table_sexp,
    SEXP data_sexp,
    SEXP cols_sexp) {
  int mode = Rf_asInteger(mode_sexp);
  if ((mode < 1 || mode > 3) || Rf_isObject(rows_sexp) || TYPEOF(table_sexp) != REALSXP ||
      TYPEOF(data_sexp) != REALSXP || TYPEOF(base_sexp) != REALSXP || XLENGTH(base_sexp) != 2) {
    return R_NilValue;
  }
  int tr = TYPEOF(rows_sexp);
  if (tr != INTSXP && tr != REALSXP) return R_NilValue;
  R_xlen_t m_xl = XLENGTH(rows_sexp);
  R_xlen_t n_xl = XLENGTH(table_sexp);
  if (m_xl > (R_xlen_t) INT_MAX || n_xl > (R_xlen_t) INT_MAX) return R_NilValue;
  int m = (int) m_xl;
  int n_tab = (int) n_xl;

  double base_off = REAL(base_sexp)[0];
  double N_max = REAL(base_sexp)[1];
  if (!R_FINITE(base_off) || !R_FINITE(N_max)) return R_NilValue;

  const double *tab = REAL(table_sexp);
  const void *vmax = vmaxget();
  int *hits = (m > 0) ? (int *) R_alloc((size_t) m, sizeof(int)) : NULL;
  int hint = -1;

  if (tr == INTSXP) {
    const int *rv = INTEGER(rows_sexp);
    for (int i = 0; i < m; ++i) {
      int r = rv[i];
      if (r == NA_INTEGER || r < 1 || (double) r > N_max) { vmaxset(vmax); return R_NilValue; }
      int pos = cdid_sel_bsearch_dbl(base_off + (double) r, tab, n_tab, &hint);
      if (pos < 0) { vmaxset(vmax); return R_NilValue; }
      hits[i] = pos;
    }
  } else {
    const double *rv = REAL(rows_sexp);
    for (int i = 0; i < m; ++i) {
      double r = rv[i];
      if (!R_FINITE(r) || r != floor(r) || r < 1.0 || r > N_max) { vmaxset(vmax); return R_NilValue; }
      int pos = cdid_sel_bsearch_dbl(base_off + r, tab, n_tab, &hint);
      if (pos < 0) { vmaxset(vmax); return R_NilValue; }
      hits[i] = pos;
    }
  }

  if (mode == 1) {
    if (XLENGTH(data_sexp) != (R_xlen_t) n_tab) { vmaxset(vmax); return R_NilValue; }
    const double *w = REAL(data_sexp);
    SEXP out_s = PROTECT(Rf_allocVector(REALSXP, m));
    double *out = REAL(out_s);
    for (int i = 0; i < m; ++i) out[i] = w[hits[i]];
    UNPROTECT(1);
    vmaxset(vmax);
    return out_s;
  }

  SEXP z_dim = Rf_getAttrib(data_sexp, R_DimSymbol);
  if (TYPEOF(z_dim) != INTSXP || XLENGTH(z_dim) != 2 || INTEGER(z_dim)[0] != n_tab) {
    vmaxset(vmax);
    return R_NilValue;
  }
  int p_full = INTEGER(z_dim)[1];
  const double *Z = REAL(data_sexp);

  if (mode == 2) {
    SEXP out_s = PROTECT(Rf_allocMatrix(REALSXP, m, p_full));
    double *out = REAL(out_s);
    for (int j = 0; j < p_full; ++j) {
      const double *zj = Z + (size_t) j * (size_t) n_tab;
      double *oj = out + (size_t) j * (size_t) m;
      for (int i = 0; i < m; ++i) oj[i] = zj[hits[i]];
    }
    SEXP zdn = Rf_getAttrib(data_sexp, R_DimNamesSymbol);
    if (TYPEOF(zdn) == VECSXP && XLENGTH(zdn) == 2 && VECTOR_ELT(zdn, 1) != R_NilValue) {
      SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(dn, 0, R_NilValue);
      SET_VECTOR_ELT(dn, 1, VECTOR_ELT(zdn, 1));
      Rf_setAttrib(out_s, R_DimNamesSymbol, dn);
      UNPROTECT(1);
    }
    UNPROTECT(1);
    vmaxset(vmax);
    return out_s;
  }

  /* mode == 3 */
  if (TYPEOF(cols_sexp) != INTSXP) { vmaxset(vmax); return R_NilValue; }
  int k_cols = (int) XLENGTH(cols_sexp);
  const int *cols = INTEGER(cols_sexp);
  for (int j = 0; j < k_cols; ++j) {
    if (cols[j] == NA_INTEGER || cols[j] < 1 || cols[j] > p_full) { vmaxset(vmax); return R_NilValue; }
  }
  SEXP out_s = PROTECT(Rf_allocMatrix(REALSXP, m, k_cols));
  double *out = REAL(out_s);
  for (int j = 0; j < k_cols; ++j) {
    const double *zj = Z + (size_t) (cols[j] - 1) * (size_t) n_tab;
    double *oj = out + (size_t) j * (size_t) m;
    for (int i = 0; i < m; ++i) oj[i] = zj[hits[i]];
  }
  SEXP zdn = Rf_getAttrib(data_sexp, R_DimNamesSymbol);
  if (TYPEOF(zdn) == VECSXP && XLENGTH(zdn) == 2) {
    SEXP cn = VECTOR_ELT(zdn, 1);
    if (TYPEOF(cn) == STRSXP && XLENGTH(cn) == p_full) {
      SEXP sub_cn = PROTECT(Rf_allocVector(STRSXP, k_cols));
      for (int j = 0; j < k_cols; ++j) SET_STRING_ELT(sub_cn, j, STRING_ELT(cn, cols[j] - 1));
      SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(dn, 0, R_NilValue);
      SET_VECTOR_ELT(dn, 1, sub_cn);
      Rf_setAttrib(out_s, R_DimNamesSymbol, dn);
      UNPROTECT(2);
    }
  }
  UNPROTECT(1);
  vmaxset(vmax);
  return out_s;
}

/* ---------------------------------------------------------------------------
   Internal panel / reference / date / pair helpers
   --------------------------------------------------------------------------- */
typedef struct {
  int N;
  int n_A;
  const int *a_index;      /* length N */
  const int *cohort_index; /* length N */
  const int *cohort_order; /* length N */
  SEXP cohort_start_s;     /* length J (INTSXP or REALSXP) */
  SEXP cohort_end_s;       /* length J (INTSXP or REALSXP) */
  int J;
  SEXP id_sexp;            /* INTSXP or REALSXP or STRSXP */
  const double *Z_static;  /* N x d_Z or NULL */
  int d_Z;
  const double *w_static;  /* length N */
  const int *U_logical;    /* length N (LGLSXP) */
  const double *rho;       /* length N */
  const double *a_vec;     /* length J */
  const double *mass_U;    /* length J */
  double N_dbl;

  /* Panel postings */
  int n_dates;
  const int *date_start;
  const int *date_end;
  const int *date_count;
  const int *post_a;
  const int *post_obs;
  R_xlen_t n_post;
  const double *y_obs;
  R_xlen_t n_y;

  /* Singleton lookup (for independent_period_logit) */
  int n_singletons;
  int sing_numeric_ok;
  const double *unit_id_key;
  int n_panel_units;
  const int *units_n_obs;
  SEXP units_obs_start_s;
  const int *obs_cal_idx;
} cdid_sel_reader_t;

static int cdid_sel_init_reader(SEXP reader_sexp, cdid_sel_reader_t *rd, int need_singletons) {
  memset(rd, 0, sizeof(*rd));
  SEXP st = cdid_sel_elt(reader_sexp, "state");
  SEXP panel_ref = cdid_sel_elt(reader_sexp, "panel_ref");
  if (TYPEOF(st) != VECSXP || TYPEOF(panel_ref) != ENVSXP) return 0;
  int panel_bound = 0;
  SEXP panel = cdid_frame_value(panel_ref, Rf_install("panel_owner"), &panel_bound);
  if (!panel_bound || TYPEOF(panel) != VECSXP) return 0;

  SEXP a_idx_s = cdid_sel_elt(st, "a_index");
  SEXP ci_s = cdid_sel_elt(st, "cohort_index");
  SEXP co_s = cdid_sel_elt(st, "cohort_order");
  SEXP cs_s = cdid_sel_elt(st, "cohort_start");
  SEXP ce_s = cdid_sel_elt(st, "cohort_end");
  SEXP U_s = cdid_sel_elt(st, "U");
  SEXP rho_s = cdid_sel_elt(st, "rho");
  SEXP w_s = cdid_sel_elt(st, "w");
  SEXP a_s = cdid_sel_elt(st, "a");
  SEXP mu_s = cdid_sel_elt(st, "mass_U");
  SEXP Z_s = cdid_sel_elt(st, "Z");
  SEXP N_s = cdid_sel_elt(st, "N");
  SEXP nA_s = cdid_sel_elt(st, "n_A");
  if (TYPEOF(a_idx_s) != INTSXP || TYPEOF(ci_s) != INTSXP || TYPEOF(co_s) != INTSXP ||
      (TYPEOF(cs_s) != INTSXP && TYPEOF(cs_s) != REALSXP) ||
      (TYPEOF(ce_s) != INTSXP && TYPEOF(ce_s) != REALSXP) || TYPEOF(U_s) != LGLSXP ||
      TYPEOF(rho_s) != REALSXP || TYPEOF(w_s) != REALSXP || TYPEOF(a_s) != REALSXP ||
      TYPEOF(mu_s) != REALSXP) return 0;

  int N = (int) XLENGTH(a_idx_s);
  int J = (int) XLENGTH(cs_s);
  if (N < 1 || J < 1 || XLENGTH(ci_s) != N || XLENGTH(co_s) != N || XLENGTH(U_s) != N ||
      XLENGTH(rho_s) != N || XLENGTH(w_s) != N || XLENGTH(ce_s) != J ||
      XLENGTH(a_s) != J || XLENGTH(mu_s) != J) return 0;

  rd->N = N;
  rd->n_A = Rf_asInteger(nA_s);
  rd->N_dbl = Rf_asReal(N_s);
  if (!R_FINITE(rd->N_dbl) || rd->N_dbl <= 0.0) return 0;
  rd->a_index = INTEGER(a_idx_s);
  rd->cohort_index = INTEGER(ci_s);
  rd->cohort_order = INTEGER(co_s);
  rd->cohort_start_s = cs_s;
  rd->cohort_end_s = ce_s;
  rd->J = J;
  rd->id_sexp = cdid_sel_elt(st, "id");
  rd->w_static = REAL(w_s);
  rd->U_logical = LOGICAL(U_s);
  rd->rho = REAL(rho_s);
  rd->a_vec = REAL(a_s);
  rd->mass_U = REAL(mu_s);

  if (TYPEOF(Z_s) == REALSXP) {
    SEXP zdim = Rf_getAttrib(Z_s, R_DimSymbol);
    if (TYPEOF(zdim) == INTSXP && XLENGTH(zdim) == 2 && INTEGER(zdim)[0] == N) {
      rd->Z_static = REAL(Z_s);
      rd->d_Z = INTEGER(zdim)[1];
    }
  }

  SEXP post = cdid_sel_elt(panel, "postings");
  SEXP obs = cdid_sel_elt(panel, "observations");
  if (TYPEOF(post) != VECSXP || TYPEOF(obs) != VECSXP) return 0;
  SEXP ds_s = cdid_sel_elt(post, "date_start");
  SEXP de_s = cdid_sel_elt(post, "date_end");
  SEXP dc_s = cdid_sel_elt(post, "date_count");
  SEXP pa_s = cdid_sel_elt(post, "post_a_index");
  SEXP po_s = cdid_sel_elt(post, "post_obs_index");
  SEXP y_s = cdid_sel_elt(obs, "y");
  if (TYPEOF(ds_s) != INTSXP || TYPEOF(de_s) != INTSXP || TYPEOF(dc_s) != INTSXP ||
      TYPEOF(pa_s) != INTSXP || Rf_isObject(pa_s) || TYPEOF(po_s) != INTSXP || Rf_isObject(po_s) ||
      TYPEOF(y_s) != REALSXP || Rf_isObject(y_s)) return 0;

  rd->n_dates = (int) XLENGTH(dc_s);
  rd->date_start = INTEGER(ds_s);
  rd->date_end = INTEGER(de_s);
  rd->date_count = INTEGER(dc_s);
  rd->post_a = INTEGER(pa_s);
  rd->post_obs = INTEGER(po_s);
  rd->n_post = XLENGTH(pa_s);
  rd->y_obs = REAL(y_s);
  rd->n_y = XLENGTH(y_s);
  if (XLENGTH(po_s) != rd->n_post) return 0;

  SEXP n_sing_s = cdid_sel_elt(panel, "n_singleton_units");
  rd->n_singletons = (n_sing_s != R_NilValue) ? Rf_asInteger(n_sing_s) : 0;
  rd->sing_numeric_ok = 1;
  if (rd->n_singletons > 0 && need_singletons) {
    SEXP id_type_s = cdid_sel_elt(panel, "id_type");
    SEXP ukey_s = cdid_sel_elt(panel, "unit_id_key");
    SEXP units_s = cdid_sel_elt(panel, "units");
    SEXP cal_s = cdid_sel_elt(obs, "calendar_index");
    if (TYPEOF(id_type_s) != STRSXP || XLENGTH(id_type_s) != 1 ||
        strcmp(CHAR(STRING_ELT(id_type_s, 0)), "numeric") != 0 ||
        TYPEOF(ukey_s) != REALSXP || TYPEOF(units_s) != VECSXP || TYPEOF(cal_s) != INTSXP ||
        (TYPEOF(rd->id_sexp) != REALSXP && TYPEOF(rd->id_sexp) != INTSXP) || Rf_isObject(rd->id_sexp)) {
      rd->sing_numeric_ok = 0;
      return 0;
    }
    SEXP nobs_s = cdid_sel_elt(units_s, "n_obs");
    SEXP ostart_s = cdid_sel_elt(units_s, "obs_start");
    if (TYPEOF(nobs_s) != INTSXP || (TYPEOF(ostart_s) != INTSXP && TYPEOF(ostart_s) != REALSXP)) return 0;
    rd->unit_id_key = REAL(ukey_s);
    rd->n_panel_units = (int) XLENGTH(ukey_s);
    rd->units_n_obs = INTEGER(nobs_s);
    rd->units_obs_start_s = ostart_s;
    rd->obs_cal_idx = INTEGER(cal_s);
  }
  return 1;
}

typedef struct {
  int has_ref;
  double N1;
  int n_recipes;
  int J;
  SEXP design_idx_s;
  SEXP weight_idx_s;
  const double *w_keys;
  int n_w_keys;
  const double *w_vals;

  int x_is_dated;
  const double *x_keys;
  int n_x_keys;
  const double *x_Z;
  int n_x_rows;
  int p_x;

  int s_is_dated;
  const double *s_keys;
  int n_s_keys;
  const double *s_Z;
  int n_s_rows;
  int p_s;

  const double *stat_a;      /* n_recipes x J or NULL */
  const double *stat_mass_U; /* n_recipes x J or NULL */
} cdid_sel_ref_t;

static int cdid_sel_init_ref(SEXP ref_sexp, cdid_sel_ref_t *rf) {
  memset(rf, 0, sizeof(*rf));
  if (ref_sexp == R_NilValue) return 1;
  if (TYPEOF(ref_sexp) != VECSXP) return 0;
  rf->has_ref = 1;
  double N_val = Rf_asReal(cdid_sel_elt(ref_sexp, "N"));
  if (!R_FINITE(N_val) || N_val < 1.0) return 0;
  rf->N1 = N_val + 1.0;

  SEXP recipes = cdid_sel_elt(ref_sexp, "recipes");
  if (TYPEOF(recipes) != VECSXP) return 0;
  rf->design_idx_s = cdid_sel_elt(recipes, "design_index");
  rf->weight_idx_s = cdid_sel_elt(recipes, "weight_index");
  if (rf->design_idx_s == R_NilValue || rf->weight_idx_s == R_NilValue) return 0;
  rf->n_recipes = (int) XLENGTH(rf->design_idx_s);

  SEXP wk_s = cdid_sel_elt(ref_sexp, "weight_keys");
  SEXP wv_s = cdid_sel_elt(ref_sexp, "w");
  if (TYPEOF(wk_s) != REALSXP || TYPEOF(wv_s) != REALSXP || XLENGTH(wk_s) != XLENGTH(wv_s)) return 0;
  rf->w_keys = REAL(wk_s);
  rf->n_w_keys = (int) XLENGTH(wk_s);
  rf->w_vals = REAL(wv_s);

  SEXP role_slot = cdid_sel_elt(ref_sexp, "role_slot");
  SEXP designs = cdid_sel_elt(ref_sexp, "designs");
  if (TYPEOF(role_slot) != STRSXP || TYPEOF(designs) != VECSXP) return 0;

  /* x design */
  SEXP x_slot = cdid_sel_elt(designs, "x");
  if (TYPEOF(x_slot) == VECSXP) {
    SEXP kind_s = cdid_sel_elt(x_slot, "kind");
    if (TYPEOF(kind_s) == STRSXP && XLENGTH(kind_s) == 1 && strcmp(CHAR(STRING_ELT(kind_s, 0)), "dated") == 0) {
      rf->x_is_dated = 1;
      SEXP k_s = cdid_sel_elt(x_slot, "keys");
      SEXP Z_s = cdid_sel_elt(x_slot, "Z");
      SEXP zdim = Rf_getAttrib(Z_s, R_DimSymbol);
      if (TYPEOF(k_s) != REALSXP || TYPEOF(Z_s) != REALSXP || TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return 0;
      rf->x_keys = REAL(k_s);
      rf->n_x_keys = (int) XLENGTH(k_s);
      rf->x_Z = REAL(Z_s);
      rf->n_x_rows = INTEGER(zdim)[0];
      rf->p_x = INTEGER(zdim)[1];
    } else {
      SEXP Z_s = cdid_sel_elt(x_slot, "static");
      SEXP zdim = Rf_getAttrib(Z_s, R_DimSymbol);
      if (TYPEOF(Z_s) != REALSXP || TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return 0;
      rf->x_is_dated = 0;
      rf->x_Z = REAL(Z_s);
      rf->n_x_rows = INTEGER(zdim)[0];
      rf->p_x = INTEGER(zdim)[1];
    }
  }

  /* s design (may share x slot if role_slot["s"] == "x") */
  const char *s_nm = NULL;
  SEXP rs_nms = Rf_getAttrib(role_slot, R_NamesSymbol);
  if (TYPEOF(rs_nms) == STRSXP) {
    for (R_xlen_t i = 0; i < XLENGTH(role_slot); ++i) {
      if (strcmp(CHAR(STRING_ELT(rs_nms, i)), "s") == 0) {
        s_nm = CHAR(STRING_ELT(role_slot, i));
        break;
      }
    }
  }
  if (s_nm != NULL) {
    SEXP s_slot = cdid_sel_elt(designs, s_nm);
    if (TYPEOF(s_slot) != VECSXP) return 0;
    SEXP kind_s = cdid_sel_elt(s_slot, "kind");
    if (TYPEOF(kind_s) == STRSXP && XLENGTH(kind_s) == 1 && strcmp(CHAR(STRING_ELT(kind_s, 0)), "dated") == 0) {
      rf->s_is_dated = 1;
      SEXP k_s = cdid_sel_elt(s_slot, "keys");
      SEXP Z_s = cdid_sel_elt(s_slot, "Z");
      SEXP zdim = Rf_getAttrib(Z_s, R_DimSymbol);
      if (TYPEOF(k_s) != REALSXP || TYPEOF(Z_s) != REALSXP || TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return 0;
      rf->s_keys = REAL(k_s);
      rf->n_s_keys = (int) XLENGTH(k_s);
      rf->s_Z = REAL(Z_s);
      rf->n_s_rows = INTEGER(zdim)[0];
      rf->p_s = INTEGER(zdim)[1];
    } else {
      SEXP Z_s = cdid_sel_elt(s_slot, "static");
      SEXP zdim = Rf_getAttrib(Z_s, R_DimSymbol);
      if (TYPEOF(Z_s) != REALSXP || TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return 0;
      rf->s_is_dated = 0;
      rf->s_Z = REAL(Z_s);
      rf->n_s_rows = INTEGER(zdim)[0];
      rf->p_s = INTEGER(zdim)[1];
    }
  }

  SEXP stats = cdid_sel_elt(ref_sexp, "stats");
  if (TYPEOF(stats) == VECSXP) {
    SEXP a_m = cdid_sel_elt(stats, "a");
    SEXP m_m = cdid_sel_elt(stats, "mass_U");
    SEXP adim = Rf_getAttrib(a_m, R_DimSymbol);
    if (TYPEOF(a_m) == REALSXP && TYPEOF(m_m) == REALSXP && TYPEOF(adim) == INTSXP && XLENGTH(adim) == 2 &&
        INTEGER(adim)[0] == rf->n_recipes) {
      rf->stat_a = REAL(a_m);
      rf->stat_mass_U = REAL(m_m);
      rf->J = INTEGER(adim)[1];
    }
  }
  return 1;
}

static inline int cdid_sel_ref_w_one(const cdid_sel_ref_t *rf, int recipe1, int row1, int *hint, double *w_out) {
  if (recipe1 < 1 || recipe1 > rf->n_recipes) return 0;
  double a_idx = cdid_sel_num_at(rf->weight_idx_s, recipe1 - 1);
  int pos = cdid_sel_bsearch_dbl(a_idx * rf->N1 + (double) row1, rf->w_keys, rf->n_w_keys, hint);
  if (pos < 0) return 0;
  *w_out = rf->w_vals[pos];
  return 1;
}

static inline int cdid_sel_ref_s_row_idx(const cdid_sel_ref_t *rf, int recipe1, int row1, int *hint) {
  if (!rf->s_is_dated) return (row1 >= 1 && row1 <= rf->n_s_rows) ? (row1 - 1) : -1;
  if (recipe1 < 1 || recipe1 > rf->n_recipes) return -1;
  double d_idx = cdid_sel_num_at(rf->design_idx_s, recipe1 - 1);
  return cdid_sel_bsearch_dbl(d_idx * rf->N1 + (double) row1, rf->s_keys, rf->n_s_keys, hint);
}

static inline int cdid_sel_ref_x_row_idx(const cdid_sel_ref_t *rf, int recipe1, int row1, int *hint) {
  if (!rf->x_is_dated) return (row1 >= 1 && row1 <= rf->n_x_rows) ? (row1 - 1) : -1;
  if (recipe1 < 1 || recipe1 > rf->n_recipes) return -1;
  double d_idx = cdid_sel_num_at(rf->design_idx_s, recipe1 - 1);
  return cdid_sel_bsearch_dbl(d_idx * rf->N1 + (double) row1, rf->x_keys, rf->n_x_keys, hint);
}

static inline int cdid_sel_ref_stat_one(const cdid_sel_ref_t *rf, int recipe1, int cohort1, int is_mass, double *out) {
  if (rf->stat_a == NULL || recipe1 < 1 || recipe1 > rf->n_recipes || cohort1 < 1 || cohort1 > rf->J) return 0;
  const double *mat = is_mass ? rf->stat_mass_U : rf->stat_a;
  double val = mat[(recipe1 - 1) + (cohort1 - 1) * rf->n_recipes];
  if (ISNAN(val)) return 0;
  *out = val;
  return 1;
}

static inline double cdid_sel_date_one(const cdid_sel_reader_t *rd, int row1, int dt1, int *hint_a, int *hint_u) {
  int ai = rd->a_index[row1 - 1];
  if (ai > 0) {
    if (dt1 >= 1 && dt1 <= rd->n_dates && rd->date_count[dt1 - 1] > 0) {
      int b0 = rd->date_start[dt1 - 1] - 1;
      int b1 = rd->date_end[dt1 - 1] - 1;
      int pos = cdid_sel_bsearch_int_range(ai, rd->post_a, b0, b1, hint_a);
      return (pos >= 0) ? 1.0 : 0.0;
    }
    return 0.0;
  }
  if (rd->n_singletons > 0) {
    double id_d = (TYPEOF(rd->id_sexp) == INTSXP) ? (double) INTEGER(rd->id_sexp)[row1 - 1] : REAL(rd->id_sexp)[row1 - 1];
    if (id_d == 0.0) id_d = 0.0;
    int u = cdid_sel_bsearch_dbl(id_d, rd->unit_id_key, rd->n_panel_units, hint_u);
    if (u >= 0 && rd->units_n_obs[u] == 1) {
      int opos = cdid_sel_int_at(rd->units_obs_start_s, u) - 1;
      if (opos >= 0 && rd->obs_cal_idx[opos] == dt1) return 1.0;
    }
  }
  return 0.0;
}

static inline int cdid_sel_pair_one(
    const cdid_sel_reader_t *rd,
    int row1,
    int ds1,
    int dt1,
    int read_dy,
    int *hint_s,
    int *hint_t,
    double *S_out,
    double *dy_out) {
  *S_out = 0.0;
  if (read_dy) *dy_out = 0.0;
  int ai = rd->a_index[row1 - 1];
  if (ai <= 0) return 1;
  if (ds1 < 1 || ds1 > rd->n_dates || dt1 < 1 || dt1 > rd->n_dates) return 1;
  if (rd->date_count[ds1 - 1] <= 0 || rd->date_count[dt1 - 1] <= 0) return 1;
  int ps = cdid_sel_bsearch_int_range(ai, rd->post_a, rd->date_start[ds1 - 1] - 1, rd->date_end[ds1 - 1] - 1, hint_s);
  if (ps < 0) return 1;
  int pt = cdid_sel_bsearch_int_range(ai, rd->post_a, rd->date_start[dt1 - 1] - 1, rd->date_end[dt1 - 1] - 1, hint_t);
  if (pt < 0) return 1;
  *S_out = 1.0;
  if (read_dy) {
    int os = rd->post_obs[ps] - 1;
    int ot = rd->post_obs[pt] - 1;
    if (os < 0 || (R_xlen_t) os >= rd->n_y || ot < 0 || (R_xlen_t) ot >= rd->n_y) return 0;
    double d = rd->y_obs[ot] - rd->y_obs[os];
    if (!R_FINITE(d)) return 0;
    *dy_out = d;
  }
  return 1;
}

/* ---------------------------------------------------------------------------
   3. C_cdid_sel_date_rows(rows_sexp, a_index_sexp, post_a_sexp, bounds2_sexp, sing_sexp, dt_sexp)
   --------------------------------------------------------------------------- */
SEXP C_cdid_sel_date_rows(
    SEXP rows_sexp,
    SEXP a_index_sexp,
    SEXP post_a_sexp,
    SEXP bounds2_sexp,
    SEXP sing_sexp,
    SEXP dt_sexp) {
  if (Rf_isObject(rows_sexp) || TYPEOF(a_index_sexp) != INTSXP || TYPEOF(post_a_sexp) != INTSXP) {
    return R_NilValue;
  }
  int tr = TYPEOF(rows_sexp);
  if (tr != INTSXP && tr != REALSXP) return R_NilValue;
  R_xlen_t m_xl = XLENGTH(rows_sexp);
  if (m_xl < 1 || m_xl > (R_xlen_t) INT_MAX) return R_NilValue;
  int m = (int) m_xl;
  int N = (int) XLENGTH(a_index_sexp);
  const int *a_idx = INTEGER(a_index_sexp);
  const int *post_a = INTEGER(post_a_sexp);
  int n_post = (int) XLENGTH(post_a_sexp);
  int dt1 = Rf_asInteger(dt_sexp);

  int has_bounds = 0, b0 = 0, b1 = -1;
  if (TYPEOF(bounds2_sexp) == INTSXP && XLENGTH(bounds2_sexp) == 2) {
    b0 = INTEGER(bounds2_sexp)[0] - 1;
    b1 = INTEGER(bounds2_sexp)[1] - 1;
    if (b0 >= 0 && b1 >= b0 && b1 < n_post) has_bounds = 1;
  }

  int has_sing = 0, n_panel_units = 0;
  SEXP id_s = R_NilValue, os_s = R_NilValue;
  const double *ukey = NULL;
  const int *u_nobs = NULL, *o_cal = NULL;
  if (sing_sexp != R_NilValue) {
    if (TYPEOF(sing_sexp) != VECSXP || XLENGTH(sing_sexp) != 5) return R_NilValue;
    id_s = VECTOR_ELT(sing_sexp, 0);
    SEXP uk_s = VECTOR_ELT(sing_sexp, 1);
    SEXP no_s = VECTOR_ELT(sing_sexp, 2);
    os_s = VECTOR_ELT(sing_sexp, 3);
    SEXP ci_s = VECTOR_ELT(sing_sexp, 4);
    if ((TYPEOF(id_s) != REALSXP && TYPEOF(id_s) != INTSXP) || Rf_isObject(id_s) || XLENGTH(id_s) != N ||
        TYPEOF(uk_s) != REALSXP || TYPEOF(no_s) != INTSXP ||
        (TYPEOF(os_s) != INTSXP && TYPEOF(os_s) != REALSXP) || TYPEOF(ci_s) != INTSXP) {
      return R_NilValue;
    }
    has_sing = 1;
    ukey = REAL(uk_s);
    n_panel_units = (int) XLENGTH(uk_s);
    u_nobs = INTEGER(no_s);
    o_cal = INTEGER(ci_s);
  }

  SEXP R_s = PROTECT(Rf_allocVector(REALSXP, m));
  double *R_out = REAL(R_s);
  int hint_a = -1, hint_u = -1;

  for (int i = 0; i < m; ++i) {
    int r1 = cdid_sel_int_at(rows_sexp, i);
    if (r1 == NA_INTEGER || r1 < 1 || r1 > N) {
      UNPROTECT(1);
      return R_NilValue;
    }
    int ai = a_idx[r1 - 1];
    double val = 0.0;
    if (ai > 0) {
      if (has_bounds && cdid_sel_bsearch_int_range(ai, post_a, b0, b1, &hint_a) >= 0) val = 1.0;
    } else if (has_sing) {
      double id_d = (TYPEOF(id_s) == INTSXP) ? (double) INTEGER(id_s)[r1 - 1] : REAL(id_s)[r1 - 1];
      if (id_d == 0.0) id_d = 0.0;
      int u = cdid_sel_bsearch_dbl(id_d, ukey, n_panel_units, &hint_u);
      if (u >= 0 && u_nobs[u] == 1) {
        int opos = cdid_sel_int_at(os_s, u) - 1;
        if (opos >= 0 && o_cal[opos] == dt1) val = 1.0;
      }
    }
    R_out[i] = val;
  }
  UNPROTECT(1);
  return R_s;
}

/* ---------------------------------------------------------------------------
   4. C_cdid_sel_fit_block(reader_sexp, cells_sexp, cids_sexp, dspec_sexp, ref_sexp, params_sexp)
      params_sexp: list(mode_int, is_date_int, eps_dbl, maxit_int)
      mode 1: pooled/date selection block (fit_selection_block_provider)
      mode 2: cp04_separate_pair block (B_m == 1 in fit_pooled_selection_block)
      mode 3: standalone cp04_selection_fit (cells_sexp is single spec list(pair, cohort, [reference]))
   --------------------------------------------------------------------------- */
SEXP C_cdid_sel_fit_block(
    SEXP reader_sexp,
    SEXP cells_sexp,
    SEXP cids_sexp,
    SEXP dspec_sexp,
    SEXP ref_sexp,
    SEXP params_sexp) {
  int mode = 0, is_date = 0, maxit = 0;
  double eps = 0.0;
  if (TYPEOF(params_sexp) == REALSXP && XLENGTH(params_sexp) == 4) {
    mode = (int) REAL(params_sexp)[0];
    is_date = (int) REAL(params_sexp)[1];
    eps = REAL(params_sexp)[2];
    maxit = (int) REAL(params_sexp)[3];
  } else if (TYPEOF(params_sexp) == VECSXP && XLENGTH(params_sexp) == 4) {
    mode = Rf_asInteger(VECTOR_ELT(params_sexp, 0));
    is_date = Rf_asInteger(VECTOR_ELT(params_sexp, 1));
    eps = Rf_asReal(VECTOR_ELT(params_sexp, 2));
    maxit = Rf_asInteger(VECTOR_ELT(params_sexp, 3));
  } else {
    return R_NilValue;
  }
  if ((mode < 1 || mode > 3) || !R_FINITE(eps) || eps <= 0.0 || maxit < 1) return R_NilValue;

  cdid_sel_reader_t rd;
  if (!cdid_sel_init_reader(reader_sexp, &rd, is_date)) return R_NilValue;
  cdid_sel_ref_t rf;
  if (!cdid_sel_init_ref(ref_sexp, &rf)) return R_NilValue;

  const void *vmax = vmaxget();

  /* Mode 2 or 3: single-cell full-design logit + H_chol */
  if (mode == 2 || mode == 3) {
    SEXP cell = cells_sexp;
    if (TYPEOF(cell) != VECSXP) { vmaxset(vmax); return R_NilValue; }
    int k = Rf_asInteger(cdid_sel_elt(cell, "cohort"));
    if (k < 1 || k > rd.J) { vmaxset(vmax); return R_NilValue; }
    SEXP ref_r_s = cdid_sel_elt(cell, "reference");
    int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
    SEXP pair = cdid_sel_elt(cell, "pair");
    int ds1 = Rf_asInteger(cdid_sel_elt(pair, "s_calendar_index"));
    int dt1 = Rf_asInteger(cdid_sel_elt(pair, "t_calendar_index"));
    if (ds1 == NA_INTEGER || dt1 == NA_INTEGER) { vmaxset(vmax); return R_NilValue; }

    double a_k = 0.0, mass_U_k = 0.0;
    if (!rf.has_ref) {
      a_k = rd.a_vec[k - 1];
      mass_U_k = rd.mass_U[k - 1];
    } else {
      if (!cdid_sel_ref_stat_one(&rf, r_id, k, 0, &a_k) ||
          !cdid_sel_ref_stat_one(&rf, r_id, k, 1, &mass_U_k)) {
        vmaxset(vmax); return R_NilValue;
      }
    }
    if (!R_FINITE(a_k) || a_k <= 0.0 || !R_FINITE(mass_U_k) || mass_U_k <= 0.0) {
      vmaxset(vmax); return R_NilValue;
    }

    int d = rf.has_ref ? rf.p_s : rd.d_Z;
    if (d < 1 || d > 128) { vmaxset(vmax); return R_NilValue; }

    int c0 = cdid_sel_int_at(rd.cohort_start_s, k - 1) - 1;
    int c1 = cdid_sel_int_at(rd.cohort_end_s, k - 1) - 1;
    int n_fit = 0;
    for (int idx = c0; idx <= c1; ++idx) {
      int row1 = rd.cohort_order[idx];
      if (rd.U_logical[row1 - 1] == 1) ++n_fit;
    }
    if (n_fit < d) { vmaxset(vmax); return R_NilValue; }

    double *X = (double *) R_alloc((size_t) n_fit * (size_t) d, sizeof(double));
    double *y = (double *) R_alloc((size_t) n_fit, sizeof(double));
    double *w = (double *) R_alloc((size_t) n_fit, sizeof(double));

    int pos = 0, hs = -1, ht = -1, hw = -1, hz = -1;
    long double wy_ld = 0.0L, sum_w_ld = 0.0L;
    int n_pos_w = 0;
    for (int idx = c0; idx <= c1; ++idx) {
      int row1 = rd.cohort_order[idx];
      if (rd.U_logical[row1 - 1] != 1) continue;
      double S_val = 0.0, dummy = 0.0;
      if (!cdid_sel_pair_one(&rd, row1, ds1, dt1, 0, &hs, &ht, &S_val, &dummy)) {
        vmaxset(vmax); return R_NilValue;
      }
      double wi = 0.0;
      if (!rf.has_ref) wi = rd.w_static[row1 - 1];
      else if (!cdid_sel_ref_w_one(&rf, r_id, row1, &hw, &wi)) { vmaxset(vmax); return R_NilValue; }
      if (!R_FINITE(wi) || wi <= 0.0) { vmaxset(vmax); return R_NilValue; }

      if (!rf.has_ref) {
        for (int j = 0; j < d; ++j) {
          double v = rd.Z_static[(row1 - 1) + (size_t) j * (size_t) rd.N];
          if (!R_FINITE(v)) { vmaxset(vmax); return R_NilValue; }
          X[pos + (size_t) j * (size_t) n_fit] = v;
        }
      } else {
        int zpos = cdid_sel_ref_s_row_idx(&rf, r_id, row1, &hz);
        if (zpos < 0) { vmaxset(vmax); return R_NilValue; }
        for (int j = 0; j < d; ++j) {
          double v = rf.s_Z[zpos + (size_t) j * (size_t) rf.n_s_rows];
          if (!R_FINITE(v)) { vmaxset(vmax); return R_NilValue; }
          X[pos + (size_t) j * (size_t) n_fit] = v;
        }
      }
      y[pos] = S_val;
      w[pos] = wi;
      wy_ld += (long double) (wi * S_val);
      sum_w_ld += (long double) wi;
      ++n_pos_w;
      ++pos;
    }
    double wy = (double) wy_ld;
    double sum_w = (double) sum_w_ld;
    if (n_pos_w < d || !R_FINITE(sum_w) || sum_w <= 0.0) { vmaxset(vmax); return R_NilValue; }

    /* Signed-coordinate separation check */
    for (int j = 0; j < d; ++j) {
      const double *xj = X + (size_t) j * (size_t) n_fit;
      int has_pos = 0, has_neg = 0;
      for (int i = 0; i < n_fit; ++i) {
        double v = xj[i];
        if (v != 0.0) {
          double sgn = (y[i] == 1.0) ? v : -v;
          if (sgn > 0.0) has_pos = 1;
          else if (sgn < 0.0) has_neg = 1;
        }
      }
      if ((has_pos && !has_neg) || (has_neg && !has_pos)) { vmaxset(vmax); return R_NilValue; }
    }

    /* Check exact column aliasing up front (if near-collinear, QR safe_tol catches it too) */
    double *eta = (double *) R_alloc((size_t) n_fit, sizeof(double));
    double *mu = (double *) R_alloc((size_t) n_fit, sizeof(double));
    double *wt = (double *) R_alloc((size_t) n_fit, sizeof(double));
    double *wz = (double *) R_alloc((size_t) n_fit, sizeof(double));
    double *W_X = (double *) R_alloc((size_t) n_fit * (size_t) d, sizeof(double));
    double *H = (double *) R_alloc((size_t) d * (size_t) d, sizeof(double));
    double *R_chol = (double *) R_alloc((size_t) d * (size_t) d, sizeof(double));
    double *beta = (double *) R_alloc((size_t) d, sizeof(double));

    long double devold_ld = 0.0L;
    const double ln2_val = log(2.0);
    for (int i = 0; i < n_fit; ++i) {
      eta[i] = 0.0;
      mu[i] = 0.5;
      devold_ld += (long double) (2.0 * w[i] * ln2_val);
    }
    memset(beta, 0, (size_t) d * sizeof(double));
    double devold = (double) devold_ld;
    int conv = 0, iter_used = 0;
    const double bnd_eps = 10.0 * DBL_EPSILON;
    const double qr_tol = (eps / 1000.0 < 1e-7) ? (eps / 1000.0) : 1e-7;
    BLAS_INT nfit_b = (BLAS_INT) n_fit, d_b = (BLAS_INT) d, one_b = 1;
    double one_d = 1.0, zero_d = 0.0;

    for (int iter = 1; iter <= maxit; ++iter) {
      iter_used = iter;
      for (int i = 0; i < n_fit; ++i) {
        double mui = mu[i];
        double varmu = mui * (1.0 - mui);
        double opexp = exp(eta[i]);
        double mu_eta = opexp / ((1.0 + opexp) * (1.0 + opexp));
        if (!(varmu > 0.0) || !(mu_eta > 0.0) || !R_FINITE(mu_eta)) { vmaxset(vmax); return R_NilValue; }
        double w_i = sqrt((w[i] * (mu_eta * mu_eta)) / varmu);
        double z_i = eta[i] + (y[i] - mui) / mu_eta;
        wt[i] = w_i;
        wz[i] = z_i * w_i;
      }
      for (int j = 0; j < d; ++j) {
        const double *xj = X + (size_t) j * (size_t) n_fit;
        double *wxj = W_X + (size_t) j * (size_t) n_fit;
        for (int i = 0; i < n_fit; ++i) wxj[i] = xj[i] * wt[i];
      }
      if (!cdid_sel_qr_wls_solve(W_X, wz, R_chol, beta, n_fit, d, qr_tol)) {
        vmaxset(vmax); return R_NilValue;
      }
      F77_CALL(dgemm)("N", "N", &nfit_b, &one_b, &d_b, &one_d,
                      X, &nfit_b, beta, &d_b, &zero_d, eta, &nfit_b FCONE FCONE);
      long double dev_ld = 0.0L;
      for (int i = 0; i < n_fit; ++i) {
        double ei = eta[i];
        if (!R_FINITE(ei) || fabs(ei) >= 30.0) { vmaxset(vmax); return R_NilValue; }
        double opexp = exp(ei);
        double mui = opexp / (1.0 + opexp);
        if (mui <= bnd_eps || mui >= 1.0 - bnd_eps) { vmaxset(vmax); return R_NilValue; }
        mu[i] = mui;
        double loss = (y[i] == 1.0) ? log(1.0 / mui) : log(1.0 / (1.0 - mui));
        dev_ld += (long double) (2.0 * w[i] * loss);
      }
      double dev = (double) dev_ld;
      if (!R_FINITE(dev)) { vmaxset(vmax); return R_NilValue; }
      if (fabs(dev - devold) / (0.1 + fabs(dev)) < eps) { conv = 1; break; }
      devold = dev;
    }
    if (!conv) { vmaxset(vmax); return R_NilValue; }

    long double q_sum_ld = 0.0L;
    for (int i = 0; i < n_fit; ++i) {
      double od = exp(eta[i]);
      double pr = od / (1.0 + od);
      if (pr <= bnd_eps || pr >= 1.0 - bnd_eps || !(pr * (1.0 - pr) > 0.0)) { vmaxset(vmax); return R_NilValue; }
      double qi = w[i] * pr * (1.0 - pr);
      wt[i] = qi;
      q_sum_ld += (long double) qi;
    }
    double q_sum = (double) q_sum_ld;
    for (int j = 0; j < d; ++j) {
      const double *xj = X + (size_t) j * (size_t) n_fit;
      double *wxj = W_X + (size_t) j * (size_t) n_fit;
      for (int i = 0; i < n_fit; ++i) wxj[i] = wt[i] * xj[i];
    }
    F77_CALL(dgemm)("T", "N", &d_b, &d_b, &nfit_b, &one_d,
                    X, &nfit_b, W_X, &nfit_b, &zero_d, H, &d_b FCONE FCONE);
    for (int idx = 0; idx < d * d; ++idx) {
      H[idx] /= rd.N_dbl;
      if (!R_FINITE(H[idx])) { vmaxset(vmax); return R_NilValue; }
    }
    if (!cdid_sel_chol_upper(H, R_chol, d, 0.0)) { vmaxset(vmax); return R_NilValue; }

    double pi_val = wy / (mass_U_k * rd.N_dbl);
    double r_val = pi_val / a_k;
    double m_b = q_sum / rd.N_dbl;

    /* Return list(beta, R_chol, iter_used, pi, r, m_b, z_bar_vec_or_null) */
    SEXP z_bar_s = R_NilValue;
    if (mode == 2 && TYPEOF(dspec_sexp) == VECSXP) {
      int has_int = Rf_asLogical(cdid_sel_elt(dspec_sexp, "has_intercept")) == 1;
      SEXP sc_s = cdid_sel_elt(dspec_sexp, "slope_cols");
      int d_s = (TYPEOF(sc_s) == INTSXP) ? (int) XLENGTH(sc_s) : 0;
      if (has_int && d_s > 0) {
        const int *sc = INTEGER(sc_s);
        z_bar_s = PROTECT(Rf_allocMatrix(REALSXP, 1, d_s));
        for (int j = 0; j < d_s; ++j) {
          const double *xj = X + (size_t) (sc[j] - 1) * (size_t) n_fit;
          long double zs = 0.0L;
          for (int i = 0; i < n_fit; ++i) zs += (long double) (wt[i] * xj[i]);
          REAL(z_bar_s)[j] = ((double) zs) / q_sum;
        }
      } else {
        z_bar_s = PROTECT(Rf_allocMatrix(REALSXP, 1, 0));
      }
    } else {
      z_bar_s = PROTECT(Rf_allocVector(REALSXP, 0));
    }

    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 7));
    SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, d));
    memcpy(REAL(beta_s), beta, (size_t) d * sizeof(double));
    SEXP R_s = PROTECT(Rf_allocMatrix(REALSXP, d, d));
    memcpy(REAL(R_s), R_chol, (size_t) d * (size_t) d * sizeof(double));
    SET_VECTOR_ELT(out_s, 0, beta_s);
    SET_VECTOR_ELT(out_s, 1, R_s);
    SET_VECTOR_ELT(out_s, 2, Rf_ScalarInteger(iter_used));
    SET_VECTOR_ELT(out_s, 3, Rf_ScalarReal(pi_val));
    SET_VECTOR_ELT(out_s, 4, Rf_ScalarReal(r_val));
    SET_VECTOR_ELT(out_s, 5, Rf_ScalarReal(m_b));
    SET_VECTOR_ELT(out_s, 6, z_bar_s);
    UNPROTECT(4);
    vmaxset(vmax);
    return out_s;
  }

  /* Mode 1: general pooled / date selection block (B_m >= 1) */
  if (TYPEOF(cells_sexp) != VECSXP || TYPEOF(cids_sexp) != INTSXP || TYPEOF(dspec_sexp) != VECSXP) {
    vmaxset(vmax); return R_NilValue;
  }
  int B_m = (int) XLENGTH(cids_sexp);
  const int *cids = INTEGER(cids_sexp);
  int has_int = Rf_asLogical(cdid_sel_elt(dspec_sexp, "has_intercept")) == 1;
  SEXP sc_s = cdid_sel_elt(dspec_sexp, "slope_cols");
  if (TYPEOF(sc_s) != INTSXP) { vmaxset(vmax); return R_NilValue; }
  int d_s = (int) XLENGTH(sc_s);
  const int *sc = INTEGER(sc_s);
  int p_m = (has_int ? B_m : 0) + d_s;
  if (B_m < 1 || (!has_int && d_s == 0) || d_s > 128) { vmaxset(vmax); return R_NilValue; }

  int max_z_cols = rf.has_ref ? rf.p_s : rd.d_Z;
  for (int j = 0; j < d_s; ++j) {
    if (sc[j] < 1 || sc[j] > max_z_cols) { vmaxset(vmax); return R_NilValue; }
  }

  int *cell_start = (int *) R_alloc((size_t) (B_m + 1), sizeof(int));
  cell_start[0] = 0;
  int n_cells_total = (int) XLENGTH(cells_sexp);
  for (int b = 0; b < B_m; ++b) {
    int cid = cids[b];
    if (cid < 1 || cid > n_cells_total) { vmaxset(vmax); return R_NilValue; }
    SEXP cell = VECTOR_ELT(cells_sexp, cid - 1);
    int k = Rf_asInteger(cdid_sel_elt(cell, "cohort"));
    if (k < 1 || k > rd.J) { vmaxset(vmax); return R_NilValue; }
    int c0 = cdid_sel_int_at(rd.cohort_start_s, k - 1) - 1;
    int c1 = cdid_sel_int_at(rd.cohort_end_s, k - 1) - 1;
    int nb = 0;
    for (int idx = c0; idx <= c1; ++idx) {
      if (rd.U_logical[rd.cohort_order[idx] - 1] == 1) ++nb;
    }
    if (nb < 1) { vmaxset(vmax); return R_NilValue; }
    cell_start[b + 1] = cell_start[b] + nb;
  }
  int n_fit = cell_start[B_m];
  if (n_fit < p_m) { vmaxset(vmax); return R_NilValue; }

  double *D_mat = (d_s > 0) ? (double *) R_alloc((size_t) n_fit * (size_t) d_s, sizeof(double)) : NULL;
  double *y_all = (double *) R_alloc((size_t) n_fit, sizeof(double));
  double *w_all = (double *) R_alloc((size_t) n_fit, sizeof(double));
  double *pi_vec = (double *) R_alloc((size_t) B_m, sizeof(double));
  double *r_vec = (double *) R_alloc((size_t) B_m, sizeof(double));

  int *pos_sgn = (d_s > 0) ? (int *) R_alloc((size_t) d_s, sizeof(int)) : NULL;
  int *neg_sgn = (d_s > 0) ? (int *) R_alloc((size_t) d_s, sizeof(int)) : NULL;
  if (d_s > 0) {
    memset(pos_sgn, 0, (size_t) d_s * sizeof(int));
    memset(neg_sgn, 0, (size_t) d_s * sizeof(int));
  }

  double sum_w_total = 0.0;
  for (int b = 0; b < B_m; ++b) {
    SEXP cell = VECTOR_ELT(cells_sexp, cids[b] - 1);
    int k = Rf_asInteger(cdid_sel_elt(cell, "cohort"));
    SEXP ref_r_s = cdid_sel_elt(cell, "reference");
    int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

    double a_k = 0.0, mass_U_k = 0.0;
    if (!rf.has_ref) {
      a_k = rd.a_vec[k - 1];
      mass_U_k = rd.mass_U[k - 1];
    } else {
      if (!cdid_sel_ref_stat_one(&rf, r_id, k, 0, &a_k) ||
          !cdid_sel_ref_stat_one(&rf, r_id, k, 1, &mass_U_k)) {
        vmaxset(vmax); return R_NilValue;
      }
    }
    if (!R_FINITE(a_k) || a_k <= 0.0 || !R_FINITE(mass_U_k) || mass_U_k <= 0.0) {
      vmaxset(vmax); return R_NilValue;
    }

    int ds1 = 0, dt1 = 0;
    if (is_date) {
      SEXP dspec = cdid_sel_elt(cell, "date_spec");
      dt1 = Rf_asInteger(cdid_sel_elt(dspec, "calendar_index"));
      if (dt1 == NA_INTEGER || dt1 < 1 || dt1 > rd.n_dates) { vmaxset(vmax); return R_NilValue; }
    } else {
      SEXP pair = cdid_sel_elt(cell, "pair");
      ds1 = Rf_asInteger(cdid_sel_elt(pair, "s_calendar_index"));
      dt1 = Rf_asInteger(cdid_sel_elt(pair, "t_calendar_index"));
      if (ds1 == NA_INTEGER || dt1 == NA_INTEGER) { vmaxset(vmax); return R_NilValue; }
    }

    int c0 = cdid_sel_int_at(rd.cohort_start_s, k - 1) - 1;
    int c1 = cdid_sel_int_at(rd.cohort_end_s, k - 1) - 1;
    int pos = cell_start[b];
    int hs = -1, ht = -1, hu = -1, hw = -1, hz = -1;
    double wy = 0.0;
    int all_one = 1, all_zero = 1;

    for (int idx = c0; idx <= c1; ++idx) {
      int row1 = rd.cohort_order[idx];
      if (rd.U_logical[row1 - 1] != 1) continue;
      double yi = 0.0, dummy = 0.0;
      if (is_date) {
        yi = cdid_sel_date_one(&rd, row1, dt1, &hs, &hu);
      } else {
        if (!cdid_sel_pair_one(&rd, row1, ds1, dt1, 0, &hs, &ht, &yi, &dummy)) {
          vmaxset(vmax); return R_NilValue;
        }
      }
      if (yi == 1.0) all_zero = 0;
      else all_one = 0;

      double wi = 0.0;
      if (!rf.has_ref) wi = rd.w_static[row1 - 1];
      else if (!cdid_sel_ref_w_one(&rf, r_id, row1, &hw, &wi)) { vmaxset(vmax); return R_NilValue; }
      if (!R_FINITE(wi) || wi <= 0.0) { vmaxset(vmax); return R_NilValue; }

      if (d_s > 0) {
        int zpos = rf.has_ref ? cdid_sel_ref_s_row_idx(&rf, r_id, row1, &hz) : (row1 - 1);
        if (zpos < 0) { vmaxset(vmax); return R_NilValue; }
        const double *Z_base = rf.has_ref ? rf.s_Z : rd.Z_static;
        int z_ld = rf.has_ref ? rf.n_s_rows : rd.N;
        for (int j = 0; j < d_s; ++j) {
          double v = Z_base[zpos + (size_t) (sc[j] - 1) * (size_t) z_ld];
          if (!R_FINITE(v)) { vmaxset(vmax); return R_NilValue; }
          D_mat[pos + (size_t) j * (size_t) n_fit] = v;
          double sgn = (yi == 1.0) ? v : -v;
          if (sgn > 0.0) pos_sgn[j] = 1;
          else if (sgn < 0.0) neg_sgn[j] = 1;
        }
      }

      y_all[pos] = yi;
      w_all[pos] = wi;
      wy += wi * yi;
      sum_w_total += wi;
      ++pos;
    }

    if (has_int && (all_one || all_zero)) {
      /* Separation witness: fall back to R for exact structured diagnostic */
      vmaxset(vmax); return R_NilValue;
    }
    pi_vec[b] = wy / (mass_U_k * rd.N_dbl);
    r_vec[b] = pi_vec[b] / a_k;
  }

  for (int j = 0; j < d_s; ++j) {
    if ((pos_sgn[j] && !neg_sgn[j]) || (neg_sgn[j] && !pos_sgn[j])) {
      vmaxset(vmax); return R_NilValue;
    }
  }

  double rank_tol = eps / 1000.0;
  if (rank_tol > 1e-7) rank_tol = 1e-7;
  if (rank_tol < 1e-12) rank_tol = 1e-12;
  double safe_tol = (rank_tol > 1e-6) ? rank_tol : 1e-6;
  const double eps_prob = 10.0 * DBL_EPSILON;

  double *alpha = has_int ? (double *) R_alloc((size_t) B_m, sizeof(double)) : NULL;
  double *beta = (d_s > 0) ? (double *) R_alloc((size_t) d_s, sizeof(double)) : NULL;
  if (has_int) memset(alpha, 0, (size_t) B_m * sizeof(double));
  if (d_s > 0) memset(beta, 0, (size_t) d_s * sizeof(double));

  double *eta = (double *) R_alloc((size_t) n_fit, sizeof(double));
  double *q_vec = (double *) R_alloc((size_t) n_fit, sizeof(double));
  double *z_work = (double *) R_alloc((size_t) n_fit, sizeof(double));
  double *W_D = (d_s > 0) ? (double *) R_alloc((size_t) n_fit * (size_t) d_s, sizeof(double)) : NULL;
  double *wz = (d_s > 0) ? (double *) R_alloc((size_t) n_fit, sizeof(double)) : NULL;
  double *R_mat = (d_s > 0) ? (double *) R_alloc((size_t) d_s * (size_t) d_s, sizeof(double)) : NULL;
  double *Q_cell = (double *) R_alloc((size_t) B_m, sizeof(double));
  double *ybar_cell = has_int ? (double *) R_alloc((size_t) B_m, sizeof(double)) : NULL;
  double *zbar_cell = (has_int && d_s > 0) ? (double *) R_alloc((size_t) B_m * (size_t) d_s, sizeof(double)) : NULL;

  memset(eta, 0, (size_t) n_fit * sizeof(double));
  double dev_old = 2.0 * M_LN2 * sum_w_total;
  int conv = 0, iter_used = 0;

  for (int iter = 1; iter <= maxit; ++iter) {
    iter_used = iter;
    for (int b = 0; b < B_m; ++b) {
      int s0 = cell_start[b], s1 = cell_start[b + 1];
      double qb_hi = 0.0, qb_lo = 0.0, qz_hi = 0.0, qz_lo = 0.0;
      for (int i = s0; i < s1; ++i) {
        double ei = eta[i];
        if (!R_FINITE(ei) || fabs(ei) > 30.0) { vmaxset(vmax); return R_NilValue; }
        double odds = exp(ei);
        double mu = (ei >= 0.0) ? (1.0 / (1.0 + exp(-ei))) : (odds / (1.0 + odds));
        if (mu <= eps_prob || mu >= 1.0 - eps_prob) { vmaxset(vmax); return R_NilValue; }
        double mu_eta = mu * (1.0 - mu);
        double qi = w_all[i] * mu_eta;
        double zi = ei + (y_all[i] - mu) / mu_eta;
        if (!R_FINITE(qi) || !(qi > 0.0) || !R_FINITE(zi)) { vmaxset(vmax); return R_NilValue; }
        q_vec[i] = qi;
        z_work[i] = zi;
        cdid_dd_add(&qb_hi, &qb_lo, qi);
        if (has_int) cdid_dd_prod_add(&qz_hi, &qz_lo, qi, zi);
      }
      double qb_sum = qb_hi + qb_lo;
      if (!(qb_sum > 0.0) || !R_FINITE(qb_sum)) { vmaxset(vmax); return R_NilValue; }
      Q_cell[b] = qb_sum;
      if (has_int) {
        ybar_cell[b] = (qz_hi + qz_lo) / qb_sum;
        if (d_s > 0) {
          for (int j = 0; j < d_s; ++j) {
            const double *dj = D_mat + (size_t) j * (size_t) n_fit;
            double s_hi = 0.0, s_lo = 0.0;
            for (int i = s0; i < s1; ++i) cdid_dd_prod_add(&s_hi, &s_lo, q_vec[i], dj[i]);
            zbar_cell[b + (size_t) j * (size_t) B_m] = (s_hi + s_lo) / qb_sum;
          }
        }
      }
    }

    if (d_s == 0) {
      for (int b = 0; b < B_m; ++b) {
        if (!R_FINITE(ybar_cell[b])) { vmaxset(vmax); return R_NilValue; }
        alpha[b] = ybar_cell[b];
      }
    } else if (has_int) {
      for (int b = 0; b < B_m; ++b) {
        int s0 = cell_start[b], s1 = cell_start[b + 1];
        double yb = ybar_cell[b];
        for (int i = s0; i < s1; ++i) {
          double sq = sqrt(q_vec[i]);
          wz[i] = sq * (z_work[i] - yb);
        }
        for (int j = 0; j < d_s; ++j) {
          double zbj = zbar_cell[b + (size_t) j * (size_t) B_m];
          const double *dj = D_mat + (size_t) j * (size_t) n_fit;
          double *wdj = W_D + (size_t) j * (size_t) n_fit;
          for (int i = s0; i < s1; ++i) wdj[i] = sqrt(q_vec[i]) * (dj[i] - zbj);
        }
      }
      if (!cdid_sel_qr_wls_solve_dd(W_D, wz, R_mat, beta, n_fit, d_s, safe_tol)) {
        vmaxset(vmax); return R_NilValue;
      }
      for (int b = 0; b < B_m; ++b) {
        double ab = ybar_cell[b];
        for (int j = 0; j < d_s; ++j) ab -= zbar_cell[b + (size_t) j * (size_t) B_m] * beta[j];
        if (!R_FINITE(ab)) { vmaxset(vmax); return R_NilValue; }
        alpha[b] = ab;
      }
    } else {
      for (int i = 0; i < n_fit; ++i) wz[i] = sqrt(q_vec[i]) * z_work[i];
      for (int j = 0; j < d_s; ++j) {
        const double *dj = D_mat + (size_t) j * (size_t) n_fit;
        double *wdj = W_D + (size_t) j * (size_t) n_fit;
        for (int i = 0; i < n_fit; ++i) wdj[i] = sqrt(q_vec[i]) * dj[i];
      }
      if (!cdid_sel_qr_wls_solve_dd(W_D, wz, R_mat, beta, n_fit, d_s, safe_tol)) {
        vmaxset(vmax); return R_NilValue;
      }
    }

    /* Update eta and compute deviance */
    for (int b = 0; b < B_m; ++b) {
      int s0 = cell_start[b], s1 = cell_start[b + 1];
      double ab = has_int ? alpha[b] : 0.0;
      for (int i = s0; i < s1; ++i) eta[i] = ab;
    }
    for (int j = 0; j < d_s; ++j) {
      double bj = beta[j];
      const double *dj = D_mat + (size_t) j * (size_t) n_fit;
      for (int i = 0; i < n_fit; ++i) eta[i] += dj[i] * bj;
    }
    double dev_hi = 0.0, dev_lo = 0.0;
    for (int i = 0; i < n_fit; ++i) {
      double ei = eta[i];
      if (!R_FINITE(ei) || fabs(ei) > 30.0) { vmaxset(vmax); return R_NilValue; }
      double odds = exp(ei);
      double mu = (ei >= 0.0) ? (1.0 / (1.0 + exp(-ei))) : (odds / (1.0 + odds));
      if (mu <= eps_prob || mu >= 1.0 - eps_prob) { vmaxset(vmax); return R_NilValue; }
      double loss = (y_all[i] == 1.0) ? -log(mu) : -log1p(-mu);
      cdid_dd_add(&dev_hi, &dev_lo, 2.0 * w_all[i] * loss);
    }
    double dev = dev_hi + dev_lo;
    if (!R_FINITE(dev)) { vmaxset(vmax); return R_NilValue; }
    if (fabs(dev - dev_old) / (0.1 + fabs(dev)) < eps) { conv = 1; break; }
    dev_old = dev;
  }

  if (!conv) { vmaxset(vmax); return R_NilValue; }

  /* Final Fisher info evaluation (final = TRUE) */
  for (int b = 0; b < B_m; ++b) {
    int s0 = cell_start[b], s1 = cell_start[b + 1];
    double qb_hi = 0.0, qb_lo = 0.0;
    for (int i = s0; i < s1; ++i) {
      double ei = eta[i];
      double odds = exp(ei);
      double mu = odds / (1.0 + odds);
      if (!R_FINITE(ei) || !R_FINITE(odds) || !(odds > 0.0) || !R_FINITE(mu) ||
          mu <= eps_prob || mu >= 1.0 - eps_prob || !(mu * (1.0 - mu) > 0.0)) {
        vmaxset(vmax); return R_NilValue;
      }
      double qi = w_all[i] * mu * (1.0 - mu);
      if (!R_FINITE(qi) || !(qi > 0.0)) { vmaxset(vmax); return R_NilValue; }
      q_vec[i] = qi;
      cdid_dd_add(&qb_hi, &qb_lo, qi);
    }
    double qb_sum = qb_hi + qb_lo;
    if (!R_FINITE(qb_sum) || !(qb_sum > 0.0)) { vmaxset(vmax); return R_NilValue; }
    Q_cell[b] = qb_sum;
    if (has_int && d_s > 0) {
      for (int j = 0; j < d_s; ++j) {
        const double *dj = D_mat + (size_t) j * (size_t) n_fit;
        double s_hi = 0.0, s_lo = 0.0;
        for (int i = s0; i < s1; ++i) cdid_dd_prod_add(&s_hi, &s_lo, q_vec[i], dj[i]);
        double zb = (s_hi + s_lo) / qb_sum;
        if (!R_FINITE(zb)) { vmaxset(vmax); return R_NilValue; }
        zbar_cell[b + (size_t) j * (size_t) B_m] = zb;
      }
    }
  }

  double cutoff = 0.0;
  if (d_s > 0) {
    if (has_int) {
      for (int b = 0; b < B_m; ++b) {
        int s0 = cell_start[b], s1 = cell_start[b + 1];
        for (int j = 0; j < d_s; ++j) {
          double zbj = zbar_cell[b + (size_t) j * (size_t) B_m];
          const double *dj = D_mat + (size_t) j * (size_t) n_fit;
          double *wdj = W_D + (size_t) j * (size_t) n_fit;
          for (int i = s0; i < s1; ++i) wdj[i] = sqrt(q_vec[i]) * (dj[i] - zbj);
        }
      }
    } else {
      for (int j = 0; j < d_s; ++j) {
        const double *dj = D_mat + (size_t) j * (size_t) n_fit;
        double *wdj = W_D + (size_t) j * (size_t) n_fit;
        for (int i = 0; i < n_fit; ++i) wdj[i] = sqrt(q_vec[i]) * dj[i];
      }
    }
    if (!cdid_sel_qr_wls_solve_dd(W_D, NULL, R_mat, NULL, n_fit, d_s, safe_tol)) {
      vmaxset(vmax); return R_NilValue;
    }
    double max_diag = 0.0;
    for (int j = 0; j < d_s; ++j) {
      double ad = fabs(R_mat[j + j * d_s]);
      if (ad > max_diag) max_diag = ad;
    }
    cutoff = rank_tol * max_diag;
  }

  /* Pack result: list(alpha, beta, m_b, z_bar, R_or_null, pi, r, iter_used, n_fit, rank_tol, cutoff) */
  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 11));
  SEXP alpha_s = PROTECT(Rf_allocVector(REALSXP, has_int ? B_m : 0));
  if (has_int) memcpy(REAL(alpha_s), alpha, (size_t) B_m * sizeof(double));
  SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, d_s));
  if (d_s > 0) memcpy(REAL(beta_s), beta, (size_t) d_s * sizeof(double));
  SEXP mb_s = PROTECT(Rf_allocVector(REALSXP, B_m));
  for (int b = 0; b < B_m; ++b) REAL(mb_s)[b] = Q_cell[b] / rd.N_dbl;
  SEXP zb_s = PROTECT(Rf_allocMatrix(REALSXP, B_m, (has_int && d_s > 0) ? d_s : 0));
  if (has_int && d_s > 0) memcpy(REAL(zb_s), zbar_cell, (size_t) B_m * (size_t) d_s * sizeof(double));
  SEXP R_s = R_NilValue;
  int n_prot = 5;
  if (d_s > 0) {
    R_s = PROTECT(Rf_allocMatrix(REALSXP, d_s, d_s));
    memcpy(REAL(R_s), R_mat, (size_t) d_s * (size_t) d_s * sizeof(double));
    ++n_prot;
  }
  SEXP pi_s = PROTECT(Rf_allocVector(REALSXP, B_m));
  memcpy(REAL(pi_s), pi_vec, (size_t) B_m * sizeof(double));
  SEXP r_s = PROTECT(Rf_allocVector(REALSXP, B_m));
  memcpy(REAL(r_s), r_vec, (size_t) B_m * sizeof(double));
  n_prot += 2;

  SET_VECTOR_ELT(out_s, 0, alpha_s);
  SET_VECTOR_ELT(out_s, 1, beta_s);
  SET_VECTOR_ELT(out_s, 2, mb_s);
  SET_VECTOR_ELT(out_s, 3, zb_s);
  SET_VECTOR_ELT(out_s, 4, R_s);
  SET_VECTOR_ELT(out_s, 5, pi_s);
  SET_VECTOR_ELT(out_s, 6, r_s);
  SET_VECTOR_ELT(out_s, 7, Rf_ScalarInteger(iter_used));
  SET_VECTOR_ELT(out_s, 8, Rf_ScalarReal((double) n_fit));
  SET_VECTOR_ELT(out_s, 9, Rf_ScalarReal(rank_tol));
  SET_VECTOR_ELT(out_s, 10, Rf_ScalarReal(cutoff));
  UNPROTECT(n_prot);
  vmaxset(vmax);
  return out_s;
}

/* ---------------------------------------------------------------------------
   Internal Selection Prediction & Upstream Evaluation Engine
   --------------------------------------------------------------------------- */
typedef struct {
  int is_date_model;
  int d_full;
  int d_s;
  const int *slope_cols; /* 1-based length d_s */
  SEXP cells;            /* VECSXP */
  SEXP blocks;           /* VECSXP */
  SEXP block_specs;      /* VECSXP */
  SEXP pair_cohorts;     /* VECSXP */
  int n_cells;
  int n_blocks;
  int n_pc;
} cdid_sel_model_t;

static int cdid_sel_init_model(SEXP sel_fit_sexp, cdid_sel_model_t *sm) {
  memset(sm, 0, sizeof(*sm));
  if (TYPEOF(sel_fit_sexp) != VECSXP) return 0;
  SEXP sm_s = cdid_sel_elt(sel_fit_sexp, "selection_model");
  if (TYPEOF(sm_s) != STRSXP || XLENGTH(sm_s) != 1) return 0;
  sm->is_date_model = (strcmp(CHAR(STRING_ELT(sm_s, 0)), "independent_period_logit") == 0) ? 1 : 0;
  SEXP dspec = cdid_sel_elt(sel_fit_sexp, "design_spec");
  if (TYPEOF(dspec) != VECSXP) return 0;
  sm->d_full = Rf_asInteger(cdid_sel_elt(dspec, "d"));
  SEXP sc_s = cdid_sel_elt(dspec, "slope_cols");
  if (TYPEOF(sc_s) != INTSXP) return 0;
  sm->d_s = (int) XLENGTH(sc_s);
  sm->slope_cols = INTEGER(sc_s);
  sm->cells = cdid_sel_elt(sel_fit_sexp, "cells");
  sm->blocks = cdid_sel_elt(sel_fit_sexp, "blocks");
  sm->block_specs = cdid_sel_elt(sel_fit_sexp, "block_specs");
  sm->pair_cohorts = cdid_sel_elt(sel_fit_sexp, "pair_cohorts");
  if (TYPEOF(sm->cells) != VECSXP || TYPEOF(sm->blocks) != VECSXP ||
      TYPEOF(sm->block_specs) != VECSXP || TYPEOF(sm->pair_cohorts) != VECSXP) return 0;
  sm->n_cells = (int) XLENGTH(sm->cells);
  sm->n_blocks = (int) XLENGTH(sm->blocks);
  sm->n_pc = (int) XLENGTH(sm->pair_cohorts);
  return 1;
}

static int cdid_sel_predict_cell_one(
    const cdid_sel_model_t *sm,
    const cdid_sel_reader_t *rd,
    const cdid_sel_ref_t *rf,
    int cell_id1,
    int row1,
    int *hint_z,
    double *prob_out) {
  if (cell_id1 < 1 || cell_id1 > sm->n_cells) return 0;
  SEXP cell = VECTOR_ELT(sm->cells, cell_id1 - 1);
  int block_id1 = Rf_asInteger(cdid_sel_elt(cell, "block_id"));
  int local_b1 = Rf_asInteger(cdid_sel_elt(cell, "local_b"));
  if (block_id1 < 1 || block_id1 > sm->n_blocks) return 0;
  SEXP block = VECTOR_ELT(sm->blocks, block_id1 - 1);
  SEXP branch_s = cdid_sel_elt(block, "branch");
  if (TYPEOF(branch_s) != STRSXP || XLENGTH(branch_s) != 1) return 0;
  const char *br = CHAR(STRING_ELT(branch_s, 0));

  if (strcmp(br, "zero_col") == 0) {
    *prob_out = 0.5;
    return 1;
  }

  SEXP ref_r_s = cdid_sel_elt(cell, "reference");
  int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
  int zpos = rf->has_ref ? cdid_sel_ref_s_row_idx(rf, r_id, row1, hint_z) : (row1 - 1);
  if (zpos < 0) return 0;
  const double *Z_base = rf->has_ref ? rf->s_Z : rd->Z_static;
  int z_ld = rf->has_ref ? rf->n_s_rows : rd->N;
  int z_cols = rf->has_ref ? rf->p_s : rd->d_Z;

  double eta = 0.0;
  if (strcmp(br, "cp04_separate_pair") == 0) {
    SEXP coef_s = cdid_sel_elt(block, "coefficients");
    if (TYPEOF(coef_s) != REALSXP || XLENGTH(coef_s) != z_cols) return 0;
    const double *coef = REAL(coef_s);
    for (int j = 0; j < z_cols; ++j) eta += Z_base[zpos + (size_t) j * (size_t) z_ld] * coef[j];
  } else {
    int has_int = Rf_asLogical(cdid_sel_elt(block, "has_intercept")) == 1;
    double alpha_val = 0.0;
    if (has_int) {
      SEXP al_s = cdid_sel_elt(block, "alpha");
      if (TYPEOF(al_s) != REALSXP || local_b1 < 1 || local_b1 > XLENGTH(al_s)) return 0;
      alpha_val = REAL(al_s)[local_b1 - 1];
    }
    double dot_zb = 0.0;
    if (sm->d_s > 0) {
      SEXP be_s = cdid_sel_elt(block, "beta");
      if (TYPEOF(be_s) != REALSXP || XLENGTH(be_s) != sm->d_s) return 0;
      const double *be = REAL(be_s);
      for (int j = 0; j < sm->d_s; ++j) {
        dot_zb += Z_base[zpos + (size_t) (sm->slope_cols[j] - 1) * (size_t) z_ld] * be[j];
      }
    }
    eta = alpha_val + dot_zb;
  }

  double p = 0.0, od = 0.0, act = 0.0;
  cdid_sel_clip_logit_scalar(eta, &p, &od, &act);
  if (!R_FINITE(p) || !(p > 0.0) || p >= 1.0) return 0;
  *prob_out = p;
  return 1;
}

static int cdid_sel_predict_pair_one(
    const cdid_sel_model_t *sm,
    const cdid_sel_reader_t *rd,
    const cdid_sel_ref_t *rf,
    int q1,
    int row1,
    int *hint_z1,
    int *hint_z2,
    int *hint_w,
    double *prob_out,
    double *ps_out,
    double *pt_out,
    double *w_sel_out,
    double *a_pc_out) {
  if (q1 < 1 || q1 > sm->n_pc) return 0;
  SEXP pc = VECTOR_ELT(sm->pair_cohorts, q1 - 1);
  int coh1 = Rf_asInteger(cdid_sel_elt(pc, "cohort"));
  SEXP ref_r_s = cdid_sel_elt(pc, "reference");
  int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

  double p_pair = 0.0, ps = 0.0, pt = 0.0;
  if (!sm->is_date_model) {
    int cid1 = Rf_asInteger(cdid_sel_elt(pc, "pair_cell_id"));
    if (!cdid_sel_predict_cell_one(sm, rd, rf, cid1, row1, hint_z1, &p_pair)) return 0;
  } else {
    int scid1 = Rf_asInteger(cdid_sel_elt(pc, "s_cell_id"));
    int tcid1 = Rf_asInteger(cdid_sel_elt(pc, "t_cell_id"));
    if (!cdid_sel_predict_cell_one(sm, rd, rf, scid1, row1, hint_z1, &ps)) return 0;
    if (!cdid_sel_predict_cell_one(sm, rd, rf, tcid1, row1, hint_z2, &pt)) return 0;
    p_pair = ps * pt;
    if (!R_FINITE(p_pair) || !(p_pair > 0.0) || p_pair >= 1.0) return 0;
  }

  double a_pc = 0.0, w_row = 0.0;
  if (!rf->has_ref) {
    a_pc = rd->a_vec[coh1 - 1];
    w_row = rd->w_static[row1 - 1];
  } else {
    if (!cdid_sel_ref_stat_one(rf, r_id, coh1, 0, &a_pc)) return 0;
    if (!cdid_sel_ref_w_one(rf, r_id, row1, hint_w, &w_row)) return 0;
  }
  double p_weight = (p_pair > CDID_SEL_CLIP_EPS) ? p_pair : CDID_SEL_CLIP_EPS;
  double rho_i = rd->rho[row1 - 1];
  double w_sel = w_row * a_pc / (rho_i * p_weight);
  if (!R_FINITE(w_sel) || !(w_sel > 0.0)) return 0;

  *prob_out = p_pair;
  if (ps_out) *ps_out = ps;
  if (pt_out) *pt_out = pt;
  if (w_sel_out) *w_sel_out = w_sel;
  if (a_pc_out) *a_pc_out = a_pc;
  return 1;
}

static int cdid_sel_predict_cell_batch(
    const cdid_sel_model_t *sm,
    const cdid_sel_reader_t *rd,
    const cdid_sel_ref_t *rf,
    int cell_id1,
    const int *rows1_sub,
    int n_sub,
    double *prob_out) {
  if (n_sub < 0 || cell_id1 < 1 || cell_id1 > sm->n_cells) return 0;
  if (n_sub == 0) return 1;
  SEXP cell = VECTOR_ELT(sm->cells, cell_id1 - 1);
  int block_id1 = Rf_asInteger(cdid_sel_elt(cell, "block_id"));
  int local_b1 = Rf_asInteger(cdid_sel_elt(cell, "local_b"));
  if (block_id1 < 1 || block_id1 > sm->n_blocks) return 0;
  SEXP block = VECTOR_ELT(sm->blocks, block_id1 - 1);
  SEXP branch_s = cdid_sel_elt(block, "branch");
  if (TYPEOF(branch_s) != STRSXP || XLENGTH(branch_s) != 1) return 0;
  const char *br = CHAR(STRING_ELT(branch_s, 0));

  if (strcmp(br, "zero_col") == 0) {
    for (int i = 0; i < n_sub; ++i) prob_out[i] = 0.5;
    return 1;
  }

  SEXP ref_r_s = cdid_sel_elt(cell, "reference");
  int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
  const double *Z_base = rf->has_ref ? rf->s_Z : rd->Z_static;
  int z_ld = rf->has_ref ? rf->n_s_rows : rd->N;
  int z_cols = rf->has_ref ? rf->p_s : rd->d_Z;
  if (Z_base == NULL) return 0;

  const void *vm_cb = vmaxget();
  double *eta_sub = (double *) R_alloc((size_t) n_sub, sizeof(double));
  int hz = -1;

  if (strcmp(br, "cp04_separate_pair") == 0) {
    SEXP coef_s = cdid_sel_elt(block, "coefficients");
    if (TYPEOF(coef_s) != REALSXP || XLENGTH(coef_s) != z_cols || z_cols < 1) {
      vmaxset(vm_cb);
      return 0;
    }
    const double *coef = REAL(coef_s);
    double *Z_sub = (double *) R_alloc((size_t) n_sub * (size_t) z_cols, sizeof(double));
    for (int i = 0; i < n_sub; ++i) {
      int r1 = rows1_sub[i];
      int zpos = rf->has_ref ? cdid_sel_ref_s_row_idx(rf, r_id, r1, &hz) : (r1 - 1);
      if (zpos < 0) { vmaxset(vm_cb); return 0; }
      for (int j = 0; j < z_cols; ++j) {
        double v = Z_base[zpos + (size_t) j * (size_t) z_ld];
        if (!R_FINITE(v)) { vmaxset(vm_cb); return 0; }
        Z_sub[i + (size_t) j * (size_t) n_sub] = v;
      }
    }
    BLAS_INT m_b = (BLAS_INT) n_sub, n_b = 1, k_b = (BLAS_INT) z_cols;
    double one_d = 1.0, zero_d = 0.0;
    F77_CALL(dgemm)("N", "N", &m_b, &n_b, &k_b, &one_d, Z_sub, &m_b, coef, &k_b, &zero_d, eta_sub, &m_b FCONE FCONE);
  } else {
    int has_int = Rf_asLogical(cdid_sel_elt(block, "has_intercept")) == 1;
    double alpha_val = 0.0;
    if (has_int) {
      SEXP al_s = cdid_sel_elt(block, "alpha");
      if (TYPEOF(al_s) != REALSXP || local_b1 < 1 || local_b1 > XLENGTH(al_s)) {
        vmaxset(vm_cb);
        return 0;
      }
      alpha_val = REAL(al_s)[local_b1 - 1];
    }
    if (sm->d_s > 0) {
      SEXP be_s = cdid_sel_elt(block, "beta");
      if (TYPEOF(be_s) != REALSXP || XLENGTH(be_s) != sm->d_s) {
        vmaxset(vm_cb);
        return 0;
      }
      const double *be = REAL(be_s);
      double *Z_sub = (double *) R_alloc((size_t) n_sub * (size_t) sm->d_s, sizeof(double));
      for (int i = 0; i < n_sub; ++i) {
        int r1 = rows1_sub[i];
        int zpos = rf->has_ref ? cdid_sel_ref_s_row_idx(rf, r_id, r1, &hz) : (r1 - 1);
        if (zpos < 0) { vmaxset(vm_cb); return 0; }
        for (int j = 0; j < sm->d_s; ++j) {
          double v = Z_base[zpos + (size_t) (sm->slope_cols[j] - 1) * (size_t) z_ld];
          if (!R_FINITE(v)) { vmaxset(vm_cb); return 0; }
          Z_sub[i + (size_t) j * (size_t) n_sub] = v;
        }
      }
      BLAS_INT m_b = (BLAS_INT) n_sub, n_b = 1, k_b = (BLAS_INT) sm->d_s;
      double one_d = 1.0, zero_d = 0.0;
      F77_CALL(dgemm)("N", "N", &m_b, &n_b, &k_b, &one_d, Z_sub, &m_b, be, &k_b, &zero_d, eta_sub, &m_b FCONE FCONE);
      if (has_int) {
        for (int i = 0; i < n_sub; ++i) eta_sub[i] = alpha_val + eta_sub[i];
      }
    } else {
      for (int i = 0; i < n_sub; ++i) eta_sub[i] = alpha_val;
    }
  }

  for (int i = 0; i < n_sub; ++i) {
    double p = 0.0, od = 0.0, act = 0.0;
    cdid_sel_clip_logit_scalar(eta_sub[i], &p, &od, &act);
    if (!R_FINITE(p) || !(p > 0.0) || p >= 1.0) {
      vmaxset(vm_cb);
      return 0;
    }
    prob_out[i] = p;
  }
  vmaxset(vm_cb);
  return 1;
}

static int cdid_sel_predict_pair_batch(
    const cdid_sel_model_t *sm,
    const cdid_sel_reader_t *rd,
    const cdid_sel_ref_t *rf,
    int q1,
    const int *rows1_sub,
    int n_sub,
    double *prob_out,
    double *ps_out,
    double *pt_out,
    double *w_sel_out,
    double *a_pc_out) {
  if (n_sub < 0 || q1 < 1 || q1 > sm->n_pc) return 0;
  SEXP pc = VECTOR_ELT(sm->pair_cohorts, q1 - 1);
  int coh1 = Rf_asInteger(cdid_sel_elt(pc, "cohort"));
  SEXP ref_r_s = cdid_sel_elt(pc, "reference");
  int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
  double a_pc = 0.0;
  if (!rf->has_ref) {
    a_pc = rd->a_vec[coh1 - 1];
  } else {
    if (!cdid_sel_ref_stat_one(rf, r_id, coh1, 0, &a_pc)) return 0;
  }
  if (a_pc_out) *a_pc_out = a_pc;
  if (n_sub == 0) return 1;

  const void *vm_pb = vmaxget();
  if (!sm->is_date_model) {
    int cid1 = Rf_asInteger(cdid_sel_elt(pc, "pair_cell_id"));
    if (!cdid_sel_predict_cell_batch(sm, rd, rf, cid1, rows1_sub, n_sub, prob_out)) {
      vmaxset(vm_pb);
      return 0;
    }
    if (ps_out) memset(ps_out, 0, (size_t) n_sub * sizeof(double));
    if (pt_out) memset(pt_out, 0, (size_t) n_sub * sizeof(double));
  } else {
    int scid1 = Rf_asInteger(cdid_sel_elt(pc, "s_cell_id"));
    int tcid1 = Rf_asInteger(cdid_sel_elt(pc, "t_cell_id"));
    double *ps_buf = ps_out ? ps_out : (double *) R_alloc((size_t) n_sub, sizeof(double));
    double *pt_buf = pt_out ? pt_out : (double *) R_alloc((size_t) n_sub, sizeof(double));
    if (!cdid_sel_predict_cell_batch(sm, rd, rf, scid1, rows1_sub, n_sub, ps_buf) ||
        !cdid_sel_predict_cell_batch(sm, rd, rf, tcid1, rows1_sub, n_sub, pt_buf)) {
      vmaxset(vm_pb);
      return 0;
    }
    for (int i = 0; i < n_sub; ++i) {
      double pp = ps_buf[i] * pt_buf[i];
      if (!R_FINITE(pp) || !(pp > 0.0) || pp >= 1.0) {
        vmaxset(vm_pb);
        return 0;
      }
      prob_out[i] = pp;
    }
  }

  if (w_sel_out) {
    int hw = -1;
    for (int i = 0; i < n_sub; ++i) {
      int r1 = rows1_sub[i];
      double w_row = 0.0;
      if (!rf->has_ref) {
        w_row = rd->w_static[r1 - 1];
      } else {
        if (!cdid_sel_ref_w_one(rf, r_id, r1, &hw, &w_row)) {
          vmaxset(vm_pb);
          return 0;
        }
      }
      double p_weight = (prob_out[i] > CDID_SEL_CLIP_EPS) ? prob_out[i] : CDID_SEL_CLIP_EPS;
      double w_sel = (w_row * a_pc) / (rd->rho[r1 - 1] * p_weight);
      if (!R_FINITE(w_sel) || !(w_sel > 0.0)) {
        vmaxset(vm_pb);
        return 0;
      }
      w_sel_out[i] = w_sel;
    }
  }
  vmaxset(vm_pb);
  return 1;
}

static int cdid_sel_eval_upstream_internal(
    const cdid_sel_model_t *sm,
    SEXP up_map_sexp,
    const cdid_sel_reader_t *rd,
    const cdid_sel_ref_t *rf,
    const int *rows1,
    int B,
    const int *uids,
    int nu,
    SEXP prep_idx_sexp,
    double *up_out,
    int *n_act_cells_out,
    int *n_act_blocks_out) {
  memset(up_out, 0, (size_t) B * (size_t) nu * sizeof(double));
  if (TYPEOF(up_map_sexp) != VECSXP) return 0;
  int n_a = Rf_asInteger(cdid_sel_elt(up_map_sexp, "n_a"));
  SEXP coh_a_s = cdid_sel_elt(up_map_sexp, "cohorts_a");
  SEXP rref_s = cdid_sel_elt(up_map_sexp, "rates_ref");
  SEXP tc_list = cdid_sel_elt(up_map_sexp, "theta_coords");
  if (TYPEOF(coh_a_s) != INTSXP || TYPEOF(tc_list) != VECSXP) return 0;
  const int *cohorts_a = INTEGER(coh_a_s);
  const int *rates_ref = (TYPEOF(rref_s) == INTSXP) ? INTEGER(rref_s) : NULL;

  /* 1. Rate columns */
  for (int j = 0; j < nu; ++j) {
    int uid = uids[j];
    if (uid < 1 || uid > n_a) continue;
    int k = cohorts_a[uid - 1];
    int r_id = rates_ref ? rates_ref[uid - 1] : 0;
    double a_k = 0.0, m_k = 0.0;
    int have_stats = 0;
    int hw = -1;
    double *col_j = up_out + (size_t) j * (size_t) B;
    for (int i = 0; i < B; ++i) {
      int r1 = rows1[i];
      if (rd->cohort_index[r1 - 1] == k && rd->U_logical[r1 - 1] == 1) {
        if (!have_stats) {
          if (!rf->has_ref) {
            a_k = rd->a_vec[k - 1];
            m_k = rd->mass_U[k - 1];
          } else {
            if (!cdid_sel_ref_stat_one(rf, r_id, k, 0, &a_k) ||
                !cdid_sel_ref_stat_one(rf, r_id, k, 1, &m_k)) return 0;
          }
          have_stats = 1;
        }
        double wi = 0.0;
        if (!rf->has_ref) wi = rd->w_static[r1 - 1];
        else if (!cdid_sel_ref_w_one(rf, r_id, r1, &hw, &wi)) return 0;
        double Ai = (rd->a_index[r1 - 1] > 0) ? 1.0 : 0.0;
        col_j[i] = wi * (Ai - a_k) / m_k;
      }
    }
  }

  /* 2. Determine active blocks and U rows */
  int n_U = 0;
  int *u_rows_local = (int *) R_alloc((size_t) B, sizeof(int)); /* local chunk index 0..B-1 of U rows */
  int *u_pos_of_i = (int *) R_alloc((size_t) B, sizeof(int));   /* -1 if !U, else 0..n_U-1 */
  for (int i = 0; i < B; ++i) {
    if (rd->U_logical[rows1[i] - 1] == 1) {
      u_pos_of_i[i] = n_U;
      u_rows_local[n_U] = i;
      ++n_U;
    } else {
      u_pos_of_i[i] = -1;
    }
  }

  int D_max = Rf_asInteger(cdid_sel_elt(up_map_sexp, "D"));
  int *uid_to_col = NULL;
  if (prep_idx_sexp == R_NilValue && D_max > 0) {
    uid_to_col = (int *) R_alloc((size_t) (D_max + 1), sizeof(int));
    memset(uid_to_col, 0, (size_t) (D_max + 1) * sizeof(int));
    for (int j = 0; j < nu; ++j) {
      if (uids[j] >= 1 && uids[j] <= D_max) uid_to_col[uids[j]] = j + 1; /* 1-based */
    }
  }

  double *Z_slope_static = NULL;
  if (!rf->has_ref && sm->d_s > 0 && n_U > 0) {
    Z_slope_static = (double *) R_alloc((size_t) n_U * (size_t) sm->d_s, sizeof(double));
    for (int j = 0; j < sm->d_s; ++j) {
      const double *zj = rd->Z_static + (size_t) (sm->slope_cols[j] - 1) * (size_t) rd->N;
      double *dst = Z_slope_static + (size_t) j * (size_t) n_U;
      for (int u = 0; u < n_U; ++u) {
        int r1 = rows1[u_rows_local[u]];
        dst[u] = zj[r1 - 1];
      }
    }
  }

  int *cell_seen = (sm->n_cells > 0) ? (int *) R_alloc((size_t) (sm->n_cells + 1), sizeof(int)) : NULL;
  if (cell_seen) memset(cell_seen, 0, (size_t) (sm->n_cells + 1) * sizeof(int));
  int n_act_cells = 0, n_act_blocks = 0;

  SEXP prep_act_blocks = (prep_idx_sexp != R_NilValue) ? cdid_sel_elt(prep_idx_sexp, "active_blocks") : R_NilValue;
  SEXP prep_keep_list = (prep_idx_sexp != R_NilValue) ? cdid_sel_elt(prep_idx_sexp, "keep") : R_NilValue;
  SEXP prep_dest_list = (prep_idx_sexp != R_NilValue) ? cdid_sel_elt(prep_idx_sexp, "destination") : R_NilValue;
  SEXP prep_act_cells = (prep_idx_sexp != R_NilValue) ? cdid_sel_elt(prep_idx_sexp, "active_cells") : R_NilValue;

  int n_loop = (prep_act_blocks != R_NilValue) ? (int) XLENGTH(prep_act_blocks) : sm->n_blocks;

  for (int idx_m = 0; idx_m < n_loop; ++idx_m) {
    int m1 = (prep_act_blocks != R_NilValue) ? INTEGER(prep_act_blocks)[idx_m] : (idx_m + 1);
    if (m1 < 1 || m1 > sm->n_blocks) return 0;
    SEXP tc_s = VECTOR_ELT(tc_list, m1 - 1);
    int n_tc = (int) XLENGTH(tc_s);
    const int *tc = (n_tc > 0) ? INTEGER(tc_s) : NULL;

    int n_keep = 0;
    const int *keep1 = NULL;
    const int *dest1 = NULL;
    int *local_keep = NULL, *local_dest = NULL;

    if (prep_idx_sexp != R_NilValue) {
      SEXP ks = VECTOR_ELT(prep_keep_list, m1 - 1);
      SEXP ds = VECTOR_ELT(prep_dest_list, m1 - 1);
      n_keep = (ks != R_NilValue) ? (int) XLENGTH(ks) : 0;
      if (n_keep > 0) {
        keep1 = INTEGER(ks);
        dest1 = INTEGER(ds);
      }
    } else {
      if (n_tc == 0) continue;
      local_keep = (int *) R_alloc((size_t) n_tc, sizeof(int));
      local_dest = (int *) R_alloc((size_t) n_tc, sizeof(int));
      for (int t = 0; t < n_tc; ++t) {
        int u_coord = tc[t];
        if (u_coord >= 1 && u_coord <= D_max && uid_to_col[u_coord] > 0) {
          local_keep[n_keep] = t + 1;
          local_dest[n_keep] = uid_to_col[u_coord];
          ++n_keep;
        }
      }
      if (n_keep == 0) continue;
      keep1 = local_keep;
      dest1 = local_dest;
    }

    ++n_act_blocks;
    SEXP bl = VECTOR_ELT(sm->blocks, m1 - 1);
    SEXP cids_s = cdid_sel_elt(bl, "cell_ids");
    int B_m = Rf_asInteger(cdid_sel_elt(bl, "B_m"));
    const int *cids = INTEGER(cids_s);
    for (int b = 0; b < B_m; ++b) {
      int cid = cids[b];
      if (cell_seen && cid >= 1 && cid <= sm->n_cells && !cell_seen[cid]) {
        cell_seen[cid] = 1;
        ++n_act_cells;
      }
    }

    if (n_U == 0 || n_keep == 0) continue;
    SEXP br_s = cdid_sel_elt(bl, "branch");
    const char *br = CHAR(STRING_ELT(br_s, 0));
    if (strcmp(br, "zero_col") == 0) continue;

    double *cell_scores_U = (double *) R_alloc((size_t) n_U * (size_t) B_m, sizeof(double));
    memset(cell_scores_U, 0, (size_t) n_U * (size_t) B_m * sizeof(double));

    double *Zb = Z_slope_static;
    if (rf->has_ref && sm->d_s > 0) {
      Zb = (double *) R_alloc((size_t) n_U * (size_t) sm->d_s, sizeof(double));
      memset(Zb, 0, (size_t) n_U * (size_t) sm->d_s * sizeof(double));
      SEXP bspec = VECTOR_ELT(sm->block_specs, m1 - 1);
      int rb = Rf_asInteger(cdid_sel_elt(bspec, "reference"));
      int *coh_in_blk = (int *) R_alloc((size_t) (rd->J + 1), sizeof(int));
      memset(coh_in_blk, 0, (size_t) (rd->J + 1) * sizeof(int));
      for (int b = 0; b < B_m; ++b) {
        int ck = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(sm->cells, cids[b] - 1), "cohort"));
        if (ck >= 1 && ck <= rd->J) coh_in_blk[ck] = 1;
      }
      int hz = -1;
      for (int u = 0; u < n_U; ++u) {
        int r1 = rows1[u_rows_local[u]];
        int ck = rd->cohort_index[r1 - 1];
        if (ck >= 1 && ck <= rd->J && coh_in_blk[ck]) {
          int zpos = cdid_sel_ref_s_row_idx(rf, rb, r1, &hz);
          if (zpos < 0) return 0;
          for (int j = 0; j < sm->d_s; ++j) {
            Zb[u + (size_t) j * (size_t) n_U] = rf->s_Z[zpos + (size_t) (sm->slope_cols[j] - 1) * (size_t) rf->n_s_rows];
          }
        }
      }
    }

    int *ck_u = (int *) R_alloc((size_t) n_U, sizeof(int));
    int *ck_rows1 = (int *) R_alloc((size_t) n_U, sizeof(int));
    double *ck_pr = (double *) R_alloc((size_t) n_U, sizeof(double));
    for (int b = 0; b < B_m; ++b) {
      int cid = cids[b];
      SEXP cell = VECTOR_ELT(sm->cells, cid - 1);
      int ck = Rf_asInteger(cdid_sel_elt(cell, "cohort"));
      SEXP ref_r_s = cdid_sel_elt(cell, "reference");
      int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
      int ds1 = 0, dt1 = 0;
      if (sm->is_date_model) {
        dt1 = Rf_asInteger(cdid_sel_elt(cdid_sel_elt(cell, "date_spec"), "calendar_index"));
      } else {
        SEXP pair = cdid_sel_elt(cell, "pair");
        ds1 = Rf_asInteger(cdid_sel_elt(pair, "s_calendar_index"));
        dt1 = Rf_asInteger(cdid_sel_elt(pair, "t_calendar_index"));
      }
      int n_ck = 0;
      for (int u = 0; u < n_U; ++u) {
        int r1 = rows1[u_rows_local[u]];
        if (rd->cohort_index[r1 - 1] == ck) {
          ck_u[n_ck] = u;
          ck_rows1[n_ck] = r1;
          ++n_ck;
        }
      }
      if (n_ck == 0) continue;
      if (!cdid_sel_predict_cell_batch(sm, rd, rf, cid, ck_rows1, n_ck, ck_pr)) return 0;
      int hs = -1, ht = -1, hu = -1, hw = -1;
      double *cs_b = cell_scores_U + (size_t) b * (size_t) n_U;
      for (int idx = 0; idx < n_ck; ++idx) {
        int u = ck_u[idx];
        int r1 = ck_rows1[idx];
        double Y = 0.0, dummy = 0.0;
        if (sm->is_date_model) {
          Y = cdid_sel_date_one(rd, r1, dt1, &hs, &hu);
        } else {
          if (!cdid_sel_pair_one(rd, r1, ds1, dt1, 0, &hs, &ht, &Y, &dummy)) return 0;
        }
        double wi = 0.0;
        if (!rf->has_ref) wi = rd->w_static[r1 - 1];
        else if (!cdid_sel_ref_w_one(rf, r_id, r1, &hw, &wi)) return 0;
        cs_b[u] = wi * (Y - ck_pr[idx]);
      }
    }

    int has_int = Rf_asLogical(cdid_sel_elt(bl, "has_intercept")) == 1;
    int d_s = sm->d_s;
    int p_theta = (has_int ? B_m : 0) + d_s;
    double *xi_theta_U = (double *) R_alloc((size_t) n_U * (size_t) p_theta, sizeof(double));
    double *tmp_p = (double *) R_alloc((size_t) (p_theta > 0 ? p_theta : 1), sizeof(double));

    if (strcmp(br, "cp04_separate_pair") == 0) {
      SEXP H_s = cdid_sel_elt(bl, "H_chol");
      SEXP coord_s = cdid_sel_elt(bl, "coordinate_order");
      if (TYPEOF(H_s) != REALSXP || TYPEOF(coord_s) != INTSXP || XLENGTH(coord_s) != p_theta) return 0;
      const int *coord = INTEGER(coord_s);
      double *score_raw = (double *) R_alloc((size_t) n_U * (size_t) p_theta, sizeof(double));
      double *xi_raw = (double *) R_alloc((size_t) n_U * (size_t) p_theta, sizeof(double));
      for (int c = 0; c < p_theta; ++c) {
        int raw_col = coord[c] - 1;
        if (raw_col < 0 || raw_col >= p_theta) return 0;
        double *sr = score_raw + (size_t) raw_col * (size_t) n_U;
        if (has_int && c == 0) {
          memcpy(sr, cell_scores_U, (size_t) n_U * sizeof(double));
        } else {
          int sj = has_int ? (c - 1) : c;
          const double *zcol = Zb + (size_t) sj * (size_t) n_U;
          for (int u = 0; u < n_U; ++u) sr[u] = cell_scores_U[u] * zcol[u];
        }
      }
      if (!cdid_sel_rtr_solve_rows(REAL(H_s), NULL, score_raw, xi_raw, n_U, p_theta, 1.0, tmp_p)) return 0;
      for (int c = 0; c < p_theta; ++c) {
        int raw_col = coord[c] - 1;
        memcpy(xi_theta_U + (size_t) c * (size_t) n_U,
               xi_raw + (size_t) raw_col * (size_t) n_U,
               (size_t) n_U * sizeof(double));
      }
    } else {
      const double *mb = REAL(cdid_sel_elt(bl, "m_b"));
      double *xi_beta_U = xi_theta_U + (size_t) (has_int ? B_m : 0) * (size_t) n_U;
      BLAS_INT nU_b = (BLAS_INT) n_U, Bm_b = (BLAS_INT) B_m, ds_b = (BLAS_INT) d_s;
      double one_d = 1.0, zero_d = 0.0;
      if (d_s > 0) {
        SEXP sf_s = cdid_sel_elt(bl, "slope_fit");
        SEXP R_s = cdid_sel_elt(sf_s, "R");
        SEXP piv_s = cdid_sel_elt(sf_s, "pivot");
        if (TYPEOF(R_s) != REALSXP || TYPEOF(piv_s) != INTSXP) return 0;
        double *s_tilde = (double *) R_alloc((size_t) n_U * (size_t) d_s, sizeof(double));
        for (int u = 0; u < n_U; ++u) {
          double sum_s = 0.0;
          for (int b = 0; b < B_m; ++b) sum_s += cell_scores_U[u + (size_t) b * (size_t) n_U];
          for (int j = 0; j < d_s; ++j) {
            s_tilde[u + (size_t) j * (size_t) n_U] = sum_s * Zb[u + (size_t) j * (size_t) n_U];
          }
        }
        if (has_int && B_m > 0) {
          const double *zbar = REAL(cdid_sel_elt(bl, "z_bar"));
          double *cs_zb = (double *) R_alloc((size_t) n_U * (size_t) d_s, sizeof(double));
          F77_CALL(dgemm)("N", "N", &nU_b, &ds_b, &Bm_b, &one_d,
                          cell_scores_U, &nU_b, zbar, &Bm_b, &zero_d, cs_zb, &nU_b FCONE FCONE);
          size_t tot_st = (size_t) n_U * (size_t) d_s;
          for (size_t idx = 0; idx < tot_st; ++idx) s_tilde[idx] -= cs_zb[idx];
        }
        if (!cdid_sel_rtr_solve_rows(REAL(R_s), INTEGER(piv_s), s_tilde, xi_beta_U, n_U, d_s, rd->N_dbl, tmp_p)) return 0;
      }
      if (has_int) {
        const double *zbar = (d_s > 0) ? REAL(cdid_sel_elt(bl, "z_bar")) : NULL;
        for (int b = 0; b < B_m; ++b) {
          double mb_b = mb[b];
          const double *cs_b = cell_scores_U + (size_t) b * (size_t) n_U;
          double *xa_b = xi_theta_U + (size_t) b * (size_t) n_U;
          for (int u = 0; u < n_U; ++u) xa_b[u] = cs_b[u] / mb_b;
        }
        if (d_s > 0 && B_m > 0) {
          double *xb_zb = (double *) R_alloc((size_t) n_U * (size_t) B_m, sizeof(double));
          F77_CALL(dgemm)("N", "T", &nU_b, &Bm_b, &ds_b, &one_d,
                          xi_beta_U, &nU_b, zbar, &Bm_b, &zero_d, xb_zb, &nU_b FCONE FCONE);
          size_t tot_xa = (size_t) n_U * (size_t) B_m;
          for (size_t idx = 0; idx < tot_xa; ++idx) xi_theta_U[idx] -= xb_zb[idx];
        }
      }
    }

    for (int t = 0; t < n_keep; ++t) {
      int src_col = keep1[t] - 1;
      int dst_col = dest1[t] - 1;
      if (src_col < 0 || src_col >= p_theta || dst_col < 0 || dst_col >= nu) return 0;
      const double *src = xi_theta_U + (size_t) src_col * (size_t) n_U;
      double *dst = up_out + (size_t) dst_col * (size_t) B;
      for (int u = 0; u < n_U; ++u) dst[u_rows_local[u]] = src[u];
    }
  }

  if (prep_act_cells != R_NilValue) n_act_cells = (int) XLENGTH(prep_act_cells);
  if (n_act_cells_out) *n_act_cells_out = n_act_cells;
  if (n_act_blocks_out) *n_act_blocks_out = n_act_blocks;
  return 1;
}

/* ---------------------------------------------------------------------------
   5. C_cdid_sel_upstream_chunk(sel_fit_sexp, up_map_sexp, reader_sexp, ref_sexp, rows_sexp, uids_sexp, prep_idx_sexp)
   --------------------------------------------------------------------------- */
SEXP C_cdid_sel_upstream_chunk(
    SEXP sel_fit_sexp,
    SEXP up_map_sexp,
    SEXP reader_sexp,
    SEXP ref_sexp,
    SEXP rows_sexp,
    SEXP uids_sexp,
    SEXP prep_idx_sexp) {
  cdid_sel_model_t sm;
  if (!cdid_sel_init_model(sel_fit_sexp, &sm)) return R_NilValue;
  cdid_sel_reader_t rd;
  if (!cdid_sel_init_reader(reader_sexp, &rd, sm.is_date_model)) return R_NilValue;
  cdid_sel_ref_t rf;
  if (!cdid_sel_init_ref(ref_sexp, &rf)) return R_NilValue;

  if (Rf_isObject(rows_sexp) || (TYPEOF(rows_sexp) != INTSXP && TYPEOF(rows_sexp) != REALSXP) ||
      (TYPEOF(uids_sexp) != INTSXP && TYPEOF(uids_sexp) != REALSXP)) return R_NilValue;
  int B = (int) XLENGTH(rows_sexp);
  int nu = (int) XLENGTH(uids_sexp);
  if (B < 1) return R_NilValue;

  const void *vmax = vmaxget();
  int *rows1 = (int *) R_alloc((size_t) B, sizeof(int));
  for (int i = 0; i < B; ++i) {
    int r1 = cdid_sel_int_at(rows_sexp, i);
    if (r1 == NA_INTEGER || r1 < 1 || r1 > rd.N) { vmaxset(vmax); return R_NilValue; }
    rows1[i] = r1;
  }
  int *uids = (nu > 0) ? (int *) R_alloc((size_t) nu, sizeof(int)) : NULL;
  for (int j = 0; j < nu; ++j) {
    int u = cdid_sel_int_at(uids_sexp, j);
    if (u == NA_INTEGER) { vmaxset(vmax); return R_NilValue; }
    uids[j] = u;
  }

  SEXP up_s = PROTECT(Rf_allocMatrix(REALSXP, B, nu));
  int n_act_cells = 0, n_act_blocks = 0;
  if (!cdid_sel_eval_upstream_internal(&sm, up_map_sexp, &rd, &rf, rows1, B, uids, nu,
                                       prep_idx_sexp, REAL(up_s), &n_act_cells, &n_act_blocks)) {
    UNPROTECT(1);
    vmaxset(vmax);
    return R_NilValue;
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_s, 0, up_s);
  SET_VECTOR_ELT(out_s, 1, Rf_ScalarInteger(n_act_cells));
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarInteger(n_act_blocks));
  UNPROTECT(2);
  vmaxset(vmax);
  return out_s;
}

/* ---------------------------------------------------------------------------
   Internal Helpers for Treatment Fits, Support Collection & IF Chunk Kernels
   --------------------------------------------------------------------------- */

/* Build O(1) lookup maps obs_s[1..n_A], obs_t[1..n_A] (1-based obs_index, 0 if not observed at BOTH dates)
   via O(count_s + count_t) two-pointer merge on sorted post_a slices. */
static int cdid_sel_build_pair_map(
    const cdid_sel_reader_t *rd,
    int ds1,
    int dt1,
    int *obs_s_out,
    int *obs_t_out) {
  memset(obs_s_out, 0, (size_t) (rd->n_A + 1) * sizeof(int));
  if (obs_t_out) memset(obs_t_out, 0, (size_t) (rd->n_A + 1) * sizeof(int));
  if (ds1 < 1 || ds1 > rd->n_dates || dt1 < 1 || dt1 > rd->n_dates) return 0;
  if (rd->date_count[ds1 - 1] <= 0 || rd->date_count[dt1 - 1] <= 0) return 1;

  int is0 = rd->date_start[ds1 - 1] - 1;
  int is1 = rd->date_end[ds1 - 1] - 1;
  int it0 = rd->date_start[dt1 - 1] - 1;
  int it1 = rd->date_end[dt1 - 1] - 1;
  if (is0 < 0 || is1 >= rd->n_post || it0 < 0 || it1 >= rd->n_post) return 0;

  int ps = is0, pt = it0;
  while (ps <= is1 && pt <= it1) {
    int as = rd->post_a[ps];
    int at = rd->post_a[pt];
    if (as == at) {
      if (as < 1 || as > rd->n_A) return 0;
      int os = rd->post_obs[ps];
      int ot = rd->post_obs[pt];
      if (os < 1 || (R_xlen_t) os > rd->n_y || ot < 1 || (R_xlen_t) ot > rd->n_y) return 0;
      obs_s_out[as] = os;
      if (obs_t_out) obs_t_out[as] = ot;
      ++ps;
      ++pt;
    } else if (as < at) {
      ++ps;
    } else {
      ++pt;
    }
  }
  return 1;
}

/* Predict legacy CP04 selection probability for row1 (1-based) */
static inline int cdid_sel_predict_legacy_one(
    const cdid_sel_reader_t *rd,
    const cdid_sel_ref_t *rf,
    int r_id,
    const double *coef,
    int d_cols,
    int row1,
    int *hint_z,
    double *prob_out) {
  int zpos = rf->has_ref ? cdid_sel_ref_s_row_idx(rf, r_id, row1, hint_z) : (row1 - 1);
  if (zpos < 0) return 0;
  const double *Z_base = rf->has_ref ? rf->s_Z : rd->Z_static;
  int z_ld = rf->has_ref ? rf->n_s_rows : rd->N;
  double eta = 0.0;
  for (int j = 0; j < d_cols; ++j) {
    eta += Z_base[zpos + (size_t) j * (size_t) z_ld] * coef[j];
  }
  double p = 0.0, od = 0.0, act = 0.0;
  cdid_sel_clip_logit_scalar(eta, &p, &od, &act);
  if (!R_FINITE(p) || !(p > 0.0) || p >= 1.0) return 0;
  *prob_out = p;
  return 1;
}

/* Compute T_logit row vector (length dq_expected) for pair_cohort q1 on row1 given precomputed (prob, ps, pt) */
static int cdid_sel_deriv_tlogit_one(
    const cdid_sel_model_t *sm,
    const cdid_sel_reader_t *rd,
    const cdid_sel_ref_t *rf,
    int q1,
    int row1,
    double prob,
    double ps,
    double pt,
    int *hint_z,
    double *tlogit_out,
    int dq_expected) {
  if (dq_expected == 0) return 1;
  if (q1 < 1 || q1 > sm->n_pc) return 0;
  SEXP pc = VECTOR_ELT(sm->pair_cohorts, q1 - 1);
  SEXP ref_r_s = cdid_sel_elt(pc, "reference");
  int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

  double act = (prob > CDID_SEL_CLIP_EPS) ? 1.0 : 0.0;
  int d_s = sm->d_s;
  const double *Z_base = rf->has_ref ? rf->s_Z : rd->Z_static;
  int z_ld = rf->has_ref ? rf->n_s_rows : rd->N;
  int zpos = 0;
  if (d_s > 0) {
    zpos = rf->has_ref ? cdid_sel_ref_s_row_idx(rf, r_id, row1, hint_z) : (row1 - 1);
    if (zpos < 0) return 0;
  }

  int pos = 0;
  if (!sm->is_date_model) {
    int m1 = Rf_asInteger(cdid_sel_elt(pc, "block_id"));
    SEXP bl = VECTOR_ELT(sm->blocks, m1 - 1);
    const char *br = CHAR(STRING_ELT(cdid_sel_elt(bl, "branch"), 0));
    if (strcmp(br, "zero_col") == 0) return (dq_expected == 0) ? 1 : 0;
    int has_int = Rf_asLogical(cdid_sel_elt(bl, "has_intercept")) == 1;
    double w_p = -((1.0 - prob) * act);
    if (has_int) tlogit_out[pos++] = w_p;
    for (int j = 0; j < d_s; ++j) {
      tlogit_out[pos++] = w_p * Z_base[zpos + (size_t) (sm->slope_cols[j] - 1) * (size_t) z_ld];
    }
  } else {
    int ms1 = Rf_asInteger(cdid_sel_elt(pc, "s_block_id"));
    int mt1 = Rf_asInteger(cdid_sel_elt(pc, "t_block_id"));
    SEXP bl_s = VECTOR_ELT(sm->blocks, ms1 - 1);
    SEXP bl_t = VECTOR_ELT(sm->blocks, mt1 - 1);
    const char *br_s = CHAR(STRING_ELT(cdid_sel_elt(bl_s, "branch"), 0));
    if (strcmp(br_s, "zero_col") == 0) return (dq_expected == 0) ? 1 : 0;
    double w_s = -((1.0 - ps) * act);
    double w_t = -((1.0 - pt) * act);
    int hi_s = Rf_asLogical(cdid_sel_elt(bl_s, "has_intercept")) == 1;
    int hi_t = Rf_asLogical(cdid_sel_elt(bl_t, "has_intercept")) == 1;
    if (ms1 != mt1) {
      if (hi_s) tlogit_out[pos++] = w_s;
      for (int j = 0; j < d_s; ++j) {
        tlogit_out[pos++] = w_s * Z_base[zpos + (size_t) (sm->slope_cols[j] - 1) * (size_t) z_ld];
      }
      if (hi_t) tlogit_out[pos++] = w_t;
      for (int j = 0; j < d_s; ++j) {
        tlogit_out[pos++] = w_t * Z_base[zpos + (size_t) (sm->slope_cols[j] - 1) * (size_t) z_ld];
      }
    } else {
      if (hi_s) {
        tlogit_out[pos++] = w_s;
        tlogit_out[pos++] = w_t;
      }
      double w_st = w_s + w_t;
      for (int j = 0; j < d_s; ++j) {
        tlogit_out[pos++] = w_st * Z_base[zpos + (size_t) (sm->slope_cols[j] - 1) * (size_t) z_ld];
      }
    }
  }
  return (pos == dq_expected) ? 1 : 0;
}

/* Internal IRLS solver matching C_cdid_cp04_logit_fast exactly */
static int cdid_sel_logit_fit_core(
    const double *X,
    const double *y,
    const double *w,
    int n,
    int p,
    double eps,
    int maxit,
    double N_val,
    double *beta_out,
    double *R_chol_out,
    int *iter_out) {
  if (p < 1 || p > 128 || n < p) return 0;
  int n_pos_w = 0;
  long double sum_w_ld = 0.0L;
  long double devold_ld = 0.0L;
  const double ln2_val = log(2.0);
  for (int i = 0; i < n; ++i) {
    double wi = w[i];
    if (!R_FINITE(wi) || !(wi > 0.0)) return 0;
    ++n_pos_w;
    sum_w_ld += (long double) wi;
    devold_ld += (long double) (2.0 * wi * ln2_val);
  }
  double sum_w = (double) sum_w_ld;
  if (n_pos_w < p || !R_FINITE(sum_w) || sum_w <= 0.0) return 0;

  for (int j = 0; j < p; ++j) {
    const double *xj = X + (size_t) j * (size_t) n;
    int has_pos = 0, has_neg = 0;
    for (int i = 0; i < n; ++i) {
      double v = xj[i];
      if (!R_FINITE(v)) return 0;
      if (v != 0.0) {
        double sgn = (y[i] == 1.0) ? v : -v;
        if (sgn > 0.0) has_pos = 1;
        else if (sgn < 0.0) has_neg = 1;
      }
    }
    if ((has_pos && !has_neg) || (has_neg && !has_pos)) return 0;
  }

  double *eta = (double *) R_alloc((size_t) n, sizeof(double));
  double *mu = (double *) R_alloc((size_t) n, sizeof(double));
  double *wt = (double *) R_alloc((size_t) n, sizeof(double));
  double *wz = (double *) R_alloc((size_t) n, sizeof(double));
  double *W_X = (double *) R_alloc((size_t) n * (size_t) p, sizeof(double));
  double *H = (double *) R_alloc((size_t) p * (size_t) p, sizeof(double));

  for (int i = 0; i < n; ++i) { eta[i] = 0.0; mu[i] = 0.5; }
  memset(beta_out, 0, (size_t) p * sizeof(double));

  double devold = (double) devold_ld;
  int conv = 0, iter_used = 0;
  const double bnd_eps = 10.0 * DBL_EPSILON;
  const double qr_tol = (eps / 1000.0 < 1e-7) ? (eps / 1000.0) : 1e-7;
  BLAS_INT n_b = (BLAS_INT) n, p_b = (BLAS_INT) p, one_b = 1;
  double one_d = 1.0, zero_d = 0.0;

  for (int iter = 1; iter <= maxit; ++iter) {
    iter_used = iter;
    for (int i = 0; i < n; ++i) {
      double mui = mu[i];
      double varmu = mui * (1.0 - mui);
      double opexp = exp(eta[i]);
      double mu_eta = opexp / ((1.0 + opexp) * (1.0 + opexp));
      if (!(varmu > 0.0) || !(mu_eta > 0.0) || !R_FINITE(mu_eta)) return 0;
      double w_i = sqrt((w[i] * (mu_eta * mu_eta)) / varmu);
      double z_i = eta[i] + (y[i] - mui) / mu_eta;
      wt[i] = w_i;
      wz[i] = z_i * w_i;
    }
    for (int j = 0; j < p; ++j) {
      const double *xj = X + (size_t) j * (size_t) n;
      double *wxj = W_X + (size_t) j * (size_t) n;
      for (int i = 0; i < n; ++i) wxj[i] = xj[i] * wt[i];
    }
    if (!cdid_sel_qr_wls_solve(W_X, wz, R_chol_out, beta_out, n, p, qr_tol)) return 0;

    F77_CALL(dgemm)("N", "N", &n_b, &one_b, &p_b, &one_d,
                    X, &n_b, beta_out, &p_b, &zero_d, eta, &n_b FCONE FCONE);

    long double dev_ld = 0.0L;
    for (int i = 0; i < n; ++i) {
      double ei = eta[i];
      if (!R_FINITE(ei) || fabs(ei) >= 30.0) return 0;
      double opexp = exp(ei);
      double mui = opexp / (1.0 + opexp);
      if (mui <= bnd_eps || mui >= 1.0 - bnd_eps) return 0;
      mu[i] = mui;
      double loss = (y[i] == 1.0) ? log(1.0 / mui) : log(1.0 / (1.0 - mui));
      dev_ld += (long double) (2.0 * w[i] * loss);
    }
    double dev = (double) dev_ld;
    if (!R_FINITE(dev)) return 0;
    if (fabs(dev - devold) / (0.1 + fabs(dev)) < eps) { conv = 1; break; }
    devold = dev;
  }
  if (!conv) return 0;

  for (int i = 0; i < n; ++i) {
    double od = exp(eta[i]);
    double pr = od / (1.0 + od);
    if (pr <= bnd_eps || pr >= 1.0 - bnd_eps || !(pr * (1.0 - pr) > 0.0)) return 0;
    wt[i] = w[i] * pr * (1.0 - pr);
  }
  for (int j = 0; j < p; ++j) {
    const double *xj = X + (size_t) j * (size_t) n;
    double *wxj = W_X + (size_t) j * (size_t) n;
    for (int i = 0; i < n; ++i) wxj[i] = wt[i] * xj[i];
  }
  F77_CALL(dgemm)("T", "N", &p_b, &p_b, &n_b, &one_d,
                  X, &n_b, W_X, &n_b, &zero_d, H, &p_b FCONE FCONE);
  for (int idx = 0; idx < p * p; ++idx) {
    H[idx] /= N_val;
    if (!R_FINITE(H[idx])) return 0;
  }
  if (!cdid_sel_chol_upper(H, R_chol_out, p, 0.0)) return 0;
  *iter_out = iter_used;
  return 1;
}

/* Build ipw_clip_summary named list:
   if has_prob: (prob, q, controls=(y==0), weights) -> 13 elements
   else:        (q, weights)                        -> 7 elements */
static SEXP cdid_sel_make_clip_summary(
    const double *prob,
    const double *q,
    const double *y,
    const double *weights,
    int n,
    int has_prob) {
  int len = has_prob ? 13 : 7;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, len));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, len));
  int idx = 0;
  #define SET_CLIP_ELT(nm_str, val_dbl) do { \
    SET_STRING_ELT(nms, idx, Rf_mkChar(nm_str)); \
    SET_VECTOR_ELT(out, idx, Rf_ScalarReal(val_dbl)); \
    ++idx; \
  } while (0)

  SET_CLIP_ELT("eps", CDID_SEL_CLIP_EPS);
  if (has_prob) {
    double p_min = (n > 0) ? prob[0] : R_NaReal;
    double p_max = (n > 0) ? prob[0] : R_NaReal;
    double p_capped = 0.0, p_capped_ctrl = 0.0;
    for (int i = 0; i < n; ++i) {
      double pv = prob[i];
      if (pv < p_min) p_min = pv;
      if (pv > p_max) p_max = pv;
      if (!(pv < CDID_SEL_CLIP_CAP)) {
        p_capped += 1.0;
        if (y[i] == 0.0) p_capped_ctrl += 1.0;
      }
    }
    double p_clip_max = (n > 0) ? ((p_max < CDID_SEL_CLIP_CAP) ? p_max : CDID_SEL_CLIP_CAP) : R_NaReal;
    SET_CLIP_ELT("p_min", p_min);
    SET_CLIP_ELT("p_max", p_max);
    SET_CLIP_ELT("p_clip_max", p_clip_max);
    SET_CLIP_ELT("p_n", (double) n);
    SET_CLIP_ELT("p_capped", p_capped);
    SET_CLIP_ELT("p_capped_ctrl", p_capped_ctrl);
  }
  {
    double q_min = (n > 0) ? q[0] : R_NaReal;
    double q_max = (n > 0) ? q[0] : R_NaReal;
    double q_floored = 0.0;
    for (int i = 0; i < n; ++i) {
      double qv = q[i];
      if (qv < q_min) q_min = qv;
      if (qv > q_max) q_max = qv;
      if (qv <= CDID_SEL_CLIP_EPS) q_floored += 1.0;
    }
    double q_floor_min = (n > 0) ? ((q_min > CDID_SEL_CLIP_EPS) ? q_min : CDID_SEL_CLIP_EPS) : R_NaReal;
    SET_CLIP_ELT("q_min", q_min);
    SET_CLIP_ELT("q_max", q_max);
    SET_CLIP_ELT("q_floor_min", q_floor_min);
    SET_CLIP_ELT("q_n", (double) n);
    SET_CLIP_ELT("q_floored", q_floored);
  }
  {
    double w_max = (n > 0) ? weights[0] : R_NaReal;
    for (int i = 1; i < n; ++i) if (weights[i] > w_max) w_max = weights[i];
    SET_CLIP_ELT("w_max", w_max);
  }
  #undef SET_CLIP_ELT
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}

/* ---------------------------------------------------------------------------
   6. C_cdid_sel_trt_fit(mode_sexp, reader_sexp, spec_sexp, sel_sexp, up_map_sexp,
                         ref_sexp, or_info_sexp, cols_sexp, ctrl_sexp)
   modes:
     1 = cp04_treatment_fit (legacy CP04 IPW)
     2 = cp04_selection_treatment_fit (general selection IPW, static & D60)
     3 = dr_selection_treatment_fit (general selection DR, static & D60)
     4 = cp04_custom_callback_fit pre-collector (returns rows, kk, fw, LD, y0, y1)
   --------------------------------------------------------------------------- */
SEXP C_cdid_sel_trt_fit(
    SEXP mode_sexp,
    SEXP reader_sexp,
    SEXP spec_sexp,
    SEXP sel_sexp,
    SEXP up_map_sexp,
    SEXP ref_sexp,
    SEXP or_info_sexp,
    SEXP cols_sexp,
    SEXP ctrl_sexp) {
  int mode = Rf_asInteger(mode_sexp);
  if (mode < 1 || mode > 4 || TYPEOF(spec_sexp) != VECSXP) return R_NilValue;

  cdid_sel_model_t sm;
  memset(&sm, 0, sizeof(sm));
  if (mode >= 2) {
    if (!cdid_sel_init_model(sel_sexp, &sm)) return R_NilValue;
  }
  cdid_sel_reader_t rd;
  if (!cdid_sel_init_reader(reader_sexp, &rd, (mode >= 2) ? sm.is_date_model : 0)) return R_NilValue;
  cdid_sel_ref_t rf;
  if (!cdid_sel_init_ref(ref_sexp, &rf)) return R_NilValue;

  SEXP pair_s = cdid_sel_elt(spec_sexp, "pair");
  SEXP coh_s = cdid_sel_elt(spec_sexp, "cohorts");
  SEXP sids_s = cdid_sel_elt(spec_sexp, "selection_ids");
  if (TYPEOF(pair_s) != VECSXP || (TYPEOF(coh_s) != INTSXP && TYPEOF(coh_s) != REALSXP) ||
      (TYPEOF(sids_s) != INTSXP && TYPEOF(sids_s) != REALSXP)) return R_NilValue;

  int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
  int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
  int m = (int) XLENGTH(coh_s);
  int n_exp = Rf_asInteger(cdid_sel_elt(spec_sexp, "n_observers"));
  SEXP ref_r_s = cdid_sel_elt(spec_sexp, "reference");
  int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
  if (m < 1 || XLENGTH(sids_s) != m || n_exp < 1) return R_NilValue;

  const void *vmax = vmaxget();

  /* Cohort map: coh_pos[k] = 1..m (0 if cohort k is not in spec$cohorts) */
  int *coh_pos = (int *) R_alloc((size_t) (rd.J + 1), sizeof(int));
  memset(coh_pos, 0, (size_t) (rd.J + 1) * sizeof(int));
  int *spec_cohorts = (int *) R_alloc((size_t) m, sizeof(int));
  int *spec_sids = (int *) R_alloc((size_t) m, sizeof(int));
  for (int k = 0; k < m; ++k) {
    int ck = cdid_sel_int_at(coh_s, k);
    int sk = cdid_sel_int_at(sids_s, k);
    if (ck == NA_INTEGER || ck < 1 || ck > rd.J || sk == NA_INTEGER || sk < 1) {
      vmaxset(vmax); return R_NilValue;
    }
    coh_pos[ck] = k + 1;
    spec_cohorts[k] = ck;
    spec_sids[k] = sk;
  }

  int *obs_s_map = (int *) R_alloc((size_t) (rd.n_A + 1), sizeof(int));
  int *obs_t_map = (int *) R_alloc((size_t) (rd.n_A + 1), sizeof(int));
  if (!cdid_sel_build_pair_map(&rd, ds1, dt1, obs_s_map, obs_t_map)) {
    vmaxset(vmax); return R_NilValue;
  }

  /* Mode 4: custom selection callback fit pre-collector (canonical a_index order) */
  if (mode == 4) {
    if (n_exp < 2) { vmaxset(vmax); return R_NilValue; }
    if (TYPEOF(up_map_sexp) != VECSXP) { vmaxset(vmax); return R_NilValue; }
    int D_up = Rf_asInteger(cdid_sel_elt(up_map_sexp, "D"));
    SEXP pc_ch_list = cdid_sel_elt(up_map_sexp, "pair_channels");
    if (D_up < 0 || TYPEOF(pc_ch_list) != VECSXP) { vmaxset(vmax); return R_NilValue; }

    /* Collect observed rows across spec$cohorts, then sort by a_index[row] */
    int *a_to_row = (int *) R_alloc((size_t) (rd.n_A + 1), sizeof(int));
    int *a_to_kk = (int *) R_alloc((size_t) (rd.n_A + 1), sizeof(int));
    memset(a_to_row, 0, (size_t) (rd.n_A + 1) * sizeof(int));
    int n_found = 0;
    for (int r1 = 1; r1 <= rd.N; ++r1) {
      int ck = rd.cohort_index[r1 - 1];
      if (ck < 1 || ck > rd.J) continue;
      int kp = coh_pos[ck];
      if (kp == 0) continue;
      int ai = rd.a_index[r1 - 1];
      if (ai > 0 && ai <= rd.n_A && obs_s_map[ai] > 0) {
        a_to_row[ai] = r1;
        a_to_kk[ai] = kp;
        ++n_found;
      }
    }
    if (n_found != n_exp) { vmaxset(vmax); return R_NilValue; }

    /* Build sorted unique u_ids */
    int *u_pos_map = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
    memset(u_pos_map, 0, (size_t) (D_up + 1) * sizeof(int));
    for (int k = 0; k < m; ++k) {
      int q1 = spec_sids[k];
      if (q1 < 1 || q1 > XLENGTH(pc_ch_list)) { vmaxset(vmax); return R_NilValue; }
      SEXP ch_s = VECTOR_ELT(pc_ch_list, q1 - 1);
      if (TYPEOF(ch_s) != INTSXP) { vmaxset(vmax); return R_NilValue; }
      const int *ch = INTEGER(ch_s);
      for (R_xlen_t j = 0; j < XLENGTH(ch_s); ++j) {
        if (ch[j] < 1 || ch[j] > D_up) { vmaxset(vmax); return R_NilValue; }
        u_pos_map[ch[j]] = 1;
      }
    }
    int u_len = 0;
    for (int u = 1; u <= D_up; ++u) if (u_pos_map[u]) u_pos_map[u] = ++u_len;
    SEXP uids_out = PROTECT(Rf_allocVector(INTSXP, u_len));
    int u_fill = 0;
    for (int u = 1; u <= D_up; ++u) if (u_pos_map[u]) INTEGER(uids_out)[u_fill++] = u;

    SEXP rows_out = PROTECT(Rf_allocVector(INTSXP, n_exp));
    SEXP kk_out = PROTECT(Rf_allocVector(INTSXP, n_exp));
    SEXP fw_out = PROTECT(Rf_allocVector(REALSXP, n_exp));
    SEXP LD_out = PROTECT(Rf_allocMatrix(REALSXP, n_exp, u_len));
    SEXP y0_out = PROTECT(Rf_allocVector(REALSXP, n_exp));
    SEXP y1_out = PROTECT(Rf_allocVector(REALSXP, n_exp));
    memset(REAL(LD_out), 0, (size_t) n_exp * (size_t) u_len * sizeof(double));

    int *rows_p = INTEGER(rows_out);
    int *kk_p = INTEGER(kk_out);
    double *fw_p = REAL(fw_out);
    double *LD_p = REAL(LD_out);
    double *y0_p = REAL(y0_out);
    double *y1_p = REAL(y1_out);

    double *tlogit_buf = (double *) R_alloc((size_t) (u_len + 2), sizeof(double));
    int hz1 = -1, hz2 = -1, hw = -1, hzt = -1;
    int pos = 0;
    for (int ai = 1; ai <= rd.n_A; ++ai) {
      int r1 = a_to_row[ai];
      if (r1 == 0) continue;
      int kp = a_to_kk[ai];
      int q1 = spec_sids[kp - 1];
      double y0v = rd.y_obs[obs_s_map[ai] - 1];
      double y1v = rd.y_obs[obs_t_map[ai] - 1];
      if (!R_FINITE(y0v) || !R_FINITE(y1v)) { UNPROTECT(7); vmaxset(vmax); return R_NilValue; }

      double prob = 0.0, ps = 0.0, pt = 0.0, w_sel = 0.0, a_pc = 0.0;
      if (!cdid_sel_predict_pair_one(&sm, &rd, &rf, q1, r1, &hz1, &hz2, &hw, &prob, &ps, &pt, &w_sel, &a_pc)) {
        UNPROTECT(7); vmaxset(vmax); return R_NilValue;
      }
      SEXP ch_s = VECTOR_ELT(pc_ch_list, q1 - 1);
      int n_ch = (int) XLENGTH(ch_s);
      const int *ch = INTEGER(ch_s);
      if (n_ch < 1 || !cdid_sel_deriv_tlogit_one(&sm, &rd, &rf, q1, r1, prob, ps, pt, &hzt, tlogit_buf, n_ch - 1)) {
        UNPROTECT(7); vmaxset(vmax); return R_NilValue;
      }
      LD_p[pos + (size_t) (u_pos_map[ch[0]] - 1) * (size_t) n_exp] = 1.0 / a_pc;
      for (int j = 1; j < n_ch; ++j) {
        LD_p[pos + (size_t) (u_pos_map[ch[j]] - 1) * (size_t) n_exp] = tlogit_buf[j - 1];
      }
      rows_p[pos] = r1;
      kk_p[pos] = kp;
      fw_p[pos] = w_sel;
      y0_p[pos] = y0v;
      y1_p[pos] = y1v;
      ++pos;
    }

    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 7));
    SET_VECTOR_ELT(out_s, 0, rows_out);
    SET_VECTOR_ELT(out_s, 1, kk_out);
    SET_VECTOR_ELT(out_s, 2, fw_out);
    SET_VECTOR_ELT(out_s, 3, LD_out);
    SET_VECTOR_ELT(out_s, 4, y0_out);
    SET_VECTOR_ELT(out_s, 5, y1_out);
    SET_VECTOR_ELT(out_s, 6, uids_out);
    UNPROTECT(8);
    vmaxset(vmax);
    return out_s;
  }

  /* Modes 1, 2, 3: treatment logit fit + cross derivatives */
  if (TYPEOF(ctrl_sexp) != REALSXP || XLENGTH(ctrl_sexp) < 2 || TYPEOF(cols_sexp) != VECSXP) {
    vmaxset(vmax); return R_NilValue;
  }
  double eps = REAL(ctrl_sexp)[0];
  int maxit = (int) REAL(ctrl_sexp)[1];
  if (!R_FINITE(eps) || eps <= 0.0 || maxit < 1) { vmaxset(vmax); return R_NilValue; }

  SEXP Za_s = cdid_sel_elt(cols_sexp, "Z_a");
  SEXP trt_idx_s = cdid_sel_elt(cols_sexp, "trt_idx");
  SEXP or_idx_s = cdid_sel_elt(cols_sexp, "or_idx");

  const double *Za_ptr = NULL;
  int n_Za_rows = 0, n_Za_cols = 0;
  if (!rf.has_ref) {
    SEXP zdim = Rf_getAttrib(Za_s, R_DimSymbol);
    if (TYPEOF(Za_s) != REALSXP || TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) {
      vmaxset(vmax); return R_NilValue;
    }
    Za_ptr = REAL(Za_s);
    n_Za_rows = INTEGER(zdim)[0];
    n_Za_cols = INTEGER(zdim)[1];
  } else {
    Za_ptr = rf.x_Z;
    n_Za_rows = rf.n_x_rows;
    n_Za_cols = rf.p_x;
  }

  int p_trt = (trt_idx_s != R_NilValue && TYPEOF(trt_idx_s) == INTSXP) ? (int) XLENGTH(trt_idx_s) : n_Za_cols;
  const int *trt_idx1 = (trt_idx_s != R_NilValue && TYPEOF(trt_idx_s) == INTSXP) ? INTEGER(trt_idx_s) : NULL;
  int p_or = (mode == 3 && or_idx_s != R_NilValue && TYPEOF(or_idx_s) == INTSXP) ? (int) XLENGTH(or_idx_s) : 0;
  const int *or_idx1 = (p_or > 0) ? INTEGER(or_idx_s) : NULL;
  if (p_trt < 1 || p_trt > 128 || n_exp < p_trt) { vmaxset(vmax); return R_NilValue; }
  if (trt_idx1) {
    for (int j = 0; j < p_trt; ++j) if (trt_idx1[j] < 1 || trt_idx1[j] > n_Za_cols) { vmaxset(vmax); return R_NilValue; }
  }
  if (or_idx1) {
    for (int j = 0; j < p_or; ++j) if (or_idx1[j] < 1 || or_idx1[j] > n_Za_cols) { vmaxset(vmax); return R_NilValue; }
  }

  /* Legacy selection coefs for mode 1 */
  const double **leg_coefs = NULL;
  if (mode == 1) {
    if (TYPEOF(sel_sexp) != VECSXP || rd.Z_static == NULL || rd.d_Z < 1) { vmaxset(vmax); return R_NilValue; }
    leg_coefs = (const double **) R_alloc((size_t) m, sizeof(const double *));
    for (int k = 0; k < m; ++k) {
      int sj = spec_sids[k];
      if (sj < 1 || sj > XLENGTH(sel_sexp)) { vmaxset(vmax); return R_NilValue; }
      SEXP c_s = cdid_sel_elt(VECTOR_ELT(sel_sexp, sj - 1), "coefficients");
      if (TYPEOF(c_s) != REALSXP || XLENGTH(c_s) != rd.d_Z) { vmaxset(vmax); return R_NilValue; }
      leg_coefs[k] = REAL(c_s);
    }
  }

  double *X_trt = (double *) R_alloc((size_t) n_exp * (size_t) p_trt, sizeof(double));
  double *Z_or = (p_or > 0) ? (double *) R_alloc((size_t) n_exp * (size_t) p_or, sizeof(double)) : NULL;
  double *y_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
  double *w_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
  double *pp_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
  double *ps_vec = (mode >= 2 && sm.is_date_model) ? (double *) R_alloc((size_t) n_exp, sizeof(double)) : NULL;
  double *pt_vec = (mode >= 2 && sm.is_date_model) ? (double *) R_alloc((size_t) n_exp, sizeof(double)) : NULL;
  double *dy_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
  int *kk_vec = (int *) R_alloc((size_t) n_exp, sizeof(int));
  int *obs_rows = (int *) R_alloc((size_t) n_exp, sizeof(int));

  int pos = 0, hx = -1, hz1 = -1, hz2 = -1, hw = -1;
  for (int r1 = 1; r1 <= rd.N; ++r1) {
    int ck = rd.cohort_index[r1 - 1];
    if (ck < 1 || ck > rd.J) continue;
    int kp = coh_pos[ck];
    if (kp == 0) continue;
    int ai = rd.a_index[r1 - 1];
    if (ai <= 0 || ai > rd.n_A || obs_s_map[ai] == 0) continue;
    if (pos >= n_exp) { vmaxset(vmax); return R_NilValue; }

    double dy_i = rd.y_obs[obs_t_map[ai] - 1] - rd.y_obs[obs_s_map[ai] - 1];
    if (!R_FINITE(dy_i)) { vmaxset(vmax); return R_NilValue; }

    int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (ai - 1);
    if (xpos < 0 || xpos >= n_Za_rows) { vmaxset(vmax); return R_NilValue; }
    for (int j = 0; j < p_trt; ++j) {
      int cj = trt_idx1 ? (trt_idx1[j] - 1) : j;
      X_trt[pos + (size_t) j * (size_t) n_exp] = Za_ptr[xpos + (size_t) cj * (size_t) n_Za_rows];
    }
    for (int j = 0; j < p_or; ++j) {
      int cj = or_idx1[j] - 1;
      Z_or[pos + (size_t) j * (size_t) n_exp] = Za_ptr[xpos + (size_t) cj * (size_t) n_Za_rows];
    }

    double pp_i = 0.0, w_i = 0.0, ps_i = 0.0, pt_i = 0.0, a_pc = 0.0;
    if (mode == 1) {
      if (!cdid_sel_predict_legacy_one(&rd, &rf, 0, leg_coefs[kp - 1], rd.d_Z, r1, &hz1, &pp_i)) {
        vmaxset(vmax); return R_NilValue;
      }
      double p_fl = (pp_i > CDID_SEL_CLIP_EPS) ? pp_i : CDID_SEL_CLIP_EPS;
      w_i = rd.w_static[r1 - 1] * rd.a_vec[ck - 1] / (rd.rho[r1 - 1] * p_fl);
      if (!R_FINITE(w_i) || !(w_i > 0.0)) { vmaxset(vmax); return R_NilValue; }
    } else {
      if (!cdid_sel_predict_pair_one(&sm, &rd, &rf, spec_sids[kp - 1], r1, &hz1, &hz2, &hw,
                                     &pp_i, &ps_i, &pt_i, &w_i, &a_pc)) {
        vmaxset(vmax); return R_NilValue;
      }
      if (ps_vec) { ps_vec[pos] = ps_i; pt_vec[pos] = pt_i; }
    }

    y_vec[pos] = (kp == 1) ? 1.0 : 0.0;
    w_vec[pos] = w_i;
    pp_vec[pos] = pp_i;
    dy_vec[pos] = dy_i;
    kk_vec[pos] = kp;
    obs_rows[pos] = r1;
    ++pos;
  }
  if (pos != n_exp) { vmaxset(vmax); return R_NilValue; }

  double *beta_trt = (double *) R_alloc((size_t) p_trt, sizeof(double));
  double *R_chol = (double *) R_alloc((size_t) p_trt * (size_t) p_trt, sizeof(double));
  int iter_used = 0;
  if (!cdid_sel_logit_fit_core(X_trt, y_vec, w_vec, n_exp, p_trt, eps, maxit, rd.N_dbl,
                               beta_trt, R_chol, &iter_used)) {
    vmaxset(vmax); return R_NilValue;
  }

  /* Replay treatment logit values, normalizers, means */
  double *eta_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
  double *lv_prob = (double *) R_alloc((size_t) n_exp, sizeof(double));
  double *lv_act = (double *) R_alloc((size_t) n_exp, sizeof(double));
  double *wp_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
  double *wm_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
  BLAS_INT nexp_b = (BLAS_INT) n_exp, ptrt_b2 = (BLAS_INT) p_trt, one_b2 = 1;
  double one_d2 = 1.0, zero_d2 = 0.0;
  F77_CALL(dgemm)("N", "N", &nexp_b, &one_b2, &ptrt_b2, &one_d2,
                  X_trt, &nexp_b, beta_trt, &ptrt_b2, &zero_d2, eta_vec, &nexp_b FCONE FCONE);
  long double sum_wp_ld = 0.0L, sum_wm_ld = 0.0L, sum_wp_dy_ld = 0.0L, sum_wm_dy_ld = 0.0L;

  for (int i = 0; i < n_exp; ++i) {
    double eta = eta_vec[i];
    double pr = 0.0, od = 0.0, act = 0.0;
    cdid_sel_clip_logit_scalar(eta, &pr, &od, &act);
    if (!R_FINITE(eta) || !R_FINITE(pr) || !R_FINITE(od) || !(od > 0.0)) {
      vmaxset(vmax); return R_NilValue;
    }
    lv_prob[i] = pr;
    lv_act[i] = act;
    double wpi = w_vec[i] * y_vec[i];
    double wmi = w_vec[i] * (1.0 - y_vec[i]) * od;
    wp_vec[i] = wpi;
    wm_vec[i] = wmi;
    sum_wp_ld += (long double) wpi;
    sum_wm_ld += (long double) wmi;
    sum_wp_dy_ld += (long double) (wpi * dy_vec[i]);
    sum_wm_dy_ld += (long double) (wmi * dy_vec[i]);
  }
  double sum_wp = (double) sum_wp_ld;
  double sum_wm = (double) sum_wm_ld;
  double dp = sum_wp / rd.N_dbl;
  double dm = sum_wm / rd.N_dbl;
  if (!R_FINITE(dp) || !(dp > 0.0) || !R_FINITE(dm) || !(dm > 0.0)) {
    vmaxset(vmax); return R_NilValue;
  }
  double T_ipw = ((double) sum_wp_dy_ld) / sum_wp;
  double M_ipw = ((double) sum_wm_dy_ld) / sum_wm;
  if (!R_FINITE(T_ipw) || !R_FINITE(M_ipw)) { vmaxset(vmax); return R_NilValue; }

  /* DR residualized quantities when mode == 3 */
  double T_dr = T_ipw, M_dr = M_ipw;
  double *r_obs = NULL;
  double *d_bar_plus = NULL, *d_bar_minus = NULL, *d_bar_diff = NULL;
  if (mode == 3) {
    if (TYPEOF(or_info_sexp) != VECSXP) { vmaxset(vmax); return R_NilValue; }
    double reg_L = Rf_asReal(cdid_sel_elt(or_info_sexp, "L"));
    double reg_alpha = Rf_asReal(cdid_sel_elt(or_info_sexp, "alpha"));
    int has_icpt = Rf_asLogical(cdid_sel_elt(or_info_sexp, "has_icpt")) == 1;
    SEXP sh_beta_s = cdid_sel_elt(or_info_sexp, "beta");
    const double *sh_beta = (p_or > 0 && TYPEOF(sh_beta_s) == REALSXP) ? REAL(sh_beta_s) : NULL;
    if (p_or > 0 && (sh_beta == NULL || XLENGTH(sh_beta_s) != p_or)) {
      vmaxset(vmax); return R_NilValue;
    }
    r_obs = (double *) R_alloc((size_t) n_exp, sizeof(double));
    d_bar_plus = (p_or > 0) ? (double *) R_alloc((size_t) p_or, sizeof(double)) : NULL;
    d_bar_minus = (p_or > 0) ? (double *) R_alloc((size_t) p_or, sizeof(double)) : NULL;
    d_bar_diff = (p_or > 0) ? (double *) R_alloc((size_t) p_or, sizeof(double)) : NULL;
    double *lin_or_vec = (double *) R_alloc((size_t) n_exp, sizeof(double));
    long double *dbp_ld = (p_or > 0) ? R_allocLD((size_t) p_or) : NULL;
    long double *dbm_ld = (p_or > 0) ? R_allocLD((size_t) p_or) : NULL;
    if (p_or > 0) {
      double *D_or = (double *) R_alloc((size_t) n_exp * (size_t) p_or, sizeof(double));
      for (size_t idx = 0; idx < (size_t) n_exp * (size_t) p_or; ++idx) D_or[idx] = reg_L * Z_or[idx];
      BLAS_INT por_b = (BLAS_INT) p_or;
      F77_CALL(dgemm)("N", "N", &nexp_b, &one_b2, &por_b, &one_d2,
                      D_or, &nexp_b, sh_beta, &por_b, &zero_d2, lin_or_vec, &nexp_b FCONE FCONE);
      for (int j = 0; j < p_or; ++j) { dbp_ld[j] = 0.0L; dbm_ld[j] = 0.0L; }
    } else {
      memset(lin_or_vec, 0, (size_t) n_exp * sizeof(double));
    }
    long double sum_wp_r_ld = 0.0L, sum_wm_r_ld = 0.0L;
    for (int i = 0; i < n_exp; ++i) {
      double m_hat = (has_icpt ? reg_alpha : 0.0) + lin_or_vec[i];
      for (int j = 0; j < p_or; ++j) {
        double dij = reg_L * Z_or[i + (size_t) j * (size_t) n_exp];
        dbp_ld[j] += (long double) (wp_vec[i] * dij);
        dbm_ld[j] += (long double) (wm_vec[i] * dij);
      }
      double ri = dy_vec[i] - m_hat;
      r_obs[i] = ri;
      sum_wp_r_ld += (long double) (wp_vec[i] * ri);
      sum_wm_r_ld += (long double) (wm_vec[i] * ri);
    }
    T_dr = ((double) sum_wp_r_ld) / sum_wp;
    M_dr = ((double) sum_wm_r_ld) / sum_wm;
    if (!R_FINITE(T_dr) || !R_FINITE(M_dr)) { vmaxset(vmax); return R_NilValue; }
    for (int j = 0; j < p_or; ++j) {
      d_bar_plus[j] = ((double) dbp_ld[j]) / sum_wp;
      d_bar_minus[j] = ((double) dbm_ld[j]) / sum_wm;
      d_bar_diff[j] = d_bar_plus[j] - d_bar_minus[j];
      if (!R_FINITE(d_bar_plus[j]) || !R_FINITE(d_bar_minus[j])) { vmaxset(vmax); return R_NilValue; }
    }
  }

  /* Upstream channel mapping */
  int u_len = 0;
  int *u_pos_map = NULL;
  SEXP uids_out = R_NilValue;
  SEXP pc_ch_list = R_NilValue;
  int n_prot = 0;
  if (mode == 1) {
    u_len = m * (rd.d_Z + 1);
  } else {
    if (TYPEOF(up_map_sexp) != VECSXP) { vmaxset(vmax); return R_NilValue; }
    int D_up = Rf_asInteger(cdid_sel_elt(up_map_sexp, "D"));
    pc_ch_list = cdid_sel_elt(up_map_sexp, "pair_channels");
    if (D_up < 0 || TYPEOF(pc_ch_list) != VECSXP) { vmaxset(vmax); return R_NilValue; }
    u_pos_map = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
    memset(u_pos_map, 0, (size_t) (D_up + 1) * sizeof(int));
    for (int k = 0; k < m; ++k) {
      int q1 = spec_sids[k];
      if (q1 < 1 || q1 > XLENGTH(pc_ch_list)) { vmaxset(vmax); return R_NilValue; }
      SEXP ch_s = VECTOR_ELT(pc_ch_list, q1 - 1);
      if (TYPEOF(ch_s) != INTSXP) { vmaxset(vmax); return R_NilValue; }
      const int *ch = INTEGER(ch_s);
      for (R_xlen_t j = 0; j < XLENGTH(ch_s); ++j) {
        if (ch[j] < 1 || ch[j] > D_up) { vmaxset(vmax); return R_NilValue; }
        u_pos_map[ch[j]] = 1;
      }
    }
    for (int u = 1; u <= D_up; ++u) if (u_pos_map[u]) u_pos_map[u] = ++u_len;
    uids_out = PROTECT(Rf_allocVector(INTSXP, u_len));
    ++n_prot;
    int uf = 0;
    for (int u = 1; u <= D_up; ++u) if (u_pos_map[u]) INTEGER(uids_out)[uf++] = u;
  }

  SEXP C_s = PROTECT(Rf_allocMatrix(REALSXP, p_trt, u_len)); ++n_prot;
  SEXP U_main_s = PROTECT(Rf_allocMatrix(REALSXP, 4, u_len)); ++n_prot;
  SEXP V_main_s = PROTECT(Rf_allocMatrix(REALSXP, 4, p_trt)); ++n_prot;
  SEXP U_ipw_s = R_NilValue, V_ipw_s = R_NilValue;
  if (mode == 3) {
    U_ipw_s = PROTECT(Rf_allocMatrix(REALSXP, 4, u_len)); ++n_prot;
    V_ipw_s = PROTECT(Rf_allocMatrix(REALSXP, 4, p_trt)); ++n_prot;
    memset(REAL(U_ipw_s), 0, (size_t) 4 * (size_t) u_len * sizeof(double));
    memset(REAL(V_ipw_s), 0, (size_t) 4 * (size_t) p_trt * sizeof(double));
  }
  double *C_p = REAL(C_s);
  double *Um_p = REAL(U_main_s);
  double *Vm_p = REAL(V_main_s);
  double *Ui_p = (mode == 3) ? REAL(U_ipw_s) : NULL;
  double *Vi_p = (mode == 3) ? REAL(V_ipw_s) : NULL;
  memset(C_p, 0, (size_t) p_trt * (size_t) u_len * sizeof(double));
  memset(Um_p, 0, (size_t) 4 * (size_t) u_len * sizeof(double));
  memset(Vm_p, 0, (size_t) 4 * (size_t) p_trt * sizeof(double));

  int max_ch = u_len + rd.d_Z + 4;
  double *ld_buf = (double *) R_alloc((size_t) max_ch, sizeof(double));
  double *c_acc = (double *) R_alloc((size_t) p_trt * (size_t) (max_ch > 0 ? max_ch : 1), sizeof(double));
  double *um_acc = (double *) R_alloc((size_t) 4 * (size_t) (max_ch > 0 ? max_ch : 1), sizeof(double));
  double *ui_acc = (mode == 3) ? (double *) R_alloc((size_t) 4 * (size_t) (max_ch > 0 ? max_ch : 1), sizeof(double)) : NULL;
  int hzt = -1;

  for (int kp = 1; kp <= m; ++kp) {
    int n_rr = 0;
    for (int i = 0; i < n_exp; ++i) {
      if (kk_vec[i] == kp) ++n_rr;
    }
    if (n_rr == 0) continue;

    int n_ch = 0;
    const int *ch_ids = NULL;
    double a_pc = 0.0;
    int q1 = 0;
    if (mode == 1) {
      int ck = spec_cohorts[kp - 1];
      n_ch = rd.d_Z + 1;
      a_pc = rd.a_vec[ck - 1];
    } else {
      q1 = spec_sids[kp - 1];
      SEXP pc = VECTOR_ELT(sm.pair_cohorts, q1 - 1);
      int ck = Rf_asInteger(cdid_sel_elt(pc, "cohort"));
      SEXP ref_q_s = cdid_sel_elt(pc, "reference");
      int rq = (ref_q_s != R_NilValue) ? Rf_asInteger(ref_q_s) : 0;
      if (!rf.has_ref) a_pc = rd.a_vec[ck - 1];
      else if (!cdid_sel_ref_stat_one(&rf, rq, ck, 0, &a_pc)) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
      SEXP ch_s = VECTOR_ELT(pc_ch_list, q1 - 1);
      n_ch = (int) XLENGTH(ch_s);
      ch_ids = INTEGER(ch_s);
      if (n_ch < 1) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
    }

    const void *vm_kp = vmaxget();
    double *ld_mat = (double *) R_alloc((size_t) n_rr * (size_t) n_ch, sizeof(double));
    double *sc_sub = (double *) R_alloc((size_t) n_rr * (size_t) p_trt, sizeof(double));
    double *rm_sub = (double *) R_alloc((size_t) n_rr * 4, sizeof(double));
    double *ri_sub = (mode == 3) ? (double *) R_alloc((size_t) n_rr * 4, sizeof(double)) : NULL;

    int t_idx = 0;
    for (int i = 0; i < n_exp; ++i) {
      if (kk_vec[i] != kp) continue;
      int r1 = obs_rows[i];
      if (mode == 1) {
        ld_buf[0] = 1.0 / a_pc;
        double act_q = (pp_vec[i] > CDID_SEL_CLIP_EPS) ? 1.0 : 0.0;
        double fac = -((1.0 - pp_vec[i]) * act_q);
        for (int j = 0; j < rd.d_Z; ++j) {
          ld_buf[1 + j] = fac * rd.Z_static[(r1 - 1) + (size_t) j * (size_t) rd.N];
        }
      } else {
        double ps_i = ps_vec ? ps_vec[i] : 0.0;
        double pt_i = pt_vec ? pt_vec[i] : 0.0;
        if (!cdid_sel_deriv_tlogit_one(&sm, &rd, &rf, q1, r1, pp_vec[i], ps_i, pt_i,
                                       &hzt, ld_buf + 1, n_ch - 1)) {
          UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue;
        }
        ld_buf[0] = 1.0 / a_pc;
      }

      double wpi = wp_vec[i], wmi = wm_vec[i];
      double res_ipw0 = wpi * (dy_vec[i] - T_ipw) / dp;
      double res_ipw1 = wmi * (dy_vec[i] - M_ipw) / dm;
      double res_m0 = (mode == 3) ? (wpi * (r_obs[i] - T_dr) / dp) : res_ipw0;
      double res_m1 = (mode == 3) ? (wmi * (r_obs[i] - M_dr) / dm) : res_ipw1;
      double sc_fac = w_vec[i] * (y_vec[i] - lv_prob[i]);

      for (int c = 0; c < n_ch; ++c) {
        ld_mat[t_idx + (size_t) c * (size_t) n_rr] = ld_buf[c];
      }
      for (int j = 0; j < p_trt; ++j) {
        sc_sub[t_idx + (size_t) j * (size_t) n_rr] = sc_fac * X_trt[i + (size_t) j * (size_t) n_exp];
      }
      rm_sub[t_idx + 0 * (size_t) n_rr] = res_m0;
      rm_sub[t_idx + 1 * (size_t) n_rr] = res_m1;
      rm_sub[t_idx + 2 * (size_t) n_rr] = wpi;
      rm_sub[t_idx + 3 * (size_t) n_rr] = wmi;
      if (mode == 3) {
        ri_sub[t_idx + 0 * (size_t) n_rr] = res_ipw0;
        ri_sub[t_idx + 1 * (size_t) n_rr] = res_ipw1;
        ri_sub[t_idx + 2 * (size_t) n_rr] = wpi;
        ri_sub[t_idx + 3 * (size_t) n_rr] = wmi;
      }
      ++t_idx;
    }

    if (!cdid_sel_all_finite(ld_mat, (R_xlen_t) n_rr * n_ch) ||
        !cdid_sel_all_finite(sc_sub, (R_xlen_t) n_rr * p_trt) ||
        !cdid_sel_all_finite(rm_sub, (R_xlen_t) n_rr * 4) ||
        (mode == 3 && !cdid_sel_all_finite(ri_sub, (R_xlen_t) n_rr * 4))) {
      UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue;
    }

    BLAS_INT k_b = (BLAS_INT) n_rr;
    BLAS_INT nch_b = (BLAS_INT) n_ch;
    BLAS_INT ptrt_b = (BLAS_INT) p_trt;
    BLAS_INT four_b = 4;
    double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)("T", "N", &ptrt_b, &nch_b, &k_b, &one, sc_sub, &k_b, ld_mat, &k_b, &zero, c_acc, &ptrt_b FCONE FCONE);
    F77_CALL(dgemm)("T", "N", &four_b, &nch_b, &k_b, &one, rm_sub, &k_b, ld_mat, &k_b, &zero, um_acc, &four_b FCONE FCONE);
    if (mode == 3) {
      F77_CALL(dgemm)("T", "N", &four_b, &nch_b, &k_b, &one, ri_sub, &k_b, ld_mat, &k_b, &zero, ui_acc, &four_b FCONE FCONE);
    }
    vmaxset(vm_kp);

    for (int c = 0; c < n_ch; ++c) {
      int dst_u = (mode == 1) ? ((kp - 1) * (rd.d_Z + 1) + c) : (u_pos_map[ch_ids[c]] - 1);
      for (int j = 0; j < p_trt; ++j) {
        C_p[j + (size_t) dst_u * (size_t) p_trt] += c_acc[j + (size_t) c * (size_t) p_trt] / rd.N_dbl;
      }
      Um_p[0 + 4 * dst_u] += um_acc[0 + 4 * c] / rd.N_dbl;
      Um_p[1 + 4 * dst_u] += um_acc[1 + 4 * c] / rd.N_dbl;
      Um_p[2 + 4 * dst_u] += um_acc[2 + 4 * c] / rd.N_dbl;
      Um_p[3 + 4 * dst_u] += um_acc[3 + 4 * c] / rd.N_dbl;
      if (mode == 3) {
        Ui_p[0 + 4 * dst_u] += ui_acc[0 + 4 * c] / rd.N_dbl;
        Ui_p[1 + 4 * dst_u] += ui_acc[1 + 4 * c] / rd.N_dbl;
        Ui_p[2 + 4 * dst_u] += ui_acc[2 + 4 * c] / rd.N_dbl;
        Ui_p[3 + 4 * dst_u] += ui_acc[3 + 4 * c] / rd.N_dbl;
      }
    }
  }

  double ndm = rd.N_dbl * dm;
  for (int j = 0; j < p_trt; ++j) {
    const double *xj = X_trt + (size_t) j * (size_t) n_exp;
    long double s_v2m = 0.0L, s_v4 = 0.0L, s_v2i = 0.0L;
    for (int i = 0; i < n_exp; ++i) {
      double wm_act = wm_vec[i] * lv_act[i];
      double xij = xj[i];
      double diff_m = (mode == 3) ? (r_obs[i] - M_dr) : (dy_vec[i] - M_ipw);
      s_v2m += (long double) ((wm_act * diff_m) * xij);
      s_v4 += (long double) (wm_act * xij);
      if (mode == 3) s_v2i += (long double) ((wm_act * (dy_vec[i] - M_ipw)) * xij);
    }
    Vm_p[1 + 4 * j] = ((double) s_v2m) / ndm;
    Vm_p[3 + 4 * j] = ((double) s_v4) / rd.N_dbl;
    if (mode == 3) {
      Vi_p[1 + 4 * j] = ((double) s_v2i) / ndm;
      Vi_p[3 + 4 * j] = Vm_p[3 + 4 * j];
    }
  }

  for (int idx = 0; idx < p_trt * u_len; ++idx) if (!R_FINITE(C_p[idx])) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
  for (int idx = 0; idx < 4 * u_len; ++idx) if (!R_FINITE(Um_p[idx])) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
  for (int idx = 0; idx < 4 * p_trt; ++idx) if (!R_FINITE(Vm_p[idx])) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
  if (mode == 3) {
    for (int idx = 0; idx < 4 * u_len; ++idx) if (!R_FINITE(Ui_p[idx])) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
    for (int idx = 0; idx < 4 * p_trt; ++idx) if (!R_FINITE(Vi_p[idx])) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
  }

  SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, p_trt)); ++n_prot;
  memcpy(REAL(beta_s), beta_trt, (size_t) p_trt * sizeof(double));
  SEXP Hchol_s = PROTECT(Rf_allocMatrix(REALSXP, p_trt, p_trt)); ++n_prot;
  memcpy(REAL(Hchol_s), R_chol, (size_t) p_trt * (size_t) p_trt * sizeof(double));
  SEXP clip_s = PROTECT(cdid_sel_make_clip_summary(lv_prob, pp_vec, y_vec, w_vec, n_exp, 1)); ++n_prot;

  SEXP means_s = R_NilValue, dbp_s = R_NilValue, dbm_s = R_NilValue, dbd_s = R_NilValue;
  if (mode != 3) {
    means_s = PROTECT(Rf_allocVector(REALSXP, 4)); ++n_prot;
    REAL(means_s)[0] = T_ipw; REAL(means_s)[1] = M_ipw; REAL(means_s)[2] = dp; REAL(means_s)[3] = dm;
  } else {
    means_s = PROTECT(Rf_allocVector(REALSXP, 12)); ++n_prot;
    double mu_dr = T_dr - M_dr;
    double mu_ipw = T_ipw - M_ipw;
    double dot_dbd = 0.0;
    if (p_or > 0) {
      const double *sh_beta = REAL(cdid_sel_elt(or_info_sexp, "beta"));
      for (int j = 0; j < p_or; ++j) dot_dbd += d_bar_diff[j] * sh_beta[j];
    }
    double mu_dr_decomp = mu_ipw - dot_dbd;
    if (!R_FINITE(mu_dr_decomp)) { UNPROTECT(n_prot); vmaxset(vmax); return R_NilValue; }
    REAL(means_s)[0] = T_dr; REAL(means_s)[1] = M_dr; REAL(means_s)[2] = dp; REAL(means_s)[3] = dm;
    REAL(means_s)[4] = T_dr; REAL(means_s)[5] = M_dr; REAL(means_s)[6] = mu_dr;
    REAL(means_s)[7] = T_ipw; REAL(means_s)[8] = M_ipw; REAL(means_s)[9] = mu_ipw;
    REAL(means_s)[10] = T_dr; REAL(means_s)[11] = mu_dr_decomp;

    dbp_s = PROTECT(Rf_allocVector(REALSXP, p_or)); ++n_prot;
    dbm_s = PROTECT(Rf_allocVector(REALSXP, p_or)); ++n_prot;
    dbd_s = PROTECT(Rf_allocVector(REALSXP, p_or)); ++n_prot;
    if (p_or > 0) {
      memcpy(REAL(dbp_s), d_bar_plus, (size_t) p_or * sizeof(double));
      memcpy(REAL(dbm_s), d_bar_minus, (size_t) p_or * sizeof(double));
      memcpy(REAL(dbd_s), d_bar_diff, (size_t) p_or * sizeof(double));
    }
  }

  /* Return list(beta, H_chol, iter, C, U_main, V_main, means, clip, uids, U_ipw, V_ipw, dbp, dbm, dbd) */
  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 14)); ++n_prot;
  SET_VECTOR_ELT(out_s, 0, beta_s);
  SET_VECTOR_ELT(out_s, 1, Hchol_s);
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarInteger(iter_used));
  SET_VECTOR_ELT(out_s, 3, C_s);
  SET_VECTOR_ELT(out_s, 4, U_main_s);
  SET_VECTOR_ELT(out_s, 5, V_main_s);
  SET_VECTOR_ELT(out_s, 6, means_s);
  SET_VECTOR_ELT(out_s, 7, clip_s);
  SET_VECTOR_ELT(out_s, 8, uids_out);
  SET_VECTOR_ELT(out_s, 9, U_ipw_s);
  SET_VECTOR_ELT(out_s, 10, V_ipw_s);
  SET_VECTOR_ELT(out_s, 11, dbp_s);
  SET_VECTOR_ELT(out_s, 12, dbm_s);
  SET_VECTOR_ELT(out_s, 13, dbd_s);

  UNPROTECT(n_prot);
  vmaxset(vmax);
  return out_s;
}

/* ---------------------------------------------------------------------------
   7. C_cdid_sel_support_collect(mode_sexp, reader_sexp, spec_sexp, sel_sexp,
                                 up_map_sexp, q_id_sexp, ref_sexp, cols_sexp)
   modes:
     1 = or_corollary1_support (legacy CP04 OR)
     2 = or_selection_support  (general selection OR, static & D60)
     3 = dr_selection_support  (general selection DR, static & D60)
   --------------------------------------------------------------------------- */
static int cdid_sel_qr_compress_buf(
    double *A,
    int m,
    int k,
    double *dst_fac,
    int *jpvt,
    double *vn) {
  (void) jpvt;
  (void) vn;
  if (m <= 0 || k <= 0) return 0;
  R_xlen_t tot = (R_xlen_t) m * (R_xlen_t) k;
  for (R_xlen_t idx = 0; idx < tot; ++idx) {
    if (!R_FINITE(A[idx])) return 0;
  }
  const void *vm = vmaxget();
  int r = (m < k) ? m : k;
  BLAS_INT mb = (BLAS_INT) m;
  BLAS_INT kb = (BLAS_INT) k;
  BLAS_INT *piv = (BLAS_INT *) R_alloc((size_t) k, sizeof(BLAS_INT));
  for (int j = 0; j < k; ++j) piv[j] = 0;
  double *tau = (double *) R_alloc((size_t) r, sizeof(double));
  BLAS_INT lwork = -1, info = 0;
  double wkopt = 0.0;
  F77_CALL(dgeqp3)(&mb, &kb, A, &mb, piv, tau, &wkopt, &lwork, &info);
  if (info != 0 || !R_FINITE(wkopt)) {
    vmaxset(vm);
    return 0;
  }
  lwork = (BLAS_INT) wkopt;
  if (lwork < 3 * kb + 1) lwork = 3 * kb + 1;
  double *work = (double *) R_alloc((size_t) lwork, sizeof(double));
  F77_CALL(dgeqp3)(&mb, &kb, A, &mb, piv, tau, work, &lwork, &info);
  if (info != 0) {
    vmaxset(vm);
    return 0;
  }
  memset(dst_fac, 0, (size_t) r * (size_t) k * sizeof(double));
  for (int col = 0; col < k; ++col) {
    int orig_col = (int) piv[col] - 1;
    if (orig_col < 0 || orig_col >= k) {
      vmaxset(vm);
      return 0;
    }
    double *dst = dst_fac + (size_t) orig_col * (size_t) r;
    const double *sc = A + (size_t) col * (size_t) m;
    int rlim = (col + 1 < r) ? (col + 1) : r;
    for (int i = 0; i < rlim; ++i) {
      if (!R_FINITE(sc[i])) {
        vmaxset(vm);
        return 0;
      }
      dst[i] = sc[i];
    }
  }
  vmaxset(vm);
  return 1;
}

SEXP C_cdid_sel_support_collect(
    SEXP mode_sexp,
    SEXP reader_sexp,
    SEXP spec_sexp,
    SEXP sel_sexp,
    SEXP up_map_sexp,
    SEXP q_id_sexp,
    SEXP ref_sexp,
    SEXP cols_sexp) {
  int mode = Rf_asInteger(mode_sexp);
  if (mode < 1 || mode > 3 || TYPEOF(spec_sexp) != VECSXP || TYPEOF(cols_sexp) != VECSXP) return R_NilValue;

  cdid_sel_model_t sm;
  memset(&sm, 0, sizeof(sm));
  if (mode >= 2) {
    if (!cdid_sel_init_model(sel_sexp, &sm)) return R_NilValue;
  }
  cdid_sel_reader_t rd;
  if (!cdid_sel_init_reader(reader_sexp, &rd, (mode >= 2) ? sm.is_date_model : 0)) return R_NilValue;
  cdid_sel_ref_t rf;
  if (!cdid_sel_init_ref(ref_sexp, &rf)) return R_NilValue;

  int ck = Rf_asInteger(cdid_sel_elt(spec_sexp, "cohort"));
  SEXP pair_s = cdid_sel_elt(spec_sexp, "pair");
  if (ck < 1 || ck > rd.J || TYPEOF(pair_s) != VECSXP) return R_NilValue;
  int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
  int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
  SEXP ref_r_s = cdid_sel_elt(spec_sexp, "reference");
  int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

  int first0 = cdid_sel_int_at(rd.cohort_start_s, ck - 1) - 1;
  int last0 = cdid_sel_int_at(rd.cohort_end_s, ck - 1) - 1;
  int coh_len = last0 - first0 + 1;
  int max_chunk = Rf_asInteger(cdid_sel_elt(cols_sexp, "max_chunk_rows"));
  if (coh_len < 1 || max_chunk < 1) return R_NilValue;

  SEXP Za_s = cdid_sel_elt(cols_sexp, "Z_a");
  SEXP z_idx_s = cdid_sel_elt(cols_sexp, "z_idx");
  SEXP slope_idx_s = cdid_sel_elt(cols_sexp, "slope_idx");
  int icpt_flag = Rf_asInteger(cdid_sel_elt(cols_sexp, "icpt_flag")) ? 1 : 0;
  int icpt_col = Rf_asInteger(cdid_sel_elt(cols_sexp, "icpt_col"));
  double L_val = Rf_asReal(cdid_sel_elt(cols_sexp, "L_val"));
  if (TYPEOF(slope_idx_s) != INTSXP || !R_FINITE(L_val)) return R_NilValue;
  int p_slope = (int) XLENGTH(slope_idx_s);
  const int *slope_idx1 = INTEGER(slope_idx_s);
  if ((!icpt_flag && p_slope < 1) || p_slope > 256) return R_NilValue;

  const double *Z_src = NULL;
  int n_Z_rows = 0, n_Z_cols = 0;
  if (!rf.has_ref) {
    SEXP zdim = Rf_getAttrib(Za_s, R_DimSymbol);
    if (TYPEOF(Za_s) != REALSXP || TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return R_NilValue;
    Z_src = REAL(Za_s);
    n_Z_rows = INTEGER(zdim)[0];
    n_Z_cols = INTEGER(zdim)[1];
  } else {
    Z_src = rf.x_Z;
    n_Z_rows = rf.n_x_rows;
    n_Z_cols = rf.p_x;
  }

  int p_Z_out = (z_idx_s != R_NilValue && TYPEOF(z_idx_s) == INTSXP) ? (int) XLENGTH(z_idx_s) : n_Z_cols;
  const int *z_idx1 = (z_idx_s != R_NilValue && TYPEOF(z_idx_s) == INTSXP) ? INTEGER(z_idx_s) : NULL;
  if (z_idx1) {
    for (int j = 0; j < p_Z_out; ++j) if (z_idx1[j] < 1 || z_idx1[j] > n_Z_cols) return R_NilValue;
  }
  if (icpt_flag && (icpt_col < 1 || icpt_col > p_Z_out)) return R_NilValue;
  for (int j = 0; j < p_slope; ++j) if (slope_idx1[j] < 1 || slope_idx1[j] > p_Z_out) return R_NilValue;

  double a_val = 0.0;
  if (!rf.has_ref) a_val = rd.a_vec[ck - 1];
  else if (!cdid_sel_ref_stat_one(&rf, r_id, ck, 0, &a_val)) return R_NilValue;
  if (!R_FINITE(a_val) || !(a_val > 0.0)) return R_NilValue;

  int dq = 0, q1 = 0;
  const double *leg_coef = NULL;
  if (mode == 1) {
    SEXP c_s = cdid_sel_elt(sel_sexp, "coefficients");
    if (TYPEOF(c_s) != REALSXP || rd.Z_static == NULL || XLENGTH(c_s) != rd.d_Z) return R_NilValue;
    leg_coef = REAL(c_s);
    dq = rd.d_Z;
  } else {
    q1 = Rf_asInteger(q_id_sexp);
    SEXP pc_ch_list = cdid_sel_elt(up_map_sexp, "pair_channels");
    if (TYPEOF(pc_ch_list) != VECSXP || q1 < 1 || q1 > XLENGTH(pc_ch_list)) return R_NilValue;
    SEXP ch_s = VECTOR_ELT(pc_ch_list, q1 - 1);
    if (TYPEOF(ch_s) != INTSXP || XLENGTH(ch_s) < 1) return R_NilValue;
    dq = (int) XLENGTH(ch_s) - 1;
  }

  const void *vmax = vmaxget();
  int *obs_s_map = (int *) R_alloc((size_t) (rd.n_A + 1), sizeof(int));
  int *obs_t_map = (int *) R_alloc((size_t) (rd.n_A + 1), sizeof(int));
  if (!cdid_sel_build_pair_map(&rd, ds1, dt1, obs_s_map, obs_t_map)) {
    vmaxset(vmax); return R_NilValue;
  }

  int k_aug = p_slope + 1 + icpt_flag;
  int max_buf = (max_chunk < coh_len) ? max_chunk : coh_len;
  double *buf_D = (p_slope > 0) ? (double *) R_alloc((size_t) max_buf * (size_t) p_slope, sizeof(double)) : NULL;
  double *buf_dy = (double *) R_alloc((size_t) max_buf, sizeof(double));
  double *buf_w = (double *) R_alloc((size_t) max_buf, sizeof(double));
  double *aug_work = (double *) R_alloc((size_t) (max_buf + k_aug) * (size_t) k_aug, sizeof(double));
  double *fac_cur = (double *) R_alloc((size_t) k_aug * (size_t) k_aug, sizeof(double));
  double *fac_chunk = (double *) R_alloc((size_t) k_aug * (size_t) k_aug, sizeof(double));
  double *fac_merged = (double *) R_alloc((size_t) k_aug * (size_t) k_aug, sizeof(double));
  double *stack_work = (double *) R_alloc((size_t) (2 * k_aug) * (size_t) k_aug, sizeof(double));
  int *jpvt = (int *) R_alloc((size_t) k_aug, sizeof(int));
  double *vn = (double *) R_alloc((size_t) k_aug, sizeof(double));

  double *sumZ = (p_slope > 0) ? (double *) R_alloc((size_t) p_slope, sizeof(double)) : NULL;
  long double *ch_sumZ = (p_slope > 0) ? R_allocLD((size_t) p_slope) : NULL;
  double *colmax = (p_slope > 0) ? (double *) R_alloc((size_t) p_slope, sizeof(double)) : NULL;
  double *anchor = (p_slope > 0) ? (double *) R_alloc((size_t) p_slope, sizeof(double)) : NULL;
  double *ch_anchor = (p_slope > 0) ? (double *) R_alloc((size_t) p_slope, sizeof(double)) : NULL;
  int *constant = (p_slope > 0) ? (int *) R_alloc((size_t) p_slope, sizeof(int)) : NULL;
  int *ch_const = (p_slope > 0) ? (int *) R_alloc((size_t) p_slope, sizeof(int)) : NULL;
  double *z_row = (p_Z_out > 0) ? (double *) R_alloc((size_t) p_Z_out, sizeof(double)) : NULL;
  int *ch_r1 = (int *) R_alloc((size_t) max_buf, sizeof(int));
  double *buf_Z_sl = (p_slope > 0) ? (double *) R_alloc((size_t) max_buf * (size_t) p_slope, sizeof(double)) : NULL;
  double *ch_pr = (double *) R_alloc((size_t) max_buf, sizeof(double));
  double *ch_ps = (mode >= 2) ? (double *) R_alloc((size_t) max_buf, sizeof(double)) : NULL;
  double *ch_pt = (mode >= 2) ? (double *) R_alloc((size_t) max_buf, sizeof(double)) : NULL;
  double *prob_buf = (mode <= 2) ? (double *) R_alloc((size_t) coh_len, sizeof(double)) : NULL;
  double *q_all_buf = (mode <= 2) ? (double *) R_alloc((size_t) coh_len, sizeof(double)) : NULL;

  int m_dim = p_slope + 2;
  size_t tot_mom = (size_t) m_dim * (size_t) m_dim * (size_t) dq;
  SEXP mom_s = PROTECT(Rf_alloc3DArray(REALSXP, m_dim, m_dim, dq));
  double *mom_p = REAL(mom_s);
  memset(mom_p, 0, tot_mom * sizeof(double));
  double *ch_mom = (dq > 0) ? (double *) R_alloc(tot_mom, sizeof(double)) : NULL;
  double *x_chunk = (dq > 0) ? (double *) R_alloc((size_t) max_buf * (size_t) m_dim, sizeof(double)) : NULL;
  double *wx_chunk = (dq > 0) ? (double *) R_alloc((size_t) max_buf * (size_t) m_dim, sizeof(double)) : NULL;
  double *tl_chunk = (dq > 0) ? (double *) R_alloc((size_t) max_buf * (size_t) dq, sizeof(double)) : NULL;

  double *tlogit_buf = (double *) R_alloc((size_t) (dq + 2), sizeof(double));
  int hx = -1, hzt = -1;

  for (int j = 0; j < p_slope; ++j) {
    sumZ[j] = 0.0;
    colmax[j] = 0.0;
    constant[j] = 1;
    anchor[j] = 0.0;
  }
  double weight_mass = 0.0, sumY = 0.0;
  int n_obs_tot = 0, fac_rows = 0, acc_n_qr = 0, max_merge_rows = 0, n_comp = 0;

  for (int start0 = first0; start0 <= last0; start0 += max_chunk) {
    int end0 = start0 + max_chunk - 1;
    if (end0 > last0) end0 = last0;

    long double ch_w_mass = 0.0L, ch_sumY = 0.0L;
    for (int j = 0; j < p_slope; ++j) ch_sumZ[j] = 0.0L;
    if (dq > 0) memset(ch_mom, 0, tot_mom * sizeof(double));
    int filled = 0;

    for (int idx = start0; idx <= end0; ++idx) {
      int r1 = rd.cohort_order[idx];
      int ai = rd.a_index[r1 - 1];
      if (ai <= 0 || ai > rd.n_A || obs_s_map[ai] <= 0) continue;

      double dy = rd.y_obs[obs_t_map[ai] - 1] - rd.y_obs[obs_s_map[ai] - 1];
      if (!R_FINITE(dy)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

      int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (ai - 1);
      if (xpos < 0 || xpos >= n_Z_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      for (int j = 0; j < p_Z_out; ++j) {
        int cj = z_idx1 ? (z_idx1[j] - 1) : j;
        double zv = Z_src[xpos + (size_t) cj * (size_t) n_Z_rows];
        if (!R_FINITE(zv)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        z_row[j] = zv;
      }
      if (icpt_flag && z_row[icpt_col - 1] != 1.0) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      for (int j = 0; j < p_slope; ++j) {
        double zj = z_row[slope_idx1[j] - 1];
        double dj = L_val * zj;
        if (!R_FINITE(dj)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        buf_Z_sl[filled + (size_t) j * (size_t) max_buf] = zj;
        buf_D[filled + (size_t) j * (size_t) max_buf] = dj;
      }
      ch_r1[filled] = r1;
      buf_dy[filled] = dy;
      ++filled;
    }

    if (filled > 0) {
      int m_fl = filled;
      if (mode == 1) {
        const void *vm_m1 = vmaxget();
        double *Z_leg = (double *) R_alloc((size_t) m_fl * (size_t) dq, sizeof(double));
        double *eta_leg = (double *) R_alloc((size_t) m_fl, sizeof(double));
        for (int o = 0; o < m_fl; ++o) {
          int r1 = ch_r1[o];
          if (rd.U_logical[r1 - 1] != 1) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          for (int j = 0; j < dq; ++j) {
            double v = rd.Z_static[(r1 - 1) + (size_t) j * (size_t) rd.N];
            if (!R_FINITE(v)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            Z_leg[o + (size_t) j * (size_t) m_fl] = v;
          }
        }
        BLAS_INT m_b = (BLAS_INT) m_fl, n_b = 1, k_b = (BLAS_INT) dq;
        double one_d = 1.0, zero_d = 0.0;
        F77_CALL(dgemm)("N", "N", &m_b, &n_b, &k_b, &one_d, Z_leg, &m_b, leg_coef, &k_b, &zero_d, eta_leg, &m_b FCONE FCONE);
        for (int o = 0; o < m_fl; ++o) {
          int r1 = ch_r1[o];
          double pr = 0.0, od = 0.0, act_l = 0.0;
          cdid_sel_clip_logit_scalar(eta_leg[o], &pr, &od, &act_l);
          if (!R_FINITE(pr) || !(pr > 0.0) || pr >= 1.0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          double p_fl = (pr > CDID_SEL_CLIP_EPS) ? pr : CDID_SEL_CLIP_EPS;
          double q_w = rd.w_static[r1 - 1] * (a_val / (rd.rho[r1 - 1] * p_fl));
          if (!R_FINITE(q_w) || !(q_w > 0.0)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          ch_pr[o] = pr;
          buf_w[o] = q_w;
          double act = (pr > CDID_SEL_CLIP_EPS) ? 1.0 : 0.0;
          double fac = -((1.0 - pr) * act);
          for (int l = 0; l < dq; ++l) {
            tl_chunk[o + (size_t) l * (size_t) max_buf] = fac * Z_leg[o + (size_t) l * (size_t) m_fl];
          }
        }
        vmaxset(vm_m1);
      } else {
        double a_pc_dummy = 0.0;
        if (!cdid_sel_predict_pair_batch(&sm, &rd, &rf, q1, ch_r1, m_fl, ch_pr, ch_ps, ch_pt, buf_w, &a_pc_dummy)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        for (int o = 0; o < m_fl; ++o) {
          int r1 = ch_r1[o];
          if (!cdid_sel_deriv_tlogit_one(&sm, &rd, &rf, q1, r1, ch_pr[o], ch_ps[o], ch_pt[o], &hzt, tlogit_buf, dq)) {
            UNPROTECT(1); vmaxset(vmax); return R_NilValue;
          }
          for (int l = 0; l < dq; ++l) {
            tl_chunk[o + (size_t) l * (size_t) max_buf] = tlogit_buf[l];
          }
        }
      }

      for (int o = 0; o < m_fl; ++o) {
        double q_w = buf_w[o];
        double dy = buf_dy[o];
        double sw = sqrt(q_w);
        if (!R_FINITE(sw) || !R_FINITE(dy * sw)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        for (int j = 0; j < p_slope; ++j) {
          double zj = buf_Z_sl[o + (size_t) j * (size_t) max_buf];
          double dj = buf_D[o + (size_t) j * (size_t) max_buf];
          if (!R_FINITE(dj * sw)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          ch_sumZ[j] += (long double) (q_w * zj);
          double az = fabs(zj);
          if (az > colmax[j]) colmax[j] = az;
        }
        if (mode <= 2) {
          prob_buf[n_obs_tot + o] = ch_pr[o];
          q_all_buf[n_obs_tot + o] = q_w;
        }
        ch_w_mass += (long double) q_w;
        ch_sumY += (long double) (q_w * dy);
        if (dq > 0) {
          x_chunk[o + 0 * (size_t) max_buf] = 1.0;
          for (int j = 0; j < p_slope; ++j) {
            x_chunk[o + (size_t) (1 + j) * (size_t) max_buf] = buf_Z_sl[o + (size_t) j * (size_t) max_buf];
          }
          x_chunk[o + (size_t) (1 + p_slope) * (size_t) max_buf] = dy;
        }
      }

      weight_mass += (double) ch_w_mass;
      n_obs_tot += m_fl;
      sumY += (double) ch_sumY;
      for (int j = 0; j < p_slope; ++j) sumZ[j] += (double) ch_sumZ[j];
      if (dq > 0) {
        BLAS_INT md_b = (BLAS_INT) m_dim;
        BLAS_INT fl_b = (BLAS_INT) m_fl;
        BLAS_INT ld_b = (BLAS_INT) max_buf;
        double one = 1.0, zero = 0.0;
        for (int l = 0; l < dq; ++l) {
          const double *tl_col = tl_chunk + (size_t) l * (size_t) max_buf;
          for (int c = 0; c < m_dim; ++c) {
            const double *xc = x_chunk + (size_t) c * (size_t) max_buf;
            double *wxc = wx_chunk + (size_t) c * (size_t) max_buf;
            for (int r = 0; r < m_fl; ++r) {
              double v = (buf_w[r] * tl_col[r]) * xc[r];
              if (!R_FINITE(v)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
              wxc[r] = v;
            }
          }
          double *Ml = ch_mom + (size_t) l * (size_t) m_dim * (size_t) m_dim;
          F77_CALL(dgemm)("T", "N", &md_b, &md_b, &fl_b, &one, x_chunk, &ld_b, wx_chunk, &ld_b, &zero, Ml, &md_b FCONE FCONE);
        }
        for (size_t u = 0; u < tot_mom; ++u) mom_p[u] += ch_mom[u];
      }

      if (icpt_flag) {
        for (int j = 0; j < p_slope; ++j) {
          const double *dc = buf_D + (size_t) j * (size_t) max_buf;
          ch_anchor[j] = dc[0];
          int cst = 1;
          for (int i = 1; i < m_fl; ++i) if (dc[i] != ch_anchor[j]) { cst = 0; break; }
          ch_const[j] = cst;
        }
      }
      for (int i = 0; i < m_fl; ++i) {
        double swi = sqrt(buf_w[i]);
        if (icpt_flag) aug_work[i] = swi;
        for (int j = 0; j < p_slope; ++j) {
          aug_work[i + (size_t) (j + icpt_flag) * (size_t) m_fl] =
              buf_D[i + (size_t) j * (size_t) max_buf] * swi;
        }
        aug_work[i + (size_t) (p_slope + icpt_flag) * (size_t) m_fl] = buf_dy[i] * swi;
      }
      int r_ch = (m_fl < k_aug) ? m_fl : k_aug;
      if (!cdid_sel_qr_compress_buf(aug_work, m_fl, k_aug, fac_chunk, jpvt, vn)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      int stacked = fac_rows + r_ch;
      if (stacked > max_merge_rows) max_merge_rows = stacked;
      if (icpt_flag) {
        if (acc_n_qr == 0) {
          for (int j = 0; j < p_slope; ++j) { anchor[j] = ch_anchor[j]; constant[j] = ch_const[j]; }
        } else {
          for (int j = 0; j < p_slope; ++j) {
            constant[j] = (constant[j] && ch_const[j] && (anchor[j] == ch_anchor[j])) ? 1 : 0;
          }
        }
      }
      for (int col = 0; col < k_aug; ++col) {
        double *dst_c = stack_work + (size_t) col * (size_t) stacked;
        for (int i = 0; i < fac_rows; ++i) dst_c[i] = fac_cur[i + (size_t) col * (size_t) fac_rows];
        for (int i = 0; i < r_ch; ++i) dst_c[fac_rows + i] = fac_chunk[i + (size_t) col * (size_t) r_ch];
      }
      int r_new = (stacked < k_aug) ? stacked : k_aug;
      if (!cdid_sel_qr_compress_buf(stack_work, stacked, k_aug, fac_merged, jpvt, vn)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      memcpy(fac_cur, fac_merged, (size_t) r_new * (size_t) k_aug * sizeof(double));
      fac_rows = r_new;
      acc_n_qr += m_fl;
      n_comp += 2;
    }
  }

  if (n_obs_tot < 1 || !R_FINITE(weight_mass) || weight_mass <= 0.0 ||
      (weight_mass / rd.N_dbl) <= 0.0 || !R_FINITE(sumY)) {
    UNPROTECT(1); vmaxset(vmax); return R_NilValue;
  }
  for (int j = 0; j < p_slope; ++j) {
    if (!R_FINITE(sumZ[j])) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
  }
  for (size_t idx = 0; idx < (size_t) m_dim * (size_t) m_dim * (size_t) dq; ++idx) {
    if (!R_FINITE(mom_p[idx])) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
  }

  SEXP fac_s = PROTECT(Rf_allocMatrix(REALSXP, fac_rows, k_aug));
  if (fac_rows > 0) memcpy(REAL(fac_s), fac_cur, (size_t) fac_rows * (size_t) k_aug * sizeof(double));
  SEXP anchor_s = PROTECT(Rf_allocVector(REALSXP, p_slope));
  if (p_slope > 0) memcpy(REAL(anchor_s), anchor, (size_t) p_slope * sizeof(double));
  SEXP const_s = PROTECT(Rf_allocVector(LGLSXP, p_slope));
  for (int j = 0; j < p_slope; ++j) LOGICAL(const_s)[j] = constant[j];
  SEXP sumZ_s = PROTECT(Rf_allocVector(REALSXP, p_slope));
  if (p_slope > 0) memcpy(REAL(sumZ_s), sumZ, (size_t) p_slope * sizeof(double));
  SEXP colmax_s = PROTECT(Rf_allocVector(REALSXP, p_slope));
  if (p_slope > 0) memcpy(REAL(colmax_s), colmax, (size_t) p_slope * sizeof(double));
  SEXP scalars_s = PROTECT(Rf_allocVector(REALSXP, 5));
  REAL(scalars_s)[0] = (double) n_obs_tot;
  REAL(scalars_s)[1] = (double) max_merge_rows;
  REAL(scalars_s)[2] = (double) n_comp;
  REAL(scalars_s)[3] = weight_mass;
  REAL(scalars_s)[4] = sumY;
  SEXP clip_s = (mode <= 2) ? PROTECT(cdid_sel_make_clip_summary(NULL, prob_buf, NULL, q_all_buf, n_obs_tot, 0)) : R_NilValue;
  int n_prot = (mode <= 2) ? 8 : 7;

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_s, 0, fac_s);
  SET_VECTOR_ELT(out_s, 1, anchor_s);
  SET_VECTOR_ELT(out_s, 2, const_s);
  SET_VECTOR_ELT(out_s, 3, sumZ_s);
  SET_VECTOR_ELT(out_s, 4, colmax_s);
  SET_VECTOR_ELT(out_s, 5, scalars_s);
  SET_VECTOR_ELT(out_s, 6, mom_s);
  SET_VECTOR_ELT(out_s, 7, clip_s);
  UNPROTECT(n_prot + 1);
  vmaxset(vmax);
  return out_s;
}

/* ---------------------------------------------------------------------------
   8. C_cdid_sel_if_chunk(mode_sexp, backend_sexp, reader_sexp, frame_sexp,
                          ref_sexp, rows_sexp, links_sexp, pre_sexp, prep_sexp)
   Non-detailed (details == FALSE) complete IF chunk evaluator:
     1 = evaluate_cp04_if_chunk (legacy CP04 IPW)
     2 = evaluate_cp04_selection_if_chunk (general selection IPW, static & D60)
     3 = evaluate_dr_separate_if_chunk (conditional or D60 separate DR)
     4 = evaluate_dr_corollary1_if_chunk (pooled DR, static & D60)
     5 = evaluate_or_corollary1_if_chunk & evaluate_or_selection_if_chunk (separate & pooled OR)
     6 = evaluate_cp04_custom_if_chunk (custom selection IPW)
   --------------------------------------------------------------------------- */
typedef struct {
  int n_obs;
  int *obs;      /* 0-based chunk positions in 0..B-1 */
  double *dy;    /* length n_obs */
  double *fw;    /* length n_obs: effective weight * factor (w * v) */
  double *z;     /* n_obs x p_z col-major */
} cdid_sel_hit_t;

SEXP C_cdid_sel_if_chunk(
    SEXP mode_sexp,
    SEXP backend_sexp,
    SEXP reader_sexp,
    SEXP frame_sexp,
    SEXP ref_sexp,
    SEXP rows_sexp,
    SEXP links_sexp,
    SEXP pre_sexp,
    SEXP prep_sexp) {
  int mode = Rf_asInteger(mode_sexp);
  if (mode < 1 || mode > 6 || TYPEOF(backend_sexp) != VECSXP || TYPEOF(pre_sexp) != VECSXP) {
    return R_NilValue;
  }
  if (Rf_isObject(rows_sexp) || Rf_isObject(links_sexp) ||
      (TYPEOF(rows_sexp) != INTSXP && TYPEOF(rows_sexp) != REALSXP) ||
      (TYPEOF(links_sexp) != INTSXP && TYPEOF(links_sexp) != REALSXP)) {
    return R_NilValue;
  }
  int B = (int) XLENGTH(rows_sexp);
  int n_links = (int) XLENGTH(links_sexp);
  if (B < 1 || n_links < 1) return R_NilValue;

  int uses_sel = Rf_asLogical(cdid_sel_elt(backend_sexp, "uses_selection_models")) == 1;
  cdid_sel_model_t sm;
  memset(&sm, 0, sizeof(sm));
  SEXP sel_fit_s = R_NilValue, up_map_s = R_NilValue;
  if (mode == 2 || mode == 4 || mode == 6 || ((mode == 3 || mode == 5) && uses_sel)) {
    sel_fit_s = cdid_sel_elt(backend_sexp, "selection_backend");
    up_map_s = cdid_sel_elt(backend_sexp, "up_map");
    if (!cdid_sel_init_model(sel_fit_s, &sm) || TYPEOF(up_map_s) != VECSXP) return R_NilValue;
  }

  cdid_sel_reader_t rd;
  if (!cdid_sel_init_reader(reader_sexp, &rd, sm.is_date_model)) return R_NilValue;
  cdid_sel_ref_t rf;
  if (!cdid_sel_init_ref(ref_sexp, &rf)) return R_NilValue;

  const double *Za_ptr = NULL;
  int n_Za_rows = 0, n_Za_cols = 0;
  if (!rf.has_ref && mode != 6) {
    SEXP zdim = Rf_getAttrib(frame_sexp, R_DimSymbol);
    if (TYPEOF(frame_sexp) != REALSXP || TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return R_NilValue;
    Za_ptr = REAL(frame_sexp);
    n_Za_rows = INTEGER(zdim)[0];
    n_Za_cols = INTEGER(zdim)[1];
  } else if (rf.has_ref) {
    Za_ptr = rf.x_Z;
    n_Za_rows = rf.n_x_rows;
    n_Za_cols = rf.p_x;
  }

  SEXP all_links = cdid_sel_elt(backend_sexp, "links");
  if (TYPEOF(all_links) != VECSXP) return R_NilValue;
  int total_links = (int) XLENGTH(all_links);

  const void *vmax = vmaxget();
  int *rows1 = (int *) R_alloc((size_t) B, sizeof(int));
  int *ci_chunk = (int *) R_alloc((size_t) B, sizeof(int));
  for (int i = 0; i < B; ++i) {
    int r1 = cdid_sel_int_at(rows_sexp, i);
    if (r1 == NA_INTEGER || r1 < 1 || r1 > rd.N) { vmaxset(vmax); return R_NilValue; }
    rows1[i] = r1;
    ci_chunk[i] = rd.cohort_index[r1 - 1];
  }
  int *links1 = (int *) R_alloc((size_t) n_links, sizeof(int));
  for (int j = 0; j < n_links; ++j) {
    int lj = cdid_sel_int_at(links_sexp, j);
    if (lj == NA_INTEGER || lj < 1 || lj > total_links) { vmaxset(vmax); return R_NilValue; }
    links1[j] = lj;
  }

  SEXP phi_s = PROTECT(Rf_allocMatrix(REALSXP, B, n_links));
  double *phi = REAL(phi_s);
  memset(phi, 0, (size_t) B * (size_t) n_links * sizeof(double));
  int n_act_cells = 0, n_act_blocks = 0;

  /* =========================================================================
     MODE 1: Legacy CP04 IPW (evaluate_cp04_if_chunk)
     ========================================================================= */
  if (mode == 1) {
    SEXP sids_s = cdid_sel_elt(pre_sexp, "sids");
    SEXP tids_s = cdid_sel_elt(pre_sexp, "tids");
    SEXP sel_specs = cdid_sel_elt(backend_sexp, "selection_specs");
    SEXP sels = cdid_sel_elt(backend_sexp, "selections");
    SEXP trt_specs = cdid_sel_elt(backend_sexp, "treatment_specs");
    SEXP trts = cdid_sel_elt(backend_sexp, "treatments");
    int d = rd.d_Z, p = n_Za_cols;
    if (rd.Z_static == NULL || d < 1 || p < 1 || TYPEOF(sel_specs) != VECSXP || TYPEOF(sels) != VECSXP ||
        TYPEOF(trt_specs) != VECSXP || TYPEOF(trts) != VECSXP) {
      UNPROTECT(1); vmaxset(vmax); return R_NilValue;
    }

    int n_sids = (int) XLENGTH(sids_s);
    int n_tids = (int) XLENGTH(tids_s);
    int total_sels = (int) XLENGTH(sels);
    int total_trts = (int) XLENGTH(trts);

    double *a_if = (double *) R_alloc((size_t) B * (size_t) rd.J, sizeof(double));
    memset(a_if, 0, (size_t) B * (size_t) rd.J * sizeof(double));
    for (int i = 0; i < B; ++i) {
      int r1 = rows1[i];
      int ck = ci_chunk[i];
      if (ck >= 1 && ck <= rd.J && rd.U_logical[r1 - 1] == 1) {
        double Ai = (rd.a_index[r1 - 1] > 0) ? 1.0 : 0.0;
        a_if[i + (size_t) (ck - 1) * (size_t) B] =
          rd.w_static[r1 - 1] * (Ai - rd.a_vec[ck - 1]) / rd.mass_U[ck - 1];
      }
    }

    double **alpha_xi = (double **) R_alloc((size_t) (total_sels + 1), sizeof(double *));
    double **sel_pr = (double **) R_alloc((size_t) (total_sels + 1), sizeof(double *));
    double **sel_S = (double **) R_alloc((size_t) (total_sels + 1), sizeof(double *));
    memset(alpha_xi, 0, (size_t) (total_sels + 1) * sizeof(double *));

    int max_dim = (d > p) ? d : p;
    double *tmp_p = (double *) R_alloc((size_t) (max_dim + 2), sizeof(double));
    int *rr_loc = (int *) R_alloc((size_t) B, sizeof(int));

    for (int j = 0; j < n_sids; ++j) {
      int sj = cdid_sel_int_at(sids_s, j);
      if (sj < 1 || sj > total_sels) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      SEXP s_sp = VECTOR_ELT(sel_specs, sj - 1);
      SEXP f_sel = VECTOR_ELT(sels, sj - 1);
      int ck = Rf_asInteger(cdid_sel_elt(s_sp, "cohort"));
      SEXP pair_s = cdid_sel_elt(s_sp, "pair");
      int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
      int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
      const double *coef = REAL(cdid_sel_elt(f_sel, "coefficients"));
      const double *H_chol = REAL(cdid_sel_elt(f_sel, "H_chol"));

      double *xi = (double *) R_alloc((size_t) B * (size_t) d, sizeof(double));
      double *pr = (double *) R_alloc((size_t) B, sizeof(double));
      double *Sv = (double *) R_alloc((size_t) B, sizeof(double));
      memset(xi, 0, (size_t) B * (size_t) d * sizeof(double));
      memset(pr, 0, (size_t) B * sizeof(double));
      memset(Sv, 0, (size_t) B * sizeof(double));

      int n_rr = 0, hs = -1, ht = -1, hz = -1;
      for (int i = 0; i < B; ++i) {
        int r1 = rows1[i];
        if (ci_chunk[i] == ck && rd.U_logical[r1 - 1] == 1) rr_loc[n_rr++] = i;
      }
      if (n_rr > 0) {
        double *sc_sub = (double *) R_alloc((size_t) n_rr * (size_t) d, sizeof(double));
        double *xi_sub = (double *) R_alloc((size_t) n_rr * (size_t) d, sizeof(double));
        for (int idx = 0; idx < n_rr; ++idx) {
          int i = rr_loc[idx];
          int r1 = rows1[i];
          double pri = 0.0, Si = 0.0, dummy = 0.0;
          if (!cdid_sel_predict_legacy_one(&rd, &rf, 0, coef, d, r1, &hz, &pri)) {
            UNPROTECT(1); vmaxset(vmax); return R_NilValue;
          }
          if (!cdid_sel_pair_one(&rd, r1, ds1, dt1, 0, &hs, &ht, &Si, &dummy)) {
            UNPROTECT(1); vmaxset(vmax); return R_NilValue;
          }
          pr[i] = pri;
          Sv[i] = Si;
          double wf = rd.w_static[r1 - 1] * (Si - pri);
          for (int c = 0; c < d; ++c) {
            sc_sub[idx + (size_t) c * (size_t) n_rr] = wf * rd.Z_static[(r1 - 1) + (size_t) c * (size_t) rd.N];
          }
        }
        if (!cdid_sel_rtr_solve_rows(H_chol, NULL, sc_sub, xi_sub, n_rr, d, 1.0, tmp_p)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        for (int c = 0; c < d; ++c) {
          double *dst = xi + (size_t) c * (size_t) B;
          const double *src = xi_sub + (size_t) c * (size_t) n_rr;
          for (int idx = 0; idx < n_rr; ++idx) dst[rr_loc[idx]] = src[idx];
        }
      }
      alpha_xi[sj] = xi;
      sel_pr[sj] = pr;
      sel_S[sj] = Sv;
    }

    double **mod_phi = (double **) R_alloc((size_t) (total_trts + 1), sizeof(double *));
    memset(mod_phi, 0, (size_t) (total_trts + 1) * sizeof(double *));
    double *rhs_trt = (double *) R_alloc((size_t) B * (size_t) p, sizeof(double));
    double *xb = (double *) R_alloc((size_t) B * (size_t) p, sizeof(double));
    double *dir_diff = (double *) R_alloc((size_t) B, sizeof(double));
    double *S_m = (double *) R_alloc((size_t) B, sizeof(double));
    double *pp_m = (double *) R_alloc((size_t) B, sizeof(double));

    for (int t = 0; t < n_tids; ++t) {
      int tj = cdid_sel_int_at(tids_s, t);
      if (tj < 1 || tj > total_trts) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      SEXP sp = VECTOR_ELT(trt_specs, tj - 1);
      SEXP f_t = VECTOR_ELT(trts, tj - 1);
      SEXP coh_s = cdid_sel_elt(sp, "cohorts");
      SEXP sids_m = cdid_sel_elt(sp, "selection_ids");
      SEXP pair_s = cdid_sel_elt(sp, "pair");
      int m_c = (int) XLENGTH(coh_s);
      int u_len = m_c * (d + 1);
      int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
      int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));

      const double *coef_t = REAL(cdid_sel_elt(f_t, "coefficients"));
      const double *H_t = REAL(cdid_sel_elt(f_t, "H_chol"));
      const double *C_t = REAL(cdid_sel_elt(f_t, "cross"));
      SEXP parts = cdid_sel_elt(f_t, "partials");
      const double *U_t = REAL(cdid_sel_elt(parts, "upstream"));
      const double *V_t = REAL(cdid_sel_elt(parts, "beta"));
      const double *means = REAL(cdid_sel_elt(f_t, "means"));
      if (!cdid_sel_all_finite(coef_t, p) || !cdid_sel_all_finite(H_t, (R_xlen_t) p * p) ||
          !cdid_sel_all_finite(C_t, (R_xlen_t) p * u_len) || !cdid_sel_all_finite(U_t, (R_xlen_t) 4 * u_len) ||
          !cdid_sel_all_finite(V_t, (R_xlen_t) 4 * p) || !cdid_sel_all_finite(means, 4)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      double T_val = means[0], M_val = means[1], dp = means[2], dm = means[3];

      double *up_m = (double *) R_alloc((size_t) B * (size_t) u_len, sizeof(double));
      memset(S_m, 0, (size_t) B * sizeof(double));
      memset(pp_m, 0, (size_t) B * sizeof(double));

      int target_coh = cdid_sel_int_at(coh_s, 0);
      for (int k = 0; k < m_c; ++k) {
        int ck = cdid_sel_int_at(coh_s, k);
        int sj = cdid_sel_int_at(sids_m, k);
        if (sj < 1 || sj > total_sels || alpha_xi[sj] == NULL) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        int base_c = k * (d + 1);
        memcpy(up_m + (size_t) base_c * (size_t) B,
               a_if + (size_t) (ck - 1) * (size_t) B,
               (size_t) B * sizeof(double));
        memcpy(up_m + (size_t) (base_c + 1) * (size_t) B,
               alpha_xi[sj],
               (size_t) B * (size_t) d * sizeof(double));
        const double *Sv = sel_S[sj], *prv = sel_pr[sj];
        for (int i = 0; i < B; ++i) {
          if (ci_chunk[i] == ck) { S_m[i] = Sv[i]; pp_m[i] = prv[i]; }
        }
      }

      memset(rhs_trt, 0, (size_t) B * (size_t) p * sizeof(double));
      memset(dir_diff, 0, (size_t) B * sizeof(double));
      int hs = -1, ht = -1;
      for (int i = 0; i < B; ++i) {
        if (S_m[i] != 1.0) continue;
        int r1 = rows1[i];
        int ai = rd.a_index[r1 - 1];
        if (ai < 1 || ai > n_Za_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        double Si = 0.0, dyi = 0.0;
        if (!cdid_sel_pair_one(&rd, r1, ds1, dt1, 1, &hs, &ht, &Si, &dyi) || Si != 1.0) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        double eta = 0.0;
        for (int c = 0; c < p; ++c) eta += Za_ptr[(ai - 1) + (size_t) c * (size_t) n_Za_rows] * coef_t[c];
        double pr = 0.0, od = 0.0, act = 0.0;
        cdid_sel_clip_logit_scalar(eta, &pr, &od, &act);
        if (!R_FINITE(eta) || !R_FINITE(pr) || !R_FINITE(od)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        double D_ind = (ci_chunk[i] == target_coh) ? 1.0 : 0.0;
        double p_fl = (pp_m[i] > CDID_SEL_CLIP_EPS) ? pp_m[i] : CDID_SEL_CLIP_EPS;
        double fw = rd.w_static[r1 - 1] * rd.a_vec[ci_chunk[i] - 1] / (rd.rho[r1 - 1] * p_fl);
        double wp = fw * D_ind;
        double wm = fw * (1.0 - D_ind) * od;
        double scf = fw * (D_ind - pr);
        for (int c = 0; c < p; ++c) {
          rhs_trt[i + (size_t) c * (size_t) B] = scf * Za_ptr[(ai - 1) + (size_t) c * (size_t) n_Za_rows];
        }
        dir_diff[i] = wp * (dyi - T_val) / dp - wm * (dyi - M_val) / dm;
      }

      /* rhs_trt += up_m %*% t(C_t) where C_t is p x u_len */
      for (int u = 0; u < u_len; ++u) {
        const double *up_u = up_m + (size_t) u * (size_t) B;
        for (int c = 0; c < p; ++c) {
          double cu = C_t[c + (size_t) u * (size_t) p];
          if (cu == 0.0) continue;
          double *rc = rhs_trt + (size_t) c * (size_t) B;
          for (int i = 0; i < B; ++i) rc[i] += up_u[i] * cu;
        }
      }
      if (!cdid_sel_rtr_solve_rows(H_t, NULL, rhs_trt, xb, B, p, 1.0, tmp_p)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }

      double *mphi = (double *) R_alloc((size_t) B, sizeof(double));
      memcpy(mphi, dir_diff, (size_t) B * sizeof(double));
      for (int u = 0; u < u_len; ++u) {
        double du = U_t[0 + 4 * u] - U_t[1 + 4 * u];
        if (du == 0.0) continue;
        const double *up_u = up_m + (size_t) u * (size_t) B;
        for (int i = 0; i < B; ++i) mphi[i] += up_u[i] * du;
      }
      for (int c = 0; c < p; ++c) {
        double dv = U_t ? (V_t[0 + 4 * c] - V_t[1 + 4 * c]) : 0.0;
        if (dv == 0.0) continue;
        const double *xbc = xb + (size_t) c * (size_t) B;
        for (int i = 0; i < B; ++i) mphi[i] += xbc[i] * dv;
      }
      for (int i = 0; i < B; ++i) {
        if (!R_FINITE(mphi[i])) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      }
      mod_phi[tj] = mphi;
    }

    for (int j = 0; j < n_links; ++j) {
      int m_id = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(all_links, links1[j] - 1), "model"));
      if (m_id < 1 || m_id > total_trts || mod_phi[m_id] == NULL) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      memcpy(phi + (size_t) j * (size_t) B, mod_phi[m_id], (size_t) B * sizeof(double));
    }
  }

  /* =========================================================================
     MODE 2: General Selection IPW (evaluate_cp04_selection_if_chunk)
     ========================================================================= */
  else if (mode == 2) {
    SEXP uids_s = cdid_sel_elt(pre_sexp, "uids");
    SEXP tids_s = cdid_sel_elt(pre_sexp, "tids");
    SEXP trt_specs = cdid_sel_elt(backend_sexp, "treatment_specs");
    SEXP trts = cdid_sel_elt(backend_sexp, "treatments");
    int nu = (int) XLENGTH(uids_s);
    int n_tids = (int) XLENGTH(tids_s);
    int total_trts = (int) XLENGTH(trts);
    int D_up = Rf_asInteger(cdid_sel_elt(up_map_s, "D"));
    int p = n_Za_cols;
    if (p < 1 || D_up < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

    int *uids = (nu > 0) ? (int *) R_alloc((size_t) nu, sizeof(int)) : NULL;
    int *uid_col = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
    for (int u = 0; u <= D_up; ++u) uid_col[u] = -1;
    for (int j = 0; j < nu; ++j) {
      int u = cdid_sel_int_at(uids_s, j);
      if (u == NA_INTEGER || u < 1 || u > D_up) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      uids[j] = u;
      uid_col[u] = j;
    }

    double *up = (double *) R_alloc((size_t) B * (size_t) (nu > 0 ? nu : 1), sizeof(double));
    if (!cdid_sel_eval_upstream_internal(&sm, up_map_s, &rd, &rf, rows1, B, uids, nu,
                                         R_NilValue, up, &n_act_cells, &n_act_blocks)) {
      UNPROTECT(1); vmaxset(vmax); return R_NilValue;
    }

    double **mod_phi = (double **) R_alloc((size_t) (total_trts + 1), sizeof(double *));
    memset(mod_phi, 0, (size_t) (total_trts + 1) * sizeof(double *));
    double *rhs_trt = (double *) R_alloc((size_t) B * (size_t) p, sizeof(double));
    double *xb = (double *) R_alloc((size_t) B * (size_t) p, sizeof(double));
    double *dir_diff = (double *) R_alloc((size_t) B, sizeof(double));
    double *tmp_p = (double *) R_alloc((size_t) (p + 2), sizeof(double));
    int *coh_map = (int *) R_alloc((size_t) (rd.J + 1), sizeof(int));

    for (int t = 0; t < n_tids; ++t) {
      int tj = cdid_sel_int_at(tids_s, t);
      if (tj < 1 || tj > total_trts) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      SEXP sp = VECTOR_ELT(trt_specs, tj - 1);
      SEXP f_t = VECTOR_ELT(trts, tj - 1);
      SEXP coh_s = cdid_sel_elt(sp, "cohorts");
      SEXP sids_m = cdid_sel_elt(sp, "selection_ids");
      SEXP pair_s = cdid_sel_elt(sp, "pair");
      SEXP ref_r_s = cdid_sel_elt(sp, "reference");
      int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
      int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
      int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
      int m_c = (int) XLENGTH(coh_s);

      memset(coh_map, 0, (size_t) (rd.J + 1) * sizeof(int));
      for (int k = 0; k < m_c; ++k) {
        int ck = cdid_sel_int_at(coh_s, k);
        if (ck >= 1 && ck <= rd.J) coh_map[ck] = k + 1;
      }

      const double *coef_t = REAL(cdid_sel_elt(f_t, "coefficients"));
      const double *H_t = REAL(cdid_sel_elt(f_t, "H_chol"));
      const double *C_t = REAL(cdid_sel_elt(f_t, "cross"));
      SEXP f_uids_s = cdid_sel_elt(f_t, "upstream_ids");
      int u_len = (int) XLENGTH(f_uids_s);
      const int *f_uids = INTEGER(f_uids_s);
      SEXP parts = cdid_sel_elt(f_t, "partials");
      const double *U_t = REAL(cdid_sel_elt(parts, "upstream"));
      const double *V_t = REAL(cdid_sel_elt(parts, "beta"));
      const double *means = REAL(cdid_sel_elt(f_t, "means"));
      if (!cdid_sel_all_finite(coef_t, p) || !cdid_sel_all_finite(H_t, (R_xlen_t) p * p) ||
          !cdid_sel_all_finite(C_t, (R_xlen_t) p * u_len) || !cdid_sel_all_finite(U_t, (R_xlen_t) 4 * u_len) ||
          !cdid_sel_all_finite(V_t, (R_xlen_t) 4 * p) || !cdid_sel_all_finite(means, 4)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      double T_val = means[0], M_val = means[1], dp = means[2], dm = means[3];

      memset(rhs_trt, 0, (size_t) B * (size_t) p * sizeof(double));
      memset(dir_diff, 0, (size_t) B * sizeof(double));
      int hs = -1, ht = -1, hx = -1, hz1 = -1, hz2 = -1, hw = -1;

      for (int i = 0; i < B; ++i) {
        int ck = ci_chunk[i];
        if (ck < 1 || ck > rd.J || coh_map[ck] == 0) continue;
        int kp = coh_map[ck];
        int r1 = rows1[i];
        double Si = 0.0, dyi = 0.0;
        if (!cdid_sel_pair_one(&rd, r1, ds1, dt1, 1, &hs, &ht, &Si, &dyi)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        if (Si != 1.0) continue;
        int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (rd.a_index[r1 - 1] - 1);
        if (xpos < 0 || xpos >= n_Za_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        double eta = 0.0;
        for (int c = 0; c < p; ++c) eta += Za_ptr[xpos + (size_t) c * (size_t) n_Za_rows] * coef_t[c];
        double pr = 0.0, od = 0.0, act = 0.0;
        cdid_sel_clip_logit_scalar(eta, &pr, &od, &act);
        if (!R_FINITE(eta) || !R_FINITE(pr) || !R_FINITE(od)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

        double pp_i = 0.0, fw = 0.0;
        if (!cdid_sel_predict_pair_one(&sm, &rd, &rf, cdid_sel_int_at(sids_m, kp - 1), r1,
                                       &hz1, &hz2, &hw, &pp_i, NULL, NULL, &fw, NULL)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        double D_ind = (kp == 1) ? 1.0 : 0.0;
        double wp = fw * D_ind;
        double wm = fw * (1.0 - D_ind) * od;
        double scf = fw * (D_ind - pr);
        for (int c = 0; c < p; ++c) {
          rhs_trt[i + (size_t) c * (size_t) B] = scf * Za_ptr[xpos + (size_t) c * (size_t) n_Za_rows];
        }
        dir_diff[i] = wp * (dyi - T_val) / dp - wm * (dyi - M_val) / dm;
      }

      for (int u = 0; u < u_len; ++u) {
        int uc = uid_col[f_uids[u]];
        if (uc < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        const double *up_u = up + (size_t) uc * (size_t) B;
        for (int c = 0; c < p; ++c) {
          double cu = C_t[c + (size_t) u * (size_t) p];
          if (cu == 0.0) continue;
          double *rc = rhs_trt + (size_t) c * (size_t) B;
          for (int i = 0; i < B; ++i) rc[i] += up_u[i] * cu;
        }
      }
      if (!cdid_sel_rtr_solve_rows(H_t, NULL, rhs_trt, xb, B, p, 1.0, tmp_p)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }

      double *mphi = (double *) R_alloc((size_t) B, sizeof(double));
      memcpy(mphi, dir_diff, (size_t) B * sizeof(double));
      for (int u = 0; u < u_len; ++u) {
        double du = U_t[0 + 4 * u] - U_t[1 + 4 * u];
        if (du == 0.0) continue;
        const double *up_u = up + (size_t) uid_col[f_uids[u]] * (size_t) B;
        for (int i = 0; i < B; ++i) mphi[i] += up_u[i] * du;
      }
      for (int c = 0; c < p; ++c) {
        double dv = V_t[0 + 4 * c] - V_t[1 + 4 * c];
        if (dv == 0.0) continue;
        const double *xbc = xb + (size_t) c * (size_t) B;
        for (int i = 0; i < B; ++i) mphi[i] += xbc[i] * dv;
      }
      for (int i = 0; i < B; ++i) {
        if (!R_FINITE(mphi[i])) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      }
      mod_phi[tj] = mphi;
    }

    for (int j = 0; j < n_links; ++j) {
      int m_id = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(all_links, links1[j] - 1), "model"));
      if (m_id < 1 || m_id > total_trts || mod_phi[m_id] == NULL) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      memcpy(phi + (size_t) j * (size_t) B, mod_phi[m_id], (size_t) B * sizeof(double));
    }
  }

  /* =========================================================================
     MODE 3: Separate DR (evaluate_dr_separate_if_chunk: conditional or D60)
     ========================================================================= */
  else if (mode == 3) {
    SEXP sids_s = cdid_sel_elt(pre_sexp, "sids");
    SEXP tids_s = cdid_sel_elt(pre_sexp, "tids");
    SEXP gids_s = cdid_sel_elt(pre_sexp, "gids");
    SEXP uids_s = cdid_sel_elt(pre_sexp, "uids");

    SEXP sup_specs = cdid_sel_elt(backend_sexp, "support_specs");
    SEXP supports = cdid_sel_elt(backend_sexp, "supports");
    SEXP q_of_s_sexp = cdid_sel_elt(backend_sexp, "q_of_s");
    SEXP og = cdid_sel_elt(backend_sexp, "outcome_groups");
    SEXP out_rec = cdid_sel_elt(backend_sexp, "outcome_recipe");
    SEXP trt_rec = cdid_sel_elt(backend_sexp, "treatment_recipe");
    SEXP or_fidx_s = cdid_sel_elt(backend_sexp, "or_frame_index");
    SEXP regs = cdid_sel_elt(backend_sexp, "regressions");
    SEXP reg_specs = cdid_sel_elt(backend_sexp, "regression_specs");
    SEXP trt_specs = cdid_sel_elt(backend_sexp, "treatment_specs");
    SEXP trts = cdid_sel_elt(backend_sexp, "treatments");

    if (TYPEOF(og) != VECSXP || TYPEOF(out_rec) != VECSXP || TYPEOF(trt_rec) != VECSXP ||
        TYPEOF(or_fidx_s) != INTSXP || TYPEOF(sup_specs) != VECSXP || TYPEOF(supports) != VECSXP ||
        TYPEOF(regs) != VECSXP || TYPEOF(reg_specs) != VECSXP || TYPEOF(trt_specs) != VECSXP ||
        TYPEOF(trts) != VECSXP) {
      UNPROTECT(1); vmaxset(vmax); return R_NilValue;
    }

    int icpt = Rf_asLogical(cdid_sel_elt(out_rec, "intercept")) == 1;
    SEXP slope_idx_s = cdid_sel_elt(out_rec, "slope_index");
    SEXP trt_idx_s = cdid_sel_elt(trt_rec, "index");
    if (TYPEOF(slope_idx_s) != INTSXP || TYPEOF(trt_idx_s) != INTSXP) {
      UNPROTECT(1); vmaxset(vmax); return R_NilValue;
    }
    int p_all_or = (int) XLENGTH(slope_idx_s);
    int p_trt = (int) XLENGTH(trt_idx_s);
    const int *trt_cols1 = INTEGER(trt_idx_s);
    const int *or_fidx1 = INTEGER(or_fidx_s);
    const int *slope_idx1 = INTEGER(slope_idx_s);
    int *all_or_cols1 = (p_all_or > 0) ? (int *) R_alloc((size_t) p_all_or, sizeof(int)) : NULL;
    for (int c = 0; c < p_all_or; ++c) {
      int oc = or_fidx1[slope_idx1[c] - 1];
      if (oc < 1 || oc > n_Za_cols) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      all_or_cols1[c] = oc;
    }
    for (int c = 0; c < p_trt; ++c) {
      if (trt_cols1[c] < 1 || trt_cols1[c] > n_Za_cols) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
    }

    int n_sids = (int) XLENGTH(sids_s);
    int n_tids = (int) XLENGTH(tids_s);
    int n_gids = (int) XLENGTH(gids_s);
    int total_sups = (int) XLENGTH(supports);
    int total_trts = (int) XLENGTH(trts);
    int total_regs = (int) XLENGTH(regs);

    /* 1. Upstream evaluation (conditional) or allocation (D60 baseline) */
    int nu = uses_sel ? (int) XLENGTH(uids_s) : n_sids;
    int D_up = uses_sel ? Rf_asInteger(cdid_sel_elt(up_map_s, "D")) : 0;
    int *uid_col = NULL;
    double *up = (double *) R_alloc((size_t) B * (size_t) (nu > 0 ? nu : 1), sizeof(double));
    if (uses_sel) {
      int *uids = (nu > 0) ? (int *) R_alloc((size_t) nu, sizeof(int)) : NULL;
      uid_col = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
      for (int u = 0; u <= D_up; ++u) uid_col[u] = -1;
      for (int j = 0; j < nu; ++j) {
        int u = cdid_sel_int_at(uids_s, j);
        if (u == NA_INTEGER || u < 1 || u > D_up) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        uids[j] = u;
        uid_col[u] = j;
      }
      SEXP prep_up_idx = (TYPEOF(prep_sexp) == VECSXP) ? cdid_sel_elt(prep_sexp, "up_idx") : R_NilValue;
      if (!cdid_sel_eval_upstream_internal(&sm, up_map_s, &rd, &rf, rows1, B, uids, nu,
                                           prep_up_idx, up, &n_act_cells, &n_act_blocks)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
    } else {
      memset(up, 0, (size_t) B * (size_t) (nu > 0 ? nu : 1) * sizeof(double));
    }

    /* 2. Sparse support hits */
    int *sid_to_hpos = (int *) R_alloc((size_t) (total_sups + 1), sizeof(int));
    for (int s = 0; s <= total_sups; ++s) sid_to_hpos[s] = -1;
    cdid_sel_hit_t *hits = (cdid_sel_hit_t *) R_alloc((size_t) (n_sids > 0 ? n_sids : 1), sizeof(cdid_sel_hit_t));
    const int *q_of_s = (uses_sel && TYPEOF(q_of_s_sexp) == INTSXP) ? INTEGER(q_of_s_sexp) : NULL;
    int *tmp_obs = (int *) R_alloc((size_t) B, sizeof(int));
    double *tmp_dy = (double *) R_alloc((size_t) B, sizeof(double));

    for (int j = 0; j < n_sids; ++j) {
      int sid = cdid_sel_int_at(sids_s, j);
      if (sid < 1 || sid > total_sups) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      sid_to_hpos[sid] = j;
      SEXP sp = VECTOR_ELT(sup_specs, sid - 1);
      SEXP u_sup = VECTOR_ELT(supports, sid - 1);
      int ck = Rf_asInteger(cdid_sel_elt(sp, "cohort"));
      SEXP pair_s = cdid_sel_elt(sp, "pair");
      int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
      int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
      SEXP ref_r_s = cdid_sel_elt(sp, "reference");
      int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

      double u_rate = (!uses_sel) ? Rf_asReal(cdid_sel_elt(u_sup, "rate")) : 0.0;
      double u_massA = (!uses_sel) ? Rf_asReal(cdid_sel_elt(u_sup, "massA")) : 1.0;
      double u_v = (!uses_sel) ? Rf_asReal(cdid_sel_elt(u_sup, "v")) : 0.0;
      double a_sup = 0.0;
      if (uses_sel) {
        if (!rf.has_ref) a_sup = rd.a_vec[ck - 1];
        else if (!cdid_sel_ref_stat_one(&rf, r_id, ck, 0, &a_sup)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      }

      int n_obs = 0, hs = -1, ht = -1, hwa = -1;
      double *up_j = (!uses_sel) ? (up + (size_t) j * (size_t) B) : NULL;
      for (int i = 0; i < B; ++i) {
        if (ci_chunk[i] != ck) continue;
        int r1 = rows1[i];
        double Si = 0.0, dyi = 0.0;
        if (!cdid_sel_pair_one(&rd, r1, ds1, dt1, 1, &hs, &ht, &Si, &dyi)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        if (!uses_sel) {
          if (!rf.has_ref) {
            double Ai = (rd.a_index[r1 - 1] > 0) ? 1.0 : 0.0;
            up_j[i] = rd.w_static[r1 - 1] * (Si - Ai * u_rate) / u_massA;
          } else if (rd.a_index[r1 - 1] > 0) {
            double wa = 0.0;
            if (!cdid_sel_ref_w_one(&rf, r_id, r1, &hwa, &wa)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            up_j[i] = wa * (Si - u_rate) / u_massA;
          }
        }
        if (Si == 1.0) {
          tmp_obs[n_obs] = i;
          tmp_dy[n_obs] = dyi;
          ++n_obs;
        }
      }

      hits[j].n_obs = n_obs;
      if (n_obs > 0) {
        hits[j].obs = (int *) R_alloc((size_t) n_obs, sizeof(int));
        hits[j].dy = (double *) R_alloc((size_t) n_obs, sizeof(double));
        hits[j].fw = (double *) R_alloc((size_t) n_obs, sizeof(double));
        hits[j].z = (p_all_or > 0) ? (double *) R_alloc((size_t) n_obs * (size_t) p_all_or, sizeof(double)) : NULL;
        memcpy(hits[j].obs, tmp_obs, (size_t) n_obs * sizeof(int));
        memcpy(hits[j].dy, tmp_dy, (size_t) n_obs * sizeof(double));

        int hx = -1, hz1 = -1, hz2 = -1, hw = -1;
        for (int o = 0; o < n_obs; ++o) {
          int i = tmp_obs[o];
          int r1 = rows1[i];
          int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (rd.a_index[r1 - 1] - 1);
          if (xpos < 0 || xpos >= n_Za_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          for (int c = 0; c < p_all_or; ++c) {
            hits[j].z[o + (size_t) c * (size_t) n_obs] = Za_ptr[xpos + (size_t) (all_or_cols1[c] - 1) * (size_t) n_Za_rows];
          }
          double wi = 0.0;
          if (!rf.has_ref) wi = rd.w_static[r1 - 1];
          else if (!cdid_sel_ref_w_one(&rf, r_id, r1, &hw, &wi)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          double vi = u_v;
          if (uses_sel) {
            double pri = 0.0;
            if (!cdid_sel_predict_pair_one(&sm, &rd, &rf, q_of_s[sid - 1], r1, &hz1, &hz2, &hw,
                                           &pri, NULL, NULL, NULL, NULL)) {
              UNPROTECT(1); vmaxset(vmax); return R_NilValue;
            }
            double p_fl = (pri > CDID_SEL_CLIP_EPS) ? pri : CDID_SEL_CLIP_EPS;
            vi = a_sup / (rd.rho[r1 - 1] * p_fl);
          }
          hits[j].fw[o] = wi * vi;
        }
      } else {
        hits[j].obs = NULL; hits[j].dy = NULL; hits[j].fw = NULL; hits[j].z = NULL;
      }
    }

    /* 3. Group requested treatment models by block gid and map links -> models */
    int *gid_to_kb = (int *) R_alloc((size_t) (total_regs + 1), sizeof(int));
    for (int b = 0; b <= total_regs; ++b) gid_to_kb[b] = -1;
    for (int kb = 0; kb < n_gids; ++kb) {
      int b1 = cdid_sel_int_at(gids_s, kb);
      if (b1 < 1 || b1 > total_regs) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      gid_to_kb[b1] = kb;
    }

    int *kb_head = (int *) R_alloc((size_t) (n_gids > 0 ? n_gids : 1), sizeof(int));
    for (int kb = 0; kb < n_gids; ++kb) kb_head[kb] = -1;
    int *tid_next = (int *) R_alloc((size_t) (n_tids > 0 ? n_tids : 1), sizeof(int));
    for (int t = n_tids - 1; t >= 0; --t) {
      int tid1 = cdid_sel_int_at(tids_s, t);
      if (tid1 < 1 || tid1 > total_trts) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      int reg_b = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(trt_specs, tid1 - 1), "regression"));
      if (reg_b < 1 || reg_b > total_regs || gid_to_kb[reg_b] < 0) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      int kb = gid_to_kb[reg_b];
      tid_next[t] = kb_head[kb];
      kb_head[kb] = t;
    }

    SEXP og_bmap = cdid_sel_elt(og, "block_map");
    SEXP og_groups = cdid_sel_elt(og, "groups");
    double og_N = Rf_asReal(cdid_sel_elt(og, "N"));
    const int *bmap_gid = INTEGER(cdid_sel_elt(og_bmap, "group_id"));
    const int *bmap_lid = INTEGER(cdid_sel_elt(og_bmap, "local_block_id"));

    int max_p = (p_all_or > p_trt) ? p_all_or : p_trt;
    double *tmp_p = (double *) R_alloc((size_t) (max_p + 2), sizeof(double));
    double *hb = (p_all_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_all_or, sizeof(double)) : NULL;
    double *bf_beta = (p_all_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_all_or, sizeof(double)) : NULL;
    double *rhs_trt = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
    double *gamma = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
    double *pd_vec = (double *) R_alloc((size_t) B, sizeof(double));
    double **mod_phi = (double **) R_alloc((size_t) (total_trts + 1), sizeof(double *));
    memset(mod_phi, 0, (size_t) (total_trts + 1) * sizeof(double *));

    for (int kb = 0; kb < n_gids; ++kb) {
      int b1 = cdid_sel_int_at(gids_s, kb);
      int gid1 = bmap_gid[b1 - 1];
      int lid1 = bmap_lid[b1 - 1];
      SEXP g_obj = VECTOR_ELT(og_groups, gid1 - 1);
      SEXP r_loc = VECTOR_ELT(cdid_sel_elt(g_obj, "regressions"), lid1 - 1);
      SEXP keep_s = cdid_sel_elt(g_obj, "keep");
      int of_p = (int) XLENGTH(keep_s);
      const int *of_keep1 = (of_p > 0) ? INTEGER(keep_s) : NULL;
      const double *of_beta = (of_p > 0) ? REAL(cdid_sel_elt(g_obj, "beta")) : NULL;
      double of_L = Rf_asReal(cdid_sel_elt(r_loc, "L"));
      double of_alpha = Rf_asReal(cdid_sel_elt(r_loc, "alpha"));
      const double *of_mean_d = (icpt && of_p > 0) ? REAL(cdid_sel_elt(r_loc, "mean_d")) : NULL;

      /* Compute bf_beta ONLY when of_p > 0 (bf_alpha is unused when details == FALSE!) */
      if (of_p > 0) {
        SEXP rec_b = VECTOR_ELT(regs, b1 - 1);
        SEXP b_sids_s = cdid_sel_elt(VECTOR_ELT(reg_specs, b1 - 1), "support_ids");
        int n_bsids = (int) XLENGTH(b_sids_s);
        const int *b_sids = INTEGER(b_sids_s);
        memset(hb, 0, (size_t) B * (size_t) of_p * sizeof(double));

        for (int si = 0; si < n_bsids; ++si) {
          int hp = sid_to_hpos[b_sids[si]];
          if (hp < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          const cdid_sel_hit_t *h = &hits[hp];
          for (int o = 0; o < h->n_obs; ++o) {
            int i = h->obs[o];
            double dot_db = 0.0;
            for (int c = 0; c < of_p; ++c) {
              double dc = of_L * h->z[o + (size_t) (of_keep1[c] - 1) * (size_t) h->n_obs];
              tmp_p[c] = dc;
              dot_db += dc * of_beta[c];
            }
            double qr = h->fw[o] * (h->dy[o] - of_alpha - dot_db);
            for (int c = 0; c < of_p; ++c) {
              double d_centered = icpt ? (tmp_p[c] - of_mean_d[c]) : tmp_p[c];
              hb[i + (size_t) c * (size_t) B] += qr * d_centered;
            }
          }
        }

        const double *smb = REAL(cdid_sel_elt(rec_b, "score_mean_beta"));
        for (int c = 0; c < of_p; ++c) {
          double smc = smb[c];
          double *hbc = hb + (size_t) c * (size_t) B;
          for (int i = 0; i < B; ++i) hbc[i] -= smc;
        }

        const double *fb = REAL(cdid_sel_elt(rec_b, "forcing_beta"));
        if (uses_sel) {
          SEXP ch_s = cdid_sel_elt(rec_b, "channel_ids");
          int n_ch = (int) XLENGTH(ch_s);
          const int *ch_ids = INTEGER(ch_s);
          for (int u = 0; u < n_ch; ++u) {
            int uc = uid_col[ch_ids[u]];
            if (uc < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            const double *up_u = up + (size_t) uc * (size_t) B;
            for (int c = 0; c < of_p; ++c) {
              double fcu = fb[c + (size_t) u * (size_t) of_p];
              if (fcu == 0.0) continue;
              double *hbc = hb + (size_t) c * (size_t) B;
              for (int i = 0; i < B; ++i) hbc[i] += up_u[i] * fcu;
            }
          }
        } else {
          for (int si = 0; si < n_bsids; ++si) {
            int hp = sid_to_hpos[b_sids[si]];
            const double *up_u = up + (size_t) hp * (size_t) B;
            for (int c = 0; c < of_p; ++c) {
              double fcu = fb[c + (size_t) si * (size_t) of_p];
              if (fcu == 0.0) continue;
              double *hbc = hb + (size_t) c * (size_t) B;
              for (int i = 0; i < B; ++i) hbc[i] += up_u[i] * fcu;
            }
          }
        }

        SEXP fit_obj = cdid_sel_elt(g_obj, "fit");
        SEXP sf_obj = icpt ? cdid_sel_elt(fit_obj, "slope_fit") : fit_obj;
        const double *R_fit = REAL(cdid_sel_elt(sf_obj, "R"));
        const int *piv_fit = INTEGER(cdid_sel_elt(sf_obj, "pivot"));
        if (!cdid_sel_rtr_solve_rows(R_fit, piv_fit, hb, bf_beta, B, of_p, og_N, tmp_p)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
      }

      /* Evaluate every requested treatment model belonging to block b1 */
      for (int t = kb_head[kb]; t >= 0; t = tid_next[t]) {
        int tid1 = cdid_sel_int_at(tids_s, t);
        SEXP sp = VECTOR_ELT(trt_specs, tid1 - 1);
        SEXP f_t = VECTOR_ELT(trts, tid1 - 1);
        SEXP m_sids_s = cdid_sel_elt(sp, "support_ids");
        int n_msids = (int) XLENGTH(m_sids_s);
        const int *m_sids = INTEGER(m_sids_s);
        SEXP ref_r_s = cdid_sel_elt(sp, "reference");
        int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

        const double *coef_t = REAL(cdid_sel_elt(f_t, "coefficients"));
        const double *H_t = REAL(cdid_sel_elt(f_t, "H_chol"));
        const double *C_t = REAL(cdid_sel_elt(f_t, "cross"));
        SEXP parts = cdid_sel_elt(f_t, "partials");
        const double *U_t = REAL(cdid_sel_elt(parts, "upstream"));
        const double *V_t = REAL(cdid_sel_elt(parts, "beta"));
        const double *means = REAL(cdid_sel_elt(f_t, "means"));
        int u_len = 0;
        const int *f_uids = NULL;
        if (uses_sel) {
          SEXP fu_s = cdid_sel_elt(f_t, "upstream_ids");
          u_len = (int) XLENGTH(fu_s);
          f_uids = INTEGER(fu_s);
        } else {
          u_len = n_msids;
        }
        const double *dbd = (of_p > 0) ? REAL(cdid_sel_elt(f_t, "d_bar_diff")) : NULL;
        if (!cdid_sel_all_finite(coef_t, p_trt) || !cdid_sel_all_finite(H_t, (R_xlen_t) p_trt * p_trt) ||
            !cdid_sel_all_finite(C_t, (R_xlen_t) p_trt * u_len) || !cdid_sel_all_finite(U_t, (R_xlen_t) 4 * u_len) ||
            !cdid_sel_all_finite(V_t, (R_xlen_t) 4 * p_trt) || !cdid_sel_all_finite(means, 12) ||
            (of_p > 0 && !cdid_sel_all_finite(dbd, of_p))) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        double T_dr = means[4], M_dr = means[5], dp = means[2], dm = means[3];

        memset(rhs_trt, 0, (size_t) B * (size_t) p_trt * sizeof(double));
        memset(pd_vec, 0, (size_t) B * sizeof(double));
        int hx = -1;

        for (int k = 0; k < n_msids; ++k) {
          int hp = sid_to_hpos[m_sids[k]];
          if (hp < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          const cdid_sel_hit_t *h = &hits[hp];
          double D_ind = (k == 0) ? 1.0 : 0.0;
          for (int o = 0; o < h->n_obs; ++o) {
            int i = h->obs[o];
            int r1 = rows1[i];
            int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (rd.a_index[r1 - 1] - 1);
            if (xpos < 0 || xpos >= n_Za_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

            double eta = 0.0;
            for (int c = 0; c < p_trt; ++c) {
              double xc = Za_ptr[xpos + (size_t) (trt_cols1[c] - 1) * (size_t) n_Za_rows];
              tmp_p[c] = xc;
              eta += xc * coef_t[c];
            }
            double pr = 0.0, od = 0.0, act = 0.0;
            cdid_sel_clip_logit_scalar(eta, &pr, &od, &act);
            if (!R_FINITE(eta) || !R_FINITE(pr) || !R_FINITE(od)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

            double dot_db = 0.0;
            for (int c = 0; c < of_p; ++c) {
              dot_db += (of_L * h->z[o + (size_t) (of_keep1[c] - 1) * (size_t) h->n_obs]) * of_beta[c];
            }
            double resid = h->dy[o] - of_alpha - dot_db;
            double fw = h->fw[o];
            double wp = fw * D_ind;
            double wm = fw * (1.0 - D_ind) * od;
            double scf = fw * (D_ind - pr);
            for (int c = 0; c < p_trt; ++c) {
              rhs_trt[i + (size_t) c * (size_t) B] = scf * tmp_p[c];
            }
            pd_vec[i] = wp * (resid - T_dr) / dp - wm * (resid - M_dr) / dm;
          }
        }

        for (int u = 0; u < u_len; ++u) {
          int uc = uses_sel ? uid_col[f_uids[u]] : sid_to_hpos[m_sids[u]];
          if (uc < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          const double *up_u = up + (size_t) uc * (size_t) B;
          for (int c = 0; c < p_trt; ++c) {
            double cu = C_t[c + (size_t) u * (size_t) p_trt];
            if (cu == 0.0) continue;
            double *rc = rhs_trt + (size_t) c * (size_t) B;
            for (int i = 0; i < B; ++i) rc[i] += up_u[i] * cu;
          }
        }
        if (!cdid_sel_rtr_solve_rows(H_t, NULL, rhs_trt, gamma, B, p_trt, 1.0, tmp_p)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }

        double *mphi = (double *) R_alloc((size_t) B, sizeof(double));
        for (int i = 0; i < B; ++i) {
          double po_s = 0.0, pp_s = 0.0, ps_s = 0.0;
          for (int c = 0; c < of_p; ++c) po_s += bf_beta[i + (size_t) c * (size_t) B] * dbd[c];
          for (int c = 0; c < p_trt; ++c) pp_s += gamma[i + (size_t) c * (size_t) B] * V_t[1 + 4 * c];
          for (int u = 0; u < u_len; ++u) {
            int uc = uses_sel ? uid_col[f_uids[u]] : sid_to_hpos[m_sids[u]];
            ps_s += up[i + (size_t) uc * (size_t) B] * (U_t[0 + 4 * u] - U_t[1 + 4 * u]);
          }
          double po_i = (of_p > 0) ? (-po_s) : 0.0;
          double pp_i = -pp_s;
          double val = ((pd_vec[i] + po_i) + pp_i) + ps_s;
          if (!R_FINITE(val)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          mphi[i] = val;
        }
        mod_phi[tid1] = mphi;
      }
    }

    for (int j = 0; j < n_links; ++j) {
      int m_id = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(all_links, links1[j] - 1), "model"));
      if (m_id < 1 || m_id > total_trts || mod_phi[m_id] == NULL) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      memcpy(phi + (size_t) j * (size_t) B, mod_phi[m_id], (size_t) B * sizeof(double));
    }
  }

  /* =========================================================================
     MODE 4: Pooled DR (evaluate_dr_corollary1_if_chunk)
     ========================================================================= */
  else if (mode == 4) {
    SEXP sids_s = cdid_sel_elt(pre_sexp, "sids");
    SEXP tids_s = cdid_sel_elt(pre_sexp, "tids");
    SEXP bgids_s = cdid_sel_elt(pre_sexp, "beta_gids");
    SEXP uids_s = cdid_sel_elt(pre_sexp, "uids");

    SEXP sh_obj = cdid_sel_elt(backend_sexp, "shared");
    SEXP out_rec = cdid_sel_elt(backend_sexp, "outcome_recipe");
    SEXP trt_rec = cdid_sel_elt(backend_sexp, "treatment_recipe");
    SEXP or_fidx_s = cdid_sel_elt(backend_sexp, "or_frame_index");
    SEXP sup_specs = cdid_sel_elt(backend_sexp, "support_specs");
    SEXP supports = cdid_sel_elt(backend_sexp, "supports");
    SEXP q_of_s_sexp = cdid_sel_elt(backend_sexp, "q_of_s");
    SEXP regs = cdid_sel_elt(backend_sexp, "regressions");
    SEXP reg_specs = cdid_sel_elt(backend_sexp, "regression_specs");
    SEXP trt_specs = cdid_sel_elt(backend_sexp, "treatment_specs");
    SEXP trts = cdid_sel_elt(backend_sexp, "treatments");

    if (TYPEOF(sh_obj) != VECSXP || TYPEOF(out_rec) != VECSXP || TYPEOF(trt_rec) != VECSXP ||
        TYPEOF(or_fidx_s) != INTSXP || TYPEOF(q_of_s_sexp) != INTSXP) {
      UNPROTECT(1); vmaxset(vmax); return R_NilValue;
    }

    int p_or = Rf_asInteger(cdid_sel_elt(sh_obj, "p"));
    int p_trt = Rf_asInteger(cdid_sel_elt(trt_rec, "p"));
    int icpt = Rf_asLogical(cdid_sel_elt(out_rec, "intercept")) == 1;
    const int *trt_cols1 = INTEGER(cdid_sel_elt(trt_rec, "index"));
    int *kept_or_cols1 = (p_or > 0) ? (int *) R_alloc((size_t) p_or, sizeof(int)) : NULL;
    if (p_or > 0) {
      const int *or_fidx1 = INTEGER(or_fidx_s);
      const int *sl_idx1 = INTEGER(cdid_sel_elt(out_rec, "slope_index"));
      const int *keep1 = INTEGER(cdid_sel_elt(sh_obj, "keep"));
      for (int c = 0; c < p_or; ++c) {
        int kc = or_fidx1[sl_idx1[keep1[c] - 1] - 1];
        if (kc < 1 || kc > n_Za_cols) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        kept_or_cols1[c] = kc;
      }
    }
    for (int c = 0; c < p_trt; ++c) {
      if (trt_cols1[c] < 1 || trt_cols1[c] > n_Za_cols) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
    }

    int nu = (int) XLENGTH(uids_s);
    int D_up = Rf_asInteger(cdid_sel_elt(up_map_s, "D"));
    int *uids = (nu > 0) ? (int *) R_alloc((size_t) nu, sizeof(int)) : NULL;
    int *uid_col = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
    for (int u = 0; u <= D_up; ++u) uid_col[u] = -1;
    for (int j = 0; j < nu; ++j) {
      int u = cdid_sel_int_at(uids_s, j);
      if (u == NA_INTEGER || u < 1 || u > D_up) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      uids[j] = u;
      uid_col[u] = j;
    }

    double *up = (double *) R_alloc((size_t) B * (size_t) (nu > 0 ? nu : 1), sizeof(double));
    if (!cdid_sel_eval_upstream_internal(&sm, up_map_s, &rd, &rf, rows1, B, uids, nu,
                                         R_NilValue, up, &n_act_cells, &n_act_blocks)) {
      UNPROTECT(1); vmaxset(vmax); return R_NilValue;
    }

    int max_p = (p_or > p_trt) ? p_or : p_trt;
    double *tmp_p = (double *) R_alloc((size_t) (max_p + 2), sizeof(double));
    double *xi_beta = (p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL;
    const double *sh_beta = (p_or > 0) ? REAL(cdid_sel_elt(sh_obj, "beta")) : NULL;

    /* Shared beta-only OR core when p_or > 0 */
    if (p_or > 0) {
      int n_sids = (int) XLENGTH(sids_s);
      int total_sups = (int) XLENGTH(supports);
      const int *q_of_s = INTEGER(q_of_s_sexp);
      int *sid_to_hpos = (int *) R_alloc((size_t) (total_sups + 1), sizeof(int));
      for (int s = 0; s <= total_sups; ++s) sid_to_hpos[s] = -1;
      cdid_sel_hit_t *hits = (cdid_sel_hit_t *) R_alloc((size_t) (n_sids > 0 ? n_sids : 1), sizeof(cdid_sel_hit_t));
      int *tmp_obs = (int *) R_alloc((size_t) B, sizeof(int));
      double *tmp_dy = (double *) R_alloc((size_t) B, sizeof(double));

      for (int j = 0; j < n_sids; ++j) {
        int sid = cdid_sel_int_at(sids_s, j);
        if (sid < 1 || sid > total_sups) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        sid_to_hpos[sid] = j;
        SEXP sp = VECTOR_ELT(sup_specs, sid - 1);
        int ck = Rf_asInteger(cdid_sel_elt(sp, "cohort"));
        SEXP pair_s = cdid_sel_elt(sp, "pair");
        int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
        int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
        SEXP ref_r_s = cdid_sel_elt(sp, "reference");
        int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

        double a_sup = 0.0;
        if (!rf.has_ref) a_sup = rd.a_vec[ck - 1];
        else if (!cdid_sel_ref_stat_one(&rf, r_id, ck, 0, &a_sup)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

        int n_obs = 0, hs = -1, ht = -1;
        for (int i = 0; i < B; ++i) {
          int r1 = rows1[i];
          if (ci_chunk[i] != ck || rd.U_logical[r1 - 1] != 1) continue;
          double Si = 0.0, dyi = 0.0;
          if (!cdid_sel_pair_one(&rd, r1, ds1, dt1, 1, &hs, &ht, &Si, &dyi)) {
            UNPROTECT(1); vmaxset(vmax); return R_NilValue;
          }
          if (Si == 1.0) { tmp_obs[n_obs] = i; tmp_dy[n_obs] = dyi; ++n_obs; }
        }
        hits[j].n_obs = n_obs;
        if (n_obs > 0) {
          hits[j].obs = (int *) R_alloc((size_t) n_obs, sizeof(int));
          hits[j].dy = (double *) R_alloc((size_t) n_obs, sizeof(double));
          hits[j].fw = (double *) R_alloc((size_t) n_obs, sizeof(double));
          hits[j].z = (double *) R_alloc((size_t) n_obs * (size_t) p_or, sizeof(double));
          memcpy(hits[j].obs, tmp_obs, (size_t) n_obs * sizeof(int));
          memcpy(hits[j].dy, tmp_dy, (size_t) n_obs * sizeof(double));

          int hx = -1, hw = -1;
          int *obs_rows1 = (int *) R_alloc((size_t) n_obs, sizeof(int));
          double *pr_obs = (double *) R_alloc((size_t) n_obs, sizeof(double));
          for (int o = 0; o < n_obs; ++o) obs_rows1[o] = rows1[tmp_obs[o]];
          if (!cdid_sel_predict_pair_batch(&sm, &rd, &rf, q_of_s[sid - 1], obs_rows1, n_obs,
                                           pr_obs, NULL, NULL, NULL, NULL)) {
            UNPROTECT(1); vmaxset(vmax); return R_NilValue;
          }
          for (int o = 0; o < n_obs; ++o) {
            int r1 = obs_rows1[o];
            int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (rd.a_index[r1 - 1] - 1);
            if (xpos < 0 || xpos >= n_Za_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            for (int c = 0; c < p_or; ++c) {
              hits[j].z[o + (size_t) c * (size_t) n_obs] = Za_ptr[xpos + (size_t) (kept_or_cols1[c] - 1) * (size_t) n_Za_rows];
            }
            double wi = 0.0;
            if (!rf.has_ref) wi = rd.w_static[r1 - 1];
            else if (!cdid_sel_ref_w_one(&rf, r_id, r1, &hw, &wi)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            double p_fl = (pr_obs[o] > CDID_SEL_CLIP_EPS) ? pr_obs[o] : CDID_SEL_CLIP_EPS;
            hits[j].fw[o] = wi * (a_sup / (rd.rho[r1 - 1] * p_fl));
          }
        } else {
          hits[j].obs = NULL; hits[j].dy = NULL; hits[j].fw = NULL; hits[j].z = NULL;
        }
      }

      double *h_beta = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
      double *hb_blk = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
      memset(h_beta, 0, (size_t) B * (size_t) p_or * sizeof(double));
      int n_bgids = (int) XLENGTH(bgids_s);
      BLAS_INT por_b = (BLAS_INT) p_or, one_b = 1, B_b = (BLAS_INT) B;
      double one_d = 1.0, zero_d = 0.0;
      for (int k = 0; k < n_bgids; ++k) {
        int b1 = cdid_sel_int_at(bgids_s, k);
        SEXP blk = VECTOR_ELT(regs, b1 - 1);
        double blk_L = Rf_asReal(cdid_sel_elt(blk, "L"));
        double blk_alpha = Rf_asReal(cdid_sel_elt(blk, "alpha"));
        const double *blk_md = icpt ? REAL(cdid_sel_elt(blk, "mean_d")) : NULL;
        SEXP bsids_s = cdid_sel_elt(VECTOR_ELT(reg_specs, b1 - 1), "support_ids");
        int n_bs = (int) XLENGTH(bsids_s);
        const int *bsids = INTEGER(bsids_s);
        memset(hb_blk, 0, (size_t) B * (size_t) p_or * sizeof(double));

        for (int si = 0; si < n_bs; ++si) {
          int hp = sid_to_hpos[bsids[si]];
          if (hp < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          const cdid_sel_hit_t *h = &hits[hp];
          int n_o = h->n_obs;
          if (n_o <= 0) continue;
          const void *vm_h = vmaxget();
          double *d_mat = (double *) R_alloc((size_t) n_o * (size_t) p_or, sizeof(double));
          double *lin_o = (double *) R_alloc((size_t) n_o, sizeof(double));
          for (size_t idx = 0; idx < (size_t) n_o * (size_t) p_or; ++idx) d_mat[idx] = blk_L * h->z[idx];
          BLAS_INT no_b = (BLAS_INT) n_o;
          F77_CALL(dgemm)("N", "N", &no_b, &one_b, &por_b, &one_d,
                          d_mat, &no_b, sh_beta, &por_b, &zero_d, lin_o, &no_b FCONE FCONE);
          for (int o = 0; o < n_o; ++o) {
            int i = h->obs[o];
            double qr = h->fw[o] * (h->dy[o] - blk_alpha - lin_o[o]);
            for (int c = 0; c < p_or; ++c) {
              double dc = d_mat[o + (size_t) c * (size_t) n_o];
              if (icpt) dc -= blk_md[c];
              hb_blk[i + (size_t) c * (size_t) B] += qr * dc;
            }
          }
          vmaxset(vm_h);
        }
        for (size_t idx = 0; idx < (size_t) B * (size_t) p_or; ++idx) h_beta[idx] += hb_blk[idx];
      }

      /* ctrl_uids = sort(unique(unlist(supports[ctrl_sids]$channels))) */
      SEXP ctrl_sids_s = cdid_sel_elt(sh_obj, "control_support_ids");
      int n_ctrl = (int) XLENGTH(ctrl_sids_s);
      const int *ctrl_sids = INTEGER(ctrl_sids_s);
      int *ctrl_mark = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
      memset(ctrl_mark, 0, (size_t) (D_up + 1) * sizeof(int));
      for (int k = 0; k < n_ctrl; ++k) {
        SEXP ch_s = cdid_sel_elt(VECTOR_ELT(supports, ctrl_sids[k] - 1), "channels");
        const int *ch = INTEGER(ch_s);
        for (R_xlen_t c = 0; c < XLENGTH(ch_s); ++c) ctrl_mark[ch[c]] = 1;
      }
      const double *fb_sh = REAL(cdid_sel_elt(sh_obj, "forcing_beta")); /* p_or x D_up */
      int n_cu = 0;
      for (int u = 1; u <= D_up; ++u) if (ctrl_mark[u]) ++n_cu;
      if (n_cu > 0) {
        const void *vm_cu = vmaxget();
        double *up_ctrl = (double *) R_alloc((size_t) B * (size_t) n_cu, sizeof(double));
        double *fb_ctrl = (double *) R_alloc((size_t) p_or * (size_t) n_cu, sizeof(double));
        int cu = 0;
        for (int u = 1; u <= D_up; ++u) {
          if (!ctrl_mark[u]) continue;
          int uc = uid_col[u];
          if (uc < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          memcpy(up_ctrl + (size_t) cu * (size_t) B, up + (size_t) uc * (size_t) B, (size_t) B * sizeof(double));
          memcpy(fb_ctrl + (size_t) cu * (size_t) p_or, fb_sh + (size_t) (u - 1) * (size_t) p_or, (size_t) p_or * sizeof(double));
          ++cu;
        }
        BLAS_INT ncu_b = (BLAS_INT) n_cu;
        F77_CALL(dgemm)("N", "T", &B_b, &por_b, &ncu_b, &one_d,
                        up_ctrl, &B_b, fb_ctrl, &por_b, &zero_d, hb_blk, &B_b FCONE FCONE);
        for (size_t idx = 0; idx < (size_t) B * (size_t) p_or; ++idx) h_beta[idx] += hb_blk[idx];
        vmaxset(vm_cu);
      }

      SEXP fit_sh = cdid_sel_elt(sh_obj, "fit");
      SEXP sf_sh = icpt ? cdid_sel_elt(fit_sh, "slope_fit") : fit_sh;
      const double *R_sh = REAL(cdid_sel_elt(sf_sh, "R"));
      const int *piv_sh = INTEGER(cdid_sel_elt(sf_sh, "pivot"));
      if (!cdid_sel_rtr_solve_rows(R_sh, piv_sh, h_beta, xi_beta, B, p_or, rd.N_dbl, tmp_p)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
    }

    /* Treatment models */
    int n_tids = (int) XLENGTH(tids_s);
    int total_trts = (int) XLENGTH(trts);
    double **mod_phi = (double **) R_alloc((size_t) (total_trts + 1), sizeof(double *));
    memset(mod_phi, 0, (size_t) (total_trts + 1) * sizeof(double *));
    double *rhs_trt = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
    double *xi_gamma = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
    double *pd_vec = (double *) R_alloc((size_t) B, sizeof(double));
    double *po_vec = (double *) R_alloc((size_t) B, sizeof(double));
    double *pp_vec_m = (double *) R_alloc((size_t) B, sizeof(double));
    double *ps_vec_m = (double *) R_alloc((size_t) B, sizeof(double));
    int *rr_idx = (int *) R_alloc((size_t) B, sizeof(int));
    double *dy_rr = (double *) R_alloc((size_t) B, sizeof(double));
    int *coh_map = (int *) R_alloc((size_t) (rd.J + 1), sizeof(int));

    for (int t = 0; t < n_tids; ++t) {
      int tj = cdid_sel_int_at(tids_s, t);
      SEXP sp = VECTOR_ELT(trt_specs, tj - 1);
      SEXP f_t = VECTOR_ELT(trts, tj - 1);
      SEXP reg_obj = VECTOR_ELT(regs, Rf_asInteger(cdid_sel_elt(sp, "regression")) - 1);
      double reg_L = Rf_asReal(cdid_sel_elt(reg_obj, "L"));
      double reg_alpha = Rf_asReal(cdid_sel_elt(reg_obj, "alpha"));

      SEXP coh_s = cdid_sel_elt(sp, "cohorts");
      SEXP sids_m = cdid_sel_elt(sp, "selection_ids");
      SEXP pair_s = cdid_sel_elt(sp, "pair");
      SEXP ref_r_s = cdid_sel_elt(sp, "reference");
      int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;
      int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
      int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
      int m_c = (int) XLENGTH(coh_s);

      memset(coh_map, 0, (size_t) (rd.J + 1) * sizeof(int));
      for (int k = 0; k < m_c; ++k) {
        int ck = cdid_sel_int_at(coh_s, k);
        if (ck >= 1 && ck <= rd.J) coh_map[ck] = k + 1;
      }

      const double *coef_t = REAL(cdid_sel_elt(f_t, "coefficients"));
      const double *H_t = REAL(cdid_sel_elt(f_t, "H_chol"));
      const double *C_t = REAL(cdid_sel_elt(f_t, "cross"));
      SEXP fu_s = cdid_sel_elt(f_t, "upstream_ids");
      int u_len = (int) XLENGTH(fu_s);
      const int *f_uids = INTEGER(fu_s);
      SEXP parts = cdid_sel_elt(f_t, "partials");
      const double *U_t = REAL(cdid_sel_elt(parts, "upstream"));
      const double *V_t = REAL(cdid_sel_elt(parts, "beta"));
      const double *means = REAL(cdid_sel_elt(f_t, "means"));
      const double *dbd = (p_or > 0) ? REAL(cdid_sel_elt(f_t, "d_bar_diff")) : NULL;
      if (!cdid_sel_all_finite(coef_t, p_trt) || !cdid_sel_all_finite(H_t, (R_xlen_t) p_trt * p_trt) ||
          !cdid_sel_all_finite(C_t, (R_xlen_t) p_trt * u_len) || !cdid_sel_all_finite(U_t, (R_xlen_t) 4 * u_len) ||
          !cdid_sel_all_finite(V_t, (R_xlen_t) 4 * p_trt) || !cdid_sel_all_finite(means, 12) ||
          (p_or > 0 && !cdid_sel_all_finite(dbd, p_or))) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      double T_dr = means[4], M_dr = means[5], dp = means[2], dm = means[3];

      memset(rhs_trt, 0, (size_t) B * (size_t) p_trt * sizeof(double));
      memset(pd_vec, 0, (size_t) B * sizeof(double));
      int hs = -1, ht = -1, hx = -1;
      int n_rr = 0;

      for (int i = 0; i < B; ++i) {
        int ck = ci_chunk[i];
        if (ck < 1 || ck > rd.J || coh_map[ck] == 0) continue;
        int r1 = rows1[i];
        double Si = 0.0, dyi = 0.0;
        if (!cdid_sel_pair_one(&rd, r1, ds1, dt1, 1, &hs, &ht, &Si, &dyi)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        if (Si != 1.0) continue;
        rr_idx[n_rr] = i;
        dy_rr[n_rr] = dyi;
        ++n_rr;
      }

      const void *vm_t = vmaxget();
      BLAS_INT B_b = (BLAS_INT) B, ptrt_b = (BLAS_INT) p_trt, one_b = 1;
      double one_d = 1.0, zero_d = 0.0;

      if (n_rr > 0) {
        double *X_trt_rr = (double *) R_alloc((size_t) n_rr * (size_t) p_trt, sizeof(double));
        double *eta_rr = (double *) R_alloc((size_t) n_rr, sizeof(double));
        double *dot_or_rr = (p_or > 0) ? (double *) R_alloc((size_t) n_rr, sizeof(double)) : NULL;
        double *D_or_rr = (p_or > 0) ? (double *) R_alloc((size_t) n_rr * (size_t) p_or, sizeof(double)) : NULL;
        double *fw_rr = (double *) R_alloc((size_t) n_rr, sizeof(double));
        int *sub_idx = (int *) R_alloc((size_t) n_rr, sizeof(int));
        int *sub_rows1 = (int *) R_alloc((size_t) n_rr, sizeof(int));
        double *sub_pr = (double *) R_alloc((size_t) n_rr, sizeof(double));
        double *sub_fw = (double *) R_alloc((size_t) n_rr, sizeof(double));

        for (int o = 0; o < n_rr; ++o) {
          int i = rr_idx[o];
          int r1 = rows1[i];
          int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (rd.a_index[r1 - 1] - 1);
          if (xpos < 0 || xpos >= n_Za_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          for (int c = 0; c < p_trt; ++c) {
            double xc = Za_ptr[xpos + (size_t) (trt_cols1[c] - 1) * (size_t) n_Za_rows];
            if (!R_FINITE(xc)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            X_trt_rr[o + (size_t) c * (size_t) n_rr] = xc;
          }
          for (int c = 0; c < p_or; ++c) {
            double zc = Za_ptr[xpos + (size_t) (kept_or_cols1[c] - 1) * (size_t) n_Za_rows];
            if (!R_FINITE(zc)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            D_or_rr[o + (size_t) c * (size_t) n_rr] = reg_L * zc;
          }
        }

        BLAS_INT nrr_b = (BLAS_INT) n_rr;
        F77_CALL(dgemm)("N", "N", &nrr_b, &one_b, &ptrt_b, &one_d,
                        X_trt_rr, &nrr_b, coef_t, &ptrt_b, &zero_d, eta_rr, &nrr_b FCONE FCONE);
        if (p_or > 0) {
          BLAS_INT por_b = (BLAS_INT) p_or;
          F77_CALL(dgemm)("N", "N", &nrr_b, &one_b, &por_b, &one_d,
                          D_or_rr, &nrr_b, sh_beta, &por_b, &zero_d, dot_or_rr, &nrr_b FCONE FCONE);
        }

        for (int kp = 1; kp <= m_c; ++kp) {
          int n_sk = 0;
          for (int o = 0; o < n_rr; ++o) {
            int i = rr_idx[o];
            if (coh_map[ci_chunk[i]] == kp) {
              sub_idx[n_sk] = o;
              sub_rows1[n_sk] = rows1[i];
              ++n_sk;
            }
          }
          if (n_sk > 0) {
            if (!cdid_sel_predict_pair_batch(&sm, &rd, &rf, cdid_sel_int_at(sids_m, kp - 1),
                                             sub_rows1, n_sk, sub_pr, NULL, NULL, sub_fw, NULL)) {
              UNPROTECT(1); vmaxset(vmax); return R_NilValue;
            }
            for (int s = 0; s < n_sk; ++s) fw_rr[sub_idx[s]] = sub_fw[s];
          }
        }

        for (int o = 0; o < n_rr; ++o) {
          int i = rr_idx[o];
          int kp = coh_map[ci_chunk[i]];
          double pr = 0.0, od = 0.0, act = 0.0;
          cdid_sel_clip_logit_scalar(eta_rr[o], &pr, &od, &act);
          if (!R_FINITE(eta_rr[o]) || !R_FINITE(pr) || !R_FINITE(od)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          double m_hat = (icpt ? reg_alpha : 0.0) + ((p_or > 0) ? dot_or_rr[o] : 0.0);
          double r_rr = dy_rr[o] - m_hat;
          double fw = fw_rr[o];
          double D_ind = (kp == 1) ? 1.0 : 0.0;
          double wp = fw * D_ind;
          double wm = fw * (1.0 - D_ind) * od;
          double scf = fw * (D_ind - pr);
          for (int c = 0; c < p_trt; ++c) {
            rhs_trt[i + (size_t) c * (size_t) B] = scf * X_trt_rr[o + (size_t) c * (size_t) n_rr];
          }
          pd_vec[i] = wp * (r_rr - T_dr) / dp - wm * (r_rr - M_dr) / dm;
        }
      }

      double *upstream_m = (u_len > 0) ? (double *) R_alloc((size_t) B * (size_t) u_len, sizeof(double)) : NULL;
      for (int u = 0; u < u_len; ++u) {
        int uc = uid_col[f_uids[u]];
        if (uc < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        memcpy(upstream_m + (size_t) u * (size_t) B, up + (size_t) uc * (size_t) B, (size_t) B * sizeof(double));
      }
      if (u_len > 0) {
        double *up_cross = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
        BLAS_INT ulen_b = (BLAS_INT) u_len;
        F77_CALL(dgemm)("N", "T", &B_b, &ptrt_b, &ulen_b, &one_d,
                        upstream_m, &B_b, C_t, &ptrt_b, &zero_d, up_cross, &B_b FCONE FCONE);
        for (size_t idx = 0; idx < (size_t) B * (size_t) p_trt; ++idx) rhs_trt[idx] += up_cross[idx];
      }
      if (!cdid_sel_rtr_solve_rows(H_t, NULL, rhs_trt, xi_gamma, B, p_trt, 1.0, tmp_p)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }

      if (p_or > 0) {
        BLAS_INT por_b = (BLAS_INT) p_or;
        F77_CALL(dgemm)("N", "N", &B_b, &one_b, &por_b, &one_d,
                        xi_beta, &B_b, dbd, &por_b, &zero_d, po_vec, &B_b FCONE FCONE);
        for (int i = 0; i < B; ++i) po_vec[i] = -po_vec[i];
      } else {
        memset(po_vec, 0, (size_t) B * sizeof(double));
      }
      {
        double *v2_vec = (double *) R_alloc((size_t) p_trt, sizeof(double));
        for (int c = 0; c < p_trt; ++c) v2_vec[c] = V_t[1 + 4 * c];
        F77_CALL(dgemm)("N", "N", &B_b, &one_b, &ptrt_b, &one_d,
                        xi_gamma, &B_b, v2_vec, &ptrt_b, &zero_d, pp_vec_m, &B_b FCONE FCONE);
        for (int i = 0; i < B; ++i) pp_vec_m[i] = -pp_vec_m[i];
      }
      if (u_len > 0) {
        double *udiff = (double *) R_alloc((size_t) u_len, sizeof(double));
        for (int u = 0; u < u_len; ++u) udiff[u] = U_t[0 + 4 * u] - U_t[1 + 4 * u];
        BLAS_INT ulen_b = (BLAS_INT) u_len;
        F77_CALL(dgemm)("N", "N", &B_b, &one_b, &ulen_b, &one_d,
                        upstream_m, &B_b, udiff, &ulen_b, &zero_d, ps_vec_m, &B_b FCONE FCONE);
      } else {
        memset(ps_vec_m, 0, (size_t) B * sizeof(double));
      }
      vmaxset(vm_t);

      double *mphi = (double *) R_alloc((size_t) B, sizeof(double));
      for (int i = 0; i < B; ++i) {
        double val = ((pd_vec[i] + po_vec[i]) + pp_vec_m[i]) + ps_vec_m[i];
        if (!R_FINITE(val)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        mphi[i] = val;
      }
      mod_phi[tj] = mphi;
    }

    for (int j = 0; j < n_links; ++j) {
      int m_id = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(all_links, links1[j] - 1), "model"));
      if (m_id < 1 || m_id > total_trts || mod_phi[m_id] == NULL) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      memcpy(phi + (size_t) j * (size_t) B, mod_phi[m_id], (size_t) B * sizeof(double));
    }
  }

  /* =========================================================================
     MODE 5: Legacy & General Selection OR (evaluate_or_corollary1_if_chunk &
             evaluate_or_selection_if_chunk, both separate & pooled)
     ========================================================================= */
  else if (mode == 5) {
    SEXP sids_s = cdid_sel_elt(pre_sexp, "sids");
    SEXP gids_s = cdid_sel_elt(pre_sexp, "gids");
    SEXP bgids_s = cdid_sel_elt(pre_sexp, "beta_gids");
    SEXP uids_s = cdid_sel_elt(pre_sexp, "uids");

    SEXP op_s = cdid_sel_elt(backend_sexp, "outcome_pooling");
    int separate = (TYPEOF(op_s) == STRSXP && XLENGTH(op_s) == 1 &&
                    strcmp(CHAR(STRING_ELT(op_s, 0)), "separate") == 0) ? 1 : 0;

    SEXP out_rec = cdid_sel_elt(backend_sexp, "outcome_recipe");
    SEXP sup_specs = cdid_sel_elt(backend_sexp, "support_specs");
    SEXP supports = cdid_sel_elt(backend_sexp, "supports");
    SEXP regs = cdid_sel_elt(backend_sexp, "regressions");
    SEXP reg_specs = cdid_sel_elt(backend_sexp, "regression_specs");
    if (TYPEOF(out_rec) != VECSXP || TYPEOF(sup_specs) != VECSXP || TYPEOF(supports) != VECSXP ||
        TYPEOF(regs) != VECSXP || TYPEOF(reg_specs) != VECSXP) {
      UNPROTECT(1); vmaxset(vmax); return R_NilValue;
    }

    int icpt = Rf_asLogical(cdid_sel_elt(out_rec, "intercept")) == 1;
    SEXP sl_idx_s = cdid_sel_elt(out_rec, "slope_index");
    const int *sl_idx1 = INTEGER(sl_idx_s);
    SEXP sh_obj = separate ? R_NilValue : cdid_sel_elt(backend_sexp, "shared");
    int p_kept = 0;
    int *kept_cols1 = NULL;
    SEXP ref_or_idx_s = rf.has_ref ? cdid_sel_elt(ref_sexp, "or_index") : R_NilValue;
    const int *ref_or_idx1 = (rf.has_ref && TYPEOF(ref_or_idx_s) == INTSXP) ? INTEGER(ref_or_idx_s) : NULL;
    if (rf.has_ref && ref_or_idx1 == NULL) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

    if (separate) {
      p_kept = (int) XLENGTH(sl_idx_s);
      kept_cols1 = (p_kept > 0) ? (int *) R_alloc((size_t) p_kept, sizeof(int)) : NULL;
      for (int c = 0; c < p_kept; ++c) {
        int kc = rf.has_ref ? ref_or_idx1[sl_idx1[c] - 1] : sl_idx1[c];
        if (kc < 1 || kc > n_Za_cols) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        kept_cols1[c] = kc;
      }
    } else {
      p_kept = Rf_asInteger(cdid_sel_elt(sh_obj, "p"));
      kept_cols1 = (p_kept > 0) ? (int *) R_alloc((size_t) p_kept, sizeof(int)) : NULL;
      if (p_kept > 0) {
        const int *keep1 = INTEGER(cdid_sel_elt(sh_obj, "keep"));
        for (int c = 0; c < p_kept; ++c) {
          int slc = sl_idx1[keep1[c] - 1];
          int kc = rf.has_ref ? ref_or_idx1[slc - 1] : slc;
          if (kc < 1 || kc > n_Za_cols) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          kept_cols1[c] = kc;
        }
      }
    }

    int nu = (int) XLENGTH(uids_s);
    int D_up = uses_sel ? Rf_asInteger(cdid_sel_elt(up_map_s, "D")) :
                          Rf_asInteger(cdid_sel_elt(cdid_sel_elt(backend_sexp, "upstream"), "D"));
    int *uids = (nu > 0) ? (int *) R_alloc((size_t) nu, sizeof(int)) : NULL;
    int *uid_col = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
    for (int u = 0; u <= D_up; ++u) uid_col[u] = -1;
    for (int j = 0; j < nu; ++j) {
      int u = cdid_sel_int_at(uids_s, j);
      if (u == NA_INTEGER || u < 1 || u > D_up) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      uids[j] = u;
      uid_col[u] = j;
    }

    double *up = (double *) R_alloc((size_t) B * (size_t) (nu > 0 ? nu : 1), sizeof(double));
    if (uses_sel) {
      if (!cdid_sel_eval_upstream_internal(&sm, up_map_s, &rd, &rf, rows1, B, uids, nu,
                                           R_NilValue, up, &n_act_cells, &n_act_blocks)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
    } else {
      memset(up, 0, (size_t) B * (size_t) (nu > 0 ? nu : 1) * sizeof(double));
      SEXP up_info = cdid_sel_elt(backend_sexp, "upstream");
      int n_a = Rf_asInteger(cdid_sel_elt(up_info, "n_a"));
      const int *up_coh = INTEGER(cdid_sel_elt(up_info, "cohort"));
      for (int j = 0; j < nu; ++j) {
        if (uids[j] <= n_a) {
          int ck = up_coh[uids[j] - 1];
          double *up_j = up + (size_t) j * (size_t) B;
          for (int i = 0; i < B; ++i) {
            int r1 = rows1[i];
            if (ci_chunk[i] == ck && rd.U_logical[r1 - 1] == 1) {
              double Ai = (rd.a_index[r1 - 1] > 0) ? 1.0 : 0.0;
              up_j[i] = rd.w_static[r1 - 1] * (Ai - rd.a_vec[ck - 1]) / rd.mass_U[ck - 1];
            }
          }
        }
      }
    }

    /* Sparse support hits (plus legacy selection score solve when !uses_sel) */
    int n_sids = (int) XLENGTH(sids_s);
    int total_sups = (int) XLENGTH(supports);
    int *sid_to_hpos = (int *) R_alloc((size_t) (total_sups + 1), sizeof(int));
    for (int s = 0; s <= total_sups; ++s) sid_to_hpos[s] = -1;
    cdid_sel_hit_t *hits = (cdid_sel_hit_t *) R_alloc((size_t) (n_sids > 0 ? n_sids : 1), sizeof(cdid_sel_hit_t));
    const int *q_of_s = uses_sel ? INTEGER(cdid_sel_elt(backend_sexp, "q_of_s")) : NULL;
    SEXP leg_sels = (!uses_sel) ? cdid_sel_elt(backend_sexp, "selections") : R_NilValue;

    int max_p = (p_kept > rd.d_Z) ? p_kept : rd.d_Z;
    double *tmp_p = (double *) R_alloc((size_t) (max_p + 4), sizeof(double));
    int *ix_buf = (int *) R_alloc((size_t) B, sizeof(int));
    int *tmp_obs = (int *) R_alloc((size_t) B, sizeof(int));
    double *tmp_dy = (double *) R_alloc((size_t) B, sizeof(double));
    double *tmp_pr = (double *) R_alloc((size_t) B, sizeof(double));
    double *sc_leg = (!uses_sel && rd.d_Z > 0) ? (double *) R_alloc((size_t) B * (size_t) rd.d_Z, sizeof(double)) : NULL;
    double *sc_pack = (!uses_sel && rd.d_Z > 0) ? (double *) R_alloc((size_t) B * (size_t) rd.d_Z, sizeof(double)) : NULL;
    double *xi_pack = (!uses_sel && rd.d_Z > 0) ? (double *) R_alloc((size_t) B * (size_t) rd.d_Z, sizeof(double)) : NULL;
    int *coh_cnt = (int *) R_alloc((size_t) (rd.J + 1), sizeof(int));
    memset(coh_cnt, 0, (size_t) (rd.J + 1) * sizeof(int));
    for (int i = 0; i < B; ++i) {
      int ck_i = ci_chunk[i];
      if (ck_i >= 1 && ck_i <= rd.J && rd.U_logical[rows1[i] - 1] == 1) ++coh_cnt[ck_i];
    }

    for (int j = 0; j < n_sids; ++j) {
      int sid = cdid_sel_int_at(sids_s, j);
      if (sid < 1 || sid > total_sups) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      sid_to_hpos[sid] = j;
      SEXP sp = VECTOR_ELT(sup_specs, sid - 1);
      int ck = Rf_asInteger(cdid_sel_elt(sp, "cohort"));
      if (ck < 1 || ck > rd.J || coh_cnt[ck] == 0) {
        hits[j].n_obs = 0;
        hits[j].obs = NULL; hits[j].dy = NULL; hits[j].fw = NULL; hits[j].z = NULL;
        continue;
      }
      SEXP pair_s = cdid_sel_elt(sp, "pair");
      int ds1 = Rf_asInteger(cdid_sel_elt(pair_s, "s_calendar_index"));
      int dt1 = Rf_asInteger(cdid_sel_elt(pair_s, "t_calendar_index"));
      SEXP ref_r_s = cdid_sel_elt(sp, "reference");
      int r_id = (ref_r_s != R_NilValue) ? Rf_asInteger(ref_r_s) : 0;

      double a_sup = 0.0;
      if (!rf.has_ref) a_sup = rd.a_vec[ck - 1];
      else if (!cdid_sel_ref_stat_one(&rf, r_id, ck, 0, &a_sup)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }

      int n_ix = 0, n_obs = 0, hs = -1, ht = -1, hz = -1;
      const double *leg_coef = NULL, *leg_H = NULL;
      if (!uses_sel) {
        SEXP f_sel = VECTOR_ELT(leg_sels, sid - 1);
        leg_coef = REAL(cdid_sel_elt(f_sel, "coefficients"));
        leg_H = REAL(cdid_sel_elt(f_sel, "H_chol"));
      }

      for (int i = 0; i < B; ++i) {
        int r1 = rows1[i];
        if (ci_chunk[i] != ck || rd.U_logical[r1 - 1] != 1) continue;
        double Si = 0.0, dyi = 0.0;
        if (!cdid_sel_pair_one(&rd, r1, ds1, dt1, 1, &hs, &ht, &Si, &dyi)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        if (!uses_sel) {
          double pri = 0.0;
          if (!cdid_sel_predict_legacy_one(&rd, &rf, 0, leg_coef, rd.d_Z, r1, &hz, &pri)) {
            UNPROTECT(1); vmaxset(vmax); return R_NilValue;
          }
          double wf = rd.w_static[r1 - 1] * (Si - pri);
          for (int c = 0; c < rd.d_Z; ++c) {
            sc_leg[n_ix + (size_t) c * (size_t) B] = wf * rd.Z_static[(r1 - 1) + (size_t) c * (size_t) rd.N];
          }
          ix_buf[n_ix++] = i;
          if (Si == 1.0) { tmp_obs[n_obs] = i; tmp_dy[n_obs] = dyi; tmp_pr[n_obs] = pri; ++n_obs; }
        } else if (Si == 1.0) {
          tmp_obs[n_obs] = i; tmp_dy[n_obs] = dyi; ++n_obs;
        }
      }

      if (!uses_sel && n_ix > 0) {
        for (int c = 0; c < rd.d_Z; ++c) {
          memcpy(sc_pack + (size_t) c * (size_t) n_ix, sc_leg + (size_t) c * (size_t) B, (size_t) n_ix * sizeof(double));
        }
        if (!cdid_sel_rtr_solve_rows(leg_H, NULL, sc_pack, xi_pack, n_ix, rd.d_Z, 1.0, tmp_p)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
        SEXP ch_s = cdid_sel_elt(VECTOR_ELT(supports, sid - 1), "channels");
        const int *ch = INTEGER(ch_s);
        for (int c = 0; c < rd.d_Z; ++c) {
          int uc = uid_col[ch[1 + c]];
          if (uc < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          double *dst = up + (size_t) uc * (size_t) B;
          const double *src = xi_pack + (size_t) c * (size_t) n_ix;
          for (int idx = 0; idx < n_ix; ++idx) dst[ix_buf[idx]] = src[idx];
        }
      }

      hits[j].n_obs = n_obs;
      if (n_obs > 0) {
        hits[j].obs = (int *) R_alloc((size_t) n_obs, sizeof(int));
        hits[j].dy = (double *) R_alloc((size_t) n_obs, sizeof(double));
        hits[j].fw = (double *) R_alloc((size_t) n_obs, sizeof(double));
        hits[j].z = (p_kept > 0) ? (double *) R_alloc((size_t) n_obs * (size_t) p_kept, sizeof(double)) : NULL;
        memcpy(hits[j].obs, tmp_obs, (size_t) n_obs * sizeof(int));
        memcpy(hits[j].dy, tmp_dy, (size_t) n_obs * sizeof(double));

        int hx = -1, hz1 = -1, hz2 = -1, hw = -1;
        for (int o = 0; o < n_obs; ++o) {
          int r1 = rows1[tmp_obs[o]];
          int xpos = rf.has_ref ? cdid_sel_ref_x_row_idx(&rf, r_id, r1, &hx) : (rd.a_index[r1 - 1] - 1);
          if (xpos < 0 || xpos >= n_Za_rows) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          for (int c = 0; c < p_kept; ++c) {
            hits[j].z[o + (size_t) c * (size_t) n_obs] = Za_ptr[xpos + (size_t) (kept_cols1[c] - 1) * (size_t) n_Za_rows];
          }
          double wi = 0.0;
          if (!rf.has_ref) wi = rd.w_static[r1 - 1];
          else if (!cdid_sel_ref_w_one(&rf, r_id, r1, &hw, &wi)) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          double pri = (!uses_sel) ? tmp_pr[o] : 0.0;
          if (uses_sel) {
            if (!cdid_sel_predict_pair_one(&sm, &rd, &rf, q_of_s[sid - 1], r1, &hz1, &hz2, &hw,
                                           &pri, NULL, NULL, NULL, NULL)) {
              UNPROTECT(1); vmaxset(vmax); return R_NilValue;
            }
          }
          double p_fl = (pri > CDID_SEL_CLIP_EPS) ? pri : CDID_SEL_CLIP_EPS;
          hits[j].fw[o] = wi * (a_sup / (rd.rho[r1 - 1] * p_fl));
        }
      } else {
        hits[j].obs = NULL; hits[j].dy = NULL; hits[j].fw = NULL; hits[j].z = NULL;
      }
    }

    int n_gids = (int) XLENGTH(gids_s);
    int total_regs = (int) XLENGTH(regs);

    if (separate) {
      /* Separate OR IF core */
      int *gid_to_kb = (int *) R_alloc((size_t) (total_regs + 1), sizeof(int));
      for (int b = 0; b <= total_regs; ++b) gid_to_kb[b] = -1;
      for (int kb = 0; kb < n_gids; ++kb) gid_to_kb[cdid_sel_int_at(gids_s, kb)] = kb;

      int *kb_lhead = (int *) R_alloc((size_t) (n_gids > 0 ? n_gids : 1), sizeof(int));
      for (int kb = 0; kb < n_gids; ++kb) kb_lhead[kb] = -1;
      int *link_next = (int *) R_alloc((size_t) n_links, sizeof(int));
      for (int j = n_links - 1; j >= 0; --j) {
        int rb = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(all_links, links1[j] - 1), "regression"));
        int kb = gid_to_kb[rb];
        link_next[j] = kb_lhead[kb];
        kb_lhead[kb] = j;
      }

      SEXP og = cdid_sel_elt(backend_sexp, "outcome_groups");
      SEXP og_bmap = cdid_sel_elt(og, "block_map");
      SEXP og_groups = cdid_sel_elt(og, "groups");
      double og_N = Rf_asReal(cdid_sel_elt(og, "N"));
      const int *bmap_gid = INTEGER(cdid_sel_elt(og_bmap, "group_id"));
      const int *bmap_lid = INTEGER(cdid_sel_elt(og_bmap, "local_block_id"));

      double *ha = (double *) R_alloc((size_t) B, sizeof(double));
      double *xa = (double *) R_alloc((size_t) B, sizeof(double));
      double *hb = (p_kept > 0) ? (double *) R_alloc((size_t) B * (size_t) p_kept, sizeof(double)) : NULL;
      double *xb = (p_kept > 0) ? (double *) R_alloc((size_t) B * (size_t) p_kept, sizeof(double)) : NULL;

      for (int kb = 0; kb < n_gids; ++kb) {
        int b1 = cdid_sel_int_at(gids_s, kb);
        int gid1 = bmap_gid[b1 - 1];
        int lid1 = bmap_lid[b1 - 1];
        SEXP g_obj = VECTOR_ELT(og_groups, gid1 - 1);
        SEXP r_loc = VECTOR_ELT(cdid_sel_elt(g_obj, "regressions"), lid1 - 1);
        SEXP keep_s = cdid_sel_elt(g_obj, "keep");
        int of_p = (int) XLENGTH(keep_s);
        const int *of_keep1 = (of_p > 0) ? INTEGER(keep_s) : NULL;
        const double *of_beta = (of_p > 0) ? REAL(cdid_sel_elt(g_obj, "beta")) : NULL;
        double of_L = Rf_asReal(cdid_sel_elt(r_loc, "L"));
        double of_alpha = Rf_asReal(cdid_sel_elt(r_loc, "alpha"));
        const double *of_md = (icpt && of_p > 0) ? REAL(cdid_sel_elt(r_loc, "mean_d")) : NULL;

        SEXP rec_b = VECTOR_ELT(regs, b1 - 1);
        SEXP b_sids_s = cdid_sel_elt(VECTOR_ELT(reg_specs, b1 - 1), "support_ids");
        int n_bsids = (int) XLENGTH(b_sids_s);
        const int *b_sids = INTEGER(b_sids_s);

        memset(ha, 0, (size_t) B * sizeof(double));
        if (of_p > 0) memset(hb, 0, (size_t) B * (size_t) of_p * sizeof(double));

        for (int si = 0; si < n_bsids; ++si) {
          const cdid_sel_hit_t *h = &hits[sid_to_hpos[b_sids[si]]];
          for (int o = 0; o < h->n_obs; ++o) {
            int i = h->obs[o];
            double dot_db = 0.0;
            for (int c = 0; c < of_p; ++c) {
              double dc = of_L * h->z[o + (size_t) (of_keep1[c] - 1) * (size_t) h->n_obs];
              tmp_p[c] = dc;
              dot_db += dc * of_beta[c];
            }
            double qr = h->fw[o] * (h->dy[o] - of_alpha - dot_db);
            ha[i] += qr;
            for (int c = 0; c < of_p; ++c) {
              double dc = icpt ? (tmp_p[c] - of_md[c]) : tmp_p[c];
              hb[i + (size_t) c * (size_t) B] += qr * dc;
            }
          }
        }

        SEXP ch_s = cdid_sel_elt(rec_b, "channel_ids");
        int n_ch = (int) XLENGTH(ch_s);
        const int *ch_ids = INTEGER(ch_s);

        if (of_p > 0) {
          const double *smb = REAL(cdid_sel_elt(rec_b, "score_mean_beta"));
          for (int c = 0; c < of_p; ++c) {
            double smc = smb[c];
            double *hbc = hb + (size_t) c * (size_t) B;
            for (int i = 0; i < B; ++i) hbc[i] -= smc;
          }
          const double *fb = REAL(cdid_sel_elt(rec_b, "forcing_beta"));
          for (int u = 0; u < n_ch; ++u) {
            int uc = uid_col[ch_ids[u]];
            if (uc < 0) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
            const double *up_u = up + (size_t) uc * (size_t) B;
            for (int c = 0; c < of_p; ++c) {
              double fcu = fb[c + (size_t) u * (size_t) of_p];
              if (fcu == 0.0) continue;
              double *hbc = hb + (size_t) c * (size_t) B;
              for (int i = 0; i < B; ++i) hbc[i] += up_u[i] * fcu;
            }
          }
          SEXP fit_obj = cdid_sel_elt(g_obj, "fit");
          SEXP sf_obj = icpt ? cdid_sel_elt(fit_obj, "slope_fit") : fit_obj;
          if (!cdid_sel_rtr_solve_rows(REAL(cdid_sel_elt(sf_obj, "R")), INTEGER(cdid_sel_elt(sf_obj, "pivot")),
                                       hb, xb, B, of_p, og_N, tmp_p)) {
            UNPROTECT(1); vmaxset(vmax); return R_NilValue;
          }
        }

        if (icpt) {
          double sma = Rf_asReal(cdid_sel_elt(rec_b, "score_mean_alpha"));
          double inv_mass = 1.0 / Rf_asReal(cdid_sel_elt(rec_b, "mass"));
          const double *fa = REAL(cdid_sel_elt(rec_b, "forcing_alpha"));
          const double *rec_md = (of_p > 0) ? REAL(cdid_sel_elt(rec_b, "mean_d")) : NULL;
          for (int i = 0; i < B; ++i) xa[i] = ha[i] - sma;
          for (int u = 0; u < n_ch; ++u) {
            double fau = fa[u];
            if (fau == 0.0) continue;
            const double *up_u = up + (size_t) uid_col[ch_ids[u]] * (size_t) B;
            for (int i = 0; i < B; ++i) xa[i] += up_u[i] * fau;
          }
          for (int i = 0; i < B; ++i) xa[i] *= inv_mass;
          for (int c = 0; c < of_p; ++c) {
            double mdc = rec_md[c];
            const double *xbc = xb + (size_t) c * (size_t) B;
            for (int i = 0; i < B; ++i) xa[i] -= xbc[i] * mdc;
          }
        } else {
          memset(xa, 0, (size_t) B * sizeof(double));
        }

        for (int j = kb_lhead[kb]; j >= 0; j = link_next[j]) {
          SEXP lk = VECTOR_ELT(all_links, links1[j] - 1);
          int s_id = Rf_asInteger(cdid_sel_elt(lk, "support_id"));
          double l_mu = Rf_asReal(cdid_sel_elt(lk, "mu"));
          double l_d = Rf_asReal(cdid_sel_elt(lk, "d"));
          const double *l_c = (of_p > 0) ? REAL(cdid_sel_elt(lk, "c")) : NULL;
          const double *rl = REAL(cdid_sel_elt(lk, "rate_loading"));
          SEXP l_ch_s = cdid_sel_elt(VECTOR_ELT(supports, s_id - 1), "channels");
          int nl_ch = (int) XLENGTH(l_ch_s);
          const int *l_ch = INTEGER(l_ch_s);

          double *phi_j = phi + (size_t) j * (size_t) B;
          const cdid_sel_hit_t *h = &hits[sid_to_hpos[s_id]];
          for (int o = 0; o < h->n_obs; ++o) {
            int i = h->obs[o];
            double dot_db = 0.0;
            for (int c = 0; c < of_p; ++c) {
              dot_db += (of_L * h->z[o + (size_t) (of_keep1[c] - 1) * (size_t) h->n_obs]) * of_beta[c];
            }
            double r_val = h->dy[o] - of_alpha - dot_db;
            phi_j[i] = h->fw[o] * (r_val - l_mu) / l_d;
          }
          for (int i = 0; i < B; ++i) phi_j[i] -= xa[i];
          for (int c = 0; c < of_p; ++c) {
            double cc = l_c[c];
            if (cc == 0.0) continue;
            const double *xbc = xb + (size_t) c * (size_t) B;
            for (int i = 0; i < B; ++i) phi_j[i] -= xbc[i] * cc;
          }
          for (int u = 0; u < nl_ch; ++u) {
            double rlu = rl[u];
            if (rlu == 0.0) continue;
            const double *up_u = up + (size_t) uid_col[l_ch[u]] * (size_t) B;
            for (int i = 0; i < B; ++i) phi_j[i] += up_u[i] * rlu;
          }
          for (int i = 0; i < B; ++i) {
            if (!R_FINITE(phi_j[i])) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
          }
        }
      }
    } else {
      /* Pooled OR IF core */
      const double *sh_beta = (p_kept > 0) ? REAL(cdid_sel_elt(sh_obj, "beta")) : NULL;
      double *xi_beta = (p_kept > 0) ? (double *) R_alloc((size_t) B * (size_t) p_kept, sizeof(double)) : NULL;
      int *u_mark = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));

      if (p_kept > 0) {
        double *h_beta = (double *) R_alloc((size_t) B * (size_t) p_kept, sizeof(double));
        memset(h_beta, 0, (size_t) B * (size_t) p_kept * sizeof(double));
        int n_bgids = (int) XLENGTH(bgids_s);
        for (int k = 0; k < n_bgids; ++k) {
          int b1 = cdid_sel_int_at(bgids_s, k);
          SEXP blk = VECTOR_ELT(regs, b1 - 1);
          double blk_L = Rf_asReal(cdid_sel_elt(blk, "L"));
          double blk_alpha = Rf_asReal(cdid_sel_elt(blk, "alpha"));
          const double *blk_md = icpt ? REAL(cdid_sel_elt(blk, "mean_d")) : NULL;
          SEXP bsids_s = cdid_sel_elt(VECTOR_ELT(reg_specs, b1 - 1), "support_ids");
          int n_bs = (int) XLENGTH(bsids_s);
          const int *bsids = INTEGER(bsids_s);
          for (int si = 0; si < n_bs; ++si) {
            const cdid_sel_hit_t *h = &hits[sid_to_hpos[bsids[si]]];
            for (int o = 0; o < h->n_obs; ++o) {
              int i = h->obs[o];
              double dot_db = 0.0;
              for (int c = 0; c < p_kept; ++c) {
                double dc = blk_L * h->z[o + (size_t) c * (size_t) h->n_obs];
                tmp_p[c] = dc;
                dot_db += dc * sh_beta[c];
              }
              double qr = h->fw[o] * (h->dy[o] - blk_alpha - dot_db);
              for (int c = 0; c < p_kept; ++c) {
                double dc = icpt ? (tmp_p[c] - blk_md[c]) : tmp_p[c];
                h_beta[i + (size_t) c * (size_t) B] += qr * dc;
              }
            }
          }
        }

        SEXP ctrl_sids_s = cdid_sel_elt(sh_obj, "control_support_ids");
        int n_ctrl = (int) XLENGTH(ctrl_sids_s);
        const int *ctrl_sids = INTEGER(ctrl_sids_s);
        memset(u_mark, 0, (size_t) (D_up + 1) * sizeof(int));
        for (int k = 0; k < n_ctrl; ++k) {
          SEXP ch_s = cdid_sel_elt(VECTOR_ELT(supports, ctrl_sids[k] - 1), "channels");
          const int *ch = INTEGER(ch_s);
          for (R_xlen_t c = 0; c < XLENGTH(ch_s); ++c) u_mark[ch[c]] = 1;
        }
        const double *fb_sh = REAL(cdid_sel_elt(sh_obj, "forcing_beta"));
        for (int u = 1; u <= D_up; ++u) {
          if (!u_mark[u]) continue;
          const double *up_u = up + (size_t) uid_col[u] * (size_t) B;
          for (int c = 0; c < p_kept; ++c) {
            double fcu = fb_sh[c + (size_t) (u - 1) * (size_t) p_kept];
            if (fcu == 0.0) continue;
            double *hbc = h_beta + (size_t) c * (size_t) B;
            for (int i = 0; i < B; ++i) hbc[i] += up_u[i] * fcu;
          }
        }
        SEXP fit_sh = cdid_sel_elt(sh_obj, "fit");
        SEXP sf_sh = icpt ? cdid_sel_elt(fit_sh, "slope_fit") : fit_sh;
        if (!cdid_sel_rtr_solve_rows(REAL(cdid_sel_elt(sf_sh, "R")), INTEGER(cdid_sel_elt(sf_sh, "pivot")),
                                     h_beta, xi_beta, B, p_kept, rd.N_dbl, tmp_p)) {
          UNPROTECT(1); vmaxset(vmax); return R_NilValue;
        }
      }

      int *gid_to_kb = (int *) R_alloc((size_t) (total_regs + 1), sizeof(int));
      for (int b = 0; b <= total_regs; ++b) gid_to_kb[b] = -1;
      double *alpha_mat = icpt ? (double *) R_alloc((size_t) B * (size_t) (n_gids > 0 ? n_gids : 1), sizeof(double)) : NULL;
      if (icpt) {
        memset(alpha_mat, 0, (size_t) B * (size_t) (n_gids > 0 ? n_gids : 1) * sizeof(double));
        for (int kb = 0; kb < n_gids; ++kb) {
          int b1 = cdid_sel_int_at(gids_s, kb);
          gid_to_kb[b1] = kb;
          SEXP blk = VECTOR_ELT(regs, b1 - 1);
          double blk_L = Rf_asReal(cdid_sel_elt(blk, "L"));
          double blk_alpha = Rf_asReal(cdid_sel_elt(blk, "alpha"));
          double inv_mass = 1.0 / Rf_asReal(cdid_sel_elt(blk, "mass"));
          const double *blk_md = (p_kept > 0) ? REAL(cdid_sel_elt(blk, "mean_d")) : NULL;
          const double *fa = REAL(cdid_sel_elt(blk, "forcing_alpha"));
          SEXP bsids_s = cdid_sel_elt(VECTOR_ELT(reg_specs, b1 - 1), "support_ids");
          int n_bs = (int) XLENGTH(bsids_s);
          const int *bsids = INTEGER(bsids_s);

          double *ak = alpha_mat + (size_t) kb * (size_t) B;
          memset(u_mark, 0, (size_t) (D_up + 1) * sizeof(int));
          for (int si = 0; si < n_bs; ++si) {
            int s_id = bsids[si];
            const cdid_sel_hit_t *h = &hits[sid_to_hpos[s_id]];
            for (int o = 0; o < h->n_obs; ++o) {
              double dot_db = 0.0;
              for (int c = 0; c < p_kept; ++c) {
                dot_db += (blk_L * h->z[o + (size_t) c * (size_t) h->n_obs]) * sh_beta[c];
              }
              ak[h->obs[o]] += h->fw[o] * (h->dy[o] - blk_alpha - dot_db);
            }
            SEXP ch_s = cdid_sel_elt(VECTOR_ELT(supports, s_id - 1), "channels");
            const int *ch = INTEGER(ch_s);
            for (R_xlen_t c = 0; c < XLENGTH(ch_s); ++c) u_mark[ch[c]] = 1;
          }
          int u_idx = 0;
          for (int u = 1; u <= D_up; ++u) {
            if (!u_mark[u]) continue;
            double fau = fa[u_idx++];
            if (fau != 0.0) {
              const double *up_u = up + (size_t) uid_col[u] * (size_t) B;
              for (int i = 0; i < B; ++i) ak[i] += up_u[i] * fau;
            }
          }
          for (int i = 0; i < B; ++i) ak[i] *= inv_mass;
          for (int c = 0; c < p_kept; ++c) {
            double mdc = blk_md[c];
            const double *xbc = xi_beta + (size_t) c * (size_t) B;
            for (int i = 0; i < B; ++i) ak[i] -= xbc[i] * mdc;
          }
        }
      }

      for (int j = 0; j < n_links; ++j) {
        SEXP lk = VECTOR_ELT(all_links, links1[j] - 1);
        int b1 = Rf_asInteger(cdid_sel_elt(lk, "regression"));
        int s_id = Rf_asInteger(cdid_sel_elt(lk, "support_id"));
        double l_mu = Rf_asReal(cdid_sel_elt(lk, "mu"));
        double l_d = Rf_asReal(cdid_sel_elt(lk, "d"));
        const double *l_c = (p_kept > 0) ? REAL(cdid_sel_elt(lk, "c")) : NULL;
        const double *rl = REAL(cdid_sel_elt(lk, "rate_loading"));
        SEXP blk = VECTOR_ELT(regs, b1 - 1);
        double blk_L = Rf_asReal(cdid_sel_elt(blk, "L"));
        double blk_alpha = Rf_asReal(cdid_sel_elt(blk, "alpha"));

        double *phi_j = phi + (size_t) j * (size_t) B;
        const cdid_sel_hit_t *h = &hits[sid_to_hpos[s_id]];
        for (int o = 0; o < h->n_obs; ++o) {
          double dot_db = 0.0;
          for (int c = 0; c < p_kept; ++c) {
            dot_db += (blk_L * h->z[o + (size_t) c * (size_t) h->n_obs]) * sh_beta[c];
          }
          phi_j[h->obs[o]] = h->fw[o] * ((h->dy[o] - blk_alpha - dot_db) - l_mu) / l_d;
        }
        if (icpt) {
          const double *ak = alpha_mat + (size_t) gid_to_kb[b1] * (size_t) B;
          for (int i = 0; i < B; ++i) phi_j[i] -= ak[i];
        }
        for (int c = 0; c < p_kept; ++c) {
          double cc = l_c[c];
          if (cc == 0.0) continue;
          const double *xbc = xi_beta + (size_t) c * (size_t) B;
          for (int i = 0; i < B; ++i) phi_j[i] -= xbc[i] * cc;
        }
        /* In or_pooled_if_core, ucol(s) uses sort(unique(supports[[s]]$channels)), whereas rate_loading
           was built in supports[[s]]$channels order. Check if channels are strictly sorted; if not, fallback! */
        SEXP ch_s = cdid_sel_elt(VECTOR_ELT(supports, s_id - 1), "channels");
        int nl_ch = (int) XLENGTH(ch_s);
        const int *l_ch = INTEGER(ch_s);
        memset(u_mark, 0, (size_t) (D_up + 1) * sizeof(int));
        for (int u = 0; u < nl_ch; ++u) u_mark[l_ch[u]] = 1;
        int rl_idx = 0;
        for (int u = 1; u <= D_up; ++u) {
          if (!u_mark[u]) continue;
          double rlu = rl[rl_idx++];
          if (rlu != 0.0) {
            const double *up_u = up + (size_t) uid_col[u] * (size_t) B;
            for (int i = 0; i < B; ++i) phi_j[i] += up_u[i] * rlu;
          }
        }
        for (int i = 0; i < B; ++i) {
          if (!R_FINITE(phi_j[i])) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
        }
      }
    }
  }

  /* =========================================================================
     MODE 6: Custom Selection IPW (evaluate_cp04_custom_if_chunk)
     ========================================================================= */
  else if (mode == 6) {
    SEXP tids_s = cdid_sel_elt(pre_sexp, "tids");
    SEXP uids_s = cdid_sel_elt(pre_sexp, "uids");
    SEXP trts = cdid_sel_elt(backend_sexp, "treatments");
    int nu = (int) XLENGTH(uids_s);
    int n_tids = (int) XLENGTH(tids_s);
    int total_trts = (int) XLENGTH(trts);
    int D_up = Rf_asInteger(cdid_sel_elt(up_map_s, "D"));

    int *uids = (nu > 0) ? (int *) R_alloc((size_t) nu, sizeof(int)) : NULL;
    int *uid_col = (int *) R_alloc((size_t) (D_up + 1), sizeof(int));
    for (int u = 0; u <= D_up; ++u) uid_col[u] = -1;
    for (int j = 0; j < nu; ++j) {
      int u = cdid_sel_int_at(uids_s, j);
      if (u == NA_INTEGER || u < 1 || u > D_up) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      uids[j] = u;
      uid_col[u] = j;
    }

    double *up = (double *) R_alloc((size_t) B * (size_t) (nu > 0 ? nu : 1), sizeof(double));
    if (nu > 0) {
      if (!cdid_sel_eval_upstream_internal(&sm, up_map_s, &rd, &rf, rows1, B, uids, nu,
                                           R_NilValue, up, &n_act_cells, &n_act_blocks)) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
    } else {
      memset(up, 0, (size_t) B * sizeof(double));
    }

    int *row_to_cpos = (int *) R_alloc((size_t) (rd.N + 1), sizeof(int));
    memset(row_to_cpos, 0, (size_t) (rd.N + 1) * sizeof(int));
    for (int i = 0; i < B; ++i) row_to_cpos[rows1[i]] = i + 1;

    double **mod_phi = (double **) R_alloc((size_t) (total_trts + 1), sizeof(double *));
    memset(mod_phi, 0, (size_t) (total_trts + 1) * sizeof(double *));

    for (int t = 0; t < n_tids; ++t) {
      int tj = cdid_sel_int_at(tids_s, t);
      SEXP f_t = VECTOR_ELT(trts, tj - 1);
      SEXP fu_s = cdid_sel_elt(f_t, "upstream_ids");
      int u_len = (int) XLENGTH(fu_s);
      const int *f_uids = INTEGER(fu_s);
      const double *b_vec = REAL(cdid_sel_elt(f_t, "b"));
      const int *obs_rows = INTEGER(cdid_sel_elt(f_t, "obs_rows"));
      const double *psi_obs = REAL(cdid_sel_elt(f_t, "psi_obs"));
      int n_risk = Rf_asInteger(cdid_sel_elt(f_t, "n_risk"));
      double scale_r = rd.N_dbl / (double) n_risk;

      double *mphi = (double *) R_alloc((size_t) B, sizeof(double));
      memset(mphi, 0, (size_t) B * sizeof(double));
      for (int o = 0; o < n_risk; ++o) {
        int r1 = obs_rows[o];
        if (r1 >= 1 && r1 <= rd.N) {
          int cp = row_to_cpos[r1];
          if (cp > 0) mphi[cp - 1] = scale_r * psi_obs[o];
        }
      }
      for (int u = 0; u < u_len; ++u) {
        double bu = b_vec[u];
        if (bu == 0.0) continue;
        const double *up_u = up + (size_t) uid_col[f_uids[u]] * (size_t) B;
        for (int i = 0; i < B; ++i) mphi[i] += up_u[i] * bu;
      }
      for (int i = 0; i < B; ++i) {
        if (!R_FINITE(mphi[i])) { UNPROTECT(1); vmaxset(vmax); return R_NilValue; }
      }
      mod_phi[tj] = mphi;
    }

    for (int j = 0; j < n_links; ++j) {
      int m_id = Rf_asInteger(cdid_sel_elt(VECTOR_ELT(all_links, links1[j] - 1), "model"));
      if (m_id < 1 || m_id > total_trts || mod_phi[m_id] == NULL) {
        UNPROTECT(1); vmaxset(vmax); return R_NilValue;
      }
      memcpy(phi + (size_t) j * (size_t) B, mod_phi[m_id], (size_t) B * sizeof(double));
    }
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_s, 0, phi_s);
  SET_VECTOR_ELT(out_s, 1, Rf_ScalarInteger(n_act_cells));
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarInteger(n_act_blocks));
  UNPROTECT(2);
  vmaxset(vmax);
  return out_s;
}
