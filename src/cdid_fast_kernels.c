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
#define R_NO_REMAP
#include <Rmath.h>
#undef R_NO_REMAP
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CDID_MAX_EXACT_INT 9007199254740991ULL /* 2^53 - 1 */

static inline double cdid_fp_add(double x, double y) {
  volatile double vx = x;
  volatile double vy = y;
  volatile double vz = vx + vy;
  return vz;
}

static inline double cdid_fp_sub(double x, double y) {
  volatile double vx = x;
  volatile double vy = y;
  volatile double vz = vx - vy;
  return vz;
}

static inline double cdid_fp_mul(double x, double y) {
  volatile double vx = x;
  volatile double vy = y;
  volatile double vz = vx * vy;
  return vz;
}

static inline uint64_t cdid_dbl_bits(double x) {
  uint64_t u = 0;
  memcpy(&u, &x, sizeof(uint64_t));
  return u;
}

/* ---------------------------------------------------------------------------
   1. Fast sorted binary lookup: ipw_owner_match_sorted_double
   --------------------------------------------------------------------------- */
SEXP C_cdid_match_sorted_double(
    SEXP x_sexp,
    SEXP sorted_sexp,
    SEXP left_sexp,
    SEXP right_sexp) {
  if (Rf_isObject(x_sexp) || Rf_isObject(sorted_sexp)) return R_NilValue;
  int tx = TYPEOF(x_sexp);
  int ts = TYPEOF(sorted_sexp);
  if ((tx != INTSXP && tx != REALSXP) || (ts != INTSXP && ts != REALSXP)) return R_NilValue;

  R_xlen_t m = XLENGTH(x_sexp);
  R_xlen_t n = XLENGTH(sorted_sexp);
  if (n > (R_xlen_t) INT_MAX) return R_NilValue;

  if ((TYPEOF(left_sexp) != INTSXP && TYPEOF(left_sexp) != REALSXP) || XLENGTH(left_sexp) != 1 ||
      (TYPEOF(right_sexp) != INTSXP && TYPEOF(right_sexp) != REALSXP) || XLENGTH(right_sexp) != 1) {
    return R_NilValue;
  }
  double left_d = (TYPEOF(left_sexp) == INTSXP)
                      ? (INTEGER(left_sexp)[0] == NA_INTEGER ? R_NaN : (double) INTEGER(left_sexp)[0])
                      : REAL(left_sexp)[0];
  double right_d = (TYPEOF(right_sexp) == INTSXP)
                       ? (INTEGER(right_sexp)[0] == NA_INTEGER ? R_NaN : (double) INTEGER(right_sexp)[0])
                       : REAL(right_sexp)[0];
  if (!R_FINITE(left_d) || !R_FINITE(right_d)) return R_NilValue;

  SEXP out_sexp = PROTECT(Rf_allocVector(INTSXP, m));
  int *out = INTEGER(out_sexp);
  if (m == 0 || left_d > right_d) {
    for (R_xlen_t i = 0; i < m; ++i) out[i] = NA_INTEGER;
    UNPROTECT(1);
    return out_sexp;
  }
  if (left_d < 1.0 || right_d > (double) n || left_d != floor(left_d) || right_d != floor(right_d)) {
    UNPROTECT(1);
    return R_NilValue;
  }
  int left_i = (int) left_d;
  int right_i = (int) right_d;

  if (ts == INTSXP) {
    const int *sv = INTEGER(sorted_sexp);
    for (int k = left_i - 1; k < right_i; ++k) {
      if (sv[k] == NA_INTEGER) {
        UNPROTECT(1);
        return R_NilValue;
      }
    }
    if (tx == INTSXP) {
      const int *xv = INTEGER(x_sexp);
      for (R_xlen_t i = 0; i < m; ++i) {
        int target = xv[i];
        if (target == NA_INTEGER) {
          UNPROTECT(1);
          return R_NilValue;
        }
        int lo = left_i, hi = right_i, found = NA_INTEGER;
        while (lo <= hi) {
          int mid = lo + ((hi - lo) >> 1);
          int val = sv[mid - 1];
          if (val == target) {
            found = mid;
            break;
          }
          if (val < target) lo = mid + 1;
          else hi = mid - 1;
        }
        out[i] = found;
      }
    } else {
      const double *xv = REAL(x_sexp);
      for (R_xlen_t i = 0; i < m; ++i) {
        double target = xv[i];
        if (ISNAN(target)) {
          UNPROTECT(1);
          return R_NilValue;
        }
        int lo = left_i, hi = right_i, found = NA_INTEGER;
        while (lo <= hi) {
          int mid = lo + ((hi - lo) >> 1);
          double val = (double) sv[mid - 1];
          if (val == target) {
            found = mid;
            break;
          }
          if (val < target) lo = mid + 1;
          else hi = mid - 1;
        }
        out[i] = found;
      }
    }
  } else {
    const double *sv = REAL(sorted_sexp);
    for (int k = left_i - 1; k < right_i; ++k) {
      if (ISNAN(sv[k])) {
        UNPROTECT(1);
        return R_NilValue;
      }
    }
    const int *x_i = (tx == INTSXP) ? INTEGER(x_sexp) : NULL;
    const double *x_d = (tx == REALSXP) ? REAL(x_sexp) : NULL;
    for (R_xlen_t i = 0; i < m; ++i) {
      double target;
      if (tx == INTSXP) {
        if (x_i[i] == NA_INTEGER) {
          UNPROTECT(1);
          return R_NilValue;
        }
        target = (double) x_i[i];
      } else {
        target = x_d[i];
        if (ISNAN(target)) {
          UNPROTECT(1);
          return R_NilValue;
        }
      }
      int lo = left_i, hi = right_i, found = NA_INTEGER;
      while (lo <= hi) {
        int mid = lo + ((hi - lo) >> 1);
        double val = sv[mid - 1];
        if (val == target) {
          found = mid;
          break;
        }
        if (val < target) lo = mid + 1;
        else hi = mid - 1;
      }
      out[i] = found;
    }
  }

  UNPROTECT(1);
  return out_sexp;
}

/* ---------------------------------------------------------------------------
   2. Bitwise-exact chunked pair support scan:
      cp04_pair_support_scan & ipw_owner_chunked_pair_support
   --------------------------------------------------------------------------- */
SEXP C_cdid_pair_support_scan(
    SEXP post_a_sexp,
    SEXP bounds4_sexp,
    SEXP a_cohort_idx_sexp,
    SEXP w_a_sexp,
    SEXP cohort_map_sexp,
    SEXP dims2_sexp,
    SEXP step_sexp,
    SEXP check_fn_sexp,
    SEXP rho_sexp) {
  if (TYPEOF(post_a_sexp) != INTSXP || Rf_isObject(post_a_sexp) ||
      TYPEOF(a_cohort_idx_sexp) != INTSXP || Rf_isObject(a_cohort_idx_sexp) ||
      TYPEOF(w_a_sexp) != REALSXP || Rf_isObject(w_a_sexp) ||
      TYPEOF(dims2_sexp) != INTSXP || XLENGTH(dims2_sexp) != 2) {
    return R_NilValue;
  }
  int m = INTEGER(dims2_sexp)[0];
  int n_a = INTEGER(dims2_sexp)[1];
  if (m < 0 || n_a < 0 || XLENGTH(a_cohort_idx_sexp) < (R_xlen_t) n_a ||
      XLENGTH(w_a_sexp) < (R_xlen_t) n_a) {
    return R_NilValue;
  }
  double step = Rf_asReal(step_sexp);
  if (!R_FINITE(step) || step < 1.0 || step != floor(step)) return R_NilValue;

  const int *cmap = NULL;
  int cmap_len = 0;
  if (cohort_map_sexp != R_NilValue) {
    if (TYPEOF(cohort_map_sexp) != INTSXP) return R_NilValue;
    cmap = INTEGER(cohort_map_sexp);
    cmap_len = (int) XLENGTH(cohort_map_sexp);
  }

  SEXP pair_count_sexp = PROTECT(Rf_allocVector(INTSXP, m));
  SEXP pair_mass_sexp = PROTECT(Rf_allocVector(REALSXP, m));
  int *pair_count = INTEGER(pair_count_sexp);
  double *pair_mass = REAL(pair_mass_sexp);
  for (int j = 0; j < m; ++j) {
    pair_count[j] = 0;
    pair_mass[j] = 0.0;
  }

  int n_observed = 0;
  double chunks = 0.0;
  double ops = 0.0;

  if (bounds4_sexp == R_NilValue) {
    chunks = (n_a == 0) ? 0.0 : ceil((double) n_a / step);
    ops = (double) n_a;
  } else {
    if (TYPEOF(bounds4_sexp) != INTSXP || XLENGTH(bounds4_sexp) != 4) {
      UNPROTECT(2);
      return R_NilValue;
    }
    const int *b4 = INTEGER(bounds4_sexp);
    R_xlen_t pa_len = XLENGTH(post_a_sexp);
    int s_start = b4[0], s_end = b4[1], t_start = b4[2], t_end = b4[3];
    if (s_start < 1 || s_end < s_start || (R_xlen_t) s_end > pa_len ||
        t_start < 1 || t_end < t_start || (R_xlen_t) t_end > pa_len) {
      UNPROTECT(2);
      return R_NilValue;
    }

    const void *vmax = vmaxget();
    long double *chunk_sum = (m > 0) ? R_allocLD((size_t) m) : NULL;
    int *chunk_cnt = (m > 0) ? (int *) R_alloc((size_t) m, sizeof(int)) : NULL;
    int *active = (m > 0) ? (int *) R_alloc((size_t) m, sizeof(int)) : NULL;
    for (int j = 0; j < m; ++j) chunk_cnt[j] = 0;

    SEXP call_sexp = R_NilValue;
    int has_cb = (check_fn_sexp != R_NilValue && Rf_isFunction(check_fn_sexp) && TYPEOF(rho_sexp) == ENVSXP);
    if (has_cb) {
      call_sexp = PROTECT(Rf_lang1(check_fn_sexp));
    }

    const int *pa = INTEGER(post_a_sexp);
    const int *aci = INTEGER(a_cohort_idx_sexp);
    const double *w_a = REAL(w_a_sexp);

    int ps = s_start - 1;
    int ps_lim = s_end;
    int pt = t_start - 1;
    int pt_lim = t_end;

    double first_d = 1.0;
    while (first_d <= (double) n_a) {
      if (has_cb) {
        Rf_eval(call_sexp, rho_sexp);
      } else {
        R_CheckUserInterrupt();
      }
      int first = (int) first_d;
      double end_cand = first_d + step - 1.0;
      int last = (end_cand < (double) n_a) ? (int) end_cand : n_a;

      while (ps < ps_lim && pa[ps] < first) ++ps;
      while (pt < pt_lim && pa[pt] < first) ++pt;

      int n_act = 0;
      while (ps < ps_lim && pt < pt_lim) {
        int us = pa[ps];
        int ut = pa[pt];
        if (us > last || ut > last) break;
        if (us == ut) {
          ++n_observed;
          int oc = aci[us - 1];
          int j = (cmap != NULL) ? ((oc >= 1 && oc <= cmap_len) ? cmap[oc - 1] : 0) : oc;
          if (j >= 1 && j <= m) {
            int j0 = j - 1;
            if (chunk_cnt[j0] == 0) {
              active[n_act++] = j0;
              chunk_sum[j0] = 0.0L;
            }
            chunk_cnt[j0] += 1;
            chunk_sum[j0] += (long double) w_a[us - 1];
          }
          ++ps;
          ++pt;
        } else if (us < ut) {
          ++ps;
        } else {
          ++pt;
        }
      }
      while (ps < ps_lim && pa[ps] <= last) ++ps;
      while (pt < pt_lim && pa[pt] <= last) ++pt;

      for (int a = 0; a < n_act; ++a) {
        int j0 = active[a];
        pair_count[j0] += chunk_cnt[j0];
        pair_mass[j0] = cdid_fp_add(pair_mass[j0], (double) chunk_sum[j0]);
        chunk_cnt[j0] = 0;
      }

      chunks += 1.0;
      ops += 4.0 * (double) (last - first + 1);
      first_d = (double) last + 1.0;
    }

    if (has_cb) UNPROTECT(1);
    vmaxset(vmax);
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP nms_sexp = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(nms_sexp, 0, Rf_mkChar("pair_count"));
  SET_STRING_ELT(nms_sexp, 1, Rf_mkChar("pair_mass"));
  SET_STRING_ELT(nms_sexp, 2, Rf_mkChar("n_observed"));
  SET_STRING_ELT(nms_sexp, 3, Rf_mkChar("chunks"));
  SET_STRING_ELT(nms_sexp, 4, Rf_mkChar("ops"));
  SET_VECTOR_ELT(out_sexp, 0, pair_count_sexp);
  SET_VECTOR_ELT(out_sexp, 1, pair_mass_sexp);
  SET_VECTOR_ELT(out_sexp, 2, Rf_ScalarInteger(n_observed));
  SET_VECTOR_ELT(out_sexp, 3, Rf_ScalarReal(chunks));
  SET_VECTOR_ELT(out_sexp, 4, Rf_ScalarReal(ops));
  Rf_setAttrib(out_sexp, R_NamesSymbol, nms_sexp);
  UNPROTECT(4);
  return out_sexp;
}

/* ---------------------------------------------------------------------------
   3. Exact sparse occupancy accumulation & CSR row mapping:
      pipeline_occupancy_accumulate & pipeline_metadata_resolve_catalog_mappings
   --------------------------------------------------------------------------- */
typedef struct {
  int64_t k1;
  int k2;
  int idx;
} cdid_sort_pair_t;

static void cdid_mergesort_pairs(cdid_sort_pair_t *arr, cdid_sort_pair_t *tmp, int n) {
  for (int width = 1; width < n; width <<= 1) {
    for (int i = 0; i < n; i += (width << 1)) {
      int left = i;
      int mid = (i + width < n) ? (i + width) : n;
      int right = (i + (width << 1) < n) ? (i + (width << 1)) : n;
      int p = left, q = mid, r = left;
      while (p < mid && q < right) {
        if (arr[p].k1 < arr[q].k1 || (arr[p].k1 == arr[q].k1 && arr[p].k2 <= arr[q].k2)) {
          tmp[r++] = arr[p++];
        } else {
          tmp[r++] = arr[q++];
        }
      }
      while (p < mid) tmp[r++] = arr[p++];
      while (q < right) tmp[r++] = arr[q++];
    }
    for (int i = 0; i < n; ++i) arr[i] = tmp[i];
  }
}

SEXP C_cdid_occupancy_accumulate(
    SEXP a_units_sexp,
    SEXP unit_n_obs_sexp,
    SEXP unit_obs_start_sexp,
    SEXP a_cohort_idx_sexp,
    SEXP cal_idx_sexp,
    SEXP params_sexp) {
  if (TYPEOF(a_units_sexp) != INTSXP ||
      (TYPEOF(unit_n_obs_sexp) != INTSXP && TYPEOF(unit_n_obs_sexp) != REALSXP) ||
      (TYPEOF(unit_obs_start_sexp) != INTSXP && TYPEOF(unit_obs_start_sexp) != REALSXP) ||
      TYPEOF(a_cohort_idx_sexp) != INTSXP ||
      TYPEOF(cal_idx_sexp) != INTSXP ||
      TYPEOF(params_sexp) != VECSXP || XLENGTH(params_sexp) < 8) {
    return R_NilValue;
  }

  int n_a = Rf_asInteger(VECTOR_ELT(params_sexp, 0));
  int n_cal = Rf_asInteger(VECTOR_ELT(params_sexp, 1));
  double max_occ_cells = Rf_asReal(VECTOR_ELT(params_sexp, 2));
  double slice_pairs = Rf_asReal(VECTOR_ELT(params_sexp, 3));
  double acc_cap = Rf_asReal(VECTOR_ELT(params_sexp, 4));
  double grp_cap = Rf_asReal(VECTOR_ELT(params_sexp, 5));
  SEXP check_fn_sexp = VECTOR_ELT(params_sexp, 6);
  SEXP rho_sexp = VECTOR_ELT(params_sexp, 7);

  if (n_a < 1 || n_cal < 2 || XLENGTH(a_units_sexp) != (R_xlen_t) n_a ||
      XLENGTH(a_cohort_idx_sexp) != (R_xlen_t) n_a ||
      !R_FINITE(max_occ_cells) || !R_FINITE(slice_pairs) || slice_pairs < 1.0) {
    return R_NilValue;
  }

  const int *a_units = INTEGER(a_units_sexp);
  const int *aci = INTEGER(a_cohort_idx_sexp);
  const int *cal_idx = INTEGER(cal_idx_sexp);
  R_xlen_t cal_len = XLENGTH(cal_idx_sexp);
  R_xlen_t n_units_tot = XLENGTH(unit_n_obs_sexp);
  if (XLENGTH(unit_obs_start_sexp) != n_units_tot) return R_NilValue;

  const void *vmax = vmaxget();
  int has_cb = (check_fn_sexp != R_NilValue && Rf_isFunction(check_fn_sexp) && TYPEOF(rho_sexp) == ENVSXP);
  SEXP call_sexp = R_NilValue;
  if (has_cb) call_sexp = PROTECT(Rf_lang1(check_fn_sexp));

  cdid_sort_pair_t *sp = (cdid_sort_pair_t *) R_alloc((size_t) n_a, sizeof(cdid_sort_pair_t));
  cdid_sort_pair_t *sp_tmp = (cdid_sort_pair_t *) R_alloc((size_t) n_a, sizeof(cdid_sort_pair_t));
  int *nl = (int *) R_alloc((size_t) n_a, sizeof(int));
  int *ks = (int *) R_alloc((size_t) n_a, sizeof(int));
  int64_t *sts = (int64_t *) R_alloc((size_t) n_a, sizeof(int64_t));

  const int *uno_i = (TYPEOF(unit_n_obs_sexp) == INTSXP) ? INTEGER(unit_n_obs_sexp) : NULL;
  const double *uno_d = (TYPEOF(unit_n_obs_sexp) == REALSXP) ? REAL(unit_n_obs_sexp) : NULL;
  const int *uos_i = (TYPEOF(unit_obs_start_sexp) == INTSXP) ? INTEGER(unit_obs_start_sexp) : NULL;
  const double *uos_d = (TYPEOF(unit_obs_start_sexp) == REALSXP) ? REAL(unit_obs_start_sexp) : NULL;

  int *orig_nl = (int *) R_alloc((size_t) n_a, sizeof(int));
  int64_t *orig_st = (int64_t *) R_alloc((size_t) n_a, sizeof(int64_t));

  for (int i = 0; i < n_a; ++i) {
    int u = a_units[i];
    if (u < 1 || (R_xlen_t) u > n_units_tot) {
      if (has_cb) UNPROTECT(1);
      vmaxset(vmax);
      return R_NilValue;
    }
    int len_u = (uno_i != NULL) ? uno_i[u - 1] : (int) uno_d[u - 1];
    int64_t st_u = (uos_i != NULL) ? (int64_t) uos_i[u - 1] : (int64_t) uos_d[u - 1];
    if (len_u < 2 || len_u > n_cal || st_u < 1 || st_u + (int64_t) len_u - 1 > (int64_t) cal_len) {
      if (has_cb) UNPROTECT(1);
      vmaxset(vmax);
      return R_NilValue;
    }
    orig_nl[i] = len_u;
    orig_st[i] = st_u;
    sp[i].k1 = -(int64_t) len_u;
    sp[i].k2 = aci[i];
    sp[i].idx = i;
  }

  cdid_mergesort_pairs(sp, sp_tmp, n_a);
  for (int i = 0; i < n_a; ++i) {
    int oi = sp[i].idx;
    nl[i] = orig_nl[oi];
    ks[i] = aci[oi];
    sts[i] = orig_st[oi];
  }

  int64_t *grp = (int64_t *) R_alloc((size_t) n_a, sizeof(int64_t));
  grp[0] = 1;
  for (int i = 1; i < n_a; ++i) {
    grp[i] = grp[i - 1] + ((nl[i] != nl[i - 1] || ks[i] != ks[i - 1]) ? 1 : 0);
  }
  int64_t next_id = grp[n_a - 1];
  int max_len = nl[0];

  /* nl is sorted descending: cnt_ge[j - 1] is the prefix length with nl[i] >= j */
  int *cnt_ge = (int *) R_alloc((size_t) max_len, sizeof(int));
  int cur_p = n_a;
  for (int j = max_len; j >= 1; --j) {
    while (cur_p > 0 && nl[cur_p - 1] < j) --cur_p;
    cnt_ge[j - 1] = cur_p;
  }
  /* Fix: reverse logic - at j=1 cur_p is largest, at j=max_len cur_p is smallest! */
  cur_p = 0;
  for (int j = max_len; j >= 1; --j) {
    while (cur_p < n_a && nl[cur_p] >= j) ++cur_p;
    cnt_ge[j - 1] = cur_p;
  }

  double steps = 0.0;
  double visits = 0.0;
  for (int j = 1; j <= max_len; ++j) {
    int p = cnt_ge[j - 1];
    if (p < 2) break;
    if (has_cb) {
      Rf_eval(call_sexp, rho_sexp);
    } else {
      R_CheckUserInterrupt();
    }
    for (int i = 0; i < p; ++i) {
      sp[i].k1 = grp[i];
      sp[i].k2 = cal_idx[(R_xlen_t) (sts[i] - 1 + (int64_t) (j - 1))];
      sp[i].idx = i;
    }
    cdid_mergesort_pairs(sp, sp_tmp, p);
    int64_t cur = next_id + 1;
    grp[sp[0].idx] = cur;
    for (int i = 1; i < p; ++i) {
      if (sp[i].k1 != sp[i - 1].k1 || sp[i].k2 != sp[i - 1].k2) ++cur;
      grp[sp[i].idx] = cur;
    }
    steps += 1.0;
    visits += (double) p;
    if (cur - next_id == (int64_t) p) break;
    next_id = cur;
  }

  for (int i = 0; i < n_a; ++i) {
    sp[i].k1 = grp[i];
    sp[i].k2 = 0;
    sp[i].idx = i;
  }
  cdid_mergesort_pairs(sp, sp_tmp, n_a);

  int n_hist = 1;
  for (int i = 1; i < n_a; ++i) {
    if (sp[i].k1 != sp[i - 1].k1) ++n_hist;
  }

  double *w_hist = (double *) R_alloc((size_t) n_hist, sizeof(double));
  int *h_len = (int *) R_alloc((size_t) n_hist, sizeof(int));
  int *h_k = (int *) R_alloc((size_t) n_hist, sizeof(int));
  int64_t *h_st = (int64_t *) R_alloc((size_t) n_hist, sizeof(int64_t));

  double max_mult = 0.0;
  double p_hist = 0.0;
  double n_rep = 0.0;
  int64_t q_post_64 = 0;

  int h_idx = 0;
  int run_start = 0;
  for (int i = 1; i <= n_a; ++i) {
    if (i == n_a || sp[i].k1 != sp[run_start].k1) {
      int mult = i - run_start;
      int rep = sp[run_start].idx;
      double wm = (double) mult;
      if (wm > max_mult) max_mult = wm;
      w_hist[h_idx] = wm;
      int hl_i = nl[rep];
      h_len[h_idx] = hl_i;
      h_k[h_idx] = ks[rep];
      h_st[h_idx] = sts[rep];
      double hl_d = (double) hl_i;
      double ph_d = (hl_i % 2 == 0) ? ((hl_d / 2.0) * (hl_d - 1.0)) : (hl_d * ((hl_d - 1.0) / 2.0));
      p_hist += ph_d;
      n_rep += wm * ph_d;
      q_post_64 += (int64_t) (hl_i - 1);
      if (p_hist > (double) CDID_MAX_EXACT_INT || n_rep > (double) CDID_MAX_EXACT_INT ||
          q_post_64 > (int64_t) INT_MAX) {
        if (has_cb) UNPROTECT(1);
        vmaxset(vmax);
        return R_NilValue;
      }
      ++h_idx;
      run_start = i;
    }
  }

  double q_post = (double) q_post_64;
  double cap_cells_d = (max_occ_cells < p_hist) ? max_occ_cells : p_hist;
  if (cap_cells_d > (double) INT_MAX) {
    if (has_cb) UNPROTECT(1);
    vmaxset(vmax);
    return R_NilValue;
  }
  int cap_cells = (int) cap_cells_d;
  double grouping_entries = 20.0 * (double) n_a + 2.0 * (double) max_len + 256.0;
  double plan_entries = 4.0 * (double) n_hist + 11.0 * q_post + 4.0 * cap_cells_d +
                        4.0 * (double) n_cal + 14.0 * slice_pairs + 256.0;
  if (plan_entries > acc_cap || grouping_entries > grp_cap) {
    if (has_cb) UNPROTECT(1);
    vmaxset(vmax);
    return R_NilValue;
  }

  int n_post = (int) q_post_64;
  int64_t *raw_from = (int64_t *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(int64_t));
  int *raw_len = (int *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(int));
  double *raw_w = (double *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(double));
  cdid_sort_pair_t *pp = (cdid_sort_pair_t *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(cdid_sort_pair_t));
  cdid_sort_pair_t *pp_tmp = (cdid_sort_pair_t *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(cdid_sort_pair_t));

  int pos = 0;
  for (int h = 0; h < n_hist; ++h) {
    int len_h = h_len[h];
    int64_t st_h = h_st[h];
    double wh = w_hist[h];
    int kh = h_k[h];
    for (int pi = 1; pi < len_h; ++pi) {
      int64_t pf = st_h + (int64_t) pi;
      raw_from[pos] = pf;
      raw_len[pos] = len_h - pi;
      raw_w[pos] = wh;
      pp[pos].k1 = (int64_t) kh;
      pp[pos].k2 = cal_idx[(R_xlen_t) (pf - 2)];
      pp[pos].idx = pos;
      ++pos;
    }
  }

  if (n_post > 1) cdid_mergesort_pairs(pp, pp_tmp, n_post);

  int64_t *p_from = (int64_t *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(int64_t));
  int *p_len = (int *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(int));
  double *p_w = (double *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(double));
  int *pk = (int *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(int));
  int *ps = (int *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(int));
  double *cum_len = (double *) R_alloc((size_t) (n_post > 0 ? n_post : 1), sizeof(double));

  double csum = 0.0;
  for (int i = 0; i < n_post; ++i) {
    int oi = pp[i].idx;
    pk[i] = (int) pp[i].k1;
    ps[i] = pp[i].k2;
    p_from[i] = raw_from[oi];
    p_len[i] = raw_len[oi];
    p_w[i] = raw_w[oi];
    csum += (double) p_len[i];
    cum_len[i] = csum;
  }

  int n_rows = 0;
  if (n_post > 0) {
    n_rows = 1;
    for (int i = 1; i < n_post; ++i) {
      if (pk[i] != pk[i - 1] || ps[i] != ps[i - 1]) ++n_rows;
    }
  }

  double *buf = (double *) R_alloc((size_t) n_cal, sizeof(double));
  for (int t = 0; t < n_cal; ++t) buf[t] = 0.0;

  int *out_k = (int *) R_alloc((size_t) (cap_cells > 0 ? cap_cells : 1), sizeof(int));
  int *out_s = (int *) R_alloc((size_t) (cap_cells > 0 ? cap_cells : 1), sizeof(int));
  int *out_t = (int *) R_alloc((size_t) (cap_cells > 0 ? cap_cells : 1), sizeof(int));
  double *out_c = (double *) R_alloc((size_t) (cap_cells > 0 ? cap_cells : 1), sizeof(double));

  int n_out = 0;
  double n_slices = 0.0;
  double max_slice = 0.0;

  int r_start = 0;
  while (r_start < n_post) {
    int r_end = r_start;
    while (r_end + 1 < n_post && pk[r_end + 1] == pk[r_start] && ps[r_end + 1] == ps[r_start]) {
      ++r_end;
    }
    if (has_cb) {
      Rf_eval(call_sexp, rho_sexp);
    } else {
      R_CheckUserInterrupt();
    }

    int a = r_start;
    int b = r_end;
    for (;;) {
      double base = (a > 0) ? cum_len[a - 1] : 0.0;
      int e;
      if (cum_len[b] - base <= slice_pairs) {
        e = b;
      } else {
        double target = base + slice_pairs;
        int lo = a, hi = b, ans = a - 1;
        while (lo <= hi) {
          int mid = lo + ((hi - lo) >> 1);
          if (cum_len[mid] <= target) {
            ans = mid;
            lo = mid + 1;
          } else {
            hi = mid - 1;
          }
        }
        e = ans;
        if (e < a) {
          if (has_cb) UNPROTECT(1);
          vmaxset(vmax);
          return R_NilValue;
        }
      }

      double n_tt = cum_len[e] - base;
      for (int idx = a; idx <= e; ++idx) {
        int len_i = p_len[idx];
        int64_t f0 = p_from[idx] - 1;
        double wi = p_w[idx];
        for (int off = 0; off < len_i; ++off) {
          int t_val = cal_idx[(R_xlen_t) (f0 + (int64_t) off)];
          buf[t_val - 1] += wi;
        }
      }
      n_slices += 1.0;
      if (n_tt > max_slice) max_slice = n_tt;
      if (e >= b) break;
      a = e + 1;
    }

    int pk_r = pk[r_start];
    int ps_r = ps[r_start];
    for (int t = ps_r; t < n_cal; ++t) {
      if (buf[t] > 0.0) {
        if (n_out >= cap_cells) {
          if (has_cb) UNPROTECT(1);
          vmaxset(vmax);
          return R_NilValue;
        }
        out_k[n_out] = pk_r;
        out_s[n_out] = ps_r;
        out_t[n_out] = t + 1;
        out_c[n_out] = buf[t];
        buf[t] = 0.0;
        ++n_out;
      }
    }
    r_start = r_end + 1;
  }

  int nprot = has_cb ? 1 : 0;
  SEXP ck_sexp = PROTECT(Rf_allocVector(INTSXP, n_out));
  SEXP si_sexp = PROTECT(Rf_allocVector(INTSXP, n_out));
  SEXP ti_sexp = PROTECT(Rf_allocVector(INTSXP, n_out));
  SEXP ct_sexp = PROTECT(Rf_allocVector(REALSXP, n_out));
  nprot += 4;
  if (n_out > 0) {
    memcpy(INTEGER(ck_sexp), out_k, (size_t) n_out * sizeof(int));
    memcpy(INTEGER(si_sexp), out_s, (size_t) n_out * sizeof(int));
    memcpy(INTEGER(ti_sexp), out_t, (size_t) n_out * sizeof(int));
    memcpy(REAL(ct_sexp), out_c, (size_t) n_out * sizeof(double));
  }

  double peak_phase = grouping_entries;
  if (plan_entries > peak_phase) peak_phase = plan_entries;
  double out_phase = 9.0 * (double) n_out + 256.0;
  if (out_phase > peak_phase) peak_phase = out_phase;
  double actual_work = 2.0 * (double) n_a + visits + 2.0 * q_post + p_hist +
                       (double) n_rows + 2.0 * (double) n_out;

  SEXP diag_sexp = PROTECT(Rf_allocVector(VECSXP, 13));
  SEXP dnms_sexp = PROTECT(Rf_allocVector(STRSXP, 13));
  nprot += 2;
  const char *dnames[13] = {
    "n_distinct_histories", "max_history_multiplicity", "n_refinement_positions",
    "refinement_unit_visits", "n_postings", "n_rows", "n_slices", "max_slice_pairs",
    "n_history_pair_visits", "n_pairs_represented", "fill_buffer_cells",
    "peak_phase_entries", "actual_work_ops"
  };
  for (int i = 0; i < 13; ++i) SET_STRING_ELT(dnms_sexp, i, Rf_mkChar(dnames[i]));
  SET_VECTOR_ELT(diag_sexp, 0, Rf_ScalarInteger(n_hist));
  SET_VECTOR_ELT(diag_sexp, 1, Rf_ScalarReal(max_mult));
  SET_VECTOR_ELT(diag_sexp, 2, Rf_ScalarReal(steps));
  SET_VECTOR_ELT(diag_sexp, 3, Rf_ScalarReal(visits));
  SET_VECTOR_ELT(diag_sexp, 4, Rf_ScalarReal(q_post));
  SET_VECTOR_ELT(diag_sexp, 5, Rf_ScalarInteger(n_rows));
  SET_VECTOR_ELT(diag_sexp, 6, Rf_ScalarReal(n_slices));
  SET_VECTOR_ELT(diag_sexp, 7, Rf_ScalarReal(max_slice));
  SET_VECTOR_ELT(diag_sexp, 8, Rf_ScalarReal(p_hist));
  SET_VECTOR_ELT(diag_sexp, 9, Rf_ScalarReal(n_rep));
  SET_VECTOR_ELT(diag_sexp, 10, Rf_ScalarReal(cap_cells_d));
  SET_VECTOR_ELT(diag_sexp, 11, Rf_ScalarReal(peak_phase));
  SET_VECTOR_ELT(diag_sexp, 12, Rf_ScalarReal(actual_work));
  Rf_setAttrib(diag_sexp, R_NamesSymbol, dnms_sexp);

  SEXP res_sexp = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP rnms_sexp = PROTECT(Rf_allocVector(STRSXP, 5));
  nprot += 2;
  SET_STRING_ELT(rnms_sexp, 0, Rf_mkChar("cohort_a_index"));
  SET_STRING_ELT(rnms_sexp, 1, Rf_mkChar("s_index"));
  SET_STRING_ELT(rnms_sexp, 2, Rf_mkChar("t_index"));
  SET_STRING_ELT(rnms_sexp, 3, Rf_mkChar("count"));
  SET_STRING_ELT(rnms_sexp, 4, Rf_mkChar("diagnostics"));
  SET_VECTOR_ELT(res_sexp, 0, ck_sexp);
  SET_VECTOR_ELT(res_sexp, 1, si_sexp);
  SET_VECTOR_ELT(res_sexp, 2, ti_sexp);
  SET_VECTOR_ELT(res_sexp, 3, ct_sexp);
  SET_VECTOR_ELT(res_sexp, 4, diag_sexp);
  Rf_setAttrib(res_sexp, R_NamesSymbol, rnms_sexp);

  UNPROTECT(nprot);
  vmaxset(vmax);
  return res_sexp;
}

SEXP C_cdid_catalog_row_mappings(
    SEXP csr_sexp,
    SEXP pools_sexp,
    SEXP rows_sexp) {
  if (TYPEOF(csr_sexp) != VECSXP || XLENGTH(csr_sexp) < 7 ||
      TYPEOF(pools_sexp) != VECSXP || XLENGTH(pools_sexp) < 3 ||
      TYPEOF(rows_sexp) != VECSXP || XLENGTH(rows_sexp) < 4) {
    return R_NilValue;
  }

  SEXP ps_key_s = VECTOR_ELT(csr_sexp, 0);
  SEXP pt_key_s = VECTOR_ELT(csr_sexp, 1);
  SEXP p_starts_s = VECTOR_ELT(csr_sexp, 2);
  SEXP p_ends_s = VECTOR_ELT(csr_sexp, 3);
  SEXP p_ki_s = VECTOR_ELT(csr_sexp, 4);
  SEXP p_cnt_s = VECTOR_ELT(csr_sexp, 5);
  SEXP p_cid_s = VECTOR_ELT(csr_sexp, 6);

  SEXP pl_entries_s = VECTOR_ELT(pools_sexp, 0);
  SEXP pl_start_s = VECTOR_ELT(pools_sexp, 1);
  SEXP pl_size_s = VECTOR_ELT(pools_sexp, 2);

  SEXP r_gi_s = VECTOR_ELT(rows_sexp, 0);
  SEXP r_si_s = VECTOR_ELT(rows_sexp, 1);
  SEXP r_ti_s = VECTOR_ELT(rows_sexp, 2);
  SEXP r_pid_s = VECTOR_ELT(rows_sexp, 3);

  if (TYPEOF(ps_key_s) != INTSXP || TYPEOF(pt_key_s) != INTSXP ||
      TYPEOF(p_starts_s) != INTSXP || TYPEOF(p_ends_s) != INTSXP ||
      TYPEOF(p_ki_s) != INTSXP || TYPEOF(p_cnt_s) != REALSXP || TYPEOF(p_cid_s) != INTSXP ||
      TYPEOF(pl_entries_s) != INTSXP || TYPEOF(pl_start_s) != INTSXP || TYPEOF(pl_size_s) != INTSXP ||
      TYPEOF(r_gi_s) != INTSXP || TYPEOF(r_si_s) != INTSXP ||
      TYPEOF(r_ti_s) != INTSXP || TYPEOF(r_pid_s) != INTSXP) {
    return R_NilValue;
  }

  int n_pairs = (int) XLENGTH(ps_key_s);
  int n_occ = (int) XLENGTH(p_ki_s);
  int n_pools = (int) XLENGTH(pl_start_s);
  int n_pl_entries = (int) XLENGTH(pl_entries_s);
  int n_rows = (int) XLENGTH(r_gi_s);

  const int *ps_key = INTEGER(ps_key_s);
  const int *pt_key = INTEGER(pt_key_s);
  const int *p_starts = INTEGER(p_starts_s);
  const int *p_ends = INTEGER(p_ends_s);
  const int *p_ki = INTEGER(p_ki_s);
  const double *p_cnt = REAL(p_cnt_s);
  const int *p_cid = INTEGER(p_cid_s);

  const int *pl_entries = INTEGER(pl_entries_s);
  const int *pl_start = INTEGER(pl_start_s);
  const int *pl_size = INTEGER(pl_size_s);

  const int *r_gi = INTEGER(r_gi_s);
  const int *r_si = INTEGER(r_si_s);
  const int *r_ti = INTEGER(r_ti_s);
  const int *r_pid = INTEGER(r_pid_s);

  SEXP target_cid_s = PROTECT(Rf_allocVector(INTSXP, n_rows));
  SEXP target_cnt_s = PROTECT(Rf_allocVector(REALSXP, n_rows));
  SEXP ctrl_cnt_s = PROTECT(Rf_allocVector(REALSXP, n_rows));
  SEXP risk_cnt_s = PROTECT(Rf_allocVector(REALSXP, n_rows));
  SEXP all_cnt_s = PROTECT(Rf_allocVector(REALSXP, n_rows));
  SEXP n_sup_s = PROTECT(Rf_allocVector(INTSXP, n_rows));

  int *target_cid = INTEGER(target_cid_s);
  double *target_cnt = REAL(target_cnt_s);
  double *ctrl_cnt = REAL(ctrl_cnt_s);
  double *risk_cnt = REAL(risk_cnt_s);
  double *all_cnt = REAL(all_cnt_s);
  int *n_sup = INTEGER(n_sup_s);

  for (int r = 0; r < n_rows; ++r) {
    int gi_r = r_gi[r];
    int si_r = r_si[r];
    int ti_r = r_ti[r];
    int pid_r = r_pid[r];
    if (pid_r < 1 || pid_r > n_pools) {
      UNPROTECT(6);
      return R_NilValue;
    }

    int lo = 0, hi = n_pairs - 1, p_pos = -1;
    while (lo <= hi) {
      int mid = lo + ((hi - lo) >> 1);
      int sm = ps_key[mid], tm = pt_key[mid];
      if (sm == si_r && tm == ti_r) {
        p_pos = mid;
        break;
      }
      if (sm < si_r || (sm == si_r && tm < ti_r)) lo = mid + 1;
      else hi = mid - 1;
    }
    if (p_pos < 0) {
      UNPROTECT(6);
      return R_NilValue;
    }

    int s0 = p_starts[p_pos] - 1;
    int s1 = p_ends[p_pos] - 1;
    if (s0 < 0 || s1 < s0 || s1 >= n_occ) {
      UNPROTECT(6);
      return R_NilValue;
    }

    lo = s0;
    hi = s1;
    int t_hit = -1;
    while (lo <= hi) {
      int mid = lo + ((hi - lo) >> 1);
      int kv = p_ki[mid];
      if (kv == gi_r) {
        t_hit = mid;
        break;
      }
      if (kv < gi_r) lo = mid + 1;
      else hi = mid - 1;
    }
    if (t_hit < 0 || !(p_cnt[t_hit] >= 1.0)) {
      UNPROTECT(6);
      return R_NilValue;
    }
    target_cid[r] = p_cid[t_hit];
    target_cnt[r] = p_cnt[t_hit];

    int pst = pl_start[pid_r - 1] - 1;
    int psz = pl_size[pid_r - 1];
    if (pst < 0 || psz < 0 || pst + psz > n_pl_entries) {
      UNPROTECT(6);
      return R_NilValue;
    }
    double c_sum = 0.0;
    for (int k = 0; k < psz; ++k) {
      int pk_val = pl_entries[pst + k];
      lo = s0;
      hi = s1;
      int c_hit = -1;
      while (lo <= hi) {
        int mid = lo + ((hi - lo) >> 1);
        int kv = p_ki[mid];
        if (kv == pk_val) {
          c_hit = mid;
          break;
        }
        if (kv < pk_val) lo = mid + 1;
        else hi = mid - 1;
      }
      if (c_hit < 0 || !(p_cnt[c_hit] >= 1.0)) {
        UNPROTECT(6);
        return R_NilValue;
      }
      c_sum += p_cnt[c_hit];
      if (c_sum > (double) CDID_MAX_EXACT_INT) {
        UNPROTECT(6);
        return R_NilValue;
      }
    }

    double a_sum = 0.0;
    for (int idx = s0; idx <= s1; ++idx) {
      a_sum += p_cnt[idx];
      if (a_sum > (double) CDID_MAX_EXACT_INT) {
        UNPROTECT(6);
        return R_NilValue;
      }
    }

    double rk_sum = target_cnt[r] + c_sum;
    if (rk_sum > (double) CDID_MAX_EXACT_INT) {
      UNPROTECT(6);
      return R_NilValue;
    }

    n_sup[r] = psz;
    ctrl_cnt[r] = c_sum;
    risk_cnt[r] = rk_sum;
    all_cnt[r] = a_sum;
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("target_occ_cell_id"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("target_pair_count"));
  SET_STRING_ELT(nms_s, 2, Rf_mkChar("control_pool_pair_count"));
  SET_STRING_ELT(nms_s, 3, Rf_mkChar("risk_pair_count"));
  SET_STRING_ELT(nms_s, 4, Rf_mkChar("all_cohort_pair_count"));
  SET_STRING_ELT(nms_s, 5, Rf_mkChar("n_supported_controls"));
  SET_VECTOR_ELT(out_s, 0, target_cid_s);
  SET_VECTOR_ELT(out_s, 1, target_cnt_s);
  SET_VECTOR_ELT(out_s, 2, ctrl_cnt_s);
  SET_VECTOR_ELT(out_s, 3, risk_cnt_s);
  SET_VECTOR_ELT(out_s, 4, all_cnt_s);
  SET_VECTOR_ELT(out_s, 5, n_sup_s);
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
  UNPROTECT(8);
  return out_s;
}

/* ---------------------------------------------------------------------------
   4. Fast payload budget gate: ipw_point_view_payload_budget
   --------------------------------------------------------------------------- */
typedef struct {
  double entries;
  double work;
  double bytes;
  double max_bytes;
  double cap_me;
  double cap_wo;
  int ok;
} cdid_pb_acc_t;

static void cdid_pb_walk(SEXP x, int depth, int scan_attrs, cdid_pb_acc_t *acc);

typedef struct {
  uint64_t n_attrs;
  int bail;
} cdid_pb_attr_cnt_t;

static SEXP cdid_pb_count_attr_cb(SEXP tag, SEXP value, void *data) {
  cdid_pb_attr_cnt_t *st = (cdid_pb_attr_cnt_t *) data;
  if (TYPEOF(tag) != SYMSXP) {
    st->bail = 1;
    return R_NilValue;
  }
  if (tag == R_DimSymbol && TYPEOF(value) == INTSXP && XLENGTH(value) == 1) {
    /* 1D arrays can synthesize names from dimnames in do_attributes() */
    st->bail = 1;
    return R_NilValue;
  }
  st->n_attrs += 1ULL;
  return NULL;
}

typedef struct {
  int depth;
  cdid_pb_acc_t *acc;
} cdid_pb_attr_walk_t;

static SEXP cdid_pb_visit_attr_cb(SEXP tag, SEXP value, void *data) {
  cdid_pb_attr_walk_t *w = (cdid_pb_attr_walk_t *) data;
  cdid_pb_acc_t *acc = w->acc;
  if (!acc->ok) return R_NilValue;
  if (tag == R_RowNamesSymbol && TYPEOF(value) == INTSXP && !Rf_inherits(value, "factor") &&
      XLENGTH(value) == 2 && INTEGER(value)[0] == NA_INTEGER) {
    /* Compact c(NA, -n) row.names: attributes(df)$row.names is integer(|n|) at depth + 2 */
    if (w->depth + 2 > 12) {
      acc->ok = 0;
      return R_NilValue;
    }
    int n = INTEGER(value)[1];
    double len = 0.0;
    if (n != NA_INTEGER && n != 0) len = (double) (n < 0 ? -(int64_t) n : (int64_t) n);
    double ent_new = acc->entries + len + 8.0;
    double wrk_new = acc->work + 64.0 * (len + 8.0);
    if (ent_new > acc->cap_me || wrk_new > acc->cap_wo) {
      acc->ok = 0;
      return R_NilValue;
    }
    acc->entries = ent_new;
    acc->work = wrk_new;
    return NULL;
  }
  cdid_pb_walk(value, w->depth + 2, 1, acc);
  return acc->ok ? NULL : R_NilValue;
}

static int cdid_pb_safe_object(SEXP x) {
  if (!Rf_isObject(x)) return 1;
  if (Rf_isS4(x)) return 0;
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  if (TYPEOF(cls) != STRSXP || XLENGTH(cls) < 1) return 0;
  R_xlen_t nc = XLENGTH(cls);
  for (R_xlen_t i = 0; i < nc; ++i) {
    SEXP s = STRING_ELT(cls, i);
    if (s == NA_STRING) return 0;
    const char *c = CHAR(s);
    if (strcmp(c, "data.frame") != 0 && strncmp(c, "cdid_", 5) != 0) return 0;
  }
  return 1;
}

static void cdid_pb_walk(SEXP x, int depth, int scan_attrs, cdid_pb_acc_t *acc) {
  if (!acc->ok) return;
  if (x == R_NilValue || TYPEOF(x) == ENVSXP) return;
  if (depth > 12 || ALTREP(x) || !cdid_pb_safe_object(x)) {
    acc->ok = 0;
    return;
  }
  int t = TYPEOF(x);
  if (t != VECSXP && t != STRSXP && t != REALSXP && t != INTSXP &&
      t != LGLSXP && t != RAWSXP && t != CPLXSXP) {
    acc->ok = 0;
    return;
  }

  double k = (double) XLENGTH(x);
  double ent_new = acc->entries + k + 8.0;
  double wrk_new = acc->work + 64.0 * (k + 8.0);
  if (ent_new > acc->cap_me || wrk_new > acc->cap_wo) {
    acc->ok = 0;
    return;
  }
  acc->entries = ent_new;
  acc->work = wrk_new;

  if (scan_attrs && cdid_has_attrib(x)) {
    if (depth + 1 > 12) {
      acc->ok = 0;
      return;
    }
    cdid_pb_attr_cnt_t cnt = {0ULL, 0};
    cdid_each_attrib(x, cdid_pb_count_attr_cb, &cnt);
    if (cnt.bail) {
      acc->ok = 0;
      return;
    }
    if (cnt.n_attrs > 0ULL) {
      double ka = (double) cnt.n_attrs;
      double ea = acc->entries + ka + 8.0;
      double wa = acc->work + 64.0 * (ka + 8.0);
      if (ea > acc->cap_me || wa > acc->cap_wo) {
        acc->ok = 0;
        return;
      }
      acc->entries = ea;
      acc->work = wa;
      cdid_pb_attr_walk_t aw = {depth, acc};
      cdid_each_attrib(x, cdid_pb_visit_attr_cb, &aw);
      if (!acc->ok) return;
    }
  }

  if (t == VECSXP) {
    R_xlen_t n = XLENGTH(x);
    for (R_xlen_t j = 0; j < n; ++j) {
      cdid_pb_walk(VECTOR_ELT(x, j), depth + 1, 1, acc);
      if (!acc->ok) return;
    }
  } else if (t == STRSXP) {
    R_xlen_t n = XLENGTH(x);
    for (R_xlen_t j = 0; j < n; ++j) {
      SEXP s = STRING_ELT(x, j);
      double nb = (s == NA_STRING) ? 0.0 : (double) LENGTH(s);
      double b_new = acc->bytes + nb;
      if (b_new > (double) CDID_MAX_EXACT_INT) {
        acc->ok = 0;
        return;
      }
      acc->bytes = b_new;
      if (nb > acc->max_bytes) acc->max_bytes = nb;
      double half = (s == NA_STRING) ? 0.0 : (double) (((uint64_t) LENGTH(s) + 1ULL) / 2ULL);
      double e2 = acc->entries + half;
      double w2 = acc->work + 16.0 * nb;
      if (e2 > acc->cap_me || w2 > acc->cap_wo) {
        acc->ok = 0;
        return;
      }
      acc->entries = e2;
      acc->work = w2;
    }
  }
}

SEXP C_cdid_point_payload_budget(SEXP objects_sexp, SEXP cap_me_sexp, SEXP cap_wo_sexp) {
  double cap_me = Rf_asReal(cap_me_sexp);
  double cap_wo = Rf_asReal(cap_wo_sexp);
  if (!R_FINITE(cap_me) || !R_FINITE(cap_wo)) return R_NilValue;
  const double limit = (double) CDID_MAX_EXACT_INT;
  if (cap_me > limit) cap_me = limit;
  if (cap_wo > limit) cap_wo = limit;

  cdid_pb_acc_t acc;
  acc.entries = 0.0;
  acc.work = 128.0;
  acc.bytes = 0.0;
  acc.max_bytes = 0.0;
  acc.cap_me = cap_me;
  acc.cap_wo = cap_wo;
  acc.ok = 1;

  cdid_pb_walk(objects_sexp, 0, 1, &acc);
  if (!acc.ok) return R_NilValue;

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("entries"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("work"));
  SET_STRING_ELT(nms_s, 2, Rf_mkChar("bytes"));
  SET_STRING_ELT(nms_s, 3, Rf_mkChar("max_bytes"));
  SET_VECTOR_ELT(out_s, 0, Rf_ScalarReal(acc.entries));
  SET_VECTOR_ELT(out_s, 1, Rf_ScalarReal(acc.work));
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarReal(acc.bytes));
  SET_VECTOR_ELT(out_s, 3, Rf_ScalarReal(acc.max_bytes));
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
  UNPROTECT(2);
  return out_s;
}

/* ---------------------------------------------------------------------------
   5. Numeric unit-ID hash/bit tracker & weighted unit-batch IF accumulator:
      project_gmm_if_stream
   --------------------------------------------------------------------------- */
static uint64_t cdid_mix64(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

static int cdid_is_plain_numeric_vec(SEXP x, R_xlen_t expected_len) {
  if (Rf_isObject(x) || Rf_getAttrib(x, R_DimSymbol) != R_NilValue) return 0;
  int t = TYPEOF(x);
  if (t != REALSXP && t != INTSXP) return 0;
  return XLENGTH(x) == expected_len;
}

static SEXP cdid_format_double_vec(const uint64_t *bits, int n) {
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  char buf[64];
  for (int i = 0; i < n; ++i) {
    uint64_t u = bits[i];
    double d = 0.0;
    memcpy(&d, &u, sizeof(double));
    /* Fast integer-to-ASCII for exact integers |d| <= 2^53 - 1, excluding -0.0 */
    if (u != 0x8000000000000000ULL && fabs(d) <= (double) CDID_MAX_EXACT_INT && d == floor(d)) {
      int64_t val = (int64_t) d;
      char *p = buf + sizeof(buf);
      uint64_t mag = (val < 0) ? (uint64_t) (-(val + 1)) + 1ULL : (uint64_t) val;
      do {
        *(--p) = (char) ('0' + (int) (mag % 10ULL));
        mag /= 10ULL;
      } while (mag > 0ULL);
      if (val < 0) *(--p) = '-';
      SET_STRING_ELT(out, i, Rf_mkCharLen(p, (int) ((buf + sizeof(buf)) - p)));
    } else {
      int len = snprintf(buf, sizeof(buf), "%.17g", d);
      SET_STRING_ELT(out, i, Rf_mkCharLen(buf, len));
    }
  }
  UNPROTECT(1);
  return out;
}

/*
 * op = 0: init tracker for `ni` rows -> list(bits_raw, hash_slots_int, state_int)
 * op = 1: pass 1 batch check+insert -> logical(B) closed_flags (or NULL if invalid/non-numeric)
 * op = 2: pass 2 batch check against stored bits -> logical(2) c(ok_unit, ok_meta_w) or NULL
 * op = 3: compare expected numeric/character IDs against stored bits -> ScalarLogical
 * op = 4: format first `len` stored bits as `%.17g` STRSXP
 * op = 5: single-call check of length-`ni` numeric unit IDs `x_sexp` (finite + no duplicates)
 *         and optional `%.17g` STRSXP formatting when `extra_sexp` is TRUE
 */
SEXP C_cdid_num_id_tracker(
    SEXP op_sexp,
    SEXP tracker_sexp,
    SEXP x_sexp,
    SEXP offset_len_sexp,
    SEXP extra_sexp) {
  int op = Rf_asInteger(op_sexp);
  if (op == 0) {
    int ni = Rf_asInteger(offset_len_sexp);
    if (ni < 1) return R_NilValue;
    int cap = 16;
    while (cap < 2 * ni && cap < (1 << 29)) cap <<= 1;
    SEXP bits_raw = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t) ni * 8));
    SEXP hslots = PROTECT(Rf_allocVector(INTSXP, cap));
    SEXP state = PROTECT(Rf_allocVector(INTSXP, 2)); /* [0] = ni, [1] = has_dup */
    int *hs = INTEGER(hslots);
    for (int i = 0; i < cap; ++i) hs[i] = 0;
    INTEGER(state)[0] = ni;
    INTEGER(state)[1] = 0;
    SEXP tr = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(tr, 0, bits_raw);
    SET_VECTOR_ELT(tr, 1, hslots);
    SET_VECTOR_ELT(tr, 2, state);
    UNPROTECT(4);
    return tr;
  }

  if (op == 5) {
    int ni = Rf_asInteger(offset_len_sexp);
    if (ni < 1 || !cdid_is_plain_numeric_vec(x_sexp, (R_xlen_t) ni)) return R_NilValue;
    int want_fmt = (TYPEOF(extra_sexp) == LGLSXP && XLENGTH(extra_sexp) == 1 && LOGICAL(extra_sexp)[0] == 1);
    const int *xi = (TYPEOF(x_sexp) == INTSXP) ? INTEGER(x_sexp) : NULL;
    const double *xd = (TYPEOF(x_sexp) == REALSXP) ? REAL(x_sexp) : NULL;

    /* Fast path: check finiteness and whether strictly monotonic (common case for a_id / f_id) */
    int strictly_inc = 1;
    double prev_d = 0.0;
    for (int i = 0; i < ni; ++i) {
      double d;
      if (xi != NULL) {
        if (xi[i] == NA_INTEGER) return R_NilValue;
        d = (double) xi[i];
      } else {
        if (!R_FINITE(xd[i])) return R_NilValue;
        d = xd[i];
      }
      if (i > 0 && !(d > prev_d)) strictly_inc = 0;
      prev_d = d;
    }

    const void *vmax = vmaxget();
    uint64_t *bits = NULL;
    if (!strictly_inc || want_fmt) {
      bits = (uint64_t *) R_alloc((size_t) ni, sizeof(uint64_t));
      for (int i = 0; i < ni; ++i) {
        double d = (xi != NULL) ? (double) xi[i] : xd[i];
        bits[i] = cdid_dbl_bits(d);
      }
    }
    if (!strictly_inc) {
      int cap = 16;
      while (cap < 2 * ni && cap < (1 << 29)) cap <<= 1;
      int mask = cap - 1;
      int *hs = (int *) R_alloc((size_t) cap, sizeof(int));
      memset(hs, 0, (size_t) cap * sizeof(int));
      for (int i = 0; i < ni; ++i) {
        uint64_t u = bits[i];
        int slot = (int) (cdid_mix64(u) & (uint64_t) mask);
        for (;;) {
          int prev = hs[slot];
          if (prev == 0) {
            hs[slot] = i + 1;
            break;
          }
          if (bits[prev - 1] == u) {
            vmaxset(vmax);
            return R_NilValue;
          }
          slot = (slot + 1) & mask;
        }
      }
    }
    SEXP out = want_fmt ? cdid_format_double_vec(bits, ni) : Rf_ScalarLogical(1);
    vmaxset(vmax);
    return out;
  }

  if (TYPEOF(tracker_sexp) != VECSXP || XLENGTH(tracker_sexp) != 3) return R_NilValue;
  SEXP bits_raw = VECTOR_ELT(tracker_sexp, 0);
  SEXP hslots = VECTOR_ELT(tracker_sexp, 1);
  SEXP state = VECTOR_ELT(tracker_sexp, 2);
  uint64_t *bits = (uint64_t *) RAW(bits_raw);
  int *hs = INTEGER(hslots);
  int mask = (int) XLENGTH(hslots) - 1;
  int ni = INTEGER(state)[0];

  if (op == 1) {
    if (TYPEOF(offset_len_sexp) != INTSXP || XLENGTH(offset_len_sexp) < 2) return R_NilValue;
    int nin = INTEGER(offset_len_sexp)[0];
    int B = INTEGER(offset_len_sexp)[1];
    if (nin < 0 || B < 1 || nin + B > ni) return R_NilValue;
    if (!cdid_is_plain_numeric_vec(x_sexp, (R_xlen_t) B)) return R_NilValue;

    /* Optional cluster_id check when input_row_level == "cluster" */
    if (extra_sexp != R_NilValue) {
      if (!cdid_is_plain_numeric_vec(extra_sexp, (R_xlen_t) B)) return R_NilValue;
    }

    const int *xi = (TYPEOF(x_sexp) == INTSXP) ? INTEGER(x_sexp) : NULL;
    const double *xd = (TYPEOF(x_sexp) == REALSXP) ? REAL(x_sexp) : NULL;
    for (int i = 0; i < B; ++i) {
      if (xi != NULL) {
        if (xi[i] == NA_INTEGER) return R_NilValue;
      } else {
        if (!R_FINITE(xd[i])) return R_NilValue;
      }
    }
    if (extra_sexp != R_NilValue) {
      const int *ci = (TYPEOF(extra_sexp) == INTSXP) ? INTEGER(extra_sexp) : NULL;
      const double *cd = (TYPEOF(extra_sexp) == REALSXP) ? REAL(extra_sexp) : NULL;
      for (int i = 0; i < B; ++i) {
        double du = (xi != NULL) ? (double) xi[i] : xd[i];
        double dc;
        if (ci != NULL) {
          if (ci[i] == NA_INTEGER) return R_NilValue;
          dc = (double) ci[i];
        } else {
          if (!R_FINITE(cd[i])) return R_NilValue;
          dc = cd[i];
        }
        if (cdid_dbl_bits(du) != cdid_dbl_bits(dc)) return R_NilValue;
      }
    }

    SEXP flags_s = PROTECT(Rf_allocVector(LGLSXP, B));
    int *flags = LOGICAL(flags_s);
    for (int i = 0; i < B; ++i) {
      double d = (xi != NULL) ? (double) xi[i] : xd[i];
      uint64_t u = cdid_dbl_bits(d);
      int row_pos = nin + i;
      bits[row_pos] = u;

      int slot = (int) (cdid_mix64(u) & (uint64_t) mask);
      int dup = 0;
      for (;;) {
        int prev = hs[slot];
        if (prev == 0) {
          hs[slot] = row_pos + 1;
          break;
        }
        if (bits[prev - 1] == u) {
          dup = 1;
          INTEGER(state)[1] = 1;
          break;
        }
        slot = (slot + 1) & mask;
      }
      flags[i] = dup;
    }
    UNPROTECT(1);
    return flags_s;
  }

  if (op == 2) {
    if (TYPEOF(offset_len_sexp) != INTSXP || XLENGTH(offset_len_sexp) < 2) return R_NilValue;
    int nin = INTEGER(offset_len_sexp)[0];
    int B = INTEGER(offset_len_sexp)[1];
    if (nin < 0 || B < 1 || nin + B > ni) return R_NilValue;
    if (!cdid_is_plain_numeric_vec(x_sexp, (R_xlen_t) B)) return R_NilValue;
    if (extra_sexp != R_NilValue && !cdid_is_plain_numeric_vec(extra_sexp, (R_xlen_t) B)) return R_NilValue;

    const int *xi = (TYPEOF(x_sexp) == INTSXP) ? INTEGER(x_sexp) : NULL;
    const double *xd = (TYPEOF(x_sexp) == REALSXP) ? REAL(x_sexp) : NULL;
    for (int i = 0; i < B; ++i) {
      double d;
      if (xi != NULL) {
        if (xi[i] == NA_INTEGER) return R_NilValue;
        d = (double) xi[i];
      } else {
        if (!R_FINITE(xd[i])) return R_NilValue;
        d = xd[i];
      }
      if (cdid_dbl_bits(d) != bits[nin + i]) {
        return Rf_ScalarLogical(0);
      }
    }
    if (extra_sexp != R_NilValue) {
      const int *ci = (TYPEOF(extra_sexp) == INTSXP) ? INTEGER(extra_sexp) : NULL;
      const double *cd = (TYPEOF(extra_sexp) == REALSXP) ? REAL(extra_sexp) : NULL;
      for (int i = 0; i < B; ++i) {
        double dc;
        if (ci != NULL) {
          if (ci[i] == NA_INTEGER) return R_NilValue;
          dc = (double) ci[i];
        } else {
          if (!R_FINITE(cd[i])) return R_NilValue;
          dc = cd[i];
        }
        if (cdid_dbl_bits(dc) != bits[nin + i]) return R_NilValue;
      }
    }
    return Rf_ScalarLogical(1);
  }

  if (op == 3) {
    int n_exp = Rf_asInteger(offset_len_sexp);
    if (n_exp < 0 || n_exp > ni) return R_NilValue;
    if (!cdid_is_plain_numeric_vec(x_sexp, (R_xlen_t) n_exp)) return R_NilValue;
    const int *xi = (TYPEOF(x_sexp) == INTSXP) ? INTEGER(x_sexp) : NULL;
    const double *xd = (TYPEOF(x_sexp) == REALSXP) ? REAL(x_sexp) : NULL;
    for (int i = 0; i < n_exp; ++i) {
      double d;
      if (xi != NULL) {
        if (xi[i] == NA_INTEGER) return R_NilValue;
        d = (double) xi[i];
      } else {
        if (!R_FINITE(xd[i])) return R_NilValue;
        d = xd[i];
      }
      if (cdid_dbl_bits(d) != bits[i]) return Rf_ScalarLogical(0);
    }
    return Rf_ScalarLogical(1);
  }

  if (op == 4) {
    int n_fmt = Rf_asInteger(offset_len_sexp);
    if (n_fmt < 0 || n_fmt > ni) return R_NilValue;
    return cdid_format_double_vec(bits, n_fmt);
  }

  return R_NilValue;
}

/* Weighted unit-row accumulator matching C_gmm_if_unit_batch 0-ULP arithmetic.
   pass == 3: single-call unweighted two-pass unit accumulator over full `N x J` `psi_sexp`
              (with `w_sexp` holding `N x m_R` `phi_full` for finiteness verification and
              `flags_sexp` holding `ScalarLogical(centered)`). */
SEXP C_cdid_if_unit_batch_weighted(
    SEXP pass_sexp,
    SEXP psi_sexp,
    SEXP w_sexp,
    SEXP flags_sexp,
    SEXP counts_sexp,
    SEXP weighted_sum_sexp,
    SEXP raw_sum_sexp,
    SEXP psi_mean_sexp,
    SEXP S_sexp) {
  int pass = Rf_asInteger(pass_sexp);
  if (pass == 3) {
    if (TYPEOF(psi_sexp) != REALSXP || Rf_isObject(psi_sexp) || !Rf_isMatrix(psi_sexp) ||
        TYPEOF(w_sexp) != REALSXP || !Rf_isMatrix(w_sexp) ||
        TYPEOF(flags_sexp) != LGLSXP || XLENGTH(flags_sexp) != 1 || LOGICAL(flags_sexp)[0] == NA_LOGICAL) {
      return R_NilValue;
    }
    int N = Rf_nrows(w_sexp);
    int m_R = Rf_ncols(w_sexp);
    int J = Rf_ncols(psi_sexp);
    int compute_psi = 0;
    int B_step = N;
    if (TYPEOF(counts_sexp) == LGLSXP && XLENGTH(counts_sexp) == 1 && LOGICAL(counts_sexp)[0] == 1) {
      compute_psi = 1;
      B_step = N;
    } else if (TYPEOF(counts_sexp) == INTSXP && XLENGTH(counts_sexp) == 1 && INTEGER(counts_sexp)[0] >= 1) {
      compute_psi = 1;
      B_step = INTEGER(counts_sexp)[0];
    }
    if (N < 1 || J < 1 || m_R < 1) return R_NilValue;
    if (compute_psi ? (Rf_nrows(psi_sexp) != m_R) : (Rf_nrows(psi_sexp) != N)) return R_NilValue;
    int centered = (LOGICAL(flags_sexp)[0] != 0);

    const double *phi = REAL(w_sexp);
    R_xlen_t tot_phi = (R_xlen_t) N * (R_xlen_t) m_R;
    for (R_xlen_t idx = 0; idx < tot_phi; ++idx) {
      if (!R_FINITE(phi[idx])) return R_NilValue;
    }

    int nprot = 0;
    SEXP psi_out_s = R_NilValue;
    const double *psi = NULL;
    if (compute_psi) {
      const double *lam = REAL(psi_sexp);
      R_xlen_t tot_lam = (R_xlen_t) m_R * (R_xlen_t) J;
      for (R_xlen_t idx = 0; idx < tot_lam; ++idx) {
        if (!R_FINITE(lam[idx])) return R_NilValue;
      }
      psi_out_s = PROTECT(Rf_allocMatrix(REALSXP, N, J));
      ++nprot;
      double *psi_mat = REAL(psi_out_s);
      BLAS_INT N_b = (BLAS_INT) N;
      BLAS_INT J_b = (BLAS_INT) J;
      BLAS_INT mR_b = (BLAS_INT) m_R;
      double one_d = 1.0, zero_d = 0.0;
      for (int s0 = 0; s0 < N; s0 += B_step) {
        BLAS_INT b_cur = (BLAS_INT) ((s0 + B_step <= N) ? B_step : (N - s0));
        F77_CALL(dgemm)("N", "N", &b_cur, &J_b, &mR_b, &one_d,
                        phi + s0, &N_b, lam, &mR_b, &zero_d, psi_mat + s0, &N_b FCONE FCONE);
      }
      psi = psi_mat;
    } else {
      psi = REAL(psi_sexp);
    }

    R_xlen_t tot_psi = (R_xlen_t) N * (R_xlen_t) J;
    for (R_xlen_t idx = 0; idx < tot_psi; ++idx) {
      if (!R_FINITE(psi[idx])) {
        if (nprot > 0) UNPROTECT(nprot);
        return R_NilValue;
      }
    }

    SEXP ws_out_s = PROTECT(Rf_allocVector(REALSXP, J));
    SEXP rs_out_s = PROTECT(Rf_allocVector(REALSXP, J));
    SEXP pm_out_s = PROTECT(Rf_allocVector(REALSXP, J));
    SEXP S_out_s = PROTECT(Rf_allocMatrix(REALSXP, J, J));
    nprot += 4;
    double *ws_out = REAL(ws_out_s);
    double *rs_out = REAL(rs_out_s);
    double *pm_out = REAL(pm_out_s);
    double *S_out = REAL(S_out_s);
    memset(S_out, 0, (size_t) J * (size_t) J * sizeof(double));

    double sum_w = (double) N;
    for (int j = 0; j < J; ++j) {
      const double *pcol = psi + (size_t) j * (size_t) N;
      double s = 0.0;
      for (int i = 0; i < N; ++i) {
        s += pcol[i];
      }
      ws_out[j] = s;
      rs_out[j] = s;
      pm_out[j] = centered ? (s / sum_w) : 0.0;
      if (!R_FINITE(s) || !R_FINITE(pm_out[j])) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
    }

    const void *vmax = vmaxget();
    double *z = (double *) R_alloc((size_t) J, sizeof(double));
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < J; ++j) {
        z[j] = psi[i + (R_xlen_t) j * N] - pm_out[j];
      }
      for (int k = 0; k < J; ++k) {
        double zk = z[k];
        double *S_col = S_out + (R_xlen_t) k * J;
        for (int j = 0; j < J; ++j) {
          S_col[j] += z[j] * zk;
        }
      }
    }
    vmaxset(vmax);

    R_xlen_t tot_S = (R_xlen_t) J * (R_xlen_t) J;
    for (R_xlen_t idx = 0; idx < tot_S; ++idx) {
      if (!R_FINITE(S_out[idx])) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
    }

    int n_out = compute_psi ? 6 : 5;
    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, n_out));
    SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, n_out));
    nprot += 2;
    SET_STRING_ELT(nms_s, 0, Rf_mkChar("sum_w"));
    SET_STRING_ELT(nms_s, 1, Rf_mkChar("weighted_sum"));
    SET_STRING_ELT(nms_s, 2, Rf_mkChar("raw_sum"));
    SET_STRING_ELT(nms_s, 3, Rf_mkChar("psi_mean"));
    SET_STRING_ELT(nms_s, 4, Rf_mkChar("S"));
    SET_VECTOR_ELT(out_s, 0, Rf_ScalarReal(sum_w));
    SET_VECTOR_ELT(out_s, 1, ws_out_s);
    SET_VECTOR_ELT(out_s, 2, rs_out_s);
    SET_VECTOR_ELT(out_s, 3, pm_out_s);
    SET_VECTOR_ELT(out_s, 4, S_out_s);
    if (compute_psi) {
      SET_STRING_ELT(nms_s, 5, Rf_mkChar("psi"));
      SET_VECTOR_ELT(out_s, 5, psi_out_s);
    }
    Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
    UNPROTECT(nprot);
    return out_s;
  }
  if ((pass != 1 && pass != 2) || TYPEOF(psi_sexp) != REALSXP || Rf_isObject(psi_sexp) ||
      TYPEOF(w_sexp) != REALSXP || TYPEOF(flags_sexp) != LGLSXP ||
      TYPEOF(counts_sexp) != REALSXP || XLENGTH(counts_sexp) != 4) {
    return R_NilValue;
  }
  SEXP dim = Rf_getAttrib(psi_sexp, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) return R_NilValue;
  int B = INTEGER(dim)[0];
  int J = INTEGER(dim)[1];
  if (B < 1 || J < 1 || XLENGTH(w_sexp) != (R_xlen_t) B || XLENGTH(flags_sexp) != (R_xlen_t) B) {
    return R_NilValue;
  }

  const double *psi = REAL(psi_sexp);
  const double *w = REAL(w_sexp);
  const int *flags = LOGICAL(flags_sexp);
  const double *counts_in = REAL(counts_sexp);
  int64_t nout = (int64_t) counts_in[0];
  int64_t no = (int64_t) counts_in[1];
  int64_t npos = (int64_t) counts_in[2];
  double sum_w = counts_in[3];

  for (int i = 0; i < B; ++i) {
    if (flags[i] != 0) return R_NilValue;
    for (int j = 0; j < J; ++j) {
      if (!R_FINITE(psi[i + (R_xlen_t) j * B])) return R_NilValue;
    }
  }
  if (nout + B > no) return R_NilValue;

  if (pass == 1) {
    if (TYPEOF(weighted_sum_sexp) != REALSXP || XLENGTH(weighted_sum_sexp) != (R_xlen_t) J ||
        TYPEOF(raw_sum_sexp) != REALSXP || XLENGTH(raw_sum_sexp) != (R_xlen_t) J) {
      return R_NilValue;
    }
    SEXP counts_out_s = PROTECT(Rf_allocVector(REALSXP, 4));
    SEXP ws_out_s = PROTECT(Rf_allocVector(REALSXP, J));
    SEXP rs_out_s = PROTECT(Rf_allocVector(REALSXP, J));
    double *ws_out = REAL(ws_out_s);
    double *rs_out = REAL(rs_out_s);
    memcpy(ws_out, REAL(weighted_sum_sexp), (size_t) J * sizeof(double));
    memcpy(rs_out, REAL(raw_sum_sexp), (size_t) J * sizeof(double));

    for (int i = 0; i < B; ++i) {
      ++nout;
      double wi = w[i];
      if (wi > 0.0) ++npos;
      sum_w = cdid_fp_add(sum_w, wi);
      if (!R_FINITE(sum_w)) {
        UNPROTECT(3);
        return R_NilValue;
      }
      for (int j = 0; j < J; ++j) {
        double v = psi[i + (R_xlen_t) j * B];
        double ws = cdid_fp_add(ws_out[j], cdid_fp_mul(wi, v));
        double rs = cdid_fp_add(rs_out[j], v);
        if (!R_FINITE(ws) || !R_FINITE(rs)) {
          UNPROTECT(3);
          return R_NilValue;
        }
        ws_out[j] = ws;
        rs_out[j] = rs;
      }
    }
    double *co = REAL(counts_out_s);
    co[0] = (double) nout;
    co[1] = (double) no;
    co[2] = (double) npos;
    co[3] = sum_w;

    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(nms_s, 0, Rf_mkChar("counts"));
    SET_STRING_ELT(nms_s, 1, Rf_mkChar("weighted_sum"));
    SET_STRING_ELT(nms_s, 2, Rf_mkChar("raw_sum"));
    SET_VECTOR_ELT(out_s, 0, counts_out_s);
    SET_VECTOR_ELT(out_s, 1, ws_out_s);
    SET_VECTOR_ELT(out_s, 2, rs_out_s);
    Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
    UNPROTECT(5);
    return out_s;
  } else {
    if (TYPEOF(psi_mean_sexp) != REALSXP || XLENGTH(psi_mean_sexp) != (R_xlen_t) J ||
        TYPEOF(S_sexp) != REALSXP || XLENGTH(S_sexp) != (R_xlen_t) J * (R_xlen_t) J) {
      return R_NilValue;
    }
    const void *vmax = vmaxget();
    double *z = (double *) R_alloc((size_t) J, sizeof(double));
    const double *pm = REAL(psi_mean_sexp);

    SEXP counts_out_s = PROTECT(Rf_allocVector(REALSXP, 4));
    SEXP S_out_s = PROTECT(Rf_allocMatrix(REALSXP, J, J));
    SEXP s_dn = Rf_getAttrib(S_sexp, R_DimNamesSymbol);
    if (s_dn != R_NilValue) Rf_setAttrib(S_out_s, R_DimNamesSymbol, s_dn);
    double *S_out = REAL(S_out_s);
    memcpy(S_out, REAL(S_sexp), (size_t) J * (size_t) J * sizeof(double));

    for (int i = 0; i < B; ++i) {
      ++nout;
      double wi = w[i];
      for (int j = 0; j < J; ++j) {
        z[j] = cdid_fp_sub(psi[i + (R_xlen_t) j * B], pm[j]);
      }
      for (int k = 0; k < J; ++k) {
        double zk = z[k];
        double *S_col = S_out + (R_xlen_t) k * J;
        for (int j = 0; j < J; ++j) {
          double zjk = cdid_fp_mul(z[j], zk);
          double wzjk = cdid_fp_mul(wi, zjk);
          double s_new = cdid_fp_add(S_col[j], wzjk);
          if (!R_FINITE(s_new)) {
            UNPROTECT(2);
            vmaxset(vmax);
            return R_NilValue;
          }
          S_col[j] = s_new;
        }
      }
    }
    double *co = REAL(counts_out_s);
    co[0] = (double) nout;
    co[1] = (double) no;
    co[2] = (double) npos;
    co[3] = sum_w;

    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms_s, 0, Rf_mkChar("counts"));
    SET_STRING_ELT(nms_s, 1, Rf_mkChar("S"));
    SET_VECTOR_ELT(out_s, 0, counts_out_s);
    SET_VECTOR_ELT(out_s, 1, S_out_s);
    Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
    UNPROTECT(4);
    vmaxset(vmax);
    return out_s;
  }
}

/* ---------------------------------------------------------------------------
   6. Column-pivoted Householder QR compression & pooled WLS solver:
      or_wls_compress / cdid_fast_or_baseline_pooled_fit
   --------------------------------------------------------------------------- */
static int cdid_wls_compress_buf(const double *src, int m, int k, double *dst, int *r_out) {
  if (m <= 0 || k <= 0) return 0;
  size_t tot = (size_t) m * (size_t) k;
  for (size_t idx = 0; idx < tot; ++idx) {
    if (!R_FINITE(src[idx])) return 0;
  }
  const void *vmax = vmaxget();
  double *A = (double *) R_alloc(tot, sizeof(double));
  memcpy(A, src, tot * sizeof(double));
  int r = (m < k) ? m : k;
  BLAS_INT mb = (BLAS_INT) m;
  BLAS_INT kb = (BLAS_INT) k;
  BLAS_INT *jpvt = (BLAS_INT *) R_alloc((size_t) k, sizeof(BLAS_INT));
  for (int j = 0; j < k; ++j) jpvt[j] = 0;
  double *tau = (double *) R_alloc((size_t) r, sizeof(double));
  BLAS_INT lwork = -1, info = 0;
  double wkopt = 0.0;
  F77_CALL(dgeqp3)(&mb, &kb, A, &mb, jpvt, tau, &wkopt, &lwork, &info);
  if (info != 0 || !R_FINITE(wkopt)) {
    vmaxset(vmax);
    return 0;
  }
  lwork = (BLAS_INT) wkopt;
  if (lwork < 3 * kb + 1) lwork = 3 * kb + 1;
  double *work = (double *) R_alloc((size_t) lwork, sizeof(double));
  F77_CALL(dgeqp3)(&mb, &kb, A, &mb, jpvt, tau, work, &lwork, &info);
  if (info != 0) {
    vmaxset(vmax);
    return 0;
  }
  memset(dst, 0, (size_t) r * (size_t) k * sizeof(double));
  for (int col = 0; col < k; ++col) {
    int orig_col = (int) jpvt[col] - 1;
    if (orig_col < 0 || orig_col >= k) {
      vmaxset(vmax);
      return 0;
    }
    double *dcol = dst + (size_t) orig_col * (size_t) r;
    const double *sc = A + (size_t) col * (size_t) m;
    int rlim = (col + 1 < r) ? (col + 1) : r;
    for (int i = 0; i < rlim; ++i) {
      if (!R_FINITE(sc[i])) {
        vmaxset(vmax);
        return 0;
      }
      dcol[i] = sc[i];
    }
  }
  vmaxset(vmax);
  *r_out = r;
  return 1;
}

static SEXP cdid_wls_pooled_fit_list(SEXP req_sexp) {
  if (XLENGTH(req_sexp) != 4) return R_NilValue;
  SEXP fac_list = VECTOR_ELT(req_sexp, 0);
  SEXP const_list = VECTOR_ELT(req_sexp, 1);
  SEXP p_s = VECTOR_ELT(req_sexp, 2);
  SEXP tol_s = VECTOR_ELT(req_sexp, 3);
  if (TYPEOF(fac_list) != VECSXP || TYPEOF(const_list) != VECSXP ||
      TYPEOF(p_s) != INTSXP || XLENGTH(p_s) != 1 ||
      TYPEOF(tol_s) != REALSXP || XLENGTH(tol_s) != 1) {
    return R_NilValue;
  }
  int n_b = (int) XLENGTH(fac_list);
  int p = INTEGER(p_s)[0];
  double rank_tol = REAL(tol_s)[0];
  if (n_b < 1 || XLENGTH(const_list) != n_b || p < 0 || p == NA_INTEGER ||
      !R_FINITE(rank_tol) || rank_tol <= 0.0) {
    return R_NilValue;
  }
  int p1 = p + 1;
  int p2 = p + 2;
  const void *vmax = vmaxget();

  double *acc_buf = (double *) R_alloc((size_t) p1 * (size_t) p1, sizeof(double));
  double *res_buf = (double *) R_alloc((size_t) p1 * (size_t) p1, sizeof(double));
  double *stk_buf = (double *) R_alloc((size_t) (2 * p1) * (size_t) p1, sizeof(double));
  int m_acc = 0;
  int max_aug = 0, max_res = 0, max_res_merge = 0, extra_comp = 0;

  SEXP r11_s = PROTECT(Rf_allocVector(REALSXP, n_b));
  SEXP r1d_s = PROTECT(Rf_allocMatrix(REALSXP, p, n_b));
  SEXP r1y_s = PROTECT(Rf_allocVector(REALSXP, n_b));
  SEXP res_rows_s = PROTECT(Rf_allocVector(INTSXP, n_b));
  double *r11_v = REAL(r11_s);
  double *r1d_m = REAL(r1d_s);
  double *r1y_v = REAL(r1y_s);
  int *res_rows_v = INTEGER(res_rows_s);

  BLAS_INT one_b = 1;
  BLAS_INT p1_b = (BLAS_INT) p1;

  for (int b = 0; b < n_b; ++b) {
    SEXP fb = VECTOR_ELT(fac_list, b);
    SEXP cb = VECTOR_ELT(const_list, b);
    if (TYPEOF(fb) != REALSXP || Rf_isObject(fb) || TYPEOF(cb) != LGLSXP || XLENGTH(cb) != p) {
      UNPROTECT(4);
      vmaxset(vmax);
      return R_NilValue;
    }
    SEXP dim_b = Rf_getAttrib(fb, R_DimSymbol);
    if (TYPEOF(dim_b) != INTSXP || XLENGTH(dim_b) != 2) {
      UNPROTECT(4);
      vmaxset(vmax);
      return R_NilValue;
    }
    int m_fac = INTEGER(dim_b)[0];
    int ncol_b = INTEGER(dim_b)[1];
    if (m_fac < 1 || ncol_b != p2) {
      UNPROTECT(4);
      vmaxset(vmax);
      return R_NilValue;
    }
    const double *f_ptr = REAL(fb);
    size_t ftot = (size_t) m_fac * (size_t) p2;
    for (size_t idx = 0; idx < ftot; ++idx) {
      if (!R_FINITE(f_ptr[idx])) {
        UNPROTECT(4);
        vmaxset(vmax);
        return R_NilValue;
      }
    }
    const void *b_vmax = vmaxget();
    double *A_int = (double *) R_alloc((size_t) m_fac, sizeof(double));
    memcpy(A_int, f_ptr, (size_t) m_fac * sizeof(double));
    BLAS_INT mb = (BLAS_INT) m_fac;
    BLAS_INT jpvt_int = 0;
    double tau_int = 0.0;
    BLAS_INT lwork = -1, info = 0;
    double wkopt = 0.0;
    F77_CALL(dgeqp3)(&mb, &one_b, A_int, &mb, &jpvt_int, &tau_int, &wkopt, &lwork, &info);
    if (info != 0 || !R_FINITE(wkopt)) {
      UNPROTECT(4);
      vmaxset(vmax);
      return R_NilValue;
    }
    lwork = (BLAS_INT) wkopt;
    if (lwork < 4) lwork = 4;
    double *work = (double *) R_alloc((size_t) lwork, sizeof(double));
    F77_CALL(dgeqp3)(&mb, &one_b, A_int, &mb, &jpvt_int, &tau_int, work, &lwork, &info);
    double r11 = A_int[0];
    if (info != 0 || !R_FINITE(r11) || fabs(r11) <= 0.0) {
      UNPROTECT(4);
      vmaxset(vmax);
      return R_NilValue;
    }

    size_t qty_tot = (size_t) m_fac * (size_t) p1;
    double *qty = (double *) R_alloc(qty_tot, sizeof(double));
    memcpy(qty, f_ptr + (size_t) m_fac, qty_tot * sizeof(double));
    lwork = -1;
    wkopt = 0.0;
    F77_CALL(dormqr)("L", "T", &mb, &p1_b, &one_b, A_int, &mb, &tau_int, qty, &mb, &wkopt, &lwork, &info FCONE FCONE);
    if (info != 0 || !R_FINITE(wkopt)) {
      UNPROTECT(4);
      vmaxset(vmax);
      return R_NilValue;
    }
    lwork = (BLAS_INT) wkopt;
    if (lwork < p1_b) lwork = p1_b;
    double *work_q = (double *) R_alloc((size_t) lwork, sizeof(double));
    F77_CALL(dormqr)("L", "T", &mb, &p1_b, &one_b, A_int, &mb, &tau_int, qty, &mb, work_q, &lwork, &info FCONE FCONE);
    if (info != 0) {
      UNPROTECT(4);
      vmaxset(vmax);
      return R_NilValue;
    }
    for (size_t idx = 0; idx < qty_tot; ++idx) {
      if (!R_FINITE(qty[idx])) {
        UNPROTECT(4);
        vmaxset(vmax);
        return R_NilValue;
      }
    }

    r11_v[b] = r11;
    for (int j = 0; j < p; ++j) {
      r1d_m[j + (size_t) b * (size_t) p] = qty[(size_t) j * (size_t) m_fac];
    }
    r1y_v[b] = qty[(size_t) p * (size_t) m_fac];

    int m_raw = m_fac - 1;
    int m_res = 0;
    const int *c_const = LOGICAL(cb);
    if (m_raw > 0) {
      double *res_raw = (double *) R_alloc((size_t) m_raw * (size_t) p1, sizeof(double));
      for (int j = 0; j < p1; ++j) {
        double *dst_c = res_raw + (size_t) j * (size_t) m_raw;
        if (j < p && c_const[j] != 0 && c_const[j] != NA_LOGICAL) {
          memset(dst_c, 0, (size_t) m_raw * sizeof(double));
        } else {
          memcpy(dst_c, qty + (size_t) j * (size_t) m_fac + 1, (size_t) m_raw * sizeof(double));
        }
      }
      if (!cdid_wls_compress_buf(res_raw, m_raw, p1, res_buf, &m_res)) {
        UNPROTECT(4);
        vmaxset(vmax);
        return R_NilValue;
      }
      extra_comp += 1;
    }
    res_rows_v[b] = m_res;

    if (b == 0) {
      m_acc = m_res;
      if (m_res > 0) {
        memcpy(acc_buf, res_buf, (size_t) m_res * (size_t) p1 * sizeof(double));
      }
      max_aug = m_fac;
      max_res = m_res;
      max_res_merge = m_res;
    } else {
      int stacked = m_acc + m_res;
      if (stacked > max_res_merge) max_res_merge = stacked;
      if (m_fac > max_aug) max_aug = m_fac;
      if (m_res > max_res) max_res = m_res;
      if (stacked > 0) {
        for (int j = 0; j < p1; ++j) {
          double *scol = stk_buf + (size_t) j * (size_t) stacked;
          if (m_acc > 0) memcpy(scol, acc_buf + (size_t) j * (size_t) m_acc, (size_t) m_acc * sizeof(double));
          if (m_res > 0) memcpy(scol + m_acc, res_buf + (size_t) j * (size_t) m_res, (size_t) m_res * sizeof(double));
        }
        if (!cdid_wls_compress_buf(stk_buf, stacked, p1, acc_buf, &m_acc)) {
          UNPROTECT(4);
          vmaxset(vmax);
          return R_NilValue;
        }
        extra_comp += 1;
      } else {
        m_acc = 0;
      }
    }
    vmaxset(b_vmax);
  }

  SEXP coef_s = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP R_s = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP piv_s = PROTECT(Rf_allocVector(INTSXP, p));
  double *coef_v = REAL(coef_s);
  double *R_v = REAL(R_s);
  int *piv_v = INTEGER(piv_s);
  double rss = 0.0;
  double cutoff = 0.0;

  if (p > 0) {
    if (m_acc < p) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }
    double *X_qr = (double *) R_alloc((size_t) m_acc * (size_t) p, sizeof(double));
    memcpy(X_qr, acc_buf, (size_t) m_acc * (size_t) p * sizeof(double));
    BLAS_INT m_acc_b = (BLAS_INT) m_acc;
    BLAS_INT pb = (BLAS_INT) p;
    BLAS_INT *jpvt = (BLAS_INT *) R_alloc((size_t) p, sizeof(BLAS_INT));
    for (int j = 0; j < p; ++j) jpvt[j] = 0;
    double *tau = (double *) R_alloc((size_t) p, sizeof(double));
    BLAS_INT lwork = -1, info = 0;
    double wkopt = 0.0;
    F77_CALL(dgeqp3)(&m_acc_b, &pb, X_qr, &m_acc_b, jpvt, tau, &wkopt, &lwork, &info);
    if (info != 0 || !R_FINITE(wkopt)) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }
    lwork = (BLAS_INT) wkopt;
    if (lwork < 3 * pb + 1) lwork = 3 * pb + 1;
    double *work = (double *) R_alloc((size_t) lwork, sizeof(double));
    F77_CALL(dgeqp3)(&m_acc_b, &pb, X_qr, &m_acc_b, jpvt, tau, work, &lwork, &info);
    if (info != 0) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }

    memset(R_v, 0, (size_t) p * (size_t) p * sizeof(double));
    double max_diag = 0.0;
    double min_diag = R_PosInf;
    for (int col = 0; col < p; ++col) {
      int pv = (int) jpvt[col];
      if (pv < 1 || pv > p) {
        UNPROTECT(7);
        vmaxset(vmax);
        return R_NilValue;
      }
      piv_v[col] = pv;
      for (int i = 0; i <= col; ++i) {
        double rv = X_qr[i + (size_t) col * (size_t) m_acc];
        if (!R_FINITE(rv)) {
          UNPROTECT(7);
          vmaxset(vmax);
          return R_NilValue;
        }
        R_v[i + (size_t) col * (size_t) p] = rv;
      }
      double dg = fabs(R_v[col + (size_t) col * (size_t) p]);
      if (dg > max_diag) max_diag = dg;
      if (dg < min_diag) min_diag = dg;
    }
    cutoff = rank_tol * max_diag;
    int rank = 0;
    if (cutoff > 0.0) {
      for (int col = 0; col < p; ++col) {
        if (fabs(R_v[col + (size_t) col * (size_t) p]) > cutoff) ++rank;
      }
    }
    if (rank < p || !(min_diag > 16.0 * (double) p * cutoff)) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }

    double *qty_s = (double *) R_alloc((size_t) m_acc, sizeof(double));
    memcpy(qty_s, acc_buf + (size_t) p * (size_t) m_acc, (size_t) m_acc * sizeof(double));
    lwork = -1;
    wkopt = 0.0;
    F77_CALL(dormqr)("L", "T", &m_acc_b, &one_b, &pb, X_qr, &m_acc_b, tau, qty_s, &m_acc_b, &wkopt, &lwork, &info FCONE FCONE);
    if (info != 0 || !R_FINITE(wkopt)) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }
    lwork = (BLAS_INT) wkopt;
    if (lwork < 1) lwork = 1;
    double *work_s = (double *) R_alloc((size_t) lwork, sizeof(double));
    F77_CALL(dormqr)("L", "T", &m_acc_b, &one_b, &pb, X_qr, &m_acc_b, tau, qty_s, &m_acc_b, work_s, &lwork, &info FCONE FCONE);
    if (info != 0) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }

    double *sol = (double *) R_alloc((size_t) p, sizeof(double));
    memcpy(sol, qty_s, (size_t) p * sizeof(double));
    double one_d = 1.0;
    F77_CALL(dtrsm)("L", "U", "N", "N", &pb, &one_b, &one_d, R_v, &pb, sol, &pb FCONE FCONE FCONE FCONE);
    for (int col = 0; col < p; ++col) {
      if (!R_FINITE(sol[col])) {
        UNPROTECT(7);
        vmaxset(vmax);
        return R_NilValue;
      }
      coef_v[piv_v[col] - 1] = sol[col];
    }
    rss = 0.0;
    for (int i = p; i < m_acc; ++i) {
      rss = cdid_fp_add(rss, qty_s[i] * qty_s[i]);
    }
    if (!R_FINITE(rss)) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }
  } else {
    rss = 0.0;
    for (int i = 0; i < m_acc; ++i) {
      rss = cdid_fp_add(rss, acc_buf[i] * acc_buf[i]);
    }
    if (!R_FINITE(rss)) {
      UNPROTECT(7);
      vmaxset(vmax);
      return R_NilValue;
    }
  }

  SEXP stats_s = PROTECT(Rf_allocVector(REALSXP, 2 + n_b));
  double *stats_v = REAL(stats_s);
  stats_v[0] = rss;
  stats_v[1] = cutoff;
  double one_d_icpt = 1.0;
  for (int b = 0; b < n_b; ++b) {
    long double dot_ld = 0.0L;
    for (int j = 0; j < p; ++j) {
      volatile double prod = r1d_m[j + (size_t) b * (size_t) p] * coef_v[j];
      dot_ld += (long double) prod;
    }
    double rhs = cdid_fp_sub(r1y_v[b], (double) dot_ld);
    double alpha_b = rhs;
    double r11_b = r11_v[b];
    F77_CALL(dtrsm)("L", "U", "N", "N", &one_b, &one_b, &one_d_icpt, &r11_b, &one_b, &alpha_b, &one_b FCONE FCONE FCONE FCONE);
    if (!R_FINITE(alpha_b)) {
      UNPROTECT(8);
      vmaxset(vmax);
      return R_NilValue;
    }
    stats_v[2 + b] = alpha_b;
  }

  SEXP acc_fac_s = PROTECT(Rf_allocMatrix(REALSXP, m_acc, p1));
  if (m_acc > 0) {
    memcpy(REAL(acc_fac_s), acc_buf, (size_t) m_acc * (size_t) p1 * sizeof(double));
  }
  SEXP acct_s = PROTECT(Rf_allocVector(INTSXP, 4));
  INTEGER(acct_s)[0] = max_aug;
  INTEGER(acct_s)[1] = max_res;
  INTEGER(acct_s)[2] = max_res_merge;
  INTEGER(acct_s)[3] = extra_comp;

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 10));
  SET_VECTOR_ELT(out_s, 0, acc_fac_s);
  SET_VECTOR_ELT(out_s, 1, r11_s);
  SET_VECTOR_ELT(out_s, 2, r1d_s);
  SET_VECTOR_ELT(out_s, 3, r1y_s);
  SET_VECTOR_ELT(out_s, 4, res_rows_s);
  SET_VECTOR_ELT(out_s, 5, acct_s);
  SET_VECTOR_ELT(out_s, 6, coef_s);
  SET_VECTOR_ELT(out_s, 7, R_s);
  SET_VECTOR_ELT(out_s, 8, piv_s);
  SET_VECTOR_ELT(out_s, 9, stats_s);

  UNPROTECT(11);
  vmaxset(vmax);
  return out_s;
}

SEXP C_cdid_wls_compress(SEXP aug_sexp) {
  if (TYPEOF(aug_sexp) == VECSXP && !Rf_isObject(aug_sexp)) {
    return cdid_wls_pooled_fit_list(aug_sexp);
  }
  if (TYPEOF(aug_sexp) != REALSXP || Rf_isObject(aug_sexp)) return R_NilValue;
  SEXP dim = Rf_getAttrib(aug_sexp, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) return R_NilValue;
  int m = INTEGER(dim)[0];
  int k = INTEGER(dim)[1];
  if (m <= 0 || k <= 0) return R_NilValue;
  int r = (m < k) ? m : k;
  SEXP factor_s = PROTECT(Rf_allocMatrix(REALSXP, r, k));
  int r_out = 0;
  if (!cdid_wls_compress_buf(REAL(aug_sexp), m, k, REAL(factor_s), &r_out)) {
    UNPROTECT(1);
    return R_NilValue;
  }
  UNPROTECT(1);
  return factor_s;
}

/* ---------------------------------------------------------------------------
   7. Fast logit kernels:
      C_cdid_logit_exact_alias, C_cdid_clip_logit, C_cdid_cp04_logit_fast
   --------------------------------------------------------------------------- */
SEXP C_cdid_logit_exact_alias(SEXP X_sexp) {
  if (TYPEOF(X_sexp) != REALSXP || Rf_isObject(X_sexp)) return R_NilValue;
  SEXP dim = Rf_getAttrib(X_sexp, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) return R_NilValue;
  int n = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];
  if (n < 0 || p < 0) return R_NilValue;
  if (p == 0) return Rf_ScalarLogical(0);
  if (n == 0) return Rf_ScalarLogical(1);

  const double *X = REAL(X_sexp);
  for (int j = 0; j < p; ++j) {
    const double *xj = X + (size_t) j * (size_t) n;
    int all_zero = 1;
    for (int i = 0; i < n; ++i) {
      double v = xj[i];
      if (ISNAN(v)) return R_NilValue;
      if (v != 0.0) {
        all_zero = 0;
        break;
      }
    }
    if (all_zero) return Rf_ScalarLogical(1);

    for (int k = 0; k < j; ++k) {
      const double *xk = X + (size_t) k * (size_t) n;
      int eq = 1;
      for (int i = 0; i < n; ++i) {
        if (xj[i] != xk[i]) {
          eq = 0;
          break;
        }
      }
      if (eq) return Rf_ScalarLogical(1);
    }
  }
  return Rf_ScalarLogical(0);
}

SEXP C_cdid_clip_logit(SEXP X_or_eta_sexp, SEXP beta_sexp) {
  const double clip_cap = 1.0 - 1e-6;
  const double clip_cap_odds = clip_cap / (1.0 - clip_cap);

  SEXP eta_s = R_NilValue;
  SEXP rn_s = R_NilValue;
  int n = 0;
  int nprot = 0;

  if (beta_sexp == R_NilValue) {
    if (TYPEOF(X_or_eta_sexp) != REALSXP || Rf_isObject(X_or_eta_sexp) || cdid_has_attrib(X_or_eta_sexp)) {
      return R_NilValue;
    }
    R_xlen_t nx = XLENGTH(X_or_eta_sexp);
    if (nx > (R_xlen_t) INT_MAX) return R_NilValue;
    n = (int) nx;
    const double *ein = REAL(X_or_eta_sexp);
    for (int i = 0; i < n; ++i) {
      if (!R_FINITE(ein[i])) return R_NilValue;
    }
    eta_s = X_or_eta_sexp;
  } else {
    if (TYPEOF(X_or_eta_sexp) != REALSXP || Rf_isObject(X_or_eta_sexp) ||
        TYPEOF(beta_sexp) != REALSXP || Rf_isObject(beta_sexp)) {
      return R_NilValue;
    }
    SEXP dim = Rf_getAttrib(X_or_eta_sexp, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) return R_NilValue;
    n = INTEGER(dim)[0];
    int p = INTEGER(dim)[1];
    if (n <= 1 || p < 0 || XLENGTH(beta_sexp) != (R_xlen_t) p) return R_NilValue;

    const double *X = REAL(X_or_eta_sexp);
    const double *beta = REAL(beta_sexp);
    for (int j = 0; j < p; ++j) {
      if (!R_FINITE(beta[j])) return R_NilValue;
    }
    R_xlen_t tot_x = (R_xlen_t) n * (R_xlen_t) p;
    for (R_xlen_t idx = 0; idx < tot_x; ++idx) {
      if (!R_FINITE(X[idx])) return R_NilValue;
    }
    eta_s = PROTECT(Rf_allocVector(REALSXP, n));
    ++nprot;
    double *eta = REAL(eta_s);
    if (p == 0) {
      memset(eta, 0, (size_t) n * sizeof(double));
    } else {
      BLAS_INT n_b = (BLAS_INT) n, p_b = (BLAS_INT) p, one_b = 1;
      double one_d = 1.0, zero_d = 0.0;
      F77_CALL(dgemm)("N", "N", &n_b, &one_b, &p_b, &one_d,
                      X, &n_b, beta, &p_b, &zero_d, eta, &n_b FCONE FCONE);
    }
    for (int i = 0; i < n; ++i) {
      if (!R_FINITE(eta[i])) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
    }
    SEXP dn = Rf_getAttrib(X_or_eta_sexp, R_DimNamesSymbol);
    if (dn != R_NilValue && TYPEOF(dn) == VECSXP && XLENGTH(dn) == 2) {
      SEXP r0 = VECTOR_ELT(dn, 0);
      if (r0 != R_NilValue && TYPEOF(r0) == STRSXP && XLENGTH(r0) == (R_xlen_t) n) {
        rn_s = r0;
        Rf_setAttrib(eta_s, R_NamesSymbol, rn_s);
      }
    }
  }

  SEXP odds_s = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP prob_s = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP act_s = PROTECT(Rf_allocVector(REALSXP, n));
  nprot += 3;
  if (rn_s != R_NilValue) {
    Rf_setAttrib(odds_s, R_NamesSymbol, rn_s);
    Rf_setAttrib(prob_s, R_NamesSymbol, rn_s);
  }

  const double *eta = REAL(eta_s);
  double *odds = REAL(odds_s);
  double *prob = REAL(prob_s);
  double *act = REAL(act_s);

  for (int i = 0; i < n; ++i) {
    double e = eta[i];
    double od, pr;
    if (e > 40.0) {
      pr = 1.0 / (1.0 + exp(-e));
      od = 1.0;
    } else {
      od = exp(e);
      pr = od / (1.0 + od);
    }
    if (pr < clip_cap) {
      act[i] = 1.0;
      odds[i] = od;
    } else {
      act[i] = 0.0;
      odds[i] = clip_cap_odds;
    }
    prob[i] = pr;
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 4));
  nprot += 2;
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("eta"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("odds"));
  SET_STRING_ELT(nms_s, 2, Rf_mkChar("prob"));
  SET_STRING_ELT(nms_s, 3, Rf_mkChar("active"));
  SET_VECTOR_ELT(out_s, 0, eta_s);
  SET_VECTOR_ELT(out_s, 1, odds_s);
  SET_VECTOR_ELT(out_s, 2, prob_s);
  SET_VECTOR_ELT(out_s, 3, act_s);
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
  UNPROTECT(nprot);
  return out_s;
}

/* Upper-triangular Cholesky R (R^T R = H) stored column-major in p x p R_out.
   Returns 1 if PD and well-conditioned (min_diag > rel_tol * max_diag), else 0. */
static int cdid_chol_upper(const double *H, double *R_out, int p, double rel_tol) {
  if (p < 1) return 0;
  memset(R_out, 0, (size_t) p * (size_t) p * sizeof(double));
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i <= j; ++i) {
      double hij = H[i + j * p];
      if (!R_FINITE(hij)) return 0;
      R_out[i + j * p] = hij;
    }
  }
  BLAS_INT pb = (BLAS_INT) p;
  BLAS_INT info = 0;
  F77_CALL(dpotrf)("U", &pb, R_out, &pb, &info FCONE);
  if (info != 0) return 0;
  double max_d = 0.0;
  double min_d = DBL_MAX;
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i <= j; ++i) {
      if (!R_FINITE(R_out[i + j * p])) return 0;
    }
    double rd = R_out[j + j * p];
    if (!(rd > 0.0)) return 0;
    if (rd > max_d) max_d = rd;
    if (rd < min_d) min_d = rd;
  }
  if (!(max_d > 0.0) || !(min_d > rel_tol * max_d)) return 0;
  return 1;
}

static void cdid_chol_solve(const double *R_up, const double *rhs, double *sol, int p) {
  /* Solve R^T y = rhs, then R sol = y */
  for (int i = 0; i < p; ++i) {
    double s = rhs[i];
    for (int k = 0; k < i; ++k) s -= R_up[k + i * p] * sol[k];
    sol[i] = s / R_up[i + i * p];
  }
  for (int i = p - 1; i >= 0; --i) {
    double s = sol[i];
    for (int k = i + 1; k < p; ++k) s -= R_up[i + k * p] * sol[k];
    sol[i] = s / R_up[i + i * p];
  }
}

/* Unpivoted Householder QR WLS solver matching stats::glm.fit's Cdqrls on
   full-rank weighted designs (W_X is n x p column-major, w_z is length n).
   Writes p x p R_out and solves R_out * beta_out = Q^T w_z.
   Returns 1 if full-rank and well-conditioned, 0 otherwise (caller falls back). */
static int cdid_qr_wls_solve(
    double *W_X,
    double *w_z,
    double *R_out,
    double *beta_out,
    int n,
    int p,
    double dqr_tol,
    double rel_tol,
    int keep_lower) {
  if (n < p || p < 1 || p > 256) return 0;
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
  double max_d = 0.0, min_d = DBL_MAX;
  for (int j = 0; j < p; ++j) {
    const double *col_j = W_X + (size_t) j * (size_t) n;
    int rlim = keep_lower ? p : (j + 1);
    for (int i = 0; i < rlim; ++i) {
      if (!R_FINITE(col_j[i])) return 0;
      R_out[i + j * p] = col_j[i];
    }
    double ad = fabs(col_j[j]);
    if (ad > max_d) max_d = ad;
    if (ad < min_d) min_d = ad;
  }
  if (!(max_d > 0.0) || !(min_d > rel_tol * max_d)) return 0;
  return 1;
}

SEXP C_cdid_cp04_logit_fast(
    SEXP X_sexp,
    SEXP y_sexp,
    SEXP w_sexp,
    SEXP eps_sexp,
    SEXP maxit_sexp,
    SEXP N_sexp) {
  if (TYPEOF(X_sexp) != REALSXP || Rf_isObject(X_sexp) ||
      (TYPEOF(y_sexp) != REALSXP && TYPEOF(y_sexp) != INTSXP) || Rf_isObject(y_sexp) ||
      TYPEOF(w_sexp) != REALSXP || Rf_isObject(w_sexp)) {
    return R_NilValue;
  }
  SEXP dim = Rf_getAttrib(X_sexp, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) return R_NilValue;
  int n = INTEGER(dim)[0];
  int p = INTEGER(dim)[1];
  if (p < 1 || p > 128 || n < p || XLENGTH(y_sexp) != (R_xlen_t) n || XLENGTH(w_sexp) != (R_xlen_t) n) {
    return R_NilValue;
  }

  double eps = Rf_asReal(eps_sexp);
  int maxit = Rf_asInteger(maxit_sexp);
  double N_val = Rf_asReal(N_sexp);
  if (!R_FINITE(eps) || eps <= 0.0 || maxit < 1 || !R_FINITE(N_val) || N_val <= 0.0) {
    return R_NilValue;
  }

  const double *X = REAL(X_sexp);
  const double *w = REAL(w_sexp);
  const void *vmax = vmaxget();
  double *y = (double *) R_alloc((size_t) n, sizeof(double));
  if (TYPEOF(y_sexp) == INTSXP) {
    const int *yi = INTEGER(y_sexp);
    for (int i = 0; i < n; ++i) {
      if (yi[i] != 0 && yi[i] != 1) {
        vmaxset(vmax);
        return R_NilValue;
      }
      y[i] = (double) yi[i];
    }
  } else {
    const double *yd = REAL(y_sexp);
    for (int i = 0; i < n; ++i) {
      if (yd[i] != 0.0 && yd[i] != 1.0) {
        vmaxset(vmax);
        return R_NilValue;
      }
      y[i] = yd[i];
    }
  }

  int n_pos_w = 0;
  long double sum_w_ld = 0.0L;
  long double devold_ld = 0.0L;
  const double ln2_val = log(2.0);
  for (int i = 0; i < n; ++i) {
    double wi = w[i];
    if (!R_FINITE(wi) || wi < 0.0) {
      vmaxset(vmax);
      return R_NilValue;
    }
    if (wi > 0.0) {
      ++n_pos_w;
      sum_w_ld += (long double) wi;
      devold_ld += (long double) (2.0 * wi * ln2_val);
    }
  }
  double sum_w = (double) sum_w_ld;
  if (n_pos_w < p || !R_FINITE(sum_w) || sum_w <= 0.0) {
    vmaxset(vmax);
    return R_NilValue;
  }

  /* Validate finite X and check signed-coordinate separation up front */
  for (int j = 0; j < p; ++j) {
    const double *xj = X + (size_t) j * (size_t) n;
    int has_pos = 0, has_neg = 0;
    for (int i = 0; i < n; ++i) {
      double v = xj[i];
      if (!R_FINITE(v)) {
        vmaxset(vmax);
        return R_NilValue;
      }
      if (w[i] > 0.0 && v != 0.0) {
        double sgn = (y[i] == 1.0) ? v : -v;
        if (sgn > 0.0) has_pos = 1;
        else if (sgn < 0.0) has_neg = 1;
      }
    }
    if ((has_pos && !has_neg) || (has_neg && !has_pos)) {
      vmaxset(vmax);
      return R_NilValue;
    }
  }

  double *eta = (double *) R_alloc((size_t) n, sizeof(double));
  double *mu = (double *) R_alloc((size_t) n, sizeof(double));
  double *wt = (double *) R_alloc((size_t) n, sizeof(double));
  double *wz = (double *) R_alloc((size_t) n, sizeof(double));
  double *W_X = (double *) R_alloc((size_t) n * (size_t) p, sizeof(double));
  double *H = (double *) R_alloc((size_t) p * (size_t) p, sizeof(double));
  double *R_chol = (double *) R_alloc((size_t) p * (size_t) p, sizeof(double));
  double *beta = (double *) R_alloc((size_t) p, sizeof(double));

  for (int i = 0; i < n; ++i) {
    eta[i] = 0.0;
    mu[i] = 0.5;
  }
  for (int j = 0; j < p; ++j) beta[j] = 0.0;

  double devold = (double) devold_ld;
  int conv = 0;
  int iter_used = 0;
  const double bnd_eps = 10.0 * DBL_EPSILON;
  const double qr_tol = (eps / 1000.0 < 1e-7) ? (eps / 1000.0) : 1e-7;
  int n_b = n, p_b = p, one_b = 1;
  double one_d = 1.0, zero_d = 0.0;

  for (int iter = 1; iter <= maxit; ++iter) {
    iter_used = iter;
    for (int i = 0; i < n; ++i) {
      double wi = w[i];
      if (wi > 0.0) {
        double mui = mu[i];
        double varmu = mui * (1.0 - mui);
        double opexp = exp(eta[i]);
        double mu_eta = opexp / ((1.0 + opexp) * (1.0 + opexp));
        if (!(varmu > 0.0) || !(mu_eta > 0.0) || !R_FINITE(mu_eta)) {
          vmaxset(vmax);
          return R_NilValue;
        }
        double w_i = sqrt((wi * (mu_eta * mu_eta)) / varmu);
        double z_i = eta[i] + (y[i] - mui) / mu_eta;
        wt[i] = w_i;
        wz[i] = z_i * w_i;
      } else {
        wt[i] = 0.0;
        wz[i] = 0.0;
      }
    }

    for (int j = 0; j < p; ++j) {
      const double *xj = X + (size_t) j * (size_t) n;
      double *wxj = W_X + (size_t) j * (size_t) n;
      for (int i = 0; i < n; ++i) wxj[i] = xj[i] * wt[i];
    }

    if (!cdid_qr_wls_solve(W_X, wz, R_chol, beta, n, p, qr_tol, 0.0, 0)) {
      vmaxset(vmax);
      return R_NilValue;
    }

    F77_CALL(dgemm)("N", "N", &n_b, &one_b, &p_b, &one_d,
                    X, &n_b, beta, &p_b, &zero_d, eta, &n_b FCONE FCONE);

    long double dev_ld = 0.0L;
    for (int i = 0; i < n; ++i) {
      double ei = eta[i];
      /* Stay strictly inside logit_linkinv / boundary thresholds */
      if (!R_FINITE(ei) || fabs(ei) >= 30.0) {
        vmaxset(vmax);
        return R_NilValue;
      }
      double opexp = exp(ei);
      double mui = opexp / (1.0 + opexp);
      if (mui <= bnd_eps || mui >= 1.0 - bnd_eps) {
        vmaxset(vmax);
        return R_NilValue;
      }
      mu[i] = mui;
      double wi = w[i];
      if (wi > 0.0) {
        double loss = (y[i] == 1.0) ? log(1.0 / mui) : log(1.0 / (1.0 - mui));
        dev_ld += (long double) (2.0 * wi * loss);
      }
    }
    double dev = (double) dev_ld;
    if (!R_FINITE(dev)) {
      vmaxset(vmax);
      return R_NilValue;
    }

    if (fabs(dev - devold) / (0.1 + fabs(dev)) < eps) {
      conv = 1;
      break;
    }
    devold = dev;
  }

  if (!conv) {
    vmaxset(vmax);
    return R_NilValue;
  }

  /* Final information matrix H_info = crossprod(X, wt * X) / N */
  for (int i = 0; i < n; ++i) {
    double ei = eta[i];
    double od = exp(ei);
    double pr = od / (1.0 + od);
    if (pr <= bnd_eps || pr >= 1.0 - bnd_eps || !(pr * (1.0 - pr) > 0.0)) {
      vmaxset(vmax);
      return R_NilValue;
    }
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
    if (!R_FINITE(H[idx])) {
      vmaxset(vmax);
      return R_NilValue;
    }
  }
  if (!cdid_chol_upper(H, R_chol, p, 0.0)) {
    vmaxset(vmax);
    return R_NilValue;
  }

  SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, p));
  memcpy(REAL(beta_s), beta, (size_t) p * sizeof(double));

  SEXP R_s = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  memcpy(REAL(R_s), R_chol, (size_t) p * (size_t) p * sizeof(double));
  SEXP xdn = Rf_getAttrib(X_sexp, R_DimNamesSymbol);
  if (xdn != R_NilValue && TYPEOF(xdn) == VECSXP && XLENGTH(xdn) == 2) {
    SEXP cn = VECTOR_ELT(xdn, 1);
    if (cn != R_NilValue) {
      SEXP rdn = PROTECT(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(rdn, 0, cn);
      SET_VECTOR_ELT(rdn, 1, cn);
      Rf_setAttrib(R_s, R_DimNamesSymbol, rdn);
      UNPROTECT(1);
    }
  }

  SEXP diag_s = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP dnms_s = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(dnms_s, 0, Rf_mkChar("converged"));
  SET_STRING_ELT(dnms_s, 1, Rf_mkChar("boundary"));
  SET_STRING_ELT(dnms_s, 2, Rf_mkChar("rank"));
  SET_STRING_ELT(dnms_s, 3, Rf_mkChar("iter"));
  SET_STRING_ELT(dnms_s, 4, Rf_mkChar("warnings"));
  SET_VECTOR_ELT(diag_s, 0, Rf_ScalarLogical(1));
  SET_VECTOR_ELT(diag_s, 1, Rf_ScalarLogical(0));
  SET_VECTOR_ELT(diag_s, 2, Rf_ScalarInteger(p));
  SET_VECTOR_ELT(diag_s, 3, Rf_ScalarInteger(iter_used));
  SET_VECTOR_ELT(diag_s, 4, Rf_allocVector(STRSXP, 0));
  Rf_setAttrib(diag_s, R_NamesSymbol, dnms_s);

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP onms_s = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(onms_s, 0, Rf_mkChar("coefficients"));
  SET_STRING_ELT(onms_s, 1, Rf_mkChar("H_chol"));
  SET_STRING_ELT(onms_s, 2, Rf_mkChar("diagnostics"));
  SET_VECTOR_ELT(out_s, 0, beta_s);
  SET_VECTOR_ELT(out_s, 1, R_s);
  SET_VECTOR_ELT(out_s, 2, diag_s);
  Rf_setAttrib(out_s, R_NamesSymbol, onms_s);

  UNPROTECT(6);
  vmaxset(vmax);
  return out_s;
}

/* ---------------------------------------------------------------------------
   8. L1-tiled centered Gram chunk kernel: cp04_stream_covariance
   --------------------------------------------------------------------------- */
SEXP C_cdid_center_gram(
    SEXP phi_sexp,
    SEXP means_sexp,
    SEXP map_a_sexp,
    SEXP map_b_sexp,
    SEXP sym_sexp,
    SEXP denom_sexp) {
  if (TYPEOF(phi_sexp) != REALSXP || Rf_isObject(phi_sexp) || TYPEOF(means_sexp) != REALSXP) {
    return R_NilValue;
  }
  SEXP dim = Rf_getAttrib(phi_sexp, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) return R_NilValue;
  int B = INTEGER(dim)[0];
  int K = INTEGER(dim)[1];
  if (B < 1 || K < 1 || XLENGTH(means_sexp) != (R_xlen_t) K) return R_NilValue;

  int sym = Rf_asLogical(sym_sexp);
  double denom = Rf_asReal(denom_sexp);
  if (sym == NA_LOGICAL || !R_FINITE(denom) || denom <= 0.0) return R_NilValue;

  const double *phi = REAL(phi_sexp);
  const double *means = REAL(means_sexp);
  int B_b = B;
  double one_d = 1.0, zero_d = 0.0;

  if (sym) {
    SEXP gram_s = PROTECT(Rf_allocMatrix(REALSXP, K, K));
    double *gram = REAL(gram_s);
    memset(gram, 0, (size_t) K * (size_t) K * sizeof(double));

    const void *vmax = vmaxget();
    double *xt = (double *) R_alloc((size_t) B * (size_t) K, sizeof(double));
    for (int c = 0; c < K; ++c) {
      const double *src = phi + (size_t) c * (size_t) B;
      double *dst = xt + (size_t) c * (size_t) B;
      double mc = means[c];
      for (int i = 0; i < B; ++i) dst[i] = src[i] - mc;
    }
    int K_b = K;
    F77_CALL(dsyrk)("U", "T", &K_b, &B_b, &one_d, xt, &B_b, &zero_d, gram, &K_b FCONE FCONE);
    for (int c2 = 0; c2 < K; ++c2) {
      for (int c1 = 0; c1 <= c2; ++c1) {
        double v = gram[c1 + (size_t) c2 * (size_t) K] / denom;
        gram[c1 + (size_t) c2 * (size_t) K] = v;
        gram[c2 + (size_t) c1 * (size_t) K] = v;
      }
    }
    vmaxset(vmax);
    UNPROTECT(1);
    return gram_s;
  } else {
    if (TYPEOF(map_a_sexp) != INTSXP || TYPEOF(map_b_sexp) != INTSXP) return R_NilValue;
    int na = (int) XLENGTH(map_a_sexp);
    int nb = (int) XLENGTH(map_b_sexp);
    if (na < 1 || nb < 1) return R_NilValue;
    const int *ma = INTEGER(map_a_sexp);
    const int *mb = INTEGER(map_b_sexp);
    for (int i = 0; i < na; ++i) if (ma[i] < 1 || ma[i] > K) return R_NilValue;
    for (int i = 0; i < nb; ++i) if (mb[i] < 1 || mb[i] > K) return R_NilValue;

    SEXP gram_s = PROTECT(Rf_allocMatrix(REALSXP, na, nb));
    double *gram = REAL(gram_s);
    memset(gram, 0, (size_t) na * (size_t) nb * sizeof(double));

    const void *vmax = vmaxget();
    double *xa = (double *) R_alloc((size_t) B * (size_t) na, sizeof(double));
    double *xb = (double *) R_alloc((size_t) B * (size_t) nb, sizeof(double));
    for (int c = 0; c < na; ++c) {
      int col = ma[c] - 1;
      const double *src = phi + (size_t) col * (size_t) B;
      double *dst = xa + (size_t) c * (size_t) B;
      double mc = means[col];
      for (int i = 0; i < B; ++i) dst[i] = src[i] - mc;
    }
    for (int c = 0; c < nb; ++c) {
      int col = mb[c] - 1;
      const double *src = phi + (size_t) col * (size_t) B;
      double *dst = xb + (size_t) c * (size_t) B;
      double mc = means[col];
      for (int i = 0; i < B; ++i) dst[i] = src[i] - mc;
    }
    int na_b = na, nb_b = nb;
    F77_CALL(dgemm)("T", "N", &na_b, &nb_b, &B_b, &one_d,
                    xa, &B_b, xb, &B_b, &zero_d, gram, &na_b FCONE FCONE);
    for (size_t idx = 0; idx < (size_t) na * (size_t) nb; ++idx) {
      gram[idx] /= denom;
    }
    vmaxset(vmax);
    UNPROTECT(1);
    return gram_s;
  }
}

/* ===========================================================================
   WAVE 2 KERNELS: Custom IPW, Raw IPW, Incomplete Baseline OR & DR
   =========================================================================== */

/* Helper: galloping / binary search for target `a` in strictly increasing
   1-based subarray post_a[lo0 .. hi0] (0-based indices lo0 <= hi0).
   Updates *hint_lo so monotone queries run in O(1) amortized time.
   Returns 0-based index if found, else -1. */
static R_INLINE int cdid_find_posting(
    const int *post_a,
    int lo0,
    int hi0,
    int target,
    int *hint_lo) {
  if (lo0 > hi0) return -1;
  int lo = *hint_lo;
  if (lo < lo0 || lo > hi0 || post_a[lo] > target) {
    lo = lo0;
  } else if (post_a[lo] == target) {
    *hint_lo = lo;
    return lo;
  } else {
    int step = 1;
    int cur = lo;
    while (cur + step <= hi0 && post_a[cur + step] < target) {
      cur += step;
      step <<= 1;
    }
    lo = cur + 1;
    int hi_cand = cur + step;
    if (hi_cand < hi0) hi0 = hi_cand;
  }
  int hi = hi0;
  while (lo <= hi) {
    int mid = lo + ((hi - lo) >> 1);
    int v = post_a[mid];
    if (v == target) {
      *hint_lo = mid;
      return mid;
    }
    if (v < target) lo = mid + 1;
    else hi = mid - 1;
  }
  if (lo <= hi0) *hint_lo = lo;
  return -1;
}

/* Helper: in-place column-pivoted Householder QR compression of an `m x k`
   column-major matrix `A` into `r x k` (`r = min(m, k)`) column-major `dst_fac`
   using caller-supplied scratch `jpvt[k]` and `vn[k]`. Returns 1 if all finite, 0 else. */
static int cdid_qr_compress_buf(
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

/* 9. Custom-IPW chunked observer & endpoint collector */
SEXP C_cdid_custom_ipw_collect(
    SEXP post_a_sexp,
    SEXP post_obs_sexp,
    SEXP bounds4_sexp,
    SEXP a_cohort_idx_sexp,
    SEXP w_a_sexp,
    SEXP Z_a_sexp,
    SEXP y_sexp,
    SEXP cmap_sexp,
    SEXP params_sexp) {
  if (TYPEOF(post_a_sexp) != INTSXP || Rf_isObject(post_a_sexp) ||
      TYPEOF(post_obs_sexp) != INTSXP || Rf_isObject(post_obs_sexp) ||
      TYPEOF(bounds4_sexp) != INTSXP || XLENGTH(bounds4_sexp) != 4 ||
      TYPEOF(a_cohort_idx_sexp) != INTSXP || Rf_isObject(a_cohort_idx_sexp) ||
      TYPEOF(w_a_sexp) != REALSXP || Rf_isObject(w_a_sexp) ||
      TYPEOF(Z_a_sexp) != REALSXP || Rf_isObject(Z_a_sexp) ||
      TYPEOF(y_sexp) != REALSXP || Rf_isObject(y_sexp) ||
      TYPEOF(cmap_sexp) != INTSXP ||
      TYPEOF(params_sexp) != INTSXP || XLENGTH(params_sexp) != 4) {
    return R_NilValue;
  }
  SEXP zdim = Rf_getAttrib(Z_a_sexp, R_DimSymbol);
  if (TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return R_NilValue;
  int n_a = INTEGER(params_sexp)[0];
  int n_obs = INTEGER(params_sexp)[1];
  int step = INTEGER(params_sexp)[2];
  int m = INTEGER(params_sexp)[3];
  int p = INTEGER(zdim)[1];
  int c_own = (int) XLENGTH(cmap_sexp);
  R_xlen_t n_post = XLENGTH(post_a_sexp);
  R_xlen_t n_y = XLENGTH(y_sexp);

  if (n_a < 1 || n_obs < 1 || step < 1 || m < 1 || p < 1 || c_own < 1 ||
      INTEGER(zdim)[0] != n_a || XLENGTH(a_cohort_idx_sexp) != (R_xlen_t) n_a ||
      XLENGTH(w_a_sexp) != (R_xlen_t) n_a || XLENGTH(post_obs_sexp) != n_post) {
    return R_NilValue;
  }

  const int *b4 = INTEGER(bounds4_sexp);
  R_xlen_t s_lo = (R_xlen_t) b4[0] - 1, s_hi = (R_xlen_t) b4[1] - 1;
  R_xlen_t t_lo = (R_xlen_t) b4[2] - 1, t_hi = (R_xlen_t) b4[3] - 1;
  if (s_lo < 0 || s_hi < s_lo || s_hi >= n_post || t_lo < 0 || t_hi < t_lo || t_hi >= n_post) {
    return R_NilValue;
  }

  const int *post_a = INTEGER(post_a_sexp);
  const int *post_obs = INTEGER(post_obs_sexp);
  const int *a_coh = INTEGER(a_cohort_idx_sexp);
  const double *w_a = REAL(w_a_sexp);
  const double *Z_a = REAL(Z_a_sexp);
  const double *y = REAL(y_sexp);
  const int *cmap = INTEGER(cmap_sexp);

  SEXP a_idx_s = PROTECT(Rf_allocVector(INTSXP, n_obs));
  SEXP gi_s = PROTECT(Rf_allocVector(INTSXP, n_obs));
  SEXP w_s = PROTECT(Rf_allocVector(REALSXP, n_obs));
  SEXP y0_s = PROTECT(Rf_allocVector(REALSXP, n_obs));
  SEXP y1_s = PROTECT(Rf_allocVector(REALSXP, n_obs));
  SEXP Z_use_s = PROTECT(Rf_allocMatrix(REALSXP, n_obs, p));
  SEXP zdn = Rf_getAttrib(Z_a_sexp, R_DimNamesSymbol);
  if (zdn != R_NilValue && TYPEOF(zdn) == VECSXP && XLENGTH(zdn) == 2) {
    SEXP cn = VECTOR_ELT(zdn, 1);
    if (cn != R_NilValue) {
      SEXP udn = PROTECT(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(udn, 0, R_NilValue);
      SET_VECTOR_ELT(udn, 1, cn);
      Rf_setAttrib(Z_use_s, R_DimNamesSymbol, udn);
      UNPROTECT(1);
    }
  }

  int *out_a = INTEGER(a_idx_s);
  int *out_gi = INTEGER(gi_s);
  double *out_w = REAL(w_s);
  double *out_y0 = REAL(y0_s);
  double *out_y1 = REAL(y1_s);
  double *out_Z = REAL(Z_use_s);

  R_xlen_t ps = s_lo, pt = t_lo;
  int write_pos = 0;
  int last_a = 0;
  int cur_s = 1;
  while (cur_s <= n_a) {
    R_CheckUserInterrupt();
    int64_t ce64 = (int64_t) cur_s + (int64_t) step - 1;
    int cur_e = (ce64 > (int64_t) n_a) ? n_a : (int) ce64;

    while (ps <= s_hi && post_a[ps] < cur_s) ++ps;
    while (pt <= t_hi && post_a[pt] < cur_s) ++pt;

    R_xlen_t is = ps, it = pt;
    while (is <= s_hi && it <= t_hi) {
      int vs = post_a[is], vt = post_a[it];
      if (vs > cur_e || vt > cur_e) break;
      if (vs == vt) {
        int a = vs;
        if (a <= last_a || a < cur_s || a > cur_e) {
          UNPROTECT(6);
          return R_NilValue;
        }
        int ci = a_coh[a - 1];
        if (ci >= 1 && ci <= c_own) {
          int gi = cmap[ci - 1];
          if (gi >= 1) {
            if (gi > m || write_pos >= n_obs) {
              UNPROTECT(6);
              return R_NilValue;
            }
            int os = post_obs[is], ot = post_obs[it];
            if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y) {
              UNPROTECT(6);
              return R_NilValue;
            }
            double yl = y[os - 1];
            double yr = y[ot - 1];
            double dy = yr - yl;
            if (!R_FINITE(yl) || !R_FINITE(yr) || !R_FINITE(dy)) {
              UNPROTECT(6);
              return R_NilValue;
            }
            out_a[write_pos] = a;
            out_gi[write_pos] = gi;
            out_w[write_pos] = w_a[a - 1];
            out_y0[write_pos] = yl;
            out_y1[write_pos] = yr;
            for (int col = 0; col < p; ++col) {
              out_Z[write_pos + (size_t) col * (size_t) n_obs] =
                  Z_a[(a - 1) + (size_t) col * (size_t) n_a];
            }
            ++write_pos;
            last_a = a;
          }
        }
        ++is;
        ++it;
      } else if (vs < vt) {
        ++is;
      } else {
        ++it;
      }
    }
    ps = is;
    pt = it;
    if (cur_e == n_a) break;
    cur_s = cur_e + 1;
  }

  if (write_pos != n_obs) {
    UNPROTECT(6);
    return R_NilValue;
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("risk_a_idx"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("risk_gi"));
  SET_STRING_ELT(nms_s, 2, Rf_mkChar("risk_w"));
  SET_STRING_ELT(nms_s, 3, Rf_mkChar("y0"));
  SET_STRING_ELT(nms_s, 4, Rf_mkChar("y1"));
  SET_STRING_ELT(nms_s, 5, Rf_mkChar("Z_use"));
  SET_VECTOR_ELT(out_s, 0, a_idx_s);
  SET_VECTOR_ELT(out_s, 1, gi_s);
  SET_VECTOR_ELT(out_s, 2, w_s);
  SET_VECTOR_ELT(out_s, 3, y0_s);
  SET_VECTOR_ELT(out_s, 4, y1_s);
  SET_VECTOR_ELT(out_s, 5, Z_use_s);
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
  UNPROTECT(8);
  return out_s;
}

/* 10. Custom-IPW multi-link IF chunk evaluator */
SEXP C_cdid_custom_ipw_if_chunk(
    SEXP a_range_sexp,
    SEXP u_w_sexp,
    SEXP c_chunk_sexp,
    SEXP cmap_mat_sexp,
    SEXP obs_idx_list_sexp,
    SEXP psi_list_sexp,
    SEXP rate_par_list_sexp,
    SEXP link_scales_sexp) {
  if (TYPEOF(a_range_sexp) != INTSXP || XLENGTH(a_range_sexp) != 2 ||
      TYPEOF(u_w_sexp) != REALSXP || TYPEOF(c_chunk_sexp) != INTSXP ||
      TYPEOF(cmap_mat_sexp) != INTSXP || TYPEOF(obs_idx_list_sexp) != VECSXP ||
      TYPEOF(psi_list_sexp) != VECSXP || TYPEOF(rate_par_list_sexp) != VECSXP ||
      TYPEOF(link_scales_sexp) != REALSXP) {
    return R_NilValue;
  }
  int a_s = INTEGER(a_range_sexp)[0];
  int a_e = INTEGER(a_range_sexp)[1];
  if (a_s < 1 || a_e < a_s) return R_NilValue;
  int b_len = a_e - a_s + 1;
  if (XLENGTH(u_w_sexp) != (R_xlen_t) b_len || XLENGTH(c_chunk_sexp) != (R_xlen_t) b_len) {
    return R_NilValue;
  }

  SEXP cdim = Rf_getAttrib(cmap_mat_sexp, R_DimSymbol);
  if (TYPEOF(cdim) != INTSXP || XLENGTH(cdim) != 2) return R_NilValue;
  int c_own = INTEGER(cdim)[0];
  int L = INTEGER(cdim)[1];
  if (c_own < 1 || L < 1 || XLENGTH(obs_idx_list_sexp) != (R_xlen_t) L ||
      XLENGTH(psi_list_sexp) != (R_xlen_t) L || XLENGTH(rate_par_list_sexp) != (R_xlen_t) L ||
      XLENGTH(link_scales_sexp) != (R_xlen_t) L) {
    return R_NilValue;
  }

  const double *u_w = REAL(u_w_sexp);
  const int *c_chunk = INTEGER(c_chunk_sexp);
  const int *cmap_mat = INTEGER(cmap_mat_sexp);
  const double *scales = REAL(link_scales_sexp);

  SEXP phi_s = PROTECT(Rf_allocMatrix(REALSXP, b_len, L));
  double *phi = REAL(phi_s);

  for (int ell = 0; ell < L; ++ell) {
    SEXP obs_s = VECTOR_ELT(obs_idx_list_sexp, ell);
    SEXP psi_s = VECTOR_ELT(psi_list_sexp, ell);
    SEXP rp_s = VECTOR_ELT(rate_par_list_sexp, ell);
    if (TYPEOF(obs_s) != INTSXP || TYPEOF(psi_s) != REALSXP || TYPEOF(rp_s) != REALSXP) {
      UNPROTECT(1);
      return R_NilValue;
    }
    int n_obs = (int) XLENGTH(obs_s);
    int rp_len = (int) XLENGTH(rp_s);
    if (XLENGTH(psi_s) != (R_xlen_t) n_obs || rp_len < 3 || (rp_len % 3) != 0) {
      UNPROTECT(1);
      return R_NilValue;
    }
    int m_ell = rp_len / 3;
    const int *obs_idx = INTEGER(obs_s);
    const double *psi_obs = REAL(psi_s);
    const double *b_rate = REAL(rp_s);
    const double *q_vec = b_rate + m_ell;
    const double *rates_vec = q_vec + m_ell;
    const int *cmap_col = cmap_mat + (size_t) ell * (size_t) c_own;
    double scale_ell = scales[ell];
    double *phi_col = phi + (size_t) ell * (size_t) b_len;

    /* Binary search first obs_idx >= a_s */
    int lo = 0, hi = n_obs - 1, ptr = n_obs;
    while (lo <= hi) {
      int mid = lo + ((hi - lo) >> 1);
      if (obs_idx[mid] >= a_s) {
        ptr = mid;
        hi = mid - 1;
      } else {
        lo = mid + 1;
      }
    }
    /* Verify strict sort on active window */
    for (int k = ptr; k < n_obs && obs_idx[k] <= a_e; ++k) {
      if (k > ptr && obs_idx[k] <= obs_idx[k - 1]) {
        UNPROTECT(1);
        return R_NilValue;
      }
    }

    for (int i = 0; i < b_len; ++i) {
      int a = a_s + i;
      double S_i = 0.0, psi_i = 0.0;
      if (ptr < n_obs && obs_idx[ptr] == a) {
        S_i = 1.0;
        psi_i = psi_obs[ptr];
        ++ptr;
      }
      int ci = c_chunk[i];
      if (ci >= 1 && ci <= c_own) {
        int gk = cmap_col[ci - 1];
        if (gk >= 1 && gk <= m_ell) {
          double base_val = scale_ell * psi_i;
          double rho_val = u_w[i] * (S_i - rates_vec[gk - 1]) / q_vec[gk - 1];
          double rate_val = rho_val * b_rate[gk - 1];
          double v = base_val + rate_val;
          if (!R_FINITE(v)) {
            UNPROTECT(1);
            return R_NilValue;
          }
          phi_col[i] = v;
        } else {
          phi_col[i] = 0.0;
        }
      } else {
        phi_col[i] = 0.0;
      }
    }
  }

  UNPROTECT(1);
  return phi_s;
}

/* 11. Raw-IPW pair summary + 4-double mask fingerprint scanner */
SEXP C_cdid_ipw_raw_pair_summary(
    SEXP post_a_sexp,
    SEXP bounds4_sexp,
    SEXP a_cohort_idx_sexp,
    SEXP w_a_sexp,
    SEXP params_sexp) {
  if (TYPEOF(post_a_sexp) != INTSXP || Rf_isObject(post_a_sexp) ||
      TYPEOF(a_cohort_idx_sexp) != INTSXP || Rf_isObject(a_cohort_idx_sexp) ||
      TYPEOF(w_a_sexp) != REALSXP || Rf_isObject(w_a_sexp) ||
      TYPEOF(params_sexp) != INTSXP || XLENGTH(params_sexp) != 4) {
    return R_NilValue;
  }
  int c_a = INTEGER(params_sexp)[0];
  int n_a = INTEGER(params_sexp)[1];
  int step = INTEGER(params_sexp)[2];
  int complete_pair = INTEGER(params_sexp)[3];
  if (c_a < 1 || n_a < 1 || step < 1 ||
      XLENGTH(a_cohort_idx_sexp) != (R_xlen_t) n_a || XLENGTH(w_a_sexp) != (R_xlen_t) n_a) {
    return R_NilValue;
  }

  int has_pair = 0;
  R_xlen_t s_lo = 0, s_hi = -1, t_lo = 0, t_hi = -1;
  R_xlen_t n_post = XLENGTH(post_a_sexp);
  if (bounds4_sexp != R_NilValue) {
    if (TYPEOF(bounds4_sexp) != INTSXP || XLENGTH(bounds4_sexp) != 4) return R_NilValue;
    const int *b4 = INTEGER(bounds4_sexp);
    s_lo = (R_xlen_t) b4[0] - 1;
    s_hi = (R_xlen_t) b4[1] - 1;
    t_lo = (R_xlen_t) b4[2] - 1;
    t_hi = (R_xlen_t) b4[3] - 1;
    if (s_lo < 0 || s_hi < s_lo || s_hi >= n_post || t_lo < 0 || t_hi < t_lo || t_hi >= n_post) {
      return R_NilValue;
    }
    has_pair = 1;
  }

  SEXP cnt_s = PROTECT(Rf_allocVector(INTSXP, c_a));
  SEXP mass_s = PROTECT(Rf_allocVector(REALSXP, c_a));
  SEXP fp1_s = PROTECT(Rf_allocVector(REALSXP, c_a));
  SEXP fp2_s = PROTECT(Rf_allocVector(REALSXP, c_a));
  SEXP fp3_s = PROTECT(Rf_allocVector(REALSXP, c_a));
  int *cnt = INTEGER(cnt_s);
  double *mass = REAL(mass_s);
  double *fp1 = REAL(fp1_s);
  double *fp2 = REAL(fp2_s);
  double *fp3 = REAL(fp3_s);
  memset(cnt, 0, (size_t) c_a * sizeof(int));
  memset(mass, 0, (size_t) c_a * sizeof(double));
  memset(fp1, 0, (size_t) c_a * sizeof(double));
  memset(fp2, 0, (size_t) c_a * sizeof(double));
  memset(fp3, 0, (size_t) c_a * sizeof(double));

  int n_obs_all = 0;
  int chunks = 0;
  double ops = 0.0;
  const int64_t P = 67108859LL;

  if (has_pair) {
    const void *vmax = vmaxget();
    long double *ch_mass = R_allocLD((size_t) c_a);
    int *ch_touched = (int *) R_alloc((size_t) c_a, sizeof(int));
    int64_t *s1 = (int64_t *) R_alloc((size_t) c_a, sizeof(int64_t));
    int64_t *s2 = (int64_t *) R_alloc((size_t) c_a, sizeof(int64_t));
    int64_t *s3 = (int64_t *) R_alloc((size_t) c_a, sizeof(int64_t));
    memset(s1, 0, (size_t) c_a * sizeof(int64_t));
    memset(s2, 0, (size_t) c_a * sizeof(int64_t));
    memset(s3, 0, (size_t) c_a * sizeof(int64_t));

    const int *post_a = INTEGER(post_a_sexp);
    const int *a_coh = INTEGER(a_cohort_idx_sexp);
    const double *w_a = REAL(w_a_sexp);

    if (complete_pair) {
      if ((s_hi - s_lo + 1) < (R_xlen_t) n_a || (t_hi - t_lo + 1) < (R_xlen_t) n_a) {
        vmaxset(vmax);
        UNPROTECT(5);
        return R_NilValue;
      }
    }

    R_xlen_t ps = s_lo, pt = t_lo;
    int a_start = 1;
    while (a_start <= n_a) {
      R_CheckUserInterrupt();
      int64_t ae64 = (int64_t) a_start + (int64_t) step - 1;
      int a_end = (ae64 > (int64_t) n_a) ? n_a : (int) ae64;
      for (int j = 0; j < c_a; ++j) {
        ch_mass[j] = 0.0L;
        ch_touched[j] = 0;
      }

      if (complete_pair) {
        for (int a = a_start; a <= a_end; ++a) {
          if (post_a[s_lo + (a - 1)] != a || post_a[t_lo + (a - 1)] != a) {
            vmaxset(vmax);
            UNPROTECT(5);
            return R_NilValue;
          }
          ++n_obs_all;
          int j = a_coh[a - 1];
          if (j >= 1 && j <= c_a) {
            int j0 = j - 1;
            ++cnt[j0];
            ch_mass[j0] += (long double) w_a[a - 1];
            ch_touched[j0] = 1;
            int64_t u = ((int64_t) a) % P;
            int64_t u2 = (u * u) % P;
            int64_t u3 = (u2 * u) % P;
            s1[j0] = (s1[j0] + u) % P;
            s2[j0] = (s2[j0] + u2) % P;
            s3[j0] = (s3[j0] + u3) % P;
          }
        }
      } else {
        while (ps <= s_hi && post_a[ps] < a_start) ++ps;
        while (pt <= t_hi && post_a[pt] < a_start) ++pt;
        R_xlen_t is = ps, it = pt;
        while (is <= s_hi && it <= t_hi) {
          int vs = post_a[is], vt = post_a[it];
          if (vs > a_end || vt > a_end) break;
          if (vs == vt) {
            int a = vs;
            ++n_obs_all;
            int j = a_coh[a - 1];
            if (j >= 1 && j <= c_a) {
              int j0 = j - 1;
              ++cnt[j0];
              ch_mass[j0] += (long double) w_a[a - 1];
              ch_touched[j0] = 1;
              int64_t u = ((int64_t) a) % P;
              int64_t u2 = (u * u) % P;
              int64_t u3 = (u2 * u) % P;
              s1[j0] = (s1[j0] + u) % P;
              s2[j0] = (s2[j0] + u2) % P;
              s3[j0] = (s3[j0] + u3) % P;
            }
            ++is;
            ++it;
          } else if (vs < vt) {
            ++is;
          } else {
            ++it;
          }
        }
        ps = is;
        pt = it;
      }

      for (int j0 = 0; j0 < c_a; ++j0) {
        if (ch_touched[j0]) mass[j0] += (double) ch_mass[j0];
      }
      ++chunks;
      ops += 4.0 * (double) (a_end - a_start + 1);
      if (a_end == n_a) break;
      a_start = a_end + 1;
    }

    for (int j0 = 0; j0 < c_a; ++j0) {
      fp1[j0] = (double) s1[j0];
      fp2[j0] = (double) s2[j0];
      fp3[j0] = (double) s3[j0];
    }
    vmaxset(vmax);
  } else {
    chunks = (int) ceil((double) n_a / (double) step);
    ops = (double) n_a;
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 8));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 8));
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("n_observed_pair_all"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("cohort_pair_count"));
  SET_STRING_ELT(nms_s, 2, Rf_mkChar("cohort_pair_mass"));
  SET_STRING_ELT(nms_s, 3, Rf_mkChar("cohort_fp_s1"));
  SET_STRING_ELT(nms_s, 4, Rf_mkChar("cohort_fp_s2"));
  SET_STRING_ELT(nms_s, 5, Rf_mkChar("cohort_fp_s3"));
  SET_STRING_ELT(nms_s, 6, Rf_mkChar("chunks"));
  SET_STRING_ELT(nms_s, 7, Rf_mkChar("ops"));
  SET_VECTOR_ELT(out_s, 0, Rf_ScalarInteger(n_obs_all));
  SET_VECTOR_ELT(out_s, 1, cnt_s);
  SET_VECTOR_ELT(out_s, 2, mass_s);
  SET_VECTOR_ELT(out_s, 3, fp1_s);
  SET_VECTOR_ELT(out_s, 4, fp2_s);
  SET_VECTOR_ELT(out_s, 5, fp3_s);
  SET_VECTOR_ELT(out_s, 6, Rf_ScalarInteger(chunks));
  SET_VECTOR_ELT(out_s, 7, Rf_ScalarReal(ops));
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);
  UNPROTECT(7);
  return out_s;
}

/* 12. Raw-IPW mask equality check (op=1) and shared-mask endpoint collector (op=2) */
SEXP C_cdid_ipw_raw_masks_ops(
    SEXP op_sexp,
    SEXP post_a_sexp,
    SEXP post_obs_sexp,
    SEXP bounds_sexp,
    SEXP a_cohort_idx_sexp,
    SEXP in_coh_sexp,
    SEXP extra_sexp) {
  int op = Rf_asInteger(op_sexp);
  if ((op != 1 && op != 2) || TYPEOF(post_a_sexp) != INTSXP || Rf_isObject(post_a_sexp) ||
      TYPEOF(bounds_sexp) != INTSXP || TYPEOF(a_cohort_idx_sexp) != INTSXP || Rf_isObject(a_cohort_idx_sexp) ||
      TYPEOF(in_coh_sexp) != INTSXP) {
    return R_NilValue;
  }
  const int *post_a = INTEGER(post_a_sexp);
  const int *a_coh = INTEGER(a_cohort_idx_sexp);
  const int *in_coh = INTEGER(in_coh_sexp);
  R_xlen_t n_post = XLENGTH(post_a_sexp);
  int n_a = (int) XLENGTH(a_cohort_idx_sexp);
  int c_a = (int) XLENGTH(in_coh_sexp);
  if (n_a < 1 || c_a < 1) return R_NilValue;

  if (op == 1) {
    if (XLENGTH(bounds_sexp) != 8 || TYPEOF(extra_sexp) != INTSXP || XLENGTH(extra_sexp) != 2) {
      return R_NilValue;
    }
    const int *b = INTEGER(bounds_sexp);
    R_xlen_t s1 = (R_xlen_t) b[0] - 1, e1 = (R_xlen_t) b[1] - 1;
    R_xlen_t s2 = (R_xlen_t) b[2] - 1, e2 = (R_xlen_t) b[3] - 1;
    R_xlen_t s3 = (R_xlen_t) b[4] - 1, e3 = (R_xlen_t) b[5] - 1;
    R_xlen_t s4 = (R_xlen_t) b[6] - 1, e4 = (R_xlen_t) b[7] - 1;

    int ok1 = (b[0] >= 1 && e1 >= s1 && e1 < n_post && b[2] >= 1 && e2 >= s2 && e2 < n_post);
    int ok2 = (b[4] >= 1 && e3 >= s3 && e3 < n_post && b[6] >= 1 && e4 >= s4 && e4 < n_post);
    if (!ok1) { s1 = 0; e1 = -1; s2 = 0; e2 = -1; }
    if (!ok2) { s3 = 0; e3 = -1; s4 = 0; e4 = -1; }

    R_xlen_t i1 = s1, j1 = s2, i2 = s3, j2 = s4;
    for (;;) {
      int a1 = 0;
      while (i1 <= e1 && j1 <= e2) {
        int u = post_a[i1], v = post_a[j1];
        if (u == v) {
          if (u >= 1 && u <= n_a) {
            int c = a_coh[u - 1];
            if (c >= 1 && c <= c_a && in_coh[c - 1]) {
              a1 = u;
              ++i1;
              ++j1;
              break;
            }
          }
          ++i1;
          ++j1;
        } else if (u < v) ++i1;
        else ++j1;
      }

      int a2 = 0;
      while (i2 <= e3 && j2 <= e4) {
        int u = post_a[i2], v = post_a[j2];
        if (u == v) {
          if (u >= 1 && u <= n_a) {
            int c = a_coh[u - 1];
            if (c >= 1 && c <= c_a && in_coh[c - 1]) {
              a2 = u;
              ++i2;
              ++j2;
              break;
            }
          }
          ++i2;
          ++j2;
        } else if (u < v) ++i2;
        else ++j2;
      }

      if (a1 != a2) return Rf_ScalarLogical(0);
      if (a1 == 0) break;
    }
    return Rf_ScalarLogical(1);
  }

  /* op == 2: collect pair endpoints matching expected_a_idx */
  if (TYPEOF(post_obs_sexp) != INTSXP || Rf_isObject(post_obs_sexp) || XLENGTH(post_obs_sexp) != n_post ||
      XLENGTH(bounds_sexp) != 4 || TYPEOF(extra_sexp) != INTSXP || Rf_isObject(extra_sexp)) {
    return R_NilValue;
  }
  const int *b4 = INTEGER(bounds_sexp);
  R_xlen_t s_lo = (R_xlen_t) b4[0] - 1, s_hi = (R_xlen_t) b4[1] - 1;
  R_xlen_t t_lo = (R_xlen_t) b4[2] - 1, t_hi = (R_xlen_t) b4[3] - 1;
  if (b4[0] < 1 || s_hi < s_lo || s_hi >= n_post || b4[2] < 1 || t_hi < t_lo || t_hi >= n_post) {
    return R_NilValue;
  }
  int n_obs = (int) XLENGTH(extra_sexp);
  const int *expected_a = INTEGER(extra_sexp);
  const int *post_obs = INTEGER(post_obs_sexp);

  SEXP obs_s_s = PROTECT(Rf_allocVector(INTSXP, n_obs));
  SEXP obs_t_s = PROTECT(Rf_allocVector(INTSXP, n_obs));
  int *out_s = INTEGER(obs_s_s);
  int *out_t = INTEGER(obs_t_s);

  int pos = 0;
  R_xlen_t is = s_lo, it = t_lo;
  while (is <= s_hi && it <= t_hi) {
    int vs = post_a[is], vt = post_a[it];
    if (vs == vt) {
      if (vs >= 1 && vs <= n_a) {
        int c = a_coh[vs - 1];
        if (c >= 1 && c <= c_a && in_coh[c - 1]) {
          if (pos >= n_obs || expected_a[pos] != vs) {
            UNPROTECT(2);
            return R_NilValue;
          }
          out_s[pos] = post_obs[is];
          out_t[pos] = post_obs[it];
          ++pos;
        }
      }
      ++is;
      ++it;
    } else if (vs < vt) ++is;
    else ++it;
  }

  if (pos != n_obs) {
    UNPROTECT(2);
    return R_NilValue;
  }

  SEXP res_s = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("obs_s"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("obs_t"));
  SET_VECTOR_ELT(res_s, 0, obs_s_s);
  SET_VECTOR_ELT(res_s, 1, obs_t_s);
  Rf_setAttrib(res_s, R_NamesSymbol, nms_s);
  UNPROTECT(4);
  return res_s;
}

static R_INLINE SEXP cdid_list_get_named(SEXP lst, const char *target_name) {
  if (TYPEOF(lst) != VECSXP) return R_UnboundValue;
  SEXP nms = Rf_getAttrib(lst, R_NamesSymbol);
  if (TYPEOF(nms) != STRSXP) return R_UnboundValue;
  R_xlen_t n = XLENGTH(lst);
  if (XLENGTH(nms) != n) return R_UnboundValue;
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP s = STRING_ELT(nms, i);
    if (s != NA_STRING && strcmp(CHAR(s), target_name) == 0) {
      return VECTOR_ELT(lst, i);
    }
  }
  return R_UnboundValue;
}

/* 13. Fast weighted binomial IRLS logit for ipw_nuisance_observer_fit_worker */
SEXP C_cdid_ipw_glm_fit(
    SEXP Z_sexp,
    SEXP D_sexp,
    SEXP w_sexp,
    SEXP ctrl_sexp) {
  if (TYPEOF(Z_sexp) != REALSXP || Rf_isObject(Z_sexp)) return R_NilValue;
  SEXP zdim = Rf_getAttrib(Z_sexp, R_DimSymbol);
  if (TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return R_NilValue;
  int n = INTEGER(zdim)[0];
  int p = INTEGER(zdim)[1];
  if (p < 1 || p > 128 || n < p) return R_NilValue;

  int full_worker = (TYPEOF(D_sexp) == INTSXP && TYPEOF(ctrl_sexp) == VECSXP && XLENGTH(ctrl_sexp) == 3);
  int nprot = 0;
  SEXP sup_s = R_NilValue, risk_a_idx_s = R_NilValue;
  SEXP n_s = R_NilValue, m_coh_s = R_NilValue, group_mass_s = R_NilValue;
  SEXP pair_count_s = R_NilValue, pair_mass_s = R_NilValue, rates_s = R_NilValue;
  SEXP D_use_s = R_NilValue, fit_w_s = R_NilValue;
  int m_cohorts = 0;
  double n_val = 0.0;
  double eps = 0.0;
  int maxit = 0;
  const double *Z = REAL(Z_sexp);
  const double *D = NULL;
  const double *w = NULL;
  int n_pos_w = 0;
  long double sum_w_ld = 0.0L;
  long double devold_ld = 0.0L;
  const double ln2_val = log(2.0);
  double sum_w = 0.0;

  if (full_worker) {
    if (Rf_isObject(D_sexp) || TYPEOF(w_sexp) != REALSXP || Rf_isObject(w_sexp) ||
        XLENGTH(D_sexp) != (R_xlen_t) n || XLENGTH(w_sexp) != (R_xlen_t) n) {
      return R_NilValue;
    }
    sup_s = VECTOR_ELT(ctrl_sexp, 0);
    risk_a_idx_s = VECTOR_ELT(ctrl_sexp, 1);
    SEXP cpar_s = VECTOR_ELT(ctrl_sexp, 2);
    if (TYPEOF(sup_s) != VECSXP || TYPEOF(cpar_s) != REALSXP || XLENGTH(cpar_s) != 2) {
      return R_NilValue;
    }
    eps = REAL(cpar_s)[0];
    maxit = (int) REAL(cpar_s)[1];
    if (!R_FINITE(eps) || eps <= 0.0 || maxit < 1) return R_NilValue;

    n_s = cdid_list_get_named(sup_s, "n");
    m_coh_s = cdid_list_get_named(sup_s, "m_cohorts");
    group_mass_s = cdid_list_get_named(sup_s, "group_mass");
    if (n_s == R_UnboundValue || m_coh_s == R_UnboundValue || group_mass_s == R_UnboundValue ||
        TYPEOF(group_mass_s) != REALSXP) {
      return R_NilValue;
    }
    n_val = Rf_asReal(n_s);
    m_cohorts = Rf_asInteger(m_coh_s);
    if (!R_FINITE(n_val) || n_val <= 0.0 || m_cohorts < 2 || XLENGTH(group_mass_s) != (R_xlen_t) m_cohorts) {
      return R_NilValue;
    }
    const double *gmass = REAL(group_mass_s);
    const int *rgi = INTEGER(D_sexp);
    const double *rw = REAL(w_sexp);

    pair_count_s = PROTECT(Rf_allocVector(INTSXP, m_cohorts));
    pair_mass_s = PROTECT(Rf_allocVector(REALSXP, m_cohorts));
    rates_s = PROTECT(Rf_allocVector(REALSXP, m_cohorts));
    nprot += 3;
    int *pcount = INTEGER(pair_count_s);
    double *pmass = REAL(pair_mass_s);
    double *rates = REAL(rates_s);
    memset(pcount, 0, (size_t) m_cohorts * sizeof(int));
    const void *vm_pm = vmaxget();
    long double *pmass_ld = R_allocLD((size_t) m_cohorts);
    for (int j = 0; j < m_cohorts; ++j) pmass_ld[j] = 0.0L;
    for (int i = 0; i < n; ++i) {
      int g = rgi[i];
      double rwi = rw[i];
      if (g < 1 || g > m_cohorts || !R_FINITE(rwi)) {
        vmaxset(vm_pm);
        UNPROTECT(nprot);
        return R_NilValue;
      }
      pcount[g - 1] += 1;
      pmass_ld[g - 1] += (long double) rwi;
    }
    for (int j = 0; j < m_cohorts; ++j) {
      pmass[j] = (double) pmass_ld[j];
    }
    vmaxset(vm_pm);
    for (int j = 0; j < m_cohorts; ++j) {
      if (pcount[j] < 1 || !R_FINITE(pmass[j]) || pmass[j] <= 0.0 ||
          !R_FINITE(gmass[j]) || gmass[j] <= 0.0) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
      double rj = pmass[j] / gmass[j];
      if (!R_FINITE(rj) || rj <= 0.0 || rj > 1.0) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
      rates[j] = rj;
    }
    SEXP gm_nms = Rf_getAttrib(group_mass_s, R_NamesSymbol);
    if (gm_nms != R_NilValue) Rf_setAttrib(rates_s, R_NamesSymbol, gm_nms);

    D_use_s = PROTECT(Rf_allocVector(REALSXP, n));
    fit_w_s = PROTECT(Rf_allocVector(REALSXP, n));
    nprot += 2;
    double *D_buf = REAL(D_use_s);
    double *w_buf = REAL(fit_w_s);
    for (int i = 0; i < n; ++i) {
      int g = rgi[i];
      double rc = (rates[g - 1] > 1e-6) ? rates[g - 1] : 1e-6;
      double wi = rw[i] / rc;
      if (!R_FINITE(wi) || wi <= 0.0) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
      D_buf[i] = (g == 1) ? 1.0 : 0.0;
      w_buf[i] = wi;
      sum_w_ld += (long double) wi;
      devold_ld += (long double) (2.0 * wi * ln2_val);
    }
    sum_w = (double) sum_w_ld;
    if (!R_FINITE(sum_w) || sum_w <= 0.0) {
      UNPROTECT(nprot);
      return R_NilValue;
    }
    n_pos_w = n;
    D = D_buf;
    w = w_buf;

    /* Check finite Z, exact alias, and coordinate separation up front */
    for (int j = 0; j < p; ++j) {
      const double *zj = Z + (size_t) j * (size_t) n;
      int all_zero = 1;
      int has_pos = 0, has_neg = 0;
      for (int i = 0; i < n; ++i) {
        double v = zj[i];
        if (!R_FINITE(v)) {
          UNPROTECT(nprot);
          return R_NilValue;
        }
        if (v != 0.0) {
          all_zero = 0;
          double sgn = (D[i] == 1.0) ? v : -v;
          if (sgn > 0.0) has_pos = 1;
          else if (sgn < 0.0) has_neg = 1;
        }
      }
      if (all_zero || (has_pos && !has_neg) || (has_neg && !has_pos)) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
      for (int k = 0; k < j; ++k) {
        const double *zk = Z + (size_t) k * (size_t) n;
        int eq = 1;
        for (int i = 0; i < n; ++i) {
          if (zj[i] != zk[i]) {
            eq = 0;
            break;
          }
        }
        if (eq) {
          UNPROTECT(nprot);
          return R_NilValue;
        }
      }
    }
  } else {
    if (TYPEOF(D_sexp) != REALSXP || Rf_isObject(D_sexp) ||
        TYPEOF(w_sexp) != REALSXP || Rf_isObject(w_sexp) ||
        TYPEOF(ctrl_sexp) != REALSXP || XLENGTH(ctrl_sexp) != 2 ||
        XLENGTH(D_sexp) != (R_xlen_t) n || XLENGTH(w_sexp) != (R_xlen_t) n) {
      return R_NilValue;
    }
    eps = REAL(ctrl_sexp)[0];
    maxit = (int) REAL(ctrl_sexp)[1];
    if (!R_FINITE(eps) || eps <= 0.0 || maxit < 1) return R_NilValue;

    D = REAL(D_sexp);
    w = REAL(w_sexp);
    for (int i = 0; i < n; ++i) {
      double wi = w[i], di = D[i];
      if (!R_FINITE(wi) || wi <= 0.0 || (di != 0.0 && di != 1.0)) return R_NilValue;
      ++n_pos_w;
      sum_w_ld += (long double) wi;
      devold_ld += (long double) (2.0 * wi * ln2_val);
    }
    sum_w = (double) sum_w_ld;
    if (n_pos_w < p || !R_FINITE(sum_w) || sum_w <= 0.0) return R_NilValue;

    /* Check finite Z and coordinate separation up front */
    for (int j = 0; j < p; ++j) {
      const double *zj = Z + (size_t) j * (size_t) n;
      int has_pos = 0, has_neg = 0;
      for (int i = 0; i < n; ++i) {
        double v = zj[i];
        if (!R_FINITE(v)) return R_NilValue;
        if (v != 0.0) {
          double sgn = (D[i] == 1.0) ? v : -v;
          if (sgn > 0.0) has_pos = 1;
          else if (sgn < 0.0) has_neg = 1;
        }
      }
      if ((has_pos && !has_neg) || (has_neg && !has_pos)) return R_NilValue;
    }
  }

  const void *vmax = vmaxget();
  double *eta = (double *) R_alloc((size_t) n, sizeof(double));
  double *mu = (double *) R_alloc((size_t) n, sizeof(double));
  double *wt = (double *) R_alloc((size_t) n, sizeof(double));
  double *wz = (double *) R_alloc((size_t) n, sizeof(double));
  double *W_Z = (double *) R_alloc((size_t) n * (size_t) p, sizeof(double));
  double *R_qr = (double *) R_alloc((size_t) p * (size_t) p, sizeof(double));
  double *beta = (double *) R_alloc((size_t) p, sizeof(double));

  /* Match ipw_glm_fit_internal: start = rep(0, p), intercept = FALSE */
  for (int i = 0; i < n; ++i) {
    eta[i] = 0.0;
    mu[i] = 0.5;
  }
  for (int j = 0; j < p; ++j) beta[j] = 0.0;

  double devold = (double) devold_ld;
  int conv = 0;
  int iter_used = 0;
  double dev = devold;
  const double bnd_eps = 10.0 * DBL_EPSILON;
  const double qr_tol = (eps / 1000.0 < 1e-7) ? (eps / 1000.0) : 1e-7;
  int n_b = n, p_b = p, one_b = 1;
  double one_d = 1.0, zero_d = 0.0;

  for (int iter = 1; iter <= maxit; ++iter) {
    iter_used = iter;
    for (int i = 0; i < n; ++i) {
      double wi = w[i];
      double mui = mu[i];
      double varmu = mui * (1.0 - mui);
      double opexp = exp(eta[i]);
      double mu_eta = opexp / ((1.0 + opexp) * (1.0 + opexp));
      if (!(varmu > 0.0) || !(mu_eta > 0.0) || !R_FINITE(mu_eta)) {
        vmaxset(vmax);
        if (nprot > 0) UNPROTECT(nprot);
        return R_NilValue;
      }
      double w_i = sqrt((wi * (mu_eta * mu_eta)) / varmu);
      double z_i = eta[i] + (D[i] - mui) / mu_eta;
      wt[i] = w_i;
      wz[i] = z_i * w_i;
    }

    for (int j = 0; j < p; ++j) {
      const double *zj = Z + (size_t) j * (size_t) n;
      double *wzj = W_Z + (size_t) j * (size_t) n;
      for (int i = 0; i < n; ++i) wzj[i] = zj[i] * wt[i];
    }

    if (!cdid_qr_wls_solve(W_Z, wz, R_qr, beta, n, p, qr_tol, 0.0, 1)) {
      vmaxset(vmax);
      if (nprot > 0) UNPROTECT(nprot);
      return R_NilValue;
    }

    F77_CALL(dgemm)("N", "N", &n_b, &one_b, &p_b, &one_d,
                    Z, &n_b, beta, &p_b, &zero_d, eta, &n_b FCONE FCONE);

    long double dev_ld = 0.0L;
    for (int i = 0; i < n; ++i) {
      double ei = eta[i];
      if (!R_FINITE(ei) || fabs(ei) >= 30.0) {
        vmaxset(vmax);
        if (nprot > 0) UNPROTECT(nprot);
        return R_NilValue;
      }
      double opexp = exp(ei);
      double mui = opexp / (1.0 + opexp);
      if (mui <= bnd_eps || mui >= 1.0 - bnd_eps) {
        vmaxset(vmax);
        if (nprot > 0) UNPROTECT(nprot);
        return R_NilValue;
      }
      mu[i] = mui;
      double loss = (D[i] == 1.0) ? log(1.0 / mui) : log(1.0 / (1.0 - mui));
      dev_ld += (long double) (2.0 * w[i] * loss);
    }
    dev = (double) dev_ld;
    if (!R_FINITE(dev)) {
      vmaxset(vmax);
      if (nprot > 0) UNPROTECT(nprot);
      return R_NilValue;
    }
    if (fabs(dev - devold) / (0.1 + fabs(dev)) < eps) {
      conv = 1;
      break;
    }
    devold = dev;
  }

  if (!conv) {
    vmaxset(vmax);
    if (nprot > 0) UNPROTECT(nprot);
    return R_NilValue;
  }

  /* intercept = FALSE -> wtdmu = linkinv(0) = 0.5 -> null_deviance = sum(2 * w * log(2)) */
  double null_dev = (double) devold_ld;

  long double aic_ll_ld = 0.0L;
  for (int i = 0; i < n; ++i) {
    if (w[i] > 0.0) {
      double wr = w[i] / w[i];
      double term = wr * Rf_dbinom(Rf_fround(w[i] * D[i], 0.0), Rf_fround(w[i], 0.0), mu[i], 1);
      aic_ll_ld += (long double) term;
    }
  }
  double aic = -2.0 * ((double) aic_ll_ld) + 2.0 * (double) p;

  SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, p));
  memcpy(REAL(beta_s), beta, (size_t) p * sizeof(double));
  SEXP qr_R_s = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  memcpy(REAL(qr_R_s), R_qr, (size_t) p * (size_t) p * sizeof(double));
  nprot += 2;
  SEXP zdn = Rf_getAttrib(Z_sexp, R_DimNamesSymbol);
  SEXP rn_s = R_NilValue;
  if (zdn != R_NilValue && TYPEOF(zdn) == VECSXP && XLENGTH(zdn) == 2) {
    SEXP r0 = VECTOR_ELT(zdn, 0);
    SEXP cn = VECTOR_ELT(zdn, 1);
    if (r0 != R_NilValue && TYPEOF(r0) == STRSXP && XLENGTH(r0) == (R_xlen_t) n) {
      rn_s = r0;
    }
    if (cn != R_NilValue) {
      Rf_setAttrib(beta_s, R_NamesSymbol, cn);
    }
    if (rn_s != R_NilValue || cn != R_NilValue) {
      SEXP rdn = PROTECT(Rf_allocVector(VECSXP, 2));
      if (rn_s != R_NilValue) {
        SEXP rn_p = PROTECT(Rf_allocVector(STRSXP, p));
        for (int j = 0; j < p; ++j) SET_STRING_ELT(rn_p, j, STRING_ELT(rn_s, j));
        SET_VECTOR_ELT(rdn, 0, rn_p);
        UNPROTECT(1);
      } else {
        SET_VECTOR_ELT(rdn, 0, R_NilValue);
      }
      SET_VECTOR_ELT(rdn, 1, cn);
      Rf_setAttrib(qr_R_s, R_DimNamesSymbol, rdn);
      UNPROTECT(1);
    }
  }

  if (full_worker) {
    const double clip_cap = 1.0 - 1e-6;
    const double clip_cap_odds = clip_cap / (1.0 - clip_cap);
    int n_t = INTEGER(pair_count_s)[0];
    int n_c = n - n_t;
    SEXP odds_s = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP prob_s = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP is_t_s = PROTECT(Rf_allocVector(LGLSXP, n));
    SEXP is_c_s = PROTECT(Rf_allocVector(LGLSXP, n));
    SEXP a_obs_s = PROTECT(Rf_allocVector(REALSXP, n_t));
    SEXP b_obs_s = PROTECT(Rf_allocVector(REALSXP, n_c));
    nprot += 6;
    if (rn_s != R_NilValue) {
      Rf_setAttrib(odds_s, R_NamesSymbol, rn_s);
      Rf_setAttrib(prob_s, R_NamesSymbol, rn_s);
    }
    double *odds = REAL(odds_s);
    double *prob = REAL(prob_s);
    int *is_t = LOGICAL(is_t_s);
    int *is_c = LOGICAL(is_c_s);
    double *a_obs = REAL(a_obs_s);
    double *b_obs = REAL(b_obs_s);
    const int *rgi = INTEGER(D_sexp);

    int it = 0, ic = 0;
    long double sum_a_ld = 0.0L, sum_b_ld = 0.0L;
    double p_min = R_PosInf, p_max = R_NegInf, w_max = R_NegInf;
    int p_capped = 0, p_capped_ctrl = 0;

    for (int i = 0; i < n; ++i) {
      double e = eta[i];
      double od = exp(e > 40.0 ? 0.0 : e);
      double pr = (e > 40.0) ? (1.0 / (1.0 + exp(-e))) : (od / (1.0 + od));
      int capped = !(pr < clip_cap);
      if (capped) {
        od = clip_cap_odds;
        ++p_capped;
      }
      if (!R_FINITE(od) || !(od > 0.0) || !R_FINITE(pr) ||
          pr <= bnd_eps || pr >= 1.0 - bnd_eps || !(pr * (1.0 - pr) > 0.0)) {
        UNPROTECT(nprot);
        vmaxset(vmax);
        return R_NilValue;
      }
      odds[i] = od;
      prob[i] = pr;
      if (pr < p_min) p_min = pr;
      if (pr > p_max) p_max = pr;
      if (w[i] > w_max) w_max = w[i];

      if (rgi[i] == 1) {
        is_t[i] = 1;
        is_c[i] = 0;
        a_obs[it++] = w[i];
        sum_a_ld += (long double) w[i];
      } else {
        is_t[i] = 0;
        is_c[i] = 1;
        double bi = w[i] * od;
        b_obs[ic++] = bi;
        sum_b_ld += (long double) bi;
        if (capped) ++p_capped_ctrl;
      }
    }
    double sum_a = (double) sum_a_ld;
    double sum_b = (double) sum_b_ld;
    double mean_a = sum_a / n_val;
    double mean_b = sum_b / n_val;
    if (!R_FINITE(sum_a) || !(sum_a > 0.0) || !R_FINITE(mean_a) || !(mean_a > 0.0) ||
        !R_FINITE(sum_b) || !(sum_b > 0.0) || !R_FINITE(mean_b) || !(mean_b > 0.0)) {
      UNPROTECT(nprot);
      vmaxset(vmax);
      return R_NilValue;
    }

    const double *rates = REAL(rates_s);
    double r_min = rates[0], r_max = rates[0];
    int r_floored = 0;
    for (int j = 0; j < m_cohorts; ++j) {
      if (rates[j] < r_min) r_min = rates[j];
      if (rates[j] > r_max) r_max = rates[j];
      if (rates[j] <= 1e-6) ++r_floored;
    }
    double p_clip_max = (p_max < clip_cap) ? p_max : clip_cap;
    double r_floor_min = (r_min > 1e-6) ? r_min : 1e-6;

    SEXP clip_s = PROTECT(Rf_allocVector(VECSXP, 13));
    SEXP clip_nms = PROTECT(Rf_allocVector(STRSXP, 13));
    nprot += 2;
    const char *cnms[13] = {
      "eps", "p_min", "p_max", "p_clip_max", "p_n", "p_capped", "p_capped_ctrl",
      "r_min", "r_max", "r_floor_min", "r_n", "r_floored", "w_max"
    };
    double cvals[13] = {
      1e-6, p_min, p_max, p_clip_max, (double) n, (double) p_capped, (double) p_capped_ctrl,
      r_min, r_max, r_floor_min, (double) m_cohorts, (double) r_floored, w_max
    };
    for (int k = 0; k < 13; ++k) {
      SET_STRING_ELT(clip_nms, k, Rf_mkChar(cnms[k]));
      SET_VECTOR_ELT(clip_s, k, Rf_ScalarReal(cvals[k]));
    }
    Rf_setAttrib(clip_s, R_NamesSymbol, clip_nms);

    SEXP piv_s = PROTECT(Rf_allocVector(INTSXP, p));
    ++nprot;
    for (int j = 0; j < p; ++j) INTEGER(piv_s)[j] = j + 1;

    SEXP finfo_s = PROTECT(Rf_allocVector(VECSXP, 12));
    SEXP finfo_nms = PROTECT(Rf_allocVector(STRSXP, 12));
    nprot += 2;
    const char *fnms[12] = {
      "qr_R", "qr_pivot", "qr_tol", "rank", "converged", "boundary",
      "iter", "df_residual", "deviance", "null_deviance", "aic", "clipping"
    };
    for (int k = 0; k < 12; ++k) SET_STRING_ELT(finfo_nms, k, Rf_mkChar(fnms[k]));
    SET_VECTOR_ELT(finfo_s, 0, qr_R_s);
    SET_VECTOR_ELT(finfo_s, 1, piv_s);
    SET_VECTOR_ELT(finfo_s, 2, Rf_ScalarReal(qr_tol));
    SET_VECTOR_ELT(finfo_s, 3, Rf_ScalarInteger(p));
    SET_VECTOR_ELT(finfo_s, 4, Rf_ScalarLogical(1));
    SET_VECTOR_ELT(finfo_s, 5, Rf_ScalarLogical(0));
    SET_VECTOR_ELT(finfo_s, 6, Rf_ScalarInteger(iter_used));
    SET_VECTOR_ELT(finfo_s, 7, Rf_ScalarReal((double) (n_pos_w - p)));
    SET_VECTOR_ELT(finfo_s, 8, Rf_ScalarReal(dev));
    SET_VECTOR_ELT(finfo_s, 9, Rf_ScalarReal(null_dev));
    SET_VECTOR_ELT(finfo_s, 10, Rf_ScalarReal(aic));
    SET_VECTOR_ELT(finfo_s, 11, clip_s);
    Rf_setAttrib(finfo_s, R_NamesSymbol, finfo_nms);

    SEXP st_s = PROTECT(Rf_allocVector(VECSXP, 24));
    SEXP st_nms = PROTECT(Rf_allocVector(STRSXP, 24));
    nprot += 2;
    const char *snms[24] = {
      "n", "m_cohorts", "p", "risk_a_idx", "risk_gi", "risk_w",
      "pair_count", "pair_mass", "rates", "Z_use", "D_use", "fit_weights",
      "beta", "odds_use", "prob_use", "is_t", "is_c", "a_obs", "b_obs",
      "sum_a", "mean_a", "sum_b", "mean_b", "fit_info"
    };
    for (int k = 0; k < 24; ++k) SET_STRING_ELT(st_nms, k, Rf_mkChar(snms[k]));
    SET_VECTOR_ELT(st_s, 0, n_s);
    SET_VECTOR_ELT(st_s, 1, m_coh_s);
    SET_VECTOR_ELT(st_s, 2, Rf_ScalarInteger(p));
    SET_VECTOR_ELT(st_s, 3, risk_a_idx_s);
    SET_VECTOR_ELT(st_s, 4, D_sexp);
    SET_VECTOR_ELT(st_s, 5, w_sexp);
    SET_VECTOR_ELT(st_s, 6, pair_count_s);
    SET_VECTOR_ELT(st_s, 7, pair_mass_s);
    SET_VECTOR_ELT(st_s, 8, rates_s);
    SET_VECTOR_ELT(st_s, 9, Z_sexp);
    SET_VECTOR_ELT(st_s, 10, D_use_s);
    SET_VECTOR_ELT(st_s, 11, fit_w_s);
    SET_VECTOR_ELT(st_s, 12, beta_s);
    SET_VECTOR_ELT(st_s, 13, odds_s);
    SET_VECTOR_ELT(st_s, 14, prob_s);
    SET_VECTOR_ELT(st_s, 15, is_t_s);
    SET_VECTOR_ELT(st_s, 16, is_c_s);
    SET_VECTOR_ELT(st_s, 17, a_obs_s);
    SET_VECTOR_ELT(st_s, 18, b_obs_s);
    SET_VECTOR_ELT(st_s, 19, Rf_ScalarReal(sum_a));
    SET_VECTOR_ELT(st_s, 20, Rf_ScalarReal(mean_a));
    SET_VECTOR_ELT(st_s, 21, Rf_ScalarReal(sum_b));
    SET_VECTOR_ELT(st_s, 22, Rf_ScalarReal(mean_b));
    SET_VECTOR_ELT(st_s, 23, finfo_s);
    Rf_setAttrib(st_s, R_NamesSymbol, st_nms);

    UNPROTECT(nprot);
    vmaxset(vmax);
    return st_s;
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 7));
  nprot += 2;
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("beta"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("qr_R"));
  SET_STRING_ELT(nms_s, 2, Rf_mkChar("iter"));
  SET_STRING_ELT(nms_s, 3, Rf_mkChar("df_residual"));
  SET_STRING_ELT(nms_s, 4, Rf_mkChar("deviance"));
  SET_STRING_ELT(nms_s, 5, Rf_mkChar("null_deviance"));
  SET_STRING_ELT(nms_s, 6, Rf_mkChar("aic"));
  SET_VECTOR_ELT(out_s, 0, beta_s);
  SET_VECTOR_ELT(out_s, 1, qr_R_s);
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarInteger(iter_used));
  SET_VECTOR_ELT(out_s, 3, Rf_ScalarReal((double) (n_pos_w - p)));
  SET_VECTOR_ELT(out_s, 4, Rf_ScalarReal(dev));
  SET_VECTOR_ELT(out_s, 5, Rf_ScalarReal(null_dev));
  SET_VECTOR_ELT(out_s, 6, Rf_ScalarReal(aic));
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);

  UNPROTECT(nprot);
  vmaxset(vmax);
  return out_s;
}

/* 14. Outcome-free Hessian, Cholesky & all-cohort rate Jacobian A_beta_r (mode 1: Z_sexp is REALSXP),
       OR outcome means & replay loadings ipw_link_means + ipw_link_partials (mode 2: Z_sexp is VECSXP st). */
SEXP C_cdid_ipw_hessian(
    SEXP Z_sexp,
    SEXP fit_w_sexp,
    SEXP prob_sexp,
    SEXP risk_gi_sexp,
    SEXP risk_w_sexp,
    SEXP rates_sexp,
    SEXP n_m_sexp) {
  if (TYPEOF(Z_sexp) == VECSXP) {
    /* Mode 2: Z_sexp = st, fit_w_sexp = hs, prob_sexp = risk_dy (or y_vec),
               risk_gi_sexp = obs_s (or NULL), risk_w_sexp = obs_t (or NULL),
               rates_sexp = q, n_m_sexp = cohort_labels */
    if (TYPEOF(fit_w_sexp) != VECSXP || TYPEOF(prob_sexp) != REALSXP ||
        TYPEOF(rates_sexp) != REALSXP || TYPEOF(n_m_sexp) != STRSXP) {
      return R_NilValue;
    }
    SEXP Z_use_s = cdid_list_get_named(Z_sexp, "Z_use");
    SEXP is_t_s = cdid_list_get_named(Z_sexp, "is_t");
    SEXP is_c_s = cdid_list_get_named(Z_sexp, "is_c");
    SEXP a_obs_s = cdid_list_get_named(Z_sexp, "a_obs");
    SEXP b_obs_s = cdid_list_get_named(Z_sexp, "b_obs");
    SEXP sum_a_s = cdid_list_get_named(Z_sexp, "sum_a");
    SEXP sum_b_s = cdid_list_get_named(Z_sexp, "sum_b");
    SEXP prob_u_s = cdid_list_get_named(Z_sexp, "prob_use");
    SEXP odds_u_s = cdid_list_get_named(Z_sexp, "odds_use");
    SEXP rgi_s = cdid_list_get_named(Z_sexp, "risk_gi");
    SEXP rw_s = cdid_list_get_named(Z_sexp, "risk_w");
    SEXP st_rates_s = cdid_list_get_named(Z_sexp, "rates");
    SEXP H_chol_s = cdid_list_get_named(fit_w_sexp, "H_chol");
    SEXP A_br_s = cdid_list_get_named(fit_w_sexp, "A_beta_r");

    if (TYPEOF(Z_use_s) != REALSXP || !Rf_isMatrix(Z_use_s) ||
        TYPEOF(is_t_s) != LGLSXP || TYPEOF(is_c_s) != LGLSXP ||
        TYPEOF(a_obs_s) != REALSXP || TYPEOF(b_obs_s) != REALSXP ||
        TYPEOF(prob_u_s) != REALSXP || TYPEOF(odds_u_s) != REALSXP ||
        TYPEOF(rgi_s) != INTSXP || TYPEOF(rw_s) != REALSXP ||
        TYPEOF(st_rates_s) != REALSXP ||
        TYPEOF(H_chol_s) != REALSXP || !Rf_isMatrix(H_chol_s) ||
        TYPEOF(A_br_s) != REALSXP || !Rf_isMatrix(A_br_s)) {
      return R_NilValue;
    }

    int n_obs = Rf_nrows(Z_use_s);
    int p = Rf_ncols(Z_use_s);
    int m_cohorts = (int) XLENGTH(st_rates_s);
    int n_t = (int) XLENGTH(a_obs_s);
    int n_c = (int) XLENGTH(b_obs_s);
    double sum_a = Rf_asReal(sum_a_s);
    double sum_b = Rf_asReal(sum_b_s);

    if (n_obs < 1 || p < 1 || m_cohorts < 2 || n_t + n_c != n_obs ||
        XLENGTH(is_t_s) != (R_xlen_t) n_obs || XLENGTH(is_c_s) != (R_xlen_t) n_obs ||
        XLENGTH(prob_u_s) != (R_xlen_t) n_obs || XLENGTH(odds_u_s) != (R_xlen_t) n_obs ||
        XLENGTH(rgi_s) != (R_xlen_t) n_obs || XLENGTH(rw_s) != (R_xlen_t) n_obs ||
        XLENGTH(rates_sexp) != (R_xlen_t) m_cohorts || XLENGTH(n_m_sexp) != (R_xlen_t) m_cohorts ||
        Rf_nrows(H_chol_s) != p || Rf_ncols(H_chol_s) != p ||
        Rf_nrows(A_br_s) != p || Rf_ncols(A_br_s) != m_cohorts ||
        !R_FINITE(sum_a) || sum_a == 0.0 || !R_FINITE(sum_b) || sum_b == 0.0) {
      return R_NilValue;
    }

    const double *A_br = REAL(A_br_s);
    R_xlen_t tot_abr = (R_xlen_t) p * (R_xlen_t) m_cohorts;
    for (R_xlen_t k = 0; k < tot_abr; ++k) {
      if (!R_FINITE(A_br[k])) return R_NilValue;
    }

    const void *vmax = vmaxget();
    const double *dy = NULL;
    if (risk_gi_sexp != R_NilValue) {
      if (TYPEOF(risk_gi_sexp) != INTSXP || TYPEOF(risk_w_sexp) != INTSXP ||
          XLENGTH(risk_gi_sexp) != (R_xlen_t) n_obs || XLENGTH(risk_w_sexp) != (R_xlen_t) n_obs) {
        vmaxset(vmax);
        return R_NilValue;
      }
      const int *obs_s = INTEGER(risk_gi_sexp);
      const int *obs_t = INTEGER(risk_w_sexp);
      const double *y_vec = REAL(prob_sexp);
      R_xlen_t ny = XLENGTH(prob_sexp);
      double *dy_buf = (double *) R_alloc((size_t) n_obs, sizeof(double));
      for (int i = 0; i < n_obs; ++i) {
        int ps = obs_s[i], pt = obs_t[i];
        if (ps < 1 || (R_xlen_t) ps > ny || pt < 1 || (R_xlen_t) pt > ny) {
          vmaxset(vmax);
          return R_NilValue;
        }
        double ys = y_vec[ps - 1], yt = y_vec[pt - 1];
        if (!R_FINITE(ys) || !R_FINITE(yt)) {
          vmaxset(vmax);
          return R_NilValue;
        }
        double d = yt - ys;
        if (!R_FINITE(d)) {
          vmaxset(vmax);
          return R_NilValue;
        }
        dy_buf[i] = d;
      }
      dy = dy_buf;
    } else {
      if (XLENGTH(prob_sexp) != (R_xlen_t) n_obs) {
        vmaxset(vmax);
        return R_NilValue;
      }
      dy = REAL(prob_sexp);
      for (int i = 0; i < n_obs; ++i) {
        if (!R_FINITE(dy[i])) {
          vmaxset(vmax);
          return R_NilValue;
        }
      }
    }

    const int *is_t = LOGICAL(is_t_s);
    const int *is_c = LOGICAL(is_c_s);
    const double *a_obs = REAL(a_obs_s);
    const double *b_obs = REAL(b_obs_s);
    const double *prob_u = REAL(prob_u_s);
    const double *odds_u = REAL(odds_u_s);
    const int *rgi = INTEGER(rgi_s);
    const double *rw = REAL(rw_s);
    const double *rates = REAL(st_rates_s);
    const double *q = REAL(rates_sexp);
    const double *Z_use = REAL(Z_use_s);
    const double *H_chol = REAL(H_chol_s);

    int *c_idx = (int *) R_alloc((size_t) (n_c > 0 ? n_c : 1), sizeof(int));
    double *wt_c = (double *) R_alloc((size_t) (n_c > 0 ? n_c : 1), sizeof(double));
    int it = 0, ic = 0;
    long double s_t_ld = 0.0L, s_c_ld = 0.0L;
    for (int i = 0; i < n_obs; ++i) {
      if (is_t[i]) {
        if (it >= n_t) { vmaxset(vmax); return R_NilValue; }
        s_t_ld += (long double) (a_obs[it] * dy[i]);
        ++it;
      }
      if (is_c[i]) {
        if (ic >= n_c) { vmaxset(vmax); return R_NilValue; }
        s_c_ld += (long double) (b_obs[ic] * dy[i]);
        c_idx[ic] = i;
        ++ic;
      }
    }
    if (it != n_t || ic != n_c) { vmaxset(vmax); return R_NilValue; }
    double Tmean = ((double) s_t_ld) / sum_a;
    double Cmean = ((double) s_c_ld) / sum_b;
    double mu = Tmean - Cmean;
    if (!R_FINITE(Tmean) || !R_FINITE(Cmean) || !R_FINITE(mu)) {
      vmaxset(vmax);
      return R_NilValue;
    }

    const double clip_cap = 1.0 - 1e-6;
    for (int k = 0; k < n_c; ++k) {
      int i = c_idx[k];
      double act = (prob_u[i] < clip_cap) ? 1.0 : 0.0;
      wt_c[k] = act * b_obs[k] * (dy[i] - Cmean);
    }

    SEXP G_beta_s = PROTECT(Rf_allocVector(REALSXP, p));
    SEXP G_r_s = PROTECT(Rf_allocVector(REALSXP, m_cohorts));
    SEXP l_beta_s = PROTECT(Rf_allocVector(REALSXP, p));
    SEXP l_rate_s = PROTECT(Rf_allocVector(REALSXP, m_cohorts));
    SEXP rc_s0_s = PROTECT(Rf_allocVector(REALSXP, m_cohorts));
    SEXP rc_s1_s = PROTECT(Rf_allocVector(REALSXP, m_cohorts));
    double *G_beta = REAL(G_beta_s);
    double *G_r = REAL(G_r_s);
    double *l_beta = REAL(l_beta_s);
    double *l_rate = REAL(l_rate_s);
    double *rc_s0 = REAL(rc_s0_s);
    double *rc_s1 = REAL(rc_s1_s);

    for (int col = 0; col < p; ++col) {
      const double *zcol = Z_use + (size_t) col * (size_t) n_obs;
      long double s_ld = 0.0L;
      for (int k = 0; k < n_c; ++k) s_ld += (long double) (zcol[c_idx[k]] * wt_c[k]);
      double gb = -((double) s_ld) / sum_b;
      if (!R_FINITE(gb)) {
        UNPROTECT(6);
        vmaxset(vmax);
        return R_NilValue;
      }
      G_beta[col] = gb;
    }

    G_r[0] = 0.0;
    for (int j = 1; j < m_cohorts; ++j) {
      if (rates[j] <= 1e-6) {
        G_r[j] = 0.0;
      } else {
        double rc = rates[j];
        int gj = j + 1;
        long double s_ld = 0.0L;
        for (int i = 0; i < n_obs; ++i) {
          if (rgi[i] == gj) {
            double term = (rw[i] * odds_u[i]) * (dy[i] - Cmean);
            s_ld += (long double) term;
          }
        }
        double gr = ((double) s_ld) / (rc * rc * sum_b);
        if (!R_FINITE(gr)) {
          UNPROTECT(6);
          vmaxset(vmax);
          return R_NilValue;
        }
        G_r[j] = gr;
      }
    }

    /* Forward solve t(H_chol) * y_tmp = G_beta, then backward solve H_chol * l_beta = y_tmp via dtrsm */
    for (int c = 0; c < p; ++c) {
      double r_cc = H_chol[c + (size_t) c * (size_t) p];
      if (!R_FINITE(r_cc) || r_cc == 0.0) {
        UNPROTECT(6);
        vmaxset(vmax);
        return R_NilValue;
      }
      l_beta[c] = G_beta[c];
    }
    double *L_low = (double *) R_alloc((size_t) p * (size_t) p, sizeof(double));
    for (int c = 0; c < p; ++c) {
      for (int r = 0; r < p; ++r) {
        L_low[r + (size_t) c * (size_t) p] = (r >= c) ? H_chol[c + (size_t) r * (size_t) p] : 0.0;
      }
    }
    int p_b = p, one_b = 1, m_b = m_cohorts;
    double one_d = 1.0, zero_d = 0.0;
    F77_CALL(dtrsm)("L", "L", "N", "N", &p_b, &one_b, &one_d, L_low, &p_b, l_beta, &p_b FCONE FCONE FCONE FCONE);
    F77_CALL(dtrsm)("L", "U", "N", "N", &p_b, &one_b, &one_d, H_chol, &p_b, l_beta, &p_b FCONE FCONE FCONE FCONE);
    for (int c = 0; c < p; ++c) {
      if (!R_FINITE(l_beta[c])) {
        UNPROTECT(6);
        vmaxset(vmax);
        return R_NilValue;
      }
    }

    double *abr_dot = (double *) R_alloc((size_t) m_cohorts, sizeof(double));
    F77_CALL(dgemm)("T", "N", &m_b, &one_b, &p_b, &one_d,
                    A_br, &p_b, l_beta, &p_b, &zero_d, abr_dot, &m_b FCONE FCONE);

    for (int j = 0; j < m_cohorts; ++j) {
      double lr = G_r[j] + abr_dot[j];
      double rj = rates[j], qj = q[j];
      double r0 = -(rj / qj) * lr;
      double r1 = ((1.0 - rj) / qj) * lr;
      if (!R_FINITE(lr) || !R_FINITE(r0) || !R_FINITE(r1)) {
        UNPROTECT(6);
        vmaxset(vmax);
        return R_NilValue;
      }
      l_rate[j] = lr;
      rc_s0[j] = r0;
      rc_s1[j] = r1;
    }

    SEXP zdn = Rf_getAttrib(Z_use_s, R_DimNamesSymbol);
    if (zdn != R_NilValue && TYPEOF(zdn) == VECSXP && XLENGTH(zdn) == 2) {
      SEXP cn = VECTOR_ELT(zdn, 1);
      if (cn != R_NilValue) {
        Rf_setAttrib(G_beta_s, R_NamesSymbol, cn);
        Rf_setAttrib(l_beta_s, R_NamesSymbol, cn);
      }
    }
    Rf_setAttrib(G_r_s, R_NamesSymbol, n_m_sexp);
    Rf_setAttrib(l_rate_s, R_NamesSymbol, n_m_sexp);
    Rf_setAttrib(rc_s0_s, R_NamesSymbol, n_m_sexp);
    Rf_setAttrib(rc_s1_s, R_NamesSymbol, n_m_sexp);

    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 9));
    SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 9));
    const char *onms[9] = {
      "Tmean", "Cmean", "mu", "G_beta", "G_r",
      "loading_beta", "loading_rate", "rate_coef_s0", "rate_coef_s1"
    };
    for (int k = 0; k < 9; ++k) SET_STRING_ELT(nms_s, k, Rf_mkChar(onms[k]));
    SET_VECTOR_ELT(out_s, 0, Rf_ScalarReal(Tmean));
    SET_VECTOR_ELT(out_s, 1, Rf_ScalarReal(Cmean));
    SET_VECTOR_ELT(out_s, 2, Rf_ScalarReal(mu));
    SET_VECTOR_ELT(out_s, 3, G_beta_s);
    SET_VECTOR_ELT(out_s, 4, G_r_s);
    SET_VECTOR_ELT(out_s, 5, l_beta_s);
    SET_VECTOR_ELT(out_s, 6, l_rate_s);
    SET_VECTOR_ELT(out_s, 7, rc_s0_s);
    SET_VECTOR_ELT(out_s, 8, rc_s1_s);
    Rf_setAttrib(out_s, R_NamesSymbol, nms_s);

    UNPROTECT(8);
    vmaxset(vmax);
    return out_s;
  }

  if (TYPEOF(Z_sexp) != REALSXP || Rf_isObject(Z_sexp) ||
      TYPEOF(fit_w_sexp) != REALSXP || TYPEOF(prob_sexp) != REALSXP ||
      TYPEOF(risk_gi_sexp) != INTSXP || TYPEOF(risk_w_sexp) != REALSXP ||
      TYPEOF(rates_sexp) != REALSXP || TYPEOF(n_m_sexp) != REALSXP || XLENGTH(n_m_sexp) != 2) {
    return R_NilValue;
  }
  SEXP zdim = Rf_getAttrib(Z_sexp, R_DimSymbol);
  if (TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return R_NilValue;
  int n_obs = INTEGER(zdim)[0];
  int p = INTEGER(zdim)[1];
  double n_val = REAL(n_m_sexp)[0];
  int m_cohorts = (int) REAL(n_m_sexp)[1];
  if (n_obs < 1 || p < 1 || p > 256 || m_cohorts < 1 || !R_FINITE(n_val) || n_val <= 0.0 ||
      XLENGTH(fit_w_sexp) != (R_xlen_t) n_obs || XLENGTH(prob_sexp) != (R_xlen_t) n_obs ||
      XLENGTH(risk_gi_sexp) != (R_xlen_t) n_obs || XLENGTH(risk_w_sexp) != (R_xlen_t) n_obs ||
      XLENGTH(rates_sexp) != (R_xlen_t) m_cohorts) {
    return R_NilValue;
  }

  const double *Z = REAL(Z_sexp);
  const double *fit_w = REAL(fit_w_sexp);
  const double *prob = REAL(prob_sexp);
  const int *risk_gi = INTEGER(risk_gi_sexp);
  const double *risk_w = REAL(risk_w_sexp);
  const double *rates = REAL(rates_sexp);

  const void *vmax = vmaxget();
  double *w_hess = (double *) R_alloc((size_t) n_obs, sizeof(double));
  double *W_Z = (double *) R_alloc((size_t) n_obs * (size_t) p, sizeof(double));
  for (int i = 0; i < n_obs; ++i) {
    w_hess[i] = fit_w[i] * prob[i] * (1.0 - prob[i]);
  }
  for (int j = 0; j < p; ++j) {
    const double *zj = Z + (size_t) j * (size_t) n_obs;
    double *wzj = W_Z + (size_t) j * (size_t) n_obs;
    for (int i = 0; i < n_obs; ++i) wzj[i] = w_hess[i] * zj[i];
  }

  SEXP H_s = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP R_s = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  SEXP A_s = PROTECT(Rf_allocMatrix(REALSXP, p, m_cohorts));
  double *H = REAL(H_s);
  double *R_chol = REAL(R_s);
  double *A_br = REAL(A_s);
  memset(A_br, 0, (size_t) p * (size_t) m_cohorts * sizeof(double));

  int n_b = n_obs, p_b = p;
  double one_d = 1.0, zero_d = 0.0;
  F77_CALL(dgemm)("T", "N", &p_b, &p_b, &n_b, &one_d,
                  Z, &n_b, W_Z, &n_b, &zero_d, H, &p_b FCONE FCONE);
  for (int idx = 0; idx < p * p; ++idx) {
    H[idx] /= n_val;
    if (!R_FINITE(H[idx])) {
      vmaxset(vmax);
      UNPROTECT(3);
      return R_NilValue;
    }
  }

  if (!cdid_chol_upper(H, R_chol, p, 0.0)) {
    vmaxset(vmax);
    UNPROTECT(3);
    return R_NilValue;
  }

  const double clip_eps = 1e-6;
  double *rc2 = (double *) R_alloc((size_t) m_cohorts, sizeof(double));
  long double *A_ld = R_allocLD((size_t) p * (size_t) m_cohorts);
  for (int j = 0; j < m_cohorts; ++j) {
    if (rates[j] > clip_eps) {
      double rc = (rates[j] > clip_eps) ? rates[j] : clip_eps;
      rc2[j] = rc * rc;
    } else {
      rc2[j] = 0.0;
    }
    for (int c = 0; c < p; ++c) A_ld[c + (size_t) j * (size_t) p] = 0.0L;
  }

  for (int i = 0; i < n_obs; ++i) {
    int g = risk_gi[i];
    if (g >= 1 && g <= m_cohorts && rc2[g - 1] > 0.0) {
      double d_i = (g == 1) ? 1.0 : 0.0;
      double wt_i = (risk_w[i] * (d_i - prob[i])) / rc2[g - 1];
      long double *col_a = A_ld + (size_t) (g - 1) * (size_t) p;
      for (int c = 0; c < p; ++c) {
        double prod = Z[i + (size_t) c * (size_t) n_obs] * wt_i;
        col_a[c] += (long double) prod;
      }
    }
  }
  for (int j = 0; j < m_cohorts; ++j) {
    double *col_a = A_br + (size_t) j * (size_t) p;
    const long double *col_ld = A_ld + (size_t) j * (size_t) p;
    for (int c = 0; c < p; ++c) {
      col_a[c] = -((double) col_ld[c]) / n_val;
      if (!R_FINITE(col_a[c])) {
        vmaxset(vmax);
        UNPROTECT(3);
        return R_NilValue;
      }
    }
  }

  SEXP zdn = Rf_getAttrib(Z_sexp, R_DimNamesSymbol);
  if (zdn != R_NilValue && TYPEOF(zdn) == VECSXP && XLENGTH(zdn) == 2) {
    SEXP cn = VECTOR_ELT(zdn, 1);
    if (cn != R_NilValue) {
      SEXP hdn = PROTECT(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(hdn, 0, cn);
      SET_VECTOR_ELT(hdn, 1, cn);
      Rf_setAttrib(H_s, R_DimNamesSymbol, hdn);
      Rf_setAttrib(R_s, R_DimNamesSymbol, hdn);
      UNPROTECT(1);
    }
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("H"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("H_chol"));
  SET_STRING_ELT(nms_s, 2, Rf_mkChar("A_beta_r"));
  SET_VECTOR_ELT(out_s, 0, H_s);
  SET_VECTOR_ELT(out_s, 1, R_s);
  SET_VECTOR_ELT(out_s, 2, A_s);
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);

  UNPROTECT(5);
  vmaxset(vmax);
  return out_s;
}

/* 15. Raw-IPW chunk IF evaluator (!want_details) */
SEXP C_cdid_ipw_raw_if_chunk(
    SEXP rows_sexp,
    SEXP post_a_sexp,
    SEXP post_obs_sexp,
    SEXP y_sexp,
    SEXP a_cohort_sexp,
    SEXP w_a_sexp,
    SEXP Z_a_sexp,
    SEXP pair_tab_sexp,
    SEXP pair_active_coh_sexp,
    SEXP trt_int_sexp,
    SEXP trt_dbl_list_sexp) {
  if (TYPEOF(rows_sexp) != INTSXP || TYPEOF(post_a_sexp) != INTSXP || Rf_isObject(post_a_sexp) ||
      TYPEOF(post_obs_sexp) != INTSXP || Rf_isObject(post_obs_sexp) ||
      TYPEOF(y_sexp) != REALSXP || Rf_isObject(y_sexp) ||
      TYPEOF(a_cohort_sexp) != INTSXP || Rf_isObject(a_cohort_sexp) ||
      TYPEOF(w_a_sexp) != REALSXP || Rf_isObject(w_a_sexp) ||
      TYPEOF(Z_a_sexp) != REALSXP || Rf_isObject(Z_a_sexp) ||
      TYPEOF(pair_tab_sexp) != INTSXP || TYPEOF(pair_active_coh_sexp) != INTSXP ||
      TYPEOF(trt_int_sexp) != INTSXP || TYPEOF(trt_dbl_list_sexp) != VECSXP) {
    return R_NilValue;
  }

  int B = (int) XLENGTH(rows_sexp);
  int n_a = (int) XLENGTH(a_cohort_sexp);
  R_xlen_t n_post = XLENGTH(post_a_sexp);
  R_xlen_t n_y = XLENGTH(y_sexp);
  SEXP zdim = Rf_getAttrib(Z_a_sexp, R_DimSymbol);
  SEXP pdim = Rf_getAttrib(pair_tab_sexp, R_DimSymbol);
  SEXP adim = Rf_getAttrib(pair_active_coh_sexp, R_DimSymbol);
  SEXP tdim = Rf_getAttrib(trt_int_sexp, R_DimSymbol);
  if (B < 1 || n_a < 1 || XLENGTH(post_obs_sexp) != n_post || XLENGTH(w_a_sexp) != (R_xlen_t) n_a ||
      TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2 || INTEGER(zdim)[0] != n_a ||
      TYPEOF(pdim) != INTSXP || XLENGTH(pdim) != 2 || INTEGER(pdim)[0] != 5 ||
      TYPEOF(adim) != INTSXP || XLENGTH(adim) != 2 ||
      TYPEOF(tdim) != INTSXP || XLENGTH(tdim) != 2) {
    return R_NilValue;
  }

  int p = INTEGER(zdim)[1];
  int n_pids = INTEGER(pdim)[1];
  int c_a = INTEGER(adim)[0];
  int n_tids = INTEGER(tdim)[1];
  if (p < 1 || n_pids < 1 || c_a < 1 || n_tids < 1 ||
      INTEGER(adim)[1] != n_pids || INTEGER(tdim)[0] != (4 + c_a) ||
      XLENGTH(trt_dbl_list_sexp) != (R_xlen_t) (n_tids + 1)) {
    return R_NilValue;
  }

  SEXP mpos_s = VECTOR_ELT(trt_dbl_list_sexp, n_tids);
  if (TYPEOF(mpos_s) != INTSXP || XLENGTH(mpos_s) < 1) return R_NilValue;
  int L_req = (int) XLENGTH(mpos_s);
  const int *model_pos = INTEGER(mpos_s);
  for (int j = 0; j < L_req; ++j) {
    if (model_pos[j] < 1 || model_pos[j] > n_tids) return R_NilValue;
  }

  const int *rows = INTEGER(rows_sexp);
  const int *post_a = INTEGER(post_a_sexp);
  const int *post_obs = INTEGER(post_obs_sexp);
  const double *y = REAL(y_sexp);
  const int *a_coh = INTEGER(a_cohort_sexp);
  const double *w_a = REAL(w_a_sexp);
  const double *Z_a = REAL(Z_a_sexp);
  const int *pair_tab = INTEGER(pair_tab_sexp);
  const int *pair_act = INTEGER(pair_active_coh_sexp);
  const int *trt_int = INTEGER(trt_int_sexp);

  int a_s = rows[0], a_e = rows[B - 1];
  int is_contig = ((a_e - a_s + 1) == B);
  for (int i = 0; i < B; ++i) {
    int r = rows[i];
    if (r < 1 || r > n_a) return R_NilValue;
    if (is_contig && r != a_s + i) is_contig = 0;
  }

  const void *vmax = vmaxget();
  int *u_coh = (int *) R_alloc((size_t) B, sizeof(int));
  double *u_w = (double *) R_alloc((size_t) B, sizeof(double));
  for (int i = 0; i < B; ++i) {
    u_coh[i] = a_coh[rows[i] - 1];
    u_w[i] = w_a[rows[i] - 1];
  }

  uint8_t *obs_mat = (uint8_t *) R_alloc((size_t) B * (size_t) n_pids, sizeof(uint8_t));
  double *dy_mat = (double *) R_alloc((size_t) B * (size_t) n_pids, sizeof(double));
  memset(obs_mat, 0, (size_t) B * (size_t) n_pids * sizeof(uint8_t));
  memset(dy_mat, 0, (size_t) B * (size_t) n_pids * sizeof(double));
  int *pos_s_buf = (int *) R_alloc((size_t) B, sizeof(int));
  int *pos_t_buf = (int *) R_alloc((size_t) B, sizeof(int));

  /* 1. Query each active pair once on `rows` */
  for (int pi = 0; pi < n_pids; ++pi) {
    const int *ptab = pair_tab + (size_t) pi * 5;
    int s_start = ptab[0], s_end = ptab[1], t_start = ptab[2], t_end = ptab[3], comp_flag = ptab[4];
    if (s_start <= 0 || t_start <= 0) continue;
    R_xlen_t s_lo = (R_xlen_t) s_start - 1, s_hi = (R_xlen_t) s_end - 1;
    R_xlen_t t_lo = (R_xlen_t) t_start - 1, t_hi = (R_xlen_t) t_end - 1;
    if (s_hi < s_lo || s_hi >= n_post || t_hi < t_lo || t_hi >= n_post) {
      vmaxset(vmax);
      return R_NilValue;
    }

    uint8_t *obs_col = obs_mat + (size_t) pi * (size_t) B;
    double *dy_col = dy_mat + (size_t) pi * (size_t) B;
    const int *act_col = pair_act + (size_t) pi * (size_t) c_a;

    if (is_contig && comp_flag) {
      if ((s_hi - s_lo + 1) < (R_xlen_t) a_e || (t_hi - t_lo + 1) < (R_xlen_t) a_e) {
        vmaxset(vmax);
        return R_NilValue;
      }
      for (int i = 0; i < B; ++i) {
        int a = a_s + i;
        R_xlen_t ps = s_lo + (a - 1), pt = t_lo + (a - 1);
        if (post_a[ps] != a || post_a[pt] != a) {
          vmaxset(vmax);
          return R_NilValue;
        }
        obs_col[i] = 1;
        pos_s_buf[i] = post_obs[ps];
        pos_t_buf[i] = post_obs[pt];
      }
    } else if (is_contig) {
      int hs = (int) s_lo, ht = (int) t_lo;
      /* Binary search first >= a_s */
      int lo = (int) s_lo, hi = (int) s_hi;
      hs = hi + 1;
      while (lo <= hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (post_a[mid] >= a_s) { hs = mid; hi = mid - 1; } else lo = mid + 1;
      }
      lo = (int) t_lo; hi = (int) t_hi;
      ht = hi + 1;
      while (lo <= hi) {
        int mid = lo + ((hi - lo) >> 1);
        if (post_a[mid] >= a_s) { ht = mid; hi = mid - 1; } else lo = mid + 1;
      }
      while (hs <= (int) s_hi && ht <= (int) t_hi) {
        int vs = post_a[hs], vt = post_a[ht];
        if (vs > a_e || vt > a_e) break;
        if (vs == vt) {
          int loc = vs - a_s;
          obs_col[loc] = 1;
          pos_s_buf[loc] = post_obs[hs];
          pos_t_buf[loc] = post_obs[ht];
          ++hs;
          ++ht;
        } else if (vs < vt) ++hs;
        else ++ht;
      }
    } else {
      int hs_hint = (int) s_lo, ht_hint = (int) t_lo;
      for (int i = 0; i < B; ++i) {
        int a = rows[i];
        int hs = cdid_find_posting(post_a, (int) s_lo, (int) s_hi, a, &hs_hint);
        if (hs >= 0) {
          int ht = cdid_find_posting(post_a, (int) t_lo, (int) t_hi, a, &ht_hint);
          if (ht >= 0) {
            obs_col[i] = 1;
            pos_s_buf[i] = post_obs[hs];
            pos_t_buf[i] = post_obs[ht];
          }
        }
      }
    }

    for (int i = 0; i < B; ++i) {
      if (obs_col[i]) {
        int c = u_coh[i];
        if (c >= 1 && c <= c_a && act_col[c - 1]) {
          int os = pos_s_buf[i], ot = pos_t_buf[i];
          if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y) {
            vmaxset(vmax);
            return R_NilValue;
          }
          double yl = y[os - 1], yr = y[ot - 1];
          double dy = yr - yl;
          if (!R_FINITE(yl) || !R_FINITE(yr) || !R_FINITE(dy)) {
            vmaxset(vmax);
            return R_NilValue;
          }
          dy_col[i] = dy;
        }
      }
    }
  }

  /* 3. Evaluate distinct treatment models reusing shared-nuisance buffer */
  double *model_phi = (double *) R_alloc((size_t) B * (size_t) n_tids, sizeof(double));
  int *in_risk = (int *) R_alloc((size_t) B, sizeof(int));
  int *gi_chunk = (int *) R_alloc((size_t) B, sizeof(int));
  int *pos_buf = (int *) R_alloc((size_t) B, sizeof(int));
  int *gj_buf = (int *) R_alloc((size_t) B, sizeof(int));
  double *Z_buf = (double *) R_alloc((size_t) B * (size_t) p, sizeof(double));
  double *w_buf = (double *) R_alloc((size_t) B, sizeof(double));
  double *r_buf = (double *) R_alloc((size_t) B, sizeof(double));
  double *odds_buf = (double *) R_alloc((size_t) B, sizeof(double));
  double *score_base = (double *) R_alloc((size_t) B, sizeof(double));
  double *eta_buf = (double *) R_alloc((size_t) B, sizeof(double));
  double *zload_buf = (double *) R_alloc((size_t) B, sizeof(double));

  double *col_rate = (double *) R_alloc((size_t) B, sizeof(double));
  double *col_score = (double *) R_alloc((size_t) B, sizeof(double));
  double *col_treated = (double *) R_alloc((size_t) B, sizeof(double));
  double *col_control = (double *) R_alloc((size_t) B, sizeof(double));

  int current_nuisance = NA_INTEGER;
  int n_in_risk = 0, n_pos = 0;
  int nuisance_evaluations = 0;
  const double clip_cap = 1.0 - 1e-6;
  const double clip_cap_odds = clip_cap / (1.0 - clip_cap);

  for (int ti = 0; ti < n_tids; ++ti) {
    const int *tcol = trt_int + (size_t) ti * (size_t) (4 + c_a);
    int pi = tcol[0] - 1;
    int nid = tcol[1];
    int m_coh = tcol[2];
    const int *cmap = tcol + 4;
    if (pi < 0 || pi >= n_pids || m_coh < 1) {
      vmaxset(vmax);
      return R_NilValue;
    }

    SEXP td_s = VECTOR_ELT(trt_dbl_list_sexp, ti);
    if (TYPEOF(td_s) != REALSXP || XLENGTH(td_s) != (R_xlen_t) (4 + 2 * p + 3 * m_coh)) {
      vmaxset(vmax);
      return R_NilValue;
    }
    const double *td = REAL(td_s);
    double mean_a = td[0], trt_mean = td[1], mean_b = td[2], ctrl_mean = td[3];
    const double *beta = td + 4;
    const double *load_beta = beta + p;
    const double *rates_clip = load_beta + p;
    const double *rc_s0 = rates_clip + m_coh;
    const double *rc_s1 = rc_s0 + m_coh;

    const uint8_t *obs_col = obs_mat + (size_t) pi * (size_t) B;
    const double *dy_col = dy_mat + (size_t) pi * (size_t) B;

    if (nid != current_nuisance || nid == NA_INTEGER) {
      n_in_risk = 0;
      n_pos = 0;
      for (int i = 0; i < B; ++i) {
        int c = u_coh[i];
        int g = (c >= 1 && c <= c_a) ? cmap[c - 1] : 0;
        gi_chunk[i] = g;
        if (g >= 1 && g <= m_coh) {
          in_risk[n_in_risk++] = i;
          if (obs_col[i]) {
            pos_buf[n_pos] = i;
            gj_buf[n_pos] = g;
            w_buf[n_pos] = u_w[i];
            r_buf[n_pos] = rates_clip[g - 1];
            int row0 = rows[i] - 1;
            for (int col = 0; col < p; ++col) {
              Z_buf[n_pos + (size_t) col * (size_t) B] = Z_a[row0 + (size_t) col * (size_t) n_a];
            }
            ++n_pos;
          }
        }
      }
      for (int col = 0; col < p; ++col) {
        if (!R_FINITE(beta[col])) {
          vmaxset(vmax);
          return R_NilValue;
        }
        const double *zc = Z_buf + (size_t) col * (size_t) B;
        for (int k = 0; k < n_pos; ++k) {
          if (!R_FINITE(zc[k])) {
            vmaxset(vmax);
            return R_NilValue;
          }
        }
      }
      if (n_pos > 0) {
        BLAS_INT npos_b = (BLAS_INT) n_pos, B_lda = (BLAS_INT) B, p_b = (BLAS_INT) p, one_b = 1;
        double one_d = 1.0, zero_d = 0.0;
        F77_CALL(dgemm)("N", "N", &npos_b, &one_b, &p_b, &one_d,
                        Z_buf, &B_lda, beta, &p_b, &zero_d, eta_buf, &npos_b FCONE FCONE);
      }
      for (int k = 0; k < n_pos; ++k) {
        double e = eta_buf[k];
        if (!R_FINITE(e)) {
          vmaxset(vmax);
          return R_NilValue;
        }
        double od, pr;
        if (e > 40.0) {
          pr = 1.0 / (1.0 + exp(-e));
          od = 1.0;
        } else {
          od = exp(e);
          pr = od / (1.0 + od);
        }
        odds_buf[k] = (pr < clip_cap) ? od : clip_cap_odds;
        double d_ind = (gj_buf[k] == 1) ? 1.0 : 0.0;
        score_base[k] = (w_buf[k] / r_buf[k]) * (d_ind - pr);
      }
      current_nuisance = nid;
      ++nuisance_evaluations;
    } else {
      /* Verify shared-nuisance pair mask identity on this chunk */
      int chk_pos = 0;
      for (int r_i = 0; r_i < n_in_risk; ++r_i) {
        int idx = in_risk[r_i];
        if (obs_col[idx]) {
          if (chk_pos >= n_pos || pos_buf[chk_pos] != idx) {
            vmaxset(vmax);
            return R_NilValue;
          }
          ++chk_pos;
        }
      }
      if (chk_pos != n_pos) {
        vmaxset(vmax);
        return R_NilValue;
      }
    }

    memset(col_rate, 0, (size_t) B * sizeof(double));
    memset(col_score, 0, (size_t) B * sizeof(double));
    memset(col_treated, 0, (size_t) B * sizeof(double));
    memset(col_control, 0, (size_t) B * sizeof(double));

    for (int r_i = 0; r_i < n_in_risk; ++r_i) {
      int idx = in_risk[r_i];
      int g = gi_chunk[idx];
      col_rate[idx] = u_w[idx] * rc_s0[g - 1];
    }

    for (int col = 0; col < p; ++col) {
      if (!R_FINITE(load_beta[col])) {
        vmaxset(vmax);
        return R_NilValue;
      }
    }
    if (n_pos > 0) {
      BLAS_INT npos_b = (BLAS_INT) n_pos, B_lda = (BLAS_INT) B, p_b = (BLAS_INT) p, one_b = 1;
      double one_d = 1.0, zero_d = 0.0;
      F77_CALL(dgemm)("N", "N", &npos_b, &one_b, &p_b, &one_d,
                      Z_buf, &B_lda, load_beta, &p_b, &zero_d, zload_buf, &npos_b FCONE FCONE);
    }

    for (int k = 0; k < n_pos; ++k) {
      int idx = pos_buf[k];
      int g = gj_buf[k];
      double dy_k = dy_col[idx];
      double wk = w_buf[k], rk = r_buf[k];
      col_rate[idx] = wk * rc_s1[g - 1];
      col_score[idx] = score_base[k] * zload_buf[k];

      if (g == 1) {
        col_treated[idx] = (wk / rk) / mean_a * (dy_k - trt_mean);
      } else {
        col_control[idx] = -(wk * odds_buf[k] / rk) / mean_b * (dy_k - ctrl_mean);
      }
    }

    double *mphi = model_phi + (size_t) ti * (size_t) B;
    for (int i = 0; i < B; ++i) {
      double v = col_treated[i] + col_control[i] + col_score[i] + col_rate[i];
      if (!R_FINITE(v)) {
        vmaxset(vmax);
        return R_NilValue;
      }
      mphi[i] = v;
    }
  }

  SEXP phi_s = PROTECT(Rf_allocMatrix(REALSXP, B, L_req));
  double *phi_out = REAL(phi_s);
  for (int j = 0; j < L_req; ++j) {
    int m_idx = model_pos[j] - 1;
    memcpy(phi_out + (size_t) j * (size_t) B,
           model_phi + (size_t) m_idx * (size_t) B,
           (size_t) B * sizeof(double));
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms_s = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms_s, 0, Rf_mkChar("phi"));
  SET_STRING_ELT(nms_s, 1, Rf_mkChar("nuisance_evaluations"));
  SET_VECTOR_ELT(out_s, 0, phi_s);
  SET_VECTOR_ELT(out_s, 1, Rf_ScalarInteger(nuisance_evaluations));
  Rf_setAttrib(out_s, R_NamesSymbol, nms_s);

  UNPROTECT(3);
  vmaxset(vmax);
  return out_s;
}

/* 16. Batched incomplete-panel OR/DR support scanner + Householder QR compressor */
SEXP C_cdid_incomplete_or_supports(
    SEXP y_sexp,
    SEXP post_a_sexp,
    SEXP post_obs_sexp,
    SEXP date_bounds_sexp,
    SEXP cohort_order_sexp,
    SEXP a_index_sexp,
    SEXP w_sexp,
    SEXP Z_a_sexp,
    SEXP spec_int_sexp,
    SEXP spec_L_sexp,
    SEXP recipe_params_sexp) {
  if (TYPEOF(y_sexp) != REALSXP || Rf_isObject(y_sexp) ||
      TYPEOF(post_a_sexp) != INTSXP || Rf_isObject(post_a_sexp) ||
      TYPEOF(post_obs_sexp) != INTSXP || Rf_isObject(post_obs_sexp) ||
      TYPEOF(date_bounds_sexp) != INTSXP ||
      TYPEOF(cohort_order_sexp) != INTSXP || Rf_isObject(cohort_order_sexp) ||
      TYPEOF(a_index_sexp) != INTSXP || Rf_isObject(a_index_sexp) ||
      TYPEOF(w_sexp) != REALSXP || Rf_isObject(w_sexp) ||
      TYPEOF(Z_a_sexp) != REALSXP || Rf_isObject(Z_a_sexp) ||
      TYPEOF(spec_int_sexp) != INTSXP || TYPEOF(spec_L_sexp) != REALSXP ||
      TYPEOF(recipe_params_sexp) != VECSXP || XLENGTH(recipe_params_sexp) != 3) {
    return R_NilValue;
  }

  SEXP db_dim = Rf_getAttrib(date_bounds_sexp, R_DimSymbol);
  SEXP z_dim = Rf_getAttrib(Z_a_sexp, R_DimSymbol);
  SEXP sp_dim = Rf_getAttrib(spec_int_sexp, R_DimSymbol);
  if (TYPEOF(db_dim) != INTSXP || XLENGTH(db_dim) != 2 || INTEGER(db_dim)[0] != 3 ||
      TYPEOF(z_dim) != INTSXP || XLENGTH(z_dim) != 2 ||
      TYPEOF(sp_dim) != INTSXP || XLENGTH(sp_dim) != 2 || INTEGER(sp_dim)[0] != 4) {
    return R_NilValue;
  }

  int n_dates = INTEGER(db_dim)[1];
  int n_a = INTEGER(z_dim)[0];
  int p_Z = INTEGER(z_dim)[1];
  int S = INTEGER(sp_dim)[1];
  int N = (int) XLENGTH(cohort_order_sexp);
  R_xlen_t n_post = XLENGTH(post_a_sexp);
  R_xlen_t n_y = XLENGTH(y_sexp);

  if (n_dates < 1 || n_a < 1 || p_Z < 1 || S < 1 || N < 1 ||
      XLENGTH(post_obs_sexp) != n_post || XLENGTH(a_index_sexp) != (R_xlen_t) N ||
      XLENGTH(w_sexp) != (R_xlen_t) N || XLENGTH(spec_L_sexp) != (R_xlen_t) S) {
    return R_NilValue;
  }

  SEXP rp_int_s = VECTOR_ELT(recipe_params_sexp, 0);
  SEXP slope_cols_s = VECTOR_ELT(recipe_params_sexp, 1);
  SEXP N_num_s = VECTOR_ELT(recipe_params_sexp, 2);
  if (TYPEOF(rp_int_s) != INTSXP || XLENGTH(rp_int_s) != 5 ||
      TYPEOF(slope_cols_s) != INTSXP || TYPEOF(N_num_s) != REALSXP || XLENGTH(N_num_s) != 1) {
    return R_NilValue;
  }

  int icpt_flag = INTEGER(rp_int_s)[0] ? 1 : 0;
  int icpt_col = INTEGER(rp_int_s)[1];
  int chunk_rows = INTEGER(rp_int_s)[2];
  int qr_buf_rows = INTEGER(rp_int_s)[3];
  int use_qr_buf = INTEGER(rp_int_s)[4] ? 1 : 0;
  int p = (int) XLENGTH(slope_cols_s);
  double N_val = REAL(N_num_s)[0];

  if ((!icpt_flag && p < 1) || p > 256 || chunk_rows < 1 || qr_buf_rows < 1 || !R_FINITE(N_val) || N_val <= 0.0) {
    return R_NilValue;
  }
  if (icpt_flag && (icpt_col < 1 || icpt_col > p_Z)) return R_NilValue;
  const int *slope_cols = INTEGER(slope_cols_s);
  for (int j = 0; j < p; ++j) {
    if (slope_cols[j] < 1 || slope_cols[j] > p_Z) return R_NilValue;
  }

  const double *y = REAL(y_sexp);
  const int *post_a = INTEGER(post_a_sexp);
  const int *post_obs = INTEGER(post_obs_sexp);
  const int *db = INTEGER(date_bounds_sexp);
  const int *coh_ord = INTEGER(cohort_order_sexp);
  const int *a_idx = INTEGER(a_index_sexp);
  const double *w = REAL(w_sexp);
  const double *Z_a = REAL(Z_a_sexp);
  const int *spec_int = INTEGER(spec_int_sexp);
  const double *spec_L = REAL(spec_L_sexp);

  int k_aug = p + 1 + icpt_flag;
  int max_buf = use_qr_buf ? qr_buf_rows : chunk_rows;
  if (max_buf > N) max_buf = N;

  const void *vmax = vmaxget();
  double *buf_D = (p > 0) ? (double *) R_alloc((size_t) max_buf * (size_t) p, sizeof(double)) : NULL;
  double *buf_dy = (double *) R_alloc((size_t) max_buf, sizeof(double));
  double *buf_w = (double *) R_alloc((size_t) max_buf, sizeof(double));
  double *aug_work = (double *) R_alloc((size_t) (max_buf + k_aug) * (size_t) k_aug, sizeof(double));
  double *fac_cur = (double *) R_alloc((size_t) k_aug * (size_t) k_aug, sizeof(double));
  double *fac_chunk = (double *) R_alloc((size_t) k_aug * (size_t) k_aug, sizeof(double));
  double *fac_merged = (double *) R_alloc((size_t) k_aug * (size_t) k_aug, sizeof(double));
  double *stack_work = (double *) R_alloc((size_t) (2 * k_aug) * (size_t) k_aug, sizeof(double));
  int *jpvt = (int *) R_alloc((size_t) k_aug, sizeof(int));
  double *vn = (double *) R_alloc((size_t) k_aug, sizeof(double));

  double *sumZ = (p > 0) ? (double *) R_alloc((size_t) p, sizeof(double)) : NULL;
  long double *ch_sumZ = (p > 0) ? R_allocLD((size_t) p) : NULL;
  double *colmax = (p > 0) ? (double *) R_alloc((size_t) p, sizeof(double)) : NULL;
  double *anchor = (p > 0) ? (double *) R_alloc((size_t) p, sizeof(double)) : NULL;
  double *ch_anchor = (p > 0) ? (double *) R_alloc((size_t) p, sizeof(double)) : NULL;
  int *constant = (p > 0) ? (int *) R_alloc((size_t) p, sizeof(int)) : NULL;
  int *ch_const = (p > 0) ? (int *) R_alloc((size_t) p, sizeof(int)) : NULL;

  SEXP out_list_s = PROTECT(Rf_allocVector(VECSXP, S));

  for (int s = 0; s < S; ++s) {
    const int *sp = spec_int + (size_t) s * 4;
    int first = sp[0], last = sp[1], ds = sp[2], dt = sp[3];
    double L_s = spec_L[s];
    if (first < 1 || last < first || last > N || ds < 1 || ds > n_dates || dt < 1 || dt > n_dates || !R_FINITE(L_s)) {
      UNPROTECT(1);
      vmaxset(vmax);
      return R_NilValue;
    }

    int ds_cnt = db[(ds - 1) * 3 + 2], dt_cnt = db[(dt - 1) * 3 + 2];
    int s_lo = db[(ds - 1) * 3 + 0] - 1, s_hi = db[(ds - 1) * 3 + 1] - 1;
    int t_lo = db[(dt - 1) * 3 + 0] - 1, t_hi = db[(dt - 1) * 3 + 1] - 1;
    int has_pair = (ds_cnt > 0 && dt_cnt > 0);
    if (has_pair && (s_lo < 0 || s_hi < s_lo || (R_xlen_t) s_hi >= n_post ||
                     t_lo < 0 || t_hi < t_lo || (R_xlen_t) t_hi >= n_post)) {
      UNPROTECT(1);
      vmaxset(vmax);
      return R_NilValue;
    }

    double repeat_mass = 0.0, weight_mass = 0.0, sumY = 0.0;
    int n_repeat = 0, n_obs_tot = 0;
    for (int j = 0; j < p; ++j) {
      sumZ[j] = 0.0;
      colmax[j] = 0.0;
      constant[j] = 1;
      anchor[j] = 0.0;
    }
    int fac_rows = 0;
    int acc_n_qr = 0;
    int max_merge_rows = 0;
    int n_comp = 0;
    int filled = 0;
    int hs_hint = s_lo, ht_hint = t_lo;

    for (int start = first; start <= last; start += chunk_rows) {
      int64_t e64 = (int64_t) start + (int64_t) chunk_rows - 1;
      int end = (e64 > (int64_t) last) ? last : (int) e64;

      long double ch_rep_mass = 0.0L, ch_w_mass = 0.0L, ch_sumY = 0.0L;
      int ch_rep_cnt = 0, ch_obs_cnt = 0;
      for (int j = 0; j < p; ++j) ch_sumZ[j] = 0.0L;

      for (int pos = start; pos <= end; ++pos) {
        int r = coh_ord[pos - 1] - 1;
        if (r < 0 || r >= N) {
          UNPROTECT(1);
          vmaxset(vmax);
          return R_NilValue;
        }
        int ai = a_idx[r];
        if (ai > 0) {
          if (ai > n_a) {
            UNPROTECT(1);
            vmaxset(vmax);
            return R_NilValue;
          }
          double wr = w[r];
          ch_rep_mass += (long double) wr;
          ++ch_rep_cnt;

          if (has_pair) {
            int hs = cdid_find_posting(post_a, s_lo, s_hi, ai, &hs_hint);
            if (hs >= 0) {
              int ht = cdid_find_posting(post_a, t_lo, t_hi, ai, &ht_hint);
              if (ht >= 0) {
                int os = post_obs[hs], ot = post_obs[ht];
                if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y) {
                  UNPROTECT(1);
                  vmaxset(vmax);
                  return R_NilValue;
                }
                double yl = y[os - 1], yr = y[ot - 1];
                double dy = yr - yl;
                if (!R_FINITE(yl) || !R_FINITE(yr) || !R_FINITE(dy) || !R_FINITE(wr) || wr <= 0.0) {
                  UNPROTECT(1);
                  vmaxset(vmax);
                  return R_NilValue;
                }
                double sw = sqrt(wr);
                if (!R_FINITE(sw) || !R_FINITE(dy * sw)) {
                  UNPROTECT(1);
                  vmaxset(vmax);
                  return R_NilValue;
                }
                if (icpt_flag) {
                  double z_icpt = Z_a[(ai - 1) + (size_t) (icpt_col - 1) * (size_t) n_a];
                  if (z_icpt != 1.0) {
                    UNPROTECT(1);
                    vmaxset(vmax);
                    return R_NilValue;
                  }
                }
                for (int j = 0; j < p; ++j) {
                  double zj = Z_a[(ai - 1) + (size_t) (slope_cols[j] - 1) * (size_t) n_a];
                  double dj = L_s * zj;
                  if (!R_FINITE(zj) || !R_FINITE(dj) || !R_FINITE(dj * sw)) {
                    UNPROTECT(1);
                    vmaxset(vmax);
                    return R_NilValue;
                  }
                  ch_sumZ[j] += (long double) (wr * zj);
                  double az = fabs(zj);
                  if (az > colmax[j]) colmax[j] = az;
                  buf_D[filled + (size_t) j * (size_t) max_buf] = dj;
                }
                buf_dy[filled] = dy;
                buf_w[filled] = wr;
                ++filled;
                ++ch_obs_cnt;
                ch_w_mass += (long double) wr;
                ch_sumY += (long double) (wr * dy);

                if (use_qr_buf && filled == max_buf) {
                  /* Flush full QR buffer */
                  int m_fl = filled;
                  if (icpt_flag) {
                    for (int j = 0; j < p; ++j) {
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
                    for (int j = 0; j < p; ++j) {
                      aug_work[i + (size_t) (j + icpt_flag) * (size_t) m_fl] =
                          buf_D[i + (size_t) j * (size_t) max_buf] * swi;
                    }
                    aug_work[i + (size_t) (p + icpt_flag) * (size_t) m_fl] = buf_dy[i] * swi;
                  }
                  int r_ch = (m_fl < k_aug) ? m_fl : k_aug;
                  if (!cdid_qr_compress_buf(aug_work, m_fl, k_aug, fac_chunk, jpvt, vn)) {
                    UNPROTECT(1);
                    vmaxset(vmax);
                    return R_NilValue;
                  }
                  int stacked = fac_rows + r_ch;
                  if (stacked > max_merge_rows) max_merge_rows = stacked;
                  if (icpt_flag) {
                    if (acc_n_qr == 0) {
                      for (int j = 0; j < p; ++j) { anchor[j] = ch_anchor[j]; constant[j] = ch_const[j]; }
                    } else {
                      for (int j = 0; j < p; ++j) {
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
                  if (!cdid_qr_compress_buf(stack_work, stacked, k_aug, fac_merged, jpvt, vn)) {
                    UNPROTECT(1);
                    vmaxset(vmax);
                    return R_NilValue;
                  }
                  memcpy(fac_cur, fac_merged, (size_t) r_new * (size_t) k_aug * sizeof(double));
                  fac_rows = r_new;
                  acc_n_qr += m_fl;
                  n_comp += 2;
                  filled = 0;
                }
              }
            }
          }
        }
      }

      repeat_mass += (double) ch_rep_mass;
      n_repeat += ch_rep_cnt;
      if (ch_obs_cnt > 0) {
        weight_mass += (double) ch_w_mass;
        n_obs_tot += ch_obs_cnt;
        sumY += (double) ch_sumY;
        for (int j = 0; j < p; ++j) sumZ[j] += (double) ch_sumZ[j];

        if (!use_qr_buf && filled > 0) {
          int m_fl = filled;
          if (icpt_flag) {
            for (int j = 0; j < p; ++j) {
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
            for (int j = 0; j < p; ++j) {
              aug_work[i + (size_t) (j + icpt_flag) * (size_t) m_fl] =
                  buf_D[i + (size_t) j * (size_t) max_buf] * swi;
            }
            aug_work[i + (size_t) (p + icpt_flag) * (size_t) m_fl] = buf_dy[i] * swi;
          }
          int r_ch = (m_fl < k_aug) ? m_fl : k_aug;
          if (!cdid_qr_compress_buf(aug_work, m_fl, k_aug, fac_chunk, jpvt, vn)) {
            UNPROTECT(1);
            vmaxset(vmax);
            return R_NilValue;
          }
          int stacked = fac_rows + r_ch;
          if (stacked > max_merge_rows) max_merge_rows = stacked;
          if (icpt_flag) {
            if (acc_n_qr == 0) {
              for (int j = 0; j < p; ++j) { anchor[j] = ch_anchor[j]; constant[j] = ch_const[j]; }
            } else {
              for (int j = 0; j < p; ++j) {
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
          if (!cdid_qr_compress_buf(stack_work, stacked, k_aug, fac_merged, jpvt, vn)) {
            UNPROTECT(1);
            vmaxset(vmax);
            return R_NilValue;
          }
          memcpy(fac_cur, fac_merged, (size_t) r_new * (size_t) k_aug * sizeof(double));
          fac_rows = r_new;
          acc_n_qr += m_fl;
          n_comp += 2;
          filled = 0;
        }
      }
    }

    if (use_qr_buf && filled > 0) {
      int m_fl = filled;
      if (icpt_flag) {
        for (int j = 0; j < p; ++j) {
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
        for (int j = 0; j < p; ++j) {
          aug_work[i + (size_t) (j + icpt_flag) * (size_t) m_fl] =
              buf_D[i + (size_t) j * (size_t) max_buf] * swi;
        }
        aug_work[i + (size_t) (p + icpt_flag) * (size_t) m_fl] = buf_dy[i] * swi;
      }
      int r_ch = (m_fl < k_aug) ? m_fl : k_aug;
      if (!cdid_qr_compress_buf(aug_work, m_fl, k_aug, fac_chunk, jpvt, vn)) {
        UNPROTECT(1);
        vmaxset(vmax);
        return R_NilValue;
      }
      int stacked = fac_rows + r_ch;
      if (stacked > max_merge_rows) max_merge_rows = stacked;
      if (icpt_flag) {
        if (acc_n_qr == 0) {
          for (int j = 0; j < p; ++j) { anchor[j] = ch_anchor[j]; constant[j] = ch_const[j]; }
        } else {
          for (int j = 0; j < p; ++j) {
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
      if (!cdid_qr_compress_buf(stack_work, stacked, k_aug, fac_merged, jpvt, vn)) {
        UNPROTECT(1);
        vmaxset(vmax);
        return R_NilValue;
      }
      memcpy(fac_cur, fac_merged, (size_t) r_new * (size_t) k_aug * sizeof(double));
      fac_rows = r_new;
      acc_n_qr += m_fl;
      n_comp += 2;
      filled = 0;
    }

    if (!R_FINITE(repeat_mass) || !R_FINITE(weight_mass) || repeat_mass <= 0.0 ||
        weight_mass <= 0.0 || (repeat_mass / N_val) <= 0.0 || (weight_mass / repeat_mass) <= 0.0 ||
        n_obs_tot < 1 || !R_FINITE(sumY)) {
      UNPROTECT(1);
      vmaxset(vmax);
      return R_NilValue;
    }
    double rate = weight_mass / repeat_mass;
    if (!R_FINITE(rate)) {
      UNPROTECT(1);
      vmaxset(vmax);
      return R_NilValue;
    }
    for (int j = 0; j < p; ++j) {
      if (!R_FINITE(sumZ[j])) {
        UNPROTECT(1);
        vmaxset(vmax);
        return R_NilValue;
      }
    }

    SEXP fac_s = PROTECT(Rf_allocMatrix(REALSXP, fac_rows, k_aug));
    if (fac_rows > 0) memcpy(REAL(fac_s), fac_cur, (size_t) fac_rows * (size_t) k_aug * sizeof(double));
    SEXP anchor_s = PROTECT(Rf_allocVector(REALSXP, (n_obs_tot > 0) ? p : 0));
    if (n_obs_tot > 0 && p > 0) memcpy(REAL(anchor_s), anchor, (size_t) p * sizeof(double));
    SEXP const_s = PROTECT(Rf_allocVector(LGLSXP, p));
    for (int j = 0; j < p; ++j) LOGICAL(const_s)[j] = constant[j];
    SEXP sumZ_s = PROTECT(Rf_allocVector(REALSXP, p));
    if (p > 0) memcpy(REAL(sumZ_s), sumZ, (size_t) p * sizeof(double));
    SEXP colmax_s = PROTECT(Rf_allocVector(REALSXP, p));
    if (p > 0) memcpy(REAL(colmax_s), colmax, (size_t) p * sizeof(double));
    SEXP scalars_s = PROTECT(Rf_allocVector(REALSXP, 8));
    REAL(scalars_s)[0] = (double) n_obs_tot;
    REAL(scalars_s)[1] = (double) max_merge_rows;
    REAL(scalars_s)[2] = (double) n_comp;
    REAL(scalars_s)[3] = weight_mass;
    REAL(scalars_s)[4] = sumY;
    REAL(scalars_s)[5] = repeat_mass;
    REAL(scalars_s)[6] = (double) n_repeat;
    REAL(scalars_s)[7] = rate;

    SEXP item_s = PROTECT(Rf_allocVector(VECSXP, 6));
    SET_VECTOR_ELT(item_s, 0, fac_s);
    SET_VECTOR_ELT(item_s, 1, anchor_s);
    SET_VECTOR_ELT(item_s, 2, const_s);
    SET_VECTOR_ELT(item_s, 3, sumZ_s);
    SET_VECTOR_ELT(item_s, 4, colmax_s);
    SET_VECTOR_ELT(item_s, 5, scalars_s);
    SET_VECTOR_ELT(out_list_s, s, item_s);
    UNPROTECT(7);
  }

  UNPROTECT(1);
  vmaxset(vmax);
  return out_list_s;
}

/* 17. Incomplete-panel DR baseline treatment gather + IRLS logit + DR/IPW moments */
SEXP C_cdid_incomplete_dr_trt_fit(
    SEXP y_sexp,
    SEXP post_a_sexp,
    SEXP post_obs_sexp,
    SEXP bounds4_sexp,
    SEXP cohort_idx_sexp,
    SEXP a_index_sexp,
    SEXP w_sexp,
    SEXP Z_a_sexp,
    SEXP cohorts_sexp,
    SEXP support_rates_v_sexp,
    SEXP cols_list_sexp,
    SEXP reg_par_sexp,
    SEXP fit_par_sexp) {
  if (TYPEOF(y_sexp) != REALSXP || Rf_isObject(y_sexp) ||
      TYPEOF(post_a_sexp) != INTSXP || Rf_isObject(post_a_sexp) ||
      TYPEOF(post_obs_sexp) != INTSXP || Rf_isObject(post_obs_sexp) ||
      TYPEOF(bounds4_sexp) != INTSXP || XLENGTH(bounds4_sexp) != 4 ||
      TYPEOF(cohort_idx_sexp) != INTSXP || Rf_isObject(cohort_idx_sexp) ||
      TYPEOF(a_index_sexp) != INTSXP || Rf_isObject(a_index_sexp) ||
      TYPEOF(w_sexp) != REALSXP || Rf_isObject(w_sexp) ||
      TYPEOF(Z_a_sexp) != REALSXP || Rf_isObject(Z_a_sexp) ||
      TYPEOF(cohorts_sexp) != INTSXP || TYPEOF(support_rates_v_sexp) != REALSXP ||
      TYPEOF(cols_list_sexp) != VECSXP || XLENGTH(cols_list_sexp) != 2 ||
      TYPEOF(reg_par_sexp) != VECSXP || XLENGTH(reg_par_sexp) != 2 ||
      TYPEOF(fit_par_sexp) != REALSXP || XLENGTH(fit_par_sexp) != 5) {
    return R_NilValue;
  }

  SEXP zdim = Rf_getAttrib(Z_a_sexp, R_DimSymbol);
  if (TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return R_NilValue;
  int n_a = INTEGER(zdim)[0];
  int p_Z = INTEGER(zdim)[1];
  int N = (int) XLENGTH(cohort_idx_sexp);
  int m = (int) XLENGTH(cohorts_sexp);
  R_xlen_t n_post = XLENGTH(post_a_sexp);
  R_xlen_t n_y = XLENGTH(y_sexp);

  SEXP trt_cols_s = VECTOR_ELT(cols_list_sexp, 0);
  SEXP or_cols_s = VECTOR_ELT(cols_list_sexp, 1);
  SEXP reg_la_s = VECTOR_ELT(reg_par_sexp, 0);
  SEXP reg_beta_s = VECTOR_ELT(reg_par_sexp, 1);
  if (TYPEOF(trt_cols_s) != INTSXP || TYPEOF(or_cols_s) != INTSXP ||
      TYPEOF(reg_la_s) != REALSXP || XLENGTH(reg_la_s) != 2 || TYPEOF(reg_beta_s) != REALSXP) {
    return R_NilValue;
  }

  int p_trt = (int) XLENGTH(trt_cols_s);
  int p_or = (int) XLENGTH(or_cols_s);
  double N_val = REAL(fit_par_sexp)[0];
  int n_obs_exp = (int) REAL(fit_par_sexp)[2];
  double eps = REAL(fit_par_sexp)[3];
  int maxit = (int) REAL(fit_par_sexp)[4];

  if (n_a < 1 || p_Z < 1 || N < 1 || m < 2 || p_trt < 1 || p_trt > 128 || p_or < 0 ||
      n_obs_exp < p_trt || !R_FINITE(N_val) || N_val <= 0.0 || !R_FINITE(eps) || eps <= 0.0 || maxit < 1 ||
      XLENGTH(a_index_sexp) != (R_xlen_t) N || XLENGTH(w_sexp) != (R_xlen_t) N ||
      XLENGTH(post_obs_sexp) != n_post || XLENGTH(support_rates_v_sexp) != (R_xlen_t) (2 * m) ||
      XLENGTH(reg_beta_s) != (R_xlen_t) p_or) {
    return R_NilValue;
  }

  const int *b4 = INTEGER(bounds4_sexp);
  int s_lo = b4[0] - 1, s_hi = b4[1] - 1, t_lo = b4[2] - 1, t_hi = b4[3] - 1;
  if (b4[0] < 1 || s_hi < s_lo || (R_xlen_t) s_hi >= n_post ||
      b4[2] < 1 || t_hi < t_lo || (R_xlen_t) t_hi >= n_post) {
    return R_NilValue;
  }

  const int *trt_cols = INTEGER(trt_cols_s);
  const int *or_cols = INTEGER(or_cols_s);
  for (int j = 0; j < p_trt; ++j) if (trt_cols[j] < 1 || trt_cols[j] > p_Z) return R_NilValue;
  for (int j = 0; j < p_or; ++j) if (or_cols[j] < 1 || or_cols[j] > p_Z) return R_NilValue;

  const double *y_obs_vec = REAL(y_sexp);
  const int *post_a = INTEGER(post_a_sexp);
  const int *post_obs = INTEGER(post_obs_sexp);
  const int *coh_idx = INTEGER(cohort_idx_sexp);
  const int *a_idx = INTEGER(a_index_sexp);
  const double *w_frame = REAL(w_sexp);
  const double *Z_a = REAL(Z_a_sexp);
  const int *cohorts = INTEGER(cohorts_sexp);
  const double *rates_vec = REAL(support_rates_v_sexp);
  const double *v_vec = rates_vec + m;
  double reg_L = REAL(reg_la_s)[0];
  double reg_alpha = REAL(reg_la_s)[1];
  const double *or_beta = REAL(reg_beta_s);

  int max_coh = 0;
  for (int k = 0; k < m; ++k) {
    if (cohorts[k] < 1) return R_NilValue;
    if (cohorts[k] > max_coh) max_coh = cohorts[k];
  }
  if (max_coh > 1000000) return R_NilValue;

  const void *vmax = vmaxget();
  int *cmap = (int *) R_alloc((size_t) max_coh, sizeof(int));
  memset(cmap, 0, (size_t) max_coh * sizeof(int));
  for (int k = 0; k < m; ++k) cmap[cohorts[k] - 1] = k + 1;

  int n = n_obs_exp;
  double *X_trt = (double *) R_alloc((size_t) n * (size_t) p_trt, sizeof(double));
  double *Z_or = (p_or > 0) ? (double *) R_alloc((size_t) n * (size_t) p_or, sizeof(double)) : NULL;
  double *y_bin = (double *) R_alloc((size_t) n, sizeof(double));
  double *w_fit = (double *) R_alloc((size_t) n, sizeof(double));
  double *dy_vec = (double *) R_alloc((size_t) n, sizeof(double));
  int *kk_vec = (int *) R_alloc((size_t) n, sizeof(int));

  int pos = 0;
  int hs_hint = s_lo, ht_hint = t_lo;
  for (int r = 0; r < N; ++r) {
    int c = coh_idx[r];
    if (c >= 1 && c <= max_coh) {
      int ci = cmap[c - 1];
      if (ci > 0) {
        int ai = a_idx[r];
        if (ai >= 1 && ai <= n_a) {
          int hs = cdid_find_posting(post_a, s_lo, s_hi, ai, &hs_hint);
          if (hs >= 0) {
            int ht = cdid_find_posting(post_a, t_lo, t_hi, ai, &ht_hint);
            if (ht >= 0) {
              if (pos >= n) {
                vmaxset(vmax);
                return R_NilValue;
              }
              int os = post_obs[hs], ot = post_obs[ht];
              if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y) {
                vmaxset(vmax);
                return R_NilValue;
              }
              double yl = y_obs_vec[os - 1], yr = y_obs_vec[ot - 1];
              double dy = yr - yl;
              if (!R_FINITE(yl) || !R_FINITE(yr) || !R_FINITE(dy)) {
                vmaxset(vmax);
                return R_NilValue;
              }
              for (int j = 0; j < p_trt; ++j) {
                X_trt[pos + (size_t) j * (size_t) n] =
                    Z_a[(ai - 1) + (size_t) (trt_cols[j] - 1) * (size_t) n_a];
              }
              for (int j = 0; j < p_or; ++j) {
                Z_or[pos + (size_t) j * (size_t) n] =
                    Z_a[(ai - 1) + (size_t) (or_cols[j] - 1) * (size_t) n_a];
              }
              y_bin[pos] = (ci == 1) ? 1.0 : 0.0;
              kk_vec[pos] = ci;
              dy_vec[pos] = dy;
              w_fit[pos] = w_frame[r] * v_vec[ci - 1];
              ++pos;
            }
          }
        }
      }
    }
  }
  if (pos != n) {
    vmaxset(vmax);
    return R_NilValue;
  }

  /* Validate weights & X_trt + separation check */
  int n_pos_w = 0;
  long double sum_w_ld = 0.0L;
  long double devold_ld = 0.0L;
  const double ln2_val = log(2.0);
  for (int i = 0; i < n; ++i) {
    double wi = w_fit[i];
    if (!R_FINITE(wi) || wi < 0.0) {
      vmaxset(vmax);
      return R_NilValue;
    }
    if (wi > 0.0) {
      ++n_pos_w;
      sum_w_ld += (long double) wi;
      devold_ld += (long double) (2.0 * wi * ln2_val);
    }
  }
  double sum_w = (double) sum_w_ld;
  if (n_pos_w < p_trt || !R_FINITE(sum_w) || sum_w <= 0.0) {
    vmaxset(vmax);
    return R_NilValue;
  }
  for (int j = 0; j < p_trt; ++j) {
    const double *xj = X_trt + (size_t) j * (size_t) n;
    int has_pos = 0, has_neg = 0;
    for (int i = 0; i < n; ++i) {
      double v = xj[i];
      if (!R_FINITE(v)) {
        vmaxset(vmax);
        return R_NilValue;
      }
      if (w_fit[i] > 0.0 && v != 0.0) {
        double sgn = (y_bin[i] == 1.0) ? v : -v;
        if (sgn > 0.0) has_pos = 1;
        else if (sgn < 0.0) has_neg = 1;
      }
    }
    if ((has_pos && !has_neg) || (has_neg && !has_pos)) {
      vmaxset(vmax);
      return R_NilValue;
    }
  }

  /* IRLS logit matching C_cdid_cp04_logit_fast */
  double *eta = (double *) R_alloc((size_t) n, sizeof(double));
  double *mu = (double *) R_alloc((size_t) n, sizeof(double));
  double *wt = (double *) R_alloc((size_t) n, sizeof(double));
  double *wz = (double *) R_alloc((size_t) n, sizeof(double));
  double *W_X = (double *) R_alloc((size_t) n * (size_t) p_trt, sizeof(double));
  double *H = (double *) R_alloc((size_t) p_trt * (size_t) p_trt, sizeof(double));
  double *R_chol = (double *) R_alloc((size_t) p_trt * (size_t) p_trt, sizeof(double));
  double *beta = (double *) R_alloc((size_t) p_trt, sizeof(double));

  for (int i = 0; i < n; ++i) { eta[i] = 0.0; mu[i] = 0.5; }
  for (int j = 0; j < p_trt; ++j) beta[j] = 0.0;

  double devold = (double) devold_ld;
  int conv = 0, iter_used = 0;
  const double bnd_eps = 10.0 * DBL_EPSILON;
  const double qr_tol = (eps / 1000.0 < 1e-7) ? (eps / 1000.0) : 1e-7;
  int n_b = n, ptrt_b = p_trt, one_b = 1;
  double one_d = 1.0, zero_d = 0.0;

  for (int iter = 1; iter <= maxit; ++iter) {
    iter_used = iter;
    for (int i = 0; i < n; ++i) {
      double wi = w_fit[i];
      if (wi > 0.0) {
        double mui = mu[i];
        double varmu = mui * (1.0 - mui);
        double opexp = exp(eta[i]);
        double mu_eta = opexp / ((1.0 + opexp) * (1.0 + opexp));
        if (!(varmu > 0.0) || !(mu_eta > 0.0) || !R_FINITE(mu_eta)) { vmaxset(vmax); return R_NilValue; }
        double w_i = sqrt((wi * (mu_eta * mu_eta)) / varmu);
        double z_i = eta[i] + (y_bin[i] - mui) / mu_eta;
        wt[i] = w_i;
        wz[i] = z_i * w_i;
      } else {
        wt[i] = 0.0;
        wz[i] = 0.0;
      }
    }
    for (int j = 0; j < p_trt; ++j) {
      const double *xj = X_trt + (size_t) j * (size_t) n;
      double *wxj = W_X + (size_t) j * (size_t) n;
      for (int i = 0; i < n; ++i) wxj[i] = xj[i] * wt[i];
    }
    if (!cdid_qr_wls_solve(W_X, wz, R_chol, beta, n, p_trt, qr_tol, 0.0, 0)) { vmaxset(vmax); return R_NilValue; }

    F77_CALL(dgemm)("N", "N", &n_b, &one_b, &ptrt_b, &one_d,
                    X_trt, &n_b, beta, &ptrt_b, &zero_d, eta, &n_b FCONE FCONE);

    long double dev_ld = 0.0L;
    for (int i = 0; i < n; ++i) {
      double ei = eta[i];
      if (!R_FINITE(ei) || fabs(ei) >= 30.0) { vmaxset(vmax); return R_NilValue; }
      double opexp = exp(ei);
      double mui = opexp / (1.0 + opexp);
      if (mui <= bnd_eps || mui >= 1.0 - bnd_eps) { vmaxset(vmax); return R_NilValue; }
      mu[i] = mui;
      if (w_fit[i] > 0.0) {
        double loss = (y_bin[i] == 1.0) ? log(1.0 / mui) : log(1.0 / (1.0 - mui));
        dev_ld += (long double) (2.0 * w_fit[i] * loss);
      }
    }
    double dev = (double) dev_ld;
    if (!R_FINITE(dev)) { vmaxset(vmax); return R_NilValue; }
    if (fabs(dev - devold) / (0.1 + fabs(dev)) < eps) { conv = 1; break; }
    devold = dev;
  }
  if (!conv) { vmaxset(vmax); return R_NilValue; }

  /* Final information matrix H_info and clip_logit values */
  double *prob = mu;
  double *odds = wz;
  double *act = eta; /* reuse after reading eta[i] */
  const double clip_cap = 1.0 - 1e-6;
  const double clip_cap_odds = clip_cap / (1.0 - clip_cap);
  double p_min = DBL_MAX, p_max = -DBL_MAX, w_max = -DBL_MAX;
  int p_capped = 0, p_capped_ctrl = 0;

  for (int i = 0; i < n; ++i) {
    double ei = eta[i];
    double od, pr;
    if (ei > 40.0) {
      pr = 1.0 / (1.0 + exp(-ei));
      od = 1.0;
    } else {
      od = exp(ei);
      pr = od / (1.0 + od);
    }
    if (pr <= bnd_eps || pr >= 1.0 - bnd_eps || !(pr * (1.0 - pr) > 0.0)) {
      vmaxset(vmax);
      return R_NilValue;
    }
    wt[i] = w_fit[i] * pr * (1.0 - pr);
    if (pr < clip_cap) {
      act[i] = 1.0;
      odds[i] = od;
    } else {
      act[i] = 0.0;
      odds[i] = clip_cap_odds;
      ++p_capped;
      if (y_bin[i] == 0.0) ++p_capped_ctrl;
    }
    prob[i] = pr;
    if (pr < p_min) p_min = pr;
    if (pr > p_max) p_max = pr;
    if (w_fit[i] > w_max) w_max = w_fit[i];
  }

  for (int j = 0; j < p_trt; ++j) {
    const double *xj = X_trt + (size_t) j * (size_t) n;
    double *wxj = W_X + (size_t) j * (size_t) n;
    for (int i = 0; i < n; ++i) wxj[i] = wt[i] * xj[i];
  }
  F77_CALL(dgemm)("T", "N", &ptrt_b, &ptrt_b, &n_b, &one_d,
                  X_trt, &n_b, W_X, &n_b, &zero_d, H, &ptrt_b FCONE FCONE);
  for (int idx = 0; idx < p_trt * p_trt; ++idx) {
    H[idx] /= N_val;
    if (!R_FINITE(H[idx])) { vmaxset(vmax); return R_NilValue; }
  }
  if (!cdid_chol_upper(H, R_chol, p_trt, 0.0)) { vmaxset(vmax); return R_NilValue; }

  double *wp = (double *) R_alloc((size_t) n, sizeof(double));
  double *wm = (double *) R_alloc((size_t) n, sizeof(double));
  double *r_obs = (double *) R_alloc((size_t) n, sizeof(double));
  double *lin_or_vec = (double *) R_alloc((size_t) n, sizeof(double));
  if (p_or > 0) {
    double *D_or = (double *) R_alloc((size_t) n * (size_t) p_or, sizeof(double));
    for (size_t idx = 0; idx < (size_t) n * (size_t) p_or; ++idx) D_or[idx] = reg_L * Z_or[idx];
    int por_b = p_or;
    F77_CALL(dgemm)("N", "N", &n_b, &one_b, &por_b, &one_d,
                    D_or, &n_b, or_beta, &por_b, &zero_d, lin_or_vec, &n_b FCONE FCONE);
  } else {
    memset(lin_or_vec, 0, (size_t) n * sizeof(double));
  }

  long double sum_wp_ld = 0.0L, sum_wm_ld = 0.0L;
  long double tipw_ld = 0.0L, mipw_ld = 0.0L, tdr_ld = 0.0L, mdr_ld = 0.0L;

#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#endif
  for (int i = 0; i < n; ++i) {
    double wpi = w_fit[i] * y_bin[i];
    double wmi = w_fit[i] * (1.0 - y_bin[i]) * odds[i];
    wp[i] = wpi;
    wm[i] = wmi;
    sum_wp_ld += (long double) wpi;
    sum_wm_ld += (long double) wmi;

    double mhat = reg_alpha + lin_or_vec[i];
    double ri = dy_vec[i] - mhat;
    r_obs[i] = ri;
    tipw_ld += (long double) (wpi * dy_vec[i]);
    mipw_ld += (long double) (wmi * dy_vec[i]);
    tdr_ld += (long double) (wpi * ri);
    mdr_ld += (long double) (wmi * ri);
  }
  double sum_wp = (double) sum_wp_ld;
  double sum_wm = (double) sum_wm_ld;

  double dp = sum_wp / N_val;
  double dm = sum_wm / N_val;
  if (!R_FINITE(dp) || dp <= 0.0 || !R_FINITE(dm) || dm <= 0.0) {
    vmaxset(vmax);
    return R_NilValue;
  }

  double T_ipw = ((double) tipw_ld) / sum_wp;
  double M_ipw = ((double) mipw_ld) / sum_wm;
  double mu_ipw = T_ipw - M_ipw;
  double T_dr = ((double) tdr_ld) / sum_wp;
  double M_dr = ((double) mdr_ld) / sum_wm;
  double mu_dr = T_dr - M_dr;

  SEXP dbp_s = PROTECT(Rf_allocVector(REALSXP, p_or));
  SEXP dbm_s = PROTECT(Rf_allocVector(REALSXP, p_or));
  SEXP dbd_s = PROTECT(Rf_allocVector(REALSXP, p_or));
  double *dbp = REAL(dbp_s), *dbm = REAL(dbm_s), *dbd = REAL(dbd_s);
  long double dot_diff_ld = 0.0L;
  for (int j = 0; j < p_or; ++j) {
    const double *zj = Z_or + (size_t) j * (size_t) n;
    long double sp_j = 0.0L, sm_j = 0.0L;
#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#endif
    for (int i = 0; i < n; ++i) {
      double dij = reg_L * zj[i];
      sp_j += (long double) (wp[i] * dij);
      sm_j += (long double) (wm[i] * dij);
    }
    dbp[j] = ((double) sp_j) / sum_wp;
    dbm[j] = ((double) sm_j) / sum_wm;
    dbd[j] = dbp[j] - dbm[j];
    dot_diff_ld += (long double) (dbd[j] * or_beta[j]);
    if (!R_FINITE(dbp[j]) || !R_FINITE(dbm[j])) {
      UNPROTECT(3);
      vmaxset(vmax);
      return R_NilValue;
    }
  }
  double mu_dr_decomp = mu_ipw - ((double) dot_diff_ld);
  if (!R_FINITE(T_ipw) || !R_FINITE(M_ipw) || !R_FINITE(mu_ipw) ||
      !R_FINITE(T_dr) || !R_FINITE(M_dr) || !R_FINITE(mu_dr) || !R_FINITE(mu_dr_decomp)) {
    UNPROTECT(3);
    vmaxset(vmax);
    return R_NilValue;
  }

  SEXP C_s = PROTECT(Rf_allocMatrix(REALSXP, p_trt, m));
  SEXP Udr_s = PROTECT(Rf_allocMatrix(REALSXP, 4, m));
  SEXP Uipw_s = PROTECT(Rf_allocMatrix(REALSXP, 4, m));
  SEXP Vdr_s = PROTECT(Rf_allocMatrix(REALSXP, 4, p_trt));
  SEXP Vipw_s = PROTECT(Rf_allocMatrix(REALSXP, 4, p_trt));
  double *C_mat = REAL(C_s);
  double *Udr = REAL(Udr_s);
  double *Uipw = REAL(Uipw_s);
  double *Vdr = REAL(Vdr_s);
  double *Vipw = REAL(Vipw_s);
  memset(C_mat, 0, (size_t) p_trt * (size_t) m * sizeof(double));
  memset(Udr, 0, (size_t) 4 * (size_t) m * sizeof(double));
  memset(Uipw, 0, (size_t) 4 * (size_t) m * sizeof(double));
  memset(Vdr, 0, (size_t) 4 * (size_t) p_trt * sizeof(double));
  memset(Vipw, 0, (size_t) 4 * (size_t) p_trt * sizeof(double));

  long double *C_ld = R_allocLD((size_t) p_trt * (size_t) m);
  long double *Udr_ld = R_allocLD((size_t) 4 * (size_t) m);
  long double *Uipw_ld = R_allocLD((size_t) 4 * (size_t) m);
  for (int idx = 0; idx < p_trt * m; ++idx) C_ld[idx] = 0.0L;
  for (int idx = 0; idx < 4 * m; ++idx) { Udr_ld[idx] = 0.0L; Uipw_ld[idx] = 0.0L; }

  const double clip_eps = 1e-6;
#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#endif
  for (int i = 0; i < n; ++i) {
    int k0 = kk_vec[i] - 1;
    double sc_base = w_fit[i] * (y_bin[i] - prob[i]);
    long double *ccol = C_ld + (size_t) k0 * (size_t) p_trt;
    for (int j = 0; j < p_trt; ++j) {
      ccol[j] += (long double) (sc_base * X_trt[i + (size_t) j * (size_t) n]);
    }
    long double *ud = Udr_ld + (size_t) k0 * 4;
    ud[0] += (long double) (wp[i] * (r_obs[i] - T_dr) / dp);
    ud[1] += (long double) (wm[i] * (r_obs[i] - M_dr) / dm);
    ud[2] += (long double) wp[i];
    ud[3] += (long double) wm[i];

    long double *ui = Uipw_ld + (size_t) k0 * 4;
    ui[0] += (long double) (wp[i] * (dy_vec[i] - T_ipw) / dp);
    ui[1] += (long double) (wm[i] * (dy_vec[i] - M_ipw) / dm);
    ui[2] += (long double) wp[i];
    ui[3] += (long double) wm[i];
  }

  for (int k0 = 0; k0 < m; ++k0) {
    double rk = rates_vec[k0];
    int active_k = (rk > clip_eps);
    double scale_k = active_k ? (-1.0 / (N_val * rk)) : 0.0;
    double *ccol = C_mat + (size_t) k0 * (size_t) p_trt;
    const long double *cld = C_ld + (size_t) k0 * (size_t) p_trt;
    for (int j = 0; j < p_trt; ++j) {
      ccol[j] = active_k ? (((double) cld[j]) * scale_k) : 0.0;
      if (!R_FINITE(ccol[j])) { UNPROTECT(8); vmaxset(vmax); return R_NilValue; }
    }
    double *ud = Udr + (size_t) k0 * 4;
    double *ui = Uipw + (size_t) k0 * 4;
    const long double *udl = Udr_ld + (size_t) k0 * 4;
    const long double *uil = Uipw_ld + (size_t) k0 * 4;
    for (int r = 0; r < 4; ++r) {
      ud[r] = active_k ? (((double) udl[r]) * scale_k) : 0.0;
      ui[r] = active_k ? (((double) uil[r]) * scale_k) : 0.0;
      if (!R_FINITE(ud[r]) || !R_FINITE(ui[r])) { UNPROTECT(8); vmaxset(vmax); return R_NilValue; }
    }
  }

  double ndm = N_val * dm;
  for (int j = 0; j < p_trt; ++j) {
    const double *xj = X_trt + (size_t) j * (size_t) n;
    long double s_wmd = 0.0L, s_dr2 = 0.0L, s_ipw2 = 0.0L;
#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#endif
    for (int i = 0; i < n; ++i) {
      double wm_act_i = wm[i] * act[i];
      s_wmd += (long double) (wm_act_i * xj[i]);
      s_dr2 += (long double) ((wm_act_i * (r_obs[i] - M_dr)) * xj[i]);
      s_ipw2 += (long double) ((wm_act_i * (dy_vec[i] - M_ipw)) * xj[i]);
    }
    double v4 = ((double) s_wmd) / N_val;
    double vdr2 = ((double) s_dr2) / ndm;
    double vipw2 = ((double) s_ipw2) / ndm;
    if (!R_FINITE(v4) || !R_FINITE(vdr2) || !R_FINITE(vipw2)) {
      UNPROTECT(8);
      vmaxset(vmax);
      return R_NilValue;
    }
    Vdr[1 + j * 4] = vdr2;
    Vdr[3 + j * 4] = v4;
    Vipw[1 + j * 4] = vipw2;
    Vipw[3 + j * 4] = v4;
  }

  SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, p_trt));
  memcpy(REAL(beta_s), beta, (size_t) p_trt * sizeof(double));
  SEXP R_s = PROTECT(Rf_allocMatrix(REALSXP, p_trt, p_trt));
  memcpy(REAL(R_s), R_chol, (size_t) p_trt * (size_t) p_trt * sizeof(double));
  SEXP means_s = PROTECT(Rf_allocVector(REALSXP, 12));
  double *mv = REAL(means_s);
  mv[0] = T_dr; mv[1] = M_dr; mv[2] = dp; mv[3] = dm;
  mv[4] = T_dr; mv[5] = M_dr; mv[6] = mu_dr;
  mv[7] = T_ipw; mv[8] = M_ipw; mv[9] = mu_ipw;
  mv[10] = T_dr; mv[11] = mu_dr_decomp;

  SEXP clip_s = PROTECT(Rf_allocVector(REALSXP, 6));
  REAL(clip_s)[0] = p_min;
  REAL(clip_s)[1] = p_max;
  REAL(clip_s)[2] = (double) n;
  REAL(clip_s)[3] = (double) p_capped;
  REAL(clip_s)[4] = (double) p_capped_ctrl;
  REAL(clip_s)[5] = w_max;

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 13));
  SET_VECTOR_ELT(out_s, 0, beta_s);
  SET_VECTOR_ELT(out_s, 1, R_s);
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarInteger(iter_used));
  SET_VECTOR_ELT(out_s, 3, C_s);
  SET_VECTOR_ELT(out_s, 4, Udr_s);
  SET_VECTOR_ELT(out_s, 5, Vdr_s);
  SET_VECTOR_ELT(out_s, 6, Uipw_s);
  SET_VECTOR_ELT(out_s, 7, Vipw_s);
  SET_VECTOR_ELT(out_s, 8, means_s);
  SET_VECTOR_ELT(out_s, 9, dbp_s);
  SET_VECTOR_ELT(out_s, 10, dbm_s);
  SET_VECTOR_ELT(out_s, 11, dbd_s);
  SET_VECTOR_ELT(out_s, 12, clip_s);

  UNPROTECT(13);
  vmaxset(vmax);
  return out_s;
}

/* ===========================================================================
   WAVE 2.5 KERNELS: Fast Seals, Pointer/Seal Verifiers, Retained-IF Gram,
                     Batch DR Cholesky Solve & Diagnostic Summaries
   =========================================================================== */

static R_INLINE SEXP cdid_env_get_checked_cached(SEXP env, SEXP sym, int *bound_out,
                                                int *locked_out, SEXP cache) {
  if (TYPEOF(env) != ENVSXP) {
    *bound_out = 0;
    *locked_out = 0;
    return R_UnboundValue;
  }
  int bound = 0;
  SEXP val = cdid_frame_value_cached(env, sym, &bound, cache);
  *bound_out = bound;
  if (!bound || val == R_UnboundValue) {
    *locked_out = 0;
    return R_UnboundValue;
  }
  *locked_out = (R_BindingIsLocked(sym, env) != 0);
  return val;
}

static R_INLINE SEXP cdid_env_get_checked(SEXP env, SEXP sym, int *bound_out, int *locked_out) {
  return cdid_env_get_checked_cached(env, sym, bound_out, locked_out, NULL);
}

static SEXP cdid_default_fns_match(SEXP names, SEXP owner, SEXP snapshot, SEXP cache) {
  SEXP trace_class = PRINTNAME(Rf_install("functionWithTrace"));
  for (R_xlen_t i = 0; i < XLENGTH(names); ++i) {
    SEXP s = STRING_ELT(names, i);
    if (s == NA_STRING) return Rf_ScalarLogical(0);
    SEXP sym = Rf_installChar(s);
    int bound = 0, locked = 0;
    SEXP live = cdid_env_get_checked_cached(owner, sym, &bound, &locked, cache);
    if (!bound || live == R_UnboundValue || TYPEOF(live) == PROMSXP) return R_NilValue;
    if (!Rf_isFunction(live)) return Rf_ScalarLogical(0);
    if (Rf_isObject(live)) {
      SEXP cls = Rf_getAttrib(live, R_ClassSymbol);
      if (TYPEOF(cls) == STRSXP) {
        for (R_xlen_t j = 0; j < XLENGTH(cls); ++j) {
          if (STRING_ELT(cls, j) == trace_class) return Rf_ScalarLogical(0);
        }
      }
    }
#if R_VERSION < R_Version(4, 6, 0)
    PROTECT(live);
#endif
    SEXP expected = cdid_env_get_checked(snapshot, sym, &bound, &locked);
    int matches = bound && expected != R_UnboundValue && TYPEOF(expected) != PROMSXP &&
                  live == expected;
#if R_VERSION < R_Version(4, 6, 0)
    UNPROTECT(1);
#endif
    if (!matches) return R_NilValue;
  }
  return Rf_ScalarLogical(1);
}

static R_INLINE int cdid_c_seal_ok(SEXP obj, SEXP tag_sym) {
  if (TYPEOF(obj) != VECSXP) return 0;
  SEXP ref = cdid_list_get_named(obj, "ref");
  if (ref == R_UnboundValue || TYPEOF(ref) != ENVSXP) return 0;
  if (cdid_parent_env(ref) != R_EmptyEnv || !R_EnvironmentIsLocked(ref)) return 0;
  int bound = 0, locked = 0;
  SEXP val = cdid_env_get_checked(ref, tag_sym, &bound, &locked);
  return (bound && locked && val == obj) ? 1 : 0;
}

static R_INLINE int cdid_scalar_str_eq(SEXP x, const char *expected) {
  if (TYPEOF(x) != STRSXP || XLENGTH(x) != 1 || Rf_isObject(x)) return 0;
  SEXP s = STRING_ELT(x, 0);
  return (s != NA_STRING && strcmp(CHAR(s), expected) == 0) ? 1 : 0;
}

static R_INLINE int cdid_strsxp_eq(SEXP a, SEXP b) {
  if (a == b) return (TYPEOF(a) == STRSXP);
  if (TYPEOF(a) != STRSXP || TYPEOF(b) != STRSXP) return 0;
  if (XLENGTH(a) != XLENGTH(b)) return 0;
  R_xlen_t n = XLENGTH(a);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (STRING_ELT(a, i) != STRING_ELT(b, i)) return 0;
  }
  return 1;
}

static R_INLINE int cdid_as_scalar_int(SEXP x, int *out) {
  if (Rf_isObject(x)) return 0;
  if (TYPEOF(x) == INTSXP && XLENGTH(x) == 1) {
    int v = INTEGER(x)[0];
    if (v == NA_INTEGER) return 0;
    *out = v;
    return 1;
  }
  if (TYPEOF(x) == REALSXP && XLENGTH(x) == 1) {
    double d = REAL(x)[0];
    if (!R_FINITE(d) || d != floor(d) || d < (double) INT_MIN || d > (double) INT_MAX) return 0;
    *out = (int) d;
    return 1;
  }
  return 0;
}

SEXP C_cdid_fast_seal(SEXP obj, SEXP tag_sexp) {
  if (TYPEOF(obj) != VECSXP || TYPEOF(tag_sexp) != STRSXP || XLENGTH(tag_sexp) != 1) {
    return R_NilValue;
  }
  SEXP tag_chr = STRING_ELT(tag_sexp, 0);
  if (tag_chr == NA_STRING) return R_NilValue;
  SEXP nms = Rf_getAttrib(obj, R_NamesSymbol);
  if (TYPEOF(nms) != STRSXP || XLENGTH(nms) != XLENGTH(obj)) return R_NilValue;
  R_xlen_t ref_idx = -1;
  R_xlen_t n = XLENGTH(obj);
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP s = STRING_ELT(nms, i);
    if (s != NA_STRING && strcmp(CHAR(s), "ref") == 0) {
      ref_idx = i;
      break;
    }
  }
  int nprot = 0;
  if (ref_idx < 0) {
    SEXP new_obj = PROTECT(Rf_allocVector(VECSXP, n + 1));
    SEXP new_nms = PROTECT(Rf_allocVector(STRSXP, n + 1));
    nprot += 2;
    for (R_xlen_t i = 0; i < n; ++i) {
      SET_VECTOR_ELT(new_obj, i, VECTOR_ELT(obj, i));
      SET_STRING_ELT(new_nms, i, STRING_ELT(nms, i));
    }
    SET_STRING_ELT(new_nms, n, Rf_mkChar("ref"));
    Rf_setAttrib(new_obj, R_NamesSymbol, new_nms);
    SEXP cls = Rf_getAttrib(obj, R_ClassSymbol);
    if (cls != R_NilValue) Rf_setAttrib(new_obj, R_ClassSymbol, cls);
    obj = new_obj;
    ref_idx = n;
  } else if (MAYBE_SHARED(obj)) {
    obj = PROTECT(Rf_shallow_duplicate(obj));
    ++nprot;
  }
  SEXP ref = PROTECT(R_NewEnv(R_EmptyEnv, FALSE, 0));
  ++nprot;
  SET_VECTOR_ELT(obj, ref_idx, ref);
  SEXP sym = Rf_installChar(tag_chr);
  Rf_defineVar(sym, obj, ref);
  R_LockEnvironment(ref, TRUE);
  UNPROTECT(nprot);
  return obj;
}

static int cdid_is_already_clean(SEXP x) {
  if (x == R_NilValue || Rf_isObject(x)) return 0;
  SEXPTYPE t = TYPEOF(x);
  if (t == LGLSXP || t == INTSXP || t == REALSXP || t == CPLXSXP || t == STRSXP || t == RAWSXP) {
    if (!cdid_has_attrib(x)) return 1;
    const SEXP allowed[3] = {R_NamesSymbol, R_DimSymbol, R_DimNamesSymbol};
    return cdid_attrib_only(x, allowed, 3);
  }
  if (t == VECSXP) {
    if (cdid_has_attrib(x)) {
      const SEXP allowed[1] = {R_NamesSymbol};
      if (!cdid_attrib_only(x, allowed, 1)) return 0;
    }
    R_xlen_t n = XLENGTH(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      if (!cdid_is_already_clean(VECTOR_ELT(x, i))) return 0;
    }
    return 1;
  }
  return 0;
}

static SEXP cdid_scalar_state_c(SEXP x, int *ok) {
  if (x == R_NilValue) {
    *ok = 1;
    return R_NilValue;
  }
  if (cdid_is_already_clean(x)) {
    *ok = 1;
    return x;
  }
  SEXPTYPE t = TYPEOF(x);
  if (t == ENVSXP || t == EXTPTRSXP || t == CLOSXP || t == BUILTINSXP || t == SPECIALSXP || Rf_isFunction(x)) {
    *ok = 1;
    return R_NilValue;
  }
  if (t == VECSXP && !Rf_inherits(x, "data.frame")) {
    R_xlen_t n = XLENGTH(x);
    SEXP tmp = PROTECT(Rf_allocVector(VECSXP, n));
    R_xlen_t kept = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
      int sub_ok = 0;
      SEXP c_elt = cdid_scalar_state_c(VECTOR_ELT(x, i), &sub_ok);
      if (!sub_ok) {
        UNPROTECT(1);
        *ok = 0;
        return R_NilValue;
      }
      if (c_elt != R_NilValue) {
        SET_VECTOR_ELT(tmp, i, c_elt);
        ++kept;
      }
    }
    SEXP nms = Rf_getAttrib(x, R_NamesSymbol);
    int has_nms = (nms != R_NilValue && TYPEOF(nms) == STRSXP && XLENGTH(nms) == n);
    SEXP out = PROTECT(Rf_allocVector(VECSXP, kept));
    SEXP out_nms = has_nms ? PROTECT(Rf_allocVector(STRSXP, kept)) : R_NilValue;
    R_xlen_t pos = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
      SEXP c_elt = VECTOR_ELT(tmp, i);
      if (c_elt != R_NilValue) {
        SET_VECTOR_ELT(out, pos, c_elt);
        if (has_nms) SET_STRING_ELT(out_nms, pos, STRING_ELT(nms, i));
        ++pos;
      }
    }
    if (has_nms) {
      Rf_setAttrib(out, R_NamesSymbol, out_nms);
      UNPROTECT(3);
    } else {
      UNPROTECT(2);
    }
    *ok = 1;
    return out;
  }
  *ok = 0;
  return R_NilValue;
}

/* Multi-mode fast C pointer/environment-lock verifier.
   Returns ScalarLogical(1) / list when 100% verified; returns ScalarLogical(0) / R_NilValue
   whenever any condition does not match so the caller runs the full R validator for exact
   refusal semantics. */
SEXP C_cdid_verify_fast(SEXP mode_sexp, SEXP a1, SEXP a2, SEXP a3, SEXP a4) {
  int mode = 0;
  if (!cdid_as_scalar_int(mode_sexp, &mode)) return R_NilValue;

  if (mode == 1) {
    /* custom_ipw_if_plan_verify(plan = a1, control = a2) */
    SEXP plan = a1;
    SEXP control = a2;
    if (TYPEOF(plan) != VECSXP) return Rf_ScalarLogical(0);
    SEXP kind = cdid_list_get_named(plan, "kind");
    if (!cdid_scalar_str_eq(kind, "cdid_custom_ipw_if_plan")) return Rf_ScalarLogical(0);
    if (!cdid_c_seal_ok(plan, Rf_install("custom_ipw_if_plan"))) return Rf_ScalarLogical(0);

    SEXP ctrl_in = cdid_list_get_named(plan, "control_in");
    SEXP ctrl_pl = cdid_list_get_named(plan, "control");
    if (control != ctrl_in && control != ctrl_pl) return Rf_ScalarLogical(0);

    SEXP norm = cdid_list_get_named(plan, "norm");
    SEXP fo = cdid_list_get_named(plan, "frame_owner");
    SEXP po = cdid_list_get_named(plan, "panel_owner");
    SEXP L_s = cdid_list_get_named(plan, "n_links");
    int L = 0, norm_L = 0;
    if (TYPEOF(norm) != VECSXP || TYPEOF(fo) != VECSXP || TYPEOF(po) != VECSXP) return Rf_ScalarLogical(0);
    if (!cdid_as_scalar_int(L_s, &L) || L < 1) return Rf_ScalarLogical(0);
    if (!cdid_as_scalar_int(cdid_list_get_named(norm, "n_links"), &norm_L) || norm_L != L) return Rf_ScalarLogical(0);
    if (cdid_list_get_named(norm, "frame_owner") != fo || cdid_list_get_named(norm, "panel_owner") != po) {
      return Rf_ScalarLogical(0);
    }
    SEXP plan_lids = cdid_list_get_named(plan, "link_ids");
    SEXP norm_lids = cdid_list_get_named(norm, "link_ids");
    if (!cdid_strsxp_eq(norm_lids, plan_lids) || XLENGTH(plan_lids) != L) return Rf_ScalarLogical(0);

    SEXP pls = cdid_list_get_named(plan, "point_links");
    SEXP lvs = cdid_list_get_named(plan, "link_views");
    SEXP nvs = cdid_list_get_named(plan, "nuisance_views");
    SEXP frs = cdid_list_get_named(norm, "frozen");
    SEXP nlinks = cdid_list_get_named(norm, "links");
    SEXP obs_refs = cdid_list_get_named(plan, "obs_a_idx_refs");
    SEXP psi_refs = cdid_list_get_named(plan, "psi_obs_refs");
    SEXP br_refs = cdid_list_get_named(plan, "b_rate_refs");
    if (TYPEOF(pls) != VECSXP || XLENGTH(pls) != L ||
        TYPEOF(lvs) != VECSXP || XLENGTH(lvs) != L ||
        TYPEOF(nvs) != VECSXP || XLENGTH(nvs) != L ||
        TYPEOF(frs) != VECSXP || XLENGTH(frs) != L ||
        TYPEOF(nlinks) != VECSXP || XLENGTH(nlinks) != L ||
        TYPEOF(obs_refs) != VECSXP || XLENGTH(obs_refs) != L ||
        TYPEOF(psi_refs) != VECSXP || XLENGTH(psi_refs) != L ||
        TYPEOF(br_refs) != VECSXP || XLENGTH(br_refs) != L) {
      return Rf_ScalarLogical(0);
    }

    SEXP sym_fo = Rf_install("frame_owner");
    SEXP sym_po = Rf_install("panel_owner");
    SEXP sym_cil = Rf_install("custom_ipw_link");
    SEXP sym_cif = Rf_install("custom_ipw_if_frozen");
    SEXP sym_lv = Rf_install("link_view");
    SEXP sym_nv = Rf_install("nuisance_view");

    if (!cdid_c_seal_ok(fo, sym_fo) || !cdid_c_seal_ok(po, sym_po)) return Rf_ScalarLogical(0);
    SEXP fo_ref = cdid_list_get_named(fo, "ref");
    SEXP fo_pref = cdid_list_get_named(fo, "panel_ref");
    SEXP po_ref = cdid_list_get_named(po, "ref");
    if (fo_pref != po_ref) return Rf_ScalarLogical(0);

    for (int ell = 0; ell < L; ++ell) {
      SEXP pl = VECTOR_ELT(pls, ell);
      SEXP lv = VECTOR_ELT(lvs, ell);
      SEXP nv = VECTOR_ELT(nvs, ell);
      SEXP fr = VECTOR_ELT(frs, ell);
      if (!cdid_c_seal_ok(pl, sym_cil) || !cdid_c_seal_ok(fr, sym_cif)) return Rf_ScalarLogical(0);
      SEXP pl_ref = cdid_list_get_named(pl, "ref");
      if (cdid_list_get_named(fr, "source_ref") != pl_ref) return Rf_ScalarLogical(0);
      if (VECTOR_ELT(nlinks, ell) != pl) return Rf_ScalarLogical(0);
      if (!cdid_c_seal_ok(lv, sym_lv)) return Rf_ScalarLogical(0);
      SEXP lv_ref = cdid_list_get_named(lv, "ref");
      SEXP lv_nref = cdid_list_get_named(lv, "nuisance_ref");
      if (cdid_list_get_named(pl, "link_view_ref") != lv_ref ||
          cdid_list_get_named(pl, "link_ref") != lv_ref ||
          cdid_list_get_named(pl, "frame_ref") != fo_ref ||
          cdid_list_get_named(pl, "panel_ref") != fo_pref ||
          cdid_list_get_named(lv, "frame_ref") != fo_ref ||
          cdid_list_get_named(lv, "panel_ref") != fo_pref ||
          cdid_list_get_named(pl, "nuisance_view_ref") != lv_nref ||
          cdid_list_get_named(pl, "nuisance_ref") != lv_nref) {
        return Rf_ScalarLogical(0);
      }
      if (nv == R_NilValue) {
        if (lv_nref != R_NilValue) return Rf_ScalarLogical(0);
      } else {
        if (!cdid_c_seal_ok(nv, sym_nv)) return Rf_ScalarLogical(0);
        if (lv_nref != cdid_list_get_named(nv, "ref") ||
            cdid_list_get_named(nv, "frame_ref") != fo_ref ||
            cdid_list_get_named(nv, "panel_ref") != fo_pref) {
          return Rf_ScalarLogical(0);
        }
      }
      if (cdid_list_get_named(pl, "obs_a_idx") != VECTOR_ELT(obs_refs, ell) ||
          cdid_list_get_named(pl, "psi_obs") != VECTOR_ELT(psi_refs, ell) ||
          cdid_list_get_named(pl, "b_rate") != VECTOR_ELT(br_refs, ell)) {
        return Rf_ScalarLogical(0);
      }
      SEXP pl_lid = cdid_list_get_named(pl, "link_id");
      if (pl_lid != R_NilValue && pl_lid != R_UnboundValue) {
        if (TYPEOF(pl_lid) != STRSXP || XLENGTH(pl_lid) != 1) return Rf_ScalarLogical(0);
        SEXP s = STRING_ELT(pl_lid, 0);
        if (s != NA_STRING && CHAR(s)[0] != '\0' && s != STRING_ELT(plan_lids, ell) &&
            strcmp(CHAR(s), CHAR(STRING_ELT(plan_lids, ell))) != 0) {
          return Rf_ScalarLogical(0);
        }
      }
    }
    return Rf_ScalarLogical(1);
  }

  if (mode == 2) {
    /* ipw_point_if_retained_check(carrier = a1, plan = a2, expected_state = a3) */
    SEXP carrier = a1;
    SEXP plan = a2;
    SEXP exp_st = a3;
    if (TYPEOF(carrier) != ENVSXP || cdid_parent_env(carrier) != R_EmptyEnv) return Rf_ScalarLogical(0);
    int bnd = 0, lkd = 0;
    SEXP c_kind = cdid_env_get_checked(carrier, Rf_install("kind"), &bnd, &lkd);
    if (!bnd || !cdid_scalar_str_eq(c_kind, "ipw_point_if_retained")) return Rf_ScalarLogical(0);
    SEXP c_state = cdid_env_get_checked(carrier, Rf_install("state"), &bnd, &lkd);
    if (!bnd || lkd || TYPEOF(c_state) != STRSXP || XLENGTH(c_state) != 1 || STRING_ELT(c_state, 0) == NA_STRING) {
      return Rf_ScalarLogical(0);
    }
    const char *st_str = CHAR(STRING_ELT(c_state, 0));
    int st_code = -1;
    if (strcmp(st_str, "reserved") == 0) st_code = 0;
    else if (strcmp(st_str, "filling") == 0) st_code = 1;
    else if (strcmp(st_str, "filled") == 0) st_code = 2;
    else if (strcmp(st_str, "slicing") == 0) st_code = 3;
    if (st_code < 0) return Rf_ScalarLogical(0);
    if (exp_st != R_NilValue) {
      if (TYPEOF(exp_st) != STRSXP || XLENGTH(exp_st) != 1 || STRING_ELT(exp_st, 0) == NA_STRING ||
          strcmp(st_str, CHAR(STRING_ELT(exp_st, 0))) != 0) {
        return Rf_ScalarLogical(0);
      }
    }

    SEXP seal = cdid_env_get_checked(carrier, Rf_install("seal"), &bnd, &lkd);
    if (!bnd || !lkd || TYPEOF(seal) != ENVSXP || cdid_parent_env(seal) != R_EmptyEnv ||
        !R_EnvironmentIsLocked(seal)) {
      return Rf_ScalarLogical(0);
    }

    SEXP snames = PROTECT(R_lsInternal3(seal, TRUE, FALSE));
    R_xlen_t nsn = XLENGTH(snames);
    for (R_xlen_t i = 0; i < nsn; ++i) {
      SEXP sym = Rf_installChar(STRING_ELT(snames, i));
      int cb = 0, cl = 0, sb = 0, sl = 0;
      SEXP cur = cdid_env_get_checked(carrier, sym, &cb, &cl);
      SEXP sld = cdid_env_get_checked(seal, sym, &sb, &sl);
      if (!cb || !cl || !sb || !sl || cur != sld) {
        UNPROTECT(1);
        return Rf_ScalarLogical(0);
      }
    }
    UNPROTECT(1);

    int phi_b = 0, phi_l = 0, raw_b = 0, raw_l = 0, pp_b = 0, pp_l = 0;
    SEXP phi = cdid_env_get_checked(carrier, Rf_install("phi"), &phi_b, &phi_l);
    SEXP raw_ref = cdid_env_get_checked(carrier, Rf_install("raw_ref"), &raw_b, &raw_l);
    SEXP pp = cdid_env_get_checked(carrier, Rf_install("producer_passes"), &pp_b, &pp_l);
    if (st_code <= 1) {
      int pp_v = -1;
      if (phi_b || raw_b || !pp_b || !cdid_as_scalar_int(pp, &pp_v) || pp_v != 0) return Rf_ScalarLogical(0);
    }

    if (TYPEOF(plan) != VECSXP) return Rf_ScalarLogical(0);
    int ps_b = 0, ps_l = 0;
    SEXP psnap = cdid_env_get_checked(carrier, Rf_install("plan_snapshot"), &ps_b, &ps_l);
    if (!ps_b || ps_l || psnap != plan) return Rf_ScalarLogical(0);

    int ver_b = 0, ver_l = 0;
    SEXP ver_s = cdid_env_get_checked(carrier, Rf_install("verifications"), &ver_b, &ver_l);
    if (!ver_b || ver_l || TYPEOF(ver_s) != INTSXP || XLENGTH(ver_s) != 1) return Rf_ScalarLogical(0);

    int cb1 = 0, cl1 = 0;
    SEXP c_pkind = cdid_env_get_checked(carrier, Rf_install("plan_kind"), &cb1, &cl1);
    SEXP c_pref = cdid_env_get_checked(carrier, Rf_install("plan_ref"), &cb1, &cl1);
    SEXP c_fref = cdid_env_get_checked(carrier, Rf_install("frame_ref"), &cb1, &cl1);
    SEXP c_panref = cdid_env_get_checked(carrier, Rf_install("panel_ref"), &cb1, &cl1);
    SEXP c_N = cdid_env_get_checked(carrier, Rf_install("N"), &cb1, &cl1);
    SEXP c_L = cdid_env_get_checked(carrier, Rf_install("L"), &cb1, &cl1);
    SEXP c_B = cdid_env_get_checked(carrier, Rf_install("B"), &cb1, &cl1);
    SEXP c_lids = cdid_env_get_checked(carrier, Rf_install("link_ids"), &cb1, &cl1);

    SEXP p_kind = cdid_list_get_named(plan, "kind");
    if (!cdid_strsxp_eq(p_kind, c_pkind) || XLENGTH(p_kind) != 1) return Rf_ScalarLogical(0);
    const char *pk_str = CHAR(STRING_ELT(p_kind, 0));
    int plan_L = 0;
    int same_refs = 0;
    if (strcmp(pk_str, "cdid_ipw_point_if_batch_plan") == 0) {
      SEXP batch = cdid_list_get_named(plan, "batch");
      SEXP fv = cdid_list_get_named(plan, "fitted_view");
      if (TYPEOF(batch) != VECSXP || TYPEOF(fv) != VECSXP) return Rf_ScalarLogical(0);
      if (!cdid_as_scalar_int(cdid_list_get_named(batch, "n_links"), &plan_L)) return Rf_ScalarLogical(0);
      SEXP c_fitref = cdid_env_get_checked(carrier, Rf_install("fitted_ref"), &cb1, &cl1);
      same_refs = (cdid_list_get_named(fv, "ref") == c_fitref);
    } else if (strcmp(pk_str, "cdid_ipw_point_if_general_plan") == 0) {
      if (!cdid_as_scalar_int(cdid_list_get_named(plan, "n_links"), &plan_L)) return Rf_ScalarLogical(0);
      SEXP c_frefs = cdid_env_get_checked(carrier, Rf_install("fitted_refs"), &cb1, &cl1);
      SEXP c_nrefs = cdid_env_get_checked(carrier, Rf_install("nuisance_refs"), &cb1, &cl1);
      same_refs = (cdid_list_get_named(plan, "fitted_refs") == c_frefs &&
                   cdid_list_get_named(plan, "nuisance_refs") == c_nrefs);
    } else if (strcmp(pk_str, "cdid_custom_ipw_if_plan") == 0) {
      if (!cdid_as_scalar_int(cdid_list_get_named(plan, "n_links"), &plan_L)) return Rf_ScalarLogical(0);
      same_refs = 1;
    } else {
      return Rf_ScalarLogical(0);
    }
    if (!same_refs) return Rf_ScalarLogical(0);

    SEXP p_fo = cdid_list_get_named(plan, "frame_owner");
    int N_v = 0, L_v = 0, B_v = 0, pN = 0, pB = 0;
    if (TYPEOF(p_fo) != VECSXP ||
        cdid_list_get_named(plan, "ref") != c_pref ||
        cdid_list_get_named(p_fo, "ref") != c_fref ||
        cdid_list_get_named(p_fo, "panel_ref") != c_panref ||
        !cdid_as_scalar_int(c_N, &N_v) || !cdid_as_scalar_int(c_L, &L_v) || !cdid_as_scalar_int(c_B, &B_v) ||
        !cdid_as_scalar_int(cdid_list_get_named(plan, "n_a_units"), &pN) ||
        !cdid_as_scalar_int(cdid_list_get_named(plan, "step"), &pB) ||
        N_v != pN || L_v != plan_L || B_v != pB ||
        cdid_list_get_named(plan, "link_ids") != c_lids) {
      return Rf_ScalarLogical(0);
    }

    if (st_code >= 2) {
      if (!phi_b || phi_l || !raw_b || raw_l || TYPEOF(raw_ref) != ENVSXP ||
          cdid_parent_env(raw_ref) != R_EmptyEnv || !R_EnvironmentIsLocked(raw_ref)) {
        return Rf_ScalarLogical(0);
      }
      int rpb = 0, rpl = 0;
      SEXP rphi = cdid_env_get_checked(raw_ref, Rf_install("phi"), &rpb, &rpl);
      if (!rpb || rpl || rphi != phi) return Rf_ScalarLogical(0);
      if (TYPEOF(phi) != REALSXP || !Rf_isMatrix(phi) || Rf_nrows(phi) != N_v || Rf_ncols(phi) != L_v) {
        return Rf_ScalarLogical(0);
      }
      SEXP dn = Rf_getAttrib(phi, R_DimNamesSymbol);
      if (TYPEOF(dn) != VECSXP || XLENGTH(dn) != 2 || !cdid_strsxp_eq(VECTOR_ELT(dn, 1), c_lids)) {
        return Rf_ScalarLogical(0);
      }
    }
    return Rf_ScalarLogical(1);
  }

  if (mode == 3) {
    /* cp04_retained_if_check(store = a1, backend = a2, links = a3, args_list = a4 = list(B, scale, check_slice)) */
    SEXP store = a1;
    SEXP backend = a2;
    SEXP links = a3;
    SEXP extra = a4;
    if (TYPEOF(store) != ENVSXP || TYPEOF(backend) != VECSXP || TYPEOF(extra) != VECSXP || XLENGTH(extra) < 3) {
      return R_NilValue;
    }
    int bnd = 0, lkd = 0;
    SEXP s_kind = cdid_env_get_checked(store, Rf_install("kind"), &bnd, &lkd);
    if (!bnd || !cdid_scalar_str_eq(s_kind, "cp04_retained_if")) return R_NilValue;
    SEXP s_state = cdid_env_get_checked(store, Rf_install("state"), &bnd, &lkd);
    if (!bnd || TYPEOF(s_state) != STRSXP || XLENGTH(s_state) != 1 || STRING_ELT(s_state, 0) == NA_STRING) {
      return R_NilValue;
    }
    const char *st_str = CHAR(STRING_ELT(s_state, 0));
    if (strcmp(st_str, "released") == 0) return R_NilValue;

    SEXP s_bsnap = cdid_env_get_checked(store, Rf_install("backend_snapshot"), &bnd, &lkd);
    if (!bnd || backend != s_bsnap) return R_NilValue;
    SEXP s_bref = cdid_env_get_checked(store, Rf_install("backend_ref"), &bnd, &lkd);
    if (!bnd || cdid_list_get_named(backend, "ref") != s_bref) return R_NilValue;

    /* Determine backend class and verify exact seal */
    const char *seal_tag = NULL;
    int is_raw = 0;
    if (Rf_inherits(backend, "cdid_ipw_raw_backend")) { seal_tag = "ipw_raw_backend"; is_raw = 1; }
    else if (Rf_inherits(backend, "cdid_dr_baseline")) { seal_tag = "dr_baseline_backend"; }
    else if (Rf_inherits(backend, "cdid_or_baseline")) { seal_tag = "or_baseline_backend"; }
    else if (Rf_inherits(backend, "cdid_cp04_corollary1")) { seal_tag = "cp04_backend"; }
    else if (Rf_inherits(backend, "cdid_or_corollary1")) { seal_tag = "or_corollary1_backend"; }
    else if (Rf_inherits(backend, "cdid_dr_corollary1")) { seal_tag = "dr_corollary1_backend"; }
    else if (Rf_inherits(backend, "cdid_imputation_backend")) { seal_tag = "imputation_backend"; }
    else if (Rf_inherits(backend, "cdid_augmented_backend")) { seal_tag = "augmented_backend"; }
    else if (Rf_inherits(backend, "cdid_cp04_custom_backend")) { seal_tag = "cp04_custom_backend"; }
    else if (Rf_inherits(backend, "cdid_reference_ipw_backend")) { seal_tag = "reference_ipw_backend"; }
    else return R_NilValue;

    if (!cdid_c_seal_ok(backend, Rf_install(seal_tag))) return R_NilValue;

    SEXP s_fsnap = cdid_env_get_checked(store, Rf_install("frame_snapshot"), &bnd, &lkd);
    SEXP s_psnap = cdid_env_get_checked(store, Rf_install("panel_snapshot"), &bnd, &lkd);
    SEXP s_fref = cdid_env_get_checked(store, Rf_install("frame_ref"), &bnd, &lkd);
    SEXP s_pref = cdid_env_get_checked(store, Rf_install("panel_ref"), &bnd, &lkd);
    SEXP s_rsnap = cdid_env_get_checked(store, Rf_install("reader_snapshot"), &bnd, &lkd);
    SEXP s_stsnap = cdid_env_get_checked(store, Rf_install("state_snapshot"), &bnd, &lkd);
    SEXP s_rref = cdid_env_get_checked(store, Rf_install("reader_ref"), &bnd, &lkd);
    if (s_fsnap == R_UnboundValue || s_psnap == R_UnboundValue ||
        s_fref == R_UnboundValue || s_pref == R_UnboundValue ||
        s_rsnap == R_UnboundValue || s_stsnap == R_UnboundValue || s_rref == R_UnboundValue) {
      return R_NilValue;
    }

    SEXP sym_fo = Rf_install("frame_owner");
    SEXP sym_po = Rf_install("panel_owner");
    int fb = 0, fl = 0, pb = 0, pl = 0;
    SEXP frame = cdid_env_get_checked(s_fref, sym_fo, &fb, &fl);
    SEXP panel = cdid_env_get_checked(s_pref, sym_po, &pb, &pl);
    if (!fb || !fl || frame != s_fsnap || !R_EnvironmentIsLocked(s_fref) || cdid_parent_env(s_fref) != R_EmptyEnv) return R_NilValue;
    if (!pb || !pl || panel != s_psnap || !R_EnvironmentIsLocked(s_pref) || cdid_parent_env(s_pref) != R_EmptyEnv) return R_NilValue;
    if (!cdid_c_seal_ok(frame, sym_fo) || !cdid_c_seal_ok(panel, sym_po)) return R_NilValue;
    if (cdid_list_get_named(frame, "panel_ref") != s_pref) return R_NilValue;

    if (is_raw) {
      if (s_rsnap != R_NilValue || s_stsnap != R_NilValue || s_rref != R_NilValue) return R_NilValue;
      if (cdid_list_get_named(backend, "frame_ref") != s_fref ||
          cdid_list_get_named(backend, "panel_ref") != s_pref) return R_NilValue;
      SEXP osnap = cdid_list_get_named(backend, "owner_snapshot");
      if (TYPEOF(osnap) != VECSXP ||
          cdid_list_get_named(osnap, "frame") != frame ||
          cdid_list_get_named(osnap, "panel") != panel) return R_NilValue;
    } else {
      SEXP b_rref = cdid_list_get_named(backend, "reader_ref");
      if (b_rref != s_rref || TYPEOF(b_rref) != ENVSXP || !R_EnvironmentIsLocked(b_rref) ||
          cdid_parent_env(b_rref) != R_EmptyEnv) return R_NilValue;
      SEXP sym_rd = Rf_install("cp04_reader");
      int rb = 0, rl = 0;
      SEXP reader = cdid_env_get_checked(b_rref, sym_rd, &rb, &rl);
      if (!rb || !rl || reader != s_rsnap || !cdid_c_seal_ok(reader, sym_rd)) return R_NilValue;
      if (cdid_list_get_named(reader, "state") != s_stsnap ||
          cdid_list_get_named(reader, "frame_ref") != s_fref ||
          cdid_list_get_named(reader, "panel_ref") != s_pref) return R_NilValue;
      /* Decline fast check if complete_support is attached so R checks complete_panel */
      SEXP cs = cdid_list_get_named(backend, "complete_support");
      if (cs != R_UnboundValue && cs != R_NilValue) return R_NilValue;
    }

    /* Optional attached references */
    SEXP b_refref = cdid_list_get_named(backend, "reference_ref");
    if (b_refref == R_UnboundValue) b_refref = R_NilValue;
    SEXP s_refref = cdid_env_get_checked(store, Rf_install("reference_ref"), &bnd, &lkd);
    if (!bnd || b_refref != s_refref) return R_NilValue;
    if (s_refref != R_NilValue) {
      SEXP s_refsnap = cdid_env_get_checked(store, Rf_install("reference_snapshot"), &bnd, &lkd);
      int rfb = 0, rfl = 0;
      SEXP ref_obj = cdid_env_get_checked(s_refref, Rf_install("reference_reader"), &rfb, &rfl);
      if (!bnd || !rfb || !rfl || ref_obj != s_refsnap || !cdid_c_seal_ok(ref_obj, Rf_install("reference_reader"))) {
        return R_NilValue;
      }
    }

    SEXP b_sel = cdid_list_get_named(backend, "selection_backend");
    if (b_sel == R_UnboundValue) b_sel = R_NilValue;
    SEXP s_selref = cdid_env_get_checked(store, Rf_install("selection_ref"), &bnd, &lkd);
    if (!bnd) return R_NilValue;
    if (b_sel == R_NilValue) {
      if (s_selref != R_NilValue) return R_NilValue;
    } else {
      SEXP s_selsnap = cdid_env_get_checked(store, Rf_install("selection_snapshot"), &bnd, &lkd);
      SEXP s_selrsnap = cdid_env_get_checked(store, Rf_install("selection_reader_snapshot"), &bnd, &lkd);
      if (cdid_list_get_named(b_sel, "ref") != s_selref || b_sel != s_selsnap ||
          !cdid_c_seal_ok(b_sel, Rf_install("selection_models"))) return R_NilValue;
      SEXP sr_ref = cdid_list_get_named(b_sel, "reader_ref");
      int srb = 0, srl = 0;
      SEXP sr_obj = cdid_env_get_checked(sr_ref, Rf_install("cp04_reader"), &srb, &srl);
      if (!srb || !srl || sr_obj != s_selrsnap) return R_NilValue;
    }

    SEXP b_dref = cdid_list_get_named(backend, "design_ref");
    if (b_dref == R_UnboundValue) b_dref = R_NilValue;
    SEXP s_dref = cdid_env_get_checked(store, Rf_install("design_ref"), &bnd, &lkd);
    if (!bnd || b_dref != s_dref) return R_NilValue;
    if (s_dref != R_NilValue) {
      SEXP s_dsnap = cdid_env_get_checked(store, Rf_install("design_snapshot"), &bnd, &lkd);
      int dsb = 0, dsl = 0;
      SEXP d_obj = cdid_env_get_checked(s_dref, Rf_install("imputation_design_owner"), &dsb, &dsl);
      if (!bnd || !dsb || !dsl || d_obj != s_dsnap || !cdid_c_seal_ok(d_obj, Rf_install("imputation_design_owner"))) {
        return R_NilValue;
      }
    }

    SEXP b_lay = cdid_list_get_named(backend, "layout");
    if (b_lay == R_UnboundValue) b_lay = R_NilValue;
    SEXP s_laysnap = cdid_env_get_checked(store, Rf_install("layout_snapshot"), &bnd, &lkd);
    if (!bnd || ((b_lay == R_NilValue) != (s_laysnap == R_NilValue))) return R_NilValue;
    if (s_laysnap != R_NilValue) {
      if (b_lay != s_laysnap) return R_NilValue;
      SEXP b_lst = cdid_list_get_named(b_lay, "state");
      if (b_lst == R_UnboundValue) b_lst = R_NilValue;
      SEXP s_lstsnap = cdid_env_get_checked(store, Rf_install("layout_state_snapshot"), &bnd, &lkd);
      if (!bnd || b_lst != s_lstsnap) return R_NilValue;
    }

    SEXP b_xf = cdid_list_get_named(backend, "X_F");
    if (b_xf == R_UnboundValue) b_xf = R_NilValue;
    SEXP s_xfsnap = cdid_env_get_checked(store, Rf_install("X_F_snapshot"), &bnd, &lkd);
    if (!bnd || b_xf != s_xfsnap) return R_NilValue;

    SEXP b_po = cdid_list_get_named(backend, "pair_outcomes");
    if (b_po == R_UnboundValue) b_po = R_NilValue;
    SEXP s_posnap = cdid_env_get_checked(store, Rf_install("pair_outcomes_snapshot"), &bnd, &lkd);
    if (!bnd || b_po != s_posnap) return R_NilValue;

    /* Partition, links, scale */
    int B_req = 0, s_B = 0;
    SEXP s_B_s = cdid_env_get_checked(store, Rf_install("B"), &bnd, &lkd);
    if (!bnd || !cdid_as_scalar_int(VECTOR_ELT(extra, 0), &B_req) || !cdid_as_scalar_int(s_B_s, &s_B)) {
      return R_NilValue;
    }
    int part_ok = (B_req == s_B);
    if (!part_ok && strcmp(st_str, "filled") == 0) {
      SEXP tr = cdid_env_get_checked(store, Rf_install("gmm_transfer"), &bnd, &lkd);
      if (bnd && TYPEOF(tr) == VECSXP) {
        int tr_B = 0;
        if (cdid_as_scalar_int(cdid_list_get_named(tr, "chunk_rows"), &tr_B) && B_req == tr_B) {
          part_ok = 1;
        }
      }
    }
    if (!part_ok) return R_NilValue;

    SEXP s_scale = cdid_env_get_checked(store, Rf_install("scale"), &bnd, &lkd);
    if (!bnd || !cdid_strsxp_eq(VECTOR_ELT(extra, 1), s_scale)) return R_NilValue;

    SEXP s_links = cdid_env_get_checked(store, Rf_install("links"), &bnd, &lkd);
    if (!bnd || TYPEOF(s_links) != INTSXP || Rf_isObject(links)) return R_NilValue;
    R_xlen_t nl = XLENGTH(s_links);
    if (XLENGTH(links) != nl) return R_NilValue;
    const int *slp = INTEGER(s_links);
    if (TYPEOF(links) == INTSXP) {
      const int *lp = INTEGER(links);
      if (lp != slp && memcmp(lp, slp, (size_t) nl * sizeof(int)) != 0) return R_NilValue;
    } else if (TYPEOF(links) == REALSXP) {
      const double *lp = REAL(links);
      for (R_xlen_t i = 0; i < nl; ++i) {
        if (!R_FINITE(lp[i]) || (int) lp[i] != slp[i] || lp[i] != (double) slp[i]) return R_NilValue;
      }
    } else {
      return R_NilValue;
    }

    int check_slice = 0;
    if (cdid_as_scalar_int(VECTOR_ELT(extra, 2), &check_slice) && check_slice != 0) {
      if (strcmp(st_str, "filled") != 0) return R_NilValue;
      int ph_b = 0, ph_l = 0, rw_b = 0, rw_l = 0;
      SEXP full = cdid_env_get_checked(store, Rf_install("phi"), &ph_b, &ph_l);
      SEXP raw_ref = cdid_env_get_checked(store, Rf_install("raw_ref"), &rw_b, &rw_l);
      if (!ph_b || ph_l || !rw_b || rw_l || TYPEOF(raw_ref) != ENVSXP ||
          cdid_parent_env(raw_ref) != R_EmptyEnv || !R_EnvironmentIsLocked(raw_ref)) {
        return R_NilValue;
      }
      int rpb = 0, rpl = 0;
      SEXP rphi = cdid_env_get_checked(raw_ref, Rf_install("phi"), &rpb, &rpl);
      int sN = 0, sL = 0;
      SEXP s_N = cdid_env_get_checked(store, Rf_install("N"), &bnd, &lkd);
      SEXP s_L = cdid_env_get_checked(store, Rf_install("L"), &bnd, &lkd);
      if (!rpb || rpl || rphi != full || !cdid_as_scalar_int(s_N, &sN) || !cdid_as_scalar_int(s_L, &sL) ||
          TYPEOF(full) != REALSXP || !Rf_isMatrix(full) || Rf_nrows(full) != sN || Rf_ncols(full) != sL) {
        return R_NilValue;
      }
    }

    SEXP s_id = cdid_env_get_checked(store, Rf_install("id"), &bnd, &lkd);
    SEXP s_N = cdid_env_get_checked(store, Rf_install("N"), &bnd, &lkd);
    SEXP s_nA = cdid_env_get_checked(store, Rf_install("n_A"), &bnd, &lkd);
    SEXP s_lam = cdid_env_get_checked(store, Rf_install("lambda"), &bnd, &lkd);
    if (s_id == R_UnboundValue || s_N == R_UnboundValue || s_nA == R_UnboundValue || s_lam == R_UnboundValue) {
      return R_NilValue;
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, Rf_mkChar("id"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("N"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("n_A"));
    SET_STRING_ELT(nms, 3, Rf_mkChar("lambda"));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    SET_VECTOR_ELT(out, 0, s_id);
    SET_VECTOR_ELT(out, 1, s_N);
    SET_VECTOR_ELT(out, 2, s_nA);
    SET_VECTOR_ELT(out, 3, s_lam);
    UNPROTECT(2);
    return out;
  }

  if (mode == 4) {
    /* Re-verify a previously-validated backend snapshot `rec = a2` against live `backend = a1` */
    SEXP backend = a1;
    SEXP rec = a2;
    if (TYPEOF(backend) != VECSXP || TYPEOF(rec) != VECSXP || XLENGTH(rec) < 6) return Rf_ScalarLogical(0);
    if (VECTOR_ELT(rec, 0) != backend) return Rf_ScalarLogical(0);
    SEXP seal_tag_s = VECTOR_ELT(rec, 1);
    if (TYPEOF(seal_tag_s) != STRSXP || XLENGTH(seal_tag_s) != 1) return Rf_ScalarLogical(0);
    if (!cdid_c_seal_ok(backend, Rf_installChar(STRING_ELT(seal_tag_s, 0)))) return Rf_ScalarLogical(0);

    SEXP exp_reader = VECTOR_ELT(rec, 2);
    SEXP exp_frame = VECTOR_ELT(rec, 3);
    SEXP exp_panel = VECTOR_ELT(rec, 4);
    SEXP exp_fref = VECTOR_ELT(rec, 5);

    SEXP sym_fo = Rf_install("frame_owner");
    SEXP sym_po = Rf_install("panel_owner");
    int fb = 0, fl = 0, pb = 0, pl = 0;

    if (exp_reader == R_NilValue) {
      /* Raw IPW backend */
      SEXP fref = cdid_list_get_named(backend, "frame_ref");
      SEXP pref = cdid_list_get_named(backend, "panel_ref");
      if (fref != exp_fref || TYPEOF(fref) != ENVSXP || TYPEOF(pref) != ENVSXP) return Rf_ScalarLogical(0);
      if (!R_EnvironmentIsLocked(fref) || !R_EnvironmentIsLocked(pref) ||
          cdid_parent_env(fref) != R_EmptyEnv || cdid_parent_env(pref) != R_EmptyEnv) return Rf_ScalarLogical(0);
      SEXP f_obj = cdid_env_get_checked(fref, sym_fo, &fb, &fl);
      SEXP p_obj = cdid_env_get_checked(pref, sym_po, &pb, &pl);
      if (!fb || !fl || f_obj != exp_frame || !pb || !pl || p_obj != exp_panel) return Rf_ScalarLogical(0);
      if (!cdid_c_seal_ok(f_obj, sym_fo) || !cdid_c_seal_ok(p_obj, sym_po)) return Rf_ScalarLogical(0);
      if (cdid_list_get_named(f_obj, "panel_ref") != pref) return Rf_ScalarLogical(0);
    } else {
      SEXP rref = cdid_list_get_named(backend, "reader_ref");
      if (TYPEOF(rref) != ENVSXP || !R_EnvironmentIsLocked(rref) || cdid_parent_env(rref) != R_EmptyEnv) {
        return Rf_ScalarLogical(0);
      }
      int rb = 0, rl = 0;
      SEXP sym_rd = Rf_install("cp04_reader");
      SEXP reader = cdid_env_get_checked(rref, sym_rd, &rb, &rl);
      if (!rb || !rl || reader != exp_reader || !cdid_c_seal_ok(reader, sym_rd)) return Rf_ScalarLogical(0);
      SEXP fref = cdid_list_get_named(reader, "frame_ref");
      SEXP pref = cdid_list_get_named(reader, "panel_ref");
      if (fref != exp_fref || TYPEOF(fref) != ENVSXP || TYPEOF(pref) != ENVSXP) return Rf_ScalarLogical(0);
      if (!R_EnvironmentIsLocked(fref) || !R_EnvironmentIsLocked(pref) ||
          cdid_parent_env(fref) != R_EmptyEnv || cdid_parent_env(pref) != R_EmptyEnv) return Rf_ScalarLogical(0);
      SEXP f_obj = cdid_env_get_checked(fref, sym_fo, &fb, &fl);
      SEXP p_obj = cdid_env_get_checked(pref, sym_po, &pb, &pl);
      if (!fb || !fl || f_obj != exp_frame || !pb || !pl || p_obj != exp_panel) return Rf_ScalarLogical(0);
      if (!cdid_c_seal_ok(f_obj, sym_fo) || !cdid_c_seal_ok(p_obj, sym_po)) return Rf_ScalarLogical(0);
      if (cdid_list_get_named(f_obj, "panel_ref") != pref) return Rf_ScalarLogical(0);
    }
    return Rf_ScalarLogical(1);
  }

  if (mode == 5) {
    /* Fast batch check of link_views list `a1` against an already-validated `frame_owner = a2` */
    SEXP lvs = a1;
    SEXP frame = a2;
    if (TYPEOF(lvs) != VECSXP || TYPEOF(frame) != VECSXP) return Rf_ScalarLogical(0);
    SEXP sym_fo = Rf_install("frame_owner");
    SEXP sym_lv = Rf_install("link_view");
    if (!cdid_c_seal_ok(frame, sym_fo)) return Rf_ScalarLogical(0);
    SEXP f_ref = cdid_list_get_named(frame, "ref");
    SEXP p_ref = cdid_list_get_named(frame, "panel_ref");
    R_xlen_t nl = XLENGTH(lvs);
    for (R_xlen_t i = 0; i < nl; ++i) {
      SEXP v = VECTOR_ELT(lvs, i);
      if (TYPEOF(v) != VECSXP || !Rf_inherits(v, "cdid_ipw_link_view") || !cdid_c_seal_ok(v, sym_lv)) {
        return Rf_ScalarLogical(0);
      }
      if (cdid_list_get_named(v, "frame_ref") != f_ref ||
          cdid_list_get_named(v, "panel_ref") != p_ref ||
          cdid_list_get_named(v, "nuisance_ref") != R_NilValue) {
        return Rf_ScalarLogical(0);
      }
    }
    return Rf_ScalarLogical(1);
  }

  if (mode == 6) {
    if (TYPEOF(a1) != STRSXP || TYPEOF(a2) != ENVSXP || TYPEOF(a3) != ENVSXP)
      return R_NilValue;
#if R_VERSION < R_Version(4, 6, 0)
    SEXP cache = PROTECT(cdid_legacy_binding_cache(a2));
    SEXP matches = cdid_default_fns_match(a1, a2, a3, cache);
    UNPROTECT(1);
    return matches;
#else
    return cdid_default_fns_match(a1, a2, a3, R_NilValue);
#endif
  }

  if (mode == 7) {
    /* Fast recursive `ipw_pipeline_scalar_state(x)` for non-data.frame objects. */
    int ok = 0;
    SEXP res = cdid_scalar_state_c(a1, &ok);
    return ok ? res : Rf_mkString("CDID_FALLBACK_DF");
  }

  if (mode == 8) {
    /* Fast C `ipw_owner_seal_ok(obj = a1, binding = a2)` */
    if (TYPEOF(a1) != VECSXP || TYPEOF(a2) != STRSXP || XLENGTH(a2) != 1 || STRING_ELT(a2, 0) == NA_STRING) {
      return Rf_ScalarLogical(0);
    }
    return Rf_ScalarLogical(cdid_c_seal_ok(a1, Rf_installChar(STRING_ELT(a2, 0))));
  }

  if (mode == 9) {
    /* Fast SEXP pointer-identity match for lists of environments/SEXPs:
     * a1 = VECSXP x, a2 = VECSXP table (or R_NilValue for self-match match(x, x)).
     * Returns 1-based INTSXP index of first identical SEXP pointer in table (0 if not found).
     */
    if (TYPEOF(a1) != VECSXP) return R_NilValue;
    R_xlen_t nx = XLENGTH(a1);
    SEXP tbl = (a2 == R_NilValue) ? a1 : a2;
    if (TYPEOF(tbl) != VECSXP) return R_NilValue;
    R_xlen_t nt = XLENGTH(tbl);
    SEXP out = PROTECT(Rf_allocVector(INTSXP, nx));
    int *pout = INTEGER(out);
    /* Both inputs are lists. Use the read-only accessor available before R 4.5. */
    const SEXP *px = DATAPTR_RO(a1);
    const SEXP *pt = DATAPTR_RO(tbl);
    for (R_xlen_t i = 0; i < nx; ++i) {
      SEXP xi = px[i];
      R_xlen_t lim = (a2 == R_NilValue) ? (i + 1) : nt;
      int found = 0;
      for (R_xlen_t j = 0; j < lim; ++j) {
        if (pt[j] == xi) {
          found = (int)(j + 1);
          break;
        }
      }
      pout[i] = found;
    }
    UNPROTECT(1);
    return out;
  }

  if (mode == 10) {
    /* Fast C intake for `chained_request_prepare_inputs` and direct covariate agreement:
     * Submode A (`a1 == R_NilValue`, `a3` = `list(data_covs, frame_covs, row_frame_pos)`):
     *   Checks each numeric covariate column `k`: no NA in `dv`, and `dv[i] == fv[row_frame_pos[i] - 1]`.
     *   Returns ScalarInteger(1) if complete & agree, ScalarInteger(-1) if NA in `dv`, ScalarInteger(0) if disagree,
     *   or R_NilValue if types are unsupported.
     * Submode B (`a1` = `list(id, tt, y, g)`, `a2` = `fid` or `R_NilValue`):
     *   Checks that `(id, tt)` is already sorted by strictly increasing numeric `id` across units and
     *   strictly increasing finite `tt` within each unit, with finite nonnegative `g` invariant within unit,
     *   and `y` having no NaN/Inf. Returns canonical intake vectors + `fid_eq`.
     */
    if (a1 == R_NilValue) {
      if (TYPEOF(a3) != VECSXP || XLENGTH(a3) < 3) return R_NilValue;
      SEXP dcovs = VECTOR_ELT(a3, 0);
      SEXP fcovs = VECTOR_ELT(a3, 1);
      SEXP rfp_s = VECTOR_ELT(a3, 2);
      if (TYPEOF(dcovs) != VECSXP || TYPEOF(fcovs) != VECSXP ||
          XLENGTH(dcovs) != XLENGTH(fcovs) || TYPEOF(rfp_s) != INTSXP) {
        return R_NilValue;
      }
      R_xlen_t K = XLENGTH(dcovs);
      R_xlen_t n_rows = XLENGTH(rfp_s);
      const int *rfp = INTEGER(rfp_s);
      int agree = 1;
      for (R_xlen_t k = 0; k < K; ++k) {
        SEXP dv = VECTOR_ELT(dcovs, k);
        SEXP fv = VECTOR_ELT(fcovs, k);
        if (Rf_isObject(dv) || Rf_isObject(fv) || XLENGTH(dv) != n_rows) return R_NilValue;
        R_xlen_t nf = XLENGTH(fv);
        if (TYPEOF(dv) == REALSXP && TYPEOF(fv) == REALSXP) {
          const double *pd = REAL(dv);
          const double *pf = REAL(fv);
          for (R_xlen_t i = 0; i < n_rows; ++i) {
            double d = pd[i];
            if (ISNAN(d)) return Rf_ScalarInteger(-1);
            int idx = rfp[i];
            if (idx < 1 || idx > nf) return R_NilValue;
            if (agree && d != pf[idx - 1]) agree = 0;
          }
        } else if (TYPEOF(dv) == INTSXP && TYPEOF(fv) == INTSXP) {
          const int *pd = INTEGER(dv);
          const int *pf = INTEGER(fv);
          for (R_xlen_t i = 0; i < n_rows; ++i) {
            int d = pd[i];
            if (d == NA_INTEGER) return Rf_ScalarInteger(-1);
            int idx = rfp[i];
            if (idx < 1 || idx > nf) return R_NilValue;
            if (agree && d != pf[idx - 1]) agree = 0;
          }
        } else {
          return R_NilValue;
        }
      }
      return Rf_ScalarInteger(agree);
    }

    if (TYPEOF(a1) != VECSXP || XLENGTH(a1) < 4) return R_NilValue;
    SEXP id_s = VECTOR_ELT(a1, 0);
    SEXP tt_s = VECTOR_ELT(a1, 1);
    SEXP y_s = VECTOR_ELT(a1, 2);
    SEXP g_s = VECTOR_ELT(a1, 3);
    if (Rf_isObject(id_s) || Rf_isObject(tt_s) || Rf_isObject(y_s) || Rf_isObject(g_s)) return R_NilValue;
    int id_t = TYPEOF(id_s), tt_t = TYPEOF(tt_s), y_t = TYPEOF(y_s), g_t = TYPEOF(g_s);
    if ((id_t != INTSXP && id_t != REALSXP) ||
        (tt_t != INTSXP && tt_t != REALSXP) ||
        (y_t != INTSXP && y_t != REALSXP) ||
        (g_t != INTSXP && g_t != REALSXP)) {
      return R_NilValue;
    }
    R_xlen_t n_rows = XLENGTH(id_s);
    if (n_rows < 1 || n_rows > INT_MAX ||
        XLENGTH(tt_s) != n_rows || XLENGTH(y_s) != n_rows || XLENGTH(g_s) != n_rows) {
      return R_NilValue;
    }
    const int *id_i = (id_t == INTSXP) ? INTEGER(id_s) : NULL;
    const double *id_d = (id_t == REALSXP) ? REAL(id_s) : NULL;
    const int *tt_i = (tt_t == INTSXP) ? INTEGER(tt_s) : NULL;
    const double *tt_d_in = (tt_t == REALSXP) ? REAL(tt_s) : NULL;
    const int *y_i = (y_t == INTSXP) ? INTEGER(y_s) : NULL;
    const double *y_d_in = (y_t == REALSXP) ? REAL(y_s) : NULL;
    const int *g_i = (g_t == INTSXP) ? INTEGER(g_s) : NULL;
    const double *g_d_in = (g_t == REALSXP) ? REAL(g_s) : NULL;
    int allow_inf_g = (a4 != R_NilValue && Rf_asInteger(a4) == 1);

    int n_units = 0;
    double cur_id = 0.0, cur_t = 0.0, cur_g = 0.0;
    double cal_buf[256];
    int n_cal = 0;
    int cal_overflow = 0;

    for (R_xlen_t i = 0; i < n_rows; ++i) {
      double u;
      if (id_i != NULL) {
        if (id_i[i] == NA_INTEGER) return R_NilValue;
        u = (double) id_i[i];
      } else {
        u = id_d[i];
        if (!R_FINITE(u)) return R_NilValue;
        if (u == 0.0) u = 0.0;
      }
      double t;
      if (tt_i != NULL) {
        if (tt_i[i] == NA_INTEGER) return R_NilValue;
        t = (double) tt_i[i];
      } else {
        t = tt_d_in[i];
        if (!R_FINITE(t)) return R_NilValue;
        if (t == 0.0) t = 0.0;
      }
      double gv;
      if (g_i != NULL) {
        if (g_i[i] == NA_INTEGER || (!allow_inf_g && g_i[i] < 0)) return R_NilValue;
        gv = (double) g_i[i];
      } else {
        gv = g_d_in[i];
        if (allow_inf_g) {
          if (ISNAN(gv) || gv == R_NegInf) return R_NilValue;
        } else {
          if (!R_FINITE(gv) || gv < 0.0) return R_NilValue;
        }
        if (gv == 0.0) gv = 0.0;
      }
      if (y_i != NULL) {
        /* INTSXP has only finite ints or NA_INTEGER */
      } else {
        double yv = y_d_in[i];
        if (!R_FINITE(yv) && !R_IsNA(yv)) return R_NilValue;
      }

      if (i == 0) {
        n_units = 1;
        cur_id = u;
        cur_t = t;
        cur_g = gv;
      } else if (u == cur_id) {
        if (!(t > cur_t) || gv != cur_g) return R_NilValue;
        cur_t = t;
      } else if (u > cur_id) {
        n_units++;
        cur_id = u;
        cur_t = t;
        cur_g = gv;
      } else {
        return R_NilValue;
      }

      if (!cal_overflow) {
        int lo = 0, hi = n_cal - 1, found = 0;
        while (lo <= hi) {
          int mid = (lo + hi) >> 1;
          double cm = cal_buf[mid];
          if (cm < t) lo = mid + 1;
          else if (cm > t) hi = mid - 1;
          else { found = 1; break; }
        }
        if (!found) {
          if (n_cal >= 256) {
            cal_overflow = 1;
          } else {
            for (int k = n_cal; k > lo; --k) cal_buf[k] = cal_buf[k - 1];
            cal_buf[lo] = t;
            n_cal++;
          }
        }
      }
    }
    if (cal_overflow) return R_NilValue;

    SEXP first_s = PROTECT(Rf_allocVector(INTSXP, n_units));
    SEXP unit_id_s = PROTECT(Rf_allocVector(id_t, n_units));
    SEXP unit_key_s = PROTECT(Rf_allocVector(REALSXP, n_units));
    SEXP pos_unit_s = PROTECT(Rf_allocVector(INTSXP, n_rows));
    SEXP g_unit_s = PROTECT(Rf_allocVector(g_t, n_units));
    SEXP g_unit_int_s = PROTECT(Rf_allocVector(REALSXP, n_units));
    SEXP n_obs_unit_s = PROTECT(Rf_allocVector(INTSXP, n_units));
    SEXP cal_s = PROTECT(Rf_allocVector(REALSXP, n_cal));
    SEXP tt_d_s = PROTECT(Rf_allocVector(REALSXP, n_rows));
    SEXP y_ord_s = PROTECT(Rf_allocVector(REALSXP, n_rows));
    SEXP g_ord_s = PROTECT(Rf_allocVector(REALSXP, n_rows));

    memcpy(REAL(cal_s), cal_buf, (size_t) n_cal * sizeof(double));
    int *p_first = INTEGER(first_s);
    int *p_uid_i = (id_t == INTSXP) ? INTEGER(unit_id_s) : NULL;
    double *p_uid_d = (id_t == REALSXP) ? REAL(unit_id_s) : NULL;
    double *p_ukey = REAL(unit_key_s);
    int *p_pos = INTEGER(pos_unit_s);
    int *p_gu_i = (g_t == INTSXP) ? INTEGER(g_unit_s) : NULL;
    double *p_gu_d = (g_t == REALSXP) ? REAL(g_unit_s) : NULL;
    double *p_gu_int = REAL(g_unit_int_s);
    int *p_nobs = INTEGER(n_obs_unit_s);
    double *p_tt_d = REAL(tt_d_s);
    double *p_y_ord = REAL(y_ord_s);
    double *p_g_ord = REAL(g_ord_s);
    memset(p_nobs, 0, (size_t) n_units * sizeof(int));

    int u_idx = -1;
    int n_recoded = 0;
    cur_id = 0.0;
    double cur_g_int = 0.0;
    for (R_xlen_t i = 0; i < n_rows; ++i) {
      double u = (id_i != NULL) ? (double) id_i[i] : id_d[i];
      if (u == 0.0) u = 0.0;
      double t = (tt_i != NULL) ? (double) tt_i[i] : tt_d_in[i];
      if (t == 0.0) t = 0.0;
      double yv;
      int y_obs;
      if (y_i != NULL) {
        if (y_i[i] == NA_INTEGER) { yv = NA_REAL; y_obs = 0; }
        else { yv = (double) y_i[i]; y_obs = 1; }
      } else {
        yv = y_d_in[i];
        y_obs = !ISNAN(yv);
      }
      if (i == 0 || u != cur_id) {
        u_idx++;
        cur_id = u;
        p_first[u_idx] = (int) (i + 1);
        if (p_uid_i != NULL) p_uid_i[u_idx] = id_i[i];
        else p_uid_d[u_idx] = id_d[i];
        p_ukey[u_idx] = u;
        double gv;
        if (p_gu_i != NULL) {
          p_gu_i[u_idx] = g_i[i];
          gv = (double) g_i[i];
        } else {
          p_gu_d[u_idx] = g_d_in[i];
          gv = g_d_in[i];
        }
        if (!allow_inf_g && gv == 0.0) {
          cur_g_int = R_PosInf;
          n_recoded++;
        } else {
          cur_g_int = (gv == 0.0) ? 0.0 : gv;
        }
        p_gu_int[u_idx] = cur_g_int;
      }
      p_pos[i] = u_idx + 1;
      if (y_obs) p_nobs[u_idx]++;
      p_tt_d[i] = t;
      p_y_ord[i] = yv;
      p_g_ord[i] = cur_g_int;
    }

    int fid_eq = 0;
    if (a2 != R_NilValue && !Rf_isObject(a2) && TYPEOF(a2) == id_t && XLENGTH(a2) == n_units) {
      fid_eq = 1;
      if (id_t == INTSXP) {
        const int *pf = INTEGER(a2);
        for (int u = 0; u < n_units; ++u) {
          if (pf[u] != p_uid_i[u]) { fid_eq = 0; break; }
        }
      } else {
        const double *pf = REAL(a2);
        for (int u = 0; u < n_units; ++u) {
          if (ISNAN(pf[u]) || pf[u] != p_uid_d[u]) { fid_eq = 0; break; }
        }
      }
    }

    const char *nms[] = {
      "first", "unit_id", "unit_key", "pos_unit", "g_unit", "g_unit_int",
      "n_obs_unit", "calendar", "tt_d", "y_ord", "g_ord", "n_recoded", "fid_eq", ""
    };
    SEXP out = PROTECT(Rf_mkNamed(VECSXP, nms));
    SET_VECTOR_ELT(out, 0, first_s);
    SET_VECTOR_ELT(out, 1, unit_id_s);
    SET_VECTOR_ELT(out, 2, unit_key_s);
    SET_VECTOR_ELT(out, 3, pos_unit_s);
    SET_VECTOR_ELT(out, 4, g_unit_s);
    SET_VECTOR_ELT(out, 5, g_unit_int_s);
    SET_VECTOR_ELT(out, 6, n_obs_unit_s);
    SET_VECTOR_ELT(out, 7, cal_s);
    SET_VECTOR_ELT(out, 8, tt_d_s);
    SET_VECTOR_ELT(out, 9, y_ord_s);
    SET_VECTOR_ELT(out, 10, g_ord_s);
    SET_VECTOR_ELT(out, 11, Rf_ScalarInteger(n_recoded));
    SET_VECTOR_ELT(out, 12, Rf_ScalarLogical(fid_eq));
    UNPROTECT(12);
    return out;
  }

  if (mode == 11) {
    /* Fast C assembly for `ipw_panel_owner_assemble`:
     * a1 = ui (INTSXP, length n, sorted 1..nu)
     * a2 = time_ord (REALSXP, length n, strictly increasing within unit)
     * a3 = y_ord (REALSXP, length n)
     * a4 = list(g_ord, nu_int, all_times)
     */
    if (TYPEOF(a1) != INTSXP || TYPEOF(a2) != REALSXP || TYPEOF(a3) != REALSXP ||
        TYPEOF(a4) != VECSXP || XLENGTH(a4) < 3) {
      return R_NilValue;
    }
    R_xlen_t n = XLENGTH(a1);
    if (n < 1 || n > INT_MAX || XLENGTH(a2) != n || XLENGTH(a3) != n) return R_NilValue;
    SEXP g_ord_s = VECTOR_ELT(a4, 0);
    int nu = 0;
    SEXP all_times_s = VECTOR_ELT(a4, 2);
    if (TYPEOF(g_ord_s) != REALSXP || XLENGTH(g_ord_s) != n ||
        !cdid_as_scalar_int(VECTOR_ELT(a4, 1), &nu) || nu < 1 ||
        TYPEOF(all_times_s) != REALSXP) {
      return R_NilValue;
    }
    int n_cal = (int) XLENGTH(all_times_s);
    if (n_cal < 1) return R_NilValue;

    const int *ui = INTEGER(a1);
    const double *time_ord = REAL(a2);
    const double *y_ord = REAL(a3);
    const double *g_ord = REAL(g_ord_s);
    const double *all_times = REAL(all_times_s);

    SEXP nr_s = PROTECT(Rf_allocVector(INTSXP, nu));
    SEXP no_s = PROTECT(Rf_allocVector(INTSXP, nu));
    SEXP n_na_u_s = PROTECT(Rf_allocVector(INTSXP, nu));
    SEXP ug_s = PROTECT(Rf_allocVector(REALSXP, nu));
    int *nr = INTEGER(nr_s);
    int *no = INTEGER(no_s);
    int *n_na_u = INTEGER(n_na_u_s);
    double *ug = REAL(ug_s);
    memset(nr, 0, (size_t) nu * sizeof(int));
    memset(no, 0, (size_t) nu * sizeof(int));
    for (int u = 0; u < nu; ++u) ug[u] = NA_REAL;

    int n_obs = 0, n_na = 0;
    int prev_u = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
      int u = ui[i];
      if (u < 1 || u > nu || u < prev_u) {
        UNPROTECT(4);
        return R_NilValue;
      }
      if (u != prev_u) {
        ug[u - 1] = g_ord[i];
        prev_u = u;
      }
      nr[u - 1]++;
      if (!ISNAN(y_ord[i])) {
        no[u - 1]++;
        n_obs++;
      } else {
        n_na++;
      }
    }

    SEXP by_unit_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP obs_start_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP obs_end_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP na_start_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP na_end_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP is_unobs_s = PROTECT(Rf_allocVector(LGLSXP, nu));
    SEXP is_single_s = PROTECT(Rf_allocVector(LGLSXP, nu));
    SEXP pair_elig_s = PROTECT(Rf_allocVector(LGLSXP, nu));
    SEXP unit_to_a_s = PROTECT(Rf_allocVector(INTSXP, nu));

    double *by_unit = REAL(by_unit_s);
    double *obs_start = REAL(obs_start_s);
    double *obs_end = REAL(obs_end_s);
    double *na_start = REAL(na_start_s);
    double *na_end = REAL(na_end_s);
    int *is_unobs = LOGICAL(is_unobs_s);
    int *is_single = LOGICAL(is_single_s);
    int *pair_elig = LOGICAL(pair_elig_s);
    int *unit_to_a = INTEGER(unit_to_a_s);

    double total_pairs = 0.0, cum_obs = 0.0, cum_na = 0.0;
    int n_single = 0, n_unobs = 0, n_a = 0;
    for (int u = 0; u < nu; ++u) {
      int no_u = no[u];
      int na_u_cnt = nr[u] - no_u;
      n_na_u[u] = na_u_cnt;
      double pairs_u = (no_u >= 2) ? ((double) no_u * (double) (no_u - 1) * 0.5) : 0.0;
      by_unit[u] = pairs_u;
      total_pairs += pairs_u;
      cum_obs += (double) no_u;
      obs_end[u] = cum_obs;
      obs_start[u] = cum_obs - (double) no_u + 1.0;
      cum_na += (double) na_u_cnt;
      na_end[u] = cum_na;
      na_start[u] = cum_na - (double) na_u_cnt + 1.0;
      is_unobs[u] = (no_u == 0);
      is_single[u] = (no_u == 1);
      pair_elig[u] = (no_u >= 2);
      if (no_u == 0) {
        n_unobs++;
        unit_to_a[u] = 0;
      } else if (no_u == 1) {
        n_single++;
        unit_to_a[u] = 0;
      } else {
        n_a++;
        unit_to_a[u] = n_a;
      }
    }

    SEXP single_s = PROTECT(Rf_allocVector(INTSXP, n_single));
    SEXP unobs_s = PROTECT(Rf_allocVector(INTSXP, n_unobs));
    SEXP a_unit_index_s = PROTECT(Rf_allocVector(INTSXP, n_a));
    int *single = INTEGER(single_s);
    int *unobs = INTEGER(unobs_s);
    int *a_unit_index = INTEGER(a_unit_index_s);
    int is_idx = 0, iu_idx = 0, ia_idx = 0;
    for (int u = 0; u < nu; ++u) {
      int no_u = no[u];
      if (no_u == 0) unobs[iu_idx++] = u + 1;
      else if (no_u == 1) single[is_idx++] = u + 1;
      else a_unit_index[ia_idx++] = u + 1;
    }

    SEXP ou_s = PROTECT(Rf_allocVector(INTSXP, n_obs));
    SEXP ot_s = PROTECT(Rf_allocVector(REALSXP, n_obs));
    SEXP oy_s = PROTECT(Rf_allocVector(REALSXP, n_obs));
    SEXP ocal_s = PROTECT(Rf_allocVector(INTSXP, n_obs));
    SEXP obs_a_index_s = PROTECT(Rf_allocVector(INTSXP, n_obs));
    SEXP na_u_s = PROTECT(Rf_allocVector(INTSXP, n_na));
    SEXP na_t_s = PROTECT(Rf_allocVector(REALSXP, n_na));
    SEXP na_cal_s = PROTECT(Rf_allocVector(INTSXP, n_na));
    SEXP date_count_s = PROTECT(Rf_allocVector(INTSXP, n_cal));

    int *ou = INTEGER(ou_s);
    double *ot = REAL(ot_s);
    double *oy = REAL(oy_s);
    int *ocal = INTEGER(ocal_s);
    int *obs_a_index = INTEGER(obs_a_index_s);
    int *na_u = INTEGER(na_u_s);
    double *na_t = REAL(na_t_s);
    int *na_cal = INTEGER(na_cal_s);
    int *date_count = INTEGER(date_count_s);
    memset(date_count, 0, (size_t) n_cal * sizeof(int));

    int io = 0, ina = 0, q_a = 0;
    prev_u = 0;
    int c_idx = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
      int u = ui[i];
      if (u != prev_u) {
        c_idx = 0;
        prev_u = u;
      }
      double t = time_ord[i];
      while (c_idx < n_cal && all_times[c_idx] < t) c_idx++;
      if (c_idx >= n_cal || all_times[c_idx] != t) {
        int lo = 0, hi = n_cal - 1, found = -1;
        while (lo <= hi) {
          int mid = (lo + hi) >> 1;
          if (all_times[mid] < t) lo = mid + 1;
          else if (all_times[mid] > t) hi = mid - 1;
          else { found = mid; break; }
        }
        if (found < 0) {
          UNPROTECT(25);
          return R_NilValue;
        }
        c_idx = found;
      }
      int cal_1 = c_idx + 1;
      double yv = y_ord[i];
      if (!ISNAN(yv)) {
        ou[io] = u;
        ot[io] = t;
        oy[io] = yv;
        ocal[io] = cal_1;
        int a_idx = unit_to_a[u - 1];
        obs_a_index[io] = a_idx;
        if (a_idx > 0) {
          date_count[c_idx]++;
          q_a++;
        }
        io++;
      } else {
        na_u[ina] = u;
        na_t[ina] = t;
        na_cal[ina] = cal_1;
        ina++;
      }
    }

    SEXP first_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP last_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP sy_s = PROTECT(Rf_allocVector(REALSXP, nu));
    SEXP stime_s = PROTECT(Rf_allocVector(REALSXP, nu));
    double *first = REAL(first_s);
    double *last = REAL(last_s);
    double *sy = REAL(sy_s);
    double *stime = REAL(stime_s);
    for (int u = 0; u < nu; ++u) {
      int no_u = no[u];
      if (no_u > 0) {
        int s_idx = (int) obs_start[u] - 1;
        int e_idx = (int) obs_end[u] - 1;
        first[u] = ot[s_idx];
        last[u] = ot[e_idx];
        if (no_u == 1) {
          sy[u] = oy[s_idx];
          stime[u] = first[u];
        } else {
          sy[u] = NA_REAL;
          stime[u] = NA_REAL;
        }
      } else {
        first[u] = NA_REAL;
        last[u] = NA_REAL;
        sy[u] = NA_REAL;
        stime[u] = NA_REAL;
      }
    }

    SEXP date_start_s = PROTECT(Rf_allocVector(INTSXP, n_cal));
    SEXP date_end_s = PROTECT(Rf_allocVector(INTSXP, n_cal));
    SEXP post_a_s = PROTECT(Rf_allocVector(INTSXP, q_a));
    SEXP post_u_s = PROTECT(Rf_allocVector(INTSXP, q_a));
    SEXP post_o_s = PROTECT(Rf_allocVector(INTSXP, q_a));
    int *date_start = INTEGER(date_start_s);
    int *date_end = INTEGER(date_end_s);
    int *post_a = INTEGER(post_a_s);
    int *post_u = INTEGER(post_u_s);
    int *post_o = INTEGER(post_o_s);

    const void *vmax = vmaxget();
    int *cursor = (int *) R_alloc((size_t) n_cal, sizeof(int));
    int cum_d = 0;
    for (int c = 0; c < n_cal; ++c) {
      cursor[c] = cum_d;
      date_start[c] = cum_d + 1;
      cum_d += date_count[c];
      date_end[c] = cum_d;
    }
    for (int k = 0; k < n_obs; ++k) {
      int a_idx = obs_a_index[k];
      if (a_idx > 0) {
        int pos = cursor[ocal[k] - 1]++;
        post_a[pos] = a_idx;
        post_u[pos] = ou[k];
        post_o[pos] = k + 1;
      }
    }
    vmaxset(vmax);

    const char *nms[] = {
      "nr", "no", "n_na_u", "ug", "by_unit", "total_pairs",
      "obs_start", "obs_end", "na_start", "na_end",
      "is_unobserved", "is_singleton", "pair_eligible", "unit_to_a",
      "single", "unobs", "a_unit_index",
      "ou", "ot", "oy", "ocal", "obs_a_index",
      "na_u", "na_t", "na_cal",
      "first", "last", "sy", "singleton_time",
      "date_count", "date_start", "date_end",
      "post_a_index", "post_unit_index", "post_obs_index", "q_a", ""
    };
    SEXP out = PROTECT(Rf_mkNamed(VECSXP, nms));
    SET_VECTOR_ELT(out, 0, nr_s);
    SET_VECTOR_ELT(out, 1, no_s);
    SET_VECTOR_ELT(out, 2, n_na_u_s);
    SET_VECTOR_ELT(out, 3, ug_s);
    SET_VECTOR_ELT(out, 4, by_unit_s);
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(total_pairs));
    SET_VECTOR_ELT(out, 6, obs_start_s);
    SET_VECTOR_ELT(out, 7, obs_end_s);
    SET_VECTOR_ELT(out, 8, na_start_s);
    SET_VECTOR_ELT(out, 9, na_end_s);
    SET_VECTOR_ELT(out, 10, is_unobs_s);
    SET_VECTOR_ELT(out, 11, is_single_s);
    SET_VECTOR_ELT(out, 12, pair_elig_s);
    SET_VECTOR_ELT(out, 13, unit_to_a_s);
    SET_VECTOR_ELT(out, 14, single_s);
    SET_VECTOR_ELT(out, 15, unobs_s);
    SET_VECTOR_ELT(out, 16, a_unit_index_s);
    SET_VECTOR_ELT(out, 17, ou_s);
    SET_VECTOR_ELT(out, 18, ot_s);
    SET_VECTOR_ELT(out, 19, oy_s);
    SET_VECTOR_ELT(out, 20, ocal_s);
    SET_VECTOR_ELT(out, 21, obs_a_index_s);
    SET_VECTOR_ELT(out, 22, na_u_s);
    SET_VECTOR_ELT(out, 23, na_t_s);
    SET_VECTOR_ELT(out, 24, na_cal_s);
    SET_VECTOR_ELT(out, 25, first_s);
    SET_VECTOR_ELT(out, 26, last_s);
    SET_VECTOR_ELT(out, 27, sy_s);
    SET_VECTOR_ELT(out, 28, stime_s);
    SET_VECTOR_ELT(out, 29, date_count_s);
    SET_VECTOR_ELT(out, 30, date_start_s);
    SET_VECTOR_ELT(out, 31, date_end_s);
    SET_VECTOR_ELT(out, 32, post_a_s);
    SET_VECTOR_ELT(out, 33, post_u_s);
    SET_VECTOR_ELT(out, 34, post_o_s);
    SET_VECTOR_ELT(out, 35, Rf_ScalarInteger(q_a));
    UNPROTECT(35);
    return out;
  }

  return R_NilValue;
}

/* Check all entries of `src_mat` finite AND copy `B_r x L` into `dst_mat[row_start .. row_start+B_r-1, ]`
   in place (caller owns `dst_mat` as a fresh local matrix before sealing).
   Also supports multi-column signed scatter when `row_start_sexp` is `list(rows_or_NULL, col_dst, col_src, signs_or_NULL)`. */
SEXP C_cdid_matrix_set_rows(SEXP dst_mat, SEXP row_start_sexp, SEXP src_mat) {
  if (TYPEOF(dst_mat) != REALSXP || !Rf_isMatrix(dst_mat) ||
      TYPEOF(src_mat) != REALSXP || !Rf_isMatrix(src_mat)) {
    return Rf_ScalarLogical(0);
  }
  int N = Rf_nrows(dst_mat);
  int L = Rf_ncols(dst_mat);
  int B_r = Rf_nrows(src_mat);
  int L_src = Rf_ncols(src_mat);
  if (TYPEOF(row_start_sexp) == VECSXP) {
    if (XLENGTH(row_start_sexp) < 3 || N < 1 || B_r < 1) return Rf_ScalarLogical(0);
    SEXP r_sexp = VECTOR_ELT(row_start_sexp, 0);
    SEXP cd_sexp = VECTOR_ELT(row_start_sexp, 1);
    SEXP cs_sexp = VECTOR_ELT(row_start_sexp, 2);
    SEXP sgn_sexp = (XLENGTH(row_start_sexp) >= 4) ? VECTOR_ELT(row_start_sexp, 3) : R_NilValue;
    if (TYPEOF(cd_sexp) != INTSXP || TYPEOF(cs_sexp) != INTSXP || XLENGTH(cd_sexp) != XLENGTH(cs_sexp)) {
      return Rf_ScalarLogical(0);
    }
    R_xlen_t K = XLENGTH(cd_sexp);
    if (sgn_sexp != R_NilValue && (TYPEOF(sgn_sexp) != REALSXP || XLENGTH(sgn_sexp) != K)) {
      return Rf_ScalarLogical(0);
    }
    const int *r_ptr = NULL;
    int contig_r0 = -1;
    if (r_sexp == R_NilValue) {
      if (B_r != N) return Rf_ScalarLogical(0);
      contig_r0 = 0;
    } else if (TYPEOF(r_sexp) == INTSXP && XLENGTH(r_sexp) == B_r) {
      r_ptr = INTEGER(r_sexp);
      int is_contig = (r_ptr[0] >= 1 && r_ptr[0] - 1 + B_r <= N);
      for (int i = 0; i < B_r; ++i) {
        if (r_ptr[i] < 1 || r_ptr[i] > N) return Rf_ScalarLogical(0);
        if (is_contig && r_ptr[i] != r_ptr[0] + i) is_contig = 0;
      }
      if (is_contig) contig_r0 = r_ptr[0] - 1;
    } else {
      return Rf_ScalarLogical(0);
    }
    const int *cd = INTEGER(cd_sexp);
    const int *cs = INTEGER(cs_sexp);
    const double *sgn = (sgn_sexp != R_NilValue) ? REAL(sgn_sexp) : NULL;
    for (R_xlen_t k = 0; k < K; ++k) {
      if (cd[k] < 1 || cd[k] > L || cs[k] < 1 || cs[k] > L_src) return Rf_ScalarLogical(0);
    }
    double *dst = REAL(dst_mat);
    const double *src = REAL(src_mat);
    for (R_xlen_t k = 0; k < K; ++k) {
      double *dcol = dst + (size_t) (cd[k] - 1) * (size_t) N;
      const double *scol = src + (size_t) (cs[k] - 1) * (size_t) B_r;
      double sk = (sgn != NULL) ? sgn[k] : 1.0;
      if (contig_r0 >= 0) {
        double *dptr = dcol + (size_t) contig_r0;
        if (sk == 1.0) {
          memcpy(dptr, scol, (size_t) B_r * sizeof(double));
        } else if (sk == -1.0) {
          for (int i = 0; i < B_r; ++i) dptr[i] = -scol[i];
        } else {
          for (int i = 0; i < B_r; ++i) dptr[i] = sk * scol[i];
        }
      } else {
        for (int i = 0; i < B_r; ++i) dcol[r_ptr[i] - 1] = sk * scol[i];
      }
    }
    return Rf_ScalarLogical(1);
  }
  int s1 = 0;
  if (!cdid_as_scalar_int(row_start_sexp, &s1)) {
    return Rf_ScalarLogical(0);
  }
  if (L_src != L || s1 < 1 || B_r < 1 || s1 - 1 + B_r > N) {
    return Rf_ScalarLogical(0);
  }
  const double *src = REAL(src_mat);
  R_xlen_t tot = (R_xlen_t) B_r * (R_xlen_t) L;
  for (R_xlen_t i = 0; i < tot; ++i) {
    if (!R_FINITE(src[i])) return Rf_ScalarLogical(0);
  }
  double *dst = REAL(dst_mat);
  size_t r0 = (size_t) (s1 - 1);
  for (int c = 0; c < L; ++c) {
    memcpy(dst + (size_t) c * (size_t) N + r0,
           src + (size_t) c * (size_t) B_r,
           (size_t) B_r * sizeof(double));
  }
  return Rf_ScalarLogical(1);
}

/* Single-call centered Gram matrix over retained `N x L` `phi_full`.
 * When `mode == 1`, `means_sexp` may be `R_NilValue`; `means` is computed in C
 * using the exact chunked `colSums(phi[s0 : min(N, s0+step-1), ]) / N` loop order
 * and `list(means = means_s, V = V_s)` is returned.
 */
SEXP C_cdid_retained_cov_gram(SEXP phi_full, SEXP means_sexp, SEXP step_sexp, SEXP denom_sexp, SEXP mode_sexp) {
  int step = 0, mode = 0;
  if (TYPEOF(phi_full) != REALSXP || !Rf_isMatrix(phi_full) ||
      !cdid_as_scalar_int(step_sexp, &step) || step < 1 ||
      TYPEOF(denom_sexp) != REALSXP || XLENGTH(denom_sexp) != 1 ||
      !cdid_as_scalar_int(mode_sexp, &mode)) {
    return R_NilValue;
  }
  int N = Rf_nrows(phi_full);
  int L = Rf_ncols(phi_full);
  if (N < 1 || L < 1) return R_NilValue;
  double denom = REAL(denom_sexp)[0];
  if (!R_FINITE(denom) || denom <= 0.0) return R_NilValue;
  const double *phi = REAL(phi_full);
  R_xlen_t tot_phi = (R_xlen_t) N * (R_xlen_t) L;
  for (R_xlen_t idx = 0; idx < tot_phi; ++idx) {
    if (!R_FINITE(phi[idx])) return R_NilValue;
  }

  int nprot = 0;
  SEXP means_s = R_NilValue;
  const double *means = NULL;
  if (mode == 2) {
    means_s = PROTECT(Rf_allocVector(REALSXP, L));
    nprot++;
    double *pmeans = REAL(means_s);
    double dN = (double) N;
    for (int c = 0; c < L; ++c) {
      const double *pcol = phi + (size_t) c * (size_t) N;
      double m_acc = 0.0;
      for (int s0 = 0; s0 < N; s0 += step) {
        int rb = (s0 + step <= N) ? step : (N - s0);
        long double cs = 0.0;
        for (int i = 0; i < rb; ++i) cs += pcol[s0 + i];
        m_acc += ((double) cs) / dN;
      }
      if (!R_FINITE(m_acc)) {
        UNPROTECT(nprot);
        return R_NilValue;
      }
      pmeans[c] = m_acc;
    }
    means = pmeans;
  } else {
    if (TYPEOF(means_sexp) != REALSXP || XLENGTH(means_sexp) != L) return R_NilValue;
    means = REAL(means_sexp);
    for (int c = 0; c < L; ++c) {
      if (!R_FINITE(means[c])) return R_NilValue;
    }
  }
  SEXP V_s = PROTECT(Rf_allocMatrix(REALSXP, L, L));
  nprot++;
  double *V = REAL(V_s);
  memset(V, 0, (size_t) L * (size_t) L * sizeof(double));

  const void *vmax = vmaxget();
  int max_rb = (N < step) ? N : step;
  double *xt = (double *) R_alloc((size_t) max_rb * (size_t) L, sizeof(double));
  double *chunk_g = (double *) R_alloc((size_t) L * (size_t) L, sizeof(double));
  int L_b = L;
  double one_d = 1.0, zero_d = 0.0;

  for (int r0 = 0; r0 < N; r0 += step) {
    int rb = (r0 + step <= N) ? step : (N - r0);
    for (int c = 0; c < L; ++c) {
      const double *src = phi + (size_t) c * (size_t) N + (size_t) r0;
      double *dst = xt + (size_t) c * (size_t) rb;
      double mc = means[c];
      for (int i = 0; i < rb; ++i) dst[i] = src[i] - mc;
    }
    int rb_b = rb;
    memset(chunk_g, 0, (size_t) L * (size_t) L * sizeof(double));
    F77_CALL(dsyrk)("U", "T", &L_b, &rb_b, &one_d, xt, &rb_b, &zero_d, chunk_g, &L_b FCONE FCONE);
    if (mode == 1) {
      for (int c2 = 0; c2 < L; ++c2) {
        double *vcol = V + (size_t) c2 * (size_t) L;
        const double *gcol = chunk_g + (size_t) c2 * (size_t) L;
        for (int c1 = 0; c1 <= c2; ++c1) vcol[c1] += gcol[c1];
      }
    } else {
      for (int c2 = 0; c2 < L; ++c2) {
        double *vcol = V + (size_t) c2 * (size_t) L;
        const double *gcol = chunk_g + (size_t) c2 * (size_t) L;
        for (int c1 = 0; c1 <= c2; ++c1) vcol[c1] += gcol[c1] / denom;
      }
    }
  }

  for (int c2 = 0; c2 < L; ++c2) {
    for (int c1 = 0; c1 <= c2; ++c1) {
      double v = V[c1 + (size_t) c2 * (size_t) L];
      if (mode == 1) v /= denom;
      if (!R_FINITE(v)) {
        vmaxset(vmax);
        UNPROTECT(nprot);
        return R_NilValue;
      }
      V[c1 + (size_t) c2 * (size_t) L] = v;
      V[c2 + (size_t) c1 * (size_t) L] = v;
    }
  }

  vmaxset(vmax);
  if (mode == 2) {
    const char *nms[] = {"means", "V", ""};
    SEXP out = PROTECT(Rf_mkNamed(VECSXP, nms));
    SET_VECTOR_ELT(out, 0, means_s);
    SET_VECTOR_ELT(out, 1, V_s);
    UNPROTECT(nprot + 1);
    return out;
  }
  UNPROTECT(nprot);
  return V_s;
}

/* Batch Cholesky solve for `dr_sparse_native_if`:
   Computes `xi[[j]] = t(backsolve(R_j, forwardsolve(t(R_j), t(rhs[, (j-1)*pt + 1:pt]))))`
   for `j = 1 .. n_ok` directly on column-major `B x pt` blocks with zero transposes. */
SEXP C_cdid_dr_batch_chol_solve(SEXP trt_list, SEXP tids_sexp, SEXP rhs_mat, SEXP dims_sexp) {
  if (TYPEOF(trt_list) != VECSXP || TYPEOF(tids_sexp) != INTSXP ||
      TYPEOF(rhs_mat) != REALSXP || !Rf_isMatrix(rhs_mat) ||
      TYPEOF(dims_sexp) != INTSXP || XLENGTH(dims_sexp) < 2) {
    return R_NilValue;
  }
  int pt = INTEGER(dims_sexp)[0];
  int n_ok = INTEGER(dims_sexp)[1];
  if (pt < 0 || n_ok < 0) return R_NilValue;
  int B = Rf_nrows(rhs_mat);
  int rhs_cols = Rf_ncols(rhs_mat);
  if (n_ok > XLENGTH(tids_sexp) || (n_ok > 0 && (long long) n_ok * (long long) pt > rhs_cols)) {
    return R_NilValue;
  }

  SEXP xi_list = PROTECT(Rf_allocVector(VECSXP, n_ok));
  if (n_ok == 0) {
    UNPROTECT(1);
    return xi_list;
  }

  const int *tids = INTEGER(tids_sexp);
  const double *rhs = REAL(rhs_mat);
  R_xlen_t n_trts = XLENGTH(trt_list);

  for (int j = 0; j < n_ok; ++j) {
    int tid = tids[j];
    if (tid < 1 || tid > n_trts) {
      UNPROTECT(1);
      return R_NilValue;
    }
    SEXP trt = VECTOR_ELT(trt_list, tid - 1);
    SEXP R_s = cdid_list_get_named(trt, "H_chol");
    if (R_s == R_UnboundValue || TYPEOF(R_s) != REALSXP || !Rf_isMatrix(R_s) ||
        Rf_nrows(R_s) != pt || Rf_ncols(R_s) != pt) {
      UNPROTECT(1);
      return R_NilValue;
    }
    const double *R = REAL(R_s);
    for (int k = 0; k < pt; ++k) {
      double d = R[k + (size_t) k * (size_t) pt];
      if (!R_FINITE(d) || d == 0.0) {
        UNPROTECT(1);
        return R_NilValue;
      }
    }

    SEXP x_s = PROTECT(Rf_allocMatrix(REALSXP, B, pt));
    double *X = REAL(x_s);
    const double *Y = rhs + (size_t) j * (size_t) pt * (size_t) B;

    /* Step 1: Z = Y * R^{-1} (exact column-vectorized equivalent of forwardsolve(t(R), t(Y))) */
    for (int c = 0; c < pt; ++c) {
      double *z_c = X + (size_t) c * (size_t) B;
      const double *y_c = Y + (size_t) c * (size_t) B;
      const double *r_col = R + (size_t) c * (size_t) pt;
      memcpy(z_c, y_c, (size_t) B * sizeof(double));
      for (int k = 0; k < c; ++k) {
        double r_kc = r_col[k];
        const double *z_k = X + (size_t) k * (size_t) B;
        for (int i = 0; i < B; ++i) {
          z_c[i] -= z_k[i] * r_kc;
        }
      }
      double r_cc = r_col[c];
      for (int i = 0; i < B; ++i) {
        z_c[i] /= r_cc;
      }
    }

    /* Step 2: X = Z * R^{-T} (exact column-vectorized equivalent of backsolve(R, t(Z))) */
    for (int c = pt - 1; c >= 0; --c) {
      double *x_c = X + (size_t) c * (size_t) B;
      for (int k = c + 1; k < pt; ++k) {
        double r_ck = R[c + (size_t) k * (size_t) pt];
        const double *x_k = X + (size_t) k * (size_t) B;
        for (int i = 0; i < B; ++i) {
          x_c[i] -= x_k[i] * r_ck;
        }
      }
      double r_cc = R[c + (size_t) c * (size_t) pt];
      for (int i = 0; i < B; ++i) {
        x_c[i] /= r_cc;
      }
    }

    SET_VECTOR_ELT(xi_list, j, x_s);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return xi_list;
}

static SEXP cdid_make_diag_unavailable(int kish, const char *reason) {
  int n_elt = kish ? 3 : 4;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n_elt));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, n_elt));
  SET_STRING_ELT(nms, 0, Rf_mkChar("status"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("reason"));
  SET_VECTOR_ELT(out, 0, Rf_mkString("unavailable"));
  SET_VECTOR_ELT(out, 1, Rf_mkString(reason));
  if (kish) {
    SET_STRING_ELT(nms, 2, Rf_mkChar("value"));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(NA_REAL));
  } else {
    SET_STRING_ELT(nms, 2, Rf_mkChar("min"));
    SET_STRING_ELT(nms, 3, Rf_mkChar("max"));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(NA_REAL));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(NA_REAL));
  }
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}

/* Exact chunk-preserving C implementation of `ipw_fitted_diagnostic_summary`. */
SEXP C_cdid_diag_summary(SEXP v_sexp, SEXP chunk_rows_sexp, SEXP kish_sexp) {
  if (TYPEOF(v_sexp) != REALSXP || Rf_isObject(v_sexp) ||
      TYPEOF(kish_sexp) != LGLSXP || XLENGTH(kish_sexp) != 1 || LOGICAL(kish_sexp)[0] == NA_LOGICAL) {
    return R_NilValue;
  }
  int kish = (LOGICAL(kish_sexp)[0] != 0);
  double chunk_d = 0.0;
  if (TYPEOF(chunk_rows_sexp) == INTSXP && XLENGTH(chunk_rows_sexp) == 1 && INTEGER(chunk_rows_sexp)[0] != NA_INTEGER) {
    chunk_d = (double) INTEGER(chunk_rows_sexp)[0];
  } else if (TYPEOF(chunk_rows_sexp) == REALSXP && XLENGTH(chunk_rows_sexp) == 1) {
    chunk_d = REAL(chunk_rows_sexp)[0];
  } else {
    return R_NilValue;
  }
  if (!R_FINITE(chunk_d) || chunk_d < 1.0 || chunk_d != floor(chunk_d)) return R_NilValue;

  R_xlen_t n = XLENGTH(v_sexp);
  if (n == 0) return cdid_make_diag_unavailable(kish, "empty_vector");

  const double *v = REAL(v_sexp);
  R_xlen_t step = (R_xlen_t) chunk_d;
  double lo = R_PosInf, hi = R_NegInf;

  for (R_xlen_t start = 0; start < n; start += step) {
    R_xlen_t end = (n - start > step) ? (start + step) : n;
    for (R_xlen_t i = start; i < end; ++i) {
      if (!R_FINITE(v[i])) return cdid_make_diag_unavailable(kish, "nonfinite_vector");
    }
    double c_lo = v[start], c_hi = v[start];
    for (R_xlen_t i = start + 1; i < end; ++i) {
      if (v[i] < c_lo) c_lo = v[i];
      if (v[i] > c_hi) c_hi = v[i];
    }
    if (c_lo < lo) lo = c_lo;
    if (c_hi > hi) hi = c_hi;
  }

  if (!kish) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, Rf_mkChar("status"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("reason"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("min"));
    SET_STRING_ELT(nms, 3, Rf_mkChar("max"));
    SET_VECTOR_ELT(out, 0, Rf_mkString("available"));
    SET_VECTOR_ELT(out, 1, Rf_ScalarString(NA_STRING));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(lo));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(hi));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(2);
    return out;
  }

  if (lo <= 0.0) return cdid_make_diag_unavailable(1, "nonpositive_weights");

  double total = 0.0, squares = 0.0;
  for (R_xlen_t start = 0; start < n; start += step) {
    R_xlen_t end = (n - start > step) ? (start + step) : n;
    double s1 = 0.0, s2 = 0.0;
    for (R_xlen_t i = start; i < end; ++i) {
      double u = v[i] / hi;
      s1 += u;
      s2 += u * u;
    }
    total += s1;
    squares += s2;
  }
  double val = (total * total) / squares;
  if (!R_FINITE(val)) return cdid_make_diag_unavailable(1, "nonfinite_diagnostic_arithmetic");

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(nms, 0, Rf_mkChar("status"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("reason"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("value"));
  SET_VECTOR_ELT(out, 0, Rf_mkString("available"));
  SET_VECTOR_ELT(out, 1, Rf_ScalarString(NA_STRING));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(val));
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}


