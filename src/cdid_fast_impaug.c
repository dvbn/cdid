#pragma STDC FP_CONTRACT OFF

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "cdid_r_api.h"

#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>

static SEXP cdid_ia_get_elt(SEXP list, const char *name) {
  if (TYPEOF(list) != VECSXP) return R_NilValue;
  SEXP nms = Rf_getAttrib(list, R_NamesSymbol);
  if (TYPEOF(nms) != STRSXP) return R_NilValue;
  R_xlen_t n = XLENGTH(list);
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP s = STRING_ELT(nms, i);
    if (s != NA_STRING && strcmp(CHAR(s), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

static inline void cdid_ia_clip_logit(double e, double *prob_out, double *odds_out) {
  const double clip_cap = 1.0 - 1e-6;
  const double clip_cap_odds = clip_cap / (1.0 - clip_cap);
  double od, pr;
  if (e > 40.0) {
    pr = 1.0 / (1.0 + exp(-e));
    od = 1.0;
  } else {
    od = exp(e);
    pr = od / (1.0 + od);
  }
  *prob_out = pr;
  *odds_out = (pr < clip_cap) ? od : clip_cap_odds;
}

static inline int cdid_ia_find_posting(
    const int *post_a,
    int lo,
    int hi,
    int target,
    int *hint) {
  if (lo > hi) return -1;
  int h = *hint;
  if (h >= lo && h <= hi) {
    int vh = post_a[h];
    if (vh == target) return h;
    if (vh < target) {
      if (h + 1 <= hi && post_a[h + 1] == target) {
        *hint = h + 1;
        return h + 1;
      }
      lo = h + 1;
    } else {
      hi = h - 1;
    }
  }
  while (lo <= hi) {
    int mid = lo + ((hi - lo) >> 1);
    int vm = post_a[mid];
    if (vm == target) {
      *hint = mid;
      return mid;
    }
    if (vm < target) lo = mid + 1;
    else hi = mid - 1;
  }
  if (lo >= 0) *hint = lo;
  return -1;
}

static int cdid_ia_get_wls_fit(SEXP fit_s, const double **R_out, const int **piv_out, int p) {
  if (p <= 0) {
    *R_out = NULL;
    *piv_out = NULL;
    return 1;
  }
  if (TYPEOF(fit_s) != VECSXP) return 0;
  SEXP sf_s = cdid_ia_get_elt(fit_s, "slope_fit");
  if (sf_s != R_NilValue && TYPEOF(sf_s) == VECSXP) {
    fit_s = sf_s;
  }
  SEXP R_s = cdid_ia_get_elt(fit_s, "R");
  SEXP piv_s = cdid_ia_get_elt(fit_s, "pivot");
  if (TYPEOF(R_s) != REALSXP || XLENGTH(R_s) != (R_xlen_t) p * (R_xlen_t) p ||
      TYPEOF(piv_s) != INTSXP || XLENGTH(piv_s) < (R_xlen_t) p) {
    return 0;
  }
  const double *R = REAL(R_s);
  const int *piv = INTEGER(piv_s);
  for (int j = 0; j < p; ++j) {
    double d = R[j + j * p];
    if (!R_FINITE(d) || d == 0.0 || piv[j] < 1 || piv[j] > p) return 0;
  }
  *R_out = R;
  *piv_out = piv;
  return 1;
}

static int cdid_ia_gram_solve_rows(
    const double *R,
    const int *piv,
    int p,
    const double *rhs_mat,
    int B,
    double normalizer,
    double *out_mat) {
  if (p <= 0 || B <= 0) return 1;
  double v[256], y[256], z[256], inv_diag[256];
  if (p > 256) return 0;
  for (int k = 0; k < p; ++k) {
    inv_diag[k] = 1.0 / R[k + (size_t) k * (size_t) p];
  }
  for (int r = 0; r < B; ++r) {
    for (int k = 0; k < p; ++k) {
      double val = rhs_mat[r + (size_t) (piv[k] - 1) * (size_t) B];
      if (!R_FINITE(val)) return 0;
      v[k] = val;
    }
    for (int k = 0; k < p; ++k) {
      double s = v[k];
      const double *col_k = R + (size_t) k * (size_t) p;
      for (int j = 0; j < k; ++j) {
        s = fma(-col_k[j], y[j], s);
      }
      y[k] = s * inv_diag[k];
    }
    for (int k = p - 1; k >= 0; --k) {
      double s = y[k];
      for (int j = k + 1; j < p; ++j) {
        s = fma(-R[k + (size_t) j * (size_t) p], z[j], s);
      }
      double zk = s * inv_diag[k];
      if (!R_FINITE(zk)) return 0;
      z[k] = zk;
    }
    for (int k = 0; k < p; ++k) {
      double val = z[k] * normalizer;
      if (!R_FINITE(val)) return 0;
      out_mat[r + (size_t) (piv[k] - 1) * (size_t) B] = val;
    }
  }
  return 1;
}

static int cdid_ia_chol_solve_rows(
    const double *R,
    int p,
    const double *rhs_mat,
    int B,
    double *out_mat) {
  if (p <= 0 || B <= 0) return 1;
  double y[256], z[256], inv_diag[256];
  if (p > 256) return 0;
  for (int j = 0; j < p; ++j) {
    double d = R[j + j * p];
    if (!R_FINITE(d) || d == 0.0) return 0;
    inv_diag[j] = 1.0 / d;
  }
  for (int r = 0; r < B; ++r) {
    for (int k = 0; k < p; ++k) {
      double s = rhs_mat[r + (size_t) k * (size_t) B];
      if (!R_FINITE(s)) return 0;
      const double *col_k = R + (size_t) k * (size_t) p;
      for (int j = 0; j < k; ++j) {
        s = fma(-col_k[j], y[j], s);
      }
      y[k] = s * inv_diag[k];
    }
    for (int k = p - 1; k >= 0; --k) {
      double s = y[k];
      for (int j = k + 1; j < p; ++j) {
        s = fma(-R[k + (size_t) j * (size_t) p], z[j], s);
      }
      double zk = s * inv_diag[k];
      if (!R_FINITE(zk)) return 0;
      z[k] = zk;
      out_mat[r + (size_t) k * (size_t) B] = zk;
    }
  }
  return 1;
}

static int cdid_ia_resolve_pair_dates(
    SEXP pair_s,
    const double *all_times,
    int n_dates,
    int *ds_out,
    int *dt_out) {
  if (TYPEOF(pair_s) != VECSXP) return 0;
  SEXP ds_s = cdid_ia_get_elt(pair_s, "s_calendar_index");
  SEXP dt_s = cdid_ia_get_elt(pair_s, "t_calendar_index");
  int ds = -1, dt = -1;
  if ((TYPEOF(ds_s) == INTSXP || TYPEOF(ds_s) == REALSXP) && XLENGTH(ds_s) == 1) {
    ds = Rf_asInteger(ds_s);
  }
  if ((TYPEOF(dt_s) == INTSXP || TYPEOF(dt_s) == REALSXP) && XLENGTH(dt_s) == 1) {
    dt = Rf_asInteger(dt_s);
  }
  if (ds < 1 || ds > n_dates) {
    SEXP s_s = cdid_ia_get_elt(pair_s, "s");
    if (XLENGTH(s_s) == 1) {
      double sv = Rf_asReal(s_s);
      for (int i = 0; i < n_dates; ++i) {
        if (all_times[i] == sv) { ds = i + 1; break; }
      }
    }
  }
  if (dt < 1 || dt > n_dates) {
    SEXP t_s = cdid_ia_get_elt(pair_s, "t");
    if (XLENGTH(t_s) == 1) {
      double tv = Rf_asReal(t_s);
      for (int i = 0; i < n_dates; ++i) {
        if (all_times[i] == tv) { dt = i + 1; break; }
      }
    }
  }
  if (ds < 1 || ds > n_dates || dt < 1 || dt > n_dates) return 0;
  *ds_out = ds;
  *dt_out = dt;
  return 1;
}

/* 1. Fast C Design Matrix Finiteness + Repeat-Unit Match Check */
SEXP C_cdid_ia_check_design_matrix(
    SEXP X_s,
    SEXP Z_a_s,
    SEXP a_index_s,
    SEXP x_match_s) {
  if (TYPEOF(X_s) != REALSXP || Rf_isObject(X_s)) return R_NilValue;
  SEXP dim = Rf_getAttrib(X_s, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2) return R_NilValue;
  int N = INTEGER(dim)[0];
  int d = INTEGER(dim)[1];
  const double *X = REAL(X_s);
  R_xlen_t total = (R_xlen_t) N * (R_xlen_t) d;
  for (R_xlen_t i = 0; i < total; ++i) {
    if (!R_FINITE(X[i])) return Rf_ScalarInteger(1);
  }
  if (Z_a_s != R_NilValue && a_index_s != R_NilValue && x_match_s != R_NilValue) {
    if (TYPEOF(Z_a_s) != REALSXP || TYPEOF(a_index_s) != INTSXP || TYPEOF(x_match_s) != INTSXP) {
      return R_NilValue;
    }
    SEXP zdim = Rf_getAttrib(Z_a_s, R_DimSymbol);
    if (TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2) return R_NilValue;
    int n_a = INTEGER(zdim)[0];
    int p_z = INTEGER(zdim)[1];
    if (XLENGTH(a_index_s) != (R_xlen_t) N || XLENGTH(x_match_s) != (R_xlen_t) p_z) {
      return R_NilValue;
    }
    const double *Z_a = REAL(Z_a_s);
    const int *ai = INTEGER(a_index_s);
    const int *xm = INTEGER(x_match_s);
    for (int j = 0; j < p_z; ++j) {
      int col_x = xm[j] - 1;
      if (col_x < 0 || col_x >= d) return Rf_ScalarInteger(2);
      const double *x_col = X + (size_t) col_x * (size_t) N;
      const double *z_col = Z_a + (size_t) j * (size_t) n_a;
      for (int r = 0; r < N; ++r) {
        int a = ai[r];
        if (a > 0) {
          if (a > n_a || x_col[r] != z_col[a - 1]) {
            return Rf_ScalarInteger(2);
          }
        }
      }
    }
  }
  return Rf_ScalarInteger(0);
}

/* 2. Fast C Multiplier Bootstrap Helpers */
SEXP C_cdid_ia_multiplier_acc_if(
    SEXP S_mat_s,
    SEXP counts_s,
    SEXP if_mat_s,
    SEXP u_pos_s,
    SEXP unit_cluster_idx_s) {
  if (TYPEOF(S_mat_s) != REALSXP || TYPEOF(counts_s) != INTSXP ||
      TYPEOF(if_mat_s) != REALSXP || TYPEOF(u_pos_s) != INTSXP ||
      TYPEOF(unit_cluster_idx_s) != INTSXP) {
    return R_NilValue;
  }
  SEXP sdim = Rf_getAttrib(S_mat_s, R_DimSymbol);
  SEXP ifdim = Rf_getAttrib(if_mat_s, R_DimSymbol);
  if (TYPEOF(sdim) != INTSXP || XLENGTH(sdim) != 2 ||
      TYPEOF(ifdim) != INTSXP || XLENGTH(ifdim) != 2) {
    return R_NilValue;
  }
  int n_cl = INTEGER(sdim)[0];
  int K = INTEGER(sdim)[1];
  int B = INTEGER(ifdim)[0];
  int K_if = INTEGER(ifdim)[1];
  if (K != K_if || XLENGTH(counts_s) != (R_xlen_t) n_cl || XLENGTH(u_pos_s) != (R_xlen_t) B) {
    return R_NilValue;
  }
  double *S_mat = REAL(S_mat_s);
  int *counts = INTEGER(counts_s);
  const double *if_mat = REAL(if_mat_s);
  const int *u_pos = INTEGER(u_pos_s);
  const int *ucl = INTEGER(unit_cluster_idx_s);
  int n_units = (int) XLENGTH(unit_cluster_idx_s);

  for (int i = 0; i < B; ++i) {
    int up = u_pos[i];
    if (up == NA_INTEGER || up < 1 || up > n_units) continue;
    int cl = ucl[up - 1];
    if (cl == NA_INTEGER || cl < 1 || cl > n_cl) continue;
    int c0 = cl - 1;
    counts[c0] += 1;
    for (int k = 0; k < K; ++k) {
      double v = if_mat[i + (size_t) k * (size_t) B];
      if (!R_FINITE(v)) return Rf_ScalarLogical(0);
      S_mat[c0 + (size_t) k * (size_t) n_cl] += v;
    }
  }
  return Rf_ScalarLogical(1);
}

SEXP C_cdid_ia_check_sign_matrix(SEXP smat_s) {
  if (TYPEOF(smat_s) != REALSXP && TYPEOF(smat_s) != INTSXP) return Rf_ScalarLogical(0);
  R_xlen_t n = XLENGTH(smat_s);
  if (TYPEOF(smat_s) == REALSXP) {
    const double *v = REAL(smat_s);
    for (R_xlen_t i = 0; i < n; ++i) {
      double x = v[i];
      if (x != 1.0 && x != -1.0) return Rf_ScalarLogical(0);
    }
  } else {
    const int *v = INTEGER(smat_s);
    for (R_xlen_t i = 0; i < n; ++i) {
      int x = v[i];
      if (x != 1 && x != -1) return Rf_ScalarLogical(0);
    }
  }
  return Rf_ScalarLogical(1);
}

/* 3. Fast C Pure-Imputation IF Chunk Kernel (non-detailed) */
SEXP C_cdid_ia_imp_if_chunk(
    SEXP backend_s,
    SEXP st_s,
    SEXP panel_s,
    SEXP X_F_s,
    SEXP rows_s,
    SEXP links_s,
    SEXP sets_s,
    SEXP scale_s) {
  if (TYPEOF(backend_s) != VECSXP || TYPEOF(st_s) != VECSXP ||
      TYPEOF(panel_s) != VECSXP || TYPEOF(X_F_s) != REALSXP ||
      TYPEOF(rows_s) != INTSXP || TYPEOF(links_s) != INTSXP ||
      TYPEOF(sets_s) != VECSXP || TYPEOF(scale_s) != STRSXP || XLENGTH(scale_s) != 1) {
    return R_NilValue;
  }
  const void *vmax = vmaxget();

  int B = (int) XLENGTH(rows_s);
  int nl = (int) XLENGTH(links_s);
  const int *rows = INTEGER(rows_s);
  const int *links = INTEGER(links_s);
  int is_core_scale = (strcmp(CHAR(STRING_ELT(scale_s, 0)), "core") == 0);

  int N = Rf_asInteger(cdid_ia_get_elt(st_s, "N"));
  double lambda = Rf_asReal(cdid_ia_get_elt(st_s, "lambda"));
  SEXP cohorts_s = cdid_ia_get_elt(st_s, "cohorts");
  int G = (int) XLENGTH(cohorts_s);
  const int *cohort_index = INTEGER(cdid_ia_get_elt(st_s, "cohort_index"));
  const double *w_all = REAL(cdid_ia_get_elt(st_s, "w"));
  const int *a_index = INTEGER(cdid_ia_get_elt(st_s, "a_index"));
  const double *a_vec = REAL(cdid_ia_get_elt(st_s, "a"));
  const double *mass_F = REAL(cdid_ia_get_elt(st_s, "mass_F"));

  SEXP xdim = Rf_getAttrib(X_F_s, R_DimSymbol);
  if (TYPEOF(xdim) != INTSXP || XLENGTH(xdim) != 2 || INTEGER(xdim)[0] != N) {
    vmaxset(vmax);
    return R_NilValue;
  }
  int ncol_X = INTEGER(xdim)[1];
  const double *X_F = REAL(X_F_s);

  SEXP all_times_s = cdid_ia_get_elt(panel_s, "all_times");
  int n_dates = (int) XLENGTH(all_times_s);
  const double *all_times = REAL(all_times_s);
  SEXP post_s = cdid_ia_get_elt(panel_s, "postings");
  const int *post_a = INTEGER(cdid_ia_get_elt(post_s, "post_a_index"));
  const int *post_obs = INTEGER(cdid_ia_get_elt(post_s, "post_obs_index"));
  const int *date_start = INTEGER(cdid_ia_get_elt(post_s, "date_start"));
  const int *date_end = INTEGER(cdid_ia_get_elt(post_s, "date_end"));
  const int *date_count = INTEGER(cdid_ia_get_elt(post_s, "date_count"));
  SEXP y_obs_s = cdid_ia_get_elt(cdid_ia_get_elt(panel_s, "observations"), "y");
  const double *y_obs = REAL(y_obs_s);
  R_xlen_t n_y = XLENGTH(y_obs_s);
  R_xlen_t n_post = XLENGTH(cdid_ia_get_elt(post_s, "post_a_index"));

  SEXP recipe_s = cdid_ia_get_elt(backend_s, "outcome_recipe");
  int intercept = (Rf_asLogical(cdid_ia_get_elt(recipe_s, "intercept")) == 1);
  SEXP trt_recipe_s = cdid_ia_get_elt(backend_s, "treatment_recipe");
  const char *est_str = CHAR(STRING_ELT(cdid_ia_get_elt(backend_s, "est_method"), 0));
  int mode_est = (strcmp(est_str, "ipw") == 0) ? 0 : ((strcmp(est_str, "reg") == 0) ? 1 : 2);

  SEXP cells_s = cdid_ia_get_elt(backend_s, "cells");
  SEXP groups_s = cdid_ia_get_elt(backend_s, "groups");
  SEXP shared_or_s = cdid_ia_get_elt(backend_s, "shared_or");
  SEXP treatments_s = cdid_ia_get_elt(backend_s, "treatments");
  SEXP links_all_s = cdid_ia_get_elt(backend_s, "links");
  int n_cells = (int) XLENGTH(cells_s);
  int n_groups = (int) XLENGTH(groups_s);

  int *ci = (int *) R_alloc((size_t) B, sizeof(int));
  double *w = (double *) R_alloc((size_t) B, sizeof(double));
  int *ai = (int *) R_alloc((size_t) B, sizeof(int));
  for (int r = 0; r < B; ++r) {
    int row0 = rows[r] - 1;
    ci[r] = cohort_index[row0];
    w[r] = w_all[row0];
    ai[r] = a_index[row0];
  }

  /* 1. Repeat-rate IFs xi_a (B x G) */
  double *xi_a = (double *) R_alloc((size_t) B * (size_t) G, sizeof(double));
  memset(xi_a, 0, (size_t) B * (size_t) G * sizeof(double));
  for (int r = 0; r < B; ++r) {
    int k = ci[r];
    if (k >= 1 && k <= G) {
      double mf = mass_F[k - 1];
      if (mf > 0.0) {
        double A_r = (ai[r] > 0) ? 1.0 : 0.0;
        xi_a[r + (size_t) (k - 1) * (size_t) B] = w[r] * (A_r - a_vec[k - 1]) / mf;
      }
    }
  }

  /* 2. Cell scores & group WLS forward IFs (alpha_if and xi_beta_imp) */
  SEXP cids_s = cdid_ia_get_elt(sets_s, "cids");
  SEXP gids_s = cdid_ia_get_elt(sets_s, "gids");
  SEXP bids_s = cdid_ia_get_elt(sets_s, "bids");
  SEXP tids_s = cdid_ia_get_elt(sets_s, "tids");
  int n_cids = (int) XLENGTH(cids_s);
  int n_gids = (int) XLENGTH(gids_s);
  int n_bids = (int) XLENGTH(bids_s);
  int n_tids = (int) XLENGTH(tids_s);
  const int *cids = INTEGER(cids_s);
  const int *gids = INTEGER(gids_s);
  const int *bids = INTEGER(bids_s);
  const int *tids = INTEGER(tids_s);

  double **alpha_if = (double **) R_alloc((size_t) (n_cells + 1), sizeof(double *));
  for (int c = 0; c <= n_cells; ++c) alpha_if[c] = NULL;
  double **h_beta = (double **) R_alloc((size_t) (n_groups + 1), sizeof(double *));
  double **xi_beta = (double **) R_alloc((size_t) (n_groups + 1), sizeof(double *));
  int *p_group = (int *) R_alloc((size_t) (n_groups + 1), sizeof(int));
  for (int g = 0; g <= n_groups; ++g) {
    h_beta[g] = NULL;
    xi_beta[g] = NULL;
    p_group[g] = 0;
  }

  for (int ig = 0; ig < n_gids; ++ig) {
    int gid = gids[ig];
    if (gid < 1 || gid > n_groups) { vmaxset(vmax); return R_NilValue; }
    SEXP gf = VECTOR_ELT(groups_s, gid - 1);
    int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
    p_group[gid] = pg;
    if (pg > 0) {
      h_beta[gid] = (double *) R_alloc((size_t) B * (size_t) pg, sizeof(double));
      xi_beta[gid] = (double *) R_alloc((size_t) B * (size_t) pg, sizeof(double));
      memset(h_beta[gid], 0, (size_t) B * (size_t) pg * sizeof(double));
      memset(xi_beta[gid], 0, (size_t) B * (size_t) pg * sizeof(double));
    }
  }

  double d_buf[256];
  for (int ic = 0; ic < n_cids; ++ic) {
    int cid = cids[ic];
    if (cid < 1 || cid > n_cells) { vmaxset(vmax); return R_NilValue; }
    SEXP cf = VECTOR_ELT(cells_s, cid - 1);
    int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
    int c_coh = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
    double L_val = Rf_asReal(cdid_ia_get_elt(cf, "L"));
    double alpha_c = intercept ? Rf_asReal(cdid_ia_get_elt(cf, "alpha")) : 0.0;
    double mass_c = Rf_asReal(cdid_ia_get_elt(cf, "mass"));
    SEXP gf = VECTOR_ELT(groups_s, gid - 1);
    int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
    if (pg > 256) { vmaxset(vmax); return R_NilValue; }
    const int *kept_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf, "kept_cols_F")) : NULL;
    const double *beta_g = (pg > 0) ? REAL(cdid_ia_get_elt(gf, "beta")) : NULL;
    const double *mean_d = (pg > 0 && intercept) ? REAL(cdid_ia_get_elt(cf, "mean_d")) : NULL;

    int ds = 0, dt = 0;
    if (!cdid_ia_resolve_pair_dates(cdid_ia_get_elt(cf, "pair"), all_times, n_dates, &ds, &dt)) {
      vmaxset(vmax);
      return R_NilValue;
    }
    int has_pair = (date_count[ds - 1] > 0 && date_count[dt - 1] > 0);
    int s_lo = date_start[ds - 1] - 1, s_hi = date_end[ds - 1] - 1;
    int t_lo = date_start[dt - 1] - 1, t_hi = date_end[dt - 1] - 1;
    if (has_pair && (s_lo < 0 || s_hi < s_lo || (R_xlen_t) s_hi >= n_post ||
                     t_lo < 0 || t_hi < t_lo || (R_xlen_t) t_hi >= n_post)) {
      vmaxset(vmax);
      return R_NilValue;
    }

    double *a_if_c = (double *) R_alloc((size_t) B, sizeof(double));
    memset(a_if_c, 0, (size_t) B * sizeof(double));
    alpha_if[cid] = a_if_c;
    double *hb = (gid >= 1 && gid <= n_groups) ? h_beta[gid] : NULL;

    if (has_pair) {
      int hs_hint = s_lo, ht_hint = t_lo;
      for (int r = 0; r < B; ++r) {
        if (ci[r] != c_coh || ai[r] <= 0) continue;
        int hs = cdid_ia_find_posting(post_a, s_lo, s_hi, ai[r], &hs_hint);
        if (hs < 0) continue;
        int ht = cdid_ia_find_posting(post_a, t_lo, t_hi, ai[r], &ht_hint);
        if (ht < 0) continue;
        int os = post_obs[hs], ot = post_obs[ht];
        if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y) {
          vmaxset(vmax);
          return R_NilValue;
        }
        double dy = y_obs[ot - 1] - y_obs[os - 1];
        if (!R_FINITE(dy)) { vmaxset(vmax); return R_NilValue; }
        int row0 = rows[r] - 1;
        double jhat_sl = 0.0;
        for (int m = 0; m < pg; ++m) {
          int col0 = kept_F[m] - 1;
          double d_m = L_val * X_F[row0 + (size_t) col0 * (size_t) N];
          d_buf[m] = d_m;
          jhat_sl += d_m * beta_g[m];
        }
        double jhat = alpha_c + jhat_sl;
        double s_rr = w[r] * (dy - jhat);
        if (!R_FINITE(s_rr)) { vmaxset(vmax); return R_NilValue; }
        if (intercept) {
          a_if_c[r] = s_rr / mass_c;
        }
        if (pg > 0 && hb != NULL) {
          for (int m = 0; m < pg; ++m) {
            double dt_m = intercept ? (d_buf[m] - mean_d[m]) : d_buf[m];
            hb[r + (size_t) m * (size_t) B] += s_rr * dt_m;
          }
        }
      }
    }
  }

  for (int ig = 0; ig < n_gids; ++ig) {
    int gid = gids[ig];
    int pg = p_group[gid];
    if (pg > 0) {
      SEXP gf = VECTOR_ELT(groups_s, gid - 1);
      const double *R_g = NULL;
      const int *piv_g = NULL;
      if (!cdid_ia_get_wls_fit(cdid_ia_get_elt(gf, "fit"), &R_g, &piv_g, pg)) {
        vmaxset(vmax);
        return R_NilValue;
      }
      if (!cdid_ia_gram_solve_rows(R_g, piv_g, pg, h_beta[gid], B, (double) N, xi_beta[gid])) {
        vmaxset(vmax);
        return R_NilValue;
      }
      if (intercept) {
        SEXP g_cids_s = cdid_ia_get_elt(gf, "cell_ids");
        int n_gc = (int) XLENGTH(g_cids_s);
        const int *g_cids = INTEGER(g_cids_s);
        const double *xb = xi_beta[gid];
        for (int k = 0; k < n_gc; ++k) {
          int cid = g_cids[k];
          if (cid >= 1 && cid <= n_cells && alpha_if[cid] != NULL) {
            double *a_if_c = alpha_if[cid];
            const double *mean_d = REAL(cdid_ia_get_elt(VECTOR_ELT(cells_s, cid - 1), "mean_d"));
            for (int m = 0; m < pg; ++m) {
              double md = mean_d[m];
              if (md != 0.0) {
                const double *xb_m = xb + (size_t) m * (size_t) B;
                for (int r = 0; r < B; ++r) a_if_c[r] -= xb_m[r] * md;
              }
            }
          }
        }
      }
    }
  }

  /* 3. Shared / Separate Counterfactual OR IFs (direct, a, j channels) */
  int has_or = (shared_or_s != R_NilValue && n_bids > 0);
  int separate_or = 0;
  int n_or_blocks = 0;
  SEXP or_blocks_s = R_NilValue;
  if (has_or) {
    SEXP op_s = cdid_ia_get_elt(shared_or_s, "outcome_pooling");
    if (TYPEOF(op_s) == STRSXP && XLENGTH(op_s) == 1 &&
        strcmp(CHAR(STRING_ELT(op_s, 0)), "separate") == 0) {
      separate_or = 1;
    }
    or_blocks_s = cdid_ia_get_elt(shared_or_s, "blocks");
    n_or_blocks = (int) XLENGTH(or_blocks_s);
  }

  double **or_alpha_dir = NULL, **or_alpha_a = NULL, **or_alpha_j = NULL;
  double **or_beta_dir = NULL, **or_beta_a = NULL, **or_beta_j = NULL;
  int *p_or_by_bid = NULL;
  int **kept_or_by_bid = NULL;
  const double **beta_or_by_bid = NULL;

  if (has_or) {
    or_alpha_dir = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_alpha_a   = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_alpha_j   = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_beta_dir  = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_beta_a    = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_beta_j    = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    p_or_by_bid  = (int *) R_alloc((size_t) (n_or_blocks + 1), sizeof(int));
    kept_or_by_bid = (int **) R_alloc((size_t) (n_or_blocks + 1), sizeof(int *));
    beta_or_by_bid = (const double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(const double *));
    for (int b = 0; b <= n_or_blocks; ++b) {
      or_alpha_dir[b] = or_alpha_a[b] = or_alpha_j[b] = NULL;
      or_beta_dir[b] = or_beta_a[b] = or_beta_j[b] = NULL;
      p_or_by_bid[b] = 0;
      kept_or_by_bid[b] = NULL;
      beta_or_by_bid[b] = NULL;
    }

    const int *slope_cols_F_all = separate_or ? INTEGER(cdid_ia_get_elt(shared_or_s, "slope_cols_F")) : NULL;
    SEXP og_groups_s = R_NilValue;
    const int *bm_gid = NULL;
    if (separate_or) {
      SEXP og_s = cdid_ia_get_elt(shared_or_s, "outcome_groups");
      og_groups_s = cdid_ia_get_elt(og_s, "groups");
      bm_gid = INTEGER(cdid_ia_get_elt(cdid_ia_get_elt(og_s, "block_map"), "group_id"));
    }

    int p_or_pool = separate_or ? 0 : Rf_asInteger(cdid_ia_get_elt(shared_or_s, "p"));
    if (p_or_pool > 256) { vmaxset(vmax); return R_NilValue; }
    int *kept_pool = (!separate_or && p_or_pool > 0) ? INTEGER(cdid_ia_get_elt(shared_or_s, "kept_cols_F")) : NULL;
    const double *beta_pool = (!separate_or && p_or_pool > 0) ? REAL(cdid_ia_get_elt(shared_or_s, "beta")) : NULL;

    double *pool_rhs_b_dir = (!separate_or && p_or_pool > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double)) : NULL;
    double *pool_rhs_b_a   = (!separate_or && p_or_pool > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double)) : NULL;
    double *pool_rhs_b_j   = (!separate_or && p_or_pool > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double)) : NULL;
    if (!separate_or && p_or_pool > 0) {
      memset(pool_rhs_b_dir, 0, (size_t) B * (size_t) p_or_pool * sizeof(double));
      memset(pool_rhs_b_a,   0, (size_t) B * (size_t) p_or_pool * sizeof(double));
      memset(pool_rhs_b_j,   0, (size_t) B * (size_t) p_or_pool * sizeof(double));
    }

    /* Also make sure every link's or_block_id has alpha allocated in pooled mode */
    int *active_bid_mask = (int *) R_alloc((size_t) (n_or_blocks + 1), sizeof(int));
    memset(active_bid_mask, 0, (size_t) (n_or_blocks + 1) * sizeof(int));
    for (int ib = 0; ib < n_bids; ++ib) {
      int bid = bids[ib];
      if (bid >= 1 && bid <= n_or_blocks) active_bid_mask[bid] = 1;
    }
    for (int j = 0; j < nl; ++j) {
      SEXP lf = VECTOR_ELT(links_all_s, links[j] - 1);
      int bid = Rf_asInteger(cdid_ia_get_elt(lf, "or_block_id"));
      if (bid >= 1 && bid <= n_or_blocks) active_bid_mask[bid] = 1;
    }

    for (int bid = 1; bid <= n_or_blocks; ++bid) {
      if (!active_bid_mask[bid]) continue;
      SEXP bf = VECTOR_ELT(or_blocks_s, bid - 1);
      int p_or = p_or_pool;
      int *kept_or_F = kept_pool;
      const double *beta_or = beta_pool;
      SEXP fit_or_s = R_NilValue;
      if (separate_or) {
        SEXP og_g = VECTOR_ELT(og_groups_s, bm_gid[bid - 1] - 1);
        SEXP keep_s = cdid_ia_get_elt(og_g, "keep");
        p_or = (int) XLENGTH(keep_s);
        if (p_or > 256) { vmaxset(vmax); return R_NilValue; }
        const int *keep_idx = INTEGER(keep_s);
        kept_or_F = (p_or > 0) ? (int *) R_alloc((size_t) p_or, sizeof(int)) : NULL;
        for (int m = 0; m < p_or; ++m) kept_or_F[m] = slope_cols_F_all[keep_idx[m] - 1];
        beta_or = (p_or > 0) ? REAL(cdid_ia_get_elt(og_g, "beta")) : NULL;
        fit_or_s = cdid_ia_get_elt(og_g, "fit");
      }
      p_or_by_bid[bid] = p_or;
      kept_or_by_bid[bid] = kept_or_F;
      beta_or_by_bid[bid] = beta_or;

      double *ra_dir = (double *) R_alloc((size_t) B, sizeof(double));
      double *ra_a   = (double *) R_alloc((size_t) B, sizeof(double));
      double *ra_j   = (double *) R_alloc((size_t) B, sizeof(double));
      memset(ra_dir, 0, (size_t) B * sizeof(double));
      memset(ra_a,   0, (size_t) B * sizeof(double));
      memset(ra_j,   0, (size_t) B * sizeof(double));
      or_alpha_dir[bid] = ra_dir;
      or_alpha_a[bid]   = ra_a;
      or_alpha_j[bid]   = ra_j;

      double *rb_dir = separate_or ? ((p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL) : pool_rhs_b_dir;
      double *rb_a   = separate_or ? ((p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL) : pool_rhs_b_a;
      double *rb_j   = separate_or ? ((p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL) : pool_rhs_b_j;
      if (separate_or && p_or > 0) {
        memset(rb_dir, 0, (size_t) B * (size_t) p_or * sizeof(double));
        memset(rb_a,   0, (size_t) B * (size_t) p_or * sizeof(double));
        memset(rb_j,   0, (size_t) B * (size_t) p_or * sizeof(double));
      }

      /* Check if bid is in sets$bids */
      int in_wanted = 0;
      for (int ib = 0; ib < n_bids; ++ib) {
        if (bids[ib] == bid) { in_wanted = 1; break; }
      }
      if (!in_wanted) continue;

      double L_b = Rf_asReal(cdid_ia_get_elt(bf, "L"));
      double alpha_b = intercept ? Rf_asReal(cdid_ia_get_elt(bf, "alpha")) : 0.0;
      const double *mean_d_b = (p_or > 0 && intercept) ? REAL(cdid_ia_get_elt(bf, "mean_d")) : NULL;
      SEXP ctrl_cids_s = cdid_ia_get_elt(bf, "control_cell_ids");
      int n_kb = (int) XLENGTH(ctrl_cids_s);
      const int *ctrl_cids = INTEGER(ctrl_cids_s);
      const double *a_force_al = REAL(cdid_ia_get_elt(bf, "a_force_alpha"));
      const double *a_force_be = (p_or > 0) ? REAL(cdid_ia_get_elt(bf, "a_force_beta")) : NULL;
      const double *j_force_aa = REAL(cdid_ia_get_elt(bf, "j_force_alpha_alpha"));
      SEXP j_force_ab_s = cdid_ia_get_elt(bf, "j_force_alpha_beta");
      const double *j_force_ba = (p_or > 0) ? REAL(cdid_ia_get_elt(bf, "j_force_beta_alpha")) : NULL;
      SEXP j_force_bb_s = cdid_ia_get_elt(bf, "j_force_beta_beta");

      for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
        int cid = ctrl_cids[idx_k];
        SEXP cf = VECTOR_ELT(cells_s, cid - 1);
        int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
        SEXP gf = VECTOR_ELT(groups_s, gid - 1);
        int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
        const int *kept_imp_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf, "kept_cols_F")) : NULL;
        const double *beta_imp = (pg > 0) ? REAL(cdid_ia_get_elt(gf, "beta")) : NULL;
        double alpha_c = Rf_asReal(cdid_ia_get_elt(cf, "alpha"));
        int ci_k = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
        double a_k = a_vec[ci_k - 1];

        for (int r = 0; r < B; ++r) {
          if (ci[r] != ci_k) continue;
          int row0 = rows[r] - 1;
          double w_syn = w[r] * a_k;
          double jhat_sl = 0.0;
          for (int m = 0; m < pg; ++m) {
            jhat_sl += (L_b * X_F[row0 + (size_t) (kept_imp_F[m] - 1) * (size_t) N]) * beta_imp[m];
          }
          double jhat = alpha_c + jhat_sl;
          double mhat_sl = 0.0;
          for (int m = 0; m < p_or; ++m) {
            double d_or_m = L_b * X_F[row0 + (size_t) (kept_or_F[m] - 1) * (size_t) N];
            d_buf[m] = d_or_m;
            mhat_sl += d_or_m * beta_or[m];
          }
          double mhat = alpha_b + mhat_sl;
          double wr = w_syn * (jhat - mhat);
          if (intercept) ra_dir[r] += wr;
          for (int m = 0; m < p_or; ++m) {
            double dt_m = intercept ? (d_buf[m] - mean_d_b[m]) : d_buf[m];
            rb_dir[r + (size_t) m * (size_t) B] += wr * dt_m;
          }
        }

        const double *xi_ak = xi_a + (size_t) (ci_k - 1) * (size_t) B;
        if (intercept) {
          double afa = a_force_al[idx_k];
          if (afa != 0.0) {
            for (int r = 0; r < B; ++r) ra_a[r] += afa * xi_ak[r];
          }
        }
        if (p_or > 0) {
          const double *afb = a_force_be + (size_t) idx_k * (size_t) p_or;
          for (int m = 0; m < p_or; ++m) {
            double afb_m = afb[m];
            if (afb_m != 0.0) {
              double *dst = rb_a + (size_t) m * (size_t) B;
              for (int r = 0; r < B; ++r) dst[r] += xi_ak[r] * afb_m;
            }
          }
        }

        const double *xi_ac = alpha_if[cid];
        const double *xi_br = xi_beta[gid];
        if (intercept) {
          double jfaa = j_force_aa[idx_k];
          if (jfaa != 0.0 && xi_ac != NULL) {
            for (int r = 0; r < B; ++r) ra_j[r] += jfaa * xi_ac[r];
          }
          if (pg > 0 && xi_br != NULL) {
            const double *jfab = REAL(VECTOR_ELT(j_force_ab_s, idx_k));
            for (int m = 0; m < pg; ++m) {
              double v = jfab[m];
              if (v != 0.0) {
                const double *xcol = xi_br + (size_t) m * (size_t) B;
                for (int r = 0; r < B; ++r) ra_j[r] += xcol[r] * v;
              }
            }
          }
        }
        if (p_or > 0) {
          if (intercept && xi_ac != NULL) {
            const double *jfba = j_force_ba + (size_t) idx_k * (size_t) p_or;
            for (int q = 0; q < p_or; ++q) {
              double v = jfba[q];
              if (v != 0.0) {
                double *dst = rb_j + (size_t) q * (size_t) B;
                for (int r = 0; r < B; ++r) dst[r] += xi_ac[r] * v;
              }
            }
          }
          if (pg > 0 && xi_br != NULL) {
            const double *jfbb = REAL(VECTOR_ELT(j_force_bb_s, idx_k));
            for (int q = 0; q < p_or; ++q) {
              double *dst = rb_j + (size_t) q * (size_t) B;
              for (int m = 0; m < pg; ++m) {
                double v = jfbb[q + (size_t) m * (size_t) p_or];
                if (v != 0.0) {
                  const double *xcol = xi_br + (size_t) m * (size_t) B;
                  for (int r = 0; r < B; ++r) dst[r] += xcol[r] * v;
                }
              }
            }
          }
        }
      }

      if (separate_or) {
        if (intercept) {
          double sma = Rf_asReal(cdid_ia_get_elt(bf, "score_mean_alpha"));
          for (int r = 0; r < B; ++r) ra_dir[r] -= sma;
        }
        if (p_or > 0) {
          const double *smb = REAL(cdid_ia_get_elt(bf, "score_mean_beta"));
          for (int m = 0; m < p_or; ++m) {
            double smbm = smb[m];
            double *dst = rb_dir + (size_t) m * (size_t) B;
            for (int r = 0; r < B; ++r) dst[r] -= smbm;
          }
          double *xb_dir = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
          double *xb_a   = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
          double *xb_j   = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
          const double *R_or = NULL;
          const int *piv_or = NULL;
          if (!cdid_ia_get_wls_fit(fit_or_s, &R_or, &piv_or, p_or) ||
              !cdid_ia_gram_solve_rows(R_or, piv_or, p_or, rb_dir, B, (double) N, xb_dir) ||
              !cdid_ia_gram_solve_rows(R_or, piv_or, p_or, rb_a,   B, (double) N, xb_a)   ||
              !cdid_ia_gram_solve_rows(R_or, piv_or, p_or, rb_j,   B, (double) N, xb_j)) {
            vmaxset(vmax);
            return R_NilValue;
          }
          or_beta_dir[bid] = xb_dir;
          or_beta_a[bid]   = xb_a;
          or_beta_j[bid]   = xb_j;
        }
        if (intercept) {
          double mass_b = Rf_asReal(cdid_ia_get_elt(bf, "mass"));
          for (int r = 0; r < B; ++r) {
            ra_dir[r] /= mass_b;
            ra_a[r]   /= mass_b;
            ra_j[r]   /= mass_b;
          }
          if (p_or > 0) {
            for (int m = 0; m < p_or; ++m) {
              double md = mean_d_b[m];
              if (md != 0.0) {
                const double *bd = or_beta_dir[bid] + (size_t) m * (size_t) B;
                const double *ba = or_beta_a[bid]   + (size_t) m * (size_t) B;
                const double *bj = or_beta_j[bid]   + (size_t) m * (size_t) B;
                for (int r = 0; r < B; ++r) {
                  ra_dir[r] -= bd[r] * md;
                  ra_a[r]   -= ba[r] * md;
                  ra_j[r]   -= bj[r] * md;
                }
              }
            }
          }
        }
      }
    }

    if (!separate_or) {
      double *xb_dir = NULL, *xb_a = NULL, *xb_j = NULL;
      if (p_or_pool > 0) {
        xb_dir = (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double));
        xb_a   = (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double));
        xb_j   = (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double));
        const double *R_or = NULL;
        const int *piv_or = NULL;
        if (!cdid_ia_get_wls_fit(cdid_ia_get_elt(shared_or_s, "fit"), &R_or, &piv_or, p_or_pool) ||
            !cdid_ia_gram_solve_rows(R_or, piv_or, p_or_pool, pool_rhs_b_dir, B, (double) N, xb_dir) ||
            !cdid_ia_gram_solve_rows(R_or, piv_or, p_or_pool, pool_rhs_b_a,   B, (double) N, xb_a)   ||
            !cdid_ia_gram_solve_rows(R_or, piv_or, p_or_pool, pool_rhs_b_j,   B, (double) N, xb_j)) {
          vmaxset(vmax);
          return R_NilValue;
        }
      }
      for (int bid = 1; bid <= n_or_blocks; ++bid) {
        if (!active_bid_mask[bid]) continue;
        or_beta_dir[bid] = xb_dir;
        or_beta_a[bid]   = xb_a;
        or_beta_j[bid]   = xb_j;
        if (intercept) {
          SEXP bf = VECTOR_ELT(or_blocks_s, bid - 1);
          double mass_b = Rf_asReal(cdid_ia_get_elt(bf, "mass"));
          double *ra_dir = or_alpha_dir[bid];
          double *ra_a   = or_alpha_a[bid];
          double *ra_j   = or_alpha_j[bid];
          for (int r = 0; r < B; ++r) {
            ra_dir[r] /= mass_b;
            ra_a[r]   /= mass_b;
            ra_j[r]   /= mass_b;
          }
          if (p_or_pool > 0) {
            const double *mean_d_b = REAL(cdid_ia_get_elt(bf, "mean_d"));
            for (int m = 0; m < p_or_pool; ++m) {
              double md = mean_d_b[m];
              if (md != 0.0) {
                const double *bd = xb_dir + (size_t) m * (size_t) B;
                const double *ba = xb_a   + (size_t) m * (size_t) B;
                const double *bj = xb_j   + (size_t) m * (size_t) B;
                for (int r = 0; r < B; ++r) {
                  ra_dir[r] -= bd[r] * md;
                  ra_a[r]   -= ba[r] * md;
                  ra_j[r]   -= bj[r] * md;
                }
              }
            }
          }
        }
      }
    }
  }

  /* 4. Synthetic-weight treatment propensity logit IFs (direct and a channels) + cached odds */
  int has_trt = (treatments_s != R_NilValue && n_tids > 0);
  int n_trt_tot = (treatments_s != R_NilValue) ? (int) XLENGTH(treatments_s) : 0;
  int p_trt = (trt_recipe_s != R_NilValue) ? Rf_asInteger(cdid_ia_get_elt(trt_recipe_s, "p")) : 0;
  if (p_trt > 256) { vmaxset(vmax); return R_NilValue; }
  const int *trt_cols_F = (p_trt > 0) ? INTEGER(cdid_ia_get_elt(trt_recipe_s, "index")) : NULL;

  double **trt_dir = NULL, **trt_a = NULL, **trt_odds = NULL;
  if (has_trt) {
    trt_dir  = (double **) R_alloc((size_t) (n_trt_tot + 1), sizeof(double *));
    trt_a    = (double **) R_alloc((size_t) (n_trt_tot + 1), sizeof(double *));
    trt_odds = (double **) R_alloc((size_t) (n_trt_tot + 1), sizeof(double *));
    for (int t = 0; t <= n_trt_tot; ++t) {
      trt_dir[t] = trt_a[t] = trt_odds[t] = NULL;
    }
    double *score_dir = (p_trt > 0) ? (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double)) : NULL;
    double *score_a   = (p_trt > 0) ? (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double)) : NULL;
    int *in_all_coh   = (int *) R_alloc((size_t) (G + 1), sizeof(int));

    for (int it = 0; it < n_tids; ++it) {
      int tid = tids[it];
      if (tid < 1 || tid > n_trt_tot) { vmaxset(vmax); return R_NilValue; }
      SEXP tf = VECTOR_ELT(treatments_s, tid - 1);
      int target_coh = Rf_asInteger(cdid_ia_get_elt(tf, "target_cohort_index"));
      SEXP all_coh_s = cdid_ia_get_elt(tf, "all_cohort_indices");
      int n_ac = (int) XLENGTH(all_coh_s);
      const int *all_coh = INTEGER(all_coh_s);
      const double *coef_trt = (p_trt > 0) ? REAL(cdid_ia_get_elt(tf, "coefficients")) : NULL;

      /* Check representative reuse among already-evaluated tids */
      int rep_tid = 0;
      for (int it2 = 0; it2 < it; ++it2) {
        int tid2 = tids[it2];
        SEXP tf2 = VECTOR_ELT(treatments_s, tid2 - 1);
        if (Rf_asInteger(cdid_ia_get_elt(tf2, "target_cohort_index")) != target_coh) continue;
        SEXP ac2_s = cdid_ia_get_elt(tf2, "all_cohort_indices");
        if ((int) XLENGTH(ac2_s) != n_ac) continue;
        if (memcmp(INTEGER(ac2_s), all_coh, (size_t) n_ac * sizeof(int)) != 0) continue;
        if (p_trt > 0 && memcmp(REAL(cdid_ia_get_elt(tf2, "coefficients")), coef_trt, (size_t) p_trt * sizeof(double)) != 0) continue;
        rep_tid = tid2;
        break;
      }
      if (rep_tid > 0) {
        trt_dir[tid]  = trt_dir[rep_tid];
        trt_a[tid]    = trt_a[rep_tid];
        trt_odds[tid] = trt_odds[rep_tid];
        continue;
      }

      memset(in_all_coh, 0, (size_t) (G + 1) * sizeof(int));
      for (int k = 0; k < n_ac; ++k) {
        int c = all_coh[k];
        if (c >= 1 && c <= G) in_all_coh[c] = 1;
      }

      double *odds_vec = (double *) R_alloc((size_t) B, sizeof(double));
      memset(odds_vec, 0, (size_t) B * sizeof(double));
      trt_odds[tid] = odds_vec;

      if (p_trt > 0) {
        memset(score_dir, 0, (size_t) B * (size_t) p_trt * sizeof(double));
        memset(score_a,   0, (size_t) B * (size_t) p_trt * sizeof(double));
        for (int r = 0; r < B; ++r) {
          int c = ci[r];
          if (c < 1 || c > G || !in_all_coh[c]) continue;
          int row0 = rows[r] - 1;
          double eta = 0.0;
          for (int m = 0; m < p_trt; ++m) {
            double x_rm = X_F[row0 + (size_t) (trt_cols_F[m] - 1) * (size_t) N];
            d_buf[m] = x_rm;
            eta = fma(x_rm, coef_trt[m], eta);
          }
          if (!R_FINITE(eta)) { vmaxset(vmax); return R_NilValue; }
          double pr = 0.0, od = 0.0;
          cdid_ia_clip_logit(eta, &pr, &od);
          odds_vec[r] = od;
          double y_r = (c == target_coh) ? 1.0 : 0.0;
          double sc = (w[r] * a_vec[c - 1]) * (y_r - pr);
          for (int m = 0; m < p_trt; ++m) {
            score_dir[r + (size_t) m * (size_t) B] = sc * d_buf[m];
          }
        }
        const double *a_forc = REAL(cdid_ia_get_elt(tf, "a_forcing"));
        for (int k = 0; k < n_ac; ++k) {
          int c = all_coh[k];
          const double *xi_ac = xi_a + (size_t) (c - 1) * (size_t) B;
          for (int m = 0; m < p_trt; ++m) {
            double af = a_forc[m + (size_t) k * (size_t) p_trt];
            double *dst = score_a + (size_t) m * (size_t) B;
            for (int r = 0; r < B; ++r) dst[r] = fma(xi_ac[r], af, dst[r]);
          }
        }
        double *xg_dir = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
        double *xg_a   = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
        const double *H_chol = REAL(cdid_ia_get_elt(tf, "H_chol"));
        if (!cdid_ia_chol_solve_rows(H_chol, p_trt, score_dir, B, xg_dir) ||
            !cdid_ia_chol_solve_rows(H_chol, p_trt, score_a,   B, xg_a)) {
          vmaxset(vmax);
          return R_NilValue;
        }
        trt_dir[tid] = xg_dir;
        trt_a[tid]   = xg_a;
      } else {
        for (int r = 0; r < B; ++r) odds_vec[r] = 1.0;
      }
    }
  }

  /* 5. Link IFs assembly */
  SEXP phi_s = PROTECT(Rf_allocMatrix(REALSXP, B, nl));
  SEXP link_names_s = PROTECT(Rf_allocVector(STRSXP, nl));
  double *phi = REAL(phi_s);

  double *w_plus   = (double *) R_alloc((size_t) B, sizeof(double));
  double *j_g_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  double *m_g_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  double *j_T_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  double *w_minus  = (double *) R_alloc((size_t) B, sizeof(double));
  double *j_k_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  double *m_k_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  double *a_C_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  double *j_C_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  double *a_Bo_vec = (double *) R_alloc((size_t) B, sizeof(double));
  int *ctrl_coh_map = (int *) R_alloc((size_t) (G + 1), sizeof(int));

  for (int j = 0; j < nl; ++j) {
    int ell = links[j];
    SEXP lf = VECTOR_ELT(links_all_s, ell - 1);
    SET_STRING_ELT(link_names_s, j, STRING_ELT(cdid_ia_get_elt(lf, "link_id"), 0));

    int bid = Rf_asInteger(cdid_ia_get_elt(lf, "or_block_id"));
    int t_cid = Rf_asInteger(cdid_ia_get_elt(lf, "target_cell_id"));
    int ci_g = Rf_asInteger(cdid_ia_get_elt(lf, "target_cohort_index"));
    double a_g = a_vec[ci_g - 1];
    double L_lk = Rf_asReal(cdid_ia_get_elt(lf, "L"));
    double d_plus = Rf_asReal(cdid_ia_get_elt(lf, "d_plus"));
    double T_gl = Rf_asReal(cdid_ia_get_elt(lf, "T_gl"));
    double a_coeff_T_g = Rf_asReal(cdid_ia_get_elt(lf, "a_coeff_T_g"));

    SEXP cf_g = VECTOR_ELT(cells_s, t_cid - 1);
    int gid_g = Rf_asInteger(cdid_ia_get_elt(cf_g, "group_id"));
    SEXP gf_g = VECTOR_ELT(groups_s, gid_g - 1);
    int pg_g = Rf_asInteger(cdid_ia_get_elt(gf_g, "p"));
    const int *kept_g_F = (pg_g > 0) ? INTEGER(cdid_ia_get_elt(gf_g, "kept_cols_F")) : NULL;
    const double *beta_g = (pg_g > 0) ? REAL(cdid_ia_get_elt(gf_g, "beta")) : NULL;
    double alpha_cg = Rf_asReal(cdid_ia_get_elt(cf_g, "alpha"));

    int p_or = (has_or && bid >= 1 && bid <= n_or_blocks) ? p_or_by_bid[bid] : 0;
    const int *kept_or_F = (has_or && bid >= 1 && bid <= n_or_blocks) ? kept_or_by_bid[bid] : NULL;
    const double *beta_or = (has_or && bid >= 1 && bid <= n_or_blocks) ? beta_or_by_bid[bid] : NULL;
    double alpha_bf = (has_or && intercept) ? Rf_asReal(cdid_ia_get_elt(VECTOR_ELT(or_blocks_s, bid - 1), "alpha")) : 0.0;

    memset(w_plus,  0, (size_t) B * sizeof(double));
    memset(j_g_vec, 0, (size_t) B * sizeof(double));
    memset(m_g_vec, 0, (size_t) B * sizeof(double));
    for (int r = 0; r < B; ++r) {
      if (ci[r] != ci_g) continue;
      int row0 = rows[r] - 1;
      w_plus[r] = w[r] * a_g;
      double jv_sl = 0.0;
      for (int m = 0; m < pg_g; ++m) {
        jv_sl += (L_lk * X_F[row0 + (size_t) (kept_g_F[m] - 1) * (size_t) N]) * beta_g[m];
      }
      j_g_vec[r] = alpha_cg + jv_sl;
      if (has_or) {
        double mv_sl = 0.0;
        for (int m = 0; m < p_or; ++m) {
          mv_sl += (L_lk * X_F[row0 + (size_t) (kept_or_F[m] - 1) * (size_t) N]) * beta_or[m];
        }
        m_g_vec[r] = alpha_bf + mv_sl;
      }
    }

    if (intercept && alpha_if[t_cid] != NULL) {
      memcpy(j_T_vec, alpha_if[t_cid], (size_t) B * sizeof(double));
    } else {
      memset(j_T_vec, 0, (size_t) B * sizeof(double));
    }
    if (pg_g > 0 && xi_beta[gid_g] != NULL) {
      const double *dbpj = REAL(cdid_ia_get_elt(lf, "d_bar_plus_j"));
      const double *xb = xi_beta[gid_g];
      for (int m = 0; m < pg_g; ++m) {
        double v = dbpj[m];
        if (v != 0.0) {
          const double *xcol = xb + (size_t) m * (size_t) B;
          for (int r = 0; r < B; ++r) j_T_vec[r] += xcol[r] * v;
        }
      }
    }

    int tid = has_trt ? Rf_asInteger(cdid_ia_get_elt(lf, "trt_id")) : 0;
    double d_minus = has_trt ? Rf_asReal(cdid_ia_get_elt(lf, "d_minus")) : 1.0;
    double C_ol = has_trt ? Rf_asReal(cdid_ia_get_elt(lf, "C_ol")) : 0.0;
    double B_ol = (has_trt && has_or) ? Rf_asReal(cdid_ia_get_elt(lf, "B_ol")) : 0.0;
    double B_gl = has_or ? Rf_asReal(cdid_ia_get_elt(lf, "B_gl")) : 0.0;
    double a_coeff_Bg_g = has_or ? Rf_asReal(cdid_ia_get_elt(lf, "a_coeff_Bg_g")) : 0.0;
    const double *d_bar_plus_or = (has_or && p_or > 0) ? REAL(cdid_ia_get_elt(lf, "d_bar_plus_or")) : NULL;
    const double *d_bar_minus_or = (has_trt && has_or && p_or > 0) ? REAL(cdid_ia_get_elt(lf, "d_bar_minus_or")) : NULL;

    if (has_trt) {
      memset(w_minus,  0, (size_t) B * sizeof(double));
      memset(j_k_vec,  0, (size_t) B * sizeof(double));
      memset(m_k_vec,  0, (size_t) B * sizeof(double));
      memset(a_C_vec,  0, (size_t) B * sizeof(double));
      memset(j_C_vec,  0, (size_t) B * sizeof(double));
      memset(a_Bo_vec, 0, (size_t) B * sizeof(double));

      SEXP ctrl_cohs_s = cdid_ia_get_elt(lf, "control_cohort_indices");
      SEXP ctrl_cids_s = cdid_ia_get_elt(lf, "control_cell_ids");
      int n_ctrl = (int) XLENGTH(ctrl_cohs_s);
      const int *ctrl_cohs = INTEGER(ctrl_cohs_s);
      const int *ctrl_cids = INTEGER(ctrl_cids_s);
      const double *a_coeff_C_k = REAL(cdid_ia_get_elt(lf, "a_coeff_C_k"));
      const double *m_minus_k = REAL(cdid_ia_get_elt(lf, "m_minus_k"));
      SEXP dbmj_s = cdid_ia_get_elt(lf, "d_bar_minus_j_k");
      const double *a_coeff_Bo_k = has_or ? REAL(cdid_ia_get_elt(lf, "a_coeff_Bo_k")) : NULL;
      const double *odds_t = trt_odds[tid];

      memset(ctrl_coh_map, 0, (size_t) (G + 1) * sizeof(int));
      for (int idx_k = 0; idx_k < n_ctrl; ++idx_k) {
        int ci_k = ctrl_cohs[idx_k];
        if (ci_k >= 1 && ci_k <= G) ctrl_coh_map[ci_k] = idx_k + 1;
      }

      for (int r = 0; r < B; ++r) {
        int c = ci[r];
        int idx1 = (c >= 1 && c <= G) ? ctrl_coh_map[c] : 0;
        if (idx1 <= 0) continue;
        int idx_k = idx1 - 1;
        int cid_k = ctrl_cids[idx_k];
        SEXP cf_k = VECTOR_ELT(cells_s, cid_k - 1);
        int gid_k = Rf_asInteger(cdid_ia_get_elt(cf_k, "group_id"));
        SEXP gf_k = VECTOR_ELT(groups_s, gid_k - 1);
        int pg_k = Rf_asInteger(cdid_ia_get_elt(gf_k, "p"));
        const int *kept_k_F = (pg_k > 0) ? INTEGER(cdid_ia_get_elt(gf_k, "kept_cols_F")) : NULL;
        const double *beta_k = (pg_k > 0) ? REAL(cdid_ia_get_elt(gf_k, "beta")) : NULL;
        double alpha_ck = Rf_asReal(cdid_ia_get_elt(cf_k, "alpha"));

        int row0 = rows[r] - 1;
        w_minus[r] = (w[r] * a_vec[c - 1]) * odds_t[r];
        double jv_sl = 0.0;
        for (int m = 0; m < pg_k; ++m) {
          jv_sl += (L_lk * X_F[row0 + (size_t) (kept_k_F[m] - 1) * (size_t) N]) * beta_k[m];
        }
        j_k_vec[r] = alpha_ck + jv_sl;
        if (has_or) {
          double mv_sl = 0.0;
          for (int m = 0; m < p_or; ++m) {
            mv_sl += (L_lk * X_F[row0 + (size_t) (kept_or_F[m] - 1) * (size_t) N]) * beta_or[m];
          }
          m_k_vec[r] = alpha_bf + mv_sl;
        }
      }

      for (int idx_k = 0; idx_k < n_ctrl; ++idx_k) {
        int ci_k = ctrl_cohs[idx_k];
        int cid_k = ctrl_cids[idx_k];
        SEXP cf_k = VECTOR_ELT(cells_s, cid_k - 1);
        int gid_k = Rf_asInteger(cdid_ia_get_elt(cf_k, "group_id"));
        int pg_k = p_group[gid_k];
        const double *xi_ak = xi_a + (size_t) (ci_k - 1) * (size_t) B;
        double ack = a_coeff_C_k[idx_k];
        if (ack != 0.0) {
          for (int r = 0; r < B; ++r) a_C_vec[r] += ack * xi_ak[r];
        }
        if (has_or) {
          double abok = a_coeff_Bo_k[idx_k];
          if (abok != 0.0) {
            for (int r = 0; r < B; ++r) a_Bo_vec[r] += abok * xi_ak[r];
          }
        }
        if (intercept && alpha_if[cid_k] != NULL) {
          double mmk = m_minus_k[idx_k];
          const double *aic = alpha_if[cid_k];
          if (mmk != 0.0) {
            for (int r = 0; r < B; ++r) j_C_vec[r] += mmk * aic[r];
          }
        }
        if (pg_k > 0 && xi_beta[gid_k] != NULL) {
          const double *dbm = REAL(VECTOR_ELT(dbmj_s, idx_k));
          const double *xb = xi_beta[gid_k];
          for (int m = 0; m < pg_k; ++m) {
            double v = dbm[m];
            if (v != 0.0) {
              const double *xcol = xb + (size_t) m * (size_t) B;
              for (int r = 0; r < B; ++r) j_C_vec[r] += xcol[r] * v;
            }
          }
        }
      }
    }

    const double *xi_ag = xi_a + (size_t) (ci_g - 1) * (size_t) B;
    const double *g_C_gam  = (has_trt && p_trt > 0) ? REAL(cdid_ia_get_elt(lf, "g_C_gamma")) : NULL;
    const double *g_Bo_gam = (has_trt && has_or && p_trt > 0) ? REAL(cdid_ia_get_elt(lf, "g_Bo_gamma")) : NULL;
    const double *t_dir = (has_trt && p_trt > 0) ? trt_dir[tid] : NULL;
    const double *t_a   = (has_trt && p_trt > 0) ? trt_a[tid] : NULL;
    const double *oa_dir = (has_or && intercept) ? or_alpha_dir[bid] : NULL;
    const double *oa_a   = (has_or && intercept) ? or_alpha_a[bid] : NULL;
    const double *oa_j   = (has_or && intercept) ? or_alpha_j[bid] : NULL;
    const double *ob_dir = (has_or && p_or > 0) ? or_beta_dir[bid] : NULL;
    const double *ob_a   = (has_or && p_or > 0) ? or_beta_a[bid] : NULL;
    const double *ob_j   = (has_or && p_or > 0) ? or_beta_j[bid] : NULL;

    double *phi_col = phi + (size_t) j * (size_t) B;
    for (int r = 0; r < B; ++r) {
      double direct_T = w_plus[r] * (j_g_vec[r] - T_gl) / d_plus;
      double a_T = a_coeff_T_g * xi_ag[r];
      double j_T = j_T_vec[r];

      double direct_Bg = 0.0, a_Bg = 0.0, or_Bg_dir = 0.0, or_Bg_a = 0.0, or_Bg_j = 0.0;
      if (has_or) {
        direct_Bg = w_plus[r] * (m_g_vec[r] - B_gl) / d_plus;
        a_Bg = a_coeff_Bg_g * xi_ag[r];
        for (int m = 0; m < p_or; ++m) {
          double dbp = d_bar_plus_or[m];
          or_Bg_dir = fma(ob_dir[r + (size_t) m * (size_t) B], dbp, or_Bg_dir);
          or_Bg_a   = fma(ob_a[r   + (size_t) m * (size_t) B], dbp, or_Bg_a);
          or_Bg_j   = fma(ob_j[r   + (size_t) m * (size_t) B], dbp, or_Bg_j);
        }
        if (intercept) {
          or_Bg_dir = oa_dir[r] + or_Bg_dir;
          or_Bg_a   = oa_a[r]   + or_Bg_a;
          or_Bg_j   = oa_j[r]   + or_Bg_j;
        }
      }

      double direct_C = 0.0, a_C = 0.0, j_C = 0.0, trt_C_dir = 0.0, trt_C_a = 0.0;
      double direct_Bo = 0.0, a_Bo = 0.0, trt_Bo_dir = 0.0, trt_Bo_a = 0.0;
      double or_Bo_dir = 0.0, or_Bo_a = 0.0, or_Bo_j = 0.0;
      if (has_trt) {
        direct_C = w_minus[r] * (j_k_vec[r] - C_ol) / d_minus;
        a_C = a_C_vec[r];
        j_C = j_C_vec[r];
        for (int m = 0; m < p_trt; ++m) {
          double gcg = g_C_gam[m];
          trt_C_dir = fma(t_dir[r + (size_t) m * (size_t) B], gcg, trt_C_dir);
          trt_C_a   = fma(t_a[r   + (size_t) m * (size_t) B], gcg, trt_C_a);
        }
        if (has_or) {
          direct_Bo = w_minus[r] * (m_k_vec[r] - B_ol) / d_minus;
          a_Bo = a_Bo_vec[r];
          for (int m = 0; m < p_trt; ++m) {
            double gbog = g_Bo_gam[m];
            trt_Bo_dir = fma(t_dir[r + (size_t) m * (size_t) B], gbog, trt_Bo_dir);
            trt_Bo_a   = fma(t_a[r   + (size_t) m * (size_t) B], gbog, trt_Bo_a);
          }
          for (int m = 0; m < p_or; ++m) {
            double dbm = d_bar_minus_or[m];
            or_Bo_dir = fma(ob_dir[r + (size_t) m * (size_t) B], dbm, or_Bo_dir);
            or_Bo_a   = fma(ob_a[r   + (size_t) m * (size_t) B], dbm, or_Bo_a);
            or_Bo_j   = fma(ob_j[r   + (size_t) m * (size_t) B], dbm, or_Bo_j);
          }
          if (intercept) {
            or_Bo_dir = oa_dir[r] + or_Bo_dir;
            or_Bo_a   = oa_a[r]   + or_Bo_a;
            or_Bo_j   = oa_j[r]   + or_Bo_j;
          }
        }
      }

      double val = 0.0;
      if (mode_est == 0) {
        val = (direct_T - direct_C) + (a_T - a_C) + (j_T - j_C) + (-trt_C_dir) + (-trt_C_a);
      } else if (mode_est == 1) {
        val = (direct_T - direct_Bg) + (a_T - a_Bg) + j_T + (-or_Bg_dir) + (-or_Bg_a) + (-or_Bg_j);
      } else {
        val = ((direct_T - direct_Bg) - (direct_C - direct_Bo)) +
              ((a_T - a_Bg) - (a_C - a_Bo)) +
              (j_T - j_C) +
              (-(trt_C_dir - trt_Bo_dir)) +
              (-(trt_C_a - trt_Bo_a)) +
              (-or_Bg_dir + or_Bo_dir) +
              (-or_Bg_a + or_Bo_a) +
              (-or_Bg_j + or_Bo_j);
      }
      if (!R_FINITE(val)) {
        UNPROTECT(2);
        vmaxset(vmax);
        return R_NilValue;
      }
      phi_col[r] = is_core_scale ? (lambda * val) : val;
    }
  }

  SEXP dn_s = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn_s, 0, R_NilValue);
  SET_VECTOR_ELT(dn_s, 1, link_names_s);
  Rf_setAttrib(phi_s, R_DimNamesSymbol, dn_s);

  UNPROTECT(3);
  vmaxset(vmax);
  return phi_s;
}

/* Helper: evaluate fitted selection probability for one selection cell_id at 0-based frame row0 */
static int cdid_ia_eval_sel_cell_prob(
    SEXP sel_cells_s,
    SEXP sel_blocks_s,
    const double *Z_sel,
    int N,
    int d_full,
    const int *slope_cols,
    int cell_id,
    int row0,
    double *prob_out) {
  if (cell_id < 1 || cell_id > (int) XLENGTH(sel_cells_s)) return 0;
  SEXP cell_s = VECTOR_ELT(sel_cells_s, cell_id - 1);
  int block_id = Rf_asInteger(cdid_ia_get_elt(cell_s, "block_id"));
  int local_b = Rf_asInteger(cdid_ia_get_elt(cell_s, "local_b"));
  if (block_id < 1 || block_id > (int) XLENGTH(sel_blocks_s)) return 0;
  SEXP bl_s = VECTOR_ELT(sel_blocks_s, block_id - 1);
  const char *branch = CHAR(STRING_ELT(cdid_ia_get_elt(bl_s, "branch"), 0));
  if (strcmp(branch, "zero_col") == 0) {
    *prob_out = 0.5;
    return 1;
  }
  double eta = 0.0;
  if (strcmp(branch, "cp04_separate_pair") == 0) {
    const double *coef = REAL(cdid_ia_get_elt(bl_s, "coefficients"));
    for (int col = 0; col < d_full; ++col) {
      eta += Z_sel[row0 + (size_t) col * (size_t) N] * coef[col];
    }
  } else {
    int has_int = (Rf_asLogical(cdid_ia_get_elt(bl_s, "has_intercept")) == 1);
    int d_s = Rf_asInteger(cdid_ia_get_elt(bl_s, "d_s"));
    double al = has_int ? REAL(cdid_ia_get_elt(bl_s, "alpha"))[local_b - 1] : 0.0;
    double sl = 0.0;
    if (d_s > 0) {
      const double *beta = REAL(cdid_ia_get_elt(bl_s, "beta"));
      for (int j = 0; j < d_s; ++j) {
        int col0 = slope_cols[j] - 1;
        sl += Z_sel[row0 + (size_t) col0 * (size_t) N] * beta[j];
      }
    }
    eta = al + sl;
  }
  if (!R_FINITE(eta)) return 0;
  double pr = 0.0, od = 0.0;
  cdid_ia_clip_logit(eta, &pr, &od);
  if (!R_FINITE(pr) || pr <= 0.0 || pr >= 1.0) return 0;
  *prob_out = pr;
  return 1;
}

static int cdid_ia_eval_sel_pair_prob(
    SEXP sel_pcs_s,
    SEXP sel_cells_s,
    SEXP sel_blocks_s,
    int is_indep_period,
    const double *Z_sel,
    int N,
    int d_full,
    const int *slope_cols,
    int q_id,
    int row0,
    double *prob_out) {
  if (q_id < 1 || q_id > (int) XLENGTH(sel_pcs_s)) return 0;
  SEXP pc_s = VECTOR_ELT(sel_pcs_s, q_id - 1);
  if (is_indep_period) {
    int s_cid = Rf_asInteger(cdid_ia_get_elt(pc_s, "s_cell_id"));
    int t_cid = Rf_asInteger(cdid_ia_get_elt(pc_s, "t_cell_id"));
    double ps = 0.0, pt = 0.0;
    if (!cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N, d_full, slope_cols, s_cid, row0, &ps) ||
        !cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N, d_full, slope_cols, t_cid, row0, &pt)) {
      return 0;
    }
    double pp = ps * pt;
    if (!R_FINITE(pp) || pp <= 0.0 || pp >= 1.0) return 0;
    *prob_out = pp;
    return 1;
  } else {
    int p_cid = Rf_asInteger(cdid_ia_get_elt(pc_s, "pair_cell_id"));
    return cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N, d_full, slope_cols, p_cid, row0, prob_out);
  }
}

/* 4. Fast C Builder for augmented_cell_replay_store */
SEXP C_cdid_ia_cell_replay_build(
    SEXP backend_s,
    SEXP st_s,
    SEXP panel_s,
    SEXP X_F_s,
    SEXP ids_s,
    SEXP cids_s,
    SEXP zids_s,
    SEXP B_build_s) {
  if (TYPEOF(backend_s) != VECSXP || TYPEOF(st_s) != VECSXP ||
      TYPEOF(panel_s) != VECSXP || TYPEOF(X_F_s) != REALSXP ||
      TYPEOF(ids_s) != INTSXP || TYPEOF(cids_s) != INTSXP ||
      TYPEOF(zids_s) != INTSXP) {
    return R_NilValue;
  }
  const void *vmax = vmaxget();

  int N = Rf_asInteger(cdid_ia_get_elt(st_s, "N"));
  const int *cohort_order = INTEGER(cdid_ia_get_elt(st_s, "cohort_order"));
  SEXP c_start_s          = cdid_ia_get_elt(st_s, "cohort_start");
  SEXP c_end_s            = cdid_ia_get_elt(st_s, "cohort_end");
  const int *a_index      = INTEGER(cdid_ia_get_elt(st_s, "a_index"));
  const double *w_all     = REAL(cdid_ia_get_elt(st_s, "w"));
  SEXP Z_sel_s            = cdid_ia_get_elt(st_s, "Z");
  const double *Z_sel     = REAL(Z_sel_s);
  int G = (int) XLENGTH(cdid_ia_get_elt(st_s, "cohorts"));
  int B_build = Rf_asInteger(B_build_s);
  if (B_build < 1) B_build = N;

  const double *X_F = REAL(X_F_s);
  SEXP all_times_s = cdid_ia_get_elt(panel_s, "all_times");
  int n_dates = (int) XLENGTH(all_times_s);
  const double *all_times = REAL(all_times_s);
  SEXP post_s = cdid_ia_get_elt(panel_s, "postings");
  const int *post_a = INTEGER(cdid_ia_get_elt(post_s, "post_a_index"));
  const int *post_obs = INTEGER(cdid_ia_get_elt(post_s, "post_obs_index"));
  const int *date_start = INTEGER(cdid_ia_get_elt(post_s, "date_start"));
  const int *date_end = INTEGER(cdid_ia_get_elt(post_s, "date_end"));
  const int *date_count = INTEGER(cdid_ia_get_elt(post_s, "date_count"));
  SEXP y_obs_s = cdid_ia_get_elt(cdid_ia_get_elt(panel_s, "observations"), "y");
  const double *y_obs = REAL(y_obs_s);
  R_xlen_t n_y = XLENGTH(y_obs_s);
  R_xlen_t n_post = XLENGTH(cdid_ia_get_elt(post_s, "post_a_index"));

  SEXP cells_s = cdid_ia_get_elt(backend_s, "cells");
  SEXP groups_s = cdid_ia_get_elt(backend_s, "groups");
  SEXP sel_fit_s = cdid_ia_get_elt(backend_s, "selection_backend");
  SEXP sel_cells_s = cdid_ia_get_elt(sel_fit_s, "cells");
  SEXP sel_blocks_s = cdid_ia_get_elt(sel_fit_s, "blocks");
  SEXP sel_pcs_s = cdid_ia_get_elt(sel_fit_s, "pair_cohorts");
  SEXP dspec_s = cdid_ia_get_elt(sel_fit_s, "design_spec");
  int d_full = Rf_asInteger(cdid_ia_get_elt(dspec_s, "d"));
  int d_s = Rf_asInteger(cdid_ia_get_elt(dspec_s, "d_s"));
  const int *slope_cols = (d_s > 0) ? INTEGER(cdid_ia_get_elt(dspec_s, "slope_cols")) : NULL;
  int is_indep_period = (strcmp(CHAR(STRING_ELT(cdid_ia_get_elt(sel_fit_s, "selection_model"), 0)), "independent_period_logit") == 0);

  int nc = (int) XLENGTH(cells_s);
  int n_ids = (int) XLENGTH(ids_s);
  int n_cids = (int) XLENGTH(cids_s);
  int n_zids = (int) XLENGTH(zids_s);
  const int *ids = INTEGER(ids_s);
  const int *cids = INTEGER(cids_s);
  const int *zids = INTEGER(zids_s);

  int *need_score_mask = (int *) R_alloc((size_t) (nc + 1), sizeof(int));
  int *need_Z_mask     = (int *) R_alloc((size_t) (nc + 1), sizeof(int));
  int *coh_seen        = (int *) R_alloc((size_t) (G + 1), sizeof(int));
  memset(need_score_mask, 0, (size_t) (nc + 1) * sizeof(int));
  memset(need_Z_mask,     0, (size_t) (nc + 1) * sizeof(int));
  memset(coh_seen,        0, (size_t) (G + 1) * sizeof(int));

  for (int i = 0; i < n_cids; ++i) {
    int c = cids[i];
    if (c >= 1 && c <= nc) need_score_mask[c] = 1;
  }
  for (int i = 0; i < n_zids; ++i) {
    int c = zids[i];
    if (c >= 1 && c <= nc) need_Z_mask[c] = 1;
  }

  SEXP records_s = PROTECT(Rf_allocVector(VECSXP, nc));
  SEXP rank_s = PROTECT(Rf_allocVector(INTSXP, N));
  int *rank_vec = INTEGER(rank_s);
  memset(rank_vec, 0, (size_t) N * sizeof(int));

  int batches = 0;
  double d_buf[256];

  for (int idx = 0; idx < n_ids; ++idx) {
    int cid = ids[idx];
    if (cid < 1 || cid > nc) { UNPROTECT(2); vmaxset(vmax); return R_NilValue; }
    SEXP cf = VECTOR_ELT(cells_s, cid - 1);
    int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
    int cohort = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
    int np = Rf_asInteger(cdid_ia_get_elt(cf, "n_pair"));
    int q_id = Rf_asInteger(cdid_ia_get_elt(cf, "q_id"));
    double L_val = Rf_asReal(cdid_ia_get_elt(cf, "L"));
    double alpha_c = Rf_asReal(cdid_ia_get_elt(cf, "alpha"));

    SEXP gf = VECTOR_ELT(groups_s, gid - 1);
    int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
    if (pg > 256) { UNPROTECT(2); vmaxset(vmax); return R_NilValue; }
    const int *kept_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf, "kept_cols_F")) : NULL;
    const double *beta_g = (pg > 0) ? REAL(cdid_ia_get_elt(gf, "beta")) : NULL;

    int c_start = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[cohort - 1] : (int) REAL(c_start_s)[cohort - 1];
    int c_end   = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[cohort - 1]   : (int) REAL(c_end_s)[cohort - 1];
    int n_c     = c_end - c_start + 1;
    if (!coh_seen[cohort]) {
      coh_seen[cohort] = 1;
      for (int k = 0; k < n_c; ++k) {
        int r1 = cohort_order[c_start - 1 + k];
        if (r1 >= 1 && r1 <= N) rank_vec[r1 - 1] = k + 1;
      }
    }

    int need_score = need_score_mask[cid];
    int need_Z     = need_Z_mask[cid];

    int ds = 0, dt = 0;
    if (!cdid_ia_resolve_pair_dates(cdid_ia_get_elt(cf, "pair"), all_times, n_dates, &ds, &dt)) {
      UNPROTECT(2);
      vmaxset(vmax);
      return R_NilValue;
    }
    int has_pair = (date_count[ds - 1] > 0 && date_count[dt - 1] > 0);
    int s_lo = date_start[ds - 1] - 1, s_hi = date_end[ds - 1] - 1;
    int t_lo = date_start[dt - 1] - 1, t_hi = date_end[dt - 1] - 1;

    SEXP rec_s = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP rec_nms = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(rec_nms, 0, Rf_mkChar("rows"));
    SET_STRING_ELT(rec_nms, 1, Rf_mkChar("d"));
    SET_STRING_ELT(rec_nms, 2, Rf_mkChar("score"));
    SET_STRING_ELT(rec_nms, 3, Rf_mkChar("Z"));
    SET_STRING_ELT(rec_nms, 4, Rf_mkChar("cohort"));
    Rf_setAttrib(rec_s, R_NamesSymbol, rec_nms);

    SEXP sr_s = need_score ? PROTECT(Rf_allocVector(INTSXP, np)) : R_NilValue;
    SEXP d_s_mat = need_score ? PROTECT(Rf_allocMatrix(REALSXP, np, pg)) : R_NilValue;
    SEXP sc_s = need_score ? PROTECT(Rf_allocVector(REALSXP, np)) : R_NilValue;
    SEXP z_s = need_Z ? PROTECT(Rf_allocVector(REALSXP, n_c)) : R_NilValue;
    int nprot_rec = 2 + (need_score ? 3 : 0) + (need_Z ? 1 : 0);

    if (need_score && pg > 0) {
      SEXP kc_s = cdid_ia_get_elt(gf, "kept_columns");
      if (TYPEOF(kc_s) == STRSXP) {
        SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dn, 0, R_NilValue);
        SET_VECTOR_ELT(dn, 1, kc_s);
        Rf_setAttrib(d_s_mat, R_DimNamesSymbol, dn);
        UNPROTECT(1);
      }
    }

    int *score_rows = need_score ? INTEGER(sr_s) : NULL;
    double *d_mat   = need_score ? REAL(d_s_mat) : NULL;
    double *score   = need_score ? REAL(sc_s)    : NULL;
    double *z_vec   = need_Z     ? REAL(z_s)     : NULL;

    int filled = 0;
    int hs_hint = s_lo, ht_hint = t_lo;
    batches += (n_c + B_build - 1) / B_build;

    for (int k = 0; k < n_c; ++k) {
      int r1 = cohort_order[c_start - 1 + k];
      int row0 = r1 - 1;
      double jhat_sl = 0.0;
      for (int m = 0; m < pg; ++m) {
        double d_m = L_val * X_F[row0 + (size_t) (kept_F[m] - 1) * (size_t) N];
        d_buf[m] = d_m;
        jhat_sl += d_m * beta_g[m];
      }
      double jhat = alpha_c + jhat_sl;
      int a = a_index[row0];
      int is_obs = 0;
      double dy = 0.0;
      if (has_pair && a > 0) {
        int hs = cdid_ia_find_posting(post_a, s_lo, s_hi, a, &hs_hint);
        if (hs >= 0) {
          int ht = cdid_ia_find_posting(post_a, t_lo, t_hi, a, &ht_hint);
          if (ht >= 0) {
            int os = post_obs[hs], ot = post_obs[ht];
            if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y) {
              UNPROTECT(nprot_rec + 2);
              vmaxset(vmax);
              return R_NilValue;
            }
            dy = y_obs[ot - 1] - y_obs[os - 1];
            if (!R_FINITE(dy)) {
              UNPROTECT(nprot_rec + 2);
              vmaxset(vmax);
              return R_NilValue;
            }
            is_obs = 1;
          }
        }
      }
      if (need_Z) {
        double zv = jhat;
        if (is_obs) {
          double q = 0.0;
          if (!cdid_ia_eval_sel_pair_prob(sel_pcs_s, sel_cells_s, sel_blocks_s, is_indep_period,
                                          Z_sel, N, d_full, slope_cols, q_id, row0, &q)) {
            UNPROTECT(nprot_rec + 2);
            vmaxset(vmax);
            return R_NilValue;
          }
          double q_fl = (q < 1e-6) ? 1e-6 : q;
          double pi_r = 1.0 / q_fl;
          zv = jhat + pi_r * (dy - jhat);
        }
        z_vec[k] = zv;
      }
      if (need_score && is_obs) {
        if (filled >= np) {
          UNPROTECT(nprot_rec + 2);
          vmaxset(vmax);
          return R_NilValue;
        }
        score_rows[filled] = r1;
        for (int m = 0; m < pg; ++m) {
          d_mat[filled + (size_t) m * (size_t) np] = d_buf[m];
        }
        score[filled] = w_all[row0] * (dy - jhat);
        ++filled;
      }
    }

    if (need_score && filled != np) {
      UNPROTECT(nprot_rec + 2);
      vmaxset(vmax);
      return R_NilValue;
    }

    SET_VECTOR_ELT(rec_s, 0, sr_s);
    SET_VECTOR_ELT(rec_s, 1, d_s_mat);
    SET_VECTOR_ELT(rec_s, 2, sc_s);
    SET_VECTOR_ELT(rec_s, 3, z_s);
    SET_VECTOR_ELT(rec_s, 4, Rf_ScalarInteger(cohort));
    SET_VECTOR_ELT(records_s, cid - 1, rec_s);
    UNPROTECT(nprot_rec);
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP out_nms = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(out_nms, 0, Rf_mkChar("records"));
  SET_STRING_ELT(out_nms, 1, Rf_mkChar("rank"));
  SET_STRING_ELT(out_nms, 2, Rf_mkChar("batches"));
  SET_VECTOR_ELT(out_s, 0, records_s);
  SET_VECTOR_ELT(out_s, 1, rank_s);
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarInteger(batches));
  Rf_setAttrib(out_s, R_NamesSymbol, out_nms);

  UNPROTECT(4);
  vmaxset(vmax);
  return out_s;
}

/* Helper: contract an adjoint (t_lam_by_cid, t_nu_by_gid) with chunk observed cell scores */
static void cdid_ia_contract_adjoint_chunk(
    const double *t_lam_by_cid,
    double * const *t_nu_by_gid,
    const int *cids,
    int n_cids,
    const int *cell_gid,
    const int *p_group,
    int intercept,
    const int *cell_n_obs,
    int * const *cell_obs_rr,
    double * const *cell_s_obs,
    double * const *cell_dt_obs,
    int B,
    double *out_vec) {
  memset(out_vec, 0, (size_t) B * sizeof(double));
  for (int ic = 0; ic < n_cids; ++ic) {
    int cid = cids[ic];
    int n_ob = cell_n_obs[cid];
    if (n_ob <= 0) continue;
    int gid = cell_gid[cid];
    int pg = p_group[gid];
    double t_lam = (intercept && t_lam_by_cid != NULL) ? t_lam_by_cid[cid] : 0.0;
    const double *t_nu = (pg > 0 && t_nu_by_gid != NULL) ? t_nu_by_gid[gid] : NULL;
    const int *rr = cell_obs_rr[cid];
    const double *s_rr = cell_s_obs[cid];
    const double *dt = cell_dt_obs[cid];
    for (int k = 0; k < n_ob; ++k) {
      double val = t_lam;
      if (pg > 0 && t_nu != NULL) {
        double dot = 0.0;
        for (int m = 0; m < pg; ++m) {
          dot += dt[k + (size_t) m * (size_t) n_ob] * t_nu[m];
        }
        val += dot;
      }
      out_vec[rr[k]] += s_rr[k] * val;
    }
  }
}

/* Helper: unpack SEXP adjoint (dense or compact_reconstruction_adjoint) and contract on chunk */
static int cdid_ia_eval_sexp_adjoint_chunk(
    SEXP adj_s,
    int n_cells,
    int n_groups,
    const int *cids,
    int n_cids,
    const int *cell_gid,
    const int *p_group,
    int intercept,
    const int *cell_n_obs,
    int * const *cell_obs_rr,
    double * const *cell_s_obs,
    double * const *cell_dt_obs,
    int B,
    double *t_lam_tmp,
    double **t_nu_tmp,
    const double *zero_pg_buf,
    double *out_vec) {
  if (adj_s == R_NilValue || TYPEOF(adj_s) != VECSXP) {
    memset(out_vec, 0, (size_t) B * sizeof(double));
    return 1;
  }
  memset(t_lam_tmp, 0, (size_t) (n_cells + 1) * sizeof(double));
  for (int g = 0; g <= n_groups; ++g) t_nu_tmp[g] = (double *) zero_pg_buf;
  SEXP kind_s = cdid_ia_get_elt(adj_s, "kind");
  int is_compact = (kind_s != R_NilValue && TYPEOF(kind_s) == STRSXP && XLENGTH(kind_s) == 1 &&
                    strcmp(CHAR(STRING_ELT(kind_s, 0)), "compact_reconstruction_adjoint") == 0);
  SEXP tl_s = cdid_ia_get_elt(adj_s, "tilde_lambda");
  SEXP tn_s = cdid_ia_get_elt(adj_s, "tilde_nu");
  if (is_compact) {
    SEXP ac_s = cdid_ia_get_elt(adj_s, "cell_ids");
    SEXP ag_s = cdid_ia_get_elt(adj_s, "group_ids");
    int nac = (int) XLENGTH(ac_s);
    int nag = (int) XLENGTH(ag_s);
    const int *ac = INTEGER(ac_s);
    const int *ag = INTEGER(ag_s);
    const double *tl = REAL(tl_s);
    for (int k = 0; k < nac; ++k) {
      int c = ac[k];
      if (c >= 1 && c <= n_cells) t_lam_tmp[c] = tl[k];
    }
    for (int k = 0; k < nag; ++k) {
      int g = ag[k];
      if (g >= 1 && g <= n_groups) t_nu_tmp[g] = REAL(VECTOR_ELT(tn_s, k));
    }
  } else {
    const double *tl = REAL(tl_s);
    for (int c = 1; c <= n_cells && c <= (int) XLENGTH(tl_s); ++c) {
      t_lam_tmp[c] = tl[c - 1];
    }
    for (int g = 1; g <= n_groups && g <= (int) XLENGTH(tn_s); ++g) {
      t_nu_tmp[g] = REAL(VECTOR_ELT(tn_s, g - 1));
    }
  }
  cdid_ia_contract_adjoint_chunk(
      t_lam_tmp, t_nu_tmp, cids, n_cids, cell_gid, p_group, intercept,
      cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs, B, out_vec);
  return 1;
}

/* Helper: solve scalar group adjoints for (lam_by_cid, nu_by_cid) and contract on chunk */
static int cdid_ia_eval_imp_linear_chunk(
    const double *lam_by_cid,
    double * const *nu_by_cid,
    SEXP cells_s,
    SEXP groups_s,
    int n_cells,
    int n_groups,
    int N,
    const int *cids,
    int n_cids,
    const int *cell_gid,
    const int *p_group,
    int intercept,
    const int *cell_n_obs,
    int * const *cell_obs_rr,
    double * const *cell_s_obs,
    double * const *cell_dt_obs,
    int B,
    double *t_lam_tmp,
    double **t_nu_alloc,
    double *out_vec) {
  memset(t_lam_tmp, 0, (size_t) (n_cells + 1) * sizeof(double));
  double rhs_g[256];
  for (int gid = 1; gid <= n_groups; ++gid) {
    SEXP gf = VECTOR_ELT(groups_s, gid - 1);
    int pg = p_group[gid];
    SEXP g_cids_s = cdid_ia_get_elt(gf, "cell_ids");
    int n_gc = (int) XLENGTH(g_cids_s);
    const int *g_cids = INTEGER(g_cids_s);
    for (int m = 0; m < pg; ++m) rhs_g[m] = 0.0;
    for (int k = 0; k < n_gc; ++k) {
      int cid = g_cids[k];
      SEXP cf = VECTOR_ELT(cells_s, cid - 1);
      double lam = (intercept && lam_by_cid != NULL) ? lam_by_cid[cid] : 0.0;
      if (intercept) t_lam_tmp[cid] = lam / Rf_asReal(cdid_ia_get_elt(cf, "mass"));
      if (pg > 0) {
        const double *nu_c = (nu_by_cid != NULL) ? nu_by_cid[cid] : NULL;
        const double *md_c = intercept ? REAL(cdid_ia_get_elt(cf, "mean_d")) : NULL;
        for (int m = 0; m < pg; ++m) {
          double nc_m = (nu_c != NULL) ? nu_c[m] : 0.0;
          rhs_g[m] += nc_m - (intercept ? (lam * md_c[m]) : 0.0);
        }
      }
    }
    if (pg > 0) {
      const double *R_g = NULL;
      const int *piv_g = NULL;
      if (!cdid_ia_get_wls_fit(cdid_ia_get_elt(gf, "fit"), &R_g, &piv_g, pg) ||
          !cdid_ia_gram_solve_rows(R_g, piv_g, pg, rhs_g, 1, (double) N, t_nu_alloc[gid])) {
        return 0;
      }
    }
  }
  cdid_ia_contract_adjoint_chunk(
      t_lam_tmp, t_nu_alloc, cids, n_cids, cell_gid, p_group, intercept,
      cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs, B, out_vec);
  return 1;
}

/* 5. Fast C Augmented IF Chunk Kernel (non-detailed) */
SEXP C_cdid_ia_aug_if_chunk(
    SEXP backend_s,
    SEXP st_s,
    SEXP panel_s,
    SEXP X_F_s,
    SEXP rows_s,
    SEXP links_s,
    SEXP sets_s,
    SEXP scale_s,
    SEXP up_in_s) {
  if (TYPEOF(backend_s) != VECSXP || TYPEOF(st_s) != VECSXP ||
      TYPEOF(panel_s) != VECSXP || TYPEOF(X_F_s) != REALSXP ||
      TYPEOF(rows_s) != INTSXP || TYPEOF(links_s) != INTSXP ||
      TYPEOF(sets_s) != VECSXP || TYPEOF(scale_s) != STRSXP || XLENGTH(scale_s) != 1) {
    return R_NilValue;
  }
  const void *vmax = vmaxget();

  int B = (int) XLENGTH(rows_s);
  int nl = (int) XLENGTH(links_s);
  const int *rows = INTEGER(rows_s);
  const int *links = INTEGER(links_s);
  int is_core_scale = (strcmp(CHAR(STRING_ELT(scale_s, 0)), "core") == 0);

  int N = Rf_asInteger(cdid_ia_get_elt(st_s, "N"));
  double lambda = Rf_asReal(cdid_ia_get_elt(st_s, "lambda"));
  SEXP cohorts_s = cdid_ia_get_elt(st_s, "cohorts");
  int G = (int) XLENGTH(cohorts_s);
  const int *cohort_index = INTEGER(cdid_ia_get_elt(st_s, "cohort_index"));
  const double *w_all = REAL(cdid_ia_get_elt(st_s, "w"));
  const int *a_index = INTEGER(cdid_ia_get_elt(st_s, "a_index"));
  const double *a_vec = REAL(cdid_ia_get_elt(st_s, "a"));
  const double *mass_F = REAL(cdid_ia_get_elt(st_s, "mass_F"));
  const double *mass_U = REAL(cdid_ia_get_elt(st_s, "mass_U"));
  const int *U_all = LOGICAL(cdid_ia_get_elt(st_s, "U"));
  const double *Z_sel = REAL(cdid_ia_get_elt(st_s, "Z"));

  const double *X_F = REAL(X_F_s);
  SEXP all_times_s = cdid_ia_get_elt(panel_s, "all_times");
  int n_dates = (int) XLENGTH(all_times_s);
  const double *all_times = REAL(all_times_s);
  SEXP post_s = cdid_ia_get_elt(panel_s, "postings");
  const int *post_a = INTEGER(cdid_ia_get_elt(post_s, "post_a_index"));
  const int *post_obs = INTEGER(cdid_ia_get_elt(post_s, "post_obs_index"));
  const int *date_start = INTEGER(cdid_ia_get_elt(post_s, "date_start"));
  const int *date_end = INTEGER(cdid_ia_get_elt(post_s, "date_end"));
  const int *date_count = INTEGER(cdid_ia_get_elt(post_s, "date_count"));
  SEXP y_obs_s = cdid_ia_get_elt(cdid_ia_get_elt(panel_s, "observations"), "y");
  const double *y_obs = REAL(y_obs_s);
  R_xlen_t n_y = XLENGTH(y_obs_s);
  R_xlen_t n_post = XLENGTH(cdid_ia_get_elt(post_s, "post_a_index"));

  SEXP recipe_s = cdid_ia_get_elt(backend_s, "outcome_recipe");
  int intercept = (Rf_asLogical(cdid_ia_get_elt(recipe_s, "intercept")) == 1);
  SEXP trt_recipe_s = cdid_ia_get_elt(backend_s, "treatment_recipe");
  const char *est_str = CHAR(STRING_ELT(cdid_ia_get_elt(backend_s, "est_method"), 0));
  int mode_est = (strcmp(est_str, "ipw") == 0) ? 0 : ((strcmp(est_str, "reg") == 0) ? 1 : 2);

  SEXP cells_s = cdid_ia_get_elt(backend_s, "cells");
  SEXP groups_s = cdid_ia_get_elt(backend_s, "groups");
  SEXP shared_or_s = cdid_ia_get_elt(backend_s, "shared_or");
  SEXP treatments_s = cdid_ia_get_elt(backend_s, "treatments");
  SEXP links_all_s = cdid_ia_get_elt(backend_s, "links");
  SEXP sel_fit_s = cdid_ia_get_elt(backend_s, "selection_backend");
  SEXP up_map_s = cdid_ia_get_elt(backend_s, "up_map");
  SEXP sel_cells_s = cdid_ia_get_elt(sel_fit_s, "cells");
  SEXP sel_blocks_s = cdid_ia_get_elt(sel_fit_s, "blocks");
  SEXP sel_pcs_s = cdid_ia_get_elt(sel_fit_s, "pair_cohorts");
  SEXP dspec_s = cdid_ia_get_elt(sel_fit_s, "design_spec");
  int d_full = Rf_asInteger(cdid_ia_get_elt(dspec_s, "d"));
  int d_s = Rf_asInteger(cdid_ia_get_elt(dspec_s, "d_s"));
  if (d_full > 256 || d_s > 256) { vmaxset(vmax); return R_NilValue; }
  const int *slope_cols = (d_s > 0) ? INTEGER(cdid_ia_get_elt(dspec_s, "slope_cols")) : NULL;
  int is_indep_period = (strcmp(CHAR(STRING_ELT(cdid_ia_get_elt(sel_fit_s, "selection_model"), 0)), "independent_period_logit") == 0);

  int n_cells = (int) XLENGTH(cells_s);
  int n_groups = (int) XLENGTH(groups_s);

  int *ci = (int *) R_alloc((size_t) B, sizeof(int));
  double *w = (double *) R_alloc((size_t) B, sizeof(double));
  int *ai = (int *) R_alloc((size_t) B, sizeof(int));
  int *U_rows = (int *) R_alloc((size_t) B, sizeof(int));
  for (int r = 0; r < B; ++r) {
    int row0 = rows[r] - 1;
    ci[r] = cohort_index[row0];
    w[r] = w_all[row0];
    ai[r] = a_index[row0];
    U_rows[r] = (U_all[row0] == 1);
  }

  /* 1. Repeat-rate IFs xi_a (B x G) */
  double *xi_a = (double *) R_alloc((size_t) B * (size_t) G, sizeof(double));
  memset(xi_a, 0, (size_t) B * (size_t) G * sizeof(double));
  for (int r = 0; r < B; ++r) {
    int k = ci[r];
    if (k >= 1 && k <= G) {
      double mf = mass_F[k - 1];
      if (mf > 0.0) {
        double A_r = (ai[r] > 0) ? 1.0 : 0.0;
        xi_a[r + (size_t) (k - 1) * (size_t) B] = w[r] * (A_r - a_vec[k - 1]) / mf;
      }
    }
  }

  SEXP cids_s = cdid_ia_get_elt(sets_s, "cids");
  SEXP zids_s = cdid_ia_get_elt(sets_s, "zids");
  SEXP gids_s = cdid_ia_get_elt(sets_s, "gids");
  SEXP bids_s = cdid_ia_get_elt(sets_s, "bids");
  SEXP tids_s = cdid_ia_get_elt(sets_s, "tids");
  SEXP uids_s = cdid_ia_get_elt(sets_s, "uids");
  int n_cids = (int) XLENGTH(cids_s);
  int n_zids = (int) XLENGTH(zids_s);
  int n_gids = (int) XLENGTH(gids_s);
  int n_bids = (int) XLENGTH(bids_s);
  int n_tids = (int) XLENGTH(tids_s);
  int n_uids = (int) XLENGTH(uids_s);
  const int *cids = INTEGER(cids_s);
  const int *zids = INTEGER(zids_s);
  const int *gids = INTEGER(gids_s);
  const int *bids = INTEGER(bids_s);
  const int *tids = INTEGER(tids_s);
  const int *uids = INTEGER(uids_s);

  /* 2. Upstream selection IFs up (B x n_uids) */
  const double *up = NULL;
  int max_uid = 0;
  for (int j = 0; j < n_uids; ++j) {
    if (uids[j] > max_uid) max_uid = uids[j];
  }
  int n_up_tot = Rf_asInteger(cdid_ia_get_elt(up_map_s, "n_total"));
  if (n_up_tot > max_uid) max_uid = n_up_tot;
  int *uid_pos = (int *) R_alloc((size_t) (max_uid + 2), sizeof(int));
  memset(uid_pos, 0, (size_t) (max_uid + 2) * sizeof(int));
  for (int j = 0; j < n_uids; ++j) {
    int u = uids[j];
    if (u >= 1 && u <= max_uid) uid_pos[u] = j + 1;
  }

  if (n_uids > 0) {
    if (up_in_s != R_NilValue) {
      if (TYPEOF(up_in_s) != REALSXP || XLENGTH(up_in_s) != (R_xlen_t) B * (R_xlen_t) n_uids) {
        vmaxset(vmax);
        return R_NilValue;
      }
      up = REAL(up_in_s);
    } else {
      int n_sing = 0, n_panel_units = 0;
      SEXP id_s_st = R_NilValue, obs_start_s_u = R_NilValue;
      const double *ukey_p = NULL;
      const int *u_nobs_p = NULL, *o_cal_p = NULL;
      if (is_indep_period) {
        SEXP ns_s = cdid_ia_get_elt(panel_s, "n_singleton_units");
        n_sing = (ns_s != R_NilValue) ? Rf_asInteger(ns_s) : 0;
        if (n_sing > 0) {
          SEXP id_type_s = cdid_ia_get_elt(panel_s, "id_type");
          SEXP ukey_s = cdid_ia_get_elt(panel_s, "unit_id_key");
          SEXP units_s = cdid_ia_get_elt(panel_s, "units");
          SEXP cal_s = cdid_ia_get_elt(cdid_ia_get_elt(panel_s, "observations"), "calendar_index");
          id_s_st = cdid_ia_get_elt(st_s, "id");
          if (TYPEOF(id_type_s) != STRSXP || XLENGTH(id_type_s) != 1 ||
              strcmp(CHAR(STRING_ELT(id_type_s, 0)), "numeric") != 0 ||
              TYPEOF(ukey_s) != REALSXP || TYPEOF(units_s) != VECSXP || TYPEOF(cal_s) != INTSXP ||
              (TYPEOF(id_s_st) != REALSXP && TYPEOF(id_s_st) != INTSXP) || Rf_isObject(id_s_st)) {
            vmaxset(vmax);
            return R_NilValue;
          }
          SEXP nobs_s = cdid_ia_get_elt(units_s, "n_obs");
          obs_start_s_u = cdid_ia_get_elt(units_s, "obs_start");
          if (TYPEOF(nobs_s) != INTSXP ||
              (TYPEOF(obs_start_s_u) != INTSXP && TYPEOF(obs_start_s_u) != REALSXP)) {
            vmaxset(vmax);
            return R_NilValue;
          }
          ukey_p = REAL(ukey_s);
          n_panel_units = (int) XLENGTH(ukey_s);
          u_nobs_p = INTEGER(nobs_s);
          o_cal_p = INTEGER(cal_s);
        }
      }
      double *up_buf = (double *) R_alloc((size_t) B * (size_t) n_uids, sizeof(double));
      memset(up_buf, 0, (size_t) B * (size_t) n_uids * sizeof(double));
      int n_a_up = Rf_asInteger(cdid_ia_get_elt(up_map_s, "n_a"));
      const int *cohorts_a = INTEGER(cdid_ia_get_elt(up_map_s, "cohorts_a"));
      for (int j = 0; j < n_uids; ++j) {
        int u = uids[j];
        if (u <= n_a_up) {
          int k = cohorts_a[u - 1];
          double mk = mass_U[k - 1];
          double ak = a_vec[k - 1];
          double *dst = up_buf + (size_t) j * (size_t) B;
          for (int r = 0; r < B; ++r) {
            if (ci[r] == k && U_rows[r]) {
              double A_r = (ai[r] > 0) ? 1.0 : 0.0;
              dst[r] = w[r] * (A_r - ak) / mk;
            }
          }
        }
      }

      SEXP theta_coords_s = cdid_ia_get_elt(up_map_s, "theta_coords");
      int M_sel = (int) XLENGTH(sel_blocks_s);
      for (int m = 0; m < M_sel; ++m) {
        SEXP tc_s = VECTOR_ELT(theta_coords_s, m);
        int n_tc = (int) XLENGTH(tc_s);
        if (n_tc == 0) continue;
        const int *tc = INTEGER(tc_s);
        int any_active = 0;
        for (int k = 0; k < n_tc; ++k) {
          if (tc[k] >= 1 && tc[k] <= max_uid && uid_pos[tc[k]] > 0) { any_active = 1; break; }
        }
        if (!any_active) continue;

        SEXP bl_s = VECTOR_ELT(sel_blocks_s, m);
        const char *branch = CHAR(STRING_ELT(cdid_ia_get_elt(bl_s, "branch"), 0));
        if (strcmp(branch, "zero_col") == 0) continue;
        int B_m = Rf_asInteger(cdid_ia_get_elt(bl_s, "B_m"));
        int ds_bl = Rf_asInteger(cdid_ia_get_elt(bl_s, "d_s"));
        int has_int = (Rf_asLogical(cdid_ia_get_elt(bl_s, "has_intercept")) == 1);
        const int *bl_cids = INTEGER(cdid_ia_get_elt(bl_s, "cell_ids"));

        double *c_scores = (double *) R_alloc((size_t) B * (size_t) B_m, sizeof(double));
        memset(c_scores, 0, (size_t) B * (size_t) B_m * sizeof(double));
        for (int b = 0; b < B_m; ++b) {
          int scid = bl_cids[b];
          SEXP scell = VECTOR_ELT(sel_cells_s, scid - 1);
          int scoh = Rf_asInteger(cdid_ia_get_elt(scell, "cohort"));
          double *cs_b = c_scores + (size_t) b * (size_t) B;

          if (is_indep_period) {
            SEXP dspec_c = cdid_ia_get_elt(scell, "date_spec");
            int dt = Rf_asInteger(cdid_ia_get_elt(dspec_c != R_NilValue ? dspec_c : scell, "calendar_index"));
            if (dt == NA_INTEGER || dt < 1 || dt > n_dates) {
              vmaxset(vmax);
              return R_NilValue;
            }
            int has_dt = (date_count[dt - 1] > 0);
            int t_lo = date_start[dt - 1] - 1, t_hi = date_end[dt - 1] - 1;
            int ht_hint = t_lo, hu_hint = 0;
            for (int r = 0; r < B; ++r) {
              if (ci[r] != scoh || !U_rows[r]) continue;
              int row0 = rows[r] - 1;
              double pr = 0.0;
              if (!cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N, d_full, slope_cols, scid, row0, &pr)) {
                vmaxset(vmax);
                return R_NilValue;
              }
              double Y = 0.0;
              if (ai[r] > 0) {
                if (has_dt && cdid_ia_find_posting(post_a, t_lo, t_hi, ai[r], &ht_hint) >= 0) Y = 1.0;
              } else if (n_sing > 0) {
                double id_d = (TYPEOF(id_s_st) == INTSXP) ? (double) INTEGER(id_s_st)[row0] : REAL(id_s_st)[row0];
                if (id_d == 0.0) id_d = 0.0;
                int u = -1;
                if (hu_hint >= 0 && hu_hint + 1 < n_panel_units && ukey_p[hu_hint + 1] == id_d) {
                  u = ++hu_hint;
                } else {
                  int lo = 0, hi = n_panel_units - 1;
                  while (lo <= hi) {
                    int mid = lo + ((hi - lo) >> 1);
                    double v = ukey_p[mid];
                    if (v == id_d) { u = hu_hint = mid; break; }
                    if (v < id_d) lo = mid + 1;
                    else hi = mid - 1;
                  }
                }
                if (u >= 0 && u_nobs_p[u] == 1) {
                  int opos = (TYPEOF(obs_start_s_u) == INTSXP)
                    ? (INTEGER(obs_start_s_u)[u] - 1)
                    : ((int) REAL(obs_start_s_u)[u] - 1);
                  if (opos >= 0 && o_cal_p[opos] == dt) Y = 1.0;
                }
              }
              cs_b[r] = w[r] * (Y - pr);
            }
          } else {
            int ds = 0, dt = 0;
            if (!cdid_ia_resolve_pair_dates(cdid_ia_get_elt(scell, "pair"), all_times, n_dates, &ds, &dt)) {
              vmaxset(vmax);
              return R_NilValue;
            }
            int has_pair = (date_count[ds - 1] > 0 && date_count[dt - 1] > 0);
            int s_lo = date_start[ds - 1] - 1, s_hi = date_end[ds - 1] - 1;
            int t_lo = date_start[dt - 1] - 1, t_hi = date_end[dt - 1] - 1;
            int hs_hint = s_lo, ht_hint = t_lo;

            for (int r = 0; r < B; ++r) {
              if (ci[r] != scoh || !U_rows[r]) continue;
              int row0 = rows[r] - 1;
              double pr = 0.0;
              if (!cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N, d_full, slope_cols, scid, row0, &pr)) {
                vmaxset(vmax);
                return R_NilValue;
              }
              double Y = 0.0;
              if (has_pair && ai[r] > 0) {
                int hs = cdid_ia_find_posting(post_a, s_lo, s_hi, ai[r], &hs_hint);
                if (hs >= 0) {
                  int ht = cdid_ia_find_posting(post_a, t_lo, t_hi, ai[r], &ht_hint);
                  if (ht >= 0) Y = 1.0;
                }
              }
              cs_b[r] = w[r] * (Y - pr);
            }
          }
        }

        if (strcmp(branch, "cp04_separate_pair") == 0) {
          const int *coord = INTEGER(cdid_ia_get_elt(bl_s, "coordinate_order"));
          double *score_mat = (double *) R_alloc((size_t) B * (size_t) d_full, sizeof(double));
          double *xi_raw    = (double *) R_alloc((size_t) B * (size_t) d_full, sizeof(double));
          memset(score_mat, 0, (size_t) B * (size_t) d_full * sizeof(double));
          for (int r = 0; r < B; ++r) {
            if (!U_rows[r]) continue;
            double s1 = c_scores[r];
            if (s1 == 0.0) continue;
            int row0 = rows[r] - 1;
            for (int k = 0; k < d_full; ++k) {
              double can_k = (has_int && k == 0) ? 1.0 : Z_sel[row0 + (size_t) (slope_cols[k - (has_int ? 1 : 0)] - 1) * (size_t) N];
              int dst_col = coord[k] - 1;
              score_mat[r + (size_t) dst_col * (size_t) B] = s1 * can_k;
            }
          }
          const double *H_chol = REAL(cdid_ia_get_elt(bl_s, "H_chol"));
          if (!cdid_ia_chol_solve_rows(H_chol, d_full, score_mat, B, xi_raw)) {
            vmaxset(vmax);
            return R_NilValue;
          }
          for (int k = 0; k < d_full && k < n_tc; ++k) {
            int pos1 = (tc[k] >= 1 && tc[k] <= max_uid) ? uid_pos[tc[k]] : 0;
            if (pos1 > 0) {
              int src_col = coord[k] - 1;
              const double *src = xi_raw + (size_t) src_col * (size_t) B;
              double *dst = up_buf + (size_t) (pos1 - 1) * (size_t) B;
              for (int r = 0; r < B; ++r) dst[r] = U_rows[r] ? src[r] : 0.0;
            }
          }
        } else {
          double *xi_beta_bl = (ds_bl > 0) ? (double *) R_alloc((size_t) B * (size_t) ds_bl, sizeof(double)) : NULL;
          if (ds_bl > 0) {
            double *s_tilde = (double *) R_alloc((size_t) B * (size_t) ds_bl, sizeof(double));
            memset(s_tilde, 0, (size_t) B * (size_t) ds_bl * sizeof(double));
            const double *z_bar = has_int ? REAL(cdid_ia_get_elt(bl_s, "z_bar")) : NULL;
            for (int r = 0; r < B; ++r) {
              if (!U_rows[r]) continue;
              double sum_s = 0.0;
              for (int b = 0; b < B_m; ++b) sum_s += c_scores[r + (size_t) b * (size_t) B];
              int row0 = rows[r] - 1;
              for (int col = 0; col < ds_bl; ++col) {
                double val = sum_s * Z_sel[row0 + (size_t) (slope_cols[col] - 1) * (size_t) N];
                if (has_int) {
                  for (int b = 0; b < B_m; ++b) {
                    val -= c_scores[r + (size_t) b * (size_t) B] * z_bar[b + (size_t) col * (size_t) B_m];
                  }
                }
                s_tilde[r + (size_t) col * (size_t) B] = val;
              }
            }
            const double *R_sl = NULL;
            const int *piv_sl = NULL;
            if (!cdid_ia_get_wls_fit(cdid_ia_get_elt(bl_s, "slope_fit"), &R_sl, &piv_sl, ds_bl) ||
                !cdid_ia_gram_solve_rows(R_sl, piv_sl, ds_bl, s_tilde, B, (double) N, xi_beta_bl)) {
              vmaxset(vmax);
              return R_NilValue;
            }
          }
          int n_al = has_int ? B_m : 0;
          if (has_int) {
            const double *m_b_vec = REAL(cdid_ia_get_elt(bl_s, "m_b"));
            const double *z_bar = (ds_bl > 0) ? REAL(cdid_ia_get_elt(bl_s, "z_bar")) : NULL;
            for (int b = 0; b < B_m; ++b) {
              int pos1 = (tc[b] >= 1 && tc[b] <= max_uid) ? uid_pos[tc[b]] : 0;
              if (pos1 <= 0) continue;
              double *dst = up_buf + (size_t) (pos1 - 1) * (size_t) B;
              double mb = m_b_vec[b];
              const double *cs_b = c_scores + (size_t) b * (size_t) B;
              for (int r = 0; r < B; ++r) {
                if (!U_rows[r]) { dst[r] = 0.0; continue; }
                double v = cs_b[r] / mb;
                for (int col = 0; col < ds_bl; ++col) {
                  v -= xi_beta_bl[r + (size_t) col * (size_t) B] * z_bar[b + (size_t) col * (size_t) B_m];
                }
                dst[r] = v;
              }
            }
          }
          for (int col = 0; col < ds_bl; ++col) {
            int tc_u = tc[n_al + col];
            int pos1 = (tc_u >= 1 && tc_u <= max_uid) ? uid_pos[tc_u] : 0;
            if (pos1 > 0) {
              double *dst = up_buf + (size_t) (pos1 - 1) * (size_t) B;
              const double *src = xi_beta_bl + (size_t) col * (size_t) B;
              for (int r = 0; r < B; ++r) dst[r] = U_rows[r] ? src[r] : 0.0;
            }
          }
        }
      }
      up = up_buf;
    }
  }

  /* 3 & 4. Reconstructed cell responses Z_c (for zids) and imputation forward/adjoint IFs (for cids/gids) */
  SEXP s1f_s = cdid_ia_get_elt(sets_s, "stage1_forward");
  int use_imp_forward = (s1f_s != R_NilValue && Rf_asLogical(cdid_ia_get_elt(s1f_s, "enabled")) == TRUE);

  int *need_c = (int *) R_alloc((size_t) (n_cells + 1), sizeof(int));
  int *need_z = (int *) R_alloc((size_t) (n_cells + 1), sizeof(int));
  int *cell_gid = (int *) R_alloc((size_t) (n_cells + 1), sizeof(int));
  memset(need_c, 0, (size_t) (n_cells + 1) * sizeof(int));
  memset(need_z, 0, (size_t) (n_cells + 1) * sizeof(int));
  memset(cell_gid, 0, (size_t) (n_cells + 1) * sizeof(int));
  for (int c = 1; c <= n_cells; ++c) {
    cell_gid[c] = Rf_asInteger(cdid_ia_get_elt(VECTOR_ELT(cells_s, c - 1), "group_id"));
  }
  for (int i = 0; i < n_cids; ++i) {
    int c = cids[i];
    if (c >= 1 && c <= n_cells) need_c[c] = 1;
  }
  for (int i = 0; i < n_zids; ++i) {
    int c = zids[i];
    if (c >= 1 && c <= n_cells) need_z[c] = 1;
  }

  double **Z_cell   = (double **) R_alloc((size_t) (n_cells + 1), sizeof(double *));
  double **alpha_if = (double **) R_alloc((size_t) (n_cells + 1), sizeof(double *));
  for (int c = 0; c <= n_cells; ++c) { Z_cell[c] = NULL; alpha_if[c] = NULL; }

  int *cell_n_obs      = (!use_imp_forward) ? (int *) R_alloc((size_t) (n_cells + 1), sizeof(int)) : NULL;
  int **cell_obs_rr    = (!use_imp_forward) ? (int **) R_alloc((size_t) (n_cells + 1), sizeof(int *)) : NULL;
  double **cell_s_obs  = (!use_imp_forward) ? (double **) R_alloc((size_t) (n_cells + 1), sizeof(double *)) : NULL;
  double **cell_dt_obs = (!use_imp_forward) ? (double **) R_alloc((size_t) (n_cells + 1), sizeof(double *)) : NULL;
  int *tmp_rr_buf      = (!use_imp_forward) ? (int *) R_alloc((size_t) B, sizeof(int)) : NULL;
  double *tmp_s_buf    = (!use_imp_forward) ? (double *) R_alloc((size_t) B, sizeof(double)) : NULL;
  double *tmp_dt_buf   = (!use_imp_forward) ? (double *) R_alloc((size_t) B * 256, sizeof(double)) : NULL;
  double *t_lam_tmp    = (!use_imp_forward) ? (double *) R_alloc((size_t) (n_cells + 1), sizeof(double)) : NULL;
  double **t_nu_tmp    = (!use_imp_forward) ? (double **) R_alloc((size_t) (n_groups + 1), sizeof(double *)) : NULL;
  double **t_nu_alloc  = (!use_imp_forward) ? (double **) R_alloc((size_t) (n_groups + 1), sizeof(double *)) : NULL;
  double *zero_pg_buf  = (!use_imp_forward) ? (double *) R_alloc(256, sizeof(double)) : NULL;
  if (!use_imp_forward) {
    memset(cell_n_obs, 0, (size_t) (n_cells + 1) * sizeof(int));
    memset(zero_pg_buf, 0, 256 * sizeof(double));
    for (int c = 0; c <= n_cells; ++c) {
      cell_obs_rr[c] = NULL;
      cell_s_obs[c] = NULL;
      cell_dt_obs[c] = NULL;
    }
    for (int g = 0; g <= n_groups; ++g) {
      t_nu_tmp[g] = NULL;
      t_nu_alloc[g] = NULL;
    }
  }

  double **h_beta  = (double **) R_alloc((size_t) (n_groups + 1), sizeof(double *));
  double **xi_beta = (double **) R_alloc((size_t) (n_groups + 1), sizeof(double *));
  int *p_group     = (int *) R_alloc((size_t) (n_groups + 1), sizeof(int));
  for (int g = 0; g <= n_groups; ++g) { h_beta[g] = NULL; xi_beta[g] = NULL; p_group[g] = 0; }
  for (int gid = 1; gid <= n_groups; ++gid) {
    SEXP gf = VECTOR_ELT(groups_s, gid - 1);
    int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
    if (pg > 256) { vmaxset(vmax); return R_NilValue; }
    p_group[gid] = pg;
    if (!use_imp_forward && pg > 0) {
      t_nu_alloc[gid] = (double *) R_alloc((size_t) pg, sizeof(double));
    }
  }

  if (use_imp_forward) {
    for (int ig = 0; ig < n_gids; ++ig) {
      int gid = gids[ig];
      if (gid < 1 || gid > n_groups) { vmaxset(vmax); return R_NilValue; }
      int pg = p_group[gid];
      if (pg > 0) {
        h_beta[gid]  = (double *) R_alloc((size_t) B * (size_t) pg, sizeof(double));
        xi_beta[gid] = (double *) R_alloc((size_t) B * (size_t) pg, sizeof(double));
        memset(h_beta[gid],  0, (size_t) B * (size_t) pg * sizeof(double));
        memset(xi_beta[gid], 0, (size_t) B * (size_t) pg * sizeof(double));
      }
    }
  }

  double d_buf[256];
  for (int cid = 1; cid <= n_cells; ++cid) {
    int do_c = need_c[cid];
    int do_z = need_z[cid];
    if (!do_c && !do_z) continue;

    SEXP cf = VECTOR_ELT(cells_s, cid - 1);
    int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
    int c_coh = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
    int q_id = Rf_asInteger(cdid_ia_get_elt(cf, "q_id"));
    double L_val = Rf_asReal(cdid_ia_get_elt(cf, "L"));
    double alpha_c = Rf_asReal(cdid_ia_get_elt(cf, "alpha"));
    double mass_c = Rf_asReal(cdid_ia_get_elt(cf, "mass"));
    SEXP gf = VECTOR_ELT(groups_s, gid - 1);
    int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
    if (pg > 256) { vmaxset(vmax); return R_NilValue; }
    const int *kept_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf, "kept_cols_F")) : NULL;
    const double *beta_g = (pg > 0) ? REAL(cdid_ia_get_elt(gf, "beta")) : NULL;
    const double *mean_d = (pg > 0 && intercept) ? REAL(cdid_ia_get_elt(cf, "mean_d")) : NULL;

    int ds = 0, dt = 0;
    if (!cdid_ia_resolve_pair_dates(cdid_ia_get_elt(cf, "pair"), all_times, n_dates, &ds, &dt)) {
      vmaxset(vmax);
      return R_NilValue;
    }
    int has_pair = (date_count[ds - 1] > 0 && date_count[dt - 1] > 0);
    int s_lo = date_start[ds - 1] - 1, s_hi = date_end[ds - 1] - 1;
    int t_lo = date_start[dt - 1] - 1, t_hi = date_end[dt - 1] - 1;

    double *z_c = do_z ? (double *) R_alloc((size_t) B, sizeof(double)) : NULL;
    double *a_if_c = (do_c && use_imp_forward) ? (double *) R_alloc((size_t) B, sizeof(double)) : NULL;
    if (do_z) { memset(z_c, 0, (size_t) B * sizeof(double)); Z_cell[cid] = z_c; }
    if (do_c && use_imp_forward) { memset(a_if_c, 0, (size_t) B * sizeof(double)); alpha_if[cid] = a_if_c; }
    double *hb = (do_c && use_imp_forward && gid >= 1 && gid <= n_groups) ? h_beta[gid] : NULL;
    int n_ob_c = 0;

    int hs_hint = s_lo, ht_hint = t_lo;
    for (int r = 0; r < B; ++r) {
      if (ci[r] != c_coh) continue;
      int row0 = rows[r] - 1;
      double jhat_sl = 0.0;
      for (int m = 0; m < pg; ++m) {
        double d_m = L_val * X_F[row0 + (size_t) (kept_F[m] - 1) * (size_t) N];
        d_buf[m] = d_m;
        jhat_sl += d_m * beta_g[m];
      }
      double jhat = alpha_c + jhat_sl;
      if (do_z) z_c[r] = jhat;

      if (has_pair && ai[r] > 0) {
        int hs = cdid_ia_find_posting(post_a, s_lo, s_hi, ai[r], &hs_hint);
        if (hs >= 0) {
          int ht = cdid_ia_find_posting(post_a, t_lo, t_hi, ai[r], &ht_hint);
          if (ht >= 0) {
            int os = post_obs[hs], ot = post_obs[ht];
            if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y) {
              vmaxset(vmax);
              return R_NilValue;
            }
            double dy = y_obs[ot - 1] - y_obs[os - 1];
            if (!R_FINITE(dy)) { vmaxset(vmax); return R_NilValue; }
            double e_r = dy - jhat;
            if (do_z) {
              double q = 0.0;
              if (!cdid_ia_eval_sel_pair_prob(sel_pcs_s, sel_cells_s, sel_blocks_s, is_indep_period,
                                              Z_sel, N, d_full, slope_cols, q_id, row0, &q)) {
                vmaxset(vmax);
                return R_NilValue;
              }
              double q_fl = (q < 1e-6) ? 1e-6 : q;
              double pi_r = 1.0 / q_fl;
              z_c[r] = jhat + pi_r * e_r;
            }
            if (do_c) {
              double s_rr = w[r] * e_r;
              if (use_imp_forward) {
                if (intercept) a_if_c[r] = s_rr / mass_c;
                if (pg > 0 && hb != NULL) {
                  for (int m = 0; m < pg; ++m) {
                    double dt_m = intercept ? (d_buf[m] - mean_d[m]) : d_buf[m];
                    hb[r + (size_t) m * (size_t) B] += s_rr * dt_m;
                  }
                }
              } else {
                tmp_rr_buf[n_ob_c] = r;
                tmp_s_buf[n_ob_c] = s_rr;
                for (int m = 0; m < pg; ++m) {
                  tmp_dt_buf[n_ob_c + (size_t) m * (size_t) B] = intercept ? (d_buf[m] - mean_d[m]) : d_buf[m];
                }
                n_ob_c++;
              }
            }
          }
        }
      }
    }
    if (do_c && !use_imp_forward && n_ob_c > 0) {
      cell_n_obs[cid] = n_ob_c;
      cell_obs_rr[cid] = (int *) R_alloc((size_t) n_ob_c, sizeof(int));
      cell_s_obs[cid]  = (double *) R_alloc((size_t) n_ob_c, sizeof(double));
      memcpy(cell_obs_rr[cid], tmp_rr_buf, (size_t) n_ob_c * sizeof(int));
      memcpy(cell_s_obs[cid],  tmp_s_buf,  (size_t) n_ob_c * sizeof(double));
      if (pg > 0) {
        double *dt_c = (double *) R_alloc((size_t) n_ob_c * (size_t) pg, sizeof(double));
        for (int m = 0; m < pg; ++m) {
          memcpy(dt_c + (size_t) m * (size_t) n_ob_c, tmp_dt_buf + (size_t) m * (size_t) B, (size_t) n_ob_c * sizeof(double));
        }
        cell_dt_obs[cid] = dt_c;
      }
    }
  }

  if (use_imp_forward) {
    for (int ig = 0; ig < n_gids; ++ig) {
      int gid = gids[ig];
      int pg = p_group[gid];
      if (pg > 0) {
        SEXP gf = VECTOR_ELT(groups_s, gid - 1);
        const double *R_g = NULL;
        const int *piv_g = NULL;
        if (!cdid_ia_get_wls_fit(cdid_ia_get_elt(gf, "fit"), &R_g, &piv_g, pg) ||
            !cdid_ia_gram_solve_rows(R_g, piv_g, pg, h_beta[gid], B, (double) N, xi_beta[gid])) {
          vmaxset(vmax);
          return R_NilValue;
        }
        if (intercept) {
          SEXP g_cids_s = cdid_ia_get_elt(gf, "cell_ids");
          int n_gc = (int) XLENGTH(g_cids_s);
          const int *g_cids = INTEGER(g_cids_s);
          const double *xb = xi_beta[gid];
          for (int k = 0; k < n_gc; ++k) {
            int cid = g_cids[k];
            if (cid >= 1 && cid <= n_cells && alpha_if[cid] != NULL) {
              double *a_if_c = alpha_if[cid];
              const double *mean_d = REAL(cdid_ia_get_elt(VECTOR_ELT(cells_s, cid - 1), "mean_d"));
              for (int r = 0; r < B; ++r) {
                double dot = 0.0;
                for (int m = 0; m < pg; ++m) {
                  dot = fma(xb[r + (size_t) m * (size_t) B], mean_d[m], dot);
                }
                a_if_c[r] -= dot;
              }
            }
          }
        }
      }
    }
  }

  /* 5. Shared / Separate Counterfactual OR IFs (direct, a, j, q channels) */
  int has_or = (shared_or_s != R_NilValue && n_bids > 0);
  int separate_or = 0;
  int n_or_blocks = 0;
  SEXP or_blocks_s = R_NilValue;
  if (has_or) {
    SEXP op_s = cdid_ia_get_elt(shared_or_s, "outcome_pooling");
    if (TYPEOF(op_s) == STRSXP && XLENGTH(op_s) == 1 &&
        strcmp(CHAR(STRING_ELT(op_s, 0)), "separate") == 0) {
      separate_or = 1;
    }
    or_blocks_s = cdid_ia_get_elt(shared_or_s, "blocks");
    n_or_blocks = (int) XLENGTH(or_blocks_s);
  }

  double **or_alpha_dir = NULL, **or_alpha_a = NULL, **or_alpha_j = NULL, **or_alpha_q = NULL;
  double **or_beta_dir  = NULL, **or_beta_a  = NULL, **or_beta_j  = NULL, **or_beta_q  = NULL;
  int *p_or_by_bid = NULL;
  int **kept_or_by_bid = NULL;
  const double **beta_or_by_bid = NULL;

  if (has_or) {
    or_alpha_dir = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_alpha_a   = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_alpha_j   = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_alpha_q   = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_beta_dir  = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_beta_a    = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_beta_j    = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    or_beta_q    = (double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(double *));
    p_or_by_bid  = (int *) R_alloc((size_t) (n_or_blocks + 1), sizeof(int));
    kept_or_by_bid = (int **) R_alloc((size_t) (n_or_blocks + 1), sizeof(int *));
    beta_or_by_bid = (const double **) R_alloc((size_t) (n_or_blocks + 1), sizeof(const double *));
    for (int b = 0; b <= n_or_blocks; ++b) {
      or_alpha_dir[b] = or_alpha_a[b] = or_alpha_j[b] = or_alpha_q[b] = NULL;
      or_beta_dir[b]  = or_beta_a[b]  = or_beta_j[b]  = or_beta_q[b]  = NULL;
      p_or_by_bid[b] = 0;
      kept_or_by_bid[b] = NULL;
      beta_or_by_bid[b] = NULL;
    }

    const int *slope_cols_F_all = separate_or ? INTEGER(cdid_ia_get_elt(shared_or_s, "slope_cols_F")) : NULL;
    SEXP og_groups_s = R_NilValue;
    const int *bm_gid = NULL;
    if (separate_or) {
      SEXP og_s = cdid_ia_get_elt(shared_or_s, "outcome_groups");
      og_groups_s = cdid_ia_get_elt(og_s, "groups");
      bm_gid = INTEGER(cdid_ia_get_elt(cdid_ia_get_elt(og_s, "block_map"), "group_id"));
    }

    int p_or_pool = separate_or ? 0 : Rf_asInteger(cdid_ia_get_elt(shared_or_s, "p"));
    if (p_or_pool > 256) { vmaxset(vmax); return R_NilValue; }
    int *kept_pool = (!separate_or && p_or_pool > 0) ? INTEGER(cdid_ia_get_elt(shared_or_s, "kept_cols_F")) : NULL;
    const double *beta_pool = (!separate_or && p_or_pool > 0) ? REAL(cdid_ia_get_elt(shared_or_s, "beta")) : NULL;

    double *pool_rhs_b_dir = (!separate_or && p_or_pool > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double)) : NULL;
    double *pool_rhs_b_a   = (!separate_or && p_or_pool > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double)) : NULL;
    double *pool_rhs_b_j   = (!separate_or && p_or_pool > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double)) : NULL;
    double *pool_rhs_b_q   = (!separate_or && p_or_pool > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double)) : NULL;
    if (!separate_or && p_or_pool > 0) {
      memset(pool_rhs_b_dir, 0, (size_t) B * (size_t) p_or_pool * sizeof(double));
      memset(pool_rhs_b_a,   0, (size_t) B * (size_t) p_or_pool * sizeof(double));
      memset(pool_rhs_b_j,   0, (size_t) B * (size_t) p_or_pool * sizeof(double));
      memset(pool_rhs_b_q,   0, (size_t) B * (size_t) p_or_pool * sizeof(double));
      SEXP q_b_ids_s = cdid_ia_get_elt(shared_or_s, "q_beta_ids");
      int n_qb = (int) XLENGTH(q_b_ids_s);
      if (n_qb > 0) {
        const int *q_b_ids = INTEGER(q_b_ids_s);
        const double *q_b = REAL(cdid_ia_get_elt(shared_or_s, "q_beta"));
        for (int k = 0; k < n_qb; ++k) {
          int u = q_b_ids[k];
          int pos1 = (u >= 1 && u <= max_uid) ? uid_pos[u] : 0;
          if (pos1 <= 0 || up == NULL) { vmaxset(vmax); return R_NilValue; }
          const double *up_col = up + (size_t) (pos1 - 1) * (size_t) B;
          for (int m = 0; m < p_or_pool; ++m) {
            double v = q_b[m + (size_t) k * (size_t) p_or_pool];
            double *dst = pool_rhs_b_q + (size_t) m * (size_t) B;
            for (int r = 0; r < B; ++r) dst[r] = fma(up_col[r], v, dst[r]);
          }
        }
      }
    }

    int *active_bid_mask = (int *) R_alloc((size_t) (n_or_blocks + 1), sizeof(int));
    memset(active_bid_mask, 0, (size_t) (n_or_blocks + 1) * sizeof(int));
    for (int ib = 0; ib < n_bids; ++ib) {
      int bid = bids[ib];
      if (bid >= 1 && bid <= n_or_blocks) active_bid_mask[bid] = 1;
    }
    for (int j = 0; j < nl; ++j) {
      SEXP lf = VECTOR_ELT(links_all_s, links[j] - 1);
      int bid = Rf_asInteger(cdid_ia_get_elt(lf, "or_block_id"));
      if (bid >= 1 && bid <= n_or_blocks) active_bid_mask[bid] = 1;
    }

    for (int bid = 1; bid <= n_or_blocks; ++bid) {
      if (!active_bid_mask[bid]) continue;
      SEXP bf = VECTOR_ELT(or_blocks_s, bid - 1);
      int p_or = p_or_pool;
      int *kept_or_F = kept_pool;
      const double *beta_or = beta_pool;
      SEXP fit_or_s = R_NilValue;
      if (separate_or) {
        SEXP og_g = VECTOR_ELT(og_groups_s, bm_gid[bid - 1] - 1);
        SEXP keep_s = cdid_ia_get_elt(og_g, "keep");
        p_or = (int) XLENGTH(keep_s);
        if (p_or > 256) { vmaxset(vmax); return R_NilValue; }
        const int *keep_idx = INTEGER(keep_s);
        kept_or_F = (p_or > 0) ? (int *) R_alloc((size_t) p_or, sizeof(int)) : NULL;
        for (int m = 0; m < p_or; ++m) kept_or_F[m] = slope_cols_F_all[keep_idx[m] - 1];
        beta_or = (p_or > 0) ? REAL(cdid_ia_get_elt(og_g, "beta")) : NULL;
        fit_or_s = cdid_ia_get_elt(og_g, "fit");
      }
      p_or_by_bid[bid] = p_or;
      kept_or_by_bid[bid] = kept_or_F;
      beta_or_by_bid[bid] = beta_or;

      double *ra_dir = (double *) R_alloc((size_t) B, sizeof(double));
      double *ra_a   = (double *) R_alloc((size_t) B, sizeof(double));
      double *ra_j   = (double *) R_alloc((size_t) B, sizeof(double));
      double *ra_q   = (double *) R_alloc((size_t) B, sizeof(double));
      memset(ra_dir, 0, (size_t) B * sizeof(double));
      memset(ra_a,   0, (size_t) B * sizeof(double));
      memset(ra_j,   0, (size_t) B * sizeof(double));
      memset(ra_q,   0, (size_t) B * sizeof(double));
      or_alpha_dir[bid] = ra_dir;
      or_alpha_a[bid]   = ra_a;
      or_alpha_j[bid]   = ra_j;
      or_alpha_q[bid]   = ra_q;

      double *rb_dir = separate_or ? ((p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL) : pool_rhs_b_dir;
      double *rb_a   = separate_or ? ((p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL) : pool_rhs_b_a;
      double *rb_j   = separate_or ? ((p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL) : pool_rhs_b_j;
      double *rb_q   = separate_or ? ((p_or > 0) ? (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double)) : NULL) : pool_rhs_b_q;
      if (separate_or && p_or > 0) {
        memset(rb_dir, 0, (size_t) B * (size_t) p_or * sizeof(double));
        memset(rb_a,   0, (size_t) B * (size_t) p_or * sizeof(double));
        memset(rb_j,   0, (size_t) B * (size_t) p_or * sizeof(double));
        memset(rb_q,   0, (size_t) B * (size_t) p_or * sizeof(double));
      }

      int in_wanted = 0;
      for (int ib = 0; ib < n_bids; ++ib) {
        if (bids[ib] == bid) { in_wanted = 1; break; }
      }
      if (!in_wanted) continue;

      double L_b = Rf_asReal(cdid_ia_get_elt(bf, "L"));
      double alpha_b = intercept ? Rf_asReal(cdid_ia_get_elt(bf, "alpha")) : 0.0;
      const double *mean_d_b = (p_or > 0 && intercept) ? REAL(cdid_ia_get_elt(bf, "mean_d")) : NULL;
      SEXP ctrl_cids_s = cdid_ia_get_elt(bf, "control_cell_ids");
      int n_kb = (int) XLENGTH(ctrl_cids_s);
      const int *ctrl_cids = INTEGER(ctrl_cids_s);
      const double *a_force_al = REAL(cdid_ia_get_elt(bf, "a_force_alpha"));
      const double *a_force_be = (p_or > 0) ? REAL(cdid_ia_get_elt(bf, "a_force_beta")) : NULL;
      const double *j_force_aa = REAL(cdid_ia_get_elt(bf, "j_force_alpha_alpha"));
      SEXP j_force_ab_s = cdid_ia_get_elt(bf, "j_force_alpha_beta");
      const double *j_force_ba = (p_or > 0) ? REAL(cdid_ia_get_elt(bf, "j_force_beta_alpha")) : NULL;
      SEXP j_force_bb_s = cdid_ia_get_elt(bf, "j_force_beta_beta");

      for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
        int cid = ctrl_cids[idx_k];
        SEXP cf = VECTOR_ELT(cells_s, cid - 1);
        int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
        int pg = p_group[gid];
        int ci_k = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
        double a_k = a_vec[ci_k - 1];
        const double *z_c = Z_cell[cid];
        if (z_c == NULL) { vmaxset(vmax); return R_NilValue; }

        for (int r = 0; r < B; ++r) {
          if (ci[r] != ci_k) continue;
          int row0 = rows[r] - 1;
          double w_syn = w[r] * a_k;
          double mhat_sl = 0.0;
          for (int m = 0; m < p_or; ++m) {
            double d_or_m = L_b * X_F[row0 + (size_t) (kept_or_F[m] - 1) * (size_t) N];
            d_buf[m] = d_or_m;
            mhat_sl += d_or_m * beta_or[m];
          }
          double mhat = alpha_b + mhat_sl;
          double wr = w_syn * (z_c[r] - mhat);
          if (intercept) ra_dir[r] += wr;
          for (int m = 0; m < p_or; ++m) {
            double dt_m = intercept ? (d_buf[m] - mean_d_b[m]) : d_buf[m];
            rb_dir[r + (size_t) m * (size_t) B] += wr * dt_m;
          }
        }

        const double *xi_ak = xi_a + (size_t) (ci_k - 1) * (size_t) B;
        if (intercept) {
          double afa = a_force_al[idx_k];
          if (afa != 0.0) {
            for (int r = 0; r < B; ++r) ra_a[r] += afa * xi_ak[r];
          }
        }
        if (p_or > 0) {
          const double *afb = a_force_be + (size_t) idx_k * (size_t) p_or;
          for (int m = 0; m < p_or; ++m) {
            double afb_m = afb[m];
            if (afb_m != 0.0) {
              double *dst = rb_a + (size_t) m * (size_t) B;
              for (int r = 0; r < B; ++r) dst[r] += xi_ak[r] * afb_m;
            }
          }
        }

        if (use_imp_forward) {
          const double *xi_ac = alpha_if[cid];
          const double *xi_br = xi_beta[gid];
          if (intercept) {
            double jfaa = j_force_aa[idx_k];
            if (jfaa != 0.0 && xi_ac != NULL) {
              for (int r = 0; r < B; ++r) ra_j[r] += jfaa * xi_ac[r];
            }
            if (pg > 0 && xi_br != NULL) {
              const double *jfab = REAL(VECTOR_ELT(j_force_ab_s, idx_k));
              for (int r = 0; r < B; ++r) {
                double dot = 0.0;
                for (int m = 0; m < pg; ++m) {
                  dot = fma(xi_br[r + (size_t) m * (size_t) B], jfab[m], dot);
                }
                ra_j[r] += dot;
              }
            }
          }
          if (p_or > 0) {
            if (intercept && xi_ac != NULL) {
              const double *jfba = j_force_ba + (size_t) idx_k * (size_t) p_or;
              for (int q = 0; q < p_or; ++q) {
                double v = jfba[q];
                if (v != 0.0) {
                  double *dst = rb_j + (size_t) q * (size_t) B;
                  for (int r = 0; r < B; ++r) dst[r] += xi_ac[r] * v;
                }
              }
            }
            if (pg > 0 && xi_br != NULL) {
              const double *jfbb = REAL(VECTOR_ELT(j_force_bb_s, idx_k));
              for (int q = 0; q < p_or; ++q) {
                double *dst = rb_j + (size_t) q * (size_t) B;
                for (int r = 0; r < B; ++r) {
                  double dot = 0.0;
                  for (int m = 0; m < pg; ++m) {
                    dot = fma(xi_br[r + (size_t) m * (size_t) B], jfbb[q + (size_t) m * (size_t) p_or], dot);
                  }
                  dst[r] += dot;
                }
              }
            }
          }
        }
      }

      if (!use_imp_forward && separate_or) {
        double *lam_lin = (double *) R_alloc((size_t) (n_cells + 1), sizeof(double));
        double **nu_lin = (double **) R_alloc((size_t) (n_cells + 1), sizeof(double *));
        for (int c = 1; c <= n_cells; ++c) {
          int pg_c = p_group[cell_gid[c]];
          nu_lin[c] = (pg_c > 0) ? (double *) R_alloc((size_t) pg_c, sizeof(double)) : NULL;
        }
        for (int q = 0; q < p_or; ++q) {
          memset(lam_lin, 0, (size_t) (n_cells + 1) * sizeof(double));
          for (int c = 1; c <= n_cells; ++c) {
            int pg_c = p_group[cell_gid[c]];
            if (pg_c > 0) memset(nu_lin[c], 0, (size_t) pg_c * sizeof(double));
          }
          for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
            int cid = ctrl_cids[idx_k];
            int pg_c = p_group[cell_gid[cid]];
            if (intercept) lam_lin[cid] += j_force_ba[q + (size_t) idx_k * (size_t) p_or];
            if (pg_c > 0) {
              const double *jfbb = REAL(VECTOR_ELT(j_force_bb_s, idx_k));
              for (int m = 0; m < pg_c; ++m) nu_lin[cid][m] += jfbb[q + (size_t) m * (size_t) p_or];
            }
          }
          if (!cdid_ia_eval_imp_linear_chunk(
                lam_lin, nu_lin, cells_s, groups_s, n_cells, n_groups, N,
                cids, n_cids, cell_gid, p_group, intercept,
                cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs, B,
                t_lam_tmp, t_nu_alloc, rb_j + (size_t) q * (size_t) B)) {
            vmaxset(vmax);
            return R_NilValue;
          }
        }
        if (intercept) {
          memset(lam_lin, 0, (size_t) (n_cells + 1) * sizeof(double));
          for (int c = 1; c <= n_cells; ++c) {
            int pg_c = p_group[cell_gid[c]];
            if (pg_c > 0) memset(nu_lin[c], 0, (size_t) pg_c * sizeof(double));
          }
          for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
            int cid = ctrl_cids[idx_k];
            int pg_c = p_group[cell_gid[cid]];
            lam_lin[cid] += j_force_aa[idx_k];
            if (pg_c > 0) {
              const double *jfab = REAL(VECTOR_ELT(j_force_ab_s, idx_k));
              for (int m = 0; m < pg_c; ++m) nu_lin[cid][m] += jfab[m];
            }
          }
          if (!cdid_ia_eval_imp_linear_chunk(
                lam_lin, nu_lin, cells_s, groups_s, n_cells, n_groups, N,
                cids, n_cids, cell_gid, p_group, intercept,
                cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs, B,
                t_lam_tmp, t_nu_alloc, ra_j)) {
            vmaxset(vmax);
            return R_NilValue;
          }
        }
      }

      if (intercept) {
        SEXP q_a_ids_s = cdid_ia_get_elt(bf, "q_alpha_ids");
        int n_qa = (int) XLENGTH(q_a_ids_s);
        if (n_qa > 0) {
          const int *q_a_ids = INTEGER(q_a_ids_s);
          const double *q_a = REAL(cdid_ia_get_elt(bf, "q_alpha"));
          for (int k = 0; k < n_qa; ++k) {
            int u = q_a_ids[k];
            int pos1 = (u >= 1 && u <= max_uid) ? uid_pos[u] : 0;
            if (pos1 <= 0 || up == NULL) { vmaxset(vmax); return R_NilValue; }
            double v = q_a[k];
            const double *up_col = up + (size_t) (pos1 - 1) * (size_t) B;
            for (int r = 0; r < B; ++r) ra_q[r] = fma(up_col[r], v, ra_q[r]);
          }
        }
      }

      if (separate_or) {
        if (p_or > 0) {
          SEXP q_b_ids_s = cdid_ia_get_elt(bf, "q_beta_ids");
          int n_qb = (int) XLENGTH(q_b_ids_s);
          if (n_qb > 0) {
            const int *q_b_ids = INTEGER(q_b_ids_s);
            const double *q_b = REAL(cdid_ia_get_elt(bf, "q_beta"));
            for (int k = 0; k < n_qb; ++k) {
              int u = q_b_ids[k];
              int pos1 = (u >= 1 && u <= max_uid) ? uid_pos[u] : 0;
              if (pos1 <= 0 || up == NULL) { vmaxset(vmax); return R_NilValue; }
              const double *up_col = up + (size_t) (pos1 - 1) * (size_t) B;
              for (int m = 0; m < p_or; ++m) {
                double v = q_b[m + (size_t) k * (size_t) p_or];
                double *dst = rb_q + (size_t) m * (size_t) B;
                for (int r = 0; r < B; ++r) dst[r] = fma(up_col[r], v, dst[r]);
              }
            }
          }
        }
        if (intercept) {
          double sma = Rf_asReal(cdid_ia_get_elt(bf, "score_mean_alpha"));
          for (int r = 0; r < B; ++r) ra_dir[r] -= sma;
        }
        if (p_or > 0) {
          const double *smb = REAL(cdid_ia_get_elt(bf, "score_mean_beta"));
          for (int m = 0; m < p_or; ++m) {
            double smbm = smb[m];
            double *dst = rb_dir + (size_t) m * (size_t) B;
            for (int r = 0; r < B; ++r) dst[r] -= smbm;
          }
          double *xb_dir = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
          double *xb_a   = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
          double *xb_j   = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
          double *xb_q   = (double *) R_alloc((size_t) B * (size_t) p_or, sizeof(double));
          const double *R_or = NULL;
          const int *piv_or = NULL;
          if (!cdid_ia_get_wls_fit(fit_or_s, &R_or, &piv_or, p_or) ||
              !cdid_ia_gram_solve_rows(R_or, piv_or, p_or, rb_dir, B, (double) N, xb_dir) ||
              !cdid_ia_gram_solve_rows(R_or, piv_or, p_or, rb_a,   B, (double) N, xb_a)   ||
              !cdid_ia_gram_solve_rows(R_or, piv_or, p_or, rb_j,   B, (double) N, xb_j)   ||
              !cdid_ia_gram_solve_rows(R_or, piv_or, p_or, rb_q,   B, (double) N, xb_q)) {
            vmaxset(vmax);
            return R_NilValue;
          }
          or_beta_dir[bid] = xb_dir;
          or_beta_a[bid]   = xb_a;
          or_beta_j[bid]   = xb_j;
          or_beta_q[bid]   = xb_q;
        }
        if (intercept) {
          double mass_b = Rf_asReal(cdid_ia_get_elt(bf, "mass"));
          for (int r = 0; r < B; ++r) {
            double dot_dir = 0.0, dot_a = 0.0, dot_j = 0.0, dot_q = 0.0;
            for (int m = 0; m < p_or; ++m) {
              double md = mean_d_b[m];
              dot_dir = fma(or_beta_dir[bid][r + (size_t) m * (size_t) B], md, dot_dir);
              dot_a   = fma(or_beta_a[bid][r   + (size_t) m * (size_t) B], md, dot_a);
              dot_j   = fma(or_beta_j[bid][r   + (size_t) m * (size_t) B], md, dot_j);
              dot_q   = fma(or_beta_q[bid][r   + (size_t) m * (size_t) B], md, dot_q);
            }
            ra_dir[r] = (ra_dir[r] / mass_b) - dot_dir;
            ra_a[r]   = (ra_a[r]   / mass_b) - dot_a;
            ra_j[r]   = (ra_j[r]   / mass_b) - dot_j;
            ra_q[r]   = (ra_q[r]   / mass_b) - dot_q;
          }
        }
      }
    }

    if (!separate_or) {
      if (!use_imp_forward) {
        double *lam_lin = (double *) R_alloc((size_t) (n_cells + 1), sizeof(double));
        double **nu_lin = (double **) R_alloc((size_t) (n_cells + 1), sizeof(double *));
        for (int c = 1; c <= n_cells; ++c) {
          int pg_c = p_group[cell_gid[c]];
          nu_lin[c] = (pg_c > 0) ? (double *) R_alloc((size_t) pg_c, sizeof(double)) : NULL;
        }
        for (int q = 0; q < p_or_pool; ++q) {
          memset(lam_lin, 0, (size_t) (n_cells + 1) * sizeof(double));
          for (int c = 1; c <= n_cells; ++c) {
            int pg_c = p_group[cell_gid[c]];
            if (pg_c > 0) memset(nu_lin[c], 0, (size_t) pg_c * sizeof(double));
          }
          for (int ib = 0; ib < n_bids; ++ib) {
            int bid = bids[ib];
            SEXP bf = VECTOR_ELT(or_blocks_s, bid - 1);
            SEXP ctrl_cids_s = cdid_ia_get_elt(bf, "control_cell_ids");
            int n_kb = (int) XLENGTH(ctrl_cids_s);
            const int *ctrl_cids = INTEGER(ctrl_cids_s);
            const double *j_force_ba = REAL(cdid_ia_get_elt(bf, "j_force_beta_alpha"));
            SEXP j_force_bb_s = cdid_ia_get_elt(bf, "j_force_beta_beta");
            for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
              int cid = ctrl_cids[idx_k];
              int pg_c = p_group[cell_gid[cid]];
              if (intercept) lam_lin[cid] += j_force_ba[q + (size_t) idx_k * (size_t) p_or_pool];
              if (pg_c > 0) {
                const double *jfbb = REAL(VECTOR_ELT(j_force_bb_s, idx_k));
                for (int m = 0; m < pg_c; ++m) nu_lin[cid][m] += jfbb[q + (size_t) m * (size_t) p_or_pool];
              }
            }
          }
          if (!cdid_ia_eval_imp_linear_chunk(
                lam_lin, nu_lin, cells_s, groups_s, n_cells, n_groups, N,
                cids, n_cids, cell_gid, p_group, intercept,
                cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs, B,
                t_lam_tmp, t_nu_alloc, pool_rhs_b_j + (size_t) q * (size_t) B)) {
            vmaxset(vmax);
            return R_NilValue;
          }
        }
        if (intercept) {
          for (int bid = 1; bid <= n_or_blocks; ++bid) {
            if (!active_bid_mask[bid]) continue;
            int in_wanted = 0;
            for (int ib = 0; ib < n_bids; ++ib) {
              if (bids[ib] == bid) { in_wanted = 1; break; }
            }
            if (!in_wanted) continue;
            SEXP bf = VECTOR_ELT(or_blocks_s, bid - 1);
            SEXP ctrl_cids_s = cdid_ia_get_elt(bf, "control_cell_ids");
            int n_kb = (int) XLENGTH(ctrl_cids_s);
            const int *ctrl_cids = INTEGER(ctrl_cids_s);
            const double *j_force_aa = REAL(cdid_ia_get_elt(bf, "j_force_alpha_alpha"));
            SEXP j_force_ab_s = cdid_ia_get_elt(bf, "j_force_alpha_beta");
            memset(lam_lin, 0, (size_t) (n_cells + 1) * sizeof(double));
            for (int c = 1; c <= n_cells; ++c) {
              int pg_c = p_group[cell_gid[c]];
              if (pg_c > 0) memset(nu_lin[c], 0, (size_t) pg_c * sizeof(double));
            }
            for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
              int cid = ctrl_cids[idx_k];
              int pg_c = p_group[cell_gid[cid]];
              lam_lin[cid] += j_force_aa[idx_k];
              if (pg_c > 0) {
                const double *jfab = REAL(VECTOR_ELT(j_force_ab_s, idx_k));
                for (int m = 0; m < pg_c; ++m) nu_lin[cid][m] += jfab[m];
              }
            }
            if (!cdid_ia_eval_imp_linear_chunk(
                  lam_lin, nu_lin, cells_s, groups_s, n_cells, n_groups, N,
                  cids, n_cids, cell_gid, p_group, intercept,
                  cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs, B,
                  t_lam_tmp, t_nu_alloc, or_alpha_j[bid])) {
              vmaxset(vmax);
              return R_NilValue;
            }
          }
        }
      }
      double *xb_dir = NULL, *xb_a = NULL, *xb_j = NULL, *xb_q = NULL;
      if (p_or_pool > 0) {
        xb_dir = (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double));
        xb_a   = (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double));
        xb_j   = (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double));
        xb_q   = (double *) R_alloc((size_t) B * (size_t) p_or_pool, sizeof(double));
        const double *R_or = NULL;
        const int *piv_or = NULL;
        if (!cdid_ia_get_wls_fit(cdid_ia_get_elt(shared_or_s, "fit"), &R_or, &piv_or, p_or_pool) ||
            !cdid_ia_gram_solve_rows(R_or, piv_or, p_or_pool, pool_rhs_b_dir, B, (double) N, xb_dir) ||
            !cdid_ia_gram_solve_rows(R_or, piv_or, p_or_pool, pool_rhs_b_a,   B, (double) N, xb_a)   ||
            !cdid_ia_gram_solve_rows(R_or, piv_or, p_or_pool, pool_rhs_b_j,   B, (double) N, xb_j)   ||
            !cdid_ia_gram_solve_rows(R_or, piv_or, p_or_pool, pool_rhs_b_q,   B, (double) N, xb_q)) {
          vmaxset(vmax);
          return R_NilValue;
        }
      }
      for (int bid = 1; bid <= n_or_blocks; ++bid) {
        if (!active_bid_mask[bid]) continue;
        or_beta_dir[bid] = xb_dir;
        or_beta_a[bid]   = xb_a;
        or_beta_j[bid]   = xb_j;
        or_beta_q[bid]   = xb_q;
        if (intercept) {
          SEXP bf = VECTOR_ELT(or_blocks_s, bid - 1);
          double mass_b = Rf_asReal(cdid_ia_get_elt(bf, "mass"));
          double *ra_dir = or_alpha_dir[bid];
          double *ra_a   = or_alpha_a[bid];
          double *ra_j   = or_alpha_j[bid];
          double *ra_q   = or_alpha_q[bid];
          const double *mean_d_b = (p_or_pool > 0) ? REAL(cdid_ia_get_elt(bf, "mean_d")) : NULL;
          for (int r = 0; r < B; ++r) {
            double dot_dir = 0.0, dot_a = 0.0, dot_j = 0.0, dot_q = 0.0;
            for (int m = 0; m < p_or_pool; ++m) {
              double md = mean_d_b[m];
              dot_dir = fma(xb_dir[r + (size_t) m * (size_t) B], md, dot_dir);
              dot_a   = fma(xb_a[r   + (size_t) m * (size_t) B], md, dot_a);
              dot_j   = fma(xb_j[r   + (size_t) m * (size_t) B], md, dot_j);
              dot_q   = fma(xb_q[r   + (size_t) m * (size_t) B], md, dot_q);
            }
            ra_dir[r] = (ra_dir[r] / mass_b) - dot_dir;
            ra_a[r]   = (ra_a[r]   / mass_b) - dot_a;
            ra_j[r]   = (ra_j[r]   / mass_b) - dot_j;
            ra_q[r]   = (ra_q[r]   / mass_b) - dot_q;
          }
        }
      }
    }
  }

  /* 6. Roster treatment logit IFs (direct and a channels) + cached odds */
  int has_trt = (treatments_s != R_NilValue && n_tids > 0);
  int n_trt_tot = (treatments_s != R_NilValue) ? (int) XLENGTH(treatments_s) : 0;
  int p_trt = (trt_recipe_s != R_NilValue) ? Rf_asInteger(cdid_ia_get_elt(trt_recipe_s, "p")) : 0;
  if (p_trt > 256) { vmaxset(vmax); return R_NilValue; }
  const int *trt_cols_F = (p_trt > 0) ? INTEGER(cdid_ia_get_elt(trt_recipe_s, "index")) : NULL;

  double **trt_dir = NULL, **trt_a = NULL, **trt_odds = NULL;
  if (has_trt) {
    trt_dir  = (double **) R_alloc((size_t) (n_trt_tot + 1), sizeof(double *));
    trt_a    = (double **) R_alloc((size_t) (n_trt_tot + 1), sizeof(double *));
    trt_odds = (double **) R_alloc((size_t) (n_trt_tot + 1), sizeof(double *));
    for (int t = 0; t <= n_trt_tot; ++t) {
      trt_dir[t] = trt_a[t] = trt_odds[t] = NULL;
    }
    double *score_dir = (p_trt > 0) ? (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double)) : NULL;
    double *score_a   = (p_trt > 0) ? (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double)) : NULL;
    int *in_all_coh   = (int *) R_alloc((size_t) (G + 1), sizeof(int));

    for (int it = 0; it < n_tids; ++it) {
      int tid = tids[it];
      if (tid < 1 || tid > n_trt_tot) { vmaxset(vmax); return R_NilValue; }
      SEXP tf = VECTOR_ELT(treatments_s, tid - 1);
      int target_coh = Rf_asInteger(cdid_ia_get_elt(tf, "target_cohort_index"));
      SEXP all_coh_s = cdid_ia_get_elt(tf, "all_cohort_indices");
      int n_ac = (int) XLENGTH(all_coh_s);
      const int *all_coh = INTEGER(all_coh_s);
      const double *coef_trt = (p_trt > 0) ? REAL(cdid_ia_get_elt(tf, "coefficients")) : NULL;

      int rep_tid = 0;
      for (int it2 = 0; it2 < it; ++it2) {
        int tid2 = tids[it2];
        SEXP tf2 = VECTOR_ELT(treatments_s, tid2 - 1);
        if (Rf_asInteger(cdid_ia_get_elt(tf2, "target_cohort_index")) != target_coh) continue;
        SEXP ac2_s = cdid_ia_get_elt(tf2, "all_cohort_indices");
        if ((int) XLENGTH(ac2_s) != n_ac) continue;
        if (memcmp(INTEGER(ac2_s), all_coh, (size_t) n_ac * sizeof(int)) != 0) continue;
        if (p_trt > 0 && memcmp(REAL(cdid_ia_get_elt(tf2, "coefficients")), coef_trt, (size_t) p_trt * sizeof(double)) != 0) continue;
        rep_tid = tid2;
        break;
      }
      if (rep_tid > 0) {
        trt_dir[tid]  = trt_dir[rep_tid];
        trt_a[tid]    = trt_a[rep_tid];
        trt_odds[tid] = trt_odds[rep_tid];
        continue;
      }

      memset(in_all_coh, 0, (size_t) (G + 1) * sizeof(int));
      for (int k = 0; k < n_ac; ++k) {
        int c = all_coh[k];
        if (c >= 1 && c <= G) in_all_coh[c] = 1;
      }

      double *odds_vec = (double *) R_alloc((size_t) B, sizeof(double));
      memset(odds_vec, 0, (size_t) B * sizeof(double));
      trt_odds[tid] = odds_vec;

      if (p_trt > 0) {
        memset(score_dir, 0, (size_t) B * (size_t) p_trt * sizeof(double));
        memset(score_a,   0, (size_t) B * (size_t) p_trt * sizeof(double));
        for (int r = 0; r < B; ++r) {
          int c = ci[r];
          if (c < 1 || c > G || !in_all_coh[c]) continue;
          int row0 = rows[r] - 1;
          double eta = 0.0;
          for (int m = 0; m < p_trt; ++m) {
            double x_rm = X_F[row0 + (size_t) (trt_cols_F[m] - 1) * (size_t) N];
            d_buf[m] = x_rm;
            eta = fma(x_rm, coef_trt[m], eta);
          }
          if (!R_FINITE(eta)) { vmaxset(vmax); return R_NilValue; }
          double pr = 0.0, od = 0.0;
          cdid_ia_clip_logit(eta, &pr, &od);
          odds_vec[r] = od;
          double y_r = (c == target_coh) ? 1.0 : 0.0;
          double sc = (w[r] * a_vec[c - 1]) * (y_r - pr);
          for (int m = 0; m < p_trt; ++m) {
            score_dir[r + (size_t) m * (size_t) B] = sc * d_buf[m];
          }
        }
        const double *a_forc = REAL(cdid_ia_get_elt(tf, "a_forcing"));
        for (int k = 0; k < n_ac; ++k) {
          int c = all_coh[k];
          const double *xi_ac = xi_a + (size_t) (c - 1) * (size_t) B;
          for (int m = 0; m < p_trt; ++m) {
            double af = a_forc[m + (size_t) k * (size_t) p_trt];
            double *dst = score_a + (size_t) m * (size_t) B;
            for (int r = 0; r < B; ++r) dst[r] = fma(xi_ac[r], af, dst[r]);
          }
        }
        double *xg_dir = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
        double *xg_a   = (double *) R_alloc((size_t) B * (size_t) p_trt, sizeof(double));
        const double *H_chol = REAL(cdid_ia_get_elt(tf, "H_chol"));
        if (!cdid_ia_chol_solve_rows(H_chol, p_trt, score_dir, B, xg_dir) ||
            !cdid_ia_chol_solve_rows(H_chol, p_trt, score_a,   B, xg_a)) {
          vmaxset(vmax);
          return R_NilValue;
        }
        trt_dir[tid] = xg_dir;
        trt_a[tid]   = xg_a;
      } else {
        for (int r = 0; r < B; ++r) odds_vec[r] = 1.0;
      }
    }
  }

  /* 7. Link IFs assembly */
  SEXP phi_s = PROTECT(Rf_allocMatrix(REALSXP, B, nl));
  SEXP link_names_s = PROTECT(Rf_allocVector(STRSXP, nl));
  double *phi = REAL(phi_s);

  double *w_plus    = (double *) R_alloc((size_t) B, sizeof(double));
  double *Z_g_vec   = (double *) R_alloc((size_t) B, sizeof(double));
  double *m_g_vec   = (double *) R_alloc((size_t) B, sizeof(double));
  double *j_T_vec   = (double *) R_alloc((size_t) B, sizeof(double));
  double *sel_T_vec = (double *) R_alloc((size_t) B, sizeof(double));
  double *w_minus   = (double *) R_alloc((size_t) B, sizeof(double));
  double *Z_k_vec   = (double *) R_alloc((size_t) B, sizeof(double));
  double *m_k_vec   = (double *) R_alloc((size_t) B, sizeof(double));
  double *a_C_vec   = (double *) R_alloc((size_t) B, sizeof(double));
  double *j_C_vec   = (double *) R_alloc((size_t) B, sizeof(double));
  double *sel_C_vec = (double *) R_alloc((size_t) B, sizeof(double));
  double *a_Bo_vec  = (double *) R_alloc((size_t) B, sizeof(double));
  int *ctrl_coh_map = (int *) R_alloc((size_t) (G + 1), sizeof(int));

  for (int j = 0; j < nl; ++j) {
    int ell = links[j];
    SEXP lf = VECTOR_ELT(links_all_s, ell - 1);
    SET_STRING_ELT(link_names_s, j, STRING_ELT(cdid_ia_get_elt(lf, "link_id"), 0));

    int bid = Rf_asInteger(cdid_ia_get_elt(lf, "or_block_id"));
    int t_cid = Rf_asInteger(cdid_ia_get_elt(lf, "target_cell_id"));
    int ci_g = Rf_asInteger(cdid_ia_get_elt(lf, "target_cohort_index"));
    double a_g = a_vec[ci_g - 1];
    double L_lk = Rf_asReal(cdid_ia_get_elt(lf, "L"));
    double d_plus = Rf_asReal(cdid_ia_get_elt(lf, "d_plus"));
    double T_gl = Rf_asReal(cdid_ia_get_elt(lf, "T_gl"));
    double a_coeff_T_g = Rf_asReal(cdid_ia_get_elt(lf, "a_coeff_T_g"));

    SEXP cf_g = VECTOR_ELT(cells_s, t_cid - 1);
    int gid_g = Rf_asInteger(cdid_ia_get_elt(cf_g, "group_id"));
    int pg_g = p_group[gid_g];

    int p_or = (has_or && bid >= 1 && bid <= n_or_blocks) ? p_or_by_bid[bid] : 0;
    const int *kept_or_F = (has_or && bid >= 1 && bid <= n_or_blocks) ? kept_or_by_bid[bid] : NULL;
    const double *beta_or = (has_or && bid >= 1 && bid <= n_or_blocks) ? beta_or_by_bid[bid] : NULL;
    double alpha_bf = (has_or && intercept) ? Rf_asReal(cdid_ia_get_elt(VECTOR_ELT(or_blocks_s, bid - 1), "alpha")) : 0.0;

    const double *zg_src = Z_cell[t_cid];
    if (zg_src == NULL) { UNPROTECT(2); vmaxset(vmax); return R_NilValue; }

    memset(w_plus,  0, (size_t) B * sizeof(double));
    memset(Z_g_vec, 0, (size_t) B * sizeof(double));
    memset(m_g_vec, 0, (size_t) B * sizeof(double));
    for (int r = 0; r < B; ++r) {
      if (ci[r] != ci_g) continue;
      int row0 = rows[r] - 1;
      w_plus[r] = w[r] * a_g;
      Z_g_vec[r] = zg_src[r];
      if (has_or) {
        double mv_sl = 0.0;
        for (int m = 0; m < p_or; ++m) {
          mv_sl += (L_lk * X_F[row0 + (size_t) (kept_or_F[m] - 1) * (size_t) N]) * beta_or[m];
        }
        m_g_vec[r] = alpha_bf + mv_sl;
      }
    }

    if (use_imp_forward) {
      memset(j_T_vec, 0, (size_t) B * sizeof(double));
      double lam_T = Rf_asReal(cdid_ia_get_elt(lf, "lam_T"));
      if (intercept && lam_T != 0.0 && alpha_if[t_cid] != NULL) {
        const double *ait = alpha_if[t_cid];
        for (int r = 0; r < B; ++r) j_T_vec[r] = lam_T * ait[r];
      }
      SEXP nu_T_s = cdid_ia_get_elt(lf, "nu_T");
      if (pg_g > 0 && XLENGTH(nu_T_s) == (R_xlen_t) pg_g && xi_beta[gid_g] != NULL) {
        const double *nu_T = REAL(nu_T_s);
        const double *xb = xi_beta[gid_g];
        for (int r = 0; r < B; ++r) {
          double dot = 0.0;
          for (int m = 0; m < pg_g; ++m) {
            dot = fma(xb[r + (size_t) m * (size_t) B], nu_T[m], dot);
          }
          j_T_vec[r] += dot;
        }
      }
    } else {
      if (!cdid_ia_eval_sexp_adjoint_chunk(
            cdid_ia_get_elt(lf, "adjoint_T"), n_cells, n_groups, cids, n_cids,
            cell_gid, p_group, intercept, cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs,
            B, t_lam_tmp, t_nu_tmp, zero_pg_buf, j_T_vec)) {
        UNPROTECT(2); vmaxset(vmax); return R_NilValue;
      }
    }

    memset(sel_T_vec, 0, (size_t) B * sizeof(double));
    SEXP sel_T_ids_s = cdid_ia_get_elt(lf, "sel_T_ids");
    int n_st = (int) XLENGTH(sel_T_ids_s);
    if (n_st > 0) {
      const int *st_ids = INTEGER(sel_T_ids_s);
      const double *st_coef = REAL(cdid_ia_get_elt(lf, "sel_T"));
      for (int k = 0; k < n_st; ++k) {
        int u = st_ids[k];
        int pos1 = (u >= 1 && u <= max_uid) ? uid_pos[u] : 0;
        if (pos1 <= 0 || up == NULL) { UNPROTECT(2); vmaxset(vmax); return R_NilValue; }
        double v = st_coef[k];
        const double *up_col = up + (size_t) (pos1 - 1) * (size_t) B;
        for (int r = 0; r < B; ++r) sel_T_vec[r] = fma(up_col[r], v, sel_T_vec[r]);
      }
    }

    int tid = has_trt ? Rf_asInteger(cdid_ia_get_elt(lf, "trt_id")) : 0;
    double d_minus = has_trt ? Rf_asReal(cdid_ia_get_elt(lf, "d_minus")) : 1.0;
    double C_ol = has_trt ? Rf_asReal(cdid_ia_get_elt(lf, "C_ol")) : 0.0;
    double B_ol = (has_trt && has_or) ? Rf_asReal(cdid_ia_get_elt(lf, "B_ol")) : 0.0;
    double B_gl = has_or ? Rf_asReal(cdid_ia_get_elt(lf, "B_gl")) : 0.0;
    double a_coeff_Bg_g = has_or ? Rf_asReal(cdid_ia_get_elt(lf, "a_coeff_Bg_g")) : 0.0;
    const double *d_bar_plus_or  = (has_or && p_or > 0) ? REAL(cdid_ia_get_elt(lf, "d_bar_plus_or")) : NULL;
    const double *d_bar_minus_or = (has_trt && has_or && p_or > 0) ? REAL(cdid_ia_get_elt(lf, "d_bar_minus_or")) : NULL;

    if (has_trt) {
      memset(w_minus,   0, (size_t) B * sizeof(double));
      memset(Z_k_vec,   0, (size_t) B * sizeof(double));
      memset(m_k_vec,   0, (size_t) B * sizeof(double));
      memset(a_C_vec,   0, (size_t) B * sizeof(double));
      memset(j_C_vec,   0, (size_t) B * sizeof(double));
      memset(sel_C_vec, 0, (size_t) B * sizeof(double));
      memset(a_Bo_vec,  0, (size_t) B * sizeof(double));

      SEXP ctrl_cohs_s = cdid_ia_get_elt(lf, "control_cohort_indices");
      SEXP ctrl_cids_s = cdid_ia_get_elt(lf, "control_cell_ids");
      int n_ctrl = (int) XLENGTH(ctrl_cohs_s);
      const int *ctrl_cohs = INTEGER(ctrl_cohs_s);
      const int *ctrl_cids = INTEGER(ctrl_cids_s);
      const double *a_coeff_C_k = REAL(cdid_ia_get_elt(lf, "a_coeff_C_k"));
      const double *lam_C_k = REAL(cdid_ia_get_elt(lf, "lam_C_k"));
      SEXP nu_C_k_s = cdid_ia_get_elt(lf, "nu_C_k");
      const double *a_coeff_Bo_k = has_or ? REAL(cdid_ia_get_elt(lf, "a_coeff_Bo_k")) : NULL;
      const double *odds_t = trt_odds[tid];

      memset(ctrl_coh_map, 0, (size_t) (G + 1) * sizeof(int));
      for (int idx_k = 0; idx_k < n_ctrl; ++idx_k) {
        int ci_k = ctrl_cohs[idx_k];
        if (ci_k >= 1 && ci_k <= G) ctrl_coh_map[ci_k] = idx_k + 1;
      }

      for (int r = 0; r < B; ++r) {
        int c = ci[r];
        int idx1 = (c >= 1 && c <= G) ? ctrl_coh_map[c] : 0;
        if (idx1 <= 0) continue;
        int idx_k = idx1 - 1;
        int cid_k = ctrl_cids[idx_k];
        const double *zk_src = Z_cell[cid_k];
        if (zk_src == NULL) { UNPROTECT(2); vmaxset(vmax); return R_NilValue; }
        int row0 = rows[r] - 1;
        w_minus[r] = (w[r] * a_vec[c - 1]) * odds_t[r];
        Z_k_vec[r] = zk_src[r];
        if (has_or) {
          double mv_sl = 0.0;
          for (int m = 0; m < p_or; ++m) {
            mv_sl += (L_lk * X_F[row0 + (size_t) (kept_or_F[m] - 1) * (size_t) N]) * beta_or[m];
          }
          m_k_vec[r] = alpha_bf + mv_sl;
        }
      }

      for (int idx_k = 0; idx_k < n_ctrl; ++idx_k) {
        int ci_k = ctrl_cohs[idx_k];
        int cid_k = ctrl_cids[idx_k];
        SEXP cf_k = VECTOR_ELT(cells_s, cid_k - 1);
        int gid_k = Rf_asInteger(cdid_ia_get_elt(cf_k, "group_id"));
        int pg_k = p_group[gid_k];
        const double *xi_ak = xi_a + (size_t) (ci_k - 1) * (size_t) B;
        double ack = a_coeff_C_k[idx_k];
        if (ack != 0.0) {
          for (int r = 0; r < B; ++r) a_C_vec[r] += ack * xi_ak[r];
        }
        if (has_or) {
          double abok = a_coeff_Bo_k[idx_k];
          if (abok != 0.0) {
            for (int r = 0; r < B; ++r) a_Bo_vec[r] += abok * xi_ak[r];
          }
        }
        if (use_imp_forward) {
          if (intercept && alpha_if[cid_k] != NULL) {
            double lck = lam_C_k[idx_k];
            const double *aic = alpha_if[cid_k];
            if (lck != 0.0) {
              for (int r = 0; r < B; ++r) j_C_vec[r] += lck * aic[r];
            }
          }
          if (pg_k > 0 && xi_beta[gid_k] != NULL) {
            const double *nuck = REAL(VECTOR_ELT(nu_C_k_s, idx_k));
            const double *xb = xi_beta[gid_k];
            for (int r = 0; r < B; ++r) {
              double dot = 0.0;
              for (int m = 0; m < pg_k; ++m) {
                dot = fma(xb[r + (size_t) m * (size_t) B], nuck[m], dot);
              }
              j_C_vec[r] += dot;
            }
          }
        }
      }
      if (!use_imp_forward) {
        if (!cdid_ia_eval_sexp_adjoint_chunk(
              cdid_ia_get_elt(lf, "adjoint_C"), n_cells, n_groups, cids, n_cids,
              cell_gid, p_group, intercept, cell_n_obs, cell_obs_rr, cell_s_obs, cell_dt_obs,
              B, t_lam_tmp, t_nu_tmp, zero_pg_buf, j_C_vec)) {
          UNPROTECT(2); vmaxset(vmax); return R_NilValue;
        }
      }

      SEXP sel_C_ids_s = cdid_ia_get_elt(lf, "sel_C_ids");
      int n_sc = (int) XLENGTH(sel_C_ids_s);
      if (n_sc > 0) {
        const int *sc_ids = INTEGER(sel_C_ids_s);
        const double *sc_coef = REAL(cdid_ia_get_elt(lf, "sel_C"));
        for (int k = 0; k < n_sc; ++k) {
          int u = sc_ids[k];
          int pos1 = (u >= 1 && u <= max_uid) ? uid_pos[u] : 0;
          if (pos1 <= 0 || up == NULL) { UNPROTECT(2); vmaxset(vmax); return R_NilValue; }
          double v = sc_coef[k];
          const double *up_col = up + (size_t) (pos1 - 1) * (size_t) B;
          for (int r = 0; r < B; ++r) sel_C_vec[r] = fma(up_col[r], v, sel_C_vec[r]);
        }
      }
    }

    const double *xi_ag = xi_a + (size_t) (ci_g - 1) * (size_t) B;
    const double *g_C_gam  = (has_trt && p_trt > 0) ? REAL(cdid_ia_get_elt(lf, "g_C_gamma")) : NULL;
    const double *g_Bo_gam = (has_trt && has_or && p_trt > 0) ? REAL(cdid_ia_get_elt(lf, "g_Bo_gamma")) : NULL;
    const double *t_dir = (has_trt && p_trt > 0) ? trt_dir[tid] : NULL;
    const double *t_a   = (has_trt && p_trt > 0) ? trt_a[tid] : NULL;
    const double *oa_dir = (has_or && intercept) ? or_alpha_dir[bid] : NULL;
    const double *oa_a   = (has_or && intercept) ? or_alpha_a[bid] : NULL;
    const double *oa_j   = (has_or && intercept) ? or_alpha_j[bid] : NULL;
    const double *oa_q   = (has_or && intercept) ? or_alpha_q[bid] : NULL;
    const double *ob_dir = (has_or && p_or > 0) ? or_beta_dir[bid] : NULL;
    const double *ob_a   = (has_or && p_or > 0) ? or_beta_a[bid] : NULL;
    const double *ob_j   = (has_or && p_or > 0) ? or_beta_j[bid] : NULL;
    const double *ob_q   = (has_or && p_or > 0) ? or_beta_q[bid] : NULL;

    double *phi_col = phi + (size_t) j * (size_t) B;
    for (int r = 0; r < B; ++r) {
      double direct_T = w_plus[r] * (Z_g_vec[r] - T_gl) / d_plus;
      double a_T = a_coeff_T_g * xi_ag[r];
      double j_T = j_T_vec[r];
      double sel_T = sel_T_vec[r];

      double direct_Bg = 0.0, a_Bg = 0.0, or_Bg_dir = 0.0, or_Bg_a = 0.0, or_Bg_j = 0.0, or_Bg_q = 0.0;
      if (has_or) {
        direct_Bg = w_plus[r] * (m_g_vec[r] - B_gl) / d_plus;
        a_Bg = a_coeff_Bg_g * xi_ag[r];
        for (int m = 0; m < p_or; ++m) {
          double dbp = d_bar_plus_or[m];
          or_Bg_dir = fma(ob_dir[r + (size_t) m * (size_t) B], dbp, or_Bg_dir);
          or_Bg_a   = fma(ob_a[r   + (size_t) m * (size_t) B], dbp, or_Bg_a);
          or_Bg_j   = fma(ob_j[r   + (size_t) m * (size_t) B], dbp, or_Bg_j);
          or_Bg_q   = fma(ob_q[r   + (size_t) m * (size_t) B], dbp, or_Bg_q);
        }
        if (intercept) {
          or_Bg_dir = oa_dir[r] + or_Bg_dir;
          or_Bg_a   = oa_a[r]   + or_Bg_a;
          or_Bg_j   = oa_j[r]   + or_Bg_j;
          or_Bg_q   = oa_q[r]   + or_Bg_q;
        }
      }

      double direct_C = 0.0, a_C = 0.0, j_C = 0.0, sel_C = 0.0, trt_C_dir = 0.0, trt_C_a = 0.0;
      double direct_Bo = 0.0, a_Bo = 0.0, trt_Bo_dir = 0.0, trt_Bo_a = 0.0;
      double or_Bo_dir = 0.0, or_Bo_a = 0.0, or_Bo_j = 0.0, or_Bo_q = 0.0;
      if (has_trt) {
        direct_C = w_minus[r] * (Z_k_vec[r] - C_ol) / d_minus;
        a_C = a_C_vec[r];
        j_C = j_C_vec[r];
        sel_C = sel_C_vec[r];
        for (int m = 0; m < p_trt; ++m) {
          double gcg = g_C_gam[m];
          trt_C_dir = fma(t_dir[r + (size_t) m * (size_t) B], gcg, trt_C_dir);
          trt_C_a   = fma(t_a[r   + (size_t) m * (size_t) B], gcg, trt_C_a);
        }
        if (has_or) {
          direct_Bo = w_minus[r] * (m_k_vec[r] - B_ol) / d_minus;
          a_Bo = a_Bo_vec[r];
          for (int m = 0; m < p_trt; ++m) {
            double gbog = g_Bo_gam[m];
            trt_Bo_dir = fma(t_dir[r + (size_t) m * (size_t) B], gbog, trt_Bo_dir);
            trt_Bo_a   = fma(t_a[r   + (size_t) m * (size_t) B], gbog, trt_Bo_a);
          }
          for (int m = 0; m < p_or; ++m) {
            double dbm = d_bar_minus_or[m];
            or_Bo_dir = fma(ob_dir[r + (size_t) m * (size_t) B], dbm, or_Bo_dir);
            or_Bo_a   = fma(ob_a[r   + (size_t) m * (size_t) B], dbm, or_Bo_a);
            or_Bo_j   = fma(ob_j[r   + (size_t) m * (size_t) B], dbm, or_Bo_j);
            or_Bo_q   = fma(ob_q[r   + (size_t) m * (size_t) B], dbm, or_Bo_q);
          }
          if (intercept) {
            or_Bo_dir = oa_dir[r] + or_Bo_dir;
            or_Bo_a   = oa_a[r]   + or_Bo_a;
            or_Bo_j   = oa_j[r]   + or_Bo_j;
            or_Bo_q   = oa_q[r]   + or_Bo_q;
          }
        }
      }

      double val = 0.0;
      if (mode_est == 0) {
        val = (direct_T - direct_C) + (a_T - a_C) + (j_T - j_C) + (sel_T - sel_C) + (-trt_C_dir) + (-trt_C_a);
      } else if (mode_est == 1) {
        val = (direct_T - direct_Bg) + (a_T - a_Bg) + j_T + sel_T + (-or_Bg_dir) + (-or_Bg_a) + (-or_Bg_j) + (-or_Bg_q);
      } else {
        val = ((direct_T - direct_Bg) - (direct_C - direct_Bo)) +
              ((a_T - a_Bg) - (a_C - a_Bo)) +
              (j_T - j_C) +
              (sel_T - sel_C) +
              (-(trt_C_dir - trt_Bo_dir)) +
              (-(trt_C_a - trt_Bo_a)) +
              (-or_Bg_dir + or_Bo_dir) +
              (-or_Bg_a + or_Bo_a) +
              (-or_Bg_j + or_Bo_j) +
              (-or_Bg_q + or_Bo_q);
      }
      if (!R_FINITE(val)) {
        UNPROTECT(2);
        vmaxset(vmax);
        return R_NilValue;
      }
      phi_col[r] = is_core_scale ? (lambda * val) : val;
    }
  }

  SEXP dn_s = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn_s, 0, R_NilValue);
  SET_VECTOR_ELT(dn_s, 1, link_names_s);
  Rf_setAttrib(phi_s, R_DimNamesSymbol, dn_s);

  UNPROTECT(3);
  vmaxset(vmax);
  return phi_s;
}


/* ===========================================================================
   6. Fast C Kernels for Imputation & Augmented Fitting Stages
   =========================================================================== */

static int cdid_ia_wls_compress_buf(double *A, int m, int k, double *fac_out) {
  if (m <= 0 || k <= 0) return 1;
  size_t tot = (size_t) m * (size_t) k;
  for (size_t idx = 0; idx < tot; ++idx) {
    if (!R_FINITE(A[idx])) return 0;
  }
  int jpvt_stack[64];
  double vn_stack[64];
  int *jpvt = (k <= 64) ? jpvt_stack : (int *) R_alloc((size_t) k, sizeof(int));
  double *vn = (k <= 64) ? vn_stack : (double *) R_alloc((size_t) k, sizeof(double));

  for (int j = 0; j < k; ++j) {
    jpvt[j] = j;
    const double *col = A + (size_t) j * (size_t) m;
    double s = 0.0;
    for (int i = 0; i < m; ++i) s += col[i] * col[i];
    vn[j] = s;
  }

  int r = (m < k) ? m : k;
  for (int step = 0; step < r; ++step) {
    int pvt = step;
    double best = vn[step];
    for (int j = step + 1; j < k; ++j) {
      if (vn[j] > best) {
        best = vn[j];
        pvt = j;
      }
    }
    if (pvt != step) {
      double *c1 = A + (size_t) step * (size_t) m;
      double *c2 = A + (size_t) pvt * (size_t) m;
      for (int i = 0; i < m; ++i) {
        double tmp = c1[i];
        c1[i] = c2[i];
        c2[i] = tmp;
      }
      int ip = jpvt[step];
      jpvt[step] = jpvt[pvt];
      jpvt[pvt] = ip;
      vn[pvt] = vn[step];
    }

    double *col_s = A + (size_t) step * (size_t) m;
    double alpha = col_s[step];
    if (step + 1 < m) {
      double scale = 0.0, ssq = 1.0;
      for (int i = step + 1; i < m; ++i) {
        double ax = fabs(col_s[i]);
        if (ax > 0.0) {
          if (scale < ax) {
            double ratio = scale / ax;
            ssq = 1.0 + ssq * ratio * ratio;
            scale = ax;
          } else {
            double ratio = ax / scale;
            ssq += ratio * ratio;
          }
        }
      }
      double xnorm = scale * sqrt(ssq);
      if (xnorm > 0.0) {
        double norm_full = hypot(alpha, xnorm);
        double beta = (alpha >= 0.0) ? -norm_full : norm_full;
        double tau = (beta - alpha) / beta;
        double inv_denom = 1.0 / (alpha - beta);
        for (int i = step + 1; i < m; ++i) col_s[i] *= inv_denom;
        col_s[step] = beta;

        for (int j = step + 1; j < k; ++j) {
          double *col_j = A + (size_t) j * (size_t) m;
          double dot = col_j[step];
          for (int i = step + 1; i < m; ++i) dot += col_s[i] * col_j[i];
          double tdot = tau * dot;
          col_j[step] -= tdot;
          double ns = 0.0;
          for (int i = step + 1; i < m; ++i) {
            double v = col_j[i] - tdot * col_s[i];
            col_j[i] = v;
            ns += v * v;
          }
          vn[j] = ns;
        }
      } else {
        for (int j = step + 1; j < k; ++j) {
          double *col_j = A + (size_t) j * (size_t) m;
          double ns = 0.0;
          for (int i = step + 1; i < m; ++i) ns += col_j[i] * col_j[i];
          vn[j] = ns;
        }
      }
    }
  }

  memset(fac_out, 0, (size_t) r * (size_t) k * sizeof(double));
  for (int col = 0; col < k; ++col) {
    int orig_col = jpvt[col];
    double *dst = fac_out + (size_t) orig_col * (size_t) r;
    const double *sc = A + (size_t) col * (size_t) m;
    int rlim = (col + 1 < r) ? (col + 1) : r;
    for (int i = 0; i < rlim; ++i) {
      if (!R_FINITE(sc[i])) return 0;
      dst[i] = sc[i];
    }
  }
  return 1;
}

typedef struct {
  int p;
  int icpt;
  int ncols;
  double L_val;
  int acc_n;
  int acc_nrow;
  int max_chunk_merge_rows;
  int n_compressions;
  double weight_mass;
  double sumY;
  double *sumZ;
  double *colmax;
  double *anchor;
  int *constant;
  double *acc_fac;
} cdid_ia_sup_acc_t;

static void cdid_ia_sup_acc_init(cdid_ia_sup_acc_t *acc, int p, int icpt, double L_val) {
  acc->p = p;
  acc->icpt = icpt ? 1 : 0;
  acc->ncols = p + 1 + acc->icpt;
  acc->L_val = L_val;
  acc->acc_n = 0;
  acc->acc_nrow = 0;
  acc->max_chunk_merge_rows = 0;
  acc->n_compressions = 0;
  acc->weight_mass = 0.0;
  acc->sumY = 0.0;
  acc->sumZ = (p > 0) ? (double *) R_alloc((size_t) p, sizeof(double)) : NULL;
  acc->colmax = (p > 0) ? (double *) R_alloc((size_t) p, sizeof(double)) : NULL;
  acc->anchor = (p > 0) ? (double *) R_alloc((size_t) p, sizeof(double)) : NULL;
  acc->constant = (p > 0) ? (int *) R_alloc((size_t) p, sizeof(int)) : NULL;
  acc->acc_fac = (double *) R_alloc((size_t) acc->ncols * (size_t) acc->ncols, sizeof(double));
  for (int j = 0; j < p; ++j) {
    acc->sumZ[j] = 0.0;
    acc->colmax[j] = 0.0;
    acc->anchor[j] = 0.0;
    acc->constant[j] = 1;
  }
}

static int cdid_ia_sup_acc_add_chunk(
    cdid_ia_sup_acc_t *acc,
    int n_obs,
    const int *obs_rows0,
    const double *dy_obs,
    const double *w_obs,
    const double *X_F,
    int N,
    int intercept_col0,
    const int *slope_cols0,
    double *work_A,
    double *chunk_fac,
    double *stack_A,
    double *new_fac) {
  if (n_obs <= 0) return 1;
  int p = acc->p;
  int icpt = acc->icpt;
  int ncols = acc->ncols;
  double L_val = acc->L_val;

  double chunk_w = 0.0;
  double chunk_wy = 0.0;
  double chunk_anchor[64];
  int chunk_const[64];
  for (int j = 0; j < p; ++j) chunk_const[j] = 1;

  for (int i = 0; i < n_obs; ++i) {
    int r0 = obs_rows0[i];
    double w = w_obs[i];
    double y = dy_obs[i];
    if (!R_FINITE(w) || w <= 0.0 || !R_FINITE(y)) return 0;
    if (icpt && intercept_col0 >= 0) {
      if (X_F[r0 + (size_t) intercept_col0 * (size_t) N] != 1.0) return 0;
    }
    double sw = sqrt(w);
    if (icpt) {
      work_A[i] = sw;
    }
    for (int j = 0; j < p; ++j) {
      double z = X_F[r0 + (size_t) slope_cols0[j] * (size_t) N];
      if (!R_FINITE(z)) return 0;
      double d = L_val * z;
      if (!R_FINITE(d)) return 0;
      if (i == 0) {
        chunk_anchor[j] = d;
      } else if (d != chunk_anchor[j]) {
        chunk_const[j] = 0;
      }
      work_A[i + (size_t) (icpt + j) * (size_t) n_obs] = d * sw;
      double az = fabs(z);
      if (az > acc->colmax[j]) acc->colmax[j] = az;
    }
    work_A[i + (size_t) (icpt + p) * (size_t) n_obs] = y * sw;
    chunk_w += w;
    chunk_wy += w * y;
  }

  for (int j = 0; j < p; ++j) {
    double swz = 0.0;
    int col0 = slope_cols0[j];
    for (int i = 0; i < n_obs; ++i) {
      int r0 = obs_rows0[i];
      swz += w_obs[i] * X_F[r0 + (size_t) col0 * (size_t) N];
    }
    acc->sumZ[j] += swz;
  }

  int r_chunk = (n_obs < ncols) ? n_obs : ncols;
  if (!cdid_ia_wls_compress_buf(work_A, n_obs, ncols, chunk_fac)) return 0;

  if (icpt) {
    int stacked_rows = acc->acc_nrow + r_chunk;
    if (stacked_rows > acc->max_chunk_merge_rows) {
      acc->max_chunk_merge_rows = stacked_rows;
    }
    if (acc->acc_n == 0) {
      for (int j = 0; j < p; ++j) {
        acc->anchor[j] = chunk_anchor[j];
        acc->constant[j] = chunk_const[j];
      }
    } else {
      for (int j = 0; j < p; ++j) {
        acc->constant[j] = acc->constant[j] && chunk_const[j] && (acc->anchor[j] == chunk_anchor[j]);
      }
    }
    acc->n_compressions += 1 + ((stacked_rows > 0) ? 1 : 0);
  }

  int m_stack = acc->acc_nrow + r_chunk;
  int new_nrow = (m_stack < ncols) ? m_stack : ncols;
  for (int col = 0; col < ncols; ++col) {
    double *dst = stack_A + (size_t) col * (size_t) m_stack;
    for (int i = 0; i < acc->acc_nrow; ++i) {
      dst[i] = acc->acc_fac[i + (size_t) col * (size_t) acc->acc_nrow];
    }
    for (int i = 0; i < r_chunk; ++i) {
      dst[acc->acc_nrow + i] = chunk_fac[i + (size_t) col * (size_t) r_chunk];
    }
  }
  if (!cdid_ia_wls_compress_buf(stack_A, m_stack, ncols, new_fac)) return 0;
  memcpy(acc->acc_fac, new_fac, (size_t) new_nrow * (size_t) ncols * sizeof(double));
  acc->acc_nrow = new_nrow;
  acc->acc_n += n_obs;
  acc->weight_mass += chunk_w;
  acc->sumY += chunk_wy;
  return 1;
}

static SEXP cdid_ia_sup_acc_to_sexp(const cdid_ia_sup_acc_t *acc, SEXP slope_columns_s) {
  int p = acc->p;
  int icpt = acc->icpt;
  int ncols = acc->ncols;
  int nr = acc->acc_nrow;

  SEXP fac_s = PROTECT(Rf_allocMatrix(REALSXP, nr, ncols));
  if (nr > 0) {
    memcpy(REAL(fac_s), acc->acc_fac, (size_t) nr * (size_t) ncols * sizeof(double));
  }
  SEXP raw_s = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP raw_nms = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(raw_nms, 0, Rf_mkChar("factor"));
  SET_STRING_ELT(raw_nms, 1, Rf_mkChar("n"));
  SET_STRING_ELT(raw_nms, 2, Rf_mkChar("columns"));
  Rf_setAttrib(raw_s, R_NamesSymbol, raw_nms);
  SET_VECTOR_ELT(raw_s, 0, fac_s);
  SET_VECTOR_ELT(raw_s, 1, Rf_ScalarInteger(acc->acc_n));
  SET_VECTOR_ELT(raw_s, 2, icpt ? R_NilValue : ((p > 0) ? slope_columns_s : R_NilValue));
  SEXP raw_cls = PROTECT(Rf_mkString("cdid_or_wls_summary"));
  Rf_setAttrib(raw_s, R_ClassSymbol, raw_cls);

  SEXP summary_s = R_NilValue;
  int nprot = 4;
  if (icpt) {
    summary_s = PROTECT(Rf_allocVector(VECSXP, 7));
    SEXP sum_nms = PROTECT(Rf_allocVector(STRSXP, 7));
    SET_STRING_ELT(sum_nms, 0, Rf_mkChar("raw"));
    SET_STRING_ELT(sum_nms, 1, Rf_mkChar("n"));
    SET_STRING_ELT(sum_nms, 2, Rf_mkChar("p"));
    SET_STRING_ELT(sum_nms, 3, Rf_mkChar("columns"));
    SET_STRING_ELT(sum_nms, 4, Rf_mkChar("anchor"));
    SET_STRING_ELT(sum_nms, 5, Rf_mkChar("constant"));
    SET_STRING_ELT(sum_nms, 6, Rf_mkChar("accounting"));
    Rf_setAttrib(summary_s, R_NamesSymbol, sum_nms);

    int anch_len = (acc->acc_n > 0) ? p : 0;
    SEXP anchor_s = PROTECT(Rf_allocVector(REALSXP, anch_len));
    for (int j = 0; j < anch_len; ++j) REAL(anchor_s)[j] = acc->anchor[j];
    SEXP const_s = PROTECT(Rf_allocVector(LGLSXP, p));
    for (int j = 0; j < p; ++j) LOGICAL(const_s)[j] = acc->constant[j] ? TRUE : FALSE;

    SEXP acct_s = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP acct_nms = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(acct_nms, 0, Rf_mkChar("input_rows"));
    SET_STRING_ELT(acct_nms, 1, Rf_mkChar("augmented_cols"));
    SET_STRING_ELT(acct_nms, 2, Rf_mkChar("augmented_rows"));
    SET_STRING_ELT(acct_nms, 3, Rf_mkChar("max_chunk_merge_rows"));
    SET_STRING_ELT(acct_nms, 4, Rf_mkChar("n_compressions"));
    Rf_setAttrib(acct_s, R_NamesSymbol, acct_nms);
    SET_VECTOR_ELT(acct_s, 0, Rf_ScalarInteger(acc->acc_n));
    SET_VECTOR_ELT(acct_s, 1, Rf_ScalarInteger(p + 2));
    SET_VECTOR_ELT(acct_s, 2, Rf_ScalarInteger(nr));
    SET_VECTOR_ELT(acct_s, 3, Rf_ScalarInteger(acc->max_chunk_merge_rows));
    SET_VECTOR_ELT(acct_s, 4, Rf_ScalarInteger(acc->n_compressions));

    SET_VECTOR_ELT(summary_s, 0, raw_s);
    SET_VECTOR_ELT(summary_s, 1, Rf_ScalarInteger(acc->acc_n));
    SET_VECTOR_ELT(summary_s, 2, Rf_ScalarInteger(p));
    SET_VECTOR_ELT(summary_s, 3, (p > 0) ? slope_columns_s : R_NilValue);
    SET_VECTOR_ELT(summary_s, 4, anchor_s);
    SET_VECTOR_ELT(summary_s, 5, const_s);
    SET_VECTOR_ELT(summary_s, 6, acct_s);

    SEXP sum_cls = PROTECT(Rf_mkString("cdid_or_pooled_wls_chunk"));
    Rf_setAttrib(summary_s, R_ClassSymbol, sum_cls);
    nprot += 7;
  } else {
    summary_s = raw_s;
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP out_nms = PROTECT(Rf_allocVector(STRSXP, 7));
  SET_STRING_ELT(out_nms, 0, Rf_mkChar("L"));
  SET_STRING_ELT(out_nms, 1, Rf_mkChar("summary"));
  SET_STRING_ELT(out_nms, 2, Rf_mkChar("weight_mass"));
  SET_STRING_ELT(out_nms, 3, Rf_mkChar("n"));
  SET_STRING_ELT(out_nms, 4, Rf_mkChar("sumZ"));
  SET_STRING_ELT(out_nms, 5, Rf_mkChar("sumY"));
  SET_STRING_ELT(out_nms, 6, Rf_mkChar("colmax"));
  Rf_setAttrib(out_s, R_NamesSymbol, out_nms);

  SEXP sumZ_s = PROTECT(Rf_allocVector(REALSXP, p));
  SEXP colmax_s = PROTECT(Rf_allocVector(REALSXP, p));
  for (int j = 0; j < p; ++j) {
    REAL(sumZ_s)[j] = acc->sumZ[j];
    REAL(colmax_s)[j] = acc->colmax[j];
  }
  SET_VECTOR_ELT(out_s, 0, Rf_ScalarReal(acc->L_val));
  SET_VECTOR_ELT(out_s, 1, summary_s);
  SET_VECTOR_ELT(out_s, 2, Rf_ScalarReal(acc->weight_mass));
  SET_VECTOR_ELT(out_s, 3, Rf_ScalarInteger(acc->acc_n));
  SET_VECTOR_ELT(out_s, 4, sumZ_s);
  SET_VECTOR_ELT(out_s, 5, Rf_ScalarReal(acc->sumY));
  SET_VECTOR_ELT(out_s, 6, colmax_s);

  UNPROTECT(nprot + 4);
  return out_s;
}

/* Kernel 1: Pure Imputation Stage 1 (cohort_colmax_F + cell_supports) */
SEXP C_cdid_ia_imp_stage1(
    SEXP st_s,
    SEXP panel_s,
    SEXP X_F_s,
    SEXP cat_s,
    SEXP recipe_s,
    SEXP or_frame_index_s,
    SEXP B_chunk_s) {
  if (TYPEOF(st_s) != VECSXP || TYPEOF(panel_s) != VECSXP ||
      TYPEOF(X_F_s) != REALSXP || TYPEOF(cat_s) != VECSXP ||
      TYPEOF(recipe_s) != VECSXP) {
    return R_NilValue;
  }
  const void *vmax = vmaxget();

  int N = Rf_asInteger(cdid_ia_get_elt(st_s, "N"));
  int N_rows = Rf_nrows(X_F_s);
  const int *cohort_order = INTEGER(cdid_ia_get_elt(st_s, "cohort_order"));
  SEXP c_start_s = cdid_ia_get_elt(st_s, "cohort_start");
  SEXP c_end_s   = cdid_ia_get_elt(st_s, "cohort_end");
  const int *a_index = INTEGER(cdid_ia_get_elt(st_s, "a_index"));
  const double *w_all = REAL(cdid_ia_get_elt(st_s, "w"));
  const double *mass_F = REAL(cdid_ia_get_elt(st_s, "mass_F"));
  const double *a_vec  = REAL(cdid_ia_get_elt(st_s, "a"));
  int G = (int) XLENGTH(cdid_ia_get_elt(st_s, "cohorts"));
  int B_chunk = Rf_asInteger(B_chunk_s);
  if (B_chunk < 1) B_chunk = N_rows;

  const double *X_F = REAL(X_F_s);
  int p = Rf_asInteger(cdid_ia_get_elt(recipe_s, "p"));
  if (p > 60) { vmaxset(vmax); return R_NilValue; }
  int icpt = Rf_asLogical(cdid_ia_get_elt(recipe_s, "intercept")) == TRUE;
  int icpt_idx = Rf_asInteger(cdid_ia_get_elt(recipe_s, "intercept_index"));
  SEXP slope_idx_s = cdid_ia_get_elt(recipe_s, "slope_index");
  SEXP slope_cols_nms = cdid_ia_get_elt(recipe_s, "slope_columns");

  int n_or_fi = (int) XLENGTH(or_frame_index_s);
  int *or_fi = (int *) R_alloc((size_t) (n_or_fi > 0 ? n_or_fi : 1), sizeof(int));
  if (TYPEOF(or_frame_index_s) == INTSXP) {
    for (int i = 0; i < n_or_fi; ++i) or_fi[i] = INTEGER(or_frame_index_s)[i];
  } else {
    for (int i = 0; i < n_or_fi; ++i) or_fi[i] = (int) REAL(or_frame_index_s)[i];
  }
  int intercept_col0 = (icpt && icpt_idx >= 1 && icpt_idx <= n_or_fi) ? (or_fi[icpt_idx - 1] - 1) : -1;
  int *slope_cols0 = (p > 0) ? (int *) R_alloc((size_t) p, sizeof(int)) : NULL;
  if (p > 0) {
    if (TYPEOF(slope_idx_s) == INTSXP) {
      for (int j = 0; j < p; ++j) slope_cols0[j] = or_fi[INTEGER(slope_idx_s)[j] - 1] - 1;
    } else {
      for (int j = 0; j < p; ++j) slope_cols0[j] = or_fi[(int) REAL(slope_idx_s)[j] - 1] - 1;
    }
  }

  SEXP active_coh_s = cdid_ia_get_elt(cat_s, "active_cohorts");
  int n_act = (int) XLENGTH(active_coh_s);
  for (int i = 0; i < n_act; ++i) {
    int ci = (TYPEOF(active_coh_s) == INTSXP) ? INTEGER(active_coh_s)[i] : (int) REAL(active_coh_s)[i];
    if (ci < 1 || ci > G || !R_FINITE(mass_F[ci - 1]) || mass_F[ci - 1] <= 0.0 ||
        !R_FINITE(a_vec[ci - 1]) || a_vec[ci - 1] <= 0.0) {
      vmaxset(vmax);
      return R_NilValue;
    }
  }

  SEXP colmax_list_s = PROTECT(Rf_allocVector(VECSXP, G));
  for (int i = 0; i < n_act; ++i) {
    int ci = (TYPEOF(active_coh_s) == INTSXP) ? INTEGER(active_coh_s)[i] : (int) REAL(active_coh_s)[i];
    int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci - 1] : (int) REAL(c_start_s)[ci - 1];
    int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci - 1]   : (int) REAL(c_end_s)[ci - 1];
    SEXP cm_s = PROTECT(Rf_allocVector(REALSXP, p));
    double *cm = REAL(cm_s);
    for (int j = 0; j < p; ++j) cm[j] = 0.0;
    if (p > 0) {
      for (int idx = first; idx <= last; ++idx) {
        int r0 = cohort_order[idx - 1] - 1;
        for (int j = 0; j < p; ++j) {
          double v = fabs(X_F[r0 + (size_t) slope_cols0[j] * (size_t) N_rows]);
          if (v > cm[j]) cm[j] = v;
        }
      }
    }
    SET_VECTOR_ELT(colmax_list_s, ci - 1, cm_s);
    UNPROTECT(1);
  }

  SEXP all_times_s = cdid_ia_get_elt(panel_s, "all_times");
  int n_dates = (int) XLENGTH(all_times_s);
  const double *all_times = REAL(all_times_s);
  SEXP post_s = cdid_ia_get_elt(panel_s, "postings");
  const int *post_a = INTEGER(cdid_ia_get_elt(post_s, "post_a_index"));
  const int *post_obs = INTEGER(cdid_ia_get_elt(post_s, "post_obs_index"));
  const int *date_start = INTEGER(cdid_ia_get_elt(post_s, "date_start"));
  const int *date_end = INTEGER(cdid_ia_get_elt(post_s, "date_end"));
  const int *date_count = INTEGER(cdid_ia_get_elt(post_s, "date_count"));
  SEXP y_obs_s = cdid_ia_get_elt(cdid_ia_get_elt(panel_s, "observations"), "y");
  const double *y_obs = REAL(y_obs_s);
  R_xlen_t n_y = XLENGTH(y_obs_s);
  R_xlen_t n_post = XLENGTH(cdid_ia_get_elt(post_s, "post_a_index"));

  SEXP cell_specs_s = cdid_ia_get_elt(cat_s, "cell_specs");
  int n_cells = (int) XLENGTH(cell_specs_s);
  SEXP cell_supports_s = PROTECT(Rf_allocVector(VECSXP, n_cells));

  int ncols = p + 1 + (icpt ? 1 : 0);
  int *obs_rows0 = (int *) R_alloc((size_t) B_chunk, sizeof(int));
  double *dy_obs = (double *) R_alloc((size_t) B_chunk, sizeof(double));
  double *w_obs  = (double *) R_alloc((size_t) B_chunk, sizeof(double));
  double *work_A = (double *) R_alloc((size_t) B_chunk * (size_t) ncols, sizeof(double));
  double *chunk_fac = (double *) R_alloc((size_t) ncols * (size_t) ncols, sizeof(double));
  double *stack_A   = (double *) R_alloc((size_t) (2 * ncols) * (size_t) ncols, sizeof(double));
  double *new_fac   = (double *) R_alloc((size_t) ncols * (size_t) ncols, sizeof(double));

  for (int cid = 0; cid < n_cells; ++cid) {
    SEXP cs = VECTOR_ELT(cell_specs_s, cid);
    int ci = Rf_asInteger(cdid_ia_get_elt(cs, "cohort_index"));
    double L_val = Rf_asReal(cdid_ia_get_elt(cs, "L"));
    int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci - 1] : (int) REAL(c_start_s)[ci - 1];
    int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci - 1]   : (int) REAL(c_end_s)[ci - 1];

    int ds = 0, dt = 0;
    if (!cdid_ia_resolve_pair_dates(cdid_ia_get_elt(cs, "pair"), all_times, n_dates, &ds, &dt)) {
      UNPROTECT(2); vmaxset(vmax); return R_NilValue;
    }
    int has_pair = (date_count[ds - 1] > 0 && date_count[dt - 1] > 0);
    int s_lo = date_start[ds - 1] - 1, s_hi = date_end[ds - 1] - 1;
    int t_lo = date_start[dt - 1] - 1, t_hi = date_end[dt - 1] - 1;

    cdid_ia_sup_acc_t acc;
    cdid_ia_sup_acc_init(&acc, p, icpt, L_val);
    int hs_hint = s_lo, ht_hint = t_lo;

    for (int start = first; start <= last; start += B_chunk) {
      int end = start + B_chunk - 1;
      if (end > last) end = last;
      int n_obs = 0;
      if (has_pair) {
        if (s_lo < 0 || s_hi < s_lo || (R_xlen_t) s_hi >= n_post ||
            t_lo < 0 || t_hi < t_lo || (R_xlen_t) t_hi >= n_post) {
          UNPROTECT(2); vmaxset(vmax); return R_NilValue;
        }
        for (int idx = start; idx <= end; ++idx) {
          int r0 = cohort_order[idx - 1] - 1;
          int u = a_index[r0];
          if (u != NA_INTEGER && u > 0) {
            int hs = cdid_ia_find_posting(post_a, s_lo, s_hi, u, &hs_hint);
            if (hs >= 0) {
              int ht = cdid_ia_find_posting(post_a, t_lo, t_hi, u, &ht_hint);
              if (ht >= 0) {
                int obs_s = post_obs[hs], obs_t = post_obs[ht];
                if (obs_s < 1 || (R_xlen_t) obs_s > n_y || obs_t < 1 || (R_xlen_t) obs_t > n_y) {
                  UNPROTECT(2); vmaxset(vmax); return R_NilValue;
                }
                double dy = y_obs[obs_t - 1] - y_obs[obs_s - 1];
                if (!R_FINITE(dy)) { UNPROTECT(2); vmaxset(vmax); return R_NilValue; }
                obs_rows0[n_obs] = r0;
                dy_obs[n_obs] = dy;
                w_obs[n_obs] = w_all[r0];
                n_obs++;
              }
            }
          }
        }
      }
      if (n_obs > 0) {
        if (!cdid_ia_sup_acc_add_chunk(
              &acc, n_obs, obs_rows0, dy_obs, w_obs, X_F, N_rows,
              intercept_col0, slope_cols0, work_A, chunk_fac, stack_A, new_fac)) {
          UNPROTECT(2); vmaxset(vmax); return R_NilValue;
        }
      }
    }
    if (acc.acc_n < 1 || !R_FINITE(acc.weight_mass) || acc.weight_mass <= 0.0 ||
        (acc.weight_mass / (double) N) <= 0.0) {
      UNPROTECT(2); vmaxset(vmax); return R_NilValue;
    }
    SET_VECTOR_ELT(cell_supports_s, cid, cdid_ia_sup_acc_to_sexp(&acc, slope_cols_nms));
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP out_nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(out_nms, 0, Rf_mkChar("cohort_colmax_F"));
  SET_STRING_ELT(out_nms, 1, Rf_mkChar("cell_supports"));
  Rf_setAttrib(out_s, R_NamesSymbol, out_nms);
  SET_VECTOR_ELT(out_s, 0, colmax_list_s);
  SET_VECTOR_ELT(out_s, 1, cell_supports_s);

  UNPROTECT(4);
  vmaxset(vmax);
  return out_s;
}

/* Kernel 2: Pure Imputation Cell Moments + Stage 2 Control OR Support Summaries */
SEXP C_cdid_ia_imp_stage2_cells(
    SEXP st_s,
    SEXP X_F_s,
    SEXP cat_s,
    SEXP recipe_s,
    SEXP or_frame_index_s,
    SEXP B_chunk_s,
    SEXP cell_meta_s,
    SEXP used_ctrl_s) {
  if (TYPEOF(st_s) != VECSXP || TYPEOF(X_F_s) != REALSXP ||
      TYPEOF(cat_s) != VECSXP || TYPEOF(recipe_s) != VECSXP ||
      TYPEOF(cell_meta_s) != VECSXP) {
    return R_NilValue;
  }
  const void *vmax = vmaxget();

  int N = Rf_asInteger(cdid_ia_get_elt(st_s, "N"));
  int N_rows = Rf_nrows(X_F_s);
  const int *cohort_order = INTEGER(cdid_ia_get_elt(st_s, "cohort_order"));
  SEXP c_start_s = cdid_ia_get_elt(st_s, "cohort_start");
  SEXP c_end_s   = cdid_ia_get_elt(st_s, "cohort_end");
  const double *w_all = REAL(cdid_ia_get_elt(st_s, "w"));
  const double *a_vec = REAL(cdid_ia_get_elt(st_s, "a"));
  int B_chunk = Rf_asInteger(B_chunk_s);
  if (B_chunk < 1) B_chunk = N_rows;

  const double *X_F = REAL(X_F_s);
  int p = Rf_asInteger(cdid_ia_get_elt(recipe_s, "p"));
  if (p > 60) { vmaxset(vmax); return R_NilValue; }
  int icpt = Rf_asLogical(cdid_ia_get_elt(recipe_s, "intercept")) == TRUE;
  int icpt_idx = Rf_asInteger(cdid_ia_get_elt(recipe_s, "intercept_index"));
  SEXP slope_idx_s = cdid_ia_get_elt(recipe_s, "slope_index");
  SEXP slope_cols_nms = cdid_ia_get_elt(recipe_s, "slope_columns");

  int n_or_fi = (int) XLENGTH(or_frame_index_s);
  int *or_fi = (int *) R_alloc((size_t) (n_or_fi > 0 ? n_or_fi : 1), sizeof(int));
  if (TYPEOF(or_frame_index_s) == INTSXP) {
    for (int i = 0; i < n_or_fi; ++i) or_fi[i] = INTEGER(or_frame_index_s)[i];
  } else {
    for (int i = 0; i < n_or_fi; ++i) or_fi[i] = (int) REAL(or_frame_index_s)[i];
  }
  int intercept_col0 = (icpt && icpt_idx >= 1 && icpt_idx <= n_or_fi) ? (or_fi[icpt_idx - 1] - 1) : -1;
  int *slope_cols0 = (p > 0) ? (int *) R_alloc((size_t) p, sizeof(int)) : NULL;
  if (p > 0) {
    if (TYPEOF(slope_idx_s) == INTSXP) {
      for (int j = 0; j < p; ++j) slope_cols0[j] = or_fi[INTEGER(slope_idx_s)[j] - 1] - 1;
    } else {
      for (int j = 0; j < p; ++j) slope_cols0[j] = or_fi[(int) REAL(slope_idx_s)[j] - 1] - 1;
    }
  }

  SEXP cell_specs_s = cdid_ia_get_elt(cat_s, "cell_specs");
  int n_cells = (int) XLENGTH(cell_specs_s);
  int *need_or = (int *) R_alloc((size_t) (n_cells + 1), sizeof(int));
  memset(need_or, 0, (size_t) (n_cells + 1) * sizeof(int));
  if (used_ctrl_s != R_NilValue) {
    int n_u = (int) XLENGTH(used_ctrl_s);
    for (int i = 0; i < n_u; ++i) {
      int cid = (TYPEOF(used_ctrl_s) == INTSXP) ? INTEGER(used_ctrl_s)[i] : (int) REAL(used_ctrl_s)[i];
      if (cid >= 1 && cid <= n_cells) need_or[cid] = 1;
    }
  }

  SEXP moments_list_s = PROTECT(Rf_allocVector(VECSXP, n_cells));
  SEXP or_acc_list_s  = PROTECT(Rf_allocVector(VECSXP, n_cells));

  int ncols = p + 1 + (icpt ? 1 : 0);
  int *chunk_rows0 = (int *) R_alloc((size_t) B_chunk, sizeof(int));
  double *jhat_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
  double *wsyn_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
  double *work_A   = (double *) R_alloc((size_t) B_chunk * (size_t) ncols, sizeof(double));
  double *chunk_fac = (double *) R_alloc((size_t) ncols * (size_t) ncols, sizeof(double));
  double *stack_A   = (double *) R_alloc((size_t) (2 * ncols) * (size_t) ncols, sizeof(double));
  double *new_fac   = (double *) R_alloc((size_t) ncols * (size_t) ncols, sizeof(double));

  for (int cid = 0; cid < n_cells; ++cid) {
    SEXP cs = VECTOR_ELT(cell_specs_s, cid);
    SEXP cm = VECTOR_ELT(cell_meta_s, cid);
    int ci = Rf_asInteger(cdid_ia_get_elt(cs, "cohort_index"));
    double L_val = Rf_asReal(cdid_ia_get_elt(cs, "L"));
    double a_k = a_vec[ci - 1];
    double alpha_c = Rf_asReal(cdid_ia_get_elt(cm, "alpha"));
    SEXP kept_cols_F_s = cdid_ia_get_elt(cm, "kept_cols_F");
    SEXP beta_r_s      = cdid_ia_get_elt(cm, "beta");
    SEXP keep_r_s      = cdid_ia_get_elt(cm, "keep");
    int p_r = (int) XLENGTH(beta_r_s);
    const int *kept_F = (p_r > 0) ? INTEGER(kept_cols_F_s) : NULL;
    const double *beta_r = (p_r > 0) ? REAL(beta_r_s) : NULL;
    const int *keep_r = (p_r > 0) ? INTEGER(keep_r_s) : NULL;

    int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci - 1] : (int) REAL(c_start_s)[ci - 1];
    int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci - 1]   : (int) REAL(c_end_s)[ci - 1];

    double M_one = 0.0;
    double Y_one = 0.0;
    double M_z[64], Y_z[64];
    for (int j = 0; j < p; ++j) { M_z[j] = 0.0; Y_z[j] = 0.0; }

    int do_or = need_or[cid + 1];
    cdid_ia_sup_acc_t acc_or;
    if (do_or) cdid_ia_sup_acc_init(&acc_or, p, icpt, L_val);

    for (int start = first; start <= last; start += B_chunk) {
      int end = start + B_chunk - 1;
      if (end > last) end = last;
      int n_r = end - start + 1;

      double ch_M1 = 0.0, ch_Y1 = 0.0;
      for (int i = 0; i < n_r; ++i) {
        int r0 = cohort_order[start - 1 + i] - 1;
        chunk_rows0[i] = r0;
        double ws = w_all[r0] * a_k;
        wsyn_buf[i] = ws;
        double jh_sl = 0.0;
        for (int m = 0; m < p_r; ++m) {
          int col0 = kept_F[m] - 1;
          jh_sl += (L_val * X_F[r0 + (size_t) col0 * (size_t) N_rows]) * beta_r[m];
        }
        double jh = alpha_c + jh_sl;
        jhat_buf[i] = jh;
        ch_M1 += ws;
        ch_Y1 += ws * jh;
      }
      M_one += ch_M1;
      Y_one += ch_Y1;

      for (int j = 0; j < p; ++j) {
        int col0 = slope_cols0[j];
        double ch_Mz = 0.0, ch_Yz = 0.0;
        for (int i = 0; i < n_r; ++i) {
          int r0 = chunk_rows0[i];
          double z = X_F[r0 + (size_t) col0 * (size_t) N_rows];
          double ws = wsyn_buf[i];
          ch_Mz += ws * z;
          ch_Yz += (ws * jhat_buf[i]) * z;
        }
        M_z[j] += ch_Mz;
        Y_z[j] += ch_Yz;
      }

      if (do_or) {
        if (!cdid_ia_sup_acc_add_chunk(
              &acc_or, n_r, chunk_rows0, jhat_buf, wsyn_buf, X_F, N_rows,
              intercept_col0, slope_cols0, work_A, chunk_fac, stack_A, new_fac)) {
          UNPROTECT(2); vmaxset(vmax); return R_NilValue;
        }
      }
    }

    M_one /= (double) N;
    Y_one /= (double) N;
    if (!R_FINITE(M_one) || !R_FINITE(Y_one)) {
      UNPROTECT(2); vmaxset(vmax); return R_NilValue;
    }

    SEXP M_z_s = PROTECT(Rf_allocVector(REALSXP, p));
    SEXP Y_z_s = PROTECT(Rf_allocVector(REALSXP, p));
    SEXP M_d_kept_s = PROTECT(Rf_allocVector(REALSXP, p_r));
    for (int j = 0; j < p; ++j) {
      double mz = M_z[j] / (double) N;
      double yz = Y_z[j] / (double) N;
      if (!R_FINITE(mz) || !R_FINITE(yz)) {
        UNPROTECT(5); vmaxset(vmax); return R_NilValue;
      }
      REAL(M_z_s)[j] = mz;
      REAL(Y_z_s)[j] = yz;
    }
    if (TYPEOF(slope_cols_nms) == STRSXP) {
      Rf_setAttrib(M_z_s, R_NamesSymbol, slope_cols_nms);
      Rf_setAttrib(Y_z_s, R_NamesSymbol, slope_cols_nms);
    }
    for (int m = 0; m < p_r; ++m) {
      REAL(M_d_kept_s)[m] = L_val * REAL(M_z_s)[keep_r[m] - 1];
    }
    if (p_r > 0 && TYPEOF(slope_cols_nms) == STRSXP) {
      SEXP kn_s = PROTECT(Rf_allocVector(STRSXP, p_r));
      for (int m = 0; m < p_r; ++m) {
        SET_STRING_ELT(kn_s, m, STRING_ELT(slope_cols_nms, keep_r[m] - 1));
      }
      Rf_setAttrib(M_d_kept_s, R_NamesSymbol, kn_s);
      UNPROTECT(1);
    }

    SEXP mom_s = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP mom_nms = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(mom_nms, 0, Rf_mkChar("M_one"));
    SET_STRING_ELT(mom_nms, 1, Rf_mkChar("Y_one"));
    SET_STRING_ELT(mom_nms, 2, Rf_mkChar("M_z"));
    SET_STRING_ELT(mom_nms, 3, Rf_mkChar("Y_z"));
    SET_STRING_ELT(mom_nms, 4, Rf_mkChar("M_d_kept"));
    Rf_setAttrib(mom_s, R_NamesSymbol, mom_nms);
    SET_VECTOR_ELT(mom_s, 0, Rf_ScalarReal(M_one));
    SET_VECTOR_ELT(mom_s, 1, Rf_ScalarReal(Y_one));
    SET_VECTOR_ELT(mom_s, 2, M_z_s);
    SET_VECTOR_ELT(mom_s, 3, Y_z_s);
    SET_VECTOR_ELT(mom_s, 4, M_d_kept_s);
    SET_VECTOR_ELT(moments_list_s, cid, mom_s);
    UNPROTECT(5);

    if (do_or) {
      SET_VECTOR_ELT(or_acc_list_s, cid, cdid_ia_sup_acc_to_sexp(&acc_or, slope_cols_nms));
    }
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP out_nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(out_nms, 0, Rf_mkChar("moments"));
  SET_STRING_ELT(out_nms, 1, Rf_mkChar("ctrl_cell_or_acc"));
  Rf_setAttrib(out_s, R_NamesSymbol, out_nms);
  SET_VECTOR_ELT(out_s, 0, moments_list_s);
  SET_VECTOR_ELT(out_s, 1, or_acc_list_s);

  UNPROTECT(4);
  vmaxset(vmax);
  return out_s;
}


/* Kernel 3: Pure Imputation Stage 2 OR Block Forcing + Stage 4 Link Control Scans */
SEXP C_cdid_ia_imp_stage24_forcing_links(
    SEXP st_s,
    SEXP X_F_s,
    SEXP cat_s,
    SEXP recipe_s,
    SEXP trt_recipe_s,
    SEXP B_chunk_s,
    SEXP cells_fitted_s,
    SEXP group_fits_s,
    SEXP extra_s) {
  if (TYPEOF(st_s) != VECSXP || TYPEOF(X_F_s) != REALSXP ||
      TYPEOF(cat_s) != VECSXP || TYPEOF(recipe_s) != VECSXP ||
      TYPEOF(cells_fitted_s) != VECSXP || TYPEOF(group_fits_s) != VECSXP ||
      TYPEOF(extra_s) != VECSXP) {
    return R_NilValue;
  }
  const void *vmax = vmaxget();

  int N = Rf_asInteger(cdid_ia_get_elt(st_s, "N"));
  int N_rows = Rf_nrows(X_F_s);
  const int *cohort_order = INTEGER(cdid_ia_get_elt(st_s, "cohort_order"));
  SEXP c_start_s = cdid_ia_get_elt(st_s, "cohort_start");
  SEXP c_end_s   = cdid_ia_get_elt(st_s, "cohort_end");
  const double *w_all = REAL(cdid_ia_get_elt(st_s, "w"));
  const double *a_vec = REAL(cdid_ia_get_elt(st_s, "a"));
  int B_chunk = Rf_asInteger(B_chunk_s);
  if (B_chunk < 1) B_chunk = N_rows;

  const double *X_F = REAL(X_F_s);
  int icpt = Rf_asLogical(cdid_ia_get_elt(recipe_s, "intercept")) == TRUE;

  SEXP block_specs_s = cdid_ia_get_elt(extra_s, "block_specs");
  SEXP trt_fits_s    = cdid_ia_get_elt(extra_s, "trt_fits");
  SEXP do_bf_s       = cdid_ia_get_elt(extra_s, "do_block_forcing");
  int do_bf = (do_bf_s == R_NilValue) ? 1 : (Rf_asLogical(do_bf_s) == TRUE);

  /* Part A: Stage 2 OR Block Forcing */
  SEXP block_forcing_s = R_NilValue;
  int nprot = 0;
  if (do_bf && block_specs_s != R_NilValue && TYPEOF(block_specs_s) == VECSXP) {
    SEXP or_blocks_s = cdid_ia_get_elt(cat_s, "or_blocks");
    int n_blocks = (int) XLENGTH(or_blocks_s);
    block_forcing_s = PROTECT(Rf_allocVector(VECSXP, n_blocks));
    nprot++;

    double *r_or_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *w_raw_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *w_syn_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *d_tilde_buf = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));
    double *d_imp_buf   = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));

    for (int bid = 0; bid < n_blocks; ++bid) {
      SEXP b_spec = VECTOR_ELT(or_blocks_s, bid);
      SEXP bs = VECTOR_ELT(block_specs_s, bid);
      double L_b = Rf_asReal(cdid_ia_get_elt(b_spec, "L"));
      SEXP cids_b_s = cdid_ia_get_elt(b_spec, "control_cell_ids");
      int n_kb = (int) XLENGTH(cids_b_s);

      SEXP kept_or_F_s = cdid_ia_get_elt(bs, "kept_or_cols_F");
      SEXP beta_or_s   = cdid_ia_get_elt(bs, "beta_or");
      double alpha_b   = Rf_asReal(cdid_ia_get_elt(bs, "alpha_b"));
      SEXP mean_d_b_s  = cdid_ia_get_elt(bs, "mean_d_b");
      int p_or = (int) XLENGTH(beta_or_s);
      if (p_or > 60) { UNPROTECT(nprot); vmaxset(vmax); return R_NilValue; }
      const int *kept_or_F = (p_or > 0) ? INTEGER(kept_or_F_s) : NULL;
      const double *beta_or = (p_or > 0) ? REAL(beta_or_s) : NULL;
      const double *mean_d_b = (p_or > 0) ? REAL(mean_d_b_s) : NULL;

      SEXP a_fa_s = PROTECT(Rf_allocVector(REALSXP, n_kb));
      SEXP a_fb_s = PROTECT(Rf_allocMatrix(REALSXP, p_or, n_kb));
      SEXP j_faa_s = PROTECT(Rf_allocVector(REALSXP, n_kb));
      SEXP j_fab_s = PROTECT(Rf_allocVector(VECSXP, n_kb));
      SEXP j_fba_s = PROTECT(Rf_allocMatrix(REALSXP, p_or, n_kb));
      SEXP j_fbb_s = PROTECT(Rf_allocVector(VECSXP, n_kb));

      double *a_fa = REAL(a_fa_s);
      double *a_fb = REAL(a_fb_s);
      double *j_faa = REAL(j_faa_s);
      double *j_fba = REAL(j_fba_s);
      memset(a_fb, 0, (size_t) p_or * (size_t) n_kb * sizeof(double));
      memset(j_fba, 0, (size_t) p_or * (size_t) n_kb * sizeof(double));

      for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
        int cid = (TYPEOF(cids_b_s) == INTSXP) ? INTEGER(cids_b_s)[idx_k] : (int) REAL(cids_b_s)[idx_k];
        SEXP cf = VECTOR_ELT(cells_fitted_s, cid - 1);
        int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
        SEXP gf = VECTOR_ELT(group_fits_s, gid - 1);
        int ci = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
        double a_k = a_vec[ci - 1];
        double alpha_c = Rf_asReal(cdid_ia_get_elt(cf, "alpha"));
        int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
        if (pg > 60) { UNPROTECT(nprot + 6); vmaxset(vmax); return R_NilValue; }
        const int *kept_imp_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf, "kept_cols_F")) : NULL;
        const double *beta_imp = (pg > 0) ? REAL(cdid_ia_get_elt(gf, "beta")) : NULL;

        int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci - 1] : (int) REAL(c_start_s)[ci - 1];
        int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci - 1]   : (int) REAL(c_end_s)[ci - 1];

        double sum_w_r = 0.0;
        double sum_w_r_dt[64];
        double sum_wsyn_dt[64];
        for (int m = 0; m < p_or; ++m) { sum_w_r_dt[m] = 0.0; sum_wsyn_dt[m] = 0.0; }

        SEXP fbb_k_s = PROTECT(Rf_allocMatrix(REALSXP, p_or, pg));
        double *fbb_k = REAL(fbb_k_s);
        memset(fbb_k, 0, (size_t) p_or * (size_t) pg * sizeof(double));

        for (int start = first; start <= last; start += B_chunk) {
          int end = start + B_chunk - 1;
          if (end > last) end = last;
          int n_r = end - start + 1;

          double ch_wr = 0.0;
          for (int i = 0; i < n_r; ++i) {
            int r0 = cohort_order[start - 1 + i] - 1;
            double wr = w_all[r0];
            double ws = wr * a_k;
            w_raw_buf[i] = wr;
            w_syn_buf[i] = ws;
            double jh_sl = 0.0;
            for (int m = 0; m < pg; ++m) {
              double di = L_b * X_F[r0 + (size_t) (kept_imp_F[m] - 1) * (size_t) N_rows];
              d_imp_buf[i + (size_t) m * (size_t) n_r] = di;
              jh_sl += di * beta_imp[m];
            }
            double jh = alpha_c + jh_sl;
            double mh_sl = 0.0;
            for (int m = 0; m < p_or; ++m) {
              double dor = L_b * X_F[r0 + (size_t) (kept_or_F[m] - 1) * (size_t) N_rows];
              double dt = icpt ? (dor - mean_d_b[m]) : dor;
              d_tilde_buf[i + (size_t) m * (size_t) n_r] = dt;
              mh_sl += dor * beta_or[m];
            }
            double mh = alpha_b + mh_sl;
            double r_or = jh - mh;
            r_or_buf[i] = r_or;
            ch_wr += wr * r_or;
          }
          sum_w_r += ch_wr;

          if (p_or > 0) {
            for (int m = 0; m < p_or; ++m) {
              const double *dt_col = d_tilde_buf + (size_t) m * (size_t) n_r;
              double ch_wrdt = 0.0, ch_wsdt = 0.0;
              for (int i = 0; i < n_r; ++i) {
                ch_wrdt += (w_raw_buf[i] * r_or_buf[i]) * dt_col[i];
                ch_wsdt += w_syn_buf[i] * dt_col[i];
              }
              sum_w_r_dt[m] += ch_wrdt;
              sum_wsyn_dt[m] += ch_wsdt;
              for (int g = 0; g < pg; ++g) {
                const double *di_col = d_imp_buf + (size_t) g * (size_t) n_r;
                double ch_dd = 0.0;
                for (int i = 0; i < n_r; ++i) {
                  ch_dd += (w_syn_buf[i] * dt_col[i]) * di_col[i];
                }
                fbb_k[m + (size_t) g * (size_t) p_or] += ch_dd;
              }
            }
          }
        }

        a_fa[idx_k] = sum_w_r / (double) N;
        for (int m = 0; m < p_or; ++m) {
          a_fb[m + (size_t) idx_k * (size_t) p_or] = sum_w_r_dt[m] / (double) N;
          j_fba[m + (size_t) idx_k * (size_t) p_or] = icpt ? (sum_wsyn_dt[m] / (double) N) : 0.0;
        }
        j_faa[idx_k] = icpt ? Rf_asReal(cdid_ia_get_elt(cf, "M_one")) : 0.0;

        SEXP mdk_s = cdid_ia_get_elt(cf, "M_d_kept");
        SEXP fab_k_s = PROTECT(Rf_allocVector(REALSXP, pg));
        for (int g = 0; g < pg; ++g) REAL(fab_k_s)[g] = REAL(mdk_s)[g];
        SET_VECTOR_ELT(j_fab_s, idx_k, fab_k_s);

        for (int idx = 0; idx < p_or * pg; ++idx) fbb_k[idx] /= (double) N;
        SEXP kc_or_s = cdid_ia_get_elt(bs, "kept_or_names");
        SEXP kc_imp_s = cdid_ia_get_elt(gf, "kept_columns");
        if (p_or > 0 && pg > 0 && TYPEOF(kc_or_s) == STRSXP && TYPEOF(kc_imp_s) == STRSXP) {
          SEXP dn_s = PROTECT(Rf_allocVector(VECSXP, 2));
          SET_VECTOR_ELT(dn_s, 0, kc_or_s);
          SET_VECTOR_ELT(dn_s, 1, kc_imp_s);
          Rf_setAttrib(fbb_k_s, R_DimNamesSymbol, dn_s);
          UNPROTECT(1);
        }
        SET_VECTOR_ELT(j_fbb_s, idx_k, fbb_k_s);
        UNPROTECT(2);
      }

      SEXP bf_s = PROTECT(Rf_allocVector(VECSXP, 6));
      SEXP bf_nms = PROTECT(Rf_allocVector(STRSXP, 6));
      SET_STRING_ELT(bf_nms, 0, Rf_mkChar("a_force_alpha"));
      SET_STRING_ELT(bf_nms, 1, Rf_mkChar("a_force_beta"));
      SET_STRING_ELT(bf_nms, 2, Rf_mkChar("j_force_alpha_alpha"));
      SET_STRING_ELT(bf_nms, 3, Rf_mkChar("j_force_alpha_beta"));
      SET_STRING_ELT(bf_nms, 4, Rf_mkChar("j_force_beta_alpha"));
      SET_STRING_ELT(bf_nms, 5, Rf_mkChar("j_force_beta_beta"));
      Rf_setAttrib(bf_s, R_NamesSymbol, bf_nms);
      SET_VECTOR_ELT(bf_s, 0, a_fa_s);
      SET_VECTOR_ELT(bf_s, 1, a_fb_s);
      SET_VECTOR_ELT(bf_s, 2, j_faa_s);
      SET_VECTOR_ELT(bf_s, 3, j_fab_s);
      SET_VECTOR_ELT(bf_s, 4, j_fba_s);
      SET_VECTOR_ELT(bf_s, 5, j_fbb_s);
      SET_VECTOR_ELT(block_forcing_s, bid, bf_s);
      UNPROTECT(8);
    }
  }

  /* Part B: Stage 4 Link Control Scans */
  SEXP link_ctrl_s = R_NilValue;
  if (trt_fits_s != R_NilValue && TYPEOF(trt_fits_s) == VECSXP) {
    SEXP links_s = cdid_ia_get_elt(cat_s, "links");
    int L_links = (int) XLENGTH(links_s);
    int p_trt = Rf_asInteger(cdid_ia_get_elt(trt_recipe_s, "p"));
    if (p_trt > 60) { UNPROTECT(nprot); vmaxset(vmax); return R_NilValue; }
    const int *trt_cols_F = (p_trt > 0) ? INTEGER(cdid_ia_get_elt(trt_recipe_s, "index")) : NULL;

    link_ctrl_s = PROTECT(Rf_allocVector(VECSXP, L_links));
    nprot++;

    double *wm_buf   = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *wdot_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *wodd_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *jhat_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *mhat_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *dimp_buf = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));
    double *dor_buf  = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));
    double *xtrt_buf = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));

    for (int ell = 0; ell < L_links; ++ell) {
      SEXP lk = VECTOR_ELT(links_s, ell);
      int trt_id = Rf_asInteger(cdid_ia_get_elt(lk, "trt_id"));
      SEXP tfit = VECTOR_ELT(trt_fits_s, trt_id - 1);
      const double *gamma_t = (p_trt > 0) ? REAL(cdid_ia_get_elt(tfit, "coefficients")) : NULL;
      double L_lk = Rf_asReal(cdid_ia_get_elt(lk, "L"));

      SEXP ctrl_coh_s = cdid_ia_get_elt(lk, "control_cohort_indices");
      SEXP ctrl_cid_s = cdid_ia_get_elt(lk, "control_cell_ids");
      int n_ctrl = (int) XLENGTH(ctrl_coh_s);

      int has_or = (block_specs_s != R_NilValue);
      int p_or = 0;
      double alpha_b = 0.0;
      const int *kept_or_F = NULL;
      const double *beta_or = NULL;
      if (has_or) {
        int bid = Rf_asInteger(cdid_ia_get_elt(lk, "or_block_id"));
        SEXP bs = VECTOR_ELT(block_specs_s, bid - 1);
        SEXP beta_or_s = cdid_ia_get_elt(bs, "beta_or");
        p_or = (int) XLENGTH(beta_or_s);
        kept_or_F = (p_or > 0) ? INTEGER(cdid_ia_get_elt(bs, "kept_or_cols_F")) : NULL;
        beta_or = (p_or > 0) ? REAL(beta_or_s) : NULL;
        alpha_b = icpt ? Rf_asReal(cdid_ia_get_elt(bs, "alpha_b")) : 0.0;
      }

      double sum_wm = 0.0, sum_wm_j = 0.0, sum_wm_m = 0.0;
      double sum_wm_dor[64];
      for (int m = 0; m < p_or; ++m) sum_wm_dor[m] = 0.0;
      double sum_wdot_j_X[64], sum_wdot_X[64], sum_wdot_m_X[64];
      for (int m = 0; m < p_trt; ++m) {
        sum_wdot_j_X[m] = 0.0;
        sum_wdot_X[m] = 0.0;
        sum_wdot_m_X[m] = 0.0;
      }

      SEXP wm_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
      SEXP dimp_k_s = PROTECT(Rf_allocVector(VECSXP, n_ctrl));
      SEXP wodd_j_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
      SEXP wodd_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
      SEXP wodd_m_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));

      double *wm_k = REAL(wm_k_s);
      double *wodd_j_k = REAL(wodd_j_k_s);
      double *wodd_k = REAL(wodd_k_s);
      double *wodd_m_k = REAL(wodd_m_k_s);

      for (int idx_k = 0; idx_k < n_ctrl; ++idx_k) {
        int ci_k = (TYPEOF(ctrl_coh_s) == INTSXP) ? INTEGER(ctrl_coh_s)[idx_k] : (int) REAL(ctrl_coh_s)[idx_k];
        int cid_k = (TYPEOF(ctrl_cid_s) == INTSXP) ? INTEGER(ctrl_cid_s)[idx_k] : (int) REAL(ctrl_cid_s)[idx_k];
        SEXP cf_k = VECTOR_ELT(cells_fitted_s, cid_k - 1);
        int gid_k = Rf_asInteger(cdid_ia_get_elt(cf_k, "group_id"));
        SEXP gf_k = VECTOR_ELT(group_fits_s, gid_k - 1);
        double a_k = a_vec[ci_k - 1];
        double alpha_c = Rf_asReal(cdid_ia_get_elt(cf_k, "alpha"));
        int pg = Rf_asInteger(cdid_ia_get_elt(gf_k, "p"));
        const int *kept_imp_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf_k, "kept_cols_F")) : NULL;
        const double *beta_imp = (pg > 0) ? REAL(cdid_ia_get_elt(gf_k, "beta")) : NULL;

        int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci_k - 1] : (int) REAL(c_start_s)[ci_k - 1];
        int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci_k - 1]   : (int) REAL(c_end_s)[ci_k - 1];

        double s_wm_k = 0.0, s_wodd_j_k = 0.0, s_wodd_k = 0.0, s_wodd_m_k = 0.0;
        double s_dimp_k[64];
        for (int g = 0; g < pg; ++g) s_dimp_k[g] = 0.0;

        for (int start = first; start <= last; start += B_chunk) {
          int end = start + B_chunk - 1;
          if (end > last) end = last;
          int n_r = end - start + 1;

          double ch_wm = 0.0, ch_wm_j = 0.0, ch_wodd_j = 0.0, ch_wodd = 0.0;
          for (int i = 0; i < n_r; ++i) {
            int r0 = cohort_order[start - 1 + i] - 1;
            double wr = w_all[r0];
            double ws = wr * a_k;
            double eta = 0.0;
            for (int m = 0; m < p_trt; ++m) {
              double xv = X_F[r0 + (size_t) (trt_cols_F[m] - 1) * (size_t) N_rows];
              xtrt_buf[i + (size_t) m * (size_t) n_r] = xv;
              eta += xv * gamma_t[m];
            }
            double pr = 0.0, odds = 0.0;
            cdid_ia_clip_logit(eta, &pr, &odds);
            double act = (pr < (1.0 - 1e-6)) ? 1.0 : 0.0;
            double wm = ws * odds;
            double wdot = wm * act;
            double wodd = wr * odds;
            wm_buf[i] = wm;
            wdot_buf[i] = wdot;
            wodd_buf[i] = wodd;

            double jh_sl = 0.0;
            for (int g = 0; g < pg; ++g) {
              double di = L_lk * X_F[r0 + (size_t) (kept_imp_F[g] - 1) * (size_t) N_rows];
              dimp_buf[i + (size_t) g * (size_t) n_r] = di;
              jh_sl += di * beta_imp[g];
            }
            double jh = alpha_c + jh_sl;
            jhat_buf[i] = jh;

            ch_wm += wm;
            ch_wm_j += wm * jh;
            ch_wodd_j += wodd * jh;
            ch_wodd += wodd;
          }
          sum_wm += ch_wm;
          sum_wm_j += ch_wm_j;
          s_wm_k += ch_wm;
          s_wodd_j_k += ch_wodd_j;
          s_wodd_k += ch_wodd;

          for (int g = 0; g < pg; ++g) {
            const double *di_col = dimp_buf + (size_t) g * (size_t) n_r;
            double ch_di = 0.0;
            for (int i = 0; i < n_r; ++i) ch_di += wm_buf[i] * di_col[i];
            s_dimp_k[g] += ch_di;
          }
          for (int m = 0; m < p_trt; ++m) {
            const double *xt_col = xtrt_buf + (size_t) m * (size_t) n_r;
            double ch_wdj = 0.0, ch_wd = 0.0;
            for (int i = 0; i < n_r; ++i) {
              ch_wdj += (wdot_buf[i] * jhat_buf[i]) * xt_col[i];
              ch_wd  += wdot_buf[i] * xt_col[i];
            }
            sum_wdot_j_X[m] += ch_wdj;
            sum_wdot_X[m]   += ch_wd;
          }

          if (has_or) {
            double ch_wm_m = 0.0, ch_wodd_m = 0.0;
            for (int i = 0; i < n_r; ++i) {
              int r0 = cohort_order[start - 1 + i] - 1;
              double mh_sl = 0.0;
              for (int m = 0; m < p_or; ++m) {
                double dor = L_lk * X_F[r0 + (size_t) (kept_or_F[m] - 1) * (size_t) N_rows];
                dor_buf[i + (size_t) m * (size_t) n_r] = dor;
                mh_sl += dor * beta_or[m];
              }
              double mh = alpha_b + mh_sl;
              mhat_buf[i] = mh;
              ch_wm_m += wm_buf[i] * mh;
              ch_wodd_m += wodd_buf[i] * mh;
            }
            sum_wm_m += ch_wm_m;
            s_wodd_m_k += ch_wodd_m;
            for (int m = 0; m < p_or; ++m) {
              const double *do_col = dor_buf + (size_t) m * (size_t) n_r;
              double ch_do = 0.0;
              for (int i = 0; i < n_r; ++i) ch_do += wm_buf[i] * do_col[i];
              sum_wm_dor[m] += ch_do;
            }
            for (int m = 0; m < p_trt; ++m) {
              const double *xt_col = xtrt_buf + (size_t) m * (size_t) n_r;
              double ch_wdm = 0.0;
              for (int i = 0; i < n_r; ++i) ch_wdm += (wdot_buf[i] * mhat_buf[i]) * xt_col[i];
              sum_wdot_m_X[m] += ch_wdm;
            }
          }
        }

        wm_k[idx_k] = s_wm_k;
        wodd_j_k[idx_k] = s_wodd_j_k;
        wodd_k[idx_k] = s_wodd_k;
        wodd_m_k[idx_k] = s_wodd_m_k;

        SEXP dk_s = PROTECT(Rf_allocVector(REALSXP, pg));
        for (int g = 0; g < pg; ++g) REAL(dk_s)[g] = s_dimp_k[g];
        SET_VECTOR_ELT(dimp_k_s, idx_k, dk_s);
        UNPROTECT(1);
      }

      SEXP dor_s = PROTECT(Rf_allocVector(REALSXP, has_or ? p_or : 0));
      for (int m = 0; m < (has_or ? p_or : 0); ++m) REAL(dor_s)[m] = sum_wm_dor[m];
      SEXP wdj_s = PROTECT(Rf_allocVector(REALSXP, p_trt));
      SEXP wd_s  = PROTECT(Rf_allocVector(REALSXP, p_trt));
      SEXP wdm_s = PROTECT(Rf_allocVector(REALSXP, p_trt));
      for (int m = 0; m < p_trt; ++m) {
        REAL(wdj_s)[m] = sum_wdot_j_X[m];
        REAL(wd_s)[m]  = sum_wdot_X[m];
        REAL(wdm_s)[m] = sum_wdot_m_X[m];
      }

      SEXP lc_s = PROTECT(Rf_allocVector(VECSXP, 12));
      SEXP lc_nms = PROTECT(Rf_allocVector(STRSXP, 12));
      SET_STRING_ELT(lc_nms, 0, Rf_mkChar("sum_wm"));
      SET_STRING_ELT(lc_nms, 1, Rf_mkChar("sum_wm_j"));
      SET_STRING_ELT(lc_nms, 2, Rf_mkChar("sum_wm_m"));
      SET_STRING_ELT(lc_nms, 3, Rf_mkChar("sum_wm_dor"));
      SET_STRING_ELT(lc_nms, 4, Rf_mkChar("sum_wm_by_k"));
      SET_STRING_ELT(lc_nms, 5, Rf_mkChar("sum_wm_dimp_by_k"));
      SET_STRING_ELT(lc_nms, 6, Rf_mkChar("sum_w_odds_j_by_k"));
      SET_STRING_ELT(lc_nms, 7, Rf_mkChar("sum_w_odds_by_k"));
      SET_STRING_ELT(lc_nms, 8, Rf_mkChar("sum_w_odds_m_by_k"));
      SET_STRING_ELT(lc_nms, 9, Rf_mkChar("sum_wdot_j_X"));
      SET_STRING_ELT(lc_nms, 10, Rf_mkChar("sum_wdot_X"));
      SET_STRING_ELT(lc_nms, 11, Rf_mkChar("sum_wdot_m_X"));
      Rf_setAttrib(lc_s, R_NamesSymbol, lc_nms);

      SET_VECTOR_ELT(lc_s, 0, Rf_ScalarReal(sum_wm));
      SET_VECTOR_ELT(lc_s, 1, Rf_ScalarReal(sum_wm_j));
      SET_VECTOR_ELT(lc_s, 2, Rf_ScalarReal(sum_wm_m));
      SET_VECTOR_ELT(lc_s, 3, dor_s);
      SET_VECTOR_ELT(lc_s, 4, wm_k_s);
      SET_VECTOR_ELT(lc_s, 5, dimp_k_s);
      SET_VECTOR_ELT(lc_s, 6, wodd_j_k_s);
      SET_VECTOR_ELT(lc_s, 7, wodd_k_s);
      SET_VECTOR_ELT(lc_s, 8, wodd_m_k_s);
      SET_VECTOR_ELT(lc_s, 9, wdj_s);
      SET_VECTOR_ELT(lc_s, 10, wd_s);
      SET_VECTOR_ELT(lc_s, 11, wdm_s);
      SET_VECTOR_ELT(link_ctrl_s, ell, lc_s);
      UNPROTECT(11);
    }
  }

  SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP out_nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(out_nms, 0, Rf_mkChar("block_forcing"));
  SET_STRING_ELT(out_nms, 1, Rf_mkChar("link_ctrl"));
  Rf_setAttrib(out_s, R_NamesSymbol, out_nms);
  SET_VECTOR_ELT(out_s, 0, block_forcing_s);
  SET_VECTOR_ELT(out_s, 1, link_ctrl_s);

  UNPROTECT(nprot + 2);
  vmaxset(vmax);
  return out_s;
}


/* Kernel 4: Augmented Backend S4, S5a, S5b, S6 in C */
SEXP C_cdid_ia_aug_stage456(
    SEXP mode_s,
    SEXP st_s,
    SEXP panel_s,
    SEXP X_F_s,
    SEXP cat_s,
    SEXP recipe_s,
    SEXP or_frame_index_s,
    SEXP trt_recipe_s,
    SEXP B_chunk_s,
    SEXP cells_s,
    SEXP extra_s) {
  if (TYPEOF(st_s) != VECSXP || TYPEOF(X_F_s) != REALSXP ||
      TYPEOF(cat_s) != VECSXP || TYPEOF(recipe_s) != VECSXP ||
      TYPEOF(cells_s) != VECSXP || TYPEOF(extra_s) != VECSXP) {
    return R_NilValue;
  }
  const void *vmax = vmaxget();

  int mode = Rf_asInteger(mode_s);
  int N = Rf_asInteger(cdid_ia_get_elt(st_s, "N"));
  int N_rows = Rf_nrows(X_F_s);
  const int *cohort_order = INTEGER(cdid_ia_get_elt(st_s, "cohort_order"));
  SEXP c_start_s = cdid_ia_get_elt(st_s, "cohort_start");
  SEXP c_end_s   = cdid_ia_get_elt(st_s, "cohort_end");
  const double *w_all = REAL(cdid_ia_get_elt(st_s, "w"));
  const double *a_vec = REAL(cdid_ia_get_elt(st_s, "a"));
  int G = (int) XLENGTH(cdid_ia_get_elt(st_s, "cohorts"));
  int B_chunk = Rf_asInteger(B_chunk_s);
  if (B_chunk < 1) B_chunk = N_rows;

  const double *X_F = REAL(X_F_s);
  int p = Rf_asInteger(cdid_ia_get_elt(recipe_s, "p"));
  if (p > 60) { vmaxset(vmax); return R_NilValue; }
  int icpt = Rf_asLogical(cdid_ia_get_elt(recipe_s, "intercept")) == TRUE;
  int icpt_idx = Rf_asInteger(cdid_ia_get_elt(recipe_s, "intercept_index"));
  SEXP slope_idx_s = cdid_ia_get_elt(recipe_s, "slope_index");
  SEXP slope_cols_nms = cdid_ia_get_elt(recipe_s, "slope_columns");

  int n_or_fi = (int) XLENGTH(or_frame_index_s);
  int *or_fi = (int *) R_alloc((size_t) (n_or_fi > 0 ? n_or_fi : 1), sizeof(int));
  if (TYPEOF(or_frame_index_s) == INTSXP) {
    for (int i = 0; i < n_or_fi; ++i) or_fi[i] = INTEGER(or_frame_index_s)[i];
  } else {
    for (int i = 0; i < n_or_fi; ++i) or_fi[i] = (int) REAL(or_frame_index_s)[i];
  }
  int intercept_col0 = (icpt && icpt_idx >= 1 && icpt_idx <= n_or_fi) ? (or_fi[icpt_idx - 1] - 1) : -1;
  int *slope_cols0 = (p > 0) ? (int *) R_alloc((size_t) p, sizeof(int)) : NULL;
  if (p > 0) {
    if (TYPEOF(slope_idx_s) == INTSXP) {
      for (int j = 0; j < p; ++j) slope_cols0[j] = or_fi[INTEGER(slope_idx_s)[j] - 1] - 1;
    } else {
      for (int j = 0; j < p; ++j) slope_cols0[j] = or_fi[(int) REAL(slope_idx_s)[j] - 1] - 1;
    }
  }

  int n_cells = (int) XLENGTH(cells_s);
  SEXP groups_s = cdid_ia_get_elt(extra_s, "groups");

  if (mode == 1) {
    if (TYPEOF(panel_s) != VECSXP) { vmaxset(vmax); return R_NilValue; }
    const int *a_index = INTEGER(cdid_ia_get_elt(st_s, "a_index"));
    SEXP Z_sel_s = cdid_ia_get_elt(st_s, "Z");
    const double *Z_sel = REAL(Z_sel_s);
    int N_z_rows = Rf_nrows(Z_sel_s);

    SEXP all_times_s = cdid_ia_get_elt(panel_s, "all_times");
    int n_dates = (int) XLENGTH(all_times_s);
    const double *all_times = REAL(all_times_s);
    SEXP post_s = cdid_ia_get_elt(panel_s, "postings");
    const int *post_a = INTEGER(cdid_ia_get_elt(post_s, "post_a_index"));
    const int *post_obs = INTEGER(cdid_ia_get_elt(post_s, "post_obs_index"));
    const int *date_start = INTEGER(cdid_ia_get_elt(post_s, "date_start"));
    const int *date_end = INTEGER(cdid_ia_get_elt(post_s, "date_end"));
    const int *date_count = INTEGER(cdid_ia_get_elt(post_s, "date_count"));
    SEXP y_obs_s = cdid_ia_get_elt(cdid_ia_get_elt(panel_s, "observations"), "y");
    const double *y_obs = REAL(y_obs_s);
    R_xlen_t n_y = XLENGTH(y_obs_s);
    R_xlen_t n_post = XLENGTH(cdid_ia_get_elt(post_s, "post_a_index"));

    SEXP sel_fit_s = cdid_ia_get_elt(extra_s, "sel_fit");
    SEXP q_of_cell_s = cdid_ia_get_elt(extra_s, "q_of_cell");
    SEXP theta_ids_s = cdid_ia_get_elt(extra_s, "theta_ids");
    SEXP used_ctrl_s = cdid_ia_get_elt(extra_s, "used_ctrl");

    SEXP sel_cells_s = cdid_ia_get_elt(sel_fit_s, "cells");
    SEXP sel_blocks_s = cdid_ia_get_elt(sel_fit_s, "blocks");
    SEXP sel_pcs_s = cdid_ia_get_elt(sel_fit_s, "pair_cohorts");
    SEXP dspec_s = cdid_ia_get_elt(sel_fit_s, "design_spec");
    int d_full = Rf_asInteger(cdid_ia_get_elt(dspec_s, "d"));
    int d_s = Rf_asInteger(cdid_ia_get_elt(dspec_s, "d_s"));
    if (d_s > 60) { vmaxset(vmax); return R_NilValue; }
    const int *sel_slope_cols = (d_s > 0) ? INTEGER(cdid_ia_get_elt(dspec_s, "slope_cols")) : NULL;
    int is_indep_period = (strcmp(CHAR(STRING_ELT(cdid_ia_get_elt(sel_fit_s, "selection_model"), 0)), "independent_period_logit") == 0);

    int *need_ctrl = (int *) R_alloc((size_t) (n_cells + 1), sizeof(int));
    memset(need_ctrl, 0, (size_t) (n_cells + 1) * sizeof(int));
    if (used_ctrl_s != R_NilValue) {
      int nu = (int) XLENGTH(used_ctrl_s);
      for (int i = 0; i < nu; ++i) {
        int c = (TYPEOF(used_ctrl_s) == INTSXP) ? INTEGER(used_ctrl_s)[i] : (int) REAL(used_ctrl_s)[i];
        if (c >= 1 && c <= n_cells) need_ctrl[c] = 1;
      }
    }

    SEXP colmax_list_s = PROTECT(Rf_allocVector(VECSXP, G));
    SEXP or_blocks_s = cdid_ia_get_elt(cat_s, "or_blocks");
    if (XLENGTH(or_blocks_s) > 0) {
      SEXP act_coh_s = cdid_ia_get_elt(cat_s, "active_cohorts");
      int n_act = (int) XLENGTH(act_coh_s);
      for (int i = 0; i < n_act; ++i) {
        int ci = (TYPEOF(act_coh_s) == INTSXP) ? INTEGER(act_coh_s)[i] : (int) REAL(act_coh_s)[i];
        int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci - 1] : (int) REAL(c_start_s)[ci - 1];
        int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci - 1]   : (int) REAL(c_end_s)[ci - 1];
        SEXP cm_s = PROTECT(Rf_allocVector(REALSXP, p));
        double *cm = REAL(cm_s);
        for (int j = 0; j < p; ++j) cm[j] = 0.0;
        for (int idx = first; idx <= last; ++idx) {
          int r0 = cohort_order[idx - 1] - 1;
          for (int j = 0; j < p; ++j) {
            double v = fabs(X_F[r0 + (size_t) slope_cols0[j] * (size_t) N_rows]);
            if (v > cm[j]) cm[j] = v;
          }
        }
        SET_VECTOR_ELT(colmax_list_s, ci - 1, cm_s);
        UNPROTECT(1);
      }
    }

    SEXP cell_aug_s = PROTECT(Rf_allocVector(VECSXP, n_cells));
    SEXP ctrl_acc_s = PROTECT(Rf_allocVector(VECSXP, n_cells));
    SEXP ev_cache_s = PROTECT(Rf_allocVector(VECSXP, n_cells));

    int ncols = p + 1 + (icpt ? 1 : 0);
    int *chunk_rows0 = (int *) R_alloc((size_t) B_chunk, sizeof(int));
    double *z_chunk  = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *ws_chunk = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *work_A   = (double *) R_alloc((size_t) B_chunk * (size_t) ncols, sizeof(double));
    double *chunk_fac = (double *) R_alloc((size_t) ncols * (size_t) ncols, sizeof(double));
    double *stack_A   = (double *) R_alloc((size_t) (2 * ncols) * (size_t) ncols, sizeof(double));
    double *new_fac   = (double *) R_alloc((size_t) ncols * (size_t) ncols, sizeof(double));
    double *wpi_obs_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    double *wpie_obs_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
    int *ob_rows0_buf = (int *) R_alloc((size_t) B_chunk, sizeof(int));

    for (int cid = 0; cid < n_cells; ++cid) {
      SEXP cf = VECTOR_ELT(cells_s, cid);
      int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
      SEXP gf = VECTOR_ELT(groups_s, gid - 1);
      int ci = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
      int np_expected = Rf_asInteger(cdid_ia_get_elt(cf, "n_pair"));
      double L_val = Rf_asReal(cdid_ia_get_elt(cf, "L"));
      double alpha_c = Rf_asReal(cdid_ia_get_elt(cf, "alpha"));
      double a_k = a_vec[ci - 1];
      int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
      if (pg > 60) { UNPROTECT(4); vmaxset(vmax); return R_NilValue; }
      const int *kept_imp_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf, "kept_cols_F")) : NULL;
      const double *beta_imp = (pg > 0) ? REAL(cdid_ia_get_elt(gf, "beta")) : NULL;

      int q_id = (TYPEOF(q_of_cell_s) == INTSXP) ? INTEGER(q_of_cell_s)[cid] : (int) REAL(q_of_cell_s)[cid];
      SEXP th_s = VECTOR_ELT(theta_ids_s, cid);
      int dq = (int) XLENGTH(th_s);
      if (dq > 128) { UNPROTECT(4); vmaxset(vmax); return R_NilValue; }

      SEXP pc_s = VECTOR_ELT(sel_pcs_s, q_id - 1);
      int p_cid = 0, s_cid = 0, t_cid = 0, ms = 0, mt = 0;
      int is_zero_col = 0, has_icpt_1 = 0, has_icpt_2 = 0;
      if (!is_indep_period) {
        p_cid = Rf_asInteger(cdid_ia_get_elt(pc_s, "pair_cell_id"));
        int m_bl = Rf_asInteger(cdid_ia_get_elt(pc_s, "block_id"));
        SEXP bl = VECTOR_ELT(sel_blocks_s, m_bl - 1);
        is_zero_col = (strcmp(CHAR(STRING_ELT(cdid_ia_get_elt(bl, "branch"), 0)), "zero_col") == 0);
        has_icpt_1 = Rf_asLogical(cdid_ia_get_elt(bl, "has_intercept")) == TRUE;
        int exp_dq = is_zero_col ? 0 : ((has_icpt_1 ? 1 : 0) + d_s);
        if (exp_dq != dq) { UNPROTECT(4); vmaxset(vmax); return R_NilValue; }
      } else {
        s_cid = Rf_asInteger(cdid_ia_get_elt(pc_s, "s_cell_id"));
        t_cid = Rf_asInteger(cdid_ia_get_elt(pc_s, "t_cell_id"));
        ms = Rf_asInteger(cdid_ia_get_elt(pc_s, "s_block_id"));
        mt = Rf_asInteger(cdid_ia_get_elt(pc_s, "t_block_id"));
        SEXP bl_s = VECTOR_ELT(sel_blocks_s, ms - 1);
        SEXP bl_t = VECTOR_ELT(sel_blocks_s, mt - 1);
        is_zero_col = (strcmp(CHAR(STRING_ELT(cdid_ia_get_elt(bl_s, "branch"), 0)), "zero_col") == 0);
        has_icpt_1 = Rf_asLogical(cdid_ia_get_elt(bl_s, "has_intercept")) == TRUE;
        has_icpt_2 = Rf_asLogical(cdid_ia_get_elt(bl_t, "has_intercept")) == TRUE;
        int exp_dq = 0;
        if (!is_zero_col) {
          if (ms != mt) exp_dq = (has_icpt_1 ? 1 : 0) + d_s + (has_icpt_2 ? 1 : 0) + d_s;
          else exp_dq = (has_icpt_1 ? 2 : 0) + d_s;
        }
        if (exp_dq != dq) { UNPROTECT(4); vmaxset(vmax); return R_NilValue; }
      }

      int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci - 1] : (int) REAL(c_start_s)[ci - 1];
      int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci - 1]   : (int) REAL(c_end_s)[ci - 1];
      int n_c   = last - first + 1;

      int ds = 0, dt = 0;
      if (!cdid_ia_resolve_pair_dates(cdid_ia_get_elt(cf, "pair"), all_times, n_dates, &ds, &dt)) {
        UNPROTECT(4); vmaxset(vmax); return R_NilValue;
      }
      int has_pair = (date_count[ds - 1] > 0 && date_count[dt - 1] > 0);
      int s_lo = date_start[ds - 1] - 1, s_hi = date_end[ds - 1] - 1;
      int t_lo = date_start[dt - 1] - 1, t_hi = date_end[dt - 1] - 1;
      int hs_hint = s_lo, ht_hint = t_lo;

      SEXP Z_cache_s = PROTECT(Rf_allocVector(REALSXP, n_c));
      SEXP ob_k_s    = PROTECT(Rf_allocVector(INTSXP, np_expected));
      SEXP pi_ob_s   = PROTECT(Rf_allocVector(REALSXP, np_expected));
      SEXP e_ob_s    = PROTECT(Rf_allocVector(REALSXP, np_expected));
      SEXP T_ob_s    = PROTECT(Rf_allocMatrix(REALSXP, np_expected, dq));
      double *Z_cache = REAL(Z_cache_s);
      int *ob_k_cache = INTEGER(ob_k_s);
      double *pi_cache = REAL(pi_ob_s);
      double *e_cache  = REAL(e_ob_s);
      double *T_cache  = REAL(T_ob_s);

      double Y_aug = 0.0, MW_one = 0.0;
      double Y_aug_z[64], MW_d[64], Rq_one[128];
      for (int j = 0; j < p; ++j) Y_aug_z[j] = 0.0;
      for (int g = 0; g < pg; ++g) MW_d[g] = 0.0;
      for (int m = 0; m < dq; ++m) Rq_one[m] = 0.0;

      double q_min = R_PosInf, q_max = R_NegInf, w_max = R_NegInf;
      int q_floored = 0, n_obs_tot = 0;

      int do_ctrl = need_ctrl[cid + 1];
      cdid_ia_sup_acc_t acc_ctrl;
      if (do_ctrl) cdid_ia_sup_acc_init(&acc_ctrl, p, icpt, L_val);

      for (int start = first; start <= last; start += B_chunk) {
        int end = start + B_chunk - 1;
        if (end > last) end = last;
        int n_r = end - start + 1;
        int n_ob_ch = 0;
        int ob_start_pos = n_obs_tot;

        double ch_Yaug = 0.0, ch_MW1 = 0.0;
        for (int i = 0; i < n_r; ++i) {
          int k_idx = (start - first) + i;
          int r0 = cohort_order[start - 1 + i] - 1;
          chunk_rows0[i] = r0;
          double ws = w_all[r0] * a_k;
          ws_chunk[i] = ws;
          double jh_sl = 0.0;
          for (int g = 0; g < pg; ++g) {
            jh_sl += (L_val * X_F[r0 + (size_t) (kept_imp_F[g] - 1) * (size_t) N_rows]) * beta_imp[g];
          }
          double jh = alpha_c + jh_sl;
          double zv = jh;
          int u = a_index[r0];
          if (has_pair && u != NA_INTEGER && u > 0) {
            int hs = cdid_ia_find_posting(post_a, s_lo, s_hi, u, &hs_hint);
            if (hs >= 0) {
              int ht = cdid_ia_find_posting(post_a, t_lo, t_hi, u, &ht_hint);
              if (ht >= 0) {
                int os = post_obs[hs], ot = post_obs[ht];
                if (os < 1 || (R_xlen_t) os > n_y || ot < 1 || (R_xlen_t) ot > n_y ||
                    n_obs_tot >= np_expected) {
                  UNPROTECT(9); vmaxset(vmax); return R_NilValue;
                }
                double dy = y_obs[ot - 1] - y_obs[os - 1];
                if (!R_FINITE(dy)) { UNPROTECT(9); vmaxset(vmax); return R_NilValue; }

                double q = 0.0, ps = 0.0, pt = 0.0;
                if (!is_indep_period) {
                  if (!cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N_z_rows, d_full, sel_slope_cols, p_cid, r0, &q)) {
                    UNPROTECT(9); vmaxset(vmax); return R_NilValue;
                  }
                } else {
                  if (!cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N_z_rows, d_full, sel_slope_cols, s_cid, r0, &ps) ||
                      !cdid_ia_eval_sel_cell_prob(sel_cells_s, sel_blocks_s, Z_sel, N_z_rows, d_full, sel_slope_cols, t_cid, r0, &pt)) {
                    UNPROTECT(9); vmaxset(vmax); return R_NilValue;
                  }
                  q = ps * pt;
                  if (!R_FINITE(q) || q <= 0.0 || q >= 1.0) {
                    UNPROTECT(9); vmaxset(vmax); return R_NilValue;
                  }
                }
                double act = (q > 1e-6) ? 1.0 : 0.0;
                double q_fl = (q > 1e-6) ? q : 1e-6;
                double pi = 1.0 / q_fl;
                double e = dy - jh;
                zv = jh + pi * e;
                double w_pi = ws * pi;

                if (q < q_min) q_min = q;
                if (q > q_max) q_max = q;
                if (q <= 1e-6) q_floored++;
                if (w_pi > w_max) w_max = w_pi;

                if (!is_zero_col && dq > 0) {
                  int col = 0;
                  if (!is_indep_period) {
                    double wp = -((1.0 - q) * act);
                    if (has_icpt_1) T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = wp;
                    for (int j = 0; j < d_s; ++j) {
                      double zs = Z_sel[r0 + (size_t) (sel_slope_cols[j] - 1) * (size_t) N_z_rows];
                      T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = wp * zs;
                    }
                  } else {
                    double w_s = -((1.0 - ps) * act);
                    double w_t = -((1.0 - pt) * act);
                    if (ms != mt) {
                      if (has_icpt_1) T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = w_s;
                      for (int j = 0; j < d_s; ++j) {
                        double zs = Z_sel[r0 + (size_t) (sel_slope_cols[j] - 1) * (size_t) N_z_rows];
                        T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = w_s * zs;
                      }
                      if (has_icpt_2) T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = w_t;
                      for (int j = 0; j < d_s; ++j) {
                        double zs = Z_sel[r0 + (size_t) (sel_slope_cols[j] - 1) * (size_t) N_z_rows];
                        T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = w_t * zs;
                      }
                    } else {
                      if (has_icpt_1) {
                        T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = w_s;
                        T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = w_t;
                      }
                      double wst = w_s + w_t;
                      for (int j = 0; j < d_s; ++j) {
                        double zs = Z_sel[r0 + (size_t) (sel_slope_cols[j] - 1) * (size_t) N_z_rows];
                        T_cache[n_obs_tot + (size_t) (col++) * (size_t) np_expected] = wst * zs;
                      }
                    }
                  }
                }

                ob_k_cache[n_obs_tot] = k_idx;
                pi_cache[n_obs_tot] = pi;
                e_cache[n_obs_tot] = e;
                n_obs_tot++;

                wpi_obs_buf[n_ob_ch] = w_pi;
                wpie_obs_buf[n_ob_ch] = w_pi * e;
                ob_rows0_buf[n_ob_ch] = r0;
                n_ob_ch++;
                ch_MW1 += w_pi;
              }
            }
          }
          Z_cache[k_idx] = zv;
          z_chunk[i] = zv;
          ch_Yaug += ws * zv;
        }
        Y_aug += ch_Yaug;
        MW_one += ch_MW1;

        for (int j = 0; j < p; ++j) {
          int col0 = slope_cols0[j];
          double ch_Yz = 0.0;
          for (int i = 0; i < n_r; ++i) {
            ch_Yz += (ws_chunk[i] * z_chunk[i]) * X_F[chunk_rows0[i] + (size_t) col0 * (size_t) N_rows];
          }
          Y_aug_z[j] += ch_Yz;
        }
        if (n_ob_ch > 0) {
          for (int g = 0; g < pg; ++g) {
            int col0 = kept_imp_F[g] - 1;
            double ch_mwd = 0.0;
            for (int ii = 0; ii < n_ob_ch; ++ii) {
              ch_mwd += wpi_obs_buf[ii] * (L_val * X_F[ob_rows0_buf[ii] + (size_t) col0 * (size_t) N_rows]);
            }
            MW_d[g] += ch_mwd;
          }
          for (int m = 0; m < dq; ++m) {
            const double *t_col = T_cache + (size_t) m * (size_t) np_expected + ob_start_pos;
            double ch_rq = 0.0;
            for (int ii = 0; ii < n_ob_ch; ++ii) {
              ch_rq += wpie_obs_buf[ii] * t_col[ii];
            }
            Rq_one[m] -= ch_rq;
          }
        }
        if (do_ctrl) {
          if (!cdid_ia_sup_acc_add_chunk(
                &acc_ctrl, n_r, chunk_rows0, z_chunk, ws_chunk, X_F, N_rows,
                intercept_col0, slope_cols0, work_A, chunk_fac, stack_A, new_fac)) {
            UNPROTECT(9); vmaxset(vmax); return R_NilValue;
          }
        }
      }

      if (n_obs_tot != np_expected) {
        UNPROTECT(9); vmaxset(vmax); return R_NilValue;
      }

      Y_aug /= (double) N;
      MW_one /= (double) N;
      double M_one = Rf_asReal(cdid_ia_get_elt(cf, "M_one"));
      double MD_one = M_one - MW_one;
      if (!R_FINITE(Y_aug) || !R_FINITE(MW_one) || !R_FINITE(MD_one)) {
        UNPROTECT(9); vmaxset(vmax); return R_NilValue;
      }

      SEXP Yaug_z_s = PROTECT(Rf_allocVector(REALSXP, p));
      for (int j = 0; j < p; ++j) {
        double v = Y_aug_z[j] / (double) N;
        if (!R_FINITE(v)) { UNPROTECT(10); vmaxset(vmax); return R_NilValue; }
        REAL(Yaug_z_s)[j] = v;
      }
      if (TYPEOF(slope_cols_nms) == STRSXP) Rf_setAttrib(Yaug_z_s, R_NamesSymbol, slope_cols_nms);

      SEXP MW_d_s = PROTECT(Rf_allocVector(REALSXP, pg));
      SEXP MD_d_s = PROTECT(Rf_allocVector(REALSXP, pg));
      const double *mdk = REAL(cdid_ia_get_elt(cf, "M_d_kept"));
      for (int g = 0; g < pg; ++g) {
        double mwd = MW_d[g] / (double) N;
        double mdd = mdk[g] - mwd;
        if (!R_FINITE(mwd) || !R_FINITE(mdd)) { UNPROTECT(12); vmaxset(vmax); return R_NilValue; }
        REAL(MW_d_s)[g] = mwd;
        REAL(MD_d_s)[g] = mdd;
      }
      SEXP kc_imp_s = cdid_ia_get_elt(gf, "kept_columns");
      if (TYPEOF(kc_imp_s) == STRSXP) {
        Rf_setAttrib(MW_d_s, R_NamesSymbol, kc_imp_s);
        Rf_setAttrib(MD_d_s, R_NamesSymbol, kc_imp_s);
      }

      SEXP Rq_one_s = PROTECT(Rf_allocVector(REALSXP, dq));
      for (int m = 0; m < dq; ++m) {
        double v = Rq_one[m] / (double) N;
        if (!R_FINITE(v)) { UNPROTECT(13); vmaxset(vmax); return R_NilValue; }
        REAL(Rq_one_s)[m] = v;
      }

      SEXP clip_s = PROTECT(Rf_allocVector(VECSXP, 7));
      SEXP clip_nms = PROTECT(Rf_allocVector(STRSXP, 7));
      SET_STRING_ELT(clip_nms, 0, Rf_mkChar("eps"));
      SET_STRING_ELT(clip_nms, 1, Rf_mkChar("q_min"));
      SET_STRING_ELT(clip_nms, 2, Rf_mkChar("q_max"));
      SET_STRING_ELT(clip_nms, 3, Rf_mkChar("q_floor_min"));
      SET_STRING_ELT(clip_nms, 4, Rf_mkChar("q_n"));
      SET_STRING_ELT(clip_nms, 5, Rf_mkChar("q_floored"));
      SET_STRING_ELT(clip_nms, 6, Rf_mkChar("w_max"));
      Rf_setAttrib(clip_s, R_NamesSymbol, clip_nms);
      SET_VECTOR_ELT(clip_s, 0, Rf_ScalarReal(1e-6));
      SET_VECTOR_ELT(clip_s, 1, Rf_ScalarReal((n_obs_tot > 0) ? q_min : NA_REAL));
      SET_VECTOR_ELT(clip_s, 2, Rf_ScalarReal((n_obs_tot > 0) ? q_max : NA_REAL));
      SET_VECTOR_ELT(clip_s, 3, Rf_ScalarReal((n_obs_tot > 0) ? ((q_min > 1e-6) ? q_min : 1e-6) : NA_REAL));
      SET_VECTOR_ELT(clip_s, 4, Rf_ScalarReal((double) n_obs_tot));
      SET_VECTOR_ELT(clip_s, 5, Rf_ScalarReal((double) q_floored));
      SET_VECTOR_ELT(clip_s, 6, Rf_ScalarReal((n_obs_tot > 0) ? w_max : NA_REAL));

      SEXP caug_s = PROTECT(Rf_allocVector(VECSXP, 8));
      SEXP caug_nms = PROTECT(Rf_allocVector(STRSXP, 8));
      SET_STRING_ELT(caug_nms, 0, Rf_mkChar("Yaug_one"));
      SET_STRING_ELT(caug_nms, 1, Rf_mkChar("Yaug_z"));
      SET_STRING_ELT(caug_nms, 2, Rf_mkChar("MW_one"));
      SET_STRING_ELT(caug_nms, 3, Rf_mkChar("MW_d"));
      SET_STRING_ELT(caug_nms, 4, Rf_mkChar("MD_one"));
      SET_STRING_ELT(caug_nms, 5, Rf_mkChar("MD_d"));
      SET_STRING_ELT(caug_nms, 6, Rf_mkChar("Rq_one"));
      SET_STRING_ELT(caug_nms, 7, Rf_mkChar("clipping"));
      Rf_setAttrib(caug_s, R_NamesSymbol, caug_nms);
      SET_VECTOR_ELT(caug_s, 0, Rf_ScalarReal(Y_aug));
      SET_VECTOR_ELT(caug_s, 1, Yaug_z_s);
      SET_VECTOR_ELT(caug_s, 2, Rf_ScalarReal(MW_one));
      SET_VECTOR_ELT(caug_s, 3, MW_d_s);
      SET_VECTOR_ELT(caug_s, 4, Rf_ScalarReal(MD_one));
      SET_VECTOR_ELT(caug_s, 5, MD_d_s);
      SET_VECTOR_ELT(caug_s, 6, Rq_one_s);
      SET_VECTOR_ELT(caug_s, 7, clip_s);
      SET_VECTOR_ELT(cell_aug_s, cid, caug_s);

      SEXP cev_s = PROTECT(Rf_allocVector(VECSXP, 5));
      SET_VECTOR_ELT(cev_s, 0, Z_cache_s);
      SET_VECTOR_ELT(cev_s, 1, ob_k_s);
      SET_VECTOR_ELT(cev_s, 2, pi_ob_s);
      SET_VECTOR_ELT(cev_s, 3, e_ob_s);
      SET_VECTOR_ELT(cev_s, 4, T_ob_s);
      SET_VECTOR_ELT(ev_cache_s, cid, cev_s);
      UNPROTECT(14);

      if (do_ctrl) {
        SET_VECTOR_ELT(ctrl_acc_s, cid, cdid_ia_sup_acc_to_sexp(&acc_ctrl, slope_cols_nms));
      }
    }

    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP out_nms = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(out_nms, 0, Rf_mkChar("cell_aug"));
    SET_STRING_ELT(out_nms, 1, Rf_mkChar("ctrl_acc"));
    SET_STRING_ELT(out_nms, 2, Rf_mkChar("cohort_colmax_F"));
    SET_STRING_ELT(out_nms, 3, Rf_mkChar("cell_ev_cache"));
    Rf_setAttrib(out_s, R_NamesSymbol, out_nms);
    SET_VECTOR_ELT(out_s, 0, cell_aug_s);
    SET_VECTOR_ELT(out_s, 1, ctrl_acc_s);
    SET_VECTOR_ELT(out_s, 2, colmax_list_s);
    SET_VECTOR_ELT(out_s, 3, ev_cache_s);
    UNPROTECT(6);
    vmaxset(vmax);
    return out_s;
  }

  /* Mode 2: S5b OR Block Forcing + S6 Link Control Scans from cell_ev_cache */
  if (mode == 2) {
    SEXP ev_cache_s    = cdid_ia_get_elt(extra_s, "cell_ev_cache");
    SEXP theta_ids_s   = cdid_ia_get_elt(extra_s, "theta_ids");
    SEXP block_specs_s = cdid_ia_get_elt(extra_s, "block_specs");
    SEXP treatments_s  = cdid_ia_get_elt(extra_s, "treatments");
    int nprot = 0;

    SEXP block_forcing_s = R_NilValue;
    if (block_specs_s != R_NilValue && TYPEOF(block_specs_s) == VECSXP) {
      SEXP or_blocks_s = cdid_ia_get_elt(cat_s, "or_blocks");
      int n_blocks = (int) XLENGTH(or_blocks_s);
      block_forcing_s = PROTECT(Rf_allocVector(VECSXP, n_blocks));
      nprot++;

      double *r_or_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *w_raw_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *w_syn_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *w_D_buf   = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *d_tilde_buf = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));
      double *d_imp_buf   = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));
      double *wpie_ob_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      int *ob_loc_buf     = (int *) R_alloc((size_t) B_chunk, sizeof(int));

      for (int bid = 0; bid < n_blocks; ++bid) {
        SEXP b_spec = VECTOR_ELT(or_blocks_s, bid);
        SEXP bs = VECTOR_ELT(block_specs_s, bid);
        double L_b = Rf_asReal(cdid_ia_get_elt(b_spec, "L"));
        SEXP cids_b_s = cdid_ia_get_elt(b_spec, "control_cell_ids");
        int n_kb = (int) XLENGTH(cids_b_s);

        SEXP kept_or_F_s = cdid_ia_get_elt(bs, "kept_or_cols_F");
        SEXP beta_or_s   = cdid_ia_get_elt(bs, "beta_or");
        double alpha_b   = Rf_asReal(cdid_ia_get_elt(bs, "alpha_b"));
        SEXP mean_d_b_s  = cdid_ia_get_elt(bs, "mean_d_b");
        int p_or = (int) XLENGTH(beta_or_s);
        if (p_or > 60) { UNPROTECT(nprot); vmaxset(vmax); return R_NilValue; }
        const int *kept_or_F = (p_or > 0) ? INTEGER(kept_or_F_s) : NULL;
        const double *beta_or = (p_or > 0) ? REAL(beta_or_s) : NULL;
        const double *mean_d_b = (p_or > 0) ? REAL(mean_d_b_s) : NULL;

        SEXP a_fa_s = PROTECT(Rf_allocVector(REALSXP, n_kb));
        SEXP a_fb_s = PROTECT(Rf_allocMatrix(REALSXP, p_or, n_kb));
        SEXP j_faa_s = PROTECT(Rf_allocVector(REALSXP, n_kb));
        SEXP j_fab_s = PROTECT(Rf_allocVector(VECSXP, n_kb));
        SEXP j_fba_s = PROTECT(Rf_allocMatrix(REALSXP, p_or, n_kb));
        SEXP j_fbb_s = PROTECT(Rf_allocVector(VECSXP, n_kb));
        SEXP rq_dt_s = PROTECT(Rf_allocVector(VECSXP, n_kb));

        double *a_fa = REAL(a_fa_s);
        double *a_fb = REAL(a_fb_s);
        double *j_faa = REAL(j_faa_s);
        double *j_fba = REAL(j_fba_s);
        memset(a_fb, 0, (size_t) p_or * (size_t) n_kb * sizeof(double));
        memset(j_fba, 0, (size_t) p_or * (size_t) n_kb * sizeof(double));

        for (int idx_k = 0; idx_k < n_kb; ++idx_k) {
          int cid = (TYPEOF(cids_b_s) == INTSXP) ? INTEGER(cids_b_s)[idx_k] : (int) REAL(cids_b_s)[idx_k];
          SEXP cf = VECTOR_ELT(cells_s, cid - 1);
          int gid = Rf_asInteger(cdid_ia_get_elt(cf, "group_id"));
          SEXP gf = VECTOR_ELT(groups_s, gid - 1);
          int ci = Rf_asInteger(cdid_ia_get_elt(cf, "cohort_index"));
          double a_k = a_vec[ci - 1];
          int pg = Rf_asInteger(cdid_ia_get_elt(gf, "p"));
          const int *kept_imp_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf, "kept_cols_F")) : NULL;
          int dq = (int) XLENGTH(VECTOR_ELT(theta_ids_s, cid - 1));

          SEXP cev = VECTOR_ELT(ev_cache_s, cid - 1);
          const double *Z_cache = REAL(VECTOR_ELT(cev, 0));
          const int *ob_k_cache = INTEGER(VECTOR_ELT(cev, 1));
          const double *pi_cache = REAL(VECTOR_ELT(cev, 2));
          const double *e_cache  = REAL(VECTOR_ELT(cev, 3));
          const double *T_cache  = REAL(VECTOR_ELT(cev, 4));
          int np = (int) XLENGTH(VECTOR_ELT(cev, 1));

          int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci - 1] : (int) REAL(c_start_s)[ci - 1];
          int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci - 1]   : (int) REAL(c_end_s)[ci - 1];

          double sum_w_r = 0.0;
          double sum_w_r_dt[64], sum_wD_dt[64];
          for (int m = 0; m < p_or; ++m) { sum_w_r_dt[m] = 0.0; sum_wD_dt[m] = 0.0; }

          SEXP fbb_k_s = PROTECT(Rf_allocMatrix(REALSXP, p_or, pg));
          double *fbb_k = REAL(fbb_k_s);
          memset(fbb_k, 0, (size_t) p_or * (size_t) pg * sizeof(double));

          SEXP rq_k_s = PROTECT(Rf_allocMatrix(REALSXP, p_or, dq));
          double *rq_k = REAL(rq_k_s);
          memset(rq_k, 0, (size_t) p_or * (size_t) dq * sizeof(double));

          int ob_cursor = 0;
          for (int start = first; start <= last; start += B_chunk) {
            int end = start + B_chunk - 1;
            if (end > last) end = last;
            int n_r = end - start + 1;
            int k_base = start - first;

            double ch_wr = 0.0;
            for (int i = 0; i < n_r; ++i) {
              int k_idx = k_base + i;
              int r0 = cohort_order[start - 1 + i] - 1;
              double wr = w_all[r0];
              double ws = wr * a_k;
              w_raw_buf[i] = wr;
              w_syn_buf[i] = ws;
              w_D_buf[i] = ws;

              double mh_sl = 0.0;
              for (int m = 0; m < p_or; ++m) {
                double dor = L_b * X_F[r0 + (size_t) (kept_or_F[m] - 1) * (size_t) N_rows];
                double dt = icpt ? (dor - mean_d_b[m]) : dor;
                d_tilde_buf[i + (size_t) m * (size_t) n_r] = dt;
                mh_sl += dor * beta_or[m];
              }
              double mh = alpha_b + mh_sl;
              double r_or = Z_cache[k_idx] - mh;
              r_or_buf[i] = r_or;
              ch_wr += wr * r_or;
            }
            sum_w_r += ch_wr;

            int ob_start_pos = ob_cursor;
            int n_ob_ch = 0;
            while (ob_cursor < np && ob_k_cache[ob_cursor] < k_base + n_r) {
              int loc = ob_k_cache[ob_cursor] - k_base;
              double pi = pi_cache[ob_cursor];
              double e  = e_cache[ob_cursor];
              w_D_buf[loc] = w_syn_buf[loc] * (1.0 - pi);
              wpie_ob_buf[n_ob_ch] = w_syn_buf[loc] * pi * e;
              ob_loc_buf[n_ob_ch] = loc;
              n_ob_ch++;
              ob_cursor++;
            }

            if (p_or > 0) {
              if (pg > 0) {
                for (int g = 0; g < pg; ++g) {
                  int col0 = kept_imp_F[g] - 1;
                  double *di_col = d_imp_buf + (size_t) g * (size_t) n_r;
                  for (int i = 0; i < n_r; ++i) {
                    int r0 = cohort_order[start - 1 + i] - 1;
                    di_col[i] = L_b * X_F[r0 + (size_t) col0 * (size_t) N_rows];
                  }
                }
              }
              for (int m = 0; m < p_or; ++m) {
                const double *dt_col = d_tilde_buf + (size_t) m * (size_t) n_r;
                double ch_wrdt = 0.0, ch_wDdt = 0.0;
                for (int i = 0; i < n_r; ++i) {
                  ch_wrdt += (w_raw_buf[i] * r_or_buf[i]) * dt_col[i];
                  ch_wDdt += w_D_buf[i] * dt_col[i];
                }
                sum_w_r_dt[m] += ch_wrdt;
                sum_wD_dt[m]  += ch_wDdt;
                for (int g = 0; g < pg; ++g) {
                  const double *di_col = d_imp_buf + (size_t) g * (size_t) n_r;
                  double ch_dd = 0.0;
                  for (int i = 0; i < n_r; ++i) {
                    ch_dd += (w_D_buf[i] * dt_col[i]) * di_col[i];
                  }
                  fbb_k[m + (size_t) g * (size_t) p_or] += ch_dd;
                }
                if (n_ob_ch > 0 && dq > 0) {
                  for (int col = 0; col < dq; ++col) {
                    const double *t_col = T_cache + (size_t) col * (size_t) np + ob_start_pos;
                    double ch_rq = 0.0;
                    for (int ii = 0; ii < n_ob_ch; ++ii) {
                      ch_rq += (dt_col[ob_loc_buf[ii]] * wpie_ob_buf[ii]) * t_col[ii];
                    }
                    rq_k[m + (size_t) col * (size_t) p_or] -= ch_rq;
                  }
                }
              }
            }
          }

          a_fa[idx_k] = sum_w_r / (double) N;
          for (int m = 0; m < p_or; ++m) {
            a_fb[m + (size_t) idx_k * (size_t) p_or] = sum_w_r_dt[m] / (double) N;
            j_fba[m + (size_t) idx_k * (size_t) p_or] = icpt ? (sum_wD_dt[m] / (double) N) : 0.0;
          }
          j_faa[idx_k] = icpt ? Rf_asReal(cdid_ia_get_elt(cf, "MD_one")) : 0.0;

          SEXP mdd_s = cdid_ia_get_elt(cf, "MD_d");
          SEXP fab_k_s = PROTECT(Rf_allocVector(REALSXP, pg));
          for (int g = 0; g < pg; ++g) REAL(fab_k_s)[g] = REAL(mdd_s)[g];
          SET_VECTOR_ELT(j_fab_s, idx_k, fab_k_s);

          for (int idx = 0; idx < p_or * pg; ++idx) fbb_k[idx] /= (double) N;
          SEXP kc_or_s = cdid_ia_get_elt(bs, "kept_or_names");
          SEXP kc_imp_s = cdid_ia_get_elt(gf, "kept_columns");
          if (p_or > 0 && pg > 0 && TYPEOF(kc_or_s) == STRSXP && TYPEOF(kc_imp_s) == STRSXP) {
            SEXP dn_s = PROTECT(Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(dn_s, 0, kc_or_s);
            SET_VECTOR_ELT(dn_s, 1, kc_imp_s);
            Rf_setAttrib(fbb_k_s, R_DimNamesSymbol, dn_s);
            UNPROTECT(1);
          }
          SET_VECTOR_ELT(j_fbb_s, idx_k, fbb_k_s);
          SET_VECTOR_ELT(rq_dt_s, idx_k, rq_k_s);
          UNPROTECT(3);
        }

        SEXP bf_s = PROTECT(Rf_allocVector(VECSXP, 7));
        SEXP bf_nms = PROTECT(Rf_allocVector(STRSXP, 7));
        SET_STRING_ELT(bf_nms, 0, Rf_mkChar("a_force_alpha"));
        SET_STRING_ELT(bf_nms, 1, Rf_mkChar("a_force_beta"));
        SET_STRING_ELT(bf_nms, 2, Rf_mkChar("j_force_alpha_alpha"));
        SET_STRING_ELT(bf_nms, 3, Rf_mkChar("j_force_alpha_beta"));
        SET_STRING_ELT(bf_nms, 4, Rf_mkChar("j_force_beta_alpha"));
        SET_STRING_ELT(bf_nms, 5, Rf_mkChar("j_force_beta_beta"));
        SET_STRING_ELT(bf_nms, 6, Rf_mkChar("Rq_dtilde_list"));
        Rf_setAttrib(bf_s, R_NamesSymbol, bf_nms);
        SET_VECTOR_ELT(bf_s, 0, a_fa_s);
        SET_VECTOR_ELT(bf_s, 1, a_fb_s);
        SET_VECTOR_ELT(bf_s, 2, j_faa_s);
        SET_VECTOR_ELT(bf_s, 3, j_fab_s);
        SET_VECTOR_ELT(bf_s, 4, j_fba_s);
        SET_VECTOR_ELT(bf_s, 5, j_fbb_s);
        SET_VECTOR_ELT(bf_s, 6, rq_dt_s);
        SET_VECTOR_ELT(block_forcing_s, bid, bf_s);
        UNPROTECT(9);
      }
    }

    /* Part B: S6 Link Control Scans */
    SEXP link_ctrl_s = R_NilValue;
    if (treatments_s != R_NilValue && TYPEOF(treatments_s) == VECSXP) {
      SEXP links_s = cdid_ia_get_elt(cat_s, "links");
      int L_links = (int) XLENGTH(links_s);
      int p_trt = Rf_asInteger(cdid_ia_get_elt(trt_recipe_s, "p"));
      if (p_trt > 60) { UNPROTECT(nprot); vmaxset(vmax); return R_NilValue; }
      const int *trt_cols_F = (p_trt > 0) ? INTEGER(cdid_ia_get_elt(trt_recipe_s, "index")) : NULL;

      link_ctrl_s = PROTECT(Rf_allocVector(VECSXP, L_links));
      nprot++;

      double *wm_buf   = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *wmD_buf  = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *wdot_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *wodd_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *z_buf    = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *mhat_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));
      double *dor_buf  = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));
      double *xtrt_buf = (double *) R_alloc((size_t) B_chunk * 64, sizeof(double));
      double *wmpie_ob_buf = (double *) R_alloc((size_t) B_chunk, sizeof(double));

      for (int ell = 0; ell < L_links; ++ell) {
        SEXP lk = VECTOR_ELT(links_s, ell);
        int trt_id = Rf_asInteger(cdid_ia_get_elt(lk, "trt_id"));
        SEXP tfit = VECTOR_ELT(treatments_s, trt_id - 1);
        const double *gamma_t = (p_trt > 0) ? REAL(cdid_ia_get_elt(tfit, "coefficients")) : NULL;
        double L_lk = Rf_asReal(cdid_ia_get_elt(lk, "L"));

        SEXP ctrl_coh_s = cdid_ia_get_elt(lk, "control_cohort_indices");
        SEXP ctrl_cid_s = cdid_ia_get_elt(lk, "control_cell_ids");
        int n_ctrl = (int) XLENGTH(ctrl_coh_s);

        int has_or = (block_specs_s != R_NilValue);
        int p_or = 0;
        double alpha_b = 0.0;
        const int *kept_or_F = NULL;
        const double *beta_or = NULL;
        if (has_or) {
          int bid = Rf_asInteger(cdid_ia_get_elt(lk, "or_block_id"));
          SEXP bs = VECTOR_ELT(block_specs_s, bid - 1);
          SEXP beta_or_s = cdid_ia_get_elt(bs, "beta_or");
          p_or = (int) XLENGTH(beta_or_s);
          kept_or_F = (p_or > 0) ? INTEGER(cdid_ia_get_elt(bs, "kept_or_cols_F")) : NULL;
          beta_or = (p_or > 0) ? REAL(beta_or_s) : NULL;
          alpha_b = icpt ? Rf_asReal(cdid_ia_get_elt(bs, "alpha_b")) : 0.0;
        }

        double sum_wm = 0.0, sum_wm_Z = 0.0, sum_wm_m = 0.0;
        double sum_wm_dor[64];
        for (int m = 0; m < p_or; ++m) sum_wm_dor[m] = 0.0;
        double sum_wdot_Z_X[64], sum_wdot_X[64], sum_wdot_m_X[64];
        for (int m = 0; m < p_trt; ++m) {
          sum_wdot_Z_X[m] = 0.0;
          sum_wdot_X[m] = 0.0;
          sum_wdot_m_X[m] = 0.0;
        }

        SEXP wm_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
        SEXP lam_C_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
        SEXP nu_C_k_s = PROTECT(Rf_allocVector(VECSXP, n_ctrl));
        SEXP wodd_Z_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
        SEXP wodd_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
        SEXP wodd_m_k_s = PROTECT(Rf_allocVector(REALSXP, n_ctrl));
        SEXP rq_o_k_s = PROTECT(Rf_allocVector(VECSXP, n_ctrl));

        double *wm_k = REAL(wm_k_s);
        double *lam_C_k = REAL(lam_C_k_s);
        double *wodd_Z_k = REAL(wodd_Z_k_s);
        double *wodd_k = REAL(wodd_k_s);
        double *wodd_m_k = REAL(wodd_m_k_s);

        for (int idx_k = 0; idx_k < n_ctrl; ++idx_k) {
          int ci_k = (TYPEOF(ctrl_coh_s) == INTSXP) ? INTEGER(ctrl_coh_s)[idx_k] : (int) REAL(ctrl_coh_s)[idx_k];
          int cid_k = (TYPEOF(ctrl_cid_s) == INTSXP) ? INTEGER(ctrl_cid_s)[idx_k] : (int) REAL(ctrl_cid_s)[idx_k];
          SEXP cf_k = VECTOR_ELT(cells_s, cid_k - 1);
          int gid_k = Rf_asInteger(cdid_ia_get_elt(cf_k, "group_id"));
          SEXP gf_k = VECTOR_ELT(groups_s, gid_k - 1);
          double a_k = a_vec[ci_k - 1];
          int pg = Rf_asInteger(cdid_ia_get_elt(gf_k, "p"));
          const int *kept_imp_F = (pg > 0) ? INTEGER(cdid_ia_get_elt(gf_k, "kept_cols_F")) : NULL;
          int dq = (int) XLENGTH(VECTOR_ELT(theta_ids_s, cid_k - 1));

          SEXP cev = VECTOR_ELT(ev_cache_s, cid_k - 1);
          const double *Z_cache = REAL(VECTOR_ELT(cev, 0));
          const int *ob_k_cache = INTEGER(VECTOR_ELT(cev, 1));
          const double *pi_cache = REAL(VECTOR_ELT(cev, 2));
          const double *e_cache  = REAL(VECTOR_ELT(cev, 3));
          const double *T_cache  = REAL(VECTOR_ELT(cev, 4));
          int np = (int) XLENGTH(VECTOR_ELT(cev, 1));

          int first = (TYPEOF(c_start_s) == INTSXP) ? INTEGER(c_start_s)[ci_k - 1] : (int) REAL(c_start_s)[ci_k - 1];
          int last  = (TYPEOF(c_end_s) == INTSXP)   ? INTEGER(c_end_s)[ci_k - 1]   : (int) REAL(c_end_s)[ci_k - 1];

          double s_wm_k = 0.0, s_wmD_k = 0.0, s_wodd_Z_k = 0.0, s_wodd_k = 0.0, s_wodd_m_k = 0.0;
          double s_dimpD_k[64], s_rq_o_k[128];
          for (int g = 0; g < pg; ++g) s_dimpD_k[g] = 0.0;
          for (int m = 0; m < dq; ++m) s_rq_o_k[m] = 0.0;

          int ob_cursor = 0;
          for (int start = first; start <= last; start += B_chunk) {
            int end = start + B_chunk - 1;
            if (end > last) end = last;
            int n_r = end - start + 1;
            int k_base = start - first;

            double ch_wm = 0.0, ch_wm_Z = 0.0, ch_wodd_Z = 0.0, ch_wodd = 0.0;
            for (int i = 0; i < n_r; ++i) {
              int k_idx = k_base + i;
              int r0 = cohort_order[start - 1 + i] - 1;
              double wr = w_all[r0];
              double ws = wr * a_k;
              double eta = 0.0;
              for (int m = 0; m < p_trt; ++m) {
                double xv = X_F[r0 + (size_t) (trt_cols_F[m] - 1) * (size_t) N_rows];
                xtrt_buf[i + (size_t) m * (size_t) n_r] = xv;
                eta += xv * gamma_t[m];
              }
              double pr = 0.0, odds = 0.0;
              cdid_ia_clip_logit(eta, &pr, &odds);
              double act = (pr < (1.0 - 1e-6)) ? 1.0 : 0.0;
              double wm = ws * odds;
              double wdot = wm * act;
              double wodd = wr * odds;
              double zv = Z_cache[k_idx];

              wm_buf[i] = wm;
              wmD_buf[i] = wm;
              wdot_buf[i] = wdot;
              wodd_buf[i] = wodd;
              z_buf[i] = zv;

              ch_wm += wm;
              ch_wm_Z += wm * zv;
              ch_wodd_Z += wodd * zv;
              ch_wodd += wodd;
            }
            sum_wm += ch_wm;
            sum_wm_Z += ch_wm_Z;
            s_wm_k += ch_wm;
            s_wodd_Z_k += ch_wodd_Z;
            s_wodd_k += ch_wodd;

            int ob_start_pos = ob_cursor;
            int n_ob_ch = 0;
            while (ob_cursor < np && ob_k_cache[ob_cursor] < k_base + n_r) {
              int loc = ob_k_cache[ob_cursor] - k_base;
              double pi = pi_cache[ob_cursor];
              double e  = e_cache[ob_cursor];
              wmD_buf[loc] = wm_buf[loc] * (1.0 - pi);
              wmpie_ob_buf[n_ob_ch] = wm_buf[loc] * pi * e;
              n_ob_ch++;
              ob_cursor++;
            }

            double ch_wmD = 0.0;
            for (int i = 0; i < n_r; ++i) ch_wmD += wmD_buf[i];
            s_wmD_k += ch_wmD;

            for (int g = 0; g < pg; ++g) {
              int col0 = kept_imp_F[g] - 1;
              double ch_diD = 0.0;
              for (int i = 0; i < n_r; ++i) {
                int r0 = cohort_order[start - 1 + i] - 1;
                ch_diD += wmD_buf[i] * (L_lk * X_F[r0 + (size_t) col0 * (size_t) N_rows]);
              }
              s_dimpD_k[g] += ch_diD;
            }
            for (int m = 0; m < p_trt; ++m) {
              const double *xt_col = xtrt_buf + (size_t) m * (size_t) n_r;
              double ch_wdZ = 0.0, ch_wd = 0.0;
              for (int i = 0; i < n_r; ++i) {
                ch_wdZ += (wdot_buf[i] * z_buf[i]) * xt_col[i];
                ch_wd  += wdot_buf[i] * xt_col[i];
              }
              sum_wdot_Z_X[m] += ch_wdZ;
              sum_wdot_X[m]   += ch_wd;
            }

            if (has_or) {
              double ch_wm_m = 0.0, ch_wodd_m = 0.0;
              for (int i = 0; i < n_r; ++i) {
                int r0 = cohort_order[start - 1 + i] - 1;
                double mh_sl = 0.0;
                for (int m = 0; m < p_or; ++m) {
                  double dor = L_lk * X_F[r0 + (size_t) (kept_or_F[m] - 1) * (size_t) N_rows];
                  dor_buf[i + (size_t) m * (size_t) n_r] = dor;
                  mh_sl += dor * beta_or[m];
                }
                double mh = alpha_b + mh_sl;
                mhat_buf[i] = mh;
                ch_wm_m += wm_buf[i] * mh;
                ch_wodd_m += wodd_buf[i] * mh;
              }
              sum_wm_m += ch_wm_m;
              s_wodd_m_k += ch_wodd_m;
              for (int m = 0; m < p_or; ++m) {
                const double *do_col = dor_buf + (size_t) m * (size_t) n_r;
                double ch_do = 0.0;
                for (int i = 0; i < n_r; ++i) ch_do += wm_buf[i] * do_col[i];
                sum_wm_dor[m] += ch_do;
              }
              for (int m = 0; m < p_trt; ++m) {
                const double *xt_col = xtrt_buf + (size_t) m * (size_t) n_r;
                double ch_wdm = 0.0;
                for (int i = 0; i < n_r; ++i) ch_wdm += (wdot_buf[i] * mhat_buf[i]) * xt_col[i];
                sum_wdot_m_X[m] += ch_wdm;
              }
            }

            if (n_ob_ch > 0 && dq > 0) {
              for (int m = 0; m < dq; ++m) {
                const double *t_col = T_cache + (size_t) m * (size_t) np + ob_start_pos;
                double ch_rq = 0.0;
                for (int ii = 0; ii < n_ob_ch; ++ii) {
                  ch_rq += wmpie_ob_buf[ii] * t_col[ii];
                }
                s_rq_o_k[m] -= ch_rq;
              }
            }
          }

          wm_k[idx_k] = s_wm_k;
          lam_C_k[idx_k] = s_wmD_k;
          wodd_Z_k[idx_k] = s_wodd_Z_k;
          wodd_k[idx_k] = s_wodd_k;
          wodd_m_k[idx_k] = s_wodd_m_k;

          SEXP nuk_s = PROTECT(Rf_allocVector(REALSXP, pg));
          for (int g = 0; g < pg; ++g) REAL(nuk_s)[g] = s_dimpD_k[g];
          SET_VECTOR_ELT(nu_C_k_s, idx_k, nuk_s);

          SEXP rqk_s = PROTECT(Rf_allocVector(REALSXP, dq));
          for (int m = 0; m < dq; ++m) REAL(rqk_s)[m] = s_rq_o_k[m];
          SET_VECTOR_ELT(rq_o_k_s, idx_k, rqk_s);
          UNPROTECT(2);
        }

        SEXP dor_s = PROTECT(Rf_allocVector(REALSXP, has_or ? p_or : 0));
        for (int m = 0; m < (has_or ? p_or : 0); ++m) REAL(dor_s)[m] = sum_wm_dor[m];
        SEXP wdZ_s = PROTECT(Rf_allocVector(REALSXP, p_trt));
        SEXP wd_s  = PROTECT(Rf_allocVector(REALSXP, p_trt));
        SEXP wdm_s = PROTECT(Rf_allocVector(REALSXP, p_trt));
        for (int m = 0; m < p_trt; ++m) {
          REAL(wdZ_s)[m] = sum_wdot_Z_X[m];
          REAL(wd_s)[m]  = sum_wdot_X[m];
          REAL(wdm_s)[m] = sum_wdot_m_X[m];
        }

        SEXP lc_s = PROTECT(Rf_allocVector(VECSXP, 14));
        SEXP lc_nms = PROTECT(Rf_allocVector(STRSXP, 14));
        SET_STRING_ELT(lc_nms, 0, Rf_mkChar("sum_wm"));
        SET_STRING_ELT(lc_nms, 1, Rf_mkChar("sum_wm_Z"));
        SET_STRING_ELT(lc_nms, 2, Rf_mkChar("sum_wm_m"));
        SET_STRING_ELT(lc_nms, 3, Rf_mkChar("sum_wm_dor"));
        SET_STRING_ELT(lc_nms, 4, Rf_mkChar("sum_wm_by_k"));
        SET_STRING_ELT(lc_nms, 5, Rf_mkChar("lam_C_k"));
        SET_STRING_ELT(lc_nms, 6, Rf_mkChar("nu_C_k"));
        SET_STRING_ELT(lc_nms, 7, Rf_mkChar("sum_w_odds_Z_by_k"));
        SET_STRING_ELT(lc_nms, 8, Rf_mkChar("sum_w_odds_by_k"));
        SET_STRING_ELT(lc_nms, 9, Rf_mkChar("sum_w_odds_m_by_k"));
        SET_STRING_ELT(lc_nms, 10, Rf_mkChar("sum_wdot_Z_X"));
        SET_STRING_ELT(lc_nms, 11, Rf_mkChar("sum_wdot_X"));
        SET_STRING_ELT(lc_nms, 12, Rf_mkChar("sum_wdot_m_X"));
        SET_STRING_ELT(lc_nms, 13, Rf_mkChar("Rq_o_k"));
        Rf_setAttrib(lc_s, R_NamesSymbol, lc_nms);

        SET_VECTOR_ELT(lc_s, 0, Rf_ScalarReal(sum_wm));
        SET_VECTOR_ELT(lc_s, 1, Rf_ScalarReal(sum_wm_Z));
        SET_VECTOR_ELT(lc_s, 2, Rf_ScalarReal(sum_wm_m));
        SET_VECTOR_ELT(lc_s, 3, dor_s);
        SET_VECTOR_ELT(lc_s, 4, wm_k_s);
        SET_VECTOR_ELT(lc_s, 5, lam_C_k_s);
        SET_VECTOR_ELT(lc_s, 6, nu_C_k_s);
        SET_VECTOR_ELT(lc_s, 7, wodd_Z_k_s);
        SET_VECTOR_ELT(lc_s, 8, wodd_k_s);
        SET_VECTOR_ELT(lc_s, 9, wodd_m_k_s);
        SET_VECTOR_ELT(lc_s, 10, wdZ_s);
        SET_VECTOR_ELT(lc_s, 11, wd_s);
        SET_VECTOR_ELT(lc_s, 12, wdm_s);
        SET_VECTOR_ELT(lc_s, 13, rq_o_k_s);
        SET_VECTOR_ELT(link_ctrl_s, ell, lc_s);
        UNPROTECT(13);
      }
    }

    SEXP out_s = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP out_nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(out_nms, 0, Rf_mkChar("block_forcing"));
    SET_STRING_ELT(out_nms, 1, Rf_mkChar("link_ctrl"));
    Rf_setAttrib(out_s, R_NamesSymbol, out_nms);
    SET_VECTOR_ELT(out_s, 0, block_forcing_s);
    SET_VECTOR_ELT(out_s, 1, link_ctrl_s);
    UNPROTECT(nprot + 2);
    vmaxset(vmax);
    return out_s;
  }

  vmaxset(vmax);
  return R_NilValue;
}

static inline void cdid_ia_dd_add(double *hi, double *lo, double x) {
  double s = *hi + x;
  double v = s - *hi;
  double e = (*hi - (s - v)) + (x - v);
  *hi = s;
  *lo += e;
}

static inline void cdid_ia_dd_prod_add(double *hi, double *lo, double a, double b) {
  double p = a * b;
  double err = fma(a, b, -p);
  cdid_ia_dd_add(hi, lo, p);
  *lo += err;
}

static int cdid_ia_qr_wls_solve_dd(
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
          cdid_ia_dd_prod_add(&s_hi, &s_lo, u, u);
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
        for (int i = step + 1; i < n; ++i) cdid_ia_dd_prod_add(&d_hi, &d_lo, col_s[i], col_j[i]);
        double tdot = tau * (d_hi + d_lo);
        col_j[step] -= tdot;
        for (int i = step + 1; i < n; ++i) col_j[i] = fma(-tdot, col_s[i], col_j[i]);
      }

      if (w_z != NULL) {
        double dz_hi = w_z[step], dz_lo = 0.0;
        for (int i = step + 1; i < n; ++i) cdid_ia_dd_prod_add(&dz_hi, &dz_lo, col_s[i], w_z[i]);
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
      for (int k = i + 1; k < p; ++k) cdid_ia_dd_prod_add(&s_hi, &s_lo, -R_out[i + k * p], beta_out[k]);
      double b_i = (s_hi + s_lo) / R_out[i + i * p];
      if (!R_FINITE(b_i)) return 0;
      beta_out[i] = b_i;
    }
  }
  return 1;
}

static int cdid_ia_chol_upper(
    const double *H,
    double *R_out,
    int p,
    double rel_tol) {
  if (p < 1) return 0;
  memset(R_out, 0, (size_t) p * (size_t) p * sizeof(double));
  double max_diag = 0.0;
  for (int j = 0; j < p; ++j) {
    double d = fabs(H[j + j * p]);
    if (d > max_diag) max_diag = d;
  }
  if (!R_FINITE(max_diag) || !(max_diag > 0.0)) return 0;
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      double s = H[i + j * p];
      for (int k = 0; k < i; ++k) {
        s -= R_out[k + i * p] * R_out[k + j * p];
      }
      double r_ii = R_out[i + i * p];
      if (!(r_ii > 0.0)) return 0;
      R_out[i + j * p] = s / r_ii;
    }
    double d = H[j + j * p];
    for (int k = 0; k < j; ++k) {
      double r_kj = R_out[k + j * p];
      d -= r_kj * r_kj;
    }
    if (!R_FINITE(d) || !(d > 0.0)) return 0;
    double r_jj = sqrt(d);
    if (!R_FINITE(r_jj) || !(r_jj > 0.0)) return 0;
    R_out[j + j * p] = r_jj;
  }
  double r_max = 0.0, r_min = DBL_MAX;
  for (int j = 0; j < p; ++j) {
    double r_jj = R_out[j + j * p];
    if (r_jj > r_max) r_max = r_jj;
    if (r_jj < r_min) r_min = r_jj;
  }
  if (!(r_max > 0.0) || !(r_min > rel_tol * r_max)) return 0;
  return 1;
}

/* ---------------------------------------------------------------------------
   9. C_cdid_ia_fit_sel_blocks
      Batch-fits all M selection blocks in a single .Call for fit_augmented_backend.
   --------------------------------------------------------------------------- */
SEXP C_cdid_ia_fit_sel_blocks(
    SEXP st_s,
    SEXP panel_s,
    SEXP cells_s,
    SEXP block_specs_s,
    SEXP dspec_s,
    SEXP params_s) {
  if (TYPEOF(st_s) != VECSXP || TYPEOF(panel_s) != VECSXP ||
      TYPEOF(cells_s) != VECSXP || TYPEOF(block_specs_s) != VECSXP ||
      TYPEOF(dspec_s) != VECSXP || TYPEOF(params_s) != REALSXP || XLENGTH(params_s) != 4) {
    return R_NilValue;
  }
  int is_date = (int) REAL(params_s)[0];
  int is_sep_pair = (int) REAL(params_s)[1];
  double eps = REAL(params_s)[2];
  int maxit = (int) REAL(params_s)[3];
  if (!R_FINITE(eps) || eps <= 0.0 || maxit < 1) return R_NilValue;

  SEXP a_idx_s = cdid_ia_get_elt(st_s, "a_index");
  SEXP co_s = cdid_ia_get_elt(st_s, "cohort_order");
  SEXP cs_s = cdid_ia_get_elt(st_s, "cohort_start");
  SEXP ce_s = cdid_ia_get_elt(st_s, "cohort_end");
  SEXP U_s = cdid_ia_get_elt(st_s, "U");
  SEXP w_s = cdid_ia_get_elt(st_s, "w");
  SEXP a_s = cdid_ia_get_elt(st_s, "a");
  SEXP mu_s = cdid_ia_get_elt(st_s, "mass_U");
  SEXP Z_s = cdid_ia_get_elt(st_s, "Z");
  SEXP N_s = cdid_ia_get_elt(st_s, "N");
  if (TYPEOF(a_idx_s) != INTSXP || TYPEOF(co_s) != INTSXP ||
      (TYPEOF(cs_s) != INTSXP && TYPEOF(cs_s) != REALSXP) ||
      (TYPEOF(ce_s) != INTSXP && TYPEOF(ce_s) != REALSXP) ||
      TYPEOF(U_s) != LGLSXP || TYPEOF(w_s) != REALSXP ||
      TYPEOF(a_s) != REALSXP || TYPEOF(mu_s) != REALSXP || TYPEOF(Z_s) != REALSXP) {
    return R_NilValue;
  }
  int N = (int) XLENGTH(a_idx_s);
  int G = (int) XLENGTH(cs_s);
  double N_dbl = Rf_asReal(N_s);
  if (N < 1 || G < 1 || !R_FINITE(N_dbl) || N_dbl <= 0.0) return R_NilValue;
  SEXP zdim = Rf_getAttrib(Z_s, R_DimSymbol);
  if (TYPEOF(zdim) != INTSXP || XLENGTH(zdim) != 2 || INTEGER(zdim)[0] != N) return R_NilValue;
  int d_Z = INTEGER(zdim)[1];
  const double *Z_static = REAL(Z_s);
  const int *a_index = INTEGER(a_idx_s);
  const int *cohort_order = INTEGER(co_s);
  const int *U_logical = LOGICAL(U_s);
  const double *w_static = REAL(w_s);
  const double *a_vec = REAL(a_s);
  const double *mass_U = REAL(mu_s);

  SEXP post_s = cdid_ia_get_elt(panel_s, "postings");
  SEXP obs_s = cdid_ia_get_elt(panel_s, "observations");
  if (TYPEOF(post_s) != VECSXP || TYPEOF(obs_s) != VECSXP) return R_NilValue;
  const int *date_start = INTEGER(cdid_ia_get_elt(post_s, "date_start"));
  const int *date_end   = INTEGER(cdid_ia_get_elt(post_s, "date_end"));
  const int *date_count = INTEGER(cdid_ia_get_elt(post_s, "date_count"));
  const int *post_a     = INTEGER(cdid_ia_get_elt(post_s, "post_a_index"));
  int n_dates = (int) XLENGTH(cdid_ia_get_elt(post_s, "date_count"));

  int n_sing = 0, n_panel_units = 0;
  SEXP id_s_st = R_NilValue, obs_start_s_u = R_NilValue;
  const double *ukey_p = NULL;
  const int *u_nobs_p = NULL, *o_cal_p = NULL;
  if (is_date) {
    SEXP ns_s = cdid_ia_get_elt(panel_s, "n_singleton_units");
    n_sing = (ns_s != R_NilValue) ? Rf_asInteger(ns_s) : 0;
    if (n_sing > 0) {
      SEXP id_type_s = cdid_ia_get_elt(panel_s, "id_type");
      SEXP ukey_s = cdid_ia_get_elt(panel_s, "unit_id_key");
      SEXP units_s = cdid_ia_get_elt(panel_s, "units");
      SEXP cal_s = cdid_ia_get_elt(obs_s, "calendar_index");
      id_s_st = cdid_ia_get_elt(st_s, "id");
      if (TYPEOF(id_type_s) != STRSXP || XLENGTH(id_type_s) != 1 ||
          strcmp(CHAR(STRING_ELT(id_type_s, 0)), "numeric") != 0 ||
          TYPEOF(ukey_s) != REALSXP || TYPEOF(units_s) != VECSXP || TYPEOF(cal_s) != INTSXP ||
          (TYPEOF(id_s_st) != REALSXP && TYPEOF(id_s_st) != INTSXP) || Rf_isObject(id_s_st)) {
        return R_NilValue;
      }
      SEXP nobs_s = cdid_ia_get_elt(units_s, "n_obs");
      obs_start_s_u = cdid_ia_get_elt(units_s, "obs_start");
      if (TYPEOF(nobs_s) != INTSXP ||
          (TYPEOF(obs_start_s_u) != INTSXP && TYPEOF(obs_start_s_u) != REALSXP)) {
        return R_NilValue;
      }
      ukey_p = REAL(ukey_s);
      n_panel_units = (int) XLENGTH(ukey_s);
      u_nobs_p = INTEGER(nobs_s);
      o_cal_p = INTEGER(cal_s);
    }
  }

  int has_int = (Rf_asLogical(cdid_ia_get_elt(dspec_s, "has_intercept")) == 1);
  SEXP sc_s = cdid_ia_get_elt(dspec_s, "slope_cols");
  if (TYPEOF(sc_s) != INTSXP) return R_NilValue;
  int d_s = (int) XLENGTH(sc_s);
  const int *sc = INTEGER(sc_s);
  if ((!has_int && d_s == 0) || d_Z > 128 || d_s > 128) return R_NilValue;
  for (int j = 0; j < d_s; ++j) {
    if (sc[j] < 1 || sc[j] > d_Z) return R_NilValue;
  }

  int M = (int) XLENGTH(block_specs_s);
  int n_cells_total = (int) XLENGTH(cells_s);
  if (M < 1 || n_cells_total < 1) return R_NilValue;

  SEXP blocks_out_s = PROTECT(Rf_allocVector(VECSXP, M));
  const void *vmax_outer = vmaxget();

  double rank_tol = eps / 1000.0;
  if (rank_tol > 1e-7) rank_tol = 1e-7;
  if (rank_tol < 1e-12) rank_tol = 1e-12;
  double safe_tol = (rank_tol > 1e-6) ? rank_tol : 1e-6;
  const double eps_prob = 10.0 * DBL_EPSILON;
  const double inv_N = 1.0 / N_dbl;

  for (int m_idx = 0; m_idx < M; ++m_idx) {
    const void *vmax = vmaxget();
    SEXP bs_s = VECTOR_ELT(block_specs_s, m_idx);
    SEXP cids_s = cdid_ia_get_elt(bs_s, "cell_ids");
    if (TYPEOF(cids_s) != INTSXP) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
    int B_m = (int) XLENGTH(cids_s);
    const int *cids = INTEGER(cids_s);
    if (B_m < 1) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }

    if (is_sep_pair && B_m == 1 && d_Z > 0) {
      int cid = cids[0];
      if (cid < 1 || cid > n_cells_total) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
      SEXP cell = VECTOR_ELT(cells_s, cid - 1);
      int k = Rf_asInteger(cdid_ia_get_elt(cell, "cohort"));
      if (k < 1 || k > G) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
      SEXP pair = cdid_ia_get_elt(cell, "pair");
      int ds1 = Rf_asInteger(cdid_ia_get_elt(pair, "s_calendar_index"));
      int dt1 = Rf_asInteger(cdid_ia_get_elt(pair, "t_calendar_index"));
      if (ds1 == NA_INTEGER || dt1 == NA_INTEGER || ds1 < 1 || ds1 > n_dates || dt1 < 1 || dt1 > n_dates) {
        vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
      }
      double a_k = a_vec[k - 1];
      double mass_U_k = mass_U[k - 1];
      if (!R_FINITE(a_k) || a_k <= 0.0 || !R_FINITE(mass_U_k) || mass_U_k <= 0.0) {
        vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
      }

      int c0 = ((TYPEOF(cs_s) == INTSXP) ? INTEGER(cs_s)[k - 1] : (int) REAL(cs_s)[k - 1]) - 1;
      int c1 = ((TYPEOF(ce_s) == INTSXP) ? INTEGER(ce_s)[k - 1] : (int) REAL(ce_s)[k - 1]) - 1;
      int n_fit = 0;
      for (int idx = c0; idx <= c1; ++idx) {
        if (U_logical[cohort_order[idx] - 1] == 1) ++n_fit;
      }
      if (n_fit < d_Z) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }

      double *X = (double *) R_alloc((size_t) n_fit * (size_t) d_Z, sizeof(double));
      double *y = (double *) R_alloc((size_t) n_fit, sizeof(double));
      double *w = (double *) R_alloc((size_t) n_fit, sizeof(double));

      int has_pair = (date_count[ds1 - 1] > 0 && date_count[dt1 - 1] > 0);
      int s_lo = date_start[ds1 - 1] - 1, s_hi = date_end[ds1 - 1] - 1;
      int t_lo = date_start[dt1 - 1] - 1, t_hi = date_end[dt1 - 1] - 1;
      int hs_hint = s_lo, ht_hint = t_lo;
      int pos = 0;
      double wy = 0.0, sum_w = 0.0;
      for (int idx = c0; idx <= c1; ++idx) {
        int row0 = cohort_order[idx] - 1;
        if (U_logical[row0] != 1) continue;
        double S_val = 0.0;
        int ai = a_index[row0];
        if (has_pair && ai > 0) {
          if (cdid_ia_find_posting(post_a, s_lo, s_hi, ai, &hs_hint) >= 0 &&
              cdid_ia_find_posting(post_a, t_lo, t_hi, ai, &ht_hint) >= 0) {
            S_val = 1.0;
          }
        }
        double wi = w_static[row0];
        if (!R_FINITE(wi) || wi <= 0.0) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
        for (int j = 0; j < d_Z; ++j) {
          double v = Z_static[row0 + (size_t) j * (size_t) N];
          if (!R_FINITE(v)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          X[pos + (size_t) j * (size_t) n_fit] = v;
        }
        y[pos] = S_val;
        w[pos] = wi;
        wy += wi * S_val;
        sum_w += wi;
        ++pos;
      }
      if (!R_FINITE(sum_w) || sum_w <= 0.0) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }

      for (int j = 0; j < d_Z; ++j) {
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
        if ((has_pos && !has_neg) || (has_neg && !has_pos)) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
      }

      double *eta = (double *) R_alloc((size_t) n_fit, sizeof(double));
      double *mu = (double *) R_alloc((size_t) n_fit, sizeof(double));
      double *wt = (double *) R_alloc((size_t) n_fit, sizeof(double));
      double *wz = (double *) R_alloc((size_t) n_fit, sizeof(double));
      double *W_X = (double *) R_alloc((size_t) n_fit * (size_t) d_Z, sizeof(double));
      double *H = (double *) R_alloc((size_t) d_Z * (size_t) d_Z, sizeof(double));
      double *R_chol = (double *) R_alloc((size_t) d_Z * (size_t) d_Z, sizeof(double));
      double *beta = (double *) R_alloc((size_t) d_Z, sizeof(double));

      for (int i = 0; i < n_fit; ++i) { eta[i] = 0.0; mu[i] = 0.5; }
      memset(beta, 0, (size_t) d_Z * sizeof(double));
      double devold = 2.0 * M_LN2 * sum_w;
      int conv = 0, iter_used = 0;

      for (int iter = 1; iter <= maxit; ++iter) {
        iter_used = iter;
        for (int i = 0; i < n_fit; ++i) {
          double mui = mu[i];
          double varmu = mui * (1.0 - mui);
          double opexp = exp(eta[i]);
          double mu_eta = opexp / ((1.0 + opexp) * (1.0 + opexp));
          if (!(varmu > 0.0) || !(mu_eta > 0.0) || !R_FINITE(mu_eta)) {
            vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
          }
          double w_i = sqrt((w[i] * mu_eta * mu_eta) / varmu);
          double z_i = eta[i] + (y[i] - mui) / mu_eta;
          wt[i] = w_i;
          wz[i] = z_i * w_i;
        }
        for (int j = 0; j < d_Z; ++j) {
          const double *xj = X + (size_t) j * (size_t) n_fit;
          double *wxj = W_X + (size_t) j * (size_t) n_fit;
          for (int i = 0; i < n_fit; ++i) wxj[i] = xj[i] * wt[i];
        }
        if (!cdid_ia_qr_wls_solve_dd(W_X, wz, R_chol, beta, n_fit, d_Z, safe_tol)) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
        memset(eta, 0, (size_t) n_fit * sizeof(double));
        for (int j = 0; j < d_Z; ++j) {
          double bj = beta[j];
          const double *xj = X + (size_t) j * (size_t) n_fit;
          for (int i = 0; i < n_fit; ++i) eta[i] += xj[i] * bj;
        }
        double dev = 0.0;
        for (int i = 0; i < n_fit; ++i) {
          double ei = eta[i];
          if (!R_FINITE(ei) || fabs(ei) >= 30.0) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          double mui = 1.0 / (1.0 + exp(-ei));
          if (mui <= eps_prob || mui >= 1.0 - eps_prob) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          mu[i] = mui;
          double loss = (y[i] == 1.0) ? -log(mui) : -log1p(-mui);
          dev += 2.0 * w[i] * loss;
        }
        if (!R_FINITE(dev)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
        if (fabs(dev - devold) / (0.1 + fabs(dev)) < eps) { conv = 1; break; }
        devold = dev;
      }
      if (!conv) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }

      double q_sum = 0.0;
      for (int i = 0; i < n_fit; ++i) {
        double od = exp(eta[i]);
        double pr = od / (1.0 + od);
        if (pr <= eps_prob || pr >= 1.0 - eps_prob || !(pr * (1.0 - pr) > 0.0)) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
        double qi = w[i] * pr * (1.0 - pr);
        wt[i] = qi;
        q_sum += qi;
      }
      for (int j = 0; j < d_Z; ++j) {
        const double *xj = X + (size_t) j * (size_t) n_fit;
        for (int k2 = 0; k2 <= j; ++k2) {
          const double *xk = X + (size_t) k2 * (size_t) n_fit;
          double s = 0.0;
          for (int i = 0; i < n_fit; ++i) s += wt[i] * xk[i] * xj[i];
          double hs_val = s * inv_N;
          if (!R_FINITE(hs_val)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          H[k2 + j * d_Z] = hs_val;
          H[j + k2 * d_Z] = hs_val;
        }
      }
      if (!cdid_ia_chol_upper(H, R_chol, d_Z, 1e-8)) {
        vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
      }

      double pi_val = wy / (mass_U_k * N_dbl);
      double r_val = pi_val / a_k;
      double m_b = q_sum * inv_N;

      SEXP z_bar_s = PROTECT(Rf_allocMatrix(REALSXP, 1, (has_int && d_s > 0) ? d_s : 0));
      if (has_int && d_s > 0) {
        for (int j = 0; j < d_s; ++j) {
          const double *xj = X + (size_t) (sc[j] - 1) * (size_t) n_fit;
          double zs = 0.0;
          for (int i = 0; i < n_fit; ++i) zs += wt[i] * xj[i];
          REAL(z_bar_s)[j] = zs / q_sum;
        }
      }
      SEXP b_out = PROTECT(Rf_allocVector(VECSXP, 7));
      SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, d_Z));
      memcpy(REAL(beta_s), beta, (size_t) d_Z * sizeof(double));
      SEXP R_s = PROTECT(Rf_allocMatrix(REALSXP, d_Z, d_Z));
      memcpy(REAL(R_s), R_chol, (size_t) d_Z * (size_t) d_Z * sizeof(double));
      SET_VECTOR_ELT(b_out, 0, beta_s);
      SET_VECTOR_ELT(b_out, 1, R_s);
      SET_VECTOR_ELT(b_out, 2, Rf_ScalarInteger(iter_used));
      SET_VECTOR_ELT(b_out, 3, Rf_ScalarReal(pi_val));
      SET_VECTOR_ELT(b_out, 4, Rf_ScalarReal(r_val));
      SET_VECTOR_ELT(b_out, 5, Rf_ScalarReal(m_b));
      SET_VECTOR_ELT(b_out, 6, z_bar_s);
      SET_VECTOR_ELT(blocks_out_s, m_idx, b_out);
      UNPROTECT(4);
      vmaxset(vmax);
      continue;
    }

    /* Mode 1: general pooled / date selection block */
    int p_m = (has_int ? B_m : 0) + d_s;
    int *cell_start = (int *) R_alloc((size_t) (B_m + 1), sizeof(int));
    cell_start[0] = 0;
    for (int b = 0; b < B_m; ++b) {
      int cid = cids[b];
      if (cid < 1 || cid > n_cells_total) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
      SEXP cell = VECTOR_ELT(cells_s, cid - 1);
      int k = Rf_asInteger(cdid_ia_get_elt(cell, "cohort"));
      if (k < 1 || k > G) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
      int c0 = ((TYPEOF(cs_s) == INTSXP) ? INTEGER(cs_s)[k - 1] : (int) REAL(cs_s)[k - 1]) - 1;
      int c1 = ((TYPEOF(ce_s) == INTSXP) ? INTEGER(ce_s)[k - 1] : (int) REAL(ce_s)[k - 1]) - 1;
      int nb = 0;
      for (int idx = c0; idx <= c1; ++idx) {
        if (U_logical[cohort_order[idx] - 1] == 1) ++nb;
      }
      if (nb < 1) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
      cell_start[b + 1] = cell_start[b] + nb;
    }
    int n_fit = cell_start[B_m];
    if (n_fit < p_m) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }

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
      SEXP cell = VECTOR_ELT(cells_s, cids[b] - 1);
      int k = Rf_asInteger(cdid_ia_get_elt(cell, "cohort"));
      double a_k = a_vec[k - 1];
      double mass_U_k = mass_U[k - 1];
      if (!R_FINITE(a_k) || a_k <= 0.0 || !R_FINITE(mass_U_k) || mass_U_k <= 0.0) {
        vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
      }

      int ds1 = 0, dt1 = 0;
      if (is_date) {
        SEXP dspec_c = cdid_ia_get_elt(cell, "date_spec");
        dt1 = Rf_asInteger(cdid_ia_get_elt(dspec_c != R_NilValue ? dspec_c : cell, "calendar_index"));
        if (dt1 == NA_INTEGER || dt1 < 1 || dt1 > n_dates) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
      } else {
        SEXP pair = cdid_ia_get_elt(cell, "pair");
        ds1 = Rf_asInteger(cdid_ia_get_elt(pair, "s_calendar_index"));
        dt1 = Rf_asInteger(cdid_ia_get_elt(pair, "t_calendar_index"));
        if (ds1 == NA_INTEGER || dt1 == NA_INTEGER || ds1 < 1 || ds1 > n_dates || dt1 < 1 || dt1 > n_dates) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
      }

      int c0 = ((TYPEOF(cs_s) == INTSXP) ? INTEGER(cs_s)[k - 1] : (int) REAL(cs_s)[k - 1]) - 1;
      int c1 = ((TYPEOF(ce_s) == INTSXP) ? INTEGER(ce_s)[k - 1] : (int) REAL(ce_s)[k - 1]) - 1;
      int pos = cell_start[b];
      int has_dt = (date_count[dt1 - 1] > 0);
      int t_lo = date_start[dt1 - 1] - 1, t_hi = date_end[dt1 - 1] - 1;
      int has_ds = (!is_date && date_count[ds1 - 1] > 0);
      int s_lo = !is_date ? (date_start[ds1 - 1] - 1) : 0;
      int s_hi = !is_date ? (date_end[ds1 - 1] - 1) : -1;
      int hs_hint = s_lo, ht_hint = t_lo, hu_hint = 0;
      double wy = 0.0;
      int all_one = 1, all_zero = 1;

      for (int idx = c0; idx <= c1; ++idx) {
        int row0 = cohort_order[idx] - 1;
        if (U_logical[row0] != 1) continue;
        double yi = 0.0;
        int ai = a_index[row0];
        if (is_date) {
          if (ai > 0) {
            if (has_dt && cdid_ia_find_posting(post_a, t_lo, t_hi, ai, &ht_hint) >= 0) yi = 1.0;
          } else if (n_sing > 0) {
            double id_d = (TYPEOF(id_s_st) == INTSXP) ? (double) INTEGER(id_s_st)[row0] : REAL(id_s_st)[row0];
            if (id_d == 0.0) id_d = 0.0;
            int u = -1;
            if (hu_hint >= 0 && hu_hint + 1 < n_panel_units && ukey_p[hu_hint + 1] == id_d) {
              u = ++hu_hint;
            } else {
              int lo = 0, hi = n_panel_units - 1;
              while (lo <= hi) {
                int mid = lo + ((hi - lo) >> 1);
                double v = ukey_p[mid];
                if (v == id_d) { u = hu_hint = mid; break; }
                if (v < id_d) lo = mid + 1;
                else hi = mid - 1;
              }
            }
            if (u >= 0 && u_nobs_p[u] == 1) {
              int opos = (TYPEOF(obs_start_s_u) == INTSXP)
                ? (INTEGER(obs_start_s_u)[u] - 1)
                : ((int) REAL(obs_start_s_u)[u] - 1);
              if (opos >= 0 && o_cal_p[opos] == dt1) yi = 1.0;
            }
          }
        } else {
          if (has_ds && has_dt && ai > 0) {
            if (cdid_ia_find_posting(post_a, s_lo, s_hi, ai, &hs_hint) >= 0 &&
                cdid_ia_find_posting(post_a, t_lo, t_hi, ai, &ht_hint) >= 0) {
              yi = 1.0;
            }
          }
        }
        if (yi == 1.0) all_zero = 0;
        else all_one = 0;

        double wi = w_static[row0];
        if (!R_FINITE(wi) || wi <= 0.0) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
        for (int j = 0; j < d_s; ++j) {
          double v = Z_static[row0 + (size_t) (sc[j] - 1) * (size_t) N];
          if (!R_FINITE(v)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          D_mat[pos + (size_t) j * (size_t) n_fit] = v;
          double sgn = (yi == 1.0) ? v : -v;
          if (sgn > 0.0) pos_sgn[j] = 1;
          else if (sgn < 0.0) neg_sgn[j] = 1;
        }
        y_all[pos] = yi;
        w_all[pos] = wi;
        wy += wi * yi;
        sum_w_total += wi;
        ++pos;
      }
      if (has_int && (all_one || all_zero)) {
        vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
      }
      pi_vec[b] = wy / (mass_U_k * N_dbl);
      r_vec[b] = pi_vec[b] / a_k;
    }

    for (int j = 0; j < d_s; ++j) {
      if ((pos_sgn[j] && !neg_sgn[j]) || (neg_sgn[j] && !pos_sgn[j])) {
        vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
      }
    }

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
          if (!R_FINITE(ei) || fabs(ei) > 30.0) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          double odds = exp(ei);
          double mu = (ei >= 0.0) ? (1.0 / (1.0 + exp(-ei))) : (odds / (1.0 + odds));
          if (mu <= eps_prob || mu >= 1.0 - eps_prob) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          double mu_eta = mu * (1.0 - mu);
          double qi = w_all[i] * mu_eta;
          double zi = ei + (y_all[i] - mu) / mu_eta;
          if (!R_FINITE(qi) || !(qi > 0.0) || !R_FINITE(zi)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          q_vec[i] = qi;
          z_work[i] = zi;
          cdid_ia_dd_add(&qb_hi, &qb_lo, qi);
          if (has_int) cdid_ia_dd_prod_add(&qz_hi, &qz_lo, qi, zi);
        }
        double qb_sum = qb_hi + qb_lo;
        if (!(qb_sum > 0.0) || !R_FINITE(qb_sum)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
        Q_cell[b] = qb_sum;
        if (has_int) {
          ybar_cell[b] = (qz_hi + qz_lo) / qb_sum;
          if (d_s > 0) {
            for (int j = 0; j < d_s; ++j) {
              const double *dj = D_mat + (size_t) j * (size_t) n_fit;
              double s_hi = 0.0, s_lo = 0.0;
              for (int i = s0; i < s1; ++i) cdid_ia_dd_prod_add(&s_hi, &s_lo, q_vec[i], dj[i]);
              zbar_cell[b + (size_t) j * (size_t) B_m] = (s_hi + s_lo) / qb_sum;
            }
          }
        }
      }

      if (d_s == 0) {
        for (int b = 0; b < B_m; ++b) {
          if (!R_FINITE(ybar_cell[b])) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
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
        if (!cdid_ia_qr_wls_solve_dd(W_D, wz, R_mat, beta, n_fit, d_s, safe_tol)) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
        for (int b = 0; b < B_m; ++b) {
          double ab = ybar_cell[b];
          for (int j = 0; j < d_s; ++j) ab -= zbar_cell[b + (size_t) j * (size_t) B_m] * beta[j];
          if (!R_FINITE(ab)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
          alpha[b] = ab;
        }
      } else {
        for (int i = 0; i < n_fit; ++i) wz[i] = sqrt(q_vec[i]) * z_work[i];
        for (int j = 0; j < d_s; ++j) {
          const double *dj = D_mat + (size_t) j * (size_t) n_fit;
          double *wdj = W_D + (size_t) j * (size_t) n_fit;
          for (int i = 0; i < n_fit; ++i) wdj[i] = sqrt(q_vec[i]) * dj[i];
        }
        if (!cdid_ia_qr_wls_solve_dd(W_D, wz, R_mat, beta, n_fit, d_s, safe_tol)) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
      }

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
        if (!R_FINITE(ei) || fabs(ei) > 30.0) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
        double odds = exp(ei);
        double mu = (ei >= 0.0) ? (1.0 / (1.0 + exp(-ei))) : (odds / (1.0 + odds));
        if (mu <= eps_prob || mu >= 1.0 - eps_prob) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
        double loss = (y_all[i] == 1.0) ? -log(mu) : -log1p(-mu);
        cdid_ia_dd_add(&dev_hi, &dev_lo, 2.0 * w_all[i] * loss);
      }
      double dev = dev_hi + dev_lo;
      if (!R_FINITE(dev)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
      if (fabs(dev - dev_old) / (0.1 + fabs(dev)) < eps) { conv = 1; break; }
      dev_old = dev;
    }
    if (!conv) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }

    for (int b = 0; b < B_m; ++b) {
      int s0 = cell_start[b], s1 = cell_start[b + 1];
      double qb_hi = 0.0, qb_lo = 0.0;
      for (int i = s0; i < s1; ++i) {
        double ei = eta[i];
        double odds = exp(ei);
        double mu = odds / (1.0 + odds);
        if (!R_FINITE(ei) || !R_FINITE(odds) || !(odds > 0.0) || !R_FINITE(mu) ||
            mu <= eps_prob || mu >= 1.0 - eps_prob || !(mu * (1.0 - mu) > 0.0)) {
          vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
        }
        double qi = w_all[i] * mu * (1.0 - mu);
        if (!R_FINITE(qi) || !(qi > 0.0)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
        q_vec[i] = qi;
        cdid_ia_dd_add(&qb_hi, &qb_lo, qi);
      }
      double qb_sum = qb_hi + qb_lo;
      if (!R_FINITE(qb_sum) || !(qb_sum > 0.0)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
      Q_cell[b] = qb_sum;
      if (has_int && d_s > 0) {
        for (int j = 0; j < d_s; ++j) {
          const double *dj = D_mat + (size_t) j * (size_t) n_fit;
          double s_hi = 0.0, s_lo = 0.0;
          for (int i = s0; i < s1; ++i) cdid_ia_dd_prod_add(&s_hi, &s_lo, q_vec[i], dj[i]);
          double zb = (s_hi + s_lo) / qb_sum;
          if (!R_FINITE(zb)) { vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue; }
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
      if (!cdid_ia_qr_wls_solve_dd(W_D, NULL, R_mat, NULL, n_fit, d_s, safe_tol)) {
        vmaxset(vmax_outer); UNPROTECT(1); return R_NilValue;
      }
      double max_diag = 0.0;
      for (int j = 0; j < d_s; ++j) {
        double ad = fabs(R_mat[j + j * d_s]);
        if (ad > max_diag) max_diag = ad;
      }
      cutoff = rank_tol * max_diag;
    }

    SEXP b_out = PROTECT(Rf_allocVector(VECSXP, 11));
    SEXP alpha_s = PROTECT(Rf_allocVector(REALSXP, has_int ? B_m : 0));
    if (has_int) memcpy(REAL(alpha_s), alpha, (size_t) B_m * sizeof(double));
    SEXP beta_s = PROTECT(Rf_allocVector(REALSXP, d_s));
    if (d_s > 0) memcpy(REAL(beta_s), beta, (size_t) d_s * sizeof(double));
    SEXP mb_s = PROTECT(Rf_allocVector(REALSXP, B_m));
    for (int b = 0; b < B_m; ++b) REAL(mb_s)[b] = Q_cell[b] * inv_N;
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

    SET_VECTOR_ELT(b_out, 0, alpha_s);
    SET_VECTOR_ELT(b_out, 1, beta_s);
    SET_VECTOR_ELT(b_out, 2, mb_s);
    SET_VECTOR_ELT(b_out, 3, zb_s);
    SET_VECTOR_ELT(b_out, 4, R_s);
    SET_VECTOR_ELT(b_out, 5, pi_s);
    SET_VECTOR_ELT(b_out, 6, r_s);
    SET_VECTOR_ELT(b_out, 7, Rf_ScalarInteger(iter_used));
    SET_VECTOR_ELT(b_out, 8, Rf_ScalarReal((double) n_fit));
    SET_VECTOR_ELT(b_out, 9, Rf_ScalarReal(rank_tol));
    SET_VECTOR_ELT(b_out, 10, Rf_ScalarReal(cutoff));
    SET_VECTOR_ELT(blocks_out_s, m_idx, b_out);
    UNPROTECT(n_prot);
    vmaxset(vmax);
  }

  UNPROTECT(1);
  vmaxset(vmax_outer);
  return blocks_out_s;
}
