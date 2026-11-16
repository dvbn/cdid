#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("fp-contract=off")
#endif
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#pragma STDC FP_CONTRACT OFF
#endif

#include <R.h>
#include <Rinternals.h>
#include "cdid_r_api.h"
#include <R_ext/RS.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)) || \
    (defined(__ASSOC_MATH__) && (__ASSOC_MATH__ > 0)) || \
    (defined(__RECIPROCAL_MATH__) && (__RECIPROCAL_MATH__ > 0))
#error "Fast-math / associative-math / finite-math-only compilation is forbidden."
#endif

static inline uint64_t cpt_canon_bits(double x) {
  if (x == 0.0) return 0ULL;
  uint64_t u = 0ULL;
  memcpy(&u, &x, sizeof(uint64_t));
  return u;
}

static inline uint64_t cpt_mix64(uint64_t z) {
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static inline size_t cpt_next_pow2(size_t n) {
  size_t cap = 16;
  while (cap < n) {
    if (cap > (SIZE_MAX >> 1)) return 0;
    cap <<= 1;
  }
  return cap;
}

static inline int cpt_ctz64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctzll(x);
#else
  int c = 0;
  while ((x & 1ULL) == 0ULL) {
    x >>= 1;
    ++c;
  }
  return c;
#endif
}

/* ---------------------------------------------------------------------------
   1. Canonical IEEE-754 binary64 (-0.0 == +0.0) distinct grid deduplication
   --------------------------------------------------------------------------- */
SEXP C_cdid_cpt_distinct_grid(SEXP Z_sexp) {
  if (Rf_isObject(Z_sexp) || !Rf_isMatrix(Z_sexp) || TYPEOF(Z_sexp) != REALSXP) {
    return R_NilValue;
  }
  int N = Rf_nrows(Z_sexp);
  int p = Rf_ncols(Z_sexp);
  if (N < 1 || p < 0) return R_NilValue;

  const double *Z = REAL(Z_sexp);
  R_xlen_t total = (R_xlen_t) N * (R_xlen_t) p;
  for (R_xlen_t idx = 0; idx < total; ++idx) {
    if (!R_FINITE(Z[idx])) return R_NilValue;
  }

  if (p == 0) {
    SEXP z_out = PROTECT(Rf_allocMatrix(REALSXP, 1, 0));
    SEXP first_sexp = PROTECT(Rf_allocVector(INTSXP, 1));
    SEXP map_sexp = PROTECT(Rf_allocVector(INTSXP, N));
    SEXP mult_sexp = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(first_sexp)[0] = 1;
    int *map_p = INTEGER(map_sexp);
    for (int i = 0; i < N; ++i) map_p[i] = 1;
    INTEGER(mult_sexp)[0] = N;

    SEXP res = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_STRING_ELT(nms, 0, Rf_mkChar("z"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("first"));
    SET_STRING_ELT(nms, 2, Rf_mkChar("map"));
    SET_STRING_ELT(nms, 3, Rf_mkChar("mult"));
    SET_STRING_ELT(nms, 4, Rf_mkChar("M"));
    SET_STRING_ELT(nms, 5, Rf_mkChar("n_distinct"));
    Rf_setAttrib(res, R_NamesSymbol, nms);
    SET_VECTOR_ELT(res, 0, z_out);
    SET_VECTOR_ELT(res, 1, first_sexp);
    SET_VECTOR_ELT(res, 2, map_sexp);
    SET_VECTOR_ELT(res, 3, mult_sexp);
    SET_VECTOR_ELT(res, 4, Rf_ScalarInteger(N));
    SET_VECTOR_ELT(res, 5, Rf_ScalarInteger(1));
    UNPROTECT(6);
    return res;
  }

  const void *vmax = vmaxget();
  size_t cap = cpt_next_pow2((size_t) N * 2u + 8u);
  if (cap == 0) return R_NilValue;
  size_t mask = cap - 1u;

  uint64_t *ubits = (uint64_t *) R_alloc((size_t) N * (size_t) p, sizeof(uint64_t));
  for (int j = 0; j < p; ++j) {
    const double *col = Z + (size_t) j * (size_t) N;
    for (int i = 0; i < N; ++i) {
      ubits[(size_t) i * (size_t) p + (size_t) j] = cpt_canon_bits(col[i]);
    }
  }

  int *ht = (int *) R_alloc(cap, sizeof(int));
  memset(ht, 0, cap * sizeof(int));
  int *first_buf = (int *) R_alloc((size_t) N, sizeof(int));
  int *mult_buf = (int *) R_alloc((size_t) N, sizeof(int));

  SEXP map_sexp = PROTECT(Rf_allocVector(INTSXP, N));
  int *map_p = INTEGER(map_sexp);
  int Md = 0;
  size_t row_bytes = (size_t) p * sizeof(uint64_t);

  for (int i = 0; i < N; ++i) {
    const uint64_t *ri = ubits + (size_t) i * (size_t) p;
    uint64_t h = 0x9e3779b97f4a7c15ULL * (uint64_t) p;
    for (int j = 0; j < p; ++j) {
      h = cpt_mix64(h ^ (ri[j] + 0x9e3779b97f4a7c15ULL * (uint64_t) (j + 1)));
    }
    size_t slot = (size_t) h & mask;
    for (;;) {
      int d = ht[slot];
      if (d == 0) {
        d = ++Md;
        ht[slot] = d;
        first_buf[d - 1] = i + 1;
        mult_buf[d - 1] = 1;
        map_p[i] = d;
        break;
      }
      int r0 = first_buf[d - 1] - 1;
      const uint64_t *r0_bits = ubits + (size_t) r0 * (size_t) p;
      if (memcmp(ri, r0_bits, row_bytes) == 0) {
        map_p[i] = d;
        mult_buf[d - 1]++;
        break;
      }
      slot = (slot + 1u) & mask;
    }
  }

  SEXP z_out = PROTECT(Rf_allocMatrix(REALSXP, Md, p));
  double *z_dst = REAL(z_out);
  for (int j = 0; j < p; ++j) {
    const double *src_col = Z + (size_t) j * (size_t) N;
    double *dst_col = z_dst + (size_t) j * (size_t) Md;
    for (int d = 0; d < Md; ++d) {
      dst_col[d] = src_col[first_buf[d] - 1];
    }
  }

  SEXP first_sexp = PROTECT(Rf_allocVector(INTSXP, Md));
  SEXP mult_sexp = PROTECT(Rf_allocVector(INTSXP, Md));
  memcpy(INTEGER(first_sexp), first_buf, (size_t) Md * sizeof(int));
  memcpy(INTEGER(mult_sexp), mult_buf, (size_t) Md * sizeof(int));

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(nms, 0, Rf_mkChar("z"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("first"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("map"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("mult"));
  SET_STRING_ELT(nms, 4, Rf_mkChar("M"));
  SET_STRING_ELT(nms, 5, Rf_mkChar("n_distinct"));
  Rf_setAttrib(res, R_NamesSymbol, nms);
  SET_VECTOR_ELT(res, 0, z_out);
  SET_VECTOR_ELT(res, 1, first_sexp);
  SET_VECTOR_ELT(res, 2, map_sexp);
  SET_VECTOR_ELT(res, 3, mult_sexp);
  SET_VECTOR_ELT(res, 4, Rf_ScalarInteger(N));
  SET_VECTOR_ELT(res, 5, Rf_ScalarInteger(Md));

  vmaxset(vmax);
  UNPROTECT(6);
  return res;
}

/* ---------------------------------------------------------------------------
   Shared Helper: Orthant Dominance Contraction
   Computes Out[m * K + c] = sum_{a: zd[src_idx[a], .] <= zd[m, .]} S_act[a * K + c]
   for m = 0..Md-1, c = 0..K-1.
   Uses O(Md * K) for p_nc=0, O(Md log Md + Md * K) prefix sum for p_nc=1,
   O(Md log Md * K) 2D sweep + Fenwick tree for p_nc=2, and 64-target bitmask
   tiling for p_nc >= 3.
   --------------------------------------------------------------------------- */
typedef struct {
  double v;
  int idx;
} cpt_val_idx_t;

static int cpt_cmp_val_idx(const void *a, const void *b) {
  const cpt_val_idx_t *pa = (const cpt_val_idx_t *) a;
  const cpt_val_idx_t *pb = (const cpt_val_idx_t *) b;
  if (pa->v < pb->v) return -1;
  if (pa->v > pb->v) return 1;
  return (pa->idx < pb->idx) ? -1 : ((pa->idx > pb->idx) ? 1 : 0);
}

static void cpt_orthant_sum(
    const double *zd,
    int Md,
    int p,
    const int *src_idx,
    int M_act,
    const double *S_act,
    int K,
    double *Out) {
  memset(Out, 0, (size_t) Md * (size_t) K * sizeof(double));
  if (Md <= 0 || M_act <= 0 || K <= 0) return;

  const void *vmax = vmaxget();
  int *nc_cols = (p > 0) ? (int *) R_alloc((size_t) p, sizeof(int)) : NULL;
  int p_nc = 0;
  for (int j = 0; j < p; ++j) {
    const double *col = zd + (size_t) j * (size_t) Md;
    double zmin = col[0];
    double zmax = col[0];
    for (int m = 1; m < Md; ++m) {
      double v = col[m];
      if (v < zmin) zmin = v;
      if (v > zmax) zmax = v;
    }
    if (zmin < zmax) {
      nc_cols[p_nc++] = j;
    }
  }

  if (p_nc == 0) {
    double *tot = (double *) R_alloc((size_t) K, sizeof(double));
    memset(tot, 0, (size_t) K * sizeof(double));
    for (int a = 0; a < M_act; ++a) {
      const double *sv = S_act + (size_t) a * (size_t) K;
      for (int c = 0; c < K; ++c) tot[c] += sv[c];
    }
    for (int m = 0; m < Md; ++m) {
      memcpy(Out + (size_t) m * (size_t) K, tot, (size_t) K * sizeof(double));
    }
    vmaxset(vmax);
    return;
  }

  if (p_nc == 1) {
    const double *xcol = zd + (size_t) nc_cols[0] * (size_t) Md;
    int *is_act = (int *) R_alloc((size_t) Md, sizeof(int));
    for (int m = 0; m < Md; ++m) is_act[m] = -1;
    for (int a = 0; a < M_act; ++a) is_act[src_idx[a]] = a;

    cpt_val_idx_t *px = (cpt_val_idx_t *) R_alloc((size_t) Md, sizeof(cpt_val_idx_t));
    for (int m = 0; m < Md; ++m) {
      px[m].v = (xcol[m] == 0.0) ? 0.0 : xcol[m];
      px[m].idx = m;
    }
    qsort(px, (size_t) Md, sizeof(cpt_val_idx_t), cpt_cmp_val_idx);

    double *running = (double *) R_alloc((size_t) K, sizeof(double));
    memset(running, 0, (size_t) K * sizeof(double));

    int g0 = 0;
    while (g0 < Md) {
      double xv = px[g0].v;
      int g1 = g0 + 1;
      while (g1 < Md && px[g1].v == xv) ++g1;
      for (int t = g0; t < g1; ++t) {
        int a = is_act[px[t].idx];
        if (a >= 0) {
          const double *sv = S_act + (size_t) a * (size_t) K;
          for (int c = 0; c < K; ++c) running[c] += sv[c];
        }
      }
      for (int t = g0; t < g1; ++t) {
        memcpy(Out + (size_t) px[t].idx * (size_t) K, running, (size_t) K * sizeof(double));
      }
      g0 = g1;
    }
    vmaxset(vmax);
    return;
  }

  if (p_nc == 2) {
    const double *xcol = zd + (size_t) nc_cols[0] * (size_t) Md;
    const double *ycol = zd + (size_t) nc_cols[1] * (size_t) Md;

    int *is_act = (int *) R_alloc((size_t) Md, sizeof(int));
    for (int m = 0; m < Md; ++m) is_act[m] = -1;
    for (int a = 0; a < M_act; ++a) is_act[src_idx[a]] = a;

    cpt_val_idx_t *py = (cpt_val_idx_t *) R_alloc((size_t) Md, sizeof(cpt_val_idx_t));
    for (int m = 0; m < Md; ++m) {
      py[m].v = (ycol[m] == 0.0) ? 0.0 : ycol[m];
      py[m].idx = m;
    }
    qsort(py, (size_t) Md, sizeof(cpt_val_idx_t), cpt_cmp_val_idx);

    int *ry = (int *) R_alloc((size_t) Md, sizeof(int));
    int Ky = 1;
    ry[py[0].idx] = 1;
    for (int t = 1; t < Md; ++t) {
      if (py[t].v > py[t - 1].v) ++Ky;
      ry[py[t].idx] = Ky;
    }

    cpt_val_idx_t *px = py; /* reuse buffer */
    for (int m = 0; m < Md; ++m) {
      px[m].v = (xcol[m] == 0.0) ? 0.0 : xcol[m];
      px[m].idx = m;
    }
    qsort(px, (size_t) Md, sizeof(cpt_val_idx_t), cpt_cmp_val_idx);

    double *bit = (double *) R_alloc((size_t) (Ky + 1) * (size_t) K, sizeof(double));
    memset(bit, 0, (size_t) (Ky + 1) * (size_t) K * sizeof(double));

    int g0 = 0;
    while (g0 < Md) {
      double xv = px[g0].v;
      int g1 = g0 + 1;
      while (g1 < Md && px[g1].v == xv) ++g1;
      for (int t = g0; t < g1; ++t) {
        int i = px[t].idx;
        int a = is_act[i];
        if (a >= 0) {
          const double *sv = S_act + (size_t) a * (size_t) K;
          for (int idx = ry[i]; idx <= Ky; idx += (idx & -idx)) {
            double *bv = bit + (size_t) idx * (size_t) K;
            for (int c = 0; c < K; ++c) bv[c] += sv[c];
          }
        }
      }
      for (int t = g0; t < g1; ++t) {
        int m = px[t].idx;
        double *ov = Out + (size_t) m * (size_t) K;
        for (int idx = ry[m]; idx > 0; idx -= (idx & -idx)) {
          const double *bv = bit + (size_t) idx * (size_t) K;
          for (int c = 0; c < K; ++c) ov[c] += bv[c];
        }
      }
      g0 = g1;
    }
    vmaxset(vmax);
    return;
  }

  /* General p_nc >= 3: 64-target bitmask tiling over active sources */
  double *z_act = (double *) R_alloc((size_t) M_act * (size_t) p_nc, sizeof(double));
  for (int a = 0; a < M_act; ++a) {
    int d = src_idx[a];
    double *zr = z_act + (size_t) a * (size_t) p_nc;
    for (int j = 0; j < p_nc; ++j) {
      zr[j] = zd[(size_t) d + (size_t) nc_cols[j] * (size_t) Md];
    }
  }

  double *z_blk = (double *) R_alloc((size_t) p_nc * 64u, sizeof(double));
  for (int m0 = 0; m0 < Md; m0 += 64) {
    int nb = (Md - m0 < 64) ? (Md - m0) : 64;
    for (int j = 0; j < p_nc; ++j) {
      const double *col = zd + (size_t) nc_cols[j] * (size_t) Md + m0;
      memcpy(z_blk + (size_t) j * 64u, col, (size_t) nb * sizeof(double));
    }
    uint64_t init_mask = (nb == 64) ? ~0ULL : ((1ULL << nb) - 1ULL);
    for (int a = 0; a < M_act; ++a) {
      const double *za = z_act + (size_t) a * (size_t) p_nc;
      uint64_t mask = init_mask;
      for (int j = 0; j < p_nc; ++j) {
        double zav = za[j];
        const double *zt = z_blk + (size_t) j * 64u;
        uint64_t mj = 0ULL;
        for (int u = 0; u < nb; ++u) {
          mj |= ((uint64_t) (zav <= zt[u])) << u;
        }
        mask &= mj;
        if (!mask) break;
      }
      if (mask) {
        const double *sv = S_act + (size_t) a * (size_t) K;
        while (mask) {
          int u = cpt_ctz64(mask);
          mask &= mask - 1ULL;
          double *ov = Out + (size_t) (m0 + u) * (size_t) K;
          for (int c = 0; c < K; ++c) ov[c] += sv[c];
        }
      }
    }
  }
  vmaxset(vmax);
}

/* ---------------------------------------------------------------------------
   2. Native Cell Accumulator (init, piece add, fused two-arm/OR add, finalize)
   --------------------------------------------------------------------------- */
typedef struct {
  int Md;
  int R;
  int Dl;
  int C;
  int N_total;
  int keep_rows;
  double n_rows_seen;
  double *Vd;       /* Md x R row-major: Vd[d * R + r] */
  double *Pd;       /* Md x Dl row-major: Pd[d * Dl + c] */
  double *Bl;       /* Dl x R col-major: Bl[c + r * Dl] */
  double *AggX;     /* C x Dl col-major: AggX[k + c * C] */
  double *AggKu;    /* C x R col-major: AggKu[k + r * C] */
  double *sum_nu;   /* R */
  double *unit_A;   /* N_total */
  double *rows_nu;  /* N_total x R col-major (only if keep_rows) */
  double *rows_xi;  /* N_total x Dl col-major (only if keep_rows) */
} cpt_acc_t;

static void cpt_acc_free_struct(cpt_acc_t *st) {
  if (!st) return;
  if (st->Vd) R_Free(st->Vd);
  if (st->Pd) R_Free(st->Pd);
  if (st->Bl) R_Free(st->Bl);
  if (st->AggX) R_Free(st->AggX);
  if (st->AggKu) R_Free(st->AggKu);
  if (st->sum_nu) R_Free(st->sum_nu);
  if (st->unit_A) R_Free(st->unit_A);
  if (st->rows_nu) R_Free(st->rows_nu);
  if (st->rows_xi) R_Free(st->rows_xi);
  R_Free(st);
}

static void cpt_acc_finalizer(SEXP ptr) {
  cpt_acc_t *st = (cpt_acc_t *) R_ExternalPtrAddr(ptr);
  if (st) {
    cpt_acc_free_struct(st);
    R_ClearExternalPtr(ptr);
  }
}

static cpt_acc_t *cpt_acc_unwrap(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || R_ExternalPtrTag(ptr) != Rf_install("cdid_cpt_acc")) {
    return NULL;
  }
  return (cpt_acc_t *) R_ExternalPtrAddr(ptr);
}

SEXP C_cdid_cpt_acc_init(SEXP dims_sexp) {
  if (TYPEOF(dims_sexp) != INTSXP || XLENGTH(dims_sexp) != 6) return R_NilValue;
  const int *dm = INTEGER(dims_sexp);
  int Md = dm[0], R = dm[1], Dl = dm[2], C = dm[3], N_total = dm[4], keep_rows = dm[5];
  if (Md < 1 || R < 1 || Dl < 0 || C < 1 || N_total < 1 || (keep_rows != 0 && keep_rows != 1)) {
    return R_NilValue;
  }

  SEXP ptr = PROTECT(R_MakeExternalPtr(NULL, Rf_install("cdid_cpt_acc"), R_NilValue));
  R_RegisterCFinalizerEx(ptr, cpt_acc_finalizer, TRUE);

  cpt_acc_t *st = R_Calloc(1, cpt_acc_t);
  R_SetExternalPtrAddr(ptr, st);
  st->Md = Md;
  st->R = R;
  st->Dl = Dl;
  st->C = C;
  st->N_total = N_total;
  st->keep_rows = keep_rows;
  st->n_rows_seen = 0.0;

  st->Vd = R_Calloc((size_t) Md * (size_t) R, double);
  if (Dl > 0) {
    st->Pd = R_Calloc((size_t) Md * (size_t) Dl, double);
    st->Bl = R_Calloc((size_t) Dl * (size_t) R, double);
    st->AggX = R_Calloc((size_t) C * (size_t) Dl, double);
  }
  st->AggKu = R_Calloc((size_t) C * (size_t) R, double);
  st->sum_nu = R_Calloc((size_t) R, double);
  st->unit_A = R_Calloc((size_t) N_total, double);
  if (keep_rows) {
    st->rows_nu = R_Calloc((size_t) N_total * (size_t) R, double);
    if (Dl > 0) {
      st->rows_xi = R_Calloc((size_t) N_total * (size_t) Dl, double);
    }
  }

  UNPROTECT(1);
  return ptr;
}

static int cpt_check_finite_real(SEXP x, R_xlen_t expected_len) {
  if (Rf_isObject(x) || TYPEOF(x) != REALSXP || XLENGTH(x) != expected_len) return 0;
  const double *p = REAL(x);
  for (R_xlen_t i = 0; i < expected_len; ++i) {
    if (!R_FINITE(p[i])) return 0;
  }
  return 1;
}

static int cpt_check_finite_mat(SEXP x, int nr, int nc) {
  if (Rf_isObject(x) || !Rf_isMatrix(x) || TYPEOF(x) != REALSXP) return 0;
  if (Rf_nrows(x) != nr || Rf_ncols(x) != nc) return 0;
  const double *p = REAL(x);
  R_xlen_t tot = (R_xlen_t) nr * (R_xlen_t) nc;
  for (R_xlen_t i = 0; i < tot; ++i) {
    if (!R_FINITE(p[i])) return 0;
  }
  return 1;
}

SEXP C_cdid_cpt_acc_add_piece(
    SEXP ptr_sexp,
    SEXP rows_sexp,
    SEXP nu_sexp,
    SEXP rho_sexp,
    SEXP P_sexp,
    SEXP xi_sexp,
    SEXP blam_sexp,
    SEXP sD_sexp,
    SEXP gmap_sexp,
    SEXP cmap_sexp) {
  cpt_acc_t *st = cpt_acc_unwrap(ptr_sexp);
  if (!st) return R_NilValue;
  if (TYPEOF(rows_sexp) != INTSXP || TYPEOF(gmap_sexp) != INTSXP || TYPEOF(cmap_sexp) != INTSXP) {
    return R_NilValue;
  }
  int n = (int) XLENGTH(rows_sexp);
  if (n == 0) return Rf_ScalarLogical(1);
  if (XLENGTH(gmap_sexp) != st->N_total || XLENGTH(cmap_sexp) != st->N_total) return R_NilValue;

  int R = st->R, Dl = st->Dl, C = st->C, N_total = st->N_total, Md = st->Md;
  if (!cpt_check_finite_mat(nu_sexp, n, R) ||
      !cpt_check_finite_mat(rho_sexp, n, R) ||
      !cpt_check_finite_mat(P_sexp, n, Dl) ||
      !cpt_check_finite_mat(xi_sexp, n, Dl) ||
      !cpt_check_finite_mat(blam_sexp, Dl, R) ||
      !cpt_check_finite_real(sD_sexp, R)) {
    return Rf_ScalarLogical(0);
  }

  const int *rows = INTEGER(rows_sexp);
  const int *gmap = INTEGER(gmap_sexp);
  const int *cmap = INTEGER(cmap_sexp);
  for (int i = 0; i < n; ++i) {
    int r0 = rows[i] - 1;
    if (r0 < 0 || r0 >= N_total) return R_NilValue;
    int d0 = gmap[r0] - 1;
    int k0 = cmap[r0] - 1;
    if (d0 < 0 || d0 >= Md || k0 < 0 || k0 >= C) return R_NilValue;
  }

  const double *nu = REAL(nu_sexp);
  const double *rho = REAL(rho_sexp);
  const double *P = REAL(P_sexp);
  const double *xi = REAL(xi_sexp);
  const double *blam = REAL(blam_sexp);
  const double *sD = REAL(sD_sexp);

  if (Dl > 0) {
    size_t bl_len = (size_t) Dl * (size_t) R;
    for (size_t idx = 0; idx < bl_len; ++idx) st->Bl[idx] += blam[idx];
  }

  for (int r = 0; r < R; ++r) {
    const double *nu_col = nu + (size_t) r * (size_t) n;
    const double *rho_col = rho + (size_t) r * (size_t) n;
    double *agg_ku_col = st->AggKu + (size_t) r * (size_t) C;
    double *rnu_col = st->keep_rows ? (st->rows_nu + (size_t) r * (size_t) N_total) : NULL;
    double sDr = sD[r];
    double s_nu = 0.0;
    for (int i = 0; i < n; ++i) {
      int r0 = rows[i] - 1;
      int d0 = gmap[r0] - 1;
      int k0 = cmap[r0] - 1;
      double nval = nu_col[i];
      double nr = nval * rho_col[i];
      st->Vd[(size_t) d0 * (size_t) R + (size_t) r] += nr;
      agg_ku_col[k0] += nval;
      s_nu += nval;
      st->unit_A[r0] += nr * sDr;
      if (rnu_col) rnu_col[r0] += nval;
    }
    st->sum_nu[r] += s_nu;
  }

  for (int c = 0; c < Dl; ++c) {
    const double *p_col = P + (size_t) c * (size_t) n;
    const double *xi_col = xi + (size_t) c * (size_t) n;
    double *agg_x_col = st->AggX + (size_t) c * (size_t) C;
    double *rxi_col = st->keep_rows ? (st->rows_xi + (size_t) c * (size_t) N_total) : NULL;
    for (int i = 0; i < n; ++i) {
      int r0 = rows[i] - 1;
      double pv = p_col[i];
      if (pv != 0.0) st->Pd[(size_t) (gmap[r0] - 1) * (size_t) Dl + (size_t) c] += pv;
      double xv = xi_col[i];
      if (xv != 0.0) agg_x_col[cmap[r0] - 1] += xv;
      if (rxi_col) rxi_col[r0] += xv;
    }
  }

  st->n_rows_seen += (double) n;
  return Rf_ScalarLogical(1);
}

/* Fused zero-allocation piece + scatter-accumulator for two-arm (mode=1) and OR (mode=2) */
SEXP C_cdid_cpt_acc_add_fast(
    SEXP ptr_sexp,
    SEXP mode_sexp,
    SEXP rows_sexp,
    SEXP rr_sexp,
    SEXP obs_vecs_sexp,
    SEXP obs_mats_sexp,
    SEXP xi_blocks_sexp,
    SEXP params_sexp,
    SEXP gmap_sexp,
    SEXP cmap_sexp) {
  cpt_acc_t *st = cpt_acc_unwrap(ptr_sexp);
  if (!st) return R_NilValue;
  if (TYPEOF(mode_sexp) != INTSXP || XLENGTH(mode_sexp) != 1 ||
      TYPEOF(rows_sexp) != INTSXP || TYPEOF(rr_sexp) != INTSXP ||
      TYPEOF(obs_vecs_sexp) != VECSXP || TYPEOF(obs_mats_sexp) != VECSXP ||
      TYPEOF(xi_blocks_sexp) != VECSXP || TYPEOF(params_sexp) != REALSXP ||
      TYPEOF(gmap_sexp) != INTSXP || TYPEOF(cmap_sexp) != INTSXP) {
    return R_NilValue;
  }

  int mode = INTEGER(mode_sexp)[0];
  int n = (int) XLENGTH(rows_sexp);
  int m = (int) XLENGTH(rr_sexp);
  int Dl = st->Dl, C = st->C, N_total = st->N_total, Md = st->Md;
  if (n == 0) return Rf_ScalarLogical(1);
  if (XLENGTH(gmap_sexp) != N_total || XLENGTH(cmap_sexp) != N_total) return R_NilValue;

  const int *rows = INTEGER(rows_sexp);
  const int *rr = INTEGER(rr_sexp);
  const int *gmap = INTEGER(gmap_sexp);
  const int *cmap = INTEGER(cmap_sexp);

  for (int i = 0; i < n; ++i) {
    int r0 = rows[i] - 1;
    if (r0 < 0 || r0 >= N_total) return R_NilValue;
    int d0 = gmap[r0] - 1;
    int k0 = cmap[r0] - 1;
    if (d0 < 0 || d0 >= Md || k0 < 0 || k0 >= C) return R_NilValue;
  }
  for (int u = 0; u < m; ++u) {
    if (rr[u] < 1 || rr[u] > n) return R_NilValue;
  }

  /* Validate xi_blocks first (width must sum to Dl and all entries finite) */
  int n_xb = (int) XLENGTH(xi_blocks_sexp);
  int xi_cols_seen = 0;
  for (int b = 0; b < n_xb; ++b) {
    SEXP blk = VECTOR_ELT(xi_blocks_sexp, b);
    if (blk == R_NilValue) continue;
    if (Rf_isObject(blk) || TYPEOF(blk) != REALSXP) return Rf_ScalarLogical(0);
    int wb = 0;
    if (Rf_isMatrix(blk)) {
      if (Rf_nrows(blk) != n) return Rf_ScalarLogical(0);
      wb = Rf_ncols(blk);
    } else {
      if (XLENGTH(blk) != n) return Rf_ScalarLogical(0);
      wb = 1;
    }
    const double *bp = REAL(blk);
    R_xlen_t tot = (R_xlen_t) n * (R_xlen_t) wb;
    for (R_xlen_t idx = 0; idx < tot; ++idx) {
      if (!R_FINITE(bp[idx])) return Rf_ScalarLogical(0);
    }
    xi_cols_seen += wb;
  }
  if (xi_cols_seen != Dl) return Rf_ScalarLogical(0);

  if (mode == 1) {
    /* Two-arm (IPW / DR): R must be 2 */
    if (st->R != 2 || XLENGTH(obs_vecs_sexp) != 5 || XLENGTH(obs_mats_sexp) != 3 || XLENGTH(params_sexp) != 5) {
      return Rf_ScalarLogical(0);
    }
    const double *prm = REAL(params_sexp);
    double D_T = prm[0], D_M = prm[1];
    int n_up = (int) prm[2], p_trt = (int) prm[3], or_icpt = (int) prm[4];
    if (!R_FINITE(D_T) || !R_FINITE(D_M) || D_T == 0.0 || D_M == 0.0) return Rf_ScalarLogical(0);
    double sD0 = 1.0 / D_T;
    double sD1 = -1.0 / D_M;
    if (!R_FINITE(sD0) || !R_FINITE(sD1)) return Rf_ScalarLogical(0);

    SEXP D_s = VECTOR_ELT(obs_vecs_sexp, 0);
    SEXP fw_s = VECTOR_ELT(obs_vecs_sexp, 1);
    SEXP odds_s = VECTOR_ELT(obs_vecs_sexp, 2);
    SEXP act_s = VECTOR_ELT(obs_vecs_sexp, 3);
    SEXP rho_s = VECTOR_ELT(obs_vecs_sexp, 4);
    SEXP lam_up_s = VECTOR_ELT(obs_mats_sexp, 0);
    SEXP X_s = VECTOR_ELT(obs_mats_sexp, 1);
    SEXP d_or_s = VECTOR_ELT(obs_mats_sexp, 2);

    int p_or = (d_or_s == R_NilValue) ? 0 : (Rf_isMatrix(d_or_s) ? Rf_ncols(d_or_s) : -1);
    if (p_or < 0 || n_up + p_trt + (or_icpt ? 1 : 0) + p_or != Dl) return Rf_ScalarLogical(0);

    if (m > 0) {
      if (!cpt_check_finite_real(D_s, m) ||
          !cpt_check_finite_real(fw_s, m) ||
          !cpt_check_finite_real(odds_s, m) ||
          !cpt_check_finite_real(act_s, m) ||
          !cpt_check_finite_real(rho_s, m)) return Rf_ScalarLogical(0);
      if (n_up > 0 && !cpt_check_finite_mat(lam_up_s, m, n_up)) return Rf_ScalarLogical(0);
      if (p_trt > 0 && !cpt_check_finite_mat(X_s, m, p_trt)) return Rf_ScalarLogical(0);
      if (p_or > 0 && !cpt_check_finite_mat(d_or_s, m, p_or)) return Rf_ScalarLogical(0);

      const double *Dv = REAL(D_s);
      const double *fw = REAL(fw_s);
      const double *odds = REAL(odds_s);
      const double *act = REAL(act_s);
      const double *rho = REAL(rho_s);
      const double *lam_up = (n_up > 0) ? REAL(lam_up_s) : NULL;
      const double *X = (p_trt > 0) ? REAL(X_s) : NULL;
      const double *d_or = (p_or > 0) ? REAL(d_or_s) : NULL;

      /* Pre-check finiteness of computed products before mutating state */
      for (int u = 0; u < m; ++u) {
        double nu0 = fw[u] * Dv[u];
        double nu1 = fw[u] * (1.0 - Dv[u]) * odds[u];
        double nr0 = nu0 * rho[u];
        double nr1 = nu1 * rho[u];
        if (!R_FINITE(nu0) || !R_FINITE(nu1) || !R_FINITE(nr0) || !R_FINITE(nr1)) {
          return Rf_ScalarLogical(0);
        }
      }

      double *agg_ku0 = st->AggKu;
      double *agg_ku1 = st->AggKu + (size_t) C;
      double *Bl0 = st->Bl;
      double *Bl1 = st->Bl + (size_t) Dl;
      double s_nu0 = 0.0, s_nu1 = 0.0;

      for (int u = 0; u < m; ++u) {
        int r0 = rows[rr[u] - 1] - 1;
        int d0 = gmap[r0] - 1;
        int k0 = cmap[r0] - 1;

        double nu0 = fw[u] * Dv[u];
        double nu1 = fw[u] * (1.0 - Dv[u]) * odds[u];
        double rh = rho[u];
        double nr0 = nu0 * rh;
        double nr1 = nu1 * rh;

        st->Vd[(size_t) d0 * 2u + 0u] += nr0;
        st->Vd[(size_t) d0 * 2u + 1u] += nr1;
        agg_ku0[k0] += nu0;
        agg_ku1[k0] += nu1;
        s_nu0 += nu0;
        s_nu1 += nu1;

        double wnr0 = sD0 * nr0;
        double wnr1 = sD1 * nr1;
        st->unit_A[r0] += nr0 * sD0 + nr1 * sD1;

        if (st->keep_rows) {
          st->rows_nu[r0] += nu0;
          st->rows_nu[r0 + N_total] += nu1;
        }

        if (Dl > 0) {
          double *pd_row = st->Pd + (size_t) d0 * (size_t) Dl;
          for (int c = 0; c < n_up; ++c) {
            double lu = lam_up[(size_t) u + (size_t) c * (size_t) m];
            if (lu != 0.0) {
              pd_row[c] += wnr0 * lu + wnr1 * lu;
              Bl0[c] += nu0 * lu;
              Bl1[c] += nu1 * lu;
            }
          }
          double au = act[u];
          if (au != 0.0) {
            for (int c = 0; c < p_trt; ++c) {
              double lp = au * X[(size_t) u + (size_t) c * (size_t) m];
              if (lp != 0.0) {
                int col = n_up + c;
                pd_row[col] += wnr1 * lp;
                Bl1[col] += nu1 * lp;
              }
            }
          }
          int off = n_up + p_trt;
          double wnu0 = sD0 * nu0;
          double wnu1 = sD1 * nu1;
          if (or_icpt) {
            pd_row[off] += -wnu0 - wnu1;
            ++off;
          }
          for (int c = 0; c < p_or; ++c) {
            double dv = -d_or[(size_t) u + (size_t) c * (size_t) m];
            if (dv != 0.0) {
              pd_row[off + c] += wnu0 * dv + wnu1 * dv;
            }
          }
        }
      }
      st->sum_nu[0] += s_nu0;
      st->sum_nu[1] += s_nu1;
    }
  } else if (mode == 2) {
    /* Control-only OR: R must be 1 */
    if (st->R != 1 || XLENGTH(obs_vecs_sexp) != 2 || XLENGTH(obs_mats_sexp) != 2 || XLENGTH(params_sexp) != 3) {
      return Rf_ScalarLogical(0);
    }
    const double *prm = REAL(params_sexp);
    double d_denom = prm[0];
    int n_ch = (int) prm[1], or_icpt = (int) prm[2];
    if (!R_FINITE(d_denom) || d_denom == 0.0) return Rf_ScalarLogical(0);
    double sD0 = 1.0 / d_denom;
    if (!R_FINITE(sD0)) return Rf_ScalarLogical(0);

    SEXP q_s = VECTOR_ELT(obs_vecs_sexp, 0);
    SEXP rho_s = VECTOR_ELT(obs_vecs_sexp, 1);
    SEXP lam_s = VECTOR_ELT(obs_mats_sexp, 0);
    SEXP d_or_s = VECTOR_ELT(obs_mats_sexp, 1);

    int p_or = (d_or_s == R_NilValue) ? 0 : (Rf_isMatrix(d_or_s) ? Rf_ncols(d_or_s) : -1);
    if (p_or < 0 || n_ch + (or_icpt ? 1 : 0) + p_or != Dl) return Rf_ScalarLogical(0);

    if (m > 0) {
      if (!cpt_check_finite_real(q_s, m) || !cpt_check_finite_real(rho_s, m)) return Rf_ScalarLogical(0);
      if (n_ch > 0 && !cpt_check_finite_mat(lam_s, m, n_ch)) return Rf_ScalarLogical(0);
      if (p_or > 0 && !cpt_check_finite_mat(d_or_s, m, p_or)) return Rf_ScalarLogical(0);

      const double *qv = REAL(q_s);
      const double *rho = REAL(rho_s);
      const double *lam = (n_ch > 0) ? REAL(lam_s) : NULL;
      const double *d_or = (p_or > 0) ? REAL(d_or_s) : NULL;

      for (int u = 0; u < m; ++u) {
        if (!R_FINITE(qv[u] * rho[u])) return Rf_ScalarLogical(0);
      }

      double s_nu0 = 0.0;
      for (int u = 0; u < m; ++u) {
        int r0 = rows[rr[u] - 1] - 1;
        int d0 = gmap[r0] - 1;
        int k0 = cmap[r0] - 1;

        double nu0 = qv[u];
        double nr0 = nu0 * rho[u];
        st->Vd[d0] += nr0;
        st->AggKu[k0] += nu0;
        s_nu0 += nu0;
        double wnr0 = sD0 * nr0;
        st->unit_A[r0] += nr0 * sD0;
        if (st->keep_rows) st->rows_nu[r0] += nu0;

        if (Dl > 0) {
          double *pd_row = st->Pd + (size_t) d0 * (size_t) Dl;
          for (int c = 0; c < n_ch; ++c) {
            double lc = lam[(size_t) u + (size_t) c * (size_t) m];
            if (lc != 0.0) {
              pd_row[c] += wnr0 * lc;
              st->Bl[c] += nu0 * lc;
            }
          }
          int off = n_ch;
          double wnu0 = sD0 * nu0;
          if (or_icpt) {
            pd_row[off] += -wnu0;
            ++off;
          }
          for (int c = 0; c < p_or; ++c) {
            double dv = -d_or[(size_t) u + (size_t) c * (size_t) m];
            if (dv != 0.0) pd_row[off + c] += wnu0 * dv;
          }
        }
      }
      st->sum_nu[0] += s_nu0;
    }
  } else {
    return R_NilValue;
  }

  /* Scatter-add xi_blocks across all n chunk rows */
  int col_off = 0;
  for (int b = 0; b < n_xb; ++b) {
    SEXP blk = VECTOR_ELT(xi_blocks_sexp, b);
    if (blk == R_NilValue) continue;
    int wb = Rf_isMatrix(blk) ? Rf_ncols(blk) : 1;
    const double *bp = REAL(blk);
    for (int cb = 0; cb < wb; ++cb) {
      int c = col_off + cb;
      const double *xcol = bp + (size_t) cb * (size_t) n;
      double *agg_col = st->AggX + (size_t) c * (size_t) C;
      double *rxi_col = st->keep_rows ? (st->rows_xi + (size_t) c * (size_t) N_total) : NULL;
      for (int i = 0; i < n; ++i) {
        int r0 = rows[i] - 1;
        double xv = xcol[i];
        if (xv != 0.0) agg_col[cmap[r0] - 1] += xv;
        if (rxi_col) rxi_col[r0] += xv;
      }
    }
    col_off += wb;
  }

  st->n_rows_seen += (double) n;
  return Rf_ScalarLogical(1);
}

typedef struct {
  int64_t key;
  int k;
  int d;
  double val;
} cpt_kd_entry_t;

static int cpt_cmp_kd(const void *a, const void *b) {
  const cpt_kd_entry_t *pa = (const cpt_kd_entry_t *) a;
  const cpt_kd_entry_t *pb = (const cpt_kd_entry_t *) b;
  return (pa->key < pb->key) ? -1 : ((pa->key > pb->key) ? 1 : 0);
}

SEXP C_cdid_cpt_acc_finalize(
    SEXP ptr_sexp,
    SEXP zd_sexp,
    SEXP mult_sexp,
    SEXP gmap_sexp,
    SEXP cmap_sexp,
    SEXP spec_num_sexp) {
  cpt_acc_t *st = cpt_acc_unwrap(ptr_sexp);
  if (!st) return R_NilValue;
  if (!Rf_isMatrix(zd_sexp) || TYPEOF(zd_sexp) != REALSXP ||
      TYPEOF(mult_sexp) != INTSXP || TYPEOF(gmap_sexp) != INTSXP || TYPEOF(cmap_sexp) != INTSXP ||
      TYPEOF(spec_num_sexp) != VECSXP || XLENGTH(spec_num_sexp) != 5) {
    return R_NilValue;
  }

  int Md = st->Md, R = st->R, Dl = st->Dl, C = st->C, N_total = st->N_total;
  int p = Rf_ncols(zd_sexp);
  if (Rf_nrows(zd_sexp) != Md || XLENGTH(mult_sexp) != Md ||
      XLENGTH(gmap_sexp) != N_total || XLENGTH(cmap_sexp) != N_total) {
    return R_NilValue;
  }

  double sigma = Rf_asReal(VECTOR_ELT(spec_num_sexp, 0));
  SEXP s_sexp = VECTOR_ELT(spec_num_sexp, 1);
  SEXP D_sexp = VECTOR_ELT(spec_num_sexp, 2);
  int M_grid = Rf_asInteger(VECTOR_ELT(spec_num_sexp, 3));
  double N_dbl = Rf_asReal(VECTOR_ELT(spec_num_sexp, 4));
  if (TYPEOF(s_sexp) != REALSXP || XLENGTH(s_sexp) != R ||
      TYPEOF(D_sexp) != REALSXP || XLENGTH(D_sexp) != R || N_dbl <= 0.0) {
    return R_NilValue;
  }
  const double *s_vec = REAL(s_sexp);
  const double *D_vec = REAL(D_sexp);

  const void *vmax = vmaxget();
  double *sD = (double *) R_alloc((size_t) R, sizeof(double));
  double *invD = (double *) R_alloc((size_t) R, sizeof(double));
  for (int r = 0; r < R; ++r) {
    invD[r] = 1.0 / D_vec[r];
    sD[r] = s_vec[r] / D_vec[r];
  }

  /* Pack active rows of [Vn = Vd / D, Pd] into row-major S_act (M_act x K_tot) */
  int K_tot = R + Dl;
  int *src_idx = (int *) R_alloc((size_t) Md, sizeof(int));
  int M_act = 0;
  double *tot_Vn = (double *) R_alloc((size_t) R, sizeof(double));
  memset(tot_Vn, 0, (size_t) R * sizeof(double));
  double *tot_Pd = (Dl > 0) ? (double *) R_alloc((size_t) Dl, sizeof(double)) : NULL;
  if (Dl > 0) memset(tot_Pd, 0, (size_t) Dl * sizeof(double));

  for (int d = 0; d < Md; ++d) {
    int nz = 0;
    const double *vd_r = st->Vd + (size_t) d * (size_t) R;
    for (int r = 0; r < R; ++r) {
      if (vd_r[r] != 0.0) {
        nz = 1;
        tot_Vn[r] += vd_r[r] * invD[r];
      }
    }
    if (Dl > 0) {
      const double *pd_r = st->Pd + (size_t) d * (size_t) Dl;
      for (int c = 0; c < Dl; ++c) {
        if (pd_r[c] != 0.0) {
          nz = 1;
          tot_Pd[c] += pd_r[c];
        }
      }
    }
    if (nz) src_idx[M_act++] = d;
  }

  double *S_act = (M_act > 0) ? (double *) R_alloc((size_t) M_act * (size_t) K_tot, sizeof(double)) : NULL;
  for (int a = 0; a < M_act; ++a) {
    int d = src_idx[a];
    double *dst = S_act + (size_t) a * (size_t) K_tot;
    const double *vd_r = st->Vd + (size_t) d * (size_t) R;
    for (int r = 0; r < R; ++r) dst[r] = vd_r[r] * invD[r];
    if (Dl > 0) {
      memcpy(dst + R, st->Pd + (size_t) d * (size_t) Dl, (size_t) Dl * sizeof(double));
    }
  }

  double *Out = (double *) R_alloc((size_t) Md * (size_t) K_tot, sizeof(double));
  cpt_orthant_sum(REAL(zd_sexp), Md, p, src_idx, M_act, S_act, K_tot, Out);

  /* Allocate R outputs */
  SEXP Rmat_sexp = PROTECT(Rf_allocMatrix(REALSXP, Md, R));
  SEXP Lmat_sexp = PROTECT(Rf_allocMatrix(REALSXP, Dl, Md));
  SEXP J_sexp = PROTECT(Rf_allocVector(REALSXP, Md));
  SEXP R1_sexp = PROTECT(Rf_allocVector(REALSXP, R));
  SEXP L1_sexp = PROTECT(Rf_allocVector(REALSXP, Dl));
  SEXP mult_dbl_sexp = PROTECT(Rf_allocVector(REALSXP, Md));
  SEXP AggX_sexp = PROTECT(Rf_allocMatrix(REALSXP, C, Dl));
  SEXP AggK_sexp = PROTECT(Rf_allocMatrix(REALSXP, C, R));
  SEXP accum_sexp = PROTECT(Rf_allocVector(REALSXP, R));

  double *Rmat = REAL(Rmat_sexp);
  double *Lmat = REAL(Lmat_sexp);
  double *J = REAL(J_sexp);
  double invN = 1.0 / N_dbl;

  for (int m = 0; m < Md; ++m) {
    const double *om = Out + (size_t) m * (size_t) K_tot;
    double jm = 0.0;
    for (int r = 0; r < R; ++r) {
      double rc_mr = om[r] / N_dbl;
      Rmat[(size_t) m + (size_t) r * (size_t) Md] = rc_mr;
      jm += rc_mr * s_vec[r];
    }
    J[m] = sigma * jm;
    if (Dl > 0) {
      double *lm_col = Lmat + (size_t) m * (size_t) Dl;
      const double *ph_m = om + R;
      for (int c = 0; c < Dl; ++c) {
        double sub = 0.0;
        for (int r = 0; r < R; ++r) {
          sub += st->Bl[(size_t) c + (size_t) r * (size_t) Dl] *
                 (Rmat[(size_t) m + (size_t) r * (size_t) Md] * sD[r]);
        }
        lm_col[c] = (ph_m[c] - sub) * invN;
      }
    }
  }

  double *R1 = REAL(R1_sexp);
  double mu1 = 0.0;
  for (int r = 0; r < R; ++r) {
    R1[r] = tot_Vn[r] / N_dbl;
    mu1 += s_vec[r] * R1[r];
    REAL(accum_sexp)[r] = st->sum_nu[r] / N_dbl;
  }

  if (Dl > 0) {
    double *L1 = REAL(L1_sexp);
    for (int c = 0; c < Dl; ++c) {
      double sub = 0.0;
      for (int r = 0; r < R; ++r) {
        sub += st->Bl[(size_t) c + (size_t) r * (size_t) Dl] * (R1[r] * sD[r]);
      }
      L1[c] = (tot_Pd[c] - sub) * invN;
    }
    memcpy(REAL(AggX_sexp), st->AggX, (size_t) C * (size_t) Dl * sizeof(double));
  }

  double *AggK = REAL(AggK_sexp);
  for (int r = 0; r < R; ++r) {
    const double *src = st->AggKu + (size_t) r * (size_t) C;
    double *dst = AggK + (size_t) r * (size_t) C;
    double sDr = sD[r];
    for (int k = 0; k < C; ++k) dst[k] = src[k] * sDr;
  }

  const int *mult_int = INTEGER(mult_sexp);
  double *mult_dbl = REAL(mult_dbl_sexp);
  for (int m = 0; m < Md; ++m) mult_dbl[m] = (double) mult_int[m];

  /* Group nonzero unit_A[0..N_total-1] by (k = cmap[i], d = gmap[i]) */
  const int *gmap = INTEGER(gmap_sexp);
  const int *cmap = INTEGER(cmap_sexp);
  int nnz_u = 0;
  for (int i = 0; i < N_total; ++i) {
    if (st->unit_A[i] != 0.0) ++nnz_u;
  }

  int n_kd_keep = 0;
  cpt_kd_entry_t *entries = NULL;
  if (nnz_u > 0) {
    size_t hcap = cpt_next_pow2((size_t) nnz_u * 2u + 8u);
    size_t hmask = hcap - 1u;
    int *ht = (int *) R_alloc(hcap, sizeof(int));
    memset(ht, 0, hcap * sizeof(int));
    entries = (cpt_kd_entry_t *) R_alloc((size_t) nnz_u, sizeof(cpt_kd_entry_t));
    int n_kd = 0;

    for (int i = 0; i < N_total; ++i) {
      double av = st->unit_A[i];
      if (av == 0.0) continue;
      int k = cmap[i];
      int d = gmap[i];
      int64_t key = (int64_t) (k - 1) * (int64_t) Md + (int64_t) d;
      uint64_t h = cpt_mix64((uint64_t) key);
      size_t slot = (size_t) h & hmask;
      for (;;) {
        int idx1 = ht[slot];
        if (idx1 == 0) {
          entries[n_kd].key = key;
          entries[n_kd].k = k;
          entries[n_kd].d = d;
          entries[n_kd].val = av;
          ++n_kd;
          ht[slot] = n_kd;
          break;
        }
        if (entries[idx1 - 1].key == key) {
          entries[idx1 - 1].val += av;
          break;
        }
        slot = (slot + 1u) & hmask;
      }
    }

    for (int e = 0; e < n_kd; ++e) {
      if (entries[e].val != 0.0) {
        if (n_kd_keep != e) entries[n_kd_keep] = entries[e];
        ++n_kd_keep;
      }
    }
    if (n_kd_keep > 1) {
      qsort(entries, (size_t) n_kd_keep, sizeof(cpt_kd_entry_t), cpt_cmp_kd);
    }
  }

  SEXP Ak_sexp = PROTECT(Rf_allocVector(INTSXP, n_kd_keep));
  SEXP Ad_sexp = PROTECT(Rf_allocVector(INTSXP, n_kd_keep));
  SEXP Av_sexp = PROTECT(Rf_allocVector(REALSXP, n_kd_keep));
  int *Ak = INTEGER(Ak_sexp);
  int *Ad = INTEGER(Ad_sexp);
  double *Av = REAL(Av_sexp);
  for (int e = 0; e < n_kd_keep; ++e) {
    Ak[e] = entries[e].k;
    Ad[e] = entries[e].d;
    Av[e] = entries[e].val;
  }

  SEXP rows_list = R_NilValue;
  int nprot = 12;
  if (st->keep_rows) {
    rows_list = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP rnms = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(rnms, 0, Rf_mkChar("A"));
    SET_STRING_ELT(rnms, 1, Rf_mkChar("kappa"));
    SET_STRING_ELT(rnms, 2, Rf_mkChar("xi"));
    SET_STRING_ELT(rnms, 3, Rf_mkChar("map"));
    Rf_setAttrib(rows_list, R_NamesSymbol, rnms);

    SEXP rA_sexp = PROTECT(Rf_allocVector(REALSXP, N_total));
    SEXP rK_sexp = PROTECT(Rf_allocMatrix(REALSXP, N_total, R));
    SEXP rX_sexp = PROTECT(Rf_allocMatrix(REALSXP, N_total, Dl));
    memcpy(REAL(rA_sexp), st->unit_A, (size_t) N_total * sizeof(double));
    double *rK = REAL(rK_sexp);
    for (int r = 0; r < R; ++r) {
      const double *src = st->rows_nu + (size_t) r * (size_t) N_total;
      double *dst = rK + (size_t) r * (size_t) N_total;
      double sDr = sD[r];
      for (int i = 0; i < N_total; ++i) dst[i] = src[i] * sDr;
    }
    if (Dl > 0) {
      memcpy(REAL(rX_sexp), st->rows_xi, (size_t) N_total * (size_t) Dl * sizeof(double));
    }
    SET_VECTOR_ELT(rows_list, 0, rA_sexp);
    SET_VECTOR_ELT(rows_list, 1, rK_sexp);
    SET_VECTOR_ELT(rows_list, 2, rX_sexp);
    SET_VECTOR_ELT(rows_list, 3, gmap_sexp);
    nprot += 5;
  }

  /* Release C accumulator memory immediately */
  cpt_acc_free_struct(st);
  R_ClearExternalPtr(ptr_sexp);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 15));
  SEXP onms = PROTECT(Rf_allocVector(STRSXP, 15));
  const char *field_names[15] = {
    "J", "mult", "M", "Rmat", "Lmat", "AggX", "AggK",
    "A_k", "A_d", "A_v", "mu1", "L1", "R1", "accumulated", "rows"
  };
  for (int f = 0; f < 15; ++f) SET_STRING_ELT(onms, f, Rf_mkChar(field_names[f]));
  Rf_setAttrib(out, R_NamesSymbol, onms);
  SET_VECTOR_ELT(out, 0, J_sexp);
  SET_VECTOR_ELT(out, 1, mult_dbl_sexp);
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(M_grid));
  SET_VECTOR_ELT(out, 3, Rmat_sexp);
  SET_VECTOR_ELT(out, 4, Lmat_sexp);
  SET_VECTOR_ELT(out, 5, AggX_sexp);
  SET_VECTOR_ELT(out, 6, AggK_sexp);
  SET_VECTOR_ELT(out, 7, Ak_sexp);
  SET_VECTOR_ELT(out, 8, Ad_sexp);
  SET_VECTOR_ELT(out, 9, Av_sexp);
  SET_VECTOR_ELT(out, 10, Rf_ScalarReal(mu1));
  SET_VECTOR_ELT(out, 11, L1_sexp);
  SET_VECTOR_ELT(out, 12, R1_sexp);
  SET_VECTOR_ELT(out, 13, accum_sexp);
  SET_VECTOR_ELT(out, 14, rows_list);
  nprot += 2;

  vmaxset(vmax);
  UNPROTECT(nprot);
  return out;
}

/* Helper: lookup named element in an R list */
static SEXP cpt_list_elt(SEXP list, const char *name) {
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if (TYPEOF(list) != VECSXP || TYPEOF(names) != STRSXP) return R_NilValue;
  R_xlen_t n = XLENGTH(list);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

/* ---------------------------------------------------------------------------
   3. Fast Multiplier Bootstrap Batch Kernel: cpt_batch_cvmb (serial & fork)
   --------------------------------------------------------------------------- */
SEXP C_cdid_cpt_batch_cvmb(SEXP E_sexp, SEXP cells_sexp, SEXP N_sexp) {
  if (Rf_isObject(E_sexp) || !Rf_isMatrix(E_sexp) || TYPEOF(E_sexp) != REALSXP ||
      TYPEOF(cells_sexp) != VECSXP) {
    return R_NilValue;
  }
  int b = Rf_nrows(E_sexp);
  int C = Rf_ncols(E_sexp);
  double N_dbl = Rf_asReal(N_sexp);
  if (b < 1 || C < 1 || !R_FINITE(N_dbl) || N_dbl <= 0.0) return R_NilValue;

  const double *E = REAL(E_sexp);
  int n_cells = (int) XLENGTH(cells_sexp);

  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, b));
  double *out = REAL(out_sexp);
  memset(out, 0, (size_t) b * sizeof(double));
  double invN = 1.0 / N_dbl;

  for (int ci = 0; ci < n_cells; ++ci) {
    const void *vmax = vmaxget();
    SEXP cl = VECTOR_ELT(cells_sexp, ci);
    SEXP zd_s = cpt_list_elt(cl, "zd");
    SEXP mult_s = cpt_list_elt(cl, "mult");
    SEXP M_s = cpt_list_elt(cl, "M");
    SEXP Rmat_s = cpt_list_elt(cl, "Rmat");
    SEXP Lmat_s = cpt_list_elt(cl, "Lmat");
    SEXP AggX_s = cpt_list_elt(cl, "AggX");
    SEXP AggK_s = cpt_list_elt(cl, "AggK");
    SEXP Ak_s = cpt_list_elt(cl, "A_k");
    SEXP Ad_s = cpt_list_elt(cl, "A_d");
    SEXP Av_s = cpt_list_elt(cl, "A_v");

    if (!Rf_isMatrix(zd_s) || TYPEOF(zd_s) != REALSXP ||
        TYPEOF(mult_s) != REALSXP ||
        !Rf_isMatrix(Rmat_s) || TYPEOF(Rmat_s) != REALSXP ||
        !Rf_isMatrix(Lmat_s) || TYPEOF(Lmat_s) != REALSXP ||
        !Rf_isMatrix(AggX_s) || TYPEOF(AggX_s) != REALSXP ||
        !Rf_isMatrix(AggK_s) || TYPEOF(AggK_s) != REALSXP ||
        TYPEOF(Ak_s) != INTSXP || TYPEOF(Ad_s) != INTSXP || TYPEOF(Av_s) != REALSXP) {
      vmaxset(vmax);
      UNPROTECT(1);
      return R_NilValue;
    }

    int Md = Rf_nrows(zd_s);
    int p = Rf_ncols(zd_s);
    int R = Rf_ncols(Rmat_s);
    int Dl = Rf_ncols(AggX_s);
    double M_dbl = Rf_asReal(M_s);
    R_xlen_t ne = XLENGTH(Av_s);

    if (Md < 1 || Rf_nrows(AggX_s) != C || Rf_nrows(AggK_s) != C || Rf_ncols(AggK_s) != R ||
        Rf_nrows(Rmat_s) != Md || Rf_nrows(Lmat_s) != Dl || Rf_ncols(Lmat_s) != Md ||
        XLENGTH(mult_s) != Md || XLENGTH(Ak_s) != ne || XLENGTH(Ad_s) != ne ||
        !R_FINITE(M_dbl) || M_dbl <= 0.0) {
      vmaxset(vmax);
      UNPROTECT(1);
      return R_NilValue;
    }

    const double *AggX = REAL(AggX_s);
    const double *AggK = REAL(AggK_s);
    const double *Rmat = REAL(Rmat_s);
    const double *Lmat = REAL(Lmat_s);
    const double *mult = REAL(mult_s);
    const int *Ak = INTEGER(Ak_s);
    const int *Ad = INTEGER(Ad_s);
    const double *Av = REAL(Av_s);

    /* 1. EX = E %*% AggX (b x Dl col-major) and EK = E %*% AggK (b x R col-major) */
    double *EX = (Dl > 0) ? (double *) R_alloc((size_t) b * (size_t) Dl, sizeof(double)) : NULL;
    if (Dl > 0) {
      memset(EX, 0, (size_t) b * (size_t) Dl * sizeof(double));
      for (int c = 0; c < Dl; ++c) {
        const double *ax_col = AggX + (size_t) c * (size_t) C;
        double *ex_col = EX + (size_t) c * (size_t) b;
        for (int k = 0; k < C; ++k) {
          double ax = ax_col[k];
          if (ax != 0.0) {
            const double *ek = E + (size_t) k * (size_t) b;
            for (int r = 0; r < b; ++r) ex_col[r] += ek[r] * ax;
          }
        }
      }
    }

    double *EK = (R > 0) ? (double *) R_alloc((size_t) b * (size_t) R, sizeof(double)) : NULL;
    if (R > 0) {
      memset(EK, 0, (size_t) b * (size_t) R * sizeof(double));
      for (int c = 0; c < R; ++c) {
        const double *ak_col = AggK + (size_t) c * (size_t) C;
        double *ek_col = EK + (size_t) c * (size_t) b;
        for (int k = 0; k < C; ++k) {
          double ak = ak_col[k];
          if (ak != 0.0) {
            const double *ek = E + (size_t) k * (size_t) b;
            for (int r = 0; r < b; ++r) ek_col[r] += ek[r] * ak;
          }
        }
      }
    }

    /* 2. Accumulate Q_act (n_qd x b, contiguous b-vector per active grid point d) */
    int *d_map = (int *) R_alloc((size_t) Md, sizeof(int));
    for (int d = 0; d < Md; ++d) d_map[d] = -1;
    int n_cand = 0;
    for (R_xlen_t e = 0; e < ne; ++e) {
      int d0 = Ad[e] - 1;
      int k0 = Ak[e] - 1;
      if (d0 < 0 || d0 >= Md || k0 < 0 || k0 >= C) {
        vmaxset(vmax);
        UNPROTECT(1);
        return R_NilValue;
      }
      if (Av[e] != 0.0 && d_map[d0] < 0) {
        d_map[d0] = 1;
        ++n_cand;
      }
    }

    int *cand_idx = (n_cand > 0) ? (int *) R_alloc((size_t) n_cand, sizeof(int)) : NULL;
    int pos = 0;
    for (int d = 0; d < Md; ++d) {
      if (d_map[d] > 0) {
        d_map[d] = pos;
        cand_idx[pos++] = d;
      }
    }

    int n_qd = 0;
    int *qd_idx = NULL;
    double *Q_qd = NULL;
    if (n_cand > 0) {
      double *Q_cand = (double *) R_alloc((size_t) n_cand * (size_t) b, sizeof(double));
      memset(Q_cand, 0, (size_t) n_cand * (size_t) b * sizeof(double));
      for (R_xlen_t e = 0; e < ne; ++e) {
        double v = Av[e];
        if (v == 0.0) continue;
        int a = d_map[Ad[e] - 1];
        const double *ek = E + (size_t) (Ak[e] - 1) * (size_t) b;
        double *qa = Q_cand + (size_t) a * (size_t) b;
        for (int r = 0; r < b; ++r) qa[r] += ek[r] * v;
      }

      qd_idx = (int *) R_alloc((size_t) n_cand, sizeof(int));
      for (int a = 0; a < n_cand; ++a) {
        const double *qa = Q_cand + (size_t) a * (size_t) b;
        int nz = 0;
        for (int r = 0; r < b; ++r) {
          if (qa[r] != 0.0) {
            nz = 1;
            break;
          }
        }
        if (nz) {
          qd_idx[n_qd] = cand_idx[a];
          if (n_qd != a) {
            memcpy(Q_cand + (size_t) n_qd * (size_t) b, qa, (size_t) b * sizeof(double));
          }
          ++n_qd;
        }
      }
      Q_qd = Q_cand;
    }

    /* 3. Orthant contraction T1 (Md x b, contiguous b-vector per target m) */
    double *T1 = (double *) R_alloc((size_t) Md * (size_t) b, sizeof(double));
    if (n_qd > 0) {
      cpt_orthant_sum(REAL(zd_s), Md, p, qd_idx, n_qd, Q_qd, b, T1);
    } else {
      memset(T1, 0, (size_t) Md * (size_t) b * sizeof(double));
    }

    /* 4. Combine Js and accumulate (N / M) * sum_m mult[m] * Js[r, m]^2 */
    double scale_NM = N_dbl / M_dbl;
    double *js_buf = (double *) R_alloc((size_t) b, sizeof(double));
    for (int m = 0; m < Md; ++m) {
      double wm = mult[m];
      memcpy(js_buf, T1 + (size_t) m * (size_t) b, (size_t) b * sizeof(double));
      for (int c = 0; c < R; ++c) {
        double rm = Rmat[(size_t) m + (size_t) c * (size_t) Md];
        if (rm != 0.0) {
          const double *ek_col = EK + (size_t) c * (size_t) b;
          for (int r = 0; r < b; ++r) js_buf[r] -= ek_col[r] * rm;
        }
      }
      if (Dl > 0) {
        const double *lm_col = Lmat + (size_t) m * (size_t) Dl;
        for (int c = 0; c < Dl; ++c) {
          double lm = lm_col[c];
          if (lm != 0.0) {
            const double *ex_col = EX + (size_t) c * (size_t) b;
            for (int r = 0; r < b; ++r) js_buf[r] += ex_col[r] * lm;
          }
        }
      }
      if (wm != 0.0) {
        double fac = scale_NM * wm;
        for (int r = 0; r < b; ++r) {
          double jv = js_buf[r] * invN;
          out[r] += fac * (jv * jv);
        }
      }
    }

    vmaxset(vmax);
  }

  UNPROTECT(1);
  return out_sexp;
}

/* ---------------------------------------------------------------------------
   4. Early-Exit Cluster-Process Degeneracy Check: cpt_process_nonzero
   --------------------------------------------------------------------------- */
SEXP C_cdid_cpt_process_nonzero(SEXP cells_sexp, SEXP C_sexp) {
  if (TYPEOF(cells_sexp) != VECSXP) return R_NilValue;
  int C = Rf_asInteger(C_sexp);
  if (C < 1) return R_NilValue;

  int n_cells = (int) XLENGTH(cells_sexp);
  const void *vmax = vmaxget();
  double *q_col = (double *) R_alloc((size_t) C, sizeof(double));

  for (int ci = 0; ci < n_cells; ++ci) {
    SEXP cl = VECTOR_ELT(cells_sexp, ci);
    SEXP zd_s = cpt_list_elt(cl, "zd");
    SEXP Rmat_s = cpt_list_elt(cl, "Rmat");
    SEXP Lmat_s = cpt_list_elt(cl, "Lmat");
    SEXP AggX_s = cpt_list_elt(cl, "AggX");
    SEXP AggK_s = cpt_list_elt(cl, "AggK");
    SEXP Ak_s = cpt_list_elt(cl, "A_k");
    SEXP Ad_s = cpt_list_elt(cl, "A_d");
    SEXP Av_s = cpt_list_elt(cl, "A_v");

    if (!Rf_isMatrix(zd_s) || TYPEOF(zd_s) != REALSXP ||
        !Rf_isMatrix(Rmat_s) || TYPEOF(Rmat_s) != REALSXP ||
        !Rf_isMatrix(Lmat_s) || TYPEOF(Lmat_s) != REALSXP ||
        !Rf_isMatrix(AggX_s) || TYPEOF(AggX_s) != REALSXP ||
        !Rf_isMatrix(AggK_s) || TYPEOF(AggK_s) != REALSXP ||
        TYPEOF(Ak_s) != INTSXP || TYPEOF(Ad_s) != INTSXP || TYPEOF(Av_s) != REALSXP) {
      vmaxset(vmax);
      return R_NilValue;
    }

    int Md = Rf_nrows(zd_s);
    int p = Rf_ncols(zd_s);
    int R = Rf_ncols(Rmat_s);
    int Dl = Rf_ncols(AggX_s);
    R_xlen_t ne = XLENGTH(Av_s);

    if (Md < 1 || Rf_nrows(AggX_s) != C || Rf_nrows(AggK_s) != C || Rf_ncols(AggK_s) != R ||
        Rf_nrows(Rmat_s) != Md || Rf_nrows(Lmat_s) != Dl || Rf_ncols(Lmat_s) != Md ||
        XLENGTH(Ak_s) != ne || XLENGTH(Ad_s) != ne) {
      vmaxset(vmax);
      return R_NilValue;
    }

    const double *zd = REAL(zd_s);
    const double *Rmat = REAL(Rmat_s);
    const double *Lmat = REAL(Lmat_s);
    const double *AggX = REAL(AggX_s);
    const double *AggK = REAL(AggK_s);
    const int *Ak = INTEGER(Ak_s);
    const int *Ad = INTEGER(Ad_s);
    const double *Av = REAL(Av_s);

    for (int m = 0; m < Md; ++m) {
      memset(q_col, 0, (size_t) C * sizeof(double));
      for (R_xlen_t e = 0; e < ne; ++e) {
        double v = Av[e];
        if (v == 0.0) continue;
        int d0 = Ad[e] - 1;
        int k0 = Ak[e] - 1;
        if (d0 < 0 || d0 >= Md || k0 < 0 || k0 >= C) {
          vmaxset(vmax);
          return R_NilValue;
        }
        int le = 1;
        for (int j = 0; j < p; ++j) {
          if (zd[(size_t) d0 + (size_t) j * (size_t) Md] >
              zd[(size_t) m + (size_t) j * (size_t) Md]) {
            le = 0;
            break;
          }
        }
        if (le) q_col[k0] += v;
      }

      const double *lm_col = Lmat + (size_t) m * (size_t) Dl;
      for (int k = 0; k < C; ++k) {
        double val = q_col[k];
        for (int r = 0; r < R; ++r) {
          val -= AggK[(size_t) k + (size_t) r * (size_t) C] *
                 Rmat[(size_t) m + (size_t) r * (size_t) Md];
        }
        for (int c = 0; c < Dl; ++c) {
          val += AggX[(size_t) k + (size_t) c * (size_t) C] * lm_col[c];
        }
        if (val != 0.0) {
          vmaxset(vmax);
          return Rf_ScalarLogical(1);
        }
      }
    }
  }

  vmaxset(vmax);
  return Rf_ScalarLogical(0);
}
