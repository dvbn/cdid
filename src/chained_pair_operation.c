/* Lazy, operation-owned sparse posting POSITIONS only. Outcome vector reads,
   endpoint validation and right-minus-left remain in the original R helper. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>

#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)) || \
    (defined(__ASSOC_MATH__) && (__ASSOC_MATH__ > 0)) || \
    (defined(__RECIPROCAL_MATH__) && (__RECIPROCAL_MATH__ > 0))
#error "Unsafe floating-point compilation is forbidden for checked sparse pair operations."
#endif

/* Narrow bridges to the existing arithmetic, capacity and exact search rules. */
extern void cdid_pair_operation_arithmetic(void);
extern double cdid_pair_operation_cap(SEXP);
extern void cdid_pair_operation_error(const char *, const char *);
extern int cdid_pair_operation_search(SEXP, int, int, int);

typedef struct { int key, slot; } pair49_row;
typedef struct {
  int B, T;
  double reservation;
  pair49_row *rows;
  int *positions;
  uint64_t searches, hits, queries;
} pair49_state;

static int pair49_plain(SEXP x, int type) {
  return TYPEOF(x) == type && !OBJECT(x) && !ALTREP(x) && ATTRIB(x) == R_NilValue;
}
static int pair49_row_compare(const void *x, const void *y) {
  int a = ((const pair49_row *)x)->key, b = ((const pair49_row *)y)->key;
  return a < b ? -1 : (a > b ? 1 : 0);
}
static void pair49_release(SEXP pointer) {
  pair49_state *state = (pair49_state *)R_ExternalPtrAddr(pointer);
  if (state == NULL) return;
  if (state->rows != NULL) R_Free(state->rows);
  if (state->positions != NULL) R_Free(state->positions);
  R_Free(state);
  R_ClearExternalPtr(pointer);
}
static pair49_state *pair49_get(SEXP pointer) {
  if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrTag(pointer) != Rf_install("cdid_pair_position_operation") ||
      TYPEOF(R_ExternalPtrProtected(pointer)) != VECSXP || XLENGTH(R_ExternalPtrProtected(pointer)) != 8 ||
      R_ExternalPtrAddr(pointer) == NULL)
    cdid_pair_operation_error("Sparse pair operation lacks its live owned storage.", "cdid_gmm_bound_error");
  return (pair49_state *)R_ExternalPtrAddr(pointer);
}
static int pair49_slot(pair49_state *state, int row) {
  int lo = 0, hi = state->B - 1;
  while (lo <= hi) {
    int mid = lo + (hi - lo) / 2;
    if (state->rows[mid].key == row) return state->rows[mid].slot;
    if (state->rows[mid].key < row) lo = mid + 1; else hi = mid - 1;
  }
  return -1;
}
static int pair49_index(SEXP rows, R_xlen_t i, int N) {
  if (TYPEOF(rows) == INTSXP) {
    int row = INTEGER(rows)[i];
    return row == NA_INTEGER || row < 1 || row > N ? 0 : row;
  }
  double row = REAL(rows)[i];
  if (!R_FINITE(row) || row < 1 || row > N || row != (double)(int)row) return 0;
  return (int)row;
}

SEXP C_gmm_pair_operation_create(SEXP rows, SEXP a, SEXP pa, SEXP po,
    SEXP starts, SEXP ends, SEXP counts, SEXP y, SEXP reservation) {
  /* Unsupported sources decline before touching their elements or allocating
     owned storage. Rows are the wrapper's already charged plain B-size copy. */
  if (!pair49_plain(rows, INTSXP) || !pair49_plain(a, INTSXP) ||
      !pair49_plain(pa, INTSXP) || !pair49_plain(po, INTSXP) ||
      !pair49_plain(starts, INTSXP) || !pair49_plain(ends, INTSXP) || !pair49_plain(counts, INTSXP) ||
      !(pair49_plain(y, REALSXP) || pair49_plain(y, INTSXP))) return R_NilValue;
  R_xlen_t B = XLENGTH(rows), T = XLENGTH(starts), N = XLENGTH(a), np = XLENGTH(pa);
  if (B < 1 || B > INT_MAX || T < 1 || T > INT_MAX || N < 1 || N > INT_MAX ||
      np < 1 || np > INT_MAX || XLENGTH(po) != np || XLENGTH(ends) != T || XLENGTH(counts) != T)
    return R_NilValue;
  cdid_pair_operation_arithmetic();
  double cap = cdid_pair_operation_cap(reservation);
  uint64_t bt = (uint64_t)B * (uint64_t)T;
  uint64_t other = 4096ULL + 64ULL * (uint64_t)B + 8ULL * (uint64_t)T;
  if (bt > (UINT64_MAX - other) / 3ULL) return R_NilValue;
  uint64_t required = other + 3ULL * bt;
  if (!R_FINITE(cap) || cap < (double)required || bt > SIZE_MAX / sizeof(int) ||
      (uint64_t)B > SIZE_MAX / sizeof(pair49_row)) return R_NilValue;
  for (R_xlen_t i = 0; i < B; ++i)
    if (INTEGER(rows)[i] == NA_INTEGER || INTEGER(rows)[i] < 1 || INTEGER(rows)[i] > N) return R_NilValue;

  SEXP sources = PROTECT(Rf_allocVector(VECSXP, 8));
  SEXP inputs[8] = {rows, a, pa, po, starts, ends, counts, y};
  for (int j = 0; j < 8; ++j) SET_VECTOR_ELT(sources, j, inputs[j]);
  SEXP pointer = PROTECT(R_MakeExternalPtr(NULL, Rf_install("cdid_pair_position_operation"), sources));
  pair49_state *state = (pair49_state *)R_Calloc(1, pair49_state);
  R_SetExternalPtrAddr(pointer, state);
  R_RegisterCFinalizerEx(pointer, pair49_release, TRUE);
  state->B = (int)B; state->T = (int)T; state->reservation = cap;
  state->rows = (pair49_row *)R_Calloc((size_t)B, pair49_row);
  state->positions = (int *)R_Calloc((size_t)bt, int);
  for (R_xlen_t i = 0; i < B; ++i) { state->rows[i].key = INTEGER(rows)[i]; state->rows[i].slot = (int)i; }
  qsort(state->rows, (size_t)B, sizeof(pair49_row), pair49_row_compare);
  for (R_xlen_t i = 1; i < B; ++i) if (state->rows[i].key == state->rows[i - 1].key) {
    pair49_release(pointer); UNPROTECT(2); return R_NilValue;
  }
  for (uint64_t i = 0; i < bt; ++i) state->positions[i] = -1;
  UNPROTECT(2); return pointer;
}

static int pair49_position(pair49_state *state, SEXP pa, int slot, int date,
    int a, int lo, int hi) {
  R_xlen_t at = (R_xlen_t)date * state->B + slot;
  int position = state->positions[at];
  if (position >= 0) { ++state->hits; return position; }
  ++state->searches;
  position = cdid_pair_operation_search(pa, lo, hi, a);
  state->positions[at] = position;
  return position;
}

SEXP C_gmm_pair_operation_query(SEXP pointer, SEXP rows, SEXP dates,
    SEXP read, SEXP reservation, SEXP action) {
  if (!pair49_plain(action, INTSXP) || XLENGTH(action) != 1 || INTEGER(action)[0] < 0 || INTEGER(action)[0] > 2)
    cdid_pair_operation_error("Invalid sparse pair operation action.", "cdid_gmm_bound_error");
  int what = INTEGER(action)[0];
  if (what == 2) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrTag(pointer) != Rf_install("cdid_pair_position_operation"))
      cdid_pair_operation_error("Cannot release foreign sparse pair storage.", "cdid_gmm_bound_error");
    pair49_release(pointer); return R_NilValue;
  }
  pair49_state *state = pair49_get(pointer);
  if (what == 1) {
    cdid_pair_operation_arithmetic();
    SEXP out = PROTECT(Rf_allocVector(REALSXP, 3));
    REAL(out)[0] = (double)state->searches; REAL(out)[1] = (double)state->hits; REAL(out)[2] = (double)state->queries;
    UNPROTECT(1); return out;
  }
  /* A decline has no cache update and preserves the original R/scalar path. */
  if (!(pair49_plain(rows, INTSXP) || pair49_plain(rows, REALSXP)) ||
      !pair49_plain(dates, INTSXP) || XLENGTH(dates) != 2 ||
      !pair49_plain(read, LGLSXP) || XLENGTH(read) != 1 || LOGICAL(read)[0] == NA_LOGICAL ||
      XLENGTH(rows) > state->B) return R_NilValue;
  int ds = INTEGER(dates)[0], dt = INTEGER(dates)[1];
  if (ds == NA_INTEGER || dt == NA_INTEGER || ds < 1 || dt < 1 || ds > state->T || dt > state->T) return R_NilValue;
  SEXP source = R_ExternalPtrProtected(pointer);
  SEXP a = VECTOR_ELT(source, 1), pa = VECTOR_ELT(source, 2), po = VECTOR_ELT(source, 3);
  const int *starts = INTEGER(VECTOR_ELT(source, 4)), *ends = INTEGER(VECTOR_ELT(source, 5));
  const int *counts = INTEGER(VECTOR_ELT(source, 6));
  R_xlen_t n = XLENGTH(rows); int N = (int)XLENGTH(a), np = (int)XLENGTH(pa);
  int positive = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    int row = pair49_index(rows, i, N);
    if (!row || pair49_slot(state, row) < 0) return R_NilValue;
    int ai = INTEGER(a)[row - 1];
    /* native_ok in the original R wrapper declines these rows; its scalar
       keep=which(ai>0) path, not the native malformed guard, decides them. */
    if (ai == NA_INTEGER || ai < 0) return R_NilValue;
    if (ai > 0) positive = 1;
  }
  if (counts[ds - 1] == NA_INTEGER || counts[dt - 1] == NA_INTEGER) return R_NilValue;
  int lookup = positive && counts[ds - 1] > 0 && counts[dt - 1] > 0;
  int bounds[4] = {starts[ds - 1], ends[ds - 1], starts[dt - 1], ends[dt - 1]};
  if (lookup) for (int j = 0; j < 4; ++j)
    if (bounds[j] == NA_INTEGER || bounds[j] < 1 || bounds[j] > np) return R_NilValue;
  if (lookup && (bounds[0] > bounds[1] || bounds[2] > bounds[3])) return R_NilValue;
  cdid_pair_operation_arithmetic();
  double cap = cdid_pair_operation_cap(reservation);
  uint64_t required = 4096ULL + 64ULL * (uint64_t)state->B + 8ULL * (uint64_t)state->T +
    3ULL * (uint64_t)state->B * (uint64_t)state->T;
  if (cap != state->reservation || cap < (double)required)
    cdid_pair_operation_error("Sparse pair operation live storage exceeds its admitted reservation.", "cdid_gmm_bound_error");
  size_t size = (size_t)(n > 0 ? n : 1);
  int *obs = (int *)R_alloc(size, sizeof(int));
  int *left = (int *)R_alloc(size, sizeof(int));
  int *right = (int *)R_alloc(size, sizeof(int));
  R_xlen_t k = 0;
  ++state->queries;
  if (lookup) for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && i % 4096 == 0) R_CheckUserInterrupt();
    int row = pair49_index(rows, i, N), ai = INTEGER(a)[row - 1];
    if (ai <= 0) continue;
    int slot = pair49_slot(state, row);
    int ls = pair49_position(state, pa, slot, ds - 1, ai, bounds[0], bounds[1]);
    if (!ls) continue;
    int rt = pair49_position(state, pa, slot, dt - 1, ai, bounds[2], bounds[3]);
    if (!rt) continue;
    obs[k] = (int)i + 1; left[k] = ls; right[k] = rt; ++k;
  }
  int need = LOGICAL(read)[0];
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SEXP oo = PROTECT(Rf_allocVector(INTSXP, k));
  SEXP ss = PROTECT(need ? Rf_allocVector(INTSXP, k) : R_NilValue);
  SEXP tt = PROTECT(need ? Rf_allocVector(INTSXP, k) : R_NilValue);
  for (R_xlen_t i = 0; i < k; ++i) {
    INTEGER(oo)[i] = obs[i];
    if (need) { INTEGER(ss)[i] = INTEGER(po)[left[i] - 1]; INTEGER(tt)[i] = INTEGER(po)[right[i] - 1]; }
  }
  SET_VECTOR_ELT(out, 0, oo); SET_VECTOR_ELT(out, 1, ss); SET_VECTOR_ELT(out, 2, tt);
  SET_STRING_ELT(names, 0, Rf_mkChar("obs")); SET_STRING_ELT(names, 1, Rf_mkChar("pos_s")); SET_STRING_ELT(names, 2, Rf_mkChar("pos_t"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(5); return out;
}
