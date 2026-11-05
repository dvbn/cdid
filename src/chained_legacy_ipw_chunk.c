/* Ordinary legacy CP04 full-coordinate replay. R retains its logit, clipping,
   matrix products and two-triangular-solve primitives. No eager pair reads. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)) || \
    (defined(__ASSOC_MATH__) && (__ASSOC_MATH__ > 0)) || \
    (defined(__RECIPROCAL_MATH__) && (__RECIPROCAL_MATH__ > 0))
#error "Unsafe floating-point compilation is forbidden for CP04 replay."
#endif

extern void cdid_pair_operation_arithmetic(void);
static SEXP item(SEXP x, const char *name) {
  if (TYPEOF(x) != VECSXP) return R_NilValue;
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (TYPEOF(names) != STRSXP) return R_NilValue;
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(x, i);
  return R_NilValue;
}
static int plain(SEXP x, int type) {
  return TYPEOF(x) == type && !OBJECT(x) && !ALTREP(x) && ATTRIB(x) == R_NilValue;
}
static int list_ok(SEXP x) {
  if (TYPEOF(x) != VECSXP || OBJECT(x) || ALTREP(x)) return 0;
  return ATTRIB(x) == R_NilValue || (TAG(ATTRIB(x)) == R_NamesSymbol && CDR(ATTRIB(x)) == R_NilValue);
}
static int matrix_ok(SEXP x, int nr, int nc) {
  if (TYPEOF(x) != REALSXP || OBJECT(x) || ALTREP(x)) return 0;
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 || INTEGER(dim)[0] != nr || INTEGER(dim)[1] != nc) return 0;
  for (SEXP a = ATTRIB(x); a != R_NilValue; a = CDR(a))
    if (TAG(a) != R_DimSymbol && TAG(a) != R_DimNamesSymbol) return 0;
  return XLENGTH(x) == (R_xlen_t)nr * nc;
}
/* Ordinary full-F selection designs retain model.matrix metadata in the
   immutable public/standalone readers. These attributes do not participate
   in the numeric subsetting/product; no metadata value is read or copied.
   Fitted H/cross/partial matrices keep the stricter matrix_ok above. */
static int selection_design_header_ok(SEXP x, int nr, int nc) {
  if (TYPEOF(x) != REALSXP || OBJECT(x) || ALTREP(x)) return 0;
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
      INTEGER(dim)[0] != nr || INTEGER(dim)[1] != nc) return 0;
  for (SEXP a = ATTRIB(x); a != R_NilValue; a = CDR(a)) {
    SEXP tag = TAG(a);
    if (tag != R_DimSymbol && tag != R_DimNamesSymbol &&
        tag != Rf_install("intercept") && tag != Rf_install("has_intercept") &&
        tag != Rf_install("assign") && tag != Rf_install("contrasts")) return 0;
  }
  return XLENGTH(x) == (R_xlen_t)nr * nc;
}
/* The canonical frame design may be a base REAL ALTREP wrapper. Inspect
   only its ordinary matrix header here: never Length/REAL on that borrowed
   wrapper. The original R row-subset primitive reads it at its old position. */
static int frame_matrix_header_ok(SEXP x, int nr, int nc) {
  if (TYPEOF(x) != REALSXP || OBJECT(x) || nr < 1 || nc < 1) return 0;
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || OBJECT(dim) || ALTREP(dim) || XLENGTH(dim) != 2 ||
      INTEGER(dim)[0] != nr || INTEGER(dim)[1] != nc) return 0;
  for (SEXP a = ATTRIB(x); a != R_NilValue; a = CDR(a))
    if (TAG(a) != R_DimSymbol && TAG(a) != R_DimNamesSymbol) return 0;
  return 1;
}
static int nrow(SEXP x) {
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  return TYPEOF(d) == INTSXP && XLENGTH(d) == 2 ? INTEGER(d)[0] : -1;
}
static int ncol(SEXP x) {
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  return TYPEOF(d) == INTSXP && XLENGTH(d) == 2 ? INTEGER(d)[1] : -1;
}
/* Force a binary64 store at each original vector arithmetic operation. */
static double mul(double a, double b) { volatile double z = a * b; return z; }
static double add(double a, double b) { volatile double z = a + b; return z; }
static double sub(double a, double b) { volatile double z = a - b; return z; }
static double quo(double a, double b) { volatile double z = a / b; return z; }
static SEXP zero_matrix(int nr, int nc) {
  SEXP x = Rf_allocMatrix(REALSXP, nr, nc);
  memset(REAL(x), 0, (size_t)nr * nc * sizeof(double));
  return x;
}
static SEXP zero_vector(int n) {
  SEXP x = Rf_allocVector(REALSXP, n);
  memset(REAL(x), 0, (size_t)n * sizeof(double));
  return x;
}
static SEXP eval1(SEXP f, SEXP x) {
  SEXP call = PROTECT(Rf_lang2(f, x));
  SEXP out = Rf_eval(call, R_BaseEnv);
  UNPROTECT(1); return out;
}
static SEXP eval2(SEXP f, SEXP x, SEXP y) {
  SEXP call = PROTECT(Rf_lang3(f, x, y));
  SEXP out = Rf_eval(call, R_BaseEnv);
  UNPROTECT(1); return out;
}
static SEXP subset_rows(SEXP rows, const int *rr, int n) {
  SEXP out = Rf_allocVector(TYPEOF(rows), n);
  for (int i = 0; i < n; ++i) {
    if (TYPEOF(rows) == INTSXP) INTEGER(out)[i] = INTEGER(rows)[rr[i]];
    else REAL(out)[i] = REAL(rows)[rr[i]];
  }
  return out;
}
static int row_at(SEXP rows, int i) {
  return (TYPEOF(rows) == INTSXP ? INTEGER(rows)[i] : (int)REAL(rows)[i]) - 1;
}
static SEXP pair_call(SEXP pair, SEXP rows, SEXP descriptor, int read) {
  SEXP flag = PROTECT(Rf_ScalarLogical(read));
  SEXP call = PROTECT(Rf_lang4(pair, rows, descriptor, flag));
  SEXP out = Rf_eval(call, R_BaseEnv);
  UNPROTECT(2); return out;
}
static void finite_call(SEXP finite, SEXP x, const char *label) {
  SEXP text = PROTECT(Rf_mkString(label));
  (void)eval2(finite, x, text);
  UNPROTECT(1);
}

/* All structural declines precede the first arithmetic or helper call. They
   inspect headers only, never an unvisited model's numerical values. */
SEXP C_gmm_legacy_ipw_chunk(SEXP reader, SEXP frame, SEXP backend, SEXP rows,
                          SEXP tids, SEXP sids, SEXP pair, SEXP helpers, SEXP funding) {
  SEXP st = item(reader, "state"), Z = item(st, "Z"), Xall = item(frame, "Z_a");
  SEXP ci = item(st, "cohort_index"), ai = item(st, "a_index"), w = item(st, "w");
  SEXP U = item(st, "U"), rho = item(st, "rho"), a = item(st, "a"), mass = item(st, "mass_U");
  SEXP ss = item(backend, "selection_specs"), sf = item(backend, "selections");
  SEXP ts = item(backend, "treatment_specs"), tf = item(backend, "treatments");
  int B = (int)XLENGTH(rows), N = (int)XLENGTH(w), K = (int)XLENGTH(a);
  int d = ncol(Z), p = ncol(Xall), nA = nrow(Xall);
  if (B < 1 || XLENGTH(rows) > INT_MAX || N < B || K < 1 || d < 1 || p < 1 ||
      !selection_design_header_ok(Z, N, d) || !frame_matrix_header_ok(Xall, nA, p) ||
      !plain(ci, INTSXP) || XLENGTH(ci) != N || !plain(ai, INTSXP) || XLENGTH(ai) != N ||
      !plain(w, REALSXP) || !plain(U, LGLSXP) || XLENGTH(U) != N ||
      !plain(rho, REALSXP) || XLENGTH(rho) != N || !plain(a, REALSXP) ||
      !plain(mass, REALSXP) || XLENGTH(mass) != K || !plain(tids, INTSXP) || !plain(sids, INTSXP) ||
      !list_ok(st) || !list_ok(ss) || !list_ok(sf) || !list_ok(ts) || !list_ok(tf) ||
      !Rf_isFunction(pair)) return R_NilValue;
  if ((TYPEOF(rows) != INTSXP && TYPEOF(rows) != REALSXP) || OBJECT(rows) || ATTRIB(rows) != R_NilValue) return R_NilValue;
  for (int i = 0; i < B; ++i) {
    double value = TYPEOF(rows) == INTSXP ? INTEGER(rows)[i] : REAL(rows)[i];
    if (!R_FINITE(value) || value != floor(value) || value < 1 || value > N) return R_NilValue;
    int r = (int)value - 1;
    if (INTEGER(ci)[r] < 1 || INTEGER(ci)[r] > K || INTEGER(ai)[r] < 0 || INTEGER(ai)[r] > nA ||
        LOGICAL(U)[r] == NA_LOGICAL) return R_NilValue;
  }
  for (R_xlen_t j = 0; j < XLENGTH(sids); ++j) {
    int s = INTEGER(sids)[j] - 1;
    if (s < 0 || s >= XLENGTH(sf) || s >= XLENGTH(ss)) return R_NilValue;
    SEXP spec = VECTOR_ELT(ss, s), fit = VECTOR_ELT(sf, s), cohort = item(spec, "cohort");
    if (!list_ok(spec) || !list_ok(fit) || !plain(cohort, INTSXP) || XLENGTH(cohort) != 1 || INTEGER(cohort)[0] < 1 || INTEGER(cohort)[0] > K ||
        !plain(item(fit, "coefficients"), REALSXP) || XLENGTH(item(fit, "coefficients")) != d ||
        !matrix_ok(item(fit, "H_chol"), d, d)) return R_NilValue;
  }
  double treatment_peak = 0;
  for (R_xlen_t j = 0; j < XLENGTH(tids); ++j) {
    int t = INTEGER(tids)[j] - 1;
    if (t < 0 || t >= XLENGTH(tf) || t >= XLENGTH(ts)) return R_NilValue;
    SEXP spec = VECTOR_ELT(ts, t), fit = VECTOR_ELT(tf, t);
    SEXP cohorts = item(spec, "cohorts"), ids = item(spec, "selection_ids"), part = item(fit, "partials");
    int m = (int)XLENGTH(cohorts);
    if (d == INT_MAX || m > INT_MAX / (d + 1)) return R_NilValue;
    int u = m * (d + 1);
    SEXP means = item(fit, "means");
    double peak = (double)B*(u + 10.0*p + 64) + (double)p*u + 4.0*u + 4.0*p +
      512.0*(p + d + 16) + 8.0*(p + d + 1)*(p + d + 1);
    if (peak > treatment_peak) treatment_peak = peak;
    if (!list_ok(spec) || !list_ok(fit) || !list_ok(part) || !plain(cohorts, INTSXP) || !plain(ids, INTSXP) || XLENGTH(ids) != m || m < 1 ||
        !plain(item(fit, "coefficients"), REALSXP) || XLENGTH(item(fit, "coefficients")) != p ||
        !matrix_ok(item(fit, "H_chol"), p, p) || !matrix_ok(item(fit, "cross"), p, u) ||
        !matrix_ok(item(part, "upstream"), 4, u) || !matrix_ok(item(part, "beta"), 4, p) ||
        TYPEOF(means) != REALSXP || OBJECT(means) || ALTREP(means) || XLENGTH(means) != 4 ||
        ATTRIB(means) == R_NilValue || TAG(ATTRIB(means)) != R_NamesSymbol || CDR(ATTRIB(means)) != R_NilValue ||
        (TYPEOF(Rf_getAttrib(means, R_NamesSymbol)) != STRSXP ||
          XLENGTH(Rf_getAttrib(means, R_NamesSymbol)) != 4)) return R_NilValue;
    SEXP names = Rf_getAttrib(means, R_NamesSymbol);
    const char *expected[] = {"T", "M", "dp", "dm"};
    for (int k = 0; k < 4; ++k) if (strcmp(CHAR(STRING_ELT(names, k)), expected[k]) != 0) return R_NilValue;
    for (int k = 0; k < m; ++k) {
      int s = INTEGER(ids)[k];
      if (INTEGER(cohorts)[k] < 1 || INTEGER(cohorts)[k] > K || s < 1 || s > XLENGTH(sf)) return R_NilValue;
      int found = 0;
      for (R_xlen_t q = 0; q < XLENGTH(sids); ++q) if (INTEGER(sids)[q] == s) found = 1;
      if (!found) return R_NilValue;
    }
  }
  /* Reserve all persistent model results plus the current and previous
     primitive workspaces, headers/calls and replacement overlap. The pair
     operation has its separate eightfold reservation; this fits only the
     unchanged 24-fold remainder, fixed-dimension scratch and metadata. */
  double persistent = (double)B*(K + (double)XLENGTH(sids)*(d + 2) + XLENGTH(tids)) +
    32.0*(3.0*XLENGTH(sf) + XLENGTH(tf) + K + XLENGTH(sids) + XLENGTH(tids) + 32);
  double selection_peak = (double)B*(6.0*d + 32) + 2.0*d*d + 512.0*(d + 8);
  double peak = persistent + (selection_peak > treatment_peak ? selection_peak : treatment_peak);
  if (TYPEOF(funding) != REALSXP || XLENGTH(funding) != 1 || !R_FINITE(REAL(funding)[0]) ||
      !R_FINITE(peak) || peak > REAL(funding)[0]) return R_NilValue;
  /* Checked outer evaluator already validates rows; source/header admission
     above makes no numerical substitution or retry possible beyond here. */
  cdid_pair_operation_arithmetic();
  int *rr = (int *)R_alloc((size_t)B, sizeof(int));
  SEXP predict = item(helpers, "cp04_predict_selection"), solve = item(helpers, "cp04_solve");
  SEXP logit = item(helpers, "ipw_logit_values"), clip_floor = item(helpers, "ipw_clip_floor");
  SEXP finite = item(helpers, "cp04_finite"), transpose = Rf_findFun(Rf_install("t"), R_BaseEnv);
  SEXP product = Rf_findFun(Rf_install("%*%"), R_BaseEnv);
  SEXP models = PROTECT(Rf_allocVector(VECSXP, XLENGTH(tf)));
  SEXP alpha = PROTECT(Rf_allocVector(VECSXP, XLENGTH(sf)));
  SEXP probs = PROTECT(Rf_allocVector(VECSXP, XLENGTH(sf)));
  SEXP responses = PROTECT(Rf_allocVector(VECSXP, XLENGTH(sf)));
  SEXP aif = PROTECT(zero_matrix(B, K));
  for (int i = 0; i < B; ++i) {
    int r = row_at(rows, i), k = INTEGER(ci)[r] - 1;
    if (LOGICAL(U)[r]) REAL(aif)[i + B*k] = quo(mul(REAL(w)[r], sub(INTEGER(ai)[r] > 0 ? 1.0 : 0.0, REAL(a)[k])), REAL(mass)[k]);
  }
  for (R_xlen_t jj = 0; jj < XLENGTH(sids); ++jj) {
    int sj = INTEGER(sids)[jj] - 1;
    SEXP spec = VECTOR_ELT(ss, sj), fit = VECTOR_ELT(sf, sj);
    int k = INTEGER(item(spec, "cohort"))[0], n = 0;
    for (int i = 0; i < B; ++i) { int r = row_at(rows, i); if (INTEGER(ci)[r] == k && LOGICAL(U)[r]) rr[n++] = i; }
    SEXP xi = PROTECT(zero_matrix(B, d)), pr = PROTECT(zero_vector(B)), S = PROTECT(zero_vector(B));
    SET_VECTOR_ELT(alpha, sj, xi); SET_VECTOR_ELT(probs, sj, pr); SET_VECTOR_ELT(responses, sj, S);
    if (n) {
      SEXP rs = PROTECT(subset_rows(rows, rr, n));
      SEXP call = PROTECT(Rf_lang5(predict, reader, rs, spec, fit));
      SEXP prob = PROTECT(Rf_eval(call, R_BaseEnv));
      SEXP hit = PROTECT(pair_call(pair, rs, item(spec, "pair"), 0));
      SEXP resp = item(hit, "S"), score = PROTECT(zero_matrix(n, d));
      for (int i = 0; i < n; ++i) {
        int b = rr[i], r = row_at(rows, b);
        REAL(pr)[b] = REAL(prob)[i]; REAL(S)[b] = REAL(resp)[i];
        double factor = mul(REAL(w)[r], sub(REAL(resp)[i], REAL(prob)[i]));
        for (int h = 0; h < d; ++h) REAL(score)[i + n*h] = mul(factor, REAL(Z)[r + N*h]);
      }
      SEXP score_t = PROTECT(eval1(transpose, score));
      SEXP solved = PROTECT(eval2(solve, item(fit, "H_chol"), score_t));
      SEXP value = PROTECT(eval1(transpose, solved));
      for (int h = 0; h < d; ++h) for (int i = 0; i < n; ++i) REAL(xi)[rr[i] + B*h] = REAL(value)[i + n*h];
      UNPROTECT(8);
    }
    UNPROTECT(3);
  }
  for (R_xlen_t jj = 0; jj < XLENGTH(tids); ++jj) {
    int tj = INTEGER(tids)[jj] - 1;
    SEXP spec = VECTOR_ELT(ts, tj), fit = VECTOR_ELT(tf, tj), cohorts = item(spec, "cohorts");
    SEXP ids = item(spec, "selection_ids"), means = item(fit, "means"), part = item(fit, "partials");
    int m = (int)XLENGTH(cohorts), u = m * (d + 1), n = 0;
    SEXP upstream = PROTECT(zero_matrix(B, u)), S = PROTECT(zero_vector(B)), pp = PROTECT(zero_vector(B));
    for (int k = 0; k < m; ++k) {
      int c = INTEGER(cohorts)[k] - 1, sj = INTEGER(ids)[k] - 1, ix = k * (d + 1);
      SEXP xi = VECTOR_ELT(alpha, sj), pr = VECTOR_ELT(probs, sj), resp = VECTOR_ELT(responses, sj);
      for (int i = 0; i < B; ++i) {
        REAL(upstream)[i + B*ix] = REAL(aif)[i + B*c];
        for (int h = 0; h < d; ++h) REAL(upstream)[i + B*(ix + h + 1)] = REAL(xi)[i + B*h];
        if (INTEGER(ci)[row_at(rows, i)] == c + 1) { REAL(S)[i] = REAL(resp)[i]; REAL(pp)[i] = REAL(pr)[i]; }
      }
    }
    for (int i = 0; i < B; ++i) if (REAL(S)[i] == 1) rr[n++] = i;
    SEXP beta = PROTECT(zero_matrix(B, p)), direct = PROTECT(zero_matrix(B, 4));
    for (int i = 0; i < B; ++i) { REAL(direct)[i + 2*B] = -REAL(means)[2]; REAL(direct)[i + 3*B] = -REAL(means)[3]; }
    if (n) {
      /* Exactly frame$Z_a[st$a_index[rows[rr]], , drop = FALSE], before
         the treatment pair query. R owns this bounded source read, including
         legitimate base wrappers; no full-frame data pointer is requested. */
      SEXP ar = PROTECT(Rf_allocVector(INTSXP, n));
      for (int i = 0; i < n; ++i) INTEGER(ar)[i] = INTEGER(ai)[row_at(rows, rr[i])];
      SEXP drop = PROTECT(Rf_ScalarLogical(0));
      SEXP subset = PROTECT(Rf_lang5(Rf_findFun(Rf_install("["), R_BaseEnv),
        Xall, ar, R_MissingArg, drop));
      SET_TAG(CDR(CDR(CDR(CDR(subset)))), Rf_install("drop"));
      SEXP X = PROTECT(Rf_eval(subset, R_BaseEnv));
      SEXP rs = PROTECT(subset_rows(rows, rr, n));
      SEXP hit = PROTECT(pair_call(pair, rs, item(spec, "pair"), 1));
      SEXP lv = PROTECT(eval2(logit, X, item(fit, "coefficients")));
      SEXP values = PROTECT(Rf_allocVector(REALSXP, 3*n));
      const char *fields[] = {"prob", "eta", "odds"};
      for (int h = 0; h < 3; ++h) memcpy(REAL(values) + h*n, REAL(item(lv, fields[h])), (size_t)n * sizeof(double));
      finite_call(finite, values, "treatment replay");
      SEXP pp_sub = PROTECT(Rf_allocVector(REALSXP, n));
      for (int i = 0; i < n; ++i) REAL(pp_sub)[i] = REAL(pp)[rr[i]];
      SEXP clipped = PROTECT(eval1(clip_floor, pp_sub));
      SEXP prob = item(lv, "prob"), odds = item(lv, "odds"), dy = item(hit, "dy");
      for (int i = 0; i < n; ++i) {
        int b = rr[i], r = row_at(rows, b), c = INTEGER(ci)[r] - 1;
        double D = c + 1 == INTEGER(cohorts)[0] ? 1.0 : 0.0;
        double fw = quo(mul(REAL(w)[r], REAL(a)[c]), mul(REAL(rho)[r], REAL(clipped)[i]));
        double wp = mul(fw, D), wm = mul(mul(fw, sub(1.0, D)), REAL(odds)[i]);
        double factor = mul(fw, sub(D, REAL(prob)[i]));
        for (int h = 0; h < p; ++h) REAL(beta)[b + B*h] = mul(factor, REAL(X)[i + n*h]);
        REAL(direct)[b] = quo(mul(wp, sub(REAL(dy)[i], REAL(means)[0])), REAL(means)[2]);
        REAL(direct)[b + B] = quo(mul(wm, sub(REAL(dy)[i], REAL(means)[1])), REAL(means)[3]);
        REAL(direct)[b + 2*B] = add(REAL(direct)[b + 2*B], wp);
        REAL(direct)[b + 3*B] = add(REAL(direct)[b + 3*B], wm);
      }
      UNPROTECT(10);
    }
    SEXP cross_t = PROTECT(eval1(transpose, item(fit, "cross")));
    SEXP cross_product = PROTECT(eval2(product, upstream, cross_t));
    /* Keep the original matrix RHS and its two additions, with no projection. */
    SEXP rhs = PROTECT(zero_matrix(B, p));
    for (R_xlen_t i = 0; i < XLENGTH(rhs); ++i) REAL(rhs)[i] = add(REAL(beta)[i], REAL(cross_product)[i]);
    SEXP rhs_t = PROTECT(eval1(transpose, rhs));
    SEXP solved = PROTECT(eval2(solve, item(fit, "H_chol"), rhs_t));
    SEXP xb = PROTECT(eval1(transpose, solved));
    SEXP up_t = PROTECT(eval1(transpose, item(part, "upstream")));
    SEXP up_product = PROTECT(eval2(product, upstream, up_t));
    SEXP coord = PROTECT(zero_matrix(B, 4));
    for (R_xlen_t i = 0; i < XLENGTH(coord); ++i) REAL(coord)[i] = add(REAL(direct)[i], REAL(up_product)[i]);
    SEXP beta_t = PROTECT(eval1(transpose, item(part, "beta")));
    SEXP beta_product = PROTECT(eval2(product, xb, beta_t));
    for (R_xlen_t i = 0; i < XLENGTH(coord); ++i) REAL(coord)[i] = add(REAL(coord)[i], REAL(beta_product)[i]);
    SEXP values = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)B*(p + 4)));
    memcpy(REAL(values), REAL(xb), (size_t)B*p*sizeof(double));
    memcpy(REAL(values) + (R_xlen_t)B*p, REAL(coord), (size_t)B*4*sizeof(double));
    finite_call(finite, values, "complete frame IF");
    SEXP phi = PROTECT(Rf_allocVector(REALSXP, B)), model = PROTECT(Rf_allocVector(VECSXP, 1));
    for (int i = 0; i < B; ++i) REAL(phi)[i] = sub(REAL(coord)[i], REAL(coord)[i + B]);
    SET_VECTOR_ELT(model, 0, phi);
    SEXP names = PROTECT(Rf_mkString("phi")); Rf_setAttrib(model, R_NamesSymbol, names);
    SET_VECTOR_ELT(models, tj, model);
    UNPROTECT(20);
  }
  UNPROTECT(5);
  return models;
}
