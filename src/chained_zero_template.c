#include <R.h>
#include <Rinternals.h>
#include "cdid_r_api.h"
#include <limits.h>
#include <string.h>

/* Metadata only. Numeric bits are copied; no score, RHS or solve is evaluated. */
static int zr_plain(SEXP x, SEXPTYPE type) {
  return TYPEOF(x) == type && !Rf_isObject(x) && !ALTREP(x) && !cdid_has_attrib(x);
}

static int zr_named_list(SEXP x) {
  if (TYPEOF(x) != VECSXP || Rf_isObject(x) || ALTREP(x)) return 0;
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  return cdid_single_attrib(x, R_NamesSymbol) &&
    TYPEOF(names) == STRSXP && !Rf_isObject(names) && !ALTREP(names) &&
    !cdid_has_attrib(names) && XLENGTH(names) == XLENGTH(x);
}

static SEXP zr_field(SEXP x, const char *name) {
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
    SEXP s = STRING_ELT(names, i);
    if (s != NA_STRING && strcmp(CHAR(s), name) == 0) return VECTOR_ELT(x, i);
  }
  return R_NilValue;
}

static int zr_active(SEXP active, int gid) {
  for (R_xlen_t i = 0; i < XLENGTH(active); ++i)
    if (INTEGER(active)[i] == gid) return 1;
  return 0;
}

/* Nonalias identical() operands can read ALTREP data. Decline those private
   operands before a suffix can run callbacks between checks and copies.
   Source aliases need no traversal: R_compute_identical returns immediately
   for them. This is admission only; it never replaces an equality check. */
static int zr_materialized(SEXP x, int depth, int *left);
typedef struct { int depth, *left, ok; } zr_attr_walk;
static SEXP zr_materialized_attr(SEXP tag, SEXP value, void *data) {
  zr_attr_walk *walk = (zr_attr_walk *) data;
  (void) tag;
  /* Preserve the old attribute pairlist's node/depth charges: one cons
     node followed by its value, then the next cons at depth + 1. */
  if (walk->depth > 32 || --*walk->left < 0 ||
      !zr_materialized(value, walk->depth + 1, walk->left)) {
    walk->ok = 0;
    return R_NilValue;
  }
  ++walk->depth;
  return NULL;
}
static int zr_materialized(SEXP x, int depth, int *left) {
  if (x == R_NilValue) return 1;
  if (depth > 32 || --*left < 0 || ALTREP(x)) return 0;
  switch (TYPEOF(x)) {
  case PROMSXP: case DOTSXP: case BCODESXP: return 0;
  default: break;
  }
  if (cdid_has_attrib(x)) {
    zr_attr_walk attrs = {depth + 1, left, 1};
    cdid_each_attrib(x, zr_materialized_attr, &attrs);
    if (!attrs.ok) return 0;
  }
  switch (TYPEOF(x)) {
  case VECSXP: case EXPRSXP:
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
      if (!zr_materialized(VECTOR_ELT(x, i), depth + 1, left)) return 0;
    break;
  case LISTSXP: case LANGSXP:
    if (!zr_materialized(CAR(x), depth + 1, left) ||
        !zr_materialized(CDR(x), depth + 1, left)) return 0;
    break;
  case CLOSXP:
    if (!zr_materialized(cdid_closure_formals(x), depth + 1, left) ||
        !zr_materialized(cdid_closure_body(x), depth + 1, left)) return 0;
    break;
  default: break; /* Environments are compared by identity, never inspected. */
  }
  return 1;
}

static int zr_ready(SEXP x, SEXP y) {
  if (x == y) return 1;
  int left = 4096;
  return zr_materialized(x, 0, &left) && zr_materialized(y, 0, &left);
}

/* A failed/unsupported group is not included. The caller's scalar helper,
   RHS, solve and template update process that group before any later check. */
static int zr_hit(SEXP templates, SEXP groups, SEXP cells, SEXP active,
                  SEXP N, int gid, SEXP *ids, SEXP *lambda, SEXP *nu) {
  if (zr_active(active, gid)) return 0;
  SEXP group = VECTOR_ELT(groups, gid - 1), x = VECTOR_ELT(templates, gid - 1);
  if (!zr_named_list(group) || !zr_named_list(x)) return 0;
  SEXP cids = zr_field(group, "cell_ids"), p = zr_field(group, "p");
  if (!zr_plain(cids, INTSXP) || !XLENGTH(cids) || !zr_plain(p, INTSXP) ||
      XLENGTH(p) != 1 || INTEGER(p)[0] < 0) return 0;
  for (R_xlen_t i = 0; i < XLENGTH(cids); ++i)
    if (INTEGER(cids)[i] < 1 || INTEGER(cids)[i] > XLENGTH(cells)) return 0;
  SEXP lam = zr_field(x, "lambda"), v = zr_field(x, "nu"), old_N = zr_field(x, "N");
  if (!zr_plain(lam, REALSXP) || XLENGTH(lam) != XLENGTH(cids) ||
      !zr_plain(v, REALSXP) || XLENGTH(v) != INTEGER(p)[0] ||
      !zr_named_list(zr_field(x, "group")) || !zr_plain(zr_field(x, "cells"), VECSXP) ||
      !((zr_plain(old_N, INTSXP) || zr_plain(old_N, REALSXP)) && XLENGTH(old_N) == 1)) return 0;
  if (!zr_ready(group, zr_field(x, "group"))) return 0;
  /* R's identical defaults use IDENT_USE_CLOENV (16), including closure
     environment equality, signed zero and the ordinary NA/NaN conventions. */
  if (!R_compute_identical(group, zr_field(x, "group"), IDENT_USE_CLOENV)) return 0;
  SEXP old_cells = zr_field(x, "cells");
  if (XLENGTH(old_cells) == XLENGTH(cids)) for (R_xlen_t i = 0; i < XLENGTH(cids); ++i)
    if (!zr_ready(VECTOR_ELT(cells, INTEGER(cids)[i] - 1), VECTOR_ELT(old_cells, i))) return 0;
  SEXP subset = PROTECT(Rf_allocVector(VECSXP, XLENGTH(cids)));
  for (R_xlen_t i = 0; i < XLENGTH(cids); ++i)
    SET_VECTOR_ELT(subset, i, VECTOR_ELT(cells, INTEGER(cids)[i] - 1));
  int same_cells = R_compute_identical(subset, zr_field(x, "cells"), IDENT_USE_CLOENV);
  UNPROTECT(1);
  if (!same_cells || !R_compute_identical(N, zr_field(x, "N"), IDENT_USE_CLOENV)) return 0;
  *ids = cids; *lambda = lam; *nu = v;
  return 1;
}

SEXP C_gmm_reconstruction_zero_run(SEXP templates, SEXP groups, SEXP cells,
                                  SEXP active, SEXP N, SEXP cursor) {
  if (!zr_plain(templates, VECSXP) || !zr_plain(groups, VECSXP) || !zr_plain(cells, VECSXP) ||
      !zr_plain(active, INTSXP) || !zr_plain(cursor, INTSXP) || XLENGTH(cursor) != 1 ||
      XLENGTH(groups) > INT_MAX || XLENGTH(cells) > INT_MAX || XLENGTH(templates) != XLENGTH(groups) ||
      !((zr_plain(N, INTSXP) || zr_plain(N, REALSXP)) && XLENGTH(N) == 1)) return R_NilValue;
  int start = INTEGER(cursor)[0], ng = (int) XLENGTH(groups);
  if (start < 1 || start > ng) return R_NilValue;
  for (R_xlen_t i = 0; i < XLENGTH(active); ++i)
    if (INTEGER(active)[i] < 1 || INTEGER(active)[i] > ng) return R_NilValue;
  int stop = start; R_xlen_t nc = 0, np = 0;
  SEXP ids = R_NilValue, lambda = R_NilValue, nu = R_NilValue;
  for (; stop <= ng; ++stop) {
    if (!zr_hit(templates, groups, cells, active, N, stop, &ids, &lambda, &nu)) break;
    /* The existing forcing/map slabs fund at most nc cell slots and the
       issued group widths. Oversized private incidences keep scalar replay. */
    if (XLENGTH(ids) > XLENGTH(cells) - nc || XLENGTH(nu) > INT_MAX - np) break;
    nc += XLENGTH(ids); np += XLENGTH(nu);
    if (stop == INT_MAX) return R_NilValue;
  }
  int count = stop - start;
  if (!count) return R_NilValue;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
  const char *labels[] = {"cell_ids", "lambda", "group_ids", "nu", "next_gid"};
  for (int i = 0; i < 5; ++i) SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_setAttrib(out, R_NamesSymbol, names);
  SEXP ci = PROTECT(Rf_allocVector(INTSXP, nc));
  SEXP lam = PROTECT(Rf_allocVector(REALSXP, nc));
  SEXP gi = PROTECT(Rf_allocVector(INTSXP, count));
  SEXP nus = PROTECT(Rf_allocVector(VECSXP, count));
  R_xlen_t offset = 0;
  /* No R evaluation/callback occurs between successful checks and these
     fresh copies. Only the already accepted contiguous prefix is visited. */
  for (int gid = start; gid < stop; ++gid) {
    SEXP group = VECTOR_ELT(groups, gid - 1), x = VECTOR_ELT(templates, gid - 1);
    ids = zr_field(group, "cell_ids"); lambda = zr_field(x, "lambda"); nu = zr_field(x, "nu");
    if (XLENGTH(ids)) {
      memcpy(INTEGER(ci) + offset, INTEGER(ids), (size_t) XLENGTH(ids) * sizeof(int));
      memcpy(REAL(lam) + offset, REAL(lambda), (size_t) XLENGTH(ids) * sizeof(double));
    }
    offset += XLENGTH(ids); INTEGER(gi)[gid - start] = gid;
    SEXP v = PROTECT(Rf_allocVector(REALSXP, XLENGTH(nu)));
    if (XLENGTH(nu)) memcpy(REAL(v), REAL(nu), (size_t) XLENGTH(nu) * sizeof(double));
    SET_VECTOR_ELT(nus, gid - start, v); UNPROTECT(1);
  }
  SET_VECTOR_ELT(out, 0, ci); SET_VECTOR_ELT(out, 1, lam); SET_VECTOR_ELT(out, 2, gi); SET_VECTOR_ELT(out, 3, nus);
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(stop));
  UNPROTECT(6); return out;
}

/* V53 bridges preserve the existing V50 decision and materialization rules. */
int cdid_reconstruction_zero_hit(SEXP templates, SEXP groups, SEXP cells, SEXP active,
                                  SEXP N, int gid, SEXP *ids, SEXP *lambda, SEXP *nu) {
  return zr_hit(templates, groups, cells, active, N, gid, ids, lambda, nu);
}
int cdid_reconstruction_materialized(SEXP x) {
  int left = 65536;
  return zr_materialized(x, 0, &left);
}
