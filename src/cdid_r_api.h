#ifndef CDID_R_API_H
#define CDID_R_API_H

#include <Rinternals.h>
#include <Rversion.h>
#include <string.h>

/* Raw attributes, without attributes()'s allocation or names/row.names
   transformations. NULL continues a visit; any SEXP (even R_NilValue) stops. */
typedef SEXP (*cdid_attrib_visitor)(SEXP, SEXP, void *);
typedef struct { cdid_attrib_visitor visit; void *data; } cdid_attrib_call;

static inline int cdid_has_attrib(SEXP x) {
#if R_VERSION >= R_Version(4, 5, 0)
  return ANY_ATTRIB(x) != 0;
#else
  return ATTRIB(x) != R_NilValue;
#endif
}

static inline SEXP cdid_attrib_call_visit(SEXP tag, SEXP value, void *data) {
  cdid_attrib_call *call = (cdid_attrib_call *) data;
  PROTECT(tag); PROTECT(value);
  SEXP out = call->visit(tag, value, call->data);
  UNPROTECT(2);
  return out;
}

static inline SEXP cdid_each_attrib(SEXP x, cdid_attrib_visitor visit, void *data) {
  cdid_attrib_call call = {visit, data};
  PROTECT(x);
#if R_VERSION >= R_Version(4, 6, 0)
  SEXP out = R_mapAttrib(x, cdid_attrib_call_visit, &call);
#else
  SEXP out = NULL;
  for (SEXP a = ATTRIB(x); a != R_NilValue; a = CDR(a)) {
    PROTECT(a);
    out = cdid_attrib_call_visit(TAG(a), CAR(a), &call);
    UNPROTECT(1);
    if (out != NULL) break;
  }
#endif
  UNPROTECT(1);
  return out;
}

typedef struct { const SEXP *tags; int n, seen, single; } cdid_attrib_tags;
static inline SEXP cdid_attrib_tags_visit(SEXP tag, SEXP value, void *data) {
  cdid_attrib_tags *tags = (cdid_attrib_tags *) data;
  (void) value;
  if (tags->single && tags->seen++) return R_NilValue;
  for (int i = 0; i < tags->n; ++i) if (tag == tags->tags[i]) return NULL;
  return R_NilValue;
}

/* Exact singleton, or an unrestricted number of the permitted raw tags. */
static inline int cdid_single_attrib(SEXP x, SEXP tag) {
  cdid_attrib_tags tags = {&tag, 1, 0, 1};
  return cdid_each_attrib(x, cdid_attrib_tags_visit, &tags) == NULL && tags.seen == 1;
}
static inline int cdid_attrib_only(SEXP x, const SEXP *allowed, int n) {
  cdid_attrib_tags tags = {allowed, n, 0, 0};
  return cdid_each_attrib(x, cdid_attrib_tags_visit, &tags) == NULL;
}

static inline SEXP cdid_closure_formals(SEXP x) {
#if R_VERSION >= R_Version(4, 5, 0)
  return R_ClosureFormals(x);
#else
  return FORMALS(x);
#endif
}
static inline SEXP cdid_closure_body(SEXP x) {
#if R_VERSION >= R_Version(4, 5, 0)
  return R_ClosureBody(x);
#else
  return BODY(x);
#endif
}
static inline SEXP cdid_closure_env(SEXP x) {
#if R_VERSION >= R_Version(4, 5, 0)
  return R_ClosureEnv(x);
#else
  return CLOENV(x);
#endif
}
static inline SEXP cdid_parent_env(SEXP x) {
#if R_VERSION >= R_Version(4, 5, 0)
  return R_ParentEnv(x);
#else
  return ENCLOS(x);
#endif
}

/* Distinguish an absent binding from a present but unavailable value: only
   absence permits base-frame fallback. Never evaluate active/delayed bindings. */
static inline SEXP cdid_frame_value(SEXP env, SEXP key, int *bound) {
#if R_VERSION >= R_Version(4, 6, 0)
  R_BindingType_t type = R_GetBindingType(key, env);
  *bound = type != R_BindingTypeUnbound;
  if (type != R_BindingTypeValue && type != R_BindingTypeForced) return R_UnboundValue;
  return R_getVarEx(key, env, FALSE, R_UnboundValue);
#else
  /* Existence must precede R_BindingIsActive (which errors for an absent
     binding) and the legacy lookup (which itself invokes active bindings).
     R 4.1 lacks R_existsVarInFrame; listing names does not read their values. */
#if R_VERSION >= R_Version(4, 2, 0)
  *bound = R_existsVarInFrame(env, key);
#else
  SEXP names = PROTECT(R_lsInternal3(env, TRUE, FALSE));
  *bound = 0;
  for (R_xlen_t i = 0; i < XLENGTH(names); ++i)
    if (!strcmp(CHAR(STRING_ELT(names, i)), CHAR(PRINTNAME(key)))) { *bound = 1; break; }
  UNPROTECT(1);
#endif
  if (!*bound || R_BindingIsActive(key, env)) return R_UnboundValue;
  SEXP x = Rf_findVarInFrame(env, key);
  if (x == R_MissingArg) return R_UnboundValue;
  return TYPEOF(x) == PROMSXP ? PRVALUE(x) : x;
#endif
}

#endif
