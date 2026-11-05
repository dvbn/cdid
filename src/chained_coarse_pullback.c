/* V53: original scalar group assembly plus score contraction, one coarse call.
   Original Gram helpers and V50 first-hit/suffix decisions remain authoritative. */
#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ > 0)) || \
    (defined(__ASSOC_MATH__) && (__ASSOC_MATH__ > 0)) || \
    (defined(__RECIPROCAL_MATH__) && (__RECIPROCAL_MATH__ > 0))
#error "Unsafe arithmetic is forbidden for coarse pullbacks."
#endif
extern void cdid_pair_operation_arithmetic(void);
extern void cdid_pair_operation_error(const char *, const char *);
extern int cdid_reconstruction_materialized(SEXP);
extern int cdid_reconstruction_zero_hit(SEXP, SEXP, SEXP, SEXP, SEXP, int, SEXP *, SEXP *, SEXP *);
extern SEXP C_gmm_reconstruction_score_pullback(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

typedef struct { int B, nc, ng, peak, intercept, Ncount; double total_p, cap; } coarse53;
static int cr_plain(SEXP x, int type) { return TYPEOF(x)==type && !OBJECT(x) && !ALTREP(x) && ATTRIB(x)==R_NilValue; }
static int cr_list(SEXP x) {
  if(TYPEOF(x)!=VECSXP || OBJECT(x) || ALTREP(x)) return 0;
  SEXP a=ATTRIB(x), n=Rf_getAttrib(x,R_NamesSymbol);
  return a!=R_NilValue && TAG(a)==R_NamesSymbol && CDR(a)==R_NilValue && cr_plain(n,STRSXP) && XLENGTH(n)==XLENGTH(x);
}
static SEXP cr_field(SEXP x,const char *key) {
  SEXP n=Rf_getAttrib(x,R_NamesSymbol);
  for(R_xlen_t i=0;i<XLENGTH(x);++i) if(STRING_ELT(n,i)!=NA_STRING && !strcmp(CHAR(STRING_ELT(n,i)),key)) return VECTOR_ELT(x,i);
  return R_NilValue;
}
static int cr_int(SEXP x,int min,int *v) {
  if(!((cr_plain(x,INTSXP)||cr_plain(x,REALSXP)) && XLENGTH(x)==1)) return 0;
  double d=TYPEOF(x)==INTSXP?INTEGER(x)[0]:REAL(x)[0];
  if(!R_FINITE(d)||d<min||d>INT_MAX||d!=floor(d)) return 0; *v=(int)d; return 1;
}
static int cr_ids(SEXP x,int upper,int ordered) {
  if(!cr_plain(x,INTSXP)) return 0; int old=0;
  for(R_xlen_t i=0;i<XLENGTH(x);++i) {int n=INTEGER(x)[i];if(n<1||n>upper||(ordered&&n<=old))return 0;old=n;} return 1;
}
static int cr_num(SEXP x,int n,SEXP columns) {
  if(TYPEOF(x)!=REALSXP||OBJECT(x)||ALTREP(x)||XLENGTH(x)!=n)return 0;
  if(ATTRIB(x)==R_NilValue)return 1;
  SEXP a=ATTRIB(x), nm=Rf_getAttrib(x,R_NamesSymbol);
  return TAG(a)==R_NamesSymbol&&CDR(a)==R_NilValue&&cr_plain(nm,STRSXP)&&cr_plain(columns,STRSXP)&&
    XLENGTH(nm)==n&&XLENGTH(columns)==n&&R_compute_identical(nm,columns,IDENT_USE_CLOENV);
}
static int cr_matrix(SEXP x,int n,int p) {
  if(TYPEOF(x)!=REALSXP||OBJECT(x)||ALTREP(x))return 0;
  for(SEXP a=ATTRIB(x);a!=R_NilValue;a=CDR(a))if(TAG(a)!=R_DimSymbol&&TAG(a)!=R_DimNamesSymbol)return 0;
  SEXP d=Rf_getAttrib(x,R_DimSymbol);
  return cr_plain(d,INTSXP)&&XLENGTH(d)==2&&INTEGER(d)[0]==n&&INTEGER(d)[1]==p&&XLENGTH(x)==(R_xlen_t)n*p;
}
static int cr_functions(SEXP f) {
  if(!cr_list(f))return 0; SEXP owner=cr_field(f,".owner"),names=Rf_getAttrib(f,R_NamesSymbol);
  if(TYPEOF(owner)!=ENVSXP)return 0;
  for(R_xlen_t i=0;i<XLENGTH(f);++i) {
    const char *name=CHAR(STRING_ELT(names,i));if(!strcmp(name,".owner")||!strcmp(name,".live"))continue;
    SEXP sym=Rf_install(name),actual=Rf_findVarInFrame(owner,sym);
    if(actual!=R_UnboundValue&&R_BindingIsActive(sym,owner))return 0;
    if(actual==R_UnboundValue)actual=Rf_findVarInFrame(R_BaseEnv,sym);
    if(TYPEOF(actual)==PROMSXP)actual=PRVALUE(actual);
    if(actual!=VECTOR_ELT(f,i))return 0;
  }
  return 1;
}

/* Forced values only: a private/unforced/active binding never runs a callback
   to obtain permission. Canonical factory/core formals have already been used. */
static int cr_bound(SEXP env,SEXP key) {
#if R_VERSION >= R_Version(4, 5, 0)
  return R_existsVarInFrame(env,key);
#else
  /* This optional method uses an active-safe existence API absent in R4.1.
     Older supported R keeps the unchanged V57/full recognizer path. */
  (void)env;(void)key;return 0;
#endif
}
static SEXP cr_value(SEXP env,const char *name) {
  if(TYPEOF(env)!=ENVSXP||OBJECT(env))return R_UnboundValue;
  SEXP key=Rf_install(name);
  if(!cr_bound(env,key)||R_BindingIsActive(key,env))return R_UnboundValue;
  SEXP x=Rf_findVarInFrame(env,key);
  if(TYPEOF(x)==PROMSXP)x=PRVALUE(x);
  return TYPEOF(x)==PROMSXP?R_UnboundValue:x;
}
static const char *cr_factory_names[]={"read","plan","n_cells","position","cache","seen","reader","view"};
static const char *cr_live_primitives[]={"typeof","is.object","is.environment","parent.env","environment","sort","ls",
 "bindingIsActive","get","attr","identical","body","quote","is.null","is.call","length","as.name",
 "is.integer","is.logical","is.list","attributes","c","!","!=","==","<","||","&&","[","[[","$","attr<-","isTRUE"};
#define CR_LIVE_PRIMITIVES ((int)(sizeof(cr_live_primitives)/sizeof(cr_live_primitives[0])))
static const char *cr_live_helpers[]={"reconstruction_score_cache_view","reconstruction_cache_reader_body","reconstruction_sparse_adjoint_ids"};
static SEXP cr_inherited_value(SEXP owner,const char *name) {
  return cr_bound(owner,Rf_install(name))?cr_value(owner,name):cr_value(R_BaseEnv,name);
}
static int cr_factory_names_ok(SEXP env) {
  if(TYPEOF(env)!=ENVSXP||OBJECT(env))return 0;
  SEXP n=PROTECT(R_lsInternal3(env,TRUE,FALSE));
  int ok=TYPEOF(n)==STRSXP&&XLENGTH(n)==8;
  if(ok)for(int j=0;j<8;++j){int found=0;for(int k=0;k<8;++k)
    if(STRING_ELT(n,k)!=NA_STRING&&!strcmp(CHAR(STRING_ELT(n,k)),cr_factory_names[j]))++found;
    if(found!=1){ok=0;break;}}
  UNPROTECT(1);return ok;
}
static int cr_same(SEXP a,SEXP b){return R_compute_identical(a,b,IDENT_USE_CLOENV);}
static int cr_source_tag(SEXP tag){return tag==Rf_install("srcref")||tag==Rf_install("srcfile")||tag==Rf_install("wholeSrcref");}
static int cr_code_attrs(SEXP x,SEXP y,int strip) {
  SEXP a=ATTRIB(x),b=ATTRIB(y);
  for(;;){while(strip&&a!=R_NilValue&&cr_source_tag(TAG(a)))a=CDR(a);
    while(strip&&b!=R_NilValue&&cr_source_tag(TAG(b)))b=CDR(b);
    if(a==R_NilValue||b==R_NilValue)return a==b;
    if(TAG(a)!=TAG(b)||!cr_same(CAR(a),CAR(b)))return 0;a=CDR(a);b=CDR(b);}
}
/* Exactly the two annotation locations of the existing R normalizer. Other
   executable nodes, pair tags and all other attributes remain exact. */
static int cr_reader_code(SEXP x,SEXP y,int outer) {
  if(TYPEOF(x)!=LANGSXP||TYPEOF(y)!=LANGSXP||!cr_code_attrs(x,y,1))return 0;
  SEXP a=x,b=y;int count=0;
  for(;a!=R_NilValue&&b!=R_NilValue;a=CDR(a),b=CDR(b),++count){
    if(TYPEOF(a)!=LANGSXP&&TYPEOF(a)!=LISTSXP)return 0;
    if(TYPEOF(b)!=LANGSXP&&TYPEOF(b)!=LISTSXP)return 0;
    if(count>4||TAG(a)!=TAG(b)||(count&&!cr_code_attrs(a,b,0)))return 0;
    if(outer&&count==3){SEXP u=CAR(a),v=CAR(b);
      if(TYPEOF(u)!=LANGSXP||TYPEOF(v)!=LANGSXP||!cr_code_attrs(u,v,0))return 0;
      SEXP up=u,vp=v;
      for(int k=0;k<3;++k){if(up==R_NilValue||vp==R_NilValue||TAG(up)!=TAG(vp)||(k&&!cr_code_attrs(up,vp,0)))return 0;
        if(k==2){if(!cr_reader_code(CAR(up),CAR(vp),0))return 0;}
        else if(!cr_same(CAR(up),CAR(vp)))return 0;up=CDR(up);vp=CDR(vp);}
      if(up!=R_NilValue||vp!=R_NilValue)return 0;
    }else if(!cr_same(CAR(a),CAR(b)))return 0;
  }
  return a==R_NilValue&&b==R_NilValue&&count==(outer?5:4);
}
/* Detached canonical code copies are capped before allocation. Symbols and
   source environments stay borrowed; normalized code has no attributes. */
static int cr_code_cost(SEXP x,int depth,double *cost) {
  if(depth>64||*cost>640)return 0;
  if(x==R_NilValue||TYPEOF(x)==SYMSXP)return 1;
  if(ATTRIB(x)!=R_NilValue)return 0;
  if(TYPEOF(x)==LANGSXP||TYPEOF(x)==LISTSXP){*cost+=8;
    return cr_code_cost(CAR(x),depth+1,cost)&&cr_code_cost(CDR(x),depth+1,cost);}
  if((TYPEOF(x)==INTSXP||TYPEOF(x)==REALSXP||TYPEOF(x)==LGLSXP)&&!ALTREP(x)){
    *cost+=8+XLENGTH(x);return *cost<=640;}
  return 0;
}
static SEXP cr_live_create(SEXP view,SEXP state,SEXP functions) {
  SEXP meta=cr_field(functions,".live"),reader=cr_field(view,"reader"),read=cr_field(view,"read");
  if(!cr_list(meta)||TYPEOF(reader)!=CLOSXP||TYPEOF(read)!=CLOSXP)return R_NilValue;
  SEXP env=CLOENV(reader),core=CLOENV(read),owner=cr_field(functions,".owner");
  if(TYPEOF(env)!=ENVSXP||OBJECT(env)||TYPEOF(core)!=ENVSXP||OBJECT(core)||ENCLOS(env)!=owner||ENCLOS(core)!=owner||!cr_factory_names_ok(env))return R_NilValue;
  SEXP vf=cr_value(env,"view"),witness=cr_value(core,"complete_cache"),rb=cr_field(meta,"reader_body"),vb=cr_field(meta,"view_body"),helpers=cr_field(meta,"functions");
  double cost=0;
  if(TYPEOF(vf)!=CLOSXP||!cr_list(witness)||!cr_list(helpers)||XLENGTH(helpers)!=3||
    !cr_reader_code(BODY_EXPR(reader),rb,1)||!cr_same(BODY_EXPR(vf),vb)||
    !cr_code_cost(rb,0,&cost)||!cr_code_cost(vb,0,&cost)||cost>640||
    !cr_same(cr_field(witness,"view"),view)||!cr_same(cr_field(witness,"state"),state))return R_NilValue;
  for(int j=0;j<3;++j)if(cr_value(owner,cr_live_helpers[j])!=cr_field(helpers,cr_live_helpers[j]))return R_NilValue;
  SEXP out=PROTECT(Rf_allocVector(VECSXP,12)),values=PROTECT(Rf_allocVector(VECSXP,8)),prims=PROTECT(Rf_allocVector(VECSXP,CR_LIVE_PRIMITIVES));
  SET_VECTOR_ELT(out,0,env);SET_VECTOR_ELT(out,1,core);SET_VECTOR_ELT(out,2,values);
  for(int j=0;j<8;++j){SEXP x=cr_value(env,cr_factory_names[j]);if(x==R_UnboundValue){UNPROTECT(3);return R_NilValue;}SET_VECTOR_ELT(values,j,x);}
  for(int j=0;j<CR_LIVE_PRIMITIVES;++j){SEXP x=cr_inherited_value(owner,cr_live_primitives[j]);
    if(x==R_UnboundValue){UNPROTECT(3);return R_NilValue;}SET_VECTOR_ELT(prims,j,x);}
  SET_VECTOR_ELT(out,3,Rf_duplicate(rb));SET_VECTOR_ELT(out,4,Rf_duplicate(vb));SET_VECTOR_ELT(out,5,helpers);
  SET_VECTOR_ELT(out,6,prims);SET_VECTOR_ELT(out,7,witness);SET_VECTOR_ELT(out,8,BODY_EXPR(read));
  SET_VECTOR_ELT(out,9,ATTRIB(reader));SET_VECTOR_ELT(out,10,ATTRIB(vf));SET_VECTOR_ELT(out,11,ATTRIB(read));
  UNPROTECT(3);return out;
}
static int cr_live_ok(SEXP live,SEXP reader,SEXP view,SEXP state,SEXP functions) {
  if(TYPEOF(live)!=VECSXP||XLENGTH(live)!=12||(TYPEOF(reader)!=CLOSXP||OBJECT(reader)))return 0;
  SEXP env=VECTOR_ELT(live,0),core=VECTOR_ELT(live,1),values=VECTOR_ELT(live,2),owner=cr_field(functions,".owner");
  if(OBJECT(env)||OBJECT(core)||ENCLOS(env)!=owner||ENCLOS(core)!=owner||!cr_factory_names_ok(env)||
    reader!=VECTOR_ELT(values,6)||CLOENV(reader)!=env)return 0;
  for(int j=0;j<8;++j)if(cr_value(env,cr_factory_names[j])!=VECTOR_ELT(values,j))return 0;
  SEXP read=VECTOR_ELT(values,0),vf=VECTOR_ELT(values,7),helpers=VECTOR_ELT(live,5);
  if(TYPEOF(read)!=CLOSXP||OBJECT(read)||TYPEOF(vf)!=CLOSXP||OBJECT(vf)||CLOENV(read)!=core||CLOENV(vf)!=env||
    !cr_reader_code(BODY_EXPR(reader),VECTOR_ELT(live,3),1)||!cr_same(BODY_EXPR(vf),VECTOR_ELT(live,4))||
    BODY_EXPR(read)!=VECTOR_ELT(live,8)||!cr_same(ATTRIB(reader),VECTOR_ELT(live,9))||
    !cr_same(ATTRIB(vf),VECTOR_ELT(live,10))||!cr_same(ATTRIB(read),VECTOR_ELT(live,11)))return 0;
  if(cr_value(core,"read_cell_score")!=reader||cr_value(core,"complete_cache")!=VECTOR_ELT(live,7)||
    !cr_same(cr_field(VECTOR_ELT(live,7),"view"),view)||!cr_same(cr_field(VECTOR_ELT(live,7),"state"),state))return 0;
  const char *fields[]={"backend","cells","groups","N","B"};const char *statekeys[]={"core","cells","groups","N","B"};
  for(int j=0;j<5;++j)if(cr_value(core,fields[j])!=cr_field(state,statekeys[j]))return 0;
  SEXP recipe=cr_value(core,"recipe"),sets=cr_value(core,"sets");
  if(recipe!=cr_field(cr_field(state,"core"),"outcome_recipe")||!cr_list(sets)||!cr_same(cr_field(recipe,"intercept"),cr_field(state,"intercept"))||
    cr_field(sets,"cids")!=cr_field(state,"cids"))return 0;
  for(int j=0;j<3;++j)if(cr_value(owner,cr_live_helpers[j])!=cr_field(helpers,cr_live_helpers[j]))return 0;
  SEXP prims=VECTOR_ELT(live,6);
  for(int j=0;j<CR_LIVE_PRIMITIVES;++j){SEXP x=cr_inherited_value(owner,cr_live_primitives[j]);
    if(x!=VECTOR_ELT(prims,j))return 0;}
  SEXP fnames=Rf_getAttrib(functions,R_NamesSymbol);
  for(R_xlen_t j=0;j<XLENGTH(functions);++j){const char *name=CHAR(STRING_ELT(fnames,j));
    if(!strcmp(name,".owner")||!strcmp(name,".live"))continue;
    if(cr_bound(owner,Rf_install(name))&&R_BindingIsActive(Rf_install(name),owner))return 0;}
  return cr_functions(functions);
}
static void cr_release(SEXP ptr) {
  coarse53 *p=(coarse53*)R_ExternalPtrAddr(ptr);if(p){R_Free(p);R_ClearExternalPtr(ptr);R_SetExternalPtrProtected(ptr,R_NilValue);}
}
static SEXP cr_source(SEXP ptr,int slot){return VECTOR_ELT(R_ExternalPtrProtected(ptr),slot);}
static double cr_add(double x,double y){volatile double z=x+y;return z;}
static double cr_sub(double x,double y){volatile double z=x-y;return z;}
static double cr_mul(double x,double y){volatile double z=x*y;return z;}
static double cr_div(double x,double y){volatile double z=x/y;return z;}
static int cr_active(SEXP active,int gid){for(R_xlen_t i=0;i<XLENGTH(active);++i)if(INTEGER(active)[i]==gid)return 1;return 0;}
static SEXP cr_at(SEXP call,int k){SEXP p=CDR(call);for(int i=1;i<k;++i)p=CDR(p);return p;}
static void cr_put(SEXP call,int k,SEXP value){SETCAR(cr_at(call,k),value);}

SEXP C_gmm_reconstruction_coarse_create(SEXP view,SEXP state,SEXP centered,SEXP functions,SEXP cap) {
  if(!cr_list(view)||!cr_list(state)||!cr_functions(functions)||
    !((cr_plain(cap,REALSXP)||cr_plain(cap,INTSXP))&&XLENGTH(cap)==1))return R_NilValue;
  double limit=TYPEOF(cap)==REALSXP?REAL(cap)[0]:INTEGER(cap)[0];if(!R_FINITE(limit)||limit<=0)return R_NilValue;
  SEXP cells=cr_field(state,"cells"),groups=cr_field(state,"groups"),ids=cr_field(state,"cids"),
    cache=cr_field(view,"cache"),pos=cr_field(view,"position"),intercept=cr_field(state,"intercept");
  SEXP seen=cr_field(view,"seen");
  int B=0;
  if(!cr_plain(cache,VECSXP)||!cr_plain(seen,LGLSXP)||XLENGTH(seen)!=XLENGTH(cache)||!cr_plain(cells,VECSXP)||!cr_plain(groups,VECSXP)||XLENGTH(cells)>INT_MAX||XLENGTH(groups)>INT_MAX||
    !cr_int(cr_field(state,"B"),1,&B)||!cr_plain(cache,VECSXP)||!cr_plain(pos,INTSXP)||XLENGTH(pos)!=XLENGTH(cells)||
    !cr_ids(ids,(int)XLENGTH(cells),1)||!cr_plain(intercept,LGLSXP)||XLENGTH(intercept)!=1||LOGICAL(intercept)[0]==NA_LOGICAL||
    (centered!=R_NilValue&&(!cr_plain(centered,VECSXP)||XLENGTH(centered)!=XLENGTH(cache))))return R_NilValue;
  SEXP N=cr_field(state,"N");int countN=0;if(!cr_int(N,1,&countN))return R_NilValue;
  /* Canonical complete plan caches exactly its ordered needed cells. This
     proves O(1) visit membership from position, without a new full-cell map. */
  if(XLENGTH(cache)!=XLENGTH(ids))return R_NilValue;
  R_xlen_t positive=0;
  for(R_xlen_t j=0;j<XLENGTH(pos);++j){int slot=INTEGER(pos)[j];if(slot<0||slot>XLENGTH(cache))return R_NilValue;if(slot>0)++positive;}
  if(positive!=XLENGTH(ids))return R_NilValue;
  for(R_xlen_t j=0;j<XLENGTH(ids);++j)if(INTEGER(pos)[INTEGER(ids)[j]-1]!=j+1)return R_NilValue;
  int peak=0;double sum_p=0,incidence=0;
  for(R_xlen_t g=0;g<XLENGTH(groups);++g) {
    SEXP gf=VECTOR_ELT(groups,g);int p=0;
    if(!cr_list(gf)||!cr_int(cr_field(gf,"p"),0,&p)||!cr_ids(cr_field(gf,"cell_ids"),(int)XLENGTH(cells),0)||
      !cdid_reconstruction_materialized(cr_field(gf,"fit")))return R_NilValue;
    SEXP memberids=cr_field(gf,"cell_ids");
    for(R_xlen_t k=0;k<XLENGTH(memberids);++k){SEXP cf=VECTOR_ELT(cells,INTEGER(memberids)[k]-1);int membergroup=0;
      if(!cr_list(cf)||!cr_int(cr_field(cf,"group_id"),1,&membergroup)||membergroup!=g+1)return R_NilValue;}
    if(p>peak)peak=p;sum_p+=p;incidence+=XLENGTH(memberids);
  }
  /* Every fitted cell belongs to one group; a wider private incidence keeps R. */
  if(incidence!=XLENGTH(cells))return R_NilValue;
  for(R_xlen_t c=0;c<XLENGTH(cells);++c) {
    SEXP cf=VECTOR_ELT(cells,c);int gid=0,p=0;
    if(!cr_list(cf)||!cr_int(cr_field(cf,"group_id"),1,&gid)||gid>XLENGTH(groups))return R_NilValue;
    SEXP gf=VECTOR_ELT(groups,gid-1);cr_int(cr_field(gf,"p"),0,&p);
    if(!cr_num(cr_field(cf,"mean_d"),p,cr_field(gf,"kept_columns"))||
      !((cr_plain(cr_field(cf,"mass"),REALSXP)||cr_plain(cr_field(cf,"mass"),INTSXP))&&XLENGTH(cr_field(cf,"mass"))==1))return R_NilValue;
    if(TYPEOF(cr_field(cf,"mass"))==INTSXP&&INTEGER(cr_field(cf,"mass"))[0]==NA_INTEGER)return R_NilValue;
  }
  for(R_xlen_t j=0;j<XLENGTH(ids);++j) {
    int cid=INTEGER(ids)[j],slot=INTEGER(pos)[cid-1];if(slot<1||slot>XLENGTH(cache)||LOGICAL(seen)[slot-1]!=TRUE)return R_NilValue;
    SEXP rec=VECTOR_ELT(cache,slot-1);if(rec==R_NilValue)continue;
    SEXP cf=VECTOR_ELT(cells,cid-1);int gid=0,p=0;cr_int(cr_field(cf,"group_id"),1,&gid);cr_int(cr_field(VECTOR_ELT(groups,gid-1),"p"),0,&p);
    if(!cr_list(rec))return R_NilValue;SEXP rr=cr_field(rec,"rr");
    if(!cr_ids(rr,B,1)||XLENGTH(rr)>B||!cr_plain(cr_field(rec,"s_rr"),REALSXP)||XLENGTH(cr_field(rec,"s_rr"))!=XLENGTH(rr)||
      !cr_matrix(cr_field(rec,"d_rr"),(int)XLENGTH(rr),p))return R_NilValue;
    if(centered!=R_NilValue&&p&&LOGICAL(intercept)[0]&&
      !cr_matrix(VECTOR_ELT(centered,slot-1),(int)XLENGTH(rr),p))return R_NilValue;
  }
  double retained=2048.0,peak_entries=4096.0+4.0*XLENGTH(cells)+128.0*XLENGTH(groups)+4.0*sum_p+
    32.0*peak*peak+(double)B*(4.0*peak+8.0);
  /* Existing forcing+dense solved+first-inactive records and scalar adjoint
     compartments; no new numeric cache. Explicit configured cap unchanged. */
  if(!R_FINITE(peak_entries)||retained+peak_entries>limit||(double)B*peak>R_XLEN_T_MAX)return R_NilValue;
  cdid_pair_operation_arithmetic();
  SEXP anchors=PROTECT(Rf_allocVector(VECSXP,6));SET_VECTOR_ELT(anchors,0,view);SET_VECTOR_ELT(anchors,1,state);
  SET_VECTOR_ELT(anchors,2,centered);SET_VECTOR_ELT(anchors,3,functions);SET_VECTOR_ELT(anchors,4,cap);
  SET_VECTOR_ELT(anchors,5,cr_live_create(view,state,functions));
  SEXP ptr=PROTECT(R_MakeExternalPtr(NULL,Rf_install("cdid_reconstruction_coarse_operation"),anchors));
  R_RegisterCFinalizerEx(ptr,cr_release,TRUE);coarse53 *d=R_Calloc(1,coarse53);d->B=B;d->nc=(int)XLENGTH(cells);
  d->ng=(int)XLENGTH(groups);d->peak=peak;d->Ncount=countN;d->intercept=LOGICAL(intercept)[0];d->total_p=sum_p;d->cap=limit;R_SetExternalPtrAddr(ptr,d);
  UNPROTECT(2);return ptr;
}

SEXP C_gmm_reconstruction_coarse_pullback(SEXP ptr,SEXP view,SEXP state,SEXP lambda,SEXP nu,
                                         SEXP active,SEXP visit,SEXP templates) {
  if(TYPEOF(ptr)!=EXTPTRSXP||R_ExternalPtrTag(ptr)!=Rf_install("cdid_reconstruction_coarse_operation"))return R_NilValue;
  if(view==R_NilValue&&state==R_NilValue&&lambda==R_NilValue&&nu==R_NilValue&&active==R_NilValue&&visit==R_NilValue&&templates==R_NilValue){cr_release(ptr);return R_NilValue;}
  coarse53 *d=(coarse53*)R_ExternalPtrAddr(ptr);if(!d)return R_NilValue;
  if(TYPEOF(view)==CLOSXP&&state==R_NilValue&&nu==R_NilValue&&active==R_NilValue&&templates==R_NilValue){
    SEXP oldview=cr_source(ptr,0),oldstate=cr_source(ptr,1),f=cr_source(ptr,3),live=cr_source(ptr,5);
    if(!cr_plain(visit,INTSXP))return R_NilValue;
    if(!cr_live_ok(live,view,oldview,oldstate,f))return R_NilValue;
    SEXP helper=cr_field(VECTOR_ELT(live,5),"reconstruction_sparse_adjoint_ids");
    SEXP call=PROTECT(Rf_lang3(helper,VECTOR_ELT(live,7),visit));
    SEXP ids=PROTECT(Rf_eval(call,R_BaseEnv));
    if(ids==R_NilValue){UNPROTECT(2);return R_NilValue;}
    if(!cr_live_ok(live,view,oldview,oldstate,f)){UNPROTECT(2);return R_NilValue;}
    SEXP out=PROTECT(C_gmm_reconstruction_score_pullback(oldview,oldstate,lambda,ids,cr_source(ptr,2),cr_source(ptr,4)));
    UNPROTECT(3);return out;
  }
  if(!cr_list(view)||!cr_list(state))return R_NilValue;
  SEXP oldv=cr_source(ptr,0),olds=cr_source(ptr,1),functions=cr_source(ptr,3),
    cells=cr_field(state,"cells"),groups=cr_field(state,"groups"),N=cr_field(state,"N");
  int currentB=0,currentN=0;SEXP hi=cr_field(state,"intercept"),oldids=cr_field(olds,"cids"),currentids=cr_field(state,"cids");
  if(!cr_int(cr_field(state,"B"),1,&currentB)||currentB!=d->B||!cr_int(N,1,&currentN)||currentN!=d->Ncount||
    !cr_plain(hi,LGLSXP)||XLENGTH(hi)!=1||LOGICAL(hi)[0]!=d->intercept||
    !cr_ids(currentids,d->nc,0)||XLENGTH(currentids)!=XLENGTH(oldids))return R_NilValue;
  for(R_xlen_t j=0;j<XLENGTH(oldids);++j)if(INTEGER(currentids)[j]!=INTEGER(oldids)[j])return R_NilValue;
  if(!cr_functions(functions)||cells!=cr_field(olds,"cells")||groups!=cr_field(olds,"groups")||
    cr_field(view,"cache")!=cr_field(oldv,"cache")||cr_field(view,"position")!=cr_field(oldv,"position")||
    cr_field(view,"seen")!=cr_field(oldv,"seen")||
    !R_compute_identical(N,cr_field(olds,"N"),IDENT_USE_CLOENV))return R_NilValue;
  /* Same admitted chunk/source operation, contraction of an already solved
     dense or structural compact record; no new loading/solve arithmetic. */
  if(nu==R_NilValue&&active==R_NilValue&&templates==R_NilValue){
    if(!cr_ids(visit,d->nc,0))return R_NilValue;
    SEXP pos=cr_field(oldv,"position");
    for(R_xlen_t j=0;j<XLENGTH(visit);++j)if(INTEGER(pos)[INTEGER(visit)[j]-1]==0)return R_NilValue;
    return C_gmm_reconstruction_score_pullback(view,state,lambda,visit,cr_source(ptr,2),cr_source(ptr,4));
  }
  if(!cr_plain(lambda,REALSXP)||XLENGTH(lambda)!=d->nc||!cr_plain(nu,VECSXP)||XLENGTH(nu)!=d->nc||
    !cr_ids(active,d->ng,0)||!cr_ids(visit,d->nc,0)||!cr_plain(templates,VECSXP)||XLENGTH(templates)!=d->ng)return R_NilValue;
  SEXP issued_position=cr_field(oldv,"position");
  for(R_xlen_t j=0;j<XLENGTH(visit);++j)if(INTEGER(issued_position)[INTEGER(visit)[j]-1]==0)return R_NilValue;
  for(int cid=0;cid<d->nc;++cid){SEXP cf=VECTOR_ELT(cells,cid);int gid=0,p=0;cr_int(cr_field(cf,"group_id"),1,&gid);
    SEXP gf=VECTOR_ELT(groups,gid-1);cr_int(cr_field(gf,"p"),0,&p);if(!cr_num(VECTOR_ELT(nu,cid),p,cr_field(gf,"kept_columns")))return R_NilValue;}
  for(int g=0;g<d->ng;++g){SEXP x=VECTOR_ELT(templates,g),gf=VECTOR_ELT(groups,g);int p=0;cr_int(cr_field(gf,"p"),0,&p);
    if(x!=R_NilValue&&(!cr_list(x)||!cr_plain(cr_field(x,"lambda"),REALSXP)||XLENGTH(cr_field(x,"lambda"))!=XLENGTH(cr_field(gf,"cell_ids"))||
      !cr_plain(cr_field(x,"nu"),REALSXP)||XLENGTH(cr_field(x,"nu"))!=p))return R_NilValue;
    if(x!=R_NilValue){SEXP oldcells=cr_field(x,"cells"),gids=cr_field(gf,"cell_ids");
      SEXP oldN=cr_field(x,"N");
      if(cr_field(x,"group")!=gf||!cr_plain(oldcells,VECSXP)||XLENGTH(oldcells)!=XLENGTH(gids)||
        !((cr_plain(oldN,INTSXP)||cr_plain(oldN,REALSXP))&&XLENGTH(oldN)==1)||
        !R_compute_identical(oldN,N,IDENT_USE_CLOENV))return R_NilValue;
      for(R_xlen_t k=0;k<XLENGTH(gids);++k)if(VECTOR_ELT(oldcells,k)!=VECTOR_ELT(cells,INTEGER(gids)[k]-1))return R_NilValue;}}

  cdid_pair_operation_arithmetic();
  int np=0;SEXP tl=PROTECT(Rf_allocVector(REALSXP,d->nc));++np;for(int i=0;i<d->nc;++i)REAL(tl)[i]=0;
  SEXP tn=PROTECT(Rf_allocVector(VECSXP,d->ng));++np;
  SEXP work=PROTECT(Rf_shallow_duplicate(templates));++np;
  SEXP gid_s=PROTECT(Rf_allocVector(INTSXP,1)),active_s=PROTECT(Rf_allocVector(LGLSXP,1));np+=2;
  SEXP zero_args=PROTECT(Rf_allocList(6));++np;SEXP zero_call=PROTECT(Rf_lcons(cr_field(functions,"reconstruction_zero_template"),zero_args));++np;
  cr_put(zero_call,1,work);cr_put(zero_call,2,gid_s);cr_put(zero_call,3,active_s);cr_put(zero_call,5,cells);cr_put(zero_call,6,N);
  SEXP solve_call=PROTECT(Rf_lang4(R_NilValue,R_NilValue,R_NilValue,N));++np;SET_TAG(cr_at(solve_call,3),Rf_install("normalizer"));
  SEXP double_call=PROTECT(Rf_lang2(cr_field(functions,"as.double"),R_NilValue));++np;
  int intercept=LOGICAL(cr_field(state,"intercept"))[0],g=1,solves=0,hits=0,suffix=0;
  while(g<=d->ng){if(!(g%256))R_CheckUserInterrupt();SEXP gf=VECTOR_ELT(groups,g-1),ids=cr_field(gf,"cell_ids");int p=0;cr_int(cr_field(gf,"p"),0,&p);
    INTEGER(gid_s)[0]=g;LOGICAL(active_s)[0]=cr_active(active,g);cr_put(zero_call,4,gf);
    SEXP cached=PROTECT(Rf_eval(zero_call,R_BaseEnv));
    if(cached!=R_NilValue){++hits;SEXP lam=cr_field(cached,"lambda"),v=cr_field(cached,"nu");
      for(R_xlen_t k=0;k<XLENGTH(ids);++k)REAL(tl)[INTEGER(ids)[k]-1]=REAL(lam)[k];SET_VECTOR_ELT(tn,g-1,v);UNPROTECT(1);++g;
      /* Same V50 contiguous prefix and count limits; direct final copies
         eliminate its temporary run record and per-group nu allocation. */
      R_xlen_t nc=0,nv=0;SEXP cids=R_NilValue,cl=R_NilValue,cn=R_NilValue;
      if(g<d->ng&&!cr_active(active,g)) while(g<=d->ng&&cdid_reconstruction_zero_hit(work,groups,cells,active,N,g,&cids,&cl,&cn)){
        if(XLENGTH(cids)>d->nc-nc||XLENGTH(cn)>INT_MAX-nv)break;nc+=XLENGTH(cids);nv+=XLENGTH(cn);
        for(R_xlen_t k=0;k<XLENGTH(cids);++k)REAL(tl)[INTEGER(cids)[k]-1]=REAL(cl)[k];SET_VECTOR_ELT(tn,g-1,cn);++suffix;++g;}
      continue;
    }UNPROTECT(1);
    SEXP rhs=PROTECT(Rf_allocVector(REALSXP,p));for(int j=0;j<p;++j)REAL(rhs)[j]=0;
    for(R_xlen_t k=0;k<XLENGTH(ids);++k){int cid=INTEGER(ids)[k]-1;SEXP cf=VECTOR_ELT(cells,cid),mean=cr_field(cf,"mean_d"),nv=VECTOR_ELT(nu,cid);
      double lam=intercept?REAL(lambda)[cid]:0.0;SEXP mass=cr_field(cf,"mass");
      if(intercept)REAL(tl)[cid]=cr_div(lam,TYPEOF(mass)==REALSXP?REAL(mass)[0]:INTEGER(mass)[0]);
      for(int j=0;j<p;++j)REAL(rhs)[j]=cr_sub(cr_add(REAL(rhs)[j],REAL(nv)[j]),cr_mul(lam,REAL(mean)[j]));
      if(p&&Rf_getAttrib(rhs,R_NamesSymbol)==R_NilValue){SEXP nm=Rf_getAttrib(nv,R_NamesSymbol);if(nm==R_NilValue)nm=Rf_getAttrib(mean,R_NamesSymbol);if(nm!=R_NilValue)Rf_setAttrib(rhs,R_NamesSymbol,nm);}}
    SEXP solved=R_NilValue;
    if(p){++solves;SETCAR(solve_call,cr_field(functions,intercept?"or_pooled_wls_gram_solve":"or_wls_gram_solve"));cr_put(solve_call,1,cr_field(gf,"fit"));cr_put(solve_call,2,rhs);
      SEXP answer=PROTECT(Rf_eval(solve_call,R_BaseEnv));SETCADR(double_call,answer);solved=PROTECT(Rf_eval(double_call,R_BaseEnv));
      SET_VECTOR_ELT(tn,g-1,solved);SETCADR(double_call,R_NilValue);cr_put(solve_call,1,R_NilValue);cr_put(solve_call,2,R_NilValue);UNPROTECT(2);
    }else{solved=PROTECT(Rf_allocVector(REALSXP,0));SET_VECTOR_ELT(tn,g-1,solved);UNPROTECT(1);}
    if(!cr_active(active,g)){
      SEXP record=PROTECT(Rf_allocVector(VECSXP,5)),names=PROTECT(Rf_allocVector(STRSXP,5)),lam=PROTECT(Rf_allocVector(REALSXP,XLENGTH(ids))),cs=PROTECT(Rf_allocVector(VECSXP,XLENGTH(ids)));
      const char *labels[]={"lambda","nu","group","cells","N"};for(int j=0;j<5;++j)SET_STRING_ELT(names,j,Rf_mkChar(labels[j]));
      for(R_xlen_t k=0;k<XLENGTH(ids);++k){REAL(lam)[k]=REAL(tl)[INTEGER(ids)[k]-1];SET_VECTOR_ELT(cs,k,VECTOR_ELT(cells,INTEGER(ids)[k]-1));}
      SET_VECTOR_ELT(record,0,lam);SET_VECTOR_ELT(record,1,VECTOR_ELT(tn,g-1));SET_VECTOR_ELT(record,2,gf);SET_VECTOR_ELT(record,3,cs);SET_VECTOR_ELT(record,4,N);Rf_setAttrib(record,R_NamesSymbol,names);SET_VECTOR_ELT(work,g-1,record);UNPROTECT(4);}
    UNPROTECT(1);++g;
  }
  cr_put(zero_call,1,R_NilValue);cr_put(zero_call,4,R_NilValue);
  SEXP adj=PROTECT(Rf_allocVector(VECSXP,2)),an=PROTECT(Rf_allocVector(STRSXP,2));np+=2;
  SET_STRING_ELT(an,0,Rf_mkChar("tilde_lambda"));SET_STRING_ELT(an,1,Rf_mkChar("tilde_nu"));SET_VECTOR_ELT(adj,0,tl);SET_VECTOR_ELT(adj,1,tn);Rf_setAttrib(adj,R_NamesSymbol,an);
  SEXP phi=PROTECT(C_gmm_reconstruction_score_pullback(view,state,adj,visit,cr_source(ptr,2),cr_source(ptr,4)));++np;
  if(phi==R_NilValue)cdid_pair_operation_error("Coarse pullback's admitted score shape changed during execution.","cdid_gmm_bound_error");
  SEXP out=PROTECT(Rf_allocVector(VECSXP,4)),on=PROTECT(Rf_allocVector(STRSXP,4)),counts=PROTECT(Rf_allocVector(INTSXP,3));np+=3;
  const char *labels[]={"phi","templates","adjoint","counts"};for(int j=0;j<4;++j)SET_STRING_ELT(on,j,Rf_mkChar(labels[j]));
  INTEGER(counts)[0]=solves;INTEGER(counts)[1]=hits;INTEGER(counts)[2]=suffix;
  SET_VECTOR_ELT(out,0,phi);SET_VECTOR_ELT(out,1,work);SET_VECTOR_ELT(out,2,adj);SET_VECTOR_ELT(out,3,counts);Rf_setAttrib(out,R_NamesSymbol,on);
  UNPROTECT(np);return out;
}
