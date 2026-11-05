# Coarse ordinary legacy CP04 replay. Borrow the existing admitted lazy pair
# operation, with one package-helper/source admission before numerical work.
# General selection models, private/reference/layout outcomes and details keep
# their existing evaluator. Original R logit, clipping and solve helpers remain.
cp04_legacy_ipw_native <- local({
  helpers <- list(cp04_predict_selection = cp04_predict_selection,
    cp04_solve = cp04_solve, cp04_finite = cp04_finite,
    ipw_logit_values = ipw_logit_values, ipw_clip_logit = ipw_clip_logit,
    ipw_clip_floor = ipw_clip_floor, cp04_abort = cp04_abort)
  constants <- list(ipw_clip_eps = ipw_clip_eps, ipw_clip_cap = ipw_clip_cap,
    ipw_clip_cap_odds = ipw_clip_cap_odds)
  run <- function(backend, reader, frame, rows, tids, sids, pair_op, funding_entries) {
    if (is.null(pair_op) || !identical(class(backend), "cdid_cp04_corollary1") ||
        !identical(class(reader), "cdid_cp04_frame_reader") ||
        !identical(class(frame), c("cdid_ipw_frame_owner", "cdid_ipw_context")) || isTRUE(backend$uses_selection_models) ||
        !is.null(backend$reference_ref) || !is.null(backend$layout) ||
        !is.null(backend$pair_outcomes)) return(NULL)
    owner <- environment(helpers$cp04_solve)
    if (!isNamespace(owner) || !all(vapply(c(names(helpers), names(constants)), function(name) {
      expected <- if (name %in% names(helpers)) helpers[[name]] else constants[[name]]
      exists(name, owner, inherits = FALSE) && bindingIsLocked(name, owner) &&
        !bindingIsActive(name, owner) && identical(get(name, owner, inherits = FALSE), expected)
    }, logical(1L)))) return(NULL)
    syms <- gmm_require_native_dll()
    dll <- getNamespaceInfo(owner, "DLLs")[["cdid"]]
    si <- base::getNativeSymbolInfo("C_gmm_legacy_ipw_chunk", PACKAGE = dll)
    if (!inherits(si, "NativeSymbolInfo") || is.null(si$address) ||
        !identical(as.integer(si$numParameters), 9L))
      cdid_gmm_abort("Native symbol registration/forcing verification failed for `C_gmm_legacy_ipw_chunk`.",
        c("cdid_gmm_native_error", "cdid_gmm_accuracy_error"))
    attr(syms, "calls")$C_gmm_legacy_ipw_chunk(reader, frame, backend, rows, tids, sids, pair_op$pair, helpers, funding_entries)
  }
  lockEnvironment(environment(), bindings = TRUE)
  run
})
