# Original source proof, once per owned complete chunk, never per loading.
reconstruction_coarse_source <- function(view) {
  read <- if (is.list(view) && !is.object(view)) .subset2(view, "read") else NULL
  if (typeof(read) != "closure" || is.object(read)) return(FALSE)
  core <- environment(read)
  if (!is.environment(core) || identical(core, emptyenv())) return(FALSE)
  owner <- parent.env(core)
  if (!isNamespace(owner) || !identical(owner, environment(reconstruction_coarse_source)) || sys.nframe() > 128L) return(FALSE)
  frames <- sys.frames(); parents <- sys.parents()
  for (j in seq_along(frames)) if (identical(frames[[j]], core)) {
    p <- parents[j]
    return(p > 0L && ((identical(sys.function(j), get("imputation_evaluate_if_chunk_core", owner, inherits = FALSE)) &&
      identical(sys.function(p), get("evaluate_imputation_if_chunk", owner, inherits = FALSE))) ||
      (identical(sys.function(j), get("augmented_evaluate_if_chunk_core", owner, inherits = FALSE)) &&
      identical(sys.function(p), get("evaluate_augmented_if_chunk", owner, inherits = FALSE)))))
  }
  FALSE
}
reconstruction_coarse_kernels <- function() {
  gmm_require_native_dll()
  owner <- environment(sys.function()); dll <- getNamespaceInfo(owner, "DLLs")[["cdid"]]
  names <- c("C_gmm_reconstruction_coarse_create", "C_gmm_reconstruction_coarse_pullback",
    "C_gmm_reconstruction_score_pullback")
  arities <- c(5L, 8L, 6L)
  lapply(seq_along(names), function(j) {
    si <- base::getNativeSymbolInfo(names[j], PACKAGE = dll)
    if (!inherits(si, "NativeSymbolInfo") || is.null(si$address) || !identical(as.integer(si$numParameters), arities[j]))
      cdid_gmm_abort("Native coarse-pullback registration/forcing verification failed.", c("cdid_gmm_native_error", "cdid_gmm_accuracy_error"))
    si$address
  })
}
reconstruction_coarse_operation <- local({
  helper_names <- c("reconstruction_zero_template", "reconstruction_zero_template_run",
    "or_pooled_wls_gram_solve", "or_wls_gram_solve", "reconstruction_adjoint_lambda",
    "reconstruction_adjoint_nu", "reconstruction_adjoint_position", "reconstruction_complete_cache_get",
    "reconstruction_complete_cache_centered")
  defaults <- lapply(helper_names, function(nm) get(nm, environment(), inherits = TRUE)); names(defaults) <- helper_names
  primitive_names <- c("as.double", "sweep", "%*%", "as.vector", "numeric", "+", "-", "*", "/", "[", "[[", "[<-", "%in%", "length", "identical")
  defaults <- c(defaults, setNames(lapply(primitive_names, function(nm) get(nm, baseenv())), primitive_names))
  live_names <- c("reconstruction_score_cache_view", "reconstruction_cache_reader_body", "reconstruction_sparse_adjoint_ids")
  live_defaults <- lapply(live_names, function(nm) get(nm, environment(), inherits = TRUE)); names(live_defaults) <- live_names
  function(view, state, centered = NULL) {
    if (!is.null(state$core$reference_ref) || !is.null(state$core$layout) || !is.null(state$core$pair_outcomes)) return(NULL)
    # Production eligibility is the actual ordinary backend and completed
    # score witness. Private call stacks, bodies and namespace mutations are
    # development concerns, not per-comparison input/model checks.
    if (!is.list(view) || is.object(view) || typeof(view$read) != "closure" ||
        !is.environment(environment(view$read))) return(NULL)
    frame <- environment(view$read)
    if (!exists("complete_cache", frame, inherits = FALSE) ||
        bindingIsActive("complete_cache", frame)) return(NULL)
    complete <- get("complete_cache", frame, inherits = FALSE)
    if (!is.list(complete) || is.object(complete) ||
        !identical(complete$view, view) || !identical(complete$state, state)) return(NULL)
    owner <- environment(reconstruction_coarse_kernels)
    functions <- lapply(names(defaults), function(nm) get(nm, owner, inherits = TRUE))
    names(functions) <- names(defaults)
    functions$.owner <- owner
    kernels <- reconstruction_coarse_kernels()
    ptr <- .Call(kernels[[1L]], view, state, centered, functions, state$core$control$max_matrix_entries)
    if (is.null(ptr)) return(NULL)
    issued_view <- view; issued_state <- state; issued_centered <- centered
    issued_reader <- view$reader; cap <- state$core$control$max_matrix_entries
    rm(view, state, centered, functions)
    active <- TRUE
    pullback <- function(view, state, lambda, nu, groups, ids, templates) {
      if (!active) return(NULL)
      .Call(kernels[[2L]], ptr, view, state, lambda, nu, groups, ids, templates)
    }
    contract <- function(view, state, adjoint, ids) {
      if (!active || !identical(view, issued_view) || !identical(state, issued_state)) return(NULL)
      .Call(kernels[[3L]], view, state, adjoint, ids, issued_centered, cap)
    }
    contract_current <- function(current_reader, adjoint, structural_cids) {
      if (!active || !identical(current_reader, issued_reader)) return(NULL)
      ids <- reconstruction_sparse_adjoint_ids(complete, structural_cids)
      if (is.null(ids)) return(NULL)
      .Call(kernels[[3L]], issued_view, issued_state, adjoint, ids, issued_centered, cap)
    }
    abort <- function() {
      active <<- FALSE
      .Call(kernels[[2L]], ptr, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      invisible(NULL)
    }
    finish <- function() abort()
    list(pullback = pullback, contract = contract, contract_current = contract_current, finish = finish, abort = abort)
  }
})
for (.name in c("defaults", "helper_names", "primitive_names", "live_names", "live_defaults")) lockBinding(.name, environment(reconstruction_coarse_operation))
rm(.name)
