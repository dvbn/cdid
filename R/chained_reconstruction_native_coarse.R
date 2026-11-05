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
  names <- c("C_gmm_reconstruction_coarse_create", "C_gmm_reconstruction_coarse_pullback")
  arities <- c(5L, 8L)
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
  function(view, state, centered = NULL) {
    if (!is.null(state$core$reference_ref) || !is.null(state$core$layout) || !is.null(state$core$pair_outcomes)) return(NULL)
    if (!reconstruction_coarse_source(view)) return(NULL)
    owner <- environment(reconstruction_coarse_source)
    current <- vector("list", length(defaults)); names(current) <- names(defaults)
    for (nm in names(defaults)) {
      if ((exists(nm, owner, inherits = FALSE) && bindingIsActive(nm, owner)) ||
          exists(nm, environment(view$read), inherits = FALSE)) return(NULL)
      current[[nm]] <- get(nm, owner, inherits = TRUE)
      if (!identical(current[[nm]], defaults[[nm]])) return(NULL)
    }
    # Lazy-load serialization may duplicate an otherwise identical closure.
    # Bind the ACTUAL current objects just verified, preserving subsequent
    # native pointer checks against any changed namespace binding.
    functions <- c(current, list(.owner = owner))
    rm(current)
    kernels <- reconstruction_coarse_kernels()
    ptr <- .Call(kernels[[1L]], view, state, centered, functions, state$core$control$max_matrix_entries)
    if (is.null(ptr)) return(NULL)
    frame <- environment(view$read); read <- view$read
    augmented <- inherits(state$core, "cdid_augmented_backend")
    core_name <- if (augmented) "augmented_evaluate_if_chunk_core" else "imputation_evaluate_if_chunk_core"
    wrapper_name <- if (augmented) "evaluate_augmented_if_chunk" else "evaluate_imputation_if_chunk"
    core_fun <- get(core_name, owner); wrapper_fun <- get(wrapper_name, owner)
    rm(view, state, centered, functions)
    active <- TRUE
    check <- function(current_read) !bindingIsActive(core_name, owner) && !bindingIsActive(wrapper_name, owner) &&
      identical(environment(current_read), frame) && identical(get(core_name, owner), core_fun) && identical(get(wrapper_name, owner), wrapper_fun)
    pullback <- function(view, state, lambda, nu, groups, ids, templates) {
      if (!active || !check(view$read)) return(NULL)
      .Call(kernels[[2L]], ptr, view, state, lambda, nu, groups, ids, templates)
    }
    abort <- function() {
      active <<- FALSE
      .Call(kernels[[2L]], ptr, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      invisible(NULL)
    }
    finish <- function() { check(read); abort() }
    list(pullback = pullback, finish = finish, abort = abort)
  }
})
for (.name in c("defaults", "helper_names", "primitive_names")) lockBinding(.name, environment(reconstruction_coarse_operation))
rm(.name)
