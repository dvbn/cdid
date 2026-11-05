# Shared ordinary, callback-free sparse POSITION reads. Source/helper/DLL
# permission belongs to this issued chunk only; standalone reads stay intact.
cp04_pair_operation_kernels <- function() {
  owner <- environment(cp04_pair_rows)
  if (!isNamespace(owner)) return(NULL)
  dll <- getNamespaceInfo(owner, "DLLs")[["cdid"]]
  arities <- c(C_gmm_pair_operation_create = 9L, C_gmm_pair_operation_query = 6L)
  lapply(names(arities), function(name) {
    si <- base::getNativeSymbolInfo(name, PACKAGE = dll)
    if (!inherits(si, "NativeSymbolInfo") || is.null(si$address) ||
        !identical(as.integer(si$numParameters), arities[[name]]))
      cdid_gmm_abort(paste0("Native symbol registration/forcing verification failed for `", name, "`."),
        c("cdid_gmm_native_error", "cdid_gmm_accuracy_error"))
    si$address
  })
}

cp04_pair_read_operation <- local({
  helpers <- list(cp04_pair_rows = cp04_pair_rows, cp04_pair_rows_native_ok = cp04_pair_rows_native_ok,
    cp04_pair_rows_native_kernel = cp04_pair_rows_native_kernel,
    ipw_owner_read_y_positions = ipw_owner_read_y_positions,
    ipw_owner_match_sorted_double = ipw_owner_match_sorted_double, cp04_finite = cp04_finite,
    cp04_pair_operation_kernels = cp04_pair_operation_kernels)
  issue <- function(reader, rows, funding_entries) {
    force(reader); force(rows); force(funding_entries)
    if (!is.numeric(rows) || is.object(rows) || !is.null(attributes(rows)) || !length(rows) ||
        !gmm_scalar_numeric(funding_entries) || !is.finite(funding_entries) || funding_entries <= 0) return(NULL)
    panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
    B <- length(rows); T <- length(panel$postings$date_start)
    reservation <- 4096 + 64 * B + 8 * T + 3 * B * T
    if (!is.finite(reservation) || reservation > funding_entries) return(NULL)
    # Package seq.int batches may be compact ALTREP indices. This bounded
    # owned copy precedes all source/DLL permission, never coerces source Y or
    # posting arrays, and its copy/replacement peak is charged above.
    owned_rows <- integer(B); owned_rows[] <- rows
    validate_cp04_reader(reader)
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
    helper_check <- function() {
      owner <- environment(helpers$cp04_pair_rows)
      isNamespace(owner) && all(vapply(names(helpers), function(name)
        exists(name, owner, inherits = FALSE) && bindingIsLocked(name, owner) && !bindingIsActive(name, owner) &&
          identical(get(name, envir = owner, inherits = FALSE), helpers[[name]]), logical(1L)))
    }
    if (!helper_check()) return(NULL)
    syms <- gmm_require_native_dll()
    kernels <- cp04_pair_operation_kernels()
    if (is.null(kernels)) return(NULL)
    post <- panel$postings
    state <- .Call(kernels[[1L]], owned_rows, reader$state$a_index,
      post$post_a_index, post$post_obs_index, post$date_start, post$date_end,
      post$date_count, panel$observations$y, reservation)
    if (is.null(state)) return(NULL)
    active <- TRUE; final_counts <- NULL; fallbacks <- 0L
    index_fields <- c("s", "t", "s_calendar_index", "t_calendar_index")
    owner_fields <- c(index_fields, "time_gap", "time_gap_finite")
    catalog_fields <- c(owner_fields, "L", "pair_label", "pair_id")
    index_ok <- function(x) is.numeric(x) && !is.object(x) && is.null(attributes(x)) &&
      length(x) == 1L && is.finite(x) && x == floor(x) && x >= 1 && x <= length(post$date_start)
    pair <- function(r_sub, descriptor, read_outcomes = FALSE) {
      if (!active) cp04_abort("Sparse pair operation is closed.", "integrity")
      plain_pair <- is.list(descriptor) && !is.object(descriptor) &&
        identical(names(attributes(descriptor)), "names") &&
        (identical(names(descriptor), index_fields) || identical(names(descriptor), owner_fields) ||
          identical(names(descriptor), catalog_fields))
      legacy_pair <- is.list(descriptor) && identical(class(descriptor), "cdid_ipw_pair_descriptor") &&
        identical(sort(names(attributes(descriptor))), c("class", "names")) &&
        identical(names(descriptor), index_fields)
      if (!plain_pair && !legacy_pair) {
        fallbacks <<- fallbacks + 1L
        return(cp04_pair_rows(reader, r_sub, descriptor, read_outcomes))
      }
      ds <- descriptor$s_calendar_index; dt <- descriptor$t_calendar_index
      # Owner sorted-match positions are ordinary numeric scalars. Convert
      # only exact, in-range calendar INDICES, never a numerical source array.
      if (!index_ok(ds) || !index_ok(dt)) {
        fallbacks <<- fallbacks + 1L
        return(cp04_pair_rows(reader, r_sub, descriptor, read_outcomes))
      }
      hit <- .Call(kernels[[2L]], state, r_sub, as.integer(c(ds, dt)), read_outcomes, reservation, 0L)
      if (is.null(hit)) {
        fallbacks <<- fallbacks + 1L
        return(cp04_pair_rows(reader, r_sub, descriptor, read_outcomes))
      }
      n <- length(r_sub); out <- list(S = numeric(n))
      if (read_outcomes) out$dy <- numeric(n)
      obs <- hit$obs
      if (length(obs)) {
        out$S[obs] <- 1
        if (read_outcomes) {
          yy <- ipw_owner_read_y_positions(panel, hit$pos_s, hit$pos_t)
          out$dy[obs] <- yy$y_right_obs - yy$y_left_obs
        }
      }
      if (read_outcomes) cp04_finite(out$dy, "observed pair differences")
      out
    }
    finish <- function() {
      if (!active) cp04_abort("Sparse pair operation is closed.", "integrity")
      current <- gmm_require_native_dll()
      if (!isTRUE(.Call(current$C_gmm_sexp_aliases, syms, current)) || !helper_check())
        cp04_abort("Sparse pair operation DLL or private reader helper changed.", "integrity")
      for (source in list(reader, frame, panel)) {
        binding <- if (inherits(source, "cdid_cp04_frame_reader")) "cp04_reader" else
          if (inherits(source, "cdid_ipw_frame_owner")) "frame_owner" else "panel_owner"
        ref <- source$ref
        if (!is.environment(ref) || !environmentIsLocked(ref) || !exists(binding, ref, inherits = FALSE) ||
            !bindingIsLocked(binding, ref) || bindingIsActive(binding, ref) ||
            !isTRUE(.Call(current$C_gmm_sexp_aliases, source, get(binding, ref, inherits = FALSE))))
          cp04_abort("Sparse pair operation source changed before its successful boundary.", "integrity")
      }
      final_counts <<- c(.Call(kernels[[2L]], state, NULL, NULL, FALSE, reservation, 1L), fallbacks)
      names(final_counts) <<- c("searches", "cache_hits", "queries", "fallbacks")
      active <<- FALSE
      invisible(NULL)
    }
    abort <- function() {
      active <<- FALSE
      .Call(kernels[[2L]], state, NULL, NULL, FALSE, reservation, 2L)
      invisible(NULL)
    }
    counts <- function() {
      if (!is.null(final_counts)) return(final_counts)
      if (!active) return(NULL)
      out <- c(.Call(kernels[[2L]], state, NULL, NULL, FALSE, reservation, 1L), fallbacks)
      stats::setNames(out, c("searches", "cache_hits", "queries", "fallbacks"))
    }
    list(pair = pair, finish = finish, abort = abort, counts = counts,
      reservation = reservation, funding_entries = funding_entries, chunk_rows = B, n_dates = T)
  }
  lockEnvironment(environment(), bindings = TRUE)
  issue
})
