# R --with-keep.source annotates exactly the outer cached-reader block and
# its !seen block. Strip ONLY those three standard source annotations at
# these two guarded nodes; all executable nodes/tags/other attributes remain.
# This is fixed-size work, not a recursive AST walk on every replay arm.
reconstruction_cache_reader_body <- function(x) {
  if (!is.call(x) || length(x) != 5L || !identical(x[[1L]], as.name("{")) ||
      !is.call(x[[4L]]) || length(x[[4L]]) != 3L || !identical(x[[4L]][[1L]], as.name("if")) ||
      !is.call(x[[4L]][[3L]]) || length(x[[4L]][[3L]]) != 4L ||
      !identical(x[[4L]][[3L]][[1L]], as.name("{"))) return(NULL)
  inner <- x[[4L]][[3L]]
  for (nm in c("srcref", "srcfile", "wholeSrcref")) {
    attr(x, nm) <- NULL
    attr(inner, nm) <- NULL
  }
  branch <- x[[4L]]; branch[[3L]] <- inner; x[[4L]] <- branch
  x
}

# One bounded cache-state view, tied to the actual factory reader. No caller
# assertion of completeness is accepted, and no raw reader is invoked here.
reconstruction_score_cache_view <- function(reader) {
  if (typeof(reader) != "closure" || is.object(reader)) return(NULL)
  e <- environment(reader)
  if (!is.environment(e) || !identical(parent.env(e), environment(reconstruction_cached_cell_reader))) return(NULL)
  bindings <- c("read", "plan", "n_cells", "position", "cache", "seen", "reader", "view")
  if (!identical(sort(ls(e, all.names = TRUE)), sort(bindings))) return(NULL)
  for (nm in bindings) if (bindingIsActive(nm, e)) return(NULL)
  view <- get("view", e, inherits = FALSE)
  if (!identical(reader, get("reader", e, inherits = FALSE)) || typeof(view) != "closure" ||
      !identical(environment(view), e) || !identical(attr(reader, "reconstruction_cache_view", exact = TRUE), view) ||
      !identical(body(view), quote(list(position = position, seen = seen, cache = cache, read = read)))) return(NULL)
  expected <- quote({
    slot <- position[cid]
    if (slot == 0L) return(read(cid))
    if (!seen[slot]) {
      record <- read(cid)
      cache[slot] <<- list(record)
      seen[slot] <<- TRUE
    }
    cache[[slot]]
  })
  actual_body <- reconstruction_cache_reader_body(body(reader))
  expected_body <- reconstruction_cache_reader_body(expected)
  if (is.null(actual_body) || is.null(expected_body) || !identical(actual_body, expected_body)) return(NULL)
  x <- view()
  if (!is.integer(x$position) || !is.null(attributes(x$position)) ||
      !is.logical(x$seen) || !is.null(attributes(x$seen)) ||
      !is.list(x$cache) || !is.null(attributes(x$cache)) || length(x$cache) != length(x$seen)) return(NULL)
  c(x, list(reader = reader))
}

# This check is called only AFTER the first unchanged complete traversal. Every
# needed record must now be an actual seen slot, so the validation pass cannot
# invoke read or another user callback. Validate final values, not stale early
# witnesses: a later first-pass read may have changed an earlier record/mean.
reconstruction_complete_cache_witness <- function(view, state) {
  if (is.null(view) || !is.list(state$cells) || is.object(state$cells) ||
      !is.list(state$groups) || is.object(state$groups) ||
      !is.integer(state$cids) || !is.null(attributes(state$cids)) || !length(state$cids) ||
      !identical(state$intercept, isTRUE(state$intercept)) ||
      !is.numeric(state$N) || is.object(state$N) || !is.null(attributes(state$N)) ||
      length(state$N) != 1L || !is.finite(state$N) || state$N <= 0 ||
      !is.numeric(state$B) || is.object(state$B) || !is.null(attributes(state$B)) ||
      length(state$B) != 1L || !is.finite(state$B) || state$B < 1 ||
      length(view$position) != length(state$cells)) return(NULL)
  ids <- state$cids
  if (anyNA(ids) || any(ids < 1L | ids > length(state$cells)) || anyDuplicated(ids)) return(NULL)
  slots <- view$position[ids]
  if (anyNA(slots) || any(slots < 1L | slots > length(view$cache)) || anyDuplicated(slots) ||
      anyNA(view$seen[slots]) || !all(view$seen[slots])) return(NULL)
  plain <- function(x) is.numeric(x) && !is.object(x)
  checked <- logical(length(state$groups)); membership <- integer(length(state$cells))
  # Pure ordinary callers keep all existing per-arm source/state checks. The
  # exact centered value is optional, never a new validity permission. Its
  # list/seal headers fit the existing 1024-slot cache header allowance; each
  # record must also fund source + derived + both construction/replacement
  # overlaps in its unchanged 4*B*(p+2)+256+8*column_bytes reservation.
  core_obj <- if (is.list(state) && !is.object(state)) .subset2(state, "core") else NULL
  retain_centered <- is.list(view) && !is.object(view) && is.list(core_obj) && inherits(core_obj, "cdid_imputation_backend") &&
    is.null(.subset2(core_obj, "reference_ref")) && is.null(.subset2(core_obj, "layout")) &&
    is.null(.subset2(core_obj, "pair_outcomes")) && state$intercept && length(ids) == length(view$cache) &&
    reconstruction_centered_cache_source(view)
  centered_values <- if (retain_centered) vector("list", length(view$cache)) else NULL
  any_centered <- FALSE
  for (j in seq_along(ids)) {
    cf <- state$cells[[ids[j]]]
    if (!is.list(cf) || is.object(cf) || length(cf$group_id) != 1L ||
        !is.numeric(cf$group_id) || is.object(cf$group_id) || !is.null(attributes(cf$group_id)) ||
        !is.finite(cf$group_id) || cf$group_id != floor(cf$group_id) ||
        cf$group_id < 1L || cf$group_id > length(state$groups)) return(NULL)
    gf <- state$groups[[cf$group_id]]
    if (!is.list(gf) || is.object(gf) || !is.numeric(gf$p) || is.object(gf$p) ||
        !is.null(attributes(gf$p)) || length(gf$p) != 1L ||
        !is.finite(gf$p) || gf$p < 0 || gf$p != floor(gf$p) ||
        !is.integer(gf$cell_ids) || !is.null(attributes(gf$cell_ids))) return(NULL)
    if (!checked[cf$group_id]) {
      gcids <- gf$cell_ids
      if (anyNA(gcids) || any(gcids < 1L | gcids > length(state$cells)) ||
          anyDuplicated(gcids) || any(membership[gcids] != 0L)) return(NULL)
      membership[gcids] <- cf$group_id; checked[cf$group_id] <- TRUE
    }
    if (membership[ids[j]] != cf$group_id) return(NULL)
    rec <- view$cache[[slots[j]]]
    if (is.null(rec)) next
    if (!is.list(rec) || is.object(rec) || !identical(names(rec), c("rr", "d_rr", "s_rr")) ||
        !is.integer(rec$rr) || !is.null(attributes(rec$rr)) || anyNA(rec$rr) ||
        any(rec$rr < 1L | rec$rr > state$B) || anyDuplicated(rec$rr) ||
        !plain(rec$s_rr) || !is.null(dim(rec$s_rr)) || length(rec$s_rr) != length(rec$rr) ||
        any(!names(attributes(rec$s_rr)) %in% "names") ||
        !plain(rec$d_rr) || !is.matrix(rec$d_rr) || any(!names(attributes(rec$d_rr)) %in% c("dim", "dimnames")) ||
        nrow(rec$d_rr) != length(rec$rr) ||
        ncol(rec$d_rr) != gf$p || !all(is.finite(rec$s_rr))) return(NULL)
    if (gf$p) {
      if (!plain(cf$mean_d) || length(cf$mean_d) != gf$p || !is.null(dim(cf$mean_d))) return(NULL)
      centered <- if (state$intercept) sweep(rec$d_rr, 2L, cf$mean_d, FUN = "-") else rec$d_rr
      if (!all(is.finite(centered))) return(NULL)
      if (retain_centered) {
        columns <- gf$kept_columns
        peak <- reconstruction_centered_record_peak(rec, centered, columns)
        if (is.finite(peak)) {
          slab <- 4 * state$B * (gf$p + 2) + 256 + 8 * sum(nchar(columns, type = "bytes"))
          # Count the original record and centered matrix twice: retained
          # values plus the current raw/centered replacement construction.
          # One additional list slot per record is explicit, not free aliasing.
          if (peak > slab) {
            retain_centered <- FALSE; centered_values <- NULL
          } else {
            centered_values[slots[j]] <- list(centered)
            any_centered <- TRUE
          }
        } else {
          retain_centered <- FALSE; centered_values <- NULL
        }
      }
    }
  }
  position <- integer(length(state$cells)); position[ids] <- seq_along(ids)
  out <- list(view = view, state = state, needed_position = position)
  if (retain_centered && any_centered) {
    derived <- list(values = centered_values, ref = NULL)
    # The empty-list skeleton excludes the separately funded per-record list
    # slots and matrices. Failed optional funding retains the old score cache.
    # The named two-field empty-list skeleton is 42 entries; another 64
    # covers output/list/seal replacement headers. Per-record slots are in
    # peak above. All 106 entries fit the unchanged 1024-slot header reserve.
    out$centered <- ipw_owner_seal(derived, "reconstruction_centered_cache")
  }
  out
}

# The caller has just passed the existing complete view/state equality. A
# changed derived component declines reuse; the old sweep remains available.
reconstruction_complete_cache_centered <- function(snapshot) {
  x <- snapshot$centered
  ref <- if (is.list(x)) x$ref else NULL
  if (!is.environment(ref) || !identical(parent.env(ref), emptyenv()) || !environmentIsLocked(ref) ||
      !exists("reconstruction_centered_cache", ref, inherits = FALSE) ||
      !bindingIsLocked("reconstruction_centered_cache", ref) || bindingIsActive("reconstruction_centered_cache", ref) ||
      !identical(x, get("reconstruction_centered_cache", ref, inherits = FALSE)) ||
      !is.list(x$values) || length(x$values) != length(snapshot$view$cache)) return(NULL)
  x$values
}

# A bounded existing-frame guard for this numerical-value reuse only. It
# grants no validation permission: every factory recognizer/state check stays.
# Cloned private cores and arbitrary readers retain their original sweep.
reconstruction_centered_cache_source <- function(view) {
  if (!is.list(view) || is.object(view)) return(FALSE)
  read <- .subset2(view, "read")
  if (typeof(read) != "closure" || is.object(read)) return(FALSE)
  core <- environment(read)
  if (!is.environment(core) || identical(core, emptyenv())) return(FALSE)
  owner <- parent.env(core)
  if (!isNamespace(owner) || !identical(owner, environment(reconstruction_complete_cache_witness)) ||
      exists("sweep", core, inherits = FALSE) ||
      (exists("sweep", owner, inherits = FALSE) && bindingIsActive("sweep", owner)) ||
      !identical(get("sweep", owner, inherits = TRUE), base::sweep)) return(FALSE)
  # Inspect existing call-frame identities, never force a private reader
  # promise. The bounded pointer/index lists fit the old cache header reserve.
  if (sys.nframe() > 128L) return(FALSE)
  frames <- sys.frames(); parents <- sys.parents()
  for (j in seq_along(frames)) if (identical(frames[[j]], core)) {
    p <- parents[j]
    return(p > 0L && identical(sys.function(j), get("imputation_evaluate_if_chunk_core", owner, inherits = FALSE)) &&
      identical(sys.function(p), get("evaluate_imputation_if_chunk", owner, inherits = FALSE)))
  }
  FALSE
}

# Exact cp04_payload counts for canonical score records and their already
# computed centered matrix. Unusual private decorations decline retention.
# Matrix: n*p + p + 46 + sum(ceil(column_bytes/2)); record adds 2*n+52.
reconstruction_centered_record_peak <- function(record, centered, columns) {
  if (!is.character(columns) || is.object(columns) || anyNA(columns) || !is.null(attributes(columns)) ||
      !identical(names(attributes(record)), "names") || !is.null(attributes(record$s_rr)) ||
      !identical(names(attributes(record$d_rr)), c("dim", "dimnames")) ||
      !is.null(rownames(record$d_rr)) || !identical(colnames(record$d_rr), columns) ||
      !identical(attributes(centered), attributes(record$d_rr))) return(Inf)
  p <- length(columns); n <- length(record$rr)
  matrix_entries <- n * p + p + 46 + sum(ceiling(nchar(columns, type = "bytes") / 2))
  record_entries <- matrix_entries + 2 * n + 52
  2 * record_entries + 2 * matrix_entries + 1
}

# Structural cells plus EVERY member of their active groups, in the original
# needed order. No per-arm R scan over all cells and no thresholding of values.
reconstruction_sparse_adjoint_ids <- function(snapshot, structural_cids) {
  st <- snapshot$state
  if (!is.integer(structural_cids) || !is.null(attributes(structural_cids)) ||
      anyNA(structural_cids) || any(structural_cids < 1L | structural_cids > length(st$cells))) return(NULL)
  if (any(snapshot$needed_position[structural_cids] == 0L)) return(NULL)
  # The complete witness validates every needed cell's group and membership
  # before its p>0 branch. With one group, any nonempty valid arm's closure
  # therefore contains all needed cells, already in their required order.
  # Empty private arms retain their original empty result below.
  if (length(structural_cids) && length(st$groups) == 1L) return(st$cids)
  gids <- unique(vapply(st$cells[structural_cids], `[[`, integer(1L), "group_id"))
  closure <- unique(c(structural_cids, unlist(lapply(st$groups[gids], `[[`, "cell_ids"))))
  at <- snapshot$needed_position[closure]
  st$cids[sort(unique(at[at > 0L]))]
}
reconstruction_complete_cache_get <- function(view, cid) view$cache[[view$position[cid]]]
