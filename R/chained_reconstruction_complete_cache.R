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
    }
  }
  position <- integer(length(state$cells)); position[ids] <- seq_along(ids)
  list(view = view, state = state, needed_position = position)
}

# Structural cells plus EVERY member of their active groups, in the original
# needed order. No per-arm R scan over all cells and no thresholding of values.
reconstruction_sparse_adjoint_ids <- function(snapshot, structural_cids) {
  st <- snapshot$state
  if (!is.integer(structural_cids) || !is.null(attributes(structural_cids)) ||
      anyNA(structural_cids) || any(structural_cids < 1L | structural_cids > length(st$cells))) return(NULL)
  if (any(snapshot$needed_position[structural_cids] == 0L)) return(NULL)
  gids <- unique(vapply(st$cells[structural_cids], `[[`, integer(1L), "group_id"))
  closure <- unique(c(structural_cids, unlist(lapply(st$groups[gids], `[[`, "cell_ids"))))
  at <- snapshot$needed_position[closure]
  st$cids[sort(unique(at[at > 0L]))]
}
reconstruction_complete_cache_get <- function(view, cid) view$cache[[view$position[cid]]]
