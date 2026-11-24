# ---------------------------------------------------------------------------
# Owner storage counted once in the pipeline compositions.
#
# When the fitted nuisance registry returns, the pipeline's persistent prefix
# absorbs the backend's owned entries and, on CP04 routes, the reader's owned
# entries (`combined_initial`, run_chained_ipw_pipeline). The prefix already
# contains the frame and panel owner storage (`owners_live`). Every later
# composition that admits a producer, a covariance stream, a retained moment IF
# phase, a `within` phase or an auxiliary replay adds the backend's owned plus
# borrowed entries to a retained value that extends this prefix. The backend's
# borrowed entries are the reader's owned entries, the frame and panel owner
# storage, the design owners and the external prefix X (CP04), or the frame and
# panel owner storage alone (compact raw IPW).
#
# The record below states, once and from verified identities, which of those
# entries the prefix already contains:
#   * the backend's owned entries, always (the prefix adds them explicitly);
#   * the reader's owned entries, when the backend's borrowed entries are
#     exactly the reader, owner, design and external terms;
#   * the frame and panel owner storage, when in addition the backend borrows
#     the pipeline's own frame and panel owners and the prefix's owner term
#     covers that storage.
# The external prefix X stays with cp04_registry_retained_entries(). Design
# owners, which the prefix does not contain, stay charged. Shared objects and
# strings keep their payload convention: only whole duplicated terms go.
# ---------------------------------------------------------------------------

cp04_owner_overlap_record <- function(registry, prefix_entries, frame_owner, panel_owner, owners_in_prefix) {
  raw <- is.null(registry[["cp04"]])
  b <- if (raw) registry[["raw_ipw"]] else registry[["cp04"]]
  if (!is.list(b) || !is.list(b$accounting) || !is.environment(b$ref)) return(NULL)
  prefix_entries <- ipw_owner_int_bound(prefix_entries, "pipeline prefix with backend storage", min_val = 0)
  owners_in_prefix <- ipw_owner_int_bound(owners_in_prefix, "owner storage in the pipeline prefix", min_val = 0)
  acc <- b$accounting
  backend_owned <- ipw_owner_int_bound(acc$owned_entries, "backend owned entries", min_val = 0)
  borrowed <- ipw_owner_int_bound(acc$borrowed_entries, "backend borrowed entries", min_val = 0)
  same_owner <- function(ref, owner) is.environment(ref) && is.list(owner) && identical(ref, owner$ref)
  reader_owned <- 0; owner_entries <- 0; reader_ref <- NULL
  if (raw) {
    o <- ipw_point_view_owner_entries(frame_owner, panel_owner)
    if (same_owner(b$frame_ref, frame_owner) && same_owner(b$panel_ref, panel_owner) &&
        borrowed == o && o <= owners_in_prefix) owner_entries <- o
  } else {
    reader <- registry[["cp04_reader"]]
    if (is.list(reader) && is.list(reader$accounting) && identical(b$reader_ref, reader$ref)) {
      racc <- reader$accounting
      num <- function(x) if (is.null(x)) 0 else as.double(x)
      o <- num(racc$borrowed_owner_entries)
      parts <- c(num(racc$owned_entries), o, num(acc$design_owned_entries), num(acc$design_borrowed_entries),
        num(acc$retained_external_entries))
      if (all(is.finite(parts)) && all(parts >= 0) && borrowed == sum(parts)) {
        reader_owned <- parts[1L]
        reader_ref <- reader$ref
        if (same_owner(reader$frame_ref, frame_owner) && same_owner(reader$panel_ref, panel_owner) &&
            o <= owners_in_prefix) owner_entries <- o
      }
    }
  }
  overlap <- ipw_owner_safe_sum(c(backend_owned, reader_owned, owner_entries), label = "owner storage in the prefix")
  if (overlap > prefix_entries)
    cp04_abort("Owner storage exceeds the pipeline prefix that contains it.", "integrity")
  rec <- ipw_owner_seal(list(kind = "cp04_owner_overlap", backend_ref = b$ref, reader_ref = reader_ref,
    backend_owned_entries = backend_owned, backend_borrowed_entries = borrowed, reader_owned_entries = reader_owned,
    owner_entries = owner_entries, overlap_entries = overlap, prefix_entries = prefix_entries),
    "cp04_owner_overlap")
  # Kept in a locked environment: payload walks of the registry (the
  # post-pipeline diagnostic copy) charge environments nothing, so the record
  # adds no counted entries of its own.
  holder <- new.env(parent = emptyenv())
  assign("record", rec, envir = holder)
  lockEnvironment(holder, bindings = TRUE)
  holder
}

# The sealed record of `registry`, checked against its backend and reader;
# NULL when the registry has none.
cp04_owner_overlap_checked <- function(registry) {
  holder <- registry[["prefix_owner_storage"]]
  if (is.null(holder)) return(NULL)
  rec <- if (is.environment(holder)) get0("record", envir = holder, inherits = FALSE) else NULL
  b <- if (is.null(registry[["cp04"]])) registry[["raw_ipw"]] else registry[["cp04"]]
  if (!ipw_owner_seal_ok(rec, "cp04_owner_overlap") || !identical(rec$kind, "cp04_owner_overlap") ||
      !is.list(b) || !identical(b$ref, rec$backend_ref) ||
      b$accounting$owned_entries != rec$backend_owned_entries ||
      b$accounting$borrowed_entries != rec$backend_borrowed_entries ||
      (!is.null(rec$reader_ref) && (!is.list(registry[["cp04_reader"]]) ||
        !identical(registry[["cp04_reader"]]$ref, rec$reader_ref) ||
        registry[["cp04_reader"]]$accounting$owned_entries != rec$reader_owned_entries)))
    cp04_abort("Owner storage record belongs to a different backend or reader.", "integrity")
  rec
}

# Retained entries for a composition that adds the backend's owned plus
# borrowed entries: the recorded external prefix and the owner storage the
# prefix already contains are removed once. `retained` must extend the prefix.
cp04_composed_retained_entries <- function(registry, retained) {
  rec <- cp04_owner_overlap_checked(registry)
  if (!is.null(rec) && ipw_owner_int_bound(retained, "pipeline retained entries", min_val = 0) < rec$prefix_entries)
    cp04_abort("Pipeline retained reservation does not extend its recorded owner prefix.", "integrity")
  out <- cp04_registry_retained_entries(registry, retained)
  if (is.null(rec)) out else out - rec$overlap_entries
}

# Owner storage the prefix contains, as a number (0 without a record).
cp04_owner_overlap_entries <- function(registry, part = c("overlap", "owner")) {
  part <- match.arg(part)
  rec <- cp04_owner_overlap_checked(registry)
  if (is.null(rec)) 0 else if (identical(part, "owner")) rec$owner_entries else rec$overlap_entries
}

# The returned result as its payload walks saw it before the record existed:
# the diagnostic copy checks charge the result exactly as before.
cp04_without_owner_record <- function(res) {
  if (is.list(res) && is.list(res[["nuisance_registry"]]) &&
      !is.null(res[["nuisance_registry"]][["prefix_owner_storage"]]))
    res[["nuisance_registry"]][["prefix_owner_storage"]] <- NULL
  res
}
