# Representation only: callers have already performed EVERY scalar group
# solve and its original checks. Support is structural, never a zero test.
reconstruction_pack_adjoint <- function(tilde_lambda, tilde_nu, cell_ids, cells) {
  dense <- list(tilde_lambda = tilde_lambda, tilde_nu = tilde_nu)
  cids <- sort(unique(as.integer(cell_ids)))
  gids <- sort(unique(vapply(cells[cids], `[[`, integer(1L), "group_id")))
  compact <- list(kind = "compact_reconstruction_adjoint", cell_ids = cids,
    tilde_lambda = tilde_lambda[cids], group_ids = gids, tilde_nu = unname(tilde_nu[gids]))
  # Tiny catalogs may need fewer headers in their original representation.
  # Exact payload comparison changes storage only, never a fitted coefficient.
  if (unname(cp04_payload(compact)["entries"]) < unname(cp04_payload(dense)["entries"])) compact else dense
}

# Only the fitted constructors call this: lambda/list are unnamed, and every
# nu is as.double(scalar solve) or numeric(0). No decorated caller shortcut.
reconstruction_pack_plain_adjoint <- function(tilde_lambda, tilde_nu, cell_ids, cells) {
  dense <- list(tilde_lambda = tilde_lambda, tilde_nu = tilde_nu)
  entries <- length(tilde_lambda) + 9 * length(tilde_nu) + sum(lengths(tilde_nu)) + 55
  reconstruction_choose_adjoint(dense, cell_ids, cells, entries)
}

reconstruction_choose_adjoint <- function(dense, cell_ids, cells, dense_entries) {
  tilde_lambda <- dense$tilde_lambda; tilde_nu <- dense$tilde_nu
  cids <- sort(unique(as.integer(cell_ids)))
  gids <- sort(unique(vapply(cells[cids], `[[`, integer(1L), "group_id")))
  compact <- list(kind = "compact_reconstruction_adjoint", cell_ids = cids,
    tilde_lambda = tilde_lambda[cids], group_ids = gids, tilde_nu = unname(tilde_nu[gids]))
  # Tiny catalogs may need fewer headers in their original representation.
  # Exact payload comparison changes storage only, never a fitted coefficient.
  if (unname(cp04_payload(compact)["entries"]) < dense_entries) compact else dense
}

# Sorted structural integer IDs permit bounded lookup without a dense map or
# per-call hash table. At most31 comparisons for admitted integer catalogs.
reconstruction_adjoint_position <- function(ids, id) {
  lo <- 1; hi <- as.double(length(ids))
  while (lo <= hi) {
    mid <- floor((lo + hi) / 2)
    if (ids[mid] == id) return(mid)
    if (ids[mid] < id) lo <- mid + 1 else hi <- mid - 1
  }
  0L
}
reconstruction_adjoint_lambda <- function(adj, cid) {
  if (!identical(adj$kind, "compact_reconstruction_adjoint")) return(adj$tilde_lambda[cid])
  pos <- reconstruction_adjoint_position(adj$cell_ids, cid)
  if (pos) adj$tilde_lambda[pos] else 0
}
reconstruction_adjoint_nu <- function(adj, gid, p) {
  if (!identical(adj$kind, "compact_reconstruction_adjoint")) return(adj$tilde_nu[[gid]])
  pos <- reconstruction_adjoint_position(adj$group_ids, gid)
  # A missing group still has its original p coordinates and participates in
  # the unchanged score multiplication; zero times overflow is not skipped.
  if (pos) adj$tilde_nu[[pos]] else numeric(p)
}

# Test/reference adapter for the original dense record; production replay uses
# scalar lookup and never allocates a cell-by-link or group-by-link table.
reconstruction_expand_adjoint <- function(adj, n_cells, groups) {
  if (!identical(adj$kind, "compact_reconstruction_adjoint")) return(adj)
  lambda <- numeric(n_cells); lambda[adj$cell_ids] <- adj$tilde_lambda
  nu <- lapply(groups, function(g) numeric(g$p))
  nu[adj$group_ids] <- adj$tilde_nu
  list(tilde_lambda = lambda, tilde_nu = nu)
}

# Dense input loadings, dense solved record, compact candidate, scalar RHS/
# factors and replacement/list/header copies coexist during packing.
# ng,k<=nc makes this no larger than the existing128nc(p+2)^2+8192 slab.
reconstruction_adjoint_construction_scratch <- function(nc, ng, p, k) {
  ipw_owner_safe_sum(c(8192, ipw_owner_safe_product(4, nc),
    ipw_owner_safe_product(ng, 8 * p + 128), ipw_owner_safe_product(8,
    ipw_owner_safe_sum(c(ipw_owner_safe_product(nc, p + 12),
      ipw_owner_safe_product(ng, p + 12), ipw_owner_safe_product(k, p + 12),
      ipw_owner_safe_product(p, p))))))
}

# A template is populated only AFTER the group's original inactive scalar
# arithmetic/solve succeeds. Active groups never reuse it, even for zero data.
reconstruction_zero_template <- function(templates, gid, active, group, cells, N) {
  if (active) return(NULL)
  x <- templates[[gid]]
  if (is.null(x) || !identical(group, x$group) ||
      !identical(cells[group$cell_ids], x$cells) || !identical(N, x$N)) return(NULL)
  x
}
