# Private local stage-2 outcome access. Stage-1 group fits and treatment roots
# keep their own identities and dimensions; no fitted-weight deduplication.
reconstruction_separate <- function(sh) !is.null(sh) && identical(sh$outcome_pooling, "separate")

reconstruction_or_block <- function(sh, bid) {
  if (!reconstruction_separate(sh)) return(sh)
  f <- or_outcome_block_fit(sh$outcome_groups, bid)
  list(keep = f$keep, p = f$p, kept_columns = f$basis$kept_columns,
    kept_cols_F = sh$slope_cols_F[f$keep], beta = f$beta, fit = f$fit,
    summary = sh$outcome_groups$groups[[f$group_id]]$summary, basis = f$basis,
    n_blocks = sh$n_blocks)
}

# Replace only the pooled stage-2 fields in the existing compact payload model.
# Stage1, treatment and link adjoints remain their actual distinct allocations.
reconstruction_retained_bound <- function(cat, recipe, pooled, theta_ids = NULL) {
  p <- recipe$p; G <- length(cat$or_blocks); incidence <- domains <- forcing <- 0
  for (b in seq_along(cat$or_blocks)) {
    ctrl <- cat$or_blocks[[b]]$control_cell_ids
    target <- vapply(cat$links[vapply(cat$links, function(l) l$or_block_id == b, logical(1))], `[[`, integer(1), "target_cell_id")
    incidence <- incidence + length(ctrl); domains <- domains + length(union(ctrl, target))
    if (!is.null(theta_ids)) forcing <- forcing + length(unique(unlist(theta_ids[ctrl])))
  }
  group_fields <- G * (8192 + 128 * (p + 2)^2 + 64 * sum(nchar(recipe$columns, type = "bytes"))) +
    64 * (incidence + domains + length(cat$links)) * (p + 2)^2 + 32 * forcing * (p + 4)
  # shared aliases shared_or in the current source and is charged twice.
  ipw_owner_safe_sum(c(pooled[!names(pooled) %in% c("shared_or", "shared")], 2 * group_fields))
}

reconstruction_fit_extra <- function(cat, recipe, st) {
  p <- recipe$p; peak <- aliases <- 0
  for (b in seq_along(cat$or_blocks)) {
    ctrl <- cat$or_blocks[[b]]$control_cohort_indices
    target <- vapply(cat$links[vapply(cat$links, function(l) l$or_block_id == b, logical(1))], `[[`, integer(1), "target_cohort_index")
    ids <- union(ctrl, target); peak <- max(peak, length(ids))
    aliases <- aliases + sum(st$cohort_end[ids] - st$cohort_start[ids] + 1)
  }
  list(scratch = 64 * (peak + 2) * (p + 2)^2 + 32 * length(cat$links),
    work = 256 * (length(cat$or_blocks) * (p + 2)^3 + aliases * (p + 2)^2 *
      (1 + ceiling(log2(st$n_A + 1)))))
}

# Exact score records are local to one core invocation: the reader closes over
# its fixed rows, cell/group fits, X and weights. Only declared needed cells
# persist, including NULL nonobserver records. Named-row designs retain the
# original path rather than replicating unbounded row-label metadata.
reconstruction_stage1_forward_plan <- function(backend, B, gids, X_F, n_links, details) {
  declined <- list(enabled = FALSE, gids = integer(0L), entries = 0)
  if (details || !reconstruction_separate(backend$shared_or) || !length(gids) ||
      !is.null(rownames(X_F))) return(declined)
  rows <- min(B, backend$control$max_chunk_rows)
  p <- vapply(backend$groups[gids], `[[`, integer(1L), "p")
  column_bytes <- sum(vapply(backend$groups[gids], function(g)
    sum(nchar(g$kept_columns, type = "bytes")), numeric(1L)))
  peak <- max(c(0, p))
  # Retained group rows, names and the peak score/RHS/transposed/pivoted
  # triangular solves, scaling, reorder and returned transpose coexist.
  # Fund ALL needed groups or leave the exact V6 scalar adjoint route intact.
  entries <- 1024 + 4 * rows * sum(p) + 16 * rows * peak +
    4 * peak^2 + 256 * length(gids) + 8 * column_bytes
  if (!is.finite(entries) || entries > 32 * rows * n_links) return(declined)
  list(enabled = TRUE, gids = gids, entries = entries)
}

reconstruction_score_cache_plan <- function(backend, B, cids, X_F, n_links, reserved_entries = 0, funding_entries = NULL) {
  declined <- list(enabled = FALSE, cids = integer(0L), entries = 0)
  if (!length(cids) || !is.null(rownames(X_F)))
    return(declined)
  rows <- min(B, backend$control$max_chunk_rows)
  # The unchanged envelope reserves 16*B*(4*nl) = 64*B*nl for row/link
  # work. Pure replay retains phi, three psi matrices and a scaled phi;
  # augmented replay retains phi/scaled phi. Layout storage/aggregation adds
  # at most three such matrices. Leave 32*B*nl for all of these; detailed
  # per-link records and per-link vectors keep their separate old charges.
  # Cache records, including names/headers, fit entirely in the other half.
  # Generic IPW/pooled envelopes retain larger row/link terms; the same slice
  # fits there. Their unchanged global-cell metadata funds the full position
  # map, selected-slot bookkeeping and planning vectors before any reads.
  # Use the core-subchunk cap even when the envelope's B is layout-expanded.
  # Generic callers fund their larger, unchanged 64*B*32nl or64*B*48nl
  # row/link slabs, leaving32*B*nl for all old phi/psi/layout buffers. Separate
  # callers keep their original32*B*nl cache allowance and forward sharing.
  funding <- (if (is.null(funding_entries)) 32 * rows * n_links else funding_entries) - reserved_entries
  if (funding <= 1024) return(declined)
  p <- vapply(backend$cells[cids], function(cf) backend$groups[[cf$group_id]]$p, integer(1))
  column_bytes <- vapply(backend$cells[cids], function(cf)
    sum(nchar(backend$groups[[cf$group_id]]$kept_columns, type = "bytes")), numeric(1))
  record_entries <- 4 * rows * (p + 2) + 256 + 8 * column_bytes
  entries <- 1024; selected <- logical(length(cids))
  for (j in seq_along(cids)) if (record_entries[j] <= funding - entries) {
    selected[j] <- TRUE
    entries <- entries + record_entries[j]
  }
  if (!any(selected)) return(declined)
  list(enabled = TRUE, cids = cids[selected], entries = entries)
}

reconstruction_cached_cell_reader <- function(read, plan, n_cells) {
  # Force the original reader before its caller rebinds read_cell_score to
  # this wrapper; retaining an unforced promise would recursively call itself.
  force(read)
  if (!isTRUE(plan$enabled)) return(read)
  position <- integer(n_cells); position[plan$cids] <- seq_along(plan$cids)
  cache <- vector("list", length(plan$cids)); seen <- logical(length(plan$cids))
  reader <- function(cid) {
    slot <- position[cid]
    if (slot == 0L) return(read(cid))
    if (!seen[slot]) {
      record <- read(cid)
      cache[slot] <<- list(record)
      seen[slot] <<- TRUE
    }
    cache[[slot]]
  }
  view <- function() list(position = position, seen = seen, cache = cache, read = read)
  attr(reader, "reconstruction_cache_view") <- view
  reader
}

# Both reconstructed routes retain the same local OR and stage1/treatment row
# shapes. Selection adds only its existing upstream admission and q channels.
reconstruction_if_envelope <- function(backend, B, links, details, augmented = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- if (!augmented && !is.null(backend$layout)) backend$layout$state else reader$state
  if (!is.null(backend$layout)) B <- B * max(1, backend$layout$by_orig$count)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  nc_all <- length(backend$cells); ng_all <- length(backend$groups)
  nb_all <- length(backend$or_blocks); nt_all <- length(backend$treatments)
  incidence <- 0
  for (b in backend$shared_or$blocks) incidence <- incidence + length(b$control_cell_ids)
  for (j in links) incidence <- incidence + length(backend$links[[j]]$control_cell_ids)
  for (g in backend$groups) incidence <- incidence + length(g$cell_ids)
  if (augmented) {
    for (c in backend$cells) incidence <- incidence + length(c$theta_ids)
    for (ids in backend$up_map$theta_coords) incidence <- incidence + length(ids)
  }
  maps <- 4096 + 64 * (nc_all + ng_all + nb_all + nt_all + length(links) + incidence)
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    maps, c(backend$accounting$work_ops_bound, 64 * maps), 0, "reconstructed OR dependency metadata")
  sets <- if (augmented) augmented_if_needed_sets(backend, links, details) else
    imputation_if_needed_sets(backend, links, backend$est_method, details)
  nc <- length(sets$cids); ng <- length(sets$gids); nb <- length(sets$bids); nt <- length(sets$tids)
  nl <- length(links); J <- length(st$cohorts); p <- backend$outcome_recipe$p; pt <- backend$treatment_recipe$p
  local_p <- vapply(backend$shared_or$outcome_groups$groups[sets$bids], function(g) length(g$keep), integer(1))
  imp_p <- vapply(backend$groups, `[[`, integer(1), "p")
  # The augmented OR multiplication materializes one local forcing transpose.
  operator_peak <- if (augmented) max(c(0, vapply(backend$shared_or$blocks[sets$bids],
    function(b) as.double(length(b$q_beta)), numeric(1)))) else 0
  nu <- if (augmented) length(sets$uids) else 0L
  nz <- if (augmented) length(sets$zids) else 0L
  name_bytes <- sum(nchar(c(backend$outcome_recipe$columns, backend$treatment_recipe$columns), type = "bytes"))
  link_bytes <- sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  metadata <- 8192 + 8 * link_bytes + B * frame$id_byte_metadata$max +
    if (details) 2048 * (nc + ng + nb + nt + nl + 1) +
      32 * name_bytes * (ng + nb + nt + nl + 1) + 32 * (nc_all + nt_all + nu) else 256 * nl
  channels <- if (augmented) 7 else 6
  width <- nl + 8 + if (details) 3 * J + 2 * nc_all + sum(imp_p) +
    channels * sum(local_p + 1) + nt * 5 * pt + 64 * nl + nu + 8 else 0
  output <- metadata + B * width
  X_F <- if (augmented) augmented_backend_X_F(backend) else imputation_backend_X_F(backend, reader)
  stage1_forward <- reconstruction_stage1_forward_plan(backend, B, sets$gids, X_F, nl, details)
  score_cache <- reconstruction_score_cache_plan(backend, B, sets$cids, X_F, nl,
    reserved_entries = stage1_forward$entries)
  # xi_a/tau/eta, requested treatment rows and X_chunk exist in ordinary replay;
  # stage1 forcing uses scalar adjoints and one local OR operator at a time.
  scratch <- maps + output + 2 * operator_peak + 16 * B * (3 * J + length(backend$design_columns) +
    5 * nt * pt + 4 * nl + 12 * (max(c(0, local_p)) + 1) + nu + nz + 24) +
    64 * nc_all * (p + 2) + 128 * ((p + 2)^2 + (pt + 2)^2)
  # A scalar adjoint still solves every stage1 group in the existing kernel;
  # its score replay reads only cids. Charge that visible global solve loop.
  adjoint <- sum((imp_p + 1)^3) + nc_all * (p + 2)
  work <- 128 * maps + 128 * B * ((nc + 1) * (p + 2)^2 *
    (1 + ceiling(log2(st$n_A + 1))) + (sum(local_p + 1) + 2 * nl) * (adjoint + nc * (p + 2)) +
    nt * (pt + 2) * (pt + J + 4) + J + nl * (p + pt + 8) + (nz + nb + 2 * nl) * (nu + 4))
  if (augmented) {
    sp <- selection_if_budget(backend$selection_backend, B, sets$q_ids, details = details)
    scratch <- scratch + sp$output + sp$scratch; work <- work + sp$work
  }
  c(sets, list(output = output, scratch = scratch, work = work, score_cache = score_cache,
    stage1_forward = stage1_forward))
}
