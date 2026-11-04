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
  c(sets, list(output = output, scratch = scratch, work = work))
}
