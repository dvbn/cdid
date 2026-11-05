# Issued ordinary conditional separate DR metadata only. No fitted/numeric
# operator is prepared; original chunk arithmetic and source operations remain.
dr_conditional_replay_maps <- function(backend, links) {
  incidence <- 0
  for (r in backend$regression_specs) incidence <- incidence + length(r$support_ids)
  for (s in backend$treatment_specs) incidence <- incidence + length(s$support_ids)
  for (u in backend$supports) incidence <- incidence + length(u$channels)
  for (m in backend$treatments) incidence <- incidence + length(m$upstream_ids)
  4096 + 32 * (length(backend$regression_specs) + length(backend$supports) + backend$upstream$D +
    length(backend$treatments) + length(links) + incidence)
}

dr_conditional_replay_preparation <- function(backend, links, B) {
  validate_dr_corollary1_backend(backend)
  if (!dr_separate(backend) || !isTRUE(backend$uses_selection_models)) return(NULL)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) ||
      B < 1 || B != floor(B) || B > backend$control$max_chunk_rows)
    dr_corollary1_abort("Invalid conditional DR preparation link batch or row partition.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N
  if (N <= B) return(NULL)
  maps <- dr_conditional_replay_maps(backend, links)
  # The original map admission precedes even the temporary old envelope maps.
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    maps, c(backend$accounting$work_ops_bound, 64 * maps), 0, "separate DR dependency metadata")
  full <- dr_separate_if_envelope(backend, B, links, FALSE)
  tail_B <- N %% B; if (!tail_B) tail_B <- B
  sel <- backend$selection_backend
  G <- length(backend$regression_specs); S <- length(backend$supports)
  Tn <- length(backend$treatments); M <- length(sel$blocks); C <- length(sel$cells)
  ng <- length(full$gids); nt <- length(full$tids); ns <- length(full$sids); nu <- length(full$uids); nl <- length(links)
  up_incidence <- 0
  for (t in full$tids) up_incidence <- up_incidence + length(backend$treatments[[t]]$upstream_ids)
  active_n <- theta_n <- cell_n <- 0
  for (m in seq_len(M)) if (any(backend$up_map$theta_coords[[m]] %in% full$uids)) {
    active_n <- active_n + 1
    theta_n <- theta_n + length(backend$up_map$theta_coords[[m]])
    cell_n <- cell_n + length(sel$blocks[[m]]$cell_ids)
  }
  # Two exact full/tail dependency records, three split lists with numeric
  # names, global support/model/selection slots, and requested index vectors.
  # Leaf budget cells use the q/block union even for theta-free blocks,
  # independently of the upstream active-cell index below.
  # No fitted record, design or probability is copied into this plan.
  retained <- 4096 + 2 * unname(cp04_payload(links)["entries"]) +
    4 * (nl + ns + nu + ng + nt + S + Tn + 3 * M + up_incidence) +
    24 * (2 * ng + nt) + 48 * active_n + 8 * (theta_n + cell_n) + 4 * C
  construction <- 4096 + 4 * (nl + ns + nu + M + C + backend$upstream$D) + 16 * (ng + nt)
  reservation <- ipw_owner_safe_sum(c(retained, construction))
  # Half the original row-independent map slab funds preparation + its
  # construction. The other half and every V38 extra metadata term stay free
  # for original per-chunk descriptors/current core; no caller prefix grows.
  if (!is.finite(reservation) || reservation > floor(maps / 2)) return(NULL)
  tail <- if (tail_B == B) full else dr_separate_if_envelope(backend, tail_B, links, FALSE)
  qids <- unique(backend$q_of_s[full$sids])
  selection_full <- selection_if_budget(sel, B, qids, details = FALSE)
  selection_tail <- if (tail_B == B) selection_full else selection_if_budget(sel, tail_B, qids, details = FALSE)
  blocks <- vapply(links, function(j) backend$links[[j]]$regression, integer(1))
  models <- vapply(links, function(j) backend$links[[j]]$model, integer(1))
  by_block <- split(seq_along(links), factor(blocks, levels = full$gids))
  by_model <- split(seq_along(links), factor(models, levels = full$tids)); names(by_model) <- as.character(full$tids)
  treatment_blocks <- vapply(full$tids, function(t) backend$treatment_specs[[t]]$regression, integer(1))
  models_by_block <- split(full$tids, factor(treatment_blocks, levels = full$gids))
  support_pos <- integer(S); support_pos[full$sids] <- seq_along(full$sids)
  model_up <- vector("list", Tn)
  for (t in full$tids) model_up[[t]] <- match(backend$treatments[[t]]$upstream_ids, full$uids)
  active <- which(vapply(backend$up_map$theta_coords, function(tc) any(tc %in% full$uids), logical(1)))
  keep <- destination <- cohorts <- vector("list", M)
  active_cells <- integer(0L)
  for (m in active) {
    bl <- sel$blocks[[m]]; tc <- backend$up_map$theta_coords[[m]]
    keep[[m]] <- which(tc %in% full$uids)
    destination[[m]] <- match(tc[keep[[m]]], full$uids)
    cohorts[[m]] <- vapply(sel$cells[bl$cell_ids], `[[`, integer(1), "cohort")
    active_cells <- c(active_cells, bl$cell_ids)
  }
  index <- list(rate_positions = which(full$uids <= backend$up_map$n_a), active_blocks = active,
    active_cells = unique(active_cells), keep = keep, destination = destination, cohorts = cohorts)
  plan <- list(kind = "dr_conditional_replay_preparation", backend_ref = backend$ref,
    reader_ref = backend$reader_ref, selection_ref = sel$ref, reference_ref = backend$reference_ref,
    frame_ref = reader$frame_ref, panel_ref = reader$panel_ref, links = links,
    N = N, B = B, tail_B = tail_B, maps = maps, full = full, tail = tail,
    selection_full = selection_full, selection_tail = selection_tail,
    by_block = by_block, by_model = by_model, models_by_block = models_by_block,
    support_pos = support_pos, model_up = model_up, upstream_index = index,
    active_count = active_n, theta_count = theta_n, cell_count = cell_n,
    retained_bound = retained, construction_peak = construction, reservation = reservation, entries = 0)
  plan$entries <- unname(cp04_payload(plan)["entries"]) + 64
  if (plan$entries > retained)
    dr_corollary1_abort("Conditional DR preparation payload bound underestimated.", "bound")
  ipw_owner_seal(plan, "dr_conditional_replay_preparation")
}

dr_conditional_replay_check <- function(plan, backend, links, B) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  if (!is.list(plan) || !identical(plan$kind, "dr_conditional_replay_preparation") ||
      !ipw_owner_seal_ok(plan, "dr_conditional_replay_preparation") ||
      !identical(plan$backend_ref, backend$ref) || !identical(plan$reader_ref, backend$reader_ref) ||
      !identical(plan$selection_ref, backend$selection_backend$ref) ||
      !identical(plan$reference_ref, backend$reference_ref) || !identical(plan$frame_ref, reader$frame_ref) ||
      !identical(plan$panel_ref, reader$panel_ref) || !identical(plan$links, links) ||
      !isTRUE(plan$N == reader$state$N) || !isTRUE(B == plan$B || B == plan$tail_B))
    dr_corollary1_abort("Conditional DR preparation does not belong to this source, link batch and partition.", "integrity")
  full <- B == plan$B
  pre <- if (full) plan$full else plan$tail
  sp <- if (full) plan$selection_full else plan$selection_tail
  # Every original per-chunk resource guard remains; only its already computed
  # metadata value is reused. Selection operation still validates fresh below.
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    plan$maps, c(backend$accounting$work_ops_bound, 64 * plan$maps), 0, "separate DR dependency metadata")
  sel <- backend$selection_backend
  cp04_budget(sel$control, c(sel$accounting$owned_entries, sel$accounting$borrowed_entries),
    sp$scratch, sp$work, sp$output, "selection IF chunk")
  pre
}

dr_conditional_replay_upstream_index <- function(plan, sel, up_map, uids, details) {
  if (!is.list(plan) || !identical(plan$kind, "dr_conditional_replay_preparation") || details ||
      !ipw_owner_seal_ok(plan, "dr_conditional_replay_preparation") ||
      !identical(plan$selection_ref, sel$ref) || !identical(plan$full$uids, uids) ||
      !is.environment(plan$backend_ref) || !environmentIsLocked(plan$backend_ref) ||
      !exists("dr_corollary1_backend", plan$backend_ref, inherits = FALSE) ||
      !bindingIsLocked("dr_corollary1_backend", plan$backend_ref) || bindingIsActive("dr_corollary1_backend", plan$backend_ref))
    cp04_abort("Conditional upstream preparation has different source coordinates.", "integrity")
  original <- get("dr_corollary1_backend", plan$backend_ref, inherits = FALSE)
  if (!identical(original$up_map, up_map) || !identical(original$selection_backend$ref, sel$ref))
    cp04_abort("Conditional upstream preparation has different source coordinates.", "integrity")
  plan$upstream_index
}
