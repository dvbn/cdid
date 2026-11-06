# One ordinary producer owns these row-independent, originally solved OR
# pullbacks. Scores, centered designs and IF rows never enter this plan.
augmented_adjoint_replay_plan <- function(backend, links, B, phase_store = NULL) {
  validate_augmented_backend(backend)
  if (!backend$outcome_pooling %in% c("pooled", "separate", "inactive") ||
      !backend$imputation_pooling %in% c("separate", "within_cohort", "across_cohorts") ||
      !is.null(backend$layout) || !is.null(backend$reference_ref) ||
      !is.null(backend$pair_outcomes)) return(NULL)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) ||
      B < 1 || B != floor(B) || B > backend$control$max_chunk_rows)
    augmented_abort("Invalid explicit OR adjoint producer partition.")
  pre <- augmented_if_budget(backend, B, links, FALSE)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N
  separate_or <- reconstruction_separate(backend$shared_or)
  has_pooled_or <- identical(backend$outcome_pooling, "pooled") && !is.null(backend$shared_or)
  if (B > N || (!identical(backend$outcome_pooling, "inactive") && !length(pre$bids))) return(NULL)
  if (separate_or) {
    sep <- reconstruction_if_envelope(backend, B, links, FALSE, augmented = TRUE)
    if (!identical(pre$scratch, sep$scratch)) return(NULL)
    rm(sep)
    if (B != 1L) {
      one <- augmented_if_budget(backend, 1L, links, FALSE)
      sep <- reconstruction_if_envelope(backend, 1L, links, FALSE, augmented = TRUE)
      if (!identical(one$scratch, sep$scratch)) return(NULL)
      rm(one, sep)
    }
  } else {
    # Fund only the actual V44 inventory, not an assumed legacy compartment.
    # Both affine endpoints also cover every issued full/tail row count.
    sp <- selection_if_budget(backend$selection_backend, B, pre$q_ids, details = FALSE)
    inv <- reconstruction_ordinary_inventory(backend, B, links, pre, sp)
    if (!identical(pre$scratch, inv$scratch)) return(NULL)
    rm(sp, inv)
    if (B != 1L) {
      one <- augmented_if_budget(backend, 1L, links, FALSE)
      sp <- selection_if_budget(backend$selection_backend, 1L, one$q_ids, details = FALSE)
      inv <- reconstruction_ordinary_inventory(backend, 1L, links, one, sp)
      if (!identical(one$scratch, inv$scratch)) return(NULL)
      rm(one, sp, inv)
    }
  }
  cells <- backend$cells; groups <- backend$groups
  nc <- length(cells); ng <- length(groups)
  pg <- vapply(groups, `[[`, integer(1L), "p")
  P <- sum(pg); r <- max(c(0, pg))
  q <- if (has_pooled_or) backend$shared_or$p else 0L
  alpha <- if (has_pooled_or && backend$outcome_recipe$intercept)
    sort(unique(vapply(backend$links[links], `[[`, integer(1L), "or_block_id"))) else integer(0L)
  dense <- nc + 9 * ng + P + 55
  alpha_bounds <- if (length(alpha)) vapply(alpha, function(bid) {
    cids <- backend$shared_or$blocks[[bid]]$control_cell_ids
    gids <- unique(vapply(cells[cids], `[[`, integer(1L), "group_id"))
    # Existing structural compact representation, including numerical zeros.
    min(dense, 512 + 4 * length(cids) + 32 * length(gids) + 4 * sum(pg[gids]))
  }, numeric(1L)) else numeric(0L)
  seed_bound <- if (has_pooled_or) 2048 + 256 * ng + 4 * nc + 4 * P else 0
  cell_columns <- sum(vapply(cells, function(cf) pg[cf$group_id], numeric(1L)))
  if (separate_or) {
    retained <- 2048 + 16 * (length(links) + length(pre$bids))
    peak <- 2048 + 16 * (nc + ng) + 32 * r^2
    funded <- 4096 + 64 * (nc + ng + length(backend$or_blocks) + length(links) + length(pre$uids)) +
      64 * nc * (backend$outcome_recipe$p + 2) + 128 * (backend$outcome_recipe$p + 2)^2
  } else {
    names_bytes <- sum(nchar(c(backend$design_columns, backend$outcome_recipe$columns,
      backend$treatment_recipe$columns, colnames(reader$state$Z)), type = "bytes"))
    retained <- 4096 + 32 * (q + length(alpha)) + q * (dense + 128) +
      sum(alpha_bounds + 128) + seed_bound
    # Original input loadings, dense solved/replacement record, old/new inactive
    # templates, compact packing, scalar solves and the guarded C score peak.
    peak <- 8192 + 16 * nc + 192 * ng + 4 * cell_columns + 4 * P +
      32 * r^2 + 64 * names_bytes + (if (has_pooled_or) 2 * dense else 0)
    funded <- 16384 + 320 * nc + 512 * ng + 8 * cell_columns + 8 * P +
      32 * r^2 + 80 * names_bytes
  }
  tail <- N %% B; min_B <- if (tail) min(B, tail) else B
  # The constant compartments fund retention even between chunks. Score
  # output/centering uses only the original adjoint row compartment.
  if (!is.finite(retained + peak) || retained + peak > funded ||
      1024 + min_B * (4 * r + 8) + 8 * r > 8192 + min_B * (12 * r + 32)) return(NULL)
  # Pair IDs belong to catalog metadata, not these fixed numerical models.
  # Compare only original fitted operands; keep the first requested model.
  tids <- pre$tids
  representatives <- tids
  fields <- c("target_cohort_index", "all_cohort_indices", "coefficients",
    "H_chol", "a_forcing", "column_names")
  for (pos in seq_along(tids)) {
    tf <- backend$treatments[[tids[pos]]]
    if (pos > 1L) for (prior in seq_len(pos - 1L)) {
      if (representatives[prior] != tids[prior]) next
      other <- backend$treatments[[tids[prior]]]
      same <- TRUE
      for (nm in fields) if (!identical(tf[[nm]], other[[nm]])) { same <- FALSE; break }
      if (same) { representatives[pos] <- tids[prior]; break }
    }
  }
  # Two ordered integer maps and scalar comparison/list headers use only the
  # existing fixed treatment metadata, even between chunks or at B=1.
  treatment_retained <- 256 + 8 * length(tids)
  treatment_peak <- 256 + 8 * length(tids)
  treatment_funded <- 512 + 64 * (length(tids) + length(backend$treatments))
  if (treatment_retained + treatment_peak > treatment_funded) return(NULL)
  # Keep the exact issued full-row preflight, including nonlinear cache
  # selection. No extrapolation to tail rows or different link orders.
  replay_preflight <- NULL; preflight_retained <- 0
  preflight_cost <- 128 + unname(cp04_payload(pre)["entries"])
  if (is.null(attributes(links)) && is.finite(preflight_cost) && retained + peak + preflight_cost <= funded) {
    replay_preflight <- pre
    preflight_retained <- preflight_cost
  }
  store <- new.env(parent = emptyenv())
  store$entries <- vector("list", q + length(alpha)); store$seed <- NULL
  store$active <- TRUE
  store$cell_replay <- augmented_cell_replay_store(backend, pre, phase_store)
  plan <- list(kind = "augmented_adjoint_replay_plan", backend_ref = backend$ref,
    shared = backend$shared_or, links = as.integer(links), bids = pre$bids,
    alpha = alpha, q = q, B = B, min_B = min_B, N = N, store = store,
    treatment_ids = tids, treatment_representatives = representatives,
    treatment_retained_bound = treatment_retained, treatment_construction_peak = treatment_peak,
    treatment_funded_entries = treatment_funded,
    replay_preflight = replay_preflight, replay_preflight_retained_bound = preflight_retained,
    entry_bounds = c(rep(dense, q), alpha_bounds), retained_bound = retained + treatment_retained + preflight_retained,
    construction_peak = peak + treatment_peak, funded_entries = funded + treatment_funded)
  ipw_owner_seal(plan, "augmented_adjoint_replay_plan")
}

# The public wrapper has already validated the actual backend, row/link
# indices and ordinary route. A miss uses its original complete calculation.
augmented_adjoint_replay_preflight <- function(plan, backend, links, B) {
  if (!is.null(attributes(links)) || is.null(plan) || !is.list(plan) ||
      !identical(plan$kind, "augmented_adjoint_replay_plan") ||
      !ipw_owner_seal_ok(plan, "augmented_adjoint_replay_plan") ||
      !identical(plan$backend_ref, backend$ref) || !identical(plan$shared, backend$shared_or) ||
      !identical(plan$links, as.integer(links)) || !isTRUE(B == plan$B) ||
      !isTRUE(plan$store$active)) return(NULL)
  pre <- plan$replay_preflight
  if (is.null(pre)) return(NULL)
  # Cheap actual arithmetic remains at EVERY use with current owner/control
  # counts; only invariant metadata/inventory construction is reused.
  cp04_budget(backend$control,
    c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    pre$scratch, c(backend$accounting$work_ops_bound, pre$work), pre$output,
    "augmented observation-correction IF chunk")
  pre
}

augmented_adjoint_replay_access <- local({
  helpers <- c("or_pooled_wls_gram_solve", "or_wls_gram_solve", "reconstruction_zero_template",
    "reconstruction_zero_template_run", "reconstruction_pack_plain_adjoint",
    "reconstruction_pack_adjoint", "reconstruction_adjoint_lambda", "reconstruction_adjoint_nu",
    "augmented_evaluate_if_chunk_core", "evaluate_augmented_if_chunk", "reconstruction_coarse_operation",
    "reconstruction_score_cache_view", "reconstruction_complete_cache_witness")
  defaults <- lapply(helpers, function(nm) get(nm, environment(), inherits = TRUE))
  names(defaults) <- helpers
  primitives <- c("as.double", "sweep", "%*%", "as.vector", "numeric", "+", "-", "*", "/", "[", "[[", "[<-", "%in%", "length", "identical")
  defaults <- c(defaults, setNames(lapply(primitives, function(nm) get(nm, baseenv())), primitives))
  function(plan, backend, links, B, pre) {
    if (is.null(plan)) return(NULL)
    if (!is.list(plan) || !identical(plan$kind, "augmented_adjoint_replay_plan") ||
        !ipw_owner_seal_ok(plan, "augmented_adjoint_replay_plan") ||
        !identical(plan$backend_ref, backend$ref) || !identical(plan$shared, backend$shared_or) ||
        !identical(plan$links, as.integer(links)) || !identical(plan$bids, pre$bids) || !identical(plan$treatment_ids, pre$tids) ||
        !isTRUE(B == plan$B || B == plan$min_B))
      augmented_abort("OR adjoint replay differs from its sealed producer and partition.", "integrity")
    owner <- environment(augmented_adjoint_replay_plan)
    for (nm in names(defaults)) if ((exists(nm, owner, inherits = FALSE) && bindingIsActive(nm, owner)) ||
        !identical(get(nm, owner, inherits = TRUE), defaults[[nm]])) return(NULL)
    store <- plan$store
    if (!isTRUE(store$active)) return(NULL)
    seed <- store$seed
    if (!is.null(seed) && !ipw_owner_seal_ok(seed, "augmented_adjoint_seed")) {
      store$active <- FALSE; return(NULL)
    }
    slot <- function(kind, id) if (identical(kind, "beta")) {
      if (id >= 1L && id <= plan$q) id else 0L
    } else {
      pos <- match(id, plan$alpha)
      if (is.na(pos)) 0L else plan$q + pos
    }
    get <- function(kind, id) {
      k <- slot(kind, id)
      if (!store$active || !k) return(NULL)
      x <- store$entries[[k]]
      if (!is.null(x) && !ipw_owner_seal_ok(x, "augmented_or_adjoint")) {
        store$active <- FALSE; return(NULL)
      }
      if (is.null(x)) NULL else x$adjoint
    }
    put <- function(kind, id, adjoint, templates) {
      k <- slot(kind, id)
      if (!store$active || !k || !is.null(store$entries[[k]])) return(invisible(NULL))
      if (unname(cp04_payload(adjoint)["entries"]) > plan$entry_bounds[k])
        augmented_abort("OR adjoint replay retained bound underestimated.", "bound")
      store$entries[[k]] <- ipw_owner_seal(list(adjoint = adjoint), "augmented_or_adjoint")
      if (identical(kind, "beta") && is.null(store$seed))
        store$seed <- ipw_owner_seal(list(templates = templates), "augmented_adjoint_seed")
      invisible(NULL)
    }
    list(get = get, put = put, seed = if (is.null(seed)) NULL else seed$templates,
      treatment_ids = plan$treatment_ids, treatment_representatives = plan$treatment_representatives,
      cell_replay = store$cell_replay)
  }
})
lockBinding("defaults", environment(augmented_adjoint_replay_access))
lockBinding("helpers", environment(augmented_adjoint_replay_access))
lockBinding("primitives", environment(augmented_adjoint_replay_access))

# Invoked only at the existing issued producer's return/error boundary.
augmented_adjoint_replay_release <- function(plan) {
  if (is.list(plan) && identical(plan$kind, "augmented_adjoint_replay_plan") &&
      is.environment(plan$store)) {
    plan$store$active <- FALSE
    if (!is.null(plan$store$cell_replay)) plan$store$cell_replay$release()
    plan$store$cell_replay <- NULL
    plan$store$entries <- NULL
    plan$store$seed <- NULL
  }
  invisible(NULL)
}
