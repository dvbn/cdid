# Private CP02 fitted-pipeline metadata slice:
#   1. Unexported owner-to-sparse-occupancy adapter (`derive_owner_sparse_occupancy`).
#   2. Thin unexported owner-to-existing-comparison-catalog adapter
#      (`build_owner_comparison_catalog`).
#
# Scope and boundaries:
#   * Uses sealed `cdid_ipw_panel_owner` (and optional sealed `cdid_ipw_frame_owner`)
#     metadata. Occupancy is accumulated exactly from the owner's sealed per-unit
#     calendar indices (CP08 P2): identical (cohort, observation-date sequence)
#     histories are grouped, then each distinct history's pairs are added row by
#     row into one dense calendar buffer per (cohort, s). No outcome, ID value or
#     per-pair reader batch is touched.
#   * Derives exact `(cohort, s, t, count)` rows for canonical A-units (`n_obs >= 2`),
#     while preserving singleton (`n_obs == 1`), unobserved (`n_obs == 0`), explicit-NA,
#     and optional external frame-only metadata in the referenced owner without copying.
#   * Preserves literal stored dates (`panel_owner$all_times`), cohorts
#     (`panel_owner$cohorts`), and anticipation (`anticipation`).
#   * Identifies cohorts and dates exclusively by exact integer dictionary indices
#     (`1..n_cohorts` and `1..n_calendar_dates`); never uses formatted floating-point
#     strings as identity keys.
#   * Passes `targets`, `control_group`, `conventions`, `requested_post`, and explicit
#     `catalog_bounds` directly to the unchanged `build_comparison_catalog`, retaining
#     all approved D10 original and D13 fixed-horizon comparisons.
#   * Counts grouping, accumulation, temporary dictionary, output and catalog
#     payload storage before relevant allocations/callbacks under explicit caller
#     controls (no internal defaults).
#   * Never materializes a compulsory all-unit-pair table or N-by-moment cache;
#     reads zero outcome endpoints (`endpoints_read == 0L`) and performs zero outcome
#     subtractions, fits, IF/covariance/GMM solves, Monte Carlo, or public wrapping.

cdid_pipeline_metadata_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_pipeline_metadata_refusal",
      "cdid_pipeline_metadata_error",
      "cdid_ipw_owner_error",
      "cdid_catalog_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

# Private helpers this file calls, grouped by the construction stage that
# actually calls them (CP08 stage-local dependency demand). Every entry of this
# file demands only the groups its own stage uses through
# `pipeline_metadata_require_stage`; `pipeline_metadata_require_deps()` (no
# arguments) remains the exhaustive standalone check over every group. Both use
# the existing locked-namespace memo protocol per group: a genuine installed
# namespace resolves each group's function-type check once per session, while
# sourced or cloned function environments re-check on every call. In a
# lazy-loading namespace a function-type check forces the helper's load, so a
# stage never forces the view constructors or the design bridge it does not call.
pipeline_metadata_dependency_groups <- list(
  owners = c(
    "cdid_abort", "gmm_sexp_aliases", "gmm_scalar_numeric", "gmm_scalar_logical", "gmm_scalar_character",
    "cdid_ipw_owner_abort", "ipw_owner_seal", "ipw_owner_seal_ok", "ipw_owner_refuse_decorative_attrs",
    "ipw_owner_int_bound", "ipw_owner_safe_product", "ipw_owner_safe_sum", "ipw_owner_control",
    "validate_ipw_owner_control", "validate_ipw_panel_owner", "validate_ipw_frame_owner", "ipw_owner_reset_trace",
    "catalog_abort", "validate_catalog_bounds"
  ),
  occupancy = c(
    "pair_count_safe", "pair_scalar_vector", "logit_cert_count_output", "ipw_owner_check_interrupt",
    "ipw_owner_match_sorted_double", "ipw_owner_upper_bound", "validate_catalog_occupancy"
  ),
  catalog = c(
    "design_vector", "exact_cmp_sum", "validate_calendar_spec", "validate_catalog_requests",
    "catalog_size", "catalog_product", "build_comparison_catalog"
  ),
  design = c(
    "catalog_adapter_preflight", "catalog_adapter_dimensions", "catalog_design_envelope",
    "comparison_catalog_design_links", "build_comparison_catalog_design",
    "catalog_design_handoff", "build_catalog_handoff_design"
  ),
  views = c("construct_ipw_nuisance_view", "construct_ipw_link_view"),
  # Owner pair readers stay in the exhaustive check; no stage of this file reads
  # pair positions (every construction keeps `endpoints_read = 0L`).
  readers = c("traverse_ipw_panel_pair_positions", "query_ipw_panel_pair_positions")
)

pipeline_metadata_require_stage <- function(groups) {
  ns_env <- environment(pipeline_metadata_require_stage)
  known <- names(pipeline_metadata_dependency_groups)
  if (!is.character(groups) || length(groups) < 1L || anyNA(groups) || !all(groups %in% known)) {
    cdid_pipeline_metadata_abort(
      "Internal dependency stage request names an unknown private helper group.",
      "cdid_pipeline_metadata_dependency_error"
    )
  }
  for (grp in unique(groups)) {
    key <- paste0("pipeline_metadata:", grp)
    if (cdid_dependency_check_done(key, ns_env)) next
    for (fn in pipeline_metadata_dependency_groups[[grp]]) {
      if (!exists(fn, envir = ns_env, mode = "function", inherits = TRUE)) {
        cdid_pipeline_metadata_abort(
          paste0("Required private helper `", fn, "` is not loaded in the namespace."),
          "cdid_pipeline_metadata_dependency_error"
        )
      }
    }
    cdid_dependency_check_record(key, ns_env)
  }
  invisible(TRUE)
}

pipeline_metadata_require_deps <- function() {
  pipeline_metadata_require_stage(names(pipeline_metadata_dependency_groups))
  cdid_dependency_check_record("pipeline_metadata", environment(pipeline_metadata_require_deps))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Explicit private resource controls (NO internal defaults)
# ---------------------------------------------------------------------------

pipeline_metadata_control <- function(owner_control, catalog_bounds) {
  pipeline_metadata_require_stage("owners")
  if (missing(owner_control) || missing(catalog_bounds)) {
    cdid_pipeline_metadata_abort(
      "Both `owner_control` and `catalog_bounds` must be supplied explicitly by the caller (no internal defaults).",
      c("cdid_pipeline_metadata_control_error", "cdid_pipeline_metadata_bound_error",
        "cdid_ipw_owner_control_error", "cdid_catalog_bound_error")
    )
  }
  own <- validate_ipw_owner_control(owner_control)
  cat_b <- validate_catalog_bounds(catalog_bounds)
  structure(
    list(
      owner_control = own,
      catalog_bounds = cat_b
    ),
    class = "cdid_pipeline_metadata_control"
  )
}

validate_pipeline_metadata_control <- function(control, catalog_bounds = NULL) {
  pipeline_metadata_require_stage("owners")
  if (missing(control) || is.null(control) || !is.list(control)) {
    cdid_pipeline_metadata_abort(
      "Explicit `control` is required; no internal defaults exist.",
      c("cdid_pipeline_metadata_control_error", "cdid_pipeline_metadata_bound_error",
        "cdid_ipw_owner_control_error", "cdid_catalog_bound_error")
    )
  }
  if (inherits(control, "cdid_pipeline_metadata_control") ||
      (is.null( dim(control) ) && !is.null(names(control)) &&
       setequal(names(control), c("owner_control", "catalog_bounds")))) {
    if (!is.null(catalog_bounds)) {
      cdid_pipeline_metadata_abort(
        "Supply `catalog_bounds` either inside `pipeline_metadata_control` or as the `catalog_bounds` argument, not both.",
        c("cdid_pipeline_metadata_control_error", "cdid_pipeline_metadata_bound_error")
      )
    }
    return(pipeline_metadata_control(control$owner_control, control$catalog_bounds))
  }
  own <- validate_ipw_owner_control(control)
  if (is.null(catalog_bounds)) {
    cdid_pipeline_metadata_abort(
      "Explicit `catalog_bounds` naming max_candidates, max_rows, max_pool_entries, max_scratch and max_work is required.",
      c("cdid_catalog_bound_error", "cdid_pipeline_metadata_control_error", "cdid_pipeline_metadata_bound_error")
    )
  }
  cat_b <- validate_catalog_bounds(catalog_bounds)
  pipeline_metadata_control(own, cat_b)
}

# ---------------------------------------------------------------------------
# Owner resolution and sealed identity verification
# ---------------------------------------------------------------------------

pipeline_metadata_resolve_owners <- function(panel_owner, frame_owner = NULL) {
  pipeline_metadata_require_stage("owners")
  if (missing(panel_owner) || is.null(panel_owner)) {
    cdid_pipeline_metadata_abort(
      "A sealed `panel_owner` (or `frame_owner`) must be supplied explicitly.",
      c("cdid_pipeline_metadata_owner_error", "cdid_pipeline_metadata_input_error", "cdid_ipw_owner_integrity_error")
    )
  }
  ipw_owner_refuse_decorative_attrs(panel_owner, "panel_owner")
  if (!is.null(frame_owner)) {
    ipw_owner_refuse_decorative_attrs(frame_owner, "frame_owner")
  }

  if (inherits(panel_owner, "cdid_ipw_frame_owner")) {
    validate_ipw_frame_owner(panel_owner)
    f_obj <- panel_owner
    p_obj <- get("panel_owner", envir = f_obj$panel_ref, inherits = FALSE)
    validate_ipw_panel_owner(p_obj)
    if (!is.null(frame_owner)) {
      validate_ipw_frame_owner(frame_owner, panel_owner = p_obj)
      if (!identical(frame_owner$ref, f_obj$ref)) {
        cdid_pipeline_metadata_abort(
          "Supplied `frame_owner` disagrees with the frame owner passed as `panel_owner`.",
          c("cdid_pipeline_metadata_owner_mismatch_error", "cdid_ipw_owner_mismatch_error")
        )
      }
    }
    return(list(panel_owner = p_obj, frame_owner = f_obj))
  }

  validate_ipw_panel_owner(panel_owner)
  f_obj <- NULL
  if (!is.null(frame_owner)) {
    validate_ipw_frame_owner(frame_owner, panel_owner = panel_owner)
    f_obj <- frame_owner
  }
  list(panel_owner = panel_owner, frame_owner = f_obj)
}

pipeline_metadata_owner_summary <- function(panel_owner, frame_owner = NULL) {
  list(
    owner_id = panel_owner$owner_id,
    frame_id = if (!is.null(frame_owner)) frame_owner$frame_id else NULL,
    id_type = panel_owner$id_type,
    id_storage = panel_owner$id_storage,
    id_byte_metadata = panel_owner$id_byte_metadata,
    n_input_rows = panel_owner$counts$n_input_rows,
    n_observed_rows = panel_owner$n_observed_rows,
    n_a_observed_rows = panel_owner$n_a_observed_rows,
    n_explicit_na_rows = panel_owner$n_explicit_na_rows,
    n_total_units = panel_owner$n_total_units,
    n_a_units = panel_owner$n_a_units,
    n_singleton_units = panel_owner$n_singleton_units,
    n_unobserved_units = panel_owner$n_unobserved_units,
    n_frame_only_units = if (!is.null(frame_owner)) frame_owner$n_frame_only_units else 0L,
    n_calendar_dates = panel_owner$n_calendar_dates,
    n_observed_dates = panel_owner$n_observed_dates,
    n_cohorts = length(panel_owner$cohorts),
    n_a_cohorts = length(panel_owner$a_cohorts),
    n_pairs = panel_owner$counts$n_pairs,
    has_frame_owner = !is.null(frame_owner)
  )
}

# ---------------------------------------------------------------------------
# Component 1: Owner-to-Sparse-Occupancy Adapter & Preflight
# ---------------------------------------------------------------------------

preflight_owner_sparse_occupancy <- function(
  panel_owner,
  control,
  frame_owner = NULL
) {
  pipeline_metadata_require_stage(c("owners", "occupancy"))
  owners <- pipeline_metadata_resolve_owners(panel_owner, frame_owner = frame_owner)
  p_own <- owners$panel_owner
  ctrl <- validate_ipw_owner_control(control)

  n_a <- p_own$n_a_units
  n_u <- p_own$n_total_units
  n_cal <- p_own$n_calendar_dates
  n_coh <- length(p_own$cohorts)
  n_a_coh <- length(p_own$a_cohorts)
  n_pairs <- p_own$counts$n_pairs
  q_a <- p_own$n_a_observed_rows

  if (n_cal >= .Machine$integer.max || n_coh >= .Machine$integer.max || n_a_coh >= .Machine$integer.max) {
    cdid_pipeline_metadata_abort(
      "Panel calendar or cohort dictionary count exceeds supported integer index range.",
      c("cdid_pipeline_metadata_bound_error", "cdid_ipw_owner_bound_error", "cdid_catalog_bound_error")
    )
  }

  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  max_date_pairs <- if (n_cal >= 2L) pair_count_safe(n_cal)$total else 0
  max_grid_cells <- if (n_a_coh > 0L && max_date_pairs > 0) {
    ipw_owner_safe_product(n_a_coh, max_date_pairs, label = "cohort-by-date-pair grid upper bound")
  } else {
    0
  }
  max_occ_cells <- min(n_pairs, max_grid_cells)
  if (max_occ_cells >= .Machine$integer.max) {
    cdid_pipeline_metadata_abort(
      "Occupancy cell upper bound exceeds supported data.frame integer row range.",
      c("cdid_pipeline_metadata_bound_error", "cdid_ipw_owner_bound_error", "cdid_catalog_bound_error")
    )
  }

  # 1. Retained and output entry bounds (4 columns in `occupancy`, 5 columns in `cell_indices`, plus scalars)
  retained_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(9, max_occ_cells, max_ret, "retained occupancy and cell_indices columns"),
      128
    ),
    max_bound = max_ret,
    label = "retained occupancy adapter entries"
  )
  output_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(9, max_occ_cells, max_out, "output occupancy and cell_indices columns"),
      128
    ),
    max_bound = max_out,
    label = "occupancy adapter output entries"
  )

  # 2. Exact accumulation plan, fixed from owner counts and caps before any
  #    allocation. Phase 1 labels A-units by (history length, cohort) and then
  #    refines by the calendar index at each position; the units still carrying
  #    a position form a prefix of the length-descending order, so refinement
  #    visits each A-observation at most once. Phase 2 expands each distinct
  #    history's pairs, in slices of at most `slice_pairs` that never split one
  #    (cohort, s) posting (its t-suffix has at most T - 1 dates), into one dense
  #    length-T row buffer. Upper bounds: distinct histories <= n_a, postings
  #    <= n_a_observed_rows, distinct-history pair visits <= n_pairs and output
  #    cells <= max_occupancy_cells. Phase-local scratch is released between phases.
  header_entries <- ipw_owner_safe_sum(
    c(n_a_coh, n_coh, n_cal),
    max_bound = max_me,
    label = "occupancy cohort/calendar dictionary headers"
  )
  grouping_scratch_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(20, n_a, max_me, "history grouping unit order, labels and refinement keys"),
      ipw_owner_safe_product(2, n_cal, max_me, "history length tabulation"),
      256
    ),
    max_bound = max_me,
    label = "occupancy history grouping scratch entries"
  )
  accumulation_fixed_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(4, n_a, max_me, "distinct history summaries"),
      ipw_owner_safe_product(11, q_a, max_me, "history (cohort, s) postings, order and permutation copies"),
      ipw_owner_safe_product(4, max_occ_cells, max_me, "preallocated occupancy cell fill buffers"),
      ipw_owner_safe_product(4, n_cal, max_me, "dense calendar row buffer and row emission"),
      256
    ),
    max_bound = max_me,
    label = "occupancy accumulation fixed scratch entries"
  )
  t_span <- max(1, n_cal - 1)
  slice_room <- max_me - (retained_vector_entries + header_entries + accumulation_fixed_entries)
  slice_pairs <- max(t_span, min(65536, max(1, n_pairs), if (slice_room > 0) floor(slice_room / 14) else 0))
  accumulation_scratch_entries <- ipw_owner_safe_sum(
    c(
      accumulation_fixed_entries,
      ipw_owner_safe_product(14, slice_pairs, max_me, "pair slice expansion and radix collapse scratch")
    ),
    max_bound = max_me,
    label = "occupancy accumulation scratch entries"
  )
  assembly_scratch_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(9, max_occ_cells, max_me, "dictionary truncation and table assembly copies"),
      256
    ),
    max_bound = max_me,
    label = "occupancy table assembly scratch entries"
  )

  # 3. Concurrent payload bound (retained tables, headers and the largest phase)
  concurrent_entries <- ipw_owner_safe_sum(
    c(
      retained_vector_entries,
      header_entries,
      max(grouping_scratch_entries, accumulation_scratch_entries, assembly_scratch_entries)
    ),
    max_bound = max_me,
    label = "concurrent grouping/accumulation/dictionary/output entries"
  )

  # 4. Work operations bound for the algorithm actually executed
  sort_factor <- function(x) 1 + ceiling(log2(x + 1))
  grouping_work_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(64 * sort_factor(n_a), n_a, max_wo, "initial history order and final group compaction"),
      ipw_owner_safe_product(32 * sort_factor(n_a), q_a, max_wo, "position-wise history refinement"),
      ipw_owner_safe_product(256, n_cal, max_wo, "refinement position dispatch"),
      256
    ),
    max_bound = max_wo,
    label = "occupancy history grouping work operations"
  )
  row_bound <- min(q_a, ipw_owner_safe_product(n_a_coh, n_cal, label = "cohort/s row upper bound"))
  accumulation_work_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(16, q_a, max_wo, "history posting construction"),
      ipw_owner_safe_product(32 * sort_factor(q_a), q_a, max_wo, "posting (cohort, s) radix order and permutation"),
      ipw_owner_safe_product(32 * sort_factor(slice_pairs), n_pairs, max_wo, "pair slice expansion, radix collapse and row-buffer updates"),
      ipw_owner_safe_product(8, n_pairs, max_wo, "multi-slice row buffer scans"),
      ipw_owner_safe_product(256, row_bound, max_wo, "cohort/s row dispatch"),
      ipw_owner_safe_product(16, max_occ_cells, max_wo, "row cell emission"),
      256
    ),
    max_bound = max_wo,
    label = "occupancy accumulation work operations"
  )
  assembly_work_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(32, max_occ_cells, max_wo, "dictionary truncation, table assembly and count check"),
      ipw_owner_safe_product(8, header_entries, max_wo, "cohort/calendar dictionary mapping work"),
      256
    ),
    max_bound = max_wo,
    label = "occupancy assembly work operations"
  )
  work_ops_bound <- ipw_owner_safe_sum(
    c(grouping_work_ops, accumulation_work_ops, assembly_work_ops),
    max_bound = max_wo,
    label = "total occupancy derivation work operations"
  )

  list(
    n_a_units = n_a,
    n_total_units = n_u,
    n_calendar_dates = n_cal,
    n_cohorts = n_coh,
    n_a_cohorts = n_a_coh,
    n_pairs = n_pairs,
    n_a_observed_rows = q_a,
    max_date_pairs = max_date_pairs,
    max_occupancy_cells = max_occ_cells,
    accumulation_method = "exact_repeated_history_rowwise_dense_buffer",
    slice_pairs = logit_cert_count_output(slice_pairs),
    max_refinement_positions = logit_cert_count_output(n_cal),
    matrix_dim_limit = as.double(.Machine$integer.max),
    vector_length_limit = 2^52,
    exact_count_limit = 2^53 - 1,
    retained_matrix_entries = 0,
    retained_vector_entries = retained_vector_entries,
    total_retained_entries = retained_vector_entries,
    output_vector_entries = output_vector_entries,
    header_entries = header_entries,
    grouping_scratch_entries = grouping_scratch_entries,
    accumulation_scratch_entries = accumulation_scratch_entries,
    assembly_scratch_entries = assembly_scratch_entries,
    concurrent_matrix_entries = 0,
    concurrent_vector_entries = concurrent_entries,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    grouping_work_ops = grouping_work_ops,
    accumulation_work_ops = accumulation_work_ops,
    assembly_work_ops = assembly_work_ops,
    work_ops_bound = work_ops_bound,
    plan = list(
      phase_1 = "exact grouping of A-units by (cohort, observed calendar-index sequence); released before phase 2",
      phase_2 = "distinct-history pairs in slices of at most slice_pairs into one dense length-T buffer per (cohort, s) row, emitted in ascending t",
      distinct_histories_upper = n_a,
      postings_upper = q_a,
      pair_visits_upper = n_pairs,
      fill_buffer_cells_upper = max_occ_cells,
      row_buffer_length = n_cal,
      owner_pair_reader_used = FALSE,
      unit_id_values_read = FALSE,
      chunk_rows_control_used = FALSE
    ),
    all_unit_pair_table_materialized = FALSE,
    endpoints_read = 0L,
    global_pair_subtraction_performed = FALSE,
    peak_rss_promised = FALSE,
    accounting_scope = "Conservative atomic vector/dictionary slot bounds for the exact grouping/accumulation algorithm before any allocation; excludes R SEXP/list headers and OS process peak RSS."
  )
}

preflight_pipeline_occupancy <- preflight_owner_sparse_occupancy

# Exact occupancy accumulation (CP08 P2). Reads only the sealed owner's unit
# observation ranges, A-cohort indices and observation calendar indices; never
# outcomes, unit IDs or a per-pair reader. Cells are returned sorted by
# (A-cohort index, s index, t index), which is the (cohort, s, t) order because
# the A-cohort dictionary is an increasing subsequence of the cohort dictionary.
pipeline_occupancy_accumulate <- function(p_own, preflight) {
  n_a <- p_own$n_a_units
  n_cal <- p_own$n_calendar_dates
  diag <- list(n_distinct_histories = 0, max_history_multiplicity = 0, n_refinement_positions = 0,
    refinement_unit_visits = 0, n_postings = 0, n_rows = 0, n_slices = 0, max_slice_pairs = 0,
    n_history_pair_visits = 0, n_pairs_represented = 0, fill_buffer_cells = 0,
    peak_phase_entries = 0, actual_work_ops = 0)
  if (n_a < 1L) {
    return(list(cohort_a_index = integer(0L), s_index = integer(0L), t_index = integer(0L),
      count = numeric(0L), diagnostics = diag))
  }
  if (pipeline_occupancy_complete_roster(p_own)) return(pipeline_occupancy_accumulate_complete(p_own, preflight, diag))
  cal_idx <- p_own$observations$calendar_index
  a_units <- p_own$a_unit_index

  # Phase 1: labels start from (history length, cohort); position j refines the
  # units of length >= j, a prefix of the length-descending order. Labels issued
  # at a position exceed every earlier label, so they never collide with the
  # final labels of shorter histories. Equal final labels <=> equal cohort,
  # length and calendar index at every position.
  n_len <- p_own$units$n_obs[a_units]
  ord <- order(-n_len, p_own$a_cohort_index, method = "radix")
  nl <- n_len[ord]
  n_len <- NULL
  ks <- p_own$a_cohort_index[ord]
  sts <- p_own$units$obs_start[a_units]
  sts <- sts[ord]
  ord <- NULL
  grp <- as.double(cumsum(c(TRUE, nl[-1L] != nl[-n_a] | ks[-1L] != ks[-n_a])))
  next_id <- grp[n_a]
  max_len <- nl[1L]
  cnt_ge <- rev(cumsum(rev(tabulate(nl, nbins = max_len))))
  steps <- 0
  visits <- 0
  for (j in seq_len(max_len)) {
    p <- cnt_ge[j]
    if (p < 2L) break
    ipw_owner_check_interrupt()
    key <- cal_idx[sts[seq_len(p)] + (j - 1)]
    gp <- grp[seq_len(p)]
    o <- order(gp, key, method = "radix")
    gs <- gp[o]
    kk <- key[o]
    new_id <- next_id + cumsum(c(TRUE, gs[-1L] != gs[-p] | kk[-1L] != kk[-p]))
    grp[o] <- new_id
    steps <- steps + 1
    visits <- visits + p
    if (new_id[p] - next_id == p) break # every remaining unit is alone: labels are final
    next_id <- new_id[p]
  }
  key <- gp <- o <- gs <- kk <- new_id <- cnt_ge <- NULL
  og <- order(grp, method = "radix")
  gsorted <- grp[og]
  first <- which(c(TRUE, gsorted[-1L] != gsorted[-n_a]))
  n_hist <- length(first)
  w_hist <- as.double(diff(c(first, n_a + 1L)))
  rep_pos <- og[first]
  h_len <- nl[rep_pos]
  h_k <- ks[rep_pos]
  h_st <- sts[rep_pos]
  # Release phase-1 scratch before any phase-2 allocation.
  grp <- og <- gsorted <- first <- rep_pos <- nl <- ks <- sts <- NULL

  hl <- as.double(h_len)
  pair_h <- ifelse(hl %% 2 == 0, (hl / 2) * (hl - 1), hl * ((hl - 1) / 2))
  p_hist <- ipw_owner_safe_sum(pair_h, label = "distinct-history pair visits")
  n_rep <- ipw_owner_safe_sum(w_hist * pair_h, label = "represented within-unit pairs")
  q_post <- sum(hl) - n_hist
  hl <- pair_h <- NULL
  cap_cells <- min(preflight$max_occupancy_cells, p_hist)
  slice_pairs <- preflight$slice_pairs
  grouping_entries <- 20 * n_a + 2 * max_len + 256
  plan_entries <- 4 * n_hist + 11 * q_post + 4 * cap_cells + 4 * n_cal + 14 * slice_pairs + 256
  if (plan_entries > preflight$accumulation_scratch_entries || grouping_entries > preflight$grouping_scratch_entries) {
    cdid_pipeline_metadata_abort(
      "Internal occupancy accumulation plan exceeds its preflight reservation.",
      "cdid_pipeline_metadata_internal_error"
    )
  }

  # Phase 2: one posting per (distinct history, position i < n_h), keyed by
  # (cohort, s); its t-suffix holds the history's later calendar indices.
  post_n <- h_len - 1L
  ph <- rep.int(seq_len(n_hist), post_n)
  pos_i <- sequence(post_n)
  post_n <- NULL
  p_from <- h_st[ph] + pos_i
  p_len <- h_len[ph] - pos_i
  p_w <- w_hist[ph]
  pk <- h_k[ph]
  ph <- pos_i <- NULL
  ps <- cal_idx[p_from - 1]
  po <- order(pk, ps, method = "radix")
  pk <- pk[po]
  ps <- ps[po]
  p_from <- p_from[po]
  p_len <- p_len[po]
  p_w <- p_w[po]
  po <- NULL
  n_post <- length(pk)
  row_start <- which(c(TRUE, pk[-1L] != pk[-n_post] | ps[-1L] != ps[-n_post]))
  row_end <- c(row_start[-1L] - 1L, n_post)
  cum_len <- cumsum(as.double(p_len))

  buf <- numeric(n_cal)
  out_k <- integer(cap_cells)
  out_s <- integer(cap_cells)
  out_t <- integer(cap_cells)
  out_c <- numeric(cap_cells)
  n_out <- 0
  n_slices <- 0
  max_slice <- 0
  for (r in seq_along(row_start)) {
    ipw_owner_check_interrupt()
    a <- row_start[r]
    b <- row_end[r]
    multi <- FALSE
    repeat {
      base <- if (a > 1L) cum_len[a - 1L] else 0
      e <- if (cum_len[b] - base <= slice_pairs) b else ipw_owner_upper_bound(cum_len, a, b, base + slice_pairs)
      idx <- seq.int(a, e)
      len <- p_len[idx]
      tt <- cal_idx[rep.int(p_from[idx], len) + (sequence(len) - 1)]
      n_tt <- length(tt)
      if (e == a) {
        # One history suffix: strictly increasing distinct calendar indices.
        ut <- tt
        buf[ut] <- buf[ut] + p_w[a]
      } else {
        o <- order(tt, method = "radix")
        ts <- tt[o]
        cs <- cumsum(rep.int(p_w[idx], len)[o])
        ends <- which(c(ts[-1L] != ts[-n_tt], TRUE))
        ut <- ts[ends]
        buf[ut] <- buf[ut] + diff(c(0, cs[ends]))
        o <- ts <- cs <- ends <- NULL
      }
      n_slices <- n_slices + 1
      max_slice <- max(max_slice, n_tt)
      tt <- idx <- len <- NULL
      if (e >= b) break
      multi <- TRUE
      a <- e + 1L
    }
    if (multi) ut <- which(buf > 0)
    m <- length(ut)
    if (n_out + m > cap_cells) {
      cdid_pipeline_metadata_abort(
        "Internal occupancy cell count exceeds its exact distinct-pair bound.",
        "cdid_pipeline_metadata_internal_error"
      )
    }
    hit <- n_out + seq_len(m)
    out_k[hit] <- pk[row_start[r]]
    out_s[hit] <- ps[row_start[r]]
    out_t[hit] <- ut
    out_c[hit] <- buf[ut]
    buf[ut] <- 0
    n_out <- n_out + m
  }
  n_rows <- length(row_start)
  pk <- ps <- p_from <- p_len <- p_w <- row_start <- row_end <- cum_len <- buf <- ut <- hit <- NULL

  keep <- seq_len(n_out)
  cohort_a_index <- out_k[keep]
  out_k <- NULL
  s_index <- out_s[keep]
  out_s <- NULL
  t_index <- out_t[keep]
  out_t <- NULL
  count <- out_c[keep]
  out_c <- NULL

  diag$n_distinct_histories <- n_hist
  diag$max_history_multiplicity <- max(w_hist)
  diag$n_refinement_positions <- steps
  diag$refinement_unit_visits <- visits
  diag$n_postings <- q_post
  diag$n_rows <- n_rows
  diag$n_slices <- n_slices
  diag$max_slice_pairs <- max_slice
  diag$n_history_pair_visits <- p_hist
  diag$n_pairs_represented <- n_rep
  diag$fill_buffer_cells <- cap_cells
  diag$peak_phase_entries <- max(grouping_entries, plan_entries, 9 * n_out + 256)
  diag$actual_work_ops <- 2 * n_a + visits + 2 * q_post + p_hist + n_rows + 2 * n_out
  list(cohort_a_index = cohort_a_index, s_index = s_index, t_index = t_index, count = count, diagnostics = diag)
}

# CP08 common boundary: exact complete-roster predicate on the sealed owner.
# Dates are unique within unit and belong to the calendar, so the A-observation
# count equals n_a * T exactly when every pair-eligible unit observes every
# calendar date; nothing weaker (for instance marginal date counts) is used.
pipeline_occupancy_complete_roster <- function(p_own) {
  n_a <- p_own$n_a_units
  n_cal <- p_own$n_calendar_dates
  n_a >= 1L && n_cal >= 2L &&
    identical(as.double(p_own$n_a_observed_rows), as.double(n_a) * as.double(n_cal))
}

# Closed form of the general accumulation under the complete-roster predicate:
# every A-unit of a cohort has the same history 1..T, so there is one distinct
# history per A-cohort, every (cohort, s, t) cell with s < t is occupied by the
# cohort's A-unit count, and the cells come out in the general path's
# (A-cohort, s, t) order. The returned counts, types and diagnostics are the
# ones the general grouping/refinement/slice algorithm produces on this input;
# only its T radix refinements and per-row buffer work are not repeated.
pipeline_occupancy_accumulate_complete <- function(p_own, preflight, diag) {
  ipw_owner_check_interrupt()
  n_a <- p_own$n_a_units
  n_cal <- p_own$n_calendar_dates
  n_hist <- length(p_own$a_cohorts)
  w_hist <- as.double(p_own$cohort_a_count)
  max_len <- n_cal
  # Refinement positions of the general phase 1 on identical histories: none
  # with a single A-unit, one when every A-unit is alone in its cohort, else T.
  steps <- if (n_a < 2L) 0 else if (n_hist == n_a) 1 else as.double(n_cal)
  visits <- steps * n_a
  hl <- rep(as.double(n_cal), n_hist)
  pair_h <- ifelse(hl %% 2 == 0, (hl / 2) * (hl - 1), hl * ((hl - 1) / 2))
  p_hist <- ipw_owner_safe_sum(pair_h, label = "distinct-history pair visits")
  n_rep <- ipw_owner_safe_sum(w_hist * pair_h, label = "represented within-unit pairs")
  q_post <- sum(hl) - n_hist
  cap_cells <- min(preflight$max_occupancy_cells, p_hist)
  slice_pairs <- preflight$slice_pairs
  grouping_entries <- 20 * n_a + 2 * max_len + 256
  plan_entries <- 4 * n_hist + 11 * q_post + 4 * cap_cells + 4 * n_cal + 14 * slice_pairs + 256
  if (plan_entries > preflight$accumulation_scratch_entries || grouping_entries > preflight$grouping_scratch_entries) {
    cdid_pipeline_metadata_abort(
      "Internal occupancy accumulation plan exceeds its preflight reservation.",
      "cdid_pipeline_metadata_internal_error"
    )
  }
  # One (cohort, s) row per history posting; its t-suffix is s + 1, ..., T.
  later <- rev(seq_len(n_cal - 1L))
  s_seq <- rep.int(seq_len(n_cal - 1L), later)
  t_seq <- s_seq + sequence(later)
  n_per <- length(s_seq)
  cohort_a_index <- rep.int(seq_len(n_hist), rep.int(n_per, n_hist))
  s_index <- rep.int(s_seq, n_hist)
  t_index <- rep.int(t_seq, n_hist)
  count <- rep(w_hist, each = n_per)
  n_out <- length(cohort_a_index)
  if (n_out > cap_cells) {
    cdid_pipeline_metadata_abort(
      "Internal occupancy cell count exceeds its exact distinct-pair bound.",
      "cdid_pipeline_metadata_internal_error"
    )
  }
  n_rows <- as.integer(n_hist * (n_cal - 1L))

  diag$n_distinct_histories <- n_hist
  diag$max_history_multiplicity <- max(w_hist)
  diag$n_refinement_positions <- steps
  diag$refinement_unit_visits <- visits
  diag$n_postings <- q_post
  diag$n_rows <- n_rows
  diag$n_slices <- as.double(n_rows)
  diag$max_slice_pairs <- as.double(n_cal - 1L)
  diag$n_history_pair_visits <- p_hist
  diag$n_pairs_represented <- n_rep
  diag$fill_buffer_cells <- cap_cells
  diag$peak_phase_entries <- max(grouping_entries, plan_entries, 9 * n_out + 256)
  diag$actual_work_ops <- 2 * n_a + visits + 2 * q_post + p_hist + n_rows + 2 * n_out
  list(cohort_a_index = cohort_a_index, s_index = s_index, t_index = t_index, count = count, diagnostics = diag)
}

derive_owner_sparse_occupancy <- function(
  panel_owner,
  control,
  frame_owner = NULL,
  trace_env = NULL,
  occupancy_id = NULL
) {
  pipeline_metadata_require_stage(c("owners", "occupancy"))
  if (missing(panel_owner) || missing(control)) {
    cdid_pipeline_metadata_abort(
      "Required arguments (`panel_owner`, `control`) must be supplied explicitly.",
      c("cdid_pipeline_metadata_input_error", "cdid_ipw_owner_input_error")
    )
  }
  if (!is.null(occupancy_id) && !gmm_scalar_character(occupancy_id)) {
    cdid_pipeline_metadata_abort(
      "`occupancy_id` must be NULL or a nonempty character scalar.",
      "cdid_pipeline_metadata_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_pipeline_metadata_abort(
      "`trace_env` must be NULL or an environment.",
      "cdid_pipeline_metadata_input_error"
    )
  }

  owners <- pipeline_metadata_resolve_owners(panel_owner, frame_owner = frame_owner)
  p_own <- owners$panel_owner
  f_own <- owners$frame_owner
  ctrl <- validate_ipw_owner_control(control)

  # 1. Preflight BEFORE any dictionary allocation, accumulation or trace reset
  preflight <- preflight_owner_sparse_occupancy(p_own, ctrl, frame_owner = f_own)
  ipw_owner_reset_trace(trace_env)

  n_cal <- p_own$n_calendar_dates
  n_coh <- length(p_own$cohorts)
  n_a_coh <- length(p_own$a_cohorts)
  n_pairs <- p_own$counts$n_pairs

  # Exact integer cohort dictionary map from `a_cohorts` to `cohorts`
  a_to_all_cohort <- if (n_a_coh > 0L) {
    as.integer(ipw_owner_match_sorted_double(p_own$a_cohorts, p_own$cohorts))
  } else {
    integer(0L)
  }
  if (anyNA(a_to_all_cohort)) {
    cdid_pipeline_metadata_abort(
      "Internal owner cohort dictionary mismatch: `a_cohorts` is not a subset of `cohorts`.",
      "cdid_pipeline_metadata_internal_error"
    )
  }

  # 2. Exact outcome-free accumulation from sealed calendar indices (no reader, no outcome or ID values)
  acc <- pipeline_occupancy_accumulate(p_own, preflight)
  acc_diag <- acc$diagnostics
  dict_ki <- a_to_all_cohort[acc$cohort_a_index]
  dict_si <- acc$s_index
  dict_ti <- acc$t_index
  dict_count <- acc$count
  acc <- NULL
  if (anyNA(dict_ki) || anyNA(dict_si) || anyNA(dict_ti) ||
      any(dict_si < 1L | dict_si >= dict_ti | dict_ti > n_cal) ||
      any(dict_ki < 1L | dict_ki > n_coh) || any(dict_count < 1)) {
    cdid_pipeline_metadata_abort(
      "Invalid dictionary index in accumulated occupancy cells.",
      "cdid_pipeline_metadata_internal_error"
    )
  }

  total_occ_count <- ipw_owner_safe_sum(dict_count, label = "total occupancy pair count")
  if (total_occ_count != n_pairs || acc_diag$n_pairs_represented != n_pairs) {
    cdid_pipeline_metadata_abort(
      "Derived occupancy total count disagrees with `panel_owner$counts$n_pairs`.",
      "cdid_pipeline_metadata_count_error"
    )
  }

  n_cells <- length(dict_ki)
  occupancy <- data.frame(
    cohort = as.double(p_own$cohorts[dict_ki]),
    s = as.double(p_own$all_times[dict_si]),
    t = as.double(p_own$all_times[dict_ti]),
    count = as.double(dict_count),
    stringsAsFactors = FALSE
  )
  # Validate with the existing catalog occupancy validator to guarantee schema compatibility
  validate_catalog_occupancy(occupancy)

  cell_indices <- data.frame(
    cell_id = seq_len(n_cells),
    cohort_index = as.integer(dict_ki),
    s_index = as.integer(dict_si),
    t_index = as.integer(dict_ti),
    count = as.double(dict_count),
    stringsAsFactors = FALSE
  )

  actual_retained <- 9 * n_cells + 32
  actual_concurrent <- actual_retained + preflight$header_entries + acc_diag$peak_phase_entries
  actual_ops <- 2 * (n_a_coh + n_coh + n_cal) + acc_diag$actual_work_ops
  stopifnot(
    n_cells <= preflight$max_occupancy_cells,
    actual_retained <= preflight$total_retained_entries,
    actual_concurrent <= preflight$concurrent_entries,
    actual_ops <= preflight$work_ops_bound
  )

  occ_obj <- structure(
    list(
      occupancy_id = occupancy_id,
      panel_ref = p_own$ref,
      frame_ref = if (!is.null(f_own)) f_own$ref else NULL,
      occupancy = occupancy,
      cell_indices = cell_indices,
      calendar = p_own$all_times,
      observed_times = p_own$times,
      cohorts = p_own$cohorts,
      a_cohorts = p_own$a_cohorts,
      owner_summary = pipeline_metadata_owner_summary(p_own, f_own),
      control = ctrl,
      preflight = preflight,
      diagnostics = list(
        endpoints_read = 0L,
        read_outcomes = FALSE,
        global_pair_subtraction_performed = FALSE,
        pair_differences_computed = FALSE,
        all_unit_pair_table_materialized = FALSE,
        n_by_moment_cache_allocated = FALSE,
        formatted_float_keys_used = FALSE,
        dictionary_index_keys_used = TRUE,
        singleton_metadata_preserved_in_owner = TRUE,
        unobserved_metadata_preserved_in_owner = TRUE,
        explicit_na_metadata_preserved_in_owner = TRUE,
        accumulation_method = preflight$accumulation_method,
        owner_pair_reader_used = FALSE,
        unit_id_values_read = FALSE,
        n_distinct_histories = logit_cert_count_output(acc_diag$n_distinct_histories),
        max_history_multiplicity = logit_cert_count_output(acc_diag$max_history_multiplicity),
        n_refinement_positions = logit_cert_count_output(acc_diag$n_refinement_positions),
        refinement_unit_visits = acc_diag$refinement_unit_visits,
        n_postings = acc_diag$n_postings,
        n_rows = logit_cert_count_output(acc_diag$n_rows),
        n_slices = logit_cert_count_output(acc_diag$n_slices),
        max_slice_pairs = logit_cert_count_output(acc_diag$max_slice_pairs),
        slice_pairs = preflight$slice_pairs,
        n_pairs_represented = acc_diag$n_pairs_represented,
        n_history_pair_visits = acc_diag$n_history_pair_visits,
        pair_visits_avoided = n_pairs - acc_diag$n_history_pair_visits,
        n_occupancy_cells = logit_cert_count_output(n_cells),
        total_occupancy_count = logit_cert_count_output(total_occ_count),
        actual_retained_entries = actual_retained,
        actual_concurrent_entries = actual_concurrent,
        actual_work_ops = actual_ops,
        temporary_dictionaries = list(
          cohort_a_to_all_map_length = n_a_coh,
          row_buffer_length = n_cal,
          fill_buffer_cells = acc_diag$fill_buffer_cells
        )
      )
    ),
    class = "cdid_pipeline_occupancy"
  )

  ipw_owner_seal(occ_obj, "pipeline_occupancy")
}

construct_owner_sparse_occupancy <- derive_owner_sparse_occupancy
adapt_owner_to_sparse_occupancy <- derive_owner_sparse_occupancy
chained_pipeline_occupancy <- derive_owner_sparse_occupancy

validate_owner_sparse_occupancy <- function(occupancy_obj, panel_owner = NULL, frame_owner = NULL) {
  pipeline_metadata_require_stage("owners")
  if (missing(occupancy_obj) || is.null(occupancy_obj) ||
      !inherits(occupancy_obj, "cdid_pipeline_occupancy") ||
      !is.list(occupancy_obj) || !ipw_owner_seal_ok(occupancy_obj, "pipeline_occupancy")) {
    cdid_pipeline_metadata_abort(
      "Occupancy object is malformed, foreign, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_pipeline_metadata_integrity_error", "cdid_ipw_owner_integrity_error")
    )
  }
  pref <- occupancy_obj$panel_ref
  if (!is.environment(pref) || !environmentIsLocked(pref) ||
      !exists("panel_owner", envir = pref, inherits = FALSE) ||
      !bindingIsLocked("panel_owner", pref)) {
    cdid_pipeline_metadata_abort(
      "Occupancy object lacks a valid sealed `panel_owner` reference.",
      c("cdid_pipeline_metadata_integrity_error", "cdid_ipw_owner_integrity_error")
    )
  }
  p_obj <- get("panel_owner", envir = pref, inherits = FALSE)
  validate_ipw_panel_owner(p_obj)

  if (!is.null(panel_owner)) {
    owners <- pipeline_metadata_resolve_owners(panel_owner, frame_owner = frame_owner)
    if (!identical(pref, owners$panel_owner$ref)) {
      cdid_pipeline_metadata_abort(
        "Occupancy object belongs to a different `cdid_ipw_panel_owner` snapshot.",
        c("cdid_pipeline_metadata_owner_mismatch_error", "cdid_ipw_owner_mismatch_error")
      )
    }
  }
  if (!is.null(occupancy_obj$frame_ref)) {
    fref <- occupancy_obj$frame_ref
    if (!is.environment(fref) || !environmentIsLocked(fref) ||
        !exists("frame_owner", envir = fref, inherits = FALSE) ||
        !bindingIsLocked("frame_owner", fref)) {
      cdid_pipeline_metadata_abort(
        "Occupancy object `frame_ref` is malformed.",
        c("cdid_pipeline_metadata_integrity_error", "cdid_ipw_owner_integrity_error")
      )
    }
    f_obj <- get("frame_owner", envir = fref, inherits = FALSE)
    validate_ipw_frame_owner(f_obj, panel_owner = p_obj)
    if (!is.null(frame_owner)) {
      validate_ipw_frame_owner(frame_owner, panel_owner = p_obj)
      if (!identical(fref, frame_owner$ref)) {
        cdid_pipeline_metadata_abort(
          "Occupancy object belongs to a different `cdid_ipw_frame_owner` snapshot.",
          c("cdid_pipeline_metadata_owner_mismatch_error", "cdid_ipw_owner_mismatch_error")
        )
      }
    }
  }
  if (!is.data.frame(occupancy_obj$occupancy) || ncol(occupancy_obj$occupancy) != 4L ||
      !is.data.frame(occupancy_obj$cell_indices) || ncol(occupancy_obj$cell_indices) != 5L ||
      nrow(occupancy_obj$occupancy) != nrow(occupancy_obj$cell_indices)) {
    cdid_pipeline_metadata_abort("Sealed occupancy headers disagree.", "cdid_pipeline_metadata_integrity_error")
  }
  if (!is.null(frame_owner)) validate_ipw_frame_owner(frame_owner, panel_owner = p_obj)
  invisible(occupancy_obj)
}

validate_pipeline_occupancy <- validate_owner_sparse_occupancy

# ---------------------------------------------------------------------------
# Component 2: Thin Owner-to-Existing-Comparison-Catalog Adapter & Preflight
# ---------------------------------------------------------------------------

pipeline_metadata_validate_catalog_inputs_early <- function(
  calendar,
  cohorts,
  targets,
  anticipation,
  control_group,
  conventions,
  requested_post,
  catalog_bounds,
  n_occ_lower = 0,
  n_occ_upper = 0,
  owner_control,
  occupancy_retained,
  work_prefix
) {
  bnd <- validate_catalog_bounds(catalog_bounds)
  nT0 <- length(calendar)
  nC0 <- length(cohorts)
  nJ0 <- length(targets)
  nQ0 <- if (is.data.frame(requested_post)) nrow(requested_post) else 0

  if (nT0 >= .Machine$integer.max || nC0 >= .Machine$integer.max ||
      nJ0 > .Machine$integer.max || n_occ_lower >= .Machine$integer.max ||
      n_occ_upper >= .Machine$integer.max || nQ0 > .Machine$integer.max) {
    catalog_abort(
      "Catalog dictionary/data-frame index header exceeds supported integer range.",
      "cdid_catalog_bound_error"
    )
  }

  H_min <- catalog_size(c(nT0, nC0, nJ0, n_occ_lower, nQ0, 1), "header dimensions")
  H_max <- catalog_size(c(nT0, nC0, nJ0, n_occ_upper, nQ0, 1), "header dimensions upper bound")
  byte_conventions <- if (is.character(conventions) && !is.object(conventions) &&
                          length(conventions) == 1L && !anyNA(conventions)) {
    as.double(nchar(conventions, type = "bytes"))
  } else {
    0
  }

  min_header_scratch <- catalog_size(
    c(
      catalog_product(64, H_min, "header payload"),
      catalog_product(16, byte_conventions, "convention bytes")
    ),
    "header scratch"
  )
  min_header_work <- catalog_size(
    c(
      catalog_product(64 * (1 + ceiling(log2(H_min + 1))), H_min, "header validation/sort work"),
      catalog_product(16, byte_conventions, "convention byte work")
    ),
    "header work"
  )
  max_header_scratch <- catalog_size(
    c(
      catalog_product(64, H_max, "header payload upper bound"),
      catalog_product(16, byte_conventions, "convention bytes")
    ),
    "header scratch upper bound"
  )
  max_header_work <- catalog_size(
    c(
      catalog_product(64 * (1 + ceiling(log2(H_max + 1))), H_max, "header validation/sort work upper bound"),
      catalog_product(16, byte_conventions, "convention byte work")
    ),
    "header work upper bound"
  )

  if (min_header_work > bnd$max_work) {
    catalog_abort("Header validation/sort work exceeds explicit `max_work` before scans.", "cdid_catalog_bound_error")
  }
  if (min_header_scratch > bnd$max_scratch) {
    catalog_abort("Scratch for header validation/sort buffers exceeds explicit `max_scratch` before allocation.", "cdid_catalog_bound_error")
  }

  # Caller-owned sealed occupancy and the early header buffers coexist.
  # Refuse the owner's smaller caps before validating/scanning calendar contents.
  ipw_owner_safe_sum(c(occupancy_retained),min(owner_control$max_retained_entries,owner_control$max_output_entries),"retained supplied occupancy before header scans")
  ipw_owner_safe_sum(c(occupancy_retained,min_header_scratch),owner_control$max_matrix_entries,"occupancy plus early header buffers before scans")
  ipw_owner_safe_sum(c(work_prefix,min_header_work),owner_control$max_work_ops,"composed early header work before scans")

  spec <- validate_calendar_spec(calendar, cohorts, anticipation)
  cal <- spec$calendar
  coh <- spec$cohorts
  nT <- length(cal)
  nC <- length(coh)
  if (as.double(nT) * nC > 2^53 - 1) {
    catalog_abort("Calendar/cohort dictionary product exceeds exact index arithmetic.", "cdid_catalog_bound_error")
  }
  if (!design_vector(control_group, 1L, "character") ||
      !control_group %in% c("nevertreated", "notyettreated")) {
    catalog_abort("`control_group` must be explicitly \"nevertreated\" or \"notyettreated\".", "cdid_catalog_input_error")
  }
  if (!design_vector(conventions, 1L, "character") || !nzchar(conventions)) {
    catalog_abort(
      "Declare the common population universe/mixture/covariate/weight conventions as one nonempty opaque string.",
      "cdid_catalog_input_error"
    )
  }
  if (!design_vector(targets, length(targets)) || length(targets) < 1L || anyDuplicated(as.double(targets))) {
    catalog_abort("`targets` must be a nonempty vector of unique finite cohorts.", "cdid_catalog_target_error")
  }
  tg <- match(as.double(targets), coh)
  if (anyNA(tg)) {
    catalog_abort("Every target must belong to the supplied cohort universe.", "cdid_catalog_target_error")
  }
  tg <- sort(tg)
  req <- validate_catalog_requests(requested_post, cal, coh, tg)

  # Pre-check exact_cmp_sum on the cutoff table before pair traversal so arithmetic overflow
  # on (calendar, cohorts, anticipation) refuses immediately before any reader callbacks!
  ci <- 1L
  for (d in seq_len(nT)) {
    repeat {
      below <- exact_cmp_sum(cal[d], spec$anticipation, coh[min(ci, nC)]) < 0L
      if (ci > nC || below) break
      ci <- ci + 1L
    }
  }

  list(
    spec = spec,
    target_indices = tg,
    requested = req,
    byte_conventions = byte_conventions,
    min_header_scratch = min_header_scratch,
    min_header_work = min_header_work,
    max_header_scratch = max_header_scratch,
    max_header_work = max_header_work
  )
}

preflight_owner_comparison_catalog <- function(
  panel_owner,
  targets,
  anticipation,
  control_group,
  conventions,
  requested_post = NULL,
  control,
  catalog_bounds = NULL,
  occupancy = NULL,
  frame_owner = NULL,
  declared_qualified = NULL,
  target_post_times = NULL,
  build_design = FALSE
) {
  pipeline_metadata_require_stage(c("owners", "occupancy", "catalog", if (isTRUE(build_design)) "design"))
  if (missing(panel_owner) || missing(targets) || missing(anticipation) ||
      missing(control_group) || missing(conventions) || missing(control)) {
    cdid_pipeline_metadata_abort(
      "Required arguments (`panel_owner`, `targets`, `anticipation`, `control_group`, `conventions`, `control`) must be supplied explicitly.",
      c("cdid_pipeline_metadata_input_error", "cdid_catalog_input_error")
    )
  }
  if (!gmm_scalar_logical(build_design)) {
    cdid_pipeline_metadata_abort("`build_design` must be a scalar logical.", "cdid_pipeline_metadata_input_error")
  }
  for (item in list(
    list(v = targets, nm = "targets"),
    list(v = anticipation, nm = "anticipation"),
    list(v = control_group, nm = "control_group"),
    list(v = conventions, nm = "conventions"),
    list(v = requested_post, nm = "requested_post"),
    list(v = declared_qualified, nm = "declared_qualified"),
    list(v = target_post_times, nm = "target_post_times")
  )) {
    ipw_owner_refuse_decorative_attrs(item$v, item$nm)
  }

  owners <- pipeline_metadata_resolve_owners(panel_owner, frame_owner = frame_owner)
  p_own <- owners$panel_owner
  f_own <- owners$frame_owner
  pm_ctrl <- validate_pipeline_metadata_control(control, catalog_bounds)
  own_ctrl <- pm_ctrl$owner_control
  cat_bnd <- pm_ctrl$catalog_bounds

  occ_supplied <- !is.null(occupancy)
  if (occ_supplied) {
    validate_owner_sparse_occupancy(occupancy, panel_owner = p_own, frame_owner = f_own)
    occ_pre <- occupancy$preflight
    n_occ_low <- nrow(occupancy$occupancy)
    n_occ_high <- n_occ_low
  } else {
    occ_pre <- preflight_owner_sparse_occupancy(p_own, own_ctrl, frame_owner = f_own)
    n_occ_low <- 0
    n_occ_high <- occ_pre$max_occupancy_cells
  }

  cat_early <- pipeline_metadata_validate_catalog_inputs_early(
    calendar = p_own$all_times,
    cohorts = p_own$cohorts,
    targets = targets,
    anticipation = anticipation,
    control_group = control_group,
    conventions = conventions,
    requested_post = requested_post,
    catalog_bounds = cat_bnd,
    n_occ_lower = n_occ_low,
    n_occ_upper = n_occ_high,
    owner_control = own_ctrl,
    occupancy_retained = occ_pre$total_retained_entries,
    work_prefix = if(occ_supplied) 0 else occ_pre$work_ops_bound
  )

  # Mirror of `catalog_adapter_preflight` before any catalog row exists: the
  # same typed envelope on lower-bound dimensions (no rows, pools or operator
  # keys yet). It is a necessary condition only; the built catalog's actual
  # dimensions are checked again before any adapter allocation.
  if (isTRUE(build_design)) {
    J_ad <- length(cat_early$target_indices)
    lower <- list(R = 0, C = as.double(length(cat_early$spec$cohorts)), T = as.double(length(cat_early$spec$calendar)),
      P = 0, NP = 0, PR = 0, J = as.double(J_ad),
      D = if (is.data.frame(declared_qualified)) as.double(nrow(declared_qualified)) else 0,
      Q = if (is.data.frame(target_post_times)) as.double(nrow(target_post_times)) else
        catalog_product(J_ad, length(target_post_times), "requested design grid"),
      O = 0, B = cat_early$byte_conventions, Bf = 0)
    env_ad <- catalog_design_envelope(lower, design = !is.null(target_post_times))
    if (env_ad$scratch > cat_bnd$max_scratch || env_ad$work > cat_bnd$max_work) {
      catalog_abort(
        "Catalog adapter output/string/graph work exceeds its explicit stored controls before allocation.",
        "cdid_catalog_bound_error"
      )
    }
  }

  list(
    occupancy_supplied = occ_supplied,
    occupancy_preflight = occ_pre,
    catalog_header_preflight = list(
      n_calendar_dates = length(p_own$all_times),
      n_cohorts = length(p_own$cohorts),
      n_targets = length(cat_early$target_indices),
      n_requested_post = length(cat_early$requested$gi),
      convention_bytes = cat_early$byte_conventions,
      n_occupancy_lower = n_occ_low,
      n_occupancy_upper = n_occ_high,
      min_header_scratch = cat_early$min_header_scratch,
      max_header_scratch = cat_early$max_header_scratch,
      min_header_work = cat_early$min_header_work,
      max_header_work = cat_early$max_header_work,
      catalog_bounds = cat_bnd
    ),
    owner_control = own_ctrl,
    catalog_bounds = cat_bnd,
    build_design = isTRUE(build_design),
    endpoints_read = 0L,
    peak_rss_promised = FALSE
  )
}

preflight_pipeline_metadata <- preflight_owner_comparison_catalog

# Resolve explicit row-to-operator/pool/occupancy and cohort support mappings
# using exact integer dictionary indices (zero unit scans, zero float formatting).
pipeline_metadata_match_pair <- function(si, ti, s_keys, t_keys) {
  lo <- 1; hi <- length(s_keys)
  while (lo <= hi) {
    mid <- floor(lo + (hi - lo) / 2)
    if (s_keys[mid] == si && t_keys[mid] == ti) return(as.integer(mid))
    if (s_keys[mid] < si || (s_keys[mid] == si && t_keys[mid] < ti)) lo <- mid + 1 else hi <- mid - 1
  }
  NA_integer_
}

# Conditional success envelope for the unchanged catalog. Exact candidate count
# uses only already sealed occupancy indices; no outcome or candidate arrays.
pipeline_metadata_catalog_reservation <- function(panel_owner, occupancy_obj, targets,
  anticipation, control_group, conventions, requested_post, owner_control, catalog_bounds,
  work_prefix) {
  ci <- occupancy_obj$cell_indices
  T <- panel_owner$n_calendar_dates; C <- length(panel_owner$cohorts)
  J <- length(targets); O <- nrow(ci); Q <- if(is.data.frame(requested_post)) nrow(requested_post) else 0
  H <- catalog_size(c(T,C,J,O,Q,1), "composed catalog header")
  bytes <- as.double(nchar(conventions,type="bytes"))
  hdr_scratch <- catalog_size(c(catalog_product(64,H,"catalog header"),catalog_product(16,bytes,"convention bytes")),"catalog header scratch")
  hdr_work <- catalog_size(c(catalog_product(64*(1+ceiling(log2(H+1))),H,"catalog header visits"),16*bytes),"catalog header work")
  occ_ret <- occupancy_obj$preflight$total_retained_entries
  ipw_owner_safe_sum(c(occ_ret,hdr_scratch),owner_control$max_matrix_entries,"occupancy plus catalog reservation header")
  ipw_owner_safe_sum(c(work_prefix,2*hdr_work),owner_control$max_work_ops,"catalog reservation/header work")
  tg <- sort(match(as.double(targets),panel_owner$cohorts))
  refs <- integer(J)
  for(j in seq_len(J)) {
    lo<-1;hi<-T;answer<-0L
    while(lo<=hi) {
      mid<-floor(lo+(hi-lo)/2)
      if(exact_cmp_sum(panel_owner$all_times[mid],anticipation,panel_owner$cohorts[tg[j]])<0L) {answer<-as.integer(mid);lo<-mid+1} else hi<-mid-1
    }
    refs[j]<-answer
  }
  hfirst <- findInterval(panel_owner$cohorts[tg],panel_owner$all_times,left.open=TRUE)+1L
  ntp<-0;ncand<-0
  for(i in seq_len(O)) {
    j<-ipw_owner_match_sorted_double(ci$cohort_index[i],tg)
    if(is.na(j)) next
    ntp<-ntp+1
    if(refs[j]>0L) ncand<-catalog_size(c(ncand,1,if(identical(control_group,"notyettreated")) max(0,T-max(hfirst[j],ci$t_index[i])+1) else 0),"exact catalog candidate count")
  }
  # Bounds condition on existing catalog's own refusals, without changing them.
  R <- min(ncand,catalog_bounds$max_rows)
  E <- min(catalog_product(R,C,"row pool envelope"),catalog_bounds$max_pool_entries)
  stage<-catalog_size(c(ncand,ntp,1),"catalog candidate stage")
  key_slots<-catalog_product(catalog_product(16,bytes,"key bytes"),ncand,"catalog key payload")
  scratch<-catalog_size(c(hdr_scratch,9*O,T,J,8*ntp,catalog_product(256,stage,"catalog candidate buffers"),key_slots,5*ncand,catalog_product(64,E,"pool buffers"),13*R,4*J,4*Q),"catalog cumulative reservation")
  cat_work<-catalog_size(c(hdr_work,catalog_product(128*(C+1+ceiling(log2(stage+1))),stage,"candidate work"),key_slots),"catalog work reservation")
  map_upper<-catalog_size(c(occ_ret,scratch,43*R,256),"composed retained catalog/mapping envelope")
  ipw_owner_safe_sum(c(map_upper),min(owner_control$max_retained_entries,owner_control$max_output_entries),"composed catalog/mapping retention and output before catalog allocation")
  ipw_owner_safe_sum(c(occ_ret,scratch),owner_control$max_matrix_entries,"occupancy plus catalog workspace")
  # Retained/output are rechecked using actual rows before mapping allocation.
  whole_work<-ipw_owner_safe_sum(c(work_prefix,2*hdr_work,cat_work),owner_control$max_work_ops,"composed occupancy/catalog work")
  list(catalog_scratch_upper=scratch,catalog_work_upper=cat_work,catalog_candidate_count=ncand,
    catalog_target_pair_count=ntp,mapping_retained_upper=map_upper,work_prefix=whole_work)
}

pipeline_metadata_resolve_catalog_mappings <- function(
  panel_owner,
  frame_owner = NULL,
  occupancy_obj,
  catalog,
  owner_control,
  work_prefix
) {
  n_coh <- length(panel_owner$cohorts)
  n_cal <- panel_owner$n_calendar_dates
  n_occ <- nrow(occupancy_obj$cell_indices)
  n_rows <- nrow(catalog$rows)
  n_ops <- nrow(catalog$operators)
  n_pools <- catalog$pools$n_pools
  n_pool_entries <- length(catalog$pools$entries)

  # Preflight mapping retention, output, scratch, and work against owner_control BEFORE allocation
  map_retained <- ipw_owner_safe_sum(
    c(
      occ_retained <- occupancy_obj$preflight$total_retained_entries,
      catalog$resources$scratch_charged,
      ipw_owner_safe_product(30, n_rows, owner_control$max_retained_entries, "retained catalog row mapping columns"),
      ipw_owner_safe_product(13, n_ops, owner_control$max_retained_entries, "retained catalog operator mapping columns"),
      256
    ),
    max_bound = owner_control$max_retained_entries,
    label = "retained pipeline metadata mapping entries"
  )
  if (map_retained > owner_control$max_output_entries) {
    cdid_pipeline_metadata_abort(
      paste0("Pipeline metadata output entries (", map_retained, ") exceed `control$max_output_entries` (", owner_control$max_output_entries, ")."),
      c("cdid_pipeline_metadata_bound_error", "cdid_ipw_owner_bound_error"),
      diagnostics = list(requested = map_retained, bound = owner_control$max_output_entries)
    )
  }
  map_scratch <- ipw_owner_safe_sum(
    c(
      map_retained,
      ipw_owner_safe_product(16,n_occ,owner_control$max_matrix_entries,"CSR occupancy sort/index scratch"),
      ipw_owner_safe_product(12,n_coh,owner_control$max_matrix_entries,"live row/pool/cohort slices"),
      ipw_owner_safe_product(8,n_pools+n_pool_entries+1,owner_control$max_matrix_entries,"pool lookup scratch"),
      256
    ),
    max_bound = owner_control$max_matrix_entries,
    label = "concurrent pipeline metadata mapping entries"
  )
  row_pair_visits <- ipw_owner_safe_product(n_rows,n_coh,owner_control$max_work_ops,"repeated row/pair/cohort visits")
  map_work <- ipw_owner_safe_sum(c(work_prefix,
    ipw_owner_safe_product(32*(1+ceiling(log2(n_occ+1))+ceiling(log2(n_coh+1))),
      ipw_owner_safe_sum(c(n_occ,row_pair_visits,n_pool_entries,n_ops,n_coh,1)),
      owner_control$max_work_ops,"CSR sort/search and repeated row-pool work"),256),
    owner_control$max_work_ops,"composed occupancy/catalog/mapping work")

  # Map panel A-unit cohort counts (and optional frame weights) to full `panel_owner$cohorts`
  coh_a_count <- integer(n_coh)
  coh_a_weight <- if (!is.null(frame_owner)) numeric(n_coh) else rep(NA_real_, n_coh)
  if (length(panel_owner$a_cohorts) > 0L) {
    a_in_all <- as.integer(ipw_owner_match_sorted_double(panel_owner$a_cohorts, panel_owner$cohorts))
    coh_a_count[a_in_all] <- panel_owner$cohort_a_count
    if (!is.null(frame_owner)) {
      coh_a_weight[a_in_all] <- frame_owner$cohort_a_weight
    }
  }

  # Per-pool A-unit count and weight mass
  pl <- catalog$pools
  pool_a_count <- numeric(n_pools)
  pool_a_weight <- if (!is.null(frame_owner)) numeric(n_pools) else rep(NA_real_, n_pools)
  for (p in seq_len(n_pools)) {
    p_idx <- pl$entries[seq.int(pl$start[p], length.out = pl$size[p])]
    pool_a_count[p] <- ipw_owner_safe_sum(coh_a_count[p_idx], label = "pool A-unit count")
    if (!is.null(frame_owner)) {
      pool_a_weight[p] <- sum(coh_a_weight[p_idx])
    }
  }

  # Build CSR index on `occupancy_obj$cell_indices` ordered by (s_index, t_index, cohort_index)
  ci_df <- occupancy_obj$cell_indices
  if (n_occ > 0L) {
    ord_pair <- order(ci_df$s_index, ci_df$t_index, ci_df$cohort_index, method = "radix")
    p_cell_id <- ci_df$cell_id[ord_pair]
    p_ki <- ci_df$cohort_index[ord_pair]
    p_si <- ci_df$s_index[ord_pair]
    p_ti <- ci_df$t_index[ord_pair]
    p_cnt <- ci_df$count[ord_pair]
    new_p <- if (n_occ > 1L) c(TRUE, p_si[-1L] != p_si[-n_occ] | p_ti[-1L] != p_ti[-n_occ]) else TRUE
    pair_starts <- which(new_p)
    pair_ends <- c(pair_starts[-1L] - 1L, n_occ)
    pair_s_key <- p_si[pair_starts]
    pair_t_key <- p_ti[pair_starts]
    # Keep two exact integer keys; their packed product could exceed2^53.
  } else {
    p_cell_id <- integer(0L)
    p_ki <- integer(0L)
    p_cnt <- numeric(0L)
    pair_starts <- integer(0L)
    pair_ends <- integer(0L)
    pair_s_key <- pair_t_key <- integer(0L)
  }

  rows <- catalog$rows
  target_occ_cell_id <- integer(n_rows)
  target_pair_count <- numeric(n_rows)
  control_pool_pair_count <- numeric(n_rows)
  risk_pair_count <- numeric(n_rows)
  all_cohort_pair_count <- numeric(n_rows)
  n_supported_controls <- integer(n_rows)

  for (r in seq_len(n_rows)) {
    gi_r <- rows$target_index[r]
    si_r <- rows$s_index[r]
    ti_r <- rows$t_index[r]
    pid_r <- rows$pool_id[r]
    p_pos <- pipeline_metadata_match_pair(si_r,ti_r,pair_s_key,pair_t_key)
    if (is.na(p_pos)) {
      cdid_pipeline_metadata_abort(
        "Internal catalog-to-occupancy invariant failed: retained catalog row pair is absent from occupancy.",
        "cdid_pipeline_metadata_internal_error"
      )
    }
    sl <- seq.int(pair_starts[p_pos], pair_ends[p_pos])
    sl_ki <- p_ki[sl]
    sl_cnt <- p_cnt[sl]
    sl_cid <- p_cell_id[sl]

    t_hit <- ipw_owner_match_sorted_double(gi_r, sl_ki)
    if (is.na(t_hit) || sl_cnt[t_hit] < 1) {
      cdid_pipeline_metadata_abort(
        "Internal catalog-to-occupancy invariant failed: retained target cohort lacks positive pair occupancy.",
        "cdid_pipeline_metadata_internal_error"
      )
    }
    target_occ_cell_id[r] <- sl_cid[t_hit]
    target_pair_count[r] <- sl_cnt[t_hit]

    pool_k <- pl$entries[seq.int(pl$start[pid_r], length.out = pl$size[pid_r])]
    c_hits <- ipw_owner_match_sorted_double(pool_k, sl_ki)
    if (anyNA(c_hits) || any(sl_cnt[c_hits] < 1)) {
      cdid_pipeline_metadata_abort(
        "Internal catalog-to-occupancy invariant failed: retained control cohort lacks positive pair occupancy.",
        "cdid_pipeline_metadata_internal_error"
      )
    }
    n_supported_controls[r] <- length(pool_k)
    control_pool_pair_count[r] <- ipw_owner_safe_sum(sl_cnt[c_hits], label = "control pool pair count")
    risk_pair_count[r] <- ipw_owner_safe_sum(
      c(target_pair_count[r], control_pool_pair_count[r]),
      label = "risk-set pair count"
    )
    all_cohort_pair_count[r] <- ipw_owner_safe_sum(sl_cnt, label = "all-cohort pair count")
  }

  target_a_unit_count_row <- as.double(coh_a_count[rows$target_index])
  control_a_unit_count_row <- if (n_rows > 0L) pool_a_count[rows$pool_id] else numeric(0L)
  risk_a_unit_count_row <- target_a_unit_count_row + control_a_unit_count_row

  row_mapping <- data.frame(
    row_id = rows$row_id,
    operator_id = rows$operator_id,
    pool_id = rows$pool_id,
    pool_size = rows$pool_size,
    target_index = rows$target_index,
    g = rows$g,
    s_index = rows$s_index,
    t_index = rows$t_index,
    s = rows$s,
    t = rows$t,
    ref_index = rows$ref_index,
    ref = panel_owner$all_times[rows$ref_index],
    functional_key = catalog$operators$functional_key[rows$operator_id],
    s_role = rows$s_role,
    t_role = rows$t_role,
    s_coordinate_id = rows$s_coordinate_id,
    t_coordinate_id = rows$t_coordinate_id,
    original = rows$original,
    n_horizons = rows$n_horizons,
    horizon_first_index = rows$horizon_first_index,
    horizon_last_index = rows$horizon_last_index,
    target_occupancy_cell_id = target_occ_cell_id,
    target_pair_count = target_pair_count,
    control_pool_pair_count = control_pool_pair_count,
    risk_pair_count = risk_pair_count,
    all_cohort_pair_count = all_cohort_pair_count,
    target_a_unit_count = target_a_unit_count_row,
    control_pool_a_unit_count = control_a_unit_count_row,
    risk_a_unit_count = risk_a_unit_count_row,
    n_supported_controls = n_supported_controls,
    stringsAsFactors = FALSE
  )

  ops <- catalog$operators
  op_row_counts <- if (n_ops > 0L) tabulate(rows$operator_id, nbins = n_ops) else integer(0L)
  op_target_a_count <- as.double(coh_a_count[ops$target_index])
  op_control_a_count <- if (n_ops > 0L) pool_a_count[ops$pool_id] else numeric(0L)
  op_risk_a_count <- op_target_a_count + op_control_a_count
  op_target_a_weight <- if (!is.null(frame_owner)) coh_a_weight[ops$target_index] else rep(NA_real_, n_ops)
  op_control_a_weight <- if (!is.null(frame_owner) && n_ops > 0L) pool_a_weight[ops$pool_id] else rep(NA_real_, n_ops)

  operator_mapping <- data.frame(
    operator_id = ops$operator_id,
    target_index = ops$target_index,
    g = ops$g,
    pool_id = ops$pool_id,
    pool_size = if (n_ops > 0L) pl$size[ops$pool_id] else integer(0L),
    never_treated_in_pool = if (n_ops > 0L) pl$never_treated[ops$pool_id] else logical(0L),
    functional_key = ops$functional_key,
    n_catalog_rows = op_row_counts,
    target_a_unit_count = op_target_a_count,
    control_pool_a_unit_count = op_control_a_count,
    risk_a_unit_count = op_risk_a_count,
    target_a_weight = op_target_a_weight,
    control_pool_a_weight = op_control_a_weight,
    stringsAsFactors = FALSE
  )

  list(
    row_mapping = row_mapping,
    operator_mapping = operator_mapping,
    accounting = list(
      retained_entries = map_retained,
      concurrent_entries = map_scratch,
      work_ops_bound = map_work
    )
  )
}

# CP08 compact basis (catalog reuse): binds one already constructed, sealed
# `cdid_pipeline_metadata` object as the catalog/occupancy/mapping source of a
# second construction over the same owners. Every predicate is a check on the
# actual sealed snapshot and the exact original argument objects: the source
# must pass the full metadata validator against the live owners, name the same
# sealed panel/frame/occupancy references, carry its own construction record
# whose arguments are `identical()` to the ones now supplied (no label, no
# formatted-number key), and have been built under the identical validated
# control. Anything else refuses with the original typed metadata errors or a
# reuse refusal; nothing is derived from the source when it is refused.
pipeline_metadata_bind_catalog_source <- function(source, panel_owner, frame_owner, pm_ctrl, occupancy,
                                                  targets, anticipation, control_group, conventions,
                                                  requested_post, pt_assumption, pt_start) {
  bad <- function(msg, extra = character()) cdid_pipeline_metadata_abort(msg,
    c(extra, "cdid_pipeline_metadata_reuse_error", "cdid_pipeline_metadata_input_error"))
  # Original seal/owner/occupancy/catalog/mapping validator and its typed refusals first.
  validate_owner_comparison_catalog(source, panel_owner = panel_owner, frame_owner = frame_owner)
  if (!identical(source$panel_ref, panel_owner$ref)) {
    bad("Catalog source metadata was built on a different sealed panel owner; reuse is refused.",
      c("cdid_pipeline_metadata_owner_mismatch_error", "cdid_ipw_owner_mismatch_error"))
  }
  f_ref <- if (is.null(frame_owner)) NULL else frame_owner$ref
  if (!identical(source$frame_ref, f_ref)) {
    bad("Catalog source metadata was built on a different sealed frame owner; reuse is refused.",
      c("cdid_pipeline_metadata_owner_mismatch_error", "cdid_ipw_owner_mismatch_error"))
  }
  if (!is.null(occupancy) && (!is.list(occupancy) || !identical(occupancy$ref, source$occupancy_ref) ||
      !identical(source$occupancy_object$ref, source$occupancy_ref))) {
    bad("Supplied occupancy is not the sealed occupancy the catalog source was built from; reuse is refused.",
      "cdid_pipeline_metadata_owner_mismatch_error")
  }
  cons <- source$construction
  cons_names <- c("targets", "anticipation", "control_group", "conventions", "requested_post", "pt_assumption", "pt_start")
  if (!is.list(cons) || is.object(cons) || !identical(names(cons), cons_names)) {
    bad("Catalog source metadata carries no construction record of its exact original arguments; reuse is refused.")
  }
  now <- list(targets = targets, anticipation = anticipation, control_group = control_group,
    conventions = conventions, requested_post = requested_post, pt_assumption = pt_assumption, pt_start = pt_start)
  for (nm in cons_names) {
    if (!identical(cons[[nm]], now[[nm]])) {
      bad(paste0("Catalog source metadata was constructed with a different `", nm,
        "`; the catalog is not reused for different arguments."), "cdid_catalog_input_error")
    }
  }
  if (!identical(source$control, pm_ctrl)) {
    bad("Catalog source metadata was constructed under a different validated control; reuse is refused.",
      "cdid_pipeline_metadata_control_error")
  }
  acc <- source$diagnostics$mapping_accounting
  if (!is.list(acc) || !gmm_scalar_numeric(acc$retained_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(acc$concurrent_entries, nonneg = TRUE) || !gmm_scalar_numeric(acc$work_ops_bound, nonneg = TRUE)) {
    bad("Catalog source metadata lacks its mapping accounting; reuse is refused.")
  }
  if (!is.data.frame(source$operator_mapping) || !inherits(source$catalog, "cdid_comparison_catalog") ||
      !is.data.frame(source$catalog$operators) || nrow(source$operator_mapping) != nrow(source$catalog$operators)) {
    bad("Catalog source metadata operator mapping is malformed; reuse is refused.", "cdid_pipeline_metadata_integrity_error")
  }
  list(
    occupancy_object = source$occupancy_object,
    catalog = source$catalog,
    mappings = list(row_mapping = source$row_mapping, operator_mapping = source$operator_mapping, accounting = acc),
    reservation = source$preflight$catalog_composition_reservation,
    source_ref = source$ref,
    binding = list(panel_ref = source$panel_ref, frame_ref = source$frame_ref, occupancy_ref = source$occupancy_ref,
      source_metadata_ref = source$ref, arguments_identical = cons_names, control_identical = TRUE)
  )
}

build_owner_comparison_catalog <- function(
  panel_owner,
  targets,
  anticipation,
  control_group,
  conventions,
  requested_post = NULL,
  control,
  catalog_bounds = NULL,
  occupancy = NULL,
  frame_owner = NULL,
  declared_qualified = NULL,
  target_post_times = NULL,
  build_design = FALSE,
  trace_env = NULL,
  metadata_id = NULL,
  pt_assumption = "post", pt_start = NULL,
  metadata_source = NULL
) {
  # Stage-local demand: owners/controls, the occupancy and catalog stages this
  # construction runs, and the design bridge only when a design is requested.
  pipeline_metadata_require_stage(c("owners", "occupancy", "catalog", if (isTRUE(build_design)) "design"))
  if (missing(panel_owner) || missing(targets) || missing(anticipation) ||
      missing(control_group) || missing(conventions) || missing(control)) {
    cdid_pipeline_metadata_abort(
      "Required arguments (`panel_owner`, `targets`, `anticipation`, `control_group`, `conventions`, `control`) must be supplied explicitly.",
      c("cdid_pipeline_metadata_input_error", "cdid_catalog_input_error")
    )
  }
  if (!is.null(metadata_id) && !gmm_scalar_character(metadata_id)) {
    cdid_pipeline_metadata_abort(
      "`metadata_id` must be NULL or a nonempty character scalar.",
      "cdid_pipeline_metadata_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_pipeline_metadata_abort(
      "`trace_env` must be NULL or an environment.",
      "cdid_pipeline_metadata_input_error"
    )
  }

  owners <- pipeline_metadata_resolve_owners(panel_owner, frame_owner = frame_owner)
  p_own <- owners$panel_owner
  validate_chained_pt_request(pt_assumption, pt_start, p_own$all_times)
  f_own <- owners$frame_owner
  pm_ctrl <- validate_pipeline_metadata_control(control, catalog_bounds)
  own_ctrl <- pm_ctrl$owner_control
  cat_bnd <- pm_ctrl$catalog_bounds

  # CP08 compact basis: an explicitly supplied sealed metadata source of the
  # same owners and exact arguments may be reused below. Its occupancy stands
  # in for an absent `occupancy` argument so the unchanged stage 1 preflight
  # validates the same object; the source itself is bound after that preflight,
  # so malformed arguments keep their original typed refusals and order.
  # Standalone calls (NULL) construct everything exactly as before.
  reuse <- NULL
  if (!is.null(metadata_source) && is.null(occupancy) && is.list(metadata_source)) {
    occupancy <- metadata_source$occupancy_object
  }

  # 1. Stage 1 Preflight BEFORE any pair-position reader traversal or callbacks
  pre <- preflight_owner_comparison_catalog(
    panel_owner = p_own,
    targets = targets,
    anticipation = anticipation,
    control_group = control_group,
    conventions = conventions,
    requested_post = requested_post,
    control = own_ctrl,
    catalog_bounds = cat_bnd,
    occupancy = occupancy,
    frame_owner = f_own,
    declared_qualified = declared_qualified,
    target_post_times = target_post_times,
    build_design = build_design
  )
  if (!is.null(metadata_source)) {
    reuse <- pipeline_metadata_bind_catalog_source(metadata_source, p_own, f_own, pm_ctrl, occupancy,
      targets, anticipation, control_group, conventions, requested_post, pt_assumption, pt_start)
  }

  # 2. Derive or validate sealed sparse `(cohort, s, t, count)` occupancy
  occ_obj <- if (is.null(occupancy)) {
    derive_owner_sparse_occupancy(
      panel_owner = p_own,
      control = own_ctrl,
      frame_owner = f_own,
      trace_env = trace_env
    )
  } else {
    ipw_owner_reset_trace(trace_env)
    validate_owner_sparse_occupancy(occupancy, panel_owner = p_own, frame_owner = f_own)
  }

  if (is.null(reuse)) {
    # Charge both already-completed prefix and the next catalog workspace before allocation.
    reservation <- pipeline_metadata_catalog_reservation(p_own,occ_obj,targets,anticipation,
      control_group,conventions,requested_post,own_ctrl,cat_bnd,
      work_prefix = (if(isTRUE(pre$occupancy_supplied)) 0 else occ_obj$preflight$work_ops_bound) +
        pre$catalog_header_preflight$max_header_work)

    # 3. Call unchanged `build_comparison_catalog` with literal stored dates, cohorts, and anticipation
    catalog <- build_comparison_catalog(
      calendar = p_own$all_times,
      cohorts = p_own$cohorts,
      targets = targets,
      anticipation = anticipation,
      occupancy = occ_obj$occupancy,
      control_group = control_group,
      conventions = conventions,
      requested_post = requested_post,
      bounds = cat_bnd
    )

    # 4. Resolve explicit row-to-operator/cohort/occupancy mappings
    mappings <- pipeline_metadata_resolve_catalog_mappings(
      panel_owner = p_own,
      frame_owner = f_own,
      occupancy_obj = occ_obj,
      catalog = catalog,
      owner_control = own_ctrl,
      work_prefix = reservation$work_prefix
    )
  } else {
    # 3./4. Reused: the sealed source's catalog and mappings (the exact aliased
    # snapshot its seal binds), its composed reservation borrowed; no catalog
    # row, CSR mapping or candidate buffer is rebuilt. Only the design step
    # below is constructed and charged here.
    if (!identical(occ_obj$ref, reuse$occupancy_object$ref)) {
      cdid_pipeline_metadata_abort("Validated occupancy differs from the catalog source's sealed occupancy.",
        c("cdid_pipeline_metadata_reuse_error", "cdid_pipeline_metadata_owner_mismatch_error"))
    }
    reservation <- reuse$reservation
    catalog <- reuse$catalog
    mappings <- reuse$mappings
  }

  # 5. Optional design-link and calendar-design bridge (only when explicitly requested).
  #    CP08 validated catalog-to-design handoff: one adapter preflight on the
  #    actual catalog object supplies the typed dimensions of both envelopes
  #    (the links stage and, with requested targets, the complete design; the
  #    links envelope does not depend on the target count), charged exactly as
  #    before. The declared-qualification packet is then validated once and its
  #    links, canonical calendar and typed pools enter the private design core,
  #    which performs the one complete link validation and the unchanged graph
  #    construction. Standalone constructors keep their own strict wrappers.
  design_links <- NULL
  design_obj <- NULL
  design_accounting <- list(owned_entries = 0, work_ops = 0)
  if (isTRUE(build_design)) {
    ad <- catalog_adapter_preflight(catalog,declared_qualified,target_post_times,design=!is.null(target_post_times))
    ad_links <- if(!is.null(target_post_times)) catalog_design_envelope(ad$dimensions, design = FALSE) else list(scratch=0,work=0)
    design_scratch <- ipw_owner_safe_sum(c(ad$scratch,ad_links$scratch))
    ipw_owner_safe_sum(c(mappings$accounting$retained_entries,design_scratch),min(own_ctrl$max_retained_entries,own_ctrl$max_output_entries),"retained catalog/mappings plus design output")
    ipw_owner_safe_sum(c(mappings$accounting$concurrent_entries,design_scratch),own_ctrl$max_matrix_entries,"catalog/mappings plus design workspace")
    design_accounting <- list(owned_entries = design_scratch,
      work_ops = ipw_owner_safe_sum(c(mappings$accounting$work_ops_bound,ad$work,ad_links$work),own_ctrl$max_work_ops,"composed metadata/design work") -
        mappings$accounting$work_ops_bound)
    handoff <- catalog_design_handoff(catalog, declared_qualified)
    if (!is.null(target_post_times)) {
      design_obj <- build_catalog_handoff_design(
        handoff,
        target_post_times = target_post_times,
        pt_assumption = pt_assumption, pt_start = pt_start
      )
      design_links <- design_obj$links
    } else {
      design_links <- handoff$links
    }
    handoff <- NULL
  }

  meta_obj <- structure(
    list(
      metadata_id = metadata_id,
      panel_ref = p_own$ref,
      frame_ref = if (!is.null(f_own)) f_own$ref else NULL,
      occupancy_ref = occ_obj$ref,
      occupancy = occ_obj$occupancy,
      cell_indices = occ_obj$cell_indices,
      occupancy_object = occ_obj,
      catalog = catalog,
      row_mapping = mappings$row_mapping,
      operator_mapping = mappings$operator_mapping,
      design_links = design_links,
      design = design_obj,
      owner_summary = pipeline_metadata_owner_summary(p_own,f_own),
      control = pm_ctrl,
      # Exact original argument objects of this construction (bound by the seal);
      # a later construction over the same owners may reuse this catalog only
      # when its own arguments are `identical()` to these.
      construction = list(targets = targets, anticipation = anticipation, control_group = control_group,
        conventions = conventions, requested_post = requested_post, pt_assumption = pt_assumption, pt_start = pt_start),
      catalog_source_ref = if (is.null(reuse)) NULL else reuse$source_ref,
      preflight = c(pre,list(catalog_composition_reservation=reservation,
        catalog_reuse_binding = if (is.null(reuse)) NULL else reuse$binding)),
      diagnostics = list(
        # Owned/borrowed accounting of this construction: a reused catalog and
        # mapping are borrowed from the sealed source (their retained entries
        # and work were charged once, there); only the design step is owned here.
        catalog_reuse = list(reused = !is.null(reuse),
          catalog_constructed = is.null(reuse), mappings_constructed = is.null(reuse),
          borrowed_catalog_mapping_entries = if (is.null(reuse)) 0 else mappings$accounting$retained_entries,
          owned_catalog_mapping_entries = if (is.null(reuse)) mappings$accounting$retained_entries else 0,
          owned_design_entries = design_accounting$owned_entries,
          owned_work_ops = if (is.null(reuse)) mappings$accounting$work_ops_bound + design_accounting$work_ops else design_accounting$work_ops,
          borrowed_work_ops = if (is.null(reuse)) 0 else mappings$accounting$work_ops_bound,
          source_metadata_ref = if (is.null(reuse)) NULL else reuse$source_ref),
        endpoints_read = 0L,
        read_outcomes = FALSE,
        global_pair_subtraction_performed = FALSE,
        pair_differences_computed = FALSE,
        all_unit_pair_table_materialized = FALSE,
        n_by_moment_cache_allocated = FALSE,
        formatted_float_keys_used = FALSE,
        dictionary_index_keys_used = TRUE,
        all_catalog_rows_retained = identical(nrow(mappings$row_mapping), nrow(catalog$rows)),
        all_catalog_operators_retained = identical(nrow(mappings$operator_mapping), nrow(catalog$operators)),
        singleton_metadata_preserved_in_owner = TRUE,
        unobserved_metadata_preserved_in_owner = TRUE,
        explicit_na_metadata_preserved_in_owner = TRUE,
        occupancy_diagnostics = occ_obj$diagnostics,
        catalog_resources = catalog$resources,
        mapping_accounting = mappings$accounting
      )
    ),
    class = "cdid_pipeline_metadata"
  )

  ipw_owner_seal(meta_obj, "pipeline_metadata")
}

adapt_owner_to_comparison_catalog <- build_owner_comparison_catalog
construct_pipeline_metadata <- build_owner_comparison_catalog
chained_pipeline_metadata <- build_owner_comparison_catalog

validate_owner_comparison_catalog <- function(metadata_obj, panel_owner = NULL, frame_owner = NULL) {
  pipeline_metadata_require_stage("owners")
  if (missing(metadata_obj) || is.null(metadata_obj) ||
      !inherits(metadata_obj, "cdid_pipeline_metadata") ||
      !is.list(metadata_obj) || !ipw_owner_seal_ok(metadata_obj, "pipeline_metadata")) {
    cdid_pipeline_metadata_abort(
      "Pipeline metadata object is malformed, foreign, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_pipeline_metadata_integrity_error", "cdid_ipw_owner_integrity_error")
    )
  }
  validate_owner_sparse_occupancy(metadata_obj$occupancy_object, panel_owner = panel_owner, frame_owner = frame_owner)
  p_obj <- get("panel_owner",envir=metadata_obj$panel_ref,inherits=FALSE)
  if(!is.null(metadata_obj$frame_ref)) {
    f_obj <- get("frame_owner",envir=metadata_obj$frame_ref,inherits=FALSE)
    validate_ipw_frame_owner(f_obj,panel_owner=p_obj)
    if(!is.null(frame_owner) && !identical(frame_owner$ref,metadata_obj$frame_ref))
      cdid_pipeline_metadata_abort("Metadata belongs to a different sealed frame.",c("cdid_pipeline_metadata_owner_mismatch_error","cdid_ipw_owner_mismatch_error"))
  } else if(!is.null(frame_owner)) validate_ipw_frame_owner(frame_owner,panel_owner=p_obj)
  if (!identical(metadata_obj$panel_ref, metadata_obj$occupancy_object$panel_ref) ||
      !identical(metadata_obj$occupancy_ref, metadata_obj$occupancy_object$ref)) {
    cdid_pipeline_metadata_abort(
      "Pipeline metadata `panel_ref` or `occupancy_ref` disagrees with its sealed occupancy object.",
      c("cdid_pipeline_metadata_integrity_error", "cdid_pipeline_metadata_owner_mismatch_error")
    )
  }
  if (!inherits(metadata_obj$catalog, "cdid_comparison_catalog") ||
      !is.data.frame(metadata_obj$catalog$rows) ||
      !is.data.frame(metadata_obj$row_mapping) ||
      nrow(metadata_obj$row_mapping) != nrow(metadata_obj$catalog$rows)) {
    cdid_pipeline_metadata_abort(
      "Pipeline metadata catalog or row mapping is malformed.",
      "cdid_pipeline_metadata_integrity_error"
    )
  }
  invisible(metadata_obj)
}

validate_pipeline_metadata <- validate_owner_comparison_catalog

# ---------------------------------------------------------------------------
# Exact dictionary-index row specification & compact outcome-free view helpers
# ---------------------------------------------------------------------------

pipeline_metadata_row_spec_preflight <- function(metadata_obj, row_id, control) {
  ctrl <- validate_ipw_owner_control(control)
  r <- ipw_owner_int_bound(row_id, "row_id", min_val=1)
  if(r > nrow(metadata_obj$catalog$rows)) cdid_pipeline_metadata_abort("Row index exceeds retained catalog.", "cdid_catalog_input_error")
  pid <- metadata_obj$catalog$rows$pool_id[r]
  k <- metadata_obj$catalog$pools$size[pid]
  slots <- ipw_owner_safe_sum(c(128,ipw_owner_safe_product(4,k)),min(ctrl$max_matrix_entries,ctrl$max_output_entries),"bounded row descriptor/pool slices")
  ipw_owner_safe_sum(c(128,ipw_owner_safe_product(16,k)),ctrl$max_work_ops,"bounded row descriptor work")
  invisible(slots)
}

extract_pipeline_catalog_row_spec <- function(metadata_obj, row_id) {
  validate_owner_comparison_catalog(metadata_obj)
  r <- ipw_owner_int_bound(row_id, "row_id", min_val = 1)
  n_rows <- nrow(metadata_obj$catalog$rows)
  if (r > n_rows) {
    cdid_pipeline_metadata_abort(
      paste0("`row_id` (", r, ") exceeds retained catalog row count (", n_rows, ")."),
      c("cdid_pipeline_metadata_bound_error", "cdid_catalog_input_error")
    )
  }
  pipeline_metadata_row_spec_preflight(metadata_obj,r,metadata_obj$control$owner_control)
  pipeline_metadata_row_spec_fields(metadata_obj, r)
}

# Descriptors of the retained catalog rows `rows`, in that order, under the
# single extractor's contract: the whole catalog is validated ONCE, then every
# row keeps its original bound check, preflight and failure order before its
# fields are read through one invocation-local descriptor source, created and
# dropped inside this call (CP08 construction; the registry no longer validates
# the whole catalog once per row).
extract_pipeline_catalog_row_specs <- function(metadata_obj, rows = seq_len(nrow(metadata_obj$catalog$rows)),
                                               retained_entries = 0, external_work_ops = 0,
                                               control = NULL, return_accounting = FALSE) {
  validate_owner_comparison_catalog(metadata_obj)
  n_rows <- nrow(metadata_obj$catalog$rows)
  source <- NULL
  accounting <- NULL
  out <- list()
  for (j in seq_along(rows)) {
    r <- ipw_owner_int_bound(rows[j], "row_id", min_val = 1)
    if (r > n_rows) {
      cdid_pipeline_metadata_abort(
        paste0("`row_id` (", r, ") exceeds retained catalog row count (", n_rows, ")."),
        c("cdid_pipeline_metadata_bound_error", "cdid_catalog_input_error")
      )
    }
    pipeline_metadata_row_spec_preflight(metadata_obj, r, metadata_obj$control$owner_control)
    if (is.null(source)) {
      source <- pipeline_metadata_descriptor_source(metadata_obj, control,
        retained_entries, external_work_ops, length(rows), include_output = TRUE)
      accounting <- source$accounting
      out <- vector("list", length(rows))
    }
    out[[j]] <- pipeline_metadata_row_spec_fields(metadata_obj, r, source)
  }
  source <- NULL
  if (isTRUE(return_accounting)) list(specs = out, accounting = accounting) else out
}

# Counts only existing dictionary lengths before any cohort vector is expanded.
# Two cached pool vectors use at most 2E slots over all pools; the 4E scratch
# allowance covers current slices/coercions beside them. Two indexed lists and
# their replacement copies use at most 8P slots; 256 covers the environment,
# bindings, scalar plan/list headers and batch wrapper. No new string keys or
# growing named-list prefixes are created. The optional 64-slot IPW memo is
# additional: the original 128/row block remains committed to its old records.
pipeline_metadata_descriptor_preflight <- function(metadata_obj, control = NULL,
                                                   retained_entries = 0, external_work_ops = 0,
                                                   n_descriptors = 0, include_output = FALSE,
                                                   memo_entries = 0) {
  ctrl <- validate_ipw_owner_control(if (is.null(control)) metadata_obj$control$owner_control else control)
  retained <- ipw_owner_int_bound(retained_entries, "retained_entries", min_val = 0)
  prefix <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
  nr <- ipw_owner_int_bound(n_descriptors, "n_descriptors", min_val = 0)
  memo <- ipw_owner_int_bound(memo_entries, "memo_entries", min_val = 0)
  P <- length(metadata_obj$catalog$pools$size)
  E <- length(metadata_obj$catalog$pools$entries)
  lim <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  context <- ipw_owner_safe_sum(c(256, memo,
    ipw_owner_safe_product(4, P, lim, "descriptor indexed pool lists"),
    ipw_owner_safe_product(2, E, lim, "descriptor cached cohort vectors")), lim,
    "descriptor context retained entries")
  scratch <- ipw_owner_safe_sum(c(256,
    ipw_owner_safe_product(4, P, ctrl$max_matrix_entries, "descriptor pool-list replacement copies"),
    ipw_owner_safe_product(4, E, ctrl$max_matrix_entries, "descriptor cohort slice and coercion scratch")),
    ctrl$max_matrix_entries, "descriptor context scratch")
  output <- if (isTRUE(include_output)) ipw_owner_safe_product(128, nr,
    min(ctrl$max_retained_entries, ctrl$max_output_entries), "descriptor output records") else 0
  # The caller's old envelope already includes its old descriptor records.
  # max avoids counting that original reservation twice; the new context joins it.
  live_base <- max(retained, output)
  ipw_owner_safe_sum(c(live_base, context), ctrl$max_retained_entries,
    "caller plus descriptor context retention")
  ipw_owner_safe_sum(c(live_base, context, scratch), ctrl$max_matrix_entries,
    "caller plus descriptor context and scratch")
  work <- ipw_owner_safe_sum(c(512,
    ipw_owner_safe_product(128, ipw_owner_safe_sum(c(nr, P, E, 1), ctrl$max_work_ops,
      "descriptor preparation dimensions"), ctrl$max_work_ops, "descriptor context preparation work"),
    ipw_owner_safe_product(16 * as.double(nr), P, ctrl$max_work_ops,
      "descriptor indexed-list replacement work")), ctrl$max_work_ops, "descriptor context whole-batch work")
  ipw_owner_safe_sum(c(prefix, work), ctrl$max_work_ops, "caller plus descriptor preparation work")
  list(retained_entries = context, scratch_entries = scratch, output_entries = output,
    work_ops_bound = work, n_pools = P, n_pool_entries = E, n_descriptors = nr)
}

# Invocation-local, pre-admitted aliases of one already validated catalog.
# Indexed pool slots avoid string keys; actual pool vectors are expanded lazily
# after the unchanged row preflight. This source establishes no trust itself.
pipeline_metadata_descriptor_source <- function(metadata_obj, control = NULL,
                                                retained_entries = 0, external_work_ops = 0,
                                                n_descriptors = nrow(metadata_obj$catalog$rows),
                                                include_output = FALSE, memo_entries = 0) {
  accounting <- pipeline_metadata_descriptor_preflight(metadata_obj, control,
    retained_entries, external_work_ops, n_descriptors, include_output, memo_entries)
  source <- new.env(parent = emptyenv())
  source$rows <- metadata_obj$catalog$rows
  source$map <- metadata_obj$row_mapping
  source$pools <- metadata_obj$catalog$pools
  source$cohorts <- metadata_obj$catalog$spec$cohorts
  source$calendar <- metadata_obj$catalog$spec$calendar
  source$K_indices <- vector("list", accounting$n_pools)
  source$K <- vector("list", accounting$n_pools)
  source$accounting <- accounting
  source
}

# Field reader shared by the validated extractor above and the construction
# loops, which validate the same catalog at their entry. It reads one row of
# an already validated, preflighted catalog from its existing column vectors
# (no one-row data frames) and establishes no trust itself. With a descriptor
# source of the same catalog, the row's pool K indices/values are the source's
# shared expansion.
pipeline_metadata_row_spec_fields <- function(metadata_obj, r, source = NULL) {
  ri <- as.integer(r)
  if (is.null(source)) {
    rows <- metadata_obj$catalog$rows
    map <- metadata_obj$row_mapping
    pl <- metadata_obj$catalog$pools
    coh <- metadata_obj$catalog$spec$cohorts
    cal <- metadata_obj$catalog$spec$calendar
    pid <- rows$pool_id[ri]
    k_indices <- pl$entries[seq.int(pl$start[pid], length.out = pl$size[pid])]
    K_vals <- coh[k_indices]
    k_indices <- as.integer(k_indices)
    K_vals <- as.double(K_vals)
  } else {
    rows <- source$rows
    map <- source$map
    coh <- source$cohorts
    cal <- source$calendar
    pid <- rows$pool_id[ri]
    k_indices <- source$K_indices[[pid]]
    if (is.null(k_indices)) {
      pl <- source$pools
      k_indices <- pl$entries[seq.int(pl$start[pid], length.out = pl$size[pid])]
      K_vals <- as.double(coh[k_indices])
      k_indices <- as.integer(k_indices)
      source$K_indices[[pid]] <- k_indices
      source$K[[pid]] <- K_vals
    } else {
      K_vals <- source$K[[pid]]
    }
  }

  list(
    row_id = ri,
    operator_id = rows$operator_id[ri],
    pool_id = pid,
    target_index = rows$target_index[ri],
    g = coh[rows$target_index[ri]],
    K_indices = k_indices,
    K = K_vals,
    s_index = rows$s_index[ri],
    t_index = rows$t_index[ri],
    s = cal[rows$s_index[ri]],
    t = cal[rows$t_index[ri]],
    ref_index = rows$ref_index[ri],
    ref = cal[rows$ref_index[ri]],
    functional_key = map$functional_key[ri],
    s_role = rows$s_role[ri],
    t_role = rows$t_role[ri],
    s_coordinate_id = rows$s_coordinate_id[ri],
    t_coordinate_id = rows$t_coordinate_id[ri],
    original = rows$original[ri],
    n_horizons = rows$n_horizons[ri],
    horizon_first_index = rows$horizon_first_index[ri],
    horizon_last_index = rows$horizon_last_index[ri],
    target_occupancy_cell_id = map$target_occupancy_cell_id[ri],
    target_pair_count = map$target_pair_count[ri],
    control_pool_pair_count = map$control_pool_pair_count[ri],
    risk_pair_count = map$risk_pair_count[ri],
    target_a_unit_count = map$target_a_unit_count[ri],
    control_pool_a_unit_count = map$control_pool_a_unit_count[ri],
    risk_a_unit_count = map$risk_a_unit_count[ri]
  )
}

construct_pipeline_catalog_row_views <- function(
  metadata_obj,
  row_id,
  control = NULL,
  frame_owner = NULL,
  nuisance_view = NULL,
  verify_nuisance_mask = TRUE
) {
  pipeline_metadata_require_stage(c("owners", "views"))
  validate_owner_comparison_catalog(metadata_obj, frame_owner = frame_owner)
  p_own <- get("panel_owner", envir = metadata_obj$panel_ref, inherits = FALSE)
  f_own <- if (!is.null(frame_owner)) {
    validate_ipw_frame_owner(frame_owner, panel_owner = p_own)
    frame_owner
  } else if (!is.null(metadata_obj$frame_ref)) {
    f_tmp <- get("frame_owner", envir = metadata_obj$frame_ref, inherits = FALSE)
    validate_ipw_frame_owner(f_tmp, panel_owner = p_own)
    f_tmp
  } else {
    cdid_pipeline_metadata_abort(
      "A sealed `frame_owner` is required to construct compact `cdid_ipw_nuisance_view` and `cdid_ipw_link_view` objects.",
      c("cdid_pipeline_metadata_owner_error", "cdid_ipw_owner_frame_error")
    )
  }

  ctrl <- if (is.null(control)) metadata_obj$control$owner_control else validate_ipw_owner_control(control)
  pipeline_metadata_row_spec_preflight(metadata_obj,row_id,ctrl)
  spec <- extract_pipeline_catalog_row_spec(metadata_obj, row_id)
  pipeline_metadata_row_views_from_spec(spec, f_own, ctrl, nuisance_view, verify_nuisance_mask)
}

# Views for one extracted row spec. Takes no metadata: the owner constructors
# validate the frame owner, control, any supplied nuisance view and the
# endpoints at their own boundaries on every call. `construction_context` is
# supplied only by the trusted registry loop (see ipw_owner_construction_context).
pipeline_metadata_row_views_from_spec <- function(spec, frame_owner, control,
                                                  nuisance_view = NULL, verify_nuisance_mask = TRUE,
                                                  construction_context = NULL) {
  pipeline_metadata_require_stage("views")
  nu_view <- if (is.null(nuisance_view)) {
    construct_ipw_nuisance_view(
      frame_owner = frame_owner,
      g = spec$g,
      K = spec$K,
      s = spec$s,
      t = spec$t,
      control = control,
      nuisance_id = paste0("nuisance:op", spec$operator_id, ":r", spec$row_id),
      require_support = TRUE,
      construction_context = construction_context
    )
  } else {
    nuisance_view
  }

  lk_view <- construct_ipw_link_view(
    frame_owner = frame_owner,
    s = spec$s,
    t = spec$t,
    g = spec$g,
    K = spec$K,
    nuisance_view = nu_view,
    control = control,
    link_id = paste0("link:r", spec$row_id),
    functional_key = spec$functional_key,
    require_support = TRUE,
    verify_nuisance_mask = verify_nuisance_mask,
    construction_context = construction_context
  )

  list(
    row_spec = spec,
    nuisance_view = nu_view,
    link_view = lk_view,
    endpoints_read = 0L
  )
}
