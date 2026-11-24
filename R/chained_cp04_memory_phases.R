# Source-bound lifetime accounting for ordinary static baseline OR/DR under
# `optimal`, and for every sealed CP04 backend under `within`, whose own
# reduction and working-covariance phases are charged in the pipeline.
# Other routes and standalone backend callers retain their existing contracts.
cp04_pipeline_phase_eligible <- function(cp04, transfer, pretest, imputation, augmented, custom, weighting, strategy) {
  if (!is.list(cp04) || isTRUE(transfer) || isTRUE(pretest) || is.function(custom) || is.null(strategy)) return(FALSE)
  if (identical(weighting, "within")) return(TRUE)
  (identical(cp04$or_method, "baseline") || identical(cp04$dr_method, "baseline")) &&
    !isTRUE(imputation) && !isTRUE(augmented) && is.null(cp04$reference_source) && identical(weighting, "optimal")
}

# The standalone raw IPW backend (`cp04 = NULL`, incomplete panels) under
# `optimal` and `within`, with the same exclusions as the CP04 backends.
ipw_raw_phase_eligible <- function(raw, transfer, pretest, custom, weighting, strategy) {
  isTRUE(raw) && !isTRUE(transfer) && !isTRUE(pretest) && !is.function(custom) && !is.null(strategy) &&
    weighting %in% c("optimal", "within")
}

# Retained-basis covariance and replay terms of the raw IPW structural
# reduction join the deferred terms of the open phase record.
ipw_raw_phase_defer <- function(phase, covariance, replay) {
  if (!identical(phase$kind, "cp04_pipeline_memory_phases") || ipw_owner_seal_ok(phase, "cp04_pipeline_memory_phases"))
    cp04_abort("Only an open pipeline phase record receives deferred terms.", "integrity")
  phase$whole_path_estimate <- ipw_owner_safe_sum(c(phase$whole_path_estimate, covariance, replay))
  phase$deferred_covariance_entries <- ipw_owner_safe_sum(c(phase$deferred_covariance_entries, covariance))
  phase$deferred_replay_entries <- ipw_owner_safe_sum(c(phase$deferred_replay_entries, replay))
  phase
}

# The fitted backend of a memory-phase pipeline: static baseline OR/DR keep
# their exact validators; every other sealed CP04 backend (only under
# `within`) is validated by the retained moment IF validator.
cp04_pipeline_phase_backend <- function(registry) {
  b <- registry$cp04; reader <- registry$cp04_reader
  if (is.null(b) && !is.null(registry$raw_ipw)) {
    validate_ipw_raw_backend(registry$raw_ipw)
    return(registry$raw_ipw)
  }
  if (inherits(b, c("cdid_or_baseline", "cdid_dr_baseline")) && is.null(b$reference_ref)) {
    if (inherits(b, "cdid_or_baseline")) validate_or_baseline_backend(b) else validate_dr_baseline_backend(b)
  } else if (is.list(b) && inherits(b, cp04_supported_backend_classes)) {
    cp04_retained_if_validate_backend(b)
  } else {
    cp04_abort("Pipeline memory phases require their sealed fitted backend.", "integrity")
  }
  validate_cp04_reader(reader)
  b
}

# Entries simultaneously reachable inside `within_mask_spanning_forest()` for
# M catalog links, n_a A units and n_cal dates, beyond the persistent prefix.
# Every term is counted as a double, including logical and integer arrays,
# and an array is counted until it becomes unreachable: date and link masks
# (n_cal + M + 8 vectors of length n_a); the completed forest blocks (sum of
# squared group sizes, at most M^2) together with either the current group's
# dropped-edge map, its cross-product, identity and sum (at most 3.25 M^2) or
# the assembled, reordered, transposed, summed and halved Gram (at most
# 5 M^2); and linear index, key, adjacency and queue buffers (64 M).
cp04_pipeline_within_forest_entries <- function(M, n_a, n_cal) {
  M <- as.double(M); n_a <- as.double(n_a); n_cal <- as.double(n_cal)
  ipw_owner_safe_sum(c((n_cal + M + 8) * n_a, 6 * M * M, 64 * M, 1024),
    label = "within reduction phase entries")
}

# Entries simultaneously reachable inside `compute_pipeline_within_covariance()`
# for L retained links, beyond the persistent prefix: date masks and unit
# vectors (n_cal + 8 vectors of length n_a); the two covariate subsets and the
# per-link treatment and outcome fit scratch (2 n_a p + 16 n_a (p + 1) +
# 4 p^2); the scaled loading matrix (L n_a); and at most ten L x L arrays
# (cohort and date sign masks with the operands and repetitions of each
# `outer()` call, the cross-product and its scaled copy, the signed product,
# its transpose, sum and symmetrized result), plus linear vectors (16 L).
cp04_pipeline_within_covariance_entries <- function(L, n_a, n_cal, p) {
  L <- as.double(L); n_a <- as.double(n_a); n_cal <- as.double(n_cal); p <- as.double(p)
  ipw_owner_safe_sum(c((n_cal + 8) * n_a, 2 * n_a * p, 16 * n_a * (p + 1), 4 * p * p,
    L * n_a, 10 * L * L, 16 * L, 1024), label = "within covariance phase entries")
}

# Entries simultaneously reachable inside `compute_same_pair_within_bundle()`
# beyond the persistent prefix and the fitted backend, apart from its loading
# log, which is charged record by record against the workspace that remains
# (`level_loading_log()`): the panel outcome and observation matrices
# (2 n_a n_cal); link descriptors and views of the M_u catalog and auxiliary
# links (4 M_u n_a); refitted outcome, selection and imputation models, one at
# a time, on at most n_a n_cal stacked rows with p columns (16 n_a n_cal (p + 1)
# + 4 p^2); per-cohort pair-mask tables and row scales (n_g (n_cal^2 + n_a));
# witness blocks of at most all n_a rows (7 M n_a + 4 n_a); the forest
# reduction with its tree-path maps and Gram (6 M^2 + 64 M); and mask keys.
cp04_pipeline_within_bundle_entries <- function(M, M_u, n_a, n_cal, n_g, p) {
  M <- as.double(M); M_u <- as.double(M_u); n_a <- as.double(n_a); n_cal <- as.double(n_cal)
  n_g <- as.double(n_g); p <- as.double(p)
  ipw_owner_safe_sum(c(2 * n_a * n_cal, 4 * M_u * n_a, 16 * n_a * n_cal * (p + 1), 4 * p * p,
    n_g * (n_cal * n_cal + n_a), 7 * M * n_a + 4 * n_a, 6 * M * M, 64 * M + 64 * M_u, 1024),
    label = "within loading phase entries")
}

# Entries reachable inside `within_cov_stream()` for L retained links, block
# rows b and n_a rows, beyond the persistent prefix and the loading log: the
# per-cohort Grams (at most L^2), the result (L^2), each block's cross-product
# (at most L^2), the symmetrization of one cohort (3 L^2), the replay and Gram
# block buffers (4 b L), record row selection (4 n_a) and column maps.
cp04_pipeline_within_stream_entries <- function(L, b, n_a, M) {
  L <- as.double(L); b <- as.double(b); n_a <- as.double(n_a); M <- as.double(M)
  ipw_owner_safe_sum(c(6 * L * L, 4 * b * L, 4 * n_a, M + 16 * L, 1024),
    label = "within loading covariance phase entries")
}

cp04_pipeline_phase_start <- function(whole, covariance, target_output, draws, replay, auxiliary, control) {
  deferred <- ipw_owner_safe_sum(c(covariance, target_output, draws, replay, auxiliary))
  # The small phase record replaces covariance/replay reservation that stays
  # retired even after reporting restores its output/draw/auxiliary terms.
  # Never introduce a new admission cost when that reservation cannot fund it.
  if (ipw_owner_safe_sum(c(covariance, replay)) < 1024) return(NULL)
  if (deferred > whole) cp04_abort("Pipeline phase terms exceed their original prefix.", "integrity")
  persistent <- ipw_owner_safe_sum(c(whole - deferred, 1024),
    min(control$max_retained_entries, control$max_matrix_entries), "persistent baseline pipeline prefix")
  list(kind = "cp04_pipeline_memory_phases", whole_path_estimate = whole,
    persistent_entries = persistent, metadata_entries = 1024,
    deferred_covariance_entries = covariance, deferred_target_output_entries = target_output,
    deferred_draw_entries = draws, deferred_replay_entries = replay, deferred_auxiliary_entries = auxiliary)
}

cp04_pipeline_phase_finish <- function(phase, registry, covariance_gmm, gmm_retained, reporting) {
  if (is.null(phase)) return(NULL)
  b <- cp04_pipeline_phase_backend(registry); reader <- registry$cp04_reader
  # The existing V5 provenance checker still removes EXACTLY one embedded X.
  # It is not a second alias deduction or a change to standalone reservations.
  cp04_registry_retained_entries(registry, phase$persistent_entries)
  phase$backend_ref <- b$ref; phase$reader_ref <- reader$ref
  phase$fit_retained_bound <- b$accounting$retained_bound
  phase$fit_scratch_bound <- b$accounting$fit_scratch_bound
  phase$fit_peak_bound <- ipw_owner_safe_sum(c(phase$fit_retained_bound, phase$fit_scratch_bound))
  phase$covariance_gmm_peak <- covariance_gmm
  phase$gmm_retained_entries <- gmm_retained
  phase$reporting_prefix_entries <- reporting
  phase$deferred_diagnostic_entries <- 5 * (4 * b$accounting$owned_entries + 4096)
  phase$rule <- paste("Persistent owners/input/registry reservations survive every phase; fit/model buffers,",
    "covariance/GMM consumers and reporting/draw/output buffers are composed at their allocation phase;",
    "final public diagnostics execute after pipeline return with retained-result and partial-public payload.")
  if (unname(cp04_payload(phase)["entries"]) + 64 > phase$metadata_entries)
    cp04_abort("Pipeline phase metadata reservation underestimated.", "bound")
  ipw_owner_seal(phase, "cp04_pipeline_memory_phases")
}

# Called AFTER the pipeline returns. Its GMM result retains factors/loadings/IF,
# all counted here. Replay/store/solver scratch is no longer retained; the
# pipeline cleared those aliases at the original release boundaries. res owns
# the panel/frame/reader/backend payload; environments are only sealed refs.
# Counting repeated aliases again is conservative and needs no equality guess.
cp04_public_diagnostic_retained <- function(res, public_object = NULL) {
  phase <- res$resource_accounting$memory_phases
  if (is.null(phase)) return(max(c(0, res$resource_accounting$simultaneous_entries,
    res$resource_accounting$initial_envelope_entries)))
  b <- res$nuisance_registry$cp04; reader <- res$nuisance_registry$cp04_reader
  if (is.null(b)) b <- res$nuisance_registry$raw_ipw
  if (!ipw_owner_seal_ok(phase, "cp04_pipeline_memory_phases") ||
      !identical(phase$kind, "cp04_pipeline_memory_phases") ||
      !identical(phase$backend_ref, b$ref) || !identical(phase$reader_ref, reader$ref))
    cp04_abort("Public diagnostic phase belongs to a different pipeline backend or reader.", "integrity")
  cp04_pipeline_phase_backend(res$nuisance_registry)
  cp04_registry_retained_entries(res$nuisance_registry, phase$persistent_entries)
  # The persistent prefix holds the frame and panel owner storage, and the walk
  # of `res` charges the same owners again; count that storage once.
  owners <- if (is.list(res$frame_owner) && is.list(res$panel_owner))
    cp04_owner_overlap_entries(res$nuisance_registry, "owner") else 0
  ipw_owner_safe_sum(c(phase$persistent_entries - owners,
    unname(cp04_payload(list(cp04_without_owner_record(res), public_object))["entries"])), label = "returned result and partial public diagnostic phase")
}

# The public result keeps scalar accounting, never the private backend/reader
# seal. Otherwise a diagnostic record would accidentally retain all owners.
cp04_pipeline_public_resource <- function(resource) {
  phase <- resource$memory_phases
  if (is.null(phase)) return(resource)
  if (!ipw_owner_seal_ok(phase, "cp04_pipeline_memory_phases"))
    cp04_abort("Malformed private pipeline phase accounting.", "integrity")
  resource$memory_phases <- phase[setdiff(names(phase), c("ref", "backend_ref", "reader_ref"))]
  resource
}
