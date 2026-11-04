# Source-bound lifetime accounting for ordinary static baseline OR/DR only.
# Other routes and standalone backend callers retain their existing contracts.
cp04_pipeline_phase_eligible <- function(cp04, transfer, pretest, imputation, augmented, custom, weighting, strategy) {
  is.list(cp04) && (identical(cp04$or_method, "baseline") || identical(cp04$dr_method, "baseline")) &&
    !isTRUE(transfer) && !isTRUE(pretest) && !isTRUE(imputation) && !isTRUE(augmented) &&
    !is.function(custom) && is.null(cp04$reference_source) && identical(weighting, "optimal") && !is.null(strategy)
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
  b <- registry$cp04; reader <- registry$cp04_reader
  if (!inherits(b, c("cdid_or_baseline", "cdid_dr_baseline")) || !is.null(b$reference_ref))
    cp04_abort("Pipeline memory phases require their static baseline backend.", "integrity")
  if (inherits(b, "cdid_or_baseline")) validate_or_baseline_backend(b) else validate_dr_baseline_backend(b)
  validate_cp04_reader(reader)
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
  if (!ipw_owner_seal_ok(phase, "cp04_pipeline_memory_phases") ||
      !identical(phase$kind, "cp04_pipeline_memory_phases") ||
      !identical(phase$backend_ref, b$ref) || !identical(phase$reader_ref, reader$ref))
    cp04_abort("Public diagnostic phase belongs to a different pipeline backend or reader.", "integrity")
  if (inherits(b, "cdid_or_baseline")) validate_or_baseline_backend(b) else
    if (inherits(b, "cdid_dr_baseline")) validate_dr_baseline_backend(b) else
      cp04_abort("Public diagnostic phases require their fitted baseline backend.", "integrity")
  validate_cp04_reader(reader)
  cp04_registry_retained_entries(res$nuisance_registry, phase$persistent_entries)
  ipw_owner_safe_sum(c(phase$persistent_entries,
    unname(cp04_payload(list(res, public_object))["entries"])), label = "returned result and partial public diagnostic phase")
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
