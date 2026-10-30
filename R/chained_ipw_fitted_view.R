# Private CP02 owner-backed finite weighted-logit nuisance fitting slice.
#
# Scope and boundaries:
#   * Consumes a sealed `cdid_ipw_nuisance_view` (and its sealed `cdid_ipw_frame_owner`
#     and `cdid_ipw_panel_owner`) via bounded outcome-free chunk queries
#     (`query_ipw_view_chunk(..., read_outcomes = FALSE, include_design = TRUE)`).
#   * Streams globally ordered risk-set pair observers (`O = R intersect S`) into a
#     single temporary observer buffer (`risk_a_idx`, `risk_gi`, `risk_w`, `Z_use`),
#     releasing each consumed query chunk before requesting the next chunk.
#   * Calls the shared observer-buffer worker (`ipw_nuisance_observer_fit_worker`) and
#     `ipw_nuisance_hessian_stage` in `R/chained_ipw_link.R`, keeping `stats::glm.fit`
#     as the sole solver without forging a dense context or substituting a streaming
#     optimizer.
#   * Releases all observer-buffer and per-observation model arrays before returning
#     the compact `cdid_ipw_fitted_nuisance_view` (and before any subsequent fit).
#   * Enforces explicit private resource caps before any scan, query, allocation, or
#     callback, charging simultaneous owner-reader chunk and observer buffers,
#     fitting/QR workspace, derivatives, character byte metadata, copies, and compact
#     retained state, while disclosing opaque library/process memory rather than
#     promising an OS RSS cap.
#   * Discloses both the observer-buffer point summaries that define the candidate fit
#     and the chunked view-support point summaries. Neither `glm.fit` convergence nor
#     object sealing certifies a mathematical root, enclosure, or statistical inference.

cdid_ipw_fitted_view_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_ipw_fitted_view_refusal",
      "cdid_ipw_fitted_view_error",
      "cdid_ipw_owner_error",
      "cdid_ipw_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

ipw_fitted_view_require_deps <- function() {
  ns_env <- environment(ipw_fitted_view_require_deps)
  if (cdid_dependency_check_done("ipw_fitted_view", ns_env)) return(invisible(TRUE))
  req_fns <- c(
    "cdid_abort", "gmm_sexp_aliases", "gmm_scalar_numeric", "gmm_scalar_logical",
    "gmm_scalar_character", "logit_cert_count_output", "ipw_owner_seal",
    "ipw_owner_seal_ok", "ipw_owner_refuse_decorative_attrs", "ipw_owner_int_bound",
    "ipw_owner_safe_product", "ipw_owner_safe_sum", "ipw_owner_check_matrix_dims",
    "ipw_owner_control", "validate_ipw_owner_control", "validate_ipw_panel_owner",
    "validate_ipw_frame_owner", "validate_ipw_nuisance_view", "validate_ipw_link_view",
    "construct_ipw_nuisance_view", "query_ipw_view_chunk",
    "verify_ipw_view_mask_compatibility", "validate_ipw_glm_control",
    "ipw_nuisance_observer_fit_worker", "ipw_nuisance_hessian_stage"
  )
  for (fn in req_fns) {
    if (!exists(fn, envir = ns_env, mode = "function", inherits = TRUE)) {
      cdid_ipw_fitted_view_abort(
        paste0("Required private helper `", fn, "` is not loaded in the namespace."),
        "cdid_ipw_fitted_view_dependency_error"
      )
    }
  }
  cdid_dependency_check_record("ipw_fitted_view", ns_env)
}

# ---------------------------------------------------------------------------
# Explicit private resource controls (NO internal defaults)
# ---------------------------------------------------------------------------

ipw_fitted_view_control <- function(
  max_chunk_rows,
  max_observer_units,
  max_retained_entries,
  max_matrix_entries,
  max_work_ops,
  max_output_entries
) {
  ipw_fitted_view_require_deps()
  if (missing(max_chunk_rows) || missing(max_observer_units) ||
      missing(max_retained_entries) || missing(max_matrix_entries) ||
      missing(max_work_ops) || missing(max_output_entries)) {
    cdid_ipw_fitted_view_abort(
      "All fields of `ipw_fitted_view_control` (`max_chunk_rows`, `max_observer_units`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) must be supplied explicitly by the caller (no internal defaults).",
      c("cdid_ipw_fitted_view_control_error", "cdid_ipw_owner_control_error", "cdid_ipw_owner_bound_error")
    )
  }
  for (item in list(
    list(v = max_chunk_rows, nm = "max_chunk_rows"),
    list(v = max_observer_units, nm = "max_observer_units"),
    list(v = max_retained_entries, nm = "max_retained_entries"),
    list(v = max_matrix_entries, nm = "max_matrix_entries"),
    list(v = max_work_ops, nm = "max_work_ops"),
    list(v = max_output_entries, nm = "max_output_entries")
  )) {
    ipw_owner_refuse_decorative_attrs(item$v, item$nm)
  }

  max_ch <- ipw_owner_int_bound(max_chunk_rows, "max_chunk_rows", min_val = 1)
  max_obs <- ipw_owner_int_bound(max_observer_units, "max_observer_units", min_val = 1)
  max_ret <- ipw_owner_int_bound(max_retained_entries, "max_retained_entries", min_val = 1)
  max_me <- ipw_owner_int_bound(max_matrix_entries, "max_matrix_entries", min_val = 1)
  max_wo <- ipw_owner_int_bound(max_work_ops, "max_work_ops", min_val = 1)
  max_out <- ipw_owner_int_bound(max_output_entries, "max_output_entries", min_val = 1)

  structure(
    list(
      max_chunk_rows = max_ch,
      max_observer_units = max_obs,
      max_retained_entries = max_ret,
      max_matrix_entries = max_me,
      max_work_ops = max_wo,
      max_output_entries = max_out
    ),
    class = "cdid_ipw_fitted_view_control"
  )
}

validate_ipw_fitted_view_control <- function(control) {
  ipw_fitted_view_require_deps()
  ipw_owner_refuse_decorative_attrs(control, "control")
  req <- c(
    "max_chunk_rows", "max_observer_units", "max_retained_entries",
    "max_matrix_entries", "max_work_ops", "max_output_entries"
  )
  if (missing(control) || is.null(control) || !is.list(control) ||
      length(control) != length(req) || is.null(names(control)) ||
      anyDuplicated(names(control)) || !setequal(names(control), req)) {
    cdid_ipw_fitted_view_abort(
      "Explicit `control` list (`max_chunk_rows`, `max_observer_units`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) is required; no internal defaults exist.",
      c("cdid_ipw_fitted_view_control_error", "cdid_ipw_owner_control_error", "cdid_ipw_owner_bound_error")
    )
  }
  ipw_fitted_view_control(
    max_chunk_rows = control$max_chunk_rows,
    max_observer_units = control$max_observer_units,
    max_retained_entries = control$max_retained_entries,
    max_matrix_entries = control$max_matrix_entries,
    max_work_ops = control$max_work_ops,
    max_output_entries = control$max_output_entries
  )
}

ipw_fitted_view_to_owner_control <- function(ctrl) {
  ipw_owner_control(
    max_chunk_rows = ctrl$max_chunk_rows,
    max_retained_entries = ctrl$max_retained_entries,
    max_matrix_entries = ctrl$max_matrix_entries,
    max_work_ops = ctrl$max_work_ops,
    max_output_entries = ctrl$max_output_entries
  )
}

ipw_fitted_view_char_bytes <- function(chr_list, control) {
  total_bytes <- 0; max_bytes <- 0; work <- 128
  for (x in chr_list) {
    if (!is.null(x) && is.character(x)) for(j in seq_along(x)) {
      # Scalar byte lengths avoid allocating a character-length numeric array.
      nb <- as.double(nchar(x[j],type="bytes"))
      work <- ipw_owner_safe_sum(c(work,nb,8),control$max_work_ops,"character metadata scan work")
      total_bytes <- ipw_owner_safe_sum(c(total_bytes,nb),label="character metadata byte count")
      max_bytes <- max(max_bytes,nb)
      ipw_owner_safe_sum(c(ceiling(total_bytes/2),32),control$max_matrix_entries,"character metadata scratch slots")
    }
  }
  slots <- ipw_owner_safe_sum(c(ceiling(total_bytes/2),32),label="character metadata charged slots")
  list(source_bytes=total_bytes,max_string_bytes=max_bytes,charged_slots=slots,
       charged_bytes=ipw_owner_safe_product(8,slots),work_ops_bound=work)
}

# Descriptive summaries of existing worker vectors only. At most one bounded
# chunk is normalized at a time; no original-scale squares or retained arrays.
# The caller reserves metadata, scratch and work in the fitted-view preflight.
ipw_fitted_diagnostic_summary <- function(v, chunk_rows, kish) {
  unavailable <- function(reason) {
    if (kish) list(status = "unavailable", reason = reason, value = NA_real_) else
      list(status = "unavailable", reason = reason, min = NA_real_, max = NA_real_)
  }
  n <- length(v)
  if (!n) return(unavailable("empty_vector"))
  lo <- Inf; hi <- -Inf; start <- 1.0
  while (start <= n) {
    end <- min(n, start + chunk_rows - 1.0)
    x <- v[seq.int(start, end)]
    if (any(!is.finite(x))) return(unavailable("nonfinite_vector"))
    lo <- min(lo, min(x)); hi <- max(hi, max(x))
    x <- NULL
    start <- end + 1.0
  }
  if (!kish) return(list(status = "available", reason = NA_character_, min = lo, max = hi))
  if (lo <= 0) return(unavailable("nonpositive_weights"))
  total <- 0.0; squares <- 0.0; start <- 1.0
  while (start <= n) {
    end <- min(n, start + chunk_rows - 1.0)
    u <- v[seq.int(start, end)] / hi
    total <- total + sum(u)
    squares <- squares + sum(u * u)
    u <- NULL
    start <- end + 1.0
  }
  value <- total^2 / squares
  if (!is.finite(value)) return(unavailable("nonfinite_diagnostic_arithmetic"))
  list(status = "available", reason = NA_character_, value = value)
}

# ---------------------------------------------------------------------------
# Preflight accounting (executed strictly BEFORE any query, buffer allocation,
# or callback invocation)
# ---------------------------------------------------------------------------

preflight_ipw_fitted_nuisance_view <- function(
  nuisance_view,
  glm_control,
  control,
  n_planned_fits = 1L
) {
  ipw_fitted_view_require_deps()
  if (missing(nuisance_view) || missing(glm_control) || missing(control)) {
    cdid_ipw_fitted_view_abort(
      "Required arguments (`nuisance_view`, `glm_control`, `control`) must be supplied explicitly.",
      "cdid_ipw_fitted_view_input_error"
    )
  }
  ipw_owner_refuse_decorative_attrs(nuisance_view, "nuisance_view")
  ipw_owner_refuse_decorative_attrs(glm_control, "glm_control")
  ipw_owner_refuse_decorative_attrs(control, "control")

  validate_ipw_nuisance_view(nuisance_view)
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  ctrl <- validate_ipw_fitted_view_control(control)
  n_fits <- ipw_owner_int_bound(n_planned_fits, "n_planned_fits", min_val = 1)

  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = nuisance_view$panel_ref, inherits = FALSE)
  validate_ipw_frame_owner(frame_owner, panel_owner = panel_owner)

  na <- ipw_owner_int_bound(frame_owner$n_a_units, "n_a_units", min_val = 1)
  n_obs <- ipw_owner_int_bound(nuisance_view$support$n_observed_risk_units, "n_observed_risk_units", min_val = 0)
  p <- ipw_owner_int_bound(frame_owner$n_covariates, "n_covariates", min_val = 1)
  m <- ipw_owner_int_bound(nuisance_view$cohorts$m_cohorts, "m_cohorts", min_val = 2)

  max_ch <- ctrl$max_chunk_rows
  max_obs <- ctrl$max_observer_units
  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  if (n_obs > max_obs) {
    cdid_ipw_fitted_view_abort(
      paste0("Observed risk-set pair unit count (", n_obs, ") exceeds `control$max_observer_units` (", max_obs, ")."),
      c("cdid_ipw_fitted_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error"),
      diagnostics = list(n_observed_risk_units = n_obs, max_observer_units = max_obs)
    )
  }

  B_eff <- min(na, max_ch)
  n_chunks <- ceiling(na / B_eff)
  # Refusals may intentionally return failed predictions; reserve that output too.
  failure_output_entries <- ipw_owner_safe_sum(c(ipw_owner_safe_product(4,n_obs),
    ipw_owner_safe_product(p,m),p,ipw_owner_safe_product(4,m),256),
    min(max_ret,max_out),"possible failed-fit diagnostic output")

  # Constant/header-only allowance before any character content traversal.
  header_slots <- ipw_owner_safe_sum(c(ipw_owner_safe_product(8,p+m+1),32),
    min(max_ret,max_out,max_me),"fitted metadata header before character scans")
  ipw_owner_safe_sum(c(ipw_owner_safe_product(16,p+m+1),128),max_wo,"fitted header work before character scans")

  # 1. Character metadata and copy charges (covariate names, cohort labels, view ID)
  char_meta <- ipw_fitted_view_char_bytes(list(
    colnames(frame_owner$Z_a),
    nuisance_view$cohorts$cohort_labels,
    nuisance_view$nuisance_id
  ),ctrl)

  # 2. Compact retained state per fit:
  # Matrices: H (p x p), H_chol (p x p), qr_R (p x p), A_beta_r (p x m)
  p2 <- ipw_owner_check_matrix_dims(p, p, max_matrix_entries = max_ret, label = "fitted nuisance p x p factor matrix")
  pm <- ipw_owner_check_matrix_dims(p, m, max_matrix_entries = max_ret, label = "fitted nuisance p x m score-rate Jacobian")
  retained_matrix_entries_one <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(3, p2, max_ret, "H, H_chol, and qr_R matrices"), pm,nuisance_view$preflight$retained_matrix_entries),
    max_bound = max_ret,
    label = "retained fitted nuisance matrix entries"
  )
  # Vectors & scalars: beta (2*p for coefficients + alias), qr_pivot (p),
  # cohort counts/masses/rates/chunked disclosures (10*m), scalar means/diagnostics/metadata + char slots
  retained_vector_entries_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(3, p, max_ret, "beta coefficients and qr_pivot"),
      ipw_owner_safe_product(12, m, max_ret, "cohort counts, masses, observer rates, and chunked view rates"),
      char_meta$charged_slots,
      nuisance_view$preflight$retained_vector_entries,
      2048, # D27 scalar/status/character diagnostic output (no observer arrays)
      128
    ),
    max_bound = max_ret,
    label = "retained fitted nuisance vector and metadata entries"
  )
  retained_entries_one <- ipw_owner_safe_sum(
    c(retained_matrix_entries_one, retained_vector_entries_one),
    max_bound = max_ret,
    label = "single fitted nuisance retained entries"
  )
  total_retained_entries <- ipw_owner_safe_product(
    n_fits, retained_entries_one, max_entries = max_ret,
    label = "total retained fitted nuisance entries across planned fits"
  )
  if (total_retained_entries > max_out) {
    cdid_ipw_fitted_view_abort(
      paste0("Total retained fitted nuisance entries (", total_retained_entries, ") exceed `control$max_output_entries` (", max_out, ")."),
      c("cdid_ipw_fitted_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error"),
      diagnostics = list(requested = total_retained_entries, max_output_entries = max_out)
    )
  }

  # 3. Single temporary observer buffer: risk_a_idx (n_obs), risk_gi (n_obs), risk_w (n_obs), Z_use (n_obs x p)
  obs_Z_entries <- ipw_owner_check_matrix_dims(n_obs, p, max_matrix_entries = max_me, label = "observer buffer Z_use matrix")
  observer_buffer_entries <- ipw_owner_safe_sum(
    c(obs_Z_entries, ipw_owner_safe_product(3, max(1, n_obs), max_me, "observer buffer index/cohort/weight vectors"), 32),
    max_bound = max_me,
    label = "temporary observer buffer entries"
  )

  # 4. Simultaneous bounded owner-reader chunk buffer (at most 1 chunk active at a time during streaming):
  owner_ctrl <- ipw_fitted_view_to_owner_control(ctrl)
  chunk_pre <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p,
    m_cohorts = m,
    control = owner_ctrl,
    read_outcomes = FALSE,
    include_design = TRUE,
    enclose_differences = FALSE
  )
  reader_chunk_concurrent_entries <- chunk_pre$concurrent_entries

  # 5. Fitting, QR, and derivative temporary workspace (active after reader chunks are released):
  # stats::glm.fit allocates working design copy (n_obs x p), QR matrix (n_obs x p),
  # working vectors (D_use, r_use, fit_weights, mu, eta, weights, residuals, effects, w, z, etc.),
  # and ipw_nuisance_hessian_stage allocates weighted design (n_obs x p), w_hess, score_factor,
  # plus p x p and p x m derivative workspaces and copy buffers.
  fitting_qr_matrix_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(8, obs_Z_entries, max_me, "glm.fit working x, QR matrix, and Hessian/Jacobian weighted Z copies"),
      ipw_owner_safe_product(8, p2, max_me, "QR R, Hessian H, Cholesky H_chol, and chol2inv workspaces"),
      ipw_owner_safe_product(4, pm, max_me, "score-rate Jacobian A_beta_r workspace")
    ),
    max_bound = max_me,
    label = "fitting, QR, and derivative matrix workspace entries"
  )
  fitting_qr_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(64, max(1, n_obs), max_me, "worker, glm.fit IRLS, and Hessian per-observer working vectors"),
      ipw_owner_safe_product(16, p, max_me, "QR qraux/pivot/coefficient working vectors"),
      ipw_owner_safe_product(16, m, max_me, "cohort rate and mass working vectors"),
      char_meta$charged_slots,
      8 * min(B_eff, max(1, n_obs)), # D27 normalized chunk and expression copies
      2048, # D27 scalar/status/character scratch
      128
    ),
    max_bound = max_me,
    label = "fitting, QR, and derivative vector workspace entries"
  )
  fitting_qr_entries <- ipw_owner_safe_sum(
    c(fitting_qr_matrix_entries, fitting_qr_vector_entries),
    max_bound = max_me,
    label = "total fitting, QR, and derivative workspace entries"
  )

  # 6. Simultaneous peak temporary entries (streaming phase vs solver/derivative phase, plus retained outputs):
  streaming_phase_entries <- ipw_owner_safe_sum(
    c(total_retained_entries, observer_buffer_entries, reader_chunk_concurrent_entries, char_meta$charged_slots),
    max_bound = max_me,
    label = "simultaneous streaming-phase entries (observer buffer + 1 reader chunk + retained)"
  )
  solver_phase_entries <- ipw_owner_safe_sum(
    c(total_retained_entries, observer_buffer_entries, fitting_qr_entries),
    max_bound = max_me,
    label = "simultaneous solver-phase entries (observer buffer + fitting/QR/derivative workspace + retained)"
  )
  peak_concurrent_entries <- max(streaming_phase_entries, solver_phase_entries)
  peak_concurrent_matrix_entries <- max(
    retained_matrix_entries_one * n_fits + obs_Z_entries + chunk_pre$output_matrix_entries,
    retained_matrix_entries_one * n_fits + obs_Z_entries + fitting_qr_matrix_entries
  )
  if (peak_concurrent_matrix_entries > max_me) {
    cdid_ipw_fitted_view_abort(
      paste0("Peak concurrent matrix entries (", peak_concurrent_matrix_entries, ") exceed `control$max_matrix_entries` (", max_me, ")."),
      c("cdid_ipw_fitted_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error")
    )
  }

  # 7. Work operations charge (chunk streaming + IRLS iterations + Hessian/Jacobian):
  maxit_d <- ipw_owner_int_bound(glm_ctrl$maxit, "glm_control$maxit", min_val = 1)
  stream_work_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, chunk_pre$work_ops_bound, max_wo, "bounded view chunk query work"),
      ipw_owner_safe_product(4, max(1, n_obs) * p, max_wo, "observer buffer copy and order verification work")
    ),
    max_bound = max_wo,
    label = "single-fit stream work operations"
  )
  irls_step_ops <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(4, max(1, n_obs) * p * p, max_wo, "QR decomposition ops per IRLS iteration"),
      ipw_owner_safe_product(16, max(1, n_obs) * p, max_wo, "linear predictor and score ops per IRLS iteration"),
      ipw_owner_safe_product(16, max(1, n_obs), max_wo, "link/weight ops per IRLS iteration"),
      64
    ),
    max_bound = max_wo,
    label = "single IRLS iteration work operations"
  )
  solver_work_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(maxit_d, irls_step_ops, max_wo, "total IRLS solver work bound"),
      ipw_owner_safe_product(4, max(1, n_obs) * p * p + max(1, n_obs) * p * m + p * p * p, max_wo, "Hessian, Cholesky, and score-rate Jacobian work"),
      ipw_owner_safe_product(128, max(1, n_obs), max_wo, "D27 three Kish summaries and fitted probability range"),
      4096,
      256
    ),
    max_bound = max_wo,
    label = "single-fit solver and derivative work operations"
  )
  work_ops_one <- ipw_owner_safe_sum(
    c(stream_work_one, solver_work_one,char_meta$work_ops_bound),
    max_bound = max_wo,
    label = "single-fit total work operations"
  )
  total_work_ops_bound <- ipw_owner_safe_product(
    n_fits, work_ops_one, max_entries = max_wo,
    label = "total planned fitted nuisance work operations"
  )

  list(
    n_a_units = na,
    n_observed_risk_units = n_obs,
    p_covariates = p,
    m_cohorts = m,
    n_planned_fits = n_fits,
    nuisance_diagnostics = list(retained_entries = 2048,
      scratch_entries = 8 * min(B_eff, max(1, n_obs)) + 2048,
      work_ops_bound = 128 * max(1, n_obs) + 4096,
      scope = "existing counted workspace/work; not process RSS or CPU time"),
    failure_output_entries = failure_output_entries,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_fit = logit_cert_count_output(n_chunks),
    max_simultaneous_query_chunks = 1L,
    max_simultaneous_observer_buffers = 1L,
    character_metadata = char_meta,
    observer_buffer_entries = observer_buffer_entries,
    reader_chunk_concurrent_entries = reader_chunk_concurrent_entries,
    fitting_qr_workspace_entries = fitting_qr_entries,
    streaming_phase_entries = streaming_phase_entries,
    solver_phase_entries = solver_phase_entries,
    peak_concurrent_entries = peak_concurrent_entries,
    peak_concurrent_matrix_entries = peak_concurrent_matrix_entries,
    peak_concurrent_bytes_8double = 8 * peak_concurrent_entries,
    referenced_nuisance_entries = nuisance_view$preflight$total_retained_entries,
    retained_matrix_entries_per_fit = retained_matrix_entries_one,
    retained_vector_entries_per_fit = retained_vector_entries_one,
    retained_entries_per_fit = retained_entries_one,
    total_retained_entries = total_retained_entries,
    total_retained_bytes_8double = 8 * total_retained_entries,
    work_ops_per_fit_bound = work_ops_one,
    total_work_ops_bound = total_work_ops_bound,
    endpoints_read = 0L,
    peak_rss_promised = FALSE,
    memory_disclosure = paste0(
      "Deterministic preflight charges simultaneous bounded owner-reader chunk buffers, ",
      "the single temporary observer buffer (n_obs x p design plus index/cohort/weight vectors), ",
      "stats::glm.fit IRLS/QR workspaces, Hessian/Cholesky/score-rate derivative workspaces, ",
      "character metadata byte slots, copy buffers, and compact retained outputs. ",
      "Opaque C-runtime, BLAS/LAPACK, and OS process RSS allocations are disclosed as outside ",
      "R-level object accounting and no process RSS cap is promised."
    )
  )
}

# Verify that a compact fitted nuisance view contains zero observer-length or unit-length arrays.
ipw_fitted_view_verify_compact_storage <- function(obj) {
  forbidden_fields <- c(
    "risk_a_idx", "risk_gi", "risk_w", "Z_use", "D_use", "fit_weights",
    "eta_use", "odds_use", "prob_use", "is_t", "is_c", "a_obs", "b_obs",
    "w_hess", "score_factor", "y_left_obs", "y_right_obs", "Delta_iv_obs",
    "observed_a_index", "fitted.values", "linear.predictors", "residuals",
    "effects", "weights", "y", "qr"
  )
  check_node <- function(x, path) {
    if (is.environment(x)) return(invisible(TRUE))
    if (is.list(x)) {
      nms <- names(x)
      if (!is.null(nms)) {
        # First offending field in `names(x)` order: the primitive `match` gives
        # exactly the element `intersect(nms, forbidden_fields)[1L]` named,
        # without the set helper's coercion, uniqueness pass and allocation.
        m_bad <- match(nms, forbidden_fields, 0L)
        if (any(m_bad > 0L)) {
          cdid_ipw_fitted_view_abort(
            paste0("Forbidden non-compact observer/model field `", nms[which(m_bad > 0L)[1L]], "` retained at `", path, "`."),
            c("cdid_ipw_fitted_view_storage_error", "cdid_ipw_owner_integrity_error")
          )
        }
      }
      for (i in seq_along(x)) {
        nm_i <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
        check_node(x[[i]], paste0(path, "$", nm_i))
      }
    }
    invisible(TRUE)
  }
  check_node(obj, "fitted_nuisance_view")
}

# ---------------------------------------------------------------------------
# Owner-backed finite weighted-logit nuisance fit with one temporary observer buffer
# ---------------------------------------------------------------------------

fit_ipw_nuisance_view <- function(
  nuisance_view,
  glm_control,
  control,
  callback = NULL,
  trace_env = NULL,
  failure_trace_env = NULL
) {
  ipw_fitted_view_require_deps()
  if (missing(nuisance_view) || missing(glm_control) || missing(control)) {
    cdid_ipw_fitted_view_abort(
      "All required arguments (`nuisance_view`, `glm_control`, `control`) must be supplied explicitly.",
      "cdid_ipw_fitted_view_input_error"
    )
  }
  if (!is.null(callback) && !is.function(callback)) {
    cdid_ipw_fitted_view_abort("`callback` must be NULL or a function.", "cdid_ipw_fitted_view_input_error")
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_fitted_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_fitted_view_input_error")
  }
  if (!is.null(failure_trace_env) && !is.environment(failure_trace_env)) {
    cdid_ipw_fitted_view_abort("`failure_trace_env` must be NULL or an environment.", "cdid_ipw_fitted_view_input_error")
  }

  # 1. Preflight strictly BEFORE any scan, query, buffer allocation, or callback!
  preflight <- preflight_ipw_fitted_nuisance_view(
    nuisance_view = nuisance_view,
    glm_control = glm_control,
    control = control,
    n_planned_fits = 1L
  )

  glm_ctrl <- validate_ipw_glm_control(glm_control)
  ctrl <- validate_ipw_fitted_view_control(control)
  owner_ctrl <- ipw_fitted_view_to_owner_control(ctrl)
  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = nuisance_view$panel_ref, inherits = FALSE)

  na <- as.integer(frame_owner$n_a_units)
  n_obs <- as.integer(nuisance_view$support$n_observed_risk_units)
  p <- as.integer(frame_owner$n_covariates)
  m <- as.integer(nuisance_view$cohorts$m_cohorts)
  step <- as.integer(preflight$effective_chunk_rows)

  # Dedicated temporary environments for the single observer buffer and failure diagnostics.
  # Scrubbed unconditionally via on.exit so no observer/model arrays survive refusal or return.
  obs_buf_env <- new.env(parent = emptyenv())
  diag_env <- new.env(parent = emptyenv())
  obs_buf_env$active <- FALSE
  obs_buf_env$active_query_chunks <- 0L
  obs_buf_env$max_simultaneous_query_chunks <- 0L
  obs_buf_env$released_query_chunks <- 0L
  obs_buf_env$observer_released <- FALSE

  cleanup_temporary_buffers <- function() {
    for (nm in c("risk_a_idx", "risk_gi", "risk_w", "Z_use", "st", "hs")) {
      if (exists(nm, envir = obs_buf_env, inherits = FALSE)) {
        assign(nm, NULL, envir = obs_buf_env)
        rm(list = nm, envir = obs_buf_env)
      }
    }
    obs_buf_env$active <- FALSE
    obs_buf_env$observer_released <- TRUE
    if (exists("error_diagnostics", envir = diag_env, inherits = FALSE)) {
      assign("error_diagnostics", NULL, envir = diag_env)
      rm(list = "error_diagnostics", envir = diag_env)
    }
    if (!is.null(trace_env)) {
      trace_env$active_observer_buffers <- 0L
      trace_env$observer_released <- TRUE
    }
  }
  on.exit(cleanup_temporary_buffers(), add = TRUE)

  if (!is.null(trace_env)) {
    trace_env$preflight_completed <- TRUE
    trace_env$queries_started <- 0L
    trace_env$chunks_released <- 0L
    trace_env$max_simultaneous_query_chunks <- 0L
    trace_env$active_observer_buffers <- 0L
    trace_env$max_simultaneous_observer_buffers <- 0L
    trace_env$observer_released <- FALSE
    trace_env$endpoints_read <- 0L
  }

  # 2. Allocate the single temporary observer buffer
  obs_buf_env$active <- TRUE
  obs_buf_env$risk_a_idx <- integer(n_obs)
  obs_buf_env$risk_gi <- integer(n_obs)
  obs_buf_env$risk_w <- numeric(n_obs)
  obs_buf_env$Z_use <- matrix(
    0.0,
    nrow = n_obs,
    ncol = p,
    dimnames = list(NULL, colnames(frame_owner$Z_a))
  )
  if (!is.null(trace_env)) {
    trace_env$active_observer_buffers <- 1L
    trace_env$max_simultaneous_observer_buffers <- max(
      as.integer(trace_env$max_simultaneous_observer_buffers),
      1L
    )
  }

  # 3. Stream bounded outcome-free view chunks into the single observer buffer,
  #    releasing each consumed query chunk before requesting the next one.
  write_pos <- 0L
  last_a_idx <- 0L
  chunks_queried <- 0L
  stream_ops <- 0
  y_trace <- new.env(parent = emptyenv())

  cur_s <- 1L
  while (cur_s <= na) {
    ipw_owner_check_interrupt()
    cur_e <- as.integer(min(as.double(na), as.double(cur_s) + step - 1))

    ch <- query_ipw_view_chunk(
      view = nuisance_view,
      a_start = cur_s,
      a_end = cur_e,
      control = owner_ctrl,
      read_outcomes = FALSE,
      include_design = TRUE,
      enclose_differences = FALSE,
      trace_env = y_trace
    )
    chunks_queried <- chunks_queried + 1L
    obs_buf_env$active_query_chunks <- obs_buf_env$active_query_chunks + 1L
    obs_buf_env$max_simultaneous_query_chunks <- max(
      obs_buf_env$max_simultaneous_query_chunks,
      obs_buf_env$active_query_chunks
    )
    if (!is.null(trace_env)) {
      trace_env$queries_started <- chunks_queried
      trace_env$max_simultaneous_query_chunks <- obs_buf_env$max_simultaneous_query_chunks
      trace_env$endpoints_read <- as.integer(trace_env$endpoints_read) +
        as.integer(ch$read_diagnostics$n_endpoints_read)
    }

    if (isTRUE(ch$read_diagnostics$read_outcomes) ||
        !identical(as.integer(ch$read_diagnostics$n_endpoints_read), 0L) ||
        !identical(as.integer(y_trace$calls), 0L) ||
        !is.null(ch$observed_risk$y_left_obs) ||
        !is.null(ch$observed_risk$y_right_obs) ||
        !is.null(ch$observed_risk$Delta_iv_obs)) {
      cdid_ipw_fitted_view_abort(
        "Outcome-free nuisance view query read or returned outcome endpoints.",
        c("cdid_ipw_fitted_view_outcome_leak_error", "cdid_ipw_owner_integrity_error")
      )
    }

    k_ch <- as.integer(ch$counts$observed_risk_row_count)
    if (k_ch > 0L) {
      if (write_pos + k_ch > n_obs) {
        cdid_ipw_fitted_view_abort(
          "Streamed risk-set observer count exceeds view support `n_observed_risk_units`.",
          c("cdid_ipw_fitted_view_count_error", "cdid_ipw_owner_integrity_error")
        )
      }
      a_ch <- as.integer(ch$observed_risk$a_index)
      gi_ch <- as.integer(ch$observed_risk$cohort_index)
      loc_ch <- as.integer(ch$observed_risk$local_pos)
      w_ch <- as.double(ch$weight[loc_ch])
      Z_ch <- ch$design$Z_obs

      if (length(a_ch) != k_ch || length(gi_ch) != k_ch || length(w_ch) != k_ch ||
          !is.matrix(Z_ch) || nrow(Z_ch) != k_ch || ncol(Z_ch) != p) {
        cdid_ipw_fitted_view_abort(
          "Malformed observed-risk slice in view chunk.",
          c("cdid_ipw_fitted_view_chunk_error", "cdid_ipw_owner_integrity_error")
        )
      }
      if (a_ch[1L] <= last_a_idx || (k_ch > 1L && any(diff(a_ch) <= 0L)) ||
          a_ch[1L] < cur_s || a_ch[k_ch] > cur_e) {
        cdid_ipw_fitted_view_abort(
          "Streamed risk-set observer `a_index` is not strictly increasing in global A-unit order.",
          c("cdid_ipw_fitted_view_order_error", "cdid_ipw_owner_order_mismatch_error")
        )
      }
      if (anyNA(gi_ch) || any(gi_ch < 1L | gi_ch > m)) {
        cdid_ipw_fitted_view_abort(
          "Streamed risk-set observer cohort index is out of `[1, m_cohorts]`.",
          c("cdid_ipw_fitted_view_cohort_error", "cdid_ipw_owner_cohort_error")
        )
      }

      dst <- seq.int(write_pos + 1L, write_pos + k_ch)
      obs_buf_env$risk_a_idx[dst] <- a_ch
      obs_buf_env$risk_gi[dst] <- gi_ch
      obs_buf_env$risk_w[dst] <- w_ch
      obs_buf_env$Z_use[dst, ] <- Z_ch
      write_pos <- write_pos + k_ch
      last_a_idx <- a_ch[k_ch]
      rm(a_ch,gi_ch,loc_ch,w_ch,Z_ch,dst)
    }
    stream_ops <- stream_ops + ch$work$actual_work_ops + k_ch * p

    # Release consumed query chunk BEFORE next query!
    ch <- NULL
    obs_buf_env$active_query_chunks <- obs_buf_env$active_query_chunks - 1L
    obs_buf_env$released_query_chunks <- obs_buf_env$released_query_chunks + 1L
    if (!is.null(trace_env)) {
      trace_env$chunks_released <- obs_buf_env$released_query_chunks
    }
    if (!is.null(callback)) {
      callback(list(
        phase = "after_chunk_release",
        chunk_index = chunks_queried,
        active_query_chunks = obs_buf_env$active_query_chunks,
        released_query_chunks = obs_buf_env$released_query_chunks,
        observer_rows_buffered = write_pos,
        endpoints_read = 0L
      ))
    }

    cur_s <- as.double(cur_e) + 1
  }

  if (write_pos != n_obs || obs_buf_env$active_query_chunks != 0L) {
    cdid_ipw_fitted_view_abort(
      "Streamed risk-set observer count disagrees with `nuisance_view$support$n_observed_risk_units`.",
      c("cdid_ipw_fitted_view_count_error", "cdid_ipw_owner_integrity_error")
    )
  }

  # 4. Build support descriptor for the shared observer-buffer worker.
  sup <- list(
    n = na,
    g = nuisance_view$cohorts$g,
    K = nuisance_view$cohorts$K,
    cohorts = nuisance_view$cohorts$cohorts,
    cohort_labels = nuisance_view$cohorts$cohort_labels,
    m_cohorts = m,
    group_count = as.integer(nuisance_view$support$by_cohort_a_count),
    group_a_count = as.integer(nuisance_view$support$by_cohort_a_count),
    group_mass = as.double(nuisance_view$support$by_cohort_a_mass),
    group_a_mass = as.double(nuisance_view$support$by_cohort_a_mass),
    q = as.double(nuisance_view$support$q),
    a_weight_total = as.double(frame_owner$sum_w_a)
  )

  # 5. Execute the shared observer-buffer weighted-logit worker and Hessian stage.
  #    On any refusal, attach failed predictions/output to the error condition's
  #    `error_diagnostics` (and optional `failure_trace_env`) and scrub the observer
  #    buffer before propagating the refusal.
  tryCatch(
    {
      obs_buf_env$st <- ipw_nuisance_observer_fit_worker(
        sup = sup,
        risk_a_idx = obs_buf_env$risk_a_idx,
        risk_gi = obs_buf_env$risk_gi,
        risk_w = obs_buf_env$risk_w,
        Z_use = obs_buf_env$Z_use,
        ctrl = glm_ctrl,
        diag_env = diag_env
      )
      obs_buf_env$hs <- tryCatch(
        ipw_nuisance_hessian_stage(obs_buf_env$st, sup$cohort_labels),
        error = function(e) {
          diag_env$error_diagnostics <- list(
            stage = "hessian_stage",
            converged = isTRUE(obs_buf_env$st$fit_info$converged),
            boundary = isTRUE(obs_buf_env$st$fit_info$boundary),
            iter = obs_buf_env$st$fit_info$iter,
            rank = obs_buf_env$st$fit_info$rank,
            failed_coefficients = obs_buf_env$st$beta,
            failed_fitted_values = obs_buf_env$st$prob_use,
            failed_odds_values = obs_buf_env$st$odds_use
          )
          stop(e)
        }
      )
      if (any(!is.finite(obs_buf_env$hs$A_beta_r))) {
        diag_env$error_diagnostics <- list(
          stage = "score_rate_jacobian_check",
          failed_coefficients = obs_buf_env$st$beta,
          failed_fitted_values = obs_buf_env$st$prob_use,
          failed_A_beta_r = obs_buf_env$hs$A_beta_r
        )
        cdid_abort(
          "Nonfinite propensity or outcome partial derivatives.",
          c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
        )
      }
    },
    error = function(cnd) {
      err_diag <- if (exists("error_diagnostics", envir = diag_env, inherits = FALSE)) {
        diag_env$error_diagnostics
      } else {
        list(stage = "unknown_fit_failure")
      }
      cnd$error_diagnostics <- err_diag
      if (!is.null(failure_trace_env)) {
        failure_trace_env$error_diagnostics <- err_diag
        failure_trace_env$error_classes <- class(cnd)
        failure_trace_env$error_message <- conditionMessage(cnd)
        failure_trace_env$observer_released_on_error <- TRUE
      }
      cleanup_temporary_buffers()
      stop(cnd)
    }
  )

  # 6. Extract ONLY compact O(p^2 + p*m + m) quantities, then immediately release
  #    the observer buffer and all per-observer model vectors BEFORE sealing/returning.
  st <- obs_buf_env$st
  hs <- obs_buf_env$hs

  n_treated_obs <- as.integer(sum(st$is_t))
  n_control_obs <- as.integer(sum(st$is_c))
  st_pair_count <- as.integer(st$pair_count)
  names(st_pair_count) <- sup$cohort_labels
  st_pair_mass <- as.double(st$pair_mass)
  names(st_pair_mass) <- sup$cohort_labels
  st_rates <- as.double(st$rates)
  names(st_rates) <- sup$cohort_labels
  st_beta <- as.double(st$beta)
  names(st_beta) <- colnames(frame_owner$Z_a)
  st_sum_a <- as.double(st$sum_a)
  st_mean_a <- as.double(st$mean_a)
  st_sum_b <- as.double(st$sum_b)
  st_mean_b <- as.double(st$mean_b)
  st_fit_info <- st$fit_info
  hs_H <- hs$H
  hs_H_chol <- hs$H_chol
  hs_A_beta_r <- hs$A_beta_r

  # Capture the actual arrays before the existing release; no reconstructed
  # weights, fitted probabilities, outcome scans or additional fit.
  st_weight_concentration <- list(
    propensity_fit_weights = ipw_fitted_diagnostic_summary(st$fit_weights, preflight$effective_chunk_rows, TRUE),
    target_outcome_weights = ipw_fitted_diagnostic_summary(st$a_obs, preflight$effective_chunk_rows, TRUE),
    pooled_control_outcome_ipw_weights = ipw_fitted_diagnostic_summary(st$b_obs, preflight$effective_chunk_rows, TRUE))
  st_probability_range <- ipw_fitted_diagnostic_summary(st$prob_use, preflight$effective_chunk_rows, FALSE)

  # Release st, hs, and all observer arrays NOW!
  st <- NULL
  hs <- NULL
  cleanup_temporary_buffers()

  if (!is.null(callback)) {
    callback(list(
      phase = "after_observer_release",
      active_query_chunks = obs_buf_env$active_query_chunks,
      observer_active = obs_buf_env$active,
      observer_released = obs_buf_env$observer_released,
      endpoints_read = 0L
    ))
  }

  view_pair_mass <- as.double(nuisance_view$support$by_cohort_pair_mass)
  names(view_pair_mass) <- sup$cohort_labels
  view_rates <- as.double(nuisance_view$support$rates)
  names(view_rates) <- sup$cohort_labels
  cohort_a_mass <- as.double(nuisance_view$support$by_cohort_a_mass)
  names(cohort_a_mass) <- sup$cohort_labels
  cohort_q <- as.double(nuisance_view$support$q)
  names(cohort_q) <- sup$cohort_labels

  actual_ops <- stream_ops + as.double(st_fit_info$iter) * (4 * n_obs * p * p + 16 * n_obs * p) +
    4 * (n_obs * p * p + n_obs * p * m) + preflight$nuisance_diagnostics$work_ops_bound
  stopifnot(actual_ops <= preflight$total_work_ops_bound)

  fitted_view <- structure(
    list(
      nuisance_id = nuisance_view$nuisance_id,
      nuisance_view_ref = nuisance_view$ref,
      frame_ref = nuisance_view$frame_ref,
      panel_ref = nuisance_view$panel_ref,
      defining_pair = nuisance_view$defining_pair,
      cohorts = list(
        g = sup$g,
        K = sup$K,
        cohorts = sup$cohorts,
        cohort_labels = sup$cohort_labels,
        n_cohorts = m,
        m_cohorts = m
      ),
      dimensions = list(
        n_a_units = na,
        n_observed_risk_units = n_obs,
        n_covariates = p,
        m_cohorts = m,
        covariate_names = colnames(frame_owner$Z_a)
      ),
      counts = list(
        n_a_units = na,
        n_risk_a_units = as.integer(nuisance_view$support$n_risk_a_units),
        n_outside_risk_a_units = as.integer(nuisance_view$support$n_outside_risk_a_units),
        n_observed_pair_units = as.integer(nuisance_view$support$n_observed_pair_all),
        n_observed_risk_units = n_obs,
        n_nonobserved_risk_units = as.integer(nuisance_view$support$n_nonobserved_risk_units),
        n_treated_obs = n_treated_obs,
        n_control_obs = n_control_obs,
        by_cohort_a = nuisance_view$support$by_cohort_a_count,
        by_cohort_pair = st_pair_count
      ),
      masses = list(
        a_weight_total = as.double(frame_owner$sum_w_a),
        cohort_a_weight = cohort_a_mass,
        cohort_pair_weight = st_pair_mass,
        q = cohort_q,
        rates = st_rates,
        view_chunked_pair_weight = view_pair_mass,
        view_chunked_rates = view_rates
      ),
      means = list(
        sum_a = st_sum_a,
        mean_a = st_mean_a,
        sum_b = st_sum_b,
        mean_b = st_mean_b
      ),
      coefficients = st_beta,
      beta = st_beta,
      factors = list(
        H = hs_H,
        H_chol = hs_H_chol,
        qr_R = st_fit_info$qr_R,
        qr_pivot = st_fit_info$qr_pivot,
        qr_tol = st_fit_info$qr_tol,
        rank = st_fit_info$rank
      ),
      partials = list(
        A_beta_r = hs_A_beta_r
      ),
      diagnostics = list(
        weight_concentration = st_weight_concentration,
        weight_concentration_definition = "Kish(v) = (sum(v/max(v)))^2 / sum((v/max(v))^2); descriptive concentration only",
        weight_concentration_sources = c(propensity_fit_weights = "st$fit_weights",
          target_outcome_weights = "st$a_obs", pooled_control_outcome_ipw_weights = "st$b_obs"),
        fitted_probability_range = st_probability_range,
        fitted_probability_scope = "observed fitted risk set",
        converged = st_fit_info$converged,
        boundary = st_fit_info$boundary,
        iter = st_fit_info$iter,
        rank = st_fit_info$rank,
        n_covariates = p,
        df_residual = st_fit_info$df_residual,
        deviance = st_fit_info$deviance,
        null_deviance = st_fit_info$null_deviance,
        aic = st_fit_info$aic,
        control = glm_ctrl,
        clipping = st_fit_info$clipping
      ),
      point_summary_disclosure = list(
        defining_rate_source = "observer_buffer_full_vector_sum_over_frame_owner_cohort_a_mass",
        observer_cohort_a_mass = cohort_a_mass,
        observer_cohort_pair_mass = st_pair_mass,
        observer_rates = st_rates,
        view_chunked_cohort_a_mass = cohort_a_mass,
        view_chunked_cohort_pair_mass = view_pair_mass,
        view_chunked_rates = view_rates,
        max_abs_pair_mass_chunk_diff = max(abs(st_pair_mass - view_pair_mass)),
        max_abs_rate_chunk_diff = max(abs(st_rates - view_rates)),
        candidate_point_summary_only = TRUE,
        certifies_root_or_inference = FALSE,
        glm_convergence_is_certificate = FALSE,
        seal_is_certificate = FALSE,
        disclosure = paste0(
          "Candidate weighted-logit fit uses finite binary64 observer-buffer point summaries ",
          "(observer_cohort_pair_mass / observer_cohort_a_mass summed in canonical A-unit order); ",
          "chunked view-support point summaries (view_chunked_cohort_pair_mass, view_chunked_rates) ",
          "are also recorded. Neither stats::glm.fit convergence nor object sealing certifies a ",
          "mathematical root, outward interval enclosure, or statistical inference."
        )
      ),
      glm_control = glm_ctrl,
      control = ctrl,
      preflight = preflight,
      accounting = list(
        retained_matrix_entries = preflight$retained_matrix_entries_per_fit,
        retained_vector_entries = preflight$retained_vector_entries_per_fit,
        total_retained_entries = preflight$retained_entries_per_fit,
        observer_buffer_released = TRUE,
        reader_chunks_released = TRUE,
        full_n_masks_retained = FALSE,
        observer_arrays_retained = FALSE,
        copied_Z_retained = FALSE,
        endpoints_read = 0L,
        peak_rss_promised = FALSE
      ),
      work = list(
        chunks_queried = logit_cert_count_output(chunks_queried),
        max_simultaneous_query_chunks = 1L,
        max_simultaneous_observer_buffers = 1L,
        observer_units_copied = logit_cert_count_output(n_obs),
        endpoints_read = 0L,
        actual_work_ops = actual_ops
      )
    ),
    class = "cdid_ipw_fitted_nuisance_view"
  )

  ipw_fitted_view_verify_compact_storage(fitted_view)
  ipw_owner_seal(fitted_view, "fitted_nuisance_view")
}

construct_ipw_fitted_nuisance_view <- fit_ipw_nuisance_view

fit_ipw_nuisance_from_owner <- function(
  frame_owner,
  g,
  K,
  s,
  t,
  glm_control,
  control,
  nuisance_id = NULL,
  callback = NULL,
  trace_env = NULL,
  failure_trace_env = NULL
) {
  ipw_fitted_view_require_deps()
  if (missing(frame_owner) || missing(g) || missing(K) || missing(s) || missing(t) ||
      missing(glm_control) || missing(control)) {
    cdid_ipw_fitted_view_abort(
      "All required arguments (`frame_owner`, `g`, `K`, `s`, `t`, `glm_control`, `control`) must be supplied explicitly.",
      "cdid_ipw_fitted_view_input_error"
    )
  }
  ctrl <- validate_ipw_fitted_view_control(control)
  owner_ctrl <- ipw_fitted_view_to_owner_control(ctrl)
  nu_view <- construct_ipw_nuisance_view(
    frame_owner = frame_owner,
    g = g,
    K = K,
    s = s,
    t = t,
    control = owner_ctrl,
    nuisance_id = nuisance_id,
    require_support = TRUE
  )
  # The newly constructed view remains reachable; charge its completed work
  # plus the fit and this separate preflight before entering the fitting stage.
  fit_pre <- preflight_ipw_fitted_nuisance_view(nu_view,glm_control,ctrl)
  ipw_owner_safe_sum(c(nu_view$preflight$work_ops_bound,fit_pre$total_work_ops_bound,
    fit_pre$character_metadata$work_ops_bound,128),ctrl$max_work_ops,"composed owner-view construction and fitting work")
  fit_ipw_nuisance_view(
    nuisance_view = nu_view,
    glm_control = glm_control,
    control = ctrl,
    callback = callback,
    trace_env = trace_env,
    failure_trace_env = failure_trace_env
  )
}

# Sequential multi-view fitting: preflights the set, releases each consumed query
# before the next one, and releases all observer/model arrays before the next
# distinct nuisance fit.
fit_ipw_nuisance_view_set <- function(
  nuisance_views,
  glm_control,
  control,
  trace_env = NULL
) {
  ipw_fitted_view_require_deps()
  if (missing(nuisance_views) || missing(glm_control) || missing(control)) {
    cdid_ipw_fitted_view_abort(
      "Required arguments (`nuisance_views`, `glm_control`, `control`) must be supplied explicitly.",
      "cdid_ipw_fitted_view_input_error"
    )
  }
  if (!is.list(nuisance_views) || length(nuisance_views) < 1L) {
    cdid_ipw_fitted_view_abort(
      "`nuisance_views` must be a non-empty list of sealed `cdid_ipw_nuisance_view` objects.",
      "cdid_ipw_fitted_view_input_error"
    )
  }
  ctrl <- validate_ipw_fitted_view_control(control)
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  n_fits <- length(nuisance_views)

  if(!is.null(trace_env) && !is.environment(trace_env))
    cdid_ipw_fitted_view_abort("`trace_env` must be NULL or an environment.","cdid_ipw_fitted_view_input_error")
  set_header <- ipw_owner_safe_sum(c(ipw_owner_safe_product(16,n_fits),128),
    min(ctrl$max_retained_entries,ctrl$max_matrix_entries,ctrl$max_output_entries),"fit-set list/trace/header payload before allocation")
  set_work <- ipw_owner_safe_sum(c(ipw_owner_safe_product(128,n_fits),128),ctrl$max_work_ops,"fit-set header work")
  sum_ret <- set_header; peak_extra <- 0
  for(i in seq_len(n_fits)) {
    pi <- preflight_ipw_fitted_nuisance_view(nuisance_views[[i]],glm_ctrl,ctrl,n_planned_fits=1L)
    sum_ret <- ipw_owner_safe_sum(c(sum_ret,pi$total_retained_entries),min(ctrl$max_retained_entries,ctrl$max_output_entries),"complete fit-set retained outputs")
    set_work <- ipw_owner_safe_sum(c(set_work,pi$total_work_ops_bound,pi$character_metadata$work_ops_bound,128),ctrl$max_work_ops,"complete fit-set work")
    peak_extra <- max(peak_extra,pi$peak_concurrent_entries-pi$total_retained_entries)
  }
  ipw_owner_safe_sum(c(sum_ret,peak_extra),ctrl$max_matrix_entries,"all fit-set outputs plus largest active workspace")
  pi <- NULL

  fits <- vector("list", n_fits)
  names(fits) <- names(nuisance_views)
  between_fit_active_buffers <- integer(n_fits)

  for (i in seq_len(n_fits)) {
    step_trace <- new.env(parent = emptyenv())
    fits[[i]] <- fit_ipw_nuisance_view(
      nuisance_view = nuisance_views[[i]],
      glm_control = glm_ctrl,
      control = ctrl,
      trace_env = step_trace
    )
    if (!identical(step_trace$active_observer_buffers, 0L) || !isTRUE(step_trace$observer_released)) {
      cdid_ipw_fitted_view_abort(
        "Observer buffer was not released before starting the next distinct nuisance fit.",
        "cdid_ipw_fitted_view_lifecycle_error"
      )
    }
    between_fit_active_buffers[i] <- as.integer(step_trace$active_observer_buffers)
  }

  if (!is.null(trace_env)) {
    if (!is.environment(trace_env)) {
      cdid_ipw_fitted_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_fitted_view_input_error")
    }
    trace_env$n_distinct_fits <- n_fits
    trace_env$between_fit_active_buffers <- between_fit_active_buffers
    trace_env$all_released_between_fits <- all(between_fit_active_buffers == 0L)
  }

  fits
}

# ---------------------------------------------------------------------------
# Seal validation, nuisance sharing verification, and RDS revalidation
# ---------------------------------------------------------------------------

validate_ipw_fitted_nuisance_view <- function(
  fitted_view,
  nuisance_view = NULL,
  frame_owner = NULL
) {
  ipw_fitted_view_require_deps()
  if (missing(fitted_view) || is.null(fitted_view) ||
      !inherits(fitted_view, "cdid_ipw_fitted_nuisance_view") ||
      !is.list(fitted_view) ||
      !ipw_owner_seal_ok(fitted_view, "fitted_nuisance_view")) {
    cdid_ipw_fitted_view_abort(
      "Fitted nuisance view is malformed, unsealed, deserialized without revalidation, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_ipw_fitted_view_integrity_error", "cdid_ipw_owner_integrity_error", "cdid_ipw_owner_view_error")
    )
  }
  ipw_owner_refuse_decorative_attrs(fitted_view, "fitted_nuisance_view")
  ipw_fitted_view_verify_compact_storage(fitted_view)

  nv_ref <- fitted_view$nuisance_view_ref
  if (!is.environment(nv_ref) || !environmentIsLocked(nv_ref) ||
      !exists("nuisance_view", envir = nv_ref, inherits = FALSE) ||
      !bindingIsLocked("nuisance_view", nv_ref)) {
    cdid_ipw_fitted_view_abort(
      "Fitted nuisance view lacks a valid sealed `nuisance_view` reference.",
      c("cdid_ipw_fitted_view_integrity_error", "cdid_ipw_owner_integrity_error")
    )
  }
  nv_obj <- get("nuisance_view", envir = nv_ref, inherits = FALSE)
  validate_ipw_nuisance_view(nv_obj)

  if (!identical(fitted_view$frame_ref, nv_obj$frame_ref) ||
      !identical(fitted_view$panel_ref, nv_obj$panel_ref)) {
    cdid_ipw_fitted_view_abort(
      "Fitted nuisance view owner references disagree with its `nuisance_view`.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_fitted_view_integrity_error")
    )
  }

  f_obj <- get("frame_owner", envir = fitted_view$frame_ref, inherits = FALSE)
  validate_ipw_frame_owner(f_obj)

  if (!is.null(nuisance_view)) {
    validate_ipw_nuisance_view(nuisance_view, frame_owner = f_obj)
    if (!identical(nv_ref, nuisance_view$ref)) {
      cdid_ipw_fitted_view_abort(
        "Fitted nuisance view belongs to a different `cdid_ipw_nuisance_view` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  }
  if (!is.null(frame_owner)) {
    validate_ipw_frame_owner(frame_owner)
    if (!identical(fitted_view$frame_ref, frame_owner$ref)) {
      cdid_ipw_fitted_view_abort(
        "Fitted nuisance view belongs to a different `cdid_ipw_frame_owner` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_context_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  }

  validate_ipw_glm_control(fitted_view$glm_control)
  validate_ipw_fitted_view_control(fitted_view$control)

  p <- as.integer(f_obj$n_covariates)
  m <- as.integer(nv_obj$cohorts$m_cohorts)
  beta <- fitted_view$coefficients
  if (!is.double(beta) || length(beta) != p || any(!is.finite(beta)) ||
      !identical(beta, fitted_view$beta) ||
      !identical(names(beta), colnames(f_obj$Z_a))) {
    cdid_ipw_fitted_view_abort(
      "Fitted nuisance view `coefficients` vector is malformed or inconsistent with `frame_owner`.",
      "cdid_ipw_fitted_view_integrity_error"
    )
  }
  H <- fitted_view$factors$H
  H_chol <- fitted_view$factors$H_chol
  qr_R <- fitted_view$factors$qr_R
  A_beta_r <- fitted_view$partials$A_beta_r
  if (!is.matrix(H) || !identical(dim(H), c(p, p)) || any(!is.finite(H)) ||
      !is.matrix(H_chol) || !identical(dim(H_chol), c(p, p)) || any(!is.finite(H_chol)) ||
      any(diag(H_chol) <= 0) ||
      !is.matrix(qr_R) || !identical(dim(qr_R), c(p, p)) || any(!is.finite(qr_R)) ||
      any(abs(diag(qr_R)) <= 0) ||
      !is.matrix(A_beta_r) || !identical(dim(A_beta_r), c(p, m)) || any(!is.finite(A_beta_r))) {
    cdid_ipw_fitted_view_abort(
      "Fitted nuisance view factor or derivative matrices are malformed.",
      "cdid_ipw_fitted_view_integrity_error"
    )
  }
  psd <- fitted_view$point_summary_disclosure
  if (!is.list(psd) || !isTRUE(psd$candidate_point_summary_only) ||
      !identical(psd$certifies_root_or_inference, FALSE) ||
      !identical(psd$glm_convergence_is_certificate, FALSE) ||
      !identical(psd$seal_is_certificate, FALSE)) {
    cdid_ipw_fitted_view_abort(
      "Fitted nuisance view point-summary non-certification disclosure is missing or altered.",
      "cdid_ipw_fitted_view_integrity_error"
    )
  }

  invisible(fitted_view)
}

# Verify that a fitted nuisance view can be shared with a target link/nuisance view:
# requires identical sealed frame/panel owners, identical ordered cohorts (g, K),
# identical glm_control, and streamed chunk-by-chunk equality of the actual
# risk-set pair-observation mask O = R intersect S. Counts or labels alone never suffice.
verify_ipw_fitted_nuisance_sharing <- function(
  fitted_view,
  target_view,
  glm_control,
  control
) {
  ipw_fitted_view_require_deps()
  if (missing(fitted_view) || missing(target_view) || missing(glm_control) || missing(control)) {
    cdid_ipw_fitted_view_abort(
      "Required arguments (`fitted_view`, `target_view`, `glm_control`, `control`) must be supplied explicitly.",
      "cdid_ipw_fitted_view_input_error"
    )
  }
  validate_ipw_fitted_nuisance_view(fitted_view)
  ctrl <- validate_ipw_fitted_view_control(control)
  owner_ctrl <- ipw_fitted_view_to_owner_control(ctrl)
  glm_ctrl <- validate_ipw_glm_control(glm_control)

  if (!identical(glm_ctrl, fitted_view$glm_control)) {
    cdid_ipw_fitted_view_abort(
      "Caller `glm_control` differs from the fitted nuisance view's `glm_control`; silent nuisance sharing across mismatched solver controls is refused.",
      c("cdid_ipw_control_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  if (inherits(target_view, "cdid_ipw_link_view")) {
    validate_ipw_link_view(target_view)
  } else if (inherits(target_view, "cdid_ipw_nuisance_view")) {
    validate_ipw_nuisance_view(target_view)
  } else {
    cdid_ipw_fitted_view_abort(
      "Expected `target_view` to be a sealed `cdid_ipw_link_view` or `cdid_ipw_nuisance_view`.",
      "cdid_ipw_owner_view_error"
    )
  }

  if (!identical(fitted_view$frame_ref, target_view$frame_ref) ||
      !identical(fitted_view$panel_ref, target_view$panel_ref)) {
    cdid_ipw_fitted_view_abort(
      "Target view belongs to a different frame or panel owner; nuisance sharing is refused.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_context_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  if (!identical(fitted_view$cohorts$g, target_view$cohorts$g) ||
      !identical(fitted_view$cohorts$K, target_view$cohorts$K) ||
      !identical(fitted_view$cohorts$cohorts, target_view$cohorts$cohorts) ||
      !identical(fitted_view$cohorts$cohort_labels, target_view$cohorts$cohort_labels)) {
    cdid_ipw_fitted_view_abort(
      "Target view `g` or ordered control cohorts `K` differ from the fitted nuisance view; sharing is refused.",
      c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  nv_obj <- get("nuisance_view", envir = fitted_view$nuisance_view_ref, inherits = FALSE)
  mask_check <- verify_ipw_view_mask_compatibility(
    nuisance_view = nv_obj,
    link_view = target_view,
    control = owner_ctrl
  )

  invisible(list(
    compatible = TRUE,
    n_a_units = mask_check$n_a_units,
    chunks_compared = mask_check$chunks_compared,
    endpoints_read = 0L
  ))
}

# Extract a compact, reference-free exact binary64 snapshot suitable for RDS serialization
# in small audits. Cannot be used directly as a sealed view without `revalidate_ipw_fitted_nuisance_view`.
extract_ipw_fitted_nuisance_snapshot <- function(fitted_view) {
  validate_ipw_fitted_nuisance_view(fitted_view)
  structure(
    list(
      nuisance_id = fitted_view$nuisance_id,
      defining_pair = fitted_view$defining_pair,
      cohorts = fitted_view$cohorts,
      dimensions = fitted_view$dimensions,
      counts = fitted_view$counts,
      masses = fitted_view$masses,
      means = fitted_view$means,
      coefficients = fitted_view$coefficients,
      beta = fitted_view$beta,
      factors = fitted_view$factors,
      partials = fitted_view$partials,
      diagnostics = fitted_view$diagnostics,
      point_summary_disclosure = fitted_view$point_summary_disclosure,
      glm_control = fitted_view$glm_control,
      control = fitted_view$control
    ),
    class = "cdid_ipw_fitted_nuisance_snapshot"
  )
}

# Bound a reference-free snapshot's reachable payload before recursive checks/refit.
# Environment references are borrowed identities and are never traversed here.
ipw_fitted_snapshot_payload <- function(object, max_entries, max_work) {
  entries <- 0; work <- 0
  visit <- function(x,depth=0L) {
    if(depth>4L) cdid_ipw_fitted_view_abort("Snapshot exceeds compact schema nesting.","cdid_ipw_fitted_view_revalidation_error")
    if(is.environment(x) || is.null(x)) return(invisible(NULL))
    k <- length(x)
    entries <<- ipw_owner_safe_sum(c(entries,k,1),max_entries,"snapshot reachable payload before content checks")
    work <<- ipw_owner_safe_sum(c(work,k,1),max_work,"snapshot traversal work")
    if(!is.null(names(x))) visit(names(x),depth+1L)
    if(is.list(x)) for(j in seq_len(k)) visit(x[[j]],depth+1L)
    else if(is.character(x)) for(j in seq_len(k)) {
      nb <- nchar(x[j],type="bytes")
      if(is.na(nb)) nb <- 0
      entries <<- ipw_owner_safe_sum(c(entries,ceiling(nb/2)),max_entries,"snapshot character payload")
      work <<- ipw_owner_safe_sum(c(work,nb),max_work,"snapshot character work")
    }
    invisible(NULL)
  }
  visit(object)
  list(entries=entries,work=work)
}

# Reconstruct and reseal a deserialized `cdid_ipw_fitted_nuisance_view` or
# `cdid_ipw_fitted_nuisance_snapshot` against a live sealed `nuisance_view`.
revalidate_ipw_fitted_nuisance_view <- function(
  object,
  nuisance_view,
  glm_control,
  control
) {
  ipw_fitted_view_require_deps()
  if (missing(object) || missing(nuisance_view) || missing(glm_control) || missing(control)) {
    cdid_ipw_fitted_view_abort(
      "Required arguments (`object`, `nuisance_view`, `glm_control`, `control`) must be supplied explicitly.",
      "cdid_ipw_fitted_view_input_error"
    )
  }
  ipw_owner_refuse_decorative_attrs(object, "object")
  if (!is.list(object) ||
      (!inherits(object, "cdid_ipw_fitted_nuisance_view") &&
       !inherits(object, "cdid_ipw_fitted_nuisance_snapshot"))) {
    cdid_ipw_fitted_view_abort(
      "Expected a `cdid_ipw_fitted_nuisance_view` or `cdid_ipw_fitted_nuisance_snapshot` for revalidation.",
      c("cdid_ipw_fitted_view_revalidation_error", "cdid_ipw_owner_integrity_error")
    )
  }
  snapshot_fields <- c("nuisance_id","defining_pair","cohorts","dimensions","counts","masses","means","coefficients","beta","factors","partials","diagnostics","point_summary_disclosure","glm_control","control")
  full_fields <- c(snapshot_fields,"nuisance_view_ref","frame_ref","panel_ref","preflight","accounting","work","ref")
  expected <- if(inherits(object,"cdid_ipw_fitted_nuisance_snapshot")) snapshot_fields else full_fields
  if(length(object)!=length(expected) || is.null(names(object)) || anyDuplicated(names(object)) || !setequal(names(object),expected))
    cdid_ipw_fitted_view_abort("Snapshot fields differ from the exact compact schema.",c("cdid_ipw_fitted_view_revalidation_error","cdid_ipw_owner_integrity_error"))
  validate_ipw_nuisance_view(nuisance_view)
  glm_ctrl <- validate_ipw_glm_control(glm_control)
  ctrl <- validate_ipw_fitted_view_control(control)

  if (!identical(object$glm_control, glm_ctrl)) {
    cdid_ipw_fitted_view_abort(
      "Deserialized nuisance snapshot `glm_control` disagrees with caller `glm_control`.",
      c("cdid_ipw_fitted_view_revalidation_error", "cdid_ipw_control_mismatch_error")
    )
  }

  if(!identical(object$control,ctrl))
    cdid_ipw_fitted_view_abort("Snapshot resource controls differ from explicit revalidation controls.",c("cdid_ipw_fitted_view_revalidation_error","cdid_ipw_control_mismatch_error"))
  pre <- preflight_ipw_fitted_nuisance_view(nuisance_view,glm_ctrl,ctrl)
  payload <- ipw_fitted_snapshot_payload(object,
    max_entries=max(0,ctrl$max_matrix_entries-pre$peak_concurrent_entries),
    max_work=max(0,ctrl$max_work_ops-pre$total_work_ops_bound-pre$character_metadata$work_ops_bound-128))
  ipw_fitted_view_verify_compact_storage(object)

  # Recompute the live owner-backed fit to verify every exact binary64 array and metadata field
  live_fit <- fit_ipw_nuisance_view(
    nuisance_view = nuisance_view,
    glm_control = glm_ctrl,
    control = ctrl
  )

  check_fields <- c(
    "nuisance_id", "defining_pair", "cohorts", "dimensions", "counts", "masses", "means",
    "coefficients", "beta", "factors", "partials", "diagnostics",
    "point_summary_disclosure", "glm_control"
  )
  for (nm in check_fields) {
    if (!identical(object[[nm]], live_fit[[nm]])) {
      cdid_ipw_fitted_view_abort(
        paste0("Deserialized nuisance field `", nm, "` differs from live owner-backed recomputation; revalidation is refused."),
        c("cdid_ipw_fitted_view_revalidation_error", "cdid_ipw_owner_integrity_error"),
        diagnostics = list(field = nm)
      )
    }
  }

  live_fit
}
