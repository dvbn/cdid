# Optional complete multioperator reuse only. Validate sources/controls outside
# the count-only capacity handler. Nothing below fits a nuisance or reads Y.
# A declined prospective bound selects the preexisting route before fitting;
# fitted/numerical/integrity failures are never caught for another producer.
ipw_pipeline_complete_raw_capacity <- function(metadata, frame, ctrl, aux,
    retained, work, early_reduction, weighting, J_bound, R_full, R_rep, n_post_bound) {
  validate_owner_comparison_catalog(metadata, frame_owner = frame)
  validate_ipw_frame_owner(frame)
  ctrl <- validate_ipw_pipeline_control(ctrl, weighting)
  panel <- get("panel_owner", metadata$panel_ref, inherits = FALSE)
  N <- as.double(frame$n_a_units); p <- as.double(frame$n_covariates)
  M <- as.double(nrow(metadata$catalog$rows)); Q <- as.double(length(aux$descriptors))
  C <- as.double(length(frame$a_cohorts)); T <- as.double(panel$n_calendar_dates)
  O <- as.double(nrow(metadata$catalog$operators)); D <- M + Q
  if (frame$n_total_units + frame$n_frame_only_units != N ||
      any(panel$postings$date_count != N) || O < 2) return(list(admitted = FALSE))
  capacity_decline <- function(e) {
    # Bound classes also include malformed dimensions and exact-count
    # overflow. Only a finite annotated count above a finite cap is capacity.
    d <- e$diagnostics
    count <- function(x) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
      length(x) == 1L && is.finite(x) && x >= 0 && x == floor(x) && x <= 2^53 - 1
    over <- function(x, cap) count(x) && count(cap) && x > cap
    genuine <- !inherits(e, c("cdid_ipw_owner_control_error", "cdid_ipw_owner_matrix_dim_error", "cdid_gmm_control_error")) &&
      is.list(d) && (over(d$requested, d$bound) || over(d$requested_entries, d$max_matrix_entries) ||
        over(d$total_payload_entries, d$max_matrix_entries))
    if (!genuine) stop(e)
    list(admitted = FALSE, reason = conditionMessage(e))
  }
  tryCatch({
    own <- ctrl$owner_control; fc <- ctrl$fitted_control; pc <- ctrl$point_control
    cap <- min(own$max_matrix_entries, own$max_retained_entries,
      fc$max_matrix_entries, pc$max_matrix_entries)
    B <- min(N, own$max_chunk_rows, pc$max_chunk_rows)
    if (N > min(fc$max_observer_units, pc$max_observer_units)) return(list(admitted = FALSE))
    chars <- ipw_raw_record_char_inventory(metadata, frame, aux, fc,
      retained = max(retained, ipw_point_view_owner_entries(frame, panel)), envelope = own, work = work)
    # Keep numeric tables and factor 4; local character roles replace the
    # all-catalog-strings charge repeated for every D record.
    tables <- ipw_owner_safe_sum(c(ipw_owner_safe_product(D,
      20480 + 16 * p^2 + 16 * p * C + 256 * (p + C + 1)),
      ipw_owner_safe_product(4, chars$total_record_slots), chars$metadata_entries))
    op_rows <- as.double(tabulate(match(metadata$catalog$rows$operator_id,
      unique(metadata$catalog$rows$operator_id))))
    masks <- 16 * (M + 1) + length(metadata$catalog$spec$cohorts)
    splits <- 8 * sum(op_rows * (op_rows - 1) / 2) +
      M * (11 + length(metadata$catalog$spec$cohorts)) + 128 * (M + 1)
    descriptor <- pipeline_metadata_descriptor_preflight(metadata, own, retained, work, M, memo_entries = 64)
    cert <- cert_work <- E <- 0
    if (early_reduction) {
      cb <- validate_ipw_certificate_bounds(ctrl$class_bounds)
      E <- min(cb$max_E_entries, M * min(M, T - 1, cb$max_class_dates - 1))
      cert <- 5 * E + 256 * M + 64 * nrow(metadata$design$vertices) + 256 +
        4 * M * length(panel$cohorts)
      cert_work <- 512 * M * (length(panel$cohorts) + 2) +
        256 * M * (T + (O + 1) * (length(panel$cohorts) + 1) + ceiling(log2(M + 1))) + 64 * E
    }
    # All complete pairs have one support mask within an actual catalog
    # operator. A forest therefore has at most T-1 rows per operator, even
    # when its date graph is disconnected. This sizes admission ONLY. The
    # existing exact class certificate/E/kept metric still runs unchanged.
    L <- if (early_reduction) min(M, O * (T - 1)) else M
    U <- min(D, O + Q)
    nu <- preflight_ipw_nuisance_view(N, p, C, own)
    lv <- preflight_ipw_link_view(N, p, C, own)
    dummy <- list(nuisance_id = c(paste0("nuisance:op", metadata$catalog$rows$operator_id, ":r", seq_len(M)),
        paste0("reporting_nuisance:", seq_len(Q))),
      cohorts = list(cohort_labels = as.character(frame$a_cohorts)), preflight = nu)
    fp <- preflight_ipw_fitted_nuisance_view_dimensions(N, N, p, C, fc, dummy, frame, ctrl$glm_control)
    pp <- preflight_ipw_point_link_view_dimensions(N, N, p, C, pc, NULL, NULL)
    leaf <- max(nu$concurrent_entries, lv$concurrent_entries, fp$concurrent_entries,
      fp$observer_buffer_entries + fp$worker_workspace_entries + pp$concurrent_entries)
    snapshot <- unname(cp04_payload(list(owner_snapshot = list(frame = NULL, panel = NULL)))["entries"])
    persistent <- ipw_owner_safe_sum(c(retained, tables, masks, splits, cert,
      descriptor$retained_entries, snapshot))
    construction <- max(ipw_owner_safe_sum(c(persistent, descriptor$scratch_entries, leaf)),
      ipw_owner_safe_sum(c(max(retained, ipw_point_view_owner_entries(frame, panel)), chars$scratch_entries)))
    # The returned backend/diagnostics are a subset of the construction tables.
    # Retain that entire conservative inventory here, then compose every
    # prospective consumer explicitly, including GMM, D5 and auxiliary output.
    # These explicit count checks give annotated capacity declines where
    # the GMM preflight's combined QR/materialization guards have no counts.
    k <- nrow(metadata$design$parameters)
    ipw_owner_safe_sum(c(ipw_owner_safe_product(L, k), k), ctrl$gmm_bounds$max_matrix_entries,
      "optional complete raw-IPW compact QR")
    ipw_owner_safe_product(N, J_bound, min(ctrl$gmm_bounds$max_matrix_entries,
      ctrl$gmm_bounds$max_materialized_if_entries), "optional complete raw-IPW target IF")
    ipw_owner_safe_sum(J_bound, ctrl$gmm_bounds$max_targets, "optional complete raw-IPW target count")
    gp <- preflight_gmm_allocations(L, nrow(metadata$design$parameters), J_bound,
      weighting, ctrl$gmm_bounds, n_out_rows = N, request_param_cov = TRUE,
      materialize_target_if = TRUE, n_input_rows = N,
      n_parameters = nrow(metadata$design$parameters),
      certify_error = !identical(ctrl$gmm_accuracy$method, "conventional_chol_qr"),
      covariance_strategy = ctrl$gmm_accuracy$covariance_strategy, n_dropped_triplets = E)
    request_L <- max(1, L, Q)
    env <- ipw_raw_if_envelope(B, request_L * C, request_L, request_L, p, C,
      P = min(request_L, T * (T - 1) / 2))
    # Retained-basis covariance prefix, all four centered Gram buffers, row
    # copies/plan indices, D5 ledger/basis and complete reporting/transport
    # arrays. Optional IF stores use their existing independent capacity gates.
    consumers <- 16 * L^2 + 32 * B * (L + p + 1) + 4 * L^2 + 8 * B * L +
      env$output + env$scratch + 4096 + 64 * (D * C + D + p) +
      gp$total_payload_entries + N * (R_full + n_post_bound) +
      8 * N * R_rep + 16 * R_rep^2
    peak <- ipw_owner_safe_sum(c(persistent, consumers))
    ipw_owner_safe_sum(c(max(construction, peak)), cap, "optional complete raw-IPW whole-path capacity")
    ipw_owner_safe_sum(c(fp$retained_entries_per_fit), min(fc$max_output_entries, fc$max_retained_entries))
    ipw_owner_safe_sum(c(pp$retained_entries_per_link), min(pc$max_output_entries, pc$max_retained_entries))
    ipw_owner_safe_sum(c(env$output), min(own$max_output_entries, pc$max_output_entries))
    # Complete masks never split an operator branch. Auxiliary comparisons may
    # add at most Q branches; all pair/mask scans still retain their original
    # counted work, including exact cross-operator comparison attempts.
    stream <- 128 * N * (1 + ceiling(log2(C + 1))) +
      128 * ceiling(N / B) * (1 + ceiling(log2(N + 1)))
    scans <- 16 * (D + T * (T - 1) / 2 + Q * U) * stream
    extra_work <- 2 * chars$work_ops_bound + scans + U * fp$work_ops_bound + D * pp$work_ops_bound + cert_work +
      descriptor$work_ops_bound + gp$covariance_strategy_work_ops +
      8 * ceiling(N / B) * env$work + 64 * N * (L^2 + J_bound^2 + 1)
    ipw_owner_safe_sum(c(work, extra_work), own$max_work_ops, "optional complete raw-IPW whole-path work")
    list(admitted = TRUE, peak_entries = max(construction, peak),
      persistent_entries = persistent, construction_entries = construction,
      consumer_entries = consumers, retained_links_upper = L, mandatory_work_upper = work + extra_work,
      character_inventory = chars, construction_table_entries_bound = tables)
  }, cdid_ipw_owner_bound_error = capacity_decline, cdid_gmm_bound_error = capacity_decline)
}
