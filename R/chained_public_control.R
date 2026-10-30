# D23 public CP02 candidate: conventional Cholesky/QR estimation and complete
# IF/inference; complete per-fit error certification is not requested. The
# fourteen-control private interface remains available for strict research.
# D27 adopted finite resource defaults and retained advanced overrides. These
# count numeric payload/work envelopes, without RSS, CPU or whole-call guarantees.
# D78 (2026-10-24) supersedes the D27 finite WORK default only: the public
# estimated-work cap is OPT-IN. `control = NULL` and an empty list carry no user
# work budget (`counted_work_limit = Inf`); an explicit positive finite exact
# whole number remains a user cap, and an explicit `Inf` also opts out. The
# memory (1 GiB), worker, chunk, numerical, identification, support and model
# policies are unchanged. Preflight work ESTIMATES stay diagnostic/conservative.

chained_public_abort <- function(message, class = character()) {
  cdid_abort(message, unique(c(class, "cdid_chained_public_control_error", "cdid_chained_request_error")))
}

chained_public_need <- function(ok, message, class = character()) {
  if (!isTRUE(ok)) chained_public_abort(message, class)
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# 1. Resource/accuracy profile
# ------------------------------------------------------------------------------

chained_public_profile_status <- "ADOPTED_D78_OPTIONAL_WORK_CAP_D27_MEMORY_DEFAULTS"

# Adopted defaults: D27 memory/worker defaults, D78 optional (disabled) work
# cap; D23 still defers universal per-fit accuracy. Historical evidence labelled
# ADOPTED_D27_FINITE_RESOURCE_DEFAULTS keeps its finite 1e11 work profile.
chained_public_profile_defaults <- function() {
  list(numeric_workspace_bytes = 1073741824, counted_work_limit = Inf, workers = 1)
}

# The private constructors keep their finite integer contracts (every count in
# [1, 2^53 - 1]). When the public work cap is disabled, the PUBLIC-to-private
# seam alone hands them the existing exact-count arithmetic ceiling 2^53 - 1,
# which every private sum/product already enforces: representability
# protection, not a replacement user resource budget. Private count arithmetic
# and array-index overflow rules are unchanged.
chained_public_work_arithmetic_ceiling <- function() 2^53 - 1

# Public `counted_work_limit`: `Inf` (no user work budget) or a finite positive
# exact whole number <= 2^53 - 1. Zero, negatives, -Inf, NA/NaN, non-integral,
# out-of-exact-range, non-numeric, object, dimensioned or non-scalar values refuse.
chained_public_work_limit_ok <- function(x) {
  is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) == 1L && !is.na(x) &&
    ((is.infinite(x) && x > 0) || (is.finite(x) && x >= 1 && x == floor(x) && x <= 2^53 - 1))
}

# Development settings carried unchanged from the accepted private example
# (inst/scripts/cp02_pipeline_dev.R): nuisance solver control, logit certificate
# precision and the GMM core's relative conditioning predicates. They are
# recorded as candidate settings; none is a tolerance search or a public policy.
chained_public_carried_settings <- function() {
  list(
    glm = list(epsilon = 1e-12, maxit = 50L),
    certificate = list(taylor_degree = 12L, max_abs_eta = 50.0, max_squarings = 16L),
    gmm_guard = list(max_orthogonality_defect = 1e-8, min_certified_lambda_rel = 1e-10, max_factor_residual_rel = 1e-10,
                     min_rcond_R = 1e-10, min_rcond_U = 1e-10, min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10,
                     max_loading_residual_rel = 1e-10),
    chunk_rows_cap = 256L
  )
}

chained_public_control_names <- c("numeric_workspace_bytes", "counted_work_limit", "workers",
                                  "chunk_rows", "condition_trigger", "ridge_relative")

# Legacy constructor fields retained for private control compatibility only.
# The conventional method does not interpret these as public error budgets.
chained_public_staging_envelope <- function() .Machine$double.xmax

chained_public_resolve_control <- function(control, gmm_weighting = "optimal", gmm_covariance = "tikhonov") {
  base <- chained_public_profile_defaults()
  carried <- chained_public_carried_settings()
  w_choice <- resolve_gmm_weighting_choice(
    gmm_weighting, arg = "gmm_weighting",
    error_fn = chained_public_abort, error_class = "cdid_chained_public_control_error"
  )
  cov_choice <- resolve_gmm_covariance_choice(
    gmm_covariance, arg = "gmm_covariance",
    error_fn = chained_public_abort, error_class = "cdid_chained_public_control_error"
  )
  if (is.null(control)) {
    tuning <- resolve_gmm_covariance_tuning(
      NULL, error_fn = chained_public_abort, error_class = "cdid_chained_public_control_error"
    )
    cov_strat <- gmm_covariance_strategy(
      method = cov_choice,
      condition_trigger = tuning$condition_trigger,
      ridge_relative = tuning$ridge_relative,
      overrides = tuning$overrides,
      weighting_requested = w_choice$requested,
      weighting_canonical = w_choice$canonical
    )
    return(list(profile = base, profile_status = chained_public_profile_status,
      source = "NULL: resolved to the defaults (optional work cap disabled; memory/worker defaults)",
      overrides = character(0L), relative_accuracy_source = "deferred_by_D23", relative_accuracy_override_recorded = FALSE,
      glm = carried$glm, certificate = carried$certificate, gmm_guard = carried$gmm_guard,
      chunk_rows = NULL, chunk_rows_cap = carried$chunk_rows_cap,
      gmm_weighting = w_choice, gmm_covariance = cov_choice,
      covariance_tuning = tuning, covariance_strategy = cov_strat))
  }
  chained_public_need(is.list(control) && !is.object(control) && !is.data.frame(control),
    paste0("`control` must be NULL or a plain named list of advanced overrides (", paste(chained_public_control_names, collapse = ", "),
           "). Objects, including the private fourteen-control constructor output, are refused here; the private development entry ",
           "point is chained_att_gt_request()."), "cdid_chained_public_control_error")
  nms <- names(control)
  chained_public_need(length(control) == 0L || (!is.null(nms) && !anyNA(nms) && all(nzchar(nms)) && !anyDuplicated(nms)),
    "`control` entries must be uniquely named.", "cdid_chained_public_control_error")
  if ("relative_accuracy" %in% nms) {
    chained_public_abort("`control$relative_accuracy` and the universal 1e-6-of-SE guarantee are not implemented; no public accuracy override is available.",
      "cdid_chained_public_accuracy_deferred")
  }
  unknown <- setdiff(nms, chained_public_control_names)
  chained_public_need(!length(unknown), paste0("Unknown `control` entr", if (length(unknown) == 1L) "y" else "ies", ": ",
    paste(unknown, collapse = ", "), ". Allowed: ", paste(chained_public_control_names, collapse = ", "), "."),
    "cdid_chained_public_control_error")
  whole <- function(x, min) is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= min && x == floor(x)
  prof <- base; overrides <- character(0L)
  if ("numeric_workspace_bytes" %in% nms) {
    chained_public_need(whole(control$numeric_workspace_bytes, 8) && control$numeric_workspace_bytes <= 2^53 - 1,
      "`control$numeric_workspace_bytes` must be a finite whole number of bytes >= 8 (exact binary64 integer).", "cdid_chained_public_control_error")
    prof$numeric_workspace_bytes <- as.double(control$numeric_workspace_bytes); overrides <- c(overrides, "numeric_workspace_bytes")
  }
  if ("counted_work_limit" %in% nms) {
    # D78: an explicit finite exact whole number is a user work cap; an explicit
    # Inf opts out (the default). Nothing else is accepted.
    chained_public_need(chained_public_work_limit_ok(control$counted_work_limit),
      paste0("`control$counted_work_limit` must be Inf (no user work budget; the default) or a finite positive whole number ",
             "of counted work units (exact binary64 integer <= 2^53 - 1)."), "cdid_chained_public_control_error")
    prof$counted_work_limit <- as.double(control$counted_work_limit); overrides <- c(overrides, "counted_work_limit")
  }
  rel_source <- "deferred_by_D23"

  # `workers` declares how many multiplier-bootstrap worker processes a request
  # may use (inherited pl/cores); estimation itself is always serial. The D27
  # default stays 1 (unchanged by D78): `pl = TRUE, cores = k` needs an explicit workers >= k.
  if ("workers" %in% nms) {
    chained_public_need(whole(control$workers, 1) && control$workers <= .Machine$integer.max,
      "`control$workers` must be a positive whole number of declared multiplier bootstrap worker processes.",
      "cdid_chained_public_control_error")
    prof$workers <- as.double(control$workers); overrides <- c(overrides, "workers")
  }
  chunk <- NULL
  if ("chunk_rows" %in% nms) {
    chained_public_need(whole(control$chunk_rows, 1) && control$chunk_rows <= .Machine$integer.max,
      "`control$chunk_rows` must be a positive whole number of units per chunk.", "cdid_chained_public_control_error")
    chunk <- as.double(control$chunk_rows); overrides <- c(overrides, "chunk_rows")
  }
  tuning <- resolve_gmm_covariance_tuning(
    control[intersect(names(control), gmm_covariance_tuning_names)],
    error_fn = chained_public_abort, error_class = "cdid_chained_public_control_error"
  )
  overrides <- c(overrides, tuning$overrides)
  # Solver, certificate and GMM guards are carried private settings, not public overrides.
  glm <- carried$glm; cert <- carried$certificate; guard <- carried$gmm_guard
  if ("condition_trigger" %in% tuning$overrides) {
    guard$min_certified_lambda_rel <- min(guard$min_certified_lambda_rel, tuning$condition_trigger)
    guard$min_rcond_R <- min(guard$min_rcond_R, sqrt(tuning$condition_trigger))
  }
  cov_strat <- gmm_covariance_strategy(
    method = cov_choice,
    condition_trigger = tuning$condition_trigger,
    ridge_relative = tuning$ridge_relative,
    overrides = tuning$overrides,
    weighting_requested = w_choice$requested,
    weighting_canonical = w_choice$canonical
  )
  list(profile = prof, profile_status = chained_public_profile_status,
    source = if (length(overrides)) "explicit advanced list over the defaults (optional work cap; memory/worker defaults)" else
      "empty advanced list: defaults (optional work cap disabled; memory/worker defaults)",
    overrides = overrides, relative_accuracy_source = rel_source,
    relative_accuracy_override_recorded = FALSE,
    glm = glm, certificate = cert, gmm_guard = guard, chunk_rows = chunk, chunk_rows_cap = carried$chunk_rows_cap,
    gmm_weighting = w_choice, gmm_covariance = cov_choice,
    covariance_tuning = tuning, covariance_strategy = cov_strat)
}

# ------------------------------------------------------------------------------
# 2. Actual input dimensions and finite derived upper bounds
# ------------------------------------------------------------------------------
# All products/sums go through the accepted exact-integer helpers, which refuse
# any count that leaves the exact binary64 integer range (overflow-safe).

chained_public_count <- function(x, name) {
  chained_public_need(is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 0 && x == floor(x) && x <= 2^53 - 1,
    paste0("Dimension `", name, "` must be a finite nonnegative whole number."), "cdid_chained_public_dimension_error")
  as.double(x)
}

# Counts -> derived bounds. Every bound is a finite upper bound implied by the
# accepted catalog/design construction, not a toy constant:
#   pairs_per_target  = T (T - 1) / 2 stored date pairs (s < t);
#   candidates_bound  = n_targets * pairs * (1 + T): the catalog enumerates one
#                       original candidate plus at most T fixed-horizon
#                       re-evaluations per supported target pair;
#   links_bound       = n_targets * pairs * n_cohorts: retained rows are
#                       deduplicated per (target, pair, actual control pool K)
#                       and K is an upper suffix of the sorted cohort universe
#                       (plus never-treated), so at most n_cohorts pools per pair;
#   Q_bound           = 2 * n_targets * T auxiliary reporting descriptors (both
#                       base-period modes), as in the pipeline envelope;
#   R_full            = declared post ledger cells (sum over targets of stored
#                       dates >= g), the pipeline's own ledger size;
#   params_bound      = 2 * links_bound + n_targets design vertices (two
#                       endpoints per link plus one reference per cohort);
#   J_bound           = params_bound + Q_bound + 2 R_full + n_targets, the
#                       pipeline's `J_bound` form with every term at its bound.
chained_public_dimensions_from_counts <- function(n_frame, n_panel_units, n_rows, n_dates, n_cohorts, n_targets,
                                                  n_covariates, biters, n_clusters, calendar = NULL, target_cohorts = NULL,
                                                  n_aggregation_types = 4) {
  n_frame <- chained_public_count(n_frame, "n_frame"); n_panel_units <- chained_public_count(n_panel_units, "n_panel_units")
  n_rows <- chained_public_count(n_rows, "n_rows"); n_dates <- chained_public_count(n_dates, "n_dates")
  n_cohorts <- chained_public_count(n_cohorts, "n_cohorts"); n_targets <- chained_public_count(n_targets, "n_targets")
  n_covariates <- chained_public_count(n_covariates, "n_covariates"); biters <- chained_public_count(biters, "biters")
  n_clusters <- chained_public_count(n_clusters, "n_clusters"); n_agg <- chained_public_count(n_aggregation_types, "n_aggregation_types")
  chained_public_need(n_frame >= 1 && n_panel_units >= 1 && n_panel_units <= n_frame && n_rows >= 1 && n_dates >= 1 &&
    n_cohorts >= 1 && n_targets >= 1 && n_targets <= n_cohorts && n_covariates >= 1 && biters >= 1 && n_clusters >= 1 && n_clusters <= n_frame,
    "Inconsistent public input dimensions.", "cdid_chained_public_dimension_error")
  P <- function(a, b, label) ipw_owner_safe_product(a, b, max_entries = NULL, label = label)
  S <- function(..., label) ipw_owner_safe_sum(c(...), max_bound = NULL, label = label)
  pairs <- P(n_dates, n_dates - 1, "stored date pairs") / 2
  R_full <- if (!is.null(calendar) && !is.null(target_cohorts)) {
    S(vapply(as.double(target_cohorts), function(g) sum(as.double(calendar) >= g), numeric(1L)), label = "declared post ledger cells")
  } else P(n_targets, n_dates, "declared post ledger cells (bound)")
  candidates <- P(P(n_targets, pairs, "target pairs"), 1 + n_dates, "catalog candidates")
  links <- P(P(n_targets, pairs, "target pairs"), n_cohorts, "retained catalog rows")
  Q <- P(2, P(n_targets, n_dates, "target-date grid"), "auxiliary reporting descriptors")
  params <- S(P(2, links, "link endpoints"), n_targets, label = "design parameters")
  J <- S(params, Q, R_full, R_full, n_targets, label = "target columns J_bound")
  list(
    n_frame = n_frame, n_panel_units = n_panel_units, n_rows = n_rows, n_dates = n_dates, n_cohorts = n_cohorts,
    n_targets = n_targets, n_covariates = n_covariates, biters = biters, n_clusters = n_clusters, n_aggregation_types = n_agg,
    pairs_per_target = pairs, R_full = R_full, candidates_bound = candidates, links_bound = links, Q_bound = Q,
    params_bound = params, J_bound = J,
    mask_entries_bound = P(links, n_frame, "mask verification entries"),
    pool_entries_bound = P(links, n_cohorts, "catalog pool entries"),
    E_entries_bound = P(links, n_dates, "reduction map entries"),
    materialized_if_bound = P(n_frame, J, "materialized target IF entries"),
    scope = "actual counts of the one prepared input and finite derived upper bounds implied by the accepted catalog/design construction; not toy constants"
  )
}

# From the prepared public input (`chained_request_prepare_inputs`) and the request.
chained_public_dimensions <- function(prep, biters, clustered = FALSE) {
  n_frame <- prep$n_frame
  n_clusters <- if (!is.null(prep$cluster_full)) length(unique(prep$cluster_full)) else n_frame
  dims <- chained_public_dimensions_from_counts(
    n_frame = n_frame, n_panel_units = n_frame - prep$n_frame_only, n_rows = nrow(prep$data_internal),
    n_dates = length(prep$calendar), n_cohorts = max(1L, length(prep$cohort_universe)),
    n_targets = length(prep$target_cohorts),
    n_covariates = max(ncol(prep$frame_Z), if (is.null(prep$cp04)) 0L else ncol(prep$cp04$Z),
      if (is.null(prep$reference_source$p_bound)) 0L else prep$reference_source$p_bound),
    biters = biters, n_clusters = n_clusters,
    calendar = prep$calendar, target_cohorts = prep$target_cohorts)
  if (!is.null(prep$reference_source) && identical(prep$reference_source$design_kind, "dated")) {
    dims$reference_design_width <- "resolved and separately charged on canonical required rows by the sealed reader"
    dims$n_covariates_meaning <- "identity-carrier width for initial chunk planning, not the scientific design width"
  }
  dims
}

# ------------------------------------------------------------------------------
# 3. Fourteen controls from dimensions and the resolved profile
# ------------------------------------------------------------------------------

# Units per streamed chunk: derived from the workspace and the link bound so that
# one chunk's charged slice (32 entries per unit per link/covariate column, as
# in the pipeline envelope) uses at most 1/64 of the workspace; never more than
# the unit frame, never more than the carried cap (256), at least 1. An explicit
# `control$chunk_rows` replaces the rule (still capped by the unit frame).
chained_public_chunk_rows <- function(dims, entries, resolved) {
  if (!is.null(resolved$chunk_rows)) return(max(1, min(dims$n_frame, resolved$chunk_rows)))
  per_unit <- ipw_owner_safe_product(2048, dims$links_bound + dims$n_covariates + 1, label = "chunk slice entries per unit")
  max(1, min(dims$n_frame, resolved$chunk_rows_cap, floor(entries / per_unit)))
}

# `workers` > 1 (an admitted pl = TRUE, cores request) is carried as the
# multiplier's declared worker count; serial budgets are unchanged.
chained_public_multiplier_budget <- function(dims, biters, chunk, entries, bytes, workers = 1) {
  biters <- chained_public_count(biters, "biters")
  sign_rows <- max(1, min(biters, chunk))
  mb <- list(max_units = dims$n_frame, max_clusters = dims$n_clusters, max_targets = dims$J_bound, max_draws = biters,
       max_if_batch_rows = chunk, max_if_batch_cells = ipw_owner_safe_product(chunk, dims$J_bound, label = "IF batch cells"),
       max_sign_batch_rows = sign_rows, max_sign_batch_cells = ipw_owner_safe_product(sign_rows, dims$n_clusters, label = "sign batch cells"),
       max_workspace_cells = entries, max_workspace_bytes = bytes)
  workers <- chained_public_count(workers, "workers")
  chained_public_need(workers >= 1, "The multiplier worker count must be at least one.", "cdid_chained_public_control_error")
  if (workers > 1) mb$workers <- workers
  mb
}

chained_public_pipeline_control <- function(dims, resolved, gmm_weighting = NULL, gmm_covariance = NULL, workers = 1) {
  prof <- resolved$profile
  bytes <- chained_public_count(prof$numeric_workspace_bytes, "numeric_workspace_bytes")
  entries <- floor(bytes / 8)
  chained_public_need(entries >= 1, "The numeric workspace must hold at least one binary64 entry.", "cdid_chained_public_control_error")
  # D78 PUBLIC-to-private seam. `user_work` is the transparent public value
  # (Inf = no user work budget). The finite private constructors receive either
  # the explicit finite cap unchanged or, when the cap is disabled, the existing
  # exact-count arithmetic ceiling 2^53 - 1 (representability protection only).
  user_work <- prof$counted_work_limit
  chained_public_need(chained_public_work_limit_ok(user_work),
    "The resolved `counted_work_limit` must be Inf or a finite positive exact whole number.", "cdid_chained_public_control_error")
  work_cap_enabled <- is.finite(user_work)
  work <- if (work_cap_enabled) chained_public_count(user_work, "counted_work_limit") else chained_public_work_arithmetic_ceiling()
  work_cap <- list(
    status = if (work_cap_enabled) "enabled" else "disabled",
    user_limit = user_work,
    source = if (work_cap_enabled) "explicit finite control$counted_work_limit" else
      if ("counted_work_limit" %in% resolved$overrides) "explicit control$counted_work_limit = Inf (opt-out)" else "default: no user work budget",
    effective_private_max_work_ops = work,
    internal_arithmetic_ceiling = chained_public_work_arithmetic_ceiling(),
    meaning = if (work_cap_enabled)
      "the user's finite estimated-work cap reaches every private stage control unchanged; exceeding it refuses the call before the stage runs" else
      paste("no user work budget; the private stage controls carry the exact-count arithmetic ceiling 2^53 - 1 so that counts stay",
            "representable, which is not a resource budget, not a runtime measurement and not a CPU/memory guarantee"))
  chunk <- chained_public_chunk_rows(dims, entries, resolved)
  staging <- chained_public_staging_envelope()
  cov_strat <- resolved$covariance_strategy
  if (!is.null(gmm_weighting) || !is.null(gmm_covariance) || is.null(cov_strat)) {
    w_choice <- resolve_gmm_weighting_choice(
      if (!is.null(gmm_weighting)) gmm_weighting else if (!is.null(resolved$gmm_weighting)) resolved$gmm_weighting$requested else "optimal",
      arg = "gmm_weighting", error_fn = chained_public_abort, error_class = "cdid_chained_public_control_error"
    )
    cov_choice <- resolve_gmm_covariance_choice(
      if (!is.null(gmm_covariance)) gmm_covariance else if (!is.null(resolved$gmm_covariance)) resolved$gmm_covariance else "tikhonov",
      arg = "gmm_covariance", error_fn = chained_public_abort, error_class = "cdid_chained_public_control_error"
    )
    tuning <- if (!is.null(resolved$covariance_tuning)) resolved$covariance_tuning else resolve_gmm_covariance_tuning(NULL)
    cov_strat <- gmm_covariance_strategy(
      method = cov_choice,
      condition_trigger = tuning$condition_trigger,
      ridge_relative = tuning$ridge_relative,
      overrides = tuning$overrides,
      weighting_requested = w_choice$requested,
      weighting_canonical = w_choice$canonical
    )
  }
  lnk <- ipw_link_enclosure_control(max_chunk_rows = chunk, max_retained_entries = entries, max_matrix_entries = entries,
                                    max_work_ops = work, max_output_entries = entries)
  chained_public_need(chained_public_count(workers, "workers") <= prof$workers,
    "The requested multiplier workers exceed the declared `control$workers`.", "cdid_chained_public_control_error")
  mb <- chained_public_multiplier_budget(dims, dims$biters, chunk, entries, bytes, workers = workers)
  ctrl <- ipw_pipeline_control(
    owner_control = ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = entries, max_matrix_entries = entries,
                                      max_work_ops = work, max_output_entries = entries),
    fitted_control = ipw_fitted_view_control(max_chunk_rows = chunk, max_observer_units = dims$n_frame, max_retained_entries = entries,
                                             max_matrix_entries = entries, max_work_ops = work, max_output_entries = entries),
    point_control = ipw_point_view_control(max_chunk_rows = chunk, max_observer_units = dims$n_frame, max_retained_entries = entries,
                                           max_matrix_entries = entries, max_work_ops = work, max_output_entries = entries),
    glm_control = list(epsilon = resolved$glm$epsilon, maxit = as.integer(resolved$glm$maxit), trace = FALSE),
    certificate_control = logit_certificate_control(taylor_degree = resolved$certificate$taylor_degree, max_abs_eta = resolved$certificate$max_abs_eta,
                                                    max_squarings = resolved$certificate$max_squarings, max_chunk_rows = chunk,
                                                    max_matrix_entries = entries, max_work_ops = work),
    problem_control = ipw_nuisance_problem_control(max_chunk_rows = chunk, max_retained_entries = entries, max_matrix_entries = entries,
                                                   max_work_ops = work),
    link_control = lnk,
    covariance_control = ipw_covariance_transfer_control(max_chunk_rows = chunk, max_operator_block_rows = chunk, max_cluster_rows = dims$n_clusters,
                                                         max_matrix_entries = entries, max_work_ops = work, max_output_entries = entries,
                                                         link_control = lnk),
    class_bounds = list(max_links = dims$links_bound, max_classes = dims$links_bound, max_class_dates = dims$n_dates,
                        max_E_entries = dims$E_entries_bound, max_mask_entries = dims$mask_entries_bound),
    gmm_bounds = gmm_bounds_control(max_matrix_entries = entries, max_batch_rows = chunk, max_targets = dims$J_bound,
                                    max_open_clusters = dims$n_clusters, allow_materialize_if = TRUE,
                                    max_materialized_if_entries = min(dims$materialized_if_bound, entries)),
    gmm_accuracy = gmm_accuracy_control(method = "conventional_chol_qr", target_error_unit = "abs",
                                        max_target_error_abs = staging, max_target_error_se = NA_real_,
                                        max_orthogonality_defect = resolved$gmm_guard$max_orthogonality_defect,
                                        min_certified_lambda_rel = resolved$gmm_guard$min_certified_lambda_rel,
                                        max_factor_residual_rel = resolved$gmm_guard$max_factor_residual_rel,
                                        min_rcond_R = resolved$gmm_guard$min_rcond_R, min_rcond_U = resolved$gmm_guard$min_rcond_U,
                                        min_sigma_rel_A = resolved$gmm_guard$min_sigma_rel_A,
                                        max_kkt_residual_rel = resolved$gmm_guard$max_kkt_residual_rel,
                                        max_loading_residual_rel = resolved$gmm_guard$max_loading_residual_rel,
                                        require_matching_scale = identical(cov_strat$weighting_canonical, "optimal"), max_loading_error_abs = staging,
                                        covariance_strategy = cov_strat),
    catalog_bounds = list(max_candidates = dims$candidates_bound, max_rows = dims$links_bound, max_pool_entries = dims$pool_entries_bound,
                          max_scratch = entries, max_work = work),
    multiplier_budget = mb,
    inference_accuracy = NULL)
  attr(ctrl, "covariance_strategy") <- cov_strat
  record <- list(
    profile = prof, profile_status = resolved$profile_status, profile_source = resolved$source, overrides = resolved$overrides,
    # `counted_work_limit` is the transparent public value (Inf when the D78
    # optional cap is disabled); the finite value the private constructors
    # received is recorded distinctly in `work_cap`, never relabelled as the user's cap.
    numeric_workspace_entries = entries, counted_work_limit = user_work, work_cap = work_cap, chunk_rows = chunk,
    chunk_rule = if (is.null(resolved$chunk_rows)) "min(n_frame, 256, floor(workspace_entries / (2048 (links_bound + p + 1))))" else "explicit control$chunk_rows",
    dimensions = dims,
    derived_bounds = data.frame(
      bound = c("owner/fitted/point/problem/link/covariance max_retained/matrix/output entries",
                if (work_cap_enabled) "work limits (max_work_ops, catalog max_work): explicit user counted_work_limit" else
                  "work limits (max_work_ops, catalog max_work): user work cap disabled; exact-count arithmetic ceiling 2^53 - 1",
                "max_observer_units", "catalog max_candidates", "catalog max_rows / class max_links / max_classes", "catalog max_pool_entries",
                "class max_class_dates", "class max_E_entries", "class max_mask_entries", "gmm max_targets / multiplier max_targets",
                "gmm max_materialized_if_entries", "multiplier max_units", "multiplier max_clusters / covariance max_cluster_rows / gmm max_open_clusters",
                "multiplier max_draws", "chunk rows (all chunked producers, IF batch rows, gmm batch rows)", "multiplier max_sign_batch_rows"),
      value = c(entries, work, dims$n_frame, dims$candidates_bound, dims$links_bound, dims$pool_entries_bound, dims$n_dates, dims$E_entries_bound,
                dims$mask_entries_bound, dims$J_bound, min(dims$materialized_if_bound, entries), dims$n_frame, dims$n_clusters, dims$biters, chunk,
                mb$max_sign_batch_rows),
      stringsAsFactors = FALSE),
    glm = resolved$glm, unused_private_certificate_settings = resolved$certificate, gmm_guard = resolved$gmm_guard,
    gmm_weighting = list(requested = cov_strat$weighting_requested, canonical = cov_strat$weighting_canonical, alias_used = cov_strat$alias_used),
    covariance_strategy = cov_strat,
    numerical_policy = list(decision = "D23", mode = "conventional",
      complete_per_fit_certification_requested = FALSE, complete_per_fit_certification_performed = FALSE,
      relative_accuracy_status = "deferred_by_D23", override_recorded = FALSE,
      rule = "stable Cholesky/QR with existing conditioning and factor residual policies; complete IF/inference; residual ratios diagnostic only; no universal per-fit numerical-error guarantee"),
    admission = list(
      numeric_payload_bytes = bytes, numeric_payload_entries = entries, counted_work_limit = user_work,
      work_cap_status = work_cap$status, private_max_work_ops = work,
      scope = paste("bounds the reachable numeric payload the accepted preflights count; the estimated-work cap is optional and,",
        "when set, bounds their conservative work ESTIMATES, which are not measured runtime; it is not process RSS, R heap overhead,",
        "interpreter allocations or CPU time, and no whole-call resource guarantee is made"),
      whole_process_RSS_guarantee = FALSE, measurement_status = "no whole-call resource guarantee or performance claim"),
    execution = list(declared_workers = prof$workers, multiplier_workers = as.double(workers),
      scope = paste("estimation is serial; only multiplier-bootstrap raw-draw products may use worker processes,",
        "whose whole buffers are charged against the same numeric workspace")),
    witness = "not requested/performed on the conventional-computation public path",
    retention = "no observation, moment or record is dropped on refusal; the accepted pipeline's typed bound errors stop the call")
  list(control = ctrl, record = record)
}

# The accepted automatic nuisance witness (identical to the accepted private
# example and tests): the fitted Cholesky metric factor and the one-radius request.
chained_public_automatic_witness <- function(fv) {
  list(F = unname(fv$factors$H_chol), r = logit_cert_automatic_radius_request())
}
