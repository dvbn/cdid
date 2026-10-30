# Private CP02 compact point-link, bounded complete point-IF replay, and
# fitted-candidate-to-certificate root/link/IF binding module.
#
# Scope and boundaries:
#   * Consumes sealed `cdid_ipw_fitted_nuisance_view`, `cdid_ipw_link_view`,
#     `cdid_ipw_nuisance_view`, `cdid_ipw_frame_owner`, and `cdid_ipw_panel_owner`
#     objects without refitting the weighted propensity logit (`stats::glm.fit` is
#     never called here) and without forging a dense context.
#   * Requires verified common panel/frame ownership, identical ordered cohorts
#     `(g, K)`, identical `glm.control`, and actual streamed observer-mask equality
#     (`O = R intersect S`) before reading any link outcome endpoints.
#   * Streams observed-risk rows in canonical A-unit order into at most one bounded
#     temporary observer buffer for the current outcome link, reconstructs the
#     worker's point arrays from exact fitted coefficients, rates, and normalizers,
#     reuses the unchanged `ipw_link_means` and `ipw_link_partials` helpers and
#     existing Cholesky factor `H_chol`, and releases every query slice alias and
#     the temporary observer buffer before return or the next link.
#   * Produces bounded B-by-L complete point-IF chunks and covariance-transfer
#     `candidate_rows` callbacks on the canonical A-frame using all four signed
#     contributions (`treated + control + beta_score + rate`), including risk
#     nonobserver rate terms and outside-risk zeros, with risk-first endpoint
#     filtering and optional bounded one-nuisance chunk caching.
#   * Binds actual fitted point candidates (`beta0`, `H_chol`, `loading_beta`, and
#     bounded `phi` chunks) to the accepted compact raw-root/problem, link, and IF
#     validators under explicit caller `F`, `r`, and resource/certificate controls,
#     keeping point-candidate status, raw-root certification, link/IF enclosure,
#     covariance transfer, stored-array GMM acceptance, and complete statistical
#     inference strictly distinct.

# Preserve exact finite candidate bits, including the sign of zero.
ipw_point_view_identical <- function(x, y) identical(x, y, num.eq = FALSE)

cdid_ipw_point_view_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_ipw_point_view_refusal",
      "cdid_ipw_point_view_error",
      "cdid_ipw_owner_error",
      "cdid_ipw_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

ipw_point_view_require_deps <- function() {
  ns_env <- environment(ipw_point_view_require_deps)
  if (cdid_dependency_check_done("ipw_point_view", ns_env)) return(invisible(TRUE))
  req_fns <- c(
    "cdid_abort", "gmm_sexp_aliases", "gmm_scalar_numeric", "gmm_scalar_logical",
    "gmm_scalar_character", "logit_cert_count_output", "ipw_owner_seal",
    "ipw_owner_seal_ok", "ipw_owner_refuse_decorative_attrs", "ipw_owner_int_bound",
    "ipw_owner_safe_product", "ipw_owner_safe_sum", "ipw_owner_check_matrix_dims",
    "ipw_owner_check_interrupt", "ipw_owner_control", "validate_ipw_owner_control",
    "validate_ipw_panel_owner", "validate_ipw_frame_owner", "validate_ipw_nuisance_view",
    "validate_ipw_link_view", "construct_ipw_link_view", "query_ipw_view_chunk",
    "preflight_ipw_range_query", "preflight_ipw_link_view", "verify_ipw_view_mask_compatibility",
    "validate_ipw_glm_control", "ipw_fitted_view_control", "validate_ipw_fitted_view_control",
    "validate_ipw_fitted_nuisance_view", "verify_ipw_fitted_nuisance_sharing",
    "ipw_pair_rates", "ipw_logit_values", "ipw_clip_floor", "ipw_obs_weights", "ipw_link_means",
    "ipw_link_partials", "validate_logit_certificate_control",
    "validate_ipw_nuisance_problem_control", "validate_ipw_link_enclosure_control",
    "construct_certified_ipw_nuisance_problem_view", "enclose_ipw_link_view",
    "enclose_ipw_link_if_view_chunk", "ipw_validate_nuisance_problem_mode",
    "ipw_validate_link_enclosure_mode", "preflight_ipw_link_if_view_chunk",
    "ipw_candidate_rows_spec", "gmm_require_native_dll", "cdid_gmm_abort",
    "ipw_owner_view_cohort_index", "ipw_owner_pair_complete", "ipw_owner_pair_masks_equal_fact",
    "ipw_owner_validated_view_record", "ipw_owner_validated_view_match", "ipw_owner_read_y_positions"
  )
  for (fn in req_fns) {
    if (!exists(fn, envir = ns_env, mode = "function", inherits = TRUE)) {
      cdid_ipw_point_view_abort(
        paste0("Required private helper `", fn, "` is not loaded in the namespace."),
        "cdid_ipw_point_view_dependency_error"
      )
    }
  }
  cdid_dependency_check_record("ipw_point_view", ns_env)
}

# ---------------------------------------------------------------------------
# Explicit private resource controls (NO internal defaults)
# ---------------------------------------------------------------------------

ipw_point_view_control <- function(
  max_chunk_rows,
  max_observer_units,
  max_retained_entries,
  max_matrix_entries,
  max_work_ops,
  max_output_entries
) {
  ipw_point_view_require_deps()
  if (missing(max_chunk_rows) || missing(max_observer_units) ||
      missing(max_retained_entries) || missing(max_matrix_entries) ||
      missing(max_work_ops) || missing(max_output_entries)) {
    cdid_ipw_point_view_abort(
      "All fields of `ipw_point_view_control` (`max_chunk_rows`, `max_observer_units`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) must be supplied explicitly by the caller (no internal defaults).",
      c("cdid_ipw_point_view_control_error", "cdid_ipw_owner_control_error", "cdid_ipw_owner_bound_error")
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
    class = "cdid_ipw_point_view_control"
  )
}

validate_ipw_point_view_control <- function(control) {
  ipw_point_view_require_deps()
  ipw_owner_refuse_decorative_attrs(control, "control")
  req <- c(
    "max_chunk_rows", "max_observer_units", "max_retained_entries",
    "max_matrix_entries", "max_work_ops", "max_output_entries"
  )
  if (missing(control) || is.null(control) || !is.list(control) ||
      length(control) != length(req) || is.null(names(control)) ||
      anyDuplicated(names(control)) || !setequal(names(control), req)) {
    cdid_ipw_point_view_abort(
      "Explicit `control` list (`max_chunk_rows`, `max_observer_units`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) is required; no internal defaults exist.",
      c("cdid_ipw_point_view_control_error", "cdid_ipw_owner_control_error", "cdid_ipw_owner_bound_error")
    )
  }
  ipw_point_view_control(
    max_chunk_rows = control$max_chunk_rows,
    max_observer_units = control$max_observer_units,
    max_retained_entries = control$max_retained_entries,
    max_matrix_entries = control$max_matrix_entries,
    max_work_ops = control$max_work_ops,
    max_output_entries = control$max_output_entries
  )
}

ipw_point_view_to_owner_control <- function(ctrl) {
  ipw_owner_control(
    max_chunk_rows = ctrl$max_chunk_rows,
    max_retained_entries = ctrl$max_retained_entries,
    max_matrix_entries = ctrl$max_matrix_entries,
    max_work_ops = ctrl$max_work_ops,
    max_output_entries = ctrl$max_output_entries
  )
}

ipw_point_view_to_fitted_control <- function(ctrl) {
  ipw_fitted_view_control(
    max_chunk_rows = ctrl$max_chunk_rows,
    max_observer_units = ctrl$max_observer_units,
    max_retained_entries = ctrl$max_retained_entries,
    max_matrix_entries = ctrl$max_matrix_entries,
    max_work_ops = ctrl$max_work_ops,
    max_output_entries = ctrl$max_output_entries
  )
}

ipw_point_view_char_bytes <- function(chr_list, control) {
  total_bytes <- 0
  max_bytes <- 0
  work <- 128
  for (x in chr_list) {
    if (!is.null(x) && is.character(x)) {
      for (j in seq_along(x)) {
        nb <- as.double(nchar(x[j], type = "bytes"))
        if (is.na(nb)) nb <- 0
        work <- ipw_owner_safe_sum(
          c(work, nb, 8),
          control$max_work_ops,
          "point-view character metadata scan work"
        )
        total_bytes <- ipw_owner_safe_sum(
          c(total_bytes, nb),
          label = "point-view character metadata byte count"
        )
        max_bytes <- max(max_bytes, nb)
        ipw_owner_safe_sum(
          c(ceiling(total_bytes / 2), 32),
          control$max_matrix_entries,
          "point-view character metadata scratch slots"
        )
      }
    }
  }
  slots <- ipw_owner_safe_sum(
    c(ceiling(total_bytes / 2), 32),
    label = "point-view character metadata charged slots"
  )
  list(
    source_bytes = total_bytes,
    max_string_bytes = max_bytes,
    charged_slots = slots,
    charged_bytes = ipw_owner_safe_product(8, slots),
    work_ops_bound = work
  )
}

ipw_point_view_id_bytes_bound <- function(panel_owner, b_len, n_chunks, control) {
  if (!ipw_point_view_identical(panel_owner$id_type, "character")) {
    return(list(
      id_type = "numeric",
      total_id_bytes = 0,
      max_id_bytes = 0,
      chunk_id_bytes_bound = 0,
      chunk_id_slots = 0,
      stream_id_work_bound = 0
    ))
  }
  bm <- panel_owner$id_byte_metadata
  tot_b <- as.double(bm$total)
  max_b <- as.double(bm$max)
  chunk_b <- min(tot_b, ipw_owner_safe_product(b_len, max_b, label = "chunk character ID byte bound"))
  chunk_slots <- ipw_owner_safe_sum(
    c(ceiling(chunk_b / 2), 32),
    control$max_matrix_entries,
    "chunk character ID key scratch slots"
  )
  stream_work <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(8, tot_b, control$max_work_ops, "whole-traversal character ID key byte work"),
      ipw_owner_safe_product(16, n_chunks, control$max_work_ops, "chunk character ID slice dispatch work")
    ),
    control$max_work_ops,
    "total character ID key work bound"
  )
  list(
    id_type = "character",
    total_id_bytes = tot_b,
    max_id_bytes = max_b,
    chunk_id_bytes_bound = chunk_b,
    chunk_id_slots = chunk_slots,
    stream_id_work_bound = stream_work
  )
}

# Budget compact, caller-reachable payload BEFORE validators/content comparisons.
# Environments are borrowed sealed references (or caller-owned trace state), never
# traversed. nchar(type="bytes") fetches stored CHARSXP byte lengths; only after
# this length gate may content scans, hashing, or equality comparisons run.
# Every term is a nonnegative integer, so each running sum is exact while it
# stays below 2^53; the checked helper is called only where a sum would leave
# its bound (or that exact range), and then raises the identical refusal with
# the same label, totals and diagnostics at the same node.
ipw_point_view_payload_budget <- function(objects, control) {
  entries <- 0; bytes <- 0; max_bytes <- 0; work <- 128
  limit <- 2^53 - 1
  max_me <- control$max_matrix_entries
  max_wo <- control$max_work_ops
  cap_me <- if (is.null(max_me)) limit else min(limit, max_me)
  cap_wo <- if (is.null(max_wo)) limit else min(limit, max_wo)
  visit <- function(x, depth = 0L, scan_attrs = TRUE) {
    if (is.null(x) || is.environment(x)) return(invisible(NULL))
    if (depth > 12L || is.function(x) || is.language(x) || isS4(x)) {
      cdid_ipw_point_view_abort("Unsupported compact payload.", "cdid_ipw_point_view_input_error")
    }
    k <- length(x)
    total <- entries + k + 8
    entries <<- if (total > cap_me) ipw_owner_safe_sum(c(entries, k, 8), max_me,
                                                       "compact payload length gate") else total
    total <- work + 64 * (k + 8)
    work <<- if (total > cap_wo) ipw_owner_safe_sum(c(work, 64 * (k + 8)), max_wo,
                                                    "compact payload traversal gate") else total
    ax <- attributes(x)
    if (scan_attrs && !is.null(ax)) visit(ax, depth + 1L, FALSE)
    if (is.list(x)) {
      for (j in seq_len(k)) visit(x[[j]], depth + 1L)
    } else if (is.character(x)) {
      for (j in seq_len(k)) {
        nb <- nchar(x[[j]], type = "bytes", allowNA = TRUE)
        if (is.na(nb)) nb <- 0
        total <- bytes + nb
        bytes <<- if (total > limit) ipw_owner_safe_sum(c(bytes, nb), label = "compact character bytes") else total
        max_bytes <<- max(max_bytes, nb)
        total <- entries + ceiling(nb / 2)
        entries <<- if (total > cap_me) ipw_owner_safe_sum(c(entries, ceiling(nb / 2)), max_me,
                                                           "compact character length gate") else total
        total <- work + 16 * nb
        work <<- if (total > cap_wo) ipw_owner_safe_sum(c(work, 16 * nb), max_wo,
                                                        "compact character comparison gate") else total
      }
    }
    invisible(NULL)
  }
  visit(objects)
  list(entries = entries, work = work, bytes = bytes, max_bytes = max_bytes)
}

# Constant-time seal/shape gate; deeper validators follow the payload budget.
ipw_point_view_seal_gate <- function(x, binding) {
  if (!is.list(x) || !ipw_owner_seal_ok(x, binding)) {
    cdid_ipw_point_view_abort("Unsealed or foreign compact input.",
      c("cdid_ipw_point_view_integrity_error", "cdid_ipw_owner_integrity_error"))
  }
  invisible(NULL)
}

ipw_point_view_owner_entries <- function(frame_owner, panel_owner) {
  ipw_owner_safe_sum(c(panel_owner$accounting$actual_retained_entries,
    frame_owner$accounting$retained_matrix_entries,
    frame_owner$accounting$retained_vector_entries), label = "borrowed owner storage once")
}

# Verify that a compact point link view contains zero observer-length or unit-length arrays.
ipw_point_view_verify_compact_storage <- function(obj, p_max = NULL, m_max = NULL) {
  forbidden_fields <- c(
    "risk_a_idx", "risk_gi", "risk_w", "risk_dy", "dy_c", "dy_t",
    "Z_use", "D_use", "fit_weights", "eta_use", "odds_use", "prob_use",
    "is_t", "is_c", "a_obs", "b_obs", "w_hess", "score_factor",
    "y_left_obs", "y_right_obs", "Delta_iv_obs", "observed_a_index",
    "fitted.values", "linear.predictors", "residuals", "effects",
    "weights", "y", "qr", "phi", "phi_matrix"
  )
  dim_cap <- if (!is.null(p_max) && !is.null(m_max)) max(as.integer(p_max), as.integer(m_max)) else NULL
  check_node <- function(x, path) {
    if (is.environment(x) || is.null(x)) return(invisible(TRUE))
    if (is.list(x)) {
      nms <- names(x)
      if (!is.null(nms)) {
        # First offending field in `names(x)` order: the primitive `match` gives
        # exactly the element `intersect(nms, forbidden_fields)[1L]` named,
        # without the set helper's coercion, uniqueness pass and allocation.
        m_bad <- match(nms, forbidden_fields, 0L)
        if (any(m_bad > 0L)) {
          cdid_ipw_point_view_abort(
            paste0("Forbidden non-compact observer/model field `", nms[which(m_bad > 0L)[1L]], "` retained at `", path, "`."),
            c("cdid_ipw_point_view_storage_error", "cdid_ipw_owner_integrity_error")
          )
        }
      }
      for (i in seq_along(x)) {
        nm_i <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
        check_node(x[[i]], paste0(path, "$", nm_i))
      }
      return(invisible(TRUE))
    }
    if (!is.null(dim_cap) && is.atomic(x) && (is.numeric(x) || is.integer(x) || is.logical(x))) {
      d <- dim(x)
      if (!is.null(d)) {
        if (any(d > dim_cap)) {
          cdid_ipw_point_view_abort(
            paste0("Non-compact matrix dimension (", paste(d, collapse = "x"), ") retained at `", path, "`."),
            c("cdid_ipw_point_view_storage_error", "cdid_ipw_owner_integrity_error")
          )
        }
      } else if (length(x) > dim_cap) {
        cdid_ipw_point_view_abort(
          paste0("Non-compact vector length (", length(x), ") retained at `", path, "`."),
          c("cdid_ipw_point_view_storage_error", "cdid_ipw_owner_integrity_error")
        )
      }
    }
    invisible(TRUE)
  }
  check_node(obj, "point_link_view")
}

# ---------------------------------------------------------------------------
# Section 1: Compact Point Link (Preflight, Fit, Set Fit, Validation, RDS)
# ---------------------------------------------------------------------------

preflight_ipw_point_link_view <- function(
  fitted_view,
  link_view,
  control,
  glm_control = NULL,
  n_planned_links = 1L
) {
  ipw_point_view_require_deps()
  if (missing(fitted_view) || missing(link_view) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`fitted_view`, `link_view`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  ipw_owner_refuse_decorative_attrs(fitted_view, "fitted_view")
  ipw_owner_refuse_decorative_attrs(link_view, "link_view")
  ipw_owner_refuse_decorative_attrs(control, "control")
  if (!is.null(glm_control)) {
    ipw_owner_refuse_decorative_attrs(glm_control, "glm_control")
  }

  ctrl <- validate_ipw_point_view_control(control)
  ipw_point_view_seal_gate(fitted_view, "fitted_nuisance_view")
  ipw_point_view_seal_gate(link_view, "link_view")
  payload <- ipw_point_view_payload_budget(list(fitted_view, link_view, glm_control), ctrl)
  validate_ipw_fitted_nuisance_view(fitted_view)
  validate_ipw_link_view(link_view)
  n_links <- ipw_owner_int_bound(n_planned_links, "n_planned_links", min_val = 1)

  if (!is.null(glm_control)) {
    glm_ctrl <- validate_ipw_glm_control(glm_control)
    if (!ipw_point_view_identical(glm_ctrl, fitted_view$glm_control)) {
      cdid_ipw_point_view_abort(
        "Caller `glm_control` differs from the fitted nuisance view's `glm_control`; point-link reuse is refused.",
        c("cdid_ipw_control_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  } else {
    glm_ctrl <- validate_ipw_glm_control(fitted_view$glm_control)
  }

  if (!ipw_point_view_identical(fitted_view$frame_ref, link_view$frame_ref) ||
      !ipw_point_view_identical(fitted_view$panel_ref, link_view$panel_ref)) {
    cdid_ipw_point_view_abort(
      "Link view belongs to a different frame or panel owner than the fitted nuisance view; reuse is refused.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_context_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  if (!ipw_point_view_identical(fitted_view$cohorts$g, link_view$cohorts$g) ||
      !ipw_point_view_identical(fitted_view$cohorts$K, link_view$cohorts$K) ||
      !ipw_point_view_identical(fitted_view$cohorts$cohorts, link_view$cohorts$cohorts) ||
      !ipw_point_view_identical(fitted_view$cohorts$cohort_labels, link_view$cohorts$cohort_labels)) {
    cdid_ipw_point_view_abort(
      "Link view `g` or ordered control cohorts `K` differ from the fitted nuisance view; reuse is refused.",
      c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  frame_owner <- get("frame_owner", envir = fitted_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = fitted_view$panel_ref, inherits = FALSE)
  nv_obj <- get("nuisance_view", envir = fitted_view$nuisance_view_ref, inherits = FALSE)
  validate_ipw_frame_owner(frame_owner, panel_owner = panel_owner)
  validate_ipw_link_view(link_view, frame_owner = frame_owner)
  if (!is.null(link_view$nuisance_ref)) {
    if (!ipw_point_view_identical(link_view$nuisance_ref, fitted_view$nuisance_view_ref)) {
      cdid_ipw_point_view_abort(
        "Link view is bound to a different `nuisance_view` snapshot than the fitted nuisance view; reuse is refused.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_shared_mask_error", "cdid_ipw_shared_incompatible_error")
      )
    }
    validate_ipw_link_view(link_view, frame_owner = frame_owner, nuisance_view = nv_obj)
  }

  na <- ipw_owner_int_bound(frame_owner$n_a_units, "n_a_units", min_val = 1)
  n_obs <- ipw_owner_int_bound(link_view$support$n_observed_risk_units, "n_observed_risk_units", min_val = 0)
  p <- ipw_owner_int_bound(frame_owner$n_covariates, "n_covariates", min_val = 1)
  m <- ipw_owner_int_bound(fitted_view$cohorts$m_cohorts, "m_cohorts", min_val = 2)

  if (!ipw_point_view_identical(as.integer(n_obs), as.integer(fitted_view$counts$n_observed_risk_units)) ||
      !ipw_point_view_identical(as.integer(link_view$support$by_cohort_pair_count), as.integer(unname(fitted_view$counts$by_cohort_pair)))) {
    cdid_ipw_point_view_abort(
      "Link view observed risk-set support counts differ from the fitted nuisance view.",
      c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  max_ch <- ctrl$max_chunk_rows
  max_obs <- ctrl$max_observer_units
  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  if (n_obs > max_obs) {
    cdid_ipw_point_view_abort(
      paste0("Observed risk-set pair unit count (", n_obs, ") exceeds `control$max_observer_units` (", max_obs, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error"),
      diagnostics = list(n_observed_risk_units = n_obs, max_observer_units = max_obs)
    )
  }
  if (n_obs < p || any(link_view$support$by_cohort_pair_count < 1L) ||
      any(!is.finite(link_view$support$by_cohort_pair_mass) | link_view$support$by_cohort_pair_mass <= 0)) {
    cdid_ipw_point_view_abort(
      "Link view lacks positive observed risk-set support across all declared cohorts.",
      c("cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }

  B_eff <- min(na, max_ch)
  n_chunks <- ceiling(na / B_eff)

  # Possible failure diagnostic output (observer dy/weights/probabilities + partials)
  failure_output_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(5, max(1, n_obs)),
      ipw_owner_safe_product(p, m),
      ipw_owner_safe_product(2, p),
      ipw_owner_safe_product(4, m),
      256
    ),
    min(max_ret, max_out),
    "possible failed point-link diagnostic output"
  )

  # Constant header allowance before any character metadata scan
  ipw_owner_safe_sum(
    c(ipw_owner_safe_product(8, p + m + 1), 32),
    min(max_ret, max_out, max_me),
    "point-link metadata header before character scans"
  )
  ipw_owner_safe_sum(
    c(ipw_owner_safe_product(16, p + m + 1), 128),
    max_wo,
    "point-link header work before character scans"
  )

  # 1. Character metadata and ID key byte charges
  char_meta <- ipw_point_view_char_bytes(list(
    colnames(frame_owner$Z_a),
    fitted_view$cohorts$cohort_labels,
    fitted_view$nuisance_id,
    link_view$link_id,
    link_view$functional_key
  ), ctrl)
  id_meta <- ipw_point_view_id_bytes_bound(panel_owner, B_eff, n_chunks, ctrl)

  # 2. Compact retained state per point link:
  p2 <- ipw_owner_check_matrix_dims(p, p, max_matrix_entries = max_ret, label = "point-link p x p factor matrix")
  pm <- ipw_owner_check_matrix_dims(p, m, max_matrix_entries = max_ret, label = "point-link p x m score-rate Jacobian")
  retained_matrix_entries_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(3, p2, max_ret, "H, H_chol, and qr_R matrices"),
      pm,
      link_view$preflight$retained_matrix_entries
    ),
    max_bound = max_ret,
    label = "retained point-link matrix entries"
  )
  retained_vector_entries_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(5, p, max_ret, "beta coefficients, qr_pivot, G_beta, and loading_beta"),
      ipw_owner_safe_product(18, m, max_ret, "G_r, loading_rate, rate_coef_s0, rate_coef_s1, counts, masses, and disclosures"),
      char_meta$charged_slots,
      link_view$preflight$retained_vector_entries,
      2048, # D27 fitted-view diagnostics carried by reference into each point link
      8192
    ),
    max_bound = max_ret,
    label = "retained point-link vector and metadata entries"
  )
  retained_entries_one <- ipw_owner_safe_sum(
    c(retained_matrix_entries_one, retained_vector_entries_one),
    max_bound = max_ret,
    label = "single point-link retained entries"
  )
  total_retained_entries <- ipw_owner_safe_product(
    n_links, retained_entries_one, max_entries = max_ret,
    label = "total retained point-link entries across planned links"
  )
  if (total_retained_entries > max_out) {
    cdid_ipw_point_view_abort(
      paste0("Total retained point-link entries (", total_retained_entries, ") exceed `control$max_output_entries` (", max_out, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error"),
      diagnostics = list(requested = total_retained_entries, max_output_entries = max_out)
    )
  }

  # 3. Single temporary observer buffer: risk_a_idx (n_obs), risk_gi (n_obs), risk_w (n_obs), risk_dy (n_obs), Z_use (n_obs x p)
  obs_Z_entries <- ipw_owner_check_matrix_dims(n_obs, p, max_matrix_entries = max_me, label = "point-link observer buffer Z_use matrix")
  observer_buffer_entries <- ipw_owner_safe_sum(
    c(
      obs_Z_entries,
      ipw_owner_safe_product(4, max(1, n_obs), max_me, "observer buffer index/cohort/weight/dy vectors"),
      32
    ),
    max_bound = max_me,
    label = "temporary point-link observer buffer entries"
  )

  # 4. Simultaneous bounded owner-reader chunk buffer (at most 1 chunk active at a time during streaming)
  owner_ctrl <- ipw_point_view_to_owner_control(ctrl)
  chunk_pre <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p,
    m_cohorts = m,
    control = owner_ctrl,
    read_outcomes = TRUE,
    include_design = TRUE,
    enclose_differences = FALSE
  )
  reader_chunk_concurrent_entries <- chunk_pre$concurrent_entries

  # 5. Point helper workspace (ipw_logit_values, ipw_obs_weights, ipw_link_means, ipw_link_partials, Cholesky solves)
  point_calc_matrix_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(4, obs_Z_entries, max_me, "control Z_c and weighted Z_c copies in ipw_link_partials"),
      ipw_owner_safe_product(4, p2, max_me, "H_chol and transpose triangular solve workspace"),
      ipw_owner_safe_product(2, pm, max_me, "A_beta_r crossprod workspace")
    ),
    max_bound = max_me,
    label = "point-link helper matrix workspace entries"
  )
  point_calc_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(64, max(1, n_obs), max_me, "model vectors, conversions, helper subsets, expression and diagnostic copies"),
      ipw_owner_safe_product(8, p, max_me, "G_beta and triangular solve working vectors"),
      ipw_owner_safe_product(12, m, max_me, "G_r, loading_rate, and rate_coef working vectors"),
      char_meta$charged_slots,
      id_meta$chunk_id_slots,
      128
    ),
    max_bound = max_me,
    label = "point-link helper vector workspace entries"
  )
  point_calc_entries <- ipw_owner_safe_sum(
    c(point_calc_matrix_entries, point_calc_vector_entries),
    max_bound = max_me,
    label = "total point-link helper workspace entries"
  )

  # 6. Simultaneous peak concurrent entries (including borrowed fitted_view and link_view retained envelopes)
  referenced_entries <- ipw_owner_safe_sum(
    c(payload$entries, ipw_point_view_owner_entries(frame_owner, panel_owner)),
    max_bound = max_me,
    label = "referenced fitted nuisance and link view retained entries"
  )
  mask_sharing_concurrent_entries <- ipw_owner_safe_sum(
    c(total_retained_entries, referenced_entries, 64 * B_eff, 32 * m, 256),
    max_bound = max_me,
    label = "simultaneous mask-sharing check entries"
  )
  streaming_phase_entries <- ipw_owner_safe_sum(
    c(
      total_retained_entries,
      referenced_entries,
      observer_buffer_entries,
      reader_chunk_concurrent_entries,
      observer_buffer_entries, # environment subassignment may copy a full buffer
      4 * B_eff * p + 32 * B_eff, # conversions, RHS slices, order/difference temporaries
      char_meta$charged_slots,
      id_meta$chunk_id_slots
    ),
    max_bound = max_me,
    label = "simultaneous point-link streaming-phase entries"
  )
  point_calc_phase_entries <- ipw_owner_safe_sum(
    c(
      total_retained_entries,
      referenced_entries,
      observer_buffer_entries,
      point_calc_entries,
      failure_output_entries
    ),
    max_bound = max_me,
    label = "simultaneous point-link helper-phase entries"
  )
  peak_concurrent_entries <- max(mask_sharing_concurrent_entries, streaming_phase_entries, point_calc_phase_entries)
  peak_concurrent_matrix_entries <- max(
    retained_matrix_entries_one * n_links + obs_Z_entries + chunk_pre$output_matrix_entries,
    retained_matrix_entries_one * n_links + obs_Z_entries + point_calc_matrix_entries
  )
  if (peak_concurrent_matrix_entries > max_me) {
    cdid_ipw_point_view_abort(
      paste0("Peak concurrent matrix entries (", peak_concurrent_matrix_entries, ") exceed `control$max_matrix_entries` (", max_me, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error")
    )
  }

  # 7. Whole-traversal work operations bound (mask sharing check + chunk streaming + point helpers + metadata)
  mask_sharing_work_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(128 * (1 + ceiling(log2(m + 1))), na, max_wo, "streamed observer-mask cohort check work"),
      ipw_owner_safe_product(128 * (1 + ceiling(log2(na + 1))), n_chunks, max_wo, "streamed observer-mask posting search work")
    ),
    max_bound = max_wo,
    label = "single-link observer-mask sharing verification work"
  )
  stream_work_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, chunk_pre$work_ops_bound, max_wo, "bounded link-view chunk query work"),
      ipw_owner_safe_product(6, max(1, n_obs) * p, max_wo, "observer buffer copy and order verification work"),
      ipw_owner_safe_product(4, max(1, n_obs), max_wo, "observer endpoint subtraction and finiteness check work")
    ),
    max_bound = max_wo,
    label = "single-link stream work operations"
  )
  point_calc_work_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(16, max(1, n_obs) * p, max_wo, "ipw_logit_values and G_beta work"),
      ipw_owner_safe_product(16, max(1, n_obs) * m, max_wo, "ipw_pair_rates and G_r work"),
      ipw_owner_safe_product(4, p * p + p * m, max_wo, "Cholesky triangular solves and A_beta_r crossprod work"),
      256
    ),
    max_bound = max_wo,
    label = "single-link point helper work operations"
  )
  work_ops_one <- ipw_owner_safe_sum(
    c(
      mask_sharing_work_one,
      stream_work_one,
      point_calc_work_one,
      char_meta$work_ops_bound,
      id_meta$stream_id_work_bound,
      8 * payload$work # repeated sealed validation and compact result construction
    ),
    max_bound = max_wo,
    label = "single-link total work operations"
  )
  total_work_ops_bound <- ipw_owner_safe_product(
    n_links, work_ops_one, max_entries = max_wo,
    label = "total planned point-link work operations"
  )

  list(
    n_a_units = na,
    n_observed_risk_units = n_obs,
    p_covariates = p,
    m_cohorts = m,
    n_planned_links = n_links,
    failure_output_entries = failure_output_entries,
    input_payload_entries = payload$entries,
    validation_work_bound = 8 * payload$work,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_link = logit_cert_count_output(n_chunks),
    max_simultaneous_query_chunks = 1L,
    max_simultaneous_observer_buffers = 1L,
    character_metadata = char_meta,
    id_byte_metadata = id_meta,
    observer_buffer_entries = observer_buffer_entries,
    reader_chunk_concurrent_entries = reader_chunk_concurrent_entries,
    point_helper_workspace_entries = point_calc_entries,
    referenced_entries = referenced_entries,
    borrowed_owner_entries = ipw_point_view_owner_entries(frame_owner, panel_owner),
    streaming_phase_entries = streaming_phase_entries,
    point_calc_phase_entries = point_calc_phase_entries,
    peak_concurrent_entries = peak_concurrent_entries,
    peak_concurrent_matrix_entries = peak_concurrent_matrix_entries,
    peak_concurrent_bytes_8double = 8 * peak_concurrent_entries,
    retained_matrix_entries_per_link = retained_matrix_entries_one,
    retained_vector_entries_per_link = retained_vector_entries_one,
    retained_entries_per_link = retained_entries_one,
    total_retained_entries = total_retained_entries,
    total_retained_bytes_8double = 8 * total_retained_entries,
    work_ops_per_link_bound = work_ops_one,
    total_work_ops_bound = total_work_ops_bound,
    expected_endpoints_read_per_link = logit_cert_count_output(2 * n_obs),
    propensity_refit_performed = FALSE,
    peak_rss_promised = FALSE,
    memory_disclosure = paste0(
      "Deterministic preflight charges simultaneous bounded owner-reader chunk buffers, ",
      "the single temporary observer buffer (n_obs x p design plus index/cohort/weight/dy vectors), ",
      "reconstructed point-array and Cholesky solve workspaces, possible failure diagnostic outputs, ",
      "character/ID key byte slots, and compact retained state before any scan, query, or callback. ",
      "Opaque C-runtime, BLAS/LAPACK, and OS process RSS memory remain outside R-level object accounting."
    )
  )
}

fit_ipw_point_link_view <- function(
  fitted_view,
  link_view,
  control,
  glm_control = NULL,
  link_id = NULL,
  functional_key = NULL,
  callback = NULL,
  trace_env = NULL,
  failure_trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(fitted_view) || missing(link_view) || missing(control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`fitted_view`, `link_view`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(callback) && !is.function(callback)) {
    cdid_ipw_point_view_abort("`callback` must be NULL or a function.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(failure_trace_env) && !is.environment(failure_trace_env)) {
    cdid_ipw_point_view_abort("`failure_trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  ctrl <- validate_ipw_point_view_control(control)
  override_budget <- ipw_point_view_payload_budget(list(link_id, functional_key), ctrl)
  if (!is.null(link_id)) {
    ipw_owner_refuse_decorative_attrs(link_id, "link_id")
    if (!gmm_scalar_character(link_id)) {
      cdid_ipw_point_view_abort("`link_id` must be NULL or a non-empty character scalar.", "cdid_ipw_point_view_input_error")
    }
  }
  if (!is.null(functional_key)) {
    ipw_owner_refuse_decorative_attrs(functional_key, "functional_key")
    if (!gmm_scalar_character(functional_key)) {
      cdid_ipw_point_view_abort("`functional_key` must be NULL or a non-empty character scalar.", "cdid_ipw_point_view_input_error")
    }
  }

  # 1. Preflight strictly BEFORE any scan, query, buffer allocation, or callback
  preflight <- preflight_ipw_point_link_view(
    fitted_view = fitted_view,
    link_view = link_view,
    control = control,
    glm_control = glm_control,
    n_planned_links = 1L
  )

  ipw_owner_safe_sum(c(preflight$peak_concurrent_entries, override_budget$entries),
                     ctrl$max_matrix_entries, "point-link override metadata concurrent payload")
  ipw_owner_safe_sum(c(preflight$total_retained_entries, override_budget$entries),
                     min(ctrl$max_retained_entries, ctrl$max_output_entries), "point-link override metadata output")
  ipw_owner_safe_sum(c(preflight$total_work_ops_bound, override_budget$work),
                     ctrl$max_work_ops, "point-link override metadata work")
  preflight$peak_concurrent_entries <- preflight$peak_concurrent_entries + override_budget$entries
  preflight$peak_concurrent_bytes_8double <- 8 * preflight$peak_concurrent_entries
  preflight$retained_entries_per_link <- preflight$retained_entries_per_link + override_budget$entries
  preflight$total_retained_entries <- preflight$total_retained_entries + override_budget$entries
  preflight$total_retained_bytes_8double <- 8 * preflight$total_retained_entries
  preflight$total_work_ops_bound <- preflight$total_work_ops_bound + override_budget$work

  if (!is.null(link_id) && !is.null(link_view$link_id) && !ipw_point_view_identical(link_id, link_view$link_id)) {
    cdid_ipw_point_view_abort(
      "Explicit `link_id` disagrees with `link_view$link_id`.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(functional_key) && !is.null(link_view$functional_key) && !ipw_point_view_identical(functional_key, link_view$functional_key)) {
    cdid_ipw_point_view_abort(
      "Explicit `functional_key` disagrees with `link_view$functional_key`.",
      "cdid_ipw_point_view_input_error"
    )
  }
  eff_link_id <- if (!is.null(link_id)) link_id else link_view$link_id
  eff_func_key <- if (!is.null(functional_key)) functional_key else link_view$functional_key

  ctrl <- validate_ipw_point_view_control(control)
  owner_ctrl <- ipw_point_view_to_owner_control(ctrl)
  fv_ctrl <- ipw_point_view_to_fitted_control(ctrl)
  glm_ctrl <- if (is.null(glm_control)) fitted_view$glm_control else validate_ipw_glm_control(glm_control)

  # 2. Verify actual risk-set observer-mask equality BEFORE reading any endpoints
  sharing_check <- verify_ipw_fitted_nuisance_sharing(
    fitted_view = fitted_view,
    target_view = link_view,
    glm_control = glm_ctrl,
    control = fv_ctrl
  )
  if (!isTRUE(sharing_check$compatible) || !ipw_point_view_identical(as.integer(sharing_check$endpoints_read), 0L)) {
    cdid_ipw_point_view_abort(
      "Actual observer-mask sharing check failed or read outcome endpoints.",
      c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error")
    )
  }

  frame_owner <- get("frame_owner", envir = fitted_view$frame_ref, inherits = FALSE)
  na <- as.integer(frame_owner$n_a_units)
  n_obs <- as.integer(link_view$support$n_observed_risk_units)
  p <- as.integer(frame_owner$n_covariates)
  m <- as.integer(fitted_view$cohorts$m_cohorts)
  step <- as.integer(preflight$effective_chunk_rows)
  # The preflight and sharing check above fully validated this link view and
  # its owners; each chunk query re-checks their identity through this record.
  query_record <- ipw_owner_validated_view_record(link_view, frame_owner,
    get("panel_owner", envir = link_view$panel_ref, inherits = FALSE), step = step)

  obs_buf_env <- new.env(parent = emptyenv())
  diag_env <- new.env(parent = emptyenv())
  obs_buf_env$active <- FALSE
  obs_buf_env$active_query_chunks <- 0L
  obs_buf_env$max_simultaneous_query_chunks <- 0L
  obs_buf_env$released_query_chunks <- 0L
  obs_buf_env$observer_released <- FALSE

  cleanup_temporary_buffers <- function() {
    if (exists("y_trace", envir = parent.env(environment()), inherits = FALSE)) ipw_owner_reset_trace(y_trace)
    for (nm in c("risk_a_idx", "risk_gi", "risk_w", "risk_dy", "Z_use",
                 "r_use", "fit_weights", "lv", "ow", "st", "hs", "means", "part")) {
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
    trace_env$mask_sharing_verified_before_outcomes <- TRUE
    trace_env$queries_started <- 0L
    trace_env$chunks_released <- 0L
    trace_env$max_simultaneous_query_chunks <- 0L
    trace_env$active_observer_buffers <- 0L
    trace_env$max_simultaneous_observer_buffers <- 0L
    trace_env$observer_released <- FALSE
    trace_env$endpoints_read <- 0L
    trace_env$propensity_solver_calls <- 0L
  }

  # 3. Allocate the single temporary observer buffer for the current outcome link
  obs_buf_env$active <- TRUE
  obs_buf_env$risk_a_idx <- integer(n_obs)
  obs_buf_env$risk_gi <- integer(n_obs)
  obs_buf_env$risk_w <- numeric(n_obs)
  obs_buf_env$risk_dy <- numeric(n_obs)
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

  # 4. Stream bounded link-view chunks in canonical A-unit order, reading only
  #    risk-filtered endpoints and releasing each query chunk before the next.
  write_pos <- 0L
  last_a_idx <- 0L
  chunks_queried <- 0L
  endpoints_read_total <- 0L
  stream_ops <- 0
  y_trace <- new.env(parent = emptyenv())

  cur_s <- 1L
  while (cur_s <= na) {
    ipw_owner_check_interrupt()
    cur_e <- as.integer(min(as.double(na), as.double(cur_s) + step - 1))

    ch <- query_ipw_view_chunk(
      view = link_view,
      a_start = cur_s,
      a_end = cur_e,
      control = owner_ctrl,
      read_outcomes = TRUE,
      include_design = TRUE,
      enclose_differences = FALSE,
      trace_env = y_trace,
      validated = query_record
    )
    chunks_queried <- chunks_queried + 1L
    obs_buf_env$active_query_chunks <- obs_buf_env$active_query_chunks + 1L
    obs_buf_env$max_simultaneous_query_chunks <- max(
      obs_buf_env$max_simultaneous_query_chunks,
      obs_buf_env$active_query_chunks
    )

    k_ch <- as.integer(ch$counts$observed_risk_row_count)
    ep_ch <- as.integer(ch$read_diagnostics$n_endpoints_read)
    endpoints_read_total <- endpoints_read_total + ep_ch
    if (!is.null(trace_env)) {
      trace_env$queries_started <- chunks_queried
      trace_env$max_simultaneous_query_chunks <- obs_buf_env$max_simultaneous_query_chunks
      trace_env$endpoints_read <- endpoints_read_total
    }

    if (!isTRUE(ch$read_diagnostics$risk_filtered_before_outcome_read) ||
        isTRUE(ch$read_diagnostics$off_risk_endpoints_read) ||
        isTRUE(ch$read_diagnostics$nonobserver_endpoints_read) ||
        !ipw_point_view_identical(ep_ch, as.integer(2L * k_ch)) ||
        !ipw_point_view_identical(length(y_trace$pos_s_index), k_ch) ||
        !ipw_point_view_identical(length(y_trace$pos_t_index), k_ch)) {
      cdid_ipw_point_view_abort(
        "Link view query violated risk-first endpoint read filtering.",
        c("cdid_ipw_point_view_outcome_leak_error", "cdid_ipw_owner_integrity_error")
      )
    }

    ipw_owner_reset_trace(y_trace)
    if (k_ch > 0L) {
      if (write_pos + k_ch > n_obs) {
        cdid_ipw_point_view_abort(
          "Streamed risk-set observer count exceeds `link_view$support$n_observed_risk_units`.",
          c("cdid_ipw_point_view_count_error", "cdid_ipw_owner_integrity_error")
        )
      }
      a_ch <- as.integer(ch$observed_risk$a_index)
      gi_ch <- as.integer(ch$observed_risk$cohort_index)
      loc_ch <- as.integer(ch$observed_risk$local_pos)
      w_ch <- as.double(ch$weight[loc_ch])
      yl_ch <- as.double(ch$observed_risk$y_left_obs)
      yr_ch <- as.double(ch$observed_risk$y_right_obs)
      Z_ch <- ch$design$Z_obs

      if (length(a_ch) != k_ch || length(gi_ch) != k_ch || length(w_ch) != k_ch ||
          length(yl_ch) != k_ch || length(yr_ch) != k_ch ||
          !is.matrix(Z_ch) || nrow(Z_ch) != k_ch || ncol(Z_ch) != p) {
        cdid_ipw_point_view_abort(
          "Malformed observed-risk slice in link view chunk.",
          c("cdid_ipw_point_view_chunk_error", "cdid_ipw_owner_integrity_error")
        )
      }
      if (a_ch[1L] <= last_a_idx || (k_ch > 1L && any(diff(a_ch) <= 0L)) ||
          a_ch[1L] < cur_s || a_ch[k_ch] > cur_e) {
        cdid_ipw_point_view_abort(
          "Streamed risk-set observer `a_index` is not strictly increasing in global A-unit order.",
          c("cdid_ipw_point_view_order_error", "cdid_ipw_owner_order_mismatch_error")
        )
      }
      if (anyNA(gi_ch) || any(gi_ch < 1L | gi_ch > m)) {
        cdid_ipw_point_view_abort(
          "Streamed risk-set observer cohort index is out of `[1, m_cohorts]`.",
          c("cdid_ipw_point_view_cohort_error", "cdid_ipw_owner_cohort_error")
        )
      }
      dy_ch <- yr_ch - yl_ch
      if (anyNA(dy_ch) || any(!is.finite(dy_ch))) {
        diag_env$error_diagnostics <- list(
          stage = "pair_difference_check",
          chunk_start = cur_s,
          chunk_end = cur_e,
          failed_a_index = a_ch,
          failed_dy = dy_ch
        )
        if (!is.null(failure_trace_env)) {
          failure_trace_env$error_diagnostics <- diag_env$error_diagnostics
          failure_trace_env$observer_released_on_error <- TRUE
        }
        rm(a_ch, gi_ch, loc_ch, w_ch, yl_ch, yr_ch, dy_ch, Z_ch)
        ch <- NULL
        cleanup_temporary_buffers()
        cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
      }

      dst <- seq.int(write_pos + 1L, write_pos + k_ch)
      obs_buf_env$risk_a_idx[dst] <- a_ch
      obs_buf_env$risk_gi[dst] <- gi_ch
      obs_buf_env$risk_w[dst] <- w_ch
      obs_buf_env$risk_dy[dst] <- dy_ch
      obs_buf_env$Z_use[dst, ] <- Z_ch
      write_pos <- write_pos + k_ch
      last_a_idx <- a_ch[k_ch]
      rm(a_ch, gi_ch, loc_ch, w_ch, yl_ch, yr_ch, dy_ch, Z_ch, dst)
    }
    stream_ops <- stream_ops + ch$work$actual_work_ops + k_ch * (p + 2)

    # Release consumed query chunk BEFORE callback and next query!
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
        endpoints_read = endpoints_read_total
      ))
    }

    cur_s <- as.double(cur_e) + 1
  }

  if (write_pos != n_obs || obs_buf_env$active_query_chunks != 0L ||
      !ipw_point_view_identical(endpoints_read_total, as.integer(2L * n_obs))) {
    cdid_ipw_point_view_abort(
      "Streamed risk-set observer count or endpoint read count disagrees with `link_view$support$n_observed_risk_units`.",
      c("cdid_ipw_point_view_count_error", "cdid_ipw_owner_integrity_error")
    )
  }

  # 5. Reconstruct the worker's point arrays from exact fitted coefficients,
  #    rates, and normalizers, and invoke unchanged `ipw_link_means` and
  #    `ipw_link_partials` without refitting the propensity model.
  sup <- list(
    n = as.double(na),
    g = fitted_view$cohorts$g,
    K = fitted_view$cohorts$K,
    cohorts = fitted_view$cohorts$cohorts,
    cohort_labels = fitted_view$cohorts$cohort_labels,
    m_cohorts = m,
    group_count = as.integer(unname(fitted_view$counts$by_cohort_a)),
    group_a_count = as.integer(unname(fitted_view$counts$by_cohort_a)),
    group_mass = as.double(unname(fitted_view$masses$cohort_a_weight)),
    group_a_mass = as.double(unname(fitted_view$masses$cohort_a_weight)),
    q = as.double(unname(fitted_view$masses$q)),
    a_weight_total = as.double(fitted_view$masses$a_weight_total)
  )

  tryCatch(
    {
      pr_chk <- ipw_pair_rates(obs_buf_env$risk_gi, obs_buf_env$risk_w, sup)
      if (!ipw_point_view_identical(as.integer(pr_chk$pair_count), as.integer(unname(fitted_view$counts$by_cohort_pair))) ||
          !ipw_point_view_identical(as.double(pr_chk$pair_mass), as.double(unname(fitted_view$masses$cohort_pair_weight))) ||
          !ipw_point_view_identical(as.double(pr_chk$rates), as.double(unname(fitted_view$masses$rates)))) {
        diag_env$error_diagnostics <- list(
          stage = "reconstructed_rates_check",
          reconstructed_rates = pr_chk$rates,
          fitted_rates = unname(fitted_view$masses$rates)
        )
        cdid_ipw_point_view_abort(
          "Reconstructed observer-buffer pair masses or rates disagree with the sealed fitted nuisance view.",
          c("cdid_ipw_point_view_integrity_error", "cdid_ipw_shared_mask_error")
        )
      }

      rates <- as.double(unname(fitted_view$masses$rates))
      beta <- as.double(fitted_view$coefficients)
      names(beta) <- colnames(frame_owner$Z_a)

      obs_buf_env$r_use <- ipw_clip_floor(rates)[obs_buf_env$risk_gi]
      obs_buf_env$fit_weights <- obs_buf_env$risk_w / obs_buf_env$r_use
      if (any(!is.finite(obs_buf_env$fit_weights) | obs_buf_env$fit_weights <= 0)) {
        diag_env$error_diagnostics <- list(
          stage = "inverse_rate_weights",
          rates = rates,
          fit_weights = obs_buf_env$fit_weights
        )
        cdid_abort(
          "Nonfinite or nonpositive inverse-rate propensity weights `w_i / r_G`.",
          c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
        )
      }

      obs_buf_env$lv <- ipw_logit_values(obs_buf_env$Z_use, beta)
      if (any(!is.finite(obs_buf_env$lv$eta)) || any(!is.finite(obs_buf_env$lv$odds)) ||
          any(obs_buf_env$lv$odds <= 0) || any(!is.finite(obs_buf_env$lv$prob))) {
        diag_env$error_diagnostics <- list(
          stage = "logit_values_check",
          failed_coefficients = beta
        )
        cdid_abort(
          "Weighted propensity logit produced nonfinite linear predictors, odds, or probabilities.",
          c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
        )
      }

      obs_buf_env$ow <- ipw_obs_weights(
        obs_buf_env$fit_weights,
        obs_buf_env$lv$odds,
        obs_buf_env$risk_gi,
        sup$n
      )
      if (!ipw_point_view_identical(as.integer(sum(obs_buf_env$ow$is_t)), as.integer(fitted_view$counts$n_treated_obs)) ||
          !ipw_point_view_identical(as.integer(sum(obs_buf_env$ow$is_c)), as.integer(fitted_view$counts$n_control_obs)) ||
          !ipw_point_view_identical(as.double(obs_buf_env$ow$sum_a), as.double(fitted_view$means$sum_a)) ||
          !ipw_point_view_identical(as.double(obs_buf_env$ow$mean_a), as.double(fitted_view$means$mean_a)) ||
          !ipw_point_view_identical(as.double(obs_buf_env$ow$sum_b), as.double(fitted_view$means$sum_b)) ||
          !ipw_point_view_identical(as.double(obs_buf_env$ow$mean_b), as.double(fitted_view$means$mean_b))) {
        diag_env$error_diagnostics <- list(
          stage = "reconstructed_normalizer_check",
          reconstructed_sum_a = obs_buf_env$ow$sum_a,
          fitted_sum_a = fitted_view$means$sum_a,
          reconstructed_sum_b = obs_buf_env$ow$sum_b,
          fitted_sum_b = fitted_view$means$sum_b
        )
        cdid_ipw_point_view_abort(
          "Reconstructed observer normalizers disagree with the sealed fitted nuisance view.",
          c("cdid_ipw_point_view_integrity_error", "cdid_ipw_shared_mask_error")
        )
      }

      obs_buf_env$st <- list(
        n = sup$n,
        m_cohorts = m,
        p = p,
        risk_a_idx = obs_buf_env$risk_a_idx,
        risk_gi = obs_buf_env$risk_gi,
        risk_w = obs_buf_env$risk_w,
        pair_count = as.integer(unname(fitted_view$counts$by_cohort_pair)),
        pair_mass = as.double(unname(fitted_view$masses$cohort_pair_weight)),
        rates = rates,
        Z_use = obs_buf_env$Z_use,
        D_use = as.double(obs_buf_env$risk_gi == 1L),
        fit_weights = obs_buf_env$fit_weights,
        beta = beta,
        odds_use = obs_buf_env$lv$odds,
        prob_use = obs_buf_env$lv$prob,
        is_t = obs_buf_env$ow$is_t,
        is_c = obs_buf_env$ow$is_c,
        a_obs = obs_buf_env$ow$a_obs,
        b_obs = obs_buf_env$ow$b_obs,
        sum_a = fitted_view$means$sum_a,
        mean_a = fitted_view$means$mean_a,
        sum_b = fitted_view$means$sum_b,
        mean_b = fitted_view$means$mean_b
      )
      obs_buf_env$hs <- list(
        H = fitted_view$factors$H,
        H_chol = fitted_view$factors$H_chol,
        A_beta_r = fitted_view$partials$A_beta_r
      )

      obs_buf_env$means <- tryCatch(
        ipw_link_means(obs_buf_env$st, obs_buf_env$risk_dy),
        error = function(e) {
          diag_env$error_diagnostics <- list(
            stage = "link_means_stage",
            sum_a = fitted_view$means$sum_a,
            sum_b = fitted_view$means$sum_b
          )
          stop(e)
        }
      )
      obs_buf_env$part <- tryCatch(
        ipw_link_partials(
          obs_buf_env$st,
          obs_buf_env$hs,
          obs_buf_env$risk_dy,
          obs_buf_env$means,
          sup$q,
          sup$cohort_labels
        ),
        error = function(e) {
          diag_env$error_diagnostics <- list(
            stage = "link_partials_stage",
            treated_mean = obs_buf_env$means$Tmean,
            control_mean = obs_buf_env$means$Cmean,
            estimate = obs_buf_env$means$mu
          )
          stop(e)
        }
      )
    },
    error = function(cnd) {
      err_diag <- if (exists("error_diagnostics", envir = diag_env, inherits = FALSE)) {
        diag_env$error_diagnostics
      } else {
        list(stage = "unknown_point_link_failure")
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

  # 6. Extract ONLY compact O(p^2 + p*m + m) quantities, then immediately scrub
  #    the observer buffer and all per-observer model vectors BEFORE sealing/returning.
  Tmean_val <- as.double(obs_buf_env$means$Tmean)
  Cmean_val <- as.double(obs_buf_env$means$Cmean)
  mu_val <- as.double(obs_buf_env$means$mu)

  G_beta_val <- obs_buf_env$part$G_beta
  G_r_val <- obs_buf_env$part$G_r
  loading_beta_val <- obs_buf_env$part$loading_beta
  loading_rate_val <- obs_buf_env$part$loading_rate
  rate_coef_s0_val <- obs_buf_env$part$rate_coef_s0
  rate_coef_s1_val <- obs_buf_env$part$rate_coef_s1

  # Scrub all observer/model arrays NOW before candidate checks, callbacks, or sealing!
  cleanup_temporary_buffers()

  if (!is.null(link_view$loading_beta_candidate) &&
      !ipw_point_view_identical(as.double(link_view$loading_beta_candidate), as.double(loading_beta_val))) {
    cdid_ipw_point_view_abort(
      "Candidate `loading_beta_candidate` stored on `link_view` disagrees with the fitted point link's `loading_beta`.",
      c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error")
    )
  }

  if (!is.null(callback)) {
    callback(list(
      phase = "after_observer_release",
      active_query_chunks = obs_buf_env$active_query_chunks,
      observer_active = obs_buf_env$active,
      observer_released = obs_buf_env$observer_released,
      endpoints_read = endpoints_read_total
    ))
  }

  link_pair_mass <- as.double(link_view$support$by_cohort_pair_mass)
  names(link_pair_mass) <- sup$cohort_labels
  link_rates <- as.double(link_view$support$rates)
  names(link_rates) <- sup$cohort_labels

  shared_across_dates <- !ipw_point_view_identical(link_view$pair$s, fitted_view$defining_pair$s) ||
    !ipw_point_view_identical(link_view$pair$t, fitted_view$defining_pair$t)

  actual_ops <- stream_ops + 16 * n_obs * p + 16 * n_obs * m + 4 * (p * p + p * m)
  stopifnot(actual_ops <= preflight$total_work_ops_bound)

  point_link <- structure(
    list(
      link_id = eff_link_id,
      functional_key = eff_func_key,
      nuisance_id = fitted_view$nuisance_id,
      fitted_nuisance_ref = fitted_view$ref,
      nuisance_view_ref = fitted_view$nuisance_view_ref,
      nuisance_ref = fitted_view$nuisance_view_ref,
      link_view_ref = link_view$ref,
      link_ref = link_view$ref,
      frame_ref = fitted_view$frame_ref,
      panel_ref = fitted_view$panel_ref,
      estimate = mu_val,
      treated_mean = Tmean_val,
      control_mean = Cmean_val,
      pair = list(
        s = link_view$pair$s,
        t = link_view$pair$t,
        s_calendar_index = link_view$pair$s_calendar_index,
        t_calendar_index = link_view$pair$t_calendar_index,
        link_id = eff_link_id,
        functional_key = eff_func_key
      ),
      defining_pair = fitted_view$defining_pair,
      shared_across_dates = shared_across_dates,
      nuisance_mask_verified = TRUE,
      cohorts = fitted_view$cohorts,
      dimensions = fitted_view$dimensions,
      counts = list(
        n_a_units = na,
        n_risk_a_units = fitted_view$counts$n_risk_a_units,
        n_outside_risk_a_units = fitted_view$counts$n_outside_risk_a_units,
        n_observed_pair_units = as.integer(link_view$support$n_observed_pair_all),
        n_observed_risk_units = n_obs,
        n_nonobserved_risk_units = fitted_view$counts$n_nonobserved_risk_units,
        n_treated_obs = fitted_view$counts$n_treated_obs,
        n_control_obs = fitted_view$counts$n_control_obs,
        by_cohort_a = fitted_view$counts$by_cohort_a,
        by_cohort_pair = fitted_view$counts$by_cohort_pair
      ),
      masses = fitted_view$masses,
      means = list(
        treated_mean = Tmean_val,
        control_mean = Cmean_val,
        estimate = mu_val,
        sum_a = fitted_view$means$sum_a,
        mean_a = fitted_view$means$mean_a,
        sum_b = fitted_view$means$sum_b,
        mean_b = fitted_view$means$mean_b
      ),
      coefficients = fitted_view$coefficients,
      beta = fitted_view$beta,
      factors = fitted_view$factors,
      partials = list(
        G_beta = G_beta_val,
        G_r = G_r_val,
        A_beta_r = fitted_view$partials$A_beta_r
      ),
      loadings = list(
        loading_beta = loading_beta_val,
        loading_rate = loading_rate_val,
        rate_coef_s0 = rate_coef_s0_val,
        rate_coef_s1 = rate_coef_s1_val
      ),
      diagnostics = fitted_view$diagnostics,
      point_summary_disclosure = list(
        defining_rate_source = fitted_view$point_summary_disclosure$defining_rate_source,
        observer_cohort_a_mass = fitted_view$point_summary_disclosure$observer_cohort_a_mass,
        observer_cohort_pair_mass = fitted_view$point_summary_disclosure$observer_cohort_pair_mass,
        observer_rates = fitted_view$point_summary_disclosure$observer_rates,
        view_chunked_cohort_a_mass = fitted_view$point_summary_disclosure$view_chunked_cohort_a_mass,
        view_chunked_cohort_pair_mass = fitted_view$point_summary_disclosure$view_chunked_cohort_pair_mass,
        view_chunked_rates = fitted_view$point_summary_disclosure$view_chunked_rates,
        link_chunked_cohort_pair_mass = link_pair_mass,
        link_chunked_rates = link_rates,
        max_abs_pair_mass_chunk_diff = max(abs(fitted_view$masses$cohort_pair_weight - link_pair_mass)),
        max_abs_rate_chunk_diff = max(abs(fitted_view$masses$rates - link_rates)),
        propensity_refit_performed = FALSE,
        candidate_point_summary_only = TRUE,
        certifies_root_or_inference = FALSE,
        glm_convergence_is_certificate = FALSE,
        seal_is_certificate = FALSE,
        disclosure = paste0(
          "Compact point link reuses exact fitted coefficients, observer-buffer rates, normalizers, ",
          "Hessian Cholesky factor, and score-rate Jacobian without refitting the propensity model. ",
          "Chunked link-view support summaries are disclosed alongside observer-buffer summaries. ",
          "Sealing records a finite binary64 point candidate only and never certifies raw normalization, ",
          "root enclosure, covariance transfer, or statistical inference."
        )
      ),
      glm_control = fitted_view$glm_control,
      control = ctrl,
      preflight = preflight,
      accounting = list(
        retained_matrix_entries = preflight$retained_matrix_entries_per_link,
        retained_vector_entries = preflight$retained_vector_entries_per_link,
        total_retained_entries = preflight$retained_entries_per_link,
        borrowed_fitted_nuisance_ref = TRUE,
        borrowed_nuisance_view_ref = TRUE,
        borrowed_link_view_ref = TRUE,
        borrowed_frame_ref = TRUE,
        borrowed_panel_ref = TRUE,
        observer_buffer_released = TRUE,
        reader_chunks_released = TRUE,
        full_n_masks_retained = FALSE,
        observer_arrays_retained = FALSE,
        copied_Z_retained = FALSE,
        per_link_endpoint_copies_retained = FALSE,
        propensity_refit_performed = FALSE,
        endpoints_read = endpoints_read_total,
        peak_rss_promised = FALSE
      ),
      work = list(
        chunks_queried = logit_cert_count_output(chunks_queried),
        mask_sharing_chunks_compared = logit_cert_count_output(sharing_check$chunks_compared),
        max_simultaneous_query_chunks = 1L,
        max_simultaneous_observer_buffers = 1L,
        observer_units_copied = logit_cert_count_output(n_obs),
        endpoints_read = endpoints_read_total,
        propensity_solver_calls = 0L,
        actual_work_ops = actual_ops
      )
    ),
    class = "cdid_ipw_point_link_view"
  )

  ipw_point_view_verify_compact_storage(point_link, p_max = p, m_max = m)
  actual_payload <- ipw_point_view_payload_budget(point_link, ctrl)
  if (actual_payload$entries > preflight$retained_entries_per_link) {
    cdid_ipw_point_view_abort("Point-link output exceeded its proved retained envelope.",
      "cdid_ipw_point_view_bound_error")
  }
  ipw_owner_seal(point_link, "point_link_view")
}

construct_ipw_point_link_view <- fit_ipw_point_link_view

fit_ipw_point_link_from_owner <- function(
  fitted_view,
  s,
  t,
  control,
  g = NULL,
  K = NULL,
  glm_control = NULL,
  link_id = NULL,
  functional_key = NULL,
  callback = NULL,
  trace_env = NULL,
  failure_trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(fitted_view) || missing(s) || missing(t) || missing(control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`fitted_view`, `s`, `t`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  ctrl <- validate_ipw_point_view_control(control)
  ipw_point_view_seal_gate(fitted_view, "fitted_nuisance_view")
  payload <- ipw_point_view_payload_budget(list(fitted_view, s, t, g, K,
    glm_control, link_id, functional_key), ctrl)
  validate_ipw_fitted_nuisance_view(fitted_view)
  owner_ctrl <- ipw_point_view_to_owner_control(ctrl)
  frame_owner <- get("frame_owner", envir = fitted_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = fitted_view$panel_ref, inherits = FALSE)
  nv_obj <- get("nuisance_view", envir = fitted_view$nuisance_view_ref, inherits = FALSE)
  n <- frame_owner$n_a_units; p <- frame_owner$n_covariates
  m <- fitted_view$cohorts$m_cohorts; no <- fitted_view$counts$n_observed_risk_units
  if (no > ctrl$max_observer_units) {
    cdid_ipw_point_view_abort("Observed risk support exceeds explicit point observer cap.",
      "cdid_ipw_owner_bound_error")
  }
  b <- min(n, ctrl$max_chunk_rows); chunks <- ceiling(n / b)
  vp <- preflight_ipw_link_view(n, p, m, owner_ctrl, has_candidates = FALSE)
  # Shape-only reservation before the constructor's posting/support/mask scans.
  # A compatible link has exactly the fitted observer count, but this reservation
  # is not used as proof of equality; the actual mask check remains mandatory.
  ipw_owner_safe_sum(c(vp$work_ops_bound, chunks * vp$work_ops_bound,
    n * 512 * (p + m + 8), 16 * payload$work, 16 * p^3, 8192),
    ctrl$max_work_ops, "from-owner construction plus point fitting whole work")
  ipw_owner_safe_sum(c(payload$entries, ipw_point_view_owner_entries(frame_owner,panel_owner),
    no * (16 * p + 128), b * (16 * p + 256), 16 * (p*p + p*m), 16384),
    ctrl$max_matrix_entries, "from-owner construction plus point fitting concurrent entries")
  ipw_owner_safe_sum(c(payload$entries, 8192 + 4*p*p + 2*p*m + 16*p + 32*m),
    min(ctrl$max_retained_entries,ctrl$max_output_entries), "from-owner compact output envelope")

  lk_view <- construct_ipw_link_view(
    frame_owner = frame_owner,
    s = s,
    t = t,
    g = g,
    K = K,
    nuisance_view = nv_obj,
    control = owner_ctrl,
    link_id = link_id,
    functional_key = functional_key,
    require_support = TRUE,
    verify_nuisance_mask = TRUE
  )

  pt_pre <- preflight_ipw_point_link_view(
    fitted_view = fitted_view,
    link_view = lk_view,
    control = ctrl,
    glm_control = glm_control,
    n_planned_links = 1L
  )
  ipw_owner_safe_sum(
    c(
      lk_view$preflight$work_ops_bound,
      pt_pre$total_work_ops_bound,
      pt_pre$character_metadata$work_ops_bound,
      128
    ),
    ctrl$max_work_ops,
    "composed link-view construction and point-link fitting work"
  )

  fit_ipw_point_link_view(
    fitted_view = fitted_view,
    link_view = lk_view,
    control = ctrl,
    glm_control = glm_control,
    link_id = link_id,
    functional_key = functional_key,
    callback = callback,
    trace_env = trace_env,
    failure_trace_env = failure_trace_env
  )
}

fit_ipw_point_link_view_set <- function(
  fitted_view,
  link_views,
  control,
  glm_control = NULL,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(fitted_view) || missing(link_views) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`fitted_view`, `link_views`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.list(link_views) || length(link_views) < 1L) {
    cdid_ipw_point_view_abort(
      "`link_views` must be a non-empty list of sealed `cdid_ipw_link_view` objects.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }

  ctrl <- validate_ipw_point_view_control(control)
  n_links <- length(link_views)
  set_payload <- ipw_point_view_payload_budget(list(link_views, names(link_views)), ctrl)

  single_fitted <- inherits(fitted_view, "cdid_ipw_fitted_nuisance_view")
  if (!single_fitted && (!is.list(fitted_view) || length(fitted_view) != n_links)) {
    cdid_ipw_point_view_abort(
      "`fitted_view` must be a sealed `cdid_ipw_fitted_nuisance_view` or a list of the same length as `link_views`.",
      "cdid_ipw_point_view_input_error"
    )
  }

  set_header <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(16, n_links), 128),
    min(ctrl$max_retained_entries, ctrl$max_matrix_entries, ctrl$max_output_entries),
    "point-link set list/trace/header payload before allocation"
  )
  set_work <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(128, n_links), 128, 2 * set_payload$work),
    ctrl$max_work_ops,
    "point-link set header work"
  )
  sum_ret <- ipw_owner_safe_sum(c(set_header, set_payload$entries),
    min(ctrl$max_retained_entries, ctrl$max_output_entries), "point-link set metadata and outputs")
  peak_extra <- 0
  for (i in seq_len(n_links)) {
    fv_i <- if (single_fitted) fitted_view else fitted_view[[i]]
    pi <- preflight_ipw_point_link_view(
      fitted_view = fv_i,
      link_view = link_views[[i]],
      control = ctrl,
      glm_control = glm_control,
      n_planned_links = 1L
    )
    sum_ret <- ipw_owner_safe_sum(
      c(sum_ret, pi$total_retained_entries),
      min(ctrl$max_retained_entries, ctrl$max_output_entries),
      "complete point-link set retained outputs"
    )
    set_work <- ipw_owner_safe_sum(
      c(set_work, pi$total_work_ops_bound, pi$character_metadata$work_ops_bound, 128),
      ctrl$max_work_ops,
      "complete point-link set work"
    )
    peak_extra <- max(peak_extra, pi$peak_concurrent_entries - pi$total_retained_entries)
  }
  ipw_owner_safe_sum(
    c(sum_ret, peak_extra),
    ctrl$max_matrix_entries,
    "all point-link set outputs plus largest active workspace"
  )
  pi <- NULL

  out <- vector("list", n_links)
  names(out) <- names(link_views)
  between_link_active_buffers <- integer(n_links)
  list_names <- names(link_views)

  for (i in seq_len(n_links)) {
    fv_i <- if (single_fitted) fitted_view else fitted_view[[i]]
    lid_i <- if (!is.null(list_names) && nzchar(list_names[i]) && is.null(link_views[[i]]$link_id)) {
      list_names[i]
    } else NULL
    step_trace <- new.env(parent = emptyenv())
    out[[i]] <- fit_ipw_point_link_view(
      fitted_view = fv_i,
      link_view = link_views[[i]],
      control = ctrl,
      glm_control = glm_control,
      link_id = lid_i,
      trace_env = step_trace
    )
    if (!ipw_point_view_identical(step_trace$active_observer_buffers, 0L) ||
        !isTRUE(step_trace$observer_released) ||
        !ipw_point_view_identical(step_trace$propensity_solver_calls, 0L)) {
      cdid_ipw_point_view_abort(
        "Temporary observer buffer was not released before starting the next point link.",
        "cdid_ipw_point_view_lifecycle_error"
      )
    }
    between_link_active_buffers[i] <- as.integer(step_trace$active_observer_buffers)
  }

  if (!is.null(trace_env)) {
    trace_env$n_links_fitted <- n_links
    trace_env$between_link_active_buffers <- between_link_active_buffers
    trace_env$all_released_between_links <- all(between_link_active_buffers == 0L)
    trace_env$propensity_solver_calls <- 0L
  }

  out
}

validate_ipw_point_link_view <- function(
  point_link,
  fitted_view = NULL,
  link_view = NULL,
  frame_owner = NULL
) {
  ipw_point_view_require_deps()
  if (missing(point_link) || is.null(point_link) ||
      !inherits(point_link, "cdid_ipw_point_link_view") ||
      !is.list(point_link) ||
      !ipw_owner_seal_ok(point_link, "point_link_view")) {
    cdid_ipw_point_view_abort(
      "Point link view is malformed, unsealed, deserialized without revalidation, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_ipw_point_view_integrity_error", "cdid_ipw_owner_integrity_error", "cdid_ipw_owner_view_error")
    )
  }
  ipw_owner_refuse_decorative_attrs(point_link, "point_link_view")

  fv_ref <- point_link$fitted_nuisance_ref
  if (!is.environment(fv_ref) || !environmentIsLocked(fv_ref) ||
      !exists("fitted_nuisance_view", envir = fv_ref, inherits = FALSE) ||
      !bindingIsLocked("fitted_nuisance_view", fv_ref)) {
    cdid_ipw_point_view_abort(
      "Point link view lacks a valid sealed `fitted_nuisance_view` reference.",
      c("cdid_ipw_point_view_integrity_error", "cdid_ipw_owner_integrity_error")
    )
  }
  fv_obj <- get("fitted_nuisance_view", envir = fv_ref, inherits = FALSE)
  validate_ipw_fitted_nuisance_view(fv_obj)

  lv_ref <- point_link$link_view_ref
  if (!is.environment(lv_ref) || !environmentIsLocked(lv_ref) ||
      !exists("link_view", envir = lv_ref, inherits = FALSE) ||
      !bindingIsLocked("link_view", lv_ref) ||
      !ipw_point_view_identical(point_link$link_ref, lv_ref)) {
    cdid_ipw_point_view_abort(
      "Point link view lacks a valid sealed `link_view` reference.",
      c("cdid_ipw_point_view_integrity_error", "cdid_ipw_owner_integrity_error")
    )
  }
  lv_obj <- get("link_view", envir = lv_ref, inherits = FALSE)
  f_obj <- get("frame_owner", envir = fv_obj$frame_ref, inherits = FALSE)
  nv_obj <- get("nuisance_view", envir = fv_obj$nuisance_view_ref, inherits = FALSE)
  validate_ipw_link_view(
    lv_obj,
    frame_owner = f_obj,
    nuisance_view = if (!is.null(lv_obj$nuisance_ref)) nv_obj else NULL
  )

  if (!ipw_point_view_identical(point_link$frame_ref, fv_obj$frame_ref) ||
      !ipw_point_view_identical(point_link$panel_ref, fv_obj$panel_ref) ||
      !ipw_point_view_identical(point_link$nuisance_view_ref, fv_obj$nuisance_view_ref) ||
      !ipw_point_view_identical(point_link$nuisance_ref, fv_obj$nuisance_view_ref) ||
      !ipw_point_view_identical(point_link$frame_ref, lv_obj$frame_ref) ||
      !ipw_point_view_identical(point_link$panel_ref, lv_obj$panel_ref) ||
      !isTRUE(point_link$nuisance_mask_verified)) {
    cdid_ipw_point_view_abort(
      "Point link view owner or nuisance references disagree with its bound views.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_point_view_integrity_error")
    )
  }

  if (!is.null(fitted_view)) {
    validate_ipw_fitted_nuisance_view(fitted_view)
    if (!ipw_point_view_identical(fv_ref, fitted_view$ref)) {
      cdid_ipw_point_view_abort(
        "Point link view belongs to a different `cdid_ipw_fitted_nuisance_view` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  }
  if (!is.null(link_view)) {
    validate_ipw_link_view(link_view)
    if (!ipw_point_view_identical(lv_ref, link_view$ref)) {
      cdid_ipw_point_view_abort(
        "Point link view belongs to a different `cdid_ipw_link_view` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  }
  if (!is.null(frame_owner)) {
    validate_ipw_frame_owner(frame_owner)
    if (!ipw_point_view_identical(point_link$frame_ref, frame_owner$ref)) {
      cdid_ipw_point_view_abort(
        "Point link view belongs to a different `cdid_ipw_frame_owner` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_context_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  }

  validate_ipw_glm_control(point_link$glm_control)
  validate_ipw_point_view_control(point_link$control)

  p <- as.integer(f_obj$n_covariates)
  m <- as.integer(nv_obj$cohorts$m_cohorts)
  ipw_point_view_verify_compact_storage(point_link, p_max = p, m_max = m)

  if (!ipw_point_view_identical(point_link$coefficients, fv_obj$coefficients) ||
      !ipw_point_view_identical(point_link$beta, fv_obj$beta) ||
      !ipw_point_view_identical(point_link$factors, fv_obj$factors) ||
      !ipw_point_view_identical(point_link$partials$A_beta_r, fv_obj$partials$A_beta_r) ||
      !ipw_point_view_identical(point_link$masses, fv_obj$masses) ||
      !ipw_point_view_identical(point_link$means$sum_a, fv_obj$means$sum_a) ||
      !ipw_point_view_identical(point_link$means$mean_a, fv_obj$means$mean_a) ||
      !ipw_point_view_identical(point_link$means$sum_b, fv_obj$means$sum_b) ||
      !ipw_point_view_identical(point_link$means$mean_b, fv_obj$means$mean_b)) {
    cdid_ipw_point_view_abort(
      "Point link view coefficients, factors, rates, or normalizers disagree with its bound `fitted_nuisance_view`.",
      "cdid_ipw_point_view_integrity_error"
    )
  }

  if (!gmm_scalar_numeric(point_link$estimate) ||
      !gmm_scalar_numeric(point_link$treated_mean) ||
      !gmm_scalar_numeric(point_link$control_mean) ||
      !ipw_point_view_identical(point_link$estimate, point_link$means$estimate) ||
      !ipw_point_view_identical(point_link$treated_mean, point_link$means$treated_mean) ||
      !ipw_point_view_identical(point_link$control_mean, point_link$means$control_mean)) {
    cdid_ipw_point_view_abort(
      "Point link view scalar means/estimate are malformed or inconsistent.",
      "cdid_ipw_point_view_integrity_error"
    )
  }

  G_beta <- point_link$partials$G_beta
  G_r <- point_link$partials$G_r
  l_beta <- point_link$loadings$loading_beta
  l_rate <- point_link$loadings$loading_rate
  rc0 <- point_link$loadings$rate_coef_s0
  rc1 <- point_link$loadings$rate_coef_s1
  if (!is.double(G_beta) || length(G_beta) != p || any(!is.finite(G_beta)) ||
      !is.double(l_beta) || length(l_beta) != p || any(!is.finite(l_beta)) ||
      !is.double(G_r) || length(G_r) != m || any(!is.finite(G_r)) ||
      !is.double(l_rate) || length(l_rate) != m || any(!is.finite(l_rate)) ||
      !is.double(rc0) || length(rc0) != m || any(!is.finite(rc0)) ||
      !is.double(rc1) || length(rc1) != m || any(!is.finite(rc1))) {
    cdid_ipw_point_view_abort(
      "Point link view outcome partials or loadings are malformed or non-finite.",
      "cdid_ipw_point_view_integrity_error"
    )
  }

  psd <- point_link$point_summary_disclosure
  if (!is.list(psd) || !isTRUE(psd$candidate_point_summary_only) ||
      !ipw_point_view_identical(psd$propensity_refit_performed, FALSE) ||
      !ipw_point_view_identical(psd$certifies_root_or_inference, FALSE) ||
      !ipw_point_view_identical(psd$glm_convergence_is_certificate, FALSE) ||
      !ipw_point_view_identical(psd$seal_is_certificate, FALSE)) {
    cdid_ipw_point_view_abort(
      "Point link view non-certification disclosure is missing or altered.",
      "cdid_ipw_point_view_integrity_error"
    )
  }

  invisible(point_link)
}

extract_ipw_point_link_snapshot <- function(point_link) {
  validate_ipw_point_link_view(point_link)
  structure(
    list(
      link_id = point_link$link_id,
      functional_key = point_link$functional_key,
      nuisance_id = point_link$nuisance_id,
      estimate = point_link$estimate,
      treated_mean = point_link$treated_mean,
      control_mean = point_link$control_mean,
      pair = point_link$pair,
      defining_pair = point_link$defining_pair,
      shared_across_dates = point_link$shared_across_dates,
      nuisance_mask_verified = point_link$nuisance_mask_verified,
      cohorts = point_link$cohorts,
      dimensions = point_link$dimensions,
      counts = point_link$counts,
      masses = point_link$masses,
      means = point_link$means,
      coefficients = point_link$coefficients,
      beta = point_link$beta,
      factors = point_link$factors,
      partials = point_link$partials,
      loadings = point_link$loadings,
      diagnostics = point_link$diagnostics,
      point_summary_disclosure = point_link$point_summary_disclosure,
      glm_control = point_link$glm_control,
      control = point_link$control
    ),
    class = "cdid_ipw_point_link_snapshot"
  )
}

ipw_point_link_snapshot_payload <- function(object, max_entries, max_work) {
  entries <- 0
  work <- 0
  visit <- function(x, depth = 0L) {
    if (depth > 4L) {
      cdid_ipw_point_view_abort("Point-link snapshot exceeds compact schema nesting.", "cdid_ipw_point_view_revalidation_error")
    }
    if (is.environment(x) || is.null(x)) return(invisible(NULL))
    k <- length(x)
    entries <<- ipw_owner_safe_sum(c(entries, k, 1), max_entries, "point-link snapshot reachable payload before content checks")
    work <<- ipw_owner_safe_sum(c(work, k, 1), max_work, "point-link snapshot traversal work")
    if (!is.null(names(x))) visit(names(x), depth + 1L)
    if (is.list(x)) {
      for (j in seq_len(k)) visit(x[[j]], depth + 1L)
    } else if (is.character(x)) {
      for (j in seq_len(k)) {
        nb <- nchar(x[j], type = "bytes")
        if (is.na(nb)) nb <- 0
        entries <<- ipw_owner_safe_sum(c(entries, ceiling(nb / 2)), max_entries, "point-link snapshot character payload")
        work <<- ipw_owner_safe_sum(c(work, nb), max_work, "point-link snapshot character work")
      }
    }
    invisible(NULL)
  }
  visit(object)
  list(entries = entries, work = work)
}

revalidate_ipw_point_link_view <- function(
  object,
  fitted_view,
  link_view,
  control,
  glm_control = NULL
) {
  ipw_point_view_require_deps()
  if (missing(object) || missing(fitted_view) || missing(link_view) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`object`, `fitted_view`, `link_view`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  ipw_owner_refuse_decorative_attrs(object, "object")
  if (!is.list(object) ||
      (!inherits(object, "cdid_ipw_point_link_view") &&
       !inherits(object, "cdid_ipw_point_link_snapshot"))) {
    cdid_ipw_point_view_abort(
      "Expected a `cdid_ipw_point_link_view` or `cdid_ipw_point_link_snapshot` for revalidation.",
      c("cdid_ipw_point_view_revalidation_error", "cdid_ipw_owner_integrity_error")
    )
  }
  ctrl <- validate_ipw_point_view_control(control)
  payload <- ipw_point_view_payload_budget(list(object), ctrl)
  ipw_owner_safe_sum(c(4 * payload$work), ctrl$max_work_ops, "snapshot validation/comparison work")
  snapshot_fields <- c(
    "link_id", "functional_key", "nuisance_id", "estimate", "treated_mean", "control_mean",
    "pair", "defining_pair", "shared_across_dates", "nuisance_mask_verified", "cohorts",
    "dimensions", "counts", "masses", "means", "coefficients", "beta", "factors",
    "partials", "loadings", "diagnostics", "point_summary_disclosure", "glm_control", "control"
  )
  full_fields <- c(
    snapshot_fields,
    "fitted_nuisance_ref", "nuisance_view_ref", "nuisance_ref",
    "link_view_ref", "link_ref", "frame_ref", "panel_ref",
    "preflight", "accounting", "work", "ref"
  )
  expected <- if (inherits(object, "cdid_ipw_point_link_snapshot")) snapshot_fields else full_fields
  if (length(object) != length(expected) || is.null(names(object)) ||
      anyDuplicated(names(object)) || !setequal(names(object), expected)) {
    cdid_ipw_point_view_abort(
      "Point-link snapshot fields differ from the exact compact schema.",
      c("cdid_ipw_point_view_revalidation_error", "cdid_ipw_owner_integrity_error")
    )
  }

  validate_ipw_fitted_nuisance_view(fitted_view)
  validate_ipw_link_view(link_view)
  ctrl <- validate_ipw_point_view_control(control)
  glm_ctrl <- if (is.null(glm_control)) fitted_view$glm_control else validate_ipw_glm_control(glm_control)

  if (!ipw_point_view_identical(object$glm_control, glm_ctrl)) {
    cdid_ipw_point_view_abort(
      "Deserialized point-link snapshot `glm_control` disagrees with caller `glm_control`.",
      c("cdid_ipw_point_view_revalidation_error", "cdid_ipw_control_mismatch_error")
    )
  }
  if (!ipw_point_view_identical(object$control, ctrl)) {
    cdid_ipw_point_view_abort(
      "Point-link snapshot resource controls differ from explicit revalidation controls.",
      c("cdid_ipw_point_view_revalidation_error", "cdid_ipw_control_mismatch_error")
    )
  }

  pre <- preflight_ipw_point_link_view(
    fitted_view = fitted_view,
    link_view = link_view,
    control = ctrl,
    glm_control = glm_ctrl,
    n_planned_links = 1L
  )
  ipw_point_link_snapshot_payload(
    object,
    max_entries = max(0, ctrl$max_matrix_entries - pre$peak_concurrent_entries),
    max_work = max(0, ctrl$max_work_ops - pre$total_work_ops_bound - pre$character_metadata$work_ops_bound - 128)
  )
  ipw_point_view_verify_compact_storage(object, p_max = pre$p_covariates, m_max = pre$m_cohorts)

  live_link <- fit_ipw_point_link_view(
    fitted_view = fitted_view,
    link_view = link_view,
    control = ctrl,
    glm_control = glm_ctrl,
    link_id = object$link_id,
    functional_key = object$functional_key
  )

  for (nm in snapshot_fields) {
    if (!ipw_point_view_identical(object[[nm]], live_link[[nm]])) {
      cdid_ipw_point_view_abort(
        paste0("Deserialized point-link field `", nm, "` differs from live owner-backed recomputation; revalidation is refused."),
        c("cdid_ipw_point_view_revalidation_error", "cdid_ipw_owner_integrity_error"),
        diagnostics = list(field = nm)
      )
    }
  }

  live_link
}

# ---------------------------------------------------------------------------
# Section 2: Bounded Complete Point IF Replay and Covariance Callback Producer
# ---------------------------------------------------------------------------

normalize_ipw_point_links_view <- function(point_links, frame_owner = NULL, control = NULL) {
  ipw_point_view_require_deps()
  if (missing(point_links) || is.null(point_links)) {
    cdid_ipw_point_view_abort(
      "`point_links` must be a sealed `cdid_ipw_point_link_view` or a non-empty list of them.",
      "cdid_ipw_point_view_input_error"
    )
  }
  single <- inherits(point_links, "cdid_ipw_point_link_view")
  if (!single) {
    ipw_owner_refuse_decorative_attrs(point_links, "point_links")
  }
  link_list <- if (single) list(point_links) else point_links
  if (!is.list(link_list) || length(link_list) < 1L) {
    cdid_ipw_point_view_abort(
      "`point_links` must be a sealed `cdid_ipw_point_link_view` or a non-empty list of them.",
      "cdid_ipw_point_view_input_error"
    )
  }
  L <- length(link_list)
  if (!is.null(control)) {
    ipw_owner_safe_product(
      64, L,
      min(control$max_retained_entries, control$max_matrix_entries, control$max_output_entries),
      "point-link list header before link validation"
    )
    ipw_owner_safe_product(
      128, L,
      control$max_work_ops,
      "point-link list header validation work"
    )
  }

  if (is.null(control)) {
    cdid_ipw_point_view_abort("Explicit control required for point-link normalization.", "cdid_ipw_point_view_control_error")
  }
  # Quadratic worst case: all L nuisance references distinct. Charge BEFORE search.
  search_work <- ipw_owner_safe_product(16 * L, L, control$max_work_ops,
                                        "distinct-reference search and link-ID set construction")
  for (lk in link_list) ipw_point_view_seal_gate(lk, "point_link_view")
  payload <- ipw_point_view_payload_budget(link_list, control)
  normalization_work <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(search_work, 1 + payload$max_bytes, control$max_work_ops,
      "worst-case link-ID set byte comparisons"), 8 * payload$work),
    control$max_work_ops, "all normalization validation and metadata work")
  list_names <- names(link_list)
  if (!is.null(list_names) && (anyNA(list_names) || any(!nzchar(list_names)) || anyDuplicated(list_names))) {
    cdid_ipw_point_view_abort(
      "If `point_links` is a named list, its names must be unique non-empty character link IDs.",
      "cdid_ipw_point_view_input_error"
    )
  }

  validate_ipw_point_link_view(link_list[[1L]], frame_owner = frame_owner)
  fref_1 <- link_list[[1L]]$frame_ref
  pref_1 <- link_list[[1L]]$panel_ref
  f_obj <- get("frame_owner", envir = fref_1, inherits = FALSE)
  p_obj <- get("panel_owner", envir = pref_1, inherits = FALSE)

  link_ids <- character(L)
  seen_fitted_refs <- list()
  seen_nv_refs <- list()
  existing_retained <- 0

  for (ell in seq_len(L)) {
    lk <- link_list[[ell]]
    if (ell > 1L) {
      validate_ipw_point_link_view(lk)
      if (!ipw_point_view_identical(lk$frame_ref, fref_1) || !ipw_point_view_identical(lk$panel_ref, pref_1)) {
        cdid_ipw_point_view_abort(
          "All point links in replay must share the exact same sealed `frame_owner` and `panel_owner`.",
          c("cdid_ipw_owner_mismatch_error", "cdid_ipw_context_mismatch_error")
        )
      }
    }
    lid <- lk$link_id
    if (!is.null(lid) && nzchar(lid)) {
      if (!is.null(list_names) && !ipw_point_view_identical(list_names[ell], lid)) {
        cdid_ipw_point_view_abort(
          paste0("`names(point_links)[", ell, "]` ('", list_names[ell], "') disagrees with `point_links[[", ell, "]]$link_id` ('", lid, "')."),
          "cdid_ipw_point_view_input_error"
        )
      }
      link_ids[ell] <- lid
    } else if (!is.null(list_names)) {
      link_ids[ell] <- list_names[ell]
    } else {
      link_ids[ell] <- paste0("link_", ell)
    }

    already_fv <- any(vapply(seen_fitted_refs, identical, logical(1L), y = lk$fitted_nuisance_ref))
    if (!already_fv) {
      seen_fitted_refs[[length(seen_fitted_refs) + 1L]] <- lk$fitted_nuisance_ref
      fv_ell <- get("fitted_nuisance_view", envir = lk$fitted_nuisance_ref, inherits = FALSE)
      existing_retained <- existing_retained + as.double(fv_ell$preflight$retained_entries_per_fit)
    }
    already_nv <- any(vapply(seen_nv_refs, identical, logical(1L), y = lk$nuisance_view_ref))
    if (!already_nv) {
      seen_nv_refs[[length(seen_nv_refs) + 1L]] <- lk$nuisance_view_ref
    }
    existing_retained <- existing_retained + as.double(lk$preflight$retained_entries_per_link)
  }

  if (anyDuplicated(link_ids)) {
    cdid_ipw_point_view_abort(
      "Point link IDs across `point_links` must be unique.",
      "cdid_ipw_point_view_input_error",
      diagnostics = list(link_ids = link_ids)
    )
  }

  list(
    links = link_list,
    link_ids = link_ids,
    single = single,
    n_links = L,
    n_unique_fitted_nuisances = length(seen_fitted_refs),
    frame_owner = f_obj,
    panel_owner = p_obj,
    existing_retained_entries = max(existing_retained, payload$entries) +
      ipw_point_view_owner_entries(f_obj, p_obj),
    borrowed_owner_entries = ipw_point_view_owner_entries(f_obj, p_obj),
    normalization_work_bound = normalization_work,
    normalization_payload_entries = payload$entries,
    validated = NULL
  )
}


# The caller has completed the whole replay allocation/work preflight before
# this helper reads any posting. Records reference existing owners/views;
# 64 slots per link cover the record, seal, two date facts and bounded keys.
# 256 slots cover the shared one-entry control/range-preflight memo.
ipw_point_view_add_query_records <- function(norm, control, step) {
  query_memo <- new.env(parent = emptyenv())
  norm$validated <- lapply(norm$links, function(lk) ipw_owner_validated_view_record(
    get("link_view", envir = lk$link_view_ref, inherits = FALSE), norm$frame_owner,
    norm$panel_owner, query_memo, step = step))
  norm
}

preflight_ipw_point_if_view <- function(
  point_links,
  control,
  include_components = FALSE,
  chunk_rows = NULL,
  return_norm = FALSE
) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`point_links`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  ipw_owner_refuse_decorative_attrs(control, "control")
  ipw_owner_refuse_decorative_attrs(include_components, "include_components")
  if (!gmm_scalar_logical(include_components)) {
    cdid_ipw_point_view_abort("`include_components` must be a scalar logical (`TRUE` or `FALSE`).", "cdid_ipw_point_view_input_error")
  }
  ctrl <- validate_ipw_point_view_control(control)
  norm <- normalize_ipw_point_links_view(point_links, control = ctrl)

  frame_owner <- norm$frame_owner
  panel_owner <- norm$panel_owner
  na <- ipw_owner_int_bound(frame_owner$n_a_units, "n_a_units", min_val = 1)
  p <- ipw_owner_int_bound(frame_owner$n_covariates, "n_covariates", min_val = 1)
  L <- norm$n_links

  max_ch <- ctrl$max_chunk_rows
  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  req_ch <- if (is.null(chunk_rows)) max_ch else ipw_owner_int_bound(chunk_rows, "chunk_rows", min_val = 1)
  if (req_ch > max_ch) {
    cdid_ipw_point_view_abort(
      paste0("Requested chunk rows (", req_ch, ") exceed `control$max_chunk_rows` (", max_ch, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error")
    )
  }
  B_eff <- min(na, req_ch)
  n_chunks <- ceiling(na / B_eff)

  # Constant header allowance before character metadata scans
  ipw_owner_safe_sum(
    c(ipw_owner_safe_product(16, L + p + 1), 64),
    min(max_ret, max_out, max_me),
    "point-IF metadata header before character scans"
  )
  ipw_owner_safe_sum(
    c(ipw_owner_safe_product(32, L + p + 1), 128),
    max_wo,
    "point-IF header work before character scans"
  )

  char_meta <- ipw_point_view_char_bytes(
    c(list(norm$link_ids, colnames(frame_owner$Z_a)),
      lapply(norm$links, function(lk) lk$cohorts$cohort_labels)),
    ctrl
  )
  id_meta <- ipw_point_view_id_bytes_bound(panel_owner, B_eff, n_chunks, ctrl)

  # 1. Chunk output matrices: B_eff x L phi (plus 4 B_eff x L component matrices when include_components = TRUE)
  n_out_matrices <- if (isTRUE(include_components)) 5L else 1L
  b_by_L <- ipw_owner_check_matrix_dims(
    B_eff, L,
    max_matrix_entries = min(max_me, max_out),
    label = "point-IF chunk B x L matrix"
  )
  chunk_output_matrix_entries <- ipw_owner_safe_product(
    n_out_matrices, b_by_L,
    max_entries = min(max_me, max_out),
    label = "point-IF chunk output matrices"
  )
  chunk_output_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(6, B_eff, min(max_me, max_out), "chunk index/ID/cohort/weight vectors"),
      L,
      char_meta$charged_slots,
      id_meta$chunk_id_slots,
      64
    ),
    max_bound = min(max_me, max_out),
    label = "point-IF chunk output vectors and metadata"
  )
  chunk_output_entries <- ipw_owner_safe_sum(
    c(chunk_output_matrix_entries, chunk_output_vector_entries),
    max_bound = min(max_me, max_out),
    label = "point-IF chunk total output entries"
  )

  # 2. Simultaneous reader chunk, local per-link columns, and bounded one-nuisance chunk cache
  owner_ctrl <- ipw_point_view_to_owner_control(ctrl)
  max_m <- max(vapply(norm$links, function(lk) as.integer(lk$cohorts$m_cohorts), integer(1L)))
  chunk_pre <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p,
    m_cohorts = max_m,
    control = owner_ctrl,
    read_outcomes = TRUE,
    include_design = TRUE,
    enclose_differences = FALSE
  )
  reader_chunk_concurrent_entries <- chunk_pre$concurrent_entries

  obs_Z_chunk_entries <- ipw_owner_check_matrix_dims(
    B_eff, p,
    max_matrix_entries = max_me,
    label = "point-IF chunk observer design Z_k"
  )
  local_link_working_entries <- ipw_owner_safe_sum(
    c(
      obs_Z_chunk_entries,
      ipw_owner_safe_product(64, B_eff, max_me, "columns, conversion copies, subsets, expression RHS and result columns"),
      ipw_owner_safe_product(4, p + max_m, max_me, "per-link coefficient/loading slices"),
      64
    ),
    max_bound = max_me,
    label = "point-IF local link working entries"
  )
  one_nuisance_cache_entries <- ipw_owner_safe_sum(
    c(
      obs_Z_chunk_entries,
      ipw_owner_safe_product(12, B_eff, max_me, "one-nuisance chunk cached observer positions, weights, rates, odds, prob, and score base"),
      64
    ),
    max_bound = max_me,
    label = "bounded one-nuisance chunk cache entries"
  )
  failure_output_entries <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(4, B_eff), L, 128),
    min(max_ret, max_out, max_me),
    "possible failed point-IF chunk diagnostic output"
  )

  existing_retained <- ipw_owner_safe_sum(
    c(norm$existing_retained_entries, 64 * L + 256),
    max_bound = min(max_ret, max_me),
    label = "existing retained fitted-nuisance and point-link entries"
  )

  peak_concurrent_entries <- ipw_owner_safe_sum(
    c(
      existing_retained,
      2 * chunk_output_entries, # allocation/column replacement and return wrapper copies
      reader_chunk_concurrent_entries,
      local_link_working_entries,
      one_nuisance_cache_entries,
      failure_output_entries
    ),
    max_bound = max_me,
    label = "simultaneous point-IF replay concurrent entries"
  )
  peak_concurrent_matrix_entries <- ipw_owner_safe_sum(
    c(
      chunk_output_matrix_entries,
      chunk_pre$output_matrix_entries,
      2 * obs_Z_chunk_entries
    ),
    max_bound = max_me,
    label = "simultaneous point-IF replay matrix entries"
  )

  # 3. Whole-traversal work operations across all chunks and all L links (NEVER per-chunk only)
  per_link_stream_work <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, chunk_pre$work_ops_bound, max_wo, "whole-stream reader chunk queries per link"),
      ipw_owner_safe_product(na, 8 * p + 24, max_wo, "whole-stream 4-term point IF evaluation per link"),
      ipw_owner_safe_product(n_chunks, 32, max_wo, "per-chunk dispatch and alias release work per link")
    ),
    max_bound = max_wo,
    label = "single-link whole-stream point-IF work"
  )
  all_links_stream_work <- ipw_owner_safe_product(
    L, per_link_stream_work,
    max_entries = max_wo,
    label = "all-links whole-stream point-IF work"
  )
  total_work_ops_bound <- ipw_owner_safe_sum(
    c(
      all_links_stream_work,
      char_meta$work_ops_bound,
      id_meta$stream_id_work_bound,
      4 * norm$normalization_work_bound, # preparation, validation and callback dispatch setup
      256
    ),
    max_bound = max_wo,
    label = "total point-IF replay work operations bound"
  )

  total_expected_endpoints_read <- sum(vapply(
    norm$links,
    function(lk) 2 * as.double(lk$counts$n_observed_risk_units),
    numeric(1L)
  ))

  pre <- list(
    n_a_units = na,
    n_links = L,
    link_ids = norm$link_ids,
    p_covariates = p,
    m_cohorts_max = max_m,
    n_unique_fitted_nuisances = norm$n_unique_fitted_nuisances,
    include_components = isTRUE(include_components),
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks = logit_cert_count_output(n_chunks),
    n_output_matrices_per_chunk = n_out_matrices,
    b_by_L_entries = b_by_L,
    chunk_output_matrix_entries = chunk_output_matrix_entries,
    chunk_output_vector_entries = chunk_output_vector_entries,
    chunk_output_entries = chunk_output_entries,
    reader_chunk_concurrent_entries = reader_chunk_concurrent_entries,
    local_link_working_entries = local_link_working_entries,
    one_nuisance_cache_entries = one_nuisance_cache_entries,
    existing_retained_entries = existing_retained,
    borrowed_owner_entries = norm$borrowed_owner_entries,
    failure_output_entries = failure_output_entries,
    character_metadata = char_meta,
    id_byte_metadata = id_meta,
    peak_concurrent_entries = peak_concurrent_entries,
    peak_concurrent_matrix_entries = peak_concurrent_matrix_entries,
    peak_concurrent_bytes_8double = 8 * peak_concurrent_entries,
    per_link_stream_work_bound = per_link_stream_work,
    total_work_ops_bound = total_work_ops_bound,
    normalization_work_bound = norm$normalization_work_bound,
    normalization_payload_entries = norm$normalization_payload_entries,
    producer_workspace_entries = peak_concurrent_entries - existing_retained,
    consumer_retention_included = FALSE,
    total_expected_endpoints_read = logit_cert_count_output(total_expected_endpoints_read),
    full_n_by_L_retained = FALSE,
    per_link_dense_masks_retained = FALSE,
    peak_rss_promised = FALSE
  )
  # Replay and covariance callers reuse this exact normalization (the same calls
  # on the same immutable inputs) instead of repeating it.
  if (isTRUE(return_norm)) {
    norm <- ipw_point_view_add_query_records(norm, ctrl, B_eff)
    list(preflight = pre, norm = norm, ctrl = ctrl)
  } else pre
}

ipw_point_if_eval_chunk_internal <- function(
  norm,
  a_start,
  a_end,
  ctrl,
  include_components = FALSE,
  callback_schema = FALSE,
  use_nuisance_cache = TRUE,
  batch_index = 1L,
  trace_env = NULL,
  state = NULL
) {
  if (!gmm_scalar_logical(include_components) || !gmm_scalar_logical(callback_schema) ||
      !gmm_scalar_logical(use_nuisance_cache) ||
      (callback_schema && include_components)) {
    cdid_ipw_point_view_abort("Invalid point IF output options.", "cdid_ipw_point_view_input_error")
  }
  batch_index <- ipw_owner_int_bound(batch_index, "batch_index", min_val = 1)
  owner_ctrl <- ipw_point_view_to_owner_control(ctrl)
  frame_owner <- norm$frame_owner
  panel_owner <- norm$panel_owner
  na <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  L <- norm$n_links
  link_ids <- norm$link_ids

  s_idx <- ipw_owner_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_owner_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > na) {
    cdid_ipw_point_view_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= n_a_units` (", na, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error")
    )
  }
  a_s <- as.integer(s_idx)
  a_e <- as.integer(e_idx)
  b_len <- a_e - a_s + 1L
  if (b_len > ctrl$max_chunk_rows) {
    cdid_ipw_point_view_abort(
      paste0("Requested chunk length (", b_len, ") exceeds `control$max_chunk_rows` (", ctrl$max_chunk_rows, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error")
    )
  }

  if (!is.null(state)) {
    if (!is.list(state) || !identical(state$kind, "cdid_ipw_point_if_general_state") ||
        !ipw_owner_seal_ok(state, "point_if_general_state") ||
        !is.environment(state$runtime) || !identical(state$runtime$status, "built") ||
        !identical(state$frame_ref, frame_owner$ref) ||
        !identical(state$panel_ref, frame_owner$panel_ref) ||
        !identical(state$n_a_units, na) ||
        !identical(state$n_links, L) ||
        !identical(state$link_ids, link_ids) ||
        !is.list(state$link_states) || length(state$link_states) != L ||
        !is.list(state$link_state_refs) || length(state$link_state_refs) != L) {
      ipw_point_view_general_abort("General point-IF replay state is malformed, foreign, released, or differs from its sealed snapshot.")
    }
    if (batch_index < 1L || batch_index > state$n_chunks) {
      ipw_point_view_general_abort("General point-IF replay state chunk index is out of bounds.")
    }
    exp_s <- as.integer(1 + (as.double(batch_index) - 1) * as.double(state$step))
    exp_e <- as.integer(min(as.double(na), as.double(exp_s) + as.double(state$step) - 1))
    if (!identical(a_s, exp_s) || !identical(a_e, exp_e)) {
      ipw_point_view_general_abort("General point-IF replay state chunk range disagrees with the sealed chunk partition.")
    }
  }

  a_seq <- seq.int(a_s, a_e)
  u_idx <- panel_owner$a_unit_index[a_seq]
  u_id <- panel_owner$a_id[a_seq]
  u_g <- panel_owner$a_g[a_seq]
  u_w <- frame_owner$w_a[a_seq]

  phi_mat <- matrix(0.0, nrow = b_len, ncol = L, dimnames = list(NULL, link_ids))
  comp_rate <- if (include_components) matrix(0.0, nrow = b_len, ncol = L, dimnames = list(NULL, link_ids)) else NULL
  comp_score <- if (include_components) matrix(0.0, nrow = b_len, ncol = L, dimnames = list(NULL, link_ids)) else NULL
  comp_treated <- if (include_components) matrix(0.0, nrow = b_len, ncol = L, dimnames = list(NULL, link_ids)) else NULL
  comp_control <- if (include_components) matrix(0.0, nrow = b_len, ncol = L, dimnames = list(NULL, link_ids)) else NULL

  nu_cache_env <- new.env(parent = emptyenv())
  nu_cache_env$active <- FALSE
  clear_nuisance_cache <- function() {
    for (nm in ls(nu_cache_env, all.names = TRUE)) {
      if (!ipw_point_view_identical(nm, "active")) {
        assign(nm, NULL, envir = nu_cache_env)
        rm(list = nm, envir = nu_cache_env)
      }
    }
    nu_cache_env$active <- FALSE
  }
  on.exit(clear_nuisance_cache(), add = TRUE)

  active_query_chunks <- 0L
  max_simultaneous_query_chunks <- 0L
  max_simultaneous_nuisance_caches <- 0L
  nuisance_cache_hits <- 0L
  endpoints_read_chunk <- 0L
  actual_ops <- 0
  y_trace <- new.env(parent = emptyenv())
  on.exit(ipw_owner_reset_trace(y_trace), add = TRUE)

  if (!is.null(state)) {
    u_a_cohort <- panel_owner$a_cohort_index[a_seq]
    Z_a <- frame_owner$Z_a
    prev_cohorts <- NULL
    gi_chunk <- NULL
    in_risk <- NULL
    in_risk_chunk <- NULL

    for (ell in seq_len(L)) {
      ipw_owner_check_interrupt()
      fit <- norm$links[[ell]]
      lv_ell <- get("link_view", envir = fit$link_view_ref, inherits = FALSE)

      can_reuse_cache <- isTRUE(use_nuisance_cache) &&
        isTRUE(nu_cache_env$active) &&
        ipw_point_view_identical(nu_cache_env$fitted_nuisance_ref, fit$fitted_nuisance_ref) &&
        ipw_point_view_identical(nu_cache_env$a_start, a_s) &&
        ipw_point_view_identical(nu_cache_env$a_end, a_e)

      if (!can_reuse_cache && isTRUE(nu_cache_env$active)) {
        clear_nuisance_cache()
      }

      sl <- ipw_point_if_general_state_link_slice(
        state = state,
        ell = ell,
        batch_index = batch_index,
        a_s = a_s,
        a_e = a_e,
        b_len = b_len,
        fit = fit,
        lv_ell = lv_ell
      )
      k_ch <- sl$k_ch
      a_k <- sl$a_index
      gj_k <- sl$cohort_index
      obs_s_k <- sl$obs_s_index
      obs_t_k <- sl$obs_t_index
      odds_k <- sl$odds
      score_base_k <- sl$score_base
      sl <- NULL

      if (is.null(prev_cohorts) || !identical(lv_ell$cohorts, prev_cohorts)) {
        cohort_index <- ipw_owner_view_cohort_index(lv_ell$cohorts, u_a_cohort)
        in_risk <- !is.na(cohort_index)
        gi_chunk <- as.integer(cohort_index)
        in_risk_chunk <- which(in_risk)
        prev_cohorts <- lv_ell$cohorts
      }

      # Indices, risk membership, endpoint positions and score finiteness were
      # checked once during construction and are protected by the locked COW
      # link seal checked by the slice helper. Do not rescan immutable payloads.
      pos_k <- if (k_ch > 0L) a_k - a_s + 1L else integer(0L)

      ipw_owner_reset_trace(y_trace)
      ep <- ipw_owner_read_y_positions(
        panel_owner = panel_owner,
        pos_s = obs_s_k,
        pos_t = obs_t_k,
        trace_env = y_trace
      )
      ep_ch <- as.integer(2L * k_ch)
      endpoints_read_chunk <- endpoints_read_chunk + ep_ch

      if (!ipw_point_view_identical(ep_ch, as.integer(2L * k_ch)) ||
          !ipw_point_view_identical(length(y_trace$pos_s_index), k_ch) ||
          !ipw_point_view_identical(length(y_trace$pos_t_index), k_ch)) {
        clear_nuisance_cache()
        cdid_ipw_point_view_abort(
          "Point IF chunk query violated risk-first endpoint read filtering.",
          c("cdid_ipw_point_view_outcome_leak_error", "cdid_ipw_owner_integrity_error")
        )
      }
      ipw_owner_reset_trace(y_trace)

      col_rate <- numeric(b_len)
      col_score <- numeric(b_len)
      col_treated <- numeric(b_len)
      col_control <- numeric(b_len)

      if (length(in_risk_chunk) > 0L) {
        gj_all <- gi_chunk[in_risk_chunk]
        col_rate[in_risk_chunk] <- u_w[in_risk_chunk] * fit$loadings$rate_coef_s0[gj_all]
        rm(gj_all)

        if (k_ch > 0L) {
          yl_k <- as.double(ep$y_left_obs)
          yr_k <- as.double(ep$y_right_obs)
          dy_k <- yr_k - yl_k
          if (anyNA(dy_k) || any(!is.finite(dy_k))) {
            clear_nuisance_cache()
            cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
          }

          if (can_reuse_cache) {
            if (!ipw_point_view_identical(pos_k, nu_cache_env$pos_k) ||
                !ipw_point_view_identical(a_k, nu_cache_env$a_k) ||
                !ipw_point_view_identical(gj_k, nu_cache_env$gj_k)) {
              clear_nuisance_cache()
              cdid_ipw_point_view_abort(
                "Shared-nuisance point link observer mask disagrees with cached nuisance chunk mask.",
                c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error")
              )
            }
            nuisance_cache_hits <- nuisance_cache_hits + 1L
            w_k <- nu_cache_env$w_k
            r_k <- nu_cache_env$r_k
            Z_k <- nu_cache_env$Z_k
          } else {
            w_k <- u_w[pos_k]
            r_k <- ipw_clip_floor(fit$masses$rates[gj_k])
            Z_k <- Z_a[a_k, , drop = FALSE]
            if (isTRUE(use_nuisance_cache)) {
              nu_cache_env$active <- TRUE
              nu_cache_env$fitted_nuisance_ref <- fit$fitted_nuisance_ref
              nu_cache_env$a_start <- a_s
              nu_cache_env$a_end <- a_e
              nu_cache_env$k_ch <- k_ch
              nu_cache_env$pos_k <- pos_k
              nu_cache_env$a_k <- a_k
              nu_cache_env$gj_k <- gj_k
              nu_cache_env$w_k <- w_k
              nu_cache_env$r_k <- r_k
              nu_cache_env$Z_k <- Z_k
              max_simultaneous_nuisance_caches <- max(max_simultaneous_nuisance_caches, 1L)
            }
          }

          col_rate[pos_k] <- w_k * fit$loadings$rate_coef_s1[gj_k]
          z_load_k <- drop(Z_k %*% fit$loadings$loading_beta)
          col_score[pos_k] <- score_base_k * z_load_k

          is_tk <- gj_k == 1L
          if (any(is_tk)) {
            pos_t <- pos_k[is_tk]
            col_treated[pos_t] <- (w_k[is_tk] / r_k[is_tk]) / fit$means$mean_a *
              (dy_k[is_tk] - fit$means$treated_mean)
            rm(pos_t)
          }
          is_ck <- !is_tk
          if (any(is_ck)) {
            pos_c <- pos_k[is_ck]
            col_control[pos_c] <- -(w_k[is_ck] * odds_k[is_ck] / r_k[is_ck]) /
              fit$means$mean_b * (dy_k[is_ck] - fit$means$control_mean)
            rm(pos_c)
          }
          rm(yl_k, yr_k, dy_k, w_k, r_k, Z_k, z_load_k, is_tk, is_ck)
        } else if (can_reuse_cache) {
          if (!ipw_point_view_identical(nu_cache_env$k_ch, 0L)) {
            clear_nuisance_cache()
            cdid_ipw_point_view_abort(
              "Shared-nuisance point link observer count disagrees with cached nuisance chunk.",
              c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error")
            )
          }
          nuisance_cache_hits <- nuisance_cache_hits + 1L
        } else if (isTRUE(use_nuisance_cache)) {
          nu_cache_env$active <- TRUE
          nu_cache_env$fitted_nuisance_ref <- fit$fitted_nuisance_ref
          nu_cache_env$a_start <- a_s
          nu_cache_env$a_end <- a_e
          nu_cache_env$k_ch <- 0L
          max_simultaneous_nuisance_caches <- max(max_simultaneous_nuisance_caches, 1L)
        }
      }

      col_phi <- col_treated + col_control + col_score + col_rate
      if (anyNA(col_phi) || any(!is.finite(col_phi))) {
        clear_nuisance_cache()
        cdid_abort("Nonfinite IF contribution encountered during replay.", c("cdid_ipw_numeric_error"))
      }
      phi_mat[, ell] <- col_phi
      if (include_components) {
        comp_rate[, ell] <- col_rate
        comp_score[, ell] <- col_score
        comp_treated[, ell] <- col_treated
        comp_control[, ell] <- col_control
      }
      actual_ops <- actual_ops + 2 * k_ch + k_ch * p + b_len * (4 * p + 12) + 6 * k_ch + 64
      rm(fit, lv_ell, a_k, gj_k, obs_s_k, obs_t_k, odds_k, score_base_k,
         pos_k, ep, col_rate, col_score, col_treated, col_control, col_phi)
    }
    actual_ops <- actual_ops + 256
    state$runtime$chunks_evaluated <- logit_cert_count_output(as.double(state$runtime$chunks_evaluated) + 1)
  } else {
    for (ell in seq_len(L)) {
    ipw_owner_check_interrupt()
    fit <- norm$links[[ell]]
    lv_ell <- get("link_view", envir = fit$link_view_ref, inherits = FALSE)

    can_reuse_cache <- isTRUE(use_nuisance_cache) &&
      isTRUE(nu_cache_env$active) &&
      ipw_point_view_identical(nu_cache_env$fitted_nuisance_ref, fit$fitted_nuisance_ref) &&
      ipw_point_view_identical(nu_cache_env$a_start, a_s) &&
      ipw_point_view_identical(nu_cache_env$a_end, a_e)

    if (!can_reuse_cache && isTRUE(nu_cache_env$active)) {
      clear_nuisance_cache()
    }

    ch <- query_ipw_view_chunk(
      view = lv_ell,
      a_start = a_s,
      a_end = a_e,
      control = owner_ctrl,
      read_outcomes = TRUE,
      include_design = !can_reuse_cache,
      enclose_differences = FALSE,
      trace_env = y_trace,
      validated = norm$validated[[ell]]
    )
    active_query_chunks <- active_query_chunks + 1L
    max_simultaneous_query_chunks <- max(max_simultaneous_query_chunks, active_query_chunks)

    k_ch <- as.integer(ch$counts$observed_risk_row_count)
    ep_ch <- as.integer(ch$read_diagnostics$n_endpoints_read)
    endpoints_read_chunk <- endpoints_read_chunk + ep_ch
    actual_ops <- actual_ops + ch$work$actual_work_ops

    if (!isTRUE(ch$read_diagnostics$risk_filtered_before_outcome_read) ||
        isTRUE(ch$read_diagnostics$off_risk_endpoints_read) ||
        isTRUE(ch$read_diagnostics$nonobserver_endpoints_read) ||
        !ipw_point_view_identical(ep_ch, as.integer(2L * k_ch)) ||
        !ipw_point_view_identical(length(y_trace$pos_s_index), k_ch) ||
        !ipw_point_view_identical(length(y_trace$pos_t_index), k_ch)) {
      ch <- NULL
      clear_nuisance_cache()
      cdid_ipw_point_view_abort(
        "Point IF chunk query violated risk-first endpoint read filtering.",
        c("cdid_ipw_point_view_outcome_leak_error", "cdid_ipw_owner_integrity_error")
      )
    }

    ipw_owner_reset_trace(y_trace)
    gi_chunk <- as.integer(ch$cohort_index)
    in_risk_chunk <- which(ch$in_risk)

    col_rate <- numeric(b_len)
    col_score <- numeric(b_len)
    col_treated <- numeric(b_len)
    col_control <- numeric(b_len)

    if (length(in_risk_chunk) > 0L) {
      gj_all <- gi_chunk[in_risk_chunk]
      col_rate[in_risk_chunk] <- u_w[in_risk_chunk] * fit$loadings$rate_coef_s0[gj_all]
      rm(gj_all)

      if (k_ch > 0L) {
        pos_k <- as.integer(ch$observed_risk$local_pos)
        a_k <- as.integer(ch$observed_risk$a_index)
        gj_k <- as.integer(ch$observed_risk$cohort_index)
        yl_k <- as.double(ch$observed_risk$y_left_obs)
        yr_k <- as.double(ch$observed_risk$y_right_obs)
        dy_k <- yr_k - yl_k
        if (anyNA(dy_k) || any(!is.finite(dy_k))) {
          rm(pos_k, a_k, gj_k, yl_k, yr_k, dy_k, gi_chunk, in_risk_chunk,
             col_rate, col_score, col_treated, col_control)
          ch <- NULL
          clear_nuisance_cache()
          cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
        }

        if (can_reuse_cache) {
          if (!ipw_point_view_identical(pos_k, nu_cache_env$pos_k) ||
              !ipw_point_view_identical(a_k, nu_cache_env$a_k) ||
              !ipw_point_view_identical(gj_k, nu_cache_env$gj_k)) {
            rm(pos_k, a_k, gj_k, yl_k, yr_k, dy_k, gi_chunk, in_risk_chunk,
               col_rate, col_score, col_treated, col_control)
            ch <- NULL
            clear_nuisance_cache()
            cdid_ipw_point_view_abort(
              "Shared-nuisance point link observer mask disagrees with cached nuisance chunk mask.",
              c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error")
            )
          }
          nuisance_cache_hits <- nuisance_cache_hits + 1L
          w_k <- nu_cache_env$w_k
          r_k <- nu_cache_env$r_k
          Z_k <- nu_cache_env$Z_k
          odds_k <- nu_cache_env$odds_k
          score_base_k <- nu_cache_env$score_base_k
        } else {
          w_k <- u_w[pos_k]
          r_k <- ipw_clip_floor(fit$masses$rates[gj_k])
          D_k <- as.double(gj_k == 1L)
          Z_k <- ch$design$Z_obs
          lv_k <- ipw_logit_values(Z_k, fit$coefficients)
          odds_k <- lv_k$odds
          prob_k <- lv_k$prob
          score_base_k <- (w_k / r_k) * (D_k - prob_k)
          rm(lv_k)
          if (isTRUE(use_nuisance_cache)) {
            nu_cache_env$active <- TRUE
            nu_cache_env$fitted_nuisance_ref <- fit$fitted_nuisance_ref
            nu_cache_env$a_start <- a_s
            nu_cache_env$a_end <- a_e
            nu_cache_env$k_ch <- k_ch
            nu_cache_env$pos_k <- pos_k
            nu_cache_env$a_k <- a_k
            nu_cache_env$gj_k <- gj_k
            nu_cache_env$w_k <- w_k
            nu_cache_env$r_k <- r_k
            nu_cache_env$Z_k <- Z_k
            nu_cache_env$odds_k <- odds_k
            nu_cache_env$score_base_k <- score_base_k
            max_simultaneous_nuisance_caches <- max(max_simultaneous_nuisance_caches, 1L)
          }
          rm(D_k, prob_k)
        }

        col_rate[pos_k] <- w_k * fit$loadings$rate_coef_s1[gj_k]
        z_load_k <- drop(Z_k %*% fit$loadings$loading_beta)
        col_score[pos_k] <- score_base_k * z_load_k

        is_tk <- gj_k == 1L
        if (any(is_tk)) {
          pos_t <- pos_k[is_tk]
          col_treated[pos_t] <- (w_k[is_tk] / r_k[is_tk]) / fit$means$mean_a *
            (dy_k[is_tk] - fit$means$treated_mean)
          rm(pos_t)
        }
        is_ck <- !is_tk
        if (any(is_ck)) {
          pos_c <- pos_k[is_ck]
          col_control[pos_c] <- -(w_k[is_ck] * odds_k[is_ck] / r_k[is_ck]) /
            fit$means$mean_b * (dy_k[is_ck] - fit$means$control_mean)
          rm(pos_c)
        }

        rm(pos_k, a_k, gj_k, yl_k, yr_k, dy_k, w_k, r_k, Z_k, odds_k, score_base_k, z_load_k, is_tk, is_ck)
      } else if (can_reuse_cache) {
        if (!ipw_point_view_identical(nu_cache_env$k_ch, 0L)) {
          rm(gi_chunk, in_risk_chunk, col_rate, col_score, col_treated, col_control)
          ch <- NULL
          clear_nuisance_cache()
          cdid_ipw_point_view_abort(
            "Shared-nuisance point link observer count disagrees with cached nuisance chunk.",
            c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error")
          )
        }
        nuisance_cache_hits <- nuisance_cache_hits + 1L
      } else if (isTRUE(use_nuisance_cache)) {
        nu_cache_env$active <- TRUE
        nu_cache_env$fitted_nuisance_ref <- fit$fitted_nuisance_ref
        nu_cache_env$a_start <- a_s
        nu_cache_env$a_end <- a_e
        nu_cache_env$k_ch <- 0L
        max_simultaneous_nuisance_caches <- max(max_simultaneous_nuisance_caches, 1L)
      }
    }

    # Release consumed query chunk BEFORE assembling column or moving to next link!
    ch <- NULL
    active_query_chunks <- active_query_chunks - 1L

    col_phi <- col_treated + col_control + col_score + col_rate
    if (anyNA(col_phi) || any(!is.finite(col_phi))) {
      rm(gi_chunk, in_risk_chunk, col_rate, col_score, col_treated, col_control, col_phi)
      clear_nuisance_cache()
      cdid_abort("Nonfinite IF contribution encountered during replay.", c("cdid_ipw_numeric_error"))
    }
    phi_mat[, ell] <- col_phi
    if (include_components) {
      comp_rate[, ell] <- col_rate
      comp_score[, ell] <- col_score
      comp_treated[, ell] <- col_treated
      comp_control[, ell] <- col_control
    }
    actual_ops <- actual_ops + b_len * (4 * p + 12)
    rm(fit, lv_ell, gi_chunk, in_risk_chunk, col_rate, col_score, col_treated, col_control, col_phi)
    }
  }

  # Release the one-nuisance chunk cache before returning the chunk batch
  clear_nuisance_cache()

  if (!is.null(trace_env)) {
    trace_env$last_chunk_start <- a_s
    trace_env$last_chunk_end <- a_e
    trace_env$active_query_chunks <- active_query_chunks
    trace_env$max_simultaneous_query_chunks <- max(
      if (!is.null(trace_env$max_simultaneous_query_chunks)) as.integer(trace_env$max_simultaneous_query_chunks) else 0L,
      max_simultaneous_query_chunks
    )
    trace_env$active_nuisance_caches <- 0L
    trace_env$max_simultaneous_nuisance_caches <- max(
      if (!is.null(trace_env$max_simultaneous_nuisance_caches)) as.integer(trace_env$max_simultaneous_nuisance_caches) else 0L,
      max_simultaneous_nuisance_caches
    )
    trace_env$nuisance_cache_hits <- (if (!is.null(trace_env$nuisance_cache_hits)) as.integer(trace_env$nuisance_cache_hits) else 0L) +
      nuisance_cache_hits
    trace_env$endpoints_read <- (if (!is.null(trace_env$endpoints_read)) as.integer(trace_env$endpoints_read) else 0L) +
      endpoints_read_chunk
  }

  if (isTRUE(callback_schema)) {
    return(list(
      phi = phi_mat,
      id = u_id,
      link_ids = link_ids,
      a_index = a_seq
    ))
  }

  structure(
    list(
      batch_index = as.integer(batch_index),
      a_index = a_seq,
      unit_index = u_idx,
      id = u_id,
      g = u_g,
      weight = u_w,
      link_ids = link_ids,
      phi = phi_mat,
      components = if (include_components) {
        list(
          rate = comp_rate,
          beta_score = comp_score,
          treated = comp_treated,
          control = comp_control
        )
      } else NULL,
      read_diagnostics = list(
        n_endpoints_read = logit_cert_count_output(endpoints_read_chunk),
        risk_filtered_before_outcome_read = TRUE,
        off_risk_endpoints_read = FALSE,
        nonobserver_endpoints_read = FALSE,
        nuisance_cache_hits = nuisance_cache_hits
      ),
      work = list(
        actual_work_ops = actual_ops,
        max_simultaneous_query_chunks = max_simultaneous_query_chunks,
        max_simultaneous_nuisance_caches = max_simultaneous_nuisance_caches
      )
    ),
    class = c("cdid_ipw_point_if_batch", "cdid_ipw_if_batch")
  )
}

evaluate_ipw_point_if_view_chunk <- function(
  point_links,
  a_start,
  a_end,
  control,
  include_components = FALSE,
  callback_schema = FALSE,
  use_nuisance_cache = TRUE,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(a_start) || missing(a_end) || missing(control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`point_links`, `a_start`, `a_end`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!gmm_scalar_logical(include_components) ||
      !gmm_scalar_logical(callback_schema) ||
      !gmm_scalar_logical(use_nuisance_cache)) {
    cdid_ipw_point_view_abort(
      "`include_components`, `callback_schema`, and `use_nuisance_cache` must be scalar logicals.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (isTRUE(callback_schema) && isTRUE(include_components)) {
    cdid_ipw_point_view_abort(
      "`callback_schema = TRUE` requires `include_components = FALSE` so the returned batch matches the 4-field covariance callback schema.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }

  ctrl <- validate_ipw_point_view_control(control)
  norm <- normalize_ipw_point_links_view(point_links, control = ctrl)
  na <- as.integer(norm$frame_owner$n_a_units)
  s_idx <- ipw_owner_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_owner_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > na) {
    cdid_ipw_point_view_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= n_a_units` (", na, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error")
    )
  }
  b_len <- as.integer(e_idx - s_idx + 1)

  # Preflight before any chunk query
  pre_norm <- preflight_ipw_point_if_view(
    point_links = norm$links,
    control = ctrl,
    include_components = include_components,
    chunk_rows = b_len,
    return_norm = TRUE
  )
  norm <- pre_norm$norm
  pre_norm <- NULL

  ipw_point_if_eval_chunk_internal(
    norm = norm,
    a_start = s_idx,
    a_end = e_idx,
    ctrl = ctrl,
    include_components = include_components,
    callback_schema = callback_schema,
    use_nuisance_cache = use_nuisance_cache,
    batch_index = 1L,
    trace_env = trace_env
  )
}

query_ipw_point_if_view_chunk <- evaluate_ipw_point_if_view_chunk

replay_ipw_point_if_view <- function(
  point_links,
  control,
  callback,
  include_components = FALSE,
  callback_schema = FALSE,
  use_nuisance_cache = TRUE,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control) || missing(callback)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`point_links`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.function(callback)) {
    cdid_ipw_point_view_abort(
      "`callback` must be a function.",
      c("cdid_ipw_point_view_input_error", "cdid_ipw_input_error")
    )
  }
  if (!gmm_scalar_logical(include_components) ||
      !gmm_scalar_logical(callback_schema) ||
      !gmm_scalar_logical(use_nuisance_cache)) {
    cdid_ipw_point_view_abort(
      "`include_components`, `callback_schema`, and `use_nuisance_cache` must be scalar logicals.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (isTRUE(callback_schema) && isTRUE(include_components)) {
    cdid_ipw_point_view_abort(
      "`callback_schema = TRUE` requires `include_components = FALSE`.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }

  # 1. Whole-traversal preflight BEFORE scanning any chunk or invoking callback!
  pre_norm <- preflight_ipw_point_if_view(
    point_links = point_links,
    control = control,
    include_components = include_components,
    return_norm = TRUE
  )
  pre <- pre_norm$preflight
  ctrl <- pre_norm$ctrl
  norm <- pre_norm$norm
  pre_norm <- NULL

  na <- as.integer(pre$n_a_units)
  L <- as.integer(pre$n_links)
  step <- as.integer(pre$effective_chunk_rows)

  if (!is.null(trace_env)) {
    trace_env$preflight_completed <- TRUE
    trace_env$batches_emitted <- 0L
    trace_env$units_emitted <- 0L
    trace_env$active_query_chunks <- 0L
    trace_env$max_simultaneous_query_chunks <- 0L
    trace_env$active_nuisance_caches <- 0L
    trace_env$max_simultaneous_nuisance_caches <- 0L
    trace_env$nuisance_cache_hits <- 0L
    trace_env$endpoints_read <- 0L
  }

  emitted_units <- 0L
  batches <- 0L
  max_batch_units <- 0L
  total_work_ops <- 0

  a_start <- 1L
  while (a_start <= na) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(na), as.double(a_start) + step - 1))
    b_len <- a_end - a_start + 1L
    batches <- batches + 1L

    batch <- ipw_point_if_eval_chunk_internal(
      norm = norm,
      a_start = a_start,
      a_end = a_end,
      ctrl = ctrl,
      include_components = include_components,
      callback_schema = callback_schema,
      use_nuisance_cache = use_nuisance_cache,
      batch_index = batches,
      trace_env = trace_env
    )
    if (!isTRUE(callback_schema) && !is.null(batch$work$actual_work_ops)) {
      total_work_ops <- total_work_ops + batch$work$actual_work_ops
    }
    callback(batch)
    batch <- NULL

    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    if (!is.null(trace_env)) {
      trace_env$batches_emitted <- batches
      trace_env$units_emitted <- emitted_units
    }
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != na) {
    cdid_ipw_point_view_abort("Point IF replay unit count mismatch.", c("cdid_ipw_count_error"))
  }

  invisible(list(
    n_a_units = emitted_units,
    n_links = L,
    link_ids = norm$link_ids,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * as.double(L),
    actual_work_ops = total_work_ops,
    preflight = pre
  ))
}

traverse_ipw_point_if_view <- replay_ipw_point_if_view

make_ipw_point_candidate_rows_spec <- function(
  point_links,
  control,
  use_nuisance_cache = TRUE,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`point_links`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!gmm_scalar_logical(use_nuisance_cache)) {
    cdid_ipw_point_view_abort("`use_nuisance_cache` must be a scalar logical.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }

  # Preflight the entire stream upfront before returning the callback spec
  pre_norm <- preflight_ipw_point_if_view(
    point_links = point_links,
    control = control,
    include_components = FALSE,
    return_norm = TRUE
  )
  pre <- pre_norm$preflight
  ctrl <- pre_norm$ctrl
  norm <- pre_norm$norm
  pre_norm <- NULL

  cb <- function(a_start, a_end) {
    ipw_point_if_eval_chunk_internal(
      norm = norm,
      a_start = a_start,
      a_end = a_end,
      ctrl = ctrl,
      include_components = FALSE,
      callback_schema = TRUE,
      use_nuisance_cache = use_nuisance_cache,
      batch_index = 1L,
      trace_env = trace_env
    )
  }

  cb_env <- list2env(list(norm = norm, ctrl = ctrl,
    use_nuisance_cache = use_nuisance_cache, trace_env = trace_env),
    parent = environment(make_ipw_point_candidate_rows_spec))
  environment(cb) <- cb_env
  lockEnvironment(cb_env, bindings = TRUE)
  spec <- ipw_candidate_rows_spec(callback = cb, link_ids = norm$link_ids)
  if (!is.null(trace_env)) {
    trace_env$candidate_rows_preflight <- pre
  }
  spec
}

stream_ipw_point_covariance_view <- function(
  point_links,
  control,
  use_nuisance_cache = TRUE,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`point_links`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  pre_norm <- preflight_ipw_point_if_view(
    point_links = point_links,
    control = control,
    include_components = FALSE,
    return_norm = TRUE
  )
  pre1 <- pre_norm$preflight
  ctrl <- pre_norm$ctrl
  norm <- pre_norm$norm
  pre_norm <- NULL
  n_a <- as.double(pre1$n_a_units)
  L <- as.integer(pre1$n_links)
  link_ids <- norm$link_ids

  cov_entries <- ipw_owner_check_matrix_dims(
    L, L,
    max_matrix_entries = min(ctrl$max_matrix_entries, ctrl$max_output_entries),
    label = "streamed L x L point IF covariance matrix"
  )
  ipw_owner_safe_sum(
    c(pre1$peak_concurrent_entries, 4 * cov_entries, 8 * L, 6 * pre1$b_by_L_entries, 128),
    ctrl$max_matrix_entries,
    "two-pass streamed point covariance concurrent entries"
  )
  ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(2, pre1$total_work_ops_bound, ctrl$max_work_ops, "two-pass point IF replay work"),
      ipw_owner_safe_product(n_a, 2 * L * L + 4 * L, ctrl$max_work_ops, "centering and crossprod accumulation work"),
      128
    ),
    ctrl$max_work_ops,
    "total two-pass streamed point covariance work"
  )

  sum_phi <- numeric(L)
  replay_ipw_point_if_view(
    point_links = norm$links,
    control = ctrl,
    include_components = FALSE,
    callback_schema = TRUE,
    use_nuisance_cache = use_nuisance_cache,
    trace_env = trace_env,
    callback = function(batch) {
      sum_phi <<- sum_phi + colSums(batch$phi)
      if (any(!is.finite(sum_phi))) {
        cdid_abort("Nonfinite IF sum during covariance replay.", "cdid_ipw_numeric_error")
      }
    }
  )
  mean_phi <- sum_phi / n_a
  names(mean_phi) <- link_ids

  cov_mat <- matrix(0.0, nrow = L, ncol = L, dimnames = list(link_ids, link_ids))
  replay_info <- replay_ipw_point_if_view(
    point_links = norm$links,
    control = ctrl,
    include_components = FALSE,
    callback_schema = TRUE,
    use_nuisance_cache = use_nuisance_cache,
    trace_env = trace_env,
    callback = function(batch) {
      centered <- sweep(batch$phi, 2L, mean_phi, FUN = "-")
      if (any(!is.finite(centered))) {
        cdid_abort("Nonfinite centered IF during covariance replay.", "cdid_ipw_numeric_error")
      }
      cov_mat <<- cov_mat + crossprod(centered)
      if (any(!is.finite(cov_mat))) {
        cdid_abort("Nonfinite IF crossproduct or sum during covariance replay.", "cdid_ipw_numeric_error")
      }
      rm(centered)
    }
  )
  cov_mat <- cov_mat / n_a

  list(
    covariance = cov_mat,
    mean_phi = mean_phi,
    estimates = stats::setNames(
      vapply(norm$links, `[[`, numeric(1L), "estimate"),
      link_ids
    ),
    n_a_units = pre1$n_a_units,
    replay = replay_info
  )
}

# ---------------------------------------------------------------------------
# Section 2a: Invocation-local general normalization record (CP08 missing-panel
# IPW performance). On the general route, `run_chained_ipw_pipeline` performs
# the point-link normalization and whole-replay preflight once
# (`preflight_ipw_point_if_view(..., return_norm = TRUE)`) and seals the exact
# normalization, its validated-view records and the preflight together with
# every object that normalization validated: the frame and panel owners, each
# distinct fitted nuisance view and nuisance view, and every link's link view
# and point link. The two covariance passes and the two GMM projection passes
# consume the record through the private `*_general` variants below instead of
# repeating the same validation of the same immutable objects. Each pass
# verifies the record at its boundary and after every callback, before the next
# endpoint read: the record's own seal, identical controls (hence the same
# chunk step) and flags, and the live locked binding of every bound object
# through the existing native alias check (`ipw_owner_seal_ok`) plus the
# existing validated-view record match for every link view. The record is not
# a certification token, is never stored beyond the pipeline invocation and
# never reaches a public argument. The standalone wrappers
# (`replay_ipw_point_if_view`, `stream_ipw_point_covariance_view`,
# `evaluate_ipw_point_if_view_chunk`, `make_ipw_point_candidate_rows_spec`)
# keep their full per-call validation unchanged; the complete batch plan route
# (Section 2b) stays distinct.
# ---------------------------------------------------------------------------

ipw_point_view_general_abort <- function(message, diagnostics = list()) {
  cdid_ipw_point_view_abort(message,
    c("cdid_ipw_point_view_general_error", "cdid_ipw_point_view_integrity_error", "cdid_ipw_owner_integrity_error"),
    diagnostics = diagnostics)
}

# One normalization and preflight (exactly `preflight_ipw_point_if_view` with
# `return_norm = TRUE`, `include_components = FALSE`, the control's chunk rows),
# then a sealed record of every validated object. The record's header and the
# verifications of its four planned passes (two covariance, two projection;
# one verification per pass boundary and one after each chunk callback) are
# charged against the caller's caps before the record exists; both are
# upper bounds on this record's own work, not on any later caller's. Nothing
# here reads an outcome or copies an owner array: every field references the
# existing sealed objects.
ipw_point_if_general_plan <- function(point_links, control) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`point_links`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  pre_norm <- preflight_ipw_point_if_view(point_links, control, include_components = FALSE, return_norm = TRUE)
  pre <- pre_norm$preflight
  ctrl <- pre_norm$ctrl
  norm <- pre_norm$norm
  pre_norm <- NULL
  L <- as.integer(norm$n_links)
  fo <- norm$frame_owner
  po <- norm$panel_owner
  if (!is.list(norm$validated) || length(norm$validated) != L ||
      !is.environment(fo$ref) || !is.environment(fo$panel_ref) || !identical(fo$panel_ref, po$ref)) {
    ipw_point_view_general_abort("General point-IF normalization returned an inconsistent validated-view record set.")
  }
  # Before reference/index allocation or distinct-reference searches, reserve
  # the header using n_nuisance <= n_fitted (each fitted view names exactly
  # one nuisance view). This covers construction's simultaneously live lists,
  # indices, seal/header and later consumer use, without traversing payloads.
  max_lim <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  n_fitted_bound <- as.double(norm$n_unique_fitted_nuisances)
  record_entries <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(8, L, max_lim, "general point-IF record per-link references and indices"),
    ipw_owner_safe_product(8, n_fitted_bound, max_lim, "general point-IF record distinct view references"),
    256), max_lim, "general point-IF normalization record header")
  ipw_owner_safe_sum(c(pre$peak_concurrent_entries, record_entries), ctrl$max_matrix_entries,
    "general point-IF normalization record plus point-IF producer concurrent entries")
  construction_work <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(16 * as.double(L), L, ctrl$max_work_ops, "general record distinct-reference searches"),
    64 * as.double(L), 256), ctrl$max_work_ops, "general record construction work")
  n_chunks <- as.double(pre$n_chunks)
  one_verification <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(64, L, ctrl$max_work_ops, "per-link seal, binding and record identity checks"),
    ipw_owner_safe_product(128, n_fitted_bound, ctrl$max_work_ops, "distinct fitted/nuisance view seal and binding checks"),
    256), ctrl$max_work_ops, "one general point-IF record verification")
  # The covariance entry check is additional to its two replay boundaries.
  verification_count <- 1 + 4 * (n_chunks + 1)
  verification_work <- ipw_owner_safe_product(verification_count, one_verification, ctrl$max_work_ops,
    "general point-IF record verifications across covariance entry and four passes")
  additional_work <- ipw_owner_safe_sum(c(construction_work, verification_work), ctrl$max_work_ops,
    "general point-IF record construction and verification work")
  ipw_owner_safe_sum(c(pre$total_work_ops_bound, additional_work), ctrl$max_work_ops,
    "general point-IF replay work plus record construction and verifications")
  # The distinct fitted nuisance views and nuisance views exactly as the
  # normalization resolved them (one locked reference each), and each link's
  # index into those lists so a later verification is one identity test per
  # link rather than a scan.
  link_views <- vector("list", L)
  link_fitted_index <- integer(L)
  link_nuisance_index <- integer(L)
  fitted_refs <- list()
  fitted_views <- list()
  nuisance_refs <- list()
  nuisance_views <- list()
  for (ell in seq_len(L)) {
    lk <- norm$links[[ell]]
    lv <- get("link_view", envir = lk$link_view_ref, inherits = FALSE)
    if (!identical(norm$validated[[ell]]$view, lv) || !identical(lk$link_ref, lk$link_view_ref)) {
      ipw_point_view_general_abort("General point-IF normalization record cannot bind a link view that its validated record does not name.")
    }
    link_views[[ell]] <- lv
    fref <- lk$fitted_nuisance_ref
    k_f <- 0L
    for (k in seq_along(fitted_refs)) if (identical(fitted_refs[[k]], fref)) { k_f <- k; break }
    if (k_f == 0L) {
      fitted_refs[[length(fitted_refs) + 1L]] <- fref
      fitted_views[[length(fitted_views) + 1L]] <- get("fitted_nuisance_view", envir = fref, inherits = FALSE)
      k_f <- length(fitted_refs)
    }
    link_fitted_index[ell] <- k_f
    nref <- lk$nuisance_view_ref
    k_n <- 0L
    for (k in seq_along(nuisance_refs)) if (identical(nuisance_refs[[k]], nref)) { k_n <- k; break }
    if (k_n == 0L) {
      nuisance_refs[[length(nuisance_refs) + 1L]] <- nref
      nuisance_views[[length(nuisance_views) + 1L]] <- get("nuisance_view", envir = nref, inherits = FALSE)
      k_n <- length(nuisance_refs)
    }
    link_nuisance_index[ell] <- k_n
  }
  n_fv <- length(fitted_views)
  n_nv <- length(nuisance_views)
  if (!identical(n_fv, as.integer(norm$n_unique_fitted_nuisances))) {
    ipw_point_view_general_abort("General point-IF normalization record resolved a different set of fitted nuisances than its normalization.")
  }
  # Each distinct fitted view names one of the distinct nuisance views.
  fitted_nuisance_index <- integer(n_fv)
  for (k in seq_len(n_fv)) {
    nref <- fitted_views[[k]]$nuisance_view_ref
    k_n <- 0L
    for (j in seq_len(n_nv)) if (identical(nuisance_refs[[j]], nref)) { k_n <- j; break }
    if (k_n == 0L) {
      ipw_point_view_general_abort("General point-IF normalization record found a fitted nuisance whose nuisance view no link names.")
    }
    fitted_nuisance_index[k] <- k_n
  }

  if (n_nv > n_fv) {
    ipw_point_view_general_abort("General point-IF normalization record resolved more nuisance views than fitted views.")
  }

  plan <- list(
    kind = "cdid_ipw_point_if_general_plan",
    control = ctrl,
    control_in = control,
    frame_owner = fo,
    panel_owner = po,
    point_links = norm$links,
    link_views = link_views,
    link_ids = norm$link_ids,
    n_links = L,
    fitted_refs = fitted_refs,
    fitted_views = fitted_views,
    link_fitted_index = link_fitted_index,
    nuisance_refs = nuisance_refs,
    nuisance_views = nuisance_views,
    link_nuisance_index = link_nuisance_index,
    fitted_nuisance_index = fitted_nuisance_index,
    norm = norm,
    preflight = pre,
    n_a_units = as.integer(fo$n_a_units),
    step = as.integer(pre$effective_chunk_rows),
    include_components = FALSE,
    general = list(
      route = "general_normalization_record",
      normalizations_performed = 1L,
      n_links = L,
      n_distinct_fitted_nuisances = n_fv,
      n_distinct_nuisance_views = n_nv,
      effective_chunk_rows = pre$effective_chunk_rows,
      n_chunks = pre$n_chunks,
      passes_planned = 4L,
      verifications_planned = verification_count,
      construction_work_ops_bound = construction_work,
      one_verification_work_ops_bound = one_verification,
      additional_work_ops_bound = additional_work,
      record_header_entries = record_entries,
      additional_live_entries = record_entries,
      verification_work_ops_bound = verification_work,
      full_n_by_L_retained = FALSE,
      per_link_dense_masks_retained = FALSE,
      peak_rss_promised = FALSE,
      lifetime = paste0("References to the existing sealed owners, views and point links persist from the pipeline's ",
        "producer preflight through the GMM projection and are dropped after the solve; per pass one chunk's ",
        "query and output are live with the consumer, then released, exactly as in the standalone replay.")
    )
  )
  ipw_owner_seal(plan, "point_if_general_plan")
}

# Record check at every pass boundary and after every callback, before the
# next endpoint read: the record's seal, the caller's controls and flags, and
# the live locked binding of every object the normalization validated. Every
# check is a seal/alias/identity test on an existing object; no payload is
# traversed again. A record that fails any check refuses; nothing falls back
# to a repeated normalization here.
ipw_point_if_general_verify <- function(plan, control, include_components = FALSE) {
  if (!is.list(plan) || !identical(plan$kind, "cdid_ipw_point_if_general_plan") ||
      !ipw_owner_seal_ok(plan, "point_if_general_plan")) {
    ipw_point_view_general_abort("General point-IF normalization record is malformed, foreign, or differs from its sealed snapshot.")
  }
  if (!identical(control, plan$control_in) && !identical(control, plan$control)) {
    cdid_ipw_point_view_abort("General point-IF normalization record controls differ from the traversal controls; reuse is refused.",
      c("cdid_ipw_control_mismatch_error", "cdid_ipw_point_view_general_error"))
  }
  if (!gmm_scalar_logical(include_components) || !identical(isTRUE(include_components), FALSE)) {
    cdid_ipw_point_view_abort(
      "General point-IF normalization record was preflighted without IF components; component output must use the standalone replay.",
      c("cdid_ipw_point_view_input_error", "cdid_ipw_point_view_general_error"))
  }
  norm <- plan$norm
  fo <- plan$frame_owner
  po <- plan$panel_owner
  L <- plan$n_links
  bound <- is.list(norm) && identical(norm$frame_owner, fo) && identical(norm$panel_owner, po) &&
    identical(as.integer(norm$n_links), L) && identical(norm$link_ids, plan$link_ids) &&
    is.list(norm$validated) && length(norm$validated) == L && length(plan$point_links) == L &&
    length(plan$link_views) == L && length(plan$link_fitted_index) == L && length(plan$link_nuisance_index) == L &&
    ipw_owner_seal_ok(fo, "frame_owner") && identical(fo$panel_ref, po$ref) && ipw_owner_seal_ok(po, "panel_owner")
  n_fv <- length(plan$fitted_views)
  n_nv <- length(plan$nuisance_views)
  bound <- bound && length(plan$fitted_refs) == n_fv && length(plan$nuisance_refs) == n_nv &&
    length(plan$fitted_nuisance_index) == n_fv
  for (k in seq_len(n_nv)) {
    if (!bound) break
    nv <- plan$nuisance_views[[k]]
    bound <- identical(nv$ref, plan$nuisance_refs[[k]]) && ipw_owner_seal_ok(nv, "nuisance_view") &&
      identical(nv$frame_ref, fo$ref) && identical(nv$panel_ref, fo$panel_ref)
  }
  for (k in seq_len(n_fv)) {
    if (!bound) break
    fv <- plan$fitted_views[[k]]
    j <- plan$fitted_nuisance_index[k]
    bound <- identical(fv$ref, plan$fitted_refs[[k]]) && ipw_owner_seal_ok(fv, "fitted_nuisance_view") &&
      identical(fv$frame_ref, fo$ref) && identical(fv$panel_ref, fo$panel_ref) &&
      j >= 1L && j <= n_nv && identical(fv$nuisance_view_ref, plan$nuisance_refs[[j]])
  }
  for (ell in seq_len(L)) {
    if (!bound) break
    pl <- plan$point_links[[ell]]
    lv <- plan$link_views[[ell]]
    k_f <- plan$link_fitted_index[ell]
    k_n <- plan$link_nuisance_index[ell]
    bound <- identical(norm$links[[ell]], pl) && ipw_owner_seal_ok(pl, "point_link_view") &&
      ipw_owner_seal_ok(lv, "link_view") && identical(pl$link_view_ref, lv$ref) && identical(pl$link_ref, lv$ref) &&
      identical(pl$frame_ref, fo$ref) && identical(pl$panel_ref, fo$panel_ref) &&
      k_f >= 1L && k_f <= n_fv && identical(pl$fitted_nuisance_ref, plan$fitted_refs[[k_f]]) &&
      k_n >= 1L && k_n <= n_nv && identical(pl$nuisance_view_ref, plan$nuisance_refs[[k_n]]) &&
      identical(lv$frame_ref, fo$ref) && identical(lv$panel_ref, fo$panel_ref) &&
      (is.null(pl$link_id) || !nzchar(pl$link_id) || identical(pl$link_id, plan$link_ids[[ell]])) &&
      ipw_owner_validated_view_match(norm$validated[[ell]], lv)
  }
  if (!bound) {
    ipw_point_view_general_abort("General point-IF normalization record owner, nuisance, link or point-link bindings changed; reuse is refused.")
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# D90 missing-panel IPW: invocation-local, outcome-free sparse replay-state
# record bound to a verified `cdid_ipw_point_if_general_plan`.
#   Let S = sum_ell S_ell (observed risk units across all L links), C = n_chunks,
#   B = effective_chunk_rows. Raw state has E = 6*S + (C + 1)*L slots (six
#   length-S_ell vectors per link: observer A-index, ordered cohort index,
#   left/right observation positions, exact odds_k, exact score_base_k, plus
#   length-(C + 1) chunk offsets). Capacity admission reserves
#   2*E + 320*L + 2048 + 16*B additional entries. The second E covers
#   construction/COW copies. Per-link headers include seven vector headers,
#   record/names/seal, scalars, and outer-list/cursor slots; the fixed reserve
#   includes reservation, binding, runtime and outer-state headers. Slice/index
#   scratch is released at each link. These are counted slots, not process RSS.
#   Work charges the original producer construction bound, (C + 16)*E plus
#   7*C*L^2 for initialization/fill, worst-case vector/outer-list replacement
#   copies and four slice passes, boundary state checks, chunk/link guards,
#   and both extra original plan verifications. Admission uses the full live/work
#   envelope and the unchanged producer/covariance/replay caps before storage
#   allocation or queries; the row partition is never reduced to admit state.
#   Typed capacity decline returns NULL before state construction so the
#   caller runs the unchanged general replay; invalid counts, controls,
#   identities, orders or numerical failures propagate.
# ---------------------------------------------------------------------------

ipw_point_if_general_state_entries <- function(S, C, L, B, max_lim = NULL) {
  six_S <- ipw_owner_safe_product(6, S, max_lim, "general point-IF replay state observer slots")
  offset_slots <- ipw_owner_safe_product(C + 1, L, max_lim, "general point-IF replay state chunk offset slots")
  E <- ipw_owner_safe_sum(c(six_S, offset_slots), max_lim, "general point-IF replay state raw entries E")
  two_E <- ipw_owner_safe_product(2, E, max_lim, "general point-IF replay state plus construction copy reserve")
  header_slots <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(320, L, max_lim, "general point-IF replay state per-link headers"), 2048),
    max_lim,
    "general point-IF replay state headers and seals"
  )
  slice_slots <- ipw_owner_safe_product(16, B, max_lim, "general point-IF replay state chunk slice and index scratch")
  entries <- ipw_owner_safe_sum(
    c(two_E, header_slots, slice_slots),
    max_lim,
    "general point-IF replay state reserved entries"
  )
  list(
    raw_entries = E,
    reserved_entries = entries,
    header_entries = header_slots,
    slice_entries = slice_slots
  )
}

ipw_point_if_general_state_work <- function(pre_work, E, C, L, max_wo = NULL, plan_check_work = 0) {
  construction_work <- ipw_owner_safe_sum(c(pre_work), max_wo, "general point-IF replay state producer construction work")
  copy_slice_work <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(C + 16, E, max_wo,
      "general point-IF replay state fill, vector replacement copies and slice work"),
    ipw_owner_safe_product(ipw_owner_safe_product(7 * C, L, max_wo,
      "general point-IF replay state outer-list replacements"), L, max_wo,
      "general point-IF replay state worst-case outer-list replacement copies")
  ), max_wo, "general point-IF replay state construction copies and slice work")
  one_state_check <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(64, L, max_wo, "general point-IF replay state per-link check"), 256),
    max_wo,
    "one general point-IF replay state check"
  )
  state_checks_planned <- 1 + 4 * (C + 1)
  state_check_work <- ipw_owner_safe_product(
    state_checks_planned,
    one_state_check,
    max_wo,
    "general point-IF replay state verification work"
  )
  chunk_guard_work <- ipw_owner_safe_product(
    4 * C, one_state_check, max_wo,
    "general point-IF replay state per-chunk and per-link slice guards"
  )
  setup_work <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(2, plan_check_work, max_wo,
        "general point-IF replay state reserve and build plan verifications"), one_state_check),
    max_wo, "general point-IF replay state count, reservation and seal setup work"
  )
  total_work <- ipw_owner_safe_sum(
    c(construction_work, copy_slice_work, state_check_work, chunk_guard_work, setup_work),
    max_wo,
    "general point-IF replay state total additional work"
  )
  list(
    construction_work = construction_work,
    copy_slice_work = copy_slice_work,
    one_state_check = one_state_check,
    state_checks_planned = state_checks_planned,
    state_check_work = state_check_work,
    chunk_guard_work = chunk_guard_work,
    setup_work = setup_work,
    total_work = total_work
  )
}

ipw_point_if_general_state_reserve <- function(
  plan,
  control,
  live = 0,
  work = 0,
  envelope_control = NULL,
  .verified_entry = FALSE
) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`plan`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!isTRUE(.verified_entry)) {
    ipw_point_if_general_verify(plan, control)
  }
  scalar_nonneg <- function(x) is.numeric(x) && !is.object(x) && length(x) == 1L && is.finite(x) && x >= 0
  scalar_pos <- function(x) scalar_nonneg(x) && x > 0
  count_pos <- function(x) scalar_nonneg(x) && x >= 1 && x == floor(x) && x <= .Machine$integer.max
  count_nonneg <- function(x) scalar_nonneg(x) && x == floor(x) && x <= .Machine$integer.max
  if (!scalar_nonneg(live) || !scalar_nonneg(work)) {
    cdid_ipw_point_view_abort(
      "`live` and `work` must be finite non-negative scalar numbers.",
      "cdid_ipw_point_view_input_error"
    )
  }
  live <- ipw_owner_int_bound(live, "general replay state external live entries", min_val = 0)
  work <- ipw_owner_int_bound(work, "general replay state external work", min_val = 0)
  if (is.null(envelope_control)) envelope_control <- plan$control
  if (!is.list(envelope_control) ||
      !scalar_pos(envelope_control$max_retained_entries) ||
      !scalar_pos(envelope_control$max_matrix_entries) ||
      !scalar_pos(envelope_control$max_output_entries) ||
      !scalar_pos(envelope_control$max_work_ops)) {
    cdid_ipw_point_view_abort(
      "`envelope_control` must provide finite positive retained/matrix/output/work caps.",
      "cdid_ipw_point_view_input_error"
    )
  }

  for (nm in c("max_retained_entries", "max_matrix_entries", "max_output_entries", "max_work_ops")) {
    ipw_owner_int_bound(envelope_control[[nm]], paste0("envelope_control$", nm))
  }
  pre <- plan$preflight
  ctrl <- plan$control
  N <- plan$n_a_units
  L <- plan$n_links
  B <- plan$step
  if (! is.list(pre) || !count_pos(N) || !count_pos(L) || !count_pos(B) ||
      B > N || B > ctrl$max_chunk_rows ||
      !identical(as.double(pre$n_a_units), as.double(N)) ||
      !identical(as.double(pre$n_links), as.double(L)) ||
      !identical(as.double(pre$effective_chunk_rows), as.double(B)) ||
      !count_pos(pre$n_chunks) ||
      !identical(as.double(pre$n_chunks), ceiling(as.double(N) / as.double(B)))) {
    cdid_ipw_point_view_abort(
      "General point-IF normalization record dimensions or chunk partition are inconsistent.",
      c("cdid_ipw_count_error", "cdid_ipw_point_view_general_error", "cdid_ipw_point_view_integrity_error")
    )
  }
  C <- as.integer(pre$n_chunks)
  S_double <- 0
  for (ell in seq_len(L)) {
    pl_ell <- plan$point_links[[ell]]
    lv_ell <- plan$link_views[[ell]]
    s_pl <- if (is.list(pl_ell$counts)) pl_ell$counts$n_observed_risk_units else NULL
    s_lv <- if (is.list(lv_ell$support)) lv_ell$support$n_observed_risk_units else NULL
    if (!count_nonneg(s_pl) || !count_nonneg(s_lv) || s_pl > N ||
        !identical(as.integer(s_pl), as.integer(s_lv))) {
      cdid_ipw_point_view_abort(
        "General point-IF replay state link observed-risk counts are invalid or inconsistent.",
        c("cdid_ipw_count_error", "cdid_ipw_point_view_general_error", "cdid_ipw_point_view_integrity_error")
      )
    }
    S_double <- ipw_owner_safe_sum(c(S_double, s_pl),
      label = "general point-IF replay state observed-risk count")
  }
  if (!scalar_nonneg(pre$total_expected_endpoints_read) ||
      !identical(2 * S_double, as.double(pre$total_expected_endpoints_read))) {
    cdid_ipw_point_view_abort(
      "General point-IF replay state total observed-risk count disagrees with preflight expected endpoint reads.",
      c("cdid_ipw_count_error", "cdid_ipw_point_view_general_error", "cdid_ipw_point_view_integrity_error")
    )
  }

  max_lim_ctrl <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  max_lim_env <- min(envelope_control$max_retained_entries, envelope_control$max_matrix_entries)
  max_lim <- min(max_lim_ctrl, max_lim_env)
  max_wo <- min(ctrl$max_work_ops, envelope_control$max_work_ops)
  N_d <- as.double(N)
  L_d <- as.double(L)
  B_d <- as.double(B)
  C_d <- as.double(C)

  # Exact-count arithmetic and malformed inputs must propagate even when a
  # separate capacity cap would also decline. Compute the inventory uncapped.
  ent_info <- ipw_point_if_general_state_entries(S_double, C_d, L_d, B_d)
  wrk_info <- ipw_point_if_general_state_work(pre$total_work_ops_bound,
    ent_info$raw_entries, C_d, L_d,
    plan_check_work = plan$general$one_verification_work_ops_bound)
  queries_expected <- ipw_owner_safe_product(C_d, L_d,
    label = "general point-IF replay state construction queries")
  admitted <- tryCatch({
    ent_info <- ipw_point_if_general_state_entries(S_double, C_d, L_d, B_d, max_lim = max_lim)
    E <- ent_info$raw_entries
    entries <- ent_info$reserved_entries
    wrk_info <- ipw_point_if_general_state_work(pre$total_work_ops_bound, E, C_d, L_d,
      max_wo = max_wo, plan_check_work = plan$general$one_verification_work_ops_bound)
    state_work <- wrk_info$total_work

    admission_total <- ipw_owner_safe_sum(
      c(live, entries),
      max_lim_env,
      "pipeline simultaneous live entries plus general point-IF replay state"
    )
    work_total <- ipw_owner_safe_sum(
      c(work, state_work),
      envelope_control$max_work_ops,
      "pipeline cumulative work plus general point-IF replay state"
    )

    ipw_owner_safe_sum(
      c(pre$peak_concurrent_entries, plan$general$additional_live_entries, entries),
      max_lim_ctrl,
      "general point-IF producer plus normalization record and replay state entries"
    )
    cov_entries <- ipw_owner_check_matrix_dims(
      L, L,
      max_matrix_entries = min(ctrl$max_matrix_entries, ctrl$max_output_entries),
      label = "streamed L x L point IF covariance matrix"
    )
    ipw_owner_safe_sum(
      c(pre$peak_concurrent_entries, 4 * cov_entries, 8 * L_d, 6 * pre$b_by_L_entries,
        plan$general$additional_live_entries, entries, 128),
      ctrl$max_matrix_entries,
      "two-pass streamed point covariance concurrent entries plus general replay state"
    )
    ipw_owner_safe_sum(
      c(
        ipw_owner_safe_product(2, pre$total_work_ops_bound, ctrl$max_work_ops, "two-pass point IF replay work"),
        ipw_owner_safe_product(1 + 2 * (C_d + 1), plan$general$one_verification_work_ops_bound, ctrl$max_work_ops,
          "general covariance entry and two replay record verifications"),
        ipw_owner_safe_product(N_d, 2 * L_d * L_d + 4 * L_d, ctrl$max_work_ops, "centering and crossprod accumulation work"),
        state_work,
        128
      ),
      ctrl$max_work_ops,
      "total two-pass streamed point covariance work plus general replay state"
    )
    ipw_owner_safe_sum(
      c(pre$total_work_ops_bound, plan$general$additional_work_ops_bound, state_work),
      ctrl$max_work_ops,
      "general point-IF plan work plus replay state construction, slicing and checks"
    )
    list(
      E = E,
      entries = entries,
      wrk_info = wrk_info,
      admission_total = admission_total,
      work_total = work_total
    )
  }, cdid_ipw_owner_bound_error = function(e) {
    d <- e$diagnostics
    if (is.list(d) && scalar_nonneg(d$requested) && scalar_pos(d$bound) &&
        d$requested > d$bound) return(NULL)
    stop(e)
  })

  if (is.null(admitted)) return(NULL)

  res <- new.env(parent = emptyenv())
  res$kind <- "cdid_ipw_point_if_general_state_reservation"
  res$state <- "reserved"
  res$plan_ref <- plan$ref
  res$frame_ref <- plan$frame_owner$ref
  res$panel_ref <- plan$frame_owner$panel_ref
  res$plan_snapshot <- plan
  lockBinding("plan_snapshot", res)
  res$N <- as.integer(N)
  res$L <- as.integer(L)
  res$B <- as.integer(B)
  res$C <- as.integer(C)
  res$S <- S_double
  res$s_vec <- vapply(plan$point_links,
    function(lk) as.integer(lk$counts$n_observed_risk_units), integer(1L), USE.NAMES = FALSE)
  res$queries_expected <- queries_expected
  res$link_ids <- plan$link_ids
  res$raw_entries <- admitted$E
  res$entries <- admitted$entries
  res$work <- admitted$wrk_info$total_work
  res$construction_work_ops_bound <- admitted$wrk_info$construction_work
  res$copy_slice_work_ops_bound <- admitted$wrk_info$copy_slice_work
  res$one_state_check_work_ops_bound <- admitted$wrk_info$one_state_check
  res$state_checks_planned <- admitted$wrk_info$state_checks_planned
  res$state_check_work_ops_bound <- admitted$wrk_info$state_check_work
  res$chunk_guard_work_ops_bound <- admitted$wrk_info$chunk_guard_work
  res$setup_work_ops_bound <- admitted$wrk_info$setup_work
  res$admission_total <- admitted$admission_total
  res$work_total <- admitted$work_total
  res$cap <- max_lim
  res$work_cap <- max_wo
  lockEnvironment(res, bindings = TRUE)
  unlockBinding("state", res)
  res
}

ipw_point_if_general_state_verify <- function(state, plan, control) {
  # Callers retain the original full plan verification at every boundary.
  # This additional check verifies only the sparse state and its plan binding.
  if (!ipw_owner_seal_ok(plan, "point_if_general_plan") ||
      (!identical(control, plan$control_in) && !identical(control, plan$control))) {
    ipw_point_view_general_abort("General point-IF replay state plan or controls changed; reuse is refused.")
  }
  if (!is.list(state) || !identical(state$kind, "cdid_ipw_point_if_general_state") ||
      !ipw_owner_seal_ok(state, "point_if_general_state") ||
      !is.environment(state$runtime) || !identical(state$runtime$status, "built") ||
      !is.environment(state$binding_ref) || !identical(parent.env(state$binding_ref), emptyenv()) ||
      !environmentIsLocked(state$binding_ref) ||
      !exists("plan_snapshot", envir = state$binding_ref, inherits = FALSE) ||
      !bindingIsLocked("plan_snapshot", state$binding_ref) ||
      !gmm_sexp_aliases(plan, get("plan_snapshot", envir = state$binding_ref, inherits = FALSE))) {
    ipw_point_view_general_abort("General point-IF replay state is malformed, foreign, released, or differs from its sealed snapshot.")
  }
  fo <- plan$frame_owner
  L <- plan$n_links
  C <- as.integer(plan$preflight$n_chunks)
  ok <- identical(state$plan_ref, plan$ref) &&
    identical(state$frame_ref, fo$ref) &&
    identical(state$panel_ref, fo$panel_ref) &&
    identical(state$n_a_units, plan$n_a_units) &&
    identical(state$n_links, L) &&
    identical(state$step, plan$step) &&
    identical(state$n_chunks, C) &&
    identical(state$link_ids, plan$link_ids) &&
    is.list(state$link_states) && length(state$link_states) == L &&
    is.list(state$link_state_refs) && length(state$link_state_refs) == L
  for (ell in seq_len(L)) {
    if (!ok) break
    ls_ell <- state$link_states[[ell]]
    pl_ell <- plan$point_links[[ell]]
    lv_ell <- plan$link_views[[ell]]
    S_ell <- as.integer(pl_ell$counts$n_observed_risk_units)
    offs <- if (is.list(ls_ell)) ls_ell$chunk_offsets else NULL
    ok <- is.list(ls_ell) &&
      identical(ls_ell$kind, "cdid_ipw_point_if_general_link_state") &&
      identical(ls_ell$ref, state$link_state_refs[[ell]]) &&
      ipw_owner_seal_ok(ls_ell, "point_if_general_link_state") &&
      identical(ls_ell$link_index, as.integer(ell)) &&
      identical(ls_ell$link_id, plan$link_ids[[ell]]) &&
      identical(ls_ell$point_link_ref, pl_ell$ref) &&
      identical(ls_ell$link_view_ref, lv_ell$ref) &&
      identical(ls_ell$fitted_nuisance_ref, pl_ell$fitted_nuisance_ref) &&
      identical(ls_ell$nuisance_view_ref, pl_ell$nuisance_view_ref) &&
      identical(ls_ell$n_observed_risk_units, S_ell) &&
      identical(ls_ell$n_chunks, C) &&
      is.integer(offs) && length(offs) == C + 1 &&
      !is.na(offs[1L]) && !is.na(offs[C + 1]) &&
      identical(offs[1L], 0L) && identical(offs[C + 1], S_ell) &&
      is.integer(ls_ell$a_index) && length(ls_ell$a_index) == S_ell &&
      is.integer(ls_ell$cohort_index) && length(ls_ell$cohort_index) == S_ell &&
      is.integer(ls_ell$obs_s_index) && length(ls_ell$obs_s_index) == S_ell &&
      is.integer(ls_ell$obs_t_index) && length(ls_ell$obs_t_index) == S_ell &&
      is.double(ls_ell$odds) && length(ls_ell$odds) == S_ell &&
      is.double(ls_ell$score_base) && length(ls_ell$score_base) == S_ell
  }
  if (!ok) {
    ipw_point_view_general_abort("General point-IF replay state link bindings, offsets, or vector storage changed; reuse is refused.")
  }
  state$runtime$verifications <- logit_cert_count_output(as.double(state$runtime$verifications) + 1)
  invisible(TRUE)
}

ipw_point_if_general_state_link_slice <- function(state, ell, batch_index, a_s, a_e, b_len, fit, lv_ell) {
  ls_ell <- state$link_states[[ell]]
  C <- state$n_chunks
  if (!is.list(ls_ell) || !identical(ls_ell$kind, "cdid_ipw_point_if_general_link_state") ||
      !identical(ls_ell$ref, state$link_state_refs[[ell]]) ||
      !ipw_owner_seal_ok(ls_ell, "point_if_general_link_state") ||
      !identical(ls_ell$link_index, as.integer(ell)) ||
      !identical(ls_ell$point_link_ref, fit$ref) ||
      !identical(ls_ell$link_view_ref, lv_ell$ref) ||
      !identical(ls_ell$fitted_nuisance_ref, fit$fitted_nuisance_ref) ||
      !identical(ls_ell$nuisance_view_ref, fit$nuisance_view_ref)) {
    ipw_point_view_general_abort("General point-IF replay state link binding is malformed or disagrees with the active point link.")
  }
  S_ell <- ls_ell$n_observed_risk_units
  offs <- ls_ell$chunk_offsets
  if (!is.integer(offs) || length(offs) != C + 1 ||
      is.na(offs[1L]) || is.na(offs[C + 1]) ||
      !identical(offs[1L], 0L) || !identical(offs[C + 1], S_ell) ||
      is.na(offs[batch_index]) || is.na(offs[batch_index + 1L])) {
    ipw_point_view_general_abort("General point-IF replay state chunk offsets are malformed.")
  }
  o_prev <- offs[batch_index]
  o_end <- offs[batch_index + 1L]
  if (o_prev < 0L || o_end < o_prev || o_end > S_ell ||
      (batch_index == 1L && o_prev != 0L) ||
      (batch_index == C && o_end != S_ell)) {
    ipw_point_view_general_abort("General point-IF replay state chunk offset range is out of bounds.")
  }
  k_ch <- as.integer(o_end - o_prev)
  if (k_ch > b_len) {
    ipw_point_view_general_abort("General point-IF replay state observer count exceeds chunk length.")
  }
  if (k_ch > 0L) {
    o_seq <- seq.int(o_prev + 1L, o_end)
    list(
      k_ch = k_ch,
      a_index = ls_ell$a_index[o_seq],
      cohort_index = ls_ell$cohort_index[o_seq],
      obs_s_index = ls_ell$obs_s_index[o_seq],
      obs_t_index = ls_ell$obs_t_index[o_seq],
      odds = ls_ell$odds[o_seq],
      score_base = ls_ell$score_base[o_seq]
    )
  } else {
    list(
      k_ch = 0L,
      a_index = integer(0L),
      cohort_index = integer(0L),
      obs_s_index = integer(0L),
      obs_t_index = integer(0L),
      odds = numeric(0L),
      score_base = numeric(0L)
    )
  }
}

ipw_point_if_general_state_build <- function(
  reservation,
  plan,
  control,
  use_nuisance_cache = TRUE
) {
  ipw_point_view_require_deps()
  if (missing(reservation) || missing(plan) || missing(control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`reservation`, `plan`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!gmm_scalar_logical(use_nuisance_cache)) {
    cdid_ipw_point_view_abort("`use_nuisance_cache` must be a scalar logical.", "cdid_ipw_point_view_input_error")
  }
  ipw_point_if_general_verify(plan, control)
  if (!is.environment(reservation) ||
      !identical(reservation$kind, "cdid_ipw_point_if_general_state_reservation") ||
      !identical(reservation$state, "reserved") ||
      !identical(parent.env(reservation), emptyenv()) || !environmentIsLocked(reservation) ||
      !all(vapply(setdiff(ls(reservation, all.names = TRUE), "state"),
        bindingIsLocked, logical(1L), env = reservation)) ||
      !exists("plan_snapshot", envir = reservation, inherits = FALSE) ||
      !bindingIsLocked("plan_snapshot", reservation) ||
      !gmm_sexp_aliases(plan, get("plan_snapshot", envir = reservation, inherits = FALSE)) ||
      !identical(plan$ref, reservation$plan_ref) ||
      !identical(plan$frame_owner$ref, reservation$frame_ref) ||
      !identical(plan$frame_owner$panel_ref, reservation$panel_ref)) {
    ipw_point_view_general_abort("General point-IF replay state reservation is invalid, already consumed, or does not match the sealed plan.")
  }

  ctrl <- plan$control
  owner_ctrl <- ipw_point_view_to_owner_control(ctrl)
  norm <- plan$norm
  fo <- plan$frame_owner
  po <- plan$panel_owner
  na <- reservation$N
  L <- reservation$L
  step <- reservation$B
  n_chunks <- reservation$C
  s_vec <- reservation$s_vec
  n_obs_y <- length(po$observations$y)

  reservation$state <- "building"
  on.exit({
    if (!identical(reservation$state, "built")) reservation$state <- "failed"
    unlockBinding("plan_snapshot", reservation)
    assign("plan_snapshot", NULL, envir = reservation)
    lockBinding("plan_snapshot", reservation)
  }, add = TRUE)

  a_index_list <- vector("list", L)
  cohort_index_list <- vector("list", L)
  obs_s_list <- vector("list", L)
  obs_t_list <- vector("list", L)
  odds_list <- vector("list", L)
  score_base_list <- vector("list", L)
  offsets_list <- vector("list", L)
  cursors <- integer(L)

  for (ell in seq_len(L)) {
    s_ell <- s_vec[ell]
    a_index_list[[ell]] <- integer(s_ell)
    cohort_index_list[[ell]] <- integer(s_ell)
    obs_s_list[[ell]] <- integer(s_ell)
    obs_t_list[[ell]] <- integer(s_ell)
    odds_list[[ell]] <- numeric(s_ell)
    score_base_list[[ell]] <- numeric(s_ell)
    offsets_list[[ell]] <- integer(n_chunks + 1)
  }

  nu_cache_env <- new.env(parent = emptyenv())
  nu_cache_env$active <- FALSE
  clear_nuisance_cache <- function() {
    for (nm in ls(nu_cache_env, all.names = TRUE)) {
      if (!ipw_point_view_identical(nm, "active")) {
        assign(nm, NULL, envir = nu_cache_env)
        rm(list = nm, envir = nu_cache_env)
      }
    }
    nu_cache_env$active <- FALSE
  }
  on.exit(clear_nuisance_cache(), add = TRUE)

  chunks_seen <- 0L
  emitted_units <- 0L
  a_start <- 1L
  while (a_start <= na) {
    ipw_owner_check_interrupt()
    a_s <- as.integer(a_start)
    a_e <- as.integer(min(as.double(na), as.double(a_start) + step - 1))
    b_len <- a_e - a_s + 1L
    chunks_seen <- chunks_seen + 1L
    if (chunks_seen > n_chunks) {
      clear_nuisance_cache()
      ipw_point_view_general_abort("General point-IF replay state construction exceeded planned chunk count.")
    }
    a_seq <- seq.int(a_s, a_e)
    u_w <- fo$w_a[a_seq]
    nu_cache_env$active <- FALSE

    for (ell in seq_len(L)) {
      ipw_owner_check_interrupt()
      fit <- norm$links[[ell]]
      lv_ell <- plan$link_views[[ell]]

      can_reuse_cache <- isTRUE(use_nuisance_cache) &&
        isTRUE(nu_cache_env$active) &&
        ipw_point_view_identical(nu_cache_env$fitted_nuisance_ref, fit$fitted_nuisance_ref) &&
        ipw_point_view_identical(nu_cache_env$a_start, a_s) &&
        ipw_point_view_identical(nu_cache_env$a_end, a_e)

      if (!can_reuse_cache && isTRUE(nu_cache_env$active)) {
        clear_nuisance_cache()
      }

      ch <- query_ipw_view_chunk(
        view = lv_ell,
        a_start = a_s,
        a_end = a_e,
        control = owner_ctrl,
        read_outcomes = FALSE,
        include_design = TRUE,
        enclose_differences = FALSE,
        trace_env = NULL,
        validated = norm$validated[[ell]]
      )

      if (isTRUE(ch$read_diagnostics$read_outcomes) ||
          !ipw_point_view_identical(as.integer(ch$read_diagnostics$n_endpoints_read), 0L) ||
          !isTRUE(ch$read_diagnostics$risk_filtered_before_outcome_read) ||
          isTRUE(ch$read_diagnostics$off_risk_endpoints_read) ||
          isTRUE(ch$read_diagnostics$nonobserver_endpoints_read)) {
        ch <- NULL
        clear_nuisance_cache()
        cdid_ipw_point_view_abort(
          "Point IF replay state construction violated outcome-free or risk-first query invariants.",
          c("cdid_ipw_point_view_outcome_leak_error", "cdid_ipw_point_view_general_error", "cdid_ipw_owner_integrity_error")
        )
      }

      k_ch <- as.integer(ch$counts$observed_risk_row_count)
      cur <- cursors[ell]
      s_ell <- s_vec[ell]
      if (is.na(k_ch) || k_ch < 0L || k_ch > b_len || cur + k_ch > s_ell) {
        ch <- NULL
        clear_nuisance_cache()
        cdid_ipw_point_view_abort(
          "General point-IF replay state chunk observer count disagrees with sealed link count.",
          c("cdid_ipw_count_error", "cdid_ipw_point_view_general_error")
        )
      }

      in_risk_chunk <- which(ch$in_risk)
      if (length(in_risk_chunk) > 0L) {
        if (k_ch > 0L) {
          pos_k <- as.integer(ch$observed_risk$local_pos)
          a_k <- as.integer(ch$observed_risk$a_index)
          gj_k <- as.integer(ch$observed_risk$cohort_index)
          obs_s_k <- as.integer(ch$observed_risk$obs_s_index)
          obs_t_k <- as.integer(ch$observed_risk$obs_t_index)
          if (length(pos_k) != k_ch || length(a_k) != k_ch || length(gj_k) != k_ch ||
              length(obs_s_k) != k_ch || length(obs_t_k) != k_ch ||
              anyNA(a_k) || a_k[1L] < a_s || a_k[k_ch] > a_e || is.unsorted(a_k, strictly = TRUE) ||
              !ipw_point_view_identical(pos_k, as.integer(a_k - a_s + 1L)) ||
              anyNA(gj_k) || any(gj_k < 1L) || any(gj_k > as.integer(lv_ell$cohorts$m_cohorts)) ||
              anyNA(obs_s_k) || anyNA(obs_t_k) ||
              any(obs_s_k < 1L) || any(obs_s_k > n_obs_y) ||
              any(obs_t_k < 1L) || any(obs_t_k > n_obs_y)) {
            ch <- NULL
            clear_nuisance_cache()
            ipw_point_view_general_abort("General point-IF replay state chunk query returned malformed observer indices or positions.")
          }

          if (!all(ch$in_risk[pos_k]) ||
              !ipw_point_view_identical(gj_k, as.integer(ch$cohort_index[pos_k]))) {
            ch <- NULL
            clear_nuisance_cache()
            cdid_ipw_point_view_abort(
              "Point IF chunk query violated risk-first endpoint read filtering.",
              c("cdid_ipw_point_view_outcome_leak_error", "cdid_ipw_point_view_general_error", "cdid_ipw_owner_integrity_error")
            )
          }

          if (can_reuse_cache) {
            if (!ipw_point_view_identical(pos_k, nu_cache_env$pos_k) ||
                !ipw_point_view_identical(a_k, nu_cache_env$a_k) ||
                !ipw_point_view_identical(gj_k, nu_cache_env$gj_k)) {
              ch <- NULL
              clear_nuisance_cache()
              cdid_ipw_point_view_abort(
                "Shared-nuisance point link observer mask disagrees with cached nuisance chunk mask.",
                c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error")
              )
            }
          }

          w_k <- u_w[pos_k]
          r_k <- ipw_clip_floor(fit$masses$rates[gj_k])
          D_k <- as.double(gj_k == 1L)
          Z_k <- ch$design$Z_obs
          lv_k <- ipw_logit_values(Z_k, fit$coefficients)
          odds_k <- lv_k$odds
          prob_k <- lv_k$prob
          score_base_k <- (w_k / r_k) * (D_k - prob_k)
          rm(lv_k)
          if (anyNA(odds_k) || any(!is.finite(odds_k)) || anyNA(score_base_k) || any(!is.finite(score_base_k))) {
            ch <- NULL
            clear_nuisance_cache()
            cdid_abort("Nonfinite IF contribution encountered during replay.", c("cdid_ipw_numeric_error"))
          }

          if (!can_reuse_cache && isTRUE(use_nuisance_cache)) {
            nu_cache_env$active <- TRUE
            nu_cache_env$fitted_nuisance_ref <- fit$fitted_nuisance_ref
            nu_cache_env$a_start <- a_s
            nu_cache_env$a_end <- a_e
            nu_cache_env$k_ch <- k_ch
            nu_cache_env$pos_k <- pos_k
            nu_cache_env$a_k <- a_k
            nu_cache_env$gj_k <- gj_k
          }

          dst <- seq.int(cur + 1L, cur + k_ch)
          a_index_list[[ell]][dst] <- a_k
          cohort_index_list[[ell]][dst] <- gj_k
          obs_s_list[[ell]][dst] <- obs_s_k
          obs_t_list[[ell]][dst] <- obs_t_k
          odds_list[[ell]][dst] <- odds_k
          score_base_list[[ell]][dst] <- score_base_k
          cursors[ell] <- cur + k_ch
        } else if (can_reuse_cache) {
          if (!ipw_point_view_identical(nu_cache_env$k_ch, 0L)) {
            ch <- NULL
            clear_nuisance_cache()
            cdid_ipw_point_view_abort(
              "Shared-nuisance point link observer count disagrees with cached nuisance chunk.",
              c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error")
            )
          }
        } else if (isTRUE(use_nuisance_cache)) {
          nu_cache_env$active <- TRUE
          nu_cache_env$fitted_nuisance_ref <- fit$fitted_nuisance_ref
          nu_cache_env$a_start <- a_s
          nu_cache_env$a_end <- a_e
          nu_cache_env$k_ch <- 0L
        }
      } else if (k_ch != 0L) {
        ch <- NULL
        clear_nuisance_cache()
        ipw_point_view_general_abort("General point-IF replay state chunk reported observed risk rows in an empty risk chunk.")
      }

      offsets_list[[ell]][chunks_seen + 1] <- cursors[ell]
      ch <- NULL
    }

    emitted_units <- emitted_units + b_len
    a_start <- as.double(a_e) + 1
  }

  clear_nuisance_cache()

  if (emitted_units != na || chunks_seen != n_chunks || !identical(cursors, s_vec)) {
    cdid_ipw_point_view_abort(
      "General point-IF replay state construction counts disagree with sealed point-link counts.",
      c("cdid_ipw_count_error", "cdid_ipw_point_view_general_error")
    )
  }

  link_states <- vector("list", L)
  link_state_refs <- vector("list", L)
  for (ell in seq_len(L)) {
    pl_ell <- plan$point_links[[ell]]
    lv_ell <- plan$link_views[[ell]]
    ls_ell <- ipw_owner_seal(list(
      kind = "cdid_ipw_point_if_general_link_state",
      link_index = as.integer(ell),
      link_id = plan$link_ids[[ell]],
      point_link_ref = pl_ell$ref,
      link_view_ref = lv_ell$ref,
      fitted_nuisance_ref = pl_ell$fitted_nuisance_ref,
      nuisance_view_ref = pl_ell$nuisance_view_ref,
      n_observed_risk_units = s_vec[ell],
      n_chunks = n_chunks,
      chunk_offsets = offsets_list[[ell]],
      a_index = a_index_list[[ell]],
      cohort_index = cohort_index_list[[ell]],
      obs_s_index = obs_s_list[[ell]],
      obs_t_index = obs_t_list[[ell]],
      odds = odds_list[[ell]],
      score_base = score_base_list[[ell]]
    ), "point_if_general_link_state")
    link_states[[ell]] <- ls_ell
    link_state_refs[[ell]] <- ls_ell$ref
  }

  binding_ref <- new.env(parent = emptyenv())
  binding_ref$plan_snapshot <- plan
  lockEnvironment(binding_ref, bindings = TRUE)

  runtime <- new.env(parent = emptyenv())
  runtime$status <- "built"
  runtime$passes_completed <- 0L
  runtime$chunks_evaluated <- 0L
  runtime$verifications <- 0L
  lockEnvironment(runtime, bindings = FALSE)

  state <- ipw_owner_seal(list(
    kind = "cdid_ipw_point_if_general_state",
    plan_ref = plan$ref,
    frame_ref = fo$ref,
    panel_ref = fo$panel_ref,
    binding_ref = binding_ref,
    runtime = runtime,
    n_a_units = na,
    n_links = L,
    step = step,
    n_chunks = n_chunks,
    total_observed_risk_units = reservation$S,
    raw_entries = reservation$raw_entries,
    entries = reservation$entries,
    work = reservation$work,
    construction_work_ops_bound = reservation$construction_work_ops_bound,
    copy_slice_work_ops_bound = reservation$copy_slice_work_ops_bound,
    one_state_check_work_ops_bound = reservation$one_state_check_work_ops_bound,
    state_checks_planned = reservation$state_checks_planned,
    state_check_work_ops_bound = reservation$state_check_work_ops_bound,
    chunk_guard_work_ops_bound = reservation$chunk_guard_work_ops_bound,
    setup_work_ops_bound = reservation$setup_work_ops_bound,
    admission_total = reservation$admission_total,
    work_total = reservation$work_total,
    cap = reservation$cap,
    work_cap = reservation$work_cap,
    queries_performed = logit_cert_count_output(reservation$queries_expected),
    endpoints_read_during_construction = 0L,
    link_ids = plan$link_ids,
    link_states = link_states,
    link_state_refs = link_state_refs
  ), "point_if_general_state")

  reservation$state <- "built"
  state
}

ipw_point_if_general_state <- function(
  plan,
  control,
  live = 0,
  work = 0,
  envelope_control = NULL,
  use_nuisance_cache = TRUE
) {
  res <- ipw_point_if_general_state_reserve(
    plan = plan,
    control = control,
    live = live,
    work = work,
    envelope_control = envelope_control
  )
  if (is.null(res)) return(NULL)
  ipw_point_if_general_state_build(
    reservation = res,
    plan = plan,
    control = control,
    use_nuisance_cache = use_nuisance_cache
  )
}

ipw_point_if_general_state_release <- function(state, plan = NULL, superseded = FALSE) {
  if (is.null(state)) {
    null_status <- if (isTRUE(superseded)) "superseded_by_retained_if" else "declined_capacity"
    if (is.list(plan) && identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
      s_sum <- sum(vapply(plan$point_links, function(lk) as.double(lk$counts$n_observed_risk_units), numeric(1L)))
      c_num <- as.double(plan$preflight$n_chunks)
      l_num <- as.double(plan$n_links)
      b_num <- as.double(plan$step)
      ent <- ipw_point_if_general_state_entries(s_sum, c_num, l_num, b_num)
      wrk <- ipw_point_if_general_state_work(plan$preflight$total_work_ops_bound,
        ent$raw_entries, c_num, l_num,
        plan_check_work = plan$general$one_verification_work_ops_bound)
      raw_e <- ent$raw_entries
      res_e <- ent$reserved_entries
      chk_w <- wrk$state_check_work
      tot_w <- wrk$total_work
      return(list(
        status = null_status,
        admitted = FALSE,
        released = FALSE,
        n_a_units = plan$n_a_units,
        n_links = plan$n_links,
        chunk_rows = plan$step,
        n_chunks = as.integer(plan$preflight$n_chunks),
        total_observed_risk_units = s_sum,
        raw_state_entries = raw_e,
        reserved_entries = res_e,
        construction_work_ops_bound = plan$preflight$total_work_ops_bound,
        copy_slice_work_ops_bound = wrk$copy_slice_work,
        one_state_check_work_ops_bound = wrk$one_state_check,
        state_checks_planned = wrk$state_checks_planned,
        state_check_work_ops_bound = chk_w,
        chunk_guard_work_ops_bound = wrk$chunk_guard_work,
        setup_work_ops_bound = wrk$setup_work,
        work_ops_bound = tot_w,
        queries_performed = 0L,
        endpoints_read_during_construction = 0L,
        passes_served = 0L,
        chunks_served = 0L,
        state_verifications = 0L,
        outcomes_retained = FALSE,
        differences_retained = FALSE,
        full_n_by_L_retained = FALSE,
        per_link_dense_masks_retained = FALSE,
        peak_rss_promised = FALSE
      ))
    }
    return(list(status = null_status, admitted = FALSE, released = FALSE))
  }
  rt <- state$runtime
  was_built <- is.environment(rt) && identical(rt$status, "built")
  rec <- list(
    status = if (was_built) "retained_and_released" else paste0("incomplete:", if (is.environment(rt)) rt$status else "unknown"),
    admitted = TRUE,
    released = TRUE,
    n_a_units = state$n_a_units,
    n_links = state$n_links,
    chunk_rows = state$step,
    n_chunks = state$n_chunks,
    total_observed_risk_units = state$total_observed_risk_units,
    raw_state_entries = state$raw_entries,
    reserved_entries = state$entries,
    admission_total_entries = state$admission_total,
    cap_entries = state$cap,
    construction_work_ops_bound = state$construction_work_ops_bound,
    copy_slice_work_ops_bound = state$copy_slice_work_ops_bound,
    one_state_check_work_ops_bound = state$one_state_check_work_ops_bound,
    state_checks_planned = state$state_checks_planned,
    state_check_work_ops_bound = state$state_check_work_ops_bound,
    chunk_guard_work_ops_bound = state$chunk_guard_work_ops_bound,
    setup_work_ops_bound = state$setup_work_ops_bound,
    work_ops_bound = state$work,
    admission_total_work_ops = state$work_total,
    cap_work_ops = state$work_cap,
    queries_performed = state$queries_performed,
    endpoints_read_during_construction = state$endpoints_read_during_construction,
    passes_served = if (is.environment(rt)) logit_cert_count_output(rt$passes_completed) else 0L,
    chunks_served = if (is.environment(rt)) logit_cert_count_output(rt$chunks_evaluated) else 0L,
    state_verifications = if (is.environment(rt)) logit_cert_count_output(rt$verifications) else 0L,
    outcomes_retained = FALSE,
    differences_retained = FALSE,
    full_n_by_L_retained = FALSE,
    per_link_dense_masks_retained = FALSE,
    peak_rss_promised = FALSE,
    lifetime = paste0(
      "Admitted against the complete simultaneous pipeline envelope before allocation; built once from ",
      "query_ipw_view_chunk(read_outcomes = FALSE, include_design = TRUE) in chunk/link order; sealed before ",
      "callbacks; sliced across the two covariance passes and two GMM projection passes while reading raw ",
      "endpoints through ipw_owner_read_y_positions; released immediately after solve_chained_gmm_core."
    )
  )
  if (is.environment(rt)) rt$status <- "released"
  if (is.list(state$link_state_refs)) {
    for (lref in state$link_state_refs) {
      if (is.environment(lref) && exists("point_if_general_link_state", envir = lref, inherits = FALSE)) {
        if (bindingIsLocked("point_if_general_link_state", lref)) unlockBinding("point_if_general_link_state", lref)
        assign("point_if_general_link_state", NULL, envir = lref)
      }
    }
  }
  if (is.environment(state$binding_ref) && exists("plan_snapshot", envir = state$binding_ref, inherits = FALSE)) {
    if (bindingIsLocked("plan_snapshot", state$binding_ref)) unlockBinding("plan_snapshot", state$binding_ref)
    assign("plan_snapshot", NULL, envir = state$binding_ref)
  }
  if (is.environment(state$ref) && exists("point_if_general_state", envir = state$ref, inherits = FALSE)) {
    if (bindingIsLocked("point_if_general_state", state$ref)) unlockBinding("point_if_general_state", state$ref)
    assign("point_if_general_state", NULL, envir = state$ref)
  }
  rec
}

# The general replay on a verified record: the same chunk partition, the same
# `ipw_point_if_eval_chunk_internal` chunk evaluation (risk-first endpoint
# reads, the four signed IF terms in their order, the one-nuisance chunk
# cache), the same callback order and refusals as `replay_ipw_point_if_view`,
# without that wrapper's per-call normalization. The record (and any admitted
# sparse replay state bound to it) is re-verified at the pass boundary and
# after every callback, including the last callback, before the next chunk's
# endpoint read or return.
replay_ipw_point_if_general <- function(
  plan,
  control,
  callback,
  callback_schema = TRUE,
  use_nuisance_cache = TRUE,
  trace_env = NULL,
  state = NULL,
  carrier = NULL
) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control) || missing(callback)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`plan`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.function(callback)) {
    cdid_ipw_point_view_abort("`callback` must be a function.", c("cdid_ipw_point_view_input_error", "cdid_ipw_input_error"))
  }
  if (!gmm_scalar_logical(callback_schema) || !gmm_scalar_logical(use_nuisance_cache)) {
    cdid_ipw_point_view_abort("`callback_schema` and `use_nuisance_cache` must be scalar logicals.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(carrier)) {
    if (!is.null(state)) {
      ipw_point_view_general_abort("General point-IF replay cannot use both sparse replay state and retained IF carrier simultaneously.")
    }
    if (!isTRUE(callback_schema)) {
      cdid_ipw_point_view_abort(
        "Retained point IF carrier supports only the callback-schema IF replay.",
        c("cdid_ipw_point_view_input_error", "cdid_ipw_point_view_general_error")
      )
    }
    return(replay_ipw_point_if_retained(carrier, plan, control, callback, trace_env = trace_env))
  }
  ipw_point_if_general_verify(plan, control)
  if (!is.null(state)) {
    ipw_point_if_general_state_verify(state, plan, control)
  }
  verifications <- 1L
  pre <- plan$preflight
  ctrl <- plan$control
  norm <- plan$norm

  na <- as.integer(pre$n_a_units)
  L <- as.integer(pre$n_links)
  step <- as.integer(pre$effective_chunk_rows)

  if (!is.null(trace_env)) {
    trace_env$preflight_completed <- TRUE
    trace_env$batches_emitted <- 0L
    trace_env$units_emitted <- 0L
    trace_env$active_query_chunks <- 0L
    trace_env$max_simultaneous_query_chunks <- 0L
    trace_env$active_nuisance_caches <- 0L
    trace_env$max_simultaneous_nuisance_caches <- 0L
    trace_env$nuisance_cache_hits <- 0L
    trace_env$endpoints_read <- 0L
  }

  emitted_units <- 0L
  batches <- 0L
  max_batch_units <- 0L
  total_work_ops <- 0

  a_start <- 1L
  while (a_start <= na) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(na), as.double(a_start) + step - 1))
    b_len <- a_end - a_start + 1L
    batches <- batches + 1L

    batch <- ipw_point_if_eval_chunk_internal(
      norm = norm,
      a_start = a_start,
      a_end = a_end,
      ctrl = ctrl,
      include_components = FALSE,
      callback_schema = callback_schema,
      use_nuisance_cache = use_nuisance_cache,
      batch_index = batches,
      trace_env = trace_env,
      state = state
    )
    if (!isTRUE(callback_schema) && !is.null(batch$work$actual_work_ops)) {
      total_work_ops <- total_work_ops + batch$work$actual_work_ops
    }
    callback(batch)
    batch <- NULL
    ipw_point_if_general_verify(plan, control)
    if (!is.null(state)) {
      ipw_point_if_general_state_verify(state, plan, control)
    }
    verifications <- verifications + 1L

    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    if (!is.null(trace_env)) {
      trace_env$batches_emitted <- batches
      trace_env$units_emitted <- emitted_units
    }
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != na) {
    cdid_ipw_point_view_abort("Point IF replay unit count mismatch.", c("cdid_ipw_count_error"))
  }
  if (!is.null(state) && is.environment(state$runtime)) {
    state$runtime$passes_completed <- logit_cert_count_output(as.double(state$runtime$passes_completed) + 1)
  }

  invisible(list(
    n_a_units = emitted_units,
    n_links = L,
    link_ids = norm$link_ids,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * as.double(L),
    actual_work_ops = total_work_ops,
    normalizations_performed = 0L,
    verifications = verifications,
    state_verifications = if (is.null(state)) 0L else verifications,
    served_from = if (is.null(state)) "general_chunk_query" else "general_sparse_replay_state",
    preflight = pre
  ))
}

# Two-pass centered point moment covariance on a verified record: the same
# reservations, passes, nonfinite refusals, centering and scaling as
# `stream_ipw_point_covariance_view`, with the record's one normalization in
# place of that wrapper's three, and optional admitted sparse replay state.
stream_ipw_point_covariance_general <- function(
  plan,
  control,
  use_nuisance_cache = TRUE,
  trace_env = NULL,
  state = NULL,
  carrier = NULL
) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`plan`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!gmm_scalar_logical(use_nuisance_cache)) {
    cdid_ipw_point_view_abort("`use_nuisance_cache` must be a scalar logical.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  ipw_point_if_general_verify(plan, control)
  if (!is.null(state) && !is.null(carrier)) {
    ipw_point_view_general_abort("General point-IF covariance stream cannot use both sparse replay state and retained IF carrier simultaneously.")
  }
  state_entries <- 0
  state_work <- 0
  if (!is.null(state)) {
    ipw_point_if_general_state_verify(state, plan, control)
    state_entries <- state$entries
    state_work <- state$work
  } else if (!is.null(carrier)) {
    ipw_point_if_retained_check(carrier, plan, control, verify_plan = FALSE)
    if (!identical(carrier$state, "reserved")) {
      ipw_point_view_general_abort("Retained point IF carrier must be reserved and unfilled for the covariance stream.")
    }
    state_entries <- carrier$entries
    state_work <- carrier$work
  }
  pre1 <- plan$preflight
  ctrl <- plan$control
  norm <- plan$norm
  n_a <- as.double(pre1$n_a_units)
  L <- as.integer(pre1$n_links)
  link_ids <- norm$link_ids

  cov_entries <- ipw_owner_check_matrix_dims(
    L, L,
    max_matrix_entries = min(ctrl$max_matrix_entries, ctrl$max_output_entries),
    label = "streamed L x L point IF covariance matrix"
  )
  ipw_owner_safe_sum(
    c(pre1$peak_concurrent_entries, 4 * cov_entries, 8 * L, 6 * pre1$b_by_L_entries,
      plan$general$additional_live_entries, state_entries, 128),
    ctrl$max_matrix_entries,
    "two-pass streamed point covariance concurrent entries plus retained normalization record"
  )
  ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(2, pre1$total_work_ops_bound, ctrl$max_work_ops, "two-pass point IF replay work"),
      ipw_owner_safe_product(1 + 2 * (as.double(pre1$n_chunks) + 1),
        plan$general$one_verification_work_ops_bound, ctrl$max_work_ops,
        "general covariance entry and two replay record verifications"),
      ipw_owner_safe_product(n_a, 2 * L * L + 4 * L, ctrl$max_work_ops, "centering and crossprod accumulation work"),
      state_work,
      128
    ),
    ctrl$max_work_ops,
    "total two-pass streamed point covariance work"
  )

  sum_phi <- numeric(L)
  mean_callback <- function(batch) {
    sum_phi <<- sum_phi + colSums(batch$phi)
    if (any(!is.finite(sum_phi))) {
      cdid_abort("Nonfinite IF sum during covariance replay.", "cdid_ipw_numeric_error")
    }
  }
  if (is.null(carrier)) {
    replay_ipw_point_if_general(
      plan = plan,
      control = ctrl,
      callback_schema = TRUE,
      use_nuisance_cache = use_nuisance_cache,
      trace_env = trace_env,
      state = state,
      callback = mean_callback
    )
  } else {
    ipw_point_if_retained_fill(
      carrier = carrier,
      plan = plan,
      control = ctrl,
      callback = mean_callback,
      use_nuisance_cache = use_nuisance_cache,
      trace_env = trace_env,
      .verified_entry = TRUE
    )
  }
  mean_phi <- sum_phi / n_a
  names(mean_phi) <- link_ids

  cov_mat <- matrix(0.0, nrow = L, ncol = L, dimnames = list(link_ids, link_ids))
  replay_info <- replay_ipw_point_if_general(
    plan = plan,
    control = ctrl,
    callback_schema = TRUE,
    use_nuisance_cache = use_nuisance_cache,
    trace_env = trace_env,
    state = state,
    carrier = carrier,
    callback = function(batch) {
      centered <- sweep(batch$phi, 2L, mean_phi, FUN = "-")
      if (any(!is.finite(centered))) {
        cdid_abort("Nonfinite centered IF during covariance replay.", "cdid_ipw_numeric_error")
      }
      cov_mat <<- cov_mat + crossprod(centered)
      if (any(!is.finite(cov_mat))) {
        cdid_abort("Nonfinite IF crossproduct or sum during covariance replay.", "cdid_ipw_numeric_error")
      }
      rm(centered)
    }
  )
  cov_mat <- cov_mat / n_a

  list(
    covariance = cov_mat,
    mean_phi = mean_phi,
    estimates = stats::setNames(
      vapply(norm$links, `[[`, numeric(1L), "estimate"),
      link_ids
    ),
    n_a_units = pre1$n_a_units,
    replay = replay_info
  )
}

# ---------------------------------------------------------------------------
# Section 2b: Trusted Batch Execution Plan (CP08 balanced performance, stage 2)
# ---------------------------------------------------------------------------
# One invocation-local route for the exact eligible baseline IPW subset, which
# `run_chained_ipw_pipeline` selects from actual facts before any retained
# point fit (`ipw_pipeline_batch_eligibility`): the outcome-free optimal
# reduction of built-in IPW under PT-post without numerical transfer, one
# fitted nuisance shared by every retained link, a single target cohort with
# the never-treated pool as its only control cohort, the intercept-only design,
# constant unit weights and complete calendar dates (every A-unit observed) at
# every required date. Every other request keeps the general route; a failure
# here propagates and is never retried through another route.
#   * `fit_ipw_point_link_view_batch` fits the retained point links from one
#     shared observer state. Each link still reads its own raw endpoints
#     (y_t - y_s, canonical observer order) and calls the unchanged
#     `ipw_link_means` and `ipw_link_partials`; there is no propensity refit.
#   * `ipw_point_if_batch_plan` performs the single normalization/preflight and
#     seals compact link descriptors. The two covariance passes and the two GMM
#     projection passes replay through it; each regenerates every IF from raw
#     endpoints with the registered kernel `C_gmm_ipw_if_chunk`, one call per
#     bounded chunk, and re-verifies the bound objects at its boundary and
#     after each callback. The plan is not a caller certification token:
#     consumers check its seal, its controls and the live owner/view bindings.

# Kernel resolution at each chunk boundary. The existing resolver first checks
# the current `cdid` DLL (identity, its nine routines and arities, lifecycle,
# live arithmetic); this kernel's registered symbol and arity are then resolved
# from that same DLL. The address is reused only while the existing resolver
# returns the very same symbol table object, so any re-resolution there
# re-resolves here. The kernel repeats the live arithmetic check on every call.
# No R fallback exists: a backend without the kernel refuses with a typed error.
ipw_point_if_native_cache <- new.env(parent = emptyenv())

ipw_point_if_native_kernel <- function() {
  syms <- gmm_require_native_dll()
  owner <- environment(sys.function())
  cache <- get0("ipw_point_if_native_cache", envir = owner, mode = "environment", inherits = FALSE)
  if (is.environment(cache) && identical(cache$owner, owner) && !is.null(cache$address) &&
      gmm_sexp_aliases(cache$syms, syms)) return(cache$address)
  dlls <- if (isNamespace(owner)) getNamespaceInfo(owner, "DLLs") else base::getLoadedDLLs()
  dll <- dlls[["cdid"]]
  si <- tryCatch(base::getNativeSymbolInfo("C_gmm_ipw_if_chunk", PACKAGE = dll), error = function(e) NULL)
  if (is.null(si) || !inherits(si, "NativeSymbolInfo") || is.null(si$address) ||
      !identical(as.integer(si$numParameters), 18L)) {
    cdid_gmm_abort(
      "Native symbol registration/forcing verification failed for `C_gmm_ipw_if_chunk`.",
      c("cdid_gmm_native_error", "cdid_gmm_accuracy_error")
    )
  }
  if (is.environment(cache)) {
    cache$owner <- owner
    cache$syms <- syms
    cache$address <- si$address
  }
  si$address
}

ipw_point_view_batch_abort <- function(message, diagnostics = list()) {
  cdid_ipw_point_view_abort(message,
    c("cdid_ipw_point_view_batch_error", "cdid_ipw_point_view_integrity_error", "cdid_ipw_owner_integrity_error"),
    diagnostics = diagnostics)
}

# Retained point links fitted from one shared observer state. Requires one
# fitted nuisance, link views bound to its owners and nuisance snapshot with
# identical cohort descriptors, and complete dates for the defining pair and
# every link pair, so every link's observer set is every risk A-unit in
# canonical order. Validation runs once for the nuisance and once per link
# before any allocation; a single preflight charges the shared observer state,
# one link's endpoint/difference/helper workspace, every retained output and
# the referenced owners. The shared state lives across this loop only and is
# released before returning; each link records that it was live at sealing.
fit_ipw_point_link_view_batch <- function(
  fitted_view,
  link_views,
  control,
  date_memo = NULL,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(fitted_view) || missing(link_views) || missing(control)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`fitted_view`, `link_views`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.list(link_views) || is.object(link_views) || length(link_views) < 1L) {
    cdid_ipw_point_view_abort(
      "`link_views` must be a non-empty list of sealed `cdid_ipw_link_view` objects.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(date_memo) && !is.environment(date_memo)) {
    cdid_ipw_point_view_abort("`date_memo` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  ctrl <- validate_ipw_point_view_control(control)
  L <- length(link_views)
  ipw_owner_safe_product(64, L, min(ctrl$max_retained_entries, ctrl$max_matrix_entries, ctrl$max_output_entries),
    "batch point-link list header before validation")
  ipw_point_view_seal_gate(fitted_view, "fitted_nuisance_view")
  for (lv in link_views) ipw_point_view_seal_gate(lv, "link_view")
  payload <- ipw_point_view_payload_budget(list(fitted_view, link_views), ctrl)

  # 1. Validate the nuisance and its owners once, then every link view.
  validate_ipw_fitted_nuisance_view(fitted_view)
  frame_owner <- get("frame_owner", envir = fitted_view$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = fitted_view$panel_ref, inherits = FALSE)
  nv_obj <- get("nuisance_view", envir = fitted_view$nuisance_view_ref, inherits = FALSE)
  validate_ipw_frame_owner(frame_owner, panel_owner = panel_owner)
  validate_ipw_glm_control(fitted_view$glm_control)
  na <- ipw_owner_int_bound(frame_owner$n_a_units, "n_a_units", min_val = 1)
  n_obs <- ipw_owner_int_bound(fitted_view$counts$n_observed_risk_units, "n_observed_risk_units", min_val = 0)
  p <- ipw_owner_int_bound(frame_owner$n_covariates, "n_covariates", min_val = 1)
  m <- ipw_owner_int_bound(fitted_view$cohorts$m_cohorts, "m_cohorts", min_val = 2)
  # Stored and helper-facing counts keep the scalar fitter's integer types.
  na_i <- as.integer(na)
  n_obs_i <- as.integer(n_obs)
  p_i <- as.integer(p)
  m_i <- as.integer(m)
  coh <- nv_obj$cohorts
  step <- min(na, ctrl$max_chunk_rows)
  memo <- if (is.environment(date_memo) && identical(date_memo$panel_owner, panel_owner)) date_memo else
    list2env(list(panel_owner = panel_owner), parent = emptyenv())
  if (!identical(panel_owner$n_a_units, frame_owner$n_a_units) ||
      !ipw_owner_pair_complete(panel_owner, nv_obj$defining_pair, step, memo)) {
    ipw_point_view_batch_abort("Batch point fitting requires the nuisance's defining pair at complete dates.")
  }
  for (ell in seq_len(L)) {
    lv <- link_views[[ell]]
    validate_ipw_link_view(lv, frame_owner = frame_owner, nuisance_view = nv_obj)
    if (!identical(lv$frame_ref, fitted_view$frame_ref) || !identical(lv$panel_ref, fitted_view$panel_ref)) {
      cdid_ipw_point_view_abort(
        "Link view belongs to a different frame or panel owner than the fitted nuisance view; reuse is refused.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_context_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
    if (!identical(lv$cohorts, coh) || !identical(fitted_view$cohorts$cohorts, coh$cohorts) ||
        !identical(fitted_view$cohorts$cohort_labels, coh$cohort_labels)) {
      cdid_ipw_point_view_abort(
        "Link view `g` or ordered control cohorts `K` differ from the fitted nuisance view; reuse is refused.",
        c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_mismatch_error", "cdid_ipw_shared_incompatible_error")
      )
    }
    if (!ipw_point_view_identical(as.integer(lv$support$n_observed_risk_units), as.integer(n_obs)) ||
        !ipw_point_view_identical(as.integer(lv$support$by_cohort_pair_count), as.integer(unname(fitted_view$counts$by_cohort_pair)))) {
      cdid_ipw_point_view_abort(
        "Link view observed risk-set support counts differ from the fitted nuisance view.",
        c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error", "cdid_ipw_shared_incompatible_error")
      )
    }
    # Mask equality by the exact posting fact (four complete dates), as the
    # existing verification decides it; this route requires it for every link.
    if (!ipw_owner_pair_complete(panel_owner, lv$pair, step, memo) ||
        !ipw_owner_pair_masks_equal_fact(panel_owner, nv_obj$defining_pair, lv$pair, step, memo)) {
      ipw_point_view_batch_abort("Batch point fitting requires every link pair at complete dates.",
        diagnostics = list(link = ell))
    }
  }
  if (n_obs > ctrl$max_observer_units) {
    cdid_ipw_point_view_abort(
      paste0("Observed risk-set pair unit count (", n_obs, ") exceeds `control$max_observer_units` (", ctrl$max_observer_units, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error"),
      diagnostics = list(n_observed_risk_units = n_obs, max_observer_units = ctrl$max_observer_units)
    )
  }
  if (n_obs < p || any(nv_obj$support$by_cohort_pair_count < 1L) ||
      any(!is.finite(nv_obj$support$by_cohort_pair_mass) | nv_obj$support$by_cohort_pair_mass <= 0)) {
    cdid_ipw_point_view_abort(
      "Link view lacks positive observed risk-set support across all declared cohorts.",
      c("cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }

  # 2. One preflight before any observer allocation or endpoint read.
  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  override_budget <- ipw_point_view_payload_budget(list(NULL, NULL), ctrl)
  char_meta <- lapply(unname(link_views), function(lv) ipw_point_view_char_bytes(list(
    colnames(frame_owner$Z_a), coh$cohort_labels, fitted_view$nuisance_id, lv$link_id, lv$functional_key), ctrl))
  id_meta <- ipw_point_view_id_bytes_bound(panel_owner, step, ceiling(na / step), ctrl)
  p2 <- ipw_owner_check_matrix_dims(p, p, max_matrix_entries = max_ret, label = "point-link p x p factor matrix")
  pm <- ipw_owner_check_matrix_dims(p, m, max_matrix_entries = max_ret, label = "point-link p x m score-rate Jacobian")
  # Per link, the scalar fitter's compact retained envelope (same terms).
  ret_matrix <- vapply(unname(link_views), function(lv) ipw_owner_safe_sum(c(
    ipw_owner_safe_product(3, p2, max_ret, "H, H_chol, and qr_R matrices"), pm,
    lv$preflight$retained_matrix_entries), max_ret, "retained point-link matrix entries"), numeric(1L))
  ret_vector <- vapply(seq_len(L), function(ell) ipw_owner_safe_sum(c(
    ipw_owner_safe_product(5, p, max_ret, "beta coefficients, qr_pivot, G_beta, and loading_beta"),
    ipw_owner_safe_product(18, m, max_ret, "G_r, loading_rate, rate_coef_s0, rate_coef_s1, counts, masses, and disclosures"),
    char_meta[[ell]]$charged_slots, link_views[[ell]]$preflight$retained_vector_entries, 2048, 8192),
    max_ret, "retained point-link vector and metadata entries"), numeric(1L))
  ret_one <- vapply(seq_len(L), function(ell) ipw_owner_safe_sum(c(ret_matrix[[ell]], ret_vector[[ell]],
    override_budget$entries), max_ret, "single point-link retained entries"), numeric(1L))
  total_retained <- ipw_owner_safe_sum(ret_one, min(max_ret, max_out), "all batch point-link retained outputs")
  obs_Z <- ipw_owner_check_matrix_dims(n_obs, p, max_matrix_entries = max_me, label = "batch shared observer design Z_use")
  shared_state <- ipw_owner_safe_sum(c(obs_Z,
    ipw_owner_safe_product(32, max(1, n_obs), max_me, "shared observer index/cohort/weight/rate/fitting-weight/logit/normalizer vectors and their temporaries"),
    256), max_me, "batch shared observer state")
  build_scratch <- ipw_owner_safe_sum(c(ipw_owner_safe_product(16, step, max_me, "chunked observer-state construction temporaries"),
    ipw_owner_safe_product(step, p, max_me, "chunk design rows"), 64), max_me, "batch observer-state construction scratch")
  link_workspace <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(4, obs_Z, max_me, "control Z_c and weighted Z_c copies in ipw_link_partials"),
    ipw_owner_safe_product(4, p2, max_me, "H_chol and transpose triangular solve workspace"),
    ipw_owner_safe_product(2, pm, max_me, "A_beta_r crossprod workspace"),
    ipw_owner_safe_product(64, max(1, n_obs), max_me, "posting slots, endpoints, raw difference buffer, helper subsets and expression copies"),
    ipw_owner_safe_product(8, p, max_me, "G_beta and triangular solve working vectors"),
    ipw_owner_safe_product(12, m, max_me, "G_r, loading_rate, and rate_coef working vectors"),
    max(vapply(char_meta, `[[`, numeric(1L), "charged_slots")), id_meta$chunk_id_slots, 256),
    max_me, "one link's endpoint, difference and helper workspace")
  failure_output <- ipw_owner_safe_sum(c(ipw_owner_safe_product(5, max(1, n_obs)), pm, ipw_owner_safe_product(2, p),
    ipw_owner_safe_product(4, m), 256), min(max_ret, max_out), "possible failed point-link diagnostic output")
  referenced <- ipw_owner_safe_sum(c(payload$entries, ipw_point_view_owner_entries(frame_owner, panel_owner)),
    max_me, "referenced fitted nuisance, link views and owners")
  peak <- ipw_owner_safe_sum(c(referenced, total_retained, shared_state,
    max(build_scratch, ipw_owner_safe_sum(c(link_workspace, failure_output)))), max_me,
    "batch point-link concurrent entries (shared state, one link workspace, all outputs)")
  peak_matrix <- ipw_owner_safe_sum(c(sum(ret_matrix), obs_Z, ipw_owner_safe_product(4, obs_Z)), max_me,
    "batch point-link concurrent matrix entries")
  shared_work <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(na, 32 + 2 * p, max_wo, "chunked shared observer-state construction work"),
    ipw_owner_safe_product(16, max(1, n_obs) * p, max_wo, "ipw_logit_values and pair-rate work"),
    ipw_owner_safe_product(16, max(1, n_obs) * m, max_wo, "ipw_pair_rates and normalizer work"), 256),
    max_wo, "batch shared observer-state work")
  link_work <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(16, max(1, n_obs), max_wo, "posting-slot checks, endpoint reads and raw differences"),
    ipw_owner_safe_product(16, max(1, n_obs) * p, max_wo, "ipw_link_partials G_beta work"),
    ipw_owner_safe_product(16, max(1, n_obs) * m, max_wo, "ipw_link_means and G_r work"),
    ipw_owner_safe_product(4, p * p + p * m, max_wo, "Cholesky triangular solves and A_beta_r crossprod work"),
    256), max_wo, "single batch link work")
  total_work <- ipw_owner_safe_sum(c(shared_work, ipw_owner_safe_product(L, link_work, max_wo, "all batch link work"),
    vapply(char_meta, `[[`, numeric(1L), "work_ops_bound"), id_meta$stream_id_work_bound,
    ipw_owner_safe_product(8, payload$work, max_wo, "validation and compact result construction work"),
    ipw_owner_safe_product(L, override_budget$work + 4096, max_wo, "per-link sealing and payload checks")),
    max_wo, "total batch point-link work operations")
  batch_record <- list(
    route = "batch_shared_observer_state",
    n_a_units = na,
    n_observed_risk_units = n_obs,
    p_covariates = p,
    m_cohorts = m,
    n_links = L,
    effective_chunk_rows = logit_cert_count_output(step),
    shared_state_entries = shared_state,
    construction_scratch_entries = build_scratch,
    link_workspace_entries = link_workspace,
    failure_output_entries = failure_output,
    referenced_entries = referenced,
    total_retained_entries = total_retained,
    peak_concurrent_entries = peak,
    peak_concurrent_matrix_entries = peak_matrix,
    total_work_ops_bound = total_work,
    propensity_refit_performed = FALSE,
    peak_rss_promised = FALSE
  )

  state <- new.env(parent = emptyenv())
  diag_env <- new.env(parent = emptyenv())
  release_state <- function() {
    for (nm in ls(state, all.names = TRUE)) {
      assign(nm, NULL, envir = state)
      rm(list = nm, envir = state)
    }
    if (exists("error_diagnostics", envir = diag_env, inherits = FALSE)) {
      rm(list = "error_diagnostics", envir = diag_env)
    }
    if (!is.null(trace_env)) {
      trace_env$shared_state_active <- FALSE
      trace_env$shared_state_released <- TRUE
    }
  }
  on.exit(release_state(), add = TRUE)
  if (!is.null(trace_env)) {
    trace_env$preflight <- batch_record
    trace_env$preflight_completed <- TRUE
    trace_env$links_fitted <- 0L
    trace_env$endpoints_read <- 0
    trace_env$propensity_solver_calls <- 0L
    trace_env$shared_state_active <- TRUE
    trace_env$shared_state_released <- FALSE
  }

  # 3. Shared observer state: every risk A-unit in canonical order, one chunk
  #    of cohort indices at a time (outcome-free), filled in place.
  risk_a_idx <- integer(n_obs)
  risk_gi <- integer(n_obs)
  risk_w <- numeric(n_obs)
  Z_use <- matrix(0.0, nrow = n_obs, ncol = p, dimnames = list(NULL, colnames(frame_owner$Z_a)))
  write_pos <- 0L
  cur_s <- 1
  while (cur_s <= na) {
    ipw_owner_check_interrupt()
    cur_e <- min(na, cur_s + step - 1)
    a_seq <- seq.int(as.integer(cur_s), as.integer(cur_e))
    gi <- ipw_owner_view_cohort_index(coh, panel_owner$a_cohort_index[a_seq])
    loc <- which(!is.na(gi))
    k_ch <- length(loc)
    if (k_ch > 0L) {
      if (write_pos + k_ch > n_obs) {
        cdid_ipw_point_view_abort(
          "Streamed risk-set observer count exceeds `link_view$support$n_observed_risk_units`.",
          c("cdid_ipw_point_view_count_error", "cdid_ipw_owner_integrity_error")
        )
      }
      dst <- seq.int(write_pos + 1L, write_pos + k_ch)
      a_ch <- a_seq[loc]
      risk_a_idx[dst] <- a_ch
      risk_gi[dst] <- gi[loc]
      risk_w[dst] <- as.double(frame_owner$w_a[a_ch])
      Z_use[dst, ] <- frame_owner$Z_a[a_ch, , drop = FALSE]
      write_pos <- write_pos + k_ch
      rm(dst, a_ch)
    }
    rm(a_seq, gi, loc)
    cur_s <- cur_e + 1
  }
  if (write_pos != n_obs) {
    cdid_ipw_point_view_abort(
      "Streamed risk-set observer count or endpoint read count disagrees with `link_view$support$n_observed_risk_units`.",
      c("cdid_ipw_point_view_count_error", "cdid_ipw_owner_integrity_error")
    )
  }
  state$risk_a_idx <- risk_a_idx
  state$risk_gi <- risk_gi
  state$risk_w <- risk_w
  state$Z_use <- Z_use
  rm(risk_a_idx, risk_gi, risk_w, Z_use)

  sup <- list(
    n = as.double(na),
    g = fitted_view$cohorts$g,
    K = fitted_view$cohorts$K,
    cohorts = fitted_view$cohorts$cohorts,
    cohort_labels = fitted_view$cohorts$cohort_labels,
    m_cohorts = m_i,
    group_count = as.integer(unname(fitted_view$counts$by_cohort_a)),
    group_a_count = as.integer(unname(fitted_view$counts$by_cohort_a)),
    group_mass = as.double(unname(fitted_view$masses$cohort_a_weight)),
    group_a_mass = as.double(unname(fitted_view$masses$cohort_a_weight)),
    q = as.double(unname(fitted_view$masses$q)),
    a_weight_total = as.double(fitted_view$masses$a_weight_total)
  )
  fail <- function(cnd) {
    err_diag <- if (exists("error_diagnostics", envir = diag_env, inherits = FALSE)) {
      diag_env$error_diagnostics
    } else {
      list(stage = "unknown_point_link_failure")
    }
    cnd$error_diagnostics <- err_diag
    release_state()
    stop(cnd)
  }
  post <- panel_owner$postings
  oy <- panel_owner$observations$y
  chunk_of <- (state$risk_a_idx - 1L) %/% as.integer(step) + 1L
  out <- vector("list", L)
  shared_ops <- 0
  endpoints_total <- 0

  for (ell in seq_len(L)) {
    ipw_owner_check_interrupt()
    lv <- link_views[[ell]]
    # 4. Raw endpoints of this link at its complete dates. Each posting slot is
    #    checked to hold its risk A-unit before any outcome is read; refusal
    #    precedence equals the chunked scalar stream (first failing chunk;
    #    within it, endpoints before differences).
    slot_s <- post$date_start[lv$pair$s_calendar_index] - 1L + state$risk_a_idx
    slot_t <- post$date_start[lv$pair$t_calendar_index] - 1L + state$risk_a_idx
    if (!isTRUE(all(post$post_a_index[slot_s] == state$risk_a_idx)) ||
        !isTRUE(all(post$post_a_index[slot_t] == state$risk_a_idx))) {
      ipw_point_view_batch_abort("Batch point fitting found a posting slot without its A-unit at a link date; no endpoint was read.",
        diagnostics = list(link = ell))
    }
    y_s <- oy[post$post_obs_index[slot_s]]
    y_t <- oy[post$post_obs_index[slot_t]]
    rm(slot_s, slot_t)
    risk_dy <- y_t - y_s
    if (anyNA(risk_dy) || any(!is.finite(risk_dy))) {
      bad_y <- !is.finite(y_s) | !is.finite(y_t)
      first_dy <- min(chunk_of[!is.finite(risk_dy)])
      endpoint_first <- any(bad_y) && min(chunk_of[bad_y]) == first_dy
      rm(y_s, y_t, risk_dy, bad_y)
      release_state()
      if (endpoint_first) {
        cdid_ipw_owner_abort(
          "Observed risk-set outcome endpoints contain NA, NaN, or Inf.",
          c("cdid_ipw_owner_endpoint_error", "cdid_ipw_enclosure_endpoint_error")
        )
      }
      cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
    }
    rm(y_s, y_t)
    endpoints_total <- endpoints_total + 2 * n_obs

    # 5. The shared reconstruction once (after the first link's reads, where the
    #    scalar fitter would first perform it), then this link's unchanged
    #    means and partials.
    tryCatch(
      {
        if (ell == 1L) {
          pr_chk <- ipw_pair_rates(state$risk_gi, state$risk_w, sup)
          if (!ipw_point_view_identical(as.integer(pr_chk$pair_count), as.integer(unname(fitted_view$counts$by_cohort_pair))) ||
              !ipw_point_view_identical(as.double(pr_chk$pair_mass), as.double(unname(fitted_view$masses$cohort_pair_weight))) ||
              !ipw_point_view_identical(as.double(pr_chk$rates), as.double(unname(fitted_view$masses$rates)))) {
            diag_env$error_diagnostics <- list(
              stage = "reconstructed_rates_check",
              reconstructed_rates = pr_chk$rates,
              fitted_rates = unname(fitted_view$masses$rates)
            )
            cdid_ipw_point_view_abort(
              "Reconstructed observer-buffer pair masses or rates disagree with the sealed fitted nuisance view.",
              c("cdid_ipw_point_view_integrity_error", "cdid_ipw_shared_mask_error")
            )
          }
          rm(pr_chk)
          rates <- as.double(unname(fitted_view$masses$rates))
          beta <- as.double(fitted_view$coefficients)
          names(beta) <- colnames(frame_owner$Z_a)
          state$r_use <- ipw_clip_floor(rates)[state$risk_gi]
          state$fit_weights <- state$risk_w / state$r_use
          if (any(!is.finite(state$fit_weights) | state$fit_weights <= 0)) {
            diag_env$error_diagnostics <- list(
              stage = "inverse_rate_weights",
              rates = rates,
              fit_weights = state$fit_weights
            )
            cdid_abort(
              "Nonfinite or nonpositive inverse-rate propensity weights `w_i / r_G`.",
              c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
            )
          }
          state$lv <- ipw_logit_values(state$Z_use, beta)
          if (any(!is.finite(state$lv$eta)) || any(!is.finite(state$lv$odds)) ||
              any(state$lv$odds <= 0) || any(!is.finite(state$lv$prob))) {
            diag_env$error_diagnostics <- list(
              stage = "logit_values_check",
              failed_coefficients = beta
            )
            cdid_abort(
              "Weighted propensity logit produced nonfinite linear predictors, odds, or probabilities.",
              c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
            )
          }
          state$ow <- ipw_obs_weights(state$fit_weights, state$lv$odds, state$risk_gi, sup$n)
          if (!ipw_point_view_identical(as.integer(sum(state$ow$is_t)), as.integer(fitted_view$counts$n_treated_obs)) ||
              !ipw_point_view_identical(as.integer(sum(state$ow$is_c)), as.integer(fitted_view$counts$n_control_obs)) ||
              !ipw_point_view_identical(as.double(state$ow$sum_a), as.double(fitted_view$means$sum_a)) ||
              !ipw_point_view_identical(as.double(state$ow$mean_a), as.double(fitted_view$means$mean_a)) ||
              !ipw_point_view_identical(as.double(state$ow$sum_b), as.double(fitted_view$means$sum_b)) ||
              !ipw_point_view_identical(as.double(state$ow$mean_b), as.double(fitted_view$means$mean_b))) {
            diag_env$error_diagnostics <- list(
              stage = "reconstructed_normalizer_check",
              reconstructed_sum_a = state$ow$sum_a,
              fitted_sum_a = fitted_view$means$sum_a,
              reconstructed_sum_b = state$ow$sum_b,
              fitted_sum_b = fitted_view$means$sum_b
            )
            cdid_ipw_point_view_abort(
              "Reconstructed observer normalizers disagree with the sealed fitted nuisance view.",
              c("cdid_ipw_point_view_integrity_error", "cdid_ipw_shared_mask_error")
            )
          }
          state$st <- list(
            n = sup$n,
            m_cohorts = m_i,
            p = p_i,
            risk_a_idx = state$risk_a_idx,
            risk_gi = state$risk_gi,
            risk_w = state$risk_w,
            pair_count = as.integer(unname(fitted_view$counts$by_cohort_pair)),
            pair_mass = as.double(unname(fitted_view$masses$cohort_pair_weight)),
            rates = rates,
            Z_use = state$Z_use,
            D_use = as.double(state$risk_gi == 1L),
            fit_weights = state$fit_weights,
            beta = beta,
            odds_use = state$lv$odds,
            prob_use = state$lv$prob,
            is_t = state$ow$is_t,
            is_c = state$ow$is_c,
            a_obs = state$ow$a_obs,
            b_obs = state$ow$b_obs,
            sum_a = fitted_view$means$sum_a,
            mean_a = fitted_view$means$mean_a,
            sum_b = fitted_view$means$sum_b,
            mean_b = fitted_view$means$mean_b
          )
          state$hs <- list(
            H = fitted_view$factors$H,
            H_chol = fitted_view$factors$H_chol,
            A_beta_r = fitted_view$partials$A_beta_r
          )
          shared_ops <- 16 * n_obs * p + 16 * n_obs * m
        }
        means <- tryCatch(
          ipw_link_means(state$st, risk_dy),
          error = function(e) {
            diag_env$error_diagnostics <- list(
              stage = "link_means_stage",
              sum_a = fitted_view$means$sum_a,
              sum_b = fitted_view$means$sum_b
            )
            stop(e)
          }
        )
        part <- tryCatch(
          ipw_link_partials(state$st, state$hs, risk_dy, means, sup$q, sup$cohort_labels),
          error = function(e) {
            diag_env$error_diagnostics <- list(
              stage = "link_partials_stage",
              treated_mean = means$Tmean,
              control_mean = means$Cmean,
              estimate = means$mu
            )
            stop(e)
          }
        )
      },
      error = fail
    )

    # 6. Compact extraction; this link's buffers are released before sealing.
    Tmean_val <- as.double(means$Tmean)
    Cmean_val <- as.double(means$Cmean)
    mu_val <- as.double(means$mu)
    G_beta_val <- part$G_beta
    G_r_val <- part$G_r
    loading_beta_val <- part$loading_beta
    loading_rate_val <- part$loading_rate
    rate_coef_s0_val <- part$rate_coef_s0
    rate_coef_s1_val <- part$rate_coef_s1
    rm(risk_dy, means, part)

    if (!is.null(lv$loading_beta_candidate) &&
        !ipw_point_view_identical(as.double(lv$loading_beta_candidate), as.double(loading_beta_val))) {
      cdid_ipw_point_view_abort(
        "Candidate `loading_beta_candidate` stored on `link_view` disagrees with the fitted point link's `loading_beta`.",
        c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error")
      )
    }

    link_pair_mass <- as.double(lv$support$by_cohort_pair_mass)
    names(link_pair_mass) <- sup$cohort_labels
    link_rates <- as.double(lv$support$rates)
    names(link_rates) <- sup$cohort_labels
    shared_across_dates <- !ipw_point_view_identical(lv$pair$s, fitted_view$defining_pair$s) ||
      !ipw_point_view_identical(lv$pair$t, fitted_view$defining_pair$t)
    actual_ops <- 6 * n_obs + 16 * n_obs * p + 16 * n_obs * m + 4 * (p * p + p * m) + (if (ell == 1L) shared_ops else 0)
    stopifnot(actual_ops <= link_work + shared_work)
    link_preflight <- list(
      route = "batch_shared_observer_state",
      n_a_units = na,
      n_observed_risk_units = n_obs,
      p_covariates = p,
      m_cohorts = m,
      n_planned_links = L,
      effective_chunk_rows = logit_cert_count_output(step),
      character_metadata = char_meta[[ell]],
      retained_matrix_entries_per_link = ret_matrix[[ell]],
      retained_vector_entries_per_link = ret_vector[[ell]],
      retained_entries_per_link = ret_one[[ell]],
      shared_state_entries = shared_state,
      link_workspace_entries = link_workspace,
      failure_output_entries = failure_output,
      batch_total_retained_entries = total_retained,
      batch_peak_concurrent_entries = peak,
      batch_total_work_ops_bound = total_work,
      work_ops_per_link_bound = link_work,
      expected_endpoints_read_per_link = logit_cert_count_output(2 * n_obs),
      propensity_refit_performed = FALSE,
      peak_rss_promised = FALSE
    )

    point_link <- structure(
      list(
        link_id = lv$link_id,
        functional_key = lv$functional_key,
        nuisance_id = fitted_view$nuisance_id,
        fitted_nuisance_ref = fitted_view$ref,
        nuisance_view_ref = fitted_view$nuisance_view_ref,
        nuisance_ref = fitted_view$nuisance_view_ref,
        link_view_ref = lv$ref,
        link_ref = lv$ref,
        frame_ref = fitted_view$frame_ref,
        panel_ref = fitted_view$panel_ref,
        estimate = mu_val,
        treated_mean = Tmean_val,
        control_mean = Cmean_val,
        pair = list(
          s = lv$pair$s,
          t = lv$pair$t,
          s_calendar_index = lv$pair$s_calendar_index,
          t_calendar_index = lv$pair$t_calendar_index,
          link_id = lv$link_id,
          functional_key = lv$functional_key
        ),
        defining_pair = fitted_view$defining_pair,
        shared_across_dates = shared_across_dates,
        nuisance_mask_verified = TRUE,
        cohorts = fitted_view$cohorts,
        dimensions = fitted_view$dimensions,
        counts = list(
          n_a_units = na_i,
          n_risk_a_units = fitted_view$counts$n_risk_a_units,
          n_outside_risk_a_units = fitted_view$counts$n_outside_risk_a_units,
          n_observed_pair_units = as.integer(lv$support$n_observed_pair_all),
          n_observed_risk_units = n_obs_i,
          n_nonobserved_risk_units = fitted_view$counts$n_nonobserved_risk_units,
          n_treated_obs = fitted_view$counts$n_treated_obs,
          n_control_obs = fitted_view$counts$n_control_obs,
          by_cohort_a = fitted_view$counts$by_cohort_a,
          by_cohort_pair = fitted_view$counts$by_cohort_pair
        ),
        masses = fitted_view$masses,
        means = list(
          treated_mean = Tmean_val,
          control_mean = Cmean_val,
          estimate = mu_val,
          sum_a = fitted_view$means$sum_a,
          mean_a = fitted_view$means$mean_a,
          sum_b = fitted_view$means$sum_b,
          mean_b = fitted_view$means$mean_b
        ),
        coefficients = fitted_view$coefficients,
        beta = fitted_view$beta,
        factors = fitted_view$factors,
        partials = list(
          G_beta = G_beta_val,
          G_r = G_r_val,
          A_beta_r = fitted_view$partials$A_beta_r
        ),
        loadings = list(
          loading_beta = loading_beta_val,
          loading_rate = loading_rate_val,
          rate_coef_s0 = rate_coef_s0_val,
          rate_coef_s1 = rate_coef_s1_val
        ),
        diagnostics = fitted_view$diagnostics,
        point_summary_disclosure = list(
          defining_rate_source = fitted_view$point_summary_disclosure$defining_rate_source,
          observer_cohort_a_mass = fitted_view$point_summary_disclosure$observer_cohort_a_mass,
          observer_cohort_pair_mass = fitted_view$point_summary_disclosure$observer_cohort_pair_mass,
          observer_rates = fitted_view$point_summary_disclosure$observer_rates,
          view_chunked_cohort_a_mass = fitted_view$point_summary_disclosure$view_chunked_cohort_a_mass,
          view_chunked_cohort_pair_mass = fitted_view$point_summary_disclosure$view_chunked_cohort_pair_mass,
          view_chunked_rates = fitted_view$point_summary_disclosure$view_chunked_rates,
          link_chunked_cohort_pair_mass = link_pair_mass,
          link_chunked_rates = link_rates,
          max_abs_pair_mass_chunk_diff = max(abs(fitted_view$masses$cohort_pair_weight - link_pair_mass)),
          max_abs_rate_chunk_diff = max(abs(fitted_view$masses$rates - link_rates)),
          propensity_refit_performed = FALSE,
          candidate_point_summary_only = TRUE,
          certifies_root_or_inference = FALSE,
          glm_convergence_is_certificate = FALSE,
          seal_is_certificate = FALSE,
          disclosure = paste0(
            "Compact point link reuses exact fitted coefficients, observer-buffer rates, normalizers, ",
            "Hessian Cholesky factor, and score-rate Jacobian without refitting the propensity model. ",
            "Chunked link-view support summaries are disclosed alongside observer-buffer summaries. ",
            "Sealing records a finite binary64 point candidate only and never certifies raw normalization, ",
            "root enclosure, covariance transfer, or statistical inference."
          )
        ),
        glm_control = fitted_view$glm_control,
        control = ctrl,
        preflight = link_preflight,
        accounting = list(
          retained_matrix_entries = ret_matrix[[ell]],
          retained_vector_entries = ret_vector[[ell]],
          total_retained_entries = ret_one[[ell]],
          borrowed_fitted_nuisance_ref = TRUE,
          borrowed_nuisance_view_ref = TRUE,
          borrowed_link_view_ref = TRUE,
          borrowed_frame_ref = TRUE,
          borrowed_panel_ref = TRUE,
          shared_observer_state = TRUE,
          observer_buffer_released = FALSE,
          observer_state_release = "after_last_link_of_batch",
          reader_chunks_released = TRUE,
          full_n_masks_retained = FALSE,
          observer_arrays_retained = FALSE,
          copied_Z_retained = FALSE,
          per_link_endpoint_copies_retained = FALSE,
          propensity_refit_performed = FALSE,
          endpoints_read = as.integer(2 * n_obs),
          peak_rss_promised = FALSE
        ),
        work = list(
          chunks_queried = 0L,
          mask_equality_decided_by_posting_fact = TRUE,
          max_simultaneous_query_chunks = 0L,
          max_simultaneous_observer_buffers = 1L,
          observer_units_copied = logit_cert_count_output(if (ell == 1L) n_obs else 0),
          endpoints_read = as.integer(2 * n_obs),
          propensity_solver_calls = 0L,
          actual_work_ops = actual_ops
        )
      ),
      class = "cdid_ipw_point_link_view"
    )
    ipw_point_view_verify_compact_storage(point_link, p_max = p, m_max = m)
    actual_payload <- ipw_point_view_payload_budget(point_link, ctrl)
    if (actual_payload$entries > ret_one[[ell]]) {
      cdid_ipw_point_view_abort("Point-link output exceeded its proved retained envelope.",
        "cdid_ipw_point_view_bound_error")
    }
    out[[ell]] <- ipw_owner_seal(point_link, "point_link_view")
    if (!is.null(trace_env)) {
      trace_env$links_fitted <- ell
      trace_env$endpoints_read <- endpoints_total
    }
  }

  release_state()
  names(out) <- vapply(out, `[[`, character(1L), "link_id")
  out
}

# Sealed invocation-local plan for the eligible retained links of one fitted
# nuisance. It performs the single normalization and the unchanged general
# producer preflight (whose bounds stay a conservative reservation), takes
# every link's complete-pair fact from the records that preflight creates, and
# stores compact descriptors: calendar indices, treated/control means, the
# normalizers, the p beta loadings (p = total propensity columns including the
# intercept) and rate_coef_s1 (2 + 4 + p + m per link; the p = 1 layout is the
# original 2 + 5 + m). The loadings are each link's original H_chol triangular
# solves and the coefficients the single original fitted propensity; nothing is
# refitted, inverted or admitted by closeness.
# The batch charge (descriptors and identity metadata retained through the GMM
# projection, plus the chunk's shared arrays and kernel output) is composed with
# that preflight before any descriptor is allocated; the pipeline then composes
# it into its simultaneous envelope.
ipw_point_if_batch_plan <- function(point_links, control, fitted_view) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control) || missing(fitted_view)) {
    cdid_ipw_point_view_abort(
      "Required arguments (`point_links`, `control`, `fitted_view`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  pre_norm <- preflight_ipw_point_if_view(point_links, control, include_components = FALSE, return_norm = TRUE)
  pre <- pre_norm$preflight
  ctrl <- pre_norm$ctrl
  norm <- pre_norm$norm
  pre_norm <- NULL
  L <- as.integer(norm$n_links)
  fo <- norm$frame_owner
  po <- norm$panel_owner
  ipw_point_view_seal_gate(fitted_view, "fitted_nuisance_view")
  nv <- get("nuisance_view", envir = fitted_view$nuisance_view_ref, inherits = FALSE)
  m <- as.integer(nv$cohorts$m_cohorts)
  # The shared fixed-unit design is the owner's n_A x p double matrix (p >= 1
  # including the intercept), structurally, with the general preflight's p; no
  # design entry is scanned here and no value closeness admits the route.
  p_ok <- is.numeric(fo$n_covariates) && length(fo$n_covariates) == 1L && !is.na(fo$n_covariates) &&
    fo$n_covariates >= 1 && fo$n_covariates == floor(fo$n_covariates)
  p_i <- if (p_ok) as.integer(fo$n_covariates) else NA_integer_
  ok <- identical(norm$n_unique_fitted_nuisances, 1L) && p_ok && identical(as.integer(pre$p_covariates), p_i) &&
    is.matrix(fo$Z_a) && is.double(fo$Z_a) && identical(dim(fo$Z_a), c(as.integer(fo$n_a_units), p_i)) &&
    is.double(fo$w_a) && isTRUE(all(fo$w_a == fo$w_a[1L])) && identical(fo$n_a_units, po$n_a_units) &&
    identical(fitted_view$frame_ref, fo$ref) && length(norm$validated) == L
  if (!ok) {
    ipw_point_view_batch_abort("Batch IF plan requires one fitted nuisance, the owner's shared fixed-unit n_A x p design, constant unit weights and complete link pairs.")
  }
  max_lim <- min(ctrl$max_matrix_entries, ctrl$max_output_entries)
  B <- as.double(pre$effective_chunk_rows)
  p <- as.double(p_i)
  plan_retained <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(2 + 4 + p + m + 16, L, label = "batch IF plan descriptors and identity metadata"),
    ipw_point_view_char_bytes(list(norm$link_ids), ctrl)$charged_slots, 512),
    min(ctrl$max_retained_entries, max_lim), "batch IF plan retained descriptors")
  chunk_shared <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(24, B, max_lim, "chunk index/ID/weight/cohort/observer/rate/odds/score arrays and temporaries"),
    ipw_owner_safe_product(B, p, max_lim, "chunk observer design rows"),
    if (p_i > 1L) ipw_owner_safe_sum(c(B, p), max_lim, "exact matrix-product result and private loading vector") else 0,
    pre$id_byte_metadata$chunk_id_slots, 128), max_lim, "batch IF chunk shared observer arrays")
  chunk_output <- ipw_owner_safe_sum(c(pre$b_by_L_entries, 64), max_lim, "batch IF kernel output, status and batch list")
  producer <- ipw_owner_safe_sum(c(chunk_shared, chunk_output), max_lim, "batch IF producer chunk workspace")
  additional_live <- ipw_owner_safe_sum(c(plan_retained, max(0, producer - pre$producer_workspace_entries)))
  ipw_owner_safe_sum(c(pre$peak_concurrent_entries, additional_live), ctrl$max_matrix_entries,
    "batch IF plan plus point-IF producer concurrent entries")
  # Traversal work: per chunk the shared observer arrays plus the B x p design
  # rows and their p-term propensity products; per unit and link the kernel's
  # posting checks, endpoint reads, raw difference, p-term loading product and
  # four-term IF. General p also reserves each exact R primitive dispatch.
  # The p = 1 terms are the original 64 * B and 40 * L.
  traversal_work <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(as.double(pre$n_chunks), 64 * B + 4 * (p - 1) * B + 4 * L + (if (p_i > 1L) 128 * L else 0) + 512, ctrl$max_work_ops,
      "chunk shared arrays, design rows and propensity products, kernel dispatch and binding checks"),
    ipw_owner_safe_product(as.double(fo$n_a_units), (40 + 4 * (p - 1)) * L, ctrl$max_work_ops,
      "kernel posting checks, endpoint reads, raw differences, p-term loading products and four-term IF"),
    256), ctrl$max_work_ops, "batch IF plan traversal work")
  link_views <- vector("list", L)
  link_dates <- matrix(0L, nrow = 2L, ncol = L)
  link_par <- matrix(0.0, nrow = 4L + p_i + m, ncol = L)
  z_names <- colnames(fo$Z_a)
  for (ell in seq_len(L)) {
    if (!ok) break
    lk <- norm$links[[ell]]
    lv <- get("link_view", envir = lk$link_view_ref, inherits = FALSE)
    # Exact descriptor shapes and column order: the link's coefficients are the
    # one fitted propensity (p entries); its beta loading has p entries in the
    # design's positional column order; names, when present, must agree. The
    # canonical validator also permits unnamed loadings. rate_coef_s1 has one
    # entry per ordered cohort.
    ok <- identical(lk$fitted_nuisance_ref, fitted_view$ref) && identical(lv$cohorts, nv$cohorts) &&
      isTRUE(norm$validated[[ell]]$complete_pair) && identical(norm$validated[[ell]]$view, lv) &&
      is.double(lk$coefficients) && length(lk$coefficients) == p_i &&
      identical(lk$coefficients, fitted_view$coefficients) &&
      is.double(lk$loadings$loading_beta) && length(lk$loadings$loading_beta) == p_i &&
      (is.null(names(lk$loadings$loading_beta)) || identical(names(lk$loadings$loading_beta), z_names)) &&
      is.double(lk$loadings$rate_coef_s1) && length(lk$loadings$rate_coef_s1) == m
    if (!ok) break
    link_views[[ell]] <- lv
    link_dates[, ell] <- c(as.integer(lv$pair$s_calendar_index), as.integer(lv$pair$t_calendar_index))
    link_par[, ell] <- c(lk$means$treated_mean, lk$means$control_mean, lk$means$mean_a, lk$means$mean_b,
      as.double(lk$loadings$loading_beta), as.double(lk$loadings$rate_coef_s1))
  }
  if (!ok || anyNA(link_dates) || !identical(dim(link_par), c(4L + p_i + m, L))) {
    ipw_point_view_batch_abort("Batch IF plan requires one fitted nuisance, the owner's shared fixed-unit n_A x p design, constant unit weights and complete link pairs.")
  }

  plan <- list(
    kind = "cdid_ipw_point_if_batch_plan",
    control = ctrl,
    control_in = control,
    frame_owner = fo,
    panel_owner = po,
    fitted_view = fitted_view,
    nuisance_view = nv,
    point_links = norm$links,
    link_views = link_views,
    link_ids = norm$link_ids,
    dimnames = list(NULL, norm$link_ids),
    cohorts = nv$cohorts,
    rates = norm$links[[1L]]$masses$rates,
    coefficients = norm$links[[1L]]$coefficients,
    link_dates = link_dates,
    link_par = link_par,
    n_a_units = as.integer(fo$n_a_units),
    step = as.integer(B),
    max_output_entries = max_lim,
    preflight = pre,
    batch = list(
      route = "batch_if_plan",
      kernel = "C_gmm_ipw_if_chunk",
      n_links = L,
      p_covariates = p_i,
      m_cohorts = m,
      effective_chunk_rows = pre$effective_chunk_rows,
      n_chunks = pre$n_chunks,
      normalizations_performed = 1L,
      plan_retained_entries = plan_retained,
      chunk_shared_entries = chunk_shared,
      chunk_output_entries = chunk_output,
      producer_workspace_entries = producer,
      general_producer_workspace_entries = pre$producer_workspace_entries,
      additional_live_entries = additional_live,
      traversal_work_ops_bound = traversal_work,
      full_n_by_L_retained = FALSE,
      per_link_dense_masks_retained = FALSE,
      peak_rss_promised = FALSE,
      lifetime = paste0("Descriptors and identity metadata persist from plan creation through the GMM projection; ",
        "per traversal one chunk's shared arrays and kernel output are live with the consumer, then released.")
    )
  )
  ipw_owner_seal(plan, "point_if_batch_plan")
}

# Plan binding check at every traversal boundary and after every callback,
# before the next endpoint read: the plan's own seal, the caller's controls,
# and the objects bound in the frame/panel/fitted-nuisance/nuisance/link/point
# environments must be exactly those the plan was built from.
ipw_point_if_plan_verify <- function(plan, control) {
  if (!is.list(plan) || !identical(plan$kind, "cdid_ipw_point_if_batch_plan") ||
      !ipw_owner_seal_ok(plan, "point_if_batch_plan")) {
    ipw_point_view_batch_abort("Batch IF plan is malformed, foreign, or differs from its sealed snapshot.")
  }
  if (!identical(control, plan$control_in) && !identical(control, plan$control)) {
    cdid_ipw_point_view_abort("Batch IF plan controls differ from the traversal controls; reuse is refused.",
      c("cdid_ipw_control_mismatch_error", "cdid_ipw_point_view_batch_error"))
  }
  fo <- plan$frame_owner
  fv <- plan$fitted_view
  bound <- ipw_owner_seal_ok(fo, "frame_owner") &&
    identical(fo$panel_ref, plan$panel_owner$ref) &&
    ipw_owner_seal_ok(plan$panel_owner, "panel_owner") &&
    ipw_owner_seal_ok(fv, "fitted_nuisance_view") &&
    identical(fv$frame_ref, fo$ref) && identical(fv$panel_ref, fo$panel_ref) &&
    identical(fv$nuisance_view_ref, plan$nuisance_view$ref) &&
    identical(plan$nuisance_view$frame_ref, fo$ref) && identical(plan$nuisance_view$panel_ref, fo$panel_ref) &&
    ipw_owner_seal_ok(plan$nuisance_view, "nuisance_view")
  for (ell in seq_along(plan$point_links)) {
    if (!bound) break
    pl <- plan$point_links[[ell]]
    lv <- plan$link_views[[ell]]
    bound <- ipw_owner_seal_ok(pl, "point_link_view") &&
      ipw_owner_seal_ok(lv, "link_view") &&
      identical(pl$frame_ref, fo$ref) && identical(pl$panel_ref, fo$panel_ref) &&
      identical(lv$frame_ref, fo$ref) && identical(lv$panel_ref, fo$panel_ref) &&
      identical(pl$link_view_ref, lv$ref) && identical(pl$link_ref, lv$ref) &&
      identical(pl$fitted_nuisance_ref, fv$ref) &&
      identical(pl$nuisance_view_ref, plan$nuisance_view$ref) &&
      identical(lv$nuisance_ref, plan$nuisance_view$ref)
  }
  if (!bound) {
    ipw_point_view_batch_abort("Batch IF plan owner, nuisance, link or point-link bindings changed; reuse is refused.")
  }
  invisible(TRUE)
}

ipw_point_if_plan_refuse <- function(status) {
  code <- status[1L]
  if (identical(code, 1L)) {
    cdid_ipw_owner_abort(
      "Observed risk-set outcome endpoints contain NA, NaN, or Inf.",
      c("cdid_ipw_owner_endpoint_error", "cdid_ipw_enclosure_endpoint_error")
    )
  }
  if (identical(code, 2L)) cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
  if (identical(code, 3L)) cdid_abort("Nonfinite IF contribution encountered during replay.", c("cdid_ipw_numeric_error"))
  ipw_point_view_batch_abort("Batch IF plan posting slot does not hold its A-unit at a link date; no endpoint was read.",
    diagnostics = list(status = status))
}

# Optional component output adds four B x L matrices to the plan's ordinary
# producer reservation. Compose the same output/copy allowances as the general
# preflight with the retained plan before any chunk arrays, reads or callbacks.
# The production IF-only path returns its original preflight unchanged.
ipw_point_if_plan_preflight <- function(plan, include_components) {
  pre <- plan$preflight
  if (!isTRUE(include_components)) return(pre)
  ctrl <- plan$control
  extra <- ipw_owner_safe_product(4, pre$b_by_L_entries,
    min(ctrl$max_matrix_entries, ctrl$max_output_entries), "batch IF optional component matrices")
  output <- ipw_owner_safe_sum(c(pre$chunk_output_entries, extra),
    min(ctrl$max_matrix_entries, ctrl$max_output_entries), "batch IF total output with components")
  ipw_owner_safe_sum(c(pre$peak_concurrent_entries, 2 * extra, plan$batch$additional_live_entries),
    ctrl$max_matrix_entries, "batch IF components plus retained plan concurrent entries")
  pre$include_components <- TRUE
  pre$n_output_matrices_per_chunk <- 5L
  pre$chunk_output_matrix_entries <- pre$chunk_output_matrix_entries + extra
  pre$chunk_output_entries <- output
  pre$peak_concurrent_entries <- pre$peak_concurrent_entries + 2 * extra
  pre$peak_concurrent_matrix_entries <- pre$peak_concurrent_matrix_entries + extra
  pre$peak_concurrent_bytes_8double <- 8 * pre$peak_concurrent_entries
  pre$producer_workspace_entries <- pre$producer_workspace_entries + 2 * extra
  pre
}

# One bounded chunk of every retained link's complete IF. The shared observer
# arrays are formed here exactly as the general replay's one-nuisance chunk
# cache forms them (same R operations on the same owner values); the kernel
# writes the B x L IF (and, on request, its four components).
ipw_point_if_plan_chunk <- function(plan, a_start, a_end, kernel, include_components = FALSE) {
  if (isTRUE(include_components)) ipw_point_if_plan_preflight(plan, TRUE)
  a_s <- as.integer(a_start)
  a_e <- as.integer(a_end)
  a_seq <- seq.int(a_s, a_e)
  b_len <- length(a_seq)
  fo <- plan$frame_owner
  po <- plan$panel_owner
  u_id <- po$a_id[a_seq]
  u_w <- fo$w_a[a_seq]
  gi <- ipw_owner_view_cohort_index(plan$cohorts, po$a_cohort_index[a_seq])
  pos_k <- which(!is.na(gi))
  gj_k <- gi[pos_k]
  w_k <- u_w[pos_k]
  r_k <- ipw_clip_floor(plan$rates[gj_k])
  D_k <- as.double(gj_k == 1L)
  Z_k <- fo$Z_a[a_seq[pos_k], , drop = FALSE]
  lv_k <- ipw_logit_values(Z_k, plan$coefficients)
  odds_k <- lv_k$odds
  prob_k <- lv_k$prob
  score_base_k <- (w_k / r_k) * (D_k - prob_k)
  rm(gi, D_k, lv_k, prob_k)
  post <- po$postings
  res <- .Call(kernel, po$observations$y, post$date_start, post$post_a_index, post$post_obs_index,
    a_s, b_len, pos_k, gj_k, w_k, r_k, odds_k, score_base_k, Z_k, plan$link_dates, plan$link_par,
    plan$dimnames, isTRUE(include_components), plan$max_output_entries)
  rm(pos_k, gj_k, w_k, r_k, odds_k, score_base_k, Z_k)
  if (!identical(res$status[1L], 0L)) ipw_point_if_plan_refuse(res$status)
  out <- list(phi = res$phi, id = u_id, link_ids = plan$link_ids, a_index = a_seq)
  if (isTRUE(include_components)) {
    out$components <- list(rate = res$rate, beta_score = res$beta_score, treated = res$treated, control = res$control)
  }
  out
}

replay_ipw_point_if_plan <- function(plan, control, callback, include_components = FALSE) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control) || missing(callback)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`plan`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.function(callback)) {
    cdid_ipw_point_view_abort("`callback` must be a function.", c("cdid_ipw_point_view_input_error", "cdid_ipw_input_error"))
  }
  if (!gmm_scalar_logical(include_components)) {
    cdid_ipw_point_view_abort("`include_components` must be a scalar logical.", "cdid_ipw_point_view_input_error")
  }
  ipw_point_if_plan_verify(plan, control)
  pre <- ipw_point_if_plan_preflight(plan, include_components)
  na <- plan$n_a_units
  step <- plan$step
  emitted <- 0L
  batches <- 0L
  max_batch_units <- 0L
  a_start <- 1
  while (a_start <= na) {
    ipw_owner_check_interrupt()
    a_end <- min(as.double(na), a_start + step - 1)
    batches <- batches + 1L
    batch <- ipw_point_if_plan_chunk(plan, a_start, a_end, ipw_point_if_native_kernel(), include_components)
    callback(batch)
    batch <- NULL
    ipw_point_if_plan_verify(plan, control)
    b_len <- as.integer(a_end - a_start + 1)
    emitted <- emitted + b_len
    max_batch_units <- max(max_batch_units, b_len)
    a_start <- a_end + 1
  }
  if (emitted != na) {
    cdid_ipw_point_view_abort("Point IF replay unit count mismatch.", c("cdid_ipw_count_error"))
  }
  invisible(list(
    n_a_units = emitted,
    n_links = length(plan$link_ids),
    link_ids = plan$link_ids,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * length(plan$link_ids),
    kernel_calls = batches,
    normalizations_performed = 0L,
    preflight = pre
  ))
}

# ---------------------------------------------------------------------------
# CP08/D92 retained point IF (invocation-local carrier). For ONE sealed batch,
# general, or static custom IF plan without numerical transfer, the FIRST
# covariance replay's exact raw N_A x L phi (every A unit in plan order, the
# plan's links in their order, the producer's complete A-unit IF under the
# plan's one normalization) is written once into a carrier and sliced, in the
# plan's own row partition, for the covariance Gram pass and both GMM
# projection passes. The producer is the unchanged `replay_ipw_point_if_plan`,
# `replay_ipw_point_if_general`, or `replay_custom_ipw_if_plan`; the carrier
# only observes its callback. Nothing centered or projected is stored; nothing
# is attached to the sealed plan, the owners, the fitted view, the returned MP
# or any global. Full-frame transport (`n_total >= n_a_units`) and cluster
# aggregation happen downstream of the A-unit GMM projection as before.
# Admission is decided from the plan's sealed facts and exact arithmetic BEFORE
# any allocation; a typed capacity decline leaves every original charge,
# partition and replay untouched, while malformed inputs, integrity violations
# and exact-count overflows propagate.
# ---------------------------------------------------------------------------

# Charged entries: the N_A x L carrier plus one possible duplicate, 2 N_A ID/
# index slots, 4 B x L slice scratch and 256 headers. The Gram-side buffers are
# already charged by the covariance plan's own check.
ipw_point_if_retained_entries <- function(N, L, B, max_lim = NULL) {
  nl <- ipw_owner_check_matrix_dims(N, L, max_matrix_entries = max_lim, label = "CP08 retained point IF")
  two_nl <- ipw_owner_safe_product(2, nl, max_lim, "CP08 retained point IF plus copy reserve")
  two_n <- ipw_owner_safe_product(2, N, max_lim, "CP08 retained point IF unit ID/index slots")
  bl <- ipw_owner_check_matrix_dims(B, L, max_matrix_entries = max_lim, label = "CP08 retained point IF slice chunk")
  four_bl <- ipw_owner_safe_product(4, bl, max_lim, "CP08 retained point IF slice scratch")
  ipw_owner_safe_sum(
    c(two_nl, two_n, four_bl, 256),
    max_lim,
    "CP08 retained point IF reserved entries"
  )
}
# Charged work: the fill copy, boundary binding checks and three slice passes.
ipw_point_if_retained_work <- function(N, L, max_wo = NULL) {
  nl <- ipw_owner_safe_product(N, L, max_wo, "CP08 retained point IF elements")
  ipw_owner_safe_product(128, nl, max_wo, "CP08 retained point IF fill and slice work")
}

# Additional counted storage/work introduced by this carrier, independent of
# the original producer's work. All arithmetic here is uncapped: an exact-count
# failure cannot be hidden by an earlier, otherwise valid capacity decline.
ipw_point_if_retained_inventory <- function(plan) {
  N <- as.double(plan$n_a_units)
  L <- as.double(length(plan$link_ids))
  B <- as.double(plan$step)
  C <- as.double(plan$preflight$n_chunks)
  P <- ipw_owner_safe_product
  S <- function(x, label) ipw_owner_safe_sum(x, label = label)
  controls <- ipw_point_view_payload_budget(list(plan$control_in, plan$control), list())
  labels <- ipw_point_view_char_bytes(list(plan$link_ids), list())
  ids <- ipw_point_view_id_bytes_bound(plan$panel_owner, B, C, list())
  one_plan_check <- if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
    plan$general$one_verification_work_ops_bound
  } else if (identical(plan$kind, "cdid_custom_ipw_if_plan")) {
    plan$custom$one_verification_work_ops_bound
  } else S(c(P(64, L), 512), "retained batch plan seal checks")
  one_check <- S(c(one_plan_check, P(2, controls$work), P(128, L), 2048),
    "retained boundary plan/control/carrier checks")
  # Covers reservation, covariance entry, fill, three slice passes, and checks
  # immediately before and after every callback, including the last callback.
  checks <- P(16, S(c(C, 1), "retained boundary count"))
  boundary_work <- P(checks, one_check, label = "retained boundary validation work")
  metadata <- S(c(P(2, controls$entries), labels$charged_slots,
    ids$chunk_id_slots, P(128, L), 4096), "retained labels, snapshots and headers")
  # Reserve the complete centered Gram workspace also in the caller's live
  # envelope. GMM target/reporting buffers are supplied by that caller in live.
  covariance <- S(c(P(4, ipw_owner_check_matrix_dims(L, L)), P(8, L),
    P(6, ipw_owner_check_matrix_dims(B, L)), 128), "retained centered Gram workspace")
  extra <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) plan$batch else
    if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) plan$general else plan$custom
  base_live <- S(c(plan$preflight$borrowed_owner_entries,
    plan$preflight$peak_concurrent_entries, extra$additional_live_entries),
    "retained borrowed owners, live producer and plan")
  base_work <- S(c(P(4, plan$preflight$total_work_ops_bound),
    if (is.null(extra$additional_work_ops_bound)) 0 else extra$additional_work_ops_bound,
    P(N, S(c(P(2, P(L, L)), P(4, L)), "retained covariance arithmetic per row")), 128),
    "retained original four-pass producer and covariance work")
  list(
    entries = S(c(ipw_point_if_retained_entries(N, L, B), metadata, covariance),
      "complete retained IF simultaneous inventory"),
    work = S(c(ipw_point_if_retained_work(N, L),
      P(C, P(N, L), label = "possible full raw matrix COW copy at every fill chunk"),
      boundary_work, controls$work,
      P(checks, labels$work_ops_bound), P(4, ids$stream_id_work_bound)),
      "complete retained IF additional work"),
    metadata_entries = metadata,
    covariance_entries = covariance,
    boundary_work = boundary_work,
    base_live = base_live,
    base_work = base_work
  )
}

ipw_point_if_retained_plan_verify <- function(plan, control) {
  if (!is.list(plan) || !is.character(plan$kind) || length(plan$kind) != 1L) {
    ipw_point_view_batch_abort("Retained point IF plan is malformed or missing its kind tag.")
  }
  if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) {
    return(ipw_point_if_plan_verify(plan, control))
  }
  if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
    return(ipw_point_if_general_verify(plan, control))
  }
  if (identical(plan$kind, "cdid_custom_ipw_if_plan")) {
    return(custom_ipw_if_plan_verify(plan, control))
  }
  ipw_point_view_batch_abort("Retained point IF carrier supports only sealed batch, general, or custom point-IF plans.")
}

# Reserve a carrier or return NULL. Plan integrity and exact-count bounds run
# first and propagate on violation; numerical transfer declines; then the exact
# arithmetic admission of the carrier inside the caller's complete simultaneous
# live envelope and work cap and the plan's own covariance-stream charges is
# probed, declining only on a typed capacity bound refusal. The plan's row
# partition is sealed in the plan (`step`), so the reservation cannot alter it.
# Nothing is allocated here.
ipw_point_if_retained_reserve <- function(plan, control, live, work, envelope_control,
                                          n_total = NULL, clustered = FALSE, numerical_transfer = FALSE) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control) || missing(live) || missing(work) || missing(envelope_control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`plan`, `control`, `live`, `work`, `envelope_control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  ipw_point_if_retained_plan_verify(plan, control)
  if (!gmm_scalar_logical(clustered) || !gmm_scalar_logical(numerical_transfer)) {
    cdid_ipw_point_view_abort(
      "`clustered` and `numerical_transfer` must be scalar logicals.",
      "cdid_ipw_point_view_input_error"
    )
  }
  scalar_nonneg <- function(x) is.numeric(x) && !is.object(x) && length(x) == 1L && is.finite(x) && x >= 0
  scalar_pos <- function(x) scalar_nonneg(x) && x > 0
  count_pos <- function(x) scalar_nonneg(x) && x >= 1 && x == floor(x) && x <= .Machine$integer.max
  if (!scalar_nonneg(live) || !scalar_nonneg(work)) {
    cdid_ipw_point_view_abort(
      "`live` and `work` must be finite non-negative scalar numbers.",
      "cdid_ipw_point_view_input_error"
    )
  }
  live <- ipw_owner_int_bound(live, "retained point IF external live entries", min_val = 0)
  work <- ipw_owner_int_bound(work, "retained point IF external work", min_val = 0)
  if (!is.list(envelope_control) ||
      !scalar_pos(envelope_control$max_retained_entries) ||
      !scalar_pos(envelope_control$max_matrix_entries) ||
      !scalar_pos(envelope_control$max_output_entries) ||
      !scalar_pos(envelope_control$max_work_ops)) {
    cdid_ipw_point_view_abort(
      "`envelope_control` must provide finite positive retained/matrix/output/work caps.",
      "cdid_ipw_point_view_input_error"
    )
  }
  for (nm in c("max_retained_entries", "max_matrix_entries", "max_output_entries", "max_work_ops")) {
    ipw_owner_int_bound(envelope_control[[nm]], paste0("envelope_control$", nm))
  }
  fo <- plan$frame_owner
  po <- plan$panel_owner
  pre <- plan$preflight
  ctrl <- plan$control
  N <- plan$n_a_units
  L <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) plan$batch$n_links else plan$n_links
  B <- plan$step
  if (is.null(n_total)) n_total <- N
  if (!count_pos(n_total) || as.double(n_total) < as.double(N) ||
      !identical(as.double(fo$n_a_units), as.double(N)) ||
      !identical(as.double(po$n_a_units), as.double(N))) {
    cdid_ipw_point_view_abort(
      "Retained point IF carrier frame/A-unit count `n_total` is invalid or smaller than `n_a_units`.",
      c("cdid_ipw_count_error", "cdid_ipw_point_view_integrity_error")
    )
  }
  if (!is.double(fo$w_a) || length(fo$w_a) != as.integer(N) || any(!is.finite(fo$w_a)) || any(fo$w_a <= 0)) {
    cdid_ipw_point_view_abort(
      "Retained point IF carrier frame owner weights `w_a` are malformed.",
      c("cdid_ipw_point_view_integrity_error")
    )
  }
  if (!is.list(pre) || !count_pos(N) || !count_pos(L) || !count_pos(B) ||
      B > N || B > ctrl$max_chunk_rows ||
      !identical(length(plan$link_ids), as.integer(L)) ||
      !identical(as.double(pre$n_a_units), as.double(N)) ||
      !identical(as.double(pre$n_links), as.double(L)) ||
      !identical(as.double(pre$effective_chunk_rows), as.double(B)) ||
      !count_pos(pre$n_chunks) ||
      !identical(as.double(pre$n_chunks), ceiling(as.double(N) / as.double(B)))) {
    cdid_ipw_point_view_abort(
      "Retained point IF plan dimensions or chunk partition are inconsistent.",
      c("cdid_ipw_count_error", "cdid_ipw_point_view_integrity_error")
    )
  }
  if (isTRUE(numerical_transfer)) return(NULL)

  N_d <- as.double(N)
  L_d <- as.double(L)
  B_d <- as.double(B)
  C_d <- as.double(pre$n_chunks)
  add_live <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) {
    plan$batch$additional_live_entries
  } else if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
    plan$general$additional_live_entries
  } else {
    plan$custom$additional_live_entries
  }
  one_verify_work <- if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
    plan$general$one_verification_work_ops_bound
  } else if (identical(plan$kind, "cdid_custom_ipw_if_plan")) {
    plan$custom$one_verification_work_ops_bound
  } else {
    0
  }
  add_plan_work <- if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
    plan$general$additional_work_ops_bound
  } else if (identical(plan$kind, "cdid_custom_ipw_if_plan")) {
    plan$custom$additional_work_ops_bound
  } else {
    0
  }

  # Exact-count arithmetic and matrix dimension limits must propagate even when
  # a separate capacity cap would also decline. Compute the inventory uncapped.
  inventory <- ipw_point_if_retained_inventory(plan)
  entries_uncapped <- inventory$entries
  cache_work_uncapped <- inventory$work
  # A direct caller may supply zero external live/work. The producer and its
  # borrowed owners still exist; the caller's envelope must dominate them.
  live <- max(live, inventory$base_live)
  work <- max(work, inventory$base_work)
  cov_uncapped <- ipw_owner_check_matrix_dims(L, L, label = "streamed L x L point IF covariance matrix")
  ipw_owner_safe_sum(c(live, entries_uncapped), label = "pipeline simultaneous live entries plus retained point IF")
  ipw_owner_safe_sum(c(work, cache_work_uncapped), label = "pipeline cumulative work plus retained point IF")
  ipw_owner_safe_sum(
    c(pre$peak_concurrent_entries, add_live, 4 * cov_uncapped, 8 * L_d, 6 * pre$b_by_L_entries, 128, entries_uncapped),
    label = "two-pass streamed point covariance plus retained point IF concurrent entries"
  )
  cov_verify_uncapped <- if (one_verify_work > 0) {
    ipw_owner_safe_product(1 + 2 * (C_d + 1), one_verify_work, label = "covariance entry and two replay record verifications")
  } else 0
  ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(2, pre$total_work_ops_bound, label = "two-pass point IF replay work"),
      cov_verify_uncapped,
      ipw_owner_safe_product(N_d, 2 * L_d * L_d + 4 * L_d, label = "centering and crossprod accumulation work"),
      128,
      cache_work_uncapped
    ),
    label = "two-pass streamed point covariance plus retained point IF work"
  )
  ipw_owner_safe_sum(
    c(pre$total_work_ops_bound, add_plan_work, cache_work_uncapped),
    label = "point IF plan work plus retained point IF work"
  )

  max_lim_ctrl <- min(ctrl$max_retained_entries, ctrl$max_matrix_entries)
  max_lim_env <- min(envelope_control$max_retained_entries, envelope_control$max_matrix_entries)
  cap <- min(max_lim_ctrl, max_lim_env)
  max_wo <- min(ctrl$max_work_ops, envelope_control$max_work_ops)

  admitted <- tryCatch({
    entries <- ipw_owner_safe_sum(entries_uncapped, cap, "complete retained IF inventory")
    cache_work <- ipw_owner_safe_sum(cache_work_uncapped, max_wo, "complete retained IF work")
    total <- ipw_owner_safe_sum(
      c(live, entries),
      max_lim_env,
      "pipeline simultaneous live entries plus retained point IF"
    )
    work_total <- ipw_owner_safe_sum(
      c(work, cache_work),
      envelope_control$max_work_ops,
      "pipeline cumulative work plus retained point IF"
    )
    ipw_owner_safe_sum(
      c(pre$peak_concurrent_entries, add_live, entries),
      max_lim_ctrl,
      "point IF producer plus plan record and retained point IF entries"
    )
    cov_entries <- ipw_owner_check_matrix_dims(
      L, L,
      min(ctrl$max_matrix_entries, ctrl$max_output_entries),
      "streamed L x L point IF covariance matrix"
    )
    ipw_owner_safe_sum(
      c(pre$peak_concurrent_entries, add_live, 4 * cov_entries, 8 * L_d,
        6 * pre$b_by_L_entries, 128, entries),
      ctrl$max_matrix_entries,
      "two-pass streamed point covariance plus retained point IF concurrent entries"
    )
    cov_verify_work <- if (one_verify_work > 0) {
      ipw_owner_safe_product(
        1 + 2 * (C_d + 1),
        one_verify_work,
        ctrl$max_work_ops,
        "covariance entry and two replay record verifications"
      )
    } else 0
    ipw_owner_safe_sum(
      c(
        ipw_owner_safe_product(2, pre$total_work_ops_bound, ctrl$max_work_ops, "two-pass point IF replay work"),
        cov_verify_work,
        ipw_owner_safe_product(N_d, 2 * L_d * L_d + 4 * L_d, ctrl$max_work_ops, "centering and crossprod accumulation work"),
        128,
        cache_work
      ),
      ctrl$max_work_ops,
      "two-pass streamed point covariance plus retained point IF work"
    )
    ipw_owner_safe_sum(
      c(pre$total_work_ops_bound, add_plan_work, cache_work),
      ctrl$max_work_ops,
      "point IF plan work plus retained point IF work"
    )
    list(entries = entries, cache_work = cache_work, total = total, work_total = work_total)
  }, cdid_ipw_owner_bound_error = function(e) {
    d <- e$diagnostics
    if (is.list(d) && scalar_nonneg(d$requested) && scalar_pos(d$bound) &&
        d$requested > d$bound) return(NULL)
    stop(e)
  })
  if (is.null(admitted)) return(NULL)

  carrier <- new.env(parent = emptyenv())
  carrier$kind <- "ipw_point_if_retained"
  carrier$plan_kind <- plan$kind
  carrier$state <- "reserved"
  carrier$plan_ref <- plan$ref
  carrier$frame_ref <- fo$ref
  carrier$panel_ref <- fo$panel_ref
  carrier$fitted_ref <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) plan$fitted_view$ref else NULL
  carrier$fitted_refs <- if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) plan$fitted_refs else NULL
  carrier$nuisance_refs <- if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) plan$nuisance_refs else NULL
  carrier$plan_snapshot <- plan
  lockBinding("plan_snapshot", carrier)
  carrier$N <- as.integer(N)
  carrier$n_total <- as.integer(n_total)
  carrier$clustered <- isTRUE(clustered)
  carrier$L <- as.integer(L)
  carrier$B <- as.integer(B)
  carrier$n_chunks <- as.integer(pre$n_chunks)
  carrier$link_ids <- plan$link_ids
  carrier$raw_entries <- N_d * L_d
  carrier$entries <- admitted$entries
  carrier$work <- admitted$cache_work
  carrier$admission_total <- admitted$total
  carrier$work_total <- admitted$work_total
  carrier$cap <- cap
  carrier$work_cap <- max_wo
  carrier$producer_passes <- 0L
  carrier$slice_passes <- 0L
  carrier$batches <- 0L
  carrier$slices_emitted <- 0L
  carrier$verifications <- 1L
  # Lock configuration bindings and retain their exact objects in a separate
  # locked environment. Unlock/rebind/relock of an equal COW copy is refused.
  seal <- new.env(parent = emptyenv())
  mutable <- c("state", "producer_passes", "slice_passes", "batches", "slices_emitted", "verifications")
  for (nm in setdiff(ls(carrier, all.names = TRUE), mutable)) {
    assign(nm, get(nm, carrier, inherits = FALSE), seal)
    lockBinding(nm, carrier)
  }
  lockEnvironment(seal, bindings = TRUE)
  carrier$seal <- seal
  lockBinding("seal", carrier)
  carrier
}

# Every boundary: the plan's own verifier (seal, controls, owner/nuisance/link/
# point-link bindings) and the reserved plan identity (the same sealed object,
# not a resealed or edited copy).
ipw_point_if_retained_check <- function(carrier, plan, control, verify_plan = TRUE, expected_state = NULL) {
  if (!is.environment(carrier) || !identical(carrier$kind, "ipw_point_if_retained") ||
      !identical(parent.env(carrier), emptyenv()) ||
      !is.character(carrier$state) || length(carrier$state) != 1L || is.na(carrier$state) ||
      !carrier$state %in% c("reserved", "filling", "filled", "slicing") ||
      (!is.null(expected_state) && !identical(carrier$state, expected_state))) {
    ipw_point_view_batch_abort("Retained point IF carrier is not available.")
  }
  seal <- carrier$seal
  if (!is.environment(seal) || !identical(parent.env(seal), emptyenv()) ||
      !environmentIsLocked(seal) || !bindingIsLocked("seal", carrier)) {
    ipw_point_view_batch_abort("Retained point IF carrier lost its configuration seal.")
  }
  for (nm in ls(seal, all.names = TRUE)) {
    if (!exists(nm, carrier, inherits = FALSE) || !bindingIsLocked(nm, carrier) ||
        !bindingIsLocked(nm, seal) ||
        !{
          current <- get(nm, carrier, inherits = FALSE)
          sealed <- get(nm, seal, inherits = FALSE)
          (is.null(current) && is.null(sealed)) || gmm_sexp_aliases(current, sealed)
        }) {
      ipw_point_view_batch_abort("Retained point IF carrier configuration changed.")
    }
  }
  if (isTRUE(verify_plan)) {
    ipw_point_if_retained_plan_verify(plan, control)
  }
  if (carrier$state %in% c("reserved", "filling") &&
      (exists("phi", carrier, inherits = FALSE) || exists("raw_ref", carrier, inherits = FALSE) ||
       !identical(carrier$producer_passes, 0L))) {
    ipw_point_view_batch_abort("Retained point IF carrier has raw storage outside its filled state.")
  }
  plan_L <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) as.integer(plan$batch$n_links) else as.integer(plan$n_links)
  same_refs <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) {
    identical(plan$fitted_view$ref, carrier$fitted_ref)
  } else if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
    identical(plan$fitted_refs, carrier$fitted_refs) && identical(plan$nuisance_refs, carrier$nuisance_refs)
  } else {
    identical(plan$kind, "cdid_custom_ipw_if_plan")
  }
  if (!exists("plan_snapshot", carrier, inherits = FALSE) || !bindingIsLocked("plan_snapshot", carrier) ||
      !gmm_sexp_aliases(plan, get("plan_snapshot", carrier, inherits = FALSE)) ||
      !identical(plan$kind, carrier$plan_kind) ||
      !identical(plan$ref, carrier$plan_ref) || !identical(plan$frame_owner$ref, carrier$frame_ref) ||
      !identical(plan$frame_owner$panel_ref, carrier$panel_ref) || !same_refs ||
      !identical(plan$n_a_units, carrier$N) || !identical(plan_L, carrier$L) ||
      !identical(plan$step, carrier$B) || !identical(plan$link_ids, carrier$link_ids)) {
    ipw_point_view_batch_abort("Retained point IF carrier does not describe this sealed plan.")
  }
  if (carrier$state %in% c("filled", "slicing")) {
    full <- carrier$phi; N <- carrier$N; L <- carrier$L
    if (!exists("phi", carrier, inherits = FALSE) || !bindingIsLocked("phi", carrier) ||
        !exists("raw_ref", carrier, inherits = FALSE) || !bindingIsLocked("raw_ref", carrier) ||
        !is.environment(carrier$raw_ref) || !identical(parent.env(carrier$raw_ref), emptyenv()) ||
        !environmentIsLocked(carrier$raw_ref) || !exists("phi", carrier$raw_ref, inherits = FALSE) ||
        !bindingIsLocked("phi", carrier$raw_ref) || !gmm_sexp_aliases(full, get("phi", carrier$raw_ref, inherits = FALSE)) ||
        !is.matrix(full) || !is.double(full) || nrow(full) != N || ncol(full) != L ||
        !identical(colnames(full), carrier$link_ids)) {
      ipw_point_view_batch_abort("Retained point IF carrier lost its shape after the producer pass.")
    }
  }
  carrier$verifications <- carrier$verifications + 1L
  invisible(TRUE)
}

# The one producer pass: the unchanged plan replay, whose every callback batch
# is stored raw (exact rows, IDs from the plan's bound panel owner, link order,
# finite values) and handed on unchanged. A locked raw-matrix binding retains
# the exact completed object; copy-on-modify value, label or order changes
# cannot alias it. The caller's covariance preflight admitted the replay work;
# the fill/slice work is charged in addition.
ipw_point_if_retained_fill <- function(
  carrier,
  plan,
  control,
  callback,
  use_nuisance_cache = TRUE,
  trace_env = NULL,
  .verified_entry = FALSE
) {
  completed <- FALSE
  owned_carrier <- is.environment(carrier) && identical(carrier$kind, "ipw_point_if_retained")
  on.exit(if (!completed && owned_carrier) ipw_point_if_retained_release(carrier), add = TRUE)
  ipw_point_if_retained_check(carrier, plan, control, verify_plan = !isTRUE(.verified_entry))
  if (!identical(carrier$state, "reserved")) {
    ipw_point_view_batch_abort("Retained point IF carrier can be filled only once, from its reserved state.")
  }
  if (!is.function(callback)) {
    cdid_ipw_point_view_abort("`callback` must be a function.", c("cdid_ipw_point_view_input_error", "cdid_ipw_input_error"))
  }
  ctrl <- plan$control; pre <- plan$preflight
  N <- carrier$N; L <- carrier$L; B <- carrier$B
  ipw_owner_safe_sum(c(pre$total_work_ops_bound, carrier$work), ctrl$max_work_ops, "CP08 retained point IF fill work")
  ipw_owner_check_matrix_dims(N, L, ctrl$max_matrix_entries, "CP08 retained point IF")
  po <- plan$panel_owner
  carrier$state <- "filling"
  phi <- matrix(0, N, L); filled <- 0L; batches <- 0L; cn <- NULL
  obs_cb <- function(batch) {
    rows <- seq.int(filled + 1L, min(N, filled + B))
    if (!is.list(batch) || !is.numeric(batch$a_index) || is.object(batch$a_index) ||
        !is.null(dim(batch$a_index)) ||
        !identical(as.double(batch$a_index), as.double(rows)) || !identical(batch$link_ids, plan$link_ids) ||
        !is.matrix(batch$phi) || !is.double(batch$phi) || nrow(batch$phi) != length(rows) || ncol(batch$phi) != L ||
        any(!is.finite(batch$phi)) || !identical(batch$id, po$a_id[rows]) || !is.null(batch$components)) {
      ipw_point_view_batch_abort("The producer chunk does not match the reserved retained point IF layout.")
    }
    if (filled == 0L) cn <<- colnames(batch$phi) else if (!identical(colnames(batch$phi), cn)) {
      ipw_point_view_batch_abort("The producer chunk link labels changed within one pass.")
    }
    ipw_point_if_retained_check(carrier, plan, control, expected_state = "filling")
    phi[rows, ] <<- batch$phi
    filled <<- filled + length(rows); batches <<- batches + 1L
    counts <- lapply(c("producer_passes", "slice_passes", "batches", "slices_emitted"),
      function(nm) get(nm, carrier, inherits = FALSE))
    callback(batch)
    ipw_point_if_retained_check(carrier, plan, control, expected_state = "filling")
    if (!identical(counts, lapply(c("producer_passes", "slice_passes", "batches", "slices_emitted"),
        function(nm) get(nm, carrier, inherits = FALSE)))) {
      ipw_point_view_batch_abort("Retained point IF carrier counters changed during the callback.")
    }
  }
  info <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) {
    replay_ipw_point_if_plan(plan, control, callback = obs_cb)
  } else if (identical(plan$kind, "cdid_ipw_point_if_general_plan")) {
    replay_ipw_point_if_general(
      plan = plan,
      control = control,
      callback = obs_cb,
      callback_schema = TRUE,
      use_nuisance_cache = use_nuisance_cache,
      trace_env = trace_env,
      state = NULL,
      carrier = NULL
    )
  } else if (identical(plan$kind, "cdid_custom_ipw_if_plan")) {
    replay_custom_ipw_if_plan(
      plan = plan,
      control = control,
      callback = obs_cb,
      callback_schema = TRUE,
      trace_env = trace_env,
      carrier = NULL
    )
  } else {
    ipw_point_view_batch_abort("Unsupported plan kind for retained point IF fill.")
  }
  # The producer pass verified `plan` after the last callback; check carrier snapshot/bindings here.
  ipw_point_if_retained_check(carrier, plan, control, verify_plan = FALSE, expected_state = "filling")
  if (filled != N || batches != as.integer(info$n_batches) || batches != carrier$n_chunks) {
    ipw_point_view_batch_abort("The producer pass did not cover every A unit exactly once.")
  }
  colnames(phi) <- cn
  raw_ref <- new.env(parent = emptyenv()); raw_ref$phi <- phi
  lockEnvironment(raw_ref, bindings = TRUE)
  carrier$phi <- phi; lockBinding("phi", carrier)
  carrier$raw_ref <- raw_ref; lockBinding("raw_ref", carrier)
  carrier$state <- "filled"
  carrier$producer_passes <- carrier$producer_passes + 1L
  carrier$batches <- batches
  ipw_point_if_retained_check(carrier, plan, control, expected_state = "filled")
  completed <- TRUE
  invisible(info)
}

# One slice pass over the plan's own row partition: the batch schema of the
# IF-only plan replay (phi, id, link_ids, a_index), the IDs from the plan's
# bound panel owner, with the plan and carrier re-verified at every boundary
# and after the last callback. No kernel call, no chunk query, no endpoint read.
replay_ipw_point_if_retained <- function(carrier, plan, control, callback, trace_env = NULL) {
  ipw_point_view_require_deps()
  if (missing(carrier) || missing(plan) || missing(control) || missing(callback)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`carrier`, `plan`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.function(callback)) {
    cdid_ipw_point_view_abort("`callback` must be a function.", c("cdid_ipw_point_view_input_error", "cdid_ipw_input_error"))
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  completed <- FALSE
  owned_carrier <- is.environment(carrier) && identical(carrier$kind, "ipw_point_if_retained")
  on.exit(if (!completed && owned_carrier) ipw_point_if_retained_release(carrier), add = TRUE)
  ipw_point_if_retained_check(carrier, plan, control, expected_state = "filled")
  if (!identical(carrier$state, "filled")) {
    ipw_point_view_batch_abort("Retained point IF carrier was not completely filled.")
  }
  ctrl <- plan$control; pre <- plan$preflight
  N <- carrier$N; L <- carrier$L; B <- carrier$B
  ipw_owner_safe_sum(c(pre$total_work_ops_bound, carrier$work), ctrl$max_work_ops, "CP08 retained point IF slice work")
  carrier$state <- "slicing"
  if (!is.null(trace_env)) {
    trace_env$preflight_completed <- TRUE
    trace_env$batches_emitted <- 0L
    trace_env$units_emitted <- 0L
    trace_env$active_query_chunks <- 0L
    trace_env$max_simultaneous_query_chunks <- 0L
    trace_env$active_nuisance_caches <- 0L
    trace_env$max_simultaneous_nuisance_caches <- 0L
    trace_env$nuisance_cache_hits <- 0L
    trace_env$endpoints_read <- 0L
  }
  verifications <- 1L
  emitted <- 0L; batches <- 0L; max_batch_units <- 0L
  a_start <- 1
  while (a_start <= N) {
    ipw_owner_check_interrupt()
    a_end <- min(as.double(N), a_start + B - 1)
    full <- carrier$phi
    a_seq <- seq.int(as.integer(a_start), as.integer(a_end))
    batches <- batches + 1L
    ipw_point_if_retained_check(carrier, plan, control, expected_state = "slicing")
    counts <- lapply(c("producer_passes", "slice_passes", "batches", "slices_emitted"),
      function(nm) get(nm, carrier, inherits = FALSE))
    callback(list(
      phi = full[a_seq, , drop = FALSE],
      id = plan$panel_owner$a_id[a_seq],
      link_ids = plan$link_ids,
      a_index = a_seq
    ))
    # Match the original replay's post-callback source check, including the
    # last callback; the filled carrier's raw identity is checked here too.
    ipw_point_if_retained_check(carrier, plan, control, expected_state = "slicing")
    if (!identical(counts, lapply(c("producer_passes", "slice_passes", "batches", "slices_emitted"),
        function(nm) get(nm, carrier, inherits = FALSE)))) {
      ipw_point_view_batch_abort("Retained point IF carrier counters changed during the callback.")
    }
    verifications <- verifications + 2L
    b_len <- length(a_seq)
    emitted <- emitted + b_len
    max_batch_units <- max(max_batch_units, b_len)
    if (!is.null(trace_env)) {
      trace_env$batches_emitted <- batches
      trace_env$units_emitted <- emitted
    }
    a_start <- a_end + 1
  }
  if (emitted != N) {
    cdid_ipw_point_view_abort("Point IF replay unit count mismatch.", c("cdid_ipw_count_error"))
  }
  carrier$slice_passes <- carrier$slice_passes + 1L
  carrier$slices_emitted <- carrier$slices_emitted + batches
  carrier$state <- "filled"
  completed <- TRUE
  invisible(list(
    n_a_units = emitted,
    n_links = L,
    link_ids = plan$link_ids,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * L,
    actual_work_ops = 0,
    kernel_calls = 0L,
    normalizations_performed = 0L,
    verifications = verifications,
    state_verifications = 0L,
    served_from = "retained_point_if",
    preflight = pre
  ))
}

# Release right after the GMM solve: the matrix, its locked reference and the
# plan snapshot are dropped; the scalar record of what happened is returned
# (NULL carrier: the plan replay served every pass).
ipw_point_if_retained_release <- function(carrier, plan = NULL) {
  if (is.null(carrier)) {
    if (is.list(plan)) {
      N_d <- as.double(plan$n_a_units)
      L_d <- if (identical(plan$kind, "cdid_ipw_point_if_batch_plan")) as.double(plan$batch$n_links) else as.double(plan$n_links)
      B_d <- as.double(plan$step)
      return(list(
        status = "streamed",
        admitted = FALSE,
        released = FALSE,
        plan_kind = plan$kind,
        n_a_units = plan$n_a_units,
        n_links = as.integer(L_d),
        chunk_rows = plan$step,
        n_chunks = as.integer(plan$preflight$n_chunks),
        raw_if_entries = N_d * L_d,
        reserved_entries = 0,
        work_ops_bound = 0,
        producer_passes = 0L,
        slice_passes = 0L
      ))
    }
    return(list(status = "streamed"))
  }
  was_filled <- identical(carrier$state, "filled")
  rec <- list(
    status = if (was_filled) "retained" else paste0("incomplete:", carrier$state),
    admitted = TRUE,
    released = TRUE,
    plan_kind = carrier$plan_kind,
    n_a_units = carrier$N,
    n_total = carrier$n_total,
    clustered = carrier$clustered,
    n_links = carrier$L,
    chunk_rows = carrier$B,
    n_batches = carrier$batches,
    producer_passes = carrier$producer_passes,
    slice_passes = carrier$slice_passes,
    slices_emitted = carrier$slices_emitted,
    verifications = carrier$verifications,
    raw_if_entries = carrier$raw_entries,
    reserved_entries = carrier$entries,
    admission_total_entries = carrier$admission_total,
    cap_entries = carrier$cap,
    work_ops_bound = carrier$work,
    admission_total_work_ops = carrier$work_total,
    cap_work_ops = carrier$work_cap,
    lifetime = paste(
      "reserved before allocation inside the complete simultaneous envelope; filled by the unchanged",
      "first covariance replay; sliced for the covariance Gram pass and both GMM projection passes in the plan's",
      "row partition; released before auxiliary reporting and final output"
    )
  )
  drop <- intersect(c("phi", "raw_ref", "plan_snapshot", "seal", "plan_ref", "frame_ref", "panel_ref",
    "link_ids", "fitted_ref", "fitted_refs", "nuisance_refs"),
    ls(carrier, all.names = TRUE))
  for (nm in drop) if (bindingIsLocked(nm, carrier)) unlockBinding(nm, carrier)
  rm(list = drop, envir = carrier)
  carrier$state <- "released"
  rec
}

# The two-pass centered covariance of `stream_ipw_point_covariance_view`, with
# its gates, column sums, centering and crossproduct accumulation unchanged,
# over the plan's normalization instead of three repeated normalizations.
# CP08 retained point IF: a reserved `carrier` is filled by the first (mean)
# replay and sliced for the second (Gram) replay; NULL keeps both replays.
stream_ipw_point_covariance_plan <- function(plan, control, carrier = NULL) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control)) {
    cdid_ipw_point_view_abort("Required arguments (`plan`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error")
  }
  ipw_point_if_plan_verify(plan, control)
  pre1 <- plan$preflight
  ctrl <- plan$control
  n_a <- as.double(pre1$n_a_units)
  L <- as.integer(pre1$n_links)
  link_ids <- plan$link_ids
  carrier_entries <- 0; carrier_work <- 0
  if (!is.null(carrier)) {
    ipw_point_if_retained_check(carrier, plan, control, verify_plan = FALSE)
    if (!identical(carrier$state, "reserved"))
      ipw_point_view_batch_abort("Retained point IF carrier must be reserved and unfilled for the covariance stream.")
    carrier_entries <- carrier$entries; carrier_work <- carrier$work
  }

  cov_entries <- ipw_owner_check_matrix_dims(
    L, L,
    max_matrix_entries = min(ctrl$max_matrix_entries, ctrl$max_output_entries),
    label = "streamed L x L point IF covariance matrix"
  )
  ipw_owner_safe_sum(
    c(pre1$peak_concurrent_entries, plan$batch$additional_live_entries, 4 * cov_entries, 8 * L,
      6 * pre1$b_by_L_entries, 128, carrier_entries),
    ctrl$max_matrix_entries,
    "two-pass streamed point covariance concurrent entries"
  )
  ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(2, pre1$total_work_ops_bound, ctrl$max_work_ops, "two-pass point IF replay work"),
      ipw_owner_safe_product(n_a, 2 * L * L + 4 * L, ctrl$max_work_ops, "centering and crossprod accumulation work"),
      128, carrier_work
    ),
    ctrl$max_work_ops,
    "total two-pass streamed point covariance work"
  )

  sum_phi <- numeric(L)
  mean_callback <- function(batch) {
    sum_phi <<- sum_phi + colSums(batch$phi)
    if (any(!is.finite(sum_phi))) {
      cdid_abort("Nonfinite IF sum during covariance replay.", "cdid_ipw_numeric_error")
    }
  }
  if (is.null(carrier)) {
    replay_ipw_point_if_plan(plan = plan, control = control, callback = mean_callback)
  } else {
    ipw_point_if_retained_fill(carrier, plan, control, callback = mean_callback, .verified_entry = TRUE)
  }
  mean_phi <- sum_phi / n_a
  names(mean_phi) <- link_ids

  cov_mat <- matrix(0.0, nrow = L, ncol = L, dimnames = list(link_ids, link_ids))
  gram_callback <- function(batch) {
    centered <- sweep(batch$phi, 2L, mean_phi, FUN = "-")
    if (any(!is.finite(centered))) {
      cdid_abort("Nonfinite centered IF during covariance replay.", "cdid_ipw_numeric_error")
    }
    cov_mat <<- cov_mat + crossprod(centered)
    if (any(!is.finite(cov_mat))) {
      cdid_abort("Nonfinite IF crossproduct or sum during covariance replay.", "cdid_ipw_numeric_error")
    }
    rm(centered)
  }
  replay_info <- if (is.null(carrier)) {
    replay_ipw_point_if_plan(plan = plan, control = control, callback = gram_callback)
  } else {
    replay_ipw_point_if_retained(carrier, plan, control, callback = gram_callback)
  }
  cov_mat <- cov_mat / n_a

  list(
    covariance = cov_mat,
    mean_phi = mean_phi,
    estimates = stats::setNames(
      vapply(plan$point_links, `[[`, numeric(1L), "estimate"),
      link_ids
    ),
    n_a_units = pre1$n_a_units,
    replay = replay_info
  )
}

# ---------------------------------------------------------------------------
# Section 3: Candidate-to-Certificate Composition (Root, Link, and IF Binding)
# ---------------------------------------------------------------------------

# Outward root, link and IF enclosures cover the unclipped functional. A
# nuisance whose empirical rates are floored or whose control odds are capped
# is refused rather than certified for a different functional.
ipw_point_view_refuse_clipped <- function(fitted_view) {
  clip <- fitted_view$diagnostics$clipping
  record_ok <- (is.list(clip) || is.numeric(clip)) && !is.object(clip) &&
    is.null(dim(clip)) && !is.null(names(clip)) && !anyNA(names(clip)) &&
    !anyDuplicated(names(clip)) &&
    all(c("r_floored", "p_capped_ctrl") %in% names(clip))
  inactive_count <- function(name) {
    if (!record_ok) return(FALSE)
    count <- clip[[name]]
    is.numeric(count) && !is.object(count) && is.null(dim(count)) &&
      length(count) == 1L && is.finite(count) && count == 0
  }
  if (!inactive_count("r_floored") || !inactive_count("p_capped_ctrl")) {
    cdid_ipw_point_view_abort(
      "Probability clipping binds for this nuisance; outward enclosures cover only the unclipped functional, so certification is refused.",
      "cdid_ipw_point_view_clipping_error"
    )
  }
  invisible(TRUE)
}

bind_fitted_nuisance_to_problem_view <- function(
  fitted_view,
  F,
  r,
  certificate_control,
  problem_control,
  reader_control,
  beta0 = NULL
) {
  ipw_point_view_require_deps()
  if (missing(fitted_view) || missing(F) || missing(r) ||
      missing(certificate_control) || missing(problem_control) || missing(reader_control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`fitted_view`, `F`, `r`, `certificate_control`, `problem_control`, `reader_control`) must be supplied explicitly (no automatic radius or default controls).",
      "cdid_ipw_point_view_input_error"
    )
  }
  rd_ctrl <- validate_ipw_owner_control(reader_control)
  ipw_point_view_payload_budget(list(fitted_view, F, r, beta0,
    attr(reader_control, "expected_a_id", exact = TRUE)), rd_ctrl)
  validate_ipw_fitted_nuisance_view(fitted_view)
  ipw_point_view_refuse_clipped(fitted_view)
  ipw_owner_refuse_decorative_attrs(F, "F")
  ipw_owner_refuse_decorative_attrs(r, "r")
  if (!is.null(beta0)) {
    ipw_owner_refuse_decorative_attrs(beta0, "beta0")
  }

  exact_beta0 <- as.double(fitted_view$coefficients)
  if (!is.null(beta0)) {
    if (!is.atomic(beta0) || !is.null(dim(beta0)) || is.object(beta0) ||
        !is.numeric(beta0) || !ipw_point_view_identical(as.double(beta0), exact_beta0)) {
      cdid_ipw_point_view_abort(
        "Supplied `beta0` disagrees bitwise with `as.double(fitted_view$coefficients)`; candidate substitution is refused.",
        c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_beta_error")
      )
    }
  }

  nv_obj <- get("nuisance_view", envir = fitted_view$nuisance_view_ref, inherits = FALSE)
  frame_owner <- get("frame_owner", envir = fitted_view$frame_ref, inherits = FALSE)
  validate_ipw_nuisance_view(nv_obj, frame_owner = frame_owner)

  problem <- construct_certified_ipw_nuisance_problem_view(
    nuisance_view = nv_obj,
    beta0 = exact_beta0,
    F = F,
    r = r,
    certificate_control = certificate_control,
    problem_control = problem_control,
    reader_control = reader_control
  )

  if (!ipw_point_view_identical(problem$raw$beta0, exact_beta0) ||
      !ipw_point_view_identical(problem$frame_ref, fitted_view$frame_ref) ||
      !ipw_point_view_identical(problem$panel_ref, fitted_view$panel_ref) ||
      !ipw_point_view_identical(problem$nuisance_ref, fitted_view$nuisance_view_ref)) {
    cdid_ipw_point_view_abort(
      "Constructed nuisance problem view failed post-binding identity or exact `beta0` verification.",
      "cdid_ipw_point_view_integrity_error"
    )
  }
  # Selected-radius identity: with the typed automatic request the stored
  # problem radius must be the certificate's one selected scalar; with a manual
  # radius it must be the caller's value. The request object is never compared
  # as a number.
  if (!ipw_point_view_identical(problem$raw$r, problem$certificate$r) ||
      !ipw_point_view_identical(problem$raw$r, problem$certificate$provenance$radius_r) ||
      (logit_cert_is_automatic_radius_request(r) &&
         !identical(problem$certificate$radius_selection$mode, .logit_cert_automatic_radius_mode)) ||
      (!logit_cert_is_automatic_radius_request(r) &&
         !ipw_point_view_identical(problem$raw$r, as.double(r) + 0.0))) {
    cdid_ipw_point_view_abort(
      "Constructed nuisance problem view failed selected-radius identity verification.",
      "cdid_ipw_point_view_integrity_error"
    )
  }

  problem
}

bind_point_link_to_enclosure_view <- function(
  problem,
  point_link,
  link_control,
  reader_control,
  loading_beta_candidate = NULL,
  include_candidates = FALSE
) {
  ipw_point_view_require_deps()
  if (missing(problem) || missing(point_link) || missing(link_control) || missing(reader_control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`problem`, `point_link`, `link_control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!gmm_scalar_logical(include_candidates)) {
    cdid_ipw_point_view_abort("`include_candidates` must be a scalar logical.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(loading_beta_candidate)) {
    ipw_owner_refuse_decorative_attrs(loading_beta_candidate, "loading_beta_candidate")
  }

  budget_ctrl <- validate_ipw_link_enclosure_control(link_control)
  ipw_point_view_payload_budget(list(problem, point_link, loading_beta_candidate,
    attr(reader_control, "expected_a_id", exact = TRUE)), budget_ctrl)
  ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE, expected_storage_mode = "compact_view")
  validate_ipw_point_link_view(point_link)
  if (is.environment(point_link$fitted_nuisance_ref))
    ipw_point_view_refuse_clipped(get("fitted_nuisance_view", envir = point_link$fitted_nuisance_ref, inherits = FALSE))

  if (!ipw_point_view_identical(problem$frame_ref, point_link$frame_ref) ||
      !ipw_point_view_identical(problem$panel_ref, point_link$panel_ref) ||
      !ipw_point_view_identical(problem$nuisance_ref, point_link$nuisance_view_ref)) {
    cdid_ipw_point_view_abort(
      "Certified nuisance problem and point link do not share the exact same sealed `frame_owner`, `panel_owner`, and `nuisance_view`.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_enclosure_problem_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  if (!ipw_point_view_identical(problem$raw$beta0, as.double(point_link$coefficients))) {
    cdid_ipw_point_view_abort(
      "Certified nuisance problem `raw$beta0` disagrees bitwise with `as.double(point_link$coefficients)`.",
      c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_beta_error")
    )
  }
  if (!ipw_point_view_identical(problem$raw$cohorts, point_link$cohorts$cohorts) ||
      !isTRUE(point_link$nuisance_mask_verified)) {
    cdid_ipw_point_view_abort(
      "Point link cohorts or observer-mask verification disagree with the nuisance problem.",
      c("cdid_ipw_owner_cohort_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  actual_loading <- as.double(point_link$loadings$loading_beta)
  if (!is.null(loading_beta_candidate)) {
    if (!is.atomic(loading_beta_candidate) || !is.null(dim(loading_beta_candidate)) ||
        is.object(loading_beta_candidate) || !is.numeric(loading_beta_candidate) ||
        !ipw_point_view_identical(as.double(loading_beta_candidate), actual_loading)) {
      cdid_ipw_point_view_abort(
        "Supplied `loading_beta_candidate` disagrees bitwise with `as.double(point_link$loadings$loading_beta)`.",
        c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error")
      )
    }
  }

  link_view <- get("link_view", envir = point_link$link_view_ref, inherits = FALSE)
  if (!is.null(link_view$loading_beta_candidate) &&
      !ipw_point_view_identical(as.double(link_view$loading_beta_candidate), actual_loading)) {
    cdid_ipw_point_view_abort(
      "Candidate `loading_beta_candidate` on `link_view` disagrees with `point_link$loadings$loading_beta`.",
      c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error")
    )
  }

  cand_arg <- if (isTRUE(include_candidates)) {
    list(
      mu = as.double(point_link$estimate),
      treated_mean = as.double(point_link$treated_mean),
      control_mean = as.double(point_link$control_mean),
      G_beta = as.double(point_link$partials$G_beta),
      G_r = as.double(point_link$partials$G_r),
      A_beta_r = unname(point_link$partials$A_beta_r),
      H_chol = unname(point_link$factors$H_chol),
      loading_beta = actual_loading,
      loading_rate = as.double(point_link$loadings$loading_rate)
    )
  } else NULL

  enclosed <- enclose_ipw_link_view(
    problem = problem,
    link_view = link_view,
    loading_beta_candidate = actual_loading,
    control = link_control,
    reader_control = reader_control,
    candidates = cand_arg
  )

  if (!ipw_point_view_identical(enclosed$loading_beta_candidate, actual_loading) ||
      !ipw_point_view_identical(enclosed$link_ref, point_link$link_view_ref) ||
      !ipw_point_view_identical(enclosed$frame_ref, point_link$frame_ref) ||
      !ipw_point_view_identical(enclosed$panel_ref, point_link$panel_ref) ||
      !ipw_point_view_identical(enclosed$nuisance_ref, point_link$nuisance_view_ref)) {
    cdid_ipw_point_view_abort(
      "Enclosed link view failed post-binding identity or loading candidate verification.",
      "cdid_ipw_point_view_integrity_error"
    )
  }

  enclosed
}

enclose_point_link_if_view_chunk <- function(
  link_enclosure,
  point_link,
  a_start,
  a_end,
  point_control,
  link_control,
  reader_control,
  include_components = FALSE,
  include_shared_xi = FALSE,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(link_enclosure) || missing(point_link) || missing(a_start) ||
      missing(a_end) || missing(point_control) || missing(link_control) ||
      missing(reader_control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`link_enclosure`, `point_link`, `a_start`, `a_end`, `point_control`, `link_control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!gmm_scalar_logical(include_components) || !gmm_scalar_logical(include_shared_xi)) {
    cdid_ipw_point_view_abort(
      "`include_components` and `include_shared_xi` must be scalar logicals.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_point_view_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }

  pt_ctrl <- validate_ipw_point_view_control(point_control)
  lk_ctrl <- validate_ipw_link_enclosure_control(link_control)
  rd_ctrl <- validate_ipw_owner_control(reader_control)
  ipw_point_view_payload_budget(list(link_enclosure, point_link,
    attr(reader_control, "expected_a_id", exact = TRUE)), pt_ctrl)
  prob <- ipw_validate_link_enclosure_mode(link_enclosure, require_certified = FALSE, expected_storage_mode = "compact_view")
  validate_ipw_point_link_view(point_link)

  if (!ipw_point_view_identical(link_enclosure$frame_ref, point_link$frame_ref) ||
      !ipw_point_view_identical(link_enclosure$panel_ref, point_link$panel_ref) ||
      !ipw_point_view_identical(link_enclosure$nuisance_ref, point_link$nuisance_view_ref) ||
      !ipw_point_view_identical(link_enclosure$link_ref, point_link$link_view_ref)) {
    cdid_ipw_point_view_abort(
      "Link enclosure and point link do not reference the exact same sealed `frame_owner`, `panel_owner`, `nuisance_view`, and `link_view`.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_enclosure_link_error", "cdid_ipw_shared_incompatible_error")
    )
  }
  if (!ipw_point_view_identical(prob$raw$beta0, as.double(point_link$coefficients)) ||
      !ipw_point_view_identical(link_enclosure$loading_beta_candidate, as.double(point_link$loadings$loading_beta))) {
    cdid_ipw_point_view_abort(
      "Link enclosure `beta0` or `loading_beta_candidate` disagrees with `point_link`.",
      c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_candidate_error")
    )
  }

  na <- as.integer(prob$raw$n_a_units)
  p <- as.integer(prob$raw$p_covariates)
  m <- length(prob$raw$cohorts)
  s_idx <- ipw_owner_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_owner_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > na) {
    cdid_ipw_point_view_abort(
      paste0("Chunk range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= ", na, "`."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  b_len <- as.integer(e_idx - s_idx + 1)

  # Preflight BOTH the point IF chunk evaluation AND the interval IF chunk enclosure BEFORE reading endpoints
  point_pre <- preflight_ipw_point_if_view(
    point_links = list(point_link),
    control = pt_ctrl,
    include_components = include_components,
    chunk_rows = b_len
  )
  panel_owner <- get("panel_owner", envir = prob$panel_ref, inherits = FALSE)
  expected_a_id <- attr(reader_control, "expected_a_id", exact = TRUE)
  ipw_point_view_payload_budget(list(expected_a_id), pt_ctrl)
  interval_pre <- preflight_ipw_link_if_view_chunk(
    b_len = b_len,
    p = p,
    m = m,
    control = lk_ctrl,
    reader_control = rd_ctrl,
    certificate_control = prob$controls$certificate_control,
    include_components = include_components,
    include_shared_xi = include_shared_xi,
    has_phi_candidate = TRUE,
    has_xi_beta = isTRUE(include_shared_xi),
    read_outcomes = TRUE,
    panel_owner = panel_owner,
    a_seq = seq.int(as.integer(s_idx), as.integer(e_idx)),
    expected_a_id = expected_a_id,
    caller_phi_entries = b_len,
    caller_xi_entries = 0
  )

  ipw_owner_safe_sum(c(point_pre$peak_concurrent_entries,
    interval_pre$combined_concurrent_entries, 8*b_len + 4*p*p),
    pt_ctrl$max_matrix_entries, "point and interval candidate/result composition entries")
  ipw_owner_safe_sum(c(point_pre$total_work_ops_bound,
    interval_pre$combined_work_ops_bound), pt_ctrl$max_work_ops,
    "point and interval candidate/result whole work")

  pt_batch <- evaluate_ipw_point_if_view_chunk(
    point_links = list(point_link),
    a_start = s_idx,
    a_end = e_idx,
    control = pt_ctrl,
    include_components = include_components,
    callback_schema = FALSE,
    use_nuisance_cache = FALSE,
    trace_env = trace_env
  )
  phi_cand <- as.double(pt_batch$phi[, 1L])
  H_chol_cand <- unname(point_link$factors$H_chol)
  pt_batch <- NULL

  ans <- enclose_ipw_link_if_view_chunk(
    link_enclosure = link_enclosure,
    a_start = s_idx,
    a_end = e_idx,
    control = link_control,
    reader_control = reader_control,
    phi_candidate = phi_cand,
    H_chol_candidate = H_chol_cand,
    include_components = include_components,
    include_shared_xi = include_shared_xi
  )
  rm(phi_cand, H_chol_cand)
  ans
}

# The diagnostic binding wrapper retains all chunk summaries. Charge their
# full count and every repeated producer/enclosure traversal before any outcomes.
# Bounds deliberately use fixed component envelopes where candidates do not yet
# exist. This is a computational cap, never covariance/inference acceptance.
preflight_ipw_point_certificate_binding <- function(
  fitted_view, link_view, point_control, certificate_control, problem_control,
  link_control, reader_control, enclose_if_chunks, include_if_components,
  include_shared_xi, F, r, problem = NULL
) {
  pt <- validate_ipw_point_view_control(point_control)
  lc <- validate_ipw_link_enclosure_control(link_control)
  rc <- validate_ipw_owner_control(reader_control)
  cc <- validate_logit_certificate_control(certificate_control)
  pc <- validate_ipw_nuisance_problem_control(problem_control)
  extra <- ipw_point_view_payload_budget(list(F, r, problem,
    attr(reader_control, "expected_a_id", exact = TRUE)), pt)
  pp <- preflight_ipw_point_link_view(fitted_view, link_view, pt)
  p <- pp$p_covariates; m <- pp$m_cohorts; n <- pp$n_a_units
  b <- min(n, pt$max_chunk_rows, lc$max_chunk_rows, rc$max_chunk_rows)
  chunks <- if (enclose_if_chunks) ceiling(n / b) else 0
  # Future compact point metadata: every value/name/attribute/byte charged.
  # payload_budget uses <=640 operations per retained charged slot; normalization
  # and wrapper repetitions below are deliberately counted separately.
  future_payload <- pp$retained_entries_per_link
  future_normalize_work <- 16 + 8 * 640 * future_payload
  future_if_work <- 4 * future_normalize_work +
    n * (256 * (p + m + 4)) + ceiling(n / b) * 8192 + 4096
  summary_entries <- ipw_owner_safe_product(1024, chunks,
    min(pt$max_retained_entries, pt$max_output_entries), "all binding IF chunk summaries")
  retained <- ipw_owner_safe_sum(c(pp$total_retained_entries, summary_entries,
    32768 + 128 * (p * p + p * m + p + m), extra$entries),
    min(pt$max_retained_entries, pt$max_output_entries), "complete binding retained output")
  ipre <- preflight_ipw_link_if_view_chunk(b, p, m, lc, rc, cc,
    include_components = include_if_components, include_shared_xi = include_shared_xi,
    has_phi_candidate = TRUE, has_xi_beta = include_shared_xi, read_outcomes = TRUE,
    caller_phi_entries = b, caller_xi_entries = 0)
  # Existing component caps bound root and link phases; add caller-held outputs
  # and both point/interval results. Base owners are already in pp's envelope.
  peak <- ipw_owner_safe_sum(c(retained, pp$peak_concurrent_entries,
    max(pc$max_matrix_entries, lc$max_matrix_entries),
    ipre$combined_concurrent_entries, b * (256 + 16 * p) + 8192),
    pt$max_matrix_entries, "complete binding concurrent output and workspaces")
  total_work <- ipw_owner_safe_sum(c(pp$total_work_ops_bound,
    if (is.null(problem)) pc$max_work_ops else 0, lc$max_work_ops,
    future_if_work, chunks * (2 * future_if_work + ipre$combined_work_ops_bound),
    4 * extra$work, 1024 * chunks), pt$max_work_ops, "complete binding whole-traversal work")
  list(n_if_chunks = chunks, chunk_rows = b, summary_entries = summary_entries,
    retained_entries = retained, peak_concurrent_entries = peak,
    future_point_if_work_bound = future_if_work,
    total_work_ops_bound = total_work, component_caps_conservatively_reserved = TRUE,
    covariance_or_complete_inference_acceptance = FALSE)
}

certify_ipw_point_link_view <- function(
  fitted_view,
  link_view,
  F,
  r,
  point_control,
  certificate_control,
  problem_control,
  link_control,
  reader_control,
  problem = NULL,
  enclose_if_chunks = TRUE,
  include_if_components = FALSE,
  include_shared_xi = FALSE,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(fitted_view) || missing(link_view) || missing(F) || missing(r) ||
      missing(point_control) || missing(certificate_control) ||
      missing(problem_control) || missing(link_control) || missing(reader_control)) {
    cdid_ipw_point_view_abort(
      "All required arguments (`fitted_view`, `link_view`, `F`, `r`, `point_control`, `certificate_control`, `problem_control`, `link_control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!gmm_scalar_logical(enclose_if_chunks) ||
      !gmm_scalar_logical(include_if_components) ||
      !gmm_scalar_logical(include_shared_xi)) {
    cdid_ipw_point_view_abort(
      "`enclose_if_chunks`, `include_if_components`, and `include_shared_xi` must be scalar logicals.",
      "cdid_ipw_point_view_input_error"
    )
  }

  binding_preflight <- preflight_ipw_point_certificate_binding(
    fitted_view, link_view, point_control, certificate_control, problem_control,
    link_control, reader_control, enclose_if_chunks, include_if_components,
    include_shared_xi, F, r, problem)

  # 2. Bind or verify the compact nuisance problem view
  if (is.null(problem)) {
    prob <- bind_fitted_nuisance_to_problem_view(
      fitted_view = fitted_view,
      F = F,
      r = r,
      certificate_control = certificate_control,
      problem_control = problem_control,
      reader_control = reader_control
    )
  } else {
    ipw_validate_nuisance_problem_mode(problem, require_certified = FALSE, expected_storage_mode = "compact_view")
    if (!ipw_point_view_identical(problem$controls$certificate_control, validate_logit_certificate_control(certificate_control)) ||
        !ipw_point_view_identical(problem$controls$problem_control, validate_ipw_nuisance_problem_control(problem_control)) ||
        !ipw_point_view_identical(problem$controls$reader_control, validate_ipw_owner_control(reader_control)) ||
        !ipw_point_view_identical(problem$frame_ref, fitted_view$frame_ref) ||
        !ipw_point_view_identical(problem$panel_ref, fitted_view$panel_ref) ||
        !ipw_point_view_identical(problem$nuisance_ref, fitted_view$nuisance_view_ref) ||
        !ipw_point_view_identical(problem$raw$beta0, as.double(fitted_view$coefficients)) ||
        !ipw_point_view_identical(problem$raw$F, matrix(as.double(F), nrow = nrow(F), ncol = ncol(F))) ||
        # Caller radius identity: a typed automatic request matches a problem
        # whose stored radius is the certificate's selected scalar; a manual
        # radius must equal the stored radius bitwise. The request object is
        # never forwarded or compared as a numeric radius.
        !(if (logit_cert_is_automatic_radius_request(r)) {
            identical(problem$certificate$radius_selection$mode, .logit_cert_automatic_radius_mode) &&
              ipw_point_view_identical(problem$raw$r, problem$certificate$r) &&
              ipw_point_view_identical(problem$raw$r, problem$certificate$provenance$radius_r)
          } else {
            gmm_scalar_numeric(r, positive = TRUE) && ipw_point_view_identical(problem$raw$r, as.double(r))
          })) {
      cdid_ipw_point_view_abort(
        "Supplied `problem` disagrees with `fitted_view`, `F`, or `r`.",
        c("cdid_ipw_point_view_candidate_mismatch_error", "cdid_ipw_enclosure_problem_mismatch_error")
      )
    }
    prob <- problem
  }

  # 1. Fit the compact point link (no propensity refit)
  point_link <- fit_ipw_point_link_view(
    fitted_view = fitted_view,
    link_view = link_view,
    control = point_control,
    trace_env = trace_env
  )

  # 3. Bind the point link's loading candidate to the compact link enclosure validator
  enclosed_link <- bind_point_link_to_enclosure_view(
    problem = prob,
    point_link = point_link,
    link_control = link_control,
    reader_control = reader_control
  )

  # 4. Prepare covariance callback spec and optionally stream complete-IF chunk enclosures
  cand_spec <- make_ipw_point_candidate_rows_spec(
    point_links = list(point_link),
    control = point_control,
    trace_env = trace_env
  )

  na <- as.integer(point_link$counts$n_a_units)
  pt_ctrl <- validate_ipw_point_view_control(point_control)
  lk_ctrl <- validate_ipw_link_enclosure_control(link_control)
  rd_ctrl <- validate_ipw_owner_control(reader_control)
  step <- as.integer(min(na, pt_ctrl$max_chunk_rows, lk_ctrl$max_chunk_rows, rd_ctrl$max_chunk_rows))

  if_chunk_summaries <- list()
  if_all_certified <- isTRUE(enclosed_link$certified)
  if_reasons <- if (!isTRUE(enclosed_link$certified)) enclosed_link$inconclusive_reasons else character(0L)
  max_eval_err <- 0.0
  max_root_err <- 0.0
  max_direct_star_err <- 0.0
  total_if_rows_certified <- 0L

  if (isTRUE(enclose_if_chunks)) {
    cur_s <- 1L
    while (cur_s <= na) {
      cur_e <- as.integer(min(as.double(na), as.double(cur_s) + step - 1))
      ans <- enclose_point_link_if_view_chunk(
        link_enclosure = enclosed_link,
        point_link = point_link,
        a_start = cur_s,
        a_end = cur_e,
        point_control = pt_ctrl,
        link_control = link_control,
        reader_control = reader_control,
        include_components = include_if_components,
        include_shared_xi = include_shared_xi,
        trace_env = trace_env
      )
      sep <- ans$separated_row_errors
      ch_ev <- if (!is.null(sep$eval_error)) max(sep$eval_error) else NA_real_
      ch_rt <- if (!is.null(sep$root_error)) max(sep$root_error) else NA_real_
      ch_ds <- if (!is.null(sep$direct_star_error)) max(sep$direct_star_error) else NA_real_
      if (isTRUE(ans$certified)) {
        total_if_rows_certified <- as.integer(total_if_rows_certified + as.integer(cur_e - cur_s + 1L))
        max_eval_err <- max(max_eval_err, ch_ev)
        max_root_err <- max(max_root_err, ch_rt)
        max_direct_star_err <- max(max_direct_star_err, ch_ds)
      } else {
        if_all_certified <- FALSE
        if_reasons <- unique(c(if_reasons, ans$inconclusive_reasons))
      }
      if_chunk_summaries[[length(if_chunk_summaries) + 1L]] <- list(
        a_start = cur_s,
        a_end = cur_e,
        certified = isTRUE(ans$certified),
        status = ans$status,
        inconclusive_reasons = ans$inconclusive_reasons,
        max_eval_error = ch_ev,
        max_root_error = ch_rt,
        max_direct_star_error = ch_ds,
        work = ans$work
      )
      rm(ans, sep)
      cur_s <- as.integer(cur_e + 1L)
    }
  } else {
    if_all_certified <- FALSE
  }

  overall_certified <- isTRUE(prob$certified) && isTRUE(enclosed_link$certified) &&
    (!isTRUE(enclose_if_chunks) || isTRUE(if_all_certified))
  all_reasons <- unique(c(
    prob$inconclusive_reasons,
    enclosed_link$inconclusive_reasons,
    if_reasons
  ))

  structure(
    list(
      certified = overall_certified,
      status = if (overall_certified) "CERTIFIED" else "INCONCLUSIVE",
      inconclusive_reasons = all_reasons,
      point_link = point_link,
      problem = prob,
      link_enclosure = enclosed_link,
      candidate_rows_spec = cand_spec,
      binding_preflight = binding_preflight,
      if_chunk_summaries = if_chunk_summaries,
      if_summary = list(
        enclose_if_chunks = isTRUE(enclose_if_chunks),
        if_chunks_certified = isTRUE(if_all_certified),
        n_chunks = length(if_chunk_summaries),
        chunk_rows = step,
        n_if_rows_certified = total_if_rows_certified,
        max_eval_error = if (overall_certified && isTRUE(enclose_if_chunks)) max_eval_err else NA_real_,
        max_root_error = if (overall_certified && isTRUE(enclose_if_chunks)) max_root_err else NA_real_,
        max_direct_star_error_not_root_radius = if (overall_certified && isTRUE(enclose_if_chunks)) max_direct_star_err else NA_real_
      ),
      status_separation = list(
        point_candidate_fitted = TRUE,
        point_candidate_is_certificate = FALSE,
        raw_root_certified = isTRUE(prob$certified),
        link_enclosure_certified = isTRUE(enclosed_link$certified),
        if_chunks_certified = isTRUE(enclose_if_chunks) && isTRUE(if_all_certified),
        covariance_transfer_certified = FALSE,
        stored_array_gmm_accepted = FALSE,
        complete_statistical_inference_certified = FALSE,
        public_numerical_policy_adopted = FALSE
      )
    ),
    class = "cdid_ipw_point_link_certificate_binding"
  )
}
