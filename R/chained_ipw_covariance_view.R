# ===========================================================================
# Chained IPW Compact-View Covariance Transfer Adapter
# ===========================================================================
#
# Private common-A compact-view adapter for IPW covariance transfer.
# Composes the immutable c45 streaming Gram and covariance-transfer kernel
# (`ipw_cov_transfer_stream_kernel`) with the immutable compact owner/view
# layer (`cdid_ipw_panel_owner`, `cdid_ipw_frame_owner`, `cdid_ipw_nuisance_view`,
# `cdid_ipw_link_view`) and the immutable compact IF chunk enclosure layer
# (`preflight_ipw_link_if_view_chunk`, `enclose_ipw_link_if_view_chunk`,
# `enclose_ipw_link_if_view_chunk_with_scratch`, `new_ipw_if_chunk_scratch`,
# `ipw_validate_if_chunk_scratch`).
#
# Strictly unexported and internal. All arguments of the entrypoint
# `certify_ipw_covariance_transfer_view` are explicit with no defaults.
# ===========================================================================

ipw_cov_transfer_view_require_deps <- function() {
  ipw_cov_transfer_require_deps()
  req_fns <- c(
    "validate_ipw_owner_control",
    "validate_ipw_panel_owner",
    "validate_ipw_frame_owner",
    "validate_ipw_nuisance_view",
    "validate_ipw_link_view",
    "ipw_owner_id_byte_metadata",
    "ipw_owner_safe_product",
    "ipw_owner_safe_sum",
    "ipw_validate_nuisance_problem_mode",
    "ipw_validate_link_enclosure_mode",
    "preflight_ipw_link_if_view_chunk",
    "enclose_ipw_link_if_view_chunk",
    "enclose_ipw_link_if_view_chunk_with_scratch",
    "new_ipw_if_chunk_scratch",
    "ipw_validate_if_chunk_scratch",
    "ipw_cov_transfer_stream_kernel",
    "validate_ipw_covariance_transfer_control",
    "ipw_cov_transfer_validate_operator_header",
    "ipw_cov_transfer_validate_operator_content"
  )
  for (fn in req_fns) {
    if (!exists(fn, mode = "function")) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Required compact-view covariance transfer dependency `", fn, "` is missing."),
        "cdid_ipw_covariance_transfer_dependency_error"
      )
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Compact Frame Owner and Link View Header/Content Validators
# ---------------------------------------------------------------------------

# Length/byte metadata only: nchar(type="bytes") reads stored CHARSXP byte
# lengths, not string contents. Reserve vector slots and ALL comparison work
# before any nzchar/unique/identical/match content validation. No key coercion.
ipw_cov_transfer_view_numeric_id_key <- function(id) {
  out <- as.double(id)
  out[out == 0.0] <- 0.0
  out
}

ipw_cov_transfer_view_key_budget <- function(vec, max_entries, max_work,
                                             base_entries, work_so_far,
                                             bytes_so_far, multiplicity) {
  if (is.null(vec) || !is.character(vec)) {
    return(list(bytes = bytes_so_far, work = work_so_far))
  }
  count <- length(vec)
  ipw_cov_transfer_safe_sum(c(base_entries, 4 * count, ceiling(2 * bytes_so_far / 8)),
    max_entries, "key metadata vector slots before byte inspection")
  work <- ipw_cov_transfer_safe_sum(c(work_so_far, 8 * count), max_work,
    "key metadata length work before byte inspection")
  bytes <- bytes_so_far
  for (i in seq_len(count)) {
    nb <- nchar(vec[[i]], type = "bytes", allowNA = TRUE, keepNA = TRUE)
    if (is.na(nb)) nb <- 0 # NA is refused by the existing content validators.
    bytes <- ipw_cov_transfer_safe_sum(c(bytes, nb), label = "key byte metadata")
    # Reserve before content comparisons, including adversarial common prefixes.
    work <- ipw_cov_transfer_safe_sum(c(work, multiplicity * nb), max_work,
      "key byte comparison work before content validation")
    ipw_cov_transfer_safe_sum(c(base_entries, 4 * count, ceiling(2 * bytes / 8)),
      max_entries, "key byte scratch before content validation")
  }
  list(bytes = bytes, work = work)
}

# Pre-header phase. The immutable validators may compare labels and seal fields;
# bound the keys on both supplied and sealed objects before invoking them.
ipw_cov_transfer_view_header_budget <- function(links, frame_owner, candidate_rows,
                                                C, op, V_candidate, ctrl) {
  L <- length(links)
  entries <- ipw_cov_transfer_safe_product(66, L, ctrl$max_matrix_entries,
    "link header vectors before header traversal")
  ref_count <- ipw_cov_transfer_safe_product(L, L - 1, label = "two distinct-reference scan counts")
  ref_work <- ipw_cov_transfer_safe_product(8, ref_count, ctrl$max_work_ops,
    "distinct-reference comparison work before header traversal")
  header_work <- ipw_cov_transfer_safe_sum(c(ref_work, 66 * L), ctrl$max_work_ops,
    "header plus distinct-reference work before traversal")
  budget <- list(bytes = 0, work = header_work)
  max_key_length <- 0
  charge <- function(x) {
    if (is.character(x)) max_key_length <<- max(max_key_length, length(x))
    budget <<- ipw_cov_transfer_view_key_budget(x, ctrl$max_matrix_entries,
      ctrl$max_work_ops, entries, budget$work, budget$bytes, max(64, 8 * (L + 1)))
  }
  obj_keys <- function(x) {
    if (is.list(x)) {
      charge(names(x))
      charge(x[["link_id"]])
    }
  }
  snapshot <- function(x, field, binding) {
    ref <- if (is.list(x)) x[[field]] else NULL
    if (is.environment(ref) && environmentIsLocked(ref) &&
        exists(binding, envir = ref, inherits = FALSE) && bindingIsLocked(binding, ref)) {
      get(binding, envir = ref, inherits = FALSE)
    } else NULL
  }
  charge(names(links))
  obj_keys(frame_owner)
  f <- snapshot(frame_owner, "ref", "frame_owner")
  obj_keys(f)
  if (is.list(f)) { charge(f[["a_id"]]); charge(f[["a_id_key"]]) }
  for (l in seq_len(L)) {
    arg <- links[[l]]; obj_keys(arg)
    lk <- snapshot(arg, "ref", "link_enclosure"); obj_keys(lk)
    pr <- snapshot(lk, "problem_ref", "nuisance_problem"); obj_keys(pr)
    nv <- snapshot(pr, "nuisance_ref", "nuisance_view"); obj_keys(nv)
    lv <- snapshot(lk, "link_ref", "link_view"); obj_keys(lv)
    # A bare callback is normalized to a candidate spec. Reserve its implicit
    # declared IDs too, so the same callback has a stable preflight envelope
    # whether supplied as a function or with the canonical explicit spec.
    if (is.function(candidate_rows) && is.list(lk)) charge(lk[["link_id"]])
  }
  if (is.list(candidate_rows)) { obj_keys(candidate_rows); charge(candidate_rows[["link_ids"]]) }
  if (is.function(candidate_rows)) charge(c("callback", "link_ids"))
  charge(colnames(C)); charge(rownames(C))
  charge(colnames(V_candidate)); charge(rownames(V_candidate))
  if (is.list(op)) {
    charge(names(op)); charge(op[["purpose"]]); charge(op[["row_level"]])
    charge(op[["cluster_id"]]); charge(op[["cluster_order"]])
  }
  list(reference_comparisons_bound = ref_count, reference_work_bound = ref_work,
       key_work_bound = budget$work - header_work, key_bytes = budget$bytes,
       entries_bound = entries + 4 * max_key_length + ceiling(2 * budget$bytes / 8))
}

ipw_cov_transfer_validate_frame_owner_header <- function(frame_owner) {
  if (missing(frame_owner) || is.null(frame_owner)) {
    cdid_ipw_covariance_transfer_abort(
      "Explicit sealed `frame_owner` (`cdid_ipw_frame_owner`) is required.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(
    frame_owner,
    "frame_owner",
    allowed_attrs = c("names", "class")
  )
  if (!inherits(frame_owner, "cdid_ipw_frame_owner") || !is.list(frame_owner)) {
    cdid_ipw_covariance_transfer_abort(
      "`frame_owner` must be a sealed `cdid_ipw_frame_owner` object; dense frame specs or unsealed lists are refused.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  validate_ipw_frame_owner(frame_owner)

  f_ref <- frame_owner$ref
  p_ref <- frame_owner$panel_ref
  if (!is.environment(f_ref) || !environmentIsLocked(f_ref) ||
      !exists("frame_owner", envir = f_ref, inherits = FALSE) ||
      !bindingIsLocked("frame_owner", f_ref) ||
      !is.environment(p_ref) || !environmentIsLocked(p_ref) ||
      !exists("panel_owner", envir = p_ref, inherits = FALSE) ||
      !bindingIsLocked("panel_owner", p_ref)) {
    cdid_ipw_covariance_transfer_abort(
      "`frame_owner` lacks locked `ref` or `panel_ref` environments.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_integrity_error")
    )
  }

  sealed_frame <- get("frame_owner", envir = f_ref, inherits = FALSE)
  sealed_panel <- get("panel_owner", envir = p_ref, inherits = FALSE)
  validate_ipw_frame_owner(sealed_frame)
  validate_ipw_panel_owner(sealed_panel)

  if (!identical(sealed_frame$ref, f_ref) ||
      !identical(sealed_frame$panel_ref, p_ref) ||
      !identical(sealed_panel$ref, p_ref)) {
    cdid_ipw_covariance_transfer_abort(
      "`frame_owner` and `panel_owner` reference seals are inconsistent.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_integrity_error")
    )
  }

  n_a <- ipw_cov_transfer_int_bound(sealed_frame$n_a_units, "frame_owner$n_a_units", min_val = 1)
  p_cov <- ipw_cov_transfer_int_bound(sealed_frame$n_covariates, "frame_owner$n_covariates", min_val = 1)
  if (!identical(as.integer(sealed_panel$n_a_units), as.integer(n_a)) ||
      length(sealed_frame$a_id) != n_a ||
      length(sealed_frame$a_id_key) != n_a ||
      length(sealed_frame$w_a) != n_a ||
      length(sealed_frame$a_g) != n_a ||
      !is.matrix(sealed_frame$Z_a) ||
      nrow(sealed_frame$Z_a) != n_a ||
      ncol(sealed_frame$Z_a) != p_cov ||
      !identical(sealed_frame$id_type, sealed_panel$id_type) ||
      !(sealed_frame$id_type %in% c("numeric", "character"))) {
    cdid_ipw_covariance_transfer_abort(
      "`frame_owner` internal canonical A-unit dimensions or ID types disagree with its sealed `panel_owner`.",
      c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_integrity_error")
    )
  }

  list(
    frame_owner = sealed_frame,
    panel_owner = sealed_panel,
    frame_ref = f_ref,
    panel_ref = p_ref,
    n = as.integer(n_a),
    normalizer_n = as.double(n_a),
    p = as.double(p_cov),
    id_type = sealed_frame$id_type,
    shared_owner_retained_entries = as.double(sealed_panel$preflight$total_retained_entries) +
      as.double(sealed_frame$preflight$total_retained_entries)
  )
}

ipw_cov_transfer_validate_links_view_headers <- function(links, frame_hdr) {
  if (missing(links) || is.null(links) || !is.list(links) || is.object(links) || length(links) < 1L) {
    cdid_ipw_covariance_transfer_abort(
      "`links` must be a non-empty plain list of sealed compact-view `cdid_ipw_link_enclosure` objects.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  ipw_cov_transfer_refuse_decorative_attrs(links, "links", allowed_attrs = "names")
  L <- length(links)
  list_names <- names(links)
  if (!is.null(list_names) && (anyNA(list_names) || any(!nzchar(list_names)) || anyDuplicated(list_names))) {
    cdid_ipw_covariance_transfer_abort(
      "If `links` is a named list, its names must be unique non-empty character link IDs.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }

  f_ref <- frame_hdr$frame_ref
  p_ref <- frame_hdr$panel_ref
  f_owner <- frame_hdr$frame_owner
  n_ref <- frame_hdr$n
  norm_ref <- frame_hdr$normalizer_n
  id_type_ref <- frame_hdr$id_type
  p_ref_dim <- frame_hdr$p

  links_sealed <- vector("list", L)
  problems <- vector("list", L)
  nuisance_views <- vector("list", L)
  link_views <- vector("list", L)
  p_vec <- numeric(L)
  m_vec <- numeric(L)
  cert_ctrls <- vector("list", L)
  link_ids <- character(L)

  seen_prob_refs <- list()
  seen_nv_refs <- list()
  first_problem <- logical(L)
  reference_comparisons <- 0
  existing_retained <- 0.0
  shared_view_retained <- 0.0

  for (l in seq_len(L)) {
    lk_arg <- links[[l]]
    ipw_cov_transfer_refuse_decorative_attrs(
      lk_arg,
      paste0("links[[", l, "]]"),
      allowed_attrs = c("names", "class")
    )
    prob_from_mode <- ipw_validate_link_enclosure_mode(
      lk_arg,
      require_certified = FALSE,
      expected_storage_mode = "compact_view"
    )
    lk <- get("link_enclosure", envir = lk_arg$ref, inherits = FALSE)
    prob_l <- get("nuisance_problem", envir = lk$problem_ref, inherits = FALSE)
    if (!identical(prob_l$ref, prob_from_mode$ref)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` locked problem reference disagrees with validated problem snapshot."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
      )
    }
    ipw_validate_nuisance_problem_mode(
      prob_l,
      require_certified = FALSE,
      expected_storage_mode = "compact_view"
    )

    # Refuse any forged dense fields on compact problem or link enclosure (use exact [[ indexing, never $ partial matching)
    ep_l <- lk[["endpoints"]]
    if (!is.null(prob_l$raw[["Z"]]) || !is.null(prob_l$raw[["id"]]) || !is.null(prob_l$raw[["id_key"]]) ||
        !is.null(prob_l$raw[["w"]]) || !is.null(prob_l$raw[["G"]]) || !is.null(prob_l$raw[["pair_observed"]]) ||
        !is.null(prob_l$masks[["in_risk"]]) || !is.null(prob_l$masks[["in_observed"]]) ||
        !is.null(prob_l$masks[["cohort_index"]]) || !is.null(prob_l$masks[["risk_indices"]]) ||
        !is.null(prob_l$masks[["observed_indices"]]) || !is.list(ep_l) ||
        !is.null(ep_l[["obs_a_index"]]) || !is.null(ep_l[["a_to_obs_pos"]]) ||
        !is.null(ep_l[["y_left_obs"]]) || !is.null(ep_l[["y_right_obs"]])) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` or its nuisance problem contains forged dense arrays in `compact_view` mode."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
      )
    }

    # Require actual identical sealed frame_ref and panel_ref across frame_owner, link, and problem.
    # Equal unit counts, labels, or arrays on a distinct owner instance never substitute for owner identity.
    if (!identical(prob_l$frame_ref, f_ref) || !identical(lk$frame_ref, f_ref) ||
        !identical(prob_l$panel_ref, p_ref) || !identical(lk$panel_ref, p_ref)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` or its nuisance problem is not bound to the exact supplied `frame_owner`/`panel_owner` reference."),
        c(
          "cdid_ipw_covariance_transfer_frame_error",
          "cdid_ipw_covariance_transfer_link_error",
          "cdid_ipw_owner_mismatch_error",
          "cdid_ipw_context_mismatch_error"
        )
      )
    }

    if (!is.environment(prob_l$nuisance_ref) || !environmentIsLocked(prob_l$nuisance_ref) ||
        !exists("nuisance_view", envir = prob_l$nuisance_ref, inherits = FALSE) ||
        !is.environment(lk$link_ref) || !environmentIsLocked(lk$link_ref) ||
        !exists("link_view", envir = lk$link_ref, inherits = FALSE)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` lacks locked `nuisance_ref` or `link_ref` environments."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
      )
    }

    nv_l <- get("nuisance_view", envir = prob_l$nuisance_ref, inherits = FALSE)
    lv_l <- get("link_view", envir = lk$link_ref, inherits = FALSE)
    validate_ipw_nuisance_view(nv_l, frame_owner = f_owner)
    validate_ipw_link_view(lv_l, frame_owner = f_owner, nuisance_view = nv_l)

    if (!identical(nv_l$frame_ref, f_ref) || !identical(nv_l$panel_ref, p_ref) ||
        !identical(lv_l$frame_ref, f_ref) || !identical(lv_l$panel_ref, p_ref) ||
        !identical(lk$nuisance_ref, nv_l$ref) || !identical(lv_l$nuisance_ref, nv_l$ref) ||
        !isTRUE(lv_l$nuisance_mask_verified)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` bound `nuisance_view` or `link_view` fails owner/mask reference identity."),
        c(
          "cdid_ipw_covariance_transfer_frame_error",
          "cdid_ipw_covariance_transfer_link_error",
          "cdid_ipw_owner_mismatch_error"
        )
      )
    }

    # Require problem$provenance$normalizer_n == frame_owner$n_a_units == problem$raw$n_a_units
    raw_l <- prob_l$raw
    if (is.null(prob_l$provenance) ||
        !identical(as.double(prob_l$provenance$normalizer_n), norm_ref) ||
        !identical(as.double(raw_l$n_a_units), norm_ref) ||
        !identical(as.integer(raw_l$n_a_units), n_ref) ||
        !identical(raw_l$id_type, id_type_ref) ||
        (!is.null(lk$provenance$normalizer_n) && !identical(as.double(lk$provenance$normalizer_n), norm_ref))) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` normalizer, `raw$n_a_units`, or `raw$id_type` disagrees with `frame_owner`."),
        c("cdid_ipw_covariance_transfer_frame_error", "cdid_ipw_covariance_transfer_link_error")
      )
    }

    p_l <- ipw_cov_transfer_int_bound(raw_l$p_covariates, paste0("links[[", l, "]] raw$p_covariates"), min_val = 1)
    m_l <- ipw_cov_transfer_int_bound(length(raw_l$cohorts), paste0("links[[", l, "]] length(raw$cohorts)"), min_val = 2)
    if (p_l != p_ref_dim || m_l != as.double(nv_l$cohorts$m_cohorts) || m_l != as.double(lv_l$cohorts$m_cohorts)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` covariate count `p` or cohort count `m` disagrees with `frame_owner` or bound views."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
      )
    }

    lid_lk <- lk$link_id
    lid_lv <- lv_l$link_id
    if (!is.null(lid_lk) && !is.null(lid_lv) && !identical(lid_lk, lid_lv)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]$link_id` ('", lid_lk, "') disagrees with bound `link_view$link_id` ('", lid_lv, "')."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    lid <- if (!is.null(lid_lk)) lid_lk else lid_lv
    if (!is.null(lid)) {
      if (!gmm_scalar_character(lid)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`links[[", l, "]]$link_id` must be a non-empty character scalar."),
          c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
      if (!is.null(list_names) && !identical(list_names[l], lid)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`names(links)[", l, "]` ('", list_names[l], "') disagrees with `links[[", l, "]]$link_id` ('", lid, "')."),
          c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
        )
      }
      link_ids[l] <- lid
    } else if (!is.null(list_names)) {
      link_ids[l] <- list_names[l]
    } else {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` must declare an explicit `link_id` (or `names(links)` must supply explicit link IDs)."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }

    links_sealed[[l]] <- lk
    problems[[l]] <- prob_l
    nuisance_views[[l]] <- nv_l
    link_views[[l]] <- lv_l
    p_vec[l] <- p_l
    m_vec[l] <- m_l
    cert_ctrls[[l]] <- validate_logit_certificate_control(prob_l$controls$certificate_control)

    # Deduplicate pre-existing nuisance problem state by actual sealed problem$ref
    reference_comparisons <- reference_comparisons + length(seen_prob_refs)
    already_seen_prob <- any(vapply(seen_prob_refs, identical, logical(1L), y = prob_l$ref))
    first_problem[l] <- !already_seen_prob
    if (!already_seen_prob) {
      seen_prob_refs[[length(seen_prob_refs) + 1L]] <- prob_l$ref
      existing_retained <- existing_retained + as.double(prob_l$preflight$total_retained_entries)
    }
    existing_retained <- existing_retained + as.double(lk$preflight$total_retained_entries)

    reference_comparisons <- reference_comparisons + length(seen_nv_refs)
    already_seen_nv <- any(vapply(seen_nv_refs, identical, logical(1L), y = nv_l$ref))
    if (!already_seen_nv) {
      seen_nv_refs[[length(seen_nv_refs) + 1L]] <- nv_l$ref
      shared_view_retained <- shared_view_retained + as.double(nv_l$preflight$total_retained_entries)
    }
    shared_view_retained <- shared_view_retained + as.double(lv_l$preflight$total_retained_entries)
  }

  if (anyDuplicated(link_ids)) {
    cdid_ipw_covariance_transfer_abort(
      "Declared `link_id` values across `links` must be unique.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error"),
      diagnostics = list(link_ids = link_ids)
    )
  }

  list(
    L = L,
    n = n_ref,
    normalizer_n = norm_ref,
    id_type = id_type_ref,
    link_ids = link_ids,
    links_sealed = links_sealed,
    problems = problems,
    nuisance_views = nuisance_views,
    link_views = link_views,
    first_problem = first_problem,
    reference_comparisons = reference_comparisons,
    n_unique_problems = length(seen_prob_refs),
    n_unique_nuisance_views = length(seen_nv_refs),
    p_vec = p_vec,
    m_vec = m_vec,
    certificate_controls = cert_ctrls,
    existing_4bd_retained_entries = existing_retained,
    shared_view_retained_entries = shared_view_retained
  )
}

ipw_cov_transfer_validate_links_view_content <- function(link_hdr, frame_hdr) {
  L <- link_hdr$L
  n <- frame_hdr$n
  norm_n <- frame_hdr$normalizer_n

  for (l in seq_len(L)) {
    ipw_cov_transfer_check_interrupt()
    lk <- link_hdr$links_sealed[[l]]
    prob_l <- link_hdr$problems[[l]]
    nv_l <- link_hdr$nuisance_views[[l]]
    lv_l <- link_hdr$link_views[[l]]
    raw_l <- prob_l$raw
    masks_l <- prob_l$masks

    dims_l <- prob_l$dimensions
    if (link_hdr$first_problem[l]) {
      if (!identical(as.double(raw_l$cohorts), as.double(nv_l$cohorts$cohorts)) ||
          !identical(as.integer(dims_l$n_risk_units), as.integer(nv_l$support$n_risk_a_units)) ||
          !identical(as.integer(dims_l$n_observed_risk_units), as.integer(nv_l$support$n_observed_risk_units)) ||
          !identical(as.integer(masks_l$by_cohort_risk_count), as.integer(nv_l$support$by_cohort_a_count)) ||
          !identical(as.integer(masks_l$by_cohort_observed_count), as.integer(nv_l$support$by_cohort_pair_count)) ||
          !identical(as.double(prob_l$provenance$normalizer_n), norm_n)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`links[[", l, "]]` compact nuisance problem cohorts or support counts disagree with bound `nuisance_view`."),
          c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
        )
      }
    }

    if (!identical(as.double(lv_l$cohorts$cohorts), as.double(raw_l$cohorts)) ||
        !identical(as.integer(lv_l$support$n_a_units), as.integer(n)) ||
        !identical(as.integer(lv_l$support$n_risk_a_units), as.integer(dims_l$n_risk_units)) ||
        !identical(as.integer(lv_l$support$n_observed_risk_units), as.integer(dims_l$n_observed_risk_units)) ||
        !identical(as.double(lk$provenance$normalizer_n), norm_n)) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`links[[", l, "]]` compact link enclosure cohorts or support counts disagree with bound `link_view`."),
        c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_integrity_error")
      )
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Compact-View Preflight Composition (executed BEFORE content scans/callbacks)
# ---------------------------------------------------------------------------

preflight_ipw_covariance_transfer_view <- function(
  n,
  L,
  d,
  row_level,
  n_clusters,
  p_vec,
  m_vec,
  certificate_controls,
  control,
  reader_control,
  id_is_character = FALSE,
  cluster_is_character = FALSE,
  existing_4bd_retained_entries = 0,
  header_key_work_ops = 0
) {
  ipw_cov_transfer_view_require_deps()
  if (missing(n) || missing(L) || missing(d) || missing(row_level) ||
      missing(n_clusters) || missing(p_vec) || missing(m_vec) ||
      missing(certificate_controls) || missing(control) || missing(reader_control)) {
    cdid_ipw_covariance_transfer_abort(
      "All required arguments of `preflight_ipw_covariance_transfer_view` must be supplied explicitly.",
      "cdid_ipw_covariance_transfer_input_error"
    )
  }
  ctrl <- validate_ipw_covariance_transfer_control(control)
  ipw_cov_transfer_refuse_decorative_attrs(
    reader_control,
    "reader_control",
    allowed_attrs = c("names", "class", "expected_a_id", "trace_env")
  )
  r_ctrl <- validate_ipw_owner_control(reader_control)

  n_d <- ipw_cov_transfer_int_bound(n, "n_a_units", min_val = 1)
  L_d <- ipw_cov_transfer_int_bound(L, "n_links", min_val = 1)
  d_d <- ipw_cov_transfer_int_bound(d, "n_moments", min_val = 1)
  if (!gmm_scalar_character(row_level) || !row_level %in% c("unit", "cluster")) {
    cdid_ipw_covariance_transfer_abort(
      "`row_level` in preflight must be `'unit'` or `'cluster'`.",
      "cdid_ipw_covariance_transfer_operator_error"
    )
  }
  N_c <- ipw_cov_transfer_int_bound(n_clusters, "n_clusters", min_val = 1)
  if (identical(row_level, "unit") && N_c != n_d) {
    cdid_ipw_covariance_transfer_abort(
      "Unit operator requires `n_clusters == n_a_units`.",
      "cdid_ipw_covariance_transfer_operator_error"
    )
  }
  if (identical(row_level, "cluster") && N_c > ctrl$max_cluster_rows) {
    cdid_ipw_covariance_transfer_abort(
      paste0("Declared cluster count (", N_c, ") exceeds `control$max_cluster_rows` (", ctrl$max_cluster_rows, ")."),
      "cdid_ipw_covariance_transfer_bound_error",
      diagnostics = list(requested_clusters = N_c, max_cluster_rows = ctrl$max_cluster_rows)
    )
  }

  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  lk_max_wo <- ctrl$link_control$max_work_ops
  r_max_wo <- r_ctrl$max_work_ops

  # Effective input chunk size B is bounded by n and every covariance/link/reader chunk cap.
  b_in <- min(
    n_d,
    ctrl$max_chunk_rows,
    ctrl$link_control$max_chunk_rows,
    r_ctrl$max_chunk_rows
  )
  b_op <- if (identical(row_level, "unit")) {
    min(b_in, ctrl$max_operator_block_rows)
  } else {
    min(N_c, ctrl$max_operator_block_rows)
  }
  n_input_chunks <- ceiling(n_d / b_in)
  n_full_chunks <- floor(n_d / b_in)
  rem_chunk_rows <- n_d - n_full_chunks * b_in
  n_operator_blocks <- if (identical(row_level, "unit")) {
    n_input_chunks * ceiling(b_in / b_op)
  } else {
    ceiling(N_c / b_op)
  }

  # 1. Per-link compact view preflights and actual-range chunk work composition
  if (length(p_vec) != L_d || length(m_vec) != L_d || length(certificate_controls) != L_d) {
    cdid_ipw_covariance_transfer_abort(
      "Link dimension/control vectors in preflight must have length `L`.",
      "cdid_ipw_covariance_transfer_input_error"
    )
  }

  max_view_combined_concurrent <- 0
  max_view_arith_cached_concurrent <- 0
  max_view_reader_concurrent <- 0
  max_view_cache_payload <- 0
  max_4bd_output_entries <- 0

  total_4bd_stream_work_ops <- 0
  total_cache_copy_work_ops <- 0
  total_reader_base_work_ops <- 0
  total_compact_view_stream_work_ops <- 0

  link_chunk_preflights <- vector("list", as.integer(L_d))

  for (l in seq_len(as.integer(L_d))) {
    pre_l_full <- preflight_ipw_link_if_view_chunk(
      b_len = b_in,
      p = p_vec[l],
      m = m_vec[l],
      control = ctrl$link_control,
      reader_control = r_ctrl,
      certificate_control = certificate_controls[[l]],
      include_components = FALSE,
      include_shared_xi = FALSE,
      has_phi_candidate = TRUE,
      has_xi_beta = FALSE,
      read_outcomes = TRUE,
      panel_owner = NULL,
      a_seq = NULL,
      expected_a_id = NULL,
      caller_phi_entries = b_in,
      caller_xi_entries = 0
    )
    pre_l_rem <- if (rem_chunk_rows > 0) {
      preflight_ipw_link_if_view_chunk(
        b_len = rem_chunk_rows,
        p = p_vec[l],
        m = m_vec[l],
        control = ctrl$link_control,
        reader_control = r_ctrl,
        certificate_control = certificate_controls[[l]],
        include_components = FALSE,
        include_shared_xi = FALSE,
        has_phi_candidate = TRUE,
        has_xi_beta = FALSE,
        read_outcomes = TRUE,
        panel_owner = NULL,
        a_seq = NULL,
        expected_a_id = NULL,
        caller_phi_entries = rem_chunk_rows,
        caller_xi_entries = 0
      )
    } else NULL

    # P_view,l is the compact preflight's COMBINED arithmetic/reader/cache payload
    if (pre_l_full$combined_concurrent_entries > max_view_combined_concurrent) {
      max_view_combined_concurrent <- pre_l_full$combined_concurrent_entries
    }
    if (pre_l_full$concurrent_entries > max_view_arith_cached_concurrent) {
      max_view_arith_cached_concurrent <- pre_l_full$concurrent_entries
    }
    if (pre_l_full$reader_concurrent_entries > max_view_reader_concurrent) {
      max_view_reader_concurrent <- pre_l_full$reader_concurrent_entries
    }
    if (pre_l_full$cache_payload_entries > max_view_cache_payload) {
      max_view_cache_payload <- pre_l_full$cache_payload_entries
    }
    if (pre_l_full$output_entries > max_4bd_output_entries) {
      max_4bd_output_entries <- pre_l_full$output_entries
    }

    # Exact actual-range decomposition across the n_full_chunks of size b_in and optional tail chunk of size rem_chunk_rows
    full_interval_ops <- pre_l_full$work_ops_bound - pre_l_full$cache_copy_work_ops_bound
    rem_interval_ops <- if (!is.null(pre_l_rem)) {
      pre_l_rem$work_ops_bound - pre_l_rem$cache_copy_work_ops_bound
    } else 0

    link_interval_ops <- ipw_cov_transfer_safe_sum(
      c(
        ipw_cov_transfer_safe_product(n_full_chunks, full_interval_ops, lk_max_wo, paste0("full-chunk interval work for link ", l)),
        rem_interval_ops
      ),
      max_bound = lk_max_wo,
      label = paste0("whole-stream interval-only IF work for link ", l)
    )
    link_cache_copy_ops <- ipw_cov_transfer_safe_sum(
      c(
        ipw_cov_transfer_safe_product(n_full_chunks, pre_l_full$cache_copy_work_ops_bound, lk_max_wo, paste0("full-chunk cache-copy work for link ", l)),
        if (!is.null(pre_l_rem)) pre_l_rem$cache_copy_work_ops_bound else 0
      ),
      max_bound = lk_max_wo,
      label = paste0("whole-stream quadratic 4B^2 cache-copy work for link ", l)
    )
    link_arith_cached_ops <- ipw_cov_transfer_safe_sum(
      c(link_interval_ops, link_cache_copy_ops),
      max_bound = lk_max_wo,
      label = paste0("whole-stream interval plus cache-copy work for link ", l)
    )
    link_reader_base_ops <- ipw_owner_safe_sum(
      c(
        ipw_owner_safe_product(n_full_chunks, pre_l_full$reader_base_work_ops_bound, r_max_wo, paste0("full-chunk reader work for link ", l)),
        if (!is.null(pre_l_rem)) pre_l_rem$reader_base_work_ops_bound else 0
      ),
      max_bound = r_max_wo,
      label = paste0("whole-stream base reader work for link ", l)
    )
    link_combined_base_ops <- ipw_cov_transfer_safe_sum(
      c(link_arith_cached_ops, link_reader_base_ops),
      max_bound = lk_max_wo,
      label = paste0("whole-stream combined compact IF/cache/reader work for link ", l)
    )

    pre_l_full$whole_stream_interval_work_ops_bound <- link_interval_ops
    pre_l_full$whole_stream_cache_copy_work_ops_bound <- link_cache_copy_ops
    pre_l_full$whole_stream_reader_base_work_ops_bound <- link_reader_base_ops
    pre_l_full$whole_stream_combined_work_ops_bound <- link_combined_base_ops
    link_chunk_preflights[[l]] <- pre_l_full

    total_4bd_stream_work_ops <- ipw_cov_transfer_safe_sum(
      c(total_4bd_stream_work_ops, link_interval_ops),
      max_bound = max_wo,
      label = "cumulative interval-only 4bd IF chunk work across all links"
    )
    total_cache_copy_work_ops <- ipw_cov_transfer_safe_sum(
      c(total_cache_copy_work_ops, link_cache_copy_ops),
      max_bound = max_wo,
      label = "cumulative quadratic 4B^2 cache-copy work across all links"
    )
    total_reader_base_work_ops <- ipw_owner_safe_sum(
      c(total_reader_base_work_ops, link_reader_base_ops),
      max_bound = r_max_wo,
      label = "cumulative base reader work across all links"
    )
    ipw_cov_transfer_safe_sum(
      c(total_reader_base_work_ops, 0),
      max_bound = max_wo,
      label = "cumulative base reader work against covariance max_work_ops"
    )
    total_compact_view_stream_work_ops <- ipw_cov_transfer_safe_sum(
      c(total_compact_view_stream_work_ops, link_combined_base_ops),
      max_bound = max_wo,
      label = "cumulative combined compact IF/cache/reader stream work across all links"
    )
  }

  # 2. Covariance component concurrent payload entries P_cov (c45 lines 715-777)
  bL <- ipw_cov_transfer_safe_product(b_in, L_d, max_me, "b x L candidate/error chunk block")
  b_by_L_entries <- ipw_cov_transfer_safe_sum(
    c(ipw_cov_transfer_safe_product(4, bL, max_me, "three live buffers plus retained callback backing for integer phi"),
      ipw_cov_transfer_safe_product(12, b_in, max_me, "candidate column/coercion, ID/range and validation temporaries")),
    max_bound = max_me,
    label = "bounded b x L IF and separated error buffers"
  )

  bd_in <- ipw_cov_transfer_safe_product(b_in, d_d, max_me, "b x d input transform block")
  bd_op <- ipw_cov_transfer_safe_product(b_op, d_d, max_me, "b_op x d operator Gram block")
  b_by_d_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(14, bd_in, max_me, "b x d Uhat, [U], residuals, and mapped error blocks"),
      ipw_cov_transfer_safe_product(12, bd_op, max_me, "five operator input slices, Yhat, |Yhat|, E_mat, and native transpose temporaries")
    ),
    max_bound = max_me,
    label = "bounded b x d transform and operator buffers"
  )

  dL <- ipw_cov_transfer_safe_product(d_d, L_d, max_me, "d x L moment map C")
  d_by_L_entries <- ipw_cov_transfer_safe_product(4, dL, max_me, "caller C, copied C, t(C), and |t(C)| matrices")

  d2 <- ipw_cov_transfer_safe_product(d_d, d_d, max_me, "d x d covariance matrix block")
  d_by_d_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(24, d2, max_me, "V_candidate, S, [L,U], dist, point Gram, D1..D3, and spectral workspace"),
      ipw_cov_transfer_safe_product(4, d_d, max_me, "global t accumulator and moment vector scratch"),
      256
    ),
    max_bound = max_me,
    label = "d x d covariance accumulators and spectral scratch"
  )

  id_membership_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, n_d, max_me, "conservative c45 bookkeeping bound for A-frame IDs and unit validation vectors"),
      if (identical(row_level, "cluster")) {
        ipw_cov_transfer_safe_sum(
          c(
            ipw_cov_transfer_safe_product(2, n_d, max_me, "cluster_id key and unit-to-cluster index map"),
            ipw_cov_transfer_safe_product(2, N_c, max_me, "cluster_order key and cluster omega vector")
          ),
          max_bound = max_me,
          label = "cluster membership and weight vectors"
        )
      } else {
        ipw_cov_transfer_safe_product(1, n_d, max_me, "unit omega vector")
      },
      ipw_cov_transfer_safe_product(66, L_d, max_me, "ordered link IDs, link headers, and retained per-link preflight scalar payload")
    ),
    max_bound = max_me,
    label = "conservative c45 bookkeeping bound plus operator membership and weight vectors"
  )

  cluster_buffer_entries <- if (identical(row_level, "cluster")) {
    Ncd <- ipw_cov_transfer_safe_product(N_c, d_d, max_me, "N_c x d cluster accumulator matrix")
    ipw_cov_transfer_safe_product(
      5, Ncd, max_me,
      "explicit N_c x d cluster accumulators (z_hat, z_lo, z_hi, ae_eval, ae_root)"
    )
  } else {
    0
  }

  p_cov_entries <- ipw_cov_transfer_safe_sum(
    c(
      b_by_L_entries,
      b_by_d_entries,
      d_by_L_entries,
      d_by_d_entries,
      id_membership_entries,
      cluster_buffer_entries
    ),
    max_bound = max_me,
    label = "covariance-only simultaneous buffers P_cov"
  )

  # Base simultaneous entries before extra key-byte temporaries: P_cov + max_l P_view,l
  # where P_view,l is the compact preflight's combined arithmetic/reader/cache payload.
  base_component_concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(p_cov_entries, max_view_combined_concurrent),
    max_bound = max_me,
    label = "simultaneous P_cov plus max_l P_view,l combined compact chunk buffers"
  )

  min_char_temp_bytes <- ipw_cov_transfer_safe_sum(
    c(
      if (isTRUE(id_is_character)) 2 * (n_d + b_in) else 0,
      if (isTRUE(cluster_is_character)) 2 * (n_d + N_c) else 0,
      2 * L_d
    ),
    label = "minimum character comparison temporary bytes"
  )
  min_char_temp_entries <- ceiling(min_char_temp_bytes / 8)

  min_reader_id_byte_scratch_entries <- if (isTRUE(id_is_character)) {
    ipw_owner_safe_product(
      64, b_in, r_ctrl$max_matrix_entries,
      "minimum reader character expected_a_id temporary slots"
    )
  } else {
    0
  }
  if (min_reader_id_byte_scratch_entries > 0) {
    ipw_owner_safe_sum(
      c(max_view_reader_concurrent, min_reader_id_byte_scratch_entries),
      max_bound = r_ctrl$max_matrix_entries,
      label = "reader concurrent entries plus minimum character expected_a_id scratch"
    )
    ipw_cov_transfer_safe_sum(
      c(max_view_combined_concurrent, min_reader_id_byte_scratch_entries),
      max_bound = ctrl$link_control$max_matrix_entries,
      label = "compact view combined concurrent entries plus minimum character expected_a_id scratch"
    )
  }

  p_extra_keys <- ipw_cov_transfer_safe_sum(
    c(min_char_temp_entries, min_reader_id_byte_scratch_entries),
    max_bound = max_me,
    label = "extra key-byte temporary entries P_extra_keys"
  )

  concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(base_component_concurrent_entries, p_extra_keys),
    max_bound = max_me,
    label = "total simultaneous payload entries P_cov + max_l P_view,l + P_extra_keys"
  )

  # 3. Output payload entry bound
  output_entries <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, d2, max_out, "enclosure [L,U], distance, V_candidate, and both retained S endpoints"),
      dL,
      ipw_cov_transfer_safe_product(4, d_d, max_out, "global t interval and moment summaries"),
      if (identical(row_level, "cluster")) {
        ipw_cov_transfer_safe_product(
          3,
          ipw_cov_transfer_safe_product(N_c, d_d, max_out, "N_c x d cluster output block"),
          max_out,
          "retained cluster point and both interval endpoint arrays"
        )
      } else 0,
      ipw_cov_transfer_safe_product(64, L_d, max_out, "retained per-link preflight scalar payload"),
      256
    ),
    max_bound = max_out,
    label = "covariance transfer output entries"
  )

  # 4. Whole-stream covariance component interval operation bound (c45 lines 826-884)
  per_input_row_transform_ops <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, dL, max_wo, "per-row Uhat residual and |E_eval| |C|', |E_root| |C|' dot ops"),
      ipw_cov_transfer_safe_product(32, d_d, max_wo, "per-row transform residual subtraction, native returned norms, interval addition, and error extraction ops"),
      16,
      if (identical(row_level, "cluster")) ipw_cov_transfer_safe_product(3, d_d, max_wo, "per-row outward cluster accumulation ops") else 0
    ),
    max_bound = max_wo,
    label = "per-input-row transform and cluster accumulation ops"
  )
  input_transform_ops <- ipw_cov_transfer_safe_product(
    n_d, per_input_row_transform_ops, max_wo, "whole-stream input transform and cluster sum ops"
  )

  per_operator_row_ops <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(8, d2, max_wo, "per-operator-row point Gram residual and D1, D2, D3 dot ops"),
      ipw_cov_transfer_safe_product(16, d_d, max_wo, "per-operator-row sqrt(omega), y_iv, E_mat, t, and nu_kind^2 ops"),
      8
    ),
    max_bound = max_wo,
    label = "per-operator-row Gram and error accumulation ops"
  )
  operator_row_ops <- ipw_cov_transfer_safe_product(
    N_c, per_operator_row_ops, max_wo, "whole-stream operator-row Gram ops"
  )
  operator_block_merge_ops <- ipw_cov_transfer_safe_product(
    n_operator_blocks,
    ipw_cov_transfer_safe_sum(
      c(ipw_cov_transfer_safe_product(32, d2, max_wo, "per-operator-block residual subtraction, native returned norms, D sum, widening, and S addition ops"), 32),
      max_bound = max_wo, label = "operator-block native norm constants"
    ),
    max_wo,
    "whole-stream operator-block merge ops"
  )

  d3 <- ipw_cov_transfer_safe_product(d2, d_d, max_wo, "d^3 spectral residual ops")
  post_gram_ops <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(6, d3, max_wo, "spectral certificate Q'Q and S*Q-Q*Lambda residual ops"),
      ipw_cov_transfer_safe_product(64, d2, max_wo, "global centering, kappa scaling, distance, returned norms, and eigen transfer ops"),
      ipw_cov_transfer_safe_product(64, d_d, max_wo, "post-Gram vector ops"),
      256
    ),
    max_bound = max_wo,
    label = "post-Gram centering, norm, and eigen transfer ops"
  )

  component_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(input_transform_ops, operator_row_ops, operator_block_merge_ops, post_gram_ops),
    max_bound = max_wo,
    label = "component covariance transfer interval work bound"
  )

  # Preserve the kernel's total_work_ops_bound as interval-only!
  total_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(total_4bd_stream_work_ops, component_work_ops_bound),
    max_bound = max_wo,
    label = "whole-stream interval-only 4bd plus covariance transfer work bound"
  )

  # 5. Combined outer work bounds: reader ID-byte work + header/operator/ID/key validation work
  min_reader_id_byte_work_ops <- if (isTRUE(id_is_character)) {
    ipw_owner_safe_product(
      256,
      ipw_owner_safe_product(n_d, L_d, r_max_wo, "n x L character expected_a_id slice comparisons"),
      r_max_wo,
      "minimum reader character expected_a_id comparison work across links"
    )
  } else {
    0
  }
  total_reader_work_ops_bound <- ipw_owner_safe_sum(
    c(total_reader_base_work_ops, min_reader_id_byte_work_ops),
    max_bound = r_max_wo,
    label = "total reader work operations bound across all links and chunks"
  )
  ipw_cov_transfer_safe_sum(
    c(total_reader_work_ops_bound, 0),
    max_bound = max_wo,
    label = "total reader work bound against covariance max_work_ops"
  )

  total_compact_view_stream_work_ops <- ipw_cov_transfer_safe_sum(
    c(total_4bd_stream_work_ops, total_cache_copy_work_ops, total_reader_work_ops_bound),
    max_bound = max_wo,
    label = "total compact view stream work bound (interval + 4B^2 cache-copy + reader)"
  )

  base_validation_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(66, L_d, max_wo, "link and owner header validation work"),
      ipw_cov_transfer_safe_product(8, ipw_cov_transfer_safe_product(L_d, L_d - 1),
        max_wo, "two distinct-reference deduplication scans"),
      header_key_work_ops,
      ipw_cov_transfer_safe_product(32, ipw_cov_transfer_safe_product(n_d, L_d),
        max_wo, "candidate finite checks, conversion, link columns and separated-error validation"),
      ipw_cov_transfer_safe_product(16, n_d, max_wo, "exact optional index and ID temporaries"),
      ipw_cov_transfer_safe_product(2048, n_input_chunks, max_wo, "callback field-name byte/metadata validation"),
      8 * (dL + d2),
      ipw_cov_transfer_safe_product(8, n_d, max_wo, "A-frame ID and operator weight validation work"),
      ipw_cov_transfer_safe_product(4, N_c, max_wo, "cluster order and weight validation work"),
      ipw_cov_transfer_safe_product(
        2,
        ipw_cov_transfer_safe_product(n_input_chunks, L_d, max_wo, "chunk x link count"),
        max_wo,
        "chunk callback link_id and range validation work"
      )
    ),
    max_bound = max_wo,
    label = "base header/operator/ID validation work bound"
  )
  min_pass_char_work_bytes <- ipw_cov_transfer_safe_sum(
    c(
      if (isTRUE(id_is_character)) 4 * n_d else 0,
      if (isTRUE(cluster_is_character)) 2 * (n_d + N_c) else 0,
      ipw_cov_transfer_safe_product(2, ipw_cov_transfer_safe_product(n_input_chunks + 1, L_d, max_wo, "pass link_id comparisons"), max_wo, "minimum pass link_id comparison bytes")
    ),
    label = "minimum whole-pass character comparison bytes"
  )
  key_validation_work_ops_bound <- ipw_cov_transfer_safe_product(
    32, min_pass_char_work_bytes, max_wo, "minimum key-byte validation work bound"
  )
  validation_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(base_validation_work_ops_bound, key_validation_work_ops_bound),
    max_bound = max_wo,
    label = "header, operator, ID, and key validation work bound"
  )

  combined_total_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(component_work_ops_bound, total_compact_view_stream_work_ops, validation_work_ops_bound),
    max_bound = max_wo,
    label = "whole-traversal combined covariance, compact IF/reader/cache, and validation work bound"
  )

  total_with_existing <- ipw_cov_transfer_safe_sum(
    c(existing_4bd_retained_entries, concurrent_entries),
    label = "pre-existing compact 4bd retained plus stage concurrent payload"
  )

  list(
    storage_mode = "compact_view",
    n_a_units = n_d,
    n_links = L_d,
    n_moments = d_d,
    row_level = row_level,
    n_clusters = N_c,
    effective_chunk_rows = logit_cert_count_output(b_in),
    effective_operator_block_rows = logit_cert_count_output(b_op),
    n_input_chunks = logit_cert_count_output(n_input_chunks),
    n_full_chunks = logit_cert_count_output(n_full_chunks),
    remainder_chunk_rows = logit_cert_count_output(rem_chunk_rows),
    n_operator_blocks = logit_cert_count_output(n_operator_blocks),
    allocation_phases = list(
      callback_validation = 3 * bL + 12 * b_in,
      candidate_conversion = 3 * bL + 12 * b_in,
      link_traversal = 4 * bL + 12 * b_in + max_view_combined_concurrent,
      kernel_transform_and_gram = p_cov_entries),
    b_by_L_entries = b_by_L_entries,
    b_by_d_entries = b_by_d_entries,
    d_by_L_entries = d_by_L_entries,
    d_by_d_entries = d_by_d_entries,
    id_and_membership_entries = id_membership_entries,
    cluster_buffer_entries = cluster_buffer_entries,
    p_cov_entries = p_cov_entries,
    max_4bd_chunk_concurrent_entries = max_view_combined_concurrent,
    max_view_combined_concurrent_entries = max_view_combined_concurrent,
    max_view_arith_cached_concurrent_entries = max_view_arith_cached_concurrent,
    max_view_reader_concurrent_entries = max_view_reader_concurrent,
    max_view_cache_payload_entries = max_view_cache_payload,
    max_4bd_chunk_output_entries = max_4bd_output_entries,
    base_component_concurrent_entries = base_component_concurrent_entries,
    character_temporary_bytes = min_char_temp_bytes,
    character_temporary_entries_8byte = min_char_temp_entries,
    reader_id_byte_scratch_entries = min_reader_id_byte_scratch_entries,
    extra_key_entries = p_extra_keys,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    output_entries = output_entries,
    total_4bd_stream_work_ops_bound = total_4bd_stream_work_ops,
    component_work_ops_bound = component_work_ops_bound,
    total_work_ops_bound = total_work_ops_bound,
    total_reader_base_work_ops_bound = total_reader_base_work_ops,
    reader_id_byte_work_ops_bound = min_reader_id_byte_work_ops,
    total_reader_work_ops_bound = total_reader_work_ops_bound,
    total_cache_copy_work_ops_bound = total_cache_copy_work_ops,
    total_compact_view_stream_work_ops_bound = total_compact_view_stream_work_ops,
    reference_comparisons_bound = L_d * (L_d - 1),
    reference_comparison_work_ops_bound = 8 * L_d * (L_d - 1),
    header_key_work_ops_bound = header_key_work_ops,
    base_validation_work_ops_bound = base_validation_work_ops_bound,
    key_validation_work_ops_bound = key_validation_work_ops_bound,
    validation_work_ops_bound = validation_work_ops_bound,
    combined_total_work_ops_bound = combined_total_work_ops_bound,
    existing_4bd_retained_entries = as.double(existing_4bd_retained_entries),
    total_with_existing_4bd_entries = total_with_existing,
    total_with_existing_4bd_bytes_8double = 8 * total_with_existing,
    link_chunk_preflights = link_chunk_preflights,
    n_by_L_cache_allocated = FALSE,
    dense_A_incidence_allocated = FALSE,
    per_block_full_n_validation = FALSE,
    zero_n_total_storage_claimed = FALSE,
    peak_rss_promised = FALSE,
    bookkeeping_bound_disclosure = paste0(
      "`id_and_membership_entries` reuses the conservative c45 bookkeeping bound (6*n + operator vectors + 66*L) ",
      "rather than claiming new dense per-link N-vector allocation; operator normalization retains one length-n ",
      "unit `omega` vector or, when clustered, `unit_to_cluster`/`cluster_order` plus five N_c x d cluster accumulators."
    ),
    storage_and_work_disclosure = paste0(
      "Simultaneous stage payload is P_cov + max_l P_view,l + P_extra_keys, where P_view,l is the compact view ",
      "preflight's combined arithmetic, reader, and 80B+2Bp cache payload (shared scratch counted once while live). ",
      "Whole-traversal work sums covariance component interval work, actual-range compact chunk bounds across ",
      "chunks/links (including reader work and quadratic 4B^2 cache-copy work), and header/operator/ID/key validation work. ",
      "Kernel interval-only bounds (`total_4bd_stream_work_ops_bound`, `component_work_ops_bound`, `total_work_ops_bound`) ",
      "and counters are preserved separately from combined reader/cache/key work. Shared owner/view retained state and ",
      "opaque user callback internals are disclosed separately; no process-RSS bound or zero-N total storage is claimed."
    )
  )
}

# ---------------------------------------------------------------------------
# Required Unexported Entry Point:
#   certify_ipw_covariance_transfer_view(
#     links, candidate_rows, frame_owner, C, op, V_candidate, control, reader_control
#   )
# ---------------------------------------------------------------------------

certify_ipw_covariance_transfer_view <- function(
  links,
  candidate_rows,
  frame_owner,
  C,
  op,
  V_candidate,
  control,
  reader_control
) {
  ipw_cov_transfer_view_require_deps()
  if (missing(links) || missing(candidate_rows) || missing(frame_owner) ||
      missing(C) || missing(op) || missing(V_candidate) ||
      missing(control) || missing(reader_control)) {
    cdid_ipw_covariance_transfer_abort(
      "All eight required arguments (`links`, `candidate_rows`, `frame_owner`, `C`, `op`, `V_candidate`, `control`, `reader_control`) must be supplied explicitly.",
      "cdid_ipw_covariance_transfer_input_error"
    )
  }

  # 1. Validate controls and lightweight owner/link headers BEFORE scanning content or calling callbacks
  ctrl <- validate_ipw_covariance_transfer_control(control)
  ipw_cov_transfer_refuse_decorative_attrs(
    reader_control,
    "reader_control",
    allowed_attrs = c("names", "class", "trace_env")
  )
  trace_env <- attr(reader_control, "trace_env", exact = TRUE)
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_ipw_owner_abort(
      "`attr(reader_control, 'trace_env')` must be NULL or an environment.",
      "cdid_ipw_owner_input_error"
    )
  }
  r_ctrl <- validate_ipw_owner_control(reader_control)

  if (!is.list(links) || is.object(links) || length(links) < 1L) {
    cdid_ipw_covariance_transfer_abort(
      "`links` must be a non-empty plain list of sealed compact-view `cdid_ipw_link_enclosure` objects.",
      c("cdid_ipw_covariance_transfer_link_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  # Bound the header vectors before allocating/scanning one header per link.
  ipw_cov_transfer_safe_product(
    66,
    length(links),
    ctrl$max_matrix_entries,
    "link header vectors and per-link summaries before compact header validation"
  )
  header_budget <- ipw_cov_transfer_view_header_budget(
    links, frame_owner, candidate_rows, C, op, V_candidate, ctrl)
  frame_hdr <- ipw_cov_transfer_validate_frame_owner_header(frame_owner)
  link_hdr <- ipw_cov_transfer_validate_links_view_headers(links, frame_hdr)

  n <- frame_hdr$n
  L <- link_hdr$L
  expected_link_ids <- link_hdr$link_ids
  f_owner <- frame_hdr$frame_owner
  p_owner <- frame_hdr$panel_owner

  ipw_cov_transfer_refuse_decorative_attrs(C, "C", allowed_attrs = c("dim", "dimnames"))
  if (is.null(C) || is.object(C) || !is.atomic(C) || !is.numeric(C) ||
      !is.matrix(C) || nrow(C) < 1L || ncol(C) != L) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`C` must be a plain numeric d x L matrix with L = ", L, " columns and d >= 1 rows."),
      c("cdid_ipw_covariance_transfer_map_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  d <- nrow(C)

  ipw_cov_transfer_refuse_decorative_attrs(V_candidate, "V_candidate", allowed_attrs = c("dim", "dimnames"))
  if (is.null(V_candidate) || is.object(V_candidate) || !is.atomic(V_candidate) ||
      !is.numeric(V_candidate) || !is.matrix(V_candidate) ||
      nrow(V_candidate) != d || ncol(V_candidate) != d) {
    cdid_ipw_covariance_transfer_abort(
      paste0("`V_candidate` must be a plain numeric d x d matrix (d = ", d, ")."),
      c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }

  op_spec <- ipw_cov_transfer_validate_operator_header(op, n = n)

  ipw_cov_transfer_refuse_decorative_attrs(
    candidate_rows,
    "candidate_rows",
    allowed_attrs = if (is.function(candidate_rows)) "srcref" else c("names", "class")
  )
  cand_spec <- if (is.function(candidate_rows)) {
    ipw_candidate_rows_spec(callback = candidate_rows, link_ids = NULL)
  } else if (is.list(candidate_rows)) {
    if (is.null(names(candidate_rows)) || any(!nzchar(names(candidate_rows))) || anyDuplicated(names(candidate_rows))) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` list must have unique non-empty field names.",
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    forbidden_cand <- c("certified", "trusted", "is_exact", "bypass_enclosure", "trust_metadata")
    if (any(names(candidate_rows) %in% forbidden_cand)) {
      cdid_ipw_covariance_transfer_abort(
        "Shortcut or certificate flags on `candidate_rows` are refused.",
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_metadata_error")
      )
    }
    extra_cand <- setdiff(names(candidate_rows), c("callback", "link_ids"))
    if (length(extra_cand) > 0L) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Unexpected fields (`", paste(extra_cand, collapse = "`, `"), "`) in `candidate_rows`."),
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
    ipw_candidate_rows_spec(callback = candidate_rows$callback, link_ids = candidate_rows$link_ids)
  } else {
    cdid_ipw_covariance_transfer_abort(
      "`candidate_rows` must be a callback function `function(a_start, a_end)` or an `ipw_candidate_rows_spec` list.",
      c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }

  # 2. Base preflight BEFORE ID metadata/range allocation, content scans, character byte scans, or callbacks
  n_clusters_hdr <- if (identical(op_spec$row_level, "cluster")) length(op_spec$cluster_order) else n
  preflight <- preflight_ipw_covariance_transfer_view(
    n = n,
    L = L,
    d = d,
    row_level = op_spec$row_level,
    n_clusters = n_clusters_hdr,
    p_vec = link_hdr$p_vec,
    m_vec = link_hdr$m_vec,
    certificate_controls = link_hdr$certificate_controls,
    control = ctrl,
    reader_control = r_ctrl,
    id_is_character = identical(link_hdr$id_type, "character"),
    cluster_is_character = !is.null(op_spec$cluster_id) && is.character(op_spec$cluster_id),
    existing_4bd_retained_entries = link_hdr$existing_4bd_retained_entries,
    header_key_work_ops = header_budget$key_work_bound
  )
  preflight$shared_owner_retained_entries <- frame_hdr$shared_owner_retained_entries
  preflight$shared_view_retained_entries <- link_hdr$shared_view_retained_entries

  b_in <- as.integer(preflight$effective_chunk_rows)
  b_op <- as.integer(preflight$effective_operator_block_rows)
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  r_max_wo <- r_ctrl$max_work_ops

  actual_validation_work_ops <- ipw_cov_transfer_safe_sum(c(
    ipw_cov_transfer_safe_product(66, L, max_wo, "completed link/owner header validation work"),
    8 * link_hdr$reference_comparisons, header_budget$key_work_bound),
    max_wo, "completed header and distinct-reference/key validation work")
  preflight$actual_reference_comparisons <- link_hdr$reference_comparisons
  preflight$header_key_bytes <- header_budget$key_bytes
  preflight$allocation_phases$header_keys <- header_budget$entries_bound
  ipw_cov_transfer_safe_sum(c(preflight$base_component_concurrent_entries,
    ceiling(2 * header_budget$key_bytes / 8)), max_me, "header key bytes plus stage buffers")

  # 3. Account for link_ids character temporary bytes and key work BEFORE comparing link_ids or C colnames
  link_char_info <- ipw_cov_transfer_account_char_bytes(
    vec = expected_link_ids,
    n_len = L,
    step = max(1L, b_in),
    label = "links$link_id",
    max_matrix_entries = max_me,
    base_concurrent_entries = preflight$base_component_concurrent_entries,
    char_bytes_so_far = 0
  )
  if (!is.null(cand_spec$link_ids)) {
    cand_lid_char <- ipw_cov_transfer_account_char_bytes(
      vec = cand_spec$link_ids,
      n_len = length(cand_spec$link_ids),
      step = max(1L, b_in),
      label = "candidate_rows$link_ids",
      max_matrix_entries = max_me,
      base_concurrent_entries = preflight$base_component_concurrent_entries,
      char_bytes_so_far = link_char_info$total_char_temp_bytes
    )
    link_char_info <- cand_lid_char
  }
  if (!is.null(colnames(C))) {
    c_col_char <- ipw_cov_transfer_account_char_bytes(
      vec = colnames(C),
      n_len = ncol(C),
      step = max(1L, b_in),
      label = "colnames(C)",
      max_matrix_entries = max_me,
      base_concurrent_entries = preflight$base_component_concurrent_entries,
      char_bytes_so_far = link_char_info$total_char_temp_bytes
    )
    link_char_info <- c_col_char
  }
  actual_validation_work_ops <- ipw_cov_transfer_safe_sum(
    c(
      actual_validation_work_ops,
      ipw_cov_transfer_safe_product(32, link_char_info$total_char_temp_bytes, max_wo, "link_id key comparison work")
    ),
    max_bound = max_wo,
    label = "validation work after link_id byte accounting"
  )

  if (!is.null(cand_spec$link_ids) && !identical(as.character(cand_spec$link_ids), expected_link_ids)) {
    cdid_ipw_covariance_transfer_abort(
      "`candidate_rows$link_ids` disagrees with the ordered `link_id` sequence of `links`.",
      c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }

  # 4. Validate C, V_candidate, compact link content, owner ID bytes, and operator content
  if (anyNA(C) || any(!is.finite(C))) {
    cdid_ipw_covariance_transfer_abort(
      "`C` must be a finite numeric d x L matrix without NA/NaN/Inf.",
      c("cdid_ipw_covariance_transfer_map_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  if (!is.null(colnames(C)) && !identical(colnames(C), expected_link_ids)) {
    cdid_ipw_covariance_transfer_abort(
      "`colnames(C)` disagrees with the ordered `link_id` sequence of `links`.",
      c("cdid_ipw_covariance_transfer_map_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  C_dbl <- matrix(as.double(C), nrow = d, ncol = L, dimnames = list(rownames(C), expected_link_ids))

  if (anyNA(V_candidate) || any(!is.finite(V_candidate))) {
    cdid_ipw_covariance_transfer_abort(
      "`V_candidate` must be a finite numeric d x d matrix without NA/NaN/Inf.",
      c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
    )
  }
  for (j in seq_len(d)) {
    if (!identical(as.double(V_candidate[, j]), as.double(V_candidate[j, ]))) {
      cdid_ipw_covariance_transfer_abort(
        "`V_candidate` must be an exactly symmetric numeric d x d matrix.",
        c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  if (!is.null(rownames(C))) {
    if ((!is.null(rownames(V_candidate)) && !identical(rownames(V_candidate), rownames(C))) ||
        (!is.null(colnames(V_candidate)) && !identical(colnames(V_candidate), rownames(C)))) {
      cdid_ipw_covariance_transfer_abort(
        "`V_candidate` row/column names disagree with `rownames(C)`.",
        c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_input_error")
      )
    }
  }
  V_cand_dbl <- matrix(as.double(V_candidate), nrow = d, ncol = d, dimnames = dimnames(V_candidate))

  ipw_cov_transfer_validate_links_view_content(link_hdr, frame_hdr)

  is_num_id <- identical(link_hdr$id_type, "numeric")
  expected_chunk_lid_char_bytes <- link_char_info$total_char_temp_bytes
  if (!is.null(cand_spec$link_ids) || !is.null(colnames(C))) {
    lid_only_info <- ipw_cov_transfer_account_char_bytes(
      vec = expected_link_ids,
      n_len = L,
      step = max(1L, b_in),
      label = "links$link_id",
      max_matrix_entries = max_me,
      base_concurrent_entries = preflight$base_component_concurrent_entries,
      char_bytes_so_far = 0
    )
    expected_chunk_lid_char_bytes <- lid_only_info$total_char_temp_bytes
  }
  owner_char_info <- list(
    total_char_temp_bytes = link_char_info$total_char_temp_bytes,
    total_char_temp_entries_8byte = link_char_info$total_char_temp_entries_8byte
  )
  actual_reader_id_byte_scratch_entries <- 0
  actual_reader_id_byte_work_ops <- 0
  max_chunk_paired_id_bytes <- 0
  total_pass_paired_id_bytes <- 0

  if (!is_num_id) {
    owner_char_info <- ipw_cov_transfer_account_char_bytes(
      vec = f_owner$a_id_key,
      n_len = n,
      step = b_in,
      label = "frame_owner$a_id",
      max_matrix_entries = max_me,
      base_concurrent_entries = preflight$base_component_concurrent_entries,
      char_bytes_so_far = link_char_info$total_char_temp_bytes
    )
    c_start <- 1L
    while (c_start <= n) {
      ipw_cov_transfer_check_interrupt()
      c_end <- as.integer(min(as.double(n), as.double(c_start) + b_in - 1))
      c_seq <- seq.int(c_start, c_end)
      ob <- ipw_owner_id_byte_metadata(p_owner$a_id_key[c_seq])
      paired_ch_bytes <- ipw_owner_safe_product(2, ob$total, label = "chunk paired expected_a_id bytes")
      if (paired_ch_bytes > max_chunk_paired_id_bytes) {
        max_chunk_paired_id_bytes <- paired_ch_bytes
      }
      total_pass_paired_id_bytes <- ipw_owner_safe_sum(
        c(total_pass_paired_id_bytes, paired_ch_bytes),
        label = "whole-pass paired expected_a_id bytes"
      )
      rm(ob, c_seq)
      c_start <- as.double(c_end) + 1
    }
    actual_reader_id_byte_scratch_entries <- ipw_owner_safe_product(
      32,
      max_chunk_paired_id_bytes,
      r_ctrl$max_matrix_entries,
      "reader character expected_a_id temporary slots"
    )
    ipw_owner_safe_sum(
      c(preflight$max_view_reader_concurrent_entries, actual_reader_id_byte_scratch_entries),
      max_bound = r_ctrl$max_matrix_entries,
      label = "reader concurrent entries including actual character expected_a_id temporaries"
    )
    ipw_cov_transfer_safe_sum(
      c(preflight$max_view_combined_concurrent_entries, actual_reader_id_byte_scratch_entries),
      max_bound = ctrl$link_control$max_matrix_entries,
      label = "compact view combined concurrent entries including actual character expected_a_id temporaries"
    )
    per_link_reader_id_work <- ipw_owner_safe_product(
      128,
      total_pass_paired_id_bytes,
      r_max_wo,
      "per-link reader character expected_a_id comparison work"
    )
    for (l in seq_len(L)) {
      ipw_cov_transfer_safe_sum(
        c(preflight$link_chunk_preflights[[l]]$whole_stream_combined_work_ops_bound, per_link_reader_id_work),
        max_bound = ctrl$link_control$max_work_ops,
        label = paste0("whole-stream combined compact work plus character ID work for link ", l)
      )
    }
    actual_reader_id_byte_work_ops <- ipw_owner_safe_product(
      L,
      per_link_reader_id_work,
      r_max_wo,
      "total reader character expected_a_id comparison work across links"
    )
  }

  op_norm <- ipw_cov_transfer_validate_operator_content(
    op_spec = op_spec,
    n = n,
    step = b_in,
    max_matrix_entries = max_me,
    base_concurrent_entries = preflight$base_component_concurrent_entries,
    char_bytes_so_far = owner_char_info$total_char_temp_bytes
  )

  # Update preflight with actual static character key bytes, expected callback chunk key bytes,
  # and reader ID-byte charges BEFORE any callback or outcome read executes.
  static_char_bytes <- op_norm$char_info$total_char_temp_bytes
  expected_max_chunk_char_bytes <- ipw_cov_transfer_safe_sum(
    c(2 * expected_chunk_lid_char_bytes, max_chunk_paired_id_bytes),
    label = "expected max chunk callback character bytes"
  )
  expected_total_pass_chunk_char_bytes <- ipw_cov_transfer_safe_sum(
    c(
      ipw_cov_transfer_safe_product(
        as.double(preflight$n_input_chunks),
        2 * expected_chunk_lid_char_bytes,
        max_wo,
        "whole-pass callback link_id and colnames(phi) bytes"
      ),
      total_pass_paired_id_bytes
    ),
    label = "expected whole-pass callback character bytes"
  )
  pre_stream_char_bytes <- ipw_cov_transfer_safe_sum(
    c(static_char_bytes, expected_max_chunk_char_bytes),
    label = "pre-stream static plus expected chunk character bytes"
  )
  pre_stream_char_entries <- ceiling(pre_stream_char_bytes / 8)
  static_extra_keys <- ipw_cov_transfer_safe_sum(
    c(pre_stream_char_entries, actual_reader_id_byte_scratch_entries),
    max_bound = max_me,
    label = "static extra key-byte temporary entries P_extra_keys"
  )
  preflight$character_temporary_bytes <- pre_stream_char_bytes
  preflight$character_temporary_entries_8byte <- pre_stream_char_entries
  preflight$reader_id_byte_scratch_entries <- actual_reader_id_byte_scratch_entries
  preflight$extra_key_entries <- static_extra_keys
  preflight$concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$base_component_concurrent_entries, static_extra_keys),
    max_bound = max_me,
    label = "concurrent entries after owner/operator character key-byte accounting"
  )
  preflight$concurrent_entries <- max(preflight$concurrent_entries, header_budget$entries_bound)
  preflight$concurrent_bytes_8double <- 8 * preflight$concurrent_entries
  preflight$total_with_existing_4bd_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$existing_4bd_retained_entries, preflight$concurrent_entries),
    label = "validated existing 4bd plus component concurrent payload"
  )
  preflight$total_with_existing_4bd_bytes_8double <- 8 * preflight$total_with_existing_4bd_entries

  preflight$reader_id_byte_work_ops_bound <- actual_reader_id_byte_work_ops
  preflight$total_reader_work_ops_bound <- ipw_owner_safe_sum(
    c(preflight$total_reader_base_work_ops_bound, actual_reader_id_byte_work_ops),
    max_bound = r_max_wo,
    label = "total reader work bound after character ID byte accounting"
  )
  ipw_cov_transfer_safe_sum(
    c(preflight$total_reader_work_ops_bound, 0),
    max_bound = max_wo,
    label = "total reader work bound against covariance max_work_ops"
  )
  preflight$total_compact_view_stream_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(
      preflight$total_4bd_stream_work_ops_bound,
      preflight$total_cache_copy_work_ops_bound,
      preflight$total_reader_work_ops_bound
    ),
    max_bound = max_wo,
    label = "total compact view stream work bound after character ID byte accounting"
  )
  preflight$key_validation_work_ops_bound <- ipw_cov_transfer_safe_product(
    32,
    ipw_cov_transfer_safe_sum(
      c(static_char_bytes, expected_total_pass_chunk_char_bytes),
      label = "pre-stream static plus expected whole-pass callback character bytes"
    ),
    max_wo,
    "static key-byte validation work bound"
  )
  preflight$validation_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(preflight$base_validation_work_ops_bound, preflight$key_validation_work_ops_bound),
    max_bound = max_wo,
    label = "validation work bound after static key-byte accounting"
  )
  preflight$combined_total_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(
      preflight$component_work_ops_bound,
      preflight$total_compact_view_stream_work_ops_bound,
      preflight$validation_work_ops_bound
    ),
    max_bound = max_wo,
    label = "combined total work bound after static key-byte accounting"
  )

  static_owner_op_char_bytes <- static_char_bytes - link_char_info$total_char_temp_bytes
  actual_validation_work_ops <- ipw_cov_transfer_safe_sum(
    c(
      actual_validation_work_ops,
      8 * (d * L + d * d),
      ipw_cov_transfer_safe_product(4, n, max_wo, "operator weight and membership validation work"),
      ipw_cov_transfer_safe_product(4, as.double(op_norm$n_clusters), max_wo, "cluster order and weight validation work"),
      ipw_cov_transfer_safe_product(32, static_owner_op_char_bytes, max_wo, "owner and operator key-byte validation work")
    ),
    max_bound = max_wo,
    label = "actual validation work before streaming"
  )

  # 5. Check whether any sealed source nuisance problem or link is INCONCLUSIVE
  uncert_links <- which(!vapply(link_hdr$links_sealed, function(lk) isTRUE(lk$certified), logical(1L)))
  uncert_probs <- which(!vapply(link_hdr$problems, function(pr) isTRUE(pr$certified), logical(1L)))
  if (length(uncert_links) > 0L || length(uncert_probs) > 0L) {
    reasons <- unique(unlist(lapply(seq_len(L), function(l) {
      c(
        if (!isTRUE(link_hdr$problems[[l]]$certified)) link_hdr$problems[[l]]$inconclusive_reasons else character(0L),
        if (!isTRUE(link_hdr$links_sealed[[l]]$certified)) link_hdr$links_sealed[[l]]$inconclusive_reasons else character(0L)
      )
    })))
    if (length(reasons) == 0L) reasons <- "source_link_or_problem_inconclusive"
    inc_obj <- structure(
      list(
        storage_mode = "compact_view",
        certified = FALSE,
        enclosure_certified = FALSE,
        pd_transferred = FALSE,
        status = "INCONCLUSIVE",
        inconclusive = TRUE,
        inconclusive_type = "source_link_inconclusive",
        inconclusive_reasons = reasons,
        separation_certified = FALSE,
        sampling_certified = FALSE,
        identification_certified = FALSE,
        raw_panel_eligibility_certified = FALSE,
        public_gate_adopted = FALSE,
        singularity_proved = FALSE,
        nonidentification_proved = FALSE,
        rank_deficiency_proved = FALSE,
        link_ids = expected_link_ids,
        n_a_units = as.double(n),
        n_links = as.double(L),
        n_moments = as.double(d),
        n_unique_nuisance_problems = as.integer(link_hdr$n_unique_problems),
        frame = list(
          storage_mode = "compact_view",
          normalizer_n = as.double(n),
          embedding = "identity",
          id_type = link_hdr$id_type,
          id_aliased_to_source = TRUE,
          frame_ref = frame_hdr$frame_ref,
          panel_ref = frame_hdr$panel_ref,
          larger_frame_transport_applied = FALSE,
          a_rows_rescaled = FALSE,
          zero_rows_invented = FALSE
        ),
        operator = list(
          purpose = op_norm$purpose,
          row_level = op_norm$row_level,
          n_clusters = as.double(op_norm$n_clusters),
          n_positive_weight_rows = as.double(op_norm$n_positive_weight_rows),
          n_zero_unit_declared_clusters = as.integer(op_norm$n_zero_unit_declared_clusters),
          centered = op_norm$centered,
          kappa = op_norm$kappa,
          case_weights_reapplied = FALSE
        ),
        C = C_dbl,
        V_candidate = V_cand_dbl,
        enclosure = NULL,
        epsilon_acc = NA_real_,
        v_upper = NA_real_,
        nu_eval = NA_real_,
        nu_root = NA_real_,
        nu_total = NA_real_,
        epsilon_V = NA_real_,
        error_budget = NULL,
        eigen_bounds = NULL,
        control = ctrl,
        reader_control = r_ctrl,
        preflight = preflight
      ),
      class = "cdid_ipw_covariance_transfer"
    )
    return(ipw_enclosure_seal(inc_obj, "covariance_transfer"))
  }

  # 6. Callback chunk fetcher that validates batch headers BEFORE copying or reading outcomes,
  #    binds expected_a_id only to the B-length owned A slice, optionally reuses one bounded
  #    scratch across compatible consecutive links on the exact same problem$ref and range,
  #    and explicitly releases every callback/link/scratch alias before its successor.
  max_chunk_char_bytes <- 0
  total_pass_chunk_char_bytes <- 0
  total_physical_4bd_ops <- 0
  total_reused_4bd_ops <- 0
  total_reader_work_ops <- 0
  total_cache_copy_work_ops <- 0
  total_compact_view_stream_work_ops <- 0
  total_endpoints_read <- 0
  scratch_allocations <- 0L
  scratch_reuse_count <- 0L

  all_cb_released <- TRUE
  all_lk_released <- TRUE
  all_sc_released <- TRUE
  all_id_bounded <- TRUE
  max_expected_a_id_len <- 0L

  fetch_chunk_from_links <- function(a_start, a_end) {
    fetch_env <- environment()
    b_len <- a_end - a_start + 1L
    a_seq <- seq.int(a_start, a_end)

    batch <- tryCatch(
      withCallingHandlers(
        cand_spec$callback(a_start, a_end),
        warning = function(w) {
          cdid_ipw_covariance_transfer_abort(
            paste0("`candidate_rows` callback emitted a warning: ", conditionMessage(w)),
            "cdid_ipw_covariance_transfer_callback_error",
            diagnostics = list(a_start = a_start, a_end = a_end)
          )
        }
      ),
      error = function(e) {
        if (inherits(e, "cdid_ipw_covariance_transfer_refusal") ||
            inherits(e, "cdid_ipw_enclosure_refusal") ||
            inherits(e, "cdid_ipw_owner_refusal")) {
          stop(e)
        }
        cdid_ipw_covariance_transfer_abort(
          paste0("`candidate_rows` callback failed: ", conditionMessage(e)),
          "cdid_ipw_covariance_transfer_callback_error",
          diagnostics = list(a_start = a_start, a_end = a_end)
        )
      }
    )

    # Refuse oversized or malformed callback outputs BEFORE any component copies or outcome reads
    if (is.null(batch) || !is.list(batch) || is.object(batch)) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` callback must return a plain named list (`phi`, `id`, `link_ids`).",
        "cdid_ipw_covariance_transfer_callback_error",
        diagnostics = list(a_start = a_start, a_end = a_end)
      )
    }
    # Batch field-name keys are bounded before any names/shortcut comparisons.
    callback_name_budget <- ipw_cov_transfer_view_key_budget(names(batch), max_me,
      max_wo, preflight$base_component_concurrent_entries,
      preflight$component_work_ops_bound + preflight$total_compact_view_stream_work_ops_bound +
        actual_validation_work_ops, 0, 64)
    callback_name_work <- callback_name_budget$work -
      preflight$component_work_ops_bound - preflight$total_compact_view_stream_work_ops_bound -
      actual_validation_work_ops
    ipw_cov_transfer_refuse_decorative_attrs(batch, "candidate_rows batch", allowed_attrs = "names")
    if (is.null(names(batch)) || any(!nzchar(names(batch))) || anyDuplicated(names(batch))) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` callback batch must have unique non-empty names.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    forbidden_batch <- c("certified", "trusted", "is_exact", "bypass_enclosure", "trust_metadata")
    if (any(names(batch) %in% forbidden_batch)) {
      cdid_ipw_covariance_transfer_abort(
        "Decorative certificate flags on `candidate_rows` callback output are refused.",
        c("cdid_ipw_covariance_transfer_callback_error", "cdid_ipw_covariance_transfer_metadata_error")
      )
    }
    if (!all(c("phi", "id", "link_ids") %in% names(batch))) {
      cdid_ipw_covariance_transfer_abort(
        "`candidate_rows` callback batch must contain `phi`, `id`, and `link_ids`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    extra_b <- setdiff(names(batch), c("phi", "id", "link_ids", "a_index"))
    if (length(extra_b) > 0L) {
      cdid_ipw_covariance_transfer_abort(
        paste0("Unexpected fields (`", paste(extra_b, collapse = "`, `"), "`) in `candidate_rows` callback batch."),
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }

    phi_raw <- batch[["phi"]]
    id_raw <- batch[["id"]]
    lids_raw <- batch[["link_ids"]]
    a_idx_raw <- batch[["a_index"]]
    ipw_cov_transfer_refuse_decorative_attrs(phi_raw, "batch$phi", allowed_attrs = c("dim", "dimnames"))
    ipw_cov_transfer_refuse_decorative_attrs(id_raw, "batch$id")
    ipw_cov_transfer_refuse_decorative_attrs(lids_raw, "batch$link_ids")

    if (is.null(phi_raw) || is.object(phi_raw) || !is.atomic(phi_raw) ||
        !is.numeric(phi_raw) || !is.matrix(phi_raw)) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$phi` from `candidate_rows` must be a plain numeric matrix.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    nr_p <- nrow(phi_raw)
    nc_p <- ncol(phi_raw)
    if (is.na(nr_p) || is.na(nc_p) ||
        nr_p > ctrl$max_chunk_rows ||
        nr_p > ctrl$link_control$max_chunk_rows ||
        nr_p > r_ctrl$max_chunk_rows ||
        as.double(nr_p) * as.double(nc_p) > max_me) {
      cdid_ipw_covariance_transfer_abort(
        "Oversized `batch$phi` returned by `candidate_rows` callback exceeds explicit chunk/matrix bounds (refused before copy or outcome read).",
        c("cdid_ipw_covariance_transfer_bound_error", "cdid_ipw_covariance_transfer_callback_error"),
        diagnostics = list(nrow = nr_p, ncol = nc_p, expected_rows = b_len, expected_cols = L)
      )
    }
    if (nr_p != b_len || nc_p != L) {
      cdid_ipw_covariance_transfer_abort(
        paste0("`batch$phi` dimensions (", nr_p, " x ", nc_p, ") disagree with requested chunk (", b_len, " x ", L, ")."),
        "cdid_ipw_covariance_transfer_callback_error",
        diagnostics = list(nrow = nr_p, ncol = nc_p, expected_rows = b_len, expected_cols = L)
      )
    }

    if (is.null(id_raw) || is.factor(id_raw) || is.object(id_raw) || !is.atomic(id_raw) ||
        !is.null(dim(id_raw)) || length(id_raw) != b_len || is.numeric(id_raw) != is_num_id ||
        !(is.numeric(id_raw) || is.character(id_raw))) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$id` from `candidate_rows` must be a plain vector of length `b_len` matching the source ID type.",
        c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_callback_error")
      )
    }
    if (is.null(lids_raw) || is.object(lids_raw) || !is.atomic(lids_raw) ||
        !is.null(dim(lids_raw)) || !is.character(lids_raw) || length(lids_raw) != L) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$link_ids` from `candidate_rows` must be a plain character vector of length `L`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    if (!is.null(a_idx_raw)) {
      ipw_cov_transfer_refuse_decorative_attrs(a_idx_raw, "batch$a_index")
      if (is.object(a_idx_raw) || !is.atomic(a_idx_raw) || !is.null(dim(a_idx_raw)) ||
          !is.numeric(a_idx_raw) || length(a_idx_raw) != b_len ||
          anyNA(a_idx_raw) || any(!is.finite(a_idx_raw)) ||
          !identical(as.double(a_idx_raw), as.double(a_seq))) {
        cdid_ipw_covariance_transfer_abort(
          "`batch$a_index` from `candidate_rows` disagrees with requested `[a_start, a_end]`.",
          "cdid_ipw_covariance_transfer_callback_error"
        )
      }
    }

    # Reserve all callback key work before content comparisons/scans below.
    callback_key_budget <- list(bytes = 0, work =
      preflight$component_work_ops_bound + preflight$total_compact_view_stream_work_ops_bound +
        actual_validation_work_ops + callback_name_work)
    for (key_vec in list(id_raw, lids_raw, colnames(phi_raw))) {
      callback_key_budget <- ipw_cov_transfer_view_key_budget(key_vec, max_me, max_wo,
        preflight$base_component_concurrent_entries + actual_reader_id_byte_scratch_entries,
        callback_key_budget$work, callback_key_budget$bytes, 64)
    }
    rm(key_vec, callback_key_budget, callback_name_budget)

    # Account for callback character temporary bytes and work BEFORE string comparisons
    base_plus_reader_id_scratch <- ipw_cov_transfer_safe_sum(
      c(preflight$base_component_concurrent_entries, actual_reader_id_byte_scratch_entries),
      max_bound = max_me,
      label = "base concurrent entries plus reader ID scratch"
    )
    chunk_char_bytes <- 0
    if (!is_num_id) {
      ch_id_info <- ipw_cov_transfer_account_char_bytes(
        vec = id_raw,
        n_len = b_len,
        step = b_len,
        label = "batch$id",
        max_matrix_entries = max_me,
        base_concurrent_entries = base_plus_reader_id_scratch,
        char_bytes_so_far = static_char_bytes
      )
      chunk_char_bytes <- ch_id_info$total_char_temp_bytes - static_char_bytes
    }
    ch_lid_info <- ipw_cov_transfer_account_char_bytes(
      vec = lids_raw,
      n_len = L,
      step = L,
      label = "batch$link_ids",
      max_matrix_entries = max_me,
      base_concurrent_entries = base_plus_reader_id_scratch,
      char_bytes_so_far = static_char_bytes + chunk_char_bytes
    )
    chunk_char_bytes <- ch_lid_info$total_char_temp_bytes - static_char_bytes
    if (!is.null(colnames(phi_raw))) {
      ch_col_info <- ipw_cov_transfer_account_char_bytes(
        vec = colnames(phi_raw),
        n_len = L,
        step = L,
        label = "colnames(batch$phi)",
        max_matrix_entries = max_me,
        base_concurrent_entries = base_plus_reader_id_scratch,
        char_bytes_so_far = static_char_bytes + chunk_char_bytes
      )
      chunk_char_bytes <- ch_col_info$total_char_temp_bytes - static_char_bytes
    }
    if (chunk_char_bytes > max_chunk_char_bytes) {
      max_chunk_char_bytes <<- chunk_char_bytes
    }
    total_pass_chunk_char_bytes <<- ipw_cov_transfer_safe_sum(
      c(total_pass_chunk_char_bytes, chunk_char_bytes),
      label = "whole-pass callback character comparison bytes"
    )

    chunk_val_work <- ipw_cov_transfer_safe_sum(
      c(
        ipw_cov_transfer_safe_product(4, b_len, max_wo, "chunk ID validation work"),
        32 * b_len * L + 16 * b_len + callback_name_work,
        ipw_cov_transfer_safe_product(2, L, max_wo, "chunk link_id validation work"),
        ipw_cov_transfer_safe_product(32, chunk_char_bytes, max_wo, "chunk callback key-byte comparison work")
      ),
      max_bound = max_wo,
      label = "chunk callback validation and key-byte work"
    )
    actual_validation_work_ops <<- ipw_cov_transfer_safe_sum(
      c(actual_validation_work_ops, chunk_val_work),
      max_bound = max_wo,
      label = "cumulative validation work including callback key-byte comparisons"
    )
    ipw_cov_transfer_safe_sum(
      c(
        preflight$component_work_ops_bound,
        preflight$total_compact_view_stream_work_ops_bound,
        actual_validation_work_ops
      ),
      max_bound = max_wo,
      label = "combined total work bound including callback key-byte work"
    )

    if (!identical(as.character(lids_raw), expected_link_ids)) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$link_ids` from `candidate_rows` disagrees with the ordered `link_id` sequence of `links`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }
    if (!is.null(colnames(phi_raw)) && !identical(colnames(phi_raw), expected_link_ids)) {
      cdid_ipw_covariance_transfer_abort(
        "`colnames(batch$phi)` disagrees with `batch$link_ids`.",
        "cdid_ipw_covariance_transfer_callback_error"
      )
    }

    batch_id_key <- if (is_num_id) {
      if (anyNA(id_raw) || any(!is.finite(id_raw))) {
        cdid_ipw_covariance_transfer_abort(
          "`batch$id` contains NA, NaN, or infinite numeric identifiers.",
          c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_callback_error")
        )
      }
      ipw_cov_transfer_view_numeric_id_key(id_raw)
    } else {
      as.character(id_raw)
    }
    if (!identical(batch_id_key, f_owner$a_id_key[a_seq])) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$id` from `candidate_rows` disagrees with `frame_owner$a_id_key[a_start:a_end]`.",
        c("cdid_ipw_covariance_transfer_id_error", "cdid_ipw_covariance_transfer_callback_error"),
        diagnostics = list(a_start = a_start, a_end = a_end)
      )
    }

    if (anyNA(phi_raw) || any(!is.finite(phi_raw))) {
      cdid_ipw_covariance_transfer_abort(
        "`batch$phi` from `candidate_rows` contains NA, NaN, or Inf.",
        c("cdid_ipw_covariance_transfer_candidate_error", "cdid_ipw_covariance_transfer_callback_error"),
        diagnostics = list(a_start = a_start, a_end = a_end)
      )
    }

    # Conversion phase: callback payload + optional integer-to-double vector +
    # matrix result <= 3BL. Neither error buffer exists until aliases are gone.
    Phi_block <- if (is.double(phi_raw)) phi_raw else
      matrix(as.double(phi_raw), nrow = b_len, ncol = L)
    expected_a_id_chunk <- f_owner$a_id[a_seq]
    if (length(expected_a_id_chunk) != b_len) {
      all_id_bounded <<- FALSE
      cdid_ipw_covariance_transfer_abort(
        "Internal reader `expected_a_id` slice length disagrees with chunk length `b_len`.",
        "cdid_ipw_covariance_transfer_integrity_error"
      )
    }
    if (length(expected_a_id_chunk) > max_expected_a_id_len) {
      max_expected_a_id_len <<- as.integer(length(expected_a_id_chunk))
    }

    # Explicitly release callback return object AND all field aliases before entering the link loop.
    rm(batch, phi_raw, id_raw, lids_raw, a_idx_raw, batch_id_key)
    if (any(c("batch", "phi_raw", "id_raw", "lids_raw", "a_idx_raw", "batch_id_key") %in% ls(fetch_env, all.names = TRUE))) {
      all_cb_released <<- FALSE
      cdid_ipw_covariance_transfer_abort(
        "Callback aliases were not released before link IF traversal.",
        "cdid_ipw_covariance_transfer_integrity_error"
      )
    }

    # Link phase: three adapter BL buffers. Reserve one additional BL for
    # integer callback backing even if its closure retains that input matrix.
    E_eval_block <- matrix(0.0, nrow = b_len, ncol = L)
    E_root_block <- matrix(0.0, nrow = b_len, ncol = L)

    chunk_4bd_ops <- 0
    max_direct_star_dist <- 0.0
    active_scratch <- NULL
    active_scratch_prob_ref <- NULL

    for (l in seq_len(L)) {
      if (any(c("ch_l", "sep_rows", "w_l", "phi_col_l", "r_ctrl_call") %in% ls(fetch_env, all.names = TRUE))) {
        all_lk_released <<- FALSE
        cdid_ipw_covariance_transfer_abort(
          "Previous link IF chunk or alias survived into the next link iteration.",
          "cdid_ipw_covariance_transfer_integrity_error"
        )
      }

      prob_l <- link_hdr$problems[[l]]
      # Release optional scratch immediately on problem change
      if (!is.null(active_scratch) && !identical(active_scratch_prob_ref, prob_l$ref)) {
        rm(active_scratch, active_scratch_prob_ref)
        active_scratch <- NULL
        active_scratch_prob_ref <- NULL
      }
      # Allocate one bounded scratch only when consecutive links share the exact same sealed problem$ref
      if (is.null(active_scratch) && l < L && identical(link_hdr$problems[[l + 1L]]$ref, prob_l$ref)) {
        active_scratch <- new_ipw_if_chunk_scratch(
          problem = prob_l,
          a_start = a_start,
          a_end = a_end,
          control = ctrl$link_control
        )
        active_scratch_prob_ref <- prob_l$ref
        scratch_allocations <<- scratch_allocations + 1L
      } else if (!is.null(active_scratch)) {
        ipw_validate_if_chunk_scratch(active_scratch, prob_l, a_start, a_end)
        if (!identical(active_scratch$problem_ref, prob_l$ref) ||
            !identical(active_scratch$a_start, as.integer(a_start)) ||
            !identical(active_scratch$a_end, as.integer(a_end))) {
          all_sc_released <<- FALSE
          cdid_ipw_covariance_transfer_abort(
            "Incompatible chunk scratch retained across distinct problem or range.",
            "cdid_ipw_covariance_transfer_integrity_error"
          )
        }
        scratch_reuse_count <<- scratch_reuse_count + 1L
      }

      phi_col_l <- as.double(Phi_block[, l])
      r_ctrl_call <- r_ctrl
      attr(r_ctrl_call, "expected_a_id") <- expected_a_id_chunk
      if (!is.null(trace_env)) {
        attr(r_ctrl_call, "trace_env") <- trace_env
      }

      ch_l <- if (!is.null(active_scratch)) {
        enclose_ipw_link_if_view_chunk_with_scratch(
          link_enclosure = link_hdr$links_sealed[[l]],
          a_start = a_start,
          a_end = a_end,
          control = ctrl$link_control,
          reader_control = r_ctrl_call,
          scratch = active_scratch,
          phi_candidate = phi_col_l,
          include_components = FALSE,
          include_shared_xi = FALSE
        )
      } else {
        enclose_ipw_link_if_view_chunk(
          link_enclosure = link_hdr$links_sealed[[l]],
          a_start = a_start,
          a_end = a_end,
          control = ctrl$link_control,
          reader_control = r_ctrl_call,
          phi_candidate = phi_col_l,
          include_components = FALSE,
          include_shared_xi = FALSE
        )
      }

      if (!isTRUE(ch_l$certified)) {
        inc_reasons_l <- unique(c("link_if_view_chunk_inconclusive", ch_l$inconclusive_reasons))
        rm(ch_l, phi_col_l, r_ctrl_call, prob_l)
        if (!is.null(active_scratch)) {
          rm(active_scratch, active_scratch_prob_ref)
        }
        rm(expected_a_id_chunk, a_seq, Phi_block, E_eval_block, E_root_block)
        return(list(
          certified = FALSE,
          inconclusive_type = "source_chunk_inconclusive",
          inconclusive_reasons = inc_reasons_l,
          diagnostics = list(link_index = l, link_id = expected_link_ids[l], a_start = a_start, a_end = a_end)
        ))
      }

      sep_rows <- ch_l$separated_row_errors
      w_l <- ch_l$work
      if (is.null(sep_rows) || is.null(w_l) ||
          !identical(ch_l$a_index, a_seq) ||
          !identical(sep_rows$phi_candidate, phi_col_l) ||
          length(sep_rows$eval_error) != b_len ||
          length(sep_rows$root_error) != b_len ||
          anyNA(sep_rows$eval_error) || any(!is.finite(sep_rows$eval_error)) || any(sep_rows$eval_error < 0.0) ||
          anyNA(sep_rows$root_error) || any(!is.finite(sep_rows$root_error)) || any(sep_rows$root_error < 0.0)) {
        cdid_ipw_covariance_transfer_abort(
          paste0("`enclose_ipw_link_if_view_chunk` output for link ", l, " lacks verified `separated_row_errors` or `a_index`."),
          "cdid_ipw_covariance_transfer_integrity_error"
        )
      }

      # Never relabel direct_star_error as root_error:
      # eval_error bounds discrepancy from phi(beta0); root_error bounds change from beta0 to beta*.
      E_eval_block[, l] <- sep_rows$eval_error
      E_root_block[, l] <- sep_rows$root_error
      if (!is.null(sep_rows$direct_star_error)) {
        max_direct_star_dist <- max(max_direct_star_dist, max(sep_rows$direct_star_error))
      }

      # Preserve chunk_4bd_ops as interval-only; accumulate reader and cache-copy work separately.
      chunk_4bd_ops <- chunk_4bd_ops + w_l$actual_interval_ops
      total_physical_4bd_ops <<- total_physical_4bd_ops + w_l$physical_interval_ops
      total_reused_4bd_ops <<- total_reused_4bd_ops + w_l$reused_shared_interval_ops
      total_reader_work_ops <<- total_reader_work_ops + w_l$actual_reader_work_ops
      total_cache_copy_work_ops <<- total_cache_copy_work_ops + w_l$cache_copy_work_ops_charged
      total_compact_view_stream_work_ops <<- total_compact_view_stream_work_ops + w_l$combined_work_ops
      total_endpoints_read <<- total_endpoints_read + w_l$endpoints_read

      # Release scratch immediately if this is the last link or if the next link has a different problem$ref
      if (!is.null(active_scratch) && (l == L || !identical(link_hdr$problems[[l + 1L]]$ref, prob_l$ref))) {
        rm(active_scratch, active_scratch_prob_ref)
        active_scratch <- NULL
        active_scratch_prob_ref <- NULL
      }

      # Explicitly release current IF chunk and all aliases before the next link iteration
      rm(ch_l, sep_rows, w_l, phi_col_l, r_ctrl_call, prob_l)
    }

    rm(active_scratch, active_scratch_prob_ref, expected_a_id_chunk, a_seq)
    if (any(c("active_scratch", "active_scratch_prob_ref", "expected_a_id_chunk", "ch_l", "sep_rows") %in% ls(fetch_env, all.names = TRUE))) {
      all_sc_released <<- FALSE
      cdid_ipw_covariance_transfer_abort(
        "Chunk scratch or link aliases were not released at chunk completion.",
        "cdid_ipw_covariance_transfer_integrity_error"
      )
    }

    list(
      certified = TRUE,
      Phi_block = Phi_block,
      E_eval_block = E_eval_block,
      E_root_block = E_root_block,
      chunk_4bd_ops = chunk_4bd_ops,
      max_direct_star_dist = max_direct_star_dist
    )
  }

  # 7. Run the unchanged c45 streaming Gram and covariance-transfer kernel
  kern <- ipw_cov_transfer_stream_kernel(
    n = n,
    L = L,
    d = d,
    b_in = b_in,
    b_op = b_op,
    fetch_chunk = fetch_chunk_from_links,
    C_dbl = C_dbl,
    op_norm = op_norm,
    V_cand_dbl = V_cand_dbl,
    max_matrix_entries = max_me,
    preflight = preflight
  )

  # 8. Finalize preflight key-byte accounting and enforce outer combined work bounds
  total_char_bytes <- static_char_bytes + max(expected_max_chunk_char_bytes, max_chunk_char_bytes)
  total_char_entries <- ceiling(total_char_bytes / 8)
  final_extra_keys <- ipw_cov_transfer_safe_sum(
    c(total_char_entries, actual_reader_id_byte_scratch_entries),
    max_bound = max_me,
    label = "final extra key-byte temporary entries P_extra_keys"
  )
  preflight$character_temporary_bytes <- total_char_bytes
  preflight$character_temporary_entries_8byte <- total_char_entries
  preflight$extra_key_entries <- final_extra_keys
  preflight$concurrent_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$base_component_concurrent_entries, final_extra_keys),
    max_bound = max_me,
    label = "final concurrent entries P_cov + max_l P_view,l + P_extra_keys"
  )
  preflight$concurrent_entries <- max(preflight$concurrent_entries, header_budget$entries_bound)
  preflight$concurrent_bytes_8double <- 8 * preflight$concurrent_entries
  preflight$total_with_existing_4bd_entries <- ipw_cov_transfer_safe_sum(
    c(preflight$existing_4bd_retained_entries, preflight$concurrent_entries),
    label = "validated existing 4bd plus component concurrent payload"
  )
  preflight$total_with_existing_4bd_bytes_8double <- 8 * preflight$total_with_existing_4bd_entries

  preflight$key_validation_work_ops_bound <- ipw_cov_transfer_safe_product(
    32,
    static_char_bytes + max(expected_total_pass_chunk_char_bytes, total_pass_chunk_char_bytes),
    max_wo,
    "final key-byte validation work bound"
  )
  preflight$validation_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(preflight$base_validation_work_ops_bound, preflight$key_validation_work_ops_bound),
    max_bound = max_wo,
    label = "final validation work bound"
  )
  preflight$combined_total_work_ops_bound <- ipw_cov_transfer_safe_sum(
    c(
      preflight$component_work_ops_bound,
      preflight$total_compact_view_stream_work_ops_bound,
      preflight$validation_work_ops_bound
    ),
    max_bound = max_wo,
    label = "final combined total work bound"
  )

  work_out <- if (!is.null(kern$work)) {
    actual_comp_ops <- kern$work$actual_component_interval_ops
    actual_4bd_ops <- kern$work$actual_4bd_interval_ops
    stopifnot(
      actual_4bd_ops + total_reader_work_ops + total_cache_copy_work_ops == total_compact_view_stream_work_ops,
      total_reader_work_ops <= preflight$total_reader_work_ops_bound,
      total_reader_work_ops <= r_max_wo,
      total_cache_copy_work_ops <= preflight$total_cache_copy_work_ops_bound,
      total_compact_view_stream_work_ops <= preflight$total_compact_view_stream_work_ops_bound,
      actual_validation_work_ops <= preflight$validation_work_ops_bound,
      actual_comp_ops + total_compact_view_stream_work_ops + actual_validation_work_ops <= preflight$combined_total_work_ops_bound,
      actual_comp_ops + total_compact_view_stream_work_ops + actual_validation_work_ops <= max_wo
    )
    c(
      kern$work,
      list(
        physical_4bd_interval_ops = total_physical_4bd_ops,
        reused_shared_4bd_interval_ops = total_reused_4bd_ops,
        actual_reader_work_ops = total_reader_work_ops,
        cache_copy_work_ops_charged = total_cache_copy_work_ops,
        actual_compact_view_stream_work_ops = total_compact_view_stream_work_ops,
        actual_validation_work_ops = actual_validation_work_ops,
        combined_total_work_ops = actual_comp_ops + total_compact_view_stream_work_ops + actual_validation_work_ops,
        endpoints_read = logit_cert_count_output(total_endpoints_read),
        scratch_allocations = logit_cert_count_output(scratch_allocations),
        scratch_reuse_count = logit_cert_count_output(scratch_reuse_count)
      )
    )
  } else {
    list(
      physical_4bd_interval_ops = total_physical_4bd_ops,
      reused_shared_4bd_interval_ops = total_reused_4bd_ops,
      actual_reader_work_ops = total_reader_work_ops,
      cache_copy_work_ops_charged = total_cache_copy_work_ops,
      actual_compact_view_stream_work_ops = total_compact_view_stream_work_ops,
      actual_validation_work_ops = actual_validation_work_ops,
      combined_total_work_ops = total_compact_view_stream_work_ops + actual_validation_work_ops,
      endpoints_read = logit_cert_count_output(total_endpoints_read),
      scratch_allocations = logit_cert_count_output(scratch_allocations),
      scratch_reuse_count = logit_cert_count_output(scratch_reuse_count)
    )
  }

  res_obj <- structure(
    list(
      storage_mode = "compact_view",
      certified = isTRUE(kern$certified),
      enclosure_certified = isTRUE(kern$enclosure_certified),
      pd_transferred = isTRUE(kern$pd_transferred),
      status = kern$status,
      inconclusive = isTRUE(kern$inconclusive),
      inconclusive_type = kern$inconclusive_type,
      inconclusive_reasons = kern$inconclusive_reasons,
      separation_certified = FALSE,
      sampling_certified = FALSE,
      identification_certified = FALSE,
      raw_panel_eligibility_certified = FALSE,
      public_gate_adopted = FALSE,
      singularity_proved = FALSE,
      nonidentification_proved = FALSE,
      rank_deficiency_proved = FALSE,
      link_ids = expected_link_ids,
      n_a_units = as.double(n),
      n_links = as.double(L),
      n_moments = as.double(d),
      n_unique_nuisance_problems = as.integer(link_hdr$n_unique_problems),
      frame = list(
        storage_mode = "compact_view",
        normalizer_n = as.double(n),
        embedding = "identity",
        id_type = link_hdr$id_type,
        id_aliased_to_source = TRUE,
        frame_ref = frame_hdr$frame_ref,
        panel_ref = frame_hdr$panel_ref,
        larger_frame_transport_applied = FALSE,
        a_rows_rescaled = FALSE,
        zero_rows_invented = FALSE
      ),
      operator = list(
        purpose = op_norm$purpose,
        row_level = op_norm$row_level,
        n_clusters = as.double(op_norm$n_clusters),
        n_positive_weight_rows = as.double(op_norm$n_positive_weight_rows),
        n_zero_unit_declared_clusters = as.integer(op_norm$n_zero_unit_declared_clusters),
        centered = op_norm$centered,
        kappa = op_norm$kappa,
        case_weights_reapplied = FALSE
      ),
      C = C_dbl,
      V_candidate = V_cand_dbl,
      enclosure = kern$enclosure,
      S_unnormalized_iv = kern$S_unnormalized_iv,
      t_iv = kern$t_iv,
      W_iv = kern$W_iv,
      epsilon_acc = if (!is.null(kern$epsilon_acc)) kern$epsilon_acc else NA_real_,
      v_candidate_norm_upper = if (!is.null(kern$v_candidate_norm_upper)) kern$v_candidate_norm_upper else NA_real_,
      v_upper = if (!is.null(kern$v_upper)) kern$v_upper else NA_real_,
      nu_eval = if (!is.null(kern$nu_eval)) kern$nu_eval else NA_real_,
      nu_root = if (!is.null(kern$nu_root)) kern$nu_root else NA_real_,
      nu_total = if (!is.null(kern$nu_total)) kern$nu_total else NA_real_,
      epsilon_V = if (!is.null(kern$epsilon_V)) kern$epsilon_V else NA_real_,
      error_budget = kern$error_budget,
      eigen_bounds = kern$eigen_bounds,
      row_error_summary = kern$row_error_summary,
      cluster_summary = kern$cluster_summary,
      arithmetic_refusal = if (!is.null(kern$arithmetic_refusal)) {
        kern$arithmetic_refusal
      } else {
        list(occurred = FALSE, reason = "none", message = NULL, condition_class = NULL, diagnostics = NULL)
      },
      chunk_diagnostics = kern$chunk_diagnostics,
      control = ctrl,
      reader_control = r_ctrl,
      preflight = preflight,
      work = work_out,
      lifetime = list(
        callback_aliases_released_before_links = all_cb_released,
        link_results_released_each_step = all_lk_released,
        scratch_released_on_problem_or_chunk_end = all_sc_released,
        expected_a_id_chunk_bounded = all_id_bounded,
        max_expected_a_id_length = max_expected_a_id_len
      ),
      provenance = list(
        storage_mode = "compact_view",
        frame_ref = frame_hdr$frame_ref,
        panel_ref = frame_hdr$panel_ref,
        source_link_refs = lapply(link_hdr$links_sealed, `[[`, "ref"),
        source_problem_refs = lapply(link_hdr$problems, `[[`, "ref"),
        source_nuisance_view_refs = lapply(link_hdr$nuisance_views, `[[`, "ref"),
        source_link_view_refs = lapply(link_hdr$link_views, `[[`, "ref"),
        sep_rows_semantics = "eval_error bounds discrepancy from phi(beta0); root_error is the Theorem B E_root bound from beta0 to beta*; direct_star_error is never relabeled as root_error.",
        conditional_scope = "Certifies the finite covariance functional and sufficient Weyl PD transfer conditional on supplied sealed compact-view 4bd link enclosures on the identical sealed cdid_ipw_frame_owner. Does not certify fitted roots/loadings, sampling, raw panel eligibility, missing-trends/PT/overlap, population identification, larger-frame transport, estimated aggregation-weight IF, GMM/multiplier transfer, coverage, or process RSS."
      )
    ),
    class = "cdid_ipw_covariance_transfer"
  )

  ipw_enclosure_seal(res_obj, "covariance_transfer")
}
