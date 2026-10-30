# Private Standalone Compact Raw-IPW Backend (D93; DIRECTION.md and CLIPPING_INTERFACE.md).
#
# Scope and boundaries:
#   * Standalone OR-free compact raw-IPW producer and batched influence-function
#     evaluator operating on the admissible pair-eligible roster A (1..n_a_units)
#     with n_A normalization (internal lambda = 1) and original full-frame
#     transport (`transport_if_to_full_frame`) afterward.
#   * Borrows support `(pair, cohort)` and treatment/link `(pair, g, ordered K)`
#     deduplication and batched IF evaluation organization from DR, while using
#     ZERO outcome regression (no OR recipe, no WLS, no OR rank/fit/score, no
#     full-frame covariate requirement, no `cp04_frame_reader`).
#   * Calls the unchanged standalone IPW numerical and D93 clipping leaves in
#     `R/chained_ipw_link.R` (`ipw_nuisance_observer_fit_worker`,
#     `ipw_nuisance_hessian_stage`, `ipw_link_means`, `ipw_link_partials`,
#     `ipw_logit_values`, `ipw_clip_floor`, `ipw_clip_odds`,
#     `ipw_owner_read_y_positions`).
#   * Carries at most ONE canonical observed-risk buffer (`risk_a_idx`,
#     `risk_gi`, `risk_w`, `Z_use`) from outcome-free nuisance fitting through
#     endpoint outcome reads and point link means/partials, releasing it before
#     the next model and before sealing the backend.

ipw_raw_backend_abort <- function(message, class = "input", extra = list()) {
  stop(structure(
    c(list(message = message, call = NULL), extra),
    class = unique(c(
      paste0("cdid_ipw_raw_backend_", class, "_error"),
      paste0("cdid_cp04_", class, "_error"),
      "cdid_ipw_raw_backend_error",
      "cdid_ipw_fit_error",
      "cdid_ipw_error",
      "cdid_error",
      "error",
      "condition"
    ))
  ))
}

# Lightweight 1-binding locked reference environment for diagnostic/provenance
# compatibility without deep recursive per-link snapshot cloning. The whole
# compact backend is sealed once via `ipw_owner_seal(backend, "ipw_raw_backend")`.
ipw_raw_fast_ref <- function(obj, tag) {
  ref <- new.env(parent = emptyenv())
  obj$ref <- ref
  assign(tag, obj, envir = ref)
  lockBinding(tag, ref)
  lockEnvironment(ref, bindings = TRUE)
  obj
}

ipw_raw_fast_ref_ok <- function(obj, tag) {
  is.list(obj) && is.environment(obj$ref) &&
    identical(parent.env(obj$ref), emptyenv()) &&
    environmentIsLocked(obj$ref) &&
    exists(tag, envir = obj$ref, inherits = FALSE) &&
    bindingIsLocked(tag, obj$ref)
}

validate_ipw_raw_backend <- function(backend, frame_owner = NULL, panel_owner = NULL) {
  # Check the data binding before any getter can execute an active binding or
  # an S3 method. No native address is carried across a caller callback.
  ref_ok <- function(ref, tag) {
    is.environment(ref) && identical(parent.env(ref), emptyenv()) &&
      environmentIsLocked(ref) && exists(tag, ref, inherits = FALSE) &&
      bindingIsLocked(tag, ref) && !bindingIsActive(tag, ref)
  }
  if (missing(backend) || typeof(backend) != "list" ||
      !identical(attr(backend, "class", exact = TRUE), "cdid_ipw_raw_backend") ||
      !ref_ok(.subset2(backend, "ref"), "ipw_raw_backend") ||
      !gmm_sexp_aliases(backend, get("ipw_raw_backend", .subset2(backend, "ref"), inherits = FALSE))) {
    ipw_raw_backend_abort(
      "Compact raw-IPW backend is malformed, foreign, or differs from its sealed snapshot.",
      "integrity"
    )
  }
  if (!is.null(backend$reader_ref) || !is.null(backend$reference_ref) ||
      !is.null(backend$selection_backend) || !is.null(backend$design_ref) ||
      !is.null(backend$layout)) {
    ipw_raw_backend_abort(
      "Compact raw-IPW backend must not carry a CP04 full-frame reader, reference, selection, or layout owner.",
      "integrity"
    )
  }
  fref <- backend$frame_ref
  pref <- backend$panel_ref
  if (!ref_ok(fref, "frame_owner") || !ref_ok(pref, "panel_owner")) {
    ipw_raw_backend_abort(
      "Compact raw-IPW backend lacks valid sealed `frame_owner` and `panel_owner` references.",
      "integrity"
    )
  }
  f_obj <- get("frame_owner", envir = fref, inherits = FALSE)
  p_obj <- get("panel_owner", envir = pref, inherits = FALSE)
  if (typeof(f_obj) != "list" || typeof(p_obj) != "list" ||
      !identical(attr(f_obj, "class", exact = TRUE), c("cdid_ipw_frame_owner", "cdid_ipw_context")) ||
      !identical(attr(p_obj, "class", exact = TRUE), "cdid_ipw_panel_owner") ||
      !ref_ok(.subset2(f_obj, "ref"), "frame_owner") ||
      !ref_ok(.subset2(p_obj, "ref"), "panel_owner")) {
    ipw_raw_backend_abort("Compact raw-IPW owner snapshot has an invalid data binding.", "integrity")
  }
  # A canonical binding can be unlocked and rebound by a caller callback.
  # Its current owner cannot certify itself through that same mutable ref.
  # These construction-time borrowed objects keep the original data/control
  # independently reachable; unchanged objects require only the alias check.
  original <- backend$owner_snapshot
  if (typeof(original) != "list" || is.object(original) ||
      !identical(names(original), c("frame", "panel"))) {
    ipw_raw_backend_abort("Compact raw-IPW backend lacks its original owner references.", "integrity")
  }
  for (nm in c("frame", "panel")) {
    current <- if (nm == "frame") f_obj else p_obj
    if (!gmm_sexp_aliases(current, original[[nm]]) &&
        !identical(current, original[[nm]], num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)) {
      cdid_ipw_owner_abort(
        "Compact raw-IPW owner data or control differs from its construction snapshot; reuse is refused.",
        "cdid_ipw_owner_integrity_error"
      )
    }
  }
  validate_ipw_frame_owner(f_obj, panel_owner = p_obj)
  if (!identical(pref, f_obj$panel_ref)) {
    ipw_raw_backend_abort(
      "Compact raw-IPW backend `panel_ref` disagrees with its `frame_owner`.",
      "integrity"
    )
  }
  if (!is.null(frame_owner)) {
    if (!gmm_sexp_aliases(frame_owner, f_obj)) validate_ipw_frame_owner(frame_owner)
    if (!identical(fref, frame_owner$ref)) {
      ipw_raw_backend_abort(
        "Compact raw-IPW backend belongs to a different `cdid_ipw_frame_owner` snapshot.",
        "integrity"
      )
    }
  }
  if (!is.null(panel_owner)) {
    if (!gmm_sexp_aliases(panel_owner, p_obj)) validate_ipw_panel_owner(panel_owner)
    if (!identical(pref, panel_owner$ref)) {
      ipw_raw_backend_abort(
        "Compact raw-IPW backend belongs to a different `cdid_ipw_panel_owner` snapshot.",
        "integrity"
      )
    }
  }
  validate_ipw_owner_control(backend$control)
  validate_ipw_glm_control(backend$glm_control)
  if (!identical(as.integer(backend$n_a_units), as.integer(f_obj$n_a_units)) ||
      !identical(as.integer(backend$p), as.integer(f_obj$n_covariates)) ||
      !is.list(backend$pairs) || !is.list(backend$support_specs) ||
      !is.list(backend$supports) || length(backend$support_specs) != length(backend$supports) ||
      !is.list(backend$treatment_specs) || !is.list(backend$treatments) ||
      length(backend$treatment_specs) != length(backend$treatments) ||
      !is.list(backend$links) || length(backend$links) < 1L) {
    ipw_raw_backend_abort(
      "Compact raw-IPW backend tables or dimensions are inconsistent with its sealed owners.",
      "integrity"
    )
  }
  invisible(backend)
}

# Scan ONE unique date pair `(s, t)` across `1..n_a_units` in bounded chunks of
# `step`, computing outcome-free per-cohort counts, weight masses, and exact
# 4-double mask power-sum fingerprints modulo P = 67108859 for every A-cohort.
ipw_raw_scan_pair_summary <- function(frame_owner, panel_owner, pair_desc, step, date_memo = NULL) {
  n_a <- as.integer(frame_owner$n_a_units)
  c_a <- length(frame_owner$a_cohorts)
  P <- ipw_owner_mask_fingerprint_modulus
  a_cohort_idx <- frame_owner$a_cohort_index
  w_a <- frame_owner$w_a

  cohort_pair_count <- integer(c_a)
  cohort_pair_mass <- numeric(c_a)
  cohort_fp_s1 <- numeric(c_a)
  cohort_fp_s2 <- numeric(c_a)
  cohort_fp_s3 <- numeric(c_a)
  n_obs_pair_all <- 0L
  chunks <- 0L
  ops <- 0

  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  complete_pair <- !is.na(ds) && !is.na(dt) &&
    ipw_owner_pair_complete(panel_owner, pair_desc, step, date_memo)

  if (!is.na(ds) && !is.na(dt) &&
      panel_owner$postings$date_count[ds] > 0L &&
      panel_owner$postings$date_count[dt] > 0L) {
    a_start <- 1L
    while (a_start <= n_a) {
      ipw_owner_check_interrupt()
      a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
      hit <- if (complete_pair) {
        ipw_owner_query_complete_pair_mask(panel_owner, ds, dt, a_start, a_end)
      } else {
        ipw_owner_query_pair_mask(
          panel_owner, pair_desc$s, pair_desc$t,
          a_start = a_start, a_end = a_end, ds = ds, dt = dt
        )
      }
      a <- hit$a_index
      nh <- length(a)
      if (nh > 0L) {
        n_obs_pair_all <- n_obs_pair_all + nh
        ci <- a_cohort_idx[a]
        w_r <- w_a[a]
        u <- as.double(a) %% P
        u2 <- (u * u) %% P
        u3 <- (u2 * u) %% P
        group_order <- order(ci, method = "radix")
        runs <- rle(ci[group_order])
        group_end <- cumsum(runs$lengths)
        group_start <- group_end - runs$lengths + 1L
        for (z in seq_along(runs$values)) {
          j <- runs$values[z]
          if (!is.na(j) && j >= 1L && j <= c_a) {
            idx_z <- group_order[seq.int(group_start[z], group_end[z])]
            cohort_pair_count[j] <- cohort_pair_count[j] + runs$lengths[z]
            cohort_pair_mass[j] <- cohort_pair_mass[j] + sum(w_r[idx_z])
            cohort_fp_s1[j] <- (cohort_fp_s1[j] + ipw_owner_exact_mod_sum(u[idx_z], P)) %% P
            cohort_fp_s2[j] <- (cohort_fp_s2[j] + ipw_owner_exact_mod_sum(u2[idx_z], P)) %% P
            cohort_fp_s3[j] <- (cohort_fp_s3[j] + ipw_owner_exact_mod_sum(u3[idx_z], P)) %% P
          }
        }
      }
      chunks <- chunks + 1L
      ops <- ops + 4 * (a_end - a_start + 1L)
      a_start <- as.double(a_end) + 1
    }
  } else {
    chunks <- as.integer(ceiling(n_a / step))
    ops <- n_a
  }

  list(
    pair = pair_desc,
    complete_pair = complete_pair,
    n_observed_pair_all = as.integer(n_obs_pair_all),
    cohort_pair_count = cohort_pair_count,
    cohort_pair_mass = cohort_pair_mass,
    cohort_fp_s1 = cohort_fp_s1,
    cohort_fp_s2 = cohort_fp_s2,
    cohort_fp_s3 = cohort_fp_s3,
    chunks = chunks,
    ops = ops
  )
}

# Extract cohort-specific pair support and 4-double mask fingerprint from a
# scanned `pair_summary` record, enforcing the exact support refusals of
# `ipw_owner_chunked_pair_support`.
ipw_raw_pair_support_for_cohorts <- function(frame_owner, pair_summary, coh_desc, require_support = TRUE) {
  n_a <- as.integer(frame_owner$n_a_units)
  m <- coh_desc$m_cohorts
  ci <- coh_desc$a_cohort_map
  pair_count <- integer(m)
  pair_mass <- numeric(m)
  present <- which(!is.na(ci))
  if (length(present) > 0L) {
    pair_count[present] <- pair_summary$cohort_pair_count[ci[present]]
    pair_mass[present] <- pair_summary$cohort_pair_mass[ci[present]]
  }
  n_obs_pair_all <- pair_summary$n_observed_pair_all
  n_obs_risk <- sum(pair_count)

  if (isTRUE(require_support)) {
    if (n_obs_pair_all < 1L || n_obs_risk < 1L) {
      cdid_ipw_owner_abort(
        "Requested date pair has zero empirical risk-set observations; requested operator is unavailable.",
        c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error")
      )
    }
    if (any(pair_count < 1L) || any(!is.finite(pair_mass) | pair_mass <= 0)) {
      cdid_ipw_owner_abort(
        "At least one declared cohort in `{g, K}` has zero empirical pair count or mass on `(s, t)`; zero-support cohorts are not silently discarded.",
        c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error")
      )
    }
  }
  rates <- ifelse(coh_desc$group_a_mass > 0, pair_mass / coh_desc$group_a_mass, NA_real_)
  if (isTRUE(require_support)) {
    rates[pair_count == coh_desc$group_a_count] <- 1
  }
  if (isTRUE(require_support) && any(!is.finite(rates) | rates <= 0 | rates > 1)) {
    cdid_ipw_owner_abort(
      "Invalid marginal pair observation rates `r_j`.",
      c("cdid_ipw_owner_pair_support_error", "cdid_ipw_rate_error")
    )
  }
  names(pair_count) <- names(pair_mass) <- names(rates) <- coh_desc$cohort_labels

  P <- ipw_owner_mask_fingerprint_modulus
  fp <- c(
    count = as.double(n_obs_risk),
    s1 = if (length(present) > 0L) sum(pair_summary$cohort_fp_s1[ci[present]]) %% P else 0,
    s2 = if (length(present) > 0L) sum(pair_summary$cohort_fp_s2[ci[present]]) %% P else 0,
    s3 = if (length(present) > 0L) sum(pair_summary$cohort_fp_s3[ci[present]]) %% P else 0
  )

  list(
    n_a_units = n_a,
    n_observed_pair_all = n_obs_pair_all,
    n_risk_a_units = sum(coh_desc$group_a_count),
    n_outside_risk_a_units = n_a - sum(coh_desc$group_a_count),
    n_observed_risk_units = as.integer(n_obs_risk),
    n_nonobserved_risk_units = sum(coh_desc$group_a_count) - as.integer(n_obs_risk),
    by_cohort_a_count = coh_desc$group_a_count,
    by_cohort_a_mass = coh_desc$group_a_mass,
    by_cohort_pair_count = pair_count,
    by_cohort_pair_mass = pair_mass,
    q = coh_desc$q,
    rates = rates,
    fingerprint = fp,
    chunks = pair_summary$chunks,
    ops = pair_summary$ops + 4 * m
  )
}

# Exact chunked zero-outcome comparison of two pair masks on `coh_desc`.
# Short-circuits immediately when the pairs are identical or both complete.
ipw_raw_verify_pair_masks_equal <- function(frame_owner, panel_owner, pair1_summary, pair2_summary, coh_desc, step) {
  p1 <- pair1_summary$pair
  p2 <- pair2_summary$pair
  if (identical(p1$s, p2$s) && identical(p1$t, p2$t) &&
      identical(p1$s_calendar_index, p2$s_calendar_index) &&
      identical(p1$t_calendar_index, p2$t_calendar_index)) {
    return(TRUE)
  }
  if (isTRUE(pair1_summary$complete_pair) && isTRUE(pair2_summary$complete_pair)) {
    return(TRUE)
  }
  n_a <- as.integer(frame_owner$n_a_units)
  a_cohort_idx <- frame_owner$a_cohort_index
  oci <- coh_desc$owner_cohort_index
  ds1 <- p1$s_calendar_index; dt1 <- p1$t_calendar_index
  ds2 <- p2$s_calendar_index; dt2 <- p2$t_calendar_index

  a_start <- 1L
  while (a_start <= n_a) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    hit1 <- ipw_owner_query_pair_mask(panel_owner, p1$s, p1$t, a_start = a_start, a_end = a_end, ds = ds1, dt = dt1)
    hit2 <- ipw_owner_query_pair_mask(panel_owner, p2$s, p2$t, a_start = a_start, a_end = a_end, ds = ds2, dt = dt2)
    a1 <- hit1$a_index
    a2 <- hit2$a_index
    if (length(a1) > 0L) a1 <- a1[!is.na(ipw_owner_match_sorted_double(a_cohort_idx[a1], oci))]
    if (length(a2) > 0L) a2 <- a2[!is.na(ipw_owner_match_sorted_double(a_cohort_idx[a2], oci))]
    if (!identical(as.integer(a1), as.integer(a2))) return(FALSE)
    a_start <- as.double(a_end) + 1
  }
  TRUE
}

# Stream ONE observer buffer (`risk_a_idx`, `risk_gi`, `risk_w`, `Z_use`,
# `obs_s`, `obs_t`) for `(pair_summary, coh_desc)` in bounded chunks of `step`.
# No outcome values (`y`) are read here.
ipw_raw_collect_observer_buffer <- function(frame_owner, panel_owner, pair_summary, coh_desc, n_obs, step) {
  n_a <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  a_cohort_idx <- frame_owner$a_cohort_index
  w_a <- frame_owner$w_a
  Z_a <- frame_owner$Z_a
  pair_desc <- pair_summary$pair
  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  complete_pair <- isTRUE(pair_summary$complete_pair)

  risk_a_idx <- integer(n_obs)
  risk_gi <- integer(n_obs)
  risk_w <- numeric(n_obs)
  Z_use <- matrix(0.0, nrow = n_obs, ncol = p, dimnames = list(NULL, colnames(Z_a)))
  obs_s <- integer(n_obs)
  obs_t <- integer(n_obs)

  pos <- 0L
  a_start <- 1L
  while (a_start <= n_a) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    hit <- if (complete_pair) {
      ipw_owner_query_complete_pair_mask(panel_owner, ds, dt, a_start, a_end)
    } else {
      ipw_owner_query_pair_mask(
        panel_owner, pair_desc$s, pair_desc$t,
        a_start = a_start, a_end = a_end, ds = ds, dt = dt
      )
    }
    a <- hit$a_index
    if (length(a) > 0L) {
      gi <- ipw_owner_view_cohort_index(coh_desc, a_cohort_idx[a])
      in_r <- which(!is.na(gi))
      k <- length(in_r)
      if (k > 0L) {
        if (pos + k > n_obs) {
          ipw_raw_backend_abort("Observed risk support count changed from its pair summary.", "integrity")
        }
        ii <- seq.int(pos + 1L, pos + k)
        a_sub <- a[in_r]
        risk_a_idx[ii] <- a_sub
        risk_gi[ii] <- as.integer(gi[in_r])
        risk_w[ii] <- w_a[a_sub]
        Z_use[ii, ] <- Z_a[a_sub, , drop = FALSE]
        obs_s[ii] <- hit$obs_s[in_r]
        obs_t[ii] <- hit$obs_t[in_r]
        pos <- pos + k
      }
    }
    a_start <- as.double(a_end) + 1
  }
  if (pos != n_obs) {
    ipw_raw_backend_abort("Observed risk support count disagreed with its pair summary.", "integrity")
  }
  list(
    risk_a_idx = risk_a_idx,
    risk_gi = risk_gi,
    risk_w = risk_w,
    Z_use = Z_use,
    obs_s = obs_s,
    obs_t = obs_t
  )
}

# Collect endpoint observation indices `(obs_s, obs_t)` on a second pair
# `pair_summary` that shares the exact same observed risk set `expected_a_idx`.
ipw_raw_collect_pair_endpoints <- function(panel_owner, pair_summary, coh_desc, expected_a_idx, a_cohort_idx, step) {
  n_a <- length(a_cohort_idx)
  n_obs <- length(expected_a_idx)
  pair_desc <- pair_summary$pair
  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  if (isTRUE(pair_summary$complete_pair)) {
    post <- panel_owner$postings
    ps <- post$date_start[ds] - 1L + expected_a_idx
    pt <- post$date_start[dt] - 1L + expected_a_idx
    if (!identical(as.integer(post$post_a_index[ps]), as.integer(expected_a_idx)) ||
        !identical(as.integer(post$post_a_index[pt]), as.integer(expected_a_idx))) {
      ipw_raw_backend_abort("Complete-date posting slot disagreed with its A-unit index.", "integrity")
    }
    return(list(obs_s = post$post_obs_index[ps], obs_t = post$post_obs_index[pt]))
  }
  obs_s <- integer(n_obs)
  obs_t <- integer(n_obs)
  a_found <- integer(n_obs)
  pos <- 0L
  a_start <- 1L
  while (a_start <= n_a) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    hit <- ipw_owner_query_pair_mask(
      panel_owner, pair_desc$s, pair_desc$t,
      a_start = a_start, a_end = a_end, ds = ds, dt = dt
    )
    a <- hit$a_index
    if (length(a) > 0L) {
      gi <- ipw_owner_view_cohort_index(coh_desc, a_cohort_idx[a])
      in_r <- which(!is.na(gi))
      k <- length(in_r)
      if (k > 0L) {
        if (pos + k > n_obs) {
          ipw_raw_backend_abort("Shared-mask pair observer count disagreed with its nuisance branch.", "integrity")
        }
        ii <- seq.int(pos + 1L, pos + k)
        a_found[ii] <- a[in_r]
        obs_s[ii] <- hit$obs_s[in_r]
        obs_t[ii] <- hit$obs_t[in_r]
        pos <- pos + k
      }
    }
    a_start <- as.double(a_end) + 1
  }
  if (pos != n_obs || !identical(as.integer(a_found), as.integer(expected_a_idx))) {
    ipw_raw_backend_abort("Shared-mask pair observer indices disagreed with its nuisance branch.", "integrity")
  }
  list(obs_s = obs_s, obs_t = obs_t)
}

# Build lightweight diagnostic records (`nuisance_view`, `fitted_view`, `link_view`, `point_link`)
# for downstream consumers (`chained_request_nuisance_diagnostics`, `reduction_summary`, etc.).
ipw_raw_make_nuisance_view_record <- function(nuisance_id, frame_owner, pair_desc, coh_desc, sup_info, own_ctrl, preflight) {
  nv <- structure(
    list(
      nuisance_id = nuisance_id,
      frame_ref = frame_owner$ref,
      panel_ref = frame_owner$panel_ref,
      defining_pair = pair_desc,
      cohorts = coh_desc,
      support = list(
        n_a_units = frame_owner$n_a_units,
        n_risk_a_units = sup_info$n_risk_a_units,
        n_outside_risk_a_units = sup_info$n_outside_risk_a_units,
        n_observed_pair_all = sup_info$n_observed_pair_all,
        n_observed_risk_units = sup_info$n_observed_risk_units,
        n_nonobserved_risk_units = sup_info$n_nonobserved_risk_units,
        by_cohort_a_count = coh_desc$group_a_count,
        by_cohort_a_mass = coh_desc$group_a_mass,
        by_cohort_pair_count = sup_info$by_cohort_pair_count,
        by_cohort_pair_mass = sup_info$by_cohort_pair_mass,
        q = coh_desc$q,
        rates = sup_info$rates
      ),
      candidates = NULL,
      control = own_ctrl,
      preflight = preflight,
      accounting = list(
        full_n_masks_retained = FALSE,
        copied_Z_retained = FALSE,
        endpoints_read = 0L
      ),
      work = list(
        support_chunks_processed = logit_cert_count_output(sup_info$chunks),
        actual_work_ops = sup_info$ops
      )
    ),
    class = "cdid_ipw_nuisance_view"
  )
  ipw_raw_fast_ref(nv, "nuisance_view")
}

ipw_raw_make_link_view_record <- function(link_id, functional_key, frame_owner, pair_desc, coh_desc, sup_info,
                                          nv_rec, own_ctrl, preflight) {
  lv <- structure(
    list(
      link_id = link_id,
      functional_key = functional_key,
      nuisance_id = nv_rec$nuisance_id,
      nuisance_ref = nv_rec$ref,
      frame_ref = frame_owner$ref,
      panel_ref = frame_owner$panel_ref,
      pair = pair_desc,
      defining_pair = nv_rec$defining_pair,
      cohorts = coh_desc,
      nuisance_mask_verified = TRUE,
      support = list(
        n_a_units = frame_owner$n_a_units,
        n_risk_a_units = sup_info$n_risk_a_units,
        n_outside_risk_a_units = sup_info$n_outside_risk_a_units,
        n_observed_pair_all = sup_info$n_observed_pair_all,
        n_observed_risk_units = sup_info$n_observed_risk_units,
        n_nonobserved_risk_units = sup_info$n_nonobserved_risk_units,
        by_cohort_a_count = coh_desc$group_a_count,
        by_cohort_a_mass = coh_desc$group_a_mass,
        by_cohort_pair_count = sup_info$by_cohort_pair_count,
        by_cohort_pair_mass = sup_info$by_cohort_pair_mass,
        q = coh_desc$q,
        rates = sup_info$rates
      ),
      control = own_ctrl,
      preflight = preflight
    ),
    class = "cdid_ipw_link_view"
  )
  ipw_raw_fast_ref(lv, "link_view")
}

ipw_raw_make_fitted_view_record <- function(nv_rec, frame_owner, st_summary, hs, glm_ctrl, fit_ctrl, preflight) {
  na <- as.integer(frame_owner$n_a_units)
  p <- as.integer(frame_owner$n_covariates)
  coh <- nv_rec$cohorts
  m <- coh$m_cohorts
  n_obs <- as.integer(nv_rec$support$n_observed_risk_units)
  cohort_labels <- coh$cohort_labels

  st_pair_count <- as.integer(st_summary$pair_count)
  names(st_pair_count) <- cohort_labels
  st_pair_mass <- as.double(st_summary$pair_mass)
  names(st_pair_mass) <- cohort_labels
  st_rates <- as.double(st_summary$rates)
  names(st_rates) <- cohort_labels
  st_beta <- as.double(st_summary$beta)
  names(st_beta) <- colnames(frame_owner$Z_a)

  view_pair_mass <- as.double(nv_rec$support$by_cohort_pair_mass)
  names(view_pair_mass) <- cohort_labels
  view_rates <- as.double(nv_rec$support$rates)
  names(view_rates) <- cohort_labels
  cohort_a_mass <- as.double(nv_rec$support$by_cohort_a_mass)
  names(cohort_a_mass) <- cohort_labels
  cohort_q <- as.double(nv_rec$support$q)
  names(cohort_q) <- cohort_labels

  st_fit_info <- st_summary$fit_info
  fv <- structure(
    list(
      nuisance_id = nv_rec$nuisance_id,
      nuisance_view_ref = nv_rec$ref,
      frame_ref = nv_rec$frame_ref,
      panel_ref = nv_rec$panel_ref,
      defining_pair = nv_rec$defining_pair,
      cohorts = list(
        g = coh$g,
        K = coh$K,
        cohorts = coh$cohorts,
        cohort_labels = cohort_labels,
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
        n_risk_a_units = as.integer(nv_rec$support$n_risk_a_units),
        n_outside_risk_a_units = as.integer(nv_rec$support$n_outside_risk_a_units),
        n_observed_pair_units = as.integer(nv_rec$support$n_observed_pair_all),
        n_observed_risk_units = n_obs,
        n_nonobserved_risk_units = as.integer(nv_rec$support$n_nonobserved_risk_units),
        n_treated_obs = st_summary$n_treated_obs,
        n_control_obs = st_summary$n_control_obs,
        by_cohort_a = nv_rec$support$by_cohort_a_count,
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
        sum_a = st_summary$sum_a,
        mean_a = st_summary$mean_a,
        sum_b = st_summary$sum_b,
        mean_b = st_summary$mean_b
      ),
      coefficients = st_beta,
      beta = st_beta,
      factors = list(
        H = hs$H,
        H_chol = hs$H_chol,
        qr_R = st_fit_info$qr_R,
        qr_pivot = st_fit_info$qr_pivot,
        qr_tol = st_fit_info$qr_tol,
        rank = st_fit_info$rank
      ),
      partials = list(
        A_beta_r = hs$A_beta_r
      ),
      diagnostics = list(
        weight_concentration = st_summary$weight_concentration,
        weight_concentration_definition = "Kish(v) = (sum(v/max(v)))^2 / sum((v/max(v))^2); descriptive concentration only",
        weight_concentration_sources = c(
          propensity_fit_weights = "st$fit_weights",
          target_outcome_weights = "st$a_obs",
          pooled_control_outcome_ipw_weights = "st$b_obs"
        ),
        fitted_probability_range = st_summary$probability_range,
        fitted_probability_scope = "observed fitted risk set",
        clipping = st_fit_info$clipping,
        converged = st_fit_info$converged,
        boundary = st_fit_info$boundary,
        iter = st_fit_info$iter,
        rank = st_fit_info$rank,
        n_covariates = p,
        df_residual = st_fit_info$df_residual,
        deviance = st_fit_info$deviance,
        null_deviance = st_fit_info$null_deviance,
        aic = st_fit_info$aic,
        control = glm_ctrl
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
      control = fit_ctrl,
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
      )
    ),
    class = "cdid_ipw_fitted_nuisance_view"
  )
  ipw_raw_fast_ref(fv, "fitted_nuisance_view")
}

ipw_raw_make_point_link_record <- function(lv_rec, fv_rec, trt_fit, cp04_index, pt_ctrl, preflight) {
  eff_link_id <- lv_rec$link_id
  eff_func_key <- lv_rec$functional_key
  na <- fv_rec$dimensions$n_a_units
  n_obs <- fv_rec$dimensions$n_observed_risk_units
  cohort_labels <- fv_rec$cohorts$cohort_labels

  link_pair_mass <- as.double(lv_rec$support$by_cohort_pair_mass)
  names(link_pair_mass) <- cohort_labels
  link_rates <- as.double(lv_rec$support$rates)
  names(link_rates) <- cohort_labels
  shared_across_dates <- !identical(lv_rec$pair$s, fv_rec$defining_pair$s) ||
    !identical(lv_rec$pair$t, fv_rec$defining_pair$t)

  pl <- structure(
    list(
      link_id = eff_link_id,
      functional_key = eff_func_key,
      nuisance_id = fv_rec$nuisance_id,
      fitted_nuisance_ref = fv_rec$ref,
      nuisance_view_ref = fv_rec$nuisance_view_ref,
      nuisance_ref = fv_rec$nuisance_view_ref,
      link_view_ref = lv_rec$ref,
      link_ref = lv_rec$ref,
      frame_ref = fv_rec$frame_ref,
      panel_ref = fv_rec$panel_ref,
      cp04_index = as.integer(cp04_index),
      estimate = trt_fit$means$estimate,
      treated_mean = trt_fit$means$treated_mean,
      control_mean = trt_fit$means$control_mean,
      pair = list(
        s = lv_rec$pair$s,
        t = lv_rec$pair$t,
        s_calendar_index = lv_rec$pair$s_calendar_index,
        t_calendar_index = lv_rec$pair$t_calendar_index,
        link_id = eff_link_id,
        functional_key = eff_func_key
      ),
      defining_pair = fv_rec$defining_pair,
      shared_across_dates = shared_across_dates,
      nuisance_mask_verified = TRUE,
      cohorts = fv_rec$cohorts,
      dimensions = fv_rec$dimensions,
      counts = list(
        n_a_units = na,
        n_risk_a_units = fv_rec$counts$n_risk_a_units,
        n_outside_risk_a_units = fv_rec$counts$n_outside_risk_a_units,
        n_observed_pair_units = as.integer(lv_rec$support$n_observed_pair_all),
        n_observed_risk_units = n_obs,
        n_nonobserved_risk_units = fv_rec$counts$n_nonobserved_risk_units,
        n_treated_obs = fv_rec$counts$n_treated_obs,
        n_control_obs = fv_rec$counts$n_control_obs,
        by_cohort_a = fv_rec$counts$by_cohort_a,
        by_cohort_pair = fv_rec$counts$by_cohort_pair
      ),
      masses = fv_rec$masses,
      means = trt_fit$means,
      coefficients = fv_rec$coefficients,
      beta = fv_rec$beta,
      factors = fv_rec$factors,
      partials = trt_fit$partials,
      loadings = trt_fit$loadings,
      diagnostics = fv_rec$diagnostics,
      point_summary_disclosure = list(
        defining_rate_source = fv_rec$point_summary_disclosure$defining_rate_source,
        observer_cohort_a_mass = fv_rec$point_summary_disclosure$observer_cohort_a_mass,
        observer_cohort_pair_mass = fv_rec$point_summary_disclosure$observer_cohort_pair_mass,
        observer_rates = fv_rec$point_summary_disclosure$observer_rates,
        view_chunked_cohort_a_mass = fv_rec$point_summary_disclosure$view_chunked_cohort_a_mass,
        view_chunked_cohort_pair_mass = fv_rec$point_summary_disclosure$view_chunked_cohort_pair_mass,
        view_chunked_rates = fv_rec$point_summary_disclosure$view_chunked_rates,
        link_chunked_cohort_pair_mass = link_pair_mass,
        link_chunked_rates = link_rates,
        max_abs_pair_mass_chunk_diff = max(abs(fv_rec$masses$cohort_pair_weight - link_pair_mass)),
        max_abs_rate_chunk_diff = max(abs(fv_rec$masses$rates - link_rates)),
        propensity_refit_performed = FALSE,
        candidate_point_summary_only = TRUE,
        certifies_root_or_inference = FALSE,
        glm_convergence_is_certificate = FALSE,
        seal_is_certificate = FALSE
      ),
      glm_control = fv_rec$glm_control,
      control = pt_ctrl,
      preflight = preflight
    ),
    class = "cdid_ipw_point_link_view"
  )
  ipw_raw_fast_ref(pl, "point_link_view")
}

# Build the standalone OR-free compact raw-IPW backend and pipeline registry.
ipw_raw_catalog_row_spec <- function(spec) {
  # The metadata extractor supplies row/operator IDs, not view labels. Match
  # pipeline_metadata_row_views_from_spec and the scalar compact descriptors.
  spec$link_id <- paste0("link:r", spec$row_id)
  spec$nuisance_id <- paste0("nuisance:op", spec$operator_id, ":r", spec$row_id)
  spec
}

build_ipw_raw_backend_registry <- function(
  metadata_obj,
  frame_owner = NULL,
  design = NULL,
  design_link_index = NULL,
  owner_control,
  fitted_control,
  point_control,
  glm_control,
  class_bounds = NULL,
  weighting = "optimal",
  reduce_structural_redundancy = TRUE,
  structural_reduction_stage = "outcome_free",
  aux_spec = NULL,
  retained_pipeline_entries = 0,
  external_work_ops = 0,
  progress = NULL
) {
  validate_owner_comparison_catalog(metadata_obj, frame_owner = frame_owner)
  own_ctrl <- validate_ipw_owner_control(owner_control)
  fit_ctrl <- validate_ipw_fitted_view_control(fitted_control)
  pt_ctrl <- validate_ipw_point_view_control(point_control)
  glm_ctrl <- validate_ipw_glm_control(glm_control)

  cat_rows <- metadata_obj$catalog$rows
  M <- nrow(cat_rows)
  if (M < 1L) {
    cdid_ipw_pipeline_abort(
      "Comparison catalog contains zero qualified links; cannot construct IPW backend.",
      c("cdid_ipw_pipeline_empty_catalog_error", "cdid_gmm_target_estimability_error")
    )
  }

  f_owner <- if (!is.null(frame_owner)) {
    frame_owner
  } else if (!is.null(metadata_obj$frame_ref)) {
    get("frame_owner", envir = metadata_obj$frame_ref, inherits = FALSE)
  } else {
    cdid_ipw_pipeline_abort(
      "`frame_owner` must be supplied either on `metadata_obj` or explicitly.",
      "cdid_ipw_pipeline_input_error"
    )
  }
  validate_ipw_frame_owner(f_owner)
  validate_owner_comparison_catalog(metadata_obj, frame_owner = f_owner)
  p_owner <- get("panel_owner", envir = metadata_obj$panel_ref, inherits = FALSE)
  validate_ipw_frame_owner(f_owner, panel_owner = p_owner)

  if (is.null(design)) design <- metadata_obj$design
  if (is.null(design_link_index)) design_link_index <- match(seq_len(M), design$links$catalog_row)

  retained_pipeline_entries <- ipw_owner_int_bound(retained_pipeline_entries, "retained_pipeline_entries", min_val = 0)
  external_work_ops <- ipw_owner_int_bound(external_work_ops, "external_work_ops", min_val = 0)
  work_prefix <- external_work_ops

  n_a <- as.integer(f_owner$n_a_units)
  p <- as.integer(f_owner$n_covariates)
  # Charge the new snapshot list, names and reference headers, without walking
  # or duplicating the borrowed owner buffers (already charged once below).
  owner_snapshot_entries <- unname(cp04_payload(
    list(owner_snapshot = list(frame = NULL, panel = NULL)))["entries"])
  # Reserve the complete construction lifetime, including the compact tables
  # and the diagnostic records returned to the pipeline. The fixed allowance
  # uses the existing point-view 8192 + 2048 metadata allowance twice (views,
  # fit and point records); matrices/vectors and character metadata are extra.
  # Owners are already in a supplied pipeline envelope; standalone callers
  # still charge them once. Leaf scratch is composed with this inventory below.
  Q_bound <- if (is.null(aux_spec)) 0L else length(aux_spec$descriptors)
  D_bound <- as.double(M) + Q_bound
  C_bound <- as.double(length(f_owner$a_cohorts))
  char_bound <- ipw_fitted_view_char_bytes(
    list(colnames(f_owner$Z_a), as.character(f_owner$a_cohorts),
         paste0("link:r", seq_len(M)), metadata_obj$row_mapping$functional_key), fit_ctrl)$charged_slots
  table_entries_bound <- ipw_owner_safe_product(D_bound,
    20480 + 16 * p * p + 16 * p * C_bound + 256 * (p + C_bound + 1) + 4 * char_bound,
    label = "raw-IPW support/model/link and diagnostic tables")
  construction_base <- max(retained_pipeline_entries, ipw_point_view_owner_entries(f_owner, p_owner))
  construction_retained <- ipw_owner_safe_sum(c(construction_base, table_entries_bound, owner_snapshot_entries),
    min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries),
    "raw-IPW construction owners, pipeline and compact tables")
  numerical_work <- 0
  retained_failure_entries <- 0
  construction_peak <- construction_retained
  C_all <- length(metadata_obj$catalog$spec$cohorts)
  op_rows <- as.double(tabulate(match(cat_rows$operator_id, unique(cat_rows$operator_id))))
  split_bound <- sum(op_rows * (op_rows - 1) / 2)
  mask_registry_entries <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(16, M + 1), C_all),
    own_ctrl$max_retained_entries,
    "invocation-local mask fingerprint registry"
  )
  split_cert_entries <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(8, split_bound, label = "mask split records"),
    ipw_owner_safe_product(M, 11 + C_all, label = "mask-equivalence certificate"),
    ipw_owner_safe_product(128, M + 1, label = "construction support, fingerprint and date fact records")
  ))
  retained_with_masks <- ipw_owner_safe_sum(
    c(retained_pipeline_entries, mask_registry_entries),
    min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries),
    "pipeline envelope plus invocation-local mask fingerprint registry"
  )
  retained_with_masks <- ipw_owner_safe_sum(c(retained_with_masks, split_cert_entries))
  construction_retained <- ipw_owner_safe_sum(c(construction_retained, mask_registry_entries, split_cert_entries),
    min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries),
    "raw-IPW construction tables and mask registry")

  stream_n <- as.double(n_a)
  stream_step <- max(1, min(stream_n, own_ctrl$max_chunk_rows))
  stream_cost <- function(m_c) {
    128 * stream_n * (1 + ceiling(log2(m_c + 1))) +
      128 * ceiling(stream_n / stream_step) * (1 + ceiling(log2(stream_n + 1)))
  }
  scan_work <- 0
  charge_scans <- function(n_streams, m_c) {
    scan_work <<- ipw_owner_safe_sum(c(scan_work, ipw_owner_safe_product(n_streams, stream_cost(m_c))))
    ipw_owner_safe_sum(
      c(external_work_ops, scan_work),
      own_ctrl$max_work_ops,
      "pipeline work plus mask registry scans (including exact collision resolution)"
    )
  }

  descriptor_src <- descriptor_charge <- dim_memo <- NULL
  fit_pre_memo <- pt_pre_memo <- NULL
  mismatch_reason <- ipw_owner_mask_mismatch_message()
  date_memo <- list2env(list(panel_owner = p_owner), parent = emptyenv())
  pair_cache <- new.env(hash = TRUE, parent = emptyenv())
  distinct_pairs <- list()
  pair_keys <- character()

  get_pair_summary <- function(pair_desc) {
    pkey <- paste0(
      if (is.na(pair_desc$s_calendar_index)) sprintf("s:%a", pair_desc$s) else paste0("ds:", pair_desc$s_calendar_index),
      "|",
      if (is.na(pair_desc$t_calendar_index)) sprintf("t:%a", pair_desc$t) else paste0("dt:", pair_desc$t_calendar_index)
    )
    pid <- get0(pkey, envir = pair_cache, inherits = FALSE)
    if (is.null(pid)) {
      charge_scans(1, C_bound)
      ps <- ipw_raw_scan_pair_summary(f_owner, p_owner, pair_desc, stream_step, date_memo)
      pid <- length(distinct_pairs) + 1L
      ps$pair_id <- as.integer(pid)
      distinct_pairs[[pid]] <<- ps
      pair_keys[pid] <<- pkey
      assign(pkey, pid, envir = pair_cache)
    }
    distinct_pairs[[pid]]
  }

  nuisance_buckets <- list()
  distinct_nuisance_views <- list()
  distinct_nuisance_specs <- list()
  distinct_nuisance_keys <- character()
  distinct_operator_ids <- integer()
  distinct_fingerprints <- vector("list", M)
  link_views <- vector("list", M)
  row_specs <- vector("list", M)
  row_pair_ids <- integer(M)
  row_sup_infos <- vector("list", M)
  link_to_nuisance_idx <- integer(M)
  mask_splits <- list()
  n_compared <- n_fp_rejected <- n_exact <- n_exact_mismatch <- n_fp_passes <- 0
  n_support_reused <- 0

  # 1. Catalog row loop: exact per-row preflight, pair/cohort support checks, and
  #    operator mask-branch grouping in catalog row order 1..M.
  for (m in seq_len(M)) {
    pipeline_metadata_row_spec_preflight(metadata_obj, m, own_ctrl)
    pipeline_metadata_row_spec_preflight(metadata_obj, m, metadata_obj$control$owner_control)
    if (is.null(descriptor_src)) {
      descriptor_src <- pipeline_metadata_descriptor_source(
        metadata_obj, own_ctrl, retained_with_masks, external_work_ops, M, memo_entries = 64
      )
      descriptor_charge <- descriptor_src$accounting
      split_cert_entries <- ipw_owner_safe_sum(c(split_cert_entries, descriptor_charge$retained_entries))
      retained_with_masks <- ipw_owner_safe_sum(c(retained_with_masks, descriptor_charge$retained_entries))
      external_work_ops <- ipw_owner_safe_sum(
        c(external_work_ops, descriptor_charge$work_ops_bound),
        own_ctrl$max_work_ops,
        "IPW registry plus descriptor preparation work"
      )
    }
    row_cohorts <- 1 + metadata_obj$catalog$pools$size[cat_rows$pool_id[m]]
    dim_inputs <- list(n_a, p, row_cohorts, own_ctrl)
    if (is.null(dim_memo) || !identical(dim_memo$inputs, dim_inputs)) {
      row_nuisance_preflight <- preflight_ipw_nuisance_view(n_a, p, row_cohorts, own_ctrl)
      row_link_preflight <- preflight_ipw_link_view(n_a, p, row_cohorts, own_ctrl)
      dim_memo <- list(
        inputs = dim_inputs,
        nuisance_preflight = row_nuisance_preflight,
        link_preflight = row_link_preflight,
        peak = max(row_nuisance_preflight$concurrent_entries, row_link_preflight$concurrent_entries),
        retained = max(row_nuisance_preflight$total_retained_entries, row_link_preflight$total_retained_entries)
      )
    }
    ipw_owner_safe_sum(
      c(mask_registry_entries, split_cert_entries, max(retained_pipeline_entries, dim_memo$retained)),
      own_ctrl$max_retained_entries,
      "mask fingerprint registry plus existing row constructor retention"
    )
    ipw_owner_safe_sum(
      c(mask_registry_entries, split_cert_entries, max(retained_pipeline_entries, dim_memo$peak),
        descriptor_charge$scratch_entries),
      own_ctrl$max_matrix_entries,
      "mask fingerprint registry plus row constructor and descriptor scratch envelope"
    )
    spec <- ipw_raw_catalog_row_spec(pipeline_metadata_row_spec_fields(metadata_obj, m, descriptor_src))
    op_id <- as.integer(cat_rows$operator_id[m])
    op_key <- paste0("op_", op_id)
    bucket <- nuisance_buckets[[op_key]]
    matched_idx <- NA_integer_

    pair_desc <- ipw_owner_make_pair_descriptor(p_owner, spec$s, spec$t)
    coh_desc <- ipw_owner_make_cohort_descriptor(f_owner, spec$g, spec$K, require_support = TRUE)
    pair_sum <- get_pair_summary(pair_desc)
    sup_info <- ipw_raw_pair_support_for_cohorts(f_owner, pair_sum, coh_desc, require_support = TRUE)
    row_fp <- sup_info$fingerprint

    if (!is.null(bucket) && length(bucket) > 0L) {
      for (cand_idx in bucket) {
        cand_spec <- distinct_nuisance_specs[[cand_idx]]
        cand_nv <- distinct_nuisance_views[[cand_idx]]
        n_compared <- n_compared + 1
        if (identical(as.double(spec$g), cand_nv$cohorts$g) &&
            identical(as.double(spec$K), cand_nv$cohorts$K)) {
          charge_scans(1, cand_nv$cohorts$m_cohorts)
          n_fp_passes <- n_fp_passes + 1
          if (is.null(distinct_fingerprints[[cand_idx]])) {
            charge_scans(1, cand_nv$cohorts$m_cohorts)
            distinct_fingerprints[[cand_idx]] <- cand_spec$sup_info$fingerprint
            n_fp_passes <- n_fp_passes + 1
          }
          if (!identical(row_fp, distinct_fingerprints[[cand_idx]])) {
            n_fp_rejected <- n_fp_rejected + 1
            mask_splits[[length(mask_splits) + 1L]] <- list(
              row_id = m,
              operator_id = op_id,
              candidate_nuisance_id = cand_nv$nuisance_id,
              reason = mismatch_reason
            )
            next
          }
        }
        n_exact <- n_exact + 1
        charge_scans(2, cand_nv$cohorts$m_cohorts)
        masks_eq <- identical(as.double(spec$g), cand_nv$cohorts$g) &&
          identical(as.double(spec$K), cand_nv$cohorts$K) &&
          ipw_raw_verify_pair_masks_equal(
            f_owner, p_owner, cand_spec$pair_summary, pair_sum, coh_desc, stream_step
          )
        if (isTRUE(masks_eq)) {
          matched_idx <- cand_idx
          if (identical(cand_spec$pair_summary$pair_id, pair_sum$pair_id) ||
              (isTRUE(cand_spec$pair_summary$complete_pair) && isTRUE(pair_sum$complete_pair))) {
            n_support_reused <- n_support_reused + 1
          }
          break
        } else {
          n_exact_mismatch <- n_exact_mismatch + 1
          mask_splits[[length(mask_splits) + 1L]] <- list(
            row_id = m,
            operator_id = op_id,
            candidate_nuisance_id = cand_nv$nuisance_id,
            reason = mismatch_reason
          )
        }
      }
    }

    if (is.na(matched_idx)) {
      charge_scans(3, row_cohorts)
      new_idx <- length(distinct_nuisance_views) + 1L
      branch_num <- if (is.null(bucket)) 1L else (length(bucket) + 1L)
      nuis_key <- paste0("op", op_id, ":mask", branch_num)
      nv_id <- spec$nuisance_id
      nv_rec <- ipw_raw_make_nuisance_view_record(
        nv_id, f_owner, pair_desc, coh_desc, sup_info, own_ctrl, dim_memo$nuisance_preflight
      )
      distinct_nuisance_views[[new_idx]] <- nv_rec
      distinct_nuisance_specs[[new_idx]] <- list(
        nuisance_id = nv_id,
        nuisance_key = nuis_key,
        pair_desc = pair_desc,
        pair_summary = pair_sum,
        coh_desc = coh_desc,
        sup_info = sup_info
      )
      distinct_operator_ids[new_idx] <- op_id
      distinct_nuisance_keys[new_idx] <- nuis_key
      distinct_fingerprints[[new_idx]] <- row_fp
      nuisance_buckets[[op_key]] <- c(bucket, new_idx)
      matched_idx <- new_idx
      n_support_reused <- n_support_reused + 1
    }

    lv_rec <- ipw_raw_make_link_view_record(
      link_id = spec$link_id, functional_key = spec$functional_key,
      frame_owner = f_owner, pair_desc = pair_desc, coh_desc = coh_desc, sup_info = sup_info,
      nv_rec = distinct_nuisance_views[[matched_idx]], own_ctrl = own_ctrl,
      preflight = dim_memo$link_preflight
    )
    row_specs[[m]] <- spec
    link_views[[m]] <- lv_rec
    row_pair_ids[m] <- pair_sum$pair_id
    row_sup_infos[[m]] <- sup_info
    link_to_nuisance_idx[m] <- matched_idx
    if (is.function(progress)) progress(m, M, length(distinct_nuisance_views))
  }
  names(link_views) <- vapply(link_views, `[[`, character(1L), "link_id")
  names(distinct_nuisance_views) <- distinct_nuisance_keys
  descriptor_src <- NULL

  # 2. Validate per-nuisance fitted-view preflight bounds before fitting any model.
  U_cat <- length(distinct_nuisance_views)
  fit_preflights <- vector("list", U_cat)
  for (u in seq_len(U_cat)) {
    nv_u <- distinct_nuisance_views[[u]]
    n_obs_u <- nv_u$support$n_observed_risk_units
    m_u <- nv_u$cohorts$m_cohorts
    if (n_obs_u > fit_ctrl$max_observer_units) {
      cdid_ipw_fitted_view_abort(
        paste0("Observed risk-set unit count (", n_obs_u, ") exceeds `control$max_observer_units` (", fit_ctrl$max_observer_units, ")."),
        c("cdid_ipw_fitted_view_bound_error", "cdid_ipw_owner_bound_error")
      )
    }
    fit_in <- list(n_a, n_obs_u, p, m_u, fit_ctrl)
    if (is.null(fit_pre_memo) || !identical(fit_pre_memo$inputs, fit_in)) {
      fit_pre_memo <- list(
        inputs = fit_in,
        preflight = preflight_ipw_fitted_nuisance_view_dimensions(
          na = n_a, n_obs = n_obs_u, p = p, m = m_u, ctrl = fit_ctrl,
          nuisance_view = nv_u, frame_owner = f_owner, glm_ctrl = glm_ctrl
        )
      )
    }
    fit_preflights[[u]] <- fit_pre_memo$preflight
  }

  # 3. Determine structural class certificate when `weighting == "optimal"`.
  #    Because class/forest/E certification depends only on the validated catalog
  #    design, cohort/pair descriptors, and verified outcome-free nuisance mask
  #    groups (`link_to_nuisance_idx`), under `early_reduction` it runs BEFORE
  #    any outcome endpoints are read so dropped rows never read outcomes.
  wt_mode <- match.arg(weighting, c("optimal", "identity"))
  early_reduction <- identical(wt_mode, "optimal") &&
    isTRUE(reduce_structural_redundancy) &&
    identical(structural_reduction_stage, "outcome_free")

  certify_raw_classes <- function(link_kind, object_refs) {
    b <- validate_ipw_certificate_bounds(class_bounds)
    if (M > b$max_links) {
      cdid_abort("Number of links exceeds `max_links`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
    }
    if (as.double(M) * as.double(n_a) > b$max_mask_entries) {
      cdid_abort("Mask verification work exceeds `max_mask_entries`.", c("cdid_ipw_bound_error", "cdid_ipw_certificate_error"))
    }
    link_g <- vapply(row_specs, function(s) as.double(s$g), numeric(1L))
    link_s <- vapply(row_specs, function(s) as.double(s$s), numeric(1L))
    link_t <- vapply(row_specs, function(s) as.double(s$t), numeric(1L))
    link_fk <- lapply(row_specs, function(s) s$functional_key)
    link_K <- lapply(row_specs, function(s) as.double(s$K))
    checks <- ipw_structural_design_checks(link_g, link_s, link_t, link_K, link_fk, design, design_link_index, M)

    nu_refs <- lapply(distinct_nuisance_views, `[[`, "ref")
    same_mask_fn <- function(ell, rep_ell) {
      if (identical(cat_rows$operator_id[ell], cat_rows$operator_id[rep_ell])) {
        return(identical(link_to_nuisance_idx[ell], link_to_nuisance_idx[rep_ell]))
      }
      u_ell <- link_to_nuisance_idx[ell]
      u_rep <- link_to_nuisance_idx[rep_ell]
      sp_ell <- distinct_nuisance_specs[[u_ell]]
      sp_rep <- distinct_nuisance_specs[[u_rep]]
      if (!identical(sp_ell$sup_info$fingerprint, sp_rep$sup_info$fingerprint)) return(FALSE)
      ipw_raw_verify_pair_masks_equal(
        f_owner, p_owner, sp_rep$pair_summary, sp_ell$pair_summary, sp_rep$coh_desc, stream_step
      )
    }
    ipw_structural_classes_assemble(
      M, link_s, link_t, checks$fk, checks$vs, checks$vt, checks$dli,
      link_to_nuisance_idx, nu_refs,
      link_cohort_g = lapply(row_specs, function(s) as.double(s$g)),
      link_cohort_K = link_K,
      link_glm = replicate(M, glm_ctrl, simplify = FALSE),
      same_mask_fn = same_mask_fn,
      b = b,
      link_kind = link_kind,
      object_refs = object_refs,
      context_ref = f_owner$ref
    )
  }

  class_cert <- NULL
  structural_preflight <- NULL
  if (early_reduction) {
    b_early <- validate_ipw_certificate_bounds(class_bounds)
    T_cal <- as.double(p_owner$n_calendar_dates)
    C_coh <- as.double(length(p_owner$cohorts))
    V_vert <- as.double(nrow(design$vertices))
    path_bound <- min(as.double(M), T_cal - 1, b_early$max_class_dates - 1)
    E_bound <- min(b_early$max_E_entries, ipw_owner_safe_product(M, path_bound, label = "early E entry upper bound"))
    linear <- ipw_owner_safe_sum(c(256 * M, 64 * V_vert, 256), min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries), "early certificate linear plan")
    cert_retained <- ipw_owner_safe_sum(c(5 * E_bound, linear, 4 * M * C_coh), min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries), "early certificate retained")
    work_early <- ipw_owner_safe_sum(c(
      512 * M * (C_coh + 2),
      256 * M * (T_cal + (U_cat + 1) * (C_coh + 1) + ceiling(log2(M + 1))),
      64 * E_bound
    ), own_ctrl$max_work_ops, "early structural mask/design/forest/expansion work")
    structural_preflight <- list(
      E_entries_upper = E_bound,
      fitted_nuisance_entries = sum(vapply(fit_preflights, `[[`, numeric(1L), "retained_entries_per_fit")),
      original_link_view_entries = 0,
      p1_registry_entries = 16 * (M + 1) + C_coh,
      mask_split_entries = 8 * length(mask_splits),
      mask_certificate_entries = M * (11 + C_coh),
      mask_streams_charged = 0,
      mask_evidence = "registry mask-equivalence certificate; only cross-operator comparisons with identical cohorts stream",
      certificate_retained_entries = cert_retained,
      additional_entries = cert_retained,
      work_ops_bound = work_early
    )
    construction_retained <- ipw_owner_safe_sum(c(construction_retained, cert_retained),
      min(own_ctrl$max_retained_entries, own_ctrl$max_matrix_entries),
      "raw-IPW construction tables and structural certificate")
    class_cert <- certify_raw_classes("shared_view_declaration", lapply(link_views, `[[`, "ref"))
    kept_idx <- class_cert$retained_index
    dropped_idx <- class_cert$dropped_index
  } else {
    kept_idx <- seq_len(M)
    dropped_idx <- integer(0L)
  }

  # 4. Auxiliary descriptor intake (1..Q): outcome-free support screening and
  #    matching to existing or new nuisance branches.
  aux_descriptors <- if (!is.null(aux_spec) && is.list(aux_spec$descriptors)) aux_spec$descriptors else list()
  Q <- length(aux_descriptors)
  raw_aux_errors <- vector("list", Q)
  raw_aux_fit_errors <- vector("list", Q)
  raw_aux_counts <- vector("list", Q)
  raw_aux_nuisance_idx <- rep(NA_integer_, Q)
  raw_aux_pair_ids <- rep(NA_integer_, Q)
  raw_aux_sup_infos <- vector("list", Q)
  raw_aux_link_views <- vector("list", Q)

  for (j in seq_len(Q)) {
    d <- aux_descriptors[[j]]
    trial <- tryCatch({
      p_d <- ipw_owner_make_pair_descriptor(p_owner, d$s, d$t)
      c_d <- ipw_owner_make_cohort_descriptor(f_owner, d$g, d$K, require_support = TRUE)
      ps_d <- get_pair_summary(p_d)
      sup_d <- ipw_raw_pair_support_for_cohorts(f_owner, ps_d, c_d, require_support = TRUE)
      list(pair_desc = p_d, coh_desc = c_d, pair_summary = ps_d, sup_info = sup_d)
    }, cdid_ipw_support_error = identity)

    if (inherits(trial, "error")) {
      raw_aux_errors[[j]] <- trial
      p_d <- ipw_owner_make_pair_descriptor(p_owner, d$s, d$t)
      c_d <- ipw_owner_make_cohort_descriptor(f_owner, d$g, d$K, require_support = FALSE)
      ps_d <- get_pair_summary(p_d)
      sup_d <- ipw_raw_pair_support_for_cohorts(f_owner, ps_d, c_d, require_support = FALSE)
      raw_aux_counts[[j]] <- list(
        counts_available = TRUE,
        n_a_units = sup_d$n_a_units,
        n_observed_pair_all = sup_d$n_observed_pair_all,
        n_observed_risk_units = sup_d$n_observed_risk_units,
        by_cohort_a_count = sup_d$by_cohort_a_count,
        by_cohort_pair_count = sup_d$by_cohort_pair_count
      )
      next
    }

    p_d <- trial$pair_desc
    c_d <- trial$coh_desc
    ps_d <- trial$pair_summary
    sup_d <- trial$sup_info
    fp_d <- sup_d$fingerprint
    matched_u <- NA_integer_

    for (u in seq_along(distinct_nuisance_specs)) {
      cand_u <- distinct_nuisance_specs[[u]]
      candidate <- identical(cand_u$coh_desc$g, c_d$g) && identical(cand_u$coh_desc$K, c_d$K) &&
        identical(cand_u$sup_info$fingerprint, fp_d)
      if (candidate) charge_scans(2, c_d$m_cohorts)
      if (candidate && ipw_raw_verify_pair_masks_equal(f_owner, p_owner, cand_u$pair_summary, ps_d, c_d, stream_step)) {
        matched_u <- u
        break
      }
    }
    fkey_j <- paste0("reporting_operator:", gsub(" ", "", paste(d$g, paste(d$K, collapse = ","), sep = ":")))
    if (is.na(matched_u)) {
      fkey_j <- paste0("reporting_operator:", j)
      new_u <- length(distinct_nuisance_views) + 1L
      nv_id <- paste0("reporting_nuisance:", j)
      nv_pre <- preflight_ipw_nuisance_view(n_a, p, c_d$m_cohorts, own_ctrl)
      nv_rec <- ipw_raw_make_nuisance_view_record(nv_id, f_owner, p_d, c_d, sup_d, own_ctrl, nv_pre)
      distinct_nuisance_views[[new_u]] <- nv_rec
      distinct_nuisance_specs[[new_u]] <- list(
        nuisance_id = nv_id,
        nuisance_key = nv_id,
        pair_desc = p_d,
        pair_summary = ps_d,
        coh_desc = c_d,
        sup_info = sup_d,
        is_auxiliary_only = TRUE,
        aux_j = j
      )
      matched_u <- new_u
    }
    lv_pre <- preflight_ipw_link_view(n_a, p, c_d$m_cohorts, own_ctrl)
    raw_aux_link_views[[j]] <- ipw_raw_make_link_view_record(
      link_id = paste0("reporting:", j), functional_key = fkey_j,
      frame_owner = f_owner, pair_desc = p_d, coh_desc = c_d, sup_info = sup_d,
      nv_rec = distinct_nuisance_views[[matched_u]], own_ctrl = own_ctrl, preflight = lv_pre
    )
    raw_aux_nuisance_idx[j] <- matched_u
    raw_aux_pair_ids[j] <- ps_d$pair_id
    raw_aux_sup_infos[[j]] <- sup_d
  }

  # 5. Build deduplicated `(pair_id, cohort)` supports across all active links.
  U_total <- length(distinct_nuisance_specs)
  support_keys <- character()
  support_specs <- list()
  supports <- list()

  get_support_id <- function(pair_id, cohort_pos, coh_desc, sup_info, k_in_coh) {
    skey <- paste0(pair_id, ":", cohort_pos)
    sid <- match(skey, support_keys)
    if (is.na(sid)) {
      sid <- length(support_specs) + 1L
      support_keys[sid] <<- skey
      pair_sum <- distinct_pairs[[pair_id]]
      rate_k <- as.double(sup_info$rates[k_in_coh])
      rate_clip_k <- as.double(ipw_clip_floor(rate_k))
      mass_a_k <- as.double(coh_desc$group_a_mass[k_in_coh])
      q_k <- as.double(coh_desc$q[k_in_coh])
      support_specs[[sid]] <<- list(
        pair = pair_sum$pair,
        pair_id = as.integer(pair_id),
        cohort = as.integer(cohort_pos),
        cohort_value = as.double(coh_desc$cohorts[k_in_coh])
      )
      supports[[sid]] <<- list(
        rate = rate_k,
        rate_clip = rate_clip_k,
        v = 1.0 / rate_clip_k,
        massA = mass_a_k,
        q = q_k,
        n_pair = as.integer(sup_info$by_cohort_pair_count[k_in_coh]),
        pair_mass = as.double(sup_info$by_cohort_pair_mass[k_in_coh])
      )
    }
    as.integer(sid)
  }

  # Map active catalog rows (`rows_to_fit`) and supported auxiliary descriptors
  # onto deduplicated `(pair_id, nuisance_branch)` treatment models.
  rows_to_fit <- if (early_reduction) sort(as.integer(kept_idx)) else seq_len(M)
  nuisance_to_rows <- split(rows_to_fit, link_to_nuisance_idx[rows_to_fit])
  sup_aux_j <- which(!is.na(raw_aux_nuisance_idx))
  nuisance_to_aux <- if (length(sup_aux_j)) split(sup_aux_j, raw_aux_nuisance_idx[sup_aux_j]) else list()

  fitted_views_all <- vector("list", U_total)
  treatment_keys <- character()
  treatment_specs <- list()
  treatments <- list()
  row_treatment_id <- rep(NA_integer_, M)
  aux_treatment_id <- rep(NA_integer_, Q)

  # 6. Single-observer-buffer loop over distinct nuisance branches 1..U_total:
  #    Outcome-free weighted logit + Hessian stage FIRST, then endpoint reads +
  #    point link means/partials on the SAME observer buffer, then immediate release.
  for (u in seq_len(U_total)) {
    sp_u <- distinct_nuisance_specs[[u]]
    nv_u <- distinct_nuisance_views[[u]]
    coh_u <- sp_u$coh_desc
    sup_u <- sp_u$sup_info
    pair_u <- sp_u$pair_summary
    n_obs_u <- sup_u$n_observed_risk_units
    m_u <- coh_u$m_cohorts
    is_aux_only <- isTRUE(sp_u$is_auxiliary_only)

    fit_pre_u <- if (u <= U_cat) fit_preflights[[u]] else {
      if (n_obs_u > fit_ctrl$max_observer_units) {
        err_b <- tryCatch(
          cdid_ipw_fitted_view_abort(
            paste0("Observed risk-set unit count (", n_obs_u, ") exceeds `control$max_observer_units` (", fit_ctrl$max_observer_units, ")."),
            c("cdid_ipw_fitted_view_bound_error", "cdid_ipw_owner_bound_error")
          ),
          error = identity
        )
        raw_aux_fit_errors[[sp_u$aux_j]] <- err_b
        next
      }
      preflight_ipw_fitted_nuisance_view_dimensions(
        na = n_a, n_obs = n_obs_u, p = p, m = m_u, ctrl = fit_ctrl,
        nuisance_view = nv_u, frame_owner = f_owner, glm_ctrl = glm_ctrl
      )
    }

    # Compose the one live observer/solver state with all retained records.
    # Endpoint copies, control-design copies, subsets, triangular solves and
    # expression temporaries can coexist with that state during point work.
    branch_rows <- nuisance_to_rows[[as.character(u)]]
    branch_aux <- nuisance_to_aux[[as.character(u)]]
    branch_pairs <- unique(c(row_pair_ids[branch_rows], raw_aux_pair_ids[branch_aux]))
    point_pre_u <- if (!length(branch_pairs)) list(concurrent_entries = 0, work_ops_bound = 0) else
      preflight_ipw_point_link_view_dimensions(
        n_a, n_obs_u, p, m_u, pt_ctrl, fitted_view = NULL, link_view = NULL)
    leaf_peak <- max(fit_pre_u$concurrent_entries,
      fit_pre_u$observer_buffer_entries + fit_pre_u$worker_workspace_entries +
        point_pre_u$concurrent_entries)
    construction_peak <- max(construction_peak, ipw_owner_safe_sum(
      c(construction_retained, descriptor_charge$retained_entries, retained_failure_entries, leaf_peak),
      min(own_ctrl$max_matrix_entries, fit_ctrl$max_matrix_entries, pt_ctrl$max_matrix_entries),
      "raw-IPW complete construction inventory plus one model buffer and leaf copies"))
    numerical_work <- ipw_owner_safe_sum(c(numerical_work, fit_pre_u$work_ops_bound,
      ipw_owner_safe_product(length(branch_pairs), point_pre_u$work_ops_bound)))
    ipw_owner_safe_sum(c(external_work_ops, scan_work, numerical_work), own_ctrl$max_work_ops,
      "raw-IPW cumulative construction scans, weighted fits, Hessians and point work")

    obs_buf <- ipw_raw_collect_observer_buffer(f_owner, p_owner, pair_u, coh_u, n_obs_u, stream_step)
    sup_worker <- list(
      n = n_a,
      g = coh_u$g,
      K = coh_u$K,
      cohorts = coh_u$cohorts,
      cohort_labels = coh_u$cohort_labels,
      m_cohorts = m_u,
      group_count = as.integer(sup_u$by_cohort_a_count),
      group_a_count = as.integer(sup_u$by_cohort_a_count),
      group_mass = as.double(sup_u$by_cohort_a_mass),
      group_a_mass = as.double(sup_u$by_cohort_a_mass),
      q = as.double(sup_u$q),
      a_weight_total = as.double(f_owner$sum_w_a)
    )

    diag_env <- new.env(parent = emptyenv())
    fit_attempt <- tryCatch({
      st_u <- ipw_nuisance_observer_fit_worker(
        sup = sup_worker,
        risk_a_idx = obs_buf$risk_a_idx,
        risk_gi = obs_buf$risk_gi,
        risk_w = obs_buf$risk_w,
        Z_use = obs_buf$Z_use,
        ctrl = glm_ctrl,
        diag_env = diag_env
      )
      hs_u <- tryCatch(
        ipw_nuisance_hessian_stage(st_u, coh_u$cohort_labels),
        error = function(e) {
          diag_env$error_diagnostics <- list(
            stage = "hessian_stage",
            converged = isTRUE(st_u$fit_info$converged),
            boundary = isTRUE(st_u$fit_info$boundary),
            iter = st_u$fit_info$iter,
            rank = st_u$fit_info$rank,
            failed_coefficients = st_u$beta,
            failed_fitted_values = st_u$prob_use,
            failed_odds_values = st_u$odds_use
          )
          stop(e)
        }
      )
      if (any(!is.finite(hs_u$A_beta_r))) {
        diag_env$error_diagnostics <- list(
          stage = "score_rate_jacobian_check",
          failed_coefficients = st_u$beta,
          failed_fitted_values = st_u$prob_use,
          failed_A_beta_r = hs_u$A_beta_r
        )
        cdid_abort(
          "Nonfinite propensity or outcome partial derivatives.",
          c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
        )
      }
      list(st = st_u, hs = hs_u)
    }, error = function(cnd) {
      cnd$error_diagnostics <- if (exists("error_diagnostics", envir = diag_env, inherits = FALSE)) {
        diag_env$error_diagnostics
      } else {
        list(stage = "unknown_fit_failure")
      }
      cnd
    })

    if (inherits(fit_attempt, "error")) {
      st_u <- hs_u <- obs_buf <- NULL
      if (is_aux_only) {
        retained_failure_entries <- ipw_owner_safe_sum(c(retained_failure_entries,
          unname(cp04_payload(fit_attempt)["entries"])))
        aux_rows_u <- nuisance_to_aux[[as.character(u)]]
        for (aj in aux_rows_u) raw_aux_fit_errors[[aj]] <- fit_attempt
        next
      }
      stop(fit_attempt)
    }

    st_u <- fit_attempt$st
    hs_u <- fit_attempt$hs
    fit_attempt <- NULL

    eff_chunk <- fit_pre_u$effective_chunk_rows
    st_summary <- list(
      n_treated_obs = as.integer(sum(st_u$is_t)),
      n_control_obs = as.integer(sum(st_u$is_c)),
      pair_count = st_u$pair_count,
      pair_mass = st_u$pair_mass,
      rates = st_u$rates,
      beta = st_u$beta,
      sum_a = as.double(st_u$sum_a),
      mean_a = as.double(st_u$mean_a),
      sum_b = as.double(st_u$sum_b),
      mean_b = as.double(st_u$mean_b),
      fit_info = st_u$fit_info,
      weight_concentration = list(
        propensity_fit_weights = ipw_fitted_diagnostic_summary(st_u$fit_weights, eff_chunk, TRUE),
        target_outcome_weights = ipw_fitted_diagnostic_summary(st_u$a_obs, eff_chunk, TRUE),
        pooled_control_outcome_ipw_weights = ipw_fitted_diagnostic_summary(st_u$b_obs, eff_chunk, TRUE)
      ),
      probability_range = ipw_fitted_diagnostic_summary(st_u$prob_use, eff_chunk, FALSE)
    )
    fv_rec <- ipw_raw_make_fitted_view_record(nv_u, f_owner, st_summary, hs_u, glm_ctrl, fit_ctrl, fit_pre_u)
    fitted_views_all[[u]] <- fv_rec

    # Fit each distinct pair on this nuisance branch using the ALREADY-LIVE
    # observer buffer `st_u` and `hs_u` before releasing `st_u`.
    cat_rows_u <- nuisance_to_rows[[as.character(u)]]
    aux_rows_u <- nuisance_to_aux[[as.character(u)]]

    fit_pair_on_branch <- function(pid, sup_info_p, is_aux_target = FALSE) {
      tkey <- paste0(pid, ":", u)
      tid <- match(tkey, treatment_keys)
      if (!is.na(tid)) return(tid)

      if (n_obs_u > pt_ctrl$max_observer_units) {
        cdid_ipw_point_view_abort(
          paste0("Observed risk-set unit count (", n_obs_u, ") exceeds `control$max_observer_units` (", pt_ctrl$max_observer_units, ")."),
          c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error")
        )
      }

      ep_pos <- if (identical(as.integer(pid), as.integer(pair_u$pair_id))) {
        list(obs_s = obs_buf$obs_s, obs_t = obs_buf$obs_t)
      } else {
        ipw_raw_collect_pair_endpoints(
          p_owner, distinct_pairs[[pid]], coh_u, obs_buf$risk_a_idx, f_owner$a_cohort_index, stream_step
        )
      }
      ep <- ipw_owner_read_y_positions(p_owner, ep_pos$obs_s, ep_pos$obs_t)
      risk_dy <- as.double(ep$y_right_obs) - as.double(ep$y_left_obs)
      if (anyNA(risk_dy) || any(!is.finite(risk_dy))) {
        cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
      }

      means <- ipw_link_means(st_u, risk_dy)
      part <- ipw_link_partials(st_u, hs_u, risk_dy, means, coh_u$q, coh_u$cohort_labels)

      sids <- integer(m_u)
      for (k in seq_len(m_u)) {
        sids[k] <- get_support_id(pid, coh_u$a_cohort_map[k], coh_u, sup_info_p, k)
      }
      tid <- length(treatments) + 1L
      treatment_keys[tid] <<- tkey
      treatment_specs[[tid]] <<- list(
        pair = distinct_pairs[[pid]]$pair,
        pair_id = as.integer(pid),
        cohorts = as.integer(coh_u$a_cohort_map),
        cohort_values = as.double(coh_u$cohorts),
        cohort_labels = coh_u$cohort_labels,
        support_ids = as.integer(sids),
        n_observers = as.integer(n_obs_u),
        nuisance_idx = as.integer(u)
      )
      treatments[[tid]] <<- list(
        coefficients = fv_rec$coefficients,
        H = hs_u$H,
        H_chol = hs_u$H_chol,
        A_beta_r = hs_u$A_beta_r,
        rates = fv_rec$masses$rates,
        rates_clip = ipw_clip_floor(as.double(fv_rec$masses$rates)),
        q = fv_rec$masses$q,
        means = list(
          treated_mean = as.double(means$Tmean),
          control_mean = as.double(means$Cmean),
          estimate = as.double(means$mu),
          sum_a = fv_rec$means$sum_a,
          mean_a = fv_rec$means$mean_a,
          sum_b = fv_rec$means$sum_b,
          mean_b = fv_rec$means$mean_b
        ),
        partials = list(
          G_beta = part$G_beta,
          G_r = part$G_r,
          A_beta_r = hs_u$A_beta_r
        ),
        loadings = list(
          loading_beta = part$loading_beta,
          loading_rate = part$loading_rate,
          rate_coef_s0 = part$rate_coef_s0,
          rate_coef_s1 = part$rate_coef_s1
        )
      )
      as.integer(tid)
    }

    if (length(cat_rows_u) > 0L) {
      for (m in cat_rows_u) {
        row_treatment_id[m] <- fit_pair_on_branch(row_pair_ids[m], row_sup_infos[[m]], FALSE)
      }
    }
    if (length(aux_rows_u) > 0L) {
      for (aj in aux_rows_u) {
        if (!is.null(raw_aux_fit_errors[[aj]])) next
        tid_try <- tryCatch(
          fit_pair_on_branch(raw_aux_pair_ids[aj], raw_aux_sup_infos[[aj]], TRUE),
          error = identity
        )
        if (inherits(tid_try, "error")) {
          raw_aux_fit_errors[[aj]] <- tid_try
        } else {
          aux_treatment_id[aj] <- tid_try
        }
      }
    }

    # Immediately scrub and release the observer buffer before the next nuisance branch!
    st_u <- hs_u <- obs_buf <- NULL
  }

  # 7. Assemble `point_links_all`, backend `links`, and auxiliary point links.
  point_links_all <- vector("list", M)
  backend_links <- list()

  for (m in rows_to_fit) {
    tid <- row_treatment_id[m]
    t_spec <- treatment_specs[[tid]]
    t_fit <- treatments[[tid]]
    lv_m <- link_views[[m]]
    fv_m <- fitted_views_all[[link_to_nuisance_idx[m]]]
    pt_in <- list(n_a, t_spec$n_observers, p, length(t_spec$cohorts), pt_ctrl)
    if (is.null(pt_pre_memo) || !identical(pt_pre_memo$inputs, pt_in)) {
      pt_pre_memo <- list(
        inputs = pt_in,
        preflight = preflight_ipw_point_link_view_dimensions(
          na = n_a, n_obs = t_spec$n_observers, p = p, m = length(t_spec$cohorts),
          ctrl = pt_ctrl, fitted_view = fv_m, link_view = lv_m
        )
      )
    }
    b_idx <- length(backend_links) + 1L
    backend_links[[b_idx]] <- list(
      link_id = lv_m$link_id,
      source_functional_key = lv_m$functional_key,
      pair_id = t_spec$pair_id,
      cohort = t_spec$cohorts[1L],
      support_id = t_spec$support_ids[1L],
      control_support_ids = t_spec$support_ids[-1L],
      model = as.integer(tid),
      regression = as.integer(tid),
      mu = as.double(t_fit$means$estimate)
    )
    point_links_all[[m]] <- ipw_raw_make_point_link_record(
      lv_m, fv_m, t_fit, b_idx, pt_ctrl, pt_pre_memo$preflight
    )
  }

  if (early_reduction && length(dropped_idx) > 0L) {
    for (m in dropped_idx) {
      lv_m <- link_views[[m]]
      fv_m <- fitted_views_all[[link_to_nuisance_idx[m]]]
      n_e <- length(class_cert$E$i)
      first <- ipw_owner_lower_bound(class_cert$E$i, 1, n_e, m)
      last <- ipw_owner_upper_bound(class_cert$E$i, first, n_e, m)
      hit <- if (last >= first) seq.int(first, last) else integer(0L)
      point_links_all[[m]] <- structure(
        list(
          link_id = lv_m$link_id,
          catalog_row = as.integer(m),
          link_view_ref = lv_m$ref,
          nuisance_view_ref = fv_m$nuisance_view_ref,
          fitted_nuisance_ref = fv_m$ref,
          expansion = list(retained_rows = class_cert$E$j[hit], signs = class_cert$E$x[hit]),
          point_fit_performed = FALSE,
          outcomes_read = FALSE,
          status = "structurally_reduced_outcome_free",
          scope = "Original comparison recovered through the certified expansion E (Prop. 5/Cor. 6); never fitted, replayed or read."
        ),
        class = "cdid_ipw_reduced_link_record"
      )
    }
  }
  names(point_links_all) <- names(link_views)

  if (!early_reduction && identical(wt_mode, "optimal")) {
    class_cert <- certify_raw_classes("shared_view_fit", lapply(point_links_all, `[[`, "ref"))
    if (isTRUE(reduce_structural_redundancy)) {
      kept_idx <- class_cert$retained_index
      dropped_idx <- class_cert$dropped_index
    } else {
      kept_idx <- seq_len(M)
      dropped_idx <- integer(0L)
    }
  }

  raw_aux_index <- rep(NA_integer_, Q)
  raw_aux_points <- vector("list", Q)
  for (j in seq_len(Q)) {
    tid <- aux_treatment_id[j]
    if (is.na(tid)) next
    t_spec <- treatment_specs[[tid]]
    t_fit <- treatments[[tid]]
    lv_j <- raw_aux_link_views[[j]]
    fv_j <- fitted_views_all[[raw_aux_nuisance_idx[j]]]
    pt_pre_j <- preflight_ipw_point_link_view_dimensions(
      na = n_a, n_obs = t_spec$n_observers, p = p, m = length(t_spec$cohorts),
      ctrl = pt_ctrl, fitted_view = fv_j, link_view = lv_j
    )
    b_idx <- length(backend_links) + 1L
    backend_links[[b_idx]] <- list(
      link_id = lv_j$link_id,
      source_functional_key = lv_j$functional_key,
      pair_id = t_spec$pair_id,
      cohort = t_spec$cohorts[1L],
      support_id = t_spec$support_ids[1L],
      control_support_ids = t_spec$support_ids[-1L],
      model = as.integer(tid),
      regression = as.integer(tid),
      mu = as.double(t_fit$means$estimate)
    )
    raw_aux_index[j] <- as.integer(b_idx)
    raw_aux_points[[j]] <- ipw_raw_make_point_link_record(
      lv_j, fv_j, t_fit, b_idx, pt_ctrl, pt_pre_j
    )
  }

  fitted_views_cat <- fitted_views_all[seq_len(U_cat)]
  names(fitted_views_cat) <- distinct_nuisance_keys
  nuisance_views_cat <- distinct_nuisance_views[seq_len(U_cat)]
  names(nuisance_views_cat) <- distinct_nuisance_keys

  mask_certificate <- ipw_seal(
    list(
      kind = "cdid_ipw_registry_mask_equivalence",
      frame_ref = f_owner$ref,
      panel_ref = f_owner$panel_ref,
      metadata_ref = metadata_obj$ref,
      n_rows = as.double(M),
      row_operator = as.integer(cat_rows$operator_id),
      row_view = as.integer(link_to_nuisance_idx),
      row_link_refs = unname(lapply(link_views, `[[`, "ref")),
      view_operator = as.integer(distinct_operator_ids),
      view_g = unname(vapply(nuisance_views_cat, function(v) as.double(v$cohorts$g), numeric(1L))),
      view_K = unname(lapply(nuisance_views_cat, function(v) as.double(v$cohorts$K))),
      view_nuisance_refs = unname(lapply(nuisance_views_cat, `[[`, "ref")),
      view_fit_refs = unname(lapply(fitted_views_cat, `[[`, "ref")),
      separations = list(fingerprint_differences = n_fp_rejected, exact_streamed_mismatches = n_exact_mismatch)
    ),
    "mask_certificate"
  )

  # 8. Assemble and seal the compact `cdid_ipw_raw_backend` object once.
  pairs_compact <- lapply(distinct_pairs, function(ps) {
    list(
      pair_id = ps$pair_id,
      pair = ps$pair,
      complete_pair = ps$complete_pair,
      n_observed_pair_all = ps$n_observed_pair_all
    )
  })
  owned_entries <- unname(cp04_payload(list(
    pairs_compact, support_specs, supports, treatment_specs, treatments, backend_links
  ))["entries"]) + 256 + owner_snapshot_entries
  borrowed_entries <- ipw_point_view_owner_entries(f_owner, p_owner)
  work_charged <- ipw_owner_safe_sum(c(scan_work, descriptor_charge$work_ops_bound, numerical_work))
  composed_work <- ipw_owner_safe_sum(
    c(work_prefix, work_charged),
    own_ctrl$max_work_ops,
    "compact raw-IPW backend construction work"
  )

  backend_ctrl <- ipw_owner_control(
    max_chunk_rows = min(own_ctrl$max_chunk_rows, pt_ctrl$max_chunk_rows),
    max_retained_entries = min(own_ctrl$max_retained_entries, pt_ctrl$max_retained_entries),
    max_matrix_entries = min(own_ctrl$max_matrix_entries, pt_ctrl$max_matrix_entries),
    max_work_ops = min(own_ctrl$max_work_ops, pt_ctrl$max_work_ops),
    max_output_entries = min(own_ctrl$max_output_entries, pt_ctrl$max_output_entries)
  )

  raw_backend <- structure(
    list(
      kind = "cdid_ipw_raw_backend",
      frame_ref = f_owner$ref,
      panel_ref = f_owner$panel_ref,
      owner_snapshot = list(frame = f_owner, panel = p_owner),
      reader_ref = NULL,
      reference_ref = NULL,
      selection_backend = NULL,
      design_ref = NULL,
      layout = NULL,
      complete_support = list(dates = integer(0L)),
      n_a_units = n_a,
      n_total_units = as.integer(f_owner$n_total_units + f_owner$n_frame_only_units),
      p = p,
      covariate_names = colnames(f_owner$Z_a),
      control = backend_ctrl,
      owner_control = own_ctrl,
      fitted_control = fit_ctrl,
      point_control = pt_ctrl,
      glm_control = glm_ctrl,
      pairs = pairs_compact,
      support_specs = support_specs,
      supports = supports,
      treatment_specs = treatment_specs,
      treatments = treatments,
      links = backend_links,
      accounting = list(
        owned_entries = owned_entries,
        borrowed_entries = borrowed_entries,
        owner_snapshot_entries = owner_snapshot_entries,
        construction_peak_entries = construction_peak,
        construction_table_entries_bound = table_entries_bound,
        simultaneous_scratch_bound = construction_peak - construction_retained,
        work_ops_bound = composed_work,
        full_n_masks_retained = FALSE,
        observer_arrays_retained = FALSE
      )
    ),
    class = "cdid_ipw_raw_backend"
  )
  raw_backend <- ipw_owner_seal(raw_backend, "ipw_raw_backend")

  list(
    raw_ipw = raw_backend,
    nuisance_views = nuisance_views_cat,
    fitted_views = fitted_views_cat,
    nuisance_keys = distinct_nuisance_keys,
    operator_ids = distinct_operator_ids,
    row_specs = row_specs,
    link_views = link_views,
    point_links = point_links_all,
    point_links_all = point_links_all,
    link_to_nuisance_idx = link_to_nuisance_idx,
    link_nuisance_keys = distinct_nuisance_keys[link_to_nuisance_idx],
    mask_splits = mask_splits,
    mask_certificate = mask_certificate,
    class_cert = class_cert,
    kept_idx = kept_idx,
    dropped_idx = dropped_idx,
    structural_preflight = structural_preflight,
    raw_aux_index = raw_aux_index,
    raw_aux_errors = raw_aux_errors,
    raw_aux_fit_errors = raw_aux_fit_errors,
    raw_aux_counts = raw_aux_counts,
    raw_aux_points = raw_aux_points,
    raw_aux_fitted_views = fitted_views_all,
    raw_aux_nuisance_views = distinct_nuisance_views,
    construction_work = list(
      n_rows = M,
      n_distinct_nuisance_views = U_cat,
      candidate_comparisons = n_compared,
      fingerprint_rejections = n_fp_rejected,
      exact_candidate_attempts = n_exact,
      exact_mask_mismatches = n_exact_mismatch,
      fingerprint_passes = n_fp_passes,
      link_view_constructions = n_exact + U_cat,
      endpoints_read = 0L,
      endpoints_read_scope = "mask_screening_only",
      work_ops_charged = work_charged,
      registry_entries_reserved = mask_registry_entries + split_cert_entries,
      split_records_bound = split_bound,
      link_support_scans_reused = n_support_reused,
      complete_pair_fingerprints_reused = 0
    ),
    raw_work_accounting = list(
      mode = "cumulative_includes_prefix",
      prefix = as.double(work_prefix),
      work_charged = as.double(work_charged),
      composed_total = as.double(composed_work)
    )
  )
}

# Dimension-only preflights for the retained raw records and numerical leaves.
# Whole-construction owners/tables and the live leaf copies are composed above.
preflight_ipw_fitted_nuisance_view_dimensions <- function(na, n_obs, p, m, ctrl, nuisance_view, frame_owner, glm_ctrl) {
  max_ch <- ctrl$max_chunk_rows
  max_obs <- ctrl$max_observer_units
  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  B_eff <- min(na, max_ch)
  n_chunks <- ceiling(na / B_eff)

  char_meta <- ipw_fitted_view_char_bytes(
    list(
      nuisance_view$nuisance_id,
      nuisance_view$cohorts$cohort_labels,
      colnames(frame_owner$Z_a)
    ),
    ctrl
  )

  ret_matrix_entries <- ipw_owner_safe_sum(
    c(3 * p * p, p * m, nuisance_view$preflight$retained_matrix_entries, 64),
    min(max_ret, max_out, max_me),
    "fitted nuisance view retained matrix entries"
  )
  ret_vector_entries <- ipw_owner_safe_sum(
    c(4 * p, 14 * m, 2176, nuisance_view$preflight$retained_vector_entries, char_meta$charged_slots),
    min(max_ret, max_out, max_me),
    "fitted nuisance view retained vector entries"
  )
  ret_total <- ipw_owner_safe_sum(
    c(ret_matrix_entries, ret_vector_entries),
    min(max_ret, max_out, max_me),
    "fitted nuisance view total retained entries"
  )

  q_pre <- preflight_ipw_range_query(
    b_len = B_eff,
    p_covariates = p,
    m_cohorts = m,
    control = ipw_fitted_view_to_owner_control(ctrl),
    read_outcomes = FALSE,
    include_design = TRUE,
    enclose_differences = FALSE
  )
  obs_buf_entries <- ipw_owner_safe_sum(
    c(n_obs * p, 6 * n_obs, 32),
    max_me,
    "observer buffer entries"
  )
  worker_temp_entries <- ipw_owner_safe_sum(
    c(6 * n_obs * p, 24 * n_obs, 6 * p * p, 12 * p, 4 * p * m, 8 * m, 256),
    max_me,
    "weighted-logit and Hessian workspace entries"
  )
  diag_scratch_entries <- ipw_owner_safe_sum(c(2 * B_eff, 64), max_me, "nuisance diagnostic summary scratch")

  streaming_phase_entries <- ipw_owner_safe_sum(
    c(ret_total, obs_buf_entries, q_pre$concurrent_entries, char_meta$charged_slots),
    max_bound = max_me,
    label = "simultaneous streaming-phase entries (observer buffer + 1 reader chunk + retained)"
  )
  solver_phase_entries <- ipw_owner_safe_sum(
    c(ret_total, obs_buf_entries, worker_temp_entries, diag_scratch_entries),
    max_bound = max_me,
    label = "simultaneous solver-phase entries (observer buffer + fitting/QR/derivative workspace + retained)"
  )
  peak_concurrent <- max(streaming_phase_entries, solver_phase_entries)

  stream_work <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, q_pre$work_ops_bound),
      ipw_owner_safe_product(n_obs, p + 8),
      char_meta$work_ops_bound
    ),
    max_wo,
    "fitted nuisance view streaming work"
  )
  fit_work <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(
        glm_ctrl$maxit + 2,
        ipw_owner_safe_sum(c(4 * n_obs * p * p, 16 * n_obs * p, 4 * p * p * p))
      ),
      4 * n_obs * p * p,
      4 * n_obs * p * m,
      4 * p * p * p,
      4 * p * p * m,
      16 * n_obs + 64 * ceiling(n_obs / B_eff) + 128
    ),
    max_wo,
    "fitted nuisance view solver and Hessian work"
  )
  total_work <- ipw_owner_safe_sum(c(stream_work, fit_work), max_wo, "fitted nuisance view total work")

  list(
    n_a_units = na,
    n_observed_risk_units = n_obs,
    p_covariates = p,
    m_cohorts = m,
    effective_chunk_rows = B_eff,
    n_chunks = n_chunks,
    retained_matrix_entries_per_fit = ret_matrix_entries,
    retained_vector_entries_per_fit = ret_vector_entries,
    retained_entries_per_fit = ret_total,
    observer_buffer_entries = obs_buf_entries,
    worker_workspace_entries = worker_temp_entries,
    concurrent_entries = peak_concurrent,
    work_ops_bound = total_work,
    total_work_ops_bound = total_work
  )
}

preflight_ipw_point_link_view_dimensions <- function(na, n_obs, p, m, ctrl, fitted_view, link_view) {
  max_ch <- ctrl$max_chunk_rows
  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  if (n_obs > ctrl$max_observer_units) {
    cdid_ipw_point_view_abort("Observed risk-set unit count exceeds `control$max_observer_units`.",
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error"))
  }

  B_eff <- min(na, max_ch)
  n_chunks <- ceiling(na / B_eff)

  ret_matrix_entries <- ipw_owner_safe_sum(c(3 * p * p, p * m,
    if (is.null(link_view)) 0 else link_view$preflight$retained_matrix_entries, 64),
    min(max_ret, max_out, max_me), "point link retained matrix entries")
  ret_vector_entries <- ipw_owner_safe_sum(
    c(5 * p, 18 * m, 10240,
      if (is.null(link_view)) 0 else link_view$preflight$retained_vector_entries),
    min(max_ret, max_out, max_me),
    "point link retained vector entries"
  )
  ret_total <- ipw_owner_safe_sum(
    c(ret_matrix_entries, ret_vector_entries),
    min(max_ret, max_out, max_me),
    "point link total retained entries"
  )
  obs_buf_entries <- ipw_owner_safe_sum(c(5 * n_obs * p, 68 * n_obs,
    4 * p * p, 2 * p * m, 8 * p, 12 * m, 256), max_me,
    "point link observer buffer, endpoint copies and partial workspace")
  total_work <- ipw_owner_safe_sum(
    c(32 * na + 32 * n_obs * (p + m + 4) + 16 * (p * p + p * m) + 512),
    max_wo,
    "point link total work"
  )

  list(
    n_a_units = na,
    n_observed_risk_units = n_obs,
    p_covariates = p,
    m_cohorts = m,
    effective_chunk_rows = B_eff,
    n_chunks = n_chunks,
    retained_matrix_entries_per_link = ret_matrix_entries,
    retained_vector_entries_per_link = ret_vector_entries,
    retained_entries_per_link = ret_total,
    observer_buffer_entries = obs_buf_entries,
    concurrent_entries = obs_buf_entries + ret_total,
    work_ops_bound = total_work,
    total_work_ops_bound = total_work
  )
}

# ---------------------------------------------------------------------------
# Batched Complete Original-Unit Raw-IPW Influence Function Evaluator
# ---------------------------------------------------------------------------

ipw_raw_if_envelope <- function(B, S, M_models, L, p, m_max, details = FALSE, P = S) {
  B <- as.double(B)
  S <- as.double(S)
  M_models <- as.double(M_models)
  L <- as.double(L)
  p <- as.double(p)
  m_max <- as.double(m_max)

  out_mult <- if (isTRUE(details)) 5 else 1
  output <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(out_mult, ipw_owner_safe_product(B, L)),
    6 * B, 2 * L, 256
  ), label = "raw-IPW IF chunk output")

  scratch <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(B, if (isTRUE(details)) S else 0),
    ipw_owner_safe_product(B, 3 * P + 6 * p + 96),
    ipw_owner_safe_product(out_mult, ipw_owner_safe_product(B, max(1, M_models))),
    16 * (S + M_models + L + p + 16)
  ), label = "raw-IPW IF chunk scratch")

  work <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(B, 16 * S + 24 * M_models * (p + m_max + 4) + 8 * L),
    256 * (S + M_models + L + 1)
  ), label = "raw-IPW IF chunk work")

  list(output = output, scratch = scratch, work = work)
}

ipw_raw_if_budget <- function(backend, B, links, details = FALSE) {
  validate_ipw_raw_backend(backend)
  ctrl <- backend$control
  n_a <- backend$n_a_units
  n_links_total <- length(backend$links)

  if (!is.numeric(B) || is.object(B) || !is.null(dim(B)) || length(B) != 1L || is.na(B) || !is.finite(B) ||
      B != floor(B) || B < 1 || B > n_a || B > ctrl$max_chunk_rows) {
    ipw_raw_backend_abort("Invalid IF row chunk size `B`.", "input")
  }
  valid_links <- is.numeric(links) && !is.object(links) && is.null(dim(links)) &&
    length(links) > 0L && all(is.finite(links)) && all(links == floor(links)) &&
    all(links >= 1 & links <= n_links_total) && !anyDuplicated(links)
  if (!valid_links || !identical(details, isTRUE(details))) {
    ipw_raw_backend_abort("Invalid explicit IF link batch or details flag.", "input")
  }

  lk <- backend$links[links]
  tids <- unique(vapply(lk, `[[`, integer(1L), "model"))
  sids <- unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids"), use.names = FALSE))
  pids <- unique(vapply(backend$treatment_specs[tids], `[[`, integer(1L), "pair_id"))
  m_max <- max(1L, vapply(backend$treatment_specs[tids], function(s) length(s$cohorts), integer(1L)))

  env <- ipw_raw_if_envelope(
    B = B,
    S = length(sids),
    M_models = length(tids),
    L = length(links),
    p = backend$p,
    m_max = m_max,
    details = details,
    P = length(pids)
  )
  cp04_budget(
    ctrl,
    c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    c(env$output, env$scratch),
    env$work,
    env$output,
    "raw-IPW IF chunk evaluation"
  )
  c(env, list(pids = as.integer(pids), sids = as.integer(sids), tids = as.integer(tids)))
}

ipw_raw_if_plan <- function(backend, links, B) {
  pre <- ipw_raw_if_budget(backend, B, links, details = FALSE)
  lk <- backend$links[links]
  ids <- vapply(lk, `[[`, character(1L), "link_id")
  model_pos <- match(vapply(lk, `[[`, integer(1L), "model"), pre$tids)
  plan <- list(
    kind = "ipw_raw_if_plan",
    backend_ref = backend$ref,
    links = as.integer(links),
    B = as.integer(B),
    pids = pre$pids,
    sids = pre$sids,
    tids = pre$tids,
    model_pos = as.integer(model_pos),
    link_ids = ids,
    pre = pre,
    entries = 0
  )
  plan$entries <- unname(cp04_payload(plan)["entries"])
  ipw_owner_seal(plan, "ipw_raw_if_plan")
}

# Query pair positions on `rows` (in `1..n_a_units`), supporting both contiguous
# chunks `[a_s, a_e]` and arbitrary distinct index subsets.
ipw_raw_query_pair_rows <- function(panel_owner, pair_rec, rows) {
  B <- length(rows)
  pair_desc <- pair_rec$pair
  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  obs <- logical(B)
  pos_s <- integer(B)
  pos_t <- integer(B)
  if (is.na(ds) || is.na(dt) ||
      panel_owner$postings$date_count[ds] < 1L ||
      panel_owner$postings$date_count[dt] < 1L) {
    return(list(obs = obs, pos_s = pos_s, pos_t = pos_t))
  }
  a_s <- as.integer(rows[1L])
  a_e <- as.integer(rows[B])
  is_contig <- (a_e - a_s + 1L == B) && identical(as.integer(rows), seq.int(a_s, a_e))
  if (is_contig) {
    hit <- if (isTRUE(pair_rec$complete_pair)) {
      ipw_owner_query_complete_pair_mask(panel_owner, ds, dt, a_s, a_e)
    } else {
      ipw_owner_query_pair_mask(
        panel_owner, pair_desc$s, pair_desc$t,
        a_start = a_s, a_end = a_e, ds = ds, dt = dt
      )
    }
    if (length(hit$a_index) > 0L) {
      loc <- hit$a_index - a_s + 1L
      obs[loc] <- TRUE
      pos_s[loc] <- hit$obs_s
      pos_t[loc] <- hit$obs_t
    }
    return(list(obs = obs, pos_s = pos_s, pos_t = pos_t))
  }
  post <- panel_owner$postings
  pa <- post$post_a_index
  r_dbl <- as.double(rows)
  hs <- ipw_owner_match_sorted_double(r_dbl, pa, post$date_start[ds], post$date_end[ds])
  ht <- ipw_owner_match_sorted_double(r_dbl, pa, post$date_start[dt], post$date_end[dt])
  ok <- which(!is.na(hs) & !is.na(ht))
  if (length(ok) > 0L) {
    obs[ok] <- TRUE
    pos_s[ok] <- post$post_obs_index[hs[ok]]
    pos_t[ok] <- post$post_obs_index[ht[ok]]
  }
  list(obs = obs, pos_s = pos_s, pos_t = pos_t)
}

evaluate_ipw_raw_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE,
                                      plan = NULL, include_components = FALSE) {
  validate_ipw_raw_backend(backend)
  n_a <- backend$n_a_units
  valid_index <- function(x, n) {
    is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
      length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) &&
      all(x >= 1 & x <= n) && !anyDuplicated(x)
  }
  want_details <- isTRUE(details) || isTRUE(include_components)
  if (!valid_index(rows, n_a) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) ||
      !gmm_scalar_logical(details) || !gmm_scalar_logical(include_components) ||
      !is.character(scale) || is.object(scale) || !is.null(dim(scale)) ||
      length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core")) {
    ipw_raw_backend_abort("Invalid explicit IF row/link batch, scale, or details flag.", "input")
  }

  B <- length(rows)
  rows_i <- as.integer(rows)
  links_i <- as.integer(links)

  if (!is.null(plan) && !want_details) {
    plan_ref <- if (typeof(plan) == "list" && !is.object(plan)) .subset2(plan, "ref") else NULL
    if (!is.environment(plan_ref) || !exists("ipw_raw_if_plan", plan_ref, inherits = FALSE) ||
        bindingIsActive("ipw_raw_if_plan", plan_ref) ||
        !identical(plan$kind, "ipw_raw_if_plan") ||
        !ipw_owner_seal_ok(plan, "ipw_raw_if_plan") ||
        !identical(plan$backend_ref, backend$ref) ||
        !identical(plan$links, links_i) ||
        B > plan$B) {
      ipw_raw_backend_abort("Sealed raw-IPW IF plan does not match the requested backend, links, or chunk bound.", "integrity")
    }
    pre <- plan$pre
    link_ids <- plan$link_ids
    model_pos <- plan$model_pos
    ipw_owner_safe_sum(c(backend$accounting$owned_entries, backend$accounting$borrowed_entries,
      pre$output, pre$scratch, plan$entries), backend$control$max_matrix_entries,
      "raw-IPW sealed plan plus IF output and scratch")
  } else {
    pre <- ipw_raw_if_budget(backend, B, links_i, details = want_details)
    lk <- backend$links[links_i]
    link_ids <- vapply(lk, `[[`, character(1L), "link_id")
    model_pos <- match(vapply(lk, `[[`, integer(1L), "model"), pre$tids)
  }

  frame <- get("frame_owner", envir = backend$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", envir = backend$panel_ref, inherits = FALSE)
  u_a_cohort <- frame$a_cohort_index[rows_i]
  u_w <- frame$w_a[rows_i]
  Z_a <- frame$Z_a

  # 1. Query each active pair once on `rows` and read outcome differences once
  #    on the union of observed active-support rows for that pair.
  n_pids <- length(pre$pids)
  pair_hits <- vector("list", n_pids)
  for (p_idx in seq_len(n_pids)) {
    pid <- pre$pids[p_idx]
    pair_rec <- backend$pairs[[pid]]
    hit_p <- ipw_raw_query_pair_rows(panel, pair_rec, rows_i)
    # Active cohorts for this pair across `pre$sids`
    active_sids_p <- pre$sids[vapply(backend$support_specs[pre$sids], function(s) s$pair_id == pid, logical(1L))]
    active_cohorts_p <- vapply(backend$support_specs[active_sids_p], `[[`, integer(1L), "cohort")
    obs_active <- which(hit_p$obs & (u_a_cohort %in% active_cohorts_p))
    dy_vec <- numeric(B)
    if (length(obs_active) > 0L) {
      ep <- ipw_owner_read_y_positions(panel, hit_p$pos_s[obs_active], hit_p$pos_t[obs_active])
      dy_sub <- as.double(ep$y_right_obs) - as.double(ep$y_left_obs)
      if (anyNA(dy_sub) || any(!is.finite(dy_sub))) {
        cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
      }
      dy_vec[obs_active] <- dy_sub
    }
    pair_hits[[p_idx]] <- list(obs = hit_p$obs, dy = dy_vec)
  }

  # 2. Rate channel (`B x |pre$sids|`) for diagnostic/support-column inspection
  up <- if (want_details) {
    up_mat <- matrix(0.0, nrow = B, ncol = length(pre$sids))
    for (j in seq_along(pre$sids)) {
      sid <- pre$sids[j]
      s <- backend$support_specs[[sid]]
      u <- backend$supports[[sid]]
      p_idx <- match(s$pair_id, pre$pids)
      ix <- which(u_a_cohort == s$cohort)
      if (length(ix) > 0L) {
        s_ind <- as.double(pair_hits[[p_idx]]$obs[ix])
        up_mat[ix, j] <- u_w[ix] * (s_ind - u$rate) / u$massA
      }
    }
    up_mat
  } else NULL

  # 3. Evaluate each distinct active treatment model `(pair, g, K)` once using
  #    the exact standalone four-channel arithmetic and addition order.
  n_tids <- length(pre$tids)
  model_phi <- matrix(0.0, nrow = B, ncol = n_tids)
  model_rate <- if (want_details) matrix(0.0, nrow = B, ncol = n_tids) else NULL
  model_score <- if (want_details) matrix(0.0, nrow = B, ncol = n_tids) else NULL
  model_treated <- if (want_details) matrix(0.0, nrow = B, ncol = n_tids) else NULL
  model_control <- if (want_details) matrix(0.0, nrow = B, ncol = n_tids) else NULL

  # Preserve the accepted shared-nuisance score reuse: only one bounded design,
  # raw probability and effective-odds buffer is live. Different pairs keep
  # their own outcomes, partials, model columns and original link ordering.
  nuisance_ids <- vapply(backend$treatment_specs[pre$tids], `[[`, integer(1L), "nuisance_idx")
  nuisance_buffer <- NULL
  current_nuisance <- NA_integer_
  nuisance_evaluations <- 0L
  for (t_idx in seq_len(n_tids)) {
    tj <- pre$tids[t_idx]
    spec <- backend$treatment_specs[[tj]]
    f <- backend$treatments[[tj]]
    p_idx <- match(spec$pair_id, pre$pids)
    ph <- pair_hits[[p_idx]]

    if (!identical(current_nuisance, nuisance_ids[t_idx])) {
      nuisance_buffer <- NULL
      nuisance_buffer <- local({
        gi <- match(u_a_cohort, spec$cohorts)
        in_risk <- which(!is.na(gi))
        pos <- in_risk[which(ph$obs[in_risk])]
        gj <- gi[pos]
        Z <- Z_a[rows_i[pos], , drop = FALSE]
        lv <- ipw_logit_values(Z, f$coefficients)
        r <- f$rates_clip[gj]
        w <- u_w[pos]
        list(gi = gi, in_risk = in_risk, pos = pos, gj = gj,
          Z = Z, r = r, w = w, odds = lv$odds,
          score_base = (w / r) * (as.double(gj == 1L) - lv$prob))
      })
      current_nuisance <- nuisance_ids[t_idx]
      nuisance_evaluations <- nuisance_evaluations + 1L
    } else if (!identical(nuisance_buffer$pos,
        nuisance_buffer$in_risk[which(ph$obs[nuisance_buffer$in_risk])])) {
      ipw_raw_backend_abort("Shared-nuisance pair mask changed during replay.", "integrity")
    }
    gi_chunk <- nuisance_buffer$gi
    in_risk_chunk <- nuisance_buffer$in_risk

    col_rate <- numeric(B)
    col_score <- numeric(B)
    col_treated <- numeric(B)
    col_control <- numeric(B)

    if (length(in_risk_chunk) > 0L) {
      gj_all <- gi_chunk[in_risk_chunk]
      col_rate[in_risk_chunk] <- u_w[in_risk_chunk] * f$loadings$rate_coef_s0[gj_all]

      pos_k <- nuisance_buffer$pos
      if (length(pos_k) > 0L) {
        gj_k <- nuisance_buffer$gj
        dy_k <- ph$dy[pos_k]
        w_k <- nuisance_buffer$w
        r_k <- nuisance_buffer$r
        odds_k <- nuisance_buffer$odds

        col_rate[pos_k] <- w_k * f$loadings$rate_coef_s1[gj_k]
        z_load_k <- drop(nuisance_buffer$Z %*% f$loadings$loading_beta)
        col_score[pos_k] <- nuisance_buffer$score_base * z_load_k

        is_tk <- gj_k == 1L
        if (any(is_tk)) {
          pos_t <- pos_k[is_tk]
          col_treated[pos_t] <- (w_k[is_tk] / r_k[is_tk]) / f$means$mean_a *
            (dy_k[is_tk] - f$means$treated_mean)
        }
        is_ck <- !is_tk
        if (any(is_ck)) {
          pos_c <- pos_k[is_ck]
          col_control[pos_c] <- -(w_k[is_ck] * odds_k[is_ck] / r_k[is_ck]) /
            f$means$mean_b * (dy_k[is_ck] - f$means$control_mean)
        }
      }
    }

    col_phi <- col_treated + col_control + col_score + col_rate
    if (anyNA(col_phi) || any(!is.finite(col_phi))) {
      cdid_abort("Nonfinite IF contribution encountered during replay.", c("cdid_ipw_numeric_error"))
    }
    model_phi[, t_idx] <- col_phi
    if (want_details) {
      model_rate[, t_idx] <- col_rate
      model_score[, t_idx] <- col_score
      model_treated[, t_idx] <- col_treated
      model_control[, t_idx] <- col_control
    }
  }
  nuisance_buffer <- NULL

  # 4. Assemble requested link columns in original request order
  phi <- model_phi[, model_pos, drop = FALSE]
  dimnames(phi) <- list(NULL, link_ids)

  comp <- if (want_details) {
    cr <- model_rate[, model_pos, drop = FALSE]
    cs <- model_score[, model_pos, drop = FALSE]
    ct <- model_treated[, model_pos, drop = FALSE]
    cc <- model_control[, model_pos, drop = FALSE]
    dimnames(cr) <- dimnames(cs) <- dimnames(ct) <- dimnames(cc) <- list(NULL, link_ids)
    list(
      rate = cr,
      beta_score = cs,
      treated = ct,
      control = cc
    )
  } else NULL

  out <- list(
    rows = rows_i,
    a_index = rows_i,
    unit_index = panel$a_unit_index[rows_i],
    id = panel$a_id[rows_i],
    g = panel$a_g[rows_i],
    weight = u_w,
    links = links_i,
    link_ids = link_ids,
    phi = phi,
    scale = scale,
    normalization_count = n_a,
    components = comp,
    details = if (want_details) list(components = comp, upstream = up, psi_ipw = phi) else NULL,
    diagnostics = list(
      n_rate_score_evaluations = length(pre$sids),
      n_regression_score_evaluations = 0L,
      n_global_beta_blocks = 0L,
      n_treatment_score_evaluations = length(pre$tids),
      n_nuisance_score_evaluations = nuisance_evaluations,
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  out
}

ipw_raw_link_estimate <- function(backend, link_index, details = FALSE) {
  validate_ipw_raw_backend(backend)
  if (!is.numeric(link_index) || is.object(link_index) || !is.null(dim(link_index)) ||
      length(link_index) != 1L || is.na(link_index) || !is.finite(link_index) ||
      link_index != floor(link_index) || link_index < 1 || link_index > length(backend$links)) {
    ipw_raw_backend_abort("Invalid `link_index`.", "input")
  }
  l <- backend$links[[as.integer(link_index)]]
  f <- backend$treatments[[l$model]]
  if (!isTRUE(details)) return(as.double(f$means$estimate))
  list(
    link_id = l$link_id,
    estimate = f$means$estimate,
    treated_mean = f$means$treated_mean,
    control_mean = f$means$control_mean,
    dp = f$means$mean_a,
    dm = f$means$mean_b
  )
}
