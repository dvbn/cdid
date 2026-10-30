# ==============================================================================
# Chained Custom-IPW Callback Bridge (`cdid/R/chained_custom_ipw.R`)
# ==============================================================================
# Implements the CP02 custom-IPW callback seam on the D23 candidate:
#   - Calls user-supplied did-style panel callback
#       f(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...)
#     once per link on the observed risk set O_ell = {i : R_{i,ell} S_{i,ell} = 1}
#     in canonical increasing A-unit order with calibrated weights
#       a_i = w_i / hat{r}_{ell, G_i}.
#   - Validates that the callback returns a list with finite scalar `ATT` and
#     finite length-m_ell vector (or m_ell x 1 matrix) `att.inf.func` (never
#     centering away an uncentered IF, never substituting built-in IPW, and
#     never regularizing/falling back).
#   - Propagates cohort-rate estimation uncertainty via
#       b_{ell,k} = -(m_ell * hat{r}_{ell,k})^{-1} sum_{j : G_j = k} psi_j
#     or, when the callback returns the optional O15 field
#     `log_weight_gradient` g (dATT/dlog i.weights, passed-row order),
#       b_{ell,k} = -hat{r}_{ell,k}^{-1} sum_{j : G_j = k} g_j,
#       phi_{i,ell} = (n / m_ell) * 1{i in O_ell} * psi_{pos(i)}
#                     + sum_{k in K_ell} b_{ell,k} * rho_{i,k},
#     where rho_{i,k} = w_i * 1{G_i = k} * (S_{i,ell} - hat{r}_{ell,k}) / q_k.
#   - Provides both the bounded dense-source (`cdid_ipw_context`) interface and
#     the owner/view-backed (`cdid_ipw_frame_owner` / `cdid_ipw_link_view`)
#     pipeline interface, sharing the exact same composition and chunk-replay
#     workers. Retains a bounded O(m_ell + |K_ell|) per-link representation
#     after releasing callback input arrays, with checked preflight before any
#     callback invocation.
# ==============================================================================

cdid_custom_ipw_abort <- function(message, class = character(), diagnostics = NULL) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_custom_ipw_error",
      "cdid_ipw_point_view_error",
      "cdid_ipw_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

custom_ipw_validate_callback_output <- function(res, n_obs, link_id = NULL) {
  where <- if (!is.null(link_id) && nzchar(link_id)) paste0(" on link `", link_id, "`") else ""
  if (is.null(res) || !is.list(res) || is.data.frame(res)) {
    cdid_custom_ipw_abort(
      paste0("Custom `est_method` callback", where, " must return a list containing `ATT` and `att.inf.func`."),
      c("cdid_custom_ipw_output_error", "cdid_ipw_fit_error")
    )
  }
  nms <- names(res)
  if (is.null(nms) || !all(c("ATT", "att.inf.func") %in% nms)) {
    cdid_custom_ipw_abort(
      paste0("Custom `est_method` callback result", where, " must contain named elements `ATT` and `att.inf.func`."),
      c("cdid_custom_ipw_output_error", "cdid_ipw_fit_error")
    )
  }
  att <- res[["ATT"]]
  if (!is.numeric(att) || !is.atomic(att) || is.object(att) || !is.null(dim(att)) ||
      length(att) != 1L || is.na(att) || !is.finite(att)) {
    cdid_custom_ipw_abort(
      paste0("Custom `est_method` callback `ATT`", where, " must be a finite numeric scalar."),
      c("cdid_custom_ipw_output_error", "cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  inf <- res[["att.inf.func"]]
  if (!is.numeric(inf) || !is.atomic(inf) || is.object(inf)) {
    cdid_custom_ipw_abort(
      paste0("Custom `est_method` callback `att.inf.func`", where, " must be a finite numeric vector or 1-column matrix."),
      c("cdid_custom_ipw_output_error", "cdid_ipw_fit_error")
    )
  }
  d_inf <- dim(inf)
  valid_shape <- (is.null(d_inf) && length(inf) == as.integer(n_obs)) ||
    (length(d_inf) == 2L && identical(as.integer(d_inf), c(as.integer(n_obs), 1L)))
  if (!valid_shape) {
    cdid_custom_ipw_abort(
      paste0(
        "Custom `est_method` callback `att.inf.func`", where,
        " must have length equal to the observed risk-set size (", as.integer(n_obs), ")."
      ),
      c("cdid_custom_ipw_dimension_error", "cdid_custom_ipw_output_error", "cdid_ipw_fit_error")
    )
  }
  psi <- as.double(inf)
  if (anyNA(psi) || any(!is.finite(psi))) {
    cdid_custom_ipw_abort(
      paste0("Custom `est_method` callback `att.inf.func`", where, " contains nonfinite (`NA`/`NaN`/`Inf`) entries."),
      c("cdid_custom_ipw_output_error", "cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  # O15 (D64): optional total log-weight gradient g_j = dATT / d log(i.weights[j])
  # in passed-row order, the callback fully re-run (internal nuisance refits
  # included); ATT units, never centred here.
  g <- res[["log_weight_gradient"]]
  if (!is.null(g)) {
    d_g <- dim(g)
    if (!is.numeric(g) || !is.atomic(g) || is.object(g) ||
        !((is.null(d_g) && length(g) == as.integer(n_obs)) ||
          (length(d_g) == 2L && identical(as.integer(d_g), c(as.integer(n_obs), 1L))))) {
      cdid_custom_ipw_abort(
        paste0("Custom `est_method` callback `log_weight_gradient`", where, " must be NULL or a numeric vector (or 1-column ",
          "matrix) of length equal to the observed risk-set size (", as.integer(n_obs), "), in the order of the passed rows."),
        c("cdid_custom_ipw_dimension_error", "cdid_custom_ipw_output_error", "cdid_ipw_fit_error")
      )
    }
    g <- as.double(g)
    if (anyNA(g) || any(!is.finite(g))) {
      cdid_custom_ipw_abort(
        paste0("Custom `est_method` callback `log_weight_gradient`", where, " contains nonfinite (`NA`/`NaN`/`Inf`) entries."),
        c("cdid_custom_ipw_output_error", "cdid_ipw_numeric_error", "cdid_ipw_fit_error")
      )
    }
  }
  list(ATT = as.double(att), psi = psi, g = g)
}

# Cohort-rate sensitivity of a baseline callback link, b_k = sum_j g_j d log a_j /
# d r_k with a_j = w_j / max(r_{G_j}, eps): the explicit total gradient when supplied (O15),
# otherwise the existing identity rule g = psi / m (the two are never added).
# A floored rate (r_k <= eps) has zero package-owned derivative.
custom_ipw_rate_sensitivity <- function(val, cohort_index, rates, m_cohorts, n_obs) {
  explicit <- !is.null(val$g)
  active <- ipw_clip_floor_active(rates)
  vapply(seq_len(m_cohorts), function(k) {
    if (active[k] == 0) 0 else
    if (explicit) -sum(val$g[cohort_index == k]) / rates[k] else
      -sum(val$psi[cohort_index == k]) / (as.double(n_obs) * rates[k])
  }, numeric(1L))
}

custom_ipw_derivative_source <- function(val) {
  if (is.null(val$g)) "identity_psi_over_m (existing baseline rule)" else "explicit_log_weight_gradient"
}

# Shared observer-buffer custom-IPW worker used by both dense-context and
# owner/view-backed link fitting.
custom_ipw_observer_fit_worker <- function(
  sup,
  risk_a_idx,
  risk_gi,
  risk_w,
  y1,
  y0,
  Z_use,
  custom_est_method,
  custom_dots = list(),
  link_id = NULL
) {
  n_obs <- length(risk_a_idx)
  if (n_obs < 2L) {
    cdid_custom_ipw_abort(
      "Custom IPW link has fewer than 2 observed risk-set units.",
      c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }
  pr <- ipw_pair_rates(risk_gi, risk_w, sup)
  rates <- as.double(unname(pr$rates))
  if (any(pr$pair_count < 1L) || any(!is.finite(rates) | rates <= 0 | rates > 1)) {
    cdid_custom_ipw_abort(
      "Active link cohort has zero observed-pair support or invalid marginal observation rate.",
      c("cdid_ipw_owner_pair_support_error", "cdid_ipw_rate_error", "cdid_ipw_support_error")
    )
  }

  r_use <- ipw_clip_floor(rates)[risk_gi]
  fit_weights <- risk_w / r_use
  if (any(!is.finite(fit_weights) | fit_weights <= 0)) {
    cdid_abort(
      "Nonfinite or nonpositive inverse-rate propensity weights `w_i / r_G`.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  D_obs <- as.numeric(risk_gi == 1L)

  cb_args <- c(
    list(
      y1 = y1,
      y0 = y0,
      D = D_obs,
      covariates = Z_use,
      i.weights = fit_weights,
      inffunc = TRUE
    ),
    custom_dots
  )
  raw_res <- do.call(custom_est_method, cb_args)
  rm(cb_args)

  val <- custom_ipw_validate_callback_output(raw_res, n_obs = n_obs, link_id = link_id)
  rm(raw_res)

  m_cohorts <- as.integer(sup$m_cohorts)
  b_rate <- custom_ipw_rate_sensitivity(val, risk_gi, rates, m_cohorts, n_obs)
  if (anyNA(b_rate) || any(!is.finite(b_rate))) {
    where <- if (!is.null(link_id) && nzchar(link_id)) paste0(" on link `", link_id, "`") else ""
    cdid_custom_ipw_abort(
      paste0("Nonfinite cohort-rate sensitivity `b_rate`", where, "."),
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  names(b_rate) <- sup$cohort_labels

  list(
    estimate = as.double(val$ATT),
    psi_obs = as.double(val$psi),
    callback_if_mean = mean(val$psi),
    derivative_source = custom_ipw_derivative_source(val),
    b_rate = b_rate,
    pair_count = as.integer(pr$pair_count),
    pair_mass = as.double(pr$pair_mass),
    rates = rates
  )
}

# Shared per-chunk IF evaluation worker used by both dense-context and
# owner/view-backed IF replay.
custom_ipw_eval_link_chunk <- function(
  a_s,
  a_e,
  u_w,
  gi_chunk,
  obs_a_idx,
  psi_obs,
  b_rate,
  q_vec,
  rates_vec,
  n_a_num,
  m_obs_num,
  include_components = FALSE
) {
  b_len <- a_e - a_s + 1L
  in_risk <- which(!is.na(gi_chunk))
  col_phi <- numeric(b_len)
  col_base <- if (include_components) numeric(b_len) else NULL
  col_rate <- if (include_components) numeric(b_len) else NULL

  if (length(in_risk) > 0L) {
    S_chunk <- numeric(b_len)
    psi_chunk <- numeric(b_len)
    lo <- findInterval(a_s - 1L, obs_a_idx) + 1L
    hi <- findInterval(a_e, obs_a_idx)
    if (lo <= hi) {
      pos <- seq.int(lo, hi)
      rel <- obs_a_idx[pos] - a_s + 1L
      S_chunk[rel] <- 1.0
      psi_chunk[rel] <- psi_obs[pos]
    }
    gk <- gi_chunk[in_risk]
    base_val <- (n_a_num / m_obs_num) * psi_chunk[in_risk]
    rho_val <- u_w[in_risk] * (S_chunk[in_risk] - rates_vec[gk]) / q_vec[gk]
    rate_val <- rho_val * b_rate[gk]
    col_phi[in_risk] <- base_val + rate_val
    if (include_components) {
      col_base[in_risk] <- base_val
      col_rate[in_risk] <- rate_val
    }
  }
  if (anyNA(col_phi) || any(!is.finite(col_phi))) {
    cdid_abort("Nonfinite IF contribution encountered during custom-IPW replay.", "cdid_ipw_numeric_error")
  }
  list(phi = col_phi, base = col_base, rate = col_rate)
}

# ==============================================================================
# 1. Bounded Dense-Source (`cdid_ipw_context`) Custom-IPW Interface
# ==============================================================================

preflight_custom_ipw_link <- function(context, g, K, s, t, max_entries = NULL) {
  validate_ipw_context(context)
  coh <- validate_ipw_cohorts(g, K)
  desc <- make_date_pair_descriptor(s, t, context)
  sup <- ipw_cohort_support(context, coh)
  sup$cohort_labels <- ipw_cohort_labels(coh)

  hit <- query_pair_mask_internal(
    context, desc$s, desc$t,
    ds = desc$s_calendar_index, dt = desc$t_calendar_index
  )
  risk <- ipw_pair_risk_set(context, coh$cohorts, hit)
  n_obs <- length(risk$risk_a_idx)
  if (n_obs < 2L) {
    cdid_custom_ipw_abort(
      "Custom IPW link has fewer than 2 observed risk-set units.",
      c("cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }
  p <- as.integer(context$n_covariates)
  m <- as.integer(sup$m_cohorts)
  limit <- 2^53 - 1
  if (as.double(n_obs) > limit / max(as.double(p) + 8, 1)) {
    cdid_abort(
      "Custom-IPW observer buffer exceeds supported integer/double dimensions.",
      c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
    )
  }
  observer_buffer_entries <- as.double(n_obs) * as.double(p) + 8 * as.double(n_obs) + 64
  retained_entries <- 2 * as.double(n_obs) + 8 * as.double(m) + 64
  peak_concurrent_entries <- observer_buffer_entries + retained_entries

  if (!is.null(max_entries)) {
    bound <- ipw_bound(max_entries, "max_entries", allow_zero = TRUE)
    if (!is.finite(peak_concurrent_entries) || peak_concurrent_entries > bound ||
        peak_concurrent_entries > .Machine$integer.max) {
      cdid_abort(
        "Custom-IPW link preflight exceeds explicit `max_entries` bound.",
        c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
      )
    }
  }

  list(
    coh = coh,
    desc = desc,
    sup = sup,
    hit = hit,
    risk = risk,
    preflight = list(
      n_a_units = context$n_a_units,
      n_observed_risk_units = n_obs,
      p_covariates = p,
      m_cohorts = m,
      observer_buffer_entries = observer_buffer_entries,
      retained_entries_per_link = retained_entries,
      peak_concurrent_entries = peak_concurrent_entries,
      peak_rss_promised = FALSE,
      memory_disclosure = paste0(
        "Deterministic preflight charges the temporary observer buffer (y1, y0, D, covariates, i.weights) ",
        "and the retained (obs_a_idx, psi_obs, b_rate) link representation before invoking the user callback. ",
        "Arbitrary workspace allocated inside the opaque user callback and OS process RSS are outside ",
        "package-owned array accounting."
      )
    )
  )
}

fit_custom_ipw_link <- function(
  context,
  g,
  K,
  s,
  t,
  custom_est_method,
  custom_dots = list(),
  link_id = NULL,
  functional_key = NULL,
  max_entries = NULL
) {
  validate_ipw_context(context)
  if (missing(custom_est_method) || !is.function(custom_est_method)) {
    cdid_custom_ipw_abort("`custom_est_method` must be a function.", "cdid_ipw_input_error")
  }
  if (!is.list(custom_dots) || is.object(custom_dots)) {
    cdid_custom_ipw_abort("`custom_dots` must be a plain list.", "cdid_ipw_input_error")
  }
  validate_ipw_link_labels(link_id, functional_key)

  # 1. Checked preflight strictly BEFORE reading outcomes or invoking the callback
  pre_bundle <- preflight_custom_ipw_link(context, g, K, s, t, max_entries = max_entries)
  coh <- pre_bundle$coh
  desc <- pre_bundle$desc
  sup <- pre_bundle$sup
  risk <- pre_bundle$risk
  pre <- pre_bundle$preflight

  # 2. Read pair outcomes for observed pair mask and slice risk-set observers in canonical A-unit order
  pq <- ipw_pair_outcomes(context, pre_bundle$hit)
  in_risk <- risk$in_risk
  risk_a_idx <- as.integer(risk$risk_a_idx)
  risk_gi <- as.integer(risk$risk_gi)
  risk_w <- as.double(context$w_a[risk_a_idx])
  y0 <- as.double(pq$y_s[in_risk])
  y1 <- as.double(pq$y_t[in_risk])
  Z_use <- context$Z_a[risk_a_idx, , drop = FALSE]
  rm(pq, pre_bundle)

  # 3. Invoke callback and compute cohort-rate sensitivity b_{ell,k}
  w_out <- custom_ipw_observer_fit_worker(
    sup = sup,
    risk_a_idx = risk_a_idx,
    risk_gi = risk_gi,
    risk_w = risk_w,
    y1 = y1,
    y0 = y0,
    Z_use = Z_use,
    custom_est_method = custom_est_method,
    custom_dots = custom_dots,
    link_id = link_id
  )
  rm(y0, y1, Z_use, risk_w, risk_gi)

  cohort_labels <- sup$cohort_labels
  group_a_count <- as.integer(sup$group_a_count)
  group_mass <- as.double(sup$group_mass)
  q_vec <- as.double(sup$q)
  pair_count <- as.integer(w_out$pair_count)
  pair_mass <- as.double(w_out$pair_mass)
  rates_named <- as.double(w_out$rates)
  names(group_a_count) <- names(pair_count) <- cohort_labels
  names(group_mass) <- names(pair_mass) <- names(q_vec) <- names(rates_named) <- cohort_labels

  n_obs <- length(risk_a_idx)
  fit_obj <- structure(
    list(
      link_id = link_id,
      functional_key = functional_key,
      nuisance_id = paste0("custom_ipw:", if (!is.null(link_id) && nzchar(link_id)) link_id else paste0(desc$s, "_", desc$t)),
      fitted_nuisance_ref = NULL,
      nuisance_view_ref = NULL,
      nuisance_ref = NULL,
      estimate = w_out$estimate,
      pair = list(
        s = desc$s,
        t = desc$t,
        s_calendar_index = desc$s_calendar_index,
        t_calendar_index = desc$t_calendar_index,
        link_id = link_id,
        functional_key = functional_key
      ),
      shared_across_dates = FALSE,
      nuisance_mask_verified = FALSE,
      cohorts = list(
        g = coh$g,
        K = coh$K,
        cohorts = coh$cohorts,
        cohort_labels = cohort_labels,
        n_cohorts = sup$m_cohorts,
        m_cohorts = sup$m_cohorts
      ),
      dimensions = list(
        n_a_units = context$n_a_units,
        n_observed_risk_units = n_obs,
        p_covariates = context$n_covariates,
        m_cohorts = sup$m_cohorts
      ),
      counts = list(
        n_a_units = context$n_a_units,
        n_risk_a_units = sum(group_a_count),
        n_outside_risk_a_units = context$n_a_units - sum(group_a_count),
        n_observed_pair_units = risk$n_obs_pair_all,
        n_observed_risk_units = n_obs,
        n_nonobserved_risk_units = sum(group_a_count) - n_obs,
        n_treated_obs = pair_count[1L],
        n_control_obs = sum(pair_count[-1L]),
        by_cohort_a = group_a_count,
        by_cohort_pair = pair_count
      ),
      masses = list(
        a_weight_total = context$sum_w_a,
        cohort_a_weight = group_mass,
        cohort_pair_weight = pair_mass,
        q = q_vec,
        rates = rates_named
      ),
      obs_a_idx = risk_a_idx,
      psi_obs = w_out$psi_obs,
      b_rate = w_out$b_rate,
      diagnostics = list(callback_if_mean = w_out$callback_if_mean, derivative_source = w_out$derivative_source),
      preflight = pre,
      accounting = list(
        retained_entries = pre$retained_entries_per_link,
        peak_callback_input_entries = pre$observer_buffer_entries,
        callback_internal_workspace = "untracked_opaque_callback",
        endpoints_read = as.integer(2L * n_obs),
        observer_released = TRUE
      ),
      context_ref = context$ref
    ),
    class = c("cdid_custom_ipw_link_fit", "cdid_custom_ipw_link")
  )
  ipw_seal(fit_obj, "custom_ipw_link")
}

estimate_custom_ipw_link <- fit_custom_ipw_link

normalize_custom_ipw_fits <- function(context, fits) {
  validate_ipw_context(context)
  single <- inherits(fits, "cdid_custom_ipw_link_fit")
  fit_list <- if (single) list(fits) else fits
  if (!is.list(fit_list) || length(fit_list) < 1L ||
      any(!vapply(fit_list, inherits, logical(1L), what = "cdid_custom_ipw_link_fit"))) {
    cdid_abort(
      "`fits` must be a `cdid_custom_ipw_link_fit` or a nonempty list of `cdid_custom_ipw_link_fit` objects.",
      c("cdid_ipw_input_error")
    )
  }
  retained_entries <- 0.0
  for (ell in seq_along(fit_list)) {
    validate_custom_ipw_link(fit_list[[ell]])
    if (!identical(fit_list[[ell]]$context_ref, context$ref)) {
      cdid_abort(
        "All fitted custom-IPW link caches must share the exact same `cdid_ipw_context` instance.",
        c("cdid_ipw_context_mismatch_error", "cdid_ipw_context_error")
      )
    }
    retained_entries <- retained_entries + as.double(fit_list[[ell]]$accounting$retained_entries)
  }
  link_ids <- vapply(seq_along(fit_list), function(ell) {
    id_ell <- fit_list[[ell]]$pair$link_id
    if (!is.null(id_ell) && nzchar(id_ell)) id_ell else paste0("link_", ell)
  }, character(1L))
  if (!is.null(names(fit_list)) && all(nzchar(names(fit_list)))) {
    link_ids <- names(fit_list)
  }
  if (anyDuplicated(link_ids)) {
    cdid_custom_ipw_abort("Duplicate `link_id` across supplied `cdid_custom_ipw_link_fit` objects.", "cdid_ipw_input_error")
  }
  list(
    fits = fit_list,
    link_ids = link_ids,
    single = single,
    retained_entries = retained_entries
  )
}

preflight_custom_ipw_replay <- function(context, fits, chunk_size, max_entries = NULL) {
  norm <- normalize_custom_ipw_fits(context, fits)
  csize <- ipw_bound(chunk_size, "chunk_size")
  n_a <- as.double(context$n_a_units)
  n_links <- as.double(length(norm$fits))
  eff_chunk <- min(csize, n_a)
  n_batches <- ceiling(n_a / eff_chunk)
  limit <- 2^53 - 1
  if (eff_chunk > limit / max(n_links, 1) || eff_chunk * n_links > .Machine$integer.max) {
    cdid_abort(
      "Single replay batch matrix exceeds supported integer/double dimensions.",
      c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
    )
  }
  batch_entries <- eff_chunk * n_links
  full_entries <- if (n_a <= limit / max(n_links, 1)) n_a * n_links else Inf

  if (!is.null(max_entries)) {
    bound <- ipw_bound(max_entries, "max_entries", allow_zero = TRUE)
    if (!is.finite(full_entries) || full_entries > bound || full_entries > .Machine$integer.max) {
      cdid_abort(
        "Custom-IPW IF materialization exceeds explicit entry bound.",
        c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
      )
    }
  }

  list(
    n_a_units = context$n_a_units,
    n_total_units = context$n_total_units,
    n_singleton_units = context$n_singleton_units,
    n_unobserved_units = context$n_unobserved_units,
    n_links = length(norm$fits),
    n_covariates = context$n_covariates,
    chunk_size = as.double(csize),
    effective_chunk_size = as.integer(eff_chunk),
    n_batches = as.integer(n_batches),
    batch_matrix_entries = as.double(batch_entries),
    full_matrix_entries = as.double(full_entries),
    retained_link_entries = norm$retained_entries,
    covariance_matrix_entries = if (n_links <= limit / n_links) n_links^2 else Inf,
    link_ids = norm$link_ids,
    peak_rss_promised = FALSE
  )
}

replay_custom_ipw_if <- function(
  context,
  fits,
  chunk_size,
  callback,
  include_components = FALSE
) {
  norm <- normalize_custom_ipw_fits(context, fits)
  pre <- preflight_custom_ipw_replay(context, norm$fits, chunk_size)
  if (!is.function(callback)) {
    cdid_abort("`callback` must be a function.", c("cdid_ipw_input_error", "cdid_pairs_data_error"))
  }
  validate_ipw_components(include_components)

  fit_list <- norm$fits
  link_ids <- norm$link_ids
  n_links <- length(fit_list)
  n_a <- context$n_a_units
  n_a_num <- as.double(n_a)
  step <- pre$effective_chunk_size

  emitted_units <- 0L
  batches <- 0L
  max_batch_units <- 0L

  a_start <- 1L
  while (a_start <= n_a) {
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    a_seq <- seq.int(a_start, a_end)
    b_len <- length(a_seq)
    u_idx <- context$a_unit_index[a_seq]
    u_id <- context$a_id[a_seq]
    u_g <- context$a_g[a_seq]
    u_w <- context$w_a[a_seq]

    phi_mat <- matrix(0.0, nrow = b_len, ncol = n_links, dimnames = list(NULL, link_ids))
    if (include_components) {
      comp_base <- matrix(0.0, nrow = b_len, ncol = n_links, dimnames = list(NULL, link_ids))
      comp_rate <- matrix(0.0, nrow = b_len, ncol = n_links, dimnames = list(NULL, link_ids))
    }

    for (ell in seq_len(n_links)) {
      fit <- fit_list[[ell]]
      gi_chunk <- match(u_g, fit$cohorts$cohorts)
      ev <- custom_ipw_eval_link_chunk(
        a_s = a_start,
        a_e = a_end,
        u_w = u_w,
        gi_chunk = gi_chunk,
        obs_a_idx = fit$obs_a_idx,
        psi_obs = fit$psi_obs,
        b_rate = as.double(fit$b_rate),
        q_vec = as.double(fit$masses$q),
        rates_vec = as.double(fit$masses$rates),
        n_a_num = n_a_num,
        m_obs_num = as.double(fit$counts$n_observed_risk_units),
        include_components = include_components
      )
      phi_mat[, ell] <- ev$phi
      if (include_components) {
        comp_base[, ell] <- ev$base
        comp_rate[, ell] <- ev$rate
      }
    }

    batches <- batches + 1L
    batch <- structure(
      list(
        batch_index = batches,
        a_index = a_seq,
        unit_index = u_idx,
        id = u_id,
        g = u_g,
        weight = u_w,
        phi = phi_mat,
        components = if (include_components) {
          list(
            base = comp_base,
            rate = comp_rate
          )
        } else NULL
      ),
      class = "cdid_ipw_if_batch"
    )
    callback(batch)
    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != n_a) {
    cdid_abort("Custom-IPW IF replay unit count mismatch.", c("cdid_ipw_count_error"))
  }
  invisible(list(
    n_a_units = emitted_units,
    n_links = n_links,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * as.double(n_links)
  ))
}

traverse_custom_ipw_if <- replay_custom_ipw_if

materialize_custom_ipw_if <- function(
  context,
  fits,
  max_entries,
  chunk_size,
  include_components = FALSE
) {
  if (missing(max_entries)) {
    cdid_abort(
      "Supply an explicit `max_entries` bound for custom-IPW IF materialization.",
      c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
    )
  }
  bound <- ipw_bound(max_entries, "max_entries", allow_zero = TRUE)
  validate_ipw_components(include_components)
  norm <- normalize_custom_ipw_fits(context, fits)
  pre <- preflight_custom_ipw_replay(context, norm$fits, chunk_size, max_entries = bound)
  pre$materialized_matrix_count <- if (include_components) 3L else 1L
  pre$materialized_matrix_entries <- pre$full_matrix_entries * pre$materialized_matrix_count
  n_a <- pre$n_a_units
  n_links <- pre$n_links
  link_ids <- pre$link_ids

  phi <- matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids))
  comp <- if (include_components) {
    list(
      base = matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids)),
      rate = matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids))
    )
  } else NULL

  replay_info <- replay_custom_ipw_if(
    context = context,
    fits = norm$fits,
    chunk_size = chunk_size,
    include_components = include_components,
    callback = function(batch) {
      idx <- batch$a_index
      phi[idx, ] <<- batch$phi
      if (include_components) {
        comp$base[idx, ] <<- batch$components$base
        comp$rate[idx, ] <<- batch$components$rate
      }
    }
  )

  structure(
    list(
      a_index = seq_len(n_a),
      unit_index = context$a_unit_index,
      id = context$a_id,
      g = context$a_g,
      weight = context$w_a,
      phi = if (norm$single) drop(phi) else phi,
      phi_matrix = phi,
      components = comp,
      estimates = stats::setNames(
        vapply(norm$fits, `[[`, numeric(1L), "estimate"),
        link_ids
      ),
      replay = replay_info,
      preflight = pre
    ),
    class = "cdid_ipw_if_materialized"
  )
}

stream_custom_ipw_covariance <- function(context, fits, chunk_size) {
  norm <- normalize_custom_ipw_fits(context, fits)
  pre <- preflight_custom_ipw_replay(context, norm$fits, chunk_size)
  n_a <- as.double(pre$n_a_units)
  n_links <- pre$n_links
  link_ids <- pre$link_ids

  if (!is.finite(pre$covariance_matrix_entries) ||
      pre$covariance_matrix_entries > .Machine$integer.max) {
    cdid_abort(
      "Covariance matrix exceeds supported integer/double dimensions.",
      c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
    )
  }
  sum_phi <- numeric(n_links)
  replay_custom_ipw_if(context, norm$fits, chunk_size, function(batch) {
    sum_phi <<- sum_phi + colSums(batch$phi)
    if (any(!is.finite(sum_phi))) {
      cdid_abort("Nonfinite IF sum during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
    }
  })
  mean_phi <- sum_phi / n_a
  names(mean_phi) <- link_ids

  cov_mat <- matrix(0.0, nrow = n_links, ncol = n_links, dimnames = list(link_ids, link_ids))
  replay_info <- replay_custom_ipw_if(context, norm$fits, chunk_size, function(batch) {
    centered <- sweep(batch$phi, 2L, mean_phi, FUN = "-")
    if (any(!is.finite(centered))) {
      cdid_abort("Nonfinite centered IF during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
    }
    cov_mat <<- cov_mat + crossprod(centered)
    if (any(!is.finite(cov_mat))) {
      cdid_abort("Nonfinite IF crossproduct or sum during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
    }
  })
  cov_mat <- cov_mat / n_a

  list(
    covariance = cov_mat,
    mean_phi = mean_phi,
    estimates = stats::setNames(
      vapply(norm$fits, `[[`, numeric(1L), "estimate"),
      link_ids
    ),
    n_a_units = pre$n_a_units,
    replay = replay_info
  )
}

# ==============================================================================
# 2. Owner/View-Backed (`cdid_ipw_frame_owner` / `cdid_ipw_link_view`) Interface
# ==============================================================================

preflight_custom_ipw_link_view <- function(
  frame_owner,
  link_view,
  control,
  n_planned_links = 1L
) {
  ipw_point_view_require_deps()
  if (missing(frame_owner) || missing(link_view) || missing(control)) {
    cdid_custom_ipw_abort(
      "Required arguments (`frame_owner`, `link_view`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  validate_ipw_frame_owner(frame_owner)
  validate_ipw_link_view(link_view, frame_owner = frame_owner)
  ctrl <- validate_ipw_point_view_control(control)
  n_links <- ipw_owner_int_bound(n_planned_links, "n_planned_links", min_val = 1)
  panel_owner <- get("panel_owner", envir = frame_owner$panel_ref, inherits = FALSE)

  na <- ipw_owner_int_bound(frame_owner$n_a_units, "n_a_units", min_val = 1)
  n_obs <- ipw_owner_int_bound(link_view$support$n_observed_risk_units, "n_observed_risk_units", min_val = 0)
  p <- ipw_owner_int_bound(frame_owner$n_covariates, "n_covariates", min_val = 1)
  m <- ipw_owner_int_bound(link_view$cohorts$m_cohorts, "m_cohorts", min_val = 2)

  max_ch <- ctrl$max_chunk_rows
  max_obs <- ctrl$max_observer_units
  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  if (n_obs < 2L) {
    cdid_custom_ipw_abort(
      "Custom IPW link view has fewer than 2 observed risk-set units.",
      c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }
  if (n_obs > max_obs) {
    cdid_custom_ipw_abort(
      paste0("Observed risk-set unit count (", n_obs, ") exceeds `control$max_observer_units` (", max_obs, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error", "cdid_ipw_bound_error")
    )
  }

  B_eff <- min(na, max_ch)
  n_chunks <- ceiling(na / B_eff)

  obs_Z_entries <- ipw_owner_check_matrix_dims(
    n_obs, p,
    max_matrix_entries = max_me,
    label = "custom-IPW observer buffer Z_use matrix"
  )
  observer_buffer_entries <- ipw_owner_safe_sum(
    c(
      obs_Z_entries,
      ipw_owner_safe_product(8, max(1, n_obs), max_me, "custom-IPW observer index/cohort/weight/y0/y1/D/a_obs/psi vectors"),
      64
    ),
    max_bound = max_me,
    label = "temporary custom-IPW observer buffer entries"
  )

  retained_entries_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(2, max(1, n_obs), min(max_ret, max_out, max_me), "custom-IPW retained obs_a_idx and psi_obs"),
      ipw_owner_safe_product(8, m, min(max_ret, max_out, max_me), "custom-IPW retained b_rate, q, rates, and cohort metadata"),
      64
    ),
    max_bound = min(max_ret, max_out, max_me),
    label = "single custom-IPW link retained entries"
  )
  total_retained_entries <- ipw_owner_safe_product(
    n_links, retained_entries_one,
    max_entries = min(max_ret, max_out, max_me),
    label = "total planned custom-IPW retained link entries"
  )

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
  borrowed_entries <- ipw_point_view_owner_entries(frame_owner, panel_owner)

  peak_concurrent_entries <- ipw_owner_safe_sum(
    c(
      total_retained_entries,
      borrowed_entries,
      2 * observer_buffer_entries,
      reader_chunk_concurrent_entries,
      4 * B_eff * p + 32 * B_eff,
      128
    ),
    max_bound = max_me,
    label = "simultaneous custom-IPW streaming and callback handoff entries"
  )
  peak_concurrent_matrix_entries <- ipw_owner_safe_sum(
    c(obs_Z_entries, chunk_pre$output_matrix_entries),
    max_bound = max_me,
    label = "simultaneous custom-IPW matrix entries"
  )

  stream_work_one <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(n_chunks, chunk_pre$work_ops_bound, max_wo, "bounded link-view chunk query work"),
      ipw_owner_safe_product(8, max(1, n_obs) * p, max_wo, "observer buffer copy and order verification work"),
      ipw_owner_safe_product(16, max(1, n_obs) * m, max_wo, "rate calibration and b_rate sensitivity work"),
      256
    ),
    max_bound = max_wo,
    label = "single custom-IPW link work operations"
  )
  total_work_ops_bound <- ipw_owner_safe_product(
    n_links, stream_work_one,
    max_entries = max_wo,
    label = "total planned custom-IPW link work operations"
  )

  list(
    n_a_units = na,
    n_observed_risk_units = n_obs,
    p_covariates = p,
    m_cohorts = m,
    n_planned_links = n_links,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_link = logit_cert_count_output(n_chunks),
    observer_buffer_entries = observer_buffer_entries,
    reader_chunk_concurrent_entries = reader_chunk_concurrent_entries,
    borrowed_owner_entries = borrowed_entries,
    peak_concurrent_entries = peak_concurrent_entries,
    peak_concurrent_matrix_entries = peak_concurrent_matrix_entries,
    peak_concurrent_bytes_8double = 8 * peak_concurrent_entries,
    retained_entries_per_link = retained_entries_one,
    total_retained_entries = total_retained_entries,
    total_retained_bytes_8double = 8 * total_retained_entries,
    work_ops_per_link_bound = stream_work_one,
    total_work_ops_bound = total_work_ops_bound,
    expected_endpoints_read_per_link = logit_cert_count_output(2 * n_obs),
    peak_rss_promised = FALSE,
    memory_disclosure = paste0(
      "Deterministic preflight charges package-owned chunk buffers, the single temporary observer buffer ",
      "(y1, y0, D, covariates, i.weights), and the retained (obs_a_idx, psi_obs, b_rate) link representation ",
      "before invoking the user callback. Arbitrary workspace allocated inside the opaque user callback and ",
      "OS process RSS are outside package-owned array accounting."
    )
  )
}

fit_custom_ipw_link_view <- function(
  frame_owner,
  link_view,
  control,
  custom_est_method,
  custom_dots = list(),
  link_id = NULL,
  functional_key = NULL
) {
  ipw_point_view_require_deps()
  if (missing(frame_owner) || missing(link_view) || missing(control) || missing(custom_est_method)) {
    cdid_custom_ipw_abort(
      "All required arguments (`frame_owner`, `link_view`, `control`, `custom_est_method`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.function(custom_est_method)) {
    cdid_custom_ipw_abort("`custom_est_method` must be a function.", "cdid_ipw_point_view_input_error")
  }
  if (!is.list(custom_dots) || is.object(custom_dots)) {
    cdid_custom_ipw_abort("`custom_dots` must be a plain list.", "cdid_ipw_point_view_input_error")
  }

  # 1. Preflight strictly BEFORE any chunk query, buffer allocation, or callback invocation
  preflight <- preflight_custom_ipw_link_view(
    frame_owner = frame_owner,
    link_view = link_view,
    control = control,
    n_planned_links = 1L
  )
  ctrl <- validate_ipw_point_view_control(control)
  owner_ctrl <- ipw_point_view_to_owner_control(ctrl)

  eff_link_id <- if (!is.null(link_id)) as.character(link_id)[1L] else link_view$link_id
  eff_func_key <- if (!is.null(functional_key)) as.character(functional_key)[1L] else link_view$functional_key

  na <- as.integer(frame_owner$n_a_units)
  n_obs <- as.integer(link_view$support$n_observed_risk_units)
  p <- as.integer(frame_owner$n_covariates)
  m <- as.integer(link_view$cohorts$m_cohorts)
  step <- as.integer(preflight$effective_chunk_rows)

  obs_buf_env <- new.env(parent = emptyenv())
  obs_buf_env$active <- FALSE
  obs_buf_env$observer_released <- FALSE

  cleanup_temporary_buffers <- function() {
    for (nm in c("risk_a_idx", "risk_gi", "risk_w", "y0", "y1", "Z_use")) {
      if (exists(nm, envir = obs_buf_env, inherits = FALSE)) {
        assign(nm, NULL, envir = obs_buf_env)
        rm(list = nm, envir = obs_buf_env)
      }
    }
    obs_buf_env$active <- FALSE
    obs_buf_env$observer_released <- TRUE
  }
  on.exit(cleanup_temporary_buffers(), add = TRUE)

  # 2. Allocate single temporary observer buffer for this link
  obs_buf_env$active <- TRUE
  obs_buf_env$risk_a_idx <- integer(n_obs)
  obs_buf_env$risk_gi <- integer(n_obs)
  obs_buf_env$risk_w <- numeric(n_obs)
  obs_buf_env$y0 <- numeric(n_obs)
  obs_buf_env$y1 <- numeric(n_obs)
  obs_buf_env$Z_use <- matrix(
    0.0,
    nrow = n_obs,
    ncol = p,
    dimnames = list(NULL, colnames(frame_owner$Z_a))
  )

  # 3. Stream bounded link-view chunks in canonical increasing A-unit order
  write_pos <- 0L
  last_a_idx <- 0L
  endpoints_read_total <- 0L
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
      trace_env = y_trace
    )
    k_ch <- as.integer(ch$counts$observed_risk_row_count)
    ep_ch <- as.integer(ch$read_diagnostics$n_endpoints_read)
    endpoints_read_total <- endpoints_read_total + ep_ch

    if (!isTRUE(ch$read_diagnostics$risk_filtered_before_outcome_read) ||
        isTRUE(ch$read_diagnostics$off_risk_endpoints_read) ||
        isTRUE(ch$read_diagnostics$nonobserver_endpoints_read) ||
        !ipw_point_view_identical(ep_ch, as.integer(2L * k_ch)) ||
        !ipw_point_view_identical(length(y_trace$pos_s_index), k_ch) ||
        !ipw_point_view_identical(length(y_trace$pos_t_index), k_ch)) {
      cdid_custom_ipw_abort(
        "Link view query violated risk-first endpoint read filtering.",
        c("cdid_ipw_point_view_outcome_leak_error", "cdid_ipw_owner_integrity_error")
      )
    }
    ipw_owner_reset_trace(y_trace)

    if (k_ch > 0L) {
      if (write_pos + k_ch > n_obs) {
        cdid_custom_ipw_abort(
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
        cdid_custom_ipw_abort(
          "Malformed observed-risk slice in link view chunk.",
          c("cdid_ipw_point_view_chunk_error", "cdid_ipw_owner_integrity_error")
        )
      }
      if (a_ch[1L] <= last_a_idx || (k_ch > 1L && any(diff(a_ch) <= 0L)) ||
          a_ch[1L] < cur_s || a_ch[k_ch] > cur_e) {
        cdid_custom_ipw_abort(
          "Streamed risk-set observer `a_index` is not strictly increasing in global A-unit order.",
          c("cdid_ipw_point_view_order_error", "cdid_ipw_owner_order_mismatch_error")
        )
      }
      if (anyNA(gi_ch) || any(gi_ch < 1L | gi_ch > m)) {
        cdid_custom_ipw_abort(
          "Streamed risk-set observer cohort index is out of `[1, m_cohorts]`.",
          c("cdid_ipw_point_view_cohort_error", "cdid_ipw_owner_cohort_error")
        )
      }
      dy_ch <- yr_ch - yl_ch
      if (anyNA(yl_ch) || any(!is.finite(yl_ch)) || anyNA(yr_ch) || any(!is.finite(yr_ch)) ||
          anyNA(dy_ch) || any(!is.finite(dy_ch))) {
        rm(a_ch, gi_ch, loc_ch, w_ch, yl_ch, yr_ch, dy_ch, Z_ch)
        ch <- NULL
        cleanup_temporary_buffers()
        cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
      }

      dst <- seq.int(write_pos + 1L, write_pos + k_ch)
      obs_buf_env$risk_a_idx[dst] <- a_ch
      obs_buf_env$risk_gi[dst] <- gi_ch
      obs_buf_env$risk_w[dst] <- w_ch
      obs_buf_env$y0[dst] <- yl_ch
      obs_buf_env$y1[dst] <- yr_ch
      obs_buf_env$Z_use[dst, ] <- Z_ch
      write_pos <- write_pos + k_ch
      last_a_idx <- a_ch[k_ch]
      rm(a_ch, gi_ch, loc_ch, w_ch, yl_ch, yr_ch, dy_ch, Z_ch, dst)
    }
    ch <- NULL
    cur_s <- as.double(cur_e) + 1
  }

  if (write_pos != n_obs || !ipw_point_view_identical(endpoints_read_total, as.integer(2L * n_obs))) {
    cdid_custom_ipw_abort(
      "Streamed risk-set observer count or endpoint read count disagrees with `link_view$support$n_observed_risk_units`.",
      c("cdid_ipw_point_view_count_error", "cdid_ipw_owner_integrity_error")
    )
  }

  # 4. Verify cohort-specific pair counts/masses/rates and run shared observer fit worker
  sup <- list(
    n = as.double(na),
    g = link_view$cohorts$g,
    K = link_view$cohorts$K,
    cohorts = link_view$cohorts$cohorts,
    cohort_labels = link_view$cohorts$cohort_labels,
    m_cohorts = m,
    group_count = as.integer(unname(link_view$support$by_cohort_a_count)),
    group_a_count = as.integer(unname(link_view$support$by_cohort_a_count)),
    group_mass = as.double(unname(link_view$support$by_cohort_a_mass)),
    group_a_mass = as.double(unname(link_view$support$by_cohort_a_mass)),
    q = as.double(unname(link_view$support$q)),
    a_weight_total = sum(frame_owner$w_a)
  )

  w_out <- custom_ipw_observer_fit_worker(
    sup = sup,
    risk_a_idx = obs_buf_env$risk_a_idx,
    risk_gi = obs_buf_env$risk_gi,
    risk_w = obs_buf_env$risk_w,
    y1 = obs_buf_env$y1,
    y0 = obs_buf_env$y0,
    Z_use = obs_buf_env$Z_use,
    custom_est_method = custom_est_method,
    custom_dots = custom_dots,
    link_id = eff_link_id
  )

  if (!ipw_point_view_identical(as.integer(w_out$pair_count), as.integer(unname(link_view$support$by_cohort_pair_count))) ||
      !ipw_point_view_identical(as.double(w_out$pair_mass), as.double(unname(link_view$support$by_cohort_pair_mass))) ||
      !ipw_point_view_identical(as.double(w_out$rates), as.double(unname(link_view$support$rates)))) {
    cdid_custom_ipw_abort(
      "Reconstructed observer-buffer pair masses or rates disagree with `link_view$support`.",
      c("cdid_ipw_point_view_integrity_error", "cdid_ipw_shared_mask_error")
    )
  }

  obs_a_idx <- as.integer(obs_buf_env$risk_a_idx)
  psi_obs <- w_out$psi_obs
  b_rate <- w_out$b_rate
  att_hat <- w_out$estimate
  callback_if_mean <- w_out$callback_if_mean
  derivative_source <- w_out$derivative_source
  rm(w_out)

  # Release y1, y0, Z_use, risk_w, risk_gi immediately
  cleanup_temporary_buffers()

  cohort_a_mass <- as.double(link_view$support$by_cohort_a_mass)
  cohort_pair_mass <- as.double(link_view$support$by_cohort_pair_mass)
  q_vec <- as.double(link_view$support$q)
  rates_named <- as.double(link_view$support$rates)
  names(cohort_a_mass) <- names(cohort_pair_mass) <- names(q_vec) <- names(rates_named) <- sup$cohort_labels

  custom_link <- structure(
    list(
      link_id = eff_link_id,
      functional_key = eff_func_key,
      nuisance_id = paste0("custom_ipw:", eff_link_id),
      fitted_nuisance_ref = NULL,
      nuisance_view_ref = link_view$nuisance_ref,
      nuisance_ref = link_view$nuisance_ref,
      link_view_ref = link_view$ref,
      link_ref = link_view$ref,
      frame_ref = link_view$frame_ref,
      panel_ref = link_view$panel_ref,
      estimate = att_hat,
      pair = list(
        s = link_view$pair$s,
        t = link_view$pair$t,
        s_calendar_index = link_view$pair$s_calendar_index,
        t_calendar_index = link_view$pair$t_calendar_index,
        link_id = eff_link_id,
        functional_key = eff_func_key
      ),
      shared_across_dates = FALSE,
      nuisance_mask_verified = isTRUE(link_view$nuisance_mask_verified),
      cohorts = link_view$cohorts,
      dimensions = list(
        n_a_units = na,
        n_observed_risk_units = n_obs,
        p_covariates = p,
        m_cohorts = m
      ),
      counts = list(
        n_a_units = na,
        n_risk_a_units = as.integer(link_view$support$n_risk_a_units),
        n_outside_risk_a_units = as.integer(link_view$support$n_outside_risk_a_units),
        n_observed_pair_units = as.integer(link_view$support$n_observed_pair_all),
        n_observed_risk_units = n_obs,
        n_nonobserved_risk_units = as.integer(link_view$support$n_nonobserved_risk_units),
        n_treated_obs = as.integer(link_view$support$by_cohort_pair_count[1L]),
        n_control_obs = as.integer(sum(link_view$support$by_cohort_pair_count[-1L])),
        by_cohort_a = link_view$support$by_cohort_a_count,
        by_cohort_pair = link_view$support$by_cohort_pair_count
      ),
      masses = list(
        a_weight_total = sup$a_weight_total,
        cohort_a_weight = cohort_a_mass,
        cohort_pair_weight = cohort_pair_mass,
        q = q_vec,
        rates = rates_named
      ),
      obs_a_idx = obs_a_idx,
      psi_obs = psi_obs,
      b_rate = b_rate,
      diagnostics = list(callback_if_mean = callback_if_mean, derivative_source = derivative_source),
      preflight = preflight,
      accounting = list(
        retained_entries = preflight$retained_entries_per_link,
        peak_callback_input_entries = preflight$observer_buffer_entries,
        callback_internal_workspace = "untracked_opaque_callback",
        endpoints_read = endpoints_read_total,
        observer_released = TRUE
      )
    ),
    class = c("cdid_custom_ipw_link", "cdid_ipw_point_link_view")
  )
  ipw_seal(custom_link, "custom_ipw_link")
}

validate_custom_ipw_link <- function(link, frame_owner = NULL) {
  if (!inherits(link, "cdid_custom_ipw_link") || !is.list(link) || !ipw_seal_ok(link, "custom_ipw_link")) {
    cdid_custom_ipw_abort(
      "Expected a sealed `cdid_custom_ipw_link` object.",
      c("cdid_custom_ipw_seal_error", "cdid_ipw_point_view_input_error")
    )
  }
  if (!is.null(frame_owner) && !identical(link$frame_ref, frame_owner$ref)) {
    cdid_custom_ipw_abort(
      "`cdid_custom_ipw_link` does not belong to the supplied `frame_owner`.",
      "cdid_ipw_context_mismatch_error"
    )
  }
  n_obs <- as.integer(link$counts$n_observed_risk_units)
  m <- as.integer(link$cohorts$m_cohorts)
  if (!is.integer(link$obs_a_idx) || length(link$obs_a_idx) != n_obs ||
      !is.double(link$psi_obs) || length(link$psi_obs) != n_obs ||
      !is.double(link$b_rate) || length(link$b_rate) != m ||
      !is.numeric(link$estimate) || length(link$estimate) != 1L || !is.finite(link$estimate)) {
    cdid_custom_ipw_abort("Corrupted `cdid_custom_ipw_link` internal arrays.", "cdid_ipw_point_view_integrity_error")
  }
  invisible(TRUE)
}

build_custom_ipw_catalog_links <- function(
  metadata_obj,
  frame_owner = NULL,
  owner_control,
  point_control,
  custom_est_method,
  custom_dots = list(),
  retained_pipeline_entries = 0
) {
  validate_owner_comparison_catalog(metadata_obj, frame_owner = frame_owner)
  own_ctrl <- validate_ipw_owner_control(owner_control)
  pt_ctrl <- validate_ipw_point_view_control(point_control)
  if (!is.function(custom_est_method)) {
    cdid_custom_ipw_abort("`custom_est_method` must be a function.", "cdid_ipw_pipeline_input_error")
  }
  if (!is.list(custom_dots) || is.object(custom_dots)) {
    cdid_custom_ipw_abort("`custom_dots` must be a plain list.", "cdid_ipw_pipeline_input_error")
  }

  cat_rows <- metadata_obj$catalog$rows
  M <- nrow(cat_rows)
  if (M < 1L) {
    cdid_ipw_pipeline_abort(
      "Comparison catalog contains zero qualified links; cannot construct IPW views.",
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

  # 1. Construct lightweight link views and run whole-catalog preflight BEFORE any callback invocation
  row_specs <- vector("list", M)
  link_views <- vector("list", M)
  link_ids <- character(M)
  total_retained <- 0.0
  total_work <- 0.0
  max_link_workspace <- 0.0
  nuisance_views <- list()
  link_to_nuisance_idx <- integer(M)

  for (m in seq_len(M)) {
    pipeline_metadata_row_spec_preflight(metadata_obj, m, own_ctrl)
    spec <- extract_pipeline_catalog_row_spec(metadata_obj, m)
    nv_idx <- NA_integer_
    for (k in seq_along(nuisance_views)) {
      nv <- nuisance_views[[k]]
      if (!identical(nv$cohorts$g, spec$g) || !identical(nv$cohorts$K, spec$K)) next
      trial <- tryCatch(construct_ipw_link_view(
        f_owner, spec$s, spec$t, spec$g, spec$K, nuisance_view = nv,
        control = own_ctrl, link_id = paste0("link:r", spec$row_id),
        functional_key = spec$functional_key, require_support = TRUE,
        verify_nuisance_mask = TRUE), cdid_ipw_shared_mask_error = identity)
      if (!inherits(trial, "error")) { nv_idx <- k; break }
    }
    if (is.na(nv_idx)) {
      nv_idx <- length(nuisance_views) + 1L
      nuisance_views[[nv_idx]] <- construct_ipw_nuisance_view(
        f_owner, spec$g, spec$K, spec$s, spec$t, own_ctrl,
        nuisance_id = paste0("custom_ipw:mask", nv_idx))
    }
    lv <- construct_ipw_link_view(
      frame_owner = f_owner,
      s = spec$s,
      t = spec$t,
      g = spec$g,
      K = spec$K,
      nuisance_view = nuisance_views[[nv_idx]],
      control = own_ctrl,
      link_id = paste0("link:r", spec$row_id),
      functional_key = spec$functional_key,
      require_support = TRUE,
      verify_nuisance_mask = TRUE
    )
    pre_m <- preflight_custom_ipw_link_view(
      frame_owner = f_owner,
      link_view = lv,
      control = pt_ctrl,
      n_planned_links = 1L
    )
    total_retained <- ipw_owner_safe_sum(
      c(total_retained, pre_m$retained_entries_per_link),
      max_bound = min(pt_ctrl$max_retained_entries, pt_ctrl$max_output_entries, pt_ctrl$max_matrix_entries),
      label = "whole-catalog custom-IPW retained link entries"
    )
    total_work <- ipw_owner_safe_sum(
      c(total_work, pre_m$work_ops_per_link_bound),
      max_bound = pt_ctrl$max_work_ops,
      label = "whole-catalog custom-IPW link fitting work operations"
    )
    max_link_workspace <- max(max_link_workspace,
      pre_m$peak_concurrent_entries - pre_m$retained_entries_per_link)
    link_to_nuisance_idx[m] <- nv_idx
    row_specs[[m]] <- spec
    link_views[[m]] <- lv
    link_ids[m] <- lv$link_id
  }

  ipw_owner_safe_sum(
    c(retained_pipeline_entries, total_retained, max_link_workspace),
    max_bound = min(pt_ctrl$max_matrix_entries, own_ctrl$max_matrix_entries, own_ctrl$max_retained_entries),
    label = "whole-catalog custom-IPW simultaneous retained links and callback input workspace"
  )

  # 2. Fit each catalog link via the custom callback
  point_links <- vector("list", M)
  for (m in seq_len(M)) {
    point_links[[m]] <- fit_custom_ipw_link_view(
      frame_owner = f_owner,
      link_view = link_views[[m]],
      control = pt_ctrl,
      custom_est_method = custom_est_method,
      custom_dots = custom_dots
    )
  }
  names(point_links) <- link_ids
  names(link_views) <- link_ids

  list(
    nuisance_views = nuisance_views,
    fitted_views = list(),
    nuisance_keys = vapply(nuisance_views, `[[`, character(1L), "nuisance_id"),
    operator_ids = as.integer(cat_rows$operator_id),
    row_specs = row_specs,
    link_views = link_views,
    point_links = point_links,
    link_to_nuisance_idx = link_to_nuisance_idx,
    link_nuisance_keys = vapply(nuisance_views, `[[`, character(1L), "nuisance_id")[link_to_nuisance_idx],
    mask_splits = list(),
    custom_retained_entries = total_retained
  )
}

normalize_custom_ipw_links_view <- function(point_links, control) {
  ctrl <- validate_ipw_point_view_control(control)
  if (inherits(point_links, "cdid_custom_ipw_link")) {
    links <- list(point_links)
  } else if (is.list(point_links) && !is.object(point_links) && length(point_links) >= 1L) {
    links <- point_links
  } else {
    cdid_custom_ipw_abort(
      "`point_links` must be a sealed `cdid_custom_ipw_link` or a nonempty plain list of them.",
      "cdid_ipw_point_view_input_error"
    )
  }
  L <- length(links)
  validate_custom_ipw_link(links[[1L]])
  frame_owner <- get("frame_owner", envir = links[[1L]]$frame_ref, inherits = FALSE)
  panel_owner <- get("panel_owner", envir = links[[1L]]$panel_ref, inherits = FALSE)
  validate_ipw_frame_owner(frame_owner, panel_owner = panel_owner)

  link_ids <- character(L)
  existing_retained <- 0.0
  for (ell in seq_len(L)) {
    lk <- links[[ell]]
    validate_custom_ipw_link(lk, frame_owner = frame_owner)
    id_ell <- if (!is.null(lk$link_id) && nzchar(lk$link_id)) lk$link_id else paste0("link_", ell)
    link_ids[ell] <- id_ell
    existing_retained <- ipw_owner_safe_sum(
      c(existing_retained, lk$accounting$retained_entries),
      max_bound = min(ctrl$max_retained_entries, ctrl$max_matrix_entries),
      label = "retained custom-IPW link entries"
    )
  }
  if (anyDuplicated(link_ids)) {
    cdid_custom_ipw_abort("Duplicate `link_id` across supplied `cdid_custom_ipw_link` objects.", "cdid_ipw_point_view_input_error")
  }
  borrowed <- ipw_point_view_owner_entries(frame_owner, panel_owner)
  list(
    links = links,
    n_links = L,
    link_ids = link_ids,
    frame_owner = frame_owner,
    panel_owner = panel_owner,
    existing_retained_entries = existing_retained,
    borrowed_owner_entries = borrowed,
    normalization_work_bound = 16 * L,
    normalization_payload_entries = 4 * L
  )
}

preflight_custom_ipw_if_view <- function(point_links, control, chunk_rows = NULL, return_norm = FALSE) {
  ipw_point_view_require_deps()
  if (!gmm_scalar_logical(return_norm)) {
    cdid_custom_ipw_abort("`return_norm` must be a scalar logical.", "cdid_ipw_point_view_input_error")
  }
  ctrl <- validate_ipw_point_view_control(control)
  norm <- normalize_custom_ipw_links_view(point_links, control = ctrl)
  frame_owner <- norm$frame_owner
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
    cdid_custom_ipw_abort(
      paste0("Requested chunk rows (", req_ch, ") exceed `control$max_chunk_rows` (", max_ch, ")."),
      c("cdid_ipw_point_view_bound_error", "cdid_ipw_owner_bound_error")
    )
  }
  B_eff <- min(na, req_ch)
  n_chunks <- ceiling(na / B_eff)
  max_m <- max(vapply(norm$links, function(lk) as.integer(lk$cohorts$m_cohorts), integer(1L)))

  b_by_L <- ipw_owner_check_matrix_dims(
    B_eff, L,
    max_matrix_entries = min(max_me, max_out),
    label = "custom-IPW IF chunk B x L matrix"
  )
  chunk_output_entries <- ipw_owner_safe_sum(
    c(b_by_L, 6 * B_eff, L, 64),
    max_bound = min(max_me, max_out),
    label = "custom-IPW IF chunk output entries"
  )
  local_working_entries <- ipw_owner_safe_sum(
    c(16 * B_eff, 4 * max_m, 64),
    max_bound = max_me,
    label = "custom-IPW IF chunk working entries"
  )
  existing_retained <- norm$existing_retained_entries
  peak_concurrent_entries <- ipw_owner_safe_sum(
    c(existing_retained, 2 * chunk_output_entries, local_working_entries),
    max_bound = max_me,
    label = "simultaneous custom-IPW IF replay concurrent entries"
  )
  per_link_stream_work <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(na, 16, max_wo, "custom-IPW IF chunk evaluation work per link"), 64 * n_chunks),
    max_bound = max_wo,
    label = "single-link custom-IPW IF replay work"
  )
  total_work_ops_bound <- ipw_owner_safe_sum(
    c(ipw_owner_safe_product(L, per_link_stream_work, max_wo, "all-links custom-IPW IF replay work"), 256),
    max_bound = max_wo,
    label = "total custom-IPW IF replay work operations bound"
  )

  pre <- list(
    n_a_units = na,
    n_links = L,
    link_ids = norm$link_ids,
    p_covariates = p,
    m_cohorts_max = max_m,
    n_unique_fitted_nuisances = L,
    include_components = FALSE,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks = logit_cert_count_output(n_chunks),
    n_output_matrices_per_chunk = 1L,
    b_by_L_entries = b_by_L,
    chunk_output_matrix_entries = b_by_L,
    chunk_output_vector_entries = chunk_output_entries - b_by_L,
    chunk_output_entries = chunk_output_entries,
    reader_chunk_concurrent_entries = 0,
    local_link_working_entries = local_working_entries,
    one_nuisance_cache_entries = 0,
    existing_retained_entries = existing_retained,
    borrowed_owner_entries = norm$borrowed_owner_entries,
    failure_output_entries = 0,
    peak_concurrent_entries = peak_concurrent_entries,
    peak_concurrent_matrix_entries = b_by_L,
    peak_concurrent_bytes_8double = 8 * peak_concurrent_entries,
    per_link_stream_work_bound = per_link_stream_work,
    total_work_ops_bound = total_work_ops_bound,
    normalization_work_bound = norm$normalization_work_bound,
    normalization_payload_entries = norm$normalization_payload_entries,
    producer_workspace_entries = peak_concurrent_entries - existing_retained,
    consumer_retention_included = FALSE,
    total_expected_endpoints_read = 0L,
    full_n_by_L_retained = FALSE,
    per_link_dense_masks_retained = FALSE,
    peak_rss_promised = FALSE
  )
  if (isTRUE(return_norm)) {
    return(list(preflight = pre, norm = norm, ctrl = ctrl))
  }
  pre
}

custom_ipw_if_eval_chunk_internal <- function(
  norm,
  a_start,
  a_end,
  callback_schema = FALSE,
  batch_index = 1L
) {
  frame_owner <- norm$frame_owner
  panel_owner <- norm$panel_owner
  na <- frame_owner$n_a_units
  a_s <- ipw_owner_int_bound(a_start, "a_start", min_val = 1)
  a_e <- ipw_owner_int_bound(a_end, "a_end", min_val = a_s)
  if (a_e > na) {
    cdid_custom_ipw_abort("`a_end` exceeds `frame_owner$n_a_units`.", "cdid_ipw_owner_range_error")
  }
  b_len <- a_e - a_s + 1L
  L <- norm$n_links
  link_ids <- norm$link_ids

  a_seq <- seq.int(a_s, a_e)
  u_idx <- panel_owner$a_to_unit[a_seq]
  u_id <- panel_owner$a_id[a_seq]
  u_g <- panel_owner$a_g[a_seq]
  u_w <- as.double(frame_owner$w_a[a_seq])
  c_chunk <- as.integer(frame_owner$a_cohort_index[a_seq])
  n_a_num <- as.double(na)

  phi_mat <- matrix(0.0, nrow = b_len, ncol = L, dimnames = list(NULL, link_ids))

  for (ell in seq_len(L)) {
    lnk <- norm$links[[ell]]
    frozen <- if (is.null(norm$frozen)) NULL else norm$frozen[[ell]]
    gi_chunk <- ipw_owner_view_cohort_index(lnk$cohorts, c_chunk)
    ev <- custom_ipw_eval_link_chunk(
      a_s = a_s,
      a_e = a_e,
      u_w = u_w,
      gi_chunk = gi_chunk,
      obs_a_idx = lnk$obs_a_idx,
      psi_obs = lnk$psi_obs,
      b_rate = if (is.null(frozen)) as.double(lnk$b_rate) else frozen$b_rate,
      q_vec = if (is.null(frozen)) as.double(lnk$masses$q) else frozen$q,
      rates_vec = if (is.null(frozen)) as.double(lnk$masses$rates) else frozen$rates,
      n_a_num = n_a_num,
      m_obs_num = as.double(lnk$counts$n_observed_risk_units),
      include_components = FALSE
    )
    phi_mat[, ell] <- ev$phi
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
      components = NULL,
      read_diagnostics = list(
        n_endpoints_read = 0L,
        risk_filtered_before_outcome_read = TRUE,
        off_risk_endpoints_read = FALSE,
        nonobserver_endpoints_read = FALSE,
        nuisance_cache_hits = 0L
      ),
      work = list(
        actual_work_ops = b_len * L * 8,
        max_simultaneous_query_chunks = 0L,
        max_simultaneous_nuisance_caches = 0L
      )
    ),
    class = c("cdid_ipw_point_if_batch", "cdid_ipw_if_batch")
  )
}

replay_custom_ipw_if_view <- function(
  point_links,
  control,
  callback,
  callback_schema = FALSE,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control) || missing(callback)) {
    cdid_custom_ipw_abort(
      "All required arguments (`point_links`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.function(callback)) {
    cdid_custom_ipw_abort("`callback` must be a function.", "cdid_ipw_point_view_input_error")
  }
  pn <- preflight_custom_ipw_if_view(point_links = point_links, control = control, return_norm = TRUE)
  pre <- pn$preflight
  ctrl <- pn$ctrl
  norm <- pn$norm

  na <- as.integer(pre$n_a_units)
  L <- as.integer(pre$n_links)
  step <- as.integer(pre$effective_chunk_rows)

  emitted_units <- 0L
  batches <- 0L
  max_batch_units <- 0L

  a_start <- 1L
  while (a_start <= na) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(na), as.double(a_start) + step - 1))
    b_len <- a_end - a_start + 1L
    batches <- batches + 1L

    batch <- custom_ipw_if_eval_chunk_internal(
      norm = norm,
      a_start = a_start,
      a_end = a_end,
      callback_schema = callback_schema,
      batch_index = batches
    )
    callback(batch)
    batch <- NULL

    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != na) {
    cdid_custom_ipw_abort("Custom-IPW IF replay unit count mismatch.", "cdid_ipw_count_error")
  }

  invisible(list(
    n_a_units = emitted_units,
    n_links = L,
    link_ids = norm$link_ids,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * as.double(L),
    actual_work_ops = as.double(na) * as.double(L) * 8.0,
    preflight = pre
  ))
}

stream_custom_ipw_covariance_view <- function(
  point_links,
  control,
  trace_env = NULL
) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control)) {
    cdid_custom_ipw_abort(
      "Required arguments (`point_links`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  pn1 <- preflight_custom_ipw_if_view(point_links = point_links, control = control, return_norm = TRUE)
  pre1 <- pn1$preflight
  ctrl <- pn1$ctrl
  norm <- pn1$norm
  n_a <- as.double(pre1$n_a_units)
  L <- as.integer(pre1$n_links)
  link_ids <- norm$link_ids

  cov_entries <- ipw_owner_check_matrix_dims(
    L, L,
    max_matrix_entries = min(ctrl$max_matrix_entries, ctrl$max_output_entries),
    label = "streamed L x L custom-IPW IF covariance matrix"
  )
  ipw_owner_safe_sum(
    c(pre1$peak_concurrent_entries, 4 * cov_entries, 8 * L, 6 * pre1$b_by_L_entries, 128),
    ctrl$max_matrix_entries,
    "two-pass streamed custom-IPW covariance concurrent entries"
  )

  sum_phi <- numeric(L)
  replay_custom_ipw_if_view(
    point_links = norm$links,
    control = ctrl,
    callback_schema = TRUE,
    trace_env = trace_env,
    callback = function(batch) {
      sum_phi <<- sum_phi + colSums(batch$phi)
      if (any(!is.finite(sum_phi))) {
        cdid_abort("Nonfinite IF sum during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
      }
    }
  )
  mean_phi <- sum_phi / n_a
  names(mean_phi) <- link_ids

  cov_mat <- matrix(0.0, nrow = L, ncol = L, dimnames = list(link_ids, link_ids))
  replay_info <- replay_custom_ipw_if_view(
    point_links = norm$links,
    control = ctrl,
    callback_schema = TRUE,
    trace_env = trace_env,
    callback = function(batch) {
      centered <- sweep(batch$phi, 2L, mean_phi, FUN = "-")
      if (any(!is.finite(centered))) {
        cdid_abort("Nonfinite centered IF during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
      }
      cov_mat <<- cov_mat + crossprod(centered)
      if (any(!is.finite(cov_mat))) {
        cdid_abort("Nonfinite IF crossproduct or sum during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
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
    n_links = L,
    link_ids = link_ids,
    n_batches = replay_info$n_batches,
    max_batch_units = replay_info$max_batch_units,
    preflight = pre1
  )
}

custom_ipw_if_plan <- function(point_links, control) {
  ipw_point_view_require_deps()
  if (missing(point_links) || missing(control)) {
    cdid_custom_ipw_abort(
      "Required arguments (`point_links`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  pn <- preflight_custom_ipw_if_view(point_links = point_links, control = control, return_norm = TRUE)
  pre <- pn$preflight
  norm <- pn$norm
  ctrl <- pn$ctrl
  fo <- norm$frame_owner
  po <- norm$panel_owner
  L <- as.integer(norm$n_links)
  na <- as.integer(pre$n_a_units)
  step <- as.integer(pre$effective_chunk_rows)
  n_chunks <- as.integer(pre$n_chunks)

  # The adapter's normalized numeric cohort arrays are new live objects. Count
  # their payload, conversion copies, lists/seals and control comparison work
  # before allocating any adapter reference or normalized array.
  control_budget <- ipw_point_view_payload_budget(list(control, ctrl), list())
  cohort_entries <- ipw_owner_safe_sum(vapply(norm$links, function(lk) {
    ipw_owner_safe_product(3, lk$cohorts$m_cohorts,
      label = "custom IF frozen b_rate, q and rates")
  }, numeric(1L)), label = "custom IF frozen cohort arrays")
  record_entries <- ipw_owner_safe_sum(
    c(pre$normalization_payload_entries, 128 * as.double(L),
      2 * cohort_entries, control_budget$entries, 512),
    min(ctrl$max_retained_entries, ctrl$max_matrix_entries),
    "custom-IPW IF plan retained adapter entries"
  )
  ipw_owner_safe_sum(c(pre$peak_concurrent_entries, record_entries),
    ctrl$max_matrix_entries, "custom-IPW IF replay plus frozen adapter")
  one_verification <- ipw_owner_safe_sum(
    c(128 * as.double(L), 2 * control_budget$work, 512),
    ctrl$max_work_ops, "one custom-IPW IF plan verification work bound")
  verification_count <- 2 + 4 * (2 * as.double(n_chunks) + 1)
  verification_work <- ipw_owner_safe_product(verification_count, one_verification,
    ctrl$max_work_ops, "custom-IPW IF plan verification work across four passes")
  construction_work <- ipw_owner_safe_sum(
    c(pre$normalization_work_bound, 4 * cohort_entries, control_budget$work),
    ctrl$max_work_ops, "custom-IPW IF frozen adapter construction work")
  additional_work <- ipw_owner_safe_sum(c(construction_work, verification_work),
    ctrl$max_work_ops, "custom-IPW IF plan construction and verification work")
  ipw_owner_safe_sum(c(pre$total_work_ops_bound, additional_work),
    ctrl$max_work_ops, "custom-IPW IF producer plus frozen adapter work")

  link_views <- vector("list", L)
  nuisance_views <- vector("list", L)
  frozen <- vector("list", L)
  for (ell in seq_len(L)) {
    lk <- norm$links[[ell]]
    lv <- if (is.environment(lk$link_view_ref)) get0("link_view", envir = lk$link_view_ref, inherits = FALSE) else NULL
    if (!is.list(lv) || !ipw_owner_seal_ok(lv, "link_view") ||
        !identical(lk$link_ref, lk$link_view_ref) ||
        !identical(lv$frame_ref, fo$ref) || !identical(lv$panel_ref, fo$panel_ref)) {
      cdid_custom_ipw_abort(
        "Custom-IPW IF plan cannot bind a malformed or foreign link view.",
        c("cdid_ipw_point_view_integrity_error", "cdid_custom_ipw_error")
      )
    }
    link_views[[ell]] <- lv
    nv <- if (is.environment(lv$nuisance_ref)) get0("nuisance_view", lv$nuisance_ref, inherits = FALSE) else NULL
    if (!is.null(lv$nuisance_ref) && (!ipw_owner_seal_ok(nv, "nuisance_view") ||
        !identical(nv$frame_ref, fo$ref) || !identical(nv$panel_ref, fo$panel_ref))) {
      cdid_custom_ipw_abort("Custom-IPW IF plan cannot bind a foreign nuisance view.",
        "cdid_ipw_point_view_integrity_error")
    }
    nuisance_views[ell] <- list(nv)
    frozen[[ell]] <- ipw_owner_seal(list(
      source_ref = lk$ref,
      b_rate = as.double(lk$b_rate),
      q = as.double(lk$masses$q),
      rates = as.double(lk$masses$rates)
    ), "custom_ipw_if_frozen")
  }
  norm$frozen <- frozen

  plan <- list(
    kind = "cdid_custom_ipw_if_plan",
    control = ctrl,
    control_in = control,
    frame_owner = fo,
    panel_owner = po,
    point_links = norm$links,
    link_views = link_views,
    nuisance_views = nuisance_views,
    obs_a_idx_refs = lapply(norm$links, `[[`, "obs_a_idx"),
    psi_obs_refs = lapply(norm$links, `[[`, "psi_obs"),
    b_rate_refs = lapply(norm$links, `[[`, "b_rate"),
    link_ids = norm$link_ids,
    n_links = L,
    norm = norm,
    preflight = pre,
    n_a_units = na,
    step = step,
    include_components = FALSE,
    custom = list(
      route = "custom_if_plan",
      normalizations_performed = 1L,
      n_links = L,
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
      lifetime = paste0(
        "References to the existing sealed owners, link views and frozen custom link outputs persist from the ",
        "pipeline's producer preflight through the GMM projection and are dropped after the solve; user callbacks ",
        "are never re-invoked during IF replay."
      )
    )
  )
  ipw_owner_seal(plan, "custom_ipw_if_plan")
}

custom_ipw_if_plan_verify <- function(plan, control) {
  if (!is.list(plan) || !identical(plan$kind, "cdid_custom_ipw_if_plan") ||
      !ipw_owner_seal_ok(plan, "custom_ipw_if_plan")) {
    cdid_custom_ipw_abort(
      "Custom-IPW IF plan is malformed, foreign, or differs from its sealed snapshot.",
      c("cdid_ipw_point_view_integrity_error", "cdid_custom_ipw_error")
    )
  }
  if (!identical(control, plan$control_in) && !identical(control, plan$control)) {
    cdid_custom_ipw_abort(
      "Custom-IPW IF plan controls differ from the traversal controls; reuse is refused.",
      c("cdid_ipw_control_mismatch_error", "cdid_custom_ipw_error")
    )
  }
  norm <- plan$norm
  fo <- plan$frame_owner
  po <- plan$panel_owner
  L <- plan$n_links
  bound <- is.list(norm) && gmm_sexp_aliases(norm$frame_owner, fo) && gmm_sexp_aliases(norm$panel_owner, po) &&
    identical(as.integer(norm$n_links), L) && identical(norm$link_ids, plan$link_ids) &&
    length(plan$point_links) == L && length(plan$link_views) == L &&
    length(plan$nuisance_views) == L && length(norm$frozen) == L &&
    length(plan$obs_a_idx_refs) == L && length(plan$psi_obs_refs) == L && length(plan$b_rate_refs) == L &&
    ipw_owner_seal_ok(fo, "frame_owner") && identical(fo$panel_ref, po$ref) && ipw_owner_seal_ok(po, "panel_owner")
  for (ell in seq_len(L)) {
    if (!bound) break
    pl <- plan$point_links[[ell]]
    lv <- plan$link_views[[ell]]
    nv <- plan$nuisance_views[[ell]]
    fr <- norm$frozen[[ell]]
    bound <- ipw_owner_seal_ok(pl, "custom_ipw_link") &&
      ipw_owner_seal_ok(fr, "custom_ipw_if_frozen") && identical(fr$source_ref, pl$ref) &&
      gmm_sexp_aliases(norm$links[[ell]], pl) &&
      ipw_owner_seal_ok(lv, "link_view") &&
      identical(pl$link_view_ref, lv$ref) && identical(pl$link_ref, lv$ref) &&
      identical(pl$frame_ref, fo$ref) && identical(pl$panel_ref, fo$panel_ref) &&
      identical(lv$frame_ref, fo$ref) && identical(lv$panel_ref, fo$panel_ref) &&
      identical(pl$nuisance_view_ref, lv$nuisance_ref) && identical(pl$nuisance_ref, lv$nuisance_ref) &&
      (if (is.null(nv)) is.null(lv$nuisance_ref) else
        ipw_owner_seal_ok(nv, "nuisance_view") && identical(lv$nuisance_ref, nv$ref) &&
          identical(nv$frame_ref, fo$ref) && identical(nv$panel_ref, fo$panel_ref)) &&
      gmm_sexp_aliases(pl$obs_a_idx, plan$obs_a_idx_refs[[ell]]) &&
      gmm_sexp_aliases(pl$psi_obs, plan$psi_obs_refs[[ell]]) &&
      gmm_sexp_aliases(pl$b_rate, plan$b_rate_refs[[ell]]) &&
      (is.null(pl$link_id) || !nzchar(pl$link_id) || identical(pl$link_id, plan$link_ids[[ell]]))
  }
  if (!bound) {
    cdid_custom_ipw_abort(
      "Custom-IPW IF plan owner, link view, or frozen custom-link bindings changed; reuse is refused.",
      c("cdid_ipw_point_view_integrity_error", "cdid_custom_ipw_error")
    )
  }
  invisible(TRUE)
}

replay_custom_ipw_if_plan <- function(
  plan,
  control,
  callback,
  callback_schema = TRUE,
  trace_env = NULL,
  carrier = NULL
) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control) || missing(callback)) {
    cdid_custom_ipw_abort(
      "All required arguments (`plan`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.function(callback)) {
    cdid_custom_ipw_abort("`callback` must be a function.", "cdid_ipw_point_view_input_error")
  }
  if (!gmm_scalar_logical(callback_schema)) {
    cdid_custom_ipw_abort("`callback_schema` must be a scalar logical.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_custom_ipw_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  if (!is.null(carrier)) {
    if (!isTRUE(callback_schema)) {
      cdid_custom_ipw_abort(
        "Retained point IF carrier supports only the callback-schema IF replay.",
        "cdid_ipw_point_view_input_error"
      )
    }
    return(replay_ipw_point_if_retained(carrier, plan, control, callback, trace_env = trace_env))
  }
  custom_ipw_if_plan_verify(plan, control)
  verifications <- 1L
  pre <- plan$preflight
  norm <- plan$norm
  na <- as.integer(pre$n_a_units)
  L <- as.integer(pre$n_links)
  step <- as.integer(pre$effective_chunk_rows)

  emitted_units <- 0L
  batches <- 0L
  max_batch_units <- 0L

  a_start <- 1L
  while (a_start <= na) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(na), as.double(a_start) + step - 1))
    b_len <- a_end - a_start + 1L
    batches <- batches + 1L

    batch <- custom_ipw_if_eval_chunk_internal(
      norm = norm,
      a_start = a_start,
      a_end = a_end,
      callback_schema = callback_schema,
      batch_index = batches
    )
    custom_ipw_if_plan_verify(plan, control)
    callback(batch)
    batch <- NULL
    custom_ipw_if_plan_verify(plan, control)
    verifications <- verifications + 2L

    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    if (!is.null(trace_env)) {
      trace_env$batches_emitted <- batches
      trace_env$units_emitted <- emitted_units
    }
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != na) {
    cdid_custom_ipw_abort("Custom-IPW IF replay unit count mismatch.", "cdid_ipw_count_error")
  }

  invisible(list(
    n_a_units = emitted_units,
    n_links = L,
    link_ids = norm$link_ids,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * as.double(L),
    actual_work_ops = as.double(na) * as.double(L) * 8.0,
    normalizations_performed = 0L,
    verifications = verifications,
    served_from = "custom_chunk_replay",
    preflight = pre
  ))
}

stream_custom_ipw_covariance_plan <- function(
  plan,
  control,
  trace_env = NULL,
  carrier = NULL
) {
  ipw_point_view_require_deps()
  if (missing(plan) || missing(control)) {
    cdid_custom_ipw_abort(
      "Required arguments (`plan`, `control`) must be supplied explicitly.",
      "cdid_ipw_point_view_input_error"
    )
  }
  if (!is.null(trace_env) && !is.environment(trace_env)) {
    cdid_custom_ipw_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_point_view_input_error")
  }
  custom_ipw_if_plan_verify(plan, control)
  carrier_entries <- 0
  carrier_work <- 0
  if (!is.null(carrier)) {
    ipw_point_if_retained_check(carrier, plan, control, verify_plan = FALSE)
    if (!identical(carrier$state, "reserved")) {
      cdid_custom_ipw_abort(
        "Retained point IF carrier must be reserved and unfilled for the covariance stream.",
        "cdid_ipw_point_view_integrity_error"
      )
    }
    carrier_entries <- carrier$entries
    carrier_work <- carrier$work
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
    label = "streamed L x L custom-IPW IF covariance matrix"
  )
  ipw_owner_safe_sum(
    c(pre1$peak_concurrent_entries, 4 * cov_entries, 8 * L, 6 * pre1$b_by_L_entries,
      plan$custom$additional_live_entries, carrier_entries, 128),
    ctrl$max_matrix_entries,
    "two-pass streamed custom-IPW covariance concurrent entries plus plan record"
  )
  ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(2, pre1$total_work_ops_bound, ctrl$max_work_ops, "two-pass custom-IPW IF replay work"),
      ipw_owner_safe_product(1 + 2 * (2 * as.double(pre1$n_chunks) + 1),
        plan$custom$one_verification_work_ops_bound, ctrl$max_work_ops,
        "custom-IPW covariance entry and two replay record verifications"),
      ipw_owner_safe_product(n_a, 2 * L * L + 4 * L, ctrl$max_work_ops, "centering and crossprod accumulation work"),
      carrier_work,
      128
    ),
    ctrl$max_work_ops,
    "total two-pass streamed custom-IPW covariance work"
  )

  sum_phi <- numeric(L)
  mean_callback <- function(batch) {
    sum_phi <<- sum_phi + colSums(batch$phi)
    if (any(!is.finite(sum_phi))) {
      cdid_abort("Nonfinite IF sum during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
    }
  }
  if (is.null(carrier)) {
    replay_custom_ipw_if_plan(
      plan = plan,
      control = control,
      callback_schema = TRUE,
      trace_env = trace_env,
      callback = mean_callback
    )
  } else {
    ipw_point_if_retained_fill(
      carrier = carrier,
      plan = plan,
      control = control,
      callback = mean_callback,
      trace_env = trace_env,
      .verified_entry = TRUE
    )
  }
  mean_phi <- sum_phi / n_a
  names(mean_phi) <- link_ids

  cov_mat <- matrix(0.0, nrow = L, ncol = L, dimnames = list(link_ids, link_ids))
  gram_callback <- function(batch) {
    centered <- sweep(batch$phi, 2L, mean_phi, FUN = "-")
    if (any(!is.finite(centered))) {
      cdid_abort("Nonfinite centered IF during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
    }
    cov_mat <<- cov_mat + crossprod(centered)
    if (any(!is.finite(cov_mat))) {
      cdid_abort("Nonfinite IF crossproduct or sum during custom-IPW covariance replay.", "cdid_ipw_numeric_error")
    }
    rm(centered)
  }
  replay_info <- replay_custom_ipw_if_plan(
    plan = plan,
    control = control,
    callback_schema = TRUE,
    trace_env = trace_env,
    carrier = carrier,
    callback = gram_callback
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
    n_links = L,
    link_ids = link_ids,
    n_batches = replay_info$n_batches,
    max_batch_units = replay_info$max_batch_units,
    preflight = pre1
  )
}
