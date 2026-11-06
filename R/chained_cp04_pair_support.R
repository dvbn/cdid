# Invocation-local, outcome-free CP04 support reuse across distinct g/K. The
# original owner scanner is unchanged and remains the declined-cache reference.
# Complete pairs keep their existing class/mask cache. Only the public registry
# opts in; no owner control, row batch, estimator or support rule is changed.

cp04_pair_support_context <- function(ctx, frame, n_views, retained, control, max_cohorts) {
  if (!is.environment(ctx) || !ipw_owner_construction_context_for(ctx, frame) ||
      !ipw_owner_seal_ok(frame, "frame_owner")) return(ctx)
  # Scalar date facts avoid allocating a calendar-sized comparison before
  # admission. Entirely complete panels retain the existing context unchanged.
  d <- 1L; n_dates <- length(frame$postings$date_count)
  while (d <= n_dates && frame$postings$date_count[d] == frame$n_a_units) d <- d + 1L
  if (d > n_dates) return(ctx)
  C <- length(frame$a_cohorts)
  # One record per view is a safe bound without reading/validating any later
  # descriptor. Each holds O(C) raw count/mass vectors, seal/descriptor wrappers
  # and bounded keys. A further block covers the all-cohort scan's O(C)
  # accumulator copies. Chunk scratch is the original constructor's 64*B.
  entries <- (256 + 64 * as.double(C)) * (as.double(n_views) + 1)
  total <- entries + ctx$link_entries
  m <- max(2, max_cohorts)
  constructor_retained <- 20 * m + frame$n_covariates + 256
  constructor_peak <- constructor_retained + 64 * min(frame$n_a_units, control$max_chunk_rows) + 6 * m + 256
  peak <- total + max(retained, constructor_peak)
  kept <- total + max(retained, constructor_retained)
  # Cache initialization/key work fits inside the unchanged constructor work
  # bounds or the optimization declines. It never obtains another full cap.
  initialization_work <- 64 * (as.double(C) + as.double(n_views) + as.double(d) + 1)
  if (!is.finite(peak) || !is.finite(kept) || max(peak, kept) > 2^53 - 1 ||
      peak > control$max_matrix_entries || kept > control$max_retained_entries ||
      initialization_work > control$max_work_ops) return(ctx)
  cache <- new.env(parent = emptyenv())
  cache$frame_ref <- frame$ref; cache$panel_ref <- frame$panel_ref
  cache$records <- new.env(parent = emptyenv())
  cache$capacity <- as.double(n_views); cache$n_records <- 0
  cache$entries <- entries; cache$scans <- 0; cache$reuses <- 0
  ctx$pair_support <- cache
  # Existing per-view composition guards include the whole extra reservation.
  ctx$link_entries <- total
  ctx
}

# Exact raw counts/masses in original A-row chunks. Each cohort's within-chunk
# weights keep their original order, as in ipw_owner_chunked_pair_support; only
# its owner cohort index replaces the requested-cohort binary lookup/filter.
cp04_pair_support_scan <- function(frame, pair, step) {
  C <- length(frame$a_cohorts); n_a <- frame$n_a_units
  pair_count <- integer(C); pair_mass <- numeric(C)
  n_observed <- 0L; chunks <- 0; ops <- 0
  ds <- pair$s_calendar_index; dt <- pair$t_calendar_index
  if (!is.na(ds) && !is.na(dt) && frame$postings$date_count[ds] > 0L && frame$postings$date_count[dt] > 0L) {
    first <- 1L
    while (first <= n_a) {
      ipw_owner_check_interrupt()
      last <- as.integer(min(as.double(n_a), as.double(first) + step - 1))
      hit <- ipw_owner_query_pair_mask(frame, pair$s, pair$t, a_start = first, a_end = last, ds = ds, dt = dt)
      nh <- length(hit$a_index)
      if (nh > 0L) {
        n_observed <- n_observed + nh
        gi <- frame$a_cohort_index[hit$a_index]; w <- frame$w_a[hit$a_index]
        group_order <- order(gi, method = "radix")
        runs <- rle(gi[group_order])
        group_end <- cumsum(runs$lengths)
        group_start <- group_end - runs$lengths + 1L
        for (z in seq_along(runs$values)) {
          j <- runs$values[z]
          pair_count[j] <- pair_count[j] + runs$lengths[z]
          pair_mass[j] <- pair_mass[j] + sum(w[group_order[seq.int(group_start[z], group_end[z])]])
        }
      }
      chunks <- chunks + 1
      ops <- ops + 4 * (last - first + 1L)
      first <- as.double(last) + 1
    }
  } else {
    chunks <- ceiling(n_a / step)
    ops <- n_a
  }
  ipw_owner_seal(list(frame_ref = frame$ref, panel_ref = frame$panel_ref, pair = pair, step = step,
    pair_count = pair_count, pair_mass = pair_mass, n_observed = n_observed, chunks = chunks, ops = ops),
    "cp04_pair_support")
}

# Apply the original scanner's support/rate rules to ONLY the requested ordered
# descriptor. An absent owner cohort remains zero under require_support=FALSE.
cp04_pair_support_subset <- function(raw, frame, cohorts, require_support) {
  m <- cohorts$m_cohorts
  pair_count <- integer(m); pair_mass <- numeric(m)
  present <- which(!is.na(cohorts$a_cohort_map))
  pair_count[present] <- raw$pair_count[cohorts$a_cohort_map[present]]
  pair_mass[present] <- raw$pair_mass[cohorts$a_cohort_map[present]]
  n_obs_pair_all <- raw$n_observed; ops <- raw$ops
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
  rates <- ifelse(cohorts$group_a_mass > 0, pair_mass / cohorts$group_a_mass, NA_real_)
  if (isTRUE(require_support)) {
    rates[pair_count == cohorts$group_a_count] <- 1
    ops <- ops + 4 * m
  }
  if (isTRUE(require_support) && any(!is.finite(rates) | rates <= 0 | rates > 1)) {
    cdid_ipw_owner_abort(
      "Invalid marginal pair observation rates `r_j`.",
      c("cdid_ipw_owner_pair_support_error", "cdid_ipw_rate_error")
    )
  }
  names(pair_count) <- names(pair_mass) <- names(rates) <- cohorts$cohort_labels
  list(n_observed_pair_all = n_obs_pair_all, n_risk_a_units = sum(cohorts$group_a_count),
    n_outside_risk_a_units = frame$n_a_units - sum(cohorts$group_a_count),
    n_observed_risk_units = n_obs_risk, n_nonobserved_risk_units = sum(cohorts$group_a_count) - n_obs_risk,
    by_cohort_pair_count = pair_count, by_cohort_pair_mass = pair_mass, rates = rates,
    chunks = raw$chunks, ops = ops)
}

cp04_pair_support <- function(ctx, frame, pair, cohorts, step, require_support) {
  cache <- ctx$pair_support
  if (!is.environment(cache) || !ipw_owner_construction_context_for(ctx, frame) ||
      !ipw_owner_seal_ok(frame, "frame_owner") ||
      !identical(cache$frame_ref, frame$ref) || !identical(cache$panel_ref, frame$panel_ref) ||
      !is.environment(cache$records) || !gmm_scalar_numeric(cache$capacity) ||
      !gmm_scalar_numeric(cache$n_records, nonneg = TRUE) ||
      cache$n_records > cache$capacity || ipw_owner_pair_counts_full(frame, pair)) return(NULL)
  key <- paste0("cp04-pair:", step, ":", pair$s_calendar_index, ":", pair$t_calendar_index)
  raw <- get0(key, envir = cache$records, inherits = FALSE)
  scanned <- is.null(raw)
  if (scanned) {
    if (cache$n_records >= cache$capacity) return(NULL)
    raw <- cp04_pair_support_scan(frame, pair, step)
    assign(key, raw, envir = cache$records)
    cache$n_records <- cache$n_records + 1
    cache$scans <- cache$scans + 1
    ctx$link_scans <- ctx$link_scans + 1
  } else {
    # A key is only a candidate. Its exact immutable payload, owner identities,
    # full pair descriptor and step must agree; stale/colliding records decline.
    if (!is.list(raw) || !ipw_owner_seal_ok(raw, "cp04_pair_support") ||
        !identical(raw$frame_ref, frame$ref) || !identical(raw$panel_ref, frame$panel_ref) ||
        !identical(raw$pair, pair) || !identical(raw$step, step)) return(NULL)
    cache$reuses <- cache$reuses + 1
    ctx$link_reuses <- ctx$link_reuses + 1
  }
  list(support = cp04_pair_support_subset(raw, frame, cohorts, require_support), scanned = scanned)
}
