# Private CP02 fixed-operator weighted-pair IPW kernel.
# Computes a declared empirical link and its complete case-weight derivative
# using one reusable calendar-posting context, compact fit caches, and batched
# IF replay. No moment catalog, public estimator, or survey-weight regime is
# defined here.

ipw_scalar_numeric <- function(x, finite = TRUE) {
  is.atomic(x) && is.null(dim(x)) && !is.object(x) && is.numeric(x) &&
    length(x) == 1L && !is.na(x) && (!finite || is.finite(x))
}

ipw_bound <- function(x, name, allow_zero = FALSE) {
  if (!ipw_scalar_numeric(x, finite = TRUE) || x != floor(x) ||
      x < (if (allow_zero) 0 else 1) || x > 2^53 - 1) {
    cdid_abort(
      paste0("`", name, "` must be an explicit finite integer bound."),
      c("cdid_ipw_bound_error", "cdid_pairs_bound_error")
    )
  }
  as.double(x)
}

validate_ipw_glm_control <- function(control) {
  if (missing(control) || is.null(control) || !is.list(control)) {
    cdid_abort(
      "Supply an explicit `glm.control` list; no default solver policy is adopted.",
      c("cdid_ipw_control_error", "cdid_ipw_input_error")
    )
  }
  req <- c("epsilon", "maxit", "trace")
  if (length(control) != length(req) || anyDuplicated(names(control)) ||
      !setequal(names(control), req) ||
      !ipw_scalar_numeric(control$epsilon, finite = TRUE) || control$epsilon <= 0 ||
      !ipw_scalar_numeric(control$maxit, finite = TRUE) || control$maxit < 1 ||
      control$maxit != floor(control$maxit) || control$maxit > .Machine$integer.max ||
      !is.atomic(control$trace) || !is.null(dim(control$trace)) ||
      is.object(control$trace) || !is.logical(control$trace) ||
      length(control$trace) != 1L || is.na(control$trace)) {
    cdid_abort(
      "Invalid `glm.control` specification; require exactly unambiguous `epsilon`, `maxit`, `trace` fields, finite positive `epsilon`, representable integer `maxit >= 1`, and logical `trace`.",
      c("cdid_ipw_control_error", "cdid_ipw_input_error")
    )
  }
  list(
    epsilon = as.double(control$epsilon),
    maxit = as.integer(control$maxit),
    trace = as.logical(control$trace)
  )
}

validate_ipw_context <- function(context) {
  if (!inherits(context, "cdid_ipw_context") || !is.list(context)) {
    cdid_abort("Expected `cdid_ipw_context`.", c("cdid_ipw_context_error", "cdid_ipw_input_error"))
  }
  ref <- context$ref
  if (!is.environment(ref) || !environmentIsLocked(ref) ||
      !exists("context", envir = ref, inherits = FALSE) ||
      !bindingIsLocked("context", ref) || !identical(context, ref$context)) {
    cdid_abort(
      "IPW context differs from its frozen snapshot; reconstruct the context and refit links after changing inputs.",
      c("cdid_ipw_context_mismatch_error", "cdid_ipw_context_error")
    )
  }
  # Unchanged lists and their arrays share the snapshot's R objects. `identical`
  # checks equality without copying inputs or retaining per-link unit arrays.
  invisible(context)
}

validate_ipw_components <- function(include_components) {
  if (!is.atomic(include_components) || !is.null(dim(include_components)) ||
      is.object(include_components) || !is.logical(include_components) ||
      length(include_components) != 1L || is.na(include_components)) {
    cdid_abort("`include_components` must be a scalar logical.", c("cdid_ipw_input_error"))
  }
  invisible(include_components)
}

# First index in vec[left..right] with vec[idx] >= target, or right + 1L.
lower_bound_int <- function(vec, left, right, target) {
  lo <- as.integer(left)
  hi <- as.integer(right)
  ans <- hi + 1L
  tval <- as.integer(target)
  while (lo <= hi) {
    mid <- lo + (hi - lo) %/% 2L
    if (vec[mid] >= tval) {
      ans <- mid
      hi <- mid - 1L
    } else {
      lo <- mid + 1L
    }
  }
  ans
}

# Last index in vec[left..right] with vec[idx] <= target, or left - 1L.
upper_bound_int <- function(vec, left, right, target) {
  lo <- as.integer(left)
  hi <- as.integer(right)
  ans <- lo - 1L
  tval <- as.integer(target)
  while (lo <= hi) {
    mid <- lo + (hi - lo) %/% 2L
    if (vec[mid] <= tval) {
      ans <- mid
      lo <- mid + 1L
    } else {
      hi <- mid - 1L
    }
  }
  ans
}

make_date_pair_descriptor <- function(s, t, context = NULL) {
  if (!ipw_scalar_numeric(s, finite = TRUE) || !ipw_scalar_numeric(t, finite = TRUE)) {
    cdid_abort(
      "Date-pair endpoints `s` and `t` must be finite numeric scalars.",
      c("cdid_ipw_pair_error", "cdid_ipw_input_error")
    )
  }
  s_num <- as.double(s)
  t_num <- as.double(t)
  if (!(s_num < t_num)) {
    cdid_abort(
      "Date-pair descriptor requires strict chronological order `s < t`.",
      c("cdid_ipw_pair_error", "cdid_ipw_input_error")
    )
  }
  ds <- dt <- NA_integer_
  if (!is.null(context)) {
    validate_ipw_context(context)
    ds <- match(s_num, context$all_times)
    dt <- match(t_num, context$all_times)
  }
  structure(
    list(s = s_num, t = t_num, s_calendar_index = ds, t_calendar_index = dt),
    class = "cdid_ipw_pair_descriptor"
  )
}

construct_ipw_context <- function(pair_index, id, Z, weights = NULL) {
  if (!inherits(pair_index, "cdid_observed_pair_index") ||
      !is.data.frame(pair_index$units) || !is.data.frame(pair_index$observations)) {
    cdid_abort(
      "Expected compact `cdid_observed_pair_index` from `construct_observed_pairs()`.",
      c("cdid_ipw_context_error", "cdid_ipw_input_error")
    )
  }
  units <- pair_index$units
  obs <- pair_index$observations
  nu <- nrow(units)
  a_unit_index <- which(units$pair_eligible)
  n_a <- length(a_unit_index)
  if (n_a < 1L) {
    cdid_abort(
      "At least one pair-eligible A-unit (`n_obs >= 2`) is required.",
      c("cdid_ipw_context_error", "cdid_ipw_support_error")
    )
  }

  if (missing(id) || is.null(id)) {
    cdid_abort(
      "Explicit unit `id` vector is required to align covariates and weights.",
      c("cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }
  if (is.factor(id)) id <- as.character(id)
  if (!is.atomic(id) || !is.null(dim(id)) || is.object(id) || length(id) < 1L || anyNA(id)) {
    cdid_abort(
      "`id` must be a nonmissing scalar vector of unit identifiers.",
      c("cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }
  id_is_num <- is.numeric(id)
  id_is_chr <- is.character(id)
  idx_is_num <- is.numeric(units$id)
  idx_is_chr <- is.character(units$id)
  if (!(id_is_num || id_is_chr) ||
      (id_is_num && (!idx_is_num || any(!is.finite(id)))) ||
      (id_is_chr && (!idx_is_chr || any(!nzchar(id))))) {
    cdid_abort(
      "`id` must match the exact type (finite numeric or nonempty character) of `pair_index$units$id` without string coercion.",
      c("cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }

  # Match numeric identifiers directly as IEEE doubles, never via formatted strings.
  id_key <- if (id_is_num) as.double(id) else id
  unit_key <- if (idx_is_num) as.double(units$id) else units$id
  if (anyDuplicated(id_key)) {
    cdid_abort(
      "Duplicate unit identifiers in `id` are invalid.",
      c("cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }
  unit_pos_of_supplied <- match(id_key, unit_key)
  if (anyNA(unit_pos_of_supplied)) {
    cdid_abort(
      "`id` contains unit identifiers absent from `pair_index$units$id`.",
      c("cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }
  a_pos_in_supplied <- match(unit_key[a_unit_index], id_key)
  if (anyNA(a_pos_in_supplied)) {
    cdid_abort(
      "Every pair-eligible A-unit (`n_obs >= 2`) must be present in `id`.",
      c("cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }

  if (missing(Z) || !is.matrix(Z) || !is.numeric(Z) || is.object(Z) ||
      nrow(Z) != length(id) || ncol(Z) < 1L || anyNA(Z) || any(!is.finite(Z))) {
    cdid_abort(
      "`Z` must be a finite numeric matrix with `nrow(Z) == length(id)` and `ncol(Z) >= 1`.",
      c("cdid_ipw_covariate_error", "cdid_ipw_context_error")
    )
  }
  Z_a <- Z[a_pos_in_supplied, , drop = FALSE]
  storage.mode(Z_a) <- "double"
  rownames(Z_a) <- NULL

  weights_supplied <- !is.null(weights)
  if (!weights_supplied) {
    w_all <- rep(1.0, length(id))
  } else {
    if (!is.atomic(weights) || !is.null(dim(weights)) || is.object(weights) ||
        !is.numeric(weights) || length(weights) != length(id) ||
        anyNA(weights) || any(!is.finite(weights) | weights <= 0)) {
      cdid_abort(
        "`weights` must be a finite, strictly positive numeric vector of length `length(id)`.",
        c("cdid_ipw_weight_error", "cdid_ipw_context_error")
      )
    }
    w_all <- as.double(weights)
  }
  w_a <- as.double(w_all[a_pos_in_supplied])
  sum_w_a <- sum(w_a)
  if (!is.finite(sum_w_a) || sum_w_a <= 0) {
    cdid_abort(
      "Total A-unit mechanical case weight must be finite and positive.",
      c("cdid_ipw_weight_error", "cdid_ipw_context_error")
    )
  }

  # Build compact calendar posting lists over A-units once: total O(Q + U + T).
  # Because observations in `pair_index$observations` are already ordered by
  # `(unit_index, time)` and each unit has at most one observation per date,
  # a stable radix sort by `calendar_index` produces strictly sorted unit lists
  # for every calendar date without allocating a U * T matrix.
  all_times <- pair_index$all_times
  n_cal <- length(all_times)
  unit_to_a <- integer(nu)
  unit_to_a[a_unit_index] <- seq_len(n_a)
  obs_unit_index <- rep.int(seq_len(nu), units$n_obs)
  obs_a_index <- unit_to_a[obs_unit_index]
  keep_a_obs <- which(obs_a_index > 0L)
  q_a <- length(keep_a_obs)
  cal_a <- obs$calendar_index[keep_a_obs]
  ord_cal <- order(cal_a, method = "radix")
  date_count <- tabulate(cal_a, nbins = n_cal)
  date_end <- cumsum(date_count)
  date_start <- date_end - date_count + 1L
  post_a_index <- obs_a_index[keep_a_obs[ord_cal]]
  post_unit_index <- obs_unit_index[keep_a_obs[ord_cal]]
  post_obs_index <- keep_a_obs[ord_cal]

  a_id <- units$id[a_unit_index]
  a_g <- as.double(units$g[a_unit_index])
  a_cohorts <- sort(unique(a_g))
  a_cohort_index <- match(a_g, a_cohorts)
  cohort_a_count <- tabulate(a_cohort_index, nbins = length(a_cohorts))
  cohort_a_weight <- vapply(
    seq_along(a_cohorts),
    function(j) sum(w_a[a_cohort_index == j]),
    numeric(1L)
  )
  if (any(!is.finite(cohort_a_weight) | cohort_a_weight <= 0)) {
    cdid_abort(
      "Cohort A-unit weight sums must be finite and positive.",
      c("cdid_ipw_weight_error", "cdid_ipw_context_error")
    )
  }
  cohort_a_share <- cohort_a_weight / n_a

  # Preserve singleton and unobserved unit metadata explicitly while keeping
  # them excluded from A-unit rate denominators.
  unit_supplied_pos <- match(unit_key, id_key)
  unit_weights_all <- rep(NA_real_, nu)
  present_units <- which(!is.na(unit_supplied_pos))
  unit_weights_all[present_units] <- w_all[unit_supplied_pos[present_units]]

  ref <- new.env(parent = emptyenv())
  ctx <- structure(
    list(
      n_a_units = n_a,
      n_total_units = nu,
      n_singleton_units = pair_index$counts$n_singleton_units,
      n_unobserved_units = pair_index$counts$n_unobserved_units,
      n_covariates = ncol(Z_a),
      n_calendar_dates = n_cal,
      n_observed_dates = length(pair_index$times),
      n_observed_rows = nrow(obs),
      n_a_observed_rows = q_a,
      a_unit_index = a_unit_index,
      unit_to_a = unit_to_a,
      a_id = a_id,
      a_g = a_g,
      a_cohorts = a_cohorts,
      a_cohort_index = a_cohort_index,
      cohort_a_count = cohort_a_count,
      cohort_a_weight = cohort_a_weight,
      cohort_a_share = cohort_a_share,
      Z_a = Z_a,
      w_a = w_a,
      sum_w_a = sum_w_a,
      weights_supplied = weights_supplied,
      all_times = all_times,
      times = pair_index$times,
      cohorts = pair_index$cohorts,
      observations = obs,
      units = units,
      singletons = pair_index$singletons,
      unit_weights_all = unit_weights_all,
      counts = pair_index$counts,
      cols = pair_index$cols,
      postings = list(
        date_start = date_start,
        date_end = date_end,
        date_count = date_count,
        post_a_index = post_a_index,
        post_unit_index = post_unit_index,
        post_obs_index = post_obs_index
      ),
      ref = ref
    ),
    class = "cdid_ipw_context"
  )
  ref$context <- ctx
  lockEnvironment(ref, bindings = TRUE)
  ctx
}

build_ipw_context <- construct_ipw_context
chained_ipw_context <- construct_ipw_context

get_calendar_posting <- function(context, time) {
  validate_ipw_context(context)
  if (!ipw_scalar_numeric(time, finite = TRUE)) {
    cdid_abort("Calendar `time` must be a finite numeric scalar.", c("cdid_ipw_pair_error", "cdid_ipw_input_error"))
  }
  d <- match(as.double(time), context$all_times)
  post <- context$postings
  if (is.na(d) || post$date_count[d] == 0L) {
    return(data.frame(
      a_index = integer(0L),
      unit_index = integer(0L),
      id = context$a_id[integer(0L)],
      obs_index = integer(0L),
      stringsAsFactors = FALSE
    ))
  }
  idx <- seq.int(post$date_start[d], post$date_end[d])
  a_idx <- post$post_a_index[idx]
  data.frame(
    a_index = a_idx,
    unit_index = post$post_unit_index[idx],
    id = context$a_id[a_idx],
    obs_index = post$post_obs_index[idx],
    stringsAsFactors = FALSE
  )
}

# Internal fast posting-list intersection for date pair (s, t), optionally
# restricted to a contiguous A-unit range [a_start, a_end] by binary bounds.
query_observed_pair_internal <- function(context, s, t,
                                         a_start = 1L,
                                         a_end = context$n_a_units,
                                         ds = NULL, dt = NULL) {
  if (is.null(ds)) ds <- match(as.double(s), context$all_times)
  if (is.null(dt)) dt <- match(as.double(t), context$all_times)
  empty <- list(
    a_index = integer(0L),
    unit_index = integer(0L),
    obs_s = integer(0L),
    obs_t = integer(0L),
    y_s = numeric(0L),
    y_t = numeric(0L),
    dy = numeric(0L)
  )
  if (is.na(ds) || is.na(dt) || a_start > a_end) return(empty)
  post <- context$postings
  if (post$date_count[ds] == 0L || post$date_count[dt] == 0L) return(empty)

  pa <- post$post_a_index
  lo_s <- lower_bound_int(pa, post$date_start[ds], post$date_end[ds], a_start)
  hi_s <- upper_bound_int(pa, post$date_start[ds], post$date_end[ds], a_end)
  if (lo_s > hi_s) return(empty)

  lo_t <- lower_bound_int(pa, post$date_start[dt], post$date_end[dt], a_start)
  hi_t <- upper_bound_int(pa, post$date_start[dt], post$date_end[dt], a_end)
  if (lo_t > hi_t) return(empty)

  # Tighten both intervals to the shared endpoint range via binary bounds.
  min_a <- max(pa[lo_s], pa[lo_t])
  max_a <- min(pa[hi_s], pa[hi_t])
  if (min_a > max_a) return(empty)
  if (pa[lo_s] < min_a) lo_s <- lower_bound_int(pa, lo_s, hi_s, min_a)
  if (pa[hi_s] > max_a) hi_s <- upper_bound_int(pa, lo_s, hi_s, max_a)
  if (lo_s > hi_s) return(empty)
  if (pa[lo_t] < min_a) lo_t <- lower_bound_int(pa, lo_t, hi_t, min_a)
  if (pa[hi_t] > max_a) hi_t <- upper_bound_int(pa, lo_t, hi_t, max_a)
  if (lo_t > hi_t) return(empty)

  idx_s <- seq.int(lo_s, hi_s)
  idx_t <- seq.int(lo_t, hi_t)
  u_s <- pa[idx_s]
  u_t <- pa[idx_t]

  if (length(u_s) <= length(u_t)) {
    pos <- findInterval(u_s, u_t)
    hit <- which(pos > 0L & u_t[pmax(pos, 1L)] == u_s)
    if (length(hit) == 0L) return(empty)
    hit_s <- idx_s[hit]
    hit_t <- idx_t[pos[hit]]
  } else {
    pos <- findInterval(u_t, u_s)
    hit <- which(pos > 0L & u_s[pmax(pos, 1L)] == u_t)
    if (length(hit) == 0L) return(empty)
    hit_s <- idx_s[pos[hit]]
    hit_t <- idx_t[hit]
  }

  obs_s <- post$post_obs_index[hit_s]
  obs_t <- post$post_obs_index[hit_t]
  oy <- context$observations$y
  y_s <- oy[obs_s]
  y_t <- oy[obs_t]
  dy <- y_t - y_s
  if (any(!is.finite(dy))) {
    cdid_abort("Pair difference overflow.", c("cdid_ipw_outcome_error", "cdid_pairs_outcome_error"))
  }

  list(
    a_index = pa[hit_s],
    unit_index = post$post_unit_index[hit_s],
    obs_s = obs_s,
    obs_t = obs_t,
    y_s = y_s,
    y_t = y_t,
    dy = dy
  )
}

query_observed_pair <- function(context, s, t, a_bounds = NULL, unit_bounds = NULL) {
  desc <- make_date_pair_descriptor(s, t, context)
  a_start <- 1L
  a_end <- context$n_a_units
  if (!is.null(a_bounds)) {
    if (!is.atomic(a_bounds) || !is.null(dim(a_bounds)) || is.object(a_bounds) ||
        !is.numeric(a_bounds) || length(a_bounds) != 2L || anyNA(a_bounds) ||
        any(!is.finite(a_bounds) | a_bounds != floor(a_bounds)) ||
        a_bounds[1L] < 1 || a_bounds[2L] > context$n_a_units || a_bounds[1L] > a_bounds[2L]) {
      cdid_abort(
        "`a_bounds` must be an integer pair `c(start, end)` within `1..n_a_units`.",
        c("cdid_ipw_bound_error", "cdid_ipw_input_error")
      )
    }
    a_start <- as.integer(a_bounds[1L])
    a_end <- as.integer(a_bounds[2L])
  }
  if (!is.null(unit_bounds)) {
    if (!is.atomic(unit_bounds) || !is.null(dim(unit_bounds)) || is.object(unit_bounds) ||
        !is.numeric(unit_bounds) || length(unit_bounds) != 2L || anyNA(unit_bounds) ||
        any(!is.finite(unit_bounds) | unit_bounds != floor(unit_bounds)) ||
        unit_bounds[1L] < 1 || unit_bounds[2L] > context$n_total_units ||
        unit_bounds[1L] > unit_bounds[2L]) {
      cdid_abort(
        "`unit_bounds` must be an integer pair `c(start, end)` within `1..n_total_units`.",
        c("cdid_ipw_bound_error", "cdid_ipw_input_error")
      )
    }
    u_lo <- lower_bound_int(context$a_unit_index, 1L, context$n_a_units, as.integer(unit_bounds[1L]))
    u_hi <- upper_bound_int(context$a_unit_index, 1L, context$n_a_units, as.integer(unit_bounds[2L]))
    a_start <- max(a_start, u_lo)
    a_end <- min(a_end, u_hi)
  }

  res <- query_observed_pair_internal(
    context, desc$s, desc$t,
    a_start = a_start, a_end = a_end,
    ds = desc$s_calendar_index, dt = desc$t_calendar_index
  )
  n_hit <- length(res$a_index)
  data.frame(
    a_index = res$a_index,
    unit_index = res$unit_index,
    id = context$a_id[res$a_index],
    g = context$a_g[res$a_index],
    s = rep(desc$s, n_hit),
    t = rep(desc$t, n_hit),
    obs_s = res$obs_s,
    obs_t = res$obs_t,
    y_s = res$y_s,
    y_t = res$y_t,
    dy = res$dy,
    weight = context$w_a[res$a_index],
    stringsAsFactors = FALSE
  )
}

extract_observed_pair <- query_observed_pair

validate_ipw_cohorts <- function(g, K) {
  if (!ipw_scalar_numeric(g, finite = TRUE)) {
    cdid_abort(
      "Target cohort `g` must be a single finite numeric value.",
      c("cdid_ipw_cohort_error", "cdid_ipw_input_error")
    )
  }
  if (missing(K) || is.null(K) || !is.atomic(K) || !is.null(dim(K)) ||
      is.object(K) || !is.numeric(K) || length(K) < 1L ||
      anyNA(K) || any(is.nan(K) | K == -Inf)) {
    cdid_abort(
      "Declared controls `K` must be an explicit nonempty numeric cohort vector (finite values or `+Inf`); control labels are not inferred.",
      c("cdid_ipw_cohort_error", "cdid_ipw_input_error")
    )
  }
  g_num <- as.double(g)
  K_num <- as.double(K)
  if (anyDuplicated(K_num) || any(K_num == g_num)) {
    cdid_abort(
      "Declared controls `K` must be distinct cohort values excluding target `g`.",
      c("cdid_ipw_cohort_error", "cdid_ipw_input_error")
    )
  }
  list(g = g_num, K = K_num, cohorts = c(g_num, K_num))
}

fit_ipw_link <- function(context, g, K, s, t, control,
                         link_id = NULL, functional_key = NULL) {
  validate_ipw_context(context)
  ctrl <- validate_ipw_glm_control(control)
  coh <- validate_ipw_cohorts(g, K)
  desc <- make_date_pair_descriptor(s, t, context)

  if (!is.null(link_id) &&
      (!is.character(link_id) || !is.null(dim(link_id)) || is.object(link_id) ||
       length(link_id) != 1L || is.na(link_id) || !nzchar(link_id))) {
    cdid_abort("`link_id` must be a nonempty scalar character string.", c("cdid_ipw_input_error"))
  }
  if (!is.null(functional_key) &&
      (!is.character(functional_key) || !is.null(dim(functional_key)) || is.object(functional_key) ||
       length(functional_key) != 1L || is.na(functional_key) || !nzchar(functional_key))) {
    cdid_abort("`functional_key` must be a nonempty scalar character string.", c("cdid_ipw_input_error"))
  }

  cohorts <- coh$cohorts
  m_cohorts <- length(cohorts)
  ci <- match(cohorts, context$a_cohorts)
  if (anyNA(ci) || any(context$cohort_a_count[ci] < 1L) ||
      any(!is.finite(context$cohort_a_weight[ci]) | context$cohort_a_weight[ci] <= 0)) {
    cdid_abort(
      "All declared cohorts `{g, K}` require positive A-unit count and weight mass; unsupported cohorts are not silently discarded.",
      c("cdid_ipw_cohort_support_error", "cdid_ipw_support_error")
    )
  }

  n <- as.double(context$n_a_units)
  group_a_count <- context$cohort_a_count[ci]
  group_mass <- context$cohort_a_weight[ci]
  q <- group_mass / n
  if (any(!is.finite(q) | q <= 0)) {
    cdid_abort(
      "Invalid A-unit cohort share normalization `q_j`.",
      c("cdid_ipw_normalization_error", "cdid_ipw_fit_error")
    )
  }

  pq <- query_observed_pair_internal(
    context, desc$s, desc$t,
    ds = desc$s_calendar_index, dt = desc$t_calendar_index
  )
  n_obs_pair_all <- length(pq$a_index)
  if (n_obs_pair_all < 1L) {
    cdid_abort(
      "Requested date pair has zero empirical A-unit observations; requested operator is unavailable.",
      c("cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }

  obs_g <- context$a_g[pq$a_index]
  obs_gi <- match(obs_g, cohorts)
  in_risk <- which(!is.na(obs_gi))
  if (length(in_risk) < 1L) {
    cdid_abort(
      "Requested date pair has zero empirical risk-set observations; requested operator is unavailable.",
      c("cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }

  risk_a_idx <- pq$a_index[in_risk]
  risk_gi <- obs_gi[in_risk]
  risk_dy <- pq$dy[in_risk]
  risk_w <- context$w_a[risk_a_idx]

  pair_count <- tabulate(risk_gi, nbins = m_cohorts)
  if (any(pair_count < 1L)) {
    cdid_abort(
      "At least one declared cohort in `{g, K}` has zero empirical pair count on `(s, t)`; zero-support cohorts are not silently discarded.",
      c("cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }

  pair_mass <- vapply(
    seq_len(m_cohorts),
    function(j) sum(risk_w[risk_gi == j]),
    numeric(1L)
  )
  if (length(pair_mass) != m_cohorts || any(!is.finite(pair_mass) | pair_mass <= 0)) {
    cdid_abort(
      "Every declared cohort requires finite positive observed-pair weight mass.",
      c("cdid_ipw_pair_support_error", "cdid_ipw_support_error")
    )
  }

  rates <- pair_mass / group_mass
  if (any(!is.finite(rates) | rates <= 0 | rates > 1)) {
    cdid_abort(
      "Invalid marginal pair observation rates `r_j`.",
      c("cdid_ipw_rate_error", "cdid_ipw_normalization_error", "cdid_ipw_fit_error")
    )
  }

  Z_use <- context$Z_a[risk_a_idx, , drop = FALSE]
  p <- ncol(Z_use)
  if (nrow(Z_use) < p) {
    cdid_abort(
      "Observed risk-set pair count is smaller than covariate dimension `ncol(Z)`.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }
  D_use <- as.double(risk_gi == 1L)
  r_use <- rates[risk_gi]
  fit_weights <- risk_w / r_use
  if (any(!is.finite(fit_weights) | fit_weights <= 0)) {
    cdid_abort(
      "Nonfinite or nonpositive inverse-rate propensity weights `w_i / r_G`.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }

  # `binomial()$initialize` warns on non-integer prior weights (`w_i / r_G`).
  # That exact known warning is muffled narrowly; convergence/boundary warnings
  # are captured and turned into typed errors below without retrying or trimming.
  glm_warnings <- character(0L)
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm.fit(
        x = Z_use,
        y = D_use,
        weights = fit_weights,
        family = stats::binomial(link = "logit"),
        control = ctrl,
        intercept = FALSE
      ),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("non-integer #successes", msg, fixed = TRUE)) {
          invokeRestart("muffleWarning")
        } else if (grepl("did not converge", msg, fixed = TRUE) ||
                   grepl("fitted probabilities numerically 0 or 1 occurred", msg, fixed = TRUE)) {
          glm_warnings <<- c(glm_warnings, msg)
          invokeRestart("muffleWarning")
        } else {
          cdid_abort(
            paste0("Propensity `glm.fit` warning: ", msg),
            c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
          )
        }
      }
    ),
    cdid_error = function(e) stop(e),
    error = function(e) {
      cdid_abort(
        paste0("Propensity `glm.fit` failed: ", conditionMessage(e)),
        c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
      )
    }
  )

  has_boundary_warn <- any(grepl("fitted probabilities numerically 0 or 1 occurred", glm_warnings, fixed = TRUE))
  has_conv_warn <- any(grepl("did not converge", glm_warnings, fixed = TRUE))

  if (is.null(fit$rank) || fit$rank < p || is.null(fit$qr$qr) || anyNA(fit$coefficients)) {
    cdid_abort(
      "Weighted propensity logit has deficient model rank.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }
  qr_R <- fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  qr_diag <- abs(diag(qr_R))
  if (any(!is.finite(qr_diag) | qr_diag <= 0)) {
    cdid_abort(
      "Weighted propensity logit QR factor is singular or nonfinite.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }

  if (isTRUE(fit$boundary) || has_boundary_warn) {
    cdid_abort(
      "Weighted propensity logit encountered boundary probabilities (0 or 1).",
      c("cdid_ipw_boundary_error", "cdid_ipw_fit_error")
    )
  }
  if (!isTRUE(fit$converged) || has_conv_warn) {
    cdid_abort(
      "Weighted propensity logit did not converge under caller-supplied `glm.control`.",
      c("cdid_ipw_convergence_error", "cdid_ipw_fit_error")
    )
  }

  beta <- as.double(fit$coefficients)
  names(beta) <- colnames(Z_use)
  if (length(beta) != p || any(!is.finite(beta))) {
    cdid_abort(
      "Weighted propensity logit produced nonfinite coefficients.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }

  eta_use <- drop(Z_use %*% beta)
  odds_use <- exp(eta_use)
  prob_use <- odds_use / (1 + odds_use)
  if (any(!is.finite(eta_use)) || any(!is.finite(odds_use)) ||
      any(odds_use <= 0) || any(!is.finite(prob_use))) {
    cdid_abort(
      "Weighted propensity logit produced nonfinite linear predictors, odds, or probabilities.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  prob_eps <- 10 * .Machine$double.eps
  if (any(prob_use <= prob_eps | prob_use >= 1 - prob_eps | prob_use * (1 - prob_use) <= 0)) {
    cdid_abort(
      "Weighted propensity logit fitted probabilities reached numerical boundary (0 or 1).",
      c("cdid_ipw_boundary_error", "cdid_ipw_fit_error")
    )
  }

  is_t <- risk_gi == 1L
  is_c <- !is_t
  a_obs <- fit_weights[is_t]
  b_obs <- fit_weights[is_c] * odds_use[is_c]
  dy_t <- risk_dy[is_t]
  dy_c <- risk_dy[is_c]

  sum_a <- sum(a_obs)
  mean_a <- sum_a / n
  sum_b <- sum(b_obs)
  mean_b <- sum_b / n
  if (!is.finite(sum_a) || sum_a <= 0 || !is.finite(mean_a) || mean_a <= 0 ||
      !is.finite(sum_b) || sum_b <= 0 || !is.finite(mean_b) || mean_b <= 0) {
    cdid_abort(
      "Invalid treated or control outcome normalization denominator.",
      c("cdid_ipw_normalization_error", "cdid_ipw_fit_error")
    )
  }

  Tmean <- sum(a_obs * dy_t) / sum_a
  Cmean <- sum(b_obs * dy_c) / sum_b
  mu <- Tmean - Cmean
  if (!is.finite(Tmean) || !is.finite(Cmean) || !is.finite(mu)) {
    cdid_abort(
      "Nonfinite treated mean, control mean, or link estimate.",
      c("cdid_ipw_normalization_error", "cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }

  w_hess <- fit_weights * prob_use * (1 - prob_use)
  H <- crossprod(Z_use, Z_use * as.numeric(w_hess)) / n
  if (any(!is.finite(H))) {
    cdid_abort(
      "Nonfinite weighted propensity Hessian `H`.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }
  H_chol <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(H_chol) || any(!is.finite(H_chol)) || any(diag(H_chol) <= 0)) {
    cdid_abort(
      "Weighted propensity Hessian `H` is not positive definite.",
      c("cdid_ipw_rank_error", "cdid_ipw_fit_error")
    )
  }

  # Score rate Jacobian columns A_{beta, r_j} = -mean(w_i 1_j S / r_j^2 Z (D - p))
  # for ALL declared cohorts in {g, K}, including the target cohort column.
  cohort_labels <- c("g", paste0("k_", seq_along(coh$K)))
  A_beta_r <- vapply(seq_len(m_cohorts), function(j) {
    idx_j <- which(risk_gi == j)
    zj <- Z_use[idx_j, , drop = FALSE]
    wt_j <- risk_w[idx_j] * (D_use[idx_j] - prob_use[idx_j]) / (rates[j]^2)
    -colSums(zj * as.numeric(wt_j)) / n
  }, numeric(p))
  if (p == 1L) {
    A_beta_r <- matrix(A_beta_r, nrow = 1L, ncol = m_cohorts)
  }
  dimnames(A_beta_r) <- list(colnames(Z_use), cohort_labels)

  # Outcome partials: G_beta = -sum(b Z (delta - Cmean)) / sum(b);
  # G_{r_g} = 0 directly; G_{r_j} = sum(w_i 1_j S odds (delta - Cmean)) / (r_j^2 sum(b)) for j in K.
  Z_c <- Z_use[is_c, , drop = FALSE]
  res_c <- dy_c - Cmean
  G_beta <- -colSums(Z_c * as.numeric(b_obs * res_c)) / sum_b
  names(G_beta) <- colnames(Z_use)

  G_r <- numeric(m_cohorts)
  G_r[1L] <- 0.0
  for (j in seq_len(m_cohorts)[-1L]) {
    idx_j <- which(risk_gi == j)
    G_r[j] <- sum(risk_w[idx_j] * odds_use[idx_j] * (risk_dy[idx_j] - Cmean)) /
      (rates[j]^2 * sum_b)
  }
  names(G_r) <- cohort_labels

  if (any(!is.finite(A_beta_r)) || any(!is.finite(G_beta)) || any(!is.finite(G_r))) {
    cdid_abort(
      "Nonfinite propensity or outcome partial derivatives.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }

  # Factor H once: loading_beta = H^{-1} G_beta via triangular solves on H_chol
  # (never forming H^{-1}), and loading_rate = G_r + A_{beta, r}' loading_beta.
  loading_beta <- drop(backsolve(H_chol, forwardsolve(t(H_chol), G_beta)))
  names(loading_beta) <- colnames(Z_use)
  loading_rate <- as.double(G_r + drop(crossprod(A_beta_r, loading_beta)))
  names(loading_rate) <- cohort_labels

  rate_coef_s0 <- -(rates / q) * loading_rate
  rate_coef_s1 <- ((1 - rates) / q) * loading_rate
  if (any(!is.finite(loading_beta)) || any(!is.finite(loading_rate)) ||
      any(!is.finite(rate_coef_s0)) || any(!is.finite(rate_coef_s1))) {
    cdid_abort(
      "Nonfinite replay loadings.",
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error")
    )
  }

  names(group_a_count) <- names(pair_count) <- cohort_labels
  names(group_mass) <- names(pair_mass) <- names(q) <- names(rates) <- cohort_labels

  structure(
    list(
      estimate = mu,
      treated_mean = Tmean,
      control_mean = Cmean,
      pair = list(
        s = desc$s,
        t = desc$t,
        s_calendar_index = desc$s_calendar_index,
        t_calendar_index = desc$t_calendar_index,
        link_id = link_id,
        functional_key = functional_key
      ),
      cohorts = list(
        g = coh$g,
        K = coh$K,
        cohorts = cohorts,
        cohort_labels = cohort_labels,
        n_cohorts = m_cohorts
      ),
      counts = list(
        n_a_units = context$n_a_units,
        n_risk_a_units = sum(group_a_count),
        n_outside_risk_a_units = context$n_a_units - sum(group_a_count),
        n_observed_pair_units = n_obs_pair_all,
        n_observed_risk_units = length(risk_a_idx),
        n_nonobserved_risk_units = sum(group_a_count) - length(risk_a_idx),
        n_treated_obs = sum(is_t),
        n_control_obs = sum(is_c),
        by_cohort_a = group_a_count,
        by_cohort_pair = pair_count
      ),
      masses = list(
        a_weight_total = context$sum_w_a,
        cohort_a_weight = group_mass,
        cohort_pair_weight = pair_mass,
        q = q,
        rates = rates
      ),
      means = list(
        treated_mean = Tmean,
        control_mean = Cmean,
        estimate = mu,
        sum_a = sum_a,
        mean_a = mean_a,
        sum_b = sum_b,
        mean_b = mean_b
      ),
      coefficients = beta,
      factors = list(
        H = H,
        H_chol = H_chol,
        qr_R = qr_R,
        qr_pivot = as.integer(fit$qr$pivot),
        qr_tol = as.double(fit$qr$tol),
        rank = as.integer(fit$rank)
      ),
      partials = list(
        G_beta = G_beta,
        G_r = G_r,
        A_beta_r = A_beta_r
      ),
      loadings = list(
        loading_beta = loading_beta,
        loading_rate = loading_rate,
        rate_coef_s0 = rate_coef_s0,
        rate_coef_s1 = rate_coef_s1
      ),
      diagnostics = list(
        converged = isTRUE(fit$converged),
        boundary = isTRUE(fit$boundary),
        iter = as.integer(fit$iter),
        rank = as.integer(fit$rank),
        n_covariates = p,
        df_residual = as.double(fit$df.residual),
        deviance = as.double(fit$deviance),
        null_deviance = as.double(fit$null.deviance),
        aic = as.double(fit$aic),
        control = ctrl
      ),
      context_ref = context$ref
    ),
    class = "cdid_ipw_link_fit"
  )
}

estimate_ipw_link <- fit_ipw_link

normalize_ipw_fits <- function(context, fits) {
  validate_ipw_context(context)
  single <- inherits(fits, "cdid_ipw_link_fit")
  fit_list <- if (single) list(fits) else fits
  if (!is.list(fit_list) || length(fit_list) < 1L ||
      any(!vapply(fit_list, inherits, logical(1L), what = "cdid_ipw_link_fit"))) {
    cdid_abort(
      "`fits` must be a `cdid_ipw_link_fit` or a nonempty list of `cdid_ipw_link_fit` objects.",
      c("cdid_ipw_input_error")
    )
  }
  for (ell in seq_along(fit_list)) {
    if (!identical(fit_list[[ell]]$context_ref, context$ref)) {
      cdid_abort(
        "All fitted link caches must share the exact same `cdid_ipw_context` instance.",
        c("cdid_ipw_context_mismatch_error", "cdid_ipw_context_error")
      )
    }
  }
  link_ids <- vapply(seq_along(fit_list), function(ell) {
    id_ell <- fit_list[[ell]]$pair$link_id
    if (!is.null(id_ell) && nzchar(id_ell)) id_ell else paste0("link_", ell)
  }, character(1L))
  if (!is.null(names(fit_list)) && all(nzchar(names(fit_list)))) {
    link_ids <- names(fit_list)
  }
  list(fits = fit_list, link_ids = link_ids, single = single)
}

preflight_ipw_replay <- function(context, fits, chunk_size, max_entries = NULL) {
  norm <- normalize_ipw_fits(context, fits)
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
        "IF materialization exceeds explicit entry bound.",
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
    covariance_matrix_entries = if (n_links <= limit / n_links) n_links^2 else Inf,
    link_ids = norm$link_ids
  )
}

replay_ipw_if <- function(context, fits, chunk_size, callback,
                          include_components = FALSE) {
  norm <- normalize_ipw_fits(context, fits)
  pre <- preflight_ipw_replay(context, norm$fits, chunk_size)
  if (!is.function(callback)) {
    cdid_abort("`callback` must be a function.", c("cdid_ipw_input_error", "cdid_pairs_data_error"))
  }
  validate_ipw_components(include_components)

  fit_list <- norm$fits
  link_ids <- norm$link_ids
  n_links <- length(fit_list)
  n_a <- context$n_a_units
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
      comp_rate <- matrix(0.0, nrow = b_len, ncol = n_links, dimnames = list(NULL, link_ids))
      comp_score <- matrix(0.0, nrow = b_len, ncol = n_links, dimnames = list(NULL, link_ids))
      comp_treated <- matrix(0.0, nrow = b_len, ncol = n_links, dimnames = list(NULL, link_ids))
      comp_control <- matrix(0.0, nrow = b_len, ncol = n_links, dimnames = list(NULL, link_ids))
    }

    for (ell in seq_len(n_links)) {
      fit <- fit_list[[ell]]
      gi_chunk <- match(u_g, fit$cohorts$cohorts)
      in_risk_chunk <- which(!is.na(gi_chunk))
      if (length(in_risk_chunk) == 0L) next

      gj_all <- gi_chunk[in_risk_chunk]
      col_rate <- numeric(b_len)
      col_rate[in_risk_chunk] <- u_w[in_risk_chunk] * fit$loadings$rate_coef_s0[gj_all]
      col_score <- numeric(b_len)
      col_treated <- numeric(b_len)
      col_control <- numeric(b_len)

      pq <- query_observed_pair_internal(
        context,
        fit$pair$s,
        fit$pair$t,
        a_start = a_start,
        a_end = a_end,
        ds = fit$pair$s_calendar_index,
        dt = fit$pair$t_calendar_index
      )
      if (length(pq$a_index) > 0L) {
        local_pos <- pq$a_index - a_start + 1L
        obs_gj <- gi_chunk[local_pos]
        obs_keep <- which(!is.na(obs_gj))
        if (length(obs_keep) > 0L) {
          pos_k <- local_pos[obs_keep]
          a_k <- pq$a_index[obs_keep]
          gj_k <- obs_gj[obs_keep]
          dy_k <- pq$dy[obs_keep]
          w_k <- u_w[pos_k]
          r_k <- fit$masses$rates[gj_k]
          D_k <- as.double(gj_k == 1L)

          Z_k <- context$Z_a[a_k, , drop = FALSE]
          eta_k <- drop(Z_k %*% fit$coefficients)
          odds_k <- exp(eta_k)
          prob_k <- odds_k / (1 + odds_k)

          col_rate[pos_k] <- w_k * fit$loadings$rate_coef_s1[gj_k]
          z_load_k <- drop(Z_k %*% fit$loadings$loading_beta)
          col_score[pos_k] <- (w_k / r_k) * (D_k - prob_k) * z_load_k

          is_tk <- gj_k == 1L
          if (any(is_tk)) {
            pos_t <- pos_k[is_tk]
            col_treated[pos_t] <- (w_k[is_tk] / r_k[is_tk]) / fit$means$mean_a *
              (dy_k[is_tk] - fit$means$treated_mean)
          }
          is_ck <- !is_tk
          if (any(is_ck)) {
            pos_c <- pos_k[is_ck]
            col_control[pos_c] <- -(w_k[is_ck] * odds_k[is_ck] / r_k[is_ck]) /
              fit$means$mean_b * (dy_k[is_ck] - fit$means$control_mean)
          }
        }
      }

      col_phi <- col_treated + col_control + col_score + col_rate
      if (any(!is.finite(col_phi))) {
        cdid_abort("Nonfinite IF contribution encountered during replay.", c("cdid_ipw_numeric_error"))
      }
      phi_mat[, ell] <- col_phi
      if (include_components) {
        comp_rate[, ell] <- col_rate
        comp_score[, ell] <- col_score
        comp_treated[, ell] <- col_treated
        comp_control[, ell] <- col_control
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
            rate = comp_rate,
            beta_score = comp_score,
            treated = comp_treated,
            control = comp_control
          )
        } else NULL
      ),
      class = "cdid_ipw_if_batch"
    )
    callback(batch) # Deliberately discard callback result; propagate callback errors.
    emitted_units <- emitted_units + b_len
    max_batch_units <- max(max_batch_units, b_len)
    a_start <- as.double(a_end) + 1
  }

  if (emitted_units != n_a) {
    cdid_abort("IF replay unit count mismatch.", c("cdid_ipw_count_error"))
  }
  invisible(list(
    n_a_units = emitted_units,
    n_links = n_links,
    n_batches = batches,
    max_batch_units = max_batch_units,
    max_batch_entries = as.double(max_batch_units) * as.double(n_links)
  ))
}

traverse_ipw_if <- replay_ipw_if

materialize_ipw_if <- function(context, fits, max_entries, chunk_size,
                               include_components = FALSE) {
  if (missing(max_entries)) {
    cdid_abort("Supply an explicit `max_entries` bound for IF materialization.",
               c("cdid_ipw_bound_error", "cdid_pairs_bound_error"))
  }
  bound <- ipw_bound(max_entries, "max_entries", allow_zero = TRUE)
  validate_ipw_components(include_components)
  norm <- normalize_ipw_fits(context, fits)
  pre <- preflight_ipw_replay(context, norm$fits, chunk_size, max_entries = bound)
  # The existing entry bound applies to EACH N * L matrix. Components add
  # four such matrices; phi/phi_matrix share one matrix until a caller changes it.
  pre$materialized_matrix_count <- if (include_components) 5L else 1L
  pre$materialized_matrix_entries <- pre$full_matrix_entries * pre$materialized_matrix_count
  n_a <- pre$n_a_units
  n_links <- pre$n_links
  link_ids <- pre$link_ids

  phi <- matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids))
  comp <- if (include_components) {
    list(
      rate = matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids)),
      beta_score = matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids)),
      treated = matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids)),
      control = matrix(NA_real_, nrow = n_a, ncol = n_links, dimnames = list(NULL, link_ids))
    )
  } else NULL

  replay_info <- replay_ipw_if(
    context = context,
    fits = norm$fits,
    chunk_size = chunk_size,
    include_components = include_components,
    callback = function(batch) {
      idx <- batch$a_index
      phi[idx, ] <<- batch$phi
      if (include_components) {
        comp$rate[idx, ] <<- batch$components$rate
        comp$beta_score[idx, ] <<- batch$components$beta_score
        comp$treated[idx, ] <<- batch$components$treated
        comp$control[idx, ] <<- batch$components$control
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

stream_ipw_covariance <- function(context, fits, chunk_size) {
  norm <- normalize_ipw_fits(context, fits)
  pre <- preflight_ipw_replay(context, norm$fits, chunk_size)
  n_a <- as.double(pre$n_a_units)
  n_links <- pre$n_links
  link_ids <- pre$link_ids

  if (!is.finite(pre$covariance_matrix_entries) ||
      pre$covariance_matrix_entries > .Machine$integer.max) {
    cdid_abort("Covariance matrix exceeds supported integer/double dimensions.",
               c("cdid_ipw_bound_error", "cdid_pairs_bound_error"))
  }
  sum_phi <- numeric(n_links)
  replay_ipw_if(context, norm$fits, chunk_size, function(batch) {
    sum_phi <<- sum_phi + colSums(batch$phi)
    if (any(!is.finite(sum_phi))) {
      cdid_abort("Nonfinite IF sum during covariance replay.", "cdid_ipw_numeric_error")
    }
  })
  mean_phi <- sum_phi / n_a
  names(mean_phi) <- link_ids

  cov_mat <- matrix(0.0, nrow = n_links, ncol = n_links, dimnames = list(link_ids, link_ids))
  replay_info <- replay_ipw_if(context, norm$fits, chunk_size, function(batch) {
    centered <- sweep(batch$phi, 2L, mean_phi, FUN = "-")
    if (any(!is.finite(centered))) {
      cdid_abort("Nonfinite centered IF during covariance replay.", "cdid_ipw_numeric_error")
    }
    cov_mat <<- cov_mat + crossprod(centered)
    if (any(!is.finite(cov_mat))) {
      cdid_abort("Nonfinite IF crossproduct or sum during covariance replay.", "cdid_ipw_numeric_error")
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
