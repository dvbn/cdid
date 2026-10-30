# Private CP02 immutable common panel/frame storage and bounded pair/nuisance readers.
#
# Scope and boundaries:
#   * Normalizes the full finite numeric calendar, exact unit IDs, invariant cohorts,
#     and long observation positions once in a sealed common panel owner
#     (`cdid_ipw_panel_owner`). Retains singleton and unobserved unit metadata and
#     explicit missing-outcome (`y = NA` row present) versus missing-row indexing.
#   * Validates raw values and observation metadata during ingestion without globally
#     subtracting every outcome pair (`construct_observed_pairs` remains unchanged).
#   * Separates long panel storage from explicit frame/design owners
#     (`cdid_ipw_frame_owner`) so risk/cohort changes can share long panel data,
#     while changed actual A-unit set, unit order, weights, or design Z require
#     distinct verified frame owners (labels alone never prove compatibility).
#   * Observing at least twice (`n_obs >= 2`) defines the canonical A-unit set;
#     it does not redefine population ATT.
#   * Compact nuisance views (`cdid_ipw_nuisance_view`) and link views
#     (`cdid_ipw_link_view`) retain only sealed owner/frame references, `(s, t)`
#     endpoint dates, ordered cohort/operator descriptors, and optional small
#     sufficient statistics/candidates --- never per-link endpoint copies, full-N
#     position maps, or full-N S/risk/cohort/observer masks.
#   * Queries bounded contiguous A-unit ranges `[a_start, a_end]`; exposes all
#     queried row identities, risk membership, and observation metadata (including
#     risk-set nonobservers `R \ O`); filters risk and pair observation positions
#     BEFORE reading any `y` values; and returns only raw stored endpoints for
#     eligible rows (`O = R intersect S`), with no cached rounded differences.
#   * Outcome-free queries read zero endpoints. No estimator, covariance, fit,
#     Monte Carlo, or public inference policy is implemented here.

cdid_ipw_owner_abort <- function(message, class = character(), diagnostics = list()) {
  stop(structure(
    class = unique(c(
      class,
      "cdid_ipw_owner_refusal",
      "cdid_ipw_owner_error",
      "cdid_ipw_enclosure_error",
      "cdid_ipw_error",
      "cdid_error",
      "error",
      "condition"
    )),
    list(message = message, call = NULL, diagnostics = diagnostics)
  ))
}

ipw_owner_check_interrupt <- function() {
  gmm_sexp_aliases(NULL, NULL)
  invisible(NULL)
}

# A completed dependency check is recorded only for a locked namespace, whose
# bindings cannot be removed after loading; sourced or cloned function
# environments (development loads, test clones) re-check on every call.
cdid_dependency_checks <- new.env(parent = emptyenv())

cdid_dependency_check_done <- function(key, ns_env) {
  identical(get0(key, envir = cdid_dependency_checks, inherits = FALSE), ns_env)
}

cdid_dependency_check_record <- function(key, ns_env) {
  if (isNamespace(ns_env) && environmentIsLocked(ns_env)) assign(key, ns_env, envir = cdid_dependency_checks)
  invisible(TRUE)
}

ipw_owner_require_deps <- function() {
  ns_env <- environment(ipw_owner_require_deps)
  if (cdid_dependency_check_done("ipw_owner", ns_env)) return(invisible(TRUE))
  req_fns <- c(
    "cdid_abort", "pair_scalar_vector", "pair_count_safe", "gmm_sexp_aliases",
    "lower_bound_int", "upper_bound_int", "query_pair_mask_internal",
    "gmm_iv", "gmm_iv_sub", "gmm_iv_add", "gmm_iv_mul", "gmm_iv_div",
    "gmm_scalar_numeric", "gmm_scalar_logical", "gmm_scalar_character",
    "logit_cert_count_output", "validate_logit_certificate_cohorts",
    "validate_logit_certificate_metric_factor"
  )
  for (fn in req_fns) {
    if (!exists(fn, envir = ns_env, mode = "function", inherits = TRUE)) {
      cdid_ipw_owner_abort(
        paste0("Required private helper `", fn, "` is not loaded in the namespace."),
        "cdid_ipw_owner_dependency_error"
      )
    }
  }
  cdid_dependency_check_record("ipw_owner", ns_env)
}

ipw_owner_seal <- function(obj, binding) {
  ref <- new.env(parent = emptyenv())
  obj$ref <- ref
  assign(binding, obj, envir = ref)
  lockEnvironment(ref, bindings = TRUE)
  obj
}

ipw_owner_seal_ok <- function(obj, binding) {
  ref <- if (is.list(obj)) obj$ref else NULL
  is.environment(ref) && identical(parent.env(ref), emptyenv()) &&
    environmentIsLocked(ref) &&
    exists(binding, envir = ref, inherits = FALSE) &&
    bindingIsLocked(binding, ref) &&
    gmm_sexp_aliases(obj, get(binding, envir = ref, inherits = FALSE))
}

ipw_owner_refuse_decorative_attrs <- function(x, nm) {
  ax <- attributes(x)
  if (!is.null(ax)) {
    bad_attr <- intersect(
      names(ax),
      c("certified", "certificate", "converged", "is_pd", "status", "rho", "dy", "pre_subtracted")
    )
    if (length(bad_attr) > 0L) {
      cdid_ipw_owner_abort(
        paste0("Decorative or foreign attribute (`", bad_attr[1L], "`) on `", nm, "` is refused."),
        c("cdid_ipw_owner_foreign_attribute_error", "cdid_ipw_owner_input_error")
      )
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Safe integer bounds and real matrix vs vector dimension limits
# ---------------------------------------------------------------------------

ipw_owner_int_bound <- function(x, name, min_val = 1) {
  if (!is.atomic(x) || !is.null(dim(x)) || is.object(x) || !is.numeric(x) ||
      length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < min_val || x > 2^53 - 1) {
    cdid_ipw_owner_abort(
      paste0("`", name, "` must be an explicit finite integer in [", min_val, ", 2^53 - 1]."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_owner_control_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  as.double(x)
}

ipw_owner_safe_product <- function(a, b, max_entries = NULL, label = "allocation") {
  da <- ipw_owner_int_bound(a, paste0(label, " dim_a"), min_val = 0)
  db <- ipw_owner_int_bound(b, paste0(label, " dim_b"), min_val = 0)
  limit <- 2^53 - 1
  if (da > 0 && db > limit / da) {
    cdid_ipw_owner_abort(
      paste0("Count product for ", label, " overflows exact IEEE binary64 integer range (2^53 - 1)."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(label = label, dim_a = da, dim_b = db, exact_count_limit = limit)
    )
  }
  prod_val <- da * db
  if (!is.finite(prod_val) || prod_val > limit) {
    cdid_ipw_owner_abort(
      paste0("Count product for ", label, " exceeds exact integer range."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error")
    )
  }
  if (!is.null(max_entries) && prod_val > max_entries) {
    cdid_ipw_owner_abort(
      paste0("Count for ", label, " (", prod_val, ") exceeds explicit bound (", max_entries, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(label = label, requested = prod_val, bound = max_entries)
    )
  }
  prod_val
}

ipw_owner_safe_sum <- function(vals, max_bound = NULL, label = "sum") {
  limit <- 2^53 - 1
  total <- 0
  for (v in vals) {
    dv <- as.double(v)
    if (!is.finite(dv) || dv < 0 || dv != floor(dv) || dv > limit - total) {
      cdid_ipw_owner_abort(
        paste0("Cumulative count for ", label, " overflows exact IEEE binary64 integer range."),
        c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
        diagnostics = list(label = label, total_so_far = total, term = dv)
      )
    }
    total <- total + dv
  }
  if (!is.null(max_bound) && total > max_bound) {
    cdid_ipw_owner_abort(
      paste0("Total ", label, " (", total, ") exceeds explicit bound (", max_bound, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(label = label, requested = total, bound = max_bound)
    )
  }
  total
}

# Real 2D R matrices store (nrow, ncol) in a 32-bit integer `dim` attribute
# (each <= .Machine$integer.max). R_XLEN_T_MAX is 2^52 on this 64-bit
# R interface; exact aggregate integer counting permits at most 2^53 - 1.
ipw_owner_check_matrix_dims <- function(nrow_val, ncol_val, max_matrix_entries = NULL, label = "matrix") {
  nr <- ipw_owner_int_bound(nrow_val, paste0(label, " nrow"), min_val = 0)
  nc <- ipw_owner_int_bound(ncol_val, paste0(label, " ncol"), min_val = 0)
  int_max <- as.double(.Machine$integer.max)
  if (nr > int_max || nc > int_max) {
    cdid_ipw_owner_abort(
      paste0("Real 2D matrix dimension for ", label, " (", nr, " x ", nc,
             ") exceeds R 32-bit matrix dimension limit (`INT_MAX` = ", int_max, ")."),
      c("cdid_ipw_owner_matrix_dim_error", "cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(label = label, nrow = nr, ncol = nc, matrix_dim_limit = int_max)
    )
  }
  entries <- ipw_owner_safe_product(nr, nc, max_entries = max_matrix_entries, label = label)
  if (entries > 2^52) cdid_ipw_owner_abort("Matrix vector exceeds R_XLEN_T_MAX (2^52).", "cdid_ipw_owner_bound_error")
  entries
}

# ---------------------------------------------------------------------------
# Explicit caller resource controls (NO internal defaults)
# ---------------------------------------------------------------------------

ipw_owner_control <- function(
  max_chunk_rows,
  max_retained_entries,
  max_matrix_entries,
  max_work_ops,
  max_output_entries
) {
  ipw_owner_require_deps()
  if (missing(max_chunk_rows) || missing(max_retained_entries) ||
      missing(max_matrix_entries) || missing(max_work_ops) ||
      missing(max_output_entries)) {
    cdid_ipw_owner_abort(
      "All fields of `ipw_owner_control` (`max_chunk_rows`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) must be supplied explicitly by the caller (no internal defaults).",
      c("cdid_ipw_owner_control_error", "cdid_ipw_owner_bound_error")
    )
  }
  max_ch <- ipw_owner_int_bound(max_chunk_rows, "max_chunk_rows", min_val = 1)
  max_ret <- ipw_owner_int_bound(max_retained_entries, "max_retained_entries", min_val = 1)
  max_me <- ipw_owner_int_bound(max_matrix_entries, "max_matrix_entries", min_val = 1)
  max_wo <- ipw_owner_int_bound(max_work_ops, "max_work_ops", min_val = 1)
  max_out <- ipw_owner_int_bound(max_output_entries, "max_output_entries", min_val = 1)

  structure(
    list(
      max_chunk_rows = max_ch,
      max_retained_entries = max_ret,
      max_matrix_entries = max_me,
      max_work_ops = max_wo,
      max_output_entries = max_out
    ),
    class = "cdid_ipw_owner_control"
  )
}

validate_ipw_owner_control <- function(control) {
  ipw_owner_require_deps()
  req <- c("max_chunk_rows", "max_retained_entries", "max_matrix_entries", "max_work_ops", "max_output_entries")
  if (missing(control) || is.null(control) || !is.list(control) ||
      length(control) != length(req) || is.null(names(control)) ||
      anyDuplicated(names(control)) || !setequal(names(control), req)) {
    cdid_ipw_owner_abort(
      "Explicit `control` list (`max_chunk_rows`, `max_retained_entries`, `max_matrix_entries`, `max_work_ops`, `max_output_entries`) is required; no internal defaults exist.",
      c("cdid_ipw_owner_control_error", "cdid_ipw_owner_bound_error")
    )
  }
  ipw_owner_control(
    max_chunk_rows = control$max_chunk_rows,
    max_retained_entries = control$max_retained_entries,
    max_matrix_entries = control$max_matrix_entries,
    max_work_ops = control$max_work_ops,
    max_output_entries = control$max_output_entries
  )
}

ipw_panel_owner_control <- ipw_owner_control
ipw_frame_owner_control <- ipw_owner_control
ipw_view_control <- ipw_owner_control
ipw_reader_control <- ipw_owner_control

# Byte metadata uses CHARSXP byte lengths, not enc2utf8/charToRaw conversions.
# The existing row preflight covers this bounded count-vector pass. Shared
# string contents remain excluded, but NEW conversion/raw/integer temporaries
# and their byte-dependent work are explicitly charged before comparisons.
ipw_owner_id_byte_metadata <- function(id) {
  if (!is.character(id)) return(list(total = 0, max = 0))
  lens <- as.double(nchar(id, type = "bytes"))
  list(total = ipw_owner_safe_sum(lens, label = "character ID source byte count"),
       max = if (length(lens)) max(lens) else 0)
}

ipw_owner_preflight_id_bytes <- function(preflight, control, metadata,
                                         owner_max_bytes, lookup_count,
                                         search_width, simultaneous) {
  ctrl <- validate_ipw_owner_control(control)
  own <- ipw_owner_int_bound(owner_max_bytes, "owner_max_ID_bytes", min_val = 0)
  n <- ipw_owner_int_bound(lookup_count, "ID_lookup_count", min_val = 0)
  width <- ipw_owner_int_bound(search_width, "ID_search_width", min_val = 0)
  # At most four UTF8 bytes per source byte. For a comparison, conversion/raw/
  # integer arrays plus slices, comparison flags and which indices fit within
  # 32 slots per source byte. This also covers transient converted CHARSXP bytes
  # conservatively as 8-byte slots, rather than calling them shared contents.
  owner_bytes <- ipw_owner_safe_product(if (simultaneous) n else min(1, n), own,
                                       label = "simultaneous owner ID byte bound")
  scratch_source_bytes <- ipw_owner_safe_sum(c(metadata$total, owner_bytes), label = "ID conversion source bytes")
  extra_scratch <- ipw_owner_safe_product(32, scratch_source_bytes, ctrl$max_matrix_entries,
                                         "UTF8/raw/integer ID temporary slots")
  work_source_bytes <- ipw_owner_safe_sum(c(metadata$total,
    ipw_owner_safe_product(n, own, label = "probed owner ID bytes")), label = "ID comparison source bytes")
  extra_work <- ipw_owner_safe_product(64 * (1 + ceiling(log2(width + 1))), work_source_bytes,
                                      ctrl$max_work_ops, "byte-dependent ID normalization/comparison work")
  preflight$concurrent_entries <- ipw_owner_safe_sum(c(preflight$concurrent_entries, extra_scratch),
    ctrl$max_matrix_entries, "concurrent payload including ID conversion temporaries")
  preflight$work_ops_bound <- ipw_owner_safe_sum(c(preflight$work_ops_bound, extra_work),
    ctrl$max_work_ops, "work including byte-dependent ID comparisons")
  preflight$concurrent_vector_entries <- preflight$concurrent_vector_entries + extra_scratch
  preflight$concurrent_bytes_8double <- 8 * preflight$concurrent_entries
  if (!is.null(preflight$constructor_scratch_entries)) {
    preflight$constructor_scratch_entries <- preflight$constructor_scratch_entries + extra_scratch
  }
  preflight$id_byte_preflight <- list(source_byte_count = metadata$total,
    max_source_ID_bytes = metadata$max, max_owner_ID_bytes = own,
    lookup_count = n, search_width = width, UTF8_expansion_bound = 4,
    additional_scratch_entries = extra_scratch, additional_work_bound = extra_work)
  preflight
}

# ---------------------------------------------------------------------------
# Binary search helpers (no hash table allocations or string coercions)
# ---------------------------------------------------------------------------

# Search only requested scalars. Never coerce or check the entire sealed vector.
# Optional bounds avoid allocating a complete unit date slice for cell queries.
# CP08 balanced performance: on a strictly increasing NA-free numeric table
# searched over its whole range, the scalar search returns the unique equal
# position or NA, which is exactly match(). That vectorized lookup is taken only
# when the table is no longer than the query, so its O(n) order check and hash
# never exceed the O(m log n) search; sub-range searches, NA or NaN keys,
# classed objects and every other table keep the scalar search below.
ipw_owner_sorted_lookup_vectorizable <- function(x, sorted_vec, left, right) {
  n <- length(sorted_vec)
  n > 0L && n <= length(x) && n <= .Machine$integer.max &&
    is.numeric(left) && length(left) == 1L && !is.na(left) && left == 1 &&
    is.numeric(right) && length(right) == 1L && !is.na(right) && right == n &&
    (is.integer(x) || is.double(x)) && (is.integer(sorted_vec) || is.double(sorted_vec)) &&
    !is.object(x) && !is.object(sorted_vec) && !anyNA(x) && !anyNA(sorted_vec) &&
    !is.unsorted(sorted_vec, strictly = TRUE)
}

ipw_owner_match_sorted_double <- function(x, sorted_vec, left = 1, right = length(sorted_vec)) {
  m <- length(x)
  if (m > 0L && ipw_owner_sorted_lookup_vectorizable(x, sorted_vec, left, right)) return(match(x, sorted_vec))
  out <- if (length(sorted_vec) <= .Machine$integer.max) rep.int(NA_integer_, m) else rep.int(NA_real_, m)
  if (m == 0L || left > right) return(out)
  for (i in seq_len(m)) {
    lo <- as.double(left); hi <- as.double(right); target <- x[i]
    while (lo <= hi) {
      mid <- lo + floor((hi - lo) / 2)
      value <- sorted_vec[mid]
      if (value == target) {
        out[i] <- if (is.integer(out)) as.integer(mid) else mid
        break
      }
      if (value < target) lo <- mid + 1 else hi <- mid - 1
    }
  }
  out
}

ipw_owner_chr_lt_radix <- function(a, b) {
  if (a == b) return(FALSE)
  ra <- as.integer(charToRaw(enc2utf8(a)))
  rb <- as.integer(charToRaw(enc2utf8(b)))
  k <- min(length(ra), length(rb))
  if (k > 0L) {
    diff_pos <- which(ra[seq_len(k)] != rb[seq_len(k)])
    if (length(diff_pos) > 0L) {
      p <- diff_pos[1L]
      return(ra[p] < rb[p])
    }
  }
  length(ra) < length(rb)
}

ipw_owner_match_sorted_chr <- function(x, sorted_chr) {
  n <- length(sorted_chr)
  m <- length(x)
  out <- rep.int(NA_integer_, m)
  if (n == 0L || m == 0L) return(out)
  for (i in seq_len(m)) {
    target <- x[i]
    lo <- 1L
    hi <- n
    while (lo <= hi) {
      mid <- lo + (hi - lo) %/% 2L
      val <- sorted_chr[mid]
      if (val == target) {
        out[i] <- mid
        break
      } else if (ipw_owner_chr_lt_radix(val, target)) {
        lo <- mid + 1L
      } else {
        hi <- mid - 1L
      }
    }
  }
  out
}

ipw_owner_match_sorted_id <- function(id_key, sorted_id_key, id_type) {
  if (identical(id_type, "numeric")) {
    ipw_owner_match_sorted_double(id_key, sorted_id_key)
  } else {
    ipw_owner_match_sorted_chr(id_key, sorted_id_key)
  }
}

# Accepted sorted-posting intersection adapted for long observation positions.
# No full posting scan, coercion or hash allocation; only chunk slices are made.
ipw_owner_lower_bound <- function(vec, left, right, target) {
  lo <- as.double(left); hi <- as.double(right); ans <- hi + 1
  while (lo <= hi) {
    mid <- lo + floor((hi - lo) / 2)
    if (vec[mid] >= target) { ans <- mid; hi <- mid - 1 } else lo <- mid + 1
  }
  ans
}
ipw_owner_upper_bound <- function(vec, left, right, target) {
  lo <- as.double(left); hi <- as.double(right); ans <- lo - 1
  while (lo <= hi) {
    mid <- lo + floor((hi - lo) / 2)
    if (vec[mid] <= target) { ans <- mid; lo <- mid + 1 } else hi <- mid - 1
  }
  ans
}
ipw_owner_query_pair_mask <- function(context, s, t,
                                     a_start = 1L,
                                     a_end = context$n_a_units,
                                     ds = NULL, dt = NULL) {
  if (is.null(ds)) ds <- match(as.double(s), context$all_times)
  if (is.null(dt)) dt <- match(as.double(t), context$all_times)
  empty <- list(
    a_index = integer(0L),
    unit_index = integer(0L),
    obs_s = integer(0L),
    obs_t = integer(0L)
  )
  if (is.na(ds) || is.na(dt) || a_start > a_end) return(empty)
  post <- context$postings
  if (post$date_count[ds] == 0L || post$date_count[dt] == 0L) return(empty)

  pa <- post$post_a_index
  lo_s <- ipw_owner_lower_bound(pa, post$date_start[ds], post$date_end[ds], a_start)
  hi_s <- ipw_owner_upper_bound(pa, post$date_start[ds], post$date_end[ds], a_end)
  if (lo_s > hi_s) return(empty)

  lo_t <- ipw_owner_lower_bound(pa, post$date_start[dt], post$date_end[dt], a_start)
  hi_t <- ipw_owner_upper_bound(pa, post$date_start[dt], post$date_end[dt], a_end)
  if (lo_t > hi_t) return(empty)

  # Tighten both intervals to the shared endpoint range via binary bounds.
  min_a <- max(pa[lo_s], pa[lo_t])
  max_a <- min(pa[hi_s], pa[hi_t])
  if (min_a > max_a) return(empty)
  if (pa[lo_s] < min_a) lo_s <- ipw_owner_lower_bound(pa, lo_s, hi_s, min_a)
  if (pa[hi_s] > max_a) hi_s <- ipw_owner_upper_bound(pa, lo_s, hi_s, max_a)
  if (lo_s > hi_s) return(empty)
  if (pa[lo_t] < min_a) lo_t <- ipw_owner_lower_bound(pa, lo_t, hi_t, min_a)
  if (pa[hi_t] > max_a) hi_t <- ipw_owner_upper_bound(pa, lo_t, hi_t, max_a)
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

  list(
    a_index = pa[hit_s],
    unit_index = post$post_unit_index[hit_s],
    obs_s = post$post_obs_index[hit_s],
    obs_t = post$post_obs_index[hit_t]
  )
}

# Outcome values for an already computed pair mask; the only outcome reader.
# ---------------------------------------------------------------------------
# Preflight accounting entrypoints (executed BEFORE scans, copies, callbacks)
# ---------------------------------------------------------------------------

preflight_ipw_panel_owner <- function(n_rows, control) {
  ipw_owner_require_deps()
  ctrl <- validate_ipw_owner_control(control)
  nr <- ipw_owner_int_bound(n_rows, "nrow(data)", min_val = 1)
  if (nr > 2^52) cdid_ipw_owner_abort("Panel vector exceeds R_XLEN_T_MAX (2^52).", "cdid_ipw_owner_bound_error")

  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  B_eff <- min(nr, ctrl$max_chunk_rows)

  # Worst-case before content scan: nu <= n_rows, q <= n_rows, n_cal <= n_rows.
  # Retained vectors: unit metadata (<= 24 * nu), observations/explicit_na (<= 6 * n_rows),
  # calendar & postings (<= 8 * n_rows), scalar counts (32).
  retained_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(64, nr, max_ret, "panel unit, observation, calendar and posting vectors"),
      256
    ),
    max_bound = max_ret,
    label = "retained panel owner entries"
  )
  retained_matrix_entries <- 0
  total_retained_entries <- retained_vector_entries
  if (total_retained_entries > max_out) {
    cdid_ipw_owner_abort(
      paste0("Panel owner output entries (", total_retained_entries, ") exceed `max_output_entries` (", max_out, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(requested = total_retained_entries, bound = max_out)
    )
  }

  constructor_scratch_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(128, nr, max_me, "panel radix sort, unit match, and posting construction scratch"),
      ipw_owner_safe_product(16, B_eff, max_me, "panel chunk validation scratch"),
      256
    ),
    max_bound = max_me,
    label = "panel constructor scratch entries"
  )
  concurrent_entries <- ipw_owner_safe_sum(
    c(total_retained_entries, constructor_scratch_entries),
    max_bound = max_me,
    label = "concurrent panel owner payload entries"
  )
  work_ops_bound <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(128 * (1 + ceiling(log2(nr + 1))), nr, max_wo, "panel validation, ordering, binary lookup and posting work charges"),
      256
    ),
    max_bound = max_wo,
    label = "panel owner construction work operations"
  )

  list(
    n_input_rows = nr,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(ceiling(nr / B_eff)),
    queried_row_count = nr,
    matrix_dim_limit = as.double(.Machine$integer.max),
    vector_length_limit = 2^52,
    exact_count_limit = 2^53 - 1,
    retained_vector_entries = retained_vector_entries,
    retained_matrix_entries = retained_matrix_entries,
    total_retained_entries = total_retained_entries,
    constructor_scratch_entries = constructor_scratch_entries,
    concurrent_matrix_entries = 0,
    concurrent_vector_entries = concurrent_entries,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    work_ops_bound = work_ops_bound,
    global_pair_subtraction_performed = FALSE,
    peak_rss_promised = FALSE,
    accounting_scope = "Logical/numeric/character vector and matrix entry bounds and conservative operation charges; counts are not OS process peak RSS and exclude R SEXP headers and environment bindings."
  )
}

preflight_ipw_frame_owner <- function(
  n_a_units,
  n_total_units,
  n_supplied_id,
  p_covariates,
  n_calendar_dates,
  n_a_observed_rows,
  control
) {
  ipw_owner_require_deps()
  ctrl <- validate_ipw_owner_control(control)

  na <- ipw_owner_int_bound(n_a_units, "n_a_units", min_val = 1)
  nu <- ipw_owner_int_bound(n_total_units, "n_total_units", min_val = 1)
  nsup <- ipw_owner_int_bound(n_supplied_id, "length(id)", min_val = 1)
  p <- ipw_owner_int_bound(p_covariates, "ncol(Z)", min_val = 1)
  ncal <- ipw_owner_int_bound(n_calendar_dates, "n_calendar_dates", min_val = 1)
  qa <- ipw_owner_int_bound(n_a_observed_rows, "n_a_observed_rows", min_val = 2)

  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  B_eff <- min(nsup, ctrl$max_chunk_rows)

  # Check real 2D matrix limits on supplied Z (nsup x p) and retained Z_a (na x p)
  supplied_Z_entries <- ipw_owner_check_matrix_dims(nsup, p, max_matrix_entries = max_me, label = "supplied Z matrix")
  retained_matrix_entries <- ipw_owner_check_matrix_dims(na, p, max_matrix_entries = max_ret, label = "retained Z_a matrix")
  if (retained_matrix_entries > max_me) {
    cdid_ipw_owner_abort(
      paste0("Retained Z_a matrix entries (", retained_matrix_entries, ") exceed `max_matrix_entries` (", max_me, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error")
    )
  }

  # Frame-owned vectors (panel arrays are aliased, not copied, but accounted transparently)
  retained_vector_entries <- ipw_owner_safe_sum(
    c(
      na, # w_a
      nu, # unit_weights_all
      ipw_owner_safe_product(2, nu, max_ret, "cohort_a_weight and cohort_a_share; global cohorts <= units"),
      ipw_owner_safe_product(6, max(0, nsup - na), max_ret, "six-column external frame-only metadata"),
      256
    ),
    max_bound = max_ret,
    label = "retained frame vector entries"
  )
  total_retained_entries <- ipw_owner_safe_sum(
    c(retained_matrix_entries, retained_vector_entries),
    max_bound = max_ret,
    label = "total retained frame owner entries"
  )
  if (total_retained_entries > max_out) {
    cdid_ipw_owner_abort(
      paste0("Frame owner output entries (", total_retained_entries, ") exceed `max_output_entries` (", max_out, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(requested = total_retained_entries, bound = max_out)
    )
  }

  chunk_Bp <- ipw_owner_check_matrix_dims(B_eff, p, max_matrix_entries = max_me, label = "frame chunk B x p validation slice")
  constructor_scratch_entries <- ipw_owner_safe_sum(
    c(
      chunk_Bp,
      ipw_owner_safe_product(16, nsup + nu, max_me, "frame normalization, alignment and stable cohort grouping scratch"),
      ipw_owner_safe_product(16, B_eff, max_me, "frame chunk vector scratch"),
      256
    ),
    max_bound = max_me,
    label = "frame constructor scratch entries"
  )
  concurrent_matrix_entries <- ipw_owner_safe_sum(
    c(retained_matrix_entries, chunk_Bp),
    max_bound = max_me,
    label = "concurrent frame matrix entries"
  )
  concurrent_entries <- ipw_owner_safe_sum(
    c(total_retained_entries, constructor_scratch_entries),
    max_bound = max_me,
    label = "concurrent frame owner payload entries"
  )
  work_ops_bound <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(128 * (1 + ceiling(log2(max(nsup, nu) + 1))), supplied_Z_entries + nsup + nu, max_wo, "frame validation, alignment and stable grouping work charges"),
      ipw_owner_safe_product(8, nsup, max_wo, "ID validation, duplicate check, and weight scan work"),
      ipw_owner_safe_product(4, nu, max_wo, "unit weight and cohort aggregation work"),
      256
    ),
    max_bound = max_wo,
    label = "frame owner construction work operations"
  )

  list(
    n_a_units = na,
    n_total_units = nu,
    n_supplied_id = nsup,
    p_covariates = p,
    n_calendar_dates = ncal,
    n_a_observed_rows = qa,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(ceiling(nsup / B_eff)),
    queried_row_count = nsup,
    matrix_dim_limit = as.double(.Machine$integer.max),
    vector_length_limit = 2^52,
    exact_count_limit = 2^53 - 1,
    supplied_Z_entries = supplied_Z_entries,
    retained_matrix_entries = retained_matrix_entries,
    retained_vector_entries = retained_vector_entries,
    total_retained_entries = total_retained_entries,
    constructor_scratch_entries = constructor_scratch_entries,
    concurrent_matrix_entries = concurrent_matrix_entries,
    concurrent_vector_entries = concurrent_entries - concurrent_matrix_entries,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    work_ops_bound = work_ops_bound,
    panel_arrays_aliased = TRUE,
    peak_rss_promised = FALSE,
    accounting_scope = "Frame-owned matrix/vector payload entries and conservative validation operation charges; panel observations and postings are aliased by reference. Counts are not OS process peak RSS."
  )
}

preflight_ipw_nuisance_view <- function(
  n_a_units,
  p_covariates,
  m_cohorts,
  control,
  has_candidates = FALSE
) {
  ipw_owner_require_deps()
  ctrl <- validate_ipw_owner_control(control)

  na <- ipw_owner_int_bound(n_a_units, "n_a_units", min_val = 1)
  p <- ipw_owner_int_bound(p_covariates, "p_covariates", min_val = 1)
  m <- ipw_owner_int_bound(m_cohorts, "m_cohorts", min_val = 2)

  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  B_eff <- min(na, ctrl$max_chunk_rows)

  p2 <- if (isTRUE(has_candidates)) {
    ipw_owner_check_matrix_dims(p, p, max_matrix_entries = max_ret, label = "nuisance view F candidate matrix")
  } else 0
  if (p2 > max_me) {
    cdid_ipw_owner_abort("Nuisance view candidate matrix exceeds `max_matrix_entries`.", "cdid_ipw_owner_bound_error")
  }

  # Compact nuisance view retains O(m) cohort descriptors/support and optional O(p^2 + p) candidates;
  # NO full-N masks, NO Z copy, NO |O|-length index vector!
  retained_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(20, m, max_ret, "nuisance view declared-cohort descriptors, counts, masses, and rates"),
      if (isTRUE(has_candidates)) ipw_owner_safe_sum(c(p, 4), max_ret, "nuisance view beta0 and r candidates") else 0,
      256
    ),
    max_bound = max_ret,
    label = "retained nuisance view vector entries"
  )
  retained_matrix_entries <- p2
  total_retained_entries <- ipw_owner_safe_sum(
    c(retained_matrix_entries, retained_vector_entries),
    max_bound = max_ret,
    label = "total retained nuisance view entries"
  )
  if (total_retained_entries > max_out) {
    cdid_ipw_owner_abort(
      paste0("Nuisance view entries (", total_retained_entries, ") exceed `max_output_entries` (", max_out, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error")
    )
  }

  # Candidate normalization/finiteness checks can allocate matrix-sized scratch
  # before support traversal; charge it as well as retained candidates.
  candidate_scratch <- ipw_owner_safe_product(4, retained_matrix_entries, max_me, "nuisance candidate validation scratch")
  # Support verification scans postings in bounded chunks of size B_eff
  # After division, 6*m accumulator scratch also covers completeness's m-logical
  # comparison, up to m replacement indices and an m-rate copy; division
  # temporaries are dead and the summary vectors are charged as retained above.
  constructor_scratch_entries <- ipw_owner_safe_sum(
    c(
      candidate_scratch,
      ipw_owner_safe_product(48, B_eff, max_me, "chunked posting intersection scratch"),
      ipw_owner_safe_product(6, m, max_me, "cohort accumulator scratch"),
      256
    ),
    max_bound = max_me,
    label = "nuisance view constructor scratch entries"
  )
  concurrent_entries <- ipw_owner_safe_sum(
    c(total_retained_entries, constructor_scratch_entries),
    max_bound = max_me,
    label = "concurrent nuisance view entries"
  )
  # The cohort support/rate term (>= 96*m) includes the added <= 4*m
  # completeness comparison/index/copy/replacement work reported by the helper.
  work_ops_bound <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(32 * (1 + ceiling(log2(m + 1))), na, max_wo, "posting support and declared-cohort binary lookup work charges"),
      ipw_owner_safe_product(64 * (1 + ceiling(log2(na + 1))), ceiling(na / B_eff), max_wo, "posting binary-search work charges"),
      ipw_owner_safe_product(32 * (1 + ceiling(log2(m + 1))), m, max_wo, "cohort support and rate check ops"),
      if (isTRUE(has_candidates)) ipw_owner_safe_product(8, p2 + p + 8, max_wo, "candidate validation work charges") else 0,
      256
    ),
    max_bound = max_wo,
    label = "nuisance view construction work operations"
  )

  list(
    n_a_units = na,
    p_covariates = p,
    m_cohorts = m,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(ceiling(na / B_eff)),
    queried_row_count = na,
    matrix_dim_limit = as.double(.Machine$integer.max),
    vector_length_limit = 2^52,
    exact_count_limit = 2^53 - 1,
    retained_matrix_entries = retained_matrix_entries,
    retained_vector_entries = retained_vector_entries,
    total_retained_entries = total_retained_entries,
    constructor_scratch_entries = constructor_scratch_entries,
    concurrent_matrix_entries = retained_matrix_entries + candidate_scratch,
    concurrent_vector_entries = concurrent_entries - retained_matrix_entries - candidate_scratch,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    work_ops_bound = work_ops_bound,
    full_n_masks_retained = FALSE,
    copied_Z_retained = FALSE,
    endpoints_read = 0L,
    peak_rss_promised = FALSE,
    accounting_scope = "Compact nuisance view descriptor entries and chunked posting support verification; zero full-N masks or Z copies retained."
  )
}

preflight_ipw_link_view <- function(
  n_a_units,
  p_covariates,
  m_cohorts,
  control,
  has_candidates = FALSE
) {
  ipw_owner_require_deps()
  ctrl <- validate_ipw_owner_control(control)

  na <- ipw_owner_int_bound(n_a_units, "n_a_units", min_val = 1)
  p <- ipw_owner_int_bound(p_covariates, "p_covariates", min_val = 1)
  m <- ipw_owner_int_bound(m_cohorts, "m_cohorts", min_val = 2)

  max_ret <- ctrl$max_retained_entries
  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries
  B_eff <- min(na, ctrl$max_chunk_rows)

  pm <- if (isTRUE(has_candidates)) {
    ipw_owner_check_matrix_dims(p, m, max_matrix_entries = max_ret, label = "link view optional L/H candidate matrix")
  } else 0
  p2 <- if (isTRUE(has_candidates)) {
    ipw_owner_check_matrix_dims(p, p, max_matrix_entries = max_ret, label = "link view optional H_chol candidate matrix")
  } else 0
  retained_matrix_entries <- ipw_owner_safe_sum(
    c(pm, pm, p2), # A_beta_r and L are separately normalized and retained
    max_bound = max_ret,
    label = "retained link view candidate matrix entries"
  )
  if (retained_matrix_entries > max_me) {
    cdid_ipw_owner_abort("Link view candidate matrices exceed `max_matrix_entries`.", "cdid_ipw_owner_bound_error")
  }

  # Compact link view retains ONLY O(p + m) descriptors/candidates;
  # NO y_left_obs / y_right_obs copies, NO a_to_obs_pos N-map, NO full-N masks!
  retained_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(20, m, max_ret, "link view declared-cohort descriptors and support summary"),
      p, # loading_beta_candidate
      if (isTRUE(has_candidates)) ipw_owner_safe_sum(c(3 * p, 3 * m, 8), max_ret, "link view optional candidate vectors") else 0,
      256
    ),
    max_bound = max_ret,
    label = "retained link view vector entries"
  )
  total_retained_entries <- ipw_owner_safe_sum(
    c(retained_matrix_entries, retained_vector_entries),
    max_bound = max_ret,
    label = "total retained link view entries"
  )
  if (total_retained_entries > max_out) {
    cdid_ipw_owner_abort(
      paste0("Link view entries (", total_retained_entries, ") exceed `max_output_entries` (", max_out, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error")
    )
  }

  candidate_scratch <- ipw_owner_safe_product(4, retained_matrix_entries, max_me, "link candidate validation scratch")
  # As in nuisance preflight, 6*m covers the post-division completeness
  # comparison/index/rate-copy scratch (<= 3*m), with summaries retained above.
  constructor_scratch_entries <- ipw_owner_safe_sum(
    c(
      candidate_scratch,
      ipw_owner_safe_product(64, B_eff, max_me, "chunked posting support and streamed mask comparison scratch"),
      ipw_owner_safe_product(6, m, max_me, "cohort accumulator scratch"),
      256
    ),
    max_bound = max_me,
    label = "link view constructor scratch entries"
  )
  concurrent_entries <- ipw_owner_safe_sum(
    c(total_retained_entries, constructor_scratch_entries),
    max_bound = max_me,
    label = "concurrent link view entries"
  )
  # The cohort rate term (>= 96*m) includes completeness's <= 4*m work.
  work_ops_bound <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(64 * (1 + ceiling(log2(m + 1))), na, max_wo, "posting support and streamed mask work charges"),
      ipw_owner_safe_product(64 * (1 + ceiling(log2(na + 1))), ceiling(na / B_eff), max_wo, "posting binary-search work charges"),
      ipw_owner_safe_product(32 * (1 + ceiling(log2(m + 1))), m, max_wo, "cohort rate check ops"),
      if (isTRUE(has_candidates)) ipw_owner_safe_product(8, retained_matrix_entries + 4 * p + 4 * m + 16, max_wo, "candidate validation work charges") else 0,
      256
    ),
    max_bound = max_wo,
    label = "link view construction work operations"
  )

  list(
    n_a_units = na,
    p_covariates = p,
    m_cohorts = m,
    effective_chunk_rows = logit_cert_count_output(B_eff),
    n_chunks_per_pass = logit_cert_count_output(ceiling(na / B_eff)),
    queried_row_count = na,
    matrix_dim_limit = as.double(.Machine$integer.max),
    vector_length_limit = 2^52,
    exact_count_limit = 2^53 - 1,
    retained_endpoint_entries = 0,
    retained_position_map_entries = 0,
    retained_matrix_entries = retained_matrix_entries,
    retained_vector_entries = retained_vector_entries,
    total_retained_entries = total_retained_entries,
    constructor_scratch_entries = constructor_scratch_entries,
    concurrent_matrix_entries = retained_matrix_entries + candidate_scratch,
    concurrent_vector_entries = concurrent_entries - retained_matrix_entries - candidate_scratch,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    work_ops_bound = work_ops_bound,
    per_link_endpoint_copies_retained = FALSE,
    full_n_position_map_retained = FALSE,
    full_n_masks_retained = FALSE,
    endpoints_read = 0L,
    peak_rss_promised = FALSE,
    accounting_scope = "Compact link view descriptor entries and chunked posting verification; zero endpoint copies, zero N-position maps, zero full-N masks."
  )
}

preflight_ipw_range_query <- function(
  b_len,
  p_covariates,
  m_cohorts,
  control,
  read_outcomes = FALSE,
  include_design = FALSE,
  enclose_differences = FALSE
) {
  ipw_owner_require_deps()
  ctrl <- validate_ipw_owner_control(control)

  B_d <- ipw_owner_int_bound(b_len, "chunk_rows", min_val = 1)
  p <- ipw_owner_int_bound(p_covariates, "p_covariates", min_val = 1)
  m <- ipw_owner_int_bound(m_cohorts, "m_cohorts", min_val = 2)

  if (B_d > ctrl$max_chunk_rows) {
    cdid_ipw_owner_abort(
      paste0("Requested contiguous A-range length (", B_d, ") exceeds `control$max_chunk_rows` (", ctrl$max_chunk_rows, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_enclosure_bound_error"),
      diagnostics = list(requested_chunk_rows = B_d, max_chunk_rows = ctrl$max_chunk_rows)
    )
  }

  max_me <- ctrl$max_matrix_entries
  max_wo <- ctrl$max_work_ops
  max_out <- ctrl$max_output_entries

  Bp <- if (isTRUE(include_design)) {
    ipw_owner_check_matrix_dims(B_d, p, max_matrix_entries = max_me, label = "query chunk B x p design matrix")
  } else 0
  diff_iv_mat <- if (isTRUE(read_outcomes) && isTRUE(enclose_differences)) {
    ipw_owner_check_matrix_dims(B_d, 2L, max_matrix_entries = max_me, label = "query chunk B x 2 outward difference interval matrix")
  } else 0

  # Output vector entries for queried range of length B_d:
  # All-chunk metadata (10 * B_d: a_index, unit_index, id, g, weight, cohort_index,
  # in_risk, is_target, in_observed, is_risk_nonobserver, chunk_to_obs_pos)
  # + observer index vectors (<= 6 * B_d) + raw endpoints if read_outcomes (<= 2 * B_d)
  output_vector_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(24, B_d, max_out, "query chunk row metadata, risk/observer flags, and chunk position map"),
      if (isTRUE(read_outcomes)) ipw_owner_safe_product(2, B_d, max_out, "raw stored y_left_obs and y_right_obs endpoints") else 0,
      256
    ),
    max_bound = max_out,
    label = "query chunk output vector entries"
  )
  output_matrix_entries <- ipw_owner_safe_sum(
    c(
      if (isTRUE(include_design)) ipw_owner_safe_product(2, Bp, max_out, "query chunk Z_chunk and Z_obs matrices") else 0,
      diff_iv_mat
    ),
    max_bound = max_out,
    label = "query chunk output matrix entries"
  )
  output_entries <- ipw_owner_safe_sum(
    c(output_vector_entries, output_matrix_entries),
    max_bound = max_out,
    label = "total query chunk output entries"
  )

  query_scratch_entries <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(64, B_d, max_me, "posting binary-bound and risk-filter scratch"),
      m,
      256
    ),
    max_bound = max_me,
    label = "query chunk scratch entries"
  )
  concurrent_entries <- ipw_owner_safe_sum(
    c(output_entries, query_scratch_entries),
    max_bound = max_me,
    label = "concurrent query chunk payload entries"
  )
  work_ops_bound <- ipw_owner_safe_sum(
    c(
      ipw_owner_safe_product(64 * (1 + ceiling(log2(m + 1))), B_d, max_wo, "chunk metadata, compact cohort lookup and posting intersection charges"),
      if (isTRUE(read_outcomes)) ipw_owner_safe_product(4, B_d, max_wo, "risk-filtered raw endpoint read and finiteness check ops") else 0,
      if (isTRUE(read_outcomes) && isTRUE(enclose_differences)) ipw_owner_safe_product(2, B_d, max_wo, "outward interval endpoint subtraction ops") else 0,
      if (isTRUE(include_design)) ipw_owner_safe_product(2, Bp, max_wo, "chunk design matrix slice ops") else 0,
      4096 # eight posting binary bounds over exact-index domain, plus scalar validation
    ),
    max_bound = max_wo,
    label = "query chunk work operations"
  )

  list(
    queried_row_count = logit_cert_count_output(B_d),
    p_covariates = p,
    m_cohorts = m,
    read_outcomes = isTRUE(read_outcomes),
    include_design = isTRUE(include_design),
    enclose_differences = isTRUE(enclose_differences),
    matrix_dim_limit = as.double(.Machine$integer.max),
    vector_length_limit = 2^52,
    exact_count_limit = 2^53 - 1,
    output_vector_entries = output_vector_entries,
    output_matrix_entries = output_matrix_entries,
    output_entries = output_entries,
    query_scratch_entries = query_scratch_entries,
    concurrent_matrix_entries = output_matrix_entries,
    concurrent_vector_entries = concurrent_entries - output_matrix_entries,
    concurrent_entries = concurrent_entries,
    concurrent_bytes_8double = 8 * concurrent_entries,
    work_ops_bound = work_ops_bound,
    full_n_scan_performed = FALSE,
    hash_temporaries_allocated = FALSE,
    risk_filtered_before_outcome_read = TRUE,
    peak_rss_promised = FALSE
  )
}

preflight_ipw_stream <- function(
  n_total_queried_rows,
  max_chunk_rows_req,
  p_covariates,
  m_cohorts,
  control,
  read_outcomes = FALSE,
  include_design = FALSE,
  enclose_differences = FALSE
) {
  ipw_owner_require_deps()
  ctrl <- validate_ipw_owner_control(control)
  ntot <- ipw_owner_int_bound(n_total_queried_rows, "n_total_queried_rows", min_val = 1)
  b_eff <- min(ntot, ipw_owner_int_bound(max_chunk_rows_req, "max_chunk_rows", min_val = 1), ctrl$max_chunk_rows)
  n_batches <- ceiling(ntot / b_eff)

  ch_pre <- preflight_ipw_range_query(
    b_len = b_eff,
    p_covariates = p_covariates,
    m_cohorts = m_cohorts,
    control = ctrl,
    read_outcomes = read_outcomes,
    include_design = include_design,
    enclose_differences = enclose_differences
  )
  total_work_ops_bound <- ipw_owner_safe_product(
    n_batches,
    ch_pre$work_ops_bound,
    max_entries = ctrl$max_work_ops,
    label = "whole-stream work operations"
  )

  list(
    queried_row_count = logit_cert_count_output(ntot),
    effective_chunk_rows = logit_cert_count_output(b_eff),
    n_batches = logit_cert_count_output(n_batches),
    chunk_preflight = ch_pre,
    concurrent_entries = ch_pre$concurrent_entries,
    concurrent_bytes_8double = ch_pre$concurrent_bytes_8double,
    total_work_ops_bound = total_work_ops_bound,
    peak_rss_promised = FALSE
  )
}

# ---------------------------------------------------------------------------
# Component 1: Immutable Common Panel Storage Owner (`cdid_ipw_panel_owner`)
# ---------------------------------------------------------------------------

construct_ipw_panel_owner <- function(
  data,
  idname,
  tname,
  yname,
  gname,
  control,
  owner_id = NULL
) {
  ipw_owner_require_deps()
  if (missing(data) || missing(idname) || missing(tname) || missing(yname) ||
      missing(gname) || missing(control)) {
    cdid_ipw_owner_abort(
      "All required arguments (`data`, `idname`, `tname`, `yname`, `gname`, `control`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }

  # 1. Validate control and lightweight headers BEFORE scanning content or copying
  hdr <- ipw_panel_owner_headers(data, idname, tname, yname, gname, control, owner_id)
  ctrl <- hdr$ctrl; cols <- hdr$cols; n <- hdr$n
  id <- hdr$id; time <- hdr$time; y <- hdr$y; g <- hdr$g

  # 2. Preflight BEFORE scanning any column values or sorting
  preflight <- preflight_ipw_panel_owner(n_rows = n, control = ctrl)
  step <- preflight$effective_chunk_rows

  # 3. Chunked raw value validation (NO global pairwise outcome subtraction!)
  id_is_num <- is.numeric(id)
  id_storage <- typeof(id)
  y_is_lgl <- is.logical(y)
  validation_chunks <- 0
  actual_ops <- 0
  id_byte_metadata <- list(total = 0, max = 0)

  r_start <- 1L
  while (r_start <= n) {
    ipw_owner_check_interrupt()
    r_end <- min(as.double(n), as.double(r_start) + step - 1)
    idx <- seq.int(r_start, r_end)
    b_len <- length(idx)

    id_ch <- id[idx]
    if (id_is_num) {
      if (anyNA(id_ch) || any(!is.finite(id_ch))) {
        cdid_ipw_owner_abort(
          "ID must be a nonmissing finite numeric or nonempty character vector.",
          c("cdid_ipw_owner_id_error", "cdid_pairs_id_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    } else {
      if (anyNA(id_ch) || any(!nzchar(id_ch))) {
        cdid_ipw_owner_abort(
          "ID must be a nonmissing finite numeric or nonempty character vector.",
          c("cdid_ipw_owner_id_error", "cdid_pairs_id_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    }

    if (!id_is_num) {
      byte_chunk <- ipw_owner_id_byte_metadata(id_ch)
      id_byte_metadata$total <- ipw_owner_safe_sum(c(id_byte_metadata$total, byte_chunk$total), label = "input ID byte count")
      id_byte_metadata$max <- max(id_byte_metadata$max, byte_chunk$max)
    }
    t_ch <- time[idx]
    if (anyNA(t_ch) || any(!is.finite(t_ch))) {
      cdid_ipw_owner_abort(
        "Time must be a finite numeric scalar column.",
        c("cdid_ipw_owner_time_error", "cdid_pairs_time_error"),
        diagnostics = list(chunk_start = r_start, chunk_end = r_end)
      )
    }

    g_ch <- g[idx]
    if (anyNA(g_ch) || any(is.nan(g_ch) | g_ch == -Inf)) {
      cdid_ipw_owner_abort(
        "Cohort must be a numeric scalar column; only +Inf is allowed for never-treated units.",
        c("cdid_ipw_owner_cohort_error", "cdid_pairs_cohort_error"),
        diagnostics = list(chunk_start = r_start, chunk_end = r_end)
      )
    }

    y_ch <- y[idx]
    if (y_is_lgl) {
      if (!all(is.na(y_ch))) {
        cdid_ipw_owner_abort(
          "Logical outcome column is only valid when all entries are NA.",
          c("cdid_ipw_owner_outcome_error", "cdid_pairs_outcome_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    } else {
      if (any(is.nan(y_ch) | is.infinite(y_ch))) {
        cdid_ipw_owner_abort(
          "Outcome must be numeric; only NA is unobserved (NaN and +/-Inf are refused).",
          c("cdid_ipw_owner_outcome_error", "cdid_pairs_outcome_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    }

    validation_chunks <- validation_chunks + 1
    actual_ops <- actual_ops + 4 * b_len
    r_start <- as.double(r_end) + 1
  }

  if (!id_is_num) {
    preflight <- ipw_owner_preflight_id_bytes(preflight, ctrl, id_byte_metadata,
      owner_max_bytes = 0, lookup_count = 0, search_width = n, simultaneous = TRUE)
  }
  # 4. Normalize canonical unit IDs, calendar, and long observation positions once
  # Numeric IDs and dates are matched and ordered directly as IEEE doubles, never via strings.
  id_key_col <- if (id_is_num) {
    d <- as.double(id)
    d[d == 0.0] <- 0.0
    d
  } else {
    as.character(id)
  }
  time_col <- as.double(time)
  time_col[time_col == 0.0] <- 0.0
  g_col <- as.double(g)
  g_col[g_col == 0.0] <- 0.0
  y_col <- if (y_is_lgl) rep(NA_real_, n) else as.double(y)

  ids_key <- sort(unique(id_key_col), method = "radix")
  nu <- length(ids_key)
  if (nu > .Machine$integer.max) cdid_ipw_owner_abort("Distinct unit count exceeds tabulate's INT_MAX bin limit.", "cdid_ipw_owner_bound_error")
  ui <- match(id_key_col, ids_key)
  ord <- order(ui, time_col, method = "radix")
  actual_ops <- actual_ops + 4 * n

  ui <- ui[ord]
  time_ord <- time_col[ord]
  y_ord <- y_col[ord]
  g_ord <- g_col[ord]

  if (n > 1L) {
    same <- ui[-1L] == ui[-n]
    if (any(same & g_ord[-1L] != g_ord[-n])) {
      cdid_ipw_owner_abort(
        "Within-unit cohort changes are invalid.",
        c("cdid_ipw_owner_cohort_change_error", "cdid_pairs_cohort_change_error", "cdid_pairs_cohort_error")
      )
    }
    if (any(same & time_ord[-1L] == time_ord[-n])) {
      cdid_ipw_owner_abort(
        "Duplicate (ID, time) rows are invalid.",
        c("cdid_ipw_owner_duplicate_error", "cdid_pairs_duplicate_error")
      )
    }
    actual_ops <- actual_ops + 2 * n
  }

  ipw_panel_owner_assemble(owner_id, ctrl, preflight, cols, n, id_is_num, id_storage, id_byte_metadata,
    ids_key, nu, ui, time_ord, y_ord, g_ord, validation_chunks, actual_ops)
}

# Control/header validation of the panel owner (step 1), shared by the
# standalone constructor and the prepared-intake constructor.
ipw_panel_owner_headers <- function(data, idname, tname, yname, gname, control, owner_id) {
  ctrl <- validate_ipw_owner_control(control)
  if (!is.null(owner_id) && !gmm_scalar_character(owner_id)) {
    cdid_ipw_owner_abort("`owner_id` must be NULL or a nonempty character scalar.", "cdid_ipw_owner_input_error")
  }
  if (!is.data.frame(data)) {
    cdid_ipw_owner_abort("`data` must be a data.frame.", c("cdid_ipw_owner_data_error", "cdid_pairs_data_error"))
  }
  ipw_owner_refuse_decorative_attrs(data, "data")

  cols <- c(idname = idname, tname = tname, yname = yname, gname = gname)
  args <- list(idname, tname, yname, gname)
  if (any(!vapply(args, function(x) is.character(x) && is.null(dim(x)) &&
                 length(x) == 1L && !is.na(x) && nzchar(x), logical(1L))) ||
      length(unique(cols)) != 4L ||
      any(vapply(cols, function(x) sum(names(data) == x, na.rm = TRUE) != 1L, logical(1L)))) {
    cdid_ipw_owner_abort(
      "Select four distinct, unambiguous existing column names (`idname`, `tname`, `yname`, `gname`).",
      c("cdid_ipw_owner_column_error", "cdid_pairs_column_error")
    )
  }

  n <- nrow(data)
  if (is.na(n) || n < 1L) {
    cdid_ipw_owner_abort("`data` must contain at least one row.", c("cdid_ipw_owner_data_error", "cdid_pairs_data_error"))
  }

  id <- data[[idname]]
  time <- data[[tname]]
  y <- data[[yname]]
  g <- data[[gname]]

  for (item in list(
    list(v = id, nm = "id"),
    list(v = time, nm = "time"),
    list(v = y, nm = "y"),
    list(v = g, nm = "g")
  )) {
    ipw_owner_refuse_decorative_attrs(item$v, item$nm)
  }

  if (!pair_scalar_vector(id, n) || is.factor(id) || is.object(id) ||
      !(is.numeric(id) || is.character(id))) {
    cdid_ipw_owner_abort(
      "ID column must be a plain non-object scalar numeric or character vector.",
      c("cdid_ipw_owner_id_error", "cdid_pairs_id_error")
    )
  }
  if (!pair_scalar_vector(time, n) || is.object(time) || !is.numeric(time)) {
    cdid_ipw_owner_abort(
      "Time column must be a plain numeric scalar vector.",
      c("cdid_ipw_owner_time_error", "cdid_pairs_time_error")
    )
  }
  if (!pair_scalar_vector(g, n) || is.object(g) || !is.numeric(g)) {
    cdid_ipw_owner_abort(
      "Cohort column must be a plain numeric scalar vector.",
      c("cdid_ipw_owner_cohort_error", "cdid_pairs_cohort_error")
    )
  }
  if (!pair_scalar_vector(y, n) || is.object(y) || !(is.numeric(y) || is.logical(y))) {
    cdid_ipw_owner_abort(
      "Outcome column must be a plain scalar numeric (or all-NA logical) vector.",
      c("cdid_ipw_owner_outcome_error", "cdid_pairs_outcome_error")
    )
  }

  list(ctrl = ctrl, cols = cols, n = n, id = id, time = time, y = y, g = g)
}

# Prepared-intake panel owner (CP08 common boundary). The public request layer
# (`chained_request_prepare_inputs`) has already established, on exactly this
# `data` object, every fact the chunked scan below re-derives: finite numeric
# identifiers and dates, cohorts without NaN/-Inf, outcomes without NaN/Inf,
# one row per (unit, date) and one cohort per unit. Its sealed record carries
# the first-appearance unit positions, the sorted owner order and the canonical
# (unit, date) row permutation. Headers and the preflight admission run exactly
# as in the standalone constructor; the record is accepted only when it is the
# sealed object describing this very `data` (object identity, never a flag),
# and the retained representation is then assembled and sealed by the same
# code as the standalone owner, with the same diagnostic charges.
construct_ipw_panel_owner_prepared <- function(data, idname, tname, yname, gname, control, prepared, owner_id = NULL) {
  ipw_owner_require_deps()
  if (missing(data) || missing(idname) || missing(tname) || missing(yname) ||
      missing(gname) || missing(control) || missing(prepared)) {
    cdid_ipw_owner_abort(
      "All required arguments (`data`, `idname`, `tname`, `yname`, `gname`, `control`, `prepared`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }
  hdr <- ipw_panel_owner_headers(data, idname, tname, yname, gname, control, owner_id)
  ctrl <- hdr$ctrl; cols <- hdr$cols; n <- hdr$n
  id <- hdr$id; time <- hdr$time; y <- hdr$y; g <- hdr$g
  if (!inherits(prepared, "cdid_request_prepared_intake") || !is.list(prepared) ||
      !ipw_owner_seal_ok(prepared, "prepared_intake") || !gmm_sexp_aliases(prepared$source, data) ||
      !identical(unname(prepared$cols), unname(cols)) || !identical(prepared$n_rows, n) ||
      !is.numeric(id) || is.logical(y) || !identical(prepared$id_storage, typeof(id))) {
    cdid_ipw_owner_abort(
      "Prepared intake record does not describe the supplied `data`.",
      "cdid_ipw_owner_input_error"
    )
  }

  # 2. Preflight BEFORE any ordering or allocation (same admission as the standalone owner)
  preflight <- preflight_ipw_panel_owner(n_rows = n, control = ctrl)
  step <- preflight$effective_chunk_rows
  id_storage <- typeof(id)
  # The same validation charges the standalone chunk scan records: the facts
  # were established by whole-vector passes in the request layer on this data.
  validation_chunks <- ceiling(n / step)
  actual_ops <- 4 * n
  id_byte_metadata <- list(total = 0, max = 0)

  # 4. Canonical unit IDs, calendar and long positions from the record's maps
  time_col <- as.double(time)
  time_col[time_col == 0.0] <- 0.0
  g_col <- as.double(g)
  g_col[g_col == 0.0] <- 0.0
  y_col <- as.double(y)

  ids_key <- prepared$unit_key
  nu <- prepared$n_units
  if (nu > .Machine$integer.max) cdid_ipw_owner_abort("Distinct unit count exceeds tabulate's INT_MAX bin limit.", "cdid_ipw_owner_bound_error")
  ui <- prepared$unit_rank[prepared$pos_unit]
  ord <- prepared$row_order
  actual_ops <- actual_ops + 4 * n

  ui <- ui[ord]
  time_ord <- time_col[ord]
  y_ord <- y_col[ord]
  g_ord <- g_col[ord]
  if (n > 1L) actual_ops <- actual_ops + 2 * n

  ipw_panel_owner_assemble(owner_id, ctrl, preflight, cols, n, TRUE, id_storage, id_byte_metadata,
    ids_key, nu, ui, time_ord, y_ord, g_ord, validation_chunks, actual_ops)
}

# Retained representation of the panel owner from its canonical ordered arrays
# (shared by both constructors; unchanged fields, counts, postings and seal).
ipw_panel_owner_assemble <- function(owner_id, ctrl, preflight, cols, n, id_is_num, id_storage, id_byte_metadata,
                                     ids_key, nu, ui, time_ord, y_ord, g_ord, validation_chunks, actual_ops) {
  # Preserve original integer storage if input ID was integer, else canonical double/character
  ids <- if (identical(id_storage, "integer")) as.integer(ids_key) else ids_key
  all_times <- sort(unique(time_ord), method = "radix")
  n_cal <- length(all_times)
  if (n_cal > .Machine$integer.max) cdid_ipw_owner_abort("Distinct calendar count exceeds tabulate's INT_MAX bin limit.", "cdid_ipw_owner_bound_error")

  observed <- !is.na(y_ord)
  ou <- ui[observed]
  ot <- time_ord[observed]
  oy <- y_ord[observed]
  ocal <- ipw_owner_match_sorted_double(ot, all_times)

  # Explicit NA rows (preserved separately from missing rows for exact status indexing)
  na_mask <- !observed
  na_u <- ui[na_mask]
  na_t <- time_ord[na_mask]
  na_cal <- ipw_owner_match_sorted_double(na_t, all_times)

  nr <- tabulate(ui, nbins = nu)
  no <- tabulate(ou, nbins = nu)
  n_na_u <- nr - no
  counts <- pair_count_safe(no)

  obs_end <- cumsum(as.double(no))
  obs_start <- obs_end - no + 1
  na_end <- cumsum(as.double(n_na_u))
  na_start <- na_end - n_na_u + 1

  ug <- g_ord[match(seq_len(nu), ui)]
  first <- last <- sy <- singleton_time <- rep(NA_real_, nu)
  nonzero <- which(no > 0L)
  single <- which(no == 1L)
  unobs <- which(no == 0L)
  a_unit_index <- which(no >= 2L)
  n_a <- length(a_unit_index)

  first[nonzero] <- ot[obs_start[nonzero]]
  last[nonzero] <- ot[obs_end[nonzero]]
  sy[single] <- oy[obs_start[single]]
  singleton_time[single] <- first[single]

  unit_to_a <- integer(nu)
  if (n_a > 0L) {
    unit_to_a[a_unit_index] <- seq_len(n_a)
  }
  obs_a_index <- unit_to_a[ou]

  # Build compact calendar postings over A-units once: O(Q_A + T)
  keep_a_obs <- which(obs_a_index > 0L)
  q_a <- length(keep_a_obs)
  cal_a <- ocal[keep_a_obs]
  ord_cal <- if (q_a > 0L) order(cal_a, method = "radix") else integer(0L)
  date_count <- tabulate(cal_a, nbins = n_cal)
  date_end <- cumsum(date_count)
  date_start <- date_end - date_count + 1L
  post_a_index <- obs_a_index[keep_a_obs[ord_cal]]
  post_unit_index <- ou[keep_a_obs[ord_cal]]
  post_obs_index <- keep_a_obs[ord_cal]
  actual_ops <- actual_ops + 2 * length(ot) + 2 * q_a

  a_id <- ids[a_unit_index]
  a_id_key <- ids_key[a_unit_index]
  a_g <- as.double(ug[a_unit_index])
  a_cohorts <- if (n_a > 0L) sort(unique(a_g), method = "radix") else numeric(0L)
  a_cohort_index <- if (n_a > 0L) ipw_owner_match_sorted_double(a_g, a_cohorts) else integer(0L)
  cohort_a_count <- if (length(a_cohorts) > 0L) tabulate(a_cohort_index, nbins = length(a_cohorts)) else integer(0L)

  units <- data.frame(
    unit_index = seq_len(nu),
    id = ids,
    g = ug,
    n_rows = nr,
    n_obs = no,
    n_na_rows = n_na_u,
    n_pairs = counts$by_unit,
    is_unobserved = no == 0L,
    is_singleton = no == 1L,
    pair_eligible = no >= 2L,
    first_obs_time = first,
    last_obs_time = last,
    singleton_time = singleton_time,
    singleton_y = sy,
    obs_start = obs_start,
    obs_end = obs_end,
    na_start = na_start,
    na_end = na_end,
    stringsAsFactors = FALSE
  )
  singletons <- data.frame(
    unit_index = single,
    id = ids[single],
    g = ug[single],
    time = first[single],
    y = sy[single],
    n_rows = nr[single],
    stringsAsFactors = FALSE
  )
  unobserved_units <- data.frame(
    unit_index = unobs,
    id = ids[unobs],
    g = ug[unobs],
    n_rows = nr[unobs],
    n_na_rows = n_na_u[unobs],
    stringsAsFactors = FALSE
  )
  observations <- data.frame(
    time = ot,
    y = oy,
    calendar_index = ocal,
    unit_index = ou,
    a_index = obs_a_index,
    stringsAsFactors = FALSE
  )
  explicit_na_rows <- data.frame(
    unit_index = na_u,
    time = na_t,
    calendar_index = na_cal,
    stringsAsFactors = FALSE
  )
  obs_times <- if (length(ot) > 0L) sort(unique(ot), method = "radix") else numeric(0L)
  all_cohorts <- sort(unique(ug), method = "radix")

  actual_retained_entries <- ipw_owner_safe_sum(
    c(
      20 * nu, # 18 unit columns plus unit_id_key and unit_to_a
      6 * length(single),
      5 * length(unobs),
      5 * length(ot), # observations table
      3 * length(na_u), # explicit_na_rows table
      5 * n_a, # a_unit_index, a_id, a_id_key, a_g, a_cohort_index
      2 * length(a_cohorts),
      n_cal + length(obs_times) + length(all_cohorts),
      3 * n_cal + 3 * q_a, # postings
      256
    ),
    max_bound = ctrl$max_retained_entries,
    label = "actual retained panel owner entries"
  )
  stopifnot(
    actual_retained_entries <= preflight$total_retained_entries,
    actual_ops <= preflight$work_ops_bound
  )

  panel_obj <- structure(
    list(
      owner_id = owner_id,
      id_type = if (id_is_num) "numeric" else "character",
      id_storage = id_storage,
      id_byte_metadata = id_byte_metadata,
      unit_id_key = ids_key,
      n_a_units = n_a,
      n_total_units = nu,
      n_singleton_units = length(single),
      n_unobserved_units = length(unobs),
      n_calendar_dates = n_cal,
      n_observed_dates = length(obs_times),
      n_observed_rows = length(ot),
      n_a_observed_rows = q_a,
      n_explicit_na_rows = length(na_u),
      a_unit_index = a_unit_index,
      unit_to_a = unit_to_a,
      a_id = a_id,
      a_id_key = a_id_key,
      a_g = a_g,
      a_cohorts = a_cohorts,
      a_cohort_index = a_cohort_index,
      cohort_a_count = cohort_a_count,
      units = units,
      singletons = singletons,
      unobserved_units = unobserved_units,
      observations = observations,
      explicit_na_rows = explicit_na_rows,
      times = obs_times,
      all_times = all_times,
      cohorts = all_cohorts,
      postings = list(
        date_start = date_start,
        date_end = date_end,
        date_count = date_count,
        post_a_index = post_a_index,
        post_unit_index = post_unit_index,
        post_obs_index = post_obs_index
      ),
      counts = list(
        n_input_rows = n,
        n_observed_rows = length(ot),
        n_na_rows = sum(!observed),
        n_units = nu,
        n_unobserved_units = length(unobs),
        n_singleton_units = length(single),
        n_pair_eligible_units = n_a,
        n_pairs = counts$total
      ),
      cols = cols,
      control = ctrl,
      preflight = preflight,
      accounting = list(
        actual_retained_vector_entries = actual_retained_entries,
        actual_retained_matrix_entries = 0,
        actual_retained_entries = actual_retained_entries,
        raw_outcome_validation_performed = TRUE,
        global_pair_subtraction_performed = FALSE,
        endpoints_read_during_ingestion = 0L # pair queries; raw y validation above is required
      ),
      work = list(
        validation_chunks_processed = logit_cert_count_output(validation_chunks),
        actual_work_ops = actual_ops
      )
    ),
    class = "cdid_ipw_panel_owner"
  )

  ipw_owner_seal(panel_obj, "panel_owner")
}

construct_common_panel_owner <- construct_ipw_panel_owner
chained_ipw_panel_owner <- construct_ipw_panel_owner

validate_ipw_panel_owner <- function(panel_owner) {
  ipw_owner_require_deps()
  if (missing(panel_owner) || is.null(panel_owner) ||
      !inherits(panel_owner, "cdid_ipw_panel_owner") ||
      !is.list(panel_owner) || !ipw_owner_seal_ok(panel_owner, "panel_owner")) {
    cdid_ipw_owner_abort(
      "Panel owner is malformed, foreign, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_panel_error")
    )
  }
  bm <- panel_owner$id_byte_metadata
  if (!is.list(bm) || !identical(names(bm), c("total", "max")) ||
      !gmm_scalar_numeric(bm$total, nonneg = TRUE) || !gmm_scalar_numeric(bm$max, nonneg = TRUE) ||
      bm$total != floor(bm$total) || bm$max != floor(bm$max) || bm$total < bm$max || bm$total > 2^53 - 1) {
    cdid_ipw_owner_abort("Owner lacks verified fixed-size character byte metadata.", "cdid_ipw_owner_integrity_error")
  }
  invisible(panel_owner)
}

validate_common_panel_owner <- validate_ipw_panel_owner

# Distinguish observed outcome, explicit missing outcome (`y = NA` row present),
# and missing row (no `(id, time)` record in input data) without reading `y`.
query_ipw_panel_cell_status <- function(panel_owner, id, time, control) {
  validate_ipw_panel_owner(panel_owner)
  ctrl <- validate_ipw_owner_control(control)
  if (missing(id) || missing(time) || is.null(id) || is.null(time) ||
      is.factor(id) || is.object(id) || is.object(time) ||
      !is.atomic(id) || !is.null(dim(id)) || !is.atomic(time) || !is.null(dim(time)) ||
      length(id) != length(time) || length(id) < 1L) {
    cdid_ipw_owner_abort(
      "`id` and `time` must be non-empty plain scalar vectors of equal length.",
      "cdid_ipw_owner_input_error"
    )
  }
  k <- length(id)
  if (k > ctrl$max_chunk_rows || 8 * k + 128 > ctrl$max_output_entries ||
      32 * k + 128 > ctrl$max_matrix_entries ||
      256 * k * (1 + ceiling(log2(panel_owner$counts$n_input_rows + 1))) > ctrl$max_work_ops) {
    cdid_ipw_owner_abort(
      "Cell status query exceeds explicit `control` bounds.",
      "cdid_ipw_owner_bound_error"
    )
  }
  id_is_num <- is.numeric(id)
  id_is_chr <- is.character(id)
  if ((identical(panel_owner$id_type, "numeric") && (!id_is_num || anyNA(id) || any(!is.finite(id)))) ||
      (identical(panel_owner$id_type, "character") && (!id_is_chr || anyNA(id) || any(!nzchar(id)))) ||
      !is.numeric(time) || anyNA(time) || any(!is.finite(time))) {
    cdid_ipw_owner_abort(
      "`id` must match `panel_owner$id_type` without coercion and `time` must be finite numeric.",
      "cdid_ipw_owner_input_error"
    )
  }

  if (id_is_chr) {
    byte_pre <- list(concurrent_entries = 32 * k + 128,
      concurrent_vector_entries = 32 * k + 128,
      work_ops_bound = 256 * k * (1 + ceiling(log2(panel_owner$counts$n_input_rows + 1))))
    byte_pre <- ipw_owner_preflight_id_bytes(byte_pre, ctrl, ipw_owner_id_byte_metadata(id),
      owner_max_bytes = panel_owner$id_byte_metadata$max, lookup_count = k,
      search_width = panel_owner$n_total_units, simultaneous = FALSE)
  }

  id_key <- if (id_is_num) {
    d <- as.double(id)
    d[d == 0.0] <- 0.0
    d
  } else {
    as.character(id)
  }
  t_dbl <- as.double(time)
  t_dbl[t_dbl == 0.0] <- 0.0

  u_pos <- ipw_owner_match_sorted_id(id_key, panel_owner$unit_id_key, panel_owner$id_type)
  cal_pos <- ipw_owner_match_sorted_double(t_dbl, panel_owner$all_times)
  status <- rep("missing_row", k)
  obs_row <- rep(NA_integer_, k)
  na_row <- rep(NA_integer_, k)

  units <- panel_owner$units
  obs_t <- panel_owner$observations$time
  na_t <- panel_owner$explicit_na_rows$time

  for (i in seq_len(k)) {
    u <- u_pos[i]
    if (is.na(u)) next
    ti <- t_dbl[i]
    if (units$n_obs[u] > 0L) {
      s_o <- units$obs_start[u]
      e_o <- units$obs_end[u]
      hit_o <- ipw_owner_match_sorted_double(ti, obs_t, s_o, e_o)
      if (!is.na(hit_o)) {
        status[i] <- "observed"
        obs_row[i] <- hit_o
        next
      }
    }
    if (units$n_na_rows[u] > 0L) {
      s_n <- units$na_start[u]
      e_n <- units$na_end[u]
      hit_n <- ipw_owner_match_sorted_double(ti, na_t, s_n, e_n)
      if (!is.na(hit_n)) {
        status[i] <- "explicit_na"
        na_row[i] <- hit_n
      }
    }
  }

  data.frame(
    id = id,
    time = t_dbl,
    unit_index = u_pos,
    calendar_index = cal_pos,
    in_known_calendar = !is.na(cal_pos),
    status = status,
    obs_index = obs_row,
    explicit_na_index = na_row,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Component 2: Explicit Frame / Design Owner (`cdid_ipw_frame_owner`)
# ---------------------------------------------------------------------------

construct_ipw_frame_owner <- function(
  panel_owner,
  id,
  Z,
  weights = NULL,
  control,
  frame_g = NULL,
  allow_frame_only = FALSE,
  frame_id = NULL,
  design_on_A_only = FALSE
) {
  ipw_owner_require_deps()
  if (missing(panel_owner) || missing(id) || missing(Z) || missing(control)) {
    cdid_ipw_owner_abort(
      "All required arguments (`panel_owner`, `id`, `Z`, `control`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }

  hdr <- ipw_frame_owner_headers(panel_owner, id, Z, weights, control, frame_g, allow_frame_only, frame_id, design_on_A_only)
  ctrl <- hdr$ctrl; n_a <- hdr$n_a; nsup <- hdr$nsup; p <- hdr$p
  id_is_num <- hdr$id_is_num; id_is_chr <- hdr$id_is_chr; weights_supplied <- hdr$weights_supplied

  # 2. Preflight BEFORE scanning id, Z, weights, frame_g
  preflight <- preflight_ipw_frame_owner(
    n_a_units = n_a,
    n_total_units = panel_owner$n_total_units,
    n_supplied_id = nsup,
    p_covariates = p,
    n_calendar_dates = panel_owner$n_calendar_dates,
    n_a_observed_rows = panel_owner$n_a_observed_rows,
    control = ctrl
  )
  step <- preflight$effective_chunk_rows

  # 3. Validate id, Z, weights, frame_g in bounded chunks
  validation_chunks <- 0
  actual_ops <- 0
  id_byte_metadata <- list(total = 0, max = 0)
  r_start <- 1L
  while (r_start <= nsup) {
    ipw_owner_check_interrupt()
    r_end <- min(as.double(nsup), as.double(r_start) + step - 1)
    idx <- seq.int(r_start, r_end)
    b_len <- length(idx)

    id_ch <- id[idx]
    if (id_is_num) {
      if (anyNA(id_ch) || any(!is.finite(id_ch))) {
        cdid_ipw_owner_abort(
          "`id` contains NA, NaN, or infinite numeric identifiers.",
          c("cdid_ipw_owner_id_error", "cdid_ipw_id_error", "cdid_ipw_context_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    } else {
      if (anyNA(id_ch) || any(!nzchar(id_ch))) {
        cdid_ipw_owner_abort(
          "`id` contains NA or empty character identifiers.",
          c("cdid_ipw_owner_id_error", "cdid_ipw_id_error", "cdid_ipw_context_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    }

    if (id_is_chr) {
      byte_chunk <- ipw_owner_id_byte_metadata(id_ch)
      id_byte_metadata$total <- ipw_owner_safe_sum(c(id_byte_metadata$total, byte_chunk$total), label = "frame ID byte count")
      id_byte_metadata$max <- max(id_byte_metadata$max, byte_chunk$max)
    }
    zi <- idx
    if (design_on_A_only) {
      # Sorted owner lookup per chunk; never rescan the complete A roster.
      up <- ipw_owner_match_sorted_id(id[idx], panel_owner$unit_id_key, panel_owner$id_type)
      hit <- which(!is.na(up)); active <- rep(FALSE, length(idx))
      active[hit] <- panel_owner$unit_to_a[up[hit]] > 0L
      zi <- idx[active]
    }
    Z_ch <- Z[zi, , drop = FALSE]
    if (anyNA(Z_ch) || any(!is.finite(Z_ch))) {
      cdid_ipw_owner_abort(
        "`Z` must be a finite numeric matrix without NA, NaN, or Inf.",
        c("cdid_ipw_owner_covariate_error", "cdid_ipw_covariate_error", "cdid_ipw_context_error"),
        diagnostics = list(chunk_start = r_start, chunk_end = r_end)
      )
    }

    if (weights_supplied) {
      w_ch <- weights[idx]
      if (anyNA(w_ch) || any(!is.finite(w_ch) | w_ch <= 0)) {
        cdid_ipw_owner_abort(
          "`weights` must be a finite, strictly positive numeric vector of length `length(id)`.",
          c("cdid_ipw_owner_weight_error", "cdid_ipw_weight_error", "cdid_ipw_context_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    }

    if (!is.null(frame_g)) {
      fg_ch <- frame_g[idx]
      if (anyNA(fg_ch) || any(is.nan(fg_ch) | fg_ch == -Inf)) {
        cdid_ipw_owner_abort(
          "`frame_g` must be finite numeric or `+Inf` (for never-treated units).",
          c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_error"),
          diagnostics = list(chunk_start = r_start, chunk_end = r_end)
        )
      }
    }

    validation_chunks <- validation_chunks + 1
    actual_ops <- actual_ops + b_len * (p + 3)
    r_start <- as.double(r_end) + 1
  }

  if (id_is_chr) {
    preflight <- ipw_owner_preflight_id_bytes(preflight, ctrl, id_byte_metadata,
      owner_max_bytes = panel_owner$id_byte_metadata$max, lookup_count = nsup,
      search_width = panel_owner$n_total_units, simultaneous = FALSE)
  }

  id_key <- if (id_is_num) {
    d <- as.double(id)
    d[d == 0.0] <- 0.0
    d
  } else {
    as.character(id)
  }
  if (anyDuplicated(id_key)) {
    cdid_ipw_owner_abort(
      "Duplicate unit identifiers in `id` are invalid.",
      c("cdid_ipw_owner_id_error", "cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }

  unit_pos_of_supplied <- ipw_owner_match_sorted_id(id_key, panel_owner$unit_id_key, panel_owner$id_type)
  absent_from_panel <- which(is.na(unit_pos_of_supplied))
  if (length(absent_from_panel) > 0L && !isTRUE(allow_frame_only)) {
    cdid_ipw_owner_abort(
      "`id` contains unit identifiers absent from `panel_owner$units$id` (set `allow_frame_only = TRUE` only when explicitly retaining external frame-only non-A metadata).",
      c("cdid_ipw_owner_id_error", "cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }

  # If frame_g is supplied, verify cohort invariance on all units present in panel_owner
  in_panel_sup <- which(!is.na(unit_pos_of_supplied))
  if (!is.null(frame_g) && length(in_panel_sup) > 0L) {
    fg_sub <- as.double(frame_g[in_panel_sup])
    fg_sub[fg_sub == 0.0] <- 0.0
    pg_sub <- panel_owner$units$g[unit_pos_of_supplied[in_panel_sup]]
    if (!identical(fg_sub, pg_sub)) {
      cdid_ipw_owner_abort(
        "Supplied `frame_g` disagrees with `panel_owner$units$g` on panel units.",
        c("cdid_ipw_owner_cohort_change_error", "cdid_ipw_owner_cohort_error")
      )
    }
  }

  # Every canonical A-unit (`n_obs >= 2`) must be present in supplied `id`
  a_pos_in_supplied <- match(panel_owner$a_id_key, id_key)
  if (anyNA(a_pos_in_supplied)) {
    cdid_ipw_owner_abort(
      "Every pair-eligible A-unit (`n_obs >= 2`) must be present in `id`.",
      c("cdid_ipw_owner_id_error", "cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }
  actual_ops <- actual_ops + 3 * nsup

  ipw_frame_owner_assemble(panel_owner, id, Z, weights, frame_g, ctrl, preflight, frame_id, n_a, nsup, p,
    weights_supplied, id_byte_metadata, unit_pos_of_supplied, absent_from_panel, a_pos_in_supplied,
    validation_chunks, actual_ops)
}

# Panel-owner, control and header validation of the frame owner (step 1), shared
# by the standalone constructor and the prepared-intake constructor.
ipw_frame_owner_headers <- function(panel_owner, id, Z, weights, control, frame_g, allow_frame_only, frame_id,
                                    design_on_A_only) {
  # 1. Validate panel owner, control, and lightweight headers BEFORE scanning content
  validate_ipw_panel_owner(panel_owner)
  ctrl <- validate_ipw_owner_control(control)
  if (!gmm_scalar_logical(design_on_A_only))
    cdid_ipw_owner_abort("`design_on_A_only` must be logical.", "cdid_ipw_owner_input_error")
  if (!gmm_scalar_logical(allow_frame_only)) {
    cdid_ipw_owner_abort("`allow_frame_only` must be a scalar logical.", "cdid_ipw_owner_input_error")
  }
  if (!is.null(frame_id) && !gmm_scalar_character(frame_id)) {
    cdid_ipw_owner_abort("`frame_id` must be NULL or a nonempty character scalar.", "cdid_ipw_owner_input_error")
  }
  for (item in list(
    list(v = id, nm = "id"),
    list(v = Z, nm = "Z"),
    list(v = weights, nm = "weights"),
    list(v = frame_g, nm = "frame_g")
  )) {
    ipw_owner_refuse_decorative_attrs(item$v, item$nm)
  }

  n_a <- panel_owner$n_a_units
  if (n_a < 1L) {
    cdid_ipw_owner_abort(
      "At least one pair-eligible A-unit (`n_obs >= 2`) is required in `panel_owner`.",
      c("cdid_ipw_owner_support_error", "cdid_ipw_support_error", "cdid_ipw_context_error")
    )
  }

  if (is.null(id) || is.factor(id) || is.object(id) || !is.atomic(id) ||
      !is.null(dim(id)) || length(id) < 1L) {
    cdid_ipw_owner_abort(
      "`id` must be a non-empty plain scalar vector of unit identifiers.",
      c("cdid_ipw_owner_id_error", "cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }
  nsup <- length(id)
  id_is_num <- is.numeric(id)
  id_is_chr <- is.character(id)
  panel_is_num <- identical(panel_owner$id_type, "numeric")
  if ((panel_is_num && !id_is_num) || (!panel_is_num && !id_is_chr)) {
    cdid_ipw_owner_abort(
      "`id` must match the exact type (numeric or character) of `panel_owner$units$id` without string coercion.",
      c("cdid_ipw_owner_id_error", "cdid_ipw_id_error", "cdid_ipw_context_error")
    )
  }

  if (is.null(Z) || !is.matrix(Z) || is.object(Z) || !is.numeric(Z) ||
      nrow(Z) != nsup || ncol(Z) < 1L) {
    cdid_ipw_owner_abort(
      "`Z` must be a plain numeric matrix with `nrow(Z) == length(id)` and `ncol(Z) >= 1`.",
      c("cdid_ipw_owner_covariate_error", "cdid_ipw_covariate_error", "cdid_ipw_context_error")
    )
  }
  p <- ncol(Z)

  weights_supplied <- !is.null(weights)
  if (weights_supplied) {
    if (!is.atomic(weights) || !is.null(dim(weights)) || is.object(weights) ||
        !is.numeric(weights) || length(weights) != nsup) {
      cdid_ipw_owner_abort(
        "`weights` must be a finite, strictly positive numeric vector of length `length(id)`.",
        c("cdid_ipw_owner_weight_error", "cdid_ipw_weight_error", "cdid_ipw_context_error")
      )
    }
  }
  if (!is.null(frame_g)) {
    if (!is.atomic(frame_g) || !is.null(dim(frame_g)) || is.object(frame_g) ||
        !is.numeric(frame_g) || length(frame_g) != nsup) {
      cdid_ipw_owner_abort(
        "`frame_g` must be NULL or a numeric vector of length `length(id)`.",
        c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_error")
      )
    }
  }

  list(ctrl = ctrl, n_a = n_a, nsup = nsup, p = p, id_is_num = id_is_num, id_is_chr = id_is_chr,
    weights_supplied = weights_supplied)
}

# Prepared-intake frame owner (CP08 common boundary). The request layer has
# validated, on exactly these arrays, the complete F = A roster (every frame
# unit is a pair-eligible panel unit, no frame-only unit), finite numeric
# identifiers, the finite design, unit weights and the cohorts that agree with
# the panel by identity; the sealed record carries each frame row's sorted
# owner position and the frame position of every sorted A-unit. Headers and
# the preflight admission run exactly as in the standalone constructor; the
# chunk scan, the sorted identity lookup, the cohort agreement comparison and
# the A-roster match are replaced by the record's maps after the arrays are
# recognised by identity (the panel owner must be the one assembled from the
# same record: its sorted keys are the record's own vector). Representation
# and seal come from the shared assembler with the same diagnostic charges.
construct_ipw_frame_owner_prepared <- function(
  panel_owner,
  id,
  Z,
  weights,
  control,
  frame_g,
  prepared,
  frame_id = NULL,
  design_on_A_only = FALSE
) {
  ipw_owner_require_deps()
  if (missing(panel_owner) || missing(id) || missing(Z) || missing(weights) || missing(control) ||
      missing(frame_g) || missing(prepared)) {
    cdid_ipw_owner_abort(
      "All required arguments (`panel_owner`, `id`, `Z`, `weights`, `control`, `frame_g`, `prepared`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }
  hdr <- ipw_frame_owner_headers(panel_owner, id, Z, weights, control, frame_g, FALSE, frame_id, design_on_A_only)
  ctrl <- hdr$ctrl; n_a <- hdr$n_a; nsup <- hdr$nsup; p <- hdr$p
  if (!inherits(prepared, "cdid_request_prepared_intake") || !is.list(prepared) ||
      !ipw_owner_seal_ok(prepared, "prepared_intake") ||
      !gmm_sexp_aliases(panel_owner$unit_id_key, prepared$unit_key) ||
      !gmm_sexp_aliases(id, prepared$frame_ids) || !gmm_sexp_aliases(Z, prepared$frame_Z) ||
      !gmm_sexp_aliases(weights, prepared$frame_w) || !gmm_sexp_aliases(frame_g, prepared$frame_g) ||
      !identical(prepared$n_frame, nsup) || !identical(n_a, panel_owner$n_total_units) || !hdr$id_is_num ||
      !isTRUE(hdr$weights_supplied)) {
    cdid_ipw_owner_abort(
      "Prepared intake record does not describe the supplied frame arrays and panel owner.",
      "cdid_ipw_owner_input_error"
    )
  }

  # 2. Preflight BEFORE any lookup or allocation (same admission as the standalone owner)
  preflight <- preflight_ipw_frame_owner(
    n_a_units = n_a,
    n_total_units = panel_owner$n_total_units,
    n_supplied_id = nsup,
    p_covariates = p,
    n_calendar_dates = panel_owner$n_calendar_dates,
    n_a_observed_rows = panel_owner$n_a_observed_rows,
    control = ctrl
  )
  step <- preflight$effective_chunk_rows
  # The same charges the standalone chunk scan and alignment record.
  validation_chunks <- ceiling(nsup / step)
  actual_ops <- nsup * (p + 3)
  id_byte_metadata <- list(total = 0, max = 0)

  # 3. Alignment maps from the sealed record: frame row -> sorted unit, sorted A-unit -> frame row
  unit_pos_of_supplied <- prepared$frame_unit_sorted
  absent_from_panel <- integer(0L)
  a_pos_in_supplied <- prepared$a_pos_in_frame
  actual_ops <- actual_ops + 3 * nsup

  ipw_frame_owner_assemble(panel_owner, id, Z, weights, frame_g, ctrl, preflight, frame_id, n_a, nsup, p,
    TRUE, id_byte_metadata, unit_pos_of_supplied, absent_from_panel, a_pos_in_supplied,
    validation_chunks, actual_ops)
}

# Retained representation of the frame owner from its alignment maps (shared by
# both constructors; unchanged fields, aliasing rules, accounting and seal).
ipw_frame_owner_assemble <- function(panel_owner, id, Z, weights, frame_g, ctrl, preflight, frame_id, n_a, nsup, p,
                                     weights_supplied, id_byte_metadata, unit_pos_of_supplied, absent_from_panel,
                                     a_pos_in_supplied, validation_chunks, actual_ops) {
  in_panel_sup <- which(!is.na(unit_pos_of_supplied))
  # Share input Z and weights under R copy-on-modify when already in canonical A order
  is_canonical_a_input <- (nsup == n_a) && identical(a_pos_in_supplied, seq_len(n_a))
  if (is_canonical_a_input && is.double(Z) && is.null(rownames(Z))) {
    Z_a <- Z
  } else {
    Z_a <- Z[a_pos_in_supplied, , drop = FALSE]
    storage.mode(Z_a) <- "double"
    rownames(Z_a) <- NULL
  }

  if (!weights_supplied) {
    w_all <- rep(1.0, nsup)
    w_a <- if (is_canonical_a_input) w_all else rep(1.0, n_a)
  } else {
    w_all <- if (is.double(weights) && is.null(attributes(weights))) weights else as.double(weights)
    w_a <- if (is_canonical_a_input) w_all else as.double(w_all[a_pos_in_supplied])
  }

  sum_w_a <- sum(w_a)
  if (!is.finite(sum_w_a) || sum_w_a <= 0) {
    cdid_ipw_owner_abort(
      "Total A-unit mechanical case weight must be finite and positive.",
      c("cdid_ipw_owner_weight_error", "cdid_ipw_weight_error", "cdid_ipw_context_error")
    )
  }

  a_cohorts <- panel_owner$a_cohorts
  a_cohort_index <- panel_owner$a_cohort_index
  # Stable grouping preserves each old sum's exact input order, avoiding N scans
  # for every global cohort (which may outnumber dates and approach N).
  cohort_order <- order(a_cohort_index, method = "radix")
  cohort_end <- cumsum(as.double(panel_owner$cohort_a_count))
  cohort_start <- cohort_end - panel_owner$cohort_a_count + 1
  cohort_a_weight <- vapply(seq_along(a_cohorts), function(j) {
    sum(w_a[cohort_order[seq.int(cohort_start[j], cohort_end[j])]])
  }, numeric(1L))
  if (any(!is.finite(cohort_a_weight) | cohort_a_weight <= 0)) {
    cdid_ipw_owner_abort(
      "Cohort A-unit weight sums must be finite and positive.",
      c("cdid_ipw_owner_weight_error", "cdid_ipw_weight_error", "cdid_ipw_context_error")
    )
  }
  cohort_a_share <- cohort_a_weight / n_a

  # Preserve singleton, panel-unobserved, and optional frame-only unit metadata
  nu <- panel_owner$n_total_units
  unit_weights_all <- rep(NA_real_, nu)
  if (length(in_panel_sup) > 0L) {
    unit_weights_all[unit_pos_of_supplied[in_panel_sup]] <- w_all[in_panel_sup]
  }
  actual_ops <- actual_ops + 2 * n_a + nu

  frame_only_units <- data.frame(
    supplied_pos = absent_from_panel,
    id = id[absent_from_panel],
    g = if (!is.null(frame_g)) as.double(frame_g[absent_from_panel]) else rep(NA_real_, length(absent_from_panel)),
    weight = w_all[absent_from_panel],
    n_rows = rep(0L, length(absent_from_panel)),
    n_obs = rep(0L, length(absent_from_panel)),
    stringsAsFactors = FALSE
  )

  stopifnot(actual_ops <= preflight$work_ops_bound)

  # Frame owner aliases panel_owner's long arrays directly (zero duplicate long storage)
  frame_obj <- structure(
    list(
      frame_id = frame_id,
      panel_ref = panel_owner$ref,
      id_type = panel_owner$id_type,
      id_storage = panel_owner$id_storage,
      id_byte_metadata = id_byte_metadata,
      n_a_units = n_a,
      n_total_units = nu,
      n_singleton_units = panel_owner$n_singleton_units,
      n_unobserved_units = panel_owner$n_unobserved_units,
      n_frame_only_units = length(absent_from_panel),
      n_covariates = p,
      n_calendar_dates = panel_owner$n_calendar_dates,
      n_observed_dates = panel_owner$n_observed_dates,
      n_observed_rows = panel_owner$n_observed_rows,
      n_a_observed_rows = panel_owner$n_a_observed_rows,
      n_explicit_na_rows = panel_owner$n_explicit_na_rows,
      a_unit_index = panel_owner$a_unit_index,
      unit_to_a = panel_owner$unit_to_a,
      a_id = panel_owner$a_id,
      a_id_key = panel_owner$a_id_key,
      a_g = panel_owner$a_g,
      a_cohorts = a_cohorts,
      a_cohort_index = a_cohort_index,
      cohort_a_count = panel_owner$cohort_a_count,
      cohort_a_weight = cohort_a_weight,
      cohort_a_share = cohort_a_share,
      Z_a = Z_a,
      w_a = w_a,
      sum_w_a = sum_w_a,
      weights_supplied = weights_supplied,
      all_times = panel_owner$all_times,
      times = panel_owner$times,
      cohorts = panel_owner$cohorts,
      observations = panel_owner$observations,
      explicit_na_rows = panel_owner$explicit_na_rows,
      units = panel_owner$units,
      singletons = panel_owner$singletons,
      unobserved_units = panel_owner$unobserved_units,
      frame_only_units = frame_only_units,
      unit_weights_all = unit_weights_all,
      counts = panel_owner$counts,
      cols = panel_owner$cols,
      postings = panel_owner$postings,
      control = ctrl,
      preflight = preflight,
      accounting = list(
        retained_matrix_entries = as.double(n_a) * as.double(p),
        retained_vector_entries = as.double(n_a) + as.double(nu) + 2 * length(a_cohorts) + 6 * length(absent_from_panel),
        Z_aliased_from_input = is_canonical_a_input && is.double(Z) && is.null(rownames(Z)),
        w_aliased_from_input = weights_supplied && is_canonical_a_input && is.double(weights) && is.null(attributes(weights)),
        panel_arrays_aliased = TRUE
      ),
      work = list(
        validation_chunks_processed = logit_cert_count_output(validation_chunks),
        actual_work_ops = actual_ops
      )
    ),
    class = c("cdid_ipw_frame_owner", "cdid_ipw_context")
  )

  ipw_owner_seal(frame_obj, "frame_owner")
}

construct_common_frame_owner <- construct_ipw_frame_owner
chained_ipw_frame_owner <- construct_ipw_frame_owner

validate_ipw_frame_owner <- function(frame_owner, panel_owner = NULL) {
  ipw_owner_require_deps()
  if (missing(frame_owner) || is.null(frame_owner) ||
      !inherits(frame_owner, "cdid_ipw_frame_owner") ||
      !is.list(frame_owner) || !ipw_owner_seal_ok(frame_owner, "frame_owner")) {
    cdid_ipw_owner_abort(
      "Frame owner is malformed, foreign, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_frame_error", "cdid_ipw_context_mismatch_error")
    )
  }
  pref <- frame_owner$panel_ref
  if (!is.environment(pref) || !environmentIsLocked(pref) ||
      !exists("panel_owner", envir = pref, inherits = FALSE) ||
      !bindingIsLocked("panel_owner", pref)) {
    cdid_ipw_owner_abort(
      "Frame owner lacks a valid sealed `panel_owner` reference.",
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_frame_error")
    )
  }
  p_obj <- get("panel_owner", envir = pref, inherits = FALSE)
  validate_ipw_panel_owner(p_obj)
  if (!is.null(panel_owner)) {
    if (!gmm_sexp_aliases(panel_owner, p_obj)) {
      validate_ipw_panel_owner(panel_owner)
    }
    if (!identical(pref, panel_owner$ref)) {
      cdid_ipw_owner_abort(
        "Frame owner belongs to a different `cdid_ipw_panel_owner` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_owner_frame_error")
      )
    }
  }
  invisible(frame_owner)
}

validate_common_frame_owner <- validate_ipw_frame_owner

# Verify compatibility of two frame owners: labels alone NEVER prove reuse.
verify_ipw_frame_compatibility <- function(frame_a, frame_b, require_same_owner = TRUE) {
  validate_ipw_frame_owner(frame_a); validate_ipw_frame_owner(frame_b)
  if (!gmm_scalar_logical(require_same_owner)) cdid_ipw_owner_abort("`require_same_owner` must be scalar logical.", "cdid_ipw_owner_input_error")
  if (identical(frame_a$ref, frame_b$ref)) return(invisible(TRUE))
  refuse <- function() cdid_ipw_owner_abort("Frame owners differ in actual owner, A-unit set, order, weights or design.", c("cdid_ipw_owner_mismatch_error", "cdid_ipw_owner_frame_error"))
  if (require_same_owner || !identical(frame_a$panel_ref, frame_b$panel_ref) ||
      !identical(frame_a$n_a_units, frame_b$n_a_units) || !identical(frame_a$n_covariates, frame_b$n_covariates) ||
      !identical(frame_a$weights_supplied, frame_b$weights_supplied)) refuse()
  if (!gmm_sexp_aliases(frame_a$a_id_key, frame_b$a_id_key) ||
      !gmm_sexp_aliases(frame_a$a_g, frame_b$a_g)) refuse()
  # Bounds were explicitly supplied by the constructors and remain sealed.
  ca <- frame_a$control; cb <- frame_b$control
  n <- frame_a$n_a_units; p <- frame_a$n_covariates
  step <- min(n, ca$max_chunk_rows, cb$max_chunk_rows)
  ipw_owner_safe_product(8 * (p + 4), n, min(ca$max_work_ops, cb$max_work_ops), "bounded frame comparison work")
  ipw_owner_safe_sum(c(4 * step * p, 16 * step, 128), min(ca$max_matrix_entries, cb$max_matrix_entries), "bounded frame comparison scratch")
  left <- 1
  while (left <= n) {
    ipw_owner_check_interrupt()
    right <- min(n, left + step - 1); idx <- seq.int(left, right)
    for (nm in "w_a") if (!identical(frame_a[[nm]][idx], frame_b[[nm]][idx], num.eq = FALSE)) refuse()
    if (!identical(frame_a$Z_a[idx,,drop=FALSE], frame_b$Z_a[idx,,drop=FALSE], num.eq=FALSE)) refuse()
    left <- right + 1
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Component 3: Compact Outcome-Free Nuisance & Link Views
# ---------------------------------------------------------------------------

ipw_owner_make_pair_descriptor <- function(panel_owner, s, t) {
  if (!gmm_scalar_numeric(s) || !gmm_scalar_numeric(t)) {
    cdid_ipw_owner_abort(
      "Date-pair endpoints `s` and `t` must be finite numeric scalars.",
      c("cdid_ipw_owner_pair_error", "cdid_ipw_pair_error", "cdid_ipw_input_error")
    )
  }
  s_num <- as.double(s)
  t_num <- as.double(t)
  s_num[s_num == 0.0] <- 0.0
  t_num[t_num == 0.0] <- 0.0
  if (!(s_num < t_num)) {
    cdid_ipw_owner_abort(
      "Date-pair descriptor requires strict chronological order `s < t`.",
      c("cdid_ipw_owner_pair_error", "cdid_ipw_pair_error", "cdid_ipw_input_error")
    )
  }
  ds <- ipw_owner_match_sorted_double(s_num, panel_owner$all_times)
  dt <- ipw_owner_match_sorted_double(t_num, panel_owner$all_times)
  gap <- t_num - s_num
  list(
    s = s_num,
    t = t_num,
    s_calendar_index = ds,
    t_calendar_index = dt,
    time_gap = gap,
    time_gap_finite = is.finite(gap)
  )
}

ipw_owner_make_cohort_descriptor <- function(frame_owner, g, K, require_support = TRUE) {
  if (!gmm_scalar_numeric(g)) {
    cdid_ipw_owner_abort(
      "Target cohort `g` must be a single finite numeric value.",
      c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_error", "cdid_ipw_input_error")
    )
  }
  if (missing(K) || is.null(K) || !is.atomic(K) || !is.null(dim(K)) ||
      is.object(K) || !is.numeric(K) || length(K) < 1L ||
      anyNA(K) || any(is.nan(K) | K == -Inf)) {
    cdid_ipw_owner_abort(
      "Declared controls `K` must be an explicit nonempty numeric cohort vector (finite values or `+Inf`).",
      c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_error", "cdid_ipw_input_error")
    )
  }
  g_num <- as.double(g)
  g_num[g_num == 0.0] <- 0.0
  K_num <- as.double(K)
  K_num[K_num == 0.0] <- 0.0
  if (anyDuplicated(K_num) || any(K_num == g_num)) {
    cdid_ipw_owner_abort(
      "Declared controls `K` must be distinct cohort values excluding target `g`.",
      c("cdid_ipw_owner_cohort_error", "cdid_ipw_cohort_error", "cdid_ipw_input_error")
    )
  }
  cohorts <- c(g_num, K_num)
  m <- length(cohorts)
  cohort_labels <- c("g", paste0("k_", seq_along(K_num)))

  ci <- ipw_owner_match_sorted_double(cohorts, frame_owner$a_cohorts)
  group_a_count <- integer(m)
  group_a_mass <- numeric(m)
  present <- which(!is.na(ci))
  if (length(present) > 0L) {
    group_a_count[present] <- frame_owner$cohort_a_count[ci[present]]
    group_a_mass[present] <- frame_owner$cohort_a_weight[ci[present]]
  }
  if (isTRUE(require_support)) {
    if (anyNA(ci) || any(group_a_count < 1L) || any(!is.finite(group_a_mass) | group_a_mass <= 0)) {
      cdid_ipw_owner_abort(
        "All declared cohorts `{g, K}` require positive A-unit count and weight mass; unsupported cohorts are not silently discarded.",
        c("cdid_ipw_owner_cohort_support_error", "cdid_ipw_cohort_support_error", "cdid_ipw_support_error")
      )
    }
  }
  q <- group_a_mass / as.double(frame_owner$n_a_units)
  if (isTRUE(require_support) && any(!is.finite(q) | q <= 0)) {
    cdid_ipw_owner_abort(
      "Invalid A-unit cohort share normalization `q_j`.",
      c("cdid_ipw_owner_cohort_support_error", "cdid_ipw_normalization_error")
    )
  }

  # Store only declared supported cohorts, never a global-G (potentially N)
  # lookup vector per view. Sorted owner indices map to the original declared order.
  lookup_order <- order(ci, na.last = NA, method = "radix")
  owner_cohort_index <- ci[lookup_order]
  view_cohort_index <- lookup_order

  names(group_a_count) <- names(group_a_mass) <- names(q) <- cohort_labels
  list(
    g = g_num,
    K = K_num,
    cohorts = cohorts,
    cohort_labels = cohort_labels,
    m_cohorts = m,
    a_cohort_map = ci,
    owner_cohort_index = owner_cohort_index,
    view_cohort_index = view_cohort_index,
    group_a_count = group_a_count,
    group_a_mass = group_a_mass,
    q = q
  )
}

ipw_owner_view_cohort_index <- function(coh_desc, owner_index) {
  pos <- ipw_owner_match_sorted_double(owner_index, coh_desc$owner_cohort_index)
  coh_desc$view_cohort_index[pos]
}

# ---------------------------------------------------------------------------
# Exact complete-history posting facts (CP08 balanced performance)
# ---------------------------------------------------------------------------
# Calendar index `d` is complete when its posting slice is exactly the A-units
# 1..n_A, i.e. every pair-eligible unit is observed at d. The slice is read in
# chunks of at most `step` postings (the scratch of one streamed chunk): each
# chunk must be strictly increasing integers continuing from the previous one,
# so all chunks together are exactly 1..n_A. One O(n_A) read of the postings
# every streamed pair query reads; no outcome is read and nothing is retained.
ipw_owner_posting_date_complete <- function(postings, n_a, d, step = NULL) {
  if (!is.numeric(d) || length(d) != 1L || is.na(d) || d != floor(d) ||
      !is.numeric(n_a) || length(n_a) != 1L || is.na(n_a) || n_a < 1 ||
      d < 1 || d > length(postings$date_count)) return(FALSE)
  cnt <- postings$date_count[d]
  s <- postings$date_start[d]
  e <- postings$date_end[d]
  pa <- postings$post_a_index
  if (!is.integer(pa) || is.na(cnt) || cnt != n_a || is.na(s) || is.na(e) ||
      s < 1 || e > length(pa) || e - s + 1 != n_a) return(FALSE)
  step <- if (is.numeric(step) && length(step) == 1L && !is.na(step) && step >= 1) floor(step) else n_a
  pos <- s
  first <- 1
  while (pos <= e) {
    hi <- min(e, pos + step - 1)
    sl <- pa[seq.int(pos, hi)]
    if (anyNA(sl) || sl[1L] != first || sl[length(sl)] != first + (hi - pos) ||
        is.unsorted(sl, strictly = TRUE)) return(FALSE)
    first <- first + (hi - pos) + 1
    pos <- hi + 1
  }
  TRUE
}

# Per-date completeness, memoized in `memo` (an environment bound by its owner
# to this one sealed postings object) when supplied.
ipw_owner_date_complete <- function(owner, d, step = NULL, memo = NULL) {
  if (!is.environment(memo) || !is.numeric(d) || length(d) != 1L || is.na(d)) {
    return(ipw_owner_posting_date_complete(owner$postings, owner$n_a_units, d, step))
  }
  key <- paste0("complete:", d)
  hit <- get0(key, envir = memo, inherits = FALSE)
  if (is.logical(hit)) return(hit)
  val <- ipw_owner_posting_date_complete(owner$postings, owner$n_a_units, d, step)
  assign(key, val, envir = memo)
  val
}

# Both dates of a pair descriptor are complete on this owner's postings.
ipw_owner_pair_complete <- function(owner, pair_desc, step = NULL, memo = NULL) {
  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  !is.null(ds) && !is.null(dt) &&
    ipw_owner_date_complete(owner, ds, step, memo) &&
    ipw_owner_date_complete(owner, dt, step, memo)
}

# Equality of two pair-observation masks on one frame decided without streaming:
# identical pair descriptors make both streamed reads the same query, and four
# complete dates make both masks every A-unit, chunk by chunk. FALSE means
# "not decided"; the caller then streams exactly as before.
ipw_owner_pair_masks_equal_fact <- function(owner, p1, p2, step = NULL, memo = NULL) {
  if (identical(p1$s, p2$s) && identical(p1$t, p2$t) &&
      identical(p1$s_calendar_index, p2$s_calendar_index) &&
      identical(p1$t_calendar_index, p2$t_calendar_index)) return(TRUE)
  ipw_owner_pair_complete(owner, p1, step, memo) && ipw_owner_pair_complete(owner, p2, step, memo)
}

# Pair hits of [a_start, a_end] when both dates are complete: the positions
# `ipw_owner_query_pair_mask` reaches by binary search, read directly
# (slice position date_start + a - 1 holds A-unit a at a complete date).
ipw_owner_query_complete_pair_mask <- function(context, ds, dt, a_start, a_end) {
  post <- context$postings
  a <- seq.int(a_start, a_end)
  ps <- post$date_start[ds] - 1L + a
  pt <- post$date_start[dt] - 1L + a
  list(
    a_index = post$post_a_index[ps],
    unit_index = post$post_unit_index[ps],
    obs_s = post$post_obs_index[ps],
    obs_t = post$post_obs_index[pt]
  )
}

# Invocation-local construction facts. Only the trusted registry loop
# (`build_pipeline_nuisance_and_link_views`) creates one, from its validated
# frame owner, and drops it before returning (CP08 stage 5: so does the CP04
# registry's view helper, see ipw_owner_link_support_context); standalone
# constructor calls never receive one and scan as before. It records the exact
# support scan of every nuisance view constructed in this call (by nuisance
# ID, bound to the view's seal reference), memoizes complete-pair mask
# fingerprints, whose value is the same for every complete pair under one
# cohort index and chunk step, and memoizes per-date completeness of the
# frame's postings.
ipw_owner_construction_context <- function(frame_owner) {
  ctx <- new.env(parent = emptyenv())
  ctx$frame_ref <- frame_owner$ref
  ctx$panel_ref <- frame_owner$panel_ref
  ctx$nuisance <- new.env(parent = emptyenv())
  ctx$fingerprints <- new.env(parent = emptyenv())
  ctx$dates <- new.env(parent = emptyenv())
  ctx$support_reuses <- 0
  ctx$fingerprint_reuses <- 0
  ctx
}

ipw_owner_construction_context_for <- function(ctx, frame_owner) {
  is.environment(ctx) && identical(ctx$frame_ref, frame_owner$ref) &&
    identical(ctx$panel_ref, frame_owner$panel_ref) &&
    is.environment(ctx$nuisance) && is.environment(ctx$fingerprints) && is.environment(ctx$dates)
}

# Records in the context that this exact sealed frame owner was validated in
# the current invocation. The validation is performed here, so the recorded
# reference is never a bare flag: it exists only after a completed
# `validate_ipw_frame_owner` of the snapshot the context is bound to.
ipw_owner_context_check_frame <- function(ctx, frame_owner) {
  validate_ipw_frame_owner(frame_owner)
  if (is.environment(ctx) && ipw_owner_construction_context_for(ctx, frame_owner)) {
    ctx$checked_frame_ref <- frame_owner$ref
  }
  invisible(frame_owner)
}

# A context marker only selects a candidate. Re-establish the positive predicates
# of the frame and panel validators on their actual aliased sealed snapshots;
# a caller-populated mutable marker cannot bypass either boundary. An undecided
# or malformed candidate returns FALSE so the original validator raises its
# original typed refusal in the original order.
ipw_owner_context_frame_checked <- function(ctx, frame_owner) {
  if (!(is.environment(ctx) && is.list(frame_owner) &&
    inherits(frame_owner, "cdid_ipw_frame_owner") && is.environment(ctx$checked_frame_ref) &&
    identical(ctx$checked_frame_ref, frame_owner$ref) && ipw_owner_construction_context_for(ctx, frame_owner) &&
    ipw_owner_seal_ok(frame_owner, "frame_owner"))) return(FALSE)
  pref <- frame_owner$panel_ref
  if (!is.environment(pref) || !environmentIsLocked(pref) ||
      !exists("panel_owner", envir = pref, inherits = FALSE) || !bindingIsLocked("panel_owner", pref)) return(FALSE)
  panel <- get("panel_owner", envir = pref, inherits = FALSE)
  if (!inherits(panel, "cdid_ipw_panel_owner") || !is.list(panel) ||
      !ipw_owner_seal_ok(panel, "panel_owner")) return(FALSE)
  bm <- panel$id_byte_metadata
  is.list(bm) && identical(names(bm), c("total", "max")) &&
    gmm_scalar_numeric(bm$total, nonneg = TRUE) && gmm_scalar_numeric(bm$max, nonneg = TRUE) &&
    bm$total == floor(bm$total) && bm$max == floor(bm$max) && bm$total >= bm$max && bm$total <= 2^53 - 1
}

# One checked nuisance/support template per view this invocation constructed:
# the record written by `construct_ipw_nuisance_view` after sealing (bound to
# the exact sealed reference, the context's validated frame and its panel).
# TRUE only for that same sealed snapshot with its seal intact; a stale copy,
# a foreign view with a repeated ID, or another frame never matches.
ipw_owner_context_nuisance_checked <- function(ctx, nuisance_view, frame_owner = NULL) {
  if (!is.environment(ctx) || !is.list(nuisance_view) || !inherits(nuisance_view, "cdid_ipw_nuisance_view") ||
      !is.environment(ctx$checked_frame_ref) || !is.environment(ctx$nuisance)) return(FALSE)
  if (!identical(ctx$checked_frame_ref, nuisance_view$frame_ref) ||
      !identical(ctx$frame_ref, nuisance_view$frame_ref) ||
      !identical(ctx$panel_ref, nuisance_view$panel_ref)) return(FALSE)
  if (is.null(frame_owner)) {
    fref <- nuisance_view$frame_ref
    if (!is.environment(fref) || !environmentIsLocked(fref) ||
        !exists("frame_owner", envir = fref, inherits = FALSE) || !bindingIsLocked("frame_owner", fref)) return(FALSE)
    frame_owner <- get("frame_owner", envir = fref, inherits = FALSE)
  }
  if (!ipw_owner_context_frame_checked(ctx, frame_owner) ||
      !identical(nuisance_view$panel_ref, frame_owner$panel_ref)) return(FALSE)
  id <- nuisance_view$nuisance_id
  if (!gmm_scalar_character(id)) return(FALSE)
  rec <- get0(paste0("n:", id), envir = ctx$nuisance, inherits = FALSE)
  is.list(rec) && is.environment(rec$ref) && identical(rec$ref, nuisance_view$ref) &&
    ipw_owner_seal_ok(nuisance_view, "nuisance_view")
}

ipw_owner_context_fingerprint_key <- function(owner_cohort_index, step) {
  # Bounded bucket key. A collision only replaces a memo entry: the caller
  # still compares the entire cohort vector exactly before any reuse.
  first <- if (length(owner_cohort_index)) owner_cohort_index[1L] else 0
  paste0("fp:", step, ":", first, ":", length(owner_cohort_index))
}

# Records the exact support scan of a nuisance view this invocation constructed.
ipw_owner_context_record_nuisance <- function(ctx, sealed_view, sup_info, step, require_support) {
  id <- sealed_view$nuisance_id
  if (!gmm_scalar_character(id)) return(invisible(FALSE))
  key <- paste0("n:", id)
  if (exists(key, envir = ctx$nuisance, inherits = FALSE)) {
    # A repeated ID is ambiguous: neither view's scan is reused.
    assign(key, NULL, envir = ctx$nuisance)
    return(invisible(FALSE))
  }
  assign(key, list(ref = sealed_view$ref, cohorts = sealed_view$cohorts, pair = sealed_view$defining_pair,
    step = step, require_support = require_support, support = sup_info), envir = ctx$nuisance)
  invisible(TRUE)
}

# A link bound to a nuisance view constructed in this invocation has that view's
# support summary whenever their pair masks are equal by an exact posting fact:
# same frame, same cohort descriptor, chunk step and support rule, hence the
# same chunk-by-chunk hits, counts and stable mass sums. NULL means scan.
ipw_owner_context_link_support <- function(ctx, frame_owner, nuisance_view, pair_desc, coh_desc,
                                           step, require_support) {
  id <- nuisance_view$nuisance_id
  if (!gmm_scalar_character(id)) return(NULL)
  rec <- get0(paste0("n:", id), envir = ctx$nuisance, inherits = FALSE)
  if (is.null(rec) || !identical(rec$ref, nuisance_view$ref) || !identical(rec$cohorts, coh_desc) ||
      !identical(rec$step, step) || !identical(rec$require_support, require_support) ||
      !ipw_owner_pair_masks_equal_fact(frame_owner, rec$pair, pair_desc, step, ctx$dates)) return(NULL)
  ctx$support_reuses <- ctx$support_reuses + 1
  rec$support
}

# CP08 stage 5: support records of nuisance-free link views. Only one CP04
# registry invocation (`cp04_public_link_views`) creates this context variant,
# after admitting its whole reservation, and drops it before fitting; the IPW
# registry context and standalone calls are unchanged. A record holds only the
# seal reference of a view this invocation scanned and sealed, its chunk step,
# support rule and two scalar scan counts: the support values stay in that
# immutable sealed view and are borrowed, never copied. Per view (at most one
# record each) 128 entries cover the record wrapper and scalars, its bounded
# key and at most two date facts/keys, as in the IPW registry; one further
# block covers the context itself. Admission also covers the largest declared
# constructor while these records stay live. Not admitted: NULL, and every view scans.
ipw_owner_link_support_context <- function(frame_owner, n_views, retained_entries, control,
                                           max_cohorts = length(frame_owner$cohorts)) {
  if (!is.list(frame_owner) || !is.environment(frame_owner$ref) || !is.environment(frame_owner$panel_ref) ||
      !is.list(control) || !gmm_scalar_numeric(control$max_retained_entries) ||
      !gmm_scalar_numeric(control$max_matrix_entries) || !gmm_scalar_numeric(n_views) || n_views < 2 ||
      !gmm_scalar_numeric(retained_entries, nonneg = TRUE)) return(NULL)
  entries <- 128 * (n_views + 1)
  # Every CP04 link has no optional candidates. The exact retained/peak
  # formulas from preflight_ipw_link_view are monotone in m; decide using the
  # largest declared class before holding any context. A per-view bypass is
  # insufficient: records of an earlier class would still be live.
  if (!gmm_scalar_numeric(max_cohorts, nonneg = TRUE) ||
      !gmm_scalar_numeric(frame_owner$n_covariates) || !gmm_scalar_numeric(frame_owner$n_a_units) ||
      !gmm_scalar_numeric(control$max_chunk_rows)) return(NULL)
  m <- max(2, max_cohorts)
  constructor_retained <- 20 * m + frame_owner$n_covariates + 256
  constructor_peak <- constructor_retained + 64 * min(frame_owner$n_a_units, control$max_chunk_rows) + 6 * m + 256
  total_retained <- entries + max(retained_entries, constructor_retained)
  total_peak <- entries + max(retained_entries, constructor_peak)
  if (!is.finite(total_retained) || !is.finite(total_peak) || max(total_retained, total_peak) > 2^53 - 1 ||
      total_retained > control$max_retained_entries || total_peak > control$max_matrix_entries) return(NULL)
  ctx <- ipw_owner_construction_context(frame_owner)
  ctx$links <- new.env(parent = emptyenv())
  ctx$link_capacity <- as.double(n_views)
  ctx$link_entries <- entries
  ctx$link_retained <- as.double(retained_entries)
  ctx$link_records <- 0
  ctx$link_scans <- 0
  ctx$link_reuses <- 0
  ctx
}

# A nuisance-free view uses a link-support context only when the context's
# whole reservation also composes with this view's own completed constructor
# preflight (the IPW registry composition: records plus the larger of the
# pipeline envelope and the constructor), decided before any context lookup
# or record. Otherwise the view scans exactly as without a context.
ipw_owner_context_link_admits <- function(ctx, frame_owner, preflight, ctrl) {
  isTRUE(ipw_owner_construction_context_for(ctx, frame_owner) && is.environment(ctx$links) &&
    gmm_scalar_numeric(ctx$link_entries) && gmm_scalar_numeric(ctx$link_retained) &&
    gmm_scalar_numeric(ctx$link_capacity) && gmm_scalar_numeric(ctx$link_records) &&
    ctx$link_entries + max(ctx$link_retained, preflight$total_retained_entries) <= ctrl$max_retained_entries &&
    ctx$link_entries + max(ctx$link_retained, preflight$concurrent_entries) <= ctrl$max_matrix_entries)
}

# Bounded scalar bucket key of a record class; never evidence of equality.
# Complete-pair records share one key per class, other records are keyed by
# their pair's calendar indices. A colliding class finds a record whose exact
# identity differs, so it scans and adds no record.
ipw_owner_context_link_key <- function(coh_desc, step, require_support, pair_desc = NULL) {
  map <- coh_desc$a_cohort_map
  key <- paste0("ls:", step, ":", require_support, ":", length(map), ":", map[1L], ":", map[2L], ":",
    map[length(map)])
  if (is.null(pair_desc)) paste0(key, ":complete") else
    paste0(key, ":pair:", pair_desc$s_calendar_index, ":", pair_desc$t_calendar_index)
}

# Both pair dates post exactly n_A observations: a count-only precondition of
# completeness, checked before any posting slice is read.
ipw_owner_pair_counts_full <- function(frame_owner, pair_desc) {
  d <- c(pair_desc$s_calendar_index, pair_desc$t_calendar_index)
  length(d) == 2L && !anyNA(d) && isTRUE(all(frame_owner$postings$date_count[d] == frame_owner$n_a_units))
}

# Support of a nuisance-free link view sealed earlier in this invocation, or
# NULL (scan). Bounded keys (this exact pair; the class's complete pairs) only
# find a candidate. It is accepted when the sealed view behind the record is
# intact, carries no nuisance view, belongs to these frame/panel snapshots,
# matches the record's scan counts and has exactly this cohort descriptor,
# chunk step and support rule, and when the existing posting fact proves the
# two pair masks equal (identical pair, or four complete dates). Equal masks
# under these give the same chunk-by-chunk hits, counts and stable mass sums,
# so the borrowed values are the original scan's own values, names and order.
ipw_owner_context_link_lookup <- function(ctx, frame_owner, pair_desc, coh_desc, step, require_support) {
  counts_full <- ipw_owner_pair_counts_full(frame_owner, pair_desc)
  for (key in c(ipw_owner_context_link_key(coh_desc, step, require_support, pair_desc),
                ipw_owner_context_link_key(coh_desc, step, require_support))) {
    rec <- get0(key, envir = ctx$links, inherits = FALSE)
    if (!is.list(rec) || !is.environment(rec$ref) || !identical(rec$step, step) ||
        !identical(rec$require_support, require_support) || !gmm_scalar_numeric(rec$chunks) ||
        !gmm_scalar_numeric(rec$ops)) next
    sv <- get0("link_view", envir = rec$ref, inherits = FALSE)
    if (!inherits(sv, "cdid_ipw_link_view") || !identical(sv$ref, rec$ref) || !ipw_owner_seal_ok(sv, "link_view") ||
        !is.null(sv$nuisance_ref) || !identical(sv$frame_ref, frame_owner$ref) ||
        !identical(sv$panel_ref, frame_owner$panel_ref) || !identical(sv$preflight$effective_chunk_rows, step) ||
        !identical(sv$work$actual_work_ops, rec$ops) ||
        !identical(sv$work$support_chunks_processed, logit_cert_count_output(rec$chunks)) ||
        !identical(sv$cohorts, coh_desc)) next
    if (!(identical(sv$pair, pair_desc) || counts_full) || !ipw_owner_pair_masks_equal_fact(frame_owner, sv$pair, pair_desc, step, ctx$dates)) next
    ctx$link_reuses <- ctx$link_reuses + 1
    s <- sv$support
    return(list(
      n_observed_pair_all = s$n_observed_pair_all,
      n_risk_a_units = s$n_risk_a_units,
      n_outside_risk_a_units = s$n_outside_risk_a_units,
      n_observed_risk_units = s$n_observed_risk_units,
      n_nonobserved_risk_units = s$n_nonobserved_risk_units,
      by_cohort_pair_count = s$by_cohort_pair_count,
      by_cohort_pair_mass = s$by_cohort_pair_mass,
      rates = s$rates,
      chunks = rec$chunks,
      ops = rec$ops
    ))
  }
  NULL
}

# Records a nuisance-free link view this invocation has just scanned and
# sealed: one record per bounded key (the first wins, so a collision adds
# none) and at most the admitted number of views. Failed scans and borrowed
# supports never reach this point.
ipw_owner_context_record_link <- function(ctx, frame_owner, sealed_link, sup_info, step, require_support) {
  if (ctx$link_records >= ctx$link_capacity) return(invisible(FALSE))
  pair <- sealed_link$pair
  complete <- ipw_owner_pair_counts_full(frame_owner, pair) &&
    ipw_owner_pair_complete(frame_owner, pair, step, ctx$dates)
  key <- ipw_owner_context_link_key(sealed_link$cohorts, step, require_support, if (complete) NULL else pair)
  if (exists(key, envir = ctx$links, inherits = FALSE)) return(invisible(FALSE))
  assign(key, list(ref = sealed_link$ref, step = step, require_support = require_support,
    chunks = sup_info$chunks, ops = sup_info$ops), envir = ctx$links)
  ctx$link_records <- ctx$link_records + 1
  invisible(TRUE)
}

# Outcome-free chunked pair-support summary using only O(B + m) scratch
ipw_owner_chunked_pair_support <- function(frame_owner, pair_desc, coh_desc, step, require_support = TRUE) {
  n_a <- frame_owner$n_a_units
  m <- coh_desc$m_cohorts
  pair_count <- integer(m)
  pair_mass <- numeric(m)
  n_obs_pair_all <- 0L
  chunks <- 0
  ops <- 0

  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  a_cohort_idx <- frame_owner$a_cohort_index
  w_a <- frame_owner$w_a

  if (!is.na(ds) && !is.na(dt) &&
      frame_owner$postings$date_count[ds] > 0L &&
      frame_owner$postings$date_count[dt] > 0L) {
    a_start <- 1L
    while (a_start <= n_a) {
      ipw_owner_check_interrupt()
      a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
      hit <- ipw_owner_query_pair_mask(
        frame_owner, pair_desc$s, pair_desc$t,
        a_start = a_start, a_end = a_end, ds = ds, dt = dt
      )
      nh <- length(hit$a_index)
      if (nh > 0L) {
        n_obs_pair_all <- n_obs_pair_all + nh
        gi <- ipw_owner_view_cohort_index(coh_desc, a_cohort_idx[hit$a_index])
        in_r <- which(!is.na(gi))
        if (length(in_r) > 0L) {
          gi_r <- gi[in_r]
          w_r <- w_a[hit$a_index[in_r]]
          # Group only this bounded chunk. Stable order preserves each old sum;
          # avoid scanning the whole chunk for each declared cohort or creating
          # an m-vector accumulator on every batch.
          group_order <- order(gi_r, method = "radix")
          runs <- rle(gi_r[group_order])
          group_end <- cumsum(runs$lengths)
          group_start <- group_end - runs$lengths + 1L
          for (z in seq_along(runs$values)) {
            j <- runs$values[z]
            pair_count[j] <- pair_count[j] + runs$lengths[z]
            pair_mass[j] <- pair_mass[j] + sum(w_r[group_order[seq.int(group_start[z], group_end[z])]])
          }
        }
      }
      chunks <- chunks + 1
      ops <- ops + 4 * (a_end - a_start + 1L)
      a_start <- as.double(a_end) + 1
    }
  } else {
    chunks <- ceiling(n_a / step)
    ops <- n_a
  }

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
    # Constructors validate unique (ID, time) rows and positive weights. Posting
    # intersections in disjoint ordered A-unit chunks are unique subsets of the
    # same roster/cohort map used for group_a_count and group_a_mass. Equal
    # positive counts thus prove complete support here, not for arbitrary arrays.
    # Its empirical rate is exactly one; retain the independent pair_mass sum.
    # require_support = FALSE keeps the original unchecked division semantics.
    rates[pair_count == coh_desc$group_a_count] <- 1
    ops <- ops + 4 * m # comparison, logical indexing, rate copy, replacement
  }
  if (isTRUE(require_support) && any(!is.finite(rates) | rates <= 0 | rates > 1)) {
    cdid_ipw_owner_abort(
      "Invalid marginal pair observation rates `r_j`.",
      c("cdid_ipw_owner_pair_support_error", "cdid_ipw_rate_error")
    )
  }
  names(pair_count) <- names(pair_mass) <- names(rates) <- coh_desc$cohort_labels

  list(
    n_observed_pair_all = n_obs_pair_all,
    n_risk_a_units = sum(coh_desc$group_a_count),
    n_outside_risk_a_units = n_a - sum(coh_desc$group_a_count),
    n_observed_risk_units = n_obs_risk,
    n_nonobserved_risk_units = sum(coh_desc$group_a_count) - n_obs_risk,
    by_cohort_pair_count = pair_count,
    by_cohort_pair_mass = pair_mass,
    rates = rates,
    chunks = chunks,
    ops = ops
  )
}

ipw_owner_refuse_dense_masks <- function(
  pair_observed = NULL,
  risk_mask = NULL,
  observed_mask = NULL,
  cohort_index = NULL,
  caller = "view"
) {
  if (!is.null(pair_observed) || !is.null(risk_mask) ||
      !is.null(observed_mask) || !is.null(cohort_index)) {
    cdid_ipw_owner_abort(
      paste0(
        "Arbitrary supplied dense masks (`pair_observed`, `risk_mask`, `observed_mask`, `cohort_index`) are outside the compact posting-derived representation in `",
        caller,
        "` (use `verify_ipw_view_dense_masks` only when explicitly bounding and disclosing external dense mask verification)."
      ),
      c("cdid_ipw_owner_dense_mask_error", "cdid_ipw_enclosure_mask_error")
    )
  }
  invisible(TRUE)
}

construct_ipw_nuisance_view <- function(
  frame_owner,
  g,
  K,
  s,
  t,
  control,
  beta0 = NULL,
  F = NULL,
  r = NULL,
  nuisance_id = NULL,
  require_support = TRUE,
  pair_observed = NULL,
  risk_mask = NULL,
  observed_mask = NULL,
  cohort_index = NULL,
  construction_context = NULL
) {
  ipw_owner_require_deps()
  if (missing(frame_owner) || missing(g) || missing(K) || missing(s) || missing(t) || missing(control)) {
    cdid_ipw_owner_abort(
      "All required arguments (`frame_owner`, `g`, `K`, `s`, `t`, `control`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }

  # The exact sealed frame owner a bound construction context validated in
  # this invocation is not re-validated per view; any other object is.
  if (!ipw_owner_context_frame_checked(construction_context, frame_owner)) validate_ipw_frame_owner(frame_owner)
  ctrl <- validate_ipw_owner_control(control)
  ipw_owner_refuse_dense_masks(
    pair_observed = pair_observed,
    risk_mask = risk_mask,
    observed_mask = observed_mask,
    cohort_index = cohort_index,
    caller = "construct_ipw_nuisance_view"
  )
  if (!gmm_scalar_logical(require_support)) {
    cdid_ipw_owner_abort("`require_support` must be a scalar logical.", "cdid_ipw_owner_input_error")
  }
  if (!is.null(nuisance_id) && !gmm_scalar_character(nuisance_id)) {
    cdid_ipw_owner_abort("`nuisance_id` must be NULL or a nonempty character scalar.", "cdid_ipw_owner_input_error")
  }
  for (item in list(
    list(v = g, nm = "g"), list(v = K, nm = "K"),
    list(v = s, nm = "s"), list(v = t, nm = "t"),
    list(v = beta0, nm = "beta0"), list(v = F, nm = "F"), list(v = r, nm = "r")
  )) {
    ipw_owner_refuse_decorative_attrs(item$v, item$nm)
  }

  panel_owner <- get("panel_owner", envir = frame_owner$panel_ref, inherits = FALSE)
  p <- frame_owner$n_covariates
  m <- 1L + (if (is.atomic(K) && is.null(dim(K))) length(K) else 0L)
  has_cand <- !is.null(beta0) || !is.null(F) || !is.null(r)

  # Preflight BEFORE scanning postings or candidate matrices
  preflight <- preflight_ipw_nuisance_view(
    n_a_units = frame_owner$n_a_units,
    p_covariates = p,
    m_cohorts = max(2L, m),
    control = ctrl,
    has_candidates = has_cand
  )
  step <- preflight$effective_chunk_rows

  pair_desc <- ipw_owner_make_pair_descriptor(panel_owner, s, t)
  coh_desc <- ipw_owner_make_cohort_descriptor(frame_owner, g, K, require_support = require_support)

  beta0_dbl <- NULL
  F_dbl <- NULL
  r_dbl <- NULL
  cand_ops <- 0
  if (!is.null(beta0)) {
    if (!is.atomic(beta0) || !is.null(dim(beta0)) || is.object(beta0) ||
        !is.numeric(beta0) || length(beta0) != p || anyNA(beta0) || any(!is.finite(beta0))) {
      cdid_ipw_owner_abort(
        paste0("`beta0` must be a finite numeric vector of length `n_covariates` (", p, ")."),
        c("cdid_ipw_owner_candidate_error", "cdid_ipw_enclosure_beta_error")
      )
    }
    beta0_dbl <- as.double(beta0)
    cand_ops <- cand_ops + p
  }
  if (!is.null(F)) {
    F_dbl <- validate_logit_certificate_metric_factor(F, p)
    cand_ops <- cand_ops + p * p
  }
  if (!is.null(r)) {
    if (!gmm_scalar_numeric(r, positive = TRUE)) {
      cdid_ipw_owner_abort(
        "`r` must be a finite positive numeric scalar radius.",
        c("cdid_ipw_owner_candidate_error", "cdid_ipw_enclosure_radius_error")
      )
    }
    r_dbl <- as.double(r)
    cand_ops <- cand_ops + 1
  }

  sup_info <- ipw_owner_chunked_pair_support(
    frame_owner = frame_owner,
    pair_desc = pair_desc,
    coh_desc = coh_desc,
    step = step,
    require_support = require_support
  )
  actual_ops <- sup_info$ops + cand_ops
  stopifnot(actual_ops <= preflight$work_ops_bound)

  nu_view <- structure(
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
      candidates = if (has_cand) list(beta0 = beta0_dbl, F = F_dbl, r = r_dbl) else NULL,
      control = ctrl,
      preflight = preflight,
      accounting = list(
        full_n_masks_retained = FALSE,
        copied_Z_retained = FALSE,
        endpoints_read = 0L
      ),
      work = list(
        support_chunks_processed = logit_cert_count_output(sup_info$chunks),
        actual_work_ops = actual_ops
      )
    ),
    class = "cdid_ipw_nuisance_view"
  )

  sealed_view <- ipw_owner_seal(nu_view, "nuisance_view")
  if (!is.null(construction_context) && ipw_owner_construction_context_for(construction_context, frame_owner)) {
    ipw_owner_context_record_nuisance(construction_context, sealed_view, sup_info, step, require_support)
  }
  sealed_view
}

declare_ipw_nuisance_view <- construct_ipw_nuisance_view

validate_ipw_nuisance_view <- function(nuisance_view, frame_owner = NULL) {
  ipw_owner_require_deps()
  if (missing(nuisance_view) || is.null(nuisance_view) ||
      !inherits(nuisance_view, "cdid_ipw_nuisance_view") ||
      !is.list(nuisance_view) || !ipw_owner_seal_ok(nuisance_view, "nuisance_view")) {
    cdid_ipw_owner_abort(
      "Nuisance view is malformed, foreign, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_view_error")
    )
  }
  fref <- nuisance_view$frame_ref
  if (!is.environment(fref) || !environmentIsLocked(fref) ||
      !exists("frame_owner", envir = fref, inherits = FALSE) ||
      !bindingIsLocked("frame_owner", fref)) {
    cdid_ipw_owner_abort(
      "Nuisance view lacks a valid sealed `frame_owner` reference.",
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_view_error")
    )
  }
  f_obj <- get("frame_owner", envir = fref, inherits = FALSE)
  validate_ipw_frame_owner(f_obj)
  if (!identical(nuisance_view$panel_ref, f_obj$panel_ref)) {
    cdid_ipw_owner_abort(
      "Nuisance view `panel_ref` disagrees with its `frame_owner`.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_owner_view_error")
    )
  }
  if (!is.null(frame_owner)) {
    if (!gmm_sexp_aliases(frame_owner, f_obj)) {
      validate_ipw_frame_owner(frame_owner)
    }
    if (!identical(fref, frame_owner$ref)) {
      cdid_ipw_owner_abort(
        "Nuisance view belongs to a different `cdid_ipw_frame_owner` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_owner_view_error", "cdid_ipw_context_mismatch_error")
      )
    }
  }
  invisible(nuisance_view)
}

# Exact reason text of the streamed mask-mismatch refusal below. A caller that
# has already proved two masks unequal (different fingerprints, see below)
# records this same text, so retained split diagnostics keep their meaning.
ipw_owner_mask_mismatch_message <- function() {
  "Link risk-set pair-observation mask differs from the nuisance view mask on chunk; reuse is refused."
}

# Streamed chunk-by-chunk verification that two views share identical owner,
# identical ordered cohorts, and identical risk-set pair-observation mask O = R intersect S.
# The exact guards of `verify_ipw_view_mask_compatibility` that precede its
# streamed comparison, in their original order: frame/panel owner identity,
# the streamed-work bound, identical ordered cohorts, then the posting fact.
# `decided = TRUE` carries the result the stream would return for masks equal
# by an exact posting fact (identical pair, or four complete dates); otherwise
# the caller streams with the returned owner, step and pair descriptors
# exactly as before. Callers validate `nuisance_view`, the link's seal (or its
# freshness within one constructor invocation) and `ctrl` before calling.
ipw_view_mask_compatibility_fact <- function(nuisance_view, link_view, ctrl, construction_context = NULL) {
  if (!identical(nuisance_view$frame_ref, link_view$frame_ref) ||
      !identical(nuisance_view$panel_ref, link_view$panel_ref)) {
    cdid_ipw_owner_abort(
      "Views belong to different frame or panel owners; reuse is refused.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  n_a <- frame_owner$n_a_units
  step <- min(n_a, ctrl$max_chunk_rows)
  if (64 * step + 32 * nuisance_view$cohorts$m_cohorts + 256 > ctrl$max_matrix_entries ||
      128 * n_a * (1 + ceiling(log2(nuisance_view$cohorts$m_cohorts + 1))) +
      128 * ceiling(n_a / step) * (1 + ceiling(log2(n_a + 1))) > ctrl$max_work_ops) {
    cdid_ipw_owner_abort(
      "Streamed mask comparison exceeds explicit `control` bounds.",
      "cdid_ipw_owner_bound_error"
    )
  }

  if (!identical(nuisance_view$cohorts$cohorts, link_view$cohorts$cohorts) ||
      !identical(nuisance_view$cohorts$cohort_labels, link_view$cohorts$cohort_labels)) {
    cdid_ipw_owner_abort(
      "Target cohort `g` and exact ordered controls `K` must be identical between nuisance and link views.",
      c("cdid_ipw_owner_cohort_error", "cdid_ipw_shared_incompatible_error")
    )
  }

  p1 <- nuisance_view$defining_pair
  p2 <- if (!is.null(link_view$pair)) link_view$pair else link_view$defining_pair

  # Both streamed reads of the comparison filter the hits of each chunk by the
  # same cohort descriptor, so when the two pair masks are equal by an exact
  # posting fact (identical pair, or four complete dates) every chunk compares
  # equal and the stream would return exactly this result. Any other case streams.
  date_memo <- if (!is.null(construction_context) &&
    ipw_owner_construction_context_for(construction_context, frame_owner)) construction_context$dates else NULL
  if (ipw_owner_pair_masks_equal_fact(frame_owner, p1, p2, step, date_memo)) {
    return(list(
      decided = TRUE,
      result = list(
        compatible = TRUE,
        n_a_units = n_a,
        chunks_compared = as.integer(ceiling(n_a / step)),
        endpoints_read = 0L
      ),
      frame_owner = frame_owner, n_a = n_a, step = step, p1 = p1, p2 = p2
    ))
  }
  list(decided = FALSE, result = NULL, frame_owner = frame_owner, n_a = n_a, step = step, p1 = p1, p2 = p2)
}

verify_ipw_view_mask_compatibility <- function(nuisance_view, link_view, control, construction_context = NULL) {
  validate_ipw_nuisance_view(nuisance_view)
  if (!inherits(link_view, "cdid_ipw_link_view") && !inherits(link_view, "cdid_ipw_nuisance_view")) {
    cdid_ipw_owner_abort(
      "Expected a sealed `cdid_ipw_link_view` or `cdid_ipw_nuisance_view`.",
      "cdid_ipw_owner_view_error"
    )
  }
  if (inherits(link_view, "cdid_ipw_link_view")) {
    if (!ipw_owner_seal_ok(link_view, "link_view")) {
      cdid_ipw_owner_abort("Link view differs from its sealed snapshot.", "cdid_ipw_owner_integrity_error")
    }
  } else {
    validate_ipw_nuisance_view(link_view)
  }
  ctrl <- validate_ipw_owner_control(control)

  fact <- ipw_view_mask_compatibility_fact(nuisance_view, link_view, ctrl, construction_context)
  if (isTRUE(fact$decided)) return(invisible(fact$result))
  frame_owner <- fact$frame_owner
  n_a <- fact$n_a
  step <- fact$step
  p1 <- fact$p1
  p2 <- fact$p2
  a_cohort_idx <- frame_owner$a_cohort_index

  chunks <- 0L
  a_start <- 1L
  while (a_start <= n_a) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    hit1 <- ipw_owner_query_pair_mask(
      frame_owner, p1$s, p1$t,
      a_start = a_start, a_end = a_end,
      ds = p1$s_calendar_index, dt = p1$t_calendar_index
    )
    hit2 <- ipw_owner_query_pair_mask(
      frame_owner, p2$s, p2$t,
      a_start = a_start, a_end = a_end,
      ds = p2$s_calendar_index, dt = p2$t_calendar_index
    )
    r1 <- hit1$a_index[!is.na(ipw_owner_view_cohort_index(nuisance_view$cohorts, a_cohort_idx[hit1$a_index]))]
    r2 <- hit2$a_index[!is.na(ipw_owner_view_cohort_index(nuisance_view$cohorts, a_cohort_idx[hit2$a_index]))]
    if (!identical(r1, r2)) {
      cdid_ipw_owner_abort(
        ipw_owner_mask_mismatch_message(),
        c("cdid_ipw_owner_mask_mismatch_error", "cdid_ipw_shared_mask_error", "cdid_ipw_shared_incompatible_error"),
        diagnostics = list(chunk_start = a_start, chunk_end = a_end)
      )
    }
    chunks <- chunks + 1L
    a_start <- as.double(a_end) + 1
  }

  invisible(list(
    compatible = TRUE,
    n_a_units = n_a,
    chunks_compared = chunks,
    endpoints_read = 0L
  ))
}

# ---------------------------------------------------------------------------
# Outcome-free pair-observation mask fingerprints (CP08 P1)
# ---------------------------------------------------------------------------
# Computed from the same chunked set that the streamed comparison above
# compares: A-unit positions observed at both pair dates whose cohort is among
# the view's declared cohorts. The fingerprint is the count and the first three
# power sums of these positions modulo the prime P = 2^26 - 5. Positions are
# reduced below 2^26 before multiplying and block sums stay below 2^52, so every
# step is exact in binary64 and the value is a deterministic function of the
# exact position set: different fingerprints prove different masks. Equal
# fingerprints prove nothing (equal-power-sum sets exist); callers must resolve
# them with `verify_ipw_view_mask_compatibility`. Zero outcome reads.

ipw_owner_mask_fingerprint_modulus <- 67108859

# Exact sum modulo P of nonnegative integer-valued doubles below P < 2^26.
ipw_owner_exact_mod_sum <- function(v, P) {
  n <- length(v)
  if (n <= 67108864) return(sum(v) %% P)
  acc <- 0
  i <- 1
  while (i <= n) {
    j <- min(n, i + 67108863)
    acc <- (acc + sum(v[i:j]) %% P) %% P
    i <- j + 1
  }
  acc
}

# Internal streaming primitive; callers validate the owner or view and charge
# the stream before calling. O(step) scratch per chunk, nothing retained.
ipw_owner_pair_mask_fingerprint <- function(frame_owner, pair_desc, owner_cohort_index, step) {
  n_a <- frame_owner$n_a_units
  P <- ipw_owner_mask_fingerprint_modulus
  a_cohort_idx <- frame_owner$a_cohort_index
  oci <- owner_cohort_index
  fp <- c(count = 0, s1 = 0, s2 = 0, s3 = 0)
  chunks <- 0
  a_start <- 1L
  while (a_start <= n_a) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    hit <- ipw_owner_query_pair_mask(
      frame_owner, pair_desc$s, pair_desc$t,
      a_start = a_start, a_end = a_end,
      ds = pair_desc$s_calendar_index, dt = pair_desc$t_calendar_index
    )
    a <- hit$a_index
    if (length(a) > 0L && length(oci) > 0L) {
      # Exact membership in the sorted distinct owner cohort indices: the same
      # set `ipw_owner_view_cohort_index` resolves by scalar binary search.
      # The sealed descriptor already supplies sorted owner indices. Reuse
      # the existing bounded binary lookup: findInterval would revalidate
      # every breakpoint in every chunk, adding an uncharged O(m) pass.
      a <- a[!is.na(ipw_owner_match_sorted_double(a_cohort_idx[a], oci))]
      if (length(a) > 0L) {
        u <- as.double(a) %% P
        u2 <- (u * u) %% P
        fp[["count"]] <- fp[["count"]] + length(a)
        fp[["s1"]] <- (fp[["s1"]] + ipw_owner_exact_mod_sum(u, P)) %% P
        fp[["s2"]] <- (fp[["s2"]] + ipw_owner_exact_mod_sum(u2, P)) %% P
        fp[["s3"]] <- (fp[["s3"]] + ipw_owner_exact_mod_sum((u2 * u) %% P, P)) %% P
      }
    }
    chunks <- chunks + 1
    a_start <- as.double(a_end) + 1
  }
  list(fingerprint = fp, chunks = chunks)
}

# Fingerprint of a sealed nuisance view's mask at its defining pair, or at an
# explicit pair `(s, t)` under the SAME frame and ordered cohorts (the mask the
# streamed comparison would read for a link at that pair). The stream is
# charged with the one-comparison envelope above before any query.
fingerprint_ipw_view_mask <- function(nuisance_view, control, s = NULL, t = NULL,
                                      retained_entries = 0, construction_context = NULL) {
  # A template the bound context constructed and checked in this invocation
  # is fingerprinted without re-validating the same sealed snapshot per row.
  if (!ipw_owner_context_nuisance_checked(construction_context, nuisance_view)) validate_ipw_nuisance_view(nuisance_view)
  ctrl <- validate_ipw_owner_control(control)
  if (is.null(s) != is.null(t)) {
    cdid_ipw_owner_abort("Supply both `s` and `t`, or neither.", "cdid_ipw_owner_input_error")
  }
  frame_owner <- get("frame_owner", envir = nuisance_view$frame_ref, inherits = FALSE)
  n_a <- frame_owner$n_a_units
  m <- nuisance_view$cohorts$m_cohorts
  step <- min(n_a, ctrl$max_chunk_rows)
  if (64 * step + 32 * m + 256 > ctrl$max_matrix_entries ||
      128 * n_a * (1 + ceiling(log2(m + 1))) +
      128 * ceiling(n_a / step) * (1 + ceiling(log2(n_a + 1))) > ctrl$max_work_ops) {
    cdid_ipw_owner_abort(
      "Streamed mask fingerprint exceeds explicit `control` bounds.",
      "cdid_ipw_owner_bound_error"
    )
  }
  retained <- ipw_owner_int_bound(retained_entries, "retained_entries", min_val = 0)
  ipw_owner_safe_sum(c(retained, 64 * step + 32 * m + 256), ctrl$max_matrix_entries,
    "concurrent mask fingerprint workspace")
  pair_desc <- if (is.null(s)) {
    nuisance_view$defining_pair
  } else {
    ipw_owner_refuse_decorative_attrs(s, "s")
    ipw_owner_refuse_decorative_attrs(t, "t")
    ipw_owner_make_pair_descriptor(get("panel_owner", envir = frame_owner$panel_ref, inherits = FALSE), s, t)
  }
  oci <- nuisance_view$cohorts$owner_cohort_index
  res <- NULL
  # Trusted registry loop only: every complete pair has the full A-roster as
  # its mask, so one streamed fingerprint per cohort index and chunk step serves
  # every complete pair of this invocation. Other pairs stream as before.
  if (!is.null(construction_context) && ipw_owner_construction_context_for(construction_context, frame_owner) &&
      ipw_owner_pair_complete(frame_owner, pair_desc, step, construction_context$dates)) {
    key <- ipw_owner_context_fingerprint_key(oci, step)
    rec <- get0(key, envir = construction_context$fingerprints, inherits = FALSE)
    if (!is.null(rec) && identical(rec$owner_cohort_index, oci) && identical(rec$step, step)) {
      ipw_owner_check_interrupt()
      construction_context$fingerprint_reuses <- construction_context$fingerprint_reuses + 1
      res <- rec$result
    } else {
      res <- ipw_owner_pair_mask_fingerprint(frame_owner, pair_desc, oci, step)
      assign(key, list(owner_cohort_index = oci, step = step, result = res), envir = construction_context$fingerprints)
    }
  }
  if (is.null(res)) res <- ipw_owner_pair_mask_fingerprint(frame_owner, pair_desc, oci, step)
  list(
    fingerprint = res$fingerprint,
    pair = pair_desc,
    n_a_units = n_a,
    chunks_streamed = res$chunks,
    endpoints_read = 0L
  )
}

construct_ipw_link_view <- function(
  frame_owner,
  s,
  t,
  g = NULL,
  K = NULL,
  nuisance_view = NULL,
  control,
  loading_beta_candidate = NULL,
  candidates = NULL,
  link_id = NULL,
  functional_key = NULL,
  require_support = TRUE,
  verify_nuisance_mask = TRUE,
  pair_observed = NULL,
  risk_mask = NULL,
  observed_mask = NULL,
  cohort_index = NULL,
  construction_context = NULL
) {
  ipw_owner_require_deps()
  if (missing(frame_owner) || missing(s) || missing(t) || missing(control)) {
    cdid_ipw_owner_abort(
      "Required arguments (`frame_owner`, `s`, `t`, `control`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }

  # The exact sealed frame owner a bound construction context validated in
  # this invocation is not re-validated per link; any other object is.
  if (!ipw_owner_context_frame_checked(construction_context, frame_owner)) validate_ipw_frame_owner(frame_owner)
  ctrl <- validate_ipw_owner_control(control)
  ipw_owner_refuse_dense_masks(
    pair_observed = pair_observed,
    risk_mask = risk_mask,
    observed_mask = observed_mask,
    cohort_index = cohort_index,
    caller = "construct_ipw_link_view"
  )
  if (!gmm_scalar_logical(require_support) || !gmm_scalar_logical(verify_nuisance_mask)) {
    cdid_ipw_owner_abort("`require_support` and `verify_nuisance_mask` must be scalar logicals.", "cdid_ipw_owner_input_error")
  }
  if (!is.null(link_id) && !gmm_scalar_character(link_id)) {
    cdid_ipw_owner_abort("`link_id` must be NULL or a nonempty character scalar.", "cdid_ipw_owner_input_error")
  }
  if (!is.null(functional_key) && !gmm_scalar_character(functional_key)) {
    cdid_ipw_owner_abort("`functional_key` must be NULL or a nonempty character scalar.", "cdid_ipw_owner_input_error")
  }

  for (item in list(
    list(v = s, nm = "s"), list(v = t, nm = "t"),
    list(v = g, nm = "g"), list(v = K, nm = "K"),
    list(v = loading_beta_candidate, nm = "loading_beta_candidate"),
    list(v = candidates, nm = "candidates")
  )) {
    ipw_owner_refuse_decorative_attrs(item$v, item$nm)
  }

  if (!is.null(nuisance_view)) {
    # One checked nuisance/support template per view this invocation built:
    # the sealed snapshot recorded by the bound context is validated once, at
    # its construction, and every link of its operator reuses that check.
    # A view the context did not construct, a stale copy or a foreign frame
    # validates in full, with the original error classes and order.
    if (!ipw_owner_context_nuisance_checked(construction_context, nuisance_view, frame_owner)) {
      validate_ipw_nuisance_view(nuisance_view, frame_owner = frame_owner)
    }
    if (is.null(g)) g <- nuisance_view$cohorts$g
    if (is.null(K)) K <- nuisance_view$cohorts$K
  } else if (is.null(g) || is.null(K)) {
    cdid_ipw_owner_abort(
      "Supply either (`g`, `K`) or a `nuisance_view` defining the ordered cohorts.",
      "cdid_ipw_owner_cohort_error"
    )
  }

  p <- frame_owner$n_covariates
  m <- 1L + (if (is.atomic(K) && is.null(dim(K))) length(K) else 0L)
  has_cand <- !is.null( candidates )

  # Preflight BEFORE scanning postings or candidate vectors
  preflight <- preflight_ipw_link_view(
    n_a_units = frame_owner$n_a_units,
    p_covariates = p,
    m_cohorts = max(2L, m),
    control = ctrl,
    has_candidates = has_cand
  )
  step <- preflight$effective_chunk_rows

  if (!is.null(nuisance_view)) {
    if (!identical(as.double(g), nuisance_view$cohorts$g) ||
        !identical(as.double(K), nuisance_view$cohorts$K)) {
      cdid_ipw_owner_abort(
        "Link view `g` and exact ordered `K` disagree with the supplied `nuisance_view`.",
        c("cdid_ipw_owner_cohort_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  }
  panel_owner <- get("panel_owner", envir = frame_owner$panel_ref, inherits = FALSE)
  pair_desc <- ipw_owner_make_pair_descriptor(panel_owner, s, t)
  coh_desc <- if (!is.null(nuisance_view)) {
    nuisance_view$cohorts
  } else {
    ipw_owner_make_cohort_descriptor(frame_owner, g, K, require_support = require_support)
  }

  tilde_lambda <- NULL
  if (!is.null(loading_beta_candidate)) {
    if (!is.atomic(loading_beta_candidate) || !is.null(dim(loading_beta_candidate)) ||
        is.object(loading_beta_candidate) || !is.numeric(loading_beta_candidate) ||
        length(loading_beta_candidate) != p || anyNA(loading_beta_candidate) ||
        any(!is.finite(loading_beta_candidate))) {
      cdid_ipw_owner_abort(
        paste0("`loading_beta_candidate` must be a finite numeric vector of length `n_covariates` (", p, ")."),
        c("cdid_ipw_owner_candidate_error", "cdid_ipw_enclosure_candidate_error")
      )
    }
    tilde_lambda <- as.double(loading_beta_candidate)
  }

  cand_norm <- if (!is.null(candidates)) {
    ipw_enclosure_validate_candidates(candidates, p = p, m = coh_desc$m_cohorts)
  } else NULL
  if (!is.null(cand_norm$loading_beta) && !is.null(tilde_lambda) &&
      !identical(cand_norm$loading_beta, tilde_lambda)) {
    cdid_ipw_owner_abort(
      "`candidates$loading_beta` disagrees with `loading_beta_candidate`.",
      c("cdid_ipw_owner_candidate_error", "cdid_ipw_enclosure_candidate_error")
    )
  }

  sup_info <- NULL
  if (!is.null(nuisance_view) && !is.null(construction_context) &&
      ipw_owner_construction_context_for(construction_context, frame_owner)) {
    sup_info <- ipw_owner_context_link_support(construction_context, frame_owner, nuisance_view,
      pair_desc, coh_desc, step, require_support)
  }
  # CP08 stage 5: a nuisance-free view of one CP04 registry invocation may
  # borrow an exact support record of that invocation (see
  # ipw_owner_context_link_lookup); every other view scans as before. Lookup
  # and record read at most this pair's two posting dates, once per call.
  link_context <- is.null(nuisance_view) && !is.null(construction_context) &&
    ipw_owner_context_link_admits(construction_context, frame_owner, preflight, ctrl)
  if (link_context) {
    sup_info <- ipw_owner_context_link_lookup(construction_context, frame_owner, pair_desc, coh_desc, step,
      require_support)
  }
  link_scanned <- is.null(sup_info)
  if (is.null(sup_info)) {
    if (link_context) construction_context$link_scans <- construction_context$link_scans + 1
    sup_info <- ipw_owner_chunked_pair_support(
      frame_owner = frame_owner,
      pair_desc = pair_desc,
      coh_desc = coh_desc,
      step = step,
      require_support = require_support
    )
  }

  link_view <- structure(
    list(
      link_id = link_id,
      functional_key = functional_key,
      frame_ref = frame_owner$ref,
      panel_ref = frame_owner$panel_ref,
      nuisance_ref = if (!is.null(nuisance_view)) nuisance_view$ref else NULL,
      nuisance_mask_verified = is.null(nuisance_view),
      pair = pair_desc,
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
      loading_beta_candidate = tilde_lambda,
      candidates = cand_norm,
      control = ctrl,
      preflight = preflight,
      accounting = list(
        per_link_endpoint_copies_retained = FALSE,
        full_n_position_map_retained = FALSE,
        full_n_masks_retained = FALSE,
        endpoints_read = 0L
      ),
      work = list(
        support_chunks_processed = logit_cert_count_output(sup_info$chunks),
        actual_work_ops = sup_info$ops
      )
    ),
    class = "cdid_ipw_link_view"
  )

  # A link bound to a validated nuisance view runs the verifier's exact
  # owner/cohort/bound guards on this fresh invocation-local object before it
  # is sealed. When the posting fact decides mask equality, that is the very
  # result `verify_ipw_view_mask_compatibility` would return for the sealed
  # link, so the completed verification is recorded and the link sealed once.
  # An undecided fact keeps the original sequence: seal, streamed
  # verification of the sealed snapshot, then the verified re-seal.
  fact_verified <- FALSE
  if (!is.null(nuisance_view) && isTRUE(verify_nuisance_mask)) {
    fact <- ipw_view_mask_compatibility_fact(nuisance_view, link_view, ctrl, construction_context)
    if (isTRUE(fact$decided)) {
      link_view$nuisance_mask_verified <- TRUE
      fact_verified <- TRUE
    }
    fact <- NULL
  }

  sealed_link <- ipw_owner_seal(link_view, "link_view")
  if (link_context && link_scanned) {
    ipw_owner_context_record_link(construction_context, frame_owner, sealed_link, sup_info, step, require_support)
  }
  if (!is.null(nuisance_view) && isTRUE(verify_nuisance_mask) && !fact_verified) {
    verify_ipw_view_mask_compatibility(nuisance_view, sealed_link, ctrl, construction_context = construction_context)
    # Record only a completed actual mask verification, bound to this immutable
    # view and nuisance snapshot. An unchecked reference alone proves no reuse.
    checked <- sealed_link
    checked$ref <- NULL
    checked$nuisance_mask_verified <- TRUE
    sealed_link <- ipw_owner_seal(checked, "link_view")
  }
  sealed_link
}

declare_ipw_link_view <- construct_ipw_link_view

validate_ipw_link_view <- function(link_view, frame_owner = NULL, nuisance_view = NULL) {
  ipw_owner_require_deps()
  ipw_link_view_validate_body(link_view, frame_owner = frame_owner, nuisance_view = nuisance_view,
    nuisance_checked = NULL)
}

# Invocation-local record of sealed nuisance references one caller has already
# validated against one sealed frame owner within the current call. Entries
# bind the exact locked reference environment and its sealed snapshot; a later
# membership test requires both the identical environment and the aliased
# snapshot, so a stale, forged or relabelled object never matches. Holds two
# slots per distinct reference; callers charge that before use. Never stored
# beyond the call and never exposed through a public argument.
ipw_link_view_nuisance_checked_set <- function(frame_owner) {
  set <- new.env(parent = emptyenv())
  set$frame_ref <- frame_owner$ref
  set$refs <- list()
  set$objs <- list()
  set
}

# A reusable reference must name the current object's own immutable aliased
# snapshot, rather than another mutable or same-label container holding it.
ipw_owner_immutable_view_snapshot <- function(ref, obj, binding) {
  is.environment(ref) && environmentIsLocked(ref) &&
    exists(binding, envir = ref, inherits = FALSE) && bindingIsLocked(binding, ref) &&
    is.list(obj) && identical(obj$ref, ref) && ipw_owner_seal_ok(obj, binding)
}

ipw_link_view_nuisance_checked_has <- function(set, nref, n_obj, fref) {
  if (!is.environment(set) || !identical(set$frame_ref, fref)) return(FALSE)
  # The set is an index, never proof of validity. Its caller has already
  # validated the link's actual frame; retain every remaining predicate of
  # validate_ipw_nuisance_view on the current sealed snapshot before reuse.
  if (!is.list(n_obj) || !inherits(n_obj, "cdid_ipw_nuisance_view") ||
      !identical(n_obj$ref, nref) || !ipw_owner_seal_ok(n_obj, "nuisance_view") ||
      !identical(n_obj$frame_ref, fref) || !is.environment(fref) || !environmentIsLocked(fref) ||
      !exists("frame_owner", envir = fref, inherits = FALSE) || !bindingIsLocked("frame_owner", fref)) return(FALSE)
  frame <- get("frame_owner", envir = fref, inherits = FALSE)
  if (!identical(n_obj$panel_ref, frame$panel_ref)) return(FALSE)
  refs <- set$refs
  if (!is.list(refs) || !is.list(set$objs) || length(refs) != length(set$objs)) return(FALSE)
  for (k in seq_along(refs)) {
    if (identical(refs[[k]], nref)) return(gmm_sexp_aliases(set$objs[[k]], n_obj))
  }
  FALSE
}

ipw_link_view_nuisance_checked_add <- function(set, nref, n_obj, fref) {
  if (!is.environment(set) || !identical(set$frame_ref, fref)) return(invisible(FALSE))
  if (!is.list(set$refs) || !is.list(set$objs) || length(set$refs) != length(set$objs)) return(invisible(FALSE))
  for (j in seq_along(set$refs)) {
    if (identical(set$refs[[j]], nref) && gmm_sexp_aliases(set$objs[[j]], n_obj)) return(invisible(TRUE))
  }
  k <- length(set$refs) + 1L
  set$refs[[k]] <- nref
  set$objs[[k]] <- n_obj
  invisible(TRUE)
}

# Full link-view validation body. Every link keeps its own seal, frame and
# panel identity, nuisance reference, binding and completed-mask checks. Only
# the recursive validation of the referenced nuisance snapshot is performed
# once per exact sealed reference when `nuisance_checked` (a set created by
# `ipw_link_view_nuisance_checked_set` for the same frame owner in the same
# call) already holds it; the public validator passes NULL and validates in full.
ipw_link_view_validate_body <- function(link_view, frame_owner = NULL, nuisance_view = NULL,
                                        nuisance_checked = NULL) {
  ipw_owner_require_deps()
  if (missing(link_view) || is.null(link_view) ||
      !inherits(link_view, "cdid_ipw_link_view") ||
      !is.list(link_view) || !ipw_owner_seal_ok(link_view, "link_view")) {
    cdid_ipw_owner_abort(
      "Link view is malformed, foreign, or differs from its sealed snapshot; reuse is refused.",
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_view_error")
    )
  }
  fref <- link_view$frame_ref
  if (!is.environment(fref) || !environmentIsLocked(fref) ||
      !exists("frame_owner", envir = fref, inherits = FALSE) ||
      !bindingIsLocked("frame_owner", fref)) {
    cdid_ipw_owner_abort(
      "Link view lacks a valid sealed `frame_owner` reference.",
      c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_view_error")
    )
  }
  f_obj <- get("frame_owner", envir = fref, inherits = FALSE)
  validate_ipw_frame_owner(f_obj)
  if (!identical(link_view$panel_ref, f_obj$panel_ref)) {
    cdid_ipw_owner_abort(
      "Link view `panel_ref` disagrees with its `frame_owner`.",
      c("cdid_ipw_owner_mismatch_error", "cdid_ipw_owner_view_error")
    )
  }
  if (!is.null(frame_owner)) {
    if (!gmm_sexp_aliases(frame_owner, f_obj)) {
      validate_ipw_frame_owner(frame_owner)
    }
    if (!identical(fref, frame_owner$ref)) {
      cdid_ipw_owner_abort(
        "Link view belongs to a different `cdid_ipw_frame_owner` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_owner_view_error", "cdid_ipw_context_mismatch_error")
      )
    }
  }
  n_obj <- NULL
  if (!is.null(link_view$nuisance_ref)) {
    nref <- link_view$nuisance_ref
    if (!is.environment(nref) || !environmentIsLocked(nref) ||
        !exists("nuisance_view", envir = nref, inherits = FALSE)) {
      cdid_ipw_owner_abort(
        "Link view `nuisance_ref` is malformed.",
        c("cdid_ipw_owner_integrity_error", "cdid_ipw_owner_view_error")
      )
    }
    n_obj <- get("nuisance_view", envir = nref, inherits = FALSE)
    if (!ipw_link_view_nuisance_checked_has(nuisance_checked, nref, n_obj, fref)) {
      validate_ipw_nuisance_view(n_obj, frame_owner = f_obj)
      ipw_link_view_nuisance_checked_add(nuisance_checked, nref, n_obj, fref)
    }
  }
  if (!is.null(nuisance_view)) {
    if (is.null(n_obj) || !gmm_sexp_aliases(nuisance_view, n_obj)) {
      validate_ipw_nuisance_view(nuisance_view, frame_owner = f_obj)
    }
    if (!identical(link_view$nuisance_ref, nuisance_view$ref)) {
      cdid_ipw_owner_abort(
        "Link view was not bound to the supplied `nuisance_view` snapshot.",
        c("cdid_ipw_owner_mismatch_error", "cdid_ipw_owner_view_error")
      )
    }
    if (!isTRUE(link_view$nuisance_mask_verified)) {
      cdid_ipw_owner_abort(
        "Nuisance reference is present but actual risk-set observer-mask equality was not verified; reuse is refused.",
        c("cdid_ipw_owner_mask_unverified_error", "cdid_ipw_owner_view_error", "cdid_ipw_shared_incompatible_error")
      )
    }
  }
  invisible(link_view)
}

# Explicitly bounded and disclosed verification of external dense O(N) masks
# against a compact posting-derived view (dense masks are never retained by the view).
verify_ipw_view_dense_masks <- function(
  view,
  control,
  pair_observed = NULL,
  risk_mask = NULL,
  observed_mask = NULL,
  cohort_index = NULL
) {
  if (inherits(view, "cdid_ipw_link_view")) {
    validate_ipw_link_view(view)
    pair_desc <- view$pair
  } else if (inherits(view, "cdid_ipw_nuisance_view")) {
    validate_ipw_nuisance_view(view)
    pair_desc <- view$defining_pair
  } else {
    cdid_ipw_owner_abort("Expected a sealed `cdid_ipw_link_view` or `cdid_ipw_nuisance_view`.", "cdid_ipw_owner_view_error")
  }
  ctrl <- validate_ipw_owner_control(control)
  frame_owner <- get("frame_owner", envir = view$frame_ref, inherits = FALSE)
  n_a <- frame_owner$n_a_units

  n_masks <- (!is.null(pair_observed)) + (!is.null(risk_mask)) + (!is.null(observed_mask)) + (!is.null(cohort_index))
  if (n_masks == 0L) {
    cdid_ipw_owner_abort("Supply at least one dense mask to `verify_ipw_view_dense_masks`.", "cdid_ipw_owner_input_error")
  }

  # Preflight O(N) caller-supplied dense mask verification BEFORE scanning mask contents
  caller_dense_entries <- ipw_owner_safe_product(n_masks, n_a, ctrl$max_matrix_entries, "disclosed caller dense mask entries")
  step <- min(n_a, ctrl$max_chunk_rows)
  work_bound <- ipw_owner_safe_sum(c(
    ipw_owner_safe_product(128 * (1 + ceiling(log2(view$cohorts$m_cohorts + 1))), n_a, ctrl$max_work_ops, "dense mask and compact cohort work"),
    ipw_owner_safe_product(128 * (1 + ceiling(log2(n_a + 1))), ceiling(n_a / step), ctrl$max_work_ops, "posting search work")
  ), ctrl$max_work_ops, "dense mask verification work")
  ipw_owner_safe_sum(c(caller_dense_entries, 64 * step, 32 * view$cohorts$m_cohorts, 256), ctrl$max_matrix_entries, "dense mask verification concurrent entries")

  for (m_item in list(
    list(x = pair_observed, nm = "pair_observed", lgl = TRUE),
    list(x = risk_mask, nm = "risk_mask", lgl = TRUE),
    list(x = observed_mask, nm = "observed_mask", lgl = TRUE),
    list(x = cohort_index, nm = "cohort_index", lgl = FALSE)
  )) {
    x <- m_item$x
    if (is.null(x)) next
    if (!is.atomic(x) || !is.null(dim(x)) || is.object(x) || length(x) != n_a ||
        (m_item$lgl && !is.logical(x)) || (!m_item$lgl && !is.numeric(x))) {
      cdid_ipw_owner_abort(
        paste0("Supplied dense `", m_item$nm, "` must be a plain atomic vector of length `n_a_units` (", n_a, ")."),
        c("cdid_ipw_owner_dense_mask_error", "cdid_ipw_enclosure_mask_error")
      )
    }
  }

  a_cohort_idx <- frame_owner$a_cohort_index
  chunks <- 0L
  a_start <- 1L
  while (a_start <= n_a) {
    ipw_owner_check_interrupt()
    a_end <- as.integer(min(as.double(n_a), as.double(a_start) + step - 1))
    idx <- seq.int(a_start, a_end)
    b_len <- length(idx)

    ci_ch <- ipw_owner_view_cohort_index(view$cohorts, a_cohort_idx[idx])
    in_r_ch <- !is.na(ci_ch)
    hit <- ipw_owner_query_pair_mask(
      frame_owner, pair_desc$s, pair_desc$t,
      a_start = a_start, a_end = a_end,
      ds = pair_desc$s_calendar_index, dt = pair_desc$t_calendar_index
    )
    po_ch <- logical(b_len)
    if (length(hit$a_index) > 0L) {
      po_ch[hit$a_index - a_start + 1L] <- TRUE
    }
    in_o_ch <- in_r_ch & po_ch

    if (!is.null(risk_mask)) {
      rm_ch <- risk_mask[idx]
      if (anyNA(rm_ch) || !identical(as.logical(rm_ch), in_r_ch)) {
        cdid_ipw_owner_abort(
          "Supplied dense `risk_mask` disagrees with posting-derived risk set R = {G in cohorts}.",
          c("cdid_ipw_owner_dense_mask_error", "cdid_ipw_enclosure_mask_error"),
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
    }
    if (!is.null(observed_mask)) {
      om_ch <- observed_mask[idx]
      if (anyNA(om_ch) || !identical(as.logical(om_ch), in_o_ch)) {
        cdid_ipw_owner_abort(
          "Supplied dense `observed_mask` disagrees with posting-derived observer set O = R intersect S.",
          c("cdid_ipw_owner_dense_mask_error", "cdid_ipw_enclosure_mask_error"),
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
    }
    if (!is.null(cohort_index)) {
      ci_arg <- cohort_index[idx]
      if (any(is.na(ci_arg) != is.na(ci_ch)) ||
          any(!is.finite(ci_arg[!is.na(ci_arg)]) | ci_arg[!is.na(ci_arg)] != floor(ci_arg[!is.na(ci_arg)])) ||
          !identical(as.integer(ci_arg), ci_ch)) {
        cdid_ipw_owner_abort(
          "Supplied dense `cohort_index` disagrees with `match(G, cohorts)`.",
          c("cdid_ipw_owner_dense_mask_error", "cdid_ipw_enclosure_mask_error"),
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
    }
    if (!is.null(pair_observed)) {
      po_arg <- pair_observed[idx]
      if (anyNA(po_arg) || !identical(as.logical(po_arg[in_r_ch]), po_ch[in_r_ch])) {
        cdid_ipw_owner_abort(
          "Supplied dense `pair_observed` disagrees with posting-derived mask on the risk set R.",
          c("cdid_ipw_owner_dense_mask_error", "cdid_ipw_enclosure_mask_error"),
          diagnostics = list(chunk_start = a_start, chunk_end = a_end)
        )
      }
    }

    chunks <- chunks + 1L
    a_start <- as.double(a_end) + 1
  }

  list(
    verified = TRUE,
    dense_mask_disclosed = TRUE,
    dense_mask_retained_by_view = FALSE,
    caller_dense_mask_entries = caller_dense_entries,
    chunks_verified = chunks,
    work_ops_bound = work_bound,
    endpoints_read = 0L
  )
}

# ---------------------------------------------------------------------------
# Component 4: Bounded Contiguous A-Range Readers (Risk-First Filtering!)
# ---------------------------------------------------------------------------

# The ONLY function that reads `panel_owner$observations$y`.
# Invoked strictly AFTER risk-set and pair-observation position filtering.
# Bounded diagnostic state: only the latest call, never a growing history.
# Callers may accumulate per-chunk read_diagnostics in their own callback.
ipw_owner_reset_trace <- function(trace_env) {
  if (!is.null(trace_env)) {
    if (!is.environment(trace_env)) cdid_ipw_owner_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_owner_input_error")
    trace_env$calls <- 0L
    trace_env$pos_s_index <- integer(0L)
    trace_env$pos_t_index <- integer(0L)
  }
  invisible(NULL)
}

ipw_owner_read_y_positions <- function(panel_owner, pos_s, pos_t, trace_env = NULL) {
  if (!is.null(trace_env)) {
    if (!is.environment(trace_env)) cdid_ipw_owner_abort("`trace_env` must be NULL or an environment.", "cdid_ipw_owner_input_error")
    trace_env$calls <- 1L
    trace_env$pos_s_index <- pos_s
    trace_env$pos_t_index <- pos_t
  }
  oy <- panel_owner$observations$y
  y_s <- oy[pos_s]
  y_t <- oy[pos_t]
  if (anyNA(y_s) || any(!is.finite(y_s)) || anyNA(y_t) || any(!is.finite(y_t))) {
    cdid_ipw_owner_abort(
      "Observed risk-set outcome endpoints contain NA, NaN, or Inf.",
      c("cdid_ipw_owner_endpoint_error", "cdid_ipw_enclosure_endpoint_error")
    )
  }
  list(y_left_obs = as.double(y_s), y_right_obs = as.double(y_t))
}

# Sealed record of one view and its owners, made by a caller immediately after
# its own full validation of that view in the same call (point-IF normalization,
# point-link fit). A chunk query given the record re-checks only that the
# record, the view and both owner bindings are still those same objects instead
# of re-running every validator on every chunk; any other view, record or
# binding takes the full validation path. The record also carries the view
# pair's complete-history fact and an exact memo of validated controls and
# range preflights (values of the unchanged functions for identical inputs).
ipw_owner_validated_view_record <- function(view, frame_owner, panel_owner, memo = NULL, step = NULL) {
  pair <- if (inherits(view, "cdid_ipw_link_view")) view$pair else view$defining_pair
  if (!is.environment(memo)) memo <- new.env(parent = emptyenv())
  ipw_owner_seal(list(
    kind = "cdid_ipw_validated_view",
    view = view,
    frame_owner = frame_owner,
    panel_owner = panel_owner,
    complete_pair = identical(frame_owner$n_a_units, panel_owner$n_a_units) &&
      ipw_owner_pair_complete(panel_owner, pair, step, memo),
    memo = memo
  ), "validated_view")
}

ipw_owner_validated_view_match <- function(validated, view) {
  ok <- is.list(validated) && identical(validated$kind, "cdid_ipw_validated_view") &&
    is.environment(validated$ref) && environmentIsLocked(validated$ref) &&
    identical(get0("validated_view", envir = validated$ref, inherits = FALSE), validated) &&
    identical(validated$view, view) && is.list(view) &&
    is.environment(view$frame_ref) && is.environment(view$panel_ref) &&
    identical(get0("frame_owner", envir = view$frame_ref, inherits = FALSE), validated$frame_owner) &&
    identical(get0("panel_owner", envir = view$panel_ref, inherits = FALSE), validated$panel_owner) &&
    is.environment(validated$memo)
  isTRUE(ok)
}

ipw_owner_validated_control <- function(validated, control) {
  memo <- validated$memo
  if (!is.null(memo$control_in) && identical(memo$control_in, control)) return(memo$control_out)
  out <- validate_ipw_owner_control(control)
  memo$control_in <- control
  memo$control_out <- out
  memo$preflight <- NULL
  out
}

ipw_owner_validated_preflight <- function(validated, ctrl, b_len, p, m, read_outcomes, include_design,
                                          enclose_differences) {
  memo <- validated$memo
  args <- list(b_len = b_len, p = p, m = m, read_outcomes = read_outcomes, include_design = include_design,
    enclose_differences = enclose_differences)
  hit <- memo$preflight
  if (!is.null(hit) && identical(hit$args, args) && identical(hit$ctrl, ctrl)) return(hit$preflight)
  pre <- preflight_ipw_range_query(
    b_len = b_len,
    p_covariates = p,
    m_cohorts = m,
    control = ctrl,
    read_outcomes = read_outcomes,
    include_design = include_design,
    enclose_differences = enclose_differences
  )
  # One previous exact request only: bounded state even across many views.
  memo$preflight <- list(args = args, ctrl = ctrl, preflight = pre)
  pre
}

query_ipw_view_chunk <- function(
  view,
  a_start,
  a_end,
  control,
  read_outcomes = FALSE,
  include_design = FALSE,
  enclose_differences = FALSE,
  expected_a_id = NULL,
  trace_env = NULL,
  pair_observed = NULL,
  risk_mask = NULL,
  observed_mask = NULL,
  cohort_index = NULL,
  validated = NULL
) {
  ipw_owner_require_deps()
  if (missing(view) || missing(a_start) || missing(a_end) || missing(control)) {
    cdid_ipw_owner_abort(
      "Required arguments (`view`, `a_start`, `a_end`, `control`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }

  # 1. Validate sealed view and owner snapshots (O(1) pointer identity under R copy-on-modify)
  trusted <- !is.null(validated) && ipw_owner_validated_view_match(validated, view)
  if (trusted) {
    pair_desc <- if (inherits(view, "cdid_ipw_link_view")) view$pair else view$defining_pair
  } else if (inherits(view, "cdid_ipw_link_view")) {
    validate_ipw_link_view(view)
    pair_desc <- view$pair
  } else if (inherits(view, "cdid_ipw_nuisance_view")) {
    validate_ipw_nuisance_view(view)
    pair_desc <- view$defining_pair
  } else {
    cdid_ipw_owner_abort(
      "Expected a sealed `cdid_ipw_link_view` or `cdid_ipw_nuisance_view`.",
      "cdid_ipw_owner_view_error"
    )
  }
  ctrl <- if (trusted) ipw_owner_validated_control(validated, control) else validate_ipw_owner_control(control)
  ipw_owner_refuse_dense_masks(
    pair_observed = pair_observed,
    risk_mask = risk_mask,
    observed_mask = observed_mask,
    cohort_index = cohort_index,
    caller = "query_ipw_view_chunk"
  )
  if (!gmm_scalar_logical(read_outcomes) || !gmm_scalar_logical(include_design) ||
      !gmm_scalar_logical(enclose_differences)) {
    cdid_ipw_owner_abort(
      "`read_outcomes`, `include_design`, and `enclose_differences` must be scalar logicals.",
      "cdid_ipw_owner_input_error"
    )
  }
  if (isTRUE(enclose_differences) && !isTRUE(read_outcomes)) {
    cdid_ipw_owner_abort(
      "`enclose_differences = TRUE` requires `read_outcomes = TRUE`.",
      "cdid_ipw_owner_input_error"
    )
  }

  frame_owner <- if (trusted) validated$frame_owner else get("frame_owner", envir = view$frame_ref, inherits = FALSE)
  panel_owner <- if (trusted) validated$panel_owner else get("panel_owner", envir = view$panel_ref, inherits = FALSE)
  n_a <- frame_owner$n_a_units
  p <- frame_owner$n_covariates
  m <- view$cohorts$m_cohorts

  s_idx <- ipw_owner_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- ipw_owner_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > n_a) {
    cdid_ipw_owner_abort(
      paste0("Range `[a_start, a_end]` must satisfy `1 <= a_start <= a_end <= n_a_units` (", n_a, ")."),
      c("cdid_ipw_owner_bound_error", "cdid_ipw_bound_error")
    )
  }
  a_s <- as.integer(s_idx)
  a_e <- as.integer(e_idx)
  b_len <- a_e - a_s + 1L

  # Header check on optional expected_a_id BEFORE preflight and content slice
  if (!is.null(expected_a_id)) {
    ipw_owner_refuse_decorative_attrs(expected_a_id, "expected_a_id")
    if (is.factor(expected_a_id) || is.object(expected_a_id) || !is.atomic(expected_a_id) ||
        !is.null(dim(expected_a_id)) || !(length(expected_a_id) %in% c(b_len, n_a))) {
      cdid_ipw_owner_abort(
        "`expected_a_id` must be NULL or a plain atomic vector of length `b_len` or `n_a_units`.",
        c("cdid_ipw_owner_id_error", "cdid_ipw_owner_order_mismatch_error")
      )
    }
  }

  # 2. Preflight BEFORE any slice, ID verification, posting lookup, or outcome read
  preflight <- if (trusted) {
    ipw_owner_validated_preflight(validated, ctrl, b_len, p, m, read_outcomes, include_design, enclose_differences)
  } else {
    preflight_ipw_range_query(
      b_len = b_len,
      p_covariates = p,
      m_cohorts = m,
      control = ctrl,
      read_outcomes = read_outcomes,
      include_design = include_design,
      enclose_differences = enclose_differences
    )
  }
  ipw_owner_check_interrupt()

  a_seq <- seq.int(a_s, a_e)
  actual_ops <- 0

  # 3. Verify optional expected_a_id on the queried chunk ONLY (no full-N scan!)
  if (!is.null(expected_a_id)) {
    exp_sub <- if (length(expected_a_id) == b_len) expected_a_id else expected_a_id[a_seq]
    if (identical(panel_owner$id_type, "character") && is.character(exp_sub) && !anyNA(exp_sub)) {
      exp_bytes <- ipw_owner_id_byte_metadata(exp_sub)
      own_bytes <- ipw_owner_id_byte_metadata(panel_owner$a_id_key[a_seq])
      merged <- list(total = ipw_owner_safe_sum(c(exp_bytes$total, own_bytes$total), label = "paired expected ID bytes"),
                     max = max(exp_bytes$max, own_bytes$max))
      preflight <- ipw_owner_preflight_id_bytes(preflight, ctrl, merged,
        owner_max_bytes = 0, lookup_count = 0, search_width = 1, simultaneous = TRUE)
    }
    if (identical(panel_owner$id_type, "numeric")) {
      if (!is.numeric(exp_sub) || anyNA(exp_sub) || any(!is.finite(exp_sub))) {
        cdid_ipw_owner_abort(
          "Queried `expected_a_id` slice contains non-numeric, NA, or non-finite values.",
          c("cdid_ipw_owner_id_error", "cdid_ipw_owner_order_mismatch_error")
        )
      }
      exp_key <- as.double(exp_sub)
      exp_key[exp_key == 0.0] <- 0.0
      if (!identical(exp_key, panel_owner$a_id_key[a_seq])) {
        cdid_ipw_owner_abort(
          "Queried `expected_a_id` slice differs in value or unit order from `frame_owner$a_id`.",
          c("cdid_ipw_owner_order_mismatch_error", "cdid_ipw_owner_id_error")
        )
      }
    } else {
      if (!is.character(exp_sub) || anyNA(exp_sub) || any(!nzchar(exp_sub)) ||
          !identical(as.character(exp_sub), panel_owner$a_id_key[a_seq])) {
        cdid_ipw_owner_abort(
          "Queried `expected_a_id` slice differs in value, type, or unit order from `frame_owner$a_id`.",
          c("cdid_ipw_owner_order_mismatch_error", "cdid_ipw_owner_id_error")
        )
      }
    }
    actual_ops <- actual_ops + 2 * b_len
  }

  ipw_owner_reset_trace(trace_env) # after base AND optional ID-byte preflight
  # 4. Chunk-local row identities and risk membership (zero hash tables, zero coercions)
  u_idx <- panel_owner$a_unit_index[a_seq]
  u_id <- panel_owner$a_id[a_seq]
  u_g <- panel_owner$a_g[a_seq]
  u_w <- frame_owner$w_a[a_seq]
  cohort_index <- ipw_owner_view_cohort_index(view$cohorts, panel_owner$a_cohort_index[a_seq])
  in_risk <- !is.na(cohort_index)
  is_target <- in_risk & (cohort_index == 1L)
  actual_ops <- actual_ops + 4 * b_len

  # 5. Outcome-free posting-list intersection on [a_s, a_e], then filter by `in_risk`
  #    BEFORE reading any `y` values! A record's complete pair reads the same
  #    hits by position.
  hit <- if (trusted && isTRUE(validated$complete_pair)) {
    ipw_owner_query_complete_pair_mask(panel_owner, pair_desc$s_calendar_index, pair_desc$t_calendar_index, a_s, a_e)
  } else {
    ipw_owner_query_pair_mask(
      panel_owner,
      pair_desc$s,
      pair_desc$t,
      a_start = a_s,
      a_end = a_e,
      ds = pair_desc$s_calendar_index,
      dt = pair_desc$t_calendar_index
    )
  }
  hit_local_pos <- if (length(hit$a_index) > 0L) hit$a_index - a_s + 1L else integer(0L)
  in_risk_hit <- if (length(hit_local_pos) > 0L) which(in_risk[hit_local_pos]) else integer(0L)

  obs_local_pos <- hit_local_pos[in_risk_hit]
  obs_a_index <- hit$a_index[in_risk_hit]
  obs_unit_index <- hit$unit_index[in_risk_hit]
  obs_cohort_index <- cohort_index[obs_local_pos]
  obs_s_index <- hit$obs_s[in_risk_hit]
  obs_t_index <- hit$obs_t[in_risk_hit]
  n_obs_risk <- length(obs_local_pos)

  in_observed <- logical(b_len)
  if (n_obs_risk > 0L) in_observed[obs_local_pos] <- TRUE
  is_risk_nonobserver <- in_risk & !in_observed

  # Chunk-sized position map (length b_len <= B, NEVER length N!)
  chunk_to_obs_pos <- integer(b_len)
  if (n_obs_risk > 0L) chunk_to_obs_pos[obs_local_pos] <- seq_len(n_obs_risk)
  actual_ops <- actual_ops + 4 * b_len

  # 6. Read raw stored endpoints ONLY when `read_outcomes = TRUE` and ONLY on `O = R intersect S`
  y_left_obs <- NULL
  y_right_obs <- NULL
  Delta_iv_obs <- NULL
  n_endpoints_read <- 0L
  read_obs_s_idx <- integer(0L)
  read_obs_t_idx <- integer(0L)

  if (isTRUE(read_outcomes)) {
    ep <- ipw_owner_read_y_positions(
      panel_owner = panel_owner,
      pos_s = obs_s_index,
      pos_t = obs_t_index,
      trace_env = trace_env
    )
    y_left_obs <- ep$y_left_obs
    y_right_obs <- ep$y_right_obs
    read_obs_s_idx <- obs_s_index
    read_obs_t_idx <- obs_t_index
    n_endpoints_read <- logit_cert_count_output(2 * as.double(n_obs_risk))
    actual_ops <- actual_ops + 2 * n_obs_risk

    if (isTRUE(enclose_differences)) {
      Delta_iv_obs <- matrix(0.0, nrow = n_obs_risk, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
      for (k in seq_len(n_obs_risk)) {
        Delta_iv_obs[k, ] <- tryCatch(
          gmm_iv_sub(gmm_iv(y_right_obs[k]), gmm_iv(y_left_obs[k])),
          cdid_gmm_accuracy_error = function(e) {
            cdid_ipw_owner_abort(
              paste0("Outward interval endpoint difference overflowed binary64 on risk-set observer ", obs_a_index[k], "."),
              c("cdid_ipw_owner_outcome_overflow_error", "cdid_gmm_accuracy_error"),
              diagnostics = list(a_index = obs_a_index[k], y_left = y_left_obs[k], y_right = y_right_obs[k])
            )
          }
        )
      }
      actual_ops <- actual_ops + n_obs_risk
    }
  }

  Z_chunk <- NULL
  Z_obs <- NULL
  if (isTRUE(include_design)) {
    Z_chunk <- frame_owner$Z_a[a_seq, , drop = FALSE]
    Z_obs <- frame_owner$Z_a[obs_a_index, , drop = FALSE]
    actual_ops <- actual_ops + (b_len + n_obs_risk) * p
  }

  stopifnot(actual_ops <= preflight$work_ops_bound)

  structure(
    list(
      a_start = a_s,
      a_end = a_e,
      chunk_rows = b_len,
      a_index = a_seq,
      unit_index = u_idx,
      id = u_id,
      g = u_g,
      weight = u_w,
      cohort_index = cohort_index,
      in_risk = in_risk,
      is_target = is_target,
      in_observed = in_observed,
      is_risk_nonobserver = is_risk_nonobserver,
      chunk_to_obs_pos = chunk_to_obs_pos,
      observed_risk = list(
        local_pos = obs_local_pos,
        a_index = obs_a_index,
        unit_index = obs_unit_index,
        cohort_index = obs_cohort_index,
        obs_s_index = obs_s_index,
        obs_t_index = obs_t_index,
        y_left_obs = y_left_obs,
        y_right_obs = y_right_obs,
        Delta_iv_obs = Delta_iv_obs
      ),
      design = if (isTRUE(include_design)) list(Z_chunk = Z_chunk, Z_obs = Z_obs) else NULL,
      counts = list(
        queried_row_count = b_len,
        risk_row_count = sum(in_risk),
        observed_risk_row_count = n_obs_risk,
        nonobserver_risk_row_count = sum(is_risk_nonobserver),
        outside_risk_row_count = sum(!in_risk),
        all_cohort_pair_hit_count = length(hit$a_index),
        outside_risk_pair_hit_count = length(hit$a_index) - n_obs_risk
      ),
      read_diagnostics = list(
        read_outcomes = isTRUE(read_outcomes),
        n_endpoints_read = n_endpoints_read,
        read_obs_s_index = read_obs_s_idx,
        read_obs_t_index = read_obs_t_idx,
        risk_filtered_before_outcome_read = TRUE,
        off_risk_endpoints_read = FALSE,
        nonobserver_endpoints_read = FALSE,
        cached_rounded_dy_returned = FALSE
      ),
      preflight = preflight,
      work = list(
        actual_work_ops = actual_ops
      )
    ),
    class = "cdid_ipw_view_chunk"
  )
}

query_ipw_link_view_chunk <- query_ipw_view_chunk
query_ipw_nuisance_view_chunk <- query_ipw_view_chunk
read_ipw_view_chunk <- query_ipw_view_chunk

# Bounded chunk streaming over a view (preflights whole stream before first chunk/callback)
stream_ipw_view_chunks <- function(
  view,
  control,
  callback,
  read_outcomes = FALSE,
  include_design = FALSE,
  enclose_differences = FALSE,
  a_start = 1L,
  a_end = NULL,
  expected_a_id = NULL,
  trace_env = NULL
) {
  ipw_owner_require_deps()
  if (missing(view) || missing(control) || missing(callback)) {
    cdid_ipw_owner_abort(
      "Required arguments (`view`, `control`, `callback`) must be supplied explicitly.",
      "cdid_ipw_owner_input_error"
    )
  }
  if (!is.function(callback)) {
    cdid_ipw_owner_abort("`callback` must be a function.", "cdid_ipw_owner_input_error")
  }
  if (inherits(view, "cdid_ipw_link_view")) {
    validate_ipw_link_view(view)
  } else if (inherits(view, "cdid_ipw_nuisance_view")) {
    validate_ipw_nuisance_view(view)
  } else {
    cdid_ipw_owner_abort("Expected a sealed `cdid_ipw_link_view` or `cdid_ipw_nuisance_view`.", "cdid_ipw_owner_view_error")
  }
  ctrl <- validate_ipw_owner_control(control)
  frame_owner <- get("frame_owner", envir = view$frame_ref, inherits = FALSE)
  n_a <- frame_owner$n_a_units

  s_idx <- ipw_owner_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- if (is.null(a_end)) as.double(n_a) else ipw_owner_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > n_a) {
    cdid_ipw_owner_abort("Invalid `[a_start, a_end]` stream range.", "cdid_ipw_owner_bound_error")
  }
  a_s <- as.integer(s_idx)
  a_e <- as.integer(e_idx)
  total_rows <- a_e - a_s + 1L

  if (!is.null(expected_a_id)) {
    ipw_owner_refuse_decorative_attrs(expected_a_id, "expected_a_id")
    if (is.factor(expected_a_id) || is.object(expected_a_id) || !is.atomic(expected_a_id) ||
        !is.null(dim(expected_a_id)) || !(length(expected_a_id) %in% c(total_rows, n_a))) {
      cdid_ipw_owner_abort(
        "`expected_a_id` must be NULL or a plain atomic vector of length `total_rows` or `n_a_units`.",
        c("cdid_ipw_owner_id_error", "cdid_ipw_owner_order_mismatch_error")
      )
    }
  }

  # Whole-stream preflight BEFORE first slice or callback
  stream_pre <- preflight_ipw_stream(
    n_total_queried_rows = total_rows,
    max_chunk_rows_req = ctrl$max_chunk_rows,
    p_covariates = frame_owner$n_covariates,
    m_cohorts = view$cohorts$m_cohorts,
    control = ctrl,
    read_outcomes = read_outcomes,
    include_design = include_design,
    enclose_differences = enclose_differences
  )
  step <- stream_pre$effective_chunk_rows

  batches <- 0L
  emitted_rows <- 0L
  max_batch_rows <- 0L
  total_risk_rows <- 0L
  total_obs_risk_rows <- 0L
  total_nonobs_risk_rows <- 0L
  total_outside_risk_rows <- 0L
  total_endpoints_read <- 0
  total_ops <- 0

  cur_s <- a_s
  while (cur_s <= a_e) {
    ipw_owner_check_interrupt()
    cur_e <- as.integer(min(as.double(a_e), as.double(cur_s) + step - 1))
    exp_slice <- if (is.null(expected_a_id)) {
      NULL
    } else if (length(expected_a_id) == n_a) {
      expected_a_id[seq.int(cur_s, cur_e)]
    } else {
      expected_a_id[seq.int(cur_s - a_s + 1L, cur_e - a_s + 1L)]
    }

    ch <- query_ipw_view_chunk(
      view = view,
      a_start = cur_s,
      a_end = cur_e,
      control = ctrl,
      read_outcomes = read_outcomes,
      include_design = include_design,
      enclose_differences = enclose_differences,
      expected_a_id = exp_slice,
      trace_env = trace_env
    )
    batches <- batches + 1L
    ch$batch_index <- batches

    emitted_rows <- emitted_rows + ch$chunk_rows
    max_batch_rows <- max(max_batch_rows, ch$chunk_rows)
    total_risk_rows <- total_risk_rows + ch$counts$risk_row_count
    total_obs_risk_rows <- total_obs_risk_rows + ch$counts$observed_risk_row_count
    total_nonobs_risk_rows <- total_nonobs_risk_rows + ch$counts$nonobserver_risk_row_count
    total_outside_risk_rows <- total_outside_risk_rows + ch$counts$outside_risk_row_count
    total_endpoints_read <- total_endpoints_read + ch$read_diagnostics$n_endpoints_read
    total_ops <- total_ops + ch$work$actual_work_ops

    callback(ch) # Discard callback return value; let errors and interrupts propagate.
    ch <- NULL
    cur_s <- as.double(cur_e) + 1
  }

  stopifnot(emitted_rows == total_rows, total_ops <= stream_pre$total_work_ops_bound)

  invisible(list(
    a_start = a_s,
    a_end = a_e,
    queried_row_count = emitted_rows,
    n_batches = batches,
    max_batch_rows = max_batch_rows,
    total_risk_rows = total_risk_rows,
    total_observed_risk_rows = total_obs_risk_rows,
    total_nonobserver_risk_rows = total_nonobs_risk_rows,
    total_outside_risk_rows = total_outside_risk_rows,
    total_endpoints_read = logit_cert_count_output(total_endpoints_read),
    total_work_ops = total_ops,
    preflight = stream_pre
  ))
}

traverse_ipw_view_chunks <- stream_ipw_view_chunks

# Outcome-free query of pair observation positions on [a_start, a_end] directly
# from panel_owner$postings (for allowed-domain position comparison without reading y).
query_ipw_panel_pair_positions <- function(
  panel_owner,
  s,
  t,
  a_start = 1L,
  a_end = NULL,
  control
) {
  validate_ipw_panel_owner(panel_owner)
  ctrl <- validate_ipw_owner_control(control)
  n_a <- panel_owner$n_a_units
  if (n_a < 1L) {
    if (128 > ctrl$max_output_entries || 128 > ctrl$max_matrix_entries || 4096 > ctrl$max_work_ops) {
      cdid_ipw_owner_abort("Empty pair-position query exceeds explicit scalar/header bounds.", "cdid_ipw_owner_bound_error")
    }
    pair_desc <- ipw_owner_make_pair_descriptor(panel_owner, s, t)
    return(data.frame(
      a_index = integer(0L),
      unit_index = integer(0L),
      id = panel_owner$a_id[integer(0L)],
      g = numeric(0L),
      s = numeric(0L),
      t = numeric(0L),
      obs_s = integer(0L),
      obs_t = integer(0L),
      stringsAsFactors = FALSE
    ))
  }
  s_idx <- ipw_owner_int_bound(a_start, "a_start", min_val = 1)
  e_idx <- if (is.null(a_end)) as.double(n_a) else ipw_owner_int_bound(a_end, "a_end", min_val = 1)
  if (s_idx > e_idx || e_idx > n_a) {
    cdid_ipw_owner_abort("Invalid `[a_start, a_end]` bounds.", "cdid_ipw_owner_bound_error")
  }
  b_len <- as.integer(e_idx - s_idx + 1)
  if (b_len > ctrl$max_chunk_rows || 8 * b_len + 128 > ctrl$max_output_entries ||
      64 * b_len + 128 > ctrl$max_matrix_entries || 64 * b_len + 4096 > ctrl$max_work_ops) {
    cdid_ipw_owner_abort("Pair position query exceeds explicit `control` bounds.", "cdid_ipw_owner_bound_error")
  }
  pair_desc <- ipw_owner_make_pair_descriptor(panel_owner, s, t)

  hit <- ipw_owner_query_pair_mask(
    panel_owner, pair_desc$s, pair_desc$t,
    a_start = as.integer(s_idx), a_end = as.integer(e_idx),
    ds = pair_desc$s_calendar_index, dt = pair_desc$t_calendar_index
  )
  nh <- length(hit$a_index)
  data.frame(
    a_index = hit$a_index,
    unit_index = hit$unit_index,
    id = panel_owner$a_id[hit$a_index],
    g = panel_owner$a_g[hit$a_index],
    s = rep(pair_desc$s, nh),
    t = rep(pair_desc$t, nh),
    obs_s = hit$obs_s,
    obs_t = hit$obs_t,
    stringsAsFactors = FALSE
  )
}

# Bounded traversal of observed pair positions across A-units (outcome-free by default;
# reads raw stored endpoints only when `read_outcomes = TRUE`, never subtracting globally).
traverse_ipw_panel_pair_positions <- function(
  panel_owner,
  control,
  callback,
  read_outcomes = FALSE,
  trace_env = NULL
) {
  validate_ipw_panel_owner(panel_owner)
  ctrl <- validate_ipw_owner_control(control)
  if (!is.function(callback)) {
    cdid_ipw_owner_abort("`callback` must be a function.", "cdid_ipw_owner_input_error")
  }
  if (!gmm_scalar_logical(read_outcomes)) {
    cdid_ipw_owner_abort("`read_outcomes` must be a scalar logical.", "cdid_ipw_owner_input_error")
  }
  total_pairs <- panel_owner$counts$n_pairs
  size <- min(ctrl$max_chunk_rows, total_pairs)
  if ((if (read_outcomes) 17 else 15) * max(1, size) + 128 > ctrl$max_output_entries ||
      64 * max(1, size) + 256 > ctrl$max_matrix_entries ||
      64 * max(1, total_pairs) + 32 * panel_owner$n_a_units + 256 > ctrl$max_work_ops) {
    cdid_ipw_owner_abort("Pair position traversal exceeds explicit `control` bounds.", "cdid_ipw_owner_bound_error")
  }
  ipw_owner_reset_trace(trace_env) # after traversal preflight
  if (size == 0) {
    return(invisible(list(n_pairs = 0, n_batches = 0L, max_batch_size = 0L, endpoints_read = 0L)))
  }

  bu <- integer(size)
  bs <- bt <- numeric(size)
  used <- 0L
  emitted <- 0
  batches <- 0L
  maximum <- 0L
  endpoints_read <- 0

  emit <- function() {
    ipw_owner_check_interrupt()
    k <- seq_len(used)
    u <- bu[k]
    a <- if (max(bs[k]) <= .Machine$integer.max) as.integer(bs[k]) else bs[k]
    b <- if (max(bt[k]) <= .Machine$integer.max) as.integer(bt[k]) else bt[k]
    obs <- panel_owner$observations
    s_t <- obs$time[a]
    t_t <- obs$time[b]
    gap <- t_t - s_t
    if (any(is.na(gap) | gap <= 0)) {
      cdid_ipw_owner_abort("Invalid pair time order.", "cdid_ipw_owner_time_error")
    }
    io_raw <- b - a - 1
    io <- if (max(io_raw) <= .Machine$integer.max) as.integer(io_raw) else io_raw
    ip <- as.integer(obs$calendar_index[b] - obs$calendar_index[a] - 1L)
    batch <- data.frame(
      pair_index = emitted + as.double(k),
      unit_index = u,
      a_index = panel_owner$unit_to_a[u],
      id = panel_owner$units$id[u],
      g = panel_owner$units$g[u],
      s = s_t,
      t = t_t,
      obs_s = a,
      obs_t = b,
      time_gap = gap,
      time_gap_finite = is.finite(gap),
      n_intervening_obs = io,
      n_intervening_panel_dates = ip,
      is_consecutive_unit = io == 0L,
      is_consecutive_panel = ip == 0L,
      stringsAsFactors = FALSE
    )
    if (isTRUE(read_outcomes)) {
      ep <- ipw_owner_read_y_positions(panel_owner, a, b, trace_env = trace_env)
      batch$y_s <- ep$y_left_obs
      batch$y_t <- ep$y_right_obs
      endpoints_read <<- endpoints_read + 2 * as.double(used)
    }
    callback(batch)
    emitted <<- emitted + used
    batches <<- batches + 1L
    maximum <<- max(maximum, used)
    used <<- 0L
  }

  for (u in panel_owner$a_unit_index) {
    begin <- panel_owner$units$obs_start[u]
    finish <- panel_owner$units$obs_end[u]
    for (a in seq.int(begin, finish - 1)) {
      b <- a + 1
      while (b <= finish) {
        take <- min(size - used, finish - b + 1)
        k <- seq.int(used + 1, used + take)
        bu[k] <- u
        bs[k] <- a
        bt[k] <- seq.int(b, length.out = take)
        used <- used + take
        b <- b + take
        if (used == size) emit()
      }
    }
  }
  if (used > 0L) emit()
  if (emitted != total_pairs) {
    cdid_ipw_owner_abort("Pair position traversal count mismatch.", "cdid_ipw_owner_count_error")
  }
  invisible(list(
    n_pairs = emitted,
    n_batches = logit_cert_count_output(batches),
    max_batch_size = logit_cert_count_output(maximum),
    endpoints_read = logit_cert_count_output(endpoints_read)
  ))
}
