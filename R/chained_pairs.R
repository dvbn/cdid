# Private CP02 pair index. No estimator or selection model is implemented here.
# Retained storage is O(input rows + units + calendar); pairs are visited lazily.

pair_scalar_vector <- function(x, n) {
  is.atomic(x) && is.null(dim(x)) && length(x) == n
}

pair_count_safe <- function(n_obs) {
  # Halve an even factor BEFORE multiplying, in double arithmetic. Counts and
  # indices must remain exact integers in IEEE doubles (at most 2^53 - 1).
  if (!is.numeric(n_obs) || !is.null(dim(n_obs)) || anyNA(n_obs) ||
      any(!is.finite(n_obs) | n_obs < 0 | n_obs != floor(n_obs))) {
    cdid_abort("Invalid observation counts.", "cdid_pairs_count_error")
  }
  n <- as.double(n_obs)
  a <- ifelse(n %% 2 == 0, n / 2, n)
  b <- ifelse(n %% 2 == 0, pmax(n - 1, 0), pmax(n - 1, 0) / 2)
  limit <- 2^53 - 1
  if (any(a > limit / pmax(b, 1))) {
    cdid_abort("Pair count exceeds exact integer arithmetic.", "cdid_pairs_count_error")
  }
  ans <- as.double(a * b)
  total <- 0
  for (v in ans) {
    if (v > limit - total) {
      cdid_abort("Total pair count exceeds exact integer arithmetic.", "cdid_pairs_count_error")
    }
    total <- total + v
  }
  list(by_unit = ans, total = total)
}

pair_bound <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || is.object(x) || !is.null(dim(x)) || length(x) != 1L ||
      is.na(x) || !is.finite(x) || x != floor(x) ||
      x < (if (allow_zero) 0 else 1) || x > 2^53 - 1) {
    cdid_abort(paste0("`", name, "` must be an explicit finite integer bound."),
               "cdid_pairs_bound_error")
  }
  as.double(x)
}

construct_observed_pairs <- function(data, idname, tname, yname, gname) {
  if (!is.data.frame(data)) cdid_abort("`data` must be a data.frame.", "cdid_pairs_data_error")
  cols <- c(idname = idname, tname = tname, yname = yname, gname = gname)
  args <- list(idname, tname, yname, gname)
  if (any(!vapply(args, function(x) is.character(x) && is.null(dim(x)) &&
                 length(x) == 1L && !is.na(x) && nzchar(x), logical(1L))) ||
      length(unique(cols)) != 4L ||
      any(vapply(cols, function(x) sum(names(data) == x, na.rm = TRUE) != 1L, logical(1L)))) {
    cdid_abort("Select four distinct, unambiguous existing column names.", "cdid_pairs_column_error")
  }
  n <- nrow(data)
  id <- data[[idname]]; time <- data[[tname]]
  y <- data[[yname]]; g <- data[[gname]]
  if (!pair_scalar_vector(id, n)) cdid_abort("ID must be a scalar column.", "cdid_pairs_id_error")
  if (is.factor(id)) id <- as.character(id)
  if (!(is.numeric(id) || is.character(id)) || is.object(id) || anyNA(id) ||
      (is.numeric(id) && any(!is.finite(id))) ||
      (is.character(id) && any(!nzchar(id)))) {
    cdid_abort("ID must be a nonmissing finite numeric or nonempty character vector.", "cdid_pairs_id_error")
  }
  if (!pair_scalar_vector(time, n) || !is.numeric(time) || is.object(time) ||
      anyNA(time) || any(!is.finite(time))) {
    cdid_abort("Time must be a finite numeric scalar column.", "cdid_pairs_time_error")
  }
  if (!pair_scalar_vector(g, n) || !is.numeric(g) || is.object(g) || anyNA(g) || any(g == -Inf)) {
    cdid_abort("Cohort must be a numeric scalar column; only +Inf is allowed.", "cdid_pairs_cohort_error")
  }
  if (!pair_scalar_vector(y, n)) cdid_abort("Outcome must be a scalar column.", "cdid_pairs_outcome_error")
  if (is.logical(y) && all(is.na(y))) y <- as.numeric(y)
  if (!is.numeric(y) || is.object(y) || any(is.nan(y) | is.infinite(y))) {
    cdid_abort("Outcome must be numeric; only NA is unobserved.", "cdid_pairs_outcome_error")
  }
  # Numeric values are matched directly, never through printed identifiers.
  ids <- sort(unique(id), method = "radix")
  nu <- length(ids)
  ui <- match(id, ids)
  ord <- order(ui, time, method = "radix")
  ui <- ui[ord]; time <- as.double(time[ord]); y <- as.double(y[ord]); g <- as.double(g[ord])
  if (n > 1L) {
    same <- ui[-1L] == ui[-n]
    if (any(same & g[-1L] != g[-n])) cdid_abort("Within-unit cohort changes are invalid.",
      c("cdid_pairs_cohort_change_error", "cdid_pairs_cohort_error"))
    if (any(same & time[-1L] == time[-n])) cdid_abort("Duplicate (ID, time) rows are invalid.", "cdid_pairs_duplicate_error")
  }
  all_times <- sort(unique(time))
  observed <- !is.na(y)
  ou <- ui[observed]; ot <- time[observed]; oy <- y[observed]
  nr <- tabulate(ui, nbins = nu); no <- tabulate(ou, nbins = nu)
  counts <- pair_count_safe(no)
  end <- cumsum(as.double(no)); start <- end - no + 1
  ug <- g[match(seq_len(nu), ui)]
  first <- last <- sy <- singleton_time <- rep(NA_real_, nu)
  nonzero <- which(no > 0L); single <- which(no == 1L)
  first[nonzero] <- ot[start[nonzero]]; last[nonzero] <- ot[end[nonzero]]
  sy[single] <- oy[start[single]]; singleton_time[single] <- first[single]
  # The greatest chronological difference is attained at an endpoint against
  # a previous extremum. This linear scan rejects ANY overflowing outcome pair
  # before a callback can consume even an earlier batch. Finite numeric dates
  # whose elapsed gap overflows double precision remain valid pairs.
  for (u in which(no >= 2L)) {
    ys <- oy[seq.int(start[u], end[u])]
    lo <- cummin(ys); hi <- cummax(ys)
    if (any(!is.finite(ys[-1L] - lo[-length(lo)])) ||
        any(!is.finite(ys[-1L] - hi[-length(hi)]))) {
      cdid_abort("A finite outcome pair has a nonfinite difference.", "cdid_pairs_outcome_error")
    }
  }
  units <- data.frame(unit_index = seq_len(nu), id = ids, g = ug,
    n_rows = nr, n_obs = no, n_pairs = counts$by_unit,
    is_unobserved = no == 0L, is_singleton = no == 1L, pair_eligible = no >= 2L,
    first_obs_time = first, last_obs_time = last,
    singleton_time = singleton_time, singleton_y = sy,
    obs_start = start, obs_end = end, stringsAsFactors = FALSE)
  singletons <- data.frame(unit_index = single, id = ids[single], g = ug[single],
    time = first[single], y = sy[single], n_rows = nr[single], stringsAsFactors = FALSE)
  structure(list(units = units, singletons = singletons,
    observations = data.frame(time = ot, y = oy, calendar_index = match(ot, all_times)),
    times = sort(unique(ot)), all_times = all_times, cohorts = sort(unique(ug)),
    counts = list(n_input_rows = n, n_observed_rows = length(ot), n_na_rows = sum(!observed),
      n_units = nu, n_unobserved_units = sum(no == 0L), n_singleton_units = length(single),
      n_pair_eligible_units = sum(no >= 2L), n_pairs = counts$total), cols = cols),
    class = "cdid_observed_pair_index")
}

empty_pair_batch <- function(index, n = 0L) {
  data.frame(pair_index = numeric(n), unit_index = integer(n),
    id = index$units$id[rep(NA_integer_, n)], g = numeric(n),
    s = numeric(n), t = numeric(n), y_s = numeric(n), y_t = numeric(n),
    dy = numeric(n), time_gap = numeric(n), time_gap_finite = logical(n),
    n_intervening_obs = integer(n),
    n_intervening_panel_dates = integer(n), is_consecutive_unit = logical(n),
    is_consecutive_panel = logical(n), stringsAsFactors = FALSE)
}

traverse_observed_pairs <- function(index, chunk_size, callback) {
  if (!inherits(index, "cdid_observed_pair_index")) cdid_abort("Expected compact pair index.", "cdid_pairs_data_error")
  size <- pair_bound(chunk_size, "chunk_size")
  if (!is.function(callback)) cdid_abort("`callback` must be a function.", "cdid_pairs_data_error")
  size <- min(size, index$counts$n_pairs)
  if (size == 0) return(invisible(list(n_pairs = 0, n_batches = 0, max_batch_size = 0)))
  # Fixed buffers. The nested (unit, start, end) traversal also splits a single
  # long unit, rather than constructing combn(m, 2) or a per-unit pair vector.
  bu <- integer(size); bs <- bt <- numeric(size)
  used <- 0L; emitted <- batches <- maximum <- 0
  emit <- function() {
    k <- seq_len(used); u <- bu[k]; a <- bs[k]; b <- bt[k]
    obs <- index$observations
    ys <- obs$y[a]; yt <- obs$y[b]; dy <- yt - ys
    gap <- obs$time[b] - obs$time[a]
    if (any(!is.finite(dy))) cdid_abort("Pair difference overflow.", "cdid_pairs_outcome_error")
    if (any(is.na(gap) | gap <= 0)) cdid_abort("Invalid pair time order.", "cdid_pairs_time_error")
    io <- as.integer(b - a - 1); ip <- as.integer(obs$calendar_index[b] - obs$calendar_index[a] - 1L)
    batch <- data.frame(pair_index = emitted + as.double(k), unit_index = u,
      id = index$units$id[u], g = index$units$g[u], s = obs$time[a], t = obs$time[b],
      y_s = ys, y_t = yt, dy = dy, time_gap = gap, time_gap_finite = is.finite(gap),
      n_intervening_obs = io, n_intervening_panel_dates = ip,
      is_consecutive_unit = io == 0L, is_consecutive_panel = ip == 0L,
      stringsAsFactors = FALSE)
    callback(batch) # Result deliberately discarded; errors propagate unchanged.
    emitted <<- emitted + used; batches <<- batches + 1
    maximum <<- max(maximum, used); used <<- 0L
  }
  for (u in which(index$units$pair_eligible)) {
    begin <- index$units$obs_start[u]; finish <- index$units$obs_end[u]
    for (a in seq.int(begin, finish - 1)) {
      b <- a + 1
      while (b <= finish) {
        take <- min(size - used, finish - b + 1)
        k <- seq.int(used + 1, used + take)
        bu[k] <- u; bs[k] <- a; bt[k] <- seq.int(b, length.out = take)
        used <- used + take; b <- b + take
        if (used == size) emit()
      }
    }
  }
  if (used > 0L) emit()
  if (emitted != index$counts$n_pairs) cdid_abort("Pair traversal count mismatch.", "cdid_pairs_count_error")
  invisible(list(n_pairs = emitted, n_batches = batches, max_batch_size = maximum))
}

materialize_observed_pairs <- function(index, max_pairs, chunk_size) {
  bound <- pair_bound(max_pairs, "max_pairs", allow_zero = TRUE)
  if (!inherits(index, "cdid_observed_pair_index")) cdid_abort("Expected compact pair index.", "cdid_pairs_data_error")
  p <- index$counts$n_pairs
  if (p > bound) cdid_abort("Pair materialization exceeds explicit preallocation bound.", "cdid_pairs_bound_error")
  pairs <- empty_pair_batch(index, p)
  traverse_observed_pairs(index, chunk_size, function(batch) {
    k <- batch$pair_index
    for (nm in names(pairs)) pairs[[nm]][k] <<- batch[[nm]]
  })
  # Exact integer-index keys for cohort/date cells; no floating-point formatting.
  cells <- data.frame(g = numeric(), s = numeric(), t = numeric(), n_pairs = numeric())
  if (p > 0) {
    ord <- order(pairs$g, pairs$s, pairs$t, method = "radix")
    gg <- pairs$g[ord]; ss <- pairs$s[ord]; tt <- pairs$t[ord]
    change <- c(TRUE, gg[-1L] != gg[-p] | ss[-1L] != ss[-p] | tt[-1L] != tt[-p])
    k <- which(change)
    cells <- data.frame(g = gg[k], s = ss[k], t = tt[k], n_pairs = as.double(diff(c(k, p + 1))))
  }
  out <- index; out$pairs <- pairs; out$cohort_pair_cells <- cells
  class(out) <- "cdid_observed_pairs"
  out
}

chained_pairs <- construct_observed_pairs
