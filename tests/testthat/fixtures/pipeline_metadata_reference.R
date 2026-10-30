# Independent brute-force reference implementations and test fixtures for the
# CP02 fitted-pipeline metadata slice (`R/chained_pipeline_metadata.R`).
#
# Zero calls to `cdid` functions are permitted in this reference file.

ref_exact_cmp_sum <- function(a, b, c) {
  a <- as.double(a)
  b <- as.double(b)
  c <- as.double(c)
  if (is.infinite(c) && c > 0) return(-1L)
  s <- a + b
  if (s < c) return(-1L)
  if (s > c) return(1L)
  if (abs(a) >= abs(b)) {
    bb <- s - a
    err <- (a - (s - bb)) + (b - bb)
  } else {
    aa <- s - b
    err <- (b - (s - aa)) + (a - aa)
  }
  if (err < 0) -1L else if (err > 0) 1L else 0L
}

ref_match_exact_double <- function(x, table) {
  # Exact binary64 value comparison (normalizing -0.0 to +0.0), never string formatting
  x_d <- as.double(x)
  x_d[x_d == 0.0] <- 0.0
  t_d <- as.double(table)
  t_d[t_d == 0.0] <- 0.0
  vapply(x_d, function(val) {
    hits <- which(t_d == val)
    if (length(hits) == 1L) hits else NA_integer_
  }, integer(1L))
}

brute_force_panel_occupancy <- function(
  data,
  idname = "id",
  tname = "time",
  yname = "y",
  gname = "g",
  weights_by_unit = NULL
) {
  id_raw <- data[[idname]]
  t_raw <- as.double(data[[tname]])
  t_raw[t_raw == 0.0] <- 0.0
  y_raw <- if (is.logical(data[[yname]])) rep(NA_real_, nrow(data)) else as.double(data[[yname]])
  g_raw <- as.double(data[[gname]])
  g_raw[g_raw == 0.0] <- 0.0

  id_is_num <- is.numeric(id_raw)
  id_key <- if (id_is_num) {
    d <- as.double(id_raw)
    d[d == 0.0] <- 0.0
    d
  } else {
    as.character(id_raw)
  }

  ids_key <- sort(unique(id_key), method = "radix")
  ids <- if (identical(typeof(id_raw), "integer")) as.integer(ids_key) else ids_key
  nu <- length(ids_key)

  ui <- if (id_is_num) ref_match_exact_double(id_key, ids_key) else match(id_key, ids_key)
  ord <- order(ui, t_raw, method = "radix")
  ui <- ui[ord]
  t_ord <- t_raw[ord]
  y_ord <- y_raw[ord]
  g_ord <- g_raw[ord]

  calendar <- sort(unique(t_ord), method = "radix")
  n_cal <- length(calendar)

  obs_mask <- !is.na(y_ord)
  ou <- ui[obs_mask]
  ot <- t_ord[obs_mask]
  oy <- y_ord[obs_mask]
  ocal <- ref_match_exact_double(ot, calendar)

  nr <- tabulate(ui, nbins = nu)
  no <- tabulate(ou, nbins = nu)
  n_na <- nr - no

  ug <- vapply(seq_len(nu), function(u) g_ord[which(ui == u)[1L]], numeric(1L))
  cohorts <- sort(unique(ug), method = "radix")
  n_coh <- length(cohorts)
  u_coh_idx <- ref_match_exact_double(ug, cohorts)

  a_units <- which(no >= 2L)
  single_units <- which(no == 1L)
  unobs_units <- which(no == 0L)
  a_cohorts <- if (length(a_units) > 0L) sort(unique(ug[a_units]), method = "radix") else numeric(0L)

  # Brute-force 3D count array over (cohort_index, s_index, t_index)
  count_arr <- array(0, dim = c(n_coh, n_cal, n_cal))
  pair_rows <- list()
  pair_counter <- 0L

  for (u in a_units) {
    obs_pos <- which(ou == u)
    m_u <- length(obs_pos)
    ki <- u_coh_idx[u]
    for (i in seq_len(m_u - 1L)) {
      for (j in seq.int(i + 1L, m_u)) {
        pos_s <- obs_pos[i]
        pos_t <- obs_pos[j]
        si <- ocal[pos_s]
        ti <- ocal[pos_t]
        count_arr[ki, si, ti] <- count_arr[ki, si, ti] + 1
        pair_counter <- pair_counter + 1L
        pair_rows[[pair_counter]] <- data.frame(
          pair_index = as.double(pair_counter),
          unit_index = u,
          cohort_index = ki,
          cohort = cohorts[ki],
          s_index = si,
          t_index = ti,
          s = calendar[si],
          t = calendar[ti],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Extract sparse (cohort, s, t, count) ordered by (cohort_index, s_index, t_index)
  sparse_list <- list()
  idx_list <- list()
  dense_list <- list()
  cell_id <- 0L
  dense_id <- 0L

  for (ki in seq_len(n_coh)) {
    for (si in seq_len(max(0L, n_cal - 1L))) {
      for (ti in seq.int(si + 1L, n_cal)) {
        cnt <- count_arr[ki, si, ti]
        dense_id <- dense_id + 1L
        dense_list[[dense_id]] <- data.frame(
          cohort = cohorts[ki],
          s = calendar[si],
          t = calendar[ti],
          count = as.double(cnt),
          stringsAsFactors = FALSE
        )
        if (cnt > 0) {
          cell_id <- cell_id + 1L
          sparse_list[[cell_id]] <- data.frame(
            cohort = cohorts[ki],
            s = calendar[si],
            t = calendar[ti],
            count = as.double(cnt),
            stringsAsFactors = FALSE
          )
          idx_list[[cell_id]] <- data.frame(
            cell_id = cell_id,
            cohort_index = ki,
            s_index = si,
            t_index = ti,
            count = as.double(cnt),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  sparse_occupancy <- if (cell_id > 0L) {
    do.call(rbind, sparse_list)
  } else {
    data.frame(cohort = numeric(0L), s = numeric(0L), t = numeric(0L), count = numeric(0L), stringsAsFactors = FALSE)
  }
  rownames(sparse_occupancy) <- NULL

  cell_indices <- if (cell_id > 0L) {
    do.call(rbind, idx_list)
  } else {
    data.frame(
      cell_id = integer(0L),
      cohort_index = integer(0L),
      s_index = integer(0L),
      t_index = integer(0L),
      count = numeric(0L),
      stringsAsFactors = FALSE
    )
  }
  rownames(cell_indices) <- NULL

  dense_occupancy <- if (dense_id > 0L) {
    do.call(rbind, dense_list)
  } else {
    data.frame(cohort = numeric(0L), s = numeric(0L), t = numeric(0L), count = numeric(0L), stringsAsFactors = FALSE)
  }
  rownames(dense_occupancy) <- NULL

  coh_a_count <- vapply(seq_len(n_coh), function(ki) sum(u_coh_idx[a_units] == ki), integer(1L))
  coh_a_weight <- if (!is.null(weights_by_unit)) {
    vapply(seq_len(n_coh), function(ki) {
      u_in_k <- a_units[u_coh_idx[a_units] == ki]
      sum(weights_by_unit[u_in_k])
    }, numeric(1L))
  } else {
    rep(NA_real_, n_coh)
  }

  list(
    calendar = calendar,
    observed_times = sort(unique(ot), method = "radix"),
    cohorts = cohorts,
    a_cohorts = a_cohorts,
    n_input_rows = nrow(data),
    n_observed_rows = length(ot),
    n_a_observed_rows = sum(no[a_units]),
    n_explicit_na_rows = sum(n_na),
    n_total_units = nu,
    n_a_units = length(a_units),
    n_singleton_units = length(single_units),
    n_unobserved_units = length(unobs_units),
    n_pairs = as.double(pair_counter),
    count_arr = count_arr,
    sparse_occupancy = sparse_occupancy,
    dense_occupancy = dense_occupancy,
    cell_indices = cell_indices,
    coh_a_count = coh_a_count,
    coh_a_weight = coh_a_weight
  )
}

brute_force_catalog_mappings <- function(ref_occ, catalog) {
  rows <- catalog$rows
  ops <- catalog$operators
  pl <- catalog$pools
  n_rows <- nrow(rows)
  n_ops <- nrow(ops)
  n_pools <- pl$n_pools

  pool_a_count <- vapply(seq_len(n_pools), function(p) {
    p_idx <- pl$entries[seq.int(pl$start[p], length.out = pl$size[p])]
    as.double(sum(ref_occ$coh_a_count[p_idx]))
  }, numeric(1L))

  target_occ_cell_id <- integer(n_rows)
  target_pair_count <- numeric(n_rows)
  control_pool_pair_count <- numeric(n_rows)
  risk_pair_count <- numeric(n_rows)
  all_cohort_pair_count <- numeric(n_rows)
  n_supported_controls <- integer(n_rows)

  ci <- ref_occ$cell_indices
  for (r in seq_len(n_rows)) {
    gi <- rows$target_index[r]
    si <- rows$s_index[r]
    ti <- rows$t_index[r]
    pid <- rows$pool_id[r]
    pool_k <- pl$entries[seq.int(pl$start[pid], length.out = pl$size[pid])]

    hit_cell <- which(ci$cohort_index == gi & ci$s_index == si & ci$t_index == ti)
    target_occ_cell_id[r] <- ci$cell_id[hit_cell]
    target_pair_count[r] <- ref_occ$count_arr[gi, si, ti]
    control_pool_pair_count[r] <- sum(ref_occ$count_arr[pool_k, si, ti])
    risk_pair_count[r] <- target_pair_count[r] + control_pool_pair_count[r]
    all_cohort_pair_count[r] <- sum(ref_occ$count_arr[, si, ti])
    n_supported_controls[r] <- length(pool_k)
  }

  target_a_unit_count <- as.double(ref_occ$coh_a_count[rows$target_index])
  control_pool_a_unit_count <- if (n_rows > 0L) pool_a_count[rows$pool_id] else numeric(0L)
  risk_a_unit_count <- target_a_unit_count + control_pool_a_unit_count

  list(
    target_occupancy_cell_id = target_occ_cell_id,
    target_pair_count = target_pair_count,
    control_pool_pair_count = control_pool_pair_count,
    risk_pair_count = risk_pair_count,
    all_cohort_pair_count = all_cohort_pair_count,
    target_a_unit_count = target_a_unit_count,
    control_pool_a_unit_count = control_pool_a_unit_count,
    risk_a_unit_count = risk_a_unit_count,
    n_supported_controls = n_supported_controls
  )
}

make_pipeline_metadata_fixtures <- function() {
  # 1. Rotating / gapped staggered panel with D10 + D13 fixed-horizon opportunities
  rotating_gapped <- rbind(
    # Cohort g = 3: observed at (1,2), (2,3), (2,4), (3,4)
    data.frame(id = 1L, time = c(1, 2, 3, 4), y = c(1.0, 2.0, NA, 4.0), g = 3, x = 0.5),
    data.frame(id = 2L, time = c(1, 2, 3, 4), y = c(NA, 1.5, 2.5, 3.5), g = 3, x = -0.5),
    # Cohort g = 4: observed at (1,2), (2,3), (2,4), (3,4)
    data.frame(id = 3L, time = c(1, 2, 3, 4), y = c(0.5, 1.5, 2.0, NA), g = 4, x = 1.0),
    data.frame(id = 4L, time = c(1, 2, 3, 4), y = c(1.2, 2.2, NA, 3.2), g = 4, x = -1.0),
    # Never-treated g = Inf: observed at (1,2), (2,3), (2,4), (3,4)
    data.frame(id = 5L, time = c(1, 2, 3, 4), y = c(2.0, 3.0, 4.0, 5.0), g = Inf, x = 0.25),
    data.frame(id = 6L, time = c(1, 2, 3, 4), y = c(NA, 2.5, 3.5, 4.5), g = Inf, x = -0.25)
  )

  # 2. Balanced panel
  balanced <- expand.grid(time = c(1, 2, 3, 4), id = 1:6, KEEP.OUT.ATTRS = FALSE)
  balanced <- balanced[, c("id", "time")]
  balanced$g <- c(3, 3, 4, 4, Inf, Inf)[balanced$id]
  balanced$y <- 1.0 + 0.1 * balanced$id + 0.25 * balanced$time
  balanced$x <- (-1)^balanced$id * 0.5

  # 3. Singleton and unobserved units mixed with A-units (including singleton-only cohort 5 and date 0.5)
  with_singletons <- rbind(
    rotating_gapped,
    # Singleton unit in cohort 5 at time 0.5 (extra calendar date & extra cohort with zero pairs)
    data.frame(id = 7L, time = c(0.5, 2, 3), y = c(9.0, NA, NA), g = 5, x = 0.1),
    # Unobserved unit in cohort 6 (all y = NA)
    data.frame(id = 8L, time = c(1, 2), y = c(NA_real_, NA_real_), g = 6, x = -0.1)
  )

  # 4. Empty-A panel (only singletons and unobserved units; n_a_units == 0, n_pairs == 0)
  empty_a <- rbind(
    data.frame(id = 10L, time = c(1, 2, 3), y = c(1.0, NA, NA), g = 3, x = 0.0),
    data.frame(id = 20L, time = c(1, 2, 3), y = c(NA, 2.0, NA), g = Inf, x = 1.0),
    data.frame(id = 30L, time = c(1, 2, 3), y = c(NA_real_, NA_real_, NA_real_), g = 4, x = -1.0)
  )

  # 5. Exact adjacent binary64 dates and cohorts (distinct in IEEE-754 double, identical in default format/as.character)
  t1 <- 1.0
  t1_adj <- 1.0 + .Machine$double.eps
  t2 <- 2.0
  t3 <- 3.0
  t3_adj <- 3.0 + 2 * .Machine$double.eps
  g3 <- 3.0
  g3_adj <- 3.0 + 2 * .Machine$double.eps

  adjacent_binary64 <- rbind(
    # Cohort g3: observed at t1, t2, t3
    data.frame(id = 1L, time = c(t1, t2, t3), y = c(1.0, 2.0, 3.0), g = g3, x = 0.5),
    data.frame(id = 2L, time = c(t1, t1_adj, t2, t3_adj), y = c(1.1, 1.2, 2.1, 3.1), g = g3, x = -0.5),
    # Adjacent cohort g3_adj: observed at t1_adj, t2, t3_adj
    data.frame(id = 3L, time = c(t1_adj, t2, t3_adj), y = c(0.5, 1.5, 2.5), g = g3_adj, x = 0.2),
    # Never-treated Inf: observed at all 5 dates
    data.frame(id = 4L, time = c(t1, t1_adj, t2, t3, t3_adj), y = c(2.0, 2.1, 3.0, 4.0, 4.1), g = Inf, x = -0.2),
    data.frame(id = 5L, time = c(t1, t1_adj, t2, t3, t3_adj), y = c(2.5, 2.6, 3.5, 4.5, 4.6), g = Inf, x = 0.8)
  )

  # 6. Character IDs (including multi-byte UTF-8 and varying byte lengths)
  char_ids <- rotating_gapped
  char_map <- c(
    "unit_alpha_01",
    "unit_\u00e9\u03b2_02",
    "u3",
    "unit_\u6771\u4eac_04_long_identifier",
    "nt_control_05",
    "nt_\u03c9_06"
  )
  char_ids$id <- char_map[char_ids$id]

  list(
    rotating_gapped = rotating_gapped,
    balanced = balanced,
    with_singletons = with_singletons,
    empty_a = empty_a,
    adjacent_binary64 = adjacent_binary64,
    char_ids = char_ids
  )
}
