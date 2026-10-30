require_cdid_namespace <- function() {
  if (!isNamespaceLoaded("cdid") && !requireNamespace("cdid", quietly = TRUE)) {
    stop("Actual installed 'cdid' namespace is required; no fallback is permitted in packaged tests.", call. = FALSE)
  }
  ns <- asNamespace("cdid")
  if (!isNamespace(ns) || !environmentIsLocked(ns)) {
    stop("'cdid' namespace must be a locked R namespace.", call. = FALSE)
  }
  ns
}

# Require installed namespace immediately at file load (never skip)
require_cdid_namespace()

# ---------------------------------------------------------------------------
# Independent Cartesian / Hand Reference Implementations (Zero cdid calls)
# ---------------------------------------------------------------------------

hand_reference_panel_owner <- function(data, idname, tname, yname, gname) {
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
  ui <- match(id_key, ids_key)
  ord <- order(ui, t_raw, method = "radix")

  ui <- ui[ord]
  t_ord <- t_raw[ord]
  y_ord <- y_raw[ord]
  g_ord <- g_raw[ord]

  all_times <- sort(unique(t_ord), method = "radix")
  obs_mask <- !is.na(y_ord)
  ou <- ui[obs_mask]
  ot <- t_ord[obs_mask]
  oy <- y_ord[obs_mask]
  ocal <- match(ot, all_times)

  na_u <- ui[!obs_mask]
  na_t <- t_ord[!obs_mask]
  na_cal <- match(na_t, all_times)

  nr <- tabulate(ui, nbins = nu)
  no <- tabulate(ou, nbins = nu)
  n_na <- nr - no
  by_unit_pairs <- ifelse(no >= 2L, as.double(no) * (as.double(no) - 1) / 2, 0)

  obs_end <- cumsum(as.double(no))
  obs_start <- obs_end - no + 1
  na_end <- cumsum(as.double(n_na))
  na_start <- na_end - n_na + 1

  ug <- vapply(seq_len(nu), function(u) g_ord[which(ui == u)[1L]], numeric(1L))
  a_unit_index <- which(no >= 2L)
  n_a <- length(a_unit_index)
  unit_to_a <- integer(nu)
  if (n_a > 0L) unit_to_a[a_unit_index] <- seq_len(n_a)

  # Cartesian enumeration of all observed within-unit pairs (a, b) with a < b
  cartesian_pairs <- list()
  p_idx <- 0
  for (u in a_unit_index) {
    idx_u <- which(ou == u)
    for (i in seq_len(length(idx_u) - 1L)) {
      for (j in seq.int(i + 1L, length(idx_u))) {
        a <- idx_u[i]
        b <- idx_u[j]
        p_idx <- p_idx + 1
        cartesian_pairs[[p_idx]] <- data.frame(
          pair_index = as.double(p_idx),
          unit_index = u,
          a_index = unit_to_a[u],
          id = ids[u],
          g = ug[u],
          s = ot[a],
          t = ot[b],
          obs_s = a,
          obs_t = b,
          y_s = oy[a],
          y_t = oy[b],
          time_gap = ot[b] - ot[a],
          time_gap_finite = is.finite(ot[b] - ot[a]),
          n_intervening_obs = as.integer(b - a - 1L),
          n_intervening_panel_dates = as.integer(ocal[b] - ocal[a] - 1L),
          is_consecutive_unit = (b - a - 1L) == 0L,
          is_consecutive_panel = (ocal[b] - ocal[a] - 1L) == 0L,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  pairs_df <- if (length(cartesian_pairs) > 0L) do.call(rbind, cartesian_pairs) else data.frame()

  list(
    ids = ids,
    ids_key = ids_key,
    nu = nu,
    all_times = all_times,
    times = sort(unique(ot), method = "radix"),
    cohorts = sort(unique(ug), method = "radix"),
    a_unit_index = a_unit_index,
    unit_to_a = unit_to_a,
    a_id = ids[a_unit_index],
    a_g = ug[a_unit_index],
    a_cohorts = if (n_a > 0L) sort(unique(ug[a_unit_index]), method = "radix") else numeric(0L),
    n_rows = nr,
    n_obs = no,
    n_na_rows = n_na,
    n_pairs = by_unit_pairs,
    obs_start = obs_start,
    obs_end = obs_end,
    na_start = na_start,
    na_end = na_end,
    observations = data.frame(
      time = ot,
      y = oy,
      calendar_index = ocal,
      unit_index = ou,
      a_index = unit_to_a[ou],
      stringsAsFactors = FALSE
    ),
    explicit_na_rows = data.frame(
      unit_index = na_u,
      time = na_t,
      calendar_index = na_cal,
      stringsAsFactors = FALSE
    ),
    cartesian_pairs = pairs_df
  )
}

hand_reference_risk_first_chunk <- function(ref_panel, Z_a, w_a, g, K, s, t, a_start, a_end, read_outcomes = FALSE) {
  cohorts <- c(as.double(g), as.double(K))
  a_seq <- seq.int(as.integer(a_start), as.integer(a_end))
  b_len <- length(a_seq)

  u_idx <- ref_panel$a_unit_index[a_seq]
  u_id <- ref_panel$a_id[a_seq]
  u_g <- ref_panel$a_g[a_seq]
  u_w <- w_a[a_seq]
  cohort_index <- match(u_g, cohorts)
  in_risk <- !is.na(cohort_index)
  is_target <- in_risk & (cohort_index == 1L)

  obs <- ref_panel$observations
  local_pos <- integer(0L)
  obs_a <- integer(0L)
  obs_u <- integer(0L)
  obs_ci <- integer(0L)
  obs_s <- integer(0L)
  obs_t <- integer(0L)

  # Risk-first check: only inspect observation endpoints on risk-set units, or intersect
  # observation positions with risk membership BEFORE reading y
  for (k in seq_len(b_len)) {
    if (!in_risk[k]) next
    ai <- a_seq[k]
    ui <- u_idx[k]
    rows_u <- which(obs$unit_index == ui)
    hit_s <- rows_u[obs$time[rows_u] == as.double(s)]
    hit_t <- rows_u[obs$time[rows_u] == as.double(t)]
    if (length(hit_s) == 1L && length(hit_t) == 1L) {
      local_pos <- c(local_pos, k)
      obs_a <- c(obs_a, ai)
      obs_u <- c(obs_u, ui)
      obs_ci <- c(obs_ci, cohort_index[k])
      obs_s <- c(obs_s, hit_s)
      obs_t <- c(obs_t, hit_t)
    }
  }

  in_observed <- logical(b_len)
  if (length(local_pos) > 0L) in_observed[local_pos] <- TRUE
  is_risk_nonobserver <- in_risk & !in_observed
  chunk_to_obs_pos <- integer(b_len)
  if (length(local_pos) > 0L) chunk_to_obs_pos[local_pos] <- seq_along(local_pos)

  y_left_obs <- if (isTRUE(read_outcomes)) obs$y[obs_s] else NULL
  y_right_obs <- if (isTRUE(read_outcomes)) obs$y[obs_t] else NULL

  list(
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
    local_pos = local_pos,
    obs_a_index = obs_a,
    obs_unit_index = obs_u,
    obs_cohort_index = obs_ci,
    obs_s_index = obs_s,
    obs_t_index = obs_t,
    y_left_obs = y_left_obs,
    y_right_obs = y_right_obs
  )
}

default_owner_test_control <- function(ns, max_chunk_rows = 4L) {
  ns$ipw_owner_control(
    max_chunk_rows = max_chunk_rows,
    max_retained_entries = 100000L,
    max_matrix_entries = 100000L,
    max_work_ops = 500000L,
    max_output_entries = 100000L
  )
}

# ---------------------------------------------------------------------------
# Test 1: Namespace bindings and unexported status
# ---------------------------------------------------------------------------

test_that("CP02 common panel/frame owner and reader entrypoints exist in locked cdid namespace and are unexported", {
  ns <- require_cdid_namespace()
  expect_true(isNamespace(ns))
  expect_true(environmentIsLocked(ns))

  priv_fns <- c(
    "ipw_owner_control",
    "validate_ipw_owner_control",
    "ipw_panel_owner_control",
    "ipw_frame_owner_control",
    "ipw_view_control",
    "ipw_reader_control",
    "preflight_ipw_panel_owner",
    "preflight_ipw_frame_owner",
    "preflight_ipw_nuisance_view",
    "preflight_ipw_link_view",
    "preflight_ipw_range_query",
    "preflight_ipw_stream",
    "construct_ipw_panel_owner",
    "validate_ipw_panel_owner",
    "query_ipw_panel_cell_status",
    "construct_ipw_frame_owner",
    "validate_ipw_frame_owner",
    "verify_ipw_frame_compatibility",
    "construct_ipw_nuisance_view",
    "validate_ipw_nuisance_view",
    "verify_ipw_view_mask_compatibility",
    "construct_ipw_link_view",
    "validate_ipw_link_view",
    "verify_ipw_view_dense_masks",
    "query_ipw_view_chunk",
    "stream_ipw_view_chunks",
    "query_ipw_panel_pair_positions",
    "traverse_ipw_panel_pair_positions"
  )
  exports <- getNamespaceExports("cdid")
  for (fn in priv_fns) {
    expect_true(exists(fn, envir = ns, mode = "function", inherits = FALSE), info = fn)
    expect_true(bindingIsLocked(fn, ns), info = fn)
    expect_identical(environment(get(fn, envir = ns, inherits = FALSE)), ns, info = fn)
    expect_false(fn %in% exports, info = fn)
  }
})

# ---------------------------------------------------------------------------
# Test 2: Independent Cartesian/hand reference across balanced, intermittent,
# reentry, gapped, exact-two, off-grid g, never-treated (+Inf), singletons,
# unobserved units, and explicit NA vs missing-row cell status
# ---------------------------------------------------------------------------

test_that("panel and frame owners match independent Cartesian/hand reference across irregular panel topologies", {
  ns <- require_cdid_namespace()
  ctrl <- default_owner_test_control(ns, max_chunk_rows = 3L)

  # Irregular panel combining:
  # - Unit 105: balanced across 4 observed dates (10.0, 20.5, 35.0, 50.25), off-grid cohort g = 25.75
  # - Unit 101: intermittent with explicit NA at 20.5 and missing row at 35.0, reentry at 50.25 (exact-two obs: 10.0, 50.25), g = 25.75
  # - Unit 109: reentry with obs at 10.0, explicit NA at 20.5, obs at 35.0, obs at 50.25, never-treated g = Inf
  # - Unit 102: exact-two obs at 10.0 and 35.0, g = 40.0
  # - Unit 108: singleton (1 obs at 20.5, explicit NA at 35.0), g = 25.75
  # - Unit 104: panel-unobserved (explicit NA at 10.0 and 99.0; note 99.0 is in all_times but NOT in times!), g = Inf
  # - Unit 107: exact-two obs at 10.0 and 35.0, g = Inf
  df <- rbind(
    data.frame(id = c(105, 105, 105, 105), time = c(50.25, 10.0, 35.0, 20.5), y = c(4.5, 1.0, 3.25, 2.0), g = rep(25.75, 4)),
    data.frame(id = c(101, 101, 101), time = c(10.0, 20.5, 50.25), y = c(10.0, NA_real_, 15.5), g = rep(25.75, 3)),
    data.frame(id = c(109, 109, 109, 109), time = c(10.0, 20.5, 35.0, 50.25), y = c(7.0, NA_real_, 8.5, 9.25), g = rep(Inf, 4)),
    data.frame(id = c(102, 102), time = c(35.0, 10.0), y = c(22.0, 20.0), g = rep(40.0, 2)),
    data.frame(id = c(108, 108), time = c(20.5, 35.0), y = c(5.0, NA_real_), g = rep(25.75, 2)),
    data.frame(id = c(104, 104), time = c(10.0, 99.0), y = c(NA_real_, NA_real_), g = rep(Inf, 2)),
    data.frame(id = c(107, 107), time = c(10.0, 35.0), y = c(30.0, 33.5), g = rep(Inf, 2))
  )

  ref <- hand_reference_panel_owner(df, "id", "time", "y", "g")
  panel <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl, owner_id = "p_irreg")
  ns$validate_ipw_panel_owner(panel)

  expect_identical(panel$accounting$global_pair_subtraction_performed, FALSE)
  expect_identical(panel$accounting$endpoints_read_during_ingestion, 0L)
  expect_identical(panel$units$id, ref$ids)
  expect_identical(panel$all_times, ref$all_times)
  expect_identical(panel$times, ref$times)
  expect_true(99.0 %in% panel$all_times)
  expect_false(99.0 %in% panel$times)
  expect_identical(panel$cohorts, ref$cohorts)
  expect_identical(panel$a_unit_index, ref$a_unit_index)
  expect_identical(panel$unit_to_a, ref$unit_to_a)
  expect_identical(panel$a_id, ref$a_id)
  expect_identical(panel$a_g, ref$a_g)
  expect_identical(panel$a_cohorts, ref$a_cohorts)
  expect_identical(panel$units$n_rows, ref$n_rows)
  expect_identical(panel$units$n_obs, ref$n_obs)
  expect_identical(panel$units$n_na_rows, ref$n_na_rows)
  expect_identical(panel$units$n_pairs, ref$n_pairs)
  expect_identical(panel$observations, ref$observations)
  expect_identical(panel$explicit_na_rows, ref$explicit_na_rows)
  expect_identical(panel$singletons$id, 108)
  expect_identical(panel$singletons$y, 5.0)
  expect_identical(panel$unobserved_units$id, 104)
  expect_identical(panel$unobserved_units$n_na_rows, 2L)

  # Verify explicit missing-outcome (`explicit_na`) vs missing-row (`missing_row`) vs `observed`
  cell_st <- ns$query_ipw_panel_cell_status(
    panel,
    id = c(101, 101, 101, 104, 999),
    time = c(10.0, 20.5, 35.0, 99.0, 10.0),
    control = default_owner_test_control(ns, max_chunk_rows = 5L)
  )
  expect_identical(cell_st$status, c("observed", "explicit_na", "missing_row", "explicit_na", "missing_row"))
  expect_identical(cell_st$in_known_calendar, c(TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_false(is.na(cell_st$obs_index[1L]))
  expect_true(is.na(cell_st$obs_index[2L]))
  expect_false(is.na(cell_st$explicit_na_index[2L]))
  expect_true(is.na(cell_st$explicit_na_index[3L]))

  # Verify outcome-free and outcome-reading pair position traversal against Cartesian reference
  tr_env <- new.env(parent = emptyenv())
  batches_no_y <- list()
  res_no_y <- ns$traverse_ipw_panel_pair_positions(
    panel, control = ctrl,
    callback = function(b) batches_no_y[[length(batches_no_y) + 1L]] <<- b,
    read_outcomes = FALSE,
    trace_env = tr_env
  )
  expect_identical(res_no_y$endpoints_read, 0L)
  expect_identical(tr_env$calls, 0L)
  df_no_y <- do.call(rbind, batches_no_y)
  expect_false("y_s" %in% names(df_no_y))
  expect_identical(df_no_y$obs_s, ref$cartesian_pairs$obs_s)
  expect_identical(df_no_y$obs_t, ref$cartesian_pairs$obs_t)
  expect_identical(df_no_y$n_intervening_obs, ref$cartesian_pairs$n_intervening_obs)
  expect_identical(df_no_y$n_intervening_panel_dates, ref$cartesian_pairs$n_intervening_panel_dates)
  expect_identical(df_no_y$is_consecutive_unit, ref$cartesian_pairs$is_consecutive_unit)
  expect_identical(df_no_y$is_consecutive_panel, ref$cartesian_pairs$is_consecutive_panel)

  batches_with_y <- list()
  res_with_y <- ns$traverse_ipw_panel_pair_positions(
    panel, control = ctrl,
    callback = function(b) batches_with_y[[length(batches_with_y) + 1L]] <<- b,
    read_outcomes = TRUE,
    trace_env = tr_env
  )
  df_with_y <- do.call(rbind, batches_with_y)
  expect_identical(res_with_y$endpoints_read, 2L * nrow(ref$cartesian_pairs))
  expect_identical(df_with_y$y_s, ref$cartesian_pairs$y_s)
  expect_identical(df_with_y$y_t, ref$cartesian_pairs$y_t)

  # Construct frame owner with scrambled ID order + singleton/unobserved + external frame-only unit (id = 777)
  sup_id <- c(777, 109, 108, 107, 105, 104, 102, 101)
  Z_sup <- cbind(1, seq_along(sup_id) * 0.5)
  w_sup <- c(2.0, 1.5, 0.8, 1.2, 0.9, 1.1, 1.3, 0.7)
  fg_sup <- c(Inf, Inf, 25.75, Inf, 25.75, Inf, 40.0, 25.75)

  # Without allow_frame_only = TRUE, external unit 777 is refused
  expect_error(
    ns$construct_ipw_frame_owner(panel, id = sup_id, Z = Z_sup, weights = w_sup, control = ctrl, allow_frame_only = FALSE),
    class = "cdid_ipw_owner_id_error"
  )

  frame <- ns$construct_ipw_frame_owner(
    panel, id = sup_id, Z = Z_sup, weights = w_sup,
    control = ctrl, frame_g = fg_sup, allow_frame_only = TRUE, frame_id = "f_irreg"
  )
  ns$validate_ipw_frame_owner(frame, panel_owner = panel)
  expect_identical(frame$n_a_units, 5L)
  expect_identical(frame$n_singleton_units, 1L)
  expect_identical(frame$n_unobserved_units, 1L)
  expect_identical(frame$n_frame_only_units, 1L)
  expect_identical(frame$frame_only_units$id, 777)
  expect_identical(frame$frame_only_units$weight, 2.0)

  # Build nuisance view and link view on (s = 10.0, t = 35.0) with g = 25.75, K = c(Inf)
  # Note: among A-units:
  # - a_index 1 (id 101, g = 25.75): observed at 10.0 and 50.25 -> in risk R, NOT observed on (10.0, 35.0) (risk non-observer!)
  # - a_index 2 (id 102, g = 40.0): observed at 10.0 and 35.0 -> outside risk R (off-risk pair observer!)
  # - a_index 3 (id 105, g = 25.75): observed at 10.0 and 35.0 -> in O = R intersect S
  # - a_index 4 (id 107, g = Inf): observed at 10.0 and 35.0 -> in O = R intersect S
  # - a_index 5 (id 109, g = Inf): observed at 10.0 and 35.0 -> in O = R intersect S
  nu_view <- ns$construct_ipw_nuisance_view(
    frame, g = 25.75, K = c(Inf), s = 10.0, t = 35.0, control = ctrl, nuisance_id = "nu_10_35"
  )
  lk_view <- ns$construct_ipw_link_view(
    frame, s = 10.0, t = 35.0, nuisance_view = nu_view, control = ctrl, link_id = "lk_10_35"
  )
  expect_identical(nu_view$support$n_risk_a_units, 4L)
  expect_identical(nu_view$support$n_outside_risk_a_units, 1L)
  expect_identical(nu_view$support$n_observed_pair_all, 4L)
  expect_identical(nu_view$support$n_observed_risk_units, 3L)
  expect_identical(nu_view$support$n_nonobserved_risk_units, 1L)

  # Test contiguous chunk queries across arbitrary boundaries (B = 1, 2, 3, 5)
  for (b_size in c(1L, 2L, 3L, 5L)) {
    c_b <- default_owner_test_control(ns, max_chunk_rows = b_size)
    a_s <- 1L
    while (a_s <= frame$n_a_units) {
      a_e <- min(frame$n_a_units, a_s + b_size - 1L)
      ch <- ns$query_ipw_view_chunk(
        lk_view, a_start = a_s, a_end = a_e, control = c_b,
        read_outcomes = TRUE, include_design = TRUE, enclose_differences = TRUE
      )
      h_ch <- hand_reference_risk_first_chunk(
        ref, frame$Z_a, frame$w_a, g = 25.75, K = c(Inf),
        s = 10.0, t = 35.0, a_start = a_s, a_end = a_e, read_outcomes = TRUE
      )
      expect_identical(ch$a_index, h_ch$a_index)
      expect_identical(ch$unit_index, h_ch$unit_index)
      expect_identical(ch$id, h_ch$id)
      expect_identical(ch$g, h_ch$g)
      expect_identical(ch$weight, h_ch$weight)
      expect_identical(ch$cohort_index, h_ch$cohort_index)
      expect_identical(ch$in_risk, h_ch$in_risk)
      expect_identical(ch$is_target, h_ch$is_target)
      expect_identical(ch$in_observed, h_ch$in_observed)
      expect_identical(ch$is_risk_nonobserver, h_ch$is_risk_nonobserver)
      expect_identical(ch$chunk_to_obs_pos, h_ch$chunk_to_obs_pos)
      expect_identical(ch$observed_risk$local_pos, h_ch$local_pos)
      expect_identical(ch$observed_risk$a_index, h_ch$obs_a_index)
      expect_identical(ch$observed_risk$obs_s_index, h_ch$obs_s_index)
      expect_identical(ch$observed_risk$obs_t_index, h_ch$obs_t_index)
      expect_identical(ch$observed_risk$y_left_obs, h_ch$y_left_obs)
      expect_identical(ch$observed_risk$y_right_obs, h_ch$y_right_obs)
      expect_identical(ch$read_diagnostics$n_endpoints_read, 2L * length(h_ch$obs_a_index))
      a_s <- a_e + 1L
    }
  }
})

# ---------------------------------------------------------------------------
# Test 3: Character IDs (radix byte order) and large exact binary64 numeric IDs
# ---------------------------------------------------------------------------

test_that("character IDs and large exact binary64 numeric IDs preserve exact distinctions without coercion", {
  ns <- require_cdid_namespace()
  ctrl <- default_owner_test_control(ns, max_chunk_rows = 4L)

  # 1. Character IDs where lexicographic radix order differs from numeric order ("id_10" < "id_2")
  chr_ids <- c("id_2", "id_10", "id_02", "A_unit", "a_unit")
  df_chr <- do.call(rbind, lapply(seq_along(chr_ids), function(i) {
    data.frame(
      id = rep(chr_ids[i], 2L),
      time = c(1, 2),
      y = c(i * 1.0, i * 2.5),
      g = rep(if (i <= 2L) 2 else Inf, 2L),
      stringsAsFactors = FALSE
    )
  }))
  p_chr <- ns$construct_ipw_panel_owner(df_chr, "id", "time", "y", "g", control = ctrl)
  expect_identical(p_chr$id_type, "character")
  expect_identical(p_chr$a_id, sort(chr_ids, method = "radix"))

  # Numeric ID query on character panel owner is refused (no string coercion)
  expect_error(
    ns$query_ipw_panel_cell_status(p_chr, id = c(2, 10), time = c(1, 1), control = ctrl),
    class = "cdid_ipw_owner_input_error"
  )
  expect_error(
    ns$construct_ipw_frame_owner(p_chr, id = seq_along(chr_ids), Z = matrix(1, 5, 1), control = ctrl),
    class = "cdid_ipw_owner_id_error"
  )

  f_chr <- ns$construct_ipw_frame_owner(
    p_chr, id = rev(chr_ids), Z = cbind(1, seq_along(chr_ids)), control = ctrl
  )
  st_chr <- ns$query_ipw_panel_cell_status(
    p_chr, id = chr_ids, time = rep(1, 5),
    control = default_owner_test_control(ns, max_chunk_rows = 5L)
  )
  expect_true(all(st_chr$status == "observed"))
  expect_identical(p_chr$units$id[st_chr$unit_index], chr_ids)

  # 2. Large exact binary64 integer IDs near 2^53 - 1 (where 32-bit int would overflow)
  big_ids <- c(2^53 - 1, 2^53 - 2, 2^52 + 7, 2^52 + 1)
  df_big <- do.call(rbind, lapply(seq_along(big_ids), function(i) {
    data.frame(
      id = rep(big_ids[i], 2L),
      time = c(2000, 2001),
      y = c(10 + i, 20 + i),
      g = rep(if (i <= 2L) 2001 else Inf, 2L)
    )
  }))
  p_big <- ns$construct_ipw_panel_owner(df_big, "id", "time", "y", "g", control = ctrl)
  expect_identical(p_big$id_type, "numeric")
  expect_identical(p_big$a_id, sort(big_ids))
  expect_identical(diff(p_big$a_id[3:4]), 1.0)

  # Character ID query on numeric panel owner is refused (no coercion)
  expect_error(
    ns$construct_ipw_frame_owner(p_big, id = as.character(big_ids), Z = matrix(1, 4, 1), control = ctrl),
    class = "cdid_ipw_owner_id_error"
  )
  f_big <- ns$construct_ipw_frame_owner(p_big, id = big_ids, Z = cbind(1, 1:4), control = ctrl)
  lk_big <- ns$construct_ipw_link_view(f_big, s = 2000, t = 2001, g = 2001, K = c(Inf), control = ctrl)
  ch_big <- ns$query_ipw_view_chunk(lk_big, 1L, 4L, control = ctrl, read_outcomes = TRUE, expected_a_id = sort(big_ids))
  expect_identical(ch_big$id, sort(big_ids))
})

# ---------------------------------------------------------------------------
# Test 4: Overflow-prone off-risk outcomes & endpoint-read instrumentation
# ---------------------------------------------------------------------------

test_that("risk-first filtering avoids reading off-risk endpoints and prevents off-risk subtraction overflow", {
  ns <- require_cdid_namespace()
  ctrl <- default_owner_test_control(ns, max_chunk_rows = 3L)

  # Units 1..2: g = 2 (target)
  # Units 3..4: g = Inf (control in K)
  # Units 5..6: g = 3 (OFF-RISK cohort not in {g = 2, K = Inf}!), with finite extreme
  # endpoints y(1) = -.Machine$double.xmax and y(2) = +.Machine$double.xmax whose
  # difference y(2) - y(1) overflows IEEE binary64 to +Inf!
  big_val <- .Machine$double.xmax
  df_ovf <- rbind(
    data.frame(id = c(1, 1, 2, 2), time = c(1, 2, 1, 2), y = c(1.0, 2.5, 1.5, 3.0), g = rep(2, 4)),
    data.frame(id = c(3, 3, 4, 4), time = c(1, 2, 1, 2), y = c(0.5, 1.0, 2.0, 2.5), g = rep(Inf, 4)),
    data.frame(id = c(5, 5, 6, 6), time = c(1, 2, 1, 2), y = c(-big_val, big_val, -big_val, big_val), g = rep(3, 4))
  )

  # 1. Existing construct_observed_pairs subtracts all pairs globally -> aborts on off-risk units 5..6!
  expect_error(
    ns$construct_observed_pairs(df_ovf, "id", "time", "y", "g"),
    class = "cdid_pairs_outcome_error"
  )

  # 2. New construct_ipw_panel_owner validates raw finite endpoints WITHOUT global subtraction
  panel <- ns$construct_ipw_panel_owner(df_ovf, "id", "time", "y", "g", control = ctrl)
  expect_identical(panel$accounting$global_pair_subtraction_performed, FALSE)
  expect_identical(panel$accounting$endpoints_read_during_ingestion, 0L)

  frame <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), control = ctrl)
  nu_view <- ns$construct_ipw_nuisance_view(frame, g = 2, K = c(Inf), s = 1, t = 2, control = ctrl)
  lk_view <- ns$construct_ipw_link_view(frame, s = 1, t = 2, nuisance_view = nu_view, control = ctrl)

  # 3. Outcome-free query reads ZERO endpoints
  tr_free <- new.env(parent = emptyenv())
  ch_free <- ns$query_ipw_view_chunk(
    lk_view, a_start = 1L, a_end = 6L, control = default_owner_test_control(ns, 6L),
    read_outcomes = FALSE, trace_env = tr_free
  )
  expect_identical(tr_free$calls, 0L)
  expect_identical(ch_free$read_diagnostics$n_endpoints_read, 0L)
  expect_null(ch_free$observed_risk$y_left_obs)
  expect_null(ch_free$observed_risk$y_right_obs)
  expect_identical(ch_free$counts$outside_risk_pair_hit_count, 2L)

  # 4. Outcome-reading query on {g = 2, K = Inf} filters risk BEFORE reading y,
  #    reads ONLY observation indices 1..8 (units 1..4), and NEVER touches indices 9..12 (units 5..6)!
  tr_read <- new.env(parent = emptyenv())
  batches <- list()
  st_res <- ns$stream_ipw_view_chunks(
    lk_view, control = ctrl,
    callback = function(ch) batches[[length(batches) + 1L]] <<- ch,
    read_outcomes = TRUE, include_design = TRUE, enclose_differences = TRUE,
    trace_env = tr_read
  )
  expect_identical(st_res$total_endpoints_read, 8L)
  read_s <- unlist(lapply(batches, function(ch) ch$read_diagnostics$read_obs_s_index))
  read_t <- unlist(lapply(batches, function(ch) ch$read_diagnostics$read_obs_t_index))
  expect_identical(read_s, c(1L, 3L, 5L, 7L))
  expect_identical(read_t, c(2L, 4L, 6L, 8L))
  expect_false(any(read_s %in% 9:12))
  expect_false(any(read_t %in% 9:12))

  # 5. If cohort 3 IS included in the risk set and enclose_differences = TRUE is requested,
  #    outward interval subtraction detects binary64 overflow on unit 5 and aborts cleanly!
  lk_with_3 <- ns$construct_ipw_link_view(frame, s = 1, t = 2, g = 2, K = c(3, Inf), control = ctrl)
  expect_error(
    ns$query_ipw_view_chunk(
      lk_with_3, a_start = 4L, a_end = 6L, control = ctrl,
      read_outcomes = TRUE, enclose_differences = TRUE
    ),
    class = "cdid_ipw_owner_outcome_overflow_error"
  )
})

# ---------------------------------------------------------------------------
# Test 5: Allowed-domain equivalence with accepted helpers (without off-risk reads)
# ---------------------------------------------------------------------------

test_that("allowed-domain metadata and positions match construct_observed_pairs, construct_ipw_context, and query_observed_pair", {
  ns <- require_cdid_namespace()
  ctrl <- default_owner_test_control(ns, max_chunk_rows = 5L)

  fixture_path <- testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R")
  ref_fix <- dget(fixture_path)$fixture1_exact_rational

  # Build long panel data frame from fixture1_exact_rational (10 A-units observed at baseline dates 0 and 10,
  # and observed on (1, 2) iff pair_observed[i] == TRUE) plus a singleton unit 901
  rows_fix <- lapply(seq_along(ref_fix$id), function(i) {
    uid <- ref_fix$id[i]
    gi <- as.double(ref_fix$G[i])
    y1 <- if (ref_fix$pair_observed[i]) ref_fix$yL_1_float[i] else NA_real_
    y2 <- if (ref_fix$pair_observed[i]) ref_fix$yR_1_float[i] else NA_real_
    data.frame(
      id = rep(uid, 4L),
      time = c(0, 1, 2, 10),
      y = c(0.0, y1, y2, 1.0),
      g = rep(gi, 4L)
    )
  })
  df_fix <- rbind(
    do.call(rbind, rows_fix),
    data.frame(id = 901L, time = 1, y = 42.0, g = 2)
  )
  sup_ids <- c(ref_fix$id, 901L)
  Z_sup <- rbind(ref_fix$Z_float, c(1, -0.5))
  w_sup <- c(ref_fix$w_float, 0.5)

  old_pairs <- ns$construct_observed_pairs(df_fix, "id", "time", "y", "g")
  old_ctx <- ns$construct_ipw_context(old_pairs, id = sup_ids, Z = Z_sup, weights = w_sup)
  old_qp <- ns$query_observed_pair(old_ctx, s = 1, t = 2)

  panel <- ns$construct_ipw_panel_owner(df_fix, "id", "time", "y", "g", control = ctrl)
  frame <- ns$construct_ipw_frame_owner(
    panel, id = sup_ids, Z = Z_sup, weights = w_sup, control = ctrl
  )

  expect_identical(panel$a_id, old_ctx$a_id)
  expect_identical(panel$a_g, old_ctx$a_g)
  expect_identical(panel$all_times, old_ctx$all_times)
  expect_identical(panel$times, old_ctx$times)
  expect_identical(panel$cohorts, old_ctx$cohorts)
  expect_identical(panel$singletons[, names(old_pairs$singletons)], old_pairs$singletons)
  expect_identical(panel$observations[, names(old_pairs$observations)], old_pairs$observations)
  expect_identical(frame$Z_a, old_ctx$Z_a)
  expect_identical(frame$w_a, old_ctx$w_a)
  expect_identical(frame$cohort_a_count, old_ctx$cohort_a_count)
  expect_equal(frame$cohort_a_weight, old_ctx$cohort_a_weight, tolerance = 0)
  expect_equal(frame$cohort_a_share, old_ctx$cohort_a_share, tolerance = 0)

  # Outcome-free pair position query matches query_observed_pair positions
  pos_df <- ns$query_ipw_panel_pair_positions(panel, s = 1, t = 2, control = default_owner_test_control(ns, 20L))
  expect_identical(pos_df$a_index, old_qp$a_index)
  expect_identical(pos_df$obs_s, old_qp$obs_s)
  expect_identical(pos_df$obs_t, old_qp$obs_t)

  # Verify cohort shares q and pair observation rates r match fixture1_exact_rational
  nu_f1 <- ns$construct_ipw_nuisance_view(frame, g = 2, K = c(3, 4), s = 1, t = 2, control = default_owner_test_control(ns, 10L))
  expect_equal(unname(nu_f1$support$q), ref_fix$q_cohort_float, tolerance = 1e-15)
  expect_equal(unname(nu_f1$support$rates), ref_fix$r_cohort_float, tolerance = 1e-15)

  # Unlike query_observed_pair (which reads y_s and y_t for off-risk unit 9 with G = 10),
  # query_ipw_view_chunk on {g = 2, K = c(3, 4)} reads ONLY risk-set observers!
  expect_true(9L %in% old_qp$a_index)
  lk <- ns$construct_ipw_link_view(frame, s = 1, t = 2, nuisance_view = nu_f1, control = default_owner_test_control(ns, 10L))
  tr <- new.env(parent = emptyenv())
  ch <- ns$query_ipw_view_chunk(lk, 1L, frame$n_a_units, control = default_owner_test_control(ns, 10L), read_outcomes = TRUE, trace_env = tr)
  u9_obs_rows <- which(panel$observations$unit_index == panel$units$unit_index[panel$units$id == 9L])
  expect_false(any(tr$pos_s_index %in% u9_obs_rows))
  expect_false(any(tr$pos_t_index %in% u9_obs_rows))
})
