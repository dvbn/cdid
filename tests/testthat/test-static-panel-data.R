# Data-preparation tests for R/static_panel_data.R.
#
# Resolves functions from the active test / package namespace environment and
# reads independent hand-computed reference data exclusively via
# `testthat::test_path()`. Never reads plaintext production R files.

resolve_static_panel_fn <- function(name) {
  if (exists(name, envir = parent.frame(), mode = "function", inherits = TRUE)) {
    return(get(name, envir = parent.frame(), mode = "function", inherits = TRUE))
  }
  if ("cdid" %in% loadedNamespaces() &&
      exists(name, envir = asNamespace("cdid"), mode = "function", inherits = FALSE)) {
    return(get(name, envir = asNamespace("cdid"), mode = "function", inherits = FALSE))
  }
  stop(sprintf("Required function `%s` was not found in the test or `cdid` namespace environment.", name))
}

static_panel_data <- resolve_static_panel_fn("static_panel_data")
static_panel_normalize <- resolve_static_panel_fn("static_panel_normalize")
static_transition_cell <- resolve_static_panel_fn("static_transition_cell")
static_transition_cells <- resolve_static_panel_fn("static_transition_cells")
static_extract_reference_design <- resolve_static_panel_fn("static_extract_reference_design")
static_reference_design <- resolve_static_panel_fn("static_reference_design")
static_reference_design_extractor <- resolve_static_panel_fn("static_reference_design_extractor")
static_observation_pair <- resolve_static_panel_fn("static_observation_pair")
static_observation_pairs <- resolve_static_panel_fn("static_observation_pairs")
static_observation_pair_accessor <- resolve_static_panel_fn("static_observation_pair_accessor")

build_hand_5date_long_df <- function(ref) {
  grid <- expand.grid(
    u_idx = seq_along(ref$unit_id),
    t_idx = seq_along(ref$dates),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n_rows <- nrow(grid)
  id_col <- ref$unit_id[grid$u_idx]
  t_col <- ref$dates[grid$t_idx]
  y_col <- numeric(n_rows)
  d_col <- integer(n_rows)
  x1_col <- numeric(n_rows)
  x2_col <- numeric(n_rows)
  w_sel_col <- numeric(n_rows)
  sector_col <- character(n_rows)

  sectors_by_unit <- c("mfg", "svc", "tech", "mfg", "svc", "tech", "mfg", "svc", "tech", "mfg")

  for (r in seq_len(n_rows)) {
    i <- grid$u_idx[r]
    s <- grid$t_idx[r]
    y_col[r] <- ref$Y_matrix[i, s]
    d_col[r] <- ref$D_matrix[i, s]
    x1_col[r] <- 10 * i + s
    x2_col[r] <- (-1)^i * (2 * s + 0.5 * i)
    w_sel_col[r] <- 100 + 3 * i - s
    sector_col[r] <- sectors_by_unit[i]
  }

  data.frame(
    unit = id_col,
    year = t_col,
    outcome = y_col,
    treated = d_col,
    x1 = x1_col,
    x2 = x2_col,
    w_sel = w_sel_col,
    sector = factor(sector_col, levels = c("mfg", "svc", "tech")),
    sector_chr = sector_col,
    stringsAsFactors = FALSE
  )
}

test_that("5-date hand panel normalizes Y/D, preserves zero/single-outcome units, and matches independent reference", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)

  panel <- static_panel_data(
    data = df_long,
    idname = "unit",
    tname = "year",
    yname = "outcome",
    dname = "treated"
  )

  expect_s3_class(panel, "cdid_static_panel_data")
  expect_identical(panel$n, 10L)
  expect_identical(panel$T_dates, 5L)
  expect_identical(panel$unit_id, ref$unit_id)
  expect_identical(panel$unit_id_storage_mode, "character")
  expect_false(panel$unit_id_is_num)
  expect_identical(panel$unit_id_source, "sorted_unique_from_data")
  expect_identical(panel$dates, ref$dates)
  expect_identical(panel$date_grid_source, "extracted_unique_sorted")
  expect_false(panel$explicit_dates)

  expect_identical(panel$D, ref$D_matrix)
  expect_equal(panel$Y, ref$Y_matrix)
  expect_identical(panel$unit_outcome_counts, ref$unit_outcome_counts)
  expect_identical(panel$summary$n_units_zero_outcomes, 1L)
  expect_identical(panel$summary$n_units_single_outcome, 1L)
  expect_identical(panel$summary$n_d_missing_on_present_rows, 1L)
  expect_identical(panel$summary$n_cells_absent, 0L)

  # Unit u10 at 2025 has present row with D = NA: explicitly recorded as missing_d, never untreated_0.
  expect_identical(panel$d_status["u10", "2025"], "missing_d")
  expect_false(panel$d_known["u10", "2025"])
  expect_true(is.na(panel$D["u10", "2025"]))
})

test_that("transition cells match hand reference across entry, exit, re-entry, absorbing, constant, and shared first-entry histories", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)
  panel <- static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated")

  req_table <- data.frame(
    t = c(2012, 2012, 2015, 2015, 2019, 2019),
    d = c(0L, 1L, 0L, 1L, 0L, 1L),
    stringsAsFactors = FALSE
  )
  cells_obj <- static_transition_cells(panel, cells = req_table)
  expect_s3_class(cells_obj, "cdid_static_transition_cells")
  expect_identical(cells_obj$n_cells, 6L)
  expect_identical(cells_obj$cell_ids, names(ref$cells))

  for (cid in names(ref$cells)) {
    got <- cells_obj$cells[[cid]]
    exp_c <- ref$cells[[cid]]
    expect_identical(got$t, exp_c$t)
    expect_identical(got$t_prev, exp_c$t_prev)
    expect_identical(got$t_idx, exp_c$t_idx)
    expect_identical(got$t_prev_idx, exp_c$t_prev_idx)
    expect_identical(got$d, exp_c$d)
    expect_identical(got$sigma, exp_c$sigma)
    expect_identical(got$R, exp_c$R)
    expect_identical(got$G, exp_c$G)
    expect_identical(got$C, exp_c$C)
    expect_identical(got$n_risk, exp_c$n_risk)
    expect_identical(got$n_switchers, exp_c$n_switchers)
    expect_identical(got$n_stayers, exp_c$n_stayers)
    expect_identical(got$support$status, exp_c$status)
    expect_true(got$support$is_supported)
  }

  # Units u03, u04, u05 share identical first-entry date (2012) and are all switchers at c = (2012, 0),
  # but diverge completely at later transitions:
  # - At exit cell c = (2015, 1): u03 is a stayer (C=1), u04 and u05 are exit switchers (G=1).
  # - At entry cell c = (2019, 0): u03 is outside risk set (R=0), u04 is a re-entry switcher (G=1), u05 is a stayer (C=1).
  c_2012_0 <- cells_obj$cells[["t=2012,d=0"]]
  c_2015_1 <- cells_obj$cells[["t=2015,d=1"]]
  c_2019_0 <- cells_obj$cells[["t=2019,d=0"]]

  expect_identical(c_2012_0$G[3:5], c(1L, 1L, 1L))
  expect_identical(c_2015_1$C[3], 1L)
  expect_identical(c_2015_1$G[4:5], c(1L, 1L))
  expect_identical(c(c_2019_0$R[3], c_2019_0$G[4], c_2019_0$C[5]), c(0L, 1L, 1L))

  # Zero-outcome unit u02 (index 2) and single-outcome unit u05 (index 5) retain their full population roles.
  expect_identical(c_2012_0$C[2], 1L)
  expect_identical(c_2012_0$G[5], 1L)
  expect_identical(c_2015_1$G[5], 1L)
  expect_identical(c_2019_0$C[5], 1L)
})

test_that("meaningful mutations of transition-role definition and gap indexing are discriminated by reference checks", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)
  panel <- static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  cells_obj <- static_transition_cells(
    panel,
    t = c(2012, 2012, 2015, 2015, 2019, 2019),
    d = c(0L, 1L, 0L, 1L, 0L, 1L)
  )

  # Mutation A (First-switch / absorbing cohort role mutation):
  # Redefine entry switcher G_c at t to require never having been treated before t-1,
  # and stayer C_c to require D_1 = ... = D_t = 0.
  mutated_first_switch_G_2015 <- as.integer(
    panel$D[, 1L] == 0L & panel$D[, 2L] == 0L & panel$D[, 3L] == 1L
  )
  mutated_first_switch_G_2019 <- as.integer(
    panel$D[, 1L] == 0L & panel$D[, 2L] == 0L & panel$D[, 3L] == 0L & panel$D[, 4L] == 1L
  )
  mutated_never_treated_C_2019 <- as.integer(
    panel$D[, 1L] == 0L & panel$D[, 2L] == 0L & panel$D[, 3L] == 0L & panel$D[, 4L] == 0L
  )

  # Reference static masks include re-entrant u06 at (2015, 0), re-entrant u04 at (2019, 0),
  # and previously treated stayers u05 and u08 at (2019, 0). The mutated masks disagree!
  expect_identical(cells_obj$cells[["t=2015,d=0"]]$G, ref$cells[["t=2015,d=0"]]$G)
  expect_false(identical(mutated_first_switch_G_2015, ref$cells[["t=2015,d=0"]]$G))
  expect_false(identical(mutated_first_switch_G_2019, ref$cells[["t=2019,d=0"]]$G))
  expect_false(identical(mutated_never_treated_C_2019, ref$cells[["t=2019,d=0"]]$C))

  # Mutation B (Outcome-observed filtering of R/G/C masks):
  # Filtering G_c by adjacent outcome availability O_{t-1} * O_t drops u04 and u05 at (2012, 0).
  O_adj_2012 <- as.integer(panel$O[, 1L] == 1L & panel$O[, 2L] == 1L)
  mutated_obs_filtered_G_2012 <- cells_obj$cells[["t=2012,d=0"]]$G * O_adj_2012
  expect_false(identical(mutated_obs_filtered_G_2012, ref$cells[["t=2012,d=0"]]$G))

  # Mutation C (Calendar numeric difference u - s instead of panel date-index gap u_idx - s_idx):
  c_2015_1 <- cells_obj$cells[["t=2015,d=1"]]
  pairs_got <- static_observation_pairs(
    panel = panel,
    mask = c_2015_1$G,
    s = c(2012, 2010, 2012, 2010),
    u = c(2015, 2015, 2019, 2025)
  )
  ref_gaps <- vapply(ref$pairs_c2015_d1_G, `[[`, integer(1L), "gap")
  got_gaps <- stats::setNames(pairs_got$summary_table$gap, pairs_got$pair_ids)
  mutated_calendar_gaps <- stats::setNames(
    as.integer(pairs_got$summary_table$u - pairs_got$summary_table$s),
    pairs_got$pair_ids
  )
  expect_identical(got_gaps, ref_gaps)
  expect_identical(unname(got_gaps), c(1L, 2L, 2L, 4L))
  expect_identical(unname(mutated_calendar_gaps), c(3L, 5L, 7L, 15L))
  expect_false(identical(mutated_calendar_gaps, ref_gaps))
})

test_that("shuffled long rows, numeric/character IDs, and explicit vs extracted date grids produce invariant aligned structures", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)

  # Deterministic permutation of long rows (reordering units and dates arbitrarily).
  perm <- c(
    42L, 7L, 19L, 1L, 50L, 33L, 12L, 28L, 5L, 47L,
    23L, 14L, 39L, 8L, 31L, 45L, 2L, 17L, 36L, 21L,
    10L, 49L, 26L, 4L, 38L, 15L, 41L, 29L, 6L, 34L,
    48L, 11L, 24L, 3L, 37L, 16L, 44L, 20L, 9L, 32L,
    27L, 13L, 46L, 22L, 35L, 18L, 40L, 25L, 30L, 43L
  )
  expect_identical(sort(perm), seq_len(nrow(df_long)))
  df_shuffled <- df_long[perm, , drop = FALSE]
  rownames(df_shuffled) <- NULL

  p_orig <- static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  p_shuf <- static_panel_data(df_shuffled, idname = "unit", tname = "year", yname = "outcome", dname = "treated")

  expect_identical(p_shuf$unit_id, p_orig$unit_id)
  expect_identical(p_shuf$dates, p_orig$dates)
  expect_identical(p_shuf$D, p_orig$D)
  expect_identical(p_shuf$Y, p_orig$Y)
  expect_identical(p_shuf$O, p_orig$O)
  expect_identical(p_shuf$d_status, p_orig$d_status)

  c_orig <- static_transition_cell(p_orig, t = 2015, d = 1)
  c_shuf <- static_transition_cell(p_shuf, t = 2015, d = 1)
  expect_identical(c_shuf$R, c_orig$R)
  expect_identical(c_shuf$G, c_orig$G)
  expect_identical(c_shuf$C, c_orig$C)

  des_orig <- static_reference_design(p_orig, c_orig, formula = ~ x1 + sector * x2)
  des_shuf <- static_reference_design(p_shuf, c_shuf, formula = ~ x1 + sector * x2)
  expect_identical(des_shuf$X_tilde, des_orig$X_tilde)
  expect_identical(des_shuf$Z, des_orig$Z)

  pairs_orig <- static_observation_pairs(p_orig, c_orig$G, s = c(2010, 2012), u = c(2015, 2019))
  pairs_shuf <- static_observation_pairs(p_shuf, c_shuf$G, s = c(2010, 2012), u = c(2015, 2019))
  expect_identical(pairs_shuf$summary_table, pairs_orig$summary_table)
  expect_identical(pairs_shuf$pairs[["(2010,2015)"]]$dy, pairs_orig$pairs[["(2010,2015)"]]$dy)
  expect_identical(pairs_shuf$pairs[["(2012,2019)"]]$O, pairs_orig$pairs[["(2012,2019)"]]$O)

  # Test numeric integer and double unit IDs and explicit caller unit_ids ordering.
  df_int_id <- df_shuffled
  df_int_id$unit_int <- as.integer(sub("^u", "", df_int_id$unit)) * 10L
  p_int <- static_panel_data(
    df_int_id,
    idname = "unit_int",
    tname = "year",
    yname = "outcome",
    dname = "treated",
    dates = c(2010L, 2012L, 2015L, 2019L, 2025L)
  )
  expect_identical(typeof(p_int$unit_id), "integer")
  expect_true(p_int$unit_id_is_num)
  expect_identical(p_int$unit_id, seq.int(10L, 100L, by = 10L))
  expect_identical(p_int$date_grid_source, "explicit_caller")
  expect_true(p_int$explicit_dates)

  # Custom explicit caller unit_ids order (reversed) preserves caller's exact unit ordering.
  rev_ids <- rev(ref$unit_id)
  p_rev <- static_panel_data(
    df_shuffled,
    idname = "unit",
    tname = "year",
    yname = "outcome",
    dname = "treated",
    unit_ids = rev_ids
  )
  expect_identical(p_rev$unit_id, rev_ids)
  expect_identical(p_rev$unit_id_source, "explicit_caller")
  expect_identical(p_rev$D, p_orig$D[rev_ids, , drop = FALSE])
  expect_identical(p_rev$Y, p_orig$Y[rev_ids, , drop = FALSE])
})

test_that("4-date hand panel supports empty switcher and empty stayer cells as diagnostics without dropping cells", {
  # Construct a 4-date panel with dates c(1, 4, 9, 16) where:
  # - At t = 4, d = 0: both switchers and stayers exist (supported).
  # - At t = 9, d = 1: all units with D_{t-1}=1 stay treated (n_switchers = 0, status = "empty_switchers").
  # - At t = 16, d = 0: all units with D_{t-1}=0 switch to 1 (n_stayers = 0, status = "empty_stayers").
  # - At t = 4, d = 1: no unit has D_1 = 1 (n_risk = 0, status = "empty_risk").
  df4 <- data.frame(
    id = rep(1:4, each = 4L),
    t = rep(c(1, 4, 9, 16), times = 4L),
    y = c(
      1, 2, 3, 4,
      2, NA, 5, 6,
      NA, NA, NA, NA,
      3, 4, 6, NA
    ),
    d = c(
      0L, 1L, 1L, 1L, # unit 1: enters at t=4, stays 1
      0L, 1L, 1L, 0L, # unit 2: enters at t=4, stays 1 at t=9, exits at t=16
      0L, 0L, 0L, 1L, # unit 3: stays 0 through t=9, enters at t=16
      0L, 0L, 0L, 1L  # unit 4: stays 0 through t=9, enters at t=16
    )
  )

  p4 <- static_panel_normalize(df4, idname = "id", tname = "t", yname = "y", dname = "d")
  req4 <- data.frame(
    t = c(4, 4, 9, 16),
    d = c(0L, 1L, 1L, 0L)
  )
  cells4 <- static_transition_cells(p4, cells = req4)

  # All 4 requested cells are retained in exact requested order; none are dropped or altered.
  expect_identical(cells4$n_cells, 4L)
  expect_identical(cells4$cell_ids, c("t=4,d=0", "t=4,d=1", "t=9,d=1", "t=16,d=0"))
  expect_identical(cells4$supported_cell_ids, "t=4,d=0")
  expect_identical(cells4$unsupported_cell_ids, c("t=4,d=1", "t=9,d=1", "t=16,d=0"))
  expect_identical(cells4$empty_risk_cell_ids, "t=4,d=1")
  expect_identical(cells4$empty_switcher_cell_ids, c("t=4,d=1", "t=9,d=1"))
  expect_identical(cells4$empty_stayer_cell_ids, c("t=4,d=1", "t=16,d=0"))

  expect_identical(cells4$cells[["t=4,d=0"]]$support$status, "supported")
  expect_identical(cells4$cells[["t=4,d=1"]]$support$status, "empty_risk")
  expect_identical(cells4$cells[["t=9,d=1"]]$support$status, "empty_switchers")
  expect_identical(cells4$cells[["t=9,d=1"]]$n_risk, 2L)
  expect_identical(cells4$cells[["t=9,d=1"]]$n_switchers, 0L)
  expect_identical(cells4$cells[["t=9,d=1"]]$n_stayers, 2L)

  expect_identical(cells4$cells[["t=16,d=0"]]$support$status, "empty_stayers")
  expect_identical(cells4$cells[["t=16,d=0"]]$n_risk, 2L)
  expect_identical(cells4$cells[["t=16,d=0"]]$n_switchers, 2L)
  expect_identical(cells4$cells[["t=16,d=0"]]$n_stayers, 0L)
})

test_that("unknown or absent treatment rows fail only when required by the requested transition", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)

  # Remove unit u09's row at year = 2025 completely (absent long-data row),
  # while unit u10 has a present row at year = 2025 with treated = NA.
  df_partial <- df_long[!(df_long$unit == "u09" & df_long$year == 2025), , drop = FALSE]
  panel <- static_panel_data(df_partial, idname = "unit", tname = "year", yname = "outcome", dname = "treated")

  expect_identical(panel$d_status["u09", "2025"], "absent_row")
  expect_identical(panel$d_status["u10", "2025"], "missing_d")
  expect_false(panel$row_present["u09", "2025"])
  expect_true(panel$row_present["u10", "2025"])
  expect_true(is.na(panel$D["u09", "2025"]))
  expect_true(is.na(panel$D["u10", "2025"]))

  # Transitions at 2012, 2015, and 2019 do NOT use date 2025; they succeed without error!
  c_2019_0 <- static_transition_cell(panel, t = 2019, d = 0)
  expect_identical(c_2019_0$R, ref$cells[["t=2019,d=0"]]$R)

  # Transition at 2025 (2019 -> 2025) requires D_2019 and D_2025 for all 10 study units.
  # It must fail with `cdid_static_missing_treatment_error` identifying both u09 (absent row) and u10 (missing D).
  err <- expect_error(
    static_transition_cell(panel, t = 2025, d = 0),
    class = "cdid_static_missing_treatment_error"
  )
  expect_identical(err$absent_curr_units, "u09")
  expect_identical(err$na_curr_units, "u10")
  expect_identical(err$missing_units, c("u09", "u10"))

  # Similarly, if a row at t_prev = 2015 is absent or has D = NA, both (2015, d) and (2019, d) fail,
  # while (2012, d) still succeeds!
  df_gap_2015 <- df_long
  df_gap_2015$treated[df_gap_2015$unit == "u01" & df_gap_2015$year == 2015] <- NA_integer_
  panel_gap <- static_panel_data(df_gap_2015, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  expect_s3_class(static_transition_cell(panel_gap, t = 2012, d = 0), "cdid_static_transition_cell")
  expect_error(static_transition_cell(panel_gap, t = 2015, d = 0), class = "cdid_static_missing_treatment_error")
  err_prev <- expect_error(static_transition_cell(panel_gap, t = 2019, d = 0), class = "cdid_static_missing_treatment_error")
  expect_identical(err_prev$na_prev_units, "u01")
})

test_that("duplicate (id, date) rows, invalid Y/D, and invalid transition dates/states are rejected", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)

  # Duplicate (unit, year) record.
  df_dup <- rbind(df_long, df_long[7L, , drop = FALSE])
  expect_error(
    static_panel_data(df_dup, idname = "unit", tname = "year", yname = "outcome", dname = "treated"),
    class = "cdid_static_panel_duplicate_error"
  )

  # Non-binary D (e.g. cohort first-treatment year 2012 or value 2 or -1).
  df_bad_d <- df_long
  df_bad_d$treated[3L] <- 2012L
  expect_error(
    static_panel_data(df_bad_d, idname = "unit", tname = "year", yname = "outcome", dname = "treated"),
    class = "cdid_static_panel_treatment_error"
  )

  # Logical D is rejected.
  df_log_d <- df_long
  df_log_d$treated <- df_log_d$treated == 1L
  expect_error(
    static_panel_data(df_log_d, idname = "unit", tname = "year", yname = "outcome", dname = "treated"),
    class = "cdid_static_panel_treatment_error"
  )

  # Invalid Y (Inf or NaN).
  df_inf_y <- df_long
  df_inf_y$outcome[1L] <- Inf
  expect_error(
    static_panel_data(df_inf_y, idname = "unit", tname = "year", yname = "outcome", dname = "treated"),
    class = "cdid_static_panel_outcome_error"
  )
  df_nan_y <- df_long
  df_nan_y$outcome[1L] <- NaN
  expect_error(
    static_panel_data(df_nan_y, idname = "unit", tname = "year", yname = "outcome", dname = "treated"),
    class = "cdid_static_panel_outcome_error"
  )

  # Invalid study dates (unsorted explicit dates or date outside explicit grid).
  expect_error(
    static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated",
                      dates = c(2010, 2015, 2012, 2019, 2025)),
    class = "cdid_static_panel_date_error"
  )
  expect_error(
    static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated",
                      dates = c(2010, 2012, 2015, 2019)),
    class = "cdid_static_panel_date_error"
  )

  # Transition cell date errors: initial date t = 2010 has no prior date; t = 2014 not in grid; wrong t_prev.
  panel <- static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  expect_error(static_transition_cell(panel, t = 2010, d = 0), class = "cdid_static_cell_date_error")
  expect_error(static_transition_cell(panel, t = 2014, d = 0), class = "cdid_static_cell_date_error")
  expect_error(static_transition_cell(panel, t = 2015, d = 0, t_prev = 2010), class = "cdid_static_cell_date_error")
  expect_error(static_transition_cell(panel, t = 2015, d = 2), class = "cdid_static_cell_state_error")
  expect_error(static_transition_cell(panel, t = 2015, d = 0, sigma = -1), class = "cdid_static_cell_state_error")
})

test_that("reference design extraction handles needed vs irrelevant X missingness, separate formulas, factor levels, interactions, and no-intercept refusal", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)

  # Cell c = (2015, 0) has reference date t_prev = 2012 (t_prev_idx = 2).
  # Its risk set R_c = 1 consists of units u01, u02, u06, u09, u10 (indices 1, 2, 6, 9, 10).
  # Units outside R_c (R_c = 0) are u03, u04, u05, u07, u08 (indices 3, 4, 5, 7, 8).
  #
  # Introduce:
  # 1. Irrelevant X missingness / invalid domain values at t_prev = 2012 on non-risk units u03, u04, u05:
  #    - x1 = NA on u03 at 2012, x1 = -999 on u04 at 2012 (would fail `log(x1)` if evaluated!), w_sel = NA on u05 at 2012.
  # 2. Historical X missingness at OTHER dates (2010, 2015, 2019, 2025) for ALL units, including risk units u01, u02, u06, u09, u10!
  # 3. Restrict `sector` on the risk units at 2012 so only "mfg" and "svc" appear in R_c,
  #    while "tech" is a known factor level in the panel.
  df_x <- df_long
  df_x$x1[df_x$year != 2012] <- NA_real_
  df_x$x2[df_x$year != 2012] <- NA_real_
  df_x$w_sel[df_x$year != 2012] <- NA_real_

  df_x$x1[df_x$unit == "u03" & df_x$year == 2012] <- NA_real_
  df_x$x1[df_x$unit == "u04" & df_x$year == 2012] <- -999.0
  df_x$w_sel[df_x$unit == "u05" & df_x$year == 2012] <- NA_real_
  df_x$sector[df_x$unit %in% c("u06", "u09") & df_x$year == 2012] <- factor("svc", levels = c("mfg", "svc", "tech"))

  panel <- static_panel_data(df_x, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  cell <- static_transition_cell(panel, t = 2015, d = 0)

  # Build lazy callable extractor with separate outcome, treatment, and observation/selection formulas.
  extractor <- static_reference_design(
    panel = panel,
    cell = NULL,
    outcome_formula = ~ log(x1) + sector * x2,
    treatment_formula = ~ log(x1) + x2,
    observation_formula = ~ log(x1) + x2 + w_sel + sector
  )
  expect_s3_class(extractor, "cdid_static_design_extractor")
  expect_true(is.function(extractor))

  designs <- extractor(cell, role = "all", population = "risk")
  expect_s3_class(designs, "cdid_static_cell_designs")
  expect_true(designs$separate_treatment_formula)
  expect_true(designs$separate_observation_formula)

  out_des <- designs$outcome
  trt_des <- designs$treatment
  obs_des <- designs$observation

  # Check dimensions, preserved known factor level "sectortech" (even though unobserved in R_c at 2012),
  # and collinear/zero column retention (no rank preprocessing).
  expect_identical(
    out_des$colnames_X_tilde,
    c("(Intercept)", "log(x1)", "sectorsvc", "sectortech", "x2", "sectorsvc:x2", "sectortech:x2")
  )
  expect_identical(
     out_des$colnames_Z,
    c("log(x1)", "sectorsvc", "sectortech", "x2", "sectorsvc:x2", "sectortech:x2")
  )
  expect_identical(trt_des$colnames_Z, c("log(x1)", "x2"))
  expect_identical(obs_des$colnames_Z, c("log(x1)", "x2", "w_sel", "sectorsvc", "sectortech"))

  req_idx <- which(cell$R == 1L)
  irr_idx <- which(cell$R == 0L)
  expect_identical(req_idx, c(1L, 2L, 6L, 9L, 10L))
  expect_identical(irr_idx, c(3L, 4L, 5L, 7L, 8L))

  # Required rows are finite; irrelevant rows stay NA without having evaluated log(-999).
  expect_false(anyNA(out_des$X_tilde[req_idx, , drop = FALSE]))
  expect_true(all(is.na(out_des$X_tilde[irr_idx, , drop = FALSE])))
  expect_true(all(out_des$X_tilde[req_idx, "sectortech"] == 0.0))
  expect_true(all(out_des$X_tilde[req_idx, "sectortech:x2"] == 0.0))
  expect_identical(out_des$diagnostics$irrelevant_rows_evaluated, FALSE)
  expect_identical(out_des$diagnostics$rank_preprocessed, FALSE)

  # Now test needed X missingness: if a unit in R_c (say u02, which has 0 outcomes!) has NA in x1 at t_prev = 2012,
  # design extraction on R_c must fail with `cdid_static_missing_covariate_error` (never silently dropped via na.omit!).
  df_need_na <- df_x
  df_need_na$x1[df_need_na$unit == "u02" & df_need_na$year == 2012] <- NA_real_
  panel_need_na <- static_panel_data(df_need_na, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  cell_need_na <- static_transition_cell(panel_need_na, t = 2015, d = 0)
  err_cov <- expect_error(
    static_reference_design(panel_need_na, cell_need_na, formula = ~ x1 + x2, population = "risk"),
    class = "cdid_static_missing_covariate_error"
  )
  expect_identical(err_cov$missing_units, "u02")
  expect_identical(err_cov$variable, "x1")

  # However, if we extract only on the switcher role `population = "switcher"` (G_c = u06, u09),
  # where u02 is NOT in G_c, then u02's missing x1 is irrelevant and extraction on G_c succeeds!
  des_sw <- static_reference_design(panel_need_na, cell_need_na, formula = ~ x1 + x2, population = "switcher")
  expect_identical(des_sw$required_idx, c(6L, 9L))
  expect_false(anyNA(des_sw$X_tilde[c(6L, 9L), , drop = FALSE]))
  expect_true(all(is.na(des_sw$X_tilde[-c(6L, 9L), , drop = FALSE])))

  # Refuse no-intercept formulas (~ 0 + x1, ~ x1 - 1, ~ -1 + x1).
  expect_error(
    static_reference_design(panel, cell, formula = ~ 0 + x1 + x2),
    class = "cdid_static_no_intercept_error"
  )
  expect_error(
    static_reference_design(panel, cell, formula = ~ x1 - 1),
    class = "cdid_static_no_intercept_error"
  )
})

test_that("observation-pair access returns exact O, dy, counts, index gaps, retains zero-count pairs, and supports lazy/chunked callbacks", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)
  panel <- static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  c_2015_1 <- static_transition_cell(panel, t = 2015, d = 1)

  pair_df <- data.frame(
    s = c(2012, 2010, 2012, 2010),
    u = c(2015, 2015, 2019, 2025)
  )
  pairs_res <- static_observation_pairs(
    panel = panel,
    mask = c_2015_1$G,
    pairs = pair_df
  )

  expect_s3_class(pairs_res, "cdid_static_observation_pairs")
  expect_identical(pairs_res$n_pairs, 4L)
  expect_identical(pairs_res$unit_id, ref$unit_id)
  expect_identical(pairs_res$mask, c_2015_1$G)
  expect_identical(pairs_res$empty_pair_ids, "(2012,2015)")
  expect_identical(pairs_res$supported_pair_ids, c("(2010,2015)", "(2012,2019)", "(2010,2025)"))

  for (pid in names(ref$pairs_c2015_d1_G)) {
    got_p <- pairs_res$pairs[[pid]]
    exp_p <- ref$pairs_c2015_d1_G[[pid]]
    expect_identical(got_p$s, exp_p$s)
    expect_identical(got_p$u, exp_p$u)
    expect_identical(got_p$s_idx, exp_p$s_idx)
    expect_identical(got_p$u_idx, exp_p$u_idx)
    expect_identical(got_p$gap, exp_p$gap)
    expect_identical(got_p$L_k, exp_p$gap)
    expect_identical(got_p$calendar_span, exp_p$calendar_span)
    expect_identical(got_p$O, exp_p$O)
    expect_identical(got_p$A_O, exp_p$A_O)
    expect_equal(got_p$dy, exp_p$dy)
    expect_identical(got_p$n_role, exp_p$n_role)
    expect_identical(got_p$n_observed_all, exp_p$n_observed_all)
    expect_identical(got_p$n_observed_role, exp_p$n_observed_role)
    expect_identical(got_p$supported_in_role, exp_p$supported_in_role)

    # Non-observed dy rows are strictly NA_real_ (not NaN, not 0).
    expect_true(all(is.na(got_p$dy[got_p$O == 0L])))
    expect_false(any(is.nan(got_p$dy)))
    expect_true(all(is.finite(got_p$dy[got_p$O == 1L])))
  }

  # Lazy accessor and callback / chunked streaming avoid materializing N-by-all-pairs at once.
  acc <- static_observation_pairs(panel, mask = c_2015_1$G, pairs = pair_df, materialize = FALSE)
  expect_s3_class(acc, "cdid_static_pair_accessor")
  expect_identical(acc$get_by_dates(2012, 2019)$dy, ref$pairs_c2015_d1_G[["(2012,2019)"]]$dy)

  stream_single <- static_observation_pairs(
    panel = panel,
    mask = c_2015_1$G,
    pairs = pair_df,
    callback = function(p, idx) {
      list(idx = idx, gap = p$gap, n_obs_role = p$n_observed_role)
    }
  )
  expect_s3_class(stream_single, "cdid_static_observation_pairs_stream")
  expect_identical(
    vapply(stream_single$callback_results, `[[`, integer(1L), "gap"),
    c("(2012,2015)" = 1L, "(2010,2015)" = 2L, "(2012,2019)" = 2L, "(2010,2025)" = 4L)
  )

  stream_chunked <- static_observation_pairs(
    panel = panel,
    mask = c_2015_1$G,
    pairs = pair_df,
    chunk_size = 2L,
    callback = function(chunk, indices) {
      list(indices = indices, ids = names(chunk))
    }
  )
  expect_identical(length(stream_chunked$callback_results), 2L)
  expect_identical(stream_chunked$callback_results[[1L]]$indices, c(1L, 2L))
  expect_identical(stream_chunked$callback_results[[2L]]$indices, c(3L, 4L))

  # Invalid pair endpoints and ordering (s >= u) are refused.
  expect_error(
    static_observation_pair(panel, c_2015_1$G, s = 2015, u = 2015),
    class = "cdid_static_pair_order_error"
  )
  expect_error(
    static_observation_pair(panel, c_2015_1$G, s = 2019, u = 2012),
    class = "cdid_static_pair_order_error"
  )
  expect_error(
    static_observation_pair(panel, c_2015_1$G, s = 2011, u = 2015),
    class = "cdid_static_pair_endpoint_error"
  )
})

test_that("additional contract edge cases: gname refusal, character/inline factors, custom mask row checks, and pair formats", {
  ref <- dget(testthat::test_path("fixtures", "static_panel_reference_fixture.R"))
  df_long <- build_hand_5date_long_df(ref)

  # Passing `gname` to `static_panel_data` is rejected (no `gname` formal exists).
  expect_error(
    static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated", gname = "first_treat")
  )

  # Factor or NA/empty unit IDs and single-date panels are rejected.
  df_fac_id <- df_long
  df_fac_id$unit <- factor(df_fac_id$unit)
  expect_error(
    static_panel_data(df_fac_id, idname = "unit", tname = "year", yname = "outcome", dname = "treated"),
    class = "cdid_static_panel_id_error"
  )
  df_one_date <- df_long[df_long$year == 2010, , drop = FALSE]
  expect_error(
    static_panel_data(df_one_date, idname = "unit", tname = "year", yname = "outcome", dname = "treated"),
    class = "cdid_static_panel_date_error"
  )

  # Foreign unit IDs not in explicit `unit_ids` are rejected.
  expect_error(
    static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated",
                      unit_ids = ref$unit_id[1:9]),
    class = "cdid_static_panel_id_error"
  )

  # Absent row at 2012 on u07: transitions touching 2012 (t = 2012 or t = 2015) fail,
  # while genuinely compatible transition c = (2019, 0) (2015 -> 2019) succeeds and can be reused.
  panel <- static_panel_data(df_long, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  df_abs_2012 <- df_long[!(df_long$unit == "u07" & df_long$year == 2012), , drop = FALSE]
  p_abs <- static_panel_data(df_abs_2012, idname = "unit", tname = "year", yname = "outcome", dname = "treated")
  expect_error(
    static_transition_cell(p_abs, t = 2012, d = 0),
    class = "cdid_static_missing_treatment_error"
  )
  c_2019_0 <- static_transition_cell(p_abs, t = 2019, d = 0)
  c_2019_0_full <- static_transition_cell(panel, t = 2019, d = 0)
  c_2015_0_full <- static_transition_cell(panel, t = 2015, d = 0)

  mask_without_u07 <- c(1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L)
  des_ok <- static_reference_design(p_abs, c_2019_0_full, formula = ~ x1 + sector_chr, mask = mask_without_u07)
  expect_identical(des_ok$colnames_Z, c("x1", "sector_chrsvc", "sector_chrtech"))

  err_row <- expect_error(
    static_reference_design(p_abs, c_2015_0_full, formula = ~ x1),
    class = "cdid_static_missing_treatment_error"
  )
  expect_identical(err_row$missing_units, "u07")

  #Mismatched unit names on `mask` in `static_observation_pairs` or `static_reference_design` are refused.
  bad_named_mask <- stats::setNames(c_2019_0$G, rev(panel$unit_id))
  expect_error(
    static_observation_pair(panel, bad_named_mask, s = 2010, u = 2015),
    class = "cdid_static_pair_mask_error"
  )
  expect_error(
    static_reference_design(panel, c_2019_0, formula = ~ x1, mask = bad_named_mask),
    class = "cdid_static_design_error"
  )

  # Matrix and list pair specifications and duplicate pair rejection.
  p_mat <- matrix(c(2010, 2012, 2015, 2019), nrow = 2L, ncol = 2L, dimnames = list(c("p1", "p2"), c("s", "u")))
  res_mat <- static_observation_pairs(panel, c_2019_0$C, pairs = p_mat)
  expect_identical(res_mat$pair_ids, c("p1", "p2"))
  expect_identical(res_mat$summary_table$gap, c(2L, 2L))

  p_list <- list(e1 = c(2010, 2012), e2 = list(s = 2012, u = 2025))
  res_list <- static_observation_pairs(panel, c_2019_0$C, pairs = p_list)
  expect_identical(res_list$pair_ids, c("e1", "e2"))
  expect_identical(res_list$summary_table$gap, c(1L, 3L))

  expect_error(
    static_observation_pairs(panel, c_2019_0$C, s = c(2010, 2010), u = c(2015, 2015)),
    class = "cdid_static_pair_input_error"
  )

  # Reused-cell context validation (P1): reject stale date indices, masks, and missing D on both reuse paths.
  a_ctx <- data.frame(
    id = rep(c("A", "B"), each = 3L),
    t = rep(c(10, 40, 100), 2L),
    D = c(0, 0, 1, 0, 0, 0),
    Y = 1:6,
    x = c(11, 44, 101, 22, 88, 202),
    stringsAsFactors = FALSE
  )
  b_ctx <- rbind(
    a_ctx,
    data.frame(id = c("A", "B"), t = c(30, 30), D = c(0, 0), Y = c(7, 8), x = c(333, 666), stringsAsFactors = FALSE)
  )
  pa_ctx <- static_panel_data(a_ctx, "id", "t", "Y", "D")
  pb_ctx <- static_panel_data(b_ctx, "id", "t", "Y", "D")
  ca_ctx <- static_transition_cell(pa_ctx, 100, 0)
  cb_ctx <- static_transition_cell(pb_ctx, 100, 0)

  expect_error(
    static_transition_cells(pb_ctx, cells = list(ca_ctx)),
    class = "cdid_static_cell_error"
  )
  expect_error(
    static_reference_design(pb_ctx, ca_ctx, formula = ~ x),
    class = "cdid_static_cell_error"
  )
  ext_pb <- static_reference_design_extractor(pb_ctx, outcome_formula = ~ x)
  expect_error(
    ext_pb(ca_ctx),
    class = "cdid_static_cell_error"
  )
  des_cb <- static_reference_design(pb_ctx, cb_ctx, formula = ~ x)
  des_fresh <- static_reference_design(pb_ctx, list(t = 100, d = 0), formula = ~ x)
  expect_identical(unname(des_cb$Z[, 1L]), c(44, 88))
  expect_identical(unname(des_cb$Z), unname(des_fresh$Z))

  inc_ctx <- a_ctx
  inc_ctx$D[inc_ctx$id == "B" & inc_ctx$t == 100] <- NA_real_
  pm_ctx <- static_panel_data(inc_ctx, "id", "t", "Y", "D")
  expect_error(
    static_transition_cell(pm_ctx, 100, 0),
    class = "cdid_static_missing_treatment_error"
  )
  expect_error(
    static_transition_cells(pm_ctx, cells = list(ca_ctx)),
    class = "cdid_static_missing_treatment_error"
  )
  expect_error(
    static_reference_design(pm_ctx, ca_ctx, formula = ~ x),
    class = "cdid_static_missing_treatment_error"
  )
  ext_pm <- static_reference_design_extractor(pm_ctx, outcome_formula = ~ x)
  expect_error(
    ext_pm(ca_ctx),
    class = "cdid_static_missing_treatment_error"
  )

  # Exact numeric unit identity collision (P2): 1e15 vs 1e15 + 1 share printed unit_key "1e+15".
  a_num <- data.frame(
    id = rep(c(1e15, 2), each = 2L),
    t = rep(c(10, 40), 2L),
    Y = 1:4,
    D = c(0, 1, 0, 0),
    x = c(11, 12, 21, 22)
  )
  b_num <- a_num
  b_num$id[b_num$id == 1e15] <- 1e15 + 1
  b_num$D <- c(0, 0, 0, 1)
  pa_num <- static_panel_data(a_num, "id", "t", "Y", "D")
  pb_num <- static_panel_data(b_num, "id", "t", "Y", "D")
  ca_num <- static_transition_cell(pa_num, 40, 0)
  expect_false(identical(pa_num$unit_id, pb_num$unit_id))
  expect_identical(pa_num$unit_key, pb_num$unit_key)
  expect_error(
    static_transition_cells(pb_num, cells = list(ca_num)),
    class = "cdid_static_panel_id_error"
  )
  expect_error(
    static_reference_design(pb_num, ca_num, formula = ~ x, population = "switcher"),
    class = "cdid_static_panel_id_error"
  )
  ext_pb_num <- static_reference_design_extractor(pb_num, outcome_formula = ~ x)
  expect_error(
    ext_pb_num(ca_num, population = "switcher"),
    class = "cdid_static_panel_id_error"
  )
})

