# Historical hand tables use explicit bounded debug materialization.
test_make_pairs <- function(...) cdid:::materialize_observed_pairs(
  cdid:::construct_observed_pairs(...), max_pairs = 100000, chunk_size = 2)

test_that("construct_observed_pairs enumerates all within-unit s < t pairs and preserves singletons", {
  # Hand-constructed unbalanced panel with irregular finite numeric dates,
  # gapped/reentry observations, an explicit NA row, a singleton, and an
  # all-NA unobserved unit.
  raw <- data.frame(
    unit = c("u2", "u2", "u2", "u2",
             "u1", "u1", "u1",
             "u3", "u3",
             "u4"),
    date = c(10.5, 2.0, 5.25, 1.0,
             1.0, 5.25, 10.5,
             2.0, 5.25,
             2.0),
    y = c(14.0, 3.5, NA_real_, 1.5,
          10.0, 16.0, 25.0,
          7.25, NA_real_,
          NA_real_),
    cohort = c(5.25, 5.25, 5.25, 5.25,
               0, 0, 0,
               10.5, 10.5,
               Inf),
    stringsAsFactors = FALSE
  )
  raw_copy <- raw

  res <- test_make_pairs(
    data = raw, idname = "unit", tname = "date", yname = "y", gname = "cohort"
  )

  # Input data must remain completely unchanged.
  expect_identical(raw, raw_copy)
  expect_s3_class(res, "cdid_observed_pairs")

  # Canonical unit order is sorted by id: u1, u2, u3, u4.
  expect_identical(res$units$unit_index, 1:4)
  expect_identical(res$units$id, c("u1", "u2", "u3", "u4"))
  expect_identical(res$units$g, c(0, 5.25, 10.5, Inf))
  expect_identical(res$units$n_rows, c(3L, 4L, 2L, 1L))
  expect_identical(res$units$n_obs, c(3L, 3L, 1L, 0L))
  expect_identical(res$units$n_pairs, c(3, 3, 0, 0))
  expect_identical(res$units$is_unobserved, c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(res$units$is_singleton, c(FALSE, FALSE, TRUE, FALSE))
  expect_identical(res$units$pair_eligible, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(res$units$singleton_time, c(NA_real_, NA_real_, 2.0, NA_real_))
  expect_identical(res$units$singleton_y, c(NA_real_, NA_real_, 7.25, NA_real_))

  # Singleton table retains u3 with its observed time and outcome.
  expect_identical(nrow(res$singletons), 1L)
  expect_identical(res$singletons$id, "u3")
  expect_identical(res$singletons$time, 2.0)
  expect_identical(res$singletons$y, 7.25)

  # Global observed dates are (1.0, 2.0, 5.25, 10.5).
  expect_identical(res$times, c(1.0, 2.0, 5.25, 10.5))

  # Hand-verified pairs for u1 (observed at 1.0, 5.25, 10.5) and
  # u2 (observed at 1.0, 2.0, 10.5; NA at 5.25).
  expect_identical(nrow(res$pairs), 6L)
  expect_identical(res$pairs$id, c("u1", "u1", "u1", "u2", "u2", "u2"))
  expect_identical(res$pairs$s, c(1.0, 1.0, 5.25, 1.0, 1.0, 2.0))
  expect_identical(res$pairs$t, c(5.25, 10.5, 10.5, 2.0, 10.5, 10.5))
  expect_identical(res$pairs$y_s, c(10.0, 10.0, 16.0, 1.5, 1.5, 3.5))
  expect_identical(res$pairs$y_t, c(16.0, 25.0, 25.0, 3.5, 14.0, 14.0))
  expect_identical(res$pairs$dy, c(6.0, 15.0, 9.0, 2.0, 12.5, 10.5))
  expect_identical(res$pairs$time_gap, c(4.25, 9.5, 5.25, 1.0, 9.5, 8.5))
  expect_identical(res$pairs$time_gap_finite, rep(TRUE, 6L))
  expect_identical(res$pairs$n_intervening_obs, c(0L, 1L, 0L, 0L, 1L, 0L))
  expect_identical(res$pairs$n_intervening_panel_dates, c(1L, 2L, 0L, 0L, 2L, 1L))
  expect_identical(res$pairs$is_consecutive_unit, c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE))
  expect_identical(res$pairs$is_consecutive_panel, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))

  # Absence of a row vs an explicit NA outcome row produces identical pairs,
  # singletons, observed times/outcomes, and unit observation metadata.
  raw_no_na_u2 <- raw[!(raw$unit == "u2" & raw$date == 5.25), ]
  res_no_na_u2 <- test_make_pairs(
    data = raw_no_na_u2, idname = "unit", tname = "date",
    yname = "y", gname = "cohort"
  )
  expect_identical(res$pairs, res_no_na_u2$pairs)
  expect_identical(res$singletons, res_no_na_u2$singletons)
  expect_identical(res$observations, res_no_na_u2$observations)
  unit_obs_cols <- setdiff(names(res$units), "n_rows")
  expect_identical(res$units[, unit_obs_cols], res_no_na_u2$units[, unit_obs_cols])
})

test_that("construct_observed_pairs is strictly invariant to row permutations", {
  dat <- data.frame(
    id = c(30L, 10L, 20L, 30L, 10L, 20L, 10L, 40L),
    time = c(2005, 2000, 2002, 2000, 2009, 2005, 2002, 2002),
    y = c(8.5, 1.0, 4.0, 2.5, 11.0, 7.0, NA_real_, 99.0),
    g = c(2005, 2009, 0, 2005, 2009, 0, 2009, Inf)
  )
  base_res <- test_make_pairs(dat, "id", "time", "y", "g")

  perms <- list(
    rev(seq_len(nrow(dat))),
    c(2L, 5L, 7L, 3L, 6L, 1L, 4L, 8L),
    c(8L, 1L, 6L, 2L, 4L, 7L, 3L, 5L)
  )
  for (p in perms) {
    perm_res <- test_make_pairs(dat[p, ], "id", "time", "y", "g")
    expect_identical(perm_res, base_res)
  }
})

test_that("construct_observed_pairs handles nonadjacent two-outcome units and panels with no fully observed unit", {
  # Four panel dates (-2.5, 0.0, 3.5, 8.0), but no unit is observed at all 4 dates.
  # Unit 1 has exactly two nonadjacent outcomes (-2.5 and 8.0).
  dat <- data.frame(
    id = c(1L, 1L,
           2L, 2L, 2L,
           3L, 3L,
           4L),
    t = c(-2.5, 8.0,
          -2.5, 0.0, 8.0,
          0.0, 3.5,
          3.5),
    y = c(1.25, 6.75,
          2.0, 3.0, 9.0,
          4.5, 2.5,
          10.0),
    g = c(3.5, 3.5,
          8.0, 8.0, 8.0,
          0.0, 0.0,
          3.5)
  )
  res <- test_make_pairs(dat, "id", "t", "y", "g")

  expect_identical(res$times, c(-2.5, 0.0, 3.5, 8.0))
  expect_true(all(res$units$n_obs < length(res$times)))

  # Unit 1 has exactly 2 nonadjacent outcomes and 1 pair spanning 2 intervening panel dates.
  u1_pair <- res$pairs[res$pairs$id == 1L, ]
  expect_identical(nrow(u1_pair), 1L)
  expect_identical(u1_pair$s, -2.5)
  expect_identical(u1_pair$t, 8.0)
  expect_identical(u1_pair$dy, 5.5)
  expect_identical(u1_pair$n_intervening_obs, 0L)
  expect_identical(u1_pair$n_intervening_panel_dates, 2L)
  expect_true(u1_pair$is_consecutive_unit)
  expect_false(u1_pair$is_consecutive_panel)

  # Total pairs: 1 (unit 1) + 3 (unit 2) + 1 (unit 3) + 0 (unit 4 singleton) = 5.
  expect_identical(res$counts$n_pairs, 5)
  expect_identical(res$counts$n_singleton_units, 1L)
  expect_identical(res$counts$n_pair_eligible_units, 3L)
})

test_that("construct_observed_pairs rejects invalid duplicates, cohort changes, missing keys, and nonfinite outcomes", {
  base_dat <- data.frame(
    id = c(1L, 1L, 2L, 2L),
    time = c(1, 2, 1, 2),
    y = c(1.0, 2.0, 3.0, 4.0),
    g = c(2, 2, 0, 0)
  )

  # Duplicate (id, time) with observed outcomes.
  dup_obs <- rbind(base_dat, data.frame(id = 1L, time = 2, y = 2.0, g = 2))
  expect_error(
    test_make_pairs(dup_obs, "id", "time", "y", "g"),
    class = "cdid_pairs_duplicate_error"
  )

  # Duplicate (id, time) where one duplicate has NA outcome is still invalid.
  dup_na <- rbind(base_dat, data.frame(id = 2L, time = 1, y = NA_real_, g = 0))
  expect_error(
    test_make_pairs(dup_na, "id", "time", "y", "g"),
    class = "cdid_pairs_duplicate_error"
  )

  # Within-unit cohort change (including on an NA-outcome row).
  bad_g <- base_dat
  bad_g$g[2L] <- 3
  expect_error(
    test_make_pairs(bad_g, "id", "time", "y", "g"),
    class = "cdid_pairs_cohort_change_error"
  )
  bad_g_na <- rbind(base_dat, data.frame(id = 1L, time = 3, y = NA_real_, g = 4))
  expect_error(
    test_make_pairs(bad_g_na, "id", "time", "y", "g"),
    class = "cdid_pairs_cohort_change_error"
  )

  # Non-finite observed outcomes (Inf, -Inf, NaN) are rejected.
  for (bad_val in c(Inf, -Inf, NaN)) {
    d_bad_y <- base_dat
    d_bad_y$y[1L] <- bad_val
    expect_error(
      test_make_pairs(d_bad_y, "id", "time", "y", "g"),
      class = "cdid_pairs_outcome_error"
    )
  }

  # Missing/nonfinite ID, time, or cohort are rejected.
  d_bad_id <- base_dat
  d_bad_id$id[1L] <- NA_integer_
  expect_error(
    test_make_pairs(d_bad_id, "id", "time", "y", "g"),
    class = "cdid_pairs_id_error"
  )

  for (bad_t in c(NA_real_, NaN, Inf, -Inf)) {
    d_bad_t <- base_dat
    d_bad_t$time[1L] <- bad_t
    expect_error(
      test_make_pairs(d_bad_t, "id", "time", "y", "g"),
      class = "cdid_pairs_time_error"
    )
  }

  for (bad_cohort in c(NA_real_, NaN, -Inf)) {
    d_bad_c <- base_dat
    d_bad_c$g[1L] <- bad_cohort
    d_bad_c$g[2L] <- bad_cohort
    expect_error(
      test_make_pairs(d_bad_c, "id", "time", "y", "g"),
      class = "cdid_pairs_cohort_error"
    )
  }
})

test_that("extreme finite numeric dates preserve outcome differences and calendar order with honest time_gap_finite metadata", {
  d_two <- data.frame(id = c(1L, 1L), time = c(-1e308, 1e308), y = c(0, 1), g = c(0, 0))
  idx_two <- cdid:::construct_observed_pairs(d_two, "id", "time", "y", "g")
  batches <- list()
  cdid:::traverse_observed_pairs(idx_two, 1, function(b) batches[[length(batches) + 1L]] <<- b)
  expect_identical(length(batches), 1L)
  expect_identical(batches[[1L]]$s, -1e308)
  expect_identical(batches[[1L]]$t, 1e308)
  expect_identical(batches[[1L]]$dy, 1)
  expect_identical(batches[[1L]]$time_gap, Inf)
  expect_identical(batches[[1L]]$time_gap_finite, FALSE)
  expect_identical(batches[[1L]]$n_intervening_obs, 0L)
  expect_identical(batches[[1L]]$n_intervening_panel_dates, 0L)
  expect_true(batches[[1L]]$is_consecutive_unit)
  expect_true(batches[[1L]]$is_consecutive_panel)
  expect_identical(cdid:::empty_pair_batch(idx_two, 0L), batches[[1L]][FALSE, ])

  d_three <- data.frame(
    id = c("u1", "u1", "u1", "u2"),
    time = c(-1e308, 0, 1e308, 0),
    y = c(2, 5, 11, NA_real_),
    g = c(0, 0, 0, Inf),
    stringsAsFactors = FALSE
  )
  res_three <- test_make_pairs(d_three, "id", "time", "y", "g")
  expect_identical(res_three$pairs$s, c(-1e308, -1e308, 0))
  expect_identical(res_three$pairs$t, c(0, 1e308, 1e308))
  expect_identical(res_three$pairs$dy, c(3, 9, 6))
  expect_identical(res_three$pairs$time_gap, c(1e308, Inf, 1e308))
  expect_identical(res_three$pairs$time_gap_finite, c(TRUE, FALSE, TRUE))
  expect_identical(res_three$pairs$n_intervening_obs, c(0L, 1L, 0L))
  expect_identical(res_three$pairs$n_intervening_panel_dates, c(0L, 1L, 0L))
  expect_identical(cdid:::empty_pair_batch(res_three, 0L), res_three$pairs[FALSE, ])
})

