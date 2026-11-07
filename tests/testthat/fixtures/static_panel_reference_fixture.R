list(
  dates = c(2010, 2012, 2015, 2019, 2025),
  unit_id = c("u01", "u02", "u03", "u04", "u05", "u06", "u07", "u08", "u09", "u10"),
  D_matrix = matrix(
    as.integer(c(
      # 2010, 2012, 2015, 2019, 2025
      0L, 0L, 0L, 0L, 0L,        # u01: constant untreated
      0L, 0L, 0L, 0L, 0L,        # u02: constant untreated, zero outcomes
      0L, 1L, 1L, 1L, 1L,        # u03: first-entry 2012, absorbing
      0L, 1L, 0L, 1L, 1L,        # u04: first-entry 2012, exit 2015, re-entry 2019
      0L, 1L, 0L, 0L, 1L,        # u05: first-entry 2012, exit 2015, re-entry 2025 (single outcome)
      1L, 0L, 1L, 0L, 0L,        # u06: initial treated, exit 2012, re-entry 2015, exit 2019
      1L, 1L, 1L, 1L, 1L,        # u07: constant always-treated
      1L, 1L, 0L, 0L, 0L,        # u08: initial treated, exit 2015
      0L, 0L, 1L, 1L, 0L,        # u09: late entry 2015, exit 2025
      0L, 0L, 0L, 0L, NA_integer_ # u10: untreated through 2019, unknown D at 2025
    )),
    nrow = 10L,
    ncol = 5L,
    byrow = TRUE,
    dimnames = list(
      unit_id = c("u01", "u02", "u03", "u04", "u05", "u06", "u07", "u08", "u09", "u10"),
      date = c("2010", "2012", "2015", "2019", "2025")
    )
  ),
  Y_matrix = matrix(
    c(
      # 2010, 2012, 2015, 2019, 2025
      10, 12, 15, 19, 25,        # u01: 5 outcomes
      NA, NA, NA, NA, NA,        # u02: 0 outcomes
      20, 25, NA, 32, 40,        # u03: 4 outcomes
      30, NA, 33, 41, 46,        # u04: 4 outcomes
      NA, NA, 50, NA, NA,        # u05: 1 outcome
      40, 38, 45, 42, NA,        # u06: 4 outcomes
      60, 63, 67, NA, 75,        # u07: 4 outcomes
      70, 74, NA, 71, NA,        # u08: 3 outcomes
      15, 18, 26, 30, 27,        # u09: 5 outcomes
      80, NA, NA, 88, 92         # u10: 3 outcomes
    ),
    nrow = 10L,
    ncol = 5L,
    byrow = TRUE,
    dimnames = list(
      unit_id = c("u01", "u02", "u03", "u04", "u05", "u06", "u07", "u08", "u09", "u10"),
      date = c("2010", "2012", "2015", "2019", "2025")
    )
  ),
  unit_outcome_counts = c(
    u01 = 5L, u02 = 0L, u03 = 4L, u04 = 4L, u05 = 1L,
    u06 = 4L, u07 = 4L, u08 = 3L, u09 = 5L, u10 = 3L
  ),
  cells = list(
    "t=2012,d=0" = list(
      t = 2012,
      t_prev = 2010,
      t_idx = 2L,
      t_prev_idx = 1L,
      d = 0L,
      sigma = 1.0,
      R = c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 1L),
      G = c(0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L),
      C = c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L),
      n_risk = 7L,
      n_switchers = 3L,
      n_stayers = 4L,
      status = "supported"
    ),
    "t=2012,d=1" = list(
      t = 2012,
      t_prev = 2010,
      t_idx = 2L,
      t_prev_idx = 1L,
      d = 1L,
      sigma = -1.0,
      R = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 0L),
      G = c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L),
      C = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L),
      n_risk = 3L,
      n_switchers = 1L,
      n_stayers = 2L,
      status = "supported"
    ),
    "t=2015,d=0" = list(
      t = 2015,
      t_prev = 2012,
      t_idx = 3L,
      t_prev_idx = 2L,
      d = 0L,
      sigma = 1.0,
      R = c(1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L),
      G = c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L),
      C = c(1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L),
      n_risk = 5L,
      n_switchers = 2L,
      n_stayers = 3L,
      status = "supported"
    ),
    "t=2015,d=1" = list(
      t = 2015,
      t_prev = 2012,
      t_idx = 3L,
      t_prev_idx = 2L,
      d = 1L,
      sigma = -1.0,
      R = c(0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L),
      G = c(0L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 0L),
      C = c(0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L),
      n_risk = 5L,
      n_switchers = 3L,
      n_stayers = 2L,
      status = "supported"
    ),
    "t=2019,d=0" = list(
      t = 2019,
      t_prev = 2015,
      t_idx = 4L,
      t_prev_idx = 3L,
      d = 0L,
      sigma = 1.0,
      R = c(1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L),
      G = c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L),
      C = c(1L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 1L),
      n_risk = 6L,
      n_switchers = 1L,
      n_stayers = 5L,
      status = "supported"
    ),
    "t=2019,d=1" = list(
      t = 2019,
      t_prev = 2015,
      t_idx = 4L,
      t_prev_idx = 3L,
      d = 1L,
      sigma = -1.0,
      R = c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 0L),
      G = c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L),
      C = c(0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 0L),
      n_risk = 4L,
      n_switchers = 1L,
      n_stayers = 3L,
      status = "supported"
    )
  ),
  # Reference pairs for cell c = (2015, 1) switcher mask G = c(0, 0, 0, 1, 1, 0, 0, 1, 0, 0)
  pairs_c2015_d1_G = list(
    "(2012,2015)" = list(
      s = 2012, u = 2015, s_idx = 2L, u_idx = 3L, gap = 1L, calendar_span = 3.0,
      O = c(1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L),
      A_O = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
      dy = c(3, NA, NA, NA, NA, 7, 4, NA, 8, NA),
      n_role = 3L, n_observed_all = 4L, n_observed_role = 0L,
      supported_in_role = FALSE
    ),
    "(2010,2015)" = list(
      s = 2010, u = 2015, s_idx = 1L, u_idx = 3L, gap = 2L, calendar_span = 5.0,
      O = c(1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L),
      A_O = c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L),
      dy = c(5, NA, NA, 3, NA, 5, 7, NA, 11, NA),
      n_role = 3L, n_observed_all = 5L, n_observed_role = 1L,
      supported_in_role = TRUE
    ),
    "(2012,2019)" = list(
      s = 2012, u = 2019, s_idx = 2L, u_idx = 4L, gap = 2L, calendar_span = 7.0,
      O = c(1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L),
      A_O = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L),
      dy = c(7, NA, 7, NA, NA, 4, NA, -3, 12, NA),
      n_role = 3L, n_observed_all = 5L, n_observed_role = 1L,
      supported_in_role = TRUE
    ),
    "(2010,2025)" = list(
      s = 2010, u = 2025, s_idx = 1L, u_idx = 5L, gap = 4L, calendar_span = 15.0,
      O = c(1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L),
      A_O = c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L),
      dy = c(15, NA, 20, 16, NA, NA, 15, NA, 12, 12),
      n_role = 3L, n_observed_all = 6L, n_observed_role = 1L,
      supported_in_role = TRUE
    )
  )
)
