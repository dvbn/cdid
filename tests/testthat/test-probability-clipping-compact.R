# Closed deterministic tests for compact probability-clipping diagnostic summaries
# (`ipw_clip_summary` and `ipw_clip_summary_merge`). No fits, RNG, or production
# helper oracles are used.

resolve_clip_helper <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) {
    return(get(name, mode = "function", inherits = TRUE))
  }
  if ("cdid" %in% loadedNamespaces()) {
    ns <- asNamespace("cdid")
    if (exists(name, envir = ns, mode = "function", inherits = FALSE)) {
      return(get(name, envir = ns, mode = "function", inherits = FALSE))
    }
  }
  stop("Required helper `", name, "` is not available.", call. = FALSE)
}

assert_compact_scalar_list <- function(x, expected_names) {
  expect_true(is.list(x))
  expect_false(is.object(x))
  expect_identical(names(x), expected_names)
  expect_identical(length(x), length(expected_names))
  atomic_lengths <- unlist(rapply(x, length, how = "unlist"))
  expect_identical(max(atomic_lengths), 1L)
  expect_true(all(atomic_lengths == 1L))
  for (nm in expected_names) {
    val <- x[[nm]]
    expect_true(is.atomic(val) && is.double(val) && is.null(dim(val)) && !is.object(val))
    expect_identical(length(val), 1L)
    expect_null(names(val))
  }
}

test_that("ipw_clip_summary returns compact scalar lists with exact values across all field sets including 13 and 18 fields", {
  clip_summary <- resolve_clip_helper("ipw_clip_summary")

  eps <- 1e-6
  cap <- 1 - 1e-6

  # Deterministic inputs exercising interior, exact-boundary, and beyond-boundary values
  prob <- c(0.125, 0.5, cap, 0.9999995, 0.9999999)
  ctrl <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  q <- c(2e-7, eps, 0.04, 0.85)
  rates <- c(5e-7, eps, 0.2, 0.75)
  weights <- c(0.75, 1.5, 12.25, 3.0)

  p_fields <- c("p_min", "p_max", "p_clip_max", "p_n", "p_capped", "p_capped_ctrl")
  q_fields <- c("q_min", "q_max", "q_floor_min", "q_n", "q_floored")
  r_fields <- c("r_min", "r_max", "r_floor_min", "r_n", "r_floored")
  w_fields <- c("w_max")

  exp_p_ctrl <- list(
    p_min = 0.125,
    p_max = 0.9999999,
    p_clip_max = cap,
    p_n = 5.0,
    p_capped = 3.0,
    p_capped_ctrl = 2.0
  )
  exp_p_no_ctrl <- modifyList(exp_p_ctrl, list(p_capped_ctrl = NA_real_))
  exp_q <- list(
    q_min = 2e-7,
    q_max = 0.85,
    q_floor_min = eps,
    q_n = 4.0,
    q_floored = 2.0
  )
  exp_r <- list(
    r_min = 5e-7,
    r_max = 0.75,
    r_floor_min = eps,
    r_n = 4.0,
    r_floored = 2.0
  )
  exp_w <- list(w_max = 12.25)

  # 1. Default 1-field set (eps only)
  s1 <- clip_summary()
  assert_compact_scalar_list(s1, "eps")
  expect_identical(s1, list(eps = eps))

  # 2. 7-field set (q + weights: OR Corollary 1 and augmented cell consumers)
  s7_qw <- clip_summary(q = q, weights = weights)
  assert_compact_scalar_list(s7_qw, c("eps", q_fields, w_fields))
  expect_identical(s7_qw, c(list(eps = eps), exp_q, exp_w))
  expect_identical(length(s7_qw), 7L)

  # 3. 13-field set A (prob + rates + controls + weights: IPW link / fitted view / point view)
  s13_pr <- clip_summary(prob = prob, rates = rates, controls = ctrl, weights = weights)
  assert_compact_scalar_list(s13_pr, c("eps", p_fields, r_fields, w_fields))
  expect_identical(s13_pr, c(list(eps = eps), exp_p_ctrl, exp_r, exp_w))
  expect_identical(length(s13_pr), 13L)

  # 4. 13-field set B (prob + q + controls + weights: CP04 selection treatment fit)
  s13_pq <- clip_summary(prob = prob, q = q, controls = ctrl, weights = weights)
  assert_compact_scalar_list(s13_pq, c("eps", p_fields, q_fields, w_fields))
  expect_identical(s13_pq, c(list(eps = eps), exp_p_ctrl, exp_q, exp_w))
  expect_identical(length(s13_pq), 13L)

  # 5. 18-field full set (prob + q + rates + controls + weights)
  s18 <- clip_summary(prob = prob, q = q, rates = rates, controls = ctrl, weights = weights)
  assert_compact_scalar_list(s18, c("eps", p_fields, q_fields, r_fields, w_fields))
  expect_identical(s18, c(list(eps = eps), exp_p_ctrl, exp_q, exp_r, exp_w))
  expect_identical(length(s18), 18L)

  # 6. Exhaustive check across all 16 input combinations (with controls supplied and NULL)
  for (use_p in c(FALSE, TRUE)) {
    for (use_q in c(FALSE, TRUE)) {
      for (use_r in c(FALSE, TRUE)) {
        for (use_w in c(FALSE, TRUE)) {
          for (use_ctrl in if (use_p) c(FALSE, TRUE) else FALSE) {
            got <- clip_summary(
              prob = if (use_p) prob else NULL,
              q = if (use_q) q else NULL,
              rates = if (use_r) rates else NULL,
              controls = if (use_ctrl) ctrl else NULL,
              weights = if (use_w) weights else NULL
            )
            exp_list <- c(
              list(eps = eps),
              if (use_p) (if (use_ctrl) exp_p_ctrl else exp_p_no_ctrl) else NULL,
              if (use_q) exp_q else NULL,
              if (use_r) exp_r else NULL,
              if (use_w) exp_w else NULL
            )
            assert_compact_scalar_list(got, names(exp_list))
            expect_identical(got, exp_list)
          }
        }
      }
    }
  }

  # 7. Interior unclipped regime and empty-vector inputs
  s_interior <- clip_summary(
    prob = c(0.2, 0.8),
    q = c(0.1, 0.9),
    rates = c(0.3, 0.6),
    controls = c(TRUE, FALSE),
    weights = c(2.0, 5.0)
  )
  assert_compact_scalar_list(s_interior, names(s18))
  expect_identical(
    s_interior,
    list(
      eps = eps,
      p_min = 0.2, p_max = 0.8, p_clip_max = 0.8, p_n = 2.0, p_capped = 0.0, p_capped_ctrl = 0.0,
      q_min = 0.1, q_max = 0.9, q_floor_min = 0.1, q_n = 2.0, q_floored = 0.0,
      r_min = 0.3, r_max = 0.6, r_floor_min = 0.3, r_n = 2.0, r_floored = 0.0,
      w_max = 5.0
    )
  )

  s_empty <- clip_summary(
    prob = numeric(0),
    q = numeric(0),
    rates = numeric(0),
    controls = logical(0),
    weights = numeric(0)
  )
  assert_compact_scalar_list(s_empty, names(s18))
  expect_identical(
    s_empty,
    list(
      eps = eps,
      p_min = NA_real_, p_max = NA_real_, p_clip_max = NA_real_, p_n = 0.0, p_capped = 0.0, p_capped_ctrl = 0.0,
      q_min = NA_real_, q_max = NA_real_, q_floor_min = NA_real_, q_n = 0.0, q_floored = 0.0,
      r_min = NA_real_, r_max = NA_real_, r_floor_min = NA_real_, r_n = 0.0, r_floored = 0.0,
      w_max = NA_real_
    )
  )
})

test_that("ipw_clip_summary_merge preserves exact extrema, counts, eps, and compact list structure across disjoint chunks", {
  clip_summary <- resolve_clip_helper("ipw_clip_summary")
  clip_merge <- resolve_clip_helper("ipw_clip_summary_merge")

  eps <- 1e-6
  cap <- 1 - 1e-6

  # Disjoint chunk 1 (no clipping binds here)
  p1 <- c(0.15, 0.65)
  c1 <- c(TRUE, FALSE)
  q1 <- c(0.25, 0.70)
  r1 <- c(0.30, 0.80)
  w1 <- c(1.25, 3.50)

  # Disjoint chunk 2 (upper probability cap and lower floor bind at and beyond eps)
  p2 <- c(0.05, cap, 0.9999998)
  c2 <- c(FALSE, TRUE, TRUE)
  q2 <- c(1e-7, eps, 0.95)
  r2 <- c(4e-7, eps, 0.50)
  w2 <- c(9.75, 2.00, 4.25)

  # Disjoint chunk 3 (interior values with new minimum p_min and new maximum w_max)
  p3 <- c(0.02, 0.40)
  c3 <- c(TRUE, TRUE)
  q3 <- c(0.10, 0.60)
  r3 <- c(0.15, 0.90)
  w3 <- c(15.5, 0.50)

  expected_18 <- list(
    eps = eps,
    p_min = 0.02,
    p_max = 0.9999998,
    p_clip_max = cap,
    p_n = 7.0,
    p_capped = 2.0,
    p_capped_ctrl = 2.0,
    q_min = 1e-7,
    q_max = 0.95,
    q_floor_min = eps,
    q_n = 7.0,
    q_floored = 2.0,
    r_min = 4e-7,
    r_max = 0.90,
    r_floor_min = eps,
    r_n = 7.0,
    r_floored = 2.0,
    w_max = 15.5
  )

  # Empty initial accumulator followed by disjoint chunks and trailing empty chunk
  acc_18 <- clip_summary(numeric(0), numeric(0), numeric(0), logical(0), numeric(0))
  ch1_18 <- clip_summary(p1, q1, r1, c1, w1)
  ch2_18 <- clip_summary(p2, q2, r2, c2, w2)
  ch3_18 <- clip_summary(p3, q3, r3, c3, w3)
  empty_18 <- clip_summary(numeric(0), numeric(0), numeric(0), logical(0), numeric(0))

  merged_18 <- clip_merge(clip_merge(clip_merge(clip_merge(acc_18, ch1_18), ch2_18), ch3_18), empty_18)
  assert_compact_scalar_list(merged_18, names(expected_18))
  expect_identical(merged_18, expected_18)
  expect_identical(
    clip_summary(c(p1, p2, p3), c(q1, q2, q3), c(r1, r2, r3), c(c1, c2, c3), c(w1, w2, w3)),
    expected_18
  )

  # 18-field merge with controls = NULL preserves p_capped_ctrl = NA_real_
  expected_18_no_ctrl <- modifyList(expected_18, list(p_capped_ctrl = NA_real_))
  merged_18_no_ctrl <- clip_merge(
    clip_summary(p1, q1, r1, controls = NULL, weights = w1),
    clip_merge(
      clip_summary(p2, q2, r2, controls = NULL, weights = w2),
      clip_summary(p3, q3, r3, controls = NULL, weights = w3)
    )
  )
  assert_compact_scalar_list(merged_18_no_ctrl, names(expected_18_no_ctrl))
  expect_identical(merged_18_no_ctrl, expected_18_no_ctrl)

  # Associativity and NULL handling
  expect_identical(clip_merge(NULL, ch1_18), ch1_18)
  expect_identical(clip_merge(ch1_18, NULL), ch1_18)
  expect_null(clip_merge(NULL, NULL))
  expect_identical(
    clip_merge(ch1_18, clip_merge(ch2_18, ch3_18)),
    clip_merge(clip_merge(ch1_18, ch2_18), ch3_18)
  )

  # 13-field (prob + rates + controls + weights) and (prob + q + controls + weights) chunk merges
  merged_13_pr <- clip_merge(
    clip_summary(p1, rates = r1, controls = c1, weights = w1),
    clip_merge(
      clip_summary(p2, rates = r2, controls = c2, weights = w2),
      clip_summary(p3, rates = r3, controls = c3, weights = w3)
    )
  )
  assert_compact_scalar_list(merged_13_pr, names(expected_18)[!grepl("^q_", names(expected_18))])
  expect_identical(merged_13_pr, expected_18[!grepl("^q_", names(expected_18))])

  merged_13_pq <- clip_merge(
    clip_summary(p1, q = q1, controls = c1, weights = w1),
    clip_merge(
      clip_summary(p2, q = q2, controls = c2, weights = w2),
      clip_summary(p3, q = q3, controls = c3, weights = w3)
    )
  )
  assert_compact_scalar_list(merged_13_pq, names(expected_18)[!grepl("^r_", names(expected_18))])
  expect_identical(merged_13_pq, expected_18[!grepl("^r_", names(expected_18))])

  # 7-field (q + weights) streaming pattern used by augmented and OR Corollary 1 consumers
  acc_7 <- clip_summary(q = numeric(0), weights = numeric(0))
  merged_7 <- clip_merge(
    clip_merge(
      clip_merge(acc_7, clip_summary(q = q1, weights = w1)),
      clip_summary(q = q2, weights = w2)
    ),
    clip_summary(q = q3, weights = w3)
  )
  expected_7 <- expected_18[c("eps", "q_min", "q_max", "q_floor_min", "q_n", "q_floored", "w_max")]
  assert_compact_scalar_list(merged_7, names(expected_7))
  expect_identical(merged_7, expected_7)

  # Merging two empty chunks preserves NA_real_ extrema and 0 counts
  merged_empty <- clip_merge(acc_18, empty_18)
  assert_compact_scalar_list(merged_empty, names(expected_18))
  expect_identical(merged_empty, acc_18)
})

test_that("compact clipping summary passes link and point-view storage guards while the numeric-vector mutant fails", {
  clip_summary <- resolve_clip_helper("ipw_clip_summary")
  verify_point_storage <- resolve_clip_helper("ipw_point_view_verify_compact_storage")
  verify_fitted_storage <- resolve_clip_helper("ipw_fitted_view_verify_compact_storage")

  s13 <- clip_summary(
    prob = c(0.2, 0.5, 1 - 1e-6),
    rates = c(0.4, 0.8, 1e-6),
    controls = c(FALSE, TRUE, TRUE),
    weights = c(1.0, 2.0, 3.0)
  )
  s18 <- clip_summary(
    prob = c(0.2, 0.5, 1 - 1e-6),
    q = c(0.3, 0.7, 1e-6),
    rates = c(0.4, 0.8, 1e-6),
    controls = c(FALSE, TRUE, TRUE),
    weights = c(1.0, 2.0, 3.0)
  )

  # Compact named list satisfies max(atomic_lengths) == 1 <= 9 and storage verifiers
  fit_compact_13 <- list(diagnostics = list(clipping = s13))
  fit_compact_18 <- list(diagnostics = list(clipping = s18))
  expect_identical(max(unlist(rapply(fit_compact_13, length, how = "unlist"))), 1L)
  expect_identical(max(unlist(rapply(fit_compact_18, length, how = "unlist"))), 1L)
  expect_lte(max(unlist(rapply(fit_compact_13, length, how = "unlist"))), 9L)
  expect_invisible(verify_point_storage(fit_compact_13, p_max = 1L, m_max = 2L))
  expect_invisible(verify_point_storage(fit_compact_18, p_max = 1L, m_max = 2L))
  expect_invisible(verify_fitted_storage(fit_compact_13))
  expect_invisible(verify_fitted_storage(fit_compact_18))

  # Original numeric-vector mutant (same names and scalar values as an atomic double vector) fails both guards
  mutant_vec_13 <- unlist(s13)
  mutant_vec_18 <- unlist(s18)
  expect_true(is.atomic(mutant_vec_13) && is.double(mutant_vec_13))
  expect_identical(length(mutant_vec_13), 13L)
  expect_identical(length(mutant_vec_18), 18L)

  fit_mutant_13 <- list(diagnostics = list(clipping = mutant_vec_13))
  fit_mutant_18 <- list(diagnostics = list(clipping = mutant_vec_18))
  expect_identical(max(unlist(rapply(fit_mutant_13, length, how = "unlist"))), 13L)
  expect_gt(max(unlist(rapply(fit_mutant_13, length, how = "unlist"))), 9L)
  expect_identical(max(unlist(rapply(fit_mutant_18, length, how = "unlist"))), 18L)

  expect_error(
    verify_point_storage(fit_mutant_13, p_max = 3L, m_max = 3L),
    regexp = "Non-compact vector length \\(13\\) retained at `point_link_view\\$diagnostics\\$clipping`\\.",
    class = "cdid_ipw_point_view_storage_error"
  )
  expect_error(
    verify_point_storage(fit_mutant_18, p_max = 3L, m_max = 3L),
    regexp = "Non-compact vector length \\(18\\) retained at `point_link_view\\$diagnostics\\$clipping`\\.",
    class = "cdid_ipw_point_view_storage_error"
  )
})

test_that("private clipping guard accepts inactive scalar lists and refuses active or malformed counts", {
  refuse_clipped <- resolve_clip_helper("ipw_point_view_refuse_clipped")
  inactive <- list(
    eps = 1e-6,
    p_min = 0.2, p_max = 0.8, p_clip_max = 0.8, p_n = 2.0,
    p_capped = 0.0, p_capped_ctrl = 0.0,
    r_min = 0.1, r_max = 0.6, r_floor_min = 0.1, r_n = 2.0,
    r_floored = 0.0, w_max = 3.0
  )
  view <- function(clip) list(diagnostics = list(clipping = clip))
  before <- inactive
  expect_invisible(refuse_clipped(view(inactive)))
  expect_identical(inactive, before)
  expect_invisible(refuse_clipped(view(unlist(inactive, use.names = TRUE))))
  integer_counts <- inactive
  integer_counts$r_floored <- integer_counts$p_capped_ctrl <- 0L
  expect_invisible(refuse_clipped(view(integer_counts)))
  # The guard covers consumed control odds; a treated-only cap is unchanged.
  treated_only <- inactive; treated_only$p_capped <- 1.0
  expect_invisible(refuse_clipped(view(treated_only)))

  for (field in c("r_floored", "p_capped_ctrl")) {
    active <- inactive; active[[field]] <- 1.0
    expect_error(refuse_clipped(view(active)), class = "cdid_ipw_point_view_clipping_error")
    expect_error(refuse_clipped(view(unlist(active, use.names = TRUE))),
      class = "cdid_ipw_point_view_clipping_error")
    missing_count <- inactive; missing_count[[field]] <- NULL
    expect_error(refuse_clipped(view(missing_count)), class = "cdid_ipw_point_view_clipping_error")
    expect_error(refuse_clipped(view(unlist(missing_count, use.names = TRUE))),
      class = "cdid_ipw_point_view_clipping_error")
    for (bad in list(NA_real_, NaN, Inf, -Inf, -1, 0.5, "0", FALSE,
                     numeric(0), c(0, 0), list(0), matrix(0, 1L, 1L),
                     structure(0, class = "malformed_count"), 0 + 0i)) {
      malformed <- inactive; malformed[[field]] <- bad
      expect_error(refuse_clipped(view(malformed)), class = "cdid_ipw_point_view_clipping_error")
    }
  }
  duplicate <- c(inactive, list(r_floored = 0.0))
  for (malformed in list(NULL, list(), unname(inactive), duplicate,
                        unlist(duplicate, use.names = TRUE),
                        structure(inactive, class = "malformed_record"),
                        matrix(0, 1L, 2L, dimnames = list(NULL, c("r_floored", "p_capped_ctrl"))))) {
    expect_error(refuse_clipped(view(malformed)), class = "cdid_ipw_point_view_clipping_error")
  }
})
