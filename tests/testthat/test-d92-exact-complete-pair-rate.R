# SPDX-License-Identifier: GPL-3.0-only
# Prepared for the independently reviewed installed candidate. Test bodies,
# fixture RNG and fits have NOT been run in this bounded adaptation workspace.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE))
  stop("Actual installed cdid required; no source fallback.", call. = FALSE)
d92_rate_api <- asNamespace("cdid")
stopifnot(isNamespace(d92_rate_api), environmentIsLocked(d92_rate_api))

# Reuse only these existing fixture declarations. Never source the other file
# or evaluate its library setup, test bodies, fit helpers or pipeline runners.
d92_rate_fixtures <- new.env(parent = baseenv())
d92_rate_fixtures$api <- d92_rate_api
for (expr in parse(testthat::test_path("test-d92-phase-retained-if.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in%
      c("d92_ctrl", "d92_glm", "d92_dated_panel", "d92_prepare_dated_augmented"))
    eval(expr, d92_rate_fixtures)
}
stopifnot(all(c("d92_ctrl", "d92_glm", "d92_dated_panel",
                "d92_prepare_dated_augmented") %in% ls(d92_rate_fixtures)))

d92_rate_owners <- function(d, B) {
  # The retained fixture's actual first_period reference weights, including
  # its original source construction; no replacement weights or fit needed.
  prep <- d92_rate_fixtures$d92_prepare_dated_augmented(d)
  ctrl <- d92_rate_fixtures$d92_ctrl(B)
  po <- d92_rate_api$construct_ipw_panel_owner(
    prep$data_internal, "id", "period", "y", "G", ctrl)
  fo <- d92_rate_api$construct_ipw_frame_owner(
    po, prep$frame_ids, prep$frame_Z, weights = prep$frame_w,
    control = ctrl, frame_g = prep$frame_g_arg,
    allow_frame_only = FALSE, design_on_A_only = TRUE)
  fctrl <- d92_rate_api$ipw_fitted_view_control(B, fo$n_a_units,
    ctrl$max_retained_entries, ctrl$max_matrix_entries, ctrl$max_work_ops,
    ctrl$max_output_entries)
  pctrl <- d92_rate_api$ipw_point_view_control(B, fo$n_a_units,
    ctrl$max_retained_entries, ctrl$max_matrix_entries, ctrl$max_work_ops,
    ctrl$max_output_entries)
  list(po = po, fo = fo, ctrl = ctrl, fctrl = fctrl, pctrl = pctrl, prep = prep)
}

# A single local language mutant removes precisely the complete-support fact;
# all original scans, sums and refusal guards remain the installed function's.
# No namespace binding is replaced. The exact-one assertions below must also
# fail if this fact is removed from the installed candidate itself.
d92_rate_without_complete <- utils::removeSource(
  d92_rate_api$ipw_owner_chunked_pair_support)
d92_rate_fact <- utils::removeSource(quote(if (isTRUE(require_support)) {
  rates[pair_count == coh_desc$group_a_count] <- 1
  ops <- ops + 4 * m
}))
d92_rate_fact_at <- which(vapply(as.list(body(d92_rate_without_complete)),
  identical, logical(1L), y = d92_rate_fact))
stopifnot(length(d92_rate_fact_at) == 1L)
body(d92_rate_without_complete)[[d92_rate_fact_at]] <- quote(NULL)

test_that("the retained seed43 fixture proves exact complete rates across partitions", {
  d <- d92_rate_fixtures$d92_dated_panel(N = 48L, seed = 43L)
  expect_true(all(is.finite(d$w) & d$w > 0))
  expect_true(all(d$w[d$period == 2L] > d$w[d$period == 1L]))
  expect_identical(sum(d$period == 1L), 48L)
  expect_identical(sum(d$period == 2L), 48L)
  for (B in c(1L, 2L, 7L, 8L, 9L, 11L, 15L, 16L, 17L, 47L, 48L, 49L)) {
    fx <- d92_rate_owners(d, B)
    nu <- d92_rate_api$construct_ipw_nuisance_view(
      fx$fo, 3, c(5, Inf), 1, 2, fx$ctrl)
    raw <- d92_rate_api$construct_ipw_nuisance_view(
      fx$fo, 3, c(5, Inf), 1, 2, fx$ctrl, require_support = FALSE)
    lk <- d92_rate_api$construct_ipw_link_view(
      fx$fo, 1, 2, g = 3, K = c(5, Inf), control = fx$ctrl)
    expect_identical(unname(nu$support$by_cohort_pair_count), rep(16L, 3L))
    expect_identical(nu$support$by_cohort_pair_count, nu$support$by_cohort_a_count)
    expect_identical(unname(nu$support$rates), rep(1, 3L))
    expect_identical(lk$support, nu$support)
    expect_identical(nu$support$by_cohort_pair_mass, raw$support$by_cohort_pair_mass)
    expect_identical(nu$support$q, raw$support$q)
    expect_identical(nu$work$actual_work_ops, raw$work$actual_work_ops + 12)
    expect_lte(nu$work$actual_work_ops, nu$preflight$work_ops_bound)
    expect_lte(lk$work$actual_work_ops, lk$preflight$work_ops_bound)

    # Scalar and the fitted/point workers' reconstructed observer buffers use
    # canonical A order, so complete cohorts already have exact scalar rates.
    sup <- list(m_cohorts = 3L, group_mass = unname(nu$support$by_cohort_a_mass))
    hit <- d92_rate_api$ipw_owner_query_pair_mask(fx$fo, 1, 2)
    risk <- d92_rate_api$ipw_pair_risk_set(fx$fo, nu$cohorts$cohorts, hit)
    scalar <- d92_rate_api$ipw_pair_rates(
      risk$risk_gi, fx$fo$w_a[risk$risk_a_idx], sup)
    chunks <- lapply(seq.int(1L, 48L, by = B), function(start)
      d92_rate_api$query_ipw_view_chunk(lk, start, min(48L, start + B - 1L),
        control = fx$ctrl, read_outcomes = FALSE, include_design = FALSE))
    gi <- unlist(lapply(chunks, function(ch) ch$observed_risk$cohort_index), use.names = FALSE)
    a <- unlist(lapply(chunks, function(ch) ch$observed_risk$a_index), use.names = FALSE)
    reconstructed <- d92_rate_api$ipw_pair_rates(gi, fx$fo$w_a[a], sup)
    expect_identical(reconstructed, scalar)
    expect_identical(scalar$rates, unname(nu$support$rates))

    if (B %in% c(9L, 11L)) {
      # Portable mathematical and raw-rounding witness at 1e-12: raw rates equal
      # chunk-accumulated pair_mass / a_mass, agree with exact 1 within 1e-12,
      # and exhibit at least one cohort rounding above 1 that the old mutant refuses.
      expect_identical(typeof(raw$support$rates), "double")
      expect_identical(names(raw$support$rates), names(nu$support$rates))
      expect_identical(raw$support$by_cohort_pair_count, nu$support$by_cohort_pair_count)
      expect_identical(raw$support$by_cohort_a_count, nu$support$by_cohort_a_count)
      expect_identical(raw$support$by_cohort_a_mass, nu$support$by_cohort_a_mass)
      expect_identical(raw$support$rates, raw$support$by_cohort_pair_mass / raw$support$by_cohort_a_mass)
      expect_true(all(is.finite(raw$support$rates)))
      expect_true(all(abs(raw$support$rates - nu$support$rates) <= 1e-12))
      expect_true(any(raw$support$rates > 1))
      expect_error(d92_rate_without_complete(fx$fo, nu$defining_pair, nu$cohorts, B),
        "Invalid marginal pair observation rates", class = "cdid_ipw_rate_error")

      # Exercise the actual downstream reconstruction check as well.
      bound <- d92_rate_api$construct_ipw_link_view(
        fx$fo, 1, 2, nuisance_view = nu, control = fx$ctrl)
      fitted <- d92_rate_api$fit_ipw_nuisance_view(nu,
        glm_control = d92_rate_fixtures$d92_glm, control = fx$fctrl)
      point <- d92_rate_api$fit_ipw_point_link_view(fitted, bound, control = fx$pctrl)
      expect_identical(unname(fitted$masses$rates), scalar$rates)
      expect_identical(fitted$masses$view_chunked_rates, nu$support$rates)
      expect_identical(fitted$masses$view_chunked_pair_weight, nu$support$by_cohort_pair_mass)
      expect_identical(point$masses$rates, fitted$masses$rates)
      expect_identical(point$point_summary_disclosure$link_chunked_rates, lk$support$rates)
      expect_identical(point$point_summary_disclosure$max_abs_rate_chunk_diff, 0)
      expect_identical(point$means$mean_a, fitted$means$mean_a)
      expect_identical(point$means$mean_b, fitted$means$mean_b)
    }
  }
})

test_that("mixed support changes only complete cohorts and leaves incomplete IF normalizers intact", {
  d <- d92_rate_fixtures$d92_dated_panel(N = 48L, seed = 43L)
  mixed <- d[!(d$id == 1L & d$period == 2L), , drop = FALSE]
  for (B in c(9L, 11L)) {
    fx <- d92_rate_owners(mixed, B)
    nu <- d92_rate_api$construct_ipw_nuisance_view(fx$fo, 3, c(5, Inf), 1, 2, fx$ctrl)
    raw <- d92_rate_api$construct_ipw_nuisance_view(
      fx$fo, 3, c(5, Inf), 1, 2, fx$ctrl, require_support = FALSE)
    lk <- d92_rate_api$construct_ipw_link_view(
      fx$fo, 1, 2, g = 3, K = c(5, Inf), control = fx$ctrl)
    expect_identical(unname(nu$support$by_cohort_pair_count), c(15L, 16L, 16L))
    expect_identical(nu$support$rates[1L], raw$support$rates[1L])
    expect_lt(nu$support$rates[1L], 1)
    expect_identical(unname(nu$support$rates[-1L]), c(1, 1))
    expect_identical(nu$support$by_cohort_pair_mass, raw$support$by_cohort_pair_mass)
    expect_identical(nu$support$q, raw$support$q)
    expect_identical(lk$support, nu$support)

    # The unmodified retained fixture's (1,4) pair is incomplete in ALL cohorts.
    # Compare to the local original-function mutant, then fit both unchanged
    # support descriptions and check the real IF normalizers and loadings.
    fx <- d92_rate_owners(d, B)
    nu <- d92_rate_api$construct_ipw_nuisance_view(fx$fo, 3, c(5, Inf), 1, 4, fx$ctrl)
    raw <- d92_rate_api$construct_ipw_nuisance_view(
      fx$fo, 3, c(5, Inf), 1, 4, fx$ctrl, require_support = FALSE)
    original <- d92_rate_without_complete(fx$fo, nu$defining_pair, nu$cohorts, B)
    expect_true(all(nu$support$by_cohort_pair_count < nu$support$by_cohort_a_count))
    expect_identical(nu$support, raw$support)
    expect_identical(nu$support$rates, original$rates)
    expect_identical(nu$support$by_cohort_pair_mass, original$by_cohort_pair_mass)
    fits <- lapply(list(nu, raw), function(view)
      d92_rate_api$fit_ipw_nuisance_view(view,
        glm_control = d92_rate_fixtures$d92_glm, control = fx$fctrl))
    expect_identical(fits[[1L]]$masses, fits[[2L]]$masses)
    expect_identical(fits[[1L]]$means, fits[[2L]]$means)
    points <- lapply(seq_along(fits), function(i) {
      lk <- d92_rate_api$construct_ipw_link_view(fx$fo, 1, 4,
        nuisance_view = list(nu, raw)[[i]], control = fx$ctrl)
      d92_rate_api$fit_ipw_point_link_view(fits[[i]], lk, control = fx$pctrl)
    })
    expect_identical(points[[1L]]$masses, points[[2L]]$masses)
    expect_identical(points[[1L]]$means, points[[2L]]$means)
    expect_identical(points[[1L]]$loadings, points[[2L]]$loadings)
  }
})

test_that("support opt-out and genuine input and owner refusals are preserved", {
  d <- d92_rate_fixtures$d92_dated_panel(N = 48L, seed = 43L)
  fx <- d92_rate_owners(d, 9L)
  # Absent cohort: zero count/mass, NA rate even though 0 == 0.
  expect_error(d92_rate_api$construct_ipw_nuisance_view(
    fx$fo, 3, c(5, Inf, 77), 1, 2, fx$ctrl), class = "cdid_ipw_cohort_support_error")
  absent <- d92_rate_api$construct_ipw_nuisance_view(
    fx$fo, 3, c(5, Inf, 77), 1, 2, fx$ctrl, require_support = FALSE)
  expect_identical(unname(absent$support$by_cohort_pair_count[4L]), 0L)
  expect_identical(unname(absent$support$by_cohort_pair_mass[4L]), 0)
  expect_true(is.na(absent$support$rates[4L]))
  expect_true(all(is.finite(absent$support$rates[1:3]) & abs(absent$support$rates[1:3] - 1) <= 1e-12))
  expect_true(any(absent$support$rates[1:3] > 1))
  absent_link <- d92_rate_api$construct_ipw_link_view(fx$fo, 1, 2,
    g = 3, K = c(5, Inf, 77), control = fx$ctrl, require_support = FALSE)
  expect_identical(absent_link$support, absent$support)
  for (constructor in c("construct_ipw_nuisance_view", "construct_ipw_link_view")) {
    make <- get(constructor, d92_rate_api)
    expect_error(make(fx$fo, g = 3, K = c(5, Inf), s = 1, t = 99,
      control = fx$ctrl), class = "cdid_ipw_pair_support_error")
    zero <- make(fx$fo, g = 3, K = c(5, Inf), s = 1, t = 99,
      control = fx$ctrl, require_support = FALSE)
    expect_identical(unname(zero$support$rates), rep(0, 3L))
    expect_error(make(fx$fo, g = 3, K = c(5, 5), s = 1, t = 2,
      control = fx$ctrl), class = "cdid_ipw_cohort_error")
    expect_error(make(fx$fo, g = 3, K = c(3, Inf), s = 1, t = 2,
      control = fx$ctrl), class = "cdid_ipw_cohort_error")
    expect_error(make(fx$fo, g = 3, K = c(5, Inf), s = 1, t = Inf,
      control = fx$ctrl), class = "cdid_ipw_pair_error")
    expect_error(make(fx$fo, g = 3, K = c(5, Inf), s = 1, t = 2,
      control = fx$ctrl, require_support = NA), class = "cdid_ipw_owner_input_error")
    expect_error(make(fx$fo, g = 3, K = c(5, Inf), s = 1, t = 2),
      class = "cdid_ipw_owner_input_error")
    altered <- fx$fo; altered$w_a[1L] <- altered$w_a[1L] * 2
    expect_error(make(altered, g = 3, K = c(5, Inf), s = 1, t = 2,
      control = fx$ctrl), class = "cdid_ipw_owner_integrity_error")
  }
  zero_d <- d[!(d$G == 5 & d$period == 2L), , drop = FALSE]
  zx <- d92_rate_owners(zero_d, 9L)
  expect_error(d92_rate_api$construct_ipw_nuisance_view(
    zx$fo, 3, c(5, Inf), 1, 2, zx$ctrl), class = "cdid_ipw_pair_support_error")
  zero <- d92_rate_api$construct_ipw_nuisance_view(
    zx$fo, 3, c(5, Inf), 1, 2, zx$ctrl, require_support = FALSE)
  expect_identical(unname(zero$support$rates[2L]), 0)
  expect_identical(unname(zero$support$by_cohort_pair_count[2L]), 0L)
  # Positive observed mass can still produce a genuinely invalid zero rate
  # by underflow in an incomplete cohort; the original strict guard must fire.
  ux <- d92_rate_owners(d[!(d$id == 1L & d$period == 2L), , drop = FALSE], 9L)
  uw <- rep(.Machine$double.xmin, 48L); uw[1L] <- 1e300
  uf <- d92_rate_api$construct_ipw_frame_owner(ux$po, ux$fo$a_id, ux$fo$Z_a,
    weights = uw, control = ux$ctrl)
  expect_error(d92_rate_api$construct_ipw_nuisance_view(
    uf, 3, c(5, Inf), 1, 2, ux$ctrl), class = "cdid_ipw_rate_error")
  underflow <- d92_rate_api$construct_ipw_nuisance_view(
    uf, 3, c(5, Inf), 1, 2, ux$ctrl, require_support = FALSE)
  expect_gt(underflow$support$by_cohort_pair_mass[1L], 0)
  expect_identical(unname(underflow$support$rates[1L]), 0)
  nu <- d92_rate_api$construct_ipw_nuisance_view(fx$fo, 3, c(5, Inf), 1, 2, fx$ctrl)
  expect_error(d92_rate_api$validate_ipw_nuisance_view(nu, frame_owner = zx$fo),
    class = "cdid_ipw_owner_mismatch_error")
  p <- fx$prep$data_internal
  expect_error(d92_rate_api$construct_ipw_panel_owner(
    rbind(p, p[1L, ]), "id", "period", "y", "G", fx$ctrl),
    class = "cdid_ipw_owner_duplicate_error")
  for (bad in c(Inf, -Inf, NaN)) {
    p_bad <- p; p_bad$y[1L] <- bad
    expect_error(d92_rate_api$construct_ipw_panel_owner(
      p_bad, "id", "period", "y", "G", fx$ctrl), class = "cdid_ipw_owner_outcome_error")
  }
  for (bad in c(0, -1, Inf, -Inf, NaN, NA_real_)) {
    w <- fx$fo$w_a; w[1L] <- bad
    expect_error(d92_rate_api$construct_ipw_frame_owner(fx$po, fx$fo$a_id,
      fx$fo$Z_a, weights = w, control = fx$ctrl), class = "cdid_ipw_owner_weight_error")
  }
  expect_error(d92_rate_api$construct_ipw_frame_owner(fx$po, fx$fo$a_id,
    fx$fo$Z_a, weights = rep(.Machine$double.xmax, 48L), control = fx$ctrl),
    class = "cdid_ipw_owner_weight_error")
  ids <- fx$fo$a_id; ids[2L] <- ids[1L]
  expect_error(d92_rate_api$construct_ipw_frame_owner(fx$po, ids, fx$fo$Z_a,
    weights = fx$fo$w_a, control = fx$ctrl), class = "cdid_ipw_owner_id_error")
})

test_that("existing finite preflight boundaries still admit exactly their charges", {
  ctrl <- d92_rate_fixtures$d92_ctrl(9L)
  for (name in c("preflight_ipw_nuisance_view", "preflight_ipw_link_view")) {
    preflight <- get(name, d92_rate_api)
    pre <- preflight(48L, 1L, 3L, ctrl)
    for (cap in c("max_retained_entries", "max_matrix_entries", "max_work_ops", "max_output_entries")) {
      charge <- switch(cap, max_retained_entries = pre$total_retained_entries,
        max_matrix_entries = pre$concurrent_entries,
        max_work_ops = pre$work_ops_bound, max_output_entries = pre$total_retained_entries)
      exact <- ctrl; exact[[cap]] <- charge
      expect_identical(preflight(48L, 1L, 3L, exact)$work_ops_bound, pre$work_ops_bound)
      below <- exact; below[[cap]] <- charge - 1
      expect_error(preflight(48L, 1L, 3L, below), class = "cdid_ipw_owner_bound_error")
    }
  }
})
