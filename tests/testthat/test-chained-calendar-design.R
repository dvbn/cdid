# Substantive tests for the private CP02 calendar/anticipation design path and
# six-role structural class validation (Proposition prop:calendar in
# docs/theory/cp02_moment_catalog.tex and Section 3 of
# docs/theory/cp02_moment_catalog_contract.md).

if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}

cal_reference_sha <- c(
  chained_design_qualified_b3db04.R = "b3db04e7255bf4b4126e4ffafcbf5c10138cbfcfc8d0085ab9d0083ebb8c55e9",
  cp02_calendar_design_reference.R = "ed66dddd0b84ebe4f0e4bb1ef64eb0866d92eab0f1be2b697762279dcdda3031"
)
cal_reference_sha[["chained_design_qualified_b3db04.R"]] <-
  "b3db04e7255bf4b4126e4ffafcbf5c10138cbfcfc8d0085ab9d00f3ebb8c55e9"

load_cal_reference <- function(file) {
  path <- testthat::test_path("reference", file)
  if (!file.exists(path)) stop("missing packaged reference: ", file, call. = FALSE)
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required to bind packaged references.", call. = FALSE)
  expect_identical(digest::digest(file = path, algo = "sha256"), cal_reference_sha[[file]])
  env <- new.env(parent = baseenv())
  if (identical(file, "chained_design_qualified_b3db04.R")) env$cdid_abort <- cdid:::cdid_abort
  sys.source(path, envir = env)
  env
}

cal_test_build <- function(...) {
  cdid:::materialize_chained_design(cdid:::build_calendar_chained_design(...), max_entries = 100000)
}

test_that("whole anticipation, explicit premises and intrinsic unlinked target bindings are enforced", {
  cal <- c(0, 3, 7, 9, 12, 14)
  coh <- c(6, 10, 14)
  expect_error(cdid:::validate_calendar_spec(cal, coh, 1.5), class = "cdid_design_anticipation_error")
  expect_error(cdid:::calendar_reference_date(6, cal, 1.5), class = "cdid_design_anticipation_error")
  expect_error(cdid:::eligible_control_cohorts(6, 3, coh, 1.5), class = "cdid_design_anticipation_error")
  # Preserve the original undeclared fixture as an active negative control.
  undeclared <- data.frame(g = 6, s = 3, t = 7, functional_key = "opaque", stringsAsFactors = FALSE)
  expect_error(cdid:::build_calendar_chained_design(undeclared, cal, c(coh, Inf), 2),
               class = "cdid_design_qualification_error")
  # TRUE qualifications are explicit conditional premises from the actual
  # baseline, not consequences of observed dates, counts or an opaque key.
  links <- data.frame(g = 6, s = c(3, 7), t = c(7, 9), functional_key = "last",
                      control_pool = "{14}", max_qualified_horizon = 9,
                      s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE)
  des <- cdid:::build_calendar_chained_design(links, cal, coh, 2, target_post_times = 14)
  expect_identical(des$post_targets$status, c("estimable", "estimable", "isolated_unsupported"))
  expect_invisible(cdid:::validate_ipw_design_incidence(des))
  # At u=14 the no-NT cutoff fails; the desired isolated coordinate still exists.
  expect_false(14 %in% cdid:::eligible_control_cohorts(6, 14, coh, 2))
  ntdes <- cdid:::build_calendar_chained_design(links, cal, c(coh, Inf), 2, target_post_times = 12)
  iso <- which(ntdes$vertices$time == 12)
  mut <- ntdes
  mut$vertices$role[iso] <- "post_nuisance"
  mut$vertices$functional_key[iso] <- "invented"
  mut$parameters$role[mut$parameters$vertex_index == iso] <- "post_nuisance"
  mut$parameters$functional_key[mut$parameters$vertex_index == iso] <- "invented"
  expect_error(cdid:::validate_ipw_design_incidence(mut), class = "cdid_ipw_class_vertex_error")
  mut <- ntdes
  mut$vertices$g[iso] <- 10
  mut$parameters$g[mut$parameters$vertex_index == iso] <- 10
  expect_error(cdid:::validate_ipw_design_incidence(mut), class = "cdid_ipw_class_vertex_error")
  mut <- ntdes
  mut$vertices$time[iso] <- 14
  mut$parameters$time[mut$parameters$vertex_index == iso] <- 14
  expect_error(cdid:::validate_ipw_design_incidence(mut), class = "cdid_ipw_class_vertex_error")
})

cal_dense_witness <- function(design, cert) {
  out <- numeric(design$incidence$n_links)
  if (!is.null(cert)) out[cert$link_index] <- cert$sign
  out
}

cal_verify_against_ref <- function(des, W_hand, ref_env) {
  expect_identical(unname(des$W), W_hand)
  rk_qr <- as.integer(qr(W_hand, tol = 1e-12)$rank)
  rk_svd <- as.integer(sum(svd(W_hand)$d > 1e-12))
  rk_rat <- ref_env$ref_rowspace_check(W_hand, numeric(ncol(W_hand)))$rank
  expect_identical(des$rank_diagnostics$rank, rk_qr)
  expect_identical(des$rank_diagnostics$rank, rk_svd)
  expect_identical(des$rank_diagnostics$rank, rk_rat)
  expect_identical(des$rank_diagnostics$nullity, ncol(W_hand) - rk_rat)
  expect_identical(des$rank_diagnostics$cycle_rank, nrow(W_hand) - rk_rat)

  B <- des$rank_diagnostics$nullspace_basis
  expect_identical(ncol(B), des$rank_diagnostics$nullity)
  if (ncol(B) > 0L) {
    expect_true(all(W_hand %*% B == 0L))
    expect_identical(qr(B)$rank, des$rank_diagnostics$nullity)
  }

  p <- ncol(W_hand)
  for (a in seq_len(p)) {
    ea <- numeric(p); ea[a] <- 1
    lev <- cdid:::check_level_estimable(des, a)
    rat_lev <- ref_env$ref_rowspace_check(W_hand, ea)$in_rowspace
    expect_identical(lev$estimable, rat_lev)
    expect_identical(des$parameters$level_estimable[a], rat_lev)
    if (lev$estimable) {
      w_vec <- cal_dense_witness(des, lev$witness_links)
      expect_equal(unname(drop(crossprod(w_vec, W_hand))), ea, tolerance = 0)
    }
    for (b in seq_len(p)) {
      cab <- numeric(p); cab[b] <- cab[b] + 1; cab[a] <- cab[a] - 1
      con <- cdid:::check_contrast_estimable(des, a, b)
      rat_con <- ref_env$ref_rowspace_check(W_hand, cab)$in_rowspace
      expect_identical(con$estimable, rat_con)
      if (con$estimable) {
        w_vec <- cal_dense_witness(des, con$witness_links)
        expect_equal(unname(drop(crossprod(w_vec, W_hand))), cab, tolerance = 0)
      }
    }
  }
}

test_that("frozen b3db04 reference is preserved bitwise on legacy calls and rejects unadapted anticipation roles", {
  b3 <- load_cal_reference("chained_design_qualified_b3db04.R")
  expect_identical(formals(cdid:::build_chained_design), formals(b3$build_chained_design))

  # 1. Unqualified legacy fixture:
  links_unqual <- data.frame(
    g = c(3, 3, 4, 4),
    s = c(1, 2, 1, 1),
    t = c(2, 3, 3, 4),
    ref = c(2, 2, 3, 3),
    functional_key = c("op:A", "op:A", "op:B1", "op:B2"),
    stringsAsFactors = FALSE
  )
  tgt <- data.frame(g = c(3, 3, 4), time = c(3, 4, 4))
  expect_identical(
    cdid:::build_chained_design(links_unqual, tgt),
    b3$build_chained_design(links_unqual, tgt)
  )
  expect_identical(
    cdid:::materialize_chained_design(cdid:::build_chained_design(links_unqual, tgt), 10000),
    b3$materialize_chained_design(b3$build_chained_design(links_unqual, tgt), 10000)
  )

  # 2. Qualified 4-role legacy fixture:
  links_qual4 <- data.frame(
    g = c(3, 3, 3),
    s = c(2, 3, 3),
    t = c(4, 4, 4),
    ref = c(2, 2, 2),
    functional_key = c("op:K4", "op:K4", "op:K5"),
    control_group = rep("notyettreated", 3L),
    control_pool = c("{4}", "{4}", "{5}"),
    max_qualified_horizon = c(3, 3, 4),
    s_qualified = c(TRUE, TRUE, TRUE),
    t_qualified = c(FALSE, FALSE, TRUE),
    s_coordinate_role = c("reference", "post_att", "post_att"),
    t_coordinate_role = c("post_nuisance", "post_nuisance", "post_att"),
    post_shareable = c(FALSE, FALSE, TRUE),
    uses_post_nuisance = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  expect_identical(
    cdid:::build_chained_design(links_qual4, c(3, 4)),
    b3$build_chained_design(links_qual4, c(3, 4))
  )

  # 3. Legacy build_chained_design without explicit calendar/cohorts/anticipation
  #    rejects anticipation roles and pre-treatment non-reference qualification:
  for (role in c("anticipation", "anticipation_att", "anticipation_nuisance")) {
    bad <- data.frame(
      g = 2009, s = 2002, t = 2005, ref = 2002, functional_key = "op:A",
      s_qualified = TRUE,
      t_qualified = role != "anticipation_nuisance",
      s_coordinate_role = "reference",
      t_coordinate_role = role,
      stringsAsFactors = FALSE
    )
    expect_error(cdid:::build_chained_design(bad), class = "cdid_design_qualification_error")
  }
})

test_that("irregular calendar {2000, 2002, 2005, 2009} with g = 2009: delta = 1 vs delta = 4 and baseline from date absent in selected links", {
  ref <- load_cal_reference("cp02_calendar_design_reference.R")
  fx <- ref$ref_irregular_2009_fixtures()
  cal <- fx$calendar
  cohorts <- c(2007, 2009, 2010, Inf)

  # 1. Exact latest pre-anticipation reference r_g = max{t in T : t < g - delta}:
  expect_identical(cdid:::calendar_reference_date(2009, cal, anticipation = 1), fx$delta1_ref)
  expect_identical(cdid:::calendar_reference_date(2009, cal, anticipation = 4), fx$delta4_ref)
  expect_identical(ref$ref_calendar_reference(2009, cal, 1), 2005)
  expect_identical(ref$ref_calendar_reference(2009, cal, 4), 2002)

  # Boundary equality delta = 4: 2009 - 4 = 2005 is NOT < 2005, so r_2009 = 2002 (not 2005),
  # and delta = 1 gives r_2009 = 2005 (not g - delta - 1 = 2007 and not ordinal shift 2002).
  links_d1 <- data.frame(
    g = c(2009, 2009, 2009),
    s = c(2000, 2002, 2005),
    t = c(2005, 2005, 2009),
    functional_key = rep("op:KInf", 3L),
    s_qualified = c(FALSE, FALSE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_d1 <- cal_test_build(links_d1, calendar = cal, cohorts = cohorts, anticipation = 1, target_post_times = 2009)
  expect_identical(des_d1$links$ref, rep(2005, 3L))
  expect_identical(des_d1$parameters$role, c("pre", "pre", "post"))
  expect_identical(
    des_d1$parameters$param_id,
    c("beta(g=2009,t=2000,fn=op:KInf)", "beta(g=2009,t=2002,fn=op:KInf)", "ATT(g=2009,t=2009)")
  )
  cal_verify_against_ref(des_d1, matrix(c(-1L, 0L, 0L, 0L, -1L, 0L, 0L, 0L, 1L), nrow = 3L, byrow = TRUE), ref)

  # Supplying g - delta - 1 = 2007 or ordinal index shift 2002 as ref under delta = 1 refuses:
  for (bad_ref in c(2007, 2002, 2000)) {
    bad <- links_d1; bad$ref <- bad_ref
    expect_error(
      cdid:::build_calendar_chained_design(bad, calendar = cal, cohorts = cohorts, anticipation = 1),
      class = "cdid_design_reference_error"
    )
  }

  # 2. Baseline from date absent in selected links:
  # Selected link is only (2002, 2009). Under delta = 1:
  # - On sub-calendar {2000, 2002, 2009}, r_2009 = 2002, so (2002, 2009) connects ref=2002 to ATT(2009,2009).
  # - On full calendar {2000, 2002, 2005, 2009}, r_2009 = 2005 (even though 2005 is absent from selected links!),
  #   so 2002 is an earlier pre-nuisance coordinate, ref(g=2009,t=2005) is isolated, and ATT(2009,2009) is NOT estimable!
  one_link <- data.frame(g = 2009, s = 2002, t = 2009, functional_key = "op:KInf", stringsAsFactors = FALSE)
  des_subcal <- cal_test_build(transform(one_link, s_qualified = TRUE, t_qualified = TRUE),
                              calendar = c(2000, 2002, 2009), cohorts = cohorts, anticipation = 1)
  expect_identical(des_subcal$links$ref, 2002)
  expect_identical(des_subcal$post_targets$level_estimable, TRUE)

  des_fullcal <- cal_test_build(transform(one_link, s_qualified = FALSE, t_qualified = TRUE),
                               calendar = cal, cohorts = cohorts, anticipation = 1)
  expect_identical(des_fullcal$links$ref, 2005)
  expect_identical(des_fullcal$vertices$time[des_fullcal$vertices$role == "reference"], 2005)
  expect_identical(des_fullcal$vertices$is_isolated[des_fullcal$vertices$role == "reference"], TRUE)
  expect_identical(des_fullcal$parameters$role, c("pre", "post"))
  expect_identical(des_fullcal$post_targets$level_estimable, FALSE)
  expect_identical(des_fullcal$post_targets$status, "disconnected_component")
  cal_verify_against_ref(des_fullcal, matrix(c(-1L, 1L), nrow = 1L), ref)

  # Claiming ref = 2002 (inferred from selected links) on the full calendar {2000, 2002, 2005, 2009} refuses:
  bad_inferred_ref <- one_link; bad_inferred_ref$ref <- 2002
  expect_error(
    cdid:::build_calendar_chained_design(bad_inferred_ref, calendar = cal, cohorts = cohorts, anticipation = 1),
    class = "cdid_design_reference_error"
  )

  # Empty pre-anticipation reference set (all calendar dates >= g - delta) refuses:
  expect_error(
    cdid:::build_calendar_chained_design(one_link, calendar = cal, cohorts = cohorts, anticipation = 9),
    class = "cdid_design_reference_error"
  )
  expect_error(
    cdid:::build_calendar_chained_design(one_link, calendar = cal, cohorts = cohorts, anticipation = 15),
    class = "cdid_design_reference_error"
  )
})

test_that("six endpoint roles, same actual date as private nuisance vs shared causal target, anticipation bridge, and post_targets restriction", {
  ref <- load_cal_reference("cp02_calendar_design_reference.R")
  fx <- ref$ref_irregular_2009_fixtures()
  cal <- fx$calendar
  cohorts <- c(2007, 2009, 2010, Inf)

  # Cohort g = 2009, delta = 4 => r_2009 = 2002, anticipation window [2005, 2009) contains 2005.
  # Control pools:
  # - K = {2007}: earlier-than-target cohort 2007 < 2009, clean at r_2009 = 2002 (2007 > 2002 + 4 = 2006),
  #   not clean at 2005 (2007 <= 2005 + 4 = 2009).
  # - K = {2010}: clean at 2002 and 2005 (2010 > 2009), not clean at 2009 (2010 <= 2013).
  # - K = {Inf}: clean at all finite dates.
  expect_identical(cdid:::eligible_control_cohorts(2009, 2002, cohorts, 4), c(2007, 2010, Inf))
  expect_identical(cdid:::eligible_control_cohorts(2009, 2005, cohorts, 4), c(2010, Inf))
  expect_identical(cdid:::eligible_control_cohorts(2009, 2009, cohorts, 4), Inf)
  expect_identical(ref$ref_eligible_controls(2009, 2002, cohorts, 4), c(2007, 2010, Inf))
  expect_identical(ref$ref_eligible_controls(2009, 2005, cohorts, 4), c(2010, Inf))
  expect_identical(ref$ref_eligible_controls(2009, 2009, cohorts, 4), Inf)

  links_6role <- data.frame(
    g = rep(2009, 5L),
    s = c(2000, 2002, 2002, 2005, 2005),
    t = c(2002, 2005, 2005, 2009, 2009),
    ref = rep(2002, 5L),
    functional_key = c("op:K2007", "op:K2007", "op:K2010", "op:K2010", "op:KInf"),
    control_group = c("notyettreated", "notyettreated", "notyettreated", "notyettreated", "nevertreated"),
    control_pool = c("{2007}", "{2007}", "{2010}", "{2010}", "{Inf}"),
    max_qualified_horizon = c(NA_real_, NA_real_, 2005, 2005, Inf),
    s_qualified = c(FALSE, TRUE, TRUE, TRUE, TRUE),
    t_qualified = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    s_coordinate_role = c("pre_nuisance", "reference", "reference", "anticipation_att", "anticipation_att"),
    t_coordinate_role = c("reference", "anticipation_nuisance", "anticipation_att", "post_nuisance", "post_att"),
    uses_post_nuisance = c(FALSE, FALSE, FALSE, TRUE, FALSE),
    uses_anticipation_nuisance = c(FALSE, TRUE, FALSE, FALSE, FALSE),
    post_shareable = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  des <- cal_test_build(links_6role, calendar = cal, cohorts = cohorts, anticipation = 4, target_post_times = 2009)
  expect_identical(
    des$vertices$role,
    c("reference", "pre", "anticipation_nuisance", "anticipation", "post_nuisance", "post")
  )
  expect_identical(des$parameters$role, fx$roles_six_role)
  expect_identical(
    des$parameters$param_id,
    c("beta(g=2009,t=2000,fn=op:K2007)",
      "beta_anticipation_nuisance(g=2009,t=2005,fn=op:K2007)",
      "ATT(g=2009,t=2005)",
      "beta_post_nuisance(g=2009,t=2009,fn=op:K2010)",
      "ATT(g=2009,t=2009)")
  )
  expect_identical(
    des$parameters$functional_key,
    c("op:K2007", "op:K2007", NA_character_, "op:K2010", NA_character_)
  )
  cal_verify_against_ref(des, fx$W_six_role, ref)

  # Post targets remain restricted to post-treatment (u >= g), never anticipation or nuisance:
  expect_identical(des$post_targets$param_id, "ATT(g=2009,t=2009)")
  expect_identical(des$post_targets$time, 2009)
  expect_identical(des$post_targets$level_estimable, TRUE)
  expect_identical(des$post_targets$status, "estimable")

  # Role counts in by_cohort and rank_diagnostics:
  bc <- des$rank_diagnostics$by_cohort
  expect_identical(bc$n_pre_parameters, 1L)
  expect_identical(bc$n_anticipation_nuisance_parameters, 1L)
  expect_identical(bc$n_anticipation_parameters, 1L)
  expect_identical(bc$n_estimable_anticipation, 1L)
  expect_identical(bc$n_post_nuisance_parameters, 1L)
  expect_identical(bc$n_post_parameters, 1L)
  expect_identical(bc$n_estimable_post, 1L)

  # Exact rational recovery:
  theta_true <- c(3 / 7, -5 / 7, 11 / 7, 13 / 7, 19 / 7)
  mu_true <- drop(fx$W_six_role %*% theta_true)
  expect_equal(unname(drop(qr.solve(des$W, mu_true))), theta_true, tolerance = 1e-14)
  cert_post <- cdid:::check_level_estimable(des, "ATT(g=2009,t=2009)")
  expect_identical(cert_post$witness_links, data.frame(link_index = c(3L, 5L), sign = c(1L, 1L)))
  expect_equal(sum(cal_dense_witness(des, cert_post$witness_links) * mu_true), 19 / 7, tolerance = 1e-15)

  # Contradictory qualification on repeated (g=2009, fn=op:K2010, date=2005) across rows 3 and 4 refuses:
  links_broken <- links_6role[, c("g", "s", "t", "ref", "functional_key", "s_qualified", "t_qualified")]
  links_broken$t_qualified[3L] <- FALSE
  expect_error(
    cal_test_build(links_broken, calendar = cal, cohorts = cohorts, anticipation = 4, target_post_times = 2009),
    class = "cdid_design_qualification_error"
  )

  # If date 2005 under op:K2010 is consistently unqualified (anticipation_nuisance), the cross-operator bridge breaks:
  links_broken$s_qualified[4L] <- FALSE
  des_broken <- cal_test_build(links_broken, calendar = cal, cohorts = cohorts, anticipation = 4, target_post_times = 2009)
  expect_identical(des_broken$post_targets$level_estimable, FALSE)
  expect_identical(des_broken$post_targets$time_graph_connected_to_ref, TRUE)
  expect_identical(des_broken$post_targets$status, "disconnected_component")
})

test_that("negative calendar dates, off-calendar g, boundary equality, no-NT latest-cohort boundary, and clean controls without PT domain", {
  ref <- load_cal_reference("cp02_calendar_design_reference.R")

  # 1. Negative calendar dates {-10, -6, -3, 0} and off-calendar target cohort g = -1 (with delta = 2):
  # g - delta = -3 => r_{-1} = -6; -10 is pre; -3 is boundary equality u = g - delta (anticipation); 0 is post.
  neg_cal <- c(-10, -6, -3, 0)
  neg_coh <- c(-1, 1, Inf)
  expect_false(-1 %in% neg_cal)
  expect_identical(cdid:::calendar_reference_date(-1, neg_cal, 2), -6)

  links_neg <- data.frame(
    g = rep(-1, 3L),
    s = c(-10, -6, -3),
    t = c(-6, -3, 0),
    functional_key = c("op:K1", "op:K1", "op:KInf"),
    control_pool = c("{1}", "{1}", "{Inf}"),
    max_qualified_horizon = c(-3, -3, Inf),
    s_qualified = c(FALSE, TRUE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_neg <- cal_test_build(links_neg, calendar = neg_cal, cohorts = neg_coh, anticipation = 2, target_post_times = 0)
  expect_identical(des_neg$parameters$role, c("pre", "anticipation", "post"))
  expect_identical(des_neg$parameters$time, c(-10, -3, 0))
  cal_verify_against_ref(des_neg, matrix(c(-1L, 0L, 0L, 0L, 1L, 0L, 0L, -1L, 1L), nrow = 3L, byrow = TRUE), ref)

  # Boundary equality k = u + delta: at u = -3 with delta = 2, k = -1 + 0 = -1 (target) or k = -3 + 2 = -1;
  # if cohorts include k = -4 (= ref + delta = -6 + 2), k = -4 is NOT > -4 so control_pool = "{-4}" refuses:
  bad_k_ref <- links_neg[1L, ]
  bad_k_ref$control_pool <- "{-4}"
  bad_k_ref$max_qualified_horizon <- NA_real_
  expect_error(
    cdid:::build_calendar_chained_design(bad_k_ref, calendar = neg_cal, cohorts = c(-4, -1, 1, Inf), anticipation = 2),
    class = "cdid_design_qualification_error"
  )

  # Boundary equality k = u + delta at anticipation endpoint u = -3 (u + delta = -1, or with u = 0, u + delta = 2):
  # For k = 2 at u = 0, delta = 2: k = u + delta = 2 is NOT > 2, so qualifying u = 0 under K = {2} refuses:
  bad_k_u <- data.frame(
    g = -1, s = -6, t = 0, functional_key = "op:K2",
    control_pool = "{2}", s_qualified = TRUE, t_qualified = TRUE,
    stringsAsFactors = FALSE
  )
  expect_error(
    cdid:::build_calendar_chained_design(bad_k_u, calendar = neg_cal, cohorts = c(-1, 2, Inf), anticipation = 2),
    class = "cdid_design_qualification_error"
  )

  # 2. No-never-treated latest-cohort boundary u + delta >= bar_g:
  # Cohorts = {2005, 2009} (no Inf, bar_g = 2009). Target g = 2005, delta = 0, cal = {2000, 2002, 2005, 2009}.
  # At u = 2005: 2005 + 0 < 2009 (can be qualified under K = {2009}).
  # At u = 2009: 2009 + 0 >= 2009 = bar_g (cannot be qualified on any link, but kept in post_targets!).
  cal_nont <- c(2000, 2002, 2005, 2009)
  coh_nont <- c(2005, 2009)
  links_nont <- data.frame(
    g = c(2005, 2005),
    s = c(2002, 2005),
    t = c(2005, 2009),
    functional_key = c("op:K2009", "op:K2009"),
    control_pool = c("{2009}", "{2009}"),
    max_qualified_horizon = c(2005, 2005),
    s_qualified = c(TRUE, TRUE),
    t_qualified = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  des_nont <- cal_test_build(links_nont, calendar = cal_nont, cohorts = coh_nont,
                             anticipation = 0, target_post_times = c(2005, 2009))
  expect_identical(
    des_nont$parameters$param_id,
    c("beta_post_nuisance(g=2005,t=2009,fn=op:K2009)", "ATT(g=2005,t=2005)", "ATT(g=2005,t=2009)")
  )
  expect_identical(des_nont$post_targets$param_id, c("ATT(g=2005,t=2005)", "ATT(g=2005,t=2009)"))
  expect_identical(des_nont$post_targets$level_estimable, c(TRUE, FALSE))
  expect_identical(des_nont$post_targets$is_isolated, c(FALSE, TRUE))
  expect_identical(des_nont$post_targets$status, c("estimable", "isolated_unsupported"))
  cal_verify_against_ref(des_nont, matrix(c(0L, 1L, 0L, 1L, -1L, 0L), nrow = 2L, byrow = TRUE), ref)

  # Falsely qualifying t = 2009 when bar_g = 2009 (even without control_pool column) refuses:
  bad_nont_qual <- data.frame(
    g = 2005, s = 2002, t = 2009, functional_key = "op:K2009",
    s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE
  )
  expect_error(
    cdid:::build_calendar_chained_design(bad_nont_qual, calendar = cal_nont, cohorts = coh_nont, anticipation = 0),
    class = "cdid_design_horizon_error"
  )
  expect_error(
    cdid:::build_calendar_chained_design(
      data.frame(g = 2005, s = 2002, t = 2009, functional_key = "op:K2009",
                 s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE),
      calendar = cal_nont, cohorts = coh_nont, anticipation = 0
    ),
    class = "cdid_design_horizon_error"
  )

  # Latest cohort itself g = 2009 = bar_g with earlier-than-target cohort 2007 clean at ref = 2005:
  links_latest_g <- data.frame(
    g = 2009, s = 2005, t = 2009, functional_key = "op:K2007",
    control_pool = "{2007}", s_qualified = TRUE, t_qualified = FALSE,
    stringsAsFactors = FALSE
  )
  des_latest_g <- cal_test_build(links_latest_g, calendar = cal_nont, cohorts = c(2007, 2009),
                                 anticipation = 0, target_post_times = 2009)
  expect_identical(des_latest_g$parameters$role, c("post_nuisance", "post"))
  expect_identical(des_latest_g$post_targets$param_id, "ATT(g=2009,t=2009)")
  expect_identical(des_latest_g$post_targets$status, "isolated_unsupported")
  cal_verify_against_ref(des_latest_g, matrix(c(1L, 0L), nrow = 1L), ref)

  # 3. Clean controls alone do NOT imply PT-domain membership:
  # (a) Supplying control_pool = "{Inf}" without qualification or PT-horizon metadata refuses:
  expect_error(
    cdid:::build_calendar_chained_design(
      data.frame(g = 2009, s = 2002, t = 2005, functional_key = "op:KInf",
                 control_pool = "{Inf}", stringsAsFactors = FALSE),
      calendar = cal_nont, cohorts = c(2009, Inf), anticipation = 4
    ),
    class = "cdid_design_qualification_error"
  )
  # (b) Clean controls K = {Inf} with explicit t_qualified = FALSE keeps t = 2005 as anticipation_nuisance
  #     and t = 2009 as post_nuisance when PT does not cover those dates:
  links_clean_unqual <- data.frame(
    g = c(2009, 2009), s = c(2002, 2002), t = c(2005, 2009),
    functional_key = c("op:KInf_noPT", "op:KInf_noPT"),
    control_pool = c("{Inf}", "{Inf}"),
    s_qualified = c(TRUE, TRUE), t_qualified = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  des_clean_unqual <- cal_test_build(links_clean_unqual, calendar = cal_nont,
                                     cohorts = c(2009, Inf), anticipation = 4, target_post_times = 2009)
  expect_identical(des_clean_unqual$parameters$role, c("anticipation_nuisance", "post_nuisance", "post"))
  expect_identical(des_clean_unqual$post_targets$level_estimable, FALSE)
})

test_that("exact error-free cutoff arithmetic prevents false certification under cancellation and refuses float overflow", {
  ref <- load_cal_reference("cp02_calendar_design_reference.R")

  # 1. Exact agreement between cdid:::exact_cmp_sum and hex-float bit-register ref_exact_cmp_sum
  #    on precision-edge cancellation/absorption triples:
  big <- 2^53 + 2 # 9007199254740994 (exact even significand at ulp = 2)
  big4 <- 2^53 + 4 # 9007199254740996 (exact even significand at ulp = 2)
  cases <- list(
    list(a = big, b = 1, c = big, expected = 1L),          # Tie rounds upward to big4; exact a+b > c.
    list(a = big, b = 0.75, c = big, expected = 1L),       # fl(a+b) == c AND fl(c-b) == a, yet a+b > c!
    list(a = -1, b = big4, c = big4, expected = -1L),      # fl(a+b) == c, yet a+b < c!
    list(a = -1.75, b = big4, c = big, expected = 1L),     # fl(a+b) == c, yet a+b = 2^53 + 2.25 > 2^53 + 2!
    list(a = -big4, b = big4, c = 1, expected = -1L),      # fl(c-b) == a, yet a+b = 0 < 1!
    list(a = -big4, b = big4, c = -1, expected = 1L),      # fl(c-b) == a, yet a+b = 0 > -1!
    list(a = 2^52 + 1, b = 2^52 + 1, c = big, expected = 0L),
    list(a = 2^52 + 1, b = 2^52 + 2, c = big, expected = 1L), # Tie rounds to big4; exact a+b = 2^53 + 3 > big.
    list(a = big, b = 1, c = Inf, expected = -1L)
  )
  for (tc in cases) {
    got <- cdid:::exact_cmp_sum(tc$a, tc$b, tc$c)
    exp_ref <- ref$ref_exact_cmp_sum(tc$a, tc$b, tc$c)
    expect_identical(got, tc$expected)
    expect_identical(got, exp_ref)
  }

  # 2. Design-level cancellation test:
  # Calendar = {-2, -1, 0, 2^53 + 4}, g = 2^53 + 4, delta = 2^53 + 4 => exact g - delta = 0.
  # Naive (-1) + (2^53 + 4) rounds to 2^53 + 4 == g, which would falsely classify t = -1 as >= g - delta!
  # Exact comparison finds r_g = -1, classifies -2 as "pre", -1 as "reference", 0 as "anticipation", and g as "post".
  cal_cancel <- c(-2, -1, 0, big4)
  coh_cancel <- c(big4, Inf)
  expect_identical(cdid:::calendar_reference_date(big4, cal_cancel, big4), -1)
  links_cancel <- data.frame(
    g = rep(big4, 3L),
    s = c(-2, -1, 0),
    t = c(-1, 0, big4),
    functional_key = rep("op:KInf", 3L),
    s_qualified = c(FALSE, TRUE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_cancel <- cal_test_build(links_cancel, calendar = cal_cancel, cohorts = coh_cancel, anticipation = big4)
  expect_identical(des_cancel$links$ref, rep(-1, 3L))
  expect_identical(des_cancel$parameters$role, c("pre", "anticipation", "post"))
  expect_identical(des_cancel$parameters$time, c(-2, 0, big4))

  # Exact clean-pool cutoff: big + 1 < big4, although the rounded sum is big4.
  expect_identical(cdid:::eligible_control_cohorts(big, big, c(big, big4, Inf), 1), c(big4, Inf))
  # Retain the original self-pool rejection control below: K contains target g.
  cal_cp <- c(2^53 - 2, 2^53, big, big4)
  bad_cp_cancel <- data.frame(
    g = big, s = 2^53 - 2, t = big, functional_key = "op:Kbig4",
    control_pool = paste0("{", sprintf("%.17g", big), "}"),
    s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE
  )
  expect_error(
    cdid:::build_calendar_chained_design(bad_cp_cancel, calendar = cal_cp, cohorts = c(big, big4, Inf), anticipation = 1),
    class = "cdid_design_qualification_error"
  )

  # 3. Overflowing cutoff arithmetic refuses explicitly rather than certifying +/-Inf:
  expect_error(cdid:::exact_cmp_sum(1e308, 1e308, 1e308), class = "cdid_design_arithmetic_error")
  expect_error(cdid:::exact_cmp_sum(1e308, 1e308, Inf), class = "cdid_design_arithmetic_error")
  expect_error(cdid:::exact_cmp_sum(0, 1e308, -1e308), class = "cdid_design_arithmetic_error")
  expect_error(
    cdid:::calendar_reference_date(-1e308, c(-1.5e308, -1e308, 0), anticipation = 1e308),
    class = "cdid_design_arithmetic_error"
  )
})

test_that("certify_ipw_structural_classes validates 6-role calendar designs intrinsically and refuses role/pool/incidence mutants", {
  ref <- load_cal_reference("cp02_calendar_design_reference.R")

  # Panel on dates 1..8, cohorts {5, 6, 9, Inf}.
  # Target g = 5, anticipation delta = 2 => r_5 = 2;
  # date 1 is pre (< 3); dates 3, 4 are anticipation ([3, 5)); dates 5, 6 are post (>= 5).
  # Control pool K1 = {6} is clean at r_5 = 2 (6 > 4) and at u = 3 (6 > 5),
  # but NOT clean at u = 4 (6 <= 6) or u = 5 (6 <= 7)!
  # Control pool K2 = {9, Inf} is clean at all dates 1..6 (9 > 6 + 2 = 8).
  n <- 96L
  G <- rep(c(5, 6, 9, Inf), each = 24L)
  x <- rep(seq(-1.2, 1.2, length.out = 24L), 4L)
  Z <- cbind(one = 1, x = x)
  w <- 0.8 + ((seq_len(n) * 11L) %% 23L) / 20
  obs <- matrix(FALSE, n, 8L)
  for (i in seq_len(n)) {
    k <- ((i - 1L) %% 24L) + 1L
    if (k <= 18L) obs[i, 1:6] <- TRUE else {
      obs[i, 7:8] <- TRUE
      obs[i, 1L + (k %% 6L)] <- TRUE
    }
  }
  Y <- outer(seq_len(n) / 19, 1:8, function(a, b) cos(a * b) + a + b / 7) + outer(x, 1:8, `*`)
  pan <- data.frame(
    id = rep(seq_len(n), each = 8L),
    time = rep(as.double(1:8), n),
    y = as.vector(t(ifelse(obs, Y, NA_real_))),
    g = rep(G, each = 8L)
  )
  idx <- cdid:::construct_observed_pairs(pan, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(idx, id = seq_len(n), Z = Z, weights = w)
  ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)
  bnds <- list(max_links = 50, max_classes = 10, max_class_dates = 10, max_E_entries = 200, max_mask_entries = 1e6)

  ses <- cdid:::construct_ipw_shared_session(ctx, max_nuisance = 4, max_links = 20, max_unit_entries = 1e5)
  nu_K1 <- cdid:::construct_ipw_nuisance(ses, 5, 6, 1, 2, ctrl)
  nu_K2 <- cdid:::construct_ipw_nuisance(ses, 5, c(9, Inf), 1, 2, ctrl)

  # Class 1 (op:K1, K = {6}, max_qualified_horizon = 3):
  #   1: (1, 2) [pre -> ref]
  #   2: (2, 3) [ref -> anticipation ATT(5,3)]
  #   3: (1, 3) [pre -> anticipation ATT(5,3)] (cycle edge in Class 1!)
  #   4: (2, 4) [ref -> anticipation_nuisance(5,4,op:K1)]
  #   5: (2, 5) [ref -> post_nuisance(5,5,op:K1)]
  # Class 2 (op:K2, K = {9, Inf}, max_qualified_horizon = 6):
  #   6: (2, 3) [ref -> anticipation ATT(5,3)]
  #   7: (3, 4) [anticipation ATT(5,3) -> anticipation ATT(5,4)]
  #   8: (2, 4) [ref -> anticipation ATT(5,4)] (cycle edge in Class 2!)
  #   9: (4, 5) [anticipation ATT(5,4) -> post ATT(5,5)]
  #  10: (5, 6) [post ATT(5,5) -> post ATT(5,6)]
  pairs_K1 <- list(c(1, 2), c(2, 3), c(1, 3), c(2, 4), c(2, 5))
  pairs_K2 <- list(c(2, 3), c(3, 4), c(2, 4), c(4, 5), c(5, 6))
  decl_K1 <- lapply(pairs_K1, function(p) {
    cdid:::declare_ipw_shared_link(ses, nu_K1, 5, 6, p[1L], p[2L], ctrl, functional_key = "op:K1")
  })
  decl_K2 <- lapply(pairs_K2, function(p) {
    cdid:::declare_ipw_shared_link(ses, nu_K2, 5, c(9, Inf), p[1L], p[2L], ctrl, functional_key = "op:K2")
  })
  decls <- c(decl_K1, decl_K2)

  links_df <- data.frame(
    g = rep(5, 10L),
    s = c(1, 2, 1, 2, 2, 2, 3, 2, 4, 5),
    t = c(2, 3, 3, 4, 5, 3, 4, 4, 5, 6),
    functional_key = c(rep("op:K1", 5L), rep("op:K2", 5L)),
    control_pool = c(rep("{6}", 5L), rep("{9,Inf}", 5L)),
    max_qualified_horizon = c(rep(3, 5L), rep(6, 5L)),
    s_qualified = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    t_qualified = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    s_coordinate_role = c("pre", "reference", "pre", "reference", "reference",
                          "reference", "anticipation", "reference", "anticipation", "post"),
    t_coordinate_role = c("reference", "anticipation", "anticipation", "anticipation_nuisance", "post_nuisance",
                          "anticipation", "anticipation", "anticipation", "post", "post"),
    stringsAsFactors = FALSE
  )
  des <- cdid:::build_calendar_chained_design(
    links_df, calendar = as.double(1:8), cohorts = c(5, 6, 9, Inf),
    anticipation = 2, target_post_times = c(5, 6)
  )
  expect_identical(
    sort(unique(des$vertices$role)),
    sort(c("reference", "pre", "anticipation_nuisance", "anticipation", "post_nuisance", "post"))
  )

  cert <- cdid:::certify_ipw_structural_classes(ses, decls, des, 1:10, "optimal", bnds)
  expect_s3_class(cert, "cdid_ipw_structural_certificate")
  expect_invisible(cdid:::verify_ipw_structural_certificate(cert, ses, decls, des, 1:10))
  expect_identical(cert$dropped_index, c(2L, 7L))
  # A different clean pool cannot be substituted for the actual nuisance pool.
  mut_pool <- des
  mut_pool$links$control_pool[1:5] <- "{9}"
  expect_error(cdid:::certify_ipw_structural_classes(ses, decls, mut_pool, 1:10, "optimal", bnds),
               class = "cdid_ipw_class_vertex_error")

  # Verify E * W_R == W and E * est_R == est exactly:
  Wfull <- cdid:::materialize_chained_design(des, 10000)$W
  expect_identical(
    cdid:::apply_ipw_reduction_map(cert, Wfull[cert$retained_index, , drop = FALSE]),
    unname(Wfull * 1.0)
  )
  all_fits <- lapply(decls, function(d) cdid:::fit_ipw_link_shared(ses, d))
  est <- vapply(all_fits, `[[`, numeric(1L), "estimate")
  expect_lt(max(abs(cdid:::apply_ipw_reduction_map(cert, est[cert$retained_index]) - est)), 1e-12)

  # Discriminating mutants on the 6-role calendar design:
  # 1. Stripping calendar_spec causes 4-role validator to reject anticipation roles:
  mut_nocal <- des; mut_nocal$calendar_spec <- NULL; attr(mut_nocal$links, "calendar_spec") <- NULL
  expect_error(
    cdid:::certify_ipw_structural_classes(ses, decls, mut_nocal, 1:10, "optimal", bnds),
    class = "cdid_ipw_design_error"
  )

  # 2. Mutating an anticipation vertex role to post or anticipation_nuisance refuses:
  v_ant3 <- which(des$vertices$role == "anticipation" & des$vertices$time == 3)
  for (bad_role in c("post", "anticipation_nuisance", "pre")) {
    mut_role <- des
    mut_role$vertices$role[v_ant3] <- bad_role
    mut_role$parameters$role[mut_role$parameters$vertex_index == v_ant3] <- bad_role
    expect_error(
      cdid:::certify_ipw_structural_classes(ses, decls, mut_role, 1:10, "optimal", bnds),
      class = "cdid_ipw_class_vertex_error"
    )
  }

  # 3. Mutating an anticipation_nuisance vertex role to anticipation refuses:
  v_anu4 <- which(des$vertices$role == "anticipation_nuisance" & des$vertices$time == 4)
  mut_anu <- des
  mut_anu$vertices$role[v_anu4] <- "anticipation"
  mut_anu$vertices$functional_key[v_anu4] <- NA_character_
  mut_anu$parameters$role[mut_anu$parameters$vertex_index == v_anu4] <- "anticipation"
  mut_anu$parameters$functional_key[mut_anu$parameters$vertex_index == v_anu4] <- NA_character_
  expect_error(
    cdid:::certify_ipw_structural_classes(ses, decls, mut_anu, 1:10, "optimal", bnds),
    class = "cdid_ipw_class_vertex_error"
  )

  # 4. Omitting control_pool from design$links to sneak a false qualification of t = 4 under K1 = {6}
  #    past build_calendar_chained_design is caught by certify_ipw_structural_classes against actual K = {6}:
  links_lie <- data.frame(
    g = 5, s = 2, t = 4, functional_key = "op:K1",
    s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE
  )
  des_lie <- cdid:::build_calendar_chained_design(links_lie, calendar = as.double(1:8),
                                                  cohorts = c(5, 6, 9, Inf), anticipation = 2)
  expect_error(
    cdid:::certify_ipw_structural_classes(ses, list(decl_K1[[4L]]), des_lie, 1L, "optimal", bnds),
    class = "cdid_ipw_class_vertex_error"
  )

  # 5. Mutating actual incidence start/end on the 6-role design refuses:
  mut_inc <- des
  mut_inc$incidence$end[4L] <- mut_inc$incidence$end[6L]
  expect_error(
    cdid:::certify_ipw_structural_classes(ses, decls, mut_inc, 1:10, "optimal", bnds),
    class = "cdid_ipw_class_vertex_error"
  )

  # 6. Isolated target_post_times vertex (t = 7) is accepted when valid, and refused when its role/date is mutated:
  des_iso <- cdid:::build_calendar_chained_design(
    links_df, calendar = as.double(1:8), cohorts = c(5, 6, 9, Inf),
    anticipation = 2, target_post_times = c(5, 6, 7)
  )
  expect_s3_class(
    cdid:::certify_ipw_structural_classes(ses, decls, des_iso, 1:10, "optimal", bnds),
    "cdid_ipw_structural_certificate"
  )
  v_iso7 <- which(des_iso$vertices$time == 7)
  mut_iso <- des_iso
  mut_iso$vertices$role[v_iso7] <- "anticipation"
  mut_iso$parameters$role[mut_iso$parameters$vertex_index == v_iso7] <- "anticipation"
  expect_error(
    cdid:::certify_ipw_structural_classes(ses, decls, mut_iso, 1:10, "optimal", bnds),
    class = "cdid_ipw_class_vertex_error"
  )
})

test_that("thin O(V+E) representation and explicit allocation bounds hold on 4,000-link 6-role calendar chains", {
  m <- 4000L
  cal <- as.double(0:m)
  coh <- c(2000, Inf)
  # g = 2000, anticipation = 500 => r_g = 1499;
  # dates 0..1498 are pre; 1499 is ref; 1500..1999 are anticipation; 2000..4000 are post.
  s_vec <- as.double(0:(m - 1L))
  t_vec <- as.double(1:m)
  q_pattern <- rep(c(TRUE, FALSE), length.out = m)
  links_thin <- data.frame(
    g = rep(2000, m),
    s = s_vec,
    t = t_vec,
    functional_key = rep("op:thin_cal", m),
    s_qualified = ifelse(s_vec == 1499, TRUE, ifelse(s_vec < 1500, FALSE, q_pattern)),
    t_qualified = ifelse(t_vec == 1499, TRUE, ifelse(t_vec < 1500, FALSE, c(q_pattern[-1L], TRUE))),
    stringsAsFactors = FALSE
  )
  des_thin <- cdid:::build_calendar_chained_design(links_thin, calendar = cal, cohorts = coh, anticipation = 500)
  expect_null(des_thin$W)
  expect_null(des_thin$rank_diagnostics$root_paths)
  expect_null(des_thin$rank_diagnostics$path_witnesses)
  expect_null(des_thin$rank_diagnostics$nullspace_basis)
  expect_lt(as.double(object.size(des_thin)), 8e6)
  expect_identical(des_thin$rank_diagnostics$rank, m)
  expect_true(all(c("reference", "pre", "anticipation_nuisance", "anticipation", "post_nuisance", "post") %in%
                  des_thin$vertices$role))
  expect_error(cdid:::materialize_chained_design(des_thin, 1000), class = "cdid_design_bound_error")
})
