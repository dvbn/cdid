# Tests for the private qualified-coordinate extension in R/chained_design.R.
# Uses independent hand-constructed incidence matrices and QR/SVD row-space
# references; dense materialization is explicitly bounded and test-only.

qual_test_build <- function(...) cdid:::materialize_chained_design(
  cdid:::build_chained_design(...), max_entries = 100000)

qual_dense_witness <- function(design, cert) {
  out <- numeric(design$incidence$n_links)
  if (!is.null(cert)) out[cert$link_index] <- cert$sign
  out
}

qual_qr_estimable <- function(W, c_vec, tol = 1e-12) {
  Wt <- t(as.matrix(W))
  qr_Wt <- qr(Wt, tol = tol)
  if (qr_Wt$rank == 0L) return(all(abs(c_vec) <= tol))
  max(abs(qr.resid(qr_Wt, as.numeric(c_vec)))) <= tol
}

qual_svd_estimable <- function(W, c_vec, tol = 1e-12) {
  sv <- svd(as.matrix(W), nu = 0L, nv = ncol(W))
  rk <- sum(sv$d > tol)
  if (rk == 0L) return(all(abs(c_vec) <= tol))
  V1 <- sv$v[, seq_len(rk), drop = FALSE]
  max(abs(c_vec - drop(V1 %*% crossprod(V1, c_vec)))) <= tol
}

qual_verify_all_functionals <- function(des, W_hand) {
  expect_identical(unname(des$W), W_hand)
  expect_identical(des$rank_diagnostics$rank, as.integer(qr(W_hand, tol = 1e-12)$rank))
  expect_identical(des$rank_diagnostics$rank, as.integer(sum(svd(W_hand)$d > 1e-12)))
  expect_identical(des$rank_diagnostics$nullity, ncol(W_hand) - des$rank_diagnostics$rank)
  expect_identical(des$rank_diagnostics$cycle_rank, nrow(W_hand) - des$rank_diagnostics$rank)

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
    expect_identical(lev$estimable, qual_qr_estimable(W_hand, ea))
    expect_identical(lev$estimable, qual_svd_estimable(W_hand, ea))
    expect_identical(des$parameters$level_estimable[a], lev$estimable)
    if (lev$estimable) {
      w_vec <- qual_dense_witness(des, lev$witness_links)
      expect_equal(unname(drop(crossprod(w_vec, W_hand))), ea, tolerance = 0)
    }
    for (b in seq_len(p)) {
      c_ab <- numeric(p); c_ab[b] <- c_ab[b] + 1; c_ab[a] <- c_ab[a] - 1
      con <- cdid:::check_contrast_estimable(des, a, b)
      expect_identical(con$estimable, qual_qr_estimable(W_hand, c_ab))
      expect_identical(con$estimable, qual_svd_estimable(W_hand, c_ab))
      if (con$estimable) {
        w_vec <- qual_dense_witness(des, con$witness_links)
        expect_equal(unname(drop(crossprod(w_vec, W_hand))), c_ab, tolerance = 0)
      }
    }
  }
}

test_that("N3 two rows [beta4, beta4 - ATT3] identify ATT3, keep ATT4 unestimable, and exclude post_nuisance from post_targets (D5)", {
  links_n3 <- data.frame(
    g = c(3, 3),
    s = c(2, 3),
    t = c(4, 4),
    ref = c(2, 2),
    functional_key = rep("op:g3:K={4}", 2L),
    control_group = rep("notyettreated", 2L),
    control_pool = rep("{4}", 2L),
    max_qualified_horizon = c(3, 3),
    s_qualified = c(TRUE, TRUE),
    t_qualified = c(FALSE, FALSE),
    s_coordinate_role = c("reference", "post_att"),
    t_coordinate_role = c("post_nuisance", "post_nuisance"),
    post_shareable = c(FALSE, FALSE),
    uses_post_nuisance = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  # 1. With target_post_times = 3: 2 parameters [beta4(K4), ATT(3,3)]
  des_2 <- qual_test_build(links_n3, target_post_times = 3)
  expect_identical(
    des_2$parameters$param_id,
    c("beta_post_nuisance(g=3,t=4,fn=op:g3:K={4})", "ATT(g=3,t=3)")
  )
  expect_identical(des_2$parameters$role, c("post_nuisance", "post"))
  W_hand_2 <- matrix(c(1L, 0L, 1L, -1L), nrow = 2L, byrow = TRUE)
  qual_verify_all_functionals(des_2, W_hand_2)
  expect_identical(des_2$post_targets$param_id, "ATT(g=3,t=3)")
  expect_identical(des_2$post_targets$status, "estimable")

  # 2. D5: With target_post_times = c(3, 4), requested ATT(3,4) is an isolated
  #    shared target while beta4(K4) is an anchored post_nuisance parameter.
  des_3 <- qual_test_build(links_n3, target_post_times = c(3, 4))
  expect_identical(
    des_3$parameters$param_id,
    c("beta_post_nuisance(g=3,t=4,fn=op:g3:K={4})", "ATT(g=3,t=3)", "ATT(g=3,t=4)")
  )
  expect_identical(des_3$parameters$role, c("post_nuisance", "post", "post"))
  expect_identical(
    des_3$parameters$functional_key,
    c("op:g3:K={4}", NA_character_, NA_character_)
  )
  W_hand_3 <- matrix(
    c(1L,  0L, 0L,
      1L, -1L, 0L),
    nrow = 2L, ncol = 3L, byrow = TRUE
  )
  qual_verify_all_functionals(des_3, W_hand_3)

  # Post targets contain ONLY shared ATT parameters, never post_nuisance:
  expect_identical(des_3$post_targets$param_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)"))
  expect_identical(des_3$post_targets$level_estimable, c(TRUE, FALSE))
  expect_identical(des_3$post_targets$is_isolated, c(FALSE, TRUE))
  expect_identical(des_3$post_targets$status, c("estimable", "isolated_unsupported"))
  # Naive time-only graph falsely connects t=4 to ref=2 via the unqualified link:
  expect_identical(des_3$post_targets$time_graph_connected_to_ref, c(TRUE, TRUE))

  # Role counts in by_cohort and rank_diagnostics do not count post_nuisance as pre or post ATT:
  bc <- des_3$rank_diagnostics$by_cohort
  expect_identical(bc$n_parameters, 3L)
  expect_identical(bc$n_pre_parameters, 0L)
  expect_identical(bc$n_post_nuisance_parameters, 1L)
  expect_identical(bc$n_post_parameters, 2L)
  expect_identical(bc$n_estimable_post, 1L)
  expect_identical(bc$n_isolated_post, 1L)
  expect_identical(bc$n_disconnected_post, 0L)
  expect_identical(des_3$rank_diagnostics$n_pre_parameters, 0L)
  expect_identical(des_3$rank_diagnostics$n_post_nuisance_parameters, 1L)
  expect_identical(des_3$rank_diagnostics$n_post_parameters, 2L)

  # Exact recovery on both 36-atom N3 fixtures (beta4 = 5 and beta4 = 1, ATT33 = 7, ATT34 = 11):
  cert_att3 <- cdid:::check_level_estimable(des_3, "ATT(g=3,t=3)")
  expect_true(cert_att3$estimable)
  expect_false(cdid:::check_level_estimable(des_3, "ATT(g=3,t=4)")$estimable)
  expect_true(cdid:::check_level_estimable(des_3, "beta_post_nuisance(g=3,t=4,fn=op:g3:K={4})")$estimable)
  w_att3 <- qual_dense_witness(des_3, cert_att3$witness_links)
  for (beta4_val in c(5, 1)) {
    mu_n3 <- c(beta4_val, beta4_val - 7)
    expect_equal(sum(w_att3 * mu_n3), 7, tolerance = 0)
  }

  # Direct converter handoff from tools/prototypes/chained_general_bridge_search.R:
  proto <- new.env(parent = globalenv())
  sys.source(testthat::test_path("reference", "chained_general_bridge_search.R"), envir = proto)
  ps_n3 <- data.frame(s = c(2, 3), t = c(4, 4))
  ps_n3$supported_cohorts <- list(c(3, 4), c(3, 4))
  ora_n3 <- proto$declare_population_overlap_oracle(c(3, 4, Inf), matrix(1, 3, 1), c(1, 1, 1))
  res_n3 <- proto$search_general_bridge_certificates(3, 1:4, c(3, 4, Inf), c(1, 1, 1),
                                                     "fixed_Q", ps_n3, c(3, 4), ora_n3, 100)
  handoff_links <- proto$general_bridge_certificates_to_design_links(res_n3)
  des_handoff <- qual_test_build(handoff_links, target_post_times = c(3, 4))
  expect_identical(unname(des_handoff$W), W_hand_3)
  expect_identical(des_handoff$post_targets$status, c("estimable", "isolated_unsupported"))
})

test_that("same post date as nuisance under K4 and qualified ATT under K5 remains distinct", {
  # Date 4 is unqualified post_nuisance under K4={4} (max_qualified_horizon = 3)
  # and qualified post_att under K5={5} (max_qualified_horizon = 4).
  # Pair (3,4) is present under BOTH K4 and K5 with distinct endpoint qualifications at t=4.
  links <- data.frame(
    g = rep(3, 3L),
    s = c(2, 3, 3),
    t = c(4, 4, 4),
    ref = rep(2, 3L),
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

  des <- qual_test_build(links, target_post_times = c(3, 4))
  expect_identical(
    des$parameters$param_id,
    c("beta_post_nuisance(g=3,t=4,fn=op:K4)", "ATT(g=3,t=3)", "ATT(g=3,t=4)")
  )
  expect_identical(des$parameters$role, c("post_nuisance", "post", "post"))

  # Columns: [beta4(K4), ATT(3,3), ATT(3,4)]
  W_hand <- matrix(
    c(1L,  0L, 0L,
      1L, -1L, 0L,
      0L, -1L, 1L),
    nrow = 3L, ncol = 3L, byrow = TRUE
  )
  qual_verify_all_functionals(des, W_hand)
  expect_identical(des$rank_diagnostics$rank, 3L)
  expect_identical(des$rank_diagnostics$nullity, 0L)
  expect_identical(des$post_targets$param_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)"))
  expect_identical(des$post_targets$status, c("estimable", "estimable"))

  # Population moments: beta4(K4)=5, ATT(3,3)=7, ATT(3,4)=11 => mu = (5, -2, 4).
  mu <- c(5, -2, 4)
  theta_hat <- drop(qr.solve(des$W, mu))
  expect_equal(unname(theta_hat), c(5, 7, 11), tolerance = 1e-12)
})

test_that("mixed operator composition joins only qualified vertices and keeps pre/post nuisances private", {
  # 1. Connected 3-segment composition across K4={4}, K5={5}, KNT={Inf},
  #    plus extra nuisance visits at early t=1 under K4 and post t=5 under K4
  #    to verify that neither early pre nor unqualified post vertices merge across operators.
  links_comp <- data.frame(
    g = rep(3, 8L),
    s = c(1, 2, 3, 3, 3, 4, 1, 1),
    t = c(2, 4, 4, 5, 5, 5, 4, 5),
    ref = rep(2, 8L),
    functional_key = c("op:K4", "op:K4", "op:K4", "op:K4", "op:K5", "op:K5", "op:KNT", "op:KNT"),
    control_pool = c("{4}", "{4}", "{4}", "{4}", "{5}", "{5}", "{Inf}", "{Inf}"),
    max_qualified_horizon = c(3, 3, 3, 3, 4, 4, Inf, Inf),
    s_qualified = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    t_qualified = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    post_shareable = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  des <- qual_test_build(links_comp, target_post_times = c(3, 4, 5))
  # Parameter order:
  # pre: (t=1, K4), (t=1, KNT)
  # post_nuisance: (t=4, K4), (t=5, K4), (t=5, K5)
  # post: ATT(3,3), ATT(3,4), ATT(3,5)
  expect_identical(
    des$parameters$param_id,
    c("beta(g=3,t=1,fn=op:K4)",
      "beta(g=3,t=1,fn=op:KNT)",
      "beta_post_nuisance(g=3,t=4,fn=op:K4)",
      "beta_post_nuisance(g=3,t=5,fn=op:K4)",
      "beta_post_nuisance(g=3,t=5,fn=op:K5)",
      "ATT(g=3,t=3)",
      "ATT(g=3,t=4)",
      "ATT(g=3,t=5)")
  )
  W_hand <- rbind(
    c(-1L,  0L, 0L, 0L, 0L,  0L,  0L, 0L), # (1,2) K4: 0 - beta1(K4)
    c( 0L,  0L, 1L, 0L, 0L,  0L,  0L, 0L), # (2,4) K4: beta4(K4) - 0
    c( 0L,  0L, 1L, 0L, 0L, -1L,  0L, 0L), # (3,4) K4: beta4(K4) - ATT3
    c( 0L,  0L, 0L, 1L, 0L, -1L,  0L, 0L), # (3,5) K4: beta5(K4) - ATT3
    c( 0L,  0L, 0L, 0L, 1L, -1L,  0L, 0L), # (3,5) K5: beta5(K5) - ATT3
    c( 0L,  0L, 0L, 0L, 1L,  0L, -1L, 0L), # (4,5) K5: beta5(K5) - ATT4
    c( 0L, -1L, 0L, 0L, 0L,  0L,  1L, 0L), # (1,4) KNT: ATT4 - beta1(KNT)
    c( 0L, -1L, 0L, 0L, 0L,  0L,  0L, 1L)  # (1,5) KNT: ATT5 - beta1(KNT)
  )
  qual_verify_all_functionals(des, W_hand)
  expect_identical(des$rank_diagnostics$rank, 8L)
  expect_identical(des$rank_diagnostics$n_pre_parameters, 2L)
  expect_identical(des$rank_diagnostics$n_post_nuisance_parameters, 3L)
  expect_identical(des$rank_diagnostics$n_post_parameters, 3L)
  expect_true(all(des$post_targets$level_estimable))

  # 2. Disconnected when two distinct operators K4a and K4b meet ONLY at unqualified t=4:
  links_disc <- data.frame(
    g = c(3, 3),
    s = c(2, 3),
    t = c(4, 4),
    ref = c(2, 2),
    functional_key = c("op:K4a", "op:K4b"),
    s_qualified = c(TRUE, TRUE),
    t_qualified = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  des_disc <- qual_test_build(links_disc, target_post_times = c(3, 4))
  W_disc_hand <- rbind(
    c(1L, 0L,  0L, 0L), # (2,4) K4a: beta4(K4a) - 0
    c(0L, 1L, -1L, 0L)  # (3,4) K4b: beta4(K4b) - ATT3
  )
  qual_verify_all_functionals(des_disc, W_disc_hand)
  expect_identical(des_disc$post_targets$level_estimable, c(FALSE, FALSE))
  expect_identical(des_disc$post_targets$time_graph_connected_to_ref, c(TRUE, TRUE))
  expect_identical(des_disc$post_targets$status, c("disconnected_component", "isolated_unsupported"))
})

test_that("lone post pair45 under KNT identifies ATT5 - ATT4 contrast while neither level is estimable (D5)", {
  links_45 <- data.frame(
    g = 3,
    s = 4,
    t = 5,
    ref = 2,
    functional_key = "op:KNT",
    control_group = "nevertreated",
    control_pool = "{Inf}",
    max_qualified_horizon = Inf,
    s_qualified = TRUE,
    t_qualified = TRUE,
    s_coordinate_role = "post_att",
    t_coordinate_role = "post_att",
    post_shareable = TRUE,
    uses_post_nuisance = FALSE,
    stringsAsFactors = FALSE
  )
  des <- qual_test_build(links_45, target_post_times = c(3, 4, 5))
  W_hand <- matrix(c(0L, -1L, 1L), nrow = 1L, ncol = 3L)
  qual_verify_all_functionals(des, W_hand)

  expect_identical(des$parameters$param_id, c("ATT(g=3,t=3)", "ATT(g=3,t=4)", "ATT(g=3,t=5)"))
  expect_identical(des$parameters$level_estimable, c(FALSE, FALSE, FALSE))
  expect_identical(
    des$post_targets$status,
    c("isolated_unsupported", "disconnected_component", "disconnected_component")
  )
  c45 <- cdid:::check_contrast_estimable(des, "ATT(g=3,t=4)", "ATT(g=3,t=5)")
  expect_true(c45$estimable)
  expect_true(c45$same_component)
  expect_false(c45$both_level_estimable)
  expect_identical(c45$witness_links, data.frame(link_index = 1L, sign = 1L))
  expect_false(cdid:::check_contrast_estimable(des, "ATT(g=3,t=3)", "ATT(g=3,t=4)")$estimable)
})

test_that("duplicate/contradictory qualifications, roles, horizons, pools, and early nonreference sharing are rejected", {
  base_n3 <- data.frame(
    g = c(3, 3),
    s = c(2, 3),
    t = c(4, 4),
    ref = c(2, 2),
    functional_key = rep("op:K4", 2L),
    control_group = rep("notyettreated", 2L),
    control_pool = rep("{4}", 2L),
    max_qualified_horizon = c(3, 3),
    s_qualified = c(TRUE, TRUE),
    t_qualified = c(FALSE, FALSE),
    s_coordinate_role = c("reference", "post_att"),
    t_coordinate_role = c("post_nuisance", "post_nuisance"),
    post_shareable = c(FALSE, FALSE),
    uses_post_nuisance = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  # 1. Duplicate link with identical qualification -> cdid_design_duplicate_link_error
  dup_same <- rbind(base_n3, base_n3[1L, ])
  expect_error(cdid:::build_chained_design(dup_same), class = "cdid_design_duplicate_link_error")

  # 2. Duplicate (g,s,t,fn) with contradictory qualification -> cdid_design_qualification_error
  dup_contra <- rbind(base_n3[, c("g", "s", "t", "ref", "functional_key", "s_qualified", "t_qualified")],
                      data.frame(g = 3, s = 2, t = 4, ref = 2, functional_key = "op:K4",
                                 s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE))
  expect_error(cdid:::build_chained_design(dup_contra), class = "cdid_design_qualification_error")

  # 3. Contradictory qualification across distinct pairs (2,4) and (3,4) sharing (g=3, fn=K4, date=4)
  bad_ep <- base_n3[, c("g", "s", "t", "ref", "functional_key", "s_qualified", "t_qualified")]
  bad_ep$t_qualified <- c(FALSE, TRUE)
  expect_error(cdid:::build_chained_design(bad_ep), class = "cdid_design_qualification_error")

  # 4. Contradictory qualification when date=4 is t on (2,4) and s on (4,5) under same (g=3, fn=K4)
  bad_st_ep <- data.frame(
    g = c(3, 3), s = c(2, 4), t = c(4, 5), ref = c(2, 2),
    functional_key = c("op:K4", "op:K4"),
    s_qualified = c(TRUE, TRUE), t_qualified = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  expect_error(cdid:::build_chained_design(bad_st_ep), class = "cdid_design_qualification_error")

  # 5. Early nonreference sharing refused under PT-post (both at s < g and at t < g)
  bad_early_s <- data.frame(
    g = 4, s = 1, t = 4, ref = 3, functional_key = "op:KNT",
    s_qualified = TRUE, t_qualified = TRUE, stringsAsFactors = FALSE
  )
  expect_error(cdid:::build_chained_design(bad_early_s), class = "cdid_design_qualification_error")
  expect_error(cdid:::build_chained_design(bad_early_s), class = "cdid_design_post_share_error")

  bad_early_t <- data.frame(
    g = 4, s = 1, t = 2, ref = 3, functional_key = "op:KNT",
    s_qualified = FALSE, t_qualified = TRUE, stringsAsFactors = FALSE
  )
  expect_error(cdid:::build_chained_design(bad_early_t), class = "cdid_design_qualification_error")

  # 6. Reference endpoint declared unqualified -> refused
  bad_ref_unqual <- data.frame(
    g = 3, s = 2, t = 3, ref = 2, functional_key = "op:K4",
    s_qualified = FALSE, t_qualified = TRUE, stringsAsFactors = FALSE
  )
  expect_error(cdid:::build_chained_design(bad_ref_unqual), class = "cdid_design_qualification_error")
  expect_error(cdid:::build_chained_design(bad_ref_unqual), class = "cdid_design_reference_error")

  # 7. Missing one of s_qualified/t_qualified, NA values, or non-logical/matrix columns
  only_s <- base_n3[, c("g", "s", "t", "ref", "functional_key", "s_qualified")]
  expect_error(cdid:::build_chained_design(only_s), class = "cdid_design_qualification_error")
  only_t <- base_n3[, c("g", "s", "t", "ref", "functional_key", "t_qualified")]
  expect_error(cdid:::build_chained_design(only_t), class = "cdid_design_qualification_error")
  na_qual <- base_n3[, c("g", "s", "t", "ref", "functional_key", "s_qualified", "t_qualified")]
  na_qual$t_qualified[1L] <- NA
  expect_error(cdid:::build_chained_design(na_qual), class = "cdid_design_qualification_error")
  mat_qual <- base_n3[, c("g", "s", "t", "ref", "functional_key", "s_qualified", "t_qualified")]
  mat_qual$s_qualified <- matrix(mat_qual$s_qualified, ncol = 1L)
  expect_error(cdid:::build_chained_design(mat_qual), class = "cdid_design_qualification_error")

  # 8. Contradictory coordinate roles, max_qualified_horizon, control_pool, post_shareable, uses_post_nuisance
  bad_role <- base_n3; bad_role$t_coordinate_role <- c("post_att", "post_att")
  expect_error(cdid:::build_chained_design(bad_role), class = "cdid_design_qualification_error")

  bad_mh1 <- base_n3; bad_mh1$max_qualified_horizon <- c(4, 4)
  expect_error(cdid:::build_chained_design(bad_mh1), class = "cdid_design_horizon_error")
  bad_mh2 <- base_n3; bad_mh2$max_qualified_horizon <- c(2, 2)
  expect_error(cdid:::build_chained_design(bad_mh2), class = "cdid_design_horizon_error")
  bad_mh3 <- base_n3; bad_mh3$max_qualified_horizon <- c(3, 4)
  expect_error(cdid:::build_chained_design(bad_mh3), class = "cdid_design_horizon_error")

  bad_pool1 <- base_n3; bad_pool1$control_pool <- c("{4}", "{4,Inf}")
  expect_error(cdid:::build_chained_design(bad_pool1), class = "cdid_design_qualification_error")
  bad_pool2 <- base_n3; bad_pool2$control_pool <- c("{3}", "{3}")
  expect_error(cdid:::build_chained_design(bad_pool2), class = "cdid_design_qualification_error")
  bad_pool3 <- base_n3; bad_pool3$control_pool <- c("{5}", "{5}")
  expect_error(cdid:::build_chained_design(bad_pool3), class = "cdid_design_qualification_error")

  bad_ps <- base_n3; bad_ps$post_shareable <- c(TRUE, FALSE)
  expect_error(cdid:::build_chained_design(bad_ps), class = "cdid_design_post_share_error")

  bad_upn <- base_n3; bad_upn$uses_post_nuisance <- c(FALSE, TRUE)
  expect_error(cdid:::build_chained_design(bad_upn), class = "cdid_design_qualification_error")

  bad_bare <- base_n3[, c("g", "s", "t", "ref", "functional_key", "s_qualified", "t_qualified")]
  bad_bare$functional_key <- rep("nyt:h=3", 2L)
  expect_error(cdid:::build_chained_design(bad_bare), class = "cdid_design_profile_merge_error")
})

test_that("exact large numeric dates/cohorts, legacy invariance, thin/bounded materialization, and forced post-sharing perturbation hold", {
  # 1. Exact large numeric dates/cohorts near 2^53 - 1 and 2^-40 separation:
  shift <- 9007199254740900
  cal <- shift + c(1, 2, 3, 4)
  links_large <- data.frame(
    g = rep(cal[3], 2L),
    s = c(cal[2], cal[3]),
    t = c(cal[4], cal[4]),
    ref = rep(cal[2], 2L),
    functional_key = rep("op:large:K4", 2L),
    control_pool = rep(paste0("{", sprintf("%.17g", cal[4]), "}"), 2L),
    max_qualified_horizon = rep(cal[3], 2L),
    s_qualified = c(TRUE, TRUE),
    t_qualified = c(FALSE, FALSE),
    post_shareable = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  des_large <- qual_test_build(links_large, target_post_times = c(cal[3], cal[4]))
  qual_verify_all_functionals(des_large, matrix(c(1L, 0L, 0L, 1L, -1L, 0L), nrow = 2L, byrow = TRUE))
  expect_identical(des_large$parameters$time, c(cal[4], cal[3], cal[4]))
  expect_identical(des_large$post_targets$time, c(cal[3], cal[4]))
  expect_identical(des_large$post_targets$level_estimable, c(TRUE, FALSE))

  # Sub-display-precision date separation (4 vs 4 + 2^-40) in qualified post_nuisance and post:
  d1 <- 4; d2 <- 4 + 2^-40
  links_eps <- data.frame(
    g = c(3, 3, 3),
    s = c(2, 3, 2),
    t = c(d1, d1, d2),
    ref = c(2, 2, 2),
    functional_key = c("op:K4", "op:K4", "op:K5"),
    s_qualified = c(TRUE, TRUE, TRUE),
    t_qualified = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  des_eps <- qual_test_build(links_eps)
  expect_identical(des_eps$parameters$time, c(d1, 3, d2))
  expect_identical(des_eps$parameters$role, c("post_nuisance", "post", "post"))

  # 2. Legacy schema and value invariance when s_qualified/t_qualified are absent:
  legacy_links <- data.frame(
    g = c(3, 3, 4),
    s = c(1, 2, 3),
    t = c(2, 3, 4),
    ref = c(2, 2, 3),
    functional_key = c("op:A", "op:A", "op:B"),
    horizon = c(3, 3, 4),
    post_shareable = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  des_leg <- cdid:::build_chained_design(legacy_links, target_post_times = data.frame(g = c(3, 3, 4), time = c(3, 4, 4)))
  expect_false("n_post_nuisance_parameters" %in% names(des_leg$rank_diagnostics$by_cohort))
  expect_false("n_post_nuisance_parameters" %in% names(des_leg$rank_diagnostics))
  expect_identical(des_leg$parameters$role, c("pre", "post", "post", "post"))

  # 3. Thin storage and explicit materialization bound on a 4,000-link qualified chain:
  m <- 4000L
  links_thin <- data.frame(
    g = rep(2, m),
    s = 0:(m - 1L),
    t = 1:m,
    ref = rep(1, m),
    functional_key = rep("op:thin", m),
    s_qualified = c(FALSE, TRUE, rep(c(TRUE, FALSE), length.out = m - 2L)),
    t_qualified = c(TRUE, rep(c(TRUE, FALSE), length.out = m - 1L)),
    stringsAsFactors = FALSE
  )
  des_thin <- cdid:::build_chained_design(links_thin)
  expect_null(des_thin$W)
  expect_null(des_thin$rank_diagnostics$root_paths)
  expect_null(des_thin$rank_diagnostics$path_witnesses)
  expect_null(des_thin$rank_diagnostics$nullspace_basis)
  expect_lt(as.double(object.size(des_thin)), 6e6)
  expect_identical(des_thin$rank_diagnostics$rank, m)
  expect_error(cdid:::materialize_chained_design(des_thin, 1000), class = "cdid_design_bound_error")

  # 4. Forced post-sharing perturbation fails independent hand/QR reference checks:
  perturbed_n3 <- links_large[, c("g", "s", "t", "ref", "functional_key", "s_qualified", "t_qualified")]
  perturbed_n3$s_qualified <- perturbed_n3$s == perturbed_n3$ref | perturbed_n3$s >= perturbed_n3$g
  perturbed_n3$t_qualified <- perturbed_n3$t == perturbed_n3$ref | perturbed_n3$t >= perturbed_n3$g
  des_pert <- qual_test_build(perturbed_n3, target_post_times = c(cal[3], cal[4]))
  W_expected_n3 <- matrix(c(1L, 0L, 0L, 1L, -1L, 0L), nrow = 2L, byrow = TRUE)
  expect_false(identical(dim(des_pert$W), dim(W_expected_n3)))
  expect_false(identical(des_pert$post_targets$level_estimable, c(TRUE, FALSE)))
})

test_that("permitted NA horizons, legacy unselected control_pool invariance, and trailing-comma pools are handled", {
  expect_typed_err <- function(df, cls) {
    expect_s3_class(expect_error(cdid:::build_chained_design(df), class = cls), "cdid_error")
  }
  x <- data.frame(
    g = 3, s = 2, t = 4, ref = 2, functional_key = "op:private",
    s_qualified = TRUE, t_qualified = FALSE, stringsAsFactors = FALSE
  )
  for (na in list(NA_real_, NA)) {
    y <- x; y$max_qualified_horizon <- na
    d <- qual_test_build(y)
    expect_identical(d$parameters$role, "post_nuisance")
    qual_verify_all_functionals(d, matrix(1L, 1L, 1L))
  }
  y <- rbind(x, x); y$functional_key <- c("op:A", "op:B"); y$max_qualified_horizon <- c(NA, 3)
  d_mix <- qual_test_build(y)
  expect_identical(nrow(d_mix$parameters), 2L)
  qual_verify_all_functionals(d_mix, diag(1L, 2L))

  y$functional_key <- "op:A"; y$s[2] <- 3; y$s_qualified[2] <- FALSE
  expect_typed_err(y, "cdid_design_horizon_error")
  y <- x; y$t_qualified <- TRUE; y$max_qualified_horizon <- NA_real_
  expect_typed_err(y, "cdid_design_horizon_error")
  for (bad in c(NaN, -Inf)) {
    y <- x; y$max_qualified_horizon <- bad
    expect_typed_err(y, "cdid_design_horizon_error")
  }

  legacy <- x[, c("g", "s", "t", "ref", "functional_key")]
  e <- new.env(parent = globalenv()); e$cdid_abort <- cdid:::cdid_abort
  sys.source(testthat::test_path("reference", "chained_design.R"), envir = e)
  base_build <- e$build_chained_design
  for (cp in list("{}", "{3}", NA_character_, factor("pool"), matrix("pool", 1, 1))) {
    y <- legacy; y$control_pool <- cp
    expect_identical(cdid:::build_chained_design(y, c(3, 4)), base_build(y, c(3, 4)))
  }
  y <- rbind(legacy, legacy); y$s <- c(2, 3); y$control_pool <- c("A", "B")
  expect_identical(cdid:::build_chained_design(y), base_build(y))
  y <- legacy; y[["extra_pool"]] <- "B"; names(y)[ncol(y)] <- "control_pool"
  y[["another_pool"]] <- "A"; names(y)[ncol(y)] <- "control_pool"
  expect_identical(cdid:::build_chained_design(y), base_build(y))

  for (cp in c("{4,}", "{4, }", "{4,,}", "{4,4}", "{3}", "{NaN}")) {
    y <- x; y$control_pool <- cp
    expect_typed_err(y, "cdid_design_qualification_error")
  }
})
