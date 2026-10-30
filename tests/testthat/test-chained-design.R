# Dense reference objects are test-only and bounded before allocation.
test_build_design <- function(...) cdid:::materialize_chained_design(
  cdid:::build_chained_design(...), max_entries = 100000)
test_unaugmented_design <- function(...) cdid:::materialize_chained_design(
  cdid:::build_unaugmented_time_design(...), max_entries = 100000)
dense_test_witness <- function(design, cert) {
  out <- numeric(nrow(design$W)); names(out) <- rownames(design$W)
  out[cert$link_index] <- cert$sign
  out
}

# Independent linear-algebraic reference check for a design matrix W:
# Determines whether a linear functional `c_vec' theta` lies in row(W)
# using QR decomposition of W' without relying on graph diagnostics.
matrix_ref_estimable <- function(W, c_vec, tol = 1e-12) {
  Wt <- t(as.matrix(W))
  qr_Wt <- qr(Wt, tol = tol)
  rk <- qr_Wt$rank
  if (rk == 0L) {
    return(all(abs(c_vec) <= tol))
  }
  # Residual of projecting c_vec onto col(W') = row(W).
  resid <- qr.resid(qr_Wt, as.numeric(c_vec))
  max(abs(resid)) <= tol
}

test_that("NT fixed comparison functional shares common free pre-vertices and anchors reference at zero", {
  # Irregular numeric dates (2000, 2002.5, 2005, 2009), cohort g = 2005, ref = 2002.5.
  links <- data.frame(
    g = rep(2005, 5L),
    s = c(2000, 2000, 2002.5, 2002.5, 2005),
    t = c(2002.5, 2005, 2005, 2009, 2009),
    ref = rep(2002.5, 5L),
    control_group = rep("nevertreated", 5L),
    stringsAsFactors = FALSE
  )
  # Declared from this fixed population operator, not generated from metadata.
  links$functional_key <- rep("nt:fixed", 5L)
  expect_identical(links$functional_key, rep("nt:fixed", 5L))

  des <- test_build_design(links, target_post_times = c(2005, 2009))
  expect_s3_class(des, "cdid_chained_design")

  # Reference vertex (2005, 2002.5) is anchored at 0; free parameters are
  # 1 common pre-vertex (t = 2000) and 2 post vertices (t = 2005, 2009).
  expect_identical(
    des$parameters$param_id,
    c("beta(g=2005,t=2000,fn=nt:fixed)",
      "ATT(g=2005,t=2005)",
      "ATT(g=2005,t=2009)")
  )
  expect_identical(des$parameters$role, c("pre", "post", "post"))

  W_expected <- matrix(
    c(-1L,  0L, 0L,
      -1L,  1L, 0L,
       0L,  1L, 0L,
       0L,  0L, 1L,
       0L, -1L, 1L),
    nrow = 5L, ncol = 3L, byrow = TRUE,
    dimnames = list(rownames(des$W), colnames(des$W))
  )
  expect_identical(des$W, W_expected)
  expect_identical(des$rank_diagnostics$rank, 3L)
  expect_identical(des$rank_diagnostics$nullity, 0L)
  expect_identical(des$rank_diagnostics$cycle_rank, 2L)
  expect_true(des$rank_diagnostics$has_full_column_rank)
  expect_identical(des$parameters$level_estimable, c(TRUE, TRUE, TRUE))

  # Path witnesses satisfy P' W = I_3 exactly.
  expect_identical(
    crossprod(des$rank_diagnostics$path_witnesses, des$W),
    matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3L, ncol = 3L,
           dimnames = list(colnames(des$W), colnames(des$W)))
  )
})

test_that("NYT counterexample: unaugmented all-pair projection gives ATT4 = 1/4 while augmented horizon-pre design gives 0", {
  # Population counterexample from cp02_nyt_identification.tex Section 1:
  # Dates 1..4, target g = 3, ref = 2, all true ATTs = 0.
  # Pairs in order: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4).
  # Horizons h = max(3, t) are (3, 3, 4, 3, 4, 4).
  # Population link vector is mu = (-1, -1, 0, 0, 0, 0).
  links_6 <- data.frame(
    g = rep(3, 6L),
    s = c(1, 1, 1, 2, 2, 3),
    t = c(2, 3, 4, 3, 4, 4),
    ref = rep(2, 6L),
    control_group = rep("notyettreated", 6L),
    horizon = c(3, 3, 4, 3, 4, 4),
    profile_id = rep("balanced", 6L),
    stringsAsFactors = FALSE
  )
  # The population counterexample supplies two distinct operators.
  links_6$functional_key <- paste0("nyt:h=", links_6$horizon, ":profile=balanced")
  mu_6 <- c(-1, -1, 0, 0, 0, 0)

  # 1. Unaugmented design collapses pre-coordinate across horizons 3 and 4,
  #    yielding (beta_1, ATT3, ATT4) = (3/4, 0, 1/4).
  des_unaug <- test_unaugmented_design(links_6)
  expect_false(des_unaug$is_valid_augmented)
  expect_identical(dim(des_unaug$W), c(6L, 3L))
  theta_unaug <- drop(qr.solve(des_unaug$W, mu_6))
  expect_equal(unname(theta_unaug), c(3 / 4, 0, 1 / 4), tolerance = 1e-12)

  # 3-link subset (1,2), (2,4), (1,4) with mu_3 = (-1, 0, 0):
  # unaugmented gives (2/3, 1/3).
  sub_idx <- c(1L, 5L, 3L)
  des_unaug_3 <- test_unaugmented_design(links_6[sub_idx, ])
  theta_unaug_3 <- drop(qr.solve(des_unaug_3$W, mu_6[sub_idx]))
  expect_equal(unname(theta_unaug_3), c(2 / 3, 1 / 3), tolerance = 1e-12)

  # 2. Valid augmented design gives distinct early-pre vertices for h=3 and h=4
  #    and shared post vertices ATT(3,3) and ATT(3,4).
  des_aug <- test_build_design(links_6)
  expect_identical(dim(des_aug$W), c(6L, 4L))
  expect_identical(
    des_aug$parameters$param_id,
    c("beta(g=3,t=1,fn=nyt:h=3:profile=balanced)",
      "beta(g=3,t=1,fn=nyt:h=4:profile=balanced)",
      "ATT(g=3,t=3)",
      "ATT(g=3,t=4)")
  )
  expect_identical(des_aug$rank_diagnostics$rank, 4L)
  expect_identical(des_aug$rank_diagnostics$nullity, 0L)
  expect_identical(des_aug$rank_diagnostics$cycle_rank, 2L)

  theta_aug <- drop(qr.solve(des_aug$W, mu_6))
  expect_equal(unname(theta_aug), c(1, 0, 0, 0), tolerance = 1e-12)
  # Exact fit with zero residuals and zero post ATTs:
  expect_equal(unname(drop(des_aug$W %*% theta_aug)), mu_6, tolerance = 1e-12)
  post_cols <- des_aug$parameters$role == "post"
  expect_equal(unname(theta_aug[post_cols]), c(0, 0), tolerance = 1e-12)

  # 3-link subset under augmented design has rank 3 and gives (1, 0, 0) -> ATT4 = 0.
  des_aug_3 <- test_build_design(links_6[sub_idx, ])
  expect_identical(dim(des_aug_3$W), c(3L, 3L))
  expect_identical(des_aug_3$rank_diagnostics$rank, 3L)
  theta_aug_3 <- drop(qr.solve(des_aug_3$W, mu_6[sub_idx]))
  expect_equal(unname(theta_aug_3), c(1, 0, 0), tolerance = 1e-12)
})

test_that("differing observation profiles at a fixed horizon cannot be falsely merged", {
  # Section 3 of cp02_nyt_identification.tex:
  # Fixed horizon h = 3, links (1,2), (2,3), (1,3).
  # Pair 12 has cohort observation rates (1, 1, 1); pairs 23 and 13 have (1, 1, 0.5).
  # Population link values are mu = (-1, 0, -4/3).
  links_fixed_h <- data.frame(
    g = c(3, 3, 3),
    s = c(1, 2, 1),
    t = c(2, 3, 3),
    ref = c(2, 2, 2),
    control_group = rep("notyettreated", 3L),
    horizon = c(3, 3, 3),
    stringsAsFactors = FALSE
  )
  rates_exact <- rbind(
    c(1.0, 1.0, 1.0),
    c(1.0, 1.0, 0.5),
    c(1.0, 1.0, 0.5)
  )
  mu_sec3 <- c(-1, 0, -4 / 3)

  # 1. A horizon supplies no population-functional declaration.
  expect_error(
    test_build_design(links_fixed_h),
    class = "cdid_design_profile_merge_error"
  )
  expect_error(
    cdid:::validate_functional_metadata(links_fixed_h),
    class = "cdid_design_input_error"
  )

  # 2. A bare horizon key cannot declare the full comparison operator.
  bad_merge_id <- links_fixed_h
  bad_merge_id$profile_id <- c("p_12", "p_23", "p_13")
  bad_merge_id$functional_key <- rep("nyt:h=3", 3L)
  expect_error(
    test_build_design(bad_merge_id),
    class = "cdid_design_profile_merge_error"
  )

  # Conflicting rate-profile descriptions under one declared key are rejected
  # as inconsistent metadata. Equal fitted rates would certify nothing.
  rates_close <- rbind(
    c(1.0, 1.0, 0.5000),
    c(1.0, 1.0, 0.5000),
    c(1.0, 1.0, 0.4999)
  )
  bad_close <- links_fixed_h
  bad_close$functional_key <- rep("nyt:h=3:merged", 3L)
  attr(bad_close, "profile_rates") <- rates_close
  expect_error(
    test_build_design(bad_close),
    class = "cdid_design_profile_merge_error"
  )

  # 3. Explicitly declared distinct keys for distinct uncorrected profiles at h = 3
  #    give separate early-pre vertices for (1,2) and (1,3), recovering ATT(3,3) = 0.
  links_distinct <- links_fixed_h
  # Equality here follows from the specified population law in the NYT note.
  links_distinct$functional_key <- c("population:A", "population:B", "population:B")
  expect_false(identical(links_distinct$functional_key[1L],
                         links_distinct$functional_key[3L]))
  expect_identical(links_distinct$functional_key[2L],
                   links_distinct$functional_key[3L])

  des_distinct <- test_build_design(links_distinct)
  expect_identical(dim(des_distinct$W), c(3L, 3L))
  expect_identical(des_distinct$rank_diagnostics$rank, 3L)
  theta_distinct <- drop(qr.solve(des_distinct$W, mu_sec3))
  expect_equal(unname(theta_distinct), c(1, 4 / 3, 0), tolerance = 1e-12)

  # By contrast, falsely collapsing into a single pre-coordinate at h = 3
  # produces the biased ATT(3,3) = -1/9.
  des_false_merge <- test_unaugmented_design(links_distinct)
  theta_false <- drop(qr.solve(des_false_merge$W, mu_sec3))
  expect_equal(unname(theta_false), c(10 / 9, -1 / 9), tolerance = 1e-12)

  # 4. After Proposition 2 inverse pair-rate correction restores the common
  #    comparison functional L_{3,3}, mu_corr = (-1, 0, -1) and the certified
  #    declared common operator shares the pre-vertex and gives ATT(3,3) = 0.
  links_corr <- links_fixed_h
  links_corr$composition_corrected <- rep(TRUE, 3L)
  # Proposition 2, including population sampling/rate assumptions, supplies
  # this operator identity. The boolean itself establishes nothing.
  links_corr$functional_key <- rep("nyt:h=3:certified_pop", 3L)
  expect_identical(links_corr$functional_key, rep("nyt:h=3:certified_pop", 3L))
  des_corr <- test_build_design(links_corr)
  expect_identical(dim(des_corr$W), c(3L, 2L))
  theta_corr <- drop(qr.solve(des_corr$W, c(-1, 0, -1)))
  expect_equal(unname(theta_corr), c(1, 0), tolerance = 1e-12)
})

test_that("enlarged graph detects isolated targets, floating components, and false time-graph connectivity", {
  # Two cohorts:
  # Cohort g = 3 (ref = 2):
  #   - link (2, 3): connects ATT(3,3) to ref=2 (estimable)
  #   - link (4, 6): connects ATT(3,4) and ATT(3,6) to each other, but
  #     disconnected from ref=2 (floating component: contrast ATT(3,6)-ATT(3,4)
  #     is estimable, levels are not)
  #   - target t = 5 requested in target_post_times with no links (isolated)
  # Cohort g = 4 (ref = 3):
  #   - link (1, 3) with fn = "nyt:h=4:pA"
  #   - link (1, 4) with fn = "nyt:h=4:pB" (different early-pre functional!)
  #   In the naive time-only graph, 3 -- 1 -- 4 connects t=4 to ref=3,
  #   but in G_aug, (1, pA) != (1, pB), so ATT(4,4) is disconnected from ref=3!
  links <- data.frame(
    g = c(3, 3, 4, 4),
    s = c(2, 4, 1, 1),
    t = c(3, 6, 3, 4),
    ref = c(2, 2, 3, 3),
    functional_key = c("nyt:h=3:p1", "nyt:h=6:p1", "nyt:h=4:pA", "nyt:h=4:pB"),
    stringsAsFactors = FALSE
  )
  req_targets <- data.frame(
    g = c(3, 3, 3, 3, 4),
    time = c(3, 4, 5, 6, 4)
  )

  des <- test_build_design(links, target_post_times = req_targets)

  # Check post_targets statuses:
  pt <- des$post_targets
  expect_identical(pt$g, c(3, 3, 3, 3, 4))
  expect_identical(pt$time, c(3, 4, 5, 6, 4))
  expect_identical(pt$is_isolated, c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(pt$is_supported, c(TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_identical(pt$level_estimable, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(
    pt$status,
    c("estimable", "disconnected_component", "isolated_unsupported",
      "disconnected_component", "disconnected_component")
  )

  # Isolated target ATT(3,5) has an all-zero column in W.
  expect_identical(unname(des$W[, "ATT(g=3,t=5)"]), rep(0L, 4L))

  # False connectivity in naive time-only graph for ATT(4,4) is flagged:
  att44_row <- pt[pt$g == 4 & pt$time == 4, ]
  expect_true(att44_row$time_graph_connected_to_ref)
  expect_false(att44_row$connected_to_ref)
  expect_false(att44_row$level_estimable)

  # Contrast estimability: ATT(3,6) - ATT(3,4) is estimable even though neither level is.
  c_46 <- cdid:::check_contrast_estimable(des, "ATT(g=3,t=4)", "ATT(g=3,t=6)")
  expect_true(c_46$estimable)
  e_diff <- numeric(ncol(des$W))
  names(e_diff) <- colnames(des$W)
  e_diff["ATT(g=3,t=6)"] <- 1
  e_diff["ATT(g=3,t=4)"] <- -1
  expect_identical(drop(crossprod(dense_test_witness(des, c_46$witness_links), des$W)), e_diff)

  # Contrast between disconnected components ATT(3,4) and ATT(3,3) is not estimable.
  c_34 <- cdid:::check_contrast_estimable(des, "ATT(g=3,t=3)", "ATT(g=3,t=4)")
  expect_false(c_34$estimable)

  # Independent matrix reference comparison for rank, nullity, nullspace, and estimability.
  W_mat <- des$W
  qr_rk <- qr(W_mat, tol = 1e-12)$rank
  gram_rk <- qr(crossprod(W_mat), tol = 1e-12)$rank
  expect_identical(des$rank_diagnostics$rank, as.integer(qr_rk))
  expect_identical(des$rank_diagnostics$rank, as.integer(gram_rk))
  expect_identical(des$rank_diagnostics$nullity, ncol(W_mat) - as.integer(qr_rk))

  # Nullspace basis columns satisfy W B = 0 and have full column rank = nullity.
  B <- des$rank_diagnostics$nullspace_basis
  expect_identical( ncol(B), des$rank_diagnostics$nullity )
  expect_true(all(W_mat %*% B == 0L))
  expect_identical(qr(B)$rank, des$rank_diagnostics$nullity)

  # Every parameter's graph level_estimable matches the independent QR row-space check.
  for (j in seq_len(ncol(W_mat))) {
    ej <- numeric(ncol(W_mat))
    ej[j] <- 1
    expect_identical(
      des$parameters$level_estimable[j],
      matrix_ref_estimable(W_mat, ej)
    )
  }

  # Every pair of parameters' contrast estimability matches the independent QR row-space check.
  for (i in seq_len(ncol(W_mat) - 1L)) {
    for (j in (i + 1L):ncol(W_mat)) {
      e_ij <- numeric(ncol(W_mat))
      e_ij[j] <- 1
      e_ij[i] <- -1
      c_ij <- cdid:::check_contrast_estimable(des, i, j)
      expect_identical(c_ij$estimable, matrix_ref_estimable(W_mat, e_ij))
    }
  }
})

test_that("build_chained_design handles multi-pre chains, link permutations, and rejects invalid link specs", {
  # Cohort g = 4, ref = 3, dates 1, 2, 3, 4 with pre-pre link (1, 2),
  # pre-ref link (2, 3), straddling link (1, 4), and ref-post link (3, 4)
  # under a declared population functional "nyt:h=4:certified_pop".
  links_mp <- data.frame(
    g = rep(4, 4L),
    s = c(1, 2, 1, 3),
    t = c(2, 3, 4, 4),
    ref = rep(3, 4L),
    functional_key = rep("nyt:h=4:certified_pop", 4L),
    stringsAsFactors = FALSE
  )
  des_mp <- test_build_design(links_mp)
  expect_identical(
    des_mp$parameters$param_id,
    c("beta(g=4,t=1,fn=nyt:h=4:certified_pop)",
      "beta(g=4,t=2,fn=nyt:h=4:certified_pop)",
      "ATT(g=4,t=4)")
  )
  expect_identical(des_mp$links$link_type, c("pre_only", "pre_only", "straddling", "post_only"))
  expect_identical(des_mp$rank_diagnostics$rank, 3L)
  expect_identical(des_mp$rank_diagnostics$cycle_rank, 1L)
  expect_true(all(des_mp$parameters$level_estimable))

  # Permuting link rows leaves parameter columns, vertices, post_targets,
  # rank diagnostics, and W'W strictly invariant.
  perm <- c(4L, 2L, 3L, 1L)
  des_perm <- test_build_design(links_mp[perm, ])
  expect_identical(des_perm$parameters, des_mp$parameters)
  expect_identical(des_perm$vertices, des_mp$vertices)
  expect_identical(des_perm$post_targets, des_mp$post_targets)
  expect_identical(des_perm$rank_diagnostics$rank, des_mp$rank_diagnostics$rank)
  expect_identical(des_perm$rank_diagnostics$cycle_rank, des_mp$rank_diagnostics$cycle_rank)
  expect_identical(crossprod(des_perm$W), crossprod(des_mp$W))

  # Invalid link specs are rejected with structured cdid_error classes:
  # 1. Bare horizon key on early-pre links
  bad_bare <- links_mp
  bad_bare$functional_key <- rep("nyt:h=4", 4L)
  expect_error(
    test_build_design(bad_bare),
    class = "cdid_design_profile_merge_error"
  )

  # 2. s >= t
  bad_st <- links_mp
  bad_st$s[1L] <- 2
  bad_st$t[1L] <- 2
  expect_error(
    test_build_design(bad_st),
    class = "cdid_design_endpoint_error"
  )

  # 3. Inconsistent ref within cohort or ref >= g
  bad_ref <- links_mp
  bad_ref$ref[2L] <- 2
  expect_error(
    test_build_design(bad_ref),
    class = "cdid_design_reference_error"
  )
  bad_ref2 <- links_mp
  bad_ref2$ref <- rep(4, 4L)
  expect_error(
    test_build_design(bad_ref2),
    class = "cdid_design_reference_error"
  )

  # 4. Control horizon < max(g, t) on post link
  bad_hor <- links_mp
  bad_hor$horizon <- c(4, 4, 3, 4)
  expect_error(
    test_build_design(bad_hor),
    class = "cdid_design_horizon_error"
  )

  # 5. Duplicate (g, s, t, functional_key)
  bad_dup <- rbind(links_mp, links_mp[1L, ])
  expect_error(
    test_build_design(bad_dup),
    class = "cdid_design_duplicate_link_error"
  )
})

