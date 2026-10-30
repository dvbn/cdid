# SPDX-License-Identifier: GPL-3.0-only
#
# Self-contained packaged tests for the CP02 private complete-target
# influence-function cluster covariance and supplied multiplier component
# (`R/chained_multiplier.R`).
#
# Uses only base R, `testthat`, and `digest`.
# In packaged execution (`testthat::test_check("cdid")`), resolves private
# functions directly from the actual installed `cdid` namespace (`asNamespace("cdid")`)
# and verifies non-export. The separate source-only harness never substitutes
# a sourced environment for these packaged namespace tests.

resolve_cdid_private_env <- function() {
  ns <- base::asNamespace("cdid")
  if (!base::isNamespace(ns)) stop("'cdid' is not a genuine R namespace.")
  list(env = ns, mode = "installed_namespace")
}

get_private_fn <- function(name) {
  ctx <- resolve_cdid_private_env()
  if (!base::exists(name, envir = ctx$env, inherits = FALSE)) {
    stop(sprintf("Private symbol '%s' not found in target environment.", name))
  }
  fn <- base::get(name, envir = ctx$env, inherits = FALSE)
  if (!base::is.function(fn) || !base::identical(base::environment(fn), ctx$env)) {
    stop(sprintf("Private symbol '%s' is not a closure in the actual cdid namespace.", name))
  }
  fn
}

make_explicit_budget <- function(max_units = 64,
                                 max_clusters = 32,
                                 max_targets = 16,
                                 max_draws = 128,
                                 max_if_batch_rows = 64,
                                 max_if_batch_cells = 1024,
                                 max_sign_batch_rows = 128,
                                 max_sign_batch_cells = 4096,
                                 max_workspace_cells = 16384,
                                 max_workspace_bytes = 131072) {
  list(
    max_units = max_units,
    max_clusters = max_clusters,
    max_targets = max_targets,
    max_draws = max_draws,
    max_if_batch_rows = max_if_batch_rows,
    max_if_batch_cells = max_if_batch_cells,
    max_sign_batch_rows = max_sign_batch_rows,
    max_sign_batch_cells = max_sign_batch_cells,
    max_workspace_cells = max_workspace_cells,
    max_workspace_bytes = max_workspace_bytes
  )
}

full_sign_hypercube <- function(n_clusters) {
  grids <- replicate(n_clusters, c(-1.0, 1.0), simplify = FALSE)
  mat <- as.matrix(do.call(expand.grid, grids))
  dimnames(mat) <- NULL
  storage.mode(mat) <- "double"
  mat
}

# Independent reference implementation from first principles for test assertions
independent_exact_reference <- function(unit_if_matrix,
                                        unit_ids,
                                        unit_cluster_map,
                                        cluster_ids,
                                        target_ids,
                                        n_denom,
                                        center,
                                        sign_matrix,
                                        alpha) {
  C <- length(cluster_ids)
  J <- length(target_ids)
  B <- nrow(sign_matrix)

  S <- matrix(0.0, nrow = C, ncol = J, dimnames = list(cluster_ids, target_ids))
  map_pos <- match(unit_ids, unit_cluster_map$unit_ids)
  u_clusters <- unit_cluster_map$cluster_ids[map_pos]
  for (i in seq_along(unit_ids)) {
    c_idx <- match(u_clusters[i], cluster_ids)
    S[c_idx, ] <- S[c_idx, ] + unit_if_matrix[i, ]
  }

  S_bar <- colSums(S) / C
  S_tilde <- if (center) {
    S - matrix(rep(S_bar, each = C), nrow = C, ncol = J)
  } else {
    S
  }

  G <- matrix(0.0, nrow = J, ncol = J, dimnames = list(target_ids, target_ids))
  for (c_idx in seq_len(C)) {
    row_c <- S_tilde[c_idx, ]
    G <- G + outer(row_c, row_c)
  }
  cov_mat <- G / (n_denom^2)
  omega_cl <- G / n_denom
  analytic_var <- diag(cov_mat)
  analytic_se <- sqrt(analytic_var)

  bres <- (sign_matrix %*% S_tilde) / sqrt(C)
  iqr_norm <- stats::qnorm(0.75) - stats::qnorm(0.25)
  k25 <- as.integer(ceiling(0.25 * B))
  k75 <- as.integer(ceiling(0.75 * B))
  kcrit <- as.integer(ceiling((1.0 - alpha) * B))

  q25 <- numeric(J)
  q75 <- numeric(J)
  iqr_raw <- numeric(J)
  bSigma <- numeric(J)
  names(q25) <- names(q75) <- names(iqr_raw) <- names(bSigma) <- target_ids
  for (j in seq_len(J)) {
    ord <- sort.int(bres[, j], method = "quick")
    q25[j] <- ord[k25]
    q75[j] <- ord[k75]
    iqr_raw[j] <- q75[j] - q25[j]
    bSigma[j] <- iqr_raw[j] / iqr_norm
  }
  boot_se <- bSigma * sqrt(C) / n_denom

  std_abs <- abs(bres) / matrix(rep(bSigma, each = B), nrow = B, ncol = J)
  max_stat <- apply(std_abs, 1L, max)
  crit_val <- sort.int(max_stat, method = "quick")[kcrit]

  list(
    S = S,
    S_bar = S_bar,
    S_tilde = S_tilde,
    cov_mat = cov_mat,
    omega_cl = omega_cl,
    analytic_var = analytic_var,
    analytic_se = analytic_se,
    bres = bres,
    q25 = q25,
    q75 = q75,
    iqr_raw = iqr_raw,
    bSigma = bSigma,
    boot_se = boot_se,
    max_stat = max_stat,
    crit_val = crit_val
  )
}

build_fixture_c4 <- function() {
  target_ids <- c("att_g2_t2", "att_g2_t3")
  cluster_ids <- c("cl_d", "cl_a", "cl_c", "cl_b")
  expected_unit_ids <- c("u1", "u2", "u3", "u4", "u5", "u6")
  unit_cluster_map <- list(
    unit_ids = c("u6", "u2", "u4", "u1", "u5", "u3"),
    cluster_ids = c("cl_c", "cl_a", "cl_d", "cl_d", "cl_c", "cl_c")
  )
  if_rows <- rbind(
    u1 = c(2.0, 1.0),
    u2 = c(-1.0, 8.0),
    u3 = c(4.0, -5.0),
    u4 = c(3.0, 3.0),
    u5 = c(-1.0, 1.0),
    u6 = c(-1.0, 2.0)
  )
  colnames(if_rows) <- target_ids
  if_batches <- list(
    list(
      unit_ids = c("u5", "u1"),
      target_ids = target_ids,
      if_matrix = unname(if_rows[c("u5", "u1"), , drop = FALSE])
    ),
    list(
      unit_ids = c("u3", "u2", "u6"),
      target_ids = target_ids,
      if_matrix = unname(if_rows[c("u3", "u2", "u6"), , drop = FALSE])
    ),
    list(
      unit_ids = c("u4"),
      target_ids = target_ids,
      if_matrix = unname(if_rows["u4", , drop = FALSE])
    )
  )
  signs <- full_sign_hypercube(4L)
  draw_ids <- sprintf("draw_%02d", seq_len(nrow(signs)))
  sign_batches <- list(
    list(
      draw_ids = draw_ids[1:6],
      cluster_ids = cluster_ids,
      sign_matrix = signs[1:6, , drop = FALSE]
    ),
    list(
      draw_ids = draw_ids[7:16],
      cluster_ids = cluster_ids,
      sign_matrix = signs[7:16, , drop = FALSE]
    )
  )
  list(
    target_ids = target_ids,
    expected_unit_ids = expected_unit_ids,
    unit_cluster_map = unit_cluster_map,
    cluster_ids = cluster_ids,
    if_rows = if_rows,
    if_batches = if_batches,
    signs = signs,
    draw_ids = draw_ids,
    sign_batches = sign_batches,
    n_denom = 8.0,
    alpha = 0.10,
    budget = make_explicit_budget()
  )
}

test_that("private multiplier symbols reside in target namespace, have no public exports, and require all 15 explicit arguments", {
  ctx <- resolve_cdid_private_env()
  private_symbols <- c(
    "chained_cluster_multiplier",
    "compute_multiplier_workspace_estimate",
    "validate_multiplier_budget",
    "type1_order_index"
  )
  for (sym in private_symbols) {
    expect_true(exists(sym, envir = ctx$env, inherits = FALSE), info = sym)
  }
  if (identical(ctx$mode, "installed_namespace")) {
    ns <- asNamespace("cdid")
    exports <- getNamespaceExports(ns)
    for (sym in private_symbols) {
      expect_false(sym %in% exports, info = sym)
    }
  }

  fn <- get_private_fn("chained_cluster_multiplier")
  fmls <- formals(fn)
  expect_identical(
    names(fmls),
    c(
      "target_ids",
      "expected_unit_ids",
      "unit_cluster_map",
      "cluster_ids",
      "n_denom",
      "center",
      "if_batch_fn",
      "n_if_batches",
      "sign_batch_fn",
      "n_sign_batches",
      "expected_draw_ids",
      "n_draws",
      "alpha",
      "return_raw_draws",
      "budget"
    )
  )
  # Every formal argument must be a missing symbol (no defaults).
  for (nm in names(fmls)) {
    expect_true(is.symbol(fmls[[nm]]) && identical(as.character(fmls[[nm]]), ""), info = nm)
  }

  fx <- build_fixture_c4()
  base_args <- list(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) fx$if_batches[[i]],
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = nrow(fx$signs),
    alpha = fx$alpha,
    return_raw_draws = FALSE,
    budget = fx$budget
  )
  for (nm in names(base_args)) {
    incomplete_args <- base_args[setdiff(names(base_args), nm)]
    expect_error(
      do.call(fn, incomplete_args),
      class = "cdid_multiplier_missing_arg_error",
      info = paste("Missing argument:", nm)
    )
  }
})

test_that("independent exact small examples with all 2^C signs for C = 2, 3, 4, 5, 6 match exact covariance, scale, raw draws, and criticals", {
  fn <- get_private_fn("chained_cluster_multiplier")
  tol_transcendental <- 1e-15

  # --- Case 1: C = 2 (B = 4), J = 2, n = 4 ---
  t2 <- c("t1", "t2")
  c2 <- c("c2", "c1")
  u2 <- c("u1", "u2", "u3")
  map2 <- list(unit_ids = c("u3", "u1", "u2"), cluster_ids = c("c1", "c2", "c2"))
  if2 <- matrix(c(3.0, -1.0, 1.0, 5.0, -2.0, 2.0), nrow = 3, ncol = 2, byrow = TRUE)
  signs2 <- full_sign_hypercube(2L)
  draws2 <- sprintf("d%d", seq_len(nrow(signs2)))

  for (cent in c(TRUE, FALSE)) {
    res2 <- fn(
      target_ids = t2,
      expected_unit_ids = u2,
      unit_cluster_map = map2,
      cluster_ids = c2,
      n_denom = 4.0,
      center = cent,
      if_batch_fn = function(i) list(unit_ids = u2, target_ids = t2, if_matrix = if2),
      n_if_batches = 1L,
      sign_batch_fn = function(i) list(draw_ids = draws2, cluster_ids = c2, sign_matrix = signs2),
      n_sign_batches = 1L,
      expected_draw_ids = draws2,
      n_draws = 4L,
      alpha = 0.25,
      return_raw_draws = TRUE,
      budget = make_explicit_budget()
    )
    ref2 <- independent_exact_reference(if2, u2, map2, c2, t2, 4.0, cent, signs2, 0.25)
    expect_true(res2$accepted)
    expect_equal(res2$cluster_sums, ref2$S, tolerance = 0)
    expect_equal(res2$cluster_sums_centered, ref2$S_tilde, tolerance = 0)
    expect_equal(res2$cov_mat, ref2$cov_mat, tolerance = 0)
    expect_equal(unname(res2$raw_draws), unname(ref2$bres), tolerance = 0)
    expect_equal(res2$q25, ref2$q25, tolerance = 0)
    expect_equal(res2$q75, ref2$q75, tolerance = 0)
    expect_equal(res2$raw_bootstrap_iqr, ref2$iqr_raw, tolerance = 0)
    expect_equal(res2$analytic_se, ref2$analytic_se, tolerance = tol_transcendental)
    expect_equal(res2$bootstrap_se, ref2$boot_se, tolerance = tol_transcendental)
    expect_equal(res2$crit_val, ref2$crit_val, tolerance = tol_transcendental)
  }

  # --- Case 2: C = 3 (B = 8), J = 2, n = 8, exact rank-1 centering identity ---
  t3 <- c("tau_a", "tau_b")
  c3 <- c("k3", "k1", "k2")
  u3 <- c("u1", "u2", "u3", "u4", "u5")
  map3 <- list(
    unit_ids = u3,
    cluster_ids = c("k3", "k1", "k3", "k2", "k1")
  )
  if3 <- matrix(
    c(
      4.0, -2.0,
      1.0,  3.0,
      2.0,  6.0,
     -3.0,  1.0,
      5.0, -2.0
    ),
    nrow = 5, ncol = 2, byrow = TRUE
  )
  signs3 <- full_sign_hypercube(3L)
  draws3 <- sprintf("b%d", seq_len(nrow(signs3)))
  res3_cent <- fn(
    target_ids = t3,
    expected_unit_ids = u3,
    unit_cluster_map = map3,
    cluster_ids = c3,
    n_denom = 8.0,
    center = TRUE,
    if_batch_fn = function(i) list(unit_ids = u3, target_ids = t3, if_matrix = if3),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(draw_ids = draws3, cluster_ids = c3, sign_matrix = signs3),
    n_sign_batches = 1L,
    expected_draw_ids = draws3,
    n_draws = 8L,
    alpha = 0.10,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )
  res3_uncent <- fn(
    target_ids = t3,
    expected_unit_ids = u3,
    unit_cluster_map = map3,
    cluster_ids = c3,
    n_denom = 8.0,
    center = FALSE,
    if_batch_fn = function(i) list(unit_ids = u3, target_ids = t3, if_matrix = if3),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(draw_ids = draws3, cluster_ids = c3, sign_matrix = signs3),
    n_sign_batches = 1L,
    expected_draw_ids = draws3,
    n_draws = 8L,
    alpha = 0.10,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )
  # Verify Theorem 5(c) identity: Omega^0 - Omega^cent = (n * C)^{-1} t t'
  t_total <- colSums(res3_cent$cluster_sums)
  expect_equal(
    res3_uncent$omega_cl_hat - res3_cent$omega_cl_hat,
    outer(t_total, t_total) / (8.0 * 3.0),
    tolerance = 1e-15
  )

  # --- Case 3: C = 4 (B = 16), exact dyadic rationals ---
  fx4 <- build_fixture_c4()
  res4 <- fn(
    target_ids = fx4$target_ids,
    expected_unit_ids = fx4$expected_unit_ids,
    unit_cluster_map = fx4$unit_cluster_map,
    cluster_ids = fx4$cluster_ids,
    n_denom = fx4$n_denom,
    center = TRUE,
    if_batch_fn = function(i) fx4$if_batches[[i]],
    n_if_batches = length(fx4$if_batches),
    sign_batch_fn = function(i) fx4$sign_batches[[i]],
    n_sign_batches = length(fx4$sign_batches),
    expected_draw_ids = fx4$draw_ids,
    n_draws = 16L,
    alpha = fx4$alpha,
    return_raw_draws = TRUE,
    budget = fx4$budget
  )
  # Hand-verified exact dyadic values:
  # S = [[5, 4], [-1, 8], [2, -2], [0, 0]], mean = [1.5, 2.5]
  # S_tilde = [[3.5, 1.5], [-2.5, 5.5], [0.5, -4.5], [-1.5, -2.5]]
  # G = [[21, -7], [-7, 59]], cov_mat = G / 64
  expected_S4 <- matrix(
    c(5, 4, -1, 8, 2, -2, 0, 0),
    nrow = 4, ncol = 2, byrow = TRUE,
    dimnames = list(fx4$cluster_ids, fx4$target_ids)
  )
  expected_S4_tilde <- matrix(
    c(3.5, 1.5, -2.5, 5.5, 0.5, -4.5, -1.5, -2.5),
    nrow = 4, ncol = 2, byrow = TRUE,
    dimnames = list(fx4$cluster_ids, fx4$target_ids)
  )
  expected_cov4 <- matrix(
    c(21 / 64, -7 / 64, -7 / 64, 59 / 64),
    nrow = 2, ncol = 2, byrow = TRUE,
    dimnames = list(fx4$target_ids, fx4$target_ids)
  )
  expect_equal(res4$cluster_sums, expected_S4, tolerance = 0)
  expect_equal(res4$cluster_sums_centered, expected_S4_tilde, tolerance = 0)
  expect_equal(res4$cov_mat, expected_cov4, tolerance = 0)
  expect_equal(res4$q25, c(att_g2_t2 = -2.0, att_g2_t3 = -3.0), tolerance = 0)
  expect_equal(res4$q75, c(att_g2_t2 = 1.5, att_g2_t3 = 2.5), tolerance = 0)
  expect_equal(res4$raw_bootstrap_iqr, c(att_g2_t2 = 3.5, att_g2_t3 = 5.5), tolerance = 0)

  # --- Case 4: C = 5 (B = 32), J = 3, n = 10 ---
  t5 <- c("g1", "g2", "g3")
  c5 <- sprintf("cl_%d", 1:5)
  u5 <- sprintf("unit_%02d", 1:8)
  map5 <- list(
    unit_ids = u5,
    cluster_ids = c("cl_1", "cl_2", "cl_1", "cl_3", "cl_4", "cl_3", "cl_5", "cl_2")
  )
  if5 <- matrix(
    c(
       2.0, -1.0,  0.5,
      -3.0,  4.0,  1.5,
       1.0,  2.0, -2.0,
       0.5, -2.5,  3.0,
      -4.0,  1.0, -1.0,
       2.5, -0.5,  1.0,
       3.0, -3.0, -2.5,
      -1.0,  1.0,  0.5
    ),
    nrow = 8, ncol = 3, byrow = TRUE
  )
  signs5 <- full_sign_hypercube(5L)
  draws5 <- sprintf("d5_%02d", 1:32)
  res5 <- fn(
    target_ids = t5,
    expected_unit_ids = u5,
    unit_cluster_map = map5,
    cluster_ids = c5,
    n_denom = 10.0,
    center = TRUE,
    if_batch_fn = function(i) {
      idx <- if (i == 1L) c(8L, 2L, 5L) else if (i == 2L) c(1L, 6L, 3L) else c(7L, 4L)
      list(unit_ids = u5[idx], target_ids = t5, if_matrix = if5[idx, , drop = FALSE])
    },
    n_if_batches = 3L,
    sign_batch_fn = function(i) {
      idx <- if (i == 1L) 1:16 else 17:32
      list(draw_ids = draws5[idx], cluster_ids = c5, sign_matrix = signs5[idx, , drop = FALSE])
    },
    n_sign_batches = 2L,
    expected_draw_ids = draws5,
    n_draws = 32L,
    alpha = 0.05,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )
  ref5 <- independent_exact_reference(if5, u5, map5, c5, t5, 10.0, TRUE, signs5, 0.05)
  expect_true(res5$accepted)
  expect_equal(res5$cluster_sums, ref5$S, tolerance = 0)
  expect_equal(res5$cov_mat, ref5$cov_mat, tolerance = tol_transcendental)
  expect_equal(res5$bootstrap_se, ref5$boot_se, tolerance = tol_transcendental)
  expect_equal(res5$crit_val, ref5$crit_val, tolerance = tol_transcendental)

  # --- Case 5: C = 6 (B = 64), J = 2, n = 16, alpha = 0.05 ---
  t6 <- c("target_1", "target_2")
  c6 <- sprintf("c%d", 1:6)
  u6 <- sprintf("u%02d", 1:9)
  map6 <- list(
    unit_ids = u6,
    cluster_ids = c("c1", "c2", "c3", "c4", "c5", "c6", "c1", "c3", "c6")
  )
  if6 <- matrix(
    c(
       5.0,  0.5,
      -4.0,  3.0,
       1.0, -2.0,
      -1.0,  2.0,
       3.0, -4.0,
      -5.0,  1.0,
       3.0,  0.5,
       1.0, -3.0,
      -3.0,  2.0
    ),
    nrow = 9, ncol = 2, byrow = TRUE
  )
  signs6 <- full_sign_hypercube(6L)
  draws6 <- sprintf("d6_%02d", 1:64)
  res6 <- fn(
    target_ids = t6,
    expected_unit_ids = u6,
    unit_cluster_map = map6,
    cluster_ids = c6,
    n_denom = 16.0,
    center = TRUE,
    if_batch_fn = function(i) list(unit_ids = u6, target_ids = t6, if_matrix = if6),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(draw_ids = draws6, cluster_ids = c6, sign_matrix = signs6),
    n_sign_batches = 1L,
    expected_draw_ids = draws6,
    n_draws = 64L,
    alpha = 0.05,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )
  ref6 <- independent_exact_reference(if6, u6, map6, c6, t6, 16.0, TRUE, signs6, 0.05)
  expect_true(res6$accepted)
  expect_equal(res6$cluster_sums, ref6$S, tolerance = 0)
  expect_equal(res6$cov_mat, ref6$cov_mat, tolerance = tol_transcendental)
  expect_equal(res6$bootstrap_se, ref6$boot_se, tolerance = tol_transcendental)
  expect_equal(res6$crit_val, ref6$crit_val, tolerance = tol_transcendental)
})

test_that("arbitrary unit order, fragmented cluster batches, and integer IDs join by identity while keeping zero-unit and all-zero clusters", {
  fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  # Add a 7th unit in a 5th cluster whose IF row is explicitly (0, 0), while cl_b has 0 units.
  cluster_ids <- c("cl_d", "cl_a", "cl_zero_row", "cl_c", "cl_b")
  expected_unit_ids <- c("u1", "u2", "u3", "u4", "u5", "u6", "u7")
  if_rows <- rbind(fx$if_rows, u7 = c(0.0, 0.0))
  map_canonical <- list(
    unit_ids = expected_unit_ids,
    cluster_ids = c("cl_d", "cl_a", "cl_c", "cl_d", "cl_c", "cl_c", "cl_zero_row")
  )
  signs5 <- full_sign_hypercube(5L)
  draws5 <- sprintf("d_%02d", 1:32)

  # Single-batch canonical order
  res_canon <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = expected_unit_ids,
    unit_cluster_map = map_canonical,
    cluster_ids = cluster_ids,
    n_denom = 8.0,
    center = TRUE,
    if_batch_fn = function(i) list(
      unit_ids = expected_unit_ids,
      cluster_ids = map_canonical$cluster_ids,
      target_ids = fx$target_ids,
      if_matrix = unname(if_rows)
    ),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(
      draw_ids = draws5,
      cluster_ids = cluster_ids,
      sign_matrix = signs5
    ),
    n_sign_batches = 1L,
    expected_draw_ids = draws5,
    n_draws = 32L,
    alpha = 0.10,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )

  # Multi-batch scrambled unit order and scrambled unit_cluster_map with late entrant u1 in final batch
  perm_map <- c(7L, 5L, 2L, 6L, 4L, 3L, 1L)
  b1 <- c(6L, 2L)
  b2 <- c(7L, 4L)
  b3 <- c(5L, 3L)
  b4 <- c(1L) # late entrant belonging to cl_d
  batches <- list(b1, b2, b3, b4)

  res_scrambled <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = expected_unit_ids,
    unit_cluster_map = list(
      unit_ids = map_canonical$unit_ids[perm_map],
      cluster_ids = map_canonical$cluster_ids[perm_map]
    ),
    cluster_ids = cluster_ids,
    n_denom = 8.0,
    center = TRUE,
    if_batch_fn = function(i) {
      idx <- batches[[i]]
      list(
        unit_ids = expected_unit_ids[idx],
        cluster_ids = map_canonical$cluster_ids[idx],
        target_ids = fx$target_ids,
        if_matrix = unname(if_rows[idx, , drop = FALSE])
      )
    },
    n_if_batches = 4L,
    sign_batch_fn = function(i) {
      idx <- if (i == 1L) 1:11 else if (i == 2L) 12:23 else 24:32
      s_int <- signs5[idx, , drop = FALSE]
      storage.mode(s_int) <- "integer"
      list(
        draw_ids = draws5[idx],
        cluster_ids = cluster_ids,
        sign_matrix = s_int
      )
    },
    n_sign_batches = 3L,
    expected_draw_ids = draws5,
    n_draws = 32L,
    alpha = 0.10,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )

  expect_identical(res_scrambled$cluster_unit_counts, c(cl_d = 2L, cl_a = 1L, cl_zero_row = 1L, cl_c = 3L, cl_b = 0L))
  expect_identical(res_scrambled$metadata$n_zero_unit_clusters, 1L)
  expect_identical(res_scrambled$metadata$n_zero_sum_clusters, 2L)
  expect_equal(res_scrambled$cluster_sums, res_canon$cluster_sums, tolerance = 0)
  expect_equal(res_scrambled$cov_mat, res_canon$cov_mat, tolerance = 0)
  expect_equal(res_scrambled$analytic_se, res_canon$analytic_se, tolerance = 0)
  expect_equal(res_scrambled$bootstrap_se, res_canon$bootstrap_se, tolerance = 0)
  expect_equal(res_scrambled$crit_val, res_canon$crit_val, tolerance = 0)
  expect_equal(res_scrambled$raw_draws, res_canon$raw_draws, tolerance = 0)
})

test_that("mean-versus-sum, unit-vs-cluster centering, double-weighting, denominator n, and C-vs-n scaling conventions are strictly enforced", {
  fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  res_base <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = 8.0,
    center = TRUE,
    if_batch_fn = function(i) fx$if_batches[[i]],
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = TRUE,
    budget = fx$budget
  )

  # 1. Within-cluster sum vs within-cluster mean:
  # Cluster cl_d has 2 units summing to (5, 4) -> mean would be (2.5, 2.0).
  # Cluster cl_c has 3 units summing to (2, -2) -> mean would be (2/3, -2/3).
  wrong_within_cluster_mean <- rbind(
    cl_d = c(2.5, 2.0),
    cl_a = c(-1.0, 8.0),
    cl_c = c(2.0 / 3.0, -2.0 / 3.0),
    cl_b = c(0.0, 0.0)
  )
  expect_gt(max(abs(res_base$cluster_sums - wrong_within_cluster_mean)), 1.0)

  # 2. Cluster-level centering vs unit-level centering before cluster summation:
  # Unit mean is colSums(if_rows) / 6 = (6/6, 10/6) = (1, 5/3).
  # Subtracting unit mean per unit subtracts m_c * unit_mean from cluster c!
  unit_mean <- colMeans(fx$if_rows)
  m_c <- c(2, 1, 3, 0)
  wrong_unit_centered_S <- res_base$cluster_sums - outer(m_c, unit_mean)
  expect_gt(max(abs(res_base$cluster_sums_centered - wrong_unit_centered_S)), 0.5)
  wrong_unit_centered_cov <- crossprod(wrong_unit_centered_S) / 64.0
  expect_gt(max(abs(res_base$cov_mat - wrong_unit_centered_cov)), 0.05)

  # 3. Double-weighting refusal on unit_cluster_map and on IF batches:
  bad_map <- fx$unit_cluster_map
  bad_map$weights <- rep(1.0, length(fx$expected_unit_ids))
  expect_error(
    fn(
      target_ids = fx$target_ids,
      expected_unit_ids = fx$expected_unit_ids,
      unit_cluster_map = bad_map,
      cluster_ids = fx$cluster_ids,
      n_denom = 8.0,
      center = TRUE,
      if_batch_fn = function(i) fx$if_batches[[i]],
      n_if_batches = length(fx$if_batches),
      sign_batch_fn = function(i) fx$sign_batches[[i]],
      n_sign_batches = length(fx$sign_batches),
      expected_draw_ids = fx$draw_ids,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = FALSE,
      budget = fx$budget
    ),
    class = "cdid_multiplier_double_weight_error"
  )

  expect_error(
    fn(
      target_ids = fx$target_ids,
      expected_unit_ids = fx$expected_unit_ids,
      unit_cluster_map = fx$unit_cluster_map,
      cluster_ids = fx$cluster_ids,
      n_denom = 8.0,
      center = TRUE,
      if_batch_fn = function(i) {
        b <- fx$if_batches[[i]]
        b$unit_weights <- rep(1.0, length(b$unit_ids))
        b
      },
      n_if_batches = length(fx$if_batches),
      sign_batch_fn = function(i) fx$sign_batches[[i]],
      n_sign_batches = length(fx$sign_batches),
      expected_draw_ids = fx$draw_ids,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = FALSE,
      budget = fx$budget
    ),
    class = "cdid_multiplier_double_weight_error"
  )

  # 4. Normalization denominator n_denom scaling (n_denom = 16 vs n_denom = 8, C = 4, n_units = 6):
  res_n16 <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = 16.0,
    center = TRUE,
    if_batch_fn = function(i) fx$if_batches[[i]],
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = TRUE,
    budget = fx$budget
  )
  expect_equal(res_n16$cov_mat, res_base$cov_mat / 4.0, tolerance = 0)
  expect_equal(res_n16$omega_cl_hat, res_base$omega_cl_hat / 2.0, tolerance = 0)
  expect_equal(res_n16$analytic_se, res_base$analytic_se / 2.0, tolerance = 1e-15)
  expect_equal(res_n16$raw_draws, res_base$raw_draws, tolerance = 0)
  expect_equal(res_n16$bootstrap_iqr_scale, res_base$bootstrap_iqr_scale, tolerance = 0)
  expect_equal(res_n16$bootstrap_se, res_base$bootstrap_se / 2.0, tolerance = 1e-15)
  expect_equal(res_n16$crit_val, res_base$crit_val, tolerance = 0)
})

test_that("type-1 empirical quantiles for IQR and critical value are used and strictly distinguished from R default type-7", {
  fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  res <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) fx$if_batches[[i]],
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = 0.10,
    return_raw_draws = TRUE,
    budget = fx$budget
  )

  # On fixture C=4, B=16:
  # Type-1 IQR is (3.5, 5.5), whereas R default type-7 IQR is (3.25, 5.25).
  type7_iqr <- apply(res$raw_draws, 2L, function(col) {
    unname(diff(stats::quantile(col, probs = c(0.25, 0.75), type = 7)))
  })
  expect_equal(unname(res$raw_bootstrap_iqr), c(3.5, 5.5), tolerance = 0)
  expect_equal(unname(type7_iqr), c(3.25, 5.25), tolerance = 0)
  expect_gt(min(abs(res$raw_bootstrap_iqr - type7_iqr)), 0.24)

  # Type-1 critical value at alpha = 0.10 (index ceil(0.90 * 16) = 15) is ~1.716883,
  # whereas R default type-7 critical value on the same max_stat is ~1.629287.
  type7_crit <- unname(stats::quantile(res$max_stat, probs = 0.90, type = 7))
  expect_gt(abs(res$crit_val - type7_crit), 0.08)
})

test_that("scaling all IF rows by 1e-8 or mixed 1e-8 retains all target dimensions and normalized criticals without upstream absolute screen dropping", {
  fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  res_unit <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) fx$if_batches[[i]],
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = TRUE,
    budget = fx$budget
  )

  scale_factor <- 1e-8
  res_tiny <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      b$if_matrix <- b$if_matrix * scale_factor
      b
    },
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = TRUE,
    budget = fx$budget
  )

  # Confirm upstream mboot.R screen `colSums(bres^2) > sqrt(.Machine$double.eps) * 10`
  # would have dropped both targets here:
  expect_true(all(colSums(res_tiny$raw_draws^2) < sqrt(.Machine$double.eps) * 10))

  expect_true(res_tiny$accepted)
  expect_false(res_tiny$refused)
  expect_identical(res_tiny$status, "ok")
  expect_identical(res_tiny$target_status, c(att_g2_t2 = "ok", att_g2_t3 = "ok"))
  expect_identical(length(res_tiny$bootstrap_se), 2L)
  expect_identical(dim(res_tiny$cov_mat), c(2L, 2L))
  expect_equal(res_tiny$cov_mat, res_unit$cov_mat * (scale_factor^2), tolerance = 1e-15)
  expect_equal(res_tiny$analytic_se, res_unit$analytic_se * scale_factor, tolerance = 1e-15)
  expect_equal(res_tiny$bootstrap_se, res_unit$bootstrap_se * scale_factor, tolerance = 1e-15)
  expect_equal(res_tiny$max_stat, res_unit$max_stat, tolerance = 1e-14)
  expect_equal(res_tiny$crit_val, res_unit$crit_val, tolerance = 1e-14)

  # Mixed scaling: target 1 at scale 1, target 2 at scale 1e-9 (so sum(bres[,2]^2) < 10 * .Machine$double.eps)
  res_mixed <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      b$if_matrix[, 2L] <- b$if_matrix[, 2L] * 1e-9
      b
    },
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = TRUE,
    budget = fx$budget
  )
  expect_true(sum(res_mixed$raw_draws[, 2L]^2) < 10 * .Machine$double.eps)
  expect_true(res_mixed$accepted)
  expect_equal(res_mixed$max_stat, res_unit$max_stat, tolerance = 1e-14)
  expect_equal(res_mixed$crit_val, res_unit$crit_val, tolerance = 1e-14)
})

test_that("singular joint covariance with positive marginal variances is accepted as valid", {
  fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  # 3 targets where target 3 == target 1 + target 2 (rank(cov_mat) == 2 < 3), all marginal variances > 0
  t3 <- c("t1", "t2", "t1_plus_t2")
  if3 <- cbind(fx$if_rows, fx$if_rows[, 1L] + fx$if_rows[, 2L])
  colnames(if3) <- t3

  res <- fn(
    target_ids = t3,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) list(
      unit_ids = fx$expected_unit_ids,
      target_ids = t3,
      if_matrix = unname(if3)
    ),
    n_if_batches = 1L,
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = TRUE,
    budget = fx$budget
  )

  expect_true(res$accepted)
  expect_identical(res$status, "ok")
  expect_identical(res$target_status, c(t1 = "ok", t2 = "ok", t1_plus_t2 = "ok"))
  expect_true(all(res$analytic_var > 0))
  expect_true(all(is.finite(res$analytic_se) & res$analytic_se > 0))
  expect_true(all(is.finite(res$bootstrap_se) & res$bootstrap_se > 0))
  expect_lt(qr(res$cov_mat)$rank, 3L)
  expect_true(is.finite(res$crit_val) && res$crit_val > 0)
})

test_that("zero/nonfinite bootstrap IQR scale and zero/nonfinite analytical variance refuse complete target set without dropping dimensions or substituting analytic SE", {
  fn <- get_private_fn("chained_cluster_multiplier")

  # 1. Positive analytical variance (sum(S_c^2) > 0 and sum(bres^2) > 0) but ZERO bootstrap IQR on target 2:
  # C = 4 clusters, target 1 has S = (4, -2, 1, -3), target 2 has S = (1, -1, 0, 0) (centered).
  # Across 8 supplied sign draws where 6 draws have eta_1 == eta_2 and 2 draws have eta_1 == -eta_2:
  # target 2 raw draws sorted are (-1, 0, 0, 0, 0, 0, 0, 1), so q25 (index 2) == 0 and q75 (index 6) == 0!
  c4 <- c("c1", "c2", "c3", "c4")
  u4 <- c("u1", "u2", "u3", "u4")
  t2 <- c("valid_target", "zero_iqr_target")
  map4 <- list(unit_ids = u4, cluster_ids = c4)
  if_zero_iqr <- matrix(
    c(
       4.0,  1.0,
      -2.0, -1.0,
       1.0,  0.0,
      -3.0,  0.0
    ),
    nrow = 4, ncol = 2, byrow = TRUE
  )
  grid4 <- full_sign_hypercube(4L)
  # Select 6 rows where col1 == col2 and 2 rows where col1 != col2:
  same_idx <- which(grid4[, 1L] == grid4[, 2L])[1:6]
  diff_idx <- which(grid4[, 1L] != grid4[, 2L])[1:2]
  signs8 <- grid4[c(same_idx, diff_idx), , drop = FALSE]
  draws8 <- sprintf("d%d", 1:8)

  res_ziqr <- fn(
    target_ids = t2,
    expected_unit_ids = u4,
    unit_cluster_map = map4,
    cluster_ids = c4,
    n_denom = 4.0,
    center = TRUE,
    if_batch_fn = function(i) list(unit_ids = u4, target_ids = t2, if_matrix = if_zero_iqr),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(draw_ids = draws8, cluster_ids = c4, sign_matrix = signs8),
    n_sign_batches = 1L,
    expected_draw_ids = draws8,
    n_draws = 8L,
    alpha = 0.05,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )

  expect_false(res_ziqr$accepted)
  expect_true(res_ziqr$refused)
  expect_identical(res_ziqr$status, "refused")
  expect_identical(res_ziqr$target_status, c(valid_target = "ok", zero_iqr_target = "zero_bootstrap_scale"))
  expect_identical(res_ziqr$target_accepted, c(valid_target = TRUE, zero_iqr_target = FALSE))
  expect_identical(res_ziqr$diagnostics$analytic_variance_ok, c(valid_target = TRUE, zero_iqr_target = TRUE))
  expect_identical(res_ziqr$diagnostics$bootstrap_scale_ok, c(valid_target = TRUE, zero_iqr_target = FALSE))
  # Dimension 2 is preserved everywhere; analytic_se is positive on both, but bootstrap_se and crit_val are refused (NA_real_):
  expect_identical(length(res_ziqr$analytic_se), 2L)
  expect_true(all(res_ziqr$analytic_se > 0))
  expect_identical(length(res_ziqr$bootstrap_se), 2L)
  expect_true(all(is.na(res_ziqr$bootstrap_se)))
  expect_identical(res_ziqr$candidate_bootstrap_se[["zero_iqr_target"]], 0.0)
  expect_true(is.na(res_ziqr$crit_val))
  expect_identical(dim(res_ziqr$raw_draws), c(8L, 2L))

  # 2. Zero analytical variance (anti-correlated units within each cluster on target 2):
  u8 <- sprintf("u%d", 1:8)
  map8 <- list(unit_ids = u8, cluster_ids = rep(c4, each = 2L))
  if_zero_var <- matrix(
    c(
       2.0,  1.0,
       2.0, -1.0,
      -1.0,  3.0,
      -1.0, -3.0,
       0.5, -2.0,
       0.5,  2.0,
      -1.5,  0.5,
      -1.5, -0.5
    ),
    nrow = 8, ncol = 2, byrow = TRUE
  )
  draws16 <- sprintf("d%02d", 1:16)
  res_zvar <- fn(
    target_ids = c("t_pos", "t_anticorr"),
    expected_unit_ids = u8,
    unit_cluster_map = map8,
    cluster_ids = c4,
    n_denom = 8.0,
    center = TRUE,
    if_batch_fn = function(i) list(unit_ids = u8, target_ids = c("t_pos", "t_anticorr"), if_matrix = if_zero_var),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(draw_ids = draws16, cluster_ids = c4, sign_matrix = grid4),
    n_sign_batches = 1L,
    expected_draw_ids = draws16,
    n_draws = 16L,
    alpha = 0.05,
    return_raw_draws = TRUE,
    budget = make_explicit_budget()
  )
  expect_false(res_zvar$accepted)
  expect_true(res_zvar$refused)
  expect_identical(
    res_zvar$target_status,
    c(t_pos = "ok", t_anticorr = "zero_analytical_and_bootstrap_scale")
  )
  expect_true(is.finite(res_zvar$analytic_se[["t_pos"]]) && res_zvar$analytic_se[["t_pos"]] > 0)
  expect_true(is.na(res_zvar$analytic_se[["t_anticorr"]]))
  expect_true(all(is.na(res_zvar$bootstrap_se)))
  expect_true(is.na(res_zvar$crit_val))
  expect_identical(dim(res_zvar$cov_mat), c(2L, 2L))

  # 3. Non-finite analytical variance from squaring finite 1e160 cluster sums:
  if_overflow_cov <- matrix(
    c(
       1e160,  2.0,
      -1e160, -1.0,
       1e160,  0.5,
      -1e160, -1.5
    ),
    nrow = 4, ncol = 2, byrow = TRUE
  )
  res_inf_cov <- fn(
    target_ids = c("t_overflow", "t_normal"),
    expected_unit_ids = u4,
    unit_cluster_map = map4,
    cluster_ids = c4,
    n_denom = 1.0,
    center = FALSE,
    if_batch_fn = function(i) list(unit_ids = u4, target_ids = c("t_overflow", "t_normal"), if_matrix = if_overflow_cov),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(draw_ids = draws16, cluster_ids = c4, sign_matrix = grid4),
    n_sign_batches = 1L,
    expected_draw_ids = draws16,
    n_draws = 16L,
    alpha = 0.05,
    return_raw_draws = FALSE,
    budget = make_explicit_budget()
  )
  expect_false(res_inf_cov$accepted)
  expect_identical(res_inf_cov$diagnostics$analytic_status[["t_overflow"]], "nonfinite_analytical_variance")
  expect_true(is.na(res_inf_cov$crit_val))
})

test_that("duplicate, missing, unknown, and misaligned unit/cluster/target/draw inputs and invalid signs/IFs are refused", {
  fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  call_with <- function(overrides) {
    args <- list(
      target_ids = fx$target_ids,
      expected_unit_ids = fx$expected_unit_ids,
      unit_cluster_map = fx$unit_cluster_map,
      cluster_ids = fx$cluster_ids,
      n_denom = fx$n_denom,
      center = TRUE,
      if_batch_fn = function(i) fx$if_batches[[i]],
      n_if_batches = length(fx$if_batches),
      sign_batch_fn = function(i) fx$sign_batches[[i]],
      n_sign_batches = length(fx$sign_batches),
      expected_draw_ids = fx$draw_ids,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = FALSE,
      budget = fx$budget
    )
    for (nm in names(overrides)) {
      args[[nm]] <- overrides[[nm]]
    }
    do.call(fn, args)
  }

  # Duplicate declared IDs
  expect_error(call_with(list(target_ids = c("t1", "t1"))), class = "cdid_multiplier_duplicate_target_error")
  expect_error(call_with(list(expected_unit_ids = c("u1", "u1", "u3", "u4", "u5", "u6"))), class = "cdid_multiplier_duplicate_unit_error")
  expect_error(call_with(list(cluster_ids = c("c1", "c1", "c2", "c3"))), class = "cdid_multiplier_duplicate_cluster_error")
  expect_error(call_with(list(expected_draw_ids = c("d1", "d1", fx$draw_ids[3:16]))), class = "cdid_multiplier_duplicate_draw_error")

  # unit_cluster_map unknown/missing/duplicate unit or unknown cluster
  expect_error(
    call_with(list(unit_cluster_map = list(
      unit_ids = c("u_unknown", fx$expected_unit_ids[-1L]),
      cluster_ids = fx$unit_cluster_map$cluster_ids
    ))),
    class = "cdid_multiplier_unknown_unit_error"
  )
  expect_error(
    call_with(list(unit_cluster_map = list(
      unit_ids = fx$expected_unit_ids[-1L],
      cluster_ids = fx$unit_cluster_map$cluster_ids[-1L]
    ))),
    class = "cdid_multiplier_missing_unit_error"
  )
  expect_error(
    call_with(list(unit_cluster_map = list(
      unit_ids = fx$unit_cluster_map$unit_ids,
      cluster_ids = c("cl_unknown", fx$unit_cluster_map$cluster_ids[-1L])
    ))),
    class = "cdid_multiplier_unknown_cluster_error"
  )

  # IF batch: wrong target order, unknown target, duplicate unit within batch, duplicate unit across batches, unknown unit, missing unit
  expect_error(
    call_with(list(if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      if (i == 1L) b$target_ids <- rev(b$target_ids)
      b
    })),
    class = "cdid_multiplier_target_order_error"
  )
  expect_error(
    call_with(list(if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      if (i == 1L) b$target_ids <- c("att_g2_t2", "unknown_target")
      b
    })),
    class = "cdid_multiplier_unknown_target_error"
  )
  expect_error(
    call_with(list(if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      if (i == 1L) b$unit_ids <- c("u1", "u1")
      b
    })),
    class = "cdid_multiplier_duplicate_unit_error"
  )
  expect_error(
    call_with(list(if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      if (i == 2L) b$unit_ids[1L] <- "u1" # already seen in batch 1
      b
    })),
    class = "cdid_multiplier_duplicate_unit_error"
  )
  expect_error(
    call_with(list(if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      if (i == 3L) b$unit_ids <- "u_late_unknown"
      b
    })),
    class = "cdid_multiplier_unknown_unit_error"
  )
  expect_error(
    call_with(list(n_if_batches = 2L)), # omits batch 3 (missing unit u4)
    class = "cdid_multiplier_missing_unit_error"
  )
  expect_error(
    call_with(list(if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      if (i == 1L) b$cluster_ids <- c("cl_a", "cl_a") # disagrees with unit_cluster_map
      b
    })),
    class = "cdid_multiplier_cluster_alignment_error"
  )
  expect_error(
    call_with(list(if_batch_fn = function(i) {
      b <- fx$if_batches[[i]]
      if (i == 2L) b$if_matrix[1L, 1L] <- NA_real_
      b
    })),
    class = "cdid_multiplier_nonfinite_if_error"
  )

  # Sign batch: incomplete clusters, duplicate cluster, unknown cluster, wrong cluster order,
  # duplicate draw, unknown draw, wrong draw order, incomplete draws, invalid sign value
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 1L) {
        sb$cluster_ids <- sb$cluster_ids[1:3]
        sb$sign_matrix <- sb$sign_matrix[, 1:3, drop = FALSE]
      }
      sb
    })),
    class = "cdid_multiplier_incomplete_cluster_error"
  )
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 1L) sb$cluster_ids <- c("cl_d", "cl_d", "cl_c", "cl_b")
      sb
    })),
    class = "cdid_multiplier_duplicate_cluster_error"
  )
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 1L) sb$cluster_ids <- c("cl_d", "cl_a", "cl_c", "cl_unknown")
      sb
    })),
    class = "cdid_multiplier_unknown_cluster_error"
  )
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 1L) sb$cluster_ids <- rev(sb$cluster_ids)
      sb
    })),
    class = "cdid_multiplier_cluster_order_error"
  )
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 1L) sb$draw_ids[2L] <- sb$draw_ids[1L]
      sb
    })),
    class = "cdid_multiplier_duplicate_draw_error"
  )
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 1L) sb$draw_ids[1L] <- "draw_unknown"
      sb
    })),
    class = "cdid_multiplier_unknown_draw_error"
  )
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 1L) sb$draw_ids <- rev(sb$draw_ids)
      sb
    })),
    class = "cdid_multiplier_draw_order_error"
  )
  expect_error(
    call_with(list(n_sign_batches = 1L)), # only 6 of 16 draws supplied
    class = "cdid_multiplier_incomplete_draws_error"
  )
  expect_error(
    call_with(list(sign_batch_fn = function(i) {
      sb <- fx$sign_batches[[i]]
      if (i == 2L) sb$sign_matrix[1L, 1L] <- 0.0
      sb
    })),
    class = "cdid_multiplier_invalid_sign_error"
  )
})

test_that("budget preflight, internal B x J budgeting, >INT_MAX cell budgets, and pre-scan alignment order hold", {
  fn <- get_private_fn("chained_cluster_multiplier")
  est_fn <- get_private_fn("compute_multiplier_workspace_estimate")
  fx <- build_fixture_c4()

  ws <- est_fn(n_units = 6L, n_clusters = 4L, n_targets = 2L, n_draws = 16L,
               max_if_batch_rows = 64, max_sign_batch_rows = 128)
  # Independently count 2CJ + 6J^2 + 2BJ + 6n + 4C + 8B + 32J
  # plus 2*min(n,64)*J + 4*min(B,128)*C: 628 payload slots.
  expect_identical(ws$total_persistent_cells, 628.0)

  # Top-level budget failures must happen BEFORE invoking any callback:
  if_calls <- 0L
  sign_calls <- 0L
  expect_error(
    fn(
      target_ids = fx$target_ids,
      expected_unit_ids = fx$expected_unit_ids,
      unit_cluster_map = fx$unit_cluster_map,
      cluster_ids = fx$cluster_ids,
      n_denom = fx$n_denom,
      center = TRUE,
      if_batch_fn = function(i) {
        if_calls <<- if_calls + 1L
        fx$if_batches[[i]]
      },
      n_if_batches = length(fx$if_batches),
      sign_batch_fn = function(i) {
        sign_calls <<- sign_calls + 1L
        fx$sign_batches[[i]]
      },
      n_sign_batches = length(fx$sign_batches),
      expected_draw_ids = fx$draw_ids,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = FALSE, # internal B x J is STILL budgeted!
      budget = make_explicit_budget(max_workspace_cells = 627)
    ),
    class = "cdid_multiplier_budget_error"
  )
  expect_identical(if_calls, 0L)
  expect_identical(sign_calls, 0L)

  # Exactly 628 workspace slots suffice; scalar budgets above INT_MAX are allowed.
  res_ok <- fn(
    target_ids = fx$target_ids,
    expected_unit_ids = fx$expected_unit_ids,
    unit_cluster_map = fx$unit_cluster_map,
    cluster_ids = fx$cluster_ids,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) fx$if_batches[[i]],
    n_if_batches = length(fx$if_batches),
    sign_batch_fn = function(i) fx$sign_batches[[i]],
    n_sign_batches = length(fx$sign_batches),
    expected_draw_ids = fx$draw_ids,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = FALSE,
    budget = make_explicit_budget(
      max_workspace_cells = 628,
      max_workspace_bytes = 5e11 # > .Machine$integer.max, no invented INT_MAX cap
    )
  )
  expect_true(res_ok$accepted)
  expect_null(res_ok$raw_draws)

  # Pre-scan validation order: a sign batch with BOTH wrong cluster order AND invalid sign 0
  # must fail on cluster order (before scanning sign_matrix content):
  expect_error(
    fn(
      target_ids = fx$target_ids,
      expected_unit_ids = fx$expected_unit_ids,
      unit_cluster_map = fx$unit_cluster_map,
      cluster_ids = fx$cluster_ids,
      n_denom = fx$n_denom,
      center = TRUE,
      if_batch_fn = function(i) fx$if_batches[[i]],
      n_if_batches = length(fx$if_batches),
      sign_batch_fn = function(i) {
        sb <- fx$sign_batches[[i]]
        sb$cluster_ids <- rev(sb$cluster_ids)
        sb$sign_matrix[1L, 1L] <- 0.0
        sb
      },
      n_sign_batches = length(fx$sign_batches),
      expected_draw_ids = fx$draw_ids,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = FALSE,
      budget = fx$budget
    ),
    class = "cdid_multiplier_cluster_order_error"
  )

  # A sign batch exceeding max_sign_batch_rows AND containing invalid sign 0
  # must fail with cdid_multiplier_budget_error (before scanning sign_matrix content):
  expect_error(
    fn(
      target_ids = fx$target_ids,
      expected_unit_ids = fx$expected_unit_ids,
      unit_cluster_map = fx$unit_cluster_map,
      cluster_ids = fx$cluster_ids,
      n_denom = fx$n_denom,
      center = TRUE,
      if_batch_fn = function(i) fx$if_batches[[i]],
      n_if_batches = length(fx$if_batches),
      sign_batch_fn = function(i) {
        sb <- fx$sign_batches[[i]]
        sb$sign_matrix[1L, 1L] <- 0.0
        sb
      },
      n_sign_batches = length(fx$sign_batches),
      expected_draw_ids = fx$draw_ids,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = FALSE,
      budget = make_explicit_budget(max_sign_batch_rows = 5)
    ),
    class = "cdid_multiplier_budget_error"
  )
})

test_that("integer IDs, strict ID type matching without coercion, ambient state isolation, and digest reference pins hold", {
  fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  # Ambient state snapshot before execution
  wd_before <- getwd()
  opts_before <- options()
  has_seed_before <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  seed_before <- if (has_seed_before) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL

  # Integer IDs end-to-end
  int_targets <- c(10L, 20L)
  int_clusters <- c(400L, 100L, 300L, 200L)
  int_units <- c(1L, 2L, 3L, 4L, 5L, 6L)
  int_draws <- seq_len(16L)
  int_map <- list(
    unit_ids = c(6L, 2L, 4L, 1L, 5L, 3L),
    cluster_ids = c(300L, 100L, 400L, 400L, 300L, 300L)
  )

  res_int <- fn(
    target_ids = int_targets,
    expected_unit_ids = int_units,
    unit_cluster_map = int_map,
    cluster_ids = int_clusters,
    n_denom = fx$n_denom,
    center = TRUE,
    if_batch_fn = function(i) list(
      unit_ids = int_units,
      target_ids = int_targets,
      if_matrix = unname(fx$if_rows)
    ),
    n_if_batches = 1L,
    sign_batch_fn = function(i) list(
      draw_ids = int_draws,
      cluster_ids = int_clusters,
      sign_matrix = fx$signs
    ),
    n_sign_batches = 1L,
    expected_draw_ids = int_draws,
    n_draws = 16L,
    alpha = fx$alpha,
    return_raw_draws = TRUE,
    budget = fx$budget
  )

  expect_true(res_int$accepted)
  expect_identical(res_int$metadata$caller_premises_only, TRUE)
  expect_identical(res_int$metadata$sampling_provenance_certified, FALSE)
  expect_identical(res_int$metadata$complete_if_or_nuisance_error_certified, FALSE)
  expect_identical(res_int$metadata$finite_sample_or_asymptotic_validity_certified, FALSE)
  expect_identical(res_int$metadata$entry_budget_is_rss_or_cpu_bound, FALSE)
  expect_identical(res_int$metadata$candidate_floating_arithmetic_only, TRUE)

  # Digest check of exact dyadic quantities
  expect_true(requireNamespace("digest", quietly = TRUE))
  {
    dyadic_payload <- list(
      cluster_sums = unname(res_int$cluster_sums),
      cluster_means = unname(res_int$cluster_means),
      cluster_sums_centered = unname(res_int$cluster_sums_centered),
      cov_mat = unname(res_int$cov_mat),
      omega_cl_hat = unname(res_int$omega_cl_hat),
      q25 = unname(res_int$q25),
      q75 = unname(res_int$q75),
      raw_bootstrap_iqr = unname(res_int$raw_bootstrap_iqr),
      raw_draws = unname(res_int$raw_draws)
    )
    expect_identical(
      digest::digest(dyadic_payload, algo = "sha256"),
      digest::digest(
        list(
          cluster_sums = matrix(c(5, -1, 2, 0, 4, 8, -2, 0), nrow = 4, ncol = 2),
          cluster_means = c(1.5, 2.5),
          cluster_sums_centered = matrix(c(3.5, -2.5, 0.5, -1.5, 1.5, 5.5, -4.5, -2.5), nrow = 4, ncol = 2),
          cov_mat = matrix(c(21 / 64, -7 / 64, -7 / 64, 59 / 64), nrow = 2, ncol = 2),
          omega_cl_hat = matrix(c(21 / 8, -7 / 8, -7 / 8, 59 / 8), nrow = 2, ncol = 2),
          q25 = c(-2.0, -3.0),
          q75 = c(1.5, 2.5),
          raw_bootstrap_iqr = c(3.5, 5.5),
          raw_draws = (fx$signs %*% matrix(c(3.5, -2.5, 0.5, -1.5, 1.5, 5.5, -4.5, -2.5), nrow = 4, ncol = 2)) / 2.0
        ),
        algo = "sha256"
      )
    )
  }

  # Type mismatch between character and integer IDs must be refused without coercion:
  expect_error(
    fn(
      target_ids = int_targets,
      expected_unit_ids = int_units,
      unit_cluster_map = list(
        unit_ids = as.character(int_map$unit_ids),
        cluster_ids = int_map$cluster_ids
      ),
      cluster_ids = int_clusters,
      n_denom = fx$n_denom,
      center = TRUE,
      if_batch_fn = function(i) list(unit_ids = int_units, target_ids = int_targets, if_matrix = unname(fx$if_rows)),
      n_if_batches = 1L,
      sign_batch_fn = function(i) list(draw_ids = int_draws, cluster_ids = int_clusters, sign_matrix = fx$signs),
      n_sign_batches = 1L,
      expected_draw_ids = int_draws,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = FALSE,
      budget = fx$budget
    ),
    class = "cdid_multiplier_unit_alignment_error"
  )

  # Verify ambient state was not touched
  expect_identical(getwd(), wd_before)
  expect_identical(options(), opts_before)
  expect_identical(exists(".Random.seed", envir = globalenv(), inherits = FALSE), has_seed_before)
  if (has_seed_before) {
    expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE), seed_before)
  }
})

test_that("actual source mutants of chained_cluster_multiplier closure are actively detected and killed", {
  base_fn <- get_private_fn("chained_cluster_multiplier")
  fx <- build_fixture_c4()

  run_on_fn <- function(fn_under_test, scale_if = 1.0) {
    fn_under_test(
      target_ids = fx$target_ids,
      expected_unit_ids = fx$expected_unit_ids,
      unit_cluster_map = fx$unit_cluster_map,
      cluster_ids = fx$cluster_ids,
      n_denom = fx$n_denom,
      center = TRUE,
      if_batch_fn = function(i) {
        b <- fx$if_batches[[i]]
        b$if_matrix <- b$if_matrix * scale_if
        b
      },
      n_if_batches = length(fx$if_batches),
      sign_batch_fn = function(i) fx$sign_batches[[i]],
      n_sign_batches = length(fx$sign_batches),
      expected_draw_ids = fx$draw_ids,
      n_draws = 16L,
      alpha = fx$alpha,
      return_raw_draws = TRUE,
      budget = fx$budget
    )
  }

  ref_out <- run_on_fn(base_fn)

  mutate_closure <- function(pattern, replacement) {
    src <- deparse(base_fn, width.cutoff = 500L)
    hits <- grep(pattern, src, fixed = TRUE)
    expect_true(length(hits) >= 1L, info = paste("Pattern not found in closure:", pattern))
    src[hits[1L]] <- sub(pattern, replacement, src[hits[1L]], fixed = TRUE)
    mut_env <- new.env(parent = environment(base_fn))
    eval(parse(text = paste(src, collapse = "\n")), envir = mut_env)
  }

  mutants <- list(
    M01_ignore_unit_id_order = list(
      pattern = "u_pos <- match(b_uids, expected_unit_ids)",
      replacement = "u_pos <- seq.int(total_if_rows_seen + 1L, total_if_rows_seen + r_b)"
    ),
    M02_sort_clusters_alphabetically = list(
      pattern = "unit_cluster_idx <- match(ordered_unit_clusters, cluster_ids)",
      replacement = "unit_cluster_idx <- match(ordered_unit_clusters, sort(cluster_ids))"
    ),
    M03_within_cluster_mean_instead_of_sum = list(
      pattern = "cluster_means <- colMeans(S_mat)",
      replacement = "S_mat <- S_mat / pmax(1L, streamed_cluster_unit_counts); cluster_means <- colMeans(S_mat)"
    ),
    M04_unit_level_centering = list(
      pattern = "S_centered <- sweep(S_mat, 2L, cluster_means, FUN = \"-\")",
      replacement = "S_centered <- S_mat - outer(as.double(streamed_cluster_unit_counts), colSums(S_mat) / as.double(n_units))"
    ),
    M05_center_by_nonempty_clusters_Cplus = list(
      pattern = "cluster_means <- colMeans(S_mat)",
      replacement = "cluster_means <- colSums(S_mat) / sum(streamed_cluster_unit_counts > 0L)"
    ),
    M06_cov_denom_n_instead_of_n2 = list(
      pattern = "cov_mat <- cluster_crossprod/(n_denom_num * n_denom_num)",
      replacement = "cov_mat <- cluster_crossprod/n_denom_num"
    ),
    M07_cov_denom_C2_instead_of_n2 = list(
      pattern = "cov_mat <- cluster_crossprod/(n_denom_num * n_denom_num)",
      replacement = "cov_mat <- cluster_crossprod/(as.double(n_clusters) * as.double(n_clusters))"
    ),
    M08_cov_denom_n_units_instead_of_n_denom = list(
      pattern = "cov_mat <- cluster_crossprod/(n_denom_num * n_denom_num)",
      replacement = "cov_mat <- cluster_crossprod/(as.double(n_units) * as.double(n_units))"
    ),
    M09_bres_scale_sqrt_n_instead_of_sqrt_C = list(
      pattern = "inv_sqrt_C <- 1/sqrt(as.double(n_clusters))",
      replacement = "inv_sqrt_C <- 1/sqrt(n_denom_num)"
    ),
    M10_se_scale_1_over_sqrt_n = list(
      pattern = "se_scale_factor <- sqrt(as.double(n_clusters))/n_denom_num",
      replacement = "se_scale_factor <- 1/sqrt(n_denom_num)"
    ),
    M11_se_scale_1_over_sqrt_C = list(
      pattern = "se_scale_factor <- sqrt(as.double(n_clusters))/n_denom_num",
      replacement = "se_scale_factor <- 1/sqrt(as.double(n_clusters))"
    ),
    M12_type7_iqr_quantile = list(
      pattern = "iqr_j <- q75_j - q25_j",
      replacement = "iqr_j <- unname(diff(stats::quantile(col_j, probs = c(0.25, 0.75), type = 7)))"
    ),
    M13_type7_critical_quantile = list(
      pattern = "crit_val <- unname(sorted_max[k_crit])",
      replacement = "crit_val <- unname(stats::quantile(max_stat, probs = 1 - alpha_num, type = 7))"
    ),
    M14_double_weight_by_cluster_count = list(
      pattern = "S_mat[c_idx, ] <- S_mat[c_idx, ] + if_mat[r, ]",
      replacement = "S_mat[c_idx, ] <- S_mat[c_idx, ] + if_mat[r, ] * as.double(declared_cluster_unit_counts[c_idx])"
    )
  )

  for (mname in names(mutants)) {
    spec <- mutants[[mname]]
    mut_fn <- mutate_closure(spec$pattern, spec$replacement)
    mut_out <- run_on_fn(mut_fn)
    detected <- (
      !isTRUE(all.equal(mut_out$cluster_sums, ref_out$cluster_sums, tolerance = 1e-12)) ||
      !isTRUE(all.equal(mut_out$cluster_sums_centered, ref_out$cluster_sums_centered, tolerance = 1e-12)) ||
      !isTRUE(all.equal(mut_out$cov_mat, ref_out$cov_mat, tolerance = 1e-12)) ||
      !isTRUE(all.equal(mut_out$raw_draws, ref_out$raw_draws, tolerance = 1e-12)) ||
      !isTRUE(all.equal(mut_out$bootstrap_se, ref_out$bootstrap_se, tolerance = 1e-12)) ||
      !isTRUE(all.equal(mut_out$crit_val, ref_out$crit_val, tolerance = 1e-12))
    )
    expect_true(detected, info = paste("Mutant survived:", mname))
  }

  # M15: upstream absolute magnitude screen `colSums(bres^2) <= sqrt(.Machine$double.eps) * 10`
  mut_m15 <- mutate_closure(
    "col_j <- bres[, j]",
    "col_j <- if (sum(bres[, j]^2) <= sqrt(.Machine$double.eps) * 10) rep(0.0, nrow(bres)) else bres[, j]"
  )
  out_m15_tiny <- run_on_fn(mut_m15, scale_if = 1e-8)
  expect_false(identical(out_m15_tiny$accepted, TRUE), info = "M15 upstream screen mutant was not detected on 1e-8 scale")
})

