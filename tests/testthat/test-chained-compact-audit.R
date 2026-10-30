# Independent Cartesian self-join with exact numeric IDs and the full known
# calendar. Small bounded fixtures only; production traversal stores no chunks.
audit_cartesian_pairs <- function(d) {
  observed <- d[!is.na(d$y), c("id", "time", "y", "g")]
  z <- merge(observed, observed, by = "id", suffixes = c("_s", "_t"))
  z <- z[z$time_s < z$time_t, , drop = FALSE]
  ids <- sort(unique(d$id), method = "radix")
  ui <- match(z$id, ids)
  ord <- order(ui, z$time_s, z$time_t, method = "radix")
  z <- z[ord, , drop = FALSE]; ui <- ui[ord]
  io <- vapply(seq_len(nrow(z)), function(j) sum(observed$id == z$id[j] &
    observed$time > z$time_s[j] & observed$time < z$time_t[j]), integer(1L))
  calendar <- sort(unique(d$time))
  ip <- vapply(seq_len(nrow(z)), function(j) sum(calendar > z$time_s[j] & calendar < z$time_t[j]), integer(1L))
  gap <- as.double(z$time_t) - as.double(z$time_s)
  data.frame(pair_index = as.double(seq_len(nrow(z))), unit_index = ui, id = z$id,
    # The constructor contract uses double dates/cohorts/outcomes; converting
    # R integers is exact and changes neither join identity nor any predicate.
    g = as.double(z$g_s), s = as.double(z$time_s), t = as.double(z$time_t),
    y_s = as.double(z$y_s), y_t = as.double(z$y_t),
    dy = as.double(z$y_t) - as.double(z$y_s), time_gap = gap,
    time_gap_finite = is.finite(gap),
    n_intervening_obs = io, n_intervening_panel_dates = ip,
    is_consecutive_unit = io == 0L, is_consecutive_panel = ip == 0L)
}

audit_collect <- function(index, size) {
  # Explicitly bounded test accumulator; not used in production traversal.
  stopifnot(index$counts$n_pairs <= 100000)
  chunks <- list()
  result <- cdid:::traverse_observed_pairs(index, size, function(b) chunks[[length(chunks) + 1L]] <<- b)
  out <- if (length(chunks)) do.call(rbind, chunks) else cdid:::empty_pair_batch(index)
  rownames(out) <- NULL
  list(pairs = out, traversal = result)
}

audit_svd_estimable <- function(W, c) {
  v <- svd(W, nu = 0L, nv = ncol(W))
  rank <- sum(v$d > 1e-12) # Reference tolerance only; graph code uses none.
  if (!rank) return(all(c == 0))
  basis <- v$v[, seq_len(rank), drop = FALSE]
  max(abs(c - basis %*% crossprod(basis, c))) < 1e-12
}

audit_dense_certificate <- function(d, cert) {
  out <- numeric(d$incidence$n_links)
  if (!is.null(cert)) out[cert$link_index] <- cert$sign
  out
}

test_that("bounded traversal exactly matches a Cartesian join at sizes 1, 2 and above P", {
  set.seed(8102)
  sparse <- expand.grid(id = c(1, 1 + 2^-40, 3, 4, 5, 6), time = c(-2, 0, 1, 1 + 2^-40, 5))
  sparse$g <- ifelse(sparse$id <= 3, 2, Inf)
  sparse$y <- rnorm(nrow(sparse)); sparse$y[sample(nrow(sparse), 12)] <- NA_real_
  sparse$y[sparse$id == 5] <- NA_real_
  sparse$y[sparse$id == 6 & sparse$time != 0] <- NA_real_
  long <- data.frame(id = 1, time = seq_len(31), y = seq_len(31)^2, g = 10)
  gap <- data.frame(id = c(1, 1, 1, 2), time = c(1, 2, 3, 2), y = c(0, NA, 1, NA), g = c(3, 3, 3, 0))
  for (dat in list(sparse, long, gap)) {
    snapshot <- dat
    index <- cdid:::construct_observed_pairs(dat, "id", "time", "y", "g")
    oracle <- audit_cartesian_pairs(dat)
    expect_identical(dat, snapshot)
    expect_null(index$pairs)
    expect_null(index$unit_obs_times)
    for (size in c(1, 2, nrow(oracle) + 5)) {
      got <- audit_collect(index, size)
      expect_identical(got$pairs, oracle)
      expect_identical(got$traversal$n_pairs, as.double(nrow(oracle)))
      expect_lte(got$traversal$max_batch_size, size)
      expect_identical(cdid:::materialize_observed_pairs(index, max_pairs = nrow(oracle), chunk_size = size)$pairs, oracle)
    }
    for (j in 1:10) {
      shuffled <- cdid:::construct_observed_pairs(dat[sample(nrow(dat)), ], "id", "time", "y", "g")
      expect_identical(shuffled, index)
    }
  }
})

test_that("scalar selected columns, ambiguous names and invalid bounds are rejected", {
  d <- data.frame(id = c(1, 1), time = c(1, 2), y = c(2, 3), g = c(2, 2))
  for (nm in names(d)) {
    for (wide in c(FALSE, TRUE)) {
      bad <- d; bad[[nm]] <- if (wide) cbind(d[[nm]], d[[nm]]) else matrix(d[[nm]], ncol = 1)
      expect_error(cdid:::construct_observed_pairs(bad, "id", "time", "y", "g"), class = "cdid_error")
    }
    bad <- d; bad[[paste0(nm, "2")]] <- bad[[nm]]; names(bad)[5] <- nm
    expect_error(cdid:::construct_observed_pairs(bad, "id", "time", "y", "g"), class = "cdid_pairs_column_error")
  }
  index <- cdid:::construct_observed_pairs(d, "id", "time", "y", "g")
  for (size in list(0, -1, 1.5, Inf, NA_real_, matrix(2), 2^54)) {
    expect_error(cdid:::traverse_observed_pairs(index, size, identity), class = "cdid_pairs_bound_error")
  }
  expect_error(cdid:::materialize_observed_pairs(index, max_pairs = 0, chunk_size = 1), class = "cdid_pairs_bound_error")
  expect_error(cdid:::materialize_observed_pairs(index, 2^54, 1), class = "cdid_pairs_bound_error")
  expect_error(cdid:::materialize_observed_pairs(index, chunk_size = 1)) # no implicit bound
  marker <- structure(list(message = "callback marker", call = NULL), class = c("audit_callback_error", "error", "condition"))
  err <- tryCatch(cdid:::traverse_observed_pairs(index, 1, function(b) stop(marker)), error = identity)
  expect_identical(err, marker)
  empty <- cdid:::construct_observed_pairs(d[FALSE, ], "id", "time", "y", "g")
  expect_identical(cdid:::materialize_observed_pairs(empty, 0, 1)$pairs, audit_cartesian_pairs(d[FALSE, ]))
  expect_identical(cdid:::empty_pair_batch(empty, 0L), audit_cartesian_pairs(d[FALSE, ]))
  nonempty_batch <- audit_collect(index, 1)$pairs
  expect_identical(cdid:::empty_pair_batch(index, 0L), nonempty_batch[FALSE, ])
  d_chr <- d; d_chr$id <- c("u1", "u1")
  idx_chr <- cdid:::construct_observed_pairs(d_chr, "id", "time", "y", "g")
  expect_identical(cdid:::empty_pair_batch(idx_chr, 0L), audit_collect(idx_chr, 1)$pairs[FALSE, ])
  expect_identical(empty$units$singleton_time, numeric())
  expect_identical(empty$units$n_pairs, numeric())
  calls <- 0L
  cdid:::traverse_observed_pairs(empty, 1, function(b) calls <<- calls + 1L)
  expect_identical(calls, 0L)
})

test_that("pair counts avoid integer overflow, extreme finite dates preserve finite dy with honest gap metadata, and outcome overflow precedes callbacks", {
  expected <- 65536 * 65535 / 2
  expect_identical(cdid:::pair_count_safe(65536L)$total, expected)
  expect_identical(cdid:::pair_count_safe(c(65536L, 65537L))$total, 4294967296)
  expect_error(cdid:::pair_count_safe(2^28), class = "cdid_pairs_count_error")
  expect_error(cdid:::pair_count_safe(rep(2^26, 5)), class = "cdid_pairs_count_error")
  d <- data.frame(id = 1, time = 1:65537, y = 0, g = Inf)
  index <- cdid:::construct_observed_pairs(d, "id", "time", "y", "g")
  expect_identical(index$counts$n_pairs, 2147516416)
  expect_lt(as.double(object.size(index)), 5e6)
  expect_error(cdid:::materialize_observed_pairs(index, 100000, 2), class = "cdid_pairs_bound_error")
  # Finite numeric dates -1e308, +1e308 with finite outcome differences are retained;
  # overflowing elapsed time gap is recorded honestly as Inf with time_gap_finite = FALSE.
  d_ext <- data.frame(
    id = c(1, 1, 1, 2, 3),
    time = c(-1e308, 0, 1e308, -1e308, 0),
    y = c(0, 0.25, 1, 5, NA_real_),
    g = c(0, 0, 0, Inf, Inf)
  )
  idx_ext <- cdid:::construct_observed_pairs(d_ext, "id", "time", "y", "g")
  oracle_ext <- audit_cartesian_pairs(d_ext)
  for (sz in c(1, 2, 5)) {
    got_ext <- audit_collect(idx_ext, sz)
    expect_identical(got_ext$pairs, oracle_ext)
    expect_identical(cdid:::materialize_observed_pairs(idx_ext, 10, sz)$pairs, oracle_ext)
  }
  expect_identical(oracle_ext$dy, c(0.25, 1.0, 0.75))
  expect_identical(oracle_ext$time_gap, c(1e308, Inf, 1e308))
  expect_identical(oracle_ext$time_gap_finite, c(TRUE, FALSE, TRUE))
  expect_identical(cdid:::empty_pair_batch(idx_ext, 0L), oracle_ext[FALSE, ])
  # Outcome-difference overflow is still rejected before any callback can run,
  # even when a preceding unit has valid outcomes and extreme finite dates.
  d_bad_y <- rbind(
    data.frame(id = 1, time = c(-1e308, 1e308), y = c(0, 1), g = 0),
    data.frame(id = 2, time = c(1, 2, 3), y = c(0, -.Machine$double.xmax, .Machine$double.xmax), g = 2)
  )
  expect_error(cdid:::construct_observed_pairs(d_bad_y, "id", "time", "y", "g"), class = "cdid_pairs_outcome_error")
  d_bad_ext_y <- data.frame(id = 1, time = c(-1e308, 1e308), y = c(-.Machine$double.xmax, .Machine$double.xmax), g = 0)
  expect_error(cdid:::construct_observed_pairs(d_bad_ext_y, "id", "time", "y", "g"), class = "cdid_pairs_outcome_error")
})

test_that("numeric dates, cohorts and IDs remain distinct despite identical display rounding", {
  a <- 2; b <- 2 + 2^-40
  expect_identical(format(a, digits = 7), format(b, digits = 7))
  d <- data.frame(id = c(a, a, b, b), time = c(1, a, 1, b), y = 1:4, g = c(a, a, b, b))
  i <- cdid:::construct_observed_pairs(d, "id", "time", "y", "g")
  expect_identical(i$units$id, c(a, b)); expect_identical(i$units$g, c(a, b))
  expect_identical(audit_collect(i, 1)$pairs, audit_cartesian_pairs(d))
  links <- data.frame(g = c(a, b), s = 1, t = c(a, b), ref = 1, functional_key = "opaque")
  des <- cdid:::build_chained_design(links, target_post_times = c(b, b + 1))
  expect_identical(des$parameters$g, c(a, a, a, b, b))
  expect_identical(des$rank_diagnostics$rank, 2L)
  same_cohort <- links; same_cohort$g <- a
  got <- cdid:::build_chained_design(same_cohort)
  expect_identical(got$parameters$time, c(a, b))
  expect_identical(unname(cdid:::materialize_chained_design(got, 100)$W), diag(1L, 2))
})

test_that("thin incidence, lazy certificates and graph estimability match hand matrices and SVD", {
  links <- data.frame(g = c(3, 3, 4, 4, 5), s = c(2, 4, 1, 1, 4), t = c(3, 6, 3, 4, 5),
    ref = c(2, 2, 3, 3, 4), functional_key = c("a", "b", "A", "B", "c"))
  d <- cdid:::build_chained_design(links, data.frame(g = c(3, 3, 3, 3, 4, 5), time = c(3, 4, 5, 6, 4, 5)))
  # Column order: ATT33, ATT34, ATT35, ATT36, pre4A, pre4B, ATT44, ATT55.
  W <- rbind(c(1,0,0,0,0,0,0,0), c(0,-1,0,1,0,0,0,0),
    c(0,0,0,0,-1,0,0,0), c(0,0,0,0,0,-1,1,0), c(0,0,0,0,0,0,0,1))
  dense <- cdid:::materialize_chained_design(d, 1000)
  expect_identical(unname(dense$W), matrix(as.integer(W), nrow = 5))
  expect_null(d$W); expect_null(d$rank_diagnostics$root_paths)
  expect_null(d$rank_diagnostics$nullspace_basis); expect_null(d$rank_diagnostics$path_witnesses)
  expect_identical(d$rank_diagnostics$rank, as.integer(qr(W)$rank))
  expect_identical(d$rank_diagnostics$rank, as.integer(sum(svd(W)$d > 1e-12)))
  expect_identical(d$rank_diagnostics$nullity, ncol(W) - d$rank_diagnostics$rank)
  X <- matrix(seq_len(ncol(W) * 3), nrow = ncol(W))
  M <- matrix(seq_len(nrow(W) * 2), nrow = nrow(W))
  expect_equal(unname(cdid:::incidence_multiply(d, X)), W %*% X, tolerance = 0)
  expect_equal(unname(cdid:::incidence_multiply(d, M, TRUE)), crossprod(W, M), tolerance = 0)
  expect_equal(unname(cdid:::incidence_multiply(d, X[, 1])), drop(W %*% X[, 1]), tolerance = 0)
  expect_equal(unname(cdid:::incidence_multiply(d, M[, 1], TRUE)), drop(crossprod(W, M[, 1])), tolerance = 0)
  B <- dense$rank_diagnostics$nullspace_basis
  expect_true(all(W %*% B == 0)); expect_identical(qr(B)$rank, d$rank_diagnostics$nullity)
  for (a in seq_len(ncol(W))) {
    e <- numeric(ncol(W)); e[a] <- 1
    level <- cdid:::check_level_estimable(d, a)
    expect_identical(level$estimable, audit_svd_estimable(W, e))
    if (level$estimable) expect_equal(drop(crossprod(audit_dense_certificate(d, level$witness_links), W)), e, tolerance = 0)
    for (b in seq_len(ncol(W))) {
      c <- numeric(ncol(W)); c[b] <- c[b] + 1; c[a] <- c[a] - 1
      contrast <- cdid:::check_contrast_estimable(d, a, b)
      expect_identical(contrast$estimable, audit_svd_estimable(W, c))
      if (contrast$estimable) expect_equal(drop(crossprod(audit_dense_certificate(d, contrast$witness_links), W)), c, tolerance = 0)
    }
  }
  expect_true(cdid:::check_contrast_estimable(d, 1, 8)$estimable) # separately anchored cohorts
  expect_true(cdid:::check_contrast_estimable(d, 2, 4)$estimable) # same floating component
  expect_false(cdid:::check_contrast_estimable(d, 2, 7)$estimable) # distinct floating components
  for (x in list(1.5, 0, -1, Inf, NA_real_, matrix(1), 2^40)) {
    expect_error(cdid:::check_contrast_estimable(d, x, 1), class = "cdid_design_input_error")
  }
  for (nm in c("g", "s", "t", "ref", "functional_key")) {
    bad <- links; bad[[nm]] <- matrix(bad[[nm]], ncol = 1L)
    expect_error(cdid:::build_chained_design(bad), class = "cdid_error")
    bad <- links; bad$duplicate <- bad[[nm]]; names(bad)[ncol(bad)] <- nm
    expect_error(cdid:::build_chained_design(bad), class = "cdid_design_input_error")
  }
  for (j in 1:10) {
    perm <- sample(nrow(links)); pd <- cdid:::build_chained_design(links[perm, ], d$post_targets[, c("g", "time")])
    expect_identical(pd$parameters, d$parameters); expect_identical(pd$vertices, d$vertices)
    expect_identical(pd$post_targets, d$post_targets)
    expect_equal(unname(cdid:::materialize_chained_design(pd, 1000)$W), W[perm, ], tolerance = 0)
  }
  for (targets in list(matrix(c(3, 4), ncol = 1),
                       data.frame(g = I(matrix(3)), time = 3))) {
    expect_error(cdid:::build_chained_design(links, targets), class = "cdid_error")
  }
  bad_target <- data.frame(g = 3, time = 3, duplicate = 4)
  names(bad_target)[3] <- "time"
  expect_error(cdid:::build_chained_design(links, bad_target), class = "cdid_design_input_error")
  for (nm in c("horizon", "post_shareable", "composition_corrected")) {
    bad <- links
    bad[[nm]] <- matrix(if (nm == "horizon") rep(6, nrow(links)) else rep(TRUE, nrow(links)), ncol = 1)
    expect_error(cdid:::build_chained_design(bad), class = "cdid_error")
  }
  # Finite operands that overflow a product are errors, never a fallback.
  edge <- cdid:::build_chained_design(data.frame(g = 1, s = 1, t = 2, ref = 0, functional_key = "A"))
  expect_error(cdid:::incidence_multiply(edge, c(-.Machine$double.xmax, .Machine$double.xmax)), class = "cdid_design_product_error")
  fan <- cdid:::build_chained_design(data.frame(g = 1, s = c(0,1), t = c(2,2), ref = 0, functional_key = "A"))
  expect_error(cdid:::incidence_multiply(fan, rep(.Machine$double.xmax, 2), TRUE), class = "cdid_design_product_error")
  expect_error(cdid:::materialize_chained_design(d, 1), class = "cdid_design_bound_error")
  expect_error(cdid:::materialize_chained_design(d))
  bad <- data.frame(g = 2, s = 1, t = 2, ref = 1, functional_key = "a", post_shareable = FALSE)
  expect_error(cdid:::build_chained_design(bad), class = "cdid_design_post_share_error")
})

test_that("metadata never creates or certifies population-operator keys", {
  l <- data.frame(g = 3, s = c(1, 2, 1), t = c(2, 3, 3), ref = 2,
    control_group = "notyettreated", horizon = 3, profile_id = "same fitted label", composition_corrected = TRUE)
  attr(l, "profile_rates") <- matrix(1, 3, 2)
  expect_error(cdid:::build_chained_design(l), class = "cdid_design_functional_error")
  expect_error(cdid:::validate_functional_metadata(l), class = "cdid_design_input_error")
  l$functional_key <- c("declared:A", "declared:B", "declared:C")
  before <- l
  expect_true(cdid:::validate_functional_metadata(l))
  expect_identical(l, before)
  d <- cdid:::build_chained_design(l)
  expect_identical(sum(d$parameters$role == "pre"), 2L)
  # Equal fitted rates and boolean correction did not merge A and C.
  expect_identical(d$parameters$functional_key[d$parameters$role == "pre"], c("declared:A", "declared:C"))
})

test_that("bounded N/T smoke retains compact indices and splits one long unit", {
  # 1000 sparse units on a 1200-date calendar, plus a fully observed long unit.
  # This finite fixture checks storage/batch bounds, not universal performance.
  n <- 1000L; tt <- 1200L
  sparse <- data.frame(id = rep(seq_len(n), each = 4L),
    time = rep(c(1, 401, 800, tt), n), y = rep(c(0, 1, 3, 4), n), g = Inf)
  long <- data.frame(id = n + 1L, time = seq_len(tt), y = seq_len(tt), g = Inf)
  index <- cdid:::construct_observed_pairs(rbind(sparse, long), "id", "time", "y", "g")
  expect_lt(as.double(object.size(index)), 2e6)
  visited <- 0; maximum <- 0; long_batches <- 0
  got <- cdid:::traverse_observed_pairs(index, 257, function(b) {
    visited <<- visited + nrow(b); maximum <<- max(maximum, nrow(b))
    if (any(b$id == n + 1L)) long_batches <<- long_batches + 1L
  })
  expect_identical(visited, n * 6 + tt * (tt - 1) / 2)
  expect_identical(got$n_pairs, visited); expect_identical(maximum, 257)
  expect_gt(long_batches, 1L)
  # 4000-link chain: a former dense/path representation would exceed 16M
  # entries per matrix. No such matrix is retained or constructed here.
  m <- 4000L
  l <- data.frame(g = 1, s = 0:(m - 1L), t = 1:m, ref = 0, functional_key = "declared")
  d <- cdid:::build_chained_design(l)
  expect_lt(as.double(object.size(d)), 6e6)
  expect_identical(d$rank_diagnostics$rank, m)
  expect_null(d$W); expect_null(d$rank_diagnostics$root_paths)
  expect_error(cdid:::materialize_chained_design(d, 1000), class = "cdid_design_bound_error")
  cert <- cdid:::check_level_estimable(d, m)$witness_links
  expect_identical(nrow(cert), m)
})
