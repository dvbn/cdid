# CP08 Target 1: native bounded unit-row batch accumulator (`C_gmm_if_unit_batch`)
# for the unweighted nonclustered unit-row branch of `project_gmm_if_stream`.
sys.source(file.path("reference", "gmm_audit_support.R"), envir = environment())

cp08_pure_r_env <- local({
  ref_path <- file.path("reference", "chained_gmm_core_pure_r_reference.R")
  e <- new.env(parent = globalenv())
  for (fn in c("construct_chained_design", "root_path_certificate")) {
    if (exists(fn, envir = parent.env(environment()), inherits = TRUE)) {
      assign(fn, get(fn, envir = parent.env(environment()), inherits = TRUE), envir = e)
    }
  }
  sys.source(ref_path, envir = e)
  e
})

cp08_raw_bits <- function(x) {
  writeBin(as.double(x), raw(), size = 8L, endian = "little")
}

cp08_unit_replay <- function(Phi, batch_indices = list(seq_len(nrow(Phi))),
                             centered = TRUE, weighted = FALSE,
                             weight_vec = NULL, weight_type = if (weighted) "custom_nonnegative" else "uniform",
                             norm_factor = 1 / nrow(Phi),
                             unit_ids = sprintf("u%02d", seq_len(nrow(Phi))),
                             expected_unit_ids = unit_ids,
                             nuisance_fit_id = "frozen_assumed",
                             mutate_pass2 = NULL) {
  n <- nrow(Phi)
  pass_no <- 0L
  cb <- function(consume) {
    pass_no <<- pass_no + 1L
    splits <- if (!is.null(mutate_pass2) && pass_no == 2L && !is.null(mutate_pass2$batch_indices)) {
      mutate_pass2$batch_indices
    } else {
      batch_indices
    }
    for (b_idx in seq_along(splits)) {
      idx <- splits[[b_idx]]
      sub_phi <- Phi[idx, , drop = FALSE]
      sub_u <- unit_ids[idx]
      sub_w <- if (!is.null(weight_vec)) weight_vec[idx] else NULL
      if (!is.null(mutate_pass2) && pass_no == 2L && is.function(mutate_pass2$batch_fn)) {
        mut <- mutate_pass2$batch_fn(list(unit_id = sub_u, phi = sub_phi, weight = sub_w), b_idx)
        sub_u <- mut$unit_id
        sub_phi <- mut$phi
        sub_w <- mut$weight
      }
      batch <- list(unit_id = sub_u, phi = sub_phi, nuisance_fit_id = nuisance_fit_id)
      if (!is.null(sub_w)) batch$weight <- sub_w
      consume(batch)
    }
  }
  n_pos <- if (weighted && !is.null(weight_vec)) sum(weight_vec > 0) else n
  gmm_if_replay_spec(
    callback = cb,
    input_row_level = "unit",
    aggregate_to_clusters = FALSE,
    n_input_rows = n,
    n_output_rows = n,
    n_positive_contributing_rows = n_pos,
    centered = centered,
    weighted = weighted,
    weight_type = weight_type,
    normalization_factor = norm_factor,
    scale_id = "canonical_n_inv",
    nuisance_fit_id = nuisance_fit_id,
    expected_unit_ids = expected_unit_ids
  )
}

test_that("CP08 unit batch accumulator matches pure-R reference bit-for-bit across chunk sizes, centering, names, and materialization", {
  fx <- audit_fixture()
  pure_r <- cp08_pure_r_env
  b <- audit_bounds()
  p <- audit_provenance()
  Lambda_named <- matrix(
    c(0.75, -0.25, 0.5, 0.125,
      -0.5, 1.25, -0.75, 0.375,
      0.25, 0.625, -0.125, -0.875),
    nrow = 4L, ncol = 3L,
    dimnames = list(paste0("m", 1:4), c("t1", "t2", "t3"))
  )
  Lambda_unnamed <- unname(Lambda_named)

  split_patterns <- list(
    as.list(seq_len(8L)),
    list(1:2, 3:4, 5:6, 7:8),
    list(1:3, 4:6, 7:8),
    list(1L, 2:4, 5L, 6:8),
    list(1:5, 6:8),
    list(1:8)
  )

  for (centered in c(TRUE, FALSE)) {
    for (mat_if in c(FALSE, TRUE)) {
      for (use_names in c(TRUE, FALSE)) {
        Lam <- if (use_names) Lambda_named else Lambda_unnamed
        ref_rp <- cp08_unit_replay(fx$Phi, batch_indices = list(1:8), centered = centered)
        ref_out <- pure_r$project_gmm_if_stream(
          ref_rp, Lam, 4L, b, p, materialize_target_if = mat_if
        )
        for (sp in split_patterns) {
          rp <- cp08_unit_replay(fx$Phi, batch_indices = sp, centered = centered)
          expect_true(gmm_if_unit_batch_eligible(rp, ncol(Lam)))
          out <- project_gmm_if_stream(
            rp, Lam, 4L, b, p, materialize_target_if = mat_if
          )
          expect_identical(out$n_input_rows, ref_out$n_input_rows)
          expect_identical(out$n_output_rows, ref_out$n_output_rows)
          expect_identical(out$n_positive_contributing_rows, 8L)
          expect_identical(out$n_batches, length(sp))
          expect_identical(out$max_batch_rows_seen, max(vapply(sp, length, integer(1L))))
          expect_identical(out$total_weight, ref_out$total_weight)
          expect_identical(cp08_raw_bits(out$psi_weighted_mean), cp08_raw_bits(ref_out$psi_weighted_mean))
          expect_identical(names(out$psi_weighted_mean), names(ref_out$psi_weighted_mean))
          expect_identical(cp08_raw_bits(out$psi_raw_sum), cp08_raw_bits(ref_out$psi_raw_sum))
          expect_identical(names(out$psi_raw_sum), names(ref_out$psi_raw_sum))
          expect_identical(cp08_raw_bits(out$sandwich_cov), cp08_raw_bits(ref_out$sandwich_cov))
          expect_identical(dimnames(out$sandwich_cov), dimnames(ref_out$sandwich_cov))
          expect_identical(cp08_raw_bits(out$se_sandwich), cp08_raw_bits(ref_out$se_sandwich))
          expect_identical(names(out$se_sandwich), names(ref_out$se_sandwich))
          expect_identical(out$replay_verified_scope, ref_out$replay_verified_scope)
          expect_identical(out$sampling_provenance_certified, FALSE)
          if (mat_if) {
            expect_identical(cp08_raw_bits(out$materialized_target_if), cp08_raw_bits(ref_out$materialized_target_if))
            expect_identical(dimnames(out$materialized_target_if), dimnames(ref_out$materialized_target_if))
          } else {
            expect_null(out$materialized_target_if)
          }
        }
      }
    }
  }
})

test_that("CP08 unit batch accumulator preserves invariants across row/batch permutations and executes both producer passes", {
  pure_r <- cp08_pure_r_env
  b <- audit_bounds()
  p <- audit_provenance()
  # Dyadic coordinates retain exact raw sums. Centering uses 1/6, so
  # covariance and SE follow the frozen reference in each actual row order.
  Phi_dy <- matrix(
    c(0.5, -0.25, 1.0, 0.125,
      -1.5, 0.75, -0.5, 0.25,
      0.25, 0.5, -1.25, -0.5,
      1.0, -1.0, 0.75, -0.125,
      -0.75, 0.25, 0.5, 0.375,
      0.125, -0.5, -0.25, 0.75),
    nrow = 6L, ncol = 4L, byrow = TRUE
  )
  Lambda_dy <- matrix(
    c(1.0, -0.5,
      0.5, 0.25,
      -0.25, 1.0,
      0.75, -0.125),
    nrow = 4L, ncol = 2L, byrow = TRUE,
    dimnames = list(NULL, c("tau1", "tau2"))
  )
  uids <- sprintf("unit_%02d", seq_len(6L))
  base_rp <- cp08_unit_replay(Phi_dy, batch_indices = list(1:2, 3:4, 5:6), unit_ids = uids)
  base_out <- project_gmm_if_stream(base_rp, Lambda_dy, 4L, b, p, materialize_target_if = TRUE)

  perms <- list(
    c(6L, 5L, 4L, 3L, 2L, 1L),
    c(2L, 4L, 6L, 1L, 3L, 5L),
    c(3L, 1L, 5L, 2L, 6L, 4L)
  )
  for (perm in perms) {
    for (splits in list(as.list(1:6), list(1:3, 4:6), list(1:6))) {
      pass_calls <- integer(0)
      rp_perm <- cp08_unit_replay(
        Phi_dy[perm, , drop = FALSE],
        batch_indices = splits,
        unit_ids = uids[perm],
        expected_unit_ids = uids[perm]
      )
      orig_cb <- rp_perm$callback
      rp_perm$callback <- function(consume) {
        pass_calls <<- c(pass_calls, 1L)
        orig_cb(consume)
      }
      out_perm <- project_gmm_if_stream(rp_perm, Lambda_dy, 4L, b, p, materialize_target_if = TRUE)
      ref_perm <- pure_r$project_gmm_if_stream(rp_perm, Lambda_dy, 4L, b, p, materialize_target_if = TRUE)
      # Both producer passes executed on every call
      expect_identical(length(pass_calls), 4L)
      expect_identical(out_perm$n_positive_contributing_rows, 6L)
      expect_identical(out_perm$total_weight, base_out$total_weight)
      expect_equal(out_perm$psi_weighted_mean, base_out$psi_weighted_mean, tolerance = 0)
      expect_equal(out_perm$psi_raw_sum, base_out$psi_raw_sum, tolerance = 0)
      expect_equal(out_perm$sandwich_cov, base_out$sandwich_cov, tolerance = 1e-12)
      expect_equal(out_perm$se_sandwich, base_out$se_sandwich, tolerance = 1e-12)
      expect_equal(
        out_perm$materialized_target_if[uids, , drop = FALSE],
        base_out$materialized_target_if,
        tolerance = 0
      )
      expect_identical(cp08_raw_bits(out_perm$sandwich_cov), cp08_raw_bits(ref_perm$sandwich_cov))
      expect_identical(cp08_raw_bits(out_perm$se_sandwich), cp08_raw_bits(ref_perm$se_sandwich))
    }
  }
})

test_that("CP08 unit batch accumulator refuses zero/negative/nonfinite weights, nonfinite rows, and overflow in exact row order", {
  fx <- audit_fixture()
  b <- audit_bounds()
  p <- audit_provenance()
  Lam <- diag(4L)

  # Unweighted stream rejects non-unit weights (0, negative, >1, NA, Inf)
  for (bad_w in list(
    c(0, rep(1, 7)),
    c(-1, rep(1, 7)),
    c(2, rep(1, 7)),
    c(NA_real_, rep(1, 7)),
    c(Inf, rep(1, 7))
  )) {
    rp_bad_w <- cp08_unit_replay(fx$Phi, weight_vec = bad_w, weighted = FALSE)
    expect_error(
      project_gmm_if_stream(rp_bad_w, Lam, 4L, b, p),
      class = "cdid_gmm_replay_error"
    )
  }

  # Nonfinite phi and projected IF overflow in stream
  for (bad_val in c(NA_real_, NaN, Inf, -Inf)) {
    Phi_bad <- fx$Phi
    Phi_bad[3L, 2L] <- bad_val
    rp_bad <- cp08_unit_replay(Phi_bad, batch_indices = list(1:4, 5:8))
    expect_error(
      project_gmm_if_stream(rp_bad, Lam, 4L, b, p),
      class = "cdid_gmm_replay_error"
    )
  }
  Phi_proj_ovf <- matrix(1e308, nrow = 4L, ncol = 4L)
  rp_proj_ovf <- cp08_unit_replay(Phi_proj_ovf)
  expect_error(
    project_gmm_if_stream(rp_proj_ovf, matrix(10, 4L, 2L), 4L, b, p),
    regexp = "Projected IF overflowed",
    class = "cdid_gmm_replay_error"
  )

  # Pass 1 mean/sum accumulation overflow across rows
  Phi_sum_ovf <- matrix(1e308, nrow = 3L, ncol = 2L)
  rp_sum_ovf <- cp08_unit_replay(Phi_sum_ovf, batch_indices = list(1:3))
  expect_error(
    project_gmm_if_stream(rp_sum_ovf, diag(2L), 2L, b, p),
    regexp = "Replay mean/weight accumulation overflowed",
    class = "cdid_gmm_replay_error"
  )

  # Pass 2 covariance overflow (`z_j * z_k` overflows while Pass 1 sums stay finite)
  Phi_cov_ovf <- matrix(c(1e200, -1e200, 1e200, -1e200), nrow = 2L, ncol = 2L, byrow = TRUE)
  rp_cov_ovf <- cp08_unit_replay(Phi_cov_ovf, batch_indices = list(1:2), centered = FALSE)
  expect_error(
    project_gmm_if_stream(rp_cov_ovf, diag(2L), 2L, b, p),
    regexp = "Replay covariance overflowed",
    class = "cdid_gmm_replay_error"
  )

  # Direct native batch checks: zero weight vs negative/nonfinite weights and exact row refusal order
  psi_ok <- matrix(c(1.5, -2.0, 0.5, 3.0), nrow = 2L, ncol = 2L)
  res_zero_w <- gmm_if_unit_batch(
    1L, psi_ok, c(0.0, 2.0), c(FALSE, FALSE), c(0.0, 2.0, 0.0, 0.0),
    weighted_sum = c(0.0, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 1000
  )
  expect_identical(res_zero_w$counts, c(2.0, 2.0, 1.0, 2.0))
  expect_identical(res_zero_w$weighted_sum, c(2.0 * (-2.0), 2.0 * 3.0))
  expect_identical(res_zero_w$raw_sum, c(1.5 + (-2.0), 0.5 + 3.0))

  for (bad_w_vec in list(c(-0.1, 1.0), c(NA_real_, 1.0), c(NaN, 1.0), c(Inf, 1.0))) {
    expect_error(
      gmm_if_unit_batch(
        1L, psi_ok, bad_w_vec, c(FALSE, FALSE), c(0.0, 2.0, 0.0, 0.0),
        weighted_sum = c(0.0, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 1000
      ),
      regexp = "Invalid replay weights",
      class = "cdid_gmm_replay_error"
    )
  }

  # Row-by-row failure order inside a single batch:
  # 1. `nout >= no` checked before nonfinite `psi` on the exceeding row
  psi_nonfin <- matrix(c(Inf, 1.0, 1.0, 1.0), nrow = 2L, ncol = 2L)
  expect_error(
    gmm_if_unit_batch(
      1L, psi_nonfin, c(1.0, 1.0), c(TRUE, TRUE), c(2.0, 2.0, 2.0, 2.0),
      weighted_sum = c(0.0, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 1000
    ),
    regexp = "Replay exceeds declared output count",
    class = "cdid_gmm_replay_error"
  )
  # 2. Nonfinite `psi[1, ]` checked before `flags[1] == TRUE`
  expect_error(
    gmm_if_unit_batch(
      1L, psi_nonfin, c(1.0, 1.0), c(TRUE, FALSE), c(0.0, 2.0, 0.0, 0.0),
      weighted_sum = c(0.0, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 1000
    ),
    regexp = "Cluster sum overflowed",
    class = "cdid_gmm_replay_error"
  )
  # 3. Row 1 sum overflow triggers before Row 2 duplicate flag (`flags[2] == TRUE`)
  psi_r1_ovf <- matrix(c(1e308, 1.0, 1.0, 1.0), nrow = 2L, ncol = 2L)
  expect_error(
    gmm_if_unit_batch(
      1L, psi_r1_ovf, c(1.0, 1.0), c(FALSE, TRUE), c(0.0, 2.0, 0.0, 0.0),
      weighted_sum = c(1e308, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 1000
    ),
    regexp = "Replay mean/weight accumulation overflowed",
    class = "cdid_gmm_replay_error"
  )
  # 4. Row 1 duplicate flag (`flags[1] == TRUE`) triggers before Row 2 sum overflow
  psi_r2_ovf <- matrix(c(1.0, 1e308, 1.0, 1.0), nrow = 2L, ncol = 2L)
  expect_error(
    gmm_if_unit_batch(
      1L, psi_r2_ovf, c(1.0, 1.0), c(TRUE, FALSE), c(0.0, 2.0, 0.0, 0.0),
      weighted_sum = c(1e308, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 1000
    ),
    regexp = "Closed cluster reappeared",
    class = "cdid_gmm_replay_error"
  )
  # 5. Pass 2: Row 1 metadata flag (`flags[1] == TRUE`) triggers before Row 2 covariance overflow
  psi_p2 <- matrix(c(1.0, 1e200, 1.0, 1.0), nrow = 2L, ncol = 2L)
  expect_error(
    gmm_if_unit_batch(
      2L, psi_p2, c(1.0, 1.0), c(TRUE, FALSE), c(0.0, 2.0, 2.0, 2.0),
      psi_mean = c(0.0, 0.0), S = matrix(0.0, 2L, 2L), max_matrix_entries = 1000
    ),
    regexp = "Replay output metadata changed",
    class = "cdid_gmm_replay_error"
  )
  # 6. Pass 2: Row 1 covariance overflow triggers before Row 2 metadata flag (`flags[2] == TRUE`)
  psi_p2_rev <- matrix(c(1e200, 1.0, 1.0, 1.0), nrow = 2L, ncol = 2L)
  expect_error(
    gmm_if_unit_batch(
      2L, psi_p2_rev, c(1.0, 1.0), c(FALSE, TRUE), c(0.0, 2.0, 2.0, 2.0),
      psi_mean = c(0.0, 0.0), S = matrix(0.0, 2L, 2L), max_matrix_entries = 1000
    ),
    regexp = "Replay covariance overflowed",
    class = "cdid_gmm_replay_error"
  )
})

test_that("CP08 unit batch accumulator refuses duplicate IDs, wrong counts, and mutated second-pass rows/weights/IDs/batches", {
  fx <- audit_fixture()
  b <- audit_bounds()
  p <- audit_provenance()
  Lam <- diag(4L)

  # Duplicate unit IDs within the same batch
  uids_dup_in <- sprintf("u%02d", seq_len(8L))
  uids_dup_in[3L] <- uids_dup_in[1L]
  rp_dup_in <- cp08_unit_replay(fx$Phi, batch_indices = list(1:4, 5:8), unit_ids = uids_dup_in)
  expect_error(
    project_gmm_if_stream(rp_dup_in, Lam, 4L, b, p),
    regexp = "Closed cluster reappeared",
    class = "cdid_gmm_replay_error"
  )

  # Duplicate unit IDs across batches
  uids_dup_cross <- sprintf("u%02d", seq_len(8L))
  uids_dup_cross[6L] <- uids_dup_cross[2L]
  rp_dup_cross <- cp08_unit_replay(fx$Phi, batch_indices = list(1:4, 5:8), unit_ids = uids_dup_cross)
  expect_error(
    project_gmm_if_stream(rp_dup_cross, Lam, 4L, b, p),
    regexp = "Closed cluster reappeared",
    class = "cdid_gmm_replay_error"
  )

  # Wrong declared output/input count or callback emitting fewer rows
  expect_error(
    gmm_if_replay_spec(
      callback = function(consume) NULL,
      input_row_level = "unit", aggregate_to_clusters = FALSE,
      n_input_rows = 8L, n_output_rows = 7L, n_positive_contributing_rows = 8L,
      centered = TRUE, weighted = FALSE, weight_type = "uniform",
      normalization_factor = 1 / 8, scale_id = "s", nuisance_fit_id = "frozen_assumed"
    ),
    class = "cdid_gmm_replay_error"
  )
  rp_short <- cp08_unit_replay(fx$Phi, batch_indices = list(1:4, 5:7))
  expect_error(
    project_gmm_if_stream(rp_short, Lam, 4L, b, p),
    regexp = "Replay counts/batches differ",
    class = "cdid_gmm_replay_error"
  )
  rp_long <- cp08_unit_replay(fx$Phi, batch_indices = list(1:4, 5:8, 1L))
  expect_error(
    project_gmm_if_stream(rp_long, Lam, 4L, b, p),
    regexp = "Replay exceeds declared input count",
    class = "cdid_gmm_replay_error"
  )

  # Expected unit ID mismatch and nuisance label mismatch
  uids_exp_bad <- sprintf("u%02d", seq_len(8L))
  uids_exp_bad[8L] <- "u_wrong"
  rp_exp_bad <- cp08_unit_replay(fx$Phi, expected_unit_ids = uids_exp_bad)
  expect_error(
    project_gmm_if_stream(rp_exp_bad, Lam, 4L, b, p),
    regexp = "Expected unit sequence differs",
    class = "cdid_gmm_replay_error"
  )
  rp_nuis_bad <- cp08_unit_replay(fx$Phi, nuisance_fit_id = "other_nuisance")
  expect_error(
    project_gmm_if_stream(rp_nuis_bad, Lam, 4L, b, p),
    regexp = "Replay nuisance label differs",
    class = "cdid_gmm_replay_error"
  )

  # Second-pass mutations: projected rows, weights, unit IDs, batch sizes, and batch counts
  rp_mut_row <- cp08_unit_replay(
    fx$Phi, batch_indices = list(1:4, 5:8),
    mutate_pass2 = list(batch_fn = function(batch, b_idx) {
      if (b_idx == 2L) batch$phi[1L, ] <- batch$phi[1L, ] + c(1, -1, 0, 0)
      batch
    })
  )
  expect_error(
    project_gmm_if_stream(rp_mut_row, Lam, 4L, b, p),
    regexp = "IF replay changed projected rows, weights, IDs or batch structure",
    class = "cdid_gmm_replay_error"
  )

  rp_mut_id <- cp08_unit_replay(
    fx$Phi, batch_indices = list(1:4, 5:8),
    mutate_pass2 = list(batch_fn = function(batch, b_idx) {
      if (b_idx == 1L) batch$unit_id[2L] <- "u_mutated"
      batch
    })
  )
  expect_error(
    project_gmm_if_stream(rp_mut_id, Lam, 4L, b, p),
    regexp = "IF replay changed projected rows, weights, IDs or batch structure",
    class = "cdid_gmm_replay_error"
  )

  rp_mut_w <- cp08_unit_replay(
    fx$Phi, batch_indices = list(1:4, 5:8),
    mutate_pass2 = list(batch_fn = function(batch, b_idx) {
      if (b_idx == 1L) batch$weight <- c(2, 1, 1, 1)
      batch
    })
  )
  expect_error(
    project_gmm_if_stream(rp_mut_w, Lam, 4L, b, p),
    class = "cdid_gmm_replay_error"
  )

  rp_mut_bsize <- cp08_unit_replay(
    fx$Phi, batch_indices = list(1:4, 5:8),
    mutate_pass2 = list(batch_indices = list(1:3, 4:8))
  )
  expect_error(
    project_gmm_if_stream(rp_mut_bsize, Lam, 4L, b, p),
    regexp = "IF replay changed projected rows, weights, IDs or batch structure",
    class = "cdid_gmm_replay_error"
  )

  rp_mut_bcount <- cp08_unit_replay(
    fx$Phi, batch_indices = list(1:4, 5:8),
    mutate_pass2 = list(batch_indices = list(1:8))
  )
  expect_error(
    project_gmm_if_stream(rp_mut_bcount, Lam, 4L, b, p),
    regexp = "IF replay changed projected rows, weights, IDs or batch structure",
    class = "cdid_gmm_replay_error"
  )
})

test_that("CP08 unit batch accumulator preserves input SEXPs and handles read-only aliased inputs without mutation", {
  psi <- matrix(c(1.25, -0.5, 2.0, 0.75, -1.5, 0.25), nrow = 3L, ncol = 2L,
                dimnames = list(c("r1", "r2", "r3"), c("c1", "c2")))
  w <- c(1.0, 1.0, 1.0)
  flags <- c(FALSE, FALSE, FALSE)
  counts1 <- c(0.0, 3.0, 0.0, 0.0)
  ws_in <- c(c1 = 0.5, c2 = -0.25)
  rs_in <- c(c1 = 1.0, c2 = 2.0)

  psi_snap <- cp08_raw_bits(psi)
  w_snap <- cp08_raw_bits(w)
  counts1_snap <- cp08_raw_bits(counts1)
  ws_snap <- cp08_raw_bits(ws_in)
  rs_snap <- cp08_raw_bits(rs_in)

  r1 <- gmm_if_unit_batch(
    1L, psi, w, flags, counts1,
    weighted_sum = ws_in, raw_sum = rs_in, max_matrix_entries = 1000
  )
  expect_identical(cp08_raw_bits(psi), psi_snap)
  expect_identical(cp08_raw_bits(w), w_snap)
  expect_identical(cp08_raw_bits(counts1), counts1_snap)
  expect_identical(cp08_raw_bits(ws_in), ws_snap)
  expect_identical(cp08_raw_bits(rs_in), rs_snap)
  expect_false(gmm_sexp_aliases(r1$counts, counts1))
  expect_false(gmm_sexp_aliases(r1$weighted_sum, ws_in))
  expect_false(gmm_sexp_aliases(r1$raw_sum, rs_in))
  expect_false(gmm_sexp_aliases(r1$weighted_sum, r1$raw_sum))

  # Read-only aliased inputs (`weighted_sum` and `raw_sum` sharing one SEXP)
  shared_sum <- c(c1 = 0.25, c2 = -0.75)
  shared_snap <- cp08_raw_bits(shared_sum)
  expect_true(gmm_sexp_aliases(shared_sum, shared_sum))
  r1_alias <- gmm_if_unit_batch(
    1L, psi, w, flags, counts1,
    weighted_sum = shared_sum, raw_sum = shared_sum, max_matrix_entries = 1000
  )
  expect_identical(cp08_raw_bits(shared_sum), shared_snap)
  expect_false(gmm_sexp_aliases(r1_alias$weighted_sum, shared_sum))
  expect_false(gmm_sexp_aliases(r1_alias$raw_sum, shared_sum))
  expect_false(gmm_sexp_aliases(r1_alias$weighted_sum, r1_alias$raw_sum))

  # Pass 2 input immutability and output non-aliasing
  counts2 <- c(0.0, 3.0, 3.0, 3.0)
  pm_in <- c(c1 = 0.25, c2 = 0.5)
  S_in <- matrix(c(1.0, 0.25, 0.25, 2.0), nrow = 2L, ncol = 2L,
                 dimnames = list(c("c1", "c2"), c("c1", "c2")))
  counts2_snap <- cp08_raw_bits(counts2)
  pm_snap <- cp08_raw_bits(pm_in)
  S_snap <- cp08_raw_bits(S_in)

  r2 <- gmm_if_unit_batch(
    2L, psi, w, flags, counts2,
    psi_mean = pm_in, S = S_in, max_matrix_entries = 1000
  )
  expect_identical(cp08_raw_bits(counts2), counts2_snap)
  expect_identical(cp08_raw_bits(pm_in), pm_snap)
  expect_identical(cp08_raw_bits(S_in), S_snap)
  expect_false(gmm_sexp_aliases(r2$counts, counts2))
  expect_false(gmm_sexp_aliases(r2$S, S_in))
  expect_identical(dimnames(r2$S), dimnames(S_in))
})

test_that("CP08 weighted and clustered replay branches stay on the generic path before computation", {
  fx <- audit_fixture()
  pure_r <- cp08_pure_r_env
  b <- audit_bounds()
  p <- audit_provenance()
  Lam <- matrix(c(1, 0, -1, 0.5, 0, 1, 0.25, -0.5), nrow = 4L, ncol = 2L)

  # Weighted unit-row replay (even with unit weights) is ineligible for unit-batch fast path
  rp_weighted <- cp08_unit_replay(
    fx$Phi, batch_indices = list(1:4, 5:8),
    weighted = TRUE, weight_vec = c(0.5, 1.5, 0.0, 2.0, 1.0, 0.75, 1.25, 1.0),
    weight_type = "custom_nonnegative"
  )
  expect_false(gmm_if_unit_batch_eligible(rp_weighted, ncol(Lam)))

  rp_weighted_ones <- cp08_unit_replay(
    fx$Phi, batch_indices = list(1:4, 5:8),
    weighted = TRUE, weight_vec = rep(1.0, 8L),
    weight_type = "custom_nonnegative"
  )
  expect_false(gmm_if_unit_batch_eligible(rp_weighted_ones, ncol(Lam)))

  # Clustered replays (unit-to-cluster aggregation and pre-aggregated cluster input) are ineligible
  rp_cluster_agg <- audit_replay(fx$Phi, cluster = TRUE, batch_size = 3L)
  expect_false(gmm_if_unit_batch_eligible(rp_cluster_agg, ncol(Lam)))

  rp_cluster_in <- gmm_if_replay_spec(
    callback = function(consume) {
      consume(list(
        unit_id = paste0("c", 1:4),
        cluster_id = paste0("c", 1:4),
        phi = fx$Phi[1:4, , drop = FALSE],
        nuisance_fit_id = "frozen_assumed"
      ))
    },
    input_row_level = "cluster",
    aggregate_to_clusters = FALSE,
    n_input_rows = 4L,
    n_output_rows = 4L,
    n_positive_contributing_rows = 4L,
    centered = TRUE,
    weighted = FALSE,
    weight_type = "uniform",
    normalization_factor = 1 / 4,
    scale_id = "canonical_n_inv",
    nuisance_fit_id = "frozen_assumed"
  )
  expect_false(gmm_if_unit_batch_eligible(rp_cluster_in, ncol(Lam)))

  # A trap in the new resolver proves these branches never request the new
  # kernel, without mutating the established generic nine-symbol cache.
  generic_env <- new.env(parent = environment(project_gmm_if_stream))
  generic_project <- project_gmm_if_stream
  environment(generic_project) <- generic_env
  generic_env$gmm_require_unit_batch_symbol <- function(...) stop("Ineligible branch requested unit kernel")

  for (rp_generic in list(rp_weighted, rp_weighted_ones, rp_cluster_agg, rp_cluster_in)) {
    out_gen <- generic_project(rp_generic, Lam, 4L, b, p, materialize_target_if = TRUE)
    ref_gen <- pure_r$project_gmm_if_stream(rp_generic, Lam, 4L, b, p, materialize_target_if = TRUE)
    expect_identical(cp08_raw_bits(out_gen$sandwich_cov), cp08_raw_bits(ref_gen$sandwich_cov))
    expect_identical(cp08_raw_bits(out_gen$psi_weighted_mean), cp08_raw_bits(ref_gen$psi_weighted_mean))
    expect_identical(cp08_raw_bits(out_gen$psi_raw_sum), cp08_raw_bits(ref_gen$psi_raw_sum))
    expect_identical(cp08_raw_bits(out_gen$materialized_target_if), cp08_raw_bits(ref_gen$materialized_target_if))
  }
})

test_that("CP08 native symbol registration, arity, forgery refusal, and simultaneous buffer bounds are enforced", {
  dll <- base::getLoadedDLLs()[["cdid"]]
  expect_false(is.null(dll))
  reg <- base::getDLLRegisteredRoutines(dll, addNames = TRUE)[[".Call"]]
  expect_true("C_gmm_if_unit_batch" %in% names(reg))
  expect_identical(as.integer(reg[["C_gmm_if_unit_batch"]]$numParameters), 10L)
  si <- base::getNativeSymbolInfo("C_gmm_if_unit_batch", PACKAGE = dll)
  expect_s3_class(si, "NativeSymbolInfo")
  expect_identical(as.integer(si$numParameters), 10L)

  syms <- gmm_require_native_dll()
  expect_length(syms, 9L)
  expect_false("C_gmm_if_unit_batch" %in% names(syms))
  unit_sym <- gmm_require_unit_batch_symbol(syms)
  expect_identical(unit_sym, si$address)
  expect_identical(gmm_require_unit_batch_symbol(unit_sym), si$address)

  # Forged symbol inputs (NULL, non-pointer, wrong registered symbol, duplicate-aliased list)
  psi_small <- matrix(c(1.0, 2.0, 3.0, 4.0), nrow = 2L, ncol = 2L)
  for (forged in list(
    list(C_gmm_if_unit_batch = NULL),
    list(C_gmm_if_unit_batch = "C_gmm_if_unit_batch"),
    list(C_gmm_if_unit_batch = 10L),
    list(C_gmm_if_unit_batch = syms$C_gmm_next),
    list(C_gmm_if_unit_batch = unit_sym, C_gmm_other = unit_sym),
    syms$C_gmm_next
  )) {
    expect_error(
      gmm_require_unit_batch_symbol(forged),
      class = "cdid_gmm_native_error"
    )
    expect_error(
      gmm_if_unit_batch(
        1L, psi_small, c(1.0, 1.0), c(FALSE, FALSE), c(0.0, 2.0, 0.0, 0.0),
        weighted_sum = c(0.0, 0.0), raw_sum = c(0.0, 0.0),
        max_matrix_entries = 1000, syms = forged
      ),
      class = "cdid_gmm_native_error"
    )
  }

  # A supplied unit-symbol entry never becomes proof of registration. Keep
  # the original generic cache untouched and refuse forged per-call tables.
  for (bad_cached in list(NULL, "forged", syms$C_gmm_next)) {
    supplied <- syms
    supplied["C_gmm_if_unit_batch"] <- list(bad_cached)
    expect_error(gmm_require_unit_batch_symbol(supplied), class = "cdid_gmm_native_error")
  }
  expect_identical(gmm_require_native_dll(), syms)

  # Exact simultaneous buffer bounds in Pass 1 (`B*J + 2*B + 4*J + 8`) and Pass 2 (`B*J + 2*B + 2*J + 2*J^2 + 8`)
  # For B = 3, J = 2: Pass 1 exact bound = 6 + 6 + 8 + 8 = 28; Pass 2 exact bound = 6 + 6 + 4 + 8 + 8 = 32
  psi_b3 <- matrix(as.double(1:6), nrow = 3L, ncol = 2L)
  w_b3 <- c(1.0, 1.0, 1.0)
  fl_b3 <- c(FALSE, FALSE, FALSE)
  expect_error(
    gmm_if_unit_batch(
      1L, psi_b3, w_b3, fl_b3, c(0.0, 3.0, 0.0, 0.0),
      weighted_sum = c(0.0, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 27
    ),
    class = "cdid_gmm_bound_error"
  )
  expect_error(
    .Call(
      unit_sym, 1L, psi_b3, w_b3, fl_b3, c(0.0, 3.0, 0.0, 0.0),
      c(0.0, 0.0), c(0.0, 0.0), NULL, NULL, 27
    ),
    class = "cdid_gmm_bound_error"
  )
  ok_p1 <- gmm_if_unit_batch(
    1L, psi_b3, w_b3, fl_b3, c(0.0, 3.0, 0.0, 0.0),
    weighted_sum = c(0.0, 0.0), raw_sum = c(0.0, 0.0), max_matrix_entries = 28
  )
  expect_identical(ok_p1$counts, c(3.0, 3.0, 3.0, 3.0))

  expect_error(
    gmm_if_unit_batch(
      2L, psi_b3, w_b3, fl_b3, c(0.0, 3.0, 3.0, 3.0),
      psi_mean = c(0.0, 0.0), S = matrix(0.0, 2L, 2L), max_matrix_entries = 31
    ),
    class = "cdid_gmm_bound_error"
  )
  expect_error(
    .Call(
      unit_sym, 2L, psi_b3, w_b3, fl_b3, c(0.0, 3.0, 3.0, 3.0),
      NULL, NULL, c(0.0, 0.0), matrix(0.0, 2L, 2L), 31
    ),
    class = "cdid_gmm_bound_error"
  )
  ok_p2 <- gmm_if_unit_batch(
    2L, psi_b3, w_b3, fl_b3, c(0.0, 3.0, 3.0, 3.0),
    psi_mean = c(0.0, 0.0), S = matrix(0.0, 2L, 2L), max_matrix_entries = 32
  )
  expect_identical(ok_p2$counts, c(3.0, 3.0, 3.0, 3.0))

  # Stream preflight refuses insufficient `max_matrix_entries` before invoking callback
  fx <- audit_fixture()
  p <- audit_provenance()
  b_small <- audit_bounds()
  b_small$max_matrix_entries <- 35
  cb_called <- FALSE
  rp_pre <- cp08_unit_replay(fx$Phi)
  orig_cb <- rp_pre$callback
  rp_pre$callback <- function(consume) {
    cb_called <<- TRUE
    orig_cb(consume)
  }
  expect_error(
    project_gmm_if_stream(rp_pre, matrix(1, 4L, 3L), 4L, b_small, p),
    class = "cdid_gmm_bound_error"
  )
  expect_false(cb_called)
})

test_that("CP08 changed-accumulation mutant disagrees with ordered row-by-row unit-batch accumulation and fails bitwise check", {
  pure_r <- cp08_pure_r_env
  b <- audit_bounds()
  p <- audit_provenance()

  # Four-row stream where sequential row order `( (2^53 + 1) - 2^53 ) + 1 == 1`
  # differs from reordered summation `(2^53 - 2^53) + (1 + 1) == 2` and from
  # reordered/uncentered batch Gram `crossprod(psi) - n * tcrossprod(psi_mean)`.
  Phi_cancel <- matrix(
    c(2^53, 1.0,
      1.0, 2^27,
      -2^53, 1.0,
      1.0, -2^27),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )
  Lam_id <- diag(2L)
  rp <- cp08_unit_replay(Phi_cancel, batch_indices = list(1:4), centered = TRUE, norm_factor = 1.0)

  out_native <- project_gmm_if_stream(rp, Lam_id, 2L, b, p)
  out_pure_r <- pure_r$project_gmm_if_stream(rp, Lam_id, 2L, b, p)

  # Native batch accumulator matches original ordered pure-R implementation bit-for-bit
  expect_identical(cp08_raw_bits(out_native$psi_raw_sum), cp08_raw_bits(out_pure_r$psi_raw_sum))
  expect_identical(cp08_raw_bits(out_native$psi_weighted_mean), cp08_raw_bits(out_pure_r$psi_weighted_mean))
  expect_identical(cp08_raw_bits(out_native$sandwich_cov), cp08_raw_bits(out_pure_r$sandwich_cov))

  # Changed-accumulation mutant 1: reordered row summation (grouping large cancelling rows first)
  perm_mutant <- c(1L, 3L, 2L, 4L)
  mutant_raw_sum <- colSums(Phi_cancel[perm_mutant, , drop = FALSE])
  mutant_mean <- mutant_raw_sum / nrow(Phi_cancel)
  Z_mutant <- sweep(Phi_cancel[perm_mutant, , drop = FALSE], 2L, mutant_mean, `-`)
  S_mutant_reordered <- matrix(0.0, 2L, 2L)
  for (i in seq_len(nrow(Z_mutant))) {
    S_mutant_reordered <- S_mutant_reordered + tcrossprod(Z_mutant[i, ])
  }
  expect_false(identical(cp08_raw_bits(mutant_raw_sum), cp08_raw_bits(out_native$psi_raw_sum)))
  expect_false(identical(cp08_raw_bits(mutant_mean), cp08_raw_bits(out_native$psi_weighted_mean)))
  expect_false(identical(cp08_raw_bits(S_mutant_reordered), cp08_raw_bits(out_native$sandwich_cov)))

  # Changed-accumulation mutant 2: uncentered one-pass batch Gram `crossprod(psi) - n * tcrossprod(mean)`
  S_mutant_onepass <- crossprod(Phi_cancel) - nrow(Phi_cancel) * tcrossprod(out_native$psi_weighted_mean)
  expect_false(identical(cp08_raw_bits(S_mutant_onepass), cp08_raw_bits(out_native$sandwich_cov)))
})
