require_cdid_namespace <- function() {
  if (!isNamespaceLoaded("cdid") && !requireNamespace("cdid", quietly = TRUE)) {
    stop("Actual installed 'cdid' namespace is required; no fallback is permitted in packaged tests.", call. = FALSE)
  }
  ns <- asNamespace("cdid")
  if (!isNamespace(ns) || !environmentIsLocked(ns)) {
    stop("'cdid' namespace must be a locked R namespace.", call. = FALSE)
  }
  ns
}

# Require installed namespace immediately at file load (never skip)
require_cdid_namespace()

audit_owner_ctrl <- function(
  ns,
  max_chunk_rows = 4L,
  max_retained_entries = 100000L,
  max_matrix_entries = 100000L,
  max_work_ops = 500000L,
  max_output_entries = 100000L
) {
  ns$ipw_owner_control(
    max_chunk_rows = max_chunk_rows,
    max_retained_entries = max_retained_entries,
    max_matrix_entries = max_matrix_entries,
    max_work_ops = max_work_ops,
    max_output_entries = max_output_entries
  )
}

make_audit_panel_df <- function() {
  rbind(
    data.frame(id = c(1, 1, 1), time = c(1, 2, 3), y = c(1.0, 2.0, 3.5), g = c(2, 2, 2)),
    data.frame(id = c(2, 2, 2), time = c(1, 2, 3), y = c(1.5, NA_real_, 4.0), g = c(2, 2, 2)),
    data.frame(id = c(3, 3, 3), time = c(1, 2, 3), y = c(0.5, 1.5, 2.5), g = c(3, 3, 3)),
    data.frame(id = c(4, 4, 4), time = c(1, 2, 3), y = c(2.0, 3.0, 4.5), g = c(3, 3, 3)),
    data.frame(id = c(5, 5, 5), time = c(1, 2, 3), y = c(1.2, 2.2, 3.2), g = c(Inf, Inf, Inf)),
    data.frame(id = c(6, 6, 6), time = c(1, 2, 3), y = c(0.8, 1.8, 2.8), g = c(Inf, Inf, Inf))
  )
}

# ---------------------------------------------------------------------------
# Audit 1: Explicit caller controls, tiny-budget refusal BEFORE content scans,
# and real 2D matrix dimension limit (INT_MAX) vs 1D vector length limit
# ---------------------------------------------------------------------------

test_that("explicit control bounds are mandatory and preflight refuses tiny budgets and 2D matrix limits before scanning content", {
  ns <- require_cdid_namespace()
  ctrl <- audit_owner_ctrl(ns)
  df <- make_audit_panel_df()

  # Missing control fields or NULL control refused
  expect_error(ns$ipw_owner_control(max_chunk_rows = 4L), class = "cdid_ipw_owner_control_error")
  expect_error(ns$validate_ipw_owner_control(NULL), class = "cdid_ipw_owner_control_error")
  expect_error(ns$validate_ipw_owner_control(list(max_chunk_rows = 4L)), class = "cdid_ipw_owner_control_error")

  # Real 2D matrix dimension limit INT_MAX vs R_XLEN_T_MAX=2^52; exact count limit=2^53-1.
  big_ctrl <- ns$ipw_owner_control(
    max_chunk_rows = 2^40,
    max_retained_entries = 2^50,
    max_matrix_entries = 2^50,
    max_work_ops = 2^52,
    max_output_entries = 2^50
  )
  # 1D panel owner preflight with n_rows > INT_MAX (e.g. 2^31) is allowed as 1D vectors (matrix_entries = 0)
  pre_1d <- ns$preflight_ipw_panel_owner(n_rows = as.double(.Machine$integer.max) + 10, control = big_ctrl)
  expect_identical(pre_1d$retained_matrix_entries, 0)
  expect_identical(pre_1d$matrix_dim_limit, as.double(.Machine$integer.max))
  expect_identical(pre_1d$vector_length_limit, 2^52) # R 4.5.2 Rinternals.h R_XLEN_T_MAX; count range remains 2^53-1

  # 2D frame owner preflight with n_a_units > INT_MAX (real 2D matrix Z_a) is refused with cdid_ipw_owner_matrix_dim_error
  expect_error(
    ns$preflight_ipw_frame_owner(
      n_a_units = as.double(.Machine$integer.max) + 1,
      n_total_units = as.double(.Machine$integer.max) + 1,
      n_supplied_id = as.double(.Machine$integer.max) + 1,
      p_covariates = 2L,
      n_calendar_dates = 4L,
      n_a_observed_rows = as.double(.Machine$integer.max) + 2,
      control = big_ctrl
    ),
    class = "cdid_ipw_owner_matrix_dim_error"
  )

  # Poisoned data frame (contains NaN in y and time) refused by tiny budget BEFORE scanning content
  df_poison <- df
  df_poison$y[1L] <- NaN
  df_poison$time[2L] <- Inf
  tiny_ret <- audit_owner_ctrl(ns, max_retained_entries = 10L)
  expect_error(
    ns$construct_ipw_panel_owner(df_poison, "id", "time", "y", "g", control = tiny_ret),
    class = "cdid_ipw_owner_bound_error"
  )
  tiny_wo <- audit_owner_ctrl(ns, max_work_ops = 5L)
  expect_error(
    ns$construct_ipw_panel_owner(df_poison, "id", "time", "y", "g", control = tiny_wo),
    class = "cdid_ipw_owner_bound_error"
  )

  # Build valid panel and frame, then test preflight before content scan on frame, views, and queries
  panel <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl)
  Z_poison <- matrix(NaN, nrow = 6L, ncol = 2L)
  expect_error(
    ns$construct_ipw_frame_owner(panel, id = 1:6, Z = Z_poison, control = tiny_ret),
    class = "cdid_ipw_owner_bound_error"
  )

  frame <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), control = ctrl)
  # Nuisance view with poisoned beta0 = NaN refused by tiny budget BEFORE beta0 scan
  expect_error(
    ns$construct_ipw_nuisance_view(
      frame, g = 2, K = c(3, Inf), s = 1, t = 2,
      control = tiny_ret, beta0 = c(NaN, NaN)
    ),
    class = "cdid_ipw_owner_bound_error"
  )
  # Link view with poisoned loading_beta_candidate = NaN refused by tiny budget BEFORE candidate scan
  expect_error(
    ns$construct_ipw_link_view(
      frame, s = 1, t = 2, g = 2, K = c(3, Inf),
      control = tiny_ret, loading_beta_candidate = c(NaN, NaN)
    ),
    class = "cdid_ipw_owner_bound_error"
  )

  nu <- ns$construct_ipw_nuisance_view(frame, g = 2, K = c(3, Inf), s = 1, t = 2, control = ctrl)
  lk <- ns$construct_ipw_link_view(frame, s = 1, t = 2, nuisance_view = nu, control = ctrl)

  # Range query with poisoned expected_a_id = NaN refused by tiny chunk/output budget BEFORE expected_a_id content scan
  tiny_ch <- audit_owner_ctrl(ns, max_chunk_rows = 2L)
  expect_error(
    ns$query_ipw_view_chunk(lk, a_start = 1L, a_end = 4L, control = tiny_ch, expected_a_id = rep(NaN, 4L)),
    class = "cdid_ipw_owner_bound_error"
  )

  # Stream query refuses tiny work_ops budget BEFORE invoking callback even once
  cb_called <- 0L
  expect_error(
    ns$stream_ipw_view_chunks(
      lk, control = audit_owner_ctrl(ns, max_work_ops = 10L),
      callback = function(ch) cb_called <<- cb_called + 1L
    ),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(cb_called, 0L)
})

# ---------------------------------------------------------------------------
# Audit 2: Duplicates, cohort invariance, decorative attributes, and support refusals
# ---------------------------------------------------------------------------

test_that("duplicate (ID, time), within-unit cohort changes, decorative attributes, and zero support are refused", {
  ns <- require_cdid_namespace()
  ctrl <- audit_owner_ctrl(ns)
  df <- make_audit_panel_df()

  # Duplicate (ID, time) row
  df_dup <- rbind(df, data.frame(id = 1, time = 2, y = 9.0, g = 2))
  expect_error(
    ns$construct_ipw_panel_owner(df_dup, "id", "time", "y", "g", control = ctrl),
    class = "cdid_ipw_owner_duplicate_error"
  )

  # Within-unit cohort change
  df_chg <- df
  df_chg$g[df_chg$id == 1 & df_chg$time == 3] <- 3
  expect_error(
    ns$construct_ipw_panel_owner(df_chg, "id", "time", "y", "g", control = ctrl),
    class = "cdid_ipw_owner_cohort_change_error"
  )

  # Decorative/foreign attribute on input data or column
  df_attr <- df
  attr(df_attr, "certified") <- TRUE
  expect_error(
    ns$construct_ipw_panel_owner(df_attr, "id", "time", "y", "g", control = ctrl),
    class = "cdid_ipw_owner_foreign_attribute_error"
  )

  panel <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl)
  # Frame cohort disagreement
  expect_error(
    ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), frame_g = rep(2, 6), control = ctrl),
    class = "cdid_ipw_owner_cohort_change_error"
  )
  # Missing A-unit in frame id
  expect_error(
    ns$construct_ipw_frame_owner(panel, id = 1:5, Z = cbind(1, 1:5), control = ctrl),
    class = "cdid_ipw_owner_id_error"
  )
  # Duplicate id in frame
  expect_error(
    ns$construct_ipw_frame_owner(panel, id = c(1:5, 5), Z = cbind(1, 1:6), control = ctrl),
    class = "cdid_ipw_owner_id_error"
  )
  # Non-positive weight
  expect_error(
    ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), weights = c(1, 1, 0, 1, 1, 1), control = ctrl),
    class = "cdid_ipw_owner_weight_error"
  )

  frame <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), control = ctrl)
  # Unsupported cohort in K (e.g. cohort 99 has 0 units)
  expect_error(
    ns$construct_ipw_nuisance_view(frame, g = 2, K = c(99, Inf), s = 1, t = 2, control = ctrl),
    class = "cdid_ipw_owner_cohort_support_error"
  )
  # Unsupported date pair (e.g. s = 1, t = 99)
  expect_error(
    ns$construct_ipw_nuisance_view(frame, g = 2, K = c(3, Inf), s = 1, t = 99, control = ctrl),
    class = "cdid_ipw_owner_pair_support_error"
  )
})

# ---------------------------------------------------------------------------
# Audit 3: Sealed copy-on-modify snapshots, mutation refusals, foreign/order-mismatched
# owners/views, and label-only non-equivalence
# ---------------------------------------------------------------------------

test_that("sealed copy-on-modify snapshots reject mutations, foreign objects, order mismatches, and label-only equivalence", {
  ns <- require_cdid_namespace()
  ctrl <- audit_owner_ctrl(ns, max_chunk_rows = 3L)
  df <- make_audit_panel_df()

  panel <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl, owner_id = "shared_lbl")
  frame <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), control = ctrl, frame_id = "shared_lbl")
  nu <- ns$construct_ipw_nuisance_view(frame, g = 2, K = c(3, Inf), s = 1, t = 3, control = ctrl, nuisance_id = "shared_lbl")
  lk <- ns$construct_ipw_link_view(frame, s = 1, t = 3, nuisance_view = nu, control = ctrl, link_id = "shared_lbl")

  # Locked reference environments cannot be modified in place
  expect_true(environmentIsLocked(panel$ref))
  expect_true(bindingIsLocked("panel_owner", panel$ref))
  expect_true(environmentIsLocked(frame$ref))
  expect_true(bindingIsLocked("frame_owner", frame$ref))
  expect_true(environmentIsLocked(nu$ref))
  expect_true(bindingIsLocked("nuisance_view", nu$ref))
  expect_true(environmentIsLocked(lk$ref))
  expect_true(bindingIsLocked("link_view", lk$ref))

  # Copy-on-modify mutation of panel_owner is refused
  panel_mut <- panel
  panel_mut$observations$y[1L] <- panel_mut$observations$y[1L] + 1e-12
  expect_error(ns$validate_ipw_panel_owner(panel_mut), class = "cdid_ipw_owner_integrity_error")

  # Even -0.0 vs +0.0 mutation is caught by identical(..., num.eq = FALSE)
  panel_z <- panel
  panel_z$all_times[1L] <- -0.0
  expect_error(ns$validate_ipw_panel_owner(panel_z), class = "cdid_ipw_owner_integrity_error")

  # Copy-on-modify mutation of frame_owner, nuisance_view, and link_view is refused
  frame_mut <- frame
  frame_mut$w_a[2L] <- 99.0
  expect_error(ns$validate_ipw_frame_owner(frame_mut), class = "cdid_ipw_owner_integrity_error")

  nu_mut <- nu
  nu_mut$cohorts$K <- c(Inf, 3)
  expect_error(ns$validate_ipw_nuisance_view(nu_mut), class = "cdid_ipw_owner_integrity_error")

  lk_mut <- lk
  lk_mut$pair$t <- 2.0
  expect_error(ns$validate_ipw_link_view(lk_mut), class = "cdid_ipw_owner_integrity_error")

  # Label-only equivalence NEVER proves compatibility!
  # Two distinct frame owners with the exact same `frame_id = "shared_lbl"`:
  frame_same_lbl <- ns$construct_ipw_frame_owner(
    panel, id = 1:6, Z = cbind(1, (1:6) * 2), control = ctrl, frame_id = "shared_lbl"
  )
  expect_identical(frame$frame_id, frame_same_lbl$frame_id)
  expect_error(
    ns$verify_ipw_frame_compatibility(frame, frame_same_lbl, require_same_owner = TRUE),
    class = "cdid_ipw_owner_mismatch_error"
  )
  expect_error(
    ns$verify_ipw_frame_compatibility(frame, frame_same_lbl, require_same_owner = FALSE),
    class = "cdid_ipw_owner_mismatch_error"
  )
  expect_error(
    ns$validate_ipw_nuisance_view(nu, frame_owner = frame_same_lbl),
    class = "cdid_ipw_owner_mismatch_error"
  )

  # Views with reordered controls K = c(3, Inf) vs K = c(Inf, 3) or different masks are refused
  nu_reord <- ns$construct_ipw_nuisance_view(frame, g = 2, K = c(Inf, 3), s = 1, t = 3, control = ctrl)
  expect_error(
    ns$verify_ipw_view_mask_compatibility(nu_reord, lk, control = ctrl),
    class = "cdid_ipw_owner_cohort_error"
  )
  # Note: unit 2 is missing at t = 2, so on (1, 2) unit 2 is NOT observed, whereas on (1, 3) unit 2 IS observed!
  lk_12_no_verify <- ns$construct_ipw_link_view(
    frame, s = 1, t = 2, nuisance_view = nu, control = ctrl, verify_nuisance_mask = FALSE
  )
  expect_error(
    ns$verify_ipw_view_mask_compatibility(nu, lk_12_no_verify, control = ctrl),
    class = "cdid_ipw_owner_mask_mismatch_error"
  )
  expect_error(
    ns$construct_ipw_link_view(frame, s = 1, t = 2, nuisance_view = nu, control = ctrl, verify_nuisance_mask = TRUE),
    class = "cdid_ipw_owner_mask_mismatch_error"
  )

  # Order-mismatched expected_a_id on chunk query is refused without scanning outside the chunk
  expect_error(
    ns$query_ipw_view_chunk(lk, a_start = 1L, a_end = 3L, control = ctrl, expected_a_id = c(2, 1, 3)),
    class = "cdid_ipw_owner_order_mismatch_error"
  )
})

# ---------------------------------------------------------------------------
# Audit 4: Dense mask refusal in compact views/queries vs disclosed verification
# ---------------------------------------------------------------------------

test_that("arbitrary dense masks are refused by compact views/queries and verified only via disclosed bounded entrypoint", {
  ns <- require_cdid_namespace()
  ctrl <- audit_owner_ctrl(ns, max_chunk_rows = 3L)
  df <- make_audit_panel_df()

  panel <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl)
  frame <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), control = ctrl)

  # Passing dense masks directly to construct_ipw_nuisance_view, construct_ipw_link_view, or query_ipw_view_chunk is refused
  expect_error(
    ns$construct_ipw_nuisance_view(frame, g = 2, K = c(Inf), s = 1, t = 2, control = ctrl, risk_mask = rep(TRUE, 6L)),
    class = "cdid_ipw_owner_dense_mask_error"
  )
  expect_error(
    ns$construct_ipw_link_view(frame, s = 1, t = 2, g = 2, K = c(Inf), control = ctrl, pair_observed = rep(TRUE, 6L)),
    class = "cdid_ipw_owner_dense_mask_error"
  )

  lk <- ns$construct_ipw_link_view(frame, s = 1, t = 2, g = 2, K = c(Inf), control = ctrl)
  expect_error(
    ns$query_ipw_view_chunk(lk, 1L, 3L, control = ctrl, observed_mask = c(TRUE, FALSE, TRUE)),
    class = "cdid_ipw_owner_dense_mask_error"
  )

  # Disclosed bounded verification succeeds when external dense masks match posting-derived masks
  # Units 1..2 have g = 2 (unit 2 missing at t = 2); units 3..4 have g = 3 (outside risk {2, Inf}); units 5..6 have g = Inf
  valid_risk <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
  valid_obs <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
  valid_ci <- c(1L, 1L, NA_integer_, NA_integer_, 2L, 2L)
  valid_po <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)

  v_res <- ns$verify_ipw_view_dense_masks(
    lk, control = ctrl,
    pair_observed = valid_po,
    risk_mask = valid_risk,
    observed_mask = valid_obs,
    cohort_index = valid_ci
  )
  expect_true(v_res$verified)
  expect_true(v_res$dense_mask_disclosed)
  expect_false(v_res$dense_mask_retained_by_view)
  expect_identical(v_res$caller_dense_mask_entries, 24)
  expect_identical(v_res$endpoints_read, 0L)

  # Flipping a risk-set entry in observed_mask or pair_observed is caught and refused
  bad_obs <- valid_obs
  bad_obs[2L] <- TRUE
  expect_error(
    ns$verify_ipw_view_dense_masks(lk, control = ctrl, observed_mask = bad_obs),
    class = "cdid_ipw_owner_dense_mask_error"
  )
})

# ---------------------------------------------------------------------------
# Audit 5: Risk/read-order mutation discrimination via endpoint instrumentation
# ---------------------------------------------------------------------------

test_that("endpoint-read instrumentation discriminates risk-first reading from pair-first reading and mutated risk sets", {
  ns <- require_cdid_namespace()
  ctrl <- audit_owner_ctrl(ns, max_chunk_rows = 3L)
  df <- make_audit_panel_df()

  panel <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl)
  frame <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = cbind(1, 1:6), control = ctrl)

  # Target g = 2, control K = c(Inf) on pair (s = 1, t = 2).
  # Risk set R = {1, 2, 5, 6}. Pair observers across all cohorts S = {1, 3, 4, 5, 6}.
  # Eligible risk-set observers O = R intersect S = {1, 5, 6}.
  # Risk non-observer R \ O = {2}. Off-risk pair observers S \ R = {3, 4}.
  lk <- ns$construct_ipw_link_view(frame, s = 1, t = 2, g = 2, K = c(Inf), control = ctrl)

  tr_actual <- new.env(parent = emptyenv())
  actual_s <- actual_t <- integer(0L)
  ns$stream_ipw_view_chunks(
    lk, control = ctrl,
    callback = function(ch) {
      actual_s <<- c(actual_s, ch$read_diagnostics$read_obs_s_index)
      actual_t <<- c(actual_t, ch$read_diagnostics$read_obs_t_index)
    },
    read_outcomes = TRUE,
    trace_env = tr_actual
  )

  # Mutated Reader 1 (Pair-First / Filter-After-Read): reads all pair hits S = {1, 3, 4, 5, 6} before filtering by R
  tr_pair_first <- new.env(parent = emptyenv())
  all_pos <- ns$query_ipw_panel_pair_positions(panel, s = 1, t = 2, control = audit_owner_ctrl(ns, 6L))
  ns$ipw_owner_read_y_positions(panel, all_pos$obs_s, all_pos$obs_t, trace_env = tr_pair_first)
  pair_first_s <- tr_pair_first$pos_s_index

  # Mutated Reader 2 (Mutated Risk Set R' = {2, 3} instead of {2, Inf}): reads {1, 3, 4}
  lk_mut_risk <- ns$construct_ipw_link_view(frame, s = 1, t = 2, g = 2, K = c(3), control = ctrl)
  tr_mut_risk <- new.env(parent = emptyenv())
  mut_risk_s <- integer(0L)
  ns$stream_ipw_view_chunks(
    lk_mut_risk, control = ctrl,
    callback = function(ch) mut_risk_s <<- c(mut_risk_s, ch$read_diagnostics$read_obs_s_index),
    read_outcomes = TRUE,
    trace_env = tr_mut_risk
  )

  # Exact discrimination:
  # Unit 1 obs_s = 1; Unit 3 obs_s = 6; Unit 4 obs_s = 9; Unit 5 obs_s = 12; Unit 6 obs_s = 15
  expect_identical(actual_s, c(1L, 12L, 15L))
  expect_identical(actual_t, c(2L, 13L, 16L))
  expect_identical(pair_first_s, c(1L, 6L, 9L, 12L, 15L))
  expect_identical(mut_risk_s, c(1L, 6L, 9L))
  expect_false(identical(actual_s, pair_first_s))
  expect_false(identical(actual_s, mut_risk_s))
})

# ---------------------------------------------------------------------------
# Audit 6: Alias sharing, GC/interrupt propagation, and N/L retained storage ladder
# ---------------------------------------------------------------------------

test_that("owners share long arrays by alias, propagate callback errors/interrupts, and keep view storage invariant across N/L ladders", {
  ns <- require_cdid_namespace()
  ctrl <- audit_owner_ctrl(ns, max_chunk_rows = 16L, max_retained_entries = 500000L, max_matrix_entries = 500000L, max_work_ops = 2000000L, max_output_entries = 500000L)

  df <- make_audit_panel_df()
  panel <- ns$construct_ipw_panel_owner(df, "id", "time", "y", "g", control = ctrl)

  Z_in <- matrix(as.double(1:12), nrow = 6L, ncol = 2L)
  w_in <- as.double(c(1, 1.2, 0.8, 1.1, 0.9, 1.0))
  f1 <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = Z_in, weights = w_in, control = ctrl)
  f2 <- ns$construct_ipw_frame_owner(panel, id = 1:6, Z = Z_in, weights = w_in, control = ctrl)

  # Verify R copy-on-modify alias sharing of long panel tables and canonical Z/weights
  expect_true(f1$accounting$panel_arrays_aliased)
  expect_true(f1$accounting$Z_aliased_from_input)
  expect_true(f1$accounting$w_aliased_from_input)
  expect_identical(f1$panel_ref, panel$ref)
  expect_identical(f2$panel_ref, panel$ref)

  # Verify error and interrupt propagation in stream_ipw_view_chunks and traverse_ipw_panel_pair_positions
  lk <- ns$construct_ipw_link_view(f1, s = 1, t = 3, g = 2, K = c(3, Inf), control = ctrl)
  expect_error(
    ns$stream_ipw_view_chunks(
      lk, control = audit_owner_ctrl(ns, max_chunk_rows = 2L),
      callback = function(ch) {
        if (ch$batch_index == 2L) stop("custom_callback_error")
      }
    ),
    regexp = "custom_callback_error"
  )
  int_cond <- structure(list(message = "user_interrupt"), class = c("interrupt", "condition"))
  int_caught <- tryCatch({
    ns$stream_ipw_view_chunks(
      lk, control = audit_owner_ctrl(ns, max_chunk_rows = 2L),
      callback = function(ch) signalCondition(int_cond)
    )
    FALSE
  }, interrupt = function(e) TRUE)
  expect_true(int_caught)

  # Prespecified N/L ladder test: N in {36, 144}, L in {1, 4, 8}
  build_ladder_instance <- function(n_units, n_links) {
    times <- seq_len(9L)
    cohorts <- c(4, 6, Inf)
    u_g <- rep(cohorts, length.out = n_units)
    rows <- lapply(seq_len(n_units), function(i) {
      data.frame(
        id = rep(i, 9L),
        time = times,
        y = as.double(i * 0.1 + times * 0.25),
        g = rep(u_g[i], 9L)
      )
    })
    df_n <- do.call(rbind, rows)
    Z_n <- cbind(1.0, sin(seq_len(n_units)), cos(seq_len(n_units)))
    p_n <- ns$construct_ipw_panel_owner(df_n, "id", "time", "y", "g", control = ctrl)
    f_n <- ns$construct_ipw_frame_owner(p_n, id = seq_len(n_units), Z = Z_n, control = ctrl)
    nu_n <- ns$construct_ipw_nuisance_view(f_n, g = 4, K = c(6, Inf), s = 1, t = 2, control = ctrl)
    pairs_l <- list(
      c(1, 2), c(2, 3), c(3, 4), c(4, 5),
      c(5, 6), c(6, 7), c(7, 8), c(8, 9)
    )[seq_len(n_links)]
    links <- lapply(pairs_l, function(st) {
      ns$construct_ipw_link_view(
        f_n, s = st[1L], t = st[2L], nuisance_view = nu_n,
        loading_beta_candidate = c(0.1, -0.2, 0.05), control = ctrl
      )
    })
    gc(verbose = FALSE)
    list(panel = p_n, frame = f_n, nuisance = nu_n, links = links)
  }

  inst_36_1 <- build_ladder_instance(36L, 1L)
  inst_36_8 <- build_ladder_instance(36L, 8L)
  inst_144_1 <- build_ladder_instance(144L, 1L)
  inst_144_8 <- build_ladder_instance(144L, 8L)

  # 1. Retained entries and object.size of a single nuisance_view or link_view are O(p + m),
  #    strictly invariant to N (36 vs 144)!
  expect_identical(
    inst_36_1$nuisance$preflight$total_retained_entries,
    inst_144_1$nuisance$preflight$total_retained_entries
  )
  expect_identical(
    inst_36_1$links[[1L]]$preflight$total_retained_entries,
    inst_144_1$links[[1L]]$preflight$total_retained_entries
  )
  expect_identical(
    as.numeric(object.size(inst_36_1$links[[1L]][c("pair", "cohorts", "support", "loading_beta_candidate", "accounting")])),
    as.numeric(object.size(inst_144_1$links[[1L]][c("pair", "cohorts", "support", "loading_beta_candidate", "accounting")]))
  )

  # 2. Marginal retained entries when adding 7 links (L = 1 -> L = 8) are identical at N = 36 and N = 144
  #    (proving zero O(NL) per-link retained arrays in the new views!)
  marginal_links_36 <- sum(vapply(inst_36_8$links, function(x) x$preflight$total_retained_entries, numeric(1L))) -
    inst_36_1$links[[1L]]$preflight$total_retained_entries
  marginal_links_144 <- sum(vapply(inst_144_8$links, function(x) x$preflight$total_retained_entries, numeric(1L))) -
    inst_144_1$links[[1L]]$preflight$total_retained_entries
  expect_identical(marginal_links_36, marginal_links_144)
})
