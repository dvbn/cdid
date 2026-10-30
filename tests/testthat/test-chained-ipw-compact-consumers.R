# Independent testthat suite for private CP02 compact-view numerical consumers
# (`R/chained_logit_certificate.R` and `R/chained_ipw_enclosure.R`).

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

# ---------------------------------------------------------------------------
# Deterministic Fixture-to-Panel Adapter (Zero modification of frozen fixtures)
# ---------------------------------------------------------------------------
# Builds finite panel records with identical A-membership, canonical unit order,
# pair-observer sets S_l, and in-domain stored endpoints/design/weights:
# - For in-risk pair observers (u in O_l = R intersect S_l): exact stored
#   y_left[u] and y_right[u] at dates (s_l, t_l).
# - For off-risk pair observers (u in S_l \ R): if the dense fixture used
#   nonfinite sentinels (NA/NaN/Inf/-Inf) or off-domain placeholders, prescribes
#   finite extreme off-risk values (-off_domain_mag, +off_domain_mag) at (s_l, t_l)
#   whose difference would overflow IEEE binary64 if ever read/subtracted,
#   preserving u in S_l and outside_risk_pair_hit_count while proving 0 reads.
# - For pair non-observers (u not in S_l): places auxiliary finite observations
#   at auxiliary dates (aux_s, aux_t) = (100.0, 200.0) so u has >= 2 panel
#   observations and belongs to A in the exact same canonical order while
#   remaining unobserved on (s_l, t_l).

adapt_dense_fixture_to_panel_records <- function(
  id,
  Z,
  w,
  G,
  cohorts,
  links,
  aux_s = 100.0,
  aux_t = 200.0,
  off_domain_mag = 1e308
) {
  n <- nrow(Z)
  stopifnot(length(id) == n, length(w) == n, length(G) == n)
  in_risk <- !is.na(match(as.double(G), as.double(cohorts)))

  rows <- list()
  in_domain_records <- list()
  substituted_off_risk_records <- list()
  auxiliary_nonobserver_records <- list()

  for (i in seq_len(n)) {
    uid <- id[i]
    ug <- as.double(G[i])
    obs_dates_for_u <- 0L

    for (l_idx in seq_along(links)) {
      lk <- links[[l_idx]]
      po_i <- isTRUE(lk$pair_observed[i])
      if (po_i) {
        s_val <- as.double(lk$s)
        t_val <- as.double(lk$t)
        yl_raw <- if (!is.null(lk$y_left)) as.double(lk$y_left[i]) else 0.0
        yr_raw <- if (!is.null(lk$y_right)) as.double(lk$y_right[i]) else 0.0

        if (in_risk[i]) {
          # In-domain observer u in O_l = R intersect S_l: preserve exact stored values
          stopifnot(is.finite(yl_raw), is.finite(yr_raw))
          yl_use <- yl_raw
          yr_use <- yr_raw
          in_domain_records[[length(in_domain_records) + 1L]] <- data.frame(
            a_index = i, link_index = l_idx, s = s_val, t = t_val,
            y_left = yl_use, y_right = yr_use, stringsAsFactors = FALSE
          )
        } else {
          # Off-risk pair observer u in S_l \ R: prescribe finite extreme endpoints
          # if sentinel or off-domain test requests it
          use_sub <- !is.finite(yl_raw) || !is.finite(yr_raw) || isTRUE(lk$substitute_off_risk)
          yl_use <- if (use_sub) -as.double(off_domain_mag) else yl_raw
          yr_use <- if (use_sub) +as.double(off_domain_mag) else yr_raw
          substituted_off_risk_records[[length(substituted_off_risk_records) + 1L]] <- data.frame(
            a_index = i, link_index = l_idx, s = s_val, t = t_val,
            raw_y_left = yl_raw, raw_y_right = yr_raw,
            prescribed_y_left = yl_use, prescribed_y_right = yr_use,
            stringsAsFactors = FALSE
          )
        }

        rows[[length(rows) + 1L]] <- data.frame(
          id = uid, time = s_val, y = yl_use, g = ug, stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1L]] <- data.frame(
          id = uid, time = t_val, y = yr_use, g = ug, stringsAsFactors = FALSE
        )
        obs_dates_for_u <- obs_dates_for_u + 2L
      }
    }

    if (obs_dates_for_u < 2L) {
      # Non-observer across all target links: add auxiliary pair at (aux_s, aux_t)
      # so unit i belongs to A with >= 2 panel observations
      rows[[length(rows) + 1L]] <- data.frame(
        id = uid, time = as.double(aux_s), y = -as.double(off_domain_mag), g = ug,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        id = uid, time = as.double(aux_t), y = +as.double(off_domain_mag), g = ug,
        stringsAsFactors = FALSE
      )
      auxiliary_nonobserver_records[[length(auxiliary_nonobserver_records) + 1L]] <- data.frame(
        a_index = i, aux_s = as.double(aux_s), aux_t = as.double(aux_t),
        y_aux_s = -as.double(off_domain_mag), y_aux_t = +as.double(off_domain_mag),
        stringsAsFactors = FALSE
      )
    }
  }

  panel_df <- do.call(rbind, rows)
  rownames(panel_df) <- NULL

  list(
    panel_df = panel_df,
    id = id,
    Z = Z,
    w = as.double(w),
    G = as.double(G),
    cohorts = as.double(cohorts),
    in_risk = in_risk,
    in_domain_records = if (length(in_domain_records) > 0L) do.call(rbind, in_domain_records) else data.frame(),
    substituted_off_risk_records = if (length(substituted_off_risk_records) > 0L) do.call(rbind, substituted_off_risk_records) else data.frame(),
    auxiliary_nonobserver_records = if (length(auxiliary_nonobserver_records) > 0L) do.call(rbind, auxiliary_nonobserver_records) else data.frame()
  )
}

# Test-owned instrumentation keeps the full read-set assertion outside the core.
# Production diagnostics retain only the latest bounded query's positions.
record_endpoint_calls <- function(ns, run) {
  original <- get("ipw_owner_read_y_positions", envir = ns, inherits = FALSE)
  all_s <- all_t <- integer()
  instrumented <- function(panel_owner, pos_s, pos_t, trace_env = NULL) {
    all_s <<- c(all_s, pos_s); all_t <<- c(all_t, pos_t)
    original(panel_owner, pos_s, pos_t, trace_env)
  }
  unlockBinding("ipw_owner_read_y_positions", ns)
  assign("ipw_owner_read_y_positions", instrumented, envir = ns)
  lockBinding("ipw_owner_read_y_positions", ns)
  on.exit({
    unlockBinding("ipw_owner_read_y_positions", ns)
    assign("ipw_owner_read_y_positions", original, envir = ns)
    lockBinding("ipw_owner_read_y_positions", ns)
  }, add = TRUE)
  value <- run()
  list(value = value, pos_s = all_s, pos_t = all_t)
}

# Helper to verify that no retained N-length or N x p arrays exist on a compact object
assert_zero_n_retention <- function(obj, n_a, allowed_shared_envs = list()) {
  seen_envs <- list()
  env_in_list <- function(e, lst) {
    for (item in lst) if (identical(e, item)) return(TRUE)
    FALSE
  }
  walk <- function(x, path) {
    if (is.null(x)) return(invisible(NULL))
    if (is.environment(x)) {
      if (env_in_list(x, allowed_shared_envs) || env_in_list(x, seen_envs) ||
          identical(x, baseenv()) || identical(x, emptyenv()) ||
          identical(x, .GlobalEnv) || isNamespace(x)) {
        return(invisible(NULL))
      }
      seen_envs[[length(seen_envs) + 1L]] <<- x
      for (nm in ls(x, all.names = TRUE)) {
        walk(get(nm, envir = x, inherits = FALSE), paste0(path, "$", nm))
      }
      return(invisible(NULL))
    }
    if (is.function(x)) {
      walk(environment(x), paste0("env(", path, ")"))
      return(invisible(NULL))
    }
    if (is.atomic(x)) {
      if (!is.null(dim(x))) {
        if (any(dim(x) >= n_a)) {
          stop(sprintf("Retained array `%s` has dimension %s >= n_a (%d).",
                       path, paste(dim(x), collapse = "x"), n_a), call. = FALSE)
        }
      } else if (length(x) >= n_a) {
        stop(sprintf("Retained atomic vector `%s` has length %d >= n_a (%d).",
                     path, length(x), n_a), call. = FALSE)
      }
      return(invisible(NULL))
    }
    if (is.list(x)) {
      nms <- names(x)
      for (k in seq_along(x)) {
        lbl <- if (!is.null(nms) && nzchar(nms[k])) nms[k] else as.character(k)
        walk(x[[k]], paste0(path, "$", lbl))
      }
    }
    invisible(NULL)
  }
  walk(obj, "root")
  TRUE
}

# ---------------------------------------------------------------------------
# Test 1: Locked namespace bindings and unexported status of compact consumers
# ---------------------------------------------------------------------------

test_that("compact-view consumer entrypoints exist in locked cdid namespace, are unexported, and refuse unsealed contexts", {
  ns <- require_cdid_namespace()
  expect_true(isNamespace(ns))
  expect_true(environmentIsLocked(ns))

  compact_fns <- c(
    "preflight_logit_certificate_view",
    "certify_chained_logit_root_view",
    "chained_logit_certificate_view",
    "certify_weighted_logit_root_view",
    "preflight_ipw_nuisance_problem_view",
    "preflight_ipw_link_enclosure_view",
    "preflight_ipw_link_if_view_chunk",
    "construct_certified_ipw_nuisance_problem_view",
    "enclose_ipw_link_view",
    "enclose_ipw_shared_xi_view_chunk",
    "apply_ipw_shared_xi_to_link_view",
    "enclose_ipw_link_if_view_chunk",
    "stream_ipw_link_if_view_enclosure"
  )
  exports <- getNamespaceExports("cdid")
  for (fn in compact_fns) {
    expect_true(exists(fn, envir = ns, mode = "function", inherits = FALSE), info = fn)
    expect_true(bindingIsLocked(fn, ns), info = fn)
    expect_identical(environment(get(fn, envir = ns, inherits = FALSE)), ns, info = fn)
    expect_false(fn %in% exports, info = fn)
  }

  # Internal shared workers must refuse unsealed or forged context objects
  expect_error(
    ns$logit_cert_run_with_context(list(mode = "dense")),
    class = "cdid_logit_certificate_integrity_error"
  )
  forged_ctx <- structure(list(mode = "dense", ref = new.env(parent = emptyenv())), class = "cdid_logit_access_context")
  lockEnvironment(forged_ctx$ref, bindings = TRUE)
  expect_error(
    ns$logit_cert_run_with_context(forged_ctx),
    class = "cdid_logit_certificate_integrity_error"
  )
  expect_error(
    ns$ipw_enclosure_link_run_with_context(list(mode = "compact_view")),
    class = "cdid_ipw_enclosure_integrity_error"
  )
  expect_error(
    ns$ipw_enclosure_shared_xi_run_with_context(list(mode = "compact_view")),
    class = "cdid_ipw_enclosure_integrity_error"
  )
  expect_error(
    ns$ipw_enclosure_if_chunk_run_with_context(list(mode = "compact_view")),
    class = "cdid_ipw_enclosure_integrity_error"
  )
})

# ---------------------------------------------------------------------------
# Test 2: Frozen fixture SHA256 verification, deterministic adapter invariants,
# and bitwise dense vs compact equality across chunk sizes {1, 2, 3, 7, 64, N}
# ---------------------------------------------------------------------------

test_that("fixture-to-panel adapter preserves in-domain data, substitutes off-domain sentinels with 0 reads, and matches dense bitwise across chunk sizes", {
  ns <- require_cdid_namespace()

  logit_fx_path <- testthat::test_path("fixtures", "cp02_logit_certificate_fixtures.R")
  ipw_fx_path <- testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R")
  expect_identical(
    digest::digest(file = logit_fx_path, algo = "sha256"),
    "008ac83235b01f77793dc286d6171f8b119eec13776546a97dc1913f1134d759"
  )
  expect_identical(
    digest::digest(file = ipw_fx_path, algo = "sha256"),
    "2a418bee8c504ae443683f7476f7d6ee59eedc45cd87550a1865f8384e3d2f0c"
  )

  logit_env <- new.env(parent = ns)
  sys.source(logit_fx_path, envir = logit_env)
  ipw_ref <- dget(ipw_fx_path)

  # 2A. Test all raw logit certificate fixtures across chunk sizes {1, 2, 3, 7, 64, N}
  raw_wr_fx <- logit_env$make_raw_vs_cached_wr_fixture()
  raw_wr_fx$r <- 0.01
  logit_fixtures <- list(
    beta0_zero = logit_env$make_exact_beta0_zero_fixture(),
    prop_16 = logit_env$make_propagation_16unit_fixture(),
    ov_8 = logit_env$make_finite_logit_ov_fixture(),
    unsat_3cell = logit_env$make_unsaturated_3cell_fixture(),
    raw_vs_round = raw_wr_fx
  )

  for (fx_nm in names(logit_fixtures)) {
    fx <- logit_fixtures[[fx_nm]]
    n_fx <- nrow(fx$Z)
    id_fx <- seq_len(n_fx)
    adapted <- adapt_dense_fixture_to_panel_records(
      id = id_fx, Z = fx$Z, w = fx$w, G = fx$G, cohorts = fx$cohorts,
      links = list(list(
        s = 10.0, t = 20.0, pair_observed = fx$pair_observed,
        y_left = rep(0.0, n_fx), y_right = rep(0.0, n_fx),
        substitute_off_risk = TRUE
      ))
    )
    own_ctrl <- ns$ipw_owner_control(max(1L, n_fx), 100000L, 100000L, 500000L, 100000L)
    panel <- ns$construct_ipw_panel_owner(adapted$panel_df, "id", "time", "y", "g", control = own_ctrl)
    frame <- ns$construct_ipw_frame_owner(panel, id = id_fx, Z = fx$Z, weights = fx$w, control = own_ctrl)
    nu_view <- ns$construct_ipw_nuisance_view(
      frame, g = fx$cohorts[1L], K = fx$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl
    )

    # Verify adapter preserved exact A-membership, order, design, weights, and pair counts
    expect_identical(panel$a_id, id_fx)
    expect_identical(frame$Z_a, fx$Z)
    expect_identical(frame$w_a, as.double(fx$w))
    expect_identical(nu_view$support$n_observed_pair_all, sum(fx$pair_observed))

    for (b_sz in unique(c(1L, 2L, 3L, 7L, 64L, n_fx))) {
      cert_ctrl <- ns$logit_certificate_control(18L, 700.0, 11L, b_sz, 100000L, 1000000L)
      tr_env <- new.env(parent = emptyenv())
      r_ctrl <- ns$ipw_reader_control(b_sz, 100000L, 100000L, 1000000L, 100000L)
      attr(r_ctrl, "trace_env") <- tr_env

      c_dense <- ns$certify_chained_logit_root(
        fx$Z, fx$w, fx$G, fx$pair_observed, fx$cohorts, fx$beta0, fx$F, fx$r, cert_ctrl
      )
      c_view <- ns$certify_chained_logit_root_view(
        nu_view, fx$beta0, fx$F, fx$r, cert_ctrl, r_ctrl
      )

      expect_identical(c_view$certified, c_dense$certified, info = paste(fx_nm, b_sz))
      expect_identical(c_view$status, c_dense$status, info = paste(fx_nm, b_sz))
      expect_identical(c_view$enclosures, c_dense$enclosures, info = paste(fx_nm, b_sz))
      expect_identical(c_view$krawczyk, c_dense$krawczyk, info = paste(fx_nm, b_sz))
      expect_identical(c_view$dimensions, c_dense$dimensions, info = paste(fx_nm, b_sz))
      expect_identical(c_view$diagnostics$max_abs_eta_0, c_dense$diagnostics$max_abs_eta_0)
      expect_identical(c_view$diagnostics$max_abs_eta_ball, c_dense$diagnostics$max_abs_eta_ball)
      expect_identical(c_view$provenance$storage_mode, "compact_view")
      expect_identical(c_view$read_diagnostics$endpoints_read, 0L)
      expect_identical(tr_env$calls, 0L)
    }
  }

  # 2B. Test IPW link enclosure fixtures (f1 exact rational 2-link, f2 nonzero beta rootball, f3 real endpoint subtraction)
  f1 <- ipw_ref$fixture1_exact_rational
  yL1_sentinels <- f1$yL_1_float
  yR1_sentinels <- f1$yR_1_float
  # Off-observed nonfinite sentinels as in dense fixture test
  yL1_sentinels[c(3L, 6L, 8L, 9L, 10L)] <- c(NA_real_, NaN, Inf, -Inf, NA_real_)
  yR1_sentinels[c(3L, 6L, 8L, 9L, 10L)] <- c(Inf, -Inf, NA_real_, NaN, Inf)

  adapted_f1 <- adapt_dense_fixture_to_panel_records(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G, cohorts = f1$cohorts,
    links = list(
      list(s = 10.0, t = 20.0, pair_observed = f1$pair_observed, y_left = yL1_sentinels, y_right = yR1_sentinels, substitute_off_risk = TRUE),
      list(s = 30.0, t = 40.0, pair_observed = f1$pair_observed, y_left = f1$yL_2_float, y_right = f1$yR_2_float, substitute_off_risk = TRUE)
    )
  )
  # Verify off-risk observer u=9 had nonfinite sentinel substituted by +-1e308
  expect_true(9L %in% adapted_f1$substituted_off_risk_records$a_index)
  # Verify non-observers u=3, 6, 8, 10 received auxiliary observations at (100, 200)
  expect_identical(adapted_f1$auxiliary_nonobserver_records$a_index, c(3L, 6L, 8L, 10L))

  own_ctrl10 <- ns$ipw_owner_control(10L, 100000L, 100000L, 1000000L, 100000L)
  panel_f1 <- ns$construct_ipw_panel_owner(adapted_f1$panel_df, "id", "time", "y", "g", control = own_ctrl10, owner_id = "panel_f1")
  frame_f1 <- ns$construct_ipw_frame_owner(panel_f1, id = f1$id, Z = f1$Z_float, weights = f1$w_float, control = own_ctrl10, frame_id = "frame_f1")
  nu_f1 <- ns$construct_ipw_nuisance_view(frame_f1, g = f1$cohorts[1L], K = f1$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl10, nuisance_id = "nu_f1")
  lk1_f1 <- ns$construct_ipw_link_view(frame_f1, s = 10.0, t = 20.0, nuisance_view = nu_f1, control = own_ctrl10, link_id = "link1")
  lk2_f1 <- ns$construct_ipw_link_view(frame_f1, s = 30.0, t = 40.0, nuisance_view = nu_f1, control = own_ctrl10, link_id = "link2")

  # Reference dense objects
  cert_ctrl_ref <- ns$logit_certificate_control(18L, 40.0, 16L, 4L, 50000L, 1000000L)
  prob_ctrl_ref <- ns$ipw_nuisance_problem_control(4L, 50000L, 50000L, 1000000L)
  link_ctrl_ref <- ns$ipw_link_enclosure_control(10L, 50000L, 50000L, 1000000L, 50000L)

  prob_dense_f1 <- ns$construct_certified_ipw_nuisance_problem(
    f1$id, f1$Z_float, f1$w_float, f1$G, f1$pair_observed, f1$cohorts,
    f1$beta_star_float, f1$F_float, 0.05, cert_ctrl_ref, prob_ctrl_ref
  )
  enc1_dense_f1 <- ns$enclose_ipw_link(
    prob_dense_f1, f1$id, yL1_sentinels, yR1_sentinels,
    f1$link1$lambda_beta_float, link_ctrl_ref, link_id = "link1"
  )
  enc2_dense_f1 <- ns$enclose_ipw_link(
    prob_dense_f1, f1$id, f1$yL_2_float, f1$yR_2_float,
    f1$link2$lambda_beta_float, link_ctrl_ref, link_id = "link2"
  )
  R_H_f1 <- chol(f1$H_float)
  xi_dense_f1 <- ns$enclose_ipw_shared_xi_chunk(prob_dense_f1, 1L, 10L, link_ctrl_ref, H_chol_candidate = R_H_f1)
  ch1_dense_f1 <- ns$enclose_ipw_link_if_chunk(
    enc1_dense_f1, 1L, 10L, link_ctrl_ref,
    phi_candidate = f1$link1$phi_float, xi_beta_candidate = f1$link1$xi_beta_float,
    include_components = TRUE, include_shared_xi = TRUE
  )
  st1_dense_f1 <- ns$stream_ipw_link_if_enclosure(
    enc1_dense_f1, link_ctrl_ref, callback = function(b) NULL,
    phi_candidate = f1$link1$phi_float,
    H_chol_candidate = R_H_f1, include_components = TRUE, include_shared_xi = TRUE
  )

  # Expected observation indices in panel_f1 for link1 (dates 10, 20) and link2 (dates 30, 40) on O = {1, 2, 4, 5, 7}
  obs_tbl <- panel_f1$observations
  allowed_obs_lk1 <- which(obs_tbl$a_index %in% c(1L, 2L, 4L, 5L, 7L) & obs_tbl$time %in% c(10.0, 20.0))
  allowed_obs_lk2 <- which(obs_tbl$a_index %in% c(1L, 2L, 4L, 5L, 7L) & obs_tbl$time %in% c(30.0, 40.0))
  forbidden_obs <- setdiff(seq_len(nrow(obs_tbl)), c(allowed_obs_lk1, allowed_obs_lk2))
  expect_true(length(forbidden_obs) >= 10L) # includes off-risk u=9 at (10,20),(30,40) and non-observers at (100,200)

  for (b_sz in c(1L, 2L, 3L, 7L, 10L, 64L)) {
    c_ctrl <- ns$logit_certificate_control(18L, 40.0, 16L, b_sz, 50000L, 1000000L)
    p_ctrl <- ns$ipw_nuisance_problem_control(b_sz, 50000L, 50000L, 1000000L)
    l_ctrl <- ns$ipw_link_enclosure_control(b_sz, 50000L, 50000L, 1000000L, 50000L)
    tr_nuis <- new.env(parent = emptyenv())
    r_nuis <- ns$ipw_reader_control(b_sz, 50000L, 50000L, 1000000L, 50000L)
    attr(r_nuis, "trace_env") <- tr_nuis

    prob_v <- ns$construct_certified_ipw_nuisance_problem_view(
      nu_f1, f1$beta_star_float, f1$F_float, 0.05, c_ctrl, p_ctrl, r_nuis
    )
    ns$ipw_validate_nuisance_problem_mode(prob_v, require_certified = TRUE, expected_storage_mode = "compact_view")
    expect_identical(tr_nuis$calls, 0L)
    expect_identical(prob_v$read_diagnostics$endpoints_read, 0L)
    expect_identical(prob_v$dimensions, prob_dense_f1$dimensions)
    expect_identical(prob_v$nuisance_enclosures, prob_dense_f1$nuisance_enclosures)
    expect_identical(prob_v$certificate$enclosures, prob_dense_f1$certificate$enclosures)

    # Zero retained N-arrays on prob_v
    expect_null(prob_v$raw$Z)
    expect_null(prob_v$raw$id)
    expect_null(prob_v$raw$w)
    expect_null(prob_v$raw$G)
    expect_null(prob_v$raw$pair_observed)
    expect_null(prob_v$masks$in_risk)
    expect_null(prob_v$masks$in_observed)
    expect_true(assert_zero_n_retention(
      prob_v, n_a = 10L,
      allowed_shared_envs = list(panel_f1$ref, frame_f1$ref, nu_f1$ref)
    ))

    # Enclose link 1 and link 2 across chunk size b_sz
    tr_lk1 <- new.env(parent = emptyenv())
    r_lk1 <- ns$ipw_reader_control(b_sz, 50000L, 50000L, 1000000L, 50000L)
    attr(r_lk1, "trace_env") <- tr_lk1
    capture1 <- record_endpoint_calls(ns, function() ns$enclose_ipw_link_view(prob_v, lk1_f1, f1$link1$lambda_beta_float, l_ctrl, r_lk1))
    enc1_v <- capture1$value
    ns$ipw_validate_link_enclosure_mode(enc1_v, require_certified = TRUE, expected_storage_mode = "compact_view")

    expect_identical(enc1_v$enclosures_0, enc1_dense_f1$enclosures_0)
    expect_identical(enc1_v$enclosures_star, enc1_dense_f1$enclosures_star)
    expect_identical(enc1_v$root_perturbation_bounds, enc1_dense_f1$root_perturbation_bounds)
    expect_identical(enc1_v$endpoint_bounds$omega_target, enc1_dense_f1$endpoint_bounds$omega_target)
    expect_identical(enc1_v$endpoint_bounds$omega_control, enc1_dense_f1$endpoint_bounds$omega_control)
    expect_identical(enc1_v$endpoint_bounds$omega, enc1_dense_f1$endpoint_bounds$omega)

    # Verify only in-risk observed endpoints for link1 were read (2 passes * 5 units * 2 endpoints = 20)
    expect_identical(enc1_v$read_diagnostics$endpoints_read, 20L)
    read_idx_lk1 <- unique(c(capture1$pos_s, capture1$pos_t))
    expect_lte(length(tr_lk1$obs_s_read), b_sz)
    expect_lte(length(tr_lk1$obs_t_read), b_sz)
    expect_identical(sort(read_idx_lk1), sort(allowed_obs_lk1))
    expect_length(intersect(read_idx_lk1, forbidden_obs), 0L)

    # Zero retained endpoint copies or N-position maps on enc1_v
    expect_null(enc1_v$endpoints$obs_a_index)
    expect_null(enc1_v$endpoints$a_to_obs_pos)
    expect_null(enc1_v$endpoints$y_left_obs)
    expect_null(enc1_v$endpoints$y_right_obs)
    expect_true(assert_zero_n_retention(
      enc1_v, n_a = 10L,
      allowed_shared_envs = list(panel_f1$ref, frame_f1$ref, nu_f1$ref, lk1_f1$ref, prob_v$ref)
    ))

    tr_lk2 <- new.env(parent = emptyenv())
    r_lk2 <- ns$ipw_reader_control(b_sz, 50000L, 50000L, 1000000L, 50000L)
    attr(r_lk2, "trace_env") <- tr_lk2
    capture2 <- record_endpoint_calls(ns, function() ns$enclose_ipw_link_view(prob_v, lk2_f1, f1$link2$lambda_beta_float, l_ctrl, r_lk2))
    enc2_v <- capture2$value
    expect_identical(enc2_v$enclosures_0, enc2_dense_f1$enclosures_0)
    expect_identical(enc2_v$enclosures_star, enc2_dense_f1$enclosures_star)
    read_idx_lk2 <- unique(c(capture2$pos_s, capture2$pos_t))
    expect_lte(length(tr_lk2$obs_s_read), b_sz)
    expect_lte(length(tr_lk2$obs_t_read), b_sz)
    expect_identical(sort(read_idx_lk2), sort(allowed_obs_lk2))
    expect_length(intersect(read_idx_lk2, forbidden_obs), 0L)

    # Streamed IF enclosure across chunk size b_sz matches full-chunk dense reference row-by-row and sum-by-sum
    tr_st <- new.env(parent = emptyenv())
    r_st <- ns$ipw_reader_control(b_sz, 50000L, 50000L, 1000000L, 50000L)
    attr(r_st, "trace_env") <- tr_st
    streamed_phi_star <- matrix(NA_real_, nrow = 10L, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    streamed_treated <- matrix(NA_real_, nrow = 10L, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    streamed_control <- matrix(NA_real_, nrow = 10L, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    streamed_score <- matrix(NA_real_, nrow = 10L, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    streamed_rate <- matrix(NA_real_, nrow = 10L, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))
    streamed_phi_shared <- matrix(NA_real_, nrow = 10L, ncol = 2L, dimnames = list(NULL, c("lo", "hi")))

    st_v <- ns$stream_ipw_link_if_view_enclosure(
      link_enclosure = enc1_v, control = l_ctrl,
      callback = function(batch) {
        idx <- batch$a_index
        streamed_phi_star[idx, ] <<- batch$phi_star_intersected
        streamed_treated[idx, ] <<- batch$components_star$treated
        streamed_control[idx, ] <<- batch$components_star$control
        streamed_score[idx, ] <<- batch$components_star$score
        streamed_rate[idx, ] <<- batch$components_star$rate
        streamed_phi_shared[idx, ] <<- batch$shared_xi$phi_shared_star
      },
      reader_control = r_st,
      phi_candidate = f1$link1$phi_float,
      xi_beta_candidate = f1$link1$xi_beta_float,
      include_components = TRUE,
      include_shared_xi = TRUE
    )
    expect_identical(streamed_phi_star, ch1_dense_f1$phi_star_intersected)
    expect_identical(streamed_treated, ch1_dense_f1$components_star$treated)
    expect_identical(streamed_control, ch1_dense_f1$components_star$control)
    expect_identical(streamed_score, ch1_dense_f1$components_star$score)
    expect_identical(streamed_rate, ch1_dense_f1$components_star$rate)
    expect_identical(streamed_phi_shared, ch1_dense_f1$shared_xi$phi_shared_star)
    expect_identical(st_v$phi_0_sum, st1_dense_f1$phi_0_sum)
    expect_identical(st_v$phi_star_sum, st1_dense_f1$phi_star_sum)
    expect_identical(st_v$phi_star_sum_contains_zero, TRUE)
    expect_identical(st_v$max_E_root, st1_dense_f1$max_E_root)
    expect_identical(st_v$read_diagnostics$endpoints_read, 10L)
    expect_length(intersect(unique(c(tr_st$obs_s_read, tr_st$obs_t_read)), forbidden_obs), 0L)
  }

  # 2C. Nonzero beta root-ball fixture (fixture2) and exact real endpoint subtraction fixture (fixture3)
  f2 <- ipw_ref$fixture2_nonzero_beta_rootball
  z_ref <- f2$at_beta0_perturbed
  s_ref <- f2$at_beta_star
  adapted_f2 <- adapt_dense_fixture_to_panel_records(
    id = f2$id, Z = f2$Z_float, w = f2$w_float, G = f2$G, cohorts = f2$cohorts,
    links = list(list(
      s = 10.0, t = 20.0, pair_observed = f2$pair_observed,
      y_left = f2$yL_float, y_right = f2$yR_float, substitute_off_risk = TRUE
    ))
  )
  panel_f2 <- ns$construct_ipw_panel_owner(adapted_f2$panel_df, "id", "time", "y", "g", control = own_ctrl10)
  frame_f2 <- ns$construct_ipw_frame_owner(panel_f2, id = f2$id, Z = f2$Z_float, weights = f2$w_float, control = own_ctrl10)
  nu_f2 <- ns$construct_ipw_nuisance_view(frame_f2, g = f2$cohorts[1L], K = f2$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl10)
  lk_f2 <- ns$construct_ipw_link_view(frame_f2, s = 10.0, t = 20.0, nuisance_view = nu_f2, control = own_ctrl10)

  r_ctrl10 <- ns$ipw_reader_control(10L, 50000L, 50000L, 1000000L, 50000L)
  prob_d2 <- ns$construct_certified_ipw_nuisance_problem(
    f2$id, f2$Z_float, f2$w_float, f2$G, f2$pair_observed, f2$cohorts,
    f2$beta0_perturbed_float, f2$F_float, f2$r_ball, cert_ctrl_ref, prob_ctrl_ref
  )
  prob_v2 <- ns$construct_certified_ipw_nuisance_problem_view(
    nu_f2, f2$beta0_perturbed_float, f2$F_float, f2$r_ball, cert_ctrl_ref, prob_ctrl_ref, r_ctrl10
  )
  expect_identical(prob_v2$nuisance_enclosures, prob_d2$nuisance_enclosures)

  cands_f2 <- list(
    treated_mean = z_ref$T, control_mean = z_ref$C, mu = z_ref$mu,
    v = z_ref$v, gamma = z_ref$gamma, L = z_ref$L,
    loading_beta = z_ref$lambda_beta, loading_rate = z_ref$lambda_r
  )
  enc_d2 <- ns$enclose_ipw_link(prob_d2, f2$id, f2$yL_float, f2$yR_float, z_ref$lambda_beta, link_ctrl_ref, candidates = cands_f2)
  enc_v2 <- ns$enclose_ipw_link_view(prob_v2, lk_f2, z_ref$lambda_beta, link_ctrl_ref, r_ctrl10, candidates = cands_f2)
  expect_identical(enc_v2$enclosures_0, enc_d2$enclosures_0)
  expect_identical(enc_v2$enclosures_star, enc_d2$enclosures_star)
  expect_identical(enc_v2$root_perturbation_bounds, enc_d2$root_perturbation_bounds)
  expect_identical(enc_v2$separated_errors, enc_d2$separated_errors)
  expect_true(all(enc_v2$enclosures_0$lambda_r[, 1L] <= z_ref$lambda_r & z_ref$lambda_r <= enc_v2$enclosures_0$lambda_r[, 2L]))
  expect_true(all(enc_v2$enclosures_star$lambda_r[, 1L] <= s_ref$lambda_r & s_ref$lambda_r <= enc_v2$enclosures_star$lambda_r[, 2L]))

  ch_d2 <- ns$enclose_ipw_link_if_chunk(
    enc_d2, 1L, 10L, link_ctrl_ref, phi_candidate = z_ref$phi,
    H_chol_candidate = chol(z_ref$H), include_components = TRUE, include_shared_xi = TRUE
  )
  ch_v2 <- ns$enclose_ipw_link_if_view_chunk(
    enc_v2, 1L, 10L, link_ctrl_ref, r_ctrl10, phi_candidate = z_ref$phi,
    H_chol_candidate = chol(z_ref$H), include_components = TRUE, include_shared_xi = TRUE
  )
  expect_identical(ch_v2$phi_0, ch_d2$phi_0)
  expect_identical(ch_v2$phi_star, ch_d2$phi_star)
  expect_identical(ch_v2$phi_star_via_root, ch_d2$phi_star_via_root)
  expect_identical(ch_v2$phi_star_intersected, ch_d2$phi_star_intersected)
  expect_identical(ch_v2$E_root, ch_d2$E_root)
  expect_identical(ch_v2$components_0, ch_d2$components_0)
  expect_identical(ch_v2$components_star, ch_d2$components_star)
  expect_identical(ch_v2$components_E_root, ch_d2$components_E_root)
  expect_identical(ch_v2$separated_row_errors, ch_d2$separated_row_errors)
  expect_identical(ch_v2$shared_xi$phi_shared_star, ch_d2$shared_xi$phi_shared_star)
  expect_true(all(ch_v2$phi_0[, 1L] <= z_ref$phi & z_ref$phi <= ch_v2$phi_0[, 2L]))
  expect_true(all(ch_v2$phi_star_intersected[, 1L] <= s_ref$phi & s_ref$phi <= ch_v2$phi_star_intersected[, 2L]))

  # Fixture 3: exact real endpoint subtraction in compact view
  f3 <- ipw_ref$fixture3_real_endpoint_subtraction
  yL_sub <- rep(0.0, 10L); yR_sub <- rep(0.0, 10L)
  yL_sub[c(1L, 2L)] <- f3$yL_target_float
  yR_sub[c(1L, 2L)] <- f3$yR_target_float
  yR_sub[c(4L, 5L, 7L)] <- f3$yR_control_float
  yL_sub[4L] <- f3$yL_control_u4_float
  yL_sub[c(5L, 7L)] <- f3$yL_control_u5_u7_float

  adapted_f3 <- adapt_dense_fixture_to_panel_records(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G, cohorts = f1$cohorts,
    links = list(list(
      s = 10.0, t = 20.0, pair_observed = f1$pair_observed,
      y_left = yL_sub, y_right = yR_sub, substitute_off_risk = TRUE
    ))
  )
  panel_f3 <- ns$construct_ipw_panel_owner(adapted_f3$panel_df, "id", "time", "y", "g", control = own_ctrl10)
  frame_f3 <- ns$construct_ipw_frame_owner(panel_f3, id = f1$id, Z = f1$Z_float, weights = f1$w_float, control = own_ctrl10)
  nu_f3 <- ns$construct_ipw_nuisance_view(frame_f3, g = f1$cohorts[1L], K = f1$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl10)
  lk_f3 <- ns$construct_ipw_link_view(frame_f3, s = 10.0, t = 20.0, nuisance_view = nu_f3, control = own_ctrl10)
  prob_v3 <- ns$construct_certified_ipw_nuisance_problem_view(
    nu_f3, f1$beta_star_float, f1$F_float, 0.05, cert_ctrl_ref, prob_ctrl_ref, r_ctrl10
  )
  enc_v3 <- ns$enclose_ipw_link_view(prob_v3, lk_f3, c(0.0, 0.0), link_ctrl_ref, r_ctrl10)
  expect_gte(enc_v3$enclosures_star$T[2L], f3$expected_target_hi_min)
  expect_gt(enc_v3$enclosures_star$T[2L], 1.0)
  expect_lte(enc_v3$enclosures_star$omega[1L], f3$exact_omega_control_float)
  expect_gte(enc_v3$enclosures_star$omega[2L], f3$exact_omega_control_float)
})

# ---------------------------------------------------------------------------
# Test 3: Genuine physical vs logical shared-action & within-chunk scratch reuse
# ---------------------------------------------------------------------------

test_that("shared-xi chunk and within-chunk scratch genuinely reuse outcome-free primitives across compatible links", {
  ns <- require_cdid_namespace()
  ipw_ref <- dget(testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R"))
  f1 <- ipw_ref$fixture1_exact_rational

  adapted <- adapt_dense_fixture_to_panel_records(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G, cohorts = f1$cohorts,
    links = list(
      list(s = 10.0, t = 20.0, pair_observed = f1$pair_observed, y_left = f1$yL_1_float, y_right = f1$yR_1_float, substitute_off_risk = TRUE),
      list(s = 30.0, t = 40.0, pair_observed = f1$pair_observed, y_left = f1$yL_2_float, y_right = f1$yR_2_float, substitute_off_risk = TRUE)
    )
  )
  own_ctrl <- ns$ipw_owner_control(10L, 50000L, 50000L, 500000L, 50000L)
  panel <- ns$construct_ipw_panel_owner(adapted$panel_df, "id", "time", "y", "g", control = own_ctrl)
  frame <- ns$construct_ipw_frame_owner(panel, id = f1$id, Z = f1$Z_float, weights = f1$w_float, control = own_ctrl)
  nu_v <- ns$construct_ipw_nuisance_view(frame, g = f1$cohorts[1L], K = f1$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl)
  lk1_v <- ns$construct_ipw_link_view(frame, s = 10.0, t = 20.0, nuisance_view = nu_v, control = own_ctrl)
  lk2_v <- ns$construct_ipw_link_view(frame, s = 30.0, t = 40.0, nuisance_view = nu_v, control = own_ctrl)

  cert_ctrl <- ns$logit_certificate_control(18L, 40.0, 16L, 10L, 50000L, 500000L)
  prob_ctrl <- ns$ipw_nuisance_problem_control(10L, 50000L, 50000L, 500000L)
  link_ctrl <- ns$ipw_link_enclosure_control(10L, 50000L, 50000L, 500000L, 50000L)
  r_ctrl <- ns$ipw_reader_control(10L, 50000L, 50000L, 500000L, 50000L)

  prob_v <- ns$construct_certified_ipw_nuisance_problem_view(
    nu_v, f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, prob_ctrl, r_ctrl
  )
  enc1_v <- ns$enclose_ipw_link_view(prob_v, lk1_v, f1$link1$lambda_beta_float, link_ctrl, r_ctrl)
  enc2_v <- ns$enclose_ipw_link_view(prob_v, lk2_v, f1$link2$lambda_beta_float, link_ctrl, r_ctrl)

  R_H <- chol(f1$H_float)
  xi_v <- ns$enclose_ipw_shared_xi_view_chunk(prob_v, 1L, 10L, link_ctrl, r_ctrl, H_chol_candidate = R_H)
  expect_true(xi_v$certified)
  expect_identical(xi_v$read_diagnostics$endpoints_read, 0L)
  expect_gt(xi_v$work$physical_triangular_solves, 0L)
  expect_gt(xi_v$work$physical_exp_endpoint_evaluations, 0L)

  # Apply shared xi_v to link 1 and link 2: outcome-free row primitives (dual norm solves and exp evaluations)
  # and shared_xi_chunk are reused physically without recomputation!
  app1 <- ns$apply_ipw_shared_xi_to_link_view(enc1_v, xi_v, link_ctrl, r_ctrl, phi_candidate = f1$link1$phi_float)
  app2 <- ns$apply_ipw_shared_xi_to_link_view(enc2_v, xi_v, link_ctrl, r_ctrl, phi_candidate = f1$link2$phi_float)

  for (app in list(app1, app2)) {
    expect_identical(app$work$physical_triangular_solves, 0L)
    expect_identical(app$work$physical_exp_endpoint_evaluations, 0L)
    expect_gt(app$work$reused_triangular_solves, 0L)
    expect_gt(app$work$reused_exp_endpoint_evaluations, 0L)
    expect_gt(app$work$reused_shared_interval_ops, 0L)
    expect_identical(
      app$work$logical_interval_ops,
      app$work$physical_interval_ops + app$work$reused_shared_interval_ops
    )
    expect_identical(app$work$actual_interval_ops, app$work$logical_interval_ops)
    expect_lt(app$work$physical_interval_ops, app$work$logical_interval_ops)
  }

  # Also verify direct within-chunk scratch reuse across 2 links without shared_xi
  scratch <- ns$new_ipw_if_chunk_scratch(prob_v, 1L, 10L, control = link_ctrl)
  ch1_sc <- ns$enclose_ipw_link_if_view_chunk_with_scratch(
    enc1_v, 1L, 10L, link_ctrl, r_ctrl, scratch = scratch, phi_candidate = f1$link1$phi_float
  )
  ch2_sc <- ns$enclose_ipw_link_if_view_chunk_with_scratch(
    enc2_v, 1L, 10L, link_ctrl, r_ctrl, scratch = scratch, phi_candidate = f1$link2$phi_float
  )
  # First link populates scratch; second link reuses all triangular solves and exp evaluations
  expect_gt(ch1_sc$work$physical_triangular_solves, 0L)
  expect_identical(ch1_sc$work$reused_triangular_solves, 0L)
  expect_identical(ch2_sc$work$physical_triangular_solves, 0L)
  expect_identical(ch2_sc$work$reused_triangular_solves, ch1_sc$work$physical_triangular_solves)
  expect_identical(ch2_sc$work$physical_exp_endpoint_evaluations, 0L)
  expect_identical(ch2_sc$work$reused_exp_endpoint_evaluations, ch1_sc$work$physical_exp_endpoint_evaluations)
})

# ---------------------------------------------------------------------------
# Test 4: Compatibility, preflight budgets, character-ID byte charges,
# refusal precedence, and on-risk endpoint overflow with active scratch
# ---------------------------------------------------------------------------

test_that("compact consumers enforce frame/cohort/mask compatibility, combined preflight budgets, char-ID byte charges, and refusal precedence", {
  ns <- require_cdid_namespace()
  ipw_ref <- dget(testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R"))
  f1 <- ipw_ref$fixture1_exact_rational

  # Build panel with character IDs to test optional character-ID byte charges
  char_ids <- sprintf("unit_%02d", f1$id)
  adapted_char <- adapt_dense_fixture_to_panel_records(
    id = char_ids, Z = f1$Z_float, w = f1$w_float, G = f1$G, cohorts = f1$cohorts,
    links = list(
      list(s = 10.0, t = 20.0, pair_observed = f1$pair_observed, y_left = f1$yL_1_float, y_right = f1$yR_1_float, substitute_off_risk = TRUE),
      list(s = 30.0, t = 40.0, pair_observed = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), y_left = f1$yL_2_float, y_right = f1$yR_2_float, substitute_off_risk = TRUE)
    )
  )
  own_ctrl <- ns$ipw_owner_control(4L, 50000L, 50000L, 500000L, 50000L)
  panel_c <- ns$construct_ipw_panel_owner(adapted_char$panel_df, "id", "time", "y", "g", control = own_ctrl)
  frame_c <- ns$construct_ipw_frame_owner(panel_c, id = char_ids, Z = f1$Z_float, weights = f1$w_float, control = own_ctrl)
  nu_c <- ns$construct_ipw_nuisance_view(frame_c, g = f1$cohorts[1L], K = f1$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl)
  lk_c <- ns$construct_ipw_link_view(frame_c, s = 10.0, t = 20.0, nuisance_view = nu_c, control = own_ctrl)

  cert_ctrl <- ns$logit_certificate_control(18L, 40.0, 16L, 4L, 50000L, 500000L)
  prob_ctrl <- ns$ipw_nuisance_problem_control(4L, 50000L, 50000L, 500000L)
  link_ctrl <- ns$ipw_link_enclosure_control(4L, 50000L, 50000L, 500000L, 50000L)
  r_ctrl <- ns$ipw_reader_control(4L, 50000L, 50000L, 500000L, 50000L)

  # 4A. Optional character-ID expected_a_id byte charges in preflight
  pre_no_id <- ns$preflight_logit_certificate_view(nu_c, cert_ctrl, r_ctrl)
  pre_with_id <- ns$preflight_logit_certificate_view(nu_c, cert_ctrl, r_ctrl, expected_a_id = char_ids)
  # Each character ID "unit_XX" has 7 bytes -> 2 * 70 = 140 bytes per pass * 3 passes = 420 bytes -> 128 * 420 = 53760 extra work ops
  expect_identical(pre_with_id$id_comparison_bytes_per_pass, 140)
  expect_identical(
    pre_with_id$combined_work_ops_bound - pre_no_id$combined_work_ops_bound,
    3 * 128 * 140
  )

  # Tight reader work budget that allows outcome-free pass without ID check but refuses when char-ID byte charge is included
  tight_r_ctrl <- ns$ipw_reader_control(4L, 50000L, 50000L, pre_no_id$reader_work_ops_bound + 100L, 50000L)
  expect_s3_class(
    ns$certify_chained_logit_root_view(nu_c, f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, tight_r_ctrl),
    "cdid_logit_certificate"
  )
  tight_r_ctrl_id <- tight_r_ctrl
  attr(tight_r_ctrl_id, "expected_a_id") <- char_ids
  expect_error(
    ns$certify_chained_logit_root_view(nu_c, f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, tight_r_ctrl_id),
    class = "cdid_ipw_owner_bound_error"
  )

  # 4B. Preflight budget refusal happens BEFORE any endpoint reads
  tr_pref <- new.env(parent = emptyenv())
  prob_c <- ns$construct_certified_ipw_nuisance_problem_view(
    nu_c, f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, prob_ctrl, r_ctrl
  )
  tiny_link_ctrl <- ns$ipw_link_enclosure_control(4L, 50000L, 50000L, 100L, 50000L)
  r_tr <- ns$ipw_reader_control(4L, 50000L, 50000L, 500000L, 50000L)
  attr(r_tr, "trace_env") <- tr_pref
  expect_error(
    ns$enclose_ipw_link_view(prob_c, lk_c, f1$link1$lambda_beta_float, tiny_link_ctrl, r_tr),
    class = "cdid_ipw_enclosure_bound_error"
  )
  expect_null(tr_pref$calls)

  # 4C. Compatibility refusals: different frame, permuted cohorts, unverified mask, mismatched mask
  frame_c2 <- ns$construct_ipw_frame_owner(panel_c, id = char_ids, Z = f1$Z_float, weights = f1$w_float, control = own_ctrl, frame_id = "other_fr")
  lk_other_frame <- ns$construct_ipw_link_view(frame_c2, s = 10.0, t = 20.0, g = f1$cohorts[1L], K = f1$cohorts[-1L], control = own_ctrl)
  expect_error(
    ns$enclose_ipw_link_view(prob_c, lk_other_frame, c(0, 0), link_ctrl, r_ctrl),
    class = "cdid_ipw_owner_mismatch_error"
  )

  expect_error(
    ns$construct_ipw_link_view(frame_c, s = 10.0, t = 20.0, g = f1$cohorts[1L], K = rev(f1$cohorts[-1L]), nuisance_view = nu_c, control = own_ctrl),
    class = "cdid_ipw_owner_cohort_error"
  )

  lk_unverified <- ns$construct_ipw_link_view(
    frame_c, s = 10.0, t = 20.0, nuisance_view = nu_c, control = own_ctrl, verify_nuisance_mask = FALSE
  )
  expect_error(
    ns$enclose_ipw_link_view(prob_c, lk_unverified, c(0, 0), link_ctrl, r_ctrl),
    class = "cdid_ipw_shared_incompatible_error"
  )

  # Incompatible mask on (30, 40) where u=2 is not observed
  expect_error(
    ns$construct_ipw_link_view(frame_c, s = 30.0, t = 40.0, nuisance_view = nu_c, control = own_ctrl, verify_nuisance_mask = TRUE),
    class = "cdid_ipw_owner_mask_mismatch_error"
  )

  # Candidate mismatch against stored view candidates
  nu_with_cand <- ns$construct_ipw_nuisance_view(
    frame_c, g = f1$cohorts[1L], K = f1$cohorts[-1L], s = 10.0, t = 20.0,
    control = own_ctrl, beta0 = c(0.0, 0.0), F = f1$F_float, r = 0.05
  )
  expect_error(
    ns$certify_chained_logit_root_view(nu_with_cand, c(0.1, 0.0), f1$F_float, 0.05, cert_ctrl, r_ctrl),
    class = "cdid_logit_certificate_candidate_mismatch_error"
  )
  expect_error(
    ns$construct_certified_ipw_nuisance_problem_view(nu_with_cand, c(0.1, 0.0), f1$F_float, 0.05, cert_ctrl, prob_ctrl, r_ctrl),
    class = "cdid_ipw_enclosure_candidate_mismatch_error"
  )

  # Unsupported cohort on view constructed with require_support = FALSE
  nu_unsup <- ns$construct_ipw_nuisance_view(
    frame_c, g = f1$cohorts[1L], K = c(f1$cohorts[-1L], 999.0), s = 10.0, t = 20.0,
    control = own_ctrl, require_support = FALSE
  )
  expect_error(
    ns$certify_chained_logit_root_view(nu_unsup, c(0.0, 0.0), f1$F_float, 0.05, cert_ctrl, r_ctrl),
    class = "cdid_logit_certificate_cohort_support_error"
  )

  # 4D. On-risk endpoint difference overflow produces INCONCLUSIVE even when reusing chunk_scratch,
  # while structural expected_a_id mismatch on a later chunk takes precedence over chunk 1 arithmetic refusal!
  yL_ovf <- f1$yL_1_float
  yR_ovf <- f1$yR_1_float
  yL_ovf[1L] <- -1e308
  yR_ovf[1L] <- +1e308
  adapted_ovf <- adapt_dense_fixture_to_panel_records(
    id = char_ids, Z = f1$Z_float, w = f1$w_float, G = f1$G, cohorts = f1$cohorts,
    links = list(
      list(s = 10.0, t = 20.0, pair_observed = f1$pair_observed, y_left = f1$yL_1_float, y_right = f1$yR_1_float, substitute_off_risk = TRUE),
      list(s = 30.0, t = 40.0, pair_observed = f1$pair_observed, y_left = yL_ovf, y_right = yR_ovf, substitute_off_risk = TRUE)
    )
  )
  panel_ovf <- ns$construct_ipw_panel_owner(adapted_ovf$panel_df, "id", "time", "y", "g", control = own_ctrl)
  frame_ovf <- ns$construct_ipw_frame_owner(panel_ovf, id = char_ids, Z = f1$Z_float, weights = f1$w_float, control = own_ctrl)
  nu_ovf <- ns$construct_ipw_nuisance_view(frame_ovf, g = f1$cohorts[1L], K = f1$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl)
  lk1_ovf <- ns$construct_ipw_link_view(frame_ovf, s = 10.0, t = 20.0, nuisance_view = nu_ovf, control = own_ctrl)
  lk2_ovf <- ns$construct_ipw_link_view(frame_ovf, s = 30.0, t = 40.0, nuisance_view = nu_ovf, control = own_ctrl)
  prob_ovf <- ns$construct_certified_ipw_nuisance_problem_view(
    nu_ovf, f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, prob_ctrl, r_ctrl
  )

  # Link 2 has on-risk endpoint overflow on unit 1 (chunk 1): returns INCONCLUSIVE when IDs are valid
  enc2_ovf <- ns$enclose_ipw_link_view(prob_ovf, lk2_ovf, c(0.0, 0.0), link_ctrl, r_ctrl)
  expect_false(enc2_ovf$certified)
  expect_identical(enc2_ovf$status, "INCONCLUSIVE")
  expect_identical(enc2_ovf$inconclusive_reasons, "link_arithmetic_refusal")

  # But if expected_a_id has a structural mismatch on chunk 2 (unit 6), structural error takes precedence over chunk 1 overflow!
  bad_ids <- char_ids
  bad_ids[6L] <- "wrong_id_on_chunk_2"
  r_ctrl_bad_id <- r_ctrl
  attr(r_ctrl_bad_id, "expected_a_id") <- bad_ids
  expect_error(
    ns$enclose_ipw_link_view(prob_ovf, lk2_ovf, c(0.0, 0.0), link_ctrl, r_ctrl_bad_id),
    class = "cdid_ipw_owner_order_mismatch_error"
  )
})

# ---------------------------------------------------------------------------
# Test 5: Inconclusive propagation, IF chunk overflow with active scratch,
# compact seal integrity, interrupt checks, and zero global/options leaks
# ---------------------------------------------------------------------------

test_that("inconclusive propagation, IF chunk overflow with active scratch, compact seal integrity, and zero global leaks hold", {
  ns <- require_cdid_namespace()
  ipw_ref <- dget(testthat::test_path("fixtures", "ipw_enclosure_reference_fixtures.R"))
  f1 <- ipw_ref$fixture1_exact_rational

  global_before <- ls(envir = globalenv(), all.names = TRUE)
  opts_before <- options()

  adapted <- adapt_dense_fixture_to_panel_records(
    id = f1$id, Z = f1$Z_float, w = f1$w_float, G = f1$G, cohorts = f1$cohorts,
    links = list(list(
      s = 10.0, t = 20.0, pair_observed = f1$pair_observed,
      y_left = f1$yL_1_float, y_right = f1$yR_1_float, substitute_off_risk = TRUE
    ))
  )
  own_ctrl <- ns$ipw_owner_control(4L, 50000L, 50000L, 500000L, 50000L)
  panel <- ns$construct_ipw_panel_owner(adapted$panel_df, "id", "time", "y", "g", control = own_ctrl)
  frame <- ns$construct_ipw_frame_owner(panel, id = f1$id, Z = f1$Z_float, weights = f1$w_float, control = own_ctrl)
  nu_v <- ns$construct_ipw_nuisance_view(frame, g = f1$cohorts[1L], K = f1$cohorts[-1L], s = 10.0, t = 20.0, control = own_ctrl)
  lk_v <- ns$construct_ipw_link_view(frame, s = 10.0, t = 20.0, nuisance_view = nu_v, control = own_ctrl)

  cert_ctrl <- ns$logit_certificate_control(18L, 40.0, 16L, 4L, 50000L, 500000L)
  prob_ctrl <- ns$ipw_nuisance_problem_control(4L, 50000L, 50000L, 500000L)
  link_ctrl <- ns$ipw_link_enclosure_control(10L, 50000L, 50000L, 500000L, 50000L)
  r_ctrl <- ns$ipw_reader_control(10L, 50000L, 50000L, 500000L, 50000L)

  # 5A. Domain/contraction failure returns structured INCONCLUSIVE and propagates with 0 endpoint reads
  prob_inc <- ns$construct_certified_ipw_nuisance_problem_view(
    nu_v, c(3.0, -3.0), f1$F_float, 0.01, cert_ctrl, prob_ctrl, r_ctrl
  )
  expect_false(prob_inc$certified)
  expect_identical(prob_inc$status, "INCONCLUSIVE")
  expect_false(prob_inc$separation_certified)
  expect_false(prob_inc$sampling_certified)
  expect_false(prob_inc$public_gate_adopted)
  ns$ipw_validate_nuisance_problem_mode(prob_inc, require_certified = FALSE, expected_storage_mode = "compact_view")
  expect_error(
    ns$ipw_validate_nuisance_problem_mode(prob_inc, require_certified = TRUE),
    class = "cdid_ipw_enclosure_inconclusive"
  )

  tr_inc <- new.env(parent = emptyenv())
  r_tr_inc <- ns$ipw_reader_control(10L, 50000L, 50000L, 500000L, 50000L)
  attr(r_tr_inc, "trace_env") <- tr_inc
  enc_inc <- ns$enclose_ipw_link_view(prob_inc, lk_v, c(0.0, 0.0), link_ctrl, r_tr_inc)
  expect_false(enc_inc$certified)
  expect_identical(enc_inc$status, "INCONCLUSIVE")
  expect_identical(enc_inc$read_diagnostics$endpoints_read, 0L)
  expect_identical(tr_inc$calls, 0L)

  xi_inc <- ns$enclose_ipw_shared_xi_view_chunk(prob_inc, 1L, 10L, link_ctrl, r_tr_inc, H_chol_candidate = chol(f1$H_float))
  expect_false(xi_inc$certified)
  expect_identical(xi_inc$status, "INCONCLUSIVE")

  ch_inc <- ns$enclose_ipw_link_if_view_chunk(enc_inc, 1L, 10L, link_ctrl, r_tr_inc)
  expect_false(ch_inc$certified)
  expect_identical(ch_inc$status, "INCONCLUSIVE")
  expect_identical(tr_inc$calls, 0L)

  # 5B. Tampering with sealed compact problem, link enclosure, shared_xi_chunk, or scratch is refused
  prob_ok <- ns$construct_certified_ipw_nuisance_problem_view(
    nu_v, f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, prob_ctrl, r_ctrl
  )
  enc_ok <- ns$enclose_ipw_link_view(prob_ok, lk_v, f1$link1$lambda_beta_float, link_ctrl, r_ctrl)
  xi_ok <- ns$enclose_ipw_shared_xi_view_chunk(prob_ok, 1L, 10L, link_ctrl, r_ctrl, H_chol_candidate = chol(f1$H_float))

  prob_tamp <- prob_ok
  prob_tamp$nuisance_enclosures$kappa_star_lower <- 0.999
  expect_error(
    ns$ipw_validate_nuisance_problem_mode(prob_tamp),
    class = "cdid_ipw_enclosure_integrity_error"
  )

  enc_tamp <- enc_ok
  enc_tamp$enclosures_star$mu <- c(-999.0, 999.0)
  expect_error(
    ns$ipw_validate_link_enclosure_mode(enc_tamp),
    class = "cdid_ipw_enclosure_integrity_error"
  )

  xi_tamp <- xi_ok
  xi_tamp$xi_rate[1L, 1L, 1L] <- -999.0
  expect_error(
    ns$apply_ipw_shared_xi_to_link_view(enc_ok, xi_tamp, link_ctrl, r_ctrl),
    class = "cdid_ipw_enclosure_integrity_error"
  )

  # Passing a dense problem to compact view entrypoint (or vice versa) is refused
  prob_dense <- ns$construct_certified_ipw_nuisance_problem(
    f1$id, f1$Z_float, f1$w_float, f1$G, f1$pair_observed, f1$cohorts,
    f1$beta_star_float, f1$F_float, 0.05, cert_ctrl, prob_ctrl
  )
  expect_error(
    ns$enclose_ipw_link_view(prob_dense, lk_v, c(0.0, 0.0), link_ctrl, r_ctrl),
    class = "cdid_ipw_enclosure_integrity_error"
  )
  expect_error(
    ns$enclose_ipw_link(prob_ok, f1$id, f1$yL_1_float, f1$yR_1_float, c(0.0, 0.0), link_ctrl),
    class = "cdid_ipw_enclosure_integrity_error"
  )

  global_after <- ls(envir = globalenv(), all.names = TRUE)
  expect_identical(global_after, global_before)
  expect_identical(options(), opts_before)
})
