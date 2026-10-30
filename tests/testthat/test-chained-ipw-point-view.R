# Positive cases use an explicit 1e12 operation/5e6 entry diagnostic profile
# for corrected whole-traversal and component-composition accounting. Numerical
# inputs, tolerances, reference targets, and low-cap negative predicates are unchanged.
# Focused tests for the private CP02 compact point-link, bounded complete
# point-IF replay, covariance candidate_rows callback producer, and
# candidate-to-certificate root/link/IF binding (`R/chained_ipw_point_view.R`).

if(!requireNamespace("cdid",quietly=TRUE)) stop("Actual installed cdid required; no source fallback.",call.=FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api),environmentIsLocked(api))
cp02_pv_fn <- function(name) {
  f <- get(name, envir = api, inherits = FALSE)
  if (!is.function(f) || !identical(environment(f), api)) {
    stop("Foreign function binding: ", name, call. = FALSE)
  }
  f
}

cp02_pv_locate_file <- function(...) {
  parts <- c(...)
  if (identical(parts, c("reference", "panel24.csv"))) {
    return(testthat::test_path("fixtures", "fitted_panel24.csv"))
  }
  if (identical(parts, c("tests", "testthat", "reference", "cp02_shared_ipw_audit_newton_reference.R"))) {
    return(testthat::test_path("reference", "cp02_shared_ipw_audit_newton_reference.R"))
  }
  stop("Unlisted test-local reference path.", call. = FALSE)
}

cp02_pv_panel24_fixture <- function(weighted = TRUE) {
  p24_path <- cp02_pv_locate_file("reference", "panel24.csv")
  p24 <- utils::read.csv(p24_path, stringsAsFactors = FALSE)
  u24 <- p24[!duplicated(p24$id), ]
  Z24 <- cbind(one = 1, sign = as.double(u24$sign))
  w24 <- if (weighted) 0.8 + ((seq_len(24L) * 7L) %% 11L) / 20 else rep(1, 24L)
  idx24 <- cp02_pv_fn("construct_observed_pairs")(p24, "id", "time", "y", "g")
  ctx24 <- cp02_pv_fn("construct_ipw_context")(idx24, id = u24$id, Z = Z24, weights = w24)
  list(panel = p24, units = u24, Z = Z24, w = w24, index = idx24, context = ctx24)
}

cp02_pv_weighted96_fixture <- function(weighted = TRUE) {
  n <- 96L
  G <- rep(c(4, 7, Inf, 9), each = 24L)
  x <- rep(seq(-1.4, 1.4, length.out = 24L), 4L)
  Z <- cbind(one = 1, x = x)
  w <- if (weighted) 0.7 + ((seq_len(n) * 13L) %% 29L) / 20 else rep(1, n)
  obs <- matrix(FALSE, n, 8L)
  for (i in seq_len(n)) {
    k <- ((i - 1L) %% 24L) + 1L
    if (G[i] != 9 && k <= 18L) {
      obs[i, 1:6] <- TRUE
    } else {
      obs[i, 7:8] <- TRUE
      obs[i, 1L + (k %% 6L)] <- TRUE
    }
  }
  Y <- outer(seq_len(n) / 17, 1:8, function(a, b) sin(a * b) + a + b^2 / 13) +
    outer(x, 1:8, function(a, b) a * b / 3)
  pan <- data.frame(
    id = rep(seq_len(n), each = 8L),
    time = rep(as.double(1:8), n),
    y = as.vector(t(ifelse(obs, Y, NA_real_))),
    g = rep(G, each = 8L)
  )
  idx <- cp02_pv_fn("construct_observed_pairs")(pan, "id", "time", "y", "g")
  ctx <- cp02_pv_fn("construct_ipw_context")(idx, id = seq_len(n), Z = Z, weights = w)
  list(n = n, G = G, Z = Z, w = w, obs = obs, Y = Y, panel = pan, index = idx, context = ctx)
}

cp02_pv_newton_ref <- function() {
  ref_path <- cp02_pv_locate_file("tests", "testthat", "reference", "cp02_shared_ipw_audit_newton_reference.R")
  if (!identical(
    digest::digest(file = ref_path, algo = "sha256"),
    "2e5647c6b0530e83f7701aae00684fddbf5f72ac37e6000c71c6a5bfffd94b72"
  )) {
    stop("packaged reference bytes differ: cp02_shared_ipw_audit_newton_reference.R", call. = FALSE)
  }
  env <- new.env(parent = baseenv())
  sys.source(ref_path, envir = env)
  env$independent_fit
}

# Local instrumentation clones exact installed functions without altering namespace bindings.
pv_local_api <- function() {
  e <- new.env(parent = api)
  for (nm in c(
    "ipw_owner_control", "construct_ipw_panel_owner", "construct_ipw_frame_owner",
    "construct_ipw_nuisance_view", "construct_ipw_link_view", "query_ipw_view_chunk",
    "ipw_fitted_view_control", "fit_ipw_nuisance_view", "fit_ipw_nuisance_from_owner",
    "validate_ipw_fitted_nuisance_view", "verify_ipw_fitted_nuisance_sharing",
    "ipw_glm_fit_internal", "ipw_link_means", "ipw_link_partials", "ipw_nuisance_observer_fit_worker"
  )) {
    assign(nm, get(nm, envir = api, inherits = FALSE), envir = e)
  }
  for (nm in c(
    "cdid_ipw_point_view_abort", "ipw_point_view_require_deps",
    "ipw_point_view_control", "validate_ipw_point_view_control",
    "ipw_point_view_to_owner_control", "ipw_point_view_to_fitted_control",
    "ipw_point_view_char_bytes", "ipw_point_view_id_bytes_bound",
    "ipw_point_view_identical", "ipw_point_view_payload_budget",
    "ipw_point_view_seal_gate", "ipw_point_view_owner_entries",
    "preflight_ipw_point_certificate_binding",
    "ipw_point_view_verify_compact_storage", "preflight_ipw_point_link_view",
    "fit_ipw_point_link_view", "construct_ipw_point_link_view",
    "fit_ipw_point_link_from_owner", "fit_ipw_point_link_view_set",
    "validate_ipw_point_link_view", "extract_ipw_point_link_snapshot",
    "ipw_point_link_snapshot_payload", "revalidate_ipw_point_link_view",
    "normalize_ipw_point_links_view", "preflight_ipw_point_if_view",
    "ipw_point_if_eval_chunk_internal", "evaluate_ipw_point_if_view_chunk",
    "query_ipw_point_if_view_chunk", "replay_ipw_point_if_view",
    "traverse_ipw_point_if_view", "make_ipw_point_candidate_rows_spec",
    "stream_ipw_point_covariance_view", "bind_fitted_nuisance_to_problem_view",
    "bind_point_link_to_enclosure_view", "enclose_point_link_if_view_chunk",
    "certify_ipw_point_link_view"
  )) {
    f <- get(nm, envir = api, inherits = FALSE)
    environment(f) <- e
    assign(nm, f, envir = e)
  }
  e
}

test_that("compact point link matches dense fit_ipw_link, fit_ipw_link_shared, and Newton reference bitwise across chunk sizes without solver calls", {
  construct_ipw_panel_owner <- cp02_pv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_pv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_pv_fn("construct_ipw_nuisance_view")
  construct_ipw_link_view <- cp02_pv_fn("construct_ipw_link_view")
  fit_ipw_nuisance_view <- cp02_pv_fn("fit_ipw_nuisance_view")
  fit_ipw_nuisance_from_owner <- cp02_pv_fn("fit_ipw_nuisance_from_owner")
  fit_ipw_point_link_view <- cp02_pv_fn("fit_ipw_point_link_view")
  fit_ipw_point_link_from_owner <- cp02_pv_fn("fit_ipw_point_link_from_owner")
  fit_ipw_point_link_view_set <- cp02_pv_fn("fit_ipw_point_link_view_set")
  validate_ipw_point_link_view <- cp02_pv_fn("validate_ipw_point_link_view")
  ipw_owner_control <- cp02_pv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_pv_fn("ipw_fitted_view_control")
  ipw_point_view_control <- cp02_pv_fn("ipw_point_view_control")
  fit_ipw_link <- cp02_pv_fn("fit_ipw_link")
  construct_ipw_shared_session <- cp02_pv_fn("construct_ipw_shared_session")
  construct_ipw_nuisance <- cp02_pv_fn("construct_ipw_nuisance")
  declare_ipw_shared_link <- cp02_pv_fn("declare_ipw_shared_link")
  fit_ipw_link_shared <- cp02_pv_fn("fit_ipw_link_shared")

  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  # 1. Deterministic 24-unit panel (unweighted and weighted across all 6 supported pairs)
  pairs24 <- list(c(1, 2), c(1, 4), c(2, 3), c(2, 4), c(3, 4), c(4, 5))
  for (weighted in c(FALSE, TRUE)) {
    fx24 <- cp02_pv_panel24_fixture(weighted = weighted)
    own_ctrl24 <- ipw_owner_control(7, 20000, 100000, 1000000, 20000)
    fv_ctrl24 <- ipw_fitted_view_control(7, 24, 20000, 100000, 1000000, 20000)
    po24 <- construct_ipw_panel_owner(fx24$panel, "id", "time", "y", "g", control = own_ctrl24)
    fo24 <- construct_ipw_frame_owner(po24, id = fx24$units$id, Z = fx24$Z, weights = fx24$w, control = own_ctrl24)

    for (pr in pairs24) {
      dlink <- fit_ipw_link(fx24$context, g = 3, K = Inf, s = pr[1], t = pr[2], control = glm_ctrl)
      vfit24 <- fit_ipw_nuisance_from_owner(
        fo24, g = 3, K = Inf, s = pr[1], t = pr[2],
        glm_control = glm_ctrl, control = fv_ctrl24
      )
      for (ch_sz in c(1L, 5L, 24L)) {
        pv_ctrl <- ipw_point_view_control(
          max_chunk_rows = ch_sz,
          max_observer_units = 24,
          max_retained_entries = 20000,
          max_matrix_entries = 100000,
          max_work_ops = 1000000000000,
          max_output_entries = 20000
        )
        plink <- fit_ipw_point_link_from_owner(
          vfit24, s = pr[1], t = pr[2],
          control = pv_ctrl, g = 3, K = Inf, glm_control = glm_ctrl
        )
        validate_ipw_point_link_view(plink, fitted_view = vfit24, frame_owner = fo24)
        expect_identical(plink$estimate, dlink$estimate)
        expect_identical(plink$treated_mean, dlink$treated_mean)
        expect_identical(plink$control_mean, dlink$control_mean)
        expect_identical(plink$means, dlink$means)
        expect_identical(plink$partials, dlink$partials)
        expect_identical(plink$loadings, dlink$loadings)
        expect_identical(plink$counts, dlink$counts)
      }
    }
  }

  # 2. Weighted 96-unit fixture: shared nuisance across (1,2), (2,3), (1,3), distinct (7,8),
  #    distinct control pools, and Newton reference comparison
  ind_fit <- cp02_pv_newton_ref()
  for (weighted in c(FALSE, TRUE)) {
    fx96 <- cp02_pv_weighted96_fixture(weighted = weighted)
    own_ctrl96 <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
    fv_ctrl96 <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
    po96 <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl96)
    fo96 <- construct_ipw_frame_owner(po96, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl96)

    sess <- construct_ipw_shared_session(fx96$context, max_nuisance = 8L, max_links = 16L, max_unit_entries = 10000L)
    nu12 <- construct_ipw_nuisance(sess, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl)
    nu78 <- construct_ipw_nuisance(sess, g = 4, K = c(7, Inf), s = 7, t = 8, control = glm_ctrl)
    decl12 <- declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl)
    decl23 <- declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 2, t = 3, control = glm_ctrl)
    decl13 <- declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 1, t = 3, control = glm_ctrl)
    decl78 <- declare_ipw_shared_link(sess, nu78, g = 4, K = c(7, Inf), s = 7, t = 8, control = glm_ctrl)
    d12 <- fit_ipw_link_shared(sess, decl12)
    d23 <- fit_ipw_link_shared(sess, decl23)
    d13 <- fit_ipw_link_shared(sess, decl13)
    d78 <- fit_ipw_link_shared(sess, decl78)
    d12_all_ctrl <- fit_ipw_link(fx96$context, g = 4, K = Inf, s = 1, t = 2, control = glm_ctrl)

    nv12 <- construct_ipw_nuisance_view(fo96, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl96)
    nv78 <- construct_ipw_nuisance_view(fo96, g = 4, K = c(7, Inf), s = 7, t = 8, control = own_ctrl96)
    nv12_all <- construct_ipw_nuisance_view(fo96, g = 4, K = Inf, s = 1, t = 2, control = own_ctrl96)

    vfit12 <- fit_ipw_nuisance_view(nv12, glm_control = glm_ctrl, control = fv_ctrl96)
    vfit78 <- fit_ipw_nuisance_view(nv78, glm_control = glm_ctrl, control = fv_ctrl96)
    vfit12_all <- fit_ipw_nuisance_view(nv12_all, glm_control = glm_ctrl, control = fv_ctrl96)

    lv12 <- construct_ipw_link_view(fo96, g = 4, K = c(7, Inf), s = 1, t = 2, nuisance_view = nv12, control = own_ctrl96)
    lv23 <- construct_ipw_link_view(fo96, g = 4, K = c(7, Inf), s = 2, t = 3, nuisance_view = nv12, control = own_ctrl96)
    lv13 <- construct_ipw_link_view(fo96, g = 4, K = c(7, Inf), s = 1, t = 3, nuisance_view = nv12, control = own_ctrl96)
    lv78 <- construct_ipw_link_view(fo96, g = 4, K = c(7, Inf), s = 7, t = 8, nuisance_view = nv78, control = own_ctrl96)
    lv12_all <- construct_ipw_link_view(fo96, g = 4, K = Inf, s = 1, t = 2, nuisance_view = nv12_all, control = own_ctrl96)

    for (ch_sz in c(1L, 5L, 8L, 13L, 96L)) {
      pv_ctrl <- ipw_point_view_control(ch_sz, 96, 100000, 5000000, 1000000000000, 100000)
      tr <- new.env(parent = emptyenv())
      pset <- fit_ipw_point_link_view_set(
        fitted_view = list(vfit12, vfit12, vfit12, vfit78, vfit12_all),
        link_views = list(lv12, lv23, lv13, lv78, lv12_all),
        control = pv_ctrl,
        trace_env = tr
      )
      expect_identical(tr$propensity_solver_calls, 0L)
      expect_identical(tr$n_links_fitted, 5L)
      expect_true(tr$all_released_between_links)

      ref_links <- list(d12, d23, d13, d78, d12_all_ctrl)
      for (idx_l in seq_along(pset)) {
        pl <- pset[[idx_l]]
        dl <- ref_links[[idx_l]]
        expect_identical(pl$estimate, dl$estimate)
        expect_identical(pl$treated_mean, dl$treated_mean)
        expect_identical(pl$control_mean, dl$control_mean)
        expect_identical(pl$means, dl$means)
        expect_identical(pl$partials, dl$partials)
        expect_identical(pl$loadings, dl$loadings)
        expect_identical(pl$counts, dl$counts)
      }

      # Additivity along shared observer mask: mu(1,2) + mu(2,3) == mu(1,3)
      expect_equal(pset[[1L]]$estimate + pset[[2L]]$estimate, pset[[3L]]$estimate, tolerance = 1e-15)
      # Distinct control pool produces genuinely distinct point link
      expect_gt(abs(pset[[1L]]$estimate - pset[[5L]]$estimate), 1e-3)

      # Compare against independent Newton reference
      nref12 <- ind_fit(fx96, 1, 2)
      expect_equal(pset[[1L]]$estimate, nref12$estimate, tolerance = 1e-11)
      expect_equal(unname(pset[[1L]]$loadings$loading_beta), unname(nref12$lb), tolerance = 1e-10)
      expect_equal(unname(pset[[1L]]$loadings$loading_rate), unname(nref12$lr), tolerance = 1e-10)
    }
  }

  # 3. Spy on ipw_glm_fit_internal in local instrumented API confirms zero solver calls
  local_api <- pv_local_api()
  glm_calls <- 0L
  orig_glm <- local_api$ipw_glm_fit_internal
  local_api$ipw_glm_fit_internal <- function(...) {
    glm_calls <<- glm_calls + 1L
    orig_glm(...)
  }
  fx96 <- cp02_pv_weighted96_fixture(weighted = TRUE)
  oc <- local_api$ipw_owner_control(13, 100000, 500000, 5000000, 100000)
  fc <- local_api$ipw_fitted_view_control(13, 96, 100000, 500000, 5000000, 100000)
  pc <- local_api$ipw_point_view_control(13, 96, 100000, 5000000, 1000000000000, 100000)
  po <- local_api$construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = oc)
  fo <- local_api$construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = oc)
  nv <- local_api$construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = oc)
  lv <- local_api$construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, nuisance_view = nv, control = oc)
  vf <- local_api$fit_ipw_nuisance_view(nv, glm_control = glm_ctrl, control = fc)
  glm_calls <- 0L
  pl <- local_api$fit_ipw_point_link_view(vf, lv, control = pc)
  expect_identical(glm_calls, 0L)
})

test_that("bounded complete 4-term point-IF replay and covariance stream match dense paths bitwise, including nonobserver rates, risk-first endpoint filtering, and character/unobserved A-units", {
  construct_ipw_panel_owner <- cp02_pv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_pv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_pv_fn("construct_ipw_nuisance_view")
  construct_ipw_link_view <- cp02_pv_fn("construct_ipw_link_view")
  fit_ipw_nuisance_view <- cp02_pv_fn("fit_ipw_nuisance_view")
  fit_ipw_nuisance_from_owner <- cp02_pv_fn("fit_ipw_nuisance_from_owner")
  fit_ipw_point_link_view_set <- cp02_pv_fn("fit_ipw_point_link_view_set")
  fit_ipw_point_link_from_owner <- cp02_pv_fn("fit_ipw_point_link_from_owner")
  evaluate_ipw_point_if_view_chunk <- cp02_pv_fn("evaluate_ipw_point_if_view_chunk")
  replay_ipw_point_if_view <- cp02_pv_fn("replay_ipw_point_if_view")
  stream_ipw_point_covariance_view <- cp02_pv_fn("stream_ipw_point_covariance_view")
  ipw_owner_control <- cp02_pv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_pv_fn("ipw_fitted_view_control")
  ipw_point_view_control <- cp02_pv_fn("ipw_point_view_control")
  construct_ipw_shared_session <- cp02_pv_fn("construct_ipw_shared_session")
  construct_ipw_nuisance <- cp02_pv_fn("construct_ipw_nuisance")
  declare_ipw_shared_link <- cp02_pv_fn("declare_ipw_shared_link")
  fit_ipw_link_shared <- cp02_pv_fn("fit_ipw_link_shared")
  materialize_ipw_if <- cp02_pv_fn("materialize_ipw_if")
  stream_ipw_covariance <- cp02_pv_fn("stream_ipw_covariance")

  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)
  fx96 <- cp02_pv_weighted96_fixture(weighted = TRUE)

  # Off-risk corruption test: corrupt endpoints for G == 9 (outside risk set g=4, K=c(7, Inf))
  pan_corrupt <- fx96$panel
  pan_corrupt$y[pan_corrupt$g == 9 & pan_corrupt$time %in% c(1, 2, 3)] <- 1e308

  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
  po <- construct_ipw_panel_owner(pan_corrupt, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)

  nv12 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  nv78 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 7, t = 8, control = own_ctrl)
  vfit12 <- fit_ipw_nuisance_view(nv12, glm_control = glm_ctrl, control = fv_ctrl)
  vfit78 <- fit_ipw_nuisance_view(nv78, glm_control = glm_ctrl, control = fv_ctrl)

  lv12 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, nuisance_view = nv12, control = own_ctrl)
  lv23 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 2, t = 3, nuisance_view = nv12, control = own_ctrl)
  lv13 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 1, t = 3, nuisance_view = nv12, control = own_ctrl)
  lv78 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 7, t = 8, nuisance_view = nv78, control = own_ctrl)

  pv_ctrl <- ipw_point_view_control(16, 96, 100000, 5000000, 1000000000000, 100000)
  plinks <- fit_ipw_point_link_view_set(
    list(vfit12, vfit12, vfit12, vfit78),
    list(lv12, lv23, lv13, lv78),
    control = pv_ctrl
  )

  sess <- construct_ipw_shared_session(fx96$context, max_nuisance = 8L, max_links = 16L, max_unit_entries = 10000L)
  nu12 <- construct_ipw_nuisance(sess, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl)
  nu78 <- construct_ipw_nuisance(sess, g = 4, K = c(7, Inf), s = 7, t = 8, control = glm_ctrl)
  d12 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl))
  d23 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 2, t = 3, control = glm_ctrl))
  d13 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu12, g = 4, K = c(7, Inf), s = 1, t = 3, control = glm_ctrl))
  d78 <- fit_ipw_link_shared(sess, declare_ipw_shared_link(sess, nu78, g = 4, K = c(7, Inf), s = 7, t = 8, control = glm_ctrl))
  dlinks <- list(d12, d23, d13, d78)
  dense_mat <- materialize_ipw_if(fx96$context, dlinks, max_entries = 10000, chunk_size = 16L, include_components = TRUE)

  # Verify risk nonobserver rate term is nonzero on (1,2) and outside-risk rows are exact zero
  full_12 <- evaluate_ipw_point_if_view_chunk(
    list(plinks[[1L]]), a_start = 1L, a_end = 96L,
    control = ipw_point_view_control(96, 96, 100000, 5000000, 1000000000000, 100000),
    include_components = TRUE
  )
  in_risk_12 <- fx96$G %in% c(4, 7, Inf)
  in_obs_12 <- in_risk_12 & fx96$obs[, 1] & fx96$obs[, 2]
  nonobs_risk <- which(in_risk_12 & !in_obs_12)
  outside_risk <- which(!in_risk_12)
  expect_gt(length(nonobs_risk), 0L)
  expect_true(all(abs(full_12$phi[nonobs_risk, 1L]) > 1e-6))
  expect_identical(full_12$phi[nonobs_risk, 1L], full_12$components$rate[nonobs_risk, 1L])
  expect_true(all(full_12$components$treated[nonobs_risk, 1L] == 0))
  expect_true(all(full_12$components$control[nonobs_risk, 1L] == 0))
  expect_true(all(full_12$components$beta_score[nonobs_risk, 1L] == 0))
  expect_true(all(full_12$phi[outside_risk, 1L] == 0))

  dense_replay_num_ok <- function(actual, expected, tol = 1e-12) {
    identical(typeof(actual), "double") &&
      identical(typeof(expected), "double") &&
      identical(is.object(actual), is.object(expected)) &&
      identical(attributes(actual), attributes(expected)) &&
      identical(length(actual), length(expected)) &&
      length(actual) > 0L &&
      all(is.finite(actual)) &&
      all(is.finite(expected)) &&
      all(actual[expected == 0] == 0) &&
      all(abs(actual - expected) <= tol * pmax(1, abs(actual), abs(expected)))
  }
  expect_dense_replay_num <- function(actual, expected, tol = 1e-12, info = NULL) {
    expect_true(dense_replay_num_ok(actual, expected, tol = tol), info = info)
  }
  expect_false(dense_replay_num_ok(numeric(), numeric()))
  expect_false(dense_replay_num_ok(matrix(numeric(), 0L, 4L), matrix(numeric(), 0L, 4L)))
  expect_false(dense_replay_num_ok(NA_real_, NA_real_))
  expect_false(dense_replay_num_ok(5e-13, 0))
  expect_false(dense_replay_num_ok(
    structure(as.integer(dense_mat$phi_matrix), dim = dim(dense_mat$phi_matrix), dimnames = dimnames(dense_mat$phi_matrix)),
    dense_mat$phi_matrix
  ))
  expect_false(dense_replay_num_ok(as.vector(dense_mat$phi_matrix), dense_mat$phi_matrix))
  expect_false(dense_replay_num_ok(
    dense_mat$phi_matrix + 2e-12 * pmax(1, abs(dense_mat$phi_matrix)),
    dense_mat$phi_matrix
  ))
  expect_false(dense_replay_num_ok(-dense_mat$phi_matrix, dense_mat$phi_matrix))
  expect_false(dense_replay_num_ok(
    dense_mat$phi_matrix[c(2L, 1L, 3:nrow(dense_mat$phi_matrix)), , drop = FALSE],
    dense_mat$phi_matrix
  ))

  # Compare chunked IF replay and two-pass covariance across chunk sizes 1, 5, 8, 13, 96
  for (ch_sz in c(1L, 5L, 8L, 13L, 96L)) {
    ch_ctrl <- ipw_point_view_control(ch_sz, 96, 100000, 5000000, 1000000000000, 100000)
    tr <- new.env(parent = emptyenv())
    phi_accum <- matrix(0, nrow = 0L, ncol = 4L)
    comp_accum <- list(
      rate = matrix(0, nrow = 0L, ncol = 4L),
      beta_score = matrix(0, nrow = 0L, ncol = 4L),
      treated = matrix(0, nrow = 0L, ncol = 4L),
      control = matrix(0, nrow = 0L, ncol = 4L)
    )
    replay_ipw_point_if_view(
      plinks,
      control = ch_ctrl,
      include_components = TRUE,
      trace_env = tr,
      callback = function(batch) {
        idx <- batch$a_index
        expect_dense_replay_num(batch$phi, dense_mat$phi_matrix[idx, , drop = FALSE])
        expect_identical(batch$components$rate, dense_mat$components$rate[idx, , drop = FALSE])
        expect_dense_replay_num(batch$components$beta_score, dense_mat$components$beta_score[idx, , drop = FALSE])
        expect_identical(batch$components$treated, dense_mat$components$treated[idx, , drop = FALSE])
        expect_dense_replay_num(batch$components$control, dense_mat$components$control[idx, , drop = FALSE])
        phi_accum <<- rbind(phi_accum, batch$phi)
        for (nm in names(comp_accum)) {
          comp_accum[[nm]] <<- rbind(comp_accum[[nm]], batch$components[[nm]])
        }
      }
    )
    expect_dense_replay_num(phi_accum, dense_mat$phi_matrix)
    expect_identical(comp_accum$rate, dense_mat$components$rate)
    expect_dense_replay_num(comp_accum$beta_score, dense_mat$components$beta_score)
    expect_identical(comp_accum$treated, dense_mat$components$treated)
    expect_dense_replay_num(comp_accum$control, dense_mat$components$control)
    expect_gt(tr$nuisance_cache_hits, 0L)

    dense_cov <- stream_ipw_covariance(fx96$context, dlinks, chunk_size = ch_sz)
    cov_res <- stream_ipw_point_covariance_view(plinks, control = ch_ctrl)
    expect_identical(cov_res$estimates, dense_cov$estimates)
    expect_identical(cov_res$mean_phi, dense_cov$mean_phi)
    expect_identical(cov_res$covariance, dense_cov$covariance)

    # Shared-mask additivity identity: phi(1,2) + phi(2,3) - phi(1,3) == 0
    contrast <- c(1, 1, -1, 0)
    quad <- as.double(t(contrast) %*% cov_res$covariance %*% contrast)
    expect_lt(abs(quad), 1e-14)
  }

  # Character IDs, reversed A-order, and extra singleton/unobserved A-units
  pan_char <- fx96$panel
  pan_char$id <- sprintf("U%03d", pan_char$id)
  pan_extra <- rbind(
    pan_char,
    data.frame(
      id = c("U097_SINGLE", "U098_ALL_NA", "U098_ALL_NA"),
      time = c(1, 1, 2),
      y = c(2.5, NA_real_, NA_real_),
      g = c(4, 7, 7),
      stringsAsFactors = FALSE
    )
  )
  base_ids <- c(sprintf("U%03d", seq_len(96L)), "U097_SINGLE", "U098_ALL_NA")
  a_ids_rev <- rev(base_ids)
  Z_rev <- rbind(fx96$Z, c(1, 0.2), c(1, -0.2))[match(a_ids_rev, base_ids), , drop = FALSE]
  w_rev <- c(fx96$w, 1.1, 0.9)[match(a_ids_rev, base_ids)]

  idx_rev <- cp02_pv_fn("construct_observed_pairs")(pan_extra, "id", "time", "y", "g")
  ctx_rev <- cp02_pv_fn("construct_ipw_context")(idx_rev, id = a_ids_rev, Z = Z_rev, weights = w_rev)
  dlink_rev <- cp02_pv_fn("fit_ipw_link")(ctx_rev, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl)
  dif_rev <- materialize_ipw_if(ctx_rev, list(dlink_rev), max_entries = 10000, chunk_size = 16L, include_components = TRUE)

  po_rev <- construct_ipw_panel_owner(pan_extra, "id", "time", "y", "g", control = own_ctrl)
  fo_rev <- construct_ipw_frame_owner(po_rev, id = a_ids_rev, Z = Z_rev, weights = w_rev, control = own_ctrl)
  vfit_rev <- fit_ipw_nuisance_from_owner(fo_rev, g = 4, K = c(7, Inf), s = 1, t = 2, glm_control = glm_ctrl, control = fv_ctrl)
  plink_rev <- fit_ipw_point_link_from_owner(vfit_rev, s = 1, t = 2, control = pv_ctrl, g = 4, K = c(7, Inf), glm_control = glm_ctrl)

  expect_identical(plink_rev$estimate, dlink_rev$estimate)
  expect_identical(plink_rev$loadings, dlink_rev$loadings)
  na_rev <- length(ctx_rev$a_id)
  pif_rev <- evaluate_ipw_point_if_view_chunk(
    list(plink_rev), a_start = 1L, a_end = na_rev,
    control = ipw_point_view_control(na_rev, na_rev, 100000, 5000000, 1000000000000, 100000),
    include_components = TRUE
  )
  expect_identical(pif_rev$id, dif_rev$id)
  expect_identical(pif_rev$phi, dif_rev$phi_matrix)
  expect_identical(pif_rev$components$rate, dif_rev$components$rate)
  expect_dense_replay_num(pif_rev$components$beta_score, dif_rev$components$beta_score)
  expect_identical(pif_rev$components$treated, dif_rev$components$treated)
  expect_identical(pif_rev$components$control, dif_rev$components$control)
})

test_that("candidate-to-certificate root/link/IF binding and covariance candidate_rows spec match dense reference and certify without refitting", {
  construct_ipw_panel_owner <- cp02_pv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_pv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_pv_fn("construct_ipw_nuisance_view")
  construct_ipw_link_view <- cp02_pv_fn("construct_ipw_link_view")
  fit_ipw_nuisance_view <- cp02_pv_fn("fit_ipw_nuisance_view")
  bind_fitted_nuisance_to_problem_view <- cp02_pv_fn("bind_fitted_nuisance_to_problem_view")
  bind_point_link_to_enclosure_view <- cp02_pv_fn("bind_point_link_to_enclosure_view")
  enclose_point_link_if_view_chunk <- cp02_pv_fn("enclose_point_link_if_view_chunk")
  certify_ipw_point_link_view <- cp02_pv_fn("certify_ipw_point_link_view")
  make_ipw_point_candidate_rows_spec <- cp02_pv_fn("make_ipw_point_candidate_rows_spec")
  ipw_owner_control <- cp02_pv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_pv_fn("ipw_fitted_view_control")
  ipw_point_view_control <- cp02_pv_fn("ipw_point_view_control")
  logit_certificate_control <- cp02_pv_fn("logit_certificate_control")
  ipw_nuisance_problem_control <- cp02_pv_fn("ipw_nuisance_problem_control")
  ipw_link_enclosure_control <- cp02_pv_fn("ipw_link_enclosure_control")
  fit_ipw_link <- cp02_pv_fn("fit_ipw_link")
  materialize_ipw_if <- cp02_pv_fn("materialize_ipw_if")
  construct_certified_ipw_nuisance_problem_view <- cp02_pv_fn("construct_certified_ipw_nuisance_problem_view")
  enclose_ipw_link_view <- cp02_pv_fn("enclose_ipw_link_view")
  enclose_ipw_link_if_view_chunk <- cp02_pv_fn("enclose_ipw_link_if_view_chunk")

  # Exact weighted96 setup from context/fitted_root_diagnostic.R
  fx <- cp02_pv_weighted96_fixture(weighted = TRUE)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)
  own_ctrl <- ipw_owner_control(8, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(8, 96, 100000, 500000, 5000000, 100000)
  pv_ctrl <- ipw_point_view_control(8, 96, 100000, 5000000, 1000000000000, 100000)
  cc <- logit_certificate_control(18, 40, 16, 8, 500000, 100000000)
  pc <- ipw_nuisance_problem_control(8, 100000, 500000, 100000000)
  lc <- ipw_link_enclosure_control(8, 100000, 500000, 100000000, 100000)

  po <- construct_ipw_panel_owner(fx$panel, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx$n), Z = fx$Z, weights = fx$w, control = own_ctrl)
  nv12 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl, nuisance_id = "actual_weighted96_pair12")
  vfit12 <- fit_ipw_nuisance_view(nv12, glm_control = glm_ctrl, control = fv_ctrl)
  F_mat <- unname(vfit12$factors$H_chol)
  r_val <- 1e-6

  # Direct reference via dense fit_ipw_link + materialize_ipw_if
  dlink12 <- fit_ipw_link(fx$context, g = 4, K = c(7, Inf), s = 1, t = 2, control = glm_ctrl,
                          link_id = "actual_weighted96_pair12", functional_key = "A")
  lv12 <- construct_ipw_link_view(
    fo, s = 1, t = 2, nuisance_view = nv12, control = own_ctrl,
    loading_beta_candidate = as.double(dlink12$loadings$loading_beta),
    link_id = "actual_weighted96_pair12", functional_key = "A"
  )
  ref_prob <- construct_certified_ipw_nuisance_problem_view(
    nv12, as.double(vfit12$coefficients), F_mat, r_val, cc, pc, own_ctrl
  )
  ref_enc12 <- enclose_ipw_link_view(
    ref_prob, lv12, as.double(dlink12$loadings$loading_beta), lc, own_ctrl
  )
  dense_if <- materialize_ipw_if(fx$context, dlink12, 10000, 8, FALSE)

  # Compact binding via R/chained_ipw_point_view.R (zero dense fit_ipw_link or materialize_ipw_if calls)
  tr <- new.env(parent = emptyenv())
  cert12 <- certify_ipw_point_link_view(
    fitted_view = vfit12, link_view = lv12, F = F_mat, r = r_val,
    point_control = pv_ctrl, certificate_control = cc,
    problem_control = pc, link_control = lc,
    reader_control = own_ctrl, enclose_if_chunks = TRUE,
    include_if_components = FALSE, include_shared_xi = FALSE,
    trace_env = tr
  )
  expect_true(cert12$certified)
  expect_identical(cert12$status, "CERTIFIED")
  expect_identical(tr$propensity_solver_calls, 0L)
  expect_identical(cert12$problem$raw$beta0, ref_prob$raw$beta0)
  expect_identical(cert12$problem$certificate, ref_prob$certificate)
  expect_identical(cert12$link_enclosure$enclosure, ref_enc12$enclosure)
  expect_identical(cert12$if_summary$n_if_rows_certified, 96L)
  expect_false(cert12$status_separation$point_candidate_is_certificate)

  # Verify chunk-by-chunk IF enclosure parity with dense materialize_ipw_if
  ch_if <- enclose_point_link_if_view_chunk(
    link_enclosure = cert12$link_enclosure, point_link = cert12$point_link,
    a_start = 1L, a_end = 8L, point_control = pv_ctrl,
    link_control = lc, reader_control = own_ctrl,
    include_components = FALSE, include_shared_xi = FALSE
  )
  ref_if <- enclose_ipw_link_if_view_chunk(
    ref_enc12, 1L, 8L, lc, own_ctrl,
    phi_candidate = as.double(dense_if$phi_matrix[1:8, 1L]),
    H_chol_candidate = unname(vfit12$factors$H_chol),
    include_components = FALSE, include_shared_xi = FALSE
  )
  expect_identical(ch_if$certified, ref_if$certified)
  expect_identical(ch_if$phi, ref_if$phi)
  expect_identical(ch_if$separated_row_errors, ref_if$separated_row_errors)

  # Verify make_ipw_point_candidate_rows_spec callback schema and exact batch values
  lv23 <- construct_ipw_link_view(
    fo, s = 2, t = 3, nuisance_view = nv12, control = own_ctrl,
    link_id = "actual_weighted96_pair23", functional_key = "B"
  )
  plink23 <- cp02_pv_fn("fit_ipw_point_link_view")(vfit12, lv23, control = pv_ctrl)
  cand_spec <- make_ipw_point_candidate_rows_spec(
    list(cert12$point_link, plink23),
    control = pv_ctrl
  )
  expect_s3_class(cand_spec, "cdid_ipw_candidate_rows_spec")
  cb_batch <- cand_spec$callback(9L, 16L)
  expect_false(is.object(cb_batch))
  expect_identical(names(cb_batch), c("phi", "id", "link_ids", "a_index"))
  expect_identical(cb_batch$a_index, 9:16)
  expect_identical(cb_batch$link_ids, c("actual_weighted96_pair12", "actual_weighted96_pair23"))
  expect_identical(cb_batch$phi[, 1L], dense_if$phi_matrix[9:16, 1L])
})

test_that("exact binary64 point-link RDS snapshots require and pass seal/ownership revalidation and refuse tampered values", {
  construct_ipw_panel_owner <- cp02_pv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_pv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_pv_fn("construct_ipw_nuisance_view")
  construct_ipw_link_view <- cp02_pv_fn("construct_ipw_link_view")
  fit_ipw_nuisance_view <- cp02_pv_fn("fit_ipw_nuisance_view")
  fit_ipw_point_link_view <- cp02_pv_fn("fit_ipw_point_link_view")
  validate_ipw_point_link_view <- cp02_pv_fn("validate_ipw_point_link_view")
  extract_ipw_point_link_snapshot <- cp02_pv_fn("extract_ipw_point_link_snapshot")
  revalidate_ipw_point_link_view <- cp02_pv_fn("revalidate_ipw_point_link_view")
  ipw_owner_control <- cp02_pv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_pv_fn("ipw_fitted_view_control")
  ipw_point_view_control <- cp02_pv_fn("ipw_point_view_control")

  fx96 <- cp02_pv_weighted96_fixture(weighted = TRUE)
  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
  pv_ctrl <- ipw_point_view_control(16, 96, 100000, 5000000, 1000000000000, 100000)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  po <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)
  nv12 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  lv12 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, nuisance_view = nv12, control = own_ctrl)
  lv23 <- construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 2, t = 3, nuisance_view = nv12, control = own_ctrl)
  vfit12 <- fit_ipw_nuisance_view(nv12, glm_control = glm_ctrl, control = fv_ctrl)
  plink12 <- fit_ipw_point_link_view(vfit12, lv12, control = pv_ctrl)
  snap12 <- extract_ipw_point_link_snapshot(plink12)

  rds_view <- tempfile(fileext = ".rds")
  rds_snap <- tempfile(fileext = ".rds")
  on.exit(unlink(c(rds_view, rds_snap)), add = TRUE)

  saveRDS(plink12, rds_view)
  saveRDS(snap12, rds_snap)

  loaded_view <- readRDS(rds_view)
  loaded_snap <- readRDS(rds_snap)

  expect_identical(loaded_view$estimate, plink12$estimate)
  expect_identical(loaded_view$partials, plink12$partials)
  expect_identical(loaded_view$loadings, plink12$loadings)
  expect_identical(loaded_snap$estimate, plink12$estimate)

  # Raw deserialized objects must be refused prior to revalidation
  expect_error(
    validate_ipw_point_link_view(loaded_view),
    class = "cdid_ipw_point_view_integrity_error"
  )
  expect_error(
    validate_ipw_point_link_view(loaded_snap),
    class = "cdid_ipw_point_view_integrity_error"
  )

  # Revalidation against live vfit12 + lv12 succeeds and reseals
  rev_view <- revalidate_ipw_point_link_view(loaded_view, vfit12, lv12, control = pv_ctrl)
  rev_snap <- revalidate_ipw_point_link_view(loaded_snap, vfit12, lv12, control = pv_ctrl)
  validate_ipw_point_link_view(rev_view, fitted_view = vfit12, link_view = lv12, frame_owner = fo)
  validate_ipw_point_link_view(rev_snap, fitted_view = vfit12, link_view = lv12, frame_owner = fo)
  expect_identical(rev_view$estimate, plink12$estimate)
  expect_identical(rev_snap$loadings, plink12$loadings)

  # Tampered loadings or wrong link_view must be refused
  bad_snap <- loaded_snap
  bad_snap$loadings$loading_rate[1L] <- bad_snap$loadings$loading_rate[1L] + 1e-12
  expect_error(
    revalidate_ipw_point_link_view(bad_snap, vfit12, lv12, control = pv_ctrl),
    class = "cdid_ipw_point_view_revalidation_error"
  )
  expect_error(
    revalidate_ipw_point_link_view(loaded_snap, vfit12, lv23, control = pv_ctrl),
    class = "cdid_ipw_point_view_revalidation_error"
  )
})

test_that("discriminating source-fault / formula mutants and ownership/sharing/cap/lifetime guards fail loudly", {
  local <- pv_local_api()
  fx96 <- cp02_pv_weighted96_fixture(weighted = TRUE)
  oc <- local$ipw_owner_control(11, 100000, 500000, 5000000, 100000)
  fc <- local$ipw_fitted_view_control(11, 96, 100000, 500000, 5000000, 100000)
  pc <- local$ipw_point_view_control(11, 96, 100000, 5000000, 1000000000000, 100000)
  gc <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  po <- local$construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = oc)
  fo <- local$construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = oc)
  nv12 <- local$construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = oc)
  nv78 <- local$construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 7, t = 8, control = oc)
  lv12 <- local$construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, nuisance_view = nv12, control = oc)
  lv23 <- local$construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 2, t = 3, nuisance_view = nv12, control = oc)
  lv78 <- local$construct_ipw_link_view(fo, g = 4, K = c(7, Inf), s = 7, t = 8, nuisance_view = nv78, control = oc)

  vf12 <- local$fit_ipw_nuisance_view(nv12, glm_control = gc, control = fc)
  vf78 <- local$fit_ipw_nuisance_view(nv78, glm_control = gc, control = fc)

  # 1. Lifetime check: zero local query slice aliases during query_ipw_view_chunk and callbacks
  orig_q <- local$query_ipw_view_chunk
  q_aliases <- integer()
  local$query_ipw_view_chunk <- function(...) {
    q_aliases <<- c(
      q_aliases,
      length(intersect(c("a_ch", "gi_ch", "loc_ch", "w_ch", "Z_ch", "dy_ch", "dst"), ls(parent.frame())))
    )
    orig_q(...)
  }
  cb_aliases <- integer()
  pl12 <- local$fit_ipw_point_link_view(vf12, lv12, control = pc, callback = function(ev) {
    cb_aliases <<- c(
      cb_aliases,
      length(intersect(c("a_ch", "gi_ch", "loc_ch", "w_ch", "Z_ch", "dy_ch", "dst"), ls(parent.frame())))
    )
  })
  pl23 <- local$fit_ipw_point_link_view(vf12, lv23, control = pc)
  expect_gt(length(q_aliases), 1L)
  expect_true(all(q_aliases == 0L))
  expect_true(all(cb_aliases == 0L))

  if_cb_aliases <- integer()
  local$replay_ipw_point_if_view(list(pl12, pl23), control = pc, callback = function(batch) {
    if_cb_aliases <<- c(
      if_cb_aliases,
      length(intersect(c("a_ch", "gi_ch", "loc_ch", "w_ch", "Z_ch", "dy_ch"), ls(parent.frame())))
    )
  })
  expect_gt(length(if_cb_aliases), 1L)
  expect_true(all(if_cb_aliases == 0L))

  # 2. Early preflight refusal before any query_ipw_view_chunk or callback execution
  q_calls <- 0L
  local$query_ipw_view_chunk <- function(...) {
    q_calls <<- q_calls + 1L
    orig_q(...)
  }
  low_ret <- pc
  low_ret$max_retained_entries <- 5
  expect_error(local$fit_ipw_point_link_view(vf12, lv12, control = low_ret), class = "cdid_ipw_owner_bound_error")
  expect_identical(q_calls, 0L)

  low_work <- pc
  low_work$max_work_ops <- 5
  expect_error(local$replay_ipw_point_if_view(list(pl12, pl23), control = low_work, callback = function(b) NULL), class = "cdid_ipw_owner_bound_error")
  expect_identical(q_calls, 0L)

  low_out <- pc
  low_out$max_output_entries <- 2
  expect_error(local$evaluate_ipw_point_if_view_chunk(list(pl12, pl23), 1L, 11L, control = low_out), class = "cdid_ipw_owner_bound_error")
  expect_identical(q_calls, 0L)

  # 3. Sharing and ownership mismatches refused before outcome endpoint reads
  expect_error(
    local$fit_ipw_point_link_view(vf12, lv78, control = pc),
    class = "cdid_ipw_shared_mask_error"
  )
  expect_identical(q_calls, 0L)

  # Equal-count swapped-observer mask refusal (counts and masses match, but observer indices differ)
  local$query_ipw_view_chunk <- orig_q
  pan_swap <- fx96$panel
  pan_swap$y[pan_swap$id == 1L & pan_swap$time == 3] <- NA_real_
  pan_swap$y[pan_swap$id == 2L & pan_swap$time == 2] <- NA_real_
  w_swap <- rep(1.0, fx96$n)
  Z_swap <- fx96$Z
  Z_swap[2L, ] <- Z_swap[1L, ]
  po_sw <- local$construct_ipw_panel_owner(pan_swap, "id", "time", "y", "g", control = oc)
  fo_sw <- local$construct_ipw_frame_owner(po_sw, id = seq_len(fx96$n), Z = Z_swap, weights = w_swap, control = oc)
  nv_sw12 <- local$construct_ipw_nuisance_view(fo_sw, g = 4, K = c(7, Inf), s = 1, t = 2, control = oc)
  lv_sw13 <- local$construct_ipw_link_view(fo_sw, g = 4, K = c(7, Inf), s = 1, t = 3, control = oc)
  vf_sw12 <- local$fit_ipw_nuisance_view(nv_sw12, glm_control = gc, control = fc)
  outcome_q_calls <- 0L
  local$query_ipw_view_chunk <- function(..., read_outcomes = FALSE) {
    if (isTRUE(read_outcomes)) outcome_q_calls <<- outcome_q_calls + 1L
    orig_q(..., read_outcomes = read_outcomes)
  }
  expect_error(
    local$fit_ipw_point_link_view(vf_sw12, lv_sw13, control = pc),
    class = "cdid_ipw_shared_mask_error"
  )
  expect_identical(outcome_q_calls, 0L)
  q_calls <- 0L
  local$query_ipw_view_chunk <- function(...) {
    q_calls <<- q_calls + 1L
    orig_q(...)
  }

  po2 <- local$construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = oc)
  fo2 <- local$construct_ipw_frame_owner(po2, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = oc)
  nv12_foreign <- local$construct_ipw_nuisance_view(fo2, g = 4, K = c(7, Inf), s = 1, t = 2, control = oc)
  lv12_foreign <- local$construct_ipw_link_view(fo2, g = 4, K = c(7, Inf), s = 1, t = 2, nuisance_view = nv12_foreign, control = oc)
  vf12_foreign <- local$fit_ipw_nuisance_view(nv12_foreign, glm_control = gc, control = fc)
  q_calls <- 0L
  expect_error(
    local$fit_ipw_point_link_view(vf12, lv12_foreign, control = pc),
    class = "cdid_ipw_owner_mismatch_error"
  )
  expect_identical(q_calls, 0L)

  local$query_ipw_view_chunk <- orig_q
  pl12_foreign <- local$fit_ipw_point_link_view(vf12_foreign, lv12_foreign, control = pc)
  q_calls <- 0L
  local$query_ipw_view_chunk <- function(...) {
    q_calls <<- q_calls + 1L
    orig_q(...)
  }
  expect_error(
    local$replay_ipw_point_if_view(list(pl12, pl12_foreign), control = pc, callback = function(b) NULL),
    class = "cdid_ipw_owner_mismatch_error"
  )
  expect_identical(q_calls, 0L)
  local$query_ipw_view_chunk <- orig_q

  # 4. Discriminating formula mutants for all 4 complete-IF terms and nonobserver rates
  full_if <- local$evaluate_ipw_point_if_view_chunk(
    list(pl12), 1L, 96L,
    control = local$ipw_point_view_control(96, 96, 100000, 5000000, 1000000000000, 100000),
    include_components = TRUE
  )
  phi_exact <- full_if$phi[, 1L]
  c_tr <- full_if$components$treated[, 1L]
  c_ct <- full_if$components$control[, 1L]
  c_bs <- full_if$components$beta_score[, 1L]
  c_rt <- full_if$components$rate[, 1L]

  # Omitting or sign-flipping any of the 4 terms produces large separation
  expect_gt(max(abs(phi_exact - (c_ct + c_bs + c_rt))), 1e-1)
  expect_gt(max(abs(phi_exact - (c_tr + c_bs + c_rt))), 1e-1)
  expect_gt(max(abs(phi_exact - (c_tr + c_ct + c_rt))), 1e-1)
  expect_gt(max(abs(phi_exact - (c_tr + c_ct + c_bs))), 1e-1)
  expect_gt(max(abs(phi_exact - (c_tr - c_ct + c_bs + c_rt))), 1e-1)
  expect_gt(max(abs(phi_exact - (c_tr + c_ct - c_bs + c_rt))), 1e-1)
  expect_gt(max(abs(phi_exact - (c_tr + c_ct + c_bs - c_rt))), 1e-1)

  # Omitting nonobserver rate contribution (setting rate = 0 when !in_o)
  in_risk <- fx96$G %in% c(4, 7, Inf)
  in_o <- in_risk & fx96$obs[, 1] & fx96$obs[, 2]
  mut_no_nonobs_rate <- ifelse(in_o, phi_exact, 0)
  expect_gt(max(abs(phi_exact - mut_no_nonobs_rate)), 1e-1)
})

test_that("repaired byte/work gates run before validators, outcomes, and callbacks", {
  e <- pv_local_api()
  fx <- cp02_pv_weighted96_fixture(TRUE)
  gc <- stats::glm.control(epsilon=1e-12,maxit=100,trace=FALSE)
  oc <- e$ipw_owner_control(8,100000,500000,5000000,100000)
  fc <- e$ipw_fitted_view_control(8,96,100000,500000,5000000,100000)
  # CP02_ACCOUNTING_DIAGNOSTIC_1E12_5E6 is distinct from preserved 5e6 work inputs.
  diagnostic_pc <- e$ipw_point_view_control(8,96,1000000,5000000,1e12,1000000)
  po <- e$construct_ipw_panel_owner(fx$panel,"id","time","y","g",control=oc)
  fo <- e$construct_ipw_frame_owner(po,1:96,fx$Z,fx$w,control=oc)
  nv <- e$construct_ipw_nuisance_view(fo,4,c(7,Inf),1,2,control=oc)
  fv <- e$fit_ipw_nuisance_view(nv,gc,fc)
  lv <- e$construct_ipw_link_view(fo,1,2,nuisance_view=nv,control=oc)
  reads <- callbacks <- deep <- 0L
  orig_q <- e$query_ipw_view_chunk
  orig_v <- e$validate_ipw_fitted_nuisance_view
  e$query_ipw_view_chunk <- function(...) {reads <<- reads+1L; orig_q(...)}
  e$validate_ipw_fitted_nuisance_view <- function(...) {deep <<- deep+1L; orig_v(...)}
  low <- diagnostic_pc; low$max_work_ops <- 1
  expect_error(e$preflight_ipw_point_link_view(fv,lv,low),class="cdid_ipw_owner_bound_error")
  expect_identical(c(reads,deep),c(0L,0L))
  huge <- paste(rep("a",300000),collapse="")
  bounded <- diagnostic_pc; bounded$max_output_entries <- 100000
  expect_error(e$fit_ipw_point_link_view(fv,lv,bounded,link_id=huge,
    callback=function(ev) callbacks <<- callbacks+1L),class="cdid_ipw_owner_bound_error")
  expect_identical(c(reads,callbacks),c(0L,0L))
  links <- lapply(seq_len(40),function(j) e$fit_ipw_point_link_view(fv,lv,diagnostic_pc,link_id=paste0("Q",j)))
  reads <- deep <- 0L
  low <- diagnostic_pc; low$max_work_ops <- 16*40^2-1
  expect_error(e$preflight_ipw_point_if_view(links,low),class="cdid_ipw_owner_bound_error")
  expect_identical(c(reads,deep),c(0L,0L))
  pre <- e$preflight_ipw_point_if_view(links,diagnostic_pc)
  expect_gte(pre$normalization_work_bound,16*40^2)
  expect_identical(pre$borrowed_owner_entries,
    po$accounting$actual_retained_entries+fo$accounting$retained_matrix_entries+fo$accounting$retained_vector_entries)
  reads <- callbacks <- 0L
  low <- diagnostic_pc; low$max_work_ops <- pre$total_work_ops_bound-1
  expect_error(e$replay_ipw_point_if_view(links,low,function(b) callbacks <<- callbacks+1L),
    class="cdid_ipw_owner_bound_error")
  expect_identical(c(reads,callbacks),c(0L,0L))
  # Preserve the truthful original 5e6-cap refusal rather than changing its predicate.
  original_pc <- e$ipw_point_view_control(8,96,100000,500000,5000000,100000)
  expect_error(e$make_ipw_point_candidate_rows_spec(list(links[[1]]),original_pc),
    class="cdid_ipw_owner_bound_error")
  expect_identical(reads,0L)
})

test_that("reachable point, callback, trace, query, and cache graphs respect lifetimes", {
  e <- pv_local_api()
  fx <- cp02_pv_weighted96_fixture(TRUE)
  gc <- stats::glm.control(epsilon=1e-12,maxit=100,trace=FALSE)
  oc <- e$ipw_owner_control(8,100000,500000,5000000,100000)
  fc <- e$ipw_fitted_view_control(8,96,100000,500000,5000000,100000)
  diagnostic_pc <- e$ipw_point_view_control(8,96,1000000,5000000,1e12,1000000)
  po <- e$construct_ipw_panel_owner(fx$panel,"id","time","y","g",control=oc)
  fo <- e$construct_ipw_frame_owner(po,1:96,fx$Z,fx$w,control=oc)
  nv <- e$construct_ipw_nuisance_view(fo,4,c(7,Inf),1,2,control=oc)
  fv <- e$fit_ipw_nuisance_view(nv,gc,fc)
  lv <- e$construct_ipw_link_view(fo,1,2,nuisance_view=nv,control=oc,link_id="L12")
  lv23 <- e$construct_ipw_link_view(fo,2,3,nuisance_view=nv,control=oc,link_id="L23")
  nv78 <- e$construct_ipw_nuisance_view(fo,4,c(7,Inf),7,8,control=oc)
  fv78 <- e$fit_ipw_nuisance_view(nv78,gc,fc)
  lv78 <- e$construct_ipw_link_view(fo,7,8,nuisance_view=nv78,control=oc,link_id="L78")
  query_entries <- 0L
  orig_q <- e$query_ipw_view_chunk
  e$query_ipw_view_chunk <- function(...) {
    fr <- parent.frame()
    if (exists("y_trace",fr,inherits=FALSE)) {
      tr <- get("y_trace",fr)
      stopifnot(length(tr$pos_s_index)==0L,length(tr$pos_t_index)==0L)
      if (exists("ch",fr,inherits=FALSE)) stopifnot(is.null(get("ch",fr)))
      # Detect query-shaped aliases under any name, and prior replay batches
      # in the actual producer stack, without trusting release counters.
      for (nm in ls(fr,all.names=TRUE)) {
        x <- get(nm,fr,inherits=FALSE)
        if (is.list(x)) stopifnot(!all(c("observed_risk","read_diagnostics") %in% names(x)))
      }
      for (fr0 in sys.frames()) if (exists("emitted_units",fr0,inherits=FALSE) &&
          exists("batch",fr0,inherits=FALSE)) {
        for (nm in ls(fr0,all.names=TRUE)) {
          x <- get(nm,fr0,inherits=FALSE)
          if (is.list(x)) stopifnot(!all(c("phi","a_index") %in% names(x)))
        }
      }
      for (nm in c("a_ch","dy_ch","Z_ch","pos_k","dy_k","Z_k","col_phi"))
        stopifnot(!exists(nm,fr,inherits=FALSE))
      if (exists("can_reuse_cache",fr,inherits=FALSE) && !get("can_reuse_cache",fr)) {
        ce <- get("nu_cache_env",fr)
        stopifnot(identical(ls(ce),"active"),!ce$active)
      }
      query_entries <<- query_entries+1L
    }
    orig_q(...)
  }
  callback_entries <- 0L
  pl <- e$fit_ipw_point_link_view(fv,lv,diagnostic_pc,callback=function(ev) {
    fr <- parent.frame(); tr <- get("y_trace",fr)
    expect_length(tr$pos_s_index,0L); expect_length(tr$pos_t_index,0L)
    expect_null(get("ch",fr))
    for (nm in c("a_ch","dy_ch","Z_ch")) expect_false(exists(nm,fr,inherits=FALSE))
    if (ev$phase=="after_observer_release") {
      oe <- get("obs_buf_env",fr)
      expect_false(oe$active)
      for (nm in c("risk_a_idx","risk_gi","risk_w","risk_dy","Z_use","st","lv","ow","means","part"))
        expect_false(exists(nm,oe,inherits=FALSE))
    }
    callback_entries <<- callback_entries+1L
  })
  pl23 <- e$fit_ipw_point_link_view(fv,lv23,diagnostic_pc)
  pl78 <- e$fit_ipw_point_link_view(fv78,lv78,diagnostic_pc)
  expect_gt(callback_entries,1L)
  expect_lte(e$ipw_point_view_payload_budget(pl,diagnostic_pc)$entries,pl$accounting$total_retained_entries)
  trace <- new.env(parent=emptyenv())
  spec <- e$make_ipw_point_candidate_rows_spec(list(pl,pl23,pl78),diagnostic_pc,trace_env=trace)
  # Follow all non-owner list/function/environment edges, including sealed cycles.
  # Only exact shared base-owner values and the code environment are exempted.
  aliases <- cp02_pv_fn("gmm_sexp_aliases")
  seen <- list()
  walk <- function(x) {
    if (is.null(x) || aliases(x,fo) || aliases(x,po)) return(invisible(NULL))
    if (is.function(x)) return(walk(environment(x)))
    if (is.environment(x)) {
      if (identical(x,e) || identical(x,api) || identical(x,emptyenv())) return(invisible(NULL))
      if (any(vapply(seen,identical,logical(1),y=x))) return(invisible(NULL))
      seen[[length(seen)+1L]] <<- x
      for (nm in ls(x,all.names=TRUE)) walk(get(nm,x,inherits=FALSE))
      return(walk(parent.env(x)))
    }
    if (is.list(x)) {
      for (nm in intersect(names(x),c("risk_a_idx","risk_gi","risk_w","risk_dy","Z_use",
        "pos_s_index","pos_t_index","phi_matrix","a_k","Z_k","odds_k")))
        stopifnot(length(x[[nm]])==0L)
      for (j in seq_along(x)) walk(x[[j]])
    } else if (is.numeric(x) || is.logical(x)) {
      if (is.null(dim(x))) stopifnot(length(x)<=3L) else stopifnot(all(dim(x)<=3L))
    }
    invisible(NULL)
  }
  ce <- environment(spec$callback)
  expect_setequal(ls(ce,all.names=TRUE),c("norm","ctrl","use_nuisance_cache","trace_env"))
  expect_identical(parent.env(ce),e)
  expect_identical(substitute(use_nuisance_cache,ce),TRUE)
  expect_identical(substitute(trace_env,ce),trace)
  walk(pl); walk(spec); walk(trace)
  old <- NULL
  e$replay_ipw_point_if_view(list(pl,pl23,pl78),diagnostic_pc,function(b) {
    fr <- parent.frame()
    expect_identical(get("batch",fr),b)
    expect_false(exists("prior_batch",fr,inherits=FALSE))
    expect_identical(b$a_index,seq.int(min(b$a_index),max(b$a_index)))
    old <<- max(b$a_index)
  },include_components=TRUE,trace_env=trace)
  expect_identical(old,96L)
  walk(spec); walk(trace)
  expect_gt(query_entries,36L)
  expect_identical(trace$active_query_chunks,0L)
  expect_identical(trace$active_nuisance_caches,0L)
  expect_gt(trace$nuisance_cache_hits,0L)
})

test_that("whole certification budgets and exact signed-zero witness binding precede outcomes", {
  e <- pv_local_api()
  fx <- cp02_pv_weighted96_fixture(TRUE)
  gc <- stats::glm.control(epsilon=1e-12,maxit=100,trace=FALSE)
  oc <- e$ipw_owner_control(8,100000,500000,5000000,100000)
  fc <- e$ipw_fitted_view_control(8,96,100000,500000,5000000,100000)
  diagnostic_pc <- e$ipw_point_view_control(8,96,100000,5000000,1e12,100000)
  cc <- cp02_pv_fn("logit_certificate_control")(18,40,16,8,500000,100000000)
  nc <- cp02_pv_fn("ipw_nuisance_problem_control")(8,100000,500000,100000000)
  lc <- cp02_pv_fn("ipw_link_enclosure_control")(8,100000,500000,100000000,100000)
  po <- e$construct_ipw_panel_owner(fx$panel,"id","time","y","g",control=oc)
  fo <- e$construct_ipw_frame_owner(po,1:96,fx$Z,fx$w,control=oc)
  nv <- e$construct_ipw_nuisance_view(fo,4,c(7,Inf),1,2,control=oc)
  fv <- e$fit_ipw_nuisance_view(nv,gc,fc)
  lv <- e$construct_ipw_link_view(fo,1,2,nuisance_view=nv,control=oc,link_id="L12")
  F <- unname(fv$factors$H_chol); r <- 1e-6
  prob <- e$bind_fitted_nuisance_to_problem_view(fv,F,r,cc,nc,oc)
  reads <- 0L; orig_q <- e$query_ipw_view_chunk
  e$query_ipw_view_chunk <- function(...) {reads <<- reads+1L; orig_q(...)}
  call_cert <- function(pc=diagnostic_pc,F_arg=F) e$certify_ipw_point_link_view(
    fv,lv,F_arg,r,pc,cc,nc,lc,oc,problem=prob)
  pre <- e$preflight_ipw_point_certificate_binding(fv,lv,diagnostic_pc,cc,nc,lc,oc,
    TRUE,FALSE,FALSE,F,r,prob)
  low <- diagnostic_pc; low$max_work_ops <- pre$total_work_ops_bound-1
  expect_error(call_cert(low),class="cdid_ipw_owner_bound_error")
  expect_identical(reads,0L)
  low <- diagnostic_pc; low$max_output_entries <- pre$retained_entries-1
  expect_error(call_cert(low),class="cdid_ipw_owner_bound_error")
  expect_identical(reads,0L)
  expect_equal(pre$n_if_chunks,12)
  expect_gte(pre$summary_entries,12*1024)
  F_bad <- F; F_bad[2,1] <- -0
  expect_true(identical(F,F_bad))
  expect_false(identical(F,F_bad,num.eq=FALSE))
  expect_error(call_cert(F_arg=F_bad),class="cdid_ipw_point_view_candidate_mismatch_error")
  expect_identical(reads,0L)
  pl <- e$fit_ipw_point_link_view(fv,lv,diagnostic_pc)
  reads <- 0L
  expect_error(e$evaluate_ipw_point_if_view_chunk(list(pl),1.5,8,diagnostic_pc),class="cdid_ipw_owner_bound_error")
  norm <- e$normalize_ipw_point_links_view(list(pl),control=diagnostic_pc)
  expect_error(e$ipw_point_if_eval_chunk_internal(norm,1,8,diagnostic_pc,batch_index=1.5),class="cdid_ipw_owner_bound_error")
  expect_error(e$evaluate_ipw_point_if_view_chunk(list(pl),1,8,diagnostic_pc,include_components=1),class="cdid_ipw_point_view_input_error")
  expect_identical(reads,0L)
})

test_that("reference-search work is independently gated before reference comparisons", {
  e <- pv_local_api()
  fx <- cp02_pv_weighted96_fixture(TRUE)
  gc <- stats::glm.control(epsilon=1e-12,maxit=100,trace=FALSE)
  oc <- e$ipw_owner_control(8,100000,500000,5000000,100000)
  fc <- e$ipw_fitted_view_control(8,96,100000,500000,5000000,100000)
  diagnostic_pc <- e$ipw_point_view_control(8,96,1000000,5000000,1e12,1000000)
  po <- e$construct_ipw_panel_owner(fx$panel,"id","time","y","g",control=oc)
  fo <- e$construct_ipw_frame_owner(po,1:96,fx$Z,fx$w,control=oc)
  nv <- e$construct_ipw_nuisance_view(fo,4,c(7,Inf),1,2,control=oc)
  fv <- e$fit_ipw_nuisance_view(nv,gc,fc)
  lv <- e$construct_ipw_link_view(fo,1,2,nuisance_view=nv,control=oc)
  links <- lapply(seq_len(40),function(j) e$fit_ipw_point_link_view(fv,lv,diagnostic_pc,link_id=paste0("R",j)))
  # Isolate search cost from previously tested byte/metadata gates. These exact
  # sealed links were fully validated above. No production helper is changed.
  e$ipw_point_view_payload_budget <- function(...) list(entries=0,work=0,max_bytes=0)
  ref_comparisons <- 0L
  orig_equal <- e$ipw_point_view_identical
  e$ipw_point_view_identical <- function(x,y) {
    if (is.environment(x) && is.environment(y)) ref_comparisons <<- ref_comparisons+1L
    orig_equal(x,y)
  }
  low <- diagnostic_pc; low$max_work_ops <- 16*40^2-1
  expect_error(e$normalize_ipw_point_links_view(links,control=low),class="cdid_ipw_owner_bound_error")
  expect_identical(ref_comparisons,0L)
  # Above the gate, actually execute the distinct-reference search, not only a ledger.
  e$normalize_ipw_point_links_view(links,control=diagnostic_pc)
  expect_gt(ref_comparisons,40L)
})

test_that("from-owner construction gates all posting scans and oversized override metadata", {
  e <- pv_local_api()
  fx <- cp02_pv_weighted96_fixture(TRUE)
  oc <- e$ipw_owner_control(8,100000,500000,5000000,100000)
  fc <- e$ipw_fitted_view_control(8,96,100000,500000,5000000,100000)
  diagnostic_pc <- e$ipw_point_view_control(8,96,100000,5000000,1e12,100000)
  po <- e$construct_ipw_panel_owner(fx$panel,"id","time","y","g",control=oc)
  fo <- e$construct_ipw_frame_owner(po,1:96,fx$Z,fx$w,control=oc)
  nv <- e$construct_ipw_nuisance_view(fo,4,c(7,Inf),1,2,control=oc)
  fv <- e$fit_ipw_nuisance_view(nv,stats::glm.control(epsilon=1e-12,maxit=100,trace=FALSE),fc)
  calls <- 0L; orig <- e$construct_ipw_link_view
  e$construct_ipw_link_view <- function(...) {calls <<- calls+1L; orig(...)}
  low <- diagnostic_pc; low$max_observer_units <- 1
  expect_error(e$fit_ipw_point_link_from_owner(fv,1,2,low),class="cdid_ipw_owner_bound_error")
  expect_identical(calls,0L)
  low <- diagnostic_pc; low$max_work_ops <- 1
  expect_error(e$fit_ipw_point_link_from_owner(fv,1,2,low),class="cdid_ipw_owner_bound_error")
  expect_identical(calls,0L)
  huge <- paste(rep("q",300000),collapse="")
  expect_error(e$fit_ipw_point_link_from_owner(fv,1,2,diagnostic_pc,link_id=huge),class="cdid_ipw_owner_bound_error")
  expect_identical(calls,0L)
})
