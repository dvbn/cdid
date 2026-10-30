# Focused tests for the private CP02 owner-backed finite weighted-logit
# nuisance fitting slice (`R/chained_ipw_fitted_view.R` and the shared
# observer-buffer worker in `R/chained_ipw_link.R`).

if(!requireNamespace("cdid",quietly=TRUE)) stop("Actual installed cdid required; no source fallback.",call.=FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api),environmentIsLocked(api))
cp02_fv_fn <- function(name) {
 f<-get(name,envir=api,inherits=FALSE)
 if(!is.function(f) || !identical(environment(f),api)) stop("Foreign function binding: ",name,call.=FALSE)
 f
}
cp02_fv_locate_file <- function(...) {
 parts<-c(...)
 if(identical(parts,c("reference","panel24.csv"))) return(testthat::test_path("fixtures","fitted_panel24.csv"))
 if(identical(parts,c("tests","testthat","reference","cp02_shared_ipw_audit_newton_reference.R")))
   return(testthat::test_path("reference","cp02_shared_ipw_audit_newton_reference.R"))
 stop("Unlisted test-local reference path.",call.=FALSE)
}
cp02_fv_panel24_fixture <- function(weighted = TRUE) {
  p24_path <- cp02_fv_locate_file("reference", "panel24.csv")
  p24 <- utils::read.csv(p24_path, stringsAsFactors = FALSE)
  u24 <- p24[!duplicated(p24$id), ]
  Z24 <- cbind(one = 1, sign = as.double(u24$sign))
  w24 <- if (weighted) 0.8 + ((seq_len(24L) * 7L) %% 11L) / 20 else rep(1, 24L)
  idx24 <- cp02_fv_fn("construct_observed_pairs")(p24, "id", "time", "y", "g")
  ctx24 <- cp02_fv_fn("construct_ipw_context")(idx24, id = u24$id, Z = Z24, weights = w24)
  list(panel = p24, units = u24, Z = Z24, w = w24, index = idx24, context = ctx24)
}

cp02_fv_weighted96_fixture <- function(weighted = TRUE) {
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
  idx <- cp02_fv_fn("construct_observed_pairs")(pan, "id", "time", "y", "g")
  ctx <- cp02_fv_fn("construct_ipw_context")(idx, id = seq_len(n), Z = Z, weights = w)
  list(n = n, G = G, Z = Z, w = w, obs = obs, Y = Y, panel = pan, index = idx, context = ctx)
}

cp02_fv_newton_ref <- function() {
  ref_path <- cp02_fv_locate_file("tests", "testthat", "reference", "cp02_shared_ipw_audit_newton_reference.R")
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

test_that("owner-backed fitted nuisance view matches dense path bitwise and independent Newton reference across chunk sizes", {
  construct_ipw_panel_owner <- cp02_fv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_fv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_fv_fn("construct_ipw_nuisance_view")
  fit_ipw_nuisance_view <- cp02_fv_fn("fit_ipw_nuisance_view")
  fit_ipw_nuisance_from_owner <- cp02_fv_fn("fit_ipw_nuisance_from_owner")
  validate_ipw_fitted_nuisance_view <- cp02_fv_fn("validate_ipw_fitted_nuisance_view")
  ipw_owner_control <- cp02_fv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_fv_fn("ipw_fitted_view_control")
  fit_ipw_link <- cp02_fv_fn("fit_ipw_link")
  construct_ipw_shared_session <- cp02_fv_fn("construct_ipw_shared_session")
  construct_ipw_nuisance <- cp02_fv_fn("construct_ipw_nuisance")

  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  # 1. Deterministic 24-unit panel (unweighted and weighted across all 6 supported pairs)
  pairs24 <- list(c(1, 2), c(1, 4), c(2, 3), c(2, 4), c(3, 4), c(4, 5))
  for (weighted in c(FALSE, TRUE)) {
    fx24 <- cp02_fv_panel24_fixture(weighted = weighted)
    own_ctrl24 <- ipw_owner_control(
      max_chunk_rows = 7,
      max_retained_entries = 20000,
      max_matrix_entries = 100000,
      max_work_ops = 1000000,
      max_output_entries = 20000
    )
    po24 <- construct_ipw_panel_owner(fx24$panel, "id", "time", "y", "g", control = own_ctrl24)
    fo24 <- construct_ipw_frame_owner(po24, id = fx24$units$id, Z = fx24$Z, weights = fx24$w, control = own_ctrl24)

    for (pr in pairs24) {
      dfit <- fit_ipw_link(fx24$context, g = 3, K = Inf, s = pr[1], t = pr[2], control = glm_ctrl)
      for (ch_sz in c(1L, 5L, 24L)) {
        fv_ctrl <- ipw_fitted_view_control(
          max_chunk_rows = ch_sz,
          max_observer_units = 24,
          max_retained_entries = 20000,
          max_matrix_entries = 100000,
          max_work_ops = 1000000,
          max_output_entries = 20000
        )
        vfit <- fit_ipw_nuisance_from_owner(
          fo24, g = 3, K = Inf, s = pr[1], t = pr[2],
          glm_control = glm_ctrl, control = fv_ctrl
        )
        validate_ipw_fitted_nuisance_view(vfit, frame_owner = fo24)
        expect_identical(vfit$coefficients, dfit$coefficients)
        expect_identical(vfit$beta, dfit$coefficients)
        expect_identical(vfit$factors, dfit$factors)
        expect_identical(vfit$partials$A_beta_r, dfit$partials$A_beta_r)
        expect_identical(vfit$masses$rates, dfit$masses$rates)
        expect_identical(vfit$masses$q, dfit$masses$q)
        expect_identical(vfit$masses$cohort_a_weight, dfit$masses$cohort_a_weight)
        expect_identical(vfit$masses$cohort_pair_weight, dfit$masses$cohort_pair_weight)
        expect_identical(vfit$means$sum_a, dfit$means$sum_a)
        expect_identical(vfit$means$mean_a, dfit$means$mean_a)
        expect_identical(vfit$means$sum_b, dfit$means$sum_b)
        expect_identical(vfit$means$mean_b, dfit$means$mean_b)
        # The five compact D27 diagnostics are additive. Every original
        # shared diagnostic remains bit-identical to the dense fit.
        added_diagnostics <- c("weight_concentration", "weight_concentration_definition",
          "weight_concentration_sources", "fitted_probability_range", "fitted_probability_scope")
        expect_identical(names(vfit$diagnostics), c(added_diagnostics, names(dfit$diagnostics)))
        expect_identical(vfit$diagnostics[names(dfit$diagnostics)], dfit$diagnostics)
        expect_identical(vfit$diagnostics$weight_concentration_definition,
          "Kish(v) = (sum(v/max(v)))^2 / sum((v/max(v))^2); descriptive concentration only")
        expect_identical(vfit$diagnostics$weight_concentration_sources,
          c(propensity_fit_weights = "st$fit_weights", target_outcome_weights = "st$a_obs",
            pooled_control_outcome_ipw_weights = "st$b_obs"))
        expect_identical(vfit$diagnostics$fitted_probability_scope, "observed fitted risk set")

        # Reconstruct the risk set and inverse observation-rate weights from
        # the fixed panel, without a cdid support/weight/diagnostic helper.
        left_ids <- fx24$panel$id[fx24$panel$time == pr[1] & !is.na(fx24$panel$y)]
        right_ids <- fx24$panel$id[fx24$panel$time == pr[2] & !is.na(fx24$panel$y)]
        observed <- fx24$units$id %in% intersect(left_ids, right_ids)
        cohort <- match(fx24$units$g, c(3, Inf))
        rates_ref <- vapply(1:2, function(j)
          sum(fx24$w[observed & cohort == j]) / sum(fx24$w[cohort == j]), numeric(1))
        fit_w_ref <- fx24$w[observed] / rates_ref[cohort[observed]]
        eta_ref <- drop(fx24$Z[observed, , drop = FALSE] %*% dfit$coefficients)
        probability_ref <- stats::plogis(eta_ref)
        treated <- cohort[observed] == 1L
        weights_ref <- list(propensity_fit_weights = fit_w_ref,
          target_outcome_weights = fit_w_ref[treated],
          pooled_control_outcome_ipw_weights = fit_w_ref[!treated] * exp(eta_ref[!treated]))
        expect_identical(names(vfit$diagnostics$weight_concentration), names(weights_ref))
        for (nm in names(weights_ref)) {
          summary <- vfit$diagnostics$weight_concentration[[nm]]
          expect_identical(names(summary), c("status", "reason", "value"))
          expect_identical(summary$status, "available")
          expect_identical(summary$reason, NA_character_)
          w_ref <- weights_ref[[nm]]
          expect_equal(summary$value, sum(w_ref)^2 / sum(w_ref^2), tolerance = 1e-12)
        }
        expect_equal(vfit$diagnostics$fitted_probability_range,
          list(status = "available", reason = NA_character_,
            min = min(probability_ref), max = max(probability_ref)), tolerance = 1e-14)

        # Independently form H/n_A and its exact two-by-two inverse for the
        # declared 1-norm conditioning diagnostic; this is not an rcond oracle.
        Z_ref <- fx24$Z[observed, , drop = FALSE]
        H_ref <- crossprod(Z_ref, Z_ref * as.numeric(fit_w_ref * probability_ref * (1 - probability_ref))) /
          nrow(fx24$Z)
        H_inverse_ref <- matrix(c(H_ref[2, 2], -H_ref[2, 1], -H_ref[1, 2], H_ref[1, 1]), 2) /
          (H_ref[1, 1] * H_ref[2, 2] - H_ref[1, 2] * H_ref[2, 1])
        rcond_ref <- 1 / (max(colSums(abs(H_ref))) * max(colSums(abs(H_inverse_ref))))
        conditioning <- cp02_fv_fn("chained_request_H_diagnostics")(vfit$factors$H, nrow(fx24$Z))
        expect_equal(unname(vfit$factors$H), unname(H_ref), tolerance = 1e-12)
        expect_equal(conditioning$rcond_H_1_estimate$value, rcond_ref, tolerance = 1e-12)
        expect_identical(conditioning$norm_label, "1-norm")
        expect_identical(conditioning$n_A, nrow(fx24$Z))
      }
    }
  }

  # 2. 96-unit audit_fixture (weighted and unweighted) + independent Newton reference
  independent_fit <- cp02_fv_newton_ref()
  pairs96 <- list(c(1, 2), c(2, 3), c(1, 3), c(7, 8))
  for (weighted in c(TRUE, FALSE)) {
    fx96 <- cp02_fv_weighted96_fixture(weighted = weighted)
    own_ctrl96 <- ipw_owner_control(
      max_chunk_rows = 19,
      max_retained_entries = 100000,
      max_matrix_entries = 500000,
      max_work_ops = 5000000,
      max_output_entries = 100000
    )
    po96 <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl96)
    fo96 <- construct_ipw_frame_owner(po96, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl96)
    ses96 <- construct_ipw_shared_session(fx96$context, 4, 20, 1e5)

    for (pr in pairs96) {
      dfit <- fit_ipw_link(fx96$context, g = 4, K = c(7, Inf), s = pr[1], t = pr[2], control = glm_ctrl, functional_key = "A")
      dnu <- construct_ipw_nuisance(ses96, 4, c(7, Inf), pr[1], pr[2], glm_ctrl)
      ind <- independent_fit(fx96, pr[1], pr[2])

      for (ch_sz in c(1L, 5L, 13L, 96L)) {
        fv_ctrl <- ipw_fitted_view_control(
          max_chunk_rows = ch_sz,
          max_observer_units = 96,
          max_retained_entries = 100000,
          max_matrix_entries = 500000,
          max_work_ops = 5000000,
          max_output_entries = 100000
        )
        nv <- construct_ipw_nuisance_view(
          fo96, g = 4, K = c(7, Inf), s = pr[1], t = pr[2],
          control = ipw_owner_control(ch_sz, 100000, 500000, 5000000, 100000)
        )
        vfit <- fit_ipw_nuisance_view(nv, glm_control = glm_ctrl, control = fv_ctrl)
        validate_ipw_fitted_nuisance_view(vfit, nuisance_view = nv, frame_owner = fo96)

        # Bitwise equality with dense fit_ipw_link and construct_ipw_nuisance
        expect_identical(vfit$coefficients, dfit$coefficients)
        expect_identical(vfit$factors, dfit$factors)
        expect_identical(vfit$partials$A_beta_r, dfit$partials$A_beta_r)
        expect_identical(vfit$partials$A_beta_r, dnu$hessian$A_beta_r)
        expect_identical(vfit$factors$H, dnu$hessian$H)
        expect_identical(vfit$factors$H_chol, dnu$hessian$H_chol)
        expect_identical(vfit$masses$rates, dfit$masses$rates)
        expect_identical(vfit$means$sum_a, dnu$stage$sum_a)
        expect_identical(vfit$means$mean_a, dnu$stage$mean_a)
        expect_identical(vfit$means$sum_b, dnu$stage$sum_b)
        expect_identical(vfit$means$mean_b, dnu$stage$mean_b)

        # Agreement with independent Newton reference
        expect_lt(max(abs(vfit$coefficients - ind$beta)), 1e-10)
        expect_lt(max(abs(vfit$factors$H - ind$H)), 1e-12)
        expect_lt(max(abs(vfit$partials$A_beta_r - ind$L)), 1e-11)
        expect_lt(max(abs(vfit$masses$rates - ind$r)), 1e-14)

        # Point-summary disclosure records both observer full sums and chunked view sums
        psd <- vfit$point_summary_disclosure
        expect_true(isTRUE(psd$candidate_point_summary_only))
        expect_false(psd$certifies_root_or_inference)
        expect_false(psd$glm_convergence_is_certificate)
        expect_false(psd$seal_is_certificate)
        expect_identical(psd$observer_rates, vfit$masses$rates)
        expect_identical(psd$view_chunked_rates, nv$support$rates)
        expect_lt(psd$max_abs_rate_chunk_diff, 1e-14)
      }
    }
  }
})

test_that("outcome-free reads, observer ordering, chunk release, and compact retention are enforced", {
  construct_ipw_panel_owner <- cp02_fv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_fv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_fv_fn("construct_ipw_nuisance_view")
  fit_ipw_nuisance_view <- cp02_fv_fn("fit_ipw_nuisance_view")
  fit_ipw_nuisance_view_set <- cp02_fv_fn("fit_ipw_nuisance_view_set")
  ipw_owner_control <- cp02_fv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_fv_fn("ipw_fitted_view_control")

  fx96 <- cp02_fv_weighted96_fixture(weighted = TRUE)
  # Add singleton and unobserved units to verify common owner retains them while fitted view stays compact
  extra_rows <- rbind(
    data.frame(id = 97L, time = 1, y = 1.25, g = 4),
    data.frame(id = 98L, time = 1, y = NA_real_, g = 7),
    data.frame(id = 98L, time = 2, y = NA_real_, g = 7)
  )
  pan_ext <- rbind(fx96$panel, extra_rows)

  own_ctrl <- ipw_owner_control(11, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(11, 96, 100000, 500000, 5000000, 100000)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  po <- construct_ipw_panel_owner(pan_ext, "id", "time", "y", "g", control = own_ctrl)
  expect_identical(po$n_singleton_units, 1L)
  expect_identical(po$n_unobserved_units, 1L)
  expect_identical(po$n_a_units, 96L)

  # Supply frame in reversed unit order; canonical A-unit order must yield identical fit
  rev_id <- rev(seq_len(96L))
  fo_rev <- construct_ipw_frame_owner(
    po, id = rev_id, Z = fx96$Z[rev_id, , drop = FALSE],
    weights = fx96$w[rev_id], control = own_ctrl
  )
  fo_fwd <- construct_ipw_frame_owner(
    po, id = seq_len(96L), Z = fx96$Z,
    weights = fx96$w, control = own_ctrl
  )

  nv12_rev <- construct_ipw_nuisance_view(fo_rev, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  nv12_fwd <- construct_ipw_nuisance_view(fo_fwd, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  nv78_fwd <- construct_ipw_nuisance_view(fo_fwd, g = 4, K = c(7, Inf), s = 7, t = 8, control = own_ctrl)

  cb_events <- list()
  tr <- new.env(parent = emptyenv())
  vfit_rev <- fit_ipw_nuisance_view(
    nv12_rev, glm_control = glm_ctrl, control = fv_ctrl,
    callback = function(ev) { cb_events[[length(cb_events) + 1L]] <<- ev },
    trace_env = tr
  )
  vfit_fwd <- fit_ipw_nuisance_view(nv12_fwd, glm_control = glm_ctrl, control = fv_ctrl)

  expect_identical(vfit_rev$coefficients, vfit_fwd$coefficients)
  expect_identical(vfit_rev$factors, vfit_fwd$factors)
  expect_identical(vfit_rev$partials, vfit_fwd$partials)
  expect_identical(vfit_rev$masses, vfit_fwd$masses)
  expect_identical(vfit_rev$means, vfit_fwd$means)

  # Outcome-free reads and chunk release checks
  expect_identical(tr$endpoints_read, 0L)
  expect_identical(tr$max_simultaneous_query_chunks, 1L)
  expect_identical(tr$max_simultaneous_observer_buffers, 1L)
  expect_identical(tr$active_observer_buffers, 0L)
  expect_true(isTRUE(tr$observer_released))
  expect_true(length(cb_events) >= 9L)
  chunk_events <- Filter(function(e) identical(e$phase, "after_chunk_release"), cb_events)
  expect_true(all(vapply(chunk_events, function(e) identical(e$active_query_chunks, 0L), logical(1L))))
  expect_true(all(vapply(chunk_events, function(e) identical(e$endpoints_read, 0L), logical(1L))))

  # Multi-view set releases all observer/model arrays before the next distinct nuisance fit
  set_tr <- new.env(parent = emptyenv())
  fit_set <- fit_ipw_nuisance_view_set(
    list(cell12 = nv12_fwd, cell78 = nv78_fwd),
    glm_control = glm_ctrl,
    control = fv_ctrl,
    trace_env = set_tr
  )
  expect_identical(set_tr$n_distinct_fits, 2L)
  expect_identical(set_tr$between_fit_active_buffers, c(0L, 0L))
  expect_true(isTRUE(set_tr$all_released_between_fits))
  expect_identical(fit_set$cell12$counts$n_observed_risk_units, 54L)
  expect_identical(fit_set$cell78$counts$n_observed_risk_units, 18L)

  # Verify no n_obs (54) or n_a (96) vectors/matrices exist anywhere in fitted_view
  check_no_unit_arrays <- function(x) {
    if (is.environment(x)) return(TRUE)
    if (is.atomic(x)) {
      d <- dim(x)
      if (!is.null(d)) {
        expect_false(any(d %in% c(54L, 72L, 96L)))
      } else if (!is.character(x)) {
        expect_false(length(x) %in% c(54L, 72L, 96L))
      }
      return(TRUE)
    }
    if (is.list(x)) {
      for (el in x) check_no_unit_arrays(el)
    }
    TRUE
  }
  check_no_unit_arrays(fit_set$cell12)
})

test_that("explicit insufficient resource caps refuse before any query, callback, or solver call", {
  construct_ipw_panel_owner <- cp02_fv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_fv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_fv_fn("construct_ipw_nuisance_view")
  fit_ipw_nuisance_view <- cp02_fv_fn("fit_ipw_nuisance_view")
  ipw_owner_control <- cp02_fv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_fv_fn("ipw_fitted_view_control")

  fx96 <- cp02_fv_weighted96_fixture(weighted = TRUE)
  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  po <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)
  nv <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  # Missing explicit control fields refused immediately
  expect_error(
    ipw_fitted_view_control(max_chunk_rows = 16, max_observer_units = 96),
    class = "cdid_ipw_fitted_view_control_error"
  )

  bad_controls <- list(
    obs_cap = ipw_fitted_view_control(16, max_observer_units = 53, 100000, 500000, 5000000, 100000),
    ret_cap = ipw_fitted_view_control(16, 96, max_retained_entries = 10, 500000, 5000000, 100000),
    mat_cap = ipw_fitted_view_control(16, 96, 100000, max_matrix_entries = 100, 5000000, 100000),
    ops_cap = ipw_fitted_view_control(16, 96, 100000, 500000, max_work_ops = 500, 100000),
    out_cap = ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, max_output_entries = 10)
  )

  for (bc in bad_controls) {
    cb_count <- 0L
    tr <- new.env(parent = emptyenv())
    expect_error(
      fit_ipw_nuisance_view(
        nv, glm_control = glm_ctrl, control = bc,
        callback = function(...) { cb_count <<- cb_count + 1L },
        trace_env = tr
      ),
      class = "cdid_ipw_owner_bound_error"
    )
    expect_identical(cb_count, 0L)
    expect_null(tr$queries_started)
  }
})

test_that("altered owners, controls, ordered cohorts, and mismatched masks are refused on validation and sharing", {
  construct_ipw_panel_owner <- cp02_fv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_fv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_fv_fn("construct_ipw_nuisance_view")
  construct_ipw_link_view <- cp02_fv_fn("construct_ipw_link_view")
  fit_ipw_nuisance_view <- cp02_fv_fn("fit_ipw_nuisance_view")
  validate_ipw_fitted_nuisance_view <- cp02_fv_fn("validate_ipw_fitted_nuisance_view")
  verify_ipw_fitted_nuisance_sharing <- cp02_fv_fn("verify_ipw_fitted_nuisance_sharing")
  ipw_owner_control <- cp02_fv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_fv_fn("ipw_fitted_view_control")

  fx96 <- cp02_fv_weighted96_fixture(weighted = TRUE)
  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  po <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)
  nv12 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  vfit12 <- fit_ipw_nuisance_view(nv12, glm_control = glm_ctrl, control = fv_ctrl)

  # Compatible link view on (2, 3) shares exact mask with (1, 2)
  lv23 <- construct_ipw_link_view(fo, s = 2, t = 3, nuisance_view = nv12, control = own_ctrl)
  sh_ok <- verify_ipw_fitted_nuisance_sharing(vfit12, lv23, glm_control = glm_ctrl, control = fv_ctrl)
  expect_true(isTRUE(sh_ok$compatible))
  expect_identical(sh_ok$endpoints_read, 0L)

  # Tampered fitted view breaks seal
  tampered <- vfit12
  tampered$coefficients[1L] <- tampered$coefficients[1L] + 1e-12
  expect_error(
    validate_ipw_fitted_nuisance_view(tampered),
    class = "cdid_ipw_fitted_view_integrity_error"
  )

  # Distinct frame owner with altered weights or covariates refused
  fo_alt_w <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = rep(1, fx96$n), control = own_ctrl)
  expect_error(
    validate_ipw_fitted_nuisance_view(vfit12, frame_owner = fo_alt_w),
    class = "cdid_ipw_owner_mismatch_error"
  )
  lv_alt_w <- construct_ipw_link_view(fo_alt_w, s = 1, t = 2, g = 4, K = c(7, Inf), control = own_ctrl)
  expect_error(
    verify_ipw_fitted_nuisance_sharing(vfit12, lv_alt_w, glm_control = glm_ctrl, control = fv_ctrl),
    class = "cdid_ipw_owner_mismatch_error"
  )

  # Mismatched glm.control refused
  glm_alt <- stats::glm.control(epsilon = 1e-8, maxit = 100, trace = FALSE)
  expect_error(
    verify_ipw_fitted_nuisance_sharing(vfit12, lv23, glm_control = glm_alt, control = fv_ctrl),
    class = "cdid_ipw_control_mismatch_error"
  )

  # Reordered control cohorts K = c(Inf, 7) vs c(7, Inf) refused
  lv_reordered_K <- construct_ipw_link_view(fo, s = 1, t = 2, g = 4, K = c(Inf, 7), control = own_ctrl)
  expect_error(
    verify_ipw_fitted_nuisance_sharing(vfit12, lv_reordered_K, glm_control = glm_ctrl, control = fv_ctrl),
    class = "cdid_ipw_owner_cohort_error"
  )

  # Mismatched risk-set observer mask (7, 8) vs (1, 2) refused
  lv78 <- construct_ipw_link_view(fo, s = 7, t = 8, g = 4, K = c(7, Inf), control = own_ctrl)
  expect_error(
    verify_ipw_fitted_nuisance_sharing(vfit12, lv78, glm_control = glm_ctrl, control = fv_ctrl),
    class = "cdid_ipw_owner_mask_mismatch_error"
  )
})

test_that("genuine rank deficiency, separation boundary, and non-convergence refuse without retry and keep failed predictions in error diagnostics only", {
  construct_ipw_panel_owner <- cp02_fv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_fv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_fv_fn("construct_ipw_nuisance_view")
  fit_ipw_nuisance_view <- cp02_fv_fn("fit_ipw_nuisance_view")
  ipw_owner_control <- cp02_fv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_fv_fn("ipw_fitted_view_control")

  fx96 <- cp02_fv_weighted96_fixture(weighted = TRUE)
  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  po <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)

  # 1. Genuine rank deficiency on cross-boundary pair (1, 7) of audit_fixture(TRUE)
  nv17 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 7, control = own_ctrl)
  fail_rank <- new.env(parent = emptyenv())
  tr_rank <- new.env(parent = emptyenv())
  cnd_rank <- tryCatch(
    fit_ipw_nuisance_view(nv17, glm_control = glm_ctrl, control = fv_ctrl, trace_env = tr_rank, failure_trace_env = fail_rank),
    error = function(e) e
  )
  expect_s3_class(cnd_rank, "cdid_ipw_rank_error")
  expect_true(isTRUE(fail_rank$observer_released_on_error))
  expect_identical(tr_rank$active_observer_buffers, 0L)
  expect_true(is.list(cnd_rank$error_diagnostics))
  expect_identical(cnd_rank$error_diagnostics$stage, "rank_check")
  expect_true(is.numeric(cnd_rank$error_diagnostics$failed_fitted_values))

  # 2. Genuine separation / boundary probability refusal
  Z_sep <- cbind(one = 1, sep = ifelse(fx96$G == 4, 40, -40))
  fo_sep <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = Z_sep, weights = fx96$w, control = own_ctrl)
  nv_sep <- construct_ipw_nuisance_view(fo_sep, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  fail_sep <- new.env(parent = emptyenv())
  tr_sep <- new.env(parent = emptyenv())
  cnd_sep <- tryCatch(
    fit_ipw_nuisance_view(nv_sep, glm_control = glm_ctrl, control = fv_ctrl, trace_env = tr_sep, failure_trace_env = fail_sep),
    error = function(e) e
  )
  expect_s3_class(cnd_sep, "cdid_ipw_boundary_error")
  expect_true(isTRUE(fail_sep$observer_released_on_error))
  expect_identical(tr_sep$active_observer_buffers, 0L)
  expect_true(is.list(cnd_sep$error_diagnostics))
  expect_length(cnd_sep$error_diagnostics$failed_fitted_values, 54L)

  # 3. Genuine non-convergence under tight maxit = 1 (no retry with larger maxit)
  nv12 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  glm_1iter <- stats::glm.control(epsilon = 1e-12, maxit = 1, trace = FALSE)
  fail_conv <- new.env(parent = emptyenv())
  tr_conv <- new.env(parent = emptyenv())
  cnd_conv <- tryCatch(
    fit_ipw_nuisance_view(nv12, glm_control = glm_1iter, control = fv_ctrl, trace_env = tr_conv, failure_trace_env = fail_conv),
    error = function(e) e
  )
  expect_s3_class(cnd_conv, "cdid_ipw_convergence_error")
  expect_true(isTRUE(fail_conv$observer_released_on_error))
  expect_identical(tr_conv$active_observer_buffers, 0L)
  expect_identical(cnd_conv$error_diagnostics$stage, "convergence_check")
  expect_identical(cnd_conv$error_diagnostics$iter, 1L)
  expect_length(cnd_conv$error_diagnostics$failed_coefficients, 2L)
  expect_length(cnd_conv$error_diagnostics$failed_fitted_values, 54L)
})

test_that("exact binary64 RDS snapshots require and pass seal/ownership revalidation and refuse tampered values", {
  construct_ipw_panel_owner <- cp02_fv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_fv_fn("construct_ipw_frame_owner")
  construct_ipw_nuisance_view <- cp02_fv_fn("construct_ipw_nuisance_view")
  fit_ipw_nuisance_view <- cp02_fv_fn("fit_ipw_nuisance_view")
  validate_ipw_fitted_nuisance_view <- cp02_fv_fn("validate_ipw_fitted_nuisance_view")
  extract_ipw_fitted_nuisance_snapshot <- cp02_fv_fn("extract_ipw_fitted_nuisance_snapshot")
  revalidate_ipw_fitted_nuisance_view <- cp02_fv_fn("revalidate_ipw_fitted_nuisance_view")
  ipw_owner_control <- cp02_fv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_fv_fn("ipw_fitted_view_control")

  fx96 <- cp02_fv_weighted96_fixture(weighted = TRUE)
  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  po <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)
  nv12 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 1, t = 2, control = own_ctrl)
  vfit12 <- fit_ipw_nuisance_view(nv12, glm_control = glm_ctrl, control = fv_ctrl)
  snap12 <- extract_ipw_fitted_nuisance_snapshot(vfit12)

  rds_view_path <- tempfile(fileext = ".rds")
  rds_snap_path <- tempfile(fileext = ".rds")
  on.exit(unlink(c(rds_view_path, rds_snap_path)), add = TRUE)

  saveRDS(vfit12, rds_view_path)
  saveRDS(snap12, rds_snap_path)

  loaded_view <- readRDS(rds_view_path)
  loaded_snap <- readRDS(rds_snap_path)

  # Exact binary64 arrays preserved across RDS round-trip
  expect_identical(loaded_view$coefficients, vfit12$coefficients)
  expect_identical(loaded_view$factors, vfit12$factors)
  expect_identical(loaded_view$partials, vfit12$partials)
  expect_identical(loaded_snap$coefficients, vfit12$coefficients)

  # Raw deserialized objects must be refused by validate_ipw_fitted_nuisance_view prior to revalidation
  expect_error(
    validate_ipw_fitted_nuisance_view(loaded_view),
    class = "cdid_ipw_fitted_view_integrity_error"
  )
  expect_error(
    validate_ipw_fitted_nuisance_view(loaded_snap),
    class = "cdid_ipw_fitted_view_integrity_error"
  )

  # Revalidation against the live sealed nuisance_view succeeds and reseals
  rev_from_view <- revalidate_ipw_fitted_nuisance_view(loaded_view, nv12, glm_control = glm_ctrl, control = fv_ctrl)
  rev_from_snap <- revalidate_ipw_fitted_nuisance_view(loaded_snap, nv12, glm_control = glm_ctrl, control = fv_ctrl)
  validate_ipw_fitted_nuisance_view(rev_from_view, nuisance_view = nv12, frame_owner = fo)
  validate_ipw_fitted_nuisance_view(rev_from_snap, nuisance_view = nv12, frame_owner = fo)
  expect_identical(rev_from_view$coefficients, vfit12$coefficients)
  expect_identical(rev_from_snap$factors, vfit12$factors)

  # Perturbing any stored binary64 array or revalidating against a different nuisance cell is refused
  bad_snap <- loaded_snap
  bad_snap$partials$A_beta_r[1L, 1L] <- bad_snap$partials$A_beta_r[1L, 1L] + 1e-12
  expect_error(
    revalidate_ipw_fitted_nuisance_view(bad_snap, nv12, glm_control = glm_ctrl, control = fv_ctrl),
    class = "cdid_ipw_fitted_view_revalidation_error"
  )
  nv78 <- construct_ipw_nuisance_view(fo, g = 4, K = c(7, Inf), s = 7, t = 8, control = own_ctrl)
  expect_error(
    revalidate_ipw_fitted_nuisance_view(loaded_snap, nv78, glm_control = glm_ctrl, control = fv_ctrl),
    class = "cdid_ipw_fitted_view_revalidation_error"
  )
})

test_that("discriminating formula mutants are strongly separated from fitted nuisance outputs", {
  construct_ipw_panel_owner <- cp02_fv_fn("construct_ipw_panel_owner")
  construct_ipw_frame_owner <- cp02_fv_fn("construct_ipw_frame_owner")
  fit_ipw_nuisance_from_owner <- cp02_fv_fn("fit_ipw_nuisance_from_owner")
  ipw_owner_control <- cp02_fv_fn("ipw_owner_control")
  ipw_fitted_view_control <- cp02_fv_fn("ipw_fitted_view_control")

  fx96 <- cp02_fv_weighted96_fixture(weighted = TRUE)
  own_ctrl <- ipw_owner_control(16, 100000, 500000, 5000000, 100000)
  fv_ctrl <- ipw_fitted_view_control(16, 96, 100000, 500000, 5000000, 100000)
  glm_ctrl <- stats::glm.control(epsilon = 1e-12, maxit = 100, trace = FALSE)

  po <- construct_ipw_panel_owner(fx96$panel, "id", "time", "y", "g", control = own_ctrl)
  fo <- construct_ipw_frame_owner(po, id = seq_len(fx96$n), Z = fx96$Z, weights = fx96$w, control = own_ctrl)
  vfit <- fit_ipw_nuisance_from_owner(fo, g = 4, K = c(7, Inf), s = 1, t = 2, glm_control = glm_ctrl, control = fv_ctrl)

  # Mutant 1: Omitting inverse-rate propensity weight 1 / r_G (using raw w_i in glm.fit)
  in_o <- which(fx96$G %in% c(4, 7, Inf) & fx96$obs[, 1] & fx96$obs[, 2])
  D_o <- as.double(fx96$G[in_o] == 4)
  Z_o <- fx96$Z[in_o, , drop = FALSE]
  w_o <- fx96$w[in_o]
  mut1_beta <- as.double(suppressWarnings(stats::glm.fit(Z_o, D_o, weights = w_o, family = stats::binomial(), control = glm_ctrl)$coefficients))
  expect_gt(max(abs(vfit$coefficients - mut1_beta)), 1e-2)

  # Mutant 2: Omitting 1 / n scaling in Hessian H
  mut2_H <- vfit$factors$H * fx96$n
  expect_gt(max(abs(vfit$factors$H - mut2_H)), 10)

  # Mutant 3: Using -1 / (n * r_j) instead of -1 / (n * r_j^2) in score-rate Jacobian A_beta_r
  mut3_A <- sweep(vfit$partials$A_beta_r, 2L, vfit$masses$rates, `*`)
  expect_gt(max(abs(vfit$partials$A_beta_r - mut3_A)), 5e-2)

  # Mutant 4: Omitting odds exp(eta) in control normalizer sum_b
  gi_o <- match(fx96$G[in_o], c(4, 7, Inf))
  fw_o <- w_o / vfit$masses$rates[gi_o]
  mut4_sum_b <- sum(fw_o[D_o == 0])
  expect_gt(abs(vfit$means$sum_b - mut4_sum_b), 10)
})

# Local instrumentation clones exact installed functions without altering namespace bindings.
fv_local_api <- function() {
 e<-new.env(parent=api)
 for(nm in c("ipw_owner_control","construct_ipw_panel_owner","construct_ipw_frame_owner","construct_ipw_nuisance_view","query_ipw_view_chunk")) assign(nm,get(nm,envir=api,inherits=FALSE),envir=e)
 for(nm in c("cdid_ipw_fitted_view_abort","ipw_fitted_view_require_deps","ipw_fitted_view_control","validate_ipw_fitted_view_control","ipw_fitted_view_to_owner_control","ipw_fitted_view_char_bytes","preflight_ipw_fitted_nuisance_view","ipw_fitted_view_verify_compact_storage","fit_ipw_nuisance_view","construct_ipw_fitted_nuisance_view","fit_ipw_nuisance_from_owner","fit_ipw_nuisance_view_set","validate_ipw_fitted_nuisance_view","verify_ipw_fitted_nuisance_sharing","extract_ipw_fitted_nuisance_snapshot","ipw_fitted_snapshot_payload","revalidate_ipw_fitted_nuisance_view","ipw_nuisance_observer_fit_worker")) {f<-get(nm,envir=api,inherits=FALSE);environment(f)<-e;assign(nm,f,envir=e)}
 e
}

test_that("actual local slice lifetimes, early preflight and snapshot schema are bounded", {
 local<-fv_local_api();fx<-cp02_fv_weighted96_fixture(TRUE)
 oc<-local$ipw_owner_control(11,1e5,5e5,5e6,1e5)
 fc<-local$ipw_fitted_view_control(11,96,1e5,5e5,5e6,1e5)
 gc<-stats::glm.control(epsilon=1e-12,maxit=100,trace=FALSE)
 p<-local$construct_ipw_panel_owner(fx$panel,"id","time","y","g",oc)
 f<-local$construct_ipw_frame_owner(p,1:96,fx$Z,weights=fx$w,control=oc)
 v<-local$construct_ipw_nuisance_view(f,4,c(7,Inf),1,2,control=oc)
 q<-local$query_ipw_view_chunk;aliases<-integer()
 local$query_ipw_view_chunk<-function(...) {
  aliases<<-c(aliases,length(intersect(c("a_ch","gi_ch","loc_ch","w_ch","Z_ch","dst"),ls(parent.frame()))));q(...)
 }
 callbacks<-integer()
 fit<-local$fit_ipw_nuisance_view(v,gc,fc,callback=function(e) {
  callbacks<<-c(callbacks,length(intersect(c("a_ch","gi_ch","loc_ch","w_ch","Z_ch","dst"),ls(parent.frame()))))
 })
 expect_true(length(aliases)>1L);expect_true(all(aliases==0L));expect_true(all(callbacks==0L))
 expect_gte(fit$preflight$referenced_nuisance_entries,v$preflight$total_retained_entries)
 ch<-local$ipw_fitted_view_char_bytes;char_calls<-0L
 local$ipw_fitted_view_char_bytes<-function(...) {char_calls<<-char_calls+1L;ch(...)}
 low<-fc;low$max_retained_entries<-10
 expect_error(local$fit_ipw_nuisance_view(v,gc,low),class="cdid_ipw_owner_bound_error")
 expect_identical(char_calls,0L)
 low_work<-fc;low_work$max_work_ops<-1
 expect_error(local$fit_ipw_nuisance_view(v,gc,low_work),class="cdid_ipw_owner_bound_error")
 expect_identical(char_calls,0L)
 worker<-local$ipw_nuisance_observer_fit_worker;calls<-0L
 local$ipw_nuisance_observer_fit_worker<-function(...) {calls<<-calls+1L;worker(...)}
 expect_error(local$fit_ipw_nuisance_view_set(list(v,v),gc,fc,trace_env=42),class="cdid_ipw_fitted_view_input_error")
 expect_identical(calls,0L)
 snap<-local$extract_ipw_fitted_nuisance_snapshot(fit);snap$foreign_extra<-1:10
 expect_error(local$revalidate_ipw_fitted_nuisance_view(snap,v,gc,fc),class="cdid_ipw_fitted_view_revalidation_error")
 expect_identical(calls,0L)
 snap<-local$extract_ipw_fitted_nuisance_snapshot(fit);snap$coefficients<-seq_len(5e5)
 expect_error(local$revalidate_ipw_fitted_nuisance_view(snap,v,gc,fc),class="cdid_ipw_owner_bound_error")
 expect_identical(calls,0L)
 snap<-local$extract_ipw_fitted_nuisance_snapshot(fit);snap$control$max_chunk_rows<-10
 expect_error(local$revalidate_ipw_fitted_nuisance_view(snap,v,gc,fc),class="cdid_ipw_control_mismatch_error")
 expect_identical(calls,0L)
 # Independent composition rule: all output envelopes coexist with one largest workspace.
 pre<-local$preflight_ipw_fitted_nuisance_view(v,gc,fc)
 need<-2*pre$total_retained_entries+(pre$peak_concurrent_entries-pre$total_retained_entries)+16*2+128
 tight<-fc;tight$max_matrix_entries<-need-1
 expect_error(local$fit_ipw_nuisance_view_set(list(v,v),gc,tight),class="cdid_ipw_owner_bound_error")
 expect_identical(calls,0L)
})
