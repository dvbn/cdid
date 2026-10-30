# Focused compatibility tests for `validate_ipw_glm_control` and IPW link entry points
# across canonical 3-field (`epsilon`, `maxit`, `trace`) and R >= 4.7 native-default
# 5-field (`epsilon`, `maxit`, `trace`, `tol = min(1e-7, epsilon / 1000)`, `wtol = 0`)
# `stats::glm.control` lists.

if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}

build_cp02_control_fixture_panel <- function(n, Y, observed, G, ids = seq_len(n), times = c(1, 2, 3)) {
  rows <- vector("list", n * length(times))
  k <- 0L
  for (i in seq_len(n)) {
    for (t_idx in seq_along(times)) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        id = ids[i],
        time = as.double(times[t_idx]),
        y = if (observed[i, t_idx]) as.double(Y[i, t_idx]) else NA_real_,
        g = as.double(G[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("validate_ipw_glm_control normalizes canonical 3-field and native-default 5-field controls across default and custom epsilons", {
  cases <- list(
    default_eps = list(eps = 1e-8, maxit = 25L, trace = FALSE, tol = min(1e-7, 1e-8 / 1000)),
    custom_eps_1e12 = list(eps = 1e-12, maxit = 100L, trace = FALSE, tol = min(1e-7, 1e-12 / 1000)),
    custom_eps_1e15 = list(eps = 1e-15, maxit = 1L, trace = TRUE, tol = min(1e-7, 1e-15 / 1000)),
    capped_tol_eps = list(eps = 1e-2, maxit = 40L, trace = FALSE, tol = min(1e-7, 1e-2 / 1000)),
    integer_eps = list(eps = 1L, maxit = 10L, trace = TRUE, tol = min(1e-7, 1L / 1000)),
    cap_boundary_eps = list(eps = 1e-4, maxit = 25L, trace = FALSE, tol = min(1e-7, 1e-4 / 1000)),
    underflow_eps = list(eps = .Machine$double.xmin * .Machine$double.eps,
                         maxit = 25L, trace = FALSE,
                         tol = min(1e-7, (.Machine$double.xmin * .Machine$double.eps) / 1000)),
    largest_eps_maxit = list(eps = .Machine$double.xmax, maxit = .Machine$integer.max,
                            trace = TRUE, tol = min(1e-7, .Machine$double.xmax / 1000))
  )

  for (nm in names(cases)) {
    cs <- cases[[nm]]
    expected <- list(
      epsilon = as.double(cs$eps),
      maxit = as.integer(cs$maxit),
      trace = as.logical(cs$trace)
    )

    # Canonical 3-field shapes: standard order, reordered names, and numeric/integer representations.
    c3_std <- list(epsilon = cs$eps, maxit = cs$maxit, trace = cs$trace)
    c3_reord1 <- list(trace = cs$trace, maxit = as.double(cs$maxit), epsilon = cs$eps)
    c3_reord2 <- list(maxit = as.integer(cs$maxit), epsilon = as.double(cs$eps), trace = cs$trace)

    # Native-default 5-field shapes: standard order, reordered names, and integer/double 0 wtol.
    c5_std <- list(epsilon = cs$eps, maxit = cs$maxit, trace = cs$trace, tol = cs$tol, wtol = 0)
    c5_int_wtol <- list(epsilon = as.double(cs$eps), maxit = as.double(cs$maxit), trace = cs$trace, tol = cs$tol, wtol = 0L)
    c5_reord1 <- list(wtol = 0, tol = cs$tol, trace = cs$trace, maxit = cs$maxit, epsilon = cs$eps)
    c5_reord2 <- list(tol = cs$tol, epsilon = cs$eps, wtol = 0L, maxit = as.double(cs$maxit), trace = cs$trace)
    c_runtime <- stats::glm.control(epsilon = cs$eps, maxit = cs$maxit, trace = cs$trace)

    for (cand in list(c3_std, c3_reord1, c3_reord2, c5_std, c5_int_wtol, c5_reord1, c5_reord2, c_runtime)) {
      got <- cdid:::validate_ipw_glm_control(cand)
      expect_identical(got, expected, info = nm)
      expect_identical(names(got), c("epsilon", "maxit", "trace"), info = nm)
      expect_identical(length(got), 3L, info = nm)
    }
  }
})

test_that("validate_ipw_glm_control and fit_ipw_link strictly reject nondefault tol/wtol, partials, duplicates, unknowns, and malformed values", {
  set.seed(20261030)
  n <- 60L
  Z <- cbind(intercept = 1, x = rnorm(n))
  G <- rep(c(3, 4, Inf), each = 20L)
  obs <- matrix(TRUE, n, 3L)
  Y <- matrix(rnorm(n * 3L), n, 3L)
  panel <- build_cp02_control_fixture_panel(n, Y, obs, G, times = c(1, 2, 3))
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z)

  msg_missing <- "Supply an explicit `glm.control` list; no default solver policy is adopted."
  msg_invalid <- "Invalid `glm.control` specification; require exactly unambiguous `epsilon`, `maxit`, `trace` fields, finite positive `epsilon`, representable integer `maxit >= 1`, and logical `trace`."

  err_missing <- expect_error(cdid:::validate_ipw_glm_control(), class = "cdid_ipw_control_error")
  expect_true(inherits(err_missing, "cdid_ipw_input_error"))
  expect_identical(conditionMessage(err_missing), msg_missing)

  for (non_list in list(NULL, c(epsilon = 1e-8, maxit = 25, trace = 0), "glm.control", new.env(parent = emptyenv()))) {
    err_nl <- expect_error(cdid:::validate_ipw_glm_control(non_list), class = "cdid_ipw_control_error")
    expect_true(inherits(err_nl, "cdid_ipw_input_error"))
    expect_identical(conditionMessage(err_nl), msg_missing)
  }

  tol_def <- min(1e-7, 1e-8 / 1000)
  tol_c12 <- min(1e-7, 1e-12 / 1000)
  bad_controls <- list(
    # Nondefault tol values (catches mistakenly dropping/ignoring tol or wrong formula).
    tol_uncapped_at_default_eps = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = 1e-7, wtol = 0),
    tol_equal_eps_at_default = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = 1e-8, wtol = 0),
    tol_1ulp_literal_at_default_eps = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = 1e-11, wtol = 0),
    tol_perturbed_default = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def * (1 + 1e-12), wtol = 0),
    tol_zero_at_default = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = 0, wtol = 0),
    tol_negative = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = -tol_def, wtol = 0),
    tol_mismatched_from_default_to_custom_eps = list(epsilon = 1e-12, maxit = 100, trace = FALSE, tol = tol_def, wtol = 0),
    tol_uncapped_at_custom_eps = list(epsilon = 1e-12, maxit = 100, trace = FALSE, tol = 1e-7, wtol = 0),
    tol_perturbed_custom = list(epsilon = 1e-12, maxit = 100, trace = FALSE, tol = tol_c12 * (1 + 1e-12), wtol = 0),
    tol_missing_min_cap_for_large_eps = list(epsilon = 1e-2, maxit = 25, trace = FALSE, tol = 1e-5, wtol = 0),

    # Nondefault wtol values (catches mistakenly dropping/ignoring wtol or allowing weight dropping).
    wtol_tiny_pos = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = 1e-16),
    wtol_small_pos = list(epsilon = 1e-12, maxit = 100, trace = FALSE, tol = tol_c12, wtol = 1e-7),
    wtol_one = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = 1),
    wtol_neg = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = -1e-16),

    # Malformed tol and wtol types/shapes/nonfinites.
    tol_na = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = NA_real_, wtol = 0),
    tol_nan = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = NaN, wtol = 0),
    tol_inf = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = Inf, wtol = 0),
    tol_chr = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = "1e-11", wtol = 0),
    tol_lgl = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = FALSE, wtol = 0),
    tol_len2 = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = c(tol_def, tol_def), wtol = 0),
    tol_mat = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = matrix(tol_def, 1L, 1L), wtol = 0),
    tol_obj = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = structure(tol_def, class = "num_wrapper"), wtol = 0),
    wtol_na = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = NA_real_),
    wtol_nan = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = NaN),
    wtol_inf = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = Inf),
    wtol_lgl_false = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = FALSE),
    wtol_chr = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = "0"),
    wtol_len2 = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = c(0, 0)),
    wtol_mat = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = matrix(0, 1L, 1L)),
    wtol_obj = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = structure(0, class = "num_wrapper")),

    # Partial controls (lengths 0, 1, 2, 4).
    partial_empty = list(),
    partial_1 = list(epsilon = 1e-8),
    partial_2 = list(epsilon = 1e-8, maxit = 25),
    partial_4_only_tol = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def),
    partial_4_only_wtol = list(epsilon = 1e-8, maxit = 25, trace = FALSE, wtol = 0),
    partial_4_missing_trace = list(epsilon = 1e-8, maxit = 25, tol = tol_def, wtol = 0),

    # Duplicate fields across lengths 3, 4, 5, 6.
    dup_3 = list(epsilon = 1e-8, epsilon = 1e-8, maxit = 25),
    dup_4 = list(epsilon = 1e-8, maxit = 25, trace = FALSE, epsilon = 1e-8),
    dup_5_tol = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, tol = tol_def),
    dup_5_eps = list(epsilon = 1e-8, epsilon = 1e-8, maxit = 25, trace = FALSE, wtol = 0),
    dup_6 = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = 0, tol = tol_def),

    # Unknown fields and unnamed/partially named controls.
    unknown_3 = list(epsilon = 1e-8, maxit = 25, unknown = FALSE),
    unknown_4 = list(epsilon = 1e-8, maxit = 25, trace = FALSE, unused = 0),
    unknown_5_replace_wtol = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, unused = 0),
    unknown_5_replace_tol = list(epsilon = 1e-8, maxit = 25, trace = FALSE, unknown = tol_def, wtol = 0),
    unknown_6 = list(epsilon = 1e-8, maxit = 25, trace = FALSE, tol = tol_def, wtol = 0, unused = 0),
    unnamed_3 = list(1e-8, 25, FALSE),
    unnamed_5 = list(1e-8, 25, FALSE, tol_def, 0),
    blank_name_5 = stats::setNames(list(1e-8, 25, FALSE, tol_def, 0), c("epsilon", "maxit", "trace", "tol", "")),

    # Malformed epsilon, maxit, trace in 5-field shape.
    eps_nonpos_5 = list(epsilon = 0, maxit = 25, trace = FALSE, tol = 0, wtol = 0),
    eps_na_5 = list(epsilon = NA_real_, maxit = 25, trace = FALSE, tol = tol_def, wtol = 0),
    maxit_zero_5 = list(epsilon = 1e-8, maxit = 0, trace = FALSE, tol = tol_def, wtol = 0),
    maxit_frac_5 = list(epsilon = 1e-8, maxit = 25.5, trace = FALSE, tol = tol_def, wtol = 0),
    maxit_overflow_5 = list(epsilon = 1e-8, maxit = 2^31, trace = FALSE, tol = tol_def, wtol = 0),
    trace_na_5 = list(epsilon = 1e-8, maxit = 25, trace = NA, tol = tol_def, wtol = 0),
    trace_num_5 = list(epsilon = 1e-8, maxit = 25, trace = 0L, tol = tol_def, wtol = 0)
  )

  # Exercise short-circuit validation for every field in both supported shapes.
  # Single-bracket assignment deliberately retains NULL as a malformed scalar.
  c3 <- list(epsilon = 1e-8, maxit = 25L, trace = FALSE)
  c5 <- c(c3, list(tol = tol_def, wtol = 0))
  bad_numeric <- list(
    null = NULL, empty = numeric(0), list = list(1), pairlist = pairlist(1),
    environment = new.env(parent = emptyenv()), closure = function() 1,
    expression = expression(1), raw = as.raw(0), complex = 1 + 1i,
    logical = FALSE, character = "0", na = NA_real_, nan = NaN,
    pos_inf = Inf, neg_inf = -Inf, vector = c(1, 2),
    matrix = matrix(1, 1L, 1L), array = array(1, 1L),
    object = structure(1, class = "num_wrapper")
  )
  bad_trace <- c(bad_numeric[names(bad_numeric) != "logical"],
                list(empty_logical = logical(0), na_logical = NA,
                     vector_logical = c(TRUE, FALSE), matrix_logical = matrix(FALSE),
                     object_logical = structure(FALSE, class = "logical_wrapper")))
  for (shape in list(c3, c5)) {
    for (field in names(shape)) {
      values <- if (field == "trace") bad_trace else bad_numeric
      for (kind in names(values)) {
        bad <- shape
        bad[field] <- values[kind]
        bad_controls[[paste(length(shape), field, kind, sep = "_")]] <- bad
      }
    }
    for (field in names(shape)) {
      for (bad_name in list(NA_character_, "", paste0(field, "_unknown"),
                            substr(field, 1L, 1L))) {
        bad <- shape
        names(bad)[match(field, names(bad))] <- bad_name
        bad_controls[[paste(length(shape), field, "name", length(bad_controls), sep = "_")]] <- bad
      }
    }
    for (field in c("epsilon", "maxit")) {
      for (value in c(0, -1)) {
        bad <- shape
        bad[[field]] <- value
        bad_controls[[paste(length(shape), field, value, sep = "_")]] <- bad
      }
    }
    for (value in c(0.5, 25.5, 2^31)) {
      bad <- shape
      bad$maxit <- value
      bad_controls[[paste(length(shape), "maxit", value, sep = "_")]] <- bad
    }
  }

  for (nm in names(bad_controls)) {
    bad <- bad_controls[[nm]]
    err_v <- expect_error(cdid:::validate_ipw_glm_control(bad), class = "cdid_ipw_control_error", info = nm)
    expect_true(inherits(err_v, "cdid_ipw_input_error"), info = nm)
    expect_identical(conditionMessage(err_v), msg_invalid, info = nm)

    err_f <- expect_error(
      cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2, control = bad),
      class = "cdid_ipw_control_error",
      info = nm
    )
    expect_true(inherits(err_f, "cdid_ipw_input_error"), info = nm)
    expect_identical(conditionMessage(err_f), msg_invalid, info = nm)
  }
})

test_that("canonical 3-field and native-default 5-field controls yield identical fits, diagnostics, IF materialization, replay, and covariance", {
  set.seed(41029)
  n <- 80L
  Z <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n))
  w <- runif(n, 0.7, 1.8)
  G <- sample(c(3, 4, Inf), n, replace = TRUE, prob = c(0.35, 0.35, 0.30))
  obs <- matrix(TRUE, n, 3L)
  obs[sample(seq_len(n), 15L), 2L] <- FALSE
  obs[sample(seq_len(n), 14L), 3L] <- FALSE
  for (i in which(rowSums(obs) < 2L)) obs[i, 1:2] <- TRUE
  Y <- matrix(rnorm(n * 3L), n, 3L)

  panel <- build_cp02_control_fixture_panel(n, Y, obs, G, times = c(1, 2, 3))
  pidx <- cdid:::construct_observed_pairs(panel, "id", "time", "y", "g")
  ctx <- cdid:::construct_ipw_context(pidx, id = seq_len(n), Z = Z, weights = w)

  eps_specs <- list(
    default_eps = list(eps = 1e-8, maxit = 25L),
    custom_eps = list(eps = 1e-12, maxit = 100L)
  )

  for (nm in names(eps_specs)) {
    sp <- eps_specs[[nm]]
    eps <- sp$eps
    maxit <- sp$maxit
    tol <- min(1e-7, eps / 1000)

    ctrl3 <- list(epsilon = eps, maxit = maxit, trace = FALSE)
    ctrl5 <- list(epsilon = eps, maxit = maxit, trace = FALSE, tol = tol, wtol = 0)
    ctrl5_reord <- list(wtol = 0L, tol = tol, trace = FALSE, maxit = as.double(maxit), epsilon = eps)
    ctrl_native <- stats::glm.control(epsilon = eps, maxit = maxit, trace = FALSE)

    fit12_3 <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2, control = ctrl3, link_id = "12")
    fit13_3 <- cdid:::estimate_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 3, control = ctrl3, link_id = "13")
    fits_3 <- list("12" = fit12_3, "13" = fit13_3)

    expect_identical(fit12_3$diagnostics$control, list(epsilon = as.double(eps), maxit = as.integer(maxit), trace = FALSE))

    mat_3 <- cdid:::materialize_ipw_if(ctx, fits_3, max_entries = n * 2L, chunk_size = 13L, include_components = TRUE)
    cov_3 <- cdid:::stream_ipw_covariance(ctx, fits_3, chunk_size = 13L)
    batches_3 <- list()
    info_3 <- cdid:::replay_ipw_if(ctx, fits_3, chunk_size = 13L, include_components = TRUE, callback = function(b) {
      batches_3[[length(batches_3) + 1L]] <<- b
    })

    for (alt_ctrl in list(ctrl5, ctrl5_reord, ctrl_native)) {
      fit12_alt <- cdid:::fit_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 2, control = alt_ctrl, link_id = "12")
      fit13_alt <- cdid:::estimate_ipw_link(ctx, g = 3, K = c(4, Inf), s = 1, t = 3, control = alt_ctrl, link_id = "13")
      fits_alt <- list("12" = fit12_alt, "13" = fit13_alt)

      expect_identical(fit12_alt, fit12_3, info = nm)
      expect_identical(fit13_alt, fit13_3, info = nm)

      mat_alt <- cdid:::materialize_ipw_if(ctx, fits_alt, max_entries = n * 2L, chunk_size = 13L, include_components = TRUE)
      expect_identical(mat_alt, mat_3, info = nm)

      cov_alt <- cdid:::stream_ipw_covariance(ctx, fits_alt, chunk_size = 13L)
      expect_identical(cov_alt, cov_3, info = nm)

      batches_alt <- list()
      info_alt <- cdid:::traverse_ipw_if(ctx, fits_alt, chunk_size = 13L, include_components = TRUE, callback = function(b) {
        batches_alt[[length(batches_alt) + 1L]] <<- b
      })
      expect_identical(info_alt, info_3, info = nm)
      expect_identical(batches_alt, batches_3, info = nm)
    }
  }

  # Require the same immutable fixture and loader as the existing IPW tests.
  # A missing fixture must fail rather than silently bypass this coverage.
  {
    indep_path <- normalizePath(testthat::test_path(
      "reference", "cp02_weighted_pair_if_independent_reference.R"
    ), winslash = "/", mustWork = TRUE)
    ref <- new.env(parent = globalenv())
    ref$commandArgs <- function(...) paste0("--file=", indep_path)
    ref$write.csv <- function(...) invisible(NULL)
    ref$sys.source <- ipw_reference_source
    expect_identical(digest::digest(file = indep_path, algo = "sha256"),
                     "8815f60a292c61c0391e8ee79ccac0df75a3615fb00bc7bbb9f5c9aaa24fe0e7")
    capture.output(sys.source(indep_path, envir = ref))

    panel_ref <- build_cp02_control_fixture_panel(ref$n, ref$Y, ref$observed, ref$G, times = c(1, 2, 3, 4))
    pidx_ref <- cdid:::construct_observed_pairs(panel_ref, "id", "time", "y", "g")
    ctx_ref <- cdid:::construct_ipw_context(pidx_ref, id = seq_len(ref$n), Z = ref$Z)

    for (nm in names(eps_specs)) {
      sp <- eps_specs[[nm]]
      eps <- sp$eps
      maxit <- sp$maxit
      tol <- min(1e-7, eps / 1000)
      c3 <- list(epsilon = eps, maxit = maxit, trace = FALSE)
      c5 <- list(epsilon = eps, maxit = maxit, trace = FALSE, tol = tol, wtol = 0)

      f12_3 <- cdid:::fit_ipw_link(ctx_ref, g = 3, K = c(4, Inf), s = 1, t = 2, control = c3, link_id = "12")
      f13_3 <- cdid:::fit_ipw_link(ctx_ref, g = 3, K = c(4, Inf), s = 1, t = 3, control = c3, link_id = "13")
      f12_5 <- cdid:::fit_ipw_link(ctx_ref, g = 3, K = c(4, Inf), s = 1, t = 2, control = c5, link_id = "12")
      f13_5 <- cdid:::fit_ipw_link(ctx_ref, g = 3, K = c(4, Inf), s = 1, t = 3, control = c5, link_id = "13")

      expect_identical(f12_5, f12_3, info = nm)
      expect_identical(f13_5, f13_3, info = nm)
      expect_identical(
        cdid:::materialize_ipw_if(ctx_ref, list("12" = f12_5, "13" = f13_5), max_entries = 480, chunk_size = 13, include_components = TRUE),
        cdid:::materialize_ipw_if(ctx_ref, list("12" = f12_3, "13" = f13_3), max_entries = 480, chunk_size = 13, include_components = TRUE),
        info = nm
      )
    }
  }
})
