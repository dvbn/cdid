# Independent raw nuisance refits plus the unchanged general separate route.
fop_ns <- function(name) getFromNamespace(name, "cdid")
fop_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
fop_case <- function(chunk = 7L, intercept = TRUE, slopes = TRUE, alias = FALSE,
                     bad_alias = FALSE, missing = FALSE, large_level = FALSE) {
  j <- rep(1:8, 4); k <- rep(1:4, each = 8); N <- length(j); id <- seq_len(N)
  g <- c(3, 4, 5, Inf)[k]; calendar <- 1:5
  x <- (j - 4.5) / 3; z <- sin(j) + j / 20
  w <- (1 + (j %% 3) / 5) * (1 + k * x / 40)
  Y <- vapply(calendar, function(d) j / 7 + k * d / 4 +
    (.2 * d + .08 * d^2 + .03 * k * d) * x + sin(j + k * d) / 8 +
    as.double(d >= g) * (d + x * k / 5), numeric(N))
  if (large_level) Y <- Y + 1e9 * (1 + x)
  if (missing) Y[j == 1L, 2L] <- NA_real_
  X <- matrix(numeric(), N, 0L)
  if (intercept) X <- cbind(`(Intercept)` = rep(1, N))
  if (slopes) X <- cbind(X, x = x)
  if (alias) X <- cbind(X, alias = 3 + 2 * x + if (bad_alias) as.double(k == 1L & j == 1L) / 10 else 0)
  V <- cbind(`(Intercept)` = 1, z = z)
  Xa <- cbind(`(Intercept)` = 1, x = x, z = z)
  if (alias) Xa <- cbind(Xa, alias = X[, "alias"])
  ctrl <- fop_ns("ipw_owner_control")(chunk, 1e8, 2e8, 1e13, 1e8)
  long <- data.frame(id = rep(id, each = 5L), time = rep(calendar, N),
    y = as.vector(t(Y)), g = rep(g, each = 5L))
  panel <- fop_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  frame <- fop_ns("construct_ipw_frame_owner")(panel, id, Xa, w, ctrl)
  reader <- fop_ns("construct_cp04_frame_reader")(frame, id,
    matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, ctrl)
  def <- list(c(1, 2, 3), c(1, 3, 3), c(2, 3, 3), c(2, 4, 3),
    c(1, 3, 3), c(1, 4, 4), c(3, 4, 4), c(2, 5, 5), c(1, 3, 3), c(2, 3, 3), c(1, 4, 4))
  pools <- list(c(4, 5, Inf), c(4, 5, Inf), c(4, 5, Inf), c(5, Inf),
    c(4, 5), c(5, Inf), Inf, Inf, Inf, Inf, Inf)
  views <- lapply(seq_along(def), function(i) fop_ns("construct_ipw_link_view")(frame,
    def[[i]][1L], def[[i]][2L], def[[i]][3L], pools[[i]], control = ctrl,
    link_id = paste0("fixed-operator-", i), functional_key = paste0("cmpcat:fixed-operator-", i)))
  recipe <- fop_ns("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  # OR requires exactly its own frame recipe; DR permits a distinct treatment design.
  or_frame <- fop_ns("construct_ipw_frame_owner")(panel, id, X, w, ctrl)
  or_reader <- fop_ns("construct_cp04_frame_reader")(or_frame, id,
    matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, ctrl)
  or_views <- lapply(seq_along(def), function(i) fop_ns("construct_ipw_link_view")(or_frame,
    def[[i]][1L], def[[i]][2L], def[[i]][3L], pools[[i]], control = ctrl,
    link_id = paste0("fixed-operator-", i), functional_key = paste0("cmpcat:fixed-operator-", i)))
  list(reader = reader, or_reader = or_reader, views = views, or_views = or_views,
    control = ctrl, recipe = recipe, X = X, V = V, Y = Y, g = g, w = w,
    long = long, frame = data.frame(id = id, g = ifelse(is.finite(g), g, 0), x = x, z = z, w = w))
}

fop_fit <- function(a, method, fast = TRUE) {
  fn <- fop_ns(if (method == "reg") "fit_or_baseline" else "fit_dr_baseline")
  if (!fast) {
    env <- new.env(parent = environment(fn))
    env$or_separate_complete_facts <- function(...) NULL
    environment(fn) <- env
  }
  args <- list(reader = if (method == "reg") a$or_reader else a$reader,
    link_views = if (method == "reg") a$or_views else a$views, control = a$control,
    rank_tol = 1e-10, outcome_recipe = a$recipe, outcome_pooling = "separate")
  if (method == "dr") args <- c(args, list(glm_control = fop_gc, treatment_columns = colnames(a$V)))
  do.call(fn, args)
}

# Refit every original nuisance on f_i*w_i, independently of the optimization.
fop_raw <- function(a, fit, method, frequency = rep(1, length(a$w))) {
  N <- length(frequency); w <- a$w * frequency
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  raw <- lapply(fit$support_specs, function(s) {
    ix <- which(a$g == st$cohorts[s$cohort])
    list(ix = ix, dy = a$Y[ix, s$pair$t_calendar_index] - a$Y[ix, s$pair$s_calendar_index], L = s$transitions)
  })
  theta <- lapply(seq_along(fit$regression_specs), function(b) {
    keep <- fit$outcome_groups$groups[[b]]$keep
    ss <- raw[fit$regression_specs[[b]]$support_ids]
    design <- lapply(ss, function(s) {
      D <- s$L * a$X[s$ix, a$recipe$slope_index[keep], drop = FALSE]
      if (a$recipe$intercept) cbind(1, D) else D
    })
    unname(stats::lm.wfit(do.call(rbind, design), unlist(lapply(ss, `[[`, "dy")),
      unlist(lapply(ss, function(s) w[s$ix])), tol = 1e-10)$coefficients)
  })
  residual <- function(sid, b) {
    s <- raw[[sid]]; keep <- fit$outcome_groups$groups[[b]]$keep
    D <- s$L * a$X[s$ix, a$recipe$slope_index[keep], drop = FALSE]
    if (a$recipe$intercept) D <- cbind(1, D)
    s$dy - as.vector(D %*% theta[[b]])
  }
  if (method == "reg") return(vapply(fit$links, function(l) {
    s <- raw[[l$support_id]]; weighted.mean(residual(l$support_id, l$regression), w[s$ix])
  }, numeric(1)))
  models <- vapply(fit$treatment_specs, function(t) {
    ids <- t$support_ids; rr <- unlist(lapply(raw[ids], `[[`, "ix"), use.names = FALSE)
    D <- as.double(a$g[rr] == st$cohorts[t$cohorts[1L]])
    root <- suppressWarnings(stats::glm.fit(a$V[rr, , drop = FALSE], D, weights = w[rr],
      start = rep(0, ncol(a$V)), family = stats::binomial(), control = fop_gc, intercept = FALSE))
    stopifnot(root$converged, all(is.finite(root$coefficients)))
    p <- pmin(stats::plogis(as.vector(a$V[rr, , drop = FALSE] %*% root$coefficients)), 1 - 1e-6)
    r <- unlist(lapply(ids, function(s) residual(s, t$regression)), use.names = FALSE)
    wp <- w[rr] * D; wm <- w[rr] * (1 - D) * p / (1 - p)
    weighted.mean(r, wp) - weighted.mean(r, wm)
  }, numeric(1))
  models[vapply(fit$links, `[[`, integer(1), "model")]
}

fop_eval <- function(fit, method, rows, links, details = FALSE, scale = "frame")
  fop_ns(if (method == "reg") "evaluate_or_baseline_if_chunk" else "evaluate_dr_if_chunk")(
    fit, rows, links, scale, details)

testthat::test_that("fixed operators retain every multi-pool block and original-unit derivative", {
  a <- fop_case(); N <- length(a$w)
  for (method in c("reg", "dr")) {
    fit <- fop_fit(a, method); slow <- fop_fit(a, method, FALSE); L <- length(fit$links)
    testthat::expect_null(fit$complete_support)
    testthat::expect_identical(fit$accounting$fixed_operator$route, "separate_complete")
    testthat::expect_lt(fit$accounting$fixed_operator$design_qr, fit$accounting$fixed_operator$response_columns)
    testthat::expect_equal(fit$regression_specs, slow$regression_specs)
    testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), fop_raw(a, fit, method), tolerance = 1e-9)
    testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), vapply(slow$links, `[[`, numeric(1), "mu"), tolerance = 1e-9)
    for (b in seq_along(fit$regressions)) {
      bf <- fop_ns("or_outcome_block_fit")(fit$outcome_groups, b)
      bs <- fop_ns("or_outcome_block_fit")(slow$outcome_groups, b)
      testthat::expect_identical(bf$keep, bs$keep)
      testthat::expect_equal(c(bf$alpha, bf$beta), c(bs$alpha, bs$beta), tolerance = 1e-9)
    }
    actual <- matrix(0, N, L); reference <- actual
    for (start in seq.int(1L, N, by = 7L)) {
      rows <- seq.int(start, min(N, start + 6L))
      plain <- fop_eval(fit, method, rows, seq_len(L))
      detailed <- fop_eval(fit, method, rows, seq_len(L), TRUE)
      general <- fop_eval(slow, method, rows, seq_len(L))
      actual[rows, ] <- plain$phi; reference[rows, ] <- general$phi
      testthat::expect_equal(plain$phi, detailed$phi, tolerance = 1e-10)
      testthat::expect_equal(plain$phi, general$phi, tolerance = 1e-9)
      testthat::expect_lte(unname(fop_ns("cp04_payload")(plain)["entries"]), plain$diagnostics$output_entries_bound)
    }
    numerical <- t(vapply(seq_len(N), function(i) {
      direction <- rep(-1, N); direction[i] <- N - 1
      (fop_raw(a, fit, method, 1 + 1e-6 * direction) -
        fop_raw(a, fit, method, 1 - 1e-6 * direction)) / 2e-6
    }, numeric(L)))
    testthat::expect_equal(unname(actual), unname(numerical), tolerance = 2e-6)
    testthat::expect_equal(crossprod(actual) / N, crossprod(reference) / N, tolerance = 1e-9)
    subset <- fop_eval(fit, method, 7:1, c(11L, 3L, 1L), scale = "core")
    testthat::expect_equal(unname(subset$phi), actual[7:1, c(11L, 3L, 1L)], tolerance = 1e-10)
    if (method == "dr") {
      nf <- length(unique(vapply(fit$treatment_specs, function(s) paste(s$cohorts, collapse = ":"), character(1))))
      testthat::expect_equal(fit$accounting$fixed_operator$n_treatment_logit_fits, nf)
      testthat::expect_equal(fit$accounting$fixed_operator$n_treatment_logit_reuses, length(fit$treatments) - nf)
      testthat::expect_equal(lapply(fit$treatments, `[[`, "partials"), lapply(slow$treatments, `[[`, "partials"), tolerance = 1e-9)
    }
  }
})

testthat::test_that("fixed QR preserves raw differences, aliases, intercept-only fits and chunking", {
  for (spec in list(list(large_level = TRUE), list(alias = TRUE), list(slopes = FALSE))) {
    a <- do.call(fop_case, c(list(chunk = 3L), spec))
    b <- do.call(fop_case, c(list(chunk = 11L), spec))
    for (method in c("reg", "dr")) {
      fit <- fop_fit(a, method); slow <- fop_fit(a, method, FALSE); other <- fop_fit(b, method)
      testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), vapply(slow$links, `[[`, numeric(1), "mu"), tolerance = 1e-9)
      for (rows in list(1:3, 9:7, 32:30)) {
        testthat::expect_equal(fop_eval(fit, method, rows, c(5L, 1L))$phi,
          fop_eval(slow, method, rows, c(5L, 1L))$phi, tolerance = 1e-9)
        testthat::expect_equal(fop_eval(fit, method, rows, c(5L, 1L))$phi,
          fop_eval(other, method, rows, c(5L, 1L))$phi, tolerance = 1e-9)
      }
    }
  }
  a <- fop_case()
  testthat::expect_null(fop_ns("or_separate_complete_facts")(a$reader, a$recipe,
    list(), list(), a$control, reference = list()))
  bad <- fop_case(alias = TRUE, bad_alias = TRUE)
  testthat::expect_error(fop_fit(bad, "reg"), class = "cdid_or_wls_rank_error")
  for (a in list(fop_case(missing = TRUE), fop_case(intercept = FALSE))) for (method in c("reg", "dr"))
    testthat::expect_null(fop_fit(a, method)$separate_complete)
})

fop_public <- function(a, method, law, pt, slow = FALSE, tuning = list()) {
  if (slow) testthat::local_mocked_bindings(or_separate_complete_facts = function(...) NULL, .package = "cdid")
  data <- a$long; data$g[!is.finite(data$g)] <- 0
  args <- list(yname = "y", tname = "time", idname = "id", gname = "g", data = data,
    frame = a$frame, xformla = ~x + z, weightsname = "w", selection = "none",
    est_method = method, control_group = law, pt_assumption = pt, outcome_pooling = "separate",
    bstrap = FALSE, cband = FALSE, gmm_weighting = "optimal", control = c(list(chunk_rows = 7L), tuning))
  if (pt == "all") args$pt_start <- 1
  do.call(cdid::att_gt, args)
}

testthat::test_that("fixed operators keep the full public post/all GMM criterion and inference", {
  a <- fop_case()
  cases <- list(c("reg", "nevertreated", "post"), c("reg", "notyettreated", "all"),
    c("dr", "notyettreated", "post"), c("dr", "nevertreated", "all"))
  for (cs in cases) {
    tuning <- if (cs[3L] == "all") list(condition_trigger = .99, ridge_relative = .5) else
      if (cs[1L] == "dr") list(condition_trigger = 1e-20) else list()
    fast <- fop_public(a, cs[1L], cs[2L], cs[3L], tuning = tuning)
    slow <- fop_public(a, cs[1L], cs[2L], cs[3L], TRUE, tuning)
    testthat::expect_equal(fast$att, slow$att, tolerance = 1e-9)
    testthat::expect_equal(fast$inffunc, slow$inffunc, tolerance = 1e-9)
    testthat::expect_equal(fast$V_analytical, slow$V_analytical, tolerance = 1e-9)
    testthat::expect_equal(fast$se, slow$se, tolerance = 1e-9)
    testthat::expect_identical(fast$chained$available, slow$chained$available)
    hf <- fast$chained$covariance_handling; hs <- slow$chained$covariance_handling
    testthat::expect_identical(hf$handler_invoked, hs$handler_invoked)
    testthat::expect_identical(hf$handler_used, hs$handler_used)
    testthat::expect_equal(hf$lambda, hs$lambda, tolerance = 1e-9)
    if (cs[3L] == "all") {
      testthat::expect_identical(hf$handler_used, "tikhonov")
      testthat::expect_gt(hf$lambda, 0)
    }
  }
})

testthat::test_that("many consumers reserve non-row adjoint matrices before ordinary replay", {
  n <- 32L; p <- 24L; consumers <- 400L; N <- 2L * n; id <- seq_len(N)
  j <- rep(seq_len(n), 2L); g <- rep(c(3, Inf), each = n)
  Z <- cos(outer(j - .5, seq_len(p), function(x, k) pi * x * k / n))
  colnames(Z) <- paste0("x", seq_len(p)); X <- cbind(`(Intercept)` = 1, Z)
  w <- 1 + (j %% 3L) / 10
  Y <- vapply(1:3, function(d) as.vector(Z %*% (d / (seq_len(p) + 1))) +
    sin(j + d) / 8 + as.double(g == 3) * d, numeric(N))
  control <- fop_ns("ipw_owner_control")(1L, 1e8, 2e8, 1e13, 1e8)
  data <- data.frame(id = rep(id, each = 3L), time = rep(1:3, N),
    y = as.vector(t(Y)), g = rep(g, each = 3L))
  panel <- fop_ns("construct_ipw_panel_owner")(data, "id", "time", "y", "g", control)
  frame <- fop_ns("construct_ipw_frame_owner")(panel, id, X, w, control)
  reader <- fop_ns("construct_cp04_frame_reader")(frame, id,
    matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, control)
  views <- lapply(seq_len(consumers), function(i) fop_ns("construct_ipw_link_view")(frame,
    1, 3, 3, Inf, control = control, link_id = sprintf("h%03d", i), functional_key = "cmpcat:many-consumers"))
  recipe <- fop_ns("or_pooled_outcome_recipe")(colnames(X), TRUE, colnames(Z))
  fit <- fop_ns("fit_or_baseline")(reader, views, control, 1e-10, recipe, outcome_pooling = "separate")
  testthat::expect_length(fit$regressions, 1L)
  testthat::expect_length(fit$supports, 2L)
  testthat::expect_equal(fit$outcome_groups$groups[[1L]]$keep, seq_len(p))
  pre <- fop_ns("or_baseline_if_budget")(fit, 1L, seq_len(consumers), FALSE)
  testthat::expect_equal(pre$adjoint_entries_peak, p * consumers)
  # Loadings, Gram RHS/solutions and their solve copies do not have a B
  # multiplier. For this high-p/small-B case they exceed the old row bound.
  nonrow_reservation <- 32 * p * consumers
  testthat::expect_gte(pre$scratch, nonrow_reservation)
  result <- fop_ns("evaluate_or_baseline_if_chunk")(fit, 1L, rev(seq_len(consumers)))
  testthat::expect_equal(unname(result$phi), matrix(result$phi[1L, 1L], 1L, consumers), tolerance = 1e-12)
  testthat::expect_lte(unname(fop_ns("cp04_payload")(result)["entries"]), pre$output)
  small <- fop_ns("or_baseline_if_budget")(fit, 1L, c(400L, 1L, 2L), FALSE)
  testthat::expect_equal(small$adjoint_entries_peak, 3L * p)
  testthat::expect_gte(pre$scratch - small$scratch, 32 * p * (consumers - 3L))
  # A pure envelope witness, never fitted/evaluated as a modified backend.
  # Its independent cap cannot even fund the non-row solve reservation alone.
  metadata <- unclass(fit)
  retained <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  metadata$control <- fop_ns("ipw_owner_control")(1L, 1e8,
    retained + nonrow_reservation - 1, 1e13, 1e8)
  testthat::expect_error(fop_ns("or_baseline_separate_if_budget")(
    metadata, 1L, seq_len(consumers), FALSE), class = "cdid_ipw_owner_bound_error")
})

testthat::test_that("online common QR preserves each raw support's full all-coefficient objective", {
  for (chunk in c(1L, 3L, 11L)) {
    a <- fop_case(chunk = chunk); fit <- fop_fit(a, "reg")
    prepared <- fop_ns("or_separate_complete_supports")(a$or_reader, fit$support_specs,
      a$recipe, seq_along(a$recipe$columns), a$control, fit$separate_complete,
      fit$accounting$fit_scratch_bound)
    st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
    if (chunk < 8L) testthat::expect_gt(prepared$counts["common_design_merges"], 0)
    for (sid in seq_along(fit$support_specs)) {
      sp <- fit$support_specs[[sid]]; ix <- which(a$g == st$cohorts[sp$cohort])
      dy <- a$Y[ix, sp$pair$t_calendar_index] - a$Y[ix, sp$pair$s_calendar_index]
      D <- cbind(1, sp$transitions * a$X[ix, "x"])
      F <- prepared$supports[[sid]]$design$summary$raw$factor
      for (theta in list(c(-.3, .4), c(1.7, -.2))) {
        raw_objective <- sum(a$w[ix] * (dy - as.vector(D %*% theta))^2)
        compact_objective <- sum(as.vector(F %*% c(-theta, 1))^2)
        testthat::expect_equal(compact_objective, raw_objective, tolerance = 1e-9)
      }
    }
  }
})

testthat::test_that("general separate replay groups only requested positions and preserves frozen fits", {
  a <- fop_case(missing = TRUE); N <- length(a$w)
  for (method in c("reg", "dr")) {
    fit <- fop_fit(a, method); testthat::expect_null(fit$separate_complete)
    frozen <- list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, support_specs = fit$support_specs)
    L <- length(fit$links); full <- matrix(0, N, L)
    for (start in seq.int(1L, N, by = 7L)) {
      rows <- seq.int(start, min(N, start + 6L))
      full[rows, ] <- fop_eval(fit, method, rows, seq_len(L))$phi
    }
    rows <- 7:1; links <- c(11L, 3L, 1L)
    plain <- fop_eval(fit, method, rows, links)
    detailed <- fop_eval(fit, method, rows, links, TRUE)
    testthat::expect_equal(unname(plain$phi), full[rows, links, drop = FALSE], tolerance = 1e-10)
    testthat::expect_equal(plain$phi, detailed$phi, tolerance = 1e-10)
    testthat::expect_equal(detailed$details$link_block_ids,
      vapply(fit$links[links], `[[`, integer(1), "regression"))
    testthat::expect_identical(list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, support_specs = fit$support_specs), frozen)
  }
})

# V14: explicit native coverage. Private evaluation without a plan remains
# the exact accepted scalar complete producer and independent test reference.
v14_eval <- function(fit, method, rows, links, plan = NULL, scale = "frame", details = FALSE) {
  fop_ns(if (method == "reg") "evaluate_or_baseline_if_chunk" else "evaluate_dr_baseline_if_chunk")(
    fit, rows, links, scale = scale, details = details, plan = plan)
}
v14_call <- function(fit, rows, plan, data = plan$data, cap = NULL, X = NULL) {
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  panel <- fop_ns("or_baseline_complete_panel")(reader, fit$separate_complete)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE); post <- panel$postings
  if (is.null(X)) X <- frame$Z_a
  if (is.null(cap)) cap <- fop_ns("separate_complete_native_live")(length(rows), length(plan$links),
    length(plan$dates), length(plan$gids), ncol(data$family))
  .Call(fop_ns("or_baseline_native_kernel")("C_gmm_separate_complete_if"), panel$observations$y,
    post$date_start, post$post_a_index, post$post_obs_index, as.integer(rows), reader$state$cohort_index,
    reader$state$a_index, reader$state$w, plan$dates, X, data, as.double(cap))
}

testthat::test_that("complete native plans preserve every raw moment IF and independent refit derivative", {
  a <- fop_case(chunk = 32L); N <- length(a$w)
  for (method in c("reg", "dr")) {
    fit <- fop_fit(a, method); links <- seq_along(fit$links)
    plan <- fop_ns("chained_full_if_plan")(fit, links, N)
    testthat::expect_false(is.null(plan))
    testthat::expect_identical(plan$kind, "separate_complete_native_plan")
    frozen <- list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, support_specs = fit$support_specs)
    actual <- v14_eval(fit, method, seq_len(N), links, plan)
    scalar <- v14_eval(fit, method, seq_len(N), links)
    detailed <- v14_eval(fit, method, seq_len(N), links, details = TRUE)
    testthat::expect_identical(actual$diagnostics$execution, "native_fixed_operator")
    testthat::expect_equal(actual$phi, scalar$phi, tolerance = 1e-10)
    testthat::expect_equal(actual$phi, detailed$phi, tolerance = 1e-10)
    numerical <- t(vapply(seq_len(N), function(i) {
      direction <- rep(-1, N); direction[i] <- N - 1
      (fop_raw(a, fit, method, 1 + 1e-6 * direction) -
        fop_raw(a, fit, method, 1 - 1e-6 * direction)) / 2e-6
    }, numeric(length(links))))
    testthat::expect_equal(unname(actual$phi), unname(numerical), tolerance = 2e-6)
    testthat::expect_equal(crossprod(actual$phi) / N, crossprod(scalar$phi) / N, tolerance = 1e-9)
    reordered <- c(11L, 3L, 1L)
    subset_plan <- fop_ns("chained_full_if_plan")(fit, reordered, N)
    testthat::expect_false(is.null(subset_plan))
    subset <- v14_eval(fit, method, N:1, reordered, subset_plan, scale = "core")
    testthat::expect_equal(subset$phi, scalar$phi[N:1, reordered, drop = FALSE], tolerance = 1e-10)
    testthat::expect_identical(list(groups = fit$outcome_groups, regressions = fit$regressions,
      treatments = fit$treatments, support_specs = fit$support_specs), frozen)
    testthat::expect_lte(unname(fop_ns("cp04_payload")(plan)["entries"]), plan$retained_bound)
    pre <- fop_ns(if (method == "reg") "or_baseline_if_budget" else "dr_baseline_if_budget")(fit, N, links, FALSE)
    live <- fop_ns("separate_complete_native_live")(N, length(links), length(plan$dates),
      length(plan$gids), ncol(plan$data$family))
    testthat::expect_lte(plan$entries + live, pre$scratch)
    testthat::expect_error(v14_call(fit, seq_len(N), plan, cap = live - 1), class = "cdid_gmm_bound_error")
    # This concrete loading mutation would fail the independent IF check.
    wrong <- plan$data; wrong$adj <- -wrong$adj
    changed <- v14_call(fit, seq_len(N), plan, wrong)
    testthat::expect_identical(changed$status[1L], 0L)
    testthat::expect_gt(max(abs(changed$phi - numerical)), 2e-6)
    bad <- plan; bad$data$adj[1L] <- bad$data$adj[1L] + 1
    testthat::expect_error(v14_eval(fit, method, seq_len(N), links, bad), class = "cdid_or_baseline_integrity_error")
    testthat::expect_error(v14_eval(fit, method, seq_len(N), rev(links), plan), class = "cdid_or_baseline_integrity_error")
    other <- fop_fit(a, method)
    testthat::expect_error(v14_eval(other, method, seq_len(N), links, plan), class = "cdid_or_baseline_integrity_error")
    testthat::expect_error(v14_eval(fit, method, seq_len(N - 1L), links, plan), class = "cdid_or_baseline_integrity_error")
    changed_fit <- fit
    changed_fit$outcome_groups$groups[[1L]]$beta[1L] <- changed_fit$outcome_groups$groups[[1L]]$beta[1L] + 1
    testthat::expect_error(v14_eval(changed_fit, method, seq_len(N), links, plan), class = "cdid_ipw_fit_error")
    reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
    X <- get("frame_owner", reader$frame_ref, inherits = FALSE)$Z_a
    X[reader$state$a_index[1L], 2L] <- NA_real_
    testthat::expect_error(v14_call(fit, seq_len(N), plan, X = X), class = "cdid_gmm_bound_error")
    bad_data <- plan$data; bad_data$block[3L, 1L] <- NA_integer_
    testthat::expect_error(v14_call(fit, seq_len(N), plan, bad_data), class = "cdid_gmm_bound_error")
  }
})

# Source-owned frames/views for two different retained local bases, and a
# broad common control family that can fund a genuine one-row tail.
v14_case <- function(tail = FALSE) {
  if (tail) {
    k <- c(rep(seq_len(61L), each = 4L), rep(1L, 13L))
    j <- c(rep(1:4, 61L), rep(1:4, length.out = 13L))
    cohorts <- c(3, 4:62, Inf); chunk <- 256L
  } else {
    k <- rep(1:4, each = 8L); j <- rep(1:8, 4L)
    cohorts <- c(3, 4, 5, Inf); chunk <- 32L
  }
  N <- length(j); id <- seq_len(N); g <- cohorts[k]
  x <- (j - mean(j)) / 3; z <- sin(j) + j / 20
  v <- ifelse(k %in% c(1L, length(cohorts)), 0, cos(j))
  w <- 1 + (j %% 3) / 5
  X <- cbind(`(Intercept)` = 1, x = x)
  if (!tail) X <- cbind(X, v = v)
  V <- cbind(`(Intercept)` = 1, z = z); Xa <- cbind(X, z = z)
  Y <- vapply(1:5, function(d) j / 7 + k * d / 80 + (.2 * d + .08 * d^2) * x +
    d^2 * v / 30 + sin(j + d) / 8 + as.double(d >= g) * (d + x / 5), numeric(N))
  ctrl <- fop_ns("ipw_owner_control")(chunk, 1e8, 2e8, 1e13, 1e8)
  long <- data.frame(id = rep(id, each = 5L), time = rep(1:5, N), y = as.vector(t(Y)), g = rep(g, each = 5L))
  panel <- fop_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  build <- function(Z) {
    frame <- fop_ns("construct_ipw_frame_owner")(panel, id, Z, w, ctrl)
    reader <- fop_ns("construct_cp04_frame_reader")(frame, id,
      matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, ctrl)
    def <- if (tail) list(c(1, 3, 3)) else list(c(1, 3, 3), c(1, 4, 4))
    pools <- if (tail) list(cohorts[-1L]) else list(Inf, c(5, Inf))
    views <- lapply(seq_along(def), function(i) fop_ns("construct_ipw_link_view")(frame,
      def[[i]][1L], def[[i]][2L], def[[i]][3L], pools[[i]], control = ctrl,
      link_id = paste0("native-local-", i), functional_key = paste0("cmpcat:native-local-", i)))
    list(reader = reader, views = views)
  }
  a <- build(X); b <- build(Xa)
  list(reader = b$reader, or_reader = a$reader, views = b$views, or_views = a$views,
    control = ctrl, recipe = fop_ns("or_pooled_outcome_recipe")(colnames(X), TRUE, colnames(X)[-1L]),
    X = X, V = V, Y = Y, g = g, w = w)
}

testthat::test_that("native complete replay respects local bases, zero slopes and admitted tails", {
  for (a in list(v14_case(), v14_case(TRUE), fop_case(chunk = 32L, slopes = FALSE),
                 fop_case(chunk = 32L, alias = TRUE), fop_case(chunk = 32L, large_level = TRUE))) {
    N <- length(a$w); B <- min(N, a$control$max_chunk_rows)
    for (method in c("reg", "dr")) {
      fit <- fop_fit(a, method); links <- seq_along(fit$links)
      plan <- fop_ns("chained_full_if_plan")(fit, links, B)
      testthat::expect_false(is.null(plan))
      if (N == 257L) testthat::expect_equal(plan$min_B, 1)
      if (ncol(a$X) == 3L && !"alias" %in% colnames(a$X))
        testthat::expect_gt(length(unique(vapply(fit$outcome_groups$groups, function(g) length(g$keep), integer(1)))), 1L)
      for (start in seq.int(1L, N, by = B)) {
        rows <- seq.int(start, min(N, start + B - 1L))
        actual <- v14_eval(fit, method, rows, links, plan)
        scalar <- v14_eval(fit, method, rows, links)
        testthat::expect_identical(actual$diagnostics$execution, "native_fixed_operator")
        testthat::expect_equal(actual$phi, scalar$phi, tolerance = 1e-9)
      }
    }
  }
})

testthat::test_that("optional native capacity decline keeps original admission and scalar replay", {
  a <- fop_case(chunk = 31L)
  for (method in c("reg", "dr")) {
    fit <- fop_fit(a, method); links <- seq_along(fit$links)
    testthat::expect_null(fop_ns("chained_full_if_plan")(fit, links, 31L))
    before <- fop_ns(if (method == "reg") "or_baseline_if_budget" else "dr_baseline_if_budget")(fit, 1L, links, FALSE)
    out <- fop_ns("chained_full_if_chunk")(fit, 32L, links, "frame", plan = NULL)
    testthat::expect_identical(out$diagnostics$execution, "fixed_operator")
    testthat::expect_equal(out$phi, v14_eval(fit, method, 32L, links)$phi, tolerance = 1e-10)
    after <- fop_ns(if (method == "reg") "or_baseline_if_budget" else "dr_baseline_if_budget")(fit, 1L, links, FALSE)
    testthat::expect_identical(before, after)
  }
  for (a in list(fop_case(chunk = 32L, missing = TRUE), fop_case(chunk = 32L, intercept = FALSE))) {
    fit <- fop_fit(a, "reg")
    testthat::expect_null(fop_ns("chained_full_if_plan")(fit, seq_along(fit$links), 32L))
  }
})
