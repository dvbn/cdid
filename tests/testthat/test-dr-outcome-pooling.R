# Independent raw-data nuisance/refit references for private separate DR.
drp_ns <- function(name) getFromNamespace(name, "cdid")
drp_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)

drp_case <- function(universe = "F", intercept = TRUE, slopes = TRUE, dated = FALSE, route_kind = "corollary1", extra_reference = FALSE) {
  grid <- expand.grid(R3 = 0:1, R2 = 0:1, R1 = 0:1, z = c(-1, 1),
    u = c(-1, 1), k = 1:3)
  N <- nrow(grid); id <- seq_len(N); calendar <- c(2, 5, 11)
  obs <- as.matrix(grid[c("R1", "R2", "R3")]) != 0
  A <- rowSums(obs) >= 2; U <- switch(universe, F = rep(TRUE, N),
    B = rowSums(obs) >= 1, A = A)
  z <- grid$z; k <- grid$k; g <- c(11, 23, Inf)[k]
  x <- grid$u + .25 * (grid$R1 - grid$R3) + .1 * k * z
  prob <- vapply(1:3, function(d) stats::plogis(-.6 + .3 * d + .2 * k + .35 * z), numeric(N))
  w <- apply(prob^obs * (1 - prob)^(!obs), 1L, prod) * (1 + .15 * grid$u)
  pu <- if (universe == "F") rep(1, N) else if (universe == "B")
    1 - apply(1 - prob, 1L, prod) else
      prob[, 1] * prob[, 2] + prob[, 1] * prob[, 3] + prob[, 2] * prob[, 3] -
        2 * apply(prob, 1L, prod)
  rho <- pu
  for (kk in 1:3) rho[k == kk] <- pu[k == kk] /
    (sum(w[k == kk] * U[k == kk]) / sum(w[k == kk]))
  Y <- vapply(1:3, function(d) 1 + d + k * d / 3 +
    (.3 * k * d + .09 * d^2) * x + .13 * d * x^2 + .07 * d * z * grid$u, numeric(N))
  Y[!obs] <- NA_real_
  X <- matrix(numeric(), N, 0)
  if (intercept) X <- cbind(`(Intercept)` = rep(1, N))
  if (slopes) X <- cbind(X, x = x)
  Q <- cbind(`(Intercept)` = 1, z = z)
  ctrl <- drp_ns("ipw_owner_control")(11L, 1e8, 2e8, 1e13, 1e8)
  long <- data.frame(id = rep(id, each = 3), time = rep(calendar, N),
    y = as.vector(t(Y)), g = rep(g, each = 3))
  panel <- drp_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  treatment_v <- grid$u + .3 * z
  V <- cbind(`(Intercept)` = 1, treatment_v = treatment_v)
  Xa <- if (intercept) cbind(X, treatment_v = treatment_v) else cbind(`(Intercept)` = 1, X, treatment_v = treatment_v)
  Xa[!A, ] <- NA_real_; Qa <- Q; Qa[!U, ] <- NA_real_
  frame <- drp_ns("construct_ipw_frame_owner")(panel, id, Xa, w, ctrl,
    frame_g = g, design_on_A_only = TRUE)
  reader <- drp_ns("construct_cp04_frame_reader")(frame, id, Qa, g, w, ctrl,
    U = U, rho = rho, universe = universe,
    justification = "known finite-law observation-subset transport ratios")
  defs <- list(list(s = 2, t = 11, g = 11, K = c(23, Inf)),
    list(s = 5, t = 11, g = 11, K = Inf), list(s = 2, t = 11, g = 23, K = Inf))
  views <- lapply(seq_along(defs), function(i) {
    a <- defs[[i]]
    drp_ns("construct_ipw_link_view")(frame, a$s, a$t, a$g, a$K, control = ctrl,
      link_id = paste0("selection-pooling-", i), functional_key = paste0("cmpcat:selection-pooling-", i))
  })
  recipe <- drp_ns("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  ref <- NULL; ref_X <- ref_Q <- ref_W <- ref_V <- NULL
  if (dated) {
    stopifnot(universe == "F", intercept, slopes)
    di <- rep(1:3, N)
    long$x <- rep(x, each = 3) + di * rep(z * grid$u / 7, each = 3)
    long$z <- rep(z, each = 3) * (1 + di / 20)
    long$treatment_v <- rep(treatment_v, each = 3) + di * rep(z / 30, each = 3)
    long$w <- rep(w, each = 3) * (1 + di * rep((grid$u + 1) / 10, each = 3))
    key <- match(long$time, calendar) * (N + 1) + long$id; wo <- order(key)
    route <- drp_ns("chained_reference_route")(route_kind, "F"); route$or_basis <- TRUE
    formulas <- list(x = list(formula = ~x + treatment_v))
    if (route_kind != "baseline") formulas$s <- list(formula = ~z)
    reference_vars <- c("x", "z", "treatment_v")
    if (extra_reference) {
      long$unused_a <- rep(x^2 + .2 * grid$R2, each = 3) + di / 13
      long$unused_b <- rep(z * grid$u + .3 * grid$R1, each = 3) + di * rep(grid$R3 / 17, each = 3)
      formulas$x$formula <- ~x + treatment_v + unused_a + unused_b
      reference_vars <- c(reference_vars, "unused_a", "unused_b")
    }
    src <- drp_ns("chained_reference_source")(formulas,
      data.frame(variable = reference_vars, kind = "dated", source = "data"), long,
      data.frame(id = id), FALSE, long$id, long$time, "time", id, calendar, 0,
      list(kind = "varying", name = "w", key = key[wo], value = long$w[wo], V = w,
        example = "finite-law dated reference weights"), "base_period", route)
    ref <- drp_ns("build_chained_reference_reader")(src, reader, views, length(views), ctrl)
    recipe <- drp_ns("or_pooled_outcome_recipe")(c("(Intercept)", "x"), TRUE, "x")
    ref_V <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, treatment_v = treatment_v + d * z / 30))
    ref_X <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, x = x + d * z * grid$u / 7))
    ref_Q <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, z = z * (1 + d / 20)))
    ref_W <- lapply(ref$recipes$weight_index, function(d) w * (1 + d * (grid$u + 1) / 10))
  }
  list(reader = reader, views = views, control = ctrl, recipe = recipe, reference = ref,
    X = X, V = V, Q = Q, Y = Y, obs = obs, g = g, w = w, U = U, A = A,
    rho = rho, calendar = calendar, ref_X = ref_X, ref_Q = ref_Q, ref_W = ref_W, ref_V = ref_V)
}

drp_reference <- function(a, fit, f = rep(1, nrow(a$Y)), nuisance = TRUE) {
  N <- length(f); nf <- if (nuisance) f else rep(1, N)
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  general <- isTRUE(fit$uses_selection_models)
  if (general) {
    cells <- fit$selection_backend$cells; blocks <- fit$selection_backend$block_specs
  } else {
    cells <- blocks <- list()
  }
  raw <- function(r, what) if (is.null(r)) a[[what]] else
    a[[switch(what, X = "ref_X", Q = "ref_Q", w = "ref_W")]][[r]]
  pred <- vector("list", length(cells))
  for (bl in blocks) {
    ids <- bl$cell_ids; xx <- yy <- ww <- list()
    for (j in seq_along(ids)) {
      cell <- cells[[ids[j]]]; rr <- which(a$g == st$cohorts[cell$cohort] & a$U)
      Q <- raw(cell$reference, "Q"); W <- raw(cell$reference, "w")
      E <- matrix(0, length(rr), length(ids)); E[, j] <- 1
      xx[[j]] <- cbind(E, Q[rr, 2L, drop = FALSE])
      yy[[j]] <- if (identical(cell$type, "date")) as.double(a$obs[rr, cell$calendar_index]) else
        as.double(a$obs[rr, cell$pair$s_calendar_index] & a$obs[rr, cell$pair$t_calendar_index])
      ww[[j]] <- nf[rr] * W[rr]
    }
    root <- suppressWarnings(stats::glm.fit(do.call(rbind, xx), unlist(yy), weights = unlist(ww),
      start = rep(0, length(ids) + 1L), family = stats::binomial(), control = drp_gc, intercept = FALSE))
    stopifnot(root$converged, all(is.finite(root$coefficients)))
    for (j in seq_along(ids)) {
      Q <- raw(cells[[ids[j]]]$reference, "Q")
      eta <- root$coefficients[j] + Q[, 2L] * root$coefficients[length(ids) + 1L]
      pred[[ids[j]]] <- stats::plogis(eta)
    }
  }
  ss <- lapply(seq_along(fit$support_specs), function(s) {
    sp <- fit$support_specs[[s]]; r <- sp$reference; W <- raw(r, "w"); X <- raw(r, "X")
    ik <- a$g == st$cohorts[sp$cohort]; iu <- which(ik & a$U)
    rate <- if (general) sum(nf[iu] * W[iu] * a$A[iu]) / sum(nf[iu] * W[iu]) else {
      io <- which(ik & a$obs[, sp$pair$s_calendar_index] & a$obs[, sp$pair$t_calendar_index])
      ia <- which(ik & a$A)
      sum(nf[io] * W[io]) / sum(nf[ia] * W[ia])
    }
    if (general) {
      pc <- fit$selection_backend$pair_cohorts[[fit$q_of_s[s]]]
      pp <- if (fit$selection_model == "independent_period_logit") pred[[pc$s_cell_id]] * pred[[pc$t_cell_id]] else pred[[pc$pair_cell_id]]
    } else pp <- rep(1, N)
    dy <- a$Y[, sp$pair$t_calendar_index] - a$Y[, sp$pair$s_calendar_index]
    ix <- which(ik & is.finite(dy)); q <- numeric(N)
    q[ix] <- if (general) W[ix] * rate / (a$rho[ix] * pmax(drp_ns("ipw_clip_eps"), pp[ix])) else
      W[ix] / max(drp_ns("ipw_clip_eps"), rate)
    list(X = X, dy = dy, ix = ix, q = q, L = sp$transitions)
  })
  G <- length(fit$regression_specs); icpt <- a$recipe$intercept
  separate <- fit$outcome_pooling == "separate"
  pieces <- lapply(seq_len(G), function(b) {
    keep <- if (separate) fit$outcome_groups$groups[[b]]$keep else fit$shared$keep
    xx <- yy <- ww <- list()
    for (s in fit$regression_specs[[b]]$support_ids) {
      h <- ss[[s]]; ix <- h$ix; X <- h$L * h$X[ix, a$recipe$slope_index[keep], drop = FALSE]
      if (icpt) {
        E <- matrix(0, length(ix), if (separate) 1L else G); E[, if (separate) 1L else b] <- 1
        X <- cbind(E, X)
      }
      xx[[length(xx) + 1L]] <- X; yy[[length(yy) + 1L]] <- h$dy[ix]
      ww[[length(ww) + 1L]] <- f[ix] * h$q[ix]
    }
    list(X = do.call(rbind, xx), y = unlist(yy), w = unlist(ww), keep = keep)
  })
  solve_wls <- function(p) unname(stats::lm.wfit(p$X, p$y, p$w, tol = 1e-10)$coefficients)
  if (separate) theta <- lapply(pieces, solve_wls) else {
    root <- solve_wls(list(X = do.call(rbind, lapply(pieces, `[[`, "X")),
      y = unlist(lapply(pieces, `[[`, "y")), w = unlist(lapply(pieces, `[[`, "w"))))
    theta <- lapply(seq_len(G), function(b) c(if (icpt) root[b], root[icpt * G + seq_along(pieces[[b]]$keep)]))
  }
  residual <- function(b, s) {
    h <- ss[[s]]; ix <- h$ix; D <- h$L * h$X[ix, a$recipe$slope_index[pieces[[b]]$keep], drop = FALSE]
    if (icpt) D <- cbind(1, D)
    list(ix = ix, r = h$dy[ix] - as.vector(D %*% theta[[b]]))
  }
  objective <- 0
  for (b in seq_len(G)) for (s in fit$regression_specs[[b]]$support_ids) {
    e <- residual(b, s); objective <- objective + sum(f[e$ix] * ss[[s]]$q[e$ix] * e$r^2) / N
  }
  mu <- den <- num <- numeric(length(fit$links)); design <- list()
  for (j in seq_along(fit$links)) {
    l <- fit$links[[j]]; s <- l$support_id; b <- l$regression; e <- residual(b, s)
    ix <- e$ix; qw <- f[ix] * ss[[s]]$q[ix]
    den[j] <- sum(qw) / N; num[j] <- sum(qw * e$r) / N; mu[j] <- num[j] / den[j]
    design[[j]] <- colSums(qw * ss[[s]]$L * ss[[s]]$X[ix, a$recipe$slope_index[pieces[[b]]$keep], drop = FALSE]) / N
  }
  # Refit each actual risk-set treatment root on the source-owned weights.
  gamma <- vector("list", length(fit$treatment_specs)); coords <- vector("list", length(gamma))
  dr <- ipw <- numeric(length(gamma))
  for (tid in seq_along(gamma)) {
    spec <- fit$treatment_specs[[tid]]; b <- spec$regression
    rr <- unlist(lapply(ss[spec$support_ids], `[[`, "ix"), use.names = FALSE)
    V <- if (is.null(spec$reference)) a$V else a$ref_V[[spec$reference]]
    dy <- q <- residual_y <- numeric(length(rr)); pos <- 0L
    y <- as.double(a$g[rr] == st$cohorts[spec$cohorts[1L]])
    for (sid in spec$support_ids) {
      h <- ss[[sid]]; n <- length(h$ix); ix <- seq_len(n) + pos; pos <- pos + n
      dy[ix] <- h$dy[h$ix]; q[ix] <- h$q[h$ix]
      residual_y[ix] <- residual(b, sid)$r
    }
    root <- suppressWarnings(stats::glm.fit(V[rr, , drop = FALSE], y, weights = f[rr] * q,
      start = rep(0, ncol(V)), family = stats::binomial(), control = drp_gc, intercept = FALSE))
    stopifnot(root$converged, all(is.finite(root$coefficients)))
    gamma[[tid]] <- unname(root$coefficients)
    pp <- stats::plogis(as.vector(V[rr, , drop = FALSE] %*% root$coefficients))
    odds <- pmin(pp, 1 - 1e-6) / (1 - pmin(pp, 1 - 1e-6))
    wp <- f[rr] * q * y; wm <- f[rr] * q * (1 - y) * odds
    co <- c(sum(wp * residual_y) / sum(wp), sum(wm * residual_y) / sum(wm), sum(wp) / N, sum(wm) / N)
    coords[[tid]] <- co; dr[tid] <- co[1L] - co[2L]
    ipw[tid] <- sum(wp * dy) / sum(wp) - sum(wm * dy) / sum(wm)
  }
  link_tids <- vapply(fit$links, `[[`, integer(1), "model")
  c(unlist(theta, use.names = FALSE), unlist(gamma, use.names = FALSE),
    dr[link_tids], ipw[link_tids], unlist(coords, use.names = FALSE))
}

drp_fit <- function(a, selection = "none", model = "pair_logit", pooling = "separate", mode = "separate") {
  fn <- if (selection == "none") "fit_dr_baseline" else "fit_dr_corollary1"
  args <- list(reader = a$reader, link_views = a$views, control = a$control,
    glm_control = drp_gc, rank_tol = 1e-10, outcome_recipe = a$recipe,
    treatment_columns = c("(Intercept)", "treatment_v"), reference = a$reference, outcome_pooling = mode)
  if (selection != "none") args <- c(args, list(selection_model = model, selection_pooling = pooling))
  do.call(drp_ns(fn), args)
}

drp_point <- function(fit) {
  theta <- lapply(seq_along(fit$regressions), function(b) c(
    if (fit$outcome_recipe$intercept) fit$regressions[[b]]$alpha,
    if (fit$outcome_pooling == "separate") fit$outcome_groups$groups[[b]]$beta else fit$shared$beta))
  c(unlist(theta, use.names = FALSE), unlist(lapply(fit$treatments, function(t) unname(t$coefficients)), use.names = FALSE),
    vapply(fit$links, `[[`, numeric(1), "mu"), vapply(fit$links, `[[`, numeric(1), "mu_ipw"),
    unlist(lapply(fit$treatments, function(t) unname(t$means[c("T_dr", "M_dr", "dp", "dm")])), use.names = FALSE))
}

drp_details <- function(fit) {
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  answer <- matrix(NA_real_, N, length(drp_point(fit)))
  drp_ns("replay_dr_if")(fit, 11L, length(fit$links), function(o) {
    d <- o$details
    theta <- lapply(d$blocks, function(b) cbind(if (fit$outcome_recipe$intercept) b$alpha, b$beta))
    gamma <- lapply(d$models[d$treatment_ids], `[[`, "gamma")
    coords <- lapply(d$models[d$treatment_ids], `[[`, "coordinates_full")
    answer[o$rows, ] <<- do.call(cbind, c(theta, gamma, list(o$phi, d$psi_ipw), coords))
    testthat::expect_equal(o$phi, d$psi_dir + d$psi_or + d$psi_ps + d$psi_sel, tolerance = 1e-10)
  }, details = TRUE)
  answer
}

testthat::test_that("DR points use block-local control slopes and unchanged treatment roots", {
  for (selection in c("none", "conditional_mtar")) for (spec in list(c(TRUE, TRUE), c(FALSE, TRUE), c(TRUE, FALSE))) {
    a <- drp_case(intercept = spec[1L], slopes = spec[2L])
    sep <- drp_fit(a, selection); pooled <- drp_fit(a, selection, mode = "pooled")
    for (fit in list(sep, pooled)) testthat::expect_equal(drp_point(fit), drp_reference(a, fit), tolerance = 2e-8)
    testthat::expect_null(sep$shared); testthat::expect_null(sep$complete_support)
    testthat::expect_equal(sep$accounting$n_shared_slopes, 0L)
    testthat::expect_equal(sep$accounting$n_regression_fits, length(sep$regressions))
    for (j in seq_along(sep$treatments)) {
      testthat::expect_equal(sep$treatments[[j]]$coefficients, pooled$treatments[[j]]$coefficients, tolerance = 1e-12)
      testthat::expect_equal(unname(sep$treatments[[j]]$means["mu_dr"]), unname(sep$treatments[[j]]$means["mu_dr_decomp"]), tolerance = 1e-10)
    }
    if (!spec[2L]) {
      testthat::expect_equal(vapply(sep$links, `[[`, numeric(1), "mu"), vapply(sep$links, `[[`, numeric(1), "mu_ipw"), tolerance = 1e-10)
      drp_ns("replay_dr_if")(sep, 11L, 3L, function(o)
        testthat::expect_equal(o$phi, o$details$psi_ipw, tolerance = 1e-10), details = TRUE)
    }
    fn <- if (selection == "none") "fit_dr_baseline" else "fit_dr_corollary1"
    default <- drp_ns(fn)(a$reader, a$views, a$control, drp_gc, 1e-10, a$recipe,
      treatment_columns = c("(Intercept)", "treatment_v"))
    testthat::expect_identical(default$outcome_pooling, "pooled")
    testthat::expect_identical(default$shared, pooled$shared)
    testthat::expect_identical(default$treatments, pooled$treatments)
  }
})

testthat::test_that("separate DR full-N IF refits all OR, observation and treatment roots", {
  specs <- list(list(selection = "none", model = "pair_logit", pooling = "separate", u = "F"),
    list(selection = "conditional_mtar", model = "pair_logit", pooling = "separate", u = "F"),
    list(selection = "conditional_mtar", model = "pair_logit", pooling = "within_cohort", u = "B"),
    list(selection = "conditional_mtar", model = "pair_logit", pooling = "across_cohorts", u = "A"),
    list(selection = "conditional_mtar", model = "independent_period_logit", pooling = "across_cohorts", u = "F"))
  for (spec in specs) {
    a <- drp_case(spec$u); fit <- drp_fit(a, spec$selection, spec$model, spec$pooling)
    actual <- drp_details(fit); N <- nrow(actual); h <- 1e-6
    numerical <- t(vapply(seq_len(N), function(i) {
      delta <- -rep(1, N); delta[i] <- N - 1
      (drp_reference(a, fit, 1 + h * delta) - drp_reference(a, fit, 1 - h * delta)) / (2 * h)
    }, numeric(ncol(actual))))
    testthat::expect_equal(unname(actual), numerical, tolerance = 1e-5)
    testthat::expect_lt(max(abs(colMeans(actual))), 1e-8)
    if (spec$selection != "none" && spec$u != "A")
      testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
    testthat::expect_lte(fit$accounting$owned_entries, fit$accounting$retained_bound)
    detailed <- drp_ns("evaluate_dr_if_chunk")(fit, c(3, 1, 9), c(3, 1), details = TRUE)
    ordinary <- drp_ns("evaluate_dr_if_chunk")(fit, c(3, 1, 9), c(3, 1))
    testthat::expect_equal(detailed$phi, ordinary$phi, tolerance = 1e-10)
    core <- drp_ns("evaluate_dr_if_chunk")(fit, c(3, 1, 9), c(3, 1), scale = "core")
    lambda <- sum(a$A) / N
    testthat::expect_equal(core$phi, lambda * ordinary$phi, tolerance = 1e-10)
    testthat::expect_error(drp_ns("dr_if_budget")(fit, 1e20, 1:3, TRUE))
  }
})

testthat::test_that("separate DR no-intercept and dated recipes keep complete coefficient IFs", {
  specs <- list(list(selection = "none", no_intercept = TRUE, dated = FALSE),
    list(selection = "conditional_mtar", no_intercept = TRUE, dated = FALSE),
    list(selection = "none", no_intercept = FALSE, dated = TRUE),
    list(selection = "conditional_mtar", no_intercept = FALSE, dated = TRUE))
  for (spec in specs) {
    a <- drp_case(intercept = !spec$no_intercept, dated = spec$dated,
      route_kind = if (spec$selection == "none") "baseline" else "corollary1")
    fit <- drp_fit(a, spec$selection, pooling = "within_cohort")
    actual <- drp_details(fit); h <- 1e-6; N <- nrow(actual)
    testthat::expect_equal(drp_point(fit), drp_reference(a, fit), tolerance = 2e-8)
    for (i in c(1L, 10L, 24L, 82L)) {
      delta <- -rep(1, N); delta[i] <- N - 1
      numerical <- (drp_reference(a, fit, 1 + h * delta) - drp_reference(a, fit, 1 - h * delta)) / (2 * h)
      testthat::expect_equal(unname(actual[i, ]), numerical, tolerance = 1e-5)
    }
  }
})

testthat::test_that("repeated long-ID DR links reuse fits and admit requested subset metadata", {
  for (selection in c("none", "conditional_mtar")) {
    a <- drp_case(); original <- drp_fit(a, selection)
    frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
    make_views <- function(long_ids) c(a$views, lapply(4:5, function(j) {
      id <- paste0("repeat-", j, if (long_ids) paste0("-", strrep("q", 4096L)) else "")
      drp_ns("construct_ipw_link_view")(frame, 2, 11, 11, c(23, Inf),
        control = a$control, link_id = id, functional_key = paste0("cmpcat:", id))
    }))
    short_case <- long_case <- a
    short_case$views <- make_views(FALSE); long_case$views <- make_views(TRUE)
    short <- drp_fit(short_case, selection); fit <- drp_fit(long_case, selection)
    testthat::expect_equal(length(fit$outcome_groups$groups), length(original$outcome_groups$groups))
    testthat::expect_equal(fit$accounting$n_regression_fits, original$accounting$n_regression_fits)
    testthat::expect_equal(fit$accounting$n_treatment_fits, original$accounting$n_treatment_fits)
    testthat::expect_equal(fit$accounting$n_selection_fits, original$accounting$n_selection_fits)
    for (j in 4:5) {
      testthat::expect_equal(fit$links[[j]]$regression, fit$links[[1L]]$regression)
      testthat::expect_equal(fit$links[[j]]$model, fit$links[[1L]]$model)
    }
    testthat::expect_equal(drp_point(fit), drp_reference(long_case, fit), tolerance = 2e-8)
    all_rows <- 1:11; requested_rows <- c(10L, 3L, 1L); requested <- c(5L, 2L, 4L)
    whole <- drp_ns("evaluate_dr_if_chunk")(fit, all_rows, 1:5, details = TRUE)
    before <- drp_ns("evaluate_dr_if_chunk")(original, all_rows, 1:3, details = TRUE)
    subset <- drp_ns("evaluate_dr_if_chunk")(fit, requested_rows, requested, details = TRUE)
    plain <- drp_ns("evaluate_dr_if_chunk")(fit, requested_rows, requested)
    testthat::expect_equal(whole$phi[, 1:3, drop = FALSE], before$phi, tolerance = 1e-10)
    testthat::expect_equal(subset$phi, whole$phi[match(requested_rows, all_rows), requested, drop = FALSE], tolerance = 1e-10)
    testthat::expect_equal(plain$phi, subset$phi, tolerance = 1e-10)
    long_bound <- drp_ns("dr_if_envelope")(fit, length(requested_rows), requested, TRUE)
    short_bound <- drp_ns("dr_if_envelope")(short, length(requested_rows), requested, TRUE)
    name_bytes <- function(x) sum(vapply(x$links[requested], function(l) nchar(l$link_id, type = "bytes"), integer(1)))
    reserved_extra <- 128 * (name_bytes(fit) - name_bytes(short))
    testthat::expect_gt(reserved_extra, 0)
    testthat::expect_gte(long_bound$output - short_bound$output, reserved_extra)
    testthat::expect_gte(long_bound$scratch - short_bound$scratch, reserved_extra)
  }
  # The sealed dated reader materializes all source x columns before selecting
  # the unchanged outcome/treatment recipes. Unused predictors cost scratch.
  narrow_case <- drp_case(dated = TRUE, route_kind = "baseline")
  wide_case <- drp_case(dated = TRUE, route_kind = "baseline", extra_reference = TRUE)
  narrow <- drp_fit(narrow_case); wide <- drp_fit(wide_case)
  testthat::expect_identical(wide$outcome_recipe, narrow$outcome_recipe)
  testthat::expect_identical(wide$treatment_recipe$columns, narrow$treatment_recipe$columns)
  testthat::expect_gt(length(wide_case$reference$x_columns), length(narrow_case$reference$x_columns))
  testthat::expect_equal(drp_point(wide), drp_point(narrow), tolerance = 2e-8)
  wr <- drp_ns("evaluate_dr_if_chunk")(wide, c(10L, 3L, 1L), c(3L, 1L), details = TRUE)
  nr <- drp_ns("evaluate_dr_if_chunk")(narrow, c(10L, 3L, 1L), c(3L, 1L), details = TRUE)
  testthat::expect_equal(wr$phi, nr$phi, tolerance = 1e-10)
  wb <- drp_ns("dr_if_envelope")(wide, 3L, c(3L, 1L), TRUE)
  nb <- drp_ns("dr_if_envelope")(narrow, 3L, c(3L, 1L), TRUE)
  testthat::expect_gt(wb$scratch, nb$scratch)
  testthat::expect_gt(wb$work, nb$work)
})
