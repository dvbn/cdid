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

# V32 changes the treatment fitting reservation, not its fitted equations.
testthat::test_that("conditional DR fitting counts local channel and cohort-row products rather than global D squared", {
  st <- list(cohort_start = c(1, 101, 221), cohort_end = c(100, 220, 370))
  specs <- list(list(n_observers = 50L, cohorts = c(1L, 2L), selection_ids = c(1L, 2L)),
    list(n_observers = 200L, cohorts = c(2L, 3L), selection_ids = c(2L, 3L)))
  map <- list(D = 7L, pair_channels = list(1:3, c(1L, 4L, 5L, 6L), c(2L, 7L)))
  charge <- drp_ns("dr_selection_treatment_fit_scratch")(st, specs, map, 2, 1, 2, 11, name_bytes = 17)
  # Second actual model: n=200,u=4+2=6,h=4,ell=max(120*4,150*2)=480.
  # Hand-counted buffers: input2200,chunk14432,GLM16064,derivative18272,names1088.
  testthat::expect_equal(charge[c("n", "channels", "pair_channels", "cohort_derivative_entries")],
    list(n = 200, channels = 6, pair_channels = 4, cohort_derivative_entries = 480))
  testthat::expect_equal(charge[c("input", "chunk", "solver", "derivative", "names", "scratch")],
    list(input = 2200, chunk = 14432, solver = 16064, derivative = 18272, names = 1088, scratch = 44184))
  many <- map; many$D <- 1000000L
  testthat::expect_identical(drp_ns("dr_selection_treatment_fit_scratch")(st, specs, many, 2, 1, 2, 11,
    name_bytes = 17), charge)
  wide <- map; wide$pair_channels[[2L]] <- seq_len(200L)
  wide_charge <- drp_ns("dr_selection_treatment_fit_scratch")(st, specs, wide, 2, 1, 2, 11, name_bytes = 17)
  testthat::expect_gt(wide_charge$scratch, charge$scratch)
  testthat::expect_equal(wide_charge$cohort_derivative_entries, 120 * 200)
  ref <- list(x_columns = paste0("x", 1:9), role_slot = c(s = 1L), designs = list(list(columns = paste0("s", 1:7))))
  dated <- drp_ns("dr_selection_treatment_fit_scratch")(st, specs, map, 2, 1, 2, 11, ref, 17)
  testthat::expect_equal(dated$reference_x_columns, 9)
  testthat::expect_equal(dated$reference_s_columns, 7)
  testthat::expect_equal(dated$scratch - charge$scratch, 32 * 11 * (9 + 7))
  # Exact finite capacity is retained: local width cannot escape the budget.
  own <- 123; cap <- own + wide_charge$scratch
  ctrl <- drp_ns("ipw_owner_control")(11L, cap, cap, 1e13, cap)
  testthat::expect_silent(drp_ns("cp04_budget")(ctrl, own, wide_charge$scratch, 0, 0, "V32 local-wide witness"))
  below <- ctrl; below$max_matrix_entries <- cap - 1
  testthat::expect_error(drp_ns("cp04_budget")(below, own, wide_charge$scratch, 0, 0, "V32 local-wide witness"),
    class = "cdid_ipw_owner_bound_error")
})

# Observe named live fit buffers and actual matrix/slice/result arguments.
# Borrowed sealed inputs/owners are already in the caller's retained prefix;
# only scalar maxima survive these test-local checkpoints.
drp32_fit_payload <- function(a, model, pooling, mode) {
  original <- drp_ns("dr_selection_treatment_fit"); records <- list()
  wrapped <- function(...) {
    args <- list(...); names(args) <- names(formals(original))[seq_along(args)]
    st <- args$reader$state
    d <- if (is.null(args$reference)) ncol(st$Z) else
      length(args$reference$designs[[args$reference$role_slot["s"]]]$columns)
    name_bytes <- sum(nchar(c(args$recipe$columns, args$trt_recipe$columns,
      if (is.null(args$reference)) colnames(st$Z) else
        args$reference$designs[[args$reference$role_slot["s"]]]$columns,
      if (!is.null(args$reference)) args$reference$x_columns), type = "bytes"))
    inv <- drp_ns("dr_selection_treatment_fit_scratch")(st, list(args$spec), args$up_map,
      args$trt_recipe$p, args$recipe$p, d, min(st$N, args$ctrl$max_chunk_rows), args$reference, name_bytes)
    peak <- 0; matrices <- crosses <- derivatives <- 0L
    remember <- function(env, extra) {
      owned <- as.list(env, all.names = TRUE)
      owned[intersect(names(owned), c(names(formals(original)), "st", "frame"))] <- NULL
      peak <<- max(peak, unname(drp_ns("cp04_payload")(c(owned, extra))["entries"]))
    }
    e <- new.env(parent = environment(original)); run <- original; environment(run) <- e
    e$matrix <- function(...) {
      value <- base::matrix(...); matrices <<- matrices + 1L
      remember(parent.frame(), list(new_matrix = value)); value
    }
    e$crossprod <- function(x, y = NULL) {
      value <- base::crossprod(x, y); crosses <<- crosses + 1L
      remember(parent.frame(), list(slice_x = x, slice_y = y, result = value)); value
    }
    e$selection_prediction_operation <- function(...) {
      op <- drp_ns("selection_prediction_operation")(...)
      derivative <- op$deriv
      op$deriv <- function(...) {
        value <- derivative(...); derivatives <<- derivatives + 1L
        remember(parent.frame(), list(returned_derivative = value)); value
      }
      op
    }
    value <- do.call(run, args)
    testthat::expect_lte(peak, inv$scratch)
    testthat::expect_gt(matrices, 0L); testthat::expect_gt(crosses, 0L); testthat::expect_gt(derivatives, 0L)
    records[[length(records) + 1L]] <<- list(peak = peak, inventory = inv)
    value
  }
  testthat::local_mocked_bindings(dr_selection_treatment_fit = wrapped, .package = "cdid")
  fit <- drp_fit(a, "conditional_mtar", model, pooling, mode)
  list(fit = fit, records = records)
}

testthat::test_that("conditional DR local scratch dominates actual pooled and separate pair/date fit buffers", {
  cases <- list(list(model = "pair_logit", pooling = "separate", mode = "separate", dated = FALSE),
    list(model = "pair_logit", pooling = "within_cohort", mode = "pooled", dated = FALSE),
    list(model = "independent_period_logit", pooling = "across_cohorts", mode = "separate", dated = FALSE),
    list(model = "independent_period_logit", pooling = "within_cohort", mode = "pooled", dated = TRUE))
  for (spec in cases) {
    a <- drp_case(dated = spec$dated, extra_reference = spec$dated)
    observed <- drp32_fit_payload(a, spec$model, spec$pooling, spec$mode); fit <- observed$fit
    inventory <- fit$accounting$treatment_fit_inventory
    testthat::expect_length(observed$records, length(fit$treatment_specs))
    testthat::expect_equal(inventory$scratch,
      max(vapply(observed$records, function(r) r$inventory$scratch, numeric(1))))
    for (r in observed$records) testthat::expect_lte(r$peak, inventory$scratch)
    testthat::expect_lte(inventory$scratch, fit$accounting$fit_scratch_bound)
    testthat::expect_lte(fit$accounting$owned_entries, fit$accounting$retained_bound)
    if (spec$dated) {
      testthat::expect_equal(inventory$reference_x_columns, length(a$reference$x_columns))
      testthat::expect_gt(inventory$reference_x_columns, length(a$recipe$columns))
    }
    # Same independent raw-row WLS/selection/treatment refits as the accepted
    # scientific tests; every original F row remains in ordinary/detailed IF.
    testthat::expect_equal(drp_point(fit), drp_reference(a, fit), tolerance = 2e-8)
    N <- nrow(a$Y); ordinary <- detailed <- matrix(0, N, length(fit$links))
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- seq.int(first, min(N, first + 10L))
      ordinary[rows, ] <- drp_ns("evaluate_dr_if_chunk")(fit, rows, seq_along(fit$links))$phi
      detailed[rows, ] <- drp_ns("evaluate_dr_if_chunk")(fit, rows, seq_along(fit$links), details = TRUE)$phi
    }
    testthat::expect_equal(ordinary, detailed, tolerance = 1e-10)
    testthat::expect_equal(crossprod(ordinary), crossprod(detailed), tolerance = 1e-10)
    # Exact actual fitting-stage boundary, reusing the same already fitted
    # selection nuisance so a refusal is measured before any treatment fit.
    if (spec$model == "pair_logit" && spec$pooling == "separate") {
      cap <- fit$accounting$retained_bound + fit$accounting$fit_scratch_bound
      refit <- function(limit) {
        ctrl <- a$control; ctrl$max_matrix_entries <- limit
        drp_ns("fit_dr_corollary1")(a$reader, a$views, ctrl, drp_gc, 1e-10, a$recipe,
          treatment_columns = c("(Intercept)", "treatment_v"), sel_fit = fit$selection_backend,
          outcome_pooling = spec$mode)
      }
      at <- refit(cap)
      testthat::expect_equal(drp_point(at), drp_point(fit), tolerance = 1e-12)
      refuse_before_fit <- function() {
        called <- 0L
        testthat::local_mocked_bindings(dr_selection_treatment_fit = function(...) {
          called <<- called + 1L; stop("Treatment fit ran before the required refusal.")
        }, .package = "cdid")
        testthat::expect_error(refit(cap - 1), class = "cdid_ipw_owner_bound_error")
        testthat::expect_identical(called, 0L)
      }
      refuse_before_fit()
    }
  }
})

# V33 keeps the old envelope as a test-local reference; both evaluators are
# the identical installed function, with only ordinary budgeting dispatched.
drp33_old <- local({
  e <- new.env(parent = asNamespace("cdid"))
  for (name in c("dr_corollary1_if_envelope", "dr_corollary1_if_budget", "evaluate_dr_corollary1_if_chunk")) {
    f <- drp_ns(name); environment(f) <- e; assign(name, f, e)
  }
  e$dr_corollary1_ordinary_if_envelope <- function(backend, B, links, pre, selection_pre, metadata) pre
  e
})
drp33_inventory <- function(fit, B, links) {
  old <- drp33_old$dr_corollary1_if_envelope(fit, B, links, FALSE)
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  metadata <- 16384 + 1024 * (length(old$sids) + length(old$uids) + length(fit$regression_specs) +
    length(old$tids) + length(links) + length(st$cohorts) + length(fit$selection_backend$blocks)) +
    B * frame$id_byte_metadata$max + 8 * sum(nchar(vapply(fit$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  selection <- drp_ns("selection_if_budget")(fit$selection_backend, B, unique(fit$q_of_s[old$sids]), details = FALSE)
  drp_ns("dr_corollary1_ordinary_inventory")(fit, B, links, old, selection, metadata)
}

# Source-owned checkpoints: all named current arrays, actual arithmetic
# arguments/results, global pooled residual cache, returned upstream records.
# Sealed borrowed fits/readers are excluded; only numeric maxima are retained.
drp33_payload <- function(fit, rows, links) {
  e <- new.env(parent = asNamespace("cdid")); peak <- 0; solves <- 0L; labels <- character()
  outer <- NULL
  own <- function(env, exclude) {
    x <- as.list(env, all.names = TRUE); x[intersect(names(x), exclude)] <- NULL; x
  }
  borrowed <- c("backend", "reader", "st", "frame", "sel_fit", "up_map", "ref", "spec", "f", "reg")
  remember <- function(extra = list()) {
    if (!is.null(outer)) peak <<- max(peak, unname(drp_ns("cp04_payload")(c(own(outer, borrowed), extra))["entries"]))
  }
  evaluate <- drp_ns("evaluate_dr_corollary1_if_chunk"); environment(evaluate) <- e
  e$matrix <- function(...) {
    if (is.null(outer)) outer <<- parent.frame()
    value <- base::matrix(...); remember(list(new_matrix = value)); value
  }
  e$cp04_solve <- function(R, b) {
    value <- drp_ns("cp04_solve")(R, b); solves <<- solves + 1L
    remember(list(solve_rhs = b, solve_result = value)); value
  }
  e$cp04_finite <- function(x, label, positive = FALSE) {
    if (is.null(outer)) outer <<- parent.frame()
    labels <<- c(labels, label); remember(list(finite_input = x))
    drp_ns("cp04_finite")(x, label, positive)
  }
  core <- drp_ns("or_pooled_if_core"); ce <- new.env(parent = environment(core)); environment(core) <- ce
  core_own <- function(env) own(env, c("model", "N", "w", "hits", "up", "pre", "links", "details", "beta_only", "h", "blk"))
  ce[["%*%"]] <- function(x, y) {
    value <- base::`%*%`(x, y)
    remember(list(core_live = core_own(parent.frame()), product_x = x, product_y = y, product_result = value)); value
  }
  for (name in c("or_pooled_wls_gram_solve", "or_wls_gram_solve")) local({
    original <- drp_ns(name)
    assign(name, function(...) {
      value <- original(...); args <- list(...)
      remember(list(core_live = core_own(parent.frame()), gram_rhs = args[[2L]], gram_result = value)); value
    }, ce)
  })
  e$or_pooled_if_core <- function(...) {
    if (is.null(outer)) outer <<- parent.frame()
    value <- core(...); remember(list(returned_beta_core = value)); value
  }
  out <- evaluate(fit, rows, links, "frame", FALSE)
  pre <- drp_ns("dr_corollary1_if_budget")(fit, length(rows), links, FALSE)
  testthat::expect_lte(peak, pre$scratch)
  testthat::expect_lte(unname(drp_ns("cp04_payload")(out)["entries"]), pre$output)
  testthat::expect_gt(solves, 0L)
  testthat::expect_true("complete Corollary 1 DR frame IF" %in% labels)
  list(out = out, peak = peak, pre = pre)
}

testthat::test_that("ordinary pooled conditional DR charges true outputs and global/common/local live arrays", {
  cases <- list(list(u = "F", model = "pair_logit", pooling = "separate", dated = FALSE, p0 = FALSE),
    list(u = "B", model = "pair_logit", pooling = "within_cohort", dated = FALSE, p0 = FALSE),
    list(u = "A", model = "pair_logit", pooling = "across_cohorts", dated = FALSE, p0 = FALSE),
    list(u = "F", model = "independent_period_logit", pooling = "across_cohorts", dated = TRUE, p0 = FALSE),
    list(u = "F", model = "pair_logit", pooling = "separate", dated = FALSE, p0 = TRUE))
  for (spec in cases) {
    a <- drp_case(spec$u, slopes = !spec$p0, dated = spec$dated, extra_reference = spec$dated)
    if (spec$u == "F" && !spec$dated && !spec$p0) {
      frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
      a$views <- lapply(seq_along(a$views), function(j) {
        v <- a$views[[j]]; id <- paste0("unicode-", j, strrep("\u03bb\u6f22", 240L))
        drp_ns("construct_ipw_link_view")(frame, v$pair$s, v$pair$t, v$cohorts$g, v$cohorts$K,
          control = a$control, link_id = id, functional_key = paste0("cmpcat:", id))
      })
    }
    fit <- drp_fit(a, "conditional_mtar", spec$model, spec$pooling, "pooled")
    links <- rev(seq_along(fit$links)); N <- nrow(a$Y)
    pre <- drp_ns("dr_corollary1_if_envelope")(fit, 11L, links, FALSE)
    old <- drp33_old$dr_corollary1_if_envelope(fit, 11L, links, FALSE)
    testthat::expect_lt(pre$output, old$output); testthat::expect_lte(pre$scratch, old$scratch)
    testthat::expect_identical(pre$work, old$work)
    testthat::expect_identical(drp_ns("dr_corollary1_if_envelope")(fit, 11L, links, TRUE),
      drp33_old$dr_corollary1_if_envelope(fit, 11L, links, TRUE))
    inv <- drp33_inventory(fit, 11L, links)
    W <- sum(vapply(pre$active_blocks, function(m) {
      b <- fit$selection_backend$blocks[[m]]
      (if (b$has_intercept) b$B_m else 0) + b$d_s + length(fit$up_map$theta_coords[[m]])
    }, numeric(1)))
    testthat::expect_equal(inv$selection_columns, W)
    testthat::expect_equal(inv$upstream, 11 * (length(pre$uids) + W))
    testthat::expect_equal(inv$max_model_channels,
      max(vapply(fit$treatments[pre$tids], function(f) length(f$upstream_ids), integer(1))))
    testthat::expect_equal(inv$scratch, pre$scratch); testthat::expect_equal(inv$output, pre$output)
    if (!spec$p0) {
      subset <- drp_ns("dr_corollary1_if_envelope")(fit, 11L, links[1L], FALSE)
      subset_inv <- drp33_inventory(fit, 11L, links[1L])
      # This fixture has control incidences2+1+1; the requested local block
      # cannot replace the global pooled coefficient's full source set.
      testthat::expect_identical(subset$beta_gids, seq_along(fit$regression_specs))
      testthat::expect_equal(subset_inv$beta_support_incidence, 4)
      testthat::expect_gt(subset_inv$beta_support_incidence,
        sum(vapply(fit$regression_specs[subset$gids], function(b) length(b$support_ids), integer(1))))
    } else testthat::expect_equal(inv$core_peak, 0)
    if (spec$dated) testthat::expect_equal(inv$reference_x_columns, length(a$reference$x_columns))
    testthat::expect_equal(drp_point(fit), drp_reference(a, fit), tolerance = 2e-8)
    actual <- reference <- detailed <- matrix(0, N, length(links))
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(first, min(N, first + 10L)))
      observed <- drp33_payload(fit, rows, links)$out
      frozen <- drp33_old$evaluate_dr_corollary1_if_chunk(fit, rows, links, "frame", FALSE)
      rich <- drp_ns("evaluate_dr_corollary1_if_chunk")(fit, rows, links, "frame", TRUE)
      actual[rows, ] <- observed$phi; reference[rows, ] <- frozen$phi; detailed[rows, ] <- rich$phi
    }
    testthat::expect_equal(actual, reference, tolerance = 1e-10)
    testthat::expect_equal(actual, detailed, tolerance = 1e-10)
    testthat::expect_equal(crossprod(actual), crossprod(detailed), tolerance = 1e-10)
    point_offset <- length(drp_point(fit)) - 2 * length(links) - 4 * length(fit$treatments)
    for (i in c(1L, 10L, N - 1L, N)) {
      delta <- rep(-1, N); delta[i] <- N - 1
      derivative <- (drp_reference(a, fit, 1 + 1e-6 * delta) - drp_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6
      testthat::expect_equal(actual[i, ], derivative[point_offset + links], tolerance = 1e-5)
    }
    # The actual output slope depends on requested phi/row/id width, not
    # on global nuisance dimensions or the number of treatment models.
    small <- drp_ns("dr_corollary1_if_envelope")(fit, 1L, links, FALSE)
    frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
    testthat::expect_equal(pre$output - small$output,
      10 * (length(links) + 4 + frame$id_byte_metadata$max))
  }
})

testthat::test_that("ordinary pooled conditional DR keeps adjacent caps and all four-coordinate numeric guards", {
  a <- drp_case(); fit <- drp_fit(a, "conditional_mtar", mode = "pooled"); links <- 1:3
  pre <- drp_ns("dr_corollary1_if_budget")(fit, 11L, links, FALSE)
  owned <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  reseal <- function(ctrl = fit$control, treatments = fit$treatments) {
    b <- fit; b$ref <- NULL; b$control <- ctrl; b$treatments <- treatments
    drp_ns("ipw_owner_seal")(b, "dr_corollary1_backend")
  }
  ctrl <- fit$control; ctrl$max_matrix_entries <- owned + pre$scratch
  at <- reseal(ctrl)
  testthat::expect_equal(drp_ns("evaluate_dr_corollary1_if_chunk")(at, 1:11, links)$phi,
    drp33_old$evaluate_dr_corollary1_if_chunk(fit, 1:11, links)$phi, tolerance = 1e-12)
  ctrl$max_matrix_entries <- ctrl$max_matrix_entries - 1
  testthat::expect_error(drp_ns("evaluate_dr_corollary1_if_chunk")(reseal(ctrl), 1:11, links),
    class = "cdid_ipw_owner_bound_error")
  ctrl <- fit$control; ctrl$max_output_entries <- pre$output
  testthat::expect_silent(drp_ns("evaluate_dr_corollary1_if_chunk")(reseal(ctrl), 1:11, links))
  ctrl$max_output_entries <- pre$output - 1
  testthat::expect_error(drp_ns("evaluate_dr_corollary1_if_chunk")(reseal(ctrl), 1:11, links),
    class = "cdid_ipw_owner_bound_error")
  for (kind in c("logit", "dp", "dm")) {
    tr <- fit$treatments
    if (kind == "logit") tr[[1L]]$coefficients[] <- Inf else
      tr[[1L]]$partials$upstream[if (kind == "dp") 3L else 4L, 1L] <- Inf
    b <- reseal(treatments = tr)
    actual <- tryCatch(drp_ns("evaluate_dr_corollary1_if_chunk")(b, 1:11, links), error = identity)
    reference <- tryCatch(drp33_old$evaluate_dr_corollary1_if_chunk(b, 1:11, links), error = identity)
    testthat::expect_s3_class(actual, "cdid_cp04_numeric_error")
    testthat::expect_identical(class(actual), class(reference))
    testthat::expect_identical(conditionMessage(actual), conditionMessage(reference))
  }
})
