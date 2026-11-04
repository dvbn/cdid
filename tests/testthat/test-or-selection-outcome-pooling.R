# Independent raw-data logit/WLS refits; no stored outputs or fixture edits.
osop_ns <- function(name) getFromNamespace(name, "cdid")
osop_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)

osop_case <- function(universe = "F", intercept = TRUE, slopes = TRUE, dated = FALSE) {
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
  ctrl <- osop_ns("ipw_owner_control")(11L, 1e8, 2e8, 1e13, 1e8)
  long <- data.frame(id = rep(id, each = 3), time = rep(calendar, N),
    y = as.vector(t(Y)), g = rep(g, each = 3))
  panel <- osop_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  Xa <- X; Xa[!A, ] <- NA_real_; Qa <- Q; Qa[!U, ] <- NA_real_
  frame <- osop_ns("construct_ipw_frame_owner")(panel, id, Xa, w, ctrl,
    frame_g = g, design_on_A_only = TRUE)
  reader <- osop_ns("construct_cp04_frame_reader")(frame, id, Qa, g, w, ctrl,
    U = U, rho = rho, universe = universe,
    justification = "known finite-law observation-subset transport ratios")
  defs <- list(list(s = 2, t = 11, g = 11, K = c(23, Inf)),
    list(s = 5, t = 11, g = 11, K = Inf), list(s = 2, t = 11, g = 23, K = Inf))
  views <- lapply(seq_along(defs), function(i) {
    a <- defs[[i]]
    osop_ns("construct_ipw_link_view")(frame, a$s, a$t, a$g, a$K, control = ctrl,
      link_id = paste0("selection-pooling-", i), functional_key = paste0("cmpcat:selection-pooling-", i))
  })
  recipe <- osop_ns("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  ref <- NULL; ref_X <- ref_Q <- ref_W <- NULL
  if (dated) {
    stopifnot(universe == "F", intercept, slopes)
    di <- rep(1:3, N)
    long$x <- rep(x, each = 3) + di * rep(z * grid$u / 7, each = 3)
    long$z <- rep(z, each = 3) * (1 + di / 20)
    long$w <- rep(w, each = 3) * (1 + di * rep((grid$u + 1) / 10, each = 3))
    key <- match(long$time, calendar) * (N + 1) + long$id; wo <- order(key)
    route <- osop_ns("chained_reference_route")("corollary1", "F"); route$or_basis <- TRUE
    src <- osop_ns("chained_reference_source")(list(x = list(formula = ~x), s = list(formula = ~z)),
      data.frame(variable = c("x", "z"), kind = "dated", source = "data"), long,
      data.frame(id = id), FALSE, long$id, long$time, "time", id, calendar, 0,
      list(kind = "varying", name = "w", key = key[wo], value = long$w[wo], V = w,
        example = "finite-law dated reference weights"), "base_period", route)
    ref <- osop_ns("build_chained_reference_reader")(src, reader, views, length(views), ctrl)
    recipe <- ref$outcome_recipe
    ref_X <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, x = x + d * z * grid$u / 7))
    ref_Q <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, z = z * (1 + d / 20)))
    ref_W <- lapply(ref$recipes$weight_index, function(d) w * (1 + d * (grid$u + 1) / 10))
  }
  list(reader = reader, views = views, control = ctrl, recipe = recipe, reference = ref,
    X = X, Q = Q, Y = Y, obs = obs, g = g, w = w, U = U, A = A,
    rho = rho, calendar = calendar, ref_X = ref_X, ref_Q = ref_Q, ref_W = ref_W)
}

osop_fit <- function(a, model = "pair_logit", pooling = "separate", mode = "separate", general = FALSE) {
  if (!is.null(a$reference) || general) osop_ns("fit_or_selection_models")(
    a$reader, a$views, a$control, osop_gc, 1e-10, a$recipe,
    selection_model = model, selection_pooling = pooling, reference = a$reference,
    outcome_pooling = mode) else osop_ns("fit_or_corollary1")(
    a$reader, a$views, a$control, osop_gc, 1e-10, a$recipe,
    selection_model = model, selection_pooling = pooling, outcome_pooling = mode)
}

# Only semantic cell/block IDs and frozen alias columns are borrowed. All
# nuisance roots, rates and final normalized targets are re-estimated from raw
# rows under the original-unit empirical frequencies f. Sampling w/rho stay fixed.
osop_reference <- function(a, fit, f = rep(1, nrow(a$Y)), nuisance = TRUE) {
  N <- length(f); nf <- if (nuisance) f else rep(1, N)
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  general <- isTRUE(fit$uses_selection_models)
  if (general) {
    cells <- fit$selection_backend$cells; blocks <- fit$selection_backend$block_specs
  } else {
    cells <- lapply(fit$support_specs, function(s) list(cohort = s$cohort, pair = s$pair, type = "pair"))
    blocks <- lapply(seq_along(cells), function(i) list(cell_ids = i))
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
      start = rep(0, length(ids) + 1L), family = stats::binomial(), control = osop_gc, intercept = FALSE))
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
    rate <- sum(nf[iu] * W[iu] * a$A[iu]) / sum(nf[iu] * W[iu])
    if (general) {
      pc <- fit$selection_backend$pair_cohorts[[fit$q_of_s[s]]]
      pp <- if (fit$selection_model == "independent_period_logit") pred[[pc$s_cell_id]] * pred[[pc$t_cell_id]] else pred[[pc$pair_cell_id]]
    } else pp <- pred[[s]]
    dy <- a$Y[, sp$pair$t_calendar_index] - a$Y[, sp$pair$s_calendar_index]
    ix <- which(ik & is.finite(dy)); q <- numeric(N)
    q[ix] <- W[ix] * rate / (a$rho[ix] * pmax(osop_ns("ipw_clip_eps"), pp[ix]))
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
  c(unlist(theta, use.names = FALSE), mu, den, num, unlist(design, use.names = FALSE), objective)
}

osop_point <- function(fit) {
  theta <- lapply(seq_along(fit$regressions), function(b) c(
    if (fit$outcome_recipe$intercept) fit$regressions[[b]]$alpha,
    if (fit$outcome_pooling == "separate") fit$outcome_groups$groups[[b]]$beta else fit$shared$beta))
  c(unlist(theta, use.names = FALSE), vapply(fit$links, `[[`, numeric(1), "mu"),
    vapply(fit$links, `[[`, numeric(1), "d"), vapply(fit$links, `[[`, numeric(1), "numerator"),
    unlist(lapply(fit$links, function(l) l$c * l$d), use.names = FALSE),
    if (fit$outcome_pooling == "separate") fit$outcome_groups$objective else fit$shared$objective)
}

osop_details <- function(fit) {
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  out <- matrix(NA_real_, N, length(osop_point(fit)))
  osop_ns("replay_or_corollary1_if")(fit, 11L, length(fit$links), function(o) {
    d <- o$details
    if (fit$outcome_pooling == "separate") {
      testthat::expect_identical(o$diagnostics$outcome_pooling, "separate")
      testthat::expect_identical(o$diagnostics$n_global_beta_blocks, 0L)
      testthat::expect_equal(o$diagnostics$n_local_beta_blocks, length(d$blocks))
    }
    theta <- if (fit$outcome_pooling == "separate") lapply(d$blocks, function(b)
      cbind(if (fit$outcome_recipe$intercept) b$alpha, b$beta)) else
        lapply(seq_along(fit$regressions), function(b)
          cbind(if (fit$outcome_recipe$intercept) d$alpha[, match(b, d$regression_ids)], d$beta))
    out[o$rows, ] <<- do.call(cbind, c(theta, list(o$phi, d$denominator, d$numerator), d$design, list(d$objective)))
  }, details = TRUE)
  out
}

testthat::test_that("conditional outcome modes preserve control-only WLS and pooled routes", {
  for (spec in list(c(TRUE, TRUE), c(FALSE, TRUE), c(TRUE, FALSE))) {
    a <- osop_case(intercept = spec[1L], slopes = spec[2L])
    for (mode in c("separate", "pooled")) {
      fit <- osop_fit(a, mode = mode)
      testthat::expect_equal(osop_point(fit), osop_reference(a, fit), tolerance = 2e-8)
      if (mode == "pooled") {
        default <- osop_ns("fit_or_corollary1")(a$reader, a$views, a$control, osop_gc, 1e-10, a$recipe)
        testthat::expect_identical(default$shared, fit$shared)
        testthat::expect_identical(default$links, fit$links)
        testthat::expect_identical(osop_ns("evaluate_or_corollary1_if_chunk")(default, 1:11, 1:3, details = TRUE),
          osop_ns("evaluate_or_corollary1_if_chunk")(fit, 1:11, 1:3, details = TRUE))
      } else {
        testthat::expect_null(fit$shared)
        testthat::expect_equal(fit$accounting$n_regression_fits, length(fit$regressions))
        actual <- osop_details(fit); h <- 1e-6; N <- nrow(actual)
        for (i in c(1L, 10L, 24L, 82L)) {
          delta <- -rep(1, N); delta[i] <- N - 1
          numerical <- (osop_reference(a, fit, 1 + h * delta) -
            osop_reference(a, fit, 1 - h * delta)) / (2 * h)
          testthat::expect_equal(unname(actual[i, ]), numerical, tolerance = 8e-6)
        }
      }
    }
  }
})

testthat::test_that("full-N conditional IFs refit every logit, repeat rate and OR block", {
  specs <- list(list(model = "pair_logit", pooling = "separate", u = "F"),
    list(model = "pair_logit", pooling = "within_cohort", u = "B"),
    list(model = "pair_logit", pooling = "across_cohorts", u = "A"),
    list(model = "independent_period_logit", pooling = "separate", u = "F"),
    list(model = "independent_period_logit", pooling = "within_cohort", u = "F"),
    list(model = "independent_period_logit", pooling = "across_cohorts", u = "F"))
  for (spec in specs) {
    a <- osop_case(spec$u); fit <- osop_fit(a, spec$model, spec$pooling)
    actual <- osop_details(fit); N <- nrow(actual); h <- 1e-6
    testthat::expect_equal(osop_point(fit), osop_reference(a, fit), tolerance = 2e-8)
    numerical <- t(vapply(seq_len(N), function(i) {
      delta <- -rep(1, N); delta[i] <- N - 1
      (osop_reference(a, fit, 1 + h * delta) - osop_reference(a, fit, 1 - h * delta)) / (2 * h)
    }, numeric(ncol(actual))))
    testthat::expect_equal(unname(actual), numerical, tolerance = 8e-6)
    testthat::expect_lt(max(abs(colMeans(actual))), 1e-8)
    if (spec$u != "A") testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
    testthat::expect_lte(fit$accounting$owned_entries, fit$accounting$retained_bound)
    # Independent nuisance freezing must give a different derivative. This
    # witnesses detection of omitted logit/rate forcing under identical inputs.
    delta <- -rep(1, N); delta[1L] <- N - 1
    frozen <- (osop_reference(a, fit, 1 + h * delta, FALSE) -
      osop_reference(a, fit, 1 - h * delta, FALSE)) / (2 * h)
    if (spec$u == "F") testthat::expect_gt(max(abs(actual[1L, ] - frozen)), 1e-4)
  }
})

testthat::test_that("legacy and general separate pair selection have equivalent fits and IFs", {
  a <- osop_case(); legacy <- osop_fit(a); general <- osop_fit(a, general = TRUE)
  testthat::expect_equal(osop_point(legacy), osop_point(general), tolerance = 2e-8)
  testthat::expect_equal(osop_details(legacy), osop_details(general), tolerance = 2e-8)
  # Repeated links reuse the actual block; nuisance cells are not refitted per link.
  frame <- get("frame_owner", a$reader$frame_ref, inherits = FALSE)
  extra <- osop_ns("construct_ipw_link_view")(frame, 2, 11, 11, c(23, Inf), control = a$control,
    link_id = "selection-pooling-repeat", functional_key = "cmpcat:selection-pooling-repeat")
  fit <- osop_ns("fit_or_corollary1")(a$reader, c(a$views, list(extra)), a$control, osop_gc,
    1e-10, a$recipe, outcome_pooling = "separate")
  testthat::expect_equal(fit$accounting$n_regression_fits, legacy$accounting$n_regression_fits)
  testthat::expect_equal(fit$accounting$n_selection_fits, legacy$accounting$n_selection_fits)
  testthat::expect_equal(fit$links[[4L]]$mu, fit$links[[1L]]$mu, tolerance = 1e-10)
  testthat::expect_equal(fit$links[[4L]]$regression, fit$links[[1L]]$regression)
  testthat::expect_error(osop_ns("or_corollary1_if_budget")(fit, 1e20, 1:4, TRUE))
})

testthat::test_that("dated conditional recipes retain anchored weights and full-N channels", {
  a <- osop_case(dated = TRUE)
  fit <- osop_fit(a, pooling = "across_cohorts")
  actual <- osop_details(fit); N <- nrow(actual); h <- 1e-6
  testthat::expect_equal(osop_point(fit), osop_reference(a, fit), tolerance = 2e-8)
  numerical <- t(vapply(seq_len(N), function(i) {
    delta <- -rep(1, N); delta[i] <- N - 1
    (osop_reference(a, fit, 1 + h * delta) - osop_reference(a, fit, 1 - h * delta)) / (2 * h)
  }, numeric(ncol(actual))))
  testthat::expect_equal(unname(actual), numerical, tolerance = 8e-6)
})
