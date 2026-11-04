# Independent empirical-frequency refits for stage-2 reconstructed outcomes.
rcp_ns <- function(name) getFromNamespace(name, "cdid")
rcp_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)

rcp_case <- function(universe = "F", intercept = TRUE, slopes = TRUE, dated = FALSE, route_kind = "corollary1", extra_reference = FALSE) {
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
  ctrl <- rcp_ns("ipw_owner_control")(11L, 1e8, 2e8, 1e13, 1e8)
  long <- data.frame(id = rep(id, each = 3), time = rep(calendar, N),
    y = as.vector(t(Y)), g = rep(g, each = 3))
  panel <- rcp_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  treatment_v <- grid$u + .3 * z
  V <- cbind(`(Intercept)` = 1, treatment_v = treatment_v)
  Xa <- if (intercept) cbind(X, treatment_v = treatment_v) else cbind(`(Intercept)` = 1, X, treatment_v = treatment_v)
  Xa[!A, ] <- NA_real_; Qa <- Q; Qa[!U, ] <- NA_real_
  frame <- rcp_ns("construct_ipw_frame_owner")(panel, id, Xa, w, ctrl,
    frame_g = g, design_on_A_only = TRUE)
  reader <- rcp_ns("construct_cp04_frame_reader")(frame, id, Qa, g, w, ctrl,
    U = U, rho = rho, universe = universe,
    justification = "known finite-law observation-subset transport ratios")
  defs <- list(list(s = 2, t = 11, g = 11, K = c(23, Inf)),
    list(s = 5, t = 11, g = 11, K = Inf), list(s = 2, t = 11, g = 23, K = Inf))
  views <- lapply(seq_along(defs), function(i) {
    a <- defs[[i]]
    rcp_ns("construct_ipw_link_view")(frame, a$s, a$t, a$g, a$K, control = ctrl,
      link_id = paste0("selection-pooling-", i), functional_key = paste0("cmpcat:selection-pooling-", i))
  })
  recipe <- rcp_ns("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  ref <- NULL; ref_X <- ref_Q <- ref_W <- ref_V <- NULL
  if (dated) {
    stopifnot(universe == "F", intercept, slopes)
    di <- rep(1:3, N)
    long$x <- rep(x, each = 3) + di * rep(z * grid$u / 7, each = 3)
    long$z <- rep(z, each = 3) * (1 + di / 20)
    long$treatment_v <- rep(treatment_v, each = 3) + di * rep(z / 30, each = 3)
    long$w <- rep(w, each = 3) * (1 + di * rep((grid$u + 1) / 10, each = 3))
    key <- match(long$time, calendar) * (N + 1) + long$id; wo <- order(key)
    route <- rcp_ns("chained_reference_route")(route_kind, "F"); route$or_basis <- TRUE
    formulas <- list(x = list(formula = ~x + treatment_v))
    if (!route_kind %in% c("baseline", "imputation")) formulas$s <- list(formula = ~z)
    reference_vars <- c("x", "z", "treatment_v")
    if (extra_reference) {
      long$unused_a <- rep(x^2 + .2 * grid$R2, each = 3) + di / 13
      long$unused_b <- rep(z * grid$u + .3 * grid$R1, each = 3) + di * rep(grid$R3 / 17, each = 3)
      formulas$x$formula <- ~x + treatment_v + unused_a + unused_b
      reference_vars <- c(reference_vars, "unused_a", "unused_b")
    }
    src <- rcp_ns("chained_reference_source")(formulas,
      data.frame(variable = reference_vars, kind = "dated", source = "data"), long,
      data.frame(id = id), FALSE, long$id, long$time, "time", id, calendar, 0,
      list(kind = "varying", name = "w", key = key[wo], value = long$w[wo], V = w,
        example = "finite-law dated reference weights"), "base_period", route)
    ref <- rcp_ns("build_chained_reference_reader")(src, reader, views, length(views), ctrl)
    recipe <- rcp_ns("or_pooled_outcome_recipe")(c("(Intercept)", "x"), TRUE, "x")
    ref_V <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, treatment_v = treatment_v + d * z / 30))
    ref_X <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, x = x + d * z * grid$u / 7))
    ref_Q <- lapply(ref$recipes$design_index, function(d) cbind(`(Intercept)` = 1, z = z * (1 + d / 20)))
    ref_W <- lapply(ref$recipes$weight_index, function(d) w * (1 + d * (grid$u + 1) / 10))
  }
  list(reader = reader, views = views, control = ctrl, recipe = recipe, reference = ref,
    X = X, V = V, X_full = if (intercept) cbind(X, treatment_v = treatment_v) else cbind(`(Intercept)` = 1, X, treatment_v = treatment_v), Q = Q, Y = Y, obs = obs, g = g, w = w, U = U, A = A,
    rho = rho, calendar = calendar, ref_X = ref_X, ref_Q = ref_Q, ref_W = ref_W, ref_V = ref_V)
}

rcp_fit <- function(a, method = "dr", pooling = "separate", mode = "separate") {
  rcp_ns("fit_imputation_backend")(reader = a$reader, link_views = a$views, control = a$control,
    glm_control = rcp_gc, rank_tol = 1e-10, outcome_recipe = a$recipe,
    treatment_columns = c("(Intercept)", "treatment_v"), est_method = method,
    imputation_pooling = pooling, X_frame = if (is.null(a$reference)) a$X_full else NULL,
    reference = a$reference, outcome_pooling = mode)
}

# Dense grouped WLS using only raw rows, semantic catalog IDs and frozen bases.
# The original empirical frequencies perturb every nuisance stage, not w or X.
rcp_reference <- function(a, fit, freq = rep(1, nrow(a$Y))) {
  N <- length(freq); dated <- !is.null(fit$layout)
  if (dated) {
    orig <- fit$layout$orig_row; coh <- fit$layout$state$cohort_index
    st <- fit$layout$state
    # Independent values at the recorded extraction/weight dates, then placed
    # on the semantic recipe/cohort layout. No fitted design values are read.
    X <- matrix(0, length(orig), ncol(a$X)); V <- matrix(0, length(orig), ncol(a$V)); W <- numeric(length(orig))
    for (r in seq_len(nrow(a$reference$recipes))) {
      ix <- which(fit$layout$recipe[coh] == r)
      X[ix, ] <- a$ref_X[[r]][orig[ix], , drop = FALSE]
      V[ix, ] <- a$ref_V[[r]][orig[ix], , drop = FALSE]
      W[ix] <- a$ref_W[[r]][orig[ix]]
    }
  } else {
    orig <- seq_len(N); st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
    coh <- st$cohort_index; X <- a$X; V <- a$V; W <- a$w
  }
  f <- freq[orig]; A <- a$A[orig]; n <- length(orig)
  aa <- vapply(seq_along(st$cohorts), function(k) sum(f[coh == k] * W[coh == k] * A[coh == k]) /
    sum(f[coh == k] * W[coh == k]), numeric(1))
  D <- function(rows, L, keep) {
    z <- L * X[rows, a$recipe$slope_index[keep], drop = FALSE]
    if (a$recipe$intercept) cbind(1, z) else z
  }
  jhat <- vector("list", length(fit$cells)); stage1 <- vector("list", length(fit$groups))
  for (g in seq_along(fit$groups)) {
    gf <- fit$groups[[g]]; cids <- gf$cell_ids; xx <- yy <- ww <- list()
    for (loc in seq_along(cids)) {
      cf <- fit$cells[[cids[loc]]]; rr <- which(coh == cf$cohort_index)
      dy <- a$Y[orig[rr], cf$pair$t_calendar_index] - a$Y[orig[rr], cf$pair$s_calendar_index]
      ix <- which(is.finite(dy)); rr <- rr[ix]
      z <- cf$L * X[rr, a$recipe$slope_index[gf$keep], drop = FALSE]
      if (a$recipe$intercept) {
        E <- matrix(0, length(rr), length(cids)); E[, loc] <- 1; z <- cbind(E, z)
      }
      xx[[loc]] <- z; yy[[loc]] <- dy[ix]; ww[[loc]] <- f[rr] * W[rr]
    }
    theta <- unname(stats::lm.wfit(do.call(rbind, xx), unlist(yy), unlist(ww), tol = 1e-10)$coefficients)
    stage1[[g]] <- theta
    beta <- theta[a$recipe$intercept * length(cids) + seq_along(gf$keep)]
    for (loc in seq_along(cids)) {
      cf <- fit$cells[[cids[loc]]]; rr <- which(coh == cf$cohort_index)
      alpha <- if (a$recipe$intercept) theta[loc] else 0
      val <- numeric(n); val[rr] <- alpha + as.vector((cf$L * X[rr, a$recipe$slope_index[gf$keep], drop = FALSE]) %*% beta)
      jhat[[cids[loc]]] <- val
    }
  }
  sep <- identical(fit$outcome_pooling, "separate"); Bn <- length(fit$or_blocks)
  stage2 <- list(); theta_b <- vector("list", Bn); keep_b <- vector("list", Bn)
  if (Bn) {
    pieces <- lapply(seq_len(Bn), function(b) {
      block <- fit$or_blocks[[b]]
      keep <- if (sep) fit$shared_or$outcome_groups$groups[[b]]$keep else fit$shared_or$keep
      keep_b[[b]] <<- keep; xx <- yy <- ww <- list()
      for (cid in block$control_cell_ids) {
        cf <- fit$cells[[cid]]; rr <- which(coh == cf$cohort_index)
        z <- block$L * X[rr, a$recipe$slope_index[keep], drop = FALSE]
        if (a$recipe$intercept) {
          E <- matrix(0, length(rr), if (sep) 1L else Bn); E[, if (sep) 1L else b] <- 1; z <- cbind(E, z)
        }
        xx[[length(xx) + 1L]] <- z; yy[[length(yy) + 1L]] <- jhat[[cid]][rr]
        ww[[length(ww) + 1L]] <- f[rr] * W[rr] * aa[cf$cohort_index]
      }
      list(X = do.call(rbind, xx), y = unlist(yy), w = unlist(ww))
    })
    solve_wls <- function(z) unname(stats::lm.wfit(z$X, z$y, z$w, tol = 1e-10)$coefficients)
    if (sep) {
      theta_b <- lapply(pieces, solve_wls); stage2 <- theta_b
    } else {
      root <- solve_wls(list(X = do.call(rbind, lapply(pieces, `[[`, "X")),
        y = unlist(lapply(pieces, `[[`, "y")), w = unlist(lapply(pieces, `[[`, "w"))))
      stage2 <- list(root)
      theta_b <- lapply(seq_len(Bn), function(b) c(if (a$recipe$intercept) root[b], root[a$recipe$intercept * Bn + seq_along(keep_b[[b]])]))
    }
  }
  gamma <- vector("list", length(fit$treatments))
  for (tid in seq_along(gamma)) {
    ts <- fit$treatments[[tid]]; rr <- which(coh %in% ts$all_cohort_indices)
    y <- as.double(coh[rr] == ts$target_cohort_index)
    root <- suppressWarnings(stats::glm.fit(V[rr, , drop = FALSE], y, weights = f[rr] * W[rr] * aa[coh[rr]],
      start = rep(0, ncol(V)), family = stats::binomial(), control = rcp_gc, intercept = FALSE))
    stopifnot(root$converged, all(is.finite(root$coefficients)))
    gamma[[tid]] <- unname(root$coefficients)
  }
  mu <- numeric(length(fit$links))
  for (ell in seq_along(fit$links)) {
    l <- fit$links[[ell]]; cf <- fit$cells[[l$target_cell_id]]; rr <- which(coh == cf$cohort_index)
    qw <- f[rr] * W[rr] * aa[cf$cohort_index]
    T <- sum(qw * jhat[[l$target_cell_id]][rr]) / sum(qw)
    Bg <- if (Bn) sum(qw * as.vector(D(rr, l$L, keep_b[[l$or_block_id]]) %*% theta_b[[l$or_block_id]])) / sum(qw) else 0
    if (fit$est_method == "reg") mu[ell] <- T - Bg else {
      rr <- which(coh %in% l$control_cohort_indices)
      pp <- stats::plogis(as.vector(V[rr, , drop = FALSE] %*% gamma[[l$trt_id]]))
      odds <- pmin(pp, 1 - 1e-6) / (1 - pmin(pp, 1 - 1e-6))
      qw <- f[rr] * W[rr] * aa[coh[rr]] * odds
      value <- numeric(length(rr))
      for (cid in l$control_cell_ids) {
        ix <- which(coh[rr] == fit$cells[[cid]]$cohort_index); value[ix] <- jhat[[cid]][rr[ix]]
      }
      C <- sum(qw * value) / sum(qw)
      if (fit$est_method == "ipw") mu[ell] <- T - C else {
        Bo <- sum(qw * as.vector(D(rr, l$L, keep_b[[l$or_block_id]]) %*% theta_b[[l$or_block_id]])) / sum(qw)
        mu[ell] <- (T - Bg) - (C - Bo)
      }
    }
  }
  c(unlist(stage1, use.names = FALSE), unlist(stage2, use.names = FALSE), unlist(gamma, use.names = FALSE), mu)
}

rcp_point <- function(fit) {
  sh <- fit$shared_or
  stage2 <- if (is.null(sh)) numeric() else if (identical(fit$outcome_pooling, "separate"))
    unlist(lapply(seq_len(sh$n_blocks), function(b) {
      g <- sh$outcome_groups$groups[[b]]
      c(if (fit$outcome_recipe$intercept) sh$blocks[[b]]$alpha, g$beta)
    }), use.names = FALSE) else sh$theta
  unname(c(unlist(lapply(fit$groups, function(g) c(g$alpha, g$beta)), use.names = FALSE), stage2,
    unlist(lapply(fit$treatments, `[[`, "coefficients"), use.names = FALSE), vapply(fit$links, `[[`, numeric(1), "mu")))
}

rcp_collect <- function(fit) {
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  out <- matrix(NA_real_, N, length(rcp_point(fit)))
  rcp_ns("replay_imputation_if")(fit, 11L, length(fit$links), function(o) {
    d <- o$details
    stage1 <- lapply(fit$groups, function(g) cbind(if (fit$outcome_recipe$intercept)
      d$imputation_alpha[, g$cell_ids, drop = FALSE], d$imputation_beta[[g$group_id]]))
    stage2 <- if (is.null(d$or_if)) list() else if (identical(fit$outcome_pooling, "separate"))
      lapply(d$or_if$blocks, `[[`, "theta") else list(d$or_if$theta)
    gamma <- lapply(d$trt_if, `[[`, "gamma")
    out[o$rows, ] <<- do.call(cbind, c(stage1, stage2, gamma, list(o$phi)))
  }, details = TRUE)
  out
}

testthat::test_that("reconstructed control OR modes retain stage1 fits and full finite-frequency IFs", {
  specs <- list(list(method = "reg", imp = "separate", mode = "separate"),
    list(method = "dr", imp = "within_cohort", mode = "separate"),
    list(method = "dr", imp = "across_cohorts", mode = "separate"),
    list(method = "dr", imp = "within_cohort", mode = "pooled"))
  for (spec in specs) {
    a <- rcp_case(); fit <- rcp_fit(a, spec$method, spec$imp, spec$mode)
    other <- rcp_fit(a, spec$method, spec$imp, if (spec$mode == "separate") "pooled" else "separate")
    testthat::expect_identical(fit$cells, other$cells)
    testthat::expect_identical(fit$groups, other$groups)
    testthat::expect_identical(fit$treatments, other$treatments)
    testthat::expect_equal(rcp_point(fit), rcp_reference(a, fit), tolerance = 2e-8)
    actual <- rcp_collect(fit); N <- nrow(actual); h <- 1e-6
    numerical <- t(vapply(seq_len(N), function(i) {
      delta <- -rep(1, N); delta[i] <- N - 1
      (rcp_reference(a, fit, 1 + h * delta) - rcp_reference(a, fit, 1 - h * delta)) / (2 * h)
    }, numeric(ncol(actual))))
    testthat::expect_equal(unname(actual), numerical, tolerance = 1e-5)
    testthat::expect_lt(max(abs(colMeans(actual))), 1e-8)
    testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-4)
    full <- rcp_ns("evaluate_imputation_if_chunk")(fit, 1:11, 1:3, details = TRUE)
    sub <- rcp_ns("evaluate_imputation_if_chunk")(fit, c(10L, 1L, 3L), c(3L, 1L))
    testthat::expect_equal(sub$phi, full$phi[c(10L, 1L, 3L), c(3L, 1L), drop = FALSE], tolerance = 1e-10)
  }
})

testthat::test_that("reconstruction handles no-intercept, intercept-only, and inactive IPW outcome choice", {
  for (spec in list(c(FALSE, TRUE), c(TRUE, FALSE))) {
    a <- rcp_case(intercept = spec[1L], slopes = spec[2L]); fit <- rcp_fit(a)
    testthat::expect_equal(rcp_point(fit), rcp_reference(a, fit), tolerance = 2e-8)
    actual <- rcp_collect(fit); h <- 1e-6; N <- nrow(actual)
    for (i in c(1L, 10L, 24L, 82L)) {
      delta <- -rep(1, N); delta[i] <- N - 1
      numerical <- (rcp_reference(a, fit, 1 + h * delta) - rcp_reference(a, fit, 1 - h * delta)) / (2 * h)
      testthat::expect_equal(unname(actual[i, ]), numerical, tolerance = 1e-5)
    }
    if (!spec[2L]) rcp_ns("replay_imputation_if")(fit, 11L, 3L, function(o)
      testthat::expect_equal(o$phi, o$details$psi_ipw, tolerance = 1e-10), details = TRUE)
  }
  a <- rcp_case()
  fit <- rcp_ns("fit_imputation_backend")(a$reader, a$views, a$control, rcp_gc, 1e-10, a$recipe,
    treatment_columns = c("(Intercept)", "treatment_v"), est_method = "ipw", X_frame = a$X_full,
    outcome_pooling = stop("inactive outcome choice must stay lazy"))
  testthat::expect_null(fit$shared_or)
  testthat::expect_identical(fit$outcome_pooling, "inactive")
})

testthat::test_that("separate reconstructed outcomes sum dated-layout channels by original unit", {
  a <- rcp_case(dated = TRUE, route_kind = "imputation"); fit <- rcp_fit(a, pooling = "within_cohort")
  testthat::expect_equal(rcp_point(fit), rcp_reference(a, fit), tolerance = 2e-8)
  N <- nrow(a$Y); h <- 1e-6
  wanted <- c(1L, 10L, 24L, 82L)
  actual <- matrix(NA_real_, length(wanted), length(fit$links))
  for (j in seq_along(wanted)) {
    r <- wanted[j]
    actual[j, ] <- rcp_ns("evaluate_imputation_if_chunk")(fit, r, 1:3)$phi
    delta <- -rep(1, N); delta[r] <- N - 1
    numerical <- tail((rcp_reference(a, fit, 1 + h * delta) - rcp_reference(a, fit, 1 - h * delta)) / (2 * h), 3L)
    testthat::expect_equal(unname(actual[j, ]), numerical, tolerance = 1e-5)
  }
  together <- rcp_ns("evaluate_imputation_if_chunk")(fit, wanted, c(3L, 1L))
  testthat::expect_equal(unname(together$phi), actual[, c(3L, 1L), drop = FALSE], tolerance = 1e-10)
  core <- rcp_ns("evaluate_imputation_if_chunk")(fit, wanted, c(3L, 1L), scale = "core")
  testthat::expect_equal(core$phi, (sum(a$A) / N) * together$phi, tolerance = 1e-10)
})
