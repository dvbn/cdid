# Independent deterministic baseline references. No external fixture files.
obp_ns <- function(name) getFromNamespace(name, "cdid")
obp_control <- function() obp_ns("ipw_owner_control")(7L, 1e8, 2e8, 1e12, 1e8)

obp_case <- function(intercept = TRUE, slopes = TRUE, alias = FALSE, dated = FALSE,
                     local_alias = FALSE, fail_alias = FALSE) {
  N <- 36L; j <- rep(1:12, 3); k <- rep(1:3, each = 12)
  id <- seq_len(N); g <- c(11, 23, Inf)[k]; calendar <- c(2, 5, 11, 23)
  x <- (j - 5) / 3 + k / 10; w <- 0.7 + (j %% 4) / 3 + k * (j %% 3) / 10
  Y <- vapply(1:4, function(d) 1 + j / 5 + k * d / 4 +
    ((0.15 + 0.35 * k) * d + 0.07 * d^2 * (k - 2)) * x +
    sin(j + d * k) / 10 + 0.03 * d * j^2, numeric(N))
  S <- matrix(TRUE, N, 4)
  S[j == 9, 3:4] <- FALSE; S[j == 10, c(1, 3)] <- FALSE
  S[j == 11, ] <- FALSE; S[j == 11, 2] <- TRUE; S[j == 12, ] <- FALSE
  S[j == 7 & k == 2, 3] <- FALSE; S[j == 8 & k == 3, 1] <- FALSE
  if (local_alias) {
    S[,] <- FALSE; S[j <= 5, c(1, 3)] <- TRUE
    S[j >= 6 & j <= 10, c(2, 3)] <- TRUE; S[j == 11, 2] <- TRUE
  }
  Y[!S] <- NA_real_; A <- rowSums(S) >= 2
  X <- matrix(numeric(), N, 0)
  if (intercept) X <- cbind(`(Intercept)` = rep(1, N))
  if (slopes) X <- cbind(X, x = x)
  if (alias) {
    ax <- if (local_alias) ifelse(j <= 5, 3 + 2 * x, 4 - x) else
      (if (intercept) 3 else 0) + 2 * x
    if (fail_alias) ax[j == 1 & k == 1] <- ax[j == 1 & k == 1] + 0.1
    X <- cbind(X, alias = ax)
  }
  ctrl <- obp_control()
  long <- data.frame(id = rep(id, each = 4), time = rep(calendar, N),
    y = as.vector(t(Y)), g = rep(g, each = 4))
  # One original unit is absent even from the long input; other zero and
  # singleton units remain explicitly represented there.
  panel <- obp_ns("construct_ipw_panel_owner")(long[long$id != N, ], "id", "time", "y", "g", ctrl)
  order <- c(seq.int(N, 2L, -2L), seq.int(N - 1L, 1L, -2L))
  Xa <- X; Xa[!A, ] <- NA_real_
  frame <- obp_ns("construct_ipw_frame_owner")(panel, id[order], Xa[order, , drop = FALSE],
    weights = w[order], control = ctrl, frame_g = g[order],
    allow_frame_only = TRUE, design_on_A_only = TRUE)
  reader <- obp_ns("construct_cp04_frame_reader")(frame, id[order],
    matrix(1, N, 1, dimnames = list(NULL, "selection_intercept")), g[order], w[order], ctrl)
  def <- list(list(s = 2, t = 11, g = 11, K = c(23, Inf)),
    list(s = 5, t = 11, g = 11, K = Inf), list(s = 2, t = 11, g = 23, K = Inf))
  views <- lapply(seq_along(def), function(i) {
    a <- def[[i]]
    obp_ns("construct_ipw_link_view")(frame, a$s, a$t, a$g, a$K, control = ctrl,
      link_id = paste0("private-", i), functional_key = paste0("cmpcat:private-", i))
  })
  recipe <- obp_ns("or_pooled_outcome_recipe")(colnames(X), intercept,
    setdiff(colnames(X), "(Intercept)"))
  ref <- NULL; ref_Z <- ref_W <- NULL
  if (dated) {
    stopifnot(intercept, slopes, !alias)
    long$x <- rep(x, each = 4) + rep(1:4, N) * rep(x^2 / 8, each = 4)
    long$w <- rep(w, each = 4) * (1 + rep(1:4, N) * rep((j %% 3) / 10, each = 4))
    key <- match(long$time, calendar) * (N + 1) + match(long$id, id[order])
    wo <- order(key)
    route <- obp_ns("chained_reference_route")("baseline"); route$or_basis <- TRUE
    src <- obp_ns("chained_reference_source")(list(x = list(formula = ~x)),
      data.frame(variable = "x", kind = "dated", source = "data"), long,
      data.frame(id = id[order]), FALSE, long$id, long$time, "time", id[order], calendar,
      0, list(kind = "varying", name = "w", key = key[wo], value = long$w[wo],
        V = w[order], example = "deterministic dated weights"), "base_period", route)
    ref <- obp_ns("build_chained_reference_reader")(src, reader, views, length(views), ctrl)
    recipe <- ref$outcome_recipe
    # Independent raw reference rows at the recorded design/weight dates.
    ref_Z <- lapply(ref$recipes$design_index, function(d)
      cbind(`(Intercept)` = 1, x = x[order] + d * x[order]^2 / 8))
    ref_W <- lapply(ref$recipes$weight_index, function(d)
      w[order] * (1 + d * (j[order] %% 3) / 10))
  }
  list(reader = reader, frame = frame, views = views, recipe = recipe, control = ctrl,
    reference = ref, X = X[order, , drop = FALSE], Y = Y[order, , drop = FALSE],
    g = g[order], w = w[order], A = A[order], id = id[order], calendar = calendar,
    ref_Z = ref_Z, ref_W = ref_W)
}

obp_fit <- function(a, mode = "separate", ...) obp_ns("fit_or_baseline")(
  a$reader, a$views, a$control, 1e-10, a$recipe, reference = a$reference,
  outcome_pooling = mode, ...)

obp_complete_case <- function() {
  N <- 12L; calendar <- c(2, 5, 11, 23); id <- seq_len(N)
  g <- rep(c(11, Inf), each = 6); x <- rep(c(-2, -1, 0, 1, 2, 3), 2)
  X <- cbind(`(Intercept)` = 1, x = x); w <- rep(1, N)
  Y <- vapply(1:4, function(d) x / 3 + d + (d + d^2 / 5) * x +
    cos(id + d) / 5 + as.double(g == 11) * d^2, numeric(N))
  long <- data.frame(id = rep(id, each = 4), time = rep(calendar, N),
    y = as.vector(t(Y)), g = rep(g, each = 4)); ctrl <- obp_control()
  panel <- obp_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  frame <- obp_ns("construct_ipw_frame_owner")(panel, id, X, w, ctrl)
  reader <- obp_ns("construct_cp04_frame_reader")(frame, id,
    matrix(1, N, 1, dimnames = list(NULL, "selection_intercept")), g, w, ctrl)
  views <- lapply(1:2, function(i) obp_ns("construct_ipw_link_view")(frame,
    c(2, 5)[i], c(11, 23)[i], 11, Inf, control = ctrl,
    link_id = paste0("complete-", i), functional_key = paste0("cmpcat:complete-", i)))
  list(reader = reader, frame = frame, views = views,
    recipe = obp_ns("or_pooled_outcome_recipe")(colnames(X), TRUE, "x"), control = ctrl,
    reference = NULL, X = X, Y = Y, g = g, w = w, A = rep(TRUE, N), id = id, calendar = calendar)
}

# Refit the complete raw-data functional at empirical unit frequencies f.
# Only the frozen alias column IDs and semantic catalog domains are borrowed;
# no backend coefficients, moments, forcing or IF operators enter this oracle.
obp_reference <- function(a, fit, f = rep(1, nrow(a$Y)), rate_path = TRUE) {
  N <- length(f); icpt <- a$recipe$intercept; separate <- fit$outcome_pooling == "separate"
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  G <- length(fit$regression_specs)
  ss <- lapply(fit$support_specs, function(s) {
    r <- s$reference; W <- if (is.null(r)) a$w else a$ref_W[[r]]
    Z <- if (is.null(r)) a$X else a$ref_Z[[r]]
    si <- match(s$pair$s, a$calendar); ti <- match(s$pair$t, a$calendar)
    dy <- a$Y[, ti] - a$Y[, si]; obs <- which(a$g == st$cohorts[s$cohort] & is.finite(dy))
    ik <- a$g == st$cohorts[s$cohort]
    rf <- if (rate_path) f else rep(1, N)
    rate <- sum(rf[obs] * W[obs]) / sum(rf[ik & a$A] * W[ik & a$A])
    q <- numeric(N); q[obs] <- W[obs] / max(rate, obp_ns("ipw_clip_eps"))
    list(Z = Z, dy = dy, obs = obs, rate = rate, q = q, L = ti - si)
  })
  pieces <- lapply(seq_len(G), function(b) {
    keep <- if (separate) fit$outcome_groups$groups[[b]]$keep else fit$shared$keep
    sids <- fit$regression_specs[[b]]$support_ids
    xx <- yy <- ww <- list()
    for (s in sids) {
      h <- ss[[s]]; ix <- h$obs
      d <- h$L * h$Z[ix, a$recipe$slope_index[keep], drop = FALSE]
      if (icpt) {
        E <- matrix(0, length(ix), if (separate) 1L else G)
        E[, if (separate) 1L else b] <- 1
        d <- cbind(E, d)
      }
      xx[[length(xx) + 1L]] <- d; yy[[length(yy) + 1L]] <- h$dy[ix]
      ww[[length(ww) + 1L]] <- f[ix] * h$q[ix]
    }
    list(X = do.call(rbind, xx), y = unlist(yy), w = unlist(ww), keep = keep)
  })
  solve <- function(z) unname(stats::lm.wfit(z$X, z$y, z$w, tol = 1e-10)$coefficients)
  if (separate) theta <- lapply(pieces, solve) else {
    root <- solve(list(X = do.call(rbind, lapply(pieces, `[[`, "X")),
      y = unlist(lapply(pieces, `[[`, "y")), w = unlist(lapply(pieces, `[[`, "w"))))
    theta <- lapply(seq_len(G), function(b) c(if (icpt) root[b], root[icpt * G + seq_along(pieces[[b]]$keep)]))
  }
  residual <- function(b, s) {
    h <- ss[[s]]; ix <- h$obs
    D <- h$L * h$Z[ix, a$recipe$slope_index[pieces[[b]]$keep], drop = FALSE]
    if (icpt) D <- cbind(1, D)
    list(ix = ix, r = h$dy[ix] - as.vector(D %*% theta[[b]]))
  }
  objective <- 0
  for (b in seq_len(G)) for (s in fit$regression_specs[[b]]$support_ids) {
    e <- residual(b, s); objective <- objective + sum(f[e$ix] * ss[[s]]$q[e$ix] * e$r^2) / N
  }
  mu <- den <- num <- numeric(length(fit$links)); design <- list()
  for (l in seq_along(fit$links)) {
    link <- fit$links[[l]]; s <- link$support_id; b <- link$regression
    e <- residual(b, s); ix <- e$ix; qw <- f[ix] * ss[[s]]$q[ix]
    den[l] <- sum(qw) / N; num[l] <- sum(qw * e$r) / N; mu[l] <- num[l] / den[l]
    D <- ss[[s]]$L * ss[[s]]$Z[ix, a$recipe$slope_index[pieces[[b]]$keep], drop = FALSE]
    design[[l]] <- colSums(qw * D) / N
  }
  shares <- if (is.null(a$reference)) c(
    vapply(st$cohorts, function(k) sum(f * a$w * (a$g == k)) / sum(f * a$w), numeric(1)),
    vapply(st$cohorts, function(k) sum(f * a$w * a$A * (a$g == k)) /
      sum(f * a$w * a$A), numeric(1))) else numeric()
  c(unlist(theta, use.names = FALSE), vapply(ss, `[[`, numeric(1), "rate"),
    mu, den, num, unlist(design, use.names = FALSE), objective, shares, sum(f * a$A) / N)
}

obp_point <- function(fit) {
  separate <- fit$outcome_pooling == "separate"; icpt <- fit$outcome_recipe$intercept
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  theta <- lapply(seq_along(fit$regressions), function(b) {
    beta <- if (separate) fit$outcome_groups$groups[[b]]$beta else fit$shared$beta
    c(if (icpt) fit$regressions[[b]]$alpha, beta)
  })
  c(unlist(theta, use.names = FALSE), vapply(fit$supports, `[[`, numeric(1), "rate"),
    vapply(fit$links, `[[`, numeric(1), "mu"), vapply(fit$links, `[[`, numeric(1), "d"),
    vapply(fit$links, `[[`, numeric(1), "numerator"),
    unlist(lapply(fit$links, function(l) l$c * l$d), use.names = FALSE),
    if (separate) fit$outcome_groups$objective else fit$shared$objective,
    if (is.null(fit$reference_ref)) c(st$tau, st$eta), st$lambda)
}

obp_detail_matrix <- function(fit) {
  N <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state$N
  answer <- matrix(NA_real_, N, length(obp_point(fit)))
  obp_ns("replay_or_baseline_if")(fit, 7L, length(fit$links), function(o) {
    d <- o$details
    theta <- if (fit$outcome_pooling == "separate") lapply(d$blocks, function(b)
      cbind(if (fit$outcome_recipe$intercept) b$alpha, b$beta)) else
      lapply(seq_along(fit$regressions), function(b)
        cbind(if (fit$outcome_recipe$intercept) d$alpha[, match(b, d$regression_ids)], d$beta))
    answer[o$rows, ] <<- do.call(cbind, c(theta, list(d$rates, o$phi, d$denominator,
      d$numerator), d$design, list(d$objective, d$tau, d$eta, d$lambda)))
  }, details = TRUE)
  answer
}

testthat::test_that("baseline pooling fits independent calibrated WLS on the modeled calendar", {
  for (spec in list(c(TRUE, TRUE), c(FALSE, TRUE), c(TRUE, FALSE))) {
    a <- obp_case(spec[1L], spec[2L]); sep <- obp_fit(a); pooled <- obp_fit(a, "pooled")
    testthat::expect_equal(obp_point(sep), obp_reference(a, sep), tolerance = 1e-9)
    testthat::expect_equal(obp_point(pooled), obp_reference(a, pooled), tolerance = 1e-9)
    default <- obp_ns("fit_or_baseline")(a$reader, a$views, a$control, 1e-10, a$recipe)
    testthat::expect_identical(default$outcome_pooling, "pooled")
    testthat::expect_identical(default$shared, pooled$shared)
    testthat::expect_identical(default$links, pooled$links)
    testthat::expect_identical(obp_ns("evaluate_or_baseline_if_chunk")(default, 1:7,
      seq_along(default$links), details = TRUE), obp_ns("evaluate_or_baseline_if_chunk")(
      pooled, 1:7, seq_along(pooled$links), details = TRUE))
    testthat::expect_null(sep$shared); testthat::expect_null(sep$complete_support)
    testthat::expect_equal(sep$accounting$n_regression_fits, length(sep$regressions))
    testthat::expect_lte(sep$accounting$owned_entries, sep$accounting$retained_bound)
    if (spec[2L]) testthat::expect_gt(max(abs(vapply(sep$links, `[[`, numeric(1), "mu") -
      vapply(pooled$links, `[[`, numeric(1), "mu"))), 1e-3)
  }
})

testthat::test_that("separate bypasses complete pooled fitting while pooled native IF stays available", {
  a <- obp_complete_case(); sep <- obp_fit(a); pooled <- obp_fit(a, "pooled")
  testthat::expect_null(sep$complete_support)
  testthat::expect_identical(pooled$complete_support$route, "complete_support")
  testthat::expect_equal(obp_point(sep), obp_reference(a, sep), tolerance = 1e-9)
  testthat::expect_equal(obp_point(pooled), obp_reference(a, pooled), tolerance = 1e-9)
  for (fit in list(sep, pooled)) {
    plain <- obp_ns("evaluate_or_baseline_if_chunk")(fit, 1:7, 1:2)
    detailed <- obp_ns("evaluate_or_baseline_if_chunk")(fit, 1:7, 1:2, details = TRUE)
    testthat::expect_equal(plain$phi, detailed$phi, tolerance = 1e-10)
  }
})

testthat::test_that("full original-unit derivatives include all baseline ratio and rate paths", {
  for (spec in list(c(TRUE, TRUE), c(FALSE, TRUE), c(TRUE, FALSE))) {
    a <- obp_case(spec[1L], spec[2L])
    for (mode in c("separate", "pooled")) {
      fit <- obp_fit(a, mode); actual <- obp_detail_matrix(fit); N <- nrow(actual); h <- 1e-6
      numerical <- t(vapply(seq_len(N), function(i) {
        delta <- -rep(1, N); delta[i] <- N - 1
        (obp_reference(a, fit, 1 + h * delta) - obp_reference(a, fit, 1 - h * delta)) / (2 * h)
      }, numeric(ncol(actual))))
      testthat::expect_equal(unname(actual), numerical, tolerance = 2e-6)
      testthat::expect_lt(max(abs(colMeans(actual))), 1e-9)
      # Includes singleton, explicitly zero-outcome, and frame-only rows.
      testthat::expect_equal(sum(!a$A), 6L)
      testthat::expect_gt(max(abs(actual[!a$A, , drop = FALSE])), 1e-3)
      if (mode == "separate" && spec[2L]) {
        no_rates <- t(vapply(seq_len(N), function(i) {
          delta <- -rep(1, N); delta[i] <- N - 1
          (obp_reference(a, fit, 1 + h * delta, FALSE) -
            obp_reference(a, fit, 1 - h * delta, FALSE)) / (2 * h)
        }, numeric(ncol(actual))))
        nc <- sum(vapply(fit$outcome_groups$groups,
          function(g) length(g$keep) + as.integer(a$recipe$intercept), integer(1)))
        coefficient_and_link <- c(seq_len(nc), nc + length(fit$supports) + seq_along(fit$links))
        testthat::expect_gt(max(abs(actual[, coefficient_and_link] -
          no_rates[, coefficient_and_link])), 1e-3)
      }
    }
  }
})

testthat::test_that("separate block mapping survives overlapping controls, links and chunking", {
  a <- obp_case(); fit <- obp_fit(a); N <- length(a$id); L <- length(fit$links)
  all <- matrix(NA_real_, N, L); batched <- all
  obp_ns("replay_or_baseline_if")(fit, 7L, L, function(o) all[o$rows, o$links] <<- o$phi)
  obp_ns("replay_or_baseline_if")(fit, 3L, 1L, function(o) batched[o$rows, o$links] <<- o$phi)
  testthat::expect_equal(all, batched, tolerance = 1e-12)
  core <- obp_ns("evaluate_or_baseline_if_chunk")(fit, 7:1, c(3L, 1L), "core", TRUE)
  testthat::expect_equal(unname(core$phi), get("cp04_reader", fit$reader_ref)$state$lambda *
    all[7:1, c(3L, 1L)], tolerance = 1e-12)
  testthat::expect_identical(core$id, a$id[7:1])
  testthat::expect_identical(core$details$link_block_ids,
    vapply(fit$links[c(3L, 1L)], `[[`, integer(1), "regression"))
  testthat::expect_null(core$details$beta)
  local <- obp_ns("or_baseline_if_budget")(fit, 7L, 2L)
  testthat::expect_identical(local$gids, fit$links[[2L]]$regression)
  testthat::expect_equal(length(local$beta_gids), 1L)
  for (b in core$details$blocks) {
    testthat::expect_identical(b$slope_columns,
      fit$outcome_groups$groups[[b$group_id]]$basis$kept_columns)
    testthat::expect_equal(nrow(b$beta), 7L)
  }
  testthat::expect_lte(unname(obp_ns("cp04_payload")(core)["entries"]),
    core$diagnostics$output_entries_bound)
  # A new link ID reuses the same actual comparison fit exactly once.
  duplicate <- obp_ns("construct_ipw_link_view")(a$frame, 2, 11, 11, c(23, Inf), control = a$control,
    link_id = "private-reused", functional_key = "cmpcat:private-reused")
  a$views <- c(a$views, list(duplicate)); twice <- obp_fit(a)
  testthat::expect_identical(twice$outcome_groups$groups, fit$outcome_groups$groups)
  testthat::expect_identical(twice$links[[1L]]$regression, twice$links[[4L]]$regression)
  o <- obp_ns("evaluate_or_baseline_if_chunk")(twice, 1:7, c(1L, 4L))
  testthat::expect_identical(o$phi[, 1L], o$phi[, 2L])
  # Shared original rows induce cross-block covariance.
  testthat::expect_gt(abs(crossprod(all)[1L, 3L]), 1e-3)
})

testthat::test_that("local aliases preserve required predictions and seals and budgets still refuse", {
  a <- obp_case(alias = TRUE); fit <- obp_fit(a)
  testthat::expect_true(all(vapply(fit$outcome_groups$groups,
    function(g) identical(g$basis$aliased_columns, "alias"), logical(1))))
  testthat::expect_equal(obp_point(fit), obp_reference(a, fit), tolerance = 1e-9)
  bad <- fit; bad$outcome_groups$groups[[1L]]$beta <- 0
  testthat::expect_error(obp_ns("evaluate_or_baseline_if_chunk")(bad, 1:7, 1L),
    class = "cdid_or_baseline_integrity_error")
  testthat::expect_error(obp_fit(a, NA_character_), class = "cdid_or_baseline_input_error")
  testthat::expect_error(obp_ns("or_baseline_native_plan")(fit, 1L, 7L),
    class = "cdid_or_baseline_integrity_error")
  for (field in c("max_retained_entries", "max_matrix_entries", "max_work_ops", "max_output_entries")) {
    low <- a; low$control[[field]] <- 1
    testthat::expect_error(obp_fit(low), class = "cdid_ipw_owner_bound_error")
    limited <- fit; limited$ref <- NULL; limited$control[[field]] <- 1
    limited <- obp_ns("ipw_owner_seal")(limited, "or_baseline_backend"); called <- FALSE
    testthat::expect_error(obp_ns("replay_or_baseline_if")(limited, 7L, 1L,
      function(o) called <<- TRUE), class = "cdid_ipw_owner_bound_error")
    testthat::expect_false(called)
  }
  one <- obp_ns("or_baseline_if_budget")(fit, 7L, seq_along(fit$links))
  limited <- fit; limited$ref <- NULL; limited$control$max_work_ops <- ceiling(2 * one$work)
  limited <- obp_ns("ipw_owner_seal")(limited, "or_baseline_backend"); called <- FALSE
  testthat::expect_error(obp_ns("replay_or_baseline_if")(limited, 7L, length(fit$links),
    function(o) called <<- TRUE), class = "cdid_ipw_owner_bound_error")
  testthat::expect_false(called)
  prefix <- obp_fit(a, retained_external_entries = 1000, external_work_ops = 12345)
  testthat::expect_identical(prefix$accounting$work_ops_convention, "cumulative")
  testthat::expect_equal(prefix$accounting$work_ops_bound, prefix$accounting$own_work_ops + 12345)
  testthat::expect_equal(prefix$accounting$own_work_ops, fit$accounting$own_work_ops)
  testthat::expect_equal(prefix$accounting$retained_bound, fit$accounting$retained_bound + 1000)
})

testthat::test_that("block-local alias maps differ and failing treated predictions still refuse", {
  a <- obp_case(alias = TRUE, local_alias = TRUE); fit <- obp_fit(a)
  testthat::expect_equal(obp_point(fit), obp_reference(a, fit), tolerance = 1e-9)
  # The raw fixture has alias = 3 + 2*x on pair (2, 11), and
  # alias = 4 - x on pair (5, 11). Compare each fitted block, without
  # using exact floating-point uniqueness to merge distinct control groups.
  expected_alias <- vapply(fit$regression_specs, function(b) {
    s <- fit$support_specs[[b$support_ids[1L]]]$pair$s
    stopifnot(s %in% c(2, 5))
    if (s == 2) 2 else -1
  }, numeric(1))
  testthat::expect_equal(vapply(fit$outcome_groups$groups,
    function(g) unname(g$basis$alias_map[1L, 1L]), numeric(1)), expected_alias,
    tolerance = 1e-12)
  actual <- obp_detail_matrix(fit); N <- nrow(actual); h <- 1e-6
  numerical <- t(vapply(seq_len(N), function(i) {
    delta <- -rep(1, N); delta[i] <- N - 1
    (obp_reference(a, fit, 1 + h * delta) - obp_reference(a, fit, 1 - h * delta)) / (2 * h)
  }, numeric(ncol(actual))))
  testthat::expect_equal(unname(actual), numerical, tolerance = 2e-6)
  testthat::expect_error(obp_fit(obp_case(alias = TRUE, local_alias = TRUE, fail_alias = TRUE)),
    class = "cdid_or_baseline_support_error")
  ni <- obp_case(intercept = FALSE, alias = TRUE); nf <- obp_fit(ni)
  testthat::expect_true(all(vapply(nf$outcome_groups$groups,
    function(g) identical(g$basis$aliased_columns, "alias"), logical(1))))
  testthat::expect_equal(obp_point(nf), obp_reference(ni, nf), tolerance = 1e-9)
})

testthat::test_that("dated baseline reference weights retain their full-frequency forcing", {
  a <- obp_case(dated = TRUE); fit <- obp_fit(a)
  testthat::expect_equal(obp_point(fit), obp_reference(a, fit), tolerance = 1e-9)
  actual <- obp_detail_matrix(fit); N <- nrow(actual); h <- 1e-6
  numerical <- t(vapply(seq_len(N), function(i) {
    delta <- -rep(1, N); delta[i] <- N - 1
    (obp_reference(a, fit, 1 + h * delta) - obp_reference(a, fit, 1 - h * delta)) / (2 * h)
  }, numeric(ncol(actual))))
  testthat::expect_equal(unname(actual), numerical, tolerance = 2e-6)
  testthat::expect_null(fit$complete_support)
  testthat::expect_identical(fit$reference_ref, a$reference$ref)
})
