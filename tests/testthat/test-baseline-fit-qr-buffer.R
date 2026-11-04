# Independent raw weighted objectives/refits plus the unchanged unbuffered route.
qbf <- function(name) getFromNamespace(name, "cdid")
qbf_control <- function(matrix = 2e8, chunk = 3L) qbf("ipw_owner_control")(
  chunk, 1e8, matrix, 1e13, 1e8)
qbf_case <- function(intercept = TRUE, slopes = TRUE, alias = FALSE, bad_prediction = FALSE, complete = FALSE) {
  j <- rep(1:24, 3); group <- rep(1:3, each = 24); N <- if (complete) 72L else 74L; id <- seq_len(N)
  g <- c(c(3, 4, Inf)[group], 3, Inf)[seq_len(N)]
  x <- c(ifelse(j <= 8L, 0, (j - 16.5) / 6), .1, -.2)[seq_len(N)]
  if (bad_prediction) x[g == Inf] <- 0
  z <- c(sin(j) / 2 + j / 40, .2, -.3)[seq_len(N)]
  w <- c((1 + (j %% 4L) / 10) * (1 + group * x[1:72] / 40), 1.2, 1.7)[seq_len(N)]
  Y <- vapply(1:5, function(t) c(j / 7 + group * t / 6 + (.2 * t + .03 * t^2) * x[1:72] +
    sin(j + group * t) / 9 + as.double(t >= g[1:72]) * (t + group * x[1:72] / 7),
    NA_real_, NA_real_)[seq_len(N)], numeric(N))
  if (!complete) {
    Y[which(j == 1L), 2L] <- NA_real_; Y[which(j == 3L), 4L] <- NA_real_
    Y[N, 1L] <- 4.1
  }
  X <- matrix(numeric(), N, 0L)
  if (intercept) X <- cbind(`(Intercept)` = rep(1, N))
  if (slopes) X <- cbind(X, x = x)
  if (alias) X <- cbind(X, alias = 3 + 2 * x)
  Xa <- cbind(`(Intercept)` = 1, x = x, z = z)
  if (alias) Xa <- cbind(Xa, alias = X[, "alias"])
  V <- cbind(`(Intercept)` = 1, z = z)
  control <- qbf_control()
  data <- data.frame(id = rep(id, each = 5L), t = rep(1:5, N), g = rep(g, each = 5L), y = as.vector(t(Y)))
  if (!complete) {
    data <- data[!(data$id == 2L & data$t == 3L), ]; Y[2L, 3L] <- NA_real_
  }
  panel <- qbf("construct_ipw_panel_owner")(data, "id", "t", "y", "g", control)
  make_reader <- function(design) {
    frame <- qbf("construct_ipw_frame_owner")(panel, id, design, w, control)
    qbf("construct_cp04_frame_reader")(frame, id,
      matrix(1, N, 1L, dimnames = list(NULL, "selection_intercept")), g, w, control)
  }
  readers <- list(reg = make_reader(X), dr = make_reader(Xa))
  defs <- list(list(s = 1, t = 2, g = 3, K = c(4, Inf)), list(s = 1, t = 3, g = 3, K = c(4, Inf)),
    list(s = 2, t = 4, g = 3, K = Inf), list(s = 1, t = 3, g = 4, K = Inf),
    list(s = 3, t = 5, g = 3, K = c(4, Inf)), list(s = 1, t = 3, g = 3, K = Inf))
  views <- lapply(readers, function(reader) {
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    lapply(seq_along(defs), function(i) {
      d <- defs[[i]]
      qbf("construct_ipw_link_view")(frame, d$s, d$t, d$g, d$K, control = control,
        link_id = paste0("qr-buffer-", i), functional_key = paste0("cmpcat:qr-buffer-", i))
    })
  })
  recipe <- qbf("or_pooled_outcome_recipe")(colnames(X), intercept, setdiff(colnames(X), "(Intercept)"))
  list(readers = readers, views = views, control = control, recipe = recipe, X = X, V = V,
    Y = Y, g = g, w = w, defs = defs, N = N, intercept = intercept)
}
qbf_fit <- function(a, method, automatic = FALSE, control = a$control, mode = "separate") {
  args <- list(reader = a$readers[[method]], link_views = a$views[[method]], control = control,
    rank_tol = 1e-10, outcome_recipe = a$recipe, outcome_pooling = mode,
    automatic_chunk_rows = automatic)
  if (method == "dr") args <- c(args, list(glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE),
    treatment_columns = colnames(a$V)))
  do.call(qbf(if (method == "reg") "fit_or_baseline" else "fit_dr_baseline"), args)
}
qbf_eval <- function(fit, method, rows, links) qbf(if (method == "reg")
  "evaluate_or_baseline_if_chunk" else "evaluate_dr_if_chunk")(fit, rows, links)

# Independent WLS and GLM refits of every input comparison on raw observed rows.
qbf_raw <- function(a, method, frequency = rep(1, a$N)) {
  w <- a$w * frequency; A <- rowSums(!is.na(a$Y)) >= 2L
  vapply(a$defs, function(d) {
    obs <- !is.na(a$Y[, d$s]) & !is.na(a$Y[, d$t]); cohorts <- c(d$g, d$K)
    rates <- vapply(cohorts, function(k) sum(w[a$g == k & obs]) / sum(w[a$g == k & A]), numeric(1))
    rows <- which(obs & a$g %in% cohorts); D <- as.double(a$g[rows] == d$g)
    weights <- w[rows] / rates[match(a$g[rows], cohorts)]
    dy <- a$Y[rows, d$t] - a$Y[rows, d$s]
    design <- a$X[rows, , drop = FALSE]
    if (a$intercept) design[, -1L] <- design[, -1L, drop = FALSE] * (d$t - d$s) else
      design <- design * (d$t - d$s)
    fit <- stats::lm.wfit(design[D == 0, , drop = FALSE], dy[D == 0], weights[D == 0], tol = 1e-10)
    stopifnot(all(is.finite(fit$coefficients)))
    residual <- dy - as.vector(design %*% fit$coefficients)
    if (method == "reg") return(weighted.mean(residual[D == 1], weights[D == 1]))
    trt <- suppressWarnings(stats::glm.fit(a$V[rows, , drop = FALSE], D, weights = weights,
      start = rep(0, ncol(a$V)), family = stats::binomial(),
      control = stats::glm.control(epsilon = 1e-13, maxit = 60L), intercept = FALSE))
    stopifnot(trt$converged, all(is.finite(trt$coefficients)))
    prob <- stats::plogis(as.vector(a$V[rows, , drop = FALSE] %*% trt$coefficients))
    stopifnot(all(prob > 1e-6 & prob < 1 - 1e-6))
    weighted.mean(residual, weights * D) - weighted.mean(residual, weights * (1 - D) * prob / (1 - prob))
  }, numeric(1))
}

testthat::test_that("QR buffers preserve original scalar moments, owner row order and raw weighted objectives", {
  a <- qbf_case(); fit <- qbf_fit(a, "reg"); reader <- a$readers$reg
  spec <- fit$support_specs[[1L]]
  original <- qbf("or_baseline_visit"); log <- list()
  testthat::local_mocked_bindings(or_baseline_visit = function(reader, spec, ctrl, visit, ...) {
    original(reader, spec, ctrl, function(rows, S, obs, dy, Z) {
      log[[length(log) + 1L]] <<- list(rows = rows, S = S, obs = obs, dy = dy, Z = Z)
      visit(rows, S, obs, dy, Z)
    }, ...)
  }, .package = "cdid")
  old <- qbf("or_baseline_support")(reader, spec, a$recipe, a$control); before <- log; log <- list()
  new <- qbf("or_baseline_qr_buffered_support")(reader, spec, a$recipe, a$control, 7L)
  testthat::expect_identical(log, before)
  scalar <- function(x) {
    x$design$summary <- NULL; x
  }
  testthat::expect_identical(scalar(new), scalar(old))
  testthat::expect_lt(new$design$summary$accounting$n_compressions, old$design$summary$accounting$n_compressions)
  st <- reader$state; k <- st$cohorts[spec$cohort]
  ix <- which(a$g == k & !is.na(a$Y[, spec$pair$s]) & !is.na(a$Y[, spec$pair$t]))
  dy <- a$Y[ix, spec$pair$t] - a$Y[ix, spec$pair$s]
  D <- cbind(1, spec$transitions * a$X[ix, "x"])
  F <- new$design$summary$raw$factor
  for (theta in list(c(-.3, .4), c(1.7, -.2))) testthat::expect_equal(
    sum(as.vector(F %*% c(-theta, 1))^2), sum(a$w[ix] * (dy - as.vector(D %*% theta))^2), tolerance = 1e-9)
})

testthat::test_that("automatic buffered final fits preserve aliases, rank, IFs and full original-unit covariance", {
  for (options in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE))) {
    a <- do.call(qbf_case, options)
    for (method in c("reg", "dr")) {
      old <- qbf_fit(a, method); new <- qbf_fit(a, method, TRUE); L <- length(new$links)
      testthat::expect_null(new$separate_complete)
      testthat::expect_identical(new$control, old$control)
      testthat::expect_identical(new$rank_tol, old$rank_tol)
      testthat::expect_gt(new$accounting$fit_scratch_bound, old$accounting$fit_scratch_bound)
      testthat::expect_identical(lapply(new$outcome_groups$groups, `[[`, "keep"),
        lapply(old$outcome_groups$groups, `[[`, "keep"))
      testthat::expect_equal(vapply(new$links, `[[`, numeric(1), "mu"),
        vapply(old$links, `[[`, numeric(1), "mu"), tolerance = 1e-9)
      left <- right <- matrix(0, a$N, L)
      for (first in seq.int(1L, a$N, by = 3L)) {
        rows <- seq.int(first, min(a$N, first + 2L))
        left[rows, ] <- qbf_eval(new, method, rows, seq_len(L))$phi
        right[rows, ] <- qbf_eval(old, method, rows, seq_len(L))$phi
      }
      testthat::expect_equal(left, right, tolerance = 1e-9)
      testthat::expect_equal(crossprod(left) / a$N, crossprod(right) / a$N, tolerance = 1e-9)
      testthat::expect_equal(unname(qbf_eval(new, method, c(74L, 73L, 2L), c(6L, 3L, 1L))$phi),
        left[c(74L, 73L, 2L), c(6L, 3L, 1L), drop = FALSE], tolerance = 1e-10)
    }
  }
})

testthat::test_that("buffered nuisance fits match independent raw original-unit case-weight derivatives", {
  a <- qbf_case()
  for (method in c("reg", "dr")) {
    fit <- qbf_fit(a, method, TRUE); L <- length(fit$links)
    testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), qbf_raw(a, method), tolerance = 1e-9)
    actual <- matrix(0, a$N, L)
    for (first in seq.int(1L, a$N, by = 3L)) {
      rows <- seq.int(first, min(a$N, first + 2L)); actual[rows, ] <- qbf_eval(fit, method, rows, seq_len(L))$phi
    }
    numerical <- t(vapply(seq_len(a$N), function(i) {
      direction <- rep(-1, a$N); direction[i] <- a$N - 1
      (qbf_raw(a, method, 1 + 1e-6 * direction) - qbf_raw(a, method, 1 - 1e-6 * direction)) / 2e-6
    }, numeric(L)))
    testthat::expect_equal(actual, numerical, tolerance = 2e-6)
  }
})

testthat::test_that("buffer capacity and original numerical refusals retain their existing paths", {
  p <- 1; used <- 600; fixed <- 128 * (p + 2)^2 + 4096; slope <- 64 * (p + 4)
  edge <- used + fixed + slope * 4
  testthat::expect_null(qbf("or_baseline_qr_buffer_plan")(qbf_control(edge - 1), 100, 200, 300, p, 72))
  plan <- qbf("or_baseline_qr_buffer_plan")(qbf_control(edge), 100, 200, 300, p, 72)
  testthat::expect_identical(plan$rows, 4L)
  testthat::expect_equal(used + plan$extra, edge)
  a <- qbf_case(); n <- 0L; original <- qbf("or_baseline_qr_buffered_support")
  testthat::local_mocked_bindings(or_baseline_qr_buffered_support = function(...) {
    n <<- n + 1L; original(...)
  }, .package = "cdid")
  plain <- qbf_fit(a, "reg", FALSE); testthat::expect_identical(n, 0L)
  # Original fit fits exactly, but no extra buffer scratch remains.
  cap <- ceiling(plain$accounting$retained_bound + plain$accounting$fit_scratch_bound)
  old_tight <- qbf_fit(a, "reg", FALSE, qbf_control(cap))
  new_tight <- qbf_fit(a, "reg", TRUE, qbf_control(cap))
  testthat::expect_identical(n, 0L)
  testthat::expect_identical(vapply(new_tight$links, `[[`, numeric(1), "mu"),
    vapply(old_tight$links, `[[`, numeric(1), "mu"))
  qbf_fit(a, "reg", TRUE); testthat::expect_gt(n, 0L)
  before <- n
  for (method in c("reg", "dr")) {
    at_complete <- n
    complete <- qbf_case(complete = TRUE)
    testthat::expect_false(is.null(qbf_fit(complete, method, TRUE)$separate_complete))
    testthat::expect_identical(n, at_complete)
    at_pooled <- n
    pooled <- qbf_fit(a, method, TRUE, mode = "pooled")
    testthat::expect_identical(n - at_pooled, length(pooled$support_specs))
  }
  testthat::expect_gt(n, before)
  # The OR public adapter forwards its existing automatic flag.
  before <- n
  qbf("or_public_fit")("baseline", a$readers$reg, a$views$reg, list(owner_control = a$control),
    a$readers$reg$accounting$borrowed_owner_entries, 0, TRUE, outcome_recipe = a$recipe,
    outcome_pooling = "separate")
  testthat::expect_gt(n, before)
  for (method in c("reg", "dr")) {
    old <- tryCatch(qbf_fit(a, method, FALSE, qbf_control(1)), error = identity)
    new <- tryCatch(qbf_fit(a, method, TRUE, qbf_control(1)), error = identity)
    testthat::expect_s3_class(old, "cdid_ipw_owner_bound_error")
    testthat::expect_identical(new, old)
    bad <- qbf_case(bad_prediction = TRUE)
    old <- tryCatch(qbf_fit(bad, method), error = identity)
    new <- tryCatch(qbf_fit(bad, method, TRUE), error = identity)
    testthat::expect_s3_class(old, "cdid_or_wls_rank_error")
    testthat::expect_identical(class(new), class(old))
    testthat::expect_identical(conditionMessage(new), conditionMessage(old))
  }
})

testthat::test_that("raw invalid rows are refused before buffering with original error messages", {
  recipe <- qbf("or_pooled_outcome_recipe")(c("(Intercept)", "x"), TRUE, "x")
  start <- qbf("or_pooled_support_start")(recipe, 1)
  Z <- cbind(`(Intercept)` = 1, x = c(0, 1))
  cases <- list(list(Z = Z, dy = c(Inf, 1), weights = c(1, 1)),
    list(Z = Z, dy = c(1, 2), weights = c(0, 1)),
    list(Z = Z, dy = c(1e308, 2), weights = c(1e308, 1)),
    list(Z = cbind(`(Intercept)` = 2, x = c(0, 1)), dy = c(1, 2), weights = c(1, 1)))
  for (args in cases) {
    old <- tryCatch(do.call(qbf("or_pooled_support_add"), c(list(acc = start, recipe = recipe), args,
      list(route = "baseline"))), error = identity)
    new <- tryCatch(do.call(qbf("or_baseline_qr_buffer_moments"), c(list(acc = start, recipe = recipe), args)),
      error = identity)
    testthat::expect_s3_class(old, "error")
    testthat::expect_identical(class(new), class(old))
    testthat::expect_identical(conditionMessage(new), conditionMessage(old))
  }
})


# V22 independent pooled WLS: distinct actual (pair, control set) intercepts,
# one shared scaled slope, all rate weights refitted from raw unit frequency.
qsp_raw_pooled <- function(a, frequency = rep(1, a$N)) {
  w <- a$w * frequency; A <- rowSums(!is.na(a$Y)) >= 2L
  key <- vapply(a$defs, function(d) paste(d$s, d$t, paste(sort(d$K), collapse = ":"), sep = "|"), character(1))
  keys <- unique(key); L <- length(keys)
  columns <- if (a$intercept) seq.int(2L, ncol(a$X)) else seq_len(ncol(a$X))
  if (a$intercept && ncol(a$X) == 1L) columns <- integer(0)
  blocks <- lapply(seq_along(keys), function(b) {
    d <- a$defs[[match(keys[b], key)]]
    obs <- !is.na(a$Y[, d$s]) & !is.na(a$Y[, d$t]); ix <- which(obs & a$g %in% d$K)
    rate <- vapply(d$K, function(k) sum(w[a$g == k & obs]) / sum(w[a$g == k & A]), numeric(1))
    stopifnot(all(rate > 1e-6))
    ww <- w[ix] / rate[match(a$g[ix], d$K)]
    D <- (d$t - d$s) * a$X[ix, columns, drop = FALSE]
    if (a$intercept) {
      dummy <- matrix(0, length(ix), L); dummy[, b] <- 1; D <- cbind(dummy, D)
    }
    list(X = D, dy = a$Y[ix, d$t] - a$Y[ix, d$s], w = ww)
  })
  fit <- stats::lm.wfit(do.call(rbind, lapply(blocks, `[[`, "X")),
    unlist(lapply(blocks, `[[`, "dy"), use.names = FALSE),
    unlist(lapply(blocks, `[[`, "w"), use.names = FALSE), tol = 1e-10)
  beta <- fit$coefficients; beta[is.na(beta)] <- 0
  vapply(seq_along(a$defs), function(j) {
    d <- a$defs[[j]]; ix <- which(a$g == d$g & !is.na(a$Y[, d$s]) & !is.na(a$Y[, d$t]))
    D <- (d$t - d$s) * a$X[ix, columns, drop = FALSE]
    if (a$intercept) {
      dummy <- matrix(0, length(ix), L); dummy[, match(key[j], keys)] <- 1; D <- cbind(dummy, D)
    }
    weighted.mean(a$Y[ix, d$t] - a$Y[ix, d$s] - as.vector(D %*% beta), w[ix])
  }, numeric(1))
}

qsp_design <- function(a) {
  # Existing deterministic covariate z with a changed unit of measurement:
  # its norm forces a nonidentity LAPACK pivot; the outcome/sample stay fixed.
  a$X <- cbind(`(Intercept)` = 1, x = a$X[, "x"], z = 100 * a$V[, "z"])
  reader <- a$readers$reg; panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  frame <- qbf("construct_ipw_frame_owner")(panel, seq_len(a$N), a$X, a$w, a$control)
  a$readers$reg <- qbf("construct_cp04_frame_reader")(frame, seq_len(a$N),
    matrix(1, a$N, 1L, dimnames = list(NULL, "selection_intercept")), a$g, a$w, a$control)
  a$views$reg <- lapply(seq_along(a$defs), function(j) {
    d <- a$defs[[j]]
    qbf("construct_ipw_link_view")(frame, d$s, d$t, d$g, d$K, control = a$control,
      link_id = paste0("qr-buffer-", j), functional_key = paste0("cmpcat:qr-buffer-", j))
  })
  a$recipe <- qbf("or_pooled_outcome_recipe")(colnames(a$X), TRUE, c("x", "z"))
  a
}

qsp_native <- function(fit, plan, rows, data = plan$data, cap = NULL) {
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE); st <- reader$state
  if (is.null(cap)) cap <- qbf("or_sparse_native_live")(length(rows), length(plan$links),
    length(plan$dates), length(plan$sids), ncol(data$block), length(data$beta), max(c(0L, data$block[3L, ])))
  .Call(qbf("or_baseline_native_kernel")("C_gmm_or_sparse_if"), panel$observations$y,
    panel$postings$post_a_index, panel$postings$post_obs_index, as.integer(rows), st$cohort_index,
    st$a_index, st$w, plan$bounds, frame$Z_a, data, as.double(cap))
}

testthat::test_that("pooled incomplete QR buffers retain the raw joint objective, aliases and global identification", {
  for (options in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE),
      list(bad_prediction = TRUE))) {
    a <- do.call(qbf_case, options)
    old <- qbf_fit(a, "reg", mode = "pooled"); new <- qbf_fit(a, "reg", TRUE, mode = "pooled")
    testthat::expect_null(new$complete_support)
    testthat::expect_gt(new$accounting$fit_scratch_bound, old$accounting$fit_scratch_bound)
    testthat::expect_identical(new$control, old$control)
    testthat::expect_identical(new$shared$keep, old$shared$keep)
    testthat::expect_equal(new$shared$beta, old$shared$beta, tolerance = 1e-9)
    testthat::expect_equal(vapply(new$links, `[[`, numeric(1), "mu"), qsp_raw_pooled(a), tolerance = 1e-9)
    left <- right <- matrix(0, a$N, length(new$links))
    for (first in seq.int(1L, a$N, by = 3L)) {
      rows <- seq.int(first, min(a$N, first + 2L))
      left[rows, ] <- qbf_eval(new, "reg", rows, seq_along(new$links))$phi
      right[rows, ] <- qbf_eval(old, "reg", rows, seq_along(old$links))$phi
    }
    testthat::expect_equal(left, right, tolerance = 1e-9)
    testthat::expect_equal(crossprod(left) / a$N, crossprod(right) / a$N, tolerance = 1e-9)
    if (isTRUE(options$bad_prediction)) {
      testthat::expect_true(all(a$X[a$g == Inf, "x"] == 0))
      testthat::expect_identical(length(new$shared$keep), 1L)
    }
  }
  a <- qsp_design(qbf_case())
  old <- qbf_fit(a, "reg", mode = "pooled"); new <- qbf_fit(a, "reg", TRUE, mode = "pooled")
  testthat::expect_false(identical(new$shared$fit$slope_fit$pivot, seq_len(new$shared$p)))
  testthat::expect_equal(new$shared$beta, old$shared$beta, tolerance = 1e-9)
  # The underlying mode-neutral buffer also admits a tiny, non-owner-sized flush.
  s <- old$support_specs[[1L]]
  small <- qbf("or_baseline_qr_buffered_support")(a$readers$reg, s, a$recipe, a$control, 2L)
  ordinary <- qbf("or_baseline_support")(a$readers$reg, s, a$recipe, a$control)
  for (theta in list(c(.1, -.2, .03), c(-.7, .1, -.04))) {
    F <- small$design$summary$raw$factor; G <- ordinary$design$summary$raw$factor
    testthat::expect_equal(sum(as.vector(F %*% c(-theta, 1))^2),
      sum(as.vector(G %*% c(-theta, 1))^2), tolerance = 1e-9)
  }
  # Finite old admission is preserved if the optional pooled buffer is unfunded.
  cap <- ceiling(old$accounting$retained_bound + old$accounting$fit_scratch_bound)
  plain <- qbf_fit(a, "reg", FALSE, qbf_control(cap), mode = "pooled")
  tight <- qbf_fit(a, "reg", TRUE, qbf_control(cap), mode = "pooled")
  testthat::expect_identical(tight$accounting$fit_scratch_bound, plain$accounting$fit_scratch_bound)
  testthat::expect_identical(tight$shared$beta, plain$shared$beta)
})

testthat::test_that("sparse native pooled and separate OR retain every raw frequency nuisance path", {
  a <- qbf_case(); cc <- qbf_control(chunk = a$N)
  for (mode in c("pooled", "separate")) {
    fit <- qbf_fit(a, "reg", control = cc, mode = mode); ll <- seq_along(fit$links)
    plan <- qbf("chained_full_if_plan")(fit, ll, 37L)
    testthat::expect_false(is.null(plan))
    testthat::expect_true(qbf("ipw_owner_seal_ok")(plan, "or_sparse_native_plan"))
    testthat::expect_lte(plan$entries, plan$retained_bound)
    pre <- qbf("or_baseline_if_budget")(fit, 37L, ll)
    testthat::expect_lte(plan$reservation, qbf("or_sparse_native_slice")(fit, pre, plan$min_B, length(ll)))
    actual <- scalar <- matrix(0, a$N, length(ll))
    for (rows in list(1:37, 38:74)) {
      z <- qbf("evaluate_or_baseline_if_chunk")(fit, rows, ll, plan = plan)
      testthat::expect_identical(z$diagnostics$execution, "native_sparse_fixed_operator")
      direct <- qsp_native(fit, plan, rows)
      testthat::expect_identical(direct$status[1L], 0L)
      testthat::expect_equal(unname(z$phi), direct$phi, tolerance = 1e-10)
      actual[rows, ] <- z$phi
      scalar[rows, ] <- qbf_eval(fit, "reg", rows, ll)$phi
      testthat::expect_equal(z$phi,
        qbf("evaluate_or_baseline_if_chunk")(fit, rows, ll, details = TRUE)$phi, tolerance = 1e-10)
    }
    testthat::expect_equal(actual, scalar, tolerance = 1e-10)
    raw <- if (mode == "pooled") function(freq) qsp_raw_pooled(a, freq) else function(freq) qbf_raw(a, "reg", freq)
    testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), raw(rep(1, a$N)), tolerance = 1e-9)
    numeric_if <- t(vapply(seq_len(a$N), function(i) {
      direction <- rep(-1, a$N); direction[i] <- a$N - 1
      (raw(1 + 1e-6 * direction) - raw(1 - 1e-6 * direction)) / 2e-6
    }, numeric(length(ll))))
    testthat::expect_equal(actual, numeric_if, tolerance = 2e-6)
    testthat::expect_equal(crossprod(actual) / a$N, crossprod(scalar) / a$N, tolerance = 1e-10)
    # The row set includes one zero-outcome and one singleton original F unit.
    testthat::expect_identical(rowSums(!is.na(a$Y))[73:74], c(0, 1))
    testthat::expect_equal(actual[73:74, ], scalar[73:74, ], tolerance = 1e-10)
    rows <- rev(seq_len(37L)); sub <- c(6L, 3L, 1L)
    sp <- qbf("chained_full_if_plan")(fit, sub, 37L)
    testthat::expect_false(is.null(sp))
    testthat::expect_equal(unname(qbf("evaluate_or_baseline_if_chunk")(fit, rows, sub, plan = sp)$phi),
      scalar[rows, sub, drop = FALSE], tolerance = 1e-10)
    testthat::expect_equal(qbf("evaluate_or_baseline_if_chunk")(fit, rows, sub, "core", plan = sp)$phi,
      a$readers$reg$state$lambda * qbf("evaluate_or_baseline_if_chunk")(fit, rows, sub, plan = sp)$phi,
      tolerance = 1e-10)
    if (mode == "pooled") {
      omitted <- setdiff(sp$beta_gids, sp$gids)
      testthat::expect_gt(length(omitted), 0L)
      all_gids <- sort(union(sp$gids, sp$beta_gids)); perturbed <- sp$data
      perturbed$block[6L, match(omitted, all_gids)] <- 0L
      testthat::expect_gt(max(abs(qsp_native(fit, sp, 1:37, perturbed)$phi - scalar[1:37, sub])), 1e-6)
    }
    # Both coefficient-rate forcing paths are needed on these actual masks.
    no_rate <- plan$data; no_rate$forcing_alpha[] <- 0; no_rate$forcing_beta[] <- 0
    testthat::expect_gt(max(abs(qsp_native(fit, plan, 1:37, no_rate)$phi - actual[1:37, ])), 1e-6)
    # Independent linear input-channel check: centers act on ALL F rows.
    centered <- plan$data
    center_delta <- seq_along(centered$center) / 100
    alpha_delta <- seq_len(ncol(centered$block)) / 1000
    centered$center <- centered$center + center_delta
    centered$block_par[4L, ] <- centered$block_par[4L, ] + alpha_delta
    expected <- actual
    for (j in ll) {
      b <- centered$link[1L, j]; oi <- centered$block[1L, b]
      p <- centered$owner[2L, oi]
      ix <- seq.int(centered$owner[1L, oi], length.out = p)
      ai <- seq.int(centered$link[3L, j], length.out = p)
      expected[, j] <- expected[, j] + sum(center_delta[ix] * centered$adj[ai]) +
        alpha_delta[b] / centered$block_par[3L, b]
    }
    for (rows in list(1:37, 38:74)) testthat::expect_equal(
      qsp_native(fit, plan, rows, centered)$phi, expected[rows, , drop = FALSE], tolerance = 1e-10)
    testthat::expect_gt(max(abs(expected[73:74, ] - actual[73:74, ])), 1e-6)
  }
})

testthat::test_that("sparse OR native bases, declined tails and source integrity retain old contracts", {
  for (options in list(list(slopes = FALSE), list(intercept = FALSE), list(alias = TRUE))) {
    a <- do.call(qbf_case, options)
    for (mode in c("pooled", "separate")) {
      fit <- qbf_fit(a, "reg", control = qbf_control(chunk = a$N), mode = mode)
      ll <- seq_along(fit$links); plan <- qbf("chained_full_if_plan")(fit, ll, 37L)
      testthat::expect_false(is.null(plan))
      testthat::expect_equal(qbf("evaluate_or_baseline_if_chunk")(fit, 38:74, ll, plan = plan)$phi,
        qbf_eval(fit, "reg", 38:74, ll)$phi, tolerance = 1e-10)
      # This one-row tail cannot fund descriptor lifetime; exact old stream runs.
      testthat::expect_null(qbf("chained_full_if_plan")(fit, ll, 73L))
      testthat::expect_equal(qbf("evaluate_or_baseline_if_chunk")(fit, 1:73, ll, plan = NULL)$phi,
        qbf_eval(fit, "reg", 1:73, ll)$phi, tolerance = 0)
      changed <- plan; changed$kind <- "wrong-kind"
      testthat::expect_error(qbf("evaluate_or_baseline_if_chunk")(fit, 1:37, ll, plan = changed),
        class = "cdid_or_baseline_integrity_error")
      changed <- plan; changed$data$link[1L, 1L] <- 0L
      testthat::expect_error(qbf("evaluate_or_baseline_if_chunk")(fit, 1:37, ll, plan = changed),
        class = "cdid_or_baseline_integrity_error")
      testthat::expect_error(qbf("evaluate_or_baseline_if_chunk")(fit, 1:37, rev(ll), plan = plan),
        class = "cdid_or_baseline_integrity_error")
      testthat::expect_error(qbf("evaluate_or_baseline_if_chunk")(fit, 1:36, ll, plan = plan),
        class = "cdid_or_baseline_integrity_error")
      broken <- fit; broken$links[[1L]]$mu <- broken$links[[1L]]$mu + .1
      testthat::expect_error(qbf("evaluate_or_baseline_if_chunk")(broken, 1:37, ll, plan = plan),
        class = "cdid_or_baseline_integrity_error")
      invalid <- plan$data; invalid$link[2L, 1L] <- as.integer(length(plan$sids) + 1L)
      testthat::expect_error(qsp_native(fit, plan, 1:37, invalid), class = "cdid_gmm_bound_error")
      testthat::expect_error(qsp_native(fit, plan, 1:37, cap = 37 * length(ll) - 1),
        class = "cdid_gmm_bound_error")
    }
  }
})


# V27 raw full pooled DR: one joint outcome regression and a separate actual
# treatment likelihood for each input comparison, all weights/rates refitted.
dsp27_raw <- function(a, frequency = rep(1, a$N)) {
  w <- a$w * frequency; A <- rowSums(!is.na(a$Y)) >= 2L
  keys <- vapply(a$defs, function(d) paste(d$s, d$t, paste(sort(d$K), collapse = ":"), sep = "|"), character(1))
  unique_keys <- unique(keys); ng <- length(unique_keys)
  slope <- if (a$intercept && ncol(a$X) == 1L) integer(0) else
    if (a$intercept) seq.int(2L, ncol(a$X)) else seq_len(ncol(a$X))
  design <- function(d, ix, b) {
    X <- (d$t - d$s) * a$X[ix, slope, drop = FALSE]
    if (a$intercept) {
      dummy <- matrix(0, length(ix), ng); dummy[, b] <- 1; X <- cbind(dummy, X)
    }
    X
  }
  blocks <- lapply(seq_along(unique_keys), function(b) {
    d <- a$defs[[match(unique_keys[b], keys)]]
    obs <- !is.na(a$Y[, d$s]) & !is.na(a$Y[, d$t]); ix <- which(obs & a$g %in% d$K)
    rates <- vapply(d$K, function(k) sum(w[a$g == k & obs]) / sum(w[a$g == k & A]), numeric(1))
    stopifnot(all(rates > 1e-6))
    list(X = design(d, ix, b), y = a$Y[ix, d$t] - a$Y[ix, d$s],
      w = w[ix] / rates[match(a$g[ix], d$K)])
  })
  outcome <- stats::lm.wfit(do.call(rbind, lapply(blocks, `[[`, "X")),
    unlist(lapply(blocks, `[[`, "y"), use.names = FALSE),
    unlist(lapply(blocks, `[[`, "w"), use.names = FALSE), tol = 1e-10)
  beta <- outcome$coefficients; beta[is.na(beta)] <- 0
  vapply(seq_along(a$defs), function(j) {
    d <- a$defs[[j]]; cohorts <- c(d$g, d$K)
    obs <- !is.na(a$Y[, d$s]) & !is.na(a$Y[, d$t]); ix <- which(obs & a$g %in% cohorts)
    rates <- vapply(cohorts, function(k) sum(w[a$g == k & obs]) / sum(w[a$g == k & A]), numeric(1))
    ww <- w[ix] / rates[match(a$g[ix], cohorts)]; D <- as.double(a$g[ix] == d$g)
    dy <- a$Y[ix, d$t] - a$Y[ix, d$s]
    residual <- dy - as.vector(design(d, ix, match(keys[j], unique_keys)) %*% beta)
    likelihood <- suppressWarnings(stats::glm.fit(a$V[ix, , drop = FALSE], D, weights = ww,
      start = rep(0, ncol(a$V)), family = stats::binomial(), intercept = FALSE,
      control = stats::glm.control(epsilon = 1e-13, maxit = 60L)))
    stopifnot(likelihood$converged, all(is.finite(likelihood$coefficients)))
    eta <- as.vector(a$V[ix, , drop = FALSE] %*% likelihood$coefficients)
    prob <- stats::plogis(eta)
    stopifnot(all(prob > 1e-6 & prob < 1 - 1e-6))
    weighted.mean(residual, ww * D) - weighted.mean(residual, ww * (1 - D) * prob / (1 - prob))
  }, numeric(1))
}

dsp27_prepare <- function(fit, plan, rows, data = plan$data, outcomes = NULL, cap = NULL) {
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE); st <- reader$state
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(outcomes)) outcomes <- panel$observations$y
  if (is.null(cap)) cap <- qbf("dr_sparse_native_inventory")(length(rows), length(plan$sids),
    length(plan$dates), length(plan$tids), length(plan$links), fit$shared$p,
    fit$treatment_recipe$p, plan$metadata)$prepare
  .Call(qbf("or_baseline_native_kernel")("C_gmm_dr_sparse_prepare"), outcomes,
    panel$postings$post_a_index, panel$postings$post_obs_index, as.integer(rows), st$cohort_index,
    st$a_index, st$w, plan$bounds, frame$Z_a, data, as.double(cap))
}
dsp27_solve <- function(fit, plan, prepared) {
  p <- fit$shared$p; B <- nrow(prepared$up)
  xi_beta <- if (p) t(if (fit$outcome_recipe$intercept)
    qbf("or_pooled_wls_gram_solve")(fit$shared$fit, t(prepared$or_score), normalizer = plan$N) else
    qbf("or_wls_gram_solve")(fit$shared$fit, t(prepared$or_score), normalizer = plan$N)) else matrix(0, B, 0L)
  n_ok <- if (prepared$status[1L] == 0L) length(plan$tids) else prepared$status[3L] - 1L
  xi <- lapply(seq_len(n_ok), function(j) {
    ix <- (j - 1L) * fit$treatment_recipe$p + seq_len(fit$treatment_recipe$p)
    t(qbf("cp04_solve")(fit$treatments[[plan$tids[j]]]$H_chol, t(prepared$rhs[, ix, drop = FALSE])))
  })
  list(beta = xi_beta, gamma = xi, n_ok = n_ok)
}
dsp27_assemble <- function(fit, plan, prepared, solved, data = plan$data, cap = NULL) {
  if (is.null(cap)) cap <- qbf("dr_sparse_native_inventory")(nrow(prepared$up), length(plan$sids),
    length(plan$dates), length(plan$tids), length(plan$links), fit$shared$p,
    fit$treatment_recipe$p, plan$metadata)$assembler
  .Call(qbf("or_baseline_native_kernel")("C_gmm_dr_sparse_assemble"), prepared$up,
    solved$beta, prepared$direct, solved$gamma, data, as.integer(solved$n_ok), as.double(cap))
}

testthat::test_that("automatic pooled DR QR buffers preserve joint outcome and treatment fits and raw frequency IF", {
  for (options in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE),
      list(bad_prediction = TRUE))) {
    a <- do.call(qbf_case, options)
    old <- qbf_fit(a, "dr", FALSE, mode = "pooled"); new <- qbf_fit(a, "dr", TRUE, mode = "pooled")
    testthat::expect_null(new$complete_support)
    testthat::expect_gt(new$accounting$fit_scratch_bound, old$accounting$fit_scratch_bound)
    testthat::expect_identical(new$control, old$control)
    testthat::expect_identical(new$shared$keep, old$shared$keep)
    testthat::expect_equal(new$shared$beta, old$shared$beta, tolerance = 1e-9)
    testthat::expect_equal(lapply(new$treatments, `[[`, "coefficients"),
      lapply(old$treatments, `[[`, "coefficients"), tolerance = 1e-9)
    testthat::expect_equal(vapply(new$links, `[[`, numeric(1), "mu"), dsp27_raw(a), tolerance = 1e-9)
    H <- G <- matrix(0, a$N, length(new$links))
    for (first in seq.int(1L, a$N, by = 3L)) {
      rows <- seq.int(first, min(a$N, first + 2L))
      H[rows, ] <- qbf_eval(new, "dr", rows, seq_along(new$links))$phi
      G[rows, ] <- qbf_eval(old, "dr", rows, seq_along(old$links))$phi
    }
    testthat::expect_equal(H, G, tolerance = 1e-9)
    testthat::expect_equal(crossprod(H) / a$N, crossprod(G) / a$N, tolerance = 1e-9)
  }
  a <- qbf_case(); fit <- qbf_fit(a, "dr", TRUE, mode = "pooled")
  actual <- matrix(0, a$N, length(fit$links))
  for (first in seq.int(1L, a$N, by = 3L)) {
    rows <- seq.int(first, min(a$N, first + 2L)); actual[rows, ] <- qbf_eval(fit, "dr", rows, seq_along(fit$links))$phi
  }
  numerical <- t(vapply(seq_len(a$N), function(i) {
    direction <- rep(-1, a$N); direction[i] <- a$N - 1
    (dsp27_raw(a, 1 + 1e-6 * direction) - dsp27_raw(a, 1 - 1e-6 * direction)) / 2e-6
  }, numeric(length(fit$links))))
  testthat::expect_equal(actual, numerical, tolerance = 2e-6)
  # Unfunded automatic fit buffer retains the old admitted fit, not a new refusal.
  cap <- ceiling(fit$accounting$retained_bound + fit$accounting$fit_scratch_bound)
  plain <- qbf_fit(a, "dr", FALSE, qbf_control(cap), mode = "pooled")
  cap <- ceiling(plain$accounting$retained_bound + plain$accounting$fit_scratch_bound)
  plain <- qbf_fit(a, "dr", FALSE, qbf_control(cap), mode = "pooled")
  tight <- qbf_fit(a, "dr", TRUE, qbf_control(cap), mode = "pooled")
  testthat::expect_identical(tight$accounting$fit_scratch_bound, plain$accounting$fit_scratch_bound)
  testthat::expect_identical(lapply(tight$treatments, `[[`, "coefficients"), lapply(plain$treatments, `[[`, "coefficients"))
})

testthat::test_that("native incomplete pooled DR keeps full gamma, global outcome, rate paths and full F rows", {
  for (options in list(list(), list(intercept = FALSE), list(slopes = FALSE), list(alias = TRUE))) {
    a <- do.call(qbf_case, options)
    fit <- qbf_fit(a, "dr", control = qbf_control(chunk = a$N), mode = "pooled")
    ll <- seq_along(fit$links); plan <- qbf("chained_full_if_plan")(fit, ll, 37L)
    testthat::expect_false(is.null(plan))
    testthat::expect_true(qbf("ipw_owner_seal_ok")(plan, "dr_sparse_native_plan"))
    testthat::expect_lte(plan$entries, plan$retained_bound)
    scalar <- actual <- matrix(0, a$N, length(ll))
    for (rows in list(1:37, 38:74)) {
      out <- qbf("evaluate_dr_baseline_if_chunk")(fit, rows, ll, plan = plan)
      testthat::expect_identical(out$diagnostics$execution, "native_sparse_pooled_dr")
      actual[rows, ] <- out$phi; scalar[rows, ] <- qbf_eval(fit, "dr", rows, ll)$phi
      prepared <- dsp27_prepare(fit, plan, rows); testthat::expect_identical(prepared$status[1L], 0L)
      solved <- dsp27_solve(fit, plan, prepared)
      direct <- dsp27_assemble(fit, plan, prepared, solved)
      testthat::expect_identical(direct$status[1L], 0L)
      testthat::expect_equal(unname(out$phi), direct$phi, tolerance = 1e-10)
      detailed <- qbf("evaluate_dr_baseline_if_chunk")(fit, rows, ll, details = TRUE)
      testthat::expect_equal(out$phi, detailed$phi, tolerance = 1e-10)
      for (j in seq_along(plan$tids)) testthat::expect_equal(unname(solved$gamma[[j]]),
        unname(detailed$details$models[[plan$tids[j]]]$gamma), tolerance = 1e-10)
    }
    testthat::expect_equal(actual, scalar, tolerance = 1e-10)
    testthat::expect_equal(crossprod(actual) / a$N, crossprod(scalar) / a$N, tolerance = 1e-10)
    testthat::expect_equal(actual[73:74, ], scalar[73:74, ], tolerance = 1e-10)
    rows <- rev(1:37); sub <- c(6L, 3L, 1L)
    sp <- qbf("chained_full_if_plan")(fit, sub, 37L)
    testthat::expect_false(is.null(sp))
    testthat::expect_equal(unname(qbf("evaluate_dr_baseline_if_chunk")(fit, rows, sub, plan = sp)$phi),
      scalar[rows, sub, drop = FALSE], tolerance = 1e-10)
    testthat::expect_equal(qbf("evaluate_dr_baseline_if_chunk")(fit, rows, sub, "core", plan = sp)$phi,
      a$readers$dr$state$lambda * qbf("evaluate_dr_baseline_if_chunk")(fit, rows, sub, plan = sp)$phi,
      tolerance = 1e-10)
    # This descriptor is funded even through its one-row tail; require use.
    tail <- qbf("chained_full_if_plan")(fit, ll, 73L)
    testthat::expect_false(is.null(tail)); testthat::expect_identical(tail$min_B, 1L)
    testthat::expect_identical(qbf("evaluate_dr_baseline_if_chunk")(fit, 74L, ll, plan = tail)$diagnostics$execution,
      "native_sparse_pooled_dr")
    testthat::expect_equal(unname(qbf("evaluate_dr_baseline_if_chunk")(fit, 74L, ll, plan = tail)$phi),
      scalar[74L, , drop = FALSE], tolerance = 1e-10)
    if (fit$shared$p) {
      testthat::expect_gt(length(setdiff(sp$beta_gids, sp$gids)), 0L)
      prepared <- dsp27_prepare(fit, sp, 1:37); changed <- sp$data
      omitted <- match(setdiff(sp$beta_gids, sp$gids), sp$beta_gids)
      # Leave valid nonempty spans but suppress an omitted block's OR slope
      # weight by zeroing its centered design contribution through L/mean.
      changed$block_par[, omitted] <- 0
      changed$block_par[1L, omitted] <- 1e-12
      wrong <- dsp27_prepare(fit, sp, 1:37, changed)
      wrong_solved <- dsp27_solve(fit, sp, wrong)
      testthat::expect_gt(max(abs(dsp27_assemble(fit, sp, wrong, wrong_solved)$phi - scalar[1:37, sub])), 1e-6)
      # Actual fitted rate forcing must be retained in the global score.
      changed <- plan$data; changed$forcing_beta[] <- 0
      left <- dsp27_prepare(fit, plan, 1:37, changed)
      testthat::expect_gt(max(abs(dsp27_assemble(fit, plan, left, dsp27_solve(fit, plan, left))$phi - actual[1:37, ])), 1e-6)
    }
  }
})

testthat::test_that("sparse pooled DR validates all four frame coordinates and full gamma with staged failures", {
  a <- qbf_case(slopes = FALSE)
  fit <- qbf_fit(a, "dr", control = qbf_control(chunk = a$N), mode = "pooled")
  ll <- seq_along(fit$links); plan <- qbf("chained_full_if_plan")(fit, ll, 37L)
  testthat::expect_false(is.null(plan))
  prepared <- dsp27_prepare(fit, plan, 38:74); solved <- dsp27_solve(fit, plan, prepared)
  ordinary <- dsp27_assemble(fit, plan, prepared, solved); testthat::expect_identical(ordinary$status[1L], 0L)
  # A finite projected psi is insufficient: every coordinate must be finite,
  # including dp/dm on an unobserved F row. Inject only direct-coordinate input.
  row <- 37L
  for (cc in seq_len(4L)) {
    changed <- prepared
    changed$direct[row, cc] <- Inf
    z <- dsp27_assemble(fit, plan, changed, solved)
    testthat::expect_identical(z$status[1L], 7L)
    testthat::expect_identical(z$status[3L], 1L)
  }
  # Even if every loading on a gamma direction is zero, its nonfinite value
  # itself is refused; the C boundary checks the complete root IF.
  changed <- solved; changed$gamma[[1L]][row, 1L] <- Inf
  data <- plan$data; data$partial_beta[seq.int(1L, 4L), 1L] <- 0
  testthat::expect_identical(dsp27_assemble(fit, plan, prepared, changed, data)$status[1L], 7L)
  # The first model's failed coordinate outranks a later model endpoint failure.
  reader <- a$readers$dr; panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  first <- fit$treatment_specs[[plan$tids[1L]]]$pair
  later <- which(vapply(fit$treatment_specs[plan$tids], function(s)
    !s$pair$t_calendar_index %in% c(first$s_calendar_index, first$t_calendar_index), logical(1)))[1L]
  testthat::expect_true(is.finite(later))
  date <- fit$treatment_specs[[plan$tids[later]]]$pair$t_calendar_index
  y <- panel$observations$y
  position <- panel$postings$post_obs_index[panel$postings$date_start[date]]
  y[position] <- Inf
  bad <- dsp27_prepare(fit, plan, 1:37, outcomes = y)
  testthat::expect_identical(bad$status[2L], 2L)
  testthat::expect_gt(bad$status[3L], 1L)
  prefix <- dsp27_solve(fit, plan, bad)
  bad$direct[1L, 3L] <- Inf
  testthat::expect_identical(dsp27_assemble(fit, plan, bad, prefix)$status[3L], 1L)
  # Exact source/link/partition checks and genuine original-cap refusals remain.
  changed <- plan; changed$kind <- "wrong-kind"
  testthat::expect_error(qbf("evaluate_dr_baseline_if_chunk")(fit, 1:37, ll, plan = changed),
    class = "cdid_dr_baseline_integrity_error")
  testthat::expect_error(qbf("evaluate_dr_baseline_if_chunk")(fit, 1:36, ll, plan = plan),
    class = "cdid_dr_baseline_integrity_error")
  testthat::expect_error(qbf("evaluate_dr_baseline_if_chunk")(fit, 1:37, rev(ll), plan = plan),
    class = "cdid_dr_baseline_integrity_error")
  broken <- fit; broken$treatments[[1L]]$coefficients[1L] <- broken$treatments[[1L]]$coefficients[1L] + .1
  testthat::expect_error(qbf("evaluate_dr_baseline_if_chunk")(broken, 1:37, ll, plan = plan),
    class = "cdid_dr_baseline_integrity_error")
  invalid <- plan$data; invalid$model_support[1L] <- as.integer(length(plan$sids) + 1L)
  testthat::expect_error(dsp27_prepare(fit, plan, 1:37, invalid), class = "cdid_gmm_bound_error")
  testthat::expect_error(dsp27_prepare(fit, plan, 1:37, cap = 37 * length(plan$sids) - 1),
    class = "cdid_gmm_bound_error")
  pre <- qbf("dr_baseline_if_budget")(fit, 37L, ll)
  limited <- fit; limited$ref <- NULL
  limited$control$max_matrix_entries <- limited$accounting$owned_entries + limited$accounting$borrowed_entries + pre$scratch - 1
  limited <- qbf("ipw_owner_seal")(limited, "dr_baseline_backend")
  testthat::expect_error(qbf("dr_sparse_native_plan")(limited, ll, 37L), class = "cdid_ipw_owner_bound_error")
  # The private replay status must keep the scalar numeric condition, rather
  # than reusing the unrelated complete-support status9 identity refusal.
  replay <- plan$data; replay$trt_beta[, 1L] <- .Machine$double.xmax
  invalid_replay <- dsp27_prepare(fit, plan, 1:37, replay)
  testthat::expect_identical(invalid_replay$status[1L], 27L)
  condition <- tryCatch(qbf("dr_sparse_native_refuse")(invalid_replay$status), error = identity)
  testthat::expect_s3_class(condition, "cdid_cp04_numeric_error")
  testthat::expect_identical(conditionMessage(condition), "Nonfinite baseline DR treatment replay.")
  # Independent native-input clipping check. These are deliberately supplied
  # operator values, not a claim that the baseline nuisance fit binds the cap.
  clipped <- plan$data; clipped$trt_beta[, 1L] <- c(52, rep(0, fit$treatment_recipe$p - 1L))
  z <- dsp27_prepare(fit, plan, 1:37, clipped)
  testthat::expect_identical(z$status[1L], 0L)
  spec <- fit$treatment_specs[[plan$tids[1L]]]; f <- fit$treatments[[plan$tids[1L]]]
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  rows <- 1:37; ci <- st$cohort_index[rows]
  source_cohorts <- vapply(fit$support_specs[spec$support_ids], `[[`, integer(1), "cohort")
  obs <- which(ci %in% source_cohorts & !is.na(a$Y[rows, spec$pair$s_calendar_index]) &
    !is.na(a$Y[rows, spec$pair$t_calendar_index]))
  X <- frame$Z_a[st$a_index[rows[obs]], fit$treatment_recipe$index, drop = FALSE]
  D <- as.double(ci[obs] == source_cohorts[1L])
  prob <- stats::plogis(rep(52, length(obs)))
  testthat::expect_true(any(D == 0)); testthat::expect_true(all(prob > 1 - 1e-6))
  fw <- st$w[rows[obs]] * vapply(ci[obs], function(k)
    fit$supports[[spec$support_ids[match(k, source_cohorts)]]]$v, numeric(1))
  score <- matrix(0, length(rows), fit$treatment_recipe$p)
  score[obs, ] <- fw * (D - prob) * X
  upstream <- z$up[, match(spec$support_ids, plan$sids), drop = FALSE]
  testthat::expect_equal(z$rhs[, seq_len(fit$treatment_recipe$p), drop = FALSE],
    score + upstream %*% t(f$cross), tolerance = 1e-10)
  wp <- fw * D; wm <- fw * (1 - D) * ((1 - 1e-6) / (1 - (1 - 1e-6)))
  residual <- a$Y[rows[obs], spec$pair$t_calendar_index] - a$Y[rows[obs], spec$pair$s_calendar_index] -
    clipped$model_par[2L, 1L]
  direct <- matrix(0, length(rows), 4L)
  direct[, 3L] <- -f$means["dp"]; direct[, 4L] <- -f$means["dm"]
  direct[obs, 1L] <- wp * (residual - f$means["T_dr"]) / f$means["dp"]
  direct[obs, 2L] <- wm * (residual - f$means["M_dr"]) / f$means["dm"]
  direct[obs, 3L] <- direct[obs, 3L] + wp; direct[obs, 4L] <- direct[obs, 4L] + wm
  testthat::expect_equal(z$direct[, 1:4, drop = FALSE], direct, tolerance = 1e-10)
  # Optional plan construction declines when its real inventoried lifetime is
  # one entry above the available old slice, leaving the scalar evaluator.
  original <- qbf("dr_sparse_native_inventory")
  testthat::local_mocked_bindings(dr_sparse_native_inventory = function(...) {
    z <- original(...); z$peak <- pre$scratch + 1; z
  }, .package = "cdid")
  testthat::expect_null(qbf("dr_sparse_native_plan")(fit, ll, 37L))
  testthat::expect_equal(qbf_eval(fit, "dr", 1:37, ll)$phi,
    qbf("evaluate_dr_baseline_if_chunk")(fit, 1:37, ll, plan = NULL)$phi, tolerance = 0)
})
