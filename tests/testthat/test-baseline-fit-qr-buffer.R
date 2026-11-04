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
    complete <- qbf_case(complete = TRUE)
    testthat::expect_false(is.null(qbf_fit(complete, method, TRUE)$separate_complete))
    qbf_fit(a, method, TRUE, mode = "pooled")
  }
  testthat::expect_identical(n, before)
  # The OR public adapter forwards its existing automatic flag.
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
