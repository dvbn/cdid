# Independent D93 checkpoint reference. Expected values below use only explicit
# finite-support masses, saturated raw roots and implicit differentiation of
# their score equations. The raw-H formula is checked separately at the returned
# coefficient, supplied only as an input. No cdid helper generates an expectation.
# This compact check covers one corrected IPW link, not the whole clipping section. The
# separate author tests must cover OR/DR, imputation, augmentation and dispatch.
# Reference IF tolerance 1e-8 follows the accepted case-weight convention;
# neither a production solver setting nor an existing test is changed.
d93i_fn <- function(name) getFromNamespace(name, "cdid")

d93i_fixture <- function() {
  d <- expand.grid(r1 = 0:1, r2 = 0:1, r3 = 0:1, x = 0:1, k = 1:2)
  p <- c(.3, 4e-4)[d$x + 1L]
  p3 <- c(.4, .7)[d$x + 1L]
  mass <- matrix(c(.7, 1.2, 1.5, 2e-7), 2L, 2L)
  w <- mass[cbind(d$k, d$x + 1L)] *
    ifelse(d$r1, p, 1 - p) * ifelse(d$r2, p, 1 - p) *
    ifelse(d$r3, p3, 1 - p3)
  w <- nrow(d) * w / sum(w)
  obs <- as.matrix(d[c("r1", "r2", "r3")])
  storage.mode(obs) <- "double"
  dy <- 1 + d$k + 2 * d$x + .7 * d$r3 + .6 * d$r1 * d$r3 + .3 * d$k * d$x
  Y <- cbind(0, dy, 2 * dy + d$r1)
  Y[obs == 0] <- NA_real_
  X <- cbind("(Intercept)" = 1, x = as.double(d$x))
  list(N = nrow(d), k = d$k, x = d$x, w = w, obs = obs, dy = dy,
    A = as.double(rowSums(obs) >= 2), X = X, Y = Y,
    g = ifelse(d$k == 1L, 30, Inf), rho = rep(1, nrow(d)),
    calendar = c(10, 20, 50), cluster = rep(1:8, 4L),
    links = list(list(s = 10, t = 20, g = 30, K = Inf, link_id = "independent")))
}

d93i_reference <- function(f, model, pooling, universe = "F", mutant = "", V = f$X) {
  N <- f$N; w <- f$w; X <- V; k <- f$k; x <- f$x
  U <- switch(universe, F = rep(TRUE, N), B = rowSums(f$obs) > 0, A = f$A == 1)
  S <- f$obs[, 1L] * f$obs[, 2L]
  a <- vapply(1:2, function(j) sum(w[k == j & U] * f$A[k == j & U]) /
    sum(w[k == j & U]), numeric(1))
  xa <- sapply(1:2, function(j) N * w * (k == j & U) * (f$A - a[j]) /
    sum(w[k == j & U]))
  # Explicit cell dummy columns, with slope sharing dictated by the model.
  # Both endpoint scores belong to the same original row before inversion.
  cells <- if (model == "pair_logit") data.frame(k = 1:2, date = 0L) else
    expand.grid(k = 1:2, date = 1:2)
  nc <- nrow(cells)
  group <- switch(pooling, separate = seq_len(nc), within_cohort = cells$k,
    across_cohorts = rep(1L, nc))
  ng <- max(group); ds <- nc + ng
  Hs <- matrix(0, ds, ds); ss <- matrix(0, N, ds)
  pred <- design <- vector("list", nc); slopes <- numeric(nc)
  for (c in seq_len(nc)) {
    z <- matrix(0, N, ds); z[, c] <- 1; z[, nc + group[c]] <- x
    y <- if (model == "pair_logit") S else f$obs[, cells$date[c]]
    on <- k == cells$k[c] & U
    # Cohort masses cancel within x, even after conditioning on B or A:
    # the observation law given x is identical across cohorts. Both dates in
    # this fixture also have the same law. Verify the shared-slope premise.
    px <- vapply(0:1, function(v) sum(w[on & x == v] * y[on & x == v]) /
      sum(w[on & x == v]), numeric(1))
    p <- px[x + 1L]
    slopes[c] <- diff(stats::qlogis(px))
    Hs <- Hs + crossprod(z, (w * on * p * (1 - p)) * z) / N
    ss <- ss + (w * on * (y - p)) * z
    pred[[c]] <- p; design[[c]] <- z
  }
  stopifnot(all(vapply(split(slopes, group), function(s)
    max(abs(s - s[1L])) < 1e-12, logical(1))))
  xs <- ss %*% solve(Hs)
  q <- numeric(N); dlogq <- matrix(0, N, ds)
  for (j in 1:2) {
    rr <- which(k == j); ids <- which(cells$k == j)
    q[rr] <- Reduce(`*`, lapply(pred[ids], `[`, rr))
    for (c in ids) {
      if (mutant == "drop_endpoint" && model != "pair_logit" && cells$date[c] == 2L) next
      dlogq[rr, ] <- dlogq[rr, , drop = FALSE] +
        (1 - pred[[c]][rr]) * design[[c]][rr, , drop = FALSE]
    }
  }
  qc <- pmax(q, 1e-6)
  if (mutant == "endpoint_floor" && model != "pair_logit") {
    for (j in 1:2) {
      rr <- which(k == j); ids <- which(cells$k == j)
      qc[rr] <- Reduce(`*`, lapply(pred[ids], function(p) pmax(p[rr], 1e-6)))
    }
  }
  L <- matrix(0, N, 2L + ds)
  L[cbind(seq_len(N), k)] <- if (mutant == "drop_numerator") 0 else 1 / a[k]
  L[, 2L + seq_len(ds)] <- -dlogq * if (mutant == "unmasked_q") 1 else as.double(q > 1e-6)
  xu <- cbind(xa, xs)
  h <- w * S * a[k] / (f$rho * qc)
  D <- as.double(k == 1L)
  pt <- vapply(0:1, function(v) sum(h[x == v] * D[x == v]) / sum(h[x == v]), numeric(1))
  p <- pt[x + 1L]; pc <- pmin(p, 1 - 1e-6); odds <- pc / (1 - pc)
  act <- if (mutant == "unmasked_odds") rep(1, N) else as.double(p < 1 - 1e-6)
  # Only prediction derivatives are masked. The score and H use raw p.
  score <- (h * (D - if (mutant == "clipped_score") pc else p)) * X
  H <- crossprod(X, (h * p * (1 - p)) * X) / N
  C <- crossprod(score, L) / N
  xb <- (score + xu %*% t(C)) %*% solve(H)
  wp <- h * D; wm <- h * (1 - D) * odds
  dp <- sum(wp) / N; dm <- sum(wm) / N
  T <- sum(wp * f$dy) / sum(wp); M <- sum(wm * f$dy) / sum(wm)
  residual <- cbind(wp * (f$dy - T) / dp, wm * (f$dy - M) / dm, wp, wm)
  Up <- crossprod(residual, L) / N
  V <- matrix(0, 4L, 2L)
  V[2L, ] <- colSums((wm * act * (f$dy - M)) * X) / (N * dm)
  V[4L, ] <- colSums((wm * act) * X) / N
  direct <- sweep(residual, 2L, c(0, 0, dp, dm), "-")
  coord <- direct + xu %*% t(Up) + xb %*% t(V)
  list(means = c(T, M, dp, dm), coord = unname(coord), phi = coord[, 1L] - coord[, 2L],
    H = unname(H), xb = unname(xb), q = q, p = p, h = h, odds = odds, act = act,
    L = L, xu = xu, wp = wp, wm = wm, Up = Up, V = V, S = S)
}

# Match the documented raw-probability arithmetic on the observed fitting rows,
# independently of production H/probability/clipping helpers. These fixtures
# stay in eta <= 40. A finite accepted coefficient need not give the exact
# saturated information, especially after information conditioning of V.
d93i_raw_at_coefficient <- function(f, V, h, coefficients) {
  rows <- which(f$obs[, 1L] != 0 & f$obs[, 2L] != 0)
  X <- V[rows, , drop = FALSE]
  eta <- drop(X %*% coefficients)
  stopifnot(all(is.finite(eta)), all(eta <= 40))
  raw_odds <- exp(eta)
  p <- raw_odds / (1 + raw_odds)
  H <- crossprod(X, (h[rows] * p * (1 - p)) * X) / f$N
  score <- (h[rows] * (as.double(f$k[rows] == 1L) - p)) * X
  list(rows = rows, p = p, H = unname(H), score = score)
}

testthat::test_that("independent raw-root calculus matches the complete clipped IPW coordinates", {
  original <- d93i_fixture()
  for (model in c("pair_logit", "independent_period_logit")) {
    for (pooling in c("separate", "within_cohort", "across_cohorts")) {
      for (universe in if (model == "pair_logit") c("F", "B", "A") else "F") {
        f <- original
        f$rho <- if (universe == "F") rep(1, f$N) else .8 + .1 * f$x
        if (universe == "A") {
          # The old rare control mass made the conditioned selection root
          # almost invisible to a deviance stopping rule. Balance that cell's
          # mass; compensate with known rho on treated x=1 so both corrected
          # arm masses scale together and the treatment cap stays active.
          f$w[f$k == 2L & f$x == 1L] <- 1e6 * f$w[f$k == 2L & f$x == 1L]
          f$w <- f$N * f$w / sum(f$w)
          f$rho[f$k == 1L & f$x == 1L] <- f$rho[f$k == 1L & f$x == 1L] / 1e6
        }
        U <- switch(universe, F = rep(TRUE, f$N), B = rowSums(f$obs) > 0, A = f$A == 1)
        # Condition the treatment design using independent cell information,
        # then mix both cell-root coordinates. This invertible change keeps
        # the model span and avoids ill-conditioned subtraction and comparisons
        # against cancellation zeros; no IF coordinate is discarded.
        base <- d93i_reference(f, model, pooling, universe)
        cell_H <- vapply(0:1, function(v) {
          rr <- f$x == v
          sum(base$h[rr] * base$p[rr] * (1 - base$p[rr])) / f$N
        }, numeric(1))
        v0 <- (1 - f$x) / sqrt(cell_H[1L]); v1 <- f$x / sqrt(cell_H[2L])
        V <- cbind(mix1 = v0 + v1, mix2 = v0 - 2 * v1)
        ref <- d93i_reference(f, model, pooling, universe, V = V)
        got <- d93i_fn("fit_cp04_selection_table")(Z = f$X, V = V, Y = f$Y,
          obs = f$obs, g = f$g, w = f$w, U = U, rho = f$rho, universe = universe,
          calendar = f$calendar, links = f$links, selection_model = model,
          selection_pooling = pooling, has_intercept = TRUE, chunk_rows = 7L,
          glm_control = list(epsilon = 1e-15, maxit = 60L, trace = FALSE))
        tr <- got$treatments[[1L]]
        at_beta <- d93i_raw_at_coefficient(f, V, ref$h, tr$coefficients)
        testthat::expect_equal(unname(tr$means[1:4]), ref$means, tolerance = 1e-10)
        testthat::expect_equal(unname(crossprod(tr$H_chol)), at_beta$H, tolerance = 1e-10)
        # Root accuracy is a separate check against the closed saturated root;
        # the raw score is checked in coefficient units at the IF tolerance.
        testthat::expect_equal(at_beta$p, ref$p[at_beta$rows], tolerance = 1e-10)
        testthat::expect_lt(max(abs((colSums(at_beta$score) / f$N) %*% solve(at_beta$H))), 1e-8)
        testthat::expect_equal(unname(tr$beta_if), ref$xb, tolerance = 1e-8)
        testthat::expect_equal(unname(tr$coordinates_if), ref$coord, tolerance = 1e-8)
        testthat::expect_equal(as.vector(got$phi), ref$phi, tolerance = 1e-8)
        testthat::expect_lt(max(abs(colMeans(ref$coord))), 1e-8)
        testthat::expect_true(any(ref$p > 1 - 1e-6) && any(ref$p < 1 - 1e-6))
        if (universe != "A") testthat::expect_gt(max(abs(ref$phi[f$A == 0])), 1e-6)
        else testthat::expect_true(all(ref$phi[f$A == 0] == 0))
        centered <- sweep(ref$coord, 2L, colMeans(ref$coord))
        got_centered <- sweep(tr$coordinates_if, 2L, colMeans(tr$coordinates_if))
        testthat::expect_equal(unname(crossprod(got_centered)), crossprod(centered), tolerance = 1e-8)
        testthat::expect_equal(unname(crossprod(rowsum(got_centered, f$cluster))),
          crossprod(rowsum(centered, f$cluster)), tolerance = 1e-8)
      }
    }
  }
})

testthat::test_that("active clipping witnesses distinguish the scoped omissions", {
  f <- d93i_fixture()
  r <- d93i_reference(f, "independent_period_logit", "across_cohorts")
  testthat::expect_true(all(r$q[f$x == 1] < 1e-6))
  testthat::expect_true(all(c(.3, 4e-4) > 1e-6))
  testthat::expect_true(all(r$L[f$x == 1, -(1:2)] == 0))
  testthat::expect_true(all(rowSums(abs(r$L[f$x == 1, 1:2])) > 0))
  for (m in c("unmasked_odds", "clipped_score", "unmasked_q", "drop_numerator",
              "endpoint_floor", "drop_endpoint")) {
    bad <- d93i_reference(f, "independent_period_logit", "across_cohorts", mutant = m)
    # The saturated capped stratum has zero beta loading in the scalar link.
    # Raw-score substitution is detected at the raw-root IF, a separate gate.
    gap <- if (m == "clipped_score") max(abs(bad$xb - r$xb)) else max(abs(bad$coord - r$coord))
    testthat::expect_true(gap > 1e-6, info = m)
  }
  # Equality is checked as a strict-mask contract only, never by a smooth FD.
  testthat::expect_identical(d93i_fn("ipw_clip_floor_active")(1e-6), 0)
  testthat::expect_identical(d93i_fn("ipw_clip_odds_active")(1 - 1e-6), 0)
})

testthat::test_that("sealed streamed IF and its centered Gram match the independent rows", {
  f <- d93i_fixture()
  ctrl <- d93i_fn("ipw_owner_control")(max_chunk_rows = 7L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1e9, max_output_entries = 1e7)
  long <- data.frame(id = rep(seq_len(f$N), each = 3L), time = rep(f$calendar, f$N),
    y = as.vector(t(f$Y)), g = rep(f$g, each = 3L))
  frame_only <- tail(which(rowSums(f$obs) == 0), 1L)
  long <- long[long$id != frame_only, ]
  panel <- d93i_fn("construct_ipw_panel_owner")(long, "id", "time", "y", "g", ctrl)
  frame <- d93i_fn("construct_ipw_frame_owner")(panel, seq_len(f$N), f$X,
    weights = f$w, control = ctrl, frame_g = f$g, allow_frame_only = TRUE)
  reader <- d93i_fn("construct_cp04_frame_reader")(frame, seq_len(f$N), f$X, f$g, f$w,
    ctrl, rho = f$rho, has_intercept = TRUE)
  view <- d93i_fn("construct_ipw_link_view")(frame, 10, 20, g = 30, K = Inf,
    control = ctrl, link_id = "independent", functional_key = "cmpcat:independent-D93")
  for (model in c("pair_logit", "independent_period_logit")) {
    r <- d93i_reference(f, model, "across_cohorts")
    b <- d93i_fn("fit_cp04_selection_models")(reader, list(view), ctrl,
      list(epsilon = 1e-13, maxit = 60L, trace = FALSE), selection_model = model,
      selection_pooling = "across_cohorts")
    testthat::expect_true(d93i_fn("ipw_owner_seal_ok")(b, "cp04_backend"))
    testthat::expect_equal(unname(d93i_fn("cp04_payload")(b)["entries"]), b$accounting$owned_entries)
    for (size in c(3L, 7L)) {
      phi <- numeric(f$N)
      d93i_fn("replay_cp04_if")(b, size, 1L, function(z) phi[z$rows] <<- as.vector(z$phi))
      testthat::expect_equal(phi, r$phi, tolerance = 1e-8)
      core <- numeric(f$N)
      d93i_fn("replay_cp04_if")(b, size, 1L, function(z) core[z$rows] <<- as.vector(z$phi), scale = "core")
      testthat::expect_equal(core / sum(f$A), r$phi / f$N, tolerance = 1e-8)
    }
    cv <- d93i_fn("cp04_stream_covariance")(b, 1L)
    # The package Gram is in core units with divisor n_A, not N^2.
    want <- sum(f$A) * sum((r$phi - mean(r$phi))^2) / f$N^2
    testthat::expect_equal(as.vector(cv$covariance), want, tolerance = 1e-8)
  }
})

testthat::test_that("conditional loading integrands use the same independent derivatives", {
  f <- d93i_fixture(); r <- d93i_reference(f, "independent_period_logit", "within_cohort")
  nu <- cbind(r$wp, r$wm); xi <- cbind(r$xu, r$xb); du <- ncol(r$xu)
  # Independent derivative of each unnormalized term, with explicit masks.
  LT <- cbind(r$L, matrix(0, f$N, 2L))
  LM <- cbind(r$L, r$act * f$X)
  wantP <- (r$wp * f$dy / r$means[3L]) * LT - (r$wm * f$dy / r$means[4L]) * LM
  wantB <- cbind(colSums(r$wp * LT), colSums(r$wm * LM))
  got <- d93i_fn("cpt_two_arm_piece")(seq_len(f$N), seq_len(f$N), as.double(f$k == 1L),
    r$h, f$X, r$odds, r$act, f$dy, r$L, xi, r$means[3L], r$means[4L], du, 2L)
  testthat::expect_equal(unname(got$P), unname(wantP), tolerance = 1e-10)
  testthat::expect_equal(unname(got$blam), unname(wantB), tolerance = 1e-10)
  for (h in list(rep(1, f$N), as.double(f$x == 0))) {
    # Conditional numerator, unconditional normalizer; all nuisance roots stay
    # those of the full frame. This is a case-weight derivative, not h*phi.
    R <- colSums(h * nu * f$dy) / (f$N * r$means[3:4])
    ld <- colMeans(h * wantP) - as.vector(wantB %*% (c(1, -1) * R / r$means[3:4])) / f$N
    phi <- h * (r$wp / r$means[3L] - r$wm / r$means[4L]) * f$dy -
      r$wp * R[1L] / r$means[3L] + r$wm * R[2L] / r$means[4L] + as.vector(xi %*% ld)
    gotR <- colSums(h * got$nu * got$rho) / (f$N * r$means[3:4])
    gotL <- colMeans(h * got$P) - as.vector(got$blam %*% (c(1, -1) * gotR / r$means[3:4])) / f$N
    from_piece <- h * as.vector((got$nu * got$rho) %*% (c(1, -1) / r$means[3:4])) -
      as.vector(got$nu %*% (c(1, -1) * gotR / r$means[3:4])) + as.vector(got$xi %*% gotL)
    testthat::expect_equal(from_piece, phi, tolerance = 1e-10)
    if (all(h == 1)) testthat::expect_equal(phi, r$phi, tolerance = 1e-8)
  }
})
