make_cohort_time_panel <- function(incomplete = FALSE) {
  withr::local_seed(14603L)
  d <- expand.grid(t = 1:6, id = 1:72)
  d$g <- c(0, 3, 5)[(d$id - 1) %% 3 + 1]
  d$cluster <- (d$id - 1) %/% 4
  d$y <- rnorm(72)[d$id] + sin(d$t) + rnorm(nrow(d)) +
    ifelse(d$g > 0 & d$t >= d$g, (d$t - d$g + 1) / 2, 0)
  if (incomplete) d$y[d$t > 1 & (d$id + d$t) %% 5 < 2] <- NA_real_
  d
}

fit_cohort_time_fe <- function(d, pool = "nevertreated", ..., point_only = FALSE, method = "reg") {
  cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g", data = d,
    est_method = method, base_period = "universal", control_group = pool,
    pt_assumption = if (pool == "notyettreated") "all" else "post",
    pt_start = if (pool == "notyettreated") min(d$t) else NULL,
    bstrap = FALSE, cband = FALSE, compute_inffunc = !point_only,
    control = list(preset = "cohort_time_fe"), ...
  )
}

dense_cohort_time_fe <- function(d, fit, cluster) {
  d <- d[!is.na(d$y), ]
  ids <- rownames(fit$inffunc)
  id <- match(as.character(d$id), ids)
  event <- outer(paste(d$g, d$t), paste(fit$group, fit$t), "==")
  date <- outer(d$t, fit$DIDparams$tlist[-1L], "==")
  unit <- outer(as.character(d$id), ids, "==")
  Z <- cbind(unit, date, event) * 1.0
  qr_fit <- stats::lm.fit(Z, d$y)
  stopifnot(qr_fit$rank == ncol(Z))
  influence <- t(solve(crossprod(Z), t(rowsum(Z * qr_fit$residuals, id, reorder = TRUE))))
  retained <- ncol(unit) + ncol(date) + seq_along(fit$att)
  target <- influence[, retained, drop = FALSE]
  csum <- rowsum(target, d[[cluster]][match(ids, as.character(d$id))])
  list(att = unname(qr_fit$coefficients[retained]), IF = unname(nrow(target) * target),
       vcov = crossprod(csum))
}

test_that("cohort-by-date FE matches dense OLS and original-unit clustered scores", {
  for (incomplete in c(FALSE, TRUE)) {
    for (pool in c("nevertreated", "notyettreated")) {
      d <- make_cohort_time_panel(incomplete)
      fit <- fit_cohort_time_fe(d, pool, clustervars = "cluster")
      ref <- dense_cohort_time_fe(d, fit, "cluster")
      expect_s3_class(fit, "cdid_MP")
      expect_equal(fit$att, ref$att, tolerance = 1e-12)
      expect_equal(unname(fit$inffunc), ref$IF, tolerance = 1e-12)
      expect_equal(unname(fit$V_analytical) / fit$n, ref$vcov, tolerance = 1e-12)
      expect_equal(fit$se, sqrt(diag(ref$vcov)), tolerance = 1e-12)
      expect_equal(fit_cohort_time_fe(d, pool, point_only = TRUE)$att, fit$att, tolerance = 1e-12)
      if (pool == "notyettreated") {
        expect_true(all(fit$t >= fit$group))
        expect_identical(fit$chained$pretest$status, "not_applicable")
      } else {
        expect_false(any(fit$t == fit$group - 1L))
        expect_identical(fit$chained$pretest$status, "computed")
      }
    }
  }
})

test_that("cohort-by-date FE preserves frame order and singleton score rows", {
  d <- make_cohort_time_panel(TRUE)
  singleton <- d$id == 1L & d$t > 1L
  d$y[singleton] <- NA_real_
  frame <- unique(d[c("id", "g")])
  frame <- rbind(frame[nrow(frame):1L, ], data.frame(id = 73L, g = 0))
  fit <- fit_cohort_time_fe(d, frame = frame)
  native <- fit_cohort_time_fe(d[!(d$id == 1L), ])
  expect_equal(fit$att, native$att, tolerance = 1e-12)
  expect_equal(unname(fit$V_analytical) / fit$n, unname(native$V_analytical) / native$n, tolerance = 1e-12)
  expect_identical(rownames(fit$inffunc), as.character(frame$id))
  expect_true(all(fit$inffunc[c("1", "73"), ] == 0))
})

test_that("cohort-by-date FE refuses requests with unsupported modeling options", {
  d <- make_cohort_time_panel()
  expect_error(fit_cohort_time_fe(d, xformla = ~t), class = "cdid_etwfe_scope_error")
  expect_error(fit_cohort_time_fe(d, xformla = ~offset(cluster)), class = "cdid_etwfe_scope_error")
  expect_error(fit_cohort_time_fe(d, cross_pair = TRUE), class = "cdid_etwfe_scope_error")
  expect_error(fit_cohort_time_fe(d, weightsname = "cluster"), class = "cdid_etwfe_scope_error")
})

test_that("cohort-by-date FE estimates supported effects when another cohort is disconnected", {
  d <- make_cohort_time_panel()
  d$y[d$g == 3L & d$t == 2L] <- NA_real_
  fit <- fit_cohort_time_fe(d)
  unattached <- fit$group == 3L
  expect_true(all(is.na(fit$att[unattached])))
  expect_true(all(fit$chained$support_status[unattached] == "unidentified"))
  expect_true(all(is.na(fit$inffunc[, unattached])))
  expect_true(all(is.finite(fit$att[!unattached])))
  expect_identical(fit$chained$pretest$status, "unavailable")

  empty <- fit_cohort_time_fe(transform(d, y = NA_real_))
  expect_true(all(is.na(empty$att)))
  expect_true(all(is.na(empty$inffunc)))
  expect_true(all(empty$chained$support_status == "no_support"))
})

test_that("cohort-by-date FE reuses package inference and full-frame aggregation", {
  withr::local_seed(9512)
  d <- make_cohort_time_panel()
  fit <- cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g", data = d,
    est_method = "reg", base_period = "universal", biters = 49L,
    control = list(preset = "cohort_time_fe")
  )
  analytic <- fit_cohort_time_fe(d)
  expect_equal(fit$att, analytic$att, tolerance = 1e-12)
  expect_equal(fit$V_analytical, analytic$V_analytical, tolerance = 1e-12)
  expect_true(fit$chained$inference$bstrap)
  expect_true(fit$chained$inference$cband)
  simple <- cdid::aggte(analytic, type = "simple")
  expect_equal(simple$overall.att, mean(analytic$att[analytic$t >= analytic$group]),
               tolerance = 1e-12)
})

test_that("cohort-by-date FE does not populate the GMM panel cache", {
  d <- make_cohort_time_panel()
  cache <- getFromNamespace("cdid_fast_kernel_cache", "cdid")
  prior <- cache$last_panel_prep
  withr::defer(cache$last_panel_prep <- prior)
  cache$last_panel_prep <- NULL
  fit <- fit_cohort_time_fe(d)
  expect_true(all(is.finite(fit$att)))
  expect_null(cache$last_panel_prep)
})

# Fit the conventional dummy and interaction regression. Compute the
# coefficient and cohort-average derivatives directly from its design.
cohort_time_x_oracle <- function(d, frame, cells, cluster) {
  observed <- d[!is.na(d$y), , drop = FALSE]
  ids <- frame$id
  used <- unique(observed$id)
  roster_index <- match(observed$id, ids)
  unit_dummies <- outer(observed$id, used, "==")
  dates <- outer(observed$t, sort(unique(d$t))[-1L], "==")
  events <- outer(paste(observed$g, observed$t), paste(cells$g, cells$time), "==")
  X <- cbind(unit_dummies, dates, dates * observed$x, events, events * observed$x) * 1.0
  dense <- stats::lm.fit(X, observed$y)
  stopifnot(dense$rank == ncol(X))
  n <- length(ids)
  j <- nrow(cells)
  intercept <- ncol(unit_dummies) + 2L * ncol(dates) + seq_len(j)
  slope <- intercept + j
  mean_x <- vapply(cells$g, function(g) mean(frame$x[frame$g == g]), 0.0)
  C <- matrix(0, ncol(X), j)
  C[cbind(intercept, seq_len(j))] <- 1
  C[cbind(slope, seq_len(j))] <- mean_x
  A <- solve(crossprod(X), C)
  score <- rowsum(X %*% A * dense$residuals, roster_index)
  IF <- matrix(0, n, j)
  IF[as.integer(rownames(score)), ] <- n * score
  for (k in seq_len(j)) {
    g <- frame$g == cells$g[k]
    IF[g, k] <- IF[g, k] + n / sum(g) * (frame$x[g] - mean_x[k]) * dense$coefficients[slope[k]]
  }
  csum <- rowsum(IF / n, cluster)
  list(att = unname(drop(crossprod(C, dense$coefficients))), IF = unname(IF),
       vcov = unname(crossprod(csum)))
}

test_that("covariate cohort-by-date FE matches dense OLS and full-cohort standardization", {
  withr::local_seed(4082L)
  for (missing in c(FALSE, TRUE)) {
    for (pool in c("nevertreated", "notyettreated")) {
      d <- make_cohort_time_panel(missing)
      x <- rnorm(72L)
      d$x <- x[d$id]
      d$y <- d$y + (d$t / 3 + 0.2 * (d$g > 0 & d$t >= d$g)) * d$x
      frame <- unique(d[c("id", "g", "x", "cluster")])
      frame <- rbind(frame[nrow(frame):1L, ],
                     data.frame(id = 73L, g = 3, x = 1.4, cluster = 18))
      fit <- fit_cohort_time_fe(d, pool, xformla = ~ x, frame = frame,
                                clustervars = "cluster")
      cells <- data.frame(g = fit$group, time = fit$t)
      target <- cohort_time_x_oracle(d, frame, cells, frame$cluster)
      expect_true(all(fit$chained$available))
      expect_equal(fit$att, target$att, tolerance = 1e-12)
      expect_equal(unname(fit$inffunc), target$IF, tolerance = 1e-12)
      expect_equal(unname(fit$V_analytical) / fit$n, target$vcov, tolerance = 1e-12)
      expect_true(any(fit$inffunc["73", fit$group == 3] != 0))
      expect_true(all(fit$inffunc["73", fit$group != 3] == 0))
      point <- fit_cohort_time_fe(d, pool, xformla = ~ x, frame = frame, point_only = TRUE)
      expect_equal(point$att, fit$att, tolerance = 1e-12)
      expect_null(point$inffunc)
    }
  }
})

test_that("covariate coordinates preserve the effects and isolated unsupported cells", {
  d <- make_cohort_time_panel(incomplete = TRUE)
  d$x <- (d$id %% 7L) / 16
  adjusted <- fit_cohort_time_fe(d, xformla = ~ x)
  shifted <- fit_cohort_time_fe(transform(d, x = x + 2^24), xformla = ~ x)
  expect_equal(shifted$att, adjusted$att, tolerance = 1e-12)
  expect_equal(shifted$inffunc, adjusted$inffunc, tolerance = 1e-12)

  d$x[d$g == 3L] <- 0.25
  d$y[d$g == 3L & d$t == 2L] <- NA_real_
  disconnected <- fit_cohort_time_fe(d, xformla = ~ x)
  removed <- fit_cohort_time_fe(d[d$g != 3L, ], xformla = ~ x)
  keep <- disconnected$group == 5L
  expect_true(all(is.na(disconnected$att[!keep])))
  expect_equal(disconnected$att[keep], removed$att, tolerance = 1e-12)
  expect_equal(unname(disconnected$V_analytical[keep, keep]) / disconnected$n,
               unname(removed$V_analytical) / removed$n, tolerance = 1e-12)
})

test_that("cohort-time treatment scores reduce to the installed complete two-date scores", {
  withr::local_seed(886120L)
  n <- 80L
  u <- data.frame(id = seq_len(n), x = rnorm(n))
  u$g <- 2L * rbinom(n, 1L, plogis(-0.1 + 0.6 * u$x))
  d <- merge(u, data.frame(t = 1:2), sort = FALSE)
  u0 <- rnorm(n)
  dy <- 0.27 + 0.6 * u$x + 0.4 * u$x^2 + 1.3 * (u$g == 2L) + rnorm(n)
  d$y <- u0[d$id] + (d$t == 2) * dy[d$id]
  for (method in c("ipw", "dr")) {
    reference <- cdid::att_gt("y", "t", "id", "g", xformla = ~x, data = d,
      control_group = "nevertreated", base_period = "universal", est_method = method,
      gmm_weighting = "identity", bstrap = FALSE, cband = FALSE)
    j <- which(reference$t == 2L)
    for (weighting in c("identity", "optimal")) {
      fit <- fit_cohort_time_fe(d, xformla = ~x, method = method, gmm_weighting = weighting)
      expect_equal(fit$att, reference$att[j], tolerance = 1e-12)
      expect_equal(unname(fit$inffunc), unname(reference$inffunc[, j, drop = FALSE]),
                   tolerance = 1e-10)
      expect_equal(fit$se, reference$se[j], tolerance = 1e-12)
      point <- fit_cohort_time_fe(d, xformla = ~x, method = method,
                                  gmm_weighting = weighting, point_only = TRUE)
      expect_equal(point$att, fit$att, tolerance = 1e-12)
      expect_null(point$inffunc)
    }
  }
})

cohort_time_transport_panel <- function(correct_propensity, correct_outcome) {
  g <- c(0L, 3L, 5L)
  x <- -1:1
  counts <- if (correct_propensity) rbind(c(2, 2, 2), c(1, 2, 4), c(4, 2, 1)) else
    rbind(c(10, 7, 3), c(4, 11, 5), c(2, 5, 13))
  rows <- expand.grid(x = x, mask = 1:5, g = g)
  rows <- rows[rep(seq_len(nrow(rows)), counts[cbind(match(rows$g, g), match(rows$x, x))]), ]
  rows$id <- seq_len(nrow(rows))
  d <- merge(rows, data.frame(t = 1:6), sort = FALSE)
  masks <- list(c(1, 2, 4), c(2, 3, 6), c(3, 4, 5, 6), 2, numeric())
  observed <- vapply(seq_len(nrow(d)), function(i) d$t[i] %in% masks[[d$mask[i]]], logical(1))
  tau <- (1.15 + 0.35 * d$x + 0.10 * (d$t - d$g)) * (d$g > 0 & d$t >= d$g)
  untreated <- d$id / 100 + 0.25 * d$t + 0.08 * d$t^2 * d$x
  if (!correct_outcome) untreated <- untreated + 0.16 * d$t * d$x^2
  d$y <- ifelse(observed, untreated + tau, NA_real_)
  list(data = d, frame = rows,
       true_x = drop(counts %*% x) / rowSums(counts), g = g)
}

test_that("cohort-time DR recovers cohort-standardized effects with either nuisance correct", {
  for (correct_propensity in c(TRUE, FALSE)) {
    p <- cohort_time_transport_panel(correct_propensity, correct_outcome = !correct_propensity)
    for (control in c("nevertreated", "notyettreated")) {
      fit <- fit_cohort_time_fe(p$data, control, xformla = ~x, frame = p$frame,
                                method = "dr", gmm_weighting = "identity")
      post <- fit$t >= fit$group
      truth <- 1.15 + 0.35 * p$true_x[match(fit$group[post], p$g)] +
        0.1 * (fit$t[post] - fit$group[post])
      expect_equal(fit$att[post], truth, tolerance = 1e-9)
      expect_identical(rownames(fit$inffunc), as.character(p$frame$id))
      expect_true(any(abs(fit$inffunc[p$frame$mask == 5L, post]) > 1e-8))
      expect_equal(unname(fit$V_analytical) / fit$n, crossprod(unname(fit$inffunc) / fit$n),
                   tolerance = 1e-12)
    }
  }
})

test_that("treatment weighting includes a cluster without repeat observations", {
  frame <- data.frame(id = 1:5, g = c(2L, 2L, 0L, 0L, 2L),
                      x = c(0, 1, 0, 1, 1), cluster = c(1L, 1L, 1L, 1L, 2L))
  d <- merge(frame, data.frame(t = 1:2), sort = FALSE)
  d$y <- ifelse(d$t == 2L, 2 * d$x + (d$g == 2L), 0)
  d$y[d$id == 5L & d$t == 2L] <- NA_real_
  for (method in c("ipw", "dr")) {
    fit <- fit_cohort_time_fe(d, frame = frame, xformla = ~x, method = method,
                              clustervars = "cluster", gmm_weighting = "identity")
    expect_equal(fit$att, 2 / 3, tolerance = 1e-10)
    expect_equal(unname(fit$inffunc["5", 1]), -10 / 9, tolerance = 1e-9)
    expect_equal(unname(fit$se), sqrt(2) * 2 / 9, tolerance = 1e-9)
    expect_equal(fit$chained$inference$n_zero_eligible_clusters, 1)
  }
})

test_that("point-only treatment fits do not project influence scores", {
  d <- make_cohort_time_panel(incomplete = TRUE)
  for (method in c("ipw", "dr")) {
    for (weighting in c("identity", "optimal")) {
      fit <- testthat::with_mocked_bindings(
        fit_cohort_time_fe(d, method = method, gmm_weighting = weighting, point_only = TRUE),
        project_gmm_if_stream = function(...) stop("Unexpected influence projection"),
        .package = "cdid"
      )
      expect_true(all(fit$chained$available))
      expect_null(fit$inffunc)
    }
  }
})

cohort_time_corrected_two_date_panel <- function() {
  u <- expand.grid(mask = 0:3, g = c(0L, 2L), x = -1:1, copy = 1:4)
  u$id <- seq_len(nrow(u))
  u$z <- cos(u$id)
  u$cluster <- (u$id + u$copy) %% 8L
  u$rho <- exp(0.1 * u$z)
  study <- u[rep(seq_len(nrow(u)), each = 2L), ]
  study$t <- rep(1:2, nrow(u))
  levels <- cbind(sin(u$id / 7), sin(u$id / 7) + 0.4 * u$x +
                    0.1 * u$x^2 + (u$g == 2L) * (1.4 + 0.3 * u$x) + cos(u$id / 3))
  observed <- cbind(u$mask %% 2L == 1L, u$mask >= 2L)
  levels[!observed] <- NA_real_
  study$y <- c(t(levels))
  list(units = u, study = study, observed = observed, delta = levels[, 2L] - levels[, 1L])
}

cohort_time_corrected_dense <- function(d, weights, method, correction, probability, universe) {
  u <- d$units
  U <- switch(universe, F = rep(TRUE, nrow(u)), A = rowSums(d$observed) == 2L,
              B = rowSums(d$observed) >= 1L)
  pair <- rowSums(d$observed) == 2L
  control <- stats::glm.control(epsilon = 1e-13, maxit = 100L)
  logit <- function(sub, response, design, fit_weights) {
    stats::glm.fit(design[sub, , drop = FALSE], response[sub], weights = fit_weights[sub],
                   family = stats::quasibinomial(), control = control)$coefficients
  }
  W <- cbind(1, u$x)
  C <- cbind(1, u$z)
  reconstructed <- numeric(nrow(u))
  for (g in c(0L, 2L)) {
    k <- u$g == g
    local <- numeric(sum(k))
    if (correction != "ipw") {
      fits <- stats::lm.wfit(W[k & pair, , drop = FALSE], d$delta[k & pair], weights[k & pair])
      local <- drop(W[k, , drop = FALSE] %*% fits$coefficients)
    }
    if (correction != "imputation") {
      pair_k <- k & pair
      Q <- C[pair_k, , drop = FALSE]
      if (universe == "A") {
        # With only two dates and a degenerate A logit, the fixed design
        # in the test uses 1/2, its actual declared probability.
        p <- rep(0.5, sum(pair_k))
      } else if (probability == "pair_logit") {
        p <- plogis(drop(Q %*% logit(k & U, pair, C, weights)))
      } else {
        p <- apply(d$observed, 2L, function(response)
          plogis(drop(Q %*% logit(k & U, response, C, weights))))
        p <- p[, 1L] * p[, 2L]
      }
      denominator <- if (universe == "F") p else
        p * sum(weights[k & U]) / sum(weights[k]) * u$rho[pair_k]
      local[pair[k]] <- local[pair[k]] + (d$delta[pair_k] - local[pair[k]]) / denominator
    }
    reconstructed[k] <- local
  }
  treated <- u$g == 2L
  target <- weighted.mean(reconstructed[treated], weights[treated])
  if (method == "reg") {
    control_fit <- stats::lm.wfit(W[!treated, , drop = FALSE], reconstructed[!treated], weights[!treated])
    donor_at_target <- drop(W[treated, , drop = FALSE] %*% control_fit$coefficients)
    return(target - weighted.mean(donor_at_target, weights[treated]))
  }
  beta <- logit(rep(TRUE, nrow(u)), treated, W, weights)
  odds <- exp(drop(W[!treated, , drop = FALSE] %*% beta))
  if (method == "ipw") return(target - weighted.mean(reconstructed[!treated], weights[!treated] * odds))
  fit <- stats::lm.wfit(W[!treated, , drop = FALSE], reconstructed[!treated], weights[!treated])
  mean_model <- drop(W %*% fit$coefficients)
  target - weighted.mean(mean_model[treated], weights[treated]) -
    weighted.mean(reconstructed[!treated] - mean_model[!treated], weights[!treated] * odds)
}

test_that("corrected cohort-time fits match independently refitted unit derivatives", {
  d <- cohort_time_corrected_two_date_panel()
  n <- nrow(d$units)
  cases <- data.frame(
    method = c("reg", "ipw", "dr", "dr", "reg", "dr", "ipw"),
    correction = c("imputation", "ipw", "dr", "dr", "ipw", "ipw", "imputation"),
    probability = c("pair_logit", "pair_logit", "pair_logit", "independent_period_logit",
                    "pair_logit", "pair_logit", "pair_logit"),
    universe = c("F", "F", "F", "F", "A", "B", "F")
  )
  v <- d$units$z - mean(d$units$z)
  step <- 1e-4
  for (j in seq_len(nrow(cases))) {
    spec <- as.list(cases[j, ])
    fit <- fit_cohort_time_fe(d$study, xformla = ~ x, method = spec$method, clustervars = "cluster",
      selection = "conditional_mtar", selection_est_method = spec$correction,
      selection_model = spec$probability, selection_universe = spec$universe,
      selection_xformla = if (spec$correction != "imputation")
        if (spec$universe == "A") ~ 0 else ~ z,
      selection_rho = if (spec$universe != "F") "rho",
      selection_justification = if (spec$universe != "F") "A supplied ratio for the synthetic check.",
      gmm_weighting = if (spec$method == "reg") "within" else "optimal")
    ref <- function(w) do.call(cohort_time_corrected_dense, c(list(d, w), spec))
    expect_equal(fit$att, ref(rep(1, n)), tolerance = 1e-10)
    gradient <- (ref(1 + step * v) - ref(1 - step * v)) / (2 * step)
    scores <- fit$inffunc[as.character(d$units$id), 1L]
    expect_equal(unname(mean(v * scores)), gradient, tolerance = 1e-8,
                 info = paste(unlist(spec), collapse = "/"))
    cl <- rowsum(scores, d$units$cluster)
    expect_equal(unname(fit$se), sqrt(sum((cl - mean(cl))^2)) / n, tolerance = 1e-10)
    expect_true(all(is.finite(scores)))
  }
})

test_that("corrected cohort-time projection needs only pairs used by available targets", {
  d <- make_cohort_time_panel()
  d$y[d$g == 3L & d$t == 2L] <- NA_real_
  d$y[d$g != 3L & (d$id + d$t) %% 5L < 1L] <- NA_real_
  # The disconnected g = 3 block contains all-success logits for several
  # pairs. It should not block identified effects in the other cohort.
  fit <- fit_cohort_time_fe(d, selection = "conditional_mtar", selection_xformla = ~ 1)
  expect_true(all(is.finite(fit$att[fit$group == 5L])))
  expect_true(all(is.na(fit$att[fit$group == 3L])))
  expect_identical(fit$chained$support_status[fit$group == 3L], rep("unidentified", 5L))
  obs <- fit$chained$cohort_time_fe$observation
  expect_false(any(obs$support$cohort[obs$used_pairs] == 3L))

  only_once <- transform(d, y = ifelse(t == 6L, 0, NA_real_))
  fit <- fit_cohort_time_fe(only_once, selection = "conditional_mtar", selection_xformla = ~ 1)
  expect_true(all(is.na(fit$att)))
  expect_true(all(fit$chained$support_status[fit$t == 6L] == "unidentified"))
  expect_identical(fit$chained$cohort_time_fe$observation$model_diagnostics$status, "no_required_pairs")
})

test_that("cohort-time DR uses both robustness margins with observation adjustment", {
  # Conditional date probabilities are logistic in x with odds 1/2, 1, 2.
  # Replication within each stratum fixes the date-mask frequencies exactly.
  mask_counts <- rbind(c(4, 2, 2, 1) * 4, rep(9, 4), c(1, 2, 2, 4) * 4)
  for (correct_models in c(FALSE, TRUE)) {
    count_g <- if (correct_models) rbind(c(4, 2, 1), c(2, 10, 2)) else
      rbind(c(2, 2, 2), c(1, 2, 4))
    units <- expand.grid(x = -1:1, g = c(0L, 2L), mask = 0:3)
    n <- mask_counts[cbind(units$x + 2L, units$mask + 1L)] *
      count_g[cbind(match(units$g, c(0L, 2L)), units$x + 2L)]
    units <- units[rep(seq_len(nrow(units)), n), ]
    units$id <- seq_len(nrow(units))
    panel <- merge(units, data.frame(t = 1:2), sort = FALSE)
    change <- 0.7 + panel$x * 0.3 + 0.4 * panel$x^2 * (!correct_models) +
      (1.4 + 0.6 * panel$x) * (panel$g == 2L)
    panel$y <- cos(panel$id) + ifelse(panel$t == 2L, change, 0)
    missing <- bitwAnd(panel$mask, bitwShiftL(1L, panel$t - 1L)) == 0L
    panel$y[missing] <- NA_real_
    result <- fit_cohort_time_fe(panel, xformla = ~ x, method = "dr", gmm_weighting = "identity",
      selection = "conditional_mtar", selection_est_method = "dr",
      selection_model = "independent_period_logit",
      selection_xformla = if (correct_models) ~1 else ~x)
    expect_equal(result$att, 1.4 + 0.6 * mean(units$x[units$g == 2L]), tolerance = 1e-10)
    expect_true(all(is.finite(result$inffunc)))
  }
})

test_that("corrected treatment fits ignore outcome pairs in a disconnected target component", {
  frame <- data.frame(id = 1:5, g = c(0L, 0L, 2L, 3L, 3L))
  d <- merge(frame, data.frame(t = 1:3), sort = FALSE)
  d$y <- ifelse(d$id == 1L | (d$id %in% c(3L, 4L) & d$t > 1L),
                (d$id == 4L & d$t == 3L) * 2, NA_real_)
  for (pool in c("nevertreated", "notyettreated")) {
    for (method in c("ipw", "dr")) {
      fit <- fit_cohort_time_fe(d, pool, frame = frame, method = method,
        point_only = TRUE, gmm_weighting = "identity",
        selection = "conditional_mtar", selection_xformla = ~1)
      expect_equal(fit$att[fit$group == 3L & fit$t == 3L], 2, tolerance = 1e-12)
      expect_true(all(is.na(fit$att[fit$group == 2L])))
      obs <- fit$chained$cohort_time_fe$observation
      expect_false(any(obs$support$cohort[obs$used_pairs] == 2L))
    }
  }
})

test_that("the corrected projection reserves space for its caller's live outputs", {
  n <- 100L
  cohort <- rep(c(0L, 2L), each = n / 2L)
  index <- list(unit = rep(seq_len(n), each = 2L), period = rep(1:2, n), y = rep(c(0, 1), n))
  catalog <- etwfe_corrected_catalog(index, cohort, 2L, 1e6)
  geometry <- etwfe_corrected_regression_geometry(catalog, "nevertreated", 1e6)
  obs <- list(n_parameters = 0L, read = function(cell, rows, gradient) {
    list(value = rep(1, length(rows)), coordinates = integer(), gradient = NULL)
  })
  A <- geometry$n_slots
  current <- catalog$entries + 9 * A^2 + 4 * A + nrow(catalog$cells) + 11
  expect_error(etwfe_corrected_projection(catalog, matrix(1, n), geometry, obs,
    compute_scores = TRUE, chunk_rows = 1L, entries = current, extra = n + A),
    class = "cdid_etwfe_workspace_error")
  fit <- etwfe_corrected_projection(catalog, matrix(1, n), geometry, obs,
    compute_scores = TRUE, chunk_rows = 1L, entries = current + n + A, extra = n + A)
  expect_identical(fit$block_rows, 1)
})
