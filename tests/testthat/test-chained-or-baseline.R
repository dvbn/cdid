# Host integration tests: ONLY locked installed cdid namespace functions with
# real panel/frame/reader/link owners. Expected values come from the frozen
# independent packaged reference (machine-exported CSVs) or base-R lm.wfit.
orb_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}
orb_dir <- testthat::test_path("fixtures", "cp06_pooled_or")
orb_loader <- new.env(parent = globalenv())
sys.source(file.path(orb_dir, "load_fixture.R"), envir = orb_loader)
orb_csv <- function(name) orb_loader$por_fixture_csv(orb_dir, name)
orb_raw <- orb_csv("raw.csv"); orb_raw$cohort <- as.numeric(orb_raw$cohort)
orb_cases <- orb_csv("cases.csv"); orb_calendar <- orb_csv("calendar.csv")
orb_pairs <- orb_csv("pairs.csv"); orb_blocks <- orb_csv("blocks.csv"); orb_links <- orb_csv("links.csv")
orb_cohorts <- orb_csv("cohorts.csv")
orb_cohorts$label <- as.numeric(orb_cohorts$scientific_treatment_label)
orb_order <- c(seq.int(108L, 2L, -2L), seq.int(107L, 1L, -2L))
orb_baseline <- c("NT_baseline", "NYT_baseline", "NYT_baseline_no_intercept", "NYT_baseline_zero_slope")
orb_ctrl <- orb_ns("ipw_owner_control")(max_chunk_rows = 23L, max_retained_entries = 1e7,
  max_matrix_entries = 2e7, max_work_ops = 1e10, max_output_entries = 1e7)
orb_mx <- function(x) max(c(0, abs(x)))

orb_load <- function(case, order = orb_order) orb_loader$load_pooled_or_reference(orb_dir, case, order)
orb_recipe <- function(frame, intercept, slopes) orb_ns("or_pooled_outcome_recipe")(
  colnames(frame$Z_a), intercept, slopes)
orb_fit <- function(fx, views = fx$views, control = orb_ctrl, recipe = NULL) {
  if (is.null(recipe)) recipe <- orb_recipe(fx$frame, fx$outcome_recipe$intercept,
    fx$outcome_recipe$shared_columns)
  orb_ns("fit_or_baseline")(fx$reader, views, control, rank_tol = 1e-10, outcome_recipe = recipe)
}

# Real owners for modified fixed designs, constructed exactly as load_fixture.R.
orb_build <- function(case, X, raw = orb_raw, order = orb_order) {
  N <- nrow(raw); Y <- as.matrix(raw[paste0("Y", seq_len(nrow(orb_calendar)))])
  long <- data.frame(id = rep(raw$id, each = nrow(orb_calendar)),
    time = rep(orb_calendar$time, N), y = as.vector(t(Y)), g = rep(raw$cohort, each = nrow(orb_calendar)))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  if (length(zero) > 1L) long <- long[!long$id %in% raw$id[zero[-1L]], , drop = FALSE]
  panel <- orb_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", orb_ctrl)
  X[raw$A == 0L, ] <- NA_real_
  frame <- orb_ns("construct_ipw_frame_owner")(panel, raw$id[order], X[order, , drop = FALSE],
    weights = raw$w[order], control = orb_ctrl, frame_g = raw$cohort[order],
    allow_frame_only = TRUE, design_on_A_only = TRUE)
  reader <- orb_ns("construct_cp04_frame_reader")(frame, raw$id[order],
    cbind(selection_intercept = rep(1, N))[order, , drop = FALSE], raw$cohort[order],
    raw$w[order], orb_ctrl)
  links <- orb_links[orb_links$case == case, ]; blocks <- orb_blocks[orb_blocks$case == case, ]
  views <- lapply(seq_len(nrow(links)), function(i) {
    b <- blocks[match(links$block_id[i], blocks$block_id), ]
    orb_ns("construct_ipw_link_view")(frame, links$s[i], links$t[i], g = links$target[i],
      K = as.numeric(strsplit(b$controls, "|", fixed = TRUE)[[1L]]), control = orb_ctrl,
      link_id = links$link_id[i], functional_key = links$functional_key[i])
  })
  list(panel = panel, frame = frame, reader = reader, views = views, raw = raw[order, ],
    links = links, blocks = blocks)
}
orb_design <- function(intercept, x = orb_raw$x, x2 = orb_raw$x^2, slopes = c("x", "x2")) {
  X <- cbind(x = x, x2 = x2)[, slopes, drop = FALSE]
  if (intercept) X <- cbind(`(Intercept)` = 1, X)
  X
}

# Independent stacked dense WLS on raw fixture rows: rates from raw sums,
# local block dummies, d = (t_index - s_index) z, one weight per unique block.
orb_dense <- function(case, X, intercept, slopes, raw = orb_raw) {
  links <- orb_links[orb_links$case == case, ]; blocks <- orb_blocks[orb_blocks$case == case, ]
  idx <- function(time) match(time, orb_calendar$time)
  rate <- function(s, t, k) {
    S <- raw[[paste0("R", s)]] * raw[[paste0("R", t)]]; ik <- raw$cohort == k
    sum(raw$w[ik] * S[ik]) / sum(raw$w[ik] * raw$A[ik])
  }
  XX <- NULL; YY <- WW <- numeric()
  for (b in seq_len(nrow(blocks))) {
    pr <- orb_pairs[orb_pairs$pair_id == blocks$pair_id[b], ]; s <- idx(pr$s); t <- idx(pr$t)
    for (k in as.numeric(strsplit(blocks$controls[b], "|", fixed = TRUE)[[1L]])) {
      rows <- which(raw$cohort == k & raw[[paste0("R", s)]] == 1 & raw[[paste0("R", t)]] == 1)
      E <- matrix(0, length(rows), nrow(blocks)); E[, b] <- 1
      XX <- rbind(XX, cbind(if (intercept) E, (t - s) * X[rows, slopes, drop = FALSE]))
      YY <- c(YY, raw[[paste0("Y", t)]][rows] - raw[[paste0("Y", s)]][rows])
      WW <- c(WW, raw$w[rows] / rate(s, t, k))
    }
  }
  coef <- stats::lm.wfit(XX, YY, WW, tol = 1e-12)$coefficients
  alpha <- if (intercept) coef[seq_len(nrow(blocks))] else numeric(nrow(blocks))
  beta <- coef[(if (intercept) nrow(blocks) else 0L) + seq_along(slopes)]
  mu <- vapply(seq_len(nrow(links)), function(l) {
    b <- match(links$block_id[l], blocks$block_id); s <- idx(links$s[l]); t <- idx(links$t[l])
    rows <- which(raw$cohort == links$target[l] & raw[[paste0("R", s)]] == 1 & raw[[paste0("R", t)]] == 1)
    T <- raw$w[rows] / rate(s, t, links$target[l])
    r <- raw[[paste0("Y", t)]][rows] - raw[[paste0("Y", s)]][rows] - alpha[b] -
      as.vector((t - s) * X[rows, slopes, drop = FALSE] %*% beta)
    sum(T * r) / sum(T)
  }, numeric(1))
  list(alpha = unname(alpha), beta = unname(beta), mu = mu)
}

# Map actual backend supports/blocks to frozen zero-based reference coordinates.
orb_rate_names <- function(fx, fit) vapply(fit$support_specs, function(s) {
  e <- fx$catalog$pairs$reference_pair_index[fx$catalog$pairs$s == s$pair$s &
    fx$catalog$pairs$t == s$pair$t]
  k <- orb_cohorts$reference_cohort_index[match(fx$reader$state$cohorts[s$cohort], orb_cohorts$label)]
  sprintf("r%d_%d", e, k)
}, character(1))
orb_block_ids <- function(fx, fit) vapply(fit$regression_specs, function(sp) {
  pr <- fit$support_specs[[sp$support_ids[1L]]]$pair; bl <- fx$catalog$blocks
  pid <- fx$catalog$pairs$pair_id[fx$catalog$pairs$s == pr$s & fx$catalog$pairs$t == pr$t]
  ctl <- vapply(strsplit(bl$controls, "|", fixed = TRUE), function(x)
    paste(sort(as.numeric(x)), collapse = ","), character(1))
  bl$block_id[bl$pair_id == pid & ctl == paste(fx$reader$state$cohorts[sp$cohorts], collapse = ",")]
}, character(1))

orb_collect <- function(fit, chunk = 23L, batch = length(fit$links), scale = "frame",
                        details = TRUE) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  N <- st$N; L <- length(fit$links); p <- fit$shared$p; G <- length(fit$regressions)
  na <- function(k) matrix(NA_real_, N, k)
  out <- list(phi = na(L), rates = na(length(fit$supports)), beta = na(p), alpha = na(G),
    direct = na(L), den = na(L), num = na(L), design = lapply(seq_len(L), function(j) na(p)),
    objective = rep(NA_real_, N), beta_score = na(p), alpha_score = na(G),
    block_scores = lapply(seq_len(G), function(b) na(p)), tau = na(length(st$cohorts)),
    eta = na(length(st$cohorts)), lambda = rep(NA_real_, N), support_sets = list())
  out$info <- orb_ns("replay_or_baseline_if")(fit, chunk, batch, function(b) {
    stopifnot(identical(b$id, st$id[b$rows]))
    out$phi[b$rows, b$links] <<- b$phi
    if (!details) return(invisible(NULL))
    d <- b$details
    out$support_sets[[length(out$support_sets) + 1L]] <<- d$support_ids
    out$rates[b$rows, d$support_ids] <<- d$rates
    out$beta[b$rows, ] <<- d$beta; out$beta_score[b$rows, ] <<- d$beta_score
    for (k in seq_along(d$beta_regression_ids))
      out$block_scores[[d$beta_regression_ids[k]]][b$rows, ] <<- d$beta_block_scores[[k]]
    out$alpha[b$rows, d$regression_ids] <<- d$alpha
    out$alpha_score[b$rows, d$regression_ids] <<- d$alpha_score
    out$direct[b$rows, b$links] <<- d$direct
    out$den[b$rows, b$links] <<- d$denominator; out$num[b$rows, b$links] <<- d$numerator
    for (k in seq_along(b$links)) out$design[[b$links[k]]][b$rows, ] <<- d$design[[k]]
    out$objective[b$rows] <<- d$objective
    out$tau[b$rows, ] <<- d$tau; out$eta[b$rows, ] <<- d$eta; out$lambda[b$rows] <<- d$lambda
  }, scale, details)
  out
}

# Every frozen coordinate that the actual control/target catalog fits.
orb_compare <- function(fx, fit, out) {
  pt <- setNames(fx$expectations$points$value, fx$expectations$points$name)
  ref <- fx$expectations$rows; icpt <- fx$outcome_recipe$intercept
  slopes <- fx$outcome_recipe$shared_columns; L <- length(fit$links)
  beta_names <- if (length(slopes)) paste0("beta_", slopes) else character()
  rn <- orb_rate_names(fx, fit); bid <- orb_block_ids(fx, fit)
  pe <- c(rates = orb_mx(vapply(fit$supports, `[[`, numeric(1), "rate") - pt[rn]),
    mu = orb_mx(vapply(fit$links, `[[`, numeric(1), "mu") - pt[paste0("mu_l", seq_len(L))]),
    den = orb_mx(vapply(fit$links, `[[`, numeric(1), "d") - pt[paste0("denominator_l", seq_len(L))]),
    num = orb_mx(vapply(fit$links, `[[`, numeric(1), "numerator") - pt[paste0("numerator_l", seq_len(L))]),
    objective = abs(fit$shared$objective - pt[["objective"]]),
    alpha = if (icpt) orb_mx(vapply(fit$regressions, `[[`, numeric(1), "alpha") -
      pt[paste0("alpha_", bid)]) else 0,
    beta = orb_mx(fit$shared$beta - pt[beta_names]))
  re <- c(rates = orb_mx(out$rates - ref[, rn]), beta = orb_mx(out$beta - ref[, beta_names, drop = FALSE]),
    alpha = if (icpt) orb_mx(out$alpha - ref[, paste0("alpha_", bid)]) else 0,
    phi = orb_mx(out$phi - fx$expectations$phi), mu = orb_mx(out$phi - ref[, paste0("mu_l", seq_len(L))]),
    den = orb_mx(out$den - ref[, paste0("denominator_l", seq_len(L))]),
    num = orb_mx(out$num - ref[, paste0("numerator_l", seq_len(L))]),
    objective = orb_mx(out$objective - ref[, "objective"]),
    tau = orb_mx(out$tau - ref[, paste0("tau_g", 1:4)]), eta = orb_mx(out$eta - ref[, paste0("eta_g", 1:4)]),
    lambda = orb_mx(out$lambda - ref[, "lambda"]))
  design_point <- design_row <- 0
  for (l in seq_len(L)) {
    lk <- fit$links[[l]]
    if (icpt) for (b in seq_along(bid)) {
      nm <- paste0("design_l", l, "_alpha_", bid[b])
      design_point <- max(design_point, abs(pt[[nm]] - if (b == lk$regression) lk$d else 0))
      design_row <- max(design_row, orb_mx(ref[, nm] - if (b == lk$regression) out$den[, l] else 0))
    }
    for (k in seq_along(slopes)) {
      nm <- paste0("design_l", l, "_beta_", slopes[k])
      design_point <- max(design_point, abs(pt[[nm]] - lk$c[k] * lk$d))
      design_row <- max(design_row, orb_mx(ref[, nm] - out$design[[l]][, k]))
    }
  }
  list(points = c(pe, design = design_point), rows = c(re, design = design_row), rate_names = rn,
    block_ids = bid)
}

testthat::test_that("shared baseline backend reproduces frozen points and complete original-unit rows", {
  for (case in orb_baseline) {
    fx <- orb_load(case); st <- fx$reader$state; fit <- orb_fit(fx); out <- orb_collect(fit)
    meta <- fx$catalog$case
    testthat::expect_equal(c(st$N, st$n_A, fx$panel$n_total_units, fx$panel$n_singleton_units,
      fx$frame$n_frame_only_units, sum(fx$panel$units$n_obs == 0L)), c(108, 72, 97, 24, 11, 1))
    cmp <- orb_compare(fx, fit, out)
    testthat::expect_lt(max(cmp$points), 1e-10, label = paste(case, "point error"))
    testthat::expect_lt(max(cmp$rows), 1e-9, label = paste(case, "row error"))
    testthat::expect_lt(max(abs(colMeans(cbind(out$phi, out$rates, out$beta, out$alpha)))), 2e-11)
    # Only the actual control/target supports are fitted; extra frozen columns are unused.
    testthat::expect_equal(length(fit$supports), if (meta$catalog == "NT") 6L else 10L)
    testthat::expect_true(all(cmp$rate_names %in% colnames(fx$expectations$rows)))
    # One joint fit, three unique blocks, local intercepts only when the formula has one.
    testthat::expect_equal(c(fit$accounting$n_regression_fits, fit$accounting$n_regression_blocks,
      fit$accounting$n_local_intercepts, fit$accounting$n_shared_slopes,
      fit$accounting$n_selection_fits, fit$accounting$n_treatment_fits),
      c(1, 3, if (meta$intercept == 1L) 3 else 0, meta$shared_slopes, 0, 0))
    testthat::expect_setequal(cmp$block_ids, fx$catalog$blocks$block_id)
    # Durations are modeled calendar transitions 1/2/3, never label gaps 3/24/27.
    L <- vapply(fit$regression_specs, `[[`, numeric(1), "transitions")
    gap <- vapply(fit$regression_specs, function(s) {
      pr <- fit$support_specs[[s$support_ids[1L]]]$pair; pr$t - pr$s }, numeric(1))
    testthat::expect_equal(L, unname(fx$outcome_recipe$pair_transitions[
      fx$catalog$blocks$pair_id[match(cmp$block_ids, fx$catalog$blocks$block_id)]]))
    testthat::expect_true(all(L != gap))
    testthat::expect_equal(unname(orb_ns("cp04_payload")(fit)["entries"]), fit$accounting$owned_entries)
    testthat::expect_false(fit$accounting$full_F_by_link_IF_retained)
    testthat::expect_false(fit$accounting$inverse_matrices_retained)
    testthat::expect_null(fit$supports[[1L]]$design)
    # The ENTIRE profiled objective is retained: RSS/N at the fit, larger elsewhere.
    beta <- fit$shared$beta; N <- st$N
    obj <- function(b) if (meta$intercept == 1L) orb_ns("or_pooled_wls_objective")(fit$shared$summary, b,
      vapply(fit$regressions, `[[`, numeric(1), "alpha")) / N else {
      F <- fit$shared$summary$factor; sum((F[, length(b) + 1L] - F[, seq_along(b), drop = FALSE] %*% b)^2) / N }
    testthat::expect_lt(abs(obj(beta) - fit$shared$objective), 1e-14)
    if (length(beta)) testthat::expect_gt(obj(beta + 0.01) - obj(beta), 0)
    # Full-F universe: baseline rate and link rows are exactly zero off A.
    off <- st$a_index == 0L
    testthat::expect_identical(unname(out$phi[off, , drop = FALSE]), matrix(0, sum(off), length(fit$links)))
    testthat::expect_identical(unname(out$rates[off, , drop = FALSE]), matrix(0, sum(off), length(fit$supports)))
    testthat::expect_gt(max(abs(out$tau[off, ])), .1)
  }
})

testthat::test_that("rates, shares, core scaling and original-unit overlap are complete", {
  fx <- orb_load("NYT_baseline"); st <- fx$reader$state; fit <- orb_fit(fx)
  out <- orb_collect(fit); core <- orb_collect(fit, 7L, 2L, "core")
  testthat::expect_equal(core$phi, st$lambda * out$phi, tolerance = 2e-11)
  testthat::expect_equal(core$rates, out$rates, tolerance = 2e-11)
  testthat::expect_equal(core$info$normalization_count, 72)
  ref <- fx$expectations$phi
  testthat::expect_lt(max(abs(crossprod(out$phi) / 108^2 - crossprod(ref) / 108^2)), 1e-12)
  testthat::expect_lt(max(abs(crossprod(core$phi) / 72^2 - crossprod(out$phi) / 108^2)), 1e-13)
  # Cohort 20 targets link 2 and controls block (2,5): its units enter both rows.
  testthat::expect_gt(abs(crossprod(out$phi)[1L, 2L] / 108^2), 1e-4)
  # Repeat rows missing one pair contribute only through that pair's rate.
  e1 <- which(vapply(fit$support_specs, function(s) s$pair$s == 2 && s$pair$t == 5, logical(1)))
  miss <- st$a_index > 0L & !(fx$raw$R1 == 1 & fx$raw$R2 == 1)
  testthat::expect_gt(max(abs(out$rates[miss, e1])), .1)
  # Replicated raw records keep per-record rows; clusters aggregate by original ID.
  ix <- rep(seq_len(nrow(orb_raw)), each = 2L); raw2 <- orb_raw[ix, ]; original <- raw2$id
  raw2$id <- seq_len(nrow(raw2))
  fx2 <- orb_build("NYT_baseline", orb_design(TRUE)[ix, , drop = FALSE], raw2, seq_len(nrow(raw2)))
  fit2 <- orb_ns("fit_or_baseline")(fx2$reader, fx2$views, orb_ctrl, 1e-10,
    orb_recipe(fx2$frame, TRUE, c("x", "x2")))
  out2 <- orb_collect(fit2, 19L, 2L, details = FALSE)
  expected <- ref[match(original, fx$reader$state$id), , drop = FALSE]
  testthat::expect_lt(max(abs(out2$phi - expected)), 1e-9)
  testthat::expect_equal(vapply(fit2$links, `[[`, numeric(1), "mu"),
    vapply(fit$links, `[[`, numeric(1), "mu"), tolerance = 1e-12)
  agg <- rowsum(out2$phi, original, reorder = FALSE); refagg <- rowsum(expected, original, reorder = FALSE)
  testthat::expect_lt(max(abs(crossprod(agg) - crossprod(refagg))) / nrow(raw2)^2, 1e-12)
})

testthat::test_that("duplicate consumers share blocks; order, chunks and link batches agree", {
  fx <- orb_load("NYT_baseline"); fit <- orb_fit(fx); base <- orb_collect(fit)
  testthat::expect_identical(fit$links[[1L]]$regression, fit$links[[4L]]$regression)
  testthat::expect_false(identical(fit$links[[1L]]$source_functional_key, fit$links[[4L]]$source_functional_key))
  testthat::expect_identical(base$phi[, 1L], base$phi[, 4L])
  testthat::expect_equal(length(fit$regressions), 3L)
  # A consumer counted as a second objective block would move the shared slope.
  pt <- setNames(fx$expectations$points$value, fx$expectations$points$name)
  blocks <- lapply(seq_len(nrow(fx$catalog$blocks)), function(b) {
    bl <- fx$catalog$blocks[b, ]; pr <- fx$catalog$pairs[fx$catalog$pairs$pair_id == bl$pair_id, ]
    s <- match(pr$s, orb_calendar$time); t <- match(pr$t, orb_calendar$time); D <- y <- w <- NULL
    for (k in as.numeric(strsplit(bl$controls, "|", fixed = TRUE)[[1L]])) {
      rows <- which(orb_raw$cohort == k & orb_raw[[paste0("R", s)]] == 1 & orb_raw[[paste0("R", t)]] == 1)
      kk <- orb_cohorts$reference_cohort_index[match(k, orb_cohorts$label)]
      D <- rbind(D, (t - s) * cbind(x = orb_raw$x[rows], x2 = orb_raw$x[rows]^2))
      y <- c(y, orb_raw[[paste0("Y", t)]][rows] - orb_raw[[paste0("Y", s)]][rows])
      w <- c(w, orb_raw$w[rows] / pt[[sprintf("r%d_%d", pr$reference_pair_index, kk)]])
    }
    orb_ns("or_pooled_wls_block")(D, y, w)
  })
  once <- orb_ns("or_pooled_wls_fit")(orb_ns("or_pooled_wls_merge")(blocks), 1e-10)
  twice <- orb_ns("or_pooled_wls_fit")(orb_ns("or_pooled_wls_merge")(c(blocks[1L], blocks)), 1e-10)
  testthat::expect_lt(max(abs(once$coefficients - pt[c("beta_x", "beta_x2")])), 1e-10)
  testthat::expect_gt(max(abs(twice$coefficients - pt[c("beta_x", "beta_x2")])), 3e-7)
  # Rows/link batches, explicit chunks, frame order and link-view order.
  for (cfg in list(c(1L, 1L), c(7L, 3L), c(23L, 2L))) {
    small <- orb_collect(fit, cfg[1L], cfg[2L], details = FALSE)
    testthat::expect_equal(small$phi, base$phi, tolerance = 2e-11)
  }
  rr <- c(23L, 1L, 108L, 42L); ll <- c(4L, 2L)
  got <- orb_ns("evaluate_or_baseline_if_chunk")(fit, rr, ll)
  testthat::expect_equal(unname(got$phi), unname(base$phi[rr, ll]), tolerance = 2e-11)
  plain <- orb_load("NYT_baseline", seq_len(108L)); pfit <- orb_fit(plain)
  pout <- orb_collect(pfit, 13L, 3L, details = FALSE)
  testthat::expect_equal(unname(pout$phi[match(fx$reader$state$id, plain$reader$state$id), ]),
    unname(base$phi), tolerance = 1e-10)
  rev_fit <- orb_fit(fx, rev(fx$views)); rev_out <- orb_collect(rev_fit, 11L, 2L, details = FALSE)
  testthat::expect_equal(rev_fit$shared$beta, fit$shared$beta, tolerance = 1e-12)
  testthat::expect_equal(vapply(rev_fit$regressions, `[[`, numeric(1), "alpha"),
    vapply(fit$regressions, `[[`, numeric(1), "alpha"), tolerance = 1e-12)
  testthat::expect_equal(unname(rev_out$phi[, 4:1]), unname(base$phi), tolerance = 1e-10)
  # Every link batch carries ALL global beta control dependencies (Trap 2).
  one <- orb_ns("evaluate_or_baseline_if_chunk")(fit, 1:23, 2L, details = TRUE)
  testthat::expect_true(all(fit$shared$control_support_ids %in% one$details$support_ids))
  testthat::expect_equal(one$details$beta_regression_ids, 1:3)
  testthat::expect_equal(unname(one$phi), unname(base$phi[1:23, 2L, drop = FALSE]), tolerance = 2e-11)
  zfx <- orb_load("NYT_baseline_zero_slope"); zfit <- orb_fit(zfx)
  zone <- orb_ns("evaluate_or_baseline_if_chunk")(zfit, 1:23, 2L)
  testthat::expect_equal(zone$diagnostics$n_global_beta_blocks, 0L)
  testthat::expect_lt(zone$diagnostics$n_rate_score_evaluations, length(zfit$supports))
})

testthat::test_that("omitted global beta, rate and reused-unit effects are rejected", {
  for (case in c("NT_baseline", "NYT_baseline", "NYT_baseline_no_intercept", "NYT_baseline_zero_slope")) {
    fx <- orb_load(case); fit <- orb_fit(fx); o <- orb_collect(fit)
    ref <- fx$expectations$phi; N <- fx$reader$state$N; p <- fit$shared$p
    icpt <- fx$outcome_recipe$intercept
    testthat::expect_lt(max(abs(o$phi - ref)), 1e-9)
    solve_beta <- function(rhs) if (!p) rhs else t(if (icpt)
      orb_ns("or_pooled_wls_gram_solve")(fit$shared$fit, t(rhs), N) else
      orb_ns("or_wls_gram_solve")(fit$shared$fit, t(rhs), N))
    alpha_of <- function(b, xb, rates = TRUE) {
      r <- fit$regressions[[b]]; s <- fit$regression_specs[[b]]$support_ids
      (o$alpha_score[, b] + if (rates) as.vector(o$rates[, s, drop = FALSE] %*% r$forcing_alpha) else 0) /
        r$mass - as.vector(xb %*% r$mean_d)
    }
    link_phi <- function(xb, rates = TRUE, local = NULL) sapply(seq_along(fit$links), function(j) {
      l <- fit$links[[j]]; xbl <- if (is.null(local)) xb else local(l$regression)
      o$direct[, j] - (if (icpt) alpha_of(l$regression, xbl, rates) else 0) - as.vector(xbl %*% l$c) +
        if (rates) l$rate_loading * o$rates[, l$support_id] else 0
    })
    no_rate <- link_phi(solve_beta(o$beta_score), rates = FALSE)
    testthat::expect_gt(max(abs(no_rate - ref)), 3e-7)
    if (p) {
      local_beta <- function(b) {
        s <- fit$regression_specs[[b]]$support_ids
        solve_beta(o$block_scores[[b]] + o$rates[, s, drop = FALSE] %*% t(fit$regressions[[b]]$forcing_beta))
      }
      testthat::expect_gt(max(abs(link_phi(NULL, local = local_beta) - ref)), 3e-7)
      indep <- Reduce(`+`, lapply(seq_along(fit$regressions), function(b) crossprod(local_beta(b))))
      testthat::expect_gt(max(abs(crossprod(o$beta) - indep)) / N^2, 3e-7)
      testthat::expect_lt(max(abs(Reduce(`+`, lapply(seq_along(fit$regressions), local_beta)) - o$beta)), 1e-10)
      if (fx$catalog$case$catalog == "NYT") {
        targets <- unique(vapply(fit$links, `[[`, integer(1), "cohort"))
        ctrl <- fit$shared$control_support_ids
        reused <- ctrl[vapply(fit$support_specs[ctrl], function(s) s$cohort %in% targets, logical(1))]
        keep <- setdiff(ctrl, reused); testthat::expect_gt(length(reused), 0L)
        xb <- solve_beta(o$beta_score + o$rates[, keep, drop = FALSE] %*%
          t(fit$shared$forcing_beta[, keep, drop = FALSE]))
        testthat::expect_gt(max(abs(link_phi(xb) - ref)), 3e-7)
      } else {
        # Trap 1: single-cohort intercept forcing vanishes; slope forcing does not.
        testthat::expect_lt(max(abs(unlist(lapply(fit$regressions, `[[`, "forcing_alpha")))), 1e-12)
        testthat::expect_gt(max(abs(fit$shared$forcing_beta)), 1e-4)
      }
    }
  }
})

testthat::test_that("local deficiency, combined prediction support and exact aliases follow the convention", {
  # Cohorts 30/Inf have constant x: blocks with only those controls are locally
  # deficient, the (2,5) block identifies the shared slope for every block.
  xl <- orb_raw$x; xl[orb_raw$cohort %in% c(30, Inf)] <- 1
  fl <- orb_build("NYT_baseline", orb_design(TRUE, xl, xl^2))
  fit <- orb_ns("fit_or_baseline")(fl$reader, fl$views, orb_ctrl, 1e-10, orb_recipe(fl$frame, TRUE, c("x", "x2")))
  dense <- orb_dense("NYT_baseline", orb_design(TRUE, xl, xl^2), TRUE, c("x", "x2"))
  bid <- orb_block_ids(list(catalog = list(pairs = orb_pairs, blocks = fl$blocks), reader = fl$reader), fit)
  testthat::expect_equal(fit$shared$beta, dense$beta, tolerance = 1e-9)
  testthat::expect_equal(vapply(fit$regressions, `[[`, numeric(1), "alpha"),
    dense$alpha[match(bid, fl$blocks$block_id)], tolerance = 1e-9)
  testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), dense$mu, tolerance = 1e-9)
  testthat::expect_length(fit$shared$basis$proposed_columns, 0L)
  # Every control row has x = 1 while treated rows vary: raw global X on A has
  # full rank, yet the pooled null direction cannot recover treated predictions.
  xc <- orb_raw$x; xc[orb_raw$cohort == Inf] <- 1
  fc <- orb_build("NT_baseline", orb_design(TRUE, xc, xc^2))
  testthat::expect_equal(qr(cbind(1, xc, xc^2)[orb_raw$A == 1, ])$rank, 3L)
  err <- tryCatch(orb_ns("fit_or_baseline")(fc$reader, fc$views, orb_ctrl, 1e-10,
    orb_recipe(fc$frame, TRUE, c("x", "x2"))), error = identity)
  testthat::expect_s3_class(err, "cdid_or_baseline_support_error")
  testthat::expect_s3_class(err, "cdid_or_wls_rank_error")
  testthat::expect_match(conditionMessage(err), "treated predictions")
  # Prediction-preserving exact alias x2 = 2x + 3: the x-only basis and fit.
  fa <- orb_build("NYT_baseline", orb_design(TRUE, orb_raw$x, 2 * orb_raw$x + 3))
  afit <- orb_ns("fit_or_baseline")(fa$reader, fa$views, orb_ctrl, 1e-10, orb_recipe(fa$frame, TRUE, c("x", "x2")))
  f1 <- orb_build("NYT_baseline", orb_design(TRUE, slopes = "x"))
  fit1 <- orb_ns("fit_or_baseline")(f1$reader, f1$views, orb_ctrl, 1e-10, orb_recipe(f1$frame, TRUE, "x"))
  testthat::expect_identical(afit$shared$basis$aliased_columns, "x2")
  testthat::expect_equal(afit$accounting$n_shared_slopes, 1L)
  testthat::expect_equal(unname(afit$shared$basis$alias_map[1L, 1L]), 2, tolerance = 1e-12)
  testthat::expect_equal(unname(afit$shared$basis$block_constants[, 1L]), rep(3, 3L), tolerance = 1e-12)
  testthat::expect_equal(orb_collect(afit, details = FALSE)$phi, orb_collect(fit1, details = FALSE)$phi,
    tolerance = 1e-10)
  testthat::expect_equal(vapply(fit1$links, `[[`, numeric(1), "mu"),
    orb_dense("NYT_baseline", orb_design(TRUE, slopes = "x"), TRUE, "x")$mu, tolerance = 1e-9)
  # Exactly constant slopes on all control and treated rows: intercept-only fit.
  one <- rep(1, nrow(orb_raw))
  f0 <- orb_build("NYT_baseline", orb_design(TRUE, one, one))
  fit0 <- orb_ns("fit_or_baseline")(f0$reader, f0$views, orb_ctrl, 1e-10, orb_recipe(f0$frame, TRUE, c("x", "x2")))
  zfx <- orb_load("NYT_baseline_zero_slope")
  testthat::expect_equal(fit0$shared$p, 0L)
  testthat::expect_equal(unname(orb_collect(fit0, details = FALSE)$phi[match(zfx$reader$state$id, f0$reader$state$id), ]),
    unname(zfx$expectations$phi), tolerance = 1e-9)
  # Explicit no-intercept: a column zero on every used row is a homogeneous
  # alias; zero only on controls with treated variation is refused.
  fz <- orb_build("NYT_baseline_no_intercept", orb_design(FALSE, orb_raw$x, 0 * orb_raw$x))
  zfit <- orb_ns("fit_or_baseline")(fz$reader, fz$views, orb_ctrl, 1e-10, orb_recipe(fz$frame, FALSE, c("x", "x2")))
  testthat::expect_identical(zfit$shared$basis$aliased_columns, "x2")
  testthat::expect_equal(vapply(zfit$links, `[[`, numeric(1), "mu"),
    orb_dense("NYT_baseline_no_intercept", orb_design(FALSE, slopes = "x"), FALSE, "x")$mu, tolerance = 1e-9)
  xz <- as.numeric(orb_raw$cohort == 10)
  fn <- orb_build("NYT_baseline_no_intercept", orb_design(FALSE, orb_raw$x, xz))
  testthat::expect_error(orb_ns("fit_or_baseline")(fn$reader, fn$views, orb_ctrl, 1e-10,
    orb_recipe(fn$frame, FALSE, c("x", "x2"))), class = "cdid_or_baseline_support_error")
  # The sealed frame is never mutated by alias handling.
  testthat::expect_equal(colnames(fa$frame$Z_a), c("(Intercept)", "x", "x2"))
})

testthat::test_that("explicit recipe, inputs, seals and typed support errors", {
  fx <- orb_load("NYT_baseline"); fit <- orb_fit(fx)
  rec <- orb_recipe(fx$frame, TRUE, c("x", "x2"))
  run <- function(...) orb_ns("fit_or_baseline")(fx$reader, fx$views, orb_ctrl, ...)
  testthat::expect_error(orb_ns("fit_or_baseline")(fx$reader, fx$views, orb_ctrl), "rank_tol")
  testthat::expect_error(run(0, rec), "rank_tol")
  testthat::expect_error(run(1e-10), class = "cdid_or_baseline_input_error")
  testthat::expect_error(run(1e-10, unclass(rec)), class = "cdid_or_baseline_input_error")
  testthat::expect_error(run(1e-10, orb_ns("or_pooled_outcome_recipe")(c("(Intercept)", "x", "x2"),
    FALSE, c("(Intercept)", "x", "x2"))), NA)
  testthat::expect_error(orb_ns("or_pooled_outcome_recipe")("x", FALSE, character()),
    class = "cdid_or_pooled_input_error")
  testthat::expect_error(orb_ns("or_pooled_outcome_recipe")(c("a", "b", "x"), TRUE, "x"),
    class = "cdid_or_pooled_input_error")
  testthat::expect_error(run(1e-10, orb_ns("or_pooled_outcome_recipe")(c("(Intercept)", "x"), TRUE, "x")),
    class = "cdid_or_baseline_input_error")
  X2 <- orb_design(TRUE); X2[, 1L] <- 2
  f2 <- orb_build("NYT_baseline", X2)
  testthat::expect_error(orb_ns("fit_or_baseline")(f2$reader, f2$views, orb_ctrl, 1e-10,
    orb_recipe(f2$frame, TRUE, c("x", "x2"))), class = "cdid_or_baseline_input_error")
  testthat::expect_error(run(1e-10, rec, retained_external_entries = -1), "reservation")
  testthat::expect_error(orb_ns("fit_or_baseline")(fx$reader, list(list()), orb_ctrl, 1e-10, rec),
    class = "cdid_ipw_owner_integrity_error")
  testthat::expect_error(orb_ns("fit_or_baseline")(fx$reader, rep(fx$views[1L], 2L), orb_ctrl, 1e-10, rec),
    "Distinct nonempty link IDs")
  empty <- orb_ns("construct_ipw_link_view")(fx$frame, 2, 9, g = 10, K = Inf, control = orb_ctrl,
    require_support = FALSE, link_id = "empty-pair", functional_key = "cmpcat:empty:g10:KInf")
  testthat::expect_error(orb_ns("fit_or_baseline")(fx$reader, list(empty), orb_ctrl, 1e-10, rec),
    class = "cdid_or_baseline_support_error")
  for (args in list(list(rows = 1:24, links = 1), list(rows = c(1, 1), links = 1),
    list(rows = c(1, Inf), links = 1), list(rows = 1:3, links = 5),
    list(rows = 1:3, links = 1, scale = "observer"), list(rows = 1:3, links = 1, details = NA))) {
    testthat::expect_error(do.call(orb_ns("evaluate_or_baseline_if_chunk"), c(list(backend = fit), args)),
      class = "cdid_or_baseline_input_error")
  }
  bad <- fit; bad$links[[1L]]$mu <- 0
  testthat::expect_error(orb_ns("evaluate_or_baseline_if_chunk")(bad, 1:3, 1L),
    class = "cdid_or_baseline_integrity_error")
  bad <- fit; bad$shared$beta <- bad$shared$beta + 1
  testthat::expect_error(orb_ns("replay_or_baseline_if")(bad, 7L, 1L, function(b) NULL),
    class = "cdid_or_baseline_integrity_error")
})

testthat::test_that("fit and replay charges are preflighted before any scan or callback", {
  fx <- orb_load("NYT_baseline"); fit <- orb_fit(fx)
  for (field in c("max_retained_entries", "max_matrix_entries", "max_work_ops", "max_output_entries")) {
    low <- orb_ctrl; low[[field]] <- 1
    testthat::expect_error(orb_fit(fx, control = low), class = "cdid_ipw_owner_bound_error")
    limited <- fit; limited$ref <- NULL; limited$control[[field]] <- 1
    limited <- orb_ns("ipw_owner_seal")(limited, "or_baseline_backend")
    called <- FALSE
    testthat::expect_error(orb_ns("replay_or_baseline_if")(limited, 7L, 2L,
      function(b) called <<- TRUE), class = "cdid_ipw_owner_bound_error")
    testthat::expect_false(called)
  }
  # Each chunk passes; only the cumulative replay reservation fails.
  one <- orb_ns("or_baseline_if_budget")(fit, 7L, 1:4)
  limited <- fit; limited$ref <- NULL; limited$control$max_work_ops <- ceiling(one$work * 2)
  limited <- orb_ns("ipw_owner_seal")(limited, "or_baseline_backend")
  called <- FALSE
  testthat::expect_error(orb_ns("replay_or_baseline_if")(limited, 7L, 4L,
    function(b) called <<- TRUE), class = "cdid_ipw_owner_bound_error")
  testthat::expect_false(called)
  # Global dependencies and diagnostic fields are charged, not just local blocks.
  local <- orb_ns("or_baseline_if_budget")(fit, 7L, 2L)
  detail <- orb_ns("or_baseline_if_budget")(fit, 7L, 2L, TRUE)
  testthat::expect_true(all(fit$shared$control_support_ids %in% local$sids))
  testthat::expect_gt(detail$output, local$output)
  got <- orb_ns("evaluate_or_baseline_if_chunk")(fit, 1:7, 2L, details = TRUE)
  testthat::expect_lte(unname(orb_ns("cp04_payload")(got)["entries"]), got$diagnostics$output_entries_bound)
  # The fit charge covers the verification rescan used by exact alias handling.
  fa <- orb_build("NYT_baseline", orb_design(TRUE, orb_raw$x, 2 * orb_raw$x + 3))
  afit <- orb_ns("fit_or_baseline")(fa$reader, fa$views, orb_ctrl, 1e-10, orb_recipe(fa$frame, TRUE, c("x", "x2")))
  testthat::expect_equal(afit$accounting$n_support_scans, 2 * length(afit$supports))
  testthat::expect_lte(afit$accounting$owned_entries, afit$accounting$retained_bound)
})
