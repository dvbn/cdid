# Host integration tests: ONLY locked installed cdid namespace functions with
# real panel/frame/reader/link owners and the actual CP04 selection logits.
# Expected values come from the frozen independent packaged reference
# (machine-exported CSVs) or base-R lm.wfit, never from production output.
orc_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}
orc_dir <- testthat::test_path("fixtures", "cp06_pooled_or")
orc_loader <- new.env(parent = globalenv())
sys.source(file.path(orc_dir, "load_fixture.R"), envir = orc_loader)
orc_csv <- function(name) orc_loader$por_fixture_csv(orc_dir, name)
orc_raw <- orc_csv("raw.csv"); orc_raw$cohort <- as.numeric(orc_raw$cohort)
orc_cases <- orc_csv("cases.csv"); orc_calendar <- orc_csv("calendar.csv")
orc_pairs <- orc_csv("pairs.csv"); orc_blocks <- orc_csv("blocks.csv"); orc_links <- orc_csv("links.csv")
orc_cohorts <- orc_csv("cohorts.csv")
orc_cohorts$label <- as.numeric(orc_cohorts$scientific_treatment_label)
orc_transport <- orc_csv("transport.csv")
orc_order <- c(seq.int(108L, 2L, -2L), seq.int(107L, 1L, -2L))
orc_c1 <- c("NT_F", "NT_B", "NT_A", "NYT_F", "NYT_B", "NYT_A", "NYT_F_no_intercept", "NYT_F_zero_slope")
orc_qcols <- c("selection_x_neg", "selection_x_zero", "selection_x_pos")
orc_ctrl <- orc_ns("ipw_owner_control")(max_chunk_rows = 23L, max_retained_entries = 1e7,
  max_matrix_entries = 2e7, max_work_ops = 1e10, max_output_entries = 1e7)
orc_gc <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
orc_mx <- function(x) max(c(0, abs(x)))

orc_load <- function(case, order = orc_order) orc_loader$load_pooled_or_reference(orc_dir, case, order)
orc_recipe <- function(frame, intercept, slopes) orc_ns("or_pooled_outcome_recipe")(
  colnames(frame$Z_a), intercept, slopes)
orc_fit <- function(fx, views = fx$views, control = orc_ctrl, recipe = NULL, ...) {
  if (is.null(recipe)) recipe <- orc_recipe(fx$frame, fx$outcome_recipe$intercept,
    fx$outcome_recipe$shared_columns)
  orc_ns("fit_or_corollary1")(fx$reader, views, control, orc_gc, rank_tol = 1e-10,
    outcome_recipe = recipe, ...)
}

# Real owners for modified designs or selection inputs, constructed exactly as
# load_fixture.R; sel holds U, rho and Q aligned with raw (NA off U unless replaced).
orc_selection <- function(case, ids = orc_raw$id) {
  s <- orc_csv(paste0("selection_", orc_cases$reader_universe[orc_cases$case == case], ".csv"))
  s[match(ids, s$id), , drop = FALSE]
}
orc_build <- function(case, X, raw = orc_raw, order = orc_order, sel = orc_selection(case, raw$id),
                      Q = as.matrix(sel[orc_qcols])) {
  N <- nrow(raw); Y <- as.matrix(raw[paste0("Y", seq_len(nrow(orc_calendar)))])
  long <- data.frame(id = rep(raw$id, each = nrow(orc_calendar)),
    time = rep(orc_calendar$time, N), y = as.vector(t(Y)), g = rep(raw$cohort, each = nrow(orc_calendar)))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  if (length(zero) > 1L) long <- long[!long$id %in% raw$id[zero[-1L]], , drop = FALSE]
  panel <- orc_ns("construct_ipw_panel_owner")(long, "id", "time", "y", "g", orc_ctrl)
  X[raw$A == 0L, ] <- NA_real_
  frame <- orc_ns("construct_ipw_frame_owner")(panel, raw$id[order], X[order, , drop = FALSE],
    weights = raw$w[order], control = orc_ctrl, frame_g = raw$cohort[order],
    allow_frame_only = TRUE, design_on_A_only = TRUE)
  reader <- orc_ns("construct_cp04_frame_reader")(frame, raw$id[order], Q[order, , drop = FALSE],
    raw$cohort[order], raw$w[order], orc_ctrl, U = (sel$U == 1L)[order], rho = sel$rho[order],
    universe = orc_cases$reader_universe[orc_cases$case == case],
    justification = "known transport under the frozen independent fixed atom law")
  links <- orc_links[orc_links$case == case, ]; blocks <- orc_blocks[orc_blocks$case == case, ]
  views <- lapply(seq_len(nrow(links)), function(i) {
    b <- blocks[match(links$block_id[i], blocks$block_id), ]
    orc_ns("construct_ipw_link_view")(frame, links$s[i], links$t[i], g = links$target[i],
      K = as.numeric(strsplit(b$controls, "|", fixed = TRUE)[[1L]]), control = orc_ctrl,
      link_id = links$link_id[i], functional_key = links$functional_key[i])
  })
  list(panel = panel, frame = frame, reader = reader, views = views, raw = raw[order, ],
    links = links, blocks = blocks)
}
orc_design <- function(intercept, x = orc_raw$x, x2 = orc_raw$x^2, slopes = c("x", "x2")) {
  X <- cbind(x = x, x2 = x2)[, slopes, drop = FALSE]
  if (intercept) X <- cbind(`(Intercept)` = 1, X)
  X
}
orc_kidx <- function(label) orc_cohorts$reference_cohort_index[match(label, orc_cohorts$label)]
orc_eidx <- function(pair) orc_pairs$reference_pair_index[orc_pairs$s == pair$s & orc_pairs$t == pair$t]

# Frozen fitted weights q = w a_k/(rho p) of the given aligned rows, from the
# frozen points only; the selection logits never involve the outcome design.
orc_q <- function(fx, pair_id, k, rows) {
  pt <- setNames(fx$expectations$points$value, fx$expectations$points$name)
  e <- orc_pairs$reference_pair_index[orc_pairs$pair_id == pair_id]; kk <- orc_kidx(k)
  eta <- as.vector(fx$Q[rows, , drop = FALSE] %*% pt[sprintf("alpha%d_%d_%d", e, kk, 0:2)])
  fx$raw$w[rows] * pt[[sprintf("a%d", kk)]] / (fx$rho[rows] * stats::plogis(eta))
}

# Independent stacked dense WLS on raw rows with frozen weights: local block
# dummies, d = (t_index - s_index) z, one weight per unique block.
orc_dense <- function(fx, case, X, intercept, slopes) {
  links <- orc_links[orc_links$case == case, ]; blocks <- orc_blocks[orc_blocks$case == case, ]
  raw <- fx$raw; X <- X[match(raw$id, orc_raw$id), , drop = FALSE]
  idx <- function(time) match(time, orc_calendar$time)
  rows_of <- function(k, s, t) which(raw$cohort == k & raw[[paste0("R", s)]] == 1 & raw[[paste0("R", t)]] == 1)
  XX <- NULL; YY <- WW <- numeric()
  for (b in seq_len(nrow(blocks))) {
    pr <- orc_pairs[orc_pairs$pair_id == blocks$pair_id[b], ]; s <- idx(pr$s); t <- idx(pr$t)
    for (k in as.numeric(strsplit(blocks$controls[b], "|", fixed = TRUE)[[1L]])) {
      rows <- rows_of(k, s, t); E <- matrix(0, length(rows), nrow(blocks)); E[, b] <- 1
      XX <- rbind(XX, cbind(if (intercept) E, (t - s) * X[rows, slopes, drop = FALSE]))
      YY <- c(YY, raw[[paste0("Y", t)]][rows] - raw[[paste0("Y", s)]][rows])
      WW <- c(WW, orc_q(fx, blocks$pair_id[b], k, rows))
    }
  }
  coef <- stats::lm.wfit(XX, YY, WW, tol = 1e-12)$coefficients
  alpha <- if (intercept) coef[seq_len(nrow(blocks))] else numeric(nrow(blocks))
  beta <- coef[(if (intercept) nrow(blocks) else 0L) + seq_along(slopes)]
  mu <- vapply(seq_len(nrow(links)), function(l) {
    b <- match(links$block_id[l], blocks$block_id); s <- idx(links$s[l]); t <- idx(links$t[l])
    rows <- rows_of(links$target[l], s, t); T <- orc_q(fx, links$pair_id[l], links$target[l], rows)
    r <- raw[[paste0("Y", t)]][rows] - raw[[paste0("Y", s)]][rows] - alpha[b] -
      as.vector((t - s) * X[rows, slopes, drop = FALSE] %*% beta)
    sum(T * r) / sum(T)
  }, numeric(1))
  list(alpha = unname(alpha), beta = unname(beta), mu = mu)
}

# Map actual backend coordinates/blocks to frozen zero-based reference names.
orc_up_names <- function(fx, fit) {
  up <- fit$upstream; lab <- fx$reader$state$cohorts
  vapply(seq_len(up$D), function(u) if (u <= up$n_a) sprintf("a%d", orc_kidx(lab[up$cohort[u]])) else {
    sp <- fit$support_specs[[up$support[u]]]
    sprintf("alpha%d_%d_%d", orc_eidx(sp$pair), orc_kidx(lab[sp$cohort]), up$column[u] - 1L)
  }, character(1))
}
orc_block_ids <- function(fx, fit) vapply(fit$regression_specs, function(sp) {
  pr <- fit$support_specs[[sp$support_ids[1L]]]$pair; bl <- fx$catalog$blocks
  pid <- orc_pairs$pair_id[orc_pairs$s == pr$s & orc_pairs$t == pr$t]
  ctl <- vapply(strsplit(bl$controls, "|", fixed = TRUE), function(x)
    paste(sort(as.numeric(x)), collapse = ","), character(1))
  bl$block_id[bl$pair_id == pid & ctl == paste(fx$reader$state$cohorts[sp$cohorts], collapse = ",")]
}, character(1))

orc_collect <- function(fit, chunk = 23L, batch = length(fit$links), scale = "frame",
                        details = TRUE) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  N <- st$N; L <- length(fit$links); p <- fit$shared$p; G <- length(fit$regressions)
  na <- function(k) matrix(NA_real_, N, k)
  out <- list(phi = na(L), upstream = na(fit$upstream$D), beta = na(p), alpha = na(G),
    direct = na(L), den = na(L), num = na(L), design = lapply(seq_len(L), function(j) na(p)),
    objective = rep(NA_real_, N), beta_score = na(p), alpha_score = na(G),
    block_scores = lapply(seq_len(G), function(b) na(p)), tau = na(length(st$cohorts)),
    eta = na(length(st$cohorts)), lambda = rep(NA_real_, N), upstream_sets = list())
  out$info <- orc_ns("replay_or_corollary1_if")(fit, chunk, batch, function(b) {
    stopifnot(identical(b$id, st$id[b$rows]))
    out$phi[b$rows, b$links] <<- b$phi
    if (!details) return(invisible(NULL))
    d <- b$details
    out$upstream_sets[[length(out$upstream_sets) + 1L]] <<- d$upstream_ids
    out$upstream[b$rows, d$upstream_ids] <<- d$upstream
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
orc_compare <- function(fx, fit, out) {
  pt <- setNames(fx$expectations$points$value, fx$expectations$points$name)
  ref <- fx$expectations$rows; icpt <- fx$outcome_recipe$intercept; st <- fx$reader$state
  slopes <- fx$outcome_recipe$shared_columns; L <- length(fit$links); up <- fit$upstream
  beta_names <- if (length(slopes)) paste0("beta_", slopes) else character()
  un <- orc_up_names(fx, fit); bid <- orc_block_ids(fx, fit); ia <- seq_len(up$n_a)
  pe <- c(a = orc_mx(st$a[up$cohort[ia]] - pt[un[ia]]),
    logit = orc_mx(unlist(lapply(fit$selections, `[[`, "coefficients")) - pt[un[-ia]]),
    mu = orc_mx(vapply(fit$links, `[[`, numeric(1), "mu") - pt[paste0("mu_l", seq_len(L))]),
    den = orc_mx(vapply(fit$links, `[[`, numeric(1), "d") - pt[paste0("denominator_l", seq_len(L))]),
    num = orc_mx(vapply(fit$links, `[[`, numeric(1), "numerator") - pt[paste0("numerator_l", seq_len(L))]),
    objective = abs(fit$shared$objective - pt[["objective"]]),
    alpha = if (icpt) orc_mx(vapply(fit$regressions, `[[`, numeric(1), "alpha") -
      pt[paste0("alpha_", bid)]) else 0,
    beta = orc_mx(fit$shared$beta - pt[beta_names]))
  re <- c(upstream = orc_mx(out$upstream - ref[, un]),
    beta = orc_mx(out$beta - ref[, beta_names, drop = FALSE]),
    alpha = if (icpt) orc_mx(out$alpha - ref[, paste0("alpha_", bid)]) else 0,
    phi = orc_mx(out$phi - fx$expectations$phi), mu = orc_mx(out$phi - ref[, paste0("mu_l", seq_len(L))]),
    den = orc_mx(out$den - ref[, paste0("denominator_l", seq_len(L))]),
    num = orc_mx(out$num - ref[, paste0("numerator_l", seq_len(L))]),
    objective = orc_mx(out$objective - ref[, "objective"]),
    tau = orc_mx(out$tau - ref[, paste0("tau_g", 1:4)]), eta = orc_mx(out$eta - ref[, paste0("eta_g", 1:4)]),
    lambda = orc_mx(out$lambda - ref[, "lambda"]))
  design_point <- design_row <- 0
  for (l in seq_len(L)) {
    lk <- fit$links[[l]]
    if (icpt) for (b in seq_along(bid)) {
      nm <- paste0("design_l", l, "_alpha_", bid[b])
      design_point <- max(design_point, abs(pt[[nm]] - if (b == lk$regression) lk$d else 0))
      design_row <- max(design_row, orc_mx(ref[, nm] - if (b == lk$regression) out$den[, l] else 0))
    }
    for (k in seq_along(slopes)) {
      nm <- paste0("design_l", l, "_beta_", slopes[k])
      design_point <- max(design_point, abs(pt[[nm]] - lk$c[k] * lk$d))
      design_row <- max(design_row, orc_mx(ref[, nm] - out$design[[l]][, k]))
    }
  }
  list(points = c(pe, design = design_point), rows = c(re, design = design_row), up_names = un,
    block_ids = bid)
}

testthat::test_that("Corollary 1 shared backend reproduces frozen points and complete original-unit rows", {
  for (case in orc_c1) {
    fx <- orc_load(case); st <- fx$reader$state; fit <- orc_fit(fx); out <- orc_collect(fit)
    meta <- fx$catalog$case; U <- meta$reader_universe
    testthat::expect_equal(c(st$N, st$n_A, fx$panel$n_total_units, fx$panel$n_singleton_units,
      fx$frame$n_frame_only_units, sum(fx$panel$units$n_obs == 0L)), c(108, 72, 97, 24, 11, 1))
    cmp <- orc_compare(fx, fit, out)
    testthat::expect_lt(max(cmp$points), 1e-10, label = paste(case, "point error"))
    testthat::expect_lt(max(cmp$rows), 3e-7, label = paste(case, "complete row error"))
    testthat::expect_lt(max(abs(colMeans(cbind(out$phi, out$upstream, out$beta, out$alpha)))), 2e-11)
    # Only the actual control/target supports are fitted; unused frozen coordinates need no fit.
    S <- if (meta$catalog == "NT") 6L else 10L
    testthat::expect_equal(c(length(fit$supports), length(fit$selections), fit$upstream$n_a,
      fit$upstream$D), c(S, S, 4, 4 + 3 * S))
    testthat::expect_true(all(cmp$up_names %in% colnames(fx$expectations$rows)))
    testthat::expect_false(anyDuplicated(cmp$up_names) > 0L)
    # The selection fits are exactly the actual CP04 fits, not a separately weighted logit.
    for (sid in seq_along(fit$selections))
      testthat::expect_identical(fit$selections[[sid]],
        orc_ns("cp04_selection_fit")(fx$reader, fit$selection_specs[[sid]], orc_gc, orc_ctrl))
    # One joint fit, three unique blocks, local intercepts only with the formula intercept.
    testthat::expect_equal(c(fit$accounting$n_regression_fits, fit$accounting$n_regression_blocks,
      fit$accounting$n_local_intercepts, fit$accounting$n_shared_slopes,
      fit$accounting$n_selection_fits, fit$accounting$n_rate_fits, fit$accounting$n_treatment_fits),
      c(1, 3, if (meta$intercept == 1L) 3 else 0, meta$shared_slopes, S, 4, 0))
    testthat::expect_setequal(cmp$block_ids, fx$catalog$blocks$block_id)
    testthat::expect_true(all(startsWith(vapply(fit$links, `[[`, character(1), "functional_key"),
      paste0("cmpcat:OR-Cor1-", U, "-shared-slope-WLS-a/(rho*p)|"))))
    testthat::expect_equal(unname(vapply(fit$links, `[[`, character(1), "source_functional_key")),
      unname(vapply(fx$views, `[[`, character(1), "functional_key")))
    # Durations are modeled calendar transitions 1/2/3, never label gaps 3/24/27.
    L <- vapply(fit$regression_specs, `[[`, numeric(1), "transitions")
    gap <- vapply(fit$regression_specs, function(s) {
      pr <- fit$support_specs[[s$support_ids[1L]]]$pair; pr$t - pr$s }, numeric(1))
    testthat::expect_equal(L, unname(fx$outcome_recipe$pair_transitions[
      fx$catalog$blocks$pair_id[match(cmp$block_ids, fx$catalog$blocks$block_id)]]))
    testthat::expect_true(all(L != gap))
    testthat::expect_equal(unname(orc_ns("cp04_payload")(fit)["entries"]), fit$accounting$owned_entries)
    testthat::expect_lte(fit$accounting$owned_entries, fit$accounting$retained_bound - fit$accounting$borrowed_entries)
    testthat::expect_false(fit$accounting$full_F_by_link_IF_retained)
    testthat::expect_false(fit$accounting$inverse_matrices_retained)
    testthat::expect_false(fit$accounting$dense_dummy_regression)
    testthat::expect_true(all(vapply(fit$supports, function(u) is.null(u$design) && is.null(u$moments),
      logical(1))))
    # The ENTIRE profiled objective is retained: RSS/N at the fit, larger elsewhere.
    beta <- fit$shared$beta; N <- st$N
    obj <- function(b) if (meta$intercept == 1L) orc_ns("or_pooled_wls_objective")(fit$shared$summary, b,
      vapply(fit$regressions, `[[`, numeric(1), "alpha")) / N else {
      F <- fit$shared$summary$factor; sum((F[, length(b) + 1L] - F[, seq_along(b), drop = FALSE] %*% b)^2) / N }
    testthat::expect_lt(abs(obj(beta) - fit$shared$objective), 1e-14)
    if (length(beta)) testthat::expect_gt(obj(beta + 0.01) - obj(beta), 0)
    # Known fixed transport on U (never refit or sample-normalized); Q/rho/X absent off their universes.
    tr <- orc_transport[orc_transport$route == U, ]
    testthat::expect_identical(fx$rho[fx$U],
      tr$rho[match(paste(fx$raw$cohort, fx$raw$x), paste(tr$cohort, tr$x))][fx$U])
    testthat::expect_true(all(is.na(fx$Q[!fx$U, ])) && all(is.na(fx$rho[!fx$U])))
    testthat::expect_true(all(is.na(fx$X[fx$raw$A == 0L, ])))
    # Off-U link and upstream rows are exactly zero; full-F share rows remain.
    off <- !fx$U
    testthat::expect_equal(sum(off), c(F = 0, B = 12, A = 36)[[U]])
    testthat::expect_identical(unname(out$phi[off, , drop = FALSE]), matrix(0, sum(off), length(fit$links)))
    testthat::expect_identical(unname(out$upstream[off, , drop = FALSE]), matrix(0, sum(off), fit$upstream$D))
    if (any(off)) testthat::expect_gt(max(abs(out$tau[off, ])), .1)
    if (U == "A") {
      # U = A: a_k^A = 1 with an identically zero IF (retained as a structural zero).
      testthat::expect_identical(st$a[fit$upstream$cohort[1:4]], rep(1, 4L))
      testthat::expect_identical(unname(out$upstream[, 1:4]), matrix(0, 108L, 4L))
    } else {
      # F/B: off-A rows on U carry a_k and logit nuisance rows into every link row.
      testthat::expect_gt(max(abs(out$phi[st$a_index == 0L & fx$U, ])), .1)
    }
  }
})

testthat::test_that("shared a_k and logit rows, core scaling and original-unit overlap are complete", {
  for (case in c("NYT_F", "NYT_B")) {
    fx <- orc_load(case); st <- fx$reader$state; fit <- orc_fit(fx)
    out <- orc_collect(fit); core <- orc_collect(fit, 7L, 2L, "core")
    testthat::expect_equal(core$phi, st$lambda * out$phi, tolerance = 2e-11)
    testthat::expect_equal(core$upstream, out$upstream, tolerance = 2e-11)
    testthat::expect_equal(core$info$normalization_count, 72)
    ref <- fx$expectations$phi
    testthat::expect_lt(max(abs(crossprod(out$phi) - crossprod(ref))) / 108^2, 1e-10)
    testthat::expect_lt(max(abs(crossprod(core$phi) / 72^2 - crossprod(out$phi) / 108^2)), 1e-13)
    # ONE a_k coordinate per cohort, read by that cohort's supports at every pair.
    a_of <- vapply(fit$supports, function(u) u$channels[1L], integer(1))
    k_of <- vapply(fit$support_specs, `[[`, integer(1), "cohort")
    testthat::expect_true(all(tapply(a_of, k_of, function(x) length(unique(x))) == 1L))
    testthat::expect_gt(max(table(k_of)), 1L)
    testthat::expect_identical(fit$upstream$cohort[a_of], k_of)
    testthat::expect_false(anyDuplicated(unlist(lapply(fit$supports, function(u) u$channels[-1L]))) > 0L)
    # Cohort 20 targets link 2 and controls block (2,5): its units enter both rows.
    testthat::expect_gt(abs(crossprod(out$phi)[1L, 2L] / 108^2), 1e-4)
    # Off-A identities are part of the original-unit covariance.
    keepA <- st$a_index > 0L
    testthat::expect_gt(max(abs(crossprod(out$phi[keepA, ]) - crossprod(out$phi))) / 108^2, 3e-7)
    # Replicated raw records keep per-record rows; clusters aggregate by original ID.
    ix <- rep(seq_len(nrow(orc_raw)), each = 2L); raw2 <- orc_raw[ix, ]; original <- raw2$id
    raw2$id <- seq_len(nrow(raw2))
    fx2 <- orc_build(case, orc_design(TRUE)[ix, , drop = FALSE], raw2, seq_len(nrow(raw2)),
      orc_selection(case, original))
    fit2 <- orc_fit(fx2, recipe = orc_recipe(fx2$frame, TRUE, c("x", "x2")))
    out2 <- orc_collect(fit2, 19L, 2L, details = FALSE)
    expected <- ref[match(original, st$id), , drop = FALSE]
    testthat::expect_lt(max(abs(out2$phi - expected)), 3e-7)
    testthat::expect_equal(vapply(fit2$links, `[[`, numeric(1), "mu"),
      vapply(fit$links, `[[`, numeric(1), "mu"), tolerance = 1e-12)
    agg <- rowsum(out2$phi, original, reorder = FALSE); refagg <- rowsum(expected, original, reorder = FALSE)
    testthat::expect_lt(max(abs(crossprod(agg) - crossprod(refagg))) / nrow(raw2)^2, 1e-10)
    testthat::expect_lt(max(abs(crossprod(rowsum(st$lambda * out2$phi, original, reorder = FALSE)) /
      (2 * 72)^2 - crossprod(agg) / nrow(raw2)^2)), 1e-13)
  }
})

testthat::test_that("duplicate consumers share blocks; order, chunks and link batches agree", {
  fx <- orc_load("NYT_B"); st <- fx$reader$state; fit <- orc_fit(fx); base <- orc_collect(fit)
  testthat::expect_identical(fit$links[[1L]]$regression, fit$links[[4L]]$regression)
  testthat::expect_false(identical(fit$links[[1L]]$source_functional_key, fit$links[[4L]]$source_functional_key))
  testthat::expect_identical(base$phi[, 1L], base$phi[, 4L])
  testthat::expect_equal(length(fit$regressions), 3L)
  # A consumer counted as a second objective block would move the shared slope.
  pt <- setNames(fx$expectations$points$value, fx$expectations$points$name)
  blocks <- lapply(seq_len(nrow(fx$catalog$blocks)), function(b) {
    bl <- fx$catalog$blocks[b, ]; pr <- orc_pairs[orc_pairs$pair_id == bl$pair_id, ]
    s <- match(pr$s, orc_calendar$time); t <- match(pr$t, orc_calendar$time); D <- y <- w <- NULL
    for (k in as.numeric(strsplit(bl$controls, "|", fixed = TRUE)[[1L]])) {
      rows <- which(fx$raw$cohort == k & fx$raw[[paste0("R", s)]] == 1 & fx$raw[[paste0("R", t)]] == 1)
      D <- rbind(D, (t - s) * cbind(x = fx$raw$x[rows], x2 = fx$raw$x[rows]^2))
      y <- c(y, fx$raw[[paste0("Y", t)]][rows] - fx$raw[[paste0("Y", s)]][rows])
      w <- c(w, orc_q(fx, bl$pair_id, k, rows))
    }
    orc_ns("or_pooled_wls_block")(D, y, w)
  })
  once <- orc_ns("or_pooled_wls_fit")(orc_ns("or_pooled_wls_merge")(blocks), 1e-10)
  twice <- orc_ns("or_pooled_wls_fit")(orc_ns("or_pooled_wls_merge")(c(blocks[1L], blocks)), 1e-10)
  testthat::expect_lt(max(abs(once$coefficients - pt[c("beta_x", "beta_x2")])), 1e-10)
  testthat::expect_gt(max(abs(twice$coefficients - pt[c("beta_x", "beta_x2")])), 3e-7)
  # Rows/link batches, explicit chunks, frame order and link-view order.
  for (cfg in list(c(1L, 1L), c(7L, 3L), c(23L, 2L))) {
    small <- orc_collect(fit, cfg[1L], cfg[2L], details = FALSE)
    testthat::expect_equal(small$phi, base$phi, tolerance = 2e-11)
  }
  rr <- c(23L, 1L, 108L, 42L); ll <- c(4L, 2L)
  got <- orc_ns("evaluate_or_corollary1_if_chunk")(fit, rr, ll)
  testthat::expect_equal(unname(got$phi), unname(base$phi[rr, ll]), tolerance = 2e-11)
  plain <- orc_load("NYT_B", seq_len(108L)); pfit <- orc_fit(plain)
  pout <- orc_collect(pfit, 13L, 3L, details = FALSE)
  testthat::expect_equal(unname(pout$phi[match(st$id, plain$reader$state$id), ]),
    unname(base$phi), tolerance = 1e-10)
  rev_fit <- orc_fit(fx, rev(fx$views)); rev_out <- orc_collect(rev_fit, 11L, 2L, details = FALSE)
  testthat::expect_equal(rev_fit$shared$beta, fit$shared$beta, tolerance = 1e-12)
  testthat::expect_equal(vapply(rev_fit$regressions, `[[`, numeric(1), "alpha"),
    vapply(fit$regressions, `[[`, numeric(1), "alpha"), tolerance = 1e-12)
  testthat::expect_equal(unname(rev_out$phi[, 4:1]), unname(base$phi), tolerance = 1e-10)
  # Every link batch carries ALL global beta control dependencies (Trap 2): every
  # control block, control support logit and a_k, not only its own pair/pool.
  ctrl <- fit$shared$control_support_ids
  ctrl_up <- sort(unique(unlist(lapply(fit$supports[ctrl], `[[`, "channels"))))
  one <- orc_ns("evaluate_or_corollary1_if_chunk")(fit, 1:23, 2L)
  testthat::expect_equal(one$diagnostics$n_global_beta_blocks, 3L)
  testthat::expect_gte(one$diagnostics$n_selection_score_evaluations, length(ctrl))
  testthat::expect_gte(one$diagnostics$n_upstream_coordinates, length(ctrl_up))
  testthat::expect_equal(unname(one$phi), unname(base$phi[1:23, 2L, drop = FALSE]), tolerance = 2e-11)
  dl <- orc_ns("evaluate_or_corollary1_if_chunk")(fit, 1:23, 2L, details = TRUE)
  testthat::expect_true(all(ctrl %in% dl$details$support_ids))
  testthat::expect_true(all(ctrl_up %in% dl$details$upstream_ids))
  testthat::expect_equal(dl$details$beta_regression_ids, 1:3)
  zfx <- orc_load("NYT_F_zero_slope"); zfit <- orc_fit(zfx)
  zone <- orc_ns("evaluate_or_corollary1_if_chunk")(zfit, 1:23, 2L)
  testthat::expect_equal(zone$diagnostics$n_global_beta_blocks, 0L)
  testthat::expect_lt(zone$diagnostics$n_selection_score_evaluations, length(zfit$supports))
})

testthat::test_that("omitted global beta, a_k, logit and reused-unit effects are rejected", {
  for (case in c("NT_F", "NT_B", "NYT_F", "NYT_B", "NYT_A", "NYT_F_no_intercept", "NYT_F_zero_slope")) {
    fx <- orc_load(case); fit <- orc_fit(fx); o <- orc_collect(fit)
    ref <- fx$expectations$phi; N <- fx$reader$state$N; p <- fit$shared$p; up <- fit$upstream
    icpt <- fx$outcome_recipe$intercept; ia <- seq_len(up$n_a)
    testthat::expect_lt(max(abs(o$phi - ref)), 3e-7)
    testthat::expect_false(anyNA(o$upstream))
    solve_beta <- function(rhs) if (!p) rhs else t(if (icpt)
      orc_ns("or_pooled_wls_gram_solve")(fit$shared$fit, t(rhs), N) else
      orc_ns("or_wls_gram_solve")(fit$shared$fit, t(rhs), N))
    global_beta <- function(drop = integer()) {
      Fb <- fit$shared$forcing_beta; Fb[, drop] <- 0
      solve_beta(o$beta_score + o$upstream %*% t(Fb))
    }
    alpha_of <- function(b, xb, upstream = TRUE) {
      r <- fit$regressions[[b]]
      (o$alpha_score[, b] + if (upstream) as.vector(o$upstream[, r$channel_ids, drop = FALSE] %*%
        r$forcing_alpha) else 0) / r$mass - as.vector(xb %*% r$mean_d)
    }
    link_phi <- function(xb, upstream = TRUE, link_logit = TRUE, local = NULL) sapply(seq_along(fit$links), function(j) {
      l <- fit$links[[j]]; xbl <- if (is.null(local)) xb else local(l$regression)
      ch <- fit$supports[[l$support_id]]$channels; rl <- l$rate_loading
      if (!link_logit) rl[-1L] <- 0
      o$direct[, j] - (if (icpt) alpha_of(l$regression, xbl, upstream) else 0) - as.vector(xbl %*% l$c) +
        if (upstream) as.vector(o$upstream[, ch, drop = FALSE] %*% rl) else 0
    })
    testthat::expect_lt(max(abs(link_phi(global_beta()) - ref)), 3e-7)
    testthat::expect_gt(max(abs(link_phi(solve_beta(o$beta_score), upstream = FALSE) - ref)), 3e-7)
    testthat::expect_gt(max(abs(link_phi(global_beta(), link_logit = FALSE) - ref)), 3e-7)
    if (icpt) testthat::expect_gt(max(abs(sapply(seq_along(fit$links), function(j) {
      l <- fit$links[[j]]; xb <- global_beta()
      o$direct[, j] - alpha_of(l$regression, xb, FALSE) - as.vector(xb %*% l$c) +
        as.vector(o$upstream[, fit$supports[[l$support_id]]$channels, drop = FALSE] %*% l$rate_loading)
    }) - ref)), 3e-7)
    if (p) {
      # Logit channels of EVERY control block force the shared slope.
      testthat::expect_gt(max(abs(link_phi(global_beta(seq_len(up$D)[-ia])) - ref)), 3e-7)
      local_beta <- function(b) {
        r <- fit$regressions[[b]]
        solve_beta(o$block_scores[[b]] + o$upstream[, r$channel_ids, drop = FALSE] %*% t(r$forcing_beta))
      }
      testthat::expect_gt(max(abs(link_phi(NULL, local = local_beta) - ref)), 3e-7)
      testthat::expect_lt(max(abs(Reduce(`+`, lapply(seq_along(fit$regressions), local_beta)) - o$beta)), 1e-10)
      indep <- Reduce(`+`, lapply(seq_along(fit$regressions), function(b) crossprod(local_beta(b))))
      testthat::expect_gt(max(abs(crossprod(o$beta) - indep)) / N^2, 3e-7)
      if (fx$catalog$case$catalog == "NYT") {
        # Trap 2: target cohorts reused as controls elsewhere force beta.
        targets <- unique(vapply(fit$links, `[[`, integer(1), "cohort"))
        ctrl <- fit$shared$control_support_ids
        reused <- ctrl[vapply(fit$support_specs[ctrl], function(s) s$cohort %in% targets, logical(1))]
        testthat::expect_gt(length(reused), 0L)
        drop <- unlist(lapply(fit$supports[reused], function(u) u$channels[-1L]))
        testthat::expect_gt(max(abs(link_phi(global_beta(drop)) - ref)), 3e-7)
        # Mixed-control a_k forcing remains (zero only as the structural U = A zero IF).
        if (fx$catalog$case$route != "A")
          testthat::expect_gt(max(abs(link_phi(global_beta(ia)) - ref)), 3e-7)
      } else {
        # NT: every block's sole control is cohort Inf, so a_Inf is one common
        # scalar of the WHOLE objective and its slope and intercept forcing cancel
        # by algebra; its logit (relative, row-varying) forcing does not.
        testthat::expect_lt(max(abs(fit$shared$forcing_beta[, ia])), 1e-12)
        testthat::expect_lt(max(abs(unlist(lapply(fit$regressions, function(r)
          r$forcing_alpha[r$channel_ids <= up$n_a])))), 1e-12)
        testthat::expect_gt(max(abs(fit$shared$forcing_beta[, -ia])), 1e-4)
      }
    }
  }
})

testthat::test_that("local deficiency, combined prediction support, aliases and branches follow the convention", {
  fx <- orc_load("NYT_F")
  rec <- function(f, intercept, slopes) orc_recipe(f$frame, intercept, slopes)
  bid_of <- function(f, fit) orc_block_ids(list(catalog = list(blocks = f$blocks), reader = f$reader), fit)
  # Cohorts 30/Inf have constant x: blocks with only those controls are locally
  # deficient, the (2,5) block identifies the shared slope for every block.
  xl <- orc_raw$x; xl[orc_raw$cohort %in% c(30, Inf)] <- 1
  fl <- orc_build("NYT_F", orc_design(TRUE, xl, xl^2))
  fit <- orc_fit(fl, recipe = rec(fl, TRUE, c("x", "x2")))
  dense <- orc_dense(fx, "NYT_F", orc_design(TRUE, xl, xl^2), TRUE, c("x", "x2"))
  testthat::expect_equal(fit$shared$beta, dense$beta, tolerance = 1e-9)
  testthat::expect_equal(vapply(fit$regressions, `[[`, numeric(1), "alpha"),
    dense$alpha[match(bid_of(fl, fit), fl$blocks$block_id)], tolerance = 1e-9)
  testthat::expect_equal(vapply(fit$links, `[[`, numeric(1), "mu"), dense$mu, tolerance = 1e-9)
  testthat::expect_length(fit$shared$basis$proposed_columns, 0L)
  testthat::expect_lt(max(abs(colMeans(orc_collect(fit, details = FALSE)$phi))), 2e-11)
  # Every control row has x = 1 while treated rows vary: the pooled null direction
  # cannot recover treated predictions, so it is refused with no fallback.
  xc <- orc_raw$x; xc[orc_raw$cohort == Inf] <- 1
  fc <- orc_build("NT_F", orc_design(TRUE, xc, xc^2))
  err <- tryCatch(orc_fit(fc, recipe = rec(fc, TRUE, c("x", "x2"))), error = identity)
  testthat::expect_s3_class(err, "cdid_or_corollary1_support_error")
  testthat::expect_s3_class(err, "cdid_or_wls_rank_error")
  testthat::expect_match(conditionMessage(err), "treated predictions")
  # Prediction-preserving exact alias x2 = 2x + 3: the x-only basis and fit.
  fa <- orc_build("NYT_F", orc_design(TRUE, orc_raw$x, 2 * orc_raw$x + 3))
  afit <- orc_fit(fa, recipe = rec(fa, TRUE, c("x", "x2")))
  f1 <- orc_build("NYT_F", orc_design(TRUE, slopes = "x"))
  fit1 <- orc_fit(f1, recipe = rec(f1, TRUE, "x"))
  testthat::expect_identical(afit$shared$basis$aliased_columns, "x2")
  testthat::expect_equal(afit$accounting$n_shared_slopes, 1L)
  testthat::expect_equal(unname(afit$shared$basis$alias_map[1L, 1L]), 2, tolerance = 1e-12)
  testthat::expect_equal(unname(afit$shared$basis$block_constants[, 1L]), rep(3, 3L), tolerance = 1e-12)
  testthat::expect_equal(afit$accounting$n_support_scans, 2 * length(afit$supports))
  testthat::expect_equal(orc_collect(afit, details = FALSE)$phi, orc_collect(fit1, details = FALSE)$phi,
    tolerance = 1e-10)
  testthat::expect_equal(vapply(fit1$links, `[[`, numeric(1), "mu"),
    orc_dense(fx, "NYT_F", orc_design(TRUE, slopes = "x"), TRUE, "x")$mu, tolerance = 1e-9)
  # Exactly constant slopes on all control and treated rows: intercept-only fit.
  one <- rep(1, nrow(orc_raw))
  f0 <- orc_build("NYT_F", orc_design(TRUE, one, one))
  fit0 <- orc_fit(f0, recipe = rec(f0, TRUE, c("x", "x2")))
  zfx <- orc_load("NYT_F_zero_slope")
  testthat::expect_equal(fit0$shared$p, 0L)
  testthat::expect_lt(max(abs(orc_collect(fit0, details = FALSE)$phi[match(zfx$reader$state$id,
    f0$reader$state$id), ] - zfx$expectations$phi)), 3e-7)
  # Explicit no-intercept: a column zero on every used row is a homogeneous
  # alias; zero only on controls with treated variation is refused.
  nfx <- orc_load("NYT_F_no_intercept")
  fz <- orc_build("NYT_F_no_intercept", orc_design(FALSE, orc_raw$x, 0 * orc_raw$x))
  zfit <- orc_fit(fz, recipe = rec(fz, FALSE, c("x", "x2")))
  testthat::expect_identical(zfit$shared$basis$aliased_columns, "x2")
  testthat::expect_equal(vapply(zfit$links, `[[`, numeric(1), "mu"),
    orc_dense(nfx, "NYT_F_no_intercept", orc_design(FALSE, slopes = "x"), FALSE, "x")$mu, tolerance = 1e-9)
  xz <- as.numeric(orc_raw$cohort == 10)
  fn <- orc_build("NYT_F_no_intercept", orc_design(FALSE, orc_raw$x, xz))
  testthat::expect_error(orc_fit(fn, recipe = rec(fn, FALSE, c("x", "x2"))),
    class = "cdid_or_corollary1_support_error")
  # The sealed frame is never mutated by alias handling.
  testthat::expect_equal(colnames(fa$frame$Z_a), c("(Intercept)", "x", "x2"))
})

testthat::test_that("explicit recipe, inputs, owners, seals and typed support errors", {
  fx <- orc_load("NYT_B"); fit <- orc_fit(fx); base <- orc_collect(fit, details = FALSE)
  rec <- orc_recipe(fx$frame, TRUE, c("x", "x2"))
  run <- function(...) orc_ns("fit_or_corollary1")(fx$reader, fx$views, orc_ctrl, orc_gc, ...)
  testthat::expect_error(run(), "rank_tol")
  for (tol in list(0, 1, Inf, NA_real_, c(1e-10, 1e-10), "1e-10"))
    testthat::expect_error(run(tol, rec), class = "cdid_or_corollary1_input_error")
  testthat::expect_error(run(1e-10), class = "cdid_or_corollary1_input_error")
  testthat::expect_error(run(1e-10, unclass(rec)), class = "cdid_or_corollary1_input_error")
  testthat::expect_error(run(1e-10, orc_ns("or_pooled_outcome_recipe")(c("(Intercept)", "x"), TRUE, "x")),
    class = "cdid_or_corollary1_input_error")
  testthat::expect_error(orc_ns("fit_or_corollary1")(fx$reader, fx$views, orc_ctrl,
    list(epsilon = 1e-13, maxit = 60L), 1e-10, rec), class = "cdid_ipw_control_error")
  testthat::expect_error(run(1e-10, rec, retained_external_entries = -1), "reservation")
  testthat::expect_error(run(1e-10, rec, automatic_chunk_rows = NA), class = "cdid_or_corollary1_input_error")
  X2 <- orc_design(TRUE); X2[, 1L] <- 2
  f2 <- orc_build("NYT_B", X2)
  testthat::expect_error(orc_fit(f2, recipe = orc_recipe(f2$frame, TRUE, c("x", "x2"))),
    class = "cdid_or_corollary1_input_error")
  testthat::expect_error(orc_ns("fit_or_corollary1")(fx$reader, list(list()), orc_ctrl, orc_gc, 1e-10, rec),
    class = "cdid_ipw_owner_integrity_error")
  testthat::expect_error(orc_ns("fit_or_corollary1")(fx$reader, rep(fx$views[1L], 2L), orc_ctrl, orc_gc,
    1e-10, rec), "Distinct nonempty link IDs")
  empty <- orc_ns("construct_ipw_link_view")(fx$frame, 2, 9, g = 10, K = Inf, control = orc_ctrl,
    require_support = FALSE, link_id = "empty-pair", functional_key = "cmpcat:empty:g10:KInf")
  testthat::expect_error(orc_ns("fit_or_corollary1")(fx$reader, list(empty), orc_ctrl, orc_gc, 1e-10, rec),
    class = "cdid_or_corollary1_support_error")
  # Existing CP04 selection failures are preserved, never refit or redrawn.
  fs <- orc_build("NYT_F", orc_design(TRUE), Q = cbind(one = rep(1, nrow(orc_raw)), duplicate = 1))
  testthat::expect_error(orc_fit(fs, recipe = orc_recipe(fs$frame, TRUE, c("x", "x2"))),
    class = "cdid_cp04_rank_error")
  # Off-U selection predictors and transport are never read: finite values there
  # change nothing, while views from another owner are refused.
  sel <- orc_selection("NYT_B"); tr <- orc_transport[orc_transport$route == "B", ]
  cells <- cbind(as.double(orc_raw$x == -1), as.double(orc_raw$x == 0), as.double(orc_raw$x == 1))
  colnames(cells) <- orc_qcols
  clean_sel <- sel; clean_sel$rho <- tr$rho[match(paste(orc_raw$cohort, orc_raw$x), paste(tr$cohort, tr$x))]
  testthat::expect_identical(clean_sel$rho[sel$U == 1L], sel$rho[sel$U == 1L])
  testthat::expect_true(all(is.finite(clean_sel$rho)) && sum(sel$U != 1L) > 0L)
  clean <- orc_build("NYT_B", orc_design(TRUE), sel = clean_sel, Q = cells)
  cfit <- orc_fit(clean, recipe = orc_recipe(clean$frame, TRUE, c("x", "x2")))
  testthat::expect_equal(orc_collect(cfit, details = FALSE)$phi, base$phi, tolerance = 2e-11)
  testthat::expect_error(orc_ns("fit_or_corollary1")(fx$reader, clean$views, orc_ctrl, orc_gc, 1e-10, rec),
    class = "cdid_ipw_owner_mismatch_error")
  # A further consumer of an existing block reuses its fit and selection logits.
  v <- orc_ns("construct_ipw_link_view")(fx$frame, 2, 5, g = 10, K = c(Inf, 30, 20), control = orc_ctrl,
    link_id = "another-consumer", functional_key = "cmpcat:another-consumer:g10")
  vv <- c(fx$views, list(v)); before <- vv
  dup <- orc_ns("fit_or_corollary1")(fx$reader, vv, orc_ctrl, orc_gc, 1e-10, rec)
  testthat::expect_identical(vv, before)
  testthat::expect_equal(c(dup$accounting$n_regression_blocks, dup$accounting$n_selection_fits), c(3, 10))
  testthat::expect_identical(dup$links[[5L]]$regression, dup$links[[1L]]$regression)
  testthat::expect_equal(dup$links[[5L]]$mu, dup$links[[1L]]$mu, tolerance = 0)
  testthat::expect_equal(dup$shared$beta, fit$shared$beta, tolerance = 0)
  # Invalid chunk arguments and tampered backends are typed refusals.
  for (args in list(list(rows = 1:24, links = 1), list(rows = c(1, 1), links = 1),
    list(rows = c(1, Inf), links = 1), list(rows = 1:3, links = 5),
    list(rows = 1:3, links = 1, scale = "observer"), list(rows = 1:3, links = 1, details = NA))) {
    testthat::expect_error(do.call(orc_ns("evaluate_or_corollary1_if_chunk"), c(list(backend = fit), args)),
      class = "cdid_or_corollary1_input_error")
  }
  bad <- fit; bad$links[[1L]]$mu <- 0
  testthat::expect_error(orc_ns("evaluate_or_corollary1_if_chunk")(bad, 1:3, 1L),
    class = "cdid_or_corollary1_integrity_error")
  bad <- fit; bad$shared$forcing_beta[1L, 1L] <- 99
  testthat::expect_error(orc_ns("replay_or_corollary1_if")(bad, 7L, 1L, function(b) NULL),
    class = "cdid_or_corollary1_integrity_error")
})

testthat::test_that("fit and replay charges are preflighted before any scan or callback", {
  fx <- orc_load("NYT_B"); fit <- orc_fit(fx); st <- fx$reader$state
  for (field in c("max_retained_entries", "max_matrix_entries", "max_work_ops", "max_output_entries")) {
    low <- orc_ctrl; low[[field]] <- 1
    testthat::expect_error(orc_fit(fx, control = low), class = "cdid_ipw_owner_bound_error")
    limited <- fit; limited$ref <- NULL; limited$control[[field]] <- 1
    limited <- orc_ns("ipw_owner_seal")(limited, "or_corollary1_backend")
    called <- FALSE
    testthat::expect_error(orc_ns("replay_or_corollary1_if")(limited, 7L, 2L,
      function(b) called <<- TRUE), class = "cdid_ipw_owner_bound_error")
    testthat::expect_false(called)
  }
  for (field in c("retained_external_entries", "external_work_ops"))
    testthat::expect_error(do.call(orc_fit, c(list(fx = fx), setNames(list(1e11), field))),
      class = "cdid_ipw_owner_bound_error")
  reserved <- orc_fit(fx, retained_external_entries = 1234, external_work_ops = 5678)
  testthat::expect_equal(reserved$accounting$borrowed_entries - fit$accounting$borrowed_entries, 1234)
  testthat::expect_equal(reserved$accounting$work_ops_bound - fit$accounting$work_ops_bound, 5678)
  # Each chunk passes; only the cumulative fit plus replay reservation fails.
  one <- orc_ns("or_corollary1_if_budget")(fit, 7L, 1:4)
  limited <- fit; limited$ref <- NULL
  limited$control$max_work_ops <- ceiling(fit$accounting$work_ops_bound + 2 * one$work)
  limited <- orc_ns("ipw_owner_seal")(limited, "or_corollary1_backend")
  testthat::expect_silent(orc_ns("or_corollary1_if_budget")(limited, 7L, 1:4))
  called <- FALSE
  testthat::expect_error(orc_ns("replay_or_corollary1_if")(limited, 7L, 4L,
    function(b) called <<- TRUE), class = "cdid_ipw_owner_bound_error")
  testthat::expect_false(called)
  # Global dependencies, upstream coordinates and diagnostic fields are charged.
  local <- orc_ns("or_corollary1_if_budget")(fit, 7L, 2L)
  detail <- orc_ns("or_corollary1_if_budget")(fit, 7L, 2L, TRUE)
  testthat::expect_true(all(fit$shared$control_support_ids %in% local$sids))
  testthat::expect_true(all(unlist(lapply(fit$supports[fit$shared$control_support_ids], `[[`, "channels"))
    %in% local$uids))
  testthat::expect_gt(detail$output, local$output)
  got <- orc_ns("evaluate_or_corollary1_if_chunk")(fit, 1:7, 2L, details = TRUE)
  testthat::expect_lte(unname(orc_ns("cp04_payload")(got)["entries"]), got$diagnostics$output_entries_bound)
  # Fit work counts each support's actual cohort rows, not N per support.
  rows <- vapply(fit$support_specs, function(s) st$cohort_end[s$cohort] - st$cohort_start[s$cohort] + 1,
    numeric(1))
  testthat::expect_equal(fit$accounting$support_cohort_row_visits, sum(rows), tolerance = 0)
  testthat::expect_lt(sum(rows), length(rows) * st$N)
  # The envelope is affine in the row count, as the automatic row planner assumes.
  env <- lapply(1:3, function(B) orc_ns("or_corollary1_if_envelope")(fit, B, 1:4))
  for (f in c("output", "scratch", "work"))
    testthat::expect_equal(env[[3L]][[f]] - env[[2L]][[f]], env[[2L]][[f]] - env[[1L]][[f]], tolerance = 0)
})

testthat::test_that("automatic Corollary 1 rows are planned before replay with all concurrent reservations", {
  fx <- orc_load("NYT_B"); fit <- orc_fit(fx, automatic_chunk_rows = TRUE)
  testthat::expect_true(fit$accounting$automatic_chunk_rows)
  links <- seq_along(fit$links); retained <- 12345
  env <- orc_ns("or_corollary1_if_envelope")(fit, 8L, links)
  own <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
  limited <- fit; limited$ref <- NULL
  limited$control$max_matrix_entries <- own + retained + env$output + env$scratch
  limited <- orc_ns("ipw_owner_seal")(limited, "or_corollary1_backend")
  plan <- orc_ns("cp04_preflight_links")(limited, links, retained)
  testthat::expect_equal(plan$chunk_rows, 8, tolerance = 0)
  testthat::expect_equal(plan$work_ops_bound, ceiling(fx$reader$state$N / 8) * env$work, tolerance = 0)
  testthat::expect_lte(retained + plan$peak_concurrent_entries, limited$control$max_matrix_entries)
  testthat::expect_error(orc_ns("cp04_preflight_links")(limited, links, retained, 9L),
    class = "cdid_ipw_owner_bound_error")
  # An explicit (non-automatic) backend never reduces its rows.
  fixed <- orc_fit(fx); fixed$ref <- NULL; fixed$control$max_matrix_entries <- limited$control$max_matrix_entries
  fixed <- orc_ns("ipw_owner_seal")(fixed, "or_corollary1_backend")
  called <- FALSE
  testthat::expect_error(orc_ns("cp04_replay_links")(fixed, links, function(b) called <<- TRUE,
    scale = "frame", chunk_rows = NULL), class = "cdid_ipw_owner_bound_error")
  testthat::expect_false(called)
  reference <- orc_collect(fit, chunk = 23L, details = FALSE)$phi
  actual <- matrix(NA_real_, nrow(reference), ncol(reference)); seen <- integer()
  orc_ns("cp04_replay_links")(limited, links, function(b) {
    actual[b$rows, ] <<- b$phi; seen <<- c(seen, length(b$rows))
  }, scale = "frame", chunk_rows = plan$chunk_rows)
  testthat::expect_lte(max(seen), 8L)
  testthat::expect_equal(actual, unname(reference), tolerance = 2e-11)
  # Complete core-scale covariance through the existing streaming consumer.
  V <- orc_ns("cp04_stream_covariance")(limited, links, retained)
  core_ref <- fx$reader$state$lambda * fx$expectations$phi
  Vref <- crossprod(sweep(core_ref, 2L, colMeans(core_ref))) / 72
  testthat::expect_equal(V$estimates, vapply(fit$links, `[[`, numeric(1), "mu"), tolerance = 0)
  testthat::expect_lt(max(abs(V$covariance - Vref)), 1e-10)
  # A work budget may refuse the preplanned replay (its chunk or its cumulative
  # charge); it never runs a callback and retries another partition.
  tiny <- limited; tiny$ref <- NULL
  tiny$control$max_work_ops <- ceiling(fx$reader$state$N / 8) * env$work - 1
  tiny <- orc_ns("ipw_owner_seal")(tiny, "or_corollary1_backend")
  called <- FALSE
  testthat::expect_error(orc_ns("cp04_replay_links")(tiny, links, function(b) called <<- TRUE,
    chunk_rows = 8L), class = "cdid_ipw_owner_bound_error")
  testthat::expect_false(called)
  # The original cap permits a one-row chunk but refuses the requested
  # cumulative one-row replay before any callback.
  cumulative <- limited; cumulative$ref <- NULL
  cumulative$control$max_work_ops <- fit$accounting$work_ops_bound + env$work
  cumulative <- orc_ns("ipw_owner_seal")(cumulative, "or_corollary1_backend")
  testthat::expect_silent(orc_ns("or_corollary1_if_budget")(cumulative, 1L, links))
  testthat::expect_gt(fx$reader$state$N * orc_ns("or_corollary1_if_envelope")(cumulative, 1L, links)$work,
    cumulative$control$max_work_ops)
  called <- FALSE
  testthat::expect_error(orc_ns("cp04_replay_links")(cumulative, links, function(b) called <<- TRUE,
    chunk_rows = 1L), class = "cdid_ipw_owner_bound_error")
  testthat::expect_false(called)
})
