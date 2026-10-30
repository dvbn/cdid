# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 general-route (missing-panel) baseline DR IF envelope: the generic
# branch of `dr_baseline_if_envelope`.
#
# Expected values are HAND COUNTED here from the sealed backend's own catalog
# (support specs, pools, treatment specs, links, shared control supports) and,
# for one configuration, from the literal structure of the designed fixture;
# none is produced by calling the installed envelope. Covered: unequal
# per-model support incidence m_t, block-support incidence I_beta/I_g, the
# largest number of requested supports sharing one cohort (hits), one model's
# live transients, details on/off, p_or = 0 and p_or > 0, link subsets, exact
# affine row coefficients and the UNCHANGED automatic planner at adjacent
# thresholds, the actual evaluated payload against the bound, work-cap
# adjacency of the chunk evaluator, D60 reference binding negatives,
# chunk-partition invariance of the full IF and streamed covariance against a
# test-local clone carrying the FROZEN previous envelope, explicit controls and
# the public DR call on the same missing panel.
#
# Fixtures are self-contained: test-local panel/frame/reader/link-view
# constructors over the installed owners. `api` is the actual installed
# package namespace. `ge_old_ns()` is a disclosed TEST-LOCAL clone whose
# `dr_baseline_if_envelope` is a verbatim copy of the previous (dense) generic
# branch; it is never a source fallback and is compared, not trusted.
#
# STATUS: UNEXECUTED at authoring time (no R, build, install or fit was run in
# the authoring job). Root owns sole R and measured validation.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures (test-local, self-contained)
# ---------------------------------------------------------------------------

ge_ctrl <- function(chunk = 7L, work = 1e11, matrix = 2e7, retained = 1e7, output = 1e7) {
  api$ipw_owner_control(
    max_chunk_rows = chunk,
    max_retained_entries = retained,
    max_matrix_entries = matrix,
    max_work_ops = work,
    max_output_entries = output
  )
}
ge_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
ge_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))

# Deterministic panel with cohorts 3, 4 (units 5:9) and never-treated (Inf).
# `drop` removes single (id, t) observations, so the panel is NOT complete and
# every baseline DR fit below takes the general (missing-panel) route.
ge_panel <- function(N = 41L, T = 5L, seed = 19L, cohort4 = TRUE,
                     drop = list(c(5L, 3L), c(12L, 2L), c(30L, 5L))) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:4] <- c(3, Inf, 3, Inf)
  if (cohort4) g[5:9] <- 4
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  u <- d$id
  x1 <- as.numeric((u %% 2L) == 1L)
  x2 <- cos(u * 0.37) + 0.25 * sin(u * 0.73)
  d$g <- g[u]
  d$x <- x1
  d$x2 <- x2
  eff <- ifelse(is.finite(d$g) & d$t >= d$g, 1.1 + 0.2 * (d$t - d$g), 0)
  d$y <- 0.15 * u + 0.3 * d$t + 0.35 * sin(0.61 * u * d$t) +
    0.25 * x1 * (d$t^2) - 0.18 * x2 * d$t + eff + stats::rnorm(nrow(d))
  d$cl <- (u - 1L) %/% 3L + 1L
  keep <- rep(TRUE, nrow(d))
  for (o in drop) keep[d$id == o[1L] & d$t == o[2L]] <- FALSE
  d[keep, , drop = FALSE]
}

# Seven links with UNEQUAL control sets: pairs (1,2), (1,3) and one (2,3) link
# use controls {4, Inf} (three cohorts per model, two control supports per
# block); the others use the never-treated pool only (two cohorts, one control
# support). Pair (2,3) appears twice with different pools, so one target
# support and one control support are shared by two models/blocks.
ge_links <- list(
  list(s = 1, t = 2, g = 3, K = c(4, Inf)),
  list(s = 1, t = 3, g = 3, K = c(4, Inf)),
  list(s = 2, t = 3, g = 3, K = Inf),
  list(s = 1, t = 4, g = 3, K = Inf),
  list(s = 2, t = 4, g = 3, K = Inf),
  list(s = 3, t = 5, g = 3, K = Inf),
  list(s = 2, t = 3, g = 3, K = c(4, Inf))
)
ge_link_id <- function(i, l) sprintf("L%d:g%g:%g_%g:K%s", i, l$g, l$s, l$t, paste(l$K, collapse = "+"))

# Designed structure of `ge_links` on this panel (literal hand count).
ge_structure <- list(
  n_supports = 15L,           # 6 target supports (one per pair) + 9 control supports
  n_blocks = 7L,              # one pool per link (the two (2,3) links differ in controls)
  n_models = 7L,              # one treatment model per link
  m_t_sorted = c(2, 2, 2, 2, 3, 3, 3),
  block_sizes_sorted = c(1, 1, 1, 1, 2, 2, 2),
  block_incidence = 10,       # sum over blocks of |support_ids|
  cohort_support_max = 6L,    # cohort 3 (and Inf) is the cohort of 6 requested supports
  n_control_supports = 9L,
  n_cohorts = 3L
)

ge_design <- function(N = 41L, kind = "x_bin") {
  u <- seq_len(N)
  x1 <- as.numeric((u %% 2L) == 1L)
  if (identical(kind, "icpt")) {
    # p_or = 0, p_trt = 1
    list(X = matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)")), recipe = ge_icpt, trt = "(Intercept)")
  } else if (identical(kind, "icpt_trt_x")) {
    # p_or = 0 with a nontrivial propensity design, p_trt = 2
    list(X = cbind(`(Intercept)` = 1, x = x1), recipe = ge_icpt, trt = c("(Intercept)", "x"))
  } else {
    # p_or = 1, p_trt = 2
    list(X = cbind(`(Intercept)` = 1, x = x1),
      recipe = api$or_pooled_outcome_recipe(c("(Intercept)", "x"), TRUE, "x"), trt = c("(Intercept)", "x"))
  }
}

ge_fixture <- function(d = ge_panel(), chunk = 7L, X = NULL, ctrl = ge_ctrl(chunk), links = ge_links) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  u <- d[!duplicated(d$id), c("id", "g")]
  u <- u[order(u$id), ]
  N <- nrow(u)
  Z <- if (is.null(X)) matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)")) else X[u$id, , drop = FALSE]
  w_vec <- rep(1, N)
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = w_vec, control = ctrl, frame_g = u$g,
    allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, u$id, Z, u$g, w_vec, ctrl, allow_unknown_off_A = TRUE)
  views <- vector("list", length(links))
  for (i in seq_along(links)) {
    l <- links[[i]]
    views[[i]] <- api$construct_ipw_link_view(fo, l$s, l$t, g = l$g, K = l$K, control = ctrl,
      link_id = ge_link_id(i, l),
      functional_key = sprintf("cmpcat:ge:%d:g%g:%g-%g:K%s", i, l$g, l$s, l$t, paste(l$K, collapse = "+")))
  }
  list(d = d, po = po, fo = fo, reader = reader, views = views, ctrl = ctrl, u = u, links = links)
}

ge_fit <- function(fx, des, ctrl = fx$ctrl, automatic = FALSE, ns = api) {
  fit <- ns$fit_dr_baseline(fx$reader, fx$views, ctrl, glm_control = ge_glm, rank_tol = 1e-10,
    outcome_recipe = des$recipe, treatment_columns = des$trt, automatic_chunk_rows = automatic)
  if (!is.null(fit$complete_support)) stop("fixture must take the general route", call. = FALSE)
  fit
}

ge_state <- function(fit) get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
ge_frame <- function(fit) {
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  get("frame_owner", reader$frame_ref, inherits = FALSE)
}

# ---------------------------------------------------------------------------
# Independent hand count of the general envelope from the backend catalog
# ---------------------------------------------------------------------------

# Buffer inventory per chunk row (see the test header and the generic branch):
#  support stage  rows/id/ci/w/A (5) + up (ns) + one support's read (16)
#  OR core        hits (2 + p_or) h_max + resid cache (3 + p_or)(I_beta + I_g + nl)
#                 + 8 p_or + n_ctrl + 2 ng + c_max + 2 nl + 4 (p_or + 4)
#                 [details: + p_or nb + (2 + p_or) nl + 5 p_or + 8]
#  models         retained (5 nt | sum(3 p_trt + m_t + 16)) + one model's
#                 transients m_max + 13 p_trt + 2 p_or + 52 [details: + 32]
#  assembly       3 nl [details: + 5 nl + 2 J + 6]
# These fixtures carry no reference reader, so the D60 terms are zero here
# (the negative cases below show a forged binding is refused, never charged).
ge_hand <- function(fit, B, links, details = FALSE) {
  st <- ge_state(fit); frame <- ge_frame(fit)
  stopifnot(is.null(fit$reference_ref))
  p_or <- fit$shared$p; p_trt <- fit$treatment_recipe$p
  G <- length(fit$regression_specs); J <- length(st$cohorts)
  lk <- fit$links[links]
  gids <- sort(unique(vapply(lk, `[[`, integer(1), "regression")))
  beta_gids <- if (p_or > 0L || details) seq_len(G) else integer()
  tids <- sort(unique(vapply(lk, `[[`, integer(1), "model")))
  sids <- sort(unique(c(vapply(lk, `[[`, integer(1), "support_id"),
    unlist(lapply(fit$regression_specs[union(gids, beta_gids)], `[[`, "support_ids")),
    unlist(lapply(fit$treatment_specs[tids], `[[`, "support_ids")))))
  nl <- length(links); ns <- length(sids); nt <- length(tids); nb <- length(beta_gids); ng <- length(gids)
  m_t <- vapply(fit$treatment_specs[tids], function(s) length(s$support_ids), numeric(1))
  m_max <- max(m_t)
  size <- function(g) vapply(fit$regression_specs[g], function(s) length(s$support_ids), numeric(1))
  I_beta <- sum(size(beta_gids)); I_g <- sum(size(gids)); c_max <- max(c(0, size(gids)))
  n_ctrl <- length(fit$shared$control_support_ids)
  coh <- vapply(fit$support_specs[sids], function(s) as.integer(s$cohort), integer(1))
  h_max <- max(tabulate(coh))
  id_max <- frame$id_byte_metadata$max
  link_bytes <- sum(nchar(vapply(lk, `[[`, character(1), "link_id"), type = "bytes"))
  width <- 4 * ns + 4 * (p_or + 1) * (nb + 1) + 4 * ng + sum(3 * p_trt + m_t + 16) + 8 * nl * (p_or + 4) +
    (if (details) 3 * J + 8 * nl + 16 else 8)
  metadata0 <- 16384 + 1024 * (ns + G + nt + nl + J) + 8 * link_bytes
  support_live <- 5 + ns + 16
  core_live <- (2 + p_or) * h_max + (3 + p_or) * (I_beta + I_g + nl) + 8 * p_or + n_ctrl + 2 * ng + c_max +
    2 * nl + 4 * (p_or + 4) + (if (details) p_or * nb + (2 + p_or) * nl + 5 * p_or + 8 else 0)
  trt_live <- (if (details) sum(3 * p_trt + m_t + 16) else 5 * nt) +
    m_max + 13 * p_trt + 2 * p_or + 52 + (if (details) 32 else 0)
  tail_live <- 3 * nl + (if (details) 5 * nl + 2 * J + 6 else 0)
  live <- support_live + core_live + trt_live + tail_live
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  units <- (ns + 1) * (p_or + 2)^2 * log_a + (p_or + 1) * (p_or * (nb + ng) + I_beta + I_g) +
    nt * ((p_trt + 2)^2 * log_a + p_or + 8) + sum((p_trt + 10) * (m_t + 4)) + J + nl * (p_or + p_trt + 4)
  I_all <- I_beta + I_g + nl
  lookups <- (I_all + ng + nl + 1) * (ns + ng + 1) + 2 * I_all^2 + nt * (ns + m_max + 1)
  small <- 128 * ((p_or + 2)^2 + (p_trt + 2)^2)
  metadata <- metadata0 + B * id_max
  list(
    sids = sids, tids = tids, gids = gids, beta_gids = beta_gids, m_t = m_t, ns = ns, nt = nt, G = G, J = J,
    I_beta = I_beta, I_g = I_g, c_max = c_max, n_ctrl = n_ctrl, h_max = h_max, log_a = log_a,
    width = width, live = live, units = units, lookups = lookups, metadata0 = metadata0, id_max = id_max,
    output = B * (width + 16) + metadata,
    scratch = 2 * B * live + small + metadata,
    work = 128 * B * units + 64 * lookups + 64 * metadata,
    slope = width + 16 + 2 * live + 2 * id_max,
    intercept = 2 * metadata0 + small,
    output_slope = width + 16 + id_max,
    output_intercept = metadata0
  )
}

# Disclosed TEST-LOCAL clone: the evaluator, budget and replay of the installed
# package re-bound to an environment whose `dr_baseline_if_envelope` is a
# verbatim copy of the PREVIOUS generic branch (dense nt * ns width, 64x
# scratch, dense block-by-support work). Everything else resolves to `api`.
ge_old_envelope <- function(backend, B, links, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p_or <- backend$shared$p
  p_trt <- backend$treatment_recipe$p
  G <- length(backend$regression_specs)
  J <- length(st$cohorts)
  sets <- or_pooled_if_sets(backend, links, details)
  gids <- sets$gids
  beta_gids <- sets$beta_gids
  tids <- sort(unique(vapply(backend$links[links], `[[`, integer(1), "model")))
  trt_sids <- sort(unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids"))))
  sids <- sort(unique(c(sets$sids, trt_sids)))
  nl <- length(links)
  ns <- length(sids)
  nt <- length(tids)
  if (!is.null(backend$complete_support)) {
    sets$sids <- sids
    return(dr_baseline_complete_if_envelope(backend, B, links, details, sets, tids, J, frame))
  }
  width <- 4 * ns + 4 * (p_or + 1) * (length(beta_gids) + 1) + 4 * length(gids) +
    nt * (2 * p_trt + ns + 16) + 8 * nl * (p_or + 4) +
    (if (details) 3 * J + 8 * nl + 16 else 8) + if (is.null(backend$reference_ref)) 0 else 4 * ns
  metadata <- 16384 + 1024 * (ns + G + nt + nl + J) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  output <- B * (width + 16) + metadata
  scratch <- 64 * B * (width + p_or + p_trt + 8) + 128 * ((p_or + 2)^2 + (p_trt + 2)^2) + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- 128 * B * (
    (ns + 1) * (p_or + 2)^2 * log_a +
      (length(beta_gids) + length(gids)) * (p_or + 1) * (p_or + ns) +
      nt * ((p_trt + 2)^2 * log_a + p_trt * (ns + 4)) +
      J + nl * (p_or + p_trt + 4)
  ) + 64 * metadata
  list(
    gids = gids, beta_gids = beta_gids, sids = sids, tids = tids,
    output = output, scratch = scratch, work = work
  )
}

ge_old_ns <- function() {
  e <- new.env(parent = api)
  for (nm in c("evaluate_dr_baseline_if_chunk", "dr_baseline_if_budget", "replay_dr_baseline_if")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  old <- ge_old_envelope; environment(old) <- e
  e$dr_baseline_if_envelope <- old
  e
}

# Full frame IF (and detailed terms) of one backend through an evaluator
# namespace, assembled from a replay at `chunk_rows`.
ge_full_if <- function(fit, ns = api, chunk_rows = 7L, link_batch = 3L, scale = "frame", details = FALSE) {
  st <- ge_state(fit); N <- st$N; L <- length(fit$links)
  phi <- matrix(NA_real_, N, L); terms <- NULL; batches <- 0L
  if (details) terms <- lapply(c("psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw"), function(x) matrix(NA_real_, N, L))
  info <- ns$replay_dr_baseline_if(fit, chunk_rows, link_batch, function(b) {
    stopifnot(identical(b$id, st$id[b$rows]))
    phi[b$rows, b$links] <<- b$phi
    batches <<- batches + 1L
    if (details) for (k in seq_along(terms)) {
      nm <- c("psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw")[k]
      terms[[k]][b$rows, b$links] <<- b$details[[nm]]
    }
  }, scale = scale, details = details)
  list(phi = phi, terms = terms, batches = batches, info = info)
}

# ---------------------------------------------------------------------------
# 1. Hand-counted structure, exact envelope identities, affine coefficients
# ---------------------------------------------------------------------------

test_that("general DR IF envelope equals the hand-counted incidence/liveness bound on unequal-support fixtures", {
  fx <- ge_fixture(X = ge_design(41L, "x_bin")$X)
  fits <- list(
    x_bin = ge_fit(fx, ge_design(41L, "x_bin")),
    icpt_trt_x = ge_fit(fx, ge_design(41L, "icpt_trt_x"))
  )
  fx_i <- ge_fixture()
  fits$icpt <- ge_fit(fx_i, ge_design(41L, "icpt"))

  # (a) Designed catalog structure (literal expectations), identical for every design.
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    expect_null(fit$complete_support, info = nm)
    expect_identical(length(fit$support_specs), ge_structure$n_supports, info = nm)
    expect_identical(length(fit$regression_specs), ge_structure$n_blocks, info = nm)
    expect_identical(length(fit$treatment_specs), ge_structure$n_models, info = nm)
    expect_identical(sort(vapply(fit$treatment_specs, function(s) length(s$support_ids), numeric(1))),
      ge_structure$m_t_sorted, info = nm)
    expect_identical(sort(vapply(fit$regression_specs, function(s) length(s$support_ids), numeric(1))),
      ge_structure$block_sizes_sorted, info = nm)
    expect_identical(sum(vapply(fit$regression_specs, function(s) length(s$support_ids), numeric(1))),
      ge_structure$block_incidence, info = nm)
    expect_identical(length(fit$shared$control_support_ids), ge_structure$n_control_supports, info = nm)
    expect_identical(length(ge_state(fit)$cohorts), ge_structure$n_cohorts, info = nm)
    coh <- vapply(fit$support_specs, function(s) as.integer(s$cohort), integer(1))
    expect_identical(max(tabulate(coh)), ge_structure$cohort_support_max, info = nm)
    # Each model's support_ids are one per cohort, in the model's cohort order.
    for (s in fit$treatment_specs) expect_identical(length(s$support_ids), length(s$cohorts), info = nm)
    # The two (2,3) links share a target support and a control support.
    sup23 <- vapply(fit$support_specs, function(s) s$pair$s_calendar_index == 2L && s$pair$t_calendar_index == 3L, logical(1))
    expect_identical(sum(sup23), 3L, info = nm)
  }

  # (b) Literal hand count for x_bin (p_or = 1, p_trt = 2), all seven links.
  fit <- fits$x_bin
  expect_identical(fit$shared$p, 1L); expect_identical(fit$treatment_recipe$p, 2L)
  st <- ge_state(fit); frame <- ge_frame(fit)
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  id_max <- frame$id_byte_metadata$max
  link_bytes <- sum(nchar(vapply(seq_along(ge_links), function(i) ge_link_id(i, ge_links[[i]]), character(1)), type = "bytes"))
  expect_identical(link_bytes, sum(nchar(vapply(fit$links, `[[`, character(1), "link_id"), type = "bytes")))
  lit <- list()
  # width: 4 ns + 4 (p_or + 1)(nb + 1) + 4 ng + sum(3 p_trt + m_t + 16) + 8 nl (p_or + 4) + 8
  lit$width_nd <- 4 * 15 + 4 * 2 * 8 + 4 * 7 + (3 * 2 * 7 + 17 + 16 * 7) + 8 * 7 * 5 + 8
  expect_identical(lit$width_nd, 611)
  lit$width_d <- lit$width_nd - 8 + (3 * 3 + 8 * 7 + 16)
  expect_identical(lit$width_d, 684)
  # live, non-detailed: support 5 + 15 + 16; core 3*6 + 4*(10 + 10 + 7) + 8 + 9 + 14 + 2 + 14 + 20;
  # models 5*7 + (3 + 26 + 2 + 52); assembly 21
  lit$live_nd <- (5 + 15 + 16) + (3 * 6 + 4 * 27 + 8 + 9 + 14 + 2 + 14 + 20) + (35 + 3 + 26 + 2 + 52) + 21
  expect_identical(lit$live_nd, 368)
  # live, detailed: core + (7 + 21 + 5 + 8); models retained 171 + 83 + 32; assembly 21 + 35 + 6 + 6
  lit$live_d <- (5 + 15 + 16) + (193 + 7 + 21 + 5 + 8) + (171 + 83 + 32) + (21 + 35 + 6 + 6)
  expect_identical(lit$live_d, 624)
  lit$metadata0 <- 16384 + 1024 * (15 + 7 + 7 + 7 + 3) + 8 * link_bytes
  lit$units <- 16 * 9 * log_a + 2 * (14 + 20) + 7 * (16 * log_a + 1 + 8) + 12 * (17 + 28) + 3 + 7 * 7
  lit$lookups <- (27 + 7 + 7 + 1) * (15 + 7 + 1) + 2 * 27^2 + 7 * (15 + 3 + 1)
  expect_identical(lit$lookups, 2557)
  small <- 128 * (9 + 16)
  for (B in c(1, 2, 7)) {
    for (details in c(FALSE, TRUE)) {
      width <- if (details) lit$width_d else lit$width_nd
      live <- if (details) lit$live_d else lit$live_nd
      metadata <- lit$metadata0 + B * id_max
      pre <- api$dr_baseline_if_envelope(fit, B, seq_len(7L), details)
      expect_identical(pre$output, B * (width + 16) + metadata, info = paste("B", B, details))
      expect_identical(pre$scratch, 2 * B * live + small + metadata, info = paste("B", B, details))
      expect_identical(pre$work, 128 * B * lit$units + 64 * lit$lookups + 64 * metadata, info = paste("B", B, details))
      expect_identical(pre$general$live_entries_per_row, live)
      expect_identical(pre$general$row_units, lit$units)
      expect_identical(pre$general$lookup_elements, lit$lookups)
      expect_identical(pre$general$cohort_support_max, 6L)
      expect_identical(unname(pre$general$block_support_incidence), c(10, 10))
      expect_identical(pre$general$treatment_support_incidence, 17)
      expect_identical(pre$general$route, "general")
      expect_true(isTRUE(pre$general$or_core_evaluated))
      expect_identical(pre$general$reference_columns, 0L)
      # The catalog-derived hand count agrees with the literal one.
      hc <- ge_hand(fit, B, seq_len(7L), details)
      expect_identical(hc$width, width); expect_identical(hc$live, live)
      expect_identical(hc$units, lit$units); expect_identical(hc$lookups, lit$lookups)
    }
  }

  # (c) Exact identities across designs, link subsets, details and B.
  batches <- list(seq_len(7L), c(1L, 3L), c(3L, 7L), 6L, c(2L, 4L, 5L, 6L), 7:1)
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    for (links in batches) for (details in c(FALSE, TRUE)) for (B in c(1, 2, 7)) {
      info <- paste(nm, paste(links, collapse = ","), details, B)
      hc <- ge_hand(fit, B, links, details)
      pre <- api$dr_baseline_if_envelope(fit, B, links, details)
      expect_identical(pre$sids, hc$sids, info = info)
      expect_identical(pre$tids, hc$tids, info = info)
      expect_identical(pre$gids, hc$gids, info = info)
      expect_identical(pre$beta_gids, hc$beta_gids, info = info)
      expect_identical(pre$output, hc$output, info = info)
      expect_identical(pre$scratch, hc$scratch, info = info)
      expect_identical(pre$work, hc$work, info = info)
      expect_identical(pre$general$live_entries_per_row, hc$live, info = info)
      expect_identical(pre$general$cohort_support_max, hc$h_max, info = info)
      expect_identical(unname(pre$general$block_support_incidence), c(hc$I_beta, hc$I_g), info = info)
      expect_identical(pre$general$treatment_support_incidence, sum(hc$m_t), info = info)
      # p_or = 0 without details: no global beta blocks, yet the core and its hits are charged.
      if (fit$shared$p == 0L && !details) {
        expect_identical(hc$I_beta, 0); expect_length(pre$beta_gids, 0L)
        expect_gt(hc$I_g, 0)
      }
      # Unequal supports enter per model: the model term is sum(3 p_trt + m_t + 16), never nt * ns.
      nt <- length(hc$tids); ns <- hc$ns; p_trt <- fit$treatment_recipe$p
      expect_identical(sum(3 * p_trt + hc$m_t + 16) - nt * (2 * p_trt + ns + 16),
        nt * p_trt + sum(hc$m_t) - nt * ns, info = info)
      if (length(links) > 1L) expect_true(sum(hc$m_t) < nt * ns, info = info)
      expect_lte(max(hc$m_t), ns)
      if (length(links) == 7L) {
        expect_identical(sort(unique(hc$m_t)), c(2, 3), info = info)
        expect_lt(sum(3 * p_trt + hc$m_t + 16), nt * (2 * p_trt + ns + 16))
      }
    }
    # Affine row coefficients (what cp04_plan_rows reads from B = 1 and B = 2).
    for (links in batches) {
      hc <- ge_hand(fit, 1, links, FALSE)
      one <- api$dr_baseline_if_envelope(fit, 1, links)
      two <- api$dr_baseline_if_envelope(fit, 2, links)
      seven <- api$dr_baseline_if_envelope(fit, 7, links)
      slope <- two$output + two$scratch - one$output - one$scratch
      expect_identical(slope, hc$slope, info = nm)
      expect_identical(one$output + one$scratch - slope, hc$intercept, info = nm)
      expect_identical(two$output - one$output, hc$output_slope, info = nm)
      expect_identical(one$output - (two$output - one$output), hc$output_intercept, info = nm)
      expect_identical(seven$output + seven$scratch, hc$intercept + 7 * hc$slope, info = nm)
      expect_identical(seven$output, hc$output_intercept + 7 * hc$output_slope, info = nm)
      expect_identical(seven$work - one$work, 6 * (128 * hc$units + 64 * hc$id_max), info = nm)
    }
  }

  # (d) Against the frozen previous generic branch: same sets, strictly smaller
  # concurrent charge and no larger output on these fixtures.
  old_ns <- ge_old_ns()
  for (nm in names(fits)) for (links in batches) for (details in c(FALSE, TRUE)) {
    new <- api$dr_baseline_if_envelope(fits[[nm]], 7, links, details)
    old <- old_ns$dr_baseline_if_envelope(fits[[nm]], 7, links, details)
    expect_identical(new$sids, old$sids); expect_identical(new$tids, old$tids)
    expect_identical(new$gids, old$gids); expect_identical(new$beta_gids, old$beta_gids)
    expect_lt(new$scratch, old$scratch)
    hc <- ge_hand(fits[[nm]], 7, links, details)
    expect_identical(new$output - old$output,
      7 * (length(hc$tids) * fits[[nm]]$treatment_recipe$p + sum(hc$m_t) - length(hc$tids) * hc$ns))
    if (identical(links, seq_len(7L))) expect_lt(new$output, old$output)
    expect_null(old$general)
  }
})

# ---------------------------------------------------------------------------
# 2. Unchanged automatic planner: adjacent thresholds, refusal, explicit controls
# ---------------------------------------------------------------------------

test_that("cp04_plan_rows admits the general DR row batch from the proved affine bound at adjacent thresholds", {
  des <- ge_design(41L, "x_bin")
  fx <- ge_fixture(X = des$X, ctrl = ge_ctrl(retained = 2e7))
  fit_auto <- ge_fit(fx, des, automatic = TRUE)
  fit_expl <- ge_fit(fx, des, automatic = FALSE)
  expect_true(isTRUE(fit_auto$accounting$automatic_chunk_rows))
  expect_false(isTRUE(fit_expl$accounting$automatic_chunk_rows))
  M <- fx$ctrl$max_matrix_entries
  N <- ge_state(fit_auto)$N
  expect_identical(N, 41L)
  links <- seq_len(7L)
  own <- fit_auto$accounting$owned_entries + fit_auto$accounting$borrowed_entries
  hc <- ge_hand(fit_auto, 1, links, FALSE)

  # Threshold: retained R5 leaves exactly 6 slope - 1 entries above the intercept -> B = 5.
  R5 <- M - own - hc$intercept - 6 * hc$slope + 1
  expect_gt(R5, 0)
  expect_identical(api$cp04_plan_rows(fit_auto, list(links), R5), 5)
  expect_identical(api$cp04_plan_rows(fit_auto, list(links), R5 - 1), 6)
  # Disclosed test-local one-entry undercharge: the independent threshold
  # predicate changes from B=5 to B=6 and therefore detects this mutation.
  bad <- new.env(parent = api)
  bad$dr_baseline_if_envelope <- function(...) {
    pre <- api$dr_baseline_if_envelope(...); pre$scratch <- pre$scratch - 1; pre
  }
  planner_mutant <- api$cp04_plan_rows; environment(planner_mutant) <- bad
  expect_identical(planner_mutant(fit_auto, list(links), R5), 6)
  expect_false(identical(planner_mutant(fit_auto, list(links), R5), 5))

  expect_identical(api$cp04_plan_rows(fit_auto, list(links), R5 + hc$slope), 4)
  # Several link sets: the smallest admitted batch over the sets.
  hc_b <- ge_hand(fit_auto, 1, c(1L, 3L), FALSE)
  expect_lt(hc_b$slope, hc$slope)
  expect_identical(api$cp04_plan_rows(fit_auto, list(c(1L, 3L), links), R5), 5)
  expect_identical(api$cp04_plan_rows(fit_auto, list(c(1L, 3L)), R5),
    min(N, fx$ctrl$max_chunk_rows, floor((M - R5 - own - hc_b$intercept) / hc_b$slope)))
  # Ample budget: the outer explicit row cap, never more.
  expect_identical(api$cp04_plan_rows(fit_auto, list(links), 0), 7)
  # One admitted row at the exact one-row charge; one entry less refuses with the typed bound error.
  one <- api$dr_baseline_if_envelope(fit_auto, 1, links)
  expect_identical(one$output + one$scratch, hc$intercept + hc$slope)
  R1 <- M - own - one$output - one$scratch
  expect_identical(api$cp04_plan_rows(fit_auto, list(links), R1), 1)
  expect_error(api$cp04_plan_rows(fit_auto, list(links), R1 + 1), class = "cdid_ipw_owner_bound_error")
  # Explicit controls stay respected: an explicit chunk is returned unchanged
  # (no planning), one above the cap refuses, a non-automatic fit never plans.
  expect_identical(api$cp04_plan_rows(fit_auto, list(links), R5, chunk_rows = 3L), 3)
  expect_error(api$cp04_plan_rows(fit_auto, list(links), R5, chunk_rows = 8L), class = "cdid_ipw_owner_bound_error")
  expect_identical(api$cp04_plan_rows(fit_expl, list(links), R5), 7)
  expect_identical(api$cp04_plan_rows(fit_expl, list(links), R1 + 1), 7)

  # The streamed covariance plans its rows through the same unchanged planner
  # (retained + 4 L L): the saved chunk_rows is the admitted B, and the
  # covariance does not depend on it.
  L <- length(fit_auto$links)
  cov_5 <- api$stream_dr_covariance(fit_auto, retained = R5 - 4 * L * L)
  cov_7 <- api$stream_dr_covariance(fit_expl)
  cov_a <- api$stream_dr_covariance(fit_auto)
  expect_identical(cov_5$chunk_rows, 5)
  expect_identical(cov_7$chunk_rows, 7)
  expect_identical(cov_a$chunk_rows, 7)
  expect_identical(cov_5$estimates, cov_7$estimates)
  expect_equal(cov_5$covariance, cov_7$covariance, tolerance = 1e-12)
  expect_equal(cov_5$mean_phi, cov_7$mean_phi, tolerance = 1e-12)
  expect_equal(cov_a$covariance, cov_7$covariance, tolerance = 1e-12)
  expect_true(all(is.finite(cov_5$covariance)) && max(abs(cov_5$covariance)) > 0)
  expect_error(api$stream_dr_covariance(fit_auto, retained = M), class = "cdid_ipw_owner_bound_error")
})

# ---------------------------------------------------------------------------
# 3. Actual evaluated payloads against the bound; work-cap adjacency
# ---------------------------------------------------------------------------

test_that("general DR IF chunks stay within the corrected output bound and keep every term; caps refuse before evaluation", {
  cases <- list(
    x_bin = list(fx = ge_fixture(X = ge_design(41L, "x_bin")$X), des = ge_design(41L, "x_bin")),
    icpt_trt_x = list(fx = ge_fixture(X = ge_design(41L, "icpt_trt_x")$X), des = ge_design(41L, "icpt_trt_x")),
    icpt = list(fx = ge_fixture(), des = ge_design(41L, "icpt"))
  )
  old_ns <- ge_old_ns()
  for (nm in names(cases)) {
    fit <- ge_fit(cases[[nm]]$fx, cases[[nm]]$des)
    st <- ge_state(fit); p_trt <- fit$treatment_recipe$p; p_or <- fit$shared$p
    for (links in list(seq_len(7L), c(1L, 3L), 7L)) for (details in c(FALSE, TRUE)) for (scale in c("frame", "core")) {
      for (rows in list(1:7, 8:14, 36:41, 5L)) {
        B <- length(rows); info <- paste(nm, paste(range(links), collapse = "-"), details, scale, B)
        hc <- ge_hand(fit, B, links, details)
        out <- api$evaluate_dr_baseline_if_chunk(fit, rows, links, scale, details)
        expect_identical(out$diagnostics$output_entries_bound, hc$output, info = info)
        expect_identical(out$diagnostics$concurrent_scratch_bound, hc$scratch, info = info)
        expect_identical(out$diagnostics$work_ops_bound, hc$work, info = info)
        expect_lte(unname(api$cp04_payload(out)["entries"]), hc$output, label = info)
        # D90 skips the unused p0 core only for non-detailed, non-reference chunks.
        # These fixtures have no reference reader; conservative allocation bounds stay unchanged.
        expected_core <- details || p_or > 0L
        expect_identical(out$diagnostics$n_regression_score_evaluations,
          if (expected_core) length(union(hc$gids, hc$beta_gids)) else 0L, info = info)
        expect_identical(out$diagnostics$n_global_beta_blocks, length(hc$beta_gids), info = info)
        expect_identical(out$diagnostics$n_treatment_score_evaluations, length(hc$tids), info = info)
        expect_identical(out$diagnostics$n_rate_score_evaluations, hc$ns, info = info)
        expect_identical(dim(out$phi), c(B, length(links)), info = info)
        expect_true(all(is.finite(out$phi)), info = info)
        if (details) {
          d <- out$details
          expect_identical(dim(d$rates), c(B, hc$ns), info = info)
          expect_identical(d$treatment_ids, hc$tids, info = info)
          expect_identical(length(d$models), length(fit$treatments), info = info)
          for (k in seq_along(fit$treatments)) {
            m <- d$models[[k]]
            if (!k %in% hc$tids) { expect_null(m, info = info); next }
            m_k <- length(fit$treatment_specs[[k]]$support_ids)
            expect_identical(dim(m$upstream), c(B, m_k), info = info)
            for (f in c("beta", "gamma", "beta_score")) expect_identical(dim(m[[f]]), c(B, p_trt), info = paste(info, f))
            expect_identical(dim(m$coordinates), c(B, 4L), info = info)
            expect_identical(dim(m$direct), c(B, 4L), info = info)
            expect_identical(dim(m$normalizers_A), c(B, 2L), info = info)
            for (f in c("psi_dr", "psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw"))
              expect_identical(length(m[[f]]), B, info = paste(info, f))
            expect_identical(sort(names(m)), sort(c("beta", "gamma", "coordinates", "direct", "beta_score", "upstream",
              "psi_dr", "psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw", "normalizers_A")), info = info)
            # Detailed payload per model is exactly the charged (3 p_trt + m_t + 16) B-scaled entries.
            expect_identical(sum(vapply(m, function(v) as.double(length(v)), numeric(1))), B * (3 * p_trt + m_k + 16), info = info)
          }
          for (f in c("psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw", "psi_or_estimator", "direct", "denominator", "numerator"))
            expect_identical(dim(d[[f]]), c(B, length(links)), info = paste(info, f))
          expect_identical(dim(d$tau), c(B, hc$J), info = info)
          expect_identical(dim(d$beta), c(B, p_or), info = info)
          expect_identical(dim(d$alpha), c(B, length(hc$gids)), info = info)
          expect_identical(length(d$beta_block_scores), length(hc$beta_gids), info = info)
        } else {
          expect_null(out$details, info = info)
        }
        # Same arithmetic as the frozen previous envelope's evaluator: bitwise identical IF and details.
        old <- old_ns$evaluate_dr_baseline_if_chunk(fit, rows, links, scale, details)
        expect_identical(out$phi, old$phi, info = info)
        expect_identical(out$details, old$details, info = info)
        expect_identical(out$rows, old$rows); expect_identical(out$id, old$id)
        expect_lt(out$diagnostics$concurrent_scratch_bound, old$diagnostics$concurrent_scratch_bound, label = info)
      }
    }
    # Work cap adjacency: the chunk's cumulative charge is validated before any
    # row is read; one op less refuses with the typed bound error.
    links <- seq_len(7L); rows <- 1:7
    pre <- api$dr_baseline_if_envelope(fit, 7, links, FALSE)
    W_fit <- fit$accounting$work_ops_bound
    fit_edge <- ge_fit(cases[[nm]]$fx, cases[[nm]]$des, ctrl = ge_ctrl(7L, work = W_fit + pre$work))
    expect_identical(fit_edge$accounting$work_ops_bound, W_fit)
    out_edge <- api$evaluate_dr_baseline_if_chunk(fit_edge, rows, links)
    expect_identical(out_edge$phi, api$evaluate_dr_baseline_if_chunk(fit, rows, links)$phi)
    fit_fail <- ge_fit(cases[[nm]]$fx, cases[[nm]]$des, ctrl = ge_ctrl(7L, work = W_fit + pre$work - 1))
    expect_error(api$evaluate_dr_baseline_if_chunk(fit_fail, rows, links), class = "cdid_ipw_owner_bound_error")
    # The whole replay charges ceiling(N / B) chunks of the same bound.
    expect_error(api$replay_dr_baseline_if(fit_edge, 7L, 7L, function(b) NULL), class = "cdid_ipw_owner_bound_error")
    # Invalid batches are still refused up front.
    expect_error(api$evaluate_dr_baseline_if_chunk(fit, 1:8, links), class = "cdid_dr_baseline_input_error")
    expect_error(api$evaluate_dr_baseline_if_chunk(fit, 1:7, 8L), class = "cdid_dr_baseline_input_error")
  }
})

# ---------------------------------------------------------------------------
# 4. D60 reference binding: control and negatives
# ---------------------------------------------------------------------------

test_that("the general envelope charges reference reads only through the sealed reader; forged bindings refuse", {
  des <- ge_design(41L, "x_bin")
  fx <- ge_fixture(X = des$X)
  fit <- ge_fit(fx, des, automatic = TRUE)
  links <- seq_len(7L)
  # Control: no reference reader on these fixtures, no reference term in the charge.
  expect_null(fit$reference_ref)
  expect_null(api$chained_backend_reference(fit))
  pre <- api$dr_baseline_if_envelope(fit, 7, links)
  expect_identical(pre$general$reference_columns, 0L)
  hc <- ge_hand(fit, 7, links, FALSE)
  expect_identical(pre$scratch, hc$scratch)
  # Negatives: a binding without a sealed reader, or with an unsealed fake
  # reader carrying a `x_columns` width, is refused by the envelope itself
  # (never charged from the forged width) and by the evaluator's validation.
  forged_env <- unclass(fit); attributes(forged_env) <- attributes(fit)
  forged_env$reference_ref <- new.env(parent = emptyenv())
  expect_error(api$dr_baseline_if_envelope(forged_env, 7, links), class = "cdid_reference_integrity_error")
  expect_error(api$evaluate_dr_baseline_if_chunk(forged_env, 1:7, links), class = "cdid_dr_baseline_integrity_error")
  forged_fake <- forged_env
  assign("reference_reader", list(x_columns = c("(Intercept)", "x", "z1", "z2"), N = 41L), envir = forged_fake$reference_ref)
  expect_error(api$dr_baseline_if_envelope(forged_fake, 7, links), class = "cdid_reference_integrity_error")
  expect_error(api$evaluate_dr_baseline_if_chunk(forged_fake, 1:7, links), class = "cdid_dr_baseline_integrity_error")
  expect_error(api$cp04_plan_rows(forged_fake, list(links), 0), class = "cdid_reference_integrity_error")
  expect_identical(api$dr_baseline_if_envelope(fit, 7, links)$scratch, hc$scratch)
})

# ---------------------------------------------------------------------------
# 5. Chunk-partition invariance against the frozen previous evaluator path
# ---------------------------------------------------------------------------

test_that("full general DR IF, terms and streamed covariance are invariant to the admitted chunk size and match the frozen evaluator path", {
  des <- ge_design(41L, "x_bin")
  fx7 <- ge_fixture(X = des$X, chunk = 7L)
  fx41 <- ge_fixture(X = des$X, chunk = 41L)
  fit7 <- ge_fit(fx7, des)
  fit41 <- ge_fit(fx41, des)
  old_ns <- ge_old_ns()
  N <- ge_state(fit7)$N
  ref <- ge_full_if(fit7, api, chunk_rows = 1L, link_batch = 1L, details = TRUE)
  expect_identical(ref$batches, N * 7L)
  expect_true(all(is.finite(ref$phi)))
  for (cfg in list(list(fit = fit7, B = 7L, LB = 3L), list(fit = fit7, B = 5L, LB = 7L),
                   list(fit = fit41, B = 41L, LB = 2L), list(fit = fit41, B = 13L, LB = 4L))) {
    new <- ge_full_if(cfg$fit, api, chunk_rows = cfg$B, link_batch = cfg$LB, details = TRUE)
    old <- ge_full_if(cfg$fit, old_ns, chunk_rows = cfg$B, link_batch = cfg$LB, details = TRUE)
    expect_identical(new$batches, as.integer(ceiling(N / cfg$B) * ceiling(7L / cfg$LB)))
    expect_identical(new$phi, old$phi)
    expect_identical(new$terms, old$terms)
    expect_equal(new$phi, ref$phi, tolerance = 1e-12)
    for (k in seq_along(ref$terms)) expect_equal(new$terms[[k]], ref$terms[[k]], tolerance = 1e-12)
    expect_identical(new$info$all_F_identities, TRUE)
    # The replay's cumulative work bound is the chunk count times the corrected chunk bound.
    pre <- api$dr_baseline_if_envelope(cfg$fit, min(cfg$B, N), seq_len(min(cfg$LB, 7L)), TRUE)
    expect_lte(pre$work * ceiling(N / cfg$B), new$info$work_ops_bound)
    hand_work <- cfg$fit$accounting$work_ops_bound
    for (start in seq.int(1L, 7L, by = cfg$LB)) {
      ids <- seq.int(start, min(7L, start + cfg$LB - 1L))
      hand_work <- hand_work + ceiling(N / cfg$B) * ge_hand(cfg$fit, min(cfg$B, N), ids, TRUE)$work
    }
    expect_identical(new$info$work_ops_bound, hand_work)
  }
  core <- ge_full_if(fit7, api, chunk_rows = 7L, link_batch = 7L, scale = "core", details = FALSE)
  expect_equal(core$phi, ge_state(fit7)$lambda * ref$phi, tolerance = 1e-12)
  cov7 <- api$stream_dr_covariance(fit7)
  cov41 <- api$stream_dr_covariance(fit41)
  # The genuine accepted OLD library exhibits this same ~2e-16 chunk-fit
  # rounding; use the established scientific chunk-invariance criterion.
  expect_equal(cov7$estimates, cov41$estimates, tolerance = 1e-12)
  expect_equal(cov7$covariance, cov41$covariance, tolerance = 1e-12)
  expect_identical(cov7$chunk_rows, 7); expect_identical(cov41$chunk_rows, 41)
  # Independent assembly of the covariance from the core-scale full IF.
  centered <- sweep(core$phi, 2L, colMeans(core$phi), "-")
  expect_equal(unname(crossprod(centered) / ge_state(fit7)$n_A), unname(cov7$covariance), tolerance = 1e-12)
  expect_identical(unname(cov7$estimates), vapply(fit7$links, `[[`, numeric(1), "mu_dr"))
})

# ---------------------------------------------------------------------------
# 6. Public DR call on the same missing panel: explicit versus automatic rows
# ---------------------------------------------------------------------------

test_that("public chained DR on a missing panel matches between explicit chunk rows and the automatically admitted batch", {
  skip_if_not_installed("did")
  d <- ge_panel(cohort4 = FALSE)
  d$g[d$g == Inf] <- 0
  run <- function(control) cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g", data = d, panel = TRUE,
    xformla = ~x, est_method = "dr", control_group = "nevertreated", bstrap = FALSE, cband = FALSE,
    base_period = "varying", mode = "chained", pt_assumption = "post", control = control
  )
  mp_auto <- run(NULL)
  mp_expl <- run(list(chunk_rows = 7L))
  for (mp in list(mp_auto, mp_expl)) {
    expect_s3_class(mp, "cdid_MP")
    expect_true(all(is.finite(mp$att)) && all(is.finite(mp$se)))
    expect_identical(mp$chained$estimator, "chained doubly robust GMM")
    expect_identical(mp$DIDparams$est_method, "dr")
  }
  acc_auto <- mp_auto$chained$resource_accounting$point_producer
  acc_expl <- mp_expl$chained$resource_accounting$point_producer
  n_frame <- length(unique(d$id))
  expect_identical(acc_expl$chunk_rows, 7)
  # Ample default budgets: the planner admits the public row cap (here the frame size) from the proved bound.
  expect_identical(as.integer(acc_auto$chunk_rows), as.integer(min(n_frame, mp_auto$chained$public$chunk_rows)))
  expect_lt(acc_auto$producer_workspace_entries, mp_auto$chained$public$numeric_workspace_entries)
  expect_equal(mp_auto$att, mp_expl$att, tolerance = 1e-12)
  expect_equal(mp_auto$se, mp_expl$se, tolerance = 1e-12)
  expect_equal(mp_auto$W, mp_expl$W, tolerance = 1e-12)
  expect_equal(mp_auto$Wpval, mp_expl$Wpval, tolerance = 1e-12)
  if (!is.null(mp_auto$inffunc)) expect_equal(as.matrix(mp_auto$inffunc), as.matrix(mp_expl$inffunc), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 7. Prospective arithmetic on the saved half-missing DR request (no fit, no claim)
# ---------------------------------------------------------------------------

test_that("the saved half-missing DR planner reconstruction is exact under the previous charge and admits the cap under the corrected charge (arithmetic only)", {
  # Saved request facts (result.json of h01_half_dr_residual): N = n_A = 1000,
  # 190 links/models/blocks, 380 supports, p_or = 1, p_trt = 2, J = 2,
  # 8-byte ids, 2,914 link-id bytes, own 5,022,048, workspace 134,217,728,
  # public row cap 171, admitted B = 18, producer entries 109,344,304 and
  # work 69,177,092,096 over 56 chunks. Every treatment model has two
  # supports; each block one control support; 190 supports share each cohort.
  N <- 1000; nl <- nt <- nb <- ng <- G <- 190; ns <- 380; J <- 2; p_or <- 1; p_trt <- 2; id_max <- 8; lb <- 2914
  own <- 5022048; M <- 134217728; cap <- 171
  log_a <- 1 + ceiling(log2(N + 1))
  meta <- function(B) 16384 + 1024 * (ns + G + nt + nl + J) + B * id_max + 8 * lb
  # Previous generic charge.
  w_old <- 4 * ns + 4 * (p_or + 1) * (nb + 1) + 4 * ng + nt * (2 * p_trt + ns + 16) + 8 * nl * (p_or + 4) + 8
  out_old <- function(B) B * (w_old + 16) + meta(B)
  scr_old <- function(B) 64 * B * (w_old + p_or + p_trt + 8) + 128 * ((p_or + 2)^2 + (p_trt + 2)^2) + meta(B)
  work_old <- function(B) 128 * B * ((ns + 1) * (p_or + 2)^2 * log_a + (nb + ng) * (p_or + 1) * (p_or + ns) +
    nt * ((p_trt + 2)^2 * log_a + p_trt * (ns + 4)) + J + nl * (p_or + p_trt + 4)) + 64 * meta(B)
  expect_identical(w_old, 87416)
  expect_identical(own + out_old(18) + scr_old(18), 109344304)
  expect_identical(56 * work_old(18), 69177092096)
  slope_old <- out_old(2) + scr_old(2) - out_old(1) - scr_old(1)
  icpt_old <- out_old(1) + scr_old(1) - slope_old
  expect_identical(slope_old, 5682776)
  # Corrected generic charge (the hand inventory above with m_t = 2, I = 190, h_max = 190, n_ctrl = 190, c_max = 1).
  m_t <- rep(2, nt); m_max <- 2; I_beta <- I_g <- 190; c_max <- 1; n_ctrl <- 190; h_max <- 190
  w_new <- 4 * ns + 4 * (p_or + 1) * (nb + 1) + 4 * ng + sum(3 * p_trt + m_t + 16) + 8 * nl * (p_or + 4) + 8
  live <- (5 + ns + 16) +
    ((2 + p_or) * h_max + (3 + p_or) * (I_beta + I_g + nl) + 8 * p_or + n_ctrl + 2 * ng + c_max + 2 * nl + 4 * (p_or + 4)) +
    (5 * nt + m_max + 13 * p_trt + 2 * p_or + 52) + 3 * nl
  expect_identical(w_new, 15976); expect_identical(live, 5832)
  out_new <- function(B) B * (w_new + 16) + meta(B)
  scr_new <- function(B) 2 * B * live + 128 * ((p_or + 2)^2 + (p_trt + 2)^2) + meta(B)
  slope_new <- out_new(2) + scr_new(2) - out_new(1) - scr_new(1)
  icpt_new <- out_new(1) + scr_new(1) - slope_new
  expect_identical(slope_new, w_new + 16 + 2 * live + 2 * id_max)
  expect_identical(icpt_new, icpt_old)
  plan <- function(icpt, slope, out_fn, retained) min(N, cap, floor((M - retained - own - icpt) / slope),
    floor((M - meta(0)) / (out_fn(2) - out_fn(1))))
  # Every pipeline retention consistent with the saved B = 18 admits the row cap under the corrected charge.
  lo <- M - own - icpt_old - 19 * slope_old + 1; hi <- M - own - icpt_old - 18 * slope_old
  expect_true(lo <= hi)
  for (retained in c(lo, hi)) {
    expect_identical(plan(icpt_old, slope_old, out_old, retained), 18)
    expect_identical(plan(icpt_new, slope_new, out_new, retained), cap)
    expect_lt(own + out_new(cap) + scr_new(cap), own + out_old(18) + scr_old(18))
  }
})
