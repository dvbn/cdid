# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 balanced-panel performance, stage 4: native complete-support batches
# (fit-time support preparation and the non-detailed OR/DR IF replay), the
# exact incremental catalog metadata charge and the two zero-variable
# provenance corrections. Expected values come from the unchanged general
# route of the same build, the retained R complete-support functions, the
# recursive cp04_payload walk and raw panel rows; the native kernels never
# supply their own truth.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

nb_ctrl <- function(chunk = 7L, work = 1e10, matrix = 2e7) api$ipw_owner_control(max_chunk_rows = chunk,
  max_retained_entries = 1e7, max_matrix_entries = matrix, max_work_ops = work, max_output_entries = 1e7)
nb_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
nb_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))
nb_kernels <- c("C_gmm_cp08_support_chunk", "C_gmm_cp08_pair_rows", "C_gmm_cp08_fit_pass1",
  "C_gmm_cp08_fit_pass2", "C_gmm_cp08_or_if_chunk", "C_gmm_cp08_dr_if_prepare", "C_gmm_cp08_dr_if_assemble", "C_gmm_cp08_dr_treatment_batch")

# Balanced long panel: g = 3 treated, Inf never treated (optionally a third,
# unreferenced cohort 4), 41 units, 5 dates.
nb_panel <- function(N = 41L, T = 5L, seed = 11L, extra = FALSE) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:2] <- c(3, Inf)
  if (extra) g[3:6] <- 4
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$y <- 0.2 * d$id + 0.4 * d$t + ifelse(d$t >= d$g, 1.25, 0) + stats::rnorm(nrow(d))
  d
}
nb_pairs <- function(T = 5L) {
  pr <- list()
  for (s in seq_len(T - 1L)) for (t in seq.int(s + 1L, T)) pr[[length(pr) + 1L]] <- c(s, t)
  pr
}
# Real sealed owners, the full-roster reader and link views; link IDs carry the
# link position so repeated pairs are distinct links of one block. `wt` is the
# constant unit weight; dyadic values keep the existing link-view rate sums
# exact while their products still round.
nb_fixture <- function(d = nb_panel(), chunk = 7L, perm = FALSE, pairs = nb_pairs(), char_ids = FALSE, wt = 1,
                       X = NULL, V = NULL, ctrl = nb_ctrl(chunk)) {
  if (char_ids) d$id <- sprintf("unit-%03d", d$id)
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  u <- d[!duplicated(d$id), c("id", "g")]
  u <- u[order(u$id), ]
  if (perm) { set.seed(5L); u <- u[sample(nrow(u)), ] }
  N <- nrow(u)
  ord <- match(u$id, sort(u$id))
  Z <- if (is.null(X)) matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)")) else X[ord, , drop = FALSE]
  V_mat <- if (is.null(V)) matrix(1, N, 1L, dimnames = list(NULL, "reader_intercept")) else V[ord, , drop = FALSE]
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = rep(wt, N), control = ctrl, frame_g = u$g,
    allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, u$id, V_mat,
    u$g, rep(wt, N), ctrl, allow_unknown_off_A = TRUE)
  views <- lapply(seq_along(pairs), function(k) {
    p <- pairs[[k]]
    api$construct_ipw_link_view(fo, p[1], p[2], g = 3, K = Inf, control = ctrl,
      link_id = paste0("L", k, ":", p[1], "_", p[2]), functional_key = paste0("cmpcat:nb:", k, ":", p[1], "-", p[2]))
  })
  list(d = d, po = po, fo = fo, reader = reader, views = views, ctrl = ctrl, chunk = chunk, u = u, pairs = pairs)
}
nb_general <- local({
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$or_baseline_complete_facts <- function(...) NULL
  e
})
nb_or <- function(fx, ns = api, ctrl = fx$ctrl, recipe = nb_icpt) ns$fit_or_baseline(fx$reader, fx$views, ctrl, 1e-10,
  outcome_recipe = recipe)
nb_dr <- function(fx, ns = api, ctrl = fx$ctrl, recipe = nb_icpt, treatment_columns = NULL)
  ns$fit_dr_baseline(fx$reader, fx$views, ctrl, glm_control = nb_glm,
    rank_tol = 1e-10, outcome_recipe = recipe, treatment_columns = treatment_columns)
nb_fit <- function(method, fx, ns = api, ctrl = fx$ctrl, recipe = nb_icpt, treatment_columns = NULL)
  if (method == "or") nb_or(fx, ns, ctrl, recipe) else nb_dr(fx, ns, ctrl, recipe, treatment_columns)
nb_chunk <- function(method, fit, ...) if (method == "or") api$evaluate_or_baseline_if_chunk(fit, ...) else
  api$evaluate_dr_baseline_if_chunk(fit, ...)
nb_binding <- function(method) if (method == "or") "or_baseline_backend" else "dr_baseline_backend"
nb_err <- function(method, cls) paste0(if (method == "or") "cdid_or_baseline_" else "cdid_dr_baseline_", cls, "_error")
nb_reseal <- function(x, binding, edit) { x$ref <- NULL; x <- edit(x); api$ipw_owner_seal(x, binding) }
nb_strip <- function(b) { b$diagnostics <- b$diagnostics[setdiff(names(b$diagnostics),
  c("output_entries_bound", "concurrent_scratch_bound", "work_ops_bound", "n_regression_score_evaluations"))]; b }
nb_expect_match <- function(a, b, tol = 1e-13, info = NULL) {
  # Same strict schema/leaf predicates, one reported assertion. Repeated waldo
  # assertions per scalar rebuild test metadata and obscure these whole-object checks.
  mismatch <- function(x, y, path) {
    if (!identical(class(x), class(y)) || !identical(typeof(x), typeof(y)) ||
        !identical(names(x), names(y)) || !identical(dim(x), dim(y)) ||
        !identical(dimnames(x), dimnames(y))) return(paste(path, "schema"))
    if (is.list(x)) {
      if (length(x) != length(y)) return(paste(path, "length"))
      return(unlist(lapply(seq_along(x), function(i) mismatch(x[[i]], y[[i]],
        paste0(path, "$", if (is.null(names(x))) i else names(x)[i]))), use.names = FALSE))
    }
    same <- if (typeof(x) == "double") length(waldo::compare(x, y, tolerance = tol)) == 0L else identical(x, y)
    if (same) character() else paste(path, "value")
  }
  bad <- mismatch(a, b, "root")
  expect_true(length(bad) == 0L, info = paste(c(info, head(bad, 12L)), collapse = "\n"))
}
# Counts every .Call of a CP08 kernel made by the re-bound native adapters.
nb_counted <- function() {
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline", "or_baseline_complete_supports", "or_baseline_complete_state",
               "or_pooled_fit_common", "or_pooled_wls_block_summary", "or_pooled_wls_merge_chunks",
               "or_pooled_wls_partial_block", "or_wls_summary", "or_wls_merge", "dr_baseline_treatment_fit", "dr_baseline_treatment_batch",
               "replay_or_baseline_if", "replay_dr_baseline_if", "evaluate_or_baseline_if_chunk",
               "evaluate_dr_baseline_if_chunk", "or_baseline_native_if_chunk", "dr_baseline_native_if_chunk")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  n <- new.env(parent = emptyenv())
  n$calls <- setNames(integer(length(nb_kernels)), nb_kernels); n$reads <- 0; n$slots <- 0; n$compress <- 0L
  n$pair_rows <- 0L; n$solve <- 0L
  addr <- lapply(nb_kernels, api$or_baseline_native_kernel)
  e$.Call <- function(.NAME, ...) {
    k <- which(vapply(addr, identical, logical(1), .NAME))
    res <- base::.Call(.NAME, ...)
    if (length(k)) {
      n$calls[k] <- n$calls[k] + 1L
      if (!is.null(res$counts) && nb_kernels[k] != "C_gmm_cp08_dr_if_assemble") {
        n$slots <- n$slots + res$counts[1L]; n$reads <- n$reads + res$counts[2L]
      }
    }
    res
  }
  e$or_wls_compress <- function(...) { n$compress <- n$compress + 1L; api$or_wls_compress(...) }
  e$or_baseline_complete_pair_rows <- function(...) { n$pair_rows <- n$pair_rows + 1L; api$or_baseline_complete_pair_rows(...) }
  e$cp04_solve <- function(...) { n$solve <- n$solve + 1L; api$cp04_solve(...) }
  e$n <- n
  e
}

test_that("incremental catalog metadata charges equal cp04_payload on every prefix and at the refusal edge", {
  recs <- list(
    list(pair = list(s = 1, t = 2, s_calendar_index = 1L, t_calendar_index = 2L), target = 2L, controls = 1L,
      link_id = "L1_2", source_functional_key = "cmpcat:a", functional_key = "cmpcat:OR|a"),
    list(pair = structure(list(s = 1, t = 3), class = "x_pair", note = "attribute"), target = 2L, controls = c(1L, 3L),
      link_id = "Lé_3", source_functional_key = "cmpcat:b", functional_key = NA_character_, reference = 2L),
    list(),
    list(x = c(a = 1, b = 2), y = character(0), z = NULL, w = list(NULL, "qß")))
  for (named in c(FALSE, TRUE)) for (abort in list(api$or_baseline_abort, api$dr_baseline_abort)) {
    records <- list(); st <- list(n = 0L, records = c(entries = 0, char_bytes = 0))
    for (k in seq_along(recs)) {
      if (named) records[[paste0("r", k)]] <- recs[[k]] else records[[k]] <- recs[[k]]
      st <- api$or_baseline_payload_append(st, records, abort)
      expect_identical(st$total, api$cp04_payload(records))
    }
    expect_identical(api$or_baseline_payload_verify(st, records, abort), api$cp04_payload(records))
  }
  # Meaningful negative witnesses: a skipped append and a wrong aggregate refuse.
  st0 <- list(n = 0L, records = c(entries = 0, char_bytes = 0))
  expect_error(api$or_baseline_payload_append(st0, records), class = "cdid_or_baseline_integrity_error")
  bad <- st; bad$total[["char_bytes"]] <- bad$total[["char_bytes"]] + 2
  expect_error(api$or_baseline_payload_verify(bad, records, api$dr_baseline_abort), class = "cdid_dr_baseline_integrity_error")
  # The last catalog append's work charge, recomputed from the recursive walk of
  # the records the fit builds, is admitted exactly and refused one below.
  fx <- nb_fixture(); st <- fx$reader$state; L <- length(fx$views); J <- length(st$cohorts)
  for (method in c("or", "dr")) {
    tag <- if (method == "or") "OR-baseline-shared-slope-WLS|" else "DR-baseline-shared-slope-WLS-IPW|"
    records <- lapply(fx$views, function(v) list(pair = v$pair, target = match(v$cohorts$g, st$cohorts),
      controls = match(sort(v$cohorts$K, method = "radix"), st$cohorts), link_id = v$link_id,
      source_functional_key = v$functional_key, functional_key = paste0("cmpcat:", tag, substring(v$functional_key, 8L))))
    W <- 128 * L * L * (J + 1) + 64 * unname(api$cp04_payload(records)["char_bytes"])
    label <- if (method == "or") "baseline catalog metadata work charge" else "baseline DR catalog metadata work charge"
    e <- tryCatch(nb_fit(method, fx, ctrl = nb_ctrl(7L, work = W - 1)), error = identity)
    expect_s3_class(e, "cdid_ipw_owner_bound_error")
    expect_identical(e$diagnostics$label, label)
    expect_identical(e$diagnostics$requested, W)
    e <- tryCatch(nb_fit(method, fx, ctrl = nb_ctrl(7L, work = W)), error = identity)
    expect_false(identical(e$diagnostics$label, label))
  }
})

test_that("native fit-time supports equal the retained R complete-support and general supports bit for bit", {
  for (case in list(list(chunk = 7L, perm = FALSE), list(chunk = 16L, perm = TRUE, wt = 2.5),
                    list(chunk = 41L, perm = TRUE, char = TRUE), list(chunk = 1L, perm = FALSE, wt = 0.375))) {
    fx <- nb_fixture(chunk = case$chunk, perm = case$perm, char_ids = isTRUE(case$char),
      wt = if (is.null(case$wt)) 1 else case$wt)
    for (method in c("or", "dr")) {
      fit <- nb_fit(method, fx)
      expect_false(is.null(fit$complete_support))
      rec <- fit$outcome_recipe; specs <- fit$support_specs; ctrl <- fit$control
      zc <- if (method == "or") rec$intercept_index else fit$or_frame_index[rec$intercept_index]
      ec <- nb_counted()
      native <- ec$or_baseline_complete_supports(fx$reader, specs, rec, ctrl, fit$complete_support, 1e9, zc,
        if (method == "or") api$or_baseline_abort else api$dr_baseline_abort)
      r_complete <- lapply(specs, function(s) if (method == "or")
        api$or_baseline_support(fx$reader, s, rec, ctrl, NULL, fit$complete_support) else
        api$dr_baseline_support(fx$reader, s, rec, fit$or_frame_index, ctrl, NULL, fit$complete_support))
      general <- lapply(specs, function(s) api$or_baseline_support(fx$reader, s, rec, ctrl, NULL, NULL))
      expect_identical(native, r_complete)
      expect_identical(native, general)
      n_k <- vapply(specs, function(s) fx$reader$state$cohort_end[s$cohort] - fx$reader$state$cohort_start[s$cohort] + 1,
        numeric(1))
      chunks <- as.integer(sum(ceiling(n_k / case$chunk)))
      # One native call and the two retained QR compressions per support chunk;
      # two slot checks and two endpoint reads per support row; no R pair rows.
      expect_identical(c(ec$n$calls[[1L]], ec$n$compress, ec$n$pair_rows), c(chunks, 2L * chunks, 0L))
      expect_identical(c(ec$n$slots, ec$n$reads), rep(2 * sum(n_k), 2L))
      expect_identical(fit$accounting$complete_support$n_support_chunks, as.double(chunks))
    }
  }
  # The intercept-column fact and the fit scratch reservation are enforced.
  fx <- nb_fixture(); fit <- nb_or(fx)
  expect_error(api$or_baseline_complete_supports(fx$reader, fit$support_specs, fit$outcome_recipe, fit$control,
    fit$complete_support, 2 * 7 + 63, 1L), class = "cdid_ipw_owner_bound_error")
})

test_that("native OR and DR chunks equal the general route for row sets, link orders, repeated pairs and scales", {
  pairs <- c(nb_pairs(), list(c(2, 4), c(1, 5), c(2, 4)))
  for (case in list(list(fx = nb_fixture(pairs = pairs), rows = list(1:7, c(41:35), 5L, c(3, 17, 9, 40))),
                    list(fx = nb_fixture(pairs = pairs, wt = 2.5), rows = list(1:7, c(40, 2, 23))),
                    list(fx = nb_fixture(nb_panel(extra = TRUE), chunk = 16L, perm = TRUE, char_ids = TRUE, wt = 0.375),
                      rows = list(1:16, 41L, c(16:1)))) ) {
    fx <- case$fx; L <- length(fx$views)
    for (method in c("or", "dr")) {
      fc <- nb_fit(method, fx); fg <- nb_fit(method, fx, nb_general)
      expect_false(is.null(fc$complete_support)); expect_null(fg$complete_support)
      for (rows in case$rows) for (links in list(seq_len(L), rev(seq_len(L)), c(12L, 3L, 7L), L)) for (scale in c("frame", "core")) {
        links <- links[links <= L]
        a <- nb_strip(nb_chunk(method, fc, rows, links, scale)); b <- nb_strip(nb_chunk(method, fg, rows, links, scale))
        nb_expect_match(a, b)
      }
      # Detailed requests keep the existing route and agree with the general one.
      a <- nb_chunk(method, fc, 1:7, seq_len(L), "frame", TRUE); b <- nb_chunk(method, fg, 1:7, seq_len(L), "frame", TRUE)
      nb_expect_match(a$details, b$details); nb_expect_match(a$phi, b$phi)
      cov <- function(f) if (method == "or") api$cp04_stream_covariance(f, seq_along(f$links)) else api$stream_dr_covariance(f)
      cc <- cov(fc); cg <- cov(fg)
      for (nm in c("estimates", "covariance", "mean_phi", "chunk_rows")) nb_expect_match(cc[[nm]], cg[[nm]], info = nm)
    }
  }
})

test_that("independent raw-row reference: the native OR IF columns are the intercept-only formula", {
  fx <- nb_fixture(chunk = 16L, perm = TRUE); st <- fx$reader$state; N <- st$N
  fit <- nb_or(fx)
  phi <- do.call(rbind, lapply(seq.int(1L, N, by = 16L), function(s)
    api$evaluate_or_baseline_if_chunk(fit, seq.int(s, min(N, s + 15L)), seq_along(fit$links))$phi))
  y_at <- function(id, t) fx$d$y[match(paste(id, t), paste(fx$d$id, fx$d$t))]
  tg <- st$g == 3; nv <- is.infinite(st$g)
  for (j in seq_along(fx$pairs)) {
    p <- fx$pairs[[j]]; dy <- y_at(st$id, p[2]) - y_at(st$id, p[1])
    alpha <- mean(dy[nv]); mu <- mean(dy[tg]) - alpha
    expect_lt(max(abs(phi[, j] - ifelse(tg, N / sum(tg) * (dy - alpha - mu), -N / sum(nv) * (dy - alpha)))), 1e-9)
  }
})

test_that("native batches count their calls, gathers and solves per row chunk", {
  fx <- nb_fixture(); N <- 41L; L <- length(fx$views)
  for (method in c("or", "dr")) {
    ec <- nb_counted()
    fit <- if (method == "or") ec$fit_or_baseline(fx$reader, fx$views, fx$ctrl, 1e-10, outcome_recipe = nb_icpt) else
      ec$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = nb_glm, rank_tol = 1e-10, outcome_recipe = nb_icpt)
    # Fit: pass 1 and pass 2 run once in native C; DR adds one first-model pair_rows call and one remaining-model summary batch; 0 QR at p = 0.
    expect_identical(ec$n$calls[["C_gmm_cp08_fit_pass1"]], 1L)
    expect_identical(ec$n$calls[["C_gmm_cp08_fit_pass2"]], 1L)
    expect_identical(ec$n$calls[["C_gmm_cp08_pair_rows"]], if (method == "or") 0L else 1L)
    expect_identical(ec$n$calls[["C_gmm_cp08_dr_treatment_batch"]], if (method == "or") 0L else 1L)
    expect_identical(ec$n$compress, 0L)
    expect_identical(ec$n$pair_rows, 0L)
    ec$n$calls[] <- 0L; ec$n$reads <- 0; ec$n$slots <- 0; ec$n$pair_rows <- 0L; ec$n$solve <- 0L
    replay <- if (method == "or") ec$replay_or_baseline_if else ec$replay_dr_baseline_if
    # (The OR replay hands each chunk to its callback as a promise; force it.)
    info <- replay(fit, 7L, L, function(b) force(b))
    chunks <- ceiling(N / 7)
    expect_identical(info$batches, as.integer(chunks))
    expect_identical(unname(ec$n$calls), if (method == "or") c(0L, 0L, 0L, 0L, as.integer(chunks), 0L, 0L, 0L) else
      c(0L, 0L, 0L, 0L, 0L, as.integer(chunks), as.integer(chunks), 0L))
    # Every row's five dates are gathered once per chunk (5 N reads), not twice
    # per pair support (2 x 10 N); slots likewise; one solve per DR chunk.
    expect_identical(c(ec$n$reads, ec$n$slots), c(5 * N, 5 * N))
    expect_identical(c(ec$n$pair_rows, ec$n$solve), c(0L, if (method == "or") 0L else as.integer(chunks)))
  }
})

test_that("native preflight, plan identity, mutants and malformed kernel arguments refuse", {
  fx <- nb_fixture(); L <- length(fx$views)
  for (method in c("or", "dr")) {
    fit <- nb_fit(method, fx); fg <- nb_fit(method, fx, nb_general)
    budget <- if (method == "or") api$or_baseline_if_budget else api$dr_baseline_if_budget
    planner <- if (method == "or") api$or_baseline_native_plan else api$dr_baseline_native_plan
    live <- if (method == "or") function(p, B) api$or_baseline_native_live(B, L, length(p$gids), length(p$dates)) else
      function(p, B) api$dr_baseline_native_live(B, L, length(p$tids), length(p$dates))
    # The plan and one chunk's native entries sit inside the unchanged scratch.
    for (B in c(1L, 7L)) {
      pre <- budget(fit, B, seq_len(L)); plan <- planner(fit, seq_len(L), B)
      expect_lte(plan$entries + live(plan, B), pre$scratch)
      pre$scratch <- plan$entries + live(plan, B) - 1
      expect_error(planner(fit, seq_len(L), B, pre), class = "cdid_ipw_owner_bound_error")
    }
    # One-row insufficient budget: refused before any kernel call or callback.
    pre1 <- budget(fit, 1, seq_len(L))
    edge <- fit$accounting$owned_entries + fit$accounting$borrowed_entries + pre1$scratch
    small <- nb_reseal(fit, nb_binding(method), function(x) { x$control$max_matrix_entries <- edge - 1; x })
    ec <- nb_counted(); calls <- 0L
    replay <- if (method == "or") ec$replay_or_baseline_if else ec$replay_dr_baseline_if
    expect_error(replay(small, 1L, L, function(b) calls <<- calls + 1L), class = "cdid_ipw_owner_bound_error")
    expect_identical(c(calls, sum(ec$n$calls)), c(0L, 0L))
    # Plan identity: a foreign batch, an edited unsealed plan and a larger chunk refuse.
    plan <- planner(fit, 1:3, 7L)
    expect_error(nb_chunk(method, fit, 1:7, 1:4, "frame", FALSE, plan), class = nb_err(method, "integrity"))
    edited <- plan; edited$B <- 41
    expect_error(nb_chunk(method, fit, 1:7, 1:3, "frame", FALSE, edited), class = nb_err(method, "integrity"))
    expect_error(nb_chunk(method, fit, 1:2, 1:3, "frame", FALSE, planner(fit, 1:3, 1L)), class = nb_err(method, "integrity"))
    good <- nb_strip(nb_chunk(method, fit, 1:7, 1:3, "frame", FALSE, plan))
    nb_expect_match(good, nb_strip(nb_chunk(method, fg, 1:7, 1:3, "frame")))
    # Deliberate wrong native term and wrong index (resealed mutants): the
    # general-route equality check detects both.
    kind <- if (method == "or") "or_baseline_native_plan" else "dr_baseline_native_plan"
    term <- nb_reseal(plan, kind, function(x) {
      if (method == "or") x$lnk_par[1L, 2L] <- x$lnk_par[1L, 2L] + 2^-20 else x$mdl_par[2L, 2L] <- x$mdl_par[2L, 2L] + 2^-20
      x })
    index <- nb_reseal(plan, kind, function(x) {
      if (method == "or") x$lnk[1L, 2:3] <- x$lnk[1L, 3:2] else x$lnk[2:3] <- x$lnk[3:2]
      x })
    for (mutant in list(term, index)) {
      bad <- nb_strip(nb_chunk(method, fit, 1:7, 1:3, "frame", FALSE, mutant))
      expect_false(identical(bad$phi, good$phi))
      expect_identical(bad$phi[, 1L], good$phi[, 1L])
    }
  }
  # Malformed kernel arguments abort with typed bound errors; nothing is read.
  fit <- nb_or(fx); plan <- api$or_baseline_native_plan(fit, seq_len(L), 7L)
  po <- get("panel_owner", fx$po$ref); post <- po$postings; st <- fx$reader$state
  k_or <- api$or_baseline_native_kernel("C_gmm_cp08_or_if_chunk")
  or_call <- function(...) {
    a <- list(po$observations$y, post$date_start, post$post_a_index, post$post_obs_index, 1:7, st$cohort_index,
      st$a_index, st$w, plan$dates, plan$sup, plan$sup_par, plan$blk, plan$blk_par, plan$lnk, plan$lnk_par,
      plan$dimnames, 1e7)
    edits <- list(...); for (k in names(edits)) a[[as.integer(k)]] <- edits[[k]]
    do.call(.Call, c(list(k_or), a))
  }
  expect_identical(or_call()$status, c(0L, 0L, 0L, 0L))
  bad <- list(`5` = as.double(1:7), `5` = c(1L, 99L), `9` = c(plan$dates, NA_integer_), `10` = plan$sup[1:2, ],
    `12` = c(plan$blk[-1L], 99L), `14` = plan$lnk + 1000L, `16` = list(NULL, "x"), `17` = 10,
    `8` = st$w[-1L], `11` = matrix(plan$sup_par, 1L))
  for (k in seq_along(bad)) expect_error(do.call(or_call, bad[k]), class = "cdid_gmm_bound_error", info = names(bad)[k])
  k_sup <- api$or_baseline_native_kernel("C_gmm_cp08_support_chunk")
  sup_call <- function(ctl = c(1L, 7L, 1L), dates = c(1L, 2L), ld = TRUE, z = fx$fo$Z_a) .Call(k_sup,
    po$observations$y, post$date_start, post$post_a_index, post$post_obs_index, st$cohort_order, st$a_index, st$w,
    z, ctl, dates, ld, 1e7)
  expect_identical(sup_call()$status, c(0L, 0L, 0L, 0L))
  for (args in list(list(ctl = c(1L, 99L, 1L)), list(ctl = c(1L, 7L, 2L)), list(dates = c(1L, 9L)),
                    list(ld = NA), list(z = as.vector(fx$fo$Z_a)), list(ctl = c(7L, 1L, 1L))))
    expect_error(do.call(sup_call, args), class = "cdid_gmm_bound_error")
  k_pr <- api$or_baseline_native_kernel("C_gmm_cp08_pair_rows")
  pr_call <- function(rows = 1:7, dates = c(1L, 2L), chunk = 7L) .Call(k_pr, po$observations$y, post$date_start,
    post$post_a_index, post$post_obs_index, rows, st$a_index, dates, chunk, 1e7)
  expect_identical(pr_call()$status, c(0L, 0L, 0L, 0L))
  for (args in list(list(rows = c(2L, 1L)), list(rows = c(1L, 99L)), list(dates = c(0L, 2L)), list(chunk = 0L),
                    list(rows = as.double(1:3))))
    expect_error(do.call(pr_call, args), class = "cdid_gmm_bound_error")
  expect_error(api$or_baseline_native_kernel("C_gmm_cp08_absent"), class = "cdid_gmm_native_error")
})

test_that("native and general routes raise the same first failure", {
  rebind <- function(env, nm, value) {
    old <- get(nm, envir = env)
    unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
    function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
  }
  fx <- nb_fixture(); po <- fx$po; L <- length(fx$views)
  y <- po$observations$y
  cell <- function(id, t) po$postings$post_obs_index[po$postings$date_start[t] - 1L + match(id, po$a_id)]
  bad_y <- list(endpoint = { z <- y; z[cell(4, 5)] <- NaN; z },
                difference = { z <- y; z[cell(4, 2)] <- 1e308; z[cell(4, 3)] <- -1e308; z })
  for (method in c("or", "dr")) {
    fc <- nb_fit(method, fx); fg <- nb_fit(method, fx, nb_general)
    for (nm in names(bad_y)) {
      bp <- po; bp$observations$y <- bad_y[[nm]]
      restore <- rebind(po$ref, "panel_owner", bp)
      for (what in c("replay", "fit")) {
        run <- function(f, ns) if (what == "fit") nb_fit(method, fx, ns) else
          (if (method == "or") api$replay_or_baseline_if else api$replay_dr_baseline_if)(f, 7L, L, function(b) force(b))
        ec <- tryCatch(run(fc, api), error = identity); eg <- tryCatch(run(fg, nb_general), error = identity)
        expect_s3_class(ec, "cdid_error")
        expect_identical(class(ec), class(eg), info = paste(method, nm, what))
        expect_identical(conditionMessage(ec), conditionMessage(eg), info = paste(method, nm, what))
      }
      restore()
    }
  }
})

test_that("common-design fit_pass1, fit_pass2 and IF kernels count calls, QR compressions and refuse malformed inputs", {
  N <- 41L; chunk <- 7L
  Xx <- cbind(`(Intercept)` = 1, x1 = sin(seq_len(N) * 0.41), x2 = cos(seq_len(N) * 0.59))
  rec_x <- api$or_pooled_outcome_recipe(c("(Intercept)", "x1", "x2"), TRUE, c("x1", "x2"))
  fx <- nb_fixture(chunk = chunk, X = Xx, V = Xx)
  L <- length(fx$views); st <- fx$reader$state
  ctrl_coh <- match(Inf, st$cohorts)
  n_c <- st$cohort_end[ctrl_coh] - st$cohort_start[ctrl_coh] + 1L
  n_c_chunks <- as.integer(ceiling(n_c / chunk))

  for (method in c("or", "dr")) {
    ec <- nb_counted()
    fit <- if (method == "or") ec$fit_or_baseline(fx$reader, fx$views, fx$ctrl, 1e-10, outcome_recipe = rec_x) else
      ec$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = nb_glm, rank_tol = 1e-10,
        outcome_recipe = rec_x, treatment_columns = c("(Intercept)", "x1", "x2"))
    fg <- nb_fit(method, fx, nb_general, recipe = rec_x, treatment_columns = c("(Intercept)", "x1", "x2"))

    expect_false(is.null(fit$complete_support))
    expect_identical(fit$shared$p, 2L)
    expect_identical(fit$complete_support$n_z, 3L)
    expect_identical(ec$n$calls[["C_gmm_cp08_fit_pass1"]], 1L)
    expect_identical(ec$n$calls[["C_gmm_cp08_fit_pass2"]], 1L)
    expect_identical(ec$n$calls[["C_gmm_cp08_pair_rows"]], if (method == "or") 0L else 1L)
    expect_identical(ec$n$calls[["C_gmm_cp08_dr_treatment_batch"]], if (method == "or") 0L else 1L)
    # Common-design QR compressions: n_c_chunks chunk QRs + (n_c_chunks - 1) chunk merges + 1 partial-block QR = 2 * n_c_chunks,
    # instead of 2 * n_chunks across all S = 20 supports + G = 10 block partial QRs + (G - 1) merges!
    expect_identical(ec$n$compress, 2L * n_c_chunks)
    expect_identical(fit$accounting$complete_support$n_common_qr_compressions, 2L * n_c_chunks)
    expect_identical(ec$n$pair_rows, 0L)

    # IF replay counts and gathers for p_or = 2, p_trt = 3:
    ec$n$calls[] <- 0L; ec$n$reads <- 0; ec$n$slots <- 0; ec$n$pair_rows <- 0L; ec$n$solve <- 0L
    replay <- if (method == "or") ec$replay_or_baseline_if else ec$replay_dr_baseline_if
    batches_c <- list()
    info <- replay(fit, chunk, L, function(b) { force(b); batches_c[[length(batches_c) + 1L]] <<- nb_strip(b) })
    batches_g <- list()
    invisible((if (method == "or") api$replay_or_baseline_if else api$replay_dr_baseline_if)(fg, chunk, L,
      function(b) { force(b); batches_g[[length(batches_g) + 1L]] <<- nb_strip(b) }))
    chunks_n <- as.integer(ceiling(N / chunk))
    expect_identical(info$batches, chunks_n)
    expect_identical(unname(ec$n$calls), if (method == "or") c(0L, 0L, 0L, 0L, chunks_n, 0L, 0L, 0L) else
      c(0L, 0L, 0L, 0L, 0L, chunks_n, chunks_n, 0L))
    expect_identical(c(ec$n$reads, ec$n$slots), c(5 * N, 5 * N))
    expect_identical(c(ec$n$pair_rows, ec$n$solve), c(0L, if (method == "or") 0L else chunks_n))
    nb_expect_match(batches_c, batches_g, tol = 1e-12)
  }

  # Malformed arguments to C_gmm_cp08_fit_pass1 and C_gmm_cp08_fit_pass2 abort with cdid_gmm_bound_error:
  po <- get("panel_owner", fx$po$ref); post <- po$postings; frame <- get("frame_owner", fx$reader$frame_ref)
  fit <- nb_or(fx, recipe = rec_x)
  sup_mat <- vapply(fit$support_specs, function(s) as.integer(c(st$cohort_start[s$cohort], st$cohort_end[s$cohort],
    s$pair$s_calendar_index, s$pair$t_calendar_index)), integer(4))
  L_vec <- vapply(fit$regression_specs, `[[`, numeric(1), "transitions")
  ctrl_sids <- vapply(fit$regression_specs, function(pl) pl$support_ids[1L], integer(1))
  ystar_w <- numeric(length(fit$support_specs)); ystar_w[ctrl_sids] <- L_vec / sum(L_vec^2)
  c_first <- as.integer(st$cohort_start[ctrl_coh]); c_last <- as.integer(st$cohort_end[ctrl_coh])
  k_p1 <- api$or_baseline_native_kernel("C_gmm_cp08_fit_pass1")
  p1_call <- function(...) {
    a <- list(po$observations$y, post$date_start, post$post_a_index, post$post_obs_index, st$cohort_order,
      st$a_index, st$w, frame$Z_a, sup_mat, c(2L, 3L), ystar_w, c(1L, 7L, c_first, c_last), TRUE, 1e7, NULL)
    edits <- list(...); for (k in names(edits)) a[[as.integer(k)]] <- edits[[k]]
    do.call(.Call, c(list(k_p1), a))
  }
  p1_ok <- p1_call()
  expect_identical(p1_ok$status, c(0L, 0L, 0L, 0L))
  bad_p1 <- list(list(`9` = sup_mat[1:3, ]), list(`10` = c(2L, 99L)), list(`11` = ystar_w[-1L]),
    list(`12` = c(0L, 7L, c_first, c_last)), list(`12` = c(1L, 7L, c_last, c_first)), list(`13` = NA), list(`14` = 5))
  for (ed in bad_p1) expect_error(do.call(p1_call, ed), class = "cdid_gmm_bound_error")

  k_p2 <- api$or_baseline_native_kernel("C_gmm_cp08_fit_pass2")
  blk_mat <- vapply(fit$regression_specs, function(pl) {
    pr <- fit$support_specs[[pl$support_ids[1L]]]$pair
    as.integer(c(pr$s_calendar_index, pr$t_calendar_index))
  }, integer(2))
  blk_par <- rbind(L_vec, rep(1, length(L_vec)), vapply(fit$regressions, `[[`, numeric(1), "alpha"), rep(0, length(L_vec)))
  p2_call <- function(...) {
    a <- list(po$observations$y, post$date_start, post$post_a_index, post$post_obs_index, st$cohort_order,
      st$a_index, st$w, frame$Z_a, blk_mat, blk_par, c(2L, 3L), fit$shared$beta, p1_ok$ystar, c(c_first, c_last, 7L), 1e7)
    edits <- list(...); for (k in names(edits)) a[[as.integer(k)]] <- edits[[k]]
    do.call(.Call, c(list(k_p2), a))
  }
  expect_identical(p2_call()$status, c(0L, 0L, 0L, 0L))
  # The second pass must validate positions before dereferencing either posting.
  bad_start <- post$date_start; bad_start[blk_mat[1L, 1L]] <- NA_integer_
  expect_identical(p2_call(`2` = bad_start)$status[1L], 1L)
  bad_start[blk_mat[1L, 1L]] <- 0L
  expect_identical(p2_call(`2` = bad_start)$status[1L], 1L)
  # A finite weighted column norm can have an overflowing square. Only the
  # actual original QR validator may refuse this; the squared norm is no test.
  y_large <- rep(1e155, length(po$observations$y))
  # Alternate stored levels create finite +/-1e155 raw changes.
  y_large[seq.int(1L, length(y_large), by = 2L)] <- 0
  validation_calls <- 0L
  validation <- function(sid) {
    validation_calls <<- validation_calls + 1L
    ff <- api$or_wls_compress(cbind(rep(1, 7L), rep(1e155, 7L)))
    stopifnot(all(is.finite(ff)))
  }
  expect_identical(p1_call(`1` = y_large, `15` = validation)$status, c(0L, 0L, 0L, 0L))
  expect_gt(validation_calls, 0L)
  expect_error(p1_call(`15` = 1), class = "cdid_gmm_bound_error")
  bad_p2 <- list(list(`9` = blk_mat[1L, , drop = FALSE]), list(`10` = blk_par[1:3, ]), list(`11` = c(2L, 99L)),
    list(`12` = fit$shared$beta[1L]), list(`13` = p1_ok$ystar[-1L]), list(`14` = c(c_first, c_last, 0L)), list(`15` = 5))
  for (ed in bad_p2) expect_error(do.call(p2_call, ed), class = "cdid_gmm_bound_error")
})

test_that("zero-variable provenance tables keep their columns and reach the intended routes", {
  api_classify <- api$chained_reference_classify
  empty <- api_classify(character(0), data.frame(id = 1), NULL, FALSE, 1, 1, "t", 1, 1L, "x")
  expect_identical(empty, data.frame(variable = character(0), source = character(0), kind = character(0),
    example = character(0), stringsAsFactors = FALSE))
  d <- data.frame(id = rep(1:2, each = 2), t = rep(1:2, 2), x = c(1, 1, 2, 3))
  one <- api_classify("x", d, NULL, FALSE, d$id, d$t, "t", 1:2, d$id, "x")
  expect_identical(names(one), c("variable", "source", "kind", "example"))
  expect_identical(c(one$variable, one$source, one$kind), c("x", "data", "dated"))
  # Intercept-only reg with within-unit changing weights (fix_weights =
  # base_period) prepares its reference-recipe estimation inputs; ~0 keeps the
  # original typed refusal.
  set.seed(31L)
  N <- 41L; T <- 5L
  dt <- data.frame(id = rep(seq_len(N), each = T), period = rep(seq_len(T), N))
  dt$G <- ifelse(dt$id %% 2 == 0, 3, 0)
  dt$y <- stats::rnorm(nrow(dt)) + 0.1 * dt$period
  dt$w <- 1 + (dt$id %% 3) / 4 + dt$period / 10
  resolved <- api$chained_public_resolve_control(NULL)
  prep <- api$chained_request_prepare_inputs(data = dt, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = NULL, weightsname = "w", outcome_regression = TRUE, or_resource = resolved$profile,
    fix_weights = "base_period")
  expect_identical(prep$weight_reduction, "reference_recipes")
  expect_identical(prep$reference_route, "baseline")
  expect_false(is.null(prep$reference_source))
  expect_error(api$chained_request_prepare_inputs(data = transform(dt, w = NULL), yname = "y", tname = "period",
    idname = "id", gname = "G", xformla = ~0, outcome_regression = TRUE, or_resource = resolved$profile),
    class = "cdid_chained_request_error")
  # Genuine public requests (the installed package with its pinned did).
  skip_if_not_installed("did")
  mp <- cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = dt, est_method = "reg",
    weightsname = "w", fix_weights = "base_period", gmm_weighting = "identity", control_group = "nevertreated",
    bstrap = FALSE, cband = FALSE, base_period = "varying", mode = "chained", pt_assumption = "post")
  expect_s3_class(mp, "cdid_MP")
  expect_true(all(is.finite(mp$att)))
  expect_false(is.null(mp$chained$reference_inputs))
  expect_error(cdid::att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = transform(dt, w = NULL),
    est_method = "reg", xformla = ~0, bstrap = FALSE, cband = FALSE, mode = "chained"), class = "cdid_chained_request_error")
})

test_that("the demo uses the setup library when present, with explicit overrides first", {
  script <- normalizePath(file.path("..", "..", "scripts", "positron_demo.R"), mustWork = FALSE)
  skip_if_not(file.exists(script), "source-tree demo script not present in this test location")
  exprs <- parse(script, keep.source = FALSE)
  at <- which(vapply(exprs, function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
    identical(e[[2L]], as.name("package_library")), logical(1)))
  expect_length(at, 1L)
  expect_false(any(grepl("final_package_gate_v1|/Users/", deparse(exprs[[at]]))))
  choose <- function(setup_library, env) {
    root <- tempfile("cdid-demo-root-"); dir.create(root)
    old <- setwd(root); on.exit({ setwd(old); unlink(root, recursive = TRUE) }, add = TRUE)
    if (setup_library) dir.create("validation/.lib", recursive = TRUE)
    old_env <- Sys.getenv("CDID_DEMO_PACKAGE_LIBRARY", unset = NA)
    on.exit(if (is.na(old_env)) Sys.unsetenv("CDID_DEMO_PACKAGE_LIBRARY") else
      Sys.setenv(CDID_DEMO_PACKAGE_LIBRARY = old_env), add = TRUE)
    if (is.null(env)) Sys.unsetenv("CDID_DEMO_PACKAGE_LIBRARY") else Sys.setenv(CDID_DEMO_PACKAGE_LIBRARY = env)
    e <- new.env(parent = baseenv()); eval(exprs[[at]], e); e$package_library
  }
  expect_identical(choose(TRUE, NULL), "validation/.lib")
  expect_null(choose(FALSE, NULL))
  expect_identical(choose(TRUE, "/explicit/lib"), "/explicit/lib")
  expect_identical(choose(FALSE, "/explicit/lib"), "/explicit/lib")
})
