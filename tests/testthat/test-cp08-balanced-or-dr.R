# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 balanced-panel performance, stage 3: the complete-support baseline OR/DR
# route (complete-date positions instead of endpoint searches, one shared DR
# treatment logit, one treatment replay per IF chunk, eligible-route charges),
# with the stage-4 native support preparation and non-detailed IF batches.
# Every expected value comes from the unchanged general route of the same
# build (forced by test-local instrumentation), from the accepted dense DR
# table helper, or from formulas written here over raw panel rows; the
# complete-support route never supplies its own truth.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

od_ctrl <- function(chunk = 7L, work = 1e10, matrix = 2e7) api$ipw_owner_control(max_chunk_rows = chunk,
  max_retained_entries = 1e7, max_matrix_entries = matrix, max_work_ops = work, max_output_entries = 1e7)
od_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
od_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))

# Balanced long panel: g = 3 treated, Inf never treated (optionally a second
# treated cohort 4), 41 units (partial last chunks for 7 and 16 rows), 5 dates.
od_panel <- function(N = 41L, T = 5L, seed = 11L, second = FALSE) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:2] <- c(3, Inf)
  if (second) g[3:8] <- 4
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$y <- 0.2 * d$id + 0.4 * d$t + ifelse(d$t >= d$g, 1.25, 0) + stats::rnorm(nrow(d))
  d
}
od_pairs <- function(T = 5L) {
  pr <- list()
  for (s in seq_len(T - 1L)) for (t in seq.int(s + 1L, T)) pr[[length(pr) + 1L]] <- c(s, t)
  pr
}

# Real sealed panel/frame owners, the full-roster reader and link views.
od_fixture <- function(d = od_panel(), chunk = 7L, perm = FALSE, w = NULL, X = NULL, V = NULL, pairs = od_pairs(),
                       K = Inf, target = 3, frame_only = FALSE, ctrl = od_ctrl(chunk)) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  u <- d[!duplicated(d$id), c("id", "g")]
  u <- u[order(u$id), ]
  if (frame_only) u <- rbind(u, data.frame(id = 999, g = Inf))
  if (perm) { set.seed(5L); u <- u[sample(nrow(u)), ] }
  N <- nrow(u)
  Z <- if (is.null(X)) matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)")) else X[match(u$id, sort(u$id)), , drop = FALSE]
  V_mat <- if (is.null(V)) matrix(1, N, 1L, dimnames = list(NULL, "reader_intercept")) else V[match(u$id, sort(u$id)), , drop = FALSE]
  wt <- if (is.null(w)) rep(1, N) else w[match(u$id, sort(u$id))]
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = wt, control = ctrl, frame_g = u$g,
    allow_frame_only = frame_only, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, u$id, V_mat,
    u$g, wt, ctrl, allow_unknown_off_A = TRUE)
  views <- lapply(seq_along(pairs), function(k) {
    p <- pairs[[k]]
    api$construct_ipw_link_view(fo, p[1], p[2], g = target[(k - 1L) %% length(target) + 1L], K = K, control = ctrl,
      link_id = paste0("L", p[1], "_", p[2],
        if (sum(vapply(pairs, function(q) identical(q, p), logical(1))) > 1L) paste0(":", k) else ""), functional_key = paste0("cmpcat:od:", p[1], "-", p[2]))
  })
  list(d = d, po = po, fo = fo, reader = reader, views = views, ctrl = ctrl, chunk = chunk, u = u, pairs = pairs)
}

# The unchanged general route of the same build: the fit functions re-bound to
# an environment whose route predicate declines.
od_general <- local({
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$or_baseline_complete_facts <- function(...) NULL
  e
})
od_or <- function(fx, ns = api, ctrl = fx$ctrl, recipe = od_icpt) ns$fit_or_baseline(fx$reader, fx$views, ctrl, 1e-10,
  outcome_recipe = recipe)
od_dr <- function(fx, ns = api, ctrl = fx$ctrl, recipe = od_icpt, glm = od_glm, treatment_columns = NULL)
  ns$fit_dr_baseline(fx$reader, fx$views, ctrl, glm_control = glm, rank_tol = 1e-10,
    outcome_recipe = recipe, treatment_columns = treatment_columns)

# Test-local counters around the fit and IF chunk call graph. Stage 4 + common-design:
# the native complete-support adapters are re-bound too, and every .Call of a CP08
# kernel is counted by name.
od_kernels <- c("C_gmm_cp08_support_chunk", "C_gmm_cp08_pair_rows", "C_gmm_cp08_fit_pass1",
  "C_gmm_cp08_fit_pass2", "C_gmm_cp08_or_if_chunk", "C_gmm_cp08_dr_if_prepare", "C_gmm_cp08_dr_if_assemble", "C_gmm_cp08_dr_treatment_batch")
od_counted <- function(decline = FALSE) {
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline", "or_baseline_support", "or_baseline_visit", "dr_baseline_support",
               "dr_support_visit", "dr_baseline_treatment_fit", "dr_baseline_treatment_batch", "dr_baseline_shared_logit", "replay_or_baseline_if",
               "evaluate_or_baseline_if_chunk", "replay_dr_baseline_if", "evaluate_dr_baseline_if_chunk",
               "dr_baseline_complete_replay", "or_baseline_complete_supports", "or_baseline_complete_state",
               "or_pooled_fit_common", "or_baseline_native_if_chunk", "dr_baseline_native_if_chunk")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  n <- new.env(parent = emptyenv())
  n$logit <- n$search <- n$complete <- n$replay <- n$core <- 0L
  n$native <- setNames(integer(length(od_kernels)), od_kernels)
  addr <- lapply(od_kernels, api$or_baseline_native_kernel)
  e$.Call <- function(.NAME, ...) {
    k <- which(vapply(addr, identical, logical(1), .NAME))
    if (length(k)) n$native[k] <- n$native[k] + 1L
    base::.Call(.NAME, ...)
  }
  e$cp04_logit <- function(...) { n$logit <- n$logit + 1L; api$cp04_logit(...) }
  e$or_pooled_if_core <- function(...) { n$core <- n$core + 1L; api$or_pooled_if_core(...) }
  e$cp04_pair_rows <- function(...) { n$search <- n$search + 1L; api$cp04_pair_rows(...) }
  e$or_baseline_complete_pair_rows <- function(...) { n$complete <- n$complete + 1L; api$or_baseline_complete_pair_rows(...) }
  e$ipw_logit_values <- function(...) { n$replay <- n$replay + 1L; api$ipw_logit_values(...) }
  if (decline) e$or_baseline_complete_facts <- function(...) NULL
  e$n <- n
  e
}

od_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, od_strip); attributes(x) <- a }
  x
}
od_norm_shared <- function(s) {
  if (!is.list(s)) return(s)
  norm_blocks <- function(blks) {
    if (!is.list(blks)) return(blks)
    lapply(blks, function(bk) {
      if (is.list(bk) && !is.null(bk$r11)) {
        sgn <- if (bk$r11 < 0) -1 else 1
        bk$r11 <- abs(bk$r11)
        if (!is.null(bk$r1d)) bk$r1d <- sgn * bk$r1d
        if (!is.null(bk$r1y)) bk$r1y <- sgn * bk$r1y
      }
      bk
    })
  }
  drop_acc <- c("max_chunk_merge_rows", "max_residual_merge_rows", "n_compressions")
  if (is.list(s$fit)) {
    s$fit$blocks <- norm_blocks(s$fit$blocks)
    if (is.list(s$fit$accounting)) s$fit$accounting <- s$fit$accounting[setdiff(names(s$fit$accounting), drop_acc)]
    if (is.matrix(s$fit$R) && ncol(s$fit$R) > 0L) {
      piv <- s$fit$pivot
      R_unpiv <- s$fit$R[, order(piv), drop = FALSE]
      s$fit$R <- crossprod(R_unpiv)
      s$fit$pivot <- seq_len(ncol(s$fit$R))
      if (is.list(s$fit$slope_fit)) {
        s$fit$slope_fit$R <- s$fit$R
        s$fit$slope_fit$pivot <- s$fit$pivot
      }
    }
  }
  if (is.list(s$summary)) {
    s$summary$blocks <- norm_blocks(s$summary$blocks)
    if (is.list(s$summary$accounting)) s$summary$accounting <- s$summary$accounting[setdiff(names(s$summary$accounting), drop_acc)]
    if (is.list(s$summary$residual) && is.matrix(s$summary$residual$factor)) {
      s$summary$residual$factor <- crossprod(s$summary$residual$factor)
    }
  }
  s
}
# Route-only differences: the complete-support record, route merge counters and resource bounds.
# Route work is a resource contract, not an equal-estimator output. Recompute
# its reservation from input rows, catalog dimensions and the documented compact
# payload convention; do not call the implementation's charge/payload helpers.
od_expect_work <- function(b, treatment_route = "batch") {
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st <- reader$state; N <- length(st$id); nA <- sum(st$a_index > 0L)
  S <- length(b$support_specs); G <- length(b$regression_specs); L <- length(b$links)
  p <- b$outcome_recipe$p; q <- p + 2; B <- min(N, b$control$max_chunk_rows)
  dr <- inherits(b, "cdid_dr_baseline")
  # Reconstruct only the pre-fit link descriptors; fitted moments are not
  # metadata charged by the fit preflight. Preserve their original field order.
  keys <- c("regression", "pair_id", "cohort", "link_id", "source_functional_key", "functional_key",
    if (!is.null(b$reference_ref)) "reference", "support_id",
    if (dr) c("model", "control_support_ids"))
  links <- lapply(b$links, function(l) l[keys])
  meta <- list(b$support_specs, b$regression_specs, links)
  if (dr) meta <- c(meta, list(b$treatment_specs, b$treatment_recipe))
  entries <- function(x, attrs = TRUE) {
    if (is.null(x) || is.environment(x)) return(0)
    z <- length(x) + 8.0
    if (is.list(x)) z <- z + sum(vapply(x, entries, numeric(1)))
    if (is.character(x)) z <- z + sum(ceiling(nchar(x[!is.na(x)], type = "bytes") / 2))
    if (attrs && length(attributes(x))) z <- z + entries(attributes(x), FALSE)
    z
  }
  metadata <- entries(meta); header <- 4096 + 1024 * L * (length(st$cohorts) + 1)
  common <- !is.null(b$complete_support)
  if (!common) {
    own <- 512 * S * N * q^2 * (1 + ceiling(log2(nA + 1))) +
      256 * S * N * q * (p + G + 2) + 256 * G * (S + 1) * q^3 +
      256 * (G + S) * q^2 * (G + 2) + 128 * (metadata + header + L * q)
  } else {
    nk <- vapply(b$support_specs, function(s) sum(st$g == st$cohorts[s$cohort]), numeric(1))
    nc_vec <- vapply(b$complete_support$control, function(ck) sum(st$g == st$cohorts[ck]), numeric(1))
    nc <- sum(nc_vec); c_chunks <- sum(ceiling(nc_vec / B)); rows <- sum(nk)
    multi_ctrl <- length(b$complete_support$control) > 1L
    dates <- unique(unlist(lapply(b$support_specs, function(s) c(s$pair$s_calendar_index, s$pair$t_calendar_index))))
    used <- b$outcome_recipe$columns
    if (dr) used <- union(used, b$treatment_recipe$columns)
    facts <- length(dates) * nA + (1 + length(used)) * N + 8 * L
    qr <- (if (p > 0) 256 * nc * q^2 + 256 * c_chunks * q^3 + 512 * q^3 + 256 * rows * p^2 else 0) +
      if (multi_ctrl) 256 * (G * nc) * q^2 + 256 * (G * c_chunks + G) * q^3 else 0
    own <- 256 * (rows + G * nc) * q + qr + 256 * (G + S) * q +
      128 * (metadata + header + L * q) + 64 * facts
    # These bounded comparison fixtures do not trigger raw-factor overflow
    # validation. Leave the separate large-value validation tests intact.
    expect_equal(b$accounting$complete_support$original_factor_validation_work, 0, tolerance = 0)
  }
  if (dr) {
    Tn <- length(b$treatment_specs); pt <- b$treatment_recipe$p
    obs <- vapply(b$treatment_specs, function(s) sum(st$g %in% st$cohorts[s$cohorts]), numeric(1))
    expect_equal(unname(vapply(b$treatment_specs, `[[`, integer(1), "n_observers")), obs, tolerance = 0)
    m <- max(vapply(b$treatment_specs, function(s) length(s$cohorts), integer(1)))
    models <- Tn * (pt^2 + 12 * (pt + p) + 8 * m + pt * m + 4096) + L * (8 * p + 8 * pt + 1024) + 8192
    if (!common) {
      trt <- 128 * (b$glm_control$maxit + 2) * sum(obs) * (pt + 2)^2 + 128 * Tn * N * m + 128 * models
    } else {
      trt_keys <- vapply(b$treatment_specs, function(s) paste(s$cohorts, collapse = ":"), character(1))
      trt_groups <- lapply(unique(trt_keys), function(k) which(trt_keys == k))
      n1 <- obs[1L]; rem <- max(0L, Tn - length(trt_groups))
      if (treatment_route == "batch" || rem == 0L) {
        trt <- sum(vapply(trt_groups, function(g_idx) {
          ng <- obs[g_idx[1L]]; rem_g <- max(0L, length(g_idx) - 1L)
          128 * (b$glm_control$maxit + 2) * ng * (pt + 2)^2 +
            128 * ng * ((pt + 2)^2 + p + m) + 64 * rem_g * ng * (p + pt + m + 8)
        }, numeric(1))) + 128 * models
      } else {
        # A deliberately declined summary batch reserves the old per-model
        # loop AND its attempted validation, before executing either.
        summary <- 12 + 3 * p + pt * m + 8 * (m + pt) + 16
        trt <- 128 * (b$glm_control$maxit + 2) * n1 * (pt + 2)^2 + 128 * models +
          128 * sum(obs) * ((pt + 2)^2 + p + m) + 128 * Tn * N * m +
          64 * rem * n1 * (p + pt + m + 8) + 128 * rem * summary
      }
    }
    own <- own + trt
  }
  expect_equal(b$accounting$own_work_ops, own, tolerance = 0)
  expect_equal(b$accounting$work_ops_bound, b$accounting$external_work_ops + own, tolerance = 0)
  expect_identical(b$accounting$work_ops_convention, "cumulative")
}

od_backend_view <- function(b, treatment_route = "batch") {
  od_expect_work(b, treatment_route)
  x <- unclass(b)[setdiff(names(b), c("ref", "reader_ref", "complete_support", "accounting"))]
  if (!is.null(x$shared)) x$shared <- od_norm_shared(x$shared)
  list(fields = od_strip(x), accounting = b$accounting[setdiff(names(b$accounting),
    c("complete_support", "work_ops_bound", "own_work_ops", "fit_scratch_bound", "retained_bound", "owned_entries"))])
}
od_expect_route_match <- function(a, b, tol = 1e-13, info = NULL, mixed_scalar_path = NULL) {
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
    # D75 applies only to the explicitly bound scalar in the three new cases below.
    same <- if (identical(path, mixed_scalar_path)) {
      typeof(x) == "double" && length(x) == 1L && is.finite(x) && is.finite(y) &&
        abs(x - y) <= 1e-12 * max(1, abs(x), abs(y))
    } else if (typeof(x) == "double") length(waldo::compare(x, y, tolerance = tol)) == 0L else identical(x, y)
    if (same) character() else paste(path, "value")
  }
  bad <- mismatch(a, b, "root")
  expect_true(length(bad) == 0L, info = paste(c(info, head(bad, 12L)), collapse = "\n"))
}
od_batches <- function(fit, replay, chunk, LB, scale = "frame", details = FALSE) {
  out <- list()
  info <- replay(fit, chunk, LB, function(b) {
    b$diagnostics <- b$diagnostics[setdiff(names(b$diagnostics),
      c("output_entries_bound", "concurrent_scratch_bound", "work_ops_bound"))]
    out[[length(out) + 1L]] <<- b
  }, scale, details)
  list(batches = out, info = info)
}

od_fx <- od_fixture()
od_fx_perm <- od_fixture(chunk = 16L, perm = TRUE)

test_that("eligible OR and DR fits, IF replays and covariances equal the general route", {
  for (case in list(list(fx = od_fx, chunk = 7L, LB = 4L, details = TRUE, scale = "frame"),
                    list(fx = od_fx_perm, chunk = 16L, LB = 10L, details = FALSE, scale = "core"))) {
    fx <- case$fx
    for (method in c("or", "dr")) {
      fit_c <- if (method == "or") od_or(fx) else od_dr(fx)
      fit_g <- if (method == "or") od_or(fx, od_general) else od_dr(fx, od_general)
      expect_false(is.null(fit_c$complete_support))
      expect_null(fit_g$complete_support)
      expect_null(fit_g$accounting$complete_support)
      od_expect_route_match(od_backend_view(fit_c), od_backend_view(fit_g))
      replay <- if (method == "or") api$replay_or_baseline_if else api$replay_dr_baseline_if
      bc <- od_batches(fit_c, replay, case$chunk, case$LB, case$scale, case$details)
      bg <- od_batches(fit_g, replay, case$chunk, case$LB, case$scale, case$details)
      if (method == "dr") {
        # At p_or = 0, both routes omit the unused OR core unless details = TRUE.
        n_reg <- function(b) vapply(b$batches, function(x) x$diagnostics$n_regression_score_evaluations, integer(1))
        if (!case$details) {
          expect_true(all(n_reg(bc) == 0L))
          expect_true(all(n_reg(bg) == 0L))
        } else {
          expect_true(all(n_reg(bc) > 0L))
          expect_identical(n_reg(bc), n_reg(bg))
        }
      }
      od_expect_route_match(bc$batches, bg$batches)
      expect_identical(bc$info[setdiff(names(bc$info), "work_ops_bound")], bg$info[setdiff(names(bg$info), "work_ops_bound")])
      expect_identical(bc$info$batches, as.integer(ceiling(41 / case$chunk) * ceiling(length(fx$views) / case$LB)))
      cov <- function(f) if (method == "or") api$cp04_stream_covariance(f, seq_along(f$links)) else api$stream_dr_covariance(f)
      cc <- cov(fit_c); cg <- cov(fit_g)
      for (nm in c("estimates", "covariance", "mean_phi", "chunk_rows"))
        od_expect_route_match(cc[[nm]], cg[[nm]], info = nm)
    }
  }
})

test_that("DR fits ONE logit, keeps every treatment model and replays the treatment score once per chunk", {
  fx <- od_fx; T_count <- length(fx$views)
  ec <- od_counted(); eg <- od_counted(decline = TRUE)
  fit_c <- od_dr(fx, ec); fit_g <- od_dr(fx, eg)
  expect_identical(ec$n$logit, 1L)
  expect_identical(eg$n$logit, T_count)
  # Complete positions replace every endpoint search; the general route searches.
  expect_identical(ec$n$search, 0L)
  expect_gt(eg$n$search, 0L)
  n_k <- vapply(fit_c$support_specs, function(s) {
    st <- fx$reader$state; st$cohort_end[s$cohort] - st$cohort_start[s$cohort] + 1 }, numeric(1))
  # Common-design fit passes 1 and 2 prepare the supports/moments in native C
  # and one first-model pair-row call plus one remaining-model batch read every original pair; no R pair-row query remains.
  expect_identical(ec$n$native[["C_gmm_cp08_fit_pass1"]], 1L)
  expect_identical(ec$n$native[["C_gmm_cp08_fit_pass2"]], 1L)
  expect_identical(ec$n$native[["C_gmm_cp08_pair_rows"]], 1L)
  expect_identical(ec$n$native[["C_gmm_cp08_dr_treatment_batch"]], 1L)
  expect_identical(ec$n$complete, 0L)
  expect_identical(sum(eg$n$native), 0L)
  acc <- fit_c$accounting$complete_support
  expect_identical(c(acc$n_treatment_logit_fits, acc$n_treatment_logit_reuses, acc$n_treatment_models),
    c(1L, T_count - 1L, T_count))
  expect_identical(fit_c$accounting$n_treatment_fits, T_count)
  expect_identical(acc$n_support_row_visits, sum(n_k))
  expect_identical(acc$n_support_chunks, sum(ceiling(n_k / fx$chunk)))
  expect_identical(c(acc$n_endpoint_searches, acc$n_alias_rescans), c(0, 0))
  # All T_count treatment identities and pair-specific quantities remain.
  expect_length(fit_c$treatments, T_count)
  expect_identical(length(unique(lapply(fit_c$treatments, `[[`, "support_ids"))), T_count)
  expect_identical(length(unique(lapply(fit_c$treatments, `[[`, "coefficients"))), 1L)
  od_expect_route_match(fit_c$treatments, fit_g$treatments)
  expect_identical(vapply(fit_c$links, `[[`, integer(1), "model"), seq_len(T_count))
  # One treatment replay per eligible chunk; one per model on the general route.
  # At p_or = 0 the unused OR core runs only for detailed output. Stage 4: the
  # eligible chunk is two native calls and reads no pair rows in R.
  ec$n$replay <- ec$n$core <- ec$n$complete <- 0L; ec$n$native[] <- 0L
  got_c <- ec$evaluate_dr_baseline_if_chunk(fit_c, 1:7, seq_len(T_count))
  expect_identical(c(ec$n$replay, ec$n$core, ec$n$complete), c(1L, 0L, 0L))
  expect_identical(unname(ec$n$native), c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L))
  eg$n$replay <- eg$n$core <- 0L
  got_g <- eg$evaluate_dr_baseline_if_chunk(fit_g, 1:7, seq_len(T_count))
  expect_identical(c(eg$n$replay, eg$n$core), c(T_count, 0L))
  od_expect_route_match(got_c$phi, got_g$phi)
  expect_identical(c(got_c$diagnostics$n_regression_score_evaluations, got_g$diagnostics$n_regression_score_evaluations),
    c(0L, 0L))
  ec$n$core <- eg$n$core <- 0L
  det_c <- ec$evaluate_dr_baseline_if_chunk(fit_c, 1:7, seq_len(T_count), details = TRUE)
  expect_identical(ec$n$core, 1L)
  det_g <- eg$evaluate_dr_baseline_if_chunk(fit_g, 1:7, seq_len(T_count), details = TRUE)
  expect_identical(eg$n$core, 1L)
  od_expect_route_match(det_c$details, det_g$details)
  expect_identical(c(det_c$diagnostics$n_regression_score_evaluations, det_g$diagnostics$n_regression_score_evaluations),
    rep(length(fit_g$regressions), 2L))
  # A replay whose models' inputs differ refuses instead of sharing.
  f2 <- fit_c$treatments[[2L]]; f2$coefficients <- f2$coefficients + 1
  expect_error(api$dr_baseline_complete_replay(fit_c, list(rr = 1:3, cohorts = fit_c$treatment_specs[[1L]]$cohorts,
    coefficients = fit_c$treatments[[1L]]$coefficients, v = c(1, 1)), fit_c$treatment_specs[[2L]], f2, 1:3,
    rep(1L, 3L), rep(1, 3L), matrix(1, 3L, 1L)), class = "cdid_dr_baseline_identity_error")
  share <- new.env(parent = emptyenv())
  spec <- fit_c$treatment_specs[[1L]]
  api$dr_baseline_shared_logit(share, spec, 1:4, matrix(1, 4L, 1L), c(1, 0, 1, 0), rep(1, 4L), c(1, 1), od_glm, 41L)
  expect_error(api$dr_baseline_shared_logit(share, spec, 1:4, matrix(1, 4L, 1L), c(1, 0, 1, 0), c(1, 1, 1, 2), c(1, 1),
    od_glm, 41L), class = "cdid_dr_baseline_identity_error")
  expect_identical(c(share$calls, share$reuses), c(1L, 0L))
})

test_that("independent references: raw-row intercept-only OR rows and the accepted dense DR table helper", {
  fx <- od_fx_perm; st <- fx$reader$state; N <- st$N
  fit_or <- od_or(fx); fit_dr <- od_dr(fx)
  y_at <- function(id, t) fx$d$y[match(paste(id, t), paste(fx$d$id, fx$d$t))]
  tg <- st$g == 3; nv <- is.infinite(st$g)
  phi_or <- od_batches(fit_or, api$replay_or_baseline_if, 16L, 10L)$batches
  phi_or <- do.call(rbind, lapply(phi_or, function(b) b$phi))
  for (j in seq_along(fx$pairs)) {
    p <- fx$pairs[[j]]; dy <- y_at(st$id, p[2]) - y_at(st$id, p[1])
    alpha <- mean(dy[nv]); mu <- mean(dy[tg]) - alpha
    ref <- ifelse(tg, N / sum(tg) * (dy - alpha - mu), -N / sum(nv) * (dy - alpha))
    expect_lt(abs(fit_or$regressions[[fit_or$links[[j]]$regression]]$alpha - alpha), 1e-10)
    expect_lt(abs(fit_or$links[[j]]$mu - mu), 1e-10)
    expect_lt(max(abs(phi_or[, j] - ref)), 1e-9)
  }
  Y <- sapply(seq_len(5L), function(t) y_at(st$id, t))
  links <- lapply(fx$pairs, function(p) list(s = p[1], t = p[2], g = 3, K = Inf, link_id = paste0("L", p[1], "_", p[2])))
  X1 <- matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)"))
  tb <- api$fit_dr_baseline_table(X_or = X1, V = X1, Y = Y, obs = matrix(1, N, 5L), g = st$g, w = st$w,
    calendar = 1:5, links = links, outcome_recipe = od_icpt, glm_control = od_glm, rank_tol = 1e-10)
  dr <- od_batches(fit_dr, api$replay_dr_baseline_if, 16L, 10L, details = TRUE)$batches
  col <- function(nm) do.call(rbind, lapply(dr, function(b) if (nm == "phi") b$phi else b$details[[nm]]))
  expect_lt(max(abs(vapply(fit_dr$links, `[[`, numeric(1), "mu_dr") - tb$estimates)), 2e-14)
  expect_lt(max(abs(vapply(fit_dr$links, `[[`, numeric(1), "mu_ipw") - tb$estimates_ipw)), 2e-14)
  for (nm in c("phi", "psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw"))
    expect_lt(max(abs(unname(col(nm)) - unname(tb[[nm]]))), 2e-14, label = nm)
})

test_that("route facts are actual: every broken fact keeps the general route without raising or reading", {
  d <- od_panel()
  wv <- 1 + (seq_len(41L) %% 3) / 4
  Xx <- cbind(`(Intercept)` = 1, x = seq_len(41L) %% 2)
  X_no_icpt <- matrix(seq_len(41L) %% 2 + 1, 41L, 1L, dimnames = list(NULL, "x"))
  X_bad_icpt <- cbind(`(Intercept)` = 2, x = seq_len(41L) %% 2)
  variants <- list(
    missing_endpoint = list(fx = od_fixture(d[!(d$id == 5 & d$t == 3), ])),
    varying_weights = list(fx = od_fixture(w = wv)),
    no_intercept_recipe = list(fx = od_fixture(X = X_no_icpt), recipe = api$or_pooled_outcome_recipe("x", FALSE, "x")),
    non_unit_intercept = list(fx = od_fixture(X = X_bad_icpt), recipe = api$or_pooled_outcome_recipe(c("(Intercept)", "x"), TRUE, "x")),
    frame_only_unit = list(fx = od_fixture(frame_only = TRUE)),
    other_pool = list(fx = od_fixture(od_panel(second = TRUE), K = c(4, Inf))),
    second_target = list(fx = od_fixture(od_panel(second = TRUE), target = c(3, 4))))
  for (nm in names(variants)) {
    v <- variants[[nm]]; rec <- if (is.null(v$recipe)) od_icpt else v$recipe
    for (method in c("or", "dr")) {
      ec <- od_counted()
      if (nm == "non_unit_intercept") {
        err_c <- tryCatch(if (method == "or") od_or(v$fx, ec, recipe = rec) else od_dr(v$fx, ec, recipe = rec), error = identity)
        err_g <- tryCatch(if (method == "or") od_or(v$fx, od_general, recipe = rec) else od_dr(v$fx, od_general, recipe = rec), error = identity)
        expect_s3_class(err_c, "cdid_error")
        expect_identical(class(err_c), class(err_g))
        expect_identical(conditionMessage(err_c), conditionMessage(err_g))
        expect_identical(sum(ec$n$native), 0L)
        next
      }
      fit <- if (method == "or") od_or(v$fx, ec, recipe = rec) else od_dr(v$fx, ec, recipe = rec)
      gen <- if (method == "or") od_or(v$fx, od_general, recipe = rec) else od_dr(v$fx, od_general, recipe = rec)
      cp09_eligible <- (nm == "other_pool" && method == "dr") || nm == "second_target"
      if (cp09_eligible) {
        expect_false(is.null(fit$complete_support), info = paste(nm, method))
        expect_null(gen$complete_support, info = paste(nm, method))
        expect_gt(sum(ec$n$native), 0L, label = paste(nm, method))
        if (method == "dr") {
          n_grp <- if (nm == "second_target") 2L else 1L
          expect_identical(ec$n$logit, n_grp, info = nm)
          expect_identical(ec$n$native[["C_gmm_cp08_dr_treatment_batch"]], n_grp, info = nm)
        }
        od_expect_route_match(od_backend_view(fit), od_backend_view(gen), tol = 1e-12, info = paste(nm, method))
      } else {
        expect_null(fit$complete_support, info = paste(nm, method))
        expect_identical(ec$n$complete, 0L, info = paste(nm, method))
        expect_identical(sum(ec$n$native), 0L, info = paste(nm, method))
        if (method == "dr") expect_identical(ec$n$logit, length(v$fx$views), info = nm)
        expect_identical(od_strip(unclass(fit)[names(fit) != "ref"]), od_strip(unclass(gen)[names(gen) != "ref"]),
          info = paste(nm, method))
      }
    }
  }
  # Common-design fixed X with intercept IS eligible for complete_support!
  fx_x <- od_fixture(X = Xx, V = Xx)
  rec_x <- api$or_pooled_outcome_recipe(c("(Intercept)", "x"), TRUE, "x")
  for (method in c("or", "dr")) {
    fit_x <- if (method == "or") od_or(fx_x, recipe = rec_x) else od_dr(fx_x, recipe = rec_x, treatment_columns = c("(Intercept)", "x"))
    gen_x <- if (method == "or") od_or(fx_x, od_general, recipe = rec_x) else od_dr(fx_x, od_general, recipe = rec_x, treatment_columns = c("(Intercept)", "x"))
    expect_false(is.null(fit_x$complete_support), info = method)
    expect_identical(fit_x$complete_support$p, 1L)
    od_expect_route_match(od_backend_view(fit_x), od_backend_view(gen_x), info = method)
  }
  # The predicate itself, on the eligible reader: a dated reference or an
  # out-of-range/non-intercept treatment column declines.
  ok <- od_or(od_fx)
  expect_false(is.null(ok$complete_support))
  rec_rows <- lapply(seq_along(ok$links), function(j) list(target = ok$links[[j]]$cohort,
    controls = ok$regression_specs[[ok$links[[j]]$regression]]$cohorts))
  expect_false(is.null(api$or_baseline_complete_facts(od_fx$reader, od_icpt, rec_rows, ok$support_specs, od_fx$ctrl)))
  expect_null(api$or_baseline_complete_facts(od_fx$reader, od_icpt, rec_rows, ok$support_specs, od_fx$ctrl,
    reference = list(dated = TRUE)))
  expect_null(api$or_baseline_complete_facts(od_fx$reader, od_icpt, rec_rows, ok$support_specs, od_fx$ctrl,
    treatment = list(index = 2L, p = 1L)))
})

test_that("fit, chunk, replay, producer and stream bounds refuse at the adjacent limit before any read", {
  fx <- od_fx; L <- length(fx$views)
  reseal <- function(fit, field, value, binding) {
    x <- fit; x$ref <- NULL; x$control[[field]] <- value; api$ipw_owner_seal(x, binding)
  }
  for (method in c("or", "dr")) {
    fit <- if (method == "or") od_or(fx) else od_dr(fx)
    binding <- if (method == "or") "or_baseline_backend" else "dr_baseline_backend"
    replay <- if (method == "or") api$replay_or_baseline_if else api$replay_dr_baseline_if
    budget <- if (method == "or") api$or_baseline_if_budget else api$dr_baseline_if_budget
    W <- fit$accounting$work_ops_bound
    ec <- od_counted()
    fit_at <- if (method == "or") od_or(fx, ec, od_ctrl(7L, work = W)) else od_dr(fx, ec, od_ctrl(7L, work = W))
    expect_identical(fit_at$accounting$work_ops_bound, W)
    ec <- od_counted()
    expect_error(if (method == "or") od_or(fx, ec, od_ctrl(7L, work = W - 1)) else od_dr(fx, ec, od_ctrl(7L, work = W - 1)),
      class = "cdid_ipw_owner_bound_error")
    expect_identical(c(ec$n$complete, ec$n$search, ec$n$logit, sum(ec$n$native)), c(0L, 0L, 0L, 0L))
    # Cumulative replay work: admitted at its total, refused one below before any callback.
    total <- replay(fit, 7L, L, function(b) NULL)$work_ops_bound
    calls <- 0L
    invisible(replay(reseal(fit, "max_work_ops", total, binding), 7L, L, function(b) { force(b); calls <<- calls + 1L }))
    expect_identical(calls, as.integer(ceiling(41 / 7)))
    calls <- 0L
    expect_error(replay(reseal(fit, "max_work_ops", total - 1, binding), 7L, L, function(b) { force(b); calls <<- calls + 1L }),
      class = "cdid_ipw_owner_bound_error")
    expect_identical(calls, 0L)
    # Concurrent chunk payload: retained owner entries plus the live buffers.
    pre <- budget(fit, 7L, seq_len(L))
    expect_identical(pre$complete$route, "complete_support")
    edge <- fit$accounting$owned_entries + fit$accounting$borrowed_entries + pre$scratch
    calls <- 0L
    invisible(replay(reseal(fit, "max_matrix_entries", edge, binding), 7L, L, function(b) { force(b); calls <<- calls + 1L }))
    expect_gt(calls, 0L)
    calls <- 0L
    expect_error(replay(reseal(fit, "max_matrix_entries", edge - 1, binding), 7L, L, function(b) { force(b); calls <<- calls + 1L }),
      class = "cdid_ipw_owner_bound_error")
    expect_identical(calls, 0L)
    # Producer plus downstream simultaneous entries at the public preflight.
    own <- fit$accounting$owned_entries + fit$accounting$borrowed_entries
    room <- fit$control$max_matrix_entries - own - pre$output - pre$scratch
    ok <- api$cp04_preflight_links(fit, seq_len(L), retained = room)
    expect_identical(ok$chunk_rows, 7)
    expect_identical(ok$work_ops_bound, ceiling(41 / 7) * pre$work)
    expect_error(api$cp04_preflight_links(fit, seq_len(L), retained = room + 1), class = "cdid_ipw_owner_bound_error")
    # The legacy covariance plan composes both passes of every block pair;
    # D92 cp04_stream_plan admits a single full-width block when capacity and work allow.
    legacy_plan <- api$cp04_stream_plan_legacy(fit, seq_len(L))
    blocks <- split(seq_len(L), ceiling(seq_len(L) / min(L, 7L)))
    expect_identical(length(legacy_plan$blocks), length(blocks))
    w <- 0
    for (a in seq_along(blocks)) for (b in a:length(blocks)) {
      ids <- unique(c(blocks[[a]], blocks[[b]]))
      w <- w + 2 * ceiling(41 / legacy_plan$chunk_rows) * budget(fit, legacy_plan$chunk_rows, ids)$work + 4 * 41 * length(ids)^2
    }
    expect_identical(legacy_plan$work, w)
    plan <- api$cp04_stream_plan(fit, seq_len(L))
    expect_identical(plan$blocks, split(seq_len(L), rep.int(1L, L)))
    w_full <- 2 * ceiling(41 / plan$chunk_rows) * budget(fit, plan$chunk_rows, seq_len(L))$work +
      4 * 41 * L^2 + api$cp04_retained_if_boundary_cost(fit, seq_len(L), plan$chunk_rows)$work
    expect_identical(plan$work, w_full)
    expect_lt(plan$work, legacy_plan$work)
    expect_identical(plan$chunk_rows, 7)
  }
})

test_that("complete-support bindings are re-verified at invocation and replay boundaries", {
  fx <- od_fixture()
  rebind <- function(env, nm, value) {
    old <- get(nm, envir = env)
    unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
    function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
  }
  po <- fx$po
  shifted <- po; shifted$postings$date_start[2L] <- shifted$postings$date_start[2L] + 1L
  swapped <- po; k <- po$postings$date_start[3L] - 1L + c(1L, 2L)
  swapped$postings$post_a_index[k] <- swapped$postings$post_a_index[rev(k)]
  for (method in c("or", "dr")) {
    fit <- if (method == "or") od_or(fx) else od_dr(fx)
    replay <- if (method == "or") api$replay_or_baseline_if else api$replay_dr_baseline_if
    err <- if (method == "or") "cdid_or_baseline_integrity_error" else "cdid_dr_baseline_integrity_error"
    base <- od_batches(fit, replay, 7L, 4L)$batches
    for (bad in list(shifted, swapped)) {
      restore <- rebind(po$ref, "panel_owner", bad)
      calls <- 0L
      expect_error(replay(fit, 7L, 4L, function(b) { force(b); calls <<- calls + 1L }), class = err)
      expect_identical(calls, 0L)
      restore()
    }
    # A callback that rebinds the postings stops the next chunk before it reads.
    # (The OR replay hands each chunk to its callback as a promise; force it.)
    restore <- NULL; calls <- 0L
    expect_error(replay(fit, 7L, 4L, function(b) {
      force(b)
      calls <<- calls + 1L
      if (is.null(restore)) restore <<- rebind(po$ref, "panel_owner", shifted)
    }), class = err)
    restore()
    expect_identical(calls, 1L)
    expect_identical(od_batches(fit, replay, 7L, 4L)$batches, base)
  }
})

test_that("first genuine failures keep their class and identity on both routes", {
  fx <- od_fx
  one <- list(epsilon = 1e-13, maxit = 1L, trace = FALSE)
  ec <- tryCatch(od_dr(fx, glm = one), error = identity)
  eg <- tryCatch(od_dr(fx, od_general, glm = one), error = identity)
  expect_s3_class(ec, "cdid_cp04_convergence_error")
  expect_identical(class(ec), class(eg))
  expect_identical(conditionMessage(ec), conditionMessage(eg))
  expect_identical(ec$diagnostics, eg$diagnostics)
  # Finite stored endpoints whose raw difference overflows refuse identically.
  d <- od_panel(); d$y[d$id == 4 & d$t == 2] <- 1e308; d$y[d$id == 4 & d$t == 3] <- -1e308
  bad <- od_fixture(d)
  for (method in c("or", "dr")) {
    ec <- tryCatch(if (method == "or") od_or(bad) else od_dr(bad), error = identity)
    eg <- tryCatch(if (method == "or") od_or(bad, od_general) else od_dr(bad, od_general), error = identity)
    expect_s3_class(ec, "cdid_error")
    expect_identical(class(ec), class(eg))
    expect_identical(conditionMessage(ec), conditionMessage(eg))
  }
})

test_that("complete-support fact scans are admitted before reading bounded chunks", {
  fx <- od_fixture(d = od_panel(N = 81L), pairs = list(c(1, 2)))
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline", "or_baseline_complete_facts")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  scans <- 0L
  e$ipw_owner_posting_date_complete <- function(...) {
    scans <<- scans + 1L; api$ipw_owner_posting_date_complete(...)
  }
  # Two dates, the weight/design/eligibility pass, and one link. This is the
  # already recorded fact-work term, now checked before any posting scan.
  W <- 64 * (2 * 81 + 2 * 81 + 8)
  for (method in c("or", "dr")) {
    scans <- 0L
    expect_error(if (method == "or") od_or(fx, e, od_ctrl(work = W - 1)) else
      od_dr(fx, e, od_ctrl(work = W - 1)), class = "cdid_ipw_owner_bound_error")
    expect_identical(scans, 0L)
  }
  fit <- od_or(fx)
  records <- list(list(target = fit$links[[1L]]$cohort,
    controls = fit$regression_specs[[1L]]$cohorts))
  facts <- function(ctrl, ...) e$or_baseline_complete_facts(fx$reader, od_icpt, records,
    fit$support_specs, ctrl, ...)
  scans <- 0L
  expect_false(is.null(facts(od_ctrl(work = W))))
  expect_identical(scans, 2L)
  scans <- 0L
  expect_error(facts(od_ctrl(work = W - 1)), class = "cdid_ipw_owner_bound_error")
  expect_identical(scans, 0L)
  own <- fx$reader$accounting$owned_entries + fx$reader$accounting$borrowed_owner_entries
  edge <- own + 16 * 7 + 8 * 2 + 1024
  expect_false(is.null(facts(od_ctrl(matrix = edge))))
  scans <- 0L
  expect_error(facts(od_ctrl(matrix = edge - 1)), class = "cdid_ipw_owner_bound_error")
  expect_identical(scans, 0L)
  scans <- 0L
  expect_error(facts(od_ctrl(matrix = edge), retained_external_entries = 1),
    class = "cdid_ipw_owner_bound_error")
  expect_identical(scans, 0L)
  scans <- 0L
  expect_error(facts(od_ctrl(work = W), external_work_ops = 1), class = "cdid_ipw_owner_bound_error")
  expect_identical(scans, 0L)
})

# ---------------------------------------------------------------------------
# CP08 stage 6: post-fit complete-support p = 0 structural reduction for
# baseline OR/DR under optimal GMM weighting and PT-post.
# ---------------------------------------------------------------------------

od_nonadd_panel <- function(N = 41L, T = 5L, n_g = 18L, seed = 31L, id_perm = FALSE, rows_perm = FALSE) {
  set.seed(seed)
  g <- rep(0, N)
  g[sample.int(N, n_g)] <- 3
  ids <- if (id_perm) sample(1000L + seq_len(5L * N), N) else seq_len(N)
  xi <- stats::rbinom(N, 1, 0.5)
  ue <- stats::rnorm(N)
  d <- data.frame(id = rep(ids, each = T), period = rep(seq_len(T), N))
  u <- match(d$id, ids)
  d$G <- g[u]
  d$x <- xi[u]
  eff <- ifelse(d$G > 0 & d$period >= d$G, 1 + 0.25 * (d$period - d$G), 0)
  # Nonadditive in (unit, period) with nonlinear interaction plus noise.
  d$y <- ue[u] + 0.1 * d$period + 0.35 * sin(0.7 * u * d$period) + 0.2 * d$x * (d$period^2) + eff + stats::rnorm(nrow(d))
  d$cl <- (u - 1L) %/% 4L + 1L
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

od_run_pipeline <- function(d, est = "reg", chunk = 7L, reduce = TRUE, weighting = "optimal",
                            pt_assumption = "post", pt_start = NULL, xformla = NULL,
                            weightsname = NULL, control_group = "nevertreated",
                            target_cohorts = NULL, frame = NULL, bstrap = FALSE,
                            cband = FALSE, biters = 99L, seed = 101L,
                            aggregations = c("simple", "group", "dynamic", "calendar"),
                            condition_trigger = NULL, ridge_relative = NULL, ns = api) {
  ctrl_list <- if (is.null(chunk)) list() else list(chunk_rows = chunk)
  if (!is.null(condition_trigger)) ctrl_list$condition_trigger <- condition_trigger
  if (!is.null(ridge_relative)) ctrl_list$ridge_relative <- ridge_relative
  if (!length(ctrl_list)) ctrl_list <- NULL
  resolved <- api$chained_public_resolve_control(ctrl_list, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  prep <- api$chained_request_prepare_inputs(
    data = d, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xformla, weightsname = weightsname, clustervars = "cl", anticipation = 0,
    frame = frame, outcome_regression = identical(est, "reg"),
    doubly_robust = identical(est, "dr"), or_resource = resolved$profile
  )
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = biters, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  tc <- if (is.null(target_cohorts)) prep$target_cohorts else target_cohorts
  ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = built$control, xformla = prep$xformla, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = tc, anticipation = 0L, control_group = control_group, base_period = "varying",
    weighting = weighting, reduce_structural_redundancy = reduce, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = aggregations, cband = cband, alp = 0.05, biters = biters, seed = seed,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = bstrap, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = pt_assumption, pt_start = pt_start, cp04 = prep$cp04
  )
}

od_dense_E <- function(cert) {
  E <- matrix(0, cert$n_links, length(cert$retained_index))
  E[cbind(cert$retained_index, seq_along(cert$retained_index))] <- 1
  if (length(cert$E$i)) {
    rpos <- match(cert$E$j, cert$retained_index)
    for (k in seq_along(cert$E$i)) {
      E[cert$E$i[k], rpos[k]] <- E[cert$E$i[k], rpos[k]] + cert$E$x[k]
    }
  }
  E
}

test_that("stage 6 local case-weight functional identity and complete-F IF forest relations hold on fitted OR and DR", {
  for (perm in c(FALSE, TRUE)) {
    d <- od_nonadd_panel(41L, 5L, 18L, seed = 31L, id_perm = perm, rows_perm = perm)
    chunk <- if (perm) 16L else 7L
    for (est in c("reg", "dr")) {
      res_r <- od_run_pipeline(d, est = est, chunk = chunk, reduce = TRUE)
      res_u <- od_run_pipeline(d, est = est, chunk = chunk, reduce = FALSE)
      cert <- res_r$structural_certificate
      expect_false(is.null(cert))
      expect_false(inherits(cert, "cdid_ipw_structural_certificate"))
      expect_true(api$ipw_owner_seal_ok(cert, "cp04_complete_reduction"))
      expect_identical(cert$n_links, 10L)
      expect_identical(length(cert$retained_index), 4L)
      expect_identical(length(cert$dropped_index), 6L)
      expect_identical(cert$link_kind, rep(if (est == "dr") "complete_support_DR" else "complete_support_OR", 10L))
      expect_identical(res_r$gmm_fit$provenance$reduction$mode, "sparse_class_certificate")

      # Strict IPW certificate guard is preserved on apply_ipw_reduction_map.
      expect_error(api$apply_ipw_reduction_map(cert, rep(1, 4L)), class = "cdid_ipw_certificate_error")

      # Full backend fits, point_links_all, auxiliary reporting and public diagnostics are preserved bitwise.
      b_r <- res_r$nuisance_registry$cp04
      b_u <- res_u$nuisance_registry$cp04
      expect_identical(od_backend_view(b_r), od_backend_view(b_u))
      expect_identical(res_r$point_links_all, res_u$point_links_all)
      expect_identical(res_r$nuisance_registry$cp04_points, res_u$nuisance_registry$cp04_points)
      expect_identical(res_r$nuisance_registry$cp04_aux_index, res_u$nuisance_registry$cp04_aux_index)
      expect_identical(od_strip(res_r$auxiliary_reporting), od_strip(res_u$auxiliary_reporting))
      diag_r <- api$cp04_public_diagnostics(res_r, res_r$controls$pipeline_control)
      diag_u <- api$cp04_public_diagnostics(res_u, res_u$controls$pipeline_control)
      expect_identical(od_strip(diag_r), od_strip(diag_u))

      # Independent case-weight functional and Gateaux derivative on original F unit order.
      reader <- res_r$nuisance_registry$cp04_reader
      st <- reader$state
      N <- st$N
      tg <- st$g == 3
      nv <- is.infinite(st$g)
      pi_g <- sum(tg) / N
      pi_c <- sum(nv) / N
      y_at <- function(id, t) d$y[match(paste(id, t), paste(d$id, d$period))]

      m_cw <- numeric(10L)
      H_cw <- matrix(0, N, 10L)
      H_fd <- matrix(0, N, 10L)
      eps_fd <- 1e-6
      for (m in seq_len(10L)) {
        pr <- res_r$point_links_all[[m]]$pair
        dy <- y_at(st$id, pr$t) - y_at(st$id, pr$s)
        mu_plus <- mean(dy[tg])
        mu_minus <- mean(dy[nv])
        m_cw[m] <- mu_plus - mu_minus
        H_cw[, m] <- ifelse(tg, (dy - mu_plus) / pi_g, -((dy - mu_minus) / pi_c))
        for (i in seq_len(N)) {
          w_p <- rep(1, N); w_p[i] <- 1 + eps_fd * N
          w_m <- rep(1, N); w_m[i] <- 1 - eps_fd * N
          f_p <- sum(w_p[tg] * dy[tg]) / sum(w_p[tg]) - sum(w_p[nv] * dy[nv]) / sum(w_p[nv])
          f_m <- sum(w_m[tg] * dy[tg]) / sum(w_m[tg]) - sum(w_m[nv] * dy[nv]) / sum(w_m[nv])
          H_fd[i, m] <- (f_p - f_m) / (2 * eps_fd)
        }
      }
      expect_lt(max(abs(H_cw - H_fd)), 1e-9)

      # Actual fitted moments and complete-F IFs (core scale = frame scale since n_A == N).
      m_all <- vapply(res_r$point_links_all, `[[`, numeric(1L), "estimate")
      H_all <- matrix(0, N, 10L)
      pos <- 0L
      api$cp04_replay_links(b_r, seq_len(10L), function(batch) {
        idx <- seq.int(pos + 1L, pos + nrow(batch$phi))
        H_all[idx, ] <<- batch$phi
        pos <<- pos + nrow(batch$phi)
      }, chunk_rows = chunk)

      expect_lt(max(abs(m_all - m_cw)), 1e-13)
      expect_lt(max(abs(H_all - H_cw)), 1e-12)

      if (est == "dr") {
        # Verify intercept propensity odds, zero selection/OR scores and propensity centering cancellation.
        psi_dir <- matrix(0, N, length(b_r$links))
        max_sel <- max_or <- max_ps <- 0
        api$replay_dr_baseline_if(b_r, chunk, chunk, function(b) {
          psi_dir[b$rows, b$links] <<- b$details$psi_dir
          max_sel <<- max(max_sel, max(abs(b$details$psi_sel)))
          max_or <<- max(max_or, max(abs(b$details$psi_or)))
          max_ps <<- max(max_ps, max(abs(b$details$psi_ps)))
        }, scale = "frame", details = TRUE)
        expect_identical(max_sel, 0)
        expect_identical(max_or, 0)
        expect_lt(max_ps, 1e-14)
        expect_lt(max(abs(psi_dir[, 1:10] - H_cw)), 1e-12)
        p_hat <- 1 / (1 + exp(-b_r$treatments[[1L]]$coefficients[1L]))
        expect_lt(abs(p_hat - sum(tg) / N), 1e-13)
      }

      # Exact forest algebraic relations on actual fitted quantities: m = E m_R, H = H_R E', V = E V_R E'.
      E <- od_dense_E(cert)
      m_R <- m_all[cert$retained_index]
      H_R <- H_all[, cert$retained_index, drop = FALSE]
      V_all <- crossprod(sweep(H_all, 2L, colMeans(H_all), "-")) / N
      V_R <- crossprod(sweep(H_R, 2L, colMeans(H_R), "-")) / N
      expect_lt(max(abs(m_all - drop(E %*% m_R))), 1e-13)
      expect_lt(max(abs(H_all - H_R %*% t(E))), 1e-12)
      expect_lt(max(abs(V_all - E %*% V_R %*% t(E))), 1e-12)
      expect_lt(res_r$reduction_summary$max_moment_expansion_diff, 1e-13)
      expect_equal(api$cp04_expand_reduction_map(cert, m_R), drop(E %*% m_R), tolerance = 1e-15)
      expect_equal(t(api$cp04_expand_reduction_map(cert, t(H_R))), H_R %*% t(E), tolerance = 1e-15)
    }
  }
})

test_that("stage 6 full-metric GMM, 4 aggregations, uniform bands, pretests and 190-to-19 column reduction work", {
  d <- od_nonadd_panel(41L, 5L, 18L, seed = 31L, id_perm = TRUE, rows_perm = TRUE)
  for (est in c("reg", "dr")) {
    # Both pointwise and multiplier bootstrap uniform bands with fixed seed.
    res_r <- od_run_pipeline(d, est = est, chunk = 16L, reduce = TRUE, bstrap = TRUE, cband = TRUE, biters = 99L, seed = 42L)
    res_u <- od_run_pipeline(d, est = est, chunk = 16L, reduce = FALSE, bstrap = TRUE, cband = TRUE, biters = 99L, seed = 42L)
    # On N = 41, T = 5, reduced 4 x 4 transformed covariance B = C V_R C' is well-conditioned PD ("none"),
    # while unreduced 10 x 10 covariance has 6 structural cycle zeros and triggers 1e-8 Tikhonov ridge.
    expect_identical(res_r$covariance_handling$handler_used, "none")
    expect_identical(res_r$covariance_handling$trigger_reason, "well_conditioned_pd")
    expect_identical(res_r$covariance_handling$gram_provenance, "sparse_structural_reduction_E_t_E")
    expect_identical(res_r$covariance_handling$dimensions$M_full, 10L)
    expect_identical(res_r$covariance_handling$dimensions$L_retained, 4L)
    expect_equal(res_r$covariance_handling$trace_B, res_u$covariance_handling$trace_B, tolerance = 1e-12)
    expect_equal(res_r$covariance_handling$scale, res_u$covariance_handling$scale, tolerance = 1e-12)

    # ATT, SE, full-F IF, vcov, uniform bands, all 4 aggregations and unconditional pretest match to machine precision.
    expect_equal(res_r$reporting$estimates, res_u$reporting$estimates, tolerance = 1e-13)
    expect_equal(res_r$reporting$se, res_u$reporting$se, tolerance = 1e-13)
    expect_equal(res_r$reporting$analytical_se, res_u$reporting$analytical_se, tolerance = 1e-13)
    expect_equal(res_r$reporting$if_full, res_u$reporting$if_full, tolerance = 1e-12)
    expect_equal(res_r$reporting$vcov, res_u$reporting$vcov, tolerance = 1e-13)
    expect_equal(res_r$reporting$cband_lower, res_u$reporting$cband_lower, tolerance = 1e-13)
    expect_equal(res_r$reporting$cband_upper, res_u$reporting$cband_upper, tolerance = 1e-13)
    for (ag in c("simple", "group", "dynamic", "calendar")) {
      expect_equal(res_r$aggregations[[ag]]$overall$estimate, res_u$aggregations[[ag]]$overall$estimate, tolerance = 1e-13)
      expect_equal(res_r$aggregations[[ag]]$overall$se, res_u$aggregations[[ag]]$overall$se, tolerance = 1e-13)
      expect_equal(res_r$aggregations[[ag]]$overall$if_full, res_u$aggregations[[ag]]$overall$if_full, tolerance = 1e-12)
      expect_equal(res_r$aggregations[[ag]]$estimates, res_u$aggregations[[ag]]$estimates, tolerance = 1e-13)
      expect_equal(res_r$aggregations[[ag]]$se, res_u$aggregations[[ag]]$se, tolerance = 1e-13)
      expect_equal(res_r$aggregations[[ag]]$if_full, res_u$aggregations[[ag]]$if_full, tolerance = 1e-12)
    }
    pre_r <- api$chained_pretest_record(res_r$reporting, max_matrix_entries = 1e7)
    pre_u <- api$chained_pretest_record(res_u$reporting, max_matrix_entries = 1e7)
    expect_identical(pre_r$status, "computed")
    expect_equal(pre_r$W, pre_u$W, tolerance = 1e-13)
    expect_equal(pre_r$Wpval, pre_u$Wpval, tolerance = 1e-13)
  }

  # Demonstrate actual 190-to-19 replay column reduction on a T = 20 balanced panel.
  d20 <- od_nonadd_panel(45L, 20L, 20L, seed = 77L)
  for (est in c("reg", "dr")) {
    e_cnt <- new.env(parent = api)
    for (nm in c("run_chained_ipw_pipeline", "cp04_replay_links")) {
      f <- get(nm, api); environment(f) <- e_cnt; assign(nm, f, e_cnt)
    }
    cov_fn <- get("cp04_stream_covariance", api); environment(cov_fn) <- e_cnt
    cov_cols <- integer()
    rep_cols <- integer()
    e_cnt$cp04_stream_covariance <- function(backend, links, retained = 0) {
      cov_cols <<- c(cov_cols, length(links))
      cov_fn(backend, links, retained)
    }
    e_cnt$cp04_replay_links <- function(backend, links, callback, scale = "core", chunk_rows = NULL) {
      rep_cols <<- c(rep_cols, length(links))
      api$cp04_replay_links(backend, links, callback, scale = scale, chunk_rows = chunk_rows)
    }
    res20 <- od_run_pipeline(d20, est = est, chunk = 50L, reduce = TRUE, aggregations = "simple", ns = e_cnt)
    expect_identical(res20$reduction_summary$n_links_all, 190L)
    expect_identical(res20$reduction_summary$n_links_used, 19L)
    expect_true(res20$reduction_summary$reduction_active)
    expect_identical(cov_cols, 19L)
    # 2 covariance passes (19 cols each) + 2 GMM IF projection passes (19 cols each).
    expect_identical(rep_cols, rep(19L, 4L))
  }
})

test_that("stage 6 ridge-active and zero-covariance witnesses preserve full-catalog G = E'E and refusal semantics", {
  # Ridge-active witness: N = 6 units (3 treated, 3 never-treated), T = 5 dates (M = 10, L = 4).
  # On 6 units with 2 cohorts, the 4 retained moment IFs lie in a 4-dimensional within-cohort subspace,
  # and with condition_trigger = 0.99 both reduced and unreduced systems trigger Tikhonov ridge with identical lambda.
  d_ridge <- od_nonadd_panel(12L, 5L, 6L, seed = 55L)
  for (est in c("reg", "dr")) {
    res_r <- od_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = TRUE, condition_trigger = 0.99)
    res_u <- od_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = FALSE, condition_trigger = 0.99)
    ch_r <- res_r$covariance_handling
    ch_u <- res_u$covariance_handling
    expect_identical(ch_r$handler_used, "tikhonov")
    expect_identical(ch_u$handler_used, "tikhonov")
    expect_identical(ch_r$dimensions$M_full, 10L)
    expect_identical(ch_r$dimensions$L_retained, 4L)
    expect_false(ch_r$gram_is_identity)
    expect_false(ch_r$gram_is_diagonal)
    expect_equal(ch_r$trace_B, ch_u$trace_B, tolerance = 1e-13)
    expect_equal(ch_r$scale, ch_u$scale, tolerance = 1e-13)
    expect_equal(ch_r$lambda, ch_u$lambda, tolerance = 1e-13)
    # With identical Tikhonov lambda and G = E'E, reduced and unreduced full-criterion solutions agree to machine precision!
    expect_equal(res_r$reporting$estimates, res_u$reporting$estimates, tolerance = 1e-12)
    expect_equal(res_r$reporting$se, res_u$reporting$se, tolerance = 1e-12)
    expect_equal(res_r$reporting$if_full, res_u$reporting$if_full, tolerance = 1e-12)
    expect_equal(res_r$reporting$vcov, res_u$reporting$vcov, tolerance = 1e-12)
    for (ag in c("simple", "group", "dynamic", "calendar")) {
      expect_equal(res_r$aggregations[[ag]]$overall$estimate, res_u$aggregations[[ag]]$overall$estimate, tolerance = 1e-12)
      expect_equal(res_r$aggregations[[ag]]$overall$se, res_u$aggregations[[ag]]$overall$se, tolerance = 1e-12)
    }

    # Verify that an unweighted reduced solve (G = I_L, M_full = L) would disagree strongly under ridge_relative = 0.5.
    res_r_big <- od_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = TRUE, condition_trigger = 0.99, ridge_relative = 0.5)
    res_u_big <- od_run_pipeline(d_ridge, est = est, chunk = 7L, reduce = FALSE, condition_trigger = 0.99, ridge_relative = 0.5)
    expect_equal(res_r_big$reporting$estimates, res_u_big$reporting$estimates, tolerance = 1e-12)
    expect_equal(res_r_big$reporting$se, res_u_big$reporting$se, tolerance = 1e-12)
    # Solve the unweighted reduced problem (mode = "none" on the 4 retained moments) and confirm it differs!
    nuis_id <- "stage6_unw_nuis"
    scale_id <- "stage6_unw_scale"
    prov_unw <- api$gmm_covariance_provenance(
      row_level = "unit", n_units = 12L, n_clusters = 12L, n_positive_contributing_rows = 12L,
      centered = TRUE, weighted = FALSE, weight_type = "unweighted", normalization_factor = 1 / 12,
      scale_id = scale_id, nuisance_fit_id = nuis_id,
      reduction = list(mode = "none", reduction_active = FALSE)
    )
    prov_unw$complete_unit_if_gram <- api$ipw_seal(
      list(covariance = res_r_big$covariance_stream$V, n_units = 12L, nuisance_fit_id = nuis_id),
      "complete_unit_if_gram"
    )
    gmm_unw <- api$solve_chained_gmm_core(
      design = res_r_big$design_used,
      mhat = res_r_big$covariance_stream$mhat,
      V = res_r_big$covariance_stream$V,
      provenance = prov_unw,
      targets = api$ipw_pipeline_core_target_plan(res_r_big$support$ledger, res_r_big$design_used)$target_specs,
      weighting = "optimal",
      if_replay = api$gmm_if_replay_spec(
        callback = function(consume) {
          cp04_used <- vapply(res_r_big$point_links_used, `[[`, integer(1), "cp04_index")
          api$cp04_replay_links(res_r_big$nuisance_registry$cp04, cp04_used, function(b) {
            phi <- b$phi; colnames(phi) <- res_r_big$design_used$incidence$link_ids
            consume(list(unit_id = b$id, phi = phi, weight = rep(1, nrow(phi)), nuisance_fit_id = prov_unw$nuisance_fit_id))
          }, chunk_rows = 7L)
        },
        input_row_level = "unit", aggregate_to_clusters = FALSE, n_input_rows = 12L, n_output_rows = 12L,
        n_positive_contributing_rows = 12L, centered = TRUE, weighted = FALSE, weight_type = "unweighted",
        normalization_factor = 1 / 12, scale_id = prov_unw$scale_id, nuisance_fit_id = prov_unw$nuisance_fit_id,
        expected_unit_ids = res_r_big$nuisance_registry$cp04_reader$state$id
      ),
      accuracy = res_r_big$controls$pipeline_control$gmm_accuracy,
      bounds = res_r_big$controls$pipeline_control$gmm_bounds,
      gauge_rule = "lowest", request_param_cov = TRUE, materialize_target_if = TRUE
    )
    expect_equal(res_r_big$gmm_fit$targets$se_criterion_conditional, res_u_big$gmm_fit$targets$se_criterion_conditional, tolerance = 1e-12)
    expect_gt(abs(gmm_unw$covariance_handling$trace_B - res_r_big$covariance_handling$trace_B), 1e-1)
    expect_gt(abs(gmm_unw$covariance_handling$scale - res_r_big$covariance_handling$scale), 1e-1)
    expect_gt(abs(gmm_unw$covariance_handling$lambda - res_r_big$covariance_handling$lambda), 1e-1)
    se_unw <- gmm_unw$targets$se_criterion_conditional
    se_full <- res_r_big$gmm_fit$targets$se_criterion_conditional
    expect_gt(length(se_unw), 0L)
    expect_identical(length(se_unw), length(se_full))
    expect_identical(is.finite(se_unw), is.finite(se_full))
    finite_se <- is.finite(se_unw) & is.finite(se_full)
    expect_true(any(finite_se))
    expect_gt(max(abs(se_unw[finite_se] - se_full[finite_se])), 1e-2)
  }

  # Zero-covariance witnesses:
  # (a) Identically zero outcomes -> V_R == 0 and V_all == 0 -> triggers exact "zero_covariance" branch.
  d_zero0 <- od_nonadd_panel(12L, 5L, 6L, seed = 55L)
  d_zero0$y <- 0
  # (b) Constant non-zero treatment effect with zero within-cohort variance -> refuses with cdid_ipw_pipeline_scale_error.
  d_zero1 <- od_nonadd_panel(12L, 5L, 6L, seed = 55L)
  d_zero1$y <- ifelse(d_zero1$G == 3 & d_zero1$period >= 3, 2 * (d_zero1$period - 2), 0)
  for (est in c("reg", "dr")) {
    z_r <- od_run_pipeline(d_zero0, est = est, chunk = 7L, reduce = TRUE)
    z_u <- od_run_pipeline(d_zero0, est = est, chunk = 7L, reduce = FALSE)
    expect_identical(z_r$covariance_handling$trigger_reason, "zero_covariance")
    expect_identical(z_u$covariance_handling$trigger_reason, "zero_covariance")
    expect_identical(z_r$covariance_handling$handler_used, "tikhonov")
    expect_identical(z_r$covariance_handling$lambda, 0)
    expect_equal(z_r$reporting$estimates, z_u$reporting$estimates, tolerance = 1e-13)

    err_r <- tryCatch(od_run_pipeline(d_zero1, est = est, chunk = 7L, reduce = TRUE), error = identity)
    err_u <- tryCatch(od_run_pipeline(d_zero1, est = est, chunk = 7L, reduce = FALSE), error = identity)
    expect_s3_class(err_r, "cdid_ipw_pipeline_scale_error")
    expect_identical(class(err_r), class(err_u))
  }
})

test_that("stage 6 signed-E mutation and seal tampering fail loudly", {
  d <- od_nonadd_panel(41L, 5L, 18L, seed = 31L)
  res <- od_run_pipeline(d, est = "dr", chunk = 7L, reduce = TRUE)
  cert <- res$structural_certificate
  m_all <- vapply(res$point_links_all, `[[`, numeric(1L), "estimate")
  m_R <- m_all[cert$retained_index]

  # Unsealed tampering of E$x is caught by cp04_expand_reduction_map seal verification.
  tampered <- cert
  tampered$E$x[1L] <- -tampered$E$x[1L]
  expect_error(api$cp04_expand_reduction_map(tampered, m_R), class = "cdid_ipw_certificate_error")

  # Re-sealed mutated sign produces a large moment expansion error (> 0.1) and alters G = E'E.
  mut_sealed <- tampered
  mut_sealed$ref <- NULL
  mut_sealed <- api$ipw_owner_seal(mut_sealed, "cp04_complete_reduction")
  m_mut <- api$cp04_expand_reduction_map(mut_sealed, m_R)
  expect_gt(max(abs(m_all - m_mut)), 0.1)
  g_orig <- api$assemble_gmm_gram_matrix(4L, "optimal", list(mode = "sparse_class_certificate",
    class_certificate = cert, kept_indices = res$reduction_summary$kept_indices, n_links_all = 10L, reduction_active = TRUE))
  g_mut <- api$assemble_gmm_gram_matrix(4L, "optimal", list(mode = "sparse_class_certificate",
    class_certificate = mut_sealed, kept_indices = res$reduction_summary$kept_indices, n_links_all = 10L, reduction_active = TRUE))
  expect_false(identical(g_orig$G, g_mut$G))
})

test_that("stage 6 negative eligibility cases keep the unreduced path and preserve weight-rate refusals", {
  d0 <- od_nonadd_panel(41L, 5L, 18L, seed = 31L)
  dm <- d0[!(d0$id == 5L & d0$period == 3L), ]
  dw <- d0; dw$w <- 1 + (dw$id %% 3) / 4
  d2 <- d0; d2$G[d2$id %in% c(4L, 9L, 14L, 19L)] <- 4
  fr <- d0[!duplicated(d0$id), c("id", "G", "x", "cl")]
  fr <- rbind(fr[order(fr$id), ], data.frame(id = 9999L, G = 0, x = 0, cl = 99L))

  neg_cases <- list(
    missing_date = function() od_run_pipeline(dm, est = "reg"),
    covariate_x_or = function() od_run_pipeline(d0, est = "reg", xformla = ~x),
    covariate_x_dr = function() od_run_pipeline(d0, est = "dr", xformla = ~x),
    varying_weights = function() od_run_pipeline(dw, est = "dr", weightsname = "w"),
    frame_only_unit = function() od_run_pipeline(d0, est = "reg", frame = fr),
    multiple_targets = function() od_run_pipeline(d2, est = "dr", control_group = "notyettreated"),
    identity_weighting = function() od_run_pipeline(d0, est = "dr", weighting = "identity"),
    pt_all = function() od_run_pipeline(d0, est = "dr", pt_assumption = "all", pt_start = 1)
  )
  for (nm in names(neg_cases)) {
    res <- neg_cases[[nm]]()
    if (nm %in% c("covariate_x_or", "covariate_x_dr")) {
      # Authorized shared-X extension: these same fixtures now have a map;
      # every other negative and every numerical comparison is retained.
      expect_false(is.null(res$structural_certificate), info = nm)
      expect_true(res$reduction_summary$reduction_active, info = nm)
      expect_true(res$reduction_summary$n_links_used < res$reduction_summary$n_links_all, info = nm)
    } else {
      expect_null(res$structural_certificate, info = nm)
      expect_false(res$reduction_summary$reduction_active, info = nm)
      expect_identical(res$reduction_summary$n_links_used, res$reduction_summary$n_links_all, info = nm)
    }
  }

  # Direct cp04_complete_reduction negative and bound checks.
  pos_res <- od_run_pipeline(d0, est = "dr", chunk = 7L, reduce = TRUE)
  b_dr <- pos_res$nuisance_registry$cp04
  pl_all <- pos_res$point_links_all
  des_all <- pos_res$design_all
  dli <- match(seq_len(10L), des_all$links$catalog_row)
  fo <- pos_res$frame_owner
  cb <- pos_res$controls$pipeline_control$class_bounds
  oc <- pos_res$controls$pipeline_control$owner_control

  expect_null(api$cp04_complete_reduction(NULL, pl_all, des_all, dli, fo, cb, oc))
  bad_pl <- pl_all; bad_pl[[2L]]$link_id <- "wrong_id"
  expect_null(api$cp04_complete_reduction(b_dr, bad_pl, des_all, dli, fo, cb, oc))

  cb_links <- cb; cb_links$max_links <- 5L
  expect_error(api$cp04_complete_reduction(b_dr, pl_all, des_all, dli, fo, cb_links, oc), class = "cdid_ipw_bound_error")
  cb_dates <- cb; cb_dates$max_class_dates <- 3L
  expect_error(api$cp04_complete_reduction(b_dr, pl_all, des_all, dli, fo, cb_dates, oc), class = "cdid_ipw_bound_error")
  cb_e <- cb; cb_e$max_E_entries <- 2L
  expect_error(api$cp04_complete_reduction(b_dr, pl_all, des_all, dli, fo, cb_e, oc), class = "cdid_ipw_bound_error")
  oc_work <- oc; oc_work$max_work_ops <- b_dr$accounting$work_ops_bound + 10
  expect_error(api$cp04_complete_reduction(b_dr, pl_all, des_all, dli, fo, cb, oc_work), class = "cdid_ipw_owner_bound_error")

  # Non-dyadic complete-support rates: guarded rates are exactly one, while
  # unchecked diagnostics preserve the independently summed raw quotient.
  u0 <- d0[!duplicated(d0$id), ]
  u0 <- u0[order(u0$id), ]
  tg0 <- u0$G == 3
  csum <- function(wv, mask) {
    s <- 0
    for (a0 in seq.int(1L, 41L, by = 7L)) {
      idx <- seq.int(a0, min(41L, a0 + 6L))
      h <- idx[mask[idx]]
      if (length(h)) s <- s + sum(wv[h])
    }
    s
  }
  w_nd <- NULL
  for (k in seq_len(200L)) {
    cand <- 1 + ((seq_len(41L) * k) %% 11L) / 10 + 1 / (2 * (k %% 5L) + 3)
    if (csum(cand, tg0) > sum(cand[tg0]) || csum(cand, !tg0) > sum(cand[!tg0])) {
      w_nd <- cand
      break
    }
  }
  expect_false(is.null(w_nd))
  d_nd <- d0
  d_nd$w <- w_nd[match(d_nd$id, u0$id)]
  for (est in c("reg", "dr")) {
    e_r <- od_run_pipeline(d_nd, est = est, chunk = 7L, weightsname = "w", reduce = TRUE)
    e_u <- od_run_pipeline(d_nd, est = est, chunk = 7L, weightsname = "w", reduce = FALSE)
    expect_s3_class(e_r, "cdid_ipw_pipeline_result")
    expect_s3_class(e_u, "cdid_ipw_pipeline_result")
    expect_null(e_r$structural_certificate)
    expect_true(all(vapply(e_r$nuisance_registry$cp04$supports, `[[`, numeric(1), "rate") == 1))
    expect_equal(e_r$reporting$estimates, e_u$reporting$estimates, tolerance = 1e-13)
    expect_equal(e_r$reporting$se, e_u$reporting$se, tolerance = 1e-13)
    # The unchecked reader keeps raw chunked division; it does not apply the
    # guarded rate refusal. Full support is independently proved by unit counts.
    unchecked <- api$construct_ipw_link_view(e_r$frame_owner, 1, 2, g = 3, K = Inf,
      control = e_r$controls$pipeline_control$owner_control, require_support = FALSE)
    checked <- api$construct_ipw_link_view(e_r$frame_owner, 1, 2, g = 3, K = Inf,
      control = e_r$controls$pipeline_control$owner_control, require_support = TRUE)
    expect_identical(e_r$frame_owner$a_id, u0$id)
    wa <- e_r$frame_owner$w_a
    pair_mass <- c(csum(wa, tg0), csum(wa, !tg0))
    cohort_mass <- c(sum(wa[tg0]), sum(wa[!tg0]))
    cohort_count <- as.integer(c(sum(tg0), sum(!tg0)))
    expect_identical(unname(unchecked$support$by_cohort_pair_count), cohort_count)
    expect_identical(unname(unchecked$support$by_cohort_pair_mass), pair_mass)
    expect_identical(unname(unchecked$support$rates), pair_mass / cohort_mass)
    expect_identical(unname(checked$support$by_cohort_pair_count), cohort_count)
    expect_identical(unname(checked$support$by_cohort_pair_mass), pair_mass)
    expect_identical(unname(checked$support$rates), c(1, 1))
  }
})


# Native review corrections: the construction cap must refuse before design,
# forest or path work; exact issued support facts establish eligibility.
test_that("stage 6 full map work is composed before construction and structural admission is exact", {
  d <- od_nonadd_panel(41L, 5L, 18L, seed = 31L)
  for (est in c("reg", "dr")) {
    res <- od_run_pipeline(d, est = est, chunk = 7L, reduce = TRUE)
    backend <- res$nuisance_registry$cp04
    points <- res$point_links_all
    design <- res$design_all
    M <- length(points)
    dli <- match(seq_len(M), design$links$catalog_row)
    cb <- res$controls$pipeline_control$class_bounds
    oc <- res$controls$pipeline_control$owner_control
    e <- new.env(parent = api)
    fn <- api$cp04_complete_reduction; environment(fn) <- e
    calls <- c(design = 0L, forest = 0L, path = 0L)
    e$validate_ipw_design_for_classes <- function(...) {
      calls["design"] <<- calls["design"] + 1L
      api$validate_ipw_design_for_classes(...)
    }
    e$incidence_forest <- function(...) {
      calls["forest"] <<- calls["forest"] + 1L
      api$incidence_forest(...)
    }
    e$ipw_forest_path <- function(...) {
      calls["path"] <<- calls["path"] + 1L
      api$ipw_forest_path(...)
    }
    # Independently count the issued-date path/sort upper bound.
    D <- min(length(backend$complete_support$dates), 2 * M)
    W <- 256 * M * (D + ceiling(log2(M + 1))) + 64 * M * min(M, max(0, D - 1))
    B <- backend$accounting$work_ops_bound
    invoke <- function(control, external = 0, b = backend, pl = points) {
      fn(b, pl, design, dli, res$frame_owner, cb, control,
        external_work_ops = external)
    }
    cap <- oc; cap$max_work_ops <- B + W - 1
    expect_error(invoke(cap), class = "cdid_ipw_owner_bound_error")
    expect_identical(unname(calls), rep(0L, 3L))
    cap$max_work_ops <- B + W
    expect_error(invoke(cap, external = B + 1), class = "cdid_ipw_owner_bound_error")
    expect_identical(unname(calls), rep(0L, 3L))
    cert <- invoke(cap)
    expect_identical(cert$E, res$structural_certificate$E)
    expect_gt(calls[["forest"]], 0L)
    expect_gt(calls[["path"]], 0L)

    # A deliberate removed-early-check mutant reaches construction at the same
    # one-below cap: the zero-call assertion above can distinguish this bug.
    calls[] <- 0L
    e$ipw_owner_safe_sum <- function(vals, max_bound = NULL, label = "sum") {
      if (identical(label, "pipeline/backend plus CP04 complete-support reduction work before construction"))
        api$ipw_owner_safe_sum(vals, label = label) else
          api$ipw_owner_safe_sum(vals, max_bound, label)
    }
    cap$max_work_ops <- B + W - 1
    mutant <- invoke(cap)
    expect_gt(calls[["forest"]], 0L)
    expect_gt(calls[["path"]], 0L)
    rm("ipw_owner_safe_sum", envir = e)

    # No numerical-near-equality predicate licenses a map. The copied helper's
    # abs trap would catch the removed mass-tolerance check on valid input.
    e$abs <- function(...) stop("numerical eligibility comparison reached")
    expect_identical(invoke(oc)$E, cert$E)
    tiny_pair <- points
    tiny_pair[[1L]]$pair$s <- tiny_pair[[1L]]$pair$s + .Machine$double.eps
    expect_null(invoke(oc, pl = tiny_pair))
    changed_rate <- backend
    sid <- changed_rate$links[[points[[1L]]$cp04_index]]$support_id
    changed_rate$supports[[sid]]$rate <- 1 - .Machine$double.eps
    changed_rate <- api$ipw_owner_seal(changed_rate,
      if (est == "dr") "dr_baseline_backend" else "or_baseline_backend")
    expect_null(invoke(oc, b = changed_rate))
  }
})

# ---------------------------------------------------------------------------
# Common-design OR fit compression and bounded complete-support batching
# (p = 0, 1, 2, repeated L, B = 1, case-weight derivatives, weights, chunks,
# +1e8 offsets, aliases/constants, rank refusals and wrong-formula controls)
# ---------------------------------------------------------------------------

od_common_design <- function(N = 41L, p = 2L, seed = 19L) {
  set.seed(seed)
  if (p == 0L) return(list(X = matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)")),
    recipe = od_icpt, trt_cols = "(Intercept)"))
  x1 <- sin(seq_len(N) * 0.37) + 0.25 * ((seq_len(N) %% 3L) - 1)
  if (p == 1L) {
    X <- cbind(`(Intercept)` = 1, x1 = x1)
    rec <- api$or_pooled_outcome_recipe(c("(Intercept)", "x1"), TRUE, "x1")
    return(list(X = X, recipe = rec, trt_cols = c("(Intercept)", "x1")))
  }
  x2 <- cos(seq_len(N) * 0.53) - 0.15 * ((seq_len(N) %% 5L) - 2)
  X <- cbind(`(Intercept)` = 1, x1 = x1, x2 = x2)
  rec <- api$or_pooled_outcome_recipe(c("(Intercept)", "x1", "x2"), TRUE, c("x1", "x2"))
  list(X = X, recipe = rec, trt_cols = c("(Intercept)", "x1", "x2"))
}

# Independent raw-row pooled WLS reference for any common-design fixture.
od_independent_pooled_wls <- function(fx, recipe) {
  st <- fx$reader$state; N <- st$N
  fo <- get("frame_owner", fx$reader$frame_ref, inherits = FALSE)
  Z_f <- fo$Z_a[st$a_index, , drop = FALSE]
  z_f <- if (recipe$p > 0L) Z_f[, recipe$slope_index, drop = FALSE] else matrix(0, N, 0L)
  y_at <- function(id, t) fx$d$y[match(paste(id, t), paste(fx$d$id, fx$d$t))]
  nv <- is.infinite(st$g); tg <- st$g == 3
  w_c <- st$w[nv]; z_c <- z_f[nv, , drop = FALSE]
  # Unique control blocks in link appearance order (keyed by s_t).
  pair_keys <- vapply(fx$pairs, function(p) paste0(p[1L], "_", p[2L]), character(1L))
  u_idx <- which(!duplicated(pair_keys))
  # The independent raw-row reference uses canonical numeric pair order, as
  # required by the retained block schema; expected values still use raw rows.
  u_idx <- u_idx[order(vapply(fx$pairs[u_idx], function(q) as.double(q[1L]), numeric(1)),
    vapply(fx$pairs[u_idx], function(q) as.double(q[2L]), numeric(1)))]
  u_pairs <- fx$pairs[u_idx]
  G <- length(u_pairs)
  L_vec <- vapply(u_pairs, function(p) as.double(p[2L] - p[1L]), numeric(1L))
  v_vec <- rep(1, G) # complete-support rate = 1 => v_b = 1
  D <- sum(v_vec * L_vec^2)
  m_c <- sum(w_c)
  zbar <- if (recipe$p > 0L) colSums(w_c * z_c) / m_c else numeric(0L)
  ztilde <- if (recipe$p > 0L) sweep(z_c, 2L, zbar, "-") else matrix(0, sum(nv), 0L)
  dy_c_mat <- vapply(u_pairs, function(p) y_at(st$id[nv], p[2L]) - y_at(st$id[nv], p[1L]), numeric(sum(nv)))
  if (!is.matrix(dy_c_mat)) dy_c_mat <- matrix(dy_c_mat, ncol = G)
  ybar_vec <- colSums(w_c * dy_c_mat) / m_c
  # Collapsed response ystar_i = sum_b (v_b L_b / D) dy_{i,b}
  ystar <- as.vector(dy_c_mat %*% (v_vec * L_vec / D))
  ystar_bar <- sum(w_c * ystar) / m_c
  ystar_tilde <- ystar - ystar_bar
  Gram <- if (recipe$p > 0L) D * crossprod(sqrt(w_c) * ztilde) else matrix(0, 0L, 0L)
  rhs <- if (recipe$p > 0L) D * as.vector(crossprod(w_c * ztilde, ystar_tilde)) else numeric(0L)
  beta <- if (recipe$p > 0L) as.vector(solve(Gram, rhs)) else numeric(0L)
  alpha <- ybar_vec - L_vec * as.vector(zbar %*% beta)
  Q_star <- sum(w_c * (ystar_tilde - as.vector(ztilde %*% beta))^2)
  e_mat <- dy_c_mat - outer(ystar, L_vec)
  e_bar <- colSums(w_c * e_mat) / m_c
  K_val <- sum(v_vec * colSums(w_c * sweep(e_mat, 2L, e_bar, "-")^2))
  # Direct stacked block-by-block WLS sum of squares:
  rss_direct <- 0
  for (b in seq_len(G)) {
    r_b <- dy_c_mat[, b] - alpha[b] - L_vec[b] * as.vector(z_c %*% beta)
    rss_direct <- rss_direct + v_vec[b] * sum(w_c * r_b^2)
  }
  # Full-F OR IF for each link:
  m_t <- sum(st$w[tg])
  zbar_t <- if (recipe$p > 0L) colSums(st$w[tg] * z_f[tg, , drop = FALSE]) / m_t else numeric(0L)
  # Global slope IF xi_beta (N x p):
  xi_beta <- matrix(0, N, recipe$p)
  if (recipe$p > 0L) {
    score_c <- (w_c * ztilde) * (D * (ystar_tilde - as.vector(ztilde %*% beta)))
    xi_c <- score_c %*% solve(Gram / N)
    xi_beta[nv, ] <- xi_c
  }
  phi_or <- matrix(0, N, length(fx$pairs))
  mu_or <- numeric(length(fx$pairs))
  for (j in seq_along(fx$pairs)) {
    p <- fx$pairs[[j]]
    b <- match(paste0(p[1L], "_", p[2L]), pair_keys[u_idx])
    dy <- y_at(st$id, p[2L]) - y_at(st$id, p[1L])
    L_b <- L_vec[b]
    pred <- alpha[b] + L_b * as.vector(z_f %*% beta)
    res <- dy - pred
    mu_j <- sum(st$w[tg] * res[tg]) / m_t
    mu_or[j] <- mu_j
    xi_alpha <- numeric(N)
    xi_alpha[nv] <- (st$w[nv] * res[nv]) / (m_c / N)
    if (recipe$p > 0L) xi_alpha <- xi_alpha - as.vector(xi_beta %*% (L_b * zbar))
    direct <- numeric(N)
    direct[tg] <- (st$w[tg] * (res[tg] - mu_j)) / (m_t / N)
    c_j <- if (recipe$p > 0L) L_b * zbar_t else numeric(0L)
    slope_term <- if (recipe$p > 0L) as.vector(xi_beta %*% c_j) else 0
    phi_or[, j] <- direct - xi_alpha - slope_term
  }
  list(alpha = alpha, beta = beta, D = D, Gram = Gram, K = K_val, Q_star = Q_star,
    rss = K_val + D * Q_star, rss_direct = rss_direct, mu = mu_or, phi = phi_or, xi_beta = xi_beta)
}

test_that("common-design p = 0, 1, 2, repeated L and B = 1 match generic route and independent pooled WLS references", {
  d <- od_panel(N = 41L, T = 5L, seed = 23L)
  pair_sets <- list(
    all_pairs = od_pairs(5L),
    repeated_L = list(c(1, 3), c(2, 4), c(3, 5), c(1, 2), c(2, 3), c(1, 3)),
    single_block = list(c(2, 4))
  )
  for (p in c(0L, 1L, 2L)) {
    des <- od_common_design(41L, p = p, seed = 19L + p)
    for (pname in names(pair_sets)) {
      prs <- pair_sets[[pname]]
      fx <- od_fixture(d = d, chunk = 7L, perm = (p == 2L), X = des$X, V = des$X, pairs = prs)
      ref <- od_independent_pooled_wls(fx, des$recipe)
      expect_equal(ref$rss, ref$rss_direct, tolerance = 1e-12)

      for (method in c("or", "dr")) {
        fit_c <- if (method == "or") od_or(fx, recipe = des$recipe) else
          od_dr(fx, recipe = des$recipe, treatment_columns = des$trt_cols)
        fit_g <- if (method == "or") od_or(fx, od_general, recipe = des$recipe) else
          od_dr(fx, od_general, recipe = des$recipe, treatment_columns = des$trt_cols)

        info_tag <- paste(method, "p=", p, pname)
        expect_false(is.null(fit_c$complete_support), info = info_tag)
        expect_identical(fit_c$shared$p, p, info = info_tag)
        expect_identical(fit_c$complete_support$n_z, p + 1L, info = info_tag)
        expect_null(fit_g$complete_support, info = info_tag)
        # D75: no other case, method, field or comparison uses the mixed predicate.
        approved_M_dr <- list(
          "dr:1:repeated_L" = list(index = 3L, link_id = "L3_5"),
          "dr:1:single_block" = list(index = 1L, link_id = "L2_4"),
          "dr:2:all_pairs" = list(index = 10L, link_id = "L4_5")
        )[[paste(method, p, pname, sep = ":")]]
        mixed_scalar_path <- NULL
        if (!is.null(approved_M_dr)) {
          j <- approved_M_dr$index
          expect_identical(fit_c$links[[j]]$link_id, approved_M_dr$link_id, info = info_tag)
          expect_identical(fit_g$links[[j]]$link_id, approved_M_dr$link_id, info = info_tag)
          mixed_scalar_path <- paste0("root$fields$links$", j, "$M_dr")
        }
        od_expect_route_match(od_backend_view(fit_c), od_backend_view(fit_g), tol = 1e-12, info = info_tag,
          mixed_scalar_path = mixed_scalar_path)

        # Verify alpha, beta, Q = K + D * Q_star, objective and D-scaled Gram against independent reference.
        alpha_hat <- vapply(fit_c$regressions, `[[`, numeric(1L), "alpha")
        expect_equal(unname(alpha_hat), unname(ref$alpha), tolerance = 1e-12, info = info_tag)
        expect_equal(unname(fit_c$shared$beta), unname(ref$beta), tolerance = 1e-12, info = info_tag)
        expect_equal(fit_c$shared$rss, ref$rss, tolerance = 1e-12, info = info_tag)
        expect_equal(fit_c$shared$objective, ref$rss / 41, tolerance = 1e-12, info = info_tag)
        if (p > 0L) {
          R_unpiv <- fit_c$shared$fit$R[, order(fit_c$shared$fit$pivot), drop = FALSE]
          expect_equal(crossprod(R_unpiv), unname(ref$Gram), tolerance = 1e-12, info = info_tag)
          fac <- fit_c$shared$summary$residual$factor
          expect_equal(crossprod(fac[, seq_len(p), drop = FALSE]), unname(ref$Gram), tolerance = 1e-12, info = info_tag)
        }

        # Replay full-F IF and stream covariance across both chunk/link batch sizes.
        replay <- if (method == "or") api$replay_or_baseline_if else api$replay_dr_baseline_if
        for (details in c(FALSE, TRUE)) {
          bc <- od_batches(fit_c, replay, 7L, min(3L, length(prs)), "frame", details)
          bg <- od_batches(fit_g, replay, 7L, min(3L, length(prs)), "frame", details)
          if (method == "dr" && !details) {
            drop <- function(b) { b$batches <- lapply(b$batches, function(x) {
              x$diagnostics$n_regression_score_evaluations <- NULL; x }); b }
            bc <- drop(bc); bg <- drop(bg)
          }
          od_expect_route_match(bc$batches, bg$batches, tol = 1e-12, info = paste(info_tag, details))
          if (method == "or" && !details) {
                # Reassemble link blocks into full L matrix:
            phi_full <- matrix(0, 41L, length(prs))
            for (b in bc$batches) phi_full[b$rows, b$links] <- b$phi
            expect_equal(unname(phi_full), unname(ref$phi), tolerance = 1e-12, info = info_tag)
            expect_equal(vapply(fit_c$links, `[[`, numeric(1L), "mu"), ref$mu, tolerance = 1e-12, info = info_tag)
          }
        }
        cov_fn <- if (method == "or") function(f) api$cp04_stream_covariance(f, seq_along(f$links)) else api$stream_dr_covariance
        cc <- cov_fn(fit_c); cg <- cov_fn(fit_g)
        for (nm in c("estimates", "covariance", "mean_phi", "chunk_rows"))
          od_expect_route_match(cc[[nm]], cg[[nm]], tol = 1e-12, info = paste(info_tag, nm))
      }

      # Dense DR table helper check for p = 0, 1, 2:
      st <- fx$reader$state; N <- st$N
      y_at <- function(id, t) fx$d$y[match(paste(id, t), paste(fx$d$id, fx$d$t))]
      Y <- sapply(seq_len(5L), function(t) y_at(st$id, t))
      fo <- get("frame_owner", fx$reader$frame_ref, inherits = FALSE)
      Z_f <- fo$Z_a[st$a_index, , drop = FALSE]
      links_tb <- lapply(seq_along(prs), function(k) {
        pr <- prs[[k]]
        list(s = pr[1L], t = pr[2L], g = 3, K = Inf, link_id = paste0("L", k, "_", pr[1L], "_", pr[2L]))
      })
      tb <- api$fit_dr_baseline_table(X_or = Z_f, V = Z_f, Y = Y, obs = matrix(1, N, 5L), g = st$g, w = st$w,
        calendar = 1:5, links = links_tb, outcome_recipe = des$recipe, glm_control = od_glm, rank_tol = 1e-10)
      fit_dr <- od_dr(fx, recipe = des$recipe, treatment_columns = des$trt_cols)
      dr_det <- od_batches(fit_dr, api$replay_dr_baseline_if, min(16L, fit_dr$control$max_chunk_rows), length(prs), details = TRUE)$batches
      col_dr <- function(nm) do.call(rbind, lapply(dr_det, function(b) if (nm == "phi") b$phi else b$details[[nm]]))
      expect_equal(vapply(fit_dr$links, `[[`, numeric(1L), "mu_dr"), unname(tb$estimates), tolerance = 1e-12)
      expect_equal(vapply(fit_dr$links, `[[`, numeric(1L), "mu_ipw"), unname(tb$estimates_ipw), tolerance = 1e-12)
      for (nm in c("phi", "psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw"))
        expect_equal(unname(col_dr(nm)), unname(tb[[nm]]), tolerance = 1e-12, info = paste("tb", p, pname, nm))
    }
  }

  # Full pipeline check (OR and DR with xformla = ~x and ~x + x2): public estimates, SE, IF, vcov, aggregations and diagnostics.
  dp <- od_nonadd_panel(41L, 5L, 18L, seed = 41L, id_perm = TRUE, rows_perm = TRUE)
  dp$x2 <- cos(dp$id * 0.31)
  for (xf in list(~x, ~x + x2)) {
    for (est in c("reg", "dr")) {
      res_c <- od_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, bstrap = TRUE, cband = TRUE, biters = 49L, seed = 17L)
      res_g <- od_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, bstrap = TRUE, cband = TRUE, biters = 49L, seed = 17L, ns = local({
        e <- new.env(parent = api)
        for (nm in c("run_chained_ipw_pipeline", "build_cp04_public_registry", "or_public_fit", "dr_public_fit", "fit_or_baseline", "fit_dr_baseline")) {
          f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
        }
        e$or_baseline_complete_facts <- function(...) NULL
        e
      }))
      expect_false(is.null(res_c$nuisance_registry$cp04$complete_support))
      expect_null(res_g$nuisance_registry$cp04$complete_support)
      # The shared-X extension now admits this same complete fitted object;
      # all original public/IF/covariance/aggregation checks below remain.
      expect_false(is.null(res_c$structural_certificate))
      av <- res_c$reporting$available
      expect_true(length(av) > 0L && any(av))
      expect_identical(res_c$reporting$target_meta, res_g$reporting$target_meta)
      expect_identical(av, res_g$reporting$available)
      phi_saved <- res_c$reporting$if_full[, av, drop = FALSE]
      expect_identical(nrow(phi_saved), 41L)
      expect_true(ncol(phi_saved) > 0L && all(is.finite(phi_saved)))
      original_ids <- get("cp04_reader", res_c$nuisance_registry$cp04$reader_ref, inherits = FALSE)$state$id
      input_units <- dp[!duplicated(dp$id), c("id", "cl")]
      original_clusters <- input_units$cl[match(original_ids, input_units$id)]
      expect_true(length(original_ids) == nrow(phi_saved) && !anyDuplicated(original_ids) && !anyNA(original_clusters))
      cluster_sums <- rowsum(phi_saved, original_clusters, reorder = FALSE)
      centered <- sweep(cluster_sums, 2L, colMeans(cluster_sums), "-")
      cov_saved <- crossprod(centered) / nrow(phi_saved)^2
      expect_equal(unname(cov_saved), unname(res_c$reporting$vcov[av, av, drop = FALSE]), tolerance = 1e-12)
      expect_equal(unname(sqrt(diag(cov_saved))), unname(res_c$reporting$analytical_se[av]), tolerance = 1e-12)
      expect_equal(res_c$reporting$estimates, res_g$reporting$estimates, tolerance = 1e-12)
      expect_equal(res_c$reporting$se, res_g$reporting$se, tolerance = 1e-12)
      expect_equal(res_c$reporting$analytical_se, res_g$reporting$analytical_se, tolerance = 1e-12)
      expect_equal(res_c$reporting$if_full, res_g$reporting$if_full, tolerance = 1e-12)
      expect_equal(res_c$reporting$vcov, res_g$reporting$vcov, tolerance = 1e-12)
      expect_equal(res_c$reporting$cband_lower, res_g$reporting$cband_lower, tolerance = 1e-12)
      expect_equal(res_c$reporting$cband_upper, res_g$reporting$cband_upper, tolerance = 1e-12)
      for (ag in c("simple", "group", "dynamic", "calendar")) {
        expect_equal(res_c$aggregations[[ag]]$overall$estimate, res_g$aggregations[[ag]]$overall$estimate, tolerance = 1e-12)
        expect_equal(res_c$aggregations[[ag]]$overall$se, res_g$aggregations[[ag]]$overall$se, tolerance = 1e-12)
        expect_equal(res_c$aggregations[[ag]]$estimates, res_g$aggregations[[ag]]$estimates, tolerance = 1e-12)
        expect_equal(res_c$aggregations[[ag]]$se, res_g$aggregations[[ag]]$se, tolerance = 1e-12)
      }
      diag_c <- api$cp04_public_diagnostics(res_c, res_c$controls$pipeline_control)
      diag_g <- api$cp04_public_diagnostics(res_g, res_g$controls$pipeline_control)
      expect_equal(diag_c$coefficients, diag_g$coefficients, tolerance = 1e-12)
      expect_equal(diag_c$fitted_catalog$mu, diag_g$fitted_catalog$mu, tolerance = 1e-12)
    }
  }
})

test_that("common-design case-weight Gateaux derivatives including treatment PS and global slope forcing match finite differences", {
  N <- 25L; T <- 4L
  d <- od_panel(N = N, T = T, seed = 61L)
  des <- od_common_design(N = N, p = 2L, seed = 62L)
  prs <- list(c(1, 2), c(1, 3), c(2, 4), c(1, 4))
  fx0 <- od_fixture(d = d, chunk = N, X = des$X, V = des$X, pairs = prs)
  st <- fx0$reader$state
  y_at <- function(id, t) d$y[match(paste(id, t), paste(d$id, d$t))]
  Y <- sapply(seq_len(T), function(t) y_at(st$id, t))
  fo <- get("frame_owner", fx0$reader$frame_ref, inherits = FALSE)
  Z_f <- fo$Z_a[st$a_index, , drop = FALSE]
  links_tb <- lapply(seq_along(prs), function(k) {
    pr <- prs[[k]]
    list(s = pr[1L], t = pr[2L], g = 3, K = Inf, link_id = paste0("L", k))
  })

  fit_or <- od_or(fx0, recipe = des$recipe)
  fit_dr <- od_dr(fx0, recipe = des$recipe, treatment_columns = des$trt_cols)
  phi_or <- api$evaluate_or_baseline_if_chunk(fit_or, seq_len(N), seq_along(prs), "frame", FALSE)$phi
  phi_dr <- api$evaluate_dr_baseline_if_chunk(fit_dr, seq_len(N), seq_along(prs), "frame", FALSE)$phi
  det_dr <- api$evaluate_dr_baseline_if_chunk(fit_dr, seq_len(N), seq_along(prs), "frame", TRUE)$details
  expect_equal(unname(phi_dr), unname(det_dr$psi_dir + det_dr$psi_or + det_dr$psi_ps + det_dr$psi_sel), tolerance = 1e-13)
  # Both psi_or and psi_ps are non-trivial when p_or = 2 and p_trt = 3:
  expect_gt(max(abs(det_dr$psi_or)), 1e-2)
  expect_gt(max(abs(det_dr$psi_ps)), 1e-2)

  # Finite-difference case-weight Gateaux derivatives:
  eps <- 1e-6
  fd_or <- matrix(0, N, length(prs))
  fd_dr <- matrix(0, N, length(prs))
  fd_beta <- matrix(0, N, des$recipe$p)
  for (i in seq_len(N)) {
    w_p <- rep(1, N); w_p[i] <- 1 + eps * N
    w_m <- rep(1, N); w_m[i] <- 1 - eps * N
    tb_p <- api$fit_dr_baseline_table(X_or = Z_f, V = Z_f, Y = Y, obs = matrix(1, N, T), g = st$g, w = w_p,
      calendar = seq_len(T), links = links_tb, outcome_recipe = des$recipe, glm_control = od_glm, rank_tol = 1e-10)
    tb_m <- api$fit_dr_baseline_table(X_or = Z_f, V = Z_f, Y = Y, obs = matrix(1, N, T), g = st$g, w = w_m,
      calendar = seq_len(T), links = links_tb, outcome_recipe = des$recipe, glm_control = od_glm, rank_tol = 1e-10)
    fd_dr[i, ] <- (tb_p$estimates - tb_m$estimates) / (2 * eps)
    fd_or[i, ] <- (tb_p$estimates_or - tb_m$estimates_or) / (2 * eps)
    fd_beta[i, ] <- (tb_p$shared$beta - tb_m$shared$beta) / (2 * eps)
  }
  expect_lt(max(abs(unname(phi_or) - fd_or)), 1e-7)
  expect_lt(max(abs(unname(phi_dr) - fd_dr)), 1e-7)
  ref_wls <- od_independent_pooled_wls(fx0, des$recipe)
  expect_lt(max(abs(ref_wls$xi_beta - fd_beta)), 1e-7)
})

test_that("common-design admissible weights, permutations, partial chunks, +1e8 offsets, aliases, constants, rank refusals and wrong-formula controls", {
  d0 <- od_panel(N = 41L, T = 5L, seed = 71L)
  des2 <- od_common_design(41L, p = 2L, seed = 72L)
  expect_measured_m_dr_structure <- function(x, y) {
    expect_true(typeof(x) == "double" && typeof(y) == "double" &&
      length(x) == 1L && length(y) == 1L &&
      identical(attributes(x), attributes(y), num.eq = FALSE) &&
      all(is.finite(x)) && all(is.finite(y)) &&
      identical(x == 0, y == 0) && identical(sign(x), sign(y)) &&
      identical(1 / x[x == 0], 1 / y[y == 0]))
  }

  # 1. Admissible weights, permutations, partial chunks (1, 7, 16, 41) and +1e8 outcome offsets:
  d_off <- d0; d_off$y <- d_off$y + 1e8
  fx_base <- od_fixture(d = d0, chunk = 41L, X = des2$X, V = des2$X)
  or_base <- od_or(fx_base, recipe = des2$recipe)
  dr_base <- od_dr(fx_base, recipe = des2$recipe, treatment_columns = des2$trt_cols)

  for (cfg in list(
    list(d = d_off, chunk = 7L, perm = FALSE, wt = 1),
    list(d = d0, chunk = 1L, perm = FALSE, wt = 0.375),
    list(d = d0, chunk = 16L, perm = TRUE, wt = 2.5),
    list(d = d_off, chunk = 41L, perm = TRUE, wt = 1)
  )) {
    fx <- od_fixture(d = cfg$d, chunk = cfg$chunk, perm = cfg$perm, w = rep(cfg$wt, 41L), X = des2$X, V = des2$X)
    for (method in c("or", "dr")) {
      fc <- if (method == "or") od_or(fx, recipe = des2$recipe) else od_dr(fx, recipe = des2$recipe, treatment_columns = des2$trt_cols)
      fg <- if (method == "or") od_or(fx, od_general, recipe = des2$recipe) else od_dr(fx, od_general, recipe = des2$recipe, treatment_columns = des2$trt_cols)
      expect_false(is.null(fc$complete_support))
      if (method == "dr" && identical(cfg$chunk, 1L) && !cfg$perm && identical(cfg$wt, 0.375))
        expect_measured_m_dr_structure(fc$links[[4L]]$M_dr, fg$links[[4L]]$M_dr)
      od_expect_route_match(od_backend_view(fc), od_backend_view(fg), tol = 1e-12,
        mixed_scalar_path = if (method == "dr" && identical(cfg$chunk, 1L) && !cfg$perm && identical(cfg$wt, 0.375))
          "root$fields$links$4$M_dr" else NULL)
      if (!cfg$perm && cfg$wt == 1) {
        # +1e8 offset invariance: raw endpoint subtraction before summation preserves alpha, beta and mu!
        ref_fit <- if (method == "or") or_base else dr_base
        expect_equal(fc$shared$beta, ref_fit$shared$beta, tolerance = 1e-8)
        expect_equal(vapply(fc$regressions, `[[`, numeric(1L), "alpha"),
          vapply(ref_fit$regressions, `[[`, numeric(1L), "alpha"), tolerance = 1e-8)
        expect_equal(vapply(fc$links, `[[`, numeric(1L), "mu"),
          vapply(ref_fit$links, `[[`, numeric(1L), "mu"), tolerance = 1e-8)
      }
    }
  }

  # 2. Exact aliases and constant slope columns trigger conditional original-row alias verification:
  X_alias <- cbind(des2$X, x_const = 5, x_alias = 2 * des2$X[, "x1"] - 0.5 * des2$X[, "x2"] + 3)
  rec_alias <- api$or_pooled_outcome_recipe(colnames(X_alias), TRUE, c("x1", "x2", "x_const", "x_alias"))
  fx_alias <- od_fixture(d = d0, chunk = 7L, X = X_alias, V = des2$X)
  for (method in c("or", "dr")) {
    fc <- if (method == "or") od_or(fx_alias, recipe = rec_alias) else od_dr(fx_alias, recipe = rec_alias, treatment_columns = des2$trt_cols)
    fg <- if (method == "or") od_or(fx_alias, od_general, recipe = rec_alias) else od_dr(fx_alias, od_general, recipe = rec_alias, treatment_columns = des2$trt_cols)
    expect_false(is.null(fc$complete_support))
    expect_true(fc$shared$verification_scan)
    expect_gt(fc$accounting$complete_support$n_alias_rescans, 0)
    expect_identical(fc$shared$keep, c(1L, 2L))
    expect_identical(fc$shared$basis$aliased_columns, c("x_const", "x_alias"))
    if (method == "dr") expect_measured_m_dr_structure(fc$links[[4L]]$M_dr, fg$links[[4L]]$M_dr)
    od_expect_route_match(od_backend_view(fc), od_backend_view(fg), tol = 1e-12,
      mixed_scalar_path = if (method == "dr") "root$fields$links$4$M_dr" else NULL)
  }

  # 3. Control-only null direction that does NOT reconstruct treated rows refuses identically on both routes:
  u0 <- d0[!duplicated(d0$id), c("id", "g")]; u0 <- u0[order(u0$id), ]
  x_trt_only <- ifelse(u0$g == 3, seq_len(41L), 0)
  X_bad_null <- cbind(des2$X, x_trt = x_trt_only)
  rec_bad_null <- api$or_pooled_outcome_recipe(colnames(X_bad_null), TRUE, c("x1", "x2", "x_trt"))
  fx_bad_null <- od_fixture(d = d0, chunk = 7L, X = X_bad_null, V = des2$X)
  for (method in c("or", "dr")) {
    ec <- tryCatch(if (method == "or") od_or(fx_bad_null, recipe = rec_bad_null) else
      od_dr(fx_bad_null, recipe = rec_bad_null, treatment_columns = des2$trt_cols), error = identity)
    eg <- tryCatch(if (method == "or") od_or(fx_bad_null, od_general, recipe = rec_bad_null) else
      od_dr(fx_bad_null, od_general, recipe = rec_bad_null, treatment_columns = des2$trt_cols), error = identity)
    expect_s3_class(ec, "cdid_or_wls_rank_error")
    expect_identical(class(ec), class(eg))
    expect_identical(conditionMessage(ec), conditionMessage(eg))
  }

  # 4. Deliberate wrong-formula negative controls MUST fail:
  ref <- od_independent_pooled_wls(fx_base, des2$recipe)
  phi_or <- api$evaluate_or_baseline_if_chunk(or_base, seq_len(41L), seq_along(fx_base$views), "frame", FALSE)$phi
  det_dr <- api$evaluate_dr_baseline_if_chunk(dr_base, seq_len(41L), seq_along(fx_base$views), "frame", TRUE)$details
  phi_dr <- api$evaluate_dr_baseline_if_chunk(dr_base, seq_len(41L), seq_along(fx_base$views), "frame", FALSE)$phi

  # (a) Omitted global slope IF control: dropping psi_or in DR or dropping non-local blocks from beta_gids in OR:
  phi_dr_no_slope <- det_dr$psi_dir + det_dr$psi_ps + det_dr$psi_sel
  expect_false(isTRUE(all.equal(phi_dr, phi_dr_no_slope, tolerance = 1e-12)))
  plan_or <- api$or_baseline_native_plan(or_base, seq_along(or_base$links), 41L)
  plan_local_only <- plan_or; plan_local_only$ref <- NULL
  plan_local_only$dimnames[[1L]][[6L]] <- plan_local_only$dimnames[[1L]][[6L]][1L]
  plan_local_only$dimnames[[1L]][[7L]] <- plan_local_only$dimnames[[1L]][[7L]][, 1L, drop = FALSE]
  plan_local_only <- api$ipw_owner_seal(plan_local_only, "or_baseline_native_plan")
  phi_or_local_only <- api$evaluate_or_baseline_if_chunk(or_base, seq_len(41L), seq_along(or_base$links), "frame", FALSE, plan_local_only)$phi
  expect_false(isTRUE(all.equal(phi_or, phi_or_local_only, tolerance = 1e-12)))

  # (b) Wrong D and wrong alpha_b controls:
  L_vec <- vapply(or_base$regressions, `[[`, numeric(1L), "L")
  D_wrong <- sum(L_vec) # sum(L_b) instead of sum(v_b * L_b^2)
  beta_wrong_D <- ref$beta * (ref$D / D_wrong)
  expect_gt(max(abs(or_base$shared$beta - beta_wrong_D)), 1e-1)
  alpha_no_slope <- vapply(or_base$regressions, `[[`, numeric(1L), "mean_y") # omitting -L_b * zbar' * beta
  expect_gt(max(abs(vapply(or_base$regressions, `[[`, numeric(1L), "alpha") - alpha_no_slope)), 1e-2)

  # (c) Level-averaging before differencing under large +1e15 per-unit/date offsets:
  # Subtracting raw endpoints per unit first cancels a large unit-level base offset exactly,
  # whereas summing/averaging levels across units with slightly perturbed weights before differencing amplifies roundoff.
  d_big <- d0; d_big$y <- d_big$y + 1e14
  fx_big <- od_fixture(d = d_big, chunk = 7L, X = des2$X, V = des2$X)
  fit_big <- od_or(fx_big, recipe = des2$recipe)
  expect_lt(max(abs(fit_big$shared$beta - or_base$shared$beta)), 1e-2)
})

# The reduced reporting overhead retains every schema and numeric predicate.
test_that("whole-object comparisons detect schema and numeric mutations", {
  expect_failure(od_expect_route_match(list(a = 1), list(a = 1L)))
  expect_failure(od_expect_route_match(list(a = 1), list(a = 2), tol = 1e-12))
  expect_failure(od_expect_route_match(list(a = matrix(1, 1L)), list(a = 1)))
  expect_failure(od_expect_route_match(list(a = 1), list(b = 1)))
})
