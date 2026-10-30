# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 stage 7: complete-roster DR treatment-summary batch
# (`dr_baseline_treatment_batch` and `C_gmm_cp08_dr_treatment_batch`).
# Expected values come from the unchanged per-model complete-support route,
# the unchanged general (non-complete-support) route, the independent dense
# `fit_dr_baseline_table` reference, and central finite differences of unit
# case weights. Counter clones (`sb_counted`, `sb_per_model_ns`, `sb_general_ns`)
# are test-local and disclosed separately from genuine public `cdid::att_gt`
# calls.
#
# STATUS: UNEXECUTED at authoring time (no R, build, install or fit was run in
# the authoring job). Root owns sole R and measured validation.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Fixtures and test-local counter clones (disclosed; separate from public gate)
# ---------------------------------------------------------------------------

sb_ctrl <- function(chunk = 7L, work = 1e11, matrix = 2e7, retained = 1e7, output = 1e7) {
  api$ipw_owner_control(
    max_chunk_rows = chunk,
    max_retained_entries = retained,
    max_matrix_entries = matrix,
    max_work_ops = work,
    max_output_entries = output
  )
}
sb_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
sb_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))
sb_kernels <- c(
  "C_gmm_cp08_support_chunk", "C_gmm_cp08_pair_rows", "C_gmm_cp08_fit_pass1",
  "C_gmm_cp08_fit_pass2", "C_gmm_cp08_or_if_chunk", "C_gmm_cp08_dr_if_prepare",
  "C_gmm_cp08_dr_if_assemble", "C_gmm_cp08_dr_treatment_batch"
)

sb_panel <- function(N = 41L, T = 5L, seed = 19L, extra_cohort = FALSE, second_target = FALSE) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:4] <- c(3, Inf, 3, Inf)
  if (extra_cohort || second_target) g[5:9] <- 4
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
  d
}

sb_pairs <- function(T = 5L) {
  pr <- list()
  for (s in seq_len(T - 1L)) {
    for (t in seq.int(s + 1L, T)) pr[[length(pr) + 1L]] <- c(s, t)
  }
  pr
}

sb_design_spec <- function(N = 41L, kind = "icpt") {
  u <- seq_len(N)
  x1 <- as.numeric((u %% 2L) == 1L)
  x2 <- cos(u * 0.37) + 0.25 * sin(u * 0.73)
  if (identical(kind, "icpt")) {
    X <- matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)"))
    rec <- sb_icpt
    trt_cols <- "(Intercept)"
  } else if (identical(kind, "x_bin")) {
    X <- cbind(`(Intercept)` = 1, x = x1)
    rec <- api$or_pooled_outcome_recipe(c("(Intercept)", "x"), TRUE, "x")
    trt_cols <- c("(Intercept)", "x")
  } else if (identical(kind, "x_cont")) {
    X <- cbind(`(Intercept)` = 1, x2 = x2)
    rec <- api$or_pooled_outcome_recipe(c("(Intercept)", "x2"), TRUE, "x2")
    trt_cols <- c("(Intercept)", "x2")
  } else {
    X <- cbind(`(Intercept)` = 1, x = x1, x2 = x2)
    rec <- api$or_pooled_outcome_recipe(c("(Intercept)", "x", "x2"), TRUE, c("x", "x2"))
    trt_cols <- c("(Intercept)", "x", "x2")
  }
  list(X = X, recipe = rec, trt_cols = trt_cols)
}

sb_fixture <- function(d = sb_panel(), chunk = 7L, perm = FALSE, pairs = sb_pairs(max(d$t)),
                       char_ids = FALSE, wt = 1, X = NULL, V = NULL, ctrl = sb_ctrl(chunk),
                       frame_only = FALSE, K = Inf, targets = 3) {
  if (char_ids) d$id <- sprintf("unit-%03d", d$id)
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  u <- d[!duplicated(d$id), c("id", "g")]
  u <- u[order(u$id), ]
  if (frame_only) {
    extra_id <- if (char_ids) "unit-999" else 999L
    u <- rbind(u, data.frame(id = extra_id, g = Inf))
  }
  if (perm) {
    set.seed(13L)
    u <- u[sample(nrow(u)), ]
  }
  N <- nrow(u)
  ord <- match(u$id, sort(u$id))
  Z <- if (is.null(X)) matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)")) else X[ord, , drop = FALSE]
  V_mat <- if (is.null(V)) Z else V[ord, , drop = FALSE]
  w_vec <- if (length(wt) == 1L) rep(wt, N) else wt[ord]
  fo <- api$construct_ipw_frame_owner(
    po, u$id, Z, weights = w_vec, control = ctrl, frame_g = u$g,
    allow_frame_only = frame_only, design_on_A_only = !frame_only
  )
  reader <- api$construct_cp04_frame_reader(
    fo, u$id, V_mat, u$g, w_vec, ctrl, allow_unknown_off_A = TRUE
  )
  views <- list()
  idx <- 0L
  for (g_t in targets) {
    for (k in seq_along(pairs)) {
      idx <- idx + 1L
      p <- pairs[[k]]
      views[[idx]] <- api$construct_ipw_link_view(
        fo, p[1], p[2], g = g_t, K = K, control = ctrl,
        link_id = paste0("L", idx, ":g", g_t, ":", p[1], "_", p[2]),
        functional_key = paste0("cmpcat:sb:", idx, ":g", g_t, ":", p[1], "-", p[2])
      )
    }
  }
  list(d = d, po = po, fo = fo, reader = reader, views = views, ctrl = ctrl, chunk = chunk, u = u, pairs = pairs)
}

# Disclosed test-local clone forcing the unbatched per-model complete-support loop.
sb_per_model_ns <- function() {
  e <- new.env(parent = api)
  for (nm in c("run_chained_ipw_pipeline", "build_cp04_public_registry", "dr_public_fit", "fit_dr_baseline")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  e$dr_baseline_treatment_batch <- function(...) NULL
  e
}

# Disclosed test-local clone forcing the general (non-complete-support) route.
sb_general_ns <- function() {
  e <- new.env(parent = api)
  for (nm in c("run_chained_ipw_pipeline", "build_cp04_public_registry", "or_public_fit",
               "dr_public_fit", "fit_or_baseline", "fit_dr_baseline")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  e$or_baseline_complete_facts <- function(...) NULL
  e
}

# Disclosed test-local counter clone around the OR/DR fit and IF call graph.
sb_counted <- function(mode = "batch") {
  e <- new.env(parent = api)
  for (nm in c("fit_or_baseline", "fit_dr_baseline", "or_baseline_support", "or_baseline_visit",
               "dr_baseline_support", "dr_support_visit", "dr_baseline_treatment_fit",
               "dr_baseline_treatment_batch", "dr_baseline_shared_logit",
               "or_baseline_complete_supports", "or_baseline_complete_state",
               "or_pooled_fit_common", "replay_or_baseline_if", "replay_dr_baseline_if",
               "evaluate_or_baseline_if_chunk", "evaluate_dr_baseline_if_chunk",
               "or_baseline_native_if_chunk", "dr_baseline_native_if_chunk")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  n <- new.env(parent = emptyenv())
  n$calls <- setNames(integer(length(sb_kernels)), sb_kernels)
  n$logit <- n$trt_fit <- n$trt_batch <- n$search <- n$complete <- 0L
  n$batch_slots <- n$batch_reads <- n$batch_models <- 0
  addr <- lapply(sb_kernels, api$or_baseline_native_kernel)
  e$.Call <- function(.NAME, ...) {
    k <- which(vapply(addr, identical, logical(1), .NAME))
    res <- base::.Call(.NAME, ...)
    if (length(k)) {
      n$calls[k] <- n$calls[k] + 1L
      if (identical(sb_kernels[k], "C_gmm_cp08_dr_treatment_batch") && !is.null(res$counts)) {
        n$batch_slots <- n$batch_slots + res$counts[1L]
        n$batch_reads <- n$batch_reads + res$counts[2L]
        n$batch_models <- n$batch_models + res$counts[3L]
      }
    }
    res
  }
  e$cp04_logit <- function(...) { n$logit <- n$logit + 1L; api$cp04_logit(...) }
  orig_trt_fit <- e$dr_baseline_treatment_fit
  e$dr_baseline_treatment_fit <- function(...) { n$trt_fit <- n$trt_fit + 1L; orig_trt_fit(...) }
  orig_trt_batch <- e$dr_baseline_treatment_batch
  e$dr_baseline_treatment_batch <- function(...) {
    n$trt_batch <- n$trt_batch + 1L
    if (identical(mode, "per_model")) return(NULL)
    orig_trt_batch(...)
  }
  e$cp04_pair_rows <- function(...) { n$search <- n$search + 1L; api$cp04_pair_rows(...) }
  e$or_baseline_complete_pair_rows <- function(...) { n$complete <- n$complete + 1L; api$or_baseline_complete_pair_rows(...) }
  if (identical(mode, "general")) e$or_baseline_complete_facts <- function(...) NULL
  e$n <- n
  e
}

sb_fit_dr <- function(fx, ns = api, ctrl = fx$ctrl, recipe = sb_icpt, trt_cols = NULL, glm = sb_glm) {
  ns$fit_dr_baseline(
    fx$reader, fx$views, ctrl, glm_control = glm, rank_tol = 1e-10,
    outcome_recipe = recipe, treatment_columns = trt_cols
  )
}

sb_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (is.list(x)) {
    a <- attributes(x)
    x <- lapply(x, sb_strip)
    attributes(x) <- a
  }
  x
}

sb_norm_shared <- function(s) {
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

# Route work is a resource contract, not an equal-estimator output. Recompute
# its reservation from input rows, catalog dimensions and the documented compact
# payload convention; do not call the implementation's charge/payload helpers.
sb_expect_work <- function(b, treatment_route = "batch") {
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

sb_backend_view <- function(b, treatment_route = "batch") {
  sb_expect_work(b, treatment_route)
  x <- unclass(b)[setdiff(names(b), c("ref", "reader_ref", "complete_support", "accounting"))]
  if (!is.null(x$shared)) x$shared <- sb_norm_shared(x$shared)
  list(
    fields = sb_strip(x),
    accounting = b$accounting[setdiff(
      names(b$accounting),
      c("complete_support", "work_ops_bound", "own_work_ops", "fit_scratch_bound", "retained_bound", "owned_entries")
    )]
  )
}

sb_expect_match <- function(a, b, tol = 1e-12, info = NULL) {
  mismatch <- function(x, y, path) {
    if (!identical(class(x), class(y)) || !identical(typeof(x), typeof(y)) ||
        !identical(names(x), names(y)) || !identical(dim(x), dim(y)) ||
        !identical(dimnames(x), dimnames(y))) {
      return(paste(path, "schema"))
    }
    if (is.list(x)) {
      if (length(x) != length(y)) return(paste(path, "length"))
      return(unlist(lapply(seq_along(x), function(i) {
        mismatch(x[[i]], y[[i]], paste0(path, "$", if (is.null(names(x))) i else names(x)[i]))
      }), use.names = FALSE))
    }
    same <- if (typeof(x) == "double") {
      if (tol == 0) identical(x, y) else length(waldo::compare(x, y, tolerance = tol)) == 0L
    } else {
      identical(x, y)
    }
    if (same) character() else paste(path, "value")
  }
  bad <- mismatch(a, b, "root")
  expect_true(length(bad) == 0L, info = paste(c(info, head(bad, 12L)), collapse = "\n"))
}

# The OR QR/state is byte-identical to the original per-model control.
# Generic QR factors/centering records need not have the same representation.
# Compare the generic functional and complete unit IF at the established gate.
sb_expect_generic <- function(batch, original, general, info = NULL) {
  expect_identical(batch$shared, original$shared, info = info)
  expect_identical(batch$regressions, original$regressions, info = info)
  sb_expect_match(batch$treatments, original$treatments, tol = 1e-12, info = info)
  sb_expect_match(batch$links, original$links, tol = 1e-12, info = info)
  expect_identical(vapply(batch$links, `[[`, "", "link_id"), vapply(general$links, `[[`, "", "link_id"))
  for (field in c("mu_dr", "mu_ipw", "mu_or"))
    expect_equal(vapply(batch$links, `[[`, 0.0, field), vapply(general$links, `[[`, 0.0, field), tolerance = 1e-12, info = paste(info, field))
  reader <- get("cp04_reader", batch$reader_ref); N <- reader$state$N
  for (lo in seq.int(1L, N, by = batch$control$max_chunk_rows)) {
    rows <- seq.int(lo, min(N, lo + batch$control$max_chunk_rows - 1L));links <- seq_along(batch$links)
    a <- api$evaluate_dr_baseline_if_chunk(batch, rows, links, details = TRUE)
    p <- api$evaluate_dr_baseline_if_chunk(original, rows, links, details = TRUE)
    b <- api$evaluate_dr_baseline_if_chunk(general, rows, links, details = TRUE)
    sb_expect_match(a$phi, p$phi, tol = 1e-12, info = info)
    sb_expect_match(a$details, p$details, tol = 1e-12, info = info)
    expect_equal(a$phi, b$phi, tolerance = 1e-12, info = info)
    for (field in c("psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw")) {
      expect_true(is.matrix(a$details[[field]]), info = field)
      expect_equal(a$details[[field]], b$details[[field]], tolerance = 1e-12, info = paste(info, field))
    }
  }
}

# Covariances must index the same ORIGINAL moments. The generic route may
# retain all M rows while the complete route uses the existing sealed E map.
sb_full_moment_cov <- function(res) {
  V <- res$covariance_stream$V; cert <- res$structural_certificate
  if (is.null(cert) || !isTRUE(res$reduction_summary$reduction_active)) return(unname(V))
  M <- length(res$point_links_all)
  left <- vapply(seq_len(ncol(V)), function(j) api$cp04_expand_reduction_map(cert, V[, j]), numeric(M))
  t(vapply(seq_len(M), function(i) api$cp04_expand_reduction_map(cert, left[i, ]), numeric(M)))
}

sb_run_pipeline <- function(d, chunk = 7L, reduce = TRUE, weighting = "optimal",
                            xformla = NULL, condition_trigger = NULL, clustered = TRUE,
                            ridge_relative = NULL, ns = api) {
  ctrl_list <- if (is.null(chunk)) list() else list(chunk_rows = chunk)
  if (!is.null(condition_trigger)) ctrl_list$condition_trigger <- condition_trigger
  if (!is.null(ridge_relative)) ctrl_list$ridge_relative <- ridge_relative
  if (!length(ctrl_list)) ctrl_list <- NULL
  resolved <- api$chained_public_resolve_control(ctrl_list, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  prep <- api$chained_request_prepare_inputs(
    data = d, yname = "y", tname = "t", idname = "id", gname = "g",
    xformla = xformla, weightsname = NULL, clustervars = if (clustered) "cl" else NULL, anticipation = 0,
    frame = NULL, outcome_regression = FALSE, doubly_robust = TRUE, or_resource = resolved$profile
  )
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = 99L, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "t", yname = "y", gname = "g",
    control = built$control, xformla = prep$xformla, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = prep$target_cohorts, anticipation = 0L, control_group = "nevertreated", base_period = "varying",
    weighting = weighting, reduce_structural_redundancy = reduce, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = c("simple", "group", "dynamic", "calendar"), cband = FALSE, alp = 0.05, biters = 99L, seed = 101L,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = FALSE, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = "post", pt_start = NULL, cp04 = prep$cp04
  )
}

# ---------------------------------------------------------------------------
# 1. Batch vs per-model complete-support route vs general route across ~1, ~x,
#    ~x2, ~x + x2, chunk sizes, permutations, repeated pairs and full pipeline
# ---------------------------------------------------------------------------

test_that("DR treatment summary batch matches per-model complete-support and general routes across ~1, ~x, ~x2, ~x + x2", {
  ns_pm <- sb_per_model_ns()
  ns_gen <- sb_general_ns()
  pairs_rep <- c(sb_pairs(5L), list(c(2, 4), c(1, 5), c(2, 4)))

  for (kind in c("icpt", "x_bin", "x_cont", "x_two")) {
    des <- sb_design_spec(41L, kind)
    for (case in list(
      list(chunk = 7L, perm = FALSE, char = FALSE, wt = 1, pairs = sb_pairs(5L)),
      list(chunk = 16L, perm = TRUE, char = TRUE, wt = 2.5, pairs = pairs_rep),
      list(chunk = 1L, perm = FALSE, char = FALSE, wt = 0.375, pairs = sb_pairs(5L)),
      list(chunk = 41L, perm = TRUE, char = FALSE, wt = 1, pairs = sb_pairs(5L))
    )) {
      info <- paste(kind, "chunk", case$chunk, "perm", case$perm)
      fx <- sb_fixture(
        chunk = case$chunk, perm = case$perm, pairs = case$pairs,
        char_ids = case$char, wt = case$wt, X = des$X, V = des$X
      )
      fit_b <- sb_fit_dr(fx, api, recipe = des$recipe, trt_cols = des$trt_cols)
      fit_p <- sb_fit_dr(fx, ns_pm, recipe = des$recipe, trt_cols = des$trt_cols)
      fit_g <- sb_fit_dr(fx, ns_gen, recipe = des$recipe, trt_cols = des$trt_cols)

      expect_false(is.null(fit_b$complete_support), info = info)
      expect_false(is.null(fit_p$complete_support), info = info)
      expect_null(fit_g$complete_support, info = info)

      T_count <- length(fit_b$treatment_specs)
      expect_identical(fit_b$accounting$complete_support$n_treatment_models, T_count, info = info)
      expect_identical(fit_b$accounting$complete_support$n_treatment_logit_fits, 1L, info = info)
      expect_identical(fit_b$accounting$complete_support$n_treatment_logit_reuses, T_count - 1L, info = info)
      expect_identical(fit_b$accounting$complete_support$n_treatment_summary_batches, 1L, info = info)
      expect_identical(fit_b$accounting$complete_support$n_treatment_summary_batched_models, T_count - 1L, info = info)
      expect_identical(fit_p$accounting$complete_support$n_treatment_summary_batches, 0L, info = info)
      expect_identical(fit_p$accounting$complete_support$n_treatment_summary_batched_models, 0L, info = info)

      # Established 1e-12 agreement with the per-model complete-support route on the same share$prep,
      # and route agreement with the general non-complete-support route across all 10/13 links and models.
      sb_expect_match(fit_b$treatments, fit_p$treatments, tol = 1e-12, info = paste(info, "treatments vs per-model"))
      sb_expect_match(fit_b$links, fit_p$links, tol = 1e-12, info = paste(info, "links vs per-model"))
      sb_expect_match(sb_backend_view(fit_b), sb_backend_view(fit_p, "declined"), tol = 1e-12, info = paste(info, "backend vs per-model"))
      sb_expect_generic(fit_b, fit_p, fit_g, info = info)

      if (!identical(kind, "icpt")) {
        coef1 <- fit_b$treatments[[1L]]$coefficients
        expect_gt(max(abs(coef1[-1L])), 1e-3)
        expect_identical(length(unique(lapply(fit_b$treatments, `[[`, "coefficients"))), 1L)
        for (t in seq_len(T_count)) {
          expect_identical(names(fit_b$treatments[[t]]$d_bar_plus), des$recipe$slope_columns)
          expect_identical(names(fit_b$treatments[[t]]$d_bar_minus), des$recipe$slope_columns)
          expect_identical(names(fit_b$treatments[[t]]$d_bar_diff), des$recipe$slope_columns)
          expect_identical(dimnames(fit_b$treatments[[t]]$partials$beta), list(c("T", "M", "dp", "dm"), des$trt_cols))
          expect_identical(dimnames(fit_b$treatments[[t]]$partials_ipw$beta), list(c("T", "M", "dp", "dm"), des$trt_cols))
        }
      } else {
        for (t in seq_len(T_count)) {
          expect_identical(fit_b$treatments[[t]]$d_bar_plus, numeric(0))
          expect_identical(fit_b$treatments[[t]]$d_bar_minus, numeric(0))
          expect_identical(fit_b$treatments[[t]]$d_bar_diff, numeric(0))
        }
      }

      # Full unit IF (detailed and non-detailed, frame and core scales) and streamed covariance.
      L <- length(fx$views)
      for (details in c(FALSE, TRUE)) {
        for (scale in c("frame", "core")) {
          cb <- api$evaluate_dr_baseline_if_chunk(fit_b, seq_len(min(11L, 41L, case$chunk)), seq_len(L), scale, details)
          cp <- api$evaluate_dr_baseline_if_chunk(fit_p, seq_len(min(11L, 41L, case$chunk)), seq_len(L), scale, details)
          cg <- api$evaluate_dr_baseline_if_chunk(fit_g, seq_len(min(11L, 41L, case$chunk)), seq_len(L), scale, details)
          sb_expect_match(cb$phi, cp$phi, tol = 1e-12, info = paste(info, "phi vs per-model", scale, details))
          sb_expect_match(cb$phi, cg$phi, tol = 1e-12, info = paste(info, "phi vs general", scale, details))
          if (details) {
            sb_expect_match(cb$details, cp$details, tol = 1e-12, info = paste(info, "details vs per-model", scale))
            sb_expect_match(cb$details, cg$details, tol = 1e-12, info = paste(info, "details vs general", scale))
          }
        }
      }
      cov_b <- api$stream_dr_covariance(fit_b)
      cov_p <- api$stream_dr_covariance(fit_p)
      cov_g <- api$stream_dr_covariance(fit_g)
      for (nm in c("estimates", "covariance", "mean_phi", "chunk_rows")) {
        sb_expect_match(cov_b[[nm]], cov_p[[nm]], tol = 1e-12, info = paste(info, "cov vs per-model", nm))
        sb_expect_match(cov_b[[nm]], cov_g[[nm]], tol = 1e-12, info = paste(info, "cov vs general", nm))
      }
    }
  }
})

test_that("full DR pipeline with treatment summary batch preserves global beta, nonconstant PS, full IF, E'WE metric, ridge, ATT, V and SE", {
  ns_pm <- sb_per_model_ns()
  ns_gen <- sb_general_ns()
  d <- sb_panel(41L, 5L, seed = 29L)
  d$g[d$g == Inf] <- 0

  for (xf in list(NULL, ~x, ~x2, ~x + x2)) {
    for (reduce in c(TRUE, FALSE)) {
      for (wt_mode in c("optimal", "identity")) {
        trig <- if (identical(wt_mode, "optimal")) 0.99 else NULL
        rrel <- if (identical(wt_mode, "optimal")) 1e-4 else NULL
        info <- paste(deparse(xf), "reduce", reduce, wt_mode)

        res_b <- sb_run_pipeline(d, chunk = 7L, reduce = reduce, weighting = wt_mode,
          xformla = xf, condition_trigger = trig, ridge_relative = rrel, ns = api)
        res_p <- sb_run_pipeline(d, chunk = 7L, reduce = reduce, weighting = wt_mode,
          xformla = xf, condition_trigger = trig, ridge_relative = rrel, ns = ns_pm)
        res_g <- sb_run_pipeline(d, chunk = 7L, reduce = reduce, weighting = wt_mode,
          xformla = xf, condition_trigger = trig, ridge_relative = rrel, ns = ns_gen)

        b_b <- res_b$nuisance_registry$cp04
        expect_identical(b_b$accounting$complete_support$n_treatment_summary_batches, 1L, info = info)
        expect_identical(b_b$accounting$complete_support$n_treatment_summary_batched_models, 9L, info = info)

        if (identical(wt_mode, "optimal")) {
          expect_identical(res_b$covariance_handling$handler_used, "tikhonov")
          expect_true(res_b$covariance_handling$lambda > 0)
          expect_identical(res_b$covariance_handling$dimensions$M_full, 10L)
        }

        # Retained main and auxiliary comparisons, target metadata, estimates, SE, covariance, GMM objective/ridge.
        sb_expect_match(res_b$point_links_all, res_p$point_links_all, tol = 1e-12, info = info)
        expect_identical(res_b$nuisance_registry$cp04_aux_index, res_p$nuisance_registry$cp04_aux_index, info = info)
        sb_expect_match(sb_strip(res_b$auxiliary_reporting), sb_strip(res_p$auxiliary_reporting), tol = 1e-12, info = info)
        expect_identical(res_b$reporting$target_meta, res_g$reporting$target_meta, info = info)

        for (other in list(res_p, res_g)) {
          expect_equal(res_b$reporting$estimates, other$reporting$estimates, tolerance = 1e-12, info = info)
          expect_equal(res_b$reporting$se, other$reporting$se, tolerance = 1e-12, info = info)
          expect_true(is.matrix(res_b$reporting$if_full) && nrow(res_b$reporting$if_full) == 41L)
          expect_equal(res_b$reporting$if_full, other$reporting$if_full, tolerance = 1e-12, info = info)
          expect_equal(res_b$reporting$vcov, other$reporting$vcov, tolerance = 1e-12, info = info)
          expect_equal(res_b$covariance_handling$lambda, other$covariance_handling$lambda, tolerance = 1e-12, info = info)
          expect_equal(sb_full_moment_cov(res_b), sb_full_moment_cov(other), tolerance = 1e-12, info = info)
          # Actual accepted GMM schema; no NULL-vs-NULL surrogate fields.
          for (field in c("targets", "parameters", "sandwich_cov", "raw_Z1tZ1", "raw_quadratic_residual")) {
            expect_true(length(res_b$gmm_fit[[field]]) > 0L, info = paste(info, field))
            expect_equal(res_b$gmm_fit[[field]], other$gmm_fit[[field]], tolerance = 1e-12, info = paste(info, field))
          }
          for (agg_nm in names(res_b$aggregations)) {
            expect_gt(length(res_b$aggregations[[agg_nm]]$overall$estimate), 0L)
            expect_gt(length(res_b$aggregations[[agg_nm]]$overall$se), 0L)
            expect_equal(res_b$aggregations[[agg_nm]], other$aggregations[[agg_nm]], tolerance = 1e-12, info = info)
          }
        }

        if (!is.null(xf) && reduce && identical(wt_mode, "optimal")) {
          N <- res_b$nuisance_registry$cp04_reader$state$N
          psi_or <- psi_ps <- psi_dir <- matrix(0, N, 10L)
          cat_idx <- vapply(res_b$point_links_all, `[[`, integer(1L), "cp04_index")
          api$replay_dr_baseline_if(b_b, 7L, 10L, function(batch) {
            dst <- match(batch$links, cat_idx)
            take <- which(!is.na(dst))
            if (!length(take)) return(invisible(NULL))
            psi_or[batch$rows, dst[take]] <<- batch$details$psi_or[, take, drop = FALSE]
            psi_ps[batch$rows, dst[take]] <<- batch$details$psi_ps[, take, drop = FALSE]
            psi_dir[batch$rows, dst[take]] <<- batch$details$psi_dir[, take, drop = FALSE]
          }, scale = "frame", details = TRUE)
          expect_gt(max(abs(psi_ps)), 1e-6)
          if (!identical(all.vars(xf), "x")) {
            expect_gt(max(abs(psi_or)), 1e-6)
          }
        }
      }
    }
  }
})

# ---------------------------------------------------------------------------
# 2. Independent dense table reference and case-weight finite-difference checks
# ---------------------------------------------------------------------------

test_that("independent fit_dr_baseline_table and case-weight Gateaux finite differences agree with batched summaries", {
  for (kind in c("icpt", "x_two")) {
    des <- sb_design_spec(41L, kind)
    fx <- sb_fixture(chunk = 11L, perm = TRUE, X = des$X, V = des$X)
    st <- fx$reader$state
    N <- st$N
    fit_b <- sb_fit_dr(fx, api, recipe = des$recipe, trt_cols = des$trt_cols)

    y_at <- function(id, t) fx$d$y[match(paste(id, t), paste(fx$d$id, fx$d$t))]
    Y <- sapply(seq_len(5L), function(t) y_at(st$id, t))
    links_tbl <- lapply(fx$pairs, function(p) {
      list(s = p[1], t = p[2], g = 3, K = Inf, link_id = paste0("L", p[1], "_", p[2]))
    })
    Z_ord <- fx$fo$Z_a[st$a_index, , drop = FALSE]
    tb <- api$fit_dr_baseline_table(
      X_or = Z_ord, V = Z_ord, Y = Y, obs = matrix(1, N, 5L), g = st$g, w = st$w,
      calendar = 1:5, links = links_tbl, outcome_recipe = des$recipe,
      glm_control = sb_glm, rank_tol = 1e-10
    )

    expect_lt(max(abs(vapply(fit_b$links, `[[`, numeric(1), "mu_dr") - tb$estimates)), 1e-12)
    expect_lt(max(abs(vapply(fit_b$links, `[[`, numeric(1), "mu_ipw") - tb$estimates_ipw)), 1e-12)
    if (des$recipe$p > 0L) {
      expect_lt(max(abs(unname(fit_b$shared$beta) - unname(tb$shared$beta))), 1e-12)
    }
    batches <- list()
    api$replay_dr_baseline_if(fit_b, 11L, 10L, function(b) {
      batches[[length(batches) + 1L]] <<- b
    }, details = TRUE)
    col_if <- function(nm) do.call(rbind, lapply(batches, function(b) if (nm == "phi") b$phi else b$details[[nm]]))
    for (nm in c("phi", "psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw")) {
      expect_lt(max(abs(unname(col_if(nm)) - unname(tb[[nm]]))), 1e-12, label = paste(kind, nm))
    }
  }

  # Central finite differences of unit case weights on a small deterministic panel
  # for batched models t = 2..10: verify mu_dr, mu_ipw and analytic partials.
  d_sm <- sb_panel(16L, 5L, seed = 43L)
  des_sm <- sb_design_spec(16L, "x_two")
  fx0 <- sb_fixture(d = d_sm, chunk = 5L, X = des_sm$X, V = des_sm$X)
  fit0 <- sb_fit_dr(fx0, api, recipe = des_sm$recipe, trt_cols = des_sm$trt_cols)
  N_sm <- fx0$reader$state$N
  H_sm <- matrix(0, N_sm, 10L)
  api$replay_dr_baseline_if(fit0, 5L, 10L, function(b) {
    H_sm[b$rows, b$links] <<- b$phi
  }, scale = "frame", details = FALSE)

  eps <- 1e-6
  fd_mu <- matrix(0, N_sm, 10L)
  y_at_sm <- function(id, t) d_sm$y[match(paste(id, t), paste(d_sm$id, d_sm$t))]
  st_sm <- fx0$reader$state
  Y_sm <- sapply(seq_len(5L), function(t) y_at_sm(st_sm$id, t))
  links_sm <- lapply(fx0$pairs, function(p) {
    list(s = p[1], t = p[2], g = 3, K = Inf, link_id = paste0("L", p[1], "_", p[2]))
  })
  Z_sm <- fx0$fo$Z_a[st_sm$a_index, , drop = FALSE]
  for (i in seq_len(N_sm)) {
    wp_i <- st_sm$w; wp_i[i] <- wp_i[i] + eps
    wm_i <- st_sm$w; wm_i[i] <- wm_i[i] - eps
    tp <- api$fit_dr_baseline_table(
      X_or = Z_sm, V = Z_sm, Y = Y_sm, obs = matrix(1, N_sm, 5L), g = st_sm$g, w = wp_i,
      calendar = 1:5, links = links_sm, outcome_recipe = des_sm$recipe,
      glm_control = sb_glm, rank_tol = 1e-10
    )
    tm <- api$fit_dr_baseline_table(
      X_or = Z_sm, V = Z_sm, Y = Y_sm, obs = matrix(1, N_sm, 5L), g = st_sm$g, w = wm_i,
      calendar = 1:5, links = links_sm, outcome_recipe = des_sm$recipe,
      glm_control = sb_glm, rank_tol = 1e-10
    )
    fd_mu[i, ] <- N_sm * (tp$estimates - tm$estimates) / (2 * eps)
  }
  expect_lt(max(abs(H_sm - fd_mu)), 1e-6)
})

# ---------------------------------------------------------------------------
# 3. Summary dispatch counts, single-pair boundary, OR isolation and resource bounds
# ---------------------------------------------------------------------------

test_that("treatment summary batch dispatches one C_gmm_cp08_pair_rows call and one C_gmm_cp08_dr_treatment_batch call", {
  for (T_periods in c(5L, 8L)) {
    d <- sb_panel(41L, T_periods, seed = 61L)
    des <- sb_design_spec(41L, "x_two")
    fx <- sb_fixture(d = d, chunk = 7L, X = des$X, V = des$X)
    T_count <- length(fx$views)

    ec_b <- sb_counted("batch")
    ec_p <- sb_counted("per_model")
    fit_b <- sb_fit_dr(fx, ec_b, recipe = des$recipe, trt_cols = des$trt_cols)
    fit_p <- sb_fit_dr(fx, ec_p, recipe = des$recipe, trt_cols = des$trt_cols)

    expect_identical(ec_b$n$logit, 1L)
    expect_identical(ec_b$n$trt_fit, 1L)
    expect_identical(ec_b$n$trt_batch, 1L)
    expect_identical(ec_b$n$calls[["C_gmm_cp08_fit_pass1"]], 1L)
    expect_identical(ec_b$n$calls[["C_gmm_cp08_fit_pass2"]], 1L)
    expect_identical(ec_b$n$calls[["C_gmm_cp08_pair_rows"]], 1L)
    expect_identical(ec_b$n$calls[["C_gmm_cp08_dr_treatment_batch"]], 1L)
    expect_identical(ec_b$n$search, 0L)
    expect_identical(ec_b$n$complete, 0L)
    expect_identical(as.integer(ec_b$n$batch_models), T_count - 1L)
    expect_identical(ec_b$n$batch_slots, 2 * (T_count - 1L) * 41)
    expect_identical(ec_b$n$batch_reads, 2 * (T_count - 1L) * 41)

    # Unbatched per-model clone runs dr_baseline_treatment_fit and C_gmm_cp08_pair_rows T_count times.
    expect_identical(ec_p$n$trt_fit, T_count)
    expect_identical(ec_p$n$calls[["C_gmm_cp08_pair_rows"]], T_count)
    expect_identical(ec_p$n$calls[["C_gmm_cp08_dr_treatment_batch"]], 0L)
    sb_expect_match(fit_b$treatments, fit_p$treatments, tol = 1e-12)

    # Single-model roster (T_count == 1): model 1 fits via dr_baseline_treatment_fit, 0 batch calls.
    fx1 <- fx
    fx1$views <- fx$views[1L]
    ec1 <- sb_counted("batch")
    fit1 <- sb_fit_dr(fx1, ec1, recipe = des$recipe, trt_cols = des$trt_cols)
    expect_identical(ec1$n$trt_fit, 1L)
    expect_identical(ec1$n$trt_batch, 0L)
    expect_identical(ec1$n$calls[["C_gmm_cp08_pair_rows"]], 1L)
    expect_identical(ec1$n$calls[["C_gmm_cp08_dr_treatment_batch"]], 0L)
    expect_identical(fit1$accounting$complete_support$n_treatment_summary_batches, 0L)
    expect_identical(fit1$accounting$complete_support$n_treatment_summary_batched_models, 0L)

    # OR fit never invokes treatment logit, pair_rows, or treatment summary batch.
    ec_or <- sb_counted("batch")
    fit_or <- ec_or$fit_or_baseline(fx$reader, fx$views, fx$ctrl, 1e-10, outcome_recipe = des$recipe)
    expect_false(is.null(fit_or$complete_support))
    expect_identical(ec_or$n$logit, 0L)
    expect_identical(ec_or$n$trt_fit, 0L)
    expect_identical(ec_or$n$trt_batch, 0L)
    expect_identical(ec_or$n$calls[["C_gmm_cp08_pair_rows"]], 0L)
    expect_identical(ec_or$n$calls[["C_gmm_cp08_dr_treatment_batch"]], 0L)
  }
})

# ---------------------------------------------------------------------------
# 4. Alias / rank boundaries and raw large offset / overflow
# ---------------------------------------------------------------------------

test_that("alias/rank boundaries, large common offsets and overflow conditions match per-model and general routes", {
  ns_pm <- sb_per_model_ns()
  ns_gen <- sb_general_ns()
  des <- sb_design_spec(41L, "x_two")

  # (a) Partial aliasing on OR design (2 kept out of 4 slopes) with 3-column treatment design.
  X_alias <- cbind(
    des$X[, c("(Intercept)", "x")],
    x_const = 5,
    x2 = des$X[, "x2"],
    x_alias = 2 * des$X[, "x"] - 0.5 * des$X[, "x2"]
  )
  rec_alias <- api$or_pooled_outcome_recipe(
    c("(Intercept)", "x", "x_const", "x2", "x_alias"), TRUE,
    c("x", "x_const", "x2", "x_alias")
  )
  fx_alias <- sb_fixture(X = X_alias, V = X_alias)
  fit_ab <- sb_fit_dr(fx_alias, api, recipe = rec_alias, trt_cols = c("(Intercept)", "x", "x2"))
  fit_ap <- sb_fit_dr(fx_alias, ns_pm, recipe = rec_alias, trt_cols = c("(Intercept)", "x", "x2"))
  fit_ag <- sb_fit_dr(fx_alias, ns_gen, recipe = rec_alias, trt_cols = c("(Intercept)", "x", "x2"))
  expect_true(fit_ab$shared$verification_scan)
  expect_identical(fit_ab$shared$p, 2L)
  expect_identical(fit_ab$shared$basis$kept_columns, c("x", "x2"))
  expect_identical(fit_ab$shared$basis$aliased_columns, c("x_const", "x_alias"))
  expect_identical(fit_ab$accounting$complete_support$n_treatment_summary_batches, 1L)
  sb_expect_match(fit_ab$treatments, fit_ap$treatments, tol = 1e-12)
  sb_expect_generic(fit_ab, fit_ap, fit_ag)

  # (b) Full slope aliasing (recipe$p == 1, shared$p == 0).
  X_full_alias <- cbind(`(Intercept)` = 1, x_const = 3, x = des$X[, "x"])
  rec_full_alias <- api$or_pooled_outcome_recipe(c("(Intercept)", "x_const"), TRUE, "x_const")
  fx_fa <- sb_fixture(X = X_full_alias, V = X_full_alias)
  fit_fab <- sb_fit_dr(fx_fa, api, recipe = rec_full_alias, trt_cols = c("(Intercept)", "x"))
  fit_fap <- sb_fit_dr(fx_fa, ns_pm, recipe = rec_full_alias, trt_cols = c("(Intercept)", "x"))
  fit_fag <- sb_fit_dr(fx_fa, ns_gen, recipe = rec_full_alias, trt_cols = c("(Intercept)", "x"))
  expect_identical(fit_fab$shared$p, 0L)
  expect_identical(fit_fab$treatments[[2L]]$d_bar_diff, numeric(0))
  sb_expect_generic(fit_fab, fit_fap, fit_fag)

  # (c) Treated-only covariate null-space refusal matches general route.
  d0 <- sb_panel()
  u0 <- d0[!duplicated(d0$id), c("id", "g")]
  u0 <- u0[order(u0$id), ]
  X_bad_null <- cbind(des$X, x_trt = ifelse(u0$g == 3, seq_len(41L), 0))
  rec_bad_null <- api$or_pooled_outcome_recipe(c("(Intercept)", "x", "x2", "x_trt"), TRUE, c("x", "x2", "x_trt"))
  fx_bad_null <- sb_fixture(d = d0, X = X_bad_null, V = X_bad_null)
  eb <- tryCatch(sb_fit_dr(fx_bad_null, api, recipe = rec_bad_null, trt_cols = des$trt_cols), error = identity)
  eg <- tryCatch(sb_fit_dr(fx_bad_null, ns_gen, recipe = rec_bad_null, trt_cols = des$trt_cols), error = identity)
  expect_s3_class(eb, "cdid_error")
  expect_identical(class(eb), class(eg))
  expect_identical(conditionMessage(eb), conditionMessage(eg))

  # (d) Raw large offset (y + 1e12): raw subtraction before weighted summation cancels the offset.
  d_off <- d0
  d_off$y <- d_off$y + 1e12
  fx_base <- sb_fixture(d = d0, X = des$X, V = des$X)
  fx_off <- sb_fixture(d = d_off, X = des$X, V = des$X)
  fit_base <- sb_fit_dr(fx_base, api, recipe = des$recipe, trt_cols = des$trt_cols)
  fit_off_b <- sb_fit_dr(fx_off, api, recipe = des$recipe, trt_cols = des$trt_cols)
  fit_off_p <- sb_fit_dr(fx_off, ns_pm, recipe = des$recipe, trt_cols = des$trt_cols)
  fit_off_g <- sb_fit_dr(fx_off, ns_gen, recipe = des$recipe, trt_cols = des$trt_cols)
  sb_expect_match(fit_off_b$treatments, fit_off_p$treatments, tol = 1e-12)
  sb_expect_generic(fit_off_b, fit_off_p, fit_off_g)
  expect_lt(max(abs(vapply(fit_off_b$links, `[[`, numeric(1), "mu_dr") -
                    vapply(fit_base$links, `[[`, numeric(1), "mu_dr"))), 1e-4)

  # (e) Large finite raw changes (+/- 1e155) triggering original-factor QR validation in Pass 1
  #     while remaining finite in the treatment summary batch.
  d_big <- d0
  d_big$y <- ifelse((d_big$t %% 2L) == 0L, 1e155, 0) + 1e140 * d0$y
  fx_big <- sb_fixture(d = d_big, X = des$X, V = des$X)
  fit_big_b <- sb_fit_dr(fx_big, api, recipe = des$recipe, trt_cols = des$trt_cols)
  fit_big_p <- sb_fit_dr(fx_big, ns_pm, recipe = des$recipe, trt_cols = des$trt_cols)
  fit_big_g <- sb_fit_dr(fx_big, ns_gen, recipe = des$recipe, trt_cols = des$trt_cols)
  expect_gt(fit_big_b$accounting$complete_support$original_factor_validations, 0L)
  expect_identical(fit_big_b$accounting$complete_support$n_treatment_summary_batches, 1L)
  sb_expect_match(fit_big_b$treatments, fit_big_p$treatments, tol = 1e-12)
  # This NEW extreme fixture demonstrates the inherited common-QR versus
  # generic residual-cancellation limit (saved diagnostic_original_OR.rds).
  # The batch must preserve the actual original common-fit result, not alter
  # that accepted OR fit to reproduce a different QR representation.
  expect_identical(fit_big_b$shared, fit_big_p$shared)
  expect_identical(fit_big_b$regressions, fit_big_p$regressions)
  for (lo in seq.int(1L, 41L, by = fx_big$ctrl$max_chunk_rows)) {
    rows <- seq.int(lo, min(41L, lo + fx_big$ctrl$max_chunk_rows - 1L))
    a <- api$evaluate_dr_baseline_if_chunk(fit_big_b, rows, seq_along(fit_big_b$links), details = TRUE)
    b <- api$evaluate_dr_baseline_if_chunk(fit_big_p, rows, seq_along(fit_big_p$links), details = TRUE)
    expect_true(all(is.finite(a$phi)))
    sb_expect_match(a$phi, b$phi, tol = 1e-12)
    sb_expect_match(a$details, b$details, tol = 1e-12)
  }

  # (f) First-failure parity on corrupted panel endpoints/differences and nonfinite means/partials.
  rebind <- function(env, nm, value) {
    old <- get(nm, envir = env)
    unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
    function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
  }
  po <- fx_base$po
  y_orig <- po$observations$y
  cell <- function(id, t) po$postings$post_obs_index[po$postings$date_start[t] - 1L + match(id, po$a_id)]
  bad_panels <- list(
    endpoint_nan = { z <- y_orig; z[cell(4, 5)] <- NaN; z },
    diff_overflow = { z <- y_orig; z[cell(4, 2)] <- 1e308; z[cell(4, 3)] <- -1e308; z }
  )
  for (nm in names(bad_panels)) {
    bp <- po; bp$observations$y <- bad_panels[[nm]]
    restore <- rebind(po$ref, "panel_owner", bp)
    eb <- tryCatch(sb_fit_dr(fx_base, api, recipe = des$recipe, trt_cols = des$trt_cols), error = identity)
    ep <- tryCatch(sb_fit_dr(fx_base, ns_pm, recipe = des$recipe, trt_cols = des$trt_cols), error = identity)
    eg <- tryCatch(sb_fit_dr(fx_base, ns_gen, recipe = des$recipe, trt_cols = des$trt_cols), error = identity)
    restore()
    expect_s3_class(eb, "cdid_error")
    expect_identical(class(eb), class(ep), info = nm)
    expect_identical(conditionMessage(eb), conditionMessage(ep), info = nm)
    expect_identical(class(eb), class(eg), info = nm)
    expect_identical(conditionMessage(eb), conditionMessage(eg), info = nm)
  }
})

# ---------------------------------------------------------------------------
# 5. Structural negatives, actual bindings/stale objects, kernel status/guards,
#    and adjacent memory/work caps before allocation
# ---------------------------------------------------------------------------

test_that("structural negatives, broken bindings, stale objects, kernel guards and adjacent caps refuse or fall back cleanly", {
  ns_gen <- sb_general_ns()
  d0 <- sb_panel()
  des <- sb_design_spec(41L, "x_two")

  # (a) Structural negatives decline complete_support and run 0 native batch calls.
  wv <- 1 + (seq_len(41L) %% 3L) / 4
  X_no_icpt <- des$X[, c("x", "x2"), drop = FALSE]
  rec_no_icpt <- api$or_pooled_outcome_recipe(c("x", "x2"), FALSE, c("x", "x2"))
  neg_cases <- list(
    missing_obs = list(fx = sb_fixture(d = d0[!(d0$id == 5 & d0$t == 3), ], X = des$X, V = des$X), rec = des$recipe, trt = des$trt_cols),
    varying_w = list(fx = sb_fixture(d = d0, wt = wv, X = des$X, V = des$X), rec = des$recipe, trt = des$trt_cols),
    frame_only = list(fx = sb_fixture(d = d0, frame_only = TRUE), rec = sb_icpt, trt = "(Intercept)"),
    two_controls = list(fx = sb_fixture(d = sb_panel(extra_cohort = TRUE), K = c(4, Inf)), rec = sb_icpt, trt = "(Intercept)"),
    two_targets = list(fx = sb_fixture(d = sb_panel(second_target = TRUE), targets = c(3, 4)), rec = sb_icpt, trt = "(Intercept)"),
    no_intercept = list(fx = sb_fixture(d = d0, X = X_no_icpt, V = X_no_icpt), rec = rec_no_icpt, trt = c("x", "x2"))
  )
  for (nm in names(neg_cases)) {
    cs <- neg_cases[[nm]]
    ec <- sb_counted("batch")
    fit <- sb_fit_dr(cs$fx, ec, recipe = cs$rec, trt_cols = cs$trt)
    gen <- sb_fit_dr(cs$fx, ns_gen, recipe = cs$rec, trt_cols = cs$trt)
    if (nm %in% c("two_controls", "two_targets")) {
      n_grp <- if (nm == "two_targets") 2L else 1L
      expect_false(is.null(fit$complete_support), info = nm)
      expect_null(gen$complete_support, info = nm)
      expect_identical(ec$n$trt_batch, n_grp, info = nm)
      expect_identical(ec$n$calls[["C_gmm_cp08_dr_treatment_batch"]], n_grp, info = nm)
      expect_gt(sum(ec$n$calls), 0L, label = nm)
      expect_identical(ec$n$logit, n_grp, info = nm)
      sb_expect_match(sb_backend_view(fit), sb_backend_view(gen), tol = 1e-12, info = nm)
    } else {
      expect_null(fit$complete_support, info = nm)
      expect_identical(ec$n$trt_batch, 0L, info = nm)
      expect_identical(ec$n$calls[["C_gmm_cp08_dr_treatment_batch"]], 0L, info = nm)
      expect_identical(sum(ec$n$calls), 0L, info = nm)
      expect_identical(ec$n$logit, length(cs$fx$views), info = nm)
      expect_identical(sb_strip(unclass(fit)[names(fit) != "ref"]), sb_strip(unclass(gen)[names(gen) != "ref"]), info = nm)
    }
  }

  # (b) Direct binding and stale-object negatives on dr_baseline_treatment_batch.
  fx <- sb_fixture(d = d0, X = des$X, V = des$X)
  fit <- sb_fit_dr(fx, api, recipe = des$recipe, trt_cols = des$trt_cols)
  trt_recipe <- fit$treatment_recipe
  share0 <- new.env(parent = emptyenv())
  first_trt <- api$dr_baseline_treatment_fit(
    fx$reader, fit$treatment_specs[[1L]], fit$supports,
    fit$regressions[[fit$treatment_specs[[1L]]$regression]], fit$shared,
    fit$outcome_recipe, fit$or_frame_index, trt_recipe, sb_glm, fx$ctrl,
    NULL, fit$complete_support, share0
  )
  clone_share <- function(sh = share0) {
    e <- new.env(parent = emptyenv())
    for (k in ls(sh, all.names = TRUE)) assign(k, get(k, envir = sh), envir = e)
    e
  }
  run_batch <- function(specs = fit$treatment_specs[-1L], supports = fit$supports,
                        regressions = fit$regressions, shared = fit$shared,
                        recipe = fit$outcome_recipe, or_idx = fit$or_frame_index,
                        t_rec = trt_recipe, gc = sb_glm, ctrl = fx$ctrl,
                        ref = NULL, complete = fit$complete_support, sh = clone_share(),
                        scratch = NULL, call_fn = base::.Call) {
    api$dr_baseline_treatment_batch(
      fx$reader, specs, supports, regressions, shared, recipe, or_idx, t_rec,
      gc, ctrl, ref, complete, sh, scratch = scratch, call_fn = call_fn
    )
  }
  got_ok <- run_batch()
  expect_length(got_ok, length(fit$treatment_specs) - 1L)
  sb_expect_match(got_ok, fit$treatments[-1L], tol = 1e-12)

  batch_calls <- 0L
  counting_call <- function(.NAME, ...) { batch_calls <<- batch_calls + 1L; base::.Call(.NAME, ...) }

  # Stale share counters or mutated prep bindings decline with 0 native calls:
  sh_stale <- clone_share(); sh_stale$reuses <- 1L
  expect_null(run_batch(sh = sh_stale, call_fn = counting_call))
  sh_bad_ref <- clone_share(); sh_bad_ref$prep$reader_ref <- new.env(parent = emptyenv())
  expect_null(run_batch(sh = sh_bad_ref, call_fn = counting_call))
  sh_bad_fref <- clone_share(); sh_bad_fref$prep$frame_ref <- new.env(parent = emptyenv())
  expect_null(run_batch(sh = sh_bad_fref, call_fn = counting_call))
  sh_bad_cols <- clone_share(); sh_bad_cols$prep$kept_columns <- c("x", "wrong")
  expect_null(run_batch(sh = sh_bad_cols, call_fn = counting_call))
  sh_bad_rows <- clone_share(); sh_bad_rows$inputs$rows <- rev(sh_bad_rows$inputs$rows)
  expect_null(run_batch(sh = sh_bad_rows, call_fn = counting_call))
  sh_bad_cr <- clone_share(); sh_bad_cr$prep$cohort_rows[[1L]] <- rev(sh_bad_cr$prep$cohort_rows[[1L]])
  expect_null(run_batch(sh = sh_bad_cr, call_fn = counting_call))

  # Spec / support / regression descriptor mutations decline with 0 native calls:
  sp_bad_coh <- fit$treatment_specs[-1L]; sp_bad_coh[[2L]]$cohorts <- rev(sp_bad_coh[[2L]]$cohorts)
  expect_null(run_batch(specs = sp_bad_coh, call_fn = counting_call))
  sp_bad_n <- fit$treatment_specs[-1L]; sp_bad_n[[2L]]$n_observers <- sp_bad_n[[2L]]$n_observers - 1L
  expect_null(run_batch(specs = sp_bad_n, call_fn = counting_call))
  sp_bad_ref <- fit$treatment_specs[-1L]; sp_bad_ref[[1L]]$reference <- 1L
  expect_null(run_batch(specs = sp_bad_ref, call_fn = counting_call))
  sp_bad_date <- fit$treatment_specs[-1L]; sp_bad_date[[1L]]$pair$s_calendar_index <- 99L
  expect_null(run_batch(specs = sp_bad_date, call_fn = counting_call))
  sp_bad_reg <- fit$treatment_specs[-1L]; sp_bad_reg[[1L]]$regression <- 99L
  expect_null(run_batch(specs = sp_bad_reg, call_fn = counting_call))
  sup_bad_v <- fit$supports; sup_bad_v[[fit$treatment_specs[[3L]]$support_ids[1L]]]$v <- 2
  expect_null(run_batch(supports = sup_bad_v, call_fn = counting_call))
  sup_bad_rate <- fit$supports; sup_bad_rate[[fit$treatment_specs[[3L]]$support_ids[1L]]]$rate <- 0
  expect_null(run_batch(supports = sup_bad_rate, call_fn = counting_call))
  reg_bad_L <- fit$regressions; reg_bad_L[[fit$treatment_specs[[2L]]$regression]]$L <- NaN
  expect_null(run_batch(regressions = reg_bad_L, call_fn = counting_call))
  expect_identical(batch_calls, 0L)

  # Nonfinite DR/IPW means (status 14) and nonfinite cross derivatives (status 15) raise cp04_finite errors:
  reg_huge_alpha <- fit$regressions; reg_huge_alpha[[fit$treatment_specs[[2L]]$regression]]$alpha <- 1e308
  reg_huge_L <- reg_huge_alpha; reg_huge_L[[fit$treatment_specs[[2L]]$regression]]$L <- -1e308
  expect_error(run_batch(regressions = reg_huge_L), "Nonfinite treated/control DR and IPW means",
    class = "cdid_cp04_numeric_error")
  # D93 zeros the rate log-derivative when rate <= ipw_clip_eps (1e-6); an unclipped small rate (> ipw_clip_eps) still overflows:
  sup_tiny_rate <- fit$supports; sup_tiny_rate[[fit$treatment_specs[[2L]]$support_ids[1L]]]$rate <- 1e-308
  sh_huge_ss <- clone_share(); sh_huge_ss$prep$score_sums[] <- 1e308
  before_clipped <- batch_calls
  got_clipped <- run_batch(supports = sup_tiny_rate, sh = sh_huge_ss, call_fn = counting_call)
  expect_type(got_clipped, "list")
  expect_length(got_clipped, length(fit$treatment_specs) - 1L)
  expect_identical(batch_calls, before_clipped + 1L)
  expect_true(all(got_clipped[[1L]]$cross[, 1L] == 0))
  expect_true(all(got_clipped[[1L]]$partials$upstream[, 1L] == 0))
  expect_true(all(got_clipped[[1L]]$partials_ipw$upstream[, 1L] == 0))
  sup_above_clip <- fit$supports; sup_above_clip[[fit$treatment_specs[[2L]]$support_ids[1L]]]$rate <- 2 * api$ipw_clip_eps
  # A successful batch consumes its invocation-local share counters. Use a
  # fresh share so this negative reaches the native numerical check.
  sh_above_clip <- clone_share(); sh_above_clip$prep$score_sums[] <- 1e308
  before_above <- batch_calls
  expect_error(run_batch(supports = sup_above_clip, sh = sh_above_clip, call_fn = counting_call),
    "Nonfinite baseline DR analytic cross derivatives", class = "cdid_cp04_numeric_error")
  expect_identical(batch_calls, before_above + 1L)

  # (c) Adjacent memory, scratch, and work caps refuse BEFORE allocation or native call.
  W_fit <- fit$accounting$work_ops_bound
  ec_ok <- sb_counted("batch")
  fit_at <- sb_fit_dr(fx, ec_ok, ctrl = sb_ctrl(7L, work = W_fit), recipe = des$recipe, trt_cols = des$trt_cols)
  expect_identical(fit_at$accounting$work_ops_bound, W_fit)
  ec_fail <- sb_counted("batch")
  expect_error(
    sb_fit_dr(fx, ec_fail, ctrl = sb_ctrl(7L, work = W_fit - 1), recipe = des$recipe, trt_cols = des$trt_cols),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(c(ec_fail$n$logit, ec_fail$n$trt_fit, ec_fail$n$trt_batch, sum(ec_fail$n$calls)), c(0L, 0L, 0L, 0L))
  # A structural decline cannot use the smaller batch work reservation for
  # the remaining ordinary model loop. Its first fit is already complete.
  ec_decline <- sb_counted("per_model")
  expect_error(sb_fit_dr(fx, ec_decline, ctrl = sb_ctrl(7L, work = W_fit),
    recipe = des$recipe, trt_cols = des$trt_cols), class = "cdid_ipw_owner_bound_error")
  expect_identical(ec_decline$n$trt_fit, 1L)
  expect_identical(ec_decline$n$calls[["C_gmm_cp08_dr_treatment_batch"]], 0L)

  n_obs <- share0$prep$n
  p_trt <- trt_recipe$p
  p_or <- fit$shared$p
  m_coh <- length(share0$inputs$cohorts)
  M_b <- length(fit$treatment_specs) - 1L
  share_entries <- n_obs * (p_trt + p_or + 14) + p_trt * p_trt + p_trt * (m_coh + 1) + 2 * m_coh + 128
  per_model_summary <- 12 + 3 * p_or + p_trt * m_coh + 8 * (m_coh + p_trt) + 16
  batch_native_entries <- 2 * n_obs + M_b * (per_model_summary + 4 + m_coh) + 128
  batch_work <- 64 * M_b * n_obs * (p_or + p_trt + m_coh + 8) + 128 * M_b * per_model_summary

  batch_calls <- 0L
  expect_length(run_batch(scratch = share_entries + batch_native_entries, call_fn = counting_call), M_b)
  expect_identical(batch_calls, 1L)
  batch_calls <- 0L
  expect_error(
    run_batch(scratch = share_entries + batch_native_entries - 1, call_fn = counting_call),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(batch_calls, 0L)

  expect_length(run_batch(ctrl = sb_ctrl(7L, work = batch_work), call_fn = counting_call), M_b)
  expect_identical(batch_calls, 1L)
  batch_calls <- 0L
  expect_error(
    run_batch(ctrl = sb_ctrl(7L, work = batch_work - 1), call_fn = counting_call),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(batch_calls, 0L)

  ret_base <- fx$reader$accounting$owned_entries + fx$reader$accounting$borrowed_owner_entries
  mat_edge <- ret_base + share_entries + batch_native_entries
  expect_length(run_batch(ctrl = sb_ctrl(7L, matrix = mat_edge), call_fn = counting_call), M_b)
  expect_identical(batch_calls, 1L)
  batch_calls <- 0L
  expect_error(
    run_batch(ctrl = sb_ctrl(7L, matrix = mat_edge - 1), call_fn = counting_call),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(batch_calls, 0L)

  # (d) Direct C_gmm_cp08_dr_treatment_batch kernel checks: exact max_entries edge, slot failure and malformed args.
  k_batch <- api$or_baseline_native_kernel("C_gmm_cp08_dr_treatment_batch")
  po <- get("panel_owner", fx$po$ref)
  post <- po$postings
  st <- fx$reader$state
  pairs_mat <- vapply(fit$treatment_specs[-1L], function(s) {
    as.integer(c(s$pair$s_calendar_index, s$pair$t_calendar_index))
  }, integer(2))
  model_par <- vapply(fit$treatment_specs[-1L], function(s) {
    reg <- fit$regressions[[s$regression]]
    c(reg$L, reg$alpha, vapply(fit$supports[s$support_ids], `[[`, numeric(1), "rate"))
  }, numeric(2 + m_coh))
  norm_vec <- as.double(c(share0$prep$sum_wp, share0$prep$sum_wm, share0$prep$dp, share0$prep$dm, st$N))
  ctl_vec <- as.integer(c(7L, isTRUE(unname(capabilities("long.double")))))
  c_exact_entries <- 2 * n_obs + M_b * (12 + 3 * p_or + p_trt * m_coh + 8 * (m_coh + p_trt))

  c_call <- function(...) {
    a <- list(
      po$observations$y, post$date_start, post$post_a_index, post$post_obs_index,
      share0$inputs$rows, st$a_index, share0$inputs$X, share0$prep$Z_or, fit$shared$beta,
      share0$prep$wp, share0$prep$wm, share0$prep$cohort_rows, share0$prep$score_sums,
      share0$prep$wm_design, norm_vec, pairs_mat, model_par, ctl_vec, c_exact_entries
    )
    edits <- list(...)
    for (k in names(edits)) a[[as.integer(k)]] <- edits[[k]]
    do.call(base::.Call, c(list(k_batch), a))
  }
  c_ok <- c_call()
  expect_identical(c_ok$status, c(0L, 0L, 0L, 0L))
  expect_error(c_call(`19` = c_exact_entries - 1), class = "cdid_gmm_bound_error")

  bad_pa <- post$post_a_index
  bad_pa[post$date_start[pairs_mat[1L, 1L]]] <- 999L
  expect_identical(c_call(`3` = bad_pa)$status[1L], 1L)

  bad_args <- list(
    list(`5` = rev(share0$inputs$rows)),
    list(`5` = c(share0$inputs$rows[-1L], 999L)),
    list(`7` = share0$inputs$X[-1L, , drop = FALSE]),
    list(`8` = as.vector(share0$prep$Z_or)),
    list(`9` = fit$shared$beta[-1L]),
    list(`10` = share0$prep$wp[-1L]),
    list(`12` = list(rev(share0$prep$cohort_rows[[1L]]), share0$prep$cohort_rows[[2L]])),
    list(`13` = share0$prep$score_sums[, 1L, drop = FALSE]),
    list(`14` = share0$prep$wm_design[-1L]),
    list(`15` = c(0, norm_vec[-1L])),
    list(`16` = pairs_mat + 100L),
    list(`17` = model_par[-1L, , drop = FALSE]),
    list(`18` = c(0L, 1L))
  )
  for (ed in bad_args) {
    expect_error(do.call(c_call, ed), class = "cdid_gmm_bound_error")
  }
})

# ---------------------------------------------------------------------------
# 6. Genuine public gate (no counter clone)
# ---------------------------------------------------------------------------

test_that("genuine public cdid::att_gt DR calls on ~1 and ~x succeed and match unbatched reference", {
  skip_if_not_installed("did")
  d <- sb_panel(41L, 5L, seed = 73L)
  d$g[d$g == Inf] <- 0
  ns_pm <- sb_per_model_ns()

  for (xf in list(~1, ~x, ~x + x2)) {
    mp <- cdid::att_gt(
      yname = "y", tname = "t", idname = "id", gname = "g", data = d,
      xformla = xf, est_method = "dr", control_group = "nevertreated",
      bstrap = FALSE, cband = FALSE, base_period = "varying",
      mode = "chained", pt_assumption = "post"
    )
    expect_s3_class(mp, "cdid_MP")
    expect_true(all(is.finite(mp$att)))
    expect_true(all(is.finite(mp$se)))
    ref <- sb_run_pipeline(d, chunk = NULL, reduce = TRUE, weighting = "optimal",
      xformla = if (identical(xf, ~1)) NULL else xf, clustered = FALSE, ns = ns_pm)
    expect_equal(as.numeric(mp$att), as.numeric(ref$reporting$estimates), tolerance = 1e-12)
    expect_equal(as.numeric(mp$se), as.numeric(ref$reporting$se), tolerance = 1e-12)
  }
})
