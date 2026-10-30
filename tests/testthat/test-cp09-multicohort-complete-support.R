# D90 multi-cohort baseline OR/DR. UNEXECUTED here: root must run these on the
# exact installed revision. References are the unchanged general route and the
# existing pure-table backend; no candidate formula supplies expected IFs.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
cp09_api <- asNamespace("cdid")

# Continuous, non-saturated designs; unequal cohort masses, interleaved IDs,
# irregular calendar labels and a nonadditive outcome. No RNG is needed.
cp09_panel <- function() {
  N <- 53L; tt <- c(2, 4, 9, 15, 26, 40)
  u <- seq_len(N)
  gg <- rep(c(9, 26, 0), c(13L, 17L, 23L))[order(sin(u * 0.83))]
  d <- data.frame(id = rep(1000L + rev(u), each = length(tt)), period = rep(tt, N))
  i <- rep(u, each = length(tt)); k <- rep(seq_along(tt), N)
  d$G <- gg[i]
  d$x <- sin(i * 0.71) + i / 100
  d$x2 <- cos(i * 0.37) + sin(i * 0.19)
  d$x3 <- sin(i * 0.27) + cos(i * 0.43)
  d$y <- 0.2 * i + 0.1 * k + 0.12 * d$x * k^2 + 0.09 * d$x2 * k +
    sin(i * k * 0.31) + cos(i * k * 0.17) +
    ifelse(d$G > 0 & d$period >= d$G, 1 + k / 7, 0)
  d$cl <- (i - 1L) %/% 3L + 1L
  d
}
cp09_ctrl <- function(chunk = 9L, matrix = 2e7, work = 1e13)
  cp09_api$ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = 2e7,
    max_matrix_entries = matrix, max_work_ops = work, max_output_entries = 2e7)
cp09_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
cp09_recipe <- function(columns = c("(Intercept)", "x", "x2"), intercept = TRUE)
  cp09_api$or_pooled_outcome_recipe(columns, intercept, setdiff(columns, "(Intercept)"))

# Same test-local closure rebinding as the CP08 tests; namespace stays locked.
cp09_ns <- function(general = FALSE, stream = FALSE, cap = NULL) {
  e <- new.env(parent = cp09_api)
  for (nm in c("att_gt", "chained_public_att_gt", "run_chained_ipw_pipeline", "build_cp04_public_registry", "or_public_fit", "dr_public_fit",
               "fit_or_baseline", "fit_dr_baseline")) {
    f <- get(nm, cp09_api); environment(f) <- e; assign(nm, f, e)
  }
  if (general) e$or_baseline_complete_facts <- function(...) NULL
  if (stream) e$cp04_retained_if_reserve <- function(...) NULL
  if (!is.null(cap)) e$cp04_retained_if_reserve <- function(..., control) {
    control$max_retained_entries <- cap
    cp09_api$cp04_retained_if_reserve(..., control = control)
  }
  e
}
cp09_fixture <- function(cg = "notyettreated", est = "dr", chunk = 9L,
                         perm = FALSE, wt = 1, alias = FALSE, no_intercept = FALSE) {
  d <- cp09_panel(); d$G[d$G == 0] <- Inf
  if (perm) d <- d[order(cos(seq_len(nrow(d)) * 0.7)), ]
  ctrl <- cp09_ctrl(chunk)
  po <- cp09_api$construct_ipw_panel_owner(d, "id", "period", "y", "G", ctrl)
  u <- d[!duplicated(d$id), ]; u <- u[order(sin(u$id * 0.3)), ]
  X <- cbind(`(Intercept)` = 1, x = u$x, x2 = u$x2)
  if (alias) X <- cbind(X, x_alias = 2 * u$x)
  recipe <- cp09_recipe(colnames(X))
  if (no_intercept) { X <- X[, -1L, drop = FALSE]; recipe <- cp09_recipe(colnames(X), FALSE) }
  # Private scope permits distinct frozen OR/treatment columns in one frame.
  Z <- if (est == "dr") cbind(X, x3 = u$x3) else X
  ww <- if (length(wt) == 1L) rep(wt, nrow(u)) else wt[match(u$id, sort(u$id))]
  fo <- cp09_api$construct_ipw_frame_owner(po, u$id, Z, weights = ww, control = ctrl,
    frame_g = u$G, allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- cp09_api$construct_cp04_frame_reader(fo, u$id,
    matrix(1, nrow(u), 1L, dimnames = list(NULL, "reader_intercept")), u$G, ww, ctrl,
    allow_unknown_off_A = TRUE)
  tt <- sort(unique(d$period)); views <- lks <- list()
  for (g in c(9, 26)) for (s in seq_len(length(tt) - 1L)) for (t in seq.int(s + 1L, length(tt))) {
    K <- if (cg == "nevertreated") Inf else c(9, 26, Inf)[c(9, 26, Inf) > max(g, tt[t])]
    j <- length(views) + 1L; id <- paste0("L", j)
    views[[j]] <- cp09_api$construct_ipw_link_view(fo, tt[s], tt[t], g = g, K = K, control = ctrl,
      link_id = id, functional_key = paste0("cmpcat:cp09:", j))
    lks[[j]] <- list(s = tt[s], t = tt[t], g = g, K = K, link_id = id)
  }
  if (perm) {
    # Later-cohort links arrive first, with pairs interleaved. Neither the
    # complete target vector nor group traversal may assume cohort order.
    o <- c(seq.int(length(views), 1L, by = -2L), seq.int(length(views) - 1L, 1L, by = -2L))
    views <- views[o]; lks <- lks[o]
  }
  Y <- matrix(0, nrow(u), length(tt))
  for (k in seq_along(tt)) {
    dk <- d[d$period == tt[k], ]; Y[, k] <- dk$y[match(u$id, dk$id)]
  }
  list(d = d, u = u, po = po, fo = fo, reader = reader, views = views, links = lks,
    X = X, V = if (est == "dr") Z[, c(if (!no_intercept) "(Intercept)", "x", "x3"), drop = FALSE] else X,
    Y = Y, calendar = tt, weights = ww, recipe = recipe, ctrl = ctrl, est = est)
}
cp09_fit <- function(fx, ns = cp09_api) {
  if (fx$est == "or") return(ns$fit_or_baseline(fx$reader, fx$views, fx$ctrl,
    rank_tol = 1e-10, outcome_recipe = fx$recipe))
  ns$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = cp09_glm, rank_tol = 1e-10,
    outcome_recipe = fx$recipe, treatment_columns = colnames(fx$V))
}
cp09_chunk <- function(fit, rows, links, details = FALSE, scale = "frame", plan = NULL) {
  fn <- if (inherits(fit, "cdid_dr_baseline")) cp09_api$evaluate_dr_baseline_if_chunk else
    cp09_api$evaluate_or_baseline_if_chunk
  fn(fit, rows, links, details = details, scale = scale, plan = plan)
}
cp09_phi <- function(fit, links = seq_along(fit$links), chunk = fit$control$max_chunk_rows, details = FALSE) {
  st <- get("cp04_reader", fit$reader_ref, inherits = FALSE)$state
  H <- matrix(0, st$N, length(links))
  for (first in seq.int(1L, st$N, by = chunk)) {
    rr <- seq.int(first, min(st$N, first + chunk - 1L))
    H[rr, ] <- cp09_chunk(fit, rr, links, details)$phi
  }
  H
}
cp09_run <- function(d = cp09_panel(), est = "dr", cg = "notyettreated", ns = cp09_api,
                     chunk = 32L, xf = ~x + x2, weightsname = NULL, condition_trigger = NULL,
                     ridge_relative = NULL, base_period = "varying", clustered = TRUE) {
  cc <- list(chunk_rows = chunk)
  if (!is.null(condition_trigger)) cc$condition_trigger <- condition_trigger
  if (!is.null(ridge_relative)) cc$ridge_relative <- ridge_relative
  r <- cp09_api$chained_public_resolve_control(cc, gmm_weighting = "optimal", gmm_covariance = "tikhonov")
  p <- cp09_api$chained_request_prepare_inputs(data = d, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xf, weightsname = weightsname, clustervars = if (clustered) "cl" else NULL, anticipation = 0, frame = NULL,
    outcome_regression = identical(est, "reg"), doubly_robust = identical(est, "dr"), or_resource = r$profile)
  p$cp04$automatic_chunk_rows <- FALSE
  dims <- cp09_api$chained_public_dimensions(p, biters = 99L, clustered = clustered)
  b <- cp09_api$chained_public_pipeline_control(dims, r, gmm_weighting = "optimal", gmm_covariance = "tikhonov")
  ns$run_chained_ipw_pipeline(data = p$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = b$control, xformla = p$xformla, weightsname = NULL, cluster_vec = if (clustered) p$cluster_full else NULL,
    target_cohorts = p$target_cohorts, anticipation = 0L, control_group = cg, base_period = base_period,
    weighting = "optimal", reduce_structural_redundancy = TRUE, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = c("simple", "group", "dynamic", "calendar"), cband = FALSE, alp = 0.05, biters = 99L,
    seed = 101L, frame_id = p$frame_ids, frame_Z = p$frame_Z, frame_weights = p$frame_w, frame_g = p$frame_g_arg,
    allow_frame_only = p$n_frame_only > 0L, na.rm = FALSE, bstrap = FALSE, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = "post", pt_start = NULL, cp04 = p$cp04)
}
cp09_equal <- function(a, b) {
  # Assert populated fields before parity: NULL == NULL is not evidence.
  expect_true(is.matrix(a$reporting$if_full)); expect_true(is.matrix(a$covariance_stream$V))
  for (nm in c("mhat", "V", "mean_phi"))
    expect_equal(a$covariance_stream[[nm]], b$covariance_stream[[nm]], tolerance = 1e-12)
  expect_equal(a$gmm_fit$targets, b$gmm_fit$targets, tolerance = 1e-12)
  for (nm in c("sandwich_cov", "materialized_target_if"))
    expect_equal(a$gmm_fit$if_projection[[nm]], b$gmm_fit$if_projection[[nm]], tolerance = 1e-12)
  # Signs of roundoff-sized nullspace eigenvalues need not agree. Keep the
  # spectrum, numerical rank/rule, trigger, lambda and every policy field.
  ah <- a$covariance_handling; bh <- b$covariance_handling
  costs <- c("producer_workspace_entries", "peak_concurrent_entries", "work_ops_bound")
  for (z in list(a, b)) {
    conditioning <- z$covariance_handling$conditioning
    count <- conditioning$nonpositive_eigenvalue_count
    expect_true(is.numeric(count) && length(count) == 1L && is.finite(count) && count == floor(count))
    expect_gte(count, 0L); expect_lte(count, length(conditioning$eigenvalues_B))
    expect_identical(count, as.integer(sum(conditioning$eigenvalues_B <= 0)))
    # These are route-dependent charges, still finite integral counts within
    # this run's original owner limits; the remaining preflight fields agree.
    pre <- z$auxiliary_reporting$preflight
    ctrl <- z$controls$pipeline_control$owner_control
    expect_true(all(costs %in% names(pre)))
    for (nm in costs) {
      value <- pre[[nm]]
      expect_true(is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= 0 && value == floor(value) && value <= 2^53 - 1)
      expect_lte(value, if (nm == "work_ops_bound") ctrl$max_work_ops else ctrl$max_matrix_entries)
    }
    expect_lte(pre$producer_workspace_entries, pre$peak_concurrent_entries)
  }
  ah$conditioning$nonpositive_eigenvalue_count <- NULL
  bh$conditioning$nonpositive_eigenvalue_count <- NULL
  expect_equal(ah, bh, tolerance = 1e-12)
  expect_equal(a$full_frame_transport$Phi_full, b$full_frame_transport$Phi_full, tolerance = 1e-12)
  for (nm in c("estimates", "se", "analytical_se", "if_full", "vcov"))
    expect_equal(a$reporting[[nm]], b$reporting[[nm]], tolerance = 1e-12)
  expect_identical(a$reporting$target_meta, b$reporting$target_meta)
  ax <- a$auxiliary_reporting; bx <- b$auxiliary_reporting
  ax$preflight[costs] <- NULL; bx$preflight[costs] <- NULL
  expect_equal(ax, bx, tolerance = 1e-12)
  expect_identical(a$nuisance_registry$cp04_aux_index, b$nuisance_registry$cp04_aux_index)
  expect_identical(a$reduction_summary$kept_indices, b$reduction_summary$kept_indices)
  for (ag in names(a$aggregations)) for (nm in c("estimates", "se", "if_full", "overall"))
    expect_equal(a$aggregations[[ag]][[nm]], b$aggregations[[ag]][[nm]], tolerance = 1e-12)
}

test_that("multi-cohort NT/NYT native full IFs include all pooled regression blocks and both cohort roles", {
  for (cg in c("nevertreated", "notyettreated")) {
    fx <- cp09_fixture(cg, perm = TRUE); fit <- cp09_fit(fx); gen <- cp09_fit(fx, cp09_ns(general = TRUE))
    expect_false(is.null(fit$complete_support)); expect_null(gen$complete_support)
    expect_identical(fit$complete_support$target, c(2L, 1L))
    ng <- if (cg == "nevertreated") 2L else 3L
    ac <- fit$accounting$complete_support
    expect_equal(ac$n_treatment_logit_fits, ng)
    expect_equal(ac$n_treatment_summary_batches, ng)
    expect_equal(ac$n_treatment_logit_fits + ac$n_treatment_logit_reuses, length(fit$treatments))
    st <- fx$reader$state
    for (sp in fit$treatment_specs) {
      expected <- which(st$cohort_index %in% sp$cohorts)
      expect_identical(length(expected), sp$n_observers)
      expect_identical(sp$cohorts, as.integer(c(fit$support_specs[[sp$support_ids[1L]]]$cohort,
        fit$regression_specs[[sp$regression]]$cohorts)))
      for (k in seq_along(sp$support_ids)) {
        sid <- sp$support_ids[k]
        expect_identical(fit$support_specs[[sid]]$pair, sp$pair)
        expect_identical(fit$support_specs[[sid]]$cohort, sp$cohorts[k])
        expect_equal(fit$supports[[sid]]$n_pair, sum(st$cohort_index == sp$cohorts[k]))
        expect_equal(fit$supports[[sid]]$rate, 1)
      }
    }
    H <- cp09_phi(fit); Hd <- cp09_phi(fit, details = TRUE)
    expect_equal(H, Hd, tolerance = 1e-12); expect_equal(H, cp09_phi(gen), tolerance = 1e-12)
    tab <- cp09_api$fit_dr_baseline_table(X_or = fx$X, V = fx$V, Y = fx$Y,
      obs = matrix(1, nrow(fx$Y), ncol(fx$Y)), g = fx$u$G, w = fx$weights,
      calendar = fx$calendar, links = fx$links, outcome_recipe = fx$recipe,
      glm_control = cp09_glm, rank_tol = 1e-10)
    expect_equal(H, unname(tab$phi), tolerance = 1e-12)
    expect_equal(vapply(fit$treatments, function(f) unname(f$means["mu_dr"]), numeric(1)),
      unname(tab$estimates), tolerance = 1e-12)
    expect_gt(max(abs(tab$psi_or)), 1e-6)
    expect_gt(max(abs(do.call(rbind, lapply(fit$treatments, `[[`, "d_bar_diff")))), 1e-6)
    all_links <- seq_along(fit$links)
    # A small link batch still needs the global beta score from ALL blocks.
    sub <- c(length(all_links), 1L, 7L)
    expect_equal(cp09_phi(fit, sub), H[, sub, drop = FALSE], tolerance = 1e-12)
    rr <- c(9L, 1L, 15L, 2L, 40L)
    b <- cp09_chunk(fit, rr, sub, TRUE)
    expect_identical(b$id, st$id[rr])
    expect_identical(b$details$beta_regression_ids, seq_along(fit$regression_specs))
    expect_equal(b$details$beta_score, Reduce(`+`, b$details$beta_block_scores), tolerance = 1e-12)
    local <- unique(vapply(fit$links[sub], `[[`, integer(1), "regression"))
    omitted <- setdiff(b$details$beta_regression_ids, local)
    expect_gt(length(omitted), 0L)
    expect_gt(max(abs(Reduce(`+`, b$details$beta_block_scores[
      match(omitted, b$details$beta_regression_ids)]))), 1e-6)
    expect_equal(unname(b$details$beta), unname(tab$beta_if[rr, , drop = FALSE]), tolerance = 1e-12)
    expect_equal(unname(cp09_chunk(fit, rr, sub)$phi), H[rr, sub, drop = FALSE], tolerance = 1e-12)
    expect_equal(cp09_chunk(fit, rr, sub, scale = "core")$phi,
      st$lambda * cp09_chunk(fit, rr, sub)$phi, tolerance = 1e-12)
    if (cg == "notyettreated") {
      late <- match(26, st$cohorts)
      expect_true(late %in% fit$complete_support$target && late %in% fit$complete_support$control)
      expect_true(any(vapply(fit$regression_specs, function(s) late %in% s$cohorts, logical(1))))
    }
  }
  # OR's new scope is multi-target NT only; NYT retains its general fallback.
  for (cg in c("nevertreated", "notyettreated")) {
    fx <- cp09_fixture(cg, est = "or", perm = TRUE)
    a <- cp09_fit(fx); g <- cp09_fit(fx, cp09_ns(general = TRUE))
    expect_identical(is.null(a$complete_support), cg == "notyettreated")
    expect_equal(cp09_phi(a), cp09_phi(g), tolerance = 1e-12)
    expect_equal(cp09_phi(a), cp09_phi(a, details = TRUE), tolerance = 1e-12)
  }
})

test_that("multi-cohort full inference preserves the full Gram, active ridge, auxiliary reports and aggregations", {
  for (case in list(list(est = "dr", cg = "nevertreated"), list(est = "dr", cg = "notyettreated"),
                    list(est = "reg", cg = "nevertreated"))) for (clustered in c(TRUE, FALSE)) {
    # Preserve the original clustered full-inference comparisons separately.
    # Only the unclustered request can admit the optional store; an explicit
    # 256 also admits all links in one block under the private owner row cap.
    args <- c(case, list(condition_trigger = 0.99, ridge_relative = 1e-2,
      clustered = clustered, chunk = if (clustered) 32L else 256L))
    a <- do.call(cp09_run, args)
    g <- do.call(cp09_run, c(args, list(ns = cp09_ns(general = TRUE))))
    s <- do.call(cp09_run, c(args, list(ns = cp09_ns(stream = TRUE))))
    rec <- a$resource_accounting$point_producer$retained_moment_if
    expect_identical(a$reporting$inference$clustered, clustered)
    expect_identical(rec$status, "retained")
    expect_identical(g$resource_accounting$point_producer$retained_moment_if$status, "retained")
    expect_identical(s$resource_accounting$point_producer$retained_moment_if$status, "streamed")
    expect_identical(rec$producer_passes, 1L)
    expect_identical(rec$slice_passes, 3L)
    expect_null(a$structural_certificate)
    expect_identical(a$reduction_summary$n_links_used, a$reduction_summary$n_links_all)
    expect_identical(a$covariance_handling$handler_used, "tikhonov")
    expect_gt(a$covariance_handling$lambda, 0)
    cp09_equal(a, g); cp09_equal(a, s)
    b <- a$nuisance_registry$cp04
    H <- cp09_phi(b, links = vapply(a$point_links_all, `[[`, integer(1), "cp04_index"))
    # The GMM moment Gram is the full centered original-unit Gram. Cluster
    # cross-terms in final inference are checked by cp09_equal above.
    Hc <- sweep(H, 2L, colMeans(H), "-")
    expect_equal(unname(a$covariance_stream$V), unname(crossprod(Hc) / nrow(H)), tolerance = 1e-12)
    if (!clustered) {
      # Assert the original structural, row/link partition and charged-planner
      # admission conditions, using the pipeline's recorded live reservations.
      st <- get("cp04_reader", b$reader_ref, inherits = FALSE)$state
      ids <- vapply(a$point_links_used, `[[`, integer(1), "cp04_index")
      ra <- a$resource_accounting; ctrl <- a$controls$pipeline_control$owner_control
      expect_true(inherits(b, c("cdid_or_baseline", "cdid_dr_baseline")))
      expect_false(is.null(b$complete_support)); expect_null(b$reference_ref)
      expect_equal(st$N, st$n_A, tolerance = 0); expect_true(all(st$w == 1))
      expect_false(isTRUE(b$accounting$automatic_chunk_rows))
      expect_lte(rec$chunk_rows, ctrl$max_chunk_rows)
      # Admission is one exact ordered link block at the original rows;
      # a capacity-admitted DR block can be wider than the row cap.
      expect_equal(rec$n_rows, st$N, tolerance = 0); expect_equal(rec$n_links, length(ids), tolerance = 0)
      expect_identical(rec$scale, "core")
      expect_equal(rec$chunk_rows, min(st$N, b$control$max_chunk_rows), tolerance = 0)
      expect_equal(rec$chunk_rows, ra$point_producer$chunk_rows, tolerance = 0)
      expect_equal(rec$n_batches, ceiling(st$N / rec$chunk_rows), tolerance = 0)
      b_cost <- cp09_api$cp04_retained_if_boundary_cost(b, ids, rec$chunk_rows)
      expect_equal(rec$snapshot_metadata_entries, b_cost$entries, tolerance = 0)
      expect_equal(rec$boundary_check_work_ops, b_cost$work, tolerance = 0)
      expect_equal(rec$reserved_entries,
        cp09_api$cp04_retained_if_entries(st$N, length(ids), rec$chunk_rows) + b_cost$entries, tolerance = 0)
      expect_equal(rec$work_ops_bound,
        cp09_api$cp04_retained_if_work(st$N, length(ids)) + (32 + 8 * (rec$n_batches + 1)) * b_cost$work,
        tolerance = 0)
      retained <- ra$initial_envelope_entries + ra$d5_arrays$retained_diagnostic_entries
      # Baseline DR uses the pipeline's late GMM row planner; baseline OR's
      # original GMM preflight has retained = 0.
      gmm_retained <- if (inherits(b, "cdid_dr_baseline"))
        ra$final_diagnostic_copy_entries + ra$initial_envelope_entries + ra$gmm_consumer$total_payload_entries +
          sum(unlist(ra$d5_arrays[c("ledger_if_entries", "basis_if_entries", "ledger_enclosure_entries",
            "basis_enclosure_entries", "retained_diagnostic_entries")])) else 0
      cov <- cp09_api$cp04_stream_plan(b, ids, retained)
      cov_charged <- cp09_api$cp04_stream_plan(b, ids, retained + rec$reserved_entries)
      expect_length(cov$blocks, 1L)
      expect_identical(as.integer(cov$blocks[[1L]]), seq_along(ids))
      expect_identical(cov_charged$blocks, cov$blocks)
      for (plan in list(cov, cov_charged,
        cp09_api$cp04_preflight_links(b, ids, gmm_retained),
        cp09_api$cp04_preflight_links(b, ids, gmm_retained + rec$reserved_entries)))
        expect_equal(plan$chunk_rows, rec$chunk_rows, tolerance = 0)
      live <- ra$simultaneous_entries + ra$final_diagnostic_copy_entries + 8 * st$N *
        (length(a$nuisance_registry$cp04_aux$selections[[a$controls$base_period]]) + nrow(a$support$ledger))
      expect_equal(rec$admission_total_entries, live + rec$reserved_entries, tolerance = 0)
      expect_equal(rec$cap_entries, min(ctrl$max_retained_entries, ctrl$max_matrix_entries), tolerance = 0)
      expect_lte(rec$admission_total_entries, rec$cap_entries)
      expect_lte(ra$downstream_work_ops, ctrl$max_work_ops)
      # Optional admission only: earlier fitting/resource bounds stay intact.
      at <- do.call(cp09_run, c(args, list(ns = cp09_ns(cap = rec$admission_total_entries))))
      below <- do.call(cp09_run, c(args, list(ns = cp09_ns(cap = rec$admission_total_entries - 1))))
      at_rec <- at$resource_accounting$point_producer$retained_moment_if
      expect_identical(at_rec$status, "retained")
      expect_equal(at_rec$admission_total_entries, rec$admission_total_entries, tolerance = 0)
      expect_equal(at_rec$cap_entries, rec$admission_total_entries, tolerance = 0)
      expect_identical(below$resource_accounting$point_producer$retained_moment_if$status, "streamed")
      cp09_equal(a, at); cp09_equal(a, below)
    }
  }
})

test_that("actual public att_gt fields, pretests, ordering, weights and reference requests are preserved", {
  d <- cp09_panel()
  run <- function(d, est = "dr", cg = "notyettreated", xf = ~x + x2, weightsname = NULL,
                  fix_weights = NULL, base_period = "varying", chunk = 9L, ns = NULL) {
    fn <- if (is.null(ns)) cdid::att_gt else ns$att_gt
    fn(data = d, yname = "y", tname = "period", idname = "id", gname = "G", xformla = xf,
      est_method = est, control_group = cg, bstrap = FALSE, cband = FALSE, mode = "chained",
      weightsname = weightsname, fix_weights = fix_weights, base_period = base_period,
      control = list(chunk_rows = chunk))
  }
  a <- run(d)
  general <- run(d, ns = cp09_ns(general = TRUE))
  for (nm in c("att", "se", "V_analytical", "inffunc", "W", "Wpval"))
    expect_equal(a[[nm]], general[[nm]], tolerance = 1e-12)
  expect_true(is.matrix(a$inffunc)); expect_identical(a$chained$pretest$status, "computed")
  expect_true(is.finite(a$W)); expect_true(is.finite(a$Wpval))
  p <- run(d[order(cos(seq_len(nrow(d)) * 0.7)), ], chunk = 7L)
  expect_equal(a$att, p$att, tolerance = 1e-12); expect_equal(a$se, p$se, tolerance = 1e-12)
  expect_equal(a$inffunc, p$inffunc[match(rownames(a$inffunc), rownames(p$inffunc)), , drop = FALSE], tolerance = 1e-12)
  expect_equal(a$W, p$W, tolerance = 1e-12); expect_equal(a$Wpval, p$Wpval, tolerance = 1e-12)
  for (bp in c("varying", "universal")) {
    aa <- cp09_run(base_period = bp); gg <- cp09_run(base_period = bp, ns = cp09_ns(general = TRUE))
    cp09_equal(aa, gg)
  }
  # Constant non-unit and varying weights both admit D92 retained IF; only
  # constant weights are native complete-support eligible.
  for (w in list(rep(2.5, 53L), 1 + (seq_len(53L) %% 5) / 8)) {
    dw <- d; dw$w <- w[match(dw$id, sort(unique(dw$id)))]
    aa <- cp09_run(dw, weightsname = "w"); gg <- cp09_run(dw, weightsname = "w", ns = cp09_ns(general = TRUE))
    expect_identical(aa$resource_accounting$point_producer$retained_moment_if$status, "retained")
    expect_identical(is.null(aa$nuisance_registry$cp04$complete_support), length(unique(w)) > 1L)
    cp09_equal(aa, gg)
  }
  # Equal date counts do not prove equal observer masks.
  dm <- d[!((d$id == 1053 & d$period == 4) | (d$id == 1052 & d$period == 15)), ]
  aa <- cp09_run(dm); gg <- cp09_run(dm, ns = cp09_ns(general = TRUE))
  expect_null(aa$nuisance_registry$cp04$complete_support); cp09_equal(aa, gg)
  # Dated weight recipes use the existing reference route (public OR scope).
  dw <- d; dw$w <- 1 + (dw$id %% 3) / 4 + match(dw$period, sort(unique(dw$period))) / 10
  ref <- run(dw, est = "reg", cg = "nevertreated", xf = NULL, weightsname = "w", fix_weights = "base_period")
  expect_false(is.null(ref$chained$reference_inputs))
  expect_true(all(is.finite(ref$att[ref$chained$available])))
  ref_g <- run(dw, est = "reg", cg = "nevertreated", xf = NULL, weightsname = "w",
    fix_weights = "base_period", ns = cp09_ns(general = TRUE))
  expect_equal(ref$inffunc, ref_g$inffunc, tolerance = 1e-12)
  expect_equal(ref$att, ref_g$att, tolerance = 1e-12)
  expect_error(run(d, est = "reg", xf = ~0), class = "cdid_chained_request_error")
})

test_that("group replay identity, buffer reservations, aliases and no-intercept fallback remain enforced", {
  fx <- cp09_fixture(); fit <- cp09_fit(fx); ll <- seq_along(fit$links)
  plan <- cp09_api$dr_baseline_native_plan(fit, ll, 9L)
  expect_true(cp09_api$ipw_owner_seal_ok(plan, "dr_baseline_native_plan"))
  unsealed <- plan; unsealed$ref <- NULL
  expect_equal(plan$entries, unname(cp09_api$cp04_payload(unsealed)["entries"]))
  pre <- cp09_api$dr_baseline_if_budget(fit, 9L, ll)
  expect_lte(2 * cp09_api$dr_baseline_group_plan_bound(fit, ll, pre, pre$tids), pre$scratch)
  expect_equal(pre$complete$treatment_replays_per_chunk, 3L)
  expect_equal(cp09_chunk(fit, 1:9, ll, plan = plan)$phi, cp09_chunk(fit, 1:9, ll)$phi, tolerance = 1e-12)
  broken <- fit; broken$treatments[[1L]]$coefficients[1L] <- broken$treatments[[1L]]$coefficients[1L] + 0.1
  expect_error(cp09_chunk(broken, 1:9, ll), class = "cdid_dr_baseline_integrity_error")
  # Resealing a within-group coefficient discrepancy reaches the exact replay
  # equality guard, independently of the outer owner's immutability guard.
  keys <- vapply(fit$treatment_specs, function(s) paste(s$cohorts, collapse = ":"), character(1))
  idx <- which(keys == keys[1L])[2L]
  bad <- fit; bad$ref <- NULL
  bad$treatments[[idx]]$coefficients <- bad$treatments[[idx]]$coefficients + 0.1
  bad <- cp09_api$ipw_owner_seal(bad, "dr_baseline_backend")
  expect_error(cp09_api$dr_baseline_native_plan(bad, ll, 9L), class = "cdid_dr_baseline_identity_error")
  expect_error(cp09_chunk(bad, seq_len(9L), ll, TRUE), class = "cdid_dr_baseline_identity_error")
  mutated <- plan; mutated$groups[[1L]]$cohorts <- rev(mutated$groups[[1L]]$cohorts)
  expect_error(cp09_chunk(fit, 1:9, ll, plan = mutated), class = "cdid_dr_baseline_integrity_error")
  # Budget refusal occurs at the envelope boundary, before descriptors/read.
  tight <- fit; tight$ref <- NULL
  tight$control$max_matrix_entries <- fit$accounting$owned_entries + fit$accounting$borrowed_entries + pre$scratch - 1
  tight <- cp09_api$ipw_owner_seal(tight, "dr_baseline_backend")
  expect_error(cp09_api$dr_baseline_native_plan(tight, ll, 9L), class = "cdid_ipw_owner_bound_error")
  work_tight <- fit; work_tight$ref <- NULL
  work_tight$control$max_work_ops <- fit$accounting$work_ops_bound + pre$work - 1
  work_tight <- cp09_api$ipw_owner_seal(work_tight, "dr_baseline_backend")
  expect_error(cp09_api$dr_baseline_native_plan(work_tight, ll, 9L), class = "cdid_ipw_owner_bound_error")
  for (args in list(list(alias = TRUE), list(no_intercept = TRUE))) {
    ff <- do.call(cp09_fixture, args); aa <- cp09_fit(ff); gg <- cp09_fit(ff, cp09_ns(general = TRUE))
    expect_identical(is.null(aa$complete_support), isTRUE(args$no_intercept))
    if (isTRUE(args$alias)) expect_true("x_alias" %in% aa$shared$basis$aliased_columns)
    expect_equal(cp09_phi(aa), cp09_phi(gg), tolerance = 1e-12)
  }
})
