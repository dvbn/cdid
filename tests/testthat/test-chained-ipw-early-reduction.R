# CP08 P3: outcome-free structural reduction before point fits, licensed only for
# built-in baseline IPW under optimal weighting (Prop. 5 / Cor. 6). Strict
# actual-installed namespace tests; the container source audit uses a separate
# disclosed adapter.
# >>> installed-namespace binding
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api), environmentIsLocked(api))
api_binding <- list(mode = "installed_package", env = api)
# <<< installed-namespace binding

# ---- fixtures and helpers ----------------------------------------------------
er_balanced48 <- function() {
  set.seed(20261016L)
  u_g <- rep(c(3, 4, Inf), each = 16L)
  u_x <- rep(c(-1.0, 1.0, -0.5, 0.5), times = 12L)
  u_w <- rep(c(1.0, 1.5, 0.8, 1.2), times = 12L)
  u_cl <- rep(seq_len(12L), each = 4L)
  d <- expand.grid(t = 1:4, id = 1:48, KEEP.OUT.ATTRS = FALSE)
  d$g <- u_g[d$id]; d$x <- u_x[d$id]; d$weight <- u_w[d$id]; d$cluster <- u_cl[d$id]
  eff <- ifelse(d$g == 3 & d$t >= 3, 2.0 + 0.5 * (d$t - 3), ifelse(d$g == 4 & d$t >= 4, 1.5, 0))
  d$y <- 0.25 * d$id + 0.5 * d$t + 0.3 * d$x * d$t + eff + 0.2 * stats::rnorm(nrow(d))
  d[, c("id", "t", "g", "x", "weight", "cluster", "y")]
}
er_frame48 <- list(id = c(1:48, 49:52), g = c(rep(c(3, 4, Inf), each = 16L), 3, 4, Inf, Inf),
  w = c(rep(c(1.0, 1.5, 0.8, 1.2), times = 12L), 1, 1, 1, 1), cl = c(rep(seq_len(12L), each = 4L), rep(13L, 4L)),
  Z = stats::model.matrix(~ x, data = data.frame(x = c(rep(c(-1.0, 1.0, -0.5, 0.5), times = 12L), -1, 1, -0.5, 0.5))))

er_control <- function(n_frame, n_rows, n_dates, n_cohorts, n_targets, n_cov, n_clusters, tuning = NULL,
                       weighting = "optimal") {
  dims <- api$chained_public_dimensions_from_counts(n_frame = n_frame, n_panel_units = min(n_frame, 48L),
    n_rows = n_rows, n_dates = n_dates, n_cohorts = n_cohorts, n_targets = n_targets, n_covariates = n_cov,
    biters = 99, n_clusters = n_clusters)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(tuning, gmm_weighting = weighting),
    gmm_weighting = weighting)$control
}

er_run48 <- function(df, ctrl, stage, weighting = "optimal", reduce = TRUE, custom = NULL) {
  api$run_chained_ipw_pipeline(data = df, idname = "id", tname = "t", yname = "y", gname = "g",
    xformla = ~ x, weightsname = "weight", cluster_vec = er_frame48$cl, frame_id = er_frame48$id,
    frame_Z = er_frame48$Z, frame_weights = er_frame48$w, frame_g = er_frame48$g, allow_frame_only = TRUE,
    control = ctrl, control_group = "notyettreated", base_period = "varying", weighting = weighting,
    reduce_structural_redundancy = reduce, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = c("simple", "group", "dynamic", "calendar"), cband = TRUE, biters = 99L, seed = 7L,
    custom_est_method = custom, structural_reduction_stage = stage)
}

# Strip environment identities (fresh owners/views per run) for structural comparison.
er_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (inherits(x, "cdid_ipw_structural_certificate")) {
    x <- unclass(x)
    return(lapply(x[setdiff(names(x), c("link_kind", "object_refs", "ref"))], er_strip))
  }
  if (is.list(x)) {
    a <- attributes(x)
    y <- lapply(x, er_strip)
    attributes(y) <- a
    return(y)
  }
  x
}
er_norm_cert_refs <- function(cert, n_cat) {
  if (!inherits(cert, "cdid_ipw_structural_certificate")) return(cert)
  u <- unclass(cert)
  u$nuisance_refs <- u$nuisance_refs[seq_len(n_cat)]
  class(u) <- class(cert)
  u
}
er_expect_same <- function(a, b, fields) {
  raw_early_vs_late <- !is.null(a$nuisance_registry$raw_ipw) && !is.null(b$nuisance_registry$raw_ipw) &&
    identical(a$reduction_summary$reduction_stage, "outcome_free_before_point_fits") &&
    identical(b$reduction_summary$reduction_stage, "after_point_fits")
  for (f in fields) {
    af <- a[[f]]
    bf <- b[[f]]
    if (raw_early_vs_late && identical(f, "structural_certificate")) {
      n_cat <- length(a$nuisance_registry$nuisance_keys)
      n_aux_extra <- length(b$nuisance_registry$raw_aux_nuisance_views) - n_cat
      expect_identical(names(af$nuisance_refs), a$nuisance_registry$nuisance_keys)
      expect_identical(names(bf$nuisance_refs), c(b$nuisance_registry$nuisance_keys, rep("", n_aux_extra)))
      af <- er_norm_cert_refs(af, n_cat)
      bf <- er_norm_cert_refs(bf, n_cat)
    } else if (raw_early_vs_late && identical(f, "gmm_fit")) {
      n_cat <- length(a$nuisance_registry$nuisance_keys)
      af$provenance$reduction$class_certificate <- er_norm_cert_refs(af$provenance$reduction$class_certificate, n_cat)
      bf$provenance$reduction$class_certificate <- er_norm_cert_refs(bf$provenance$reduction$class_certificate, n_cat)
      if (!is.null(af$provenance$weight_covariance$reduction$class_certificate)) {
        af$provenance$weight_covariance$reduction$class_certificate <- er_norm_cert_refs(
          af$provenance$weight_covariance$reduction$class_certificate, n_cat)
        bf$provenance$weight_covariance$reduction$class_certificate <- er_norm_cert_refs(
          bf$provenance$weight_covariance$reduction$class_certificate, n_cat)
      }
    } else if (raw_early_vs_late && identical(f, "point_links_used")) {
      kept <- a$reduction_summary$kept_indices
      expect_identical(unname(vapply(af, `[[`, integer(1L), "cp04_index")), match(kept, sort(kept)))
      expect_identical(unname(vapply(bf, `[[`, integer(1L), "cp04_index")), as.integer(kept))
      af <- lapply(af, function(pl) { pl$cp04_index <- NULL; pl })
      bf <- lapply(bf, function(pl) { pl$cp04_index <- NULL; pl })
    } else if (raw_early_vs_late && identical(f, "auxiliary_reporting") && length(af$point_links) > 0L) {
      n_drop <- b$reduction_summary$n_rows_fitted - a$reduction_summary$n_rows_fitted
      for (j in seq_along(af$point_links)) {
        if (!is.null(af$point_links[[j]])) {
          expect_identical(bf$point_links[[j]]$cp04_index - af$point_links[[j]]$cp04_index, as.integer(n_drop))
          af$point_links[[j]]$cp04_index <- bf$point_links[[j]]$cp04_index
        }
      }
      owned_diff <- b$nuisance_registry$raw_ipw$accounting$owned_entries -
        a$nuisance_registry$raw_ipw$accounting$owned_entries
      expect_identical(bf$preflight$producer_workspace_entries - af$preflight$producer_workspace_entries, owned_diff)
      expect_identical(bf$preflight$peak_concurrent_entries - af$preflight$peak_concurrent_entries, owned_diff)
      af$preflight$producer_workspace_entries <- bf$preflight$producer_workspace_entries
      af$preflight$peak_concurrent_entries <- bf$preflight$peak_concurrent_entries
    }
    expect_identical(er_strip(af), er_strip(bf), info = f)
  }
}
er_complete_fields <- c("support_status", "support", "reporting", "basis", "aggregations", "covariance_stream",
  "gmm_fit", "covariance_handling", "auxiliary_reporting", "full_frame_transport", "design_all", "design_used",
  "point_links_used", "structural_certificate", "controls", "metadata", "panel_owner", "frame_owner", "selection")
er_new_summary_fields <- c("max_moment_expansion_diff", "reduction_stage", "n_rows_fitted", "expansion_self_test")

# Count calls of namespace functions during `run`, restoring the bindings afterwards.
er_count_calls <- function(names, run) {
  counts <- stats::setNames(integer(length(names)), names)
  originals <- lapply(names, function(nm) api[[nm]])
  names(originals) <- names
  put <- function(nm, fn) {
    locked <- bindingIsLocked(nm, api)
    if (locked) unlockBinding(nm, api)
    assign(nm, fn, envir = api)
    if (locked) lockBinding(nm, api)
  }
  on.exit(for (nm in names) put(nm, originals[[nm]]), add = TRUE)
  for (nm in names) local({
    nm0 <- nm; f0 <- originals[[nm]]
    put(nm0, function(...) { counts[[nm0]] <<- counts[[nm0]] + 1L; f0(...) })
  })
  value <- run()
  list(value = value, counts = counts)
}

# Independent full-dimensional expansion matrix in the retained (design) column order.
er_expansion_matrix <- function(res) {
  cert <- res$structural_certificate
  kept <- res$reduction_summary$kept_indices
  E <- matrix(0, cert$n_links, length(kept))
  E[cbind(kept, seq_along(kept))] <- 1
  for (k in seq_along(cert$E$i)) {
    j <- match(cert$E$j[k], kept)
    E[cert$E$i[k], j] <- E[cert$E$i[k], j] + cert$E$x[k]
  }
  E
}

df48 <- er_balanced48()
ctrl48 <- er_control(52L, nrow(df48), 4L, 3L, 2L, 2L, 13L)
ridge48 <- er_control(52L, nrow(df48), 4L, 3L, 2L, 2L, 13L, tuning = list(condition_trigger = 0.99, ridge_relative = 0.01))

# Each tuning runs once per stage; the results are shared by the tests below.
er_counted <- c("fit_ipw_point_link_view", "ipw_raw_make_point_link_record", "ipw_owner_read_y_positions")
runs48 <- lapply(list(plain = ctrl48, ridge = ridge48), function(ctrl) list(ctrl = ctrl,
  early = er_count_calls(er_counted, function() er_run48(df48, ctrl, "outcome_free")),
  late = er_count_calls(er_counted, function() er_run48(df48, ctrl, "post_point_fit"))))

test_that("outcome-free reduction reproduces the accepted post-fit path's complete outputs and fits only retained rows", {
  for (run in runs48) {
    early <- run$early; late <- run$late
    e <- early$value; l <- late$value
    rs <- e$reduction_summary
    expect_true(rs$reduction_active)
    expect_identical(rs$reduction_stage, "outcome_free_before_point_fits")
    expect_identical(l$reduction_summary$reduction_stage, "after_point_fits")
    expect_identical(rs$n_rows_fitted, rs$n_links_used)
    expect_identical(l$reduction_summary$n_rows_fitted, rs$n_links_all)
    er_expect_same(e, l, er_complete_fields)
    keep <- setdiff(names(l$reduction_summary), er_new_summary_fields)
    expect_identical(er_strip(e$reduction_summary[keep]), er_strip(l$reduction_summary[keep]))
    reg_shared <- c("nuisance_views", "fitted_views", "nuisance_keys", "operator_ids", "row_specs",
      "link_views", "link_to_nuisance_idx", "link_nuisance_keys", "mask_splits", "mask_certificate",
      "kept_idx", "dropped_idx", "raw_aux_errors", "raw_aux_fit_errors", "raw_aux_counts",
      "raw_aux_fitted_views", "raw_aux_nuisance_views")
    expect_identical(er_strip(e$nuisance_registry[intersect(reg_shared, names(e$nuisance_registry))]),
      er_strip(l$nuisance_registry[intersect(reg_shared, names(l$nuisance_registry))]))
    cw_shared <- setdiff(names(e$nuisance_registry$construction_work), "work_ops_charged")
    expect_identical(e$nuisance_registry$construction_work[cw_shared], l$nuisance_registry$construction_work[cw_shared])
    expect_identical(unclass(e$structural_certificate)$link_kind, rep("shared_view_declaration", rs$n_links_all))
    # Dropped rows: outcome-free identity records sharing the right fitted nuisance; never fitted.
    drop <- rs$dropped_indices
    expect_true(length(drop) > 0L)
    for (m in drop) {
      rec <- e$point_links_all[[m]]
      expect_s3_class(rec, "cdid_ipw_reduced_link_record")
      expect_false(inherits(rec, "cdid_ipw_point_link_view"))
      expect_identical(rec$link_id, l$point_links_all[[m]]$link_id)
      expect_true(is.null(rec$estimate))
      expect_false(rec$point_fit_performed)
    }
    expect_identical(names(e$point_links_all), names(l$point_links_all))
    # Avoided work, by exact instrumented counts: one point fit per dropped row,
    # fewer outcome-endpoint reader calls; reporting fits are common to both paths.
    pt_fit_diff <- (late$counts[["fit_ipw_point_link_view"]] - early$counts[["fit_ipw_point_link_view"]]) +
      (late$counts[["ipw_raw_make_point_link_record"]] - early$counts[["ipw_raw_make_point_link_record"]])
    expect_identical(pt_fit_diff, length(drop))
    expect_lt(early$counts[["ipw_owner_read_y_positions"]], late$counts[["ipw_owner_read_y_positions"]])
    # Self-test values: the accepted path fits dropped rows; the early path cannot.
    expect_lt(l$reduction_summary$max_moment_expansion_diff, 1e-12)
    expect_true(is.na(rs$max_moment_expansion_diff))
    expect_identical(rs$expanded_m, l$reduction_summary$expanded_m)
  }
})

test_that("independent full-criterion references: moments, complete IF, covariance, Tikhonov metric and stationarity", {
  for (run in runs48) {
    ctrl <- run$ctrl
    e <- run$early$value
    l <- run$late$value
    E <- er_expansion_matrix(e)
    M <- nrow(E)
    pl_all <- l$point_links_all # every original row independently point-fitted
    m_full <- vapply(pl_all, `[[`, numeric(1L), "estimate")
    expect_lt(max(abs(m_full - drop(E %*% e$covariance_stream$mhat))), 1e-12)
    # Complete original-unit IF: all M replayed columns equal E times the retained columns.
    Phi <- matrix(NA_real_, l$covariance_stream$n_a_units, M)
    api$replay_ipw_point_if_view(pl_all, ctrl$point_control, callback_schema = TRUE,
      callback = function(b) Phi[b$a_index, ] <<- b$phi)
    Phi_R <- Phi[, e$reduction_summary$kept_indices, drop = FALSE]
    expect_lt(max(abs(Phi - Phi_R %*% t(E))) / max(1, max(abs(Phi))), 1e-12)
    # Full covariance of all original rows equals E V_R E' (accepted covariance stream).
    V_full <- api$stream_ipw_point_covariance_view(pl_all, ctrl$point_control)$covariance
    V_R <- unname(e$covariance_stream$V)
    expect_lt(max(abs(V_full - E %*% V_R %*% t(E))) / max(abs(V_full)), 1e-12)
    ch <- e$gmm_fit$covariance_handling
    expect_identical(ch$dimensions$M_full, as.integer(M))
    expect_identical(ch$dimensions$L_retained, as.integer(e$reduction_summary$n_links_used))
    expect_lt(abs(ch$trace_B - sum(diag(V_full))) / sum(diag(V_full)), 1e-12)
    if (identical(ch$handler_used, "tikhonov")) {
      # Ridge active: lambda in the original metric, and the reduced solution is
      # stationary for the FULL M-dimensional Tikhonov criterion built from the
      # independently fitted rows, full covariance and full incidence.
      lam <- ch$ridge_relative * sum(diag(V_full)) / M
      expect_lt(abs(ch$lambda - lam) / lam, 1e-12)
      inc <- e$design_all$incidence
      W <- matrix(0, inc$n_links, inc$n_parameters, dimnames = list(NULL, inc$param_ids))
      k <- which(inc$start > 0L); W[cbind(k, inc$start[k])] <- -1
      k <- which(inc$end > 0L); W[cbind(k, inc$end[k])] <- 1
      theta <- e$gmm_fit$internal_gauge$theta_gauge_coords
      expect_true(all(names(theta) %in% colnames(W)))
      ord <- e$design_all$links$catalog_row
      W_G <- W[, names(theta), drop = FALSE]
      r <- m_full[ord] - drop(W_G %*% theta)
      y <- solve(V_full[ord, ord] + diag(lam, M), r)
      grad <- drop(crossprod(W_G, y))
      expect_lt(sqrt(sum(grad^2)) / (sqrt(sum(W_G^2)) * sqrt(sum(y^2))), 1e-10)
    } else {
      expect_identical(ch$trigger_reason, "well_conditioned_pd")
    }
  }
  expect_identical(runs48$ridge$early$value$gmm_fit$covariance_handling$handler_used, "tikhonov")
})

test_that("PT-all with irregular dates and unequal observation masks: early equals the accepted path", {
  set.seed(31L)
  dates <- c(1, 2, 4.5, 5, 7)
  d <- expand.grid(t = dates, id = 1:45, KEEP.OUT.ATTRS = FALSE)
  d$g <- rep(c(4.5, 7, Inf), each = 15L)[d$id]; d$x <- rep(c(-1, 0, 1), times = 15L)[d$id]; d$weight <- 1
  d$y <- 0.1 * d$id + 0.3 * d$t + 0.2 * d$x * d$t + ifelse(d$t >= d$g, 1, 0) + 0.3 * stats::rnorm(nrow(d))
  ctrl <- er_control(45L, nrow(d), 5L, 3L, 2L, 2L, 45L)
  run <- function(df, stage, pt) api$run_chained_ipw_pipeline(data = df, idname = "id", tname = "t", yname = "y",
    gname = "g", xformla = ~ x, control = ctrl, control_group = "notyettreated", base_period = "universal",
    weighting = "optimal", run_numerical_transfer = FALSE, root_witness = NULL, aggregations = c("simple", "dynamic"),
    cband = FALSE, pt_assumption = pt, pt_start = if (identical(pt, "all")) 1 else NULL, structural_reduction_stage = stage)
  masked <- d[!(d$id %% 3L == 0L & d$t == dates[1L + (d$id %% 5L)]), ]
  for (case in list(list(df = d, pt = "all"), list(df = masked, pt = "post"), list(df = masked, pt = "all"))) {
    e <- run(case$df, "outcome_free", case$pt)
    l <- run(case$df, "post_point_fit", case$pt)
    er_expect_same(e, l, er_complete_fields)
    expect_identical(e$reduction_summary$n_rows_fitted, e$reduction_summary$n_links_used)
    expect_identical(e$reduction_summary$kept_indices, l$reduction_summary$kept_indices)
  }
})

test_that("identity weighting, unreduced requests and custom callbacks keep their accepted paths", {
  id48 <- er_control(52L, nrow(df48), 4L, 3L, 2L, 2L, 13L, weighting = "identity")
  for (args in list(list(weighting = "identity", reduce = TRUE, ctrl = id48), list(weighting = "optimal", reduce = FALSE, ctrl = ctrl48))) {
    a <- er_run48(df48, args$ctrl, "outcome_free", weighting = args$weighting, reduce = args$reduce)
    b <- er_run48(df48, args$ctrl, "post_point_fit", weighting = args$weighting, reduce = args$reduce)
    er_expect_same(a, b, c(er_complete_fields, "nuisance_registry", "point_links_all", "resource_accounting"))
    expect_false(a$reduction_summary$reduction_active)
    expect_identical(a$reduction_summary$n_rows_fitted, a$reduction_summary$n_links_all)
    expect_false(identical(a$reduction_summary$reduction_stage, "outcome_free_before_point_fits"))
  }
  custom_ipw <- function(y1, y0, D, covariates, i.weights, inffunc = TRUE, ...) {
    dy <- y1 - y0; wt <- i.weights * D; wc <- i.weights * (1 - D)
    mt <- sum(wt * dy) / sum(wt); mc <- sum(wc * dy) / sum(wc)
    list(ATT = mt - mc, att.inf.func = length(dy) * (wt * (dy - mt) / sum(wt) - wc * (dy - mc) / sum(wc)))
  }
  a <- er_run48(df48, ctrl48, "outcome_free", custom = custom_ipw)
  b <- er_run48(df48, ctrl48, "post_point_fit", custom = custom_ipw)
  er_expect_same(a, b, c(er_complete_fields, "point_links_all", "resource_accounting"))
  expect_identical(a$reduction_summary$reduction_stage, "after_point_fits")
  expect_identical(unique(unclass(a$structural_certificate)$link_kind), "custom_original_ipw")
})

test_that("CP04 OR/DR keep their accepted path and built-in IPW takes the outcome-free route through public-route preparation", {
  # Public preparation up to the pipeline call (no did object is built here);
  # the public never-treated code is 0.
  dpub <- df48
  dpub$g[!is.finite(dpub$g)] <- 0
  public_route <- function(est, stage) {
    gate <- api$chained_request_scope_gate(mode = "chained", panel = TRUE, allow_unbalanced_panel = NULL,
      est_method = est, selection = "none", selection_xformla = NULL, pt_assumption = "post", pt_start = NULL,
      gmm_weighting = "optimal", gmm_covariance = "tikhonov", fix_weights = NULL, compute_inffunc = TRUE, pl = FALSE,
      cores = 1, faster_mode = TRUE, print_details = FALSE, control_group = "notyettreated", base_period = "varying",
      anticipation = 0, alp = 0.05, bstrap = FALSE, cband = FALSE, biters = 99, dots = list(), entry = "public")
    resolved <- api$chained_public_resolve_control(NULL, gmm_weighting = gate$gmm_weighting_requested,
      gmm_covariance = gate$gmm_covariance)
    prep <- api$chained_request_prepare_inputs(data = dpub, frame = NULL, yname = "y", tname = "t", idname = "id",
      gname = "g", xformla = ~ x, weightsname = "weight", clustervars = NULL, anticipation = 0, target_cohorts = NULL,
      selection = gate$selection, outcome_regression = identical(est, "reg"), or_resource = resolved$profile,
      doubly_robust = identical(est, "dr"), selection_est_method = gate$selection_est_method,
      imputation_pooling = gate$imputation_pooling, imputation_est_method = NULL)
    if (est %in% c("reg", "dr") && !is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
    produced <- api$chained_public_pipeline_control(api$chained_public_dimensions(prep, biters = 99, clustered = FALSE),
      resolved, gmm_weighting = gate$gmm_weighting_requested, gmm_covariance = gate$gmm_covariance)
    api$run_chained_ipw_pipeline(data = prep$data_internal, idname = "id", tname = "t", yname = "y", gname = "g",
      control = produced$control, xformla = prep$xformla, weightsname = NULL, clustername = NULL,
      cluster_vec = prep$cluster_full, target_cohorts = prep$target_cohorts, anticipation = 0L,
      control_group = gate$control_group, base_period = "varying", weighting = gate$gmm_weighting_requested,
      pt_assumption = gate$pt_assumption, pt_start = gate$pt_start, run_numerical_transfer = FALSE, root_witness = NULL,
      aggregations = character(0L), cband = FALSE, alp = 0.05, biters = 99L, btype = "rademacher", seed = NULL,
      fixed_signs = NULL, cp04 = prep$cp04, frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
      frame_weights = prep$frame_w, frame_g = prep$frame_g_arg, allow_frame_only = prep$n_frame_only > 0L,
      na.rm = FALSE, bstrap = FALSE, caller_rng = FALSE, custom_est_method = NULL, custom_dots = list(),
      gmm_covariance = gate$gmm_covariance, structural_reduction_stage = stage)
  }
  for (est in c("reg", "dr", "ipw")) {
    a <- public_route(est, "outcome_free")
    b <- public_route(est, "post_point_fit")
    er_expect_same(a, b, er_complete_fields)
    if (est == "ipw") {
      expect_identical(a$reduction_summary$reduction_stage, "outcome_free_before_point_fits")
    } else {
      expect_false(a$reduction_summary$reduction_active)
      expect_identical(a$reduction_summary$reduction_stage, "none")
      expect_identical(a$reduction_summary$n_rows_fitted, a$reduction_summary$n_links_all)
      er_expect_same(a, b, c("point_links_all", "resource_accounting"))
    }
  }
})

test_that("outcome-free certificate equals the point-link certificate; tampered, unverified, fitted or oversized inputs refuse", {
  l <- runs48$plain$late$value
  reg <- api$build_pipeline_nuisance_and_link_views(metadata_obj = l$metadata, frame_owner = l$frame_owner,
    owner_control = ctrl48$owner_control, fitted_control = ctrl48$fitted_control, point_control = ctrl48$point_control,
    glm_control = ctrl48$glm_control, fit_point_links = FALSE)
  expect_null(reg$point_links)
  M <- length(reg$link_views)
  dli <- match(seq_len(M), l$design_all$links$catalog_row)
  cv <- api$certify_ipw_structural_classes_from_views(reg, l$design_all, dli, ctrl48$class_bounds, ctrl48$owner_control)
  cp <- api$certify_ipw_structural_classes_view(l$point_links_all, l$design_all, dli, "optimal", ctrl48$class_bounds, ctrl48$owner_control)
  expect_identical(er_strip(cv), er_strip(cp))
  expect_identical(unclass(cp)$link_kind, rep("shared_view_fit", M))
  expect_identical(unclass(cv)$object_refs, unname(lapply(reg$link_views, `[[`, "ref")))
  bad <- reg
  lv <- bad$link_views[[2L]]; lv$functional_key <- "tampered"; bad$link_views[[2L]] <- lv
  expect_error(api$certify_ipw_structural_classes_from_views(bad, l$design_all, dli, ctrl48$class_bounds, ctrl48$owner_control),
    class = "cdid_ipw_owner_integrity_error")
  fitted <- reg; fitted$point_links <- l$point_links_all
  expect_error(api$certify_ipw_structural_classes_from_views(fitted, l$design_all, dli, ctrl48$class_bounds, ctrl48$owner_control),
    class = "cdid_ipw_certificate_error")
  k <- reg$link_to_nuisance_idx[2L]
  nv <- reg$nuisance_views[[k]]
  unverified <- reg
  unverified$link_views[[2L]] <- api$construct_ipw_link_view(l$frame_owner, reg$link_views[[2L]]$pair$s,
    reg$link_views[[2L]]$pair$t, nuisance_view = nv, control = ctrl48$owner_control,
    link_id = reg$link_views[[2L]]$link_id, functional_key = reg$link_views[[2L]]$functional_key,
    require_support = TRUE, verify_nuisance_mask = FALSE)
  expect_error(api$certify_ipw_structural_classes_from_views(unverified, l$design_all, dli, ctrl48$class_bounds, ctrl48$owner_control),
    class = "cdid_ipw_owner_mask_unverified_error")
  small <- ctrl48$class_bounds; small$max_links <- M - 1
  expect_error(api$certify_ipw_structural_classes_from_views(reg, l$design_all, dli, small, ctrl48$owner_control),
    class = "cdid_ipw_bound_error")
  # Covariates with missing values are refused before any view, fit or certificate, on both stages.
  dna <- df48; dna$x[dna$id == 5L] <- NA_real_
  for (stage in c("outcome_free", "post_point_fit")) {
    expect_error(api$run_chained_ipw_pipeline(data = dna, idname = "id", tname = "t", yname = "y", gname = "g",
      xformla = ~ x, control = ctrl48, control_group = "notyettreated", weighting = "optimal",
      run_numerical_transfer = FALSE, root_witness = NULL, aggregations = character(0L), structural_reduction_stage = stage),
      class = "cdid_ipw_pipeline_input_error")
  }
  expect_error(er_run48(df48, ctrl48, "early"), class = "cdid_ipw_pipeline_input_error")
})

test_that("deterministic envelope arithmetic: covariance/GMM terms move from M_full to the retained basis, original rows stay charged", {
  e <- runs48$plain$early$value
  l <- runs48$plain$late$value
  rs <- e$reduction_summary
  # CP08 scaling action 3 (contract change): the exact union of both auxiliary
  # reporting plans replaces the former 2 J T = 16 bound, and the auxiliary/report
  # buffers are reserved separately from the retained covariance, identically on
  # both paths, so only the covariance/GMM/replay terms differ here. Independent
  # union for cohorts {3, 4, Inf}, dates 1:4, targets {3, 4}, not-yet-treated:
  # (3,1,2,{4,Inf}), (4,1,2,{3,Inf}), (4,2,3,{Inf}), (4,1,3,{Inf}); the varying
  # family reports three pretrend columns.
  for (r in list(e, l)) {
    ap <- r$resource_accounting$auxiliary_plan
    expect_identical(c(ap$n_descriptors_union, ap$n_reported_pretrend, ap$replaced_bound), c(4, 3, 16))
  }
  expect_identical(e$resource_accounting$auxiliary_plan, l$resource_accounting$auxiliary_plan)
  Q <- 0
  Mb <- rs$n_links_all + Q
  Mu <- rs$n_links_used + Q
  cert <- e$structural_certificate
  cert_entries <- e$resource_accounting$early_structural$additional_entries
  chunk <- ctrl48$owner_control$max_chunk_rows
  p <- 2
  raw_owned_diff <- if (is.null(e$nuisance_registry$raw_ipw)) 0 else
    l$nuisance_registry$raw_ipw$accounting$owned_entries - e$nuisance_registry$raw_ipw$accounting$owned_entries
  raw_work_diff <- if (is.null(e$nuisance_registry$raw_ipw)) 0 else
    l$nuisance_registry$construction_work$work_ops_charged - e$nuisance_registry$construction_work$work_ops_charged
  expect_identical(l$resource_accounting$initial_envelope_entries - e$resource_accounting$initial_envelope_entries,
    16 * (Mb^2 - Mu^2) + 32 * chunk * (Mb - Mu - p - 1) - cert_entries + raw_owned_diff)
  expect_identical(l$resource_accounting$downstream_work_ops - e$resource_accounting$downstream_work_ops,
    32 * 52 * (Mb^2 - Mu^2) - e$resource_accounting$early_structural$work_ops_bound + raw_work_diff)
  expect_identical(e$resource_accounting$gmm_consumer, l$resource_accounting$gmm_consumer)
})


test_that("early certificate reserves actual fitted views and all expansion copies before construction", {
  e <- runs48$plain$early$value
  reg <- e$nuisance_registry
  reg$point_links <- reg$point_links_all <- NULL
  pf <- api$preflight_ipw_early_structural_reduction(reg, e$design_all, ctrl48$class_bounds, ctrl48$owner_control)
  es <- e$resource_accounting$early_structural
  pf_shared <- setdiff(names(pf), c("original_link_view_entries", "additional_entries"))
  expect_identical(pf[pf_shared], es[pf_shared])
  expect_identical(es$original_link_view_entries, 0)
  expect_identical(es$additional_entries, pf$certificate_retained_entries)
  actual_fit <- sum(vapply(reg$fitted_views, function(fv) fv$accounting$total_retained_entries, numeric(1)))
  actual_view <- sum(vapply(reg$link_views, function(lv) lv$preflight$total_retained_entries, numeric(1)))
  expect_identical(pf$fitted_nuisance_entries, actual_fit)
  expect_identical(pf$original_link_view_entries, actual_view)
  expect_gte(pf$certificate_retained_entries, 5 * length(e$structural_certificate$E$i))
  expect_gte(pf$E_entries_upper, length(e$structural_certificate$E$i))
  expect_identical(pf$p1_registry_entries, 16 * (length(reg$link_views) + 1) + length(e$panel_owner$cohorts))
  # Refuse the composed cap before entering the accepted certificate core, even
  # though each resource control individually admits this small registry.
  dli <- match(seq_along(reg$link_views), e$design_all$links$catalog_row)
  core <- api$ipw_view_structural_classes_core
  calls <- 0L
  put <- function(fn) {
    locked <- bindingIsLocked("ipw_view_structural_classes_core", api)
    if (locked) unlockBinding("ipw_view_structural_classes_core", api)
    assign("ipw_view_structural_classes_core", fn, envir = api)
    if (locked) lockBinding("ipw_view_structural_classes_core", api)
  }
  on.exit(put(core), add = TRUE)
  put(function(...) { calls <<- calls + 1L; core(...) })
  cap <- min(ctrl48$owner_control$max_retained_entries, ctrl48$owner_control$max_matrix_entries)
  expect_error(api$certify_ipw_structural_classes_from_views(reg, e$design_all, dli,
    ctrl48$class_bounds, ctrl48$owner_control, retained_pipeline_entries = cap - pf$additional_entries + 1),
    class = "cdid_ipw_owner_bound_error")
  expect_identical(calls, 0L)
  expect_error(api$certify_ipw_structural_classes_from_views(reg, e$design_all, dli,
    ctrl48$class_bounds, ctrl48$owner_control,
    external_work_ops = ctrl48$owner_control$max_work_ops - pf$work_ops_bound + 1),
    class = "cdid_ipw_owner_bound_error")
  expect_identical(calls, 0L)
  invisible(api$certify_ipw_structural_classes_from_views(reg, e$design_all, dli,
    ctrl48$class_bounds, ctrl48$owner_control, retained_pipeline_entries = cap - pf$additional_entries,
    external_work_ops = ctrl48$owner_control$max_work_ops - pf$work_ops_bound))
  expect_identical(calls, 1L)
})

test_that("reduced records copy exactly their contiguous E row and spies restore after errors and interrupts", {
  e <- runs48$plain$early$value
  for (m in e$reduction_summary$dropped_indices) {
    rec <- e$point_links_all[[m]]
    hit <- which(e$structural_certificate$E$i == m)
    expect_identical(rec$expansion$retained_rows, e$structural_certificate$E$j[hit])
    expect_identical(rec$expansion$signs, e$structural_certificate$E$x[hit])
  }
  nm <- "fit_ipw_point_link_view"
  original <- api[[nm]]
  locked <- bindingIsLocked(nm, api)
  expect_error(er_count_calls(nm, function() stop("test spy error")), "test spy error")
  expect_identical(api[[nm]], original)
  expect_identical(bindingIsLocked(nm, api), locked)
  result <- tryCatch(er_count_calls(nm, function()
    stop(structure(list(message = "test spy interrupt", call = NULL), class = c("interrupt", "condition")))),
    interrupt = function(e) "interrupted")
  expect_identical(result, "interrupted")
  expect_identical(api[[nm]], original)
  expect_identical(bindingIsLocked(nm, api), locked)
})
