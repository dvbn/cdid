# CP08 scaling corrections (typed catalog/design envelopes, registry
# mask-equivalence certificate, exact auxiliary reporting plan, live reference
# source charge) and inherited execution options (compute_inffunc, pl/cores,
# print_details). Strict actual-installed namespace tests; did 2.5.1 is needed
# for the public att_gt()/aggte() blocks. Expectations come from independent
# enumerations written here, frozen copies of the accepted algorithms, exact
# old-formula arithmetic, or the unchanged accepted computations themselves.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api), environmentIsLocked(api))

cp8_rebind <- function(name, fn) {
  locked <- bindingIsLocked(name, api)
  if (locked) unlockBinding(name, api)
  assign(name, fn, envir = api)
  if (locked) lockBinding(name, api)
}
cp8_count <- function(names, run) {
  originals <- lapply(stats::setNames(names, names), function(nm) get(nm, envir = api, inherits = FALSE))
  counts <- stats::setNames(integer(length(names)), names)
  on.exit(for (nm in names) cp8_rebind(nm, originals[[nm]]), add = TRUE)
  for (nm in names) local({
    nm0 <- nm; f0 <- originals[[nm]]
    cp8_rebind(nm0, function(...) { counts[[nm0]] <<- counts[[nm0]] + 1L; f0(...) })
  })
  value <- run()
  list(value = value, counts = counts)
}
cp8_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (inherits(x, "cdid_ipw_structural_certificate")) {
    x <- unclass(x); return(lapply(x[setdiff(names(x), c("link_kind", "object_refs", "ref"))], cp8_strip))
  }
  if (is.list(x)) { a <- attributes(x); y <- lapply(x, cp8_strip); attributes(y) <- a; return(y) }
  x
}
cp8_grid_occ <- function(calendar, cohorts, rule) {
  out <- list()
  for (si in seq_along(calendar)) for (ti in seq_along(calendar)) if (si < ti)
    for (ki in seq_along(cohorts)) out[[length(out) + 1L]] <- data.frame(cohort = cohorts[ki], s = calendar[si],
      t = calendar[ti], count = rule(si, ti, ki))
  do.call(rbind, out)
}
cp8_caps <- function(max_work = 1e9, max_scratch = 1e8)
  list(max_candidates = 1e6, max_rows = 1e6, max_pool_entries = 1e7, max_scratch = max_scratch, max_work = max_work)
cp8_qualified <- function(cat) {
  r <- cat$rows
  s_c <- r$s_role %in% c("anticipation_conditional", "post_conditional")
  t_c <- r$t_role %in% c("anticipation_conditional", "post_conditional")
  dq <- unique(data.frame(g = c(r$g[s_c], r$g[t_c]), t = c(r$s[s_c], r$t[t_c])))
  if (!nrow(dq)) return(NULL)
  dq <- dq[order(dq$g, dq$t), , drop = FALSE]; rownames(dq) <- NULL; dq
}
# Accepted pre-CP08 adapter work formula (rows-squared graph term), for arithmetic only.
cp8_old_graph_work <- function(d) {
  H <- d$R + d$C + d$T + d$P + d$J + d$D + d$Q + 1
  512 * (d$C + 1 + ceiling(log2(H + 1))) * H * H
}

# ---- 1. typed catalog/design accounting ---------------------------------------
test_that("typed adapter envelopes are linear in actual dimensions; the D62 rows-squared charge is gone", {
  d62 <- list(R = 4950, C = 2, T = 100, P = 1, NP = 1, PR = 4950, J = 1, D = 98, Q = 98, O = 1, B = 24, Bf = 60)
  expect_identical(cp8_old_graph_work(d62), 225878024192)
  e <- api$catalog_design_envelope(d62, design = TRUE)
  expect_lt(e$work, 1e9); expect_lt(e$scratch, 134217728 / 4)
  e2 <- api$catalog_design_envelope(modifyList(d62, list(R = 9900, PR = 9900)), design = TRUE)
  expect_lt(e2$work / e$work, 2.5)                       # doubling rows does not quadruple work
  lower <- api$catalog_design_envelope(modifyList(d62, list(R = 0, P = 0, NP = 0, PR = 0, O = 0, Bf = 0)), design = TRUE)
  expect_lte(lower$work, e$work); expect_lte(lower$scratch, e$scratch)
  links <- api$catalog_design_envelope(d62, design = FALSE)
  expect_lt(links$work, e$work); expect_lt(links$scratch, e$scratch)
})

test_that("a large-row design builds under bounds the former charge refused; support equals an independent union-find", {
  cal <- 1:24; coh <- c(8, 16, Inf)
  cat <- api$build_comparison_catalog(cal, coh, c(8, 16), 0, cp8_grid_occ(cal, coh, function(si, ti, ki) 1),
    "notyettreated", "conv:cp08", NULL, cp8_caps())
  dq <- cp8_qualified(cat)
  tpt <- expand.grid(g = c(8, 16), t = cal); tpt <- tpt[tpt$t >= tpt$g, , drop = FALSE]
  pre <- api$catalog_adapter_preflight(cat, dq, tpt, design = TRUE)
  d <- pre$dimensions
  old <- cp8_old_graph_work(d)
  expect_gt(nrow(cat$rows), 400L)
  expect_lt(pre$work * 10, old)
  bounded <- cat; bounded$resources$bounds$max_work <- floor(sqrt(pre$work * old))   # between the two charges
  des <- api$build_comparison_catalog_design(bounded, dq, tpt)
  # Independent union-find over the design's own link endpoints.
  nv <- nrow(des$vertices); parent <- seq_len(nv)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
  for (k in seq_len(nrow(des$links))) {
    a <- find(des$links$start_vertex_index[k]); b <- find(des$links$end_vertex_index[k]); if (a != b) parent[a] <- b
  }
  roots <- vapply(seq_len(nv), find, integer(1))
  anchored <- unique(roots[des$vertices$role == "reference"])
  pv <- which(!is.na(des$vertices$param_index))
  expect_identical(des$post_targets$level_estimable,
    roots[pv[match(des$post_targets$param_index, des$vertices$param_index[pv])]] %in% anchored)
  expect_identical(des$rank_diagnostics$n_components, length(unique(roots)))
})

test_that("forests equal a frozen copy of the accepted per-component anchor scan, including PT-all designs", {
  accepted_forest <- function(nv, from, to, anchors, edge_order = seq_along(from)) {
    m <- length(from); anchors <- unique(as.integer(anchors))
    u <- c(from, to); v <- c(to, from); e <- rep(seq_len(m), 2L); sign <- c(rep(1L, m), rep(-1L, m))
    ord <- order(u, v, rep(edge_order, 2L), method = "radix"); degree <- tabulate(u, nbins = nv)
    v <- v[ord]; e <- e[ord]; sign <- sign[ord]
    end <- cumsum(as.double(degree)); start <- end - degree + 1
    comp <- integer(nv); parent <- parent_edge <- integer(nv); parent_sign <- integer(nv)
    queue <- integer(nv); roots <- integer(nv); nc <- 0L
    for (root in c(anchors, seq_len(nv))) {
      if (comp[root] != 0L) next
      nc <- nc + 1L; roots[nc] <- root; comp[root] <- nc; head <- tail <- 1L; queue[1L] <- root
      while (head <= tail) { a <- queue[head]; head <- head + 1L; if (degree[a] == 0L) next
        for (k in seq.int(start[a], end[a])) { b <- v[k]; if (comp[b] == 0L) { comp[b] <- nc; parent[b] <- a; parent_edge[b] <- e[k]
          parent_sign[b] <- sign[k]; tail <- tail + 1L; queue[tail] <- b } } }
    }
    if (length(anchors) > 1L && anyDuplicated(comp[anchors]) > 0L) {
      parent <- parent_edge <- parent_sign <- integer(nv); vis <- logical(nv)
      for (cid in seq_len(nc)) {
        comp_anchors <- anchors[comp[anchors] == cid]
        seeds <- if (length(comp_anchors) > 0L) comp_anchors else roots[cid]
        head <- 1L; tail <- length(seeds); queue[seq_len(tail)] <- seeds; vis[seeds] <- TRUE
        while (head <= tail) { a <- queue[head]; head <- head + 1L; if (degree[a] == 0L) next
          for (k in seq.int(start[a], end[a])) { b <- v[k]; if (!vis[b]) { vis[b] <- TRUE; parent[b] <- a; parent_edge[b] <- e[k]
            parent_sign[b] <- sign[k]; tail <- tail + 1L; queue[tail] <- b } } }
      }
    }
    list(component_id = comp, n_components = nc, roots = roots[seq_len(nc)],
         parent_vertex = parent, parent_link = parent_edge, parent_sign = parent_sign, degree = degree)
  }
  graphs <- list(list(nv = 9L, from = c(1L, 2L, 4L, 5L, 7L, 1L), to = c(2L, 3L, 5L, 6L, 8L, 3L), anchors = c(3L, 1L, 6L, 9L)),
    list(nv = 7L, from = c(1L, 1L, 2L, 3L, 5L), to = c(2L, 3L, 4L, 4L, 6L), anchors = c(4L, 2L, 7L, 1L, 6L)),
    list(nv = 6L, from = integer(), to = integer(), anchors = c(2L, 2L, 5L)))
  for (g in graphs) expect_identical(api$incidence_forest(g$nv, g$from, g$to, g$anchors), accepted_forest(g$nv, g$from, g$to, g$anchors))
  cal <- 1:6; coh <- c(4, 5, Inf)
  cat <- api$build_comparison_catalog(cal, coh, c(4, 5), 0, cp8_grid_occ(cal, coh, function(si, ti, ki) 1),
    "notyettreated", "conv:cp08", NULL, cp8_caps())
  des <- api$build_comparison_catalog_design(cat, cp8_qualified(cat), data.frame(g = c(4, 4, 5), t = c(4, 6, 6)),
    pt_assumption = "all", pt_start = 1)
  anchors <- which(des$vertices$role == "reference" | (des$vertices$role == "pre" & des$vertices$time >= 1))
  expect_gt(sum(duplicated(des$forest$component_id[anchors])), 0L)        # multi-anchor re-rooting is exercised
  expect_identical(des$forest, accepted_forest(nrow(des$vertices), des$links$start_vertex_index,
    des$links$end_vertex_index, anchors, order(order(match(des$links$g, sort(unique(des$links$g))), des$links$s, des$links$t,
      match(des$links$functional_key, sort(unique(des$links$functional_key), method = "radix")), method = "radix"))))
})

test_that("malformed catalog headers refuse before allocation and the early owner mirror evaluates the same envelope", {
  cal <- 1:6; coh <- c(3, 4, 5, Inf)
  cat <- api$build_comparison_catalog(cal, coh, c(3, 4), 0, cp8_grid_occ(cal, coh, function(si, ti, ki) as.numeric((si + ti + ki) %% 4 != 0)),
    "notyettreated", "conv:cp08", NULL, cp8_caps())
  bad <- cat; bad$rows$pool_id[1L] <- bad$pools$n_pools + 1L
  expect_error(api$catalog_adapter_preflight(bad, NULL), class = "cdid_catalog_input_error")
  bad <- cat; bad$pools$size[1L] <- bad$pools$size[1L] + 1L
  expect_error(api$catalog_adapter_preflight(bad, NULL), class = "cdid_catalog_input_error")
  small <- cat; small$resources$bounds$max_work <- 10
  expect_error(api$build_comparison_catalog_design(small, NULL, 4), class = "cdid_catalog_bound_error")
  # Early owner preflight: the mirror refuses exactly when its lower-bound envelope exceeds max_work.
  d <- data.frame(id = rep(1:12, each = 4), t = rep(1:4, 12), g = rep(c(3, 4, Inf), each = 16), y = 0)
  d$y <- d$id + d$t
  oc <- api$ipw_owner_control(max_chunk_rows = 64L, max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", oc)
  tpt <- data.frame(g = c(3, 3, 4), t = c(3, 4, 4))
  lower <- list(R = 0, C = 3, T = 4, P = 0, NP = 0, PR = 0, J = 2, D = 0, Q = 3, O = 0, B = nchar("cp08", type = "bytes"), Bf = 0)
  env <- api$catalog_design_envelope(lower, design = TRUE)
  b <- cp8_caps(max_work = env$work - 1)
  expect_error(api$preflight_owner_comparison_catalog(po, c(3, 4), 0, "notyettreated", "cp08", control = oc, catalog_bounds = b,
    target_post_times = tpt, build_design = TRUE), class = "cdid_catalog_bound_error")
  expect_silent(api$preflight_owner_comparison_catalog(po, c(3, 4), 0, "notyettreated", "cp08", control = oc,
    catalog_bounds = cp8_caps(max_work = env$work), target_post_times = tpt, build_design = TRUE))
})

test_that("design validation keeps its per-row refusals with the typed pool cache", {
  cal <- 1:6; coh <- c(3, 4, 5, Inf)
  cat <- api$build_comparison_catalog(cal, coh, c(3, 4), 0, cp8_grid_occ(cal, coh, function(si, ti, ki) 1),
    "notyettreated", "conv:cp08", NULL, cp8_caps())
  links <- api$comparison_catalog_design_links(cat, NULL)
  spec <- api$validate_calendar_spec(cal, coh, 0)
  check <- function(mut, cls, pattern) {
    l <- mut(links); attr(l, "calendar_spec") <- spec
    e <- expect_error(api$validate_design_links(l, spec), class = cls)
    expect_match(conditionMessage(e), pattern, fixed = TRUE)
  }
  check(function(l) { l$control_pool[2L] <- "{3,9}"; l }, "cdid_design_cohort_error", "must belong to the supplied cohort universe")
  check(function(l) { l$control_pool[3L] <- "{x}"; l }, "cdid_design_qualification_error", "Invalid `control_pool` cohort set syntax")
  check(function(l) { l$control_pool[1L] <- "3,4"; l }, "cdid_design_qualification_error", "explicit `{...}` cohort set")
  # The first offending row decides, as in a per-row parse.
  check(function(l) { l$control_pool[2L] <- "3,4"; l$control_pool[5L] <- "{x}"; l }, "cdid_design_qualification_error", "explicit `{...}` cohort set")
  ok <- links; attr(ok, "calendar_spec") <- spec
  expect_identical(api$validate_design_links(ok, spec)$control_pool, links$control_pool)
})

# ---- 2. registry mask-equivalence certificate ---------------------------------
cp8_panel48 <- function() {
  set.seed(20261016L)
  u_g <- rep(c(3, 4, Inf), each = 16L); u_x <- rep(c(-1.0, 1.0, -0.5, 0.5), times = 12L); u_w <- rep(c(1.0, 1.5, 0.8, 1.2), times = 12L)
  d <- expand.grid(t = 1:4, id = 1:48, KEEP.OUT.ATTRS = FALSE)
  d$g <- u_g[d$id]; d$x <- u_x[d$id]; d$weight <- u_w[d$id]
  eff <- ifelse(d$g == 3 & d$t >= 3, 2.0 + 0.5 * (d$t - 3), ifelse(d$g == 4 & d$t >= 4, 1.5, 0))
  d$y <- 0.25 * d$id + 0.5 * d$t + 0.3 * d$x * d$t + eff + 0.2 * stats::rnorm(nrow(d))
  d[(d$id %% 7 != 0) | (d$t != 2), c("id", "t", "g", "x", "weight", "y")]   # unequal observation masks
}
cp8_ctrl48 <- function(n_rows, tuning = NULL) {
  dims <- api$chained_public_dimensions_from_counts(n_frame = 48L, n_panel_units = 48L, n_rows = n_rows, n_dates = 4L,
    n_cohorts = 3L, n_targets = 2L, n_covariates = 2L, biters = 99, n_clusters = 48L)
  api$chained_public_pipeline_control(dims, api$chained_public_resolve_control(tuning, gmm_weighting = "optimal"),
    gmm_weighting = "optimal")$control
}
cp8_run48 <- function(d, ctrl, stage, ...) api$run_chained_ipw_pipeline(data = d, idname = "id", tname = "t", yname = "y",
  gname = "g", xformla = ~ x, weightsname = "weight", control = ctrl, control_group = "notyettreated", base_period = "varying",
  weighting = "optimal", run_numerical_transfer = FALSE, root_witness = NULL, structural_reduction_stage = stage, ...)
d48 <- cp8_panel48(); ctrl48 <- cp8_ctrl48(nrow(d48))
runs48 <- list(early = cp8_run48(d48, ctrl48, "outcome_free", aggregations = character(0L), cband = FALSE, bstrap = FALSE),
  late = cp8_run48(d48, ctrl48, "post_point_fit", aggregations = character(0L), cband = FALSE, bstrap = FALSE))
reg48 <- api$build_pipeline_nuisance_and_link_views(metadata_obj = runs48$early$metadata, frame_owner = runs48$early$frame_owner,
  owner_control = ctrl48$owner_control, fitted_control = ctrl48$fitted_control, point_control = ctrl48$point_control,
  glm_control = ctrl48$glm_control, fit_point_links = FALSE)
dli48 <- match(seq_along(reg48$link_views), runs48$early$design_all$links$catalog_row)

test_that("the sealed registry certificate replaces every streamed mask comparison; classes and E are unchanged", {
  M <- length(reg48$link_views)
  expect_gt(length(reg48$fitted_views), length(unique(reg48$operator_ids)))  # several mask views per operator
  with <- cp8_count("verify_ipw_view_mask_compatibility", function()
    api$certify_ipw_structural_classes_from_views(reg48, runs48$early$design_all, dli48, ctrl48$class_bounds, ctrl48$owner_control))
  noc <- reg48; noc$mask_certificate <- NULL
  without <- cp8_count("verify_ipw_view_mask_compatibility", function()
    api$certify_ipw_structural_classes_from_views(noc, runs48$early$design_all, dli48, ctrl48$class_bounds, ctrl48$owner_control))
  expect_identical(unname(with$counts), 0L)
  expect_gte(unname(without$counts), M)
  expect_identical(cp8_strip(with$value), cp8_strip(without$value))
  # Early certification precedes auxiliary intake; late certification retains
  # its additional auxiliary nuisance view. The backend ref owns the raw backend
  # itself and is never a nuisance-view ref (raw construction, sections 3/4/7).
  early_nvs <- runs48$early$nuisance_registry$nuisance_views
  late_nvs <- runs48$late$nuisance_registry$raw_aux_nuisance_views
  expect_identical(length(early_nvs), length(reg48$nuisance_keys))
  expect_identical(length(late_nvs), length(early_nvs) + 1L)
  expect_identical(names(early_nvs), reg48$nuisance_keys)
  expect_identical(names(late_nvs), c(reg48$nuisance_keys, ""))
  expect_identical(runs48$early$structural_certificate$nuisance_refs, lapply(early_nvs, `[[`, "ref"))
  expect_identical(runs48$late$structural_certificate$nuisance_refs, lapply(late_nvs, `[[`, "ref"))
  for (run in runs48) {
    raw <- run$nuisance_registry$raw_ipw
    expect_s3_class(raw, "cdid_ipw_raw_backend")
    expect_identical(ls(raw$ref, all.names = TRUE), "ipw_raw_backend")
    expect_identical(get("ipw_raw_backend", raw$ref, inherits = FALSE), raw)
    expect_false(any(vapply(run$structural_certificate$nuisance_refs, identical, logical(1L), raw$ref)))
  }
  for (j in seq_along(late_nvs)) {
    ref <- runs48$late$structural_certificate$nuisance_refs[[j]]
    expect_identical(get("nuisance_view", ref, inherits = FALSE), late_nvs[[j]])
  }
  early_cert_cmp <- runs48$early$structural_certificate
  names(early_cert_cmp$nuisance_refs) <- NULL
  late_cert_cmp <- runs48$late$structural_certificate
  late_cert_cmp$nuisance_refs <- unname(late_cert_cmp$nuisance_refs[seq_along(reg48$nuisance_keys)])
  expect_identical(cp8_strip(with$value), cp8_strip(early_cert_cmp))
  expect_identical(cp8_strip(with$value), cp8_strip(late_cert_cmp))
  expect_identical(runs48$early$reduction_summary$kept_indices, runs48$late$reduction_summary$kept_indices)
  expect_identical(runs48$early$reporting$estimates, runs48$late$reporting$estimates)
  pf <- api$preflight_ipw_early_structural_reduction(reg48, runs48$early$design_all, ctrl48$class_bounds, ctrl48$owner_control)
  pn <- api$preflight_ipw_early_structural_reduction(noc, runs48$early$design_all, ctrl48$class_bounds, ctrl48$owner_control)
  expect_identical(pf$mask_streams_charged, 0)
  expect_identical(pn$mask_streams_charged, M * (length(reg48$fitted_views) + 1))
  expect_lt(pf$work_ops_bound, pn$work_ops_bound)
  expect_match(runs48$early$resource_accounting$early_structural$mask_evidence, "certificate")
})

test_that("mismatched references, broken seals and foreign certificates refuse before any stream", {
  certify <- function(r) cp8_count("verify_ipw_view_mask_compatibility", function() tryCatch(
    api$certify_ipw_structural_classes_from_views(r, runs48$early$design_all, dli48, ctrl48$class_bounds, ctrl48$owner_control),
    error = identity))
  i <- which(reg48$link_to_nuisance_idx == reg48$link_to_nuisance_idx[1L])
  expect_gte(length(i), 2L)
  swap <- reg48; tmp <- swap$link_views[[i[1L]]]; swap$link_views[[i[1L]]] <- swap$link_views[[i[2L]]]; swap$link_views[[i[2L]]] <- tmp
  z <- certify(swap); expect_s3_class(z$value, "cdid_ipw_certificate_error"); expect_identical(unname(z$counts), 0L)
  bad <- reg48; bad$mask_certificate$view_operator[1L] <- bad$mask_certificate$view_operator[1L] + 1L
  z <- certify(bad); expect_s3_class(z$value, "cdid_ipw_shared_mask_error")
  expect_match(conditionMessage(z$value), "sealed snapshot"); expect_identical(unname(z$counts), 0L)
  foreign <- reg48; foreign$mask_certificate <- runs48$late$nuisance_registry$mask_certificate
  z <- certify(foreign); expect_s3_class(z$value, "cdid_ipw_certificate_error"); expect_identical(unname(z$counts), 0L)
})

test_that("adversarial fingerprint collisions are resolved exactly and every registry scan is charged before it runs", {
  obs <- vector("list", 24L)
  for (i in c(1, 5, 8, 12)) obs[[i]] <- c(1, 2)
  for (i in c(2, 3, 10, 11)) obs[[i]] <- c(1, 3)
  for (i in c(4, 6, 7, 9)) obs[[i]] <- c(2, 4)
  for (i in 13:24) obs[[i]] <- 1:4
  g <- ifelse(seq_len(24) %% 2 == 1, 4, Inf)
  df <- do.call(rbind, lapply(seq_along(obs), function(i) { tt <- as.double(obs[[i]])
    data.frame(id = i, t = tt, g = g[i], x = ((5 * i) %% 7) / 7, y = 0.5 * i + 0.25 * tt + ifelse(tt >= g[i], 1.5, 0) + ((7 * i + 3 * tt) %% 5) / 10) }))
  oc <- api$ipw_owner_control(max_chunk_rows = 64L, max_retained_entries = 2e6, max_matrix_entries = 2e6, max_work_ops = 5e8, max_output_entries = 5e5)
  fc <- api$ipw_fitted_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2e6, max_matrix_entries = 2e6, max_work_ops = 5e8, max_output_entries = 5e5)
  pc <- api$ipw_point_view_control(max_chunk_rows = 64L, max_observer_units = 2000L, max_retained_entries = 2e6, max_matrix_entries = 2e6, max_work_ops = 5e8, max_output_entries = 5e5)
  glm <- list(epsilon = 1e-12, maxit = 50L, trace = FALSE)
  po <- api$construct_ipw_panel_owner(df, "id", "t", "y", "g", oc)
  u <- df[!duplicated(df$id), ]; u <- u[order(u$id), ]
  fo <- api$construct_ipw_frame_owner(po, u$id, cbind(1, u$x), control = oc)
  meta <- api$build_owner_comparison_catalog(po, targets = 4, anticipation = 0, control_group = "nevertreated", conventions = "p1",
    control = oc, catalog_bounds = list(max_candidates = 1e4, max_rows = 1e4, max_pool_entries = 1e5, max_scratch = 1e7, max_work = 5e8),
    frame_owner = fo)
  reg <- api$build_pipeline_nuisance_and_link_views(meta, fo, oc, fc, pc, glm, fit_point_links = FALSE)
  expect_identical(reg$mask_certificate$separations$exact_streamed_mismatches, 1)
  expect_identical(reg$construction_work$exact_mask_mismatches, 1)
  expect_gt(reg$construction_work$work_ops_charged, 0)
  # A composed envelope with room for less than one row's streams refuses before any view is built.
  z <- cp8_count("construct_ipw_nuisance_view", function() tryCatch(api$build_pipeline_nuisance_and_link_views(meta, fo, oc, fc, pc, glm,
    fit_point_links = FALSE, external_work_ops = oc$max_work_ops - 1), error = identity))
  expect_s3_class(z$value, "cdid_ipw_owner_bound_error")
  # The descriptor context is admitted before the first registry scan; the
  # caller's explicit finite cap must already cover that work.
  expect_match(conditionMessage(z$value), "caller plus descriptor preparation work", fixed = TRUE)
  expect_identical(unname(z$counts), 0L)
})

# ---- 3. exact auxiliary reporting plan ----------------------------------------
# Independent enumeration of both families (did varying/universal pretrend cells).
cp8_aux_reference <- function(groups, cal, coh, ant, cg) {
  cal <- sort(cal); coh <- sort(coh); desc <- list(); key <- character(); sel <- list(varying = list(), universal = list())
  for (bp in c("varying", "universal")) for (g in sort(unique(groups))) {
    ref <- if (any(cal + ant < g)) max(cal[cal + ant < g]) else NA_real_
    times <- if (bp == "varying") cal[-1L] else cal
    for (t in times[times < g]) {
      if (is.na(ref)) { sel[[bp]][[length(sel[[bp]]) + 1L]] <- list(g = g, t = t, status = "unavailable", index = NA_integer_); next }
      from <- if (bp == "varying") cal[match(t, cal) - 1L] else ref
      if (bp == "universal" && t == ref) { sel[[bp]][[length(sel[[bp]]) + 1L]] <- list(g = g, t = t, status = "reference_zero", index = 0L); next }
      u <- max(t, from)
      K <- if (cg == "nevertreated") coh[is.infinite(coh)] else coh[coh != g & u + ant < coh]
      if (!length(K)) { sel[[bp]][[length(sel[[bp]]) + 1L]] <- list(g = g, t = t, status = "unavailable", index = NA_integer_); next }
      k <- paste(g, min(t, from), u, paste(K, collapse = ","))
      hit <- match(k, key)
      if (is.na(hit)) { key <- c(key, k); desc[[length(desc) + 1L]] <- list(g = g, s = min(t, from), t = u, K = K); hit <- length(desc) }
      sel[[bp]][[length(sel[[bp]]) + 1L]] <- list(g = g, t = t, status = "pending", index = hit, sign = if (t > from) 1 else -1)
    }
  }
  list(descriptors = desc, selections = sel)
}

test_that("the auxiliary union plan equals an independent enumeration of both reporting families", {
  cases <- list(list(g = c(3, 4), cal = 1:6, coh = c(3, 4, 5, Inf), ant = 0, cg = "notyettreated"),
    list(g = c(3, 4), cal = 1:6, coh = c(3, 4, 5, Inf), ant = 0, cg = "nevertreated"),
    list(g = c(4, 6), cal = 1:7, coh = c(4, 6, Inf), ant = 1, cg = "notyettreated"),
    list(g = c(2, 5), cal = c(1, 2, 4.5, 5, 7), coh = c(2, 5, 7), ant = 0, cg = "notyettreated"),
    list(g = 1, cal = 1:4, coh = c(1, 3, Inf), ant = 0, cg = "nevertreated"))
  for (z in cases) {
    ledger <- api$ipw_pipeline_request_ledger(z$g, z$cal, z$coh, z$ant, z$cg)
    got <- api$ipw_auxiliary_spec(ledger, z$cal, z$coh, z$ant, z$cg)
    ref <- cp8_aux_reference(ledger$g, z$cal, z$coh, z$ant, z$cg)
    expect_identical(lapply(got$descriptors, function(d) unname(unlist(d))), lapply(ref$descriptors, function(d) unname(unlist(d))))
    for (bp in c("varying", "universal")) {
      expect_identical(vapply(got$selections[[bp]], `[[`, character(1), "status"), vapply(ref$selections[[bp]], `[[`, character(1), "status"))
      expect_identical(vapply(got$selections[[bp]], function(s) as.integer(s$index), integer(1)),
        vapply(ref$selections[[bp]], function(s) as.integer(s$index), integer(1)))
    }
  }
})

test_that("pipelines reserve the exact plan and record it; reports keep both families' zero and unavailable cells", {
  ap <- runs48$early$resource_accounting$auxiliary_plan
  ref <- cp8_aux_reference(c(3, 4), 1:4, c(3, 4, Inf), 0, "notyettreated")
  expect_identical(ap$n_descriptors_union, as.double(length(ref$descriptors)))
  expect_identical(ap$n_reported_pretrend, as.double(length(ref$selections$varying)))
  expect_identical(ap$replaced_bound, 16)
  expect_identical(runs48$early$auxiliary_reporting$n_descriptors, length(ref$descriptors))
  u <- cp8_run48(d48, ctrl48, "outcome_free", aggregations = character(0L), cband = FALSE, bstrap = FALSE)
  uu <- api$run_chained_ipw_pipeline(data = d48, idname = "id", tname = "t", yname = "y", gname = "g", xformla = ~ x,
    weightsname = "weight", control = ctrl48, control_group = "notyettreated", base_period = "universal", weighting = "optimal",
    run_numerical_transfer = FALSE, root_witness = NULL, aggregations = character(0L), cband = FALSE, bstrap = FALSE)
  expect_identical(uu$resource_accounting$auxiliary_plan$n_descriptors_union, ap$n_descriptors_union)   # one union for both families
  expect_identical(uu$resource_accounting$auxiliary_plan$n_reference_zero, 2L)
  expect_identical(sum(uu$reporting$target_meta$support_status == "deterministic_reference_zero"), 2L)
  post <- u$reporting$target_meta$role == "post"; postu <- uu$reporting$target_meta$role == "post"
  expect_identical(unname(u$reporting$estimates[post]), unname(uu$reporting$estimates[postu]))
})

# ---- 7. live reference source before the requirement scan ---------------------
cp8_rid_panel <- function(n = 36, dates = c(1, 2, 4, 5, 7), cohorts = c(4, 5, 0), drop_y = TRUE) {
  id <- seq_len(n); g <- rep(cohorts, length.out = n)
  d <- expand.grid(time = dates, id = id)[, c("id", "time")]
  ti <- match(d$time, dates); d$G <- g[d$id]
  bx <- round(sin(1.7 * id), 2)
  d$x <- round(bx[d$id] + 0.25 * ti + 0.3 * cos(2.3 * d$id + 1.1 * ti), 2)
  d$w <- round(round(0.6 + (id %% 7) / 7, 3)[d$id] * (1 + 0.15 * ti * ((d$id %% 2) == 0)), 3)
  te <- ifelse(d$G > 0 & d$time >= d$G, 1 + 0.5 * (d$time - d$G), 0)
  d$y <- round(cos(d$id) + 0.4 * d$time + 0.8 * d$x + te + 0.5 * sin(5.1 * d$id + 2.7 * ti), 3)
  if (drop_y) d$y[(d$id * 7 + ti * 3) %% 9 == 0] <- NA
  d
}

test_that("the dated reference reader charges its live source before the requirement scan", {
  skip_if_not_installed("did")
  captured <- NULL
  orig <- get("build_chained_reference_reader", envir = api)
  on.exit(cp8_rebind("build_chained_reference_reader", orig), add = TRUE)
  cp8_rebind("build_chained_reference_reader", function(src, reader, views, n_post, ctrl, retained = 0, work = 0) {
    if (is.null(captured)) captured <<- list(src = src, reader = reader, views = views, n_post = n_post, ctrl = ctrl,
      retained = retained, work = work)
    orig(src, reader, views, n_post, ctrl, retained, work)
  })
  invisible(cdid::att_gt(yname = "y", tname = "time", idname = "id", gname = "G", data = cp8_rid_panel(), xformla = ~x,
    weightsname = "w", est_method = "reg", fix_weights = "base_period", bstrap = FALSE, cband = FALSE, gmm_weighting = "identity"))
  cp8_rebind("build_chained_reference_reader", orig)
  expect_false(is.null(captured))
  st <- captured$reader$state
  n_roles <- 1 + isTRUE(captured$src$formulas$x$dated) + (!is.null(captured$src$formulas$s) && isTRUE(captured$src$formulas$s$dated))
  n_req <- 0
  for (v in captured$views) for (k in match(c(v$cohorts$g, v$cohorts$K), st$cohorts))
    n_req <- n_req + n_roles * (st$cohort_end[k] - st$cohort_start[k] + 1)
  source <- unname(api$cp04_payload(captured$src)["entries"])
  expect_gt(source, 0)
  tight <- captured$ctrl; tight$max_retained_entries <- captured$retained + 8 * n_req + source - 1
  z <- cp8_count("chained_reference_requirements", function() tryCatch(api$build_chained_reference_reader(captured$src,
    captured$reader, captured$views, captured$n_post, tight, captured$retained, captured$work), error = identity))
  expect_s3_class(z$value, "cdid_ipw_owner_bound_error")
  expect_match(conditionMessage(z$value), "reference requirement scan", fixed = TRUE)
  expect_identical(unname(z$counts), 0L)
  roomy <- captured$ctrl; roomy$max_retained_entries <- roomy$max_retained_entries + source
  ref <- api$build_chained_reference_reader(captured$src, captured$reader, captured$views, captured$n_post, roomy,
    captured$retained, captured$work)
  expect_identical(ref$accounting$source_entries, source)
})

# ---- 4-6. inherited execution options (public, did 2.5.1) ---------------------
cp8_route_panel <- function() {   # D60 route fixture: every date observed except a deterministic pattern
  d <- cp8_rid_panel(n = 36, dates = c(1, 2, 4, 5), drop_y = FALSE)
  ti <- match(d$time, c(1, 2, 4, 5))
  d$y[((d$id * 37 + ti * 11) %% 13) < 3] <- NA
  d
}
cp8_public_panel <- function() {
  d <- cp8_panel48(); d$g[!is.finite(d$g)] <- 0; d
}
cp8_fit <- function(d, ...) cdid::att_gt(yname = "y", tname = "t", idname = "id", gname = "g", data = d, xformla = ~ x,
  weightsname = "weight", control_group = "notyettreated", ...)

test_that("compute_inffunc = FALSE returns did's point-only result with estimates equal to the full fit on every route", {
  skip_if_not_installed("did")
  dp <- cp8_public_panel(); dr <- cp8_rid_panel(); rp <- cp8_route_panel()
  routes <- list(
    list(d = dp, args = list(est_method = "ipw", control = list(condition_trigger = 0.99, ridge_relative = 0.01))),
    list(d = dp, args = list(est_method = "reg")), list(d = dp, args = list(est_method = "dr")),
    list(d = dr, args = list(est_method = "reg", fix_weights = "base_period", gmm_weighting = "identity"), dated = TRUE),
    list(d = rp, args = list(est_method = "reg", fix_weights = "base_period", gmm_weighting = "identity",
      selection = "conditional_mtar", selection_est_method = "imputation"), dated = TRUE),
    list(d = rp, args = list(est_method = "dr", fix_weights = "base_period", gmm_weighting = "identity",
      selection = "conditional_mtar", selection_est_method = "dr"), dated = TRUE))
  for (r in routes) {
    call <- function(...) if (isTRUE(r$dated)) do.call(cdid::att_gt, c(list(yname = "y", tname = "time", idname = "id", gname = "G",
      data = r$d, xformla = ~x, weightsname = "w", control_group = "notyettreated"), r$args, list(...))) else
      do.call(cp8_fit, c(list(r$d), r$args, list(...)))
    full <- call(bstrap = FALSE, cband = FALSE)
    pt <- call(compute_inffunc = FALSE)              # bstrap/cband defaults TRUE become FALSE, as in did
    expect_identical(pt$att, full$att)
    expect_identical(pt$chained$available, full$chained$available)
    expect_null(pt$inffunc); expect_null(pt$V_analytical); expect_null(pt$W); expect_null(pt$Wpval)
    expect_true(all(is.na(pt$se))); expect_equal(pt$c, stats::qnorm(0.975), tolerance = 0)
    expect_false(pt$DIDparams$bstrap); expect_false(pt$DIDparams$cband); expect_false(pt$DIDparams$compute_inffunc)
    expect_identical(pt$chained$inference$kind, "not_computed")
    expect_identical(pt$chained$covariance_handling$handler_used, full$chained$covariance_handling$handler_used)
    expect_identical(pt$chained$pretest$status, "not_computed")
    expect_error(cdid::aggte(pt, type = "simple"), class = "cdid_chained_request_no_inffunc_error")
  }
  ridge <- cp8_fit(dp, est_method = "ipw", compute_inffunc = FALSE, control = list(condition_trigger = 0.99, ridge_relative = 0.01))
  expect_identical(ridge$chained$covariance_handling$handler_used, "tikhonov")
  out <- utils::capture.output(summary(ridge))
  expect_true(any(grepl("not computed (compute_inffunc = FALSE", out, fixed = TRUE)))
})

test_that("pl/cores need declared workers; admitted workers reproduce cores = 1 draws, bands and aggte exactly", {
  skip_if_not_installed("did")
  d <- cp8_public_panel()
  expect_error(cp8_fit(d, est_method = "reg", pl = TRUE, cores = 2, biters = 49), class = "cdid_chained_request_parallel_error")
  set.seed(3); s <- cp8_fit(d, est_method = "reg", biters = 49); s_rng <- .Random.seed
  set.seed(3); p <- cp8_fit(d, est_method = "reg", biters = 49, pl = TRUE, cores = 2, control = list(workers = 2)); p_rng <- .Random.seed
  expect_identical(p_rng, s_rng)
  for (f in c("att", "se", "c", "inffunc", "V_analytical")) expect_identical(p[[f]], s[[f]], info = f)
  expect_identical(p$chained$request$execution$workers, 2)
  expect_true(p$DIDparams$pl); expect_identical(p$DIDparams$cores, 2L)
  expect_identical(p$chained$inference$multiplier$metadata$execution$workers, 2)
  set.seed(9); sa <- cdid::aggte(s, type = "dynamic"); sa_rng <- .Random.seed
  set.seed(9); pa <- cdid::aggte(p, type = "dynamic"); pa_rng <- .Random.seed
  expect_identical(pa_rng, sa_rng)
  for (f in c("overall.att", "overall.se", "att.egt", "se.egt", "crit.val.egt")) expect_identical(pa[[f]], sa[[f]], info = f)
})

test_that("multiplier worker buffers are preflighted before any draw, failed workers refuse, and the PSOCK route multiplies identically", {
  X <- rbind(c(-2, 4), c(1, 1), c(-1, 2), c(1, -3), c(1, -2), c(2, 1), c(2, -4), c(-1, 2), c(0.5, 1), c(-3.5, -3), c(0, 0), c(0, 0))
  colnames(X) <- c("a", "b"); th <- c(a = 0.5, b = -0.2); cl <- c(rep(1:5, each = 2), 6, 6)
  budget <- function(workers = NULL, cells = 1e6) { b <- list(max_units = 2000, max_clusters = 2000, max_targets = 100,
    max_draws = 2000, max_if_batch_rows = 64, max_if_batch_cells = 6400, max_sign_batch_rows = 16, max_sign_batch_cells = 64000,
    max_workspace_cells = cells, max_workspace_bytes = 8 * cells); if (!is.null(workers)) b$workers <- workers; b }
  run <- function(b) api$compute_chained_ipw_inference(th, X, cl, 0.05, cband = TRUE, biters = 200L, multiplier_budget = b,
    inference_accuracy = NULL, bstrap = TRUE, caller_rng = TRUE)
  set.seed(11); s <- run(budget()); s_rng <- .Random.seed
  set.seed(11); p <- run(budget(workers = 2)); p_rng <- .Random.seed
  expect_identical(p_rng, s_rng); expect_identical(p$multiplier$raw_draws, s$multiplier$raw_draws)
  expect_identical(p$cband_crit_val, s$cband_crit_val); expect_identical(p$se, s$se)
  expect_identical(p$multiplier$metadata$n_zero_sum_clusters, 1L)          # zero-IF cluster keeps its signs
  serial_cells <- 4 * 12 * 2 + 4 * 6 * 2 + 4 * 2 * 2 + 4 * 200 * 2 + 4 * 16 * 6
  set.seed(5); before <- .Random.seed
  expect_error(run(budget(workers = 4, cells = serial_cells)), class = "cdid_ipw_pipeline_bound_error")
  expect_identical(.Random.seed, before)
  # The PSOCK route (Windows) runs here too: two R worker processes for these products only.
  cl2 <- parallel::makePSOCKcluster(2L); on.exit(parallel::stopCluster(cl2), add = TRUE)
  S <- matrix(c(1.5, -2, 0.25, 3, 0, -1), 3, 2); sm <- list(matrix(c(1, -1, 1, 1, -1, -1), 2, 3), matrix(c(-1, 1, 1), 1, 3))
  expect_identical(api$chained_multiplier_parallel_products(sm, S, 1 / sqrt(3), 2, route = "psock", cluster = cl2),
    lapply(sm, function(m) (m %*% S) * (1 / sqrt(3))))
  failing <- list(matrix(1, 2, 3), matrix("a", 1, 3))
  expect_error(api$chained_multiplier_parallel_products(failing, matrix(1, 3, 2), 1, 2, route = "psock", cluster = cl2),
    class = "cdid_multiplier_parallel_error")
  if (.Platform$OS.type != "windows")
    expect_error(suppressWarnings(api$chained_multiplier_parallel_products(failing, matrix(1, 3, 2), 1, 2, route = "fork")),
      class = "cdid_multiplier_parallel_error")
})

test_that("print_details = TRUE reports progress messages; FALSE is quiet; estimates and RNG are unchanged", {
  skip_if_not_installed("did")
  d <- cp8_public_panel(); msgs <- character()
  set.seed(1); a <- withCallingHandlers(cp8_fit(d, est_method = "ipw", biters = 49, print_details = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }); a_rng <- .Random.seed
  quiet <- character()
  set.seed(1); b <- withCallingHandlers(cp8_fit(d, est_method = "ipw", biters = 49),
    message = function(m) { quiet <<- c(quiet, conditionMessage(m)); invokeRestart("muffleMessage") }); b_rng <- .Random.seed
  expect_gte(length(msgs), 8L); expect_length(quiet, 0L)
  expect_true(any(grepl("comparison catalog", msgs))); expect_true(any(grepl("comparisons", msgs)))
  expect_identical(a_rng, b_rng)
  for (f in c("att", "se", "c", "inffunc")) expect_identical(a[[f]], b[[f]], info = f)
  expect_true(a$DIDparams$print_details); expect_false(b$DIDparams$print_details); expect_false(a$DIDparams$faster_mode)
})
