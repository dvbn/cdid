# SPDX-License-Identifier: GPL-3.0-only
#
# Focused CP08 P1 tests for `build_pipeline_nuisance_and_link_views`: an
# outcome-free mask fingerprint rejects unequal observation-mask candidates
# before full comparison construction; only the exact streamed comparison
# accepts reuse. Strict installed namespace, no source fallback. Expectations
# come from references built here from raw panel rows, or from a plain
# per-candidate replay through the unchanged validated helper
# `construct_pipeline_catalog_row_views`; the new grouping never supplies them.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api), environmentIsLocked(api))

p1_reason <- "Link risk-set pair-observation mask differs from the nuisance view mask on chunk; reuse is refused."
p1_P <- 67108859 # 2^26 - 5, the declared fingerprint modulus

p1_controls <- function(chunk = 64L, max_retained = 2e6, max_matrix = 2e6, max_work = 5e8) {
  list(
    owner = api$ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = max_retained,
      max_matrix_entries = max_matrix, max_work_ops = max_work, max_output_entries = 5e5),
    fitted = api$ipw_fitted_view_control(max_chunk_rows = chunk, max_observer_units = 2000L,
      max_retained_entries = 2e6, max_matrix_entries = 2e6, max_work_ops = max_work, max_output_entries = 5e5),
    point = api$ipw_point_view_control(max_chunk_rows = chunk, max_observer_units = 2000L,
      max_retained_entries = 2e6, max_matrix_entries = 2e6, max_work_ops = max_work, max_output_entries = 5e5),
    glm = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    bounds = list(max_candidates = 1e4, max_rows = 1e4, max_pool_entries = 1e5, max_scratch = 1e7, max_work = max_work)
  )
}

# Long panel from explicit per-unit observed dates; deterministic, no draws.
p1_panel <- function(obs_dates, g) {
  do.call(rbind, lapply(seq_along(obs_dates), function(i) {
    tt <- as.double(obs_dates[[i]])
    data.frame(id = i, t = tt, g = g[i], x = ((5 * i) %% 7) / 7,
      y = 0.5 * i + 0.25 * tt + ifelse(tt >= g[i], 1.5, 0) + ((7 * i + 3 * tt) %% 5) / 10)
  }))
}

# Positions {1,5,8,12} and {2,3,10,11} share count and the first three power
# sums (a Prouhet-Tarry-Escott pair): masks of rows (1,2) and (1,3) collide.
# {4,6,7,9} matches their count and first power sum only. Units 13-24 are
# balanced and carry both cohorts on every pair.
p1_fixture_pte <- function() {
  obs <- vector("list", 24L)
  for (i in c(1, 5, 8, 12)) obs[[i]] <- c(1, 2)
  for (i in c(2, 3, 10, 11)) obs[[i]] <- c(1, 3)
  for (i in c(4, 6, 7, 9)) obs[[i]] <- c(2, 4)
  for (i in 13:24) obs[[i]] <- 1:4
  p1_panel(obs, ifelse(seq_len(24) %% 2 == 1, 4, Inf))
}

# Units 19 and 20 (same cohort, last A positions) alone distinguish the masks.
p1_fixture_late <- function() {
  obs <- c(rep(list(1:4), 18L), list(c(1, 2, 3)), list(c(2, 3, 4)))
  p1_panel(obs, c(ifelse(seq_len(18) %% 2 == 1, 4, Inf), 4, 4))
}

# Three cohorts and several operators; one unit per cohort misses date 1 and
# another misses date 2, giving equal-count unequal masks and shared masks.
p1_fixture_multi <- function() {
  g <- c(4, 5, Inf)[(seq_len(30) %% 3) + 1]
  obs <- lapply(seq_len(30), function(i) if (i <= 3) 2:5 else if (i <= 6) c(1, 3:5) else 1:5)
  p1_panel(obs, g)
}

p1_build <- function(df, targets = 4, control_group = "nevertreated", chunk = 64L,
                     intercept_only = FALSE, weights = NULL, ctrl = p1_controls(chunk)) {
  po <- api$construct_ipw_panel_owner(df, "id", "t", "y", "g", ctrl$owner)
  u <- df[!duplicated(df$id), ]
  u <- u[order(u$id), ]
  Z <- if (intercept_only) matrix(1, nrow(u), 1L) else cbind(1, u$x)
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = weights, control = ctrl$owner)
  meta <- api$build_owner_comparison_catalog(po, targets = targets, anticipation = 0,
    control_group = control_group, conventions = "p1", control = ctrl$owner,
    catalog_bounds = ctrl$bounds, frame_owner = fo)
  list(df = df, po = po, fo = fo, meta = meta, ctrl = ctrl)
}

p1_construct <- function(fx, ctrl = fx$ctrl) {
  api$build_pipeline_nuisance_and_link_views(fx$meta, fx$fo, ctrl$owner, ctrl$fitted, ctrl$point, ctrl$glm)
}

p1_ref_fingerprint <- function(pos) {
  u <- as.double(pos) %% p1_P
  u2 <- (u * u) %% p1_P
  c(length(pos), sum(u) %% p1_P, sum(u2) %% p1_P, sum((u2 * u) %% p1_P) %% p1_P)
}

# Reference from raw rows: A-units have >= 2 non-missing outcomes, ordered by
# id; a row mask is the sorted A positions with cohort in {g} U K observed at
# both s and t. Rows are grouped per operator in catalog order, first exact
# match wins; each earlier unequal representative is one split entry.
p1_reference <- function(df, meta) {
  obs <- df[!is.na(df$y), ]
  ids <- sort(unique(df$id))
  a_ids <- ids[tabulate(match(obs$id, ids), nbins = length(ids)) >= 2L]
  a_g <- df$g[match(a_ids, df$id)]
  seen <- function(tt) a_ids %in% obs$id[obs$t == tt]
  rows <- meta$catalog$rows
  pl <- meta$catalog$pools
  masks <- lapply(seq_len(nrow(rows)), function(r) {
    K <- meta$catalog$spec$cohorts[pl$entries[seq.int(pl$start[rows$pool_id[r]], length.out = pl$size[rows$pool_id[r]])]]
    which(a_g %in% c(rows$g[r], K) & seen(rows$s[r]) & seen(rows$t[r]))
  })
  grp_mask <- list(); grp_op <- integer(); grp_row <- integer(); keys <- character()
  l2n <- integer(nrow(rows)); splits <- list(); collisions <- 0; matches <- 0
  for (r in seq_len(nrow(rows))) {
    op <- as.integer(rows$operator_id[r])
    hit <- NA_integer_
    for (k in which(grp_op == op)) {
      if (identical(grp_mask[[k]], masks[[r]])) { hit <- k; matches <- matches + 1; break }
      if (identical(p1_ref_fingerprint(grp_mask[[k]]), p1_ref_fingerprint(masks[[r]]))) collisions <- collisions + 1
      splits[[length(splits) + 1L]] <- list(row_id = r, operator_id = op,
        candidate_nuisance_id = paste0("nuisance:op", op, ":r", grp_row[k]), reason = p1_reason)
    }
    if (is.na(hit)) {
      hit <- length(grp_mask) + 1L
      grp_mask[[hit]] <- masks[[r]]; grp_op[hit] <- op; grp_row[hit] <- r
      keys[hit] <- paste0("op", op, ":mask", sum(grp_op == op))
    }
    l2n[r] <- hit
  }
  list(masks = masks, a_g = a_g, l2n = l2n, keys = keys, ops = grp_op, splits = splits,
    collisions = collisions, matches = matches, comparisons = length(splits) + matches)
}

# Straightforward per-candidate replay through the unchanged validated helper.
p1_replay <- function(fx) {
  rows <- fx$meta$catalog$rows
  views <- list(); ops <- integer(); l2n <- integer(nrow(rows)); links <- vector("list", nrow(rows))
  for (m in seq_len(nrow(rows))) {
    hit <- NA_integer_
    for (k in which(ops == rows$operator_id[m])) {
      a <- tryCatch(api$construct_pipeline_catalog_row_views(fx$meta, m, control = fx$ctrl$owner,
        frame_owner = fx$fo, nuisance_view = views[[k]]), cdid_ipw_shared_mask_error = function(e) NULL)
      if (!is.null(a)) { hit <- k; links[[m]] <- a$link_view; break }
    }
    if (is.na(hit)) {
      a <- api$construct_pipeline_catalog_row_views(fx$meta, m, control = fx$ctrl$owner, frame_owner = fx$fo)
      views[[length(views) + 1L]] <- a$nuisance_view; ops <- c(ops, rows$operator_id[m])
      hit <- length(views); links[[m]] <- a$link_view
    }
    l2n[m] <- hit
  }
  list(views = views, links = links, l2n = l2n)
}

# Counting wrappers that delegate to the real binding and are always restored
# (the instrumentation pattern of the metadata tests).
p1_instrument <- function(fns, run, snapshot_at = NULL) {
  calls <- stats::setNames(integer(length(fns)), fns)
  snapshot <- NULL
  originals <- lapply(stats::setNames(fns, fns), function(nm) get(nm, envir = api, inherits = FALSE))
  rebind <- function(nm, f) {
    locked <- bindingIsLocked(nm, api)
    if (locked) unlockBinding(nm, api)
    assign(nm, f, envir = api)
    if (locked) lockBinding(nm, api)
  }
  on.exit(for (nm in fns) rebind(nm, originals[[nm]]), add = TRUE)
  for (nm in fns) local({
    nm <- nm
    orig <- originals[[nm]]
    rebind(nm, function(...) {
      calls[[nm]] <<- calls[[nm]] + 1L
      if (identical(nm, snapshot_at) && is.null(snapshot)) snapshot <<- calls
      orig(...)
    })
  })
  value <- run()
  list(value = value, calls = calls, snapshot = snapshot)
}

# The real interrupt check always runs first; while a fingerprint stream is
# active the wrapper then raises a simulated interrupt condition once.
p1_with_fingerprint_interrupt <- function(run) {
  fp_orig <- get("fingerprint_ipw_view_mask", envir = api, inherits = FALSE)
  ck_orig <- get("ipw_owner_check_interrupt", envir = api, inherits = FALSE)
  in_fp <- FALSE
  raised <- 0L
  rebind <- function(nm, f) {
    locked <- bindingIsLocked(nm, api)
    if (locked) unlockBinding(nm, api)
    assign(nm, f, envir = api)
    if (locked) lockBinding(nm, api)
  }
  on.exit({ rebind("fingerprint_ipw_view_mask", fp_orig); rebind("ipw_owner_check_interrupt", ck_orig) }, add = TRUE)
  rebind("fingerprint_ipw_view_mask", function(...) { in_fp <<- TRUE; on.exit(in_fp <<- FALSE); fp_orig(...) })
  rebind("ipw_owner_check_interrupt", function() {
    ck_orig()
    if (in_fp) {
      raised <<- raised + 1L
      stop(structure(class = c("cdid_p1_simulated_interrupt", "interrupt", "error", "condition"),
        list(message = "simulated user interrupt at a fingerprint chunk", call = NULL)))
    }
    invisible(NULL)
  })
  value <- run()
  list(value = value, raised = raised, restored = list(fp = fp_orig, check = ck_orig))
}

p1_strip_env <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.list(x)) { at <- attributes(x); x <- lapply(x, p1_strip_env); attributes(x) <- at }
  x
}

p1_link_if <- function(point_links, ctrl) {
  out <- list()
  api$replay_ipw_point_if_view(point_links = unname(point_links), control = ctrl, callback_schema = TRUE,
    callback = function(batch) out[[length(out) + 1L]] <<- batch$phi)
  do.call(rbind, out)
}

test_that("equal-power-sum masks collide in the fingerprint and are separated only by the exact comparison", {
  fx <- p1_build(p1_fixture_pte())
  rows <- fx$meta$catalog$rows
  M <- nrow(rows)
  ref <- p1_reference(fx$df, fx$meta)
  r12 <- which(rows$s == 1 & rows$t == 2); r13 <- which(rows$s == 1 & rows$t == 3); r24 <- which(rows$s == 2 & rows$t == 4)
  expect_identical(M, 6L)
  expect_false(identical(ref$masks[[r12]], ref$masks[[r13]]))
  expect_identical(p1_ref_fingerprint(ref$masks[[r12]]), p1_ref_fingerprint(ref$masks[[r13]]))
  expect_identical(length(ref$masks[[r24]]), length(ref$masks[[r12]]))
  expect_identical(sum(ref$masks[[r24]]), sum(ref$masks[[r12]]))
  expect_identical(ref$collisions, 1)

  run <- p1_instrument(c("construct_ipw_link_view", "construct_ipw_nuisance_view"), function() p1_construct(fx))
  res <- run$value
  expect_identical(res$link_to_nuisance_idx, ref$l2n)
  expect_identical(res$nuisance_keys, ref$keys)
  expect_identical(res$operator_ids, ref$ops)
  expect_identical(res$link_nuisance_keys, ref$keys[ref$l2n])
  expect_identical(res$mask_splits, ref$splits)
  expect_identical(unname(run$calls[["construct_ipw_link_view"]]), M + as.integer(ref$collisions))
  expect_identical(unname(run$calls[["construct_ipw_nuisance_view"]]), length(ref$keys))
  w <- res$construction_work
  expect_identical(w$exact_mask_mismatches, ref$collisions)
  expect_identical(w$exact_candidate_attempts, ref$collisions + ref$matches)
  expect_identical(w$candidate_comparisons, ref$comparisons)
  expect_identical(w$fingerprint_rejections, ref$comparisons - ref$collisions - ref$matches)
  expect_identical(w$endpoints_read, 0L)
  expect_identical(w$endpoints_read_scope, "mask_screening_only")

  # The production value is the declared function and genuinely collides.
  fp12 <- api$fingerprint_ipw_view_mask(res$nuisance_views[[ref$l2n[r12]]], fx$ctrl$owner)
  fp13 <- api$fingerprint_ipw_view_mask(res$nuisance_views[[ref$l2n[r13]]], fx$ctrl$owner)
  expect_identical(unname(fp12$fingerprint), p1_ref_fingerprint(ref$masks[[r12]]))
  expect_identical(fp12$fingerprint, fp13$fingerprint)
  expect_identical(fp12$endpoints_read, 0L)
  expect_error(api$verify_ipw_view_mask_compatibility(res$nuisance_views[[ref$l2n[r12]]],
    res$nuisance_views[[ref$l2n[r13]]], fx$ctrl$owner), class = "cdid_ipw_shared_mask_error")
})

test_that("late-chunk differences with equal per-cohort counts are rejected without building a candidate link", {
  fx <- p1_build(p1_fixture_late(), chunk = 4L)
  rows <- fx$meta$catalog$rows
  M <- nrow(rows)
  ref <- p1_reference(fx$df, fx$meta)
  r12 <- which(rows$s == 1 & rows$t == 2); r24 <- which(rows$s == 2 & rows$t == 4)
  a <- ref$masks[[r12]]; b <- ref$masks[[r24]]
  expect_identical(length(a), length(b))
  expect_identical(table(ref$a_g[a]), table(ref$a_g[b]))
  expect_true(all(c(setdiff(a, b), setdiff(b, a)) > 16L)) # only the fifth 4-row chunk differs
  expect_false(identical(ref$l2n[r12], ref$l2n[r24]))

  run <- p1_instrument(c("construct_ipw_link_view", "construct_ipw_nuisance_view"), function() p1_construct(fx))
  res <- run$value
  expect_identical(ref$collisions, 0)
  expect_identical(res$link_to_nuisance_idx, ref$l2n)
  expect_identical(res$nuisance_keys, ref$keys)
  expect_identical(res$mask_splits, ref$splits)
  expect_identical(unname(run$calls[["construct_ipw_link_view"]]), M)
  expect_identical(unname(run$calls[["construct_ipw_nuisance_view"]]), length(ref$keys))
  expect_identical(res$construction_work$fingerprint_rejections, ref$comparisons - ref$matches)

  # Discriminating count: the per-candidate procedure builds a link for every
  # failed candidate, so a regression to it fails the expectation above.
  old <- p1_instrument("construct_ipw_link_view", function() p1_replay(fx))
  expect_identical(old$value$l2n, ref$l2n)
  expect_identical(unname(old$calls[["construct_ipw_link_view"]]), M + length(ref$splits))
  expect_gt(length(ref$splits), 0L)
})

test_that("views, links, fits, point estimates and complete link IFs equal the per-candidate replay", {
  for (case in list(list(df = p1_fixture_pte(), chunk = 64L, targets = 4, cg = "nevertreated"),
                    list(df = p1_fixture_late(), chunk = 4L, targets = 4, cg = "nevertreated"),
                    list(df = p1_fixture_multi(), chunk = 5L, targets = c(4, 5), cg = "notyettreated"))) {
    fx <- p1_build(case$df, targets = case$targets, control_group = case$cg, chunk = case$chunk)
    res <- p1_construct(fx)
    rep <- p1_replay(fx)
    expect_identical(res$link_to_nuisance_idx, rep$l2n)
    expect_identical(p1_strip_env(unname(res$nuisance_views)), p1_strip_env(rep$views))
    expect_identical(p1_strip_env(unname(res$link_views)), p1_strip_env(rep$links))
    for (m in seq_along(rep$links)) {
      expect_identical(res$link_views[[m]]$nuisance_ref, res$nuisance_views[[res$link_to_nuisance_idx[m]]]$ref)
      expect_true(isTRUE(res$link_views[[m]]$nuisance_mask_verified))
    }
    fits <- api$fit_ipw_nuisance_view_set(rep$views, fx$ctrl$glm, fx$ctrl$fitted)
    pts <- lapply(seq_along(rep$links), function(m)
      api$fit_ipw_point_link_view(fits[[rep$l2n[m]]], rep$links[[m]], fx$ctrl$point))
    expect_identical(p1_strip_env(unname(res$fitted_views)), p1_strip_env(fits))
    expect_identical(vapply(unname(res$point_links), `[[`, numeric(1L), "estimate"),
      vapply(pts, `[[`, numeric(1L), "estimate"))
    expect_identical(p1_strip_env(unname(res$point_links)), p1_strip_env(pts))
    if_new <- p1_link_if(res$point_links, fx$ctrl$point)
    expect_identical(dim(if_new), c(fx$fo$n_a_units, length(rep$links)))
    expect_identical(unname(if_new), unname(p1_link_if(pts, fx$ctrl$point)))
  }
})

test_that("the trusted loop validates the catalog once per call and reads no outcome while grouping", {
  fx <- p1_build(p1_fixture_multi(), targets = c(4, 5), control_group = "notyettreated", chunk = 5L)
  ref <- p1_reference(fx$df, fx$meta)
  run <- p1_instrument(c("validate_owner_comparison_catalog", "ipw_owner_read_y_positions",
    "fingerprint_ipw_view_mask", "fit_ipw_nuisance_view_set"),
    function() p1_construct(fx), snapshot_at = "fit_ipw_nuisance_view_set")
  expect_identical(run$value$link_to_nuisance_idx, ref$l2n)
  expect_identical(run$value$mask_splits, ref$splits)
  expect_gt(length(unique(ref$ops)), 1L)
  expect_identical(unname(run$calls[["validate_owner_comparison_catalog"]]), 2L)
  expect_identical(unname(run$snapshot[["ipw_owner_read_y_positions"]]), 0L)
  expect_identical(unname(run$snapshot[["fingerprint_ipw_view_mask"]]),
    as.integer(run$value$construction_work$fingerprint_passes))
  expect_gt(unname(run$calls[["ipw_owner_read_y_positions"]]), 0L) # point fits do read endpoints
  expect_identical(run$value$construction_work$endpoints_read_scope, "mask_screening_only")
})

test_that("rows of different operators never share and cohort or owner changes are refused or kept apart", {
  fx <- p1_build(p1_fixture_multi(), targets = c(4, 5), control_group = "notyettreated", chunk = 5L)
  res <- p1_construct(fx)
  ops_of_view <- res$operator_ids[res$link_to_nuisance_idx]
  expect_identical(ops_of_view, as.integer(fx$meta$catalog$rows$operator_id))
  nv <- res$nuisance_views[[1L]]
  other <- res$nuisance_views[[which(!vapply(res$nuisance_views, function(v)
    identical(v$cohorts$owner_cohort_index, nv$cohorts$owner_cohort_index), logical(1L)))[1L]]]
  s <- nv$defining_pair$s; t <- nv$defining_pair$t
  same <- api$fingerprint_ipw_view_mask(nv, fx$ctrl$owner)$fingerprint
  expect_identical(api$fingerprint_ipw_view_mask(nv, fx$ctrl$owner, s = s, t = t)$fingerprint, same)
  expect_false(identical(api$fingerprint_ipw_view_mask(other, fx$ctrl$owner, s = s, t = t)$fingerprint, same))
  expect_error(api$construct_ipw_link_view(fx$fo, s, t, g = other$cohorts$g, K = other$cohorts$K,
    nuisance_view = nv, control = fx$ctrl$owner), class = "cdid_ipw_shared_incompatible_error")

  u <- fx$df[!duplicated(fx$df$id), ]; u <- u[order(u$id), ]
  fo2 <- api$construct_ipw_frame_owner(fx$po, u$id, cbind(1, u$x), weights = rep(2, nrow(u)), control = fx$ctrl$owner)
  fx2 <- fx; fx2$fo <- fo2
  expect_error(p1_construct(fx2), class = "cdid_ipw_owner_mismatch_error")
  expect_error(api$construct_ipw_link_view(fo2, s, t, nuisance_view = nv, control = fx$ctrl$owner),
    class = "cdid_ipw_owner_mismatch_error")
})

test_that("covariates and weights leave the mask grouping unchanged; invalid X or weights are refused first", {
  df <- p1_fixture_pte()
  ref <- NULL
  fits <- list()
  for (cfg in list(list(io = TRUE, w = NULL), list(io = FALSE, w = NULL), list(io = FALSE, w = 1 + (seq_len(24) %% 4) / 4))) {
    fx <- p1_build(df, intercept_only = cfg$io, weights = cfg$w)
    if (is.null(ref)) ref <- p1_reference(df, fx$meta)
    res <- p1_construct(fx)
    expect_identical(res$link_to_nuisance_idx, ref$l2n)
    expect_identical(res$mask_splits, ref$splits)
    expect_identical(res$construction_work$exact_mask_mismatches, ref$collisions)
    fits[[length(fits) + 1L]] <- res$fitted_views[[1L]]$beta
  }
  expect_false(identical(fits[[2L]], fits[[3L]])) # weights reach the fits, not the grouping
  ctrl <- p1_controls()
  po <- api$construct_ipw_panel_owner(df, "id", "t", "y", "g", ctrl$owner)
  u <- df[!duplicated(df$id), ]; u <- u[order(u$id), ]
  expect_error(api$construct_ipw_frame_owner(po, u$id, cbind(1, u$x), weights = c(0, rep(1, 23)), control = ctrl$owner),
    class = "cdid_ipw_owner_weight_error")
  expect_error(api$construct_ipw_frame_owner(po, u$id, cbind(1, c(NA, u$x[-1L])), control = ctrl$owner),
    class = "cdid_ipw_owner_covariate_error")
})

test_that("tampered metadata, frame owners and views are refused before any construction", {
  fx <- p1_build(p1_fixture_pte())
  bad_meta <- fx$meta
  bad_meta$catalog$rows$t[1L] <- 3
  run <- p1_instrument("construct_ipw_nuisance_view", function()
    tryCatch(api$build_pipeline_nuisance_and_link_views(bad_meta, fx$fo, fx$ctrl$owner, fx$ctrl$fitted,
      fx$ctrl$point, fx$ctrl$glm), error = identity))
  expect_s3_class(run$value, "cdid_pipeline_metadata_integrity_error")
  expect_identical(unname(run$calls[["construct_ipw_nuisance_view"]]), 0L)
  expect_error(api$construct_pipeline_catalog_row_views(bad_meta, 1L, control = fx$ctrl$owner, frame_owner = fx$fo),
    class = "cdid_pipeline_metadata_integrity_error")
  expect_error(api$extract_pipeline_catalog_row_spec(bad_meta, 1L), class = "cdid_pipeline_metadata_integrity_error")

  bad_fo <- fx$fo
  bad_fo$w_a[1L] <- 5
  expect_error(api$build_pipeline_nuisance_and_link_views(fx$meta, bad_fo, fx$ctrl$owner, fx$ctrl$fitted,
    fx$ctrl$point, fx$ctrl$glm), class = "cdid_ipw_owner_integrity_error")

  res <- p1_construct(fx)
  bad_nv <- res$nuisance_views[[1L]]
  bad_nv$cohorts$owner_cohort_index <- rev(bad_nv$cohorts$owner_cohort_index)
  expect_error(api$fingerprint_ipw_view_mask(bad_nv, fx$ctrl$owner), class = "cdid_ipw_owner_integrity_error")
  expect_error(api$fingerprint_ipw_view_mask(res$nuisance_views[[1L]], fx$ctrl$owner, s = 1),
    class = "cdid_ipw_owner_input_error")
})

test_that("resource refusals precede streaming and interrupts at fingerprint chunks propagate", {
  fx <- p1_build(p1_fixture_late(), chunk = 4L)
  res <- p1_construct(fx)
  nv <- res$nuisance_views[[1L]]
  tiny <- p1_controls(chunk = 4L, max_matrix = 300)
  run <- p1_instrument("ipw_owner_query_pair_mask", function()
    tryCatch(api$fingerprint_ipw_view_mask(nv, tiny$owner), error = identity))
  expect_s3_class(run$value, "cdid_ipw_owner_bound_error")
  expect_match(conditionMessage(run$value), "Streamed mask fingerprint exceeds", fixed = TRUE)
  expect_identical(unname(run$calls[["ipw_owner_query_pair_mask"]]), 0L)

  M <- nrow(fx$meta$catalog$rows)
  low <- p1_controls(chunk = 4L, max_retained = 16 * (M + 1) + length(fx$meta$catalog$spec$cohorts) - 1)
  run <- p1_instrument("construct_ipw_nuisance_view", function()
    tryCatch(api$build_pipeline_nuisance_and_link_views(fx$meta, fx$fo, low$owner, fx$ctrl$fitted,
      fx$ctrl$point, fx$ctrl$glm), error = identity))
  expect_s3_class(run$value, "cdid_ipw_owner_bound_error")
  expect_match(conditionMessage(run$value), "invocation-local mask fingerprint registry", fixed = TRUE)
  expect_identical(unname(run$calls[["construct_ipw_nuisance_view"]]), 0L)

  # One interrupt check per streamed chunk.
  n_chunks <- ceiling(fx$fo$n_a_units / 4)
  run <- p1_instrument("ipw_owner_check_interrupt", function() api$fingerprint_ipw_view_mask(nv, fx$ctrl$owner))
  expect_identical(unname(run$calls[["ipw_owner_check_interrupt"]]), as.integer(n_chunks))
  expect_identical(run$value$chunks_streamed, n_chunks)

  # A condition raised at a fingerprint checkpoint (after the real check runs)
  # unwinds the whole construction call; the loop only absorbs mask refusals.
  out <- p1_with_fingerprint_interrupt(function() tryCatch(p1_construct(fx), error = identity))
  expect_s3_class(out$value, "cdid_p1_simulated_interrupt")
  expect_identical(out$raised, 1L)
  expect_identical(get("fingerprint_ipw_view_mask", envir = api), out$restored$fp)
  expect_identical(get("ipw_owner_check_interrupt", envir = api), out$restored$check)
})

test_that("fingerprint scratch is composed with the live retained reservation before querying", {
  fx <- p1_build(p1_fixture_late(), chunk = 4L)
  spec <- api$extract_pipeline_catalog_row_spec(fx$meta, 1L)
  nv <- api$construct_pipeline_catalog_row_views(fx$meta, 1L,
    control = fx$ctrl$owner, frame_owner = fx$fo)$nuisance_view
  ctrl <- p1_controls(chunk = 4L, max_matrix = 2000)$owner
  # Each allocation separately fits: retained 1500 entries and at least 512
  # entries of chunk/cohort/header scratch. Their simultaneous sum does not.
  run <- p1_instrument("ipw_owner_query_pair_mask", function()
    tryCatch(api$fingerprint_ipw_view_mask(nv, ctrl, s = spec$s, t = spec$t,
      retained_entries = 1500), error = identity))
  expect_s3_class(run$value, "cdid_ipw_owner_bound_error")
  expect_match(conditionMessage(run$value), "concurrent mask fingerprint workspace", fixed = TRUE)
  expect_identical(unname(run$calls[["ipw_owner_query_pair_mask"]]), 0L)
})

test_that("the mask registry shares the existing constructor and whole-pipeline envelopes", {
  fx <- p1_build(p1_fixture_late(), chunk = 4L)
  spec <- api$extract_pipeline_catalog_row_spec(fx$meta, 1L)
  cohorts <- 1 + length(spec$K)
  nu <- api$preflight_ipw_nuisance_view(fx$fo$n_a_units, fx$fo$n_covariates, cohorts, fx$ctrl$owner)
  lk <- api$preflight_ipw_link_view(fx$fo$n_a_units, fx$fo$n_covariates, cohorts, fx$ctrl$owner)
  # Row-only caps omit the live registry and descriptor context. The composed
  # descriptor preflight now refuses before any row constructor or mask query.
  cases <- list(
    list(ctrl = p1_controls(chunk = 4L, max_matrix = max(nu$concurrent_entries, lk$concurrent_entries)),
      reason = "caller plus descriptor context and scratch"),
    list(ctrl = p1_controls(chunk = 4L, max_retained = max(nu$total_retained_entries, lk$total_retained_entries)),
      reason = "descriptor context retained entries"))
  for (case in cases) {
    ctrl <- case$ctrl
    run <- p1_instrument(c("construct_ipw_nuisance_view", "ipw_owner_query_pair_mask"), function()
      tryCatch(api$build_pipeline_nuisance_and_link_views(fx$meta, fx$fo, ctrl$owner,
        ctrl$fitted, ctrl$point, ctrl$glm), error = identity))
    expect_s3_class(run$value, "cdid_ipw_owner_bound_error")
    expect_match(conditionMessage(run$value), case$reason, fixed = TRUE)
    expect_true(all(run$calls == 0L))
  }
  run <- p1_instrument(c("construct_ipw_nuisance_view", "ipw_owner_query_pair_mask"), function()
    tryCatch(api$build_pipeline_nuisance_and_link_views(fx$meta, fx$fo, fx$ctrl$owner,
      fx$ctrl$fitted, fx$ctrl$point, fx$ctrl$glm,
      retained_pipeline_entries = fx$ctrl$owner$max_matrix_entries), error = identity))
  expect_s3_class(run$value, "cdid_ipw_owner_bound_error")
  expect_match(conditionMessage(run$value), "pipeline envelope plus invocation-local mask fingerprint registry", fixed = TRUE)
  expect_true(all(run$calls == 0L))
})

test_that("the complete pipeline retains the reference grouping and the construction counts", {
  df <- p1_fixture_multi()
  lnk <- api$ipw_link_enclosure_control(max_chunk_rows = 5L, max_retained_entries = 2e6, max_matrix_entries = 2e6,
    max_work_ops = 5e9, max_output_entries = 5e5)
  ctrl <- api$ipw_pipeline_control(
    owner_control = api$ipw_owner_control(max_chunk_rows = 5L, max_retained_entries = 2e6, max_matrix_entries = 2e6,
      max_work_ops = 5e9, max_output_entries = 5e5),
    catalog_bounds = list(max_candidates = 500L, max_rows = 500L, max_pool_entries = 5000L, max_scratch = 1e5, max_work = 5e9),
    class_bounds = list(max_links = 500L, max_classes = 50L, max_class_dates = 50L, max_E_entries = 5000L, max_mask_entries = 1e5),
    fitted_control = api$ipw_fitted_view_control(max_chunk_rows = 5L, max_observer_units = 2000L, max_retained_entries = 2e6,
      max_matrix_entries = 2e6, max_work_ops = 5e9, max_output_entries = 5e5),
    point_control = api$ipw_point_view_control(max_chunk_rows = 5L, max_observer_units = 2000L, max_retained_entries = 2e6,
      max_matrix_entries = 2e6, max_work_ops = 5e9, max_output_entries = 5e5),
    glm_control = list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
    certificate_control = api$logit_certificate_control(taylor_degree = 12L, max_abs_eta = 50, max_squarings = 16L,
      max_chunk_rows = 5L, max_matrix_entries = 2e6, max_work_ops = 5e9),
    problem_control = api$ipw_nuisance_problem_control(max_chunk_rows = 5L, max_retained_entries = 2e6,
      max_matrix_entries = 2e6, max_work_ops = 5e9),
    link_control = lnk,
    covariance_control = api$ipw_covariance_transfer_control(max_chunk_rows = 5L, max_operator_block_rows = 5L,
      max_cluster_rows = 2000L, max_matrix_entries = 2e6, max_work_ops = 5e9, max_output_entries = 1e6, link_control = lnk),
    gmm_accuracy = api$gmm_accuracy_control(method = "certified_eigen_chol_kkt", target_error_unit = "abs",
      max_target_error_abs = 1e-4, max_target_error_se = 1e-2, max_orthogonality_defect = 1e-8,
      min_certified_lambda_rel = 1e-10, max_factor_residual_rel = 1e-10, min_rcond_R = 1e-10, min_rcond_U = 1e-10,
      min_sigma_rel_A = 1e-10, max_kkt_residual_rel = 1e-10, max_loading_residual_rel = 1e-10,
      require_matching_scale = FALSE, max_loading_error_abs = 1e-4),
    gmm_bounds = api$gmm_bounds_control(max_matrix_entries = 2e6, max_batch_rows = 5L, max_targets = 100L,
      max_open_clusters = 2000L, allow_materialize_if = TRUE, max_materialized_if_entries = 1e6),
    multiplier_budget = list(max_units = 2000, max_clusters = 2000, max_targets = 100, max_draws = 2000,
      max_if_batch_rows = 64, max_if_batch_cells = 6400, max_sign_batch_rows = 32, max_sign_batch_cells = 64000,
      max_workspace_cells = 1e6, max_workspace_bytes = 8e6),
    inference_accuracy = list(max_estimate_error_abs = 1e-3, max_se_error_abs = 1e-2, max_vcov_error_abs = 1e-2,
      max_multiplier_error_abs = 1e-1, max_band_error_abs = 1e-1))
  out <- api$run_chained_ipw_pipeline(data = df, idname = "id", tname = "t", yname = "y", gname = "g",
    control = ctrl, xformla = ~ x, target_cohorts = c(4, 5), control_group = "notyettreated",
    base_period = "varying", weighting = "identity", run_numerical_transfer = FALSE, aggregations = "simple")
  reg <- out$nuisance_registry
  ref <- p1_reference(df, out$metadata)
  expect_identical(reg$link_to_nuisance_idx, ref$l2n)
  expect_identical(reg$nuisance_keys, ref$keys)
  expect_identical(reg$mask_splits, ref$splits)
  expect_identical(reg$construction_work$n_rows, nrow(out$metadata$catalog$rows))
  expect_identical(reg$construction_work$exact_mask_mismatches, ref$collisions)
  expect_true(all(is.finite(out$reporting$estimates[out$reporting$available])))
  expect_identical(nrow(out$reporting$if_full), length(unique(df$id)))
})
