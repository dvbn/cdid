# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 balanced-panel performance candidate: exact complete-history posting
# facts, the vectorized sorted lookup, construction reuse (stage 5: also by the
# nuisance-free CP04 registry views, below), validated chunk queries, the IF
# projection duplicate-key check and the payload budget. Every
# expectation comes from a reference written here (raw panel rows, a plain
# scalar binary search, plain checked sums) or from the unchanged streamed
# path of the same build; the fast paths never supply their own truth.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

fp_ctrl <- function(chunk = 5L) api$ipw_owner_control(max_chunk_rows = chunk, max_retained_entries = 1e6,
  max_matrix_entries = 1e6, max_work_ops = 1e9, max_output_entries = 1e6)

# Deterministic long panel; g = Inf codes never treated.
fp_panel <- function(N = 21L, T = 5L, drop = NULL, na = NULL) {
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- ifelse(d$id %% 3 == 0, 3, Inf)
  d$y <- 0.3 * d$id + 0.2 * d$t + ifelse(d$t >= d$g, 1, 0) + ((7 * d$id + 3 * d$t) %% 5) / 10
  if (!is.null(drop)) d <- d[!(paste(d$id, d$t) %in% drop), ]
  if (!is.null(na)) d$y[paste(d$id, d$t) %in% na] <- NA
  d
}

fp_owners <- function(d, chunk = 5L, weights = NULL) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", fp_ctrl(chunk))
  u <- d[!duplicated(d$id), ]; u <- u[order(u$id), ]
  w <- if (is.null(weights)) NULL else weights(u$id)
  fo <- api$construct_ipw_frame_owner(po, u$id, matrix(1, nrow(u), 1L), weights = w, control = fp_ctrl(chunk))
  list(d = d, po = po, fo = fo)
}

# Raw-row reference: the A-units (>= 2 observed outcomes) observed at both dates
# whose cohort is declared, as sorted A positions.
fp_ref_mask <- function(d, s, t, cohorts) {
  obs <- d[!is.na(d$y), ]
  ids <- sort(unique(d$id))
  a_ids <- ids[tabulate(match(obs$id, ids), nbins = length(ids)) >= 2L]
  g <- d$g[match(a_ids, d$id)]
  which(g %in% cohorts & a_ids %in% obs$id[obs$t == s] & a_ids %in% obs$id[obs$t == t])
}

fp_scalar_search <- function(x, v) {
  vapply(x, function(z) {
    lo <- 1; hi <- length(v)
    while (lo <= hi) {
      mid <- lo + floor((hi - lo) / 2)
      if (v[mid] == z) return(as.integer(mid))
      if (v[mid] < z) lo <- mid + 1 else hi <- mid - 1
    }
    NA_integer_
  }, integer(1L))
}

test_that("the vectorized sorted lookup returns the scalar search positions and keeps it elsewhere", {
  tables <- list(1:9, c(-2.5, 0, 1e-300, 7, Inf), c(2L, 4L, 8L, 16L), c(-Inf, -1, 3))
  for (v in tables) {
    x <- c(v, v[1], 0, -0, 5, 1e300)[c(3, 1, 2, 5, 4, 6, 7, 8, 9)[seq_len(min(9L, length(v) + 5L))]]
    x <- x[!is.na(x)]
    x <- c(x, rev(x))
    expect_true(api$ipw_owner_sorted_lookup_vectorizable(x, v, 1, length(v)))
    expect_identical(api$ipw_owner_match_sorted_double(x, v), fp_scalar_search(x, v))
  }
  # Duplicated, unsorted, NA-bearing, signed-zero or longer tables and bounded
  # searches are not vectorized; they keep the scalar search's exact answers.
  dup <- c(1L, 2L, 2L, 2L, 3L)
  expect_false(api$ipw_owner_sorted_lookup_vectorizable(c(2L, 2L, 2L, 2L, 2L, 2L), dup, 1, 5L))
  expect_identical(api$ipw_owner_match_sorted_double(rep(2L, 6L), dup), rep(3L, 6L))
  expect_false(api$ipw_owner_sorted_lookup_vectorizable(c(1, 2, 3), c(-0, 0, 1), 1, 3L))
  expect_false(api$ipw_owner_sorted_lookup_vectorizable(c(3L, 1L, 2L), c(3L, 1L, 2L), 1, 3L))
  expect_false(api$ipw_owner_sorted_lookup_vectorizable(c(1L, NA, 2L), 1:3, 1, 3L))
  expect_false(api$ipw_owner_sorted_lookup_vectorizable(2L, 1:3, 1, 3L))
  expect_false(api$ipw_owner_sorted_lookup_vectorizable(1:6, 1:6, 2, 6L))
  expect_error(api$ipw_owner_match_sorted_double(c(1L, NA, 2L), 1:3))
  expect_identical(api$ipw_owner_match_sorted_double(integer(0L), 1:3), integer(0L))
})

test_that("complete-history posting facts equal the raw rows, and forged postings are not complete", {
  cases <- list(balanced = fp_panel(), missing = fp_panel(drop = "4 2"), explicit_na = fp_panel(na = "9 4"),
    singleton = rbind(fp_panel(), data.frame(id = 22, t = 3, g = Inf, y = 1)))
  for (d in cases) {
    o <- fp_owners(d)
    obs <- d[!is.na(d$y), ]
    a_ids <- sort(unique(obs$id)); a_ids <- a_ids[tabulate(match(obs$id, a_ids)) >= 2L]
    for (k in seq_along(o$po$all_times)) {
      truth <- all(a_ids %in% obs$id[obs$t == o$po$all_times[k]])
      for (step in list(NULL, 1, 4, 100)) {
        expect_identical(api$ipw_owner_posting_date_complete(o$po$postings, o$po$n_a_units, k, step), truth)
      }
    }
  }
  o <- fp_owners(fp_panel())
  p <- o$po$postings
  i <- p$date_start[2]
  p_swap <- p; p_swap$post_a_index[c(i, i + 1L)] <- p$post_a_index[c(i + 1L, i)]
  p_dup <- p; p_dup$post_a_index[i + 1L] <- p$post_a_index[i]
  for (step in list(NULL, 1, 4)) {
    expect_false(api$ipw_owner_posting_date_complete(p_swap, o$po$n_a_units, 2L, step))
    expect_false(api$ipw_owner_posting_date_complete(p_dup, o$po$n_a_units, 2L, step))
  }
  expect_false(api$ipw_owner_posting_date_complete(p, o$po$n_a_units, NA_integer_))
  # Complete-pair hits are the binary-search hits on every chunk.
  for (a in seq_len(o$po$n_a_units)) for (b in unique(c(a, min(o$po$n_a_units, a + 3L), o$po$n_a_units))) {
    expect_identical(api$ipw_owner_query_complete_pair_mask(o$po, 2L, 5L, a, b),
      api$ipw_owner_query_pair_mask(o$po, 2, 5, a, b, ds = 2L, dt = 5L))
  }
})

test_that("mask verification decided by posting facts agrees with the raw-row masks, refusals included", {
  for (d in list(fp_panel(), fp_panel(drop = c("4 2", "8 4")), fp_panel(na = "6 1"))) {
    o <- fp_owners(d)
    tt <- o$po$all_times; pr <- utils::combn(length(tt), 2)
    nvs <- lapply(seq_len(ncol(pr)), function(k) api$construct_ipw_nuisance_view(o$fo, 3, Inf, tt[pr[1, k]], tt[pr[2, k]],
      fp_ctrl(), require_support = FALSE))
    lvs <- lapply(seq_len(ncol(pr)), function(k) api$construct_ipw_link_view(o$fo, tt[pr[1, k]], tt[pr[2, k]],
      nuisance_view = nvs[[1L]], control = fp_ctrl(), require_support = FALSE, verify_nuisance_mask = FALSE))
    for (i in seq_along(nvs)) for (j in seq_along(lvs)) {
      same <- identical(fp_ref_mask(d, tt[pr[1, i]], tt[pr[2, i]], c(3, Inf)), fp_ref_mask(d, tt[pr[1, j]], tt[pr[2, j]], c(3, Inf)))
      res <- tryCatch(api$verify_ipw_view_mask_compatibility(nvs[[i]], lvs[[j]], fp_ctrl()), cdid_ipw_shared_mask_error = function(e) e)
      expect_identical(!inherits(res, "error"), same)
      if (!inherits(res, "error")) expect_identical(res$chunks_compared, as.integer(ceiling(o$fo$n_a_units / 5)))
    }
  }
})

test_that("the trusted registry reuses exact scans only, with views equal to standalone construction", {
  for (case in list(list(d = fp_panel(), reuse_all = TRUE), list(d = fp_panel(drop = "4 2"), reuse_all = FALSE))) {
    o <- fp_owners(case$d, weights = function(id) 1 + (id %% 4) / 4)
    meta <- api$build_owner_comparison_catalog(o$po, targets = 3, anticipation = 0, control_group = "nevertreated",
      conventions = "cp08-fastpath", control = fp_ctrl(),
      catalog_bounds = list(max_candidates = 1e4, max_rows = 1e4, max_pool_entries = 1e5, max_scratch = 1e7, max_work = 1e9),
      frame_owner = o$fo)
    fc <- api$ipw_fitted_view_control(max_chunk_rows = 5L, max_observer_units = 1000L, max_retained_entries = 1e6,
      max_matrix_entries = 1e6, max_work_ops = 1e9, max_output_entries = 1e6)
    pc <- api$ipw_point_view_control(max_chunk_rows = 5L, max_observer_units = 1000L, max_retained_entries = 1e6,
      max_matrix_entries = 1e6, max_work_ops = 1e9, max_output_entries = 1e6)
    reg <- api$build_pipeline_nuisance_and_link_views(meta, o$fo, fp_ctrl(), fc, pc, list(epsilon = 1e-12, maxit = 50L, trace = FALSE),
      fit_point_links = FALSE)
    M <- nrow(meta$catalog$rows)
    w <- reg$construction_work
    if (case$reuse_all) {
      expect_identical(length(reg$nuisance_views), 1L)
      # One support scan (the view's own) and one streamed fingerprint serve
      # every row; all M fingerprint passes but the first are memo hits.
      expect_identical(w$link_support_scans_reused, as.double(M))
      expect_identical(w$fingerprint_passes, as.double(M))
      expect_identical(w$complete_pair_fingerprints_reused, as.double(M - 1L))
    } else {
      expect_gt(length(reg$nuisance_views), 1L)
      expect_lt(w$link_support_scans_reused, M)
    }
    # Standalone construction (no context) through the validated row helper.
    for (m in seq_len(M)) {
      nv <- reg$nuisance_views[[reg$link_to_nuisance_idx[m]]]
      lv <- api$construct_pipeline_catalog_row_views(meta, m, control = fp_ctrl(), frame_owner = o$fo, nuisance_view = nv)$link_view
      strip <- function(x) { if (is.environment(x)) return(NULL); if (is.list(x)) { a <- attributes(x); x <- lapply(x, strip); attributes(x) <- a }; x }
      expect_identical(strip(reg$link_views[[m]]), strip(lv))
      expect_identical(fp_ref_mask(o$d, lv$pair$s, lv$pair$t, c(3, Inf)),
        fp_ref_mask(o$d, nv$defining_pair$s, nv$defining_pair$t, c(3, Inf)))
    }
  }
})

test_that("validated chunk queries equal fully validated queries and refuse tampering as before", {
  o <- fp_owners(fp_panel())
  nv <- api$construct_ipw_nuisance_view(o$fo, 3, Inf, 1, 2, fp_ctrl())
  lv <- api$construct_ipw_link_view(o$fo, 2, 4, nuisance_view = nv, control = fp_ctrl())
  rec <- api$ipw_owner_validated_view_record(lv, o$fo, o$po)
  expect_true(rec$complete_pair)
  for (rng in list(c(1L, 5L), c(6L, 10L), c(21L, 21L))) {
    plain <- api$query_ipw_view_chunk(lv, rng[1], rng[2], fp_ctrl(), read_outcomes = TRUE, include_design = TRUE)
    expect_identical(api$query_ipw_view_chunk(lv, rng[1], rng[2], fp_ctrl(), read_outcomes = TRUE, include_design = TRUE,
      validated = rec), plain)
  }
  other <- api$construct_ipw_link_view(o$fo, 3, 5, nuisance_view = nv, control = fp_ctrl())
  expect_false(api$ipw_owner_validated_view_match(rec, other))
  expect_identical(api$query_ipw_view_chunk(other, 1L, 5L, fp_ctrl(), read_outcomes = TRUE, validated = rec),
    api$query_ipw_view_chunk(other, 1L, 5L, fp_ctrl(), read_outcomes = TRUE))
  bad <- lv; bad$pair$t <- 5
  expect_error(api$query_ipw_view_chunk(bad, 1L, 5L, fp_ctrl(), validated = rec), class = "cdid_ipw_owner_integrity_error")
  bad_rec <- rec; bad_rec$complete_pair <- FALSE
  expect_false(api$ipw_owner_validated_view_match(bad_rec, lv))
})

test_that("the IF projection consumer refuses a reappearing unit key across and within batches", {
  run <- function(batches) {
    n <- length(unlist(batches))
    spec <- api$gmm_if_replay_spec(callback = function(consume) for (b in seq_along(batches)) {
      B <- length(batches[[b]])
      consume(list(unit_id = batches[[b]], phi = matrix((seq_len(2L * B) + b) / 3, B, 2L), weight = rep(1, B)))
    }, input_row_level = "unit", aggregate_to_clusters = FALSE, n_input_rows = n, n_output_rows = n,
      n_positive_contributing_rows = n, centered = TRUE, weighted = FALSE, weight_type = "unweighted",
      normalization_factor = 1 / n, scale_id = "s", nuisance_fit_id = "nf")
    api$project_gmm_if_stream(spec, diag(2), 2L, api$gmm_bounds_control(max_matrix_entries = 1e6, max_batch_rows = 8L,
      max_targets = 4L, max_open_clusters = 50L, allow_materialize_if = TRUE, max_materialized_if_entries = 1e6),
      list(nuisance_fit_id = "nf"), materialize_target_if = TRUE)
  }
  ok <- run(list(1:4, 5:8, 9:10))
  expect_equal(ok$n_output_rows, 10)
  # Reference: the unit-level sandwich is the centered crossproduct of the rows.
  rows <- do.call(rbind, lapply(seq_along(list(1:4, 5:8, 9:10)), function(b) { B <- c(4, 4, 2)[b]; matrix((seq_len(2 * B) + b) / 3, B, 2) }))
  z <- sweep(rows, 2L, colMeans(rows))
  expect_equal(unname(ok$sandwich_cov), crossprod(z) / 10, tolerance = 1e-14)
  for (b in list(list(1:4, c(5, 2, 6)), list(c(1, 2, 1), 4:5), list(1:3, 4:6, c(7, 3)))) {
    expect_error(run(b), "Closed cluster reappeared.", fixed = TRUE, class = "cdid_gmm_replay_error")
  }
})

test_that("the payload budget equals plain checked accumulation, refusal points included", {
  ref_budget <- function(objects, ctl) {
    entries <- 0; bytes <- 0; work <- 128
    first <- NULL
    gate <- function(total, bound, label) if (is.null(first) && !is.null(bound) && total > bound) first <<- label
    visit <- function(x, depth = 0L, scan_attrs = TRUE) {
      if (is.null(x) || is.environment(x)) return(invisible(NULL))
      k <- length(x)
      entries <<- entries + k + 8; gate(entries, ctl$max_matrix_entries, "compact payload length gate")
      work <<- work + 64 * (k + 8); gate(work, ctl$max_work_ops, "compact payload traversal gate")
      if (scan_attrs && !is.null(attributes(x))) visit(attributes(x), depth + 1L, FALSE)
      if (is.list(x)) for (j in seq_len(k)) visit(x[[j]], depth + 1L) else if (is.character(x)) for (j in seq_len(k)) {
        nb <- nchar(x[[j]], type = "bytes", allowNA = TRUE); if (is.na(nb)) nb <- 0
        bytes <<- bytes + nb
        entries <<- entries + ceiling(nb / 2); gate(entries, ctl$max_matrix_entries, "compact character length gate")
        work <<- work + 16 * nb; gate(work, ctl$max_work_ops, "compact character comparison gate")
      }
    }
    visit(objects)
    list(first = first, entries = entries, work = work, bytes = bytes)
  }
  objs <- list(list(a = 1:3, b = list("x", c(u = 1.5)), m = matrix(1:4, 2)), list(letters, list(list(NA_character_, "été"))))
  for (obj in objs) for (ctl in list(list(max_matrix_entries = 1e6, max_work_ops = 1e9), list(max_matrix_entries = 20, max_work_ops = 1e9),
                                     list(max_matrix_entries = 1e6, max_work_ops = 1500))) {
    r <- ref_budget(obj, ctl)
    got <- tryCatch(api$ipw_point_view_payload_budget(obj, ctl), cdid_ipw_owner_bound_error = function(e) e)
    if (is.null(r$first)) {
      expect_identical(got$entries, r$entries); expect_identical(got$work, r$work); expect_identical(got$bytes, r$bytes)
    } else {
      expect_s3_class(got, "cdid_ipw_owner_bound_error")
      expect_match(conditionMessage(got), r$first, fixed = TRUE)
    }
  }
})

test_that("dependency checks are recorded only for the locked namespace", {
  api$ipw_owner_require_deps()
  expect_identical(api$cdid_dependency_check_done("ipw_owner", api), isNamespace(api) && environmentIsLocked(api))
  clone <- new.env(parent = api)
  f <- api$ipw_owner_require_deps; environment(f) <- clone; clone$ipw_owner_require_deps <- f
  clone$ipw_owner_require_deps()
  expect_false(api$cdid_dependency_check_done("ipw_owner", clone))
})


test_that("complete-history records are constructed only after the whole replay preflight", {
  o <- fp_owners(fp_panel())
  fc <- api$ipw_fitted_view_control(5, 1000, 1e6, 1e6, 1e9, 1e6)
  pc <- api$ipw_point_view_control(5, 1000, 1e6, 1e6, 1e9, 1e6)
  nv <- api$construct_ipw_nuisance_view(o$fo, 3, Inf, 1, 2, fp_ctrl())
  fv <- api$fit_ipw_nuisance_view(nv, stats::glm.control(epsilon = 1e-12, maxit = 50), fc)
  lv <- api$construct_ipw_link_view(o$fo, 2, 4, nuisance_view = nv, control = fp_ctrl())
  pl <- api$fit_ipw_point_link_view(fv, lv, pc)
  e <- new.env(parent = api)
  for (nm in c("normalize_ipw_point_links_view", "preflight_ipw_point_if_view", "ipw_point_view_add_query_records")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  calls <- 0L
  e$ipw_owner_validated_view_record <- function(...) { calls <<- calls + 1L; api$ipw_owner_validated_view_record(...) }
  pre <- e$preflight_ipw_point_if_view(list(pl), pc)
  expect_identical(calls, 0L)
  low <- pc; low$max_matrix_entries <- pre$peak_concurrent_entries - 1
  expect_error(e$preflight_ipw_point_if_view(list(pl), low, return_norm = TRUE), class = "cdid_ipw_owner_bound_error")
  expect_identical(calls, 0L)
  good <- e$preflight_ipw_point_if_view(list(pl), pc, return_norm = TRUE)
  expect_identical(calls, 1L)
  expect_true(good$norm$validated[[1L]]$complete_pair)
})

test_that("scalar fingerprint bucket collisions recompute and query preflight memo stays bounded", {
  o <- fp_owners(fp_panel()); ctrl <- fp_ctrl()
  nv <- api$construct_ipw_nuisance_view(o$fo, 3, Inf, 1, 2, ctrl)
  ctx <- api$ipw_owner_construction_context(o$fo)
  step <- min(o$fo$n_a_units, ctrl$max_chunk_rows)
  key <- api$ipw_owner_context_fingerprint_key(nv$cohorts$owner_cohort_index, step)
  # A colliding cached entry has the right key and step but a different complete
  # cohort vector. It must not supply the returned fingerprint.
  assign(key, list(owner_cohort_index = c(1L, 999L), step = step,
    result = list(fingerprint = c(count = -1, s1 = -1, s2 = -1, s3 = -1), chunks = 1)), ctx$fingerprints)
  expect_identical(api$fingerprint_ipw_view_mask(nv, ctrl, construction_context = ctx)$fingerprint,
    api$fingerprint_ipw_view_mask(nv, ctrl)$fingerprint)
  expect_identical(ctx$fingerprint_reuses, 0)
  expect_identical(api$ipw_owner_context_fingerprint_key(c(1, 2, 4), 5),
    api$ipw_owner_context_fingerprint_key(c(1, 3, 4), 5))
  rec <- api$ipw_owner_validated_view_record(nv, o$fo, o$po, step = 5)
  for (b in 1:5) for (design in c(FALSE, TRUE)) {
    cc <- api$ipw_owner_validated_control(rec, ctrl)
    expect_identical(api$ipw_owner_validated_preflight(rec, cc, b, 1L, 2L, TRUE, design, FALSE),
      api$preflight_ipw_range_query(b, 1L, 2L, ctrl, TRUE, design, FALSE))
    expect_named(rec$memo$preflight, c("args", "ctrl", "preflight"))
    expect_null(rec$memo$preflights)
  }
})

# ---------------------------------------------------------------------------
# Stage 5: exact invocation-local reuse of completed support scans by the
# nuisance-free link views of one CP04 registry. Expected views, support values
# and registry results come from the unchanged no-context constructor or the
# forced existing route of the same build, or from raw panel rows; the context
# never supplies its own truth. Scan, context and posting-date counts come from
# test-local instrumentation of re-bound copies.

rs_ctrl <- function(chunk = 7L, retained = 1e7, matrix = 1e7, work = 1e10) api$ipw_owner_control(
  max_chunk_rows = chunk, max_retained_entries = retained, max_matrix_entries = matrix, max_work_ops = work,
  max_output_entries = 1e7)

# Units without RNG: permuted numeric (or character) ids, cohorts cycling
# through `cohorts`, one covariate and nonconstant weights: dyadic eighths
# (exact sums, so every complete pair is supported on any platform) or
# non-dyadic tenths plus one third (rounded chunked sums).
rs_units <- function(N = 46L, char_ids = FALSE, cohorts = c(4, 6, Inf), dyadic = TRUE) {
  u <- seq_len(N)
  id <- (u * 37L) %% 1009L
  data.frame(id = if (char_ids) sprintf("id-%04d", id) else id, g = cohorts[u %% length(cohorts) + 1L],
    w = if (dyadic) 1 + (u %% 7L) / 8 else 1 + (u %% 7L) / 10 + 1 / 3, x = (u %% 5L) / 4, stringsAsFactors = FALSE)
}
# Long rows in a permuted order. `drop` removes and `na` blanks the outcome of
# "unit-position date" rows; `shift` changes every outcome value only.
rs_panel <- function(units, T = 6L, drop = NULL, na = NULL, shift = 0) {
  k <- rep(seq_len(nrow(units)), each = T)
  d <- data.frame(id = units$id[k], t = rep(seq_len(T), nrow(units)), g = units$g[k], stringsAsFactors = FALSE)
  d$y <- 0.3 * k + 0.2 * d$t + ifelse(d$t >= d$g, 1, 0) + ((7 * k + 3 * d$t) %% 5) / 10 + shift
  key <- paste(k, d$t)
  d$y[key %in% na] <- NA
  d <- d[!(key %in% drop), ]
  d[order((seq_len(nrow(d)) * 13L) %% nrow(d)), ]
}
rs_owners <- function(units, d, ctrl = rs_ctrl()) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  fo <- api$construct_ipw_frame_owner(po, units$id, cbind("(Intercept)" = 1, x = units$x), weights = units$w,
    control = ctrl)
  list(units = units, d = d, po = po, fo = fo)
}
rs_pairs <- function(T) {
  pr <- list()
  for (s in seq_len(T - 1L)) for (t in seq.int(s + 1L, T)) pr[[length(pr) + 1L]] <- c(s, t)
  pr
}
rs_specs <- function(classes, pairs) {
  out <- list()
  for (cl in classes) for (p in pairs) out[[length(out) + 1L]] <- list(s = p[1L], t = p[2L], g = cl$g, K = cl$K,
    functional_key = paste0("cmpcat:rs:g", cl$g, ":K", paste(cl$K, collapse = ","), ":", p[1L], "-", p[2L]))
  out
}

# The accepted stage-4 registry loop: no context, typed auxiliary capture.
rs_plain <- function(fo, specs, aux, ctrl) {
  views <- lapply(seq_along(specs), function(j) {
    s <- specs[[j]]
    api$construct_ipw_link_view(fo, s$s, s$t, g = s$g, K = s$K, control = ctrl, link_id = paste0("cp04:catalog:", j),
      functional_key = s$functional_key)
  })
  idx <- rep(NA_integer_, length(aux$descriptors)); errs <- vector("list", length(idx))
  for (j in seq_along(idx)) {
    s <- aux$descriptors[[j]]
    v <- tryCatch(api$construct_ipw_link_view(fo, s$s, s$t, g = s$g, K = s$K, control = ctrl,
      link_id = paste0("cp04:auxiliary:", j), functional_key = paste0("cmpcat:CP04-aux:g", s$g, ":K",
        paste(s$K, collapse = ","))), cdid_ipw_support_error = identity)
    if (inherits(v, "error")) errs[[j]] <- v else {
      views[[length(views) + 1L]] <- v; idx[j] <- length(views)
    }
  }
  list(views = views, aux_index = idx, aux_errors = errs)
}

# Re-bound copies counting actual support scans and posting-date reads and
# logging admission, context creation and scans in call order. `decline`
# forces the existing route (no context) through the same copies.
rs_bound <- c("build_cp04_public_registry", "cp04_public_link_views", "ipw_owner_link_support_context",
  "construct_ipw_link_view", "ipw_owner_context_link_lookup", "ipw_owner_context_record_link",
  "ipw_owner_pair_masks_equal_fact", "ipw_owner_pair_complete", "ipw_owner_date_complete")
rs_counted <- function(decline = FALSE) {
  e <- new.env(parent = api)
  for (nm in rs_bound) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  n <- new.env(parent = emptyenv())
  n$scans <- 0L; n$dates <- 0L; n$log <- character(); n$work <- NULL
  e$ipw_owner_chunked_pair_support <- function(...) {
    n$scans <- n$scans + 1L; n$log <- c(n$log, "scan")
    api$ipw_owner_chunked_pair_support(...)
  }
  e$ipw_owner_posting_date_complete <- function(...) { n$dates <- n$dates + 1L; api$ipw_owner_posting_date_complete(...) }
  e$ipw_owner_construction_context <- function(...) { n$log <- c(n$log, "context"); api$ipw_owner_construction_context(...) }
  admit <- e$ipw_owner_link_support_context
  e$ipw_owner_link_support_context <- function(...) {
    n$log <- c(n$log, "admission")
    ctx <- if (decline) NULL else admit(...)
    n$log <- c(n$log, if (is.null(ctx)) "refused" else "admitted")
    ctx
  }
  views_fn <- e$cp04_public_link_views
  e$cp04_public_link_views <- function(...) { out <- views_fn(...); n$work <- out$construction_work; out }
  e$rs_tally <- n
  e
}
rs_run <- function(o, specs, aux = list(descriptors = list()), ctrl = rs_ctrl(), retained = 0) {
  e <- rs_counted()
  out <- e$cp04_public_link_views(o$fo, specs, aux, ctrl, retained)
  c(out, list(counts = as.list(e$rs_tally)))
}

rs_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, rs_strip); attributes(x) <- a }
  x
}
rs_expect_same_views <- function(got, want) {
  expect_identical(length(got), length(want))
  for (j in seq_along(want)) {
    expect_identical(rs_strip(got[[j]]), rs_strip(want[[j]]))
    expect_true(identical(rs_strip(got[[j]]), rs_strip(want[[j]]), num.eq = FALSE))
  }
}
rs_work <- function(x) x$construction_work[c("support_scans", "support_reuses", "support_records")]
rs_counts <- function(scans, reuses, records) list(support_scans = scans, support_reuses = reuses,
  support_records = records)

# Raw-row reference of a link's support summary: A-units in the frame's A order
# observed (non-missing outcome rows) at both dates; per declared cohort the
# count and the weight mass accumulated chunk by chunk in A order, and the rate
# over the cohort's A mass (A order), with complete cohort support (cnt == a_cnt)
# yielding exact probability 1 when require_support is TRUE.
rs_check_support <- function(v, o, step, require_support = TRUE) {
  a_id <- o$fo$a_id
  obs <- o$d[!is.na(o$d$y), ]
  hit <- a_id %in% obs$id[obs$t == v$pair$s] & a_id %in% obs$id[obs$t == v$pair$t]
  g <- o$d$g[match(a_id, o$d$id)]
  w <- o$units$w[match(a_id, o$units$id)]
  cohorts <- c(v$cohorts$g, v$cohorts$K)
  cnt <- integer(length(cohorts)); mass <- numeric(length(cohorts)); n_a <- length(a_id)
  for (a0 in seq.int(1L, n_a, by = step)) {
    idx <- seq.int(a0, min(n_a, a0 + step - 1L))
    idx <- idx[hit[idx]]
    for (j in seq_along(cohorts)) {
      h <- idx[g[idx] == cohorts[j]]
      if (length(h)) { cnt[j] <- cnt[j] + length(h); mass[j] <- mass[j] + sum(w[h]) }
    }
  }
  a_cnt <- vapply(cohorts, function(k) sum(g == k), integer(1L))
  a_mass <- vapply(cohorts, function(k) sum(w[g == k]), numeric(1L))
  rates <- mass / a_mass
  if (isTRUE(require_support)) rates[cnt == a_cnt] <- 1
  expect_identical(unname(v$cohorts$group_a_count), a_cnt)
  expect_identical(unname(v$cohorts$group_a_mass), a_mass)
  expect_identical(unname(v$support$by_cohort_pair_count), cnt)
  expect_identical(unname(v$support$by_cohort_pair_mass), mass)
  expect_identical(unname(v$support$rates), rates)
  expect_identical(v$support$n_observed_pair_all, sum(hit))
}

test_that("registry views borrow one scan per exact class and equal the no-context constructor", {
  classes <- list(list(g = 4, K = Inf), list(g = 4, K = c(6, Inf)), list(g = 4, K = c(Inf, 6)),
    list(g = 6, K = Inf), list(g = 6, K = c(Inf, 4)))
  aux <- list(descriptors = list(list(g = 6, s = 2, t = 5, K = c(Inf, 4)), list(g = 4, s = 1, t = 2, K = c(99, Inf)),
    list(g = 6, s = 1, t = 3, K = c(4, Inf))))
  specs <- rs_specs(classes, rs_pairs(6L))
  for (char_ids in c(FALSE, TRUE)) for (chunk in c(7L, 16L)) {
    units <- rs_units(46L, char_ids)
    o <- rs_owners(units, rs_panel(units))
    ctrl <- rs_ctrl(chunk)
    got <- rs_run(o, specs, aux, ctrl)
    want <- rs_plain(o$fo, specs, aux, ctrl)
    rs_expect_same_views(got$views, want$views)
    expect_identical(got$aux_index, want$aux_index)
    expect_identical(got$aux_errors, want$aux_errors)
    expect_s3_class(got$aux_errors[[2L]], "cdid_ipw_owner_cohort_support_error")
    # Six exact classes (five catalog classes and the auxiliary (6; 4, Inf)):
    # one scan each, no scan for any later equal-mask view, each date read once.
    expect_identical(got$counts$scans, 6L)
    expect_identical(rs_work(got), rs_counts(6, length(got$views) - 6, 6))
    expect_identical(got$counts$dates, 6L)
    for (v in got$views) rs_check_support(v, o, min(46L, chunk))
    expect_true(all(vapply(got$views, function(v) identical(v$accounting$endpoints_read, 0L), logical(1L))))
  }
  # Outcome values never enter: shifted outcomes give the same views.
  units <- rs_units(46L)
  sub <- rs_specs(classes[2:3], rs_pairs(6L))
  rs_expect_same_views(rs_run(rs_owners(units, rs_panel(units, shift = 2.5)), sub)$views,
    rs_run(rs_owners(units, rs_panel(units)), sub)$views)
})

test_that("non-dyadic weights: every view and typed refusal equals the original chunked scan bitwise", {
  # Non-dyadic weights round in the chunked pair-mass sums, which the borrowed
  # values must reproduce bitwise in by_cohort_pair_mass. Under require_support = TRUE,
  # equal positive cohort counts (pair_count == group_a_count) prove complete
  # cohort coverage on the validated unique-(ID, time) roster and have exact
  # mathematical probability 1, whereas require_support = FALSE retains the
  # unchecked mass / a_mass quotient. Expected outcomes come from the no-context
  # constructor of the same build, view by view.
  classes <- list(list(g = 4, K = Inf), list(g = 4, K = c(6, Inf)), list(g = 4, K = c(Inf, 6)),
    list(g = 6, K = c(Inf, 4)))
  desc <- lapply(rs_specs(classes, rs_pairs(6L)), function(s) list(g = s$g, s = s$s, t = s$t, K = s$K))
  cls <- vapply(desc, function(x) paste(x$g, paste(x$K, collapse = ",")), "")
  for (char_ids in c(FALSE, TRUE)) for (chunk in c(7L, 16L)) for (rule in c(TRUE, FALSE)) {
    units <- rs_units(46L, char_ids, dyadic = FALSE)
    o <- rs_owners(units, rs_panel(units))
    ctrl <- rs_ctrl(chunk)
    e <- rs_counted()
    ctx <- e$ipw_owner_link_support_context(o$fo, length(desc), 0, ctrl)
    ok <- logical(length(desc))
    for (j in seq_along(desc)) {
      x <- desc[[j]]
      v <- tryCatch(e$construct_ipw_link_view(o$fo, x$s, x$t, g = x$g, K = x$K, control = ctrl, require_support = rule,
        construction_context = ctx), error = identity)
      w <- tryCatch(api$construct_ipw_link_view(o$fo, x$s, x$t, g = x$g, K = x$K, control = ctrl,
        require_support = rule), error = identity)
      ok[j] <- !inherits(w, "error")
      if (ok[j]) {
        expect_true(identical(rs_strip(v), rs_strip(w), num.eq = FALSE))
        rs_check_support(v, o, chunk, require_support = rule)
      } else {
        expect_s3_class(w, "cdid_ipw_rate_error")
        expect_identical(v, w)
      }
    }
    # One mask per class: a class borrows its one completed scan, or each of
    # its views scans again and is refused (a failed scan is never recorded).
    expect_identical(e$rs_tally$scans, sum(vapply(split(ok, cls), function(z) if (all(z)) 1L else length(z), 1L)))
    if (!rule) expect_true(all(ok))
    if (rule) {
      # The registry helper keeps the outcome: the same views or the same refusal.
      r1 <- tryCatch(rs_run(o, list(), list(descriptors = desc), ctrl)[c("views", "aux_index", "aux_errors")],
        error = identity)
      r0 <- tryCatch(rs_plain(o$fo, list(), list(descriptors = desc), ctrl), error = identity)
      expect_identical(rs_strip(r1), rs_strip(r0))
    }
  }
})

test_that("missing, NA and repeated pairs, step, rule and foreign or forged records reuse only proved-equal masks", {
  units <- rs_units(46L)
  # Unit 5 lacks date 2 and unit 11 has an NA outcome at date 5: dates 2 and 5
  # are incomplete, dates 1, 3, 4 and 6 complete.
  o <- rs_owners(units, rs_panel(units, drop = "5 2", na = "11 5"))
  ctrl <- rs_ctrl(7L)
  pairs <- rs_pairs(6L)
  c1 <- list(g = 4, K = c(6, Inf)); c2 <- list(g = 4, K = c(Inf, 6))
  specs <- c(rs_specs(list(c1), pairs), rs_specs(list(c1), list(c(2, 3), c(1, 3), c(5, 6))), rs_specs(list(c2), pairs))
  got <- rs_run(o, specs, ctrl = ctrl)
  want <- rs_plain(o$fo, specs, list(descriptors = list()), ctrl)
  rs_expect_same_views(got$views, want$views)
  for (v in got$views) rs_check_support(v, o, 7L)
  # Per ordered-K class: one scan serves every complete pair, each distinct
  # incomplete pair scans once, repeated pairs borrow; incomplete dates are
  # rejected by their counts, never read.
  incomplete <- sum(!vapply(pairs, function(p) all(p %in% c(1, 3, 4, 6)), logical(1L)))
  expect_identical(got$counts$scans, 2L * (1L + incomplete))
  expect_identical(rs_work(got), rs_counts(2 * (1 + incomplete), length(specs) - 2 * (1 + incomplete),
    2 * (1 + incomplete)))
  expect_identical(got$counts$dates, 4L)

  e <- rs_counted(); n <- e$rs_tally
  ctx <- e$ipw_owner_link_support_context(o$fo, 50, 0, ctrl)
  mk <- function(fo, s, t, K = c(6, Inf), ct = ctrl, rule = TRUE) e$construct_ipw_link_view(fo, s, t, g = 4, K = K,
    control = ct, require_support = rule, construction_context = ctx)
  plain <- function(fo, s, t, K = c(6, Inf), ct = ctrl, rule = TRUE) api$construct_ipw_link_view(fo, s, t, g = 4,
    K = K, control = ct, require_support = rule)
  v13 <- mk(o$fo, 1, 3)
  expect_identical(n$scans, 1L)
  expect_identical(rs_strip(mk(o$fo, 3, 4)), rs_strip(plain(o$fo, 3, 4)))
  expect_identical(n$scans, 1L)
  # A different support rule or chunk step is a different class.
  expect_identical(rs_strip(mk(o$fo, 3, 4, rule = FALSE)), rs_strip(plain(o$fo, 3, 4, rule = FALSE)))
  expect_identical(rs_strip(mk(o$fo, 3, 4, ct = rs_ctrl(16L))), rs_strip(plain(o$fo, 3, 4, ct = rs_ctrl(16L))))
  expect_identical(n$scans, 3L)
  # Another sealed frame of the same rows is foreign to the context; a modified
  # frame owner is refused before any scan, as without a context.
  o2 <- rs_owners(units, rs_panel(units, drop = "5 2", na = "11 5"))
  records <- ctx$link_records
  expect_identical(rs_strip(mk(o2$fo, 3, 4)), rs_strip(plain(o2$fo, 3, 4)))
  expect_identical(c(n$scans, ctx$link_records), c(4L, records))
  bad <- o$fo; bad$w_a[1L] <- bad$w_a[1L] + 1
  r1 <- tryCatch(mk(bad, 3, 4), error = identity); r0 <- tryCatch(plain(bad, 3, 4), error = identity)
  expect_s3_class(r1, "cdid_ipw_owner_integrity_error")
  expect_identical(r1, r0)
  expect_identical(n$scans, 4L)
  # Bounded keys collide for different descriptors; a record whose sealed view
  # does not match it exactly (scan counts, unsealed reference, other ordered
  # K) never grants reuse and is never replaced by a later scan.
  expect_identical(api$ipw_owner_context_link_key(list(a_cohort_map = c(1L, 2L, 5L, 4L)), 7L, TRUE),
    api$ipw_owner_context_link_key(list(a_cohort_map = c(1L, 2L, 3L, 4L)), 7L, TRUE))
  key <- api$ipw_owner_context_link_key(v13$cohorts, v13$preflight$effective_chunk_rows, TRUE)
  good <- get(key, envir = ctx$links)
  other <- mk(o$fo, 1, 3, K = c(Inf, 6))
  forged <- list(counts = utils::modifyList(good, list(ops = good$ops + 1)),
    unsealed = utils::modifyList(good, list(ref = list2env(list(link_view = v13)))),
    colliding = utils::modifyList(good, list(ref = other$ref)))
  for (f in forged) {
    assign(key, f, envir = ctx$links)
    before <- n$scans
    expect_identical(rs_strip(mk(o$fo, 4, 6)), rs_strip(plain(o$fo, 4, 6)))
    expect_identical(n$scans, before + 1L)
    expect_identical(get(key, envir = ctx$links), f)
  }
  assign(key, good, envir = ctx$links)
  before <- n$scans
  expect_identical(rs_strip(mk(o$fo, 4, 6)), rs_strip(plain(o$fo, 4, 6)))
  expect_identical(n$scans, before)
})

test_that("the context is admitted before its creation and composed with each view's constructor envelope", {
  units <- rs_units(46L)
  o <- rs_owners(units, rs_panel(units))
  specs <- rs_specs(list(list(g = 4, K = Inf)), rs_pairs(6L))
  n_views <- length(specs); entries <- 128 * (n_views + 1)
  ctrl <- rs_ctrl(7L)
  want <- rs_plain(o$fo, specs, list(descriptors = list()), ctrl)
  ok <- rs_run(o, specs, ctrl = ctrl)
  expect_identical(ok$counts$log, c("admission", "context", "admitted", "scan"))
  expect_identical(ok$construction_work$reserved_entries, entries)
  edge <- rs_run(o, specs, ctrl = ctrl, retained = 1e7 - entries)
  expect_identical(edge$counts$log, c("admission", "context", "admitted", "scan"))
  # Refused against the registry envelope: no context; every view scans.
  for (retained in c(1e7 - entries + 1, 1e7)) {
    no <- rs_run(o, specs, ctrl = ctrl, retained = retained)
    expect_identical(no$counts$log, c("admission", "refused", rep("scan", n_views)))
    expect_null(no$construction_work)
    rs_expect_same_views(no$views, want$views)
  }
  # Not composable with the largest constructor: decline before allocating
  # any context, then preserve the no-context scans and original views.
  pre <- api$preflight_ipw_link_view(o$fo$n_a_units, o$fo$n_covariates, 2L, ctrl)
  cap <- entries + pre$concurrent_entries - 1
  tight <- rs_ctrl(7L, retained = cap, matrix = cap)
  un <- rs_run(o, specs, ctrl = tight)
  expect_identical(un$counts$log, c("admission", "refused", rep("scan", n_views)))
  expect_null(un$construction_work)
  rs_expect_same_views(un$views, rs_plain(o$fo, specs, list(descriptors = list()), tight)$views)
  # The constructor's own bound refusal is unchanged and raised before any
  # scan; nothing is retried without the context.
  low <- rs_ctrl(7L, work = 10)
  e1 <- rs_counted()
  r1 <- tryCatch(e1$cp04_public_link_views(o$fo, specs, list(descriptors = list()), low, 0), error = identity)
  r0 <- tryCatch(rs_plain(o$fo, specs, list(descriptors = list()), low), error = identity)
  expect_s3_class(r1, "cdid_ipw_owner_bound_error")
  expect_identical(r1, r0)
  expect_identical(e1$rs_tally$log, c("admission", "context", "admitted"))
})

test_that("typed support failures keep their place and order and are never recorded", {
  units <- rs_units(45L)
  # Every cohort-6 unit lacks date 2.
  o <- rs_owners(units, rs_panel(units, drop = paste(which(units$g == 6), 2)))
  ctrl <- rs_ctrl(7L)
  aux <- list(descriptors = list(list(g = 4, s = 1, t = 2, K = c(6, Inf)), list(g = 4, s = 1, t = 3, K = c(6, Inf)),
    list(g = 4, s = 2, t = 4, K = c(6, Inf)), list(g = 4, s = 3, t = 4, K = c(99, Inf)),
    list(g = 4, s = 3, t = 4, K = c(6, Inf))))
  got <- rs_run(o, list(), aux, ctrl)
  want <- rs_plain(o$fo, list(), aux, ctrl)
  rs_expect_same_views(got$views, want$views)
  expect_identical(got$aux_index, c(NA, 1L, NA, NA, 2L))
  expect_identical(got$aux_index, want$aux_index)
  expect_identical(got$aux_errors, want$aux_errors)
  for (j in c(1L, 3L)) expect_s3_class(got$aux_errors[[j]], "cdid_ipw_owner_pair_support_error")
  expect_s3_class(got$aux_errors[[4L]], "cdid_ipw_owner_cohort_support_error")
  # Two failed scans and one successful scan ran; only the last view borrowed.
  expect_identical(got$counts$scans, 3L)
  expect_identical(rs_work(got), rs_counts(3, 1, 1))
  # A failing catalog view still ends the registry with its own typed refusal.
  specs <- rs_specs(list(list(g = 4, K = c(6, Inf))), list(c(1, 3), c(1, 2)))
  r1 <- tryCatch(rs_run(o, specs, ctrl = ctrl), error = identity)
  r0 <- tryCatch(rs_plain(o$fo, specs, list(descriptors = list()), ctrl), error = identity)
  expect_s3_class(r1, "cdid_ipw_owner_pair_support_error")
  expect_identical(r1, r0)
})

test_that("a 160-unit, 20-date registry borrows one scan for all 190 pairs", {
  units <- rs_units(160L, cohorts = c(10, Inf, Inf, Inf))
  o <- rs_owners(units, rs_panel(units, T = 20L))
  ctrl <- rs_ctrl(64L)
  specs <- rs_specs(list(list(g = 10, K = Inf)), rs_pairs(20L))
  expect_identical(length(specs), 190L)
  got <- rs_run(o, specs, ctrl = ctrl)
  want <- rs_plain(o$fo, specs, list(descriptors = list()), ctrl)
  rs_expect_same_views(got$views, want$views)
  expect_identical(c(got$counts$scans, got$counts$dates), c(1L, 20L))
  expect_identical(rs_work(got), rs_counts(1, 189, 1))
  for (v in got$views[c(1L, 95L, 190L)]) rs_check_support(v, o, 64L)
})

# Registry results that do not depend on fresh environments: points (link,
# estimate, cohorts, pair), auxiliary map and errors, nuisance keys and the
# complete streamed covariance of every link.
rs_registry_view <- function(reg) {
  links <- seq_along(reg$cp04_points)
  list(points = rs_strip(reg$cp04_points), aux_index = reg$cp04_aux_index, aux_errors = reg$cp04_aux_errors,
    keys = reg$nuisance_keys,
    stream = api$cp04_stream_covariance(reg$cp04, links)[c("estimates", "covariance", "mean_phi")])
}
rs_compare_registries <- function(build) {
  got_e <- rs_counted(); want_e <- rs_counted(decline = TRUE)
  got <- build(got_e); want <- build(want_e)
  expect_true(identical(rs_registry_view(got), rs_registry_view(want), num.eq = FALSE))
  expect_identical(rs_registry_view(got), rs_registry_view(want))
  # Every view reaching its support step either scanned or borrowed.
  expect_identical(want_e$rs_tally$work, NULL)
  expect_identical(got_e$rs_tally$scans, as.integer(got_e$rs_tally$work$support_scans))
  expect_identical(got_e$rs_tally$work$support_scans + got_e$rs_tally$work$support_reuses,
    as.double(want_e$rs_tally$scans))
  list(got = got, want = want, scans = c(got_e$rs_tally$scans, want_e$rs_tally$scans))
}

test_that("CP04 OR, DR and conditional registries with the context equal the forced existing route", {
  units <- rs_units(36L, cohorts = c(4, 5, Inf))
  ctrl <- rs_ctrl(7L)
  po <- api$construct_ipw_panel_owner(rs_panel(units), "id", "t", "y", "g", ctrl)
  fo <- api$construct_ipw_frame_owner(po, units$id, matrix(1, nrow(units), 1L, dimnames = list(NULL, "(Intercept)")),
    weights = units$w, control = ctrl, frame_g = units$g, allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, units$id, matrix(1, nrow(units), 1L,
    dimnames = list(NULL, "reader_intercept")), units$g, units$w, ctrl, allow_unknown_off_A = TRUE)
  rctrl <- list(owner_control = ctrl, glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE))
  for (method in c("or", "dr")) {
    cmp <- rs_compare_registries(function(e) e$build_cp04_public_registry(NULL, fo, reader, data.frame(g = c(4, 5)),
      rctrl, "notyettreated", 0L, or_method = if (method == "or") "baseline",
      dr_method = if (method == "dr") "baseline"))
    # Complete balanced dates: one scan per exact (g, ordered K) class.
    desc <- cmp$got$cp04_aux$descriptors
    expect_identical(cmp$scans, c(length(unique(lapply(desc, function(x) c(x$g, x$K)))), length(desc)))
    expect_lt(cmp$scans[1L], cmp$scans[2L])
  }
  # The accepted unbalanced conditional fixture (X, non-dyadic weights,
  # frame-only units, U-conditional selection): identical results.
  raw <- utils::read.csv(testthat::test_path("fixtures", "cp06_or_corollary1", "raw.csv"))
  raw$cohort <- as.numeric(replace(raw$cohort, raw$cohort == "NT", "Inf"))
  N <- nrow(raw); Y <- as.matrix(raw[paste0("Y", 1:4)])
  long <- data.frame(id = rep(raw$id, each = 4L), time = rep(1:4, N), y = as.vector(t(Y)), g = rep(raw$cohort, each = 4L))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  long <- long[!long$id %in% raw$id[zero[-1L]], ]
  ctrl6 <- api$ipw_owner_control(23L, 134217728, 134217728, 1e11, 134217728)
  panel <- api$construct_ipw_panel_owner(long, "id", "time", "y", "g", ctrl6)
  X <- cbind(intercept = 1, x = raw$x); X[raw$A == 0, ] <- NA_real_
  frame <- api$construct_ipw_frame_owner(panel, raw$id, X, weights = raw$w, control = ctrl6, frame_g = raw$cohort,
    allow_frame_only = TRUE, design_on_A_only = TRUE)
  reader6 <- api$construct_cp04_frame_reader(frame, raw$id, cbind(one = 1, x = raw$x), raw$cohort, raw$w, ctrl6,
    U = rep(TRUE, N), rho = rep(1, N), universe = "F", justification = "accepted deterministic finite fixture")
  cmp <- rs_compare_registries(function(e) e$build_cp04_public_registry(NULL, frame, reader6, data.frame(g = c(3, 4)),
    list(owner_control = ctrl6, glm_control = list(epsilon = 1e-13, maxit = 60L, trace = FALSE)),
    "notyettreated", 0L, or_method = "corollary1", automatic_chunk_rows = TRUE))
  expect_lte(cmp$scans[1L], cmp$scans[2L])
})


test_that("all declared classes reserve their constructor before the registry holds a context", {
  units <- rs_units(96L, cohorts = c(2:32, Inf))
  o <- rs_owners(units, rs_panel(units, T = 3L))
  specs <- rs_specs(list(list(g = 2, K = Inf), list(g = 2, K = c(3:32, Inf))), list(c(1, 2)))
  ctrl <- rs_ctrl(7L, retained = 2000, matrix = 2000)
  small <- api$preflight_ipw_link_view(96, 2, 2, ctrl)
  large <- api$preflight_ipw_link_view(96, 2, 32, ctrl)
  reservation <- 128 * (length(specs) + 1)
  expect_lte(small$concurrent_entries + reservation, ctrl$max_matrix_entries)
  expect_lte(large$concurrent_entries, ctrl$max_matrix_entries)
  expect_gt(large$concurrent_entries + reservation, ctrl$max_matrix_entries)
  got <- rs_run(o, specs, ctrl = ctrl)
  expect_identical(got$counts$log, c("admission", "refused", "scan", "scan"))
  expect_null(got$construction_work)
  rs_expect_same_views(got$views, rs_plain(o$fo, specs, list(descriptors = list()), ctrl)$views)
  # Moving the large class to the auxiliary sequence has the same live peak.
  aux <- list(descriptors = specs[2L])
  got <- rs_run(o, specs[1L], aux, ctrl)
  expect_identical(got$counts$log, c("admission", "refused", "scan", "scan"))
  rs_expect_same_views(got$views, rs_plain(o$fo, specs[1L], aux, ctrl)$views)
  # Independent comparison of the admitted arithmetic envelope with each
  # original constructor preflight, at the exact and one-below boundaries.
  for (m in c(2L, 3L, 32L)) for (chunk in c(7L, 16L)) {
    hi <- rs_ctrl(chunk)
    pre <- api$preflight_ipw_link_view(96, 2, m, hi)
    for (offset in c(-1, 0)) {
      bounded <- rs_ctrl(chunk, retained = 1e7, matrix = reservation + pre$concurrent_entries + offset)
      ctx <- api$ipw_owner_link_support_context(o$fo, 2, 0, bounded, max_cohorts = m)
      expect_identical(is.null(ctx), offset < 0)
    }
  }
})
