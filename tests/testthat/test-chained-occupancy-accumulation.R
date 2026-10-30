# CP08 P2: exact bounded occupancy accumulation. Strict actual-installed
# namespace tests; the container source audit uses a separate disclosed adapter.
# >>> installed-namespace binding
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api), environmentIsLocked(api))
api_binding <- list(mode = "installed_package", env = api)
# <<< installed-namespace binding
ref_env <- new.env(parent = baseenv())
sys.source(testthat::test_path("fixtures", "pipeline_metadata_reference.R"), envir = ref_env)

occ_ctrl <- function(max_chunk_rows = 6L, max_retained_entries = 1e7, max_matrix_entries = 1e7,
                     max_work_ops = 1e10, max_output_entries = 1e7) {
  api$ipw_owner_control(max_chunk_rows = max_chunk_rows, max_retained_entries = max_retained_entries,
    max_matrix_entries = max_matrix_entries, max_work_ops = max_work_ops, max_output_entries = max_output_entries)
}
occ_cat_bounds <- function() list(max_candidates = 1e5, max_rows = 1e5, max_pool_entries = 1e6, max_scratch = 1e7, max_work = 1e10)

# Independent reference: direct loops over units and date pairs (no cdid code).
occ_reference_counts <- function(df) {
  obs <- df[!is.na(df$y), ]
  cal <- sort(unique(as.double(df$time)))
  coh <- sort(unique(as.double(df$g)))
  out <- list()
  for (u in unique(obs$id)) {
    tu <- sort(as.double(obs$time[obs$id == u]))
    if (length(tu) < 2L) next
    gu <- as.double(obs$g[obs$id == u][1L])
    for (i in seq_len(length(tu) - 1L)) for (j in seq.int(i + 1L, length(tu))) {
      key <- paste(match(gu, coh), match(tu[i], cal), match(tu[j], cal))
      out[[key]] <- c(match(gu, coh), match(tu[i], cal), match(tu[j], cal), if (is.null(out[[key]])) 1 else out[[key]][4L] + 1)
    }
  }
  if (!length(out)) {
    return(data.frame(cohort = numeric(0L), s = numeric(0L), t = numeric(0L), count = numeric(0L)))
  }
  m <- do.call(rbind, out)
  m <- m[order(m[, 1L], m[, 2L], m[, 3L]), , drop = FALSE]
  data.frame(cohort = coh[m[, 1L]], s = cal[m[, 2L]], t = cal[m[, 3L]], count = as.double(m[, 4L]))
}

# Every observation pattern over four dates, for each of three cohorts, with
# multiplicities 1..3; empty patterns are zero-observer units (explicit NA rows
# or one NA row), single-date patterns are singletons; missing dates alternate
# between explicit NA rows and absent rows.
occ_exhaustive_panel <- function(dates, cohorts) {
  rows <- list(); uid <- 0L
  for (ci in seq_along(cohorts)) for (p in 0:15) {
    obs_set <- dates[bitwAnd(p, 2L^(0:3)) > 0L]
    mult <- 1L + ((p + ci) %% 3L)
    for (r in seq_len(mult)) {
      uid <- uid + 1L
      explicit <- (uid %% 2L) == 0L
      tt <- if (explicit || !length(obs_set)) dates else obs_set
      y <- ifelse(tt %in% obs_set, 0.5 * uid + seq_along(tt), NA_real_)
      rows[[length(rows) + 1L]] <- data.frame(id = uid, time = tt, y = y, g = cohorts[ci])
    }
  }
  do.call(rbind, rows)
}

test_that("exhaustive four-date pattern counts match two independent references for irregular and integer labels", {
  for (labels in list(c(-3.5, 0, 2.25, 10), c(2001L, 2003L, 2004L, 2010L))) {
    df <- occ_exhaustive_panel(labels, cohorts = c(2.25, 10, Inf))
    if (is.integer(labels)) df$g <- c(2003, 2010, Inf)[match(df$g, c(2.25, 10, Inf))]
    ctrl <- occ_ctrl()
    po <- api$construct_ipw_panel_owner(df, "id", "time", "y", "g", ctrl)
    tr <- new.env(parent = emptyenv())
    occ <- api$derive_owner_sparse_occupancy(po, ctrl, trace_env = tr)
    ref <- occ_reference_counts(df)
    bf <- ref_env$brute_force_panel_occupancy(df, "id", "time", "y", "g")
    expect_identical(occ$occupancy, bf$sparse_occupancy)
    expect_identical(occ$cell_indices, bf$cell_indices)
    expect_equal(unname(as.list(occ$occupancy)), unname(as.list(ref)), tolerance = 0)
    expect_identical(occ$diagnostics$total_occupancy_count, as.integer(po$counts$n_pairs))
    expect_identical(tr$calls, 0L)
    expect_identical(occ$diagnostics$endpoints_read, 0L)
    # Zero-observer and singleton units stay in the owner and contribute no cell.
    expect_identical(occ$owner_summary$n_unobserved_units, bf$n_unobserved_units)
    expect_identical(occ$owner_summary$n_singleton_units, bf$n_singleton_units)
    expect_true(bf$n_unobserved_units > 0L && bf$n_singleton_units > 0L)
    # Repeated histories are accumulated once each with exact multiplicity.
    expect_lt(occ$diagnostics$n_distinct_histories, po$n_a_units)
    expect_identical(occ$diagnostics$n_pairs_represented, po$counts$n_pairs)
    expect_identical(occ$diagnostics$pair_visits_avoided, po$counts$n_pairs - occ$diagnostics$n_history_pair_visits)
    expect_gt(occ$diagnostics$pair_visits_avoided, 0)
    # Never-treated and not-yet-treated catalogs equal those from the reference occupancy (sparse and dense).
    tg <- sort(unique(df$g[is.finite(df$g)]))[1L]
    for (cg in c("nevertreated", "notyettreated")) {
      meta <- api$build_owner_comparison_catalog(po, targets = tg, anticipation = 0, control_group = cg,
        conventions = paste0("conv:exhaustive:", cg), control = ctrl, catalog_bounds = occ_cat_bounds())
      for (occ_ref in list(bf$sparse_occupancy, bf$dense_occupancy)) {
        direct <- api$build_comparison_catalog(calendar = po$all_times, cohorts = po$cohorts, targets = tg,
          anticipation = 0, occupancy = occ_ref, control_group = cg, conventions = paste0("conv:exhaustive:", cg),
          requested_post = NULL, bounds = occ_cat_bounds())
        expect_identical(meta$catalog$rows, direct$rows)
        expect_identical(meta$catalog$operators, direct$operators)
      }
      bf_map <- ref_env$brute_force_catalog_mappings(bf, meta$catalog)
      expect_identical(meta$row_mapping$target_pair_count, bf_map$target_pair_count)
      expect_identical(meta$row_mapping$control_pool_pair_count, bf_map$control_pool_pair_count)
      expect_identical(meta$row_mapping$all_cohort_pair_count, bf_map$all_cohort_pair_count)
    }
  }
})

test_that("minimal slice budgets, frame-only units and weights leave exact counts unchanged", {
  set.seed(20261022L)
  n <- 60L; dates <- c(1, 2, 3.5, 4, 7, 8, 9.25, 12)
  rows <- lapply(seq_len(n), function(i) {
    keep <- dates[stats::runif(length(dates)) < 0.6]
    if (i %% 7L == 0L) keep <- dates
    if (!length(keep)) keep <- dates[1L]
    data.frame(id = i, time = keep, y = stats::rnorm(length(keep)), g = c(3.5, 8, Inf)[1L + (i %% 3L)],
      x = (i %% 5L) / 5)
  })
  df <- do.call(rbind, rows)
  bf <- ref_env$brute_force_panel_occupancy(df, "id", "time", "y", "g")
  base_ctrl <- occ_ctrl()
  po <- api$construct_ipw_panel_owner(df, "id", "time", "y", "g", base_ctrl)
  u <- df[!duplicated(df$id), ]
  fo <- api$construct_ipw_frame_owner(po, id = c(u$id, 900L, 901L), Z = cbind(1, c(u$x, 0, 1)),
    weights = c(seq_len(nrow(u)) / 7, 2, 3), control = base_ctrl, frame_g = c(u$g, Inf, 8),
    allow_frame_only = TRUE)
  pre <- api$preflight_owner_sparse_occupancy(po, base_ctrl)
  t_span <- po$n_calendar_dates - 1
  min_me <- pre$retained_vector_entries + pre$header_entries +
    max(pre$grouping_scratch_entries, pre$accumulation_scratch_entries - 14 * pre$slice_pairs + 14 * t_span,
        pre$assembly_scratch_entries)
  for (me in c(1e7, min_me + 14 * 3, min_me)) {
    occ <- api$derive_owner_sparse_occupancy(fo, occ_ctrl(max_matrix_entries = me))
    expect_identical(occ$occupancy, bf$sparse_occupancy, info = me)
    expect_identical(occ$cell_indices, bf$cell_indices, info = me)
    expect_identical(occ$owner_summary$n_frame_only_units, 2L)
    expect_lte(occ$diagnostics$max_slice_pairs, occ$preflight$slice_pairs)
    if (me == min_me) expect_gt(occ$diagnostics$n_slices, occ$diagnostics$n_rows)
  }
})

test_that("deterministic operation counts: a balanced T=100 panel is admitted under the D27 caps that refused the merge charge", {
  # Public D27 defaults: 2^30 bytes of numeric payload and 1e11 counted work. The
  # unchanged public chunk rule gives 6 units per chunk for this request shape
  # (links_bound = 1 target x 4950 pairs x 2 cohorts).
  entries <- 1073741824 / 8
  chunk <- min(200, 256, floor(entries / (2048 * (1 * 4950 * 2 + 1 + 1))))
  expect_identical(chunk, 6)
  ctrl <- occ_ctrl(max_chunk_rows = chunk, max_retained_entries = entries, max_matrix_entries = entries,
    max_work_ops = 1e11, max_output_entries = entries)
  df <- expand.grid(time = 1:100, id = 1:200, KEEP.OUT.ATTRS = FALSE)
  df$g <- ifelse(df$id <= 100L, 51, Inf)
  df$y <- 0
  po <- api$construct_ipw_panel_owner(df[, c("id", "time", "y", "g")], "id", "time", "y", "g", ctrl)
  n_pairs <- po$counts$n_pairs
  expect_identical(n_pairs, 200 * 4950)
  # The retired repeated-merge charge, from its accepted formula (CP07 R/chained_pipeline_metadata.R):
  m_occ <- min(n_pairs, 2 * 4950); batches <- ceiling(n_pairs / chunk)
  merge_charge <- 16 * (1 + ceiling(log2(m_occ + chunk + 1))) * (n_pairs + batches * m_occ)
  expect_gt(merge_charge, 1e11)
  pre <- api$preflight_owner_sparse_occupancy(po, ctrl)
  expect_lt(pre$work_ops_bound, 1e11)
  expect_lt(pre$work_ops_bound, merge_charge / 100)
  occ <- api$derive_owner_sparse_occupancy(po, ctrl)
  d <- occ$diagnostics
  expect_identical(nrow(occ$occupancy), 9900L)
  expect_true(all(occ$occupancy$count == 100))
  expect_identical(d$n_distinct_histories, 2L)
  expect_identical(d$max_history_multiplicity, 100L)
  expect_identical(d$n_history_pair_visits, 9900)
  expect_identical(d$pair_visits_avoided, n_pairs - 9900)
  expect_identical(d$refinement_unit_visits, 200 * 100)
  expect_identical(d$n_rows, 198L)
  expect_lte(d$actual_work_ops, pre$work_ops_bound)
  expect_lte(d$actual_concurrent_entries, pre$concurrent_entries)
})

test_that("native interrupt checks run inside accumulation and an interrupt aborts without a sealed result", {
  fx <- ref_env$make_pipeline_metadata_fixtures()
  ctrl <- occ_ctrl()
  po <- api$construct_ipw_panel_owner(fx$rotating_gapped, "id", "time", "y", "g", ctrl)
  original <- api[["ipw_owner_check_interrupt"]]
  calls <- 0L
  was_locked <- bindingIsLocked("ipw_owner_check_interrupt", api)
  swap <- function(fn) {
    if (bindingIsLocked("ipw_owner_check_interrupt", api)) unlockBinding("ipw_owner_check_interrupt", api)
    assign("ipw_owner_check_interrupt", fn, envir = api)
    if (was_locked) lockBinding("ipw_owner_check_interrupt", api)
  }
  on.exit(swap(original), add = TRUE)
  swap(function() { calls <<- calls + 1L; original() })
  occ <- api$derive_owner_sparse_occupancy(po, ctrl)
  expect_identical(calls, as.integer(occ$diagnostics$n_refinement_positions + occ$diagnostics$n_rows))
  swap(function() {
    calls <<- calls + 1L
    if (calls >= 2L) stop(structure(class = c("interrupt", "condition"), list(message = "test interrupt", call = NULL)))
    original()
  })
  calls <- 0L
  res <- tryCatch(api$derive_owner_sparse_occupancy(po, ctrl), interrupt = function(e) "interrupted")
  expect_identical(res, "interrupted")
})

test_that("tampered owners and decorated inputs are refused before any accumulation", {
  fx <- ref_env$make_pipeline_metadata_fixtures()
  ctrl <- occ_ctrl()
  po <- api$construct_ipw_panel_owner(fx$balanced, "id", "time", "y", "g", ctrl)
  bad <- po
  bad$observations$calendar_index[1L] <- bad$observations$calendar_index[2L]
  tr <- new.env(parent = emptyenv()); tr$calls <- 42L
  expect_error(api$derive_owner_sparse_occupancy(bad, ctrl, trace_env = tr), class = "cdid_ipw_owner_integrity_error")
  expect_identical(tr$calls, 42L)
  bad2 <- po
  bad2$units$n_obs[1L] <- 1L
  expect_error(api$derive_owner_sparse_occupancy(bad2, ctrl), class = "cdid_ipw_owner_integrity_error")
})
