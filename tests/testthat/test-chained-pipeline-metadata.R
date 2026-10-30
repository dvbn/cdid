# Strict actual-installed namespace tests. Source audit uses a separate disclosed adapter.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("Declared digest test dependency required.", call. = FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api), environmentIsLocked(api))
api_binding <- list(mode = "installed_package", env = api)
ref_env <- new.env(parent = baseenv())
sys.source(testthat::test_path("fixtures", "pipeline_metadata_reference.R"), envir = ref_env)
metadata_binding_names <- c("cdid_pipeline_metadata_abort","pipeline_metadata_require_deps","pipeline_metadata_control","validate_pipeline_metadata_control","pipeline_metadata_resolve_owners","pipeline_metadata_owner_summary","preflight_owner_sparse_occupancy","preflight_pipeline_occupancy","pipeline_occupancy_accumulate","derive_owner_sparse_occupancy","construct_owner_sparse_occupancy","adapt_owner_to_sparse_occupancy","chained_pipeline_occupancy","validate_owner_sparse_occupancy","validate_pipeline_occupancy","pipeline_metadata_validate_catalog_inputs_early","preflight_owner_comparison_catalog","preflight_pipeline_metadata","pipeline_metadata_match_pair","pipeline_metadata_catalog_reservation","pipeline_metadata_resolve_catalog_mappings","build_owner_comparison_catalog","adapt_owner_to_comparison_catalog","construct_pipeline_metadata","chained_pipeline_metadata","validate_owner_comparison_catalog","validate_pipeline_metadata","pipeline_metadata_row_spec_preflight","extract_pipeline_catalog_row_spec","construct_pipeline_catalog_row_views")
default_owner_ctrl <- function(
  max_chunk_rows = 50L,
  max_retained_entries = 1e6,
  max_matrix_entries = 1e6,
  max_work_ops = 1e7,
  max_output_entries = 1e6
) {
  api$ipw_owner_control(
    max_chunk_rows = max_chunk_rows,
    max_retained_entries = max_retained_entries,
    max_matrix_entries = max_matrix_entries,
    max_work_ops = max_work_ops,
    max_output_entries = max_output_entries
  )
}

default_cat_bounds <- function(
  max_candidates = 1e5,
  max_rows = 1e5,
  max_pool_entries = 1e6,
  max_scratch = 1e7,
  max_work = 1e7
) {
  list(
    max_candidates = max_candidates,
    max_rows = max_rows,
    max_pool_entries = max_pool_entries,
    max_scratch = max_scratch,
    max_work = max_work
  )
}

test_that("Execution mode and native backend primitives are authentic (no constant stubs or fake namespace)", {
  expect_identical(api_binding$mode, "installed_package")
  expect_true(isNamespace(api))
  expect_true(environmentIsLocked(api))
  expect_true(all(vapply(metadata_binding_names, function(nm) identical(environment(api[[nm]]), api), logical(1))))
  expect_false(any(metadata_binding_names %in% getNamespaceExports("cdid")))
  backend <- api$gmm_native_backend_info()
  expect_true(is.list(backend))
  x_vec <- c(1.0, 2.0, 3.0)
  y_alias <- x_vec
  z_copy <- x_vec + 0.0
  expect_true(api$gmm_sexp_aliases(x_vec, y_alias))
  expect_false(api$gmm_sexp_aliases(x_vec, z_copy))
})

test_that("Root panel24.csv and expected24.json match derived sparse occupancy and comparison catalog exactly", {
  p24_path <- testthat::test_path("fixtures", "panel24.csv")
  e24_path <- testthat::test_path("fixtures", "expected24.json")
  expect_true(file.exists(p24_path))
  expect_true(file.exists(e24_path))

  if (requireNamespace("digest", quietly = TRUE)) {
    expect_identical(
      digest::digest(file = p24_path, algo = "sha256"),
      "8d6573cf08265e0ca5f210f96f9235604d5649c468fd285ea9118f0ad292effe"
    )
    expect_identical(
      digest::digest(file = e24_path, algo = "sha256"),
      "9f4f8953df3e23490ff165c0f4c514764b6958eb38ef5c10eb1736046f734dc5"
    )
  }

  p24 <- read.csv(p24_path, stringsAsFactors = FALSE)
  exp24 <- dget(testthat::test_path("fixtures", "expected24.R"))

  ctrl <- default_owner_ctrl(max_chunk_rows = 5L)
  cat_b <- default_cat_bounds()
  tr <- new.env(parent = emptyenv())

  po <- api$construct_ipw_panel_owner(
    data = p24, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, owner_id = "panel24"
  )
  u24 <- p24[!duplicated(p24$id), ]
  fo <- api$construct_ipw_frame_owner(
    panel_owner = po, id = u24$id, Z = cbind(1, u24$pattern, u24$sign),
    control = ctrl, frame_id = "frame24"
  )

  occ_obj <- api$derive_owner_sparse_occupancy(
    panel_owner = po, control = ctrl, frame_owner = fo,
    trace_env = tr, occupancy_id = "occ24"
  )
  api$validate_owner_sparse_occupancy(occ_obj, panel_owner = po, frame_owner = fo)

  # Verify against expected24.json
  expect_identical(po$n_a_units, as.integer(exp24$n_A))
  expect_identical(po$n_observed_rows, as.integer(exp24$n_observed))
  expect_equal(po$counts$n_pairs, as.double(exp24$n_unit_pairs))
  expect_equal(occ_obj$diagnostics$total_occupancy_count, as.double(exp24$n_unit_pairs))
  expect_identical(occ_obj$diagnostics$endpoints_read, 0L)
  expect_identical(tr$calls, 0L)
  expect_identical(tr$pos_s_index, integer(0L))
  expect_identical(tr$pos_t_index, integer(0L))

  exp_occ <- data.frame(
    cohort = as.double(exp24$occupancy$cohort),
    s = as.double(exp24$occupancy$s),
    t = as.double(exp24$occupancy$t),
    count = as.double(exp24$occupancy$count),
    stringsAsFactors = FALSE
  )
  expect_identical(occ_obj$occupancy, exp_occ)

  # Brute-force reference agreement on panel24
  bf24 <- ref_env$brute_force_panel_occupancy(p24, "id", "time", "y", "g")
  expect_identical(occ_obj$occupancy, bf24$sparse_occupancy)
  expect_identical(occ_obj$cell_indices, bf24$cell_indices)

  # Comparison catalog on panel24 (never-treated and not-yet-treated both yield the 6 original rows)
  meta_nt <- api$build_owner_comparison_catalog(
    panel_owner = po,
    targets = 3,
    anticipation = 0,
    control_group = "nevertreated",
    conventions = "conv:panel24:nt",
    control = ctrl,
    catalog_bounds = cat_b,
    frame_owner = fo,
    trace_env = tr,
    metadata_id = "meta24_nt"
  )
  api$validate_owner_comparison_catalog(meta_nt, panel_owner = po, frame_owner = fo)
  expect_identical(nrow(meta_nt$catalog$rows), 6L)
  expect_identical(nrow(meta_nt$catalog$operators), 1L)
  expect_true(all(meta_nt$catalog$rows$original))
  expect_identical(meta_nt$diagnostics$endpoints_read, 0L)
  expect_identical(tr$calls, 0L)

  # Exact catalog equivalence with direct build_comparison_catalog call
  direct_cat <- api$build_comparison_catalog(
    calendar = po$all_times,
    cohorts = po$cohorts,
    targets = 3,
    anticipation = 0,
    occupancy = exp_occ,
    control_group = "nevertreated",
    conventions = "conv:panel24:nt",
    requested_post = NULL,
    bounds = cat_b
  )
  expect_identical(meta_nt$catalog, direct_cat)

  # Verify compact outcome-free row views can be built for all 6 retained rows with 0 endpoint reads
  for (r_id in seq_len(nrow(meta_nt$catalog$rows))) {
    rv <- api$construct_pipeline_catalog_row_views(meta_nt, r_id)
    expect_identical(rv$endpoints_read, 0L)
    expect_s3_class(rv$nuisance_view, "cdid_ipw_nuisance_view")
    expect_s3_class(rv$link_view, "cdid_ipw_link_view")
  }
})

test_that("Brute-force pair counts and catalog equivalence hold across rotating, gapped, balanced, singleton, and empty-A fixtures", {
  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  ctrl <- default_owner_ctrl(max_chunk_rows = 4L)
  cat_b <- default_cat_bounds()

  for (nm in c("rotating_gapped", "balanced", "with_singletons")) {
    df <- fx_list[[nm]]
    bf <- ref_env$brute_force_panel_occupancy(df, "id", "time", "y", "g")
    tr <- new.env(parent = emptyenv())

    po <- api$construct_ipw_panel_owner(
      data = df, idname = "id", tname = "time", yname = "y", gname = "g",
      control = ctrl, owner_id = paste0("po_", nm)
    )
    occ_obj <- api$derive_owner_sparse_occupancy(po, ctrl, trace_env = tr)

    expect_identical(occ_obj$occupancy, bf$sparse_occupancy, info = nm)
    expect_identical(occ_obj$cell_indices, bf$cell_indices, info = nm)
    expect_identical(occ_obj$calendar, bf$calendar, info = nm)
    expect_identical(occ_obj$observed_times, bf$observed_times, info = nm)
    expect_identical(occ_obj$cohorts, bf$cohorts, info = nm)
    expect_identical(occ_obj$a_cohorts, bf$a_cohorts, info = nm)
    expect_identical(occ_obj$owner_summary$n_singleton_units, bf$n_singleton_units, info = nm)
    expect_identical(occ_obj$owner_summary$n_unobserved_units, bf$n_unobserved_units, info = nm)
    expect_identical(occ_obj$owner_summary$n_explicit_na_rows, bf$n_explicit_na_rows, info = nm)
    expect_identical(occ_obj$diagnostics$endpoints_read, 0L, info = nm)
    expect_identical(tr$calls, 0L, info = nm)

    # Singleton and unobserved metadata remain intact in the referenced owner
    ref_po <- get("panel_owner", envir = occ_obj$panel_ref, inherits = FALSE)
    expect_true(api$gmm_sexp_aliases(po, ref_po))
    expect_identical(nrow(ref_po$singletons), bf$n_singleton_units)
    expect_identical(nrow(ref_po$unobserved_units), bf$n_unobserved_units)

    for (cg in c("nevertreated", "notyettreated")) {
      meta <- api$build_owner_comparison_catalog(
        panel_owner = po,
        targets = c(3, 4),
        anticipation = 0,
        control_group = cg,
        conventions = paste0("conv:", nm, ":", cg),
        control = ctrl,
        catalog_bounds = cat_b,
        trace_env = tr
      )
      expect_identical(meta$diagnostics$endpoints_read, 0L)
      expect_identical(tr$calls, 0L)
      expect_true(meta$diagnostics$all_catalog_rows_retained)
      expect_true(meta$diagnostics$all_catalog_operators_retained)

      # Equivalence with direct catalog from sparse occupancy AND from dense occupancy (including 0-count cells)
      cat_sparse <- api$build_comparison_catalog(
        calendar = po$all_times,
        cohorts = po$cohorts,
        targets = c(3, 4),
        anticipation = 0,
        occupancy = bf$sparse_occupancy,
        control_group = cg,
        conventions = paste0("conv:", nm, ":", cg),
        requested_post = NULL,
        bounds = cat_b
      )
      cat_dense <- api$build_comparison_catalog(
        calendar = po$all_times,
        cohorts = po$cohorts,
        targets = c(3, 4),
        anticipation = 0,
        occupancy = bf$dense_occupancy,
        control_group = cg,
        conventions = paste0("conv:", nm, ":", cg),
        requested_post = NULL,
        bounds = cat_b
      )
      expect_identical(meta$catalog$rows, cat_sparse$rows)
      expect_identical(meta$catalog$rows, cat_dense$rows)
      expect_identical(meta$catalog$operators, cat_sparse$operators)
      expect_identical(meta$catalog$operators, cat_dense$operators)
      expect_identical(meta$catalog$pools, cat_sparse$pools)
      expect_identical(meta$catalog$pools, cat_dense$pools)

      # Verify row_mapping counts against brute-force reference
      bf_map <- ref_env$brute_force_catalog_mappings(bf, meta$catalog)
      expect_identical(meta$row_mapping$target_occupancy_cell_id, bf_map$target_occupancy_cell_id)
      expect_identical(meta$row_mapping$target_pair_count, bf_map$target_pair_count)
      expect_identical(meta$row_mapping$control_pool_pair_count, bf_map$control_pool_pair_count)
      expect_identical(meta$row_mapping$risk_pair_count, bf_map$risk_pair_count)
      expect_identical(meta$row_mapping$all_cohort_pair_count, bf_map$all_cohort_pair_count)
      expect_identical(meta$row_mapping$target_a_unit_count, bf_map$target_a_unit_count)
      expect_identical(meta$row_mapping$control_pool_a_unit_count, bf_map$control_pool_a_unit_count)
      expect_identical(meta$row_mapping$risk_a_unit_count, bf_map$risk_a_unit_count)
      expect_identical(meta$row_mapping$n_supported_controls, bf_map$n_supported_controls)
    }

    # Verify D13 fixed-horizon comparison retention under notyettreated on rotating_gapped
    if (identical(nm, "rotating_gapped")) {
      meta_nyt <- api$build_owner_comparison_catalog(
        panel_owner = po,
        targets = c(3, 4),
        anticipation = 0,
        control_group = "notyettreated",
        conventions = "conv:d13",
        control = ctrl,
        catalog_bounds = cat_b
      )
      # Target g = 3 at post t = 3 has pool {4, Inf}; at post t = 4 original pool is {Inf},
      # and D13 fixed-horizon re-evaluates pair (1, 2) or (2, 3) under pool {Inf} (`original == FALSE`).
      expect_true(any(meta_nyt$catalog$rows$original))
      expect_true(any(!meta_nyt$catalog$rows$original))
    }
  }

  # Empty-A case (only singletons and unobserved units; n_a_units == 0, n_pairs == 0)
  df_empty <- fx_list$empty_a
  bf_empty <- ref_env$brute_force_panel_occupancy(df_empty, "id", "time", "y", "g")
  tr_empty <- new.env(parent = emptyenv())
  po_empty <- api$construct_ipw_panel_owner(
    data = df_empty, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl, owner_id = "po_empty_a"
  )
  expect_identical(po_empty$n_a_units, 0L)
  expect_identical(po_empty$n_singleton_units, 2L)
  expect_identical(po_empty$n_unobserved_units, 1L)
  expect_equal(po_empty$counts$n_pairs, 0)

  occ_empty <- api$derive_owner_sparse_occupancy(po_empty, ctrl, trace_env = tr_empty)
  expect_identical(nrow(occ_empty$occupancy), 0L)
  expect_identical(occ_empty$occupancy, bf_empty$sparse_occupancy)
  expect_identical(occ_empty$cell_indices, bf_empty$cell_indices)
  expect_identical(occ_empty$diagnostics$endpoints_read, 0L)
  expect_identical(tr_empty$calls, 0L)

  meta_empty <- api$build_owner_comparison_catalog(
    panel_owner = po_empty,
    targets = c(3, 4),
    anticipation = 0,
    control_group = "notyettreated",
    conventions = "conv:empty_a",
    control = ctrl,
    catalog_bounds = cat_b,
    trace_env = tr_empty
  )
  expect_identical(nrow(meta_empty$catalog$rows), 0L)
  expect_identical(nrow(meta_empty$catalog$operators), 0L)
  expect_identical(nrow(meta_empty$row_mapping), 0L)
  expect_identical(nrow(meta_empty$operator_mapping), 0L)
  expect_identical(meta_empty$diagnostics$endpoints_read, 0L)
  expect_identical(tr_empty$calls, 0L)
})

test_that("Exact adjacent IEEE binary64 dates and cohorts never collide (dictionary indices prevent float-format collisions)", {
  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  df_adj <- fx_list$adjacent_binary64
  ctrl <- default_owner_ctrl(max_chunk_rows = 2L)
  cat_b <- default_cat_bounds()

  t1 <- 1.0
  t1_adj <- 1.0 + .Machine$double.eps
  g3 <- 3.0
  g3_adj <- 3.0 + 2 * .Machine$double.eps

  # Confirm default R string formatting collides these values, while binary64 inequality holds
  expect_false(t1 == t1_adj)
  expect_identical(as.character(t1), as.character(t1_adj))
  expect_false(g3 == g3_adj)
  expect_identical(as.character(g3), as.character(g3_adj))

  bf_adj <- ref_env$brute_force_panel_occupancy(df_adj, "id", "time", "y", "g")
  po_adj <- api$construct_ipw_panel_owner(
    data = df_adj, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl
  )
  expect_identical(po_adj$n_calendar_dates, 5L)
  expect_identical(length(po_adj$cohorts), 3L)

  occ_adj <- api$derive_owner_sparse_occupancy(po_adj, ctrl)
  expect_identical(occ_adj$occupancy, bf_adj$sparse_occupancy)
  expect_identical(occ_adj$cell_indices, bf_adj$cell_indices)
  expect_false(occ_adj$diagnostics$formatted_float_keys_used)
  expect_true(occ_adj$diagnostics$dictionary_index_keys_used)

  # Both adjacent dates (s_index = 1 and s_index = 2) and both adjacent cohorts (cohort_index = 1 and 2) are distinct in cell_indices
  expect_true(all(c(1L, 2L) %in% occ_adj$cell_indices$s_index))
  expect_true(all(c(1L, 2L, 3L) %in% occ_adj$cell_indices$cohort_index))

  meta_adj <- api$build_owner_comparison_catalog(
    panel_owner = po_adj,
    targets = c(g3, g3_adj),
    anticipation = 0,
    control_group = "notyettreated",
    conventions = "conv:adj_binary64",
    control = ctrl,
    catalog_bounds = cat_b
  )
  expect_identical(nrow(meta_adj$catalog$targets), 2L)
  expect_identical(meta_adj$catalog$spec$calendar, po_adj$all_times)
  expect_identical(meta_adj$catalog$spec$cohorts, po_adj$cohorts)
  expect_true(all(c(1L, 2L) %in% meta_adj$catalog$rows$target_index))

  # Literal stored dates, cohorts, and anticipation (including non-zero anticipation) are preserved bit-for-bit
  delta_exact <- 1
  meta_delta <- api$build_owner_comparison_catalog(
    panel_owner = po_adj,
    targets = c(g3, g3_adj),
    anticipation = delta_exact,
    control_group = "nevertreated",
    conventions = "conv:adj_delta",
    control = ctrl,
    catalog_bounds = cat_b,
    occupancy = occ_adj
  )
  expect_identical(meta_delta$catalog$spec$anticipation, as.double(delta_exact))
  expect_identical(meta_delta$catalog$spec$calendar, po_adj$all_times)
  expect_identical(meta_delta$catalog$spec$cohorts, po_adj$cohorts)
})

test_that("Character IDs (including multi-byte UTF-8) are supported; exact accumulation never reads ID values, so none are charged", {
  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  df_chr <- fx_list$char_ids
  df_num <- fx_list$rotating_gapped
  ctrl <- default_owner_ctrl(max_chunk_rows = 3L)
  cat_b <- default_cat_bounds()

  po_chr <- api$construct_ipw_panel_owner(
    data = df_chr, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl
  )
  po_num <- api$construct_ipw_panel_owner(
    data = df_num, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl
  )
  expect_identical(po_chr$id_type, "character")
  expect_true(po_chr$id_byte_metadata$max > 0)
  expect_true(po_chr$id_byte_metadata$total >= po_chr$id_byte_metadata$max)

  # CP08 P2: occupancy reads unit ranges and calendar indices only. The same
  # observation structure under character and numeric IDs has identical charges.
  pre_chr <- api$preflight_owner_sparse_occupancy(po_chr, ctrl)
  pre_num <- api$preflight_owner_sparse_occupancy(po_num, ctrl)
  expect_false(pre_chr$plan$unit_id_values_read)
  expect_false(pre_chr$plan$owner_pair_reader_used)
  for (nm in c("concurrent_entries", "work_ops_bound", "total_retained_entries", "slice_pairs")) {
    expect_identical(pre_chr[[nm]], pre_num[[nm]], info = nm)
  }

  occ_chr <- api$derive_owner_sparse_occupancy(po_chr, ctrl)
  bf_chr <- ref_env$brute_force_panel_occupancy(df_chr, "id", "time", "y", "g")
  expect_identical(occ_chr$occupancy, bf_chr$sparse_occupancy)
  expect_identical(occ_chr$cell_indices, bf_chr$cell_indices)

  meta_chr <- api$build_owner_comparison_catalog(
    panel_owner = po_chr,
    targets = c(3, 4),
    anticipation = 0,
    control_group = "notyettreated",
    conventions = "conv:char_id",
    control = ctrl,
    catalog_bounds = cat_b
  )
  expect_identical(meta_chr$diagnostics$endpoints_read, 0L)
  expect_true(meta_chr$diagnostics$all_catalog_rows_retained)
})

test_that("Occupancy and comparison catalog are invariant to input row permutations, ID relabellings, and chunk sizes", {
  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  df_base <- fx_list$with_singletons
  cat_b <- default_cat_bounds()

  ctrl_1 <- default_owner_ctrl(max_chunk_rows = 1L)
  po_base <- api$construct_ipw_panel_owner(
    data = df_base, idname = "id", tname = "time", yname = "y", gname = "g",
    control = ctrl_1
  )
  occ_1 <- api$derive_owner_sparse_occupancy(po_base, ctrl_1)
  meta_1 <- api$build_owner_comparison_catalog(
    panel_owner = po_base,
    targets = c(3, 4, 5),
    anticipation = 0,
    control_group = "notyettreated",
    conventions = "conv:invariance",
    control = ctrl_1,
    catalog_bounds = cat_b,
    requested_post = data.frame(g = c(3, 4), t = c(3, 4))
  )

  # Vary max_chunk_rows across 1, 2, 3, 5, 11, 100
  for (ch_sz in c(2L, 3L, 5L, 11L, 100L)) {
    ctrl_k <- default_owner_ctrl(max_chunk_rows = ch_sz)
    occ_k <- api$derive_owner_sparse_occupancy(po_base, ctrl_k)
    expect_identical(occ_k$occupancy, occ_1$occupancy, info = paste("chunk_size =", ch_sz))
    expect_identical(occ_k$cell_indices, occ_1$cell_indices, info = paste("chunk_size =", ch_sz))

    meta_k <- api$build_owner_comparison_catalog(
      panel_owner = po_base,
      targets = c(3, 4, 5),
      anticipation = 0,
      control_group = "notyettreated",
      conventions = "conv:invariance",
      control = ctrl_k,
      catalog_bounds = cat_b,
      requested_post = data.frame(g = c(3, 4), t = c(3, 4))
    )
    expect_identical(meta_k$catalog$rows, meta_1$catalog$rows)
    expect_identical(meta_k$catalog$operators, meta_1$catalog$operators)
    expect_identical(meta_k$catalog$pools, meta_1$catalog$pools)
    expect_identical(meta_k$row_mapping, meta_1$row_mapping)
    expect_identical(meta_k$operator_mapping, meta_1$operator_mapping)
  }

  # Permute input rows and apply order-reversing unit ID bijections (both numeric and character)
  set.seed(20261015L)
  for (perm_i in 1:4) {
    perm_idx <- sample.int(nrow(df_base))
    df_perm <- df_base[perm_idx, , drop = FALSE]
    rownames(df_perm) <- NULL
    if (perm_i == 2L) {
      # Reverse numeric IDs; cohort and time observations per unit remain identical so occupancy is invariant
      df_perm$id <- 1000L - df_perm$id
    } else if (perm_i == 3L) {
      # Map to character IDs in non-monotone order
      id_map <- sprintf("id_perm_%03d", 900L - df_perm$id)
      df_perm$id <- id_map
    }
    ctrl_p <- default_owner_ctrl(max_chunk_rows = perm_i)
    po_p <- api$construct_ipw_panel_owner(
      data = df_perm, idname = "id", tname = "time", yname = "y", gname = "g",
      control = ctrl_p
    )
    occ_p <- api$derive_owner_sparse_occupancy(po_p, ctrl_p)
    expect_identical(occ_p$occupancy, occ_1$occupancy, info = paste("perm =", perm_i))
    expect_identical(occ_p$cell_indices, occ_1$cell_indices, info = paste("perm =", perm_i))

    meta_p <- api$build_owner_comparison_catalog(
      panel_owner = po_p,
      targets = c(3, 4, 5),
      anticipation = 0,
      control_group = "notyettreated",
      conventions = "conv:invariance",
      control = ctrl_p,
      catalog_bounds = cat_b,
      requested_post = data.frame(g = c(3, 4), t = c(3, 4))
    )
    expect_identical(meta_p$catalog$rows, meta_1$catalog$rows)
    expect_identical(meta_p$catalog$operators, meta_1$catalog$operators)
    expect_identical(meta_p$row_mapping, meta_1$row_mapping)
    expect_identical(meta_p$operator_mapping, meta_1$operator_mapping)
  }
})

test_that("Outcome-free traversal holds even when raw outcomes are mutated in a separate owner", {
  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  df1 <- fx_list$rotating_gapped
  df2 <- df1
  obs_idx <- which(!is.na(df2$y))
  df2$y[obs_idx] <- df2$y[obs_idx] * (-17.5) + 42.125

  ctrl <- default_owner_ctrl(max_chunk_rows = 2L)
  cat_b <- default_cat_bounds()
  tr1 <- new.env(parent = emptyenv())
  tr2 <- new.env(parent = emptyenv())

  po1 <- api$construct_ipw_panel_owner(df1, "id", "time", "y", "g", ctrl)
  po2 <- api$construct_ipw_panel_owner(df2, "id", "time", "y", "g", ctrl)

  meta1 <- api$build_owner_comparison_catalog(
    panel_owner = po1, targets = c(3, 4), anticipation = 0,
    control_group = "notyettreated", conventions = "conv:y_free",
    control = ctrl, catalog_bounds = cat_b, trace_env = tr1
  )
  meta2 <- api$build_owner_comparison_catalog(
    panel_owner = po2, targets = c(3, 4), anticipation = 0,
    control_group = "notyettreated", conventions = "conv:y_free",
    control = ctrl, catalog_bounds = cat_b, trace_env = tr2
  )

  expect_identical(tr1$calls, 0L)
  expect_identical(tr2$calls, 0L)
  expect_identical(meta1$occupancy, meta2$occupancy)
  expect_identical(meta1$cell_indices, meta2$cell_indices)
  expect_identical(meta1$catalog, meta2$catalog)
  expect_identical(meta1$row_mapping, meta2$row_mapping)
  expect_identical(meta1$operator_mapping, meta2$operator_mapping)
})

test_that("Preflight resource bounds and input validations refuse before callbacks without fallback, retry, or dropped rows", {
  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  df <- fx_list$char_ids
  ctrl_ok <- default_owner_ctrl(max_chunk_rows = 3L)
  cat_ok <- default_cat_bounds()
  po <- api$construct_ipw_panel_owner(df, "id", "time", "y", "g", ctrl_ok)

  pre_ok <- api$preflight_owner_sparse_occupancy(po, ctrl_ok)
  expect_equal(
    pre_ok$concurrent_entries,
    pre_ok$retained_vector_entries + pre_ok$header_entries +
      max(pre_ok$grouping_scratch_entries, pre_ok$accumulation_scratch_entries, pre_ok$assembly_scratch_entries)
  )
  expect_equal(pre_ok$work_ops_bound, pre_ok$grouping_work_ops + pre_ok$accumulation_work_ops + pre_ok$assembly_work_ops)

  # 1. Missing explicit controls refuse immediately
  expect_error(
    api$pipeline_metadata_control(ctrl_ok),
    class = "cdid_pipeline_metadata_control_error"
  )
  expect_error(
    api$derive_owner_sparse_occupancy(po),
    class = "cdid_pipeline_metadata_input_error"
  )
  expect_error(
    api$build_owner_comparison_catalog(
      panel_owner = po, targets = c(3, 4), anticipation = 0,
      control_group = "notyettreated", conventions = "conv:test",
      control = ctrl_ok
    ),
    class = "cdid_catalog_bound_error"
  )

  # 2. Owner control bound refusals happen in preflight BEFORE trace_env is reset by traversal
  for (field in c("max_retained_entries", "max_output_entries", "max_matrix_entries", "max_work_ops")) {
    ctrl_tight <- ctrl_ok
    ctrl_tight[[field]] <- 10
    tr_sentinel <- new.env(parent = emptyenv())
    tr_sentinel$calls <- 999L
    expect_error(
      api$derive_owner_sparse_occupancy(po, ctrl_tight, trace_env = tr_sentinel),
      class = "cdid_ipw_owner_bound_error",
      info = field
    )
    # Accumulation never started because preflight refused before the trace reset
    expect_identical(tr_sentinel$calls, 999L, info = field)
  }

  # 3. Exact caps of the executed algorithm. One work unit below its bound refuses
  #    before the trace reset; at the bound it is admitted. The smallest feasible
  #    matrix cap forces the minimal slice budget (T - 1 pairs) and still yields
  #    identical exact counts; one entry less refuses before the trace reset.
  bf_ok <- ref_env$brute_force_panel_occupancy(df, "id", "time", "y", "g")
  tr_id <- new.env(parent = emptyenv())
  tr_id$calls <- 777L
  ctrl_work <- ctrl_ok
  ctrl_work$max_work_ops <- pre_ok$work_ops_bound - 1
  expect_error(
    api$derive_owner_sparse_occupancy(po, ctrl_work, trace_env = tr_id),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(tr_id$calls, 777L)
  ctrl_work$max_work_ops <- pre_ok$work_ops_bound
  occ_work <- api$derive_owner_sparse_occupancy(po, ctrl_work, trace_env = tr_id)
  expect_identical(tr_id$calls, 0L)
  expect_identical(occ_work$occupancy, bf_ok$sparse_occupancy)
  expect_lte(occ_work$diagnostics$actual_work_ops, pre_ok$work_ops_bound)
  t_span <- max(1, po$n_calendar_dates - 1)
  min_me <- pre_ok$retained_vector_entries + pre_ok$header_entries +
    max(pre_ok$grouping_scratch_entries, pre_ok$accumulation_scratch_entries - 14 * pre_ok$slice_pairs + 14 * t_span,
        pre_ok$assembly_scratch_entries)
  ctrl_me <- ctrl_ok
  ctrl_me$max_matrix_entries <- min_me - 1
  tr_id$calls <- 777L
  expect_error(
    api$derive_owner_sparse_occupancy(po, ctrl_me, trace_env = tr_id),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(tr_id$calls, 777L)
  ctrl_me$max_matrix_entries <- min_me
  occ_me <- api$derive_owner_sparse_occupancy(po, ctrl_me)
  expect_identical(as.double(occ_me$preflight$slice_pairs), as.double(t_span))
  expect_identical(occ_me$occupancy, bf_ok$sparse_occupancy)
  expect_identical(occ_me$cell_indices, bf_ok$cell_indices)
  expect_lte(occ_me$diagnostics$actual_concurrent_entries, min_me)

  # 4. Early catalog bound and input refusals in build_owner_comparison_catalog happen BEFORE pair traversal
  for (bad_cat in list(
    list(bounds = modifyList(cat_ok, list(max_scratch = 10)), cls = "cdid_catalog_bound_error"),
    list(bounds = modifyList(cat_ok, list(max_work = 10)), cls = "cdid_catalog_bound_error"),
    list(bounds = cat_ok, targets = 999, cls = "cdid_catalog_target_error"),
    list(bounds = cat_ok, control_group = "invalid_pool", cls = "cdid_catalog_input_error"),
    list(bounds = cat_ok, conventions = "", cls = "cdid_catalog_input_error")
  )) {
    tr_cat <- new.env(parent = emptyenv())
    tr_cat$calls <- 555L
    expect_error(
      api$build_owner_comparison_catalog(
        panel_owner = po,
        targets = if (!is.null(bad_cat$targets)) bad_cat$targets else c(3, 4),
        anticipation = 0,
        control_group = if (!is.null(bad_cat$control_group)) bad_cat$control_group else "notyettreated",
        conventions = if (!is.null(bad_cat$conventions)) bad_cat$conventions else "conv:test",
        control = ctrl_ok,
        catalog_bounds = bad_cat$bounds,
        trace_env = tr_cat
      ),
      class = bad_cat$cls
    )
    expect_identical(tr_cat$calls, 555L)
  }

  # 5. Catalog candidate/row/pool_entries caps refuse cleanly through unchanged build_comparison_catalog
  for (cap_nm in c("max_candidates", "max_rows", "max_pool_entries")) {
    b_tight <- modifyList(cat_ok, setNames(list(1), cap_nm))
    expect_error(
      api$build_owner_comparison_catalog(
        panel_owner = po, targets = c(3, 4), anticipation = 0,
        control_group = "notyettreated", conventions = "conv:cap",
        control = ctrl_ok, catalog_bounds = b_tight
      ),
      class = "cdid_catalog_bound_error",
      info = cap_nm
    )
  }

  # 6. Tampered / unsealed / mismatched owner or occupancy objects are refused via native SEXP alias check
  occ_valid <- api$derive_owner_sparse_occupancy(po, ctrl_ok)
  occ_tampered <- occ_valid
  occ_tampered$occupancy$count[1L] <- occ_tampered$occupancy$count[1L] + 1
  expect_error(
    api$validate_owner_sparse_occupancy(occ_tampered),
    class = "cdid_pipeline_metadata_integrity_error"
  )

  po_other <- api$construct_ipw_panel_owner(fx_list$balanced, "id", "time", "y", "g", ctrl_ok)
  expect_error(
    api$validate_owner_sparse_occupancy(occ_valid, panel_owner = po_other),
    class = "cdid_pipeline_metadata_owner_mismatch_error"
  )
})

test_that("Discriminating source mutants in R/chained_pipeline_metadata.R are detected and killed", {
  expect_true(all(vapply(metadata_binding_names, function(nm) is.function(api[[nm]]), logical(1))))
  orig_lines <- unlist(lapply(metadata_binding_names, function(nm) {
    c(paste0(nm, " <- "), deparse(api[[nm]], width.cutoff = 500L))
  }), use.names = FALSE)

  make_mutant_env <- function(pattern, replacement) {
    mut_lines <- sub(pattern, replacement, orig_lines, fixed = TRUE)
    expect_false(identical(mut_lines, orig_lines), info = paste("Pattern not found:", pattern))
    m_env <- new.env(parent = api)
    eval(parse(text = mut_lines), envir = m_env)
    m_env
  }

  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  ctrl <- default_owner_ctrl(max_chunk_rows = 3L)
  cat_b <- default_cat_bounds()

  # Mutant 1: an outcome endpoint read inside derive_owner_sparse_occupancy is
  # visible in the caller trace, which the actual code leaves at zero reads.
  m1 <- make_mutant_env(
    "ipw_owner_reset_trace(trace_env)",
    "ipw_owner_reset_trace(trace_env); ipw_owner_read_y_positions(p_own, 1, 2, trace_env = trace_env)"
  )
  po_rot <- api$construct_ipw_panel_owner(fx_list$rotating_gapped, "id", "time", "y", "g", ctrl)
  tr_m1 <- new.env(parent = emptyenv())
  invisible(m1$derive_owner_sparse_occupancy(po_rot, ctrl, trace_env = tr_m1))
  expect_identical(tr_m1$calls, 1L)
  invisible(api$derive_owner_sparse_occupancy(po_rot, ctrl, trace_env = tr_m1))
  expect_identical(tr_m1$calls, 0L)

  # Mutant 2: replace exact calendar dictionary indices by formatted-float matching.
  # Adjacent binary64 dates (1 vs 1 + 2^-52) collapse into one index and the
  # resulting (s, s) cell is refused; the actual code matches brute force.
  m2 <- make_mutant_env(
    "cal_idx <- p_own$observations$calendar_index",
    "cal_idx <- match(as.character(p_own$observations$time), as.character(p_own$all_times))"
  )
  po_adj <- api$construct_ipw_panel_owner(fx_list$adjacent_binary64, "id", "time", "y", "g", ctrl)
  bf_adj <- ref_env$brute_force_panel_occupancy(fx_list$adjacent_binary64, "id", "time", "y", "g")
  expect_error(m2$derive_owner_sparse_occupancy(po_adj, ctrl), class = "cdid_pipeline_metadata_internal_error")
  expect_identical(api$derive_owner_sparse_occupancy(po_adj, ctrl)$cell_indices, bf_adj$cell_indices)

  # Mutant 3: Drop D13 fixed-horizon comparisons after catalog construction
  m3 <- make_mutant_env(
    "mappings <- pipeline_metadata_resolve_catalog_mappings(",
    "catalog$rows <- catalog$rows[catalog$rows$original, , drop = FALSE]; mappings <- pipeline_metadata_resolve_catalog_mappings("
  )
  meta_true <- api$build_owner_comparison_catalog(
    panel_owner = po_rot, targets = c(3, 4), anticipation = 0,
    control_group = "notyettreated", conventions = "conv:m3",
    control = ctrl, catalog_bounds = cat_b
  )
  meta_m3 <- m3$build_owner_comparison_catalog(
    panel_owner = po_rot, targets = c(3, 4), anticipation = 0,
    control_group = "notyettreated", conventions = "conv:m3",
    control = ctrl, catalog_bounds = cat_b
  )
  expect_true(nrow(meta_m3$catalog$rows) < nrow(meta_true$catalog$rows))
  expect_false(any(!meta_m3$catalog$rows$original))

  # Mutants 4 and 7 must enter the general history path. Each cohort has two
  # copies of each of the histories (1,2,4) and (1,3,4); the union calendar has
  # four dates, so no unit has a complete roster. No archived fixture changes.
  gapped <- data.frame(id = rep(1:8, each = 3L),
    time = rep(c(1, 2, 4, 1, 3, 4), 4L),
    g = rep(rep(c(3, Inf), each = 4L), each = 3L), y = as.double(1:24))
  po_gap <- api$construct_ipw_panel_owner(gapped, "id", "time", "y", "g", ctrl)
  bf_gap <- ref_env$brute_force_panel_occupancy(gapped, "id", "time", "y", "g")
  expect_true(all(po_gap$units$n_obs == 3L))
  expect_identical(po_gap$n_calendar_dates, 4L)
  expect_identical(api$derive_owner_sparse_occupancy(po_gap, ctrl)$cell_indices, bf_gap$cell_indices)

  # Mutant 4: count each distinct history once instead of its two units.
  # This represents 12 pairs instead of the panel's 24, so the total refuses.
  m4 <- make_mutant_env(
    "w_hist <- as.double(diff(c(first, n_a + 1L)))",
    "w_hist <- rep(1, length(first))"
  )
  expect_error(m4$derive_owner_sparse_occupancy(po_gap, ctrl), class = "cdid_pipeline_metadata_count_error")

  # Mutant 5: omit the pair-slice scratch from the concurrent reservation.
  m5 <- make_mutant_env(
    "ipw_owner_safe_product(14, slice_pairs, max_me, \"pair slice expansion and radix collapse scratch\")",
    "0"
  )
  pre_m5 <- m5$preflight_owner_sparse_occupancy(po_rot, ctrl)
  pre_true <- api$preflight_owner_sparse_occupancy(po_rot, ctrl)
  expect_false(identical(pre_m5$concurrent_entries, pre_true$concurrent_entries))

  # Mutant 6: label by (length, cohort) only, without positional refinement. Equal-length
  # histories of one cohort merge; the total still agrees but cells differ from brute force.
  m6 <- make_mutant_env("key <- cal_idx[sts[seq_len(p)] + (j - 1)]", "key <- rep(1L, p)")
  bf_rot <- ref_env$brute_force_panel_occupancy(fx_list$rotating_gapped, "id", "time", "y", "g")
  occ_m6 <- m6$derive_owner_sparse_occupancy(po_rot, ctrl)
  expect_false(identical(occ_m6$cell_indices, bf_rot$cell_indices))
  expect_identical(api$derive_owner_sparse_occupancy(po_rot, ctrl)$cell_indices, bf_rot$cell_indices)

  # Mutant 7: drop the cohort from the initial history label; identical histories of
  # different cohorts merge and their pairs move to one cohort.
  m7 <- make_mutant_env("nl[-1L] != nl[-n_a] | ks[-1L] != ks[-n_a]", "nl[-1L] != nl[-n_a]")
  occ_m7 <- m7$derive_owner_sparse_occupancy(po_gap, ctrl)
  expect_false(identical(occ_m7$cell_indices, bf_gap$cell_indices))
  expect_identical(api$derive_owner_sparse_occupancy(po_gap, ctrl)$cell_indices, bf_gap$cell_indices)

  # Mutant 8: never clear the dense row buffer between (cohort, s) rows.
  m8 <- make_mutant_env("buf[ut] <- 0", "invisible(NULL)")
  res_m8 <- tryCatch(m8$derive_owner_sparse_occupancy(po_rot, ctrl), error = identity)
  expect_true(inherits(res_m8, "error") || !identical(res_m8$cell_indices, bf_rot$cell_indices))
})

test_that("Frame owner resolution, design bridge, aliases, decorative attributes, and cutoff overflow preflight work as specified", {
  fx_list <- ref_env$make_pipeline_metadata_fixtures()
  df <- fx_list$balanced
  ctrl <- default_owner_ctrl(max_chunk_rows = 4L)
  cat_b <- default_cat_bounds()

  po <- api$construct_ipw_panel_owner(df, "id", "time", "y", "g", ctrl)
  u_df <- df[!duplicated(df$id), ]
  # Include an external frame-only unit (id = 99) with allow_frame_only = TRUE
  fo <- api$construct_ipw_frame_owner(
    panel_owner = po,
    id = c(u_df$id, 99L),
    Z = cbind(1, c(u_df$x, 0.0)),
    weights = c(rep(1.5, nrow(u_df)), 2.0),
    control = ctrl,
    frame_g = c(u_df$g, Inf),
    allow_frame_only = TRUE
  )

  # Passing `fo` directly as `panel_owner` resolves both panel_owner and frame_owner
  occ_from_fo <- api$adapt_owner_to_sparse_occupancy(fo, ctrl)
  api$validate_pipeline_occupancy(occ_from_fo, panel_owner = po, frame_owner = fo)
  expect_identical(occ_from_fo$owner_summary$n_frame_only_units, 1L)

  pm_ctrl <- api$pipeline_metadata_control(ctrl, cat_b)
  meta_from_fo <- api$adapt_owner_to_comparison_catalog(
    panel_owner = fo,
    targets = c(3, 4),
    anticipation = 0,
    control_group = "notyettreated",
    conventions = "conv:fo_bridge",
    control = pm_ctrl,
    declared_qualified = data.frame(g = c(3, 3, 4), t = c(3, 4, 4)),
    target_post_times = data.frame(g = c(3, 3, 4), t = c(3, 4, 4)),
    build_design = TRUE
  )
  api$validate_pipeline_metadata(meta_from_fo, panel_owner = po, frame_owner = fo)
  expect_false(is.null(meta_from_fo$design_links))
  expect_false(is.null(meta_from_fo$design))
  expect_true(all(!is.na(meta_from_fo$operator_mapping$target_a_weight)))
  expect_true(all(!is.na(meta_from_fo$operator_mapping$control_pool_a_weight)))
  expect_equal(meta_from_fo$operator_mapping$target_a_weight, rep(3.0, nrow(meta_from_fo$operator_mapping)))

  # Decorative / foreign attributes are refused before traversal
  bad_targets <- c(3, 4)
  attr(bad_targets, "certified") <- TRUE
  expect_error(
    api$build_owner_comparison_catalog(
      panel_owner = po, targets = bad_targets, anticipation = 0,
      control_group = "notyettreated", conventions = "conv:bad_attr",
      control = ctrl, catalog_bounds = cat_b
    ),
    class = "cdid_ipw_owner_foreign_attribute_error"
  )

  # Cutoff arithmetic overflow (calendar near DBL_MAX + large anticipation) refuses in preflight BEFORE pair traversal
  df_huge <- data.frame(
    id = c(1L, 1L, 2L, 2L),
    time = c(1.0, 1e308, 1.0, 1e308),
    y = c(1.0, 2.0, 3.0, 4.0),
    g = c(1e308, 1e308, Inf, Inf)
  )
  po_huge <- api$construct_ipw_panel_owner(df_huge, "id", "time", "y", "g", ctrl)
  tr_huge <- new.env(parent = emptyenv())
  tr_huge$calls <- 333L
  expect_error(
    api$build_owner_comparison_catalog(
      panel_owner = po_huge,
      targets = 1e308,
      anticipation = 1e308,
      control_group = "nevertreated",
      conventions = "conv:overflow",
      control = ctrl,
      catalog_bounds = cat_b,
      trace_env = tr_huge
    ),
    class = "cdid_design_arithmetic_error"
  )
  expect_identical(tr_huge$calls, 333L)
})


# Focused corrections: exact tuple identity, frame binding, bounded scalar queries,
# and coexistence reservations. Instrumentation is caller/test-owned and restored.
measure_metadata_calls <- function(name, run) {
  original <- api[[name]]; calls <- 0L
  was_locked <- bindingIsLocked(name,api)
  instrumented <- function(...) {calls <<- calls+1L; original(...)}
  if(was_locked) unlockBinding(name,api)
  assign(name,instrumented,envir=api)
  if(was_locked) lockBinding(name,api)
  on.exit({
    if(bindingIsLocked(name,api)) unlockBinding(name,api)
    assign(name,original,envir=api)
    if(was_locked) lockBinding(name,api)
  },add=TRUE)
  value <- run()
  list(value=value,calls=calls)
}

test_that("corrected metadata identity and resources discriminate the retained original faults", {
  d <- read.csv(testthat::test_path("fixtures","panel24.csv"),stringsAsFactors=FALSE)
  ctrl <- default_owner_ctrl(max_chunk_rows=5L);b <- default_cat_bounds()
  p <- api$construct_ipw_panel_owner(d,"id","time","y","g",ctrl)
  u <- d[!duplicated(d$id),]
  f1 <- api$construct_ipw_frame_owner(p,u$id,cbind(1,u$pattern),weights=rep(1,nrow(u)),control=ctrl,frame_id="f1")
  f2 <- api$construct_ipw_frame_owner(p,u$id,cbind(1,u$pattern),weights=rep(2,nrow(u)),control=ctrl,frame_id="f2")
  o <- api$derive_owner_sparse_occupancy(p,ctrl)
  m <- api$build_owner_comparison_catalog(p,3,0,"nevertreated","fault",control=ctrl,catalog_bounds=b,occupancy=o,frame_owner=f1)
  expect_identical(ncol(m$row_mapping),30L)
  expect_identical(ncol(m$operator_mapping),13L)
  expect_true(m$owner_summary$has_frame_owner)
  expect_error(api$validate_owner_comparison_catalog(m,frame_owner=f2),class="cdid_ipw_owner_mismatch_error")
  q <- measure_metadata_calls("validate_catalog_occupancy",function() api$extract_pipeline_catalog_row_spec(m,1))
  expect_identical(q$calls,0L)
  expect_identical(q$value$row_id,1L)
  T <- 2^27;s <- T-3
  expect_identical((s-1)*T+(T-1),(s-1)*T+T) # preserved collision witness, no large calendar allocation
  expect_identical(api$pipeline_metadata_match_pair(s,T-1,c(s,s),c(T-1,T)),1L)
  expect_identical(api$pipeline_metadata_match_pair(s,T,c(s,s),c(T-1,T)),2L)
  for(cap in c(528,758)) {
    small<-ctrl;small$max_retained_entries<-cap;small$max_output_entries<-cap
    check <- measure_metadata_calls("build_comparison_catalog",function() tryCatch(
      api$build_owner_comparison_catalog(p,3,0,"nevertreated","fault",control=small,catalog_bounds=b,occupancy=o,frame_owner=f1),error=identity))
    expect_s3_class(check$value,"cdid_ipw_owner_bound_error")
    expect_identical(check$calls,0L)
  }
  for(field in c("max_matrix_entries","max_work_ops")) {
    low<-ctrl;low[[field]]<-1
    early <- measure_metadata_calls("validate_calendar_spec",function() tryCatch(
      api$preflight_owner_comparison_catalog(p,3,0,"nevertreated","fault",control=low,catalog_bounds=b,occupancy=o),error=identity))
    expect_s3_class(early$value,"cdid_ipw_owner_bound_error")
    expect_identical(early$calls,0L)
  }
  tiny<-ctrl;tiny$max_matrix_entries<-1
  check <- measure_metadata_calls("extract_pipeline_catalog_row_spec",function() tryCatch(api$construct_pipeline_catalog_row_views(m,1,control=tiny),error=identity))
  expect_s3_class(check$value,"cdid_ipw_owner_bound_error")
  expect_identical(check$calls,0L)
  read <- measure_metadata_calls("ipw_owner_read_y_positions",function() api$derive_owner_sparse_occupancy(p,ctrl))
  expect_identical(read$calls,0L)
  expect_identical(read$value$occupancy,o$occupancy)
})
