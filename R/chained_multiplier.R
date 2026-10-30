# SPDX-License-Identifier: GPL-3.0-only
#
# Copyright (C) 2026 cdid authors
# License: GPL-3
#
# Private bounded arithmetic component for CP02 complete-target influence-
# function one-way cluster covariance and supplied multiplier statistics.
#
# Attribution and convention reference:
#   The statistical convention implemented here -- aggregating unit-level
#   influence-function rows to cluster sums S_c, scaling raw multiplier draws as
#   bres = C^{-1/2} sum_c eta_c S_c_tilde, estimating raw-draw scale via the
#   normal-consistent empirical type-1 interquartile range
#   (q_{0.75} - q_{0.25}) / (Phi^{-1}(0.75) - Phi^{-1}(0.25)) with order
#   statistic index ceil(B * p), converting to estimator standard errors via
#   sqrt(C) / n, and computing simultaneous sup-t critical values as the type-1
#   empirical (1 - alpha) quantile of max_j |bres_{b,j}| / bSigma_j -- follows
#   Callaway and Sant'Anna (2021) <doi:10.1016/j.jeconom.2020.12.001>, the
#   `did` package reference `mboot.R` (Brantly Callaway and Pedro H. C.
#   Sant'Anna), and the accepted CP02 cluster inference note
#   (`docs/theory/cp02_cluster_inference.tex`, thm:cov and thm:rad) together
#   with the CP02 numerical propagation note
#   (`docs/theory/cp02_ipw_gmm_numerical_propagation.tex`, sec:cov).
#   This file is an independent clean-room GPL-3 implementation from the
#   accepted equations with bounded streaming memory, explicit identity
#   alignment, and explicit refusal of nonpositive/nonfinite scales across the
#   complete fixed target set. Ordinary `did` delegation remains unchanged.
#
# Private signature and required fields:
#   chained_cluster_multiplier(
#     target_ids,
#     expected_unit_ids,
#     unit_cluster_map,
#     cluster_ids,
#     n_denom,
#     center,
#     if_batch_fn,
#     n_if_batches,
#     sign_batch_fn,
#     n_sign_batches,
#     expected_draw_ids,
#     n_draws,
#     alpha,
#     return_raw_draws,
#     budget
#   )
#
# All 15 arguments are required (no default values, no hidden global options,
# no RNG or seed usage, no working-directory or environment state).
#
# Required `unit_cluster_map` fields:
#   list or data.frame with components `unit_ids` and `cluster_ids` (same ID
#   storage type as `expected_unit_ids` and `cluster_ids`, bijective on
#   `expected_unit_ids`, mapping every expected unit to one declared cluster in
#   `cluster_ids`).
#
# Required `if_batch_fn(batch_index)` return fields:
#   list with components:
#     * `unit_ids`: atomic ID vector of length r_b >= 1 matching type of
#       `expected_unit_ids`
#     * `target_ids`: atomic ID vector of length J identical to `target_ids`
#     * `if_matrix`: double matrix of dimension r_b x J with complete target IF
#       rows
#   Optional component `cluster_ids` (if present, must match `unit_cluster_map`
#   for every unit in the batch). Weight fields (`weights`, `unit_weights`,
#   `cluster_weights`, `sampling_weights`) are refused to prevent double
#   weighting.
#
# Required `sign_batch_fn(batch_index)` return fields:
#   list with components:
#     * `draw_ids`: atomic ID vector of length B_b >= 1 matching the next slice
#       of `expected_draw_ids`
#     * `cluster_ids`: atomic ID vector of length C identical to `cluster_ids`
#     * `sign_matrix`: double or integer matrix of dimension B_b x C with all
#       entries in {-1, +1}
#
# Required `budget` fields (all explicit positive finite whole numbers):
#   `max_units`, `max_clusters`, `max_targets`, `max_draws`,
#   `max_if_batch_rows`, `max_if_batch_cells`, `max_sign_batch_rows`,
#   `max_sign_batch_cells`, `max_workspace_cells`, `max_workspace_bytes`.
# Optional `budget$workers` (CP08, inherited `pl`/`cores`): the declared number
# of worker processes for the raw-draw products. Absent or 1 is the serial
# path. With w > 1 the caller's sign batches are still requested one at a time,
# in draw order, from this (parent) process and validated exactly as before;
# rounds of at most w validated batches are multiplied by the same centered
# cluster sums in worker processes (fork on Unix-alikes, PSOCK on Windows,
# base `parallel` only). Each product is the serial batch product, so every
# draw, quantile and critical value equals the serial result; workers draw no
# random numbers. Their whole buffers are counted in the workspace preflight.

multiplier_abort <- function(message, subclass) {
  cdid_abort(
    message = message,
    class = unique(c(subclass, "cdid_multiplier_error"))
  )
}

validate_whole_number_limit <- function(x, name) {
  if (!is.numeric(x) || is.object(x) || !is.null(dim(x)) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x) || x > 9007199254740991) {
    multiplier_abort(
      sprintf("Budget field '%s' must be a positive finite whole number <= 2^53 - 1.", name),
      c("cdid_multiplier_budget_error", "cdid_multiplier_input_error")
    )
  }
  as.double(x)
}

validate_positive_int_scalar <- function(x, name, subclass = "cdid_multiplier_input_error",
                                         matrix_dimension = FALSE) {
  if (!is.numeric(x) || is.object(x) || !is.null(dim(x)) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x) || x > 9007199254740991) {
    multiplier_abort(
      sprintf("'%s' must be a positive finite whole number.", name),
      subclass
    )
  }
  if (matrix_dimension && x > .Machine$integer.max) {
    multiplier_abort(
      sprintf("'%s' exceeds genuine R matrix integer dimension limit (.Machine$integer.max).", name),
      c("cdid_multiplier_budget_error", subclass)
    )
  }
  if (x <= .Machine$integer.max) as.integer(x) else as.double(x)
}

validate_id_vector <- function(x, name, dup_subclass, align_subclass, scan_values = TRUE) {
  if (!(is.character(x) || is.integer(x)) || is.object(x) || !is.null(dim(x)) ||
      length(x) < 1L) {
    multiplier_abort(
      sprintf("'%s' must be a non-empty character or integer vector without NA.", name),
      c(align_subclass, "cdid_multiplier_input_error")
    )
  }
  if (!scan_values) return(x)
  if (anyNA(x)) {
    multiplier_abort(sprintf("'%s' contains NA IDs.", name), c(align_subclass, "cdid_multiplier_input_error"))
  }
  if (is.character(x) && !all(nzchar(x))) {
    multiplier_abort(
      sprintf("'%s' must not contain empty character IDs.", name),
      c(align_subclass, "cdid_multiplier_input_error")
    )
  }
  if (anyDuplicated(x) > 0L) {
    multiplier_abort(
      sprintf("'%s' contains duplicate IDs.", name),
      c(dup_subclass, "cdid_multiplier_duplicate_id_error", align_subclass)
    )
  }
  x
}

validate_multiplier_budget <- function(budget) {
  if (!is.list(budget) || is.object(budget)) {
    multiplier_abort(
      "'budget' must be a plain named list of explicit resource and dimension limits.",
      c("cdid_multiplier_budget_error", "cdid_multiplier_input_error")
    )
  }
  required_fields <- c(
    "max_units",
    "max_clusters",
    "max_targets",
    "max_draws",
    "max_if_batch_rows",
    "max_if_batch_cells",
    "max_sign_batch_rows",
    "max_sign_batch_cells",
    "max_workspace_cells",
    "max_workspace_bytes"
  )
  bnames <- names(budget)
  if (is.null(bnames) || any(!nzchar(bnames)) || anyDuplicated(bnames) > 0L ||
      !all(required_fields %in% bnames) || length(setdiff(bnames, c(required_fields, "workers"))) > 0L) {
    multiplier_abort(
      paste0(
        "'budget' must contain exactly the required fields: ",
        paste(required_fields, collapse = ", "),
        " (and optionally the declared worker count `workers`)."
      ),
      c("cdid_multiplier_budget_error", "cdid_multiplier_input_error")
    )
  }
  out <- vector("list", length(required_fields))
  names(out) <- required_fields
  for (nm in required_fields) {
    out[[nm]] <- validate_whole_number_limit(budget[[nm]], nm)
  }
  if ("workers" %in% bnames) {
    w <- validate_whole_number_limit(budget$workers, "workers")
    if (w > .Machine$integer.max) {
      multiplier_abort("Budget field 'workers' exceeds the supported process count.", "cdid_multiplier_budget_error")
    }
    out$workers <- w
  }
  out
}

# Declared multiplier worker count of a validated budget (1 = serial).
multiplier_budget_workers <- function(budget) {
  if (is.null(budget$workers)) 1 else as.double(budget$workers)
}

# Worker-process route of base `parallel` on this platform: forked workers on
# Unix-alikes (as did's mboot), a PSOCK cluster on Windows, where forking is
# unavailable. Both are documented routes of R's own `parallel` package.
chained_parallel_route <- function() {
  if (identical(.Platform$OS.type, "windows")) "psock" else "fork"
}

# Whole-buffer count of one parallel round (CP08): the parent holds w sign
# batches (b x C) and their returned products (b x J, received and stored);
# each worker holds a copy of the centered sums (C x J), its batch and its
# product, each with a serialization buffer.
multiplier_parallel_cells <- function(workers, n_clusters, n_targets, batch_rows) {
  w <- as.double(workers)
  if (w <= 1) return(0)
  b <- as.double(batch_rows); C <- as.double(n_clusters); J <- as.double(n_targets)
  w * (b * C + 2 * b * J) + w * (2 * C * J + 2 * b * C + 2 * b * J + 256)
}

# Products of one round of validated sign batches with the centered cluster
# sums, computed by `workers` processes. The worker function is the serial
# expression with no package state; a failed or malformed worker result
# refuses (no serial fallback, no regenerated signs). `cluster` is a live
# PSOCK cluster for the "psock" route (NULL on the fork route).
chained_multiplier_parallel_products <- function(smats, S_centered, inv_sqrt_C, workers,
                                                 route = chained_parallel_route(), cluster = NULL) {
  f <- function(smat, S, inv) (smat %*% S) * inv
  environment(f) <- baseenv()
  res <- if (identical(route, "fork")) {
    parallel::mclapply(smats, f, S = S_centered, inv = inv_sqrt_C, mc.cores = as.integer(workers),
      mc.preschedule = TRUE, mc.set.seed = FALSE, mc.silent = TRUE)
  } else if (identical(route, "psock") && !is.null(cluster)) {
    # parLapply re-raises worker errors in this process; keep them typed.
    tryCatch(parallel::parLapply(cluster, smats, f, S = S_centered, inv = inv_sqrt_C),
      error = function(e) multiplier_abort(paste0("A PSOCK multiplier worker failed (", conditionMessage(e),
        "); no serial fallback or sign regeneration is attempted."), c("cdid_multiplier_parallel_error", "cdid_multiplier_input_error")))
  } else {
    multiplier_abort("Unknown or unavailable parallel multiplier route.", c("cdid_multiplier_parallel_error", "cdid_multiplier_input_error"))
  }
  if (!is.list(res) || length(res) != length(smats)) {
    multiplier_abort("Parallel multiplier workers did not return one product per sign batch; no serial fallback is attempted.",
      c("cdid_multiplier_parallel_error", "cdid_multiplier_input_error"))
  }
  for (k in seq_along(res)) {
    r <- res[[k]]
    if (!is.matrix(r) || typeof(r) != "double" || !identical(dim(r), c(nrow(smats[[k]]), ncol(S_centered)))) {
      multiplier_abort(sprintf("Parallel multiplier worker %d failed or returned a malformed product (%s); no serial fallback or sign regeneration is attempted.",
        k, if (inherits(r, "try-error")) gsub("\\s+", " ", as.character(r)) else class(r)[1L]),
        c("cdid_multiplier_parallel_error", "cdid_multiplier_input_error"))
    }
  }
  res
}

compute_multiplier_workspace_estimate <- function(n_units,
                                                  n_clusters,
                                                  n_targets,
                                                  n_draws,
                                                  max_if_batch_rows,
                                                  max_sign_batch_rows,
                                                  workers = 1) {
  n_units <- validate_positive_int_scalar(n_units, "n_units")
  n_clusters <- validate_positive_int_scalar(n_clusters, "n_clusters", matrix_dimension = TRUE)
  n_targets <- validate_positive_int_scalar(n_targets, "n_targets", matrix_dimension = TRUE)
  n_draws <- validate_positive_int_scalar(n_draws, "n_draws", matrix_dimension = TRUE)
  max_if_batch_rows <- validate_whole_number_limit(max_if_batch_rows, "max_if_batch_rows")
  max_sign_batch_rows <- validate_whole_number_limit(max_sign_batch_rows, "max_sign_batch_rows")
  nu <- as.double(n_units)
  nc <- as.double(n_clusters)
  nt <- as.double(n_targets)
  nb <- as.double(n_draws)

  cluster_sum_cells <- nc * nt
  centered_cluster_sum_cells <- nc * nt
  # Three owned covariance matrices plus conservative transpose/symmetrization scratch.
  covariance_cells <- 6.0 * nt * nt
  # Internal B x J storage is always budgeted regardless of return_raw_draws.
  raw_draw_cells <- nb * nt
  draw_product_cells <- nb * nt
  unit_bookkeeping_entries <- 6.0 * nu + 4.0 * nc
  quantile_scratch_cells <- 8.0 * nb + 32.0 * nt
  if_scan_cells <- 2.0 * min(nu, max_if_batch_rows, .Machine$integer.max) * nt
  sign_scan_cells <- 4.0 * min(nb, max_sign_batch_rows) * nc
  # Declared multiplier workers (budget$workers): one round of whole worker
  # and parent buffers on top of the serial persistent payload.
  workers <- validate_whole_number_limit(workers, "workers")
  parallel_worker_cells <- multiplier_parallel_cells(workers, nc, nt, min(nb, max_sign_batch_rows))

  total_persistent_cells <- (
    cluster_sum_cells +
    centered_cluster_sum_cells +
    covariance_cells +
    raw_draw_cells +
    draw_product_cells +
    unit_bookkeeping_entries +
    quantile_scratch_cells + if_scan_cells + sign_scan_cells + parallel_worker_cells
  )
  # Conservative numeric/logical/reference slot payload. R object headers,
  # string bytes, allocator/hash/BLAS internals and callback-owned inputs are
  # excluded: this is an explicit entry estimate, never a peak-RSS bound.
  estimated_workspace_bytes <- 8.0 * total_persistent_cells
  if (!is.finite(estimated_workspace_bytes) || estimated_workspace_bytes > 9007199254740991) {
    multiplier_abort("Workspace payload exceeds the exact caller-budget integer range.", "cdid_multiplier_budget_error")
  }

  list(
    cluster_sum_cells = cluster_sum_cells,
    centered_cluster_sum_cells = centered_cluster_sum_cells,
    covariance_cells = covariance_cells,
    raw_draw_cells = raw_draw_cells,
    draw_product_cells = draw_product_cells,
    unit_bookkeeping_entries = unit_bookkeeping_entries,
    quantile_scratch_cells = quantile_scratch_cells,
    if_scan_cells = if_scan_cells,
    sign_scan_cells = sign_scan_cells,
    workers = workers,
    parallel_worker_cells = parallel_worker_cells,
    total_persistent_cells = total_persistent_cells,
    estimated_workspace_bytes = estimated_workspace_bytes
  )
}

type1_order_index <- function(n_draws, prob) {
  idx <- as.integer(ceiling(as.double(n_draws) * as.double(prob)))
  if (idx < 1L) {
    idx <- 1L
  } else if (idx > n_draws) {
    idx <- as.integer(n_draws)
  }
  idx
}

chained_cluster_multiplier <- function(target_ids,
                                       expected_unit_ids,
                                       unit_cluster_map,
                                       cluster_ids,
                                       n_denom,
                                       center,
                                       if_batch_fn,
                                       n_if_batches,
                                       sign_batch_fn,
                                       n_sign_batches,
                                       expected_draw_ids,
                                       n_draws,
                                       alpha,
                                       return_raw_draws,
                                       budget) {
  # 1. Explicit required-argument check (no defaults permitted).
  if (missing(target_ids) || missing(expected_unit_ids) ||
      missing(unit_cluster_map) || missing(cluster_ids) ||
      missing(n_denom) || missing(center) ||
      missing(if_batch_fn) || missing(n_if_batches) ||
      missing(sign_batch_fn) || missing(n_sign_batches) ||
      missing(expected_draw_ids) || missing(n_draws) ||
      missing(alpha) || missing(return_raw_draws) ||
      missing(budget)) {
    multiplier_abort(
      "All 15 arguments to 'chained_cluster_multiplier' must be explicitly supplied without defaults.",
      c("cdid_multiplier_missing_arg_error", "cdid_multiplier_input_error")
    )
  }

  # 2. Validate budget and scalar/function arguments before any workspace allocation.
  checked_budget <- validate_multiplier_budget(budget)

  if (!is.numeric(n_denom) || is.object(n_denom) || !is.null(dim(n_denom)) || length(n_denom) != 1L || is.na(n_denom) ||
      !is.finite(n_denom) || n_denom <= 0) {
    multiplier_abort(
      "'n_denom' must be a positive finite numeric scalar.",
      "cdid_multiplier_input_error"
    )
  }
  n_denom_num <- as.double(n_denom)

  if (!is.logical(center) || is.object(center) || !is.null(dim(center)) || length(center) != 1L || is.na(center)) {
    multiplier_abort(
      "'center' must be an explicit non-NA logical scalar (TRUE or FALSE).",
      "cdid_multiplier_input_error"
    )
  }
  if (!is.logical(return_raw_draws) || is.object(return_raw_draws) || !is.null(dim(return_raw_draws)) || length(return_raw_draws) != 1L ||
      is.na(return_raw_draws)) {
    multiplier_abort(
      "'return_raw_draws' must be an explicit non-NA logical scalar (TRUE or FALSE).",
      "cdid_multiplier_input_error"
    )
  }
  if (!is.numeric(alpha) || is.object(alpha) || !is.null(dim(alpha)) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    multiplier_abort(
      "'alpha' must be a finite numeric scalar strictly in (0, 1).",
      "cdid_multiplier_input_error"
    )
  }
  alpha_num <- as.double(alpha)

  n_if_batches_int <- validate_positive_int_scalar(n_if_batches, "n_if_batches")
  n_sign_batches_int <- validate_positive_int_scalar(n_sign_batches, "n_sign_batches")
  n_draws_int <- validate_positive_int_scalar(
    n_draws,
    "n_draws",
    subclass = c("cdid_multiplier_draw_alignment_error", "cdid_multiplier_input_error"),
    matrix_dimension = TRUE
  )

  if (!is.function(if_batch_fn)) {
    multiplier_abort("'if_batch_fn' must be a function.", "cdid_multiplier_input_error")
  }
  if (!is.function(sign_batch_fn)) {
    multiplier_abort("'sign_batch_fn' must be a function.", "cdid_multiplier_input_error")
  }

  # 3. Validate declared IDs and preflight dimensions/workspace against budget.
  # Shape/length checks precede content scans and duplicate/hash bookkeeping.
  for (spec in list(
    list(target_ids, "target_ids", checked_budget$max_targets, TRUE),
    list(expected_unit_ids, "expected_unit_ids", checked_budget$max_units, FALSE),
    list(cluster_ids, "cluster_ids", checked_budget$max_clusters, TRUE),
    list(expected_draw_ids, "expected_draw_ids", checked_budget$max_draws, TRUE))) {
    validate_id_vector(spec[[1L]], spec[[2L]], "cdid_multiplier_duplicate_id_error",
                       "cdid_multiplier_input_error", scan_values = FALSE)
    if (length(spec[[1L]]) > spec[[3L]] ||
        (spec[[4L]] && length(spec[[1L]]) > .Machine$integer.max)) {
      multiplier_abort(sprintf("Declared '%s' exceeds its caller bound or genuine matrix dimension limit.", spec[[2L]]),
                       "cdid_multiplier_budget_error")
    }
  }
  n_units <- length(expected_unit_ids)
  n_clusters <- length(cluster_ids)
  n_targets <- length(target_ids)

  if (n_if_batches_int > n_units || n_sign_batches_int > n_draws_int) {
    multiplier_abort("Nonempty batch counts cannot exceed the declared unit/draw counts.",
                     "cdid_multiplier_input_error")
  }

  if (length(expected_draw_ids) != n_draws_int) {
    multiplier_abort(
      "Length of 'expected_draw_ids' must equal 'n_draws'.",
      "cdid_multiplier_draw_alignment_error"
    )
  }

  if (as.double(n_units) > checked_budget$max_units) {
    multiplier_abort(
      sprintf("Number of expected units (%d) exceeds budget$max_units (%.0f).", n_units, checked_budget$max_units),
      "cdid_multiplier_budget_error"
    )
  }
  if (as.double(n_clusters) > checked_budget$max_clusters) {
    multiplier_abort(
      sprintf("Number of declared clusters (%d) exceeds budget$max_clusters (%.0f).", n_clusters, checked_budget$max_clusters),
      "cdid_multiplier_budget_error"
    )
  }
  if (as.double(n_targets) > checked_budget$max_targets) {
    multiplier_abort(
      sprintf("Number of fixed targets (%d) exceeds budget$max_targets (%.0f).", n_targets, checked_budget$max_targets),
      "cdid_multiplier_budget_error"
    )
  }
  if (as.double(n_draws_int) > checked_budget$max_draws) {
    multiplier_abort(
      sprintf("Number of draws B (%d) exceeds budget$max_draws (%.0f).", n_draws_int, checked_budget$max_draws),
      "cdid_multiplier_budget_error"
    )
  }

  workspace_estimate <- compute_multiplier_workspace_estimate(
    n_units = n_units,
    n_clusters = n_clusters,
    n_targets = n_targets,
    n_draws = n_draws_int,
    max_if_batch_rows = checked_budget$max_if_batch_rows,
    max_sign_batch_rows = checked_budget$max_sign_batch_rows,
    workers = multiplier_budget_workers(checked_budget)
  )
  if (workspace_estimate$total_persistent_cells > checked_budget$max_workspace_cells) {
    multiplier_abort(
      sprintf(
        "Estimated persistent workspace cells (%.0f) exceeds budget$max_workspace_cells (%.0f).",
        workspace_estimate$total_persistent_cells,
        checked_budget$max_workspace_cells
      ),
      "cdid_multiplier_budget_error"
    )
  }
  if (workspace_estimate$estimated_workspace_bytes > checked_budget$max_workspace_bytes) {
    multiplier_abort(
      sprintf(
        "Estimated workspace bytes (%.0f) exceeds budget$max_workspace_bytes (%.0f).",
        workspace_estimate$estimated_workspace_bytes,
        checked_budget$max_workspace_bytes
      ),
      "cdid_multiplier_budget_error"
    )
  }

  validate_id_vector(target_ids, "target_ids", "cdid_multiplier_duplicate_target_error",
                     "cdid_multiplier_target_alignment_error")
  validate_id_vector(expected_unit_ids, "expected_unit_ids", "cdid_multiplier_duplicate_unit_error",
                     "cdid_multiplier_unit_alignment_error")
  validate_id_vector(cluster_ids, "cluster_ids", "cdid_multiplier_duplicate_cluster_error",
                     "cdid_multiplier_cluster_alignment_error")
  validate_id_vector(expected_draw_ids, "expected_draw_ids", "cdid_multiplier_duplicate_draw_error",
                     "cdid_multiplier_draw_alignment_error")
  # 4. Validate unit_cluster_map and build O(n) identity lookup before callbacks.
  if (!is.list(unit_cluster_map) || (is.object(unit_cluster_map) && !identical(class(unit_cluster_map), "data.frame")) ||
      anyDuplicated(names(unit_cluster_map)) > 0L) {
    multiplier_abort(
      "'unit_cluster_map' must be a list or data.frame with 'unit_ids' and 'cluster_ids'.",
      c("cdid_multiplier_unit_alignment_error", "cdid_multiplier_input_error")
    )
  }
  forbidden_weight_fields <- c("weights", "unit_weights", "cluster_weights", "sampling_weights")
  if (any(forbidden_weight_fields %in% names(unit_cluster_map))) {
    multiplier_abort(
      "'unit_cluster_map' must not supply weights; complete IF rows are already weighted.",
      "cdid_multiplier_double_weight_error"
    )
  }
  if (!all(c("unit_ids", "cluster_ids") %in% names(unit_cluster_map))) {
    multiplier_abort(
      "'unit_cluster_map' must contain 'unit_ids' and 'cluster_ids'.",
      c("cdid_multiplier_unit_alignment_error", "cdid_multiplier_input_error")
    )
  }

  map_unit_ids <- unit_cluster_map$unit_ids
  map_cluster_ids <- unit_cluster_map$cluster_ids

  if (is.object(map_unit_ids) || !is.null(dim(map_unit_ids))) {
    multiplier_abort("Unit-map IDs must be base vectors without class or dimensions.", "cdid_multiplier_unit_alignment_error")
  }
  if (is.object(map_cluster_ids) || !is.null(dim(map_cluster_ids))) {
    multiplier_abort("Cluster-map IDs must be base vectors without class or dimensions.", "cdid_multiplier_cluster_alignment_error")
  }

  if (length(map_unit_ids) != n_units || length(map_cluster_ids) != n_units) {
    multiplier_abort("The unit-cluster map must have exactly one entry per expected unit.",
      c("cdid_multiplier_missing_unit_error", "cdid_multiplier_unit_alignment_error"))
  }

  if (typeof(map_unit_ids) != typeof(expected_unit_ids) || is.object(map_unit_ids) ||
      !is.null(dim(map_unit_ids)) || anyNA(map_unit_ids)) {
    multiplier_abort(
      "'unit_cluster_map$unit_ids' must have the same storage type as 'expected_unit_ids' without NA.",
      "cdid_multiplier_unit_alignment_error"
    )
  }
  if (is.character(map_unit_ids) && !all(nzchar(map_unit_ids))) {
    multiplier_abort(
      "'unit_cluster_map$unit_ids' must not contain empty character IDs.",
      "cdid_multiplier_unit_alignment_error"
    )
  }
  if (typeof(map_cluster_ids) != typeof(cluster_ids) || is.object(map_cluster_ids) ||
      !is.null(dim(map_cluster_ids)) || anyNA(map_cluster_ids)) {
    multiplier_abort(
      "'unit_cluster_map$cluster_ids' must have the same storage type as 'cluster_ids' without NA.",
      "cdid_multiplier_cluster_alignment_error"
    )
  }
  if (is.character(map_cluster_ids) && !all(nzchar(map_cluster_ids))) {
    multiplier_abort(
      "'unit_cluster_map$cluster_ids' must not contain empty character IDs.",
      "cdid_multiplier_cluster_alignment_error"
    )
  }
  if (length(map_unit_ids) != length(map_cluster_ids)) {
    multiplier_abort(
      "'unit_cluster_map$unit_ids' and 'unit_cluster_map$cluster_ids' must have equal length.",
      "cdid_multiplier_unit_alignment_error"
    )
  }
  if (anyDuplicated(map_unit_ids) > 0L) {
    multiplier_abort(
      "'unit_cluster_map$unit_ids' contains duplicate unit IDs.",
      c("cdid_multiplier_duplicate_unit_error", "cdid_multiplier_duplicate_id_error", "cdid_multiplier_unit_alignment_error")
    )
  }
  if (anyNA(match(map_unit_ids, expected_unit_ids))) {
    multiplier_abort(
      "'unit_cluster_map$unit_ids' contains unknown unit IDs not in 'expected_unit_ids'.",
      c("cdid_multiplier_unknown_unit_error", "cdid_multiplier_unknown_id_error", "cdid_multiplier_unit_alignment_error")
    )
  }
  if (length(map_unit_ids) != n_units || anyNA(match(expected_unit_ids, map_unit_ids))) {
    multiplier_abort(
      "'unit_cluster_map$unit_ids' is missing one or more expected units from 'expected_unit_ids'.",
      c("cdid_multiplier_missing_unit_error", "cdid_multiplier_unit_alignment_error")
    )
  }

  unit_to_map_pos <- match(expected_unit_ids, map_unit_ids)
  ordered_unit_clusters <- map_cluster_ids[unit_to_map_pos]
  unit_cluster_idx <- match(ordered_unit_clusters, cluster_ids)
  if (anyNA(unit_cluster_idx)) {
    multiplier_abort(
      "'unit_cluster_map$cluster_ids' maps one or more units to an undeclared cluster ID.",
      c("cdid_multiplier_unknown_cluster_error", "cdid_multiplier_unknown_id_error", "cdid_multiplier_cluster_alignment_error")
    )
  }

  declared_cluster_unit_counts <- tabulate(unit_cluster_idx, nbins = n_clusters)
  names(declared_cluster_unit_counts) <- as.character(cluster_ids)

  # 5. Stream complete target IF batches into C x J cluster sums S_mat.
  target_names_chr <- as.character(target_ids)
  cluster_names_chr <- as.character(cluster_ids)
  draw_names_chr <- as.character(expected_draw_ids)

  S_mat <- matrix(
    0.0,
    nrow = n_clusters,
    ncol = n_targets,
    dimnames = list(cluster_names_chr, target_names_chr)
  )
  seen_unit <- logical(n_units)
  streamed_cluster_unit_counts <- if (n_units <= .Machine$integer.max) integer(n_clusters) else numeric(n_clusters)
  total_if_rows_seen <- 0.0

  for (b_idx in seq_len(n_if_batches_int)) {
    Sys.sleep(0)
    batch <- if_batch_fn(b_idx)

    if (!is.list(batch) || is.object(batch) || anyDuplicated(names(batch)) > 0L) {
      multiplier_abort(
        sprintf("IF batch %.0f must be a plain list.", b_idx),
        c("cdid_multiplier_if_dimension_error", "cdid_multiplier_input_error")
      )
    }

    if (any(forbidden_weight_fields %in% names(batch))) {
      multiplier_abort(
        sprintf("IF batch %.0f supplies weights; complete IF rows must not be weighted again.", b_idx),
        "cdid_multiplier_double_weight_error"
      )
    }
    if (!all(c("unit_ids", "target_ids", "if_matrix") %in% names(batch))) {
      multiplier_abort(
        sprintf("IF batch %.0f must contain 'unit_ids', 'target_ids', and 'if_matrix'.", b_idx),
        c("cdid_multiplier_if_dimension_error", "cdid_multiplier_input_error")
      )
    }

    if_mat <- batch$if_matrix
    if (!is.matrix(if_mat) || typeof(if_mat) != "double" || is.object(if_mat)) {
      multiplier_abort(
        sprintf("IF batch %.0f 'if_matrix' must be a base double matrix without coercion.", b_idx),
        c("cdid_multiplier_if_dimension_error", "cdid_multiplier_input_error")
      )
    }

    r_b <- nrow(if_mat)
    c_b <- ncol(if_mat)
    if (r_b < 1L || c_b != n_targets) {
      multiplier_abort(
        sprintf(
          "IF batch %.0f 'if_matrix' has dimensions %d x %d; expected >=1 x %d.",
          b_idx, r_b, c_b, n_targets
        ),
        c("cdid_multiplier_if_dimension_error", "cdid_multiplier_target_alignment_error")
      )
    }
    if (as.double(r_b) > checked_budget$max_if_batch_rows ||
        as.double(r_b) * as.double(c_b) > checked_budget$max_if_batch_cells) {
      multiplier_abort(
        sprintf(
          "IF batch %.0f dimensions (%d x %d) exceed caller IF batch budget.",
          b_idx, r_b, c_b
        ),
        "cdid_multiplier_budget_error"
      )
    }

    if (r_b > n_units - total_if_rows_seen) {
      multiplier_abort("IF batch exceeds the remaining expected unit count.", "cdid_multiplier_unit_alignment_error")
    }
    # Validate batch target_ids alignment before scanning matrix values.
    b_tids <- batch$target_ids
    if (typeof(b_tids) != typeof(target_ids) || is.object(b_tids) ||
        !is.null(dim(b_tids)) || anyNA(b_tids)) {
      multiplier_abort(
        sprintf("IF batch %.0f 'target_ids' must match storage type of 'target_ids' without NA.", b_idx),
        "cdid_multiplier_target_alignment_error"
      )
    }
    if (anyDuplicated(b_tids) > 0L) {
      multiplier_abort(
        sprintf("IF batch %.0f 'target_ids' contains duplicate target IDs.", b_idx),
        c("cdid_multiplier_duplicate_target_error", "cdid_multiplier_duplicate_id_error", "cdid_multiplier_target_alignment_error")
      )
    }
    if (length(b_tids) != n_targets) {
      multiplier_abort(
        sprintf("IF batch %.0f 'target_ids' length (%d) does not match fixed target count (%d).", b_idx, length(b_tids), n_targets),
        "cdid_multiplier_target_alignment_error"
      )
    }
    if (anyNA(match(b_tids, target_ids))) {
      multiplier_abort(
        sprintf("IF batch %.0f 'target_ids' contains unknown target ID.", b_idx),
        c("cdid_multiplier_unknown_target_error", "cdid_multiplier_unknown_id_error", "cdid_multiplier_target_alignment_error")
      )
    }
    if (!identical(b_tids, target_ids)) {
      multiplier_abort(
        sprintf("IF batch %.0f 'target_ids' order does not match fixed 'target_ids'.", b_idx),
        c("cdid_multiplier_target_order_error", "cdid_multiplier_order_error", "cdid_multiplier_target_alignment_error")
      )
    }
    if (!is.null(colnames(if_mat)) && !identical(colnames(if_mat), target_names_chr)) {
      multiplier_abort(
        sprintf("IF batch %.0f colnames(if_matrix) do not match fixed 'target_ids'.", b_idx),
        c("cdid_multiplier_target_order_error", "cdid_multiplier_order_error", "cdid_multiplier_target_alignment_error")
      )
    }

    # Validate batch unit_ids alignment before scanning matrix values.
    b_uids <- batch$unit_ids
    if (typeof(b_uids) != typeof(expected_unit_ids) || is.object(b_uids) ||
        !is.null(dim(b_uids)) || length(b_uids) != r_b || anyNA(b_uids)) {
      multiplier_abort(
        sprintf("IF batch %.0f 'unit_ids' must have length nrow(if_matrix) and match storage type of 'expected_unit_ids' without NA.", b_idx),
        "cdid_multiplier_unit_alignment_error"
      )
    }
    if (is.character(b_uids) && !all(nzchar(b_uids))) {
      multiplier_abort(
        sprintf("IF batch %.0f 'unit_ids' contains empty string ID.", b_idx),
        "cdid_multiplier_unit_alignment_error"
      )
    }
    if (!is.null(rownames(if_mat)) && !identical(rownames(if_mat), as.character(b_uids))) {
      multiplier_abort(
        sprintf("IF batch %.0f rownames(if_matrix) do not match batch 'unit_ids'.", b_idx),
        "cdid_multiplier_unit_alignment_error"
      )
    }
    if (anyDuplicated(b_uids) > 0L) {
      multiplier_abort(
        sprintf("IF batch %.0f contains duplicate unit IDs.", b_idx),
        c("cdid_multiplier_duplicate_unit_error", "cdid_multiplier_duplicate_id_error", "cdid_multiplier_unit_alignment_error")
      )
    }
    u_pos <- match(b_uids, expected_unit_ids)
    if (anyNA(u_pos)) {
      multiplier_abort(
        sprintf("IF batch %.0f contains unknown unit ID not in 'expected_unit_ids'.", b_idx),
        c("cdid_multiplier_unknown_unit_error", "cdid_multiplier_unknown_id_error", "cdid_multiplier_unit_alignment_error")
      )
    }
    if (any(seen_unit[u_pos])) {
      multiplier_abort(
        sprintf("IF batch %.0f repeats a unit ID already seen in an earlier batch.", b_idx),
        c("cdid_multiplier_duplicate_unit_error", "cdid_multiplier_duplicate_id_error", "cdid_multiplier_unit_alignment_error")
      )
    }
    if (total_if_rows_seen + r_b > n_units) {
      multiplier_abort(
        sprintf("Total IF rows across batches exceeds expected unit count (%d).", n_units),
        "cdid_multiplier_unit_alignment_error"
      )
    }

    if ("cluster_ids" %in% names(batch) && !is.null(batch$cluster_ids)) {
      b_cids <- batch$cluster_ids
      if (typeof(b_cids) != typeof(cluster_ids) || is.object(b_cids) ||
          !is.null(dim(b_cids)) || length(b_cids) != r_b || anyNA(b_cids)) {
        multiplier_abort(
          sprintf("IF batch %.0f 'cluster_ids' must have length nrow(if_matrix) and match storage type of 'cluster_ids'.", b_idx),
          "cdid_multiplier_cluster_alignment_error"
        )
      }
      c_pos <- match(b_cids, cluster_ids)
      if (anyNA(c_pos)) {
        multiplier_abort(
          sprintf("IF batch %.0f 'cluster_ids' contains unknown cluster ID.", b_idx),
          c("cdid_multiplier_unknown_cluster_error", "cdid_multiplier_unknown_id_error", "cdid_multiplier_cluster_alignment_error")
        )
      }
      if (!identical(c_pos, unit_cluster_idx[u_pos])) {
        multiplier_abort(
          sprintf("IF batch %.0f 'cluster_ids' disagrees with declared 'unit_cluster_map'.", b_idx),
          "cdid_multiplier_cluster_alignment_error"
        )
      }
    }

    # Scan matrix values for finiteness after all dimension/budget/identity checks.
    if (!all(is.finite(if_mat))) {
      multiplier_abort(
        sprintf("IF batch %.0f 'if_matrix' contains non-finite values (NA, NaN, or Inf).", b_idx),
        c("cdid_multiplier_nonfinite_if_error", "cdid_multiplier_input_error")
      )
    }

    seen_unit[u_pos] <- TRUE
    total_if_rows_seen <- total_if_rows_seen + r_b
    for (r in seq_len(r_b)) {
      c_idx <- unit_cluster_idx[u_pos[r]]
      S_mat[c_idx, ] <- S_mat[c_idx, ] + if_mat[r, ]
      streamed_cluster_unit_counts[c_idx] <- streamed_cluster_unit_counts[c_idx] + 1L
    }
  }

  if (total_if_rows_seen != n_units || !all(seen_unit) ||
      !identical(as.double(streamed_cluster_unit_counts), as.double(declared_cluster_unit_counts))) {
    multiplier_abort(
      sprintf(
        "Incomplete IF stream: observed %.0f of %.0f expected units across %.0f batches.",
        total_if_rows_seen, n_units, n_if_batches_int
      ),
      c("cdid_multiplier_missing_unit_error", "cdid_multiplier_unit_alignment_error")
    )
  }
  names(streamed_cluster_unit_counts) <- cluster_names_chr

  # 6. Center across declared cluster set (if requested) and compute candidate covariance.
  cluster_means <- colMeans(S_mat)
  names(cluster_means) <- target_names_chr

  if (isTRUE(center)) {
    S_centered <- sweep(S_mat, 2L, cluster_means, FUN = "-")
  } else {
    S_centered <- S_mat
  }
  dimnames(S_centered) <- list(cluster_names_chr, target_names_chr)

  cluster_crossprod <- crossprod(S_centered)
  cluster_crossprod <- 0.5 * (cluster_crossprod + t(cluster_crossprod))
  dimnames(cluster_crossprod) <- list(target_names_chr, target_names_chr)

  cov_mat <- cluster_crossprod / (n_denom_num * n_denom_num)
  omega_cl_hat <- cluster_crossprod / n_denom_num
  dimnames(cov_mat) <- list(target_names_chr, target_names_chr)
  dimnames(omega_cl_hat) <- list(target_names_chr, target_names_chr)

  analytic_var <- diag(cov_mat)
  names(analytic_var) <- target_names_chr

  # 7. Stream supplied +/-1 sign matrices into bounded B x J raw draws bres.
  bres <- matrix(
    NA_real_,
    nrow = n_draws_int,
    ncol = n_targets,
    dimnames = list(draw_names_chr, target_names_chr)
  )
  seen_draw <- logical(n_draws_int)
  draw_cursor <- 0L
  inv_sqrt_C <- 1.0 / sqrt(as.double(n_clusters))
  # Declared workers: batches are still requested and validated here, one at a
  # time in draw order; only their products with S_centered run in rounds of
  # at most `workers` batches on worker processes (whole buffers preflighted).
  n_workers <- multiplier_budget_workers(checked_budget)
  parallel_route <- if (n_workers > 1) chained_parallel_route() else "serial"
  psock_cluster <- NULL
  if (identical(parallel_route, "psock")) {
    psock_cluster <- parallel::makePSOCKcluster(as.integer(n_workers))
    on.exit(parallel::stopCluster(psock_cluster), add = TRUE)
  }
  pending <- list(); n_rounds <- 0L

  for (s_idx in seq_len(n_sign_batches_int)) {
    Sys.sleep(0)
    sbatch <- sign_batch_fn(s_idx)

    if (!is.list(sbatch) || is.object(sbatch) || anyDuplicated(names(sbatch)) > 0L) {
      multiplier_abort(
        sprintf("Sign batch %.0f must be a plain list.", s_idx),
        c("cdid_multiplier_sign_dimension_error", "cdid_multiplier_input_error")
      )
    }

    if (!all(c("draw_ids", "cluster_ids", "sign_matrix") %in% names(sbatch))) {
      multiplier_abort(
        sprintf("Sign batch %.0f must contain 'draw_ids', 'cluster_ids', and 'sign_matrix'.", s_idx),
        c("cdid_multiplier_sign_dimension_error", "cdid_multiplier_input_error")
      )
    }

    smat <- sbatch$sign_matrix
    if (!is.matrix(smat) || !(typeof(smat) == "double" || typeof(smat) == "integer") ||
        is.object(smat)) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'sign_matrix' must be a base double or integer matrix without coercion.", s_idx),
        c("cdid_multiplier_sign_dimension_error", "cdid_multiplier_input_error")
      )
    }

    b_batch <- nrow(smat)
    c_batch <- ncol(smat)
    if (b_batch < 1L) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'sign_matrix' must have at least 1 row.", s_idx),
        "cdid_multiplier_sign_dimension_error"
      )
    }
    if (as.double(b_batch) > checked_budget$max_sign_batch_rows ||
        as.double(b_batch) * as.double(c_batch) > checked_budget$max_sign_batch_cells) {
      multiplier_abort(
        sprintf(
          "Sign batch %.0f dimensions (%d x %d) exceed caller sign batch budget.",
          s_idx, b_batch, c_batch
        ),
        "cdid_multiplier_budget_error"
      )
    }

    if (b_batch > n_draws_int - draw_cursor) {
      multiplier_abort("Sign batch exceeds the remaining expected draw count.", "cdid_multiplier_incomplete_draws_error")
    }
    # Validate cluster column alignment before scanning sign_matrix content.
    sb_cids <- sbatch$cluster_ids
    if (typeof(sb_cids) != typeof(cluster_ids) || is.object(sb_cids) ||
        !is.null(dim(sb_cids)) || anyNA(sb_cids)) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'cluster_ids' must match storage type of 'cluster_ids' without NA.", s_idx),
        "cdid_multiplier_cluster_alignment_error"
      )
    }
    if (is.character(sb_cids) && !all(nzchar(sb_cids))) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'cluster_ids' contains empty string ID.", s_idx),
        "cdid_multiplier_cluster_alignment_error"
      )
    }
    if (anyDuplicated(sb_cids) > 0L) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'cluster_ids' contains duplicate cluster IDs.", s_idx),
        c("cdid_multiplier_duplicate_cluster_error", "cdid_multiplier_duplicate_id_error", "cdid_multiplier_cluster_alignment_error")
      )
    }
    if (anyNA(match(sb_cids, cluster_ids))) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'cluster_ids' contains unknown cluster ID.", s_idx),
        c("cdid_multiplier_unknown_cluster_error", "cdid_multiplier_unknown_id_error", "cdid_multiplier_cluster_alignment_error")
      )
    }
    if (c_batch != n_clusters || length(sb_cids) != n_clusters ||
        anyNA(match(cluster_ids, sb_cids))) {
      multiplier_abort(
        sprintf(
          "Sign batch %.0f has incomplete/mismatched cluster columns (ncol=%d, length(cluster_ids)=%d, expected C=%d).",
          s_idx, c_batch, length(sb_cids), n_clusters
        ),
        c("cdid_multiplier_incomplete_cluster_error", "cdid_multiplier_cluster_alignment_error", "cdid_multiplier_sign_dimension_error")
      )
    }
    if (!identical(sb_cids, cluster_ids)) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'cluster_ids' order does not match declared 'cluster_ids'.", s_idx),
        c("cdid_multiplier_cluster_order_error", "cdid_multiplier_order_error", "cdid_multiplier_cluster_alignment_error")
      )
    }
    if (!is.null(colnames(smat)) && !identical(colnames(smat), cluster_names_chr)) {
      multiplier_abort(
        sprintf("Sign batch %.0f colnames(sign_matrix) do not match declared 'cluster_ids'.", s_idx),
        c("cdid_multiplier_cluster_order_error", "cdid_multiplier_order_error", "cdid_multiplier_cluster_alignment_error")
      )
    }

    # Validate draw row alignment before scanning sign_matrix content.
    sb_dids <- sbatch$draw_ids
    if (typeof(sb_dids) != typeof(expected_draw_ids) || is.object(sb_dids) ||
        !is.null(dim(sb_dids)) || length(sb_dids) != b_batch || anyNA(sb_dids)) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'draw_ids' must have length nrow(sign_matrix) and match storage type of 'expected_draw_ids' without NA.", s_idx),
        "cdid_multiplier_draw_alignment_error"
      )
    }
    if (is.character(sb_dids) && !all(nzchar(sb_dids))) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'draw_ids' contains empty string ID.", s_idx),
        "cdid_multiplier_draw_alignment_error"
      )
    }
    if (!is.null(rownames(smat)) && !identical(rownames(smat), as.character(sb_dids))) {
      multiplier_abort(
        sprintf("Sign batch %.0f rownames(sign_matrix) do not match batch 'draw_ids'.", s_idx),
        "cdid_multiplier_draw_alignment_error"
      )
    }
    if (anyDuplicated(sb_dids) > 0L) {
      multiplier_abort(
        sprintf("Sign batch %.0f contains duplicate draw IDs.", s_idx),
        c("cdid_multiplier_duplicate_draw_error", "cdid_multiplier_duplicate_id_error", "cdid_multiplier_draw_alignment_error")
      )
    }
    d_pos <- match(sb_dids, expected_draw_ids)
    if (anyNA(d_pos)) {
      multiplier_abort(
        sprintf("Sign batch %.0f contains unknown draw ID not in 'expected_draw_ids'.", s_idx),
        c("cdid_multiplier_unknown_draw_error", "cdid_multiplier_unknown_id_error", "cdid_multiplier_draw_alignment_error")
      )
    }
    if (any(seen_draw[d_pos])) {
      multiplier_abort(
        sprintf("Sign batch %.0f repeats a draw ID already seen in an earlier sign batch.", s_idx),
        c("cdid_multiplier_duplicate_draw_error", "cdid_multiplier_duplicate_id_error", "cdid_multiplier_draw_alignment_error")
      )
    }
    if (draw_cursor + b_batch > n_draws_int) {
      multiplier_abort(
        sprintf("Total draws across sign batches exceeds declared n_draws (%d).", n_draws_int),
        c("cdid_multiplier_incomplete_draws_error", "cdid_multiplier_draw_alignment_error")
      )
    }
    expected_slice <- seq.int(draw_cursor + 1L, draw_cursor + b_batch)
    if (!identical(d_pos, expected_slice)) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'draw_ids' order does not match expected_draw_ids slice.", s_idx),
        c("cdid_multiplier_draw_order_error", "cdid_multiplier_order_error", "cdid_multiplier_draw_alignment_error")
      )
    }

    # Scan sign_matrix content for exact {-1, +1} entries after alignment/budget checks.
    if (anyNA(smat) || !all(smat == 1 | smat == -1)) {
      multiplier_abort(
        sprintf("Sign batch %.0f 'sign_matrix' entries must be finite and strictly in {-1, +1}.", s_idx),
        c("cdid_multiplier_invalid_sign_error", "cdid_multiplier_input_error")
      )
    }

    seen_draw[expected_slice] <- TRUE
    if (n_workers <= 1) {
      bres[expected_slice, ] <- (smat %*% S_centered) * inv_sqrt_C
    } else {
      pending[[length(pending) + 1L]] <- list(slice = expected_slice, smat = smat)
      if (length(pending) >= n_workers || s_idx == n_sign_batches_int) {
        prods <- chained_multiplier_parallel_products(lapply(pending, `[[`, "smat"), S_centered, inv_sqrt_C,
          n_workers, parallel_route, psock_cluster)
        for (k in seq_along(pending)) bres[pending[[k]]$slice, ] <- prods[[k]]
        pending <- list(); n_rounds <- n_rounds + 1L
      }
    }
    draw_cursor <- draw_cursor + b_batch
  }

  if (draw_cursor != n_draws_int || !all(seen_draw)) {
    multiplier_abort(
      sprintf(
        "Incomplete sign draw stream: observed %d of %d expected draws across %d batches.",
        draw_cursor, n_draws_int, n_sign_batches_int
      ),
      c("cdid_multiplier_incomplete_draws_error", "cdid_multiplier_draw_alignment_error")
    )
  }

  # 8. Compute explicit TYPE 1 empirical order statistics, scales, statuses, and critical value.
  iqr_norm <- stats::qnorm(0.75) - stats::qnorm(0.25)
  se_scale_factor <- sqrt(as.double(n_clusters)) / n_denom_num
  k_q25 <- type1_order_index(n_draws_int, 0.25)
  k_q75 <- type1_order_index(n_draws_int, 0.75)
  k_crit <- type1_order_index(n_draws_int, 1.0 - alpha_num)

  q25_vec <- rep(NA_real_, n_targets)
  q75_vec <- rep(NA_real_, n_targets)
  raw_bootstrap_iqr <- rep(NA_real_, n_targets)
  raw_bootstrap_scale <- rep(NA_real_, n_targets)
  candidate_bootstrap_se <- rep(NA_real_, n_targets)
  analytic_status <- character(n_targets)
  bootstrap_status <- character(n_targets)
  target_status <- character(n_targets)
  analytic_ok <- logical(n_targets)
  bootstrap_ok <- logical(n_targets)
  target_accepted <- logical(n_targets)

  names(q25_vec) <- target_names_chr
  names(q75_vec) <- target_names_chr
  names(raw_bootstrap_iqr) <- target_names_chr
  names(raw_bootstrap_scale) <- target_names_chr
  names(candidate_bootstrap_se) <- target_names_chr
  names(analytic_status) <- target_names_chr
  names(bootstrap_status) <- target_names_chr
  names(target_status) <- target_names_chr
  names(analytic_ok) <- target_names_chr
  names(bootstrap_ok) <- target_names_chr
  names(target_accepted) <- target_names_chr

  cov_matrix_finite <- all(is.finite(cov_mat))

  for (j in seq_len(n_targets)) {
    v_j <- analytic_var[j]
    if (!is.finite(v_j) || !cov_matrix_finite) {
      analytic_ok[j] <- FALSE
      analytic_status[j] <- "nonfinite_analytical_variance"
    } else if (v_j <= 0.0) {
      analytic_ok[j] <- FALSE
      analytic_status[j] <- "zero_analytical_variance"
    } else {
      analytic_ok[j] <- TRUE
      analytic_status[j] <- "ok"
    }

    col_j <- bres[, j]
    if (!all(is.finite(col_j))) {
      bootstrap_ok[j] <- FALSE
      bootstrap_status[j] <- "nonfinite_bootstrap_scale"
    } else {
      sorted_j <- sort.int(col_j, method = "quick")
      q25_j <- sorted_j[k_q25]
      q75_j <- sorted_j[k_q75]
      iqr_j <- q75_j - q25_j
      scale_j <- iqr_j / iqr_norm
      se_j <- scale_j * se_scale_factor

      q25_vec[j] <- q25_j
      q75_vec[j] <- q75_j
      raw_bootstrap_iqr[j] <- iqr_j
      raw_bootstrap_scale[j] <- scale_j

      if (!is.finite(iqr_j) || !is.finite(scale_j) || !is.finite(se_j)) {
        bootstrap_ok[j] <- FALSE
        bootstrap_status[j] <- "nonfinite_bootstrap_scale"
      } else if (iqr_j <= 0.0 || scale_j <= 0.0 || se_j <= 0.0) {
        bootstrap_ok[j] <- FALSE
        bootstrap_status[j] <- "zero_bootstrap_scale"
        candidate_bootstrap_se[j] <- se_j
      } else {
        bootstrap_ok[j] <- TRUE
        bootstrap_status[j] <- "ok"
        candidate_bootstrap_se[j] <- se_j
      }
    }

    target_accepted[j] <- analytic_ok[j] && bootstrap_ok[j]
    if (target_accepted[j]) {
      target_status[j] <- "ok"
    } else if (!analytic_ok[j] && !bootstrap_ok[j]) {
      if (identical(analytic_status[[j]], "zero_analytical_variance") &&
          identical(bootstrap_status[[j]], "zero_bootstrap_scale")) {
        target_status[j] <- "zero_analytical_and_bootstrap_scale"
      } else {
        target_status[j] <- paste0(analytic_status[[j]], "__", bootstrap_status[[j]])
      }
    } else if (!analytic_ok[j]) {
      target_status[j] <- analytic_status[[j]]
    } else {
      target_status[j] <- bootstrap_status[[j]]
    }
  }

  analytic_se <- ifelse(analytic_ok, sqrt(pmax(0.0, analytic_var)), NA_real_)
  names(analytic_se) <- target_names_chr

  all_targets_ok <- all(target_accepted) && cov_matrix_finite
  bootstrap_iqr_scale <- rep(NA_real_, n_targets)
  bootstrap_se <- rep(NA_real_, n_targets)
  names(bootstrap_iqr_scale) <- target_names_chr
  names(bootstrap_se) <- target_names_chr

  max_stat <- rep(NA_real_, n_draws_int)
  names(max_stat) <- draw_names_chr
  crit_val <- NA_real_

  if (all_targets_ok) {
    max_work <- abs(bres[, 1L]) / raw_bootstrap_scale[1L]
    if (n_targets > 1L) {
      for (j in seq_len(n_targets)[-1L]) {
        max_work <- pmax(max_work, abs(bres[, j]) / raw_bootstrap_scale[j])
      }
    }
    if (!all(is.finite(max_work))) {
      for (j in seq_len(n_targets)) {
        if (!all(is.finite(abs(bres[, j]) / raw_bootstrap_scale[j]))) {
          bootstrap_ok[j] <- FALSE
          bootstrap_status[j] <- "nonfinite_standardized_draw"
          target_accepted[j] <- FALSE
          target_status[j] <- "nonfinite_standardized_draw"
        }
      }
      all_targets_ok <- FALSE
    } else {
      bootstrap_iqr_scale <- raw_bootstrap_scale
      bootstrap_se <- candidate_bootstrap_se
      max_stat <- max_work
      names(max_stat) <- draw_names_chr
      sorted_max <- sort.int(max_stat, method = "quick")
      crit_val <- unname(sorted_max[k_crit])
    }
  }

  accepted <- isTRUE(all_targets_ok)
  refused <- !accepted
  status <- if (accepted) "ok" else "refused"
  failed_targets <- target_ids[!target_accepted]
  refusal_reasons <- if (accepted) {
    character(0)
  } else {
    unique(c(
      analytic_status[!analytic_ok],
      bootstrap_status[!bootstrap_ok]
    ))
  }

  operator_id <- if (isTRUE(center)) {
    "declared_cluster_sum_centered_outer_product_over_n2"
  } else {
    "declared_cluster_sum_uncentered_outer_product_over_n2"
  }
  omega_operator_id <- if (isTRUE(center)) {
    "declared_cluster_sum_centered_outer_product_over_n"
  } else {
    "declared_cluster_sum_uncentered_outer_product_over_n"
  }

  metadata <- list(
    row_level = "complete_target_if_unit_rows",
    cluster_aggregation = "unit_sum_within_declared_cluster",
    operator_id = operator_id,
    omega_operator_id = omega_operator_id,
    scale_id = "estimator_covariance_1_over_n2_bootstrap_se_sqrtC_over_n",
    quantile_type = 1L,
    quantile_definition = "empirical_order_statistic_ceil_Bp",
    iqr_normalizer = iqr_norm,
    n_units = n_units,
    n_clusters = n_clusters,
    n_nonempty_clusters = sum(streamed_cluster_unit_counts > 0L),
    n_zero_unit_clusters = sum(streamed_cluster_unit_counts == 0L),
    n_zero_sum_clusters = sum(rowSums(abs(S_mat)) == 0.0),
    cluster_unit_counts = streamed_cluster_unit_counts,
    n_denom = n_denom_num,
    n_targets = n_targets,
    target_ids = target_ids,
    cluster_ids = cluster_ids,
    expected_unit_ids = expected_unit_ids,
    expected_draw_ids = expected_draw_ids,
    n_draws = n_draws_int,
    n_if_batches = n_if_batches_int,
    n_sign_batches = n_sign_batches_int,
    center = center,
    alpha = alpha_num,
    return_raw_draws = return_raw_draws,
    budget = checked_budget,
    workspace_estimate = workspace_estimate,
    execution = list(workers = n_workers, route = parallel_route, parallel_rounds = n_rounds,
      sign_order = "every sign batch requested and validated in draw order by the calling process; workers multiply only"),
    caller_premises_only = TRUE,
    sampling_provenance_certified = FALSE,
    complete_if_or_nuisance_error_certified = FALSE,
    finite_sample_or_asymptotic_validity_certified = FALSE,
    entry_budget_is_rss_or_cpu_bound = FALSE,
    candidate_floating_arithmetic_only = TRUE
  )

  diagnostics <- list(
    analytic_variance_ok = analytic_ok,
    bootstrap_scale_ok = bootstrap_ok,
    analytic_status = analytic_status,
    bootstrap_status = bootstrap_status,
    cov_matrix_finite = cov_matrix_finite,
    failed_targets = failed_targets,
    refusal_reasons = refusal_reasons,
    order_stat_indices = list(
      q25_index = k_q25,
      q75_index = k_q75,
      crit_index = k_crit
    ),
    iqr_normalizer = iqr_norm,
    se_scale_factor = se_scale_factor,
    raw_draw_scale_factor = inv_sqrt_C
  )

  structure(
    list(
      accepted = accepted,
      refused = refused,
      status = status,
      target_ids = target_ids,
      target_status = target_status,
      target_accepted = target_accepted,
      cluster_ids = cluster_ids,
      cluster_unit_counts = streamed_cluster_unit_counts,
      cluster_sums = S_mat,
      cluster_means = cluster_means,
      cluster_sums_centered = S_centered,
      cov_mat = cov_mat,
      omega_cl_hat = omega_cl_hat,
      analytic_var = analytic_var,
      analytic_se = analytic_se,
      q25 = q25_vec,
      q75 = q75_vec,
      raw_bootstrap_iqr = raw_bootstrap_iqr,
      raw_bootstrap_scale = raw_bootstrap_scale,
      candidate_bootstrap_se = candidate_bootstrap_se,
      bootstrap_iqr_scale = bootstrap_iqr_scale,
      bootstrap_se = bootstrap_se,
      max_stat = max_stat,
      crit_val = crit_val,
      raw_draws = if (isTRUE(return_raw_draws)) bres else NULL,
      diagnostics = diagnostics,
      metadata = metadata
    ),
    class = "cdid_chained_multiplier_result"
  )
}
