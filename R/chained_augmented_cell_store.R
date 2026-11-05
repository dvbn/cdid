# Ragged score/Z data for one already-admitted covariance producer. All
# selection/rate/upstream IF and final assembly retain their original B.
augmented_cell_replay_store <- function(backend, pre, phase_store) {
  if (is.null(phase_store) || !is.environment(phase_store) ||
      !identical(phase_store$kind, "cp04_retained_if") ||
      !identical(phase_store$phase, "covariance") || !identical(phase_store$state, "reserved") ||
      !identical(phase_store$backend_ref, backend$ref) ||
      exists("producer_cell_accounting", phase_store, inherits = FALSE) ||
      !identical(backend$outcome_pooling, "pooled") || !identical(backend$imputation_pooling, "separate") ||
      !is.null(backend$reference_ref) || !is.null(backend$layout) || !is.null(backend$pair_outcomes)) return(NULL)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  X <- augmented_backend_X_F(backend)
  if (!is.null(rownames(X))) return(NULL)
  ids <- sort(unique(c(pre$cids, pre$zids))); nc <- length(backend$cells)
  if (!length(ids)) return(NULL)
  cohort_n <- st$cohort_end - st$cohort_start + 1L
  width <- vapply(backend$cells[ids], function(cf) backend$groups[[cf$group_id]]$p, integer(1L))
  observed <- vapply(backend$cells[ids], `[[`, integer(1L), "n_pair")
  cohorts <- vapply(backend$cells[ids], `[[`, integer(1L), "cohort_index")
  names_bytes <- sum(vapply(backend$cells[ids], function(cf)
    sum(nchar(backend$groups[[cf$group_id]]$kept_columns, type = "bytes")), numeric(1L)))
  B_build <- min(backend$control$max_chunk_rows, max(cohort_n[cohorts]))
  numeric_entries <- sum(observed * (width + 2) + cohort_n[cohorts])
  retained <- 8192 + 8 * st$N + 1024 * nc + 64 * names_bytes + numeric_entries
  peak <- 8192 + 128 * B_build * (ncol(X) + ncol(st$Z) + max(c(0L, width)) + 32) +
    64 * max(c(0L, width))^2 + 64 * names_bytes
  # Retain a conservative potential numeric-slab duplicate as well as the
  # largest legal one-cell builder, while raw link IF and consumers coexist.
  total <- ipw_owner_safe_sum(c(phase_store$admission_total, 2 * retained, peak))
  build_work <- ipw_owner_safe_product(256, ipw_owner_safe_sum(c(numeric_entries, retained, peak)))
  work_total <- ipw_owner_safe_sum(c(phase_store$work_total, build_work))
  cap <- min(phase_store$cap, backend$control$max_retained_entries, backend$control$max_matrix_entries)
  if (total > cap || work_total > backend$control$max_work_ops) return(NULL)
  storage <- new.env(parent = emptyenv())
  storage$records <- vector("list", nc); storage$rank <- integer(st$N)
  storage$lo <- storage$hi <- storage$last <- integer(nc); storage$active <- TRUE
  for (cohort in unique(cohorts)) {
    rr <- st$cohort_order[seq.int(st$cohort_start[cohort], st$cohort_end[cohort])]
    storage$rank[rr] <- seq_along(rr)
  }
  op <- selection_prediction_operation(backend$selection_backend, source_reader = reader)
  on.exit(op$abort(), add = TRUE)
  batches <- 0L
  for (cid in ids) {
    cf <- backend$cells[[cid]]; gf <- backend$groups[[cf$group_id]]; cohort <- cf$cohort_index
    n <- cohort_n[cohort]; np <- cf$n_pair; need_score <- cid %in% pre$cids; need_Z <- cid %in% pre$zids
    global <- st$cohort_order[seq.int(st$cohort_start[cohort], st$cohort_end[cohort])]
    z <- if (need_Z) numeric(n) else NULL
    score_rows <- if (need_score) integer(np) else NULL
    score <- if (need_score) numeric(np) else NULL
    d <- if (need_score) matrix(0, np, gf$p,
      dimnames = if (gf$p) list(NULL, gf$kept_columns) else NULL) else NULL
    filled <- 0L
    for (first in seq.int(1L, n, by = B_build)) {
      local <- seq.int(first, min(n, first + B_build - 1L)); rows <- global[local]
      jhat <- if (need_Z) {
        if (gf$p) cf$alpha + as.vector((cf$L * X[rows, gf$kept_cols_F, drop = FALSE]) %*% gf$beta) else
          rep(cf$alpha, length(rows))
      } else NULL
      hit <- cp04_pair_rows(reader, rows, cf$pair, TRUE); obs <- which(hit$S == 1)
      if (need_Z) {
        value <- jhat
        if (length(obs)) {
          q <- op$pair(rows[obs], cf$q_id)$prob
          cp04_finite(q, "fitted observation probabilities on observed pairs", TRUE)
          pi <- 1 / ipw_clip_floor(q)
          e <- hit$dy[obs] - jhat[obs]
          value[obs] <- jhat[obs] + pi * e
        }
        z[local] <- value
      }
      if (need_score && length(obs)) {
        observed_rows <- rows[obs]
        d_obs <- if (gf$p) cf$L * X[observed_rows, gf$kept_cols_F, drop = FALSE] else
          matrix(0, length(obs), 0L)
        j_obs <- if (gf$p) cf$alpha + as.vector(d_obs %*% gf$beta) else rep(cf$alpha, length(obs))
        at <- seq.int(filled + 1L, filled + length(obs))
        score_rows[at] <- observed_rows; d[at, ] <- d_obs
        score[at] <- st$w[observed_rows] * (hit$dy[obs] - j_obs)
        filled <- filled + length(obs)
      }
      batches <- batches + 1L
    }
    if (need_score && filled != np)
      augmented_abort("Stored observed cell rows disagree with the original fitted support.", "integrity")
    storage$records[[cid]] <- list(rows = score_rows, d = d, score = score, Z = z, cohort = cohort)
  }
  op$finish()
  actual <- unname(cp04_payload(list(records = storage$records, rank = storage$rank,
    lo = storage$lo, hi = storage$hi, last = storage$last))["entries"])
  if (actual > retained) augmented_abort("Compact cell replay storage bound underestimated.", "bound")
  phase_store$producer_cell_peak <- total
  phase_store$producer_cell_work <- build_work
  phase_store$producer_cell_accounting <- list(retained_entries = retained, actual_entries = actual,
    duplicate_entries = retained, construction_peak = peak, admission_total = total,
    cap = cap, build_rows = B_build, build_batches = batches, n_cells = length(ids),
    observed_rows = sum(observed), cohort_rows = sum(cohort_n[cohorts]), work_ops_bound = build_work)
  augmented_cell_replay_methods(storage, st$cohort_index, st$N, phase_store$B,
    phase_store$producer_cell_accounting)
}

augmented_cell_replay_methods <- function(storage, cohort_index, N, B, accounting) {
  force(storage); force(cohort_index); force(N); force(B); force(accounting)
  matches <- function(rows) length(rows) > 0L &&
    identical(as.integer(rows), seq.int(rows[1L], min(N, rows[1L] + B - 1L))) &&
    (rows[1L] - 1L) %% B == 0L
  score <- function(cid, rows) {
    if (!storage$active) augmented_abort("Compact cell replay storage is closed.", "integrity")
    rec <- storage$records[[cid]]; first <- rows[1L]; last <- rows[length(rows)]
    if (is.null(rec$rows)) augmented_abort("Cell score outside the issued stored domain.", "integrity")
    lo <- if (first >= storage$last[cid]) max(1L, storage$lo[cid]) else 1L
    n <- length(rec$rows)
    while (lo <= n && rec$rows[lo] < first) lo <- lo + 1L
    hi <- lo - 1L
    while (hi < n && rec$rows[hi + 1L] <= last) hi <- hi + 1L
    storage$lo[cid] <- lo; storage$hi[cid] <- hi; storage$last[cid] <- first
    if (hi < lo) return(NULL)
    at <- seq.int(lo, hi); rr <- match(rec$rows[at], rows)
    list(rr = as.integer(rr), d_rr = rec$d[at, , drop = FALSE], s_rr = rec$score[at])
  }
  Z <- function(cid, rows) {
    if (!storage$active) augmented_abort("Compact cell replay storage is closed.", "integrity")
    rec <- storage$records[[cid]]
    if (is.null(rec$Z)) augmented_abort("Augmented response outside the issued stored domain.", "integrity")
    ix <- which(cohort_index[rows] == rec$cohort)
    list(ix = ix, Z = rec$Z[storage$rank[rows[ix]]])
  }
  release <- function() {
    storage$active <- FALSE
    storage$records <- storage$rank <- storage$lo <- storage$hi <- storage$last <- NULL
    invisible(NULL)
  }
  list(matches = matches, score = score, Z = Z, release = release, accounting = accounting)
}
