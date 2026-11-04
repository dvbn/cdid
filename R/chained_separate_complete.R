# Separate-only fixed-operator execution. Every original support, outcome block,
# treatment model and moment survives. The pooled complete state is not used.
or_separate_complete_facts <- function(reader, recipe, records, specs, ctrl, reference = NULL,
                                       frame_index = seq_along(recipe$columns), treatment = NULL,
                                       retained_external_entries = 0, external_work_ops = 0,
                                       metadata_entries = 0) {
  st <- reader$state
  if (!is.null(reference) || isTRUE(st$reference_carrier) || !isTRUE(recipe$intercept) ||
      !length(records) || !length(specs) || st$N != st$n_A) return(NULL)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  po <- get0("panel_owner", envir = reader$panel_ref, inherits = FALSE)
  if (!inherits(po, "cdid_ipw_panel_owner") || !identical(po$n_a_units, frame$n_a_units) ||
      st$n_A != po$n_a_units) return(NULL)
  used <- sort(unique(as.integer(c(frame_index, if (!is.null(treatment)) treatment$index))))
  dates <- sort(unique(as.integer(unlist(lapply(specs, function(s)
    c(s$pair$s_calendar_index, s$pair$t_calendar_index))))))
  B <- min(st$N, ctrl$max_chunk_rows)
  facts <- length(dates) * st$n_A + (length(used) + 2) * st$N + 8 * length(records)
  cp04_budget(ctrl, c(reader$accounting$owned_entries, reader$accounting$borrowed_owner_entries,
    retained_external_entries, metadata_entries), 16 * B * (length(used) + 4) + 32 * length(dates) + 4096,
    c(external_work_ops, 64 * facts), 4096 + 32 * length(dates), "separate fixed-operator facts")
  for (start in seq.int(1L, st$N, by = B)) {
    rows <- seq.int(start, min(st$N, start + B - 1L)); a <- st$a_index[rows]
    if (anyNA(a) || any(a < 1L | a > po$n_a_units) ||
        any(!is.finite(st$w[rows]) | st$w[rows] <= 0)) return(NULL)
    Z <- frame$Z_a[a, used, drop = FALSE]
    if (any(!is.finite(Z)) || any(frame$Z_a[a, frame_index[recipe$intercept_index]] != 1)) return(NULL)
  }
  for (d in dates) if (!isTRUE(ipw_owner_posting_date_complete(po$postings, po$n_a_units, d, B))) return(NULL)
  post <- po$postings
  list(route = "complete_support", panel_ref = reader$panel_ref, n_a = po$n_a_units,
    step = B, dates = dates, date_start = post$date_start[dates], date_end = post$date_end[dates],
    date_count = post$date_count[dates], fact_reads = facts)
}

# Compact source-owned lookup vectors. Requested batches subset these vectors;
# no row-by-block array or all-link copy is built during replay.
or_separate_complete_index <- function(complete, specs, links, treatments = NULL, families = NULL) {
  complete$support_pair <- vapply(specs, `[[`, integer(1), "pair_id")
  complete$support_cohort <- vapply(specs, `[[`, integer(1), "cohort")
  complete$link_block <- vapply(links, `[[`, integer(1), "regression")
  if (!is.null(treatments)) {
    complete$link_model <- vapply(links, `[[`, integer(1), "model")
    complete$treatment_block <- vapply(treatments, `[[`, integer(1), "regression")
    complete$treatment_family <- integer(length(treatments))
    for (j in seq_along(families)) complete$treatment_family[families[[j]]] <- j
  }
  complete
}

# Stable nonnegative combination of response-only residual coordinates. No
# squared-norm subtraction is used, including across early deficient chunks.
or_separate_complete_norm <- function(x) {
  if (!length(x)) return(0)
  if (any(!is.finite(x)))
    cdid_abort("Nonfinite compressed WLS factor.", "cdid_or_wls_numeric_error")
  scale <- max(abs(x))
  if (scale == 0) 0 else scale * sqrt(sum((x / scale)^2))
}

# One online unscaled weighted-design QR per cohort. Its row transforms serve
# bounded batches of RAW pair responses. Keep min(cumulative rows, 1+p) rows,
# never the numerical QR rank: an early zero-design response coordinate must
# survive until it actually becomes a trailing row. Pivoted predictor columns
# are restored in R, while exactly that QR's row transform is applied to h.
or_separate_complete_supports <- function(reader, specs, recipe, frame_index, ctrl, complete, scratch_bound) {
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- or_baseline_complete_panel(reader, complete)
  supports <- vector("list", length(specs))
  template <- or_pooled_support_start(recipe, 1)
  counts <- c(design_qr = 0L, common_design_merges = 0L, response_batches = 0L, response_columns = 0L)
  cohorts <- unique(vapply(specs, `[[`, integer(1), "cohort")); d <- recipe$p + 1L
  for (k in cohorts) {
    sids <- which(vapply(specs, function(s) s$cohort == k, logical(1)))
    first <- st$cohort_start[k]; last <- st$cohort_end[k]; M <- length(sids)
    B <- min(last - first + 1L, ctrl$max_chunk_rows)
    # h/e/sums and scaled-design facts coexist with one chunk's QR, up to
    # sixteen raw/projection/merge columns, small QR copies and final factors.
    # This peak must fit the already reserved S*p^2 separate-fit envelope.
    ipw_owner_safe_sum(c(8 * M * (d + 2), 16 * B * (d + 32), 128 * d^2, 8192),
      scratch_bound, "separate common-design merge scratch")
    h <- matrix(0, d, M); energy <- sumY <- numeric(M)
    R <- matrix(0, 0L, d); n_total <- 0L; mass <- 0
    sumZ <- colmax <- numeric(recipe$p); merge_rows <- 0L
    Ls <- vapply(specs[sids], `[[`, numeric(1), "transitions")
    distinct_L <- unique(Ls); L_pos <- match(Ls, distinct_L)
    anchors <- matrix(0, recipe$p, length(distinct_L))
    constant <- matrix(TRUE, recipe$p, length(distinct_L))
    for (start in seq.int(first, last, by = ctrl$max_chunk_rows)) {
      rows <- st$cohort_order[seq.int(start, min(last, start + ctrl$max_chunk_rows - 1L))]
      Z <- frame$Z_a[st$a_index[rows], frame_index, drop = FALSE]
      z <- Z[, recipe$slope_index, drop = FALSE]; w <- st$w[rows]; sw <- sqrt(w)
      X <- cbind(1, unname(z)); n <- length(rows)
      q <- qr(X * sw, LAPACK = TRUE)
      chunk_R <- qr.R(q)[, order(q$pivot), drop = FALSE]
      m_old <- nrow(R); m_chunk <- nrow(chunk_R)
      merge <- if (m_old) qr(rbind(R, chunk_R), LAPACK = TRUE) else NULL
      next_R <- if (m_old) qr.R(merge)[, order(merge$pivot), drop = FALSE] else chunk_R
      m_next <- min(n_total + n, d)
      if (nrow(next_R) != m_next || any(!is.finite(next_R)))
        cdid_abort("Nonfinite compressed WLS factor.", "cdid_or_wls_numeric_error")
      merge_rows <- max(merge_rows, m_old + m_chunk)
      counts["design_qr"] <- counts["design_qr"] + 1L
      if (m_old) counts["common_design_merges"] <- counts["common_design_merges"] + 1L
      for (lo in seq.int(1L, M, by = 16L)) {
        ids <- seq.int(lo, min(M, lo + 15L))
        raw <- matrix(0, n, length(ids))
        for (j in seq_along(ids)) raw[, j] <- or_baseline_complete_pair_rows(reader,
          panel, rows, specs[[sids[ids[j]]]]$pair, TRUE)$dy
        projected <- qr.qty(q, raw * sw)
        leading <- projected[seq_len(m_chunk), , drop = FALSE]
        transformed <- if (m_old) qr.qty(merge,
          rbind(h[seq_len(m_old), ids, drop = FALSE], leading)) else leading
        if (any(!is.finite(projected)) || any(!is.finite(transformed)))
          cdid_abort("Nonfinite compressed WLS factor.", "cdid_or_wls_numeric_error")
        h[seq_len(m_next), ids] <- transformed[seq_len(m_next), , drop = FALSE]
        for (j in seq_along(ids)) {
          tail_new <- if (n > m_chunk) projected[seq.int(m_chunk + 1L, n), j] else numeric(0)
          tail_merge <- if (m_old + m_chunk > m_next)
            transformed[seq.int(m_next + 1L, m_old + m_chunk), j] else numeric(0)
          energy[ids[j]] <- or_separate_complete_norm(c(energy[ids[j]],
            or_separate_complete_norm(tail_new), or_separate_complete_norm(tail_merge)))
          sumY[ids[j]] <- sumY[ids[j]] + sum(w * raw[, j])
        }
        counts["response_batches"] <- counts["response_batches"] + 1L
        counts["response_columns"] <- counts["response_columns"] + length(ids)
      }
      # Exact constancy is a fact of raw D=L*z, not of a rounded QR factor.
      # Equal actual L values reuse this fact without equating fitted weights.
      for (j in seq_along(distinct_L)) {
        D <- distinct_L[j] * z
        if (any(!is.finite(D))) or_baseline_abort("Nonfinite separate complete pair design.", "numeric")
        if (!n_total && recipe$p) anchors[, j] <- as.numeric(D[1L, ])
        if (recipe$p) constant[, j] <- constant[, j] & vapply(seq_len(recipe$p),
          function(col) all(D[, col] == anchors[col, j]), logical(1))
      }
      mass <- mass + sum(w); sumZ <- sumZ + unname(colSums(w * z))
      if (recipe$p) colmax <- pmax(colmax, unname(apply(abs(z), 2L, max)))
      R <- next_R; n_total <- n_total + n
    }
    for (j in seq_along(sids)) {
      # Pair duration scales predictor columns only AFTER the common merge.
      F <- cbind(sweep(R, 2L, c(1, rep(Ls[j], recipe$p)), `*`), h[seq_len(nrow(R)), j])
      if (n_total > d) F <- rbind(F, c(rep(0, d), energy[j]))
      if (any(!is.finite(F)))
        cdid_abort("Nonfinite compressed WLS factor.", "cdid_or_wls_numeric_error")
      raw_summary <- structure(list(factor = unname(F), n = n_total,
        columns = template$summary$raw$columns), class = "cdid_or_wls_summary")
      summary <- structure(list(raw = raw_summary, n = n_total, p = recipe$p,
        columns = template$summary$columns, anchor = if (recipe$p) anchors[, L_pos[j]] else numeric(0),
        constant = if (recipe$p) constant[, L_pos[j]] else logical(0),
        accounting = list(input_rows = n_total, augmented_cols = recipe$p + 2L,
          augmented_rows = nrow(F), max_chunk_merge_rows = merge_rows,
          n_compressions = 0L)), class = "cdid_or_pooled_wls_chunk")
      acc <- list(L = Ls[j], summary = summary, weight_mass = mass, n = n_total,
        sumZ = sumZ, sumY = sumY[j], colmax = colmax)
      supports[[sids[j]]] <- or_baseline_support_rate(acc, mass, n_total, st$N)
    }
  }
  list(supports = supports, counts = counts)
}

or_separate_complete_extra <- function(complete, B, p, S, G, L) {
  list(models = 8192 + 64 * (length(complete$dates) + S + G + L),
    scratch = 128 * B * (p + 20) + 128 * (p + 2)^2 + 64 * S,
    work = 128 * complete$fact_reads + 128 * (S + G + L) * (p + 2)^2)
}

# Only one bounded original-unit chunk. Raw pair responses precede all sums.
or_separate_complete_hits <- function(backend, reader, rows, sids) {
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- or_baseline_complete_panel(reader, backend$separate_complete)
  columns <- if (inherits(backend, "cdid_dr_backend")) backend$or_frame_index else
    seq_along(backend$outcome_recipe$columns)
  Z <- frame$Z_a[st$a_index[rows], columns[backend$outcome_recipe$slope_index], drop = FALSE]
  ci <- st$cohort_index[rows]; fixed <- backend$separate_complete
  cohorts <- unique(fixed$support_cohort[sids])
  cohort_rows <- split(seq_along(rows), factor(ci, levels = cohorts))
  cohort_pos <- match(fixed$support_cohort[sids], cohorts)
  pair_groups <- split(seq_along(sids), fixed$support_pair[sids])
  hits <- vector("list", length(sids))
  for (jj in pair_groups) {
    dy <- or_baseline_complete_pair_rows(reader, panel, rows,
      backend$support_specs[[sids[jj[1L]]]]$pair, TRUE)$dy
    for (j in jj) {
      ix <- cohort_rows[[cohort_pos[j]]]
      hits[[j]] <- list(obs = ix, dy = dy[ix], z = Z[ix, , drop = FALSE])
    }
  }
  hits
}

# The same coefficient IF contracted by the requested scalar loadings. Explicit
# score centering survives, and standalone OR still has its intercept IF.
or_separate_complete_adjoint <- function(backend, b, f, hits, sids, w, loadings, alpha = FALSE) {
  rec <- backend$regressions[[b]]
  B <- length(w); nl <- ncol(loadings)
  adj <- if (f$p) {
    if (backend$outcome_recipe$intercept) or_pooled_wls_gram_solve(f$fit, loadings,
      normalizer = backend$outcome_groups$N) else or_wls_gram_solve(f$fit, loadings,
      normalizer = backend$outcome_groups$N)
  } else matrix(0, 0L, nl)
  value <- matrix(0, B, nl); ha <- numeric(B)
  for (sid in backend$regression_specs[[b]]$support_ids) {
    h <- hits[[match(sid, sids)]]; ix <- h$obs
    d <- f$L * h$z[, f$keep, drop = FALSE]
    r <- h$dy - f$alpha - as.vector(d %*% f$beta)
    wr <- w[ix] * r
    if (alpha) ha[ix] <- ha[ix] + wr
    if (f$p && length(ix)) value[ix, ] <- value[ix, , drop = FALSE] +
      wr * (sweep(d, 2L, f$mean_d) %*% adj)
  }
  if (f$p) value <- sweep(value, 2L, as.vector(rec$score_mean_beta %*% adj))
  if (alpha) value <- value + (ha - rec$score_mean_alpha) / rec$mass
  value
}

or_separate_complete_if <- function(backend, reader, rows, links, scale, pre) {
  st <- reader$state; B <- length(rows); w <- st$w[rows]
  hits <- or_separate_complete_hits(backend, reader, rows, pre$sids)
  phi <- matrix(0, B, length(links))
  by_block <- split(seq_along(links), factor(backend$separate_complete$link_block[links], levels = pre$gids))
  for (k in seq_along(pre$gids)) {
    b <- pre$gids[k]; jj <- by_block[[k]]
    f <- or_outcome_block_fit(backend$outcome_groups, b)
    loads <- matrix(0, f$p, length(jj))
    if (f$p) for (j in seq_along(jj)) loads[, j] <- backend$links[[links[jj[j]]]]$c - f$mean_d
    phi[, jj] <- -or_separate_complete_adjoint(backend, b, f, hits, pre$sids, w, loads, TRUE)
    for (j in jj) {
      l <- backend$links[[links[j]]]; h <- hits[[match(l$support_id, pre$sids)]]; ix <- h$obs
      d <- f$L * h$z[, f$keep, drop = FALSE]
      r <- h$dy - f$alpha - as.vector(d %*% f$beta)
      phi[ix, j] <- phi[ix, j] + w[ix] * (r - l$mu) / l$d
    }
  }
  or_separate_complete_output(backend, st, rows, links, scale, pre, phi,
    length(pre$gids), 0L, 0L)
}

dr_separate_complete_if <- function(backend, reader, rows, links, scale, pre) {
  st <- reader$state; B <- length(rows); w <- st$w[rows]
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  hits <- or_separate_complete_hits(backend, reader, rows, pre$sids)
  phi <- matrix(0, B, length(links))
  fixed <- backend$separate_complete
  by_model <- split(seq_along(links), factor(fixed$link_model[links], levels = pre$tids))
  names(by_model) <- as.character(pre$tids)
  by_block <- split(pre$tids, factor(fixed$treatment_block[pre$tids], levels = pre$gids))
  block_of_model <- match(fixed$treatment_block[pre$tids], pre$gids)
  names(block_of_model) <- as.character(pre$tids)
  block_fits <- vector("list", length(pre$gids))
  for (k in seq_along(pre$gids)) {
    b <- pre$gids[k]; tids <- by_block[[k]]
    f <- or_outcome_block_fit(backend$outcome_groups, b)
    block_fits[[k]] <- f
    loads <- matrix(0, f$p, length(tids))
    if (f$p) for (j in seq_along(tids)) loads[, j] <- backend$treatments[[tids[j]]]$d_bar_diff
    component <- -or_separate_complete_adjoint(backend, b, f, hits, pre$sids, w, loads)
    for (j in seq_along(tids)) for (ell in by_model[[as.character(tids[j])]]) phi[, ell] <- component[, j]
  }
  families <- split(pre$tids, fixed$treatment_family[pre$tids])
  for (tids in families) {
    first <- backend$treatment_specs[[tids[1L]]]; first_fit <- backend$treatments[[tids[1L]]]
    rr <- which(st$cohort_index[rows] %in% first$cohorts)
    X <- frame$Z_a[st$a_index[rows[rr]], backend$treatment_recipe$index, drop = FALSE]
    lv <- ipw_logit_values(X, first_fit$coefficients)
    treated <- as.double(st$cohort_index[rows[rr]] == first$cohorts[1L])
    wp <- wm <- numeric(B); wp[rr] <- w[rr] * treated; wm[rr] <- w[rr] * (1 - treated) * lv$odds
    score <- matrix(0, B, backend$treatment_recipe$p)
    score[rr, ] <- w[rr] * (treated - lv$prob) * X
    gamma <- t(cp04_solve(first_fit$H_chol, t(score)))
    for (tid in tids) {
      spec <- backend$treatment_specs[[tid]]; f <- backend$treatments[[tid]]
      if (!identical(spec$cohorts, first$cohorts) || !identical(f$coefficients, first_fit$coefficients) ||
          !identical(f$H_chol, first_fit$H_chol))
        dr_baseline_abort("Separate complete propensity replay differs from its fitted operator.", "identity")
      of <- block_fits[[block_of_model[as.character(tid)]]]
      direct <- numeric(B)
      for (sid in spec$support_ids) {
        h <- hits[[match(sid, pre$sids)]]; ix <- h$obs
        d <- of$L * h$z[, of$keep, drop = FALSE]
        r <- h$dy - of$alpha - as.vector(d %*% of$beta)
        direct[ix] <- wp[ix] * (r - f$means["T_dr"]) / f$means["dp"] -
          wm[ix] * (r - f$means["M_dr"]) / f$means["dm"]
      }
      value <- direct - as.vector(gamma %*% f$partials$beta[2L, ])
      for (j in by_model[[as.character(tid)]]) phi[, j] <- phi[, j] + value
    }
  }
  or_separate_complete_output(backend, st, rows, links, scale, pre, phi,
    length(pre$gids), length(pre$tids), length(families))
}

or_separate_complete_output <- function(backend, st, rows, links, scale, pre, phi, ng, nt, nf) {
  cp04_finite(phi, "separate complete original-unit IF")
  dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id"))
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(outcome_pooling = "separate", execution = "fixed_operator",
      n_rate_score_evaluations = 0L, n_structural_zero_rate_scores = length(pre$sids),
      n_regression_score_evaluations = ng, n_global_beta_blocks = 0L, n_block_beta_blocks = ng,
      n_treatment_score_evaluations = nt, n_distinct_treatment_score_evaluations = nf,
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    or_baseline_abort("Separate complete IF output bound underestimated.", "bound")
  out
}
