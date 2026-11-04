# Optional ordinary static incomplete baseline OR replay. All original fits,
# supports and moment columns survive. Detailed/reference/complete routes stay
# on their existing evaluators; capacity decline returns their scalar stream.
or_sparse_native_live <- function(B, nl, nd, ns, ng, P, mp) {
  # Output, posting positions/endpoints/byte states (bytes charged as doubles),
  # row-local S/dy/up, owner beta scores, block alpha scores, one block score.
  1024 + B * (nl + 3 * nd) + 3 * ns + P + ng + mp
}

or_sparse_native_slice <- function(backend, pre, B, nl) {
  if (or_baseline_separate(backend))
    16 * B * (length(pre$sids) * (backend$outcome_recipe$p + 5) + nl) else
    # At most one quarter of the old 64*B*width pooled row reservation.
    16 * B * (4 * length(pre$sids) + 4 * length(pre$gids) + 4 * nl +
      4 * (backend$shared$p + 1) * (length(pre$beta_gids) + 1))
}

or_sparse_native_plan <- function(backend, links, B) {
  validate_or_baseline_backend(backend)
  if (!is.null(backend$reference_ref) || !is.null(backend$complete_support) ||
      !is.null(backend$separate_complete)) return(NULL)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) ||
      B < 1 || B != floor(B) || B > backend$control$max_chunk_rows)
    or_baseline_abort("Invalid sparse OR replay link batch or chunk rows.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  if (isTRUE(st$reference_carrier) || B > st$N || !is.double(frame$Z_a)) return(NULL)
  pre <- or_baseline_if_budget(backend, B, links, FALSE)
  if (any(vapply(pre[c("sids", "gids", "beta_gids")], function(x) !is.null(names(x)), logical(1)))) return(NULL)
  separate <- or_baseline_separate(backend); icpt <- isTRUE(backend$outcome_recipe$intercept)
  gids <- sort(union(pre$gids, pre$beta_gids)); ng <- length(gids); ns <- length(pre$sids); nl <- length(links)
  p <- if (separate) unname(vapply(backend$outcome_groups$groups[gids], function(g) length(g$keep), integer(1))) else
    rep.int(as.integer(backend$shared$p), ng)
  owner <- if (separate) seq_len(ng) else rep.int(1L, ng)
  op <- if (separate) p else as.integer(backend$shared$p)
  cn <- unname(vapply(backend$regression_specs[gids], function(s) length(s$support_ids), integer(1)))
  fb_n <- if (separate) cn else if (backend$shared$p) length(backend$shared$control_support_ids) else 0L
  link_b <- match(vapply(backend$links[links], `[[`, integer(1), "regression"), gids)
  if (anyNA(link_b)) or_baseline_abort("Sparse replay link lost its fitted outcome block.", "integrity")
  P <- sum(op); Q <- sum(p); C <- sum(cn); F <- sum(fb_n); V <- sum(op * fb_n); A <- sum(p[link_b])
  mp <- max(c(0, p)); nrhs <- tabulate(owner[link_b], nbins = length(op))
  dates <- sort(unique(as.integer(unlist(lapply(backend$support_specs[pre$sids], function(s)
    c(s$pair$s_calendar_index, s$pair$t_calendar_index)), use.names = FALSE))))
  flat <- 1 + 7 * ns + 10 * ng + 5 * length(op) + 3 * P + Q + 2 * C + F + V +
    7 * nl + A + length(pre$gids) + length(pre$beta_gids) + 3 * length(dates)
  retained <- 4096 + 4 * flat + 16 * (ns + ng + nl + length(op))
  peak <- 2048 + 16 * max(c(0, op * nrhs)) + 8 * mp^2 +
    8 * max(c(0, P, Q, C, F, V, A)) + 16 * (ns + ng + nl)
  tail <- st$N %% B; min_B <- if (tail) min(B, tail) else B
  slice <- or_sparse_native_slice(backend, pre, min_B, nl)
  live <- or_sparse_native_live(B, nl, length(dates), ns, ng, P, mp)
  # The old scalar buffer charge funds both halves. The plan is charged even
  # between chunks and at the smallest tail; no consumer admission is changed.
  construction_work <- 64 * flat + 16 * sum(op^2 * nrhs)
  if (any(!is.finite(c(retained, peak, slice, live, construction_work))) ||
      max(c(P, Q, C, F, V, A, ng, ns, nl)) >= .Machine$integer.max ||
      retained + peak > slice || live > or_sparse_native_slice(backend, pre, B, nl) ||
      retained + peak + live > pre$scratch || construction_work > pre$work / 4) return(NULL)
  start <- function(n) as.integer(cumsum(c(1, utils::head(n, -1L))))
  span <- function(s, n) if (n) seq.int(s, length.out = n) else integer(0)
  ps <- start(op); qs <- start(p); cs <- start(cn); fs <- start(fb_n); vs <- start(op * fb_n); ads <- start(p[link_b])
  data <- list(intercept = as.integer(icpt), support = matrix(0L, 3L, ns),
    support_par = matrix(0, 3L, ns), block = matrix(0L, 6L, ng), block_par = matrix(0, 4L, ng),
    owner = matrix(0L, 5L, length(op)), keep = integer(P), beta = numeric(P), center = numeric(P),
    mean_d = numeric(Q), controls = integer(C), forcing_alpha = numeric(C),
    forcing_support = integer(F), forcing_beta = numeric(V), link = matrix(0L, 3L, nl),
    link_par = matrix(0, 3L, nl), adj = numeric(A))
  for (j in seq_len(ns)) {
    s <- backend$support_specs[[pre$sids[j]]]; u <- backend$supports[[pre$sids[j]]]
    data$support[, j] <- as.integer(c(s$cohort, match(s$pair$s_calendar_index, dates), match(s$pair$t_calendar_index, dates)))
    data$support_par[, j] <- c(u$rate, u$massA, u$v)
  }
  for (kb in seq_len(ng)) {
    b <- gids[kb]; rec <- backend$regressions[[b]]; oi <- owner[kb]
    data$block[, kb] <- c(oi, qs[kb], p[kb], cs[kb], cn[kb], as.integer(b %in% pre$beta_gids))
    data$block_par[, kb] <- c(rec$L, rec$alpha, if (icpt) rec$mass else 1,
      if (separate && icpt) rec$score_mean_alpha else 0)
    data$mean_d[span(qs[kb], p[kb])] <- rec$mean_d
    data$controls[span(cs[kb], cn[kb])] <- match(backend$regression_specs[[b]]$support_ids, pre$sids)
    data$forcing_alpha[span(cs[kb], cn[kb])] <- if (icpt) rec$forcing_alpha else 0
    if (!separate) next
    f <- or_outcome_block_fit(backend$outcome_groups, b); ii <- span(ps[oi], op[oi])
    data$keep[ii] <- backend$outcome_recipe$slope_index[f$keep]
    data$beta[ii] <- f$beta; data$center[ii] <- rec$score_mean_beta
    data$forcing_support[span(fs[oi], fb_n[oi])] <- match(backend$regression_specs[[b]]$support_ids, pre$sids)
    data$forcing_beta[span(vs[oi], op[oi] * fb_n[oi])] <- as.double(rec$forcing_beta)
  }
  if (!separate) {
    data$keep <- as.integer(backend$outcome_recipe$slope_index[backend$shared$keep])
    data$beta <- as.double(backend$shared$beta)
    if (backend$shared$p) {
      data$forcing_support <- as.integer(match(backend$shared$control_support_ids, pre$sids))
      data$forcing_beta <- as.double(backend$shared$forcing_beta[, backend$shared$control_support_ids, drop = FALSE])
    }
  }
  for (oi in seq_along(op)) {
    data$owner[, oi] <- c(ps[oi], op[oi], fs[oi], fb_n[oi], vs[oi])
    jj <- which(owner[link_b] == oi)
    if (!length(jj) || !op[oi]) next
    f <- if (separate) or_outcome_block_fit(backend$outcome_groups, gids[oi]) else backend$shared
    loads <- vapply(jj, function(j) {
      l <- backend$links[[links[j]]]
      as.double(l$c - if (icpt) backend$regressions[[l$regression]]$mean_d else 0)
    }, numeric(op[oi]))
    dim(loads) <- c(op[oi], length(jj)); dimnames(loads) <- NULL
    adj <- if (icpt) or_pooled_wls_gram_solve(f$fit, loads, normalizer = st$N) else
      or_wls_gram_solve(f$fit, loads, normalizer = st$N)
    for (k in seq_along(jj)) data$adj[span(ads[jj[k]], op[oi])] <- adj[, k]
  }
  for (j in seq_len(nl)) {
    l <- backend$links[[links[j]]]
    data$link[, j] <- c(link_b[j], match(l$support_id, pre$sids), ads[j])
    data$link_par[, j] <- c(l$mu, l$d, l$rate_loading)
  }
  bounds <- rbind(panel$postings$date_start[dates], panel$postings$date_end[dates])
  empty <- panel$postings$date_count[dates] == 0L
  bounds[, empty] <- 0L; storage.mode(bounds) <- "integer"
  if (anyNA(data$controls) || anyNA(data$forcing_support) || anyNA(data$link))
    or_baseline_abort("Sparse replay descriptors lost their source supports.", "integrity")
  plan <- list(kind = "or_sparse_native_plan", backend_ref = backend$ref, links = as.integer(links),
    B = B, min_B = min_B, N = st$N, sids = pre$sids, gids = pre$gids, beta_gids = pre$beta_gids,
    dates = dates, bounds = bounds, data = data, retained_bound = retained,
    construction_peak = peak, reservation = retained + peak, entries = 0)
  plan$entries <- unname(cp04_payload(plan)["entries"]) + 64
  if (plan$entries > retained) or_baseline_abort("Sparse OR plan bound underestimated.", "bound")
  ipw_owner_seal(plan, "or_sparse_native_plan")
}

or_sparse_native_check <- function(plan, backend, links, B, pre) {
  or_baseline_native_check(plan, backend, links, B, pre, "or_sparse_native_plan", or_baseline_abort)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  if (!identical(plan$beta_gids, pre$beta_gids) || !isTRUE(plan$N == reader$state$N) ||
      !isTRUE(B == plan$B || B == plan$min_B) ||
      !isTRUE(plan$reservation <= or_sparse_native_slice(backend, pre, B, length(links))))
    or_baseline_abort("Sparse OR replay differs from its issued partition or global score blocks.", "integrity")
  invisible(plan)
}

or_sparse_native_if <- function(backend, reader, rows, links, scale, pre, plan) {
  or_sparse_native_check(plan, backend, links, length(rows), pre)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  live <- or_sparse_native_live(length(rows), length(links), length(plan$dates), length(pre$sids),
    ncol(plan$data$block), length(plan$data$beta), max(c(0L, plan$data$block[3L, ])))
  ipw_owner_safe_sum(c(plan$entries, live), pre$scratch, "sparse OR native replay")
  res <- .Call(or_baseline_native_kernel("C_gmm_or_sparse_if"), panel$observations$y,
    panel$postings$post_a_index, panel$postings$post_obs_index, as.integer(rows), st$cohort_index,
    st$a_index, st$w, plan$bounds, frame$Z_a, plan$data, as.double(live))
  if (!identical(res$status[1L], 0L)) or_baseline_native_refuse(res$status, or_baseline_abort)
  phi <- res$phi
  if (scale == "core") phi <- st$lambda * phi
  dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id"))
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_rate_score_evaluations = length(pre$sids),
      n_regression_score_evaluations = length(union(pre$gids, pre$beta_gids)),
      n_global_beta_blocks = if (or_baseline_separate(backend)) 0L else length(pre$beta_gids),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work, execution = "native_sparse_fixed_operator"))
  if (or_baseline_separate(backend)) {
    out$diagnostics$outcome_pooling <- "separate"
    out$diagnostics$n_local_beta_blocks <- length(pre$beta_gids)
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    or_baseline_abort("Sparse OR IF output bound underestimated.", "bound")
  out
}
