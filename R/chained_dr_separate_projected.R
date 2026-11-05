# Ordinary baseline, non-reference separate DR only. The caller has already
# validated and admitted its original row/support/rate-score construction.
# Detailed, conditional, reference and complete-support replay stay unchanged.

dr_separate_plan_range <- function(start, count) {
  if (count) seq.int(start, length.out = count) else integer(0L)
}

# The original count formulas, before any fixed operators are constructed.
# Their scalar/vector records fit the dependency-map reservation already
# issued by dr_separate_if_envelope; native admission also counts them while
# either construction is live.
dr_separate_baseline_replay_counts <- function(backend, links, B, pre, N) {
  count_entries <- 2048 + 32 * (length(links) + length(pre$sids) + length(pre$gids) + length(pre$tids))
  ipw_owner_safe_sum(count_entries, pre$scratch, "separate DR replay count metadata")
  nl <- length(links); ns <- length(pre$sids); ng <- length(pre$gids); nt <- length(pre$tids)
  pt <- backend$treatment_recipe$p
  p <- unname(vapply(backend$outcome_groups$groups[pre$gids], function(g) length(g$keep), integer(1L)))
  control_n <- unname(vapply(backend$regression_specs[pre$gids], function(s) length(s$support_ids), integer(1L)))
  model_b <- match(vapply(backend$treatment_specs[pre$tids], `[[`, integer(1L), "regression"), pre$gids)
  model_n <- unname(vapply(backend$treatment_specs[pre$tids], function(s) length(s$support_ids), integer(1L)))
  if (anyNA(model_b)) dr_baseline_abort("Requested model is outside its frozen outcome block.", "integrity")
  R <- sum(p[model_b]); C <- sum(control_n[model_b]); M <- sum(model_n); I <- sum(control_n)
  peak_p <- max(c(0, p)); peak_c <- max(c(0, control_n)); peak_m <- max(c(0, model_n))
  per_block <- tabulate(model_b, nbins = ng)
  peak_rhs <- max(c(0, p * per_block))
  retained <- 4096 + 2 * nl + 2 * ns + 6 * ng + 10 * nt + R + C + pt * nt + 2 * M + I
  # Complex field/subscript replacement may copy a complete flat operator.
  # Reserve old/new copies even where this R build can update it in place.
  peak <- 1024 + 4 * (nl + nt + ng) + 16 * peak_rhs + 4 * peak_p^2 +
    4 * peak_p * peak_c + 8 * pt^2 + 16 * pt + 4 * peak_m +
    2 * max(c(0, R, C, M, I, pt * nt))
  tail <- N %% B; min_B <- if (tail) min(B, tail) else B
  bound <- retained + peak
  list(nl = nl, ns = ns, ng = ng, nt = nt, pt = pt, p = p,
    control_n = control_n, model_b = model_b, model_n = model_n,
    R = R, C = C, M = M, I = I, per_block = per_block,
    retained = retained, peak = peak, min_B = min_B, bound = bound,
    count_entries = count_entries)
}

# Fixed operators of one already issued replay, never row data or fit copies.
# Count the retained descriptors AND construction peak before allocating them.
# Half of the old 32*B*nl allowance still covers phi/projected/scaled buffers;
# all other original reservations remain intact. The smallest issued tail
# protects the plan's lifetime between chunks as well as each full chunk.
dr_separate_baseline_replay_plan <- function(backend, links, B, native = FALSE) {
  validate_dr_baseline_backend(backend)
  if (!identical(native, isTRUE(native))) dr_baseline_abort("Invalid private replay funding mode.")
  if (!dr_separate(backend) || isTRUE(backend$uses_selection_models) ||
      !is.null(backend$reference_ref) || !is.null(backend$separate_complete)) return(NULL)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) ||
      B < 1 || B != floor(B) || B > backend$control$max_chunk_rows)
    dr_baseline_abort("Invalid separate DR replay link batch or chunk rows.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N
  if (B > N) return(NULL)
  pre <- dr_baseline_if_budget(backend, B, links, FALSE)
  # Descriptors carry integer identities, never replicated source labels.
  if (any(vapply(pre[c("sids", "uids", "gids", "tids")], function(x) !is.null(names(x)), logical(1L)))) return(NULL)
  counts <- dr_separate_baseline_replay_counts(backend, links, B, pre, N)
  joint <- if (native) dr_separate_sparse_prospective(backend, links, B, pre, counts) else NULL
  nl <- counts$nl; ns <- counts$ns; ng <- counts$ng; nt <- counts$nt; pt <- counts$pt
  p <- counts$p; control_n <- counts$control_n; model_b <- counts$model_b; model_n <- counts$model_n
  R <- counts$R; C <- counts$C; M <- counts$M; I <- counts$I; per_block <- counts$per_block
  retained <- counts$retained; peak <- counts$peak; min_B <- counts$min_B; bound <- counts$bound
  rm(counts)
  # A failed optional prospective check keeps the exact old R admission.
  if (is.null(joint) && (!is.finite(bound) || bound > 16 * min_B * nl)) return(NULL)
  starts <- function(n) cumsum(c(1, utils::head(n, -1L)))
  model_order <- order(model_b, seq_along(model_b))
  r_n <- p[model_b]; c_n <- control_n[model_b]
  rm(model_b)
  r_start <- c_start <- numeric(nt)
  r_start[model_order] <- starts(r_n[model_order])
  c_start[model_order] <- starts(c_n[model_order])
  rm(r_n, c_n)
  requested_models <- match(vapply(backend$links[links], `[[`, integer(1L), "model"), pre$tids)
  if (anyNA(requested_models)) dr_baseline_abort("Requested link lost its frozen treatment model.", "integrity")
  link_n <- tabulate(requested_models, nbins = nt)
  link_order <- order(requested_models, seq_along(requested_models))
  rm(requested_models)
  plan <- list(kind = "dr_separate_baseline_replay_plan", backend_ref = backend$ref,
    links = as.integer(links), B = B, min_B = min_B, N = N,
    sids = pre$sids, uids = pre$uids, gids = pre$gids, tids = pre$tids,
    p = p, model_order = model_order, block_start = starts(per_block), block_count = per_block,
    control_start = starts(control_n), control_count = control_n, control_pos = numeric(I),
    link_order = link_order, link_start = starts(link_n), link_count = link_n,
    adj_start = r_start, or_adj = numeric(R), or_center = numeric(nt),
    or_up_start = c_start, or_up = numeric(C), eta = matrix(0, pt, nt),
    support_start = starts(model_n), support_count = model_n, support_pos = numeric(M),
    prop_up = numeric(M), retained_bound = retained, construction_peak = peak,
    reservation = bound, entries = 0)
  rm(model_order, r_start, c_start, link_order, link_n, per_block, p, control_n, model_n)
  for (kb in seq_len(ng)) {
    b <- pre$gids[kb]
    models <- plan$model_order[dr_separate_plan_range(plan$block_start[kb], plan$block_count[kb])]
    tids <- pre$tids[models]
    controls <- backend$regression_specs[[b]]$support_ids
    plan$control_pos[dr_separate_plan_range(plan$control_start[kb], plan$control_count[kb])] <- match(controls, pre$sids)
    of <- or_outcome_block_fit(backend$outcome_groups, b)
    if (!of$p) next
    loads <- vapply(backend$treatments[tids], function(f) as.double(f$d_bar_diff), numeric(of$p))
    dim(loads) <- c(of$p, length(models)); dimnames(loads) <- NULL
    adj <- if (backend$outcome_recipe$intercept)
      or_pooled_wls_gram_solve(of$fit, loads, normalizer = backend$outcome_groups$N) else
      or_wls_gram_solve(of$fit, loads, normalizer = backend$outcome_groups$N)
    rec <- backend$regressions[[b]]
    for (j in seq_along(models)) {
      k <- models[j]
      plan$or_adj[dr_separate_plan_range(plan$adj_start[k], of$p)] <- adj[, j]
      plan$or_center[k] <- as.vector(rec$score_mean_beta %*% adj[, j, drop = FALSE])
      plan$or_up[dr_separate_plan_range(plan$or_up_start[k], length(controls))] <-
        as.vector(t(rec$forcing_beta) %*% adj[, j, drop = FALSE])
    }
  }
  rm(of)
  for (k in seq_len(nt)) {
    tid <- pre$tids[k]; f <- backend$treatments[[tid]]
    eta <- as.vector(cp04_solve(f$H_chol, f$partials$beta[2L, ]))
    ii <- dr_separate_plan_range(plan$support_start[k], plan$support_count[k])
    plan$eta[, k] <- eta
    plan$support_pos[ii] <- match(backend$treatment_specs[[tid]]$support_ids, pre$sids)
    plan$prop_up[ii] <- as.vector(t(f$cross) %*% eta)
  }
  if (anyNA(plan$control_pos) || anyNA(plan$support_pos))
    dr_baseline_abort("Replay support descriptors differ from the frozen model supports.", "integrity")
  plan$entries <- unname(cp04_payload(plan)["entries"]) + 64
  if (plan$entries > retained) dr_baseline_abort("Separate DR replay plan bound underestimated.", "bound")
  plan <- ipw_owner_seal(plan, "dr_separate_baseline_replay_plan")
  if (is.null(joint)) return(plan)
  # No native-only operator can escape as an R fallback. All eligibility and
  # counts were admitted before the operators above; any later underestimate
  # is an error, never a second producer after construction/execution.
  out <- dr_separate_sparse_plan(backend, links, B, plan)
  if (is.null(out) || !identical(out$kind, "dr_separate_sparse_plan"))
    dr_baseline_abort("The admitted joint separate DR construction changed its inventory.", "bound")
  out
}

dr_separate_baseline_replay_identity <- function(plan, backend, links, B, pre) {
  or_baseline_native_check(plan, backend, links, B, pre,
    "dr_separate_baseline_replay_plan", dr_baseline_abort)
  N <- get("cp04_reader", backend$reader_ref, inherits = FALSE)$state$N
  tail <- N %% plan$B; min_B <- if (tail) min(plan$B, tail) else plan$B
  if (!identical(plan$uids, pre$uids) || !isTRUE(plan$N == N) || !isTRUE(plan$min_B == min_B) ||
      !isTRUE(B == plan$B || B == plan$min_B))
    dr_baseline_abort("The separate DR replay plan differs from its issued row partition or upstream IDs.", "integrity")
  invisible(plan)
}

dr_separate_baseline_replay_check <- function(plan, backend, links, B, pre) {
  dr_separate_baseline_replay_identity(plan, backend, links, B, pre)
  if (!isTRUE(plan$reservation <= 16 * B * length(links)))
    dr_baseline_abort("The separate DR replay plan differs from its issued row partition or upstream IDs.", "integrity")
  invisible(plan)
}

# Contract the SAME frozen local coefficient IF with requested scalar model
# loadings. Keep score centering and nonzero rate forcing on every F row.
# Loadings/adjoints are p x nt; rate forcing uses one support vector per model,
# never a support-by-model matrix outside the existing conservative envelope.
dr_separate_baseline_adjoint <- function(backend, b, f, tids, hits, up, sids, w,
                                       plan = NULL, kb = NULL, model_positions = NULL) {
  B <- length(w); nt <- length(tids); p <- f$p
  value <- matrix(0, B, nt)
  if (!p) return(value)
  rec <- backend$regressions[[b]]
  if (is.null(plan)) {
    loads <- matrix(0, p, nt)
    for (j in seq_along(tids)) loads[, j] <- backend$treatments[[tids[j]]]$d_bar_diff
    adj <- if (backend$outcome_recipe$intercept)
      or_pooled_wls_gram_solve(f$fit, loads, normalizer = backend$outcome_groups$N) else
      or_wls_gram_solve(f$fit, loads, normalizer = backend$outcome_groups$N)
  } else adj <- matrix(plan$or_adj[dr_separate_plan_range(plan$adj_start[model_positions[1L]], p * nt)], p, nt)
  controls <- backend$regression_specs[[b]]$support_ids
  pos <- if (is.null(plan)) match(controls, sids) else
    plan$control_pos[dr_separate_plan_range(plan$control_start[kb], plan$control_count[kb])]
  for (k in seq_along(controls)) {
    sid <- controls[k]; h <- hits[[pos[k]]]; ix <- h$obs
    if (!length(ix)) next
    d <- f$L * h$z[, f$keep, drop = FALSE]
    residual <- h$dy - f$alpha - as.vector(d %*% f$beta)
    weighted_residual <- w[ix] * backend$supports[[sid]]$v * residual
    centered <- if (backend$outcome_recipe$intercept) sweep(d, 2L, f$mean_d) else d
    value[ix, ] <- value[ix, , drop = FALSE] + weighted_residual * (centered %*% adj)
  }
  center <- if (is.null(plan)) as.vector(rec$score_mean_beta %*% adj) else plan$or_center[model_positions]
  value <- sweep(value, 2L, center)
  us <- up[, pos, drop = FALSE]
  if (is.null(plan)) {
    forcing <- t(rec$forcing_beta)
    for (j in seq_len(nt)) value[, j] <- value[, j] +
      as.vector(us %*% (forcing %*% adj[, j, drop = FALSE]))
  } else for (j in seq_len(nt)) value[, j] <- value[, j] +
    as.vector(us %*% plan$or_up[dr_separate_plan_range(plan$or_up_start[model_positions[j]], length(controls))])
  -value
}

dr_separate_baseline_projected_if <- function(backend, reader, rows, links, scale, pre, w, hits, up, plan = NULL) {
  st <- reader$state; B <- length(rows); N <- st$N
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  requested <- backend$links[links]
  if (is.null(plan)) {
    link_models <- vapply(requested, `[[`, integer(1), "model")
    by_model <- split(seq_along(links), factor(link_models, levels = pre$tids))
    names(by_model) <- as.character(pre$tids)
    treatment_blocks <- vapply(backend$treatment_specs[pre$tids], `[[`, integer(1), "regression")
    models_by_block <- split(pre$tids, factor(treatment_blocks, levels = pre$gids))
  } else dr_separate_baseline_replay_check(plan, backend, links, B, pre)
  link_names <- vapply(requested, `[[`, character(1), "link_id")
  phi <- matrix(0, B, length(links), dimnames = list(NULL, link_names))
  for (kb in seq_along(pre$gids)) {
    b <- pre$gids[kb]
    model_positions <- if (is.null(plan)) NULL else
      plan$model_order[dr_separate_plan_range(plan$block_start[kb], plan$block_count[kb])]
    tids <- if (is.null(plan)) models_by_block[[kb]] else pre$tids[model_positions]
    of <- or_outcome_block_fit(backend$outcome_groups, b)
    projected <- dr_separate_baseline_adjoint(backend, b, of, tids, hits, up, pre$sids, w,
      plan, kb, model_positions)
    for (j_model in seq_along(tids)) {
      tid <- tids[j_model]; spec <- backend$treatment_specs[[tid]]; f <- backend$treatments[[tid]]
      ii <- if (is.null(plan)) NULL else dr_separate_plan_range(
        plan$support_start[model_positions[j_model]], plan$support_count[model_positions[j_model]])
      support_pos <- if (is.null(plan)) match(spec$support_ids, pre$sids) else plan$support_pos[ii]
      us <- up[, support_pos, drop = FALSE]
      pt <- backend$treatment_recipe$p
      score <- if (is.null(plan)) matrix(0, B, pt) else numeric(B)
      direct <- matrix(0, B, 2L)
      for (k in seq_along(spec$support_ids)) {
        sid <- spec$support_ids[k]; h <- hits[[support_pos[k]]]; rr <- h$obs
        if (!length(rr)) next
        X <- frame$Z_a[st$a_index[rows[rr]], backend$treatment_recipe$index, drop = FALSE]
        dy <- h$dy; d <- of$L * h$z[, of$keep, drop = FALSE]
        residual <- dy - of$alpha - as.vector(d %*% of$beta)
        fw <- w[rr] * backend$supports[[sid]]$v
        lv <- ipw_logit_values(X, f$coefficients); D <- as.double(k == 1L)
        wp <- fw * D; wm <- fw * (1 - D) * lv$odds
        if (is.null(plan)) score[rr, ] <- fw * (D - lv$prob) * X else
          score[rr] <- fw * (D - lv$prob) * as.vector(X %*% plan$eta[, model_positions[j_model]])
        direct[rr, 1L] <- wp * (residual - f$means["T_dr"]) / f$means["dp"]
        direct[rr, 2L] <- wm * (residual - f$means["M_dr"]) / f$means["dm"]
      }
      if (is.null(plan)) gamma <- t(cp04_solve(f$H_chol, t(score + us %*% t(f$cross))))
      pd <- direct[, 1L] - direct[, 2L]
      po <- projected[, j_model]
      pp <- if (is.null(plan)) -as.vector(gamma %*% f$partials$beta[2L, ]) else
        -score - as.vector(us %*% plan$prop_up[ii])
      ps <- as.vector(us %*% (f$partials$upstream[1L, ] - f$partials$upstream[2L, ]))
      value <- pd + po + pp + ps
      jj <- if (is.null(plan)) by_model[[as.character(tid)]] else plan$link_order[dr_separate_plan_range(
        plan$link_start[model_positions[j_model]], plan$link_count[model_positions[j_model]])]
      for (j in jj) phi[, j] <- value
    }
  }
  cp04_finite(phi, "complete separate DR original-unit IF")
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else N,
    diagnostics = list(outcome_pooling = "separate", n_rate_score_evaluations = length(pre$sids),
      n_selection_score_evaluations = 0L, n_selection_block_evaluations = 0L,
      n_regression_score_evaluations = length(pre$gids), n_global_beta_blocks = 0L,
      n_block_beta_blocks = length(pre$gids), n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    dr_route_abort("baseline", "IF output bound underestimated.", "bound")
  out
}
