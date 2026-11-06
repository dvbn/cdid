# Ordinary complete separate replay only. These are the scalar complete
# equations, with frozen scalar adjoints and a bounded native row loop. No
# level regression, shared outcome slope, moment reduction or GMM change.
separate_complete_native_live <- function(B, nl, nd, ng, nf) {
  # Fresh phi, endpoint doubles + byte states (charged as doubles), four
  # per-block row scalars, two family scalars, return headers and status.
  1024 + B * (nl + 2 * nd) + 4 * ng + 2 * nf
}

separate_complete_native_plan <- function(backend, links, B) {
  dr <- inherits(backend, "cdid_dr_baseline")
  if (dr) validate_dr_baseline_backend(backend) else validate_or_baseline_backend(backend)
  if (!identical(backend$outcome_pooling, "separate") || is.null(backend$separate_complete) ||
      isTRUE(backend$uses_selection_models) || !is.null(backend$reference_ref)) return(NULL)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) ||
      B < 1 || B != floor(B) || B > backend$control$max_chunk_rows)
    or_baseline_abort("Invalid complete separate replay link batch or chunk rows.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (B > st$N || !is.double(frame$Z_a)) return(NULL)
  pre <- if (dr) dr_baseline_if_budget(backend, B, links, FALSE) else or_baseline_if_budget(backend, B, links, FALSE)
  if (any(vapply(pre[c("sids", "gids", if (dr) "tids")], function(x) !is.null(names(x)), logical(1)))) return(NULL)
  fixed <- backend$separate_complete; ng <- length(pre$gids); ns <- length(pre$sids); nl <- length(links)
  p <- unname(vapply(backend$outcome_groups$groups[pre$gids], function(g) length(g$keep), integer(1)))
  control_n <- unname(vapply(backend$regression_specs[pre$gids], function(s) length(s$support_ids), integer(1)))
  mb <- if (dr) match(fixed$treatment_block[pre$tids], pre$gids) else match(fixed$link_block[links], pre$gids)
  if (anyNA(mb)) or_baseline_abort("Complete replay model lost its outcome block.", "integrity")
  nm <- length(mb); pt <- if (dr) backend$treatment_recipe$p else 0L
  family_ids <- if (dr) unique(fixed$treatment_family[pre$tids]) else integer(0)
  nf <- length(family_ids)
  family_first <- if (dr) pre$tids[match(family_ids, fixed$treatment_family[pre$tids])] else integer(0)
  risk_n <- if (dr) unname(vapply(backend$treatment_specs[family_first], function(s) length(s$cohorts), integer(1))) else integer(0)
  P <- sum(p); A <- sum(p[mb]); C <- sum(control_n); K <- sum(risk_n)
  if (any(!is.finite(c(P, A, C, K, pt * nm, pt * nf))) ||
      max(c(P, A, C, K, pt * nm, pt * nf, ng, nm, nl, nf)) > .Machine$integer.max) return(NULL)
  nrhs <- tabulate(mb, nbins = ng)
  peak_rhs <- max(c(0, p * nrhs)); mp <- max(c(0, p))
  # Flat arrays below, their SEXP/list headers and construction replacement
  # copies. Borrowed backend/frame/panel payload is already owned upstream.
  flat <- 6 * ng + 4 * ng + 3 * P + C + 4 * nm + 7 * nm + A + pt * nm +
    2 * nf + K + pt + pt * nf + nl + length(fixed$dates)
  retained <- 2048 + 4 * flat + 8 * (ns + ng + nm + nl + nf)
  peak <- 1024 + 16 * peak_rhs + 8 * mp^2 + 16 * pt^2 + 32 * pt +
    4 * max(c(0, A, P, C, K, pt * nm, pt * nf)) + 8 * (ng + nm + nl + nf)
  tail <- st$N %% B; min_B <- if (tail) min(B, tail) else B
  # The old envelopes reserve 32*B*(ns*(p_recipe+5)+nl), or more for DR.
  # This route replaces hits and their block row matrices. At most half of
  # that admitted slice funds the plan/construction; at least half remains
  # for fresh native output/gather/small row scratch and scaling. All other
  # old terms remain. No cap, B or whole-consumer admission is changed.
  slice <- 16 * min_B * (ns * (backend$outcome_recipe$p + 5) + nl)
  bound <- retained + peak
  live <- separate_complete_native_live(B, nl, length(fixed$dates), ng, nf)
  if (!is.finite(bound) || bound > slice ||
      live > 16 * B * (ns * (backend$outcome_recipe$p + 5) + nl)) return(NULL)
  starts <- function(n) cumsum(c(1, utils::head(n, -1L)))
  range <- function(s, n) if (n) seq.int(s, length.out = n) else integer(0)
  ps <- starts(p); cs <- starts(control_n); ads <- starts(p[mb]); fs <- starts(risk_n)
  # Native descriptor order is fixed and independently checked at the C boundary.
  data <- list(mode = as.integer(dr), block = matrix(0L, 6L, ng), block_par = matrix(0, 4L, ng),
    keep = integer(P), beta = numeric(P), mean_d = numeric(P), controls = integer(C),
    model = matrix(0L, 4L, nm), model_par = matrix(0, 7L, nm), adj = numeric(A),
    eta = numeric(pt * nm), family = matrix(0L, 2L, nf), risk = integer(K),
    trt_columns = if (dr) as.integer(backend$treatment_recipe$index) else integer(0),
    trt_beta = numeric(pt * nf), link_model = if (dr) as.integer(match(fixed$link_model[links], pre$tids)) else seq_len(nl))
  outcome_cols <- if (dr) backend$or_frame_index else seq_along(backend$outcome_recipe$columns)
  for (kb in seq_len(ng)) {
    b <- pre$gids[kb]; of <- or_outcome_block_fit(backend$outcome_groups, b)
    sids <- backend$regression_specs[[b]]$support_ids
    pairs <- vapply(backend$support_specs[sids], `[[`, integer(1), "pair_id")
    cohorts <- fixed$support_cohort[sids]
    if (length(unique(pairs)) != 1L || anyDuplicated(cohorts)) return(NULL)
    pair <- backend$support_specs[[sids[1L]]]$pair
    data$block[, kb] <- as.integer(c(match(pair$s_calendar_index, fixed$dates),
      match(pair$t_calendar_index, fixed$dates), ps[kb], p[kb], cs[kb], control_n[kb]))
    rec <- backend$regressions[[b]]
    data$block_par[, kb] <- c(of$L, of$alpha, rec$mass, rec$score_mean_alpha)
    ix <- range(ps[kb], of$p)
    data$keep[ix] <- outcome_cols[backend$outcome_recipe$slope_index[of$keep]]
    data$beta[ix] <- of$beta; data$mean_d[ix] <- of$mean_d
    data$controls[range(cs[kb], control_n[kb])] <- cohorts
    models <- which(mb == kb)
    loads <- matrix(0, of$p, length(models))
    if (of$p) for (j in seq_along(models)) {
      k <- models[j]
      loads[, j] <- if (dr) backend$treatments[[pre$tids[k]]]$d_bar_diff else
        backend$links[[links[k]]]$c - of$mean_d
    }
    adj <- if (of$p) or_pooled_wls_gram_solve(of$fit, loads,
      normalizer = backend$outcome_groups$N) else loads
    for (j in seq_along(models)) {
      k <- models[j]
      data$model[, k] <- as.integer(c(kb, ads[k], if (dr) match(fixed$treatment_family[pre$tids[k]], family_ids) else 0L,
        if (dr) backend$treatment_specs[[pre$tids[k]]]$cohorts[1L] else fixed$support_cohort[backend$links[[links[k]]]$support_id]))
      data$adj[range(ads[k], of$p)] <- adj[, j]
      data$model_par[1L, k] <- if (of$p) as.vector(rec$score_mean_beta %*% adj[, j, drop = FALSE]) else 0
      if (dr) {
        spec <- backend$treatment_specs[[pre$tids[k]]]
        ss <- backend$support_specs[spec$support_ids]
        if (any(vapply(ss, `[[`, integer(1), "pair_id") != pairs[1L]) ||
            anyDuplicated(vapply(ss, `[[`, integer(1), "cohort")) ||
            !identical(vapply(ss, `[[`, integer(1), "cohort"), spec$cohorts) ||
            !setequal(spec$cohorts[-1L], cohorts)) return(NULL)
        f <- backend$treatments[[pre$tids[k]]]
        data$model_par[4:7, k] <- f$means[c("T_dr", "M_dr", "dp", "dm")]
        data$eta[range(1 + (k - 1) * pt, pt)] <- as.vector(cp04_solve(f$H_chol, f$partials$beta[2L, ]))
      } else {
        l <- backend$links[[links[k]]]
        if (backend$support_specs[[l$support_id]]$pair_id != pairs[1L]) return(NULL)
        data$model_par[2:3, k] <- c(l$mu, l$d)
      }
    }
  }
  for (j in seq_len(nf)) {
    tid <- family_first[j]; first <- backend$treatment_specs[[tid]]; f <- backend$treatments[[tid]]
    tids <- pre$tids[fixed$treatment_family[pre$tids] == family_ids[j]]
    for (t in tids) if (!identical(backend$treatment_specs[[t]]$cohorts, first$cohorts) ||
        !identical(backend$treatments[[t]]$coefficients, f$coefficients) ||
        !identical(backend$treatments[[t]]$H_chol, f$H_chol))
      dr_baseline_abort("Separate complete propensity replay differs from its fitted operator.", "identity")
    data$family[, j] <- as.integer(c(fs[j], risk_n[j]))
    data$risk[range(fs[j], risk_n[j])] <- first$cohorts
    data$trt_beta[range(1 + (j - 1) * pt, pt)] <- f$coefficients
  }
  plan <- list(kind = "separate_complete_native_plan", backend_ref = backend$ref,
    links = as.integer(links), B = B, min_B = min_B, N = st$N,
    sids = pre$sids, gids = pre$gids, tids = pre$tids,
    dates = fixed$dates, data = data, retained_bound = retained, construction_peak = peak,
    reservation = bound, entries = 0)
  plan$entries <- unname(cp04_payload(plan)["entries"]) + 64
  if (plan$entries > retained) or_baseline_abort("Complete separate replay plan bound underestimated.", "bound")
  ipw_owner_seal(plan, "separate_complete_native_plan")
}

separate_complete_native_check <- function(plan, backend, links, B, pre) {
  or_baseline_native_check(plan, backend, links, B, pre, "separate_complete_native_plan", or_baseline_abort)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  if (!isTRUE(plan$N == reader$state$N) || !isTRUE(B == plan$B || B == plan$min_B) ||
      !identical(plan$dates, backend$separate_complete$dates) ||
      !isTRUE(plan$reservation <= 16 * B * (length(pre$sids) * (backend$outcome_recipe$p + 5) + length(links))))
    or_baseline_abort("Complete separate replay plan differs from its issued partition or source dates.", "integrity")
  invisible(plan)
}

separate_complete_native_if <- function(backend, reader, rows, links, scale, pre, plan) {
  B <- length(rows); st <- reader$state; dr <- inherits(backend, "cdid_dr_baseline")
  separate_complete_native_check(plan, backend, links, B, pre)
  ng <- length(pre$gids); nf <- ncol(plan$data$family)
  live <- separate_complete_native_live(B, length(links), length(plan$dates), ng, nf)
  ipw_owner_safe_sum(c(plan$entries, live), pre$scratch, "complete separate native replay")
  po <- or_baseline_complete_panel(reader, backend$separate_complete); post <- po$postings
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  res <- or_baseline_native_kernel("C_gmm_separate_complete_if")(po$observations$y,
    post$date_start, post$post_a_index, post$post_obs_index, as.integer(rows), st$cohort_index,
    st$a_index, st$w, plan$dates, frame$Z_a, plan$data, as.double(live))
  if (!identical(res$status[1L], 0L)) or_baseline_native_refuse(res$status,
    if (dr) dr_baseline_abort else or_baseline_abort)
  out <- or_separate_complete_output(backend, st, rows, links, scale, pre, res$phi,
    ng, if (dr) length(pre$tids) else 0L, nf)
  out$diagnostics$execution <- "native_fixed_operator"
  out
}
