# Private separate-outcome DR. Catalogs and treatment nuisance fits belong to
# chained_dr_backend.R; accepted grouped OR helpers own the local coefficient IF.
dr_separate <- function(backend) identical(backend$outcome_pooling, "separate")

# Conservative global reservation precedes dependency maps and chunk buffers.
# Shared selection groups compose their existing upstream envelope afterward.
dr_separate_if_envelope <- function(backend, B, links, details) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  G <- length(backend$regression_specs); S <- length(backend$supports)
  p <- backend$outcome_recipe$p; pt <- backend$treatment_recipe$p
  conditional <- isTRUE(backend$uses_selection_models)
  D <- if (conditional) backend$upstream$D else S
  Tn <- length(backend$treatments); nl <- length(links); J <- length(st$cohorts)
  incidence <- 0
  for (r in backend$regression_specs) incidence <- incidence + length(r$support_ids)
  for (s in backend$treatment_specs) incidence <- incidence + length(s$support_ids)
  if (conditional) {
    for (u in backend$supports) incidence <- incidence + length(u$channels)
    for (m in backend$treatments) incidence <- incidence + length(m$upstream_ids)
  }
  # Only the dependency vectors exist at this admission point. No row buffers
  # or global treatment-by-upstream matrix is constructed by this evaluator.
  maps <- 4096 + 32 * (G + S + D + Tn + nl + incidence)
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    maps, c(backend$accounting$work_ops_bound, 64 * maps), 0, "separate DR dependency metadata")
  gids <- if (details) seq_len(G) else sort(unique(vapply(backend$links[links], `[[`, integer(1), "regression")))
  tids <- sort(unique(vapply(backend$links[links], `[[`, integer(1), "model")))
  sids <- sort(unique(c(vapply(backend$links[links], `[[`, integer(1), "support_id"),
    unlist(lapply(backend$regression_specs[gids], `[[`, "support_ids")),
    unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids")))))
  uids <- if (conditional) sort(unique(c(unlist(lapply(backend$supports[sids], `[[`, "channels")),
    unlist(lapply(backend$treatments[tids], `[[`, "upstream_ids"))))) else sids
  ns <- length(sids); nu <- length(uids); ng <- length(gids); nt <- length(tids)
  ps <- vapply(backend$outcome_groups$groups[gids], function(g) length(g$keep), integer(1))
  model_u <- vapply(tids, function(t) if (conditional) length(backend$treatments[[t]]$upstream_ids) else
    length(backend$treatment_specs[[t]]$support_ids), integer(1))
  max_u <- max(c(0, model_u)); max_p <- max(c(0, ps))
  # Detailed replay visits outcome blocks outside the requested treatment
  # models. Their forcing transpose has local slope-by-channel size.
  operator_peak <- if (conditional) max(c(0, vapply(backend$regressions[gids],
    function(r) as.double(length(r$forcing_beta)), numeric(1)))) else 0
  reference <- chained_backend_reference(backend)
  p_ref <- if (is.null(reference)) 0L else length(reference$x_columns)
  link_bytes <- sum(vapply(backend$links[links], function(l) nchar(l$link_id, type = "bytes"), integer(1)))
  name_bytes <- sum(nchar(c(backend$outcome_recipe$columns, backend$treatment_recipe$columns), type = "bytes"))
  output_meta <- 8192 + 128 * link_bytes + B * frame$id_byte_metadata$max +
    if (details) 2048 * (ng + nt + nl + 1) + 32 * name_bytes * (ng + nt + nl + 1) +
      32 * (ns + nu + Tn + ng) else 128 * nl
  width <- nl + 8 + if (details) 2 * sum(ps + 1) +
    sum(3 * pt + model_u + 40) + nl * (p + 12) + nu + 2 * J + 4 else 0
  output <- output_meta + B * width
  # Hits/upstream rows survive across blocks. OR and treatment temporary rows
  # use only one local solve at a time; detailed records are output above.
  scratch <- maps + output + 2 * operator_peak + 32 * B * (ns * (p + 6) + nu +
    max_p * (max_p + max_u + nl + 8) + 4 * pt + max_u + p_ref + nl + 32) +
    128 * ((max_p + 2)^2 + (pt + 2)^2) + 64 * (ns + nu + ng + nt)
  work <- 128 * maps + 128 * B * ((ns + 1) * (p + 2)^2 *
    (1 + ceiling(log2(st$n_A + 1))) + sum((ps + 2)^2 * (nu + nl + 2)) +
    nt * (pt + 2) * (pt + max_u + 8) + (ns + nt + 1) * p_ref + J + nl * (p + pt + 8))
  if (conditional) {
    sp <- selection_if_budget(backend$selection_backend, B, unique(backend$q_of_s[sids]), details = details)
    scratch <- scratch + sp$output + sp$scratch; work <- work + sp$work
  }
  pre <- list(gids = gids, beta_gids = gids, tids = tids, sids = sids, uids = uids,
    output = output, scratch = scratch, work = work)
  if (conditional && !details) pre <- dr_conditional_separate_ordinary_if_envelope(
    backend, B, links, pre, sp, maps)
  pre
}

evaluate_dr_separate_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE, plan = NULL) {
  validate_dr_backend(backend)
  route <- if (isTRUE(backend$uses_selection_models)) "corollary1" else "baseline"
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  valid <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) > 0L &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    dr_route_abort(route, "Invalid explicit IF row/link batch, scale or details flag.")
  pre <- dr_separate_if_envelope(backend, length(rows), links, details)
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    pre$scratch, c(backend$accounting$work_ops_bound, pre$work), pre$output, "separate DR IF chunk")
  if (!details && !isTRUE(backend$uses_selection_models) && !is.null(backend$separate_complete))
    return(if (is.null(plan)) dr_separate_complete_if(backend, reader, rows, links, scale, pre) else
      separate_complete_native_if(backend, reader, rows, links, scale, pre, plan))
  if (!details && !isTRUE(backend$uses_selection_models) && is.null(backend$reference_ref) &&
      !is.null(plan) && !identical(plan$kind, "dr_separate_baseline_replay_plan"))
    return(dr_separate_sparse_if(backend, reader, rows, links, scale, pre, plan))
  B <- length(rows); N <- st$N; conditional <- isTRUE(backend$uses_selection_models)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  ref <- chained_backend_reference(backend); ci <- st$cohort_index[rows]
  w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
  A <- as.double(st$a_index[rows] > 0L)
  all_or_cols <- backend$or_frame_index[backend$outcome_recipe$slope_index]
  op <- if (conditional) selection_prediction_operation(backend$selection_backend,
    source_reader = reader, source_reference = ref) else NULL
  on.exit(if (!is.null(op)) op$abort(), add = TRUE)
  up_res <- if (conditional) op$upstream(backend$up_map, rows, pre$uids, details) else NULL
  up <- if (conditional) up_res$up else matrix(0, B, length(pre$sids))
  hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]; sp <- backend$support_specs[[sid]]; u <- backend$supports[[sid]]
    ix <- which(ci == sp$cohort)
    hit <- cp04_pair_rows(reader, rows[ix], sp$pair, TRUE)
    obs <- ix[hit$S == 1]; rr <- rows[obs]
    Z <- if (is.null(ref)) frame$Z_a[st$a_index[rr], all_or_cols, drop = FALSE] else
      chained_reference_x(ref, sp$reference, rr)[, all_or_cols, drop = FALSE]
    hits[[j]] <- list(obs = obs, dy = hit$dy[hit$S == 1], z = Z)
    if (conditional) {
      pr <- if (length(rr)) op$pair(rr, backend$q_of_s[sid])$prob else numeric(0)
      a <- if (is.null(ref)) st$a[sp$cohort] else chained_reference_stat(ref, sp$reference, sp$cohort, "a")
      hits[[j]]$v <- or_corollary1_factor(a, st$rho[rr], pr)
    } else if (is.null(ref)) {
      up[ix, j] <- w[ix] * (hit$S - A[ix] * u$rate) / u$massA
    } else {
      ia <- ix[A[ix] > 0]
      wa <- chained_reference_weight(ref, sp$reference, rows[ia])
      up[ia, j] <- wa * (hit$S[A[ix] > 0] - u$rate) / u$massA
    }
    if (!is.null(ref)) hits[[j]]$w <- chained_reference_weight(ref, sp$reference, rr)
  }
  # Ordinary baseline non-reference DR needs only each model's scalar
  # contraction of the frozen outcome coefficient IF, not detailed OR output.
  if (!details && !conditional && is.null(ref))
    return(dr_separate_baseline_projected_if(backend, reader, rows, links, scale, pre, w, hits, up, plan = plan))
  # The pure OR core needs its own target residual mean, not the DR link mean.
  # This invocation-local shallow record is never resealed or retained.
  requested <- backend$links[links]
  link_blocks <- vapply(requested, `[[`, integer(1), "regression")
  link_models <- vapply(requested, `[[`, integer(1), "model")
  by_block <- split(seq_along(links), factor(link_blocks, levels = pre$gids))
  by_model <- split(seq_along(links), factor(link_models, levels = pre$tids))
  names(by_model) <- as.character(pre$tids)
  treatment_blocks <- vapply(backend$treatment_specs[pre$tids], `[[`, integer(1), "regression")
  models_by_block <- split(pre$tids, factor(treatment_blocks, levels = pre$gids))
  or_model <- backend
  or_model$links <- lapply(requested, function(l) { l$mu <- l$mu_or; l })
  link_names <- vapply(requested, `[[`, character(1), "link_id")
  phi <- matrix(0, B, length(links), dimnames = list(NULL, link_names))
  models <- if (details) vector("list", length(backend$treatments)) else NULL
  if (details) {
    blocks <- vector("list", length(pre$gids)); names(blocks) <- as.character(pre$gids)
    den <- num <- or_phi <- direct_or <- phi; design <- vector("list", length(links)); objective <- numeric(B)
    psi_dir <- psi_or <- psi_ps <- psi_sel <- psi_ipw <- phi
  }
  for (kb in seq_along(pre$gids)) {
    b <- pre$gids[kb]; jj <- by_block[[kb]]
    of <- or_outcome_block_fit(backend$outcome_groups, b)
    local_pre <- pre; local_pre$gids <- local_pre$beta_gids <- b
    # One block's temporary coefficient IF at a time, even for ordinary chunks.
    core <- if (conditional) or_selection_separate_if_core(or_model, N, w, hits, up,
      local_pre, jj, TRUE, block_fits = list(of)) else or_baseline_separate_if_core(or_model, N, w, hits, up,
      local_pre, jj, TRUE, block_fits = list(of))
    bf <- core$blocks[[1L]]
    if (details) {
      blocks[[as.character(b)]] <- bf; objective <- objective + core$objective
      if (length(jj)) {
        den[, jj] <- core$denominator; num[, jj] <- core$numerator
        or_phi[, jj] <- core$phi; direct_or[, jj] <- core$direct
        design[jj] <- core$design
      }
    }
    tids <- models_by_block[[kb]]
    for (tid in tids) {
      spec <- backend$treatment_specs[[tid]]; f <- backend$treatments[[tid]]
      us <- up[, match(if (conditional) f$upstream_ids else spec$support_ids,
        if (conditional) pre$uids else pre$sids), drop = FALSE]
      pt <- backend$treatment_recipe$p
      score <- matrix(0, B, pt); direct <- direct_ipw <- matrix(0, B, 4L)
      direct[, 3L] <- direct_ipw[, 3L] <- -f$means["dp"]
      direct[, 4L] <- direct_ipw[, 4L] <- -f$means["dm"]
      for (k in seq_along(spec$support_ids)) {
        sid <- spec$support_ids[k]; h <- hits[[match(sid, pre$sids)]]; rr <- h$obs
        if (!length(rr)) next
        X <- if (is.null(ref)) frame$Z_a[st$a_index[rows[rr]], backend$treatment_recipe$index, drop = FALSE] else
          chained_reference_x(ref, spec$reference, rows[rr])[, backend$treatment_recipe$index, drop = FALSE]
        dy <- h$dy; d <- of$L * h$z[, of$keep, drop = FALSE]
        residual <- dy - of$alpha - as.vector(d %*% of$beta)
        fw <- (if (is.null(h$w)) w[rr] else h$w) *
          if (conditional) h$v else backend$supports[[sid]]$v
        lv <- ipw_logit_values(X, f$coefficients); D <- as.double(k == 1L)
        wp <- fw * D; wm <- fw * (1 - D) * lv$odds
        score[rr, ] <- fw * (D - lv$prob) * X
        direct[rr, 1L] <- wp * (residual - f$means["T_dr"]) / f$means["dp"]
        direct[rr, 2L] <- wm * (residual - f$means["M_dr"]) / f$means["dm"]
        direct[rr, 3L] <- direct[rr, 3L] + wp; direct[rr, 4L] <- direct[rr, 4L] + wm
        direct_ipw[rr, 1L] <- wp * (dy - f$means["T_ipw"]) / f$means["dp"]
        direct_ipw[rr, 2L] <- wm * (dy - f$means["M_ipw"]) / f$means["dm"]
        direct_ipw[rr, 3L] <- direct_ipw[rr, 3L] + wp; direct_ipw[rr, 4L] <- direct_ipw[rr, 4L] + wm
      }
      gamma <- t(cp04_solve(f$H_chol, t(score + us %*% t(f$cross))))
      colnames(gamma) <- backend$treatment_recipe$columns
      pd <- direct[, 1L] - direct[, 2L]
      po <- if (of$p) -as.vector(bf$beta %*% f$d_bar_diff) else numeric(B)
      pp <- -as.vector(gamma %*% f$partials$beta[2L, ])
      ps <- as.vector(us %*% (f$partials$upstream[1L, ] - f$partials$upstream[2L, ]))
      value <- pd + po + pp + ps
      lk <- by_model[[as.character(tid)]]
      for (j in lk) phi[, j] <- value
      if (details) {
        coord <- direct + us %*% t(f$partials$upstream) + gamma %*% t(f$partials$beta)
        coord_ipw <- direct_ipw + us %*% t(f$partials_ipw$upstream) + gamma %*% t(f$partials_ipw$beta)
        colnames(coord) <- colnames(coord_ipw) <- c("T", "M", "dp", "dm")
        full <- coord
        xa <- if (backend$outcome_recipe$intercept) bf$alpha else numeric(B)
        full[, 1L] <- full[, 1L] - xa - as.vector(bf$beta %*% f$d_bar_plus)
        full[, 2L] <- full[, 2L] - xa - as.vector(bf$beta %*% f$d_bar_minus)
        models[[tid]] <- list(beta = gamma, gamma = gamma, coordinates = coord,
          coordinates_full = full, direct = direct, beta_score = score, upstream = us,
          psi_dr = value, psi_dir = pd, psi_or = po, psi_ps = pp, psi_sel = ps,
          psi_ipw = coord_ipw[, 1L] - coord_ipw[, 2L], normalizers_A =
            coord[, 3:4, drop = FALSE] / st$lambda - outer(A - st$lambda, f$means[c("dp", "dm")] / st$lambda^2))
        for (j in lk) {
          psi_dir[, j] <- pd; psi_or[, j] <- po; psi_ps[, j] <- pp; psi_sel[, j] <- ps
          psi_ipw[, j] <- coord_ipw[, 1L] - coord_ipw[, 2L]
        }
      }
    }
    rm(core, bf)
  }
  cp04_finite(phi, "complete separate DR original-unit IF")
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else N,
    diagnostics = list(outcome_pooling = "separate", n_rate_score_evaluations = if (conditional)
      sum(pre$uids <= backend$upstream$n_a) else length(pre$sids),
      n_selection_score_evaluations = if (conditional) length(up_res$active_cells) else 0L,
      n_selection_block_evaluations = if (conditional) length(up_res$active_blocks) else 0L,
      n_regression_score_evaluations = length(pre$gids), n_global_beta_blocks = 0L,
      n_block_beta_blocks = length(pre$gids), n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (details) {
    tau <- eta <- if (is.null(ref)) matrix(0, B, length(st$cohorts)) else NULL
    if (is.null(ref)) for (k in seq_along(st$cohorts)) {
      tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
      eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
    }
    out$details <- list(outcome_pooling = "separate", blocks = blocks,
      link_block_ids = vapply(backend$links[links], `[[`, integer(1), "regression"),
      support_ids = pre$sids, upstream_ids = pre$uids, upstream = if (conditional) up else NULL,
      rates = if (!conditional) up else NULL, regression_ids = pre$gids,
      direct = direct_or, denominator = den, numerator = num, design = design, objective = objective,
      treatment_ids = pre$tids, models = models, psi_dir = psi_dir, psi_or = psi_or,
      psi_ps = psi_ps, psi_sel = psi_sel, psi_ipw = psi_ipw, psi_or_estimator = or_phi,
      tau = tau, eta = eta, lambda = A - st$lambda)
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    dr_route_abort(route, "IF output bound underestimated.", "bound")
  if (!is.null(op)) op$finish()
  out
}
