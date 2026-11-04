# Optional static incomplete pooled baseline DR. Sparse preparation and four-
# coordinate assembly are native; the original fitted-factor solves remain R.
# No selected/reference/separate/complete route or scientific kernel changes.
dr_sparse_native_inventory <- function(B, ns, nd, nt, nl, p, pt, metadata) {
  prepared <- 1024 + B * (ns + p + pt * nt + 4 * nt)
  prepare <- prepared + 3 * B * nd + p + pt
  # All prepared RHSs and direct records stay live with every full xi_gamma.
  # Count one model's RHS slice/transposes/solve results and OR solve copies.
  solve <- prepared + B * (pt * nt + 6 * pt + 8 * p) + 8 * (p^2 + pt^2) + 128 * nt + 4096
  assemble <- prepared + B * (pt * nt + p + nl + 1) + nl + 128 * nt + 2048
  output <- prepared + B * (pt * nt + p + 3 * nl) + 128 * nt + 4096
  list(prepare = prepare, peak = max(prepare, solve, assemble, output) + metadata,
    assembler = 1024 + B * (nl + 1) + nl)
}

dr_sparse_native_plan <- function(backend, links, B) {
  validate_dr_baseline_backend(backend)
  if (dr_separate(backend) || isTRUE(backend$uses_selection_models) || !is.null(backend$reference_ref) ||
      !is.null(backend$complete_support) || !is.null(backend$separate_complete)) return(NULL)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) || B < 1 ||
      B != floor(B) || B > backend$control$max_chunk_rows)
    dr_baseline_abort("Invalid sparse pooled DR replay link batch or chunk rows.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  post <- panel$postings
  if (isTRUE(st$reference_carrier) || B > st$N || st$N > .Machine$integer.max ||
      !is.double(frame$Z_a) || !is.double(panel$observations$y) ||
      !is.integer(st$a_index) || !is.integer(st$cohort_index) ||
      !is.integer(post$post_a_index) || is.object(post$post_a_index) ||
      !is.integer(post$post_obs_index) || is.object(post$post_obs_index) ||
      !is.integer(post$date_start) || !is.integer(post$date_end) ||
      length(post$post_a_index) < 1L || length(post$post_a_index) > .Machine$integer.max ||
      length(post$post_a_index) != length(post$post_obs_index)) return(NULL)
  pre <- dr_baseline_if_budget(backend, B, links, FALSE)
  if (any(vapply(pre[c("sids", "gids", "beta_gids", "tids")], function(x) !is.null(names(x)), logical(1)))) return(NULL)
  p <- backend$shared$p; pt <- backend$treatment_recipe$p
  ns <- length(pre$sids); nb <- length(pre$beta_gids); nt <- length(pre$tids); nl <- length(links)
  G <- length(backend$regression_specs); J <- length(st$cohorts)
  bc <- unname(vapply(backend$regression_specs[pre$beta_gids], function(s) length(s$support_ids), integer(1)))
  mc <- unname(vapply(backend$treatment_specs[pre$tids], function(s) length(s$support_ids), integer(1)))
  I <- sum(bc); M <- sum(mc); F <- if (p) length(backend$shared$control_support_ids) else 0L
  model <- match(vapply(backend$links[links], `[[`, integer(1), "model"), pre$tids)
  if (anyNA(model)) dr_baseline_abort("Sparse DR link lost its actual treatment model.", "integrity")
  dates <- sort(unique(as.integer(unlist(lapply(backend$support_specs[pre$sids], function(s)
    c(s$pair$s_calendar_index, s$pair$t_calendar_index)), use.names = FALSE))))
  flat <- 1 + 6 * ns + 3 * nb + (2 + p) * nb + I + 2 * p + F + p * F + pt +
    4 * nt + 6 * nt + M + pt * nt + pt * M + 4 * pt * nt + 4 * M + p * nt + nl +
    ns + length(pre$gids) + nb + nt + nl + 3 * length(dates)
  retained <- 4096 + 4 * flat + 16 * (ns + nb + nt + nl)
  construction <- 2048 + 8 * max(c(0, I, M, p * F, pt * M, 4 * pt * nt, 4 * M, p * nt)) +
    16 * (ns + nb + nt + nl)
  # New fixed maps have no replicated source labels. Existing preflight sets,
  # row identity/label copies and small list/status records remain charged.
  metadata <- 8192 + 64 * (ns + G + nt + nl + J) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  tail <- st$N %% B; min_B <- if (tail) min(B, tail) else B
  small <- if (min_B == B) pre else dr_baseline_if_budget(backend, min_B, links, FALSE)
  one <- dr_sparse_native_inventory(B, ns, length(dates), nt, nl, p, pt, metadata)
  last <- dr_sparse_native_inventory(min_B, ns, length(dates), nt, nl, p, pt, metadata)
  reservation <- retained + construction
  if (any(!is.finite(c(flat, reservation, one$peak, last$peak))) ||
      max(c(I, M, F, p * F, pt * nt, pt * M, 4 * pt * nt, 4 * M, p * nt, ns, nb, nt, nl)) >= .Machine$integer.max ||
      reservation + one$peak > pre$scratch || reservation + last$peak > small$scratch ||
      64 * flat > small$work / 4) return(NULL)
  start <- function(n) as.integer(cumsum(c(1, utils::head(n, -1L))))
  span <- function(s, n) if (n) seq.int(s, length.out = n) else integer(0)
  bs <- start(bc); ms <- start(mc); ln <- tabulate(model, nbins = nt); ls <- start(ln)
  data <- list(intercept = as.integer(backend$outcome_recipe$intercept),
    support = matrix(0L, 3L, ns), support_par = matrix(0, 3L, ns), block = matrix(0L, 3L, nb),
    block_par = matrix(0, 2L + p, nb), controls = integer(I),
    or_keep = as.integer(backend$or_frame_index[backend$outcome_recipe$slope_index[backend$shared$keep]]),
    or_beta = as.double(backend$shared$beta), forcing_support = integer(F),
    forcing_beta = numeric(p * F), trt_columns = as.integer(backend$treatment_recipe$index),
    model = matrix(0L, 4L, nt), model_par = matrix(0, 6L, nt), model_support = integer(M),
    trt_beta = matrix(0, pt, nt), cross = numeric(pt * M),
    partial_beta = matrix(0, 4L * pt, nt), partial_up = numeric(4 * M),
    d_bar_diff = matrix(0, p, nt), link_order = as.integer(order(model, seq_along(model))))
  for (j in seq_len(ns)) {
    s <- backend$support_specs[[pre$sids[j]]]; u <- backend$supports[[pre$sids[j]]]
    data$support[, j] <- as.integer(c(s$cohort, match(s$pair$s_calendar_index, dates), match(s$pair$t_calendar_index, dates)))
    data$support_par[, j] <- c(u$rate, u$massA, u$v)
  }
  for (j in seq_len(nb)) {
    b <- pre$beta_gids[j]; r <- backend$regressions[[b]]
    data$block[, j] <- c(bs[j], bc[j], as.integer(b))
    data$block_par[, j] <- c(r$L, if (backend$outcome_recipe$intercept) r$alpha else 0, r$mean_d)
    data$controls[span(bs[j], bc[j])] <- match(backend$regression_specs[[b]]$support_ids, pre$sids)
  }
  if (p) {
    data$forcing_support <- as.integer(match(backend$shared$control_support_ids, pre$sids))
    data$forcing_beta <- as.double(backend$shared$forcing_beta[, backend$shared$control_support_ids, drop = FALSE])
  }
  for (j in seq_len(nt)) {
    t <- pre$tids[j]; spec <- backend$treatment_specs[[t]]; f <- backend$treatments[[t]]
    r <- backend$regressions[[spec$regression]]; m <- mc[j]
    if (!identical(vapply(backend$support_specs[spec$support_ids], `[[`, integer(1), "cohort"), spec$cohorts) ||
        !identical(dim(f$cross), c(as.integer(pt), as.integer(m))) ||
        !identical(dim(f$partials$upstream), c(4L, as.integer(m))) ||
        !identical(dim(f$partials$beta), c(4L, as.integer(pt))))
      dr_baseline_abort("Sparse DR descriptors differ from the frozen model's risk/partial coordinates.", "integrity")
    data$model[, j] <- c(ms[j], mc[j], ls[j], ln[j])
    data$model_par[, j] <- c(r$L, if (backend$outcome_recipe$intercept) r$alpha else 0,
      f$means[c("T_dr", "M_dr", "dp", "dm")])
    data$model_support[span(ms[j], m)] <- match(spec$support_ids, pre$sids)
    data$trt_beta[, j] <- f$coefficients
    data$cross[span(1 + pt * (ms[j] - 1), pt * m)] <- as.double(f$cross)
    data$partial_beta[, j] <- as.double(f$partials$beta)
    data$partial_up[span(1 + 4 * (ms[j] - 1), 4 * m)] <- as.double(f$partials$upstream)
    if (p) data$d_bar_diff[, j] <- f$d_bar_diff
  }
  bounds <- rbind(panel$postings$date_start[dates], panel$postings$date_end[dates])
  bounds[, panel$postings$date_count[dates] == 0L] <- 0L; storage.mode(bounds) <- "integer"
  if (anyNA(data$controls) || anyNA(data$forcing_support) || anyNA(data$model_support))
    dr_baseline_abort("Sparse DR descriptors lost their source supports.", "integrity")
  plan <- list(kind = "dr_sparse_native_plan", backend_ref = backend$ref, links = as.integer(links),
    B = B, min_B = min_B, N = st$N, sids = pre$sids, gids = pre$gids, beta_gids = pre$beta_gids,
    tids = pre$tids, dates = dates, bounds = bounds, data = data, metadata = metadata,
    retained_bound = retained, construction_peak = construction, reservation = reservation, entries = 0)
  plan$entries <- unname(cp04_payload(plan)["entries"]) + 64
  if (plan$entries > retained) dr_baseline_abort("Sparse DR plan bound underestimated.", "bound")
  ipw_owner_seal(plan, "dr_sparse_native_plan")
}

dr_sparse_native_check <- function(plan, backend, links, B, pre) {
  or_baseline_native_check(plan, backend, links, B, pre, "dr_sparse_native_plan", dr_baseline_abort)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  if (!identical(plan$beta_gids, pre$beta_gids) || !isTRUE(plan$N == reader$state$N) ||
      !isTRUE(B == plan$B || B == plan$min_B))
    dr_baseline_abort("Sparse DR replay differs from its issued partition or global score blocks.", "integrity")
  one <- dr_sparse_native_inventory(B, length(pre$sids), length(plan$dates), length(pre$tids), length(links),
    backend$shared$p, backend$treatment_recipe$p, plan$metadata)
  ipw_owner_safe_sum(c(plan$reservation, one$peak), pre$scratch, "sparse DR fixed-operator replay")
  invisible(plan)
}

# A private replay status retains the scalar treatment-replay numeric error;
# existing complete-support status9 continues to mean BLAS identity refusal.
dr_sparse_native_refuse <- function(status) {
  if (identical(status[1L], 27L)) cp04_abort("Nonfinite baseline DR treatment replay.", "numeric")
  or_baseline_native_refuse(status, dr_baseline_abort)
}

dr_sparse_native_if <- function(backend, reader, rows, links, scale, pre, plan) {
  dr_sparse_native_check(plan, backend, links, length(rows), pre)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE); B <- length(rows)
  p <- backend$shared$p; pt <- backend$treatment_recipe$p; nt <- length(pre$tids)
  inventory <- dr_sparse_native_inventory(B, length(pre$sids), length(plan$dates), nt, length(links), p, pt, plan$metadata)
  a <- .Call(or_baseline_native_kernel("C_gmm_dr_sparse_prepare"), panel$observations$y,
    panel$postings$post_a_index, panel$postings$post_obs_index, as.integer(rows), st$cohort_index,
    st$a_index, st$w, plan$bounds, frame$Z_a, plan$data, as.double(inventory$prepare))
  code <- a$status
  if (!identical(code[1L], 0L) && identical(code[2L], 1L)) dr_sparse_native_refuse(code)
  # Keep the original normalized global OR solve and its finite-input refusal.
  xi_beta <- if (p) t(if (backend$outcome_recipe$intercept)
    or_pooled_wls_gram_solve(backend$shared$fit, t(a$or_score), normalizer = st$N) else
    or_wls_gram_solve(backend$shared$fit, t(a$or_score), normalizer = st$N)) else matrix(0, B, 0L)
  n_ok <- if (identical(code[1L], 0L)) nt else code[3L] - 1L
  xi <- vector("list", n_ok)
  for (j in seq_len(n_ok)) {
    ix <- (j - 1L) * pt + seq_len(pt)
    xi[[j]] <- t(cp04_solve(backend$treatments[[pre$tids[j]]]$H_chol, t(a$rhs[, ix, drop = FALSE])))
  }
  # A later model's read/replay failure is deferred until all earlier model
  # gamma/four-coordinate checks have run, matching the complete staged path.
  b <- .Call(or_baseline_native_kernel("C_gmm_dr_sparse_assemble"), a$up, xi_beta, a$direct, xi,
    plan$data, as.integer(n_ok), as.double(inventory$assembler))
  if (!identical(b$status[1L], 0L)) dr_sparse_native_refuse(b$status)
  if (!identical(code[1L], 0L)) dr_sparse_native_refuse(code)
  phi <- b$phi
  if (scale == "core") phi <- st$lambda * phi
  dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id"))
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_rate_score_evaluations = length(pre$sids),
      n_regression_score_evaluations = if (p) length(union(pre$gids, pre$beta_gids)) else 0L,
      n_global_beta_blocks = if (p) length(pre$beta_gids) else 0L, n_treatment_score_evaluations = nt,
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work, execution = "native_sparse_pooled_dr"))
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    dr_baseline_abort("Sparse DR IF output bound underestimated.", "bound")
  out
}
