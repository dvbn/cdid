# Optional native companion to the unchanged V11 projected separate-DR plan.
# It keeps V11's actual final-phi numerical contract, not pooled DR's stronger
# gamma/coordinate/replay checks. Fits, buffers, rows and stores are unchanged.
dr_separate_sparse_live <- function(B, ns, nd, nt, nl) 1024 + B * (3 * nd + ns + nt + nl) + nt + nl

dr_separate_sparse_inventory <- function(B, ns, nd, nt, nl, metadata) {
  # Keep C's date/rate/scalar-model/output buffers plus all R phi/finite-check/
  # scale copies conservatively concurrent. Preflight sets/row labels survive.
  dr_separate_sparse_live(B, ns, nd, nt, nl) + 3 * B * nl + metadata
}

dr_separate_sparse_plan <- function(backend, links, B, v11) {
  if (is.null(v11)) return(NULL)
  validate_dr_baseline_backend(backend)
  if (!dr_separate(backend) || isTRUE(backend$uses_selection_models) || !is.null(backend$reference_ref) ||
      !is.null(backend$separate_complete)) return(v11)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE); post <- panel$postings
  pre <- dr_baseline_if_budget(backend, B, links, FALSE)
  dr_separate_baseline_replay_check(v11, backend, links, B, pre)
  if (isTRUE(st$reference_carrier) || st$N > .Machine$integer.max || !is.double(frame$Z_a) ||
      !is.double(panel$observations$y) || !is.integer(st$a_index) || !is.integer(st$cohort_index) ||
      !is.integer(post$post_a_index) || is.object(post$post_a_index) ||
      !is.integer(post$post_obs_index) || is.object(post$post_obs_index) ||
      !is.integer(post$date_start) || !is.integer(post$date_end) ||
      length(post$post_a_index) < 1L || length(post$post_a_index) > .Machine$integer.max ||
      length(post$post_a_index) != length(post$post_obs_index)) return(v11)
  ns <- length(pre$sids); ng <- length(pre$gids); nt <- length(pre$tids); nl <- length(links)
  pt <- backend$treatment_recipe$p; p <- v11$p; P <- sum(p); M <- sum(v11$support_count)
  dates <- sort(unique(as.integer(unlist(lapply(backend$support_specs[pre$sids], function(s)
    c(s$pair$s_calendar_index, s$pair$t_calendar_index)), use.names = FALSE))))
  flat <- 6 * ns + 4 * ng + 3 * P + nt + pt * nt + 4 * nt + M + pt +
    ns + length(pre$uids) + ng + nt + nl + 3 * length(dates)
  retained <- 4096 + 4 * flat + 16 * (ns + ng + nt + nl)
  peak <- 2048 + 8 * max(c(0, P, M, pt * nt, 4 * nt)) + 16 * (ns + ng + nt + nl)
  metadata <- 8192 + 64 * (ns + ng + nt + nl) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  reservation <- v11$reservation + retained + peak
  small <- if (v11$min_B == B) pre else dr_baseline_if_budget(backend, v11$min_B, links, FALSE)
  if (any(!is.finite(c(flat, reservation))) || max(c(P, M, pt * nt, ng, ns, nt, nl)) >= .Machine$integer.max ||
      reservation + dr_separate_sparse_inventory(B, ns, length(dates), nt, nl, metadata) > pre$scratch ||
      reservation + dr_separate_sparse_inventory(v11$min_B, ns, length(dates), nt, nl, metadata) > small$scratch ||
      64 * flat > small$work / 4) return(v11)
  start <- as.integer(cumsum(c(1, utils::head(p, -1L))))
  span <- function(s, n) if (n) seq.int(s, length.out = n) else integer(0)
  data <- list(intercept = as.integer(backend$outcome_recipe$intercept),
    support = matrix(0L, 3L, ns), support_par = matrix(0, 3L, ns),
    block = rbind(start, as.integer(p)), block_par = matrix(0, 2L, ng),
    keep = integer(P), beta = numeric(P), mean_d = numeric(P), model_block = integer(nt),
    trt_beta = matrix(0, pt, nt), means = matrix(0, 4L, nt), selection = numeric(M),
    trt_columns = as.integer(backend$treatment_recipe$index))
  storage.mode(data$block) <- "integer"
  for (j in seq_len(ns)) {
    s <- backend$support_specs[[pre$sids[j]]]; u <- backend$supports[[pre$sids[j]]]
    data$support[, j] <- as.integer(c(s$cohort, match(s$pair$s_calendar_index, dates), match(s$pair$t_calendar_index, dates)))
    data$support_par[, j] <- c(u$rate, u$massA, u$v)
  }
  for (k in seq_len(ng)) {
    of <- or_outcome_block_fit(backend$outcome_groups, pre$gids[k]); ii <- span(start[k], of$p)
    data$block_par[, k] <- c(of$L, of$alpha)
    data$keep[ii] <- backend$or_frame_index[backend$outcome_recipe$slope_index[of$keep]]
    data$beta[ii] <- of$beta; data$mean_d[ii] <- of$mean_d
  }
  for (j in seq_len(nt)) {
    spec <- backend$treatment_specs[[pre$tids[j]]]; f <- backend$treatments[[pre$tids[j]]]
    data$model_block[j] <- match(spec$regression, pre$gids)
    data$trt_beta[, j] <- f$coefficients
    data$means[, j] <- f$means[c("T_dr", "M_dr", "dp", "dm")]
    data$selection[span(v11$support_start[j], v11$support_count[j])] <-
      f$partials$upstream[1L, ] - f$partials$upstream[2L, ]
  }
  if (anyNA(data$model_block)) dr_baseline_abort("Sparse projected model lost its local outcome block.", "integrity")
  bounds <- rbind(post$date_start[dates], post$date_end[dates])
  bounds[, post$date_count[dates] == 0L] <- 0L; storage.mode(bounds) <- "integer"
  # The original V11 plan is retained ONCE as a nested borrowed record; its
  # arrays are not copied into data. Its full reservation remains charged.
  plan <- list(kind = "dr_separate_sparse_plan", backend_ref = backend$ref,
    links = as.integer(links), B = B, min_B = v11$min_B, N = v11$N,
    sids = pre$sids, uids = pre$uids, gids = pre$gids, tids = pre$tids,
    dates = dates, bounds = bounds, data = data, v11 = v11,
    metadata = metadata, retained_bound = retained + v11$retained_bound,
    construction_peak = peak + v11$construction_peak, reservation = reservation, entries = 0)
  plan$entries <- unname(cp04_payload(plan)["entries"]) + 64
  if (plan$entries > plan$retained_bound) dr_baseline_abort("Sparse separate DR plan bound underestimated.", "bound")
  ipw_owner_seal(plan, "dr_separate_sparse_plan")
}

dr_separate_sparse_replay_plan <- function(backend, links, B) {
  v11 <- dr_separate_baseline_replay_plan(backend, links, B)
  dr_separate_sparse_plan(backend, links, B, v11)
}

dr_separate_sparse_check <- function(plan, backend, links, B, pre) {
  or_baseline_native_check(plan, backend, links, B, pre, "dr_separate_sparse_plan", dr_baseline_abort)
  dr_separate_baseline_replay_check(plan$v11, backend, links, B, pre)
  if (!identical(plan$uids, pre$uids) || !isTRUE(plan$N == plan$v11$N) || !isTRUE(plan$min_B == plan$v11$min_B))
    dr_baseline_abort("Sparse separate DR replay differs from its issued upstream/row partition.", "integrity")
  ipw_owner_safe_sum(c(plan$reservation, dr_separate_sparse_inventory(B, length(pre$sids), length(plan$dates),
    length(pre$tids), length(links), plan$metadata)), pre$scratch, "sparse separate DR full replay lifetime")
  invisible(plan)
}

dr_separate_sparse_if <- function(backend, reader, rows, links, scale, pre, plan) {
  B <- length(rows); st <- reader$state
  dr_separate_sparse_check(plan, backend, links, B, pre)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  live <- dr_separate_sparse_live(B, length(pre$sids), length(plan$dates), length(pre$tids), length(links))
  res <- .Call(or_baseline_native_kernel("C_gmm_dr_separate_sparse_if"), panel$observations$y,
    panel$postings$post_a_index, panel$postings$post_obs_index, as.integer(rows), st$cohort_index,
    st$a_index, st$w, plan$bounds, frame$Z_a, plan$data, plan$v11, as.double(live))
  if (!identical(res$status[1L], 0L)) or_baseline_native_refuse(res$status, dr_baseline_abort)
  phi <- res$phi
  # This is V11's actual final guard. No eta/probability/odds/Gamma/four-
  # coordinate check is added by this execution optimization.
  cp04_finite(phi, "complete separate DR original-unit IF")
  if (scale == "core") phi <- st$lambda * phi
  dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id"))
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(outcome_pooling = "separate", n_rate_score_evaluations = length(pre$sids),
      n_selection_score_evaluations = 0L, n_selection_block_evaluations = 0L,
      n_regression_score_evaluations = length(pre$gids), n_global_beta_blocks = 0L,
      n_block_beta_blocks = length(pre$gids), n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (unname(cp04_payload(out)["entries"]) > pre$output) dr_route_abort("baseline", "IF output bound underestimated.", "bound")
  out
}
