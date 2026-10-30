# Private Corollary 1 control-only shared-slope OR (D47). One borrowed common
# frame and CP04 reader (selection Q, U, known rho and full-F shares a_k^U),
# explicit outcome recipe and actual pairs. Each (pair, cohort) support has ONE
# actual CP04 selection logit and ONE scan with fitted weight
# q = w a_G^U / (rho_G^U max(p_{e,G}^U, eps)); ONE joint blockwise fit has a
# shared one-period slope and one intercept per unique (pair, pool) block. No
# new propensity model, probability denominator or treatment nuisance.

or_corollary1_abort <- function(message, class = "input") {
  stop(structure(list(message = message, call = NULL),
    class = c(paste0("cdid_or_corollary1_", class, "_error"),
      "cdid_ipw_fit_error", "cdid_error", "error", "condition")))
}

# Fitted factor a/(rho max(p, eps)) of observed pair rows; q = w times it in
# both the fitting scan and the IF replay.
or_corollary1_factor <- function(a, rho, prob) a / (rho * ipw_clip_floor(prob))

# Compact logit-channel moments of one observed chunk: for every selection
# column j, M_j += sum q (-(1 - p) 1{p > eps} Q_j) x x' with x = (1, z, dy) on
# all slope columns. Q_j may be negative, so no square-root factor exists.
or_corollary1_moments_add <- function(M, z, dy, q, prob, Q) {
  x <- cbind(1, z, dy); tq <- -q * ((1 - prob) * ipw_clip_floor_active(prob))
  for (j in seq_len(ncol(Q))) M[, , j] <- M[, , j] + crossprod(x, (tq * Q[, j]) * x)
  M
}

# One scan of ONE support after its selection fit. Observed pair rows satisfy
# S <= A <= U, so Z_a, Q and rho are read only there. Keeps the q-weighted
# design summary and the compact moments; both die with fitting.
or_corollary1_support <- function(reader, spec, selection, recipe, ctrl) {
  st <- reader$state; a <- st$a[spec$cohort]; m <- recipe$p + 2L
  acc <- or_pooled_support_start(recipe, spec$transitions)
  moments <- array(0, c(m, m, ncol(st$Z)))
  clip <- ipw_clip_summary(q = numeric(0), weights = numeric(0))
  or_baseline_visit(reader, spec, ctrl, function(rows, S, obs, dy, Z) {
    if (!length(obs)) return(invisible(NULL))
    rr <- rows[obs]
    prob <- cp04_predict_selection(reader, rr, spec, selection)
    q <- st$w[rr] * or_corollary1_factor(a, st$rho[rr], prob)
    cp04_finite(q, "Corollary 1 fitted pair weights", TRUE)
    clip <<- ipw_clip_summary_merge(clip, ipw_clip_summary(q = prob, weights = q))
    acc <<- or_pooled_support_add(acc, recipe, Z, dy, q, "corollary1")
    moments <<- or_corollary1_moments_add(moments, Z[, recipe$slope_index, drop = FALSE],
      dy, q, prob, st$Z[rr, , drop = FALSE])
  })
  if (!is.finite(acc$weight_mass) || acc$weight_mass <= 0 || acc$weight_mass / st$N <= 0)
    or_corollary1_abort("Positive actual-pair cohort support required.", "support")
  cp04_finite(c(acc$sumZ, acc$sumY, moments), "Corollary 1 compact support moments")
  list(design = acc, moments = moments, a = a, pair_mass = acc$weight_mass, n_pair = acc$n,
    clipping = clip)
}

# Per (block, support) forcing on the support's channels: its cohort's a_k
# (t = 1/a, q-weighted factor, F'F = M'M) and its d logit coordinates
# (t = -(1 - p) Q, compact moments). C_alpha = H(q r_b t'), C_beta =
# H(q (d - dbar_b) r_b t') (uncentered d without intercept), C_energy =
# H(q r_b^2 t'). A shared a_k coordinate sums over pairs and blocks. A
# single-cohort block's a_k C_alpha vanishes at the root; its C_beta does not.
or_corollary1_forcing <- function(recipe, fitted, supports, pools, channels, D, N) {
  keep <- fitted$keep; p <- length(keep); m <- recipe$p + 2L
  beta_total <- matrix(0, p, D); energy_total <- numeric(D)
  regressions <- fitted$regressions
  for (b in seq_along(pools)) {
    sids <- pools[[b]]$support_ids; rec <- regressions[[b]]
    fa <- fe <- numeric(D); fb <- matrix(0, p, D)
    cr <- numeric(m); cr[1L] <- -rec$alpha; cr[1L + keep] <- -rec$L * fitted$beta; cr[m] <- 1
    for (s in sids) {
      u <- supports[[s]]; ka <- channels[[s]][1L]; kl <- channels[[s]][-1L]
      g <- or_pooled_factor_moments(recipe, u$design$summary, keep, rec$alpha, fitted$beta)
      Mr <- vapply(seq_along(kl), function(j) as.vector(u$moments[, , j] %*% cr), numeric(m))
      if (recipe$intercept) {
        fa[ka] <- fa[ka] + g$g[1L] / (N * u$a)
        fb[, ka] <- fb[, ka] + (g$g[-1L] - rec$mean_d * g$g[1L]) / (N * u$a)
        fa[kl] <- fa[kl] + Mr[1L, ] / N
        fb[, kl] <- fb[, kl] + (rec$L * Mr[1L + keep, , drop = FALSE] -
          outer(rec$mean_d, Mr[1L, ])) / N
      } else {
        fb[, ka] <- fb[, ka] + g$g / (N * u$a)
        fb[, kl] <- fb[, kl] + rec$L * Mr[1L + keep, , drop = FALSE] / N
      }
      fe[ka] <- fe[ka] + g$e2 / (N * u$a)
      fe[kl] <- fe[kl] + colSums(cr * Mr) / N
    }
    ids <- sort(unique(unlist(channels[sids])))
    cp04_finite(c(fa, fb, fe), "Corollary 1 regression forcing")
    regressions[[b]]$channel_ids <- ids
    regressions[[b]]$forcing_alpha <- fa[ids]
    regressions[[b]]$forcing_beta <- fb[, ids, drop = FALSE]
    regressions[[b]]$forcing_energy <- fe[ids]
    beta_total <- beta_total + fb; energy_total <- energy_total + fe
  }
  list(regressions = regressions, beta = beta_total, energy = energy_total)
}

# ONE treated residual ratio with identical fitted weights T = q in numerator
# and denominator, and its loadings on the target support's channels:
# rate_loading = H(T (r - mu) t')/D, plus H(T t'), H(T r t') and H(T d t'). The
# a_g loading of the ratio is zero at the root and is retained.
or_corollary1_link <- function(recipe, fitted, u, rec, N) {
  keep <- fitted$keep; p <- length(keep); m <- recipe$p + 2L; dq <- dim(u$moments)[3L]
  lm <- or_pooled_link_moments(recipe, keep, u$pair_mass, u$design$sumZ, u$design$sumY,
    rec$L, rec$alpha, fitted$beta)
  d <- u$pair_mass / N; num <- lm$mu * d
  cr <- numeric(m); cr[1L] <- -rec$alpha; cr[1L + keep] <- -rec$L * fitted$beta; cr[m] <- 1
  Mr <- vapply(seq_len(dq), function(j) as.vector(u$moments[, , j] %*% cr), numeric(m))
  ans <- list(mu = lm$mu, c = lm$c, d = d, numerator = num,
    rate_loading = c(lm$residual / (N * u$a * d),
      (Mr[1L, ] - lm$mu * u$moments[1L, 1L, ]) / (N * d)),
    weight_loading = c(d / u$a, u$moments[1L, 1L, ] / N),
    numerator_loading = c(num / u$a, Mr[1L, ] / N),
    design_loading = matrix(c(lm$c * d / u$a, rec$L * u$moments[1L + keep, 1L, ] / N),
      p, 1L + dq))
  cp04_finite(unlist(ans), "Corollary 1 link moments")
  ans
}

# Charges reserved before any fit or scan. Retained: selections, compact
# supports and channel maps, block records and forcing over the global upstream
# coordinates, ONE pooled summary and fit, basis map, link moments and loadings.
# Concurrent: one selection fitting buffer, one row chunk, every support's
# summary and compact moments until forcing, merges and forcing accumulators.
# Work: selection fits and support scans over actual cohort rows, one reserved
# verification rescan, merges, proposal/map, forcing and link moments.
or_corollary1_fit_charges <- function(S, G, L, p, d, D, n_A, B, max_fit_n, visits, maxit,
                                      metadata, header, name_bytes) {
  P2 <- p + 2; d1 <- d + 1; log_a <- 1 + ceiling(log2(n_A + 1))
  models <- 16384 + 1024 * (S + G + L) + 64 * (S + G + 4) * P2^2 + 32 * G * S * P2 +
    64 * L * P2 + 64 * (2 * S + 8 * G + 16) * (P2 + name_bytes) +
    S * (d * d + d + 2048) + 64 * (S + L) * P2 * d1 + 64 * (G + 2) * P2 * (D + 1)
  scratch <- header + metadata + 128 * max_fit_n * (d + 16) + 128 * B * (P2 + d + 8) +
    64 * (2 * S + 2 * G + 8) * P2^2 + 64 * S * P2^2 * d1 + 64 * G * p * (p + 1) +
    64 * P2 * (D + 1)
  scan <- 256 * visits * (P2^2 * d1 * log_a + d1)
  work <- 256 * visits * ((maxit + 2) * d1^2 + d1 + log_a) + 256 * S * d1^3 + 2 * scan +
    256 * G * (S + 1) * P2^3 + 256 * (G + S) * P2^2 * (G + 2) +
    256 * G * S * (P2^2 * d1 + P2 * D) + 256 * L * P2^2 * d1 +
    128 * (metadata + header + L * P2)
  list(models = models, scratch = scratch, work = work)
}

# Fit-time state from finished supports: ONE joint fit, forcing over the global
# upstream coordinates (a_k once per cohort, then d logit coordinates per
# support), link moments and loadings. visit_rows(sid, visit) replays bounded
# observed Z_a chunks of support sid, used only for exact alias reconstruction.
or_corollary1_assemble <- function(recipe, supports, pools, links, channels, D, rank_tol, N,
                                   visit_rows) {
  blocks <- lapply(seq_along(pools), function(b) list(name = paste0("b", b),
    support_ids = pools[[b]]$support_ids, L = pools[[b]]$transitions))
  fitted <- or_pooled_fit_blocks(recipe, lapply(supports, function(u)
      list(summary = u$design$summary, scale = 1, colmax = u$design$colmax)),
    blocks, lapply(links, function(l) list(regression = l$regression,
      support_id = l$support_id)), rank_tol, N, visit_rows, "corollary1")
  forcing <- or_corollary1_forcing(recipe, fitted, supports, pools, channels, D, N)
  for (i in seq_along(links)) {
    l <- links[[i]]
    links[[i]] <- c(l, or_corollary1_link(recipe, fitted, supports[[l$support_id]],
      forcing$regressions[[l$regression]], N))
  }
  beta <- fitted$beta; names(beta) <- fitted$basis$kept_columns
  list(supports = lapply(seq_along(supports), function(s) {
      u <- supports[[s]]
      list(a = u$a, channels = channels[[s]], pair_mass = u$pair_mass, n_pair = u$n_pair,
        sumZ = u$design$sumZ, sumY = u$design$sumY, colmax = u$design$colmax, L = u$design$L,
        clipping = u$clipping)
    }),
    regressions = forcing$regressions, links = links,
    shared = list(p = length(fitted$keep), keep = fitted$keep, beta = fitted$beta,
      coefficients = beta, fit = fitted$fit, summary = fitted$summary, basis = fitted$basis,
      forcing_beta = forcing$beta, forcing_energy = forcing$energy,
      control_support_ids = sort(unique(unlist(lapply(pools, `[[`, "support_ids")))),
      rss = fitted$rss, objective = fitted$objective,
      verification_scan = fitted$verification_scan))
}

validate_or_corollary1_backend <- function(backend) {
  if (!inherits(backend, "cdid_or_corollary1") ||
      !ipw_owner_seal_ok(backend, "or_corollary1_backend"))
    or_corollary1_abort("Expected the immutable private Corollary 1 OR backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) &&
    exists("cp04_reader", backend$reader_ref, inherits = FALSE))
      get("cp04_reader", backend$reader_ref, inherits = FALSE) else NULL
  validate_cp04_reader(reader)
  if (!is.null(backend$reference_ref)) chained_backend_reference(backend)
  invisible(backend)
}

fit_or_corollary1 <- function(reader, link_views, control, glm_control, rank_tol, outcome_recipe,
                              retained_external_entries = 0, external_work_ops = 0,
                              automatic_chunk_rows = FALSE,
                              selection_model = "pair_logit", selection_pooling = "separate") {
  validate_cp04_reader(reader)
  if (!selection_use_legacy_consumer(reader, selection_model, selection_pooling)) {
    return(fit_or_selection_models(
      reader = reader, link_views = link_views, control = control, glm_control = glm_control,
      rank_tol = rank_tol, outcome_recipe = outcome_recipe,
      selection_model = selection_model, selection_pooling = selection_pooling,
      retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
      automatic_chunk_rows = automatic_chunk_rows
    ))
  }
  if (missing(rank_tol) || !gmm_scalar_numeric(rank_tol) || rank_tol <= 0 || rank_tol >= 1)
    or_corollary1_abort("Supply explicit finite private rank_tol in (0, 1).")
  ctrl <- validate_ipw_owner_control(control); gc <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    or_corollary1_abort("Invalid external pipeline reservation.")
  if (!identical(automatic_chunk_rows, isTRUE(automatic_chunk_rows)))
    or_corollary1_abort("Invalid private automatic row-planning flag.")
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (missing(outcome_recipe))
    or_corollary1_abort("Supply the explicit private outcome recipe (intercept flag and Z_a roles).")
  recipe <- or_pooled_check_recipe(outcome_recipe, colnames(frame$Z_a), "corollary1")
  if (!is.list(link_views) || !length(link_views))
    or_corollary1_abort("Supply actual private link views.")
  L <- length(link_views); J <- length(st$cohorts); p <- recipe$p; d <- ncol(st$Z)
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    retained_external_entries
  header <- 4096 + 1024 * L * (J + 1)
  cp04_budget(ctrl, c(base, header), header, c(external_work_ops, 128 * L * L * (J + 1)), header,
    "Corollary 1 OR metadata")
  records <- list()
  for (ell in seq_along(link_views)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:"))
      or_corollary1_abort("Corollary 1 OR requires the supplied comparison-catalog functional key.")
    key_bytes <- nchar(source, type = "bytes")
    cp04_budget(ctrl, c(base, header), header + 16 * key_bytes,
      c(external_work_ops, 128 * L * L * (J + 1) + 64 * key_bytes), header + 16 * key_bytes,
      "Corollary 1 catalog key")
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
      or_corollary1_abort("Legacy nuisance candidates are unsupported.")
    if (!is.null(v$nuisance_ref)) {
      nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair))
        or_corollary1_abort("Cross-pair or unverified nuisance sharing is unsupported.")
    }
    target <- match(v$cohorts$g, st$cohorts)
    control_cohorts <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    if (anyNA(c(target, control_cohorts)) || !isTRUE(all(v$support$by_cohort_pair_count >= 1)))
      or_corollary1_abort("Positive actual-pair cohort support required.", "support")
    or_pooled_transitions(v$pair, "corollary1")
    records[[ell]] <- list(pair = v$pair, target = target, controls = control_cohorts,
      link_id = v$link_id, source_functional_key = source,
      functional_key = paste0("cmpcat:OR-Cor1-", st$universe, "-shared-slope-WLS-a/(rho*p)|",
        substring(source, 8L)))
    cp04_budget(ctrl, c(base, header), header + cp04_payload(records)["entries"],
      c(external_work_ops, 128 * L * L * (J + 1) + 64 * cp04_payload(records)["char_bytes"]),
      header + cp04_payload(records)["entries"], "Corollary 1 catalog metadata")
  }
  ids <- vapply(records, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids))
    or_corollary1_abort("Distinct nonempty link IDs required.")
  # The canonical shared catalog: one block per unique (actual pair, sorted pool).
  catalog <- or_baseline_catalog(records, st$cohorts, recipe)
  specs <- catalog$specs; pools <- catalog$pools; links <- catalog$links
  for (i in seq_along(pools)) {
    pr <- specs[[pools[[i]]$support_ids[1L]]]$pair
    pools[[i]]$block_key <- paste0("pair:", pr$s, "-", pr$t, "|pool:",
      paste(st$cohorts[pools[[i]]$cohorts], collapse = ","), "|", recipe$key,
      "|weights:w*a_U/(rho_U*p_U)|U:", st$universe)
  }
  # Selection key: actual pair and cohort inside this reader's U/Q/rho/weights.
  selection_specs <- lapply(specs, function(s) list(pair = s$pair, cohort = s$cohort))
  S <- length(specs); G <- length(pools)
  cohorts_a <- sort(unique(vapply(specs, `[[`, integer(1), "cohort")))
  n_a <- length(cohorts_a); D <- n_a + S * d
  channels <- lapply(seq_len(S), function(s)
    c(match(specs[[s]]$cohort, cohorts_a), n_a + (s - 1L) * d + seq_len(d)))
  upstream <- list(n_a = n_a, d = d, D = D,
    cohort = c(cohorts_a, rep(vapply(specs, `[[`, integer(1), "cohort"), each = d)),
    support = c(integer(n_a), rep(seq_len(S), each = d)),
    column = c(integer(n_a), rep(seq_len(d), S)))
  metadata <- unname(cp04_payload(list(specs, pools, links, selection_specs, upstream))["entries"])
  cohort_n <- vapply(specs, function(s)
    st$cohort_end[s$cohort] - st$cohort_start[s$cohort] + 1, numeric(1))
  charges <- or_corollary1_fit_charges(S, G, L, p, d, D, st$n_A, min(st$N, ctrl$max_chunk_rows),
    max(cohort_n), sum(cohort_n), gc$maxit, metadata, header,
    sum(nchar(c(recipe$columns, colnames(st$Z)), type = "bytes")))
  scratch <- charges$scratch
  work <- ipw_owner_safe_sum(c(charges$work, external_work_ops), ctrl$max_work_ops,
    "Corollary 1 OR pipeline fit work")
  bound <- metadata + charges$models
  cp04_budget(ctrl, c(base, bound), scratch, work, bound, "Corollary 1 OR fits")
  selections <- lapply(selection_specs, function(s) cp04_selection_fit(reader, s, gc, ctrl))
  supports <- lapply(seq_len(S), function(j)
    or_corollary1_support(reader, specs[[j]], selections[[j]], recipe, ctrl))
  state <- or_corollary1_assemble(recipe, supports, pools, links, channels, D, rank_tol, st$N,
    function(sid, visit) or_baseline_visit(reader, specs[[sid]], ctrl,
      function(rows, S, obs, dy, Z) visit(Z)))
  rm(supports)
  obj <- structure(list(reader_ref = reader$ref, control = ctrl, glm_control = gc,
    rank_tol = rank_tol, outcome_recipe = recipe, support_specs = specs,
    selection_specs = selection_specs, selections = selections, upstream = upstream,
    supports = state$supports, regression_specs = pools, regressions = state$regressions,
    shared = state$shared, links = state$links,
    law = paste0("OR:Cor1:", st$universe, ":control-only:shared-one-period-slope:",
      if (recipe$intercept) "block-intercepts" else "no-intercept",
      ":actual-pair:w*a_U/(rho_U*p_U):common-frame-Z_a"),
    accounting = list(owned_entries = 0, borrowed_entries = base,
      retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
      retained_bound = base + bound, fit_scratch_bound = scratch, work_ops_bound = work,
      support_cohort_row_visits = sum(cohort_n), automatic_chunk_rows = automatic_chunk_rows,
      n_cohort_summaries = S, n_rate_fits = n_a, rates_from_reader_shares = TRUE,
      n_selection_fits = S, n_upstream_coordinates = D, n_regression_fits = 1L,
      n_regression_blocks = G, n_local_intercepts = if (recipe$intercept) G else 0L,
      n_shared_slopes = state$shared$p, n_slope_columns = p,
      n_aliased_slopes = length(state$shared$basis$aliased_columns),
      n_support_scans = S * (1L + state$shared$verification_scan), n_treatment_fits = 0L,
      full_F_masks_per_model = 0L, full_F_by_link_IF_retained = FALSE,
      inverse_matrices_retained = FALSE, dense_dummy_regression = FALSE,
      peak_rss_measured = FALSE)),
    class = "cdid_or_corollary1")
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  if (actual > bound) or_corollary1_abort("Retained Corollary 1 OR bound underestimated.", "bound")
  cp04_budget(ctrl, c(base, actual), scratch, work, actual, "retained Corollary 1 OR")
  obj$accounting$owned_entries <- actual
  ipw_owner_seal(obj, "or_corollary1_backend")
}

# Requested link blocks plus, whenever a shared slope exists (or details), EVERY
# block and control support with all of their upstream coordinates: xi_beta is
# global to each link batch. Affine in B, as the automatic row planner assumes.
# Detail matrices use backend support/regression/coordinate order; slopes follow
# the kept basis; shares follow reader$state$cohorts; details stay on the frame scale.
or_corollary1_if_envelope <- function(backend, B, links, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p <- backend$shared$p; G <- length(backend$regression_specs); J <- length(st$cohorts)
  d <- backend$upstream$d
  sets <- or_pooled_if_sets(backend, links, details)
  gids <- sets$gids; beta_gids <- sets$beta_gids; sids <- sets$sids
  uids <- sort(unique(unlist(lapply(backend$supports[sids], `[[`, "channels"))))
  nl <- length(links); ns <- length(sids); nu <- length(uids)
  width <- 4 * (nu + ns) + 4 * (p + 1) * (length(beta_gids) + 1) + 4 * length(gids) +
    4 * nl * (p + 3) + if (details) 2 * J + 8 else 8
  metadata <- 8192 + 1024 * (ns + nu + G + nl + J) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  output <- B * (width + 8) + metadata
  scratch <- 64 * B * (width + p + 2 * d + 8) + 128 * ((p + 2)^2 + (d + 1)^2) + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- 128 * B * ((ns + 1) * ((p + 2)^2 + (d + 1)^2) * log_a +
    (length(beta_gids) + length(gids)) * (p + 1) * (p + nu) + J + nl * (p + nu + 2)) +
    64 * metadata
  if (isTRUE(backend$uses_selection_models)) {
    # The global outcome slope already brings in all control supports. A
    # pooled selection slope also brings in every cell of its selection block.
    # Compose that accepted affine envelope before admission and row planning.
    selection_pre <- selection_if_budget(backend$selection_backend, B,
      unique(backend$q_of_s[sids]), details = details)
    scratch <- scratch + selection_pre$output + selection_pre$scratch
    work <- work + selection_pre$work
  }
  list(gids = gids, beta_gids = beta_gids, sids = sids, uids = uids, output = output,
    scratch = scratch, work = work)
}

or_corollary1_if_budget <- function(backend, B, links, details = FALSE) {
  pre <- or_corollary1_if_envelope(backend, B, links, details)
  cp04_budget(backend$control, c(backend$accounting$owned_entries,
    backend$accounting$borrowed_entries), pre$scratch,
    c(backend$accounting$work_ops_bound, pre$work), pre$output, "Corollary 1 OR IF chunk")
  pre
}

evaluate_or_corollary1_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_or_corollary1_backend(backend)
  if (isTRUE(backend$uses_selection_models))
    return(evaluate_or_selection_if_chunk(backend, rows, links, scale = scale, details = details))
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  valid_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) &&
    all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid_index(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    or_corollary1_abort("Invalid explicit IF row/link batch, scale or details flag.")
  pre <- or_corollary1_if_budget(backend, length(rows), links, details)
  B <- length(rows); frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  kept <- backend$outcome_recipe$slope_index[backend$shared$keep]
  ci <- st$cohort_index[rows]; w <- st$w[rows]; U <- st$U[rows]
  A <- as.double(st$a_index[rows] > 0L)
  up <- matrix(0, B, length(pre$uids))
  # Shared a_k^U rows, ONE column per cohort, on U only (structural zero for U = A).
  for (j in which(pre$uids <= backend$upstream$n_a)) {
    k <- backend$upstream$cohort[pre$uids[j]]; ix <- which(ci == k & U)
    up[ix, j] <- w[ix] * (A[ix] - st$a[k]) / st$mass_U[k]
  }
  # Logit rows and observed pairs of each support; Q and rho are read on U only.
  # Sparse hits are temporary and at most one row chunk, never retained by fit.
  hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]; sp <- backend$support_specs[[sid]]; f <- backend$selections[[sid]]
    ix <- which(ci == sp$cohort & U)
    hit <- cp04_pair_rows(reader, rows[ix], sp$pair, TRUE)
    pr <- cp04_predict_selection(reader, rows[ix], backend$selection_specs[[sid]], f)
    if (length(ix)) up[ix, match(backend$supports[[sid]]$channels[-1L], pre$uids)] <-
      t(cp04_solve(f$H_chol, t(w[ix] * (hit$S - pr) * st$Z[rows[ix], , drop = FALSE])))
    obs <- which(hit$S == 1); rr <- rows[ix[obs]]
    hits[[j]] <- list(obs = ix[obs], dy = hit$dy[obs],
      z = frame$Z_a[st$a_index[rr], kept, drop = FALSE],
      v = or_corollary1_factor(st$a[sp$cohort], st$rho[rr], pr[obs]))
  }
  core <- or_pooled_if_core(backend, st$N, w, hits, up, pre, links, details)
  phi <- core$phi
  dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id"))
  cp04_finite(phi, "complete Corollary 1 OR link IF")
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_rate_score_evaluations = sum(pre$uids <= backend$upstream$n_a),
      n_selection_score_evaluations = length(pre$sids),
      n_upstream_coordinates = length(pre$uids),
      n_regression_score_evaluations = length(union(pre$gids, pre$beta_gids)),
      n_global_beta_blocks = length(pre$beta_gids), output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (details) {
    tau <- eta <- matrix(0, B, length(st$cohorts))
    for (k in seq_along(st$cohorts)) {
      tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
      eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
    }
    out$details <- list(support_ids = pre$sids, upstream_ids = pre$uids, upstream = up,
      beta_regression_ids = pre$beta_gids, beta = core$beta, beta_score = core$beta_score,
      beta_block_scores = core$beta_block_scores, regression_ids = pre$gids,
      alpha = core$alpha, alpha_score = core$alpha_score, direct = core$direct,
      denominator = core$denominator, numerator = core$numerator, design = core$design,
      objective = core$objective, tau = tau, eta = eta, lambda = A - st$lambda)
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    or_corollary1_abort("IF output bound underestimated.", "bound")
  out
}

replay_or_corollary1_if <- function(backend, chunk_rows, link_batch_size, callback,
                                    scale = "frame", details = FALSE) {
  validate_or_corollary1_backend(backend)
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  LB <- ipw_owner_int_bound(link_batch_size, "link_batch_size")
  if (B > backend$control$max_chunk_rows || !is.function(callback) ||
      !identical(details, isTRUE(details)) || !is.character(scale) || length(scale) != 1L ||
      is.na(scale) || !scale %in% c("frame", "core"))
    or_corollary1_abort("Invalid explicit replay controls/callback.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N; L <- length(backend$links)
  total <- backend$accounting$work_ops_bound
  # Reserve ALL fit/external/replay work and each batch output before callbacks.
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1))
    pre <- or_corollary1_if_budget(backend, min(B, N), ll, details)
    total <- ipw_owner_safe_sum(c(total, ceiling(N / B) * pre$work),
      backend$control$max_work_ops, "total Corollary 1 OR pipeline replay work")
  }
  batches <- 0L
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1))
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1))
      batch <- evaluate_or_corollary1_if_chunk(backend, rr, ll, scale, details)
      callback(batch)
      rm(batch)
      batches <- batches + 1L
    }
  }
  invisible(list(n_frame_rows = N, n_links = L, batches = batches, scale = scale,
    normalization_count = if (scale == "core") reader$state$n_A else N,
    work_ops_bound = total, all_F_identities = TRUE))
}
