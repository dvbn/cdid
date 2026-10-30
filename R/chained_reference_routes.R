# CP08 D60: dated baseline IPW and custom-callback routes (selection = "none").
#
# Calls the existing static gates accept keep the accepted P1 pipeline
# unchanged. A request with a genuinely dated covariate, or with changing
# weights anchored by fix_weights, runs here on the full-roster reader with the
# sealed recipe-indexed reference reader (R/chained_reference_inputs.R):
#   * one rate support per (actual pair, cohort, recipe):
#       r = H(W_g I_k S_e) / H(W_g I_k A), read on A rows only;
#   * built-in IPW: ONE risk-set treatment logit per (pair, g, K, recipe) on the
#     observed pair rows of {g} U K with weights W_g / r_G and design z(C^r),
#     normalized treated/control means and the complete CP02 link IF
#       phi = a+ (D - T)/d+ - a- (D - M)/d- + gamma' xi_r - v' xi_beta;
#   * custom callback: the accepted did-style callback is called once per
#     (pair, g, K, recipe) on the observed risk set in canonical A order with
#     covariates z(C^r) and i.weights = W_g / r_G; its returned influence
#     function is composed with the recipe rate channel exactly as on the static
#     custom route, on the full original roster:
#       phi_i = (N/m) 1{i in O} psi_pos(i) + sum_k b_k xi_{r_k, i},
#       b_k = -sum_{j in O, G_j = k} psi_j / (m r_k)
#     (or, O15, b_k = -sum_{j in O, G_j = k} g_j / r_k with the callback's
#     explicit total log-weight gradient g).
# Every unit keeps ONE row with all channels; W and C^r are frozen data. The
# static custom contract (callback output validation, no recentering) is kept.

reference_ipw_abort <- function(message, class = "input", extra = list()) {
  stop(structure(c(list(message = message, call = NULL), extra),
    class = unique(c(paste0("cdid_reference_ipw_", class, "_error"), paste0("cdid_cp04_", class, "_error"),
      "cdid_ipw_fit_error", "cdid_error", "error", "condition"))))
}

validate_reference_ipw_backend <- function(backend) {
  if (!inherits(backend, "cdid_reference_ipw_backend") || !ipw_owner_seal_ok(backend, "reference_ipw_backend"))
    reference_ipw_abort("Expected the immutable reference-input IPW backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) && exists("cp04_reader", backend$reader_ref, inherits = FALSE))
    get("cp04_reader", backend$reader_ref, inherits = FALSE) else NULL
  validate_cp04_reader(reader)
  chained_backend_reference(backend)
  invisible(backend)
}

# Rate support of one (pair, cohort, recipe): weights at the context anchor on
# the cohort's A rows (rate denominator) and its observed pair rows.
reference_ipw_support <- function(reader, spec, ref, ctrl) {
  st <- reader$state
  a_rows <- chained_reference_cohort_rows(st, spec$cohort, "A")
  pair_mass <- repeat_mass <- 0; n_pair <- 0L
  if (length(a_rows)) for (s0 in seq.int(1L, length(a_rows), by = ctrl$max_chunk_rows)) {
    rr <- a_rows[seq.int(s0, min(length(a_rows), s0 + ctrl$max_chunk_rows - 1L))]
    S <- cp04_pair_rows(reader, rr, spec$pair)$S
    W <- chained_reference_weight(ref, spec$reference, rr)
    pair_mass <- pair_mass + sum(W * S); repeat_mass <- repeat_mass + sum(W); n_pair <- n_pair + sum(S)
  }
  if (!is.finite(pair_mass) || !is.finite(repeat_mass) || pair_mass <= 0 || repeat_mass <= 0)
    reference_ipw_abort("Positive repeat and actual-pair cohort support required.", "support")
  rate <- pair_mass / repeat_mass
  cp04_finite(rate, "reference pair observation rate", TRUE)
  list(rate = rate, massA = repeat_mass / st$N, pair_mass = pair_mass, repeat_mass = repeat_mass,
    n_pair = as.integer(n_pair))
}

# Observed risk set of one model in canonical A order: frame rows, local
# cohort position, response, outcome change, recipe design and W / r_G.
reference_ipw_risk_set <- function(reader, spec, supports, ref, ctrl) {
  st <- reader$state
  rows <- integer(); kk <- integer(); dy <- numeric()
  for (i in seq_along(spec$cohorts)) {
    a_rows <- chained_reference_cohort_rows(st, spec$cohorts[i], "A")
    if (length(a_rows)) for (s0 in seq.int(1L, length(a_rows), by = ctrl$max_chunk_rows)) {
      rr <- a_rows[seq.int(s0, min(length(a_rows), s0 + ctrl$max_chunk_rows - 1L))]
      hit <- cp04_pair_rows(reader, rr, spec$pair, TRUE)
      ob <- which(hit$S == 1)
      rows <- c(rows, rr[ob]); kk <- c(kk, rep(i, length(ob))); dy <- c(dy, hit$dy[ob])
    }
  }
  o <- order(st$a_index[rows])
  rows <- rows[o]; kk <- kk[o]; dy <- dy[o]
  if (length(rows) != spec$n_observers)
    reference_ipw_abort("Observed risk support changed from its owner view.", "identity")
  W <- chained_reference_weight(ref, spec$reference, rows)
  rate <- vapply(spec$support_ids, function(s) supports[[s]]$rate, numeric(1))
  X <- chained_reference_x(ref, spec$reference, rows)
  list(rows = rows, k = kk, y = as.double(kk == 1L), dy = dy, X = X, W = W,
    fw = W / ipw_clip_floor(rate)[kk], rate = rate)
}

# Built-in IPW: the CP02 functional and its analytic partials (the same
# formulas as the IPW component of the baseline DR backend).
reference_ipw_logit_fit <- function(reader, spec, supports, ref, glm_control, ctrl) {
  st <- reader$state; N <- st$N
  rs <- reference_ipw_risk_set(reader, spec, supports, ref, ctrl)
  X <- rs$X; y <- rs$y; w <- rs$fw; dy <- rs$dy
  fit <- cp04_logit(X, y, w, glm_control, N, spec)
  names(fit$coefficients) <- colnames(X)
  lv <- ipw_logit_values(X, fit$coefficients)
  wp <- w * y; wm <- w * (1 - y) * lv$odds
  dp <- sum(wp) / N; dm <- sum(wm) / N
  cp04_finite(c(dp, dm), "treated/control normalizers", TRUE)
  T_val <- sum(wp * dy) / sum(wp); M_val <- sum(wm * dy) / sum(wm)
  cp04_finite(c(T_val, M_val), "treated/control means")
  m <- length(spec$cohorts); p <- ncol(X)
  C <- matrix(0, p, m)
  U <- matrix(0, 4L, m, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  score <- w * (y - lv$prob) * X
  residual <- cbind(wp * (dy - T_val) / dp, wm * (dy - M_val) / dm, wp, wm)
  for (k in seq_len(m)) {
    rr <- which(rs$k == k)
    if (!length(rr) || rs$rate[k] <= ipw_clip_eps) next
    scale_k <- -1 / (N * rs$rate[k])
    C[, k] <- colSums(score[rr, , drop = FALSE]) * scale_k
    U[, k] <- colSums(residual[rr, , drop = FALSE]) * scale_k
  }
  V <- matrix(0, 4L, p, dimnames = list(c("T", "M", "dp", "dm"), colnames(X)))
  wm_act <- wm * lv$active
  V[2L, ] <- colSums(wm_act * (dy - M_val) * X) / (N * dm)
  V[4L, ] <- colSums(wm_act * X) / N
  cp04_finite(c(C, U, V), "reference IPW analytic cross derivatives")
  fit$kind <- "logit"
  fit$cross <- C
  fit$partials <- list(upstream = U, beta = V)
  fit$means <- c(T = T_val, M = M_val, dp = dp, dm = dm)
  fit$estimate <- T_val - M_val
  fit$n_risk <- length(rs$rows)
  fit$clipping <- ipw_clip_summary(prob = lv$prob, rates = rs$rate, controls = y == 0, weights = w)
  fit
}

# Custom callback on one model's observed risk set; the rate sensitivity b_k
# is that of the accepted static custom composition.
reference_ipw_custom_fit <- function(reader, spec, supports, ref, custom_est_method, custom_dots, link_id, ctrl) {
  st <- reader$state
  rs <- reference_ipw_risk_set(reader, spec, supports, ref, ctrl)
  n_obs <- length(rs$rows)
  if (n_obs < 2L) cdid_custom_ipw_abort("Custom IPW link has fewer than 2 observed risk-set units.",
    c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error"))
  if (any(!is.finite(rs$fw) | rs$fw <= 0))
    cdid_abort("Nonfinite or nonpositive inverse-rate propensity weights `w_i / r_G`.", c("cdid_ipw_numeric_error", "cdid_ipw_fit_error"))
  # Endpoint outcomes of the observed rows (both endpoints are observed).
  y01 <- reference_ipw_endpoint_outcomes(reader, rs$rows, spec$pair)
  cb_args <- c(list(y1 = y01$y1, y0 = y01$y0, D = rs$y, covariates = rs$X, i.weights = rs$fw, inffunc = TRUE),
    custom_dots)
  raw_res <- do.call(custom_est_method, cb_args)
  val <- custom_ipw_validate_callback_output(raw_res, n_obs = n_obs, link_id = link_id)
  b_rate <- custom_ipw_rate_sensitivity(val, rs$k, rs$rate, length(spec$cohorts), n_obs)
  if (anyNA(b_rate) || any(!is.finite(b_rate)))
    cdid_custom_ipw_abort(paste0("Nonfinite cohort-rate sensitivity `b_rate` on link `", link_id, "`."),
      c("cdid_ipw_numeric_error", "cdid_ipw_fit_error"))
  list(kind = "custom", estimate = as.double(val$ATT), obs_rows = rs$rows, psi_obs = as.double(val$psi),
    b_rate = b_rate, callback_if_mean = mean(val$psi), derivative_source = custom_ipw_derivative_source(val),
    n_risk = n_obs)
}

# Y_s and Y_t of observed pair rows (frame positions) from the sealed panel owner.
reference_ipw_endpoint_outcomes <- function(reader, rows, pair) {
  st <- reader$state
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  post <- panel$postings; pa <- post$post_a_index; ai <- st$a_index[rows]
  ds <- pair$s_calendar_index; dt <- pair$t_calendar_index
  ms <- ipw_owner_match_sorted_double(ai, pa, post$date_start[ds], post$date_end[ds])
  mt <- ipw_owner_match_sorted_double(ai, pa, post$date_start[dt], post$date_end[dt])
  if (anyNA(ms) || anyNA(mt)) reference_ipw_abort("Observed pair endpoint missing from the panel owner.", "identity")
  yy <- ipw_owner_read_y_positions(panel, post$post_obs_index[ms], post$post_obs_index[mt])
  cp04_finite(c(yy$y_left_obs, yy$y_right_obs), "observed pair endpoint outcomes")
  list(y0 = as.double(yy$y_left_obs), y1 = as.double(yy$y_right_obs))
}

fit_reference_ipw_backend <- function(reader, link_views, control, glm_control, reference,
                                      custom_est_method = NULL, custom_dots = list(),
                                      retained_external_entries = 0, external_work_ops = 0) {
  validate_cp04_reader(reader)
  validate_chained_reference_reader(reference)
  ctrl <- validate_ipw_owner_control(control); gc <- validate_ipw_glm_control(glm_control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    reference_ipw_abort("Invalid external pipeline reservation.")
  custom <- is.function(custom_est_method)
  if (!is.list(custom_dots) || is.object(custom_dots)) reference_ipw_abort("`custom_dots` must be a plain list.")
  st <- reader$state; N <- st$N
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (!is.list(link_views) || !length(link_views)) reference_ipw_abort("Supply actual private link views.")
  L <- length(link_views); J <- length(st$cohorts)
  if (length(reference$view_recipe) != L)
    reference_ipw_abort("The reference reader must index exactly the supplied link views.", "integrity")
  p <- length(reference$x_columns)
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries + retained_external_entries +
    reference$accounting$owned_entries
  header <- 4096 + 1024 * L * (J + 1)
  cp04_budget(ctrl, c(base, header), header, c(external_work_ops, 128 * L * L * (J + 1)), header,
    "reference IPW metadata")
  specs <- list(); skeys <- character(); models <- list(); mkeys <- character(); links <- vector("list", L)
  for (ell in seq_len(L)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:"))
      reference_ipw_abort("The reference IPW route requires the supplied comparison-catalog functional key.")
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
      reference_ipw_abort("Legacy nuisance candidates are unsupported.")
    if (!is.null(v$nuisance_ref)) reference_ipw_abort("Cross-pair or unverified nuisance sharing is unsupported.")
    r <- reference$view_recipe[ell]
    target <- match(v$cohorts$g, st$cohorts)
    controls <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    if (anyNA(c(target, controls)) || !isTRUE(all(v$support$by_cohort_pair_count >= 1)))
      reference_ipw_abort("Positive actual-pair cohort support required.", "support")
    cohorts <- c(target, controls)
    pk <- paste(v$pair$s_calendar_index, v$pair$t_calendar_index, sep = ":")
    sid <- vapply(cohorts, function(k) {
      key <- paste(pk, k, r, sep = "|"); j <- match(key, skeys)
      if (is.na(j)) {
        skeys <<- c(skeys, key)
        specs[[length(specs) + 1L]] <<- list(pair = v$pair, cohort = as.integer(k), reference = as.integer(r))
        j <- length(specs)
      }
      as.integer(j)
    }, integer(1))
    # One fit (or one callback call) per scientific (pair, g, K, recipe).
    mkey <- paste(pk, target, paste(controls, collapse = ","), r, sep = "|")
    m <- match(mkey, mkeys)
    if (is.na(m)) {
      mkeys <- c(mkeys, mkey)
      models[[length(models) + 1L]] <- list(pair = v$pair, cohorts = as.integer(cohorts), support_ids = sid,
        reference = as.integer(r), n_observers = as.integer(v$support$n_observed_risk_units), link_id = v$link_id)
      m <- length(models)
    }
    links[[ell]] <- list(link_id = v$link_id, model = as.integer(m), reference = as.integer(r),
      cohort = as.integer(target), source_functional_key = source,
      functional_key = paste0("cmpcat:IPW-baseline-", if (custom) "custom-callback" else "risk-set-logit",
        "-W_g/r|", substring(source, 8L), "|reference:", reference$recipes$key[r]))
  }
  ids <- vapply(links, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids)) reference_ipw_abort("Distinct nonempty link IDs required.")
  S <- length(specs); Tm <- length(models)
  max_n <- max(vapply(models, `[[`, integer(1), "n_observers"))
  max_m <- max(vapply(models, function(x) length(x$cohorts), integer(1)))
  metadata <- unname(cp04_payload(list(specs, models, links))["entries"])
  model_bound <- Tm * (p * p + 6 * p + 4 * max_m + p * max_m + 2048 + if (custom) 2 * max_n else 0) + 64 * S + 8192
  scratch <- 64 * max_n * (p + 16) + 32 * (p + max_m + 4)^2 + 8192
  scan_n <- sum(vapply(models, `[[`, integer(1), "n_observers"))
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- ipw_owner_safe_sum(c(external_work_ops, 128 * (gc$maxit + 2) * scan_n * (p + 2)^2 +
    128 * (S + Tm * max_m) * N * log_a + 128 * model_bound), ctrl$max_work_ops, "reference IPW fit work")
  bound <- metadata + model_bound
  cp04_budget(ctrl, c(base, bound), scratch, work, bound, "reference IPW fits")
  supports <- lapply(specs, function(s) reference_ipw_support(reader, s, reference, ctrl))
  fits <- lapply(seq_len(Tm), function(m) {
    sp <- models[[m]]
    if (custom) reference_ipw_custom_fit(reader, sp, supports, reference, custom_est_method, custom_dots,
      sp$link_id, ctrl) else reference_ipw_logit_fit(reader, sp, supports, reference, gc, ctrl)
  })
  for (ell in seq_len(L)) links[[ell]]$mu <- fits[[links[[ell]]$model]]$estimate
  obj <- structure(list(reader_ref = reader$ref, control = ctrl, glm_control = gc,
    custom = custom, support_specs = specs, supports = supports, treatment_specs = models, treatments = fits,
    links = links, universe = "A",
    law = paste0("IPW:baseline:", if (custom) "custom-callback" else "risk-set-logit",
      ":actual-pair:W_g/r:recipe-indexed-reference-design"),
    accounting = list(owned_entries = 0, borrowed_entries = base,
      retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
      retained_bound = base + bound, fit_scratch_bound = scratch, work_ops_bound = work,
      n_supports = S, n_rate_fits = S, n_treatment_fits = if (custom) 0L else Tm,
      n_custom_callbacks = if (custom) Tm else 0L, n_treatment_columns = p,
      n_selection_fits = 0L, n_regression_fits = 0L, n_reference_recipes = nrow(reference$recipes),
      full_F_masks_per_model = 0L, full_F_by_link_IF_retained = FALSE,
      inverse_matrices_retained = FALSE, peak_rss_measured = FALSE)),
    class = "cdid_reference_ipw_backend")
  obj$reference_ref <- reference$ref
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  if (actual > bound) reference_ipw_abort("Retained reference IPW bound underestimated.", "bound")
  cp04_budget(ctrl, c(base, actual), scratch, work, actual, "retained reference IPW backend")
  obj$accounting$owned_entries <- actual
  ipw_owner_seal(obj, "reference_ipw_backend")
}

reference_ipw_if_budget <- function(backend, B, links) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  tids <- sort(unique(vapply(backend$links[links], `[[`, integer(1), "model")))
  sids <- sort(unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids"))))
  p <- length(chained_backend_reference(backend)$x_columns)
  width <- 4 * length(sids) + length(tids) * (2 * p + 16) + 2 * length(links) + 8
  metadata <- 8192 + 1024 * (length(sids) + length(tids) + length(links)) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  output <- B * (width + 8) + metadata
  scratch <- 64 * B * (width + p + 8) + 128 * (p + 2)^2 + metadata
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  work <- 128 * B * ((length(sids) + length(tids)) * (p + 2)^2 * log_a + length(links)) + 64 * metadata
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    scratch, c(backend$accounting$work_ops_bound, work), output, "reference IPW IF chunk")
  list(tids = tids, sids = sids, output = output, scratch = scratch, work = work)
}

reference_ipw_if_envelope <- function(backend, B, links) reference_ipw_if_budget(backend, B, links)

# Complete original-unit rows for one bounded row chunk: the recipe rate rows
# xi_r = W_g I_k (S - r A) / H(W_g I_k A) on A rows, then each model's channels.
evaluate_reference_ipw_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_reference_ipw_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  ref <- chained_backend_reference(backend)
  valid_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) > 0L &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid_index(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    reference_ipw_abort("Invalid explicit IF row/link batch, scale or details flag.")
  pre <- reference_ipw_if_budget(backend, length(rows), links)
  B <- length(rows); N <- st$N
  ci <- st$cohort_index[rows]; A <- st$a_index[rows] > 0L
  up <- matrix(0, B, length(pre$sids)); hits <- vector("list", length(pre$sids))
  for (j in seq_along(pre$sids)) {
    s <- backend$support_specs[[pre$sids[j]]]; u <- backend$supports[[pre$sids[j]]]
    ia <- which(ci == s$cohort & A)
    if (!length(ia)) { hits[[j]] <- list(obs = integer(), dy = numeric()); next }
    hit <- cp04_pair_rows(reader, rows[ia], s$pair, TRUE)
    W <- chained_reference_weight(ref, s$reference, rows[ia])
    up[ia, j] <- W * (hit$S - u$rate) / u$massA
    ob <- which(hit$S == 1)
    hits[[j]] <- list(obs = ia[ob], dy = hit$dy[ob], W = W[ob])
  }
  models <- vector("list", length(backend$treatments))
  for (tj in pre$tids) {
    spec <- backend$treatment_specs[[tj]]; f <- backend$treatments[[tj]]
    upos <- match(spec$support_ids, pre$sids)
    upstream <- up[, upos, drop = FALSE]
    if (identical(f$kind, "custom")) {
      # (N / m) 1{i in O} psi + sum_k b_k xi_{r_k} on this chunk's original rows.
      phi <- as.vector(upstream %*% f$b_rate)
      pos <- match(rows, f$obs_rows); inO <- which(!is.na(pos))
      if (length(inO)) phi[inO] <- phi[inO] + (N / f$n_risk) * f$psi_obs[pos[inO]]
      models[[tj]] <- list(phi = phi)
      next
    }
    p <- length(f$coefficients)
    beta_score <- matrix(0, B, p); direct <- matrix(0, B, 4L)
    direct[, 3L] <- -f$means["dp"]; direct[, 4L] <- -f$means["dm"]
    for (k in seq_along(spec$cohorts)) {
      h <- hits[[upos[k]]]; rr <- h$obs
      if (!length(rr)) next
      X <- chained_reference_x(ref, spec$reference, rows[rr])
      lv <- ipw_logit_values(X, f$coefficients)
      cp04_finite(c(lv$prob, lv$eta, lv$odds), "reference IPW treatment replay")
      fw <- h$W / ipw_clip_floor(backend$supports[[spec$support_ids[k]]]$rate)
      D <- if (k == 1L) 1 else 0
      wp <- fw * D; wm <- fw * (1 - D) * lv$odds
      beta_score[rr, ] <- fw * (D - lv$prob) * X
      direct[rr, 1L] <- wp * (h$dy - f$means["T"]) / f$means["dp"]
      direct[rr, 2L] <- wm * (h$dy - f$means["M"]) / f$means["dm"]
      direct[rr, 3L] <- direct[rr, 3L] + wp
      direct[rr, 4L] <- direct[rr, 4L] + wm
    }
    xb <- t(cp04_solve(f$H_chol, t(beta_score + upstream %*% t(f$cross))))
    coord <- direct + upstream %*% t(f$partials$upstream) + xb %*% t(f$partials$beta)
    colnames(coord) <- c("T", "M", "dp", "dm")
    cp04_finite(c(xb, coord), "complete reference IPW frame IF")
    models[[tj]] <- if (details) list(phi = coord[, 1L] - coord[, 2L], beta = xb, coordinates = coord,
      direct = direct, beta_score = beta_score, upstream = upstream) else list(phi = coord[, 1L] - coord[, 2L])
  }
  phi <- matrix(0, B, length(links), dimnames = list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id")))
  for (j in seq_along(links)) phi[, j] <- models[[backend$links[[links[j]]]$model]]$phi
  cp04_finite(phi, "complete reference IPW link IF")
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_rate_score_evaluations = length(pre$sids), n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (details) out$details <- list(support_ids = pre$sids, rates = up, models = models, lambda = as.double(A) - st$lambda)
  if (unname(cp04_payload(out)["entries"]) > pre$output) reference_ipw_abort("IF output bound underestimated.", "bound")
  out
}

reference_ipw_public_fit <- function(reader, views, ctrl, retained, work, reference,
                                     custom_est_method = NULL, custom_dots = list()) {
  external <- max(0, retained - reader$accounting$borrowed_owner_entries)
  fit_reference_ipw_backend(reader, views, ctrl$owner_control, ctrl$glm_control, reference,
    custom_est_method = custom_est_method, custom_dots = custom_dots,
    retained_external_entries = external, external_work_ops = work)
}

reference_ipw_public_diagnostics <- function(res, control) {
  b <- res$nuisance_registry$cp04
  bound <- 4 * b$accounting$owned_entries + 4096
  cp04_budget(control$owner_control, c(max(c(0, res$resource_accounting$simultaneous_entries,
    res$resource_accounting$initial_envelope_entries)), bound),
    4 * bound, 8 * bound, bound, "reference IPW public diagnostic copy")
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st_cohorts <- reader$state$cohorts
  link_ids <- vapply(b$links, `[[`, character(1), "link_id")
  fitted_catalog <- data.frame(index = seq_along(b$links), link_id = link_ids,
    role = ifelse(startsWith(link_ids, "cp04:auxiliary:"), "auxiliary", "catalog"),
    cohort = st_cohorts[vapply(b$links, `[[`, integer(1), "cohort")],
    s = vapply(b$links, function(l) as.double(b$treatment_specs[[l$model]]$pair$s), numeric(1)),
    t = vapply(b$links, function(l) as.double(b$treatment_specs[[l$model]]$pair$t), numeric(1)),
    treatment_model = vapply(b$links, `[[`, integer(1), "model"),
    recipe = vapply(b$links, `[[`, integer(1), "reference"),
    mu = vapply(b$links, `[[`, numeric(1), "mu"), stringsAsFactors = FALSE)
  records <- list(
    supports = lapply(seq_along(b$supports), function(j) c(list(id = j, spec = b$support_specs[[j]]), b$supports[[j]])),
    treatment = lapply(seq_along(b$treatments), function(j) {
      f <- b$treatments[[j]]
      c(list(id = j, spec = b$treatment_specs[[j]], kind = f$kind, estimate = f$estimate, n_risk = f$n_risk),
        if (identical(f$kind, "logit")) list(coefficients = f$coefficients, means = f$means, diagnostics = f$diagnostics, clipping = f$clipping) else
          list(callback_if_mean = f$callback_if_mean, b_rate = f$b_rate))
    }))
  diag_copy <- list(records = records, fitted_catalog = fitted_catalog)
  if (unname(cp04_payload(diag_copy)["entries"]) > bound) cp04_abort("Reference IPW diagnostic copy bound underestimated.", "bound")
  list(status = "available", method = if (isTRUE(b$custom)) "custom" else "ipw",
    treatment_model = if (isTRUE(b$custom)) "custom did-style callback per (pair, g, K, recipe)" else
      "risk-set treatment propensity logit per (pair, g, K, recipe)",
    selection = "none", records = records, fitted_catalog = fitted_catalog, links = b$links,
    accounting = b$accounting,
    reference_inputs = chained_reference_record(chained_backend_reference(b)),
    frame_counts = list(n_full_frame_units = res$reporting$inference$n_total,
      n_pair_eligible_A_units = res$frame_owner$n_a_units),
    auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
    interpretation = paste("D60 recipe-indexed baseline", if (isTRUE(b$custom)) "custom-callback" else "IPW",
      "route: rates, logits/callbacks and links keyed by (pair, cohort set, recipe) and fitted once; designs",
      "z(C^r) at each context's extraction date and weights W_g at its anchor; complete F-row influence",
      "functions with the recipe rate channels. No numerical certificate."))
}
