# CP08 O15 (D64): custom did-style IPW callbacks under conditional_mtar with
# estimated observation probabilities (selection_est_method = "ipw" only).
#
# The accepted CP04 selection models (fit_selection_models, every model,
# pooling and universe of the IPW consumer) are fitted once on the union of the
# post and both auxiliary catalogs. Each scientific (pair, g, K[, recipe]) calls
# the callback once on its observed risk set in canonical A order with
#   i.weights a_j = w_j a_k / (rho_j max(p_{e,k}(X_j), eps))
# (the CP04 synthetic weights, recipe-indexed on the D60 path; d_u log a_j has
# no selection path where the floor binds). The callback
# must return the optional O15 field `log_weight_gradient`
#   g_j = d ATT / d log a_j  (total derivative, internal refits included),
# from which b = sum_j g_j d_u log a_j on the selection upstream coordinates u
# (repeat rates a_k and logit coefficients; pair_channels of each cohort), and
#   phi_i = (N/m) 1{i in O} psi_pos(i) + b' xi_{u,i}
# on every original row of F (evaluate_selection_upstream_chunk; no A-frame
# transport). The rate channel enters once (inside b); with U = A the rate IF
# is structurally zero. Without g the request is refused: the identity g = psi/m
# needs prerequisites (H1/H2) that are not verified, and g is never inferred.

validate_cp04_custom_backend <- function(backend) {
  if (!inherits(backend, "cdid_cp04_custom_backend") || !ipw_owner_seal_ok(backend, "cp04_custom_backend"))
    cp04_abort("Expected the immutable conditional_mtar custom-callback backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) && exists("cp04_reader", backend$reader_ref, inherits = FALSE))
    get("cp04_reader", backend$reader_ref, inherits = FALSE) else NULL
  validate_cp04_reader(reader)
  if (!is.null(backend$reference_ref)) chained_backend_reference(backend)
  invisible(backend)
}

# One callback call on one model's observed risk set (canonical A order).
cp04_custom_callback_fit <- function(reader, spec, sel_fit, up_map, frame, ref, custom_est_method, custom_dots, ctrl) {
  st <- reader$state
  rows <- integer(); kk <- integer()
  for (i in seq_along(spec$cohorts)) {
    a_rows <- chained_reference_cohort_rows(st, spec$cohorts[i], "A")
    if (length(a_rows)) for (s0 in seq.int(1L, length(a_rows), by = ctrl$max_chunk_rows)) {
      rr <- a_rows[seq.int(s0, min(length(a_rows), s0 + ctrl$max_chunk_rows - 1L))]
      S <- cp04_pair_rows(reader, rr, spec$pair)$S
      rows <- c(rows, rr[S == 1]); kk <- c(kk, rep(i, sum(S == 1)))
    }
  }
  o <- order(st$a_index[rows]); rows <- rows[o]; kk <- kk[o]
  n_obs <- length(rows)
  if (n_obs != spec$n_observers) cp04_abort("Observed risk support changed from its owner view.", "identity")
  if (n_obs < 2L) cdid_custom_ipw_abort("Custom IPW link has fewer than 2 observed risk-set units.",
    c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error"))
  u_ids <- sort(unique(unlist(up_map$pair_channels[spec$selection_ids])))
  fw <- numeric(n_obs); LD <- matrix(0, n_obs, length(u_ids))
  for (i in seq_along(spec$cohorts)) {
    sub <- which(kk == i)
    if (!length(sub)) next
    q <- spec$selection_ids[i]
    ld <- selection_pair_cohort_log_deriv(sel_fit, up_map, q, rows[sub])
    fw[sub] <- ld$pred$weight
    LD[sub, match(up_map$pair_channels[[q]], u_ids)] <- ld$log_deriv
  }
  if (any(!is.finite(fw) | fw <= 0))
    cdid_abort("Nonfinite or nonpositive estimated-observation weights `i.weights`.", c("cdid_ipw_numeric_error", "cdid_ipw_fit_error"))
  y01 <- reference_ipw_endpoint_outcomes(reader, rows, spec$pair)
  X <- if (is.null(ref)) frame$Z_a[st$a_index[rows], , drop = FALSE] else chained_reference_x(ref, spec$reference, rows)
  cb_args <- c(list(y1 = y01$y1, y0 = y01$y0, D = as.double(kk == 1L), covariates = X, i.weights = fw, inffunc = TRUE),
    custom_dots)
  val <- custom_ipw_validate_callback_output(do.call(custom_est_method, cb_args), n_obs = n_obs, link_id = spec$link_id)
  if (is.null(val$g))
    cdid_custom_ipw_abort(paste0("Custom IPW callbacks under selection = \"conditional_mtar\" (estimated observation ",
      "probabilities) require the explicit total log-weight gradient `log_weight_gradient` on link `",
      spec$link_id, "`: the identity g = att.inf.func / m needs prerequisites (H1/H2) that are not verified here, and ",
      "g is never inferred from att.inf.func."),
      c("cdid_custom_ipw_derivative_error", "cdid_custom_ipw_output_error", "cdid_ipw_fit_error"))
  b <- as.vector(crossprod(LD, val$g))
  cp04_finite(b, "custom callback upstream sensitivity b")
  list(kind = "custom", estimate = as.double(val$ATT), obs_rows = rows, psi_obs = as.double(val$psi),
    g = as.double(val$g), b = b, upstream_ids = as.integer(u_ids), n_risk = n_obs,
    callback_if_mean = mean(val$psi), derivative_source = custom_ipw_derivative_source(val))
}

fit_cp04_custom_backend <- function(reader, link_views, control, glm_control,
                                    retained_external_entries = 0, external_work_ops = 0,
                                    selection_model = "pair_logit", selection_pooling = "separate",
                                    reference = NULL, custom_est_method, custom_dots = list()) {
  validate_cp04_reader(reader)
  ctrl <- validate_ipw_owner_control(control); gc <- validate_ipw_glm_control(glm_control)
  if (!is.function(custom_est_method)) cp04_abort("Supply the custom callback.", "integrity")
  if (!is.list(custom_dots) || is.object(custom_dots)) cp04_abort("`custom_dots` must be a plain list.")
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) || !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    cp04_abort("Invalid external pipeline reservation.")
  st <- reader$state; N <- st$N
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (!is.list(link_views) || !length(link_views)) cp04_abort("Supply actual private link views.")
  sel_fit <- fit_selection_models(reader = reader, link_views = link_views, selection_model = selection_model,
    selection_pooling = selection_pooling, control = ctrl, glm_control = gc,
    retained_external_entries = retained_external_entries, external_work_ops = external_work_ops, reference = reference)
  up_map <- build_selection_upstream_map(sel_fit)
  p <- if (is.null(reference)) frame$n_covariates else length(reference$x_columns)
  L <- length(link_views); models <- list(); mkeys <- character(); links <- vector("list", L)
  for (i in seq_len(L)) {
    v <- link_views[[i]]
    validate_ipw_link_view(v, frame_owner = frame)
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:"))
      cp04_abort("CP04 requires the supplied comparison-catalog functional key.", "unsupported")
    pair <- list(s = v$pair$s, t = v$pair$t, s_calendar_index = v$pair$s_calendar_index,
      t_calendar_index = v$pair$t_calendar_index)
    target <- match(v$cohorts$g, st$cohorts); controls <- match(sort(v$cohorts$K, method = "radix"), st$cohorts)
    cohorts <- c(target, controls)
    if (anyNA(cohorts) || !isTRUE(all(v$support$by_cohort_pair_count >= 1)))
      cp04_abort("Positive actual-pair cohort support required.", "support")
    r <- if (is.null(reference)) NULL else reference$view_recipe[i]
    si <- vapply(cohorts, function(k) resolve_pair_cohort_id(sel_fit, if (is.null(r)) list(pair = pair, cohort = k) else
      list(pair = pair, cohort = k, reference = r)), integer(1L))
    mkey <- paste(pair$s_calendar_index, pair$t_calendar_index, target, paste(controls, collapse = ","),
      if (is.null(r)) "" else r, sep = "|")
    m <- match(mkey, mkeys)
    if (is.na(m)) {
      mkeys <- c(mkeys, mkey)
      models[[length(models) + 1L]] <- list(pair = pair, cohorts = as.integer(cohorts), selection_ids = as.integer(si),
        n_observers = as.integer(v$support$n_observed_risk_units), link_id = v$link_id,
        reference = if (is.null(r)) NULL else as.integer(r))
      m <- length(models)
    }
    links[[i]] <- list(link_id = v$link_id, model = as.integer(m), source_functional_key = source,
      functional_key = paste0("cmpcat:CP04-", st$universe, "-", selection_model, "-", selection_pooling,
        "-custom-callback-explicit-log-weight-gradient|", substring(source, 8L),
        if (is.null(r)) "" else paste0("|reference:", reference$recipes$key[r])))
  }
  ids <- vapply(links, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1L))
  if (any(!nzchar(ids)) || anyDuplicated(ids)) cp04_abort("Distinct nonempty link IDs required.")
  Tm <- length(models); max_n <- max(vapply(models, `[[`, integer(1L), "n_observers"))
  model_bound <- Tm * (2 * max_n + 2 * up_map$D + 2048) + 8192
  scratch <- 64 * max_n * (p + up_map$D + 16) + 8192
  work <- ipw_owner_safe_sum(c(sel_fit$accounting$work_ops_bound, external_work_ops,
    128 * sum(vapply(models, function(x) x$n_observers * (p + up_map$D + 4) + N * length(x$cohorts), numeric(1L)))),
    ctrl$max_work_ops, "custom conditional_mtar callback fits")
  base <- sel_fit$accounting$owned_entries + sel_fit$accounting$borrowed_entries
  cp04_budget(ctrl, c(base, model_bound), scratch, work, model_bound, "custom conditional_mtar callback links")
  fits <- lapply(models, function(sp) cp04_custom_callback_fit(reader, sp, sel_fit, up_map, frame, reference,
    custom_est_method, custom_dots, ctrl))
  for (i in seq_len(L)) links[[i]]$mu <- fits[[links[[i]]$model]]$estimate
  obj <- structure(list(uses_selection_models = TRUE, selection_model = selection_model,
    selection_pooling = selection_pooling, selection_backend = sel_fit, up_map = up_map,
    selection_specs = sel_fit$block_specs, selections = selection_blocks_diagnostic_list(sel_fit),
    treatment_specs = models, treatments = fits, links = links, reader_ref = reader$ref, control = ctrl,
    glm_control = gc, custom = TRUE, universe = st$universe, justification = st$justification,
    law = paste0("CP04:", st$universe, ":", selection_model, ":", selection_pooling,
      ":custom-callback:explicit-log-weight-gradient",
      if (is.null(reference)) ":w*a_U/(rho_U*p_U):common-frame-Z_a" else ":W_g*a_U/(rho_U*p_U):recipe-indexed-reference-design"),
    accounting = list(owned_entries = 0, borrowed_entries = reader$accounting$owned_entries +
        reader$accounting$borrowed_owner_entries + retained_external_entries,
      retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
      retained_bound = base + model_bound, fit_scratch_bound = max(sel_fit$accounting$fit_scratch_bound, scratch),
      work_ops_bound = work, n_selection_fits = length(sel_fit$blocks), n_selection_cells = length(sel_fit$cells),
      n_selection_blocks = length(sel_fit$blocks), n_selection_parameters = length(sel_fit$coefficients),
      n_treatment_fits = 0L, n_custom_callbacks = Tm, n_rate_fits = up_map$n_a,
      full_F_masks_per_model = 0L, full_F_by_link_IF_retained = FALSE, inverse_matrices_retained = FALSE,
      peak_rss_measured = FALSE)),
    class = "cdid_cp04_custom_backend")
  if (!is.null(reference)) {
    obj$reference_ref <- reference$ref
    obj$accounting$n_reference_recipes <- nrow(reference$recipes)
  }
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  obj$accounting$owned_entries <- actual
  cp04_budget(ctrl, c(obj$accounting$borrowed_entries, actual), obj$accounting$fit_scratch_bound, work, actual,
    "retained custom conditional_mtar backend")
  ipw_owner_seal(obj, "cp04_custom_backend")
}

cp04_custom_if_budget <- function(backend, B, links, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  tids <- sort(unique(vapply(backend$links[links], `[[`, integer(1L), "model")))
  uids <- sort(unique(unlist(lapply(backend$treatments[tids], `[[`, "upstream_ids"))))
  q_ids <- unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "selection_ids")))
  width <- length(uids) + 2 * length(tids) + length(links) + 8 +
    if (details) length(uids) * (1 + length(tids)) + 3 * length(tids) else 0
  metadata <- 8192 + 1024 * (length(tids) + length(links)) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1L), "link_id"), type = "bytes"))
  output <- B * (width + 8) + metadata
  scratch <- 64 * B * (width + 8) + metadata
  work <- 128 * B * (length(uids) + length(tids) + length(links)) * (1 + ceiling(log2(st$n_A + 1))) + 64 * metadata
  sp <- selection_if_budget(backend$selection_backend, B, q_ids, details = FALSE)
  scratch <- scratch + sp$output + sp$scratch
  work <- work + sp$work
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    scratch, c(backend$accounting$work_ops_bound, work), output, "custom conditional_mtar IF chunk")
  list(tids = tids, uids = uids, output = output, scratch = scratch, work = work)
}

# Complete original-unit rows: (N/m) 1{i in O} psi + b' xi_u on F.
evaluate_cp04_custom_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_cp04_custom_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  valid_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) > 0L &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid_index(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    cp04_abort("Invalid explicit IF row/link batch, scale or details flag.")
  pre <- cp04_custom_if_budget(backend, length(rows), links, details)
  B <- length(rows); N <- st$N
  up <- if (length(pre$uids)) evaluate_selection_upstream_chunk(backend$selection_backend, backend$up_map, rows,
    pre$uids)$up else matrix(0, B, 0L)
  models <- vector("list", length(backend$treatments))
  for (tj in pre$tids) {
    f <- backend$treatments[[tj]]
    xu <- up[, match(f$upstream_ids, pre$uids), drop = FALSE]
    rate_sel <- as.vector(xu %*% f$b)
    direct <- numeric(B)
    pos <- match(rows, f$obs_rows); inO <- which(!is.na(pos))
    if (length(inO)) direct[inO] <- (N / f$n_risk) * f$psi_obs[pos[inO]]
    models[[tj]] <- if (details) list(phi = direct + rate_sel, direct = direct, upstream_channel = rate_sel,
      upstream = xu) else list(phi = direct + rate_sel)
  }
  phi <- matrix(0, B, length(links), dimnames = list(NULL, vapply(backend$links[links], `[[`, character(1L), "link_id")))
  for (j in seq_along(links)) phi[, j] <- models[[backend$links[[links[j]]]$model]]$phi
  cp04_finite(phi, "complete custom conditional_mtar link IF")
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_callback_models = length(pre$tids), n_upstream_coordinates = length(pre$uids),
      output_entries_bound = pre$output, concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (details) out$details <- list(upstream_ids = pre$uids, upstream = up, models = models)
  if (unname(cp04_payload(out)["entries"]) > pre$output) cp04_abort("IF output bound underestimated.", "bound")
  out
}
