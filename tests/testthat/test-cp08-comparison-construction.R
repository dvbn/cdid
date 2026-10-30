# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 comparison construction:
# 1. Shared lightweight catalog row descriptor extraction (`extract_pipeline_catalog_row_specs`,
#    `pipeline_metadata_descriptor_preflight`, `pipeline_metadata_descriptor_source`,
#    `pipeline_metadata_row_spec_fields`) across CP04 and IPW registries, including
#    S02 pre-allocation admission, large-E pool entry bounds, and caller controls.
# 2. Indexed baseline OR/DR comparison catalog (`or_baseline_catalog`,
#    `or_baseline_catalog_key`), including S01 undated ignored-reference pool
#    identity, dated reference separation, attribute/control-order collisions,
#    and support collision fallback.
# 3. Lazy outcome-independent DR treatment preparation (`dr_baseline_shared_prepared`,
#    `dr_baseline_treatment_fit`), borrowing invariant sorted observer rows, designs,
#    weights, cohort groupings, and score/design sums across the 189 logit-reusing
#    models while preserving each pair's dy, L * Z_or, residuals, means, and partials.
#
# STATUS: UNEXECUTED at authoring time (no R, Rscript, parser, build, install,
# or fit was run in this continuation). Root owns installed numerical execution.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Frozen verbatim reference functions spliced from input/package/R
# ---------------------------------------------------------------------------
cc_frozen <- new.env(parent = api)
cc_frozen$pins <- list(
  chained_or_baseline = list(sha256 = "083b759fa8e6a0669ef6a488f7f362ea1f37c77e59c7a3ada0337c980ac0b014", lines = c(743L, 806L)),
  chained_dr_backend = list(sha256 = "5a32a76014d03e7a8d0b610c2e9fc4b37723b717a46554ca79fac56af30e5efc", lines = c(233L, 401L)),
  chained_pipeline_metadata_fields = list(sha256 = "4151a70992c8977b67ada7e93b12486a58d8e6eb1ed243281cde44f599c1c435", lines = c(1648L, 1691L)),
  chained_pipeline_metadata_extract = list(sha256 = "4151a70992c8977b67ada7e93b12486a58d8e6eb1ed243281cde44f599c1c435", lines = c(1631L, 1643L))
)

cc_frozen$or_baseline_catalog <- function(records, cohorts, recipe, reference = NULL) {
  dated <- !is.null(reference)
  pairs <- pools <- links <- list()
  for (ell in seq_along(records)) {
    x <- records[[ell]]
    same <- which(vapply(pairs, identical, logical(1), x$pair))
    if (length(same)) pi <- same[1L] else {
      pairs[[length(pairs) + 1L]] <- x$pair; pi <- length(pairs)
    }
    spec <- if (dated) list(pair_id = pi, cohorts = x$controls, reference = x$reference) else
      list(pair_id = pi, cohorts = x$controls)
    same <- which(vapply(pools, identical, logical(1), spec))
    if (length(same)) ri <- same[1L] else {
      pools[[length(pools) + 1L]] <- spec; ri <- length(pools)
    }
    links[[ell]] <- list(regression = as.integer(ri), pair_id = as.integer(pi),
      cohort = x$target, link_id = x$link_id, source_functional_key = x$source_functional_key,
      functional_key = x$functional_key)
    if (dated) links[[ell]]$reference <- x$reference
  }
  transitions <- vapply(pairs, or_pooled_transitions, numeric(1), route = "baseline")
  po <- order(vapply(pairs, function(x) as.double(x$s_calendar_index), numeric(1)),
    vapply(pairs, function(x) as.double(x$t_calendar_index), numeric(1)))
  pairs <- pairs[po]; transitions <- transitions[po]; pair_map <- match(seq_along(po), po)
  for (i in seq_along(pools)) pools[[i]]$pair_id <- as.integer(pair_map[pools[[i]]$pair_id])
  qo <- order(vapply(pools, function(s) paste(sprintf("%012d", c(s$pair_id, s$cohorts, s$reference)),
    collapse = "|"), character(1)), method = "radix")
  pools <- pools[qo]; pool_map <- match(seq_along(qo), qo)
  for (i in seq_along(links)) {
    links[[i]]$pair_id <- as.integer(pair_map[links[[i]]$pair_id])
    links[[i]]$regression <- as.integer(pool_map[links[[i]]$regression])
  }
  specs <- list()
  for (pi in seq_along(pairs)) {
    refs <- if (dated) sort(unique(c(unlist(lapply(pools, function(s) if (s$pair_id == pi) s$reference else integer())),
      unlist(lapply(links, function(s) if (s$pair_id == pi) s$reference else integer()))))) else NA_integer_
    for (r in refs) {
      in_r <- function(s) s$pair_id == pi && (!dated || s$reference == r)
      used <- sort(unique(c(unlist(lapply(pools, function(s)
        if (in_r(s)) s$cohorts else integer())),
        unlist(lapply(links, function(s) if (in_r(s)) s$cohort else integer())))))
      for (k in used) {
        sp <- list(pair = pairs[[pi]], pair_id = pi, cohort = k, transitions = transitions[pi])
        if (dated) sp$reference <- r
        specs[[length(specs) + 1L]] <- sp
      }
    }
  }
  find_support <- function(pi, k, r = NULL) which(vapply(specs, function(s)
    s$pair_id == pi && s$cohort == k && (!dated || s$reference == r), logical(1)))
  for (i in seq_along(pools)) {
    s <- pools[[i]]; pr <- pairs[[s$pair_id]]
    pools[[i]]$support_ids <- vapply(s$cohorts, function(k)
      as.integer(find_support(s$pair_id, k, s$reference)), integer(1))
    pools[[i]]$transitions <- transitions[s$pair_id]
    pools[[i]]$block_key <- if (dated) paste0("pair:", pr$s, "-", pr$t, "|pool:",
      paste(cohorts[s$cohorts], collapse = ","), "|reference-design:L_e*z(C^r)|",
      reference$recipes$key[s$reference], "|weights:W_g/r") else paste0("pair:", pr$s, "-", pr$t, "|pool:",
      paste(cohorts[s$cohorts], collapse = ","), "|", recipe$key, "|weights:w/r")
  }
  for (i in seq_along(links)) links[[i]]$support_id <- as.integer(
    find_support(links[[i]]$pair_id, links[[i]]$cohort, links[[i]]$reference))
  list(specs = specs, pools = pools, links = links)
}

cc_frozen$dr_baseline_treatment_fit <- function(reader, spec, supports, reg, shared,
                                      recipe, or_frame_index, trt_recipe,
                                      glm_control, ctrl, reference = NULL,
                                      complete = NULL, share = NULL) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- if (is.null(complete)) NULL else or_baseline_complete_panel(reader, complete, dr_baseline_abort)
  p_trt <- trt_recipe$p
  p_or <- shared$p
  kept_or_cols <- if (p_or > 0L) or_frame_index[recipe$slope_index[shared$keep]] else integer(0)
  m <- length(spec$cohorts)
  n <- spec$n_observers

  X_trt <- matrix(0, n, p_trt, dimnames = list(NULL, trt_recipe$columns))
  Z_or <- matrix(0, n, p_or, dimnames = if (p_or > 0L) list(NULL, shared$basis$kept_columns) else NULL)
  y <- w <- dy <- numeric(n)
  kk <- integer(n)
  obs_rows <- if (is.null(share)) NULL else integer(n)
  pos <- 0L

  if (!is.null(panel)) {
    # CP08 stage 4: the eligible fit's observer rows in the scan's order (every
    # frame row of the model's cohorts; complete dates observe all), with one
    # native call reading their raw pair differences under the checks and
    # errors of the scan's per-chunk or_baseline_complete_pair_rows calls; the
    # design, response, cohort and weight vectors take the scan's values. The
    # general scan below runs only without complete-date facts.
    if (p_or < 0L || !is.null(reference) || is.null(obs_rows))
      dr_baseline_abort("Complete-support treatment reads require the facts of their fit.", "integrity")
    # The frame rows of the model's cohorts in increasing order, from the
    # reader's cohort ranges (no frame-length temporary).
    rows <- sort(unlist(lapply(spec$cohorts, function(k) st$cohort_order[seq_len(st$cohort_end[k] -
      st$cohort_start[k] + 1) + st$cohort_start[k] - 1])), method = "radix")
    post <- panel$postings
    res <- .Call(or_baseline_native_kernel("C_gmm_cp08_pair_rows"), panel$observations$y, post$date_start,
      post$post_a_index, post$post_obs_index, rows, st$a_index,
      as.integer(c(spec$pair$s_calendar_index, spec$pair$t_calendar_index)), as.integer(ctrl$max_chunk_rows),
      ctrl$max_matrix_entries)
    if (!identical(res$status[1L], 0L)) or_baseline_native_refuse(res$status, dr_baseline_abort)
    pos <- length(rows)
    if (pos == n && n > 0L) {
      ci <- match(st$cohort_index[rows], spec$cohorts)
      a_idx <- st$a_index[rows]
      obs_rows[] <- rows
      X_trt[] <- frame$Z_a[a_idx, trt_recipe$index, drop = FALSE]
      if (p_or > 0L) Z_or[] <- frame$Z_a[a_idx, kept_or_cols, drop = FALSE]
      y[] <- as.double(ci == 1L)
      kk[] <- ci
      dy[] <- res$dy
      w[] <- st$w[rows] * vapply(supports[spec$support_ids], `[[`, numeric(1), "v")[ci]
    }
  }
  for (start in if (is.null(panel)) seq.int(1L, st$N, by = ctrl$max_chunk_rows)) {
    rows <- seq.int(start, min(st$N, start + ctrl$max_chunk_rows - 1L))
    ci <- match(st$cohort_index[rows], spec$cohorts)
    rows <- rows[!is.na(ci)]
    ci <- ci[!is.na(ci)]
    if (!length(rows)) next
    hit <- if (is.null(panel)) cp04_pair_rows(reader, rows, spec$pair, TRUE) else
      or_baseline_complete_pair_rows(reader, panel, rows, spec$pair, TRUE, dr_baseline_abort)
    keep <- which(hit$S == 1)
    if (!length(keep)) next
    rows <- rows[keep]
    ci <- ci[keep]
    ii <- seq.int(pos + 1L, pos + length(rows))
    pos <- pos + length(rows)
    if (!is.null(obs_rows)) obs_rows[ii] <- rows
    a_idx <- st$a_index[rows]
    if (is.null(reference)) {
      X_trt[ii, ] <- frame$Z_a[a_idx, trt_recipe$index, drop = FALSE]
      if (p_or > 0L) Z_or[ii, ] <- frame$Z_a[a_idx, kept_or_cols, drop = FALSE]
    } else {
      Xr <- chained_reference_x(reference, spec$reference, rows)
      X_trt[ii, ] <- Xr[, trt_recipe$index, drop = FALSE]
      if (p_or > 0L) Z_or[ii, ] <- Xr[, kept_or_cols, drop = FALSE]
    }
    y[ii] <- as.double(ci == 1L)
    kk[ii] <- ci
    dy[ii] <- hit$dy[keep]
    w_rows <- if (is.null(reference)) st$w[rows] else chained_reference_weight(reference, spec$reference, rows)
    for (k in unique(ci)) {
      jj <- which(ci == k)
      sid <- spec$support_ids[k]
      w[ii[jj]] <- w_rows[jj] * supports[[sid]]$v
    }
  }
  if (pos != n) dr_baseline_abort("Observed risk support changed from its owner view.", "identity")

  fit <- if (is.null(share)) cp04_logit(X_trt, y, w, glm_control, st$N, spec) else
    dr_baseline_shared_logit(share, spec, obs_rows, X_trt, y, w,
      vapply(supports[spec$support_ids], `[[`, numeric(1), "v"), glm_control, st$N)
  names(fit$coefficients) <- trt_recipe$columns
  lv <- if (is.null(share)) ipw_logit_values(X_trt, fit$coefficients) else {
    if (is.null(share$lv)) share$lv <- ipw_logit_values(X_trt, fit$coefficients)
    share$lv
  }
  wp <- w * y
  wm <- w * (1 - y) * lv$odds
  dp <- sum(wp) / st$N
  dm <- sum(wm) / st$N
  cp04_finite(c(dp, dm), "treated/control normalizers", TRUE)

  # Raw outcome IPW quantities
  T_ipw <- sum(wp * dy) / sum(wp)
  M_ipw <- sum(wm * dy) / sum(wm)
  mu_ipw <- T_ipw - M_ipw

  # Residualized DR quantities reusing fitted OR parameters (no treatment refit)
  d_obs <- reg$L * Z_or
  m_hat <- (if (recipe$intercept) reg$alpha else 0) +
    (if (p_or > 0L) as.vector(d_obs %*% shared$beta) else 0)
  r_obs <- dy - m_hat
  T_dr <- sum(wp * r_obs) / sum(wp)
  M_dr <- sum(wm * r_obs) / sum(wm)
  mu_dr <- T_dr - M_dr
  d_bar_plus <- if (p_or > 0L) as.vector(colSums(wp * d_obs) / sum(wp)) else numeric(0)
  d_bar_minus <- if (p_or > 0L) as.vector(colSums(wm * d_obs) / sum(wm)) else numeric(0)
  d_bar_diff <- d_bar_plus - d_bar_minus
  if (p_or > 0L) {
    names(d_bar_plus) <- names(d_bar_minus) <- names(d_bar_diff) <- shared$basis$kept_columns
  }
  mu_dr_decomp <- mu_ipw - (if (p_or > 0L) sum(d_bar_diff * shared$beta) else 0)
  cp04_finite(c(T_ipw, M_ipw, mu_ipw, T_dr, M_dr, mu_dr, d_bar_plus, d_bar_minus, mu_dr_decomp),
              "treated/control DR and IPW means")

  # Upstream rate cross-derivatives on spec$support_ids (log_deriv = -1 / rate_k)
  C <- matrix(0, p_trt, m)
  U_dr <- matrix(0, 4L, m, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  U_ipw <- matrix(0, 4L, m, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  score <- w * (y - lv$prob) * X_trt
  residual_dr <- cbind(wp * (r_obs - T_dr) / dp, wm * (r_obs - M_dr) / dm, wp, wm)
  residual_ipw <- cbind(wp * (dy - T_ipw) / dp, wm * (dy - M_ipw) / dm, wp, wm)

  for (k in seq_len(m)) {
    rr <- which(kk == k)
    if (!length(rr)) next
    sid <- spec$support_ids[k]
    r_k <- supports[[sid]]$rate
    scale_k <- -1 / (st$N * r_k)
    C[, k] <- colSums(score[rr, , drop = FALSE]) * scale_k
    U_dr[, k] <- colSums(residual_dr[rr, , drop = FALSE]) * scale_k
    U_ipw[, k] <- colSums(residual_ipw[rr, , drop = FALSE]) * scale_k
  }

  V_dr <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), trt_recipe$columns))
  V_dr[2L, ] <- colSums(wm * (r_obs - M_dr) * X_trt) / (st$N * dm)
  V_dr[4L, ] <- colSums(wm * X_trt) / st$N

  V_ipw <- matrix(0, 4L, p_trt, dimnames = list(c("T", "M", "dp", "dm"), trt_recipe$columns))
  V_ipw[2L, ] <- colSums(wm * (dy - M_ipw) * X_trt) / (st$N * dm)
  V_ipw[4L, ] <- V_dr[4L, ]

  cp04_finite(c(C, U_dr, U_ipw, V_dr, V_ipw), "baseline DR analytic cross derivatives")

  fit$support_ids <- as.integer(spec$support_ids)
  fit$cross <- C
  fit$partials <- list(upstream = U_dr, beta = V_dr)
  fit$partials_ipw <- list(upstream = U_ipw, beta = V_ipw)
  fit$means <- c(
    T = T_dr, M = M_dr, dp = dp, dm = dm,
    T_dr = T_dr, M_dr = M_dr, mu_dr = mu_dr,
    T_ipw = T_ipw, M_ipw = M_ipw, mu_ipw = mu_ipw,
    mu_or = T_dr, mu_dr_decomp = mu_dr_decomp
  )
  fit$d_bar_plus <- d_bar_plus
  fit$d_bar_minus <- d_bar_minus
  fit$d_bar_diff <- d_bar_diff
  fit
}

cc_frozen$pipeline_metadata_row_spec_fields <- function(metadata_obj, r) {
  ri <- as.integer(r)
  row_df <- metadata_obj$catalog$rows[ri, , drop = FALSE]
  map_df <- metadata_obj$row_mapping[ri, , drop = FALSE]
  pl <- metadata_obj$catalog$pools
  coh <- metadata_obj$catalog$spec$cohorts
  cal <- metadata_obj$catalog$spec$calendar

  pid <- row_df$pool_id
  k_indices <- pl$entries[seq.int(pl$start[pid], length.out = pl$size[pid])]
  K_vals <- coh[k_indices]

  list(
    row_id = ri,
    operator_id = row_df$operator_id,
    pool_id = pid,
    target_index = row_df$target_index,
    g = coh[row_df$target_index],
    K_indices = as.integer(k_indices),
    K = as.double(K_vals),
    s_index = row_df$s_index,
    t_index = row_df$t_index,
    s = cal[row_df$s_index],
    t = cal[row_df$t_index],
    ref_index = row_df$ref_index,
    ref = cal[row_df$ref_index],
    functional_key = map_df$functional_key,
    s_role = row_df$s_role,
    t_role = row_df$t_role,
    s_coordinate_id = row_df$s_coordinate_id,
    t_coordinate_id = row_df$t_coordinate_id,
    original = row_df$original,
    n_horizons = row_df$n_horizons,
    horizon_first_index = row_df$horizon_first_index,
    horizon_last_index = row_df$horizon_last_index,
    target_occupancy_cell_id = map_df$target_occupancy_cell_id,
    target_pair_count = map_df$target_pair_count,
    control_pool_pair_count = map_df$control_pool_pair_count,
    risk_pair_count = map_df$risk_pair_count,
    target_a_unit_count = map_df$target_a_unit_count,
    control_pool_a_unit_count = map_df$control_pool_a_unit_count,
    risk_a_unit_count = map_df$risk_a_unit_count
  )
}

cc_frozen$extract_pipeline_catalog_row_spec <- function(metadata_obj, row_id) {
  validate_owner_comparison_catalog(metadata_obj)
  r <- ipw_owner_int_bound(row_id, "row_id", min_val = 1)
  n_rows <- nrow(metadata_obj$catalog$rows)
  if (r > n_rows) {
    cdid_pipeline_metadata_abort(
      paste0("`row_id` (", r, ") exceeds retained catalog row count (", n_rows, ")."),
      c("cdid_pipeline_metadata_bound_error", "cdid_catalog_input_error")
    )
  }
  pipeline_metadata_row_spec_preflight(metadata_obj,r,metadata_obj$control$owner_control)
  pipeline_metadata_row_spec_fields(metadata_obj, r)
}

for (nm in c("or_baseline_catalog", "dr_baseline_treatment_fit",
             "pipeline_metadata_row_spec_fields", "extract_pipeline_catalog_row_spec")) {
  fn <- get(nm, envir = cc_frozen, inherits = FALSE)
  environment(fn) <- cc_frozen
  assign(nm, fn, envir = cc_frozen)
}
lockEnvironment(cc_frozen, bindings = TRUE)

# ---------------------------------------------------------------------------
# Shared fixtures and route helpers
# ---------------------------------------------------------------------------
cc_ctrl <- function(chunk = 7L, retained = 1e7, matrix = 2e7, work = 1e10, output = 1e7) {
  api$ipw_owner_control(
    max_chunk_rows = chunk,
    max_retained_entries = retained,
    max_matrix_entries = matrix,
    max_work_ops = work,
    max_output_entries = output
  )
}
cc_glm <- list(epsilon = 1e-13, maxit = 60L, trace = FALSE)
cc_icpt <- api$or_pooled_outcome_recipe("(Intercept)", TRUE, character(0))

cc_panel <- function(N = 41L, T = 5L, seed = 11L, second = FALSE, id_perm = FALSE, rows_perm = FALSE) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:2] <- c(3, Inf)
  if (second) g[3:8] <- 4
  ids <- if (id_perm) sample(1000L + seq_len(5L * N), N) else seq_len(N)
  xi <- stats::rbinom(N, 1, 0.5)
  ue <- stats::rnorm(N)
  d <- expand.grid(t = seq_len(T), id = ids, KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  u <- match(d$id, ids)
  d$g <- g[u]
  d$x <- xi[u]
  d$x2 <- cos(u * 0.31)
  d$cl <- (u - 1L) %/% 4L + 1L
  eff <- ifelse(is.finite(d$g) & d$t >= d$g, 1 + 0.25 * (d$t - d$g), 0)
  d$y <- ue[u] + 0.1 * d$t + 0.35 * sin(0.7 * u * d$t) + 0.2 * d$x * (d$t^2) + eff + stats::rnorm(nrow(d))
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

cc_d75_panel <- function(N = 41L, T = 5L, seed = 11L, second = FALSE) {
  set.seed(seed)
  g <- ifelse(stats::runif(N) < 0.45, 3, Inf)
  g[1:2] <- c(3, Inf)
  if (second) g[3:8] <- 4
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$y <- 0.2 * d$id + 0.4 * d$t + ifelse(d$t >= d$g, 1.25, 0) + stats::rnorm(nrow(d))
  d
}
cc_small_panel <- function(N = 12L, T = 5L, g_units = c(1L, 3L, 5L, 7L, 9L, 11L),
                           x_units = c(1L, 2L, 3L, 5L, 8L), seed = 55L) {
  set.seed(seed)
  g <- rep(Inf, N); g[g_units] <- 3
  x <- as.numeric(seq_len(N) %in% x_units)
  ue <- stats::rnorm(N)
  d <- data.frame(id = rep(seq_len(N), each = T), t = rep(seq_len(T), N))
  u <- d$id
  d$g <- g[u]
  d$x <- x[u]
  d$x2 <- cos(u * 0.31)
  d$cl <- (u - 1L) %/% 2L + 1L
  eff <- ifelse(is.finite(d$g) & d$t >= d$g, 1 + 0.25 * (d$t - d$g), 0)
  d$y <- ue[u] + 0.1 * d$t + 0.35 * sin(0.7 * u * d$t) + 0.2 * d$x * (d$t^2) + eff + stats::rnorm(nrow(d))
  d
}

cc_pairs <- function(T = 5L) {
  pr <- list()
  for (s in seq_len(T - 1L)) for (t in seq.int(s + 1L, T)) pr[[length(pr) + 1L]] <- c(s, t)
  pr
}

cc_meta_fixture <- function(d = cc_panel(second = TRUE), control_group = "notyettreated",
                            targets = c(3, 4), chunk = 7L, ctrl = cc_ctrl(chunk)) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  u <- d[!duplicated(d$id), c("id", "g", "x", "x2")]
  u <- u[order(u$id), ]
  Z <- cbind(`(Intercept)` = 1, x = u$x, x2 = u$x2)
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = rep(1, nrow(u)), control = ctrl,
    frame_g = u$g, design_on_A_only = TRUE)
  meta <- api$build_owner_comparison_catalog(
    po, targets = targets, anticipation = 0, control_group = control_group,
    conventions = "cp08-construction", control = ctrl,
    catalog_bounds = list(max_candidates = 1e4, max_rows = 1e4, max_pool_entries = 1e5, max_scratch = 1e7, max_work = 1e9),
    frame_owner = fo
  )
  list(d = d, po = po, fo = fo, meta = meta, ctrl = ctrl, u = u)
}

cc_od_fixture <- function(d = cc_panel(), chunk = 7L, perm = FALSE, w = NULL, X = NULL, V = NULL,
                          pairs = cc_pairs(), K = Inf, target = 3, ctrl = cc_ctrl(chunk)) {
  po <- api$construct_ipw_panel_owner(d, "id", "t", "y", "g", ctrl)
  u <- d[!duplicated(d$id), c("id", "g")]
  u <- u[order(u$id), ]
  if (perm) { set.seed(5L); u <- u[sample(nrow(u)), ] }
  N <- nrow(u)
  Z <- if (is.null(X)) matrix(1, N, 1L, dimnames = list(NULL, "(Intercept)")) else X[match(u$id, sort(u$id)), , drop = FALSE]
  V_mat <- if (is.null(V)) matrix(1, N, 1L, dimnames = list(NULL, "reader_intercept")) else V[match(u$id, sort(u$id)), , drop = FALSE]
  wt <- if (is.null(w)) rep(1, N) else w[match(u$id, sort(u$id))]
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = wt, control = ctrl, frame_g = u$g,
    allow_frame_only = FALSE, design_on_A_only = TRUE)
  reader <- api$construct_cp04_frame_reader(fo, u$id, V_mat, u$g, wt, ctrl, allow_unknown_off_A = TRUE)
  views <- lapply(seq_along(pairs), function(k) {
    p <- pairs[[k]]
    api$construct_ipw_link_view(
      fo, p[1], p[2], g = target[(k - 1L) %% length(target) + 1L], K = K, control = ctrl,
      link_id = paste0("L", p[1], "_", p[2],
        if (sum(vapply(pairs, function(q) identical(q, p), logical(1))) > 1L) paste0(":", k) else ""),
      functional_key = paste0("cmpcat:od:", p[1], "-", p[2])
    )
  })
  list(d = d, po = po, fo = fo, reader = reader, views = views, ctrl = ctrl, chunk = chunk, u = u, pairs = pairs)
}

cc_strip <- function(x) {
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return("<function>")
  if (inherits(x, "cdid_ipw_structural_certificate")) {
    x <- unclass(x)
    return(lapply(x[setdiff(names(x), c("link_kind", "object_refs", "ref"))], cc_strip))
  }
  if (is.list(x)) { a <- attributes(x); x <- lapply(x, cc_strip); attributes(x) <- a }
  x
}

cc_norm_shared <- function(s) {
  if (!is.list(s)) return(s)
  norm_blocks <- function(blks) {
    if (!is.list(blks)) return(blks)
    lapply(blks, function(bk) {
      if (is.list(bk) && !is.null(bk$r11)) {
        sgn <- if (bk$r11 < 0) -1 else 1
        bk$r11 <- abs(bk$r11)
        if (!is.null(bk$r1d)) bk$r1d <- sgn * bk$r1d
        if (!is.null(bk$r1y)) bk$r1y <- sgn * bk$r1y
      }
      bk
    })
  }
  drop_acc <- c("max_chunk_merge_rows", "max_residual_merge_rows", "n_compressions")
  if (is.list(s$fit)) {
    s$fit$blocks <- norm_blocks(s$fit$blocks)
    if (is.list(s$fit$accounting)) s$fit$accounting <- s$fit$accounting[setdiff(names(s$fit$accounting), drop_acc)]
    if (is.matrix(s$fit$R) && ncol(s$fit$R) > 0L) {
      piv <- s$fit$pivot
      R_unpiv <- s$fit$R[, order(piv), drop = FALSE]
      s$fit$R <- crossprod(R_unpiv)
      s$fit$pivot <- seq_len(ncol(s$fit$R))
      if (is.list(s$fit$slope_fit)) {
        s$fit$slope_fit$R <- s$fit$R
        s$fit$slope_fit$pivot <- s$fit$pivot
      }
    }
  }
  if (is.list(s$summary)) {
    s$summary$blocks <- norm_blocks(s$summary$blocks)
    if (is.list(s$summary$accounting)) s$summary$accounting <- s$summary$accounting[setdiff(names(s$summary$accounting), drop_acc)]
    if (is.list(s$summary$residual) && is.matrix(s$summary$residual$factor)) {
      s$summary$residual$factor <- crossprod(s$summary$residual$factor)
    }
  }
  s
}

# Route work is a resource contract, not an equal-estimator output. Recompute
# its reservation from input rows, catalog dimensions and the documented compact
# payload convention; do not call the implementation's charge/payload helpers.
cc_expect_work <- function(b, treatment_route = "batch") {
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st <- reader$state; N <- length(st$id); nA <- sum(st$a_index > 0L)
  S <- length(b$support_specs); G <- length(b$regression_specs); L <- length(b$links)
  p <- b$outcome_recipe$p; q <- p + 2; B <- min(N, b$control$max_chunk_rows)
  dr <- inherits(b, "cdid_dr_baseline")
  # Reconstruct only the pre-fit link descriptors; fitted moments are not
  # metadata charged by the fit preflight. Preserve their original field order.
  keys <- c("regression", "pair_id", "cohort", "link_id", "source_functional_key", "functional_key",
    if (!is.null(b$reference_ref)) "reference", "support_id",
    if (dr) c("model", "control_support_ids"))
  links <- lapply(b$links, function(l) l[keys])
  meta <- list(b$support_specs, b$regression_specs, links)
  if (dr) meta <- c(meta, list(b$treatment_specs, b$treatment_recipe))
  entries <- function(x, attrs = TRUE) {
    if (is.null(x) || is.environment(x)) return(0)
    z <- length(x) + 8.0
    if (is.list(x)) z <- z + sum(vapply(x, entries, numeric(1)))
    if (is.character(x)) z <- z + sum(ceiling(nchar(x[!is.na(x)], type = "bytes") / 2))
    if (attrs && length(attributes(x))) z <- z + entries(attributes(x), FALSE)
    z
  }
  metadata <- entries(meta); header <- 4096 + 1024 * L * (length(st$cohorts) + 1)
  common <- !is.null(b$complete_support)
  if (!common) {
    own <- 512 * S * N * q^2 * (1 + ceiling(log2(nA + 1))) +
      256 * S * N * q * (p + G + 2) + 256 * G * (S + 1) * q^3 +
      256 * (G + S) * q^2 * (G + 2) + 128 * (metadata + header + L * q)
  } else {
    nk <- vapply(b$support_specs, function(s) sum(st$g == st$cohorts[s$cohort]), numeric(1))
    nc <- sum(is.infinite(st$g)); rows <- sum(nk)
    dates <- unique(unlist(lapply(b$support_specs, function(s) c(s$pair$s_calendar_index, s$pair$t_calendar_index))))
    used <- b$outcome_recipe$columns
    if (dr) used <- union(used, b$treatment_recipe$columns)
    facts <- length(dates) * nA + (1 + length(used)) * N + 8 * L
    qr <- if (p > 0) 256 * nc * q^2 + 256 * ceiling(nc / B) * q^3 + 512 * q^3 + 256 * rows * p^2 else 0
    own <- 256 * (rows + G * nc) * q + qr + 256 * (G + S) * q +
      128 * (metadata + header + L * q) + 64 * facts
    # These bounded comparison fixtures do not trigger raw-factor overflow
    # validation. Leave the separate large-value validation tests intact.
    expect_equal(b$accounting$complete_support$original_factor_validation_work, 0, tolerance = 0)
  }
  if (dr) {
    Tn <- length(b$treatment_specs); pt <- b$treatment_recipe$p
    obs <- vapply(b$treatment_specs, function(s) sum(st$g %in% st$cohorts[s$cohorts]), numeric(1))
    expect_equal(unname(vapply(b$treatment_specs, `[[`, integer(1), "n_observers")), obs, tolerance = 0)
    m <- max(vapply(b$treatment_specs, function(s) length(s$cohorts), integer(1)))
    models <- Tn * (pt^2 + 12 * (pt + p) + 8 * m + pt * m + 4096) + L * (8 * p + 8 * pt + 1024) + 8192
    if (!common) {
      trt <- 128 * (b$glm_control$maxit + 2) * sum(obs) * (pt + 2)^2 + 128 * Tn * N * m + 128 * models
    } else {
      n1 <- obs[1L]; rem <- max(0L, Tn - 1L)
      trt <- 128 * (b$glm_control$maxit + 2) * n1 * (pt + 2)^2 + 128 * models
      if (treatment_route == "batch" || rem == 0L) {
        trt <- trt + 128 * n1 * ((pt + 2)^2 + p + m) + 64 * rem * n1 * (p + pt + m + 8)
      } else {
        # A deliberately declined summary batch reserves the old per-model
        # loop AND its attempted validation, before executing either.
        summary <- 12 + 3 * p + pt * m + 8 * (m + pt) + 16
        trt <- trt + 128 * sum(obs) * ((pt + 2)^2 + p + m) + 128 * Tn * N * m +
          64 * rem * n1 * (p + pt + m + 8) + 128 * rem * summary
      }
    }
    own <- own + trt
  }
  expect_equal(b$accounting$own_work_ops, own, tolerance = 0)
  expect_equal(b$accounting$work_ops_bound, b$accounting$external_work_ops + own, tolerance = 0)
  expect_identical(b$accounting$work_ops_convention, "cumulative")
}

cc_check_dr_clipping <- function(tr, info = NULL) {
  clip_names <- c("eps", "p_min", "p_max", "p_clip_max", "p_n", "p_capped", "p_capped_ctrl",
    "r_min", "r_max", "r_floor_min", "r_n", "r_floored", "w_max")
  expect_true(is.list(tr$clipping), info = info)
  expect_identical(names(tr$clipping), clip_names, info = info)
  expect_identical(tr$clipping$eps, 1e-6, info = info)
  expect_identical(as.integer(tr$clipping$r_n), length(tr$support_ids), info = info)
}

cc_strip_dr_clipping <- function(tr, info = NULL) {
  cc_check_dr_clipping(tr, info = info)
  tr[setdiff(names(tr), "clipping")]
}

cc_backend_view <- function(b, treatment_route = "batch", strip_clipping = FALSE) {
  cc_expect_work(b, treatment_route)
  x <- unclass(b)[setdiff(names(b), c("ref", "reader_ref", "complete_support", "accounting"))]
  if (!is.null(x$shared)) x$shared <- cc_norm_shared(x$shared)
  if (isTRUE(strip_clipping) && !is.null(x$treatments)) {
    x$treatments <- lapply(x$treatments, cc_strip_dr_clipping)
  }
  list(fields = cc_strip(x), accounting = b$accounting[setdiff(names(b$accounting),
    c("complete_support", "work_ops_bound", "own_work_ops", "fit_scratch_bound", "retained_bound", "owned_entries"))])
}

cc_expect_route_match <- function(a, b, tol = 1e-12, info = NULL, mixed_scalar_path = NULL) {
  mismatch <- function(x, y, path) {
    if (!identical(class(x), class(y)) || !identical(typeof(x), typeof(y)) ||
        !identical(names(x), names(y)) || !identical(dim(x), dim(y)) ||
        !identical(dimnames(x), dimnames(y))) return(paste(path, "schema"))
    if (is.list(x)) {
      if (length(x) != length(y)) return(paste(path, "length"))
      return(unlist(lapply(seq_along(x), function(i) mismatch(x[[i]], y[[i]],
        paste0(path, "$", if (is.null(names(x))) i else names(x)[i]))), use.names = FALSE))
    }
    same <- if (identical(path, mixed_scalar_path)) {
      typeof(x) == "double" && length(x) == 1L && is.finite(x) && is.finite(y) &&
        abs(x - y) <= 1e-12 * max(1, abs(x), abs(y))
    } else if (typeof(x) == "double") length(waldo::compare(x, y, tolerance = tol)) == 0L else identical(x, y)
    if (same) character() else paste(path, "value")
  }
  bad <- mismatch(a, b, "root")
  expect_true(length(bad) == 0L, info = paste(c(info, head(bad, 12L)), collapse = "\n"))
}

cc_run_pipeline <- function(d, est = "reg", chunk = 7L, reduce = TRUE, weighting = "optimal",
                            pt_assumption = "post", pt_start = NULL, xformla = NULL,
                            control_group = "nevertreated", target_cohorts = NULL,
                            bstrap = FALSE, cband = FALSE, biters = 99L, seed = 101L,
                            aggregations = c("simple", "group", "dynamic", "calendar"),
                            condition_trigger = NULL, ridge_relative = NULL,
                            structural_reduction_stage = "outcome_free", ns = api) {
  d_in <- d
  d_in$G <- ifelse(is.infinite(d_in$g), 0, d_in$g)
  d_in$period <- d_in$t
  ctrl_list <- if (is.null(chunk)) list() else list(chunk_rows = chunk)
  if (!is.null(condition_trigger)) ctrl_list$condition_trigger <- condition_trigger
  if (!is.null(ridge_relative)) ctrl_list$ridge_relative <- ridge_relative
  if (!length(ctrl_list)) ctrl_list <- NULL
  resolved <- api$chained_public_resolve_control(ctrl_list, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  prep <- api$chained_request_prepare_inputs(
    data = d_in, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xformla, weightsname = NULL, clustervars = "cl", anticipation = 0,
    frame = NULL, outcome_regression = identical(est, "reg"),
    doubly_robust = identical(est, "dr"), or_resource = resolved$profile
  )
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = biters, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = weighting, gmm_covariance = "tikhonov")
  tc <- if (is.null(target_cohorts)) prep$target_cohorts else target_cohorts
  ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = built$control, xformla = prep$xformla, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = tc, anticipation = 0L, control_group = control_group, base_period = "varying",
    weighting = weighting, reduce_structural_redundancy = reduce,
    structural_reduction_stage = structural_reduction_stage,
    run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = aggregations, cband = cband, alp = 0.05, biters = biters, seed = seed,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = bstrap, caller_rng = FALSE,
    gmm_covariance = "tikhonov", pt_assumption = pt_assumption, pt_start = pt_start, cp04 = prep$cp04
  )
}

# Namespace wrapper that forces the old frozen construction helpers (no descriptor
# source, old accumulated-list catalog, old per-model DR treatment fit) or
# declines complete-support facts altogether (`decline_complete = TRUE`).
cc_reference_ns <- function(decline_complete = FALSE, use_frozen_construction = TRUE) {
  e <- new.env(parent = api)
  for (nm in c("run_chained_ipw_pipeline", "build_pipeline_nuisance_and_link_views",
               "build_cp04_public_registry", "or_public_fit", "dr_public_fit",
               "fit_or_baseline", "fit_dr_baseline")) {
    f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
  }
  if (use_frozen_construction) {
    e$or_baseline_catalog <- cc_frozen$or_baseline_catalog
    e$dr_baseline_treatment_fit <- cc_frozen$dr_baseline_treatment_fit
    e$pipeline_metadata_row_spec_fields <- function(metadata_obj, r, source = NULL) {
      cc_frozen$pipeline_metadata_row_spec_fields(metadata_obj, r)
    }
    e$extract_pipeline_catalog_row_specs <- function(metadata_obj, rows = seq_len(nrow(metadata_obj$catalog$rows)),
                                                     retained_entries = 0, external_work_ops = 0,
                                                     control = NULL, return_accounting = FALSE) {
      specs <- lapply(rows, function(r) cc_frozen$extract_pipeline_catalog_row_spec(metadata_obj, r))
      acct <- if (length(rows)) api$pipeline_metadata_descriptor_preflight(
        metadata_obj, control, retained_entries, external_work_ops, length(rows), include_output = TRUE
      ) else NULL
      if (isTRUE(return_accounting)) list(specs = specs, accounting = acct) else specs
    }
  }
  if (decline_complete) e$or_baseline_complete_facts <- function(...) NULL
  e
}

# ---------------------------------------------------------------------------
# 1. Provenance of frozen references and descriptor extraction parity
# ---------------------------------------------------------------------------
test_that("frozen references are separate from implementation and descriptor extraction matches across orders, attributes and errors", {
  for (nm in c("or_baseline_catalog", "dr_baseline_treatment_fit",
               "pipeline_metadata_row_spec_fields", "extract_pipeline_catalog_row_spec")) {
    f_old <- get(nm, envir = cc_frozen, inherits = FALSE)
    f_new <- get(nm, envir = api, inherits = FALSE)
    expect_true(is.function(f_old) && is.function(f_new), info = nm)
    expect_identical(environment(f_old), cc_frozen, info = nm)
    expect_identical(environment(f_new), api, info = nm)
  }
  for (nm in c("or_baseline_catalog", "dr_baseline_treatment_fit", "pipeline_metadata_row_spec_fields")) {
    expect_false(identical(body(get(nm, cc_frozen)), body(get(nm, api))), info = nm)
  }
  for (nm in c("extract_pipeline_catalog_row_specs", "pipeline_metadata_descriptor_preflight",
               "pipeline_metadata_descriptor_source", "or_baseline_catalog_key",
               "dr_baseline_shared_prepared")) {
    expect_true(is.function(get(nm, envir = api, inherits = FALSE)), info = nm)
  }

  mf <- cc_meta_fixture()
  meta <- mf$meta
  M <- nrow(meta$catalog$rows)
  P <- length(meta$catalog$pools$size)
  expect_gt(M, 4L)
  expect_gt(P, 1L)

  old_all <- lapply(seq_len(M), function(r) cc_frozen$extract_pipeline_catalog_row_spec(meta, r))
  new_single <- lapply(seq_len(M), function(r) api$extract_pipeline_catalog_row_spec(meta, r))
  new_batch <- api$extract_pipeline_catalog_row_specs(meta, seq_len(M))
  new_acct <- api$extract_pipeline_catalog_row_specs(meta, seq_len(M), return_accounting = TRUE)
  expect_identical(new_single, old_all)
  expect_identical(new_batch, old_all)
  expect_identical(new_acct$specs, old_all)
  expect_identical(new_acct$accounting$n_pools, P)
  expect_identical(new_acct$accounting$n_pool_entries, length(meta$catalog$pools$entries))
  expect_identical(new_acct$accounting$n_descriptors, as.double(M))

  # Repeated, permuted, single-element, and zero-length row sequences.
  ord <- c(M, 1L, M, 2L, 1L, max(1L, M - 1L))
  expect_identical(
    api$extract_pipeline_catalog_row_specs(meta, ord),
    lapply(ord, function(r) cc_frozen$extract_pipeline_catalog_row_spec(meta, r))
  )
  expect_identical(api$extract_pipeline_catalog_row_specs(meta, integer(0)), list())
  expect_identical(
    api$extract_pipeline_catalog_row_specs(meta, integer(0), return_accounting = TRUE),
    list(specs = list(), accounting = NULL)
  )

  # Lazy pool expansion and vector sharing inside pipeline_metadata_descriptor_source.
  src <- api$pipeline_metadata_descriptor_source(meta)
  expect_identical(src$K_indices, vector("list", P))
  expect_identical(src$K, vector("list", P))
  pid1 <- meta$catalog$rows$pool_id[1L]
  s1 <- api$pipeline_metadata_row_spec_fields(meta, 1L, src)
  expect_identical(s1, old_all[[1L]])
  expect_false(is.null(src$K_indices[[pid1]]))
  expect_identical(sum(!vapply(src$K_indices, is.null, logical(1))), 1L)
  same_pool_rows <- which(meta$catalog$rows$pool_id == pid1)
  if (length(same_pool_rows) > 1L) {
    s2 <- api$pipeline_metadata_row_spec_fields(meta, same_pool_rows[2L], src)
    expect_identical(s2, old_all[[same_pool_rows[2L]]])
    expect_true(api$gmm_sexp_aliases(s1$K_indices, s2$K_indices))
    expect_true(api$gmm_sexp_aliases(s1$K, s2$K))
  }

  # Direct-vector vs one-row data-frame attribute/name parity.
  meta_attr <- unclass(meta)
  names(meta_attr$catalog$spec$cohorts) <- paste0("coh_", seq_along(meta_attr$catalog$spec$cohorts))
  names(meta_attr$catalog$spec$calendar) <- paste0("cal_", seq_along(meta_attr$catalog$spec$calendar))
  attr(meta_attr$catalog$rows$operator_id, "tag") <- "op_tag"
  attr(meta_attr$row_mapping$functional_key, "tag") <- "fk_tag"
  src_attr <- api$pipeline_metadata_descriptor_source(meta_attr)
  for (r in seq_len(M)) {
    ref_r <- cc_frozen$pipeline_metadata_row_spec_fields(meta_attr, r)
    expect_identical(api$pipeline_metadata_row_spec_fields(meta_attr, r, NULL), ref_r)
    expect_identical(api$pipeline_metadata_row_spec_fields(meta_attr, r, src_attr), ref_r)
  }

  # Invalid row bounds and preflights preserve exact failure class/message and order.
  for (bad_r in list(0L, -1L, NA_integer_, M + 1L)) {
    e_old <- tryCatch(cc_frozen$extract_pipeline_catalog_row_spec(meta, bad_r), error = identity)
    e_new1 <- tryCatch(api$extract_pipeline_catalog_row_specs(meta, c(bad_r, 1L)), error = identity)
    e_new2 <- tryCatch(api$extract_pipeline_catalog_row_specs(meta, c(1L, bad_r)), error = identity)
    expect_identical(class(e_new1), class(e_old))
    expect_identical(conditionMessage(e_new1), conditionMessage(e_old))
    expect_identical(class(e_new2), class(e_old))
    expect_identical(conditionMessage(e_new2), conditionMessage(e_old))
  }
  e_chr_old <- tryCatch(cc_frozen$extract_pipeline_catalog_row_spec(meta, "bad"), error = identity)
  e_chr_new <- tryCatch(api$extract_pipeline_catalog_row_specs(meta, c("bad", "1")), error = identity)
  expect_identical(class(e_chr_new), class(e_chr_old))
  expect_identical(conditionMessage(e_chr_new), conditionMessage(e_chr_old))
  forged_meta <- unclass(meta); forged_meta$ref <- NULL
  expect_error(api$extract_pipeline_catalog_row_specs(forged_meta), class = "cdid_pipeline_metadata_integrity_error")
})

# ---------------------------------------------------------------------------
# 2. OR catalog equality: permutation, duplicates, attributes, ordered controls,
#    dated references, S01 undated ignored reference, and support collision fallback
# ---------------------------------------------------------------------------
test_that("or_baseline_catalog equals frozen reference across permutations, duplicates, attributes, S01 undated and dated references", {
  fx <- cc_od_fixture(d = cc_panel(second = TRUE), K = c(4, Inf), target = 3)
  st <- fx$reader$state
  base_records <- lapply(seq_along(fx$views), function(ell) {
    v <- fx$views[[ell]]
    list(
      pair = v$pair,
      target = match(v$cohorts$g, st$cohorts),
      controls = match(sort(v$cohorts$K, method = "radix"), st$cohorts),
      link_id = v$link_id,
      source_functional_key = v$functional_key,
      functional_key = paste0("cmpcat:OR|", substring(v$functional_key, 8L))
    )
  })

  # 1. Standard, permuted, and duplicate-pair/pool records.
  set.seed(19L)
  perm_records <- base_records[sample.int(length(base_records))]
  dup_records <- c(base_records, lapply(base_records[1:3], function(r) {
    r$link_id <- paste0(r$link_id, "_dup"); r
  }))
  for (recs in list(base_records, perm_records, dup_records, list())) {
    expect_identical(
      api$or_baseline_catalog(recs, st$cohorts, cc_icpt, NULL),
      cc_frozen$or_baseline_catalog(recs, st$cohorts, cc_icpt, NULL)
    )
  }

  # 2. S01 exact witness: undated records (`reference = NULL`) with ignored differing
  #    `reference` values (1L vs 2L vs NULL) on the same pair and controls share ONE pool.
  r1 <- base_records[[1L]]; r1$link_id <- "L_ref1"; r1$reference <- 1L
  r2 <- base_records[[1L]]; r2$link_id <- "L_ref2"; r2$reference <- 2L
  r3 <- base_records[[1L]]; r3$link_id <- "L_ref_null"; r3$reference <- NULL
  s01_records <- list(r1, r2, r3)
  cat_s01_new <- api$or_baseline_catalog(s01_records, st$cohorts, cc_icpt, NULL)
  cat_s01_old <- cc_frozen$or_baseline_catalog(s01_records, st$cohorts, cc_icpt, NULL)
  expect_identical(cat_s01_new, cat_s01_old)
  expect_identical(length(cat_s01_new$pools), 1L)
  expect_null(cat_s01_new$pools[[1L]]$reference)

  # 3. Dated reference (`reference` non-NULL): the same records with `reference = 1L` vs `2L`
  #    stay distinct pools and distinct dated supports.
  ref_carrier <- list(recipes = data.frame(key = c("recipe_r1", "recipe_r2"), stringsAsFactors = FALSE))
  dated_records <- list(r1, r2, within(base_records[[2L]], { link_id <- "L2_r2"; reference <- 2L }),
                        within(base_records[[2L]], { link_id <- "L2_r1"; reference <- 1L }))
  cat_dated_new <- api$or_baseline_catalog(dated_records, st$cohorts, cc_icpt, ref_carrier)
  cat_dated_old <- cc_frozen$or_baseline_catalog(dated_records, st$cohorts, cc_icpt, ref_carrier)
  expect_identical(cat_dated_new, cat_dated_old)
  expect_identical(length(cat_dated_new$pools), 4L)

  # 4. Pair bucket collisions (identical s_calendar_index/t_calendar_index, distinct attributes),
  #    pool bucket collisions (identical pair/controls key, distinct control attributes),
  #    and ordered controls (`c(2L, 3L)` vs `c(3L, 2L)`).
  p_attr <- base_records[[1L]]
  p_attr$link_id <- "L_pair_attr"
  attr(p_attr$pair, "custom_tag") <- "pair_variant"
  c_attr <- base_records[[1L]]
  c_attr$link_id <- "L_ctrl_attr"
  c_attr$controls <- stats::setNames(c_attr$controls, paste0("c", seq_along(c_attr$controls)))
  c_rev <- base_records[[1L]]
  c_rev$link_id <- "L_ctrl_rev"
  c_rev$controls <- rev(c_rev$controls)
  collision_records <- list(base_records[[1L]], p_attr, c_attr, c_rev, base_records[[2L]])
  cat_col_new <- api$or_baseline_catalog(collision_records, st$cohorts, cc_icpt, NULL)
  cat_col_old <- cc_frozen$or_baseline_catalog(collision_records, st$cohorts, cc_icpt, NULL)
  expect_identical(cat_col_new, cat_col_old)
  expect_identical(length(cat_col_new$pools), 5L)

  # 5. Support-key collision fallback (`support_key` identical under `as.character()`,
  #    distinct numeric `target` values `2` vs `2 + 2 * .Machine$double.eps`):
  s_col1 <- base_records[[1L]]; s_col1$link_id <- "L_sup_1"; s_col1$target <- 2; s_col1$controls <- 1L
  s_col2 <- base_records[[1L]]; s_col2$link_id <- "L_sup_2"; s_col2$target <- 2 + 2 * .Machine$double.eps; s_col2$controls <- 1L
  expect_false(identical(s_col1$target, s_col2$target))
  expect_identical(as.character(s_col1$target), as.character(s_col2$target))
  cat_sup_new <- api$or_baseline_catalog(list(s_col1, s_col2), st$cohorts, cc_icpt, NULL)
  cat_sup_old <- cc_frozen$or_baseline_catalog(list(s_col1, s_col2), st$cohorts, cc_icpt, NULL)
  expect_identical(cat_sup_new, cat_sup_old)
  expect_identical(c(cat_sup_new$links[[1L]]$support_id, cat_sup_new$links[[2L]]$support_id), c(2L, 3L))

  # Treatment catalog constructed from the indexed OR catalog is identical too.
  tcat_new <- api$dr_build_treatment_catalog(cat_s01_new$links, cat_s01_new$pools, cat_s01_new$specs, fx$views[c(1L, 1L, 1L)])
  tcat_old <- api$dr_build_treatment_catalog(cat_s01_old$links, cat_s01_old$pools, cat_s01_old$specs, fx$views[c(1L, 1L, 1L)])
  expect_identical(tcat_new, tcat_old)
})

# ---------------------------------------------------------------------------
# 3. S02 pre-allocation resource admission: tight caps refuse before context,
#    pool expansion, or dimension memo allocation (large E + both caller controls)
# ---------------------------------------------------------------------------
test_that("tight resource caps refuse before descriptor context, pool expansion or memo allocation across large E and caller controls", {
  mf <- cc_meta_fixture()
  meta <- mf$meta
  M <- nrow(meta$catalog$rows)
  P <- length(meta$catalog$pools$size)
  E <- length(meta$catalog$pools$entries)

  # Exact preflight formulas (memo = 0, include_output = TRUE, retained = 100, prefix = 50).
  ret_in <- 100; pre_in <- 50
  pre <- api$pipeline_metadata_descriptor_preflight(
    meta, mf$ctrl, retained_entries = ret_in, external_work_ops = pre_in,
    n_descriptors = M, include_output = TRUE, memo_entries = 0
  )
  expect_identical(pre$retained_entries, 256 + 4 * P + 2 * E)
  expect_identical(pre$scratch_entries, 256 + 4 * P + 4 * E)
  expect_identical(pre$output_entries, 128 * M)
  expect_identical(pre$work_ops_bound, 512 + 128 * (M + P + E + 1) + 16 * M * P)
  live_base <- max(ret_in, pre$output_entries)

  # Instrumented environment counting context creation (new.env) and field reads.
  make_counted_md <- function() {
    e <- new.env(parent = api)
    for (nm in c("extract_pipeline_catalog_row_specs", "pipeline_metadata_descriptor_source",
                 "pipeline_metadata_descriptor_preflight", "build_cp04_public_registry")) {
      f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
    }
    cnt <- new.env(parent = emptyenv())
    cnt$env <- cnt$fields <- cnt$source <- cnt$views <- 0L
    e$new.env <- function(...) { cnt$env <- cnt$env + 1L; base::new.env(...) }
    e$pipeline_metadata_descriptor_source <- function(...) {
      cnt$source <- cnt$source + 1L
      f <- get("pipeline_metadata_descriptor_source", api); environment(f) <- e
      f(...)
    }
    e$pipeline_metadata_row_spec_fields <- function(...) {
      cnt$fields <- cnt$fields + 1L
      api$pipeline_metadata_row_spec_fields(...)
    }
    e$cp04_public_link_views <- function(...) {
      cnt$views <- cnt$views + 1L
      api$cp04_public_link_views(...)
    }
    e$cnt <- cnt
    e
  }

  # Admitted at exact edge; refused at edge - 1 before new.env or field read.
  edges <- list(
    retained = list(field = "max_retained_entries", val = live_base + pre$retained_entries),
    matrix = list(field = "max_matrix_entries", val = live_base + pre$retained_entries + pre$scratch_entries),
    work = list(field = "max_work_ops", val = pre_in + pre$work_ops_bound),
    output = list(field = "max_output_entries", val = pre$output_entries)
  )
  for (nm in names(edges)) {
    ed <- edges[[nm]]
    c_ok <- mf$ctrl; c_ok[[ed$field]] <- ed$val
    em <- make_counted_md()
    got <- em$extract_pipeline_catalog_row_specs(
      meta, seq_len(M), retained_entries = ret_in, external_work_ops = pre_in, control = c_ok
    )
    expect_length(got, M)
    expect_identical(em$cnt$env, 1L, info = nm)
    expect_identical(em$cnt$fields, M, info = nm)

    c_bad <- mf$ctrl; c_bad[[ed$field]] <- ed$val - 1
    em_bad <- make_counted_md()
    expect_error(
      em_bad$extract_pipeline_catalog_row_specs(
        meta, seq_len(M), retained_entries = ret_in, external_work_ops = pre_in, control = c_bad
      ),
      class = "cdid_ipw_owner_bound_error", info = nm
    )
    expect_identical(c(em_bad$cnt$env, em_bad$cnt$fields), c(0L, 0L), info = nm)
  }

  # Tight metadata_obj$control$owner_control refuses on row 1 preflight BEFORE descriptor_source is called.
  k1 <- meta$catalog$pools$size[meta$catalog$rows$pool_id[1L]]
  meta_tight <- meta
  meta_tight$ref <- NULL
  meta_tight$control$owner_control$max_work_ops <- 128 + 16 * k1 - 1
  meta_tight <- api$ipw_owner_seal(meta_tight, "pipeline_metadata")
  em_row1 <- make_counted_md()
  expect_error(
    em_row1$extract_pipeline_catalog_row_specs(meta_tight, seq_len(M), control = mf$ctrl),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(c(em_row1$cnt$source, em_row1$cnt$env, em_row1$cnt$fields), c(0L, 0L, 0L))

  # Large-E witness: a catalog whose pool entry dictionary E is large (e.g. E_big = 5000)
  # is refused by a cap that easily admits small E and row 1's single-pool preflight,
  # BEFORE new.env or any cohort vector slice is allocated.
  meta_bigE <- meta
  meta_bigE$ref <- NULL
  meta_bigE$catalog$pools$entries <- rep(meta$catalog$pools$entries, length.out = 5000L)
  meta_bigE <- api$ipw_owner_seal(meta_bigE, "pipeline_metadata")
  pre_bigE <- api$pipeline_metadata_descriptor_preflight(
    meta_bigE, mf$ctrl, n_descriptors = M, include_output = TRUE
  )
  expect_identical(pre_bigE$n_pool_entries, 5000L)
  expect_gt(pre_bigE$retained_entries, pre$retained_entries + 9000)
  expect_gt(pre_bigE$scratch_entries, pre$scratch_entries + 18000)
  for (field in c("max_retained_entries", "max_matrix_entries", "max_work_ops")) {
    c_big <- mf$ctrl
    c_big[[field]] <- if (field == "max_retained_entries") {
      max(pre_bigE$output_entries, 0) + pre_bigE$retained_entries - 1
    } else if (field == "max_matrix_entries") {
      max(pre_bigE$output_entries, 0) + pre_bigE$retained_entries + pre_bigE$scratch_entries - 1
    } else {
      pre_bigE$work_ops_bound - 1
    }
    # Confirm the same control admits the small-E catalog:
    expect_length(api$extract_pipeline_catalog_row_specs(meta, seq_len(M), control = c_big), M)
    # And refuses the large-E catalog before new.env or field extraction:
    em_big <- make_counted_md()
    expect_error(
      em_big$extract_pipeline_catalog_row_specs(meta_bigE, seq_len(M), control = c_big),
      class = "cdid_ipw_owner_bound_error", info = field
    )
    expect_identical(c(em_big$cnt$env, em_big$cnt$fields), c(0L, 0L), info = field)
  }

  # Caller 1 (CP04 `build_cp04_public_registry`): tight caller ctrl$owner_control refuses
  # before new.env or cp04_public_link_views; when admitted, descriptor work_ops_bound
  # is carried into external_work_ops of the subsequent fit.
  reader <- api$construct_cp04_frame_reader(mf$fo, mf$u$id, mf$fo$Z_a, mf$u$g, rep(1, nrow(mf$u)), mf$ctrl, allow_unknown_off_A = TRUE)
  ledger <- api$ipw_pipeline_request_ledger(c(3, 4), mf$po$all_times, mf$po$cohorts, 0L, "notyettreated")
  ledger <- api$ipw_pipeline_ledger_support(ledger, meta, meta$design)
  recipe_meta <- api$or_pooled_outcome_recipe(colnames(mf$fo$Z_a), TRUE, c("x", "x2"))
  cp_ctrl <- list(owner_control = mf$ctrl, glm_control = cc_glm, rank_tol = 1e-10)
  cp_tight <- cp_ctrl
  cp_tight$owner_control$max_work_ops <- pre$work_ops_bound - 1
  em_cp <- make_counted_md()
  expect_error(
    em_cp$build_cp04_public_registry(meta, mf$fo, reader, ledger, cp_tight, "notyettreated", 0, or_method = "baseline", outcome_recipe = recipe_meta),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(c(em_cp$cnt$env, em_cp$cnt$views), c(0L, 0L))
  reg_cp <- api$build_cp04_public_registry(
    meta, mf$fo, reader, ledger, cp_ctrl, "notyettreated", 0,
    external_work_ops = 25, or_method = "baseline", outcome_recipe = recipe_meta
  )
  expect_identical(reg_cp$cp04$accounting$external_work_ops, 25 + pre$work_ops_bound)

  # Caller 2 (IPW `build_pipeline_nuisance_and_link_views`): tight owner_control refuses
  # before descriptor context (new.env), dim_memo (preflight_ipw_nuisance_view), or field reads;
  # when admitted, dim_memo caches pure dimension preflights per distinct (n_A, p, row_cohorts, ctrl)
  # and construction_work$work_ops_charged includes descriptor preparation work.
  make_counted_ipw <- function() {
    e <- new.env(parent = api)
    for (nm in c("build_pipeline_nuisance_and_link_views", "pipeline_metadata_descriptor_source",
                 "pipeline_metadata_descriptor_preflight")) {
      f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
    }
    cnt <- new.env(parent = emptyenv())
    cnt$desc_env <- cnt$dim_pre <- cnt$fields <- cnt$row_pre <- 0L
    e$pipeline_metadata_row_spec_preflight <- function(...) {
      cnt$row_pre <- cnt$row_pre + 1L
      api$pipeline_metadata_row_spec_preflight(...)
    }
    e$pipeline_metadata_descriptor_source <- function(...) {
      res <- api$pipeline_metadata_descriptor_source(...)
      cnt$desc_env <- cnt$desc_env + 1L
      res
    }
    e$preflight_ipw_nuisance_view <- function(...) {
      cnt$dim_pre <- cnt$dim_pre + 1L
      api$preflight_ipw_nuisance_view(...)
    }
    e$pipeline_metadata_row_spec_fields <- function(...) {
      cnt$fields <- cnt$fields + 1L
      api$pipeline_metadata_row_spec_fields(...)
    }
    e$cnt <- cnt
    e
  }
  fit_c <- api$ipw_fitted_view_control(max_chunk_rows = 7L, max_observer_units = 1e4,
    max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
  pt_c <- api$ipw_point_view_control(max_chunk_rows = 7L, max_observer_units = 1e4,
    max_retained_entries = 1e7, max_matrix_entries = 1e7, max_work_ops = 1e10, max_output_entries = 1e7)
  pre_ipw <- api$pipeline_metadata_descriptor_preflight(meta, mf$ctrl, n_descriptors = M, include_output = FALSE, memo_entries = 64)
  expect_identical(pre_ipw$retained_entries, 256 + 64 + 4 * P + 2 * E)

  ipw_low_work <- mf$ctrl; ipw_low_work$max_work_ops <- pre_ipw$work_ops_bound - 1
  ei_bad <- make_counted_ipw()
  expect_error(
    ei_bad$build_pipeline_nuisance_and_link_views(meta, mf$fo, ipw_low_work, fit_c, pt_c, stats::glm.control(), fit_point_links = FALSE),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_identical(c(ei_bad$cnt$desc_env, ei_bad$cnt$dim_pre, ei_bad$cnt$fields), c(0L, 0L, 0L))

  ei_ok <- make_counted_ipw()
  reg_ipw <- ei_ok$build_pipeline_nuisance_and_link_views(meta, mf$fo, mf$ctrl, fit_c, pt_c, stats::glm.control(), fit_point_links = FALSE)
  expect_identical(ei_ok$cnt$desc_env, 1L)
  expect_identical(ei_ok$cnt$row_pre, 2L * M)
  expect_identical(ei_ok$cnt$fields, M)
  # Consecutive rows with identical pool size share dim_memo instead of recomputing preflight_ipw_nuisance_view M times.
  pool_sizes <- meta$catalog$pools$size[meta$catalog$rows$pool_id]
  n_transitions <- sum( c(TRUE, diff(pool_sizes) != 0) )
  expect_identical(ei_ok$cnt$dim_pre, as.integer(n_transitions))
  expect_lt(ei_ok$cnt$dim_pre, M)
  expect_gte(reg_ipw$construction_work$work_ops_charged, pre_ipw$work_ops_bound)

  # DR fit scratch cap reserves the prepared invariants (`max_trt_n * (p_or + 4) + p_trt * (max_m + 1) + max_m + 64`)
  # before any state or treatment model runs.
  fx_dr <- cc_od_fixture()
  fit_dr <- api$fit_dr_baseline(fx_dr$reader, fx_dr$views, fx_dr$ctrl, glm_control = cc_glm, rank_tol = 1e-10, outcome_recipe = cc_icpt)
  edge_mat <- fit_dr$accounting$borrowed_entries + fit_dr$accounting$retained_bound -
    fit_dr$accounting$borrowed_entries + fit_dr$accounting$fit_scratch_bound
  expect_silent(api$fit_dr_baseline(fx_dr$reader, fx_dr$views, cc_ctrl(7L, matrix = edge_mat),
    glm_control = cc_glm, rank_tol = 1e-10, outcome_recipe = cc_icpt))
  expect_error(
    api$fit_dr_baseline(fx_dr$reader, fx_dr$views, cc_ctrl(7L, matrix = edge_mat - 1),
      glm_control = cc_glm, rank_tol = 1e-10, outcome_recipe = cc_icpt),
    class = "cdid_ipw_owner_bound_error"
  )
})

# ---------------------------------------------------------------------------
# 4. Paired prepared vs frozen-reference vs forced-generic full results for
#    IPW, OR, and DR (with and without X) at 1e-12, including ATT, full-F IF,
#    V, SE, global nuisance, active ridge full E'E metric, aux IDs, and D75
# ---------------------------------------------------------------------------
test_that("prepared construction equals frozen and forced-generic routes for IPW, OR, and DR with and without X at 1e-12", {
  dp <- cc_panel(41L, 5L, seed = 41L, id_perm = TRUE, rows_perm = TRUE)
  ns_froz <- cc_reference_ns(decline_complete = FALSE, use_frozen_construction = TRUE)
  ns_gen <- cc_reference_ns(decline_complete = TRUE, use_frozen_construction = TRUE)

  # A. IPW pipeline (no-X and with-X, plus active Tikhonov ridge):
  for (xf in list(~1, ~x)) {
    ipw_prep <- cc_run_pipeline(dp, est = "ipw", chunk = 7L, xformla = xf, structural_reduction_stage = "outcome_free", ns = api)
    ipw_gen <- cc_run_pipeline(dp, est = "ipw", chunk = 7L, xformla = xf, structural_reduction_stage = "post_point_fit", ns = ns_froz)
    expect_equal(ipw_prep$reporting$estimates, ipw_gen$reporting$estimates, tolerance = 1e-12)
    expect_equal(ipw_prep$reporting$se, ipw_gen$reporting$se, tolerance = 1e-12)
    expect_equal(ipw_prep$reporting$analytical_se, ipw_gen$reporting$analytical_se, tolerance = 1e-12)
    expect_equal(ipw_prep$reporting$if_full, ipw_gen$reporting$if_full, tolerance = 1e-12)
    expect_equal(ipw_prep$reporting$vcov, ipw_gen$reporting$vcov, tolerance = 1e-12)
    expect_identical(ipw_prep$reporting$target_meta, ipw_gen$reporting$target_meta)
    expect_identical(cc_strip(ipw_prep$auxiliary_reporting), cc_strip(ipw_gen$auxiliary_reporting))
    for (ag in c("simple", "group", "dynamic", "calendar")) {
      expect_equal(ipw_prep$aggregations[[ag]]$overall$estimate, ipw_gen$aggregations[[ag]]$overall$estimate, tolerance = 1e-12)
      expect_equal(ipw_prep$aggregations[[ag]]$overall$se, ipw_gen$aggregations[[ag]]$overall$se, tolerance = 1e-12)
    }
  }
  d_small <- cc_small_panel()
  ipw_r_prep <- cc_run_pipeline(d_small, est = "ipw", chunk = 7L, xformla = ~1, condition_trigger = 0.99, ridge_relative = 0.5,
    structural_reduction_stage = "outcome_free", ns = api)
  ipw_r_gen <- cc_run_pipeline(d_small, est = "ipw", chunk = 7L, xformla = ~1, condition_trigger = 0.99, ridge_relative = 0.5,
    structural_reduction_stage = "post_point_fit", ns = ns_froz)
  expect_identical(ipw_r_prep$covariance_handling$handler_used, "tikhonov")
  expect_equal(ipw_r_prep$covariance_handling$trace_B, ipw_r_gen$covariance_handling$trace_B, tolerance = 1e-12)
  expect_equal(ipw_r_prep$covariance_handling$scale, ipw_r_gen$covariance_handling$scale, tolerance = 1e-12)
  expect_equal(ipw_r_prep$covariance_handling$lambda, ipw_r_gen$covariance_handling$lambda, tolerance = 1e-12)
  expect_equal(ipw_r_prep$reporting$estimates, ipw_r_gen$reporting$estimates, tolerance = 1e-12)
  expect_equal(ipw_r_prep$reporting$se, ipw_r_gen$reporting$se, tolerance = 1e-12)

  # B. OR and DR pipelines (no-X, ~x, ~x + x2) and active-ridge full-metric checks:
  expect_measured_m_dr_structure <- function(x, y) {
    expect_true(typeof(x) == "double" && typeof(y) == "double" &&
      length(x) == 1L && length(y) == 1L &&
      identical(attributes(x), attributes(y), num.eq = FALSE) &&
      all(is.finite(x)) && all(is.finite(y)) &&
      identical(x == 0, y == 0) && identical(sign(x), sign(y)) &&
      identical(1 / x[x == 0], 1 / y[y == 0]))
  }
  for (xf in list(NULL, ~x, ~x + x2)) {
    for (est in c("reg", "dr")) {
      info <- paste(est, if (is.null(xf)) "noX" else deparse(xf))
      res_prep <- cc_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, reduce = TRUE, ns = api)
      res_froz <- cc_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, reduce = FALSE, ns = ns_froz)
      res_gen <- cc_run_pipeline(dp, est = est, chunk = 7L, xformla = xf, reduce = TRUE, ns = ns_gen)

      b_prep <- res_prep$nuisance_registry$cp04
      b_froz <- res_froz$nuisance_registry$cp04
      b_gen <- res_gen$nuisance_registry$cp04
      expect_false(is.null(b_prep$complete_support), info = info)
      expect_null(b_gen$complete_support, info = info)

      # Prepared and frozen complete-support backends retain expect_identical semantics
      # on pre-D93 fields except measured dr ~x + x2 treatment 1-10 partials upstream/beta
      # and means, link 1-11 mu_dr/v_ps/T_dr (plus duplicate link mu aliases at positions
      # 10/18 coupled exactly to their own canonical mu_dr), and link 3/7 M_dr (elementwise
      # mixed 1e-12; exact schema and structural/signed zeros). D93 clipping is checked
      # explicitly.
      dr_prep_froz_leaf_ok <- function(actual, expected, entries, allow_tiny_sign = FALSE, tol = 1e-12) {
        identical(typeof(actual), "double") &&
          identical(typeof(expected), "double") &&
          identical(class(actual), class(expected)) &&
          identical(is.object(actual), is.object(expected)) &&
          identical(attributes(actual), attributes(expected), num.eq = FALSE) &&
          identical(length(actual), entries) &&
          identical(length(expected), entries) &&
          all(is.finite(actual)) &&
          all(is.finite(expected)) &&
          identical(actual == 0, expected == 0) &&
          identical(1 / actual[actual == 0], 1 / expected[expected == 0]) &&
          (if (allow_tiny_sign) {
            keep <- abs(actual) > tol | abs(expected) > tol
            identical(sign(actual)[keep], sign(expected)[keep])
          } else {
            identical(sign(actual), sign(expected))
          }) &&
          all(abs(actual - expected) <= tol * pmax(1, abs(actual), abs(expected)))
      }
      dr_prep_froz_has_field <- function(x, y, nm) {
        is.list(x) && is.list(y) && !is.environment(x) && !is.environment(y) &&
          !is.object(x) && !is.object(y) &&
          identical(class(x), class(y)) &&
          identical(attributes(x), attributes(y), num.eq = FALSE) &&
          identical(length(x), length(y)) &&
          identical(which(names(x) == nm), which(names(y) == nm)) &&
          identical(length(which(names(x) == nm)), 1L)
      }
      dr_prep_froz_treatments_ok <- function(trt_a, trt_b, tol = 1e-12) {
        if (!is.list(trt_a) || !is.list(trt_b) || is.environment(trt_a) || is.environment(trt_b) ||
            is.object(trt_a) || is.object(trt_b) ||
            !identical(class(trt_a), class(trt_b)) ||
            !identical(attributes(trt_a), attributes(trt_b), num.eq = FALSE) ||
            !identical(length(trt_a), length(trt_b)) ||
            !identical(length(trt_a), 10L)) {
          return(FALSE)
        }
        a_rest <- trt_a
        b_rest <- trt_b
        for (idx in c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L)) {
          ta <- a_rest[[idx]]
          tb <- b_rest[[idx]]
          if (!dr_prep_froz_has_field(ta, tb, "partials") ||
              !dr_prep_froz_has_field(ta, tb, "means")) {
            return(FALSE)
          }
          pa <- ta[["partials"]]
          pb <- tb[["partials"]]
          if (!dr_prep_froz_has_field(pa, pb, "upstream") ||
              !dr_prep_froz_has_field(pa, pb, "beta") ||
              !dr_prep_froz_leaf_ok(pa[["upstream"]], pb[["upstream"]], 8L, allow_tiny_sign = TRUE, tol = tol) ||
              !dr_prep_froz_leaf_ok(pa[["beta"]], pb[["beta"]], 12L, allow_tiny_sign = TRUE, tol = tol) ||
              !dr_prep_froz_leaf_ok(ta[["means"]], tb[["means"]], 12L, allow_tiny_sign = FALSE, tol = tol)) {
            return(FALSE)
          }
          pa[["upstream"]] <- 0
          pb[["upstream"]] <- 0
          pa[["beta"]] <- 0
          pb[["beta"]] <- 0
          ta[["partials"]] <- pa
          tb[["partials"]] <- pb
          ta[["means"]] <- 0
          tb[["means"]] <- 0
          a_rest[[idx]] <- ta
          b_rest[[idx]] <- tb
        }
        identical(a_rest, b_rest)
      }
      dr_prep_froz_backend_ok <- function(view_a, view_b, tol = 1e-12) {
        if (!dr_prep_froz_has_field(view_a, view_b, "fields")) return(FALSE)
        fa <- view_a[["fields"]]
        fb <- view_b[["fields"]]
        if (!dr_prep_froz_has_field(fa, fb, "treatments") ||
            !dr_prep_froz_has_field(fa, fb, "links") ||
            !dr_prep_froz_treatments_ok(fa[["treatments"]], fb[["treatments"]], tol = tol)) {
          return(FALSE)
        }
        la <- fa[["links"]]
        lb <- fb[["links"]]
        if (!is.list(la) || !is.list(lb) || is.environment(la) || is.environment(lb) ||
            is.object(la) || is.object(lb) ||
            !identical(class(la), class(lb)) ||
            !identical(attributes(la), attributes(lb), num.eq = FALSE) ||
            !identical(length(la), length(lb)) ||
            !identical(length(la), 11L)) {
          return(FALSE)
        }
        for (idx in c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L)) {
          lka <- la[[idx]]
          lkb <- lb[[idx]]
          # Duplicate link mu entries at positions 10/18 must exactly equal their own
          # canonical mu_dr (raw numeric/signed zero, type, class, attributes) before
          # propagating the approved canonical mu_dr cross-route comparison via local masking.
          if (!dr_prep_froz_has_field(lka, lkb, "mu_dr") ||
              !dr_prep_froz_has_field(lka, lkb, "v_ps") ||
              !identical(which(names(lka) == "mu"), c(10L, 18L)) ||
              !identical(which(names(lkb) == "mu"), c(10L, 18L)) ||
              !identical(lka[[10L]], lka[["mu_dr"]], num.eq = FALSE) ||
              !identical(lka[[18L]], lka[["mu_dr"]], num.eq = FALSE) ||
              !identical(lkb[[10L]], lkb[["mu_dr"]], num.eq = FALSE) ||
              !identical(lkb[[18L]], lkb[["mu_dr"]], num.eq = FALSE) ||
              !dr_prep_froz_leaf_ok(lka[["mu_dr"]], lkb[["mu_dr"]], 1L, allow_tiny_sign = FALSE, tol = tol) ||
              !dr_prep_froz_leaf_ok(lka[["v_ps"]], lkb[["v_ps"]], 3L, allow_tiny_sign = TRUE, tol = tol)) {
            return(FALSE)
          }
          lka[[10L]] <- 0
          lkb[[10L]] <- 0
          lka[[18L]] <- 0
          lkb[[18L]] <- 0
          lka[["mu_dr"]] <- 0
          lkb[["mu_dr"]] <- 0
          lka[["v_ps"]] <- 0
          lkb[["v_ps"]] <- 0
          if (idx %in% c(3L, 7L)) {
            if (!dr_prep_froz_has_field(lka, lkb, "M_dr") ||
                !dr_prep_froz_leaf_ok(lka[["M_dr"]], lkb[["M_dr"]], 1L, allow_tiny_sign = FALSE, tol = tol)) {
              return(FALSE)
            }
            lka[["M_dr"]] <- 0
            lkb[["M_dr"]] <- 0
          }
          if (!dr_prep_froz_has_field(lka, lkb, "T_dr") ||
              !dr_prep_froz_leaf_ok(lka[["T_dr"]], lkb[["T_dr"]], 1L, allow_tiny_sign = FALSE, tol = tol)) {
            return(FALSE)
          }
          lka[["T_dr"]] <- 0
          lkb[["T_dr"]] <- 0
          la[[idx]] <- lka
          lb[[idx]] <- lkb
        }
        fa[["treatments"]] <- 0
        fb[["treatments"]] <- 0
        fa[["links"]] <- la
        fb[["links"]] <- lb
        va_rest <- view_a
        vb_rest <- view_b
        va_rest[["fields"]] <- fa
        vb_rest[["fields"]] <- fb
        identical(va_rest, vb_rest)
      }
      if (identical(info, "dr ~x + x2")) {
        expect_true(dr_prep_froz_backend_ok(
          cc_backend_view(b_prep, strip_clipping = TRUE),
          cc_backend_view(b_froz, "declined")
        ), info = info)
      } else {
        expect_identical(cc_backend_view(b_prep, strip_clipping = TRUE), cc_backend_view(b_froz, "declined"), info = info)
      }
      if (est == "dr") {
        if (identical(info, "dr ~x + x2")) {
          expect_true(dr_prep_froz_treatments_ok(
            lapply(b_prep$treatments, cc_strip_dr_clipping, info = info),
            b_froz$treatments
          ), info = info)
        } else {
          expect_identical(lapply(b_prep$treatments, cc_strip_dr_clipping, info = info), b_froz$treatments, info = info)
        }
      }
      if (identical(info, "dr ~x + x2"))
        expect_measured_m_dr_structure(b_prep$links[[9L]]$M_dr, b_gen$links[[9L]]$M_dr)
      cc_expect_route_match(cc_backend_view(b_prep, strip_clipping = TRUE), cc_backend_view(b_gen), tol = 1e-12, info = info,
        mixed_scalar_path = if (identical(info, "dr ~x + x2")) "root$fields$links$9$M_dr" else NULL)

      # Auxiliary IDs, target metadata, and global nuisance quantities:
      expect_identical(res_prep$nuisance_registry$cp04_aux_index, res_froz$nuisance_registry$cp04_aux_index, info = info)
      expect_identical(res_prep$nuisance_registry$cp04_aux_index, res_gen$nuisance_registry$cp04_aux_index, info = info)
      aux_prep <- cc_strip(res_prep$auxiliary_reporting)
      aux_froz <- cc_strip(res_froz$auxiliary_reporting)
      if (est == "dr") {
        owned_diff <- b_prep$accounting$owned_entries - b_froz$accounting$owned_entries
        expect_equal(owned_diff, 222 * length(b_prep$treatments), tolerance = 0, info = info)
        expect_equal(aux_prep$preflight$producer_workspace_entries - aux_froz$preflight$producer_workspace_entries,
          owned_diff, tolerance = 0, info = info)
        expect_equal(aux_prep$preflight$peak_concurrent_entries - aux_froz$preflight$peak_concurrent_entries,
          owned_diff, tolerance = 0, info = info)
        aux_prep$preflight$producer_workspace_entries <- aux_froz$preflight$producer_workspace_entries
        aux_prep$preflight$peak_concurrent_entries <- aux_froz$preflight$peak_concurrent_entries
      }
      expect_identical(aux_prep, aux_froz, info = info)
      expect_equal(res_prep$auxiliary_reporting$estimates, res_gen$auxiliary_reporting$estimates, tolerance = 1e-12, info = info)
      expect_identical(res_prep$reporting$target_meta, res_gen$reporting$target_meta, info = info)
      expect_equal(b_prep$shared$beta, b_gen$shared$beta, tolerance = 1e-12, info = info)
      expect_equal(b_prep$shared$objective, b_gen$shared$objective, tolerance = 1e-12, info = info)
      if (est == "dr") {
        expect_equal(b_prep$treatments[[1L]]$coefficients, b_gen$treatments[[1L]]$coefficients, tolerance = 1e-12, info = info)
      }

      # Public ATT, full-F IF, V, SE, and all 4 aggregations match both references at 1e-12:
      for (other in list(res_froz, res_gen)) {
        expect_equal(res_prep$reporting$estimates, other$reporting$estimates, tolerance = 1e-12, info = info)
        expect_equal(res_prep$reporting$se, other$reporting$se, tolerance = 1e-12, info = info)
        expect_equal(res_prep$reporting$analytical_se, other$reporting$analytical_se, tolerance = 1e-12, info = info)
        expect_equal(res_prep$reporting$if_full, other$reporting$if_full, tolerance = 1e-12, info = info)
        expect_equal(res_prep$reporting$vcov, other$reporting$vcov, tolerance = 1e-12, info = info)
        for (ag in c("simple", "group", "dynamic", "calendar")) {
          expect_equal(res_prep$aggregations[[ag]]$overall$estimate, other$aggregations[[ag]]$overall$estimate, tolerance = 1e-12, info = info)
          expect_equal(res_prep$aggregations[[ag]]$overall$se, other$aggregations[[ag]]$overall$se, tolerance = 1e-12, info = info)
          expect_equal(res_prep$aggregations[[ag]]$estimates, other$aggregations[[ag]]$estimates, tolerance = 1e-12, info = info)
          expect_equal(res_prep$aggregations[[ag]]$se, other$aggregations[[ag]]$se, tolerance = 1e-12, info = info)
        }
      }
    }
  }

  # C. Active Tikhonov ridge full-metric E'E verification for OR and DR (no-X and ~x):
  for (xf in list(NULL, ~x)) {
    for (est in c("reg", "dr")) {
      r_prep <- cc_run_pipeline(d_small, est = est, chunk = 7L, xformla = xf, reduce = TRUE,
        condition_trigger = 0.99, ridge_relative = 0.5, ns = api)
      r_gen <- cc_run_pipeline(d_small, est = est, chunk = 7L, xformla = xf, reduce = FALSE,
        condition_trigger = 0.99, ridge_relative = 0.5, ns = ns_gen)
      ch_p <- r_prep$covariance_handling
      ch_g <- r_gen$covariance_handling
      expect_identical(ch_p$handler_used, "tikhonov")
      expect_identical(ch_g$handler_used, "tikhonov")
      expect_identical(ch_p$dimensions$M_full, 10L)
      expect_identical(ch_p$dimensions$L_retained, 4L)
      expect_identical(ch_p$gram_provenance, "sparse_structural_reduction_E_t_E")
      expect_equal(ch_p$trace_B, ch_g$trace_B, tolerance = 1e-12)
      expect_equal(ch_p$scale, ch_g$scale, tolerance = 1e-12)
      expect_equal(ch_p$lambda, ch_g$lambda, tolerance = 1e-12)
      expect_equal(r_prep$reporting$estimates, r_gen$reporting$estimates, tolerance = 1e-12)
      expect_equal(r_prep$reporting$se, r_gen$reporting$se, tolerance = 1e-12)
      expect_equal(r_prep$reporting$if_full, r_gen$reporting$if_full, tolerance = 1e-12)
    }
  }

  # D. D75's three approved mixed-1e-12 internal M_dr scalar witnesses on the seed-23 common-design fixture:
  # This exact frozen-construction predicate predates the native DR summary
  # batch. Compare its original per-model route on both sides; fc below stays
  # the genuine batch fit for every existing D75/generic scientific check.
  ns_d75_per_model <- new.env(parent = api)
  f_d75 <- get("fit_dr_baseline", api)
  environment(f_d75) <- ns_d75_per_model
  ns_d75_per_model$fit_dr_baseline <- f_d75
  ns_d75_per_model$dr_baseline_treatment_batch <- function(...) NULL
  d23 <- cc_d75_panel(N = 41L, T = 5L, seed = 23L)
  pair_sets <- list(
    all_pairs = cc_pairs(5L),
    repeated_L = list(c(1, 3), c(2, 4), c(3, 5), c(1, 2), c(2, 3), c(1, 3)),
    single_block = list(c(2, 4))
  )
  for (p in c(1L, 2L)) {
    set.seed(19L + p)
    x1 <- sin(seq_len(41L) * 0.37) + 0.25 * ((seq_len(41L) %% 3L) - 1)
    if (p == 1L) {
      X <- cbind(`(Intercept)` = 1, x1 = x1)
      rec <- api$or_pooled_outcome_recipe(c("(Intercept)", "x1"), TRUE, "x1")
      trt_cols <- c("(Intercept)", "x1")
    } else {
      x2 <- cos(seq_len(41L) * 0.53) - 0.15 * ((seq_len(41L) %% 5L) - 2)
      X <- cbind(`(Intercept)` = 1, x1 = x1, x2 = x2)
      rec <- api$or_pooled_outcome_recipe(c("(Intercept)", "x1", "x2"), TRUE, c("x1", "x2"))
      trt_cols <- c("(Intercept)", "x1", "x2")
    }
    for (pname in names(pair_sets)) {
      fx <- cc_od_fixture(d = d23, chunk = 7L, perm = (p == 2L), X = X, V = X, pairs = pair_sets[[pname]])
      fc <- api$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = cc_glm, rank_tol = 1e-10,
        outcome_recipe = rec, treatment_columns = trt_cols)
      ff <- ns_froz$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = cc_glm, rank_tol = 1e-10,
        outcome_recipe = rec, treatment_columns = trt_cols)
      fg <- ns_gen$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = cc_glm, rank_tol = 1e-10,
        outcome_recipe = rec, treatment_columns = trt_cols)
      fc_per_model <- ns_d75_per_model$fit_dr_baseline(fx$reader, fx$views, fx$ctrl,
        glm_control = cc_glm, rank_tol = 1e-10, outcome_recipe = rec, treatment_columns = trt_cols)
      expect_identical(fc_per_model$accounting$complete_support$n_treatment_summary_batches, 0L)
      expect_identical(cc_backend_view(fc_per_model, "declined", strip_clipping = TRUE), cc_backend_view(ff, "declined"))
      approved <- list(
        "1:repeated_L" = list(index = 3L, link_id = "L3_5"),
        "1:single_block" = list(index = 1L, link_id = "L2_4"),
        "2:all_pairs" = list(index = 10L, link_id = "L4_5")
      )[[paste(p, pname, sep = ":")]]
      if (!is.null(approved)) {
        expect_identical(fc$links[[approved$index]]$link_id, approved$link_id)
        expect_identical(fg$links[[approved$index]]$link_id, approved$link_id)
      }
      mixed_path <- if (!is.null(approved)) paste0("root$fields$links$", approved$index, "$M_dr") else NULL
      cc_expect_route_match(cc_backend_view(fc, strip_clipping = TRUE), cc_backend_view(fg), tol = 1e-12, mixed_scalar_path = mixed_path)
    }
  }
})

# ---------------------------------------------------------------------------
# 5. T = 20 (190 pairs / 380 supports) DR invariant preparation:
#    1 preparation, 1 fitted logit, 189 reuses, and all 190 pair contributions
# ---------------------------------------------------------------------------
test_that("T = 20 DR prepares outcome-independent invariants once, fits one logit with 189 reuses, and retains all 190 pair records", {
  d20 <- cc_panel(N = 45L, T = 20L, seed = 77L)
  for (xf in list(NULL, ~x)) {
    e <- new.env(parent = api)
    for (nm in c("run_chained_ipw_pipeline", "build_cp04_public_registry", "dr_public_fit",
                 "fit_dr_baseline", "dr_baseline_treatment_fit", "dr_baseline_shared_logit",
                 "dr_baseline_shared_prepared", "dr_baseline_treatment_batch")) {
      f <- get(nm, api); environment(f) <- e; assign(nm, f, e)
    }
    cnt <- new.env(parent = emptyenv())
    cnt$logit <- cnt$lv <- cnt$prep_calls <- cnt$prep_hits <- cnt$trt_fits <- 0L
    cnt$batches <- cnt$batched_models <- 0L
    cnt$shares <- list()
    e$cp04_logit <- function(...) { cnt$logit <- cnt$logit + 1L; api$cp04_logit(...) }
    e$ipw_logit_values <- function(...) { cnt$lv <- cnt$lv + 1L; api$ipw_logit_values(...) }
    e$dr_baseline_shared_prepared <- function(share, ...) {
      cnt$prep_calls <- cnt$prep_calls + 1L
      if (is.environment(share) && !any(vapply(cnt$shares, identical, logical(1), share))) {
        cnt$shares[[length(cnt$shares) + 1L]] <- share
      }
      ans <- api$dr_baseline_shared_prepared(share, ...)
      if (!is.null(ans)) cnt$prep_hits <- cnt$prep_hits + 1L
      ans
    }
    e$dr_baseline_treatment_batch <- function(reader, specs, ...) {
      cnt$batches <- cnt$batches + 1L
      cnt$batched_models <- cnt$batched_models + length(specs)
      f <- get("dr_baseline_treatment_batch", api); environment(f) <- e
      f(reader, specs, ...)
    }
    e$dr_baseline_treatment_fit <- function(...) {
      cnt$trt_fits <- cnt$trt_fits + 1L
      f <- get("dr_baseline_treatment_fit", api); environment(f) <- e
      f(...)
    }

    res20 <- cc_run_pipeline(d20, est = "dr", chunk = 50L, reduce = TRUE, xformla = xf, aggregations = "simple", ns = e)
    b <- res20$nuisance_registry$cp04
    expect_identical(length(b$treatments), 190L)
    expect_identical(length(b$regressions), 190L)
    n_aux <- length(res20$nuisance_registry$cp04_aux$descriptors)
    expect_identical(n_aux, 1L)
    expect_identical(length(b$links), 190L + n_aux)
    expect_identical(vapply(b$links, `[[`, character(1), "link_id"),
      c(paste0("cp04:catalog:", seq_len(190L)), paste0("cp04:auxiliary:", seq_len(n_aux))))
    expect_identical(length(b$supports), 380L)
    # One scalar seed fit, then one native summary batch for the other pairs.
    expect_identical(cnt$trt_fits, 1L)
    expect_identical(cnt$batches, 1L)
    expect_identical(cnt$batched_models, 189L)
    expect_identical(cnt$logit, 1L)
    expect_identical(cnt$lv, 1L)
    expect_identical(cnt$prep_calls, 190L)
    expect_identical(cnt$prep_hits, 189L)
    expect_identical(length(cnt$shares), 1L)
    sh <- cnt$shares[[1L]]
    expect_identical(c(sh$calls, sh$reuses, sh$prep_calls, sh$prep_reuses), c(1L, 189L, 1L, 189L))
    expect_identical(
      c(b$accounting$complete_support$n_treatment_logit_fits,
        b$accounting$complete_support$n_treatment_logit_reuses,
        b$accounting$complete_support$n_treatment_models),
      c(1L, 189L, 190L)
    )
    # Every one of the 190 treatment models matches the frozen un-prepared reference bitwise.
    frozen_share <- new.env(parent = emptyenv())
    froz_treatments <- lapply(b$treatment_specs, function(s) {
      cc_frozen$dr_baseline_treatment_fit(
        res20$nuisance_registry$cp04_reader, s, b$supports, b$regressions[[s$regression]],
        b$shared, b$outcome_recipe, b$or_frame_index, b$treatment_recipe,
        b$glm_control, b$control, NULL, b$complete_support, frozen_share
      )
    })
    expect_identical(lapply(b$treatments, cc_strip_dr_clipping), froz_treatments)
    expect_identical(c(frozen_share$calls, frozen_share$reuses), c(1L, 189L))
    expect_identical(length(unique(lapply(b$treatments, `[[`, "support_ids"))), 190L)
    expect_null(b[["share", exact = TRUE]])
    expect_null(b[["prep", exact = TRUE]])
  }
})

# ---------------------------------------------------------------------------
# 6. Malformed, equal, and unequal DR bindings preserve decline and refusal
# ---------------------------------------------------------------------------
test_that("DR shared preparation declines on malformed or unequal bindings and preserves shared-logit refusals", {
  d <- cc_panel(N = 41L, T = 5L, seed = 17L)
  u <- d[!duplicated(d$id), c("id", "g", "x", "x2")]; u <- u[order(u$id), ]
  X <- cbind(`(Intercept)` = 1, x = u$x, x2 = u$x2)
  rec1 <- api$or_pooled_outcome_recipe(c("(Intercept)", "x"), TRUE, "x")
  fx <- cc_od_fixture(d = d, chunk = 7L, X = X, V = X)
  b <- api$fit_dr_baseline(fx$reader, fx$views, fx$ctrl, glm_control = cc_glm, rank_tol = 1e-10,
    outcome_recipe = rec1, treatment_columns = c("(Intercept)", "x"))

  reader <- fx$reader
  st <- reader$state
  ts1 <- b$treatment_specs[[1L]]
  ts2 <- b$treatment_specs[[2L]]
  kept_or <- b$or_frame_index[rec1$slope_index[b$shared$keep]]

  # Populate a fresh `share` on model 1:
  new_share <- function() {
    sh <- new.env(parent = emptyenv())
    m1 <- api$dr_baseline_treatment_fit(
      reader, ts1, b$supports, b$regressions[[ts1$regression]],
      b$shared, rec1, b$or_frame_index, b$treatment_recipe,
      cc_glm, fx$ctrl, NULL, b$complete_support, sh
    )
    list(sh = sh, m1 = m1)
  }
  s0 <- new_share()
  expect_identical(c(s0$sh$calls, s0$sh$reuses, s0$sh$prep_calls, s0$sh$prep_reuses), c(1L, 0L, 1L, 0L))
  expect_true(is.list(api$dr_baseline_shared_prepared(
    s0$sh, ts2, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N
  )))

  # Equal bindings on model 2 borrow `prep` and equal frozen `dr_baseline_treatment_fit`:
  m2_prep <- api$dr_baseline_treatment_fit(
    reader, ts2, b$supports, b$regressions[[ts2$regression]],
    b$shared, rec1, b$or_frame_index, b$treatment_recipe,
    cc_glm, fx$ctrl, NULL, b$complete_support, s0$sh
  )
  m2_froz <- cc_frozen$dr_baseline_treatment_fit(
    reader, ts2, b$supports, b$regressions[[ts2$regression]],
    b$shared, rec1, b$or_frame_index, b$treatment_recipe,
    cc_glm, fx$ctrl, NULL, b$complete_support, new.env(parent = emptyenv())
  )
  expect_identical(cc_strip_dr_clipping(m2_prep), m2_froz)
  expect_identical(c(s0$sh$calls, s0$sh$reuses, s0$sh$prep_calls, s0$sh$prep_reuses), c(1L, 1L, 1L, 1L))

  # Malformed inputs to `dr_baseline_shared_prepared` cleanly return NULL without throwing:
  expect_null(api$dr_baseline_shared_prepared(NULL, ts2, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  expect_null(api$dr_baseline_shared_prepared(1L, ts2, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  expect_null(api$dr_baseline_shared_prepared(new.env(parent = emptyenv()), ts2, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  expect_null(api$dr_baseline_shared_prepared(s0$sh, NULL, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  ts_bad_sids1 <- ts2; ts_bad_sids1$support_ids <- c(1L, NA_integer_)
  ts_bad_sids2 <- ts2; ts_bad_sids2$support_ids <- c(1L, 999L)
  ts_bad_sids3 <- ts2; ts_bad_sids3$support_ids <- 1L
  expect_null(api$dr_baseline_shared_prepared(s0$sh, ts_bad_sids1, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  expect_null(api$dr_baseline_shared_prepared(s0$sh, ts_bad_sids2, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  expect_null(api$dr_baseline_shared_prepared(s0$sh, ts_bad_sids3, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  sup_bad_v <- b$supports; sup_bad_v[[ts2$support_ids[1L]]]$v <- NA_real_
  expect_null(api$dr_baseline_shared_prepared(s0$sh, ts2, sup_bad_v, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))

  # Unequal logit-relevant bindings decline `prep` and refuse with `cdid_dr_baseline_identity_error`:
  run_with_mut <- function(ts_m = ts2, sup_m = b$supports, rdr_m = reader, trt_m = b$treatment_recipe,
                           glm_m = cc_glm, ref_m = NULL) {
    s <- new_share()
    api$dr_baseline_treatment_fit(
      rdr_m, ts_m, sup_m, b$regressions[[ts2$regression]],
      b$shared, rec1, b$or_frame_index, trt_m,
      glm_m, fx$ctrl, ref_m, b$complete_support, s$sh
    )
  }
  ts_coh <- ts2; ts_coh$cohorts <- rev(ts2$cohorts)
  expect_error(run_with_mut(ts_m = ts_coh), class = "cdid_dr_baseline_identity_error")
  ts_n <- ts2; ts_n$n_observers <- ts2$n_observers - 1L
  expect_error(run_with_mut(ts_m = ts_n), class = "cdid_dr_baseline_identity_error")
  sup_v <- b$supports; sup_v[[ts2$support_ids[1L]]]$v <- 1.25
  expect_error(run_with_mut(sup_m = sup_v), class = "cdid_dr_baseline_identity_error")
  glm_alt <- cc_glm; glm_alt$maxit <- 25L
  expect_error(run_with_mut(glm_m = glm_alt), class = "cdid_dr_baseline_identity_error")
  rdr_w <- reader; rdr_w$state$w[1L] <- 2
  expect_error(run_with_mut(rdr_m = rdr_w), class = "cdid_dr_baseline_identity_error")
  trt_alt <- b$treatment_recipe; trt_alt$index <- c(1L, 3L)
  expect_error(run_with_mut(trt_m = trt_alt), class = "cdid_dr_baseline_identity_error")

  # Rebinding `frame_owner` inside `reader$frame_ref` with altered Z_a declines `prep` and refuses:
  fo_alt <- fx$fo; fo_alt$Z_a[1L, 2L] <- fo_alt$Z_a[1L, 2L] + 5
  old_fo <- get("frame_owner", envir = reader$frame_ref, inherits = FALSE)
  unlockBinding("frame_owner", reader$frame_ref)
  assign("frame_owner", fo_alt, envir = reader$frame_ref)
  lockBinding("frame_owner", reader$frame_ref)
  on.exit({
    unlockBinding("frame_owner", reader$frame_ref)
    assign("frame_owner", old_fo, envir = reader$frame_ref)
    lockBinding("frame_owner", reader$frame_ref)
  }, add = TRUE)
  expect_null(api$dr_baseline_shared_prepared(s0$sh, ts2, b$supports, reader, b$treatment_recipe, kept_or, b$shared$p, cc_glm, st$N))
  expect_error(
    api$dr_baseline_treatment_fit(
      reader, ts2, b$supports, b$regressions[[ts2$regression]],
      b$shared, rec1, b$or_frame_index, b$treatment_recipe,
      cc_glm, fx$ctrl, NULL, b$complete_support, s0$sh
    ),
    class = "cdid_dr_baseline_identity_error"
  )
  unlockBinding("frame_owner", reader$frame_ref)
  assign("frame_owner", old_fo, envir = reader$frame_ref)
  lockBinding("frame_owner", reader$frame_ref)

  # Changing OR-only bindings (`or_frame_index` / `kept_or_cols`) declines `prep`,
  # rebuilds `Z_or` from the new column while reusing the shared treatment logit,
  # and matches `cc_frozen$dr_baseline_treatment_fit` bitwise:
  s_or <- new_share()
  or_idx_alt <- c(1L, 3L)
  m2_or_alt <- api$dr_baseline_treatment_fit(
    reader, ts2, b$supports, b$regressions[[ts2$regression]],
    b$shared, rec1, or_idx_alt, b$treatment_recipe,
    cc_glm, fx$ctrl, NULL, b$complete_support, s_or$sh
  )
  m2_or_froz <- cc_frozen$dr_baseline_treatment_fit(
    reader, ts2, b$supports, b$regressions[[ts2$regression]],
    b$shared, rec1, or_idx_alt, b$treatment_recipe,
    cc_glm, fx$ctrl, NULL, b$complete_support, new.env(parent = emptyenv())
  )
  expect_identical(c(s_or$sh$calls, s_or$sh$reuses, s_or$sh$prep_calls, s_or$sh$prep_reuses), c(1L, 1L, 1L, 0L))
  expect_identical(cc_strip_dr_clipping(m2_or_alt), m2_or_froz)

  # Non-NULL `reference` with `complete` refuses with `cdid_dr_baseline_integrity_error`:
  expect_error(run_with_mut(ref_m = list(dated = TRUE)), class = "cdid_dr_baseline_integrity_error")
})
