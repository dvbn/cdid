# Private proposal for maintainer review: opt-in post-estimation output of the
# fitted nuisance models (initial vertical slice: observation-selection logits).
#
#   att_gt(..., control = list(keep_nuisance = TRUE))  -> fit$chained$nuisance_models
#   nuisance_models(fit, component = NULL, block = NULL) -> index / tables / statuses
#
# Default calls (keep_nuisance absent or FALSE) execute only the entry check in
# chained_nuisance_control(): no fit, stream, copy or retained object is added.
#
# Initial slice. selection_est_method = "ipw" with the built-in est_method =
# "ipw", "reg" or "dr" on the static (non-dated, fixed-reference) design, where
# the shared sealed `cdid_selection_models` backend of fit_selection_models()
# supplies the coefficients the estimator uses (IPW: fit_cp04_selection_models,
# OR: fit_or_selection_models, DR: fit_dr_corollary1). The R and native block
# producers fill the same block records, so the capture reads that one boundary.
# Next increment: the observation logits of the augmented observation correction
# (selection_est_method = "dr") on its ordinary static route, read from the same
# shared backend retained as cdid_augmented_backend$selection_backend/$up_map
# (chained_nuisance_augmented_scope()). With est_method = "dr" on that route the
# imputation, treatment-propensity and counterfactual OR fits are reported too
# (section 8); other augmented treatment methods keep typed statuses.
# Legacy increment: the built-in default selection_model = "pair_logit",
# selection_pooling = "separate" with a nonempty design runs the legacy
# per-(pair, cohort) consumer (selection_use_legacy_consumer()) for est_method =
# "ipw" (cdid_cp04_corollary1) and "reg" (cdid_or_corollary1). Both hold the
# same records from cp04_selection_fit()/cdid_fast_sel_cp04_fit_batch():
# selection_specs[[j]] = (pair, cohort), selections[[j]] = (coefficients over
# every reader design column, H_chol = chol(H/N), diagnostics, pi). They are
# read as one single-cell block per fit (chained_nuisance_legacy_view()) and
# mapped to the comparisons through treatment_specs$selection_ids (IPW) or the
# target support_id plus control-pool support_ids (OR).
# Every other route or component returns a typed status; nothing is refitted.
#
# Coefficient SEs reuse the estimator's own per-unit selection IF
# (selection_prediction_operation()$upstream -> evaluate_selection_block_if):
# xi_i = (H / N)^{-1} u_i with u_i the unit's weighted logit score summed over
# its cells of the block. Only the diagonal is accumulated, block by block:
#   Var(theta_j) = N^{-2} sum_g (S_gj - Sbar_j)^2,  S_gj = sum_{i in g} xi_ij,
# g = original unit (default) or the requested cluster; centered and without a
# small-sample factor, as the package's analytical ATT convention.

chained_nuisance_schema <- "cdid.nuisance_capture.v0"

chained_nuisance_selection_classes <- c("cdid_cp04_selection_backend", "cdid_or_selection_backend",
                                        "cdid_dr_selection_backend")

chained_nuisance_components <- c("selection", "outcome", "treatment", "imputation")

# ------------------------------------------------------------------------------
# 1. Public control entry (consumed before the strict normalizer)
# ------------------------------------------------------------------------------
# The public normalizer refuses unknown entries on purpose. `keep_nuisance` is
# removed here and every other entry reaches chained_public_resolve_control()
# unchanged, so all existing controls keep their meaning and refusals.
chained_nuisance_control <- function(control) {
  if (is.null(control) || !is.list(control) || is.object(control) || is.data.frame(control))
    return(list(keep = FALSE, control = control, supplied = FALSE))
  nms <- names(control)
  if (is.null(nms)) return(list(keep = FALSE, control = control, supplied = FALSE))
  hit <- which(!is.na(nms) & nms == "keep_nuisance")
  if (!length(hit)) return(list(keep = FALSE, control = control, supplied = FALSE))
  chained_public_need(length(hit) == 1L, "`control` entries must be uniquely named.",
    "cdid_chained_public_control_error")
  val <- control[[hit]]
  chained_public_need(is.logical(val) && !is.object(val) && is.null(dim(val)) && length(val) == 1L && !is.na(val),
    paste0("`control$keep_nuisance` must be TRUE or FALSE (default FALSE). TRUE keeps a compact capture of the ",
           "fitted models for nuisance_models(). Estimation and inference are unchanged."),
    "cdid_chained_public_control_error")
  list(keep = isTRUE(val), control = control[-hit], supplied = TRUE)
}

# ------------------------------------------------------------------------------
# 2. Pure streaming helpers (no owner access; shared with the pure-helper tests)
# ------------------------------------------------------------------------------
# Diagonal moment accumulator over IF row chunks. Unit level keeps two width
# vectors; cluster level keeps one n_groups x width matrix of cluster sums.
chained_nuisance_moments_new <- function(width, n_groups = NULL) {
  list(width = as.integer(width), s1 = numeric(width), s2 = numeric(width),
       gsum = if (is.null(n_groups)) NULL else matrix(0, n_groups, width), n_rows = 0)
}

chained_nuisance_moments_add <- function(acc, xi, group = NULL) {
  if (!is.matrix(xi) || ncol(xi) != acc$width || anyNA(xi) || any(!is.finite(xi)))
    cp04_abort("Selection IF chunk has the wrong width or nonfinite entries.", "integrity")
  if (is.null(acc$gsum)) {
    acc$s1 <- acc$s1 + colSums(xi)
    acc$s2 <- acc$s2 + colSums(xi * xi)
  } else {
    if (length(group) != nrow(xi) || anyNA(group))
      cp04_abort("Cluster index does not cover the IF chunk rows.", "integrity")
    rs <- rowsum(xi, group, reorder = FALSE)
    ix <- as.integer(rownames(rs))
    acc$gsum[ix, ] <- acc$gsum[ix, , drop = FALSE] + rs
  }
  acc$n_rows <- acc$n_rows + nrow(xi)
  acc
}

# Var_j = (sum_g S_gj^2 - (sum_g S_gj)^2 / G) / N^2 with G = N units (unit level)
# or the number of distinct clusters on the frame (zero-score groups included).
chained_nuisance_moments_variance <- function(acc, N, n_groups = N) {
  if (is.null(acc$gsum)) {
    s1 <- acc$s1; s2 <- acc$s2; G <- N
  } else {
    s1 <- colSums(acc$gsum); s2 <- colSums(acc$gsum * acc$gsum); G <- n_groups
  }
  v <- (s2 - s1 * s1 / G) / (as.double(N) * as.double(N))
  # Centered sums of squares are nonnegative; only rounding can cross zero.
  pmax(v, 0)
}

# Streams the U rows of one cohort range in canonical (increasing) order, in
# chunks of at most B rows, exactly as selection_owner_block_provider() does.
# `st` needs cohort_start, cohort_end, cohort_order and U only.
chained_nuisance_cohort_rows <- function(st, k, B, consume) {
  s0 <- st$cohort_start[k]; e0 <- st$cohort_end[k]
  if (!length(s0) || !is.finite(s0) || !is.finite(e0) || e0 < s0) return(invisible(NULL))
  for (start in seq.int(s0, e0, by = B)) {
    rows <- st$cohort_order[seq.int(start, min(e0, start + B - 1))]
    rows <- rows[st$U[rows]]
    if (length(rows)) consume(rows)
  }
  invisible(NULL)
}

# Block geometry from the fitted selection records (metadata only).
chained_nuisance_block_geometry <- function(sel, up_map) {
  M <- length(sel$blocks)
  if (length(up_map$theta_coords) != M)
    cp04_abort("Selection upstream map does not index the selection blocks.", "integrity")
  lapply(seq_len(M), function(m) {
    cids <- as.integer(sel$blocks[[m]]$cell_ids)
    ck <- vapply(sel$cells[cids], function(cl) as.integer(cl$cohort), integer(1))
    uids <- as.integer(up_map$theta_coords[[m]])
    list(cids = cids, ck = ck, coh = sort(unique(ck)), p = length(uids), uids = uids)
  })
}

# One bounded pass per block over the U rows of its cohorts and, for date
# models, one pass per pair-cohort for the derived product. `fns` supplies the
# estimator's own IF columns (influence(rows, m)), the fitted cell probability
# and response (cell(rows, cell_id) -> list(p, y)) and the pair probability
# (pair(rows, q)). `run = FALSE` returns the same structure filled with NA.
chained_nuisance_selection_passes <- function(sel, blk, st, B, N, fns, group = NULL, n_groups = N,
                                              do_cov = TRUE, run = TRUE) {
  C <- length(sel$cells); n_pc <- length(sel$pair_cohorts)
  is_date <- identical(sel$selection_model, "independent_period_logit")
  var_list <- lapply(blk, function(x) rep(NA_real_, x$p))
  block_groups <- rep(NA_integer_, length(blk))
  cell_stats <- matrix(NA_real_, C, 5L, dimnames = list(NULL, c("n_rows", "n_observed", "p_min", "p_max", "n_floor")))
  n_U <- rep(NA_real_, length(st$cohort_start))
  pc_stats <- matrix(NA_real_, n_pc, 2L, dimnames = list(NULL, c("p_min", "n_floor")))
  if (isTRUE(run)) for (m in seq_along(blk)) {
    info <- blk[[m]]
    acc <- chained_nuisance_moments_new(info$p, if (is.null(group)) NULL else n_groups)
    groups_seen <- if (is.null(group)) NULL else logical(n_groups)
    cs <- matrix(c(0, 0, Inf, -Inf, 0), length(info$cids), 5L, byrow = TRUE)
    for (k in info$coh) {
      nk <- 0
      chained_nuisance_cohort_rows(st, k, B, function(rows) {
        nk <<- nk + length(rows)
        if (do_cov && info$p > 0L) {
          xi <- fns$influence(rows, m)
          if (!is.matrix(xi)) xi <- matrix(xi, nrow = length(rows))
          g <- if (is.null(group)) NULL else group[rows]
          if (!is.null(g)) groups_seen[g] <<- TRUE
          acc <<- chained_nuisance_moments_add(acc, xi, g)
        }
        for (j in which(info$ck == k)) {
          py <- fns$cell(rows, info$cids[j])
          cs[j, ] <<- c(cs[j, 1L] + length(rows), cs[j, 2L] + sum(py$y), min(cs[j, 3L], py$p),
                        max(cs[j, 4L], py$p), cs[j, 5L] + sum(py$p <= ipw_clip_eps))
        }
      })
      n_U[k] <- nk
    }
    cs[!is.finite(cs[, 3L]), 3L] <- NA_real_
    cs[!is.finite(cs[, 4L]), 4L] <- NA_real_
    cell_stats[info$cids, ] <- cs
    if (do_cov && info$p > 0L) {
      if (!is.null(groups_seen)) block_groups[m] <- sum(groups_seen)
      if (is.null(groups_seen) || block_groups[m] >= 2L)
        var_list[[m]] <- chained_nuisance_moments_variance(acc, N, n_groups)
    }
  }
  for (q in seq_len(n_pc)) {
    pc <- sel$pair_cohorts[[q]]
    if (!is_date) {
      pc_stats[q, ] <- cell_stats[pc$pair_cell_id, c("p_min", "n_floor")]
      next
    }
    if (!isTRUE(run)) next
    # Date models: the pair probability is the derived product of two date
    # logits (no separate fit); its weight floor is checked on the product.
    pmin_q <- Inf; nf <- 0
    chained_nuisance_cohort_rows(st, pc$cohort, B, function(rows) {
      pr <- fns$pair(rows, q)
      pmin_q <<- min(pmin_q, pr); nf <<- nf + sum(pr <= ipw_clip_eps)
    })
    pc_stats[q, ] <- c(if (is.finite(pmin_q)) pmin_q else NA_real_, nf)
  }
  list(var_list = var_list, cell_stats = cell_stats, n_U = n_U, pc_stats = pc_stats,
       block_groups = block_groups)
}

# Coefficient table in the conventional summary.glm columns.
chained_nuisance_coef_frame <- function(block, term, role, cell, estimate, variance) {
  se <- sqrt(variance)
  z <- rep(NA_real_, length(se))
  usable <- is.finite(se) & se > 0
  z[usable] <- estimate[usable] / se[usable]
  data.frame(block = rep(block, length(term)), term = term, role = role, cell = as.integer(cell),
    estimate = as.double(estimate), std_error = as.double(se), z_value = as.double(z),
    p_value = as.double(2 * stats::pnorm(-abs(z))), stringsAsFactors = FALSE)
}

chained_nuisance_covariance_record <- function(available, code, reason, clustered, cluster_var, n_groups, N,
                                               source = NULL) {
  if (is.null(source))
    source <- paste("Per-unit IF columns of the estimator's own selection IF (selection_prediction_operation()$upstream ->",
      "evaluate_selection_block_if); no generic inverse Hessian and no second logit.")
  list(
    status = if (available) "available" else if (identical(code, "fixed_probability")) "not_applicable" else "unavailable",
    code = if (available) "ok" else code,
    reason = if (available) NULL else reason,
    type = if (clustered) "cluster" else "unit", cluster_var = cluster_var, n_groups = n_groups, N = N,
    estimator = paste0("Centered sandwich: each original unit's weighted logit score is summed over its cells ",
      "and multiplied by (H/N)^{-1}", if (clustered) ", then summed within its cluster" else "",
      ". The squared centered contributions are summed and divided by N^2, with no small-sample factor. ",
      "H is the unnormalized weighted Hessian."),
    source = source,
    diagonal_only = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Labels and equations (pure; derived from the fitted design metadata)
# ------------------------------------------------------------------------------
chained_nuisance_g <- function(g) ifelse(is.infinite(g), 0, g)

chained_nuisance_num1 <- function(x) if (length(x) == 1L && is.numeric(x)) as.double(x) else NA_real_

chained_nuisance_num <- function(x) {
  x <- as.double(x)
  out <- as.character(x)
  out[is.na(x)] <- "NA"
  out
}

chained_nuisance_cell_key <- function(cell, with_g = TRUE) {
  g <- chained_nuisance_num(chained_nuisance_g(cell$cohort_value))
  idx <- if (identical(cell$type, "date")) paste0("t=", chained_nuisance_num(cell$t)) else
    paste0("(", chained_nuisance_num(cell$pair$s), ",", chained_nuisance_num(cell$pair$t), ")")
  if (with_g) paste0("g=", g, ",", idx) else idx
}

chained_nuisance_block_label <- function(pooling, cells, cell_ids) {
  if (identical(pooling, "separate")) {
    cell <- cells[[cell_ids[1L]]]
    return(paste0(if (identical(cell$type, "date")) "date[" else "pair[", chained_nuisance_cell_key(cell), "]"))
  }
  gs <- unique(vapply(cells[cell_ids], function(cl) chained_nuisance_g(cl$cohort_value), numeric(1)))
  if (identical(pooling, "within_cohort")) paste0("within_cohort[g=", chained_nuisance_num(gs[1L]), "]") else
    paste0("across_cohorts[g in {", paste(chained_nuisance_num(sort(gs)), collapse = ","), "}]")
}

chained_nuisance_slope_text <- function(slope_names, max_terms = 6L) {
  d_s <- length(slope_names)
  if (!d_s) return(character(0))
  shown <- seq_len(min(d_s, max_terms))
  txt <- paste0("b", shown, "*", slope_names[shown])
  if (d_s > max_terms) txt <- c(txt, paste0("... (", d_s - max_terms, " more slope terms)"))
  txt
}

# One template per distinct specification (model x pooling x design branch).
chained_nuisance_selection_template <- function(selection_model, pooling, branch) {
  is_date <- identical(selection_model, "independent_period_logit")
  cell_ix <- if (is_date) "[g,t]" else "[g,(s,t)]"
  slope_ix <- switch(pooling, separate = cell_ix, within_cohort = "[g]", across_cohorts = "")
  lhs <- if (is_date) "logit Pr(R_it = 1 | X_i, G_i = g)" else "logit Pr(S_i(s,t) = 1 | X_i, G_i = g, U_i = 1)"
  rhs <- switch(branch,
    intercept_slope = paste0("a", cell_ix, " + X_i'b", slope_ix),
    intercept_only = paste0("a", cell_ix),
    no_intercept = paste0("X_i'b", slope_ix),
    zero_col = "0  (fixed: Pr = 1/2, no fitted parameters)",
    paste0("a", cell_ix, " + X_i'b", slope_ix))
  paste(lhs, "=", rhs)
}

chained_nuisance_block_equation <- function(selection_model, cells, cell_ids, has_intercept, slope_names) {
  is_date <- identical(selection_model, "independent_period_logit")
  cl <- cells[cell_ids]
  gs <- unique(vapply(cl, function(x) chained_nuisance_g(x$cohort_value), numeric(1)))
  one_g <- length(gs) == 1L
  single <- length(cell_ids) == 1L
  g_txt <- if (one_g) chained_nuisance_num(gs) else "g"
  lhs <- if (is_date) {
    t_txt <- if (single) chained_nuisance_num(cl[[1L]]$t) else "t"
    paste0("logit Pr(R_i", if (single) paste0("[t=", t_txt, "]") else "t", " = 1 | X_i, G_i = ", g_txt, ")")
  } else {
    st_txt <- if (single) paste0(chained_nuisance_num(cl[[1L]]$pair$s), ",", chained_nuisance_num(cl[[1L]]$pair$t)) else "s,t"
    paste0("logit Pr(S_i(", st_txt, ") = 1 | X_i, G_i = ", g_txt, ", U_i = 1)")
  }
  a_txt <- if (!isTRUE(has_intercept)) character(0) else if (single) "a" else if (one_g) {
    if (is_date) "a[t]" else "a[(s,t)]"
  } else if (is_date) "a[g,t]" else "a[g,(s,t)]"
  terms <- c(a_txt, chained_nuisance_slope_text(slope_names))
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "0  (fixed: Pr = 1/2, no fitted parameters)"
  paste(lhs, "=", rhs)
}

chained_nuisance_universe_text <- function(universe) {
  switch(as.character(universe),
    F = "U = F: every unit of the representative frame",
    B = "U = B: frame units with at least one observed outcome",
    A = "U = A: frame units with at least two observed outcomes",
    paste0("U = ", universe))
}

# Compact descriptive tables and the conventional coefficient table from the
# fitted block records, the pass statistics and request metadata (pure).
# `meta`: cohorts, mass_U (per N), N, design_columns, universe, formula,
# population, weights (list with text), rho_column, glm_control, route,
# covariance (record) and accounting.
chained_nuisance_selection_tables <- function(sel, blk, stats, meta) {
  M <- length(blk); C <- length(sel$cells); n_pc <- length(sel$pair_cohorts)
  ds <- sel$design_spec; zn <- meta$design_columns; N <- meta$N
  is_date <- identical(sel$selection_model, "independent_period_logit")
  cell_stats <- stats$cell_stats; n_U <- stats$n_U; pc_stats <- stats$pc_stats
  int_name <- if (isTRUE(ds$has_intercept) && length(ds$int_col) && !is.null(zn) && nzchar(zn[ds$int_col]))
    zn[ds$int_col] else "(Intercept)"
  pub_cohort <- function(k) chained_nuisance_g(meta$cohorts[k])
  links <- if (is.null(sel$links)) list() else sel$links

  block_ids <- paste0("S", seq_len(M))
  block_labels <- vapply(seq_len(M), function(m)
    chained_nuisance_block_label(sel$selection_pooling, sel$cells, blk[[m]]$cids), character(1))
  equations <- vapply(seq_len(M), function(m) {
    bl <- sel$blocks[[m]]
    if (identical(bl$branch, "zero_col"))
      return(chained_nuisance_block_equation(sel$selection_model, sel$cells, blk[[m]]$cids, FALSE, character(0)))
    chained_nuisance_block_equation(sel$selection_model, sel$cells, blk[[m]]$cids, bl$has_intercept,
      if (bl$d_s > 0L) ds$slope_names else character(0))
  }, character(1))

  # Comparisons (links) and the pair-cohorts they read.
  link_pc <- lapply(links, function(l) as.integer(l$pair_cohort_ids))
  pc_links <- tabulate(unlist(link_pc, use.names = FALSE), n_pc)
  pc_blocks <- lapply(sel$pair_cohorts, function(pc) as.integer(pc$block_ids))
  block_pcs <- split(rep(seq_len(n_pc), lengths(pc_blocks)),
                     factor(unlist(pc_blocks, use.names = FALSE), levels = seq_len(M)))
  link_blocks <- lapply(link_pc, function(ids) unique(unlist(pc_blocks[ids], use.names = FALSE)))
  block_links <- tabulate(unlist(link_blocks, use.names = FALSE), M)

  coefs <- do.call(rbind, c(list(chained_nuisance_coef_frame(character(0), character(0), character(0), integer(0),
    numeric(0), numeric(0))), lapply(seq_len(M), function(m) {
    bl <- sel$blocks[[m]]; info <- blk[[m]]
    if (identical(bl$branch, "zero_col") || info$p == 0L) return(NULL)
    est <- c(as.double(bl$alpha), as.double(bl$beta))
    n_al <- length(bl$alpha)
    if (length(est) != info$p)
      cp04_abort("Selection block coefficients disagree with their upstream coordinates.", "integrity")
    term <- c(vapply(info$cids[seq_len(n_al)], function(cid)
      paste0(int_name, "[", chained_nuisance_cell_key(sel$cells[[cid]]), "]"), character(1)),
      if (length(bl$beta)) ds$slope_names else character(0))
    role <- c(rep("local_intercept", n_al), rep(if (length(info$cids) > 1L) "shared_slope" else "slope", length(bl$beta)))
    cell <- c(info$cids[seq_len(n_al)], rep(NA_integer_, length(bl$beta)))
    chained_nuisance_coef_frame(block_ids[m], term, role, cell, est, stats$var_list[[m]])
  })))
  rownames(coefs) <- NULL

  blocks_df <- data.frame(
    block = block_ids, label = block_labels,
    cohorts = vapply(blk, function(x) paste(chained_nuisance_num(pub_cohort(x$coh)), collapse = ","), character(1)),
    n_cells = vapply(blk, function(x) length(x$cids), integer(1)),
    n_parameters = vapply(blk, `[[`, integer(1), "p"),
    n_intercepts = vapply(sel$blocks, function(bl) if (identical(bl$branch, "zero_col")) 0L else length(bl$alpha), integer(1)),
    n_slopes = vapply(sel$blocks, function(bl) if (identical(bl$branch, "zero_col")) 0L else length(bl$beta), integer(1)),
    n_fitting_rows = vapply(blk, function(x) sum(cell_stats[x$cids, "n_rows"]), numeric(1)),
    n_units = vapply(blk, function(x) sum(n_U[x$coh]), numeric(1)),
    n_fitting_clusters = stats$block_groups,
    unit_weight_sum = vapply(blk, function(x) sum(meta$mass_U[x$coh]) * N, numeric(1)),
    branch = vapply(sel$blocks, `[[`, character(1), "branch"),
    converged = vapply(sel$blocks, function(bl) isTRUE(bl$diagnostics$converged), logical(1)),
    iterations = vapply(sel$blocks, function(bl) chained_nuisance_num1(bl$diagnostics$iter), numeric(1)),
    rank = vapply(sel$blocks, function(bl) chained_nuisance_num1(bl$diagnostics$rank), numeric(1)),
    n_pair_cohorts = vapply(block_pcs, length, integer(1)),
    n_comparisons = block_links,
    equation = equations,
    stringsAsFactors = FALSE)

  cells_df <- data.frame(
    cell = seq_len(C),
    block = block_ids[vapply(sel$cells, function(cl) as.integer(cl$block_id), integer(1))],
    type = vapply(sel$cells, `[[`, character(1), "type"),
    cohort = vapply(sel$cells, function(cl) chained_nuisance_g(cl$cohort_value), numeric(1)),
    s = vapply(sel$cells, function(cl) if (identical(cl$type, "pair")) as.double(cl$pair$s) else NA_real_, numeric(1)),
    t = vapply(sel$cells, function(cl) if (identical(cl$type, "pair")) as.double(cl$pair$t) else as.double(cl$t), numeric(1)),
    key = vapply(sel$cells, chained_nuisance_cell_key, character(1)),
    n_rows = cell_stats[, "n_rows"], n_observed = cell_stats[, "n_observed"],
    observed_share_weighted = vapply(sel$cells, function(cl)
      chained_nuisance_num1(sel$blocks[[cl$block_id]]$pi[cl$local_b]), numeric(1)),
    p_min = cell_stats[, "p_min"], p_max = cell_stats[, "p_max"], n_floor = cell_stats[, "n_floor"],
    stringsAsFactors = FALSE)

  pc_df <- data.frame(
    pair_cohort = seq_len(n_pc),
    cohort = vapply(sel$pair_cohorts, function(pc) chained_nuisance_g(pc$cohort_value), numeric(1)),
    s = vapply(sel$pair_cohorts, function(pc) as.double(pc$pair$s), numeric(1)),
    t = vapply(sel$pair_cohorts, function(pc) as.double(pc$pair$t), numeric(1)),
    probability = rep(if (is_date) "derived product Pr(R_is=1|X,g) * Pr(R_it=1|X,g) of two date logits (no separate fit)" else
      "pair logit Pr(S_i(s,t)=1|X,g)", n_pc),
    cells = vapply(sel$pair_cohorts, function(pc) paste(as.integer(pc$cell_ids), collapse = ","), character(1)),
    blocks = vapply(pc_blocks, function(x) paste(block_ids[x], collapse = ","), character(1)),
    p_min = pc_stats[, "p_min"], n_floor = pc_stats[, "n_floor"],
    n_comparisons = as.integer(pc_links),
    stringsAsFactors = FALSE)

  link_df <- data.frame(
    comparison = vapply(links, function(l) as.character(l$link_id), character(1)),
    g = vapply(links, function(l) pub_cohort(l$cohorts[1L]), numeric(1)),
    s = vapply(links, function(l) as.double(l$pair$s), numeric(1)),
    t = vapply(links, function(l) as.double(l$pair$t), numeric(1)),
    controls = vapply(links, function(l) paste(chained_nuisance_num(pub_cohort(l$cohorts[-1L])), collapse = ","), character(1)),
    pair_cohorts = vapply(link_pc, function(ids) paste(ids, collapse = ","), character(1)),
    stringsAsFactors = FALSE)

  out <- list(
    status = "available",
    role = if (is.null(meta$role)) "observation-selection propensity (conditional_mtar)" else meta$role,
    route = meta$route,
    selection_model = sel$selection_model, selection_pooling = sel$selection_pooling, link = "logit",
    response = if (is_date) "R_it = 1{Y_it observed}" else "S_i(s,t) = 1{Y_is and Y_it both observed}",
    template = chained_nuisance_selection_template(sel$selection_model, sel$selection_pooling, ds$branch),
    derived = if (is_date) "p_g(s,t | X) = Pr(R_is = 1 | X, g) * Pr(R_it = 1 | X, g): a transformation of the date logits, not a second fit" else NULL,
    formula = meta$formula,
    design_columns = zn, intercept_term = int_name, slope_terms = ds$slope_names, design_branch = ds$branch,
    universe = meta$universe, universe_text = chained_nuisance_universe_text(meta$universe),
    population = meta$population, weights = meta$weights, transport = meta$transport,
    n_frame = N, glm_control = meta$glm_control, probability_floor = ipw_clip_eps,
    blocks = blocks_df, cells = cells_df, coefficients = coefs, pair_cohorts = pc_df, comparisons = link_df,
    covariance = meta$covariance, accounting = meta$accounting)
  # How a non-default route uses the fitted probability (augmented correction).
  if (!is.null(meta$use)) out$use <- meta$use
  out
}

# ------------------------------------------------------------------------------
# 4. Scope decision and typed statuses
# ------------------------------------------------------------------------------
chained_nuisance_status <- function(status, code, reason) list(status = status, code = code, reason = reason)

chained_nuisance_selection_scope <- function(b, prep, gate) {
  if (identical(gate$selection, "none"))
    return(chained_nuisance_status("not_applicable", "no_selection_model",
      'selection = "none": no observation-selection model is fitted.'))
  if (isTRUE(gate$is_custom_est))
    return(chained_nuisance_status("unavailable", "custom_callback_route",
      "Custom IPW callbacks use their own consumer of the selection fit; outside the initial slice."))
  if (identical(gate$selection_est_method, "imputation"))
    return(chained_nuisance_status("not_applicable", "pure_imputation",
      'selection_est_method = "imputation" fits no observation-selection logit.'))
  if (identical(gate$selection_est_method, "dr")) return(chained_nuisance_augmented_scope(b, prep))
  if (!is.null(prep$reference_source) || (is.list(b) && !is.null(b$reference_ref)))
    return(chained_nuisance_status("unavailable", "reference_recipe_route",
      "Dated covariates or recipe-indexed weights use reference selection models outside this reporting scope."))
  if (is.null(b))
    return(chained_nuisance_status("unavailable", "no_full_roster_backend",
      "The fit retained no full-roster selection backend."))
  if (!isTRUE(b$uses_selection_models)) return(chained_nuisance_legacy_scope(b))
  if (is.null(b$selection_backend))
    return(chained_nuisance_status("unavailable", "unrecognized_backend",
      "The backend declares the shared selection-models consumer but holds no selection backend."))
  if (!inherits(b, chained_nuisance_selection_classes) || inherits(b, "cdid_cp04_custom_backend") ||
      !is.list(b$up_map) || length(b$up_map$theta_coords) != length(b$selection_backend$blocks))
    return(chained_nuisance_status("unavailable", "unrecognized_backend",
      paste0("Backend class ", class(b)[1L], " does not expose the shared selection backend with its upstream map; ",
             "outside the initial slice.")))
  NULL
}

# Augmented observation correction (selection_est_method = "dr"). Eligible only
# on the ordinary static route, the same predicate as
# cp04_retained_if_phase_static_augmented(): no recipe/reference layout and no
# private supplied pair differences. There fit_augmented_backend() calls
# fit_selection_models() (or its native twin cdid_ia_fast_fit_selection_models)
# on the backend's own static reader with U = F and rho = 1, so the logits have
# the earlier capture's response (panel R_it / S_i(s,t)), weights st$w, rows
# (frame units of each cohort), score and normalization N. Anything else is a
# typed status; nothing is refitted.
chained_nuisance_augmented_scope <- function(b, prep) {
  if (!is.null(prep$reference_source) ||
      (is.list(b) && (!is.null(b$reference_ref) || !is.null(b$layout))))
    return(chained_nuisance_status("unavailable", "reference_recipe_route", paste0(
      "Dated covariates or recipe-indexed weights fit the augmented observation logits on the recipe-cohort layout ",
      "(weights W_g, layout rows); their score, weights and rows differ from the static capture.")))
  if (is.null(b))
    return(chained_nuisance_status("unavailable", "no_full_roster_backend",
      "The fit retained no full-roster selection backend."))
  if (!inherits(b, "cdid_augmented_backend"))
    return(chained_nuisance_status("unavailable", "unrecognized_backend",
      paste0("Backend class ", class(b)[1L], " is not the augmented observation-correction backend.")))
  if (!is.null(b$pair_outcomes))
    return(chained_nuisance_status("unavailable", "augmented_route", paste0(
      "This augmented fit carries privately supplied pair differences instead of reading them from the panel; ",
      "only the ordinary static augmented route is captured.")))
  sb <- b$selection_backend
  if (!isTRUE(b$uses_selection_models) || !identical(b$universe, "F") ||
      !inherits(sb, "cdid_selection_models") || !is.null(sb$reference_ref) ||
      !is.environment(b$reader_ref) || !identical(sb$reader_ref, b$reader_ref) ||
      !is.list(b$up_map) || length(b$up_map$theta_coords) != length(sb$blocks))
    return(chained_nuisance_status("unavailable", "unrecognized_backend", paste0(
      "The augmented backend does not hold a static shared selection backend on its own reader with a matching ",
      "upstream map.")))
  st <- if (exists("cp04_reader", b$reader_ref, inherits = FALSE))
    get("cp04_reader", b$reader_ref, inherits = FALSE)$state else NULL
  if (!is.list(st) || !identical(st$universe, "F") || !identical(st$id, prep$frame_ids) ||
      length(st$U) != length(prep$frame_ids) || !all(st$U) || !isTRUE(all(st$rho == 1)))
    return(chained_nuisance_status("unavailable", "unrecognized_backend", paste0(
      "The observation logits' reader rows are not the frame's original units on U = F with rho = 1.")))
  NULL
}

chained_nuisance_legacy_classes <- c("cdid_cp04_corollary1", "cdid_or_corollary1")

# Legacy separate pair-logit consumer (selection_use_legacy_consumer(): pair
# logit, separate pooling, nonempty design) on its static route. Eligible: the
# IPW (cdid_cp04_corollary1) and OR (cdid_or_corollary1) backends without a
# reference/recipe layout or private pair differences, whose reader holds the
# design and weights and whose selections are one finite coefficient vector over
# the reader design per (pair, cohort) spec. The comparison map is checked in
# chained_nuisance_legacy_links(). Metadata only; no row is read here.
chained_nuisance_legacy_scope <- function(b) {
  if (!inherits(b, chained_nuisance_legacy_classes) || inherits(b, "cdid_cp04_custom_backend"))
    return(chained_nuisance_status("unavailable", "legacy_separate_pair_consumer", paste0(
      "Backend class ", class(b)[1L], " does not hold the shared selection-models backend and is not one of the ",
      "captured legacy separate pair-logit consumers (cdid_cp04_corollary1 for est_method = \"ipw\", ",
      "cdid_or_corollary1 for \"reg\").")))
  if (!is.null(b$reference_ref) || !is.null(b$layout) || !is.null(b$pair_outcomes))
    return(chained_nuisance_status("unavailable", "reference_recipe_route", paste0(
      "This legacy backend carries a reference/recipe layout or private pair differences; only the static legacy ",
      "route is captured.")))
  reader <- if (is.environment(b$reader_ref) && exists("cp04_reader", b$reader_ref, inherits = FALSE))
    get("cp04_reader", b$reader_ref, inherits = FALSE) else NULL
  st <- if (is.list(reader)) reader$state else NULL
  S <- length(b$selection_specs)
  d <- if (is.list(st) && is.matrix(st$Z)) ncol(st$Z) else 0L
  fit_ok <- function(f) is.list(f) && is.numeric(f$coefficients) && length(f$coefficients) == d &&
    all(is.finite(f$coefficients))
  if (!d || isTRUE(st$reference_carrier) || !is.list(b$selection_specs) || !is.list(b$selections) || !S ||
      length(b$selections) != S || !is.list(b$links) || !length(b$links) ||
      !all(vapply(b$selections, fit_ok, logical(1))))
    return(chained_nuisance_status("unavailable", "unrecognized_backend", paste0(
      "The legacy separate pair-logit backend does not hold one finite coefficient vector over its reader design ",
      "for every (pair, cohort) fit.")))
  NULL
}

chained_nuisance_augmented_use <- paste0(
  "q_c, this fit's probability of the cell's pair (date models: the derived product), enters the augmented ",
  "response Z_c = j_c + S_e / max(q_c, eps) * (Delta_e - j_c) on observed pairs (S_e = 1); Z_c = j_c otherwise")

# Augmented route: only the observation propensity is captured. The other
# roles are fitted by the augmented backend and keep explicit statuses.
chained_nuisance_augmented_other <- function(est_label) {
  step <- paste("Not captured: beyond the observation logits, the report covers these roles only for est_method = \"dr\"",
    "on the ordinary static augmented route.")
  list(
    outcome = if (est_label %in% c("reg", "dr"))
      chained_nuisance_status("not_captured", "outside_initial_slice", paste(
        "The counterfactual outcome regression fits the reconstructed response Z_c, which depends on the",
        "fitted imputation and observation propensity; its SEs would need their derivatives and joint covariance.", step)) else
      chained_nuisance_status("not_applicable", "no_outcome_regression",
        paste0('est_method = "', est_label, '" fits no outcome regression.')),
    treatment = if (identical(est_label, "custom"))
      chained_nuisance_status("unavailable", "custom_callback",
        "A custom IPW callback supplies its own treatment weights; cdid holds no fitted treatment model for it.") else
      if (est_label %in% c("ipw", "dr"))
        chained_nuisance_status("not_captured", "outside_initial_slice",
          paste("Roster treatment propensity logits of the augmented route.", step)) else
        chained_nuisance_status("not_applicable", "no_treatment_propensity", 'est_method = "reg" fits no treatment propensity.'),
    imputation = chained_nuisance_status("not_captured", "outside_initial_slice", paste(
      "Actual-change imputation regressions j_c (the predictor in Z_c), fitted on observed pairs.", step)))
}

chained_nuisance_other_components <- function(gate, est_label) {
  if (identical(gate$selection_est_method, "dr")) return(chained_nuisance_augmented_other(est_label))
  slice <- "Not captured: the initial slice captures observation-selection logits only."
  sem <- gate$selection_est_method
  list(
    outcome = if (est_label %in% c("reg", "dr"))
      chained_nuisance_status("not_captured", "outside_initial_slice",
        paste0('Outcome-change regressions of est_method = "', est_label, '". ', slice)) else
      chained_nuisance_status("not_applicable", "no_outcome_regression",
        paste0('est_method = "', est_label, '" fits no outcome regression.')),
    treatment = if (identical(est_label, "custom"))
      chained_nuisance_status("unavailable", "custom_callback",
        "A custom IPW callback supplies its own treatment weights; cdid holds no fitted treatment model for it.") else
      if (est_label %in% c("ipw", "dr"))
        chained_nuisance_status("not_captured", "outside_initial_slice", paste("Treatment propensity logits.", slice)) else
        chained_nuisance_status("not_applicable", "no_treatment_propensity", 'est_method = "reg" fits no treatment propensity.'),
    imputation = if (sem %in% c("imputation", "dr"))
      chained_nuisance_status("not_captured", "outside_initial_slice", paste("Actual-change imputation regressions.", slice)) else
      chained_nuisance_status("not_applicable", "no_imputation",
        paste0('selection_est_method = "', sem, '" fits no actual-change imputation.')))
}

chained_nuisance_capture_assemble <- function(selection, gate, est_label, roles = NULL) {
  other <- chained_nuisance_other_components(gate, est_label)
  if (!is.null(roles)) other[names(roles)] <- roles
  structure(list(
    schema = chained_nuisance_schema, status = "captured", entry = "att_gt",
    est_method = est_label, selection_est_method = gate$selection_est_method, selection = gate$selection,
    components = c(list(selection = selection), other),
    source = paste("captured from the actual estimation objects of this call (no refit, no replay of a saved call);",
                   "compact: plain lists, data frames and vectors only")),
    class = "cdid_nuisance_capture")
}

# ------------------------------------------------------------------------------
# 5. Capture from the actual fit (opt-in only; reads the sealed owners)
# ------------------------------------------------------------------------------
chained_nuisance_capture <- function(res, prep, gate, est_label, weightsname) {
  b <- res$nuisance_registry$cp04
  scope <- chained_nuisance_selection_scope(b, prep, gate)
  selection <- if (!is.null(scope)) scope else if (isTRUE(b$uses_selection_models))
    chained_nuisance_selection_component(res, b, prep, weightsname) else
    chained_nuisance_legacy_component(res, b, prep, weightsname)
  # Augmented DR+DR: the imputation, treatment and counterfactual OR fits, after
  # (and charged on top of) the retained selection tables.
  roles <- if (identical(gate$selection_est_method, "dr") && identical(est_label, "dr"))
    chained_nuisance_augmented_roles(res, b, prep, weightsname, scope, selection$accounting)
  cap <- chained_nuisance_capture_assemble(selection, gate, est_label, roles)
  cap$entries <- unname(cp04_payload(cap)["entries"])
  cap
}

chained_nuisance_selection_component <- function(res, b, prep, weightsname) {
  sb <- b$selection_backend
  validate_selection_models(sb)
  reader <- get("cp04_reader", sb$reader_ref, inherits = FALSE)
  is_date <- identical(sb$selection_model, "independent_period_logit")
  augmented <- inherits(b, "cdid_augmented_backend")
  up_map <- b$up_map
  # Opened only after admission of the bounded pass.
  open <- function(blk) {
    op <- selection_prediction_operation(sb, source_reader = reader)
    list(abort = op$abort, finish = op$finish, fns = list(
      influence = function(rows, m) op$upstream(up_map, rows, blk[[m]]$uids)$up,
      cell = function(rows, cid) {
        cell <- sb$cells[[cid]]
        list(p = op$cell(rows, cid),
             y = if (is_date) cp04_date_rows(reader, rows, cell$date_spec)$R else cp04_pair_rows(reader, rows, cell$pair)$S)
      },
      pair = function(rows, q) {
        pc <- sb$pair_cohorts[[q]]
        op$cell(rows, pc$s_cell_id) * op$cell(rows, pc$t_cell_id)
      }))
  }
  src <- list(sel = sb, up_map = up_map, reader = reader, ctrl = sb$control, glm_control = sb$glm_control,
    open = open,
    route = if (augmented) paste0(class(b)[1L], "$selection_backend (shared cdid_selection_models backend; ",
      "ordinary static augmented route)") else paste0(class(b)[1L], " (shared cdid_selection_models backend)"),
    role = if (augmented) 'observation propensity q_c of the augmented observation correction (selection_est_method = "dr")',
    use = if (augmented) chained_nuisance_augmented_use,
    weights_use = if (augmented)
      paste("The fitted q_c enters the augmented response Z_c through 1 / max(q_c, eps) and the repeat rate a_k enters",
        "the change total a_k H_N[w I_k Z_c h]; neither is in this logit fit, which uses no other estimated model, so",
        "the reported SEs are not conditional on any estimated weight.") else
      paste("The repeat rates a_k and the transport ratio rho enter the downstream multiplier a_k / (rho p), not this",
        "logit fit, so the reported SEs are not conditional on any estimated weight."))
  chained_nuisance_selection_run(res, src, prep, weightsname)
}

# Shared admission, bounded passes and tables for one selection source `src`:
# sel (block/cell/pair-cohort/link records in the shared-backend shape), up_map,
# reader, ctrl, glm_control, open(blk) -> list(fns, finish, abort), route, role,
# use, weights_use and optional if_source / cov_code / cov_reason.
chained_nuisance_selection_run <- function(res, src, prep, weightsname, view_entries = 0, view_work = 0) {
  sb <- src$sel; up_map <- src$up_map; reader <- src$reader
  st <- reader$state
  ds <- sb$design_spec
  ctrl <- src$ctrl
  N <- st$N
  M <- length(sb$blocks); C <- length(sb$cells); n_pc <- length(sb$pair_cohorts)
  is_date <- identical(sb$selection_model, "independent_period_logit")
  cp <- prep$cp04
  n_parameters <- sum(lengths(up_map$theta_coords))
  formula_text <- if (is.null(cp$formula)) NA_character_ else paste(deparse(cp$formula), collapse = " ")
  text_entries <- function(x) sum(ceiling(nchar(x, type = "bytes", allowNA = TRUE) / 2), na.rm = TRUE)
  design_columns <- colnames(st$Z)
  # Each block repeats the shared slope labels in its equation and table.
  # Every local intercept label can repeat the intercept column name.
  label_entries <- 4 * (M + 1) * text_entries(ds$slope_names) +
    2 * (C + 1) * text_entries(design_columns[ds$int_col]) +
    4 * text_entries(c(formula_text, design_columns, weightsname, cp$rho_column,
                       cp$population, prep$clustervar))
  out_bound <- ceiling(65536 + 64 * (n_parameters + C) + label_entries +
    1024 * (M + C + n_pc + n_parameters) +
    sum(vapply(sb$links, function(l)
      512 + 64 * length(l$cohorts) + 2 * text_entries(l$link_id), numeric(1))))
  ra <- res$resource_accounting
  live <- ceiling(max(c(0, ra$simultaneous_entries, ra$initial_envelope_entries)))
  prior_work <- if (is.null(ra$downstream_work_ops)) 0 else ceiling(ra$downstream_work_ops)
  metadata_work <- 8 * out_bound
  metadata_scratch <- out_bound
  metadata_budget <- tryCatch(
    cp04_budget(ctrl, c(live, view_entries, out_bound), metadata_scratch, c(prior_work, view_work, metadata_work),
                out_bound, "nuisance selection tables"),
    cdid_ipw_owner_bound_error = identity)
  if (inherits(metadata_budget, "error"))
    return(chained_nuisance_status("unavailable", "budget",
      paste("The optional coefficient tables exceed the declared resource limits:",
            conditionMessage(metadata_budget))))

  blk <- chained_nuisance_block_geometry(sb, up_map)
  p_max <- max(c(0L, lengths(up_map$theta_coords)))
  n_range <- st$cohort_end - st$cohort_start + 1

  # ---- grouping: original units, or the requested (validated, nested) cluster --
  clustered <- !is.null(prep$cluster_full)
  group <- NULL; n_groups <- N; cov_code <- NULL; cov_reason <- NULL
  if (clustered) {
    if (length(prep$cluster_full) != N || !identical(st$id, prep$frame_ids)) {
      cov_code <- "cluster_alignment"
      cov_reason <- "Cluster labels could not be aligned with the selection fitting roster by frame position."
    } else {
      # Reserve for at most N cluster totals before constructing the map.
      n_groups <- N
    }
  }

  # ---- preflight: existing budget rule on the backend's owner control --------
  B <- as.double(min(ctrl$max_chunk_rows, reader$control$max_chunk_rows, N))
  log_b <- 1 + ceiling(log2(N + 1))
  d <- ds$d; d_s <- ds$d_s
  work_blocks <- sum(vapply(seq_len(M), function(m) {
    info <- blk[[m]]; B_m <- length(info$cids)
    rows_m <- sum(n_range[info$coh])
    per_row <- 32 * (B_m * log_b + B_m * (d_s + 2)^2 + B_m * (d + 4) + info$p + 8)
    rows_m * per_row + 64 * (info$p + 1)^2 + if (clustered) 32 * n_groups * info$p else 0
  }, numeric(1)))
  work_pairs <- if (is_date) sum(vapply(sb$pair_cohorts, function(pc)
    n_range[pc$cohort] * 32 * (2 * (d + 4) + 8), numeric(1))) else 0
  work <- ceiling(metadata_work + work_blocks + work_pairs)
  scratch <- max(metadata_scratch,
    ceiling(32 * B * (p_max + d + 8) +
      (if (clustered) 3 * as.double(N) * p_max + 4 * N else 0) + 8192))
  budget <- tryCatch(cp04_budget(ctrl, c(live, view_entries, out_bound), scratch, c(prior_work, view_work, work), out_bound,
    "nuisance selection capture"), cdid_ipw_owner_bound_error = function(e) e)
  do_pass <- !inherits(budget, "error")
  if (!do_pass) {
    cov_code <- "budget"
    cov_reason <- paste0("The bounded diagnostic/covariance pass exceeds the declared workspace or work budget: ",
      conditionMessage(budget), " Coefficients are reported; SEs and fitted-probability diagnostics are not.")
  }
  if (!n_parameters) {
    cov_code <- "fixed_probability"
    cov_reason <- "The selection probability is fixed at 1/2. There are no estimated coefficients."
  }
  if (is.null(cov_reason) && !is.null(src$cov_code)) {
    cov_code <- src$cov_code
    cov_reason <- src$cov_reason
  }
  do_cov <- do_pass && is.null(cov_reason)
  if (clustered && !do_cov) n_groups <- NA_integer_
  if (do_cov && clustered) {
    group <- match(prep$cluster_full, unique(prep$cluster_full))
    n_groups <- max(group)
    if (n_groups < 2L) {
      cov_code <- "insufficient_clusters"
      cov_reason <- "At least two clusters are needed for a clustered coefficient standard error."
      group <- NULL
      do_cov <- FALSE
    }
  }

  # ---- bounded passes through the estimator's own selection operation --------
  stats <- if (do_pass) {
    h <- src$open(blk)
    on.exit(h$abort(), add = TRUE)
    out <- chained_nuisance_selection_passes(sb, blk, st, B, N, h$fns, group, n_groups, do_cov = do_cov)
    h$finish()
    out
  } else chained_nuisance_selection_passes(sb, blk, st, B, N, NULL, run = FALSE)

  weights_text <- if (is.null(weightsname)) "w_i = 1 (no weightsname)" else
    if (identical(prep$weight_reduction, "first_period_static"))
      paste0("w_i = ", weightsname, " at the first calendar date (fix_weights = \"first_period\")") else
      paste0("w_i = ", weightsname)
  weights_use <- src$weights_use
  covariance <- chained_nuisance_covariance_record(do_cov, cov_code, cov_reason,
                                                   clustered, prep$clustervar, n_groups, N, src$if_source)
  if (do_cov && clustered) {
    admitted <- !is.na(stats$block_groups) & stats$block_groups >= 2L
    if (!all(admitted)) {
      covariance$status <- if (any(admitted)) "partial" else "unavailable"
      covariance$code <- "insufficient_clusters"
      covariance$reason <- "Each fitted block needs at least two clusters in its estimation sample."
    }
  }
  meta <- list(
    cohorts = st$cohorts, mass_U = st$mass_U, N = N, design_columns = design_columns, universe = sb$universe,
    formula = formula_text,
    population = if (is.null(cp$population)) NA_character_ else cp$population,
    weights = list(text = weights_text, column = weightsname, reduction = prep$weight_reduction, estimated = FALSE,
      note = paste("Fixed known sampling weights supplied by the caller, used as given (never estimated or self-normalized);",
        "the coefficients and their SEs are invariant to rescaling them.", weights_use)),
    transport = if (is.null(cp$rho_column)) "rho = 1" else
      paste0("known normalized rho column ", cp$rho_column, " (not in the logit fit)"),
    glm_control = src$glm_control,
    route = src$route, role = src$role, use = src$use,
    covariance = covariance,
    accounting = list(chunk_rows = B, work_ops_bound = work + view_work,
      work_ops_used_bound = view_work + if (do_pass) work else metadata_work,
      scratch_entries_bound = scratch,
      output_entries_bound = out_bound, view_entries_bound = view_entries, passes_run = do_pass,
      rule = paste("cp04_budget on the selection backend's owner control, composed with the pipeline's recorded live entries",
        "and downstream work (the public nuisance-diagnostics rule)")))
  comp <- chained_nuisance_selection_tables(sb, blk, stats, meta)
  actual <- unname(cp04_payload(comp)["entries"])
  if (actual > out_bound)
    return(chained_nuisance_status("unavailable", "budget",
      "The coefficient tables exceeded their reserved output size."))
  comp$accounting$output_entries <- actual
  comp
}

# ------------------------------------------------------------------------------
# 5b. Legacy separate pair-logit consumer: compact view of the actual records
# ------------------------------------------------------------------------------
# Direct comparison -> fitted (pair, cohort) logits, from the records the estimator
# reads. IPW: link$model -> treatment_specs[[model]], whose cohorts (target
# first, then sorted controls) align with its selection_ids. OR: the target
# support link$support_id and the control pool
# regression_specs[[link$regression]]$support_ids; support j is fitted as
# selection_specs[[j]] (fit_or_corollary1). NULL when an index is out of range
# or a mapped fit has another pair/cohort.
chained_nuisance_legacy_links <- function(b) {
  ss <- b$selection_specs; S <- length(ss)
  ok_id <- function(x, n) is.numeric(x) && length(x) >= 1L && !anyNA(x) && all(x >= 1 & x <= n & x == floor(x))
  st_pair <- function(pr) as.double(c(pr$s, pr$t))
  same <- function(j, pair, cohort) identical(st_pair(ss[[j]]$pair), st_pair(pair)) &&
    identical(as.integer(ss[[j]]$cohort), as.integer(cohort))
  is_ipw <- inherits(b, "cdid_cp04_corollary1")
  if (!is_ipw && (!is.list(b$support_specs) || length(b$support_specs) != S || !is.list(b$regression_specs)))
    return(NULL)
  out <- vector("list", length(b$links))
  for (i in seq_along(b$links)) {
    l <- b$links[[i]]
    if (!is.character(l$link_id) || length(l$link_id) != 1L || is.na(l$link_id)) return(NULL)
    if (is_ipw) {
      if (!ok_id(l$model, length(b$treatment_specs))) return(NULL)
      ts <- b$treatment_specs[[l$model]]
      sid <- as.integer(ts$selection_ids); coh <- as.integer(ts$cohorts); pair <- ts$pair
      if (!ok_id(sid, S) || length(sid) != length(coh)) return(NULL)
    } else {
      if (!ok_id(l$support_id, S) || !ok_id(l$regression, length(b$regression_specs))) return(NULL)
      ctl <- b$regression_specs[[l$regression]]$support_ids
      if (!ok_id(ctl, S)) return(NULL)
      sid <- unique(as.integer(c(l$support_id, ctl)))
      pair <- b$support_specs[[l$support_id]]$pair
      coh <- vapply(sid, function(j) as.integer(b$support_specs[[j]]$cohort), integer(1))
      if (!all(vapply(sid, function(j) identical(st_pair(b$support_specs[[j]]$pair), st_pair(pair)), logical(1))))
        return(NULL)
    }
    if (!all(vapply(seq_along(sid), function(k) same(sid[k], pair, coh[k]), logical(1)))) return(NULL)
    out[[i]] <- list(link_id = l$link_id, cohorts = coh, pair = pair, pair_cohort_ids = sid)
  }
  out
}

# The legacy fits as single-cell blocks in the shared-backend record shape, so
# the shared passes and tables apply unchanged: block = cell = pair-cohort = fit
# j, in fit order. Each coefficient vector covers every reader design column in
# Z order; the "(Intercept)" column of a design with an intercept is shown as
# the cell's local intercept, the other columns as its slopes. coord_order maps
# the Z-order IF columns to that (intercept, slopes) order.
chained_nuisance_legacy_view <- function(specs, fits, cohorts, design_columns, has_intercept, links, universe) {
  S <- length(specs); d <- length(design_columns)
  int_col <- if (isTRUE(has_intercept)) which(design_columns == "(Intercept)") else integer(0)
  slope_cols <- setdiff(seq_len(d), int_col)
  has_int <- length(int_col) == 1L
  branch <- if (has_int && length(slope_cols)) "intercept_slope" else if (has_int) "intercept_only" else "no_intercept"
  ds <- list(d = d, d_s = length(slope_cols), has_intercept = has_int, int_col = int_col,
    slope_cols = slope_cols, slope_names = design_columns[slope_cols], branch = branch)
  cells <- blocks <- pcs <- vector("list", S)
  for (j in seq_len(S)) {
    k <- as.integer(specs[[j]]$cohort); f <- fits[[j]]; cf <- as.double(f$coefficients)
    cells[[j]] <- list(type = "pair", cohort = k, cohort_value = cohorts[k], pair = specs[[j]]$pair,
      block_id = j, local_b = 1L)
    blocks[[j]] <- list(cell_ids = j, alpha = cf[int_col], beta = cf[slope_cols], branch = branch,
      has_intercept = has_int, d_s = ds$d_s, diagnostics = f$diagnostics, pi = f$pi)
    pcs[[j]] <- list(cohort = k, cohort_value = cohorts[k], pair = specs[[j]]$pair, cell_ids = j,
      block_ids = j, pair_cell_id = j)
  }
  list(sel = list(selection_model = "pair_logit", selection_pooling = "separate", universe = universe,
      design_spec = ds, cells = cells, blocks = blocks, pair_cohorts = pcs, links = links),
    up_map = list(theta_coords = lapply(seq_len(S), function(j) (j - 1L) * d + seq_len(d))),
    coord_order = c(int_col, slope_cols))
}

# The legacy consumer's own selection IF columns of one chunk of fitting rows,
# as evaluate_cp04_if_chunk(): xi = (H/N)^{-1} w (S - p) Z with the fit's stored
# N-normalized factor H_chol, in coord_order.
chained_nuisance_legacy_influence <- function(Z_rows, w, y, p, H_chol, ord) {
  score <- w * (y - p) * Z_rows
  xi <- t(cp04_solve(H_chol, t(score)))
  xi[, ord, drop = FALSE]
}

# Pass callbacks for chained_nuisance_selection_passes(). predict(rows, j) and
# respond(rows, j) give fit j's probability and pair response on its fitting
# rows; each chunk is evaluated once and shared by the IF and the cell summary.
chained_nuisance_legacy_fns <- function(Z, w, factors, ord, predict, respond) {
  last <- NULL
  at <- function(rows, j) {
    if (is.null(last) || last$j != j || !identical(last$rows, rows))
      last <<- list(j = j, rows = rows, p = predict(rows, j), y = respond(rows, j))
    last
  }
  list(
    influence = function(rows, m) {
      r <- at(rows, m)
      chained_nuisance_legacy_influence(Z[rows, , drop = FALSE], w[rows], r$y, r$p, factors[[m]], ord)
    },
    cell = function(rows, cid) {
      r <- at(rows, cid)
      list(p = r$p, y = r$y)
    },
    pair = function(rows, q) at(rows, q)$p)
}

chained_nuisance_legacy_component <- function(res, b, prep, weightsname) {
  is_ipw <- inherits(b, "cdid_cp04_corollary1")
  if (is_ipw) validate_cp04_backend(b) else validate_or_corollary1_backend(b)
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st <- reader$state
  ctrl <- b$control
  S <- length(b$selection_specs); L <- length(b$links); d <- ncol(st$Z)
  # Admit the compact view (metadata only) before building it.
  ra <- res$resource_accounting
  live <- ceiling(max(c(0, ra$simultaneous_entries, ra$initial_envelope_entries)))
  prior_work <- if (is.null(ra$downstream_work_ops)) 0 else ceiling(ra$downstream_work_ops)
  view_bound <- ceiling(4096 + 1024 * (S + L) + 64 * S * (d + 1) + 64 * L * (length(st$cohorts) + 1))
  factor_entries <- as.double(d)^2
  view_work <- ceiling(8 * view_bound + 4 * S * factor_entries)
  admitted <- tryCatch(cp04_budget(ctrl, c(live, view_bound), max(view_bound, 4 * factor_entries),
    c(prior_work, view_work), 0, "nuisance selection view"), cdid_ipw_owner_bound_error = identity)
  if (inherits(admitted, "error"))
    return(chained_nuisance_status("unavailable", "budget",
      paste("The optional coefficient tables exceed the declared resource limits:", conditionMessage(admitted))))
  links <- chained_nuisance_legacy_links(b)
  if (is.null(links))
    return(chained_nuisance_status("unavailable", "unrecognized_backend",
      "The legacy backend's comparisons could not be mapped to its fitted (pair, cohort) logits."))
  cn <- colnames(st$Z)
  hi <- attr(st$Z, "has_intercept")
  if (is.null(hi)) hi <- attr(st$Z, "intercept")
  if (is.null(hi)) hi <- "(Intercept)" %in% cn
  if (is.null(cn) || length(cn) != d || anyNA(cn) || (isTRUE(hi) && sum(cn == "(Intercept)") != 1L))
    return(chained_nuisance_status("unavailable", "unrecognized_backend",
      "The legacy selection design does not carry identifiable column labels and intercept column."))
  view <- chained_nuisance_legacy_view(b$selection_specs, b$selections, st$cohorts, cn, hi, links, st$universe)
  factors <- lapply(b$selections, `[[`, "H_chol")
  factor_ok <- all(vapply(factors, function(R) is.matrix(R) && is.numeric(R) && identical(dim(R), c(d, d)) &&
    all(is.finite(R)) && all(diag(R) > 0), logical(1)))
  open <- function(blk) list(abort = function() invisible(NULL), finish = function() invisible(NULL),
    fns = chained_nuisance_legacy_fns(st$Z, st$w, factors, view$coord_order,
      predict = function(rows, j) cp04_predict_selection(reader, rows, b$selection_specs[[j]], b$selections[[j]]),
      respond = function(rows, j) cp04_pair_rows(reader, rows, b$selection_specs[[j]]$pair)$S))
  em <- if (is_ipw) "ipw" else "reg"
  src <- list(sel = view$sel, up_map = view$up_map, reader = reader, ctrl = ctrl, glm_control = b$glm_control,
    open = open,
    route = paste0(class(b)[1L], " (legacy per-(pair, cohort) separate pair-logit consumer, est_method = \"", em, "\")"),
    role = "observation-selection propensity (conditional_mtar)",
    use = if (is_ipw)
      paste("p, this fit's pair probability, enters the observed-pair weight w a_g / (rho max(p, eps)) of each listed",
        "comparison: the fitting weights of its treatment logit and its weighted treated/control means") else
      paste("p, this fit's pair probability, enters q = w a_g / (rho max(p, eps)) on the observed pairs of its",
        "(pair, cohort) support: the control outcome-change regression weights (control pool) and the treated",
        "residual ratio (target cohort) of each listed direct comparison.",
        if (identical(b$outcome_pooling, "pooled"))
          "The control weights also enter the common slope across comparisons." else ""),
    weights_use = paste("The repeat rates a_k and the transport ratio rho enter the downstream weight",
      "a_k / (rho max(p, eps)), not this logit fit, so the reported SEs are not conditional on any estimated weight."),
    if_source = paste("Per-unit IF columns of the legacy consumer's own selection IF, (H/N)^{-1} w_i (S_i - p_i) Z_i",
      "with each fit's stored factor H_chol (as evaluate_cp04_if_chunk); no generic inverse Hessian and no second logit."),
    cov_code = if (!factor_ok) "factor_unavailable",
    cov_reason = if (!factor_ok) "A fit's stored information factor is not a finite positive-diagonal d x d Cholesky factor.")
  chained_nuisance_selection_run(res, src, prep, weightsname, view_bound, view_work)
}

# ------------------------------------------------------------------------------
# 6. Public accessor
# ------------------------------------------------------------------------------
nuisance_models <- function(fit, component = NULL, block = NULL) {
  if (!is.null(component) && !(is.character(component) && length(component) == 1L && !is.na(component) &&
      component %in% chained_nuisance_components))
    cdid_abort(paste0("`component` must be NULL or one of ", paste0("\"", chained_nuisance_components, "\"", collapse = ", "), "."),
      "cdid_nuisance_input_error")
  if (!is.null(block) && is.null(component)) component <- "selection"
  status <- chained_nuisance_fit_status(fit)
  if (!is.null(status)) return(status)
  cap <- fit$chained$nuisance_models
  if (is.null(component)) return(chained_nuisance_index(cap))
  comp <- cap$components[[component]]
  if (!identical(comp$status, "available"))
    return(structure(c(list(component = component), comp), class = "cdid_nuisance_status"))
  if (is.null(block)) return(structure(c(list(component = component), comp), class = "cdid_nuisance_component"))
  m <- chained_nuisance_block_position(comp, block)
  if (identical(component, "selection")) chained_nuisance_block_view(comp, m) else
    chained_nuisance_role_block_view(component, comp, m)
}

chained_nuisance_fit_status <- function(fit) {
  mk <- function(status, code, reason) structure(list(component = NA_character_, status = status, code = code, reason = reason),
    class = "cdid_nuisance_status")
  if (inherits(fit, c("cdid_MP", "cdid_static_att_gt"))) {
    cap <- if (is.list(fit$chained)) fit$chained$nuisance_models else NULL
    if (is.null(cap)) return(mk("not_captured", "not_requested", paste(
      "This fit carries no nuisance capture: it was estimated without control = list(keep_nuisance = TRUE)",
      "(the default) or by a cdid version without capture. Refit with control = list(keep_nuisance = TRUE);",
      "a saved call is never replayed, because changed data would silently describe a different fit.")))
    if (!is.list(cap) || !identical(cap$schema, chained_nuisance_schema))
      return(mk("unavailable", "unknown_capture_schema", "The nuisance capture has an unknown schema."))
    return(NULL)
  }
  if (inherits(fit, c("cdid_AGGTEobj", "AGGTEobj")))
    return(mk("unavailable", "aggregation_object",
      "Aggregations reuse the att_gt() fit without refitting; call nuisance_models() on the att_gt() result."))
  if (inherits(fit, "MP"))
    return(mk("unavailable", "ordinary_delegation",
      "mode = \"ordinary\" delegates to did::att_gt(); cdid captures no nuisance models for ordinary delegation."))
  cdid_abort("`fit` must be a result of cdid::att_gt().", "cdid_nuisance_input_error")
}

chained_nuisance_index <- function(cap) {
  comps <- cap$components
  pick <- function(x, nm) if (is.null(x[[nm]])) "" else as.character(x[[nm]])
  status_df <- data.frame(component = names(comps),
    status = vapply(comps, pick, character(1), "status"), code = vapply(comps, pick, character(1), "code"),
    reason = vapply(comps, pick, character(1), "reason"), stringsAsFactors = FALSE)
  rownames(status_df) <- NULL
  sel <- comps$selection
  avail <- identical(sel$status, "available")
  sel_fields <- c("selection_model", "selection_pooling", "link", "response", "template", "derived",
    "formula", "slope_terms", "universe_text", "weights", "covariance", "route")
  if (avail && !is.null(sel$use)) sel_fields <- c(sel_fields, "role", "use")
  roles <- Filter(function(x) identical(x$status, "available"), comps[setdiff(names(comps), "selection")])
  structure(list(
    est_method = cap$est_method, selection_est_method = cap$selection_est_method,
    components = status_df,
    selection = if (avail) sel[sel_fields] else NULL,
    blocks = if (avail) sel$blocks[, c("block", "label", "n_cells", "n_parameters", "n_fitting_rows", "n_units",
      "converged", "n_comparisons")] else NULL,
    roles = if (length(roles)) lapply(roles, function(r) c(r[c("role", "pooling", "template", "response", "sample",
      "weights", "covariance")], list(blocks = r$blocks[, c("block", "label", "n_parameters", "n_fitting_rows",
      "n_comparisons")]))) else NULL),
    class = "cdid_nuisance_index")
}

chained_nuisance_block_position <- function(comp, block) {
  ids <- comp$blocks$block
  m <- NA_integer_
  if (is.numeric(block) && length(block) == 1L && !is.na(block) && block == floor(block) &&
      block >= 1 && block <= length(ids)) m <- as.integer(block)
  if (is.character(block) && length(block) == 1L && !is.na(block)) {
    m <- match(block, ids)
    if (is.na(m)) m <- match(block, comp$blocks$label)
  }
  if (is.na(m)) cdid_abort(paste0("Unknown `block`. Use a block id (", paste(ids, collapse = ", "),
    "), its integer position or its label."), "cdid_nuisance_input_error")
  m
}

chained_nuisance_block_view <- function(comp, m) {
  bl <- comp$blocks[m, , drop = FALSE]
  co <- comp$coefficients[comp$coefficients$block == bl$block, , drop = FALSE]
  cm <- cbind(Estimate = co$estimate, `Std. Error` = co$std_error, `z value` = co$z_value, `Pr(>|z|)` = co$p_value)
  rownames(cm) <- co$term
  cells <- comp$cells[comp$cells$block == bl$block, , drop = FALSE]
  pcs <- comp$pair_cohorts[vapply(strsplit(comp$pair_cohorts$blocks, ",", fixed = TRUE),
    function(x) bl$block %in% x, logical(1)), , drop = FALSE]
  slopes <- co$term[co$role %in% c("shared_slope", "slope")]
  covariance <- comp$covariance
  if (!is.na(bl$n_fitting_clusters)) {
    covariance$n_fitting_clusters <- bl$n_fitting_clusters
    if (bl$n_fitting_clusters >= 2L && identical(covariance$status, "partial")) {
      covariance$status <- "available"
      covariance$code <- "ok"
      covariance$reason <- NULL
    } else if (bl$n_fitting_clusters < 2L) {
      covariance$status <- "unavailable"
      covariance$code <- "insufficient_clusters"
      covariance$reason <- "This block was fitted in fewer than two clusters."
    }
  }
  view <- list(
    block = bl$block, label = bl$label, selection_model = comp$selection_model,
    selection_pooling = comp$selection_pooling, response = comp$response, equation = bl$equation,
    intercept_cells = if (any(co$role == "local_intercept"))
      cells$key[match(co$cell[co$role == "local_intercept"], cells$cell)] else character(0),
    slope_terms = slopes,
    sample = paste0(if (is.null(bl$population))
      paste0("units with G_i in {", bl$cohorts, "} and U_i = 1 (", comp$universe_text, ")") else bl$population,
      ": ", if (!is.na(bl$n_units)) paste0(format(bl$n_units), " units, "),
      format(bl$n_fitting_rows), " fitting rows (units x cells)"),
    weights = comp$weights$text, derived = comp$derived,
    coefficients = cm, table = co, covariance = covariance,
    diagnostics = bl[, c("branch", "converged", "iterations", "rank", "n_cells", "n_parameters", "n_comparisons")],
    cells = cells, pair_cohorts = pcs, probability_floor = comp$probability_floor,
    pair_scope = if (identical(comp$route, "cohort_time_fe") || is.null(bl$population)) "cohort" else "role")
  if (!is.null(comp$use)) view[c("role", "use")] <- list(comp$role, comp$use)
  structure(view, class = "cdid_nuisance_block")
}

# ------------------------------------------------------------------------------
# 7. Printing
# ------------------------------------------------------------------------------
chained_nuisance_se_line <- function(cv) {
  if (!identical(cv$status, "available"))
    return(paste0("SE: ", cv$status, " (", cv$code, "): ", cv$reason))
  if (identical(cv$type, "cluster"))
    paste0("SE: cluster-robust by ", cv$cluster_var, " (", format(cv$n_groups), " frame clusters",
           if (!is.null(cv$n_fitting_clusters)) paste0(", ", cv$n_fitting_clusters, " in this block") else "",
           if (!is.null(cv$n_contributing_clusters))
             paste0(", ", cv$n_contributing_clusters, " with nonzero unit contributions") else "",
           "): unit scores summed within cluster, ",
           "centered, / N^2; no small-sample factor.") else
    paste0("SE: sandwich with one score contribution per original unit (N = ", format(cv$N), "), centered, / N^2; ",
           "no small-sample factor.")
}

print.cdid_nuisance_status <- function(x, ...) {
  lead <- if (is.na(x$component)) "Nuisance models" else paste0("Nuisance component \"", x$component, "\"")
  cat(lead, ": ", x$status, " (", x$code, ")\n", sep = "")
  cat(strwrap(x$reason, indent = 2L, exdent = 2L), sep = "\n")
  invisible(x)
}

print.cdid_nuisance_index <- function(x, ...) {
  cat("Nuisance models captured from att_gt() [est_method = \"", x$est_method, "\", selection_est_method = \"",
      x$selection_est_method, "\"]\n", sep = "")
  if (!is.null(x$selection)) {
    s <- x$selection
    cat("\nSelection (observation) logits: ", s$selection_model, ", ", s$selection_pooling, " pooling, link = ", s$link, "\n", sep = "")
    cat("  ", s$template, "\n", sep = "")
    if (!is.null(s$use)) cat("  role:     ", s$role, "\n  use:      ", s$use, "\n", sep = "")
    cat("  response: ", s$response, "\n", sep = "")
    if (!is.null(s$derived)) cat("  derived:  ", s$derived, "\n", sep = "")
    cat("  formula:  ", s$formula, "   slopes: ",
        if (length(s$slope_terms)) paste(s$slope_terms, collapse = ", ") else "none", "\n", sep = "")
    cat("  sample:   ", s$universe_text, "; weights: ", s$weights$text, "\n", sep = "")
    cat("  ", chained_nuisance_se_line(s$covariance), "\n\n", sep = "")
    print(x$blocks, row.names = FALSE)
  }
  for (nm_r in names(x$roles)) {
    r <- x$roles[[nm_r]]
    cat("\n", chained_nuisance_aug_title[[nm_r]], " (", r$pooling, "): ", r$role, "\n", sep = "")
    cat("  ", r$template, "\n", sep = "")
    cat("  response: ", r$response, "\n", sep = "")
    cat("  sample:   ", r$sample, "; weights: ", r$weights$text, "\n", sep = "")
    cat("  ", chained_nuisance_se_line(r$covariance), "\n\n", sep = "")
    print(r$blocks, row.names = FALSE)
  }
  other <- x$components[x$components$status != "available", , drop = FALSE]
  if (nrow(other)) {
    cat("\nOther components:\n")
    for (i in seq_len(nrow(other)))
      cat("  ", other$component[i], ": ", other$status[i], " - ", other$reason[i], "\n", sep = "")
  }
  if (!is.null(x$selection))
    cat("\nCoefficient tables: nuisance_models(fit, \"selection\") or nuisance_models(fit, \"selection\", block = \"",
        x$blocks$block[1L], "\").\n", sep = "")
  invisible(x)
}

chained_nuisance_print_block <- function(x, ...) {
  cat("Selection block ", x$block, ": ", x$label, "  (", x$selection_model, ", ", x$selection_pooling, ")\n", sep = "")
  cat("  ", x$equation, "\n", sep = "")
  cat("  response: ", x$response, "\n", sep = "")
  if (length(x$intercept_cells)) cat("  intercept cells: ", paste(x$intercept_cells, collapse = "; "), "\n", sep = "")
  cat("  slopes: ", if (length(x$slope_terms)) paste(x$slope_terms, collapse = ", ") else "none",
      if (length(x$slope_terms) && nrow(x$cells) > 1L) paste0(" (shared by the ", nrow(x$cells), " cells)") else "", "\n", sep = "")
  cat("  sample: ", x$sample, "\n", sep = "")
  cat("  weights: ", x$weights, "\n", sep = "")
  if (!is.null(x$derived)) cat("  derived: ", x$derived, "\n", sep = "")
  if (!is.null(x$use)) cat("  role: ", x$role, "\n  use: ", x$use, "\n", sep = "")
  if (nrow(x$coefficients)) {
    cat("\nCoefficients:\n")
    stats::printCoefmat(x$coefficients, has.Pvalue = TRUE, P.values = TRUE, na.print = "NA", ...)
  } else cat("\nNo fitted coefficients (zero-column design: Pr = 1/2).\n")
  cat("\n", chained_nuisance_se_line(x$covariance), "\n", sep = "")
  dg <- x$diagnostics
  cat("converged = ", dg$converged, " (", format(dg$iterations), " iterations), rank = ", format(dg$rank), "\n", sep = "")
  low <- x$cells$p_min[is.finite(x$cells$p_min)]
  high <- x$cells$p_max[is.finite(x$cells$p_max)]
  if (length(low) && length(high)) {
    cat("fitted Pr on the fitting sample in [", format(min(low), digits = 4), ", ",
        format(max(high), digits = 4), "]\n", sep = "")
  }
  if (any(is.finite(x$cells$n_floor))) {
    cat(if (identical(x$selection_model, "independent_period_logit")) {
      paste0("date probabilities at or below ", format(x$probability_floor), ": ")
    } else paste0("rows at the ", format(x$probability_floor), " floor: "),
    format(sum(x$cells$n_floor, na.rm = TRUE)), "\n", sep = "")
  }
  cat("used by ", nrow(x$pair_cohorts), " pair-", x$pair_scope, " probabilit", if (nrow(x$pair_cohorts) == 1L) "y" else "ies",
      " and ", format(dg$n_comparisons), " comparison(s)\n", sep = "")
  invisible(x)
}

print.cdid_nuisance_block <- function(x, ...) {
  if (!is.null(x$component) && !identical(x$component, "selection")) return(chained_nuisance_print_role_block(x, ...))
  chained_nuisance_print_block(x, ...)
}

print.cdid_nuisance_component <- function(x, ...) {
  if (!is.null(x$component) && !identical(x$component, "selection")) return(chained_nuisance_print_role_component(x, ...))
  cat("Selection (observation) logits: ", x$selection_model, ", ", x$selection_pooling, " pooling\n", sep = "")
  cat("  ", x$template, "\n\n", sep = "")
  for (m in seq_len(nrow(x$blocks))) {
    chained_nuisance_print_block(chained_nuisance_block_view(x, m), ...)
    cat("\n")
  }
  invisible(x)
}

# ------------------------------------------------------------------------------
# 8. Augmented DR+DR: outcome regression, treatment propensity and imputation
# ------------------------------------------------------------------------------
# Admitted scope: selection = "conditional_mtar", selection_est_method = "dr",
# est_method = "dr" on the ordinary static route (the predicate of
# chained_nuisance_augmented_scope(): no reference/layout/recipe source, no
# private pair differences, U = F, rho = 1). Read from the sealed
# cdid_augmented_backend records the estimator uses (no refit):
#   imputation  b$groups[[r]]: one WLS per imputation pooling group r of
#               Delta_e on observed pairs (S_e = 1) of its cells, weights w_i,
#               local cell intercepts (when the formula has one) and slopes on
#               L(s,t) z(X) over the group's retained basis;
#   outcome     b$shared_or: WLS of Z_c on every F row of each control cohort k
#               of each deduplicated control block, weights w_i a_k, local block
#               intercepts and L(s,t) z(X) slopes, either ONE pooled fit (common
#               retained slope, shown once) or one fit per block ("separate");
#   treatment   b$treatments: binary cohort risk-set logit 1{G = g} on every F
#               row of {g} U K, weights w_i a_{G_i}; specs with the same target
#               and risk set reuse one fit (shown once with its reuse list).
# Coefficient SEs reuse the estimator's own original-unit coefficient IF columns
# of evaluate_augmented_if_chunk(b, rows, links = 1L, details = TRUE):
# imputation_alpha/imputation_beta, or_if$theta (pooled) or
# or_if$blocks[[b]]$theta (separate; each the sum of the direct, a, j and q
# channels) and trt_if[[t]]$gamma (direct plus repeat-rate channel). details =
# TRUE evaluates every nuisance fit whatever the link batch, so one comparison
# suffices. One bounded pass over all N frame rows accumulates the diagonal with
# the selection convention (chained_nuisance_moments_*).

chained_nuisance_aug_title <- c(outcome = "Outcome regression",
  treatment = "Treatment propensity logit", imputation = "Actual-change imputation regression")

chained_nuisance_pair_txt <- function(pair) paste0("(", chained_nuisance_num(pair$s), ",", chained_nuisance_num(pair$t), ")")

chained_nuisance_gtxt <- function(v) paste(chained_nuisance_num(chained_nuisance_g(v)), collapse = ",")

# Backend shape checks for the DR+DR roles (metadata only).
chained_nuisance_aug_backend_ok <- function(b) {
  sh <- b$shared_or
  identical(b$est_method, "dr") && is.list(b$groups) && length(b$groups) > 0L &&
    is.list(b$cells) && length(b$cells) > 0L && is.list(sh) && !is.null(sh$blocks) &&
    is.list(b$treatments) && length(b$treatments) > 0L && is.list(b$links) && length(b$links) > 0L &&
    inherits(b$outcome_recipe, "cdid_or_pooled_recipe") && is.list(b$treatment_recipe) && b$treatment_recipe$p > 0L
}

# The static producer reuses the fit for the exact target and risk set, even
# when a later pair gets a distinct comparison ID.
chained_nuisance_aug_trt_reps <- function(trs) {
  key <- vapply(trs, function(tr) paste(tr$target_cohort_index, paste(tr$all_cohort_indices, collapse = ","),
    sep = "|"), character(1))
  match(key, key)
}

# Fit records of the three roles: estimates, terms, positions and the extractor
# of each fit's coefficient IF columns from the details of one chunk.
chained_nuisance_aug_fits <- function(b) {
  rc <- b$outcome_recipe; icpt <- isTRUE(rc$intercept)
  int_name <- if (icpt) rc$columns[rc$intercept_index] else NA_character_
  cell_key <- function(cf) paste0("g=", chained_nuisance_gtxt(cf$cohort), ",", chained_nuisance_pair_txt(cf$pair))
  blk_key <- function(bf) paste0(chained_nuisance_pair_txt(bf$pair), ",K={", chained_nuisance_gtxt(bf$control_cohorts), "}")
  imputation <- lapply(seq_along(b$groups), function(gid) {
    gf <- b$groups[[gid]]; cids <- as.integer(gf$cell_ids); p <- as.integer(gf$p)
    alpha <- if (icpt) as.double(gf$alpha) else numeric(0)
    if (icpt && length(alpha) != length(cids)) cp04_abort("Imputation group intercepts disagree with its cells.", "integrity")
    list(role = "imputation", id = paste0("I", gid), backend_id = gid, cids = cids,
      n_intercepts = length(alpha), n_slopes = p, estimate = c(alpha, as.double(gf$beta)),
      term = c(if (icpt) paste0(int_name, "[", vapply(b$cells[cids], cell_key, ""), "]"), as.character(gf$kept_columns)),
      coef_role = c(rep("local_intercept", length(alpha)), rep(if (length(cids) > 1L) "shared_slope" else "slope", p)),
      cell = c(if (icpt) cids, rep(NA_integer_, p)), kept = as.character(gf$kept_columns),
      aliased = as.character(gf$basis$aliased_columns),
      extract = function(d) cbind(if (icpt) d$imputation_alpha[, cids, drop = FALSE], d$imputation_beta[[gid]]))
  })
  sh <- b$shared_or; nb <- length(sh$blocks)
  outcome <- if (!reconstruction_separate(sh)) {
    alpha <- if (icpt) as.double(sh$alpha) else numeric(0)
    if (icpt && length(alpha) != nb) cp04_abort("Shared OR intercepts disagree with its blocks.", "integrity")
    list(list(role = "outcome", id = "O1", backend_id = 1L, bids = seq_len(nb), pooled = TRUE,
      n_intercepts = length(alpha), n_slopes = length(sh$beta), estimate = c(alpha, as.double(sh$beta)),
      term = c(if (icpt) paste0(int_name, "[", vapply(sh$blocks, blk_key, ""), "]"), as.character(sh$kept_columns)),
      coef_role = c(rep("local_intercept", length(alpha)), rep("shared_slope", length(sh$beta))),
      cell = c(if (icpt) seq_len(nb), rep(NA_integer_, length(sh$beta))), kept = as.character(sh$kept_columns),
      aliased = as.character(sh$basis$aliased_columns),
      extract = function(d) d$or_if$theta))
  } else lapply(seq_len(nb), function(bid) {
    loc <- reconstruction_or_block(sh, bid); bf <- sh$blocks[[bid]]
    alpha <- if (icpt) as.double(bf$alpha) else numeric(0)
    beta <- as.double(loc$beta)
    list(role = "outcome", id = paste0("O", bid), backend_id = bid, bids = bid, pooled = FALSE,
      n_intercepts = length(alpha), n_slopes = length(beta), estimate = c(alpha, beta),
      term = c(if (icpt) paste0(int_name, "[", blk_key(bf), "]"), as.character(loc$kept_columns)),
      coef_role = c(rep("local_intercept", length(alpha)), rep("slope", length(beta))),
      cell = c(if (icpt) bid, rep(NA_integer_, length(beta))), kept = as.character(loc$kept_columns),
      aliased = as.character(loc$basis$aliased_columns),
      extract = function(d) d$or_if$blocks[[as.character(bid)]]$theta)
  })
  trs <- b$treatments; rep_of <- chained_nuisance_aug_trt_reps(trs)
  reps <- which(rep_of == seq_along(trs))
  specs <- split(seq_along(trs), factor(rep_of, levels = reps))
  treatment <- lapply(seq_along(reps), function(m) {
    tid <- reps[m]; tr <- trs[[tid]]; cols <- as.character(tr$column_names)
    is_int <- if (icpt) cols == int_name else rep(FALSE, length(cols))
    list(role = "treatment", id = paste0("T", m), backend_id = tid, tids = specs[[m]],
      n_intercepts = sum(is_int), n_slopes = sum(!is_int), estimate = as.double(tr$coefficients), term = cols,
      coef_role = ifelse(is_int, "intercept", "slope"), cell = rep(NA_integer_, length(cols)), kept = cols,
      aliased = character(0), extract = function(d) d$trt_if[[tid]]$gamma)
  })
  list(outcome = outcome, treatment = treatment, imputation = imputation)
}

# One bounded pass over all N frame rows (original units, frame order). Returns
# per-fit variance vectors (NULL where unavailable) and contributing clusters.
# Containment: the R details path clears cdid_fast_kernel_cache$last_retained_if;
# that entry is restored, so the capture leaves the kernel cache as it found it.
chained_nuisance_aug_stream <- function(b, fits, N, B, group = NULL, n_groups = N) {
  cache <- get0("cdid_fast_kernel_cache", inherits = TRUE)
  if (is.environment(cache)) {
    had <- exists("last_retained_if", envir = cache, inherits = FALSE)
    old <- if (had) get("last_retained_if", envir = cache, inherits = FALSE) else NULL
    on.exit(if (had) assign("last_retained_if", old, envir = cache) else
      if (exists("last_retained_if", envir = cache, inherits = FALSE)) rm("last_retained_if", envir = cache), add = TRUE)
  }
  p <- vapply(fits, function(f) length(f$estimate), integer(1))
  W <- sum(p); end <- cumsum(p); cols <- lapply(seq_along(fits), function(f) seq_len(p[f]) + end[f] - p[f])
  acc <- chained_nuisance_moments_new(W, if (is.null(group)) NULL else n_groups)
  seen <- if (is.null(group)) NULL else matrix(FALSE, n_groups, length(fits))
  for (start in seq.int(1L, N, by = B)) {
    rr <- seq.int(start, min(N, start + B - 1L))
    d <- evaluate_augmented_if_chunk(b, rr, 1L, "frame", TRUE)$details
    xi <- matrix(0, length(rr), W)
    for (f in seq_along(fits)) {
      x <- fits[[f]]$extract(d)
      if (!is.matrix(x) || nrow(x) != length(rr) || ncol(x) != p[f])
        return(list(ok = FALSE, reason = paste0("The coefficient IF columns of fit ", fits[[f]]$id,
          " do not match its fitted coefficients.")))
      xi[, cols[[f]]] <- x
      if (!is.null(seen)) {
        hit <- rowSums(x != 0) > 0
        if (any(hit)) seen[unique(group[rr[hit]]), f] <- TRUE
      }
    }
    acc <- chained_nuisance_moments_add(acc, xi, if (is.null(group)) NULL else group[rr])
    rm(d, xi)
  }
  v <- chained_nuisance_moments_variance(acc, N, n_groups)
  list(ok = TRUE, var = lapply(cols, function(ix) v[ix]),
    clusters = if (is.null(seen)) rep(NA_integer_, length(fits)) else as.integer(colSums(seen)))
}

chained_nuisance_aug_covariance <- function(role, available, code, reason, clustered, cluster_var, n_groups, N) {
  src <- switch(role,
    outcome = paste("or_if$theta (pooled) or or_if$blocks[[b]]$theta (separate) of evaluate_augmented_if_chunk(..., details =",
      "TRUE): the direct WLS score of Z_c plus the repeat-rate a_k, imputation j and observation-propensity q channels,",
      "each pulled back through the fit's own Gram factor"),
    treatment = paste("trt_if[[t]]$gamma of evaluate_augmented_if_chunk(..., details = TRUE): the weighted logit score",
      "plus the repeat-rate a_k channel, solved with the fit's stored factor H_chol; no j or q path"),
    imputation = paste("imputation_alpha / imputation_beta of evaluate_augmented_if_chunk(..., details = TRUE): the",
      "observed-pair WLS score with fixed weights w_i through the group's own Gram factor; no estimated input"))
  accounts <- switch(role,
    outcome = "the estimated repeat rates a_k, the imputation fits j and the observation logits q that build Z_c and its weights",
    treatment = "the estimated repeat rates a_k in its weights",
    imputation = "nothing: its weights are fixed and its response is observed")
  list(status = if (available) "available" else "unavailable", code = if (available) "ok" else code,
    reason = if (available) NULL else reason, type = if (clustered) "cluster" else "unit", cluster_var = cluster_var,
    n_groups = n_groups, N = N,
    estimator = paste0("Centered sandwich: each original unit's coefficient IF contribution", if (clustered)
      ", summed within its cluster" else "", "; squared centered sums / N^2 over all N frame units (zero-contribution ",
      "units included), no small-sample factor. Marginal variances (diagonal of the joint sandwich) only."),
    source = src, accounts_for = paste("Accounts for", accounts), diagonal_only = TRUE)
}

chained_nuisance_aug_weights_text <- function(weightsname, prep) {
  if (is.null(weightsname)) return("w_i = 1 (no weightsname)")
  if (identical(prep$weight_reduction, "first_period_static"))
    paste0("w_i = ", weightsname, " at the first calendar date (fix_weights = \"first_period\")") else
    paste0("w_i = ", weightsname)
}

chained_nuisance_aug_equation <- function(f, b) {
  rc <- b$outcome_recipe; slopes <- chained_nuisance_slope_text(f$kept)
  sl <- function(L) {
    if (!length(slopes)) return(character())
    if (is.numeric(L) && L == 1) paste(slopes, collapse = " + ") else
      paste0(L, "*(", paste(slopes, collapse = " + "), ")")
  }
  join <- function(...) { x <- c(...); if (length(x)) paste(x, collapse = " + ") else "0" }
  if (identical(f$role, "imputation")) {
    cl <- b$cells[f$cids]; single <- length(cl) == 1L
    gs <- unique(vapply(cl, function(c) chained_nuisance_g(c$cohort), numeric(1)))
    g_txt <- if (length(gs) == 1L) chained_nuisance_num(gs) else "g"
    lhs <- if (single) paste0("E[Y_i", chained_nuisance_num(cl[[1L]]$pair$t), " - Y_i", chained_nuisance_num(cl[[1L]]$pair$s),
      " | X_i, G_i = ", g_txt, ", S_i", chained_nuisance_pair_txt(cl[[1L]]$pair), " = 1]") else
      paste0("E[Y_it - Y_is | X_i, G_i = ", g_txt, ", S_i(s,t) = 1]")
    a <- if (!f$n_intercepts) character(0) else if (single) "a" else if (length(gs) == 1L) "a[(s,t)]" else "a[g,(s,t)]"
    return(paste(lhs, "=", join(a, sl(if (single) cl[[1L]]$L else "L(s,t)"))))
  }
  if (identical(f$role, "outcome")) {
    bl <- b$shared_or$blocks[f$bids]
    if (length(bl) == 1L) {
      pr <- chained_nuisance_pair_txt(bl[[1L]]$pair)
      return(paste0("E[Z_i(k;", substr(pr, 2L, nchar(pr) - 1L), ") | X_i, G_i = k] = ",
        join(if (f$n_intercepts) "a", sl(bl[[1L]]$L)), ",  k in {", chained_nuisance_gtxt(bl[[1L]]$control_cohorts), "}"))
    }
    return(paste0("E[Z_i(k;s,t) | X_i, G_i = k] = ", join(if (f$n_intercepts) "a[b]", sl("L(s,t)")),
      ",  k in K_b for each control block b = ((s,t), K_b)"))
  }
  tr <- b$treatments[[f$backend_id]]
  cterm <- vapply(seq_along(f$term), function(j)
    if (f$coef_role[j] == "intercept") paste0("c", j - 1L) else paste0("c", j - 1L, "*", f$term[j]), "")
  if (length(cterm) > 7L) cterm <- c(cterm[1:6], paste0("... (", length(cterm) - 6L, " more terms)"))
  paste0("logit Pr(G_i = ", chained_nuisance_gtxt(tr$target_cohort), " | X_i, G_i in {",
    chained_nuisance_gtxt(c(tr$target_cohort, tr$control_cohorts)), "}) = ", paste(cterm, collapse = " + "))
}

chained_nuisance_aug_template <- function(role, b) {
  pool <- b$imputation_pooling
  switch(role,
    imputation = paste0("E[Y_it - Y_is | X_i, G_i = g, S_i(s,t) = 1] = ",
      if (isTRUE(b$outcome_recipe$intercept)) "a[g,(s,t)] + " else "", "L(s,t)*X_i'b",
      switch(pool, separate = "[g,(s,t)]", within_cohort = "[g]", across_cohorts = "", "[r]"),
      "  (imputation_pooling = \"", pool, "\")"),
    outcome = paste0("E[Z_i(k;s,t) | X_i, G_i = k] = ", if (isTRUE(b$outcome_recipe$intercept)) "a[b] + " else "",
      "L(s,t)*X_i'b", if (reconstruction_separate(b$shared_or)) "[b]" else "",
      ",  k in K_b  (outcome_pooling = \"", if (reconstruction_separate(b$shared_or)) "separate" else "pooled", "\")"),
    treatment = "logit Pr(G_i = g | X_i, G_i in {g} U K) = X_i'c[g,K]")
}

# Tables of one role from its fit records, variances and backend metadata (pure).
chained_nuisance_aug_role_tables <- function(role, fits, b, st, stats, meta) {
  N <- st$N; rows_k <- st$cohort_end - st$cohort_start + 1; wmass_k <- st$mass_F * N
  links <- b$links
  fit_for <- integer(switch(role, imputation = length(b$cells), outcome = length(b$shared_or$blocks),
    treatment = length(b$treatments)))
  for (m in seq_along(fits)) {
    f <- fits[[m]]
    fit_for[switch(role, imputation = f$cids, outcome = f$bids, treatment = f$tids)] <- m
  }
  or_counts <- if (role == "outcome")
    tabulate(vapply(links, `[[`, integer(1), "or_block_id"), length(fit_for))
  link_fits <- lapply(links, function(l) switch(role,
    imputation = sort(unique(fit_for[c(l$target_cell_id, l$control_cell_ids)])),
    outcome = fit_for[l$or_block_id], treatment = fit_for[l$trt_id]))
  comparison_counts <- tabulate(unlist(link_fits, use.names = FALSE), length(fits))
  ids <- vapply(fits, `[[`, "", "id")
  coefs <- do.call(rbind, c(list(chained_nuisance_coef_frame(character(0), character(0), character(0), integer(0),
    numeric(0), numeric(0))), lapply(seq_along(fits), function(m) {
      f <- fits[[m]]
      v <- if (is.null(stats$var[[m]])) rep(NA_real_, length(f$estimate)) else stats$var[[m]]
      chained_nuisance_coef_frame(f$id, f$term, f$coef_role, f$cell, f$estimate, v)
    })))
  rownames(coefs) <- NULL
  fit_rows <- function(f) switch(role,
    imputation = sum(vapply(b$cells[f$cids], function(c) as.double(c$n_pair), numeric(1))),
    outcome = sum(vapply(b$shared_or$blocks[f$bids], function(bf) sum(rows_k[bf$control_cohort_indices]), numeric(1))),
    treatment = as.double(b$treatments[[f$backend_id]]$n_risk_F))
  fit_w <- function(f) switch(role,
    imputation = sum(vapply(b$cells[f$cids], function(c) as.double(c$weight_sum), numeric(1))),
    outcome = sum(vapply(b$shared_or$blocks[f$bids], function(bf) as.double(bf$weight_sum), numeric(1))),
    treatment = { ci <- b$treatments[[f$backend_id]]$all_cohort_indices; sum(st$a[ci] * wmass_k[ci]) })
  fit_coh <- function(f) switch(role,
    imputation = chained_nuisance_gtxt(sort(unique(vapply(b$cells[f$cids], function(c) chained_nuisance_g(c$cohort), numeric(1))))),
    outcome = chained_nuisance_gtxt(sort(unique(unlist(lapply(b$shared_or$blocks[f$bids], function(bf)
      chained_nuisance_g(bf$control_cohorts)))))),
    treatment = { tr <- b$treatments[[f$backend_id]]; chained_nuisance_gtxt(c(tr$target_cohort, tr$control_cohorts)) })
  label <- function(f) switch(role,
    imputation = paste0(b$groups[[f$backend_id]]$pooling, "[", b$groups[[f$backend_id]]$name, ": ",
      paste(vapply(b$cells[f$cids], function(c) paste0("g=", chained_nuisance_gtxt(c$cohort), ",",
        chained_nuisance_pair_txt(c$pair)), ""), collapse = "; "), "]"),
    outcome = if (isTRUE(f$pooled)) paste0("pooled[", length(f$bids), " control block(s), common slope]") else
      paste0("block[", b$shared_or$blocks[[f$bids]]$name, "]"),
    treatment = { tr <- b$treatments[[f$backend_id]]; paste0("risk_set[g=", chained_nuisance_gtxt(tr$target_cohort),
      " vs K={", chained_nuisance_gtxt(tr$control_cohorts), "}]") })
  blocks_df <- data.frame(block = ids, label = vapply(fits, label, ""),
    backend_id = vapply(fits, function(f) as.integer(f$backend_id), integer(1)),
    cohorts = vapply(fits, fit_coh, ""),
    n_parameters = vapply(fits, function(f) length(f$estimate), integer(1)),
    n_intercepts = vapply(fits, function(f) as.integer(f$n_intercepts), integer(1)),
    n_slopes = vapply(fits, function(f) as.integer(f$n_slopes), integer(1)),
    kept_slopes = vapply(fits, function(f)
      paste(if (role == "treatment") f$term[f$coef_role == "slope"] else f$kept, collapse = ","), ""),
    aliased_slopes = vapply(fits, function(f) paste(f$aliased, collapse = ","), ""),
    n_fitting_rows = vapply(fits, fit_rows, numeric(1)),
    weight_sum = vapply(fits, fit_w, numeric(1)),
    n_contributing_clusters = stats$clusters,
    n_comparisons = comparison_counts,
    equation = vapply(fits, chained_nuisance_aug_equation, "", b = b),
    stringsAsFactors = FALSE)
  link_df <- data.frame(
    comparison = vapply(links, function(l) as.character(l$link_id), ""),
    g = vapply(links, function(l) chained_nuisance_g(l$target_cohort), numeric(1)),
    s = vapply(links, function(l) as.double(l$pair$s), numeric(1)),
    t = vapply(links, function(l) as.double(l$pair$t), numeric(1)),
    controls = vapply(links, function(l) chained_nuisance_gtxt(l$control_cohorts), ""),
    fits = vapply(link_fits, function(x) paste(ids[x], collapse = ","), ""),
    stringsAsFactors = FALSE)
  detail <- switch(role,
    imputation = data.frame(cell = seq_along(b$cells), block = ids[fit_for],
      cohort = vapply(b$cells, function(c) chained_nuisance_g(c$cohort), numeric(1)),
      s = vapply(b$cells, function(c) as.double(c$pair$s), numeric(1)),
      t = vapply(b$cells, function(c) as.double(c$pair$t), numeric(1)),
      L = vapply(b$cells, function(c) as.double(c$L), numeric(1)),
      n_cohort_rows = vapply(b$cells, function(c) as.double(rows_k[c$cohort_index]), numeric(1)),
      n_observed_pairs = vapply(b$cells, function(c) as.double(c$n_pair), numeric(1)),
      weight_sum = vapply(b$cells, function(c) as.double(c$weight_sum), numeric(1)),
      stringsAsFactors = FALSE),
    outcome = data.frame(control_block = seq_along(b$shared_or$blocks),
      block = ids[fit_for],
      name = vapply(b$shared_or$blocks, function(bf) as.character(bf$name), ""),
      s = vapply(b$shared_or$blocks, function(bf) as.double(bf$pair$s), numeric(1)),
      t = vapply(b$shared_or$blocks, function(bf) as.double(bf$pair$t), numeric(1)),
      L = vapply(b$shared_or$blocks, function(bf) as.double(bf$L), numeric(1)),
      controls = vapply(b$shared_or$blocks, function(bf) chained_nuisance_gtxt(bf$control_cohorts), ""),
      repeat_rates = vapply(b$shared_or$blocks, function(bf)
        paste(signif(st$a[bf$control_cohort_indices], 6), collapse = ","), ""),
      n_fitting_rows = vapply(b$shared_or$blocks, function(bf) sum(rows_k[bf$control_cohort_indices]), numeric(1)),
      weight_sum = vapply(b$shared_or$blocks, function(bf) as.double(bf$weight_sum), numeric(1)),
      n_comparisons = or_counts,
      stringsAsFactors = FALSE),
    treatment = data.frame(spec = seq_along(b$treatments),
      block = ids[fit_for],
      g = vapply(b$treatments, function(tr) chained_nuisance_g(tr$target_cohort), numeric(1)),
      controls = vapply(b$treatments, function(tr) chained_nuisance_gtxt(tr$control_cohorts), ""),
      s = vapply(b$treatments, function(tr) as.double(tr$pair$s), numeric(1)),
      t = vapply(b$treatments, function(tr) as.double(tr$pair$t), numeric(1)),
      repeat_rates = vapply(b$treatments, function(tr) paste(signif(st$a[tr$all_cohort_indices], 6), collapse = ","), ""),
      n_treated_rows = vapply(b$treatments, function(tr) as.double(rows_k[tr$target_cohort_index]), numeric(1)),
      n_comparisons = vapply(b$treatments, function(tr) length(tr$links), integer(1)),
      reuses_fit = duplicated(fit_for),
      stringsAsFactors = FALSE))
  rownames(detail) <- NULL
  rc <- b$outcome_recipe
  out <- list(status = "available", role = meta$role[[role]], route = meta$route,
    template = chained_nuisance_aug_template(role, b), response = meta$response[[role]],
    use = meta$use[[role]], formula = meta$formula,
    design_columns = if (role == "treatment") b$treatment_recipe$columns else rc$columns,
    intercept_term = if (isTRUE(rc$intercept)) rc$columns[rc$intercept_index] else NULL,
    slope_terms = rc$slope_columns,
    pooling = switch(role, imputation = b$imputation_pooling,
      outcome = if (reconstruction_separate(b$shared_or)) "separate" else "pooled", treatment = "risk_set"),
    sample = meta$sample[[role]], weights = meta$weights[[role]], n_frame = N,
    repeat_rates = stats::setNames(as.double(st$a), chained_nuisance_num(chained_nuisance_g(st$cohorts))),
    blocks = blocks_df, coefficients = coefs, comparisons = link_df,
    covariance = meta$covariance[[role]], accounting = meta$accounting)
  out[[switch(role, imputation = "cells", outcome = "control_blocks", treatment = "models")]] <- detail
  out
}

chained_nuisance_aug_meta_text <- function(weights_text) {
  list(
    role = list(
      outcome = "counterfactual control outcome-change regression on the reconstructed response Z_c (augmented DR+DR)",
      treatment = "roster treatment propensity logit (augmented DR+DR)",
      imputation = "actual-change imputation regression j_c (augmented DR+DR)"),
    response = list(
      outcome = paste("Z_i(k;s,t) = j_c + S_i(s,t) (Y_it - Y_is - j_c) / max(q_c, eps) on observed pairs, j_c otherwise,",
        "for the control cell c = (k,(s,t)) of every F unit of control cohort k"),
      treatment = "D_i = 1{G_i = g} on the risk set G_i in {g} U K",
      imputation = "Delta_i(s,t) = Y_it - Y_is on observed pairs S_i(s,t) = 1"),
    use = list(
      outcome = paste("m_b(X) = a_b + L(s,t) X'b is averaged over the target cohort (B_g, weights w a_g) and the odds-weighted",
        "controls (B_o) of each listed comparison; its estimation error is carried by the ATT IF"),
      treatment = paste("the capped treatment-propensity odds weight the controls' Z_c (C) and fitted m_b (B_o) of each listed comparison;",
        "pair identity never enters the fit, so comparisons with the same target and risk set reuse one fit"),
      imputation = paste("j_c is the predictor in Z_c = j_c + S_e / max(q_c, eps) (Delta_e - j_c) on every F row of",
        "cell c. The table lists direct comparison uses. With pooled outcome slopes, a control-cell prediction",
        "also enters the common slope.")),
    sample = list(
      outcome = "every F unit of each control cohort k of each control block (no observation requirement); rows are unit x control block",
      treatment = "every F unit of the target and control cohorts of the risk set",
      imputation = "units of the cell's cohort with both endpoints observed (S_i(s,t) = 1); rows are unit x cell"),
    weights = list(
      outcome = list(text = paste0(weights_text, " times a_k"), estimated = TRUE,
        note = "w_i a_k: the fixed weight times the estimated repeat rate of the unit's control cohort (a_k is in the SE)."),
      treatment = list(text = paste0(weights_text, " times a_{G_i}"), estimated = TRUE,
        note = "w_i a_{G_i}: the fixed weight times the estimated repeat rate of the unit's own cohort (a_k is in the SE); no q weight."),
      imputation = list(text = weights_text, estimated = FALSE,
        note = "fixed known weights on observed pairs only; no repeat rate and no observation-probability weight.")))
}

# Capture of the three DR+DR roles. `scope` is the selection-scope status
# (NULL when the ordinary static augmented route is admitted).
chained_nuisance_augmented_roles <- function(res, b, prep, weightsname, scope, selection_accounting = NULL) {
  status_all <- function(s) list(outcome = s, treatment = s, imputation = s)
  if (!is.null(scope)) return(status_all(scope))
  if (!chained_nuisance_aug_backend_ok(b))
    return(status_all(chained_nuisance_status("unavailable", "unrecognized_backend",
      "The augmented backend does not hold the DR+DR imputation, counterfactual OR and treatment fit records.")))
  validate_augmented_backend(b)
  st <- get("cp04_reader", b$reader_ref, inherits = FALSE)$state
  formula_text <- if (is.null(prep$cp04$outcome_formula)) NA_character_ else
    paste(deparse(prep$cp04$outcome_formula), collapse = " ")
  chained_nuisance_augmented_roles_run(b, st, b$control, prep$cluster_full, prep$clustervar,
    res$resource_accounting, formula_text, chained_nuisance_aug_weights_text(weightsname, prep),
    retained_entries = if (is.null(selection_accounting)) 0 else selection_accounting$output_entries,
    retained_work = if (is.null(selection_accounting)) 0 else selection_accounting$work_ops_used_bound)
}

# Shared by the capture and the low-level tests: admission, one bounded pass,
# tables. `group` borrows the original frame's cluster labels (NULL for units).
chained_nuisance_augmented_roles_run <- function(b, st, ctrl, group, cluster_var, ra, formula_text, weights_text,
                                                 retained_entries = 0, B = NULL, retained_work = 0) {
  N <- st$N
  status_all <- function(s) list(outcome = s, treatment = s, imputation = s)
  clustered <- !is.null(group)
  # Fund the fit descriptors and frame-sized grouping before constructing them.
  n_cell <- length(b$cells); n_or <- length(b$shared_or$blocks); n_trt <- length(b$treatments)
  n_fit_bound <- length(b$groups) + n_or + n_trt
  W_bound <- length(b$groups) * as.double(b$outcome_recipe$p) + n_cell +
    n_or * (as.double(b$outcome_recipe$p) + 1) + n_trt * as.double(b$treatment_recipe$p)
  text_entries <- function(x) sum(ceiling(nchar(x, type = "bytes", allowNA = TRUE) / 2), na.rm = TRUE)
  live <- ceiling(max(c(0, ra$simultaneous_entries, ra$initial_envelope_entries))) + retained_entries
  prior_work <- if (is.null(ra$downstream_work_ops)) 0 else ceiling(ra$downstream_work_ops)
  labels <- c(b$outcome_recipe$columns, b$treatment_recipe$columns, formula_text, weights_text, cluster_var)
  n_controls <- sum(vapply(b$links, function(l) length(l$control_cohort_indices), integer(1)))
  out_bound <- ceiling(3 * 65536 + 1024 * (n_fit_bound + n_cell + n_or + n_trt + 3 * length(b$links)) +
    256 * W_bound + 64 * n_controls + 16 * (n_fit_bound + n_cell + n_or + 4) * text_entries(labels) +
    24 * text_entries(vapply(b$links, function(l) as.character(l$link_id), "")))
  group_entries <- if (clustered) as.double(N) else 0
  meta_work <- 8 * out_bound + 64 * group_entries
  admitted <- tryCatch(cp04_budget(ctrl, c(live, out_bound), out_bound + 16 * group_entries,
    c(prior_work, retained_work, meta_work), c(retained_entries, out_bound),
    "nuisance augmented role tables"), cdid_ipw_owner_bound_error = identity)
  if (inherits(admitted, "error"))
    return(status_all(chained_nuisance_status("unavailable", "budget",
      paste("The optional augmented coefficient tables exceed the declared resource limits:", conditionMessage(admitted)))))

  fits <- chained_nuisance_aug_fits(b)
  flat <- unlist(fits, recursive = FALSE)
  n_fit <- length(flat)
  W <- sum(vapply(flat, function(f) length(f$estimate), integer(1)))
  code <- NULL; reason <- NULL; n_groups <- N
  if (clustered && (length(group) != N || anyNA(group))) {
    code <- "cluster_alignment"; reason <- "Cluster labels could not be aligned with the frame units by position."
    n_groups <- NA_integer_
  } else if (clustered) {
    group <- match(group, unique(group))
    n_groups <- max(group)
    if (n_groups < 2L) {
      code <- "insufficient_clusters"
      reason <- "At least two clusters are needed for a clustered coefficient standard error."
    }
  }
  B <- as.double(min(if (is.null(B)) ctrl$max_chunk_rows else B, ctrl$max_chunk_rows, N))
  work <- 0; scratch <- 0
  if (is.null(code)) {
    pre <- tryCatch(augmented_if_envelope(b, B, 1L, details = TRUE), cdid_ipw_owner_bound_error = identity,
      cdid_augmented_bound_error = identity)
    if (inherits(pre, "error")) {
      code <- "budget"
      reason <- paste("The coefficient IF replay exceeds the declared resource limits:", conditionMessage(pre))
    } else {
      n_chunks <- ceiling(N / B)
      acc_entries <- if (clustered) as.double(n_groups) * (W + n_fit) else 2 * W
      scratch <- pre$output + pre$scratch + 4 * B * W +
        (if (clustered) 2 * B * W + 3 * acc_entries else 0) + 8192
      work <- n_chunks * (pre$work + 32 * B * (W + n_fit) + 64 * (W + n_fit)) + 32 * acc_entries
      ok <- tryCatch(cp04_budget(ctrl, c(live, out_bound, acc_entries, group_entries), scratch,
        c(prior_work, retained_work, meta_work, work), c(retained_entries, out_bound),
        "nuisance augmented coefficient replay"), cdid_ipw_owner_bound_error = identity)
      if (inherits(ok, "error")) {
        code <- "budget"
        reason <- paste("The coefficient IF replay exceeds the declared workspace or work budget:", conditionMessage(ok))
      }
    }
  }
  stats <- list(var = vector("list", n_fit), clusters = rep(NA_integer_, n_fit))
  if (is.null(code)) {
    run <- tryCatch(chained_nuisance_aug_stream(b, flat, N, B, if (clustered) group else NULL, n_groups),
      cdid_ipw_owner_bound_error = identity, cdid_augmented_bound_error = identity)
    if (inherits(run, "error")) {
      code <- "budget"; reason <- paste("The coefficient IF replay was refused by a resource bound:", conditionMessage(run))
    } else if (!isTRUE(run$ok)) {
      code <- "if_layout"; reason <- run$reason
    } else {
      stats$var <- run$var; stats$clusters <- run$clusters
      if (clustered) for (j in seq_len(n_fit)) if (stats$clusters[j] < 2L) stats$var[j] <- list(NULL)
    }
  }
  if (!is.null(code)) reason <- paste(reason, "Coefficients are reported; SEs are not.")
  txt <- chained_nuisance_aug_meta_text(weights_text)
  pos <- 0L; out <- list()
  for (role in names(fits)) {
    ix <- pos + seq_along(fits[[role]]); pos <- pos + length(fits[[role]])
    cv <- chained_nuisance_aug_covariance(role, is.null(code), code, reason, clustered, cluster_var,
      if (is.null(code) || identical(code, "insufficient_clusters")) n_groups else NA_integer_, N)
    if (is.null(code) && clustered && any(stats$clusters[ix] < 2L)) {
      cv$status <- if (any(stats$clusters[ix] >= 2L)) "partial" else "unavailable"
      cv$code <- "insufficient_clusters"
      cv$reason <- "Each fit needs at least two clusters containing nonzero unit coefficient influence contributions."
    }
    meta <- c(txt, list(route = paste0(class(b)[1L], " (ordinary static augmented route, est_method = \"dr\")"),
      formula = formula_text, covariance = stats::setNames(list(cv), role),
      accounting = list(chunk_rows = B, links_replayed = 1L, replay_passes = if (is.null(code)) 1L else 0L,
        work_ops_bound = work + meta_work, scratch_entries_bound = scratch, output_entries_bound = out_bound,
        rule = paste("cp04_budget on the augmented backend's owner control with augmented_if_envelope(details = TRUE),",
          "composed with the pipeline's recorded live entries and downstream work"))))
    out[[role]] <- chained_nuisance_aug_role_tables(role, fits[[role]], b, st,
      list(var = stats$var[ix], clusters = stats$clusters[ix]), meta)
  }
  actual <- unname(cp04_payload(out)["entries"])
  if (actual > out_bound)
    return(status_all(chained_nuisance_status("unavailable", "budget",
      "The augmented coefficient tables exceeded their reserved output size.")))
  for (role in names(out)) out[[role]]$accounting$output_entries <- actual
  out
}

# Block view and printing of one role fit (equation directly above its table).
chained_nuisance_role_block_view <- function(component, comp, m) {
  bl <- comp$blocks[m, , drop = FALSE]
  co <- comp$coefficients[comp$coefficients$block == bl$block, , drop = FALSE]
  cm <- cbind(Estimate = co$estimate, `Std. Error` = co$std_error, `z value` = co$z_value, `Pr(>|z|)` = co$p_value)
  rownames(cm) <- co$term
  covariance <- comp$covariance
  if (!is.na(bl$n_contributing_clusters)) {
    covariance$n_contributing_clusters <- bl$n_contributing_clusters
    if (bl$n_contributing_clusters >= 2L && identical(covariance$status, "partial")) {
      covariance$status <- "available"; covariance$code <- "ok"; covariance$reason <- NULL
    } else if (bl$n_contributing_clusters < 2L) {
      covariance$status <- "unavailable"; covariance$code <- "insufficient_clusters"
      covariance$reason <- "This fit has fewer than two clusters containing nonzero unit coefficient influence contributions."
    }
  }
  dname <- switch(component, imputation = "cells", outcome = "control_blocks", treatment = "models")
  dt <- comp[[dname]]
  structure(list(component = component, block = bl$block, label = bl$label, pooling = comp$pooling,
    role = comp$role, use = comp$use, response = comp$response, equation = bl$equation,
    sample = if (is.null(bl$population)) {
      paste0(comp$sample, ": ", format(bl$n_fitting_rows), " fitting rows, cohorts {", bl$cohorts, "}")
    } else {
      paste0(bl$population, ": ", format(bl$n_fitting_rows), " fitting rows")
    },
    weights = comp$weights$text, weight_sum = bl$weight_sum,
    basis = list(kept = bl$kept_slopes, aliased = bl$aliased_slopes),
    coefficients = cm, table = co, covariance = covariance,
    detail = dt[dt$block == bl$block, , drop = FALSE],
    basis_transform = comp$basis_transform,
    comparisons = comp$comparisons[if ("block" %in% names(comp$comparisons)) {
      comp$comparisons$block == bl$block
    } else {
      vapply(strsplit(comp$comparisons$fits, ",", fixed = TRUE), function(x) bl$block %in% x, logical(1))
    }, , drop = FALSE]),
    class = "cdid_nuisance_block")
}

chained_nuisance_print_role_block <- function(x, ...) {
  cat(chained_nuisance_aug_title[[x$component]], " ", x$block, ": ", x$label, "  (", x$pooling, ")\n", sep = "")
  cat("  response: ", x$response, "\n", sep = "")
  cat("  sample: ", x$sample, "\n", sep = "")
  cat("  weights: ", x$weights, " (sum ", format(x$weight_sum, digits = 6), ")\n", sep = "")
  if (x$component != "treatment")
    cat("  basis: ", if (nzchar(x$basis$kept)) x$basis$kept else "none",
        if (nzchar(x$basis$aliased)) paste0("; aliased (not identified, dropped): ", x$basis$aliased) else "", "\n", sep = "")
  cat("  use: ", x$use, "\n", sep = "")
  cat("  direct comparison uses: ", nrow(x$comparisons), "\n", sep = "")
  cat("\n", x$equation, "\n", sep = "")
  if (!is.null(x$basis_transform) && nrow(x$basis_transform)) print(x$basis_transform, row.names = FALSE)
  if (nrow(x$coefficients))
    stats::printCoefmat(x$coefficients, has.Pvalue = TRUE, P.values = TRUE, na.print = "NA", ...)
  else cat("No coefficients estimated.\n")
  cat(chained_nuisance_se_line(x$covariance), "\n", sep = "")
  cat(sub("[.]$", "", x$covariance$accounts_for), ".\n", sep = "")
  invisible(x)
}

chained_nuisance_print_role_component <- function(x, ...) {
  cat(chained_nuisance_aug_title[[x$component]], ": ", x$role, "\n", sep = "")
  cat("  ", x$template, "\n", sep = "")
  cat("  ", x$covariance$accounts_for, ".\n\n", sep = "")
  for (m in seq_len(nrow(x$blocks))) {
    chained_nuisance_print_role_block(chained_nuisance_role_block_view(x$component, x, m), ...)
    cat("\n")
  }
  invisible(x)
}
