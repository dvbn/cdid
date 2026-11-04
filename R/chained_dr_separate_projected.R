# Ordinary baseline, non-reference separate DR only. The caller has already
# validated and admitted its original row/support/rate-score construction.
# Detailed, conditional, reference and complete-support replay stay unchanged.

# Contract the SAME frozen local coefficient IF with requested scalar model
# loadings. Keep score centering and nonzero rate forcing on every F row.
# Loadings/adjoints are p x nt; rate forcing uses one support vector per model,
# never a support-by-model matrix outside the existing conservative envelope.
dr_separate_baseline_adjoint <- function(backend, b, f, tids, hits, up, sids, w) {
  B <- length(w); nt <- length(tids); p <- f$p
  value <- matrix(0, B, nt)
  if (!p) return(value)
  rec <- backend$regressions[[b]]
  loads <- matrix(0, p, nt)
  for (j in seq_along(tids)) loads[, j] <- backend$treatments[[tids[j]]]$d_bar_diff
  adj <- if (backend$outcome_recipe$intercept)
    or_pooled_wls_gram_solve(f$fit, loads, normalizer = backend$outcome_groups$N) else
    or_wls_gram_solve(f$fit, loads, normalizer = backend$outcome_groups$N)
  controls <- backend$regression_specs[[b]]$support_ids
  pos <- match(controls, sids)
  for (k in seq_along(controls)) {
    sid <- controls[k]; h <- hits[[pos[k]]]; ix <- h$obs
    if (!length(ix)) next
    d <- f$L * h$z[, f$keep, drop = FALSE]
    residual <- h$dy - f$alpha - as.vector(d %*% f$beta)
    weighted_residual <- w[ix] * backend$supports[[sid]]$v * residual
    centered <- if (backend$outcome_recipe$intercept) sweep(d, 2L, f$mean_d) else d
    value[ix, ] <- value[ix, , drop = FALSE] + weighted_residual * (centered %*% adj)
  }
  value <- sweep(value, 2L, as.vector(rec$score_mean_beta %*% adj))
  us <- up[, pos, drop = FALSE]
  forcing <- t(rec$forcing_beta)
  for (j in seq_len(nt)) value[, j] <- value[, j] +
    as.vector(us %*% (forcing %*% adj[, j, drop = FALSE]))
  -value
}

dr_separate_baseline_projected_if <- function(backend, reader, rows, links, scale, pre, w, hits, up) {
  st <- reader$state; B <- length(rows); N <- st$N
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  requested <- backend$links[links]
  link_models <- vapply(requested, `[[`, integer(1), "model")
  by_model <- split(seq_along(links), factor(link_models, levels = pre$tids))
  names(by_model) <- as.character(pre$tids)
  treatment_blocks <- vapply(backend$treatment_specs[pre$tids], `[[`, integer(1), "regression")
  models_by_block <- split(pre$tids, factor(treatment_blocks, levels = pre$gids))
  link_names <- vapply(requested, `[[`, character(1), "link_id")
  phi <- matrix(0, B, length(links), dimnames = list(NULL, link_names))
  for (kb in seq_along(pre$gids)) {
    b <- pre$gids[kb]; tids <- models_by_block[[kb]]
    of <- or_outcome_block_fit(backend$outcome_groups, b)
    projected <- dr_separate_baseline_adjoint(backend, b, of, tids, hits, up, pre$sids, w)
    for (j_model in seq_along(tids)) {
      tid <- tids[j_model]; spec <- backend$treatment_specs[[tid]]; f <- backend$treatments[[tid]]
      us <- up[, match(spec$support_ids, pre$sids), drop = FALSE]
      pt <- backend$treatment_recipe$p
      score <- matrix(0, B, pt); direct <- matrix(0, B, 2L)
      for (k in seq_along(spec$support_ids)) {
        sid <- spec$support_ids[k]; h <- hits[[match(sid, pre$sids)]]; rr <- h$obs
        if (!length(rr)) next
        X <- frame$Z_a[st$a_index[rows[rr]], backend$treatment_recipe$index, drop = FALSE]
        dy <- h$dy; d <- of$L * h$z[, of$keep, drop = FALSE]
        residual <- dy - of$alpha - as.vector(d %*% of$beta)
        fw <- w[rr] * backend$supports[[sid]]$v
        lv <- ipw_logit_values(X, f$coefficients); D <- as.double(k == 1L)
        wp <- fw * D; wm <- fw * (1 - D) * lv$odds
        score[rr, ] <- fw * (D - lv$prob) * X
        direct[rr, 1L] <- wp * (residual - f$means["T_dr"]) / f$means["dp"]
        direct[rr, 2L] <- wm * (residual - f$means["M_dr"]) / f$means["dm"]
      }
      gamma <- t(cp04_solve(f$H_chol, t(score + us %*% t(f$cross))))
      pd <- direct[, 1L] - direct[, 2L]
      po <- projected[, j_model]
      pp <- -as.vector(gamma %*% f$partials$beta[2L, ])
      ps <- as.vector(us %*% (f$partials$upstream[1L, ] - f$partials$upstream[2L, ]))
      value <- pd + po + pp + ps
      for (j in by_model[[as.character(tid)]]) phi[, j] <- value
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
