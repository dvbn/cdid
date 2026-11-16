# Private baseline-only separate outcome integration. Catalog identity, owners,
# rates and reference recipes remain those of chained_or_baseline.R. No public
# dispatch or other observation/treatment route selects this implementation.

or_baseline_outcome_mode <- function(mode) {
  if (!is.character(mode) || is.object(mode) || !is.null(dim(mode)) ||
      length(mode) != 1L || is.na(mode) || !mode %in% c("pooled", "separate"))
    or_baseline_abort("Supply private outcome_pooling = pooled or separate.")
  mode
}

or_baseline_separate <- function(backend) identical(backend$outcome_pooling, "separate")

# Conservative generic envelope: reserve all G compact fits, alias replays and
# metadata before creating group layouts. No factor sharing is claimed here.
or_baseline_separate_fit_charges <- function(S, G, L, p, N, n_A, B, metadata, header, name_bytes,
                                              pools, links, support_rows) {
  P2 <- p + 2; incidence <- domains <- control_rows <- alias_rows <- 0; peak <- 0
  link_reg <- vapply(links, `[[`, integer(1L), "regression")
  link_sid <- vapply(links, `[[`, integer(1L), "support_id")
  for (b in seq_along(pools)) {
    ctrl <- pools[[b]]$support_ids
    target <- link_sid[link_reg == b]
    ids <- union(ctrl, target)
    incidence <- incidence + length(ctrl); domains <- domains + length(ids); peak <- max(peak, length(ids))
    control_rows <- control_rows + sum(support_rows[ctrl]); alias_rows <- alias_rows + sum(support_rows[ids])
  }
  # Each unique block owns one QR/summary/basis and its compact forcing. The
  # support summaries and one current group's aliases coexist during fitting.
  models <- 8192 + S * (1024 + 32 * P2^2) + G * (8192 + 128 * P2^2 + 64 * name_bytes) +
    128 * (incidence + domains + L) * P2 + L * (512 + 16 * p)
  scratch <- metadata + header + 64 * (S + peak + 2) * P2^2 +
    64 * B * (P2 + 8) + 128 * P2^2 + 32 * (domains + L + G)
  work <- 256 * (sum(support_rows) + control_rows + alias_rows) * P2^2 *
    (1 + ceiling(log2(n_A + 1))) + 256 * G * P2^3 +
    128 * (incidence + domains + L) * P2^2 + 128 * (metadata + header)
  list(models = ipw_owner_safe_sum(models), scratch = ipw_owner_safe_sum(scratch), work = ipw_owner_safe_sum(work))
}

or_baseline_separate_assemble <- function(recipe, supports, pools, links, rank_tol, N, visit_rows) {
  blocks <- lapply(seq_along(pools), function(b) list(name = paste0("b", b),
    support_ids = pools[[b]]$support_ids, L = pools[[b]]$transitions))
  fitted <- cdid_fast_or_baseline_outcome_fit(recipe, supports, blocks, links, rank_tol, N, "separate", "baseline")
  if (is.null(fitted)) {
    fitted <- or_outcome_fit_blocks(recipe, lapply(supports, function(u)
        list(summary = u$design$summary, scale = u$v, colmax = u$design$colmax)),
      blocks, lapply(links, function(l) list(regression = l$regression,
        support_id = l$support_id)), rank_tol, N, visit_rows, "separate", "baseline")
  }
  use_fast <- cdid_are_default_fns(c("or_outcome_block_fit", "or_pooled_factor_moments",
    "or_pooled_link_moments", "ipw_clip_floor", "cp04_finite"))
  icpt <- isTRUE(recipe$intercept)
  p_rec <- recipe$p
  block_fits <- if (use_fast) {
    lapply(seq_along(pools), function(b) {
      g <- fitted$groups[[b]]
      rec <- g$regressions[[1L]]
      c(list(block_id = b, group_id = b, local_block_id = 1L), rec,
        list(p = length(g$keep), coefficient_dimension = length(g$keep) + icpt,
             keep = g$keep, beta = g$beta, basis = g$basis, fit = g$fit))
    })
  } else {
    lapply(seq_along(pools), function(b) or_outcome_block_fit(fitted, b))
  }
  regressions <- lapply(seq_along(pools), function(b) {
    f <- block_fits[[b]]; sids <- pools[[b]]$support_ids
    rec <- fitted$groups[[b]]$regressions[[1L]]
    fa <- fe <- numeric(length(sids)); fb <- matrix(0, f$p, length(sids))
    # Explicit H_N(score) centering, on the frozen local alias basis.
    mean_alpha <- 0; mean_beta <- numeric(f$p)
    xi <- if (icpt) c(1L, 1L + f$keep) else f$keep
    neg_ab <- c(if (icpt) -f$alpha, -f$beta)
    y_col <- p_rec + 1L + icpt
    for (k in seq_along(sids)) {
      u <- supports[[sids[k]]]
      if (use_fast) {
        F_mat <- if (icpt) u$design$summary$raw$factor else u$design$summary$factor
        if (nrow(F_mat) > 0L) {
          F_x <- F_mat[, xi, drop = FALSE]
          rho <- F_mat[, y_col] + as.vector(F_x %*% neg_ab)
          g_vec <- as.vector(crossprod(F_x, rho))
          e2_val <- sum(rho^2)
        } else {
          g_vec <- numeric(length(xi)); e2_val <- 0
        }
        ga <- if (icpt) g_vec[1L] else 0
        gb <- if (icpt) g_vec[-1L] - f$mean_d * ga else g_vec
      } else {
        m <- or_pooled_factor_moments(recipe, u$design$summary, f$keep, f$alpha, f$beta)
        ga <- if (icpt) m$g[1L] else 0
        gb <- if (icpt) m$g[-1L] - f$mean_d * ga else m$g
        e2_val <- m$e2
      }
      mean_alpha <- mean_alpha + u$v * ga / N
      mean_beta <- mean_beta + u$v * gb / N
      scale <- if (u$rate > ipw_clip_eps) -1 / (N * u$rate^2) else 0
      fa[k] <- scale * ga; fb[, k] <- scale * gb; fe[k] <- scale * e2_val
    }
    if (!use_fast || !is.finite(sum(fa) + sum(fb) + sum(fe) + mean_alpha + sum(mean_beta))) {
      cp04_finite(c(fa, fb, fe, mean_alpha, mean_beta), "separate baseline forcing")
    }
    c(rec, list(forcing_alpha = fa, forcing_beta = fb, forcing_energy = fe,
      score_mean_alpha = mean_alpha, score_mean_beta = mean_beta))
  })
  for (i in seq_along(links)) {
    l <- links[[i]]; u <- supports[[l$support_id]]
    f <- block_fits[[l$regression]]
    if (use_fast) {
      mass <- u$pair_mass; sumZ_k <- u$design$sumZ[f$keep]; sumY <- u$design$sumY
      fitted_val <- (if (icpt) f$alpha * mass else 0) + f$L * sum(sumZ_k * f$beta)
      mu_val <- (sumY - fitted_val) / mass
      c_val <- f$L * sumZ_k / mass
      res_val <- sumY - fitted_val - mu_val * mass
      ra <- u$rate > ipw_clip_eps
      rc <- if (ra) u$rate else ipw_clip_eps
    } else {
      m <- or_pooled_link_moments(recipe, f$keep, u$pair_mass, u$design$sumZ,
        u$design$sumY, f$L, f$alpha, f$beta)
      mu_val <- m$mu; c_val <- m$c; res_val <- m$residual
      rc <- ipw_clip_floor(u$rate); ra <- u$rate > ipw_clip_eps
    }
    d <- u$pair_mass / (N * rc)
    rl <- if (ra) -res_val / (N * rc^2 * d) else 0
    links[[i]] <- c(l, list(mu = mu_val, c = c_val, d = d, numerator = mu_val * d,
      rate_loading = rl,
      log_weight_derivative = if (ra) -1 / rc else 0))
    if (!use_fast || !is.finite(mu_val + sum(c_val) + d + rl)) {
      cp04_finite(c(mu_val, c_val, d, rl), "separate baseline link moments")
    }
  }
  list(supports = lapply(supports, function(u) list(rate = u$rate, v = u$v, massA = u$massA,
      pair_mass = u$pair_mass, repeat_mass = u$repeat_mass, n_pair = u$n_pair,
      n_repeat = u$n_repeat, sumZ = u$design$sumZ, sumY = u$design$sumY,
      colmax = u$design$colmax, L = u$design$L)),
    regressions = regressions, links = links, outcome_groups = fitted)
}

# Preflight with global dimension bounds BEFORE allocating dependency sets.
# Details include every unique block, so the objective is the full objective;
# ordinary chunks need only requested blocks and their control/target supports.
or_baseline_separate_if_budget <- function(backend, B, links, details) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  G <- length(backend$regression_specs); S <- length(backend$supports)
  p <- backend$outcome_recipe$p; J <- length(st$cohorts); nl <- length(links)
  incidence <- 0
  for (r in backend$regression_specs) incidence <- incidence + length(r$support_ids)
  maps <- 4096 + 32 * (S + G + nl + incidence)
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    maps, 64 * maps, 0, "separate baseline OR dependency metadata")
  link_blocks <- if (!is.null(backend$separate_complete)) backend$separate_complete$link_block[links] else
    vapply(backend$links[links], `[[`, integer(1), "regression")
  gids <- if (details) seq_len(G) else sort(unique(link_blocks))
  sids <- sort(unique(c(vapply(backend$links[links], `[[`, integer(1), "support_id"),
    unlist(lapply(backend$regression_specs[gids], `[[`, "support_ids")))))
  ps <- vapply(backend$outcome_groups$groups[gids], function(g) length(g$keep), integer(1))
  ns <- length(sids); ng <- length(gids); mp <- max(c(0, ps))
  name_bytes <- sum(nchar(backend$outcome_recipe$columns, type = "bytes"))
  link_bytes <- sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  metadata <- 8192 + 16 * link_bytes + B * frame$id_byte_metadata$max +
    if (details) 1536 * (ng + nl + 1) + 32 * name_bytes * (ng + nl + 1) + 32 * (ns + ng) else 128 * nl
  width <- nl + 8 + if (details) 2 * sum(ps + 1) + ns + nl * (p + 3) + 2 * J + 4 else 0
  output <- B * width + metadata
  # Ordinary fixed-operator replay solves p_b-by-n_requested_links_b
  # loadings without a row dimension. Reserve all solve copies before allocation.
  adjoint_peak <- if (!details && !is.null(backend$separate_complete))
    max(c(0, ps * tabulate(match(link_blocks, gids), nbins = length(gids)))) else 0
  scratch <- maps + output + 32 * B * (ns * (p + 5) + mp * (mp + ns + 6) + nl + 16) +
    128 * (mp + 2)^2 + 32 * adjoint_peak
  work <- 128 * maps + 128 * B * ((ns + 1) * (p + 2)^2 * (1 + ceiling(log2(st$n_A + 1))) +
    sum((ps + 2)^2 * (ns + nl + 2)) + J + nl * (p + 2))
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    scratch, work, output, "separate baseline OR IF chunk")
  list(gids = gids, beta_gids = gids, sids = sids, output = output, scratch = scratch, work = work,
    adjoint_entries_peak = adjoint_peak)
}

# One bounded original-unit chunk, streamed by unique block. Only details retain
# the requested chunk's coefficient vectors. hits carry source-recipe slopes;
# each block selects its own frozen kept columns before using d = L*z.
or_baseline_separate_if_core <- function(model, N, w, hits, up, pre, links, details, block_fits = NULL) {
  B <- length(w); icpt <- model$outcome_recipe$intercept
  phi <- direct <- matrix(0, B, length(links))
  blocks <- if (details) vector("list", length(pre$gids)) else NULL
  if (details) {
    den <- num <- matrix(0, B, length(links)); design <- vector("list", length(links))
    objective <- numeric(B)
  }
  col <- function(s) match(s, pre$sids)
  # Requested positions only; do not rescan all links for every block.
  link_blocks <- vapply(model$links[links], `[[`, integer(1), "regression")
  by_block <- split(seq_along(links), factor(link_blocks, levels = pre$gids))
  for (k in seq_along(pre$gids)) {
    b <- pre$gids[k]
    f <- if (is.null(block_fits)) or_outcome_block_fit(model$outcome_groups, b) else block_fits[[k]]
    rec <- model$regressions[[b]]; sids <- model$regression_specs[[b]]$support_ids
    resid <- function(s) {
      h <- hits[[col(s)]]; d <- f$L * h$z[, f$keep, drop = FALSE]
      list(ix = h$obs, d = d, q = (if (is.null(h$w)) w[h$obs] else h$w) * model$supports[[s]]$v,
        r = h$dy - f$alpha - as.vector(d %*% f$beta))
    }
    ha <- energy <- numeric(B); hb <- matrix(0, B, f$p)
    for (s in sids) {
      e <- resid(s); ix <- e$ix
      if (length(ix)) {
        ha[ix] <- ha[ix] + e$q * e$r
        if (f$p) hb[ix, ] <- hb[ix, , drop = FALSE] + (e$q * e$r) *
          (if (icpt) sweep(e$d, 2L, f$mean_d) else e$d)
        if (details) energy[ix] <- energy[ix] + e$q * e$r^2
      }
    }
    us <- up[, col(sids), drop = FALSE]
    xb <- hb
    if (f$p) {
      rhs <- sweep(hb, 2L, rec$score_mean_beta) + us %*% t(rec$forcing_beta)
      # Reuse this block's fetched frozen fit; the operator is unchanged.
      xb <- t(if (icpt) or_pooled_wls_gram_solve(f$fit, t(rhs), normalizer = model$outcome_groups$N) else
        or_wls_gram_solve(f$fit, t(rhs), normalizer = model$outcome_groups$N))
    }
    xa <- if (icpt) (ha - rec$score_mean_alpha + as.vector(us %*% rec$forcing_alpha)) /
      rec$mass - as.vector(xb %*% rec$mean_d) else numeric(B)
    if (details) {
      colnames(xb) <- colnames(hb) <- f$basis$kept_columns
      blocks[[k]] <- list(block_id = b, group_id = f$group_id, local_block_id = f$local_block_id,
        keep = f$keep, slope_columns = f$basis$kept_columns,
        coefficient_columns = c(if (icpt) "(Intercept)", f$basis$kept_columns), beta = xb, beta_score = hb,
        alpha = if (icpt) xa else NULL, alpha_score = if (icpt) ha else NULL)
      objective <- objective + energy - model$outcome_groups$groups[[f$group_id]]$objective +
        as.vector(us %*% rec$forcing_energy)
    }
    for (j in by_block[[k]]) {
      l <- model$links[[links[j]]]; e <- resid(l$support_id); ix <- e$ix
      if (length(ix)) direct[ix, j] <- e$q * (e$r - l$mu) / l$d
      ul <- up[, col(l$support_id), drop = FALSE]
      # Standalone OR retains alpha and the treated mean loading on beta.
      phi[, j] <- direct[, j] - xa - as.vector(xb %*% l$c) + as.vector(ul %*% l$rate_loading)
      if (details) {
        tw <- tr <- numeric(B); tw[ix] <- e$q; tr[ix] <- e$q * e$r
        hd <- l$c * l$d; lo <- or_pooled_link_loadings(l, f$p)
        den[, j] <- tw - l$d + as.vector(ul %*% lo$weight)
        num[, j] <- tr - l$numerator + as.vector(ul %*% lo$numerator) -
          (if (icpt) l$d * xa else 0) - as.vector(xb %*% hd)
        dm <- matrix(0, B, f$p)
        if (length(ix) && f$p) dm[ix, ] <- e$q * e$d
        design[[j]] <- dm - matrix(hd, B, f$p, byrow = TRUE) + ul %*% t(lo$design)
        colnames(design[[j]]) <- f$basis$kept_columns
      }
    }
  }
  out <- list(phi = phi, direct = direct)
  if (details) {
    names(blocks) <- as.character(pre$gids)
    out <- c(out, list(blocks = blocks, denominator = den, numerator = num,
      design = design, objective = objective))
  }
  out
}
