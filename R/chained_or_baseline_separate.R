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
  for (b in seq_along(pools)) {
    ctrl <- pools[[b]]$support_ids
    target <- vapply(links[vapply(links, function(l) l$regression == b, logical(1))], `[[`, integer(1), "support_id")
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
  fitted <- or_outcome_fit_blocks(recipe, lapply(supports, function(u)
      list(summary = u$design$summary, scale = u$v, colmax = u$design$colmax)),
    blocks, lapply(links, function(l) list(regression = l$regression,
      support_id = l$support_id)), rank_tol, N, visit_rows, "separate", "baseline")
  regressions <- lapply(seq_along(pools), function(b) {
    f <- or_outcome_block_fit(fitted, b); sids <- pools[[b]]$support_ids
    rec <- fitted$groups[[b]]$regressions[[1L]]
    fa <- fe <- numeric(length(sids)); fb <- matrix(0, f$p, length(sids))
    # Explicit H_N(score) centering, on the frozen local alias basis.
    mean_alpha <- 0; mean_beta <- numeric(f$p)
    for (k in seq_along(sids)) {
      u <- supports[[sids[k]]]
      m <- or_pooled_factor_moments(recipe, u$design$summary, f$keep, f$alpha, f$beta)
      ga <- if (recipe$intercept) m$g[1L] else 0
      gb <- if (recipe$intercept) m$g[-1L] - f$mean_d * ga else m$g
      mean_alpha <- mean_alpha + u$v * ga / N
      mean_beta <- mean_beta + u$v * gb / N
      scale <- if (u$rate > ipw_clip_eps) -1 / (N * u$rate^2) else 0
      fa[k] <- scale * ga; fb[, k] <- scale * gb; fe[k] <- scale * m$e2
    }
    cp04_finite(c(fa, fb, fe, mean_alpha, mean_beta), "separate baseline forcing")
    c(rec, list(forcing_alpha = fa, forcing_beta = fb, forcing_energy = fe,
      score_mean_alpha = mean_alpha, score_mean_beta = mean_beta))
  })
  for (i in seq_along(links)) {
    l <- links[[i]]; u <- supports[[l$support_id]]
    f <- or_outcome_block_fit(fitted, l$regression)
    m <- or_pooled_link_moments(recipe, f$keep, u$pair_mass, u$design$sumZ,
      u$design$sumY, f$L, f$alpha, f$beta)
    rc <- ipw_clip_floor(u$rate); ra <- u$rate > ipw_clip_eps
    d <- u$pair_mass / (N * rc)
    links[[i]] <- c(l, list(mu = m$mu, c = m$c, d = d, numerator = m$mu * d,
      rate_loading = if (ra) -m$residual / (N * rc^2 * d) else 0,
      log_weight_derivative = if (ra) -1 / rc else 0))
    cp04_finite(c(m$mu, m$c, d, links[[i]]$rate_loading), "separate baseline link moments")
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
  gids <- if (details) seq_len(G) else sort(unique(vapply(backend$links[links], `[[`, integer(1), "regression")))
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
  scratch <- maps + output + 32 * B * (ns * (p + 5) + mp * (mp + ns + 6) + nl + 16) +
    128 * (mp + 2)^2
  work <- 128 * maps + 128 * B * ((ns + 1) * (p + 2)^2 * (1 + ceiling(log2(st$n_A + 1))) +
    sum((ps + 2)^2 * (ns + nl + 2)) + J + nl * (p + 2))
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    scratch, work, output, "separate baseline OR IF chunk")
  list(gids = gids, beta_gids = gids, sids = sids, output = output, scratch = scratch, work = work)
}

# One bounded original-unit chunk, streamed by unique block. Only details retain
# the requested chunk's coefficient vectors. hits carry source-recipe slopes;
# each block selects its own frozen kept columns before using d = L*z.
or_baseline_separate_if_core <- function(model, N, w, hits, up, pre, links, details) {
  B <- length(w); icpt <- model$outcome_recipe$intercept
  phi <- direct <- matrix(0, B, length(links))
  blocks <- if (details) vector("list", length(pre$gids)) else NULL
  if (details) {
    den <- num <- matrix(0, B, length(links)); design <- vector("list", length(links))
    objective <- numeric(B)
  }
  col <- function(s) match(s, pre$sids)
  for (k in seq_along(pre$gids)) {
    b <- pre$gids[k]; f <- or_outcome_block_fit(model$outcome_groups, b)
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
      xb <- t(or_outcome_block_gram_solve(model$outcome_groups, b, t(rhs)))
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
    for (j in which(vapply(model$links[links], function(l) l$regression == b, logical(1)))) {
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
