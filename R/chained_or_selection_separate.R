# Private separate counterfactual OR for accepted conditional-observation routes.
# Selection fits, cells, transport gates and upstream identities remain global.
# Only the existing unique comparison blocks receive distinct outcome slopes.

or_selection_outcome_mode <- function(mode) {
  if (!is.character(mode) || is.object(mode) || !is.null(dim(mode)) ||
      length(mode) != 1L || is.na(mode) || !mode %in% c("pooled", "separate"))
    or_corollary1_abort("Supply private outcome_pooling = pooled or separate.")
  mode
}

or_selection_separate <- function(backend) identical(backend$outcome_pooling, "separate")

# Reserve all compact group operators, forcing and possible alias rescans before
# creating any group. This conservative charge makes no factor-sharing claim.
or_selection_separate_fit_charges <- function(S, G, L, p, d, D, n_A, B, max_fit_n, visits, maxit,
                                               metadata, header, name_bytes, pools, links, channels,
                                               support_rows, selection_fitted = TRUE) {
  x <- or_baseline_separate_fit_charges(S, G, L, p, sum(support_rows), n_A, B,
    metadata, header, name_bytes, pools, links, support_rows)
  forcing <- 0
  for (b in seq_along(pools)) forcing <- forcing +
    length(unique(unlist(channels[pools[[b]]$support_ids])))
  # Global coordinates are used only by the current forcing accumulator; the
  # retained records store their actual local channel subsets.
  x$models <- x$models + 32 * (p + 3) * forcing + 32 * (p + 2) * sum(lengths(channels))
  x$scratch <- x$scratch + 64 * S * (p + 2)^2 * (d + 1) + 64 * (p + 3) * (D + 1)
  x$work <- x$work + 128 * visits * (p + 2)^2 * (d + 1) +
    128 * (p + 3) * forcing + 128 * L * (p + 2)^2 * (d + 1)
  if (!selection_fitted) {
    x$models <- x$models + S * (d * d + d + 2048)
    x$scratch <- max(x$scratch, metadata + header + 128 * max_fit_n * (d + 16))
    x$work <- x$work + 256 * visits * (maxit + 2) * (d + 1)^2 + 256 * S * (d + 1)^3
  }
  lapply(x, ipw_owner_safe_sum)
}

or_selection_separate_assemble <- function(recipe, supports, pools, links, channels,
                                           D, rank_tol, N, visit_rows,
                                           forcing_fun = or_selection_forcing,
                                           link_fun = or_selection_link) {
  blocks <- lapply(seq_along(pools), function(b) list(name = paste0("b", b),
    support_ids = pools[[b]]$support_ids, L = pools[[b]]$transitions))
  fitted <- or_outcome_fit_blocks(recipe, lapply(supports, function(u)
      list(summary = u$design$summary, scale = 1, colmax = u$design$colmax)),
    blocks, lapply(links, function(l) list(regression = l$regression,
      support_id = l$support_id)), rank_tol, N, visit_rows, "separate", "corollary1")
  regressions <- lapply(seq_along(pools), function(b) {
    f <- fitted$groups[[b]]
    rec <- forcing_fun(recipe, f, supports, pools[b], channels, D, N)$regressions[[1L]]
    mean_alpha <- 0; mean_beta <- numeric(length(f$keep))
    for (s in pools[[b]]$support_ids) {
      g <- or_pooled_factor_moments(recipe, supports[[s]]$design$summary,
        f$keep, rec$alpha, f$beta)
      ga <- if (recipe$intercept) g$g[1L] else 0
      gb <- if (recipe$intercept) g$g[-1L] - rec$mean_d * ga else g$g
      mean_alpha <- mean_alpha + ga / N
      mean_beta <- mean_beta + gb / N
    }
    c(rec, list(score_mean_alpha = mean_alpha, score_mean_beta = mean_beta))
  })
  for (i in seq_along(links)) {
    l <- links[[i]]
    links[[i]] <- c(l, link_fun(recipe, fitted$groups[[l$regression]],
      supports[[l$support_id]], regressions[[l$regression]], N))
  }
  list(supports = lapply(seq_along(supports), function(s) {
    u <- supports[[s]]
    list(a = u$a, channels = channels[[s]], pair_mass = u$pair_mass, n_pair = u$n_pair,
      sumZ = u$design$sumZ, sumY = u$design$sumY, colmax = u$design$colmax,
      L = u$design$L, clipping = u$clipping)
  }), regressions = regressions, links = links, outcome_groups = fitted)
}

or_selection_separate_record <- function(obj, state) {
  obj$outcome_groups <- state$outcome_groups
  obj$law <- sub("shared-one-period-slope", "separate-one-period-slopes", obj$law, fixed = TRUE)
  groups <- state$outcome_groups$groups
  obj$accounting$n_regression_fits <- length(groups)
  obj$accounting$n_block_slopes <- sum(vapply(groups, function(g) length(g$keep), integer(1)))
  obj$accounting$n_aliased_slopes <- sum(vapply(groups,
    function(g) length(g$basis$aliased_columns), integer(1)))
  obj$accounting$n_support_scans <- length(state$supports) +
    sum(vapply(seq_along(groups), function(b) if (groups[[b]]$verification_scan)
      length(state$outcome_groups$group_map[[b]]$support_ids) else 0L, integer(1)))
  obj
}

# Reserve global dimension bounds before constructing dependency maps. Details
# visit every outcome block for the full objective; ordinary requests visit only
# their blocks. Shared upstream selection groups still bring every needed cell.
or_selection_separate_if_envelope <- function(backend, B, links, details) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  G <- length(backend$regression_specs); S <- length(backend$supports)
  p <- backend$outcome_recipe$p; D <- backend$upstream$D; d <- backend$upstream$d
  J <- length(st$cohorts); nl <- length(links); incidence <- 0
  for (r in backend$regression_specs) incidence <- incidence + length(r$support_ids)
  for (u in backend$supports) incidence <- incidence + length(u$channels)
  maps <- 4096 + 32 * (D + S + G + nl + incidence)
  cp04_budget(backend$control, c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
    maps, c(backend$accounting$work_ops_bound, 64 * maps), 0, "separate conditional OR dependency metadata")
  gids <- if (details) seq_len(G) else sort(unique(vapply(backend$links[links], `[[`, integer(1), "regression")))
  sids <- sort(unique(c(vapply(backend$links[links], `[[`, integer(1), "support_id"),
    unlist(lapply(backend$regression_specs[gids], `[[`, "support_ids")))))
  uids <- sort(unique(unlist(lapply(backend$supports[sids], `[[`, "channels"))))
  ns <- length(sids); nu <- length(uids); ng <- length(gids)
  ps <- vapply(backend$outcome_groups$groups[gids], function(g) length(g$keep), integer(1)); mp <- max(c(0, ps))
  name_bytes <- sum(nchar(backend$outcome_recipe$columns, type = "bytes"))
  link_bytes <- sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  metadata <- 8192 + 16 * link_bytes + B * frame$id_byte_metadata$max +
    if (details) 1536 * (ng + nl + 1) + 32 * name_bytes * (ng + nl + 1) + 32 * (ns + nu + ng) else 128 * nl
  width <- nl + 8 + if (details) 2 * sum(ps + 1) + nu + nl * (p + 3) + 2 * J + 4 else 0
  output <- B * width + metadata
  scratch <- maps + output + 32 * B * (ns * (p + 6) + nu + mp * (mp + nu + 6) + d + nl + 16) +
    128 * ((mp + 2)^2 + (d + 1)^2)
  work <- 128 * maps + 128 * B * ((ns + 1) * ((p + 2)^2 + (d + 1)^2) *
    (1 + ceiling(log2(st$n_A + 1))) + sum((ps + 2)^2 * (nu + nl + 2)) + J + nl * (p + nu + 2))
  if (isTRUE(backend$uses_selection_models)) {
    pre <- selection_if_budget(backend$selection_backend, B, unique(backend$q_of_s[sids]), details = details)
    scratch <- scratch + pre$output + pre$scratch; work <- work + pre$work
  }
  list(gids = gids, beta_gids = gids, sids = sids, uids = uids, output = output, scratch = scratch, work = work)
}

# Stream one block's operator over one original-unit chunk. All selection/rate
# IFs are borrowed once in global coordinates. h$v is the source-owned accepted
# a/(rho*pi*) factor, h$w is the recipe weight when reference dates are used.
or_selection_separate_if_core <- function(model, N, w, hits, up, pre, links, details) {
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
      list(ix = h$obs, d = d, q = (if (is.null(h$w)) w[h$obs] else h$w) * h$v,
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
    us <- up[, match(rec$channel_ids, pre$uids), drop = FALSE]
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
        coefficient_columns = c(if (icpt) "(Intercept)", f$basis$kept_columns), beta = xb,
        beta_score = hb, alpha = if (icpt) xa else NULL, alpha_score = if (icpt) ha else NULL)
      objective <- objective + energy - model$outcome_groups$groups[[f$group_id]]$objective +
        as.vector(us %*% rec$forcing_energy)
    }
    for (j in which(vapply(model$links[links], function(l) l$regression == b, logical(1)))) {
      l <- model$links[[links[j]]]; e <- resid(l$support_id); ix <- e$ix
      if (length(ix)) direct[ix, j] <- e$q * (e$r - l$mu) / l$d
      ul <- up[, match(model$supports[[l$support_id]]$channels, pre$uids), drop = FALSE]
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
