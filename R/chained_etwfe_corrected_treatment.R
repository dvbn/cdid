# Corrected treatment moments are normalized on the complete cohort frame.
# The observation correction enters the response, not the treatment root.
etwfe_corrected_treatment_layout <- function(catalog, control_group, method, entries) {
  K <- ncol(catalog$Y)
  cells <- catalog$cells
  groups <- catalog$groups
  nyt <- identical(control_group, "notyettreated")
  dr <- identical(method, "dr")
  treated <- groups[groups > 0L]
  usable <- function(g) cells$cohort == g | cells$cohort == 0L |
    (nyt & cells$t < cells$cohort)
  n_links <- sum(vapply(treated, function(g) sum(usable(g)), 0))
  if (catalog$entries + 22 * n_links > entries) etwfe_tm_workspace_error()
  catalog_rows <- lapply(treated, function(g) {
    target <- cells$cohort == g
    use <- which(usable(g))
    data.frame(g = rep.int(g, length(use)), cell = use,
               donor = ifelse(target[use], NA_integer_, cells$cohort[use]))
  })
  layout <- do.call(rbind, catalog_rows)
  layout$s <- cells$s[layout$cell]
  layout$t <- cells$t[layout$cell]
  layout <- layout[order(layout$g, layout$s, layout$t, !is.na(layout$donor), layout$donor), ]
  L <- nrow(layout)
  target <- is.na(layout$donor)
  links <- data.frame(
    g = layout$g, s = layout$s, t = layout$t, ref = layout$g - 1L,
    link_role = ifelse(target, "target", "control"),
    functional_key = ifelse(target, "target", "calendar"), donor = as.integer(layout$donor)
  )
  design <- calendar_geometry <- NULL
  if (L) {
    grid <- expand.grid(g = treated, time = seq_len(K))
    ledger <- grid[grid$time >= grid$g, , drop = FALSE]
    design <- build_cross_pair_chained_design(links, seq_len(K), c(treated, Inf),
      anticipation = 0L, pt_assumption = if (nyt) "all" else "post",
      pt_start = if (nyt) 1L else NULL, ledger = ledger)
    params <- design$parameters
    observed <- catalog$observed[cbind(grid$time, match(grid$g, groups))]
    pre <- grid$time < grid$g - 1L
    requested <- which(observed & (grid$time >= grid$g | (!nyt & pre)))
    components <- vapply(requested, function(i) {
      roles <- if (pre[i]) c("pre_target", "pre") else c("post_target", "trend_post")
      a <- cross_pair_find_param(params, grid$g[i], grid$time[i], roles[1L])
      b <- cross_pair_find_param(params, grid$g[i], grid$time[i], roles[2L])
      if (is.na(a) || is.na(b)) return(NA_integer_)
      ids <- params$component_id[c(a, b)]
      if (ids[1L] == ids[2L]) ids[1L] else NA_integer_
    }, integer(1))
    keep <- design$links$component_id %in% components[!is.na(components)]
    if (!all(keep)) {
      layout <- layout[keep, , drop = FALSE]
      links <- links[keep, , drop = FALSE]
      design <- NULL
    }
    donor <- which(!is.na(layout$donor))
    if (dr && length(donor)) {
      clean <- which(cells$cohort == 0L | (nyt & cells$t < cells$cohort))
      calendar_geometry <- etwfe_corrected_projection_geometry(catalog,
        matrix(seq_len(K), K, length(groups)), clean,
        cbind(layout$s[donor], layout$t[donor]), entries)
    }
  }
  list(layout = layout, links = links, design = design, calendar = calendar_geometry,
    used_cells = sort(unique(c(layout$cell, calendar_geometry$use))))
}

etwfe_corrected_treatment <- function(catalog, cohort, features, Z, obs, control_group,
                                     method, compute_scores, glm_control, chunk_rows, entries, report = NULL,
                                     support = etwfe_corrected_treatment_layout(catalog, control_group, method, entries)) {
  n <- length(cohort)
  K <- ncol(catalog$Y)
  q <- ncol(features)
  cells <- catalog$cells
  groups <- catalog$groups
  dr <- identical(method, "dr")
  layout <- support$layout
  links <- support$links
  L <- nrow(links)
  target <- is.na(layout$donor)
  diagnostics <- list(omitted_pairs = catalog$omitted_pairs)
  if (!L) return(list(links = links, mhat = numeric(), score_rows = NULL, diagnostics = diagnostics))

  grp <- match(layout$g, groups)
  members <- catalog$units
  mu <- t(matrix(vapply(members, function(rows) colMeans(features[rows, , drop = FALSE]), numeric(q)),
                 nrow = q))
  models <- unique(layout[!target, c("g", "donor")])
  model_ix <- match(paste(layout$g, layout$donor), paste(models$g, models$donor))
  P <- obs$n_parameters
  retained_obs <- if (is.null(obs$diagnostics$budget$retained_entries)) 0 else obs$diagnostics$budget$retained_entries
  local_columns <- if (is.null(obs$local_width)) P else obs$local_width
  fixed <- as.double(catalog$entries) + retained_obs + L * (22 + 3 * q + 2 * local_columns) +
    as.double(n) * (q + 1) + nrow(models) * (q^2 + q + 8)
  fit_scratch <- 8 * as.double(n) * (q + 2)
  per_unit <- 4 * as.double(L) + 3 * P + 6 * q * K + nrow(cells) +
    nrow(models) * (q + 4) + 2 * q + 10
  if (fixed + max(fit_scratch, per_unit) > entries) etwfe_tm_workspace_error()
  block_rows <- min(chunk_rows, floor((entries - fixed) / per_unit))
  blocks <- function(ids) split(ids, (seq_along(ids) - 1L) %/% block_rows)

  calendar <- NULL
  delta <- matrix(0, L, q)
  if (dr && any(!target)) {
    geometry <- support$calendar
    calendar <- etwfe_corrected_projection(
      catalog, features, geometry, obs, compute_scores, block_rows, entries,
      extra = fixed - catalog$entries - retained_obs
    )
    fixed <- fixed + calendar$live_entries
    if (fixed + per_unit > entries) etwfe_tm_workspace_error()
    block_rows <- min(block_rows, floor((entries - fixed) / per_unit))
    if (compute_scores) block_rows <- min(block_rows, calendar$block_rows)
    delta[!target, ] <- calendar$beta[layout$t[!target], , drop = FALSE] -
      calendar$beta[layout$s[!target], , drop = FALSE]
  }

  mhat <- denom <- residual_mean <- numeric(L)
  gamma <- imbalance <- matrix(0, L, q)
  upstream <- if (compute_scores) vector("list", L)
  odds_fits <- vector("list", nrow(models))
  for (m in seq_len(nrow(models))) {
    rows <- which(cohort %in% c(models$g[m], models$donor[m]))
    odds_fits[[m]] <- cp04_logit(
      etwfe_tm_design(Z, rows), as.double(cohort[rows] == models$g[m]),
      rep.int(1, length(rows)), glm_control, n,
      list(model = "etwfe_treatment_odds", g = models$g[m], donor = models$donor[m])
    )
  }

  # Read each cohort-pair once per block. The same reconstructed change can
  # enter its own target equation and several donor-transport equations.
  by_cell <- split(seq_len(L), factor(layout$cell, levels = seq_len(nrow(cells))))
  for (c in which(lengths(by_cell) > 0L)) {
    edges <- by_cell[[c]]
    terms <- vector("list", length(edges))
    for (rows in blocks(members[[cells$group[c]]])) {
      out <- obs$read(c, rows, gradient = compute_scores)
      w <- features[rows, , drop = FALSE]
      x <- etwfe_tm_design(Z, rows)
      for (k in seq_along(edges)) {
        j <- edges[k]
        if (target[j]) {
          tilt <- rep.int(1, length(rows))
          active <- numeric(length(rows))
          e <- out$value
        } else {
          lv <- ipw_logit_values(x, odds_fits[[model_ix[j]]]$coefficients)
          tilt <- lv$odds
          active <- lv$active
          e <- if (dr) out$value - drop(w %*% delta[j, ]) else out$value
        }
        part <- terms[[k]]
        if (is.null(part)) part <- list(mass = 0, total = 0, tangent = numeric(length(out$coordinates)),
                                       d_den = numeric(q), d_num = numeric(q), x_sum = numeric(q))
        part$mass <- part$mass + sum(tilt)
        part$total <- part$total + sum(tilt * e)
        if (compute_scores) {
          if (length(out$coordinates)) part$tangent <- part$tangent + drop(crossprod(out$gradient, tilt))
          if (!target[j]) {
            derivative <- tilt * active
            part$d_den <- part$d_den + drop(crossprod(x, derivative))
            part$d_num <- part$d_num + drop(crossprod(x, derivative * e))
            if (dr) part$x_sum <- part$x_sum + drop(crossprod(w, tilt))
          }
        }
        terms[[k]] <- part
      }
    }
    for (k in seq_along(edges)) {
      j <- edges[k]
      part <- terms[[k]]
      denom[j] <- part$mass
      residual_mean[j] <- part$total / part$mass
      mhat[j] <- residual_mean[j] + sum(mu[grp[j], ] * delta[j, ])
      if (compute_scores) {
        upstream[[j]] <- list(columns = out$coordinates, load = part$tangent / part$mass)
        gamma[j, ] <- (part$d_num - residual_mean[j] * part$d_den) / part$mass
        if (dr && !target[j]) imbalance[j, ] <- mu[grp[j], ] - part$x_sum / part$mass
      }
    }
  }
  if (any(!is.finite(mhat)) || any(!is.finite(denom) | denom <= 0)) {
    etwfe_tm_fail("invalid full-frame treatment moment or normalizer.", "cdid_error_etwfe_numeric")
  }
  diagnostics$normalizer <- denom
  diagnostics$workspace <- list(fixed = fixed, score_row = per_unit)
  solver_reserve <- list(fixed = fixed, row_scratch = max(0, per_unit - 2 * L),
                         max_batch_rows = block_rows)
  result <- list(links = links, mhat = mhat, score_rows = NULL,
    diagnostics = diagnostics, solver_reserve = solver_reserve, design = support$design)
  if (!is.null(report)) {
    result$nuisance <- list(treatment = etwfe_nuisance_treatment(odds_fits, models$g, models$donor,
      cohort, Z, links, report, entries - fixed))
    if (!is.null(calendar)) {
      k <- which(!target)
      comparisons <- etwfe_nuisance_comparisons(links$g[k], links$s[k], links$t[k],
        links$donor[k], rep("calendar", length(k)), report)
      result$nuisance$outcome <- etwfe_nuisance_projection(catalog, geometry, calendar, obs, Z,
        report, entries - fixed + catalog$entries + retained_obs + calendar$live_entries,
        "corrected_calendar", comparisons)
    }
  }
  if (!compute_scores) return(result)

  score_block <- function(rows) {
    R <- length(rows)
    out <- matrix(0, R, L)
    g_of <- cohort[rows]
    w <- features[rows, , drop = FALSE]
    X <- etwfe_tm_design(Z, rows)
    nuisance <- obs$score_rows(rows)
    psi <- if (!is.null(calendar)) calendar$score(rows, nuisance)
    used_loads <- vector("list", nrow(models))
    for (m in seq_len(nrow(models))) {
      selected <- which(g_of %in% c(models$g[m], models$donor[m]))
      fit <- odds_fits[[m]]
      x <- X[selected, , drop = FALSE]
      lv <- ipw_logit_values(x, fit$coefficients)
      g <- g_of[selected] == models$g[m]
      used_loads[[m]] <- list(
        selected = selected, odds = lv$odds[!g],
        IF = t(cp04_solve(fit$H_chol, t(x * (g - lv$prob))))
      )
    }
    for (c in which(lengths(by_cell) > 0L)) {
      loc <- which(g_of == cells$cohort[c])
      z <- if (length(loc)) obs$read(c, rows[loc])$value else numeric()
      for (j in by_cell[[c]]) {
        up <- upstream[[j]]
        if (length(up$columns)) {
          out[, j] <- drop(nuisance[, up$columns, drop = FALSE] %*% up$load)
        }
        if (target[j]) {
          out[loc, j] <- out[loc, j] + n * (z - residual_mean[j]) / denom[j]
          next
        }
        a <- used_loads[[model_ix[j]]]
        residual <- z - residual_mean[j]
        if (dr) residual <- residual - drop(w[loc, , drop = FALSE] %*% delta[j, ])
        out[loc, j] <- out[loc, j] + n * a$odds * residual / denom[j]
        out[a$selected, j] <- out[a$selected, j] + drop(a$IF %*% gamma[j, ])
        if (dr) {
          rg <- which(g_of == layout$g[j])
          x_center <- sweep(w[rg, , drop = FALSE], 2L, mu[grp[j], ], "-")
          out[rg, j] <- out[rg, j] + n / length(members[[grp[j]]]) * drop(x_center %*% delta[j, ])
          s <- layout$s[j] + (seq_len(q) - 1L) * K
          t <- layout$t[j] + (seq_len(q) - 1L) * K
          out[, j] <- out[, j] + drop((psi[, t, drop = FALSE] - psi[, s, drop = FALSE]) %*% imbalance[j, ])
        }
      }
    }
    out
  }
  score_rows <- function(rows) {
    available <- floor((entries - fixed - as.double(length(rows)) * L) / per_unit)
    if (available < 1) etwfe_tm_workspace_error()
    size <- min(block_rows, available)
    result <- matrix(0, length(rows), L)
    for (pos in split(seq_along(rows), (seq_along(rows) - 1L) %/% size)) result[pos, ] <- score_block(rows[pos])
    result
  }
  result$score_rows <- score_rows
  result
}
