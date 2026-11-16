# Eliminate the unit effects before combining the slot-by-feature blocks.
# Each weighted mask product supplies one block of the same within-OLS Gram.
# Reuse its factor for the regression and joint original-unit score transfer.
etwfe_covariate_core <- function(y, unit, period, cohort, n_units, n_periods,
                                 control_group, Z, compute_inffunc = TRUE,
                                 workspace_entries, report = NULL) {
  never <- etwfe_core_contract(y, unit, period, cohort, n_units, n_periods,
                               control_group, compute_inffunc)
  if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != n_units || ncol(Z) < 1L ||
      any(!is.finite(Z))) {
    cdid_abort("ETWFE covariate model requires a finite matrix on the full unit frame.",
               "cdid_error_etwfe_contract")
  }
  workspace_error <- function() {
    cdid_abort("ETWFE covariate model exceeds the available numeric workspace.",
               "cdid_etwfe_workspace_error")
  }
  rank_error <- function() {
    cdid_abort(paste("ETWFE covariate model could not resolve full rank after",
                     "fixing the graph gauges. No estimates are returned."),
               "cdid_error_etwfe_numeric")
  }
  if (!is.numeric(workspace_entries) || length(workspace_entries) != 1L ||
      is.na(workspace_entries) || workspace_entries < 1) workspace_error()

  treated <- sort(unique(cohort[cohort > 0L]))
  n_group <- length(treated)
  group <- match(cohort, treated, nomatch = 0L) + 1L
  dates <- seq_len(n_periods)
  n_full <- as.double(n_periods) * (n_group + 1)
  if (n_full > .Machine$integer.max || 6 * n_full > workspace_entries) {
    workspace_error()
  }
  free <- if (never) outer(dates, treated - 1L, `!=`) else
    outer(dates, treated, `>=`)
  raw_slot <- matrix(dates, n_periods, n_group + 1L)
  if (n_group) {
    cell <- outer(dates, seq_len(n_group) * n_periods, `+`)
    raw_slot[, -1L] <- ifelse(free, cell, dates)
  }
  seen <- tabulate(raw_slot[cbind(period, group[unit])], n_full) > 0L
  n_slot <- sum(seen)
  slot_lookup <- integer(n_full)
  slot_lookup[seen] <- seq_len(n_slot)
  slot_by_date <- matrix(slot_lookup[raw_slot], n_periods, n_group + 1L)
  target <- which(free, arr.ind = TRUE)
  target_date <- as.integer(target[, 1L])
  target_group <- target[, 2L] + 1L
  n_target <- nrow(target)
  cell_slot <- slot_by_date[cbind(target_date, target_group)]
  date_slot <- slot_by_date[target_date]

  n_feature <- ncol(Z) + 1L
  n_coef <- as.double(n_slot) * n_feature
  fixed <- 8 * n_coef^2 + 6 * n_slot^2 + 10 * as.double(n_periods)^2 +
    8 * n_full + 4 * length(y) + n_units * (12 + 3 * n_feature) +
    n_target * (12 + 3 * n_feature)
  if (compute_inffunc) fixed <- fixed + n_target * (2 * n_units + 4 * n_coef)
  per_unit <- 12 * n_periods + 6 * n_feature + 2 * n_target * compute_inffunc
  if (n_coef > .Machine$integer.max || fixed + per_unit > workspace_entries) {
    workspace_error()
  }
  # The limit also applies to a wide full-cohort IF operation. Only the final
  # roster-aligned influence matrix holds all unit-by-target entries.
  block_size <- max(1, floor(min(workspace_entries - fixed, 2^22) / per_unit))
  blocks <- function(ids) {
    index <- seq_along(ids)
    split(ids, (index - 1L) %/% block_size)
  }
  n_obs <- tabulate(unit, n_units)
  first <- cumsum(c(1, as.double(n_obs)))[seq_len(n_units)]
  # A change of coordinate for all slot coefficients. This avoids losing
  # within-feature variation when a covariate has a large common offset.
  center <- colMeans(Z)
  scaled <- sweep(Z, 2L, center, `-`)
  amplitude <- vapply(seq_len(ncol(Z)), function(a) max(abs(scaled[, a])), 0.0)
  if (any(!is.finite(amplitude))) rank_error()
  amplitude[amplitude == 0] <- 1
  features <- cbind(1, sweep(scaled, 2L, amplitude, `/`))
  rm(scaled)
  roster <- split(seq_len(n_units), group)
  cohort_mean <- matrix(NA_real_, n_group + 1L, n_feature)
  for (g in names(roster)) {
    cohort_mean[as.integer(g), ] <- colMeans(features[roster[[g]], , drop = FALSE])
  }
  repeat_ids <- which(n_obs > 1L)
  fit_groups <- split(repeat_ids, group[repeat_ids])
  read_block <- function(ids, complete) {
    rows <- rep(first[ids] - 1, n_obs[ids]) + sequence(n_obs[ids])
    if (complete) {
      outcome <- matrix(y[rows], length(ids), n_periods, byrow = TRUE)
      return(list(within = outcome - rowMeans(outcome), mask = NULL))
    }
    coord <- cbind(rep(seq_along(ids), n_obs[ids]), period[rows])
    mask <- outcome <- matrix(0, length(ids), n_periods)
    mask[coord] <- 1
    outcome[coord] <- y[rows]
    list(within = outcome - mask * (rowSums(outcome) / n_obs[ids]), mask = mask)
  }
  gram <- matrix(0, n_coef, n_coef)
  rhs <- matrix(0, n_slot, n_feature)
  link <- matrix(FALSE, n_slot, n_slot)
  within_complete <- diag(n_periods) - 1 / n_periods
  coef_index <- matrix(seq_len(n_coef), n_slot, n_feature)
  for (g in as.integer(names(fit_groups))) {
    group_ids <- fit_groups[[as.character(g)]]
    s <- slot_by_date[, g]
    used <- which(s > 0L)
    s <- s[used]
    for (complete in c(TRUE, FALSE)) {
      for (ids in blocks(group_ids[(n_obs[group_ids] == n_periods) == complete])) {
        block <- read_block(ids, complete)
        w <- features[ids, , drop = FALSE]
        rhs[s, ] <- rhs[s, , drop = FALSE] + crossprod(block$within[, used, drop = FALSE], w)
        if (complete) {
          link[s, s] <- TRUE
          feature_gram <- crossprod(w)
        } else {
          mask <- block$mask[, used, drop = FALSE]
          link[s, s] <- link[s, s] | (crossprod(mask) > 0)
        }
        for (a in seq_len(n_feature)) {
          i <- coef_index[s, a]
          for (b in seq_len(a)) {
            j <- coef_index[s, b]
            if (complete) {
              contribution <- within_complete * feature_gram[a, b]
            } else {
              weight <- w[, a] * w[, b]
              contribution <- -crossprod(mask, mask * (weight / n_obs[ids]))
              diag(contribution) <- diag(contribution) + drop(crossprod(mask, weight))
            }
            gram[i, j] <- gram[i, j, drop = FALSE] + contribution
            if (a != b) gram[j, i] <- gram[j, i, drop = FALSE] + t(contribution)
          }
        }
      }
    }
  }

  component <- etwfe_core_components(link)
  supported <- cell_slot > 0L
  identified <- supported & date_slot > 0L
  identified[identified] <- component[cell_slot[identified]] == component[date_slot[identified]]
  id <- which(identified)
  # The coefficient system is block diagonal across slot components.
  # Components without a reported contrast do not enter any of its scores.
  kept_slot <- duplicated(component) & component %in% component[cell_slot[id]]
  keep <- as.vector(coef_index[kept_slot, , drop = FALSE])
  status <- rep("no_support", n_target)
  status[supported] <- "unidentified"
  status[id] <- "identified"
  beta <- numeric(n_coef)
  if (length(keep)) {
    G <- gram[keep, keep, drop = FALSE]
    if (any(!is.finite(G)) || any(diag(G) <= 0)) rank_error()
    scale <- sqrt(diag(G))
    G <- t(t(G / scale) / scale)
    R <- tryCatch(chol(G), error = function(e) NULL)
    rc <- if (is.null(R)) NA_real_ else rcond(R, triangular = TRUE)
    if (!is.finite(rc) || rc^2 <= .Machine$double.eps) rank_error()
    solve_scaled <- function(b) {
      backsolve(R, backsolve(R, b / scale, transpose = TRUE)) / scale
    }
    beta[keep] <- solve_scaled(as.vector(rhs)[keep])
  }
  feature_names <- if (!is.null(colnames(Z))) c("(Intercept)", colnames(Z))
  contrast <- matrix(NA_real_, n_target, n_feature, dimnames = list(NULL, feature_names))
  beta <- matrix(beta, n_slot, n_feature)
  scaled_contrast <- beta[cell_slot[id], , drop = FALSE] - beta[date_slot[id], , drop = FALSE]
  mix <- cohort_mean[target_group, , drop = FALSE]
  att <- rep(NA_real_, n_target)
  att[id] <- rowSums(scaled_contrast * mix[id, , drop = FALSE])
  raw_slopes <- sweep(scaled_contrast[, -1L, drop = FALSE], 2L, amplitude, `/`)
  contrast[id, ] <- cbind(scaled_contrast[, 1L] - drop(raw_slopes %*% center), raw_slopes)
  target_cohort <- treated[target_group - 1L]
  targets <- data.frame(cohort = as.integer(target_cohort), period = target_date,
                        is_pre = target_date < target_cohort,
                        att = att, status = status)
  result <- list(targets = targets, inffunc = NULL, contrast = contrast)
  if (!is.null(report)) {
    slots <- etwfe_nuisance_slot_codes(which(seen), n_periods, c(0L, treated))
    used_component <- component %in% component[cell_slot[id]]
    slots$gauge <- !duplicated(component) & used_component
    slot <- slot_by_date[cbind(period, group[unit])]
    score <- function(rows) etwfe_nuisance_within_score(rows, y, slot, features, beta, n_obs, first, keep, solve_scaled)
    used_units <- repeat_ids[used_component[slot[first[repeat_ids]]]]
    result$nuisance <- list(outcome = etwfe_nuisance_outcome(beta, keep, coef_index, slots, Z, center, amplitude,
      n_units, report, score, workspace_entries - fixed - length(slot), 4 * n_coef + 12 * n_periods * n_feature,
      sum(n_obs[used_units]), "levels", "Observed rows on repeat-observed units in fitted components",
      "Within-OLS score and fitted slot factor on original units."))
  }
  if (!compute_inffunc) return(result)

  influence <- matrix(0, n_units, n_target)
  influence[, !identified] <- NA_real_
  if (!length(id)) {
    result$inffunc <- influence
    return(result)
  }
  C <- transfer <- matrix(0, n_coef, length(id))
  for (a in seq_len(n_feature)) {
    C[cbind(coef_index[cell_slot[id], a], seq_along(id))] <- mix[id, a]
    C[cbind(coef_index[date_slot[id], a], seq_along(id))] <- -mix[id, a]
  }
  transfer[keep, ] <- solve_scaled(C[keep, , drop = FALSE])
  padded_beta <- rbind(0, beta)
  for (g in as.integer(names(fit_groups))) {
    s <- slot_by_date[, g]
    used <- which(s > 0L)
    group_coef <- padded_beta[s + 1L, , drop = FALSE]
    group_ids <- fit_groups[[as.character(g)]]
    for (complete in c(TRUE, FALSE)) {
      for (ids in blocks(group_ids[(n_obs[group_ids] == n_periods) == complete])) {
        block <- read_block(ids, complete)
        w <- features[ids, , drop = FALSE]
        fitted <- tcrossprod(w, group_coef)
        within_fit <- if (complete) fitted - rowMeans(fitted) else {
          fitted <- fitted * block$mask
          fitted - block$mask * (rowSums(fitted) / n_obs[ids])
        }
        residual <- (block$within - within_fit)[, used, drop = FALSE]
        score <- matrix(0, length(ids), length(id))
        for (a in seq_len(n_feature)) {
          load <- transfer[coef_index[s[used], a], , drop = FALSE]
          score <- score + (residual %*% load) * w[, a]
        }
        influence[ids, id] <- n_units * score
      }
    }
  }
  for (g in unique(target_group[id])) {
    columns <- id[target_group[id] == g]
    group_ids <- roster[[as.character(g)]]
    size_factor <- n_units / length(group_ids)
    for (ids in blocks(group_ids)) {
      centered <- sweep(features[ids, , drop = FALSE], 2L, cohort_mean[g, ], `-`)
      influence[ids, columns] <- influence[ids, columns, drop = FALSE] +
        size_factor * tcrossprod(centered, scaled_contrast[match(columns, id), , drop = FALSE])
    }
  }
  result$inffunc <- influence
  result
}
