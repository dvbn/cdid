# ------------------------------------------------------------------------------
# Chained Identification & Moment Chain Diagnostics (`chain_report()`)
# ------------------------------------------------------------------------------

build_mp_chain_info <- function(res, prep, request) {
  is_cp <- isTRUE(request$cross_pair)
  empty_links <- data.frame(
    link_index = integer(0L), link_id = character(0L),
    g = numeric(0L), s = numeric(0L), t = numeric(0L), ref = numeric(0L),
    link_role = character(0L), control_pool = character(0L),
    link_type = character(0L), functional_key = character(0L),
    mhat = numeric(0L), stringsAsFactors = FALSE
  )
  empty_all_links <- data.frame(
    g = numeric(0L), s = numeric(0L), t = numeric(0L),
    link_role = character(0L), control_pool = character(0L),
    t_qualified = logical(0L), stringsAsFactors = FALSE
  )

  # Extract same-pair direct auxiliary pre-trend records when present
  aux_pre <- list()
  if (!is_cp && !is.null(res$reporting) && is.data.frame(res$reporting$target_meta) &&
      !is.null(res$auxiliary_reporting) && is.data.frame(res$auxiliary_reporting$descriptor_status)) {
    tm_r <- res$reporting$target_meta
    ds <- res$auxiliary_reporting$descriptor_status
    bp <- res$reporting$base_period
    sels <- if (!is.null(bp) && !is.null(res$auxiliary_reporting$selections)) {
      res$auxiliary_reporting$selections[[bp]]
    } else list()
    sel_by_q <- list()
    for (sd in sels) {
      if (!is.null(sd$index) && !is.na(sd$index)) {
        sel_by_q[[as.character(sd$index)]] <- sd
      }
    }
    pre_rows <- which(tm_r$role == "pre" & isTRUE(tm_r$available) | (tm_r$role == "pre" & tm_r$available %in% TRUE))
    for (ri in pre_rows) {
      q <- tm_r$auxiliary_index[ri]
      if (!is.na(q) && q >= 1L && q <= nrow(ds)) {
        sd <- sel_by_q[[as.character(q)]]
        sgn <- if (!is.null(sd) && !is.null(sd$sign)) as.double(sd$sign) else 1.0
        est_i <- as.double(res$reporting$estimates[ri])
        aux_pre[[tm_r$target_id[ri]]] <- list(
          s = as.double(ds$s[q]),
          t = as.double(ds$t[q]),
          control_pool = as.character(ds$K[q]),
          sign = sgn,
          mhat = sgn * est_i
        )
      }
    }
  }

  if (is.null(res$design_used) || is.null(res$gmm_fit)) {
    return(list(
      cross_pair = is_cp,
      pt_assumption = request$pt_assumption,
      pt_start = request$pt_start,
      gmm_weighting = request$gmm_weighting,
      overid_df = 0L,
      n_moments = 0L,
      links = empty_links,
      all_links = empty_all_links,
      Lambda = NULL,
      targets = data.frame(
        target_id = character(0L), estimable = logical(0L),
        status = character(0L), support_reason = character(0L),
        estimate = numeric(0L), se_sandwich = numeric(0L),
        stringsAsFactors = FALSE
      ),
      witnesses = list(),
      aux_pre = aux_pre
    ))
  }

  DL <- res$design_used$links
  keep_cols <- intersect(
    c("link_index", "link_id", "g", "s", "t", "ref", "link_role", "cohort", "control_pool", "link_type", "functional_key"),
    names(DL)
  )
  links_df <- DL[, keep_cols, drop = FALSE]
  if (!("link_role" %in% names(links_df))) {
    links_df$link_role <- rep("did_pair", nrow(links_df))
  }
  links_df$mhat <- as.double(res$covariance_stream$mhat)

  DA <- if (!is.null(res$design_all) && is.data.frame(res$design_all$links)) res$design_all$links else DL
  all_links <- data.frame(
    g = as.double(DA$g),
    s = as.double(DA$s),
    t = as.double(DA$t),
    link_role = if ("link_role" %in% names(DA)) as.character(DA$link_role) else rep("did_pair", nrow(DA)),
    control_pool = as.character(DA$control_pool),
    t_qualified = if ("t_qualified" %in% names(DA)) as.logical(DA$t_qualified) else rep(TRUE, nrow(DA)),
    stringsAsFactors = FALSE
  )

  gfit <- res$gmm_fit
  t_cols <- intersect(c("target_id", "estimable", "status", "support_reason", "estimate", "se_sandwich"), names(gfit$targets))
  Lambda <- gfit$internal_gauge$Lambda
  witnesses <- gfit$witnesses
  targets_df <- gfit$targets[, t_cols, drop = FALSE]

  if (!is.null(res$aux_gmm_fit)) {
    afit <- res$aux_gmm_fit
    if (is.data.frame(afit$targets) && nrow(afit$targets) > 0L) {
      at_cols <- intersect(t_cols, names(afit$targets))
      new_t <- afit$targets[!(afit$targets$target_id %in% targets_df$target_id), at_cols, drop = FALSE]
      if (nrow(new_t) > 0L) {
        targets_df <- rbind(targets_df[, at_cols, drop = FALSE], new_t)
        rownames(targets_df) <- NULL
      }
    }
    if (length(afit$witnesses) > 0L) {
      new_w <- afit$witnesses[setdiff(names(afit$witnesses), names(witnesses))]
      witnesses <- c(witnesses, new_w)
    }
    aLam <- afit$internal_gauge$Lambda
    if (!is.null(aLam) && ncol(aLam) > 0L) {
      if (is.null(Lambda) || ncol(Lambda) == 0L) {
        Lambda <- aLam
        if (!is.null(rownames(gfit$internal_gauge$Lambda))) {
          rownames(Lambda) <- rownames(gfit$internal_gauge$Lambda)
        }
      } else {
        add_cols <- setdiff(colnames(aLam), colnames(Lambda))
        if (length(add_cols) > 0L) {
          a_sub <- aLam[, add_cols, drop = FALSE]
          rownames(a_sub) <- rownames(Lambda)
          Lambda <- cbind(Lambda, a_sub)
        }
      }
    }
  }

  list(
    cross_pair = is_cp,
    pt_assumption = request$pt_assumption,
    pt_start = request$pt_start,
    gmm_weighting = request$gmm_weighting,
    overid_df = gfit$overid_df,
    n_moments = nrow(links_df),
    links = links_df,
    all_links = all_links,
    Lambda = Lambda,
    targets = targets_df,
    witnesses = witnesses,
    aux_pre = aux_pre
  )
}

chain_report <- function(object, cohort = NULL, time = NULL, include_pre = FALSE, tol = 1e-10,
                         group = NULL, t = NULL) {
  if (!is.null(group)) {
    if (!is.null(cohort)) {
      cdid_abort("Supply either `cohort` or `group`, not both.", "cdid_chain_report_input_error")
    }
    cohort <- group
  }
  if (!is.null(t)) {
    if (!is.null(time)) {
      cdid_abort("Supply either `time` or `t`, not both.", "cdid_chain_report_input_error")
    }
    time <- t
  }
  if (!inherits(object, "cdid_MP") || is.null(object$chained) || is.null(object$chained$chain_info)) {
    cdid_abort(
      "`chain_report()` requires a chained `cdid_MP` result produced by `att_gt(mode = \"chained\", ...)`.",
      "cdid_chain_report_input_error"
    )
  }
  if (!is.numeric(tol) || is.object(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    cdid_abort("`tol` must be a finite nonnegative numeric scalar.", "cdid_chain_report_input_error")
  }
  if (!is.logical(include_pre) || is.object(include_pre) || length(include_pre) != 1L || is.na(include_pre)) {
    cdid_abort("`include_pre` must be TRUE or FALSE.", "cdid_chain_report_input_error")
  }

  ci <- object$chained$chain_info
  tm <- object$chained$target_meta
  att_vec <- as.double(object$att)
  se_vec <- as.double(object$se)
  links_df <- ci$links
  all_links <- if (!is.null(ci$all_links)) ci$all_links else links_df
  Lambda <- ci$Lambda
  aux_pre <- if (!is.null(ci$aux_pre)) ci$aux_pre else list()

  sel <- rep(TRUE, nrow(tm))
  if (!isTRUE(include_pre) && is.null(time)) {
    sel <- sel & tm$is_post
  }
  if (!is.null(cohort)) {
    if (!is.numeric(cohort) || is.object(cohort) || !length(cohort) || any(!is.finite(cohort))) {
      cdid_abort("`cohort` filter must be NULL or finite numeric values.", "cdid_chain_report_input_error")
    }
    sel <- sel & (tm$g %in% as.double(cohort))
  }
  if (!is.null(time)) {
    if (!is.numeric(time) || is.object(time) || !length(time) || any(!is.finite(time))) {
      cdid_abort("`time` filter must be NULL or finite numeric values.", "cdid_chain_report_input_error")
    }
    sel <- sel & (tm$time %in% as.double(time))
  }

  idx_cells <- which(sel)
  summary_rows <- vector("list", length(idx_cells))
  chains_out <- vector("list", length(idx_cells))
  names(chains_out) <- tm$target_id[idx_cells]

  lam_cols <- if (!is.null(Lambda) && ncol(Lambda) > 0L) colnames(Lambda) else character(0L)
  lam_col_idx <- if (length(lam_cols) > 0L) match(tm$target_id[idx_cells], lam_cols) else rep(NA_integer_, length(idx_cells))

  # Pre-index `all_links` for O(1) `has_direct` lookup per selected cell
  direct_env <- new.env(hash = TRUE, parent = emptyenv())
  if (nrow(all_links) > 0L) {
    if (ci$cross_pair) {
      d_keys <- unique(paste(
        fmt_design_num(all_links$g),
        all_links$link_role,
        fmt_design_num(all_links$s),
        fmt_design_num(all_links$t),
        sep = "|"
      ))
      for (dk in d_keys) direct_env[[dk]] <- TRUE
    } else {
      t_qual_ok <- if ("t_qualified" %in% names(all_links)) all_links$t_qualified else rep(TRUE, nrow(all_links))
      ok_al <- which(t_qual_ok)
      if (length(ok_al) > 0L) {
        d_keys <- unique(paste(
          fmt_design_num(all_links$g[ok_al]),
          fmt_design_num(all_links$s[ok_al]),
          fmt_design_num(all_links$t[ok_al]),
          sep = "|"
        ))
        for (dk in d_keys) direct_env[[dk]] <- TRUE
      }
    }
  }

  format_witness_path <- function(w_df, cross_pair_flag) {
    if (is.null(w_df) || nrow(w_df) == 0L) return("")
    parts <- vapply(seq_len(nrow(w_df)), function(r) {
      sgn <- if (w_df$sign[r] >= 0) "+" else "-"
      tag <- if (cross_pair_flag) {
        if (identical(w_df$link_role[r], "target")) paste0("g=", fmt_design_num(w_df$g[r])) else paste0("ctrl", w_df$control_pool[r])
      } else {
        paste0("K=", w_df$control_pool[r])
      }
      paste0(if (r == 1L && sgn == "+") "" else paste0(sgn, " "),
             fmt_design_num(w_df$s[r]), "->", fmt_design_num(w_df$t[r]), " [", tag, "]")
    }, character(1L))
    paste(parts, collapse = " ")
  }

  for (k in seq_along(idx_cells)) {
    i <- idx_cells[k]
    tid <- tm$target_id[i]
    gv <- tm$g[i]
    tv <- tm$time[i]
    rv <- tm$ref[i]
    bt <- if ("base_t" %in% names(tm)) tm$base_t[i] else rv
    role_i <- tm$role[i]
    av <- isTRUE(tm$available[i])

    wit_raw <- ci$witnesses[[tid]]
    w_df <- data.frame(
      target_id = character(0L), g = numeric(0L), target_time = numeric(0L),
      link_index = integer(0L), s = numeric(0L), t = numeric(0L),
      link_role = character(0L), control_pool = character(0L),
      sign = numeric(0L), link_estimate = numeric(0L), contribution = numeric(0L),
      stringsAsFactors = FALSE
    )
    gmm_df <- data.frame(
      target_id = character(0L), g = numeric(0L), target_time = numeric(0L),
      link_index = integer(0L), s = numeric(0L), t = numeric(0L),
      link_role = character(0L), control_pool = character(0L),
      gmm_weight = numeric(0L), link_estimate = numeric(0L), contribution = numeric(0L),
      stringsAsFactors = FALSE
    )

    if (!is.null(wit_raw) && nrow(wit_raw) > 0L) {
      w_col <- if ("weight" %in% names(wit_raw)) wit_raw$weight else wit_raw$sign_or_weight
      li <- wit_raw$link_index
      sub_l <- links_df[li, , drop = FALSE]
      ord_w <- order(match(sub_l$link_role, c("target", "did_pair", "control")), sub_l$s, sub_l$t)
      sub_l <- sub_l[ord_w, , drop = FALSE]
      w_col <- as.double(w_col[ord_w])
      w_df <- data.frame(
        target_id = rep(tid, nrow(sub_l)),
        g = rep(gv, nrow(sub_l)),
        target_time = rep(tv, nrow(sub_l)),
        link_index = as.integer(sub_l$link_index),
        s = as.double(sub_l$s),
        t = as.double(sub_l$t),
        link_role = as.character(sub_l$link_role),
        control_pool = as.character(sub_l$control_pool),
        sign = w_col,
        link_estimate = as.double(sub_l$mhat),
        contribution = w_col * as.double(sub_l$mhat),
        stringsAsFactors = FALSE
      )
    } else if (!ci$cross_pair && identical(role_i, "pre") && av && !is.null(aux_pre[[tid]])) {
      ap <- aux_pre[[tid]]
      w_df <- data.frame(
        target_id = tid,
        g = gv,
        target_time = tv,
        link_index = NA_integer_,
        s = ap$s,
        t = ap$t,
        link_role = "auxiliary_did_pair",
        control_pool = ap$control_pool,
        sign = ap$sign,
        link_estimate = ap$mhat,
        contribution = ap$sign * ap$mhat,
        stringsAsFactors = FALSE
      )
      gmm_df <- data.frame(
        target_id = tid,
        g = gv,
        target_time = tv,
        link_index = NA_integer_,
        s = ap$s,
        t = ap$t,
        link_role = "auxiliary_did_pair",
        control_pool = ap$control_pool,
        gmm_weight = ap$sign,
        link_estimate = ap$mhat,
        contribution = ap$sign * ap$mhat,
        stringsAsFactors = FALSE
      )
    }

    lc_k <- lam_col_idx[k]
    if (!is.na(lc_k)) {
      lam <- as.double(Lambda[, lc_k])
      act <- which(abs(lam) > tol)
      if (length(act) > 0L) {
        sub_l <- links_df[act, , drop = FALSE]
        lam_act <- lam[act]
        ord_g <- order(sub_l$g != gv, match(sub_l$link_role, c("target", "did_pair", "control")), sub_l$s, sub_l$t)
        sub_l <- sub_l[ord_g, , drop = FALSE]
        lam_act <- lam_act[ord_g]
        gmm_df <- data.frame(
          target_id = rep(tid, nrow(sub_l)),
          g = rep(gv, nrow(sub_l)),
          target_time = rep(tv, nrow(sub_l)),
          link_index = as.integer(sub_l$link_index),
          s = as.double(sub_l$s),
          t = as.double(sub_l$t),
          link_role = as.character(sub_l$link_role),
          control_pool = as.character(sub_l$control_pool),
          gmm_weight = lam_act,
          link_estimate = as.double(sub_l$mhat),
          contribution = lam_act * as.double(sub_l$mhat),
          stringsAsFactors = FALSE
        )
      }
    }

    # Check whether a direct single-pair comparison between base_t and target_time
    # exists in the unreduced comparison catalog `all_links`.
    ref_cmp <- if (!is.na(bt)) bt else rv
    lo_t <- min(ref_cmp, tv)
    hi_t <- max(ref_cmp, tv)
    has_direct <- if (identical(role_i, "reference") || is.na(ref_cmp) || ref_cmp == tv) {
      TRUE
    } else if (ci$cross_pair) {
      g_s <- fmt_design_num(gv); lo_s <- fmt_design_num(lo_t); hi_s <- fmt_design_num(hi_t)
      isTRUE(direct_env[[paste(g_s, "target", lo_s, hi_s, sep = "|")]]) &&
        isTRUE(direct_env[[paste(g_s, "control", lo_s, hi_s, sep = "|")]])
    } else if (identical(role_i, "pre")) {
      TRUE
    } else {
      isTRUE(direct_env[[paste(fmt_design_num(gv), fmt_design_num(lo_t), fmt_design_num(hi_t), sep = "|")]])
    }

    summary_rows[[k]] <- data.frame(
      target_id = tid,
      g = gv,
      time = tv,
      event_time = tv - gv,
      ref = rv,
      role = role_i,
      available = av,
      status = tm$support_status[i],
      requires_chaining = av && !has_direct,
      n_witness_links = nrow(w_df),
      n_active_gmm_links = nrow(gmm_df),
      estimate = att_vec[i],
      se = se_vec[i],
      chain_path = format_witness_path(w_df, ci$cross_pair),
      support_reason = tm$support_reason[i],
      stringsAsFactors = FALSE
    )
    chains_out[[k]] <- list(witness = w_df, gmm_weights = gmm_df)
  }

  sum_df <- if (length(summary_rows)) {
    out_s <- do.call(rbind, summary_rows)
    rownames(out_s) <- NULL
    out_s
  } else {
    data.frame(
      target_id = character(0L), g = numeric(0L), time = numeric(0L),
      event_time = numeric(0L), ref = numeric(0L), role = character(0L),
      available = logical(0L), status = character(0L), requires_chaining = logical(0L),
      n_witness_links = integer(0L), n_active_gmm_links = integer(0L),
      estimate = numeric(0L), se = numeric(0L), chain_path = character(0L),
      support_reason = character(0L), stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      summary = sum_df,
      chains = chains_out,
      meta = list(
        cross_pair = ci$cross_pair,
        pt_assumption = ci$pt_assumption,
        pt_start = ci$pt_start,
        gmm_weighting = ci$gmm_weighting,
        overid_df = ci$overid_df,
        n_moments = ci$n_moments,
        tol = tol
      )
    ),
    class = "cdid_chain_report"
  )
}

as.data.frame.cdid_chain_report <- function(x, row.names = NULL, optional = FALSE,
                                            what = c("summary", "witness", "gmm_weights"), ...) {
  if (missing(what)) what <- "summary"
  if (!is.character(what) || is.object(what) || length(what) != 1L || is.na(what) ||
      !(what %in% c("summary", "witness", "gmm_weights"))) {
    cdid_abort("`what` must be one of \"summary\", \"witness\", or \"gmm_weights\".", "cdid_chain_report_input_error")
  }
  if (identical(what, "summary")) {
    out <- x$summary
    if (!is.null(row.names)) rownames(out) <- row.names
    return(out)
  }
  parts <- lapply(x$chains, `[[`, what)
  parts <- parts[vapply(parts, function(d) is.data.frame(d) && nrow(d) > 0L, logical(1L))]
  if (!length(parts)) {
    return(if (identical(what, "witness")) {
      data.frame(target_id = character(0L), g = numeric(0L), target_time = numeric(0L),
                 link_index = integer(0L), s = numeric(0L), t = numeric(0L),
                 link_role = character(0L), control_pool = character(0L),
                 sign = numeric(0L), link_estimate = numeric(0L), contribution = numeric(0L),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(target_id = character(0L), g = numeric(0L), target_time = numeric(0L),
                 link_index = integer(0L), s = numeric(0L), t = numeric(0L),
                 link_role = character(0L), control_pool = character(0L),
                 gmm_weight = numeric(0L), link_estimate = numeric(0L), contribution = numeric(0L),
                 stringsAsFactors = FALSE)
    })
  }
  out <- do.call(rbind, parts)
  rownames(out) <- row.names
  out
}

print.cdid_chain_report <- function(x, digits = 4, show_gmm_weights = FALSE, max_cells = 30L, ...) {
  if (!is.numeric(digits) || is.object(digits) || length(digits) != 1L ||
      !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    cdid_abort("`digits` must be a finite nonnegative integer scalar.", "cdid_chain_report_input_error")
  }
  if (!is.logical(show_gmm_weights) || is.object(show_gmm_weights) ||
      length(show_gmm_weights) != 1L || is.na(show_gmm_weights)) {
    cdid_abort("`show_gmm_weights` must be TRUE or FALSE.", "cdid_chain_report_input_error")
  }
  if (!is.numeric(max_cells) || is.object(max_cells) || length(max_cells) != 1L ||
      !is.finite(max_cells) || max_cells < 1 || max_cells != floor(max_cells)) {
    cdid_abort("`max_cells` must be a finite positive integer scalar.", "cdid_chain_report_input_error")
  }
  digits <- as.integer(digits)
  m <- x$meta
  cat("Chained DiD Identification & Moment Chain Report\n")
  cat("================================================\n")
  cat(sprintf(
    "Mode: %s | PT: %s | GMM weighting: %s | Moments: %d | Overid df: %d\n",
    if (isTRUE(m$cross_pair)) "cross-pair" else "same-pair",
    m$pt_assumption, m$gmm_weighting, m$n_moments, m$overid_df
  ))
  s <- x$summary
  n_tot <- nrow(s)
  n_av <- sum(s$available)
  n_ch <- sum(s$requires_chaining, na.rm = TRUE)
  cat(sprintf(
    "Selected cells: %d (%d estimable, %d chained-only without single-cohort span, %d unavailable)\n\n",
    n_tot, n_av, n_ch, n_tot - n_av
  ))
  if (n_tot == 0L) return(invisible(x))

  n_show <- min(n_tot, as.integer(max_cells))
  for (i in seq_len(n_show)) {
    row <- s[i, ]
    est_str <- if (is.na(row$estimate)) "NA" else format(round(row$estimate, digits), nsmall = digits)
    se_str <- if (is.na(row$se)) "NA" else format(round(row$se, digits), nsmall = digits)
    tag <- if (!row$available) "UNAVAILABLE" else if (row$requires_chaining) "CHAINED" else "DIRECT+CHAIN"
    cat(sprintf(
      "[%s] %s (ref=%s, e=%s): est = %s (SE = %s)\n",
      tag, row$target_id, fmt_design_num(row$ref), fmt_design_num(row$event_time), est_str, se_str
    ))
    if (row$available && nzchar(row$chain_path)) {
      cat("  Witness chain: ", row$chain_path, "\n", sep = "")
    } else if (!row$available) {
      cat("  Reason: ", row$support_reason, "\n", sep = "")
    }
    if (isTRUE(show_gmm_weights) && row$available) {
      gw <- x$chains[[row$target_id]]$gmm_weights
      if (nrow(gw) > 0L) {
        gw_disp <- gw[, c("s", "t", "link_role", "control_pool", "gmm_weight", "link_estimate", "contribution"), drop = FALSE]
        gw_disp$gmm_weight <- round(gw_disp$gmm_weight, digits)
        gw_disp$link_estimate <- round(gw_disp$link_estimate, digits)
        gw_disp$contribution <- round(gw_disp$contribution, digits)
        print(gw_disp, row.names = FALSE)
      }
    }
  }
  if (n_tot > n_show) {
    cat(sprintf("... (%d additional cells omitted; increase `max_cells` or filter by `cohort`/`time`)\n", n_tot - n_show))
  }
  invisible(x)
}
