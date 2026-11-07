# ==============================================================================
# S3 Result Class and Methods for Static Entry / Exit / Re-Entry (`cdid_static_att_gt`)
# ==============================================================================
#
# Provides the constructor `new_cdid_static_att_gt()` and S3 methods:
#   - `print.cdid_static_att_gt`
#   - `summary.cdid_static_att_gt` / `print.summary.cdid_static_att_gt`
#   - `tidy.cdid_static_att_gt`
#   - `glance.cdid_static_att_gt`
#   - `ggdid.cdid_static_att_gt`
#   - `plot.cdid_static_att_gt`
#
# Note: `cdid_static_att_gt` does NOT inherit from `MP` or `cdid_MP` because
# static transition targets `(t, d)` are indexed by transition date `t` and
# prior treatment status `d = D_{t-1} in {0, 1}`, not absorbing adoption cohorts
# `g`. `aggte()` remains cohort-only and refuses `cdid_static_att_gt`.
# ==============================================================================

static_result_references <- function(est_method = "dr", route = "baseline") {
  em <- if (is.character(est_method) && length(est_method) == 1L && !is.na(est_method) && nzchar(est_method)) {
    toupper(est_method)
  } else {
    "DR"
  }
  rt <- if (is.character(route) && length(route) == 1L && !is.na(route) && nzchar(route)) {
    route
  } else {
    "baseline"
  }
  list(
    estimator = paste0(
      "static no-carryover entry/exit/re-entry chained difference-in-differences (",
      em, ", ", rt, " route)"
    ),
    cdid = paste0(
      "Bellego C., Benatia D. and V. Dortet-Bernadet (2025). ",
      "\"The Chained Difference-in-Differences.\" Journal of ",
      "Econometrics 248, 105783. doi:10.1016/j.jeconom.2024.105783."
    ),
    static = paste0(
      "de Chaisemartin C., D'Haultfoeuille X., Pasquier F., Sow D. and G. Vazquez-Bare (2026). ",
      "\"Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period.\" ",
      "Working paper, arXiv:2201.06898v7 (revised 2026-07-22), static binary specialization."
    )
  )
}

new_cdid_static_att_gt <- function(target_id,
                                   t,
                                   t_prev,
                                   d,
                                   direction,
                                   sigma,
                                   available,
                                   support_status,
                                   support_reason,
                                   n_risk,
                                   n_switchers,
                                   n_stayers,
                                   att,
                                   theta_hat,
                                   se,
                                   analytical_se,
                                   bootstrap_se,
                                   V_analytical,
                                   c,
                                   pointwise_crit_val,
                                   cband_crit_val,
                                   ci_lower,
                                   ci_upper,
                                   cband_lower,
                                   cband_upper,
                                   alp,
                                   bstrap,
                                   cband,
                                   biters,
                                   n,
                                   unit_id,
                                   dates,
                                   dates_num,
                                   inffunc,
                                   aggregate,
                                   graphs,
                                   request,
                                   inference,
                                   diagnostics,
                                   fit) {
  refs <- static_result_references(request$est_method, request$route)
  if (is.list(request) && is.null(request$references)) {
    request$references <- refs
  }
  if (is.list(diagnostics) && is.null(diagnostics$references)) {
    diagnostics$references <- refs
  }

  did_params <- list(
    yname = request$yname,
    tname = request$tname,
    idname = request$idname,
    dname = request$dname,
    gname = NULL,
    xformla = request$xformla,
    panel = TRUE,
    allow_unbalanced_panel = TRUE,
    alp = alp,
    bstrap = bstrap,
    cband = cband,
    biters = biters,
    clustervars = request$clustervars,
    est_method = request$est_method,
    faster_mode = request$faster_mode,
    print_details = request$print_details,
    pl = request$pl,
    cores = request$cores,
    compute_inffunc = request$compute_inffunc,
    call = request$call,
    n = n,
    references = refs
  )

  structure(
    list(
      target_id = target_id,
      t = stats::setNames(as.double(t), target_id),
      t_prev = stats::setNames(as.double(t_prev), target_id),
      d = stats::setNames(as.integer(d), target_id),
      direction = stats::setNames(as.character(direction), target_id),
      sigma = stats::setNames(as.double(sigma), target_id),
      available = stats::setNames(as.logical(available), target_id),
      support_status = stats::setNames(as.character(support_status), target_id),
      support_reason = stats::setNames(as.character(support_reason), target_id),
      n_risk = stats::setNames(as.integer(n_risk), target_id),
      n_switchers = stats::setNames(as.integer(n_switchers), target_id),
      n_stayers = stats::setNames(as.integer(n_stayers), target_id),
      att = att,
      estimate = att,
      theta_hat = theta_hat,
      se = se,
      analytical_se = analytical_se,
      bootstrap_se = bootstrap_se,
      V_analytical = V_analytical,
      c = c,
      pointwise_crit_val = pointwise_crit_val,
      cband_crit_val = cband_crit_val,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      cband_lower = cband_lower,
      cband_upper = cband_upper,
      alp = alp,
      bstrap = bstrap,
      cband = cband,
      biters = biters,
      n = n,
      unit_id = unit_id,
      dates = dates,
      dates_num = dates_num,
      inffunc = inffunc,
      aggregate = aggregate,
      graphs = graphs,
      references = refs,
      request = request,
      inference = inference,
      diagnostics = diagnostics,
      fit = fit,
      DIDparams = did_params
    ),
    class = "cdid_static_att_gt"
  )
}

if (!exists("tidy", mode = "function")) {
  tidy <- function(x, ...) UseMethod("tidy")
}
if (!exists("glance", mode = "function")) {
  glance <- function(x, ...) UseMethod("glance")
}
if (!exists("ggdid", mode = "function")) {
  ggdid <- function(object, ...) UseMethod("ggdid")
}

tidy.cdid_static_att_gt <- function(x, ...) {
  J <- length(x$target_id)
  use_cband <- isTRUE(x$cband)
  conf_low <- if (use_cband) as.double(x$cband_lower) else as.double(x$ci_lower)
  conf_high <- if (use_cband) as.double(x$cband_upper) else as.double(x$ci_upper)

  data.frame(
    term = sprintf("ATT(t=%s, %s)", as.character(x$t), as.character(x$direction)),
    target.id = as.character(x$target_id),
    time = as.double(x$t),
    t.prev = as.double(x$t_prev),
    prior.state = as.integer(x$d),
    direction = as.character(x$direction),
    sigma = as.double(x$sigma),
    estimate = as.double(x$att),
    theta.hat = as.double(x$theta_hat),
    std.error = as.double(x$se),
    analytical.se = as.double(x$analytical_se),
    bootstrap.se = as.double(x$bootstrap_se),
    point.conf.low = as.double(x$ci_lower),
    point.conf.high = as.double(x$ci_upper),
    conf.low = conf_low,
    conf.high = conf_high,
    n.risk = as.integer(x$n_risk),
    n.switchers = as.integer(x$n_switchers),
    n.stayers = as.integer(x$n_stayers),
    available = as.logical(x$available),
    support.status = as.character(x$support_status),
    support.reason = as.character(x$support_reason),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

glance.cdid_static_att_gt <- function(x, ...) {
  data.frame(
    model = "cdid_static_att_gt",
    est.method = as.character(x$request$est_method),
    route = as.character(x$request$route),
    selection = as.character(x$request$selection),
    selection.est.method = as.character(x$request$selection_est_method),
    gmm.weighting = as.character(x$request$gmm_weighting$canonical),
    nobs = as.integer(x$n),
    n.dates = as.integer(x$diagnostics$n_dates),
    n.targets = as.integer(x$diagnostics$n_declared_targets),
    n.available = as.integer(x$diagnostics$n_available_targets),
    n.clusters = as.integer(x$inference$n_clusters),
    bstrap = as.logical(x$bstrap),
    cband = as.logical(x$cband),
    aggregate.estimate = as.double(x$aggregate$estimate),
    aggregate.se = as.double(x$aggregate$se),
    aggregate.available = as.logical(x$aggregate$available),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

summary.cdid_static_att_gt <- function(object, ...) {
  td <- tidy.cdid_static_att_gt(object)
  refs <- if (is.list(object$references)) {
    object$references
  } else {
    static_result_references(object$request$est_method, object$request$route)
  }
  structure(
    list(
      call = object$request$call,
      table = td,
      aggregate = object$aggregate,
      references = refs,
      request = object$request,
      inference = object$inference,
      diagnostics = object$diagnostics,
      graphs = object$graphs,
      alp = object$alp,
      bstrap = object$bstrap,
      cband = object$cband,
      c = object$c,
      pointwise_crit_val = object$pointwise_crit_val,
      cband_crit_val = object$cband_crit_val
    ),
    class = "summary.cdid_static_att_gt"
  )
}

static_format_CellTable <- function(s) {
  td <- s$table
  conf_pct <- round(100 * (1 - s$alp))
  iv_hdr <- if (isTRUE(s$cband)) {
    sprintf("[%d%% Simult. Band]", conf_pct)
  } else {
    sprintf("[%d%% Pointwise CI]", conf_pct)
  }

  fmt_num <- function(v) {
    ifelse(is.na(v), "NA", sprintf("%8.4f", v))
  }
  ci_str <- ifelse(
    is.na(td$conf.low) | is.na(td$conf.high),
    "              [NA, NA]",
    sprintf("[%8.4f, %8.4f]", td$conf.low, td$conf.high)
  )
  sig_star <- ifelse(
    !is.na(td$conf.low) & !is.na(td$conf.high) & (td$conf.low > 0 | td$conf.high < 0),
    "*",
    " "
  )

  out_df <- data.frame(
    Cell = td$target.id,
    Time = td$time,
    Ref = td$t.prev,
    Direction = ifelse(td$prior.state == 0L, "entry (0->1)", "exit (1->0)"),
    Switchers = td$n.switchers,
    Stayers = td$n.stayers,
    ATT = fmt_num(td$estimate),
    Std.Error = fmt_num(td$std.error),
    Interval = paste0(ci_str, " ", sig_star),
    Status = td$support.status,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(out_df)[names(out_df) == "Interval"] <- iv_hdr
  out_df
}

print.summary.cdid_static_att_gt <- function(x, ...) {
  refs <- if (is.list(x$references)) {
    x$references
  } else {
    static_result_references(x$request$est_method, x$request$route)
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\nReference: Static No-Carryover Entry/Exit/Re-Entry DiD (`dname = \"", x$request$dname, "\"`)\n", sep = "")
  cat("  Chained estimators: ", refs$cdid, "\n", sep = "")
  cat("  Static binary transition design: ", refs$static, "\n", sep = "")
  cat("  Identification benchmark: static potential outcomes Y_{it}(d_{1:t}) = Y_{it}(d_t), known binary D across all study dates\n")
  cat("  Treatment comparison method: ", toupper(x$request$est_method),
      " (requested: \"", x$request$est_method_requested, "\")\n", sep = "")
  cat("  Observation route: ", x$request$route,
      " (selection = \"", x$request$selection_requested,
      "\", selection_est_method = \"", x$request$selection_est_method_requested, "\")\n", sep = "")
  if (x$request$route %in% c("weighting", "augmented")) {
    cat("  Selection model: ", x$request$selection_model,
        " (pooling = \"", x$request$selection_pooling, "\")\n", sep = "")
  }
  if (x$request$route %in% c("imputation", "augmented")) {
    cat("  Imputation model: WLS actual-change (pooling = \"", x$request$imputation_pooling,
        "\", scale = \"", x$request$imputation_scale, "\")\n", sep = "")
  }
  if (x$request$est_method %in% c("reg", "dr")) {
    cat("  Outcome regression pooling: ", x$request$outcome_pooling_mode, "\n", sep = "")
  }
  cat("  GMM weighting: ", x$request$gmm_weighting$canonical,
      " (requested: \"", x$request$gmm_weighting$requested, "\")\n", sep = "")
  cat("  Study population: n = ", x$diagnostics$n_units, " units across T = ", x$diagnostics$n_dates,
      " dates (repeat-Y: ", x$diagnostics$n_repeat_outcome_units,
      ", single-Y: ", x$diagnostics$n_single_outcome_units,
      ", zero-Y: ", x$diagnostics$n_zero_outcome_units, ")\n", sep = "")

  cat("\nTransition-Cell Treatment Effects (tau_c = sigma_c * theta_c, sigma_c = 1 - 2 * D_{t-1}):\n")
  tbl <- static_format_CellTable(x)
  print(tbl, row.names = FALSE)
  cat("---\nSignif. codes: `*` confidence interval/band does not cover 0.\n")

  cat("\nStatic Switch-Frequency Aggregate (delta_{A, s} across declared target set A):\n")
  agg <- x$aggregate
  if (isTRUE(agg$available)) {
    cat(sprintf(
      "  Estimate: %8.4f   Std. Error: %8.4f   [%d%% Pointwise CI: %8.4f, %8.4f]   (Total switchers: %d)\n",
      agg$estimate,
      agg$se,
      round(100 * (1 - x$alp)),
      agg$ci_lower,
      agg$ci_upper,
      agg$total_switchers
    ))
  } else {
    cat("  Unavailable (", agg$status, "): ", agg$reason, "\n", sep = "")
  }
  cat("  Scope: ", agg$scope, "\n", sep = "")

  cat("\nInference & Support Summary:\n")
  cat("  Standard errors: ", x$inference$se_reported, "\n", sep = "")
  cat("  Confidence intervals/bands: ", x$inference$band_reported, "\n", sep = "")
  cat("  Clustering: ", if (isTRUE(x$inference$clustered)) {
    paste0("clustered on `", x$inference$clustername, "` (", x$inference$n_clusters, " clusters)")
  } else {
    paste0("unit-level (", x$inference$n_units, " units)")
  }, "\n", sep = "")
  cat("  Role-pair graphs: ", x$diagnostics$n_total_links,
      " total positive-count edges retained across all required switcher/stayer role blocks\n", sep = "")
  cat("  Qualifications: ", x$inference$qualifications, "\n\n", sep = "")
  invisible(x)
}

print.cdid_static_att_gt <- function(x, ...) {
  print.summary.cdid_static_att_gt(summary.cdid_static_att_gt(x, ...), ...)
  invisible(x)
}

ggdid.cdid_static_att_gt <- function(object, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    static_public_abort(
      "Package `ggplot2` is required for `ggdid.cdid_static_att_gt()`.",
      "cdid_static_plot_error"
    )
  }
  td <- tidy.cdid_static_att_gt(object)
  td_avail <- td[td$available & is.finite(td$estimate), , drop = FALSE]
  td_avail$direction_label <- ifelse(
    td_avail$prior.state == 0L,
    "Entry (0 -> 1)",
    "Exit (1 -> 0)"
  )

  p <- ggplot2::ggplot(
    td_avail,
    ggplot2::aes(
      x = .data$time,
      y = .data$estimate,
       color = .data$direction_label,
      shape = .data$direction_label
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point( position = ggplot2::position_dodge(width = 0.2), size = 2.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
      width = 0.15,
      position = ggplot2::position_dodge(width = 0.2),
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = "Static Entry / Exit Transition Treatment Effects",
      subtitle = sprintf(
        "Method: %s (%s route) | Sign-oriented effect tau_c at switch date t (ref t - 1)",
        toupper(object$request$est_method),
        object$request$route
      ),
      x = "Transition Date (t)",
      y = "Treatment Effect (tau_c)",
      color = "Transition Direction",
      shape = "Transition Direction"
    ) +
    ggplot2::theme_minimal()

  p
}

plot.cdid_static_att_gt <- function(x, ...) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    return(ggdid.cdid_static_att_gt(x, ...))
  }
  td <- tidy.cdid_static_att_gt(x)
  td_avail <- td[td$available & is.finite(td$estimate), , drop = FALSE]
  if (nrow(td_avail) == 0L) {
    graphics::plot.new()
    graphics::title("No available transition cells to plot")
    return(invisible(x))
  }
  ylim <- range(c(td_avail$estimate, td_avail$conf.low, td_avail$conf.high, 0), na.rm = TRUE)
  cols <- ifelse(td_avail$prior.state == 0L, "steelblue4", "firebrick3")
  pchs <- ifelse(td_avail$prior.state == 0L, 19L, 17L)
  graphics::plot(
    td_avail$time, td_avail$estimate,
    ylim = ylim, col = cols, pch = pchs,
    xlab = "Transition Date (t)",
    ylab = "Treatment Effect (tau_c)",
    main = "Static Entry / Exit Transition Treatment Effects"
  )
  graphics::abline(h = 0, lty = 2, col = "grey50")
  for (i in seq_len(nrow(td_avail))) {
    if (is.finite(td_avail$conf.low[i]) && is.finite(td_avail$conf.high[i])) {
      graphics::segments(
        td_avail$time[i], td_avail$conf.low[i],
        td_avail$time[i], td_avail$conf.high[i],
        col = cols[i]
      )
    }
  }
  invisible(x)
}
