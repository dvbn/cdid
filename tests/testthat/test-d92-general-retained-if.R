# SPDX-License-Identifier: GPL-3.0-only
#
# D92 point-IPW performance revision: invocation-local raw full-IF carrier
# (`ipw_point_if_retained_*`) extended to the sealed general IPW plan, all
# static cohort/missingness/weight/cluster/F>=A branches, and static custom IPW
# (`custom_ipw_if_plan`).
#
# Expected numeric truth comes from forced no-carrier executions and standalone
# `replay_ipw_point_if_view` / `replay_custom_ipw_if_view` / 4-term component
# sums. The retained IF carrier never supplies its own truth.
#
# STATUS: UNEXECUTED at authoring time (no R fit or test execution in the
# authoring job). The root sequences the independent installed gate.

library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) {
  stop("Actual installed cdid required; no source fallback.", call. = FALSE)
}
api <- asNamespace("cdid")

d92_ns <- function(edits = list(), wrap = character()) {
  e <- new.env(parent = api)
  counts <- new.env(parent = emptyenv())
  genuine <- "gmm_require_native_dll"
  for (nm in ls(api, all.names = TRUE)) {
    f <- get(nm, envir = api, inherits = FALSE)
    if (is.function(f) && !is.primitive(f) && identical(environment(f), api) && !nm %in% genuine) {
      environment(f) <- e
      assign(nm, f, envir = e)
    }
  }
  for (nm in genuine) assign(nm, get(nm, envir = api, inherits = FALSE), envir = e)
  for (nm in names(edits)) assign(nm, edits[[nm]], envir = e)
  make_wrap <- function(nm, inner) {
    force(nm); force(inner)
    function(...) {
      assign(nm, get(nm, envir = counts, inherits = FALSE) + 1L, envir = counts)
      inner(...)
    }
  }
  for (nm in wrap) {
    if (!exists(nm, envir = e, inherits = FALSE)) stop("No internal function `", nm, "` to count.", call. = FALSE)
    assign(nm, 0L, envir = counts)
    assign(nm, make_wrap(nm, get(nm, envir = e, inherits = FALSE)), envir = e)
  }
  e$.counts <- counts
  e
}

d92_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0L)
d92_reset <- function(e) for (nm in ls(e$.counts, all.names = TRUE)) assign(nm, 0L, envir = e$.counts)

d92_strip_env <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.function(x)) return("<function>")
  if (is.list(x)) { a <- attributes(x); y <- lapply(x, d92_strip_env); attributes(y) <- a; return(y) }
  x
}

d92_expect_match <- function(a, b, tol = 1e-12, info = NULL) {
  mismatch <- function(x, y, path) {
    if (!identical(class(x), class(y)) || !identical(typeof(x), typeof(y)) ||
        !identical(names(x), names(y)) || !identical(dim(x), dim(y)) ||
        !identical(dimnames(x), dimnames(y))) return(paste(path, "schema"))
    if (is.list(x)) {
      if (length(x) != length(y)) return(paste(path, "length"))
      return(unlist(lapply(seq_along(x), function(i) mismatch(x[[i]], y[[i]],
        paste0(path, "$", if (is.null(names(x))) i else names(x)[i]))), use.names = FALSE))
    }
    same <- if (typeof(x) == "double") length(waldo::compare(x, y, tolerance = tol)) == 0L else identical(x, y)
    if (same) character() else paste(path, "value")
  }
  bad <- mismatch(a, b, "root")
  expect_true(length(bad) == 0L, info = paste(c(info, head(bad, 12L)), collapse = "\n"))
}

d92_rebind <- function(env, nm, value) {
  old <- get(nm, envir = env, inherits = FALSE)
  unlockBinding(nm, env); assign(nm, value, envir = env); lockBinding(nm, env)
  function() { unlockBinding(nm, env); assign(nm, old, envir = env); lockBinding(nm, env) }
}

d92_panel <- function(N = 41L, T = 6L, pattern = "half", multi_cohort = TRUE,
                      weighted = FALSE, clustered = FALSE, extra_singletons = 0L) {
  u <- seq_len(N)
  x1 <- (u %% 5L) / 4 - 0.5
  x2 <- ((3L * u) %% 7L) / 6 - 0.5
  g <- ifelse(((u %% 3L) == 0L) | (x1 > 0 & (u %% 2L) == 0L), 4, Inf)
  if (multi_cohort) g[c(5L, 11L, 23L, 29L, 35L)] <- 6
  d <- expand.grid(t = seq_len(T), id = u, KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$g <- g[d$id]
  d$x1 <- x1[d$id]
  d$x2 <- x2[d$id]
  d$y <- 0.3 * d$id + 0.2 * d$t + ifelse(d$t >= d$g, 1.0, 0.0) + 0.5 * d$x1 * d$t - 0.25 * d$x2 +
    ((7 * d$id + 3 * d$t) %% 5) / 10
  d$w <- if (weighted) 0.5 + (d$id %% 5L) * 0.35 else rep(1.0, nrow(d))
  if (clustered) d$cluster <- as.integer((d$id - 1L) %% 9L + 1L)
  keep <- rep(TRUE, nrow(d))
  if (identical(pattern, "half")) keep <- !((d$id %% 2L == 0L) & (((d$t + d$id %/% 2L) %% 2L) == 0L))
  if (identical(pattern, "mixed")) keep <- ((3L * d$id + 5L * d$t) %% 7L) != 0L
  d <- d[keep, , drop = FALSE]
  if (extra_singletons > 0L) {
    # Singleton units observed at t = 1 only: in full frame F (n_total > n_a_units) but not in A.
    s_id <- N + seq_len(extra_singletons)
    s_df <- data.frame(
      id = s_id,
      t = 1L,
      g = rep(c(4, Inf), length.out = extra_singletons),
      x1 = 0.25,
      x2 = -0.25,
      y = 0.1 * s_id,
      w = if (weighted) 1.25 else 1.0
    )
    if (clustered) s_df$cluster <- as.integer((s_id - 1L) %% 9L + 1L)
    d <- rbind(d, s_df)
  }
  rownames(d) <- NULL
  d
}

d92_formula <- function(p) switch(as.character(p), "1" = ~1, "2" = ~x1, "3" = ~x1 + x2, stop("p must be 1, 2 or 3"))

d92_ctrl <- function(d, p = 1L, chunk = 7L, weighting = "optimal") {
  gs <- unique(d$g)
  n_frame <- length(unique(d$id))
  dims <- api$chained_public_dimensions_from_counts(
    n_frame = n_frame,
    n_panel_units = n_frame,
    n_rows = nrow(d),
    n_dates = length(unique(d$t)),
    n_cohorts = length(gs),
    n_targets = sum(is.finite(gs)),
    n_covariates = p,
    biters = 99,
    n_clusters = n_frame
  )
  api$chained_public_pipeline_control(
    dims,
    api$chained_public_resolve_control(list(chunk_rows = chunk), gmm_weighting = weighting),
    gmm_weighting = weighting
  )$control
}

d92_run <- function(ns, d, ctrl, p = 1L, cg = "notyettreated", weighted = FALSE, clustered = FALSE, ...) {
  ns$run_chained_ipw_pipeline(
    data = d,
    idname = "id",
    tname = "t",
    yname = "y",
    gname = "g",
    control = ctrl,
    xformla = d92_formula(p),
    weightsname = if (weighted) "w" else NULL,
    clustername = if (clustered) "cluster" else NULL,
    control_group = cg,
    base_period = "varying",
    weighting = "optimal",
    bstrap = FALSE,
    cband = FALSE,
    run_numerical_transfer = FALSE,
    root_witness = NULL,
    structural_reduction_stage = "outcome_free",
    ...
  )
}

d92_tracked <- c(
  "ipw_point_if_retained_reserve", "ipw_point_if_retained_check",
  "ipw_point_if_retained_fill", "replay_ipw_point_if_retained",
  "ipw_point_if_retained_release", "ipw_point_if_general_state_reserve",
  "ipw_point_if_general_state_build", "replay_ipw_point_if_general",
  "stream_ipw_point_covariance_general", "query_ipw_view_chunk",
  "custom_ipw_if_plan", "custom_ipw_if_plan_verify",
  "replay_custom_ipw_if_plan", "stream_custom_ipw_covariance_plan"
)

d92_no_carrier <- function(extra_edits = list(), wrap = d92_tracked) {
  d92_ns(
    c(list(
      ipw_point_if_retained_reserve = function(...) NULL,
      ipw_point_if_general_state_reserve = function(...) NULL
    ), extra_edits),
    wrap = wrap
  )
}

d92_expect_parity <- function(a, b, info = NULL) {
  for (f in c("support", "basis", "aggregations", "gmm_fit", "structural_certificate",
              "design_all", "design_used", "covariance_stream", "covariance_handling",
              "full_frame_transport", "auxiliary_reporting")) {
    d92_expect_match(d92_strip_env(a[[f]]), d92_strip_env(b[[f]]), info = paste(info, f))
  }
  for (nm in c("mhat", "V", "mean_phi", "n_a_units")) {
    expect_identical(a$covariance_stream[[nm]], b$covariance_stream[[nm]], info = paste(info, "cov", nm))
  }
  for (nm in c("estimates", "se", "analytical_se", "z_stat", "p_value", "ci_lower", "ci_upper", "if_full", "vcov", "available")) {
    expect_identical(a$reporting[[nm]], b$reporting[[nm]], info = paste(info, "reporting", nm))
  }
  expect_identical(a$gmm_fit$if_projection$n_batches, b$gmm_fit$if_projection$n_batches, info = info)
  expect_identical(a$gmm_fit$if_projection$max_batch_rows_seen, b$gmm_fit$if_projection$max_batch_rows_seen, info = info)
}

# ---------------------------------------------------------------------------
# 1. Pipeline parity & spy verification across missing multi-cohort NYT/NT,
#    X/intercept, weights, clustering, and F > A (general + batch plans).
# ---------------------------------------------------------------------------
test_that("retained IF carrier matches forced no-carrier across missing multi-cohort NYT/NT, X/intercept, weights, clustering and F > A with 1 producer pass and 3 slice passes", {
  cases <- list(
    list(pattern = "half", multi = TRUE, p = 2L, cg = "notyettreated", wt = TRUE, cl = TRUE, extra = 3L),
    list(pattern = "mixed", multi = TRUE, p = 1L, cg = "nevertreated", wt = FALSE, cl = FALSE, extra = 0L),
    list(pattern = "half", multi = FALSE, p = 3L, cg = "notyettreated", wt = TRUE, cl = FALSE, extra = 2L),
    list(pattern = "balanced", multi = FALSE, p = 1L, cg = "nevertreated", wt = FALSE, cl = TRUE, extra = 2L)
  )
  for (cs in cases) {
    d <- d92_panel(pattern = cs$pattern, multi_cohort = cs$multi,
                   weighted = cs$wt, clustered = cs$cl, extra_singletons = cs$extra)
    ctrl <- d92_ctrl(d, p = cs$p, chunk = 7L)
    info <- paste(cs$pattern, cs$cg, "p=", cs$p, "wt=", cs$wt, "cl=", cs$cl, "extra=", cs$extra)

    e_new <- d92_ns(wrap = d92_tracked)
    e_ref <- d92_no_carrier()

    r_new <- d92_run(e_new, d, ctrl, p = cs$p, cg = cs$cg, weighted = cs$wt, clustered = cs$cl)
    r_ref <- d92_run(e_ref, d, ctrl, p = cs$p, cg = cs$cg, weighted = cs$wt, clustered = cs$cl)

    d92_expect_parity(r_new, r_ref, info = info)
    if (cs$extra > 0L) {
      expect_gt(nrow(r_new$reporting$if_full), r_new$covariance_stream$n_a_units)
    }

    # Verify exactly 1 producer pass + 3 slice passes in the unchanged row partition.
    expect_identical(d92_n(e_new, "ipw_point_if_retained_reserve"), 1L, info = info)
    expect_identical(d92_n(e_new, "ipw_point_if_retained_fill"), 1L, info = info)
    expect_identical(d92_n(e_new, "replay_ipw_point_if_retained"), 3L, info = info)
    expect_identical(d92_n(e_new, "ipw_point_if_general_state_build"), 0L, info = info)

    ret_rec <- if (!is.null(r_new$resource_accounting$point_general)) {
      pg <- r_new$resource_accounting$point_general
      expect_identical(pg$replay_state$status, "superseded_by_retained_if", info = info)
      expect_false(isTRUE(pg$replay_state$admitted), info = info)
      n_chunks <- as.integer(r_new$resource_accounting$point_producer$n_chunks)
      L_used <- length(r_new$point_links_used)
      expect_identical(
        d92_n(e_ref, "query_ipw_view_chunk") - d92_n(e_new, "query_ipw_view_chunk"),
        3L * n_chunks * L_used,
        info = info
      )
      pg$retained_if
    } else {
      r_new$resource_accounting$point_batch$retained_if
    }
    expect_identical(ret_rec$status, "retained", info = info)
    expect_true(isTRUE(ret_rec$admitted) && isTRUE(ret_rec$released), info = info)
    expect_identical(ret_rec$producer_passes, 1L, info = info)
    expect_identical(ret_rec$slice_passes, 3L, info = info)
    expect_identical(
      as.integer(ret_rec$chunk_rows),
      as.integer(r_ref$resource_accounting$point_producer$effective_chunk_rows),
      info = info
    )

    # Neither the carrier nor any sealed plan survives into the returned object.
    leaked <- FALSE
    walk <- function(x, depth = 0L) {
      if (leaked || depth > 10L) return(invisible(NULL))
      if (is.environment(x) && identical(x$kind, "ipw_point_if_retained")) leaked <<- TRUE
      if (is.list(x)) {
        if ( identical(x$kind, "cdid_ipw_point_if_general_plan") ||
             identical(x$kind, "cdid_ipw_point_if_batch_plan") ||
             identical(x$kind, "cdid_custom_ipw_if_plan") ) leaked <<- TRUE
        for (v in x) walk(v, depth + 1L)
      }
    }
    walk(r_new)
    expect_false(leaked, info = info)
  }
})

# ---------------------------------------------------------------------------
# 2. Direct general carrier pass vs standalone replay_ipw_point_if_view and
#    4-term component sums; numerical_transfer declines.
# ---------------------------------------------------------------------------
test_that("direct general retained IF carrier matches standalone replay and 4-term component sums bitwise", {
  d <- d92_panel(pattern = "half", multi_cohort = TRUE, weighted = TRUE, extra_singletons = 2L)
  ctrl <- d92_ctrl(d, p = 2L, chunk = 7L)
  res <- d92_run(d92_no_carrier(), d, ctrl, p = 2L, cg = "notyettreated", weighted = TRUE)
  links <- res$point_links_used
  pc <- ctrl$point_control
  oc <- ctrl$owner_control

  plan <- api$ipw_point_if_general_plan(links, pc)
  expect_null(api$ipw_point_if_retained_reserve(
    plan, pc, live = 0, work = 0, envelope_control = oc,
    n_total = nrow(res$reporting$if_full), clustered = TRUE, numerical_transfer = TRUE
  ))

  carrier <- api$ipw_point_if_retained_reserve(
    plan, pc, live = 0, work = 0, envelope_control = oc,
    n_total = nrow(res$reporting$if_full), clustered = TRUE, numerical_transfer = FALSE
  )
  expect_true(is.environment(carrier))
  expect_identical(carrier$state, "reserved")

  fill_batches <- list()
  api$ipw_point_if_retained_fill(carrier, plan, pc, callback = function(b) {
    fill_batches[[length(fill_batches) + 1L]] <<- b
  })
  expect_identical(carrier$state, "filled")
  expect_identical(carrier$producer_passes, 1L)

  slice_batches <- list()
  s_info <- api$replay_ipw_point_if_retained(carrier, plan, pc, callback = function(b) {
    slice_batches[[length(slice_batches) + 1L]] <<- b
  })
  expect_identical(s_info$served_from, "retained_point_if")
  expect_identical(carrier$slice_passes, 1L)

  std_batches <- list()
  api$replay_ipw_point_if_view(links, pc, callback_schema = TRUE, callback = function(b) {
    std_batches[[length(std_batches) + 1L]] <<- b
  })
  comp_batches <- list()
  api$replay_ipw_point_if_view(links, pc, include_components = TRUE, callback = function(b) {
    comp_batches[[length(comp_batches) + 1L]] <<- b
  })

  expect_identical(length(fill_batches), length(std_batches))
  expect_identical(length(slice_batches), length(std_batches))
  for (i in seq_along(std_batches)) {
    expect_identical(fill_batches[[i]], std_batches[[i]])
    expect_identical(slice_batches[[i]], std_batches[[i]])
    comp <- comp_batches[[i]]$components
    expect_identical(slice_batches[[i]]$phi, ((comp$treated + comp$control) + comp$beta_score) + comp$rate)
  }

  # Full covariance stream with a fresh reserved carrier matches standalone covariance bitwise.
  carrier_cov <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
  cg <- api$stream_ipw_point_covariance_general(plan, pc, carrier = carrier_cov)
  cv <- api$stream_ipw_point_covariance_view(links, pc)
  for (nm in c("covariance", "mean_phi", "estimates", "n_a_units")) {
    expect_identical(cg[[nm]], cv[[nm]])
  }
  expect_identical(carrier_cov$producer_passes, 1L)
  expect_identical(carrier_cov$slice_passes, 1L)
  rel <- api$ipw_point_if_retained_release(carrier_cov, plan)
  expect_identical(rel$status, "retained")
  expect_identical(carrier_cov$state, "released")
  expect_false(any(c("phi", "raw_ref", "seal", "plan_snapshot", "plan_ref", "frame_ref", "panel_ref") %in%
    ls(carrier_cov, all.names = TRUE)))
  api$ipw_point_if_retained_release(carrier, plan)
})

# ---------------------------------------------------------------------------
# 3. Static custom IPW plan adapter + retained IF carrier: user callback count
#    and order unchanged; frozen output mutations refused.
# ---------------------------------------------------------------------------
test_that("static custom IPW plan and retained IF carrier preserve user callback count/order and refuse frozen output mutation", {
  d <- d92_panel(pattern = "half", multi_cohort = TRUE, weighted = TRUE, clustered = TRUE, extra_singletons = 2L)
  ctrl <- d92_ctrl(d, p = 2L, chunk = 7L)

  mk_custom <- function(trace_vec) {
    force(trace_vec)
    function(y1, y0, D, covariates, i.weights, inffunc = TRUE) {
      dy <- y1 - y0
      w1 <- i.weights * D
      w0 <- i.weights * (1 - D)
      mu1 <- sum(w1 * dy) / sum(w1)
      mu0 <- sum(w0 * dy) / sum(w0)
      att <- mu1 - mu0
      n_obs <- length(dy)
      psi <- n_obs * (w1 * (dy - mu1) / sum(w1) - w0 * (dy - mu0) / sum(w0))
      trace_vec$calls <- trace_vec$calls + 1L
      trace_vec$log[[trace_vec$calls]] <- c(n_obs, sum(y1), sum(y0), sum(D), sum(i.weights))
      list(ATT = att, att.inf.func = psi)
    }
  }

  tr_new <- new.env(parent = emptyenv()); tr_new$calls <- 0L; tr_new$log <- list()
  tr_ref <- new.env(parent = emptyenv()); tr_ref$calls <- 0L; tr_ref$log <- list()

  e_new <- d92_ns(wrap = d92_tracked)
  # Reference forces the standalone custom view consumers without carrier or plan.
  e_ref <- d92_no_carrier(extra_edits = list(
    stream_custom_ipw_covariance_plan = function(plan, control, ...) {
      api$stream_custom_ipw_covariance_view(point_links = plan$point_links, control = control)
    },
    replay_custom_ipw_if_plan = function(plan, control, callback, callback_schema = TRUE, ...) {
      api$replay_custom_ipw_if_view(point_links = plan$point_links, control = control,
        callback = callback, callback_schema = callback_schema)
    }
  ))

  # The original decimal-weight panel has the accepted producer's strict
  # observer-mass refusal before either cache route can run. Preserve it.
  decimal_new <- tryCatch(
    d92_run(e_new, d, ctrl, p = 2L, cg = "notyettreated", weighted = TRUE, clustered = TRUE,
            custom_est_method = mk_custom(tr_new)), error = identity)
  decimal_ref <- tryCatch(
    d92_run(e_ref, d, ctrl, p = 2L, cg = "notyettreated", weighted = TRUE, clustered = TRUE,
            custom_est_method = mk_custom(tr_ref)), error = identity)
  expect_s3_class(decimal_new, "cdid_ipw_point_view_integrity_error")
  expect_s3_class(decimal_ref, "cdid_ipw_point_view_integrity_error")
  expect_identical(class(decimal_new), class(decimal_ref))
  expect_identical(conditionMessage(decimal_new),
    "Reconstructed observer-buffer pair masses or rates disagree with `link_view$support`.")
  expect_identical(conditionMessage(decimal_new), conditionMessage(decimal_ref))
  expect_true(is.integer(tr_new$calls) && length(tr_new$calls) == 1L &&
    !is.na(tr_new$calls) && tr_new$calls >= 1L && tr_new$calls <= 3L)
  expect_identical(length(tr_new$log), tr_new$calls)
  expect_identical(tr_new$calls, tr_ref$calls)
  expect_identical(tr_new$log, tr_ref$log)

  # Dyadic weights exercise successful weighted custom/F > A/cluster parity
  # without changing the panel helper, callback formulas or producer contract.
  d$w <- 1 + (d$id %% 5L) / 4
  d92_reset(e_new); d92_reset(e_ref)
  tr_new$calls <- 0L; tr_new$log <- list()
  tr_ref$calls <- 0L; tr_ref$log <- list()

  r_new <- d92_run(e_new, d, ctrl, p = 2L, cg = "notyettreated", weighted = TRUE, clustered = TRUE,
                   custom_est_method = mk_custom(tr_new))
  r_ref <- d92_run(e_ref, d, ctrl, p = 2L, cg = "notyettreated", weighted = TRUE, clustered = TRUE,
                   custom_est_method = mk_custom(tr_ref))

  expect_gt(tr_new$calls, 0L)
  expect_identical(tr_new$calls, tr_ref$calls)
  expect_identical(tr_new$log, tr_ref$log)
  d92_expect_parity(r_new, r_ref, info = "custom IPW carrier vs standalone")

  pc_rec <- r_new$resource_accounting$point_custom
  expect_false(is.null(pc_rec))
  expect_identical(pc_rec$retained_if$status, "retained")
  expect_identical(pc_rec$retained_if$producer_passes, 1L)
  expect_identical(pc_rec$retained_if$slice_passes, 3L)
  expect_identical(d92_n(e_new, "ipw_point_if_retained_fill"), 1L)
  expect_identical(d92_n(e_new, "replay_ipw_point_if_retained"), 3L)

  # Direct custom plan verification refuses mutated frozen callback outputs.
  c_links <- r_new$point_links_used
  c_plan <- api$custom_ipw_if_plan(c_links, ctrl$point_control)
  expect_true(api$custom_ipw_if_plan_verify(c_plan, ctrl$point_control))
  # The adapter preserves the standalone covariance return schema and converts
  # only the original cohort arrays, once, before replay.
  c_cov <- api$stream_custom_ipw_covariance_plan(c_plan, ctrl$point_control)
  c_ref <- api$stream_custom_ipw_covariance_view(c_links, ctrl$point_control)
  expect_identical(c_cov, c_ref)
  for (ell in seq_along(c_links)) {
    expect_identical(c_plan$norm$frozen[[ell]]$b_rate, as.double(c_links[[ell]]$b_rate))
    expect_identical(c_plan$norm$frozen[[ell]]$q, as.double(c_links[[ell]]$masses$q))
    expect_identical(c_plan$norm$frozen[[ell]]$rates, as.double(c_links[[ell]]$masses$rates))
  }
  # Equal values in a copied link do not preserve the sealed binding identity.
  pl <- c_plan$point_links[[1L]]
  pl_copy <- pl
  pl_copy$estimate <- pl$estimate
  expect_identical(pl_copy, pl)
  expect_false(api$gmm_sexp_aliases(pl_copy, pl))
  restore <- d92_rebind(pl$ref, "custom_ipw_link", pl_copy)
  expect_error(api$custom_ipw_if_plan_verify(c_plan, ctrl$point_control),
    class = "cdid_custom_ipw_error")
  restore()

  seen <- 0L
  restore <- NULL
  expect_error(api$replay_custom_ipw_if_plan(c_plan, ctrl$point_control, callback = function(b) {
    seen <<- seen + 1L
    if (seen == c_plan$preflight$n_chunks) {
      restore <<- d92_rebind(pl$ref, "custom_ipw_link", pl_copy)
    }
  }), class = "cdid_custom_ipw_error")
  expect_identical(seen, as.integer(c_plan$preflight$n_chunks))
  if (!is.null(restore)) restore()
  bad_plan <- c_plan
  bad_plan$ref <- NULL
  bad_plan$point_links[[1L]]$psi_obs <- bad_plan$point_links[[1L]]$psi_obs + 1e-6
  bad_plan <- api$ipw_owner_seal(bad_plan, "custom_ipw_if_plan")
  expect_error(
    api$replay_custom_ipw_if_plan(bad_plan, ctrl$point_control, callback = function(b) NULL),
    class = "cdid_custom_ipw_error"
  )
})

# ---------------------------------------------------------------------------
# 4. Admission boundaries: typed capacity decline falls back cleanly while
#    malformed inputs, exact-count overflows, and final-callback mutations propagate.
# ---------------------------------------------------------------------------
test_that("typed capacity bounds decline cleanly while malformed inputs, overflows and final-callback mutations propagate", {
  d <- d92_panel(pattern = "half", multi_cohort = TRUE, weighted = TRUE)
  ctrl <- d92_ctrl(d, p = 1L, chunk = 7L)
  res <- d92_run(d92_no_carrier(), d, ctrl, p = 1L, cg = "notyettreated", weighted = TRUE)
  links <- res$point_links_used
  pc <- ctrl$point_control
  oc <- ctrl$owner_control
  plan <- api$ipw_point_if_general_plan(links, pc)
  N <- as.double(plan$n_a_units)
  L <- as.double(plan$n_links)
  B <- as.double(plan$step)
  inventory <- api$ipw_point_if_retained_inventory(plan)
  need_entries <- inventory$base_live + inventory$entries
  need_work <- inventory$base_work + inventory$work
  expect_gt(inventory$entries, api$ipw_point_if_retained_entries(N, L, B))
  expect_gt(inventory$work, api$ipw_point_if_retained_work(N, L))

  # Exact boundary admits; 1 below on retained entries or work ops declines to NULL.
  oc_exact <- oc; oc_exact$max_retained_entries <- need_entries; oc_exact$max_work_ops <- need_work
  c_ok <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc_exact)
  expect_true(is.environment(c_ok))

  oc_low_e <- oc_exact; oc_low_e$max_retained_entries <- need_entries - 1
  expect_null(api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc_low_e))

  oc_low_w <- oc_exact; oc_low_w$max_work_ops <- need_work - 1
  expect_null(api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc_low_w))

  # Pipeline fallback when carrier declines: falls back to sparse state (or ordinary replay) with identical numerics.
  e_decline_raw <- d92_ns(
    edits = list(ipw_point_if_retained_reserve = function(...) NULL),
    wrap = d92_tracked
  )
  r_fallback <- d92_run(e_decline_raw, d, ctrl, p = 1L, cg = "notyettreated", weighted = TRUE)
  d92_expect_parity(r_fallback, res, info = "carrier declined -> sparse state fallback")
  expect_identical(r_fallback$resource_accounting$point_general$retained_if$status, "streamed")
  expect_true(isTRUE(r_fallback$resource_accounting$point_general$replay_state$admitted))

  # Malformed inputs and IEEE 754 exact-integer overflows (2^53 - 1) MUST propagate, never return NULL.
  expect_error(
    api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc, n_total = plan$n_a_units - 1L),
    class = "cdid_ipw_count_error"
  )
  expect_error(
    api$ipw_point_if_retained_reserve(plan, pc, live = -1, work = 0, envelope_control = oc),
    class = "cdid_ipw_point_view_input_error"
  )
  expect_error(
    api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc, clustered = NA),
    class = "cdid_ipw_point_view_input_error"
  )
  expect_error(
    api$ipw_point_if_retained_reserve(plan, pc, live = 2^53 - 2, work = 0, envelope_control = oc),
    class = "cdid_ipw_owner_bound_error"
  )
  expect_error(
    api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 2^53 - 2, envelope_control = oc),
    class = "cdid_ipw_owner_bound_error"
  )

  # Final-callback mutation refusal: mutating an owner/link binding OR the locked
  # raw carrier matrix during the LAST chunk callback is refused immediately.
  n_chunks <- as.integer(plan$preflight$n_chunks)
  expect_gt(n_chunks, 1L)
  pl1 <- plan$point_links[[1L]]
  bump <- function(x) { x$estimate <- x$estimate + 1; x }

  c_mut_fill <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
  seen <- 0L; restore <- NULL
  expect_error(
    api$ipw_point_if_retained_fill(c_mut_fill, plan, pc, callback = function(b) {
      seen <<- seen + 1L
      if (seen == n_chunks) restore <<- d92_rebind(pl1$ref, "point_link_view", bump(pl1))
    }),
    class = "cdid_ipw_point_view_general_error"
  )
  expect_identical(seen, n_chunks)
  if (!is.null(restore)) restore()
  expect_identical(c_mut_fill$state, "released")
  expect_false(exists("plan_ref", c_mut_fill, inherits = FALSE))

  # Fill a clean carrier, then mutate on the LAST callback of a slice pass:
  c_slice <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
  api$ipw_point_if_retained_fill(c_slice, plan, pc, callback = function(b) NULL)

  seen <- 0L; restore <- NULL
  expect_error(
    api$replay_ipw_point_if_retained(c_slice, plan, pc, callback = function(b) {
      seen <<- seen + 1L
      if (seen == n_chunks) restore <<- d92_rebind(pl1$ref, "point_link_view", bump(pl1))
    }),
    class = "cdid_ipw_point_view_general_error"
  )
  expect_identical(seen, n_chunks)
  if (!is.null(restore)) restore()
  expect_identical(c_slice$state, "released")
  expect_false(exists("phi", c_slice, inherits = FALSE))

  # Tampering with the carrier's locked raw matrix on the LAST slice callback is refused:
  c_raw <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
  api$ipw_point_if_retained_fill(c_raw, plan, pc, callback = function(b) NULL)
  seen <- 0L
  expect_error(
    api$replay_ipw_point_if_retained(c_raw, plan, pc, callback = function(b) {
      seen <<- seen + 1L
      if (seen == n_chunks) {
        unlockBinding("phi", c_raw)
        c_raw$phi[1L, 1L] <- c_raw$phi[1L, 1L] + 1
        lockBinding("phi", c_raw)
      }
    }),
    class = "cdid_ipw_point_view_batch_error"
  )
  expect_identical(seen, n_chunks)
  expect_identical(c_raw$state, "released")

  # State tampering used to skip the raw identity check on the final callback.
  c_state <- api$ipw_point_if_retained_reserve(plan, pc, live = 0, work = 0, envelope_control = oc)
  api$ipw_point_if_retained_fill(c_state, plan, pc, callback = function(b) NULL)
  seen <- 0L
  expect_error(
    api$replay_ipw_point_if_retained(c_state, plan, pc, callback = function(b) {
      seen <<- seen + 1L
      if (seen == n_chunks) c_state$state <- "reserved"
    }),
    class = "cdid_ipw_point_view_batch_error"
  )
  expect_identical(seen, n_chunks)
  expect_identical(c_state$state, "released")
  expect_false(exists("raw_ref", c_state, inherits = FALSE))
  api$ipw_point_if_retained_release(c_ok, plan)
})
