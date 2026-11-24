# SPDX-License-Identifier: GPL-3.0-only
#
# Covariance planning near the row cap. The certified full-width decision must
# return exactly the plan of the complete (legacy plus full-block) plan, the
# per-block-pair envelope store must not change any plan. Backends are
# captured from small public fits; the row cap and the retained charge are then
# varied so that the planner's arithmetic binds, refuses or leaves the cap.
library(testthat)
api <- asNamespace("cdid")

cert_panel <- function(N = 80L, T = 6L, seed = 5L, miss = 0.25) {
  set.seed(seed)
  g <- sample(c(4, 6, 0), N, replace = TRUE, prob = c(0.3, 0.3, 0.4))
  g[1:3] <- c(4, 6, 0)
  d <- expand.grid(t = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "t")]
  d$G <- g[d$id]
  x <- stats::rnorm(N)
  d$x <- x[d$id]
  d$xv <- round(d$x + 0.1 * d$t + stats::rnorm(nrow(d), sd = 0.3), 3)
  d$y <- 0.3 * d$x + 0.2 * d$t + ifelse(d$G > 0 & d$t >= d$G, 1, 0) + stats::rnorm(nrow(d))
  if (miss > 0) d <- d[!(stats::runif(nrow(d)) < miss & d$t > 1L), ]
  d
}

# Run a fit and keep the arguments of the first covariance plan of the pipeline.
cert_capture <- function(args, data = cert_panel()) {
  box <- new.env()
  assign(".cert_box", box, envir = globalenv())
  on.exit(rm(".cert_box", envir = globalenv()), add = TRUE)
  suppressMessages(trace("cp04_stream_plan", where = api, print = FALSE, tracer = quote(
    if (is.null(.cert_box$backend)) { .cert_box$backend <- backend; .cert_box$links <- links; .cert_box$retained <- retained })))
  on.exit(suppressMessages(untrace("cp04_stream_plan", where = api)), add = TRUE)
  base <- list(yname = "y", tname = "t", idname = "id", gname = "G", data = data, allow_unbalanced_panel = TRUE,
    mode = "chained", control_group = "nevertreated", base_period = "universal", bstrap = FALSE, cband = FALSE,
    gmm_weighting = "within")
  # A refusal by a later check (for example the diagnostic copy) still leaves
  # the captured plan arguments.
  fit <- tryCatch(do.call(cdid::att_gt, utils::modifyList(base, args)), cdid_error = function(e) e)
  if (is.null(box$backend)) stop("no covariance plan was reached: ", conditionMessage(fit))
  list(fit = fit, backend = box$backend, links = box$links, retained = box$retained)
}

# The complete plan without the certificate (the 811ed24 decision rule).
cert_complete <- function(backend, links, retained) {
  api$cp04_plan_memo_clear()
  legacy <- api$cp04_stream_plan_legacy(backend, links, retained)
  if (length(legacy$blocks) <= 1L) return(legacy)
  full <- api$cp04_stream_plan_full_block(backend, links, retained, legacy)
  if (is.null(full)) legacy else full
}

cert_try <- function(expr) tryCatch(expr, error = function(e) structure(list(class = class(e), message = conditionMessage(e)), class = "cert_error"))

# Capture the plan arguments at several workspace sizes (row caps 5, 8 and 11
# on this panel), then sweep the retained charge up to and beyond the point
# where the planner's arithmetic binds and refuses. Every certified plan must
# equal the complete plan, and the public plan must equal the complete plan or
# refuse identically.
cert_budgets <- c(8, 12, 16) * 2^20
cert_captures <- function(args, budgets = cert_budgets)
  lapply(budgets, function(bytes) cert_capture(c(args, list(control = list(numeric_workspace_bytes = bytes)))))
cert_sweep <- function(captures) {
  out <- list(certified = 0L, binding = 0L, points = 0L)
  for (cf in captures) {
    b <- cf$backend; links <- cf$links; L <- length(links)
    if (L <= b$control$max_chunk_rows) next
    room <- b$control$max_matrix_entries - b$accounting$owned_entries - b$accounting$borrowed_entries - 4 * L * L
    rs <- unique(floor(pmax(0, c(cf$retained, room * c(0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999, 1.001)))))
    for (r in rs) {
      out$points <- out$points + 1L
      complete <- cert_try(cert_complete(b, links, r))
      api$cp04_plan_memo_clear()
      certified <- cert_try(api$cp04_stream_plan_certified(b, links, r))
      if (!inherits(certified, "cert_error") && !is.null(certified)) {
        out$certified <- out$certified + 1L
        expect_identical(certified, complete)
      }
      api$cp04_plan_memo_clear()
      public <- cert_try(api$cp04_stream_plan(b, links, r))
      expect_identical(public, complete)
      # Property behind the certificate: no block pair admits fewer rows than all links.
      B_cap <- min(api$cp04_backend_domain(b)$N, b$control$max_chunk_rows)
      B_all <- cert_try(api$cp04_plan_rows(b, list(links), r + 4 * L * L))
      api$cp04_plan_memo_clear()
      legacy <- cert_try(api$cp04_stream_plan_legacy(b, links, r))
      if (!inherits(legacy, "cert_error")) {
        if (legacy$chunk_rows < B_cap) out$binding <- out$binding + 1L
        if (!inherits(B_all, "cert_error") && isTRUE(api$cp04_stream_plan_monotone_envelope(b)))
          expect_gte(legacy$chunk_rows, min(B_cap, B_all))
      }
    }
  }
  api$cp04_plan_memo_clear()
  out
}

test_that("certified plans equal the complete plan: separate DR (native envelope)", {
  cfs <- cert_captures(list(est_method = "dr", outcome_pooling = "separate", xformla = ~x))
  expect_true(api$cp04_stream_plan_monotone_envelope(cfs[[1L]]$backend))
  s <- cert_sweep(cfs)
  expect_gt(s$certified, 0L)
  expect_gt(s$binding, 0L)
})

test_that("certified plans equal the complete plan: separate DR with a time-varying covariate", {
  cfs <- cert_captures(list(est_method = "dr", outcome_pooling = "separate", xformla = ~xv), c(12, 16, 24) * 2^20)
  s <- cert_sweep(cfs)
  expect_gt(s$certified, 0L)
})

test_that("certified plans equal the complete plan: separate imputation", {
  cfs <- cert_captures(list(est_method = "dr", outcome_pooling = "separate", xformla = ~1,
    selection = "conditional_mtar", selection_est_method = "imputation", selection_xformla = ~1))
  expect_true(api$cp04_stream_plan_monotone_envelope(cfs[[1L]]$backend))
  s <- cert_sweep(cfs)
  expect_gt(s$certified, 0L)
  expect_gt(s$binding, 0L)
})

test_that("certified plans equal the complete plan: separate augmented (selection envelope)", {
  cfs <- cert_captures(list(est_method = "reg", outcome_pooling = "separate", xformla = ~1,
    selection = "conditional_mtar", selection_est_method = "dr", selection_xformla = ~1), c(12, 16, 24) * 2^20)
  expect_true(inherits(cfs[[1L]]$backend, "cdid_augmented_backend"))
  expect_true(api$cp04_stream_plan_monotone_envelope(cfs[[1L]]$backend))
  s <- cert_sweep(cfs)
  expect_gt(s$certified, 0L)
  expect_gt(s$binding, 0L)
})

test_that("backends outside the monotone list never use the certificate", {
  for (args in list(
    list(est_method = "reg", outcome_pooling = "pooled", xformla = ~1,
      selection = "conditional_mtar", selection_est_method = "imputation", selection_xformla = ~1),
    list(est_method = "dr", outcome_pooling = "pooled", xformla = ~x),
    list(est_method = "dr", outcome_pooling = "separate", xformla = ~1,
      selection = "conditional_mtar", selection_est_method = "ipw", selection_xformla = ~1))) {
    cf <- cert_capture(c(args, list(control = list(numeric_workspace_bytes = 16 * 2^20))))
    expect_false(api$cp04_stream_plan_monotone_envelope(cf$backend))
    b <- cf$backend
    expect_null(api$cp04_stream_plan_certified(b, cf$links, cf$retained))
    expect_identical(cert_try(api$cp04_stream_plan(b, cf$links, cf$retained)), cert_try(cert_complete(b, cf$links, cf$retained)))
  }
})

test_that("within one pipeline scope the per-block-pair store changes no plan and is bounded", {
  cf <- cert_capture(list(est_method = "dr", outcome_pooling = "separate", xformla = ~1,
    selection = "conditional_mtar", selection_est_method = "imputation", selection_xformla = ~1,
    control = list(numeric_workspace_bytes = 8 * 2^20)))
  b <- cf$backend
  links <- cf$links; L <- length(links); K <- ceiling(L / b$control$max_chunk_rows)
  expect_gt(K, 1)
  room <- b$control$max_matrix_entries - b$accounting$owned_entries - b$accounting$borrowed_entries - 4 * L * L
  rs <- floor(c(0, 0.9, 0.99) * room)
  fresh <- lapply(rs, function(r) cert_try(api$cp04_stream_plan_legacy(b, links, r)))
  prev <- api$cp04_plan_scope_enter()
  reused <- lapply(rs, function(r) cert_try(api$cp04_stream_plan_legacy(b, links, r)))
  m <- api$cp04_plan_scope$current$memo
  paths <- api$cp04_plan_scope$current$paths
  api$cp04_plan_scope_exit(prev)
  expect_identical(reused, fresh)
  expect_false(is.null(m))
  for (nm in c("one_output", "one_scratch", "two_output", "two_scratch", "rows", "work", "peak"))
    expect_length(m[[nm]], K * (K + 1) / 2)
  expect_identical(unname(paths[["store_pair_misses"]]), K * (K + 1) / 2)
  # At least the second build reads every pair from the store (a refused build stops early).
  expect_gte(unname(paths[["store_pair_hits"]]), K * (K + 1) / 2)
  expect_lte(unname(paths[["store_pair_hits"]]), (length(rs) - 1) * K * (K + 1) / 2)
  # A different link vector or backend starts a new store.
  prev <- api$cp04_plan_scope_enter()
  api$cp04_stream_plan_legacy(b, links, 0)
  sub <- links[-1L]
  api$cp04_stream_plan_legacy(b, sub, 0)
  expect_identical(api$cp04_plan_scope$current$memo$key[[2L]], as.integer(sub))
  b2 <- cert_capture(list(est_method = "dr", outcome_pooling = "separate", xformla = ~1,
    selection = "conditional_mtar", selection_est_method = "imputation", selection_xformla = ~1,
    control = list(numeric_workspace_bytes = 12 * 2^20)))$backend
  in_scope <- api$cp04_stream_plan_legacy(b2, links, 0)
  api$cp04_plan_scope_exit(prev)
  expect_identical(in_scope, api$cp04_stream_plan_legacy(b2, links, 0))
  expect_null(api$cp04_plan_scope$current)
})

test_that("a plan built outside a pipeline keeps nothing, so a changed envelope is always evaluated", {
  cf <- cert_capture(list(est_method = "dr", outcome_pooling = "separate", xformla = ~x,
    control = list(numeric_workspace_bytes = 12 * 2^20)))
  b <- cf$backend; links <- cf$links
  expect_null(api$cp04_plan_scope$current)
  real <- api$cp04_stream_plan_legacy(b, links, 0)
  expect_gt(length(real$blocks), 1L)
  # The same planner in an environment whose producer budget charges twice the
  # work: the plan must follow the changed budget, not a stored value.
  e <- new.env(parent = api)
  for (nm in c("cp04_stream_plan_legacy", "cp04_plan_rows")) { f <- get(nm, api); environment(f) <- e; assign(nm, f, e) }
  e$chained_full_if_budget <- function(backend, B, links) {
    out <- api$chained_full_if_budget(backend, B, links)
    out$work <- 2 * out$work
    out
  }
  inflated <- cert_try(e$cp04_stream_plan_legacy(b, links, 0))
  e2 <- new.env(parent = api)
  for (nm in c("cp04_stream_plan_legacy", "cp04_plan_rows")) { f <- get(nm, api); environment(f) <- e2; assign(nm, f, e2) }
  e2$chained_full_if_budget <- e$chained_full_if_budget
  api$cp04_plan_memo_clear()
  expect_identical(inflated, cert_try(e2$cp04_stream_plan_legacy(b, links, 0)))
  expect_false(identical(inflated, real))
})

test_that("a fit leaves no plan scope behind and records its planning paths", {
  assign("last", NULL, envir = api$cp04_plan_path_record)
  cf <- cert_capture(list(est_method = "dr", outcome_pooling = "separate", xformla = ~x,
    control = list(numeric_workspace_bytes = 12 * 2^20)))
  expect_null(api$cp04_plan_scope$current)
  p <- api$cp04_plan_path_record$last
  expect_identical(names(p), api$cp04_plan_path_names)
  expect_gt(p[["certified"]], 0)
  expect_identical(p[["complete_block_pairs"]], 0)
})
test_that("the certificate's work bounds enclose the complete plan's work", {
  for (args in list(list(est_method = "dr", outcome_pooling = "separate", xformla = ~xv),
    list(est_method = "dr", outcome_pooling = "separate", xformla = ~1,
      selection = "conditional_mtar", selection_est_method = "imputation", selection_xformla = ~1))) {
    for (cf in cert_captures(args, c(12, 16) * 2^20)) {
      b <- cf$backend; links <- cf$links
      seen <- new.env()
      cert <- api$cp04_stream_plan_certified
      e <- new.env(parent = api)
      e$cp04_stream_plan_full_block <- function(backend, links, retained, legacy) { seen$lower <- legacy$work; NULL }
      environment(cert) <- e
      expect_null(cert(b, links, cf$retained))
      api$cp04_plan_memo_clear()
      legacy <- api$cp04_stream_plan_legacy(b, links, cf$retained)
      expect_false(is.null(seen$lower))
      expect_lte(seen$lower, legacy$work)
      L <- length(links); sz <- lengths(legacy$blocks); K <- length(sz)
      gram <- 0
      for (a in seq_len(K)) for (bb in a:K) gram <- gram + 4 * api$cp04_backend_domain(b)$N * (if (a == bb) sz[a] else sz[a] + sz[bb])^2
      expect_gt(seen$lower, gram)
    }
  }
  api$cp04_plan_memo_clear()
})
