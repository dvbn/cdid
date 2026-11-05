zt_ns <- function(name) getFromNamespace(name, "cdid")
zt_fixture <- function() {
  cells <- lapply(1:8, function(i) list(cell_id = i, marker = i / 7))
  groups <- lapply(1:8, function(i) list(cell_ids = c(i), p = as.integer((i - 1L) %% 3L), marker = i))
  templates <- lapply(groups, function(g) list(lambda = as.double(g$cell_ids), nu = numeric(g$p),
    group = g, cells = cells[g$cell_ids], N = 48L))
  list(templates = templates, groups = groups, cells = cells, active = integer(), N = 48L)
}
zt_native <- function(a, cursor = 1L) zt_ns("reconstruction_zero_run_native_kernel")()(a$templates, a$groups, a$cells, a$active, a$N, cursor)
zt_reference <- function(a, cursor = 1L) {
  ids <- integer(); lambda <- numeric(); gids <- integer(); nu <- list(); next_gid <- cursor
  for (gid in seq.int(cursor, length(a$groups))) {
    g <- a$groups[[gid]]
    x <- zt_ns("reconstruction_zero_template")(a$templates, gid, gid %in% a$active, g, a$cells, a$N)
    if (is.null(x)) break
    ids <- c(ids, g$cell_ids); lambda <- c(lambda, x$lambda); gids <- c(gids, gid)
    nu <- c(nu, list(x$nu)); next_gid <- gid + 1L
  }
  if (!length(gids)) NULL else list(cell_ids = ids, lambda = lambda, group_ids = gids, nu = nu, next_gid = next_gid)
}

testthat::test_that("native inactive runs match the original three identical checks and numeric copies", {
  for (damage in c("none", "active", "miss", "group", "cells", "N", "closure")) {
    a <- zt_fixture()
    if (damage == "active") a$active <- 3L
    if (damage == "miss") a$templates[3L] <- list(NULL)
    if (damage == "group") a$groups[[3L]]$marker <- -1
    if (damage == "cells") a$cells[[3L]]$marker <- -1
    if (damage == "N") a$templates[[3L]]$N <- 49L
    if (damage == "closure") {
      f <- function(x) x
      environment(f) <- new.env(parent = baseenv())
      a$groups[[3L]]$fn <- f; a$templates[[3L]]$group <- a$groups[[3L]]
      environment(f) <- new.env(parent = baseenv()); a$groups[[3L]]$fn <- f
    }
    testthat::expect_identical(zt_native(a), zt_reference(a), info = damage)
  }
  for (N in list(0, -0, NA_real_, NaN, -NaN)) {
    a <- zt_fixture(); a$N <- N
    a$templates <- lapply(a$templates, function(x) { x$N <- N; x })
    a$templates[[1L]]$lambda <- -0
    a$templates[[2L]]$lambda <- NA_real_; a$templates[[3L]]$lambda <- NaN
    a$templates[[2L]]$nu[] <- -0; a$templates[[3L]]$nu[] <- c(NA_real_, NaN)
    got <- zt_native(a)
    testthat::expect_identical(got, zt_reference(a))
    M <- length(got$cell_ids); H <- length(got$group_ids); P <- sum(lengths(got$nu))
    testthat::expect_equal(unname(zt_ns("cp04_payload")(got)["entries"]), 93 + 2 * M + 10 * H + P)
    for (i in seq_along(got$nu)) if (length(got$nu[[i]]))
      testthat::expect_false(zt_ns("gmm_sexp_aliases")(got$nu[[i]], a$templates[[i]]$nu))
    old <- a$templates
    got$lambda[] <- 100; got$nu[[2L]][] <- 100
    testthat::expect_identical(a$templates, old)
  }
  # Equal copied metadata remains eligible; pointer identity does not stand
  # in for any of the three semantic identical() checks.
  a <- zt_fixture()
  a$groups[[3L]] <- unserialize(serialize(a$groups[[3L]], NULL))
  a$cells[[3L]] <- unserialize(serialize(a$cells[[3L]], NULL))
  testthat::expect_identical(zt_native(a), zt_reference(a))
  a <- zt_fixture(); a$N <- -0
  a$templates <- lapply(a$templates, function(x) { x$N <- 0; x })
  testthat::expect_identical(zt_native(a), zt_reference(a))
  a$N <- NaN; a$templates <- lapply(a$templates, function(x) { x$N <- -NaN; x })
  testthat::expect_identical(zt_native(a), zt_reference(a))
  a$templates[[3L]]$N <- NA_real_
  testthat::expect_identical(zt_native(a), zt_reference(a))
})

testthat::test_that("unsupported containers and output attributes stop before their group", {
  for (damage in c("named_lambda", "named_nu", "class", "empty", "bad_id", "named_ids", "double_ids")) {
    a <- zt_fixture()
    if (damage == "named_lambda") names(a$templates[[3L]]$lambda) <- "named"
    if (damage == "named_nu") names(a$templates[[3L]]$nu) <- c("a", "b")
    if (damage == "class") class(a$groups[[3L]]) <- "private_group"
    if (damage == "empty") a$groups[[3L]]$cell_ids <- integer()
    if (damage == "bad_id") a$groups[[3L]]$cell_ids <- 0L
    if (damage == "named_ids") names(a$groups[[3L]]$cell_ids) <- "named"
    if (damage == "double_ids") a$groups[[3L]]$cell_ids <- 3
    got <- zt_native(a)
    testthat::expect_identical(got$group_ids, 1:2)
    testthat::expect_identical(got$next_gid, 3L)
  }
  for (field in c("templates", "groups", "cells", "active")) {
    a <- zt_fixture(); attr(a[[field]], "private") <- TRUE
    testthat::expect_null(zt_native(a))
  }
  a <- zt_fixture(); names(a$cells) <- letters[1:8]
  testthat::expect_null(zt_native(a))
  a <- zt_fixture(); a$active <- NA_integer_
  testthat::expect_null(zt_native(a))
  # Compact integer ALTREP indices retain the scalar route, including after
  # their ordinary R subset has materialized their data.
  a <- zt_fixture(); a$groups[[1L]]$cell_ids <- 1:8
  a$templates[[1L]]$group <- a$groups[[1L]]
  a$templates[[1L]]$cells <- a$cells[1:8]; a$templates[[1L]]$lambda <- as.double(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L))
  testthat::expect_null(zt_native(a))
  a <- zt_fixture(); a$groups[[3L]]$private_altrep <- 1:100
  a$templates[[3L]]$group <- a$groups[[3L]]
  a$groups[[3L]]$marker <- -1
  testthat::expect_identical(zt_native(a)$group_ids, 1:2)
})

testthat::test_that("fresh native/scatter construction fits the unchanged forcing allowances", {
  a <- zt_fixture(); run <- zt_native(a); E <- function(x) unname(zt_ns("cp04_payload")(x)["entries"])
  nc <- length(a$cells); ng <- length(a$groups); pg <- vapply(a$groups, `[[`, integer(1), "p")
  # Borrowed source groups/cell elements remain charged by their existing
  # backend owner. Retain their actual owned record/list slots in this view.
  owned <- lapply(a$templates, function(x) list(lambda = x$lambda, nu = x$nu, group = NULL,
    cells = rep(list(NULL), length(x$cells)), N = x$N))
  testthat::expect_equal(vapply(owned, E, numeric(1)), 79 + 2 + pg)
  live <- list(incoming_lambda = numeric(nc), incoming_nu = lapply(pg, numeric),
    solved_lambda = numeric(nc), solved_nu = lapply(pg, numeric), templates = owned,
    run = run, comparison_subset = rep(list(NULL), nc),
    lambda_copy = numeric(nc), nu_list_copy = vector("list", ng),
    index_copies = list(integer(length(run$cell_ids)), integer(length(run$cell_ids))),
    nu_list_view = vector("list", length(run$group_ids)))
  forcing <- 8192 + 64 * nc + 256 * ng + 16 * sum(pg) + 16 * max(pg)^2
  testthat::expect_lte(E(live), forcing)
  # Separate-outcome replay keeps its original map/64*nc*(p+2) slab. Its
  # scalar template dimensions are the same, without the row workspace terms.
  maps <- 4096 + 64 * (nc + ng + nc)
  testthat::expect_lte(E(live), maps + 64 * nc * (max(pg) + 2))
})

testthat::test_that("each native suffix stops before an intervening scalar mutation or failure", {
  a <- zt_fixture(); a$active <- 3L
  first <- zt_native(a)
  testthat::expect_identical(first$group_ids, 1:2)
  testthat::expect_identical(first$next_gid, 3L)
  events <- integer()
  # This represents the original active scalar group/solve between runs. A
  # later group's metadata is changed here, after the first native call.
  scalar <- function(gid, fail = FALSE) {
    events <<- c(events, gid)
    a$groups[[4L]]$marker <<- -1
    if (fail) stop("intervening scalar failure")
  }
  scalar(first$next_gid)
  testthat::expect_null(zt_native(a, 4L))
  scalar(4L)
  next_run <- zt_native(a, 5L)
  testthat::expect_identical(next_run$group_ids, 5:8)
  testthat::expect_identical(events, c(3L, 4L))
  testthat::expect_error(scalar(3L, TRUE), "intervening scalar failure", fixed = TRUE)
})

testthat::test_that("run registration is checked once and altered scalar helpers keep the original path", {
  a <- zt_fixture(); scalar <- zt_ns("reconstruction_zero_template")
  run <- zt_ns("reconstruction_zero_template_run"); resolver <- zt_ns("gmm_require_native_dll")
  checks <- 0L
  testthat::local_mocked_bindings(gmm_require_native_dll = function() { checks <<- checks + 1L; resolver() }, .package = "cdid")
  call <- function(helper = scalar, cursor = 1L, active = integer()) run(
    a$templates, a$groups, a$cells, active, a$N, cursor, helper)
  testthat::expect_identical(call(), zt_reference(a))
  testthat::expect_identical(checks, 1L)
  testthat::expect_null(call(function(...) NULL))
  testthat::expect_null(call(cursor = 8L))
  testthat::expect_null(call(active = 1L))
  testthat::expect_identical(checks, 1L)
  testthat::local_mocked_bindings(reconstruction_zero_run_native_kernel = function() stop("registration failure"), .package = "cdid")
  testthat::expect_error(call(), "registration failure", fixed = TRUE)
})

zt_helpers <- function(file, wanted) {
  e <- new.env(parent = asNamespace("cdid"))
  for (x in parse(file = testthat::test_path(file))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  stopifnot(all(vapply(wanted, exists, logical(1), envir = e, inherits = FALSE)))
  e
}
zt_pure <- zt_helpers("test-reconstructed-outcome-pooling.R", c("rcp_ns", "rcp_gc", "rcp_case", "rcp_fit", "rcp_point", "rcp_reference"))
zt_aug <- zt_helpers("test-augmented-outcome-pooling.R", c("rap_ns", "rap_gc", "rap_case", "rap_fit", "rap_point", "rap_reference"))

# Count the unchanged scalar-helper calls in a test-local core body without
# replacing that helper, whose actual source identity is an admission condition.
zt_core <- function(fit, kind, rows, links, scalar = FALSE) {
  e <- new.env(parent = asNamespace("cdid")); events <- new.env(parent = emptyenv())
  events$scalar <- 0L; events$runs <- 0L; events$solves <- list(); e$zt_events <- events
  core <- zt_ns(paste0(kind, "_evaluate_if_chunk_core")); environment(core) <- e
  source <- paste(deparse(body(core), width.cutoff = 500L), collapse = "\n")
  call <- paste0("reconstruction_zero_template(zero_linear, gid, gid %in% active_groups, gf, ",
    if (kind == "augmented") "cells" else "core_fit$cells", ", N)")
  stopifnot(grepl(call, source, fixed = TRUE))
  source <- sub(call, paste0("{ zt_events$scalar <- zt_events$scalar + 1L; ", call, " }"), source, fixed = TRUE)
  body(core) <- parse(text = source)[[1L]]
  run <- zt_ns("reconstruction_zero_template_run")
  e$reconstruction_zero_template_run <- function(...) {
    out <- if (scalar) NULL else run(...)
    if (!is.null(out)) events$runs <- events$runs + 1L
    out
  }
  for (name in c("or_pooled_wls_gram_solve", "or_wls_gram_solve")) local({
    original <- zt_ns(name)
    assign(name, function(...) {
      f <- parent.frame()
      if (exists("tilde_lambda", f, inherits = FALSE) && exists("gid", f, inherits = FALSE))
        events$solves[[length(events$solves) + 1L]] <- list(gid = f$gid, rhs = f$rhs)
      original(...)
    }, e)
  })
  reader <- get("cp04_reader", fit$reader_ref, inherits = FALSE)
  pre <- zt_ns(paste0(kind, "_if_budget"))(fit, length(rows), links, FALSE)
  # Exercise this same scalar path even when optional separate-outcome forward
  # stage-one matrices would otherwise serve it. Both comparisons do this.
  pre$stage1_forward <- list(enabled = FALSE)
  args <- list(st = reader$state, X_F = if (kind == "augmented") zt_ns("augmented_backend_X_F")(fit) else
    zt_ns("imputation_backend_X_F")(fit, reader), pair_reader = function(r, pair, read_outcomes = TRUE)
      zt_ns("cp04_pair_rows")(reader, r, pair, read_outcomes), rows = rows, links = links,
    scale = "frame", details = FALSE, sets = pre)
  args <- if (kind == "augmented") c(args, list(backend = fit)) else c(args, list(core_fit = fit,
    recipe = fit$outcome_recipe, or_frame_index = fit$or_frame_index, trt_recipe = fit$treatment_recipe,
    est_method = fit$est_method))
  list(out = do.call(core, args), scalar = events$scalar, runs = events$runs, solves = events$solves)
}

testthat::test_that("native zero runs preserve scalar solve order and full F scientific references in both outcome modes", {
  specs <- list(list(kind = "imputation", method = "dr", imp = "separate", mode = "pooled", intercept = TRUE, slopes = TRUE),
    list(kind = "imputation", method = "reg", imp = "within_cohort", mode = "separate", intercept = FALSE, slopes = TRUE),
    list(kind = "imputation", method = "dr", imp = "across_cohorts", mode = "pooled", intercept = TRUE, slopes = FALSE),
    list(kind = "augmented", model = "pair_logit", imp = "separate", mode = "separate", intercept = TRUE, slopes = TRUE),
    list(kind = "augmented", model = "independent_period_logit", imp = "within_cohort", mode = "pooled", intercept = FALSE, slopes = TRUE),
    list(kind = "augmented", model = "independent_period_logit", imp = "across_cohorts", mode = "separate", intercept = TRUE, slopes = FALSE))
  total_runs <- 0L; saved_calls <- 0L
  for (s in specs) {
    h <- if (s$kind == "augmented") zt_aug else zt_pure
    a <- h[[if (s$kind == "augmented") "rap_case" else "rcp_case"]](intercept = s$intercept, slopes = s$slopes)
    fit <- if (s$kind == "augmented") h$rap_fit(a, imp = s$imp, model = s$model, mode = s$mode) else
      h$rcp_fit(a, method = s$method, pooling = s$imp, mode = s$mode)
    point <- h[[if (s$kind == "augmented") "rap_point" else "rcp_point"]]
    raw_ref <- h[[if (s$kind == "augmented") "rap_reference" else "rcp_reference"]]
    testthat::expect_equal(point(fit), raw_ref(a, fit), tolerance = 2e-8)
    N <- nrow(a$Y); links <- rev(seq_along(fit$links)); actual <- old <- matrix(0, N, length(links))
    for (start in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(start, min(N, start + 10L)))
      new <- zt_core(fit, s$kind, rows, links); ref <- zt_core(fit, s$kind, rows, links, scalar = TRUE)
      testthat::expect_identical(new$out, ref$out)
      testthat::expect_identical(new$solves, ref$solves)
      testthat::expect_lte(new$scalar, ref$scalar)
      total_runs <- total_runs + new$runs; saved_calls <- saved_calls + ref$scalar - new$scalar
      actual[rows, ] <- new$out$phi; old[rows, ] <- ref$out$phi
    }
    testthat::expect_identical(actual, old)
    testthat::expect_identical(crossprod(actual), crossprod(old))
    for (i in c(which(rowSums(!is.na(a$Y)) == 0L)[1L], which(rowSums(!is.na(a$Y)) == 1L)[1L])) {
      direction <- rep(-1, N); direction[i] <- N - 1
      numerical <- (raw_ref(a, fit, 1 + 1e-6 * direction) - raw_ref(a, fit, 1 - 1e-6 * direction)) / 2e-6
      testthat::expect_equal(unname(actual[i, ]), tail(numerical, length(links))[links], tolerance = 1e-5)
    }
  }
  testthat::expect_gt(total_runs, 0L)
  testthat::expect_gt(saved_calls, 0L)
})

testthat::test_that("the same metadata suffix preserves complete dated-layout IFs", {
  for (kind in c("imputation", "augmented")) {
    h <- if (kind == "augmented") zt_aug else zt_pure
    a <- h[[if (kind == "augmented") "rap_case" else "rcp_case"]](dated = TRUE, route_kind = kind)
    fit <- if (kind == "augmented") h$rap_fit(a, imp = "separate", mode = "separate") else
      h$rcp_fit(a, method = "dr", pooling = "separate", mode = "pooled")
    evaluate <- zt_ns(paste0("evaluate_", kind, "_if_chunk"))
    old <- function(rows, links) {
      testthat::local_mocked_bindings(reconstruction_zero_template_run = function(...) NULL, .package = "cdid")
      evaluate(fit, rows, links)
    }
    N <- nrow(a$Y); links <- rev(seq_along(fit$links)); actual <- reference <- matrix(0, N, length(links))
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- rev(seq.int(first, min(N, first + 10L)))
      got <- evaluate(fit, rows, links); ref <- old(rows, links)
      testthat::expect_identical(got, ref)
      actual[rows, ] <- got$phi; reference[rows, ] <- ref$phi
    }
    testthat::expect_identical(actual, reference)
    testthat::expect_identical(crossprod(actual), crossprod(reference))
  }
})

zt_public <- zt_helpers("test-d92-shared-if-reuse.R",
  c("d92_route_panel", "d92_counting_ns", "d92_run_public", "d92_strip", "d92_expect_parity"))
zt_public$api <- asNamespace("cdid")
zt_public$expect_equal <- testthat::expect_equal
zt_public$expect_identical <- testthat::expect_identical
testthat::test_that("public availability, full F covariance and reporting stay unchanged", {
  data <- zt_public$d92_route_panel()
  run <- function(ns) zt_public$d92_run_public(data, ns = ns,
    xformla = ~x, weightsname = "w", est_method = "dr", selection = "conditional_mtar",
    selection_est_method = "dr", selection_model = "pair_logit", selection_pooling = "within_cohort",
    imputation_pooling = "separate")
  new_ns <- zt_public$d92_counting_ns(); old_ns <- zt_public$d92_counting_ns()
  got <- run(new_ns)
  testthat::local_mocked_bindings(reconstruction_zero_template_run = function(...) NULL, .package = "cdid")
  old <- run(old_ns)
  zt_public$d92_expect_parity(new_ns$n$last_pipeline, old_ns$n$last_pipeline)
  testthat::expect_identical(new_ns$n$last_pipeline$support, old_ns$n$last_pipeline$support)
  testthat::expect_identical(new_ns$n$last_pipeline$reporting$target_meta, old_ns$n$last_pipeline$reporting$target_meta)
  testthat::expect_equal(got$att, old$att, tolerance = 1e-12)
  testthat::expect_equal(got$inffunc, old$inffunc, tolerance = 1e-12)
})
