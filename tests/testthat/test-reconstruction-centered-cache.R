# Exact witnessed-value reuse only. Factory recognition, state comparisons,
# scalar solves, products and the original finite witnesses remain active.
rcc_ns <- function(name) getFromNamespace(name, "cdid")
rcc_fixture <- local({
  e <- new.env(parent = asNamespace("cdid"))
  wanted <- c("rcp_ns", "rcp_gc", "rcp_case", "rcp_fit", "rcp_reference", "rcp_point")
  for (x in parse(file = testthat::test_path("test-reconstructed-outcome-pooling.R"))) if (is.call(x) &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) %in% wanted) eval(x, e)
  e
})
rcc_record <- function(n, columns) {
  X <- matrix(seq_len(n * length(columns)) / 7, n, length(columns), dimnames = list(NULL, columns))
  list(rr = seq_len(n), d_rr = X, s_rr = rep(1, n))
}

testthat::test_that("canonical derived-record inventory equals the independent payload walk", {
  E <- function(x) unname(rcc_ns("cp04_payload")(x)["entries"])
  peak <- rcc_ns("reconstruction_centered_record_peak")
  for (n in c(0L, 1L, 7L, 28L)) for (columns in list("x", c("x", "z"),
      c("\u03bb", paste0("long-", strrep("\u6f22", 200L))))) {
    rec <- rcc_record(n, columns)
    centered <- sweep(rec$d_rr, 2L, seq_along(columns) / 10, FUN = "-")
    testthat::expect_equal(peak(rec, centered, columns), 2 * E(rec) + 2 * E(centered) + 1)
    if (n) {
      decorated <- rec; names(decorated$s_rr) <- paste0("row-", seq_len(n))
      testthat::expect_identical(peak(decorated, centered, columns), Inf)
    }
    decorated <- rec; attr(decorated, "extra") <- TRUE
    testthat::expect_identical(peak(decorated, centered, columns), Inf)
  }
  testthat::expect_equal(E(list(values = list(), ref = NULL)), 42)
  testthat::expect_lte(42 + 64, 1024)
})

# Direct witness fixture isolates funding from production-origin selection.
# Mocking that optional selection never bypasses the old complete witness.
rcc_witness_case <- function(B = 28L, n = B) {
  cells <- list(list(group_id = 1L, mean_d = 0.25), list(group_id = 1L, mean_d = -0.5))
  groups <- list(list(p = 1L, cell_ids = 1:2, kept_columns = "x"))
  records <- list(rcc_record(n, "x"), rcc_record(n, "x"))
  factory <- rcc_ns("reconstruction_cached_cell_reader")
  reader <- factory(function(cid) records[[cid]], list(enabled = TRUE, cids = 1:2, entries = 1e6), 2L)
  reader(1L); reader(2L)
  view <- rcc_ns("reconstruction_score_cache_view")(reader)
  core <- structure(list(reference_ref = NULL, layout = NULL, pair_outcomes = NULL), class = "cdid_imputation_backend")
  state <- list(core = core, cells = cells, groups = groups, intercept = TRUE, cids = 1:2, N = B, B = B)
  list(reader = reader, view = view, state = state, records = records)
}

testthat::test_that("only derived retention declines when its old record slab cannot fund it", {
  witness <- rcc_ns("reconstruction_complete_cache_witness")
  centered <- rcc_ns("reconstruction_complete_cache_centered")
  # Arbitrary private reader state does not opt into the production reuse.
  a <- rcc_witness_case()
  plain <- witness(a$view, a$state)
  testthat::expect_null(centered(plain))
  testthat::expect_identical(plain$view, a$view)
  testthat::local_mocked_bindings(reconstruction_centered_cache_source = function(view) TRUE, .package = "cdid")
  admitted <- witness(a$view, a$state); values <- centered(admitted)
  testthat::expect_type(values, "list")
  for (j in 1:2) testthat::expect_identical(values[[a$view$position[j]]],
    sweep(a$records[[j]]$d_rr, 2L, a$state$cells[[j]]$mean_d, FUN = "-"))
  E <- function(x) unname(rcc_ns("cp04_payload")(x)["entries"])
  slab <- 1024 + 2 * (4 * a$state$B * 3 + 256 + 8)
  construction <- max(vapply(1:2, function(j) E(a$records[[j]]) + E(values[[a$view$position[j]]]), numeric(1)))
  testthat::expect_lte(E(a$view$cache) + E(admitted$centered) + construction, slab)
  small <- rcc_witness_case(1L)
  declined <- witness(small$view, small$state)
  testthat::expect_null(centered(declined))
  testthat::expect_identical(declined$view, small$view)
  testthat::expect_identical(rcc_ns("reconstruction_complete_cache_get")(declined$view, 1L), small$records[[1L]])
  bad <- admitted; bad$centered$values[[1L]][1L, 1L] <- 999
  testthat::expect_null(centered(bad))
  bad <- admitted; bad$centered$ref <- new.env()
  testthat::expect_null(centered(bad))
  # The original centered-subtraction overflow witness still fails before
  # any derived component could be used or a cache record could be skipped.
  overflow <- a; overflow$state$cells[[1L]]$mean_d <- -.Machine$double.xmax
  overflow$view$cache[[overflow$view$position[1L]]]$d_rr[] <- .Machine$double.xmax
  testthat::expect_null(witness(overflow$view, overflow$state))
  partial <- a; partial$view$seen[1L] <- FALSE
  testthat::expect_null(witness(partial$view, partial$state))
})

testthat::test_that("private source admission never forces reader promises or active bindings", {
  source <- rcc_ns("reconstruction_centered_cache_source")
  e <- new.env(parent = asNamespace("cdid"))
  delayedAssign("pair_reader", stop("private reader promise must remain unforced"), assign.env = e)
  read <- eval(quote(function(cid) NULL), e)
  testthat::expect_false(source(list(read = read)))
  e2 <- new.env(parent = asNamespace("cdid"))
  makeActiveBinding("pair_reader", function(...) stop("private active reader must remain unread"), e2)
  read2 <- eval(quote(function(cid) NULL), e2)
  testthat::expect_false(source(list(read = read2)))
})

rcc_run <- function(fit, rows, links, reuse) {
  getter <- rcc_ns("reconstruction_complete_cache_centered")
  recognizer <- rcc_ns("reconstruction_score_cache_view")
  pair <- rcc_ns("cp04_pair_rows"); formula <- rcc_ns("reconstruction_centered_record_peak")
  counts <- list(recognizer = 0L, raw = 0L, after_witness = 0L, served = 0L, record_inventory = 0L)
  witnessed <- FALSE
  testthat::local_mocked_bindings(
    reconstruction_complete_cache_centered = function(snapshot) {
      value <- getter(snapshot)
      if (!is.null(snapshot$centered)) witnessed <<- TRUE
      if (reuse && !is.null(value)) counts$served <<- counts$served + 1L
      if (reuse) value else NULL
    },
    reconstruction_score_cache_view = function(reader) {
      counts$recognizer <<- counts$recognizer + 1L; recognizer(reader)
    },
    reconstruction_centered_record_peak = function(...) {
      counts$record_inventory <<- counts$record_inventory + 1L; formula(...)
    },
    cp04_pair_rows = function(...) {
      counts$raw <<- counts$raw + 1L
      if (witnessed) counts$after_witness <<- counts$after_witness + 1L
      pair(...)
    }, .package = "cdid")
  out <- rcc_ns("evaluate_imputation_if_chunk")(fit, rows, links, scale = "frame", details = FALSE)
  list(out = out, counts = counts)
}

testthat::test_that("ordinary pure replay reuses exact centered values with every old recognition intact", {
  specs <- list(list(pool = "separate", intercept = TRUE, slopes = TRUE),
    list(pool = "within_cohort", intercept = TRUE, slopes = TRUE),
    list(pool = "across_cohorts", intercept = TRUE, slopes = TRUE),
    list(pool = "separate", intercept = FALSE, slopes = TRUE),
    list(pool = "within_cohort", intercept = TRUE, slopes = FALSE))
  for (s in specs) {
    a <- rcc_fixture$rcp_case(intercept = s$intercept, slopes = s$slopes, route_kind = "imputation")
    fit <- rcc_fixture$rcp_fit(a, method = "ipw", pooling = s$pool, mode = "pooled")
    links <- rev(seq_along(fit$links)); N <- nrow(a$Y)
    actual <- reference <- matrix(NA_real_, N, length(links))
    served <- 0L
    for (first in seq.int(1L, N, by = 11L)) {
      rows <- seq.int(first, min(N, first + 10L))
      got <- rcc_run(fit, rows, links, TRUE); old <- rcc_run(fit, rows, links, FALSE)
      testthat::expect_identical(got$out, old$out)
      testthat::expect_identical(got$counts$recognizer, old$counts$recognizer)
      testthat::expect_identical(got$counts$raw, old$counts$raw)
      testthat::expect_identical(got$counts$after_witness, 0L)
      if (got$counts$served) testthat::expect_gt(got$counts$record_inventory, 0L)
      served <- served + got$counts$served
      actual[rows, ] <- got$out$phi; reference[rows, ] <- old$out$phi
    }
    if (s$intercept && s$slopes) testthat::expect_gt(served, 0L) else testthat::expect_identical(served, 0L)
    testthat::expect_identical(actual, reference)
    testthat::expect_identical(crossprod(actual), crossprod(reference))
    testthat::expect_equal(rcc_fixture$rcp_point(fit), rcc_fixture$rcp_reference(a, fit), tolerance = 2e-8)
    for (i in c(which(rowSums(a$obs) == 0L)[1L], which(rowSums(a$obs) == 1L)[1L])) {
      delta <- rep(-1, N); delta[i] <- N - 1
      raw <- tail((rcc_fixture$rcp_reference(a, fit, 1 + 1e-6 * delta) -
        rcc_fixture$rcp_reference(a, fit, 1 - 1e-6 * delta)) / 2e-6, length(fit$links))
      testthat::expect_equal(actual[i, ], raw[links], tolerance = 1e-5)
    }
  }
})

testthat::test_that("one witnessed group returns its full needed closure without rebuilding or sorting", {
  # Literal old metadata algorithm; no numerical core or expected tensor is
  # replaced. Its structural/membership guards are copied unchanged.
  old <- function(snapshot, structural_cids) {
    st <- snapshot$state
    if (!is.integer(structural_cids) || !is.null(attributes(structural_cids)) ||
        anyNA(structural_cids) || any(structural_cids < 1L | structural_cids > length(st$cells))) return(NULL)
    if (any(snapshot$needed_position[structural_cids] == 0L)) return(NULL)
    gids <- unique(vapply(st$cells[structural_cids], `[[`, integer(1L), "group_id"))
    closure <- unique(c(structural_cids, unlist(lapply(st$groups[gids], `[[`, "cell_ids"))))
    at <- snapshot$needed_position[closure]
    st$cids[sort(unique(at[at > 0L]))]
  }
  make_snapshot <- function(p, two_groups = FALSE) {
    cells <- lapply(1:4, function(j) list(group_id = if (two_groups && j > 2L) 2L else 1L, mean_d = numeric(p)))
    groups <- if (two_groups) list(list(p = p, cell_ids = 1:2), list(p = p, cell_ids = 3:4)) else
      list(list(p = p, cell_ids = c(4L, 2L, 1L, 3L)))
    needed <- if (two_groups) c(4L, 2L, 1L, 3L) else c(3L, 1L, 4L)
    record <- rcc_record(4L, if (p) "x" else character(0L))
    reader <- rcc_ns("reconstruction_cached_cell_reader")(function(cid) record,
      list(enabled = TRUE, cids = needed, entries = 1e6), 4L)
    for (j in needed) reader(j)
    view <- rcc_ns("reconstruction_score_cache_view")(reader)
    state <- list(core = list(cells = cells, groups = groups), cells = cells, groups = groups,
      intercept = TRUE, cids = needed, N = 4L, B = 4L)
    snapshot <- rcc_ns("reconstruction_complete_cache_witness")(view, state)
    testthat::expect_type(snapshot, "list")
    snapshot
  }
  e <- new.env(parent = asNamespace("cdid")); sort_calls <- 0L
  e$sort <- function(...) { sort_calls <<- sort_calls + 1L; base::sort(...) }
  current <- rcc_ns("reconstruction_sparse_adjoint_ids"); environment(current) <- e; environment(old) <- e
  for (p in c(0L, 1L)) {
    snapshot <- make_snapshot(p)
    for (arm in list(1L, 3L, 4L, c(4L, 1L), c(1L, 1L), snapshot$state$cids)) {
      sort_calls <- 0L; expected <- old(snapshot, arm)
      testthat::expect_identical(sort_calls, 1L)
      sort_calls <- 0L
      testthat::expect_identical(current(snapshot, arm), expected)
      testthat::expect_identical(current(snapshot, arm), snapshot$state$cids)
      testthat::expect_identical(sort_calls, 0L)
    }
    changed_result <- current(snapshot, 1L); changed_result[1L] <- 99L
    testthat::expect_identical(snapshot$state$cids, c(3L, 1L, 4L))
    broken_state <- snapshot$state; broken_state$groups[[1L]]$cell_ids <- 1:3
    testthat::expect_null(rcc_ns("reconstruction_complete_cache_witness")(snapshot$view, broken_state))
    # Empty valid private arms stay empty. Unsupported, missing or foreign
    # structural IDs still stop before any closure construction or sorting.
    testthat::expect_identical(current(snapshot, integer(0L)), old(snapshot, integer(0L)))
    for (arm in list(1, NA_integer_, 0L, -1L, 5L, 2L, c(1L, 99L),
        structure(1L, names = "private"), matrix(1L, 1L, 1L))) {
      sort_calls <- 0L
      testthat::expect_null(current(snapshot, arm))
      testthat::expect_null(old(snapshot, arm))
      testthat::expect_identical(sort_calls, 0L)
    }
    multiple <- make_snapshot(p, TRUE)
    sort_calls <- 0L
    testthat::expect_identical(current(multiple, 1L), old(multiple, 1L))
    testthat::expect_identical(sort_calls, 2L)
  }
})
